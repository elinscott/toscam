!    *   ANGLE_SHIFT shifts angle ALPHA to lie between BETA and BETA+2PI.
!    * ANGLE_SHIFT_DEG shifts angle ALPHA to lie between BETA and BETA+360.
!    * ANGLE_TO_RGB returns a color on the perimeter of the color hexagon.
!    * ATAN4 computes the inverse tangent of the ratio Y / X.
!    * AXIS_LIMITS returns "nice" axis limits for a plot.
!    * BAR_CHECK computes the check digit for a barcode.
!    * BAR_CODE constructs the 113 character barcode from 11 digits.
!    * BAR_DIGIT_CODE_LEFT returns the 7 character left bar code for a digit.
!    * BAR_DIGIT_CODE_RIGHT returns the 7 character right bar code for a digit.
!    * BMI_ENGLISH computes the body mass index given English measurements.
!    * BMI_METRIC computes the body mass index given metric measurements.
!    * C4_ARGUMENT returns the argument of a C4.
!    * C4_MAGNITUDE returns the magnitude of a C4.
!    * C8_ARGUMENT returns the argument of a C8.
!    * C8_CUBE_ROOT returns the principal cube root of a C8.
!    * C8_I returns the C8 value of the aimaginary unit, i.
!    * C8_LE_L1 := X <= Y for C8 values, and the L1 norm.
!    * C8_LE_L2 := X <= Y for C8 values, and the L2 norm.
!    * C8_LE_LI := X <= Y for C8 values, and the L Infinity norm.
!    * C8_MAGNITUDE returns the magnitude of a C8.
!    * C8_NORM_L1 evaluates the L1 norm of a C8.
!    * C8_NORM_L2 evaluates the L2 norm of a C8.
!    * C8_NORM_LI evaluates the L-infinity norm of a C8.
!    * C8_NORMAL_01 returns a unit pseudonormal C8.
!    * C8_PRINT prints a C8, with an optional title.
!    * C8_SQRT returns the principal square root of a C8.
!    * C8_SWAP swaps two C8's.
!    * C8_UNIFORM_01 returns a unit pseudorandom C8.
!    * C8_ZERO returns the C8 value of 0.
!    * C8MAT_NINT rounds the entries of a C8MAT.
!    * C8MAT_PRINT prints a C8MAT.
!    * C8MAT_PRINT_OLD prints a C8MAT, with an optional title.
!    * C8MAT_PRINT_SOME prints some of a C8MAT.
!    * C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.
!    * C8VEC_INDICATOR sets a C8VEC to the indicator vector.
!    * C8VEC_NINT rounds the entries of a C8VEC.
!    * C8VEC_NORM_L2 returns the L2 norm of a C8VEC.
!    * C8VEC_PRINT prints a C8VEC, with an optional title.
!    * C8VEC_PRINT_SOME prints some of a C8VEC.
!    * C8VEC_SORT_A1 ascending sorts a C8VEC by L1 norm.
!    * C8VEC_SORT_A2 ascending sorts a C8VEC by L2 norm.
!    * C8VEC_SORT_AINF ascending sorts a C8VEC by L-infinity norm.
!    * C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.
!    * C8VEC_UNITY returns the N roots of unity.
!    * CH_IS_DIGIT is TRUE if C is a decimal digit.
!    * DEGREES_TO_RADIANS converts an angle measure from degrees to radians.
!    * E_CONSTANT returns the value of E.
!    * EULER_CONSTANT returns the value of the Euler-Mascheroni constant.
!    * FAC_DIV divides two quantities represented as prime factors.
!    * FAC_GCD finds the GCD of two products of prime factors.
!    * FAC_LCM finds the LCM of two products of prime factors.
!    * FAC_MUL multiplies two quantities represented as prime factors.
!    * FAC_PRINT prints a product of prime factors.
!    * FAC_TO_I4 converts a product of prime factors into an integer.
!    * FAC_TO_RAT converts a prime factorization into a rational value.
!    * FEET_TO_METERS converts a measurement in feet to meters.
!    * GAUSS_SUM evaluates a function that is the sum of Gaussians.
!    * GET_SEED returns a seed for the random number generator.
!    * GET_UNIT returns a free FORTRAN unit number.
!    * GRID1 finds grid points between X1 and X2 in N dimensions.
!    * GRID1N finds the I-th grid point between X1 and X2 in N dimensions.
!    * GRID2 computes grid points between X1 and X2 in N dimensions.
!    * GRID2N computes one grid point between X1 and X2 in N dimensions.
!    * GRID3 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
!    * GRID3N computes a parallelogram grid on 3 points in N dimensions.
!    * GRID4 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
!    * GRID4N computes a single point on a parallelogram grid in N space.
!    * I4_BIT_HI1 returns the position of the high 1 bit base 2 in an I4.
!    * I4_BIT_LO0 returns the position of the low 0 bit base 2 in an I4.
!    * I4_BIT_LO1 returns the position of the low 1 bit base 2 in an I4.
!    * I4_CHARACTERISTIC gives the characteristic for an I4.
!    * I4_DIV_ROUNDED computes the rounded result of I4 division.
!    * I4_DIVP returns the smallest multiple of J greater than or equal to I.
!    * I4_EVEN returns TRUE if an I4 is even.
!    * I4_GCD finds the greatest common divisor of two I4's.
!    * I4_GCDB finds the greatest common divisor of the form K**N of two I4's.
!    * I4_HUGE returns a "huge" I4.
!    * I4_HUGE_NORMALIZER returns the "normalizer" for I4_HUGE.
!    * I4_IS_POWER_OF_2 reports whether an I4 is a power of 2.
!    * I4_IS_PRIME reports whether an I4 is prime.
!    * I4_LCM computes the least common multiple of two I4's.
!    * I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
!    * I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
!    * I4_LOG_I4 returns the logarithm of an I4 to an I4 base.
!    * I4_LOG_R8 returns the logarithm of an I4 to an R8 base.
!    * I4_MANT computes the "mantissa" of a double precision number.
!    * I4_MODDIV breaks an I4 into a multiple of a divisor and remainder.
!    * I4_MODP returns the nonnegative remainder of I4 division.
!    * I4_ODD returns TRUE if an I4 is odd.
!    * I4_POWER returns the integer power of an I4.
!    * I4_SIGN evaluates the sign of an I4.
!    * I4_SWAP swaps two I4's.
!    * I4_SWAP3 swaps three I4's.
!    * I4_TO_ANGLE maps I4's to points on a circle.
!    * I4_TO_DIGITS_DECIMAL determines the last N decimal digits of an I4.
!    * I4_TO_FAC converts an I4 into a product of prime factors.
!    * I4_TO_HALTON computes one element of a leaped Halton subsequence.
!    * I4_TO_ISBN converts an I4 to an ISBN digit.
!    * I4_UNIFORM returns a scaled pseudorandom I4.
!    * I4_UNSWAP3 unswaps three I4's.
!    * I4_WALSH_1D evaluates the Walsh function.
!    * I4_WIDTH returns the "width" of an I4.
!    * I4_WRAP forces an I4 to lie between given limits by wrapping.
!    * I4_XOR calculates the exclusive OR of two I4's.
!    * I43MAT_FLIP_COLS swaps the columns of an I43MAT.
!    * I43MAT_FLIP_ROWS swaps the rows of an I43MAT.
!    * I4COL_COMPARE compares columns I and J of an I4COL.
!    * I4COL_FIND searches an I4COL for a particular column value.
!    * I4COL_FIND_ITEM searches an I4COL for a given scalar value.
!    * I4COL_FIND_PAIR_WRAP searches an I4COL for a pair of items.
!    * I4COL_SORT_A ascending sorts an I4COL.
!    * I4COL_SORT_D descending sorts an I4COL.
!    * I4COL_SORT2_A ascending sorts the elements of each column of an I4COL.
!    * I4COL_SORT2_D descending sorts elements of each column of an I4COL.
!    * I4COL_SORTED_SINGLETON_COUNT counts singletons in an I4COL.
!    * I4COL_SORTED_UNIQUE keeps unique elements in a sorted I4COL.
!    * I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
!    * I4COL_SWAP swaps columns J1 and J2 of an I4COL.
!    * I4INT_TO_R8INT maps an I4INT to an R8INT.
!    * I4MAT_ELIM carries out exact Gauss elimination on an I4MAT.
!    * I4MAT_FLIP_COLS swaps the columns of an I4MAT.
!    * I4MAT_FLIP_ROWS swaps the rows of an I4MAT.
!    * I4MAT_L1_INVERSE inverts a unit lower triangular I4MAT.
!    * I4MAT_MAX_INDEX returns the location of the maximum of an I4MAT.
!    * I4MAT_MIN_INDEX returns the location of the minimum of an I4MAT.
!    * I4MAT_MM multiplies two I4MAT's.
!    * I4MAT_PERM permutes the rows and columns of a square I4MAT.
!    * I4MAT_PERM_UNIFORM selects a random permutation of an I4MAT.
!    * I4MAT_PERM2 permutes the rows and columns of a rectangular I4MAT.
!    * I4MAT_PERM2_UNIFORM selects a random permutation of an I4MAT.
!    * I4MAT_PRINT prints an I4MAT.
!    * I4MAT_PRINT_SOME prints some of an I4MAT.
!    * I4MAT_RED divides out common factors in a row or column of an I4MAT.
!    * I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
!    * I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!    * I4MAT_U1_INVERSE inverts a unit upper triangular I4MAT.
!    * I4MAT_UNIFORM returns a scaled pseudorandom I4MAT.
!    * I4ROW_COMPARE compares two rows of an I4ROW.
!    * I4ROW_FIND_ITEM searches the rows of an I4ROW for a given value.
!    * I4ROW_FIND_PAIR_WRAP searches rows of an I4ROW for a pair of items.
!    * I4ROW_MAX returns the maximums of the rows of an I4ROW.
!    * I4ROW_MEAN returns the means of the rows of an I4ROW.
!    * I4ROW_MIN returns the minimums of the rows of an I4ROW.
!    * I4ROW_SORT_A ascending sorts the rows of an I4ROW.
!    * I4ROW_SORT_D descending sorts the rows of an I4ROW.
!    * I4ROW_SORT2_D descending sorts the elements of each row of an I4ROW.
!    * I4ROW_SORTED_UNIQUE keeps unique elements in an I4ROW.
!    * I4ROW_SORTED_UNIQUE_COUNT counts unique elements in an I4ROW.
!    * I4ROW_SUM returns the sums of the rows of an I4ROW.
!    * I4ROW_SWAP swaps two rows of an I4ROW.
!    * I4ROW_VARIANCE returns the variances of an I4ROW.
!    * I4VEC_AMAX returns the largest magnitude in an I4VEC.
!    * I4VEC_AMAX_INDEX returns the index of the largest magnitude in an I4VEC.
!    * I4VEC_AMIN returns the smallest magnitude in an I4VEC.
!    * I4VEC_AMIN_INDEX returns the index of the smallest magnitude in an I4VEC.
!    * I4VEC_AMINZ returns the smallest nonzero magnitude in an I4VEC.
!    * I4VEC_AMINZ_INDEX returns the smallest nonzero magnitude in an I4VEC.
!    * I4VEC_ASCEND_SUB computes the longest ascending subsequence of an I4VEC.
!    * I4VEC_ASCENDS determines if an I4VEC is (weakly) ascending.
!    * I4VEC_AXPY: Y(I) := Y(I) + A * X(I).
!    * I4VEC_BRACKET searches a sorted I4VEC for successive brackets of a value.
!    * I4VEC_COMPARE compares two I4VEC's.
!    * I4VEC_COPY copies an I4VEC.
!    * I4VEC_CUM computes the cumulutive sum of the entries of an I4VEC.
!    * I4VEC_DESCENDS determines if an I4VEC is decreasing.
!    * I4VEC_DIRECT_PRODUCT creates a direct product of I4VEC's.
!    * I4VEC_DIRECT_PRODUCT2 creates a direct product of I4VEC's.
!    * I4VEC_FRAC searches for the K-th smallest element in an I4VEC.
!    * I4VEC_GCD finds the greatest common divisor of an I4VEC.
!    * I4VEC_HEAP_A reorders an I4VEC into an ascending heap.
!    * I4VEC_HEAP_D reorders an I4VEC into an descending heap.
!    * I4VEC_HEAP_D_EXTRACT extracts the maximum value from a descending heap.
!    * I4VEC_HEAP_D_INSERT inserts a new I4 into a descending heap.
!    * I4VEC_HEAP_D_MAX returns the maximum value in a descending heap of integers.
!    * I4VEC_HISTOGRAM computes a histogram of the elements of an I4VEC.
!    * I4VEC_INDEX returns the first location of a given value in an I4VEC
!    * I4VEC_INDEX_DELETE_ALL deletes all occurrences of a value in an indexed sorted list.
!    * I4VEC_INDEX_DELETE_DUPES deletes duplicates from an indexed sorted I4VEC.
!    * I4VEC_INDEX_DELETE_ONE deletes one copy of an I4 from an indexed sorted list.
!    * I4VEC_INDEX_INSERT inserts an I4 into an indexed sorted I4VEC.
!    * I4VEC_INDEX_INSERT_UNIQUE inserts a unique I4 into an indexed sorted I4VEC.
!    * I4VEC_INDEX_ORDER sorts an I4VEC using an index vector.
!    * I4VEC_INDEX_SEARCH searches for an I4 in an indexed sorted I4VEC.
!    * I4VEC_INDEX_SORT_UNIQUE creates a sort index for an I4VEC.
!    * I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!    * I4VEC_INSERT inserts a value into an I4VEC.
!    * I4VEC_MASK_PRINT prints a masked I4VEC.
!    * I4VEC_MAX computes the maximum element of an I4VEC.
!    * I4VEC_MAX_INDEX computes the index of a maximum element of an I4VEC.
!    * I4VEC_MAX_INDEX_LAST returns the last maximal element location in an I4VEC
!    * I4VEC_MEAN returns the mean of an I4VEC.
!    * I4VEC_MEDIAN returns the median of an unsorted I4VEC.
!    * I4VEC_MERGE_A merges two ascending sorted I4VEC.
!    * I4VEC_MIN computes the minimum element of an I4VEC.
!    * I4VEC_MIN_INDEX computes the index of the minimum element of an I4VEC.
!    * I4VEC_NONZERO_COUNT counts the nonzero entries in an I4VEC.
!    * I4VEC_NONZERO_FIRST left-shifts all nonzeros in an I4VEC.
!    * I4VEC_ORDER_TYPE determines if an I4VEC is (non)strictly ascending/descending.
!    * I4VEC_PAIRWISE_PRIME checks whether an I4VEC's entries are pairwise prime.
!    * I4VEC_PART partitions an integer NVAL into N nearly equal parts.
!    * I4VEC_PART_QUICK_A reorders an I4VEC as part of a quick sort.
!    * I4VEC_PERMUTE permutes an I4VEC in place.
!    * I4VEC_PERMUTE_UNIFORM randomly permutes an I4VEC.
!    * I4VEC_PRINT prints an I4VEC.
!    * I4VEC_PRINT_SOME prints "some" of an I4VEC.
!    * I4VEC_PRODUCT returns the product of the entries of an I4VEC.
!    * I4VEC_RED divides out common factors in an I4VEC.
!    * I4VEC_REVERSE reverses the elements of an I4VEC.
!    * I4VEC_ROTATE rotates an I4VEC in place.
!    * I4VEC_RUN_COUNT counts runs of equal values in an I4VEC.
!    * I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC for a value.
!    * I4VEC_SEARCH_BINARY_D searches a descending sorted I4VEC for a value.
!    * I4VEC_SORT_BUBBLE_A ascending sorts an I4VEC using bubble sort.
!    * I4VEC_SORT_BUBBLE_D descending sorts an I4VEC using bubble sort.
!    * I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
!    * I4VEC_SORT_HEAP_D descending sorts an I4VEC using heap sort.
!    * I4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an I4VEC.
!    * I4VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an I4VEC.
!    * I4VEC_SORT_INSERT_A uses an ascending insertion sort on an I4VEC.
!    * I4VEC_SORT_INSERT_D uses a descending insertion sort on an I4VEC.
!    * I4VEC_SORT_QUICK_A ascending sorts an I4VEC using quick sort.
!    * I4VEC_SORT_SHELL_A ascending sorts an I4VEC using Shell's sort.
!    * I4VEC_SORTED_UNIQUE finds the unique elements in a sorted I4VEC.
!    * I4VEC_SORTED_UNIQUE_COUNT counts the unique elements in a sorted I4VEC.
!    * I4VEC_SORTED_UNIQUE_HIST histograms the unique elements of a sorted I4VEC.
!    * I4VEC_SPLIT "splits" an unsorted I4VEC based on a splitting value.
!    * I4VEC_SUM returns the sum of the entries of an I4VEC.
!    * I4VEC_SWAP swaps the entries of two I4VEC's.
!    * I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
!    * I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
!    * I4VEC_UNIQUE_COUNT counts the unique elements in an unsorted I4VEC.
!    * I4VEC_VALUE_NUM counts entries equal to a given value in an I4VEC.
!    * I4VEC_VALUE_INDEX indexesentries equal to a given value in an I4VEC.
!    * I4VEC_VARIANCE returns the variance of an I4VEC.
!    * I4VEC_WIDTH returns the "width" of an I4VEC.
!    * I4VEC_ZERO sets the entries of an I4VEC to 0.
!    * I4VEC2_COMPARE compares entries of an I4VEC2.
!    * I4VEC2_PRINT prints a pair of integer vectors.
!    * I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
!    * I4VEC2_SORT_D descending sorts a vector of pairs of integers.
!    * I4VEC2_SORTED_UNIQUE keeps the unique elements in a sorted I4VEC2.
!    * I8_HUGE returns a "huge" I8.
!    * I8_HUGE_NORMALIZER returns the "normalizer" for I8_HUGE.
!    * I8_XOR calculates the exclusive OR of two integers.
!    * I4I4_SORT_A ascending sorts a pair of integers.
!    * I4I4I4_SORT_A ascending sorts a triple of integers.
!    * IJ_NEXT returns the next matrix index.
!    * IJ_NEXT_GT returns the next matrix index, with the constraint that I < J.
!    * INDEX_BOX2_NEXT_2D produces indices on the surface of a box in 2D.
!    * INDEX_BOX2_NEXT_3D produces indices on the surface of a box in 3D.
!    * ISBN_CHECK checks an ISBN code.
!    * ISBN_FILL fills in a missing digit in an ISBN code.
!    * ISBN_TO_I4 converts an ISBN character into an integer.
!    * ISET2_COMPARE compares two I2 sets.
!    * ISET2_INDEX_INSERT_UNIQUE inserts unique I2 set values in an indexed sorted list.
!    * ISET2_INDEX_SEARCH searches for an I2 set value in an indexed sorted list.
!    * LCM_12N computes the least common multiple of the integers 1 through N.
!    * LVEC_PRINT prints an LVEC.
!    * LVEC_PRINT_SOME prints "some" of an LVEC.
!    * PAUSE_INPUT waits until an input character is entered.
!    * PERM_CHECK checks that a vector represents a permutation.
!    * PERM_CYCLE analyzes a permutation.
!    * PERM_FREE reports the number of unused items in a partial permutation.
!    * PERM_INVERSE inverts a permutation "in place".
!    * PERM_NEXT computes all of the permutations on N objects, one at a time.
!    * PERM_PRINT prints a permutation.
!    * PERM_UNIFORM selects a random permutation of N objects.
!    * POUNDS_TO_KILOGRAMS converts a measurement in pounds to kilograms.
!    * PRIME returns any of the first PRIME_MAX prime numbers.
!    * PRIME_GE returns the smallest prime greater than or equal to N.
!    * PRIMER computes the prime numbers up to a given limit.
!    * R4_EPSILON returns the R4 roundoff unit.
!    * R4_EXP computes the exponential function, avoiding overflow and underflow.
!    * R4_HUGE returns the largest legal R4.
!    * R4_TINY returns the smallest positive R4.
!    * R4_UNIFORM_01 returns a unit pseudorandom R4.
!    * R8_ABS returns the absolute value of an R8.
!    * R8_CAS returns the "casine" of an R8.
!    * R8_CEILING rounds an R8 "up" (towards +infinity) to the next integer.
!    * R8_CHOP chops an R8 to a given number of binary places.
!    * R8_CUBE_ROOT returns the cube root of an R8.
!    * R8_DIFF computes the difference of two R8's to a specified accuracy.
!    * R8_DIGIT returns a particular decimal digit of an R8.
!    * R8_EPSILON returns the R8 roundoff unit.
!    * R8_EXP computes the exponential function, avoiding overflow and underflow.
!    * R8_FLOOR rounds an R8 "down" (towards -infinity) to the next integer.
!    * R8_FRACTION returns the fraction part of an R8.
!    * R8_HUGE returns a very large R8.
!    * R8_IN_01 is TRUE if an R8 is in the range [0,1].
!    * R8_IS_INT determines if an R8 represents an integer value.
!    * R8_LOG_2 returns the logarithm base 2 of an R8.
!    * R8_LOG_10 returns the logarithm base 10 of an R8.
!    * R8_LOG_B returns the logarithm base B of an R8.
!    * R8_MANT computes the "mantissa" or "fraction part" of an R8.
!    * R8_MOD returns the remainder of R8 division.
!    * R8_MODP returns the nonnegative remainder of R8 division.
!    * R8_NINT returns the nearest integer to an R8.
!    * R8_NORMAL returns a scaled pseudonormal R8.
!    * R8_NORMAL_01 returns a unit pseudonormal R8.
!    * R8_PI returns the value of pi as an R8.
!    * R8_POWER computes the P-th power of an R8.
!    * R8_POWER_FAST computes an integer power of an R8.
!    * R8_PYTHAG computes sqrt ( A**2 + B**2 ), avoiding overflow and underflow.
!    * R8_ROUND2 rounds an R8 to a specified number of binary digits.
!    * R8_ROUNDB rounds an R8 to a given number of digits in a given base.
!    * R8_ROUNDX rounds an R8.
!    * R8_SIGN returns the sign of an R8.
!    * R8_CSQRT returns the complex square root of an R8.
!    * R8_SWAP swaps two R8's.
!    * R8_SWAP3 swaps three R8's.
!    * R8_TINY returns the smallest positive R8.
!    * R8_TO_R8_DISCRETE maps R to RD in [RMIN, RMAX] with NR possible values.
!    * R8_TO_DHMS converts decimal days into days, hours, minutes, seconds.
!    * R8_TO_I4 maps X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].
!    * R8_UNIFORM returns a scaled pseudorandom R8.
!    * R8_UNIFORM_01 returns a unit pseudorandom R8.
!    * R8_UNSWAP3 unswaps three R8's.
!    * R8_WALSH_1D evaluates the Walsh function.
!    * R82_CHEBY sets up the Chebyshev abscissas in an R8 interval.
!    * R82_DIST_L2 returns the L2 distance between a pair of R82's.
!    * R82_EQ == ( A1 == A2 ) for two R82's.
!    * R82_GE == ( A1 >= A2 ) for two R82's.
!    * R82_GT == ( A1 > A2 ) for two R82's.
!    * R82_LE == ( A1 <= A2 ) for two R82's.
!    * R82_LT == ( A1 < A2 ) for two R82's.
!    * R82_NE == ( A1 /= A2 ) for two R82's.
!    * R82_PRINT prints an R82.
!    * R82_SWAP swaps two R82 values.
!    * R82_UNIFORM returns a random R82 value in a given range.
!    * R82_UNIT_EUCLIDEAN_2D Euclidean normalizes an R82.
!    * R82POLY2_PRINT prints a second order polynomial in two variables.
!    * R82POLY2_TYPE analyzes a second order polynomial in two variables.
!    * R82POLY2_TYPE_PRINT prints the meaning of the output from R82POLY2_TYPE.
!    * R82VEC_MAX returns the maximum value in an R82VEC.
!    * R82VEC_MIN returns the minimum value in an R82VEC.
!    * R82VEC_ORDER_TYPE finds the order type of an R82VEC.
!    * R82VEC_PART_QUICK_A reorders an R82VEC as part of a quick sort.
!    * R82VEC_PERMUTE permutes an R82VEC in place.
!    * R82VEC_PRINT prints an R82VEC.
!    * R82VEC_SORT_HEAP_INDEX_A ascending index heaps an R82VEC.
!    * R82VEC_SORT_QUICK_A ascending sorts an R82VEC using quick sort.
!    * R83_CROSS computes the cross product of two R83's.
!    * R83_PRINT prints an R83.
!    * R83_STRIPLE computes the scalar triple product of three R83's.
!    * R83_SWAP swaps two R83's.
!    * R83_UNIT_EUCLIDEAN Euclidean normalizes an R83.
!    * R83_VTRIPLE computes the vector triple product of three R83's.
!    * R83VEC_UNIT_L2 makes each R83 in R83VEC have unit L2 norm.
!    * R84_UNIT_EUCLIDEAN_4D Euclidean normalizes an R84.
!    * R8BLOCK_EXPAND_LINEAR linearly interpolates new data into an R8BLOCK.
!    * R8BLOCK_PRINT prints an R8BLOCK.
!    * R8COL_COMPARE compares columns in an R8COL.
!    * R8COL_FIND seeks a column value in an R8COL.
!    * R8COL_INSERT inserts a column into an R8COL.
!    * R8COL_MAX returns the maximums in an R8COL.
!    * R8COL_MAX_INDEX returns the indices of column maximums in an R8COL.
!    * R8COL_MEAN returns the column means of an R8COL.
!    * R8COL_MIN returns the column minimums of an R8COL.
!    * R8COL_MIN_INDEX returns the indices of column minimums in an R8COL.
!    * R8COL_PART_QUICK_A reorders the columns of an R8COL.
!    * R8COL_PERMUTE permutes an R8COL in place.
!    * R8COL_SORT_HEAP_A ascending heapsorts an R8COL.
!    * R8COL_SORT_QUICK_A ascending quick sorts an R8COL.
!    * R8COL_SORTED_UNIQUE keeps unique elements in an R8COL.
!    * R8COL_SORTED_UNIQUE_COUNT counts unique elements in an R8COL.
!    * R8COL_SORTR_A ascending sorts one column of an R8COL, adjusting all columns.
!    * R8COL_SUM sums the columns of an R8COL.
!    * R8COL_SWAP swaps columns I and J of an R8COL.
!    * R8COL_TO_R8VEC converts an R8COL to an R8VEC.
!    * R8COL_VARIANCE returns the variances of an R8COL.
!    * R8R8_COMPARE compares two R8R8's.
!    * R8R8_PRINT prints an R8R8.
!    * R8R8R8_COMPARE compares two R8R8R8's.
!    * R8R8R8VEC_INDEX_INSERT_UNIQUE inserts a unique R8R8R8in an indexed sorted list.
!    * R8R8R8VEC_INDEX_SEARCH searches for an R8R8R8 value in an indexed sorted list.
!    * R8R8VEC_INDEX_INSERT_UNIQUE inserts a unique R8R8 in an indexed sorted list.
!    * R8R8VEC_INDEX_SEARCH searches for an R8R8 in an indexed sorted list.
!    * R8INT_TO_R8INT maps one R8INT to another.
!    * R8INT_TO_I4INT maps an R8INT to an integer interval.
!    * R8MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric matrix.
!    * R8MAT_CHOLESKY_SOLVE solves a Cholesky factored linear system A * x = b.
!    * R8MAT_DET computes the determinant of an R8MAT.
!    * R8MAT_DET_2D computes the determinant of a 2 by 2 R8MAT.
!    * R8MAT_DET_3D computes the determinant of a 3 by 3 R8MAT.
!    * R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
!    * R8MAT_DET_5D computes the determinant of a 5 by 5 R8MAT.
!    * R8MAT_DIAG_ADD_SCALAR adds a scalar to the diagonal of an R8MAT.
!    * R8MAT_DIAG_ADD_VECTOR adds a vector to the diagonal of an R8MAT.
!    * R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.
!    * R8MAT_DIAG_SET_SCALAR sets the diagonal of an R8MAT to a scalar value.
!    * R8MAT_DIAG_SET_VECTOR sets the diagonal of an R8MAT to a vector.
!    * R8MAT_EXPAND_LINEAR linearly interpolates new data into an R8MAT.
!    * R8MAT_EXPAND_LINEAR2 expands an R8MAT by linear interpolation.
!    * R8MAT_GIVENS_POST computes the Givens postmultiplier rotation matrix.
!    * R8MAT_GIVENS_PRE computes the Givens premultiplier rotation matrix.
!    * R8MAT_HESS approximates a Hessian matrix via finite differences.
!    * R8MAT_HOUSE_AXH computes A*H where H is a compact Householder matrix.
!    * R8MAT_HOUSE_FORM constructs a Householder matrix from its compact form.
!    * R8MAT_HOUSE_HXA computes H*A where H is a compact Householder matrix.
!    * R8MAT_HOUSE_POST computes a Householder post-multiplier matrix.
!    * R8MAT_HOUSE_PRE computes a Householder pre-multiplier matrix.
!    * R8MAT_IDENTITY stores the identity matrix in an R8MAT.
!    * R8MAT_IN_01 is TRUE if the entries of an R8MAT are in the range [0,1].
!    * R8MAT_INVERSE_2D inverts a 2 by 2 R8MAT using Cramer's rule.
!    * R8MAT_INVERSE_3D inverts a 3 by 3 R8MAT using Cramer's rule.
!    * R8MAT_INVERSE_4D inverts a 4 by 4 R8MAT using Cramer's rule.
!    * R8MAT_JAC estimates a dense jacobian matrix of the function FX.
!    * R8MAT_L_INVERSE inverts a lower triangular R8MAT.
!    * R8MAT_L_PRINT prints a lower triangular R8MAT.
!    * R8MAT_L_SOLVE solves a lower triangular linear system.
!    * R8MAT_L1_INVERSE inverts a double precision unit lower triangular R8MAT.
!    * R8MAT_LT_SOLVE solves a transposed lower triangular linear system.
!    * R8MAT_LU computes the LU factorization of a rectangular R8MAT.
!    * R8MAT_MAX returns the maximum entry of an R8MAT.
!    * R8MAT_MAX_INDEX returns the location of the maximum entry of an R8MAT.
!    * R8MAT_MAXCOL_MINROW gets the maximum column minimum row of an M by N R8MAT.
!    * R8MAT_MAXROW_MINCOL gets the maximum row minimum column of an M by N R8MAT.
!    * R8MAT_MIN returns the minimum entry of an M by N R8MAT.
!    * R8MAT_MIN_INDEX returns the location of the minimum entry of an R8MAT.
!    * R8MAT_MINCOL_MAXROW gets the minimum column maximum row of an M by N R8MAT.
!    * R8MAT_MINROW_MAXCOL gets the minimum row maximum column of an M by N R8MAT.
!    * R8MAT_MM multiplies two R8MAT's.
!    * R8MAT_MTV multiplies a transposed matrix times a vector
!    * R8MAT_MV multiplies a matrix times a vector.
!    * R8MAT_NINT rounds the entries of an R8MAT.
!    * R8MAT_NORM_EIS returns the EISPACK norm of an R8MAT.
!    * R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.
!    * R8MAT_NORM_L1 returns the matrix L1 norm of an R8MAT.
!    * R8MAT_NORM_L2 returns the matrix L2 norm of an R8MAT.
!    * R8MAT_NORM_LI returns the matrix L-infinity norm of an R8MAT.
!    * R8MAT_ORTH_UNIFORM returns a random orthogonal R8MAT.
!    * R8MAT_PLOT "plots" an R8MAT, with an optional title.
!    * R8MAT_PLOT_SYMBOL returns a symbol for a double precision number.
!    * R8MAT_POLY_CHAR computes the characteristic polynomial of an R8MAT.
!    * R8MAT_POWER computes a nonnegative power of an R8MAT.
!    * R8MAT_POWER_METHOD applies the power method to an R8MAT.
!    * R8MAT_PRINT prints an R8MAT.
!    * R8MAT_PRINT_SOME prints some of an R8MAT.
!    * R8MAT_PRINT2 prints an R8MAT.
!    * R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
!    * R8MAT_SOLVE_2D solves a 2 by 2 linear system using Cramer's rule.
!    * R8MAT_SOLVE_3D solves a 3 by 3 linear system using Cramer's rule.
!    * R8MAT_SOLVE2 computes the solution of an N by N linear system.
!    * R8MAT_SYMM_EIGEN returns a symmetric matrix with given eigensystem.
!    * R8MAT_SYMM_JACOBI applies Jacobi eigenvalue iteration to a symmetric matrix.
!    * R8MAT_TO_R8PLU factors a general R8MAT.
!    * R8MAT_TRACE computes the trace of an R8MAT.
!    * R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!    * R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!    * R8MAT_U_INVERSE inverts an upper triangular R8MAT.
!    * R8MAT_U1_INVERSE inverts a unit upper triangular R8MAT.
!    * R8MAT_UNIFORM fills an R8MAT with scaled pseudorandom numbers.
!    * R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
!    * R8MAT_VAND2 returns the N by N row Vandermonde matrix A.
!    * R8PLU_DET computes the determinant of an R8PLU matrix.
!    * R8PLU_INVERSE computes the inverse of an R8PLU matrix.
!    * R8PLU_MUL computes A * x using the PLU factors of A.
!    * R8PLU_SOL solves a linear system A*x=b from the PLU factors.
!    * R8PLU_TO_R8MAT recovers the matrix A that was factored by R8MAT_TO_R8PLU.
!    * R8POLY_DEGREE returns the degree of a polynomial.
!    * R8POLY_DERIV returns the derivative of a polynomial.
!    * R8POLY_LAGRANGE_0 evaluates the Lagrange factor at a point.
!    * R8POLY_LAGRANGE_1 evaluates the first derivative of the Lagrange factor.
!    * R8POLY_LAGRANGE_2 evaluates the second derivative of the Lagrange factor.
!    * R8POLY_LAGRANGE_COEF returns the coefficients of a Lagrange polynomial.
!    * R8POLY_LAGRANGE_FACTOR evaluates the polynomial Lagrange factor at a point.
!    * R8POLY_LAGRANGE_VAL evaluates the IPOL-th Lagrange polynomial.
!    * R8POLY_ORDER returns the order of a polynomial.
!    * R8POLY_PRINT prints out a polynomial.
!    * R8POLY_SHIFT adjusts the coefficients of a polynomial for a new argument.
!    * R8POLY_VAL_HORNER evaluates a polynomial using Horner's method.
!    * R8POLY2_EX finds the extremal point of a parabola determined by three points.
!    * R8POLY2_EX2 finds the extremal point of a parabola determined by three points.
!    * R8POLY2_ROOT returns the two roots of a quadratic polynomial.
!    * R8POLY2_RROOT returns the real parts of the roots of a quadratic polynomial.
!    * R8POLY2_VAL evaluates a parabola defined by three data values.
!    * R8POLY2_VAL2 evaluates a parabolic interpolant through tabular data.
!    * R8POLY3_ROOT returns the three roots of a cubic polynomial.
!    * R8POLY4_ROOT returns the four roots of a quartic polynomial.
!    * R8ROW_MAX returns the maximums of an R8ROW.
!    * R8ROW_MEAN returns the means of an R8ROW.
!    * R8ROW_MIN returns the minimums of an R8ROW.
!    * R8ROW_SORTED_UNIQUE_COUNT counts unique elements in an R8ROW.
!    * R8ROW_SUM returns the sums of the rows of an R8ROW.
!    * R8ROW_SWAP swaps two rows of an R8ROW.
!    * R8ROW_TO_R8VEC converts an R8ROW into an R8VEC.
!    * R8ROW_VARIANCE returns the variances of an R8ROW.
!    * R8SLMAT_PRINT prints a strict lower triangular R8MAT.
!    * R8VEC_01_TO_AB shifts and rescales an R8VEC to lie within given bounds.
!    * R8VEC_AB_TO_01 shifts and rescales an R8VEC to lie within [0,1].
!    * R8VEC_AB_TO_CD shifts and rescales an R8VEC from one interval to another.
!    * R8VEC_AMAX returns the maximum absolute value in an R8VEC.
!    * R8VEC_AMAX_INDEX returns the index of the maximum absolute value in an R8VEC.
!    * R8VEC_AMIN returns the minimum absolute value in an R8VEC.
!    * R8VEC_AMIN_INDEX returns the index of the minimum absolute value in an R8VEC.
!    * R8VEC_ASCENDS determines if an R8VEC is (weakly) ascending.
!    * R8VEC_BLEND prforms weighted interpolation of two R8VEC's.
!    * R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!    * R8VEC_BRACKET2 searches a sorted R8VEC for successive brackets of a value.
!    * R8VEC_BRACKET3 finds the interval containing or nearest a given value.
!    * R8VEC_CEILING rounds "up" (towards +infinity) entries of an R8VEC.
!    * R8VEC_CIRCULAR_VARIANCE returns the circular variance of an R8VEC.
!    * R8VEC_COMPARE compares two R8VEC's.
!    * R8VEC_CONVOLVE_CIRC returns the discrete circular convolution of two R8VEC's.
!    * R8VEC_COPY copies an R8VEC.
!    * R8VEC_CROSS_3D computes the cross product of two R8VEC's in 3D.
!    * R8VEC_CUM computes the cumulutive sum of the entries of an R8VEC.
!    * R8VEC_DIF computes coefficients for estimating the N-th derivative.
!    * R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
!    * R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
!    * R8VEC_DISTANCE returns the Euclidean distance between two R8VEC's.
!    * R8VEC_DISTINCT is true if the entries in an R8VEC are distinct.
!    * R8VEC_DOT finds the dot product of a pair of R8VEC's.
!    * R8VEC_EQ is true if two R8VECs are equal.
!    * R8VEC_EVEN returns N values, evenly spaced between ALO and AHI.
!    * R8VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
!    * R8VEC_EVEN2 linearly interpolates new numbers into an R8VEC.
!    * R8VEC_EVEN3 evenly interpolates new data into an R8VEC.
!    * R8VEC_EXPAND_LINEAR linearly interpolates new data into an R8VEC.
!    * R8VEC_FLOOR rounds "down" (towards -infinity) entries of an R8VEC.
!    * R8VEC_FRAC searches for the K-th smallest entry in an R8VEC.
!    * R8VEC_FRACTION returns the fraction parts of an R8VEC.
!    * R8VEC_GT == ( A1 > A2 ) for R8VEC's.
!    * R8VEC_HEAP_A reorders an R8VEC into an ascending heap.
!    * R8VEC_HEAP_D reorders an R8VEC into an descending heap.
!    * R8VEC_HISTOGRAM histograms an R8VEC.
!    * R8VEC_HOUSE_COLUMN defines a Householder premultiplier that "packs" a column.
!    * R8VEC_IN_01 is TRUE if the entries of an R8VEC are in the range [0,1].
!    * R8VEC_INDEX_DELETE_ALL deletes all occurrences of a value from an indexed sorted R8VEC.
!    * R8VEC_INDEX_DELETE_DUPES deletes duplicates from an indexed sorted R8VEC.
!    * R8VEC_INDEX_DELETE_ONE deletes one copy of a value from an indexed sorted R8VEC.
!    * R8VEC_INDEX_INSERT inserts a value in an indexed sorted R8VEC.
!    * R8VEC_INDEX_INSERT_UNIQUE inserts a unique value in an indexed sorted R8VEC.
!    * R8VEC_INDEX_ORDER sorts an R8VEC using an index vector.
!    * R8VEC_INDEX_SEARCH searches for a value in an indexed sorted R8VEC.
!    * R8VEC_INDEX_SORT_UNIQUE creates a sort index for an R8VEC.
!    * R8VEC_INDICATOR sets an R8VEC to the indicator vector.
!    * R8VEC_INSERT inserts a value into an R8VEC.
!    * R8VEC_IS_INT is TRUE if the entries of an R8VEC are integers.
!    * R8VEC_LENGTH returns the Euclidean length of an R8VEC.
!    * R8VEC_LT == ( A1 < A2 ) for R8VEC's.
!    * R8VEC_MASK_PRINT prints a masked R8VEC.
!    * R8VEC_MAX returns the maximum value in an R8VEC.
!    * R8VEC_MAX_INDEX returns the index of the maximum value in an R8VEC.
!    * R8VEC_MEAN returns the mean of an R8VEC.
!    * R8VEC_MEDIAN returns the median of an unsorted R8VEC.
!    * R8VEC_MIN returns the minimum value of an R8VEC.
!    * R8VEC_MIN_INDEX returns the index of the minimum value in an R8VEC.
!    * R8VEC_MIRROR_NEXT steps through all sign variations of an R8VEC.
!    * R8VEC_NINT rounds entries of an R8VEC.
!    * R8VEC_NORM_L1 returns the L1 norm of an R8VEC.
!    * R8VEC_NORM_L2 returns the L2 norm of an R8VEC.
!    * R8VEC_NORM_LI returns the L-infinity norm of an R8VEC.
!    * R8VEC_NORM_LP returns the LP norm of an R8VEC.
!    * R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!    * R8VEC_ORDER_TYPE determines if an R8VEC is (non)strictly ascending/descending.
!    * R8VEC_PART_QUICK_A reorders an R8VEC as part of a quick sort.
!    * R8VEC_PERMUTE permutes an R8VEC in place.
!    * R8VEC_PERMUTE_UNIFORM randomly permutes an R8VEC.
!    * R8VEC_POLARIZE decomposes an R8VEC into normal and parallel components.
!    * R8VEC_PRINT prints an R8VEC.
!    * R8VEC_PRINT_SOME prints "some" of an R8VEC.
!    * R8VEC_PRINT2 prints out an R8VEC.
!    * R8VEC_PRODUCT returns the product of the entries of an R8VEC.
!    * R8VEC_RANGE finds the range of Y's within a restricted X range.
!    * R8VEC_RANGE_2 updates a range to include a new array.
!    * R8VEC_READ reads an R8VEC from a file.
!    * R8VEC_READ_SIZE reads the size of an R8VEC from a file.
!    * R8VEC_REVERSE reverses the elements of an R8VEC.
!    * R8VEC_ROTATE "rotates" the entries of an R8VEC in place.
!    * R8VEC_SEARCH_BINARY_A searches an ascending sorted R8VEC.
!    * R8VEC_SORT_BUBBLE_A ascending sorts an R8VEC using bubble sort.
!    * R8VEC_SORT_BUBBLE_D descending sorts an R8VEC using bubble sort.
!    * R8VEC_SORT_HEAP_A ascending sorts an R8VEC using heap sort.
!    * R8VEC_SORT_HEAP_D descending sorts an R8VEC using heap sort.
!    * R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC.
!    * R8VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an R8VEC.
!    * R8VEC_SORT_HEAP_MASK_A: indexed heap ascending sort of a masked R8VEC.
!    * R8VEC_SORT_INSERT_A ascending sorts an R8VEC using an insertion sort.
!    * R8VEC_SORT_INSERT_INDEX_A ascending index sorts an R8VEC using insertion.
!    * R8VEC_SORT_INSERT_INDEX_D descending index sorts an R8VEC using insertion.
!    * R8VEC_SORT_QUICK_A ascending sorts an R8VEC using quick sort.
!    * R8VEC_SORT_SHELL_A ascending sorts an R8VEC using Shell's sort.
!    * R8VEC_SORT2_A ascending sorts an R8VEC and adjusts an associated R8VEC.
!    * R8VEC_SORTED_MERGE_A merges two ascending sorted R8VEC's.
!    * R8VEC_SORTED_NEAREST returns the nearest element in a sorted R8VEC.
!    * R8VEC_SORTED_SPLIT "splits" a sorted R8VEC, given a splitting value.
!    * R8VEC_SORTED_UNIQUE keeps the unique elements in a sorted R8VEC.
!    * R8VEC_SORTED_UNIQUE_COUNT counts the unique elements in a sorted R8VEC.
!    * R8VEC_SORTED_UNIQUE_HIST histograms the unique elements of a sorted R8VEC.
!    * R8VEC_SPLIT "splits" an unsorted R8VEC based on a splitting value.
!    * R8VEC_STD returns the standard deviation of an R8VEC.
!    * R8VEC_SUM returns the sum of the entries of an R8VEC.
!    * R8VEC_SWAP swaps the entries of two R8VECs.
!    * R8VEC_UNIFORM returns a scaled pseudorandom R8VEC
!    * R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!    * R8VEC_UNIQUE_COUNT counts the unique elements in an unsorted R8VEC.
!    * R8VEC_UNIT_EUCLIDEAN normalizes an R8VEC in the Euclidean norm.
!    * R8VEC_UNIT_SUM normalizes an R8VEC to have unit sum.
!    * R8VEC_VARIANCE returns the variance of an R8VEC.
!    * R8VEC_WRITE writes an R8VEC to a file.
!    * R8VEC_ZERO zeroes out an R8VEC.
!    * R8VEC2_COMPARE compares two entries in an R8VEC2.
!    * R8VEC2_PRINT prints an R8VEC2.
!    * R8VEC2_PRINT_SOME prints "some" of an R8VEC2.
!    * R8VEC2_SORT_A ascending sorts an R8VEC2.
!    * R8VEC2_SORT_D descending sorts an R8VEC2.
!    * R8VEC2_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC2.
!    * R8VEC2_SORTED_UNIQUE keeps unique elements in a sorted R8VEC2.
!    * R8VEC2_SORTED_UNIQUE_INDEX indexes unique elements in a sorted R8VEC2.
!    * R8VEC2_SUM_MAX_INDEX returns the index of the maximum sum of two R8VEC's.
!    * R8VEC3_PRINT prints an R8VEC3.
!    * RADIANS_TO_DEGREES converts an angle measure from radians to degrees.
!    * RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!    * eed = seed_input
!    * all random_seed ( )
!    * all random_seed ( size = seed_size )
!    * llocate ( seed_vector(seed_size) )
!    * f ( seed /= 0 ) then
!    * if ( debug ) then
!    * write ( *, '(a)' ) ' '
!    * write ( *, '(a)' ) 'RANDOM_INITIALIZE'
!    * write ( *, '(a,i20)' ) ' Initialize RANDOM_NUMBER, user SEED = ', seed
!    * end if
!    * lse
!    * call system_clock ( count, count_rate, count_max )
!    * seed = count
!    * if ( debug ) then
!    * write ( *, '(a)' ) ' '
!    * write ( *, '(a)' ) 'RANDOM_INITIALIZE'
!    * write ( *, '(a,i20)' ) ' Initialize RANDOM_NUMBER, arbitrary SEED = ', &
!    * seed
!    * end if
!    * nd if
!    * eed_vector(1:seed_size) = seed
!    * all random_seed ( put = seed_vector(1:seed_size) )
!    * eallocate ( seed_vector )
!    * o i = 1, warm_up
!    * call random_number ( harvest = t )
!    * nd do
!    * RAT_FACTOR factors a rational value into a product of prime factors.
!    * RICKEY evaluates Branch Rickey's baseball index.
!    * ROOTS_TO_R8POLY converts polynomial roots to polynomial coefficients.
!    * ROOTS_TO_I4POLY converts polynomial roots to polynomial coefficients.
!    * SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!    * TIMESTAMP prints the current YMDHMS date as a time stamp.
!    * TIMESTRING writes the current YMDHMS date into a string.
!    * TUPLE_NEXT2 computes the next element of an integer tuple space.
!    * TVEC_EVEN computes an evenly spaced set of angles between 0 and 2*PI.
!    * TVEC_EVEN2 computes an evenly spaced set of angles between 0 and 2*PI.
!    * TVEC_EVEN3 computes an evenly spaced set of angles between 0 and 2*PI.
!    * TVEC_EVEN_BRACKET computes evenly spaced angles between THETA1 and THETA2.
!    * TVEC_EVEN_BRACKET2 computes evenly spaced angles from THETA1 to THETA2.
!    * TVEC_EVEN_BRACKET3 computes evenly spaced angles from THETA1 to THETA2.
!    * UPC_CHECK_DIGIT returns the check digit of a UPC.


subroutine angle_shift ( alpha, beta, gamma )

!*****************************************************************************80
!
!! ANGLE_SHIFT shifts angle ALPHA to lie between BETA and BETA+2PI.
!
!  Discussion:
!
!    The input angle ALPHA is shifted by multiples of 2 * PI to lie
!    between BETA and BETA+2*PI.
!
!    The resulting angle GAMMA has all the same trigonometric function
!    values as ALPHA.
!
!  Modified:
!
!    14 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the angle to be shifted.
!
!    Input, real ( kind = 8 ) BETA, defines the lower endpoint of
!    the angle range.
!
!    Output, real ( kind = 8 ) GAMMA, the shifted angle.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  if ( alpha < beta ) then
    gamma = beta - mod ( beta - alpha, 2.0D+00 * pi ) + 2.0D+00 * pi
  else
    gamma = beta + mod ( alpha - beta, 2.0D+00 * pi )
  end if

  return
end
subroutine angle_shift_deg ( alpha, beta, gamma )

!*****************************************************************************80
!
!! ANGLE_SHIFT_DEG shifts angle ALPHA to lie between BETA and BETA+360.
!
!  Discussion:
!
!    The input angle ALPHA is shifted by multiples of 360 to lie
!    between BETA and BETA+360.
!
!    The resulting angle GAMMA has all the same trigonometric function
!    values as ALPHA.
!
!  Modified:
!
!    14 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, the angle to be shifted, in degrees.
!
!    Input, real ( kind = 8 ) BETA, defines the lower endpoint of
!    the angle range.
!
!    Output, real ( kind = 8 ) GAMMA, the shifted angle.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma

  if ( alpha < beta ) then
    gamma = beta - mod ( beta - alpha, 360.0D+00 ) + 360.0D+00
  else
    gamma = beta + mod ( alpha - beta, 360.0D+00 )
  end if

  return
end
subroutine angle_to_rgb ( angle, r, g, b )

!*****************************************************************************80
!
!! ANGLE_TO_RGB returns a color on the perimeter of the color hexagon.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle in the color hexagon.  
!    The sextants are defined by the following points:
!        0 degrees, 1, 0, 0, red;
!       60 degrees, 1, 1, 0, yellow;
!      120 degrees, 0, 1, 0, green;
!      180 degrees, 0, 1, 1, cyan;
!      240 degrees, 0, 0, 1, blue;
!      300 degrees, 1, 0, 1, magenta.
!
!    Output, real ( kind = 8 ) R, G, B, RGB specifications for the color 
!    that lies at the given angle, on the perimeter of the color hexagon.  One
!    value will be 1, and one value will be 0.
!
  implicit none

  real    ( kind = 8 ) angle
  real    ( kind = 8 ) angle2
  real    ( kind = 8 ) b
  real    ( kind = 8 ) g
  real    ( kind = 8 ), parameter :: degrees_to_radians = &
    3.141592653589793D+00 / 180.0D+00
  real    ( kind = 8 ) r

  angle = mod ( angle, 360.0D+00 )

  if ( angle < 0.0D+00 ) then
    angle = angle + 360.0D+00
  end if

  if ( angle <= 60.0D00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = 1.0D+00
    g = tan ( angle2 )
    b = 0.0D+00

  else if ( angle <= 120.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = cos ( angle2 ) / sin ( angle2 )
    g = 1.0D+00
    b = 0.0D+00

  else if ( angle <= 180.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = 1.0D+00
    b = tan ( angle2 )

  else if ( angle <= 240.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = cos ( angle2 ) / sin ( angle2 )
    b = 1.0D+00

  else if ( angle <= 300.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = tan ( angle2 )
    g = 0.0D+00
    b = 1.0D+00

  else if ( angle <= 360.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = 1.0D+00
    g = 0.0D+00
    b = cos ( angle2 ) / sin ( angle2 )

  end if

  return
end
function atan4 ( y, x )

!*****************************************************************************80
!
!! ATAN4 computes the inverse tangent of the ratio Y / X.
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
!    Input, real ( kind = 8 ) Y, X, two quantities which represent the 
!    tangent of an angle.  If Y is not zero, then the tangent is (Y/X).
!
!    Output, real ( kind = 8 ) ATAN4, an angle between 0 and 2 * PI, whose
!    tangent is (Y/X), and which lies in the appropriate quadrant so that 
!    the signs of its cosine and sine match those of X and Y.
!
  implicit none

  real    ( kind = 8 ) abs_x
  real    ( kind = 8 ) abs_y
  real    ( kind = 8 ) atan4
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) theta
  real    ( kind = 8 ) theta_0
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
!
!  Special cases:
!
  if ( x == 0.0D+00 ) then

    if ( 0.0D+00 < y ) then
      theta = pi / 2.0D+00
    else if ( y < 0.0D+00 ) then
      theta = 3.0D+00 * pi / 2.0D+00
    else if ( y == 0.0D+00 ) then
      theta = 0.0D+00
    end if

  else if ( y == 0.0D+00 ) then

    if ( 0.0D+00 < x ) then
      theta = 0.0D+00
    else if ( x < 0.0D+00 ) then
      theta = pi
    end if
!
!  We assume that ATAN2 is correct when both arguments are positive.
!
  else

    abs_y = abs ( y )
    abs_x = abs ( x )

    theta_0 = atan2 ( abs_y, abs_x )

    if ( 0.0D+00 < x .and. 0.0D+00 < y ) then
      theta = theta_0
    else if ( x < 0.0D+00 .and. 0.0D+00 < y ) then
      theta = pi - theta_0
    else if ( x < 0.0D+00 .and. y < 0.0D+00 ) then
      theta = pi + theta_0
    else if ( 0.0D+00 < x .and. y < 0.0D+00 ) then
      theta = 2.0D+00 * pi - theta_0
    end if

  end if

  atan4 = theta

  return
end
subroutine axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, nticks )

!*****************************************************************************80
!
!! AXIS_LIMITS returns "nice" axis limits for a plot.
!
!  Discussion:
!
!    The routine is given information about the range of a variable, and
!    the number of divisions desired.  It returns suggestions for
!    labeling a plotting axis for the variable, including the
!    starting and ending points, the length of a single division,
!    and a suggested tick marking for the axis.
!
!  Modified:
!
!    21 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XMIN, XMAX, the lower and upper values that
!    must be included on the axis.  XMIN must be less than XMAX.
!
!    Input, integer ( kind = 4 ) NDIVS, the number of divisions desired along
!    the axis.
!
!    Output, real ( kind = 8 ) PXMIN, PXMAX, the recommended lower and upper
!    axis bounds.  It will be the case that PXMIN <= XMIN < XMAX <= PXMAX.
!
!    Output, real ( kind = 8 ) PXDIV, the recommended size of a single division.
!
!    Output, integer ( kind = 4 ) NTICKS, a suggested number of ticks to use,
!    if subdividing each of the NDIVS divisions of the axis.
!
  implicit none

  integer ( kind = 4 ), parameter :: nsteps = 5

  real    ( kind = 8 ) best
  real    ( kind = 8 ) good
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) intlog
  integer ( kind = 4 ), dimension ( nsteps ) :: iticks = (/ 5, 4, 4, 5, 5 /)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ndivs
  integer ( kind = 4 ) nticks
  real    ( kind = 8 ) pxmax
  real    ( kind = 8 ) pxmax2
  real    ( kind = 8 ) pxmin
  real    ( kind = 8 ) pxmin2
  real    ( kind = 8 ) pxdiv
  real    ( kind = 8 ) pxdiv2
  real    ( kind = 8 ) r8_log_10
  real    ( kind = 8 ) reldif
  real    ( kind = 8 ), dimension ( nsteps ) :: steps = (/ &
    1.0D+00,  2.0D+00,  4.0D+00,  5.0D+00, 10.0D+00 /)
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) xmax
  real    ( kind = 8 ) xmin

  if ( xmin == xmax ) then
    xmin = xmin - 0.5D+00
    xmax = xmax + 0.5D+00
  else if ( xmax < xmin ) then
    temp = xmin
    xmin = xmax
    xmax = temp
  end if

  if ( ndivs <= 0 ) then
    ndivs = 5
  end if
!
!  Set RELDIF, the size of the X interval divided by the largest X.
!
  if ( xmax /= xmin ) then
    reldif = ( xmax - xmin ) / max ( abs ( xmax ), abs ( xmin ) )
  else
    reldif = 0.0D+00
  end if
!
!  If RELDIF tells us that XMIN and XMAX are extremely close,
!  do some simple things.
!
  if ( reldif < 0.00001D+00 ) then

    if ( xmax == 0.0D+00 ) then

      pxdiv = 1.0D+00

    else

      intlog = int ( r8_log_10 ( xmax ) )

      if ( intlog < 0 ) then
        intlog = intlog - 1
      end if

      pxdiv = 10.0D+00**intlog

      if ( 1.0D+00 < pxdiv ) then
        pxdiv = 1.0D+00
      end if

    end if

    nticks = 5
    pxmin = xmax - real ( ndivs / 2, kind = 8 ) * pxdiv
    pxmax = xmax + real ( ndivs - ( ndivs / 2 ), kind = 8 ) * pxdiv
!
!  But now handle the more general case, when XMIN and XMAX
!  are relatively far apart.
!
  else

    best = -999.0D+00
!
!  On second loop, increase INTLOG by 1.
!
    do j = 1, 2
!
!  Compute INTLOG, roughly the logarithm base 10 of the range
!  divided by the number of divisions.
!
      intlog = int ( r8_log_10 ( ( xmax - xmin ) &
        / real ( ndivs, kind = 8 ) ) ) + ( j - 1 )

      if ( xmax - xmin  < real ( ndivs, kind = 8 ) ) then
        intlog = intlog - 1
      end if
!
!  Now consider taking 1, 2, 4, 5 or 10 steps of size 10**INTLOG:
!
      do i = 1, nsteps
!
!  Compute the size of each step.
!
        pxdiv2 = steps(i) * 10.0D+00**intlog
!
!  Make sure NDIVS steps can reach from XMIN to XMAX, at least.
!
        if ( xmax <= xmin + ndivs * pxdiv2 ) then
!
!  Now decide where to start the axis.
!  Start the axis at PXMIN2, to the left of XMIN, and
!  representing a whole number of steps of size PXDIV2.
!
          if ( 0.0D+00 <= xmin ) then
            ival = int ( xmin / pxdiv2 )
          else
            ival = int ( xmin / pxdiv2 ) - 1
          end if

          pxmin2 = ival * pxdiv2
!
!  PXMAX2 is, of course, NDIVS steps above PXMIN2.
!
          pxmax2 = pxmin2 + ndivs * pxdiv2
!
!  Only consider going on if PXMAX2 is at least XMAX.
!
          if ( xmax <= pxmax2 ) then
!
!  Now judge this grid by the relative amount of wasted axis length.
!
            good = ( xmax - xmin ) / ( pxmax2 - pxmin2 )

            if ( best < good ) then
              best = good
              pxmax = pxmax2
              pxmin = pxmin2
              pxdiv = pxdiv2
              nticks = iticks(i)
            end if

          end if

        end if

      end do

    end do

  end if
!
!  If necessary, adjust the locations of PXMIN and PXMAX so that the
!  interval is more symmetric in containing XMIN through XMAX.
!
  do

    ilo = int ( xmin - pxmin ) / pxdiv
    ihi = int ( pxmax - xmax ) / pxdiv

    if ( ihi < ilo + 2 ) then
      exit
    end if

    pxmin = pxmin - pxdiv
    pxmax = pxmax - pxdiv

  end do

  return
end
subroutine bar_check ( digit, check )

!*****************************************************************************80
!
!! BAR_CHECK computes the check digit for a barcode.
!
!  Formula:
!
!    CHECK = SUM ( I = 1, 11, by 2's ) DIGIT(I)
!       + 3 * SUM ( I = 2, 10, by 2's ) DIGIT(I)
!
!    CHECK = MOD ( 10 - MOD ( CHECK, 10 ), 10 )
!
!  Modified:
!
!    23 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT(12), entries 1 through 11 of DIGIT 
!    contain the digits of the bar code.  Each entry must be between 0 and 9.
!    The 12th digit should be the check digit.
!
!    Output, integer ( kind = 4 ) CHECK, the correct check digit.  If the bar
!    code is correct, then DIGIT(12) should equal CHECK.
!
  implicit none

  integer ( kind = 4 ) check
  integer ( kind = 4 ) digit(12)
  integer ( kind = 4 ), parameter :: i4_ten = 10

  check = sum ( digit(1:11:2) ) + 3 * sum ( digit(2:10:2) )

  check = mod ( i4_ten - mod ( check, i4_ten ), i4_ten )

  return
end
subroutine bar_code ( digit, bar )

!*****************************************************************************80
!
!! BAR_CODE constructs the 113 character barcode from 11 digits.
!
!  Modified:
!
!    23 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) DIGIT(12).
!    On input, the first 11 entries of DIGIT contain a code to be
!    turned into a barcode.
!    On output, the 12-th entry of DIGIT is a check digit.
!
!    Output, character( len = 113 ) BAR, the bar code corresponding to the
!    digit information.
!
  implicit none

  character ( len = 113 ) bar
  integer ( kind = 4 ) check
  character ( len = 7 ) codel
  character ( len = 7 ) coder
  integer ( kind = 4 ) digit(12)
  integer ( kind = 4 ) i
!
!  9 character quiet zone.
!
  bar(1:9) = '000000000'
!
!  3 character guard pattern.
!
  bar(10:12) = '101'
!
!  7 character product category.
!
  call bar_digit_code_left ( digit(1), codel )
  bar(13:19) = codel
!
!  35 characters contain the 5 digit manufacturer code.
!
  do i = 1, 5
    call bar_digit_code_left ( digit(i+1), codel )
    bar(20+(i-1)*7:20+(i-1)*7+6) = codel
  end do
!
!  Center guard pattern.
!
  bar(55:59) = '01010'
!
!  35 characters contain the 5 digit product code.
!
  do i = 1, 5
    call bar_digit_code_right ( digit(i+6), coder )
    bar(60+(i-1)*7:60+(i-1)*7+6) = coder
  end do
!
!  Compute the check digit.
!
  call bar_check ( digit, check )
  digit(12) = check

  call bar_digit_code_right ( digit(12), coder )
  bar(95:101) = coder
!
!  Guard pattern.
!
  bar(102:104) = '101'
!
!  Quiet zone.
!
  bar(105:113) = '000000000'

  return
end
subroutine bar_digit_code_left ( digit, codel )

!*****************************************************************************80
!
!! BAR_DIGIT_CODE_LEFT returns the 7 character left bar code for a digit.
!
!  Example:
!
!    DIGIT = 3
!    CODEL = '0111101'
!
!  Modified:
!
!    26 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT, the digit, between 0 and 9.
!
!    Output, character ( len = 7 ) CODEL, the left code for the digit.
!
  implicit none

  character ( len = 7 ) codel
  integer ( kind = 4 ) digit

  if ( digit == 0 ) then
    codel = '0001101'
  else if ( digit == 1 ) then
    codel = '0011001'
  else if ( digit == 2 ) then
    codel = '0010011'
  else if ( digit == 3 ) then
    codel = '0111101'
  else if ( digit == 4 ) then
    codel = '0100011'
  else if ( digit == 5 ) then
    codel = '0110001'
  else if ( digit == 6 ) then
    codel = '0101111'
  else if ( digit == 7 ) then
    codel = '0111011'
  else if ( digit == 8 ) then
    codel = '0110111'
  else if ( digit == 9 ) then
    codel = '0001011'
  else
    codel = '???????'
  end if

  return
end
subroutine bar_digit_code_right ( digit, coder )

!*****************************************************************************80
!
!! BAR_DIGIT_CODE_RIGHT returns the 7 character right bar code for a digit.
!
!  Example:
!
!    DIGIT = 3
!    CODER = '1000010'
!
!  Modified:
!
!    26 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIGIT, the digit, between 0 and 9.
!
!    Output, character ( len = 7 ) CODER, the right code for the digit.
!
  implicit none

  character ( len = 7 ) coder
  integer ( kind = 4 ) digit

  if ( digit == 0 ) then
    coder = '1110010'
  else if ( digit == 1 ) then
    coder = '1100110'
  else if ( digit == 2 ) then
    coder = '1101100'
  else if ( digit == 3 ) then
    coder = '1000010'
  else if ( digit == 4 ) then
    coder = '1011100'
  else if ( digit == 5 ) then
    coder = '1001110'
  else if ( digit == 6 ) then
    coder = '1010000'
  else if ( digit == 7 ) then
    coder = '1000100'
  else if ( digit == 8 ) then
    coder = '1001000'
  else if ( digit == 9 ) then
    coder = '1110100'
  else
    coder = '???????'
  end if

  return
end
function bmi_english ( w_lb, h_ft, h_in )

!*****************************************************************************80
!
!! BMI_ENGLISH computes the body mass index given English measurements.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W_LB, the body weight in pounds.
!
!    Input, real ( kind = 8 ) H_FT, H_IN, the body height in feet and inches
!
!    Output, real ( kind = 8 ) BMI_ENGLISH, the body mass index.
!
  implicit none

  real    ( kind = 8 ) bmi_english
  real    ( kind = 8 ) bmi_metric
  real    ( kind = 8 ) feet_to_meters
  real    ( kind = 8 ) h_ft
  real    ( kind = 8 ) h_in
  real    ( kind = 8 ) h_m
  real    ( kind = 8 ) pounds_to_kilograms
  real    ( kind = 8 ) w_kg
  real    ( kind = 8 ) w_lb

  w_kg = pounds_to_kilograms ( w_lb )

  h_m = feet_to_meters ( h_ft + ( h_in / 12.0D+00 ) )

  bmi_english = bmi_metric ( w_kg, h_m )

  return
end
function bmi_metric ( w_kg, h_m )

!*****************************************************************************80
!
!! BMI_METRIC computes the body mass index given metric measurements.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W_KG, the body weight in kilograms.
!
!    Input, real ( kind = 8 ) H_M, the body height in meters.
!
!    Output, real ( kind = 8 ) BMI_METRIC, the body mass index.
!
  implicit none

  real    ( kind = 8 ) bmi_metric
  real    ( kind = 8 ) h_m
  real    ( kind = 8 ) w_kg

  bmi_metric = ( w_kg / h_m ) / h_m

  return
end
function c4_argument ( x )

!*****************************************************************************80
!
!! C4_ARGUMENT returns the argument of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex X, the value whose argument is desired.
!
!    Output, real C4_ARGUMENT, the argument of X.
!
  implicit none

  real c4_argument
  complex x

  if ( aimag ( x ) == 0.0E+00 .and. &
       real  ( x ) == 0.0E+00 ) then

    c4_argument = 0.0E+00

  else

    c4_argument = atan2 ( aimag ( x ), real ( x ) )

  end if

  return
end
function c4_magnitude ( x )

!*****************************************************************************80
!
!! C4_MAGNITUDE returns the magnitude of a C4.
!
!  Discussion:
!
!    A C4 is a complex ( kind = 4 ) value.
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex X, the value whose magnitude is desired.
!
!    Output, real C4_MAGNITUDE, the magnitude of X.
!
  implicit none

  real c4_magnitude
  complex x

  c4_magnitude = sqrt ( ( real ( x ) )**2 + ( aimag ( x ) )**2 )

  return
end
function c8_argument ( x )

!*****************************************************************************80
!
!! C8_ARGUMENT returns the argument of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the value whose argument is desired.
!
!    Output, real ( kind = 8 ) C8_ARGUMENT, the argument of X.
!
  implicit none

  real    ( kind = 8 ) c8_argument
  complex ( kind = 8 ) x

  if ( aimag ( x )           == 0.0D+00 .and. &
       real  ( x, kind = 8 ) == 0.0D+00 ) then

    c8_argument = 0.0D+00

  else

    c8_argument = atan2 ( aimag ( x ), real ( x ) )

  end if

  return
end
function c8_cube_root ( x )

!*****************************************************************************80
!
!! C8_CUBE_ROOT returns the principal cube root of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Modified:
!
!    08 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the number whose cube root is desired.
!
!    Output, complex ( kind = 8 ) C8_CUBE_ROOT, the cube root of X.
!
  implicit none

  real    ( kind = 8 ) argument
  real    ( kind = 8 ) c8_argument
  complex ( kind = 8 ) c8_cube_root
  real    ( kind = 8 ) c8_magnitude
  real    ( kind = 8 ) magnitude
  complex ( kind = 8 ) x

  argument = c8_argument ( x )

  magnitude = c8_magnitude ( x )

  if ( magnitude == 0.0D+00 ) then

    c8_cube_root = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  else

    c8_cube_root = magnitude**( 1.0D+00 / 3.0D+00 ) &
      * cmplx ( cos ( argument / 3.0D+00 ), &
                sin ( argument / 3.0D+00 ), kind = 8 )

  end if

  return
end
function c8_i ( )

!*****************************************************************************80
!
!! C8_I returns the C8 value of the aimaginary unit, i.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Modified:
!
!    27 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, complex ( kind = 8 ) C8_I, the value of complex i.
!
  implicit none

  complex ( kind = 8 ) c8_i

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )

  return
end
function c8_le_l1 ( x, y )

!*****************************************************************************80
!
!! C8_LE_L1 := X <= Y for C8 values, and the L1 norm.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    The L1 norm can be defined here as:
!
!      C8_NORM_L1(X) = abs ( real (X) ) + abs ( aimag (X) )
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, Y, the values to be compared.
!
!    Output, logical C8_LE_L1, is TRUE if X <= Y.
!
  implicit none

  logical c8_le_l1
  complex ( kind = 8 ) x
  complex ( kind = 8 ) y

  if ( abs ( real ( x, kind = 8 ) ) + abs ( aimag ( x ) ) <= &
       abs ( real ( y, kind = 8 ) ) + abs ( aimag ( y ) ) ) then
    c8_le_l1 = .true.
  else
    c8_le_l1 = .false.
  end if

  return
end
function c8_le_l2 ( x, y )

!*****************************************************************************80
!
!! C8_LE_L2 := X <= Y for C8 values, and the L2 norm.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    The L2 norm can be defined here as:
!
!      C8_NORM_L2(X) = sqrt ( ( real (X) )**2 + ( aimag (X) )**2 )
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, Y, the values to be compared.
!
!    Output, logical C8_LE_L2, is TRUE if X <= Y.
!
  implicit none

  logical c8_le_l2
  complex ( kind = 8 ) x
  complex ( kind = 8 ) y

  if ( ( real ( x, kind = 8 ) )**2 + ( aimag ( x ) )**2 <= &
       ( real ( y, kind = 8 ) )**2 + ( aimag ( y ) )**2 ) then
    c8_le_l2 = .true.
  else
    c8_le_l2 = .false.
  end if

  return
end
function c8_le_li ( x, y )

!*****************************************************************************80
!
!! C8_LE_LI := X <= Y for C8 values, and the L Infinity norm.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    The L Infinity norm can be defined here as:
!
!      C8_NORM_LI(X) = max ( abs ( real (X) ), abs ( aimag (X) ) )
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, Y, the values to be compared.
!
!    Output, logical C8_LE_LI, is TRUE if X <= Y.
!
  implicit none

  logical c8_le_li
  complex ( kind = 8 ) x
  complex ( kind = 8 ) y

  if ( max ( abs ( real ( x, kind = 8 ) ), abs ( aimag ( x ) ) ) <= &
       max ( abs ( real ( y, kind = 8 ) ), abs ( aimag ( y ) ) ) ) then
    c8_le_li = .true.
  else
    c8_le_li = .false.
  end if

  return
end
function c8_magnitude ( x )

!*****************************************************************************80
!
!! C8_MAGNITUDE returns the magnitude of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Modified:
!
!    21 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the value whose magnitude is desired.
!
!    Output, real ( kind = 8 ) C8_MAGNITUDE, the magnitude of X.
!
  implicit none

  real    ( kind = 8 ) c8_magnitude
  complex ( kind = 8 ) x

  c8_magnitude = sqrt ( ( real ( x, kind = 8 ) )**2 &
                   + ( aimag ( x ) )**2 )

  return
end
function c8_norm_l1 ( x )

!*****************************************************************************80
!
!! C8_NORM_L1 evaluates the L1 norm of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    Numbers of equal norm lie along diamonds centered at (0,0).
!
!    The L1 norm can be defined here as:
!
!      C8_NORM_L1(X) = abs ( real (X) ) + abs ( aimag (X) )
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the value whose norm is desired.
!
!    Output, real ( kind = 8 ) C8_NORM_L1, the norm of X.
!
  implicit none

  real    ( kind = 8 ) c8_norm_l1
  complex ( kind = 8 ) x

  c8_norm_l1 = abs ( real ( x, kind = 8 ) ) + abs ( aimag ( x ) )

  return
end
function c8_norm_l2 ( x )

!*****************************************************************************80
!
!! C8_NORM_L2 evaluates the L2 norm of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    Numbers of equal norm lie on circles centered at (0,0).
!
!    The L2 norm can be defined here as:
!
!      C8_NORM_L2(X) = sqrt ( ( real (X) )**2 + ( aimag ( X ) )**2 )
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the value whose norm is desired.
!
!    Output, real ( kind = 8 ) C8_NORM_L2, the 2-norm of X.
!
  implicit none

  real    ( kind = 8 ) c8_norm_l2
  complex ( kind = 8 ) x

  c8_norm_l2 = sqrt ( ( real ( x, kind = 8 ) )**2 &
                   + ( aimag ( x ) )**2 )

  return
end
function c8_norm_li ( x )

!*****************************************************************************80
!
!! C8_NORM_LI evaluates the L-infinity norm of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    Numbers of equal norm lie along squares whose centers are at (0,0).
!
!    The L-infinity norm can be defined here as:
!
!      C8_NORM_LI(X) = max ( abs ( real (X) ), abs ( aimag (X) ) )
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the value whose norm is desired.
!
!    Output, real ( kind = 8 ) C8_NORM_LI, the infinity norm of X.
!
  implicit none

  real    ( kind = 8 ) c8_norm_li
  complex ( kind = 8 ) x

  c8_norm_li = max ( abs ( real ( x, kind = 8 ) ), abs ( aimag ( x ) ) )

  return
end
function c8_normal_01 ( seed )

!*****************************************************************************80
!
!! C8_NORMAL_01 returns a unit pseudonormal C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, complex ( kind = 8 ) C8_NORMAL_01, a unit pseudornormal value.
!
  implicit none

  complex ( kind = 8 ) c8_normal_01
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) v1
  real    ( kind = 8 ) v2
  real    ( kind = 8 ) x_c
  real    ( kind = 8 ) x_r

  v1 = r8_uniform_01 ( seed )
  v2 = r8_uniform_01 ( seed )

  x_r = sqrt ( -2.0D+00 * log ( v1 ) ) * cos ( 2.0D+00 * pi * v2 )
  x_c = sqrt ( -2.0D+00 * log ( v1 ) ) * sin ( 2.0D+00 * pi * v2 )

  c8_normal_01 = cmplx ( x_r, x_c, kind = 8 )

  return
end
subroutine c8_print ( a, title )

!*****************************************************************************80
!
!! C8_PRINT prints a C8, with an optional title.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Modified:
!
!    28 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) A, the value to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!
  implicit none

  complex ( kind = 8 ) a
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a,2x,a,g14.6,a,g14.6,a)' ) &
      trim ( title ), '(', real ( a ), ',', aimag ( a ), ')'
  else
    write ( *, '(a,g14.6,a,g14.6,a)' ) &
      '(', real ( a ), ',', aimag ( a ), ')'
  end if

  return
end
function c8_sqrt ( x )

!*****************************************************************************80
!
!! C8_SQRT returns the principal square root of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, the number whose square root is desired.
!
!    Output, complex ( kind = 8 ) C8_SQRT, the square root of X.
!
  implicit none

  real    ( kind = 8 ) argument
  real    ( kind = 8 ) c8_argument
  real    ( kind = 8 ) c8_magnitude
  complex ( kind = 8 ) c8_sqrt
  real    ( kind = 8 ) magnitude
  complex ( kind = 8 ) x

  argument = c8_argument ( x )
  magnitude = c8_magnitude ( x )

  if ( magnitude == 0.0D+00 ) then

    c8_sqrt = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  else

    c8_sqrt = sqrt ( magnitude ) &
      * cmplx ( cos ( argument / 2.0D+00 ), &
                sin ( argument / 2.0D+00 ), kind = 8 )

  end if

  return
end
subroutine c8_swap ( x, y )

!*****************************************************************************80
!
!! C8_SWAP swaps two C8's.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, complex ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  complex ( kind = 8 ) x
  complex ( kind = 8 ) y
  complex ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
function c8_uniform_01 ( seed )

!*****************************************************************************80
!
!! C8_UNIFORM_01 returns a unit pseudorandom C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!    The angle should be uniformly distributed between 0 and 2 * PI,
!    the square root of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 8 ) C8_UNIFORM_01, a pseudorandom complex value.
!
  implicit none

  complex ( kind = 8 ) c8_uniform_01
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) r
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) theta

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if

  r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if

  theta = 2.0D+00 * pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

  c8_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

  return
end
function c8_zero ( )

!*****************************************************************************80
!
!! C8_ZERO returns the C8 value of 0.
!
!  Discussion:
!
!    A C8 is a complex ( kind = 8 ) value.
!
!  Modified:
!
!    23 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, complex ( kind = 8 ) C8_ZERO, the value of complex 0.
!
  implicit none

  complex ( kind = 8 ) c8_zero

  c8_zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  return
end
subroutine c8mat_nint ( m, n, a )

!*****************************************************************************80
!
!! C8MAT_NINT rounds the entries of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of complex ( kind = 8 ) values.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!
!    Input/output, complex ( kind = 8 ) A(M,N), the matrix to be NINT'ed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m
    do j = 1, n
      a(i,j) = cmplx ( nint ( real  ( a(i,j), kind = 8 ) ), &
                       nint ( aimag ( a(i,j) ) ), kind = 8 )
    end do
  end do

  return
end
subroutine c8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! C8MAT_PRINT prints a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of complex ( kind = 8 ) values.
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call c8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine c8mat_print_old ( m, n, a, title )

!*****************************************************************************80
!
!! C8MAT_PRINT_OLD prints a C8MAT, with an optional title.
!
!  Discussion:
!
!    A C8MAT is a matrix of complex ( kind = 8 ) values.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, complex ( kind = 8 ) A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  if ( maxval ( abs ( a(1:m,1:n) ) ) < 1000000.0D+00 ) then

    do jlo = 1, n, 3
      jhi = min ( jlo + 2, n )
      write ( *, '(a)' ) ' '
      write ( *, '(a3,3x,3(10x,i8,10x))' ) 'Col', ( j, j = jlo, jhi )
      write ( *, '(a3)' ) 'Row'
      do i = 1, m
        write ( *, '(i8,6f14.6)' ) i, a(i,jlo:jhi)
      end do
    end do

  else

    do jlo = 1, n, 3
      jhi = min ( jlo + 2, n )
      write ( *, '(a)' ) ' '
      write ( *, '(6x,3(10x,i8,10x))' ) ( j, j = jlo, jhi )
      write ( *, '(a)' ) ' '
      do i = 1, m
        write ( *, '(i8,6g14.6)' ) i, a(i,jlo:jhi)
      end do
    end do

  end if

  return
end
subroutine c8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C8MAT_PRINT_SOME prints some of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of complex ( kind = 8 ) values.
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 8 ) A(M,N), the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  character ( len = 20 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == zero ) then
          ctemp(j2) = '       0.0          '
        else if ( aimag ( a(i,j) ) == 0.0D+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = 8 )
        else
          write ( ctemp(j2), '(2g10.3)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,4a20)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine c8mat_uniform_01 ( m, n, seed, c )

!*****************************************************************************80
!
!! C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of complex ( kind = 8 ) values.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 8 ) C(M,N), the pseudorandom complex matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) j
  real    ( kind = 8 ) r
  integer ( kind = 4 ) k
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) theta

  do j = 1, n
    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge ( )
      end if

      r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge ( )
      end if

      theta = 2.0D+00 * pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

      c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

    end do

  end do

  return
end
subroutine c8vec_indicator ( n, a )

!*****************************************************************************80
!
!! C8VEC_INDICATOR sets a C8VEC to the indicator vector.
!
!  Discussion:
!
!    A C8VEC is a vector of complex ( kind = 8 ) values.
!
!    X(1:N) = ( 1-1i, 2-2i, 3-3i, 4-4i, ... )
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, complex ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = cmplx ( i, -i, kind = 8 )
  end do

  return
end
subroutine c8vec_nint ( n, a )

!*****************************************************************************80
!
!! C8VEC_NINT rounds the entries of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of complex ( kind = 8 ) values.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, complex ( kind = 8 ) A(N), the vector to be NINT'ed.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = cmplx ( nint ( real ( a(i), kind = 8 ) ), &
                   nint ( aimag ( a(i) ) ), kind = 8 )
  end do

  return
end
function c8vec_norm_l2 ( n, a )

!*****************************************************************************80
!
!! C8VEC_NORM_L2 returns the L2 norm of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of complex ( kind = 8 ) values.
!
!    The vector L2 norm is defined as:
!
!      C8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) conjg ( A(I) ) * A(I) ).
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, complex ( kind = 8 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) C8VEC_NORM_L2, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  real    ( kind = 8 ) c8vec_norm_l2

  c8vec_norm_l2 = sqrt ( sum ( conjg ( a(1:n) ) * a(1:n) ) )

  return
end
subroutine c8vec_print ( n, a, title )

!*****************************************************************************80
!
!! C8VEC_PRINT prints a C8VEC, with an optional title.
!
!  Discussion:
!
!    A C8VEC is a vector of complex ( kind = 8 ) values.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, complex ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,2g14.6)' ) i, a(i)
  end do

  return
end
subroutine c8vec_print_some ( n, x, i_lo, i_hi, title )

!*****************************************************************************80
!
!! C8VEC_PRINT_SOME prints some of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of complex ( kind = 8 ) values.
!
!  Modified:
!
!    18 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last entries
!    to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  character ( len = * ) title
  complex ( kind = 8 ) x(n)

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = max ( 1, i_lo ), min ( n, i_hi )
    write ( *, '(2x,i8,2x,2g14.6)' ) i, x(i)
  end do

  return
end
subroutine c8vec_sort_a1 ( n, x )

!*****************************************************************************80
!
!! C8VEC_SORT_A1 ascending sorts a C8VEC by L1 norm.
!
!  Discussion:
!
!    A C8VEC is a vector of complex ( kind = 8 ) values.
!
!    The L1 norm of A+Bi is abs(A) + abs(B).
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, complex ( kind = 8 ) X(N).
!    On input, an unsorted array.
!    On output, X has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  logical c8_le_l1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      call c8_swap ( x(i), x(j) )

    else if ( indx < 0 ) then

      if ( c8_le_l1 ( x(i), x(j) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c8vec_sort_a2 ( n, x )

!*****************************************************************************80
!
!! C8VEC_SORT_A2 ascending sorts a C8VEC by L2 norm.
!
!  Discussion:
!
!    A C8VEC is a vector of complex ( kind = 8 ) values.
!
!    The L2 norm of A+Bi is sqrt ( A**2 + B**2 ).
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, complex ( kind = 8 ) X(N).
!    On input, an unsorted array.
!    On output, X has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  logical c8_le_l2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      call c8_swap ( x(i), x(j) )

    else if ( indx < 0 ) then

      if ( c8_le_l2 ( x(i), x(j) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c8vec_sort_ainf ( n, x )

!*****************************************************************************80
!
!! C8VEC_SORT_AINF ascending sorts a C8VEC by L-infinity norm.
!
!  Discussion:
!
!    A C8VEC is a vector of complex ( kind = 8 ) values.
!
!    The L infinity norm of A+Bi is max ( abs ( A ), abs ( B ) ).
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, complex ( kind = 8 ) X(N).
!    On input, an unsorted array.
!    On output, X has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  logical c8_le_li
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      call c8_swap ( x(i), x(j) )

    else if ( indx < 0 ) then

      if ( c8_le_li ( x(i), x(j) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c8vec_uniform_01 ( n, seed, c )

!*****************************************************************************80
!
!! C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of complex ( kind = 8 ) values.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values to compute.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should 
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 8 ) C(N), the pseudorandom complex vector.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) r
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) theta

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge ( )
    end if

    r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge ( )
    end if

    theta = 2.0D+00 * pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

  end do

  return
end
subroutine c8vec_unity ( n, a )

!*****************************************************************************80
!
!! C8VEC_UNITY returns the N roots of unity.
!
!  Discussion:
!
!    A C8VEC is a vector of complex ( kind = 8 ) values.
!
!    X(1:N) = exp ( 2 * PI * (0:N-1) / N )
!
!    X(1:N)**N = ( (1,0), (1,0), ..., (1,0) ).
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, complex ( kind = 8 ) A(N), the N roots of unity.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) theta

  do i = 1, n
    theta = pi * real ( 2 * ( i - 1 ), kind = 8 ) / real ( n, kind = 8 )
    a(i) = cmplx ( cos ( theta ), sin ( theta ), kind = 8 )
  end do

  return
end
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT is TRUE if C is a decimal digit.
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
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, is TRUE if C is a digit.
!
  implicit none

  character c
  logical ch_is_digit

  if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
  end if

  return
end
function degrees_to_radians ( degrees )

!*****************************************************************************80
!
!! DEGREES_TO_RADIANS converts an angle measure from degrees to radians.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DEGREES, the angle measure in degrees.
!
!    Output, real ( kind = 8 ) DEGREES_TO_RADIANS, the angle measure in radians.
!
  implicit none

  real    ( kind = 8 ) degrees
  real    ( kind = 8 ) degrees_to_radians
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  degrees_to_radians = ( degrees / 180.0D+00 ) * pi

  return
end
function e_constant ( )

!*****************************************************************************80
!
!! E_CONSTANT returns the value of E.
!
!  Discussion:
!
!    "E" was named in honor of Euler, but is known as Napier's constant.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) E_CONSTANT, the base of the natural 
!    logarithm system.
!
  implicit none

  real    ( kind = 8 ) e_constant

  e_constant = 2.718281828459045D+00

  return
end
function euler_constant ( )

!*****************************************************************************80
!
!! EULER_CONSTANT returns the value of the Euler-Mascheroni constant.
!
!  Discussion:
!
!    The Euler-Mascheroni constant is often denoted by a lower-case
!    Gamma.  Gamma is defined as
!
!      Gamma = limit ( M -> Infinity ) ( Sum ( 1 <= N <= M ) 1 / N ) - Log ( M )
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) EULER_CONSTANT, the value of the 
!    Euler-Mascheroni constant.
!
  implicit none

  real    ( kind = 8 ) euler_constant

  euler_constant = 0.5772156649015328D+00

  return
end
subroutine fac_div ( prime_num, npower1, npower2, npower3 )

!*****************************************************************************80
!
!! FAC_DIV divides two quantities represented as prime factors.
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER1(PRIME_NUM), the powers of primes
!    in the representation of the first quantity.
!
!    Input, integer ( kind = 4 ) NPOWER2(PRIME_NUM), the powers of primes
!    in the representation of the second quantity.
!
!    Output, integer ( kind = 4 ) NPOWER3(PRIME_NUM), the powers of primes
!    in the representation of the quotient.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) npower1(prime_num)
  integer ( kind = 4 ) npower2(prime_num)
  integer ( kind = 4 ) npower3(prime_num)

  npower3(1:prime_num) = npower1(1:prime_num) - npower2(1:prime_num)

  return
end
subroutine fac_gcd ( prime_num, npower1, npower2, npower3 )

!*****************************************************************************80
!
!! FAC_GCD finds the GCD of two products of prime factors.
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER1(PRIME_NUM), the powers of primes
!    in the representation of the first quantity.  All the powers
!    must be nonnegative.
!
!    Input, integer ( kind = 4 ) NPOWER2(PRIME_NUM), the powers of primes
!    in the representation of the second quantity.  All the powers
!    must be nonnegative.
!
!    Output, integer ( kind = 4 ) NPOWER3(PRIME_NUM), the powers of primes
!    in the representation of the GCD.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) npower1(prime_num)
  integer ( kind = 4 ) npower2(prime_num)
  integer ( kind = 4 ) npower3(prime_num)

  do i = 1, prime_num

    if ( npower1(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAC_GCD - Fatal error!'
      write ( *, '(a)' ) '  One of the powers is negative!'
      stop
    end if

    if ( npower2(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAC_GCD - Fatal error!'
      write ( *, '(a)' ) '  One of the powers is negative!'
      stop
    end if

    npower3(i) = min ( npower1(i), npower2(i) )

  end do

  return
end
subroutine fac_lcm ( prime_num, npower1, npower2, npower3 )

!*****************************************************************************80
!
!! FAC_LCM finds the LCM of two products of prime factors.
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER1(PRIME_NUM), the powers of primes
!    in the representation of the first quantity.
!
!    Input, integer ( kind = 4 ) NPOWER2(PRIME_NUM), the powers of primes
!    in the representation of the second quantity.
!
!    Output, integer ( kind = 4 ) NPOWER3(PRIME_NUM), the powers of primes
!    in the representation of the LCM.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) npower1(prime_num)
  integer ( kind = 4 ) npower2(prime_num)
  integer ( kind = 4 ) npower3(prime_num)

  do i = 1, prime_num

    if ( npower1(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAC_LCM - Fatal error!'
      write ( *, '(a)' ) '  One of the powers is negative!'
      stop
    end if

    if ( npower2(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAC_LCM - Fatal error!'
      write ( *, '(a)' ) '  One of the powers is negative!'
      stop
    end if

    npower3(i) = max ( npower1(i), npower2(i) )

  end do

  return
end
subroutine fac_mul ( prime_num, npower1, npower2, npower3 )

!*****************************************************************************80
!
!! FAC_MUL multiplies two quantities represented as prime factors.
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER1(PRIME_NUM), the powers of primes
!    in the representation of the first quantity.
!
!    Input, integer ( kind = 4 ) NPOWER2(PRIME_NUM), the powers of primes
!    in the representation of the second quantity.
!
!    Output, integer ( kind = 4 ) NPOWER3(PRIME_NUM), the powers of primes
!    in the representation of the product.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) npower1(prime_num)
  integer ( kind = 4 ) npower2(prime_num)
  integer ( kind = 4 ) npower3(prime_num)

  npower3(1:prime_num) = npower1(1:prime_num) + npower2(1:prime_num)

  return
end
subroutine fac_print ( prime_num, npower )

!*****************************************************************************80
!
!! FAC_PRINT prints a product of prime factors.
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER(PRIME_NUM), the powers of primes
!    in the representation of the quantity.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) npower(prime_num)
  integer ( kind = 4 ) prime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Prime     Power'
  write ( *, '(a)' ) ' '
  do i = 1, prime_num
    if ( npower(i) /= 0 ) then
      write ( *, '(i8,2x,i8)' ) prime(i), npower(i)
    end if
  end do

  return
end
subroutine fac_to_i4 ( prime_num, npower, intval )

!*****************************************************************************80
!
!! FAC_TO_I4 converts a product of prime factors into an integer.
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER(PRIME_NUM), the powers of primes
!    in the representation of the quantity.  If any of these powers
!    are negative, then INTVAL will be set to 0.
!
!    Output, integer ( kind = 4 ) INTVAL, the integer represented by the 
!    product of the prime factors.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) npower(prime_num)
  integer ( kind = 4 ) prime

  intval = 1
  do i = 1, prime_num

    if ( npower(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FAC_TO_I4 - Fatal error!'
      write ( *, '(a)' ) '  One of the powers is negative!'
      stop
    end if

    intval = intval * prime(i)**npower(i)

  end do

  return
end
subroutine fac_to_rat ( prime_num, npower, top, bot )

!*****************************************************************************80
!
!! FAC_TO_RAT converts a prime factorization into a rational value.
!
!  Example:
!
!    Start with the prime factorization representation:
!
!      40/9 = 2**3 * 3**(-2) * 5
!
!    Input:
!
!      NPOWER = ( 3, -2, 1 )
!
!    Output:
!
!      TOP = 40 ( = 2**3 * 5**1 = PRIME(1)**3                 * PRIME(3)**1 )
!      BOT = 9  ( = 3**2        =               PRIME(2)**2 )
!
!  Modified:
!
!    16 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the index of the highest prime 
!    number used in the representations.
!
!    Input, integer ( kind = 4 ) NPOWER(PRIME_NUM).  NPOWER(I) is the power of
!    the I-th prime in the prime factorization.  NPOWER(I) may
!    be positive or negative.
!
!    Output, integer ( kind = 4 ) TOP, BOT, the top and bottom of a rational 
!    value.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) bot
  integer ( kind = 4 ) i
  integer ( kind = 4 ) npower(prime_num)
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) top

  top = 1
  bot = 1
  do i = 1, prime_num
    if ( 0 < npower(i) ) then
      top = top * prime(i)**npower(i)
    else if ( npower(i) < 0 ) then
      bot = bot * prime(i)**(-npower(i))
    end if
  end do

  return
end
function feet_to_meters ( ft )

!*****************************************************************************80
!
!! FEET_TO_METERS converts a measurement in feet to meters.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) FT, the length in feet.
!
!    Output, real ( kind = 8 ) FEET_TO_METERS, the corresponding 
!    length in meters.
!
  implicit none

  real    ( kind = 8 ) feet_to_meters
  real    ( kind = 8 ) ft

  feet_to_meters = 0.0254D+00 * 12.0D+00 * ft

  return
end
function gauss_sum ( dim_num, n, amplitude, center, width, x )

!*****************************************************************************80
!
!! GAUSS_SUM evaluates a function that is the sum of Gaussians.
!
!  Discussion:
!
!    Gauss_Sum(X) = Sum ( 1 <= J <= Ngauss ) Amplitude(I) * exp ( -Arg )
!
!    where
!
!      Arg = sum ( 1 <= I <= DIM_NUM ) ( ( ( X(I) - Center(I,J) ) / Width(J) )**2 )
!
!  Modified:
!
!    25 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of component Gaussian functions.
!
!    Input, real ( kind = 8 ) AMPLITUDE(N), CENTER(DIM_NUM,N), WIDTH(N),
!    the amplitude, center and width for the component Gaussian functions.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the point at which the function 
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) GAUSS_SUM, the value of the function.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n

  real    ( kind = 8 ) amplitude(n)
  real    ( kind = 8 ) arg
  real    ( kind = 8 ) center(dim_num,n)
  real    ( kind = 8 ) gauss_sum
  integer ( kind = 4 ) j
  real    ( kind = 8 ) width(n)
  real    ( kind = 8 ) x(dim_num)

  gauss_sum = 0.0D+00

  do j = 1, n

    arg = sum ( ( ( x(1:dim_num) - center(1:dim_num,j) ) / width(j) )**2 )

    gauss_sum = gauss_sum + amplitude(j) * exp ( -arg )

  end do

  return
end
subroutine get_seed ( seed )

!*****************************************************************************80
!
!! GET_SEED returns a seed for the random number generator.
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a pseudorandom seed value.
!
  implicit none

  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 ) today
  integer ( kind = 4 ) values(8)
  character ( len = 5 ) zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) /  11.0D+00
  temp = temp + real ( values(3) - 1, kind = 8 ) /  30.0D+00
  temp = temp + real ( values(5),     kind = 8 ) /  23.0D+00
  temp = temp + real ( values(6),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(7),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
  temp = temp                                    /   6.0D+00

  do while ( temp <= 0.0D+00 )
    temp = temp + 1.0D+00
  end do

  do while ( 1.0D+00 < temp )
    temp = temp - 1.0D+00
  end do

  seed = int ( real ( i4_huge ( ), kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer.
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == i4_huge ( ) ) then
    seed = seed - 1
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine grid1 ( dim_num, nstep, x1, x2, x )

!*****************************************************************************80
!
!! GRID1 finds grid points between X1 and X2 in N dimensions.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X1 and X2.
!
!    Input, integer ( kind = 4 ) NSTEP, the number of points to be generated.
!    NSTEP must be at least 2.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), the first and last
!    points, between which the equally spaced points are
!    to be computed.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,NSTEP), the set of equally spaced
!    points.  Each column of X represents one point, with X(*,1) = X1
!    and X(*,NSTEP) = X2.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep

  integer ( kind = 4 ) i
  real    ( kind = 8 ) x(dim_num,nstep)
  real    ( kind = 8 ) x1(dim_num)
  real    ( kind = 8 ) x2(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID1 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( nstep < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID1 - Fatal error!'
    write ( *, '(a)' ) '  NSTEP < 2.'
    write ( *, '(a,i8)' ) '  NSTEP = ', nstep
    stop
  end if

  do i = 1, nstep
    x(1:dim_num,i) = &
      ( real ( nstep - i,     kind = 8 ) * x1(1:dim_num)   &
      + real (         i - 1, kind = 8 ) * x2(1:dim_num) ) &
      / real ( nstep     - 1, kind = 8 )
  end do

  return
end
subroutine grid1n ( j, dim_num, nstep, x1, x2, x )

!*****************************************************************************80
!
!! GRID1N finds the I-th grid point between X1 and X2 in N dimensions.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J, the number of the desired point.
!    Normally J would be between 1 and NSTEP, but that is
!    not necessary.  Note that J = 1 returns X1 and J = NSTEP
!    returns X2.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X, 
!    X1 and X2.
!
!    Input, integer ( kind = 4 ) NSTEP, this is the number of equally
!    spaced points that are between X1 and X2.  NSTEP must
!    be at least 2, because X1 and X2 are always included
!    in the set of points.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), the first and last
!    points, between which the equally spaced points lie.
!
!    Output, real ( kind = 8 ) X(DIM_NUM), the J-th grid point between X1
!    and X2.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep

  integer ( kind = 4 ) j
  real    ( kind = 8 ) x(dim_num)
  real    ( kind = 8 ) x1(dim_num)
  real    ( kind = 8 ) x2(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID1N - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( nstep <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID1N - Fatal error!'
    write ( *, '(a)' ) '  NSTEP <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP = ', nstep
    stop
  end if

  x(1:dim_num) = ( real ( nstep - j,     kind = 8 ) * x1(1:dim_num) &
                 + real (         j - 1, kind = 8 ) * x2(1:dim_num) ) &
                 / real ( nstep     - 1, kind = 8 )

  return
end
subroutine grid2 ( j1, j2, dim_num, nstep, x1, x2, x )

!*****************************************************************************80
!
!! GRID2 computes grid points between X1 and X2 in N dimensions.
!
!  Discussion:
!
!    GRID2 computes grid points between X1 and X2 in N dimensions.
!
!    However, X1 need not be the first point computed, nor X2 the last.
!    The user must specify the steps on which X1 and X2 are passed
!    through.  These steps may even be outside the range of 1 through NSTEP.
!
!    We assume that a set of equally spaced points have
!    been drawn on the line through X1 and X2, and that
!    they have been numbered, with X1 labeled J1 and X2
!    labeled J2.  J1 or J2 may be between 1 and NSTEP,
!    in which case X1 or X2 will actually be returned in the
!    X array, but there is no requirement that J1 or J2
!    satisfy this condition.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J1, J2.  J1 specifies the step on which
!    X1 would be computed, and similarly for J2.  
!    J1 and J2 must be distinct.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X1 and X2.
!
!    Input, integer ( kind = 4 ) NSTEP, this is the number of equally
!    spaced points that are to be generated.
!    NSTEP should be at least 1.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), the points that define
!    the line along which the equally spaced points are generated, and
!    which may or may not be included in the set of computed points.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,NSTEP), the set of equally spaced
!    points.  Each column of X represents one point.
!    If 1 <= J1 <= NSTEP, then X(*,J1) = X1, and similarly for J2.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real    ( kind = 8 ) x(dim_num,nstep)
  real    ( kind = 8 ) x1(dim_num)
  real    ( kind = 8 ) x2(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID2 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( j1 == j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID2 - Fatal error!'
    write ( *, '(a)' ) '  J1 = J2, leading to zero denominator.'
    write ( *, '(a,i8)' ) '  J1 = ', j1
    write ( *, '(a,i8)' ) '  J2 = ', j2
    stop
  end if

  do j = 1, nstep
    do i = 1, dim_num
      x(i,j) = ( real ( j2 - j,      kind = 8 ) * x1(i)   &
               + real (      j - j1, kind = 8 ) * x2(i) ) &
               / real ( j2     - j1, kind = 8 )
    end do
  end do

  return
end
subroutine grid2n ( j, j1, j2, dim_num, x1, x2, x )

!*****************************************************************************80
!
!! GRID2N computes one grid point between X1 and X2 in N dimensions.
!
!  Discussion:
!
!    However, X1 need not be the first point computed, nor X2 the last.
!    The user must specify the steps on which X1 and X2 are passed through.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J, the J coordinate of the desired point.
!    Note that if J = J1, X will be returned as X1, and if
!    J = J2, X will be returned as X2.
!
!    Input, integer ( kind = 4 ) J1, J2.  J1 specifies the step on which
!    X1 would be computed, and similarly for J2.  That is,
!    we assume that a set of equally spaced points have
!    been drawn on the line through X1 and X2, and that
!    they have been numbered, with X1 labeled J1 and X2
!    labeled J2.  J1 and J2 must be distinct.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X1 and X2.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), the points that define
!    the line along which the equally spaced points are
!    generated, and which may or may not be included in the
!    set of computed points.
!
!    Output, real ( kind = 8 ) X(DIM_NUM).  X(I) is the J-th point from the
!    set of equally spaced points.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real    ( kind = 8 ) x(dim_num)
  real    ( kind = 8 ) x1(dim_num)
  real    ( kind = 8 ) x2(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID2N - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( j1 == j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID2N - Fatal error!'
    write ( *, '(a)' ) '  J1 = J2, leading to zero denominator.'
    write ( *, '(a,i8)' ) '  J1 = ', j1
    write ( *, '(a,i8)' ) '  J2 = ', j2
    stop
  end if

  do i = 1, dim_num
    x(i) = ( real ( j2 - j,      kind = 8 ) * x1(i)   &
           + real (      j - j1, kind = 8 ) * x2(i) ) &
           / real ( j2     - j1, kind = 8 )
  end do

  return
end
subroutine grid3 ( dim_num, nstep1, nstep2, x1, x2, x3, x )

!*****************************************************************************80
!
!! GRID3 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
!
!  Discussion:
!
!    The line between X1 and X2 will have NSTEP1 points generated along 
!    it, and the line between X1 and X3 will have NSTEP2 points generated
!    along it.
!
!    Fixing the second and third indices of X represents one point, with
!    the following special values:
!
!      X(*,1,1)      = X1
!      X(*,NSTEP1,1) = X2
!      X(*,1,NSTEP2) = X3.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X1, 
!    X2 and X3.
!
!    Input, integer ( kind = 4 ) NSTEP1, NSTEP2.  These are the number of
!    equally spaced points to generate in the first and second
!    directions.  NSTEP1 and NSTEP2 must be at least 2, because X1, X2 and
!    X3 are always included in the set of points.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), X3(DIM_NUM), the points
!    which define three corners of the parallelogram on
!    which the grid will be generated.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,NSTEP1,NSTEP2), the set of equally
!    spaced points.  
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep1
  integer ( kind = 4 ) nstep2

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) psi1
  real    ( kind = 8 ) psi2
  real    ( kind = 8 ) psi3
  real    ( kind = 8 ) x(dim_num,nstep1,nstep2)
  real    ( kind = 8 ) x1(dim_num)
  real    ( kind = 8 ) x2(dim_num)
  real    ( kind = 8 ) x3(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID3 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( nstep1 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID3 - Fatal error!'
    write ( *, '(a)' ) '  NSTEP1 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP1 = ', nstep1
    stop
  end if

  if ( nstep2 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID3 - Fatal error!'
    write ( *, '(a)' ) '  NSTEP2 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP2 = ', nstep2
    stop
  end if

  do j = 1, nstep1

    psi2 = real ( j - 1, kind = 8 ) &
         / real ( nstep1 - 1, kind = 8 )

    do k = 1, nstep2

      psi3 = real ( k - 1, kind = 8 ) &
           / real ( nstep2 - 1, kind = 8 )

      psi1 = 1.0D+00 - psi2 - psi3

      x(1:dim_num,j,k) = psi1 * x1(1:dim_num) &
                       + psi2 * x2(1:dim_num) &
                       + psi3 * x3(1:dim_num)

    end do
  end do

  return
end
subroutine grid3n ( j, k, dim_num, nstep1, nstep2, x1, x2, x3, x )

!*****************************************************************************80
!
!! GRID3N computes a parallelogram grid on 3 points in N dimensions.
!
!  Discussion:
!
!    The line between X1 and X2 will have NSTEP1
!    points generated along it, and the line between X1 and
!    X3 will have NSTEP2 points generated along it.
!
!    The following special values are:
!
!      J       K         X
!
!      1       1         X1
!      NSTEP1  1         X2
!      1       NSTEP2    X3
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J, K, the parallelogram coordinates
!    of the point.  J measures steps from X1 to X2, and
!    K measures steps from X1 to X3.  Normally, J would
!    be between 1 and NSTEP1, K between 1 and NSTEP2,
!    but this is not necessary.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X1, X2 and X3.
!
!    Input, integer ( kind = 4 ) NSTEP1, NSTEP2.  These are the number of
!    equally spaced points to generate in the first and second
!    directions.  NSTEP1 and NSTEP2 must be at least 2, because X1, X2 and
!    X3 are always included in the set of points.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), X3(DIM_NUM), the points
!    which define three corners of the parallelogram on
!    which the grid will be generated.
!
!    Output, real ( kind = 8 ) X(DIM_NUM), the point with coordinates (J,K)
!    from the the set of equally  spaced points.  
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep1
  integer ( kind = 4 ) nstep2

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) psi1
  real    ( kind = 8 ) psi2
  real    ( kind = 8 ) psi3
  real    ( kind = 8 ) x(dim_num)
  real    ( kind = 8 ) x1(dim_num)
  real    ( kind = 8 ) x2(dim_num)
  real    ( kind = 8 ) x3(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID3N - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( nstep1 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID3N - Fatal error!'
    write ( *, '(a)' ) '  NSTEP1 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP1 = ', nstep1
    stop
  end if

  if ( nstep2 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID3N - Fatal error!'
    write ( *, '(a)' ) '  NSTEP2 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP2 = ', nstep2
    stop
  end if

  psi2 = real ( j - 1, kind = 8 ) &
       / real ( nstep1 - 1, kind = 8 )

  psi3 = real ( k - 1, kind = 8 ) &
       / real ( nstep2 - 1, kind = 8 )

  psi1 = 1.0D+00 - psi2 - psi3

  x(1:dim_num) = psi1 * x1(1:dim_num) &
               + psi2 * x2(1:dim_num) &
               + psi3 * x3(1:dim_num)

  return
end
subroutine grid4 ( j1, j2, k1, k2, dim_num, nstep1, nstep2, x1, x2, x3, x )

!*****************************************************************************80
!
!! GRID4 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
!
!  Discussion:
!
!    Unlike GRID3, GRID4 does not necessarily place X1 at the
!    "origin" of the parallelogram, with X2 and X3 set at the
!    extreme J and K coordinates.  Instead, the user is free
!    to specify the J and K coordinates of the points, although
!    they are required to lie on a subparallelogram of the
!    larger one.
!
!    The line through X1 and X2 will have NSTEP1
!    points generated along it, and the line through X1 and
!    X3 will have NSTEP2 points generated along it.
!
!    If we aimagine that the
!    main parallelogram is drawn first, with coordinate
!    ranges 1 <= J <= NSTEP1 and 1 <= K <= NSTEP2, then
!    these indices determine the (J,K) coordinates of the
!    three points, namely:
!
!      X1 : (J1,K1)
!      X2 : (J2,K1)
!      X3 : (J1,K2)
!
!    Of course, we actually start with the points X1, X2,
!    and X3, and they define a parallelogram and a (J,K)
!    coordinate system over the plane containing them.  We
!    then are free to consider the parallelogram defined
!    by the three points (1,1), (NSTEP1,1) and (1,NSTEP2),
!    which may or may not contain any of the points X1, X2
!    and X3.
!
!    Assuming that the indices J1, J2, K1 and K2 are "within
!    bounds", the following special values will be computed:
!
!      X(*,J1,K1) = X1
!      X(*,J2,K1) = X2
!      X(*,J1,K2) = X3.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J1, J2, K1, K2, the indices.  
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X1, X2 and X3.
!
!    Input, integer ( kind = 4 ) NSTEP1, NSTEP2.  These are the number of
!    equally spaced points to generate in the first and second
!    directions.  NSTEP1 and NSTEP2 should be at least 1.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), X3(DIM_NUM), the points
!    which define three corners of the parallelogram on
!    which the grid will be generated.
!
!    Output, real ( kind = 8 ) X(DIM_NUM,NSTEP1,NSTEP2), the set of equally
!    spaced points.  Fixing the second and third indices
!    of X represents one point.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep1
  integer ( kind = 4 ) nstep2

  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real    ( kind = 8 ) psi1
  real    ( kind = 8 ) psi2
  real    ( kind = 8 ) psi3
  real    ( kind = 8 ) x(dim_num,nstep1,nstep2)
  real    ( kind = 8 ) x1(dim_num)
  real    ( kind = 8 ) x2(dim_num)
  real    ( kind = 8 ) x3(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( nstep1 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4 - Fatal error!'
    write ( *, '(a)' ) '  NSTEP1 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP1 = ', nstep1
    stop
  end if

  if ( nstep2 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4 - Fatal error!'
    write ( *, '(a)' ) '  NSTEP2 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP2 = ', nstep2
    stop
  end if

  if ( j1 == j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4 - Fatal error!'
    write ( *, '(a)' ) '  J1 = J2'
    write ( *, '(a,i8)' ) '  J1 = ', j1
    write ( *, '(a,i8)' ) '  J2 = ', j2
    stop
  end if

  if ( k1 == k2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4 - Fatal error!'
    write ( *, '(a)' ) '  K1 = K2'
    write ( *, '(a,i8)' ) '  K1 = ', k1
    write ( *, '(a,i8)' ) '  K2 = ', k2
    stop
  end if

  do j = 1, nstep1

    psi2 = real (  j - j1, kind = 8 ) &
         / real ( j2 - j1, kind = 8 )

    do k = 1, nstep2

      psi3 = real (  k - k1, kind = 8 ) &
           / real ( k2 - k1, kind = 8 )

      psi1 = 1.0D+00 - psi2 - psi3

      x(1:dim_num,j,k) = psi1 * x1(1:dim_num) &
                       + psi2 * x2(1:dim_num) &
                       + psi3 * x3(1:dim_num)

    end do
  end do

  return
end
subroutine grid4n ( j, j1, j2, k, k1, k2, dim_num, nstep1, nstep2, x1, x2, &
  x3, x )

!*****************************************************************************80
!
!! GRID4N computes a single point on a parallelogram grid in N space.
!
!  Discussion:
!
!    The computation is identical to that of GRID4, except that
!    only one point at a time is computed.
!
!    The line through X1 and X2 will have NSTEP1
!    points generated along it, and the line through X1 and
!    X3 will have NSTEP2 points generated along it.
!
!    The following special values will be computed:
!
!      J  K  X
!
!      J1 K1 X1
!      J2 K2 X2
!      J1 K2 X3
!
!    If we aimagine that the main parallelogram is drawn first, with 
!    coordinate ranges 1 <= J <= NSTEP1 and 1 <= K <= NSTEP2, then
!    the indices J and K determine the (J,K) coordinates of the
!    three points X1, X2, and X3, namely:
!
!      X1 : (J1,K1)
!      X2 : (J2,K1)
!      X3 : (J1,K2)
!
!    Of course, we actually start with the points X1, X2,
!    and X3, and they define a parallelogram and an (J,K)
!    coordinate system over the plane containing them.  We
!    then are free to consider the parallelogram defined
!    by the three points (1,1), (NSTEP1,1) and (1,NSTEP2),
!    which may or may not contain any of the points X1, X2
!    and X3.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J, the J coordinate of the point X.
!
!    Input, integer ( kind = 4 ) J1, J2.  See discussion.
!
!    Input, integer ( kind = 4 ) K, the K coordinate of the point X.
!
!    Input, integer ( kind = 4 ) K1, K2.  See discussion.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the points X, X1, X2 and X3.
!
!    Input, integer ( kind = 4 ) NSTEP1, NSTEP2.  These are the number of
!    equally spaced points generated in the first and second
!    directions.
!    NSTEP1 and NSTEP2 should be at least 1.
!
!    Input, real ( kind = 8 ) X1(DIM_NUM), X2(DIM_NUM), X3(DIM_NUM), the points
!    which define three corners of the parallelogram on
!    which the grid will be generated.
!
!    Output, real ( kind = 8 ) X(DIM_NUM), the point whose parallelogram
!    coordinates are (J,K).
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) nstep1
  integer ( kind = 4 ) nstep2

  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real    ( kind = 8 ) psi1
  real    ( kind = 8 ) psi2
  real    ( kind = 8 ) psi3
  real    ( kind = 8 ) x(dim_num)
  real    ( kind = 8 ) x1(dim_num)
  real    ( kind = 8 ) x2(dim_num)
  real    ( kind = 8 ) x3(dim_num)

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4N - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if

  if ( nstep1 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4N - Fatal error!'
    write ( *, '(a)' ) '  NSTEP1 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP1 = ', nstep1
    stop
  end if

  if ( nstep2 <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4N - Fatal error!'
    write ( *, '(a)' ) '  NSTEP2 <= 1.'
    write ( *, '(a,i8)' ) '  NSTEP2 = ', nstep2
    stop
  end if

  if ( j1 == j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4N - Fatal error!'
    write ( *, '(a)' ) '  J1 = J2'
    write ( *, '(a,i8)' ) '  J1 = ', j1
    write ( *, '(a,i8)' ) '  J2 = ', j2
    stop
  end if

  if ( k1 == k2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID4N - Fatal error!'
    write ( *, '(a)' ) '  K1 = K2'
    write ( *, '(a,i8)' ) '  K1 = ', k1
    write ( *, '(a,i8)' ) '  K2 = ', k2
    stop
  end if

  psi2 = real ( j  - j1, kind = 8 ) &
       / real ( j2 - j1, kind = 8 )

  psi3 = real ( k  - k1, kind = 8 ) &
       / real ( k2 - k1, kind = 8 )

  psi1 = 1.0D+00 - psi2 - psi3

  x(1:dim_num) = psi1 * x1(1:dim_num) &
               + psi2 * x2(1:dim_num) &
               + psi3 * x3(1:dim_num)

  return
end
function i4_bit_hi1 ( n )

!*****************************************************************************80
!
!! I4_BIT_HI1 returns the position of the high 1 bit base 2 in an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!       N    Binary    Hi 1
!    ----    --------  ----
!       0           0     0
!       1           1     1
!       2          10     2
!       3          11     2
!       4         100     3
!       5         101     3
!       6         110     3
!       7         111     3
!       8        1000     4
!       9        1001     4
!      10        1010     4
!      11        1011     4
!      12        1100     4
!      13        1101     4
!      14        1110     4
!      15        1111     4
!      16       10000     5
!      17       10001     5
!    1023  1111111111    10
!    1024 10000000000    11
!    1025 10000000001    11
!
!  Modified:
!
!    01 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be measured.
!    N should be nonnegative.  If N is nonpositive, the function
!    will always be 0.
!
!    Output, integer ( kind = 4 ) I4_BIT_HI1, the position of the highest bit.
!
  implicit none

  integer ( kind = 4 ) bit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_bit_hi1
  integer ( kind = 4 ) n

  i = n
  bit = 0

  do

    if ( i <= 0 ) then
      exit
    end if

    bit = bit + 1
    i = i / 2

  end do

  i4_bit_hi1 = bit

  return
end
function i4_bit_lo0 ( n )

!*****************************************************************************80
!
!! I4_BIT_LO0 returns the position of the low 0 bit base 2 in an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!       N    Binary    Lo 0
!    ----    --------  ----
!       0           0     1
!       1           1     2
!       2          10     1
!       3          11     3
!       4         100     1
!       5         101     2
!       6         110     1
!       7         111     4
!       8        1000     1
!       9        1001     2
!      10        1010     1
!      11        1011     3
!      12        1100     1
!      13        1101     2
!      14        1110     1
!      15        1111     5
!      16       10000     1
!      17       10001     2
!    1023  1111111111     1
!    1024 10000000000     1
!    1025 10000000001     1
!
!  Modified:
!
!    16 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be measured.
!    N should be nonnegative.
!
!    Output, integer ( kind = 4 ) I4_BIT_LO0, the position of the low 1 bit.
!
  implicit none

  integer ( kind = 4 ) bit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_bit_lo0
  integer ( kind = 4 ) n

  bit = 0
  i = n

  do

    bit = bit + 1
    i2 = i / 2

    if ( i == 2 * i2 ) then
      exit
    end if

    i = i2

  end do

  i4_bit_lo0 = bit

  return
end
function i4_bit_lo1 ( n )

!*****************************************************************************80
!
!! I4_BIT_LO1 returns the position of the low 1 bit base 2 in an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!       N    Binary    Lo 1
!    ----    --------  ----
!       0           0     0
!       1           1     1
!       2          10     2
!       3          11     1
!       4         100     3
!       5         101     1
!       6         110     2
!       7         111     1
!       8        1000     4
!       9        1001     1
!      10        1010     2
!      11        1011     1
!      12        1100     3
!      13        1101     1
!      14        1110     2
!      15        1111     1
!      16       10000     5
!      17       10001     1
!    1023  1111111111     1
!    1024 10000000000    11
!    1025 10000000001     1
!
!  Modified:
!
!    16 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be measured.
!    N should be nonnegative.
!
!    Output, integer ( kind = 4 ) I4_BIT_LO1, the position of the low 1 bit.
!
  implicit none

  integer ( kind = 4 ) bit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_bit_lo1
  integer ( kind = 4 ) n

  bit = 0
  i = n

  do

    bit = bit + 1
    i2 = i / 2

    if ( i /= 2 * i2 ) then
      exit
    end if

    i = i2

  end do

  i4_bit_lo1 = bit

  return
end
function i4_characteristic ( q )

!*****************************************************************************80
!
!! I4_CHARACTERISTIC gives the characteristic for an I4.
!
!  Discussion:
!
!    For any positive integer Q, the characteristic is:
!
!    Q, if Q is a prime;
!    P, if Q = P**N for some prime P and some integer N;
!    0, otherwise, that is, if Q is negative, 0, 1, or the product
!       of more than one distinct prime.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    05 December 2004
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738:
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, 1994, pages 494-495.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) Q, the value to be tested.
!
!    Output, integer ( kind = 4 ) I4_CHARACTERISTIC, the characteristic of Q.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_characteristic
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) q
  integer ( kind = 4 ) q_copy

  if ( q <= 1 ) then
    i4_characteristic = 0
    return
  end if
!
!  If Q is not prime, then there is at least one prime factor
!  of Q no greater than SQRT(Q)+1.
!
!  A faster code would only consider prime values of I,
!  but that entails storing a table of primes and limiting the
!  size of Q.  Simplicity and flexibility for now!
!
  i_max = int ( sqrt ( real ( q ) ) ) + 1
  q_copy = q

  do i = 2, i_max

    if ( mod ( q_copy, i ) == 0 ) then

      do while ( mod ( q_copy, i ) == 0 )
        q_copy = q_copy / i
      end do

      if ( q_copy == 1 ) then
        i4_characteristic = i
      else
        i4_characteristic = 0
      end if

      return

    end if

  end do
!
!  If no factor was found, then Q is prime.
!
  i4_characteristic = q

  return
end
function i4_div_rounded ( a, b )

!*****************************************************************************80
!
!! I4_DIV_ROUNDED computes the rounded result of I4 division.
!
!  Discussion:
!
!    This routine computes C = A / B, where A, B and C are integers
!    and C is the closest integer value to the exact real result.
!
!  Modified:
!
!    20 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the number to be divided,
!    and the divisor.
!
!    Output, integer ( kind = 4 ) I4_DIV_ROUNDED, the rounded result
!    of the division.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) a_abs
  integer ( kind = 4 ) b
  integer ( kind = 4 ) b_abs
  integer ( kind = 4 ) c
  integer ( kind = 4 ) c_abs
  integer ( kind = 4 ) c_s
  integer ( kind = 4 ) i4_div_rounded
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) i4_sign

  if ( a == 0 ) then

    c_abs = i4_huge ( )
    c_s = i4_sign ( b )

  else

    a_abs = abs ( a )
    b_abs = abs ( b )
    c_s = i4_sign ( a ) * i4_sign ( b )

    c_abs = a_abs / b_abs

    if ( ( 2 * c_abs + 1 ) * b_abs < 2 * a_abs ) then
      c_abs = c_abs + 1
    end if

  end if

  c = c_s * c_abs

  i4_div_rounded = c

  return
end
function i4_divp ( i, j )

!*****************************************************************************80
!
!! I4_DIVP returns the smallest multiple of J greater than or equal to I.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Examples:
!
!    I  J  I4_DIVP(I,J)
!
!    0  4    0
!    1  4    1
!    2  4    1
!    3  4    1
!    4  4    1
!    5  4    2
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be analyzed.
!
!    Input, integer ( kind = 4 ) J, the number, multiples of which will
!    be compared against I.  J may not be zero.
!
!    Output, integer ( kind = 4 ) I4_DIVP, the smallest multiple of J that
!    is greater than or equal to I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_divp
  integer ( kind = 4 ) j

  if ( j /= 0 ) then
    i4_divp = 1 + ( i - 1 ) / j
  else
    i4_divp = 0
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_DIVP - Fatal error!'
    write ( *, '(a)' ) '  The input value of J was zero!'
    stop
  end if

  return
end
function i4_even ( i )

!*****************************************************************************80
!
!! I4_EVEN returns TRUE if an I4 is even.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    27 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the integer to be tested.
!
!    Output, logical I4_EVEN, is TRUE if I is even.
!
  implicit none

  integer ( kind = 4 ) i
  logical i4_even
  integer ( kind = 4 ) :: i4_two = 2

  i4_even = ( mod ( i, i4_two ) == 0 )

  return
end
function i4_gcd ( i, j )

!*****************************************************************************80
!
!! I4_GCD finds the greatest common divisor of two I4's.
!
!  Discussion:
!
!    Note that only the absolute values of I and J are
!    considered, so that the result is always nonnegative.
!
!    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
!
!    If I and J have no common factor, I4_GCD is returned as 1.
!
!    Otherwise, using the Euclidean algorithm, I4_GCD is the
!    greatest common divisor of I and J.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    25 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, two numbers whose GCD is desired.
!
!    Output, integer ( kind = 4 ) I4_GCD, the greatest common divisor
!    of I and J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) :: i4_one = 1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r

  i4_gcd = 1
!
!  Return immediately if either I or J is zero.
!
  if ( i == 0 ) then
    i4_gcd = max ( i4_one, abs ( j ) )
    return
  else if ( j == 0 ) then
    i4_gcd = max ( i4_one, abs ( i ) )
    return
  end if
!
!  Set P to the larger of I and J, Q to the smaller.
!  This way, we can alter P and Q as we go.
!
  p = max ( abs ( i ), abs ( j ) )
  q = min ( abs ( i ), abs ( j ) )
!
!  Carry out the Euclidean algorithm.
!
  do

    r = mod ( p, q )

    if ( r == 0 ) then
      exit
    end if

    p = q
    q = r

  end do

  i4_gcd = q

  return
end
function i4_gcdb ( i, j, k )

!*****************************************************************************80
!
!! I4_GCDB finds the greatest common divisor of the form K**N of two I4's.
!
!  Discussion:
!
!    Note that if J is negative, I4_GCDB will also be negative.
!    This is because it is likely that the caller is forming
!    the fraction I/J, and so any minus sign should be
!    factored out of J.
!
!    If I and J are both zero, I4_GCDB is returned as 1.
!
!    If I is zero and J is not, I4_GCDB is returned as J,
!    and vice versa.
!
!    If I and J are nonzero, and have no common divisor of the
!    form K**N, I4_GCDB is returned as 1.
!
!    Otherwise, I4_GCDB is returned as the largest common divisor
!    of the form K**N shared by I and J.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, two numbers whose greatest common 
!    divisor K**N is desired.
!
!    Input, integer ( kind = 4 ) K, the possible divisor of I and J.
!
!    Output, integer ( kind = 4 ) I4_GCDB, the greatest common divisor of
!    the form K**N shared by I and J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) icopy
  integer ( kind = 4 ) i4_gcdb
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcopy
  integer ( kind = 4 ) k

  i4_gcdb = 1
!
!  If both I and J are zero, I4_GCDB is 1.
!
  if ( i == 0 .and. j == 0 ) then
    i4_gcdb = 1
    return
  end if
!
!  If just one of I and J is zero, I4_GCDB is the other one.
!
  if ( i == 0 ) then
    i4_gcdb = j
    return
  else if ( j == 0 ) then
    i4_gcdb = i
    return
  end if
!
!  Divide out K as long as you can.
!
  if ( 0 < j ) then
    i4_gcdb = 1
  else
    i4_gcdb = -1
  end if

  icopy = i
  jcopy = j

  do

    if ( mod ( icopy, k ) /= 0 .or. mod ( jcopy, k ) /= 0 ) then
      exit
    end if

    i4_gcdb = i4_gcdb * k
    icopy = icopy / k
    jcopy = jcopy / k

  end do

  return
end
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
!    bit pattern should be
!
!     01111111111111111111111111111111
!
!    In this case, its numerical value is 2147483647.
!
!    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
!    use I4_HUGE() and HUGE interchangeably.
!
!    However, when using the G95, the values returned by HUGE were
!    not equal to 2147483647, apparently, and were causing severe
!    and obscure errors in my random number generator, which needs to
!    add I4_HUGE to the seed whenever the seed is negative.  So I
!    am backing away from invoking HUGE, whereas I4_HUGE is under
!    my control.
!
!    Explanation: because under G95 the default integer type is 64 bits!
!    So HUGE ( 1 ) = a very very huge integer indeed, whereas
!    I4_HUGE ( ) = the same old 32 bit big value.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" I4.
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end
function i4_huge_normalizer ( )

!*****************************************************************************80
!
!! I4_HUGE_NORMALIZER returns the "normalizer" for I4_HUGE.
!
!  Discussion:
!
!    The value returned is 1 / ( I4_HUGE + 1 ).
!
!    For any I4, it should be the case that
!
!     -1 < I4 * I4_HUGE_NORMALIZER < 1.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    25 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) I4_HUGE_NORMALIZER, the "normalizer"
!    for I4_HUGE.
!
  implicit none

  real    ( kind = 8 ) i4_huge_normalizer

  i4_huge_normalizer = 4.656612873077392578125D-10

  return
end
function i4_is_power_of_2 ( n )

!*****************************************************************************80
!
!! I4_IS_POWER_OF_2 reports whether an I4 is a power of 2.
!
!  Discussion:
!
!    The powers of 2 are 1, 2, 4, 8, 16, and so on.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be tested.
!
!    Output, logical I4_IS_POWER_OF_2, is TRUE if N is a power of 2.
!
  implicit none

  logical i4_is_power_of_2
  integer ( kind = 4 ) :: i4_two = 2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_copy

  n_copy = n
  i4_is_power_of_2 = .false.

  if ( n_copy <= 0 ) then
    return
  end if

  do while ( n_copy /= 1 )

    if ( mod ( n_copy, i4_two ) == 1 ) then
      return
    end if

    n_copy = n_copy / 2

  end do

  i4_is_power_of_2 = .true.

  return
end
function i4_is_prime ( n )

!*****************************************************************************80
!
!! I4_IS_PRIME reports whether an I4 is prime.
!
!  Discussion:
!
!    A simple, unoptimized sieve of Erasthosthenes is used to
!    check whether N can be divided by any integer between 2
!    and SQRT(N).
!
!    Note that negative numbers, 0 and 1 are not considered prime.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    28 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the integer to be tested.
!
!    Output, logical I4_IS_PRIME, is TRUE if N is prime, and FALSE
!    otherwise.
!
  implicit none

  integer ( kind = 4 ) i
  logical i4_is_prime
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nhi

  if ( n <= 0 ) then
    i4_is_prime = .false.
    return
  end if

  if ( n == 1 ) then
    i4_is_prime = .false.
    return
  end if

  if ( n <= 3 ) then
    i4_is_prime = .true.
    return
  end if

  nhi = int ( sqrt ( real ( n ) ) )

  do i = 2, nhi
    if ( mod ( n, i ) == 0 ) then
      i4_is_prime = .false.
      return
    end if
  end do

  i4_is_prime = .true.

  return
end
function i4_lcm ( i, j )

!*****************************************************************************80
!
!! I4_LCM computes the least common multiple of two I4's.
!
!  Discussion:
!
!    The least common multiple may be defined as
!
!      LCM(I,J) = ABS( I * J ) / GCD(I,J)
!
!    where GCD(I,J) is the greatest common divisor of I and J.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the integers whose I4_LCM is desired.
!
!    Output, integer ( kind = 4 ) I4_LCM, the least common multiple of I and J.
!    I4_LCM is never negative.  I4_LCM is 0 if either I or J is zero.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) j
  integer ( kind = 4 ) i4_lcm

  i4_lcm = abs ( i * ( j / i4_gcd ( i, j ) ) )

  return
end
function i4_log_10 ( i )

!*****************************************************************************80
!
!! I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I  I4_LOG_10
!    -----  --------
!        0    0
!        1    0
!        2    0
!        9    0
!       10    1
!       11    1
!       99    1
!      100    2
!      101    2
!      999    2
!     1000    3
!     1001    3
!     9999    3
!    10000    4
!
!  Modified:
!
!    08 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 10
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_10, the integer part of the 
!    logarithm base 10 of the absolute value of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) ten_pow

  if ( i == 0 ) then

    i4_log_10 = 0

  else

    i4_log_10 = 0
    ten_pow = 10

    i_abs = abs ( i )

    do while ( ten_pow <= i_abs )
      i4_log_10 = i4_log_10 + 1
      ten_pow = ten_pow * 10
    end do

  end if

  return
end
function i4_log_2 ( i )

!*****************************************************************************80
!
!! I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
!
!  Discussion:
!
!    For positive I4_LOG_2(I), it should be true that
!      2**I4_LOG_2(X) <= |I| < 2**(I4_LOG_2(I)+1).
!    The special case of I4_LOG_2(0) returns -HUGE().
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!     I  I4_LOG_2
!
!     0  -1
!     1,  0
!     2,  1
!     3,  1
!     4,  2
!     5,  2
!     6,  2
!     7,  2
!     8,  3
!     9,  3
!    10,  3
!   127,  6
!   128,  7
!   129,  7
!
!  Modified:
!
!    13 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 2 
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_2, the integer part of the 
!    logarithm base 2 of the absolute value of I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) i4_huge

  if ( i == 0 ) then

    i4_log_2 = - i4_huge ( )

  else

    i4_log_2 = 0

    i_abs = abs ( i )

    do while ( 2 <= i_abs )
      i_abs = i_abs / 2
      i4_log_2 = i4_log_2 + 1
    end do

  end if

  return
end
function i4_log_i4 ( i4, j4 )

!*****************************************************************************80
!
!! I4_LOG_I4 returns the logarithm of an I4 to an I4 base.
!
!  Discussion:
!
!    Only the integer part of the logarithm is returned.
!
!    If 
!
!      K4 = I4_LOG_J4 ( I4, J4 ),
!
!    then we ordinarily have
!
!      J4^(K4-1) < I4 <= J4^K4.
!
!    The base J4 should be positive, and at least 2.  If J4 is negative,
!    a computation is made using the absolute value of J4.  If J4 is
!    -1, 0, or 1, the logarithm is returned as 0.
!
!    The number I4 should be positive and at least 2.  If I4 is negative,
!    a computation is made using the absolute value of I4.  If I4 is
!    -1, 0, or 1, then the logarithm is returned as 0.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!    I4  J4  K4
!
!     0   3   0
!     1   3   0
!     2   3   0
!     3   3   1
!     4   3   1
!     8   3   1
!     9   3   2
!    10   3   2
!
!  Modified:
!
!    09 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I4, the number whose logarithm is desired.
!
!    Input, integer ( kind = 4 ) J4, the base of the logarithms.
!
!    Output, integer ( kind = 4 ) I4_LOG_I4, the integer part of the logarithm
!    base abs(J4) of abs(I4).
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_abs
  integer ( kind = 4 ) i4_log_i4
  integer ( kind = 4 ) j4
  integer ( kind = 4 ) j4_abs
  integer ( kind = 4 ) value

  value = 0

  i4_abs = abs ( i4 )

  if ( 2 <= i4_abs ) then

    j4_abs = abs ( j4 )

    if ( 2 <= j4_abs ) then

      do while ( j4_abs <= i4_abs )
        i4_abs = i4_abs / j4_abs
        value = value + 1
      end do

    end if

  end if

  i4_log_i4 = value

  return
end
function i4_log_r8 ( x, b )

!*****************************************************************************80
!
!! I4_LOG_R8 returns the logarithm of an I4 to an R8 base.
!
!  Discussion:
!
!    The base B should be positive, but in any case only the absolute
!    value of B is considered.
!
!    The number X whose logarithm is desired should be positive, but
!    in any case only the absolute value of X is considered.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    An R8 is a real ( kind = 8 ) value.
!
!  Example:
!
!    If B is greater than 1, and X is positive:
!
!    if 1/B**2 <  X <= 1/B   I4_LOG_R8(X) = -1,
!    if 1/B    <  X <= 1     I4_LOG_R8(X) = 0,
!    if 1      <= X <  B,    I4_LOG_R8(X) = 0,
!    if B      <= X <  B**2  I4_LOG_R8(X) = 1,
!    if B**2   <= X <  B**3  I4_LOG_R8(X) = 2.
!
!    For positive I4_LOG_R8(X), it should be true that
!
!      ABS(B)**I4_LOG_R8(X) <= ABS(X) < ABS(B)**(I4_LOG_R8(X)+1).
!
!  Modified:
!
!    09 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the number whose logarithm base B is 
!    desired.  If X is 0, then I4_LOG_B is returned as -I4_HUGE().
!
!    Input, real ( kind = 8 ) B, the absolute value of the base of the
!    logarithms.  B must not be -1, 0, or 1.
!
!    Output, integer ( kind = 4 ) I4_LOG_R8, the integer part of the logarithm
!    base abs(B) of abs(X).
!
  implicit none

  real    ( kind = 8 ) b
  real    ( kind = 8 ) b_abs
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) i4_log_r8
  integer ( kind = 4 ) value_sign
  integer ( kind = 4 ) x
  real    ( kind = 8 ) x_abs

  if ( x == 0 ) then
    i4_log_r8 = - i4_huge ( )
    return
  end if

  b_abs = abs ( b )
  i4_log_r8 = 0

  if ( b_abs == 1.0D+00 ) then
    return
  end if

  if ( b == 0.0D+00 ) then
    return
  end if

  x_abs = abs ( real ( x ) )

  if ( b_abs < 1.0D+00 ) then
    value_sign = -1
    b_abs = 1.0D+00 / b_abs
  else
    value_sign = +1
  end if

  if ( 1.0D+00 <= x_abs .and. x_abs < b_abs ) then
    i4_log_r8 = value_sign * i4_log_r8
    return
  end if

  do while ( b_abs < x_abs )
    x_abs = x_abs / b_abs
    i4_log_r8 = i4_log_r8 + 1
  end do

  do while ( x_abs * b_abs <= 1.0D+00 )
    x_abs = x_abs * b_abs
    i4_log_r8 = i4_log_r8 - 1
  end do
!
!  If the absolute value of the base was less than 1, we inverted
!  earlier.  Now negate the logarithm to account for that.
!
  i4_log_r8 = value_sign * i4_log_r8

  return
end
subroutine i4_mant ( x, s, j, k, l )

!*****************************************************************************80
!
!! I4_MANT computes the "mantissa" of a double precision number.
!
!  Discussion:
!
!    I4_MANT computes the "mantissa" or "fraction part" of a real
!    number X, which it stores as a pair of integers, (J/K).
!
!    It also computes the sign, and the integer part of the logarithm
!    (base 2) of X.
!
!    On return:
!
!      X = S * (J/K) * 2**L
!
!    where
!
!      S is +1 or -1,
!      K is a power of 2,
!      1 <= (J/K) < 2,
!      L is an integer.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    02 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be decomposed.
!
!    Output, integer ( kind = 8 ) S, the "sign" of the number.
!    S will be -1 if X is less than 0, and +1 if X is greater
!    than or equal to zero.
!
!    Output, integer ( kind = 8 ) J, the top part of the mantissa fraction.
!
!    Output, integer ( kind = 8 ) K, the bottom part of the mantissa
!    fraction.  K is a power of 2.
!
!    Output, integer ( kind = 8 ) L, the integer part of the logarithm 
!    (base 2) of X.
!
  implicit none

  integer ( kind = 8 ) j
  integer ( kind = 8 ) k
  integer ( kind = 8 ) l
  integer ( kind = 8 ) s
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xtemp
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0D+00 ) then
    s = 1
    j = 0
    k = 1
    l = 0
    return
  end if
!
!  2: Determine the sign S.
!
  if ( 0.0D+00 < x ) then
    s = 1
    xtemp = x
  else
    s = -1
    xtemp = -x
  end if
!
!  3: Force XTEMP to lie between 1 and 2, and compute the logarithm L.
!
  l = 0

  do while ( 2.0D+00 <= xtemp )
    xtemp = xtemp / 2.0D+00
    l = l + 1
  end do

  do while ( xtemp < 1.0D+00 )
    xtemp = xtemp * 2.0D+00
    l = l - 1
  end do
!
!  4: Now strip out the mantissa as J/K.
!
  j = 0
  k = 1

  do

    j = 2 * j

    if ( 1.0D+00 <= xtemp ) then
      j = j + 1
      xtemp = xtemp - 1.0D+00
    end if

    if ( xtemp == 0.0D+00 ) then
      exit
    end if

    k = 2 * k
    xtemp = xtemp * 2.0D+00

  end do

  return
end
subroutine i4_moddiv ( n, d, m, r )

!*****************************************************************************80
!
!! I4_MODDIV breaks an I4 into a multiple of a divisor and remainder.
!
!  Discussion:
!
!    The formula used is:
!
!      N = M * D + R
!
!      0 <= || R || < || D ||
!
!    and R has the sign of N.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!    N         D       M      R
!
!   107       50      2      7
!   107      -50     -2      7
!  -107       50     -2     -7
!  -107      -50      2     -7
!
!  Modified:
!
!    01 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number to be decomposed.
!
!    Input, integer ( kind = 4 ) D, the divisor.  D may not be zero.
!
!    Output, integer ( kind = 4 ) M, the number of times N
!    is evenly divided by D.
!
!    Output, integer ( kind = 4 ) R, a remainder.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) r

  if ( d == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODDIV - Fatal error!'
    write ( *, '(a)' ) '  Input divisor D = 0'
    stop
  end if

  m = n / d
  r = n - d * m

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Examples:
!
!        I     J     MOD I4_MODP    Factorization
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
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_odd ( i )

!*****************************************************************************80
!
!! I4_ODD returns TRUE if an I4 is odd.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    27 May 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the integer to be tested.
!
!    Output, logical I4_ODD, is TRUE if I is odd.
!
  implicit none

  integer ( kind = 4 ) i
  logical i4_odd

  i4_odd = ( mod ( i + 1, 2 ) == 0 )

  return
end
function i4_power ( i, j )

!*****************************************************************************80
!
!! I4_POWER returns the integer power of an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    08 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the base and the power.  
!    J should be nonnegative.
!
!    Output, integer ( kind = 4 ) I4_POWER, the value of I^J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_power
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  if ( j < 0 ) then

    if ( i == 1 ) then
      i4_power = 1
    else if ( i == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4_POWER - Fatal error!'
      write ( *, '(a)' ) '  I^J requested, with I = 0 and J negative.'
      stop
    else
      i4_power = 0
    end if

  else if ( j == 0 ) then

    if ( i == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4_POWER - Fatal error!'
      write ( *, '(a)' ) '  I^J requested, with I = 0 and J = 0.'
      stop
    else
      i4_power = 1
    end if

  else if ( j == 1 ) then

    i4_power = i

  else

    i4_power = 1
    do k = 1, j
      i4_power = i4_power * i
    end do

  end if

  return
end
function i4_sign ( x )

!*****************************************************************************80
!
!! I4_SIGN evaluates the sign of an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the number whose sign is desired.
!
!    Output, integer ( kind = 4 ) I4_SIGN, the sign of X:
!
  implicit none

  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) x

  if ( x < 0 ) then
    i4_sign = -1
  else
    i4_sign = +1
  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
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
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
subroutine i4_swap3 ( i, j, k )

!*****************************************************************************80
!
!! I4_SWAP3 swaps three I4's.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J, K.  On output, the 
!    values of I, J, and K have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  l = i
  i = j
  j = k
  k = l

  return
end
subroutine i4_to_angle ( i, angle )

!*****************************************************************************80
!
!! I4_TO_ANGLE maps I4's to points on a circle.
!
!  Discussion:
!
!    The angles are intended to be used to select colors on a color
!    hexagon whose 6 vertices are red, yellow, green, cyan, blue,
!    magenta.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!     I   X      ANGLE
!
!     0   0/3      0
!     1   1/3    120
!     2   2/3    240
!
!     3   1/6     60
!     4   3/6    180
!     5   5/6    300
!
!     6   1/12    30
!     7   3/12    90
!     8   5/12   150
!     9   7/12   210
!    10   9/12   270
!    11  11/12   330
!
!    12   1/24    15
!    13   3/24    45
!    14   5/24    75
!    etc
!
!  Modified:
!
!    14 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the desired color.
!
!    Output, real ( kind = 8 ) ANGLE, an angle, measured in degrees, 
!    between 0 and 360.
!
  implicit none

  real    ( kind = 8 ) angle
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4

  if ( 0 <= abs ( i ) .and. abs ( i ) <= 2 ) then

    angle = 120.0D+00 * real ( abs ( i ), kind = 8 )

  else

    i1 = i4_log_2 ( abs ( i ) / 3 )
    i2 = abs ( i ) + 1 - 3 * 2**i1
    i3 = 2 * ( i2 - 1 ) + 1
    i4 = 3 * 2**( i1 + 1 )

    angle = 360.0D+00 * real ( i3, kind = 8 ) / real ( i4, kind = 8 )

  end if

  return
end
subroutine i4_to_digits_decimal ( i, n, digit )

!*****************************************************************************80
!
!! I4_TO_DIGITS_DECIMAL determines the last N decimal digits of an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
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
!    Input, integer ( kind = 4 ) I, the integer to be analyzed.
!
!    Input, integer ( kind = 4 ) N, the number of digits to determine.
!
!    Output, integer ( kind = 4 ) DIGIT(N), the last N decimal digits of I.
!    DIGIT(I) is the "coefficient" of 10**(I-1).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) digit(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_copy
  integer ( kind = 4 ), parameter :: i4_ten = 10
  integer ( kind = 4 ) j

  i_copy = i

  do j = 1, n
    digit(j) = mod ( i_copy, i4_ten )
    i_copy = ( i_copy - digit(j) ) / 10
  end do

  return
end
subroutine i4_to_fac ( intval, prime_num, npower )

!*****************************************************************************80
!
!! I4_TO_FAC converts an I4 into a product of prime factors.
!
!  Discussion:
!
!    This routine will fail if the input integer is not positive,
!    or if PRIME_NUM is too small to account for the factors of the integer.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The formula is:
!
!      INTVAL = Product ( 1 <= I <= PRIME_NUM ) PRIME(I)**NPOWER(I).
!
!  Modified:
!
!    15 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INTVAL, the integer to be factored.
!
!    Input, integer ( kind = 4 ) PRIME_NUM, the number of prime factors for
!    which storage has been allocated.
!
!    Output, integer ( kind = 4 ) NPOWER(PRIME_NUM), the powers of the primes.
!
  implicit none

  integer ( kind = 4 ) prime_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) intcopy
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) npower(prime_num)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) prime

  if ( intval <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_FAC - Fatal error!'
    write ( *, '(a)' ) '  Input integer is not positive.'
    stop
  end if
!
!  Try dividing the remainder by each prime.
!
  intcopy = intval

  do i = 1, prime_num

    npower(i) = 0

    p = prime ( i )

    do while ( mod ( intcopy, p ) == 0 )
      npower(i) = npower(i) + 1
      intcopy = intcopy / p
    end do

  end do

  return
end
subroutine i4_to_halton ( dim_num, step, seed, leap, base, r )

!*****************************************************************************80
!
!! I4_TO_HALTON computes one element of a leaped Halton subsequence.
!
!  Discussion:
!
!    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
!    sequences, each generated by a particular base.
!
!    This routine selects elements of a "leaped" subsequence of the
!    Halton sequence.  The subsequence elements are indexed by a
!    quantity called STEP, which starts at 0.  The STEP-th subsequence
!    element is simply element
!
!      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
!
!    of the original Halton sequence.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, Number 1, December 1960, pages 84-90
!
!    John Halton, GB Smith,
!    Algorithm 247:
!    Radical-Inverse Quasi-Random Point Sequence,
!    Communications of the ACM,
!    Volume 7, Number 12, December 1964, pages 701-702.
!
!    Ladislav Kocis, William Whiten,
!    Computational Investigations of Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 23, Number 2, June 1997, pages 266-294.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!    1 <= DIM_NUM is required.
!
!    Input, integer ( kind = 4 ) STEP, the index of the subsequence element.
!    0 <= STEP is required.
!
!    Input, integer ( kind = 4 ) SEED(DIM_NUM), the Halton sequence index
!    corresponding to STEP = 0.
!    0 <= SEED(1:DIM_NUM) is required.
!
!    Input, integer ( kind = 4 ) LEAP(DIM_NUM), the successive jumps in the
!    Halton sequence.  1 <= LEAP(1:DIM_NUM) is required.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the Halton bases.
!    1 < BASE(1:DIM_NUM) is required.
!
!    Output, real ( kind = 8 ) R(DIM_NUM), the STEP-th element of the leaped
!    Halton subsequence.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) base(dim_num)
  real    ( kind = 8 ) base_inv
  integer ( kind = 4 ) digit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) leap(dim_num)
  real    ( kind = 8 ) r(dim_num)
  integer ( kind = 4 ) seed(dim_num)
  integer ( kind = 4 ) seed2
  integer ( kind = 4 ) step
!
!  Check the input.
!
  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop
  end if

  if ( step < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
    write ( *, '(a)' ) ' STEP < 0.'
    stop
  end if

  if ( any ( seed(1:dim_num) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
    write ( *, '(a)' ) '  Some SEED(*) < 0.'
    stop
  end if

  if ( any ( leap(1:dim_num) < 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
    write ( *, '(a)' ) '  Some LEAP < 1.'
    stop
  end if

  if ( any ( base(1:dim_num) <= 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
    write ( *, '(a)' ) '  Some BASE <= 1.'
    stop
  end if
!
!  Calculate the data.
!
  do i = 1, dim_num

    seed2 = seed(i) + step * leap(i)

    r(i) = 0.0D+00

    base_inv = real ( 1.0D+00, kind = 8 ) / real ( base(i), kind = 8 )

    do while ( seed2 /= 0 )
      digit = mod ( seed2, base(i) )
      r(i) = r(i) + real ( digit, kind = 8 ) * base_inv
      base_inv = base_inv / real ( base(i), kind = 8 )
      seed2 = seed2 / base(i)
    end do

  end do

  return
end
function i4_to_isbn ( i )

!*****************************************************************************80
!
!! I4_TO_ISBN converts an I4 to an ISBN digit.
!
!  Discussion:
!
!    Only the integers 0 through 10 can be input.  The representation
!    of 10 is 'X'.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    16 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Book Industry Study Group,
!    The Evolution in Product Identification:
!    Sunrise 2005 and the ISBN-13,
!    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, an integer between 0 and 10.
!
!    Output, character I4_TO_ISBN, the ISBN character code of the integer.
!    If I is illegal, then I4_TO_ISBN is set to '?'.
!
  implicit none

  integer ( kind = 4 ) i
  character i4_to_isbn

       if ( i == 0 ) then
    i4_to_isbn = '0'
  else if ( i == 1 ) then
    i4_to_isbn = '1'
  else if ( i == 2 ) then
    i4_to_isbn = '2'
  else if ( i == 3 ) then
    i4_to_isbn = '3'
  else if ( i == 4 ) then
    i4_to_isbn = '4'
  else if ( i == 5 ) then
    i4_to_isbn = '5'
  else if ( i == 6 ) then
    i4_to_isbn = '6'
  else if ( i == 7 ) then
    i4_to_isbn = '7'
  else if ( i == 8 ) then
    i4_to_isbn = '8'
  else if ( i == 9 ) then
    i4_to_isbn = '9'
  else if ( i == 10 ) then
    i4_to_isbn = 'X'
  else
    i4_to_isbn = '?'
  end if

  return
end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine i4_unswap3 ( i, j, k )

!*****************************************************************************80
!
!! I4_UNSWAP3 unswaps three I4's.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J, K.  On output, the values 
!    of I, J, and K have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  l = k
  k = j
  j = i
  i = l

  return
end
function i4_walsh_1d ( x, digit )

!*****************************************************************************80
!
!! I4_WALSH_1D evaluates the Walsh function.
!
!  Discussion:
!
!    Consider the binary representation of X, and number the digits
!    in descending order, from leading to lowest, with the units digit
!    being numbered 0.
!
!    The Walsh function W(J)(X) is equal to the J-th binary digit of X.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    17 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Walsh function.
!
!    Input, integer ( kind = 4 ) DIGIT, the index of the Walsh function.
!
!    Output, integer ( kind = 4 ) I4_WALSH_1D, the value of the Walsh function.
!
  implicit none

  integer ( kind = 4 ) digit
  integer ( kind = 4 ), parameter :: i4_two = 2
  integer ( kind = 4 ) i4_walsh_1d
  integer ( kind = 4 ) n
  real    ( kind = 8 ) x
  real    ( kind = 8 ) x_copy
!
!  Hide the effect of the sign of X.
!
  x_copy = abs ( x )
!
!  If DIGIT is positive, divide by 2 DIGIT times.
!  If DIGIT is negative, multiply by 2 (-DIGIT) times.
!
  x_copy = x_copy / 2.0D+00**digit
!
!  Make it an integer.
!  Because it's positive, and we're using INT, we don't change the
!  units digit.
!
  n = int ( x_copy )
!
!  Is the units digit odd or even?
!
  if ( mod ( n, i4_two ) == 0 ) then
    i4_walsh_1d = 0
  else
    i4_walsh_1d = 1
  end if

  return
end
function i4_width ( i )

!*****************************************************************************80
!
!! I4_WIDTH returns the "width" of an I4.
!
!  Discussion:
!
!    The width of an integer is the number of characters necessary to print it.
!
!    The width of an integer can be useful when setting the appropriate output
!    format for a vector or array of values.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!        I  I4_WIDTH
!    -----  -------
!    -1234    5
!     -123    4
!      -12    3
!       -1    2
!        0    1
!        1    1
!       12    2
!      123    3
!     1234    4
!    12345    5
!
!  Modified:
!
!    04 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose width is desired.
!
!    Output, integer ( kind = 4 ) I4_WIDTH, the number of characters 
!    necessary to represent the integer in base 10, including a negative 
!    sign if necessary.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) i4_width

  if ( 0 <= i ) then
    i4_width = i4_log_10 ( i ) + 1
  else
    i4_width = i4_log_10 ( i ) + 2
  end if

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
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
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) value
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
function i4_xor ( i, j )

!*****************************************************************************80
!
!! I4_XOR calculates the exclusive OR of two I4's.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Modified:
!
!    16 February 2005
!
!  Author:
!
!   John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, two values whose exclusive OR is needed.
!
!    Output, integer ( kind = 4 ) I4_XOR, the exclusive OR of I and J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_xor
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  i1 = i
  j1 = j
  k = 0
  l = 1

  do while ( i1 /= 0 .or. j1 /= 0 )

    i2 = i1 / 2
    j2 = j1 / 2

    if ( &
      ( ( i1 == 2 * i2 ) .and. ( j1 /= 2 * j2 ) ) .or. &
      ( ( i1 /= 2 * i2 ) .and. ( j1 == 2 * j2 ) ) ) then
      k = k + l
    end if

    i1 = i2
    j1 = j2
    l = 2 * l

  end do

  i4_xor = k

  return
end
subroutine i43mat_flip_cols ( m, n, a )

!*****************************************************************************80
!
!! I43MAT_FLIP_COLS swaps the columns of an I43MAT.
!
!  Discussion:
!
!    An I43MAT is a matrix, each of whose entries is an I43, 
!    a triple of I4's.
!
!    An I43MAT can be stored as a 3 x M x N array, where M counts the "columns"
!    and N counts the "rows".
!
!  Modified:
!
!    22 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(3,M,N), the matrix whose columns 
!    are to be flipped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(3,m,n)
  integer ( kind = 4 ) b(3,m,1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n / 2
    b(1:3,1:m,    1) = a(1:3,1:m,    j)
    a(1:3,1:m,    j) = a(1:3,1:m,n+1-j)
    a(1:3,1:m,n+1-j) = b(1:3,1:m,    1)
  end do

  return
end
subroutine i43mat_flip_rows ( m, n, a )

!*****************************************************************************80
!
!! I43MAT_FLIP_ROWS swaps the rows of an I43MAT.
!
!  Discussion:
!
!    An I43MAT is a matrix, each of whose entries is an I43, 
!    a triple of I4's.
!
!    An I43MAT can be stored as a 3 x M x N array, where M counts the "columns"
!    and N counts the "rows".
!
!  Modified:
!
!    22 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(3,M,N), the matrix whose rows 
!    are to be flipped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(3,m,n)
  integer ( kind = 4 ) b(3,1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m / 2
    b(1:3,    1,1:n) = a(1:3,    i,1:n)
    a(1:3,    i,1:n) = a(1:3,m+1-i,1:n)
    a(1:3,m+1-i,1:n) = b(1:3,    1,1:n)
  end do

  return
end
subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Modified:
!
!    12 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of 
!    vectors of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Check.
!
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ICOL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  Column index I = ', i, ' is less than 1.'
    stop
  end if

  if ( n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ICOL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  N = ', n, ' is less than column index I = ', i
    stop
  end if

  if ( j < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ICOL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  Column index J = ', j, ' is less than 1.'
    stop
  end if

  if ( n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ICOL_COMPARE - Fatal error!'
    write ( *, '(a,i8,a)' ) '  N = ', n, ' is less than column index J = ', j
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4col_find ( m, n, a, ivec, col )

!*****************************************************************************80
!
!! I4COL_FIND searches an I4COL for a particular column value.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!  Example:
!
!    M = 3, N = 4,
!
!    A = (
!      1  2  3  4
!      5  6  7  8
!      9 10 11 12 )
!
!    IVEC = ( 3, 7, 11 )
!
!    COL = 3
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the table.  M is also the length of IVEC.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer ( kind = 4 ) IVEC(M), a vector to be matched with the data
!    in the array.
!
!    Output, integer ( kind = 4 ) COL, the index of the first column of 
!    the table which exactly matches every entry of IVEC, or -1 if no match
!    could be found.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) col
  integer ( kind = 4 ) ivec(m)
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    col = -1
    return
  end if

  do j = 1, n

    i = 1

    do while ( ivec(i) == a(i,j) )

      if ( i == m ) then
        col = j
        return
      end if

      i = i + 1

    end do

  end do

  col = -1

  return
end
subroutine i4col_find_item ( m, n, a, item, row, col )

!*****************************************************************************80
!
!! I4COL_FIND_ITEM searches an I4COL for a given scalar value.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the table.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer ( kind = 4 ) ITEM, the value to search for.
!
!    Output, integer ( kind = 4 ) ROW, COL, the row and column indices
!    of the first occurrence of the value ITEM.  The search
!    is conducted by columns.  If the item is not found, then
!    ROW = COL = -1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) item
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row

  do j = 1, n
    do i = 1, m
      if ( a(i,j) == item ) then
        row = i
        col = j
        return
      end if
    end do
  end do

  row = -1
  col = -1

  return
end
subroutine i4col_find_pair_wrap ( m, n, a, item1, item2, row, col )

!*****************************************************************************80
!
!! I4COL_FIND_PAIR_WRAP searches an I4COL for a pair of items.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!    The items (ITEM1, ITEM2) must occur consecutively.
!    However, wrapping is allowed, that is, if ITEM1 occurs
!    in the last row, and ITEM2 "follows" it in the first row
!    of the same column, a match is declared. 
!
!    If the pair of items is not found, then ROW = COL = -1.  
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the array.
!
!    Input, integer ( kind = 4 ) A(M,N), the array to search.
!
!    Input, integer ( kind = 4 ) ITEM1, ITEM2, the values to search for.
!
!    Output, integer ( kind = 4 ) ROW, COL, the row and column indices
!    of the first occurrence of the value ITEM1 followed immediately
!    by ITEM2.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) item1
  integer ( kind = 4 ) item2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row

  do j = 1, n
    do i = 1, m

      if ( a(i,j) == item1 ) then

        i2 = i + 1

        if ( m < i2 ) then
          i2 = 1
        end if

        if ( a(i2,j) == item2 ) then
          row = i
          col = j
          return
        end if

      end if

    end do
  end do

  row = -1
  col = -1

  return
end
subroutine i4col_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
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
    if ( 0 < indx ) then

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4col_sort_d ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_D descending sorts an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in descending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
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
    if ( 0 < indx ) then

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )
      isgn = -isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4col_sort2_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT2_A ascending sorts the elements of each column of an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!  Modified:
!
!    09 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A, and the length
!    of a vector of data.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M vectors.
!    On output, the elements of each column of A have been sorted in ascending
!    order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
!
!  Initialize.
!
  do col = 1, n

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
      if ( 0 < indx ) then

        call i4_swap ( a(i,col), a(j,col) )
!
!  Compare the I and J objects.
!
      else if ( indx < 0 ) then

        if ( a(j,col) < a(i,col) ) then
          isgn = +1
        else
          isgn = -1
        end if

      else if ( indx == 0 ) then

        exit

      end if

    end do

  end do

  return
end
subroutine i4col_sort2_d ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT2_D descending sorts elements of each column of an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!  Modified:
!
!    09 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A, and the length
!    of a vector of data.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M vectors.
!    On output, the elements of each column of A have been sorted in descending
!    order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
!
!  Initialize.
!
  do col = 1, n

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
      if ( 0 < indx ) then

        call i4_swap ( a(i,col), a(j,col) )
!
!  Compare the I and J objects.
!
      else if ( indx < 0 ) then

        if ( a(i,col) < a(j,col) ) then
          isgn = +1
        else
          isgn = -1
        end if

      else if ( indx == 0 ) then

        exit

      end if

    end do

  end do

  return
end
subroutine i4col_sorted_singleton_count ( m, n, a, singleton_num )

!*****************************************************************************80
!
!! I4COL_SORTED_SINGLETON_COUNT counts singletons in an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!    A "singleton" is an item that occurs exactly once.
!
!  Modified:
!
!    26 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Output, integer ( kind = 4 ) SINGLETON_NUM, the number of singletons.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  logical ( kind = 4 ) differ_from_next
  logical ( kind = 4 ) differ_from_previous
  integer ( kind = 4 ) j
  integer ( kind = 4 ) singleton_num

  singleton_num = 0

  if ( n <= 0 ) then
    return
  end if

  differ_from_next = .true.

  do j = 1, n

    differ_from_previous = differ_from_next

    if ( j < n ) then
      differ_from_next = any ( a(1:m,j) /= a(1:m,j+1) )
    else
      differ_from_next = .true.
    end if

    if ( differ_from_previous .and. differ_from_next ) then
      singleton_num = singleton_num + 1
    end if

  end do

  return
end
subroutine i4col_sorted_unique ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4COL_SORTED_UNIQUE keeps unique elements in a sorted I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!    The array can be sorted into ascending or descending order.
!    The important point is that identical elements must be stored
!    in adjacent positions.
!
!  Modified:
!
!    19 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the sorted array of N columns of M-vectors.
!    On output, a sorted array of columns of M-vectors.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  j1 = 1

  do j2 = 2, n

    if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
      j1 = j1 + 1
      a(1:m,j1) = a(1:m,j2)
    end if

  end do

  unique_num = j1

  return
end
subroutine i4col_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  j1 = 1

  do j2 = 2, n

    if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
      unique_num = unique_num + 1
      j1 = j2
    end if

  end do

  return
end
subroutine i4col_swap ( m, n, a, j1, j2 )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns J1 and J2 of an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, J1 = 2, J2 = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the array.
!
!    Input/output, integer ( kind = 4 ) A(M,N), an array of N columns 
!    of length M.
!
!    Input, integer ( kind = 4 ) J1, J2, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col(m)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2

  if ( j1 < 1 .or. n < j1 .or. j2 < 1 .or. n < j2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  J1 or J2 is out of bounds.'
    write ( *, '(a,i8)' ) '  J1 =    ', j1
    write ( *, '(a,i8)' ) '  J2 =    ', j2
    write ( *, '(a,i8)' ) '  N =     ', n
    stop

  end if

  if ( j1 == j2 ) then
    return
  end if

  col(1:m)  = a(1:m,j1)
  a(1:m,j1) = a(1:m,j2)
  a(1:m,j2) = col(1:m)

  return
end
subroutine i4int_to_r8int ( imin, imax, i, rmin, rmax, r )

!*****************************************************************************80
!
!! I4INT_TO_R8INT maps an I4INT to an R8INT.
!
!  Formula:
!
!    R := RMIN + ( RMAX - RMIN ) * ( I - IMIN ) / ( IMAX - IMIN )
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IMIN, IMAX, the range.
!
!    Input, integer ( kind = 4 ) I, the integer to be converted.
!
!    Input, real ( kind = 8 ) RMIN, RMAX, the range.
!
!    Output, real ( kind = 8 ) R, the corresponding value in [RMIN,RMAX].
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real    ( kind = 8 ) r
  real    ( kind = 8 ) rmax
  real    ( kind = 8 ) rmin

  if ( imax == imin ) then

    r = 0.5D+00 * ( rmin + rmax )

  else

    r = ( real ( imax - i,        kind = 8 ) * rmin   &
        + real (        i - imin, kind = 8 ) * rmax ) &
        / real ( imax     - imin, kind = 8 )

  end if

  return
end
subroutine i4mat_elim ( m, n, a )

!*****************************************************************************80
!
!! I4MAT_ELIM carries out exact Gauss elimination on an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).  On input, the M by N matrix to
!    be Gauss eliminated.  On output, the Gauss-eliminated matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) amax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol(n)
  integer ( kind = 4 ) ifact
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imult
  integer ( kind = 4 ) irow(m)
  integer ( kind = 4 ) iswap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) jmult
!
!  Initialize the swap parity counter.
!
  iswap = 1
!
!  For each column JCOL...
!
  do jcol = 1, min ( m, n )
!
!  Find the maximum element in rows JCOL through M.
!
    amax = abs ( a(jcol,jcol) )
    imax = jcol

    do i = jcol + 1, m
      if ( amax < abs ( a(i,jcol) ) ) then
        amax = abs ( a(i,jcol) )
        imax = i
      end if
    end do
!
!  If the maximum entry is nonzero, then...
!
    if ( amax /= 0 ) then
!
!  If the maximum entry does not occur in row JCOL, then swap rows.
!
      if ( imax /= jcol ) then
        iswap = -iswap
        call i4vec_swap ( n, a(jcol,1:n), a(imax,1:n) )
      end if
!
!  Eliminate all nonzero entries in column JCOL, below the diagonal entry.
!
      do i = jcol + 1, m

        if ( a(i,jcol) /= 0 ) then

          jmult = a(i,jcol)
          imult = a(jcol,jcol)
          ifact = i4_gcd ( imult, jmult )
          imult = imult / ifact
          jmult = jmult / ifact

          do j = jcol, n
            a(i,j) = jmult * a(jcol,j) - imult * a(i,j)
          end do

        end if

      end do
!
!  Remove any row or column factors.
!
      call i4mat_red ( m, n, a, irow, icol )

    end if

  end do

  return
end
subroutine i4mat_flip_cols ( m, n, a )

!*****************************************************************************80
!
!! I4MAT_FLIP_COLS swaps the columns of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is an integer matrix.
!
!    To "flip" the columns of an I4MAT is to start with something like
!
!      11 12 13 14 15
!      21 22 23 24 25
!      31 32 33 34 35
!      41 42 43 44 45
!      51 52 53 54 55
!
!    and return
!
!      15 14 13 12 11
!      25 24 23 22 21
!      35 34 33 32 31
!      45 44 43 42 41
!      55 54 53 52 51
!
!  Modified:
!
!    21 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(M,N), the matrix whose columns 
!    are to be flipped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n / 2
    b(1:m      ) = a(1:m,    j)
    a(1:m,    j) = a(1:m,n+1-j)
    a(1:m,n+1-j) = b(1:m)
  end do

  return
end
subroutine i4mat_flip_rows ( m, n, a )

!*****************************************************************************80
!
!! I4MAT_FLIP_ROWS swaps the rows of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is an integer matrix.
!
!    To "flip" the rows of an I4MAT is to start with something like
!
!      11 12 13 14 15
!      21 22 23 24 25
!      31 32 33 34 35
!      41 42 43 44 45
!      51 52 53 54 55
!
!    and return
!
!      51 52 53 54 55
!      41 42 43 44 45
!      31 32 33 34 35
!      21 22 23 24 25
!      11 12 13 14 15
!
!  Modified:
!
!    21 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(M,N), the matrix whose rows 
!    are to be flipped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m / 2
    b(      1:n) = a(    i,1:n)
    a(    i,1:n) = a(m+1-i,1:n)
    a(m+1-i,1:n) = b(      1:n)
  end do

  return
end
subroutine i4mat_l1_inverse ( n, a, b )

!*****************************************************************************80
!
!! I4MAT_L1_INVERSE inverts a unit lower triangular I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!    A unit lower triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's above the main diagonal.
!
!    The inverse of an integer unit lower triangular matrix is also
!    an integer unit lower triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call i4mat_l1_inverse ( n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Modified:
!
!    18 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, integer ( kind = 4 ) A(N,N), the unit lower triangular matrix.
!
!    Output, integer ( kind = 4 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n

    do j = 1, i - 1
      b(i,j) = -dot_product ( a(i,1:i-1), b(1:i-1,j) )
    end do

    b(i,i) = 1
    b(i,i+1:n) = 0

  end do

  return
end
subroutine i4mat_max_index ( m, n, a, i_max, j_max )

!*****************************************************************************80
!
!! I4MAT_MAX_INDEX returns the location of the maximum of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, integer ( kind = 4 ) A(M,N), the M by N matrix.
!
!    Output, integer ( kind = 4 ) I_MAX, J_MAX, the indices of the 
!    maximum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max

  i_max = -1;
  j_max = -1;

  do j = 1, n
    do i = 1, m
      if ( i == 1 .and. j == 1 ) then
        i_max = i
        j_max = j
      else if ( a(i_max,j_max) < a(i,j) ) then
        i_max = i
        j_max = j
      end if
    end do
  end do

  return
end
subroutine i4mat_min_index ( m, n, a, i_min, j_min )

!*****************************************************************************80
!
!! I4MAT_MIN_INDEX returns the location of the minimum of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Modified:
!
!    20 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, integer ( kind = 4 ) A(M,N), the M by N matrix.
!
!    Output, integer ( kind = 4 ) I_MIN, J_MIN, the indices of the
!    minimum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_min

  i_min = -1
  j_min = -1

  do j = 1, n
    do i = 1, m
      if ( i == 1 .and. j == 1 ) then
        i_min = i
        j_min = j
      else if ( a(i,j) < a(i_min,j_min) ) then
        i_min = i
        j_min = j
      end if
    end do
  end do

  return
end
subroutine i4mat_mm ( n1, n2, n3, a, b, c )

!*****************************************************************************80
!
!! I4MAT_MM multiplies two I4MAT's.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!    In FORTRAN90, this operation is more efficiently done by the
!    command:
!
!      C(1:N1,1:N3) = MATMUL ( A(1:N1,1;N2), B(1:N2,1:N3) )
!
!  Modified:
!
!    19 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, the order of the matrices.
!
!    Input, integer ( kind = 4 ) A(N1,N2), B(N2,N3), the matrices to multiply.
!
!    Output, integer ( kind = 4 ) C(N1,N3), the product matrix C = A * B.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  integer ( kind = 4 ) a(n1,n2)
  integer ( kind = 4 ) b(n2,n3)
  integer ( kind = 4 ) c(n1,n3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  do i = 1, n1
    do j = 1, n3
      c(i,j) = 0
      do k = 1, n2
        c(i,j) = c(i,j) + a(i,k) * b(k,j)
      end do
    end do
  end do

  return
end
subroutine i4mat_perm ( n, a, p )

!*****************************************************************************80
!
!! I4MAT_PERM permutes the rows and columns of a square I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Modified:
!
!    27 July 2000
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, integer ( kind = 4 ) A(N,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) is the new 
!    number of row and column I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )

  call perm_cycle ( n, p, is, nc, 1 )

  do i = 1, n

    i1 = -p(i)

    if ( 0 < i1 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( p(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( p(j1) )

              call i4_swap ( a(i1,j1), it )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( p(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:n) = abs ( p(1:n) )

  return
end
subroutine i4mat_perm_uniform ( n, a, seed )

!*****************************************************************************80
!
!! I4MAT_PERM_UNIFORM selects a random permutation of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!    The matrix is assumed to be square.  A single permutation is
!    applied to both rows and columns.
!
!  Modified:
!
!    01 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns 
!    in the array.
!
!    Input/output, integer ( kind = 4 ) A(N,N), the array to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) seed
!
!  Permute the rows and columns together.
!
  do i = 1, n

    i2 = i4_uniform ( i, n, seed )

    call i4vec_swap ( n, a(i2,1:n), a(i,1:n) )
    call i4vec_swap ( n, a(1:n,i2), a(1:n,i) )

  end do

  return
end
subroutine i4mat_perm2 ( m, n, a, p, q )

!*****************************************************************************80
!
!! I4MAT_PERM2 permutes the rows and columns of a rectangular I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Modified:
!
!    23 April 2005
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) N, number of columns in the matrix.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer ( kind = 4 ) P(M), the row permutation.  P(I) is the 
!    new number of row I.
!
!    Input, integer ( kind = 4 ) Q(N), the column permutation.  Q(I) is the 
!    new number of column I.  Note that this routine allows you to pass a 
!    single array as both P and Q.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(n)

  call perm_check ( m, p, ierror )

  call perm_cycle ( m, p, is, nc, 1 )

  call perm_check ( n, q, ierror )

  if ( 0 < q(1) ) then
    call perm_cycle ( n, q, is, nc, 1 )
  end if

  do i = 1, m

    i1 = -p(i)

    if ( 0 < i1 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( q(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( q(j1) )

              call i4_swap ( a(i1,j1), it )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                cycle
              end if

            end do

            j2 = abs ( q(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:m) = abs ( p(1:m) )

  if ( q(1) <= 0 ) then
    q(1:n) = abs ( q(1:n) )
  end if

  return
end
subroutine i4mat_perm2_uniform ( m, n, a, seed )

!*****************************************************************************80
!
!! I4MAT_PERM2_UNIFORM selects a random permutation of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!    The matrix may be rectangular.  Separate permutations are
!    applied to the rows and columns.
!
!  Modified:
!
!    01 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(M,N), the array to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) seed
!
!  Permute the rows.
!
  do i = 1, m
    i2 = i4_uniform ( i, m, seed )
    call i4vec_swap ( n, a(i2,1:n), a(i,1:n) )
  end do
!
!  Permute the columns.
!
  do j = 1, n
    j2 = i4_uniform ( j, n, seed )
    call i4vec_swap ( m, a(1:m,j2), a(1:m,j) )
  end do

  return
end
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Modified:
!
!    30 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, integer ( kind = 4 ) A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  ilo = 1
  ihi = m
  jlo = 1
  jhi = n

  call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

  return
end
subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_PRINT_SOME prints some of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Modified:
!
!    04 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 8 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8)' ) j
    end do

    write ( *, '(''  Col '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a8)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine i4mat_red ( m, n, a, row, col )

!*****************************************************************************80
!
!! I4MAT_RED divides out common factors in a row or column of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix.
!
!    Input/output, integer ( kind = 4 ) A(M,N), on input, the M by N matrix 
!    to be reduced.  On output, A has been reduced.  The greatest common 
!    factor in any row or column is 1.
!
!    Output, integer ( kind = 4 ) ROW(M), the row factors that were divided out.
!
!    Output, integer ( kind = 4 ) COL(N), the column factors that were divided
!    out.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col(n)
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row(m)

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IMAT_RED - Warning!'
    write ( *, '(a)' ) '  M must be greater than 0.'
    write ( *, '(a,i8)' ) '  Input M = ', m
    stop
  end if

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IMAT_RED - Warning!'
    write ( *, '(a)' ) '  N must be greater than 0.'
    write ( *, '(a,i8)' ) '  Input N = ', n
    stop
  end if
!
!  Remove factors common to a column.
!
  do j = 1, n
    call i4vec_red ( m, a(1:m,j), factor )
    col(j) = factor
  end do
!
!  Remove factors common to a row.
!
  do i = 1, m
    call i4vec_red ( n, a(i,1:n), factor )
    row(i) = factor
  end do

  return
end
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = * ) title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Modified:
!
!    09 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 8 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8)' ) i
    end do

    write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a8)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine i4mat_u1_inverse ( n, a, b )

!*****************************************************************************80
!
!! I4MAT_U1_INVERSE inverts a unit upper triangular I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!    A unit upper triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's below the main diagonal.
!
!    The inverse of an integer unit upper triangular matrix is also
!    an integer unit upper triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call i4mat_u1_inverse ( n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Modified:
!
!    18 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, integer ( kind = 4 ) A(N,N), the unit upper triangular matrix.
!
!    Output, integer ( kind = 4 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = n, 1, -1

    b(j+1:n,j) = 0
    b(j,j) = 1

    do i = j - 1, 1, -1
      b(i,j) = -dot_product ( a(i,i+1:j), b(i+1:j,j) )
    end do

  end do

  return
end
subroutine i4mat_uniform ( m, n, a, b, seed, x )

!*****************************************************************************80
!
!! I4MAT_UNIFORM returns a scaled pseudorandom I4MAT.
!
!  Discussion:
!
!    An I4MAT is a matrix of integer values.
!
!    The pseudorandom numbers will be scaled to be uniformly distributed
!    between A and B.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the row and column dimensions 
!    of the matrix.
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) X(M,N), a matrix of values between A and B.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x(m,n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge ( )
      end if

      r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
      r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
        +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
      value = nint ( r, kind = 4 )

      value = max ( value, min ( a, b ) )
      value = min ( value, max ( a, b ) )

      x(i,j) = value

    end do
  end do

  return
end
subroutine i4row_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4ROW_COMPARE compares two rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
!
!  Example:
!
!    Input:
!
!    M = 3, N = 4, I = 2, J = 3
!
!    A = (
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
!    14 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of M rows of vectors 
!    of length N.
!
!    Input, integer ( kind = 4 ) I, J, the rows to be compared.
!    I and J must be between 1 and M.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, row I < row J,
!     0, row I = row J,
!    +1, row J < row I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Check that I and J are legal.
!
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is less than 1.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop
  else if ( m < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
    stop
  end if

  if ( j < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is less than 1.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop
  else if ( m < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j
    write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= n )

    if ( a(i,k) < a(j,k) ) then
      isgn = -1
      return
    else if ( a(j,k) < a(i,k) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4row_find_item ( m, n, a, item, row, col )

!*****************************************************************************80
!
!! I4ROW_FIND_ITEM searches the rows of an I4ROW for a given value.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), the table to search.
!
!    Input, integer ( kind = 4 ) ITEM, the value to search for.
!
!    Output, integer ( kind = 4 ) ROW, COL, the row and column indices
!    of the first occurrence of the value ITEM.  The search
!    is conducted by rows.  If the item is not found, then
!    ROW = COL = -1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) item
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row

  row = -1
  col = -1

  do i = 1, m
    do j = 1, n
      if ( a(i,j) == item ) then
        row = i
        col = j
        return
      end if
    end do
  end do

  return
end
subroutine i4row_find_pair_wrap ( m, n, a, item1, item2, row, col )

!*****************************************************************************80
!
!! I4ROW_FIND_PAIR_WRAP searches rows of an I4ROW for a pair of items.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
!
!    The items must occur consecutively, with ITEM1 occurring
!    first.  However, wrapping is allowed.  That is, if ITEM1
!    occurs in the last column, and ITEM2 in the first, this
!    is also regarded as a match.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), the table to search.
!
!    Input, integer ( kind = 4 ) ITEM1, ITEM2, the values to search for.
!
!    Output, integer ( kind = 4 ) ROW, COL, the row and column indices
!    of the first occurrence of the value ITEM1 followed immediately
!    by ITEM2.  The search is conducted by rows.  If the pair of
!    items is not found, then ROW = COL = -1.  If COL = N,
!    the ITEM1 occurs in column N and ITEM2 occurs in column 1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) item1
  integer ( kind = 4 ) item2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) row

  row = -1
  col = -1

  do i = 1, m
    do j = 1, n

      if ( a(i,j) == item1 ) then

        if ( j < n ) then
          jp1 = j + 1
        else
          jp1 = 1
        end if

        if ( a(i,jp1) == item2 ) then
          row = i
          col = j
          return
        end if

      end if

    end do
  end do

  return
end
subroutine i4row_max ( m, n, a, amax )

!*****************************************************************************80
!
!! I4ROW_MAX returns the maximums of the rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), the array to be examined.
!
!    Output, integer ( kind = 4 ) AMAX(M), the maximums of the rows 
!    of the array.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) amax(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m

    amax(i) = a(i,1)
    do j = 2, n
      if ( amax(i) < a(i,j) ) then
        amax(i) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine i4row_mean ( m, n, a, mean )

!*****************************************************************************80
!
!! I4ROW_MEAN returns the means of the rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of data.
!
!    Output, real ( kind = 8 ) MEAN(M), the mean of each row.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) mean(m)

  do i = 1, m
    mean(i) = sum ( a(i,1:n) ) / real ( n, kind = 8 )
  end do

  return
end
subroutine i4row_min ( m, n, a, amin )

!*****************************************************************************80
!
!! I4ROW_MIN returns the minimums of the rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), the array to be examined.
!
!    Output, integer ( kind = 4 ) AMIN(M), the minimums of the rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) amin(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m

    amin(i) = a(i,1)
    do j = 2, n
      if ( a(i,j) < amin(i) ) then
        amin(i) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine i4row_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4ROW_SORT_A ascending sorts the rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
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
!      M = 5, N = 3
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
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
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
    if ( 0 < indx ) then

      call i4row_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4row_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4row_sort_d ( m, n, a )

!*****************************************************************************80
!
!! I4ROW_SORT_D descending sorts the rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows and columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in descending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
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
    if ( 0 < indx ) then

      call i4row_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4row_compare ( m, n, a, i, j, isgn )
      isgn = -isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4row_sort2_d ( m, n, a )

!*****************************************************************************80
!
!! I4ROW_SORT2_D descending sorts the elements of each row of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A, and the length
!    of a vector of data.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of M rows of N-vectors.
!    On output, the elements of each row of A have been sorted in descending
!    order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
!
!  Initialize.
!
  do row = 1, m

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
      if ( 0 < indx ) then

        call i4_swap ( a(row,i), a(row,j) )
!
!  Compare the I and J objects.
!
      else if ( indx < 0 ) then

        if ( a(row,i) < a(row,j) ) then
          isgn = +1
        else
          isgn = -1
        end if

      else if ( indx == 0 ) then

        exit

      end if

    end do

  end do

  return
end
subroutine i4row_sorted_unique ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4ROW_SORTED_UNIQUE keeps unique elements in an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
!
!    The rows of the array may be ascending or descending sorted.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(M,N), a sorted array, containing
!    M rows of data.  On output, the first UNIQUE_NUM rows
!    contain the unique rows.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  i1 = 1

  do i2 = 2, m

    if ( any ( a(i1,1:n) /= a(i2,1:n) ) ) then
      i1 = i1 + 1
      a(i1,1:n) = a(i2,1:n)
    end if

  end do

  unique_num = i1

  return
end
subroutine i4row_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4ROW_SORTED_UNIQUE_COUNT counts unique elements in an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
!
!    The rows of the array may be ascending or descending sorted.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), a sorted array, containing
!    M rows of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  i1 = 1

  do i2 = 2, m

    if ( any ( a(i1,1:n) /= a(i2,1:n) ) ) then
      unique_num = unique_num + 1
      i1 = i2
    end if

  end do

  return
end
subroutine i4row_sum ( m, n, a, rsum )

!*****************************************************************************80
!
!! I4ROW_SUM returns the sums of the rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of data.
!
!    Output, integer ( kind = 4 ) RSUM(M), the sum of the entries 
!    of each row of TABLE.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) rsum(m)

  do i = 1, m
    rsum(i) = sum ( a(i,1:n) )
  end do

  return
end
subroutine i4row_swap ( m, n, a, i1, i2 )

!*****************************************************************************80
!
!! I4ROW_SWAP swaps two rows of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(M,N), an array of data.
!
!    Input, integer ( kind = 4 ) I1, I2, the two rows to swap.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) row(n)
!
!  Check.
!
  if ( i1 < 1 .or. m < i1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I1 is out of range.'
    stop
  end if

  if ( i2 < 1 .or. m < i2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I2 is out of range.'
    stop
  end if

  if ( i1 == i2 ) then
    return
  end if

  row(1:n)  = a(i1,1:n)
  a(i1,1:n) = a(i2,1:n)
  a(i2,1:n) = row(1:n)

  return
end
subroutine i4row_variance ( m, n, a, variance )

!*****************************************************************************80
!
!! I4ROW_VARIANCE returns the variances of an I4ROW.
!
!  Discussion:
!
!    An I4ROW is an M by N array of integer values, regarded
!    as an array of M rows of length N.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), the array of data.
!
!    Output, real ( kind = 8 ) VARIANCE(M), the variance of each row.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) mean
  real    ( kind = 8 ) variance(m)

  if ( n < 2 ) then

    variance(1:m) = 0.0D+00

  else

    do i = 1, m

      mean = sum ( a(i,1:n) ) / real ( n, kind = 8 )

      variance(i) = 0.0D+00
      do j = 1, n
        variance(i) = variance(i) + ( real ( a(i,j), kind = 8 ) - mean )**2
      end do

      variance(i) = variance(i) / real ( n - 1, kind = 8 )

    end do

  end if

  return
end
subroutine i4vec_amax ( n, a, aamax )

!*****************************************************************************80
!
!! I4VEC_AMAX returns the largest magnitude in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be searched.
!
!    Output, integer ( kind = 4 ) AAMAX, the value of the entry of 
!    largest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aamax
  integer ( kind = 4 ) i

  if ( n <= 0 ) then

    aamax = 0

  else

    aamax = abs ( a(1) )

    do i = 2, n
      aamax = max ( aamax, abs ( a(i) ) )
    end do

  end if

  return
end
subroutine i4vec_amax_index ( n, a, amax_index )

!*****************************************************************************80
!
!! I4VEC_AMAX_INDEX returns the index of the largest magnitude in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be searched.
!
!    Output, integer ( kind = 4 ) AMAX_INDEX, the index of the entry 
!    of largest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aamax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) amax_index

  if ( n <= 0 ) then

    amax_index = 0

  else

    aamax = abs ( a(1) )
    amax_index = 1

    do i = 2, n

      if ( aamax < abs ( a(i) ) ) then
        aamax = abs ( a(i) )
        amax_index = i
      end if

    end do

  end if

  return
end
subroutine i4vec_amin ( n, a, aamin )

!*****************************************************************************80
!
!! I4VEC_AMIN returns the smallest magnitude in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
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
!    Input, integer ( kind = 4 ) N, the number of entries to be checked.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be checked.
!
!    Output, integer AAMIN, the value of the smallest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aamin
  integer ( kind = 4 ) i

  if ( n <= 0 ) then

    aamin = 0

  else

    aamin = abs ( a(1) )

    do i = 2, n
      aamin = min ( aamin, abs ( a(i) ) )
    end do

  end if

  return
end
subroutine i4vec_amin_index ( n, a, amin_index )

!*****************************************************************************80
!
!! I4VEC_AMIN_INDEX returns the index of the smallest magnitude in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries to be checked.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be checked.
!
!    Output, integer ( kind = 4 ) AMIN_INDEX, the entry of the smallest 
!    magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aamin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) amin_index

  if ( n <= 0 ) then

    amin_index = 0

  else

    aamin = a(1)
    amin_index = 1

    do i = 2, n

      if ( abs ( a(i) ) < aamin ) then
        aamin = abs ( a(i) )
        amin_index = i
      end if

    end do

  end if

  return
end
subroutine i4vec_aminz ( n, a, aminz )

!*****************************************************************************80
!
!! I4VEC_AMINZ returns the smallest nonzero magnitude in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
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
!    Input, integer ( kind = 4 ) N, the number of entries to be checked.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be checked.
!
!    Output, integer ( kind = 4 ) AMINZ, the value of the smallest nonzero 
!    magnitude.  If all entries are zero, AMINZ is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aminz
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iset

  aminz = 0
  iset = 0

  do i = 1, n

    if ( a(i) /= 0 ) then

      if ( iset == 0 ) then
        aminz = abs ( a(i) )
        iset = 1
      else
        aminz = min ( aminz, abs ( a(i) ) )
      end if

    end if

  end do

  return
end
subroutine i4vec_aminz_index ( n, a, aminz_index )

!*****************************************************************************80
!
!! I4VEC_AMINZ_INDEX returns the smallest nonzero magnitude in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    18 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries to be checked.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be checked.
!
!    Output, integer ( kind = 4 ) AMINZ_INDEX, the entry of the smallest
!    nonzero magnitude.  If all entries are zero, AMINZ_INDEX is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aminz
  integer ( kind = 4 ) i
  integer ( kind = 4 ) aminz_index

  aminz = 0
  aminz_index = 0

  do i = 1, n

    if ( a(i) /= 0 ) then

      if ( aminz_index == 0 .or. abs ( a(i) ) < aminz ) then
        aminz = abs ( a(i) )
        aminz_index = i
      end if

    end if

  end do

  return
end
subroutine i4vec_ascend_sub ( n, a, length, sub )

!*****************************************************************************80
!
!! I4VEC_ASCEND_SUB computes the longest ascending subsequence of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    The subsequence is required to be strictly increasing.
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
!    Input, integer ( kind = 4 ) N, the length of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be examined.
!
!    Output, integer ( kind = 4 ) LENGTH, the length of the longest 
!    increasing subsequence.
!
!    Output, integer ( kind = 4 ) SUB(N), contains in entries 1 through LENGTH
!    a longest increasing subsequence of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) length
  integer ( kind = 4 ) sub(n)
  integer ( kind = 4 ) top(n)
  integer ( kind = 4 ) top_prev(n)

  top(1:n) = 0
  top_prev(1:n) = 0
  sub(1:n) = 0

  if ( n <= 0 ) then
    length = 0
    return
  end if

  length = 0

  do i = 1, n

    k = -1

    do j = 1, length
      if ( a(i) <= a(top(j)) ) then
        k = j
        exit
      end if
    end do

    if ( k == -1 ) then
      length = length + 1
      k = length
    end if

    top(k) = i

    if ( 1 < k ) then
      top_prev(i) = top(k-1)
    else
      top_prev(i) = 0
    end if

  end do
!
!  Extract the subsequence.
!
  j = top(length)
  sub(length) = a(j)

  do i = length - 1, 1, -1
    j = top_prev(j)
    sub(i) = a(j)
  end do

  return
end
function i4vec_ascends ( n, x )

!*****************************************************************************80
!
!! I4VEC_ASCENDS determines if an I4VEC is (weakly) ascending.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Example:
!
!    X = ( -8, 1, 2, 3, 7, 7, 9 )
!
!    I4VEC_ASCENDS = TRUE
!
!    The sequence is not required to be strictly ascending.
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
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, integer ( kind = 4 ) X(N), the array to be examined.
!
!    Output, logical I4VEC_ASCENDS, is TRUE if the entries of X ascend.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical i4vec_ascends
  integer ( kind = 4 ) x(n)

  i4vec_ascends = .false.

  do i = 1, n - 1
    if ( x(i+1) < x(i) ) then
      return
    end if
  end do

  i4vec_ascends = .true.

  return
end
subroutine i4vec_axpy ( n, ia, x, incx, y, incy )

!*****************************************************************************80
!
!! I4VEC_AXPY:  Y(I) := Y(I) + A * X(I).
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    If X and Y are simple vectors, then IAXPY is equivalent to:
!
!      DO I = 1, N
!        Y(I) = Y(I) + IA * X(I)
!      END DO
!
!    However, by using the increments correctly, IAXPY can also be used
!    to manipulate rows or columns of matrices.
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of X and Y.
!
!    Input, integer ( kind = 4 ) IA, the scalar value by which each entry
!    of X is multiplied before being added to Y.
!
!    Input, integer ( kind = 4 ) X(*), the vector, a multiple of which is to be
!    added to Y.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
!    Input/output, integer ( kind = 4 ) Y(*).
!    On output, each entry of Y has been increased by
!    IA times the corresponding entry of X.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive 
!    entries of Y.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) indy
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x(*)
  integer ( kind = 4 ) y(*)

  indx = 1
  indy = 1

  do i = 1, n

    y(indy) = y(indy) + ia * x(indx)

    indx = indx + incx
    indy = indy + incy

  end do

  return
end
subroutine i4vec_bracket ( n, a, xval, left, right )

!*****************************************************************************80
!
!! I4VEC_BRACKET searches a sorted I4VEC for successive brackets of a value.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    If the values in the vector are thought of as defining intervals
!    on the number line, then this routine searches for the interval
!    containing the given value.
!
!  Modified:
!
!    25 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, integer ( kind = 4 ) A(N), an array that has been sorted 
!    into ascending order.
!
!    Input, integer ( kind = 4 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    In the most common case, 1 <= LEFT < LEFT + 1 = RIGHT <= N,
!    and A(LEFT) <= XVAL <= A(RIGHT).
!
!    Special cases:
!      Value is less than all data values:
!        LEFT = -1, RIGHT = 1, and XVAL < A(RIGHT).
!      Value is greater than all data values:
!        LEFT = N, RIGHT = -1, and A(LEFT) < XVAL.
!      Value is equal to a data value:
!        LEFT = RIGHT, and A(LEFT) = A(RIGHT) = XVAL.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) high
  integer ( kind = 4 ) left
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) right
  integer ( kind = 4 ) xval
!
!  XVAL < A(1).
!
  if ( xval < a(1) ) then
    left = -1
    right = 1
!
!  A(N) < XVAL.
!
  else if ( a(n) < xval ) then
    left = n
    right = -1
!
!  N = 1
!
  else if ( n == 1 ) then
    left = 1
    right = 1
!
!  A(1) <= XVAL <= A(N).
!
  else

    low = 1
    high = n - 1

    do

      mid = ( low + high ) / 2

      if ( high < low ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_BRACKET - Fatal error!'
        write ( *, '(a)' ) '  Algorithm or data failure.'
        stop
      end if

      if ( a(mid) == xval ) then
        left = mid
        right = mid
        exit
      else if ( a(mid+1) == xval ) then
        left = mid + 1
        right = mid + 1
        exit
      else if ( a(mid) < xval .and. xval < a(mid+1) ) then
        left = mid
        right = mid + 1
        exit
      else if ( a(mid+1) < xval ) then
        low = mid + 1
      else if ( xval < a(mid) ) then
        high = mid - 1
      end if

    end do

  end if

  return
end
subroutine i4vec_compare ( n, a1, a2, isgn )

!*****************************************************************************80
!
!! I4VEC_COMPARE compares two I4VEC's.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    Input:
!
!      A1 = ( 2, 6, 2 )
!      A2 = ( 2, 8, 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, A1 < A2,
!     0, A1 = A2,
!    +1, A2 < A1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) k

  isgn = 0

  k = 1

  do while ( k <= n )

    if ( a1(k) < a2(k) ) then
      isgn = -1
      return
    else if ( a2(k) < a1(k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4vec_copy ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC_COPY copies an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    23 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) A1(N), the vector to be copied.
!
!    Output, integer ( kind = 4 ) A2(N), a copy of A1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)

  a2(1:n) = a1(1:n)

  return
end
subroutine i4vec_cum ( n, a, a_cum )

!*****************************************************************************80
!
!! I4VEC_CUM computes the cumulutive sum of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Example:
!
!    Input:
!
!      A = (/ 1, 2, 3, 4 /)
!
!    Output:
!
!      A_CUM = (/ 0, 1, 3, 6, 10 /)
!
!  Modified:
!
!    26 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be summed.
!
!    Output, integer ( kind = 4 ) A_CUM(0:N), the cumulative sum of the 
!    entries of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) a_cum(0:n)
  integer ( kind = 4 ) i

  a_cum(0) = 0

  do i = 1, n
    a_cum(i) = a_cum(i-1) + a(i)
  end do

  return
end
function i4vec_descends ( n, x )

!*****************************************************************************80
!
!! I4VEC_DESCENDS determines if an I4VEC is decreasing.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Example:
!
!    X = ( 9, 7, 7, 3, 2, 1, -8 )
!
!    I4VEC_DESCENDS = TRUE
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, integer ( kind = 4 ) X(N), the array to be examined.
!
!    Output, logical I4VEC_DESCENDS, is TRUE if the entries of X descend.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical i4vec_descends
  integer ( kind = 4 ) x(n)

  i4vec_descends = .false.

  do i = 1, n - 1
    if ( x(i) < x(i+1) ) then
      return
    end if
  end do

  i4vec_descends = .true.

  return
end
subroutine i4vec_direct_product ( factor_index, factor_order, factor_value, &
  factor_num, point_num, x )

!*****************************************************************************80
!
!! I4VEC_DIRECT_PRODUCT creates a direct product of I4VEC's.
!
!  Discussion:
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2, 
!      ..., 
!      item JK of 1D rule K.
!
!    In particular, 
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out that task for the abscissas X.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1: 
!      Order = 4
!      X(1:4) = ( 1, 2, 3, 4 )
!
!    Rule 2:
!      Order = 3
!      X(1:3) = ( 10, 20, 30 )
!
!    Rule 3:
!      Order = 2
!      X(1:2) = ( 100, 200 )
!
!    Product Rule:
!      Order = 24
!      X(1:24) = 
!        ( 1, 10, 100 )
!        ( 2, 10, 100 )
!        ( 3, 10, 100 )
!        ( 4, 10, 100 )
!        ( 1, 20, 100 )
!        ( 2, 20, 100 )
!        ( 3, 20, 100 )
!        ( 4, 20, 100 )
!        ( 1, 30, 100 )
!        ( 2, 30, 100 )
!        ( 3, 30, 100 )
!        ( 4, 30, 100 )
!        ( 1, 10, 200 )
!        ( 2, 10, 200 )
!        ( 3, 10, 200 )
!        ( 4, 10, 200 )
!        ( 1, 20, 200 )
!        ( 2, 20, 200 )
!        ( 3, 20, 200 )
!        ( 4, 20, 200 )
!        ( 1, 30, 200 )
!        ( 2, 30, 200 )
!        ( 3, 30, 200 )
!        ( 4, 30, 200 )
!
!  Modified:
!
!    08 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, integer ( kind = 4 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the 
!    direct product.
!
!    Input/output, integer X(FACTOR_NUM,POINT_NUM), the elements of the
!    direct product, which are built up gradually.  Before the first call,
!    X might be set to 0.  After each factor has been input, X should
!    have the correct value.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of
!    values to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values 
!    to set.
!
!    Local, integer ( kind = 4 ) SKIP, the distance from the current value 
!    of START to the next location of a block of values to set.
!
!    Local, integer ( kind = 4 ) REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  integer ( kind = 4 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  integer ( kind = 4 ) x(factor_num,point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      x(factor_index,start:start+contig-1) = factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
subroutine i4vec_direct_product2 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, w )

!*****************************************************************************80
!
!! I4VEC_DIRECT_PRODUCT2 creates a direct product of I4VEC's.
!
!  Discussion:
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out the task involving the weights W.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      W(1:4) = ( 2, 3, 5, 7 )
!
!    Rule 2:
!      Order = 3
!      W(1:3) = ( 11, 13, 17 )
!
!    Rule 3:
!      Order = 2
!      W(1:2) = ( 19, 23 )
!
!    Product Rule:
!      Order = 24
!      W(1:24) =
!        ( 2 * 11 * 19 )
!        ( 3 * 11 * 19 )
!        ( 4 * 11 * 19 )
!        ( 7 * 11 * 19 )
!        ( 2 * 13 * 19 )
!        ( 3 * 13 * 19 )
!        ( 5 * 13 * 19 )
!        ( 7 * 13 * 19 )
!        ( 2 * 17 * 19 )
!        ( 3 * 17 * 19 )
!        ( 5 * 17 * 19 )
!        ( 7 * 17 * 19 )
!        ( 2 * 11 * 23 )
!        ( 3 * 11 * 23 )
!        ( 5 * 11 * 23 )
!        ( 7 * 11 * 23 )
!        ( 2 * 13 * 23 )
!        ( 3 * 13 * 23 )
!        ( 5 * 13 * 23 )
!        ( 7 * 13 * 23 )
!        ( 2 * 17 * 23 )
!        ( 3 * 17 * 23 )
!        ( 5 * 17 * 23 )
!        ( 7 * 17 * 23 )
!
!  Modified:
!
!    10 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being 
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, integer ( kind = 4 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the 
!    direct product.
!
!    Input/output, integer ( kind = 4 ) W(POINT_NUM), the elements of the
!    direct product, which are built up gradually.  Before the first call,
!    W should be set to 1.  
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of 
!    values to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values to 
!    set.
!
!    Local, integer ( kind = 4 ) SKIP, the distance from the current value 
!    of START to the next location of a block of values to set.
!
!    Local, integer ( kind = 4 ) REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  integer ( kind = 4 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  integer ( kind = 4 ) w(point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
subroutine i4vec_frac ( n, a, k, frac )

!*****************************************************************************80
!
!! I4VEC_FRAC searches for the K-th smallest element in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    Hoare's algorithm is used.
!
!  Modified:
!
!    17 July 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, integer ( kind = 4 ) A(N), array to search.  On output,
!    the elements of A have been somewhat rearranged.
!
!    Input, integer ( kind = 4 ) K, the fractile to be sought.  If K = 1, the
!    minimum entry is sought.  If K = N, the maximum is sought.
!    Other values of K search for the entry which is K-th in size.
!    K must be at least 1, and no greater than N.
!
!    Output, integer ( kind = 4 ) FRAC, the value of the K-th fractile of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) frac
  integer ( kind = 4 ) iryt
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) left

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal nonpositive value of N = ', n
    stop
  end if

  if ( k <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal nonpositive value of K = ', k
    stop
  end if

  if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal N < K, K = ', k
    stop
  end if

  left = 1
  iryt = n

  do

    if ( iryt <= left ) then
      frac = a(k)
      exit
    end if

    ix = a(k)
    i = left
    j = iryt

    do

      if ( j < i ) then

        if ( j < k ) then
          left = i
        end if

        if ( k < i ) then
          iryt = j
        end if

        exit

      end if
!
!  Find I so that IX <= A(I).
!
      do while ( a(i) < ix )
        i = i + 1
      end do
!
!  Find J so that A(J) <= IX.
!
      do while ( ix < a(j) )
        j = j - 1
      end do

      if ( i <= j ) then
        call i4_swap ( a(i), a(j) )
        i = i + 1
        j = j - 1
      end if

    end do

  end do

  return
end
function i4vec_gcd ( n, a )

!*****************************************************************************80
!
!! I4VEC_GCD finds the greatest common divisor of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
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
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be checked.
!
!    Output, integer ( kind = 4 ) I4VEC_GCD, the greatest common divisor 
!    of all entries of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) i4vec_gcd

  i4vec_gcd = maxval ( abs ( a(1:n) ) )

  do i = 1, n
    i4vec_gcd = i4_gcd ( i4vec_gcd, a(i) )
  end do

  return
end
subroutine i4vec_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_A reorders an I4VEC into an ascending heap.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    An ascending heap is an array A with the property that, for every index J,
!    A(J) <= A(2*J) and A(J) <= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
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
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the smaller of the two values,
!  and update M if necessary.
!
        if ( a(m+1) < a(m) ) then
          m = m + 1
        end if

      end if
!
!  If the small descendant is smaller than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( key <= a(m) ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_D reorders an I4VEC into an descending heap.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
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
      if ( n < m ) then
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
        if ( a(m) < a(m+1) ) then
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
subroutine i4vec_heap_d_extract ( n, a, value )

!*****************************************************************************80
!
!! I4VEC_HEAP_D_EXTRACT extracts the maximum value from a descending heap.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    In other words, the routine finds the maximum value in the
!    heap, returns that value to the user, deletes that value from
!    the heap, and restores the heap to its proper form.
!
!    This is one of three functions needed to model a priority queue.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the heap.
!
!    Input/output, integer ( kind = 4 ) A(N), the heap.
!
!    Output, integer ( kind = 4 ) VALUE, the item of maximum value, which has 
!    been removed from the heap.
!
  implicit none

  integer ( kind = 4 ) a(*)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_HEAP_D_EXTRACT - Fatal error!'
    write ( *, '(a)' ) '  The heap is empty.'
    stop
  end if
!
!  Get the maximum value.
!
  value = a(1)

  if ( n == 1 ) then
    n = 0
    return
  end if
!
!  Shift the last value down.
!
  a(1) = a(n)
!
!  Restore the heap structure.
!
  n = n - 1
  call i4vec_sort_heap_d ( n, a )

  return
end
subroutine i4vec_heap_d_insert ( n, a, value )

!*****************************************************************************80
!
!! I4VEC_HEAP_D_INSERT inserts a new I4 into a descending heap.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    This is one of three functions needed to model a priority queue.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 150.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the heap.
!
!    Input/output, integer ( kind = 4 ) A(N), the heap.
!
!    Input, integer ( kind = 4 ) VALUE, the value to be inserted.
!
  implicit none

  integer ( kind = 4 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) parent
  integer ( kind = 4 ) value

  n = n + 1
  i = n

  do while ( 1 < i )

    parent = i / 2

    if ( value <= a(parent) ) then
      exit
    end if

    a(i) = a(parent)
    i = parent

  end do

  a(i) = value

  return
end
subroutine i4vec_heap_d_max ( n, a, val_max )

!*****************************************************************************80
!
!! I4VEC_HEAP_D_MAX returns the maximum value in a descending heap of integers.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    This is one of three functions needed to model a priority queue.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 150.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the heap.
!
!    Input, integer ( kind = 4 ) A(N), the heap.
!
!    Output, integer ( kind = 4 ) VAL_MAX, the maximum value in the heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) val_max

  val_max = a(1)

  return
end
subroutine i4vec_histogram ( n, a, histo_num, histo_gram )

!*****************************************************************************80
!
!! I4VEC_HISTOGRAM computes a histogram of the elements of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    It is assumed that the entries in the vector A are nonnegative.
!    Only values between 0 and HISTO_NUM will be histogrammed.
!
!  Modified:
!
!    29 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, integer ( kind = 4 ) A(N), the array to examine.  
!
!    Input, integer ( kind = 4 ) HISTO_NUM, the maximum value for which a
!    histogram entry will be computed.
!
!    Output, integer ( kind = 4 ) HISTO_GRAM(0:HISTO_NUM), contains the 
!    number of entries of A with the values of 0 through HISTO_NUM.
!
  implicit none

  integer ( kind = 4 ) histo_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) histo_gram(0:histo_num)
  integer ( kind = 4 ) i

  histo_gram(0:histo_num) = 0

  do i = 1, n

    if ( 0 <= a(i) .and. a(i) <= histo_num ) then
      histo_gram(a(i)) = histo_gram(a(i)) + 1
    end if

  end do

  return
end
function i4vec_index ( n, a, aval )

!*****************************************************************************80
!
!! I4VEC_INDEX returns the first location of a given value in an I4VEC
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    01 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be searched.
!
!    Input, integer ( kind = 4 ) AVAL, the value to be indexed.
!
!    Output, integer ( kind = 4 ) I4VEC_INDEX, the first location in A which 
!    has the value AVAL, or -1 if no such index exists.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_index

  do i = 1, n
    if ( a(i) == aval ) then
      i4vec_index = i
      return
    end if
  end do

  i4vec_index = -1

  return
end
subroutine i4vec_index_delete_all ( n, x, indx, xval )

!*****************************************************************************80
!
!! I4VEC_INDEX_DELETE_ALL deletes all occurrences of a value in an indexed sorted list.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    16 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, integer ( kind = 4 ) X(N), the list.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, the value to be sought.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) equal1
  integer ( kind = 4 ) equal2
  integer ( kind = 4 ) get
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) put
  integer ( kind = 4 ) x(*)
  integer ( kind = 4 ) xval

  if ( n < 1 ) then
    n = 0
    return
  end if

  call i4vec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal == 0 ) then
    return
  end if

  equal1 = equal

  do

    if ( equal1 <= 1 ) then
      exit
    end if

    if ( x(indx(equal1-1)) /= xval ) then
      exit
    end if

    equal1 = equal1 - 1

  end do

  equal2 = equal

  do

    if ( n <= equal2 ) then
      exit
    end if

    if ( x(indx(equal2+1)) /= xval ) then
      exit
    end if

    equal2 = equal2 + 1

  end do
!
!  Discard certain X values.
!
  put = 0

  do get = 1, n

    if ( x(get) /= xval ) then
      put = put + 1
      x(put) = x(get)
    end if

  end do

  x(put+1:n) = 0
!
!  Adjust the INDX values.
!
  do equal = equal1, equal2
    do i = 1, n
      if ( indx(equal) < indx(i) ) then
        indx(i) = indx(i) - 1
      end if
    end do
  end do
!
!  Discard certain INDX values.
!
  indx(equal1:n+equal1-equal2-1) = indx(equal2+1:n)
  indx(n+equal1-equal2:n) = 0
!
!  Adjust N.
!
  n = put

  return
end
subroutine i4vec_index_delete_dupes ( n, x, indx, n2, x2, indx2 )

!*****************************************************************************80
!
!! I4VEC_INDEX_DELETE_DUPES deletes duplicates from an indexed sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    The output quantities N2, X2, and INDX2 are computed from the
!    input quantities by sorting, and eliminating duplicates.
!
!    The output arrays should be dimensioned of size N, unless the user
!    knows in advance what the value of N2 will be.
!
!    The output arrays may be identified with the input arrays.
!
!  Modified:
!
!    15 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer  ( kind = 4 )N, the size of the input list.
!
!    Input, integer  ( kind = 4 )X(N), the list.  
!
!    Input, integer  ( kind = 4 )INDX(N), the sort index of the list.
!
!    Output, integer  ( kind = 4 )N2, the number of unique entries in X.
!
!    Output, integer  ( kind = 4 )X2(N2), a copy of the list which has
!    been sorted, and made unique.
!
!    Output, integer  ( kind = 4 )INDX2(N2), the sort index of the new list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indx2(n)
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x2(n)
  integer ( kind = 4 ) x3(n)

  i = 0
  n3 = 0

  do

    i = i + 1

    if ( n < i ) then
      exit
    end if

    if ( 1 < i ) then
      if ( x(indx(i)) == x3(n3) ) then
        cycle
      end if
    end if

    n3 = n3 + 1
    x3(n3) = x(indx(i))

  end do
!
!  Copy data into output arrays.
!
  n2 = n3
  x2(1:n2) = x3(1:n3)
  call i4vec_indicator ( n2, indx2 )

  return
end
subroutine i4vec_index_delete_one ( n, x, indx, xval, n2, x2, indx2 )

!*****************************************************************************80
!
!! I4VEC_INDEX_DELETE_ONE deletes one copy of an I4 from an indexed sorted list.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    If the value occurs in the list more than once, only one copy is deleted.
!
!  Modified:
!
!    24 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, integer ( kind = 4 ) X(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) N2, the size of the current list.
!
!    Output, integer ( kind = 4 ) X2(N2), the list.
!
!    Output, integer ( kind = 4 ) INDX2(N2), the sort index of the list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indx2(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x2(n)
  integer ( kind = 4 ) xval

  if ( n < 1 ) then
    n2 = 0
    return
  end if

  n2 = n
  indx2(1:n2) = indx(1:n2)
  x2(1:n2) = x(1:n2)

  call i4vec_index_search ( n2, x2, indx2, xval, less, equal, more )

  if ( equal /= 0 ) then
    j = indx2(equal)
    x2(j:n2-1) = x2(j+1:n2)
    indx2(equal:n2-1) = indx2(equal+1:n2)
    do i = 1, n2 - 1
      if ( j < indx2(i) ) then
        indx2(i) = indx2(i) - 1
      end if
    end do
    n2 = n2 - 1
  end if

  return
end
subroutine i4vec_index_insert ( n, x, indx, xval )

!*****************************************************************************80
!
!! I4VEC_INDEX_INSERT inserts an I4 into an indexed sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, integer ( kind = 4 ) X(N), the list.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, the value to be sought.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) x(*)
  integer ( kind = 4 ) xval

  if ( n <= 0 ) then
    n = 1
    x(1) = xval
    indx(1) = 1
    return
  end if

  call i4vec_index_search ( n, x, indx, xval, less, equal, more )

  x(n+1) = xval
  indx(n+1:more+1:-1) = indx(n:more:-1)
  indx(more) = n + 1
  n = n + 1

  return
end
subroutine i4vec_index_insert_unique ( n, x, indx, xval )

!*****************************************************************************80
!
!! I4VEC_INDEX_INSERT_UNIQUE inserts a unique I4 into an indexed sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the size of the current list.
!    If the input value XVAL does not already occur in X, then N is increased.
!
!    Input/output, integer ( kind = 4 ) X(N), the list.
!    If the input value XVAL does not already occur in X, then it is added
!    to X.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!    If the input value XVAL does not already occur in X, then INDX is updated.
!
!    Input, integer ( kind = 4 ) XVAL, the value which will be inserted into the X
!    vector if it is not there already.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) x(*)
  integer ( kind = 4 ) xval

  if ( n <= 0 ) then
    n = 1
    x(1) = xval
    indx(1) = 1
    return
  end if
!
!  Does XVAL already occur in X?
!
  call i4vec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal == 0 ) then
    x(n+1) = xval
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1
  end if

  return
end
subroutine i4vec_index_order ( n, x, indx )

!*****************************************************************************80
!
!! I4VEC_INDEX_ORDER sorts an I4VEC using an index vector.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    The index vector itself is not modified.  Therefore, the pair
!    (X,INDX) no longer represents an index sorted vector.  If this
!    relationship is to be preserved, then simply set INDX(1:N)=(1:N).
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, integer ( kind = 4 ) X(N), the list.  On output, the list
!    has been sorted.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) y(n)

  y(1:n) = x(indx(1:n))
  x(1:n) = y(1:n)

  return
end
subroutine i4vec_index_search ( n, x, indx, xval, less, equal, more )

!*****************************************************************************80
!
!! I4VEC_INDEX_SEARCH searches for an I4 in an indexed sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, integer ( kind = 4 ) X(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, the value to be sought.
!
!    Output, integer LESS, EQUAL, MORE, the indexes in INDX of the
!    entries of X that are just less than, equal to, and just greater
!    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
!    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
!    is the greatest entry of X, then MORE is N+1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) more
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) xhi
  integer ( kind = 4 ) xlo
  integer ( kind = 4 ) xmid
  integer ( kind = 4 ) xval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n
  xlo = x(indx(lo))
  xhi = x(indx(hi))

  if ( xval < xlo ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( xval == xlo ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  if ( xhi < xval ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( xval == xhi ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))

    if ( xval == xmid ) then
      equal = mid
      less = mid - 1
      more = mid + 1
      return
    else if ( xval < xmid ) then
      hi = mid
    else if ( xmid < xval ) then
      lo = mid
    end if

  end do

  return
end
subroutine i4vec_index_sort_unique ( n, x, n2, x2, indx2 )

!*****************************************************************************80
!
!! I4VEC_INDEX_SORT_UNIQUE creates a sort index for an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    16 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, integer ( kind = 4 ) X(N), the list.  
!
!    Output, integer ( kind = 4 ) N2, the number of unique elements in X.
!
!    Output, integer ( kind = 4 ) X2(N2), a list of the unique elements of X.
!
!    Output, integer ( kind = 4 ) INDX2(N2), the sort index of the list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx2(n)
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x2(n)

  n2 = 0

  do i = 1, n
    call i4vec_index_insert_unique ( n2, x2, indx2, x(i) )
  end do

  x2(n2+1:n) = -1
  indx2(n2+1:n) = -1

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    01 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_insert ( n, a, pos, value )

!*****************************************************************************80
!
!! I4VEC_INSERT inserts a value into an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    17 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the array on input.
!
!    Input/output, integer ( kind = 4 ) A(N+1), the array.  On input, A is assumed
!    to contain N entries.  On output, A actually contains N+1 entries.
!
!    Input, integer ( kind = 4 ) POS, the position to be assigned the new entry.
!    1 <= POS <= N+1.
!
!    Input, integer ( kind = 4 ) VALUE, the value to be inserted at the given position.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) value

  if ( pos < 1 .or. n + 1 < pos ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_INSERT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal insertion position = ', pos
    stop

  else

    do i = n + 1, pos + 1, -1
      a(i) = a(i-1)
    end do

    a(pos) = value

  end if

  return
end
subroutine i4vec_mask_print ( n, a, mask_num, mask, title )

!*****************************************************************************80
!
!! I4VEC_MASK_PRINT prints a masked I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    24 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MASK_NUM, the number of masked elements.
!
!    Input, integer ( kind = 4 ) MASK(MASK_NUM), the indices of the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) mask_num
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mask(mask_num)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Masked vector printout:'

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, mask_num
    write ( *, '(2x,i8,2x,i8,2x,i10)' ) i, mask(i), a(mask(i))
  end do

  return
end
subroutine i4vec_max ( n, a, amax )

!*****************************************************************************80
!
!! I4VEC_MAX computes the maximum element of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) AMAX, the value of the largest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) amax

  amax = maxval ( a(1:n) )

  return
end
subroutine i4vec_max_index ( n, a, max_index )

!*****************************************************************************80
!
!! I4VEC_MAX_INDEX computes the index of a maximum element of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    If more than one element has the maximum value, this routine returns
!    the index of the first such element.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) MAX_INDEX, the index of the largest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) amax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_index

  if ( n <= 0 ) then

    max_index = 0

  else

    amax = a(1)
    max_index = 1

    do i = 2, n

      if ( amax < a(i) ) then
        amax = a(i)
        max_index = i
      end if

    end do

  end if

  return
end
function i4vec_max_index_last ( n, x )

!*****************************************************************************80
!
!! I4VEC_MAX_INDEX_LAST returns the last maximal element location in an I4VEC
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Example:
!
!    X = ( 5, 1, 2, 5, 0, 5, 3 )
!
!    I4VEC_MAX_INDEX_LAST = 6
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, integer ( kind = 4 ) X(N), the array to be examined.
!
!    Output, integer ( kind = 4 ) I4VEC_MAX_INDEX_LAST, the index of the last element of
!    X of maximal value.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_max_index_last
  integer ( kind = 4 ) max_last
  integer ( kind = 4 ) x(n)

  i4vec_max_index_last = 0

  do i = 1, n
    if ( i == 1 ) then
      i4vec_max_index_last = 1
      max_last = x(1)
    else if ( max_last <= x(i) ) then
      i4vec_max_index_last = i
      max_last = x(i)
    end if
  end do

  return
end
subroutine i4vec_mean ( n, a, mean )

!*****************************************************************************80
!
!! I4VEC_MEAN returns the mean of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector whose mean is desired.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  real    ( kind = 8 ) mean

  mean = real ( sum ( a(1:n) ), kind = 8 ) &
       / real ( n, kind = 8 )

  return
end
subroutine i4vec_median ( n, a, median )

!*****************************************************************************80
!
!! I4VEC_MEDIAN returns the median of an unsorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    Hoare's algorithm is used.  The values of the vector are
!    rearranged by this routine.
!
!  Modified:
!
!    18 September 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to search.  On output,
!    the order of the elements of A has been somewhat changed.
!
!    Output, integer ( kind = 4 ) MEDIAN, the value of the median of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) median

  k = ( n + 1 ) / 2

  call i4vec_frac ( n, a, k, median )

  return
end
subroutine i4vec_merge_a ( na, a, nb, b, nc, c )

!*****************************************************************************80
!
!! I4VEC_MERGE_A merges two ascending sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    The elements of A and B should be sorted in ascending order.
!
!    The elements in the output array C will also be in ascending order,
!    and unique.
!
!    The output vector C may share storage with A or B.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, integer ( kind = 4 ) A(NA), the first sorted array.
!
!    Input, integer ( kind = 4 ) NB, the dimension of B.
!
!    Input, integer ( kind = 4 ) B(NB), the second sorted array.
!
!    Output, integer ( kind = 4 ) NC, the number of elements in the output array.
!    Note that C should usually be dimensioned at least NA+NB in the
!    calling routine.
!
!    Output, integer ( kind = 4 ) C(NC), the merged unique sorted array.
!
  implicit none

  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb

  integer ( kind = 4 ) a(na)
  integer ( kind = 4 ) b(nb)
  integer ( kind = 4 ) c(na+nb)
  integer ( kind = 4 ) d(na+nb)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) na2
  integer ( kind = 4 ) nb2
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) order

  na2 = na
  nb2 = nb

  ja = 0
  jb = 0
  nc = 0

  call i4vec_order_type ( na2, a, order )

  if ( order < 0 .or. 2 < order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_MERGE - Fatal error!'
    write ( *, '(a)') '  The input array A is not ascending sorted!'
    stop
  end if

  call i4vec_order_type ( nb2, b, order )

  if ( order < 0 .or. 2 < order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_MERGE - Fatal error!'
    write ( *, '(a)' ) '  The input array B is not ascending sorted!'
    stop
  end if

  do
!
!  If we've used up all the entries of A, stick the rest of B on the end.
!
    if ( na2 <= ja ) then

      do j = 1, nb2 - jb
        jb = jb + 1
        if ( nc == 0 .or. d(nc) < b(jb) ) then
          nc = nc + 1
          d(nc) = b(jb)
        end if
      end do

      c(1:nc) = d(1:nc)

      exit
!
!  If we've used up all the entries of B, stick the rest of A on the end.
!
    else if ( nb2 <= jb ) then

      do j = 1, na2 - ja
        ja = ja + 1
        if ( nc == 0 .or. d(nc) < a(ja) ) then
          nc = nc + 1
          d(nc) = a(ja)
        end if
      end do

      c(1:nc) = d(1:nc)

      exit
!
!  Otherwise, if the next entry of A is smaller, that's our candidate.
!
    else if ( a(ja+1) <= b(jb+1) ) then

      ja = ja + 1
      if ( nc == 0 .or. d(nc) < a(ja) ) then
        nc = nc + 1
        d(nc) = a(ja)
      end if
!
!  ...or if the next entry of B is the smaller, consider that.
!
    else

      jb = jb + 1
      if ( nc == 0 .or. d(nc) < b(jb) ) then
        nc = nc + 1
        d(nc) = b(jb)
      end if
    end if

  end do

  return
end
subroutine i4vec_min ( n, a, amin )

!*****************************************************************************80
!
!! I4VEC_MIN computes the minimum element of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) AMIN, the value of the smallest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) amin

  amin = minval ( a(1:n) )

  return
end
subroutine i4vec_min_index ( n, a, imin )

!*****************************************************************************80
!
!! I4VEC_MIN_INDEX computes the index of the minimum element of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) IMIN, the index of the smallest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) amin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imin

  if ( n <= 0 ) then

    imin = 0

  else

    amin = a(1)
    imin = 1

    do i = 2, n

      if ( a(i) < amin ) then
        amin = a(i)
        imin = i
      end if

    end do

  end if

  return
end
function i4vec_nonzero_count ( n, a )

!*****************************************************************************80
!
!! I4VEC_NONZERO_COUNT counts the nonzero entries in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    01 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input, integer ( kind = 4 ) A(N), an array.
!
!    Output, integer ( kind = 4 ) I4VEC_NONZERO_COUNT, the number of nonzero entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_nonzero_count

  i4vec_nonzero_count = 0

  do i = 1, n
    if ( a(i) /= 0 ) then
      i4vec_nonzero_count = i4vec_nonzero_count + 1
    end if
  end do

  return
end
subroutine i4vec_nonzero_first ( n, x, nz, indx )

!*****************************************************************************80
!
!! I4VEC_NONZERO_FIRST left-shifts all nonzeros in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    The routine preserves the ordering of the nonzero entries.  It counts
!    the nonzeros, and returns an index vector.
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) X(N), the vector to be shifted.
!
!    Output, integer ( kind = 4 ) NZ, the number of nonzero entries in the vector.
!
!    Output, integer ( kind = 4 ) INDX(N), contains the original location of each entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) x(n)

  nz = 0

  do j = 1, n
    indx(j) = j
  end do

  j = 0

  do while ( j < n )

    j = j + 1

    if ( x(j) /= 0 ) then

      nz = nz + 1

      if ( nz /= j ) then

        x(nz) = x(j)
        x(j) = 0

        k = indx(nz)
        indx(nz) = j
        indx(j) = k

      end if
    end if
  end do

  return
end
subroutine i4vec_order_type ( n, a, order )

!*****************************************************************************80
!
!! I4VEC_ORDER_TYPE determines if an I4VEC is (non)strictly ascending/descending.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the array.
!
!    Input, integer ( kind = 4 ) A(N), the array to be checked.
!
!    Output, integer ( kind = 4 ) ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
!
!  Search for the first value not equal to A(1).
!
  i = 1

  do

    i = i + 1

    if ( n < i ) then
      order = 0
      return
    end if

    if ( a(1) < a(i) ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( a(i) < a(1) ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  Now we have a "direction".  Examine subsequent entries.
!
  do while ( i < n )

    i = i + 1

    if ( order == 1 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 3
      end if

    end if

  end do

  return
end
function i4vec_pairwise_prime ( n, a )

!*****************************************************************************80
!
!! I4VEC_PAIRWISE_PRIME checks whether an I4VEC's entries are pairwise prime.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    Two positive integers I and J are pairwise prime if they have no common
!    factor greater than 1.
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values to check.
!
!    Input, integer ( kind = 4 ) A(N), the vector of integers.
!
!    Output, logical I4VEC_PAIRWISE_PRIME, is TRUE if the vector of integers
!    is pairwise prime.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  logical i4vec_pairwise_prime
  integer ( kind = 4 ) j

  i4vec_pairwise_prime = .false.

  do i = 1, n
    do j = i + 1, n
      if ( i4_gcd ( a(i), a(j) ) /= 1 ) then
        return
      end if
    end do
  end do

  i4vec_pairwise_prime = .true.

  return
end
subroutine i4vec_part ( n, nval, a )

!*****************************************************************************80
!
!! I4VEC_PART partitions an integer NVAL into N nearly equal parts.
!
!  Example:
!
!    Input:
!
!      N = 5, NVAL = 17
!
!    Output:
!
!      A = ( 4, 4, 3, 3, 3 ).
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) NVAL, the integer to be partitioned.
!    NVAL may be positive, zero, or negative.
!
!    Output, integer ( kind = 4 ) A(N), the partition of NVAL.  The entries of
!    A add up to NVAL.  The entries of A are either all equal, or
!    differ by at most 1.  The entries of A all have the same sign
!    as NVAL, and the "largest" entries occur first.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nval

  a(1:n) = 0

  if ( 0 < nval ) then

    j = 1
    do i = 1, nval
      a(j) = a(j) + 1
      j = j + 1
      if ( n < j ) then
        j = 1
      end if
    end do

  else if ( nval < 0 ) then

    j = 1
    do i = nval, -1
      a(j) = a(j) - 1
      j = j + 1
      if ( n < j ) then
        j = 1
      end if
    end do

  end if

  return
end
subroutine i4vec_part_quick_a ( n, a, l, r )

!*****************************************************************************80
!
!! I4VEC_PART_QUICK_A reorders an I4VEC as part of a quick sort.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    The routine reorders the entries of A.  Using A(1) as a key,
!    all entries of A that are less than or equal to the key will
!    precede the key which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      N = 8
!
!      A = ( 6, 7, 3, 1, 6, 8, 2, 9 )
!
!    Output:
!
!      L = 3, R = 6
!
!      A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
!            -------        -------
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input/output, integer ( kind = 4 ) A(N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer ( kind = 4 ) L, R, the indices of A that define the three segments.
!    Let KEY = the input value of A(1).  Then
!    I <= L                 A(I) < KEY;
!         L < I < R         A(I) = KEY;
!                 R <= I    KEY < A(I).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) key
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key = a(1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( key < a(l+1) ) then
      r = r - 1
      call i4_swap ( a(r), a(l+1) )
    else if ( a(l+1) == key ) then
      m = m + 1
      call i4_swap ( a(m), a(l+1) )
      l = l + 1
    else if ( a(l+1) < key ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(i) = a(i+m)
  end do
!
!  Out of bounds here, occasionally.
!
  l = l - m

  a(l+1:l+m) = key

  return
end
subroutine i4vec_permute ( n, a, p )

!*****************************************************************************80
!
!! I4VEC_PERMUTE permutes an I4VEC in place.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    This routine permutes an array of integer "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,   4,   5,   1,   3 )
!      A = (   1,   2,   3,   4,   5 )
!
!    Output:
!
!      A    = (   2,   4,   5,   1,   3 ).
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  Pmust be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) a_temp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = -p(istart)
      cycle

    else

      a_temp = a(istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = -p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4VEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(iput) = a_temp
          exit
        end if

        a(iput) = a(iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine i4vec_permute_uniform ( n, a, seed )

!*****************************************************************************80
!
!! I4VEC_PERMUTE_UNIFORM randomly permutes an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    01 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  call perm_uniform ( n, seed, p )

  call i4vec_permute ( n, a, p )

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
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
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i12)' ) i, a(i)
  end do

  return
end
subroutine i4vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! I4VEC_PRINT_SOME prints "some" of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4 values.
!
!  Modified:
!
!    16 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices to print.
!    The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,i8)' ) i, a(i)
  end do

  return
end
function i4vec_product ( n, a )

!*****************************************************************************80
!
!! I4VEC_PRODUCT returns the product of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    In FORTRAN90, this facility is offered by the built in
!    PRODUCT function:
!
!      I4VEC_PRODUCT ( N, A ) = PRODUCT ( A(1:N) )
!
!    In MATLAB, this facility is offered by the built in
!    PROD function:
!
!      I4VEC_PRODUCT ( N, A ) = PROD ( A(1:N) )
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) I4VEC_PRODUCT, the product of the entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i4vec_product

  i4vec_product = product ( a(1:n) )

  return
end
subroutine i4vec_red ( n, a, factor )

!*****************************************************************************80
!
!! I4VEC_RED divides out common factors in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    On output, the entries of A have no common factor
!    greater than 1.
!
!  Modified:
!
!    25 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) A(*), the vector to be reduced.
!
!    Output, integer ( kind = 4 ) FACTOR, the common factor that was divided out.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
!
!  Find the smallest nonzero value.
!
  factor = 0

  do i = 1, n

    if ( a(i) /= 0 ) then

      if ( factor == 0 ) then
        factor = abs ( a(i) )
      else
        factor = min ( factor, abs ( a(i) ) )
      end if

    end if

  end do

  if ( factor == 0 ) then
    return
  end if
!
!  Find the greatest common factor of the entire vector.
!
  do i = 1, n
    factor = i4_gcd ( a(i), factor )
  end do

  if ( factor == 1 ) then
    return
  end if
!
!  Divide out the common factor.
!
  do i = 1, n
    a(i) = a(i) / factor
  end do

  return
end
subroutine i4vec_reverse ( n, a )

!*****************************************************************************80
!
!! I4VEC_REVERSE reverses the elements of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    In FORTRAN90, call I4VEC_REVERSE is equivalent to:
!
!      A(1:N) = A(N:1:-1)
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n/2
    call i4_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine i4vec_rotate ( n, m, a )

!*****************************************************************************80
!
!! I4VEC_ROTATE rotates an I4VEC in place.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Example:
!
!    Input:
!
!      N = 5, M = 2
!      A = ( 1, 2, 3, 4, 5 )
!
!    Output:
!
!      A = ( 4, 5, 1, 2, 3 ).
!
!  Modified:
!
!    07 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) M, the number of positions to the right that
!    each element should be moved.  Elements that shift pass position
!    N "wrap around" to the beginning of the array.
!
!    Input/output, integer A(N), the array to be rotated.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mcopy
  integer ( kind = 4 ) nset
  integer ( kind = 4 ) temp
!
!  Force M to be positive, between 0 and N-1.
!
  mcopy = i4_modp ( m, n )

  if ( mcopy == 0 ) then
    return
  end if

  istart = 0
  nset = 0

  do

    istart = istart + 1

    if ( n < istart ) then
      exit
    end if

    temp = a(istart)
    iget = istart
!
!  Copy the new value into the vacated entry.
!
    do

      iput = iget

      iget = iget - mcopy

      if ( iget < 1 ) then
        iget = iget + n
      end if

      if ( iget == istart ) then
        exit
      end if

      a(iput) = a(iget)
      nset = nset + 1

    end do

    a(iput) = temp
    nset = nset + 1

    if ( n <= nset ) then
      exit
    end if

  end do

  return
end
subroutine i4vec_run_count ( n, a, run_count )

!*****************************************************************************80
!
!! I4VEC_RUN_COUNT counts runs of equal values in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    A run is a sequence of equal values.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be examined.
!
!    Output, integer ( kind = 4 ) RUN_COUNT, the number of runs.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) run_count
  integer ( kind = 4 ) test

  run_count = 0

  if ( n < 1 ) then
    return
  end if

  test = 0

  do i = 1, n

    if ( i == 1 .or. a(i) /= test ) then
      run_count = run_count + 1
      test = a(i)
    end if

  end do

  return
end
subroutine i4vec_search_binary_a ( n, a, b, indx )

!*****************************************************************************80
!
!! I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC for a value.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    Binary search is used.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.9,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the array to be searched.  A must
!    be sorted in ascending order.
!
!    Input, integer ( kind = 4 ) B, the value to be searched for.
!
!    Output, integer ( kind = 4 ) INDX, the result of the search.
!    -1, B does not occur in A.
!    I, A(I) = B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) high
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid

  indx = -1

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == b ) then
      indx = mid
      exit
    else if ( a(mid) < b ) then
      low = mid + 1
    else if ( b < a(mid) ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine i4vec_search_binary_d ( n, a, b, indx )

!*****************************************************************************80
!
!! I4VEC_SEARCH_BINARY_D searches a descending sorted I4VEC for a value.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    Binary search is used.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.9,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the array to be searched.  A must
!    be sorted in descending order.
!
!    Input, integer ( kind = 4 ) B, the value to be searched for.
!
!    Output, integer ( kind = 4 ) INDX, the result of the search.
!    -1, B does not occur in A.
!    I, A(I) = B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) high
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid

  indx = -1

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == b ) then
      indx = mid
      exit
    else if ( b < a(mid) ) then
      low = mid + 1
    else if ( a(mid) < b ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine i4vec_sort_bubble_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_BUBBLE_A ascending sorts an I4VEC using bubble sort.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n - 1
    do j = i + 1, n
      if ( a(j) < a(i) ) then
        call i4_swap ( a(i), a(j) )
      end if
    end do
  end do

  return
end
subroutine i4vec_sort_bubble_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_BUBBLE_D descending sorts an I4VEC using bubble sort.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    26 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n - 1
    do j = i + 1, n
      if ( a(i) < a(j) ) then
        call i4_swap ( a(i), a(j) )
      end if
    end do
  end do

  return
end
subroutine i4vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) n1

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call i4vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  call i4_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    call i4_swap ( a(1), a(n1) )

  end do

  return
end
subroutine i4vec_sort_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_D descending sorts an I4VEC using heap sort.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) n1

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into ascending heap form.
!
  call i4vec_heap_a ( n, a )
!
!  2: Sort A.
!
!  The smallest object in the heap is in A(1).
!  Move it to position A(N).
!
  call i4_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_a ( n1, a )
!
!  Take the smallest object from A(1) and move it to A(N1).
!
    call i4_swap ( a(1), a(n1) )

  end do

  return
end
subroutine i4vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(1:N)) is sorted,
!
!    or explicitly, by the call
!
!      call i4vec_permute ( n, a, indx )
!
!    after which A(1:N) is sorted.
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine i4vec_sort_heap_index_d ( n, a, indx )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(1:N)) is sorted,
!
!    or explicitly, by the call
!
!      call i4vec_permute ( n, a, indx )
!
!    after which A(1:N) is sorted.
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j+1)) < a(indx(j)) ) then
          j = j + 1
        end if
      end if

      if ( a(indx(j)) < aval ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine i4vec_sort_insert_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_INSERT_A uses an ascending insertion sort on an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.1,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, A contains data to be sorted.
!    On output, the entries of A have been sorted in ascending order.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( a(j) <= x ) then
        exit
      end if

      a(j+1) = a(j)
      j = j - 1

    end do

    a(j+1) = x

  end do

  return
end
subroutine i4vec_sort_insert_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_INSERT_D uses a descending insertion sort on an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.1,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, A contains data to be sorted.
!    On output, the entries of A have been sorted in descending order.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( x <= a(j) ) then
        exit
      end if

      a(j+1) = a(j)
      j = j - 1

    end do

    a(j+1) = x

  end do

  return
end
subroutine i4vec_sort_quick_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_QUICK_A ascending sorts an I4VEC using quick sort.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Example:
!
!    Input:
!
!      N = 7
!
!      A = (/ 6, 7, 3, 2, 9, 1, 8 /)
!
!    Output:
!
!      A = (/ 1, 2, 3, 6, 7, 8, 9 /)
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: level_max = 25
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) l_segment
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n_segment
  integer ( kind = 4 ) rsave(level_max)
  integer ( kind = 4 ) r_segment

  if ( n <= 1 ) then
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
    call i4vec_part_quick_a ( n_segment, a(base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( level_max < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i8)' ) '  Exceeding recursion maximum of ', level_max
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

        if ( 0 < n_segment ) then
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine i4vec_sort_shell_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_SHELL_A ascending sorts an I4VEC using Shell's sort.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an array to be sorted.
!    On output, the sorted array.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) asave
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) ipow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) maxpow

  if ( n <= 1 ) then
    return
  end if
!
!  Determine the smallest MAXPOW so that
!    N <= ( 3**MAXPOW - 1 ) / 2
!
  maxpow = 1

  do while ( 3**maxpow < 2 * n + 1 )
    maxpow = maxpow + 1
  end do

  if ( 1 < maxpow ) then
    maxpow = maxpow - 1
  end if
!
!  Now sort groups of size ( 3**IPOW - 1 ) / 2.
!
  do ipow = maxpow, 1, -1

    inc = ( 3**ipow - 1 ) / 2
!
!  Sort the values with indices equal to K mod INC.
!
    do k = 1, inc
!
!  Insertion sort of the items with index
!  INC+K, 2*INC+K, 3*INC+K, ...
!
      do i = inc + k, n, inc

        asave = a(i)
        ifree = i
        j = i - inc

        do

          if ( j < 1 ) then
            exit
          end if

          if ( a(j) <= asave ) then
            exit
          end if

          ifree = j
          a(j+inc) = a(j)
          j = j - inc

        end do

        a(ifree) = asave

      end do

    end do

  end do

  return
end
subroutine i4vec_sorted_unique ( n, a, unique_num )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNIQUE finds the unique elements in a sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
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
!    Input, integer ( kind = 4 ) N, the number of elements in A.
!
!    Input/output, integer ( kind = 4 ) A(N).  On input, the sorted
!    integer array.  On output, the unique elements in A.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements in A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do itest = 2, n

    if ( a(itest) /= a(unique_num) ) then
      unique_num = unique_num + 1
      a(unique_num) = a(itest)
    end if

  end do

  return
end
subroutine i4vec_sorted_unique_count ( n, a, unique_num )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNIQUE_COUNT counts the unique elements in a sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    Because the array is sorted, this algorithm is O(N).
!
!  Modified:
!
!    29 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, integer ( kind = 4 ) A(N), the sorted array to examine.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) unique_num

  if ( n < 1 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do i = 2, n

    if ( a(i-1) /= a(i) ) then
      unique_num = unique_num + 1
    end if

  end do

  return
end
subroutine i4vec_sorted_unique_hist ( n, a, maxuniq, unique_num, auniq, acount )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNIQUE_HIST histograms the unique elements of a sorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    01 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, integer ( kind = 4 ) A(N), the array to examine.  The elements of A
!    should have been sorted.
!
!    Input, integer ( kind = 4 ) MAXUNIQ, the maximum number of unique elements
!    that can be handled.  If there are more than MAXUNIQ unique
!    elements in A, the excess will be ignored.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements of A.
!
!    Output, integer ( kind = 4 ) AUNIQ(UNIQUE_NUM), the unique elements of A.
!
!    Output, integer ( kind = 4 ) ACOUNT(UNIQUE_NUM), the number of times each element
!    of AUNIQ occurs in A.
!
  implicit none

  integer ( kind = 4 ) maxuniq
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) acount(maxuniq)
  integer ( kind = 4 ) auniq(maxuniq)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) unique_num
!
!  Start taking statistics.
!
  unique_num = 0

  do i = 1, n

    if ( i == 1 ) then

      unique_num = 1
      auniq(unique_num) = a(1)
      acount(unique_num) = 1

    else if ( a(i) == auniq(unique_num) ) then

      acount(unique_num) = acount(unique_num) + 1

    else if ( unique_num < maxuniq ) then

      unique_num = unique_num + 1
      auniq(unique_num) = a(i)
      acount(unique_num) = 1

    end if

  end do

  return
end
subroutine i4vec_split ( n, a, split, split_index )

!*****************************************************************************80
!
!! I4VEC_SPLIT "splits" an unsorted I4VEC based on a splitting value.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    If the vector is already sorted, it is simpler to do a binary search
!    on the data than to call this routine.
!
!    The vector is not assumed to be sorted before input, and is not
!    sorted during processing.  If sorting is not needed, then it is
!    more efficient to use this routine.
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to split.  On output,
!    all the entries of A that are less than or equal to SPLIT
!    are in A(1:SPLIT_INDEX).
!
!    Input, integer ( kind = 4 ) SPLIT, the value used to split the vector.
!    It is not necessary that any value of A actually equal SPLIT.
!
!    Output, integer ( kind = 4 ) SPLIT_INDEX, indicates the position of the last
!    entry of the split vector that is less than or equal to SPLIT.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) split
  integer ( kind = 4 ) split_index
!
!  Partition the vector into A1, A2, A3, where
!    A1 = A(I1:J1) holds values <= SPLIT,
!    A2 = A(I2:J2) holds untested values,
!    A3 = A(I3:J3) holds values > SPLIT.
!
  i1 = 1
  j1 = 0

  i2 = 1
  j2 = n

  i3 = n + 1
  j3 = n
!
!  Pick the next item from A2, and move it into A1 or A3.
!  Adjust indices appropriately.
!
  do i = 1, n

    if ( a(i2) <= split ) then
      i2 = i2 + 1
      j1 = j1 + 1
    else
      call i4_swap ( a(i2), a(i3-1) )
      i3 = i3 - 1
      j2 = j2 - 1
    end if

  end do

  split_index = j1

  return
end
function i4vec_sum ( n, a )

!*****************************************************************************80
!
!! I4VEC_SUM returns the sum of the entries of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    In FORTRAN90, this facility is offered by the built in
!    SUM function:
!
!      I4VEC_SUM ( N, A ) = SUM ( A(1:N) )
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) I4VEC_SUM, the sum of the entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i4vec_sum

  i4vec_sum = sum ( a(1:n) )

  return
end
subroutine i4vec_swap ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC_SWAP swaps the entries of two I4VEC's.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the arrays.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the vectors to swap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) a3(n)

  a3(1:n) = a1(1:n)
  a1(1:n) = a2(1:n)
  a2(1:n) = a3(1:n)

  return
end
subroutine i4vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Example:
!
!    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
!    TITLE = 'My vector:  '
!
!    My vector:
!
!        1    2    3    4    5
!        6    7    8    9   10
!       11
!
!  Modified:
!
!    16 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = 11 ) string
  character ( len = * ) title
  integer ( kind = 4 ) title_len

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do ilo = 1, n, 5
    ihi = min ( ilo + 5 - 1, n )
    write ( *, '(5i12)' ) a(ilo:ihi)
  end do

  return
end
subroutine i4vec_uniform ( n, a, b, seed, x )

!*****************************************************************************80
!
!! I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    The pseudorandom numbers should be scaled to be uniformly distributed
!    between A and B.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) X(N), a vector of numbers between A and B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge ( )
    end if

    r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
    r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
      +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
    value = nint ( r, kind = 4 )

    value = max ( value, min ( a, b ) )
    value = min ( value, max ( a, b ) )

    x(i) = value

  end do

  return
end
subroutine i4vec_unique_count ( n, a, unique_num )

!*****************************************************************************80
!
!! I4VEC_UNIQUE_COUNT counts the unique elements in an unsorted I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    Because the array is unsorted, this algorithm is O(N^2).
!
!  Modified:
!
!    29 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, integer ( kind = 4 ) A(N), the unsorted array to examine.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) unique_num

  unique_num = 0

  do i = 1, n

    unique_num = unique_num + 1

    do j = 1, i - 1

      if ( a(i) == a(j) ) then
        unique_num = unique_num - 1
        exit
      end if

    end do

  end do

  return
end
subroutine i4vec_value_num ( n, a, value, value_num )

!*****************************************************************************80
!
!! I4VEC_VALUE_NUM counts entries equal to a given value in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) A(N), the array to be indexed.
!
!    Input, integer ( kind = 4 ) VALUE, a value to be searched for.
!
!    Input, integer ( kind = 4 ) VALUE_NUM, the number of times the value occurs.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) value
  integer ( kind = 4 ) value_num

  value_num = 0

  do i = 1, n

    if ( a(i) == value ) then
      value_num = value_num + 1
    end if

  end do

  return
end
subroutine i4vec_value_index ( n, a, value, max_index, n_index, value_index )

!*****************************************************************************80
!
!! I4VEC_VALUE_INDEX indexesentries equal to a given value in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Example:
!
!    Input:
!
!      N = 10
!      A = (  2, 3, 1, 3, 2, 4, 2, 3, 5, 3 )
!      X_VALUE = 3
!
!    Output:
!
!      N_INDEX = 4
!      VALUE_INDEX = ( 2, 4, 8, 10 ).
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
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) A(N), the array to be indexed.
!
!    Input, integer ( kind = 4 ) VALUE, a value to be searched for.
!
!    Input, integer ( kind = 4 ) MAX_INDEX, the maximum number of indices to find.
!
!    Output, integer ( kind = 4 ) N_INDEX, the number of entries equal to VALUE.
!
!    Output, integer ( kind = 4 ) VALUE_INDEX(MAX_INDEX), the indices of entries
!    equal to VALUE.
!
  implicit none

  integer ( kind = 4 ) max_index
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n_index
  integer ( kind = 4 ) value
  integer ( kind = 4 ) value_index(max_index)

  n_index = 0

  do i = 1, n

    if ( a(i) == value ) then

      if ( max_index <= n_index ) then
        return
      end if

      n_index = n_index + 1
      value_index(n_index) = i

    end if

  end do

  return
end
subroutine i4vec_variance ( n, a, variance )

!*****************************************************************************80
!
!! I4VEC_VARIANCE returns the variance of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    01 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector whose variance is desired.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) mean
  real    ( kind = 8 ) variance

  call i4vec_mean ( n, a, mean )

  variance = 0.0D+00
  do i = 1, n
    variance = variance + ( real ( a(i), kind = 8 ) - mean )**2
  end do

  if ( 1 < n ) then
    variance = variance / real ( n - 1, kind = 8 )
  else
    variance = 0.0D+00
  end if

  return
end
function i4vec_width ( n, a )

!*****************************************************************************80
!
!! I4VEC_WIDTH returns the "width" of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    The width of an integer vector is simply the maximum of the widths of
!    its entries.
!
!    The width of a single integer is the number of characters 
!    necessary to print it.
!
!    The width of an integer vector can be useful when the vector is 
!    to be printed.
!
!  Modified:
!
!    04 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector.
!
!    Output, integer ( kind = 4 ) I4VEC_WIDTH, the width of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_width
  integer ( kind = 4 ) i4vec_width

  i4vec_width = -1

  do i = 1, n
    i4vec_width = max ( i4vec_width, i4_width ( a(i) ) )
  end do

  return
end
subroutine i4vec_zero ( n, a )

!*****************************************************************************80
!
!! I4VEC_ZERO sets the entries of an I4VEC to 0.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Modified:
!
!    22 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Output, integer ( kind = 4 ) A(N), the vector, which has been set to zero.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)

  a(1:n) = 0

  return
end
subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! I4VEC2_COMPARE compares entries of an I4VEC2.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
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
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), contain the two components of each item.
!
!    Input, integer ( kind = 4 ) I, J, the items to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item J < item I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(j) < a2(i) ) then
      isgn = +1
    end if

  else if ( a1(j) < a1(i) ) then

    isgn = +1

  end if

  return
end
subroutine i4vec2_print ( n, a, b, title )

!*****************************************************************************80
!
!! I4VEC2_PRINT prints a pair of integer vectors.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!  Modified:
!
!    09 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), B(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i10,2x,i10)' ) i, a(i), b(i)
  end do

  return
end
subroutine i4vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( n <= 1 ) then
    return
  end if
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
    if ( 0 < indx ) then

      call i4_swap ( a1(i), a1(j) )
      call i4_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4vec2_sort_d ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_D descending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( n <= 1 ) then
    return
  end if
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
    if ( 0 < indx ) then

      call i4_swap ( a1(i), a1(j) )
      call i4_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )
      isgn = -isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4vec2_sorted_unique ( n, a1, a2, unique_num )

!*****************************************************************************80
!
!! I4VEC2_SORTED_UNIQUE keeps the unique elements in a sorted I4VEC2.
!
!  Discussion:
!
!    An I4VEC2 is a pair of I4VEC's.
!
!    An I4VEC is a vector of I4's.
!
!    Entry K of an I4VEC2 is the pair of values located
!    at the K-th entries of the two I4VEC's.
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
!    Input, integer ( kind = 4 ) N, the number of items.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N).
!    On input, the array of N items.
!    On output, an array of unique items.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique items.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do itest = 2, n

    if ( a1(itest) /= a1(unique_num) .or. a2(itest) /= a2(unique_num) ) then

      unique_num = unique_num + 1

      a1(unique_num) = a1(itest)
      a2(unique_num) = a2(itest)

    end if

  end do

  return
end
function i8_huge ( )

!*****************************************************************************80
!
!! I8_HUGE returns a "huge" I8.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2**63 - 1, and its
!    bit pattern should be
!
!     0111111111111111111111111111111111111111111111111111111111111111
!
!    In this case, its numerical value is 9223372036854775807.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 8 ) I8_HUGE, a "huge" I8.
!
  implicit none

  integer ( kind = 8 ) i8
  integer ( kind = 8 ) i8_huge

  i8_huge = huge ( i8 )

  return
end
function i8_huge_normalizer ( )

!*****************************************************************************80
!
!! I8_HUGE_NORMALIZER returns the "normalizer" for I8_HUGE.
!
!  Discussion:
!
!    The value returned is 1 / ( I8_HUGE + 1 ).
!
!    For any I8, it should be the case that
!
!     -1 < I8 * I8_HUGE_NORMALIZER < 1.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) I8_HUGE_NORMALIZER, the "normalizer"
!    for I8_HUGE.
!
  implicit none

  real    ( kind = 8 ) i8_huge_normalizer

  i8_huge_normalizer = 1.084202172485504434007D-19

  return
end
function i8_xor ( i, j )

!*****************************************************************************80
!
!! I8_XOR calculates the exclusive OR of two integers.
!
!  Modified:
!
!    16 February 2005
!
!  Author:
!
!   John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) I, J, two values whose exclusive OR is needed.
!
!    Output, integer ( kind = 8 ) I8_XOR, the exclusive OR of I and J.
!
  implicit none

  integer ( kind = 8 ) i
  integer ( kind = 8 ) i1
  integer ( kind = 8 ) i2
  integer ( kind = 8 ) i8_xor
  integer ( kind = 8 ) j
  integer ( kind = 8 ) j1
  integer ( kind = 8 ) j2
  integer ( kind = 8 ) k
  integer ( kind = 8 ) l

  i1 = i
  j1 = j
  k = 0
  l = 1

  do while ( i1 /= 0 .or. j1 /= 0 )

    i2 = i1 / 2
    j2 = j1 / 2

    if ( &
      ( ( i1 == 2 * i2 ) .and. ( j1 /= 2 * j2 ) ) .or. &
      ( ( i1 /= 2 * i2 ) .and. ( j1 == 2 * j2 ) ) ) then
      k = k + l
    end if

    i1 = i2
    j1 = j2
    l = 2 * l

  end do

  i8_xor = k

  return
end
subroutine i4i4_sort_a ( i1, i2, j1, j2 )

!*****************************************************************************80
!
!! I4I4_SORT_A ascending sorts a pair of integers.
!
!  Discussion:
!
!    An I4I4 is a pair of integers, regarded as a single data item.
!
!    The program allows the reasonable call:
!
!      call i4i4_sort_a ( i1, i2, i1, i2 )
! 
!    and this will return the reasonable result.
!
!  Modified:
!
!    11 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1, I2, the values to sort.
!
!    Output, integer ( kind = 4 ) J1, J2, the sorted values.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
!
!  Copy arguments, so that the user can make "reasonable" calls like:
!
!    call i4i4_sort_a ( i1, i2, i1, i2 )
!
  k1 = i1
  k2 = i2

  j1 = min ( k1, k2 )
  j2 = max ( k1, k2 )

  return
end
subroutine i4i4i4_sort_a ( i1, i2, i3, j1, j2, j3 )

!*****************************************************************************80
!
!! I4I4I4_SORT_A ascending sorts a triple of integers.
!
!  Discussion:
!
!    An I4I4I4 is a triple of integers, regarded as a single data item.
!
!    The program allows the reasonable call:
!
!      call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
! 
!    and this will return the reasonable result.
!
!  Modified:
!
!    13 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1, I2, I3, the values to sort.
!
!    Output, integer ( kind = 4 ) J1, J2, J3, the sorted values.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
!
!  Copy arguments, so that the user can make "reasonable" calls like:
!
!    call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
!
  k1 = i1
  k2 = i2
  k3 = i3

  j1 = min ( min ( k1, k2 ), min ( k2, k3 ) )
  j2 = min ( max ( k1, k2 ), &
       min ( max ( k2, k3 ), max ( k3, k1 ) ) )
  j3 = max ( max ( k1, k2 ), max ( k2, k3 ) )

  return
end
subroutine ij_next ( i, j, n )

!*****************************************************************************80
!
!! IJ_NEXT returns the next matrix index.
!
!  Discussion:
!
!    For N = 3, the sequence of indices returned is:
!
!      (1,1), (1,2), (1,3), (2,1), (2,2), (2,3), (3,1), (3,2), (3,3), (0,0).
!
!    Note that once the value (N,N) is returned, the next value returned
!    will be (0,0).
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
!    Input/output, integer ( kind = 4 ) I, J.  On input, the current pair of indices.
!    On output, the next pair of indices.  If either index is illegal on
!    input, the output value of (I,J) will be (1,1).
!
!    Input, integer ( kind = 4 ) N, the maximum value for I and J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n

  if ( n < 1 ) then
    i = 0
    j = 0
    return
  end if

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then
    i = 1
    j = 1
    return
  end if

  if ( j < n ) then
    j = j + 1
  else if ( i < n ) then
    i = i + 1
    j = 1
  else
    i = 0
    j = 0
  end if

  return
end
subroutine ij_next_gt ( i, j, n )

!*****************************************************************************80
!
!! IJ_NEXT_GT returns the next matrix index, with the constraint that I < J.
!
!  Discussion:
!
!    For N = 5, the sequence of indices returned is:
!
!      (1,2), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5), (3,4), (3,5), (4,5).
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
!    Input/output, integer ( kind = 4 ) I, J.  On input, the current pair of indices.
!    On output, the next pair of indices.  If either index is illegal on
!    input, the output value of (I,J) will be (1,2).
!
!    Input, integer ( kind = 4 ) N, the maximum value for I and J.
!    A value of N less than 2 is nonsense.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n

  if ( n < 2 ) then
    i = 0
    j = 0
    return
  end if

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j .or. j <= i ) then
    i = 1
    j = 2
    return
  end if

  if ( j < n ) then
    j = j + 1
  else if ( i < n - 1 ) then
    i = i + 1
    j = i + 1
  else
    i = 0
    j = 0
  end if

  return
end
subroutine index_box2_next_2d ( n1, n2, ic, jc, i, j, more )

!*****************************************************************************80
!
!! INDEX_BOX2_NEXT_2D produces indices on the surface of a box in 2D.
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
!    Input, integer ( kind = 4 ) N1, N2, the half-widths of the box, that is, the
!    maximum distance allowed between (IC,JC) and (I,J).
!
!    Input, integer  ( kind = 4 )IC, JC, the central cell of the box.
!
!    Input/output, integer ( kind = 4 ) I, J.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I and J are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  logical more
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

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
  if ( jc + n2 < j ) then
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

!*****************************************************************************80
!
!! INDEX_BOX2_NEXT_3D produces indices on the surface of a box in 3D.
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
!    Input, integer ( kind = 4 ) N1, N2, N3, the "half widths" of the box, that is, the
!    maximum distances from the central cell allowed for I, J and K.
!
!    Input, integer ( kind = 4 ) IC, JC, KC, the central cell of the box.
!
!    Input/output, integer ( kind = 4 ) I, J, K.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I, J, and K are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  logical more
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

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
  if ( kc + n3 < k ) then
    k = kc - n3
    j = j + 1
  else if ( k < kc + n3 .and. &
    ( i == ic - n1 .or. &
      i == ic + n1 .or. &
      j == jc - n2 .or. &
      j == jc + n2 ) ) then
    return
  else
    k = kc + n3
    return
  end if
!
!  Check J.
!
  if ( jc + n2 < j ) then
    j = jc - n2
    i = i + 1
  else if ( j < jc + n2 .and. &
    ( i == ic - n1 .or. &
      i == ic + n1 .or. &
      k == kc - n3 .or. &
      k == kc + n3 ) ) then
    return
  else
    j = jc + n2
    return
  end if

  return
end
subroutine isbn_check ( isbn, check )

!*****************************************************************************80
!
!! ISBN_CHECK checks an ISBN code.
!
!  Discussion:
!
!    ISBN stands for International Standard Book Number.  A unique ISBN
!    is assigned to each new book.  The ISBN includes 10 digits.  There is
!    an initial digit, then a dash, then a set of digits which are a
!    code for the publisher, another digit, and then the check digit:
!
!      initial-publisher-book-check
!
!    The number of digits used for the publisher and book codes can vary,
!    but the check digit is always one digit, and the total number of
!    digits is always 10.
!
!    The check digit is interesting, because it is a way of trying to
!    make sure that an ISBN has not been incorrectly copied.  Specifically,
!    if the ISBN is correct, then its ten digits will satisfy
!
!       10 * A + 9 * B + 8 * C + 7 * D + 6 * E
!      + 5 * F * 4 * G * 3 * H + 2 * I +     J  = 0 mod 11.
!
!    Here, we've taken 'A' to represent the first digit and 'J' the
!    last (which is the check digit).  In order for the code to work,
!    the value of J must be allowed to be anything from 0 to 10.  In
!    cases where J works out to be 10, the special digit 'X' is used.
!    An 'X' digit can only occur in this last check-digit position
!    on an ISBN.
!
!  Example:
!
!    0-8493-9640-9
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Book Industry Study Group,
!    The Evolution in Product Identification:
!    Sunrise 2005 and the ISBN-13,
!    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
!
!  Parameters:
!
!    Input, character ( len = * ) ISBN, an ISBN code.
!
!    Output, integer ( kind = 4 ) CHECK, the value of the ISBN check sum.
!    If CHECK is zero, the ISBN code is legitimate.
!    If CHECK is -1, then the ISBN code is not legitimate because it does
!    not contain exactly 10 digits.  If CHECK is between 1 and 10, then
!    the ISBN code has the right number of digits, but at least one of
!    the digits is incorrect.
!
  implicit none

  character c
  logical ch_is_digit
  integer ( kind = 4 ) check
  integer ( kind = 4 ) digit(10)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_eleven = 11
  character ( len = * ) isbn
  integer ( kind = 4 ) isbn_to_i4
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) num_digit
!
!  Determine how many digits have been supplied.
!
  lenc = len_trim ( isbn )

  i = 0
  num_digit = 0

  do

    i = i + 1

    if ( lenc < i ) then
      exit
    end if

    c = isbn(i:i)

    if ( ch_is_digit ( c ) ) then

      num_digit = num_digit + 1
      digit(num_digit) = isbn_to_i4 ( c )

    else if ( ( num_digit == 9 .and. isbn(i:i) == 'X' ) .or. &
              ( num_digit == 9 .and. isbn(i:i) == 'x' ) ) then

      num_digit = num_digit + 1
      digit(num_digit) = isbn_to_i4 ( c )

    end if

    if ( 10 <= num_digit ) then
      exit
    end if

  end do
!
!  If we didn't get exactly 10 digits, return with an error.
!
  if ( num_digit /= 10 ) then
    check = -1
    return
  end if
!
!  Compute the checksum.
!
  check = 0
  do i = 1, 10
    check = check + ( 11 - i ) * digit(i)
  end do

  check = mod ( check, i4_eleven )

  return
end
subroutine isbn_fill ( isbn )

!*****************************************************************************80
!
!! ISBN_FILL fills in a missing digit in an ISBN code.
!
!  Example:
!
!    Input:
!
!      0-8493-9?40-9
!
!    Output:
!
!      0-8493-9640-9
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Book Industry Study Group,
!    The Evolution in Product Identification:
!    Sunrise 2005 and the ISBN-13,
!    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
!
!  Parameters:
!
!    Input/output, character ( len = * ) ISBN, a partial ISBN code.  On input,
!    a single digit has been replaced by the character '?', signifying
!    that that digit is missing.  The routine replaces the question
!    mark by the correct digit.
!
  implicit none

  character c
  logical ch_is_digit
  integer ( kind = 4 ) check
  integer ( kind = 4 ) digit(10)
  integer ( kind = 4 ) digit_pos
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_eleven = 11
  character i4_to_isbn
  character ( len = * ) isbn
  integer ( kind = 4 ) isbn_pos
  integer ( kind = 4 ) isbn_to_i4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) num_digit

  lenc = len_trim ( isbn )

  i = 0
  isbn_pos = -1
  digit_pos = -1
  num_digit = 0

  do

    i = i + 1

    if ( lenc < i ) then
      exit
    end if

    c = isbn(i:i)

    if ( ch_is_digit ( c ) ) then

      num_digit = num_digit + 1
      digit(num_digit) = isbn_to_i4 ( c )

    else if ( ( num_digit == 9 .and. isbn(i:i) == 'X' ) .or. &
              ( num_digit == 9 .and. isbn(i:i) == 'x' ) ) then

      num_digit = num_digit + 1
      digit(num_digit) = isbn_to_i4 ( c )

    else if ( c == '?' ) then

      if ( isbn_pos == -1 ) then

        num_digit = num_digit + 1
        digit(num_digit) = 0
        digit_pos = num_digit
        isbn_pos = i

      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ISBN_FILL - Fatal error!'
        write ( *, '(a)' ) '  Only one question mark is allowed!'
        return
      end if

    end if

    if ( 10 <= num_digit ) then
      exit
    end if

  end do

  if ( num_digit /= 10 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ISBN_FILL - Fatal error!'
    write ( *, '(a)' ) '  The input ISBN code did not have 10 digits.'
    return
  end if

  if ( isbn_pos == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ISBN_FILL - Fatal error!'
    write ( *, '(a)' ) '  A question mark is required!'
    return
  end if

  check = 0
  do i = 1, 10
    check = check + ( i4_eleven - i ) * digit(i)
  end do

  check = mod ( check, i4_eleven )

  if ( check == 0 ) then

    k = 0
!
!  Need to solve the modular equation:
!
!    A * X = B mod C
!
!  Below is a stupid way.  One day I will come back and fix this up.
!
  else

    do i = 1, 10
      j = ( i4_eleven - digit_pos ) * i + check
      if ( mod ( j, i4_eleven ) == 0 ) then
        k = i
      end if
    end do

  end if

  isbn(isbn_pos:isbn_pos) = i4_to_isbn ( k )

  return
end
function isbn_to_i4 ( c )

!*****************************************************************************80
!
!! ISBN_TO_I4 converts an ISBN character into an integer.
!
!  Discussion:
!
!    The characters '0' through '9' stand for themselves, but
!    the character 'X' or 'x' stands for 10.
!
!  Modified:
!
!    16 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Book Industry Study Group,
!    The Evolution in Product Identification:
!    Sunrise 2005 and the ISBN-13,
!    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
!
!  Parameters:
!
!    Input, character C, the ISBN character code to be converted.
!
!    Output, integer ( kind = 4 ) ISBN_TO_I4, the numeric value of the character
!    code, between 0 and 10.  This value is returned as -1 if C is
!    not a valid character code.
!
  implicit none

  character c
  integer ( kind = 4 ) isbn_to_i4

       if ( c == '0' ) then
    isbn_to_i4 = 0
  else if ( c == '1' ) then
    isbn_to_i4 = 1
  else if ( c == '2' ) then
    isbn_to_i4 = 2
  else if ( c == '3' ) then
    isbn_to_i4 = 3
  else if ( c == '4' ) then
    isbn_to_i4 = 4
  else if ( c == '5' ) then
    isbn_to_i4 = 5
  else if ( c == '6' ) then
    isbn_to_i4 = 6
  else if ( c == '7' ) then
    isbn_to_i4 = 7
  else if ( c == '8' ) then
    isbn_to_i4 = 8
  else if ( c == '9' ) then
    isbn_to_i4 = 9
  else if ( c == 'X' .or. c == 'x' ) then
    isbn_to_i4 = 10
  else
    isbn_to_i4 = -1
  end if

  return
end
function iset2_compare ( x1, y1, x2, y2 )

!*****************************************************************************80
!
!! ISET2_COMPARE compares two I2 sets.
!
!  Discussion:
!
!    The I2 set (X1,Y1) < (X2,Y2) if
!
!      min ( X1, Y1 ) < min ( X2, Y2 ) or
!      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) < max ( X2, Y2 )
!
!    The I2 set (X1,Y1) = (X2,Y2) if
!
!      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) = max ( X2, Y2 )
!
!  Modified:
!
!    18 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X1, Y1, the first I2 set.
!
!    Input, integer ( kind = 4 ) X2, Y2, the second I2 set.
!
!    Output, integer ( kind = 4 ) ISET2_COMPARE: 
!    -1, (X1,Y1) < (X2,Y2);
!     0, (X1,Y1) = (X2,Y2);
!    +1, (X1,Y1) > (X2,Y2).
!
  implicit none

  integer ( kind = 4 ) a1
  integer ( kind = 4 ) a2
  integer ( kind = 4 ) b1
  integer ( kind = 4 ) b2
  integer ( kind = 4 ) iset2_compare
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x1
  integer ( kind = 4 ) x2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  a1 = min ( x1, y1 )
  b1 = max ( x1, y1 )

  a2 = min ( x2, y2 )
  b2 = max ( x2, y2 )

  if ( a1 < a2 ) then
    value = -1
  else if ( a2 < a1 ) then
    value = +1
  else if ( b1 < b2 ) then
    value = -1
  else if ( b2 < b1 ) then
    value = +1
  else
    value = 0
  end if

  iset2_compare = value

  return
end
subroutine iset2_index_insert_unique ( n_max, n, x, y, indx, &
  xval, yval, ival, ierror )

!*****************************************************************************80
!
!! ISET2_INDEX_INSERT_UNIQUE inserts unique I2 set values in an indexed sorted list.
!
!  Discussion:
!
!    If the input value does not occur in the list, then N, X, Y and INDX
!    are updated.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_MAX, the maximum size of the list.
!
!    Input/output, integer ( kind = 4 ) N, the size of the list.
!
!    Input/output, integer ( kind = 4 ) X(N), Y(N), the list of I2 sets.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, YVAL, the value to be inserted if it is
!    not already in the list.
!
!    Output, integer ( kind = 4 ) IVAL, the index in INDX corresponding to the
!    value XVAL, YVAL.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 if an error occurred.
!
  implicit none

  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x(n_max)
  integer ( kind = 4 ) xval
  integer ( kind = 4 ) y(n_max)
  integer ( kind = 4 ) yval

  ierror = 0

  if ( n <= 0 ) then

    if ( n_max <= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ISET2_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    n = 1
    x(1) = min ( xval, yval )
    y(1) = max ( xval, yval )
    indx(1) = 1
    ival = 1
    return

  end if
!
!  Does ( XVAL, YVAL ) already occur in the list?
!
  call iset2_index_search ( n_max, n, x, y, indx, xval, yval, &
    less, equal, more )

  if ( equal == 0 ) then

    if ( n_max <= n ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ISET2_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    x(n+1) = min ( xval, yval )
    y(n+1) = max ( xval, yval )
    ival = more
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1

  else

    ival = equal

  end if

  return
end
subroutine iset2_index_search ( n_max, n, x, y, indx, xval, yval, &
  less, equal, more )

!*****************************************************************************80
!
!! ISET2_INDEX_SEARCH searches for an I2 set value in an indexed sorted list.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_MAX, the maximum size of the list.
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, integer ( kind = 4 ) X(N), Y(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, YVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) LESS, EQUAL, MORE, the indexes in INDX of the
!    list entries that are just less than, equal to, and just greater
!    than the test value.  If the test value does not occur in the list,
!    then EQUAL is zero.  If the test value is the minimum entry of the
!    list, then LESS is 0.  If the test value is the greatest entry of
!    the list, then MORE is N+1.
!
  implicit none

  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) compare
  integer ( kind = 4 ) equal
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) iset2_compare
  integer ( kind = 4 ) x(n_max)
  integer ( kind = 4 ) xhi
  integer ( kind = 4 ) xlo
  integer ( kind = 4 ) xmid
  integer ( kind = 4 ) xval
  integer ( kind = 4 ) y(n_max)
  integer ( kind = 4 ) yhi
  integer ( kind = 4 ) ylo
  integer ( kind = 4 ) ymid
  integer ( kind = 4 ) yval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n

  xlo = x(indx(lo))
  ylo = y(indx(lo))

  xhi = x(indx(hi))
  yhi = y(indx(hi))

  compare = iset2_compare ( xval, yval, xlo, ylo )

  if ( compare == -1 ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( compare == 0 ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  compare = iset2_compare ( xval, yval, xhi, yhi )

  if ( compare == +1 ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( compare == 0 ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))
    ymid = y(indx(mid))

    compare = iset2_compare ( xval, yval, xmid, ymid )

    if ( compare == 0 ) then
      equal = mid
      less = equal - 1
      more = equal + 1
      return
    else if ( compare == -1 ) then
      hi = mid
    else if ( compare == +1 ) then
      lo = mid
    end if

  end do

  return
end
function lcm_12n ( n )

!*****************************************************************************80
!
!! LCM_12N computes the least common multiple of the integers 1 through N.
!
!  Examples:
!
!    N    LCM_12N
!
!    1          1
!    2          2
!    3          3
!    4         12
!    5         60
!    6         60
!    7        420
!    8        840
!    9       2520
!   10       2520
!
!  Modified:
!
!    18 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the value of N.
!
!    Output, integer ( kind = 4 ) LCM_12N, the least common multiple of the integers 1 to N.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imult
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lcm_12n
  integer ( kind = 4 ) n

  lcm_12n = 1

  do i = 2, n

    imult = i

    do j = 1, i - 1

      if ( mod ( imult, ( i - j ) ) == 0 ) then
        imult = imult / ( i - j )
      end if

    end do

    lcm_12n = lcm_12n * imult

  end do

  return
end
subroutine lvec_print ( n, a, title )

!*****************************************************************************80
!
!! LVEC_PRINT prints an LVEC.
!
!  Modified:
!
!    12 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, logical ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  logical a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_1 = 1
  character ( len = * ) title

  call lvec_print_some ( n, a, i4_1, n, title )

  return
end
subroutine lvec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! LVEC_PRINT_SOME prints "some" of an LVEC.
!
!  Discussion:
!
!    An LVEC is a vector of logical values.
!
!  Modified:
!
!    12 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, logical A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices to 
!    print. The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  logical a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,l1)' ) i, a(i)
  end do

  return
end
function pause_input ( )

!*****************************************************************************80
!
!! PAUSE_INPUT waits until an input character is entered.
!
!  Modified:
!
!    02 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character PAUSE_INPUT, the character that was entered.
!
  implicit none

  integer ( kind = 4 ) ios
  character pause_input

  write ( *, '(a)' ) 'Press RETURN to continue.'
  read ( *, '(a)', iostat = ios ) pause_input

  return
end
subroutine perm_check ( n, p, ierror )

!*****************************************************************************80
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from 1
!    to N occurs among the N entries of the permutation.
!
!    If this is not the case, the routine stops.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) P(N), the array to check.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, the array represents a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifind
  integer ( kind = 4 ) iseek
  integer ( kind = 4 ) p(n)

  ierror = 0

  do iseek = 1, n

    ierror = iseek

    do ifind = 1, n
      if ( p(ifind) == iseek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  The input array does not represent'
      write ( *, '(a)' ) '  a proper permutation.  In particular, the'
      write ( *, '(a,i8)' ) '  array is missing the value ', ierror
      stop
    end if

  end do

  return
end
subroutine perm_cycle ( n, p, isgn, ncycle, iopt )

!*****************************************************************************80
!
!! PERM_CYCLE analyzes a permutation.
!
!  Discussion:
!
!    The routine will count cycles, find the sign of a permutation,
!    and tag a permutation.
!
!  Example:
!
!    Input:
!
!      N = 9
!      IOPT = 1
!      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      NCYCLE = 3
!      ISGN = +1
!      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
!
!  Modified:
!
!    09 July 2000
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) P(N).  On input, P describes a
!    permutation, in the sense that entry I is to be moved to P(I).
!    If IOPT = 0, then P will not be changed by this routine.
!    If IOPT = 1, then on output, P will be "tagged".  That is,
!    one element of every cycle in P will be negated.  In this way,
!    a user can traverse a cycle by starting at any entry I1 of P
!    which is negative, moving to I2 = ABS(P(I1)), then to
!    P(I2), and so on, until returning to I1.
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Output, integer ( kind = 4 ) ISGN, the "sign" of the permutation, which is
!    +1 if the permutation is even, -1 if odd.  Every permutation
!    may be produced by a certain number of pairwise switches.
!    If the number of switches is even, the permutation itself is
!    called even.
!
!    Output, integer ( kind = 4 ) NCYCLE, the number of cycles in the 
!    permutation.
!
!    Input, integer ( kind = 4 ) IOPT, requests tagging.
!    0, the permutation will not be tagged.
!    1, the permutation will be tagged.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ), parameter :: i4_two = 2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) is
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) ncycle
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )

  is = 1
  ncycle = n

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      ncycle = ncycle - 1
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    if ( iopt /= 0 ) then
      is = -i4_sign ( p(i) )
    end if

    p(i) = is * abs ( p(i) )

  end do

  isgn = 1 - 2 * mod ( n - ncycle, i4_two )

  return
end
subroutine perm_free ( npart, ipart, nfree, ifree )

!*****************************************************************************80
!
!! PERM_FREE reports the number of unused items in a partial permutation.
!
!  Discussion:
!
!    It is assumed that the N objects being permuted are the integers
!    from 1 to N, and that IPART contains a "partial" permutation, that
!    is, the NPART entries of IPART represent the beginning of a
!    permutation of all N items.
!
!    The routine returns in IFREE the items that have not been used yet.
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPART, the number of entries in IPART.  
!    NPART may be 0.
!
!    Input, integer ( kind = 4 ) IPART(NPART), the partial permutation, which 
!    should contain, at most once, some of the integers between 1 and
!    NPART+NFREE.
!
!    Input, integer ( kind = 4 ) NFREE, the number of integers that have not 
!    been used in IPART.  This is simply N - NPART.  NFREE may be zero.
!
!    Output, integer ( kind = 4 ) IFREE(NFREE), the integers between 1 and 
!    NPART+NFREE that were not used in IPART.
!
  implicit none

  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) npart

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree(nfree)
  integer ( kind = 4 ) ipart(npart)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) match
  integer ( kind = 4 ) n

  n = npart + nfree

  if ( npart < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NPART < 0.'
    write ( *, '(a,i8)' ) '  NPART = ', npart
    stop

  else if ( npart == 0 ) then

    call i4vec_indicator ( n, ifree )

  else if ( nfree < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NFREE < 0.'
    write ( *, '(a,i8)' ) '  NFREE = ', nfree
    stop

  else if ( nfree == 0 ) then

    return

  else

    k = 0

    do i = 1, n

      match = 0

      do j = 1, npart
        if ( ipart(j) == i ) then
          match = j
          exit
        end if
      end do

      if ( match == 0 ) then

        k = k + 1

        if ( nfree < k ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)'    ) 'PERM_FREE - Fatal error!'
          write ( *, '(a)'    ) '  The partial permutation is illegal.'
          write ( *, '(a)'    ) '  It should contain, at most once, some of'
          write ( *, '(a,i8)' ) '  the integers between 1 and N = ', n
          write ( *, '(a)'    ) '  The number of integers that have not'
          write ( *, '(a,i8)' ) '  been used is at least K = ', k
          write ( *, '(a,i8)' ) '  This should be exactly NFREE = ', nfree
          call i4vec_print ( npart, ipart, '  The partial permutation:' )
          stop
        end if

        ifree(k) = i

      end if

    end do

  end if

  return
end
subroutine perm_inverse ( n, p )

!*****************************************************************************80
!
!! PERM_INVERSE inverts a permutation "in place".
!
!  Modified:
!
!    02 January 2006
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard 
!    index form.  On output, P describes the inverse permutation
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ), parameter :: i4_one = 1
  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) is
  integer ( kind = 4 ) p(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    stop
  end if

  call perm_check ( n, p, ierror )

  is = 1

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    is = - i4_sign ( p(i) )
    p(i) = is * abs ( p(i) )

  end do

  do i = 1, n

    i1 = -p(i)

    if ( 0 <= i1 ) then

      i0 = i

      do

        i2 = p(i1)
        p(i1) = i0

        if ( i2 < 0 ) then
          exit
        end if

        i0 = i1
        i1 = i2

      end do

    end if

  end do

  return
end
subroutine perm_next ( n, p, more, even )

!*****************************************************************************80
!
!! PERM_NEXT computes all of the permutations on N objects, one at a time.
!
!  Discussion:
!
!    If this routine is called with MORE = TRUE, any
!    permutation in P, and EVEN = TRUE, then the successor of the input
!    permutation will be produced, unless P is the last permutation
!    on N letters, in which case P(1) will be set to 0 on return.
!
!  Modified:
!
!    12 March 2001
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N).
!    On input, P contains the previous permutation.
!    On output, P contains the next permutation.
!
!    Input/output, logical MORE.
!    On input, MORE = FALSE means this is the first call.
!    On output, MORE = FALSE means there are no more permutations.
!
!    Output, logical EVEN, is TRUE if the output permutation is even.
!
  implicit none

  integer ( kind = 4 ) n

  logical even
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ), parameter :: i4_one = 1
  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ), parameter :: i4_two = 2
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) id
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) p(n)

  if ( .not. more ) then

    call i4vec_indicator ( n, p )

    more = .true.
    even = .true.

    if ( n == 1 ) then
      more = .false.
      return
    end if

    if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, i4_two ) ) then
      return
    end if

    do i = 1, n - 3
      if ( p(i+1) /= p(i) + 1 ) then
        return
      end if
    end do

    more = .false.

  else

    if ( n == 1 ) then
      p(1) = 0
      more = .false.
      return
    end if

    if ( even ) then

      ia = p(1)
      p(1) = p(2)
      p(2) = ia
      even = .false.

      if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, i4_two ) ) then
        return
      end if

      do i = 1, n - 3
        if ( p(i+1) /= p(i) + 1 ) then
          return
        end if
      end do

      more = .false.
      return

    else

      is = 0
      more = .false.

      do i1 = 2, n

        ia = p(i1)
        i = i1 - 1
        id = 0

        do j = 1, i
          if ( i1 < p(j) ) then
            id = id + 1
          end if
        end do

        is = id + is

        if ( id /= i * mod ( is, i4_two ) ) then
          more = .true.
          exit
        end if

      end do

      if ( .not. more ) then
        p(1) = 0
        return
      end if

    end if

    m = mod ( is + i4_one, i4_two ) * ( n + 1 )

    do j = 1, i

      if ( i4_sign ( p(j) - ia ) /= i4_sign ( p(j) - m ) ) then
        m = p(j)
        l = j
      end if

    end do

    p(l) = ia
    p(i1) = m
    even = .true.

  end if

  return
end
subroutine perm_print ( n, p, title )

!*****************************************************************************80
!
!! PERM_PRINT prints a permutation.
!
!  Example:
!
!    Input:
!
!      P = 7 2 4 1 5 3 6
!
!    Printed output:
!
!      1 2 3 4 5 6 7
!      7 2 4 1 5 3 6
!
!  Modified:
!
!    16 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer which is to be partitioned.
!
!    Input, integer ( kind = 4 ) P(N), the permutation to be converted.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), parameter :: inc = 20
  integer ( kind = 4 ) p(n)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do ilo = 1, n, inc
    ihi = min ( n, ilo + inc - 1 )
    write ( *, '(a)' ) ' '
    write ( *, '(20i4)' ) ( i,i = ilo, ihi )
    write ( *, '(20i4)' ) p(ilo:ihi)
  end do

  return
end
subroutine perm_uniform ( n, seed, p )

!*****************************************************************************80
!
!! PERM_UNIFORM selects a random permutation of N objects.
!
!  Modified:
!
!    01 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) P(N), the permutation.  P(I) is the "new"
!    location of the object originally at I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  call i4vec_indicator ( n, p )

  do i = 1, n
    j = i4_uniform ( i, n, seed )
    call i4_swap ( p(i), p(j) )
  end do

  return
end
function pounds_to_kilograms ( lb )

!*****************************************************************************80
!
!! POUNDS_TO_KILOGRAMS converts a measurement in pounds to kilograms.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) LB, the weight in pounds.
!
!    Output, real ( kind = 8 ) POUNDS_TO_KILOGRAMS, the corresponding 
!    weight in kilograms.
!
  implicit none

  real    ( kind = 8 ) lb
  real    ( kind = 8 ) pounds_to_kilograms

  pounds_to_kilograms = 0.4535924D+00 * lb

  return
end
function prime ( n )

!*****************************************************************************80
!
!! PRIME returns any of the first PRIME_MAX prime numbers.
!
!  Discussion:
!
!    PRIME_MAX is 1600, and the largest prime stored is 13499.
!
!    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
!
!  Modified:
!
!    18 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964, pages 870-873.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 95-98.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the desired prime number.
!    In general, is should be true that 0 <= N <= PRIME_MAX.
!    N = -1 returns PRIME_MAX, the index of the largest prime available.
!    N = 0 is legal, returning PRIME = 1.
!
!    Output, integer ( kind = 4 ) PRIME, the N-th prime.  If N is out of range, 
!    PRIME is returned as -1.
!
  implicit none

  integer ( kind = 4 ), parameter :: prime_max = 1600

  integer ( kind = 4 ), save :: icall = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( prime_max ) :: npvec
  integer ( kind = 4 ) prime

  if ( icall == 0 ) then

    icall = 1

    npvec(1:100) = (/ &
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /)

    npvec(101:200) = (/ &
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /)

    npvec(201:300) = (/ &
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /)

    npvec(301:400) = (/ &
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /)

    npvec(401:500) = (/ &
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /)

    npvec(501:600) = (/ &
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /)

    npvec(601:700) = (/ &
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /)

    npvec(701:800) = (/ &
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /)

    npvec(801:900) = (/ &
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /)

    npvec(901:1000) = (/ &
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)

    npvec(1001:1100) = (/ &
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, &
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, &
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, &
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, &
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, &
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, &
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, &
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, &
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, &
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /)

    npvec(1101:1200) = (/ &
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, &
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, &
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, &
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, &
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, &
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, &
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, &
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, &
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, &
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /)

    npvec(1201:1300) = (/ &
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, &
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, &
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, &
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, &
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, &
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, &
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, &
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, &
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, &
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /)

    npvec(1301:1400) = (/ &
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, &
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, &
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949, &
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, &
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, &
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, &
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, &
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, &
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, &
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /)

    npvec(1401:1500) = (/ &
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777, &
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833, &
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933, &
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011, &
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109, &
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211, &
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289, &
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401, &
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487, &
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /)

   npvec(1501:1600) = (/ &
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, &
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, &
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, &
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, &
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, &
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, &
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, &
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, &
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, &
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 /)

  end if

  if ( n == -1 ) then
    prime = prime_max
  else if ( n == 0 ) then
    prime = 1
  else if ( n <= prime_max ) then
    prime = npvec(n)
  else
    prime = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRIME - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i8)' ) '  N should be between 1 and PRIME_MAX =', prime_max
    stop
  end if

  return
end
function prime_ge ( n )

!*****************************************************************************80
!
!! PRIME_GE returns the smallest prime greater than or equal to N.
!
!  Examples:
!
!    N     PRIME_GE
!
!    -10    2
!      1    2
!      2    2
!      3    3
!      4    5
!      5    5
!      6    7
!      7    7
!      8   11
!      9   11
!     10   11
!
!  Modified:
!
!    09 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number to be bounded.
!
!    Output, integer ( kind = 4 ) PRIME_GE, the smallest prime number that is 
!    greater than or equal to N.  However, if N is larger than the largest
!    prime stored, then PRIME_GE is returned as -1.
!
  implicit none

  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  integer ( kind = 4 ) i_mid
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p_hi
  integer ( kind = 4 ) p_lo
  integer ( kind = 4 ) p_mid
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) prime_ge

  if ( n <= 2 ) then
    prime_ge = 2
    return
  end if

  i_lo = 1
  p_lo = prime(i_lo)
  i_hi = prime(-1)
  p_hi = prime(i_hi)

  if ( p_hi < n ) then
    prime_ge = -p_hi
    return
  end if

  do

    if ( i_lo + 1 == i_hi ) then
      prime_ge = p_hi
      return
    end if

    i_mid = ( i_lo + i_hi ) / 2
    p_mid = prime(i_mid)

    if ( p_mid < n ) then
      i_lo = i_mid
      p_lo = p_mid
    else if ( n <= p_mid ) then
      i_hi = i_mid
      p_hi = p_mid
    end if

  end do

  return
end
subroutine primer ( n, iprime )

!*****************************************************************************80
!
!! PRIMER computes the prime numbers up to a given limit.
!
!  Discussion:
!
!    PRIMER returns the results of its computations in the vector
!    IPRIME.  IPRIME(I) is -1 if the number I is not prime, and
!    1 if I is prime.
!
!    The algorithm is a simple-minded sieve of Eristothenes, with
!    no attempt at efficiency.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of IPRIME, and the maximum
!    value that will be considered.
!
!    Output, integer ( kind = 4 ) IPRIME(N), records the results for each 
!    integer.  IPRIME(I) = -1 if I is not prime, and IPRIME(I) = 1 if I is
!    prime.  By convention, IPRIME(1) will be set to -1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iprime(n)
  integer ( kind = 4 ) next
!
!  IPRIME(I) = 0 means we don't know if I is prime.
!
  iprime(1:n) = 0
!
!  By convention, 1 is not prime.
!
  iprime(1) = -1
  next = 1
!
!  Examine the integers in order.
!
  do next = 2, n

    if ( iprime(next) == 0 ) then
      iprime(next) = 1
      do i = 2 * next, n, next
        iprime(i) = -1
      end do
    end if

  end do

  return
end
function r4_epsilon ( )

!*****************************************************************************80
!
!! R4_EPSILON returns the R4 roundoff unit.
!
!  Discussion:
!
!    The roundoff unit is a number R which is a power of 2 with the 
!    property that, to the precision of the computer's arithmetic,
!      1 < 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Modified:
!
!    08 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) R4_EPSILON, the R4 round-off unit.
!
  implicit none

  real    ( kind = 4 ) d
  real    ( kind = 4 ) d_test
  real    ( kind = 4 ) r4_epsilon

  d = 1.0E+00
  d_test = 1.0E+00 + d / 2.0E+00

  do while ( 1.0E+00 < d_test )
    d = d / 2.0E+00
    d_test = 1.0E+00 + d / 2.0E+00
  end do

  r4_epsilon = d

  return
end
function r4_exp ( x )

!*****************************************************************************80
!
!! R4_EXP computes the exponential function, avoiding overflow and underflow.
!
!  Discussion:
!
!    My experience with the G95 compiler has included many unpleasant
!    floating point exceptions when very small arguments are given to
!    the exponential function.
!
!    This routine is designed to avoid such problems.
!
!    Ideally, the rule would be:
!
!                    X <= log ( TINY ) => R4_EXP ( X ) = 0
!    log ( HUGE ) <= X                 => R4_EXP ( X ) = HUGE
!
!    However, the G95 math library seems to produce infinity for
!    EXP ( LOG ( HUGE ( X ) ), rather than HUGE ( X ), so we've
!    included a fudge factor.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X, the argument of the exponential function.
!
!    Output, real ( kind = 4 ) R4_EXP, the value of exp ( X ).
!
  implicit none

  real    ( kind = 4 ), parameter :: log_max = 88.71397E+00
  real    ( kind = 4 ), parameter :: log_min = -87.34528E+00
  real    ( kind = 4 ) r4_exp
  real    ( kind = 4 ) x

  if ( x <= log_min ) then
    r4_exp = 0.0E+00
  else if ( x < log_max ) then
    r4_exp = exp ( x )
  else
    r4_exp = huge ( x )
  end if

  return
end
function r4_huge ( )

!*****************************************************************************80
!
!! R4_HUGE returns the largest legal R4.
!
!  Discussion:
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    is more suitable for this purpose.
!
!  Modified:
!
!    08 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) R4_HUGE, a "huge" value.
!
  implicit none

  real    ( kind = 4 ) r4_huge

  r4_huge = 0.3402823466385E+39

  return
end
function r4_tiny ( )

!*****************************************************************************80
!
!! R4_TINY returns the smallest positive R4.
!
!  Discussion:
!
!    FORTRAN90 provides a built-in routine TINY ( X ) that
!    is more suitable for this purpose.
!
!  Modified:
!
!    08 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) R4_TINY, a "tiny" value.
!
  implicit none

  real    ( kind = 4 ) r4_tiny

  r4_tiny = 0.1175494350822E-37

  return
end
function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r4_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 4 ) r4_uniform_01

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
function r8_abs ( x )

!*****************************************************************************80
!
!! R8_ABS returns the absolute value of an R8.
!
!  Discussion:
!
!    FORTRAN90 supplies the ABS function, which should be used instead
!    of this function!
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose absolute value is desired.
!
!    Output, real ( kind = 8 ) R8_ABS, the absolute value of X.
!
  implicit none

  real    ( kind = 8 ) r8_abs
  real    ( kind = 8 ) x

  if ( 0.0D+00 <= x ) then
    r8_abs = x
  else
    r8_abs = -x
  end if

  return
end
function r8_cas ( x )

!*****************************************************************************80
!
!! R8_CAS returns the "casine" of an R8.
!
!  Discussion:
!
!    The "casine", used in the discrete Hartley transform, is abbreviated
!    CAS(X), and defined by:
!
!      CAS(X) = cos ( X ) + sin( X )
!             = sqrt ( 2 ) * sin ( X + pi/4 )
!             = sqrt ( 2 ) * cos ( X - pi/4 )
!
!  Modified:
!
!    06 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ralph Hartley,
!    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
!    Proceedings of the Institute of Radio Engineers,
!    Volume 30, pages 144-150, 1942.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose casine is desired.
!
!    Output, real ( kind = 8 ) R8_CAS, the casine of X, which will be between
!    plus or minus the square root of 2.
!
  implicit none

  real    ( kind = 8 ) r8_cas
  real    ( kind = 8 ) x

  r8_cas = cos ( x ) + sin ( x )

  return
end
function r8_ceiling ( r )

!*****************************************************************************80
!
!! R8_CEILING rounds an R8 "up" (towards +infinity) to the next integer.
!
!  Examples:
!
!    R     Value
!
!   -1.1  -1
!   -1.0  -1
!   -0.9   0
!    0.0   0
!    5.0   5
!    5.1   6
!    5.9   6
!    6.0   6
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the value to be rounded up.
!
!    Output, integer ( kind = 4 ) R8_CEILING, the rounded value.
!
  implicit none

  real    ( kind = 8 ) r
  integer ( kind = 4 ) r8_ceiling
  integer ( kind = 4 ) value

  value = int ( r )
  if ( real ( value, kind = 8 ) < r ) then
    value = value + 1
  end if

  r8_ceiling = value

  return
end
function r8_chop ( place, x )

!*****************************************************************************80
!
!! R8_CHOP chops an R8 to a given number of binary places.
!
!  Example:
!
!    3.875 = 2 + 1 + 1/2 + 1/4 + 1/8.
!
!    The following values would be returned for the 'chopped' value of
!    3.875:
!
!    PLACE  Value
!
!       1      2
!       2      3     = 2 + 1
!       3      3.5   = 2 + 1 + 1/2
!       4      3.75  = 2 + 1 + 1/2 + 1/4
!       5+     3.875 = 2 + 1 + 1/2 + 1/4 + 1/8
!
!  Modified:
!
!    20 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PLACE, the number of binary places to preserve.
!    PLACE = 0 means return the integer part of X.
!    PLACE = 1 means return the value of X, correct to 1/2.
!    PLACE = 2 means return the value of X, correct to 1/4.
!    PLACE = -1 means return the value of X, correct to 2.
!
!    Input, real ( kind = 8 ) X, the number to be chopped.
!
!    Output, real ( kind = 8 ) R8_CHOP, the chopped number.
!
  implicit none

  real    ( kind = 8 ) fac
  integer ( kind = 4 ) place
  real    ( kind = 8 ) r8_chop
  real    ( kind = 8 ) r8_log_2
  real    ( kind = 8 ) r8_sign
  real    ( kind = 8 ) s
  integer ( kind = 4 ) temp
  real    ( kind = 8 ) x

  s = r8_sign ( x )
  temp = int ( r8_log_2 ( abs ( x ) ) )
  fac = 2.0D+00**( temp - place + 1 )
  r8_chop = s * real ( int ( abs ( x ) / fac ), kind = 8 ) * fac

  return
end
function r8_cube_root ( x )

!*****************************************************************************80
!
!! R8_CUBE_ROOT returns the cube root of an R8.
!
!  Discussion:
!
!    This routine is designed to avoid the possible problems that can occur
!    when formulas like 0.0**(1/3) or (-1.0)**(1/3) are to be evaluated.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose cube root is desired.
!
!    Output, real ( kind = 8 ) R8_CUBE_ROOT, the cube root of X.
!
  implicit none

  real    ( kind = 8 ) r8_cube_root
  real    ( kind = 8 ) value
  real    ( kind = 8 ) x

  if ( 0.0D+00 < x ) then
    value = x**(1.0D+00/3.0D+00)
  else if ( x == 0.0D+00 ) then
    value = 0.0D+00
  else
    value = -( abs ( x ) )**(1.0D+00/3.0D+00)
  end if

  r8_cube_root = value

  return
end
function r8_diff ( x, y, n )

!*****************************************************************************80
!
!! R8_DIFF computes the difference of two R8's to a specified accuracy.
!
!  Discussion:
!
!    The user controls how many binary digits of accuracy
!    are to be used.
!
!    N determines the accuracy of the value of the result.  If N = 10,
!    for example, only 11 binary places will be used in the arithmetic.
!    In general, only N+1 binary places will be used.
!
!    N may be zero.  However, a negative value of N should
!    not be used, since this will cause both X and Y to look like 0.
!
!  Modified:
!
!    17 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the two values whose difference is desired.
!
!    Input, integer ( kind = 4 ) N, the number of binary digits to use.
!
!    Output, real ( kind = 8 ) R8_DIFF, the value of X-Y.
!
  implicit none

  real    ( kind = 8 ) cx
  real    ( kind = 8 ) cy
  integer ( kind = 4 ) n
  real    ( kind = 8 ) pow2
  real    ( kind = 8 ) r8_diff
  real    ( kind = 8 ) size
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y

  if ( x == y ) then
    r8_diff = 0.0D+00
    return
  end if

  pow2 = 2.0D+00**n
!
!  Compute the magnitude of X and Y, and take the larger of the
!  two.  At least one of the two values is not zero!
!
  size = max ( abs ( x ), abs ( y ) )
!
!  Make normalized copies of X and Y.  One of the two values will
!  actually be equal to 1.
!
  cx = x / size
  cy = y / size
!
!  Here's where rounding comes in.  We know that the larger of the
!  the two values equals 1.  We multiply both values by 2**N,
!  where N+1 is the number of binary digits of accuracy we want
!  to use, truncate the values, and divide back by 2**N.
!
  cx = real ( int ( cx * pow2 + sign ( 0.5D+00, cx ) ), kind = 8 ) / pow2
  cy = real ( int ( cy * pow2 + sign ( 0.5D+00, cy ) ), kind = 8 ) / pow2
!
!  Take the difference now.
!
  r8_diff = cx - cy
!
!  Undo the scaling.
!
  r8_diff = r8_diff * size

  return
end
subroutine r8_digit ( x, idigit, digit )

!*****************************************************************************80
!
!! R8_DIGIT returns a particular decimal digit of an R8.
!
!  Modified:
!
!    17 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose NDIG-th decimal digit
!    is desired.  If X is zero, all digits will be returned as 0.
!
!    Input, integer ( kind = 4 ) IDIGIT, the position of the desired decimal 
!    digit.  A value of 1 means the leading digit, a value of 2 the second digit
!    and so on.
!
!    Output, integer ( kind = 4 ) DIGIT, the value of the IDIGIT-th decimal 
!    digit of X.
!
  implicit none

  integer ( kind = 4 ) digit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idigit
  integer ( kind = 4 ) ival
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xcopy

  if ( x == 0.0D+00 ) then
    digit = 0
    return
  end if

  if ( idigit <= 0 ) then
    digit = 0
    return
  end if
!
!  Set XCOPY = X, and then force XCOPY to lie between 1 and 10.
!
  xcopy = abs ( x )

  do while ( xcopy < 1.0D+00 )
    xcopy = xcopy * 10.0D+00
  end do

  do while ( 10.0D+00 <= xcopy )
    xcopy = xcopy / 10.0D+00
  end do

  do i = 1, idigit
    ival = int ( xcopy )
    xcopy = ( xcopy - ival ) * 10.0D+00
  end do

  digit = ival

  return
end
function r8_epsilon ( )

!*****************************************************************************80
!
!! R8_EPSILON returns the R8 roundoff unit.
!
!  Discussion:
!
!    The roundoff unit is a number R which is a power of 2 with the 
!    property that, to the precision of the computer's arithmetic,
!      1 < 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Modified:
!
!    15 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_EPSILON, the R8 round-off unit.
!
  implicit none

  real    ( kind = 8 ) d
  real    ( kind = 8 ) d_test
  real    ( kind = 8 ) r8_epsilon

  d = 1.0D+00
  d_test = 1.0D+00 + d / 2.0D+00

  do while ( 1.0D+00 < d_test )
    d = d / 2.0D+00
    d_test = 1.0D+00 + d / 2.0D+00
  end do

  r8_epsilon = d

  return
end
function r8_exp ( x )

!*****************************************************************************80
!
!! R8_EXP computes the exponential function, avoiding overflow and underflow.
!
!  Discussion:
!
!    My experience with the G95 compiler has included many unpleasant
!    floating point exceptions when very small arguments are given to
!    the exponential function.
!
!    This routine is designed to avoid such problems.
!
!    Ideally, the rule would be:
!
!                    X <= log ( TINY ) => R8_EXP ( X ) = 0
!    log ( HUGE ) <= X                 => R8_EXP ( X ) = HUGE
!
!    However, the G95 math library seems to produce infinity for
!    EXP ( LOG ( HUGE ( X ) ), rather than HUGE ( X ), so we've
!    included a fudge factor.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the exponential function.
!
!    Output, real ( kind = 8 ) R8_EXP, the value of exp ( X ).
!
  implicit none

  real    ( kind = 8 ), parameter :: log_max = 709.711D+00
  real    ( kind = 8 ), parameter :: log_min = -708.467D+00
  real    ( kind = 8 ) r8_exp
  real    ( kind = 8 ) x

  if ( x <= log_min ) then
    r8_exp = 0.0D+00
  else if ( x < log_max ) then
    r8_exp = exp ( x )
  else
    r8_exp = huge ( x )
  end if

  return
end
function r8_floor ( r )

!*****************************************************************************80
!
!! R8_FLOOR rounds an R8 "down" (towards -infinity) to the next integer.
!
!  Examples:
!
!    R     Value
!
!   -1.1  -2
!   -1.0  -1
!   -0.9  -1
!    0.0   0
!    5.0   5
!    5.1   5
!    5.9   5
!    6.0   6
!
!  Modified:
!
!    15 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the value to be rounded down.
!
!    Output, integer ( kind = 4 ) R8_FLOOR, the rounded value.
!
  implicit none

  real    ( kind = 8 ) r
  integer ( kind = 4 ) r8_floor
  integer ( kind = 4 ) value

  value = int ( r )
  if ( r < real ( value, kind = 8 ) ) then
    value = value - 1
  end if

  r8_floor = value

  return
end
function r8_fraction ( x )

!*****************************************************************************80
!
!! R8_FRACTION returns the fraction part of an R8.
!
!  Discussion:
!
!    If we regard a real number as
!
!      R8 = SIGN * ( WHOLE + FRACTION )
!
!    where
!
!      SIGN is +1 or -1,
!      WHOLE is a nonnegative integer
!      FRACTION is a nonnegative real number strictly less than 1,
!
!    then this routine returns the value of FRACTION.
!
!  Example:
!
!     R8    R8_FRACTION
! 
!    0.00      0.00
!    1.01      0.01
!    2.02      0.02
!   19.73      0.73
!   -4.34      0.34
!    
!  Modified:
!
!    16 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) R8_FRACTION, the fraction part of X.
!
  implicit none

  real    ( kind = 8 ) r8_fraction
  real    ( kind = 8 ) x

  r8_fraction = abs ( x ) - real ( int ( abs ( x ) ), kind = 8 )

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_HUGE, a "huge" value.
!
  implicit none

  real    ( kind = 8 ) r8_huge

  r8_huge = 1.0D+30

  return
end
function r8_in_01 ( a )

!*****************************************************************************80
!
!! R8_IN_01 is TRUE if an R8 is in the range [0,1].
!
!  Modified:
!
!    06 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the value.
!
!    Output, logical R8_IN_01, is TRUE if 0 <= A <= 1.
!
  implicit none

  real    ( kind = 8 ) a
  logical r8_in_01
  logical value

  if ( a < 0.0D+00 .or. 1.0D+00 < a ) then
    value = .false.
  else
    value = .true.
  end if

  r8_in_01 = value

  return
end
function r8_is_int ( r )

!*****************************************************************************80
!
!! R8_IS_INT determines if an R8 represents an integer value.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the number to be checked.
!
!    Output, logical R8_IS_INT, is TRUE if R is an integer value.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  real    ( kind = 8 ) r
  logical r8_is_int
  logical value

  if ( real ( i4_huge ( ), kind = 8 ) < r ) then
    value = .false.
  else if ( r < - real ( i4_huge ( ), kind = 8 ) ) then
    value = .false.
  else if ( r == real ( int ( r ), kind = 8 ) ) then
    value = .true.
  else
    value = .false.
  end if

  r8_is_int = value

  return
end
function r8_log_2 ( x )

!*****************************************************************************80
!
!! R8_LOG_2 returns the logarithm base 2 of an R8.
!
!  Discussion:
!
!    value = Log ( |X| ) / Log ( 2.0 )
!
!  Modified:
!
!    27 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose base 2 logarithm is desired.
!    X should not be 0.
!
!    Output, real ( kind = 8 ) R8_LOG_2, the logarithm base 2 of the absolute
!    value of X.  It should be true that |X| = 2**D_LOG_2.
!
  implicit none

  real    ( kind = 8 ) r8_log_2
  real    ( kind = 8 ) x

  if ( x == 0.0D+00 ) then
    r8_log_2 = -huge ( x )
  else
    r8_log_2 = log ( abs ( x ) ) / log ( 2.0D+00 )
  end if

  return
end
function r8_log_10 ( x )

!*****************************************************************************80
!
!! R8_LOG_10 returns the logarithm base 10 of an R8.
!
!  Discussion:
!
!    value = Log10 ( |X| )
!
!  Modified:
!
!    27 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose base 2 logarithm is desired.
!    X should not be 0.
!
!    Output, real ( kind = 8 ) R8_LOG_10, the logarithm base 10 of the absolute
!    value of X.  It should be true that |X| = 10**R_LOG_10.
!
  implicit none

  real    ( kind = 8 ) r8_log_10
  real    ( kind = 8 ) x

  if ( x == 0.0D+00 ) then
    r8_log_10 = -huge ( x )
  else
    r8_log_10 = log10 ( abs ( x ) )
  end if

  return
end
function r8_log_b ( x, b )

!*****************************************************************************80
!
!! R8_LOG_B returns the logarithm base B of an R8.
!
!  Discussion:
!
!    value = log ( |X| ) / log ( |B| )
!
!  Modified:
!
!    27 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose base B logarithm is desired.
!    X should not be 0.
!
!    Input, real ( kind = 8 ) B, the base, which should not be 0, 1 or -1.
!
!    Output, real ( kind = 8 ) R8_LOG_B, the logarithm base B of the absolute
!    value of X.  It should be true that |X| = |B|**D_LOG_B.
!
  implicit none

  real    ( kind = 8 ) b
  real    ( kind = 8 ) r8_log_b
  real    ( kind = 8 ) x

  if ( b == 0.0D+00 .or. b == 1.0D+00 .or. b == -1.0D+00 ) then
    r8_log_b = -huge ( x )
  else if ( abs ( x ) == 0.0D+00 ) then
    r8_log_b = -huge ( x )
  else
    r8_log_b = log ( abs ( x ) ) / log ( abs ( b ) )
  end if

  return
end
subroutine r8_mant ( x, s, r, l )

!*****************************************************************************80
!
!! R8_MANT computes the "mantissa" or "fraction part" of an R8.
!
!  Formula:
!
!    X = S * R * 2**L
!
!    S is +1 or -1,
!    R is a double precision value between 1.0 and 2.0,
!    L is an integer.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be decomposed.
!
!    Output, integer ( kind = 4 ) S, the "sign" of the number.
!    S will be -1 if X is less than 0, and +1 if X is greater
!    than or equal to zero.
!
!    Output, real ( kind = 8 ) R, the mantissa of X.  R will be greater
!    than or equal to 1, and strictly less than 2.  The one
!    exception occurs if X is zero, in which case R will also
!    be zero.
!
!    Output, integer ( kind = 4 ) L, the integer part of the logarithm 
!    (base 2) of X.
!
  implicit none

  integer ( kind = 4 ) l
  real    ( kind = 8 ) r
  integer ( kind = 4 ) s
  real    ( kind = 8 ) x
!
!  Determine the sign.
!
  if ( x < 0.0D+00 ) then
    s = -1
  else
    s = 1
  end if
!
!  Set R to the absolute value of X, and L to zero.
!  Then force R to lie between 1 and 2.
!
  if ( x < 0.0D+00 ) then
    r = -x
  else
    r = x
  end if

  l = 0
!
!  Time to bail out if X is zero.
!
  if ( x == 0.0D+00 ) then
    return
  end if

  do while ( 2.0D+00 <= r )
    r = r / 2.0D+00
    l = l + 1
  end do

  do while ( r < 1.0D+00 )
    r = r * 2.0D+00
    l = l - 1
  end do

  return
end
function r8_mod ( x, y )

!*****************************************************************************80
!
!! R8_MOD returns the remainder of R8 division.
!
!  Formula:
!
!    If
!      REM = R8_MOD ( X, Y )
!      RMULT = ( X - REM ) / Y
!    then
!      X = Y * RMULT + REM
!    where REM has the same sign as X, and abs ( REM ) < Y.
!
!  Examples:
!
!        X         Y     R8_MOD  R8_MOD Factorization
!
!      107        50       7      107 =  2 *  50 + 7
!      107       -50       7      107 = -2 * -50 + 7
!     -107        50      -7     -107 = -2 *  50 - 7
!     -107       -50      -7     -107 =  2 * -50 - 7
!
!  Modified:
!
!    14 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be divided.
!
!    Input, real ( kind = 8 ) Y, the number that divides X.
!
!    Output, real ( kind = 8 ) R8_MOD, the remainder when X is divided by Y.
!
  implicit none

  real    ( kind = 8 ) r8_mod
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y

  if ( y == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_MOD - Fatal error!'
    write ( *, '(a,g14.6)' ) '  R8_MODP ( X, Y ) called with Y = ', y
    stop
  end if

  r8_mod = x - real ( int ( x / y ), kind = 8 ) * y

  if ( x < 0.0D+00 .and. 0.0D+00 < r8_mod ) then
    r8_mod = r8_mod - abs ( y )
  else if ( 0.0D+00 < x .and. r8_mod < 0.0D+00 ) then
    r8_mod = r8_mod + abs ( y )
  end if

  return
end
function r8_modp ( x, y )

!*****************************************************************************80
!
!! R8_MODP returns the nonnegative remainder of R8 division.
!
!  Formula:
!
!    If
!      REM = R8_MODP ( X, Y )
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
!    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
!
!  Examples:
!
!        X         Y     MOD R8_MODP  R8_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be divided.
!
!    Input, real ( kind = 8 ) Y, the number that divides X.
!
!    Output, real ( kind = 8 ) R8_MODP, the nonnegative remainder
!    when X is divided by Y.
!
  implicit none

  real    ( kind = 8 ) r8_modp
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y

  if ( y == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_MODP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  R8_MODP ( X, Y ) called with Y = ', y
    stop
  end if

  r8_modp = mod ( x, y )

  if ( r8_modp < 0.0D+00 ) then
    r8_modp = r8_modp + abs ( y )
  end if

  return
end
function r8_nint ( x )

!*****************************************************************************80
!
!! R8_NINT returns the nearest integer to an R8.
!
!  Examples:
!
!        X        R8_NINT
!
!      1.3         1
!      1.4         1
!      1.5         1 or 2
!      1.6         2
!      0.0         0
!     -0.7        -1
!     -1.1        -1
!     -1.6        -2
!
!  Modified:
!
!    08 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value.
!
!    Output, integer ( kind = 4 ) R8_NINT, the nearest integer to X.
!
  implicit none

  integer ( kind = 4 ) r8_nint
  integer ( kind = 4 ) s
  real    ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    s = -1
  else
    s = 1
  end if

  r8_nint = s * int ( abs ( x ) + 0.5D+00 )

  return
end
function r8_normal ( a, b, seed )

!*****************************************************************************80
!
!! R8_NORMAL returns a scaled pseudonormal R8.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the mean of the PDF.
!
!    Input, real ( kind = 8 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL, a sample of the normal PDF.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  integer ( kind = 4 ), parameter :: i4_two = 2
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) r1
  real    ( kind = 8 ) r2
  real    ( kind = 8 ) r8_normal
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save :: seed2 = 0
  integer ( kind = 4 ), save :: used = 0
  real    ( kind = 8 ) x
  real    ( kind = 8 ), save :: y = 0.0D+00
!
!  On odd numbered calls, generate two uniforms, create two normals,
!  return the first normal and its corresponding seed.
!
  if ( mod ( used, i4_two ) == 0 ) then

    r1 = r8_uniform_01 ( seed )

    if ( r1 == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_NORMAL - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    seed2 = seed
    r2 = r8_uniform_01 ( seed2 )

    x = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    y = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )
!
!  On odd calls, return the second normal and its corresponding seed.
!
  else

    seed = seed2
    x = y

  end if

  used = used + 1

  r8_normal = a + b * x

  return
end
function r8_normal_01 ( seed )

!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    Because this routine uses the Box Muller method, it requires pairs
!    of uniform random values to generate a pair of normal random values.
!    This means that on every other call, essentially, the input value of
!    SEED is ignored, since the code saves the second normal random value.
!
!    If you didn't know this, you might be confused since, usually, the
!    output of a random number generator can be completely controlled by
!    the input value of the SEED.  If I were more careful, I could rewrite
!    this routine so that it would distinguish between cases where the input
!    value of SEED is the output value from the previous call (all is well)
!    and those cases where it is not (the user has decided to do something
!    new.  Restart the uniform random number sequence.)  But I'll leave
!    that for later.
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_01, a sample of the standard 
!    normal PDF.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_2 = 2
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) r1
  real    ( kind = 8 ) r2
  real    ( kind = 8 ) r8_normal_01
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save :: seed2 = 0
  integer ( kind = 4 ), save :: used = 0
  real    ( kind = 8 ) x
  real    ( kind = 8 ), save :: y = 0.0D+00
!
!  On odd numbered calls, generate two uniforms, create two normals,
!  return the first normal and its corresponding seed.
!
  if ( mod ( used, i4_2 ) == 0 ) then

    r1 = r8_uniform_01 ( seed )

    if ( r1 == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    seed2 = seed
    r2 = r8_uniform_01 ( seed2 )

    x = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    y = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )
!
!  On odd calls, return the second normal and its corresponding seed.
!
  else

    seed = seed2
    x = y

  end if

  used = used + 1

  r8_normal_01 = x

  return
end
function r8_pi ( )

!*****************************************************************************80
!
!! R8_PI returns the value of pi as an R8.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_PI, the value of pi.
!
  implicit none

  real    ( kind = 8 ) r8_pi

  r8_pi = 3.141592653589793D+00

  return
end
function r8_power ( r, p )

!*****************************************************************************80
!
!! R8_POWER computes the P-th power of an R8.
!
!  Modified:
!
!    04 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the base.
!
!    Input, integer ( kind = 4 ) P, the power, which may be negative.
!
!    Output, real ( kind = 8 ) R8_POWER, the value of the P-th power of R.
!
  implicit none

  integer ( kind = 4 ) p
  real    ( kind = 8 ) r
  real    ( kind = 8 ) r8_power
  real    ( kind = 8 ) value
!
!  Special case.  R^0 = 1.
!
  if ( p == 0 ) then

    value = 1.0D+00
!
!  Special case.  Positive powers of 0 are 0.
!  For negative powers of 0, we go ahead and compute R**P,
!  relying on the software to complain.
!
  else if ( r == 0.0D+00 ) then

    if ( 0 < p ) then
      value = 0.0D+00
    else
      value = r**p
    end if

  else if ( 1 <= p ) then
    value = r**p
  else
    value = 1.0D+00 / r**(-p)
  end if

  r8_power = value

  return
end
subroutine r8_power_fast ( r, p, rp, mults )

!*****************************************************************************80
!
!! R8_POWER_FAST computes an integer power of an R8.
!
!  Discussion:
!
!    Obviously, R**P can be computed using P-1 multiplications.
!
!    However, R**P can also be computed using at most 2*LOG2(P) multiplications.
!    To do the calculation this way, let N = LOG2(P).
!    Compute A, A**2, A**4, ..., A**N by N-1 successive squarings.
!    Start the value of R**P at A, and each time that there is a 1 in
!    the binary expansion of P, multiply by the current result of the squarings.
!
!    This algorithm is not optimal.  For small exponents, and for special
!    cases, the result can be computed even more quickly.
!
!  Modified:
!
!    30 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the base.
!
!    Input, integer ( kind = 4 ) P, the power, which may be negative.
!
!    Output, real ( kind = 8 ) RP, the value of R**P.
!
!    Output, integer ( kind = 4 ) MULTS, the number of multiplications 
!    and divisions.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_one = 1
  integer ( kind = 4 ), parameter :: i4_two = 2
  integer ( kind = 4 ) mults
  integer ( kind = 4 ) p
  integer ( kind = 4 ) p_mag
  integer ( kind = 4 ) p_sign
  real    ( kind = 8 ) r
  real    ( kind = 8 ) r2
  real    ( kind = 8 ) rp

  mults = 0
!
!  Special bases.
!
  if ( r == 1.0D+00 ) then
    rp = 1.0D+00
    return
  end if

  if ( r == -1.0D+00 ) then

    if ( mod ( p, i4_two ) == 1 ) then
      rp = -1.0D+00
    else
      rp = 1.0D+00
    end if

    return

  end if

  if ( r == 0.0D+00 ) then

    if ( p <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_POWER_FAST - Fatal error!'
      write ( *, '(a)' ) '  Base R is zero, and exponent is negative.'
      write ( *, '(a,i8)' ) '  Exponent P = ', p
      stop
    end if

    rp = 0.0D+00
    return

  end if
!
!  Special powers.
!
  if ( p == -1 ) then
    rp = 1.0D+00 / r
    mults = mults + 1
    return
  else if ( p == 0 ) then
    rp = 1.0D+00
    return
  else if ( p == 1 ) then
    rp = r
    return
  end if
!
!  Some work to do.
!
  p_mag = abs ( p )
  p_sign = sign ( i4_one, p )

  rp = 1.0D+00
  r2 = r

  do while ( 0 < p_mag )

    if ( mod ( p_mag, i4_two ) == 1 ) then
      rp = rp * r2
      mults = mults + 1
    end if

    p_mag = p_mag / 2
    r2 = r2 * r2
    mults = mults + 1

  end do

  if ( p_sign == -1 ) then
    rp = 1.0D+00 / rp
    mults = mults + 1
  end if

  return
end
function r8_pythag ( a, b )

!*****************************************************************************80
!
!! R8_PYTHAG computes sqrt ( A**2 + B**2 ), avoiding overflow and underflow.
!
!  Modified:
!
!    17 April 2004
!
!  Reference:
!
!    The SLATEC library
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the values for which sqrt ( A**2 + B**2 )
!    is desired.
!
!    Output, real ( kind = 8 ) R8_PYTHAG, the value of sqrt ( A**2 + B**2 ).
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) a_abs
  real    ( kind = 8 ) b
  real    ( kind = 8 ) b_abs
  real    ( kind = 8 ) r8_pythag

  a_abs = abs ( a )
  b_abs = abs ( b )

  if ( b_abs < a_abs ) then
    r8_pythag = a_abs * sqrt ( 1.0D+00 + ( b_abs / a_abs ) * ( b_abs / a_abs ) )
  else if ( b_abs == 0.0D+00 ) then
    r8_pythag = 0.0D+00
  else if ( a_abs <= b_abs ) then
    r8_pythag = b_abs * sqrt ( 1.0D+00 + ( a_abs / b_abs ) * ( a_abs / b_abs ) )
  end if

  return
end
subroutine r8_round2 ( nplace, x, xround )

!*****************************************************************************80
!
!! R8_ROUND2 rounds an R8 to a specified number of binary digits.
!
!  Discussion:
!
!    Assume that the input quantity X has the form
!
!      X = S * J * 2**L
!
!    where S is plus or minus 1, L is an integer, and J is a binary
!    mantissa which is either exactly zero, or greater than or equal
!    to 0.5 and strictly less than 1.0.
!
!    Then on return, XROUND will satisfy
!
!      XROUND = S * K * 2**L
!
!    where S and L are unchanged, and K is a binary mantissa which
!    agrees with J in the first NPLACE binary digits and is zero
!    thereafter.
!
!    If NPLACE is 0, XROUND will always be zero.
!
!    If NPLACE is 1, the mantissa of XROUND will be 0 or 0.5.
!
!    If NPLACE is 2, the mantissa of XROUND will be 0, 0.25, 0.50,
!    or 0.75.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPLACE, the number of binary digits to
!    preserve.  NPLACE should be 0 or positive.
!
!    Input, real ( kind = 8 ) X, the number to be decomposed.
!
!    Output, real ( kind = 8 ) XROUND, the rounded value of X.
!
  implicit none

  integer ( kind = 4 ) iplace
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nplace
  integer ( kind = 4 ) s
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xmant
  real    ( kind = 8 ) xround
  real    ( kind = 8 ) xtemp

  xround = 0.0D+00
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0D+00 ) then
    return
  end if

  if ( nplace <= 0 ) then
    return
  end if
!
!  2: Determine the sign S.
!
  if ( 0.0D+00 < x ) then
    s = 1
    xtemp = x
  else
    s = -1
    xtemp = -x
  end if
!
!  3: Force XTEMP to lie between 1 and 2, and compute the
!  logarithm L.
!
  l = 0

  do while ( 2.0D+00 <= xtemp )
    xtemp = xtemp / 2.0D+00
    l = l + 1
  end do

  do while ( xtemp < 1.0D+00 )
    xtemp = xtemp * 2.0D+00
    l = l - 1
  end do
!
!  4: Strip out the digits of the mantissa as XMANT, and decrease L.
!
  xmant = 0.0D+00
  iplace = 0

  do

    xmant = 2.0D+00 * xmant

    if ( 1.0D+00 <= xtemp ) then
      xmant = xmant + 1.0D+00
      xtemp = xtemp - 1.0D+00
    end if

    iplace = iplace + 1

    if ( xtemp == 0.0D+00 .or. nplace <= iplace ) then
      xround = s * xmant * 2.0D+00**l
      exit
    end if

    l = l - 1
    xtemp = xtemp * 2.0D+00

  end do

  return
end
subroutine r8_roundb ( base, nplace, x, xround )

!*****************************************************************************80
!
!! R8_ROUNDB rounds an R8 to a given number of digits in a given base.
!
!  Discussion:
!
!    The code does not seem to do a good job of rounding when
!    the base is negative!
!
!    Assume that the input quantity X has the form
!
!      X = S * J * BASE**L
!
!    where S is plus or minus 1, L is an integer, and J is a
!    mantissa base BASE which is either exactly zero, or greater
!    than or equal to (1/BASE) and less than 1.0.
!
!    Then on return, XROUND will satisfy
!
!      XROUND = S * K * BASE**L
!
!    where S and L are unchanged, and K is a mantissa base BASE
!    which agrees with J in the first NPLACE digits and is zero
!    thereafter.
!
!    Note that because of rounding, for most bases, most numbers
!    with a fractional quantities cannot be stored exactly in the
!    computer, and hence will have trailing "bogus" digits.
!
!    If NPLACE is 0, XROUND will always be zero.
!
!    If NPLACE is 1, the mantissa of XROUND will be 0,
!    1/BASE, 2/BASE, ..., (BASE-1)/BASE.
!
!    If NPLACE is 2, the mantissa of XROUND will be 0,
!    BASE/BASE**2, (BASE+1)/BASE**2, ...,
!    BASE**2-2/BASE**2, BASE**2-1/BASE**2.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BASE, the base of the arithmetic.
!    BASE must not be zero.  Theoretically, BASE may be negative.
!
!    Input, integer ( kind = 4 ) NPLACE, the number of digits base BASE to
!    preserve.  NPLACE should be 0 or positive.
!
!    Input, real ( kind = 8 ) X, the number to be decomposed.
!
!    Output, real ( kind = 8 ) XROUND, the rounded value of X.
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) iplace
  integer ( kind = 4 ) is
  integer ( kind = 4 ) js
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nplace
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xmant
  real    ( kind = 8 ) xround
  real    ( kind = 8 ) xtemp

  xround = 0.0D+00
!
!  0: Error checks.
!
  if ( base == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_ROUNDB - Fatal error!'
    write ( *, '(a)' ) '  The base BASE cannot be zero.'
    stop
  end if
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0D+00 ) then
    return
  end if

  if ( nplace <= 0 ) then
    return
  end if
!
!  2: Determine the sign IS.
!
  if ( 0.0D+00 < x ) then
    is = 1
    xtemp = x
  else
    is = -1
    xtemp = -x
  end if
!
!  3: Force XTEMP to lie between 1 and ABS(BASE), and compute the
!  logarithm L.
!
  l = 0

  do while ( abs ( base ) <= abs ( xtemp ) )

    xtemp = xtemp / real ( base, kind = 8 )

    if ( xtemp < 0.0D+00 ) then
      is = -is
      xtemp = -xtemp
    end if

    l = l + 1

  end do

  do while ( abs ( xtemp ) < 1.0D+00 )

    xtemp = xtemp * base

    if ( xtemp < 0.0D+00 ) then
      is = -is
      xtemp = -xtemp
    end if

    l = l - 1

  end do
!
!  4: Now strip out the digits of the mantissa as XMANT, and
!  decrease L.
!
  xmant = 0.0D+00
  iplace = 0
  js = is

  do

    xmant = base * xmant

    if ( xmant < 0.0D+00 ) then
      js = -js
      xmant = -xmant
    end if

    if ( 1.0D+00 <= xtemp ) then
      xmant = xmant + int ( xtemp )
      xtemp = xtemp - int ( xtemp )
    end if

    iplace = iplace + 1

    if ( xtemp == 0.0D+00 .or. nplace <= iplace ) then
      xround = js * xmant * real ( base, kind = 8 )**l
      exit
    end if

    l = l - 1
    xtemp = xtemp * base

    if ( xtemp < 0.0D+00 ) then
      is = -is
      xtemp = -xtemp
    end if

  end do

  return
end
subroutine r8_roundx ( nplace, x, xround )

!*****************************************************************************80
!
!! R8_ROUNDX rounds an R8.
!
!  Discussion:
!
!    Assume that the input quantity X has the form
!
!      X = S * J * 10**L
!
!    where S is plus or minus 1, L is an integer, and J is a decimal
!    mantissa which is either exactly zero, or greater than or equal
!    to 0.1 and less than 1.0.
!
!    Then on return, XROUND will satisfy
!
!      XROUND = S * K * 10**L
!
!    where S and L are unchanged, and K is a decimal mantissa which
!    agrees with J in the first NPLACE decimal digits and is zero
!    thereafter.
!
!    Note that because of rounding, most decimal fraction quantities
!    cannot be stored exactly in the computer, and hence will have
!    trailing "bogus" digits.
!
!    If NPLACE is 0, XROUND will always be zero.
!
!    If NPLACE is 1, the mantissa of XROUND will be 0, 0.1,
!    0.2, ..., or 0.9.
!
!    If NPLACE is 2, the mantissa of XROUND will be 0, 0.01, 0.02,
!    0.03, ..., 0.98, 0.99.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPLACE, the number of decimal digits to
!    preserve.  NPLACE should be 0 or positive.
!
!    Input, real ( kind = 8 ) X, the number to be decomposed.
!
!    Output, real ( kind = 8 ) XROUND, the rounded value of X.
!
  implicit none

  integer ( kind = 4 ) iplace
  integer ( kind = 4 ) is
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nplace
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xmant
  real    ( kind = 8 ) xround
  real    ( kind = 8 ) xtemp

  xround = 0.0D+00
!
!  1: Handle the special case of 0.
!
  if ( x == 0.0D+00 ) then
    return
  end if

  if ( nplace <= 0 ) then
    return
  end if
!
!  2: Determine the sign IS.
!
  if ( 0.0D+00 < x ) then
    is = 1
    xtemp = x
  else
    is = -1
    xtemp = -x
  end if
!
!  3: Force XTEMP to lie between 1 and 10, and compute the
!  logarithm L.
!
  l = 0

  do while ( 10.0D+00 <= x )
    xtemp = xtemp / 10.0D+00
    l = l + 1
  end do

  do while ( xtemp < 1.0D+00 )
    xtemp = xtemp * 10.0D+00
    l = l - 1
  end do
!
!  4: Now strip out the digits of the mantissa as XMANT, and
!  decrease L.
!
  xmant = 0.0D+00
  iplace = 0

  do

    xmant = 10.0D+00 * xmant

    if ( 1.0D+00 <= xtemp ) then
      xmant = xmant + int ( xtemp )
      xtemp = xtemp - int ( xtemp )
    end if

    iplace = iplace + 1

    if ( xtemp == 0.0D+00 .or. nplace <= iplace ) then
      xround = is * xmant * ( 10.0D+00**l )
      exit
    end if

    l = l - 1
    xtemp = xtemp * 10.0D+00

  end do

  return
end
function r8_sign ( x )

!*****************************************************************************80
!
!! R8_SIGN returns the sign of an R8.
!
!  Discussion:
!
!    value = -1 if X < 0;
!    value =  0 if X => 0.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose sign is desired.
!
!    Output, real ( kind = 8 ) R8_SIGN, the sign of X:
!
  implicit none

  real    ( kind = 8 ) r8_sign
  real    ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    r8_sign = -1.0D+00
  else
    r8_sign = +1.0D+00
  end if

  return
end
function r8_csqrt ( x )

!*****************************************************************************80
!
!! R8_CSQRT returns the complex square root of an R8.
!
!  Modified:
!
!    23 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose square root is desired.
!
!    Output, complex ( kind = 8 ) R8_CSQRT, the square root of X:
!
  implicit none

  real    ( kind = 8 ) argument
  real    ( kind = 8 ) magnitude
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  complex ( kind = 8 ) r8_csqrt
  real    ( kind = 8 ) x

  if ( 0.0D+00 < x ) then
    magnitude = x
    argument = 0.0D+00
  else if ( 0.0D+00 == x ) then
    magnitude = 0.0D+00
    argument = 0.0D+00
  else if ( x < 0.0D+00 ) then
    magnitude = -x
    argument = pi
  end if

  magnitude = sqrt ( magnitude )
  argument = argument / 2.0D+00

  r8_csqrt = magnitude * cmplx ( cos ( argument ), sin ( argument ) )

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Modified:
!
!    22 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
  real    ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine r8_swap3 ( x, y, z )

!*****************************************************************************80
!
!! R8_SWAP3 swaps three R8's.
!
!  Example:
!
!    Input:
!
!      X = 1, Y = 2, Z = 3
!
!    Output:
!
!      X = 2, Y = 3, Z = 1
!
!  Modified:
!
!    08 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y, Z, three values to be swapped.
!
  implicit none

  real    ( kind = 8 ) w
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
  real    ( kind = 8 ) z

  w = x
  x = y
  y = z
  z = w

  return
end
function r8_tiny ( )

!*****************************************************************************80
!
!! R8_TINY returns the smallest positive R8.
!
!  Discussion:
!
!    FORTRAN90 provides a built-in routine TINY ( X ) that
!    is more suitable for this purpose.
!
!  Modified:
!
!    08 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_TINY, a "tiny" value.
!
  implicit none

  real    ( kind = 8 ) r8_tiny

  r8_tiny = 0.4450147717014D-307

  return
end
subroutine r8_to_r8_discrete ( r, rmin, rmax, nr, rd )

!*****************************************************************************80
!
!! R8_TO_R8_DISCRETE maps R to RD in [RMIN, RMAX] with NR possible values.
!
!  Formula:
!
!    if ( R < RMIN ) then
!      RD = RMIN
!    else if ( RMAX < R ) then
!      RD = RMAX
!    else
!      T = nint ( ( NR - 1 ) * ( R - RMIN ) / ( RMAX - RMIN ) )
!      RD = RMIN + T * ( RMAX - RMIN ) / real ( NR - 1 )
!
!    In the special case where NR = 1, when
!
!      XD = 0.5 * ( RMAX + RMIN )
!
!  Modified:
!
!    21 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the number to be converted.
!
!    Input, real ( kind = 8 ) RMAX, RMIN, the maximum and minimum
!    values for RD.
!
!    Input, integer ( kind = 4 ) NR, the number of allowed values for XD.
!    NR should be at least 1.
!
!    Output, real ( kind = 8 ) RD, the corresponding discrete value.
!
  implicit none

  integer ( kind = 4 ) f
  integer ( kind = 4 ) nr
  real    ( kind = 8 ) r
  real    ( kind = 8 ) rd
  real    ( kind = 8 ) rmax
  real    ( kind = 8 ) rmin
!
!  Check for errors.
!
  if ( nr < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_TO_R8_DISCRETE - Fatal error!'
    write ( *, '(a,i8)' ) '  NR = ', nr
    write ( *, '(a)' ) '  but NR must be at least 1.'
    stop
  end if

  if ( nr == 1 ) then
    rd = 0.5D+00 * ( rmin + rmax )
    return
  end if

  if ( rmax == rmin ) then
    rd = rmax
    return
  end if

  f = nint ( real ( nr, kind = 8 ) * ( rmax - r ) / ( rmax - rmin ) )
  f = max ( f, 0 )
  f = min ( f, nr )

  rd = ( real (      f, kind = 8 ) * rmin   &
       + real ( nr - f, kind = 8 ) * rmax ) &
       / real ( nr,     kind = 8 )

  return
end
subroutine r8_to_dhms ( r, d, h, m, s )

!*****************************************************************************80
!
!! R8_TO_DHMS converts decimal days into days, hours, minutes, seconds.
!
!  Modified:
!
!    08 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, a decimal number representing a time 
!    period measured in days.
!
!    Output, integer ( kind = 4 ) D, H, M, S, the equivalent number of days, 
!    hours, minutes and seconds.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  real    ( kind = 8 ) r
  real    ( kind = 8 ) r_copy
  integer ( kind = 4 ) s

  r_copy = abs ( r )

  d = int ( r_copy )

  r_copy = r_copy - d
  r_copy = 24.0D+00 * r_copy
  h = int ( r_copy )

  r_copy = r_copy - h
  r_copy = 60.0D+00 * r_copy
  m = int ( r_copy )

  r_copy = r_copy - m
  r_copy = 60.0D+00 * r_copy
  s = int ( r_copy )

  if ( r < 0.0D+00 ) then
    d = -d
    h = -h
    m = -m
    s = -s
  end if

  return
end
subroutine r8_to_i4 ( x, xmin, xmax, ixmin, ixmax, ix )

!*****************************************************************************80
!
!! R8_TO_I4 maps X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].
!
!  Formula:
!
!    IX := IXMIN + ( IXMAX - IXMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
!    IX := min ( IX, max ( IXMIN, IXMAX ) )
!    IX := max ( IX, min ( IXMIN, IXMAX ) )
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be converted.
!
!    Input, real ( kind = 8 ) XMIN, XMAX, the range.  XMAX and
!    XMIN must not be equal.  It is not necessary that XMIN be less than XMAX.
!
!    Input, integer ( kind = 4 ) IXMIN, IXMAX, the allowed range of the output
!    variable.  IXMAX corresponds to XMAX, and IXMIN to XMIN.
!    It is not necessary that IXMIN be less than IXMAX.
!
!    Output, integer ( kind = 4 ) IX, the value in the range [IXMIN,IXMAX] that
!    corresponds to X.
!
  implicit none

  integer ( kind = 4 ) ix
  integer ( kind = 4 ) ixmax
  integer ( kind = 4 ) ixmin
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xmax
  real    ( kind = 8 ) xmin

  if ( xmax == xmin ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_TO_I4 - Fatal error!'
    write ( *, '(a)' ) '  XMAX = XMIN, making a zero divisor.'
    write ( *, '(a,g14.6)' ) '  XMAX = ', xmax
    write ( *, '(a,g14.6)' ) '  XMIN = ', xmin
    stop
  end if

  temp = &
      ( ( xmax - x        ) * real ( ixmin, kind = 8 )  &
      + (        x - xmin ) * real ( ixmax, kind = 8 ) ) &
      / ( xmax     - xmin )

  if ( 0.0D+00 <= temp ) then
    temp = temp + 0.5D+00
  else
    temp = temp - 0.5D+00
  end if

  ix = int ( temp )

  return
end
function r8_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R8_UNIFORM returns a scaled pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM, a number strictly between A and B.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  real    ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if

  r8_uniform = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8_unswap3 ( x, y, z )

!*****************************************************************************80
!
!! R8_UNSWAP3 unswaps three R8's.
!
!  Example:
!
!    Input:
!
!      X = 2, Y = 3, Z = 1
!
!    Output:
!
!      X = 1, Y = 2, Z = 3
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y, Z, three values to be swapped.
!
  implicit none

  real    ( kind = 8 ) w
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
  real    ( kind = 8 ) z

  w = z
  z = y
  y = x
  x = w

  return
end
function r8_walsh_1d ( x, digit )

!*****************************************************************************80
!
!! R8_WALSH_1D evaluates the Walsh function.
!
!  Discussion:
!
!    Consider the binary representation of X, and number the digits
!    in descending order, from leading to lowest, with the units digit
!    being numbered 0.
!
!    The Walsh function W(J)(X) is equal to the J-th binary digit of X.
!
!  Modified:
!
!    17 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Walsh function.
!
!    Input, integer ( kind = 4 ) DIGIT, the index of the Walsh function.
!
!    Output, real ( kind = 8 ) R8_WALSH_1D, the value of the Walsh function.
!
  implicit none

  integer ( kind = 4 ) digit
  integer ( kind = 4 ), parameter :: i4_two = 2
  integer ( kind = 4 ) n
  real    ( kind = 8 ) r8_walsh_1d
  real    ( kind = 8 ) x
  real    ( kind = 8 ) x_copy
!
!  Hide the effect of the sign of X.
!
  x_copy = abs ( x )
!
!  If DIGIT is positive, divide by 2 DIGIT times.
!  If DIGIT is negative, multiply by 2 (-DIGIT) times.
!
  x_copy = x_copy / 2.0D+00**digit
!
!  Make it an integer.
!  Because it's positive, and we're using INT, we don't change the
!  units digit.
!
  n = int ( x_copy )
!
!  Is the units digit odd or even?
!
  if ( mod ( n, i4_two ) == 0 ) then
    r8_walsh_1d = 0.0D+00
  else
    r8_walsh_1d = 1.0D+00
  end if

  return
end
subroutine r82_cheby ( n, alo, ahi, a )

!*****************************************************************************80
!
!! R82_CHEBY sets up the Chebyshev abscissas in an R8 interval.
!
!  Discussion:
!
!    The routine sets up a vector of X values spaced between the values
!    XLO and XHI in a similar way to the spacing of the Chebyshev
!    points of the same order in the interval [-1,1].
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points to compute.
!
!    Input, real ( kind = 8 ) ALO, AHI, the range.
!
!    Output, real ( kind = 8 ) A(N), the computed X values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) ahi
  real    ( kind = 8 ) alo
  real    ( kind = 8 ) arg
  integer ( kind = 4 ) i
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else if ( 1 < n ) then

    do i = 1, n

      arg = real ( 2 * i - 1, kind = 8 ) * pi &
          / real ( 2 * n, kind = 8 )

      a(i) = 0.5D+00 * ( ( 1.0D+00 + cos ( arg ) ) * alo &
                       + ( 1.0D+00 - cos ( arg ) ) * ahi )

    end do

  end if

  return
end
function r82_dist_l2 ( a1, a2 )

!*****************************************************************************80
!
!! R82_DIST_L2 returns the L2 distance between a pair of R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    The vector L2 norm is defined as:
!
!      sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), the vectors.
!
!    Output, real ( kind = 8 ) R82_DIST_L2, the L2 norm of the distance
!    between A1 and A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a1(dim_num)
  real    ( kind = 8 ) a2(dim_num)
  real    ( kind = 8 ) r82_dist_l2

  r82_dist_l2 = sqrt ( sum ( ( a1(1:dim_num) - a2(1:dim_num) )**2 ) )

  return
end
function r82_eq ( a1, a2 )

!*****************************************************************************80
!
!! R82_EQ == ( A1 == A2 ) for two R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    The comparison is lexicographic.
!
!    A1 == A2  <=>  A1(1) == A2(1) and A1(2) == A2(2).
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), two R82 vectors to be compared.
!
!    Output, logical R82_EQ, is TRUE if and only if A1 == A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a1(dim_num)
  real    ( kind = 8 ) a2(dim_num)
  logical r82_eq

  if ( all ( a1(1:dim_num) == a2(1:dim_num) ) ) then
    r82_eq = .true.
  else
    r82_eq = .false.
  end if

  return
end
function r82_ge ( a1, a2 )

!*****************************************************************************80
!
!! R82_GE == ( A1 >= A2 ) for two R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    The comparison is lexicographic.
!
!    A1 >= A2  <=>  A1(1) > A2(1) or ( A1(1) == A2(1) and A1(2) >= A2(2) ).
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), R82 vectors to be compared.
!
!    Output, logical R92_GE, is TRUE if and only if A1 >= A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a1(dim_num)
  real    ( kind = 8 ) a2(dim_num)
  integer ( kind = 4 ) i
  logical r82_ge

  r82_ge = .true.

  do i = 1, dim_num

    if ( a2(i) < a1(i) ) then
      r82_ge = .true.
      exit
    else if ( a1(i) < a2(i) ) then
      r82_ge = .false.
      exit
    end if

  end do

  return
end
function r82_gt ( a1, a2 )

!*****************************************************************************80
!
!! R82_GT == ( A1 > A2 ) for two R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R2, with two entries.
!
!    The comparison is lexicographic.
!
!    A1 > A2  <=>  A1(1) > A2(1) or ( A1(1) == A2(1) and A1(2) > A2(2) ).
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), R82 vectors to be compared.
!
!    Output, logical R82_GT, is TRUE if and only if A1 > A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a1(dim_num)
  real    ( kind = 8 ) a2(dim_num)
  integer ( kind = 4 ) i
  logical r82_gt

  r82_gt = .false.

  do i = 1, dim_num

    if ( a2(i) < a1(i) ) then
      r82_gt = .true.
      exit
    else if ( a1(i) < a2(i) ) then
      r82_gt = .false.
      exit
    end if

  end do

  return
end
function r82_le ( a1, a2 )

!*****************************************************************************80
!
!! R82_LE == ( A1 <= A2 ) for two R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    The comparison is lexicographic.
!
!    A1 <= A2  <=>  A1(1) < A2(1) or ( A1(1) == A2(1) and A1(2) <= A2(2) ).
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), R82 vectors to be compared.
!
!    Output, logical R82_LE, is TRUE if and only if A1 <= A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a1(dim_num)
  real    ( kind = 8 ) a2(dim_num)
  integer ( kind = 4 ) i
  logical r82_le

  r82_le = .true.

  do i = 1, dim_num

    if ( a1(i) < a2(i) ) then
      r82_le = .true.
      exit
    else if ( a2(i) < a1(i) ) then
      r82_le = .false.
      exit
    end if

  end do

  return
end
function r82_lt ( a1, a2 )

!*****************************************************************************80
!
!! R82_LT == ( A1 < A2 ) for two R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    The comparison is lexicographic.
!
!    A1 < A2  <=>  A1(1) < A2(1) or ( A1(1) == A2(1) and A1(2) < A2(2) ).
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), R82 vectors to be compared.
!
!    Output, logical R82_LT, is TRUE if and only if A1 < A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a1(dim_num)
  real    ( kind = 8 ) a2(dim_num)
  integer ( kind = 4 ) i
  logical r82_lt

  r82_lt = .false.

  do i = 1, dim_num

    if ( a1(i) < a2(i) ) then
      r82_lt = .true.
      exit
    else if ( a2(i) < a1(i) ) then
      r82_lt = .false.
      exit
    end if

  end do

  return
end
function r82_ne ( a1, a2 )

!*****************************************************************************80
!
!! R82_NE == ( A1 /= A2 ) for two R82's.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    The comparison is lexicographic.
!
!    A1 /= A2  <=>  A1(1) /= A2(1) or A1(2) /= A2(2).
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), R82 vectors to be compared.
!
!    Output, logical R82_NE, is TRUE if and only if A1 /= A2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a1(dim_num)
  real    ( kind = 8 ) a2(dim_num)
  logical r82_ne

  if ( any ( a1(1:dim_num) /= a2(1:dim_num) ) ) then
    r82_ne = .true.
  else
    r82_ne = .false.
  end if

  return
end
subroutine r82_print ( a, title )

!*****************************************************************************80
!
!! R82_PRINT prints an R82.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!    A format is used which suggests a coordinate pair:
!
!  Example:
!
!    Center : ( 1.23, 7.45 )
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2), the coordinates of the vector.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  real    ( kind = 8 ) a(2)
  character ( len = * ) title

  if ( len_trim ( title ) == 0 ) then
    write ( *, '( 2x, a1, g14.6, a1, g14.6, a1 )' ) '(', a(1), ',', a(2), ')'
  else
    write ( *, '( 2x, a, a4, g14.6, a1, g14.6, a1 )' ) &
      trim ( title ), ' : (', a(1), ',', a(2), ')'
  end if

  return
end
subroutine r82_swap ( x, y )

!*****************************************************************************80
!
!! R82_SWAP swaps two R82 values.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X(2), Y(2).  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) x(dim_num)
  real    ( kind = 8 ) y(dim_num)
  real    ( kind = 8 ) z(dim_num)

  z(1:dim_num) = x(1:dim_num)
  x(1:dim_num) = y(1:dim_num)
  y(1:dim_num) = z(1:dim_num)

  return
end
subroutine r82_uniform ( b, c, seed, a )

!*****************************************************************************80
!
!! R82_UNIFORM returns a random R82 value in a given range.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) B, C, the minimum and maximum values.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A(2), the randomly chosen value.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a(dim_num)
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed

  do i = 1, dim_num
    a(i) = r8_uniform ( b, c, seed )
  end do

  return
end
subroutine r82_unit_euclidean ( a )

!*****************************************************************************80
!
!! R82_UNIT_EUCLIDEAN_2D Euclidean normalizes an R82.
!
!  Discussion:
!
!    An R82 is a vector of type R8, with two entries.
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(2), the components of the vector.
!
  implicit none

  real    ( kind = 8 ) a(2)
  real    ( kind = 8 ) norm

  norm = sqrt ( a(1) * a(1) + a(2) * a(2) )

  if ( norm /= 0.0D+00 ) then
    a(1:2) = a(1:2) / norm
  end if

  return
end
subroutine r82poly2_print ( a, b, c, d, e, f )

!*****************************************************************************80
!
!! R82POLY2_PRINT prints a second order polynomial in two variables.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, E, F, the coefficients.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  real    ( kind = 8 ) e
  real    ( kind = 8 ) f

  write ( *, &
    '( 2x, f8.4, '' * x^2 + '', f8.4, '' * y^2 + '', f8.4, '' * xy  + '' )' ) &
    a, b, c

  write ( *, &
    '( 2x, f8.4, '' * x + '', f8.4, '' * y + '', f8.4, '' = 0 '' )' ) d, e, f

  return
end
subroutine r82poly2_type ( a, b, c, d, e, f, type )

!*****************************************************************************80
!
!! R82POLY2_TYPE analyzes a second order polynomial in two variables.
!
!  Discussion:
!
!    The polynomial has the form
!
!      A x^2 + B y^2 + C xy + Dx + Ey + F = 0
!
!    The possible types of the solution set are:
!
!     1: a hyperbola;
!        9x^2 -  4y^2       -36x - 24y -  36 = 0
!     2: a parabola;
!        4x^2 +  1y^2 - 4xy + 3x -  4y +   1 = 0;
!     3: an ellipse;
!        9x^2 + 16y^2       +36x - 32y -  92 = 0;
!     4: an aimaginary ellipse (no real solutions);
!         x^2 +   y^2       - 6x - 10y + 115 = 0;
!     5: a pair of intersecting lines;
!                        xy + 3x -   y -   3 = 0
!     6: one point;
!         x^2 +  2y^2       - 2x + 16y +  33 = 0;
!     7: a pair of distinct parallel lines;
!                 y^2            -  6y +   8 = 0
!     8: a pair of aimaginary parallel lines (no real solutions);
!                 y^2            -  6y +  10 = 0
!     9: a pair of coincident lines.
!                 y^2            -  2y +   1 = 0
!    10: a single line;
!                             2x -   y +   1 = 0;
!    11; all space;
!                                          0 = 0;
!    12; no solutions;
!                                          1 = 0;
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    CRC Press, 30th Edition, 1996, pages 282-284.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, E, F, the coefficients.
!
!    Output, integer ( kind = 4 ) TYPE, indicates the type of the solution set.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  real    ( kind = 8 ) delta
  real    ( kind = 8 ) e
  real    ( kind = 8 ) f
  real    ( kind = 8 ) j
  real    ( kind = 8 ) k
  integer ( kind = 4 ) type
!
!  Handle the degenerate case.
!
  if ( a == 0.0D+00 .and. &
       b == 0.0D+00 .and. &
       c == 0.0D+00 ) then
    if ( d == 0.0D+00 .and. e == 0.0D+00 ) then
      if ( f == 0.0D+00 ) then
        type = 11
      else
        type = 12
      end if
    else
      type = 10
    end if
    return
  end if

  delta = &
      8.0D+00 * a * b * f &
    + 2.0D+00 * c * e * d &
    - 2.0D+00 * a * e * e &
    - 2.0D+00 * b * d * d &
    - 2.0D+00 * f * c * c

  j = 4.0D+00 * a * b - c * c

  if ( delta /= 0.0D+00 ) then
    if ( j < 0.0D+00 ) then
      type = 1
    else if ( j == 0.0D+00 ) then
      type = 2
    else if ( 0.0D+00 < j ) then
      if ( sign ( 1.0D+00, delta ) /= sign ( 1.0D+00, ( a + b ) ) ) then
        type = 3
      else if ( sign ( 1.0D+00, delta ) == sign ( 1.0D+00, ( a + b ) ) ) then
        type = 4
      end if
    end if
  else if ( delta == 0.0D+00 ) then
    if ( j < 0.0D+00 ) then
      type = 5
    else if ( 0.0D+00 < j ) then
      type = 6
    else if ( j == 0.0D+00 ) then

      k = 4.0D+00 * ( a + b ) * f - d * d - e * e

      if ( k < 0.0D+00 ) then
        type = 7
      else if ( 0.0D+00 < k ) then
        type = 8
      else if ( k == 0.0D+00 ) then
        type = 9
      end if

    end if
  end if

  return
end
subroutine r82poly2_type_print ( type )

!*****************************************************************************80
!
!! R82POLY2_TYPE_PRINT prints the meaning of the output from R82POLY2_TYPE.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TYPE, the type index returned by R82POLY2_TYPE.
!
  implicit none

  integer ( kind = 4 ) type

  if ( type == 1 ) then
    write ( *, '(a)' ) '  The set of solutions forms a hyperbola.'
  else if ( type == 2 ) then
    write ( *, '(a)' ) '  The set of solutions forms a parabola.'
  else if ( type == 3 ) then
    write ( *, '(a)' ) '  The set of solutions forms an ellipse.'
  else if ( type == 4 ) then
    write ( *, '(a)' ) '  The set of solutions forms an aimaginary ellipse.'
    write ( *, '(a)' ) '  (There are no real solutions).'
  else if ( type == 5 ) then
    write ( *, '(a)' ) '  The set of solutions forms a pair of intersecting lines.'
  else if ( type == 6 ) then
    write ( *, '(a)' ) '  The set of solutions is a single point.'
  else if ( type == 7 ) then
    write ( *, '(a)' ) &
      '  The set of solutions form a pair of distinct parallel lines.'
  else if ( type == 8 ) then
    write ( *, '(a)' ) &
      '  The set of solutions forms a pair of aimaginary parallel lines.'
    write ( *, '(a)' ) '  (There are no real solutions).'
  else if ( type == 9 ) then
    write ( *, '(a)' ) '  The set of solutions forms a pair of coincident lines.'
  else if ( type == 10 ) then
    write ( *, '(a)' ) '  The set of solutions forms a single line.'
  else if ( type == 11 ) then
    write ( *, '(a)' ) '  The set of solutions is all space.'
  else if ( type == 12 ) then
    write ( *, '(a)' ) '  The set of solutions is empty.'
  else
    write ( *, '(a)' ) '  This type index is unknown.'
  end if

  return
end
subroutine r82vec_max ( n, a, amax )

!*****************************************************************************80
!
!! R82VEC_MAX returns the maximum value in an R82VEC.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of R8's.
!
!  Modified:
!
!    21 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(2,N), the array.
!
!    Output, real ( kind = 8 ) AMAX(2); the largest entries in each row.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(2,n)
  real    ( kind = 8 ) amax(2)

  amax(1) = maxval ( a(1,1:n) )
  amax(2) = maxval ( a(2,1:n) )

  return
end
subroutine r82vec_min ( n, a, amin )

!*****************************************************************************80
!
!! R82VEC_MIN returns the minimum value in an R82VEC.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of R82's.
!
!  Modified:
!
!    21 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(2,N), the array.
!
!    Output, real ( kind = 8 ) AMIN(2); the smallest entries in each row.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(2,n)
  real    ( kind = 8 ) amin(2)

  amin(1) = minval ( a(1,1:n) )
  amin(2) = minval ( a(2,1:n) )

  return
end
subroutine r82vec_order_type ( n, a, order )

!*****************************************************************************80
!
!! R82VEC_ORDER_TYPE finds the order type of an R82VEC.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of double precision real values.
!
!    The dictionary or lexicographic ordering is used.
!
!    (X1,Y1) < (X2,Y2)  <=>  X1 < X2 or ( X1 = X2 and Y1 < Y2).
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the array.
!
!    Input, real ( kind = 8 ) A(2,N), the array to be checked.
!
!    Output, integer ( kind = 4 ) ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a(dim_num,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
!
!  Search for the first value not equal to A(1,1).
!
  i = 1

  do

    i = i + 1

    if ( n < i ) then
      order = 0
      return
    end if

    if ( &
         a(1,1) <  a(1,i) .or. &
       ( a(1,1) == a(1,i) .and. a(2,1) < a(2,i) ) &
    ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( &
        a(1,i) <  a(1,1)  .or. &
      ( a(1,i) == a(1,1) .and. a(2,i) < a(2,1) ) &
    ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  Now we have a "direction".  Examine subsequent entries.
!
  do

    i = i + 1
    if ( n < i ) then
      exit
    end if

    if ( order == 1 ) then

      if ( &
          a(1,i) <  a(1,i-1) .or. &
        ( a(1,i) == a(1,i-1) .and. a(2,i) < a(2,i-1) ) &
      ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( &
          a(1,i) <  a(1,i-1) .or. &
        ( a(1,i) == a(1,i-1) .and. a(2,i) < a(2,i-1) ) &
      ) then
        order = -1
        exit
      else if ( &
         a(1,i) == a(1,i-1) .and. a(2,i) == a(2,i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( &
          a(1,i-1) <  a(1,i) .or. &
        ( a(1,i-1) == a(1,i) .and. a(2,i-1) < a(2,i) ) &
      ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( &
          a(1,i-1) <  a(1,i) .or. &
        ( a(1,i-1) == a(1,i) .and. a(2,i-1) < a(2,i) ) &
      ) then
        order = -1
        exit
      else if ( a(1,i) == a(1,i-1) .and. a(2,i) == a(2,i-1) ) then
        order = 3
      end if

    end if

  end do

  return
end
subroutine r82vec_part_quick_a ( n, a, l, r )

!*****************************************************************************80
!
!! R82VEC_PART_QUICK_A reorders an R82VEC as part of a quick sort.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of R82 values.
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
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input/output, real ( kind = 8 ) A(2,N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer ( kind = 4 ) L, R, the indices of A that define the three
!    segments.  Let KEY = the input value of A(1:2,1).  Then
!    I <= L                 A(1:2,I) < KEY;
!         L < I < R         A(1:2,I) = KEY;
!                 R <= I    KEY < A(1:2,I).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a(dim_num,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) key(dim_num)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r
  logical r8vec_eq
  logical r8vec_gt
  logical r8vec_lt

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R82VEC_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i8)' ) '  N = ', n
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key(1:dim_num) = a(1:dim_num,1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( r8vec_gt ( dim_num, a(1:dim_num,l+1), key(1:dim_num) ) ) then
      r = r - 1
      call r8vec_swap ( dim_num, a(1:dim_num,r), a(1:dim_num,l+1) )
    else if ( r8vec_eq ( dim_num, a(1:dim_num,l+1), key(1:dim_num) ) ) then
      m = m + 1
      call r8vec_swap ( dim_num, a(1:dim_num,m), a(1:dim_num,l+1) )
      l = l + 1
    else if ( r8vec_lt ( dim_num, a(1:dim_num,l+1), key(1:dim_num) ) ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(1:dim_num,i) = a(1:dim_num,i+m)
  end do

  l = l - m

  do i = 1, dim_num
    a(i,l+1:l+m) = key(i)
  end do

  return
end
subroutine r82vec_permute ( n, a, p )

!*****************************************************************************80
!
!! R82VEC_PERMUTE permutes an R82VEC in place.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of double precision real values.
!
!    The same logic can be used to permute an array of objects of any 
!    arithmetic type, or an array of objects of any complexity.  The only
!    temporary storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Modified:
!
!    13 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, real ( kind = 8 ) A(2,N), the array to be permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a(dim_num,n)
  real    ( kind = 8 ) a_temp(dim_num)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = -p(istart)
      cycle

    else

      a_temp(1:dim_num) = a(1:dim_num,istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = -p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(1:dim_num,iput) = a_temp(1:dim_num)
          exit
        end if

        a(1:dim_num,iput) = a(1:dim_num,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine r82vec_print ( n, a, title )

!*****************************************************************************80
!
!! R82VEC_PRINT prints an R82VEC.
!
!  Discussion:
!
!    An R82VEC is an array of pairs of R82's.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(2,N), the R82 vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a(dim_num,n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,(5g14.6))' ) i, a(1:dim_num,i)
  end do

  return
end
subroutine r82vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R82VEC_SORT_HEAP_INDEX_A ascending index heaps an R82VEC.
!
!  Discussion:
!
!    An R82VEC is an array of R82's.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(1:2,INDX(1:N)) is sorted,
!
!    or explicitly, by the call
!
!      call r82vec_permute ( n, a, indx )
!
!    after which A(1:2,I), I = 1 to N is sorted.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(2,N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(1:2,INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a(dim_num,n)
  real    ( kind = 8 ) aval(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval(1:dim_num) = a(1:dim_num,indxt)

    else

      indxt = indx(ir)
      aval(1:dim_num) = a(1:dim_num,indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if (   a(1,indx(j)) <  a(1,indx(j+1)) .or. &
             ( a(1,indx(j)) == a(1,indx(j+1)) .and. &
               a(2,indx(j)) <  a(2,indx(j+1)) ) ) then
          j = j + 1
        end if
      end if

      if (   aval(1) <  a(1,indx(j)) .or. &
           ( aval(1) == a(1,indx(j)) .and. &
             aval(2) <  a(2,indx(j)) ) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r82vec_sort_quick_a ( n, a )

!*****************************************************************************80
!
!! R82VEC_SORT_QUICK_A ascending sorts an R82VEC using quick sort.
!
!  Discussion:
!
!    An R82VEC is an array of R82's.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(2,N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: level_max = 25
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real    ( kind = 8 ) a(dim_num,n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) l_segment
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n_segment
  integer ( kind = 4 ) rsave(level_max)
  integer ( kind = 4 ) r_segment

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R82VEC_SORT_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i8)' ) '  N = ', n
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
    call r82vec_part_quick_a ( n_segment, a(1,base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( level_max < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R82VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i8)' ) '  Exceeding recursion maximum of ', level_max
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

        if ( 0 < n_segment ) then
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine r83_cross ( v1, v2, v3 )

!*****************************************************************************80
!
!! R83_CROSS computes the cross product of two R83's.
!
!  Discussion:
!
!    An R83 is a vector of type R8, with three entries.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), the the vectors.
!
!    Output, real ( kind = 8 ) V3(3), the cross product.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real    ( kind = 8 ) v1(dim_num)
  real    ( kind = 8 ) v2(dim_num)
  real    ( kind = 8 ) v3(dim_num)

  v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

  return
end
subroutine r83_print ( x, y, z, title )

!*****************************************************************************80
!
!! R83_PRINT prints an R83.
!
!  Discussion:
!
!    An R83 is a vector of type R8, with three entries.
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
!    Input, real ( kind = 8 ) X, Y, Z, the coordinates of the vector.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  character ( len = * ) title
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
  real    ( kind = 8 ) z

  if ( len_trim ( title ) == 0 ) then
    write ( *, '( 2x, a1, g14.6, a1, g14.6, a1, g14.6, a1 )' ) &
      '(', x, ',', y, ',', z, ')'
  else
    write ( *, '( 2x, a, a4, g14.6, a1, g14.6, a1, g14.6, a1 )' ) &
      trim ( title ), ' : (', x, ',', y, ',', z, ')'
  end if

  return
end
function r83_striple ( v1, v2, v3 )

!*****************************************************************************80
!
!! R83_STRIPLE computes the scalar triple product of three R83's.
!
!  Discussion:
!
!    An R83 is a vector of type double precision, with three entries.
!
!    STRIPLE = V1 dot ( V2 x V3 ).
!
!    STRIPLE is the volume of the parallelogram whose sides are
!    formed by V1, V2 and V3.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the three vectors.
!
!    Output, real ( kind = 8 ) R83_STRIPLE, the scalar triple product.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real    ( kind = 8 ) r83_striple
  real    ( kind = 8 ) v1(dim_num)
  real    ( kind = 8 ) v2(dim_num)
  real    ( kind = 8 ) v3(dim_num)

  r83_striple =  v1(1) * ( v2(2) * v3(3) - v2(3) * v3(2) ) &
               + v1(2) * ( v2(3) * v3(1) - v2(1) * v3(3) ) &
               + v1(3) * ( v2(1) * v3(2) - v2(2) * v3(1) )

  return
end
subroutine r83_swap ( x, y )

!*****************************************************************************80
!
!! R83_SWAP swaps two R83's.
!
!  Discussion:
!
!    An R83 is a vector of type R8, with three entries.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X(3), Y(3).  On output, the values 
!    of X and Y have been interchanged.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real    ( kind = 8 ) x(dim_num)
  real    ( kind = 8 ) y(dim_num)
  real    ( kind = 8 ) z(dim_num)

  z(1:dim_num) = x(1:dim_num)
  x(1:dim_num) = y(1:dim_num)
  y(1:dim_num) = z(1:dim_num)

  return
end
subroutine r83_unit_euclidean ( x, y, z )

!*****************************************************************************80
!
!! R83_UNIT_EUCLIDEAN Euclidean normalizes an R83.
!
!  Discussion:
!
!    An R83 is a vector of type R8, with three entries.
!
!  Modified:
!
!    17 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y, Z, the components of the vector.
!
  implicit none

  real    ( kind = 8 ) norm
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
  real    ( kind = 8 ) z

  norm = sqrt ( x * x + y * y + z * z )

  if ( norm /= 0.0D+00 ) then
    x = x / norm
    y = y / norm
    z = z / norm
  end if

  return
end
subroutine r83_vtriple ( v1, v2, v3, v )

!*****************************************************************************80
!
!! R83_VTRIPLE computes the vector triple product of three R83's.
!
!  Discussion:
!
!    An R83 is a vector of type R8, with three entries.
!
!    VTRIPLE = V1 x ( V2 x V3 )
!
!    VTRIPLE is a vector perpendicular to V1, lying in the plane
!    spanned by V2 and V3.  The norm of VTRIPLE is the product
!    of the norms of V1, V2 and V3.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), three vectors.
!
!    Output, real ( kind = 8 ) V(3), the vector triple product.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real    ( kind = 8 ) v(dim_num)
  real    ( kind = 8 ) v1(dim_num)
  real    ( kind = 8 ) v2(dim_num)
  real    ( kind = 8 ) v3(dim_num)
  real    ( kind = 8 ) v4(dim_num)

  call r83_cross ( v2, v3, v4 )

  call r83_cross ( v1, v4, v )

  return
end
subroutine r83vec_unit_l2 ( n, x )

!*****************************************************************************80
!
!! R83VEC_UNIT_L2 makes each R83 in R83VEC have unit L2 norm.
!
!  Discussion:
!
!    An R83VEC is a vector of R83's.
!
!    An R83 is a vector of three R8's.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of R83 vectors.
!
!    Input/output, real ( kind = 8 ) X(3,N), the coordinates of N R83 vectors.
!    On output, the nonzero vectors have been scaled to have unit L2 norm.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) i
  real    ( kind = 8 ) norm
  real    ( kind = 8 ) x(dim_num,n)

  do i = 1, n

    norm = sqrt ( sum ( x(1:dim_num,i)**2 ) )

    if ( norm /= 0.0D+00 ) then
      x(1:dim_num,i) = x(1:dim_num,i) / norm
    end if

  end do

  return
end
subroutine r84_unit_euclidean ( v )

!*****************************************************************************80
!
!! R84_UNIT_EUCLIDEAN_4D Euclidean normalizes an R84.
!
!  Discussion:
!
!    An R84 is a vector of four R8's.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) V(4), the components of the vector.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 4

  real    ( kind = 8 ) norm
  real    ( kind = 8 ) v(dim_num)

  norm = sqrt ( sum ( v(1:dim_num)**2 ) )

  if ( norm /= 0.0D+00 ) then
    v(1:dim_num) = v(1:dim_num) / norm
  end if

  return
end
subroutine r8block_expand_linear ( l, m, n, x, lfat, mfat, nfat, xfat )

!*****************************************************************************80
!
!! R8BLOCK_EXPAND_LINEAR linearly interpolates new data into an R8BLOCK.
!
!  Discussion:
!
!    An R8BLOCK is a 3D array of R8 values.
!
!    In this routine, the expansion is specified by giving the number
!    of intermediate values to generate between each pair of original
!    data rows and columns.
!
!    The interpolation is not actually linear.  It uses the functions
!
!      1, x, y, z, xy, xz, yz, xyz.
!
!  Modified:
!
!    17 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, M, N, the dimensions of the input data.
!
!    Input, real ( kind = 8 ) X(L,M,N), the original data.
!
!    Input, integer ( kind = 4 ) LFAT, MFAT, NFAT, the number of data values
!    to interpolate original data values in the first, second and third 
!    dimensions.
!
!    Output, real ( kind = 8 ) XFAT(L2,M2,N2), the fattened data, where
!    L2 = (L-1)*(LFAT+1)+1,
!    M2 = (M-1)*(MFAT+1)+1,
!    N2 = (N-1)*(NFAT+1)+1.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) lfat
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mfat
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nfat

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iii
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jjj
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) kkk
  integer ( kind = 4 ) kp1
  real    ( kind = 8 ) r
  real    ( kind = 8 ) s
  real    ( kind = 8 ) t
  real    ( kind = 8 ) x(l,m,n)
  real    ( kind = 8 ) x000
  real    ( kind = 8 ) x001
  real    ( kind = 8 ) x010
  real    ( kind = 8 ) x011
  real    ( kind = 8 ) x100
  real    ( kind = 8 ) x101
  real    ( kind = 8 ) x110
  real    ( kind = 8 ) x111
  real    ( kind = 8 ) xfat((l-1)*(lfat+1)+1,(m-1)*(mfat+1)+1,(n-1)*(nfat+1)+1)

  do i = 1, l

    if ( i < l ) then
      ihi = lfat
    else
      ihi = 0
    end if

    do j = 1, m

      if ( j < m ) then
        jhi = mfat
      else
        jhi = 0
      end if

      do k = 1, n

        if ( k < n ) then
          khi = nfat
        else
          khi = 0
        end if

        if ( i < l ) then
          ip1 = i + 1
        else
          ip1 = i
        end if

        if ( j < m ) then
          jp1 = j + 1
        else
          jp1 = j
        end if

        if ( k < n ) then
          kp1 = k + 1
        else
          kp1 = k
        end if

        x000 = x(i,j,k)
        x001 = x(i,j,kp1)
        x100 = x(ip1,j,k)
        x101 = x(ip1,j,kp1)
        x010 = x(i,jp1,k)
        x011 = x(i,jp1,kp1)
        x110 = x(ip1,jp1,k)
        x111 = x(ip1,jp1,kp1)

        do ii = 0, ihi

          r = real ( ii,      kind = 8 ) &
            / real ( ihi + 1, kind = 8 )

          do jj = 0, jhi

            s = real ( jj,      kind = 8 ) &
              / real ( jhi + 1, kind = 8 )

            do kk = 0, khi

              t = real ( kk,      kind = 8 ) &
                / real ( khi + 1, kind = 8 )

              iii = 1 + ( i - 1 ) * ( lfat + 1 ) + ii
              jjj = 1 + ( j - 1 ) * ( mfat + 1 ) + jj
              kkk = 1 + ( k - 1 ) * ( nfat + 1 ) + kk

              xfat(iii,jjj,kkk) = &
                  x000 * ( 1.0D+00 - r ) * ( 1.0D+00 - s ) * ( 1.0D+00 - t ) &
                + x001 * ( 1.0D+00 - r ) * ( 1.0D+00 - s ) * (           t ) &
                + x010 * ( 1.0D+00 - r ) * (           s ) * ( 1.0D+00 - t ) &
                + x011 * ( 1.0D+00 - r ) * (           s ) * (           t ) &
                + x100 * (           r ) * ( 1.0D+00 - s ) * ( 1.0D+00 - t ) &
                + x101 * (           r ) * ( 1.0D+00 - s ) * (           t ) &
                + x110 * (           r ) * (           s ) * ( 1.0D+00 - t ) &
                + x111 * (           r ) * (           s ) * (           t )

            end do

          end do

        end do

      end do

    end do

  end do

  return
end
subroutine r8block_print ( l, m, n, a, title )

!*****************************************************************************80
!
!! R8BLOCK_PRINT prints an R8BLOCK.
!
!  Discussion:
!
!    An R8BLOCK is a 3D array of R8 values.
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, M, N, the dimensions of the block.
!
!    Input, real ( kind = 8 ) A(L,M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(l,m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do k = 1, n

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) ' '

    do jlo = 1, m, 5
      jhi = min ( jlo + 4, m )
      write ( *, '(a)' ) ' '
      write ( *, '(10x,5(i8,6x))' ) (j, j = jlo, jhi )
      write ( *, '(a)' ) ' '
      do i = 1, l
        write ( *, '(2x,i8,5g14.6)' ) i, a(i,jlo:jhi,k)
      end do
    end do

  end do

  return
end
subroutine r8col_compare ( m, n, a, i, j, value )

!*****************************************************************************80
!
!! R8COL_COMPARE compares columns in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!    Output:
!
!      VALUE = -1
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) VALUE, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) value
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop
  end if

  value = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      value = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      value = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine r8col_find ( m, n, a, x, col )

!*****************************************************************************80
!
!! R8COL_FIND seeks a column value in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!  Example:
!
!    Input:
!
!      M = 3,
!      N = 4,
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!      x = ( 3.,
!            7.,
!           11. )
!
!    Output:
!
!      COL = 3
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), a table of numbers, regarded as
!    N columns of vectors of length M.
!
!    Input, real ( kind = 8 ) X(M), a vector to be matched with a column of A.
!
!    Output, integer ( kind = 4 ) COL, the index of the first column of A
!    which exactly matches every entry of X, or -1 if no match
!    could be found.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x(m)

  col = -1

  do j = 1, n

    col = j

    do i = 1, m
      if ( x(i) /= a(i,j) ) then
        col = -1
        exit
      end if
    end do

    if ( col /= -1 ) then
      return
    end if

  end do

  return
end
subroutine r8col_insert ( n_max, m, n, a, x, col )

!*****************************************************************************80
!
!! R8COL_INSERT inserts a column into an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!  Example:
!
!    Input:
!
!      N_MAX = 10,
!      M = 3,
!      N = 4,
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!      X = ( 3., 4., 18. )
!
!    Output:
!
!      N = 5,
!
!      A = (
!        1.  2.  3.  3.  4.
!        5.  6.  4.  7.  8.
!        9. 10. 18. 11. 12. )
!
!      COL = 3
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_MAX, the maximum number of columns in A.
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input/output, integer ( kind = 4 ) N, the number of columns.
!    If the new column is inserted into the table, then the output
!    value of N will be increased by 1.
!
!    Input/output, real ( kind = 8 ) A(M,N_MAX), a table of numbers, regarded
!    as an array of columns.  The columns must have been sorted
!    lexicographically.
!
!    Input, real ( kind = 8 ) X(M), a vector of data which will be inserted
!    into the table if it does not already occur.
!
!    Output, integer ( kind = 4 ) COL.
!    I, X was inserted into column I.
!    -I, column I was already equal to X.
!    0, N = N_MAX.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n_max

  real    ( kind = 8 ) a(m,n_max)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) high
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) n
  real    ( kind = 8 ) x(m)
!
!  Refuse to work if N_MAX <= N.
!
  if ( n_max <= n ) then
    col = 0
    return
  end if
!
!  Stick X temporarily in column N+1, just so it's easy to use R8COL_COMPARE.
!
  a(1:m,n+1) = x(1:m)
!
!  Do a binary search.
!
  low = 1
  high = n

  do

    if ( high < low ) then
      col = low
      exit
    end if

    mid = ( low + high ) / 2

    call r8col_compare ( m, n + 1, a, mid, n + 1, isgn )

    if ( isgn == 0 ) then
      col = -mid
      return
    else if ( isgn == -1 ) then
      low = mid + 1
    else if ( isgn == +1 ) then
      high = mid - 1
    end if

  end do
!
!  Shift part of the table up to make room.
!
  do j = n, col, -1
    a(1:m,j+1) = a(1:m,j)
  end do
!
!  Insert the new column.
!
  a(1:m,col) = x(1:m)

  n = n + 1

  return
end
subroutine r8col_max ( m, n, a, amax )

!*****************************************************************************80
!
!! R8COL_MAX returns the maximums in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!  Modified:
!
!    15 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) AMAX(N), the maximums of the columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) amax(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n

    amax(j) = a(1,j)
    do i = 2, m
      if ( amax(j) < a(i,j) ) then
        amax(j) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine r8col_max_index ( m, n, a, imax )

!*****************************************************************************80
!
!! R8COL_MAX_INDEX returns the indices of column maximums in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!  Modified:
!
!    15 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, integer ( kind = 4 ) IMAX(N); IMAX(I) is the row of A in which
!    the maximum for column I occurs.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) amax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax(n)
  integer ( kind = 4 ) j

  do j = 1, n

    imax(j) = 1
    amax = a(1,j)
    do i = 2, m
      if ( amax < a(i,j) ) then
        imax(j) = i
        amax = a(i,j)
      end if
    end do

  end do

  return
end
subroutine r8col_mean ( m, n, a, mean )

!*****************************************************************************80
!
!! R8COL_MEAN returns the column means of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!  Example:
!
!    A =
!      1  2  3
!      2  6  7
!
!    MEAN =
!      1.5  4.0  5.0
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) MEAN(N), the means, or averages, of the columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  real    ( kind = 8 ) mean(n)

  do j = 1, n
    mean(j) = sum ( a(1:m,j) )
  end do

  mean(1:n) = mean(1:n) / real ( m, kind = 8  )

  return
end
subroutine r8col_min ( m, n, a, amin )

!*****************************************************************************80
!
!! R8COL_MIN returns the column minimums of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) AMIN(N), the minimums of the columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) amin(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n

    amin(j) = a(1,j)
    do i = 2, m
      if ( a(i,j) < amin(j) ) then
        amin(j) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine r8col_min_index ( m, n, a, imin )

!*****************************************************************************80
!
!! R8COL_MIN_INDEX returns the indices of column minimums in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!  Modified:
!
!    15 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, integer ( kind = 4 ) IMIN(N); IMIN(I) is the row of A in which
!    the minimum for column I occurs.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) amin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imin(n)
  integer ( kind = 4 ) j

  do j = 1, n

    imin(j) = 1
    amin = a(1,j)
    do i = 2, m
      if ( a(i,j) < amin ) then
        imin(j) = i
        amin = a(i,j)
      end if
    end do

  end do

  return
end
subroutine r8col_part_quick_a ( m, n, a, l, r )

!*****************************************************************************80
!
!! R8COL_PART_QUICK_A reorders the columns of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!    The routine reorders the columns of A.  Using A(1:M,1) as a
!    key, all entries of A that are less than or equal to the key will
!    precede the key, which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      M = 2, N = 8
!      A = ( 2  8  6  0 10 10  0  5
!            4  8  2  2  6  0  6  8 )
!
!    Output:
!
!      L = 2, R = 4
!
!      A = (  0  0  2  8  6 10 10  4
!             2  6  4  8  2  6  0  8 )
!             ----     -------------
!             LEFT KEY     RIGHT
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the row dimension of A, and the length of
!    a column.
!
!    Input, integer ( kind = 4 ) N, the column dimension of A.
!
!    Input/output, real ( kind = 8 ) A(M,N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer ( kind = 4 ) L, R, the indices of A that define the three 
!    segments.  Let KEY = the input value of A(1:M,1).  Then
!    I <= L                 A(1:M,I) < KEY;
!         L < I < R         A(1:M,I) = KEY;
!                 R <= I    KEY < A(1:M,I).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) key(m)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) r
  logical r8vec_eq
  logical r8vec_gt
  logical r8vec_lt

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    return
  end if

  if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key(1:m) = a(1:m,1)
  k = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do j = 2, n

    if ( r8vec_gt ( m, a(1:m,l+1), key(1:m) ) ) then
      r = r - 1
      call r8vec_swap ( m, a(1:m,r), a(1:m,l+1) )
    else if ( r8vec_eq ( m, a(1:m,l+1), key(1:m) ) ) then
      k = k + 1
      call r8vec_swap ( m, a(1:m,k), a(1:m,l+1) )
      l = l + 1
    else if ( r8vec_lt ( m, a(1:m,l+1), key(1:m) ) ) then
      l = l + 1
    end if

  end do
!
!  Shift small elements to the left.
!
  do j = 1, l - k
    a(1:m,j) = a(1:m,j+k)
  end do
!
!  Shift KEY elements to center.
!
  do j = l - k + 1, l
    a(1:m,j) = key(1:m)
  end do
!
!  Update L.
!
  l = l - k

  return
end
subroutine r8col_permute ( m, n, a, p )

!*****************************************************************************80
!
!! R8COL_PERMUTE permutes an R8COL in place.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!    The same logic can be used to permute an array of objects of any 
!    arithmetic type, or an array of objects of any complexity.  The only
!    temporary storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      M = 2
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of objects.
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, real ( kind = 8 ) A(M,N), the array to be permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) a_temp(m)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = -p(istart)
      cycle

    else

      a_temp(1:m) = a(1:m,istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = -p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8COL_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(1:m,iput) = a_temp(1:m)
          exit
        end if

        a(1:m,iput) = a(1:m,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine r8col_sort_heap_a ( m, n, a )

!*****************************************************************************80
!
!! R8COL_SORT_HEAP_A ascending heapsorts an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
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
    if ( 0 < indx ) then

      call r8col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call r8col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8col_sort_quick_a ( m, n, a )

!*****************************************************************************80
!
!! R8COL_SORT_QUICK_A ascending quick sorts an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the row order of A, and the length of 
!    a column.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: level_max = 25
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) l_segment
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n_segment
  integer ( kind = 4 ) rsave(level_max)
  integer ( kind = 4 ) r_segment

  if ( m <= 0 ) then
    return
  end if

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_SORT_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    write ( *, '(a,i8)' ) '  N = ', n
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
    call r8col_part_quick_a ( m, n_segment, a(1,base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( level_max < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8COL_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i8)' ) '  Exceeding recursion maximum of ', level_max
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

        if ( 0 < n_segment ) then
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine r8col_sorted_unique ( m, n, a, unique_num )

!*****************************************************************************80
!
!! R8COL_SORTED_UNIQUE keeps unique elements in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, the sorted array of N columns of M-vectors.
!    On output, a sorted array of columns of M-vectors.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  j1 = 1

  do j2 = 2, n

    if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
      j1 = j1 + 1
      a(1:m,j1) = a(1:m,j2)
    end if

  end do

  unique_num = j1

  return
end
subroutine r8col_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! R8COL_SORTED_UNIQUE_COUNT counts unique elements in an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  j1 = 1

  do j2 = 2, n

    if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
      unique_num = unique_num + 1
      j1 = j2
    end if

  end do

  return
end
subroutine r8col_sortr_a ( m, n, a, key )

!*****************************************************************************80
!
!! R8COL_SORTR_A ascending sorts one column of an R8COL, adjusting all columns.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, an unsorted M by N array.
!    On output, rows of the array have been shifted in such
!    a way that column KEY of the array is in nondecreasing order.
!
!    Input, integer ( kind = 4 ) KEY, the column in which the "key" value
!    is stored.  On output, column KEY of the array will be
!    in nondecreasing order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) key

  if ( m <= 0 ) then
    return
  end if

  if ( key < 1 .or. n < key ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_SORTR_A - Fatal error!'
    write ( *, '(a)' ) '  The value of KEY is not a legal column index.'
    write ( *, '(a,i8)' ) '  KEY = ', key
    write ( *, '(a,i8)' ) '  N = ', n
    stop
  end if
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
    if ( 0 < indx ) then

      call r8row_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      if ( a(i,key) < a(j,key) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8col_sum ( m, n, a, colsum )

!*****************************************************************************80
!
!! R8COL_SUM sums the columns of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) COLSUM(N), the sums of the columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) colsum(n)
  integer ( kind = 4 ) j

  do j = 1, n
    colsum(j) = sum ( a(1:m,j) )
  end do

  return
end
subroutine r8col_swap ( m, n, a, j1, j2 )

!*****************************************************************************80
!
!! R8COL_SWAP swaps columns I and J of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, J1 = 2, J2 = 4
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!    Output:
!
!      A = (
!        1.  4.  3.  2.
!        5.  8.  7.  6.
!        9. 12. 11. 10. )
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Input, integer ( kind = 4 ) J1, J2, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) col(m)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2

  if ( j1 < 1 .or. n < j1 .or. j2 < 1 .or. n < j2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  J1 or J2 is out of bounds.'
    write ( *, '(a,i8)' ) '  J1 =    ', j1
    write ( *, '(a,i8)' ) '  J2 =    ', j2
    write ( *, '(a,i8)' ) '  NCOL = ', n
    stop
  end if

  if ( j1 == j2 ) then
    return
  end if

  col(1:m) = a(1:m,j1)
  a(1:m,j1) = a(1:m,j2)
  a(1:m,j2) = col(1:m)

  return
end
subroutine r8col_to_r8vec ( m, n, a, x )

!*****************************************************************************80
!
!! R8COL_TO_R8VEC converts an R8COL to an R8VEC.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!    An R8VEC is a vector of double precision values.
!
!  Example:
!
!    M = 3, N = 4
!
!    A =
!      11 12 13 14
!      21 22 23 24
!      31 32 33 34
!
!    X = ( 11, 21, 31, 12, 22, 32, 13, 23, 33, 14, 24, 34 )
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Output, real ( kind = 8 ) X(M*N), a vector containing the N columns of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) x(m*n)

  k = 1
  do j = 1, n
    x(k:k+m-1) = a(1:m,j)
    k = k + m
  end do

  return
end
subroutine r8col_variance ( m, n, a, variance )

!*****************************************************************************80
!
!! R8COL_VARIANCE returns the variances of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of double precision values, regarded
!    as an array of N columns of length M.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in 
!    the array.
!
!    Input, real ( kind = 8 ) A(M,N), the array whose variances are desired.
!
!    Output, real ( kind = 8 ) VARIANCE(N), the variances of the rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) mean
  real    ( kind = 8 ) variance(n)

  do j = 1, n

    mean = sum ( a(1:m,j) ) / real ( m, kind = 8  )

    variance(j) = 0.0D+00
    do i = 1, m
      variance(j) = variance(j) + ( a(i,j) - mean )**2
    end do

    if ( 1 < m ) then
      variance(j) = variance(j) / real ( m - 1, kind = 8 )
    else
      variance(j) = 0.0D+00
    end if

  end do

  return
end
function r8r8_compare ( x1, y1, x2, y2 )

!*****************************************************************************80
!
!! R8R8_COMPARE compares two R8R8's.
!
!  Discussion:
!
!    An R8R8 is simply a pair of double precision values, stored separately.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, the first vector.
!
!    Input, real ( kind = 8 ) X2, Y2, the second vector.
!
!    Output, integer ( kind = 4 ) R8R8_COMPARE: 
!    -1, (X1,Y1) < (X2,Y2);
!     0, (X1,Y1) = (X2,Y2);
!    +1, (X1,Y1) > (X2,Y2).
!
  implicit none

  integer ( kind = 4 ) compare
  integer ( kind = 4 ) r8r8_compare
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) x2
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2

  if ( x1 < x2 ) then
    compare = -1
  else if ( x2 < x1 ) then
    compare = +1
  else if ( y1 < y2 ) then
    compare = -1
  else if ( y2 < y1 ) then
    compare = +1
  else
    compare = 0
  end if

  r8r8_compare = compare

  return
end
subroutine r8r8_print ( a1, a2, title )

!*****************************************************************************80
!
!! R8R8_PRINT prints an R8R8.
!
!  Discussion:
!
!    An R8R8 is simply a pair of R8R8's, stored separately.
!
!    A format is used which suggests a coordinate pair:
!
!  Example:
!
!    Center : ( 1.23, 7.45 )
!
!  Modified:
!
!    29 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1, A2, the coordinates of the vector.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  real    ( kind = 8 ) a1
  real    ( kind = 8 ) a2
  character ( len = * ) title

  if ( len_trim ( title ) == 0 ) then
    write ( *, '( 2x, a1, g14.6, a1, g14.6, a1 )' ) '(', a1, ',', a2, ')'
  else
    write ( *, '( 2x, a, a4, g14.6, a1, g14.6, a1 )' ) &
      trim ( title ), ' : (', a1, ',', a2, ')'
  end if

  return
end
function r8r8r8_compare ( x1, y1, z1, x2, y2, z2 )

!*****************************************************************************80
!
!! R8R8R8_COMPARE compares two R8R8R8's.
!
!  Discussion:
!
!    An R8R8R8 is simply 3 double precision values, stored as scalars.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, Z1, the first vector.
!
!    Input, real ( kind = 8 ) X2, Y2, Z2, the second vector.
!
!    Output, integer ( kind = 4 ) R8R8R8_COMPARE: 
!    -1, (X1,Y1,Z1) < (X2,Y2,Z2);
!     0, (X1,Y1,Z1) = (X2,Y2,Z2);
!    +1, (X1,Y1,Z1) > (X2,Y2,Z2).
!
  implicit none

  integer ( kind = 4 ) compare
  integer ( kind = 4 ) r8r8r8_compare
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) x2
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2
  real    ( kind = 8 ) z1
  real    ( kind = 8 ) z2

  if ( x1 < x2 ) then
    compare = -1
  else if ( x2 < x1 ) then
    compare = +1
  else if ( y1 < y2 ) then
    compare = -1
  else if ( y2 < y1 ) then
    compare = +1
  else if ( z1 < z2 ) then
    compare = -1
  else if ( z2 < z1 ) then
    compare = +1
  else
    compare = 0
  end if

  r8r8r8_compare = compare

  return
end
subroutine r8r8r8vec_index_insert_unique ( n_max, n, x, y, z, indx, &
  xval, yval, zval, ival, ierror )

!*****************************************************************************80
!
!! R8R8R8VEC_INDEX_INSERT_UNIQUE inserts a unique R8R8R8in an indexed sorted list.
!
!  Discussion:
!
!    An R8R8R8VEC is set of N R8R8R8 items.
!
!    An R8R8R8 is simply 3 double precision values, stored as scalars.
!
!    If the input value does not occur in the current list, it is added,
!    and N, X, Y, Z and INDX are updated.
!
!  Modified:
!
!    06 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_MAX, the maximum size of the list.
!
!    Input/output, integer ( kind = 4 ) N, the size of the list.
!
!    Input/output, real ( kind = 8 ) X(N), Y(N), Z(N), the R8R8R8 vector.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, YVAL, ZVAL, the value to be inserted 
!    if it is not already in the list.
!
!    Output, integer ( kind = 4 ) IVAL, the index in X, Y, Z corresponding 
!    to the value XVAL, YVAL, ZVAL.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 if an error 
!    occurred.
!
  implicit none

  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  real    ( kind = 8 ) x(n_max)
  real    ( kind = 8 ) xval
  real    ( kind = 8 ) y(n_max)
  real    ( kind = 8 ) yval
  real    ( kind = 8 ) z(n_max)
  real    ( kind = 8 ) zval

  ierror = 0

  if ( n <= 0 ) then

    if ( n_max <= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    n = 1
    x(1) = xval
    y(1) = yval
    z(1) = zval
    indx(1) = 1
    ival = 1
    return

  end if
!
!  Does ( XVAL, YVAL, ZVAL ) already occur in ( X, Y, Z)?
!
  call r8r8r8vec_index_search ( n, x, y, z, indx, xval, yval, zval, &
    less, equal, more )

  if ( equal == 0 ) then

    if ( n_max <= n ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    x(n+1) = xval
    y(n+1) = yval
    z(n+1) = zval
    ival = n + 1
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1

  else

    ival = indx(equal)

  end if

  return
end
subroutine r8r8r8vec_index_search ( n, x, y, z, indx, xval, yval, &
  zval, less, equal, more )

!*****************************************************************************80
!
!! R8R8R8VEC_INDEX_SEARCH searches for an R8R8R8 value in an indexed sorted list.
!
!  Discussion:
!
!    An R8R8R8VEC is set of N R8R8R8 items.
!
!    An R8R8R8 is simply 3 double precision values, stored as scalars.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the list.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, YVAL, ZVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) LESS, EQUAL, MORE, the indexes in INDX of the
!    entries of X that are just less than, equal to, and just greater
!    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
!    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
!    is the greatest entry of X, then MORE is N+1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) compare
  integer ( kind = 4 ) r8r8r8_compare
  integer ( kind = 4 ) equal
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) more
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo
  real    ( kind = 8 ) xmid
  real    ( kind = 8 ) xval
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) yhi
  real    ( kind = 8 ) ylo
  real    ( kind = 8 ) ymid
  real    ( kind = 8 ) yval
  real    ( kind = 8 ) z(n)
  real    ( kind = 8 ) zhi
  real    ( kind = 8 ) zlo
  real    ( kind = 8 ) zmid
  real    ( kind = 8 ) zval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n

  xlo = x(indx(lo))
  ylo = y(indx(lo))
  zlo = z(indx(lo))

  xhi = x(indx(hi))
  yhi = y(indx(hi))
  zhi = z(indx(hi))

  compare = r8r8r8_compare ( xval, yval, zval, xlo, ylo, zlo )

  if ( compare == -1 ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( compare == 0 ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  compare = r8r8r8_compare ( xval, yval, zval, xhi, yhi, zhi )

  if ( compare == 1 ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( compare == 0 ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))
    ymid = y(indx(mid))
    zmid = z(indx(mid))

    compare = r8r8r8_compare ( xval, yval, zval, xmid, ymid, zmid )

    if ( compare == 0 ) then
      equal = mid
      less = mid - 1
      more = mid + 1
      return
    else if ( compare == -1 ) then
      hi = mid
    else if ( compare == +1 ) then
      lo = mid
    end if

  end do

  return
end
subroutine r8r8vec_index_insert_unique ( n_max, n, x, y, indx, xval, yval, &
  ival, ierror )

!*****************************************************************************80
!
!! R8R8VEC_INDEX_INSERT_UNIQUE inserts a unique R8R8 in an indexed sorted list.
!
!  Discussion:
!
!    An R8R8VEC is set of N R8R8 items.
!
!    An R8R8 is simply 2 double precision values, stored as scalars.
!
!    If the input value does not occur in the current list, it is added,
!    and N, X, Y and INDX are updated.
!
!  Modified:
!
!    06 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_MAX, the maximum size of the list.
!
!    Input/output, integer ( kind = 4 ) N, the size of the list.
!
!    Input/output, real ( kind = 8 ) X(N), Y(N), the list of R8R8 vectors.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, YVAL, the value to be inserted if it is
!    not already in the list.
!
!    Output, integer ( kind = 4 ) IVAL, the index in X, Y corresponding to the
!    value XVAL, YVAL.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 if an 
!    error occurred.
!
  implicit none

  integer ( kind = 4 ) n_max

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(n_max)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  real    ( kind = 8 ) x(n_max)
  real    ( kind = 8 ) xval
  real    ( kind = 8 ) y(n_max)
  real    ( kind = 8 ) yval

  ierror = 0

  if ( n <= 0 ) then

    if ( n_max <= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    n = 1
    x(1) = xval
    y(1) = yval
    indx(1) = 1
    ival = 1
    return

  end if
!
!  Does ( XVAL, YVAL ) already occur in ( X, Y )?
!
  call r8r8vec_index_search ( n, x, y, indx, xval, yval, less, equal, more )

  if ( equal == 0 ) then

    if ( n_max <= n ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    x(n+1) = xval
    y(n+1) = yval
    ival = n + 1
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1

  else

    ival = indx(equal)

  end if

  return
end
subroutine r8r8vec_index_search ( n, x, y, indx, xval, yval, less, equal, &
  more )

!*****************************************************************************80
!
!! R8R8VEC_INDEX_SEARCH searches for an R8R8 in an indexed sorted list.
!
!  Discussion:
!
!    An R8R8VEC is set of N R8R8 items.
!
!    An R8R8 is simply 2 double precision values, stored as scalars.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, YVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) LESS, EQUAL, MORE, the indexes in INDX of the
!    entries of X that are just less than, equal to, and just greater
!    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
!    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
!    is the greatest entry of X, then MORE is N+1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) compare
  integer ( kind = 4 ) r8r8_compare
  integer ( kind = 4 ) equal
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) more
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo
  real    ( kind = 8 ) xmid
  real    ( kind = 8 ) xval
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) yhi
  real    ( kind = 8 ) ylo
  real    ( kind = 8 ) ymid
  real    ( kind = 8 ) yval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n

  xlo = x(indx(lo))
  ylo = y(indx(lo))

  xhi = x(indx(hi))
  yhi = y(indx(hi))

  compare = r8r8_compare ( xval, yval, xlo, ylo )

  if ( compare == -1 ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( compare == 0 ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  compare = r8r8_compare ( xval, yval, xhi, yhi )

  if ( compare == 1 ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( compare == 0 ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))
    ymid = y(indx(mid))

    compare = r8r8_compare ( xval, yval, xmid, ymid )

    if ( compare == 0 ) then
      equal = mid
      less = mid - 1
      more = mid + 1
      return
    else if ( compare == -1 ) then
      hi = mid
    else if ( compare == +1 ) then
      lo = mid
    end if

  end do

  return
end
subroutine r8int_to_r8int ( rmin, rmax, r, r2min, r2max, r2 )

!*****************************************************************************80
!
!! R8INT_TO_R8INT maps one R8INT to another.
!
!  Discussion:
!
!    The formula used is
!
!      R2 := R2MIN + ( R2MAX - R2MIN ) * ( R - RMIN ) / ( RMAX - RMIN )
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RMIN, RMAX, the first range.
!
!    Input, real ( kind = 8 ) R, the number to be converted.
!
!    Input, real ( kind = 8 ) R2MAX, R2MIN, the second range.
!
!    Output, real ( kind = 8 ) R2, the corresponding value in 
!    the range [R2MIN,R2MAX].
!
  implicit none

  real    ( kind = 8 ) r
  real    ( kind = 8 ) rmax
  real    ( kind = 8 ) rmin
  real    ( kind = 8 ) r2
  real    ( kind = 8 ) r2max
  real    ( kind = 8 ) r2min

  if ( rmax == rmin ) then

    r2 = ( r2max + r2min ) / 2.0D+00

  else

    r2 = ( ( ( rmax - r        ) * r2min   &
           + (        r - rmin ) * r2max ) &
           / ( rmax     - rmin ) )

  end if

  return
end
subroutine r8int_to_i4int ( rmin, rmax, r, imin, imax, i )

!*****************************************************************************80
!
!! R8INT_TO_I4INT maps an R8INT to an integer interval.
!
!  Discussion:
!
!    The formula used is
!
!      I := IMIN + ( IMAX - IMIN ) * ( R - RMIN ) / ( RMAX - RMIN )
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RMIN, RMAX, the range.
!
!    Input, real ( kind = 8 ) R, the number to be converted.
!
!    Input, integer ( kind = 4 ) IMAX, IMIN, the integer range.
!
!    Output, integer ( kind = 4 ) I, the corresponding value in the 
!    range [IMIN,IMAX].
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real    ( kind = 8 ) r
  real    ( kind = 8 ) rmax
  real    ( kind = 8 ) rmin

  if ( rmax == rmin ) then

    i = ( imax + imin ) / 2

  else

    i = nint ( &
      ( ( rmax - r        ) * real ( imin, kind = 8 )   &
      + (        r - rmin ) * real ( imax, kind = 8 ) ) &
      / ( rmax     - rmin ) )

  end if

  return
end
subroutine r8mat_cholesky_factor ( n, a, c, ierror )

!*****************************************************************************80
!
!! R8MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric matrix.
!
!  Discussion:
!
!    The matrix must be symmetric and positive semidefinite.
!
!    For a positive semidefinite symmetric matrix A, the Cholesky factorization
!    is a lower triangular matrix L such that:
!
!      A = L * L'
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of 
!    the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 8 ) C(N,N), the N by N lower triangular
!    Cholesky factor.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    1, warning, the matrix is positive semidefinite.  The factorization
!    was carried out, but the matrix is singular.
!    2, error, the matrix has at least one negative eigenvalue.  The
!    factorization could not be completed.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) c(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  real    ( kind = 8 ) sum2

  ierror = 0

  c(1:n,1:n) = a(1:n,1:n)

  do j = 1, n

    c(1:j-1,j) = 0.0D+00

    do i = j, n

      sum2 = c(j,i) - dot_product ( c(j,1:j-1), c(i,1:j-1) )

      if ( i == j ) then
        if ( sum2 < 0.0D+00 ) then
          ierror = 2
          return
        else if ( sum2 == 0.0D+00 ) then
          ierror = 1
        else
          c(i,j) = sqrt ( sum2 )
        end if
      else
        if ( c(j,j) /= 0.0D+00 ) then
          c(i,j) = sum2 / c(j,j)
        else
          c(i,j) = 0.0D+00
        end if
      end if

    end do

  end do

  return
end
subroutine r8mat_cholesky_solve ( n, a, b, x )

!*****************************************************************************80
!
!! R8MAT_CHOLESKY_SOLVE solves a Cholesky factored linear system A * x = b.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of 
!    the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N Cholesky factor of the
!    system matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) x(n)
!
!  Solve L * y = b.
!
  call r8mat_l_solve ( n, a, b, x )
!
!  Solve L' * x = y.
!
  call r8mat_lt_solve ( n, a, x, x )

  return
end
subroutine r8mat_det ( n, a, det )

!*****************************************************************************80
!
!! R8MAT_DET computes the determinant of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    07 December 2004
!
!  Reference:
!
!    Helmut Spaeth,
!    Cluster Analysis Algorithms
!    for Data Reduction and Classification of Objects,
!    Ellis Horwood, 1980, page 125-127.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n,n)
  real    ( kind = 8 ) det
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) piv(1)

  b(1:n,1:n) = a(1:n,1:n)

  det = 1.0D+00

  do k = 1, n

    piv = maxloc ( abs ( b(k:n,k) ) )

    m = piv(1) + k - 1

    if ( m /= k ) then
      det = -det
      call r8_swap ( b(m,k), b(k,k) )
    end if

    det = det * b(k,k)

    if ( b(k,k) /= 0.0D+00 ) then

      b(k+1:n,k) = -b(k+1:n,k) / b(k,k)

      do j = k + 1, n
        if ( m /= k ) then
          call r8_swap ( b(m,j), b(k,j) )
        end if
        b(k+1:n,j) = b(k+1:n,j) + b(k+1:n,k) * b(k,j)
      end do

    end if

  end do

  return
end
function r8mat_det_2d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_2D computes the determinant of a 2 by 2 R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    The formula for the determinant of a 2 by 2 matrix is
!
!      a11 * a22 - a12 * a21.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2,2), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) R8MAT_DET_2D, the determinant of the matrix.
!
  implicit none

  real    ( kind = 8 ) a(2,2)
  real    ( kind = 8 ) r8mat_det_2d

  r8mat_det_2d = a(1,1) * a(2,2) - a(1,2) * a(2,1)

  return
end
function r8mat_det_3d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_3D computes the determinant of a 3 by 3 R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    The formula for the determinant of a 3 by 3 matrix is
!
!        a11 * a22 * a33 - a11 * a23 * a32
!      + a12 * a23 * a31 - a12 * a21 * a33
!      + a13 * a21 * a32 - a13 * a22 * a31
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(3,3), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) R8MAT_DET_3D, the determinant of the matrix.
!
  implicit none

  real    ( kind = 8 ) a(3,3)
  real    ( kind = 8 ) r8mat_det_3d

  r8mat_det_3d = &
         a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
       + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
       + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )

  return
end
function r8mat_det_4d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(4,4), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) R8MAT_DET_4D, the determinant of the matrix.
!
  implicit none

  real    ( kind = 8 ) a(4,4)
  real    ( kind = 8 ) r8mat_det_4d

  r8mat_det_4d = &
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
function r8mat_det_5d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_5D computes the determinant of a 5 by 5 R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
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
!    Input, real ( kind = 8 ) A(5,5), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) R8MAT_DET_5D, the determinant of the matrix.
!
  implicit none

  real    ( kind = 8 ) a(5,5)
  real    ( kind = 8 ) b(4,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) r8mat_det_4d
  real    ( kind = 8 ) r8mat_det_5d
!
!  Expand the determinant into the sum of the determinants of the
!  five 4 by 4 matrices created by dropping row 1, and column k.
!
  r8mat_det_5d = 0.0D+00

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

    r8mat_det_5d = r8mat_det_5d + (-1)**( k + 1 ) * a(1,k) * r8mat_det_4d ( b )

  end do

  return
end
subroutine r8mat_diag_add_scalar ( n, a, s )

!*****************************************************************************80
!
!! R8MAT_DIAG_ADD_SCALAR adds a scalar to the diagonal of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N), the N by N matrix to be modified.
!
!    Input, real ( kind = 8 ) S, the value to be added to the diagonal
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) s

  do i = 1, n
    a(i,i) = a(i,i) + s
  end do

  return
end
subroutine r8mat_diag_add_vector ( n, a, v )

!*****************************************************************************80
!
!! R8MAT_DIAG_ADD_VECTOR adds a vector to the diagonal of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of 
!    the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Input, real ( kind = 8 ) V(N), the vector to be added to the diagonal of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) v(n)

  do i = 1, n
    a(i,i) = a(i,i) + v(i)
  end do

  return
end
subroutine r8mat_diag_get_vector ( n, a, v )

!*****************************************************************************80
!
!! R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 8 ) V(N), the diagonal entries
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) v(n)

  do i = 1, n
    v(i) = a(i,i)
  end do

  return
end
subroutine r8mat_diag_set_scalar ( n, a, s )

!*****************************************************************************80
!
!! R8MAT_DIAG_SET_SCALAR sets the diagonal of an R8MAT to a scalar value.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N), the N by N matrix to be modified.
!
!    Input, real ( kind = 8 ) S, the value to be assigned to the diagonal
!    of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) s

  do i = 1, n
    a(i,i) = s
  end do

  return
end
subroutine r8mat_diag_set_vector ( n, a, v )

!*****************************************************************************80
!
!! R8MAT_DIAG_SET_VECTOR sets the diagonal of an R8MAT to a vector.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Input, real ( kind = 8 ) V(N), the vector to be assigned to the
!    diagonal of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) v(n)

  do i = 1, n
    a(i,i) = v(i)
  end do

  return
end
subroutine r8mat_expand_linear ( m, n, x, mfat, nfat, xfat )

!*****************************************************************************80
!
!! R8MAT_EXPAND_LINEAR linearly interpolates new data into an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    In this routine, the expansion is specified by giving the number
!    of intermediate values to generate between each pair of original
!    data rows and columns.
!
!    The interpolation is not actually linear.  It uses the functions
!
!      1, x, y, and xy.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of 
!    input data.
!
!    Input, real ( kind = 8 ) X(M,N), the original data.
!
!    Input, integer ( kind = 4 ) MFAT, NFAT, the number of data values 
!    to interpolate between each row, and each column, of original data values.
!
!    Output, real ( kind = 8 ) XFAT(M2,N2), the fattened data, where
!    M2 = (M-1)*(MFAT+1)+1,
!    N2 = (N-1)*(NFAT+1)+1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mfat
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nfat

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iii
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jjj
  integer ( kind = 4 ) jp1
  real    ( kind = 8 ) s
  real    ( kind = 8 ) t
  real    ( kind = 8 ) x(m,n)
  real    ( kind = 8 ) x00
  real    ( kind = 8 ) x01
  real    ( kind = 8 ) x10
  real    ( kind = 8 ) x11
  real    ( kind = 8 ) xfat((m-1)*(mfat+1)+1,(n-1)*(nfat+1)+1)

  do i = 1, m

    if ( i < m ) then
      ihi = mfat
    else
      ihi = 0
    end if

    do j = 1, n

      if ( j < n ) then
        jhi = nfat
      else
        jhi = 0
      end if

      if ( i < m ) then
        ip1 = i + 1
      else
        ip1 = i
      end if

      if ( j < n ) then
        jp1 = j + 1
      else
        jp1 = j
      end if

      x00 = x(i,j)
      x10 = x(ip1,j)
      x01 = x(i,jp1)
      x11 = x(ip1,jp1)

      do ii = 0, ihi

        s = real ( ii, kind = 8 ) &
          / real ( ihi + 1, kind = 8 )

        do jj = 0, jhi

          t = real ( jj, kind = 8 ) &
            / real ( jhi + 1, kind = 8 )

          iii = 1 + ( i - 1 ) * ( mfat + 1 ) + ii
          jjj = 1 + ( j - 1 ) * ( nfat + 1 ) + jj

          xfat(iii,jjj) = &
                                            x00   &
              + s     * (       x10       - x00 ) &
              + t     * (             x01 - x00 ) &
              + s * t * ( x11 - x10 - x01 + x00 )

        end do

      end do

    end do

  end do

  return
end
subroutine r8mat_expand_linear2 ( m, n, a, m2, n2, a2 )

!*****************************************************************************80
!
!! R8MAT_EXPAND_LINEAR2 expands an R8MAT by linear interpolation.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    In this version of the routine, the expansion is indicated
!    by specifying the dimensions of the expanded array.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), a "small" M by N array.
!
!    Input, integer ( kind = 4 ) M2, N2, the number of rows and columns in A2.
!
!    Output, real ( kind = 8 ) A2(M2,N2), the expanded array, which
!    contains an interpolated version of the data in A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) a2(m2,n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real    ( kind = 8 ) r
  real    ( kind = 8 ) r1
  real    ( kind = 8 ) r2
  real    ( kind = 8 ) s
  real    ( kind = 8 ) s1
  real    ( kind = 8 ) s2

  do i = 1, m2

    if ( m2 == 1 ) then
      r = 0.5D+00
    else
      r = real ( i - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )
    end if

    i1 = 1 + int ( r * real ( m - 1, kind = 8 ) )
    i2 = i1 + 1

    if ( m < i2 ) then
      i1 = m - 1
      i2 = m
    end if

    r1 = real ( i1 - 1, kind = 8 ) &
       / real ( m - 1, kind = 8 )

    r2 = real ( i2 - 1, kind = 8 ) &
       / real ( m - 1, kind = 8 )

    do j = 1, n2

      if ( n2 == 1 ) then
        s = 0.5D+00
      else
        s = real ( j - 1, kind = 8 ) &
          / real ( n2 - 1, kind = 8 )
      end if

      j1 = 1 + int ( s * real ( n - 1, kind = 8 ) )
      j2 = j1 + 1

      if ( n < j2 ) then
        j1 = n - 1
        j2 = n
      end if

      s1 = real ( j1 - 1, kind = 8 ) &
         / real ( n - 1, kind = 8 )

      s2 = real ( j2 - 1, kind = 8 ) &
         / real ( n - 1, kind = 8 )

      a2(i,j) = &
        ( ( r2 - r ) * ( s2 - s ) * a(i1,j1) &
        + ( r - r1 ) * ( s2 - s ) * a(i2,j1) &
        + ( r2 - r ) * ( s - s1 ) * a(i1,j2) &
        + ( r - r1 ) * ( s - s1 ) * a(i2,j2) ) &
        / ( ( r2 - r1 ) * ( s2 - s1 ) )

    end do

  end do

  return
end
subroutine r8mat_givens_post ( n, a, row, col, g )

!*****************************************************************************80
!
!! R8MAT_GIVENS_POST computes the Givens postmultiplier rotation matrix.
!
!  Discussion:
!
!    The Givens post-multiplier matrix G(ROW,COL) has the property that
!    the (ROW,COL)-th entry of A*G is zero.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices A and G.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix to be operated upon.
!
!    Input, integer ( kind = 4 ) ROW, COL, the row and column of the
!    entry of A*G which is to be zeroed out.
!
!    Output, real ( kind = 8 ) G(N,N), the Givens rotation matrix.
!    G is an orthogonal matrix, that is, the inverse of
!    G is the transpose of G.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) col
  real    ( kind = 8 ) g(n,n)
  integer ( kind = 4 ) row
  real    ( kind = 8 ) theta

  call r8mat_identity ( n, g )

  theta = atan2 ( a(row,col), a(row,row) )

  g(row,row) =  cos ( theta )
  g(row,col) = -sin ( theta )
  g(col,row) =  sin ( theta )
  g(col,col) =  cos ( theta )

  return
end
subroutine r8mat_givens_pre ( n, a, row, col, g )

!*****************************************************************************80
!
!! R8MAT_GIVENS_PRE computes the Givens premultiplier rotation matrix.
!
!  Discussion:
!
!    The Givens premultiplier rotation matrix G(ROW,COL) has the
!    property that the (ROW,COL)-th entry of G*A is zero.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices A and G.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix to be operated upon.
!
!    Input, integer ( kind = 4 ) ROW, COL, the row and column of the
!    entry of the G*A which is to be zeroed out.
!
!    Output, real ( kind = 8 ) G(N,N), the Givens rotation matrix.
!    G is an orthogonal matrix, that is, the inverse of
!    G is the transpose of G.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) col
  real    ( kind = 8 ) g(n,n)
  integer ( kind = 4 ) row
  real    ( kind = 8 ) theta

  call r8mat_identity ( n, g )

  theta = atan2 ( a(row,col), a(col,col) )

  g(row,row) =  cos ( theta )
  g(row,col) = -sin ( theta )
  g(col,row) =  sin ( theta )
  g(col,col) =  cos ( theta )

  return
end
subroutine r8mat_hess ( fx, n, x, h )

!*****************************************************************************80
!
!! R8MAT_HESS approximates a Hessian matrix via finite differences.
!
!  Discussion:
!
!    H(I,J) = d2 F / d X(I) d X(J)
!
!    The values returned by this routine will be only approximate.
!    In some cases, they will be so poor that they are useless.
!    However, one of the best applications of this routine is for
!    checking your own Hessian calculations, since as Heraclitus
!    said, you'll never get the same result twice when you differentiate
!    a complicated expression by hand.
!
!    The user function routine, here called "FX", should have the form:
!
!      subroutine fx ( n, x, f )
!      integer n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FX, the name of the user function routine.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the values of the variables.
!
!    Output, real ( kind = 8 ) H(N,N), the approximated N by N Hessian matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) eps
  real    ( kind = 8 ) f00
  real    ( kind = 8 ) fmm
  real    ( kind = 8 ) fmp
  real    ( kind = 8 ) fpm
  real    ( kind = 8 ) fpp
  external fx
  real    ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) s(n)
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xi
  real    ( kind = 8 ) xj
!
!  Choose the stepsizes.
!
  eps = ( epsilon ( eps ) )**0.33D+00

  do i = 1, n
    s(i) = eps * max ( abs ( x(i) ), 1.0D+00 )
  end do
!
!  Calculate the diagonal elements.
!
  do i = 1, n

    xi = x(i)

    call fx ( n, x, f00 )

    x(i) = xi + s(i)
    call fx ( n, x, fpp )

    x(i) = xi - s(i)
    call fx ( n, x, fmm )

    h(i,i) = ( ( fpp - f00 ) + ( fmm - f00 ) ) / s(i)**2

    x(i) = xi

  end do
!
!  Calculate the off diagonal elements.
!
  do i = 1, n

    xi = x(i)

    do j = i + 1, n

      xj = x(j)

      x(i) = xi + s(i)
      x(j) = xj + s(j)
      call fx ( n, x, fpp )

      x(i) = xi + s(i)
      x(j) = xj - s(j)
      call fx ( n, x, fpm )

      x(i) = xi - s(i)
      x(j) = xj + s(j)
      call fx ( n, x, fmp )

      x(i) = xi - s(i)
      x(j) = xj - s(j)
      call fx ( n, x, fmm )

      h(j,i) = ( ( fpp - fpm ) + ( fmm - fmp ) ) / ( 4.0D+00 * s(i) * s(j) )

      h(i,j) = h(j,i)

      x(j) = xj

    end do

    x(i) = xi

  end do

  return
end
subroutine r8mat_house_axh ( n, a, v, ah )

!*****************************************************************************80
!
!! R8MAT_HOUSE_AXH computes A*H where H is a compact Householder matrix.
!
!  Discussion:
!
!    The Householder matrix H(V) is defined by
!
!      H(V) = I - 2 * v * v' / ( v' * v )
!
!    This routine is not particularly efficient.
!
!  Modified:
!
!    26 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix to be postmultiplied.
!
!    Input, real ( kind = 8 ) V(N), a vector defining a Householder matrix.
!
!    Output, real ( kind = 8 ) AH(N,N), the product A*H.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) ah(n,n)
  real    ( kind = 8 ) ah_temp(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) v(n)
  real    ( kind = 8 ) v_normsq

  v_normsq = sum ( v(1:n)**2 )
!
!  Compute A*H' = A*H
!
  do i = 1, n
    do j = 1, n
      ah_temp(i,j) = a(i,j)
      do k = 1, n
        ah_temp(i,j) = ah_temp(i,j) - 2.0D+00 * a(i,k) * v(k) * v(j) / v_normsq
      end do
    end do
  end do
!
!  Copy the temporary result into AH.
!  Doing it this way means the user can identify the input arguments A and AH.
!
  ah(1:n,1:n) = ah_temp(1:n,1:n)

  return
end
subroutine r8mat_house_form ( n, v, h )

!*****************************************************************************80
!
!! R8MAT_HOUSE_FORM constructs a Householder matrix from its compact form.
!
!  Discussion:
!
!    H(v) = I - 2 * v * v' / ( v' * v )
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) V(N), the vector defining the Householder matrix.
!
!    Output, real ( kind = 8 ) H(N,N), the Householder matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) beta
  real    ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) v(n)
!
!  Compute the L2 norm of V.
!
  beta = sum ( v(1:n)**2 )
!
!  Form the matrix H.
!
  call r8mat_identity ( n, h )

  do i = 1, n
    do j = 1, n
      h(i,j) = h(i,j) - 2.0D+00 * v(i) * v(j) / beta
    end do
  end do

  return
end
subroutine r8mat_house_hxa ( n, a, v, ha )

!*****************************************************************************80
!
!! R8MAT_HOUSE_HXA computes H*A where H is a compact Householder matrix.
!
!  Discussion:
!
!    The Householder matrix H(V) is defined by
!
!      H(V) = I - 2 * v * v' / ( v' * v )
!
!    This routine is not particularly efficient.
!
!  Modified:
!
!    26 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix to be premultiplied.
!
!    Input, real ( kind = 8 ) V(N), a vector defining a Householder matrix.
!
!    Output, real ( kind = 8 ) HA(N,N), the product H*A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) ha(n,n)
  real    ( kind = 8 ) ha_temp(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) v(n)
  real    ( kind = 8 ) v_normsq

  v_normsq = sum ( v(1:n)**2 )
!
!  Compute A*H' = A*H
!
  do i = 1, n
    do j = 1, n
      ha_temp(i,j) = a(i,j)
      do k = 1, n
        ha_temp(i,j) = ha_temp(i,j) - 2.0D+00 * v(i) * v(k) * a(k,j) / v_normsq
      end do
    end do
  end do
!
!  Copy the temporary result into HA.
!  Doing it this way means the user can identify the input arguments A and HA.
!
  ha(1:n,1:n) = ha_temp(1:n,1:n)

  return
end
subroutine r8mat_house_post ( n, a, row, col, h )

!*****************************************************************************80
!
!! R8MAT_HOUSE_POST computes a Householder post-multiplier matrix.
!
!  Discussion:
!
!    H(ROW,COL) has the property that the ROW-th column of
!    A*H(ROW,COL) is zero from entry COL+1 to the end.
!
!    In the most common case, where a QR factorization is being computed,
!    ROW = COL.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix whose Householder matrix 
!    is to be computed.
!
!    Input, integer ( kind = 4 ) ROW, COL, specify the location of the
!    entry of the matrix A which is to be preserved.  The entries in
!    the same row, but higher column, will be zeroed out if
!    A is postmultiplied by H.
!
!    Output, real ( kind = 8 ) H(N,N), the Householder matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) col
  real    ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) row
  real    ( kind = 8 ) v(n)
  real    ( kind = 8 ) w(n)
!
!  Set up the vector V.
!
  w(1:col-1) = 0.0D+00
  w(col:n) = a(row,col:n)

  call r8vec_house_column ( n, w, col, v )
!
!  Form the matrix H(V).
!
  call r8mat_house_form ( n, v, h )

  return
end
subroutine r8mat_house_pre ( n, a, row, col, h )

!*****************************************************************************80
!
!! R8MAT_HOUSE_PRE computes a Householder pre-multiplier matrix.
!
!  Discussion:
!
!    H(ROW,COL) has the property that the COL-th column of
!    H(ROW,COL)*A is zero from entry ROW+1 to the end.
!
!    In the most common case, where a QR factorization is being computed,
!    ROW = COL.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix whose Householder matrix
!    is to be computed.
!
!    Input, integer ( kind = 4 ) ROW, COL, specify the location of the
!    entry of the matrix A which is to be preserved.  The entries in
!    the same column, but higher rows, will be zeroed out if A is
!    premultiplied by H.
!
!    Output, real ( kind = 8 ) H(N,N), the Householder matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) col
  real    ( kind = 8 ) h(n,n)
  integer ( kind = 4 ) row
  real    ( kind = 8 ) v(n)
  real    ( kind = 8 ) w(n)
!
!  Set up the vector V.
!
  w(1:row-1) = 0.0D+00
  w(row:n) = a(row:n,col)

  call r8vec_house_column ( n, w, row, v )
!
!  Form the matrix H(V).
!
  call r8mat_house_form ( n, v, h )

  return
end
subroutine r8mat_identity ( n, a )

!*****************************************************************************80
!
!! R8MAT_IDENTITY stores the identity matrix in an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
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
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Output, real ( kind = 8 ) A(N,N), the N by N identity matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i

  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    a(i,i) = 1.0D+00
  end do

  return
end
function r8mat_in_01 ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_IN_01 is TRUE if the entries of an R8MAT are in the range [0,1].
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    06 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, logical R8MAT_IN_01, is TRUE if every entry of A is
!    between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  logical r8mat_in_01

  if ( any ( a(1:m,1:n) < 0.0D+00 .or. 1.0D+00 < a(1:m,1:n) ) ) then
    r8mat_in_01 = .false.
  else
    r8mat_in_01 = .true.
  end if

  return
end
subroutine r8mat_inverse_2d ( a, b, det )

!*****************************************************************************80
!
!! R8MAT_INVERSE_2D inverts a 2 by 2 R8MAT using Cramer's rule.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    If the determinant is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If the determinant is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Modified:
!
!    29 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2,2), the matrix to be inverted.
!
!    Output, real ( kind = 8 ) B(2,2), the inverse of the matrix A.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
!
  implicit none

  real    ( kind = 8 ) a(2,2)
  real    ( kind = 8 ) b(2,2)
  real    ( kind = 8 ) det
  real    ( kind = 8 ) r8mat_det_2d
!
!  Compute the determinant of A.
!
  det = r8mat_det_2d ( a )

  if ( det == 0.0D+00 ) then

    b(1:2,1:2) = 0.0D+00

  else

    b(1,1) =  a(2,2) / det
    b(1,2) = -a(1,2) / det
    b(2,1) = -a(2,1) / det
    b(2,2) =  a(1,1) / det

  end if

  return
end
subroutine r8mat_inverse_3d ( a, b, det )

!*****************************************************************************80
!
!! R8MAT_INVERSE_3D inverts a 3 by 3 R8MAT using Cramer's rule.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    If the determinant is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If the determinant is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Modified:
!
!    29 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(3,3), the matrix to be inverted.
!
!    Output, real ( kind = 8 ) B(3,3), the inverse of the matrix A.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
!
  implicit none

  real    ( kind = 8 ) a(3,3)
  real    ( kind = 8 ) b(3,3)
  real    ( kind = 8 ) det
  real    ( kind = 8 ) r8mat_det_3d
!
!  Compute the determinant of A.
!
  det = r8mat_det_3d ( a )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then
    b(1:3,1:3) = 0.0D+00
    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit
!  formula.
!
  b(1,1) =  ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
  b(1,2) = -( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) / det
  b(1,3) =  ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det

  b(2,1) = -( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) / det
  b(2,2) =  ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
  b(2,3) = -( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) / det

  b(3,1) =  ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
  b(3,2) = -( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) / det
  b(3,3) =  ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det

  return
end
subroutine r8mat_inverse_4d ( a, b, det )

!*****************************************************************************80
!
!! R8MAT_INVERSE_4D inverts a 4 by 4 R8MAT using Cramer's rule.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    If the determinant is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If the determinant is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Modified:
!
!    13 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(4,4), the matrix to be inverted.
!
!    Output, real ( kind = 8 ) B(4,4), the inverse of the matrix A.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
!
  implicit none

  real    ( kind = 8 ) a(4,4)
  real    ( kind = 8 ) b(4,4)
  real    ( kind = 8 ) det
  real    ( kind = 8 ) r8mat_det_4d
!
!  Compute the determinant of A.
!
  det = r8mat_det_4d ( a )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then

    b(1:4,1:4) = 0.0D+00

    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit formula.
!
  b(1,1) = +( &
        + a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(2,3) * ( a(3,4) * a(4,2) - a(3,2) * a(4,4) ) &
        + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        ) / det

  b(2,1) = -( &
        + a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(2,3) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
        ) / det

  b(3,1) = +( &
        + a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        + a(2,2) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(4,1) = -( &
        + a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        + a(2,2) * ( a(3,3) * a(4,1) - a(3,1) * a(4,3) ) &
        + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(1,2) = -( &
        + a(1,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(1,3) * ( a(3,4) * a(4,2) - a(3,2) * a(4,4) ) &
        + a(1,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        ) / det

  b(2,2) = +( &
        + a(1,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
        + a(1,3) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(1,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
        ) / det

  b(3,2) = -( &
        + a(1,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
        + a(1,2) * ( a(3,4) * a(4,1) - a(3,1) * a(4,4) ) &
        + a(1,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(4,2) = +( &
        + a(1,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
        + a(1,2) * ( a(3,3) * a(4,1) - a(3,1) * a(4,3) ) &
        + a(1,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) &
        ) / det

  b(1,3) = +( &
        + a(1,2) * ( a(2,3) * a(4,4) - a(2,4) * a(4,3) ) &
        + a(1,3) * ( a(2,4) * a(4,2) - a(2,2) * a(4,4) ) &
        + a(1,4) * ( a(2,2) * a(4,3) - a(2,3) * a(4,2) ) &
        ) / det

  b(2,3) = -( &
        + a(1,1) * ( a(2,3) * a(4,4) - a(2,4) * a(4,3) ) &
        + a(1,3) * ( a(2,4) * a(4,1) - a(2,1) * a(4,4) ) &
        + a(1,4) * ( a(2,1) * a(4,3) - a(2,3) * a(4,1) ) &
        ) / det

  b(3,3) = +( &
        + a(1,1) * ( a(2,2) * a(4,4) - a(2,4) * a(4,2) ) &
        + a(1,2) * ( a(2,4) * a(4,1) - a(2,1) * a(4,4) ) &
        + a(1,4) * ( a(2,1) * a(4,2) - a(2,2) * a(4,1) ) &
        ) / det

  b(4,3) = -( &
        + a(1,1) * ( a(2,2) * a(4,3) - a(2,3) * a(4,2) ) &
        + a(1,2) * ( a(2,3) * a(4,1) - a(2,1) * a(4,3) ) &
        + a(1,3) * ( a(2,1) * a(4,2) - a(2,2) * a(4,1) ) &
        ) / det

  b(1,4) = -( &
        + a(1,2) * ( a(2,3) * a(3,4) - a(2,4) * a(3,3) ) &
        + a(1,3) * ( a(2,4) * a(3,2) - a(2,2) * a(3,4) ) &
        + a(1,4) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
        ) / det

  b(2,4) = +( &
        + a(1,1) * ( a(2,3) * a(3,4) - a(2,4) * a(3,3) ) &
        + a(1,3) * ( a(2,4) * a(3,1) - a(2,1) * a(3,4) ) &
        + a(1,4) * ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) &
        ) / det

  b(3,4) = -( &
        + a(1,1) * ( a(2,2) * a(3,4) - a(2,4) * a(3,2) ) &
        + a(1,2) * ( a(2,4) * a(3,1) - a(2,1) * a(3,4) ) &
        + a(1,4) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) &
        ) / det

  b(4,4) = +( &
        + a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
        + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
        + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) &
        ) / det

  return
end
subroutine r8mat_jac ( m, n, eps, fx, x, fprime )

!*****************************************************************************80
!
!! R8MAT_JAC estimates a dense jacobian matrix of the function FX.
!
!  Discussion:
!
!    FPRIME(I,J) = d F(I) / d X(J).
!
!    The jacobian is assumed to be dense, and the LINPACK/LAPACK
!    double precision general matrix storage mode ("DGE") is used.
!
!    Forward differences are used, requiring N+1 function evaluations.
!
!    Values of EPS have typically been chosen between
!    sqrt ( EPSMCH ) and sqrt ( sqrt ( EPSMCH ) ) where EPSMCH is the
!    machine tolerance.
!
!    If EPS is too small, then F(X+EPS) will be the same as
!    F(X), and the jacobian will be full of zero entries.
!
!    If EPS is too large, the finite difference estimate will
!    be inaccurate.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) EPS, a tolerance to be used for shifting the
!    X values during the finite differencing.  No single value
!    of EPS will be reliable for all vectors X and functions FX.
!
!    Input, external FX, the name of the user written
!    routine which evaluates the function at a given point X.
!
!    FX should have the form:
!
!      subroutine fx ( m, n, x, f )
!      integer m
!      integer n
!      real ( kind = 8 ) f(m)
!      real ( kind = 8 ) x(n)
!      f(1:m) = ...
!      return
!      end
!
!    Input, real ( kind = 8 ) X(N), the point where the jacobian
!    is to be estimated.
!
!    Output, real ( kind = 8 ) FPRIME(M,N), the M by N estimated jacobian
!    matrix.
!

  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) del
  real    ( kind = 8 ) eps
  real    ( kind = 8 ) fprime(m,n)
  external fx
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xsave
  real    ( kind = 8 ) work1(m)
  real    ( kind = 8 ) work2(m)
!
!  Evaluate the function at the base point, X.
!
  call fx ( m, n, x, work2 )
!
!  Now, one by one, vary each component J of the base point X, and
!  estimate DF(I)/DX(J) = ( F(X+) - F(X) )/ DEL.
!
  do j = 1, n

    xsave = x(j)
    del = eps * ( 1.0D+00 + abs ( x(j) ) )
    x(j) = x(j) + del
    call fx ( m, n, x, work1 )
    x(j) = xsave
    fprime(1:m,j) = ( work1(1:m) - work2(1:m) ) / del

  end do

  return
end
subroutine r8mat_l_inverse ( n, a, b )

!*****************************************************************************80
!
!! R8MAT_L_INVERSE inverts a lower triangular R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    A lower triangular matrix is a matrix whose only nonzero entries
!    occur on or below the diagonal.
!
!    The inverse of a lower triangular matrix is a lower triangular matrix.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the lower triangular matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = 1, n

    do i = 1, n

      if ( i < j ) then
        b(i,j) = 0.0D+00
      else if ( j == i ) then
        b(i,j) = 1.0D+00 / a(i,j)
      else
        b(i,j) = -dot_product ( a(i,1:i-1), b(1:i-1,j) ) / a(i,i)
      end if

    end do
  end do

  return
end
subroutine r8mat_l_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_L_PRINT prints a lower triangular R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Example:
!
!    M = 5, N = 5
!    A = (/ 11, 21, 31, 41, 51, 22, 32, 42, 52, 33, 43, 53, 44, 54, 55 /)
!
!    11
!    21 22
!    31 32 33
!    41 42 43 44
!    51 52 53 54 55
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(*), the M by N matrix.  Only the lower
!    triangular elements are stored, in column major order.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(10)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) size
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  jmax = min ( n, m )

  if ( m <= n ) then
    size = ( m * ( m + 1 ) ) / 2
  else if ( n < m ) then
    size = ( n * ( n + 1 ) ) / 2 + ( m - n ) * n
  end if

  if ( all ( a(1:size) == aint ( a(1:size) ) ) ) then

    nn = 10

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(a8,10i8)' ) '  Col   ', ( j, j = jlo, jhi )
      write ( *, '(a6)' ) '  Row '
      do i = jlo, m
        jhi = min ( jlo + nn - 1, i, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2
        end do
        write ( *, '(i8,10i8)' ) i, int ( a(indx(1:jhi+1-jlo)) )
      end do
    end do

  else if ( maxval ( abs ( a(1:size) ) ) < 1000000.0D+00 ) then

    nn = 5

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m - 1, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(8x,5(i8,6x))' ) ( j, j = jlo, jhi )
      write ( *, '(a)' ) ' '
      do i = jlo, m
        jhi = min ( jlo + nn - 1, i, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2
        end do
        write ( *, '(i8,5f14.6)' ) i, a(indx(1:jhi+1-jlo))
      end do
    end do

  else

    nn = 5

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m - 1, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(8x,5(i8,6x))' ) ( j, j = jlo, jhi )
      write ( *, '(a)' ) ' '
      do i = jlo, m
        jhi = min ( jlo + nn - 1, i, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2
        end do
        write ( *, '(i8,5g14.6)' ) i, a(indx(1:jhi+1-jlo))
      end do
    end do

  end if

  return
end
subroutine r8mat_l_solve ( n, a, b, x )

!*****************************************************************************80
!
!! R8MAT_L_SOLVE solves a lower triangular linear system.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of 
!    the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N lower triangular matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) x(n)
!
!  Solve L * x = b.
!
  do i = 1, n
    x(i) = ( b(i) - dot_product ( a(i,1:i-1), x(1:i-1) ) ) / a(i,i)
  end do

  return
end
subroutine r8mat_l1_inverse ( n, a, b )

!*****************************************************************************80
!
!! R8MAT_L1_INVERSE inverts a double precision unit lower triangular R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    A unit lower triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's above the main diagonal.
!
!    The inverse of a unit lower triangular matrix is also
!    a unit lower triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call r8mat_l1_inverse ( n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the unit lower triangular matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n

    do j = 1, n

      if ( i < j ) then
        b(i,j) = 0.0D+00
      else if ( j == i ) then
        b(i,j) = 1.0D+00
      else
        b(i,j) = -dot_product ( a(i,1:i-1), b(1:i-1,j) )
      end if

    end do
  end do

  return
end
subroutine r8mat_lt_solve ( n, a, b, x )

!*****************************************************************************80
!
!! R8MAT_LT_SOLVE solves a transposed lower triangular linear system.
!
!  Discussion:
!
!    Given the lower triangular matrix A, the linear system to be solved is:
!
!      A' * x = b
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N lower triangular matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) x(n)
!
!  Solve L'*x = b.
!
  do i = n, 1, -1
    x(i) = ( b(i) - dot_product ( x(i+1:n), a(i+1:n,i) ) ) / a(i,i)
  end do

  return
end
subroutine r8mat_lu ( m, n, a, l, p, u )

!*****************************************************************************80
!
!! R8MAT_LU computes the LU factorization of a rectangular R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    The routine is given an M by N matrix A, and produces
!
!      L, an M by M unit lower triangular matrix,
!      U, an M by N upper triangular matrix, and
!      P, an M by M permutation matrix P,
!
!    so that
!
!      A = P' * L * U.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix to be factored.
!
!    Output, real ( kind = 8 ) L(M,M), the M by M unit lower triangular factor.
!
!    Output, real ( kind = 8 ) P(M,M), the M by M permutation matrix.
!
!    Output, real ( kind = 8 ) U(M,N), the M by N upper triangular factor.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipiv
  integer ( kind = 4 ) j
  real    ( kind = 8 ) l(m,m)
  real    ( kind = 8 ) p(m,m)
  real    ( kind = 8 ) pivot
  real    ( kind = 8 ) u(m,n)

!  Initialize:
!
!    U:=A
!    L:=Identity
!    P:=Identity
!
  u(1:m,1:n) = a(1:m,1:n)

  call r8mat_identity ( m, l )

  p(1:m,1:m) = l(1:m,1:m)
!
!  On step J, find the pivot row, IPIV, and the pivot value PIVOT.
!
  do j = 1, min ( m - 1, n )

    pivot = 0.0D+00
    ipiv = 0

    do i = j, m

      if ( pivot < abs ( u(i,j) ) ) then
        pivot = abs ( u(i,j) )
        ipiv = i
      end if

    end do
!
!  Unless IPIV is zero, swap rows J and IPIV.
!
    if ( ipiv /= 0 ) then

      call r8row_swap ( m, n, u, j, ipiv )

      call r8row_swap ( m, m, l, j, ipiv )

      call r8row_swap ( m, m, p, j, ipiv )
!
!  Zero out the entries in column J, from row J+1 to M.
!
      do i = j + 1, m

        if ( u(i,j) /= 0.0D+00 ) then

          l(i,j) = u(i,j) / u(j,j)

          u(i,j) = 0.0D+00

          u(i,j+1:n) = u(i,j+1:n) - l(i,j) * u(j,j+1:n)

        end if

      end do

    end if

  end do

  return
end
function r8mat_max ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_MAX returns the maximum entry of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix.
!
!    Output, real ( kind = 8 ) R8MAT_MAX, the maximum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) r8mat_max

  r8mat_max = maxval ( a(1:m,1:n) )

  return
end
subroutine r8mat_max_index ( m, n, a, i, j )

!*****************************************************************************80
!
!! R8MAT_MAX_INDEX returns the location of the maximum entry of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the maximum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj

  i = -1
  j = -1

  do jj = 1, n
    do ii = 1, m
      if ( ii == 1 .and. jj == 1 ) then
        i = ii
        j = jj
      else if ( a(i,j) < a(ii,jj) ) then
        i = ii
        j = jj
      end if
    end do
  end do

  return
end
function r8mat_maxcol_minrow ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_MAXCOL_MINROW gets the maximum column minimum row of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    R8MAT_MAXCOL_MINROW = max ( 1 <= I <= N ) ( min ( 1 <= J <= M ) A(I,J) )
!
!    For a given matrix, R8MAT_MAXCOL_MINROW <= R8MAT_MINROW_MAXCOL.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) R8MAT_MAXCOL_MINROW, the maximum column
!    minimum row entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) r8mat_maxcol_minrow
  real    ( kind = 8 ) r8mat_minrow

  r8mat_maxcol_minrow = 0.0D+00

  do i = 1, m

    r8mat_minrow = minval ( a(i,1:n) )

    if ( i == 1 ) then
      r8mat_maxcol_minrow = r8mat_minrow
    else
      r8mat_maxcol_minrow = max ( r8mat_maxcol_minrow, r8mat_minrow )
    end if

  end do

  return
end
function r8mat_maxrow_mincol ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_MAXROW_MINCOL gets the maximum row minimum column of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    R8MAT_MAXROW_MINCOL = max ( 1 <= J <= N ) ( min ( 1 <= I <= M ) A(I,J) )
!
!    For a given matrix, R8MAT_MAXROW_MINCOL <= R8MAT_MINCOL_MAXROW.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) R8MAT_MAXROW_MINCOL, the maximum row 
!    minimum column entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  real    ( kind = 8 ) r8mat_maxrow_mincol
  real    ( kind = 8 ) r8mat_mincol

  r8mat_maxrow_mincol = 0.0D+00

  do j = 1, n

    r8mat_mincol = minval ( a(1:m,j) )

    if ( j == 1 ) then
      r8mat_maxrow_mincol = r8mat_mincol
    else
      r8mat_maxrow_mincol = max ( r8mat_maxrow_mincol, r8mat_mincol )
    end if

  end do

  return
end
function r8mat_min ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_MIN returns the minimum entry of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) R8MAT_MIN, the minimum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) r8mat_min

  r8mat_min = minval ( a(1:m,1:n) )

  return
end
subroutine r8mat_min_index ( m, n, a, i, j )

!*****************************************************************************80
!
!! R8MAT_MIN_INDEX returns the location of the minimum entry of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the minimum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj

  i = -1
  j = -1

  do jj = 1, n
    do ii = 1, m
      if ( ii == 1 .and. jj == 1 ) then
        i = ii
        j = jj
      else if ( a(ii,jj) < a(i,j) ) then
        i = ii
        j = jj
      end if
    end do
  end do

  return
end
function r8mat_mincol_maxrow ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_MINCOL_MAXROW gets the minimum column maximum row of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    R8MAT_MINCOL_MAXROW = min ( 1 <= I <= N ) ( max ( 1 <= J <= M ) A(I,J) )
!
!    For a given matrix, R8MAT_MAXROW_MINCOL <= R8MAT_MINCOL_MAXROW.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) R8MAT_MINCOL_MAXROW, the minimum column
!    maximum row entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) r8mat_mincol_maxrow
  real    ( kind = 8 ) r8mat_maxrow

  r8mat_mincol_maxrow = 0.0D+00

  do i = 1, m

    r8mat_maxrow = maxval ( a(i,1:n) )

    if ( i == 1 ) then
      r8mat_mincol_maxrow = r8mat_maxrow
    else
      r8mat_mincol_maxrow = min ( r8mat_mincol_maxrow, r8mat_maxrow )
    end if

  end do

  return
end
function r8mat_minrow_maxcol ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_MINROW_MAXCOL gets the minimum row maximum column of an M by N R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    R8MAT_MINROW_MAXCOL = min ( 1 <= J <= N ) ( max ( 1 <= I <= M ) A(I,J) )
!
!    For a given matrix, R8MAT_MAXCOL_MINROW <= R8MAT_MINROW_MAXCOL.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) R8MAT_MINROW_MAXCOL, the minimum row 
!    maximum column entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  real    ( kind = 8 ) r8mat_minrow_maxcol
  real    ( kind = 8 ) r8mat_maxcol

  r8mat_minrow_maxcol = 0.0D+00

  do j = 1, n

    r8mat_maxcol = maxval ( a(1:m,j) )

    if ( j == 1 ) then
      r8mat_minrow_maxcol = r8mat_maxcol
    else
      r8mat_minrow_maxcol = min ( r8mat_minrow_maxcol, r8mat_maxcol )
    end if

  end do

  return
end
subroutine r8mat_mm ( n1, n2, n3, a, b, c )

!*****************************************************************************80
!
!! R8MAT_MM multiplies two R8MAT's.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    In FORTRAN90, this operation is more efficiently done by the
!    command:
!
!      C(1:N1,1:N3) = MATMUL ( A(1:N1,1;N2), B(1:N2,1:N3) )
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, the order of the matrices.
!
!    Input, real ( kind = 8 ) A(N1,N2), B(N2,N3), the matrices to multiply.
!
!    Output, real ( kind = 8 ) C(N1,N3), the product matrix C = A * B.
!
  implicit none

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  real    ( kind = 8 ) a(n1,n2)
  real    ( kind = 8 ) b(n2,n3)
  real    ( kind = 8 ) c(n1,n3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  do i = 1, n1
    do j = 1, n3
      c(i,j) = 0.0D+00
      do k = 1, n2
        c(i,j) = c(i,j) + a(i,k) * b(k,j)
      end do
    end do
  end do

  return
end
subroutine r8mat_mtv ( m, n, a, x, y )

!*****************************************************************************80
!
!! R8MAT_MTV multiplies a transposed matrix times a vector
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix.
!
!    Input, real ( kind = 8 ) X(M), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) Y(N), the product A'*X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) x(m)
  real    ( kind = 8 ) y(n)

  y(1:n) = matmul ( transpose ( a(1:m,1:n) ), x(1:m) )

  return
end
subroutine r8mat_mv ( m, n, a, x, y )

!*****************************************************************************80
!
!! R8MAT_MV multiplies a matrix times a vector.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    In FORTRAN90, this operation can be more efficiently carried
!    out by the command
!
!      Y(1:M) = MATMUL ( A(1:M,1:N), X(1:N) )
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) Y(M), the product A*X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) y(m)

  do i = 1, m
    y(i) = 0.0D+00
    do j = 1, n
      y(i) = y(i) + a(i,j) * x(j)
    end do
  end do

  return
end
subroutine r8mat_nint ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NINT rounds the entries of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!
!    Input/output, real ( kind = 8 ) A(M,N), the matrix to be NINT'ed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)

  a(1:m,1:n) = real ( nint ( a(1:m,1:n) ), kind = 8 )

  return
end
function r8mat_norm_eis ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_EIS returns the EISPACK norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    The EISPACK norm is defined as:
!
!      R8MAT_NORM_EIS =
!        sum ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) )
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose EISPACK norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_EIS, the EISPACK norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) r8mat_norm_eis

  r8mat_norm_eis = sum ( abs ( a(1:m,1:n) ) )

  return
end
function r8mat_norm_fro ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    The Frobenius norm is defined as
!
!      R8MAT_NORM_FRO = sqrt (
!        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)**2 )
!
!    The matrix Frobenius norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
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
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose Frobenius 
!    norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_FRO, the Frobenius norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) r8mat_norm_fro

  r8mat_norm_fro = sqrt ( sum ( a(1:m,1:n)**2 ) )

  return
end
function r8mat_norm_l1 ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_L1 returns the matrix L1 norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    The matrix L1 norm is defined as:
!
!      R8MAT_NORM_L1 = max ( 1 <= J <= N )
!        sum ( 1 <= I <= M ) abs ( A(I,J) ).
!
!    The matrix L1 norm is derived from the vector L1 norm, and
!    satisifies:
!
!      r8vec_norm_l1 ( A * x ) <= r8mat_norm_l1 ( A ) * r8vec_norm_l1 ( x ).
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose L1 norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_L1, the L1 norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) j
  real    ( kind = 8 ) r8mat_norm_l1

  r8mat_norm_l1 = 0.0D+00

  do j = 1, n
    r8mat_norm_l1 = max ( r8mat_norm_l1, sum ( abs ( a(1:m,j) ) ) )
  end do

  return
end
function r8mat_norm_l2 ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_L2 returns the matrix L2 norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    The matrix L2 norm is defined as:
!
!      R8MAT_NORM_L2 = sqrt ( max ( 1 <= I <= M ) LAMBDA(I) )
!
!    where LAMBDA contains the eigenvalues of A * A'.
!
!    The matrix L2 norm is derived from the vector L2 norm, and
!    satisifies:
!
!      r8vec_norm_l2 ( A * x ) <= r8mat_norm_l2 ( A ) * r8vec_norm_l2 ( x ).
!
!  Modified:
!
!    19 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_L2, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) b(m,m)
  real    ( kind = 8 ) diag(m)
  real    ( kind = 8 ) r8mat_norm_l2
!
!  Compute B = A * A'.
!
  b(1:m,1:m) = matmul ( a(1:m,1:n), transpose ( a(1:m,1:n) ) )
!
!  Diagonalize B.
!
  call r8mat_symm_jacobi ( m, b )
!
!  Find the maximum eigenvalue, and take its square root.
!
  call r8mat_diag_get_vector ( m, b, diag )

  r8mat_norm_l2 = sqrt ( maxval ( diag(1:m) ) )

  return
end
function r8mat_norm_li ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_NORM_LI returns the matrix L-infinity norm of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    The matrix L-infinity norm is defined as:
!
!      R8MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
!
!    The matrix L-infinity norm is derived from the vector L-infinity norm,
!    and satisifies:
!
!      r8vec_norm_li ( A * x ) <= r8mat_norm_li ( A ) * r8vec_norm_li ( x ).
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose L-infinity
!    norm is desired.
!
!    Output, real ( kind = 8 ) R8MAT_NORM_LI, the L-infinity norm of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) r8mat_norm_li

  r8mat_norm_li = 0.0D+00

  do i = 1, m
    r8mat_norm_li = max ( r8mat_norm_li, sum ( abs ( a(i,1:n) ) ) )
  end do

  return
end
subroutine r8mat_orth_uniform ( n, seed, a )

!*****************************************************************************80
!
!! R8MAT_ORTH_UNIFORM returns a random orthogonal R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    Thanks to Eugene Petrov, B I Stepanov Institute of Physics,
!    National Academy of Sciences of Belarus, for convincingly
!    pointing out the severe deficiencies of an earlier version of
!    this routine.
!
!    Essentially, the computation involves saving the Q factor of the
!    QR factorization of a matrix whose entries are normally distributed.
!    However, it is only necessary to generate this matrix a column at
!    a time, since it can be shown that when it comes time to annihilate
!    the subdiagonal elements of column K, these (transformed) elements of
!    column K are still normally distributed random values.  Hence, there
!    is no need to generate them at the beginning of the process and
!    transform them K-1 times.
!
!    For computational efficiency, the individual Householder transformations
!    could be saved, as recommended in the reference, instead of being
!    accumulated into an explicit matrix format.
!
!  Properties:
!
!    The inverse of A is equal to A'.
!
!    A * A'  = A' * A = I.
!
!    Columns and rows of A have unit Euclidean norm.
!
!    Distinct pairs of columns of A are orthogonal.
!
!    Distinct pairs of rows of A are orthogonal.
!
!    The L2 vector norm of A*x = the L2 vector norm of x for any vector x.
!
!    The L2 matrix norm of A*B = the L2 matrix norm of B for any matrix B.
!
!    The determinant of A is +1 or -1.
!
!    All the eigenvalues of A have modulus 1.
!
!    All singular values of A are 1.
!
!    All entries of A are between -1 and 1.
!
!  Modified:
!
!    04 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pete Stewart,
!    Efficient Generation of Random Orthogonal Matrices With an Application
!    to Condition Estimators,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 3, June 1980, pages 403-409.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) A(N,N), the orthogonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) v(n)
  real    ( kind = 8 ) x(n)
!
!  Start with A = the identity matrix.
!
  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 1.0D+00
      else
        a(i,j) = 0.0D+00
      end if
    end do
  end do
!
!  Now behave as though we were computing the QR factorization of
!  some other random matrix.  Generate the N elements of the first column,
!  compute the Householder matrix H1 that annihilates the subdiagonal elements,
!  and set A := A * H1' = A * H.
!
!  On the second step, generate the lower N-1 elements of the second column,
!  compute the Householder matrix H2 that annihilates them,
!  and set A := A * H2' = A * H2 = H1 * H2.
!
!  On the N-1 step, generate the lower 2 elements of column N-1,
!  compute the Householder matrix HN-1 that annihilates them, and
!  and set A := A * H(N-1)' = A * H(N-1) = H1 * H2 * ... * H(N-1).
!  This is our random orthogonal matrix.
!
  do j = 1, n-1
!
!  Set the vector that represents the J-th column to be annihilated.
!
    x(1:j-1) = 0.0D+00

    do i = j, n
      x(i) = r8_normal_01 ( seed )
    end do
!
!  Compute the vector V that defines a Householder transformation matrix
!  H(V) that annihilates the subdiagonal elements of X.
!
    call r8vec_house_column ( n, x, j, v )
!
!  Postmultiply the matrix A by H'(V) = H(V).
!
    call r8mat_house_axh ( n, a, v, a )

  end do

  return
end
subroutine r8mat_plot ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PLOT "plots" an R8MAT, with an optional title.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_ten = 10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character r8mat_plot_symbol
  character ( len = 70 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do jlo = 1, n, 70
    jhi = min ( jlo + 70-1, n )
    write ( *, '(a)' ) ' '
    write ( *, '(8x,2x,70i1)' ) ( mod ( j, i4_ten ), j = jlo, jhi )
    write ( *, '(a)' ) ' '

    do i = 1, m
      do j = jlo, jhi
        string(j+1-jlo:j+1-jlo) = r8mat_plot_symbol ( a(i,j) )
      end do
      write ( *, '(i8,2x,a)' ) i, string(1:jhi+1-jlo)
    end do
  end do

  return
end
function r8mat_plot_symbol ( r )

!*****************************************************************************80
!
!! R8MAT_PLOT_SYMBOL returns a symbol for a double precision number.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, a value whose symbol is desired.
!
!    Output, character R8MAT_PLOT_SYMBOL, is
!    '-' if R is negative,
!    '0' if R is zero,
!    '+' if R is positive.
!
  implicit none

  character r8mat_plot_symbol
  real    ( kind = 8 ) r

  if ( r < 0.0D+00 ) then
    r8mat_plot_symbol = '-'
  else if ( r == 0.0D+00 ) then
    r8mat_plot_symbol = '0'
  else if ( 0.0D+00 < r ) then
    r8mat_plot_symbol = '+'
  end if

  return
end
subroutine r8mat_poly_char ( n, a, p )

!*****************************************************************************80
!
!! R8MAT_POLY_CHAR computes the characteristic polynomial of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    15 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the N by N matrix.
!
!    Output, real ( kind = 8 ) P(0:N), the coefficients of the characteristic
!    polynomial of A.  P(N) contains the coefficient of X**N
!    (which will be 1), P(I) contains the coefficient of X**I,
!    and P(0) contains the constant term.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real    ( kind = 8 ) p(0:n)
  real    ( kind = 8 ) r8mat_trace
  real    ( kind = 8 ) trace
  real    ( kind = 8 ) work1(n,n)
  real    ( kind = 8 ) work2(n,n)
!
!  Initialize WORK1 to the identity matrix.
!
  call r8mat_identity ( n, work1 )

  p(n) = 1.0D+00

  do order = n-1, 0, -1
!
!  Work2 = A * WORK1.
!
    work2(1:n,1:n) = matmul ( a(1:n,1:n), work1(1:n,1:n) )
!
!  Take the trace.
!
    trace = r8mat_trace ( n, work2 )
!
!  P(ORDER) = -Trace ( WORK2 ) / ( N - ORDER )
!
    p(order) = -trace / real ( n - order, kind = 8 )
!
!  WORK1 := WORK2 + P(ORDER) * Identity.
!
    work1(1:n,1:n) = work2(1:n,1:n)

    do i = 1, n
      work1(i,i) = work1(i,i) + p(order)
    end do

  end do

  return
end
subroutine r8mat_power ( n, a, npow, b )

!*****************************************************************************80
!
!! R8MAT_POWER computes a nonnegative power of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    The algorithm is:
!
!      B = I
!      do NPOW times:
!        B = A * B
!      end
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix to be raised to a power.
!
!    Input, integer ( kind = 4 ) NPOW, the power to which A is to be raised.
!    NPOW must be nonnegative.
!
!    Output, real ( kind = 8 ) B(N,N), the value of A**NPOW.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) ipow
  integer ( kind = 4 ) npow

  if ( npow < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_POWER - Fatal error!'
    write ( *, '(a)' ) '  Input value of NPOW < 0.'
    write ( *, '(a,i8)' ) '  NPOW = ', npow
    stop
  end if

  call r8mat_identity ( n, b )

  do ipow = 1, npow
    b(1:n,1:n) = matmul ( a(1:n,1:n), b(1:n,1:n) )
  end do

  return
end
subroutine r8mat_power_method ( n, a, r, v )

!*****************************************************************************80
!
!! R8MAT_POWER_METHOD applies the power method to an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    If the power method has not converged, then calling the routine
!    again immediately with the output from the previous call will
!    continue the iteration.
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
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix.
!
!    Output, real ( kind = 8 ) R, the estimated eigenvalue.
!
!    Input/output, real ( kind = 8 ) V(N), on input, an estimate
!    for the eigenvector.  On output, an improved estimate for the
!    eigenvector.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) av(n)
  real    ( kind = 8 ) eps
  integer ( kind = 4 ), parameter :: i4_one = 1
  integer ( kind = 4 ) it
  real    ( kind = 8 ), parameter :: it_eps = 0.0001D+00
  integer ( kind = 4 ), parameter :: it_max = 100
  integer ( kind = 4 ), parameter :: it_min = 10
  integer ( kind = 4 ) j
  real    ( kind = 8 ) r
  real    ( kind = 8 ) r2
  real    ( kind = 8 ) r_old
  real    ( kind = 8 ) v(n)

  eps = sqrt ( epsilon ( 1.0D+00 ) )

  r = sqrt ( sum ( v(1:n)**2 ) )

  if ( r == 0.0D+00 ) then
    v(1:n) = 1.0D+00
    r = sqrt ( real ( n, kind = 8 ) )
  end if

  v(1:n) = v(1:n) / r

  do it = 1, it_max

    av(1:n) = matmul ( a(1:n,1:n), v(1:n) )

    r_old = r
    r = sqrt ( sum ( av(1:n)**2 ) )

    if ( it_min < it ) then
      if ( abs ( r - r_old ) <= it_eps * ( 1.0D+00 + abs ( r ) ) ) then
        exit
      end if
    end if

    v(1:n) = av(1:n)

    if ( r /= 0.0D+00 ) then
      v(1:n) = v(1:n) / r
    end if
!
!  Perturb V a bit, to avoid cases where the initial guess is exactly
!  the eigenvector of a smaller eigenvalue.
!
    if ( it < it_max / 2 ) then
      j = 1 + mod ( it - i4_one, n )
      v(j) = v(j) + eps * ( 1.0D+00 + abs ( v(j) ) )
      r2 = sqrt ( sum ( v(1:n)**2 ) )
      v(1:n) = v(1:n) / r2
    end if

  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    26 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine r8mat_print2 ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_PRINT2 prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    27 August 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix to be printed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) amax
  real    ( kind = 8 ) amin
  integer ( kind = 4 ) i
  character ( len = 10 ) iform
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  logical integ
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) lmax
  integer ( kind = 4 ) npline
  real    ( kind = 8 ) r8_log_10
!
!  Check if all entries are integral.
!
  integ = .true.

  do i = 1, m
    do j = 1, n

      if ( integ ) then
        if ( a(i,j) /= real ( int ( a(i,j) ), kind = 8 ) ) then
          integ = .false.
        end if
      end if

    end do
  end do
!
!  Find the maximum and minimum entries.
!
  amax = maxval ( a(1:m,1:n) )
  amin = minval ( a(1:m,1:n) )
!
!  Use the information about the maximum size of an entry to
!  compute an intelligent format for use with integer entries.
!
!  Later, we might also do this for real matrices.
!
  lmax = int ( r8_log_10 ( amax ) )

  if ( integ ) then
    npline = 79 / ( lmax + 3 )
    write ( iform, '(''('',i2,''I'',i2,'')'')' ) npline, lmax+3
  else
    npline = 5
    iform = ' '
  end if
!
!  Print a scalar quantity.
!
  if ( m == 1 .and. n == 1 ) then

    if ( integ ) then
      write ( *, iform ) int ( a(1,1) )
    else
      write ( *, '(2x,g14.6)' ) a(1,1)
    end if
!
!  Column vector of length M,
!
  else if ( n == 1 ) then

    do ilo = 1, m, npline

      ihi = min ( ilo+npline-1, m )

      if ( integ ) then
        write ( *, iform ) ( int ( a(i,1) ), i = ilo, ihi )
      else
        write ( *, '(2x,5g14.6)' ) a(ilo:ihi,1)
      end if

    end do
!
!  Row vector of length N,
!
  else if ( m == 1 ) then

    do jlo = 1, n, npline

      jhi = min ( jlo+npline-1, n )

      if ( integ ) then
        write ( *, iform ) int ( a(1,jlo:jhi) )
      else
        write ( *, '(2x,5g14.6)' ) a(1,jlo:jhi)
      end if

    end do
!
!  M by N Array
!
  else

    do jlo = 1, n, npline

      jhi = min ( jlo+npline-1, n )

      if ( npline < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8,a,i8)' ) 'Matrix columns ', jlo, ' to ', jhi
        write ( *, '(a)' ) ' '
      end if

      do i = 1, m

        if ( integ ) then
          write ( *, iform ) int ( a(i,jlo:jhi) )
        else
          write ( *, '(2x,5g14.6)' ) a(i,jlo:jhi)
        end if

      end do
    end do

  end if

  return
end
subroutine r8mat_solve ( n, rhs_num, a, info )

!*****************************************************************************80
!
!! R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) RHS_NUM, the number of right hand sides.  
!    RHS_NUM must be at least 0.
!
!    Input/output, real ( kind = 8 ) A(N,N+RHS_NUM), contains in rows and 
!    columns 1 to N the coefficient matrix, and in columns N+1 through
!    N+RHS_NUM, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) rhs_num

  real    ( kind = 8 ) a(n,n+rhs_num)
  real    ( kind = 8 ) apivot
  real    ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot
  integer ( kind = 4 ) j

  info = 0

  do j = 1, n
!
!  Choose a pivot row IPIVOT.
!
    ipivot = j
    apivot = a(j,j)

    do i = j + 1, n
      if ( abs ( apivot ) < abs ( a(i,j) ) ) then
        apivot = a(i,j)
        ipivot = i
      end if
    end do

    if ( apivot == 0.0D+00 ) then
      info = j
      return
    end if
!
!  Interchange.
!
    call r8vec_swap ( n+rhs_num, a(ipivot,1:n+rhs_num), a(j,1:n+rhs_num) )
!
!  A(J,J) becomes 1.
!
    a(j,j) = 1.0D+00
    a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
!
!  A(I,J) becomes 0.
!
    do i = 1, n

      if ( i /= j ) then

        factor = a(i,j)
        a(i,j) = 0.0D+00
        a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)

      end if

    end do

  end do

  return
end
subroutine r8mat_solve_2d ( a, b, det, x )

!*****************************************************************************80
!
!! R8MAT_SOLVE_2D solves a 2 by 2 linear system using Cramer's rule.
!
!  Discussion:
!
!    If the determinant DET is returned as zero, then the matrix A is 
!    singular, and does not have an inverse.  In that case, X is 
!    returned as the zero vector.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Modified:
!
!    16 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2,2), the matrix.
!
!    Input, real ( kind = 8 ) B(2), the right hand side.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
!
!    Output, real ( kind = 8 ) X(2), the solution of the system,
!    if DET is nonzero.
!
  implicit none

  real    ( kind = 8 ) a(2,2)
  real    ( kind = 8 ) b(2)
  real    ( kind = 8 ) det
  real    ( kind = 8 ) x(2)
!
!  Compute the determinant.
!
  det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then
    x(1:2) = 0.0D+00
    return
  end if
!
!  Compute the solution.
!
  x(1) = (  a(2,2) * b(1) - a(1,2) * b(2) ) / det
  x(2) = ( -a(2,1) * b(1) + a(1,1) * b(2) ) / det

  return
end
subroutine r8mat_solve_3d ( a, b, det, x )

!*****************************************************************************80
!
!! R8MAT_SOLVE_3D solves a 3 by 3 linear system using Cramer's rule.
!
!  Discussion:
!
!    If the determinant DET is returned as zero, then the matrix A is 
!    singular, and does not have an inverse.  In that case, X is 
!    returned as the zero vector.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(3,3), the matrix.
!
!    Input, real ( kind = 8 ) B(3), the right hand side.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
!
!    Output, real ( kind = 8 ) X(3), the solution of the system,
!    if DET is nonzero.
!
  implicit none

  real    ( kind = 8 ) a(3,3)
  real    ( kind = 8 ) b(3)
  real    ( kind = 8 ) det
  real    ( kind = 8 ) x(3)
!
!  Compute the determinant.
!
  det =  a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
       + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
       + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then
    x(1:3) = 0.0D+00
    return
  end if
!
!  Compute the solution.
!
  x(1) = (   ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) * b(1) &
           - ( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) * b(2) &
           + ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) * b(3) ) / det

  x(2) = ( - ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) * b(1) &
           + ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) * b(2) &
           - ( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) * b(3) ) / det

  x(3) = (   ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) * b(1) &
           - ( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) * b(2) &
           + ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) * b(3) ) / det

  return
end
subroutine r8mat_solve2 ( n, a, b, x, ierror )

!*****************************************************************************80
!
!! R8MAT_SOLVE2 computes the solution of an N by N linear system.
!
!  Discussion:
!
!    The linear system may be represented as
!
!      A*X = B
!
!    If the linear system is singular, but consistent, then the routine will
!    still produce a solution.
!
!  Modified:
!
!    29 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of equations.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, A is the coefficient matrix to be inverted.
!    On output, A has been overwritten.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, B is the right hand side of the system.
!    On output, B has been overwritten.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
!    Output, integer ( kind = 4 ) IERROR.
!    0, no error detected.
!    1, consistent singularity.
!    2, inconsistent singularity.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) amax
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) ipiv(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) x(n)

  ierror = 0

  ipiv(1:n) = 0
  x(1:n) = 0.0D+00
!
!  Process the matrix.
!
  do k = 1, n
!
!  In column K:
!    Seek the row IMAX with the properties that:
!      IMAX has not already been used as a pivot;
!      A(IMAX,K) is larger in magnitude than any other candidate.
!
    amax = 0.0D+00
    imax = 0
    do i = 1, n
      if ( ipiv(i) == 0 ) then
        if ( amax < abs ( a(i,k) ) ) then
          imax = i
          amax = abs ( a(i,k) )
        end if
      end if
    end do
!
!  If you found a pivot row IMAX, then,
!    eliminate the K-th entry in all rows that have not been used for pivoting.
!
    if ( imax /= 0 ) then

      ipiv(imax) = k
      a(imax,k+1:n) = a(imax,k+1:n) / a(imax,k)
      b(imax) = b(imax) / a(imax,k)
      a(imax,k) = 1.0D+00

      do i = 1, n

        if ( ipiv(i) == 0 ) then
          a(i,k+1:n) = a(i,k+1:n) - a(i,k) * a(imax,k+1:n)
          b(i) = b(i) - a(i,k) * b(imax)
          a(i,k) = 0.0D+00
        end if

      end do

    end if

  end do
!
!  Now, every row with nonzero IPIV begins with a 1, and
!  all other rows are all zero.  Begin solution.
!
  do j = n, 1, -1

    imax = 0
    do k = 1, n
      if ( ipiv(k) == j ) then
        imax = k
      end if
    end do

    if ( imax == 0 ) then

      x(j) = 0.0D+00

      if ( b(j) == 0.0D+00 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_SOLVE2 - Warning:'
        write ( *, '(a,i8)' ) '  Consistent singularity, equation = ', j
      else
        ierror = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_SOLVE2 - Error:'
        write ( *, '(a,i8)' ) '  Inconsistent singularity, equation = ', j
      end if

    else

      x(j) = b(imax)

      do i = 1, n
        if ( i /= imax ) then
          b(i) = b(i) - a(i,j) * x(j)
        end if
      end do

    end if

  end do

  return
end
subroutine r8mat_symm_eigen ( n, x, q, a )

!*****************************************************************************80
!
!! R8MAT_SYMM_EIGEN returns a symmetric matrix with given eigensystem.
!
!  Discussion:
!
!    The user must supply the desired eigenvalue vector, and the desired
!    eigenvector matrix.  The eigenvector matrix must be orthogonal.  A
!    suitable random orthogonal matrix can be generated by R8MAT_ORTH_UNIFORM. 
!
!  Modified:
!
!    18 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input, real ( kind = 8 ) X(N), the desired eigenvalues for the matrix.
!
!    Input, real ( kind = 8 ) Q(N,N), the eigenvector matrix of A.
!
!    Output, real ( kind = 8 ) A(N,N), a symmetric matrix with
!    eigenvalues X and eigenvectors the columns of Q.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) q(n,n)
  real    ( kind = 8 ) x(n)
!
!  Set A = Q * Lambda * Q'.
!
  a(1:n,1:n) = 0.0D+00

  do i = 1, n
    do j = 1, n
      do k = 1, n
        a(i,j) = a(i,j) + q(i,k) * x(k) * q(j,k)
      end do
    end do
  end do

  return
end
subroutine r8mat_symm_jacobi ( n, a )

!*****************************************************************************80
!
!! R8MAT_SYMM_JACOBI applies Jacobi eigenvalue iteration to a symmetric matrix.
!
!  Discussion:
!
!    This code was modified so that it treats as zero the off-diagonal
!    elements that are sufficiently close to, but not exactly, zero.
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of A.
!
!    Input/output, real ( kind = 8 ) A(N,N), a symmetric N by N matrix.
!    On output, the matrix has been overwritten by an approximately
!    diagonal matrix, with the eigenvalues on the diagonal.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) c
  real    ( kind = 8 ) r8mat_norm_fro
  real    ( kind = 8 ), parameter :: eps = 0.00001D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 100
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) norm_fro
  real    ( kind = 8 ) s
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) t
  real    ( kind = 8 ) t1
  real    ( kind = 8 ) t2
  real    ( kind = 8 ) u

  norm_fro = r8mat_norm_fro ( n, n, a )

  it = 0

  do

    it = it + 1

    do i = 1, n
      do j = 1, i - 1

        if ( eps * norm_fro < abs ( a(i,j) ) + abs ( a(j,i) ) ) then

          u = ( a(j,j) - a(i,i) ) / ( a(i,j) + a(j,i) )

          t = sign ( 1.0D+00, u ) / ( abs ( u ) + sqrt ( u * u + 1.0D+00 ) )
          c = 1.0D+00 / sqrt ( t * t + 1.0D+00 )
          s = t * c
!
!  A -> A * Q.
!
          do k = 1, n
            t1 = a(i,k)
            t2 = a(j,k)
            a(i,k) = t1 * c - t2 * s
            a(j,k) = t1 * s + t2 * c
          end do
!
!  A -> QT * A
!
          do k = 1, n
            t1 = a(k,i)
            t2 = a(k,j)
            a(k,i) = c * t1 - s * t2
            a(k,j) = s * t1 + c * t2
          end do

        end if
      end do
    end do
!
!  Test the size of the off-diagonal elements.
!
    sum2 = 0.0D+00
    do i = 1, n
      do j = 1, i - 1
        sum2 = sum2 + abs ( a(i,j) )
      end do
    end do

    if ( sum2 <= eps * ( norm_fro + 1.0D+00 ) ) then
      exit
    end if

    if ( it_max <= it ) then
      exit
    end if

  end do

  return
end
subroutine r8mat_to_r8plu ( n, a, pivot, lu, info )

!*****************************************************************************80
!
!! R8MAT_TO_R8PLU factors a general R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    This routine is a simplified version of the LINPACK routine DGEFA.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix to be factored.
!
!    Output, integer ( kind = 4 ) PIVOT(N), a vector of pivot indices.
!
!    Output, real ( kind = 8 ) LU(N,N), an upper triangular matrix U and
!    the multipliers L which were used to obtain it.  The factorization
!    can be written A = L * U, where L is a product of permutation and
!    unit lower triangular matrices and U is upper triangular.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real    ( kind = 8 ) lu(n,n)
  real    ( kind = 8 ) temp

  lu(1:n,1:n) = a(1:n,1:n)

  info = 0

  do k = 1, n - 1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k + 1, n
      if ( abs ( lu(l,k) ) < abs ( lu(i,k) ) ) then
        l = i
      end if
    end do

    pivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( lu(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_TO_R8PLU - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      temp    = lu(l,k)
      lu(l,k) = lu(k,k)
      lu(k,k) = temp
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    lu(k+1:n,k) = -lu(k+1:n,k) / lu(k,k)
!
!  Row elimination with column indexing.
!
    do j = k + 1, n

      if ( l /= k ) then
        temp    = lu(l,j)
        lu(l,j) = lu(k,j)
        lu(k,j) = temp
      end if

      lu(k+1:n,j) = lu(k+1:n,j) + lu(k+1:n,k) * lu(k,j)

    end do

  end do

  pivot(n) = n

  if ( lu(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_TO_R8PLU - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
  end if

  return
end
function r8mat_trace ( n, a )

!*****************************************************************************80
!
!! R8MAT_TRACE computes the trace of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    The trace of a square matrix is the sum of the diagonal elements.
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
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix whose trace is desired.
!
!    Output, real ( kind = 8 ) R8MAT_TRACE, the trace of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) r8mat_trace

  r8mat_trace = 0.0D+00
  do i = 1, n
    r8mat_trace = r8mat_trace + a(i,i)
  end do

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine r8mat_u_inverse ( n, a, b )

!*****************************************************************************80
!
!! R8MAT_U_INVERSE inverts an upper triangular R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    An upper triangular matrix is a matrix whose only nonzero entries
!    occur on or above the diagonal.
!
!    The inverse of an upper triangular matrix is an upper triangular matrix.
!
!  Modified:
!
!    11 December 2004
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the upper triangular matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = n, 1, -1

    do i = n, 1, -1

      if ( j < i ) then
        b(i,j) = 0.0D+00
      else if ( i == j ) then
        b(i,j) = 1.0D+00 / a(i,j)
      else
        b(i,j) = -dot_product ( a(i,i+1:j), b(i+1:j,j) ) / a(i,i)
      end if

    end do
  end do

  return
end
subroutine r8mat_u1_inverse ( n, a, b )

!*****************************************************************************80
!
!! R8MAT_U1_INVERSE inverts a unit upper triangular R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!    A unit upper triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's below the main diagonal.
!
!    The inverse of a unit upper triangular matrix is also
!    a unit upper triangular matrix.
!
!    This routine can invert a matrix in place, that is, with no extra
!    storage.  If the matrix is stored in A, then the call
!
!      call r8mat_u1_inverse ( n, a, a )
!
!    will result in A being overwritten by its inverse.
!
!  Modified:
!
!    11 December 2004
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of rows and columns in the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the unit upper triangular matrix.
!
!    Output, real ( kind = 8 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  real    ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do j = n, 1, -1

    do i = n, 1, -1

      if ( j < i ) then
        b(i,j) = 0.0D+00
      else if ( i == j ) then
        b(i,j) = 1.0D+00
      else
        b(i,j) = -dot_product ( a(i,i+1:j), b(i+1:j,j) )
      end if

    end do
  end do

  return
end
subroutine r8mat_uniform ( m, n, a, b, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM fills an R8MAT with scaled pseudorandom numbers.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    29 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer  ( kind = 4 )M, N, the number of rows and columns in 
!    the array.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge ( )
      end if

      r(i,j) = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in 
!    the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge ( )
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8mat_vand2 ( n, x, a )

!*****************************************************************************80
!
!! R8MAT_VAND2 returns the N by N row Vandermonde matrix A.
!
!  Discussion:
!
!    The row Vandermonde matrix returned by this routine reads "across" 
!    rather than down.  In particular, each row begins with a 1, followed by
!    some value X, followed by successive powers of X.
!
!  Formula:
!
!    A(I,J) = X(I)**(J-1)
!
!  Properties:
!
!    A is nonsingular if, and only if, the X values are distinct.
!
!    The determinant of A is
!
!      det(A) = product ( 2 <= I <= N ) (
!        product ( 1 <= J <= I-1 ) ( ( X(I) - X(J) ) ) ).
!
!    The matrix A is generally ill-conditioned.
!
!  Example:
!
!    N = 5, X = (2, 3, 4, 5, 6)
!
!    1 2  4   8   16
!    1 3  9  27   81
!    1 4 16  64  256
!    1 5 25 125  625
!    1 6 36 216 1296
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix desired.
!
!    Input, real ( kind = 8 ) X(N), the values that define A.
!
!    Output, real ( kind = 8 ) A(N,N), the N by N row Vandermonde matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x(n)

  do i = 1, n
    do j = 1, n

      if ( j == 1 .and. x(i) == 0.0D+00 ) then
        a(i,j) = 1.0D+00
      else
        a(i,j) = x(i)**(j-1)
      end if

    end do
  end do

  return
end
subroutine r8plu_det ( n, pivot, lu, det )

!*****************************************************************************80
!
!! R8PLU_DET computes the determinant of an R8PLU matrix.
!
!  Discussion:
!
!    The matrix should have been factored by R8MAT_TO_R8PLU.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector computed 
!    by R8MAT_TO_R8PLU.
!
!    Input, real ( kind = 8 ) LU(N,N), the LU factors computed 
!    by R8MAT_TO_R8PLU.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) det
  integer ( kind = 4 ) i
  real    ( kind = 8 ) lu(n,n)
  integer ( kind = 4 ) pivot(n)

  det = 1.0D+00

  do i = 1, n
    det = det * lu(i,i)
    if ( pivot(i) /= i ) then
      det = -det
    end if
  end do

  return
end
subroutine r8plu_inverse ( n, pivot, lu, a_inverse )

!*****************************************************************************80
!
!! R8PLU_INVERSE computes the inverse of an R8PLU matrix.
!
!  Discussion:
!
!    The matrix should have been factored by R8MAT_TO_R8PLU.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector from R8MAT_TO_R8PLu.
!
!    Input, real ( kind = 8 ) LU(N,N), the LU factors computed by R8MAT_TO_R8PLU.
!
!    Output, real ( kind = 8 ) A_INVERSE(N,N), the inverse of the original
!    matrix A that was factored by R8MAT_TO_R8PLU.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a_inverse(n,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) lu(n,n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) work(n)

  a_inverse(1:n,1:n) = lu(1:n,1:n)
!
!  Compute Inverse(U).
!
  do k = 1, n

    a_inverse(k,k)     = 1.0D+00 / a_inverse(k,k)
    a_inverse(1:k-1,k) = -a_inverse(1:k-1,k) * a_inverse(k,k)

    do j = k + 1, n

      temp             = a_inverse(k,j)
      a_inverse(k,j)   = 0.0D+00
      a_inverse(1:k,j) = a_inverse(1:k,j) + temp * a_inverse(1:k,k)

    end do

  end do
!
!  Form Inverse(U) * Inverse(L).
!
  do k = n - 1, 1, -1

    work(k+1:n) = a_inverse(k+1:n,k)
    a_inverse(k+1:n,k) = 0.0D+00

    do j = k + 1, n
      a_inverse(1:n,k) = a_inverse(1:n,k) + a_inverse(1:n,j) * work(j)
    end do

    if ( pivot(k) /= k ) then

      do i = 1, n
        temp                  = a_inverse(i,k)
        a_inverse(i,k)        = a_inverse(i,pivot(k))
        a_inverse(i,pivot(k)) = temp
      end do

    end if

  end do

  return
end
subroutine r8plu_mul ( n, pivot, lu, x, b )

!*****************************************************************************80
!
!! R8PLU_MUL computes A * x using the PLU factors of A.
!
!  Discussion:
!
!    It is assumed that R8MAT_TO_R8PLU has computed the PLU factors of
!    the matrix A.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector computed 
!    by R8MAT_TO_R8PLU.
!
!    Input, real ( kind = 8 ) LU(N,N), the matrix factors computed by
!    R8MAT_TO_R8PLU.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) B(N), the result of the multiplication.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) lu(n,n)
  integer ( kind = 4 ) pivot(n)
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) x(n)

  b(1:n) = x(1:n)
!
!  Y = U * X.
!
  do j = 1, n
    b(1:j-1) = b(1:j-1) + lu(1:j-1,j) * b(j)
    b(j) = lu(j,j) * b(j)
  end do
!
!  B = PL * Y = PL * U * X = A * x.
!
  do j = n - 1, 1, -1

    b(j+1:n) = b(j+1:n) - lu(j+1:n,j) * b(j)

    k = pivot(j)

    if ( k /= j ) then
      temp = b(k)
      b(k) = b(j)
      b(j) = temp
    end if

  end do

  return
end
subroutine r8plu_sol ( n, pivot, lu, b, x )

!*****************************************************************************80
!
!! R8PLU_SOL solves a linear system A*x=b from the PLU factors.
!
!  Discussion:
!
!    The PLU factors should have been computed by R8MAT_TO_R8PLU.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector from R8MAT_TO_R8PLU.
!
!    Input, real ( kind = 8 ) LU(N,N), the LU factors from R8MAT_TO_R8PLU.
!
!    Input, real ( kind = 8 ) B(N), the right hand side vector.
!
!    Output, real ( kind = 8 ) X(N), the solution vector.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) lu(n,n)
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) x(n)
!
!  Solve PL * Y = B.
!
  x(1:n) = b(1:n)

  do k = 1, n - 1

    j = pivot(k)

    if ( j /= k ) then
      temp = x(j)
      x(j) = x(k)
      x(k) = temp
    end if

    x(k+1:n) = x(k+1:n) + lu(k+1:n,k) * x(k)

  end do
!
!  Solve U * X = Y.
!
  do k = n, 1, -1
    x(k) = x(k) / lu(k,k)
    x(1:k-1) = x(1:k-1) - lu(1:k-1,k) * x(k)
  end do

  return
end
subroutine r8plu_to_r8mat ( n, pivot, lu, a )

!*****************************************************************************80
!
!! R8PLU_TO_R8MAT recovers the matrix A that was factored by R8MAT_TO_R8PLU.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector computed 
!    by R8MAT_TO_R8PLU.
!
!    Input, real ( kind = 8 ) LU(N,N), the matrix factors computed by
!    R8MAT_TO_R8PLU.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix whose factors are 
!    represented by LU and PIVOT.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) lu(n,n)
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) temp

  a(1:n,1:n) = 0.0D+00
  do i = 1, n
    a(i,i) = 1.0D+00
  end do

  do j = 1, n

    do i = 1, n
      a(1:i-1,j) = a(1:i-1,j) + lu(1:i-1,i) * a(i,j)
      a(i,j) = lu(i,i) * a(i,j)
    end do
!
!  B = PL * Y = PL * U * X = A * x.
!
    do i = n - 1, 1, -1

      a(i+1:n,j) = a(i+1:n,j) - lu(i+1:n,i) * a(i,j)

      k = pivot(i)

      if ( k /= i ) then
        temp   = a(k,j)
        a(k,j) = a(i,j)
        a(i,j) = temp
      end if

    end do

  end do

  return
end
subroutine r8poly_degree ( na, a, degree )

!*****************************************************************************80
!
!! R8POLY_DEGREE returns the degree of a polynomial.
!
!  Discussion:
!
!    The degree of a polynomial is the index of the highest power
!    of X with a nonzero coefficient.
!
!    The degree of a constant polynomial is 0.  The degree of the
!    zero polynomial is debatable, but this routine returns the
!    degree as 0.
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
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:NA), the coefficients of the polynomials.
!
!    Output, integer ( kind = 4 ) DEGREE, the degree of A.
!
  implicit none

  integer ( kind = 4 ) na

  real    ( kind = 8 ) a(0:na)
  integer ( kind = 4 ) degree

  degree = na

  do while ( 0 < degree )

    if ( a(degree) /= 0.0D+00 ) then
      return
    end if

    degree = degree - 1

  end do

  return
end
subroutine r8poly_deriv ( n, c, p, cp )

!*****************************************************************************80
!
!! R8POLY_DERIV returns the derivative of a polynomial.
!
!  Modified:
!
!    21 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Input, real ( kind = 8 ) C(0:N), the polynomial coefficients.
!    C(I) is the coefficient of X**I.
!
!    Input, integer ( kind = 4 ) P, the order of the derivative.
!    0 means no derivative is taken.
!    1 means first derivative,
!    2 means second derivative and so on.
!    Values of P less than 0 are meaningless.  Values of P greater
!    than N are meaningful, but the code will behave as though the
!    value of P was N+1.
!
!    Output, real ( kind = 8 ) CP(0:N-P), the polynomial coefficients of
!    the derivative.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) c(0:n)
  real    ( kind = 8 ) cp(0:*)
  real    ( kind = 8 ) cp_temp(0:n)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p

  if ( n < p ) then
    return
  end if

  cp_temp(0:n) = c(0:n)

  do d = 1, p
    do i = 0, n - d
      cp_temp(i) = real ( i + 1, kind = 8 ) * cp_temp(i+1)
    end do
    cp_temp(n-d+1) = 0.0D+00
  end do

  cp(0:n-p) = cp_temp(0:n-p)

  return
end
subroutine r8poly_lagrange_0 ( npol, xpol, xval, wval )

!*****************************************************************************80
!
!! R8POLY_LAGRANGE_0 evaluates the Lagrange factor at a point.
!
!  Formula:
!
!    W(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
!
!  Discussion:
!
!    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
!    polynomial L(IPOL)(X), has the property:
!
!      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
!
!    and may be expressed as:
!
!      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
!
!  Modified:
!
!    21 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas, which 
!    should be distinct.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the Lagrange 
!    factor is to be evaluated.
!
!    Output, real ( kind = 8 ) WVAL, the value of the Lagrange factor at XVAL.
!
  implicit none

  integer ( kind = 4 ) npol

  real    ( kind = 8 ) wval
  real    ( kind = 8 ) xpol(npol)
  real    ( kind = 8 ) xval

  wval = product ( xval - xpol(1:npol) )

  return
end
subroutine r8poly_lagrange_1 ( npol, xpol, xval, dwdx )

!*****************************************************************************80
!
!! R8POLY_LAGRANGE_1 evaluates the first derivative of the Lagrange factor.
!
!  Formula:
!
!    W(XPOL(1:NPOL))(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
!
!    W'(XPOL(1:NPOL))(X)
!      = Sum ( 1 <= J <= NPOL ) Product ( I /= J ) ( X - XPOL(I) )
!
!    We also have the recursion:
!
!      W'(XPOL(1:NPOL))(X) = d/dX ( ( X - XPOL(NPOL) ) * W(XPOL(1:NPOL-1))(X) )
!                    = W(XPOL(1:NPOL-1))(X)
!                    + ( X - XPOL(NPOL) ) * W'(XPOL(1:NPOL-1))(X)
!
!  Modified:
!
!    29 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOL, the number of abscissas.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas, which should 
!    be distinct.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the Lagrange 
!    factor is to be evaluated.
!
!    Output, real ( kind = 8 ) DWDX, the derivative of W with respect to X.
!
  implicit none

  integer ( kind = 4 ) npol

  real    ( kind = 8 ) dwdx
  integer ( kind = 4 ) i
  real    ( kind = 8 ) w
  real    ( kind = 8 ) xpol(npol)
  real    ( kind = 8 ) xval

  dwdx = 0.0D+00
  w = 1.0D+00

  do i = 1, npol

    dwdx = w + ( xval - xpol(i) ) * dwdx
    w = w * ( xval - xpol(i) )

  end do

  return
end
subroutine r8poly_lagrange_2 ( npol, xpol, xval, dw2dx2 )

!*****************************************************************************80
!
!! R8POLY_LAGRANGE_2 evaluates the second derivative of the Lagrange factor.
!
!  Formula:
!
!    W(X)  = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
!
!    W'(X) = Sum ( 1 <= J <= NPOL )
!            Product ( I /= J ) ( X - XPOL(I) )
!
!    W"(X) = Sum ( 1 <= K <= NPOL )
!            Sum ( J =/ K )
!            Product ( I /= K, J ) ( X - XPOL(I) )
!
!    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
!    polynomial L(IPOL)(X), has the property:
!
!      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
!
!    and may be expressed as:
!
!      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
!
!  Modified:
!
!    21 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas, which should
!    be distinct.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the Lagrange 
!    factor is to be evaluated.
!
!    Output, real ( kind = 8 ) DW2DX2, the second derivative of W 
!    with respect to XVAL.
!
  implicit none

  integer ( kind = 4 ) npol

  real    ( kind = 8 ) dw2dx2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) term
  real    ( kind = 8 ) xpol(npol)
  real    ( kind = 8 ) xval

  dw2dx2 = 0.0D+00

  do k = 1, npol

    do j = 1, npol

      if ( j /= k ) then
        term = 1.0D+00

        do i = 1, npol
          if ( i /= j .and. i /= k ) then
            term = term * ( xval - xpol(i) )
          end if
        end do

        dw2dx2 = dw2dx2 + term

      end if

    end do

  end do

  return
end
subroutine r8poly_lagrange_coef ( npol, ipol, xpol, pcof )

!*****************************************************************************80
!
!! R8POLY_LAGRANGE_COEF returns the coefficients of a Lagrange polynomial.
!
!  Discussion:
!
!    Given distinct abscissas XPOL(1:NPOL), the IPOL-th Lagrange
!    polynomial L(IPOL)(X) is defined as the polynomial of degree
!    NPOL - 1 which is 1 at XPOL(IPOL) and 0 at the NPOL - 1 other
!    abscissas.
!
!    A formal representation is:
!
!      L(IPOL)(X) = Product ( 1 <= I <= NPOL, I /= IPOL )
!       ( X - X(I) ) / ( X(IPOL) - X(I) )
!
!    However sometimes it is desirable to be able to write down
!    the standard polynomial coefficients of L(IPOL)(X).
!
!  Modified:
!
!    18 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, integer ( kind = 4 ) IPOL, the index of the polynomial to evaluate.
!    IPOL must be between 1 and NPOL.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas of the
!    Lagrange polynomials.  The entries in XPOL must be distinct.
!
!    Output, real ( kind = 8 ) PCOF(0:NPOL-1), the standard polynomial
!    coefficients of the IPOL-th Lagrange polynomial:
!      L(IPOL)(X) = SUM ( 0 <= I <= NPOL-1 ) PCOF(I) * X**I
!
  implicit none

  integer ( kind = 4 ) npol

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) ipol
  integer ( kind = 4 ) j
  real    ( kind = 8 ) pcof(0:npol-1)
  logical r8vec_distinct
  real    ( kind = 8 ) xpol(npol)
!
!  Make sure IPOL is legal.
!
  if ( ipol < 1 .or. npol < ipol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY_LAGRANGE_COEF - Fatal error!'
    write ( *, '(a)' ) '  1 <= IPOL <= NPOL is required.'
    write ( *, '(a,i8)' ) '  IPOL = ', ipol
    write ( *, '(a,i8)' ) '  NPOL = ', npol
    stop
  end if
!
!  Check that the abscissas are distinct.
!
  if ( .not. r8vec_distinct ( npol, xpol ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY_LAGRANGE_COEF - Fatal error!'
    write ( *, '(a)' ) '  Two or more entries of XPOL are equal:'
    stop
  end if

  pcof(0) = 1.0D+00
  pcof(1:npol-1) = 0.0D+00

  indx = 0

  do i = 1, npol

    if ( i /= ipol ) then

      indx = indx + 1

      do j = indx, 0, -1

        pcof(j) = -xpol(i) * pcof(j) / ( xpol(ipol) - xpol(i) )

        if ( 0 < j ) then
          pcof(j) = pcof(j) + pcof(j-1) / ( xpol(ipol) - xpol(i) )
        end if

      end do

    end if

  end do

  return
end
subroutine r8poly_lagrange_factor ( npol, xpol, xval, wval, dwdx )

!*****************************************************************************80
!
!! R8POLY_LAGRANGE_FACTOR evaluates the polynomial Lagrange factor at a point.
!
!  Formula:
!
!    W(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )
!
!  Discussion:
!
!    Suppose F(X) is at least N times continuously differentiable in the
!    interval [A,B].  Pick NPOL distinct points XPOL(I) in [A,B] and compute
!    the interpolating polynomial P(X) of order NPOL ( and degree NPOL-1)
!    which passes through all the points ( XPOL(I), F(XPOL(I)) ).
!    Then in the interval [A,B], the maximum error
!
!      abs ( F(X) - P(X) )
!
!    is bounded by:
!
!      C * FNMAX * W(X)
!
!    where
!
!      C is a constant,
!      FNMAX is the maximum value of the NPOL-th derivative of F in [A,B],
!      W(X) is the Lagrange factor.
!
!    Thus, the value of W(X) is useful as part of an estimated bound
!    for the interpolation error.
!
!    Note that the Chebyshev abscissas have the property that they minimize
!    the value of W(X) over the interval [A,B].  Hence, if the abscissas may
!    be chosen arbitrarily, the Chebyshev abscissas have this advantage over
!    other choices.
!
!    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
!    polynomial L(IPOL)(X), has the property:
!
!      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )
!
!    and may be expressed as:
!
!      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )
!
!  Modified:
!
!    18 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas, which should 
!    be distinct.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the Lagrange 
!    factor is to be evaluated.
!
!    Output, real ( kind = 8 ) WVAL, the value of the Lagrange factor at XVAL.
!
!    Output, real ( kind = 8 ) DWDX, the derivative of W with respect to XVAL.
!
  implicit none

  integer ( kind = 4 ) npol

  real    ( kind = 8 ) dwdx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) term
  real    ( kind = 8 ) wval
  real    ( kind = 8 ) xpol(npol)
  real    ( kind = 8 ) xval

  wval = product ( xval - xpol(1:npol) )

  dwdx = 0.0D+00

  do i = 1, npol

    term = 1.0D+00

    do j = 1, npol
      if ( i /= j ) then
        term = term * ( xval - xpol(j) )
      end if
    end do

    dwdx = dwdx + term

  end do

  return
end
subroutine r8poly_lagrange_val ( npol, ipol, xpol, xval, pval, dpdx )

!*****************************************************************************80
!
!! R8POLY_LAGRANGE_VAL evaluates the IPOL-th Lagrange polynomial.
!
!  Discussion:
!
!    Given NPOL distinct abscissas, XPOL(1:NPOL), the IPOL-th Lagrange
!    polynomial L(IPOL)(X) is defined as the polynomial of degree
!    NPOL - 1 which is 1 at XPOL(IPOL) and 0 at the NPOL - 1 other
!    abscissas.
!
!    A formal representation is:
!
!      L(IPOL)(X) = Product ( 1 <= I <= NPOL, I /= IPOL )
!       ( X - X(I) ) / ( X(IPOL) - X(I) )
!
!  Modified:
!
!    18 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOL, the number of abscissas.
!    NPOL must be at least 1.
!
!    Input, integer ( kind = 4 ) IPOL, the index of the polynomial to evaluate.
!    IPOL must be between 1 and NPOL.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas of the Lagrange
!    polynomials.  The entries in XPOL must be distinct.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the IPOL-th 
!    Lagrange polynomial is to be evaluated.
!
!    Output, real ( kind = 8 ) PVAL, the value of the IPOL-th Lagrange
!    polynomial at XVAL.
!
!    Output, real ( kind = 8 ) DPDX, the derivative of the IPOL-th 
!    Lagrange polynomial at XVAL.
!
  implicit none

  integer ( kind = 4 ) npol

  real    ( kind = 8 ) dpdx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipol
  integer ( kind = 4 ) j
  real    ( kind = 8 ) p2
  real    ( kind = 8 ) pval
  logical r8vec_distinct
  real    ( kind = 8 ) xpol(npol)
  real    ( kind = 8 ) xval
!
!  Make sure IPOL is legal.
!
  if ( ipol < 1 .or. npol < ipol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY_LAGRANGE_VAL - Fatal error!'
    write ( *, '(a)' ) '  1 <= IPOL <= NPOL is required.'
    write ( *, '(a,i8)' ) '  IPOL = ', ipol
    stop
  end if
!
!  Check that the abscissas are distinct.
!
  if ( .not. r8vec_distinct ( npol, xpol ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY_LAGRANGE_VAL - Fatal error!'
    write ( *, '(a)' ) '  Two or more entries of XPOL are equal:'
    stop
  end if
!
!  Evaluate the polynomial.
!
  pval = 1.0D+00

  do i = 1, npol

    if ( i /= ipol ) then

      pval = pval * ( xval - xpol(i) ) / ( xpol(ipol) - xpol(i) )

    end if

  end do
!
!  Evaluate the derivative, which can be found by summing up the result
!  of differentiating one factor at a time, successively.
!
  dpdx = 0.0D+00

  do i = 1, npol

    if ( i /= ipol ) then

      p2 = 1.0D+00
      do j = 1, npol

        if ( j == i ) then
          p2 = p2                      / ( xpol(ipol) - xpol(j) )
        else if ( j /= ipol ) then
          p2 = p2 * ( xval - xpol(j) ) / ( xpol(ipol) - xpol(j) )
        end if

      end do

      dpdx = dpdx + p2

    end if

  end do

  return
end
subroutine r8poly_order ( na, a, order )

!*****************************************************************************80
!
!! R8POLY_ORDER returns the order of a polynomial.
!
!  Discussion:
!
!    The order of a polynomial is one more than the degree.
!
!    The order of a constant polynomial is 1.  The order of the
!    zero polynomial is debatable, but this routine returns the
!    order as 1.
!
!  Modified:
!
!    19 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:NA), the coefficients of the polynomials.
!
!    Output, integer ( kind = 4 ) ORDER, the order of A.
!
  implicit none

  integer ( kind = 4 ) na

  real    ( kind = 8 ) a(0:na)
  integer ( kind = 4 ) order

  order = na + 1

  do while ( 1 < order )

    if ( a(order-1) /= 0.0D+00 ) then
      return
    end if

    order = order - 1

  end do

  return
end
subroutine r8poly_print ( n, a, title )

!*****************************************************************************80
!
!! R8POLY_PRINT prints out a polynomial.
!
!  Modified:
!
!    30 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X**N.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(0:n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) mag
  integer ( kind = 4 ) n2
  character plus_minus
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call r8poly_degree ( n, a, n2 )

  if ( n2 <= 0 ) then
    write ( *, '( ''  p(x) = 0'' )' )
    return
  end if

  if ( a(n2) < 0.0D+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n2) )

  if ( 2 <= n2 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n2
  else if ( n2 == 1 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) &
      plus_minus, mag
  else if ( n2 == 0 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n2-1, 0, -1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0D+00 ) then

      if ( 2 <= i ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''         '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
subroutine r8poly_shift ( scale, shift, n, poly_cof )

!*****************************************************************************80
!
!! R8POLY_SHIFT adjusts the coefficients of a polynomial for a new argument.
!
!  Discussion:
!
!    Assuming P(X) is a polynomial in the argument X, of the form:
!
!      P(X) =
!          C(N) * X**N
!        + ...
!        + C(1) * X
!        + C(0),
!
!    and that Z is related to X by the formula:
!
!      Z = SCALE * X + SHIFT
!
!    then this routine computes coefficients C for the polynomial Q(Z):
!
!      Q(Z) =
!          C(N) * Z**N
!        + ...
!        + C(1) * Z
!        + C(0)
!
!    so that:
!
!      Q(Z(X)) = P(X)
!
!  Example:
!
!    P(X) = 2 * X**2 - X + 6
!
!    Z = 2.0D+00 * X + 3.0D+00
!
!    Q(Z) = 0.5 *         Z**2 -  3.5 * Z + 12
!
!    Q(Z(X)) = 0.5 * ( 4.0D+00 * X**2 + 12.0D+00 * X +  9 )
!            - 3.5 * (               2.0D+00 * X +  3 )
!                                            + 12
!
!            = 2.0D+00         * X**2 -  1.0D+00 * X +  6
!
!            = P(X)
!
!  Modified:
!
!    05 October 1999
!
!  Reference:
!
!    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
!    Numerical Recipes: The Art of Scientific Computing,
!    Cambridge University Press.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) SHIFT, SCALE, the shift and scale applied to X,
!    so that Z = SCALE * X + SHIFT.
!
!    Input, integer ( kind = 4 ) N, the number of coefficients.
!
!    Input/output, real ( kind = 8 ) POLY_COF(0:N).
!    On input, the coefficient array in terms of the X variable.
!    On output, the coefficient array in terms of the Z variable.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) poly_cof(0:n)
  real    ( kind = 8 ) scale
  real    ( kind = 8 ) shift

  do i = 1, n
    poly_cof(i:n) = poly_cof(i:n) / scale
  end do

  do i = 0, n - 1
    do j = n - 1, i, -1
      poly_cof(j) = poly_cof(j) - shift * poly_cof(j+1)
    end do
  end do

  return
end
subroutine r8poly_val_horner ( n, c, x, cx )

!*****************************************************************************80
!
!! R8POLY_VAL_HORNER evaluates a polynomial using Horner's method.
!
!  Modified:
!
!    08 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of C.
!
!    Input, real ( kind = 8 ) C(0:N), the polynomial coefficients.
!    C(I) is the coefficient of X**I.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial is 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) CX, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) c(0:n)
  real    ( kind = 8 ) cx
  integer ( kind = 4 ) i
  real    ( kind = 8 ) x

  cx = c(n)
  do i = n - 1, 0, -1
    cx = cx * x + c(i)
  end do

  return
end
subroutine r8poly2_ex ( x1, y1, x2, y2, x3, y3, x, y, ierror )

!*****************************************************************************80
!
!! R8POLY2_EX finds the extremal point of a parabola determined by three points.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the coordinates of 
!    three points on the parabola.  X1, X2 and X3 must be distinct.
!
!    Output, real ( kind = 8 ) X, Y, the X coordinate of the extremal point
!    of the parabola, and the value of the parabola at that point.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    1, two of the X values are equal.
!    2, the data lies on a straight line; there is no finite extremal
!    point.
!
  implicit none

  real    ( kind = 8 ) bot
  integer ( kind = 4 ) ierror
  real    ( kind = 8 ) x
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) x2
  real    ( kind = 8 ) x3
  real    ( kind = 8 ) y
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2
  real    ( kind = 8 ) y3

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

  if ( bot == 0.0D+00 ) then
    ierror = 2
    return
  end if

  x = 0.5D+00 * ( &
          x1**2 * ( y3 - y2 ) &
        + x2**2 * ( y1 - y3 ) &
        + x3**2 * ( y2 - y1 ) ) / bot

  y = ( &
         ( x - x2 ) * ( x - x3 ) * ( x2 - x3 ) * y1 &
       - ( x - x1 ) * ( x - x3 ) * ( x1 - x3 ) * y2 &
       + ( x - x1 ) * ( x - x2 ) * ( x1 - x2 ) * y3 ) / &
       ( ( x1 - x2 ) * ( x2 - x3 ) * ( x1 - x3 ) )

  return
end
subroutine r8poly2_ex2 ( x1, y1, x2, y2, x3, y3, x, y, a, b, c, ierror )

!*****************************************************************************80
!
!! R8POLY2_EX2 finds the extremal point of a parabola determined by three points.
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
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the coordinates of 
!    three points on the parabola.  X1, X2 and X3 must be distinct.
!
!    Output, real ( kind = 8 ) X, Y, the X coordinate of the extremal 
!    point of the parabola, and the value of the parabola at that point.
!
!    Output, real ( kind = 8 ) A, B, C, the coefficients that define the
!    parabola: P(X) = A * X**2 + B * X + C.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    1, two of the X values are equal.
!    2, the data lies on a straight line; there is no finite extremal
!    point.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) det
  integer ( kind = 4 ) ierror
  real    ( kind = 8 ) v(3,3)
  real    ( kind = 8 ) w(3,3)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) x2
  real    ( kind = 8 ) x3
  real    ( kind = 8 ) y
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2
  real    ( kind = 8 ) y3

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
  v(1,1) = 1.0D+00
  v(1,2) = x1
  v(1,3) = x1 * x1

  v(2,1) = 1.0D+00
  v(2,2) = x2
  v(2,3) = x2 * x2

  v(3,1) = 1.0D+00
  v(3,2) = x3
  v(3,3) = x3 * x3
!
!  Get the inverse.
!
  call r8mat_inverse_3d ( v, w, det )
!
!  Compute the parabolic coefficients.
!
  c = w(1,1) * y1 + w(1,2) * y2 + w(1,3) * y3
  b = w(2,1) * y1 + w(2,2) * y2 + w(2,3) * y3
  a = w(3,1) * y1 + w(3,2) * y2 + w(3,3) * y3
!
!  Determine the extremal point.
!
  if ( a == 0.0D+00 ) then
    ierror = 2
    return
  end if

  x = -b / ( 2.0D+00 * a )
  y = a * x * x + b * x + c

  return
end
subroutine r8poly2_root ( a, b, c, r1, r2 )

!*****************************************************************************80
!
!! R8POLY2_ROOT returns the two roots of a quadratic polynomial.
!
!  Discussion:
!
!    The polynomial has the form:
!
!      A * X**2 + B * X + C = 0 
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the coefficients of the polynomial.
!    A must not be zero.
!
!    Output, complex ( kind = 8 ) R1, R2, the roots of the polynomial, which
!    might be real and distinct, real and equal, or complex conjugates.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  complex ( kind = 8 ) disc
  complex ( kind = 8 ) q
  complex ( kind = 8 ) r1
  complex ( kind = 8 ) r2

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_ROOT - Fatal error!'
    write ( *, '(a)' ) '  The coefficient A is zero.'
    stop
  end if

  disc = b * b - 4.0D+00 * a * c
  q = -0.5D+00 * ( b + sign ( 1.0D+00, b ) * sqrt ( disc ) )
  r1 = q / a
  r2 = c / q

  return
end
subroutine r8poly2_rroot ( a, b, c, r1, r2 )

!*****************************************************************************80
!
!! R8POLY2_RROOT returns the real parts of the roots of a quadratic polynomial.
!
!  Example:
!
!     A    B    C       roots              R1   R2
!    --   --   --     ------------------   --   --
!     1   -4    3     1          3          1    3
!     1    0    4     2*i      - 2*i        0    0
!     2   -6    5     3 +   i    3 -   i    3    3
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, the coefficients of the quadratic
!    polynomial A * X**2 + B * X + C = 0 whose roots are desired.
!    A must not be zero.
!
!    Output, real ( kind = 8 ) R1, R2, the real parts of the roots
!    of the polynomial.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) disc
  real    ( kind = 8 ) q
  real    ( kind = 8 ) r1
  real    ( kind = 8 ) r2

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_RROOT - Fatal error!'
    write ( *, '(a)' ) '  The coefficient A is zero.'
    stop
  end if

  disc = b * b - 4.0D+00 * a * c
  disc = max ( disc, 0.0D+00 )

  q = ( b + sign ( 1.0D+00, b ) * sqrt ( disc ) )
  r1 = -0.5D+00 * q / a
  r2 = -2.0D+00 * c / q

  return
end
subroutine r8poly2_val ( x1, y1, x2, y2, x3, y3, x, y, yp, ypp )

!*****************************************************************************80
!
!! R8POLY2_VAL evaluates a parabola defined by three data values.
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, three pairs of data.
!    If the X values are distinct, then all the Y values represent
!    actual values of the parabola.
!
!    Three special cases are allowed:
!
!      X1 == X2 /= X3: Y2 is the derivative at X1;
!      X1 /= X2 == X3: Y3 is the derivative at X3;
!      X1 == X2 == X3: Y2 is the derivative at X1, and
!                      Y3 is the second derivative at X1.
!
!    Input, real ( kind = 8 ) X, an abscissa at which the parabola is to be
!    evaluated.
!
!    Output, real ( kind = 8 ) Y, YP, YPP, the values of the parabola and
!    its first and second derivatives at X.
!
  implicit none

  integer ( kind = 4 ) distinct
  real    ( kind = 8 ) dif1
  real    ( kind = 8 ) dif2
  real    ( kind = 8 ) x
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) x2
  real    ( kind = 8 ) x3
  real    ( kind = 8 ) y
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2
  real    ( kind = 8 ) y3
  real    ( kind = 8 ) yp
  real    ( kind = 8 ) ypp
!
!  If any X's are equal, put them and the Y data first.
!
  if ( x1 == x2 .and. x2 == x3 ) then
    distinct = 1
  else if ( x1 == x2 ) then
    distinct = 2
  else if ( x1 == x3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_VAL - Fatal error!'
    write ( *, '(a)' ) '  X1 = X3 =/= X2.'
    write ( *, '(a,g14.6)' ) '  X1 = ', x1
    write ( *, '(a,g14.6)' ) '  X2 = ', x2
    write ( *, '(a,g14.6)' ) '  X3 = ', x3
    stop
  else if ( x2 == x3 ) then
    distinct = 2
    call r8_swap ( x1, x2 )
    call r8_swap ( x2, x3 )
    call r8_swap ( y1, y2 )
    call r8_swap ( y2, y3 )
  else
    distinct = 3
  end if
!
!  Set up the coefficients.
!
  if ( distinct == 1 ) then

    dif1 = y2
    dif2 = 0.5D+00 * y3

  else if ( distinct == 2 ) then

    dif1 = y2
    dif2 = ( ( y3 - y1 ) / ( x3 - x1 ) - y2 ) / ( x3 - x2 )

  else if ( distinct == 3 ) then

    dif1 = ( y2 - y1 ) / ( x2 - x1 )
    dif2 =  ( ( y3 - y1 ) / ( x3 - x1 ) &
            - ( y2 - y1 ) / ( x2 - x1 ) ) / ( x3 - x2 )

  end if
!
!  Evaluate.
!
  y = y1 + ( x - x1 ) * dif1 + ( x - x1 ) * ( x - x2 ) * dif2
  yp = dif1 + ( 2.0D+00 * x - x1 - x2 ) * dif2
  ypp = 2.0D+00 * dif2

  return
end
subroutine r8poly2_val2 ( dim_num, ndata, tdata, ydata, left, tval, yval )

!*****************************************************************************80
!
!! R8POLY2_VAL2 evaluates a parabolic interpolant through tabular data.
!
!  Discussion:
!
!    This routine is a utility routine used by OVERHAUSER_SPLINE_VAL.
!    It constructs the parabolic interpolant through the data in
!    3 consecutive entries of a table and evaluates this interpolant
!    at a given abscissa value.
!
!  Modified:
!
!    10 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of a single data point.
!    DIM_NUM must be at least 1.
!
!    Input, integer ( kind = 4 ) NDATA, the number of data points.
!    NDATA must be at least 3.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissas of the data points.
!    The values in TDATA must be in strictly ascending order.
!
!    Input, real ( kind = 8 ) YDATA(DIM_NUM,NDATA), the data points 
!    corresponding to the abscissas.
!
!    Input, integer ( kind = 4 ) LEFT, the location of the first of the three
!    consecutive data points through which the parabolic interpolant
!    must pass.  1 <= LEFT <= NDATA - 2.
!
!    Input, real ( kind = 8 ) TVAL, the value of T at which the parabolic
!    interpolant is to be evaluated.  Normally, TDATA(1) <= TVAL <= T(NDATA),
!    and the data will be interpolated.  For TVAL outside this range,
!    extrapolation will be used.
!
!    Output, real ( kind = 8 ) YVAL(DIM_NUM), the value of the parabolic
!    interpolant at TVAL.
!
  implicit none

  integer ( kind = 4 ) ndata
  integer ( kind = 4 ) dim_num

  real    ( kind = 8 ) dif1
  real    ( kind = 8 ) dif2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  real    ( kind = 8 ) t1
  real    ( kind = 8 ) t2
  real    ( kind = 8 ) t3
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) ydata(dim_num,ndata)
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2
  real    ( kind = 8 ) y3
  real    ( kind = 8 ) yval(dim_num)
!
!  Check.
!
  if ( left < 1 .or. ndata-2 < left ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  LEFT < 1 or NDATA-2 < LEFT.'
    write ( *, '(a,i8)' ) '  LEFT = ', left
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
    stop
  end if
!
!  Copy out the three abscissas.
!
  t1 = tdata(left)
  t2 = tdata(left+1)
  t3 = tdata(left+2)

  if ( t2 <= t1 .or. t3 <= t2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  T2 <= T1 or T3 <= T2.'
    write ( *, '(a,g14.6)' ) '  T1 = ', t1
    write ( *, '(a,g14.6)' ) '  T2 = ', t2
    write ( *, '(a,g14.6)' ) '  T3 = ', t3
    stop
  end if
!
!  Construct and evaluate a parabolic interpolant for the data
!  in each dimension.
!
  do i = 1, dim_num

    y1 = ydata(i,left)
    y2 = ydata(i,left+1)
    y3 = ydata(i,left+2)

    dif1 = ( y2 - y1 ) / ( t2 - t1 )
    dif2 = ( ( y3 - y1 ) / ( t3 - t1 ) &
           - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 )

    yval(i) = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 )

  end do

  return
end
subroutine r8poly3_root ( a, b, c, d, r1, r2, r3 )

!*****************************************************************************80
!
!! R8POLY3_ROOT returns the three roots of a cubic polynomial.
!
!  Discussion:
!
!    The polynomial has the form
!
!      A * X**3 + B * X**2 + C * X + D = 0
!
!  Modified:
!
!    10 December 2004
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, the coefficients of the polynomial.
!    A must not be zero.
!
!    Output, complex ( kind = 8 ) R1, R2, R3, the roots of the polynomial, which
!    will include at least one real root.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  complex ( kind = 8 ) i
  complex ( kind = 8 ) one
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) q
  real    ( kind = 8 ) r
  complex ( kind = 8 ) r1
  complex ( kind = 8 ) r2
  complex ( kind = 8 ) r3
  real    ( kind = 8 ) s1
  real    ( kind = 8 ) s2
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) theta

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY3_ROOT - Fatal error!'
    write ( *, '(a)' ) '  A must not be zero!'
    stop
  end if

  one = cmplx ( 1.0d+00, 0.0D+00, kind = 8 )
  i = sqrt ( -one )

  q = ( ( b / a )**2 - 3.0D+00 * ( c / a ) ) / 9.0D+00

  r = ( 2.0D+00 * ( b / a )**3 - 9.0D+00 * ( b / a ) * ( c / a ) &
      + 27.0D+00 * ( d / a ) ) / 54.0D+00

  if ( r * r < q * q * q ) then

    theta = acos ( r / sqrt ( q**3 ) )
    r1 = -2.0D+00 * sqrt ( q ) * cos (   theta                  / 3.0D+00 )
    r2 = -2.0D+00 * sqrt ( q ) * cos ( ( theta + 2.0D+00 * pi ) / 3.0D+00 )
    r3 = -2.0D+00 * sqrt ( q ) * cos ( ( theta + 4.0D+00 * pi ) / 3.0D+00 )

  else if ( q * q * q <= r * r ) then

    temp = -r + sqrt ( r**2 - q**3 )
    s1 = sign ( 1.0D+00, temp ) * ( abs ( temp ) )**(1.0D+00/3.0D+00)

    temp = -r - sqrt ( r**2 - q**3 )
    s2 = sign ( 1.0D+00, temp ) * ( abs ( temp ) )**(1.0D+00/3.0D+00)

    r1 = s1 + s2
    r2 = -0.5D+00 * ( s1 + s2 ) + i * 0.5D+00 * sqrt ( 3.0D+00 ) * ( s1 - s2 )
    r3 = -0.5D+00 * ( s1 + s2 ) - i * 0.5D+00 * sqrt ( 3.0D+00 ) * ( s1 - s2 )

  end if

  r1 = r1 - b / ( 3.0D+00 * a )
  r2 = r2 - b / ( 3.0D+00 * a )
  r3 = r3 - b / ( 3.0D+00 * a )

  return
end
subroutine r8poly4_root ( a, b, c, d, e, r1, r2, r3, r4 )

!*****************************************************************************80
!
!! R8POLY4_ROOT returns the four roots of a quartic polynomial.
!
!  Discussion:
!
!    The polynomial has the form:
!
!      A * X**4 + B * X**3 + C * X**2 + D * X + E = 0
!
!  Modified:
!
!    10 December 2004
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, the coefficients of the polynomial.
!    A must not be zero.
!
!    Output, complex ( kind = 8 ) R1, R2, R3, R4, the roots of the polynomial.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) a3
  real    ( kind = 8 ) a4
  real    ( kind = 8 ) b
  real    ( kind = 8 ) b3
  real    ( kind = 8 ) b4
  real    ( kind = 8 ) c
  real    ( kind = 8 ) c3
  real    ( kind = 8 ) c4
  real    ( kind = 8 ) d
  real    ( kind = 8 ) d3
  real    ( kind = 8 ) d4
  real    ( kind = 8 ) e
  complex ( kind = 8 ) p
  complex ( kind = 8 ) q
  complex ( kind = 8 ) r
  complex ( kind = 8 ) r1
  complex ( kind = 8 ) r2
  complex ( kind = 8 ) r3
  complex ( kind = 8 ) r4
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY4_ROOT - Fatal error!'
    write ( *, '(a)') '  A must not be zero!'
    stop
  end if

  a4 = b / a
  b4 = c / a
  c4 = d / a
  d4 = e / a
!
!  Set the coefficients of the resolvent cubic equation.
!
  a3 = 1.0D+00
  b3 = -b4
  c3 = a4 * c4 - 4.0D+00 * d4
  d3 = -a4 * a4 * d4 + 4.0D+00 * b4 * d4 - c4 * c4
!
!  Find the roots of the resolvent cubic.
!
  call r8poly3_root ( a3, b3, c3, d3, r1, r2, r3 )
!
!  Choose one root of the cubic, here R1.
!
!  Set R = sqrt ( 0.25D+00 * A4**2 - B4 + R1 )
!
  r = sqrt ( 0.25D+00 * a4**2 - b4 + r1 )

  if ( r /= zero ) then

    p = sqrt ( 0.75D+00 * a4**2 - r**2 - 2.0D+00 * b4 &
        + 0.25D+00 * ( 4.0D+00 * a4 * b4 - 8.0D+00 * c4 - a4**3 ) / r )

    q = sqrt ( 0.75D+00 * a4**2 - r**2 - 2.0D+00 * b4 &
        - 0.25D+00 * ( 4.0D+00 * a4 * b4 - 8.0D+00 * c4 - a4**3 ) / r )

  else

    p = sqrt ( 0.75D+00 * a4**2 - 2.0D+00 * b4 &
      + 2.0D+00 * sqrt ( r1**2 - 4.0D+00 * d4 ) )

    q = sqrt ( 0.75D+00 * a4**2 - 2.0D+00 * b4 &
      - 2.0D+00 * sqrt ( r1**2 - 4.0D+00 * d4 ) )

  end if
!
!  Set the roots.
!
  r1 = -0.25D+00 * a4 + 0.5D+00 * r + 0.5D+00 * p
  r2 = -0.25D+00 * a4 + 0.5D+00 * r - 0.5D+00 * p
  r3 = -0.25D+00 * a4 - 0.5D+00 * r + 0.5D+00 * q
  r4 = -0.25D+00 * a4 - 0.5D+00 * r - 0.5D+00 * q

  return
end
subroutine r8row_max ( m, n, a, amax )

!*****************************************************************************80
!
!! R8ROW_MAX returns the maximums of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of double precision values, regarded
!    as an array of M rows of length N.
!
!  Example:
!
!    A =
!      1  2  3
!      2  6  7
!
!    MAX =
!      3
!      7
!
!  Modified:
!
!    10 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the array.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) AMAX(M), the maximums of the rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) amax(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m

    amax(i) = a(i,1)
    do j = 2, n
      if ( amax(i) < a(i,j) ) then
        amax(i) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine r8row_mean ( m, n, a, mean )

!*****************************************************************************80
!
!! R8ROW_MEAN returns the means of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of double precision values, regarded
!    as an array of M rows of length N.
!
!  Example:
!
!    A =
!      1  2  3
!      2  6  7
!
!    MEAN =
!      2
!      5
!
!  Modified:
!
!    10 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in the array.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) MEAN(M), the means, or averages, of the rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) mean(m)

  do i = 1, m
    mean(i) = sum ( a(i,1:n) ) / real ( n, kind = 8 )
  end do

  return
end
subroutine r8row_min ( m, n, a, amin )

!*****************************************************************************80
!
!! R8ROW_MIN returns the minimums of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of double precision values, regarded
!    as an array of M rows of length N.
!
!  Example:
!
!    A =
!      1  2  3
!      2  6  7
!
!    MIN =
!      1
!      2
!
!  Modified:
!
!    10 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in the array.
!
!    Input, real ( kind = 8 ) A(M,N), the array to be examined.
!
!    Output, real ( kind = 8 ) AMIN(M), the minimums of the rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) amin(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, m

    amin(i) = a(i,1)
    do j = 2, n
      if ( a(i,j) < amin(i) ) then
        amin(i) = a(i,j)
      end if
    end do

  end do

  return
end
subroutine r8row_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! R8ROW_SORTED_UNIQUE_COUNT counts unique elements in an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of double precision values, regarded
!    as an array of M rows of length N.
!
!    The rows of the array may be ascending or descending sorted.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), a sorted array, containing
!    M rows of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  i1 = 1

  do i2 = 2, m

    if ( any ( a(i1,1:n) /= a(i2,1:n) ) ) then
      unique_num = unique_num + 1
      i1 = i2
    end if

  end do

  return
end
subroutine r8row_sum ( m, n, a, rowsum )

!*****************************************************************************80
!
!! R8ROW_SUM returns the sums of the rows of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of double precision values, regarded
!    as an array of M rows of length N.
!
!  Modified:
!
!    06 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Output, real ( kind = 8 ) ROWSUM(M), the sum of the entries of 
!    each row.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) rowsum(m)

  do i = 1, m
    rowsum(i) = sum ( a(i,1:n) )
  end do

  return
end
subroutine r8row_swap ( m, n, a, i1, i2 )

!*****************************************************************************80
!
!! R8ROW_SWAP swaps two rows of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of double precision values, regarded
!    as an array of M rows of length N.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Input, integer ( kind = 4 ) I1, I2, the two rows to swap.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real    ( kind = 8 ) row(n)

  if ( i1 < 1 .or. m < i1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I1 is out of range.'
    write ( *, '(a,i8)' ) '  I1 = ', i1
    stop
  end if

  if ( i2 < 1 .or. m < i2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I2 is out of range.'
    write ( *, '(a,i8)' ) '  I2 = ', i2
    stop
  end if

  if ( i1 == i2 ) then
    return
  end if

  row(1:n) = a(i1,1:n)
  a(i1,1:n) = a(i2,1:n)
  a(i2,1:n) = row(1:n)

  return
end
subroutine r8row_to_r8vec ( m, n, a, x )

!*****************************************************************************80
!
!! R8ROW_TO_R8VEC converts an R8ROW into an R8VEC.
!
!  Discussion:
!
!    An R8ROW is an M by N array of double precision values, regarded
!    as an array of M rows of length N.
!
!  Example:
!
!    M = 3, N = 4
!
!    A =
!      11 12 13 14
!      21 22 23 24
!      31 32 33 34
!
!    X = ( 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34 )
!
!  Modified:
!
!    13 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), the M by N array.
!
!    Output, real ( kind = 8 ) X(M*N), a vector containing the M rows of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x(m*n)

  j = 1
  do i = 1, m
    x(j:j+n-1) = a(i,1:n)
    j = j + n
  end do

  return
end
subroutine r8row_variance ( m, n, a, variance )

!*****************************************************************************80
!
!! R8ROW_VARIANCE returns the variances of an R8ROW.
!
!  Discussion:
!
!    An R8ROW is an M by N array of double precision values, regarded
!    as an array of M rows of length N.
!
!  Modified:
!
!    10 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in the array.
!
!    Input, real ( kind = 8 ) A(M,N), the array whose variances are desired.
!
!    Output, real ( kind = 8 ) VARIANCE(M), the variances of the rows.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) mean
  real    ( kind = 8 ) variance(m)

  do i = 1, m

    mean = sum ( a(i,1:n) ) / real ( n, kind = 8 )

    variance(i) = 0.0D+00
    do j = 1, n
      variance(i) = variance(i) + ( a(i,j) - mean )**2
    end do

    if ( 1 < n ) then
      variance(i) = variance(i) / real ( n - 1, kind = 8 )
    else
      variance(i) = 0.0D+00
    end if

  end do

  return
end
subroutine r8slmat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8SLMAT_PRINT prints a strict lower triangular R8MAT.
!
!  Example:
!
!    M = 5, N = 5
!    A = (/ 21, 31, 41, 51, 32, 42, 52, 43, 53, 54 /)
!
!    21
!    31 32
!    41 42 43
!    51 52 53 54
!
!  Modified:
!
!    30 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(*), the M by N matrix.  Only the strict
!    lower triangular elements are stored, in column major order.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(10)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) size
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  jmax = min ( n, m-1 )

  if ( m-1 <= n ) then
    size = ( m * ( m - 1 ) ) / 2
  else if ( n < m-1 ) then
    size = ( n * ( n - 1 ) ) / 2 + ( m - n - 1 ) * n
  end if

  if ( all ( a(1:size) == aint ( a(1:size) ) ) ) then

    nn = 10

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m - 1, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(a8,10i8)' ) '     Col', ( j, j = jlo, jhi )
      write ( *, '(a8)' )      '     Row'
      do i = jlo + 1, m
        jhi = min ( jlo + nn - 1, i - 1, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j + 1 ) ) / 2
        end do
        write ( *, '(2x,i8,10i8)' ) i, int ( a(indx(1:jhi+1-jlo)) )
      end do
    end do

  else if ( maxval ( abs ( a(1:size) ) ) < 1000000.0D+00 ) then

    nn = 5

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m - 1, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(a10,5(i8,6x))' ) '       Col', ( j, j = jlo, jhi )
      write ( *, '(a10)' )          '       Row'
      do i = jlo + 1, m
        jhi = min ( jlo + nn - 1, i - 1, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j + 1 ) ) / 2
        end do
        write ( *, '(2x,i8,5f14.6)' ) i, a(indx(1:jhi+1-jlo))
      end do
    end do

  else

    nn = 5

    do jlo = 1, jmax, nn
      jhi = min ( jlo + nn - 1, m - 1, jmax )
      write ( *, '(a)' ) ' '
      write ( *, '(a10,5(i8,6x))' ) '       Col', ( j, j = jlo, jhi )
      write ( *, '(a10)' ) '       Row'
      do i = jlo + 1, m
        jhi = min ( jlo + nn - 1, i - 1, jmax )
        do j = jlo, jhi
          indx(j+1-jlo) = ( j - 1 ) * m + i - ( j * ( j + 1 ) ) / 2
        end do
        write ( *, '(2x,i8,5g14.6)' ) i, a(indx(1:jhi+1-jlo))
      end do
    end do

  end if

  return
end
subroutine r8vec_01_to_ab ( n, a, amax, amin )

!*****************************************************************************80
!
!! R8VEC_01_TO_AB shifts and rescales an R8VEC to lie within given bounds.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    On input, A contains the original data, which is presumed to lie
!    between 0 and 1.  However, it is not necessary that this be so.
!
!    On output, A has been shifted and rescaled so that all entries which
!    on input lay in [0,1] now lie between AMIN and AMAX.  Other entries will
!    be mapped in a corresponding way.
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input/output, real ( kind = 8 ) A(N), the vector to be rescaled.
!
!    Input, real ( kind = 8 ) AMAX, AMIN, the maximum and minimum values 
!    allowed for A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) amax
  real    ( kind = 8 ) amax2
  real    ( kind = 8 ) amax3
  real    ( kind = 8 ) amin
  real    ( kind = 8 ) amin2
  real    ( kind = 8 ) amin3

  if ( amax == amin ) then
    a(1:n) = amin
    return
  end if

  amax2 = max ( amax, amin )
  amin2 = min ( amax, amin )

  amin3 = minval ( a(1:n) )
  amax3 = maxval ( a(1:n) )

  if ( amax3 /= amin3 ) then

    a(1:n) = ( ( amax3 - a(1:n)         ) * amin2   &
             + (         a(1:n) - amin3 ) * amax2 ) &
             / ( amax3          - amin3 )

  else

    a(1:n) = 0.5D+00 * ( amax2 + amin2 )

  end if

  return
end
subroutine r8vec_ab_to_01 ( n, a )

!*****************************************************************************80
!
!! R8VEC_AB_TO_01 shifts and rescales an R8VEC to lie within [0,1].
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    On input, A contains the original data.  On output, A has been shifted
!    and scaled so that all entries lie between 0 and 1.
!
!  Formula:
!
!    A(I) := ( A(I) - AMIN ) / ( AMAX - AMIN )
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input/output, real ( kind = 8 ) A(N), the data to be rescaled.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) amax
  real    ( kind = 8 ) amin

  amax = maxval ( a(1:n) )
  amin = minval ( a(1:n) )

  if ( amin == amax ) then
    a(1:n) = 0.5D+00
  else
    a(1:n) = ( a(1:n) - amin ) / ( amax - amin )
  end if

  return
end
subroutine r8vec_ab_to_cd ( n, a, bmin, bmax, b )

!*****************************************************************************80
!
!! R8VEC_AB_TO_CD shifts and rescales an R8VEC from one interval to another.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The mininum entry of A is mapped to BMIN, the maximum entry
!    to BMAX, and values in between are mapped linearly.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) A(N), the data to be remapped.
!
!    Input, real ( kind = 8 ) BMIN, BMAX, the values to which min(A) and max(A)
!    are to be assigned.
!
!    Output, real ( kind = 8 ) B(N), the remapped data.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) amax
  real    ( kind = 8 ) amin
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) bmax
  real    ( kind = 8 ) bmin

  if ( bmax == bmin ) then
    b(1:n) = bmin
    return
  end if

  amin = minval ( a(1:n) )
  amax = maxval ( a(1:n) )

  if ( amax == amin ) then
    b(1:n) = 0.5D+00 * ( bmax + bmin )
    return
  end if

  b(1:n) = ( ( amax - a(1:n)        ) * bmin   &
         + (          a(1:n) - amin ) * bmax ) &
           / ( amax          - amin )

  return
end
subroutine r8vec_amax ( n, a, amax )

!*****************************************************************************80
!
!! R8VEC_AMAX returns the maximum absolute value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) AMAX, the value of the entry
!    of largest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) amax

  amax = maxval ( abs ( a(1:n) ) )

  return
end
subroutine r8vec_amax_index ( n, a, amax_index )

!*****************************************************************************80
!
!! R8VEC_AMAX_INDEX returns the index of the maximum absolute value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) AMAX_INDEX, the index of the entry of largest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) amax
  integer ( kind = 4 ) amax_index
  integer ( kind = 4 ) i

  if ( n <= 0 ) then

    amax_index = -1

  else

    amax_index = 1
    amax = abs ( a(1) )

    do i = 2, n
      if ( amax < abs ( a(i) ) ) then
        amax_index = i
        amax = abs ( a(i) )
      end if
    end do

  end if

  return
end
subroutine r8vec_amin ( n, a, amin )

!*****************************************************************************80
!
!! R8VEC_AMIN returns the minimum absolute value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 )A(N), the array.
!
!    Output, real ( kind = 8 ) AMIN, the value of the entry 
!    of smallest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) amin

  amin = minval ( abs ( a(1:n) ) )

  return
end
subroutine r8vec_amin_index ( n, a, amin_index )

!*****************************************************************************80
!
!! R8VEC_AMIN_INDEX returns the index of the minimum absolute value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) AMIN_INDEX, the index of the entry of smallest magnitude.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) amin
  integer ( kind = 4 ) amin_index
  integer ( kind = 4 ) i

  if ( n <= 0 ) then

    amin_index = 0

  else

    amin_index = 1
    amin = abs ( a(1) )

    do i = 2, n
      if ( abs ( a(i) ) < amin ) then
        amin_index = i
        amin = abs ( a(i) )
      end if
    end do

  end if

  return
end
function r8vec_ascends ( n, x )

!*****************************************************************************80
!
!! R8VEC_ASCENDS determines if an R8VEC is (weakly) ascending.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    For example, if:
!
!      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
!
!    then
!
!      R8VEC_ASCENDS = TRUE
!
!    The sequence is not required to be strictly ascending.
!
!  Modified:
!
!    26 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, real ( kind = 8 ) X(N), the array to be examined.
!
!    Output, logical R8VEC_ASCENDS, is TRUE if the 
!    entries of X ascend.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical r8vec_ascends
  real    ( kind = 8 ) x(n)

  do i = 1, n - 1
    if ( x(i+1) < x(i) ) then
      r8vec_ascends = .false.
      return
    end if
  end do

  r8vec_ascends = .true.

  return
end
subroutine r8vec_blend ( n, t1, x1, t2, x2, x )

!*****************************************************************************80
!
!! R8VEC_BLEND prforms weighted interpolation of two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The formula used is:
!
!      x(i) = t * x1(i) + (1-t) * x2(i)
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in each  vector.
!
!    Input, real ( kind = 8 ) T1, the weight factor for vector 1.
!
!    Input, real ( kind = 8 ) X1(N), the first vector.
!
!    Input, real ( kind = 8 ) T2, the weight factor for vector 2.
!
!    Input, real ( kind = 8 ) X2(N), the second vector.
!
!    Output, real ( kind = 8 ) X(N), the interpolated or extrapolated value.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) t1
  real    ( kind = 8 ) t2
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) x1(n)
  real    ( kind = 8 ) x2(n)

  x(1:n) = t1 * x1(1:n) + t2 * x2(1:n)

  return
end
subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!    It is always true that RIGHT = LEFT+1.
!
!    If XVAL < X(1), then LEFT = 1, RIGHT = 2, and
!      XVAL   < X(1) < X(2);
!    If X(1) <= XVAL < X(N), then
!      X(LEFT) <= XVAL < X(RIGHT);
!    If X(N) <= XVAL, then LEFT = N-1, RIGHT = N, and
!      X(LEFT) <= X(RIGHT) <= XVAL.
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
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array that has been sorted into 
!    ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xval

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
subroutine r8vec_bracket2 ( n, x, xval, start, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET2 searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    containing the given value.
!
!    R8VEC_BRACKET2 is a variation on R8VEC_BRACKET.  It seeks to reduce
!    the search time by allowing the user to suggest an interval that
!    probably contains the value.  The routine will look in that interval
!    and the intervals to the immediate left and right.  If this does
!    not locate the point, a binary search will be carried out on
!    appropriate subportion of the sorted array.
!
!    In the most common case, 1 <= LEFT < LEFT + 1 = RIGHT <= N,
!    and X(LEFT) <= XVAL <= X(RIGHT).
!
!    Special cases:
!      Value is less than all data values:
!    LEFT = -1, RIGHT = 1, and XVAL < X(RIGHT).
!      Value is greater than all data values:
!    LEFT = N, RIGHT = -1, and X(LEFT) < XVAL.
!      Value is equal to a data value:
!    LEFT = RIGHT, and X(LEFT) = X(RIGHT) = XVAL.
!
!  Modified:
!
!    26 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of the input array.
!
!    Input, real ( kind = 8 ) X(N), an array that has been sorted into 
!    ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed by entries of X.
!
!    Input, integer ( kind = 4 ) START, between 1 and N, specifies that XVAL
!    is likely to be in the interval:
!
!      [ X(START), X(START+1) ]
!
!    or, if not in that interval, then either
!
!      [ X(START+1), X(START+2) ]
!    or
!      [ X(START-1), X(START) ].
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) high
  integer ( kind = 4 ) left
  integer ( kind = 4 ) low
  integer ( kind = 4 ) right
  integer ( kind = 4 ) start
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xval
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_BRACKET2 - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( start < 1 .or. n < start ) then
    start = ( n + 1 ) / 2
  end if
!
!  XVAL = X(START)?
!
  if ( x(start) == xval ) then

    left = start
    right = start
    return
!
!  X(START) < XVAL?
!
  else if ( x(start) < xval ) then
!
!  X(START) = X(N) < XVAL < Infinity?
!
    if ( n < start + 1 ) then

      left = start
      right = -1
      return
!
!  XVAL = X(START+1)?
!
    else if ( xval == x(start+1) ) then

      left = start + 1
      right = start + 1
      return
!
!  X(START) < XVAL < X(START+1)?
!
    else if ( xval < x(start+1) ) then

      left = start
      right = start + 1
      return
!
!  X(START+1) = X(N) < XVAL < Infinity?
!
    else if ( n < start + 2 ) then

      left = start + 1
      right = -1
      return
!
!  XVAL = X(START+2)?
!
    else if ( xval == x(start+2) ) then

      left = start + 2
      right = start + 2
      return
!
!  X(START+1) < XVAL < X(START+2)?
!
    else if ( xval < x(start+2) ) then

      left = start + 1
      right = start + 2
      return
!
!  Binary search for XVAL in [ X(START+2), X(N) ],
!  where XVAL is guaranteed to be greater than X(START+2).
!
    else

      low = start + 2
      high = n
      call r8vec_bracket ( high + 1 - low, x(low), xval, left, right )
      left = left + low - 1
      right = right + low - 1

    end if
!
!  -Infinity < XVAL < X(START) = X(1).
!
  else if ( start == 1 ) then

    left = -1
    right = start
    return
!
!  XVAL = X(START-1)?
!
  else if ( xval == x(start-1) ) then

    left = start - 1
    right = start - 1
    return
!
!  X(START-1) < XVAL < X(START)?
!
  else if ( x(start-1) <= xval ) then

    left = start - 1
    right = start
    return
!
!  Binary search for XVAL in [ X(1), X(START-1) ],
!  where XVAL is guaranteed to be less than X(START-1).
!
  else

    low = 1
    high = start - 1
    call r8vec_bracket ( high + 1 - low, x(1), xval, left, right )

  end if

  return
end
subroutine r8vec_bracket3 ( n, t, tval, left )

!*****************************************************************************80
!
!! R8VEC_BRACKET3 finds the interval containing or nearest a given value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The routine always returns the index LEFT of the sorted array
!    T with the property that either
!    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
!    *  T < T(LEFT) = T(1), or
!    *  T > T(LEFT+1) = T(N).
!
!    The routine is useful for interpolation problems, where
!    the abscissa must be located within an interval of data
!    abscissas for interpolation, or the "nearest" interval
!    to the (extreme) abscissa must be found so that extrapolation
!    can be carried out.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of the input array.
!
!    Input, real ( kind = 8 ) T(N), an array that has been sorted
!    into ascending order.
!
!    Input, real ( kind = 8 ) TVAL, a value to be bracketed by entries of T.
!
!    Input/output, integer ( kind = 4 ) LEFT.
!
!    On input, if 1 <= LEFT <= N-1, LEFT is taken as a suggestion for the
!    interval [ T(LEFT), T(LEFT+1) ] in which TVAL lies.  This interval
!    is searched first, followed by the appropriate interval to the left
!    or right.  After that, a binary search is used.
!
!    On output, LEFT is set so that the interval [ T(LEFT), T(LEFT+1) ]
!    is the closest to TVAL; it either contains TVAL, or else TVAL
!    lies outside the interval [ T(1), T(N) ].
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) high
  integer ( kind = 4 ) left
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid
  real    ( kind = 8 ) t(n)
  real    ( kind = 8 ) tval
!
!  Check the input data.
!
  if ( n < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_BRACKET3 - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 2.'
    stop
  end if
!
!  If LEFT is not between 1 and N-1, set it to the middle value.
!
  if ( left < 1 .or. n - 1 < left ) then
    left = ( n + 1 ) / 2
  end if
!
!  CASE 1: TVAL < T(LEFT):
!  Search for TVAL in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
!
  if ( tval < t(left) ) then

    if ( left == 1 ) then
      return
    else if ( left == 2 ) then
      left = 1
      return
    else if ( t(left-1) <= tval ) then
      left = left - 1
      return
    else if ( tval <= t(2) ) then
      left = 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
!
    low = 2
    high = left - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( t(mid) <= tval ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE2: T(LEFT+1) < TVAL:
!  Search for TVAL in {T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
!
  else if ( t(left+1) < tval ) then

    if ( left == n - 1 ) then
      return
    else if ( left == n - 2 ) then
      left = left + 1
      return
    else if ( tval <= t(left+2) ) then
      left = left + 1
      return
    else if ( t(n-1) <= tval ) then
      left = n - 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = LEFT+2 to N-2.
!
    low = left + 2
    high = n - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( t(mid) <= tval ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE3: T(LEFT) <= TVAL <= T(LEFT+1):
!  T is in [T(LEFT), T(LEFT+1)], as the user said it might be.
!
  else

  end if

  return
end
subroutine r8vec_ceiling ( n, r8vec, ceilingvec )

!*****************************************************************************80
!
!! R8VEC_CEILING rounds "up" (towards +infinity) entries of an R8VEC.
!
!  Examples:
!
!    R8    Value
!
!   -1.1  -1
!   -1.0  -1
!   -0.9   0
!    0.0   0
!    5.0   5
!    5.1   6
!    5.9   6
!    6.0   6
!
!  Modified:
!
!    20 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, real ( kind = 8 ) R8VEC(N), the values to be rounded up.
!
!    Output, integer ( kind = 4 ) CEILINGVEC(N), the rounded values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ceilingvec(n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) r8vec(n)
  integer ( kind = 4 ) value

  do i = 1, n

    value = int ( r8vec(i) )

    if ( real ( value, kind = 8 ) < r8vec(i) ) then
      value = value + 1
    end if

    ceilingvec(i) = value

  end do

  return
end
subroutine r8vec_circular_variance ( n, x, circular_variance )

!*****************************************************************************80
!
!! R8VEC_CIRCULAR_VARIANCE returns the circular variance of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    02 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector whose variance is desired.
!
!    Output, real ( kind = 8 ) CIRCULAR VARIANCE, the circular variance
!    of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) circular_variance
  real    ( kind = 8 ) mean
  real    ( kind = 8 ) x(n)

  call r8vec_mean ( n, x, mean )

  circular_variance = &
      ( sum ( cos ( x(1:n) - mean ) ) )**2 &
    + ( sum ( sin ( x(1:n) - mean ) ) )**2

  circular_variance = sqrt ( circular_variance ) / real ( n, kind = 8 )

  circular_variance = 1.0D+00 - circular_variance

  return
end
subroutine r8vec_compare ( n, a1, a2, isgn )

!*****************************************************************************80
!
!! R8VEC_COMPARE compares two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    Input:
!
!      A1 = ( 2.0, 6.0, 2.0 )
!      A2 = ( 2.0, 8.0, 12.0 )
!
!    Output:
!
!      ISGN = -1
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, A1 < A2,
!     0, A1 = A2,
!    +1, A1 > A2.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) k

  isgn = 0

  k = 1

  do while ( k <= n )

    if ( a1(k) < a2(k) ) then
      isgn = -1
      return
    else if ( a2(k) < a1(k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine r8vec_convolve_circ ( n, x, y, z )

!*****************************************************************************80
!
!! R8VEC_CONVOLVE_CIRC returns the discrete circular convolution of two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The formula used is:
!
!      z(1+m) = xCCy(m) = sum ( 0 <= k <= n-1 ) x(1+k) * y(1+m-k)
!
!    Here, if the index of Y becomes nonpositive, it is "wrapped around"
!    by having N added to it.
!
!    The circular convolution is equivalent to multiplication of Y by a
!    circulant matrix formed from the vector X.
!
!  Example:
!
!    Input:
!
!      X = (/ 1, 2, 3, 4 /)
!      Y = (/ 1, 2, 4, 8 /)
!
!    Output:
!
!      Circulant form:
!
!      Z = ( 1 4 3 2 )   ( 1 )
!          ( 2 1 4 3 )   ( 2 )
!          ( 3 2 1 4 ) * ( 4 )
!          ( 4 3 2 1 )   ( 8 )
!
!      The formula:
!
!      Z = (/ 1*1 + 2*8 + 3*4 + 4*2,
!             1*2 + 2*1 + 3*8 + 4*4,
!             1*4 + 2*2 + 3*1 + 4*8,
!             1*8 + 2*4 + 3*2 + 4*1 /)
!
!      Result:
!
!      Z = (/ 37, 44, 43, 26 /)
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
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the vectors to be convolved.
!
!    Output, real ( kind = 8 ) Z(N), the circular convolution of X and Y.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) m
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) z(n)

  do m = 1, n
    z(m) = dot_product ( x(1:m), y(m:1:-1) ) &
         + dot_product ( x(m+1:n), y(n:m+1:-1) )
  end do

  return
end
subroutine r8vec_copy ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_COPY copies an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    17 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), the vector to be copied.
!
!    Output, real ( kind = 8 ) A2(N), a copy of A1.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)

  a2(1:n) = a1(1:n)

  return
end
subroutine r8vec_cross_3d ( v1, v2, v3 )

!*****************************************************************************80
!
!! R8VEC_CROSS_3D computes the cross product of two R8VEC's in 3D.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
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
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), the two vectors.
!
!    Output, real ( kind = 8 ) V3(3), the cross product vector.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real    ( kind = 8 ) v1(dim_num)
  real    ( kind = 8 ) v2(dim_num)
  real    ( kind = 8 ) v3(dim_num)

  v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

  return
end
subroutine r8vec_cum ( n, a, a_cum )

!*****************************************************************************80
!
!! R8VEC_CUM computes the cumulutive sum of the entries of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Input:
!
!      A = (/ 1.0, 2.0, 3.0, 4.0 /)
!
!    Output:
!
!      A_CUM = (/ 0.0, 1.0, 3.0, 6.0, 10.0 /)
!
!  Modified:
!
!    26 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be summed.
!
!    Output, real ( kind = 8 ) A_CUM(0:N), the cumulative sum of the 
!    entries of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) a_cum(0:n)
  integer ( kind = 4 ) i

  a_cum(0) = 0.0D+00

  do i = 1, n
    a_cum(i) = a_cum(i-1) + a(i)
  end do

  return
end
subroutine r8vec_dif ( n, h, cof )

!*****************************************************************************80
!
!! R8VEC_DIF computes coefficients for estimating the N-th derivative.
!
!  Discussion:
!
!    The routine computes the N+1 coefficients for a centered finite difference
!    estimate of the N-th derivative of a function.
!
!    The estimate has the form
!
!      FDIF(N,X) = Sum (I = 0 to N) COF(I) * F ( X(I) )
!
!    To understand the computation of the coefficients, it is enough
!    to realize that the first difference approximation is
!
!      FDIF(1,X) = F(X+DX) - F(X-DX) ) / (2*DX)
!
!    and that the second difference approximation can be regarded as
!    the first difference approximation repeated:
!
!      FDIF(2,X) = FDIF(1,X+DX) - FDIF(1,X-DX) / (2*DX)
!         = F(X+2*DX) - 2 F(X) + F(X-2*DX) / (4*DX)
!
!    and so on for higher order differences.
!
!    Thus, the next thing to consider is the integer coefficients of
!    the sampled values of F, which are clearly the Pascal coefficients,
!    but with an alternating negative sign.  In particular, if we
!    consider row I of Pascal's triangle to have entries j = 0 through I,
!    then P(I,J) = P(I-1,J-1) - P(I-1,J), where P(*,-1) is taken to be 0,
!    and P(0,0) = 1.
!
!       1
!      -1  1
!       1 -2   1
!      -1  3  -3   1
!       1 -4   6  -4   1
!      -1  5 -10  10  -5  1
!       1 -6  15 -20  15 -6 1
!
!    Next, note that the denominator of the approximation for the
!    N-th derivative will be (2*DX)**N.
!
!    And finally, consider the location of the N+1 sampling
!    points for F:
!
!      X-N*DX, X-(N-2)*DX, X-(N-4)*DX, ..., X+(N-4)*DX, X+(N-2*DX), X+N*DX.
!
!    Thus, a formula for evaluating FDIF(N,X) is
!
!      fdif = 0.0
!      do i = 0, n
!        xi = x + (2*i-n) * h
!        fdif = fdif + cof(i) * f(xi)
!      end do
!
!  Modified:
!
!    17 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the derivative to be approximated.
!    N must be 0 or greater.
!
!    Input, real ( kind = 8 ) H, the half spacing between points. 
!    H must be positive.
!
!    Output, real ( kind = 8 ) COF(0:N), the coefficients needed to approximate
!    the N-th derivative of a function F.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) cof(0:n)
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  if ( n < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_DIF - Fatal error!'
    write ( *, '(a,i8)' ) '  Derivative order N = ', n
    write ( *, '(a)' ) '  but N must be at least 0.'
    stop
  end if

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_DIF - Fatal error!'
    write ( *, '(a,g14.6)' ) '  The half sampling spacing is H = ', h
    write ( *, '(a)' ) '  but H must be positive.'
    stop
  end if

  do i = 0, n

    cof(i) = 1.0D+00

    do j = i - 1, 1, -1
      cof(j) = -cof(j) + cof(j-1)
    end do

    if ( 0 < i ) then
      cof(0) = -cof(0)
    end if

  end do

  cof(0:n) = cof(0:n) / ( 2.0D+00 * h )**n

  return
end
subroutine r8vec_direct_product ( factor_index, factor_order, factor_value, &
  factor_num, point_num, x )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
!
!  Discussion:
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2, 
!      ..., 
!      item JK of 1D rule K.
!
!    In particular, 
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out that task for the abscissas X.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1: 
!      Order = 4
!      X(1:4) = ( 1, 2, 3, 4 )
!
!    Rule 2:
!      Order = 3
!      X(1:3) = ( 10, 20, 30 )
!
!    Rule 3:
!      Order = 2
!      X(1:2) = ( 100, 200 )
!
!    Product Rule:
!      Order = 24
!      X(1:24) = 
!        ( 1, 10, 100 )
!        ( 2, 10, 100 )
!        ( 3, 10, 100 )
!        ( 4, 10, 100 )
!        ( 1, 20, 100 )
!        ( 2, 20, 100 )
!        ( 3, 20, 100 )
!        ( 4, 20, 100 )
!        ( 1, 30, 100 )
!        ( 2, 30, 100 )
!        ( 3, 30, 100 )
!        ( 4, 30, 100 )
!        ( 1, 10, 200 )
!        ( 2, 10, 200 )
!        ( 3, 10, 200 )
!        ( 4, 10, 200 )
!        ( 1, 20, 200 )
!        ( 2, 20, 200 )
!        ( 3, 20, 200 )
!        ( 4, 20, 200 )
!        ( 1, 30, 200 )
!        ( 2, 30, 200 )
!        ( 3, 30, 200 )
!        ( 4, 30, 200 )
!
!  Modified:
!
!    08 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being processed.
!    The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the direct product.
!
!    Input/output, real ( kind = 8 ) X(FACTOR_NUM,POINT_NUM), the elements of the
!    direct product, which are built up gradually.  Before the first call,
!    X might be set to 0.  After each factor has been input, X should
!    have the correct value.
!
!  Local Parameters:
!
!    Local, integer START, the first location of a block of values to set.
!
!    Local, integer CONTIG, the number of consecutive values to set.
!
!    Local, integer SKIP, the distance from the current value of START
!    to the next location of a block of values to set.
!
!    Local, integer REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  real    ( kind = 8 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  real    ( kind = 8 ) x(factor_num,point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      x(factor_index,start:start+contig-1) = factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
subroutine r8vec_direct_product2 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, w )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
!
!  Discussion:
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out the task involving the weights W.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      W(1:4) = ( 2, 3, 5, 7 )
!
!    Rule 2:
!      Order = 3
!      W(1:3) = ( 11, 13, 17 )
!
!    Rule 3:
!      Order = 2
!      W(1:2) = ( 19, 23 )
!
!    Product Rule:
!      Order = 24
!      W(1:24) =
!        ( 2 * 11 * 19 )
!        ( 3 * 11 * 19 )
!        ( 4 * 11 * 19 )
!        ( 7 * 11 * 19 )
!        ( 2 * 13 * 19 )
!        ( 3 * 13 * 19 )
!        ( 5 * 13 * 19 )
!        ( 7 * 13 * 19 )
!        ( 2 * 17 * 19 )
!        ( 3 * 17 * 19 )
!        ( 5 * 17 * 19 )
!        ( 7 * 17 * 19 )
!        ( 2 * 11 * 23 )
!        ( 3 * 11 * 23 )
!        ( 5 * 11 * 23 )
!        ( 7 * 11 * 23 )
!        ( 2 * 13 * 23 )
!        ( 3 * 13 * 23 )
!        ( 5 * 13 * 23 )
!        ( 7 * 13 * 23 )
!        ( 2 * 17 * 23 )
!        ( 3 * 17 * 23 )
!        ( 5 * 17 * 23 )
!        ( 7 * 17 * 23 )
!
!  Modified:
!
!    10 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being processed.
!    The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the direct product.
!
!    Input/output, real ( kind = 8 ) W(POINT_NUM), the elements of the
!    direct product, which are built up gradually.  Before the first call,
!    W should be set to 1.  
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of values to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values to set.
!
!    Local, integer SKIP, the distance from the current value of START
!    to the next location of a block of values to set.
!
!    Local, integer REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  real    ( kind = 8 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  real    ( kind = 8 ) w(point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
function r8vec_distance ( dim_num, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_DISTANCE returns the Euclidean distance between two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    11 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) V1(DIM_NUM), V2(DIM_NUM), the vectors.
!
!    Output, real ( kind = 8 ) R8VEC_DISTANCE, the Euclidean distance
!    between the vectors.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real    ( kind = 8 ) r8vec_distance
  real    ( kind = 8 ) v1(dim_num)
  real    ( kind = 8 ) v2(dim_num)

  r8vec_distance = sqrt ( sum ( ( v1(1:dim_num) - v2(1:dim_num) )**2 ) )

  return
end
function r8vec_distinct ( n, a )

!*****************************************************************************80
!
!! R8VEC_DISTINCT is true if the entries in an R8VEC are distinct.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be checked.
!
!    Output, logical R8VEC_DISTINCT is TRUE if the elements of A are distinct.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical r8vec_distinct

  r8vec_distinct = .false.

  do i = 2, n
    do j = 1, i - 1
      if ( a(i) == a(j) ) then
        return
      end if
    end do
  end do

  r8vec_distinct = .true.

  return
end
function r8vec_dot ( dim_num, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_DOT finds the dot product of a pair of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    In FORTRAN90, the system routine DOT_PRODUCT should be called
!    directly.
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(2), V2(2), the vectors.
!
!    Output, real ( kind = 8 ) R8VEC_DOT, the dot product.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real    ( kind = 8 ) r8vec_dot
  real    ( kind = 8 ) v1(dim_num)
  real    ( kind = 8 ) v2(dim_num)

  r8vec_dot = dot_product ( v1(1:dim_num), v2(1:dim_num) )

  return
end
function r8vec_eq ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_EQ is true if two R8VECs are equal.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), two vectors to compare.
!
!    Output, logical R8VEC_EQ, is TRUE if every pair of elements A1(I) 
!    and A2(I) are equal, and FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)
  logical r8vec_eq

  r8vec_eq = ( all ( a1(1:n) == a2(1:n) ) )

  return
end
subroutine r8vec_even ( n, alo, ahi, a )

!*****************************************************************************80
!
!! R8VEC_EVEN returns N values, evenly spaced between ALO and AHI.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If N is 1, then the midpoint is returned.
!
!    Otherwise, the two endpoints are returned, and N-2 evenly
!    spaced points between them.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) ALO, AHI, the low and high values.
!
!    Output, real ( kind = 8 ) A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) ahi
  real    ( kind = 8 ) alo
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * alo   &
             + real (     i - 1, kind = 8 ) * ahi ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine r8vec_even_select ( n, xlo, xhi, ival, xval )

!*****************************************************************************80
!
!! R8VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / real ( N - 1 )
!
!    Unless N = 1, X(1) = XLO and X(N) = XHI.
!
!    If N = 1, then X(1) = 0.5*(XLO+XHI).
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) XLO, XHI, the low and high values.
!
!    Input, integer ( kind = 4 ) IVAL, the index of the desired point.
!    IVAL is normally between 1 and N, but may be any integer value.
!
!    Output, real ( kind = 8 ) XVAL, the IVAL-th of N evenly spaced values
!    between XLO and XHI.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ival
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo
  real    ( kind = 8 ) xval

  if ( n == 1 ) then

    xval = 0.5D+00 * ( xlo + xhi )

  else

    xval = ( real ( n - ival,     kind = 8 ) * xlo   &
           + real (     ival - 1, kind = 8 ) * xhi ) &
           / real ( n        - 1, kind = 8 )

  end if

  return
end
subroutine r8vec_even2 ( maxval, nfill, nold, xold, nval, xval )

!*****************************************************************************80
!
!! R8VEC_EVEN2 linearly interpolates new numbers into an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The number of values created between two old values can vary from
!    one pair of values to the next.
!
!    The interpolated values are evenly spaced.
!
!    This routine is a generalization of R8VEC_EVEN.
!
!  Modified:
!
!    30 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXVAL, the size of the XVAL array, as declared by the
!    user.  MAXVAL must be large enough to hold the NVAL values computed by
!    this routine.  In other words, MAXVAL must be at least equal to
!    NOLD + SUM (1 <= I <= NOLD-1) NFILL(I).
!
!    Input, integer ( kind = 4 ) NFILL(NOLD-1), the number of values
!    to be interpolated between XOLD(I) and XOLD(I+1).
!    NFILL(I) does not count the endpoints.  Thus, if
!    NFILL(I) is 1, there will be one new point generated
!    between XOLD(I) and XOLD(I+1).
!    NFILL(I) must be nonnegative.
!
!    Input, integer ( kind = 4 ) NOLD, the number of values XOLD,
!    between which extra values are to be interpolated.
!
!    Input, real ( kind = 8 ) XOLD(NOLD), the original vector of numbers
!    between which new values are to be interpolated.
!
!    Output, integer ( kind = 4 ) NVAL, the number of values computed
!    in the XVAL array.
!    NVAL = NOLD + SUM ( 1 <= I <= NOLD-1 ) NFILL(I)
!
!    Output, real ( kind = 8 ) XVAL(MAXVAL).  On output, XVAL contains the
!    NOLD values of XOLD, as well as the interpolated
!    values, making a total of NVAL values.
!
  implicit none

  integer ( kind = 4 ) maxval
  integer ( kind = 4 ) nold

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nadd
  integer ( kind = 4 ) nfill(nold-1)
  integer ( kind = 4 ) nval
  real    ( kind = 8 ) xold(nold)
  real    ( kind = 8 ) xval(maxval)

  nval = 1

  do i = 1, nold - 1

    if ( nfill(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_EVEN2 - Fatal error!'
      write ( *, '(a,i8)' ) '  NFILL(I) is negative for I = ', i
      write ( *, '(a,i8)' ) '  NFILL(I) = ', nfill(i)
      stop
    end if

    if ( maxval < nval + nfill(i) + 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_EVEN2 - Fatal error!'
      write ( *, '(a)' ) '  MAXVAL is not large enough.  '
      write ( *, '(a,i8)' ) '  MAXVAL = ', maxval
      write ( *, '(a)' ) '  which is exceeded by storage requirements'
      write ( *, '(a,i8)' ) '  for interpolating in interval ', i
      stop
    end if

    nadd = nfill(i) + 2

    do j = 1, nadd
      xval(nval+j-1) = ( real ( nadd - j,     kind = 8 ) * xold(i)   &
                       + real (        j - 1, kind = 8 ) * xold(i+1) ) &
                       / real ( nadd     - 1, kind = 8 )
    end do

    nval = nval + nfill(i) + 1

  end do

  return
end
subroutine r8vec_even3 ( nold, nval, xold, xval )

!*****************************************************************************80
!
!! R8VEC_EVEN3 evenly interpolates new data into an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    This routine accepts a short vector of numbers, and returns a longer
!    vector of numbers, created by interpolating new values between
!    the given values.
!
!    Between any two original values, new values are evenly interpolated.
!
!    Over the whole vector, the new numbers are interpolated in
!    such a way as to try to minimize the largest distance interval size.
!
!    The algorithm employed is not "perfect".
!
!  Modified:
!
!    29 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NOLD, the number of values XOLD, between which extra
!    values are to be interpolated.
!
!    Input, integer ( kind = 4 ) NVAL, the number of values to be computed
!    in the XVAL array.  NVAL should be at least NOLD.
!
!    Input, real ( kind = 8 ) XOLD(NOLD), the original vector of numbers
!    between which new values are to be interpolated.
!
!    Output, real ( kind = 8 ) XVAL(NVAL).  On output, XVAL contains the
!    NOLD values of XOLD, as well as interpolated
!    values, making a total of NVAL values.
!
  implicit none

  integer ( kind = 4 ) nval
  integer ( kind = 4 ) nold

  real    ( kind = 8 ) density
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nmaybe
  integer ( kind = 4 ) npts
  integer ( kind = 4 ) ntemp
  integer ( kind = 4 ) ntot
  real    ( kind = 8 ) xlen
  real    ( kind = 8 ) xleni
  real    ( kind = 8 ) xlentot
  real    ( kind = 8 ) xold(nold)
  real    ( kind = 8 ) xval(nval)

  xlen = 0.0D+00
  do i = 1, nold - 1
    xlen = xlen + abs ( xold(i+1) - xold(i) )
  end do

  ntemp = nval - nold

  density = real ( ntemp, kind = 8 ) / xlen

  ival = 1
  ntot = 0
  xlentot = 0.0D+00

  do i = 1, nold - 1

    xleni = abs ( xold(i+1) - xold(i) )
    npts = int ( density * xleni )
    ntot = ntot + npts
!
!  Determine if we have enough left-over density that it should
!  be changed into a point.  A better algorithm would agonize
!  more over where that point should go.
!
    xlentot = xlentot + xleni
    nmaybe = nint ( xlentot * density )

    if ( ntot < nmaybe ) then
      npts = npts + nmaybe - ntot
      ntot = nmaybe
    end if

    do j = 1, npts + 2
      xval(ival+j-1) = ( real ( npts+2 - j,     kind = 8 ) * xold(i)   &
                       + real (          j - 1, kind = 8 ) * xold(i+1) ) &
                       / real ( npts+2     - 1, kind = 8 )
    end do

    ival = ival + npts + 1

  end do

  return
end
subroutine r8vec_expand_linear ( n, x, fat, xfat )

!*****************************************************************************80
!
!! R8VEC_EXPAND_LINEAR linearly interpolates new data into an R8VEC.
!
!  Discussion:
!
!    This routine copies the old data, and inserts NFAT new values
!    between each pair of old data values.  This would be one way to
!    determine places to evenly sample a curve, given the (unevenly
!    spaced) points at which it was interpolated.
!
!    An R8VEC is an array of double precision real values.
!
!  Example:
!
!    N = 3
!    NFAT = 2
!
!    X(1:N)        = (/ 0.0,           6.0,             7.0 /)
!    XFAT(1:2*3+1) = (/ 0.0, 2.0, 4.0, 6.0, 6.33, 6.66, 7.0 /)
!
!  Modified:
!
!    10 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of input data values.
!
!    Input, real ( kind = 8 ) X(N), the original data.
!
!    Input, integer ( kind = 4 ) FAT, the number of data values to interpolate
!    between each pair of original data values.
!
!    Output, real ( kind = 8 ) XFAT((N-1)*(FAT+1)+1), the "fattened" data.
!
  implicit none

  integer ( kind = 4 ) fat
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xfat((n-1)*(fat+1)+1)

  k = 0

  do i = 1, n - 1

    k = k + 1
    xfat(k) = x(i)

    do j = 1, fat
      k = k + 1
      xfat(k) = ( real ( fat - j + 1, kind = 8 ) * x(i)     &
                + real (       j,     kind = 8 ) * x(i+1) ) &
                / real ( fat     + 1, kind = 8 )
    end do

  end do

  k = k + 1
  xfat(k) = x(n)

  return
end
subroutine r8vec_floor ( n, r8vec, floorvec )

!*****************************************************************************80
!
!! R8VEC_FLOOR rounds "down" (towards -infinity) entries of an R8VEC.
!
!  Examples:
!
!    R8    Value
!
!   -1.1  -2
!   -1.0  -1
!   -0.9  -1
!    0.0   0
!    5.0   5
!    5.1   5
!    5.9   5
!    6.0   6
!
!  Modified:
!
!    20 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, real ( kind = 8 ) R8VEC(N), the values to be rounded down.
!
!    Output, integer ( kind = 4 ) FLOORVEC(N), the rounded value.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) floorvec(n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) r8vec(n)
  integer ( kind = 4 ) value

  do i = 1, n

    value = int ( r8vec(i) )

    if ( r8vec(i) < real ( value, kind = 8 ) ) then
      value = value - 1
    end if

    floorvec(i) = value

  end do

  return
end
subroutine r8vec_frac ( n, a, k, frac )

!*****************************************************************************80
!
!! R8VEC_FRAC searches for the K-th smallest entry in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Hoare's algorithm is used.
!
!  Modified:
!
!    17 July 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, A is the array to search.
!    On output, the elements of A have been somewhat rearranged.
!
!    Input, integer ( kind = 4 ) K, the fractile to be sought.  If K = 1, the minimum
!    entry is sought.  If K = N, the maximum is sought.  Other values
!    of K search for the entry which is K-th in size.  K must be at
!    least 1, and no greater than N.
!
!    Output, real ( kind = 8 ) FRAC, the value of the K-th fractile of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) frac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iryt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) left
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) x

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_FRAC - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal nonpositive value of N = ', n
    stop
  end if

  if ( k <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_FRAC - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal nonpositive value of K = ', k
    stop
  end if

  if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_FRAC - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal N < K, K = ', k
    stop
  end if

  left = 1
  iryt = n

  do

    if ( iryt <= left ) then
      frac = a(k)
      exit
    end if

    x = a(k)
    i = left
    j = iryt

    do

      if ( j < i ) then
        if ( j < k ) then
          left = i
        end if
        if ( k < i ) then
          iryt = j
        end if
        exit
      end if
!
!  Find I so that X <= A(I).
!
      do while ( a(i) < x )
        i = i + 1
      end do
!
!  Find J so that A(J) <= X.
!
      do while ( x < a(j) )
        j = j - 1
      end do

      if ( i <= j ) then

        temp = a(i)
        a(i) = a(j)
        a(j) = temp

        i = i + 1
        j = j - 1
      end if

    end do

  end do

  return
end
subroutine r8vec_fraction ( n, x, fraction )

!*****************************************************************************80
!
!! R8VEC_FRACTION returns the fraction parts of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If we regard a real number as
!
!      R8 = SIGN * ( WHOLE + FRACTION )
!
!    where
!
!      SIGN is +1 or -1,
!      WHOLE is a nonnegative integer
!      FRACTION is a nonnegative real number strictly less than 1,
!
!    then this routine returns the value of FRACTION.
!
!  Example:
!
!     R8    R8_FRACTION
!
!    0.00      0.00
!    1.01      0.01
!    2.02      0.02
!   19.73      0.73
!   -4.34      0.34
!
!  Modified:
!
!    18 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(N), the arguments.
!
!    Output, real ( kind = 8 ) FRACTION(N), the fraction parts.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) fraction(n)
  real    ( kind = 8 ) x(n)

  fraction(1:n) = abs ( x(1:n) ) - real ( int ( abs ( x(1:n) ) ), kind = 8 )

  return
end
function r8vec_gt ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_GT == ( A1 > A2 ) for R8VEC's.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
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
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, logical R8VEC_GT, is TRUE if and only if A1 > A2.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)
  logical r8vec_gt
  integer ( kind = 4 ) i

  r8vec_gt = .false.

  do i = 1, n

    if ( a2(i) < a1(i) ) then
      r8vec_gt = .true.
      exit
    else if ( a1(i) < a2(i) ) then
      r8vec_gt = .false.
      exit
    end if

  end do

  return
end
subroutine r8vec_heap_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_HEAP_A reorders an R8VEC into an ascending heap.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    An ascending heap is an array A with the property that, for every index J,
!    A(J) <= A(2*J) and A(J) <= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Modified:
!
!    07 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  real    ( kind = 8 ) key
  integer ( kind = 4 ) m
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
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the smaller of the two values,
!  and update M if necessary.
!
        if ( a(m+1) < a(m) ) then
          m = m + 1
        end if

      end if
!
!  If the small descendant is smaller than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( key <= a(m) ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot.
!
    a(ifree) = key

  end do

  return
end
subroutine r8vec_heap_d ( n, a )

!*****************************************************************************80
!
!! R8VEC_HEAP_D reorders an R8VEC into an descending heap.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Modified:
!
!    07 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  real    ( kind = 8 ) key
  integer ( kind = 4 ) m
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
      if ( n < m ) then
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
        if ( a(m) < a(m+1) ) then
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
subroutine r8vec_histogram ( n, a, a_lo, a_hi, histo_num, histo_gram )

!*****************************************************************************80
!
!! R8VEC_HISTOGRAM histograms an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Values between A_LO and A_HI will be histogrammed into the bins
!    1 through HISTO_NUM.  Values below A_LO are counted in bin 0,
!    and values greater than A_HI are counted in bin HISTO_NUM+1.
!
!  Modified:
!
!    09 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), the array to examine.  
!
!    Input, real ( kind = 8 ) A_LO, A_HI, the lowest and highest
!    values to be histogrammed.  These values will also define the bins.
!
!    Input, integer ( kind = 4 ) HISTO_NUM, the number of bins to use.
!
!    Output, integer ( kind = 4 ) HISTO_GRAM(0:HISTO_NUM+1), contains the number of
!    entries of A in each bin.
!
  implicit none

  integer ( kind = 4 ) histo_num
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) a_hi
  real    ( kind = 8 ) a_lo
  real    ( kind = 8 ) delta
  integer ( kind = 4 ) histo_gram(0:histo_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  histo_gram(0:histo_num+1) = 0

  delta = ( a_hi - a_lo ) / real ( 2 * histo_num, kind = 8 )

  do i = 1, n

    if ( a(i) < a_lo ) then

      histo_gram(0) = histo_gram(0) + 1

    else if ( a(i) <= a_hi ) then

      j = nint ( &
        ( ( a_hi -           delta - a(i)        ) * real ( 1,         kind = 8 )   &
        + (      -           delta + a(i) - a_lo ) * real ( histo_num, kind = 8 ) ) &
        / ( a_hi - 2.0D+00 * delta        - a_lo ) )

      histo_gram(j) = histo_gram(j) + 1

    else if ( a_hi < a(i) ) then

      histo_gram(histo_num+1) = histo_gram(histo_num+1) + 1

    end if

  end do

  return
end
subroutine r8vec_house_column ( n, a, k, v )

!*****************************************************************************80
!
!! R8VEC_HOUSE_COLUMN defines a Householder premultiplier that "packs" a column.
!
!  Discussion:
!
!    The routine returns a vector V that defines a Householder
!    premultiplier matrix H(V) that zeros out the subdiagonal entries of
!    column K of the matrix A.
!
!       H(V) = I - 2 * v * v'
!
!  Modified:
!
!    01 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, real ( kind = 8 ) A(N), column K of the matrix A.
!
!    Input, integer ( kind = 4 ) K, the column of the matrix to be modified.
!
!    Output, real ( kind = 8 ) V(N), a vector of unit L2 norm which defines an
!    orthogonal Householder premultiplier matrix H with the property
!    that the K-th column of H*A is zero below the diagonal.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) k
  real    ( kind = 8 ) s
  real    ( kind = 8 ) v(n)

  v(1:n) = 0.0D+00

  if ( k < 1 .or. n <= k ) then
    return
  end if

  s = sqrt ( dot_product ( a(k:n), a(k:n) ) )

  if ( s == 0.0D+00 ) then
    return
  end if

  v(k) = a(k) + sign ( s, a(k) )
  v(k+1:n) = a(k+1:n)

  v(k:n) = v(k:n) / sqrt ( dot_product ( v(k:n), v(k:n) ) )

  return
end
function r8vec_in_01 ( n, a )

!*****************************************************************************80
!
!! R8VEC_IN_01 is TRUE if the entries of an R8VEC are in the range [0,1].
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    06 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, logical R8VEC_IN_01, is TRUE if every entry of A is
!    between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  logical r8vec_in_01

  if ( any ( a(1:n) < 0.0D+00 .or. 1.0D+00 < a(1:n) ) ) then
    r8vec_in_01 = .false.
  else
    r8vec_in_01 = .true.
  end if

  return
end
subroutine r8vec_index_delete_all ( n, x, indx, xval )

!*****************************************************************************80
!
!! R8VEC_INDEX_DELETE_ALL deletes all occurrences of a value from an indexed sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Note that the value of N is adjusted because of the deletions!
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, real ( kind = 8 ) X(N), the list.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, the value to be sought.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) equal1
  integer ( kind = 4 ) equal2
  integer ( kind = 4 ) get
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) put
  real    ( kind = 8 ) x(*)
  real    ( kind = 8 ) xval

  if ( n < 1 ) then
    n = 0
    return
  end if

  call r8vec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal == 0 ) then
    return
  end if

  equal1 = equal

  do

    if ( equal1 <= 1 ) then
      exit
    end if

    if ( x(indx(equal1-1)) /= xval ) then
      exit
    end if

    equal1 = equal1 - 1

  end do

  equal2 = equal

  do

    if ( n <= equal2 ) then
      exit
    end if

    if ( x(indx(equal2+1)) /= xval ) then
      exit
    end if

    equal2 = equal2 + 1

  end do
!
!  Discard certain X values.
!
  put = 0

  do get = 1, n

    if ( x(get) /= xval ) then
      put = put + 1
      x(put) = x(get)
    end if

  end do

  x(put+1:n) = 0.0D+00
!
!  Adjust the INDX values.
!
  do equal = equal1, equal2
    do i = 1, n
      if ( indx(equal) < indx(i) ) then
        indx(i) = indx(i) - 1
      end if
    end do
  end do
!
!  Discard certain INDX values.
!
  indx(equal1:n+equal1-equal2-1) = indx(equal2+1:n)
  indx(n+equal1-equal2:n) = 0
!
!  Adjust N.
!
  n = put

  return
end
subroutine r8vec_index_delete_dupes ( n, x, indx, n2, x2, indx2 )

!*****************************************************************************80
!
!! R8VEC_INDEX_DELETE_DUPES deletes duplicates from an indexed sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The output quantities N2, X2, and INDX2 are computed from the
!    input quantities by sorting, and eliminating duplicates.
!
!    The output arrays should be dimensioned of size N, unless the user
!    knows in advance what the value of N2 will be.
!
!    The output arrays may be identified with the input arrays.
!
!  Modified:
!
!    15 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input list.
!
!    Input, real ( kind = 8 ) X(N), the list.  
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Output, integer ( kind = 4 ) N2, the number of unique entries in X.
!
!    Output, real ( kind = 8 ) X2(N2), a copy of the list which has
!    been sorted, and made unique.
!
!    Output, integer ( kind = 4 ) INDX2(N2), the sort index of the new list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indx2(n)
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) x2(n)
  real    ( kind = 8 ) x3(n)

  i = 0
  n3 = 0

  do

    i = i + 1

    if ( n < i ) then
      exit
    end if

    if ( 1 < i ) then
      if ( x(indx(i)) == x3(n3) ) then
        cycle
      end if
    end if

    n3 = n3 + 1
    x3(n3) = x(indx(i))

  end do
!
!  Copy data into output arrays.
!
  n2 = n3
  x2(1:n2) = x3(1:n3)
  call i4vec_indicator ( n2, indx2 )

  return
end
subroutine r8vec_index_delete_one ( n, x, indx, xval, n2, x2, indx2 )

!*****************************************************************************80
!
!! R8VEC_INDEX_DELETE_ONE deletes one copy of a value from an indexed sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the value occurs in the list more than once, only one copy is deleted.
!
!    Note that the value of N is adjusted because of the deletions.
!
!  Modified:
!
!    24 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, real ( kind = 8 ) X(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) N2, the size of the current list.
!
!    Output, real ( kind = 8 ) X2(N2), the list.
!
!    Output, integer ( kind = 4 ) INDX2(N2), the sort index of the list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indx2(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n2
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) x2(n)
  real    ( kind = 8 ) xval

  if ( n < 1 ) then
    n2 = 0
    return
  end if

  n2 = n
  indx2(1:n2) = indx(1:n2)
  x2(1:n2) = x(1:n2)

  call r8vec_index_search ( n2, x2, indx2, xval, less, equal, more )

  if ( equal /= 0 ) then
    j = indx2(equal)
    x2(j:n2-1) = x2(j+1:n2)
    indx2(equal:n2-1) = indx2(equal+1:n2)
    do i = 1, n2-1
      if ( j < indx2(i) ) then
        indx2(i) = indx2(i) - 1
      end if
    end do
    n2 = n2 - 1
  end if

  return
end
subroutine r8vec_index_insert ( n, x, indx, xval )

!*****************************************************************************80
!
!! R8VEC_INDEX_INSERT inserts a value in an indexed sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, real ( kind = 8 ) X(N), the list.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, the value to be sought.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  real    ( kind = 8 ) x(*)
  real    ( kind = 8 ) xval

  if ( n <= 0 ) then
    n = 1
    x(1) = xval
    indx(1) = 1
    return
  end if

  call r8vec_index_search ( n, x, indx, xval, less, equal, more )

  x(n+1) = xval
  indx(n+1:more+1:-1) = indx(n:more:-1)
  indx(more) = n + 1
  n = n + 1

  return
end
subroutine r8vec_index_insert_unique ( n, x, indx, xval )

!*****************************************************************************80
!
!! R8VEC_INDEX_INSERT_UNIQUE inserts a unique value in an indexed sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the value does not occur in the list, it is included in the list,
!    and N, X and INDX are updated.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, real ( kind = 8 ) X(N), the list.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, the value to be sought.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  real    ( kind = 8 ) x(*)
  real    ( kind = 8 ) xval

  if ( n <= 0 ) then
    n = 1
    x(1) = xval
    indx(1) = 1
    return
  end if
!
!  Does XVAL already occur in X?
!
  call r8vec_index_search ( n, x, indx, xval, less, equal, more )

  if ( equal == 0 ) then
    x(n+1) = xval
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1
  end if

  return
end
subroutine r8vec_index_order ( n, x, indx )

!*****************************************************************************80
!
!! R8VEC_INDEX_ORDER sorts an R8VEC using an index vector.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The index vector itself is not modified.  Therefore, the pair
!    (X,INDX) no longer represents an index sorted vector.  If this
!    relationship is to be preserved, then simply set INDX(1:N)=(1:N).
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, real ( kind = 8 ) X(N), the list.  On output, the list
!    has been sorted.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) indx(n)
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) y(n)

  y(1:n) = x(indx(1:n))
  x(1:n) = y(1:n)

  return
end
subroutine r8vec_index_search ( n, x, indx, xval, less, equal, more )

!*****************************************************************************80
!
!! R8VEC_INDEX_SEARCH searches for a value in an indexed sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, real ( kind = 8 ) X(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) LESS, EQUAL, MORE, the indexes in INDX of the
!    entries of X that are just less than, equal to, and just greater
!    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
!    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
!    is the greatest entry of X, then MORE is N+1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) more
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xhi
  real    ( kind = 8 ) xlo
  real    ( kind = 8 ) xmid
  real    ( kind = 8 ) xval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n
  xlo = x(indx(lo))
  xhi = x(indx(hi))

  if ( xval < xlo ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( xval == xlo ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  if ( xhi < xval ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( xval == xhi ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))

    if ( xval == xmid ) then
      equal = mid
      less = equal - 1
      more = equal + 1
      return
    else if ( xval < xmid ) then
      hi = mid
    else if ( xmid < xval ) then
      lo = mid
    end if

  end do

  return
end
subroutine r8vec_index_sort_unique ( n, x, indx, n2 )

!*****************************************************************************80
!
!! R8VEC_INDEX_SORT_UNIQUE creates a sort index for an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input/output, real ( kind = 8 ) X(N), the list.  On output, X contains only
!    unique elements.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Output, integer ( kind = 4 ) N2, the number of unique elements in X.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) n2
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) y(n)

  n2 = 0

  do i = 1, n
    call r8vec_index_insert_unique ( n2, y, indx, x(i) )
  end do

  x(1:n2) = y(1:n2)

  x(n2+1:n) = 0.0D+00
  indx(n2+1:n) = 0

  return
end
subroutine r8vec_indicator ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR sets an R8VEC to the indicator vector.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, real ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = real ( i, kind = 8 )
  end do

  return
end
subroutine r8vec_insert ( n, a, pos, value )

!*****************************************************************************80
!
!! R8VEC_INSERT inserts a value into an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    17 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the array on input.
!
!    Input/output, real ( kind = 8 ) A(N+1), the array.  On input, A is 
!    assumed to contain only N entries, while on output, A actually 
!    contains N+1 entries.
!
!    Input, integer ( kind = 4 ) POS, the position to be assigned the new entry.
!    1 <= POS <= N+1.
!
!    Input, real ( kind = 8 ) VALUE, the value to be inserted.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) pos
  real    ( kind = 8 ) value

  if ( pos < 1 .or. n + 1 < pos ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_INSERT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal insertion position = ', pos
    stop

  else

    do i = n + 1, pos + 1, -1
      a(i) = a(i-1)
    end do

    a(pos) = value

  end if

  return
end
function r8vec_is_int ( n, a )

!*****************************************************************************80
!
!! R8VEC_IS_INT is TRUE if the entries of an R8VEC are integers.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    04 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, logical R8VEC_IS_INT, is TRUE if every entry of A is
!    integral.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  logical r8vec_is_int

  r8vec_is_int = all ( a(1:n) == aint ( a(1:n) ) )

  return
end
function r8vec_length ( dim_num, x )

!*****************************************************************************80
!
!! R8VEC_LENGTH returns the Euclidean length of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    08 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the vector.
!
!    Output, real ( kind = 8 ) R8VEC_LENGTH, the Euclidean length of the vector.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real    ( kind = 8 ) r8vec_length
  real    ( kind = 8 ) x(dim_num)

  r8vec_length = sqrt ( sum ( ( x(1:dim_num) )**2 ) )

  return
end
function r8vec_lt ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_LT == ( A1 < A2 ) for R8VEC's.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
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
!    05 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, logical R8VEC_LT, is TRUE if and only if A1 < A2.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)
  logical r8vec_lt
  integer ( kind = 4 ) i

  r8vec_lt = .false.

  do i = 1, n

    if ( a1(i) < a2(i) ) then
      r8vec_lt = .true.
      exit
    else if ( a2(i) < a1(i) ) then
      r8vec_lt = .false.
      exit
    end if

  end do

  return
end
subroutine r8vec_mask_print ( n, a, mask_num, mask, title )

!*****************************************************************************80
!
!! R8VEC_MASK_PRINT prints a masked R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    24 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MASK_NUM, the number of masked elements.
!
!    Input, integer ( kind = 4 ) MASK(MASK_NUM), the indices of the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) mask_num
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mask(mask_num)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Masked vector printout:'

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, mask_num
    write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) i, mask(i), a(mask(i))
  end do

  return
end
subroutine r8vec_max ( n, a, amax )

!*****************************************************************************80
!
!! R8VEC_MAX returns the maximum value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) AMAX, the value of the largest entry.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) amax

  amax = maxval ( a(1:n) )

  return
end
subroutine r8vec_max_index ( n, a, max_index )

!*****************************************************************************80
!
!! R8VEC_MAX_INDEX returns the index of the maximum value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    02 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) MAX_INDEX, the index of the largest entry.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_index

  if ( n <= 0 ) then

    max_index = -1

  else

    max_index = 1

    do i = 2, n
      if ( a(max_index) < a(i) ) then
        max_index = i
      end if
    end do

  end if

  return
end
subroutine r8vec_mean ( n, a, mean )

!*****************************************************************************80
!
!! R8VEC_MEAN returns the mean of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    02 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector whose mean is desired.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) mean

  mean = sum ( a(1:n) ) / real ( n, kind = 8 )

  return
end
subroutine r8vec_median ( n, a, median )

!*****************************************************************************80
!
!! R8VEC_MEDIAN returns the median of an unsorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Hoare's algorithm is used.  The values of the vector are
!    rearranged by this routine.
!
!  Modified:
!
!    18 September 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, real ( kind = 8 ) A(N), the array to search.  On output,
!    the order of the elements of A has been somewhat changed.
!
!    Output, real ( kind = 8 ) MEDIAN, the value of the median of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) k
  real    ( kind = 8 ) median

  k = ( n + 1 ) / 2

  call r8vec_frac ( n, a, k, median )

  return
end
subroutine r8vec_min ( n, a, amin )

!*****************************************************************************80
!
!! R8VEC_MIN returns the minimum value of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) AMIN, the value of the smallest entry.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) amin

  amin = minval ( a(1:n) )

  return
end
subroutine r8vec_min_index ( n, a, min_index )

!*****************************************************************************80
!
!! R8VEC_MIN_INDEX returns the index of the minimum value in an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    02 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) MIN_INDEX, the index of the smallest entry.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) min_index

  if ( n <= 0 ) then

    min_index = -1

  else

    min_index = 1

    do i = 2, n
      if ( a(i) < a(min_index) ) then
        min_index = i
      end if
    end do

  end if

  return
end
subroutine r8vec_mirror_next ( n, a, done )

!*****************************************************************************80
!
!! R8VEC_MIRROR_NEXT steps through all sign variations of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    In normal use, the user would set every element of A to be positive.
!    The routine will take the input value of A, and output a copy in
!    which the signs of one or more entries have been changed.  Repeatedly
!    calling the routine with the output from the previous call will generate
!    every distinct "variation" of A; that is, all possible sign variations.
!
!    When the output variable DONE is TRUE (or equal to 1), then the
!    output value of A_NEW is the last in the series.
!
!    Note that A may have some zero values.  The routine will essentially
!    ignore such entries; more exactly, it will not stupidly assume that -0
!    is a proper "variation" of 0!
!
!    Also, it is possible to call this routine with the signs of A set
!    in any way you like.  The routine will operate properly, but it
!    will nonethess terminate when it reaches the value of A in which
!    every nonzero entry has negative sign.
!
!    More efficient algorithms using the Gray code seem to require internal
!    memory in the routine, which is not one of MATLAB's strong points,
!    or the passing back and forth of a "memory array", or the use of
!    global variables, or unnatural demands on the user.  This form of
!    the routine is about as clean as I can make it.
!
!  Example:
!
!      Input         Output
!    ---------    --------------
!    A            A_NEW     DONE
!    ---------    --------  ----
!     1  2  3     -1  2  3  false
!    -1  2  3      1 -2  3  false
!     1 -2  3     -1 -2  3  false
!    -1 -2  3      1  2 -3  false
!     1  2 -3     -1  2 -3  false
!    -1  2 -3      1 -2 -3  false
!     1 -2 -3     -1 -2 -3  false
!    -1 -2 -3      1  2  3  true
!
!     1  0  3     -1  0  3  false
!    -1  0  3      1  0 -3  false
!     1  0 -3     -1  0 -3  false
!    -1  0 -3      1  0  3  true
!
!  Modified:
!
!    19 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, real ( kind = 8 ) A(N), a vector of real numbers.
!    On output, the signs of some entries have been changed.
!
!    Output, logical DONE, is TRUE if the input vector A was the last element
!    in the series (every entry was nonpositive); the output vector is reset 
!    so that all entries are nonnegative, but presumably the ride is over!
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) positive
!
!  Seek the first strictly positive entry of A.
!
  positive = 0
  do i = 1, n
    if ( 0.0D+00 < a(i) ) then
      positive = i
      exit
    end if
  end do
!
!  If there is no strictly positive entry of A, there is no successor.
!
  if ( positive == 0 ) then
    a(1:n) = -a(1:n)
    done = .true.
    return
  end if
!
!  Otherwise, negate A up to the positive entry.
!
  a(1:positive) = -a(1:positive)
  done = .false.
  
  return
end
subroutine r8vec_nint ( n, a )

!*****************************************************************************80
!
!! R8VEC_NINT rounds entries of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    10 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, real ( kind = 8 ) A(N), the vector to be NINT'ed.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)

  a(1:n) = nint ( real ( a(1:n), kind = 8 ) )

  return
end
function r8vec_norm_l1 ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM_L1 returns the L1 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The vector L1 norm is defined as:
!
!      R8VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) ).
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose L1 norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM_L1, the L1 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) r8vec_norm_l1

  r8vec_norm_l1 = sum ( abs ( a(1:n) ) )

  return
end
function r8vec_norm_l2 ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM_L2 returns the L2 norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The vector L2 norm is defined as:
!
!      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose L2 norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM_L2, the L2 norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) r8vec_norm_l2

  r8vec_norm_l2 = sqrt ( sum ( a(1:n)**2 ) )

  return
end
function r8vec_norm_li ( n, a )

!*****************************************************************************80
!
!! R8VEC_NORM_LI returns the L-infinity norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The vector L-infinity norm is defined as:
!
!      R8VEC_NORM_LI = max ( 1 <= I <= N ) abs ( A(I) ).
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose L-infinity norm is desired.
!
!    Output, real ( kind = 8 ) R8VEC_NORM_LI, the L-infinity norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) r8vec_norm_li

  r8vec_norm_li = maxval ( abs ( a(1:n) ) )

  return
end
function r8vec_norm_lp ( n, a, p )

!*****************************************************************************80
!
!! R8VEC_NORM_LP returns the LP norm of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The vector LP norm is defined as:
!
!      R8VEC_NORM_LP = ( sum ( 1 <= I <= N ) ( abs ( A(I) ) )**P )**(1/P).
!
!    Usually, the LP norms with
!      1 <= P <= Infinity 
!    are of interest.  This routine allows
!      0 < P <= Huge ( P ).  
!    If P = Huge ( P ), then the L-infinity norm is returned, which is 
!    simply the maximum of the absolute values of the vector components.
!
!  Modified:
!
!    25 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), the vector whose LP norm is desired.
!
!    Input, real ( kind = 8 ) P, the index of the norm.  
!
!    Output, real ( kind = 8 ) R8VEC_NORM_LP, the LP norm of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) p
  real    ( kind = 8 ) r8vec_norm_lp

  if ( p <= 0.0D+00 ) then
    r8vec_norm_lp = -1.0D+00
  else if ( p == huge ( p ) ) then
    r8vec_norm_lp = maxval ( abs ( a(1:n) ) )
  else if ( p == 1.0D+00 ) then
    r8vec_norm_lp = sum ( abs ( a(1:n) ) )
  else if ( p == 2.0D+00 ) then
    r8vec_norm_lp = sqrt ( sum ( a(1:n)**2 ) )
  else
    r8vec_norm_lp = ( sum ( ( abs ( a(1:n) ) )**p ) )**( 1.0D+00 / p )
  end if

  return
end
subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
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
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    In that case, it is saved for the next call.
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is 
!    negative,then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, integer ( kind = 4 ) MADE, records the number of values that have
!    been computed.  On input with negative N, this value overwrites
!    the return value of N, so the user can get an accounting of
!    how much work has been done.
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer SAVED, is 0 or 1 depending on whether there is a
!    single saved value left over from the previous call.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
!    Local, real ( kind = 8 ) Y, the value saved from the previous call, if
!    SAVED is 1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: i4_one = 1
  integer ( kind = 4 ), parameter :: i4_two = 2
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save :: made = 0
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) r(n+1)
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ), save :: saved = 0
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
  real    ( kind = 8 ), save :: y = 0.0D+00
!
!  I'd like to allow the user to reset the internal data.
!  But this won't work properly if we have a saved value Y.
!  I'm making a crock option that allows the user to signal
!  explicitly that any internal memory should be flushed,
!  by passing in a negative value for N.
!
  if ( n < 0 ) then
    n = made
    made = 0
    saved = 0
    y = 0.0D+00
    return
  else if ( n == 0 ) then
    return
  end if
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  Use up the old value, if we have it.
!
  if ( saved == 1 ) then
    x(1) = y
    saved = 0
    x_lo_index = 2
  end if
!
!  Maybe we don't need any more values.
!
  if ( x_hi_index - x_lo_index + 1 == 0 ) then
!
!  If we need just one new value, do that here to avoid null arrays.
!
  else if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * pi * r(2) )
    y =      sqrt ( -2.0D+00 * log ( r(1) ) ) * sin ( 2.0D+00 * pi * r(2) )

    saved = 1

    made = made + 2
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + i4_one, i4_two ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m:2) )

    made = made + x_hi_index - x_lo_index + 1
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * pi * r(2*m) )

    y = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * sin ( 2.0D+00 * pi * r(2*m) )

    saved = 1

    made = made + x_hi_index - x_lo_index + 2

  end if

  return
end
subroutine r8vec_order_type ( n, a, order )

!*****************************************************************************80
!
!! R8VEC_ORDER_TYPE determines if an R8VEC is (non)strictly ascending/descending.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    15 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the array.
!
!    Input, real ( kind = 8 ) A(N), the array to be checked.
!
!    Output, integer ( kind = 4 ) ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
!
!  Search for the first value not equal to A(1).
!
  i = 1

  do

    i = i + 1

    if ( n < i ) then
      order = 0
      return
    end if

    if ( a(1) < a(i) ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( a(i) < a(1) ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  Now we have a "direction".  Examine subsequent entries.
!
  do while ( i < n )

    i = i + 1

    if ( order == 1 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 3
      end if

    end if

  end do

  return
end
subroutine r8vec_part_quick_a ( n, a, l, r )

!*****************************************************************************80
!
!! R8VEC_PART_QUICK_A reorders an R8VEC as part of a quick sort.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The routine reorders the entries of A.  Using A(1) as the key,
!    all entries of A that are less than or equal to the key will
!    precede the key which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      N = 8
!
!      A = ( 6, 7, 3, 1, 6, 8, 2, 9 )
!
!    Output:
!
!      L = 3, R = 6
!
!      A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
!            -------        -------
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input/output, real ( kind = 8 ) A(N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer ( kind = 4 ) L, R, the indices of A that define 
!    the three segments.  Let KEY = the input value of A(1).  Then
!    I <= L                 A(I) < KEY;
!         L < I < R         A(I) = KEY;
!                 R <= I    KEY < A(I).
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) key
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r
  real    ( kind = 8 ) temp

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key = a(1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( key < a(l+1) ) then
      r = r - 1
      temp = a(r)
      a(r) = a(l+1)
      a(l+1) = temp
    else if ( a(l+1) == key ) then
      m = m + 1
      temp = a(m)
      a(m) = a(l+1)
      a(l+1) = temp
      l = l + 1
    else if ( a(l+1) < key ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(i) = a(i+m)
  end do
!
!  Out of bounds here, occasionally
!
  l = l - m

  a(l+1:l+m) = key

  return
end
subroutine r8vec_permute ( n, a, p )

!*****************************************************************************80
!
!! R8VEC_PERMUTE permutes an R8VEC in place.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    This routine permutes an array of real "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!    P(I) = J means that the I-th element of the output array should be 
!    the J-th element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,   4,   5,   1,   3 )
!      A = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
!
!    Output:
!
!      A    = ( 2.0, 4.0, 5.0, 1.0, 3.0 ).
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) a_temp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = -p(istart)
      cycle

    else

      a_temp = a(istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = -p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8VEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(iput) = a_temp
          exit
        end if

        a(iput) = a(iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine r8vec_permute_uniform ( n, a, seed )

!*****************************************************************************80
!
!! R8VEC_PERMUTE_UNIFORM randomly permutes an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    01 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  call perm_uniform ( n, seed, p )

  call r8vec_permute ( n, a, p )

  return
end
subroutine r8vec_polarize ( n, a, p, a_normal, a_parallel )

!*****************************************************************************80
!
!! R8VEC_POLARIZE decomposes an R8VEC into normal and parallel components.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The (nonzero) vector P defines a direction.
!
!    The vector A can be written as the sum
!
!      A = A_normal + A_parallel
!
!    where A_parallel is a linear multiple of P, and A_normal
!    is perpendicular to P.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the vector to be polarized.
!
!    Input, real ( kind = 8 ) P(N), the polarizing direction.
!
!    Output, real ( kind = 8 ) A_NORMAL(N), A_PARALLEL(N), the normal
!    and parallel components of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) a_dot_p
  real    ( kind = 8 ) a_normal(n)
  real    ( kind = 8 ) a_parallel(n)
  real    ( kind = 8 ) p(n)
  real    ( kind = 8 ) p_norm

  p_norm = sqrt ( sum ( p(1:n)**2 ) )

  if ( p_norm == 0.0D+00 ) then
    a_normal(1:n) = a(1:n)
    a_parallel(1:n) = 0.0D+00
    return
  end if

  a_dot_p = dot_product ( a(1:n), p(1:n) ) / p_norm

  a_parallel(1:n) = a_dot_p * p(1:n) / p_norm

  a_normal(1:n) = a(1:n) - a_parallel(1:n)

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end
subroutine r8vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Modified:
!
!    16 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices
!    to print.  The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,g14.8)' ) i, a(i)
  end do

  return
end
subroutine r8vec_print2 ( n, a )

!*****************************************************************************80
!
!! R8VEC_PRINT2 prints out an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    26 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) amax
  real    ( kind = 8 ) amin
  integer ( kind = 4 ) i
  character ( len = 11 ) iform
  logical integ
  integer ( kind = 4 ) lmax
  real    ( kind = 8 ) r8_log_10
!
!  Check if all entries are integral.
!
  integ = .true.

  do i = 1, n

    if ( a(i) /= real ( int ( a(i) ), kind = 8 ) ) then
      integ = .false.
      exit
    end if

  end do
!
!  Find the range of the array.
!
  amax = maxval ( abs ( a(1:n) ) )
  amin = minval ( abs ( a(1:n) ) )
!
!  Use the information about the maximum size of an entry to
!  compute an intelligent format for use with integer entries.
!
!  Later, we might also do this for real vectors.
!
  lmax = int ( r8_log_10 ( amax ) )

  if ( integ ) then
    write ( iform, '( ''(2x,i'', i2, '')'' )' ) lmax + 3
  else
    iform = ' '
  end if

  do i = 1, n

    if ( integ ) then
      write ( *, iform ) int ( a(i) )
    else
      write ( *, '(2x,g14.6)' ) a(i)
    end if

  end do

  return
end
function r8vec_product ( n, a )

!*****************************************************************************80
!
!! R8VEC_PRODUCT returns the product of the entries of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    In FORTRAN90, this facility is offered by the built in
!    PRODUCT function:
!
!      R8VEC_PRODUCT ( N, A ) = PRODUCT ( A(1:N) )
!
!    In MATLAB, this facility is offered by the built in
!    PROD function:
!
!      R8VEC_PRODUCT ( N, A ) = PROD ( A(1:N) )
!
!  Modified:
!
!    17 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) R8VEC_PRODUCT, the product of the entries.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) r8vec_product

  r8vec_product = product ( a(1:n) )

  return
end
subroutine r8vec_range ( n, x, xmin, xmax, y, ymin, ymax )

!*****************************************************************************80
!
!! R8VEC_RANGE finds the range of Y's within a restricted X range.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The routine is given a set of pairs of points (X,Y), and a range
!    XMIN to XMAX of valid X values.  Over this range, it seeks
!    YMIN and YMAX, the minimum and maximum values of Y for
!    valid X's.
!
!  Modified:
!
!    14 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) X(N), the X array.
!
!    Input, real ( kind = 8 ) XMIN, XMAX, the range of X values to check.
!
!    Input, real ( kind = 8 ) Y(N), the Y array.
!
!    Output, real ( kind = 8 ) YMIN, YMAX, the range of Y values whose
!    X value is within the X range.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xmax
  real    ( kind = 8 ) xmin
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) ymax
  real    ( kind = 8 ) ymin

  ymin =  huge ( ymin )
  ymax = -huge ( ymax )

  do i = 1, n

    if ( xmin <= x(i) .and. x(i) <= xmax ) then

      ymin = min ( ymin, y(i) )
      ymax = max ( ymax, y(i) )

    end if

  end do

  return
end
subroutine r8vec_range_2 ( n, a, amin, amax )

!*****************************************************************************80
!
!! R8VEC_RANGE_2 updates a range to include a new array.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Given a range AMIN to AMAX, and an array A, the routine will
!    decrease AMIN if necessary, or increase AMAX if necessary, so that
!    every entry of A is between AMIN and AMAX.
!
!    However, AMIN will not be increased, nor AMAX decreased.
!
!    This routine may be used to compute the maximum and minimum of a
!    collection of arrays one at a time.
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Input/output, real ( kind = 8 ) AMIN, AMAX.  On input, the
!    current legal range of values for A.  On output, AMIN and AMAX
!    are either unchanged, or else "widened" so that all entries
!    of A are within the range.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) amax
  real    ( kind = 8 ) amin

  amax = max ( amax, maxval ( a(1:n) ) )
  amin = min ( amin, minval ( a(1:n) ) )

  return
end
subroutine r8vec_read ( input_file, n, r )

!*****************************************************************************80
!
!! R8VEC_READ reads an R8VEC from a file.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file to be read.
!
!    Input, integer ( kind = 4 ) N, the size of the R8VEC.
!
!    Output, real ( kind = 8 ) R(N), the R8VEC.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = * ) input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) k
  real    ( kind = 8 ) r(n)

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    stop
  end if

  do k = 1, n

    read ( input_unit, *, iostat = ios ) r(k)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_READ - Fatal error!'
      write ( *, '(a,i8)' ) '  I/O error while reading record ', k
      stop
    end if

  end do

  close ( unit = input_unit )

  return
end
subroutine r8vec_read_size ( input_file, n )

!*****************************************************************************80
!
!! R8VEC_READ_SIZE reads the size of an R8VEC from a file.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file to
!    be read.
!
!    Output, integer ( kind = 4 ) N, the size of the R8VEC.
!
  implicit none

  character ( len = * ) input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) n
  real    ( kind = 8 ) r

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_READ_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    stop
  end if

  n = 0

  do

    read ( input_unit, *, iostat = ios ) r

    if ( ios /= 0 ) then
      exit
    end if

    n = n + 1

  end do

  close ( unit = input_unit )

  return
end
subroutine r8vec_reverse ( n, a )

!*****************************************************************************80
!
!! R8VEC_REVERSE reverses the elements of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    In FORTRAN90, calling R8VEC_REVERSE is equivalent to
!
!      A(1:N) = A(N:1:-1)
!
!  Example:
!
!    Input:
!
!      N = 5, 
!      A = ( 11.0, 12.0, 13.0, 14.0, 15.0 ).
!
!    Output:
!
!      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 ).
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n/2
    call r8_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine r8vec_rotate ( n, a, m )

!*****************************************************************************80
!
!! R8VEC_ROTATE "rotates" the entries of an R8VEC in place.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    This routine rotates an array of real "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5, M = 2
!      A    = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
!
!    Output:
!
!      A    = ( 4.0, 5.0, 1.0, 2.0, 3.0 ).
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
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) M, the number of positions to the right that
!    each element should be moved.  Elements that shift pass position
!    N "wrap around" to the beginning of the array.
!
!    Input/output, real ( kind = 8 ) A(N), the array to be rotated.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mcopy
  integer ( kind = 4 ) nset
  real    ( kind = 8 ) temp
!
!  Force M to be positive, between 0 and N-1.
!
  mcopy = i4_modp ( m, n )

  if ( mcopy == 0 ) then
    return
  end if

  istart = 0
  nset = 0

  do

    istart = istart + 1

    if ( n < istart ) then
      exit
    end if

    temp = a(istart)
    iget = istart
!
!  Copy the new value into the vacated entry.
!
    do

      iput = iget

      iget = iget - mcopy
      if ( iget < 1 ) then
        iget = iget + n
      end if

      if ( iget == istart ) then
        exit
      end if

      a(iput) = a(iget)
      nset = nset + 1

    end do

    a(iput) = temp
    nset = nset + 1

    if ( n <= nset ) then
      exit
    end if

  end do

  return
end
subroutine r8vec_search_binary_a ( n, a, aval, indx )

!*****************************************************************************80
!
!! R8VEC_SEARCH_BINARY_A searches an ascending sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Binary search is used.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.9,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 26.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the array.
!
!    Input, real ( kind = 8 ) A(N), the array to be searched.  The array must
!    be sorted in ascending order.
!
!    Input, real ( kind = 8 ) AVAL, the value to be searched for.
!
!    Output, integer ( kind = 4 ) INDX, the result of the search.
!    -1, AVAL does not occur in the array.
!    I, A(I) = AVAL.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) aval
  integer ( kind = 4 ) high
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid

  indx = -1

  low = 1
  high = n

  do while ( low <= high )

    mid = ( low + high ) / 2

    if ( a(mid) == aval ) then
      indx = mid
      exit
    else if ( a(mid) < aval ) then
      low = mid + 1
    else if ( aval < a(mid) ) then
      high = mid - 1
    end if

  end do

  return
end
subroutine r8vec_sort_bubble_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_BUBBLE_A ascending sorts an R8VEC using bubble sort.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n - 1 
    do j = i + 1, n
      if ( a(j) < a(i) ) then
        call r8_swap ( a(i), a(j) )
      end if
    end do
  end do

  return
end
subroutine r8vec_sort_bubble_d ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_BUBBLE_D descending sorts an R8VEC using bubble sort.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!
!  Modified:
!
!    31 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n - 1
    do j = i + 1, n
      if ( a(i) < a(j) ) then
        call r8_swap ( a(i), a(j) )
      end if
    end do
  end do

  return
end
subroutine r8vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_A ascending sorts an R8VEC using heap sort.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    07 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) n1
  real    ( kind = 8 ) temp

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call r8vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  temp = a(1)
  a(1) = a(n)
  a(n) = temp
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call r8vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    temp = a(1)
    a(1) = a(n1)
    a(n1) = temp

  end do

  return
end
subroutine r8vec_sort_heap_d ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_D descending sorts an R8VEC using heap sort.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    07 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) n1

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into ascending heap form.
!
  call r8vec_heap_a ( n, a )
!
!  2: Sort A.
!
!  The smallest object in the heap is in A(1).
!  Move it to position A(N).
!
  call r8_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call r8vec_heap_a ( n1, a )
!
!  Take the smallest object from A(1) and move it to A(N1).
!
    call r8_swap ( a(1), a(n1) )

  end do

  return
end
subroutine r8vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I:N)) is sorted,
!
!    or explicitly, by the call
!
!      call r8vec_permute ( n, a, indx )
!
!    after which A(1:N) is sorted.
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8vec_sort_heap_index_d ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(1:N)) is sorted,
!
!    or explicitly, by the call
!
!      call r8vec_permute ( n, a, index )
!
!    after which A(1:N) is sorted.
!
!  Modified:
!
!    21 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j+1)) < a(indx(j)) ) then
          j = j + 1
        end if
      end if

      if ( a(indx(j)) < aval ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8vec_sort_heap_mask_a ( n, a, mask_num, mask, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_MASK_A: indexed heap ascending sort of a masked R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    An array A is given.  An array MASK of indices into A is given.
!    The routine produces a vector INDX, which is a permutation of the
!    entries of MASK, so that:
!
!      A(MASK(INDX(I)) <= A(MASK(INDX(J))
!
!    whenever
!
!      I <= J
!
!    In other words, only the elements of A that are indexed by MASK
!    are to be considered, and the only thing that happens is that
!    a rearrangment of the indices in MASK is returned that orders the
!    masked elements.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Input, integer ( kind = 4 ) MASK_NUM, the number of mask elements.
!
!    Input, integer ( kind = 4 ) MASK(MASK_NUM), the mask array.  This is
!    simply a list of indices of A.  The entries of MASK should
!    be unique, and each one should be between 1 and N.
!
!    Output, integer ( kind = 4 ) INDX(MASK_NUM), the sort index.  There are 
!    MASK_NUM elements of A selected by MASK.  If we want to list those
!    elements in order, then the I-th element is A(MASK(INDX(I))).
!
  implicit none

  integer ( kind = 4 ) mask_num
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(mask_num)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mask(mask_num)

  if ( n < 1 ) then
    return
  end if

  if ( mask_num < 1 ) then
    return
  end if

  if ( mask_num == 1 ) then
    indx(1) = 1
    return
  end if

  call i4vec_indicator ( mask_num, indx )

  l = mask_num / 2 + 1
  ir = mask_num

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(mask(indxt))

    else

      indxt = indx(ir)
      aval = a(mask(indxt))
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(mask(indx(j))) < a(mask(indx(j+1))) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(mask(indx(j))) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8vec_sort_insert_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_INSERT_A ascending sorts an R8VEC using an insertion sort.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.1,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( a(j) <= x ) then
        exit
      end if

      a(j+1) = a(j)
      j = j - 1

    end do

    a(j+1) = x

  end do

  return
end
subroutine r8vec_sort_insert_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_INSERT_INDEX_A ascending index sorts an R8VEC using insertion.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.1,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(N), the array to be sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sorted indices.  The array 
!    is sorted when listed from A(INDX(1)) through A(INDX(N)).
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x

  if ( n < 1 ) then
    return
  end if

  call i4vec_indicator ( n, indx )

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( a(indx(j)) <= x ) then
        exit
      end if

      indx(j+1) = indx(j)
      j = j - 1

    end do

    indx(j+1) = i

  end do

  return
end
subroutine r8vec_sort_insert_index_d ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_INSERT_INDEX_D descending index sorts an R8VEC using insertion.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    07 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Algorithm 1.1,
!    Combinatorial Algorithms,
!    CRC Press, 1998, page 11.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the vector.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(N), the array to be sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sorted indices.  The array 
!    is sorted when listed from A(INDX(1)) through A(INDX(N)).
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x

  if ( n < 1 ) then
    return
  end if

  call i4vec_indicator ( n, indx )

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( x <= a(indx(j)) ) then
        exit
      end if

      indx(j+1) = indx(j)
      j = j - 1

    end do

    indx(j+1) = i

  end do

  return
end
subroutine r8vec_sort_quick_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_QUICK_A ascending sorts an R8VEC using quick sort.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Example:
!
!    Input:
!
!      N = 7
!      A = ( 6, 7, 3, 2, 9, 1, 8 )
!
!    Output:
!
!      A = ( 1, 2, 3, 6, 7, 8, 9 )
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: level_max = 25
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) l_segment
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n_segment
  integer ( kind = 4 ) rsave(level_max)
  integer ( kind = 4 ) r_segment

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_SORT_QUICK_A - Fatal error!'
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
    call r8vec_part_quick_a ( n_segment, a(base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( level_max < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i8)' ) '  Exceeding recursion maximum of ', level_max
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

        if ( 0 < n_segment ) then
          exit
        end if

      end do

    end if

  end do

  return
end
subroutine r8vec_sort_shell_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_SHELL_A ascending sorts an R8VEC using Shell's sort.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an array to be sorted.
!    On output, the sorted array.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) asave
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) ipow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) maxpow

  if ( n <= 1 ) then
    return
  end if
!
!  Determine the smallest MAXPOW so that
!    N <= ( 3**MAXPOW - 1 ) / 2
!
  maxpow = 1

  do while ( 3**maxpow < 2 * n + 1 )
    maxpow = maxpow + 1
  end do

  if ( 1 < maxpow ) then
    maxpow = maxpow - 1
  end if
!
!  Now sort groups of size ( 3**IPOW - 1 ) / 2.
!
  do ipow = maxpow, 1, -1

    inc = ( 3**ipow - 1 ) / 2
!
!  Sort the values with indices equal to K mod INC.
!
    do k = 1, inc
!
!  Insertion sort of the items with index
!  INC+K, 2*INC+K, 3*INC+K, ...
!
      do i = inc+k, n, inc

        asave = a(i)
        ifree = i
        j = i - inc

        do

          if ( j < 1 ) then
            exit
          end if

          if ( a(j) <= asave ) then
            exit
          end if

          ifree = j
          a(j+inc) = a(j)
          j = j - inc

        end do

        a(ifree) = asave

      end do

    end do

  end do

  return
end
subroutine r8vec_sort2_a ( n, x, y )

!*****************************************************************************80
!
!! R8VEC_SORT2_A ascending sorts an R8VEC and adjusts an associated R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The routine sorts the elements of X, and whenever
!    an element of X is moved, the corresponding element of
!    Y is moved in the same way.  This action means that after
!    the sorting, every element of X is still paired to the
!    same Y value.
!
!    If you have more than one array associated with X, or
!    an integer array, or some other complication, you may want to
!    look at doing an "indexed sort" instead.
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
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, an unsorted array.
!    On output, X has been sorted.
!
!    Input/output, real ( kind = 8 ) Y(N), an array which is to be
!    shifted corresponding to the shifts made in X.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) y(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      call r8_swap ( x(i), x(j) )
      call r8_swap ( y(i), y(j) )

    else if ( indx < 0 ) then

      if ( x(i) <= x(j) ) then
        isgn = -1
      else
        isgn = + 1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8vec_sorted_merge_a ( na, a, nb, b, nc, c )

!*****************************************************************************80
!
!! R8VEC_SORTED_MERGE_A merges two ascending sorted R8VEC's.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The elements of A and B should be sorted in ascending order.
!
!    The elements in the output array C will also be in ascending order,
!    and unique.
!
!    The output vector C may share storage with A or B.
!
!  Modified:
!
!    23 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(NA), the first sorted array.
!
!    Input, integer ( kind = 4 ) NB, the dimension of B.
!
!    Input, real ( kind = 8 ) B(NB), the second sorted array.
!
!    Output, integer ( kind = 4 ) NC, the number of elements in the output
!    array.  Note that C should usually be dimensioned at least NA+NB in the
!    calling routine.
!
!    Output, real ( kind = 8 ) C(NC), the merged unique sorted array.
!
  implicit none

  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb

  real    ( kind = 8 ) a(na)
  real    ( kind = 8 ) b(nb)
  real    ( kind = 8 ) c(na+nb)
  real    ( kind = 8 ) d(na+nb)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) na2
  integer ( kind = 4 ) nb2
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) order

  na2 = na
  nb2 = nb

  ja = 0
  jb = 0
  nc = 0

  call r8vec_order_type ( na2, a, order )

  if ( order < 0 .or. 2 < order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_SORTED_MERGE_A - Fatal error!'
    write ( *, '(a)' ) '  The input array A is not ascending sorted!'
    stop
  end if

  call r8vec_order_type ( nb2, b, order )

  if ( order < 0 .or. 2 < order ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_SORTED_MERGE_A - Fatal error!'
    write ( *, '(a)' ) '  The input array B is not ascending sorted!'
    stop
  end if

  do
!
!  If we've used up all the entries of A, stick the rest of B on the end.
!
    if ( na2 <= ja ) then

      do j = 1, nb2 - jb
        jb = jb + 1
        if ( nc == 0 ) then
          nc = nc + 1
          d(nc) = b(jb)
        else if ( d(nc) < b(jb) ) then
          nc = nc + 1
          d(nc) = b(jb)
        end if
      end do

      c(1:nc) = d(1:nc)

      exit
!
!  If we've used up all the entries of B, stick the rest of A on the end.
!
    else if ( nb2 <= jb ) then

      do j = 1, na2 - ja
        ja = ja + 1
        if ( nc == 0 ) then
          nc = nc + 1
          d(nc) = a(ja)
        else if ( d(nc) < a(ja) ) then
          nc = nc + 1
          d(nc) = a(ja)
        end if
      end do

      c(1:nc) = d(1:nc)

      exit
!
!  Otherwise, if the next entry of A is smaller, that's our candidate.
!
    else if ( a(ja+1) <= b(jb+1) ) then

      ja = ja + 1
      if ( nc == 0 ) then
        nc = nc + 1
        d(nc) = a(ja)
      else if ( d(nc) < a(ja) ) then
        nc = nc + 1
        d(nc) = a(ja)
      end if
!
!  ...or if the next entry of B is the smaller, consider that.
!
    else

      jb = jb + 1
      if ( nc == 0 ) then
        nc = nc + 1
        d(nc) = b(jb)
      else if ( d(nc) < b(jb) ) then
        nc = nc + 1
        d(nc) = b(jb)
      end if
    end if

  end do

  return
end
function r8vec_sorted_nearest ( n, a, value )

!*****************************************************************************80
!
!! R8VEC_SORTED_NEAREST returns the nearest element in a sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    29 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), a sorted vector.
!
!    Input, real ( kind = 8 ) VALUE, the value whose nearest vector
!    entry is sought.
!
!    Output, integer ( kind = 4 ) R8VEC_SORTED_NEAREST, the index of the nearest
!    entry in the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) r8vec_sorted_nearest
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  real    ( kind = 8 ) value

  if ( n < 1 ) then
    r8vec_sorted_nearest = -1
    return
  end if

  if ( n == 1 ) then
    r8vec_sorted_nearest = 1
    return
  end if

  if ( a(1) < a(n) ) then

    if ( value < a(1) ) then
      r8vec_sorted_nearest = 1
      return
    else if ( a(n) < value ) then
      r8vec_sorted_nearest = n
      return
    end if
!
!  Seek an interval containing the value.
!
    lo = 1
    hi = n

    do while ( lo < hi - 1 )

      mid = ( lo + hi ) / 2

      if ( value == a(mid) ) then
        r8vec_sorted_nearest = mid
        return
      else if ( value < a(mid) ) then
        hi = mid
      else
        lo = mid
      end if

    end do
!
!  Take the nearest.
!
    if ( abs ( value - a(lo) ) < abs ( value - a(hi) ) ) then
      r8vec_sorted_nearest = lo
    else
      r8vec_sorted_nearest = hi
    end if

    return
!
!  A descending sorted vector A.
!
  else

    if ( value < a(n) ) then
      r8vec_sorted_nearest = n
      return
    else if ( a(1) < value ) then
      r8vec_sorted_nearest = 1
      return
    end if
!
!  Seek an interval containing the value.
!
    lo = n
    hi = 1

    do while ( lo < hi - 1 )

      mid = ( lo + hi ) / 2

      if ( value == a(mid) ) then
        r8vec_sorted_nearest = mid
        return
      else if ( value < a(mid) ) then
        hi = mid
      else
        lo = mid
      end if

    end do
!
!  Take the nearest.
!
    if ( abs ( value - a(lo) ) < abs ( value - a(hi) ) ) then
      r8vec_sorted_nearest = lo
    else
      r8vec_sorted_nearest = hi
    end if

    return

  end if

  return
end
subroutine r8vec_sorted_split ( n, a, split, i_lt, i_gt )

!*****************************************************************************80
!
!! R8VEC_SORTED_SPLIT "splits" a sorted R8VEC, given a splitting value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Given a splitting value SPLIT, the routine seeks indices 
!    I_LT and I_GT so that
!
!      A(I_LT) < SPLIT < A(I_GT),
!
!    and if there are intermediate index values between I_LT and
!    I_GT, then those entries of A are exactly equal to SPLIT.
!
!  Modified:
!
!    18 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, real ( kind = 8 ) A(N), a sorted array.
!
!    Input, real ( kind = 8 ) SPLIT, a value to which the entries in A are
!    to be compared.
!
!    Output, integer ( kind = 4 ) I_LT:
!    0 if no entries are less than SPLIT;
!    N if all entries are less than SPLIT;
!    otherwise, the index of the last entry in A less than SPLIT.
!
!    Output, integer ( kind = 4 ) I_GT:
!    1 if all entries are greater than SPLIT;
!    N+1 if no entries are greater than SPLIT;
!    otherwise the index of the first entry in A greater than SPLIT.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_gt
  integer ( kind = 4 ) i_lt
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  real    ( kind = 8 ) split

  if ( n < 1 ) then
    i_lt = -1
    i_gt = -1
    return
  end if

  if ( split < a(1) ) then
    i_lt = 0
    i_gt = 1
    return
  end if

  if ( a(n) < split ) then
    i_lt = n
    i_gt = n + 1
    return
  end if

  lo = 1
  hi = n

  do

    if ( lo + 1 == hi ) then
      i_lt = lo
      exit
    end if

    mid = ( lo + hi ) / 2

    if ( split <= a(mid) ) then
      hi = mid
    else
      lo = mid
    end if

  end do

  do i = i_lt + 1, n
    if ( split < a(i) ) then
      i_gt = i
      return
    end if
  end do

  i_gt = n + 1

  return
end
subroutine r8vec_sorted_unique ( n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8VEC_SORTED_UNIQUE keeps the unique elements in a sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    29 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, the sorted array of N elements;
!    On output, the sorted unique array of UNIQUE_NUM elements.
!
!    Input, real ( kind = 8 ) TOL, a nonnegative tolerance for equality.
!    Set it to 0.0 for the strictest test.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements 
!    of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) unique_num
  real    ( kind = 8 ) tol

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do i = 2, n

    if ( tol < abs ( a(i) - a(unique_num) ) ) then
      unique_num = unique_num + 1
      a(unique_num) = a(i)
    end if

  end do

  return
end
subroutine r8vec_sorted_unique_count ( n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8VEC_SORTED_UNIQUE_COUNT counts the unique elements in a sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Because the array is sorted, this algorithm is O(N).
!
!  Modified:
!
!    29 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), the sorted array to examine.
!
!    Input, real ( kind = 8 ) TOL, a nonnegative tolerance for equality.  
!    Set it to 0.0 for the strictest test.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements 
!    of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) unique_num
  real    ( kind = 8 ) tol

  if ( n < 1 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do i = 2, n

    if ( tol < abs ( a(i-1) - a(i) ) ) then
      unique_num = unique_num + 1
    end if

  end do

  return
end
subroutine r8vec_sorted_unique_hist ( n, a, tol, maxuniq, unique_num, &
  auniq, acount )

!*****************************************************************************80
!
!! R8VEC_SORTED_UNIQUE_HIST histograms the unique elements of a sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    09 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), the array to examine.  The elements of A
!    should have been sorted.
!
!    Input, real ( kind = 8 ) TOL, a nonnegative tolerance for equality.
!    Set it to 0.0 for the strictest test.
!
!    Input, integer ( kind = 4 ) MAXUNIQ, the maximum number of unique elements
!    that can be handled.  If there are more than MAXUNIQ unique
!    elements in A, the excess will be ignored.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements
!    of A.
!
!    Output, real ( kind = 8 ) AUNIQ(UNIQUE_NUM), the unique elements of A.
!
!    Output, integer ( kind = 4 ) ACOUNT(UNIQUE_NUM), the number of times
!    each element of AUNIQ occurs in A.
!
  implicit none

  integer ( kind = 4 ) maxuniq
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) acount(maxuniq)
  real    ( kind = 8 ) auniq(maxuniq)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) unique_num
  real    ( kind = 8 ) tol
!
!  Start taking statistics.
!
  unique_num = 0

  do i = 1, n

    if ( i == 1 ) then

      unique_num = 1
      auniq(unique_num) = a(1)
      acount(unique_num) = 1

    else if ( abs ( a(i) - auniq(unique_num) ) <= tol ) then

      acount(unique_num) = acount(unique_num) + 1

    else if ( unique_num < maxuniq ) then

      unique_num = unique_num + 1
      auniq(unique_num) = a(i)
      acount(unique_num) = 1

    end if

  end do

  return
end
subroutine r8vec_split ( n, a, split, isplit )

!*****************************************************************************80
!
!! R8VEC_SPLIT "splits" an unsorted R8VEC based on a splitting value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the vector is already sorted, it is simpler to do a binary search
!    on the data than to call this routine.
!
!    The vector is not assumed to be sorted before input, and is not
!    sorted during processing.  If sorting is not needed, then it is
!    more efficient to use this routine.
!
!  Modified:
!
!    21 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, real ( kind = 8 ) A(N), the array to split.  On output,
!    all the entries of A that are less than or equal to SPLIT
!    are in A(1:ISPLIT).
!
!    Input, real ( kind = 8 ) SPLIT, the value used to split the vector.
!    It is not necessary that any value of A actually equal SPLIT.
!
!    Output, integer ( kind = 4 ) ISPLIT, indicates the position of the last
!    entry of the split vector that is less than or equal to SPLIT.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) isplit
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  real    ( kind = 8 ) split
!
!  Partition the vector into A1, A2, A3, where
!    A1 = A(I1:J1) holds values <= SPLIT,
!    A2 = A(I2:J2) holds untested values,
!    A3 = A(I3:J3) holds values > SPLIT.
!
  i1 = 1
  j1 = 0

  i2 = 1
  j2 = n

  i3 = n + 1
  j3 = n
!
!  Pick the next item from A2, and move it into A1 or A3.
!  Adjust indices appropriately.
!
  do i = 1, n

    if ( a(i2) <= split ) then
      i2 = i2 + 1
      j1 = j1 + 1
    else
      call r8_swap ( a(i2), a(i3-1) )
      i3 = i3 - 1
      j2 = j2 - 1
    end if

  end do

  isplit = j1

  return
end
subroutine r8vec_std ( n, a, std )

!*****************************************************************************80
!
!! R8VEC_STD returns the standard deviation of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The standard deviation of a vector X of length N is defined as
!
!      mean ( X(1:n) ) = sum ( X(1:n) ) / n
!
!      std ( X(1:n) ) = sqrt ( sum ( ( X(1:n) - mean )**2 ) / ( n - 1 ) )
!
!  Modified:
!
!    06 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!    N should be at least 2.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, real ( kind = 8 ) STD, the standard deviation of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) mean
  real    ( kind = 8 ) std

  if ( n < 2 ) then

    std = 0.0D+00

  else

    mean = sum ( a(1:n) ) / real ( n, kind = 8 )

    std = sum ( ( a(1:n) - mean )**2 )

    std = sqrt ( std / real ( n - 1, kind = 8 ) )

  end if

  return
end
function r8vec_sum ( n, a )

!*****************************************************************************80
!
!! R8VEC_SUM returns the sum of the entries of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    In FORTRAN90, this facility is offered by the built in
!    SUM function:
!
!      R8VEC_SUM ( N, A ) = SUM ( A(1:N) )
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) R8VEC_SUM, the sum of the entries.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) r8vec_sum

  r8vec_sum = sum ( a(1:n) )

  return
end
subroutine r8vec_swap ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_SWAP swaps the entries of two R8VECs.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    04 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the arrays.
!
!    Input/output, real ( kind = 8 ) A1(N), A2(N), the vectors to swap.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)
  real    ( kind = 8 ) a3(n)

  a3(1:n) = a1(1:n)
  a1(1:n) = a2(1:n)
  a2(1:n) = a3(1:n)

  return
end
subroutine r8vec_uniform ( n, a, b, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM returns a scaled pseudorandom R8VEC
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    29 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Philip Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) r(n)

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge ( )
    end if

    r(i) = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec_unique_count ( n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8VEC_UNIQUE_COUNT counts the unique elements in an unsorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Because the array is unsorted, this algorithm is O(N^2).
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), the unsorted array to examine.
!
!    Input, real ( kind = 8 ) TOL, a nonnegative tolerance for equality. 
!    Set it to 0.0 for the strictest test.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements
!    of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) unique_num
  real    ( kind = 8 ) tol

  unique_num = 0

  do i = 1, n

    unique_num = unique_num + 1

    do j = 1, i - 1

      if ( abs ( a(i) - a(j) ) <= tol ) then
        unique_num = unique_num - 1
        exit
      end if

    end do

  end do

  return
end
subroutine r8vec_unit_euclidean ( n, a )

!*****************************************************************************80
!
!! R8VEC_UNIT_EUCLIDEAN normalizes an R8VEC in the Euclidean norm.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The euclidean norm is also sometimes called the l2 or
!    least squares norm.
!
!  Modified:
!
!    20 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input/output, real ( kind = 8 ) A(N), the vector to be normalized.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) norm

  norm = sqrt ( sum ( a(1:n)**2 ) )

  if ( norm /= 0.0D+00 ) then
    a(1:n) = a(1:n) / norm
  end if

  return
end
subroutine r8vec_unit_sum ( n, a )

!*****************************************************************************80
!
!! R8VEC_UNIT_SUM normalizes an R8VEC to have unit sum.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    10 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, real ( kind = 8 ) A(N), the vector to be normalized.
!    On output, the entries of A should have unit sum.  However, if 
!    the input vector has zero sum, the routine halts.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) a_sum

  a_sum = sum ( a(1:n) )

  if ( a_sum == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIT_SUM - Fatal error!'
    write ( *, '(a)' ) '  The vector entries sum to 0.'
    stop
  end if

  a(1:n) = a(1:n) / a_sum

  return
end
subroutine r8vec_variance ( n, a, variance )

!*****************************************************************************80
!
!! R8VEC_VARIANCE returns the variance of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The variance of a vector X of length N is defined as
!
!      mean ( X(1:n) ) = sum ( X(1:n) ) / n
!
!      var ( X(1:n) ) = sum ( ( X(1:n) - mean )**2 ) / ( n - 1 )
!
!  Modified:
!
!    14 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!    N should be at least 2.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) mean
  real    ( kind = 8 ) variance

  if ( n < 2 ) then

    variance = 0.0D+00

  else

    mean = sum ( a(1:n) ) / real ( n, kind = 8 )

    variance = sum ( ( a(1:n) - mean )**2 )

    variance = variance / real ( n - 1, kind = 8 )

  end if

  return
end
subroutine r8vec_write ( n, r, output_file )

!*****************************************************************************80
!
!! R8VEC_WRITE writes an R8VEC to a file.
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) R(N), the vector to be written.
!
!    Input, character ( len = * ) OUTPUT_FILE, the name of the file to which
!    the information is to be written.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  character ( len = * ) output_file
  integer ( kind = 4 ) output_unit
  real    ( kind = 8 ) r(n)

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace' )

  do i = 1, n
    write ( output_unit, '(2x,g16.8)' ) r(i)
  end do

  close ( unit = output_unit )

  return
end
subroutine r8vec_zero ( n, a )

!*****************************************************************************80
!
!! R8VEC_ZERO zeroes out an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    10 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Output, real ( kind = 8 ) A(N), the vector to be zeroed.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)

  a(1:n) = 0.0D+00

  return
end
subroutine r8vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! R8VEC2_COMPARE compares two entries in an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of real values, stored
!    as two separate vectors A1 and A2.
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    A1(I) A2(I)   A1(J) A2(J)  ISGN
!    -----------   -----------  ----
!    1.0   5.0  <  1.0   6.0     -1
!    1.0   5.0  <  2.0   8.0     -1
!    1.0   5.0  <  9.0   1.0     -1
!    1.0   5.0  =  1.0   5.0      0
!    1.0   5.0  >  0.0   2.0     +1
!    1.0   5.0  >  0.0   5.0     +1
!    1.0   5.0  >  1.0   3.0     +1
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the two components of each item.
!
!    Input, integer ( kind = 4 ) I, J, the items to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item I > item J.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(j) < a2(i) ) then
      isgn = +1
    end if

  else if ( a1(j) < a1(i) ) then

    isgn = +1

  end if

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of real values, stored
!    as two separate vectors A1 and A2.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  if ( all ( a1(1:n) == aint ( a1(1:n) ) ) .and. &
       all ( a2(1:n) == aint ( a2(1:n) ) ) ) then
    do i = 1, n
      write ( *, '(i8,2i8)' ) i, int ( a1(i) ), int ( a2(i) )
    end do
  else if ( all ( abs ( a1(1:n) ) < 1000000.0D+00 ) .and. &
            all ( abs ( a2(1:n) ) < 1000000.0D+00 ) ) then
    do i = 1, n
      write ( *, '(i8,2f14.6)' ) i, a1(i), a2(i)
    end do
  else
    do i = 1, n
      write ( *, '(i8,2g14.6)' ) i, a1(i), a2(i)
    end do
  end if

  return
end
subroutine r8vec2_print_some ( n, x1, x2, max_print, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT_SOME prints "some" of an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of real values, stored
!    as two separate vectors A1 and A2.
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vectors, is no more than MAX_PRINT, then
!    the entire vectors are printed, one entry of each per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vectors.
!
!    Input, real ( kind = 8 ) X1(N), X2(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines 
!    to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title
  real    ( kind = 8 ) x1(n)
  real    ( kind = 8 ) x2(n)

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '
  end if

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print-2
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    write ( *, '(a)' ) '  ......  ..............  ..............'
    i = n
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
    end do
    i = max_print
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,a)' ) i, x1(i), x2(i), &
      '...more entries...'

  end if

  return
end
subroutine r8vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC2_SORT_A ascending sorts an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of real values, stored
!    as two separate vectors A1 and A2.
!
!    Each item to be sorted is a pair (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, real ( kind = 8 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( n <= 1 ) then
    return
  end if
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
    if ( 0 < indx ) then

      call r8_swap ( a1(i), a1(j) )
      call r8_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call r8vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8vec2_sort_d ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC2_SORT_D descending sorts an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of real values, stored
!    as two separate vectors A1 and A2.
!
!    Each item to be sorted is a pair (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, real ( kind = 8 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( n <= 1 ) then
    return
  end if
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
    if ( 0 < indx ) then

      call r8_swap ( a1(i), a1(j) )
      call r8_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!  Reverse the value of ISGN to effect a descending sort.
!
    else if ( indx < 0 ) then

      call r8vec2_compare ( n, a1, a2, i, j, isgn )

      isgn = -isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine r8vec2_sort_heap_index_a ( n, x, y, indx )

!*****************************************************************************80
!
!! R8VEC2_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of real values, stored
!    as two separate vectors A1 and A2.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    ( X(I), Y(I) ) < ( X(J), Y(J) ) if:
!
!    * X(I) < X(J), or
!
!    * X(I) = X(J), and Y(I) < Y(J).
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      ( X(INDX(1:N)), Y(INDX(1:N) ), is sorted,
!
!    or explicitly, by the call
!
!      call r8vec_permute ( n, x, indx )
!      call r8vec_permute ( n, y, indx )
!
!    after which ( X(1:N), Y(1:N) ), is sorted.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) X(N),Y(N), pairs of X, Y coordinates of points.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array has coordinates ( X(INDX(I)), Y(INDX(I) ).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xval
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) yval

  if ( n < 1 ) then
    return
  end if

  call i4vec_indicator ( n, indx )

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      xval = x(indxt)
      yval = y(indxt)

    else

      indxt = indx(ir)
      xval = x(indxt)
      yval = y(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then

        if ( x(indx(j)) < x(indx(j+1)) .or. &
          ( x(indx(j)) == x(indx(j+1)) .and. y(indx(j)) < y(indx(j+1)) ) ) then
          j = j + 1
        end if

      end if

      if ( xval < x(indx(j)) .or. &
          ( xval == x(indx(j)) .and. yval < y(indx(j)) ) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8vec2_sorted_unique ( n, a1, a2, unique_num )

!*****************************************************************************80
!
!! R8VEC2_SORTED_UNIQUE keeps unique elements in a sorted R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of real values, stored
!    as two separate vectors A1 and A2.
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
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items.
!
!    Input/output, real ( kind = 8 ) A1(N), A2(N).
!    On input, the array of N items.
!    On output, an array of unique items.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique items.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1

  do itest = 2, n

    if ( a1(itest) /= a1(unique_num) .or. a2(itest) /= a2(unique_num) ) then

      unique_num = unique_num + 1

      a1(unique_num) = a1(itest)
      a2(unique_num) = a2(itest)

    end if

  end do

  return
end
subroutine r8vec2_sorted_unique_index ( n, a1, a2, unique_num, indx )

!*****************************************************************************80
!
!! R8VEC2_SORTED_UNIQUE_INDEX indexes unique elements in a sorted R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of real values, stored
!    as two separate vectors A1 and A2.
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it should be the
!    case that equal items are stored in adjacent vector locations.
!
!  Modified:
!
!    02 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the array of N items.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique items.
!
!    Output, integer ( kind = 4 ) INDX(N), contains in entries 1 through 
!    UNIQUE_NUM an index array of the unique items.  To build new arrays 
!    with no repeated elements:
!      B1(1:UNIQUE_NUM) = A1(INDX(1:UNIQUE_NUM))
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  indx(1) = 1

  do itest = 2, n

    if ( a1(itest-1) /= a1(itest) .or. a2(itest-1) /= a2(itest) ) then

      unique_num = unique_num + 1

      indx(unique_num) = itest

    end if

  end do

  indx(unique_num+1:n) = 0

  return
end
subroutine r8vec2_sum_max_index ( n, a, b, sum_max_index )

!*****************************************************************************80
!
!! R8VEC2_SUM_MAX_INDEX returns the index of the maximum sum of two R8VEC's.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of real values, stored
!    as two separate vectors A1 and A2.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), B(N), two arrays whose sum
!    is to be examined.
!
!    Output, integer ( kind = 4 ) SUM_MAX_INDEX, the index of the largest 
!    entry in A+B.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) sum_max
  integer ( kind = 4 ) sum_max_index

  if ( n <= 0 ) then

    sum_max_index = -1

  else

    sum_max_index = 1
    sum_max = a(1) + b(1)

    do i = 2, n
      if ( sum_max < a(i) + b(i) ) then
        sum_max = a(i) + b(i)
        sum_max_index = i
      end if
    end do

  end if

  return
end
subroutine r8vec3_print ( n, a1, a2, a3, title )

!*****************************************************************************80
!
!! R8VEC3_PRINT prints an R8VEC3.
!
!  Modified:
!
!    12 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), A3(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)
  real    ( kind = 8 ) a3(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  if ( all ( a1(1:n) == aint ( a1(1:n) ) ) .and. &
       all ( a2(1:n) == aint ( a2(1:n) ) ) .and. &
       all ( a3(1:n) == aint ( a3(1:n) ) )) then
    do i = 1, n
      write ( *, '(i8,3i8)' ) i, int ( a1(i) ), int ( a2(i) ), int ( a3(i) )
    end do
  else if ( all ( abs ( a1(1:n) ) < 1000000.0D+00 )  .and. &
            all ( abs ( a2(1:n) ) < 1000000.0D+00 )  .and. &
            all ( abs ( a3(1:n) ) < 1000000.0D+00 ) ) then
    do i = 1, n
      write ( *, '(i8,3f14.6)' ) i, a1(i), a2(i), a3(i)
    end do
  else
    do i = 1, n
      write ( *, '(i8,3g14.6)' ) i, a1(i), a2(i), a3(i)
    end do
  end if

  return
end
function radians_to_degrees ( radians )

!*****************************************************************************80
!
!! RADIANS_TO_DEGREES converts an angle measure from radians to degrees.
!
!  Modified:
!
!    11 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RADIANS, the angle measure in radians.
!
!    Output, real ( kind = 8 ) RADIANS_TO_DEGREES, the angle measure in degrees.
!
  implicit none

  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) radians
  real    ( kind = 8 ) radians_to_degrees

  radians_to_degrees = ( radians / pi ) * 180.0D+00

  return
end
subroutine random_initialize ( seed_input )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!
!    The G95 compiler is currently generating internal compiler errors
!    when it tries to compile this routine.  Just one more exasperation
!    on the mountain of complications because of the ragged interface with
!    the nonstandard random number generator standard!
!
!    If you don't initialize the FORTRAN90 random number generator
!    routine RANDOM_NUMBER, which is used by calls like
!
!      call random_number ( harvest = x )
!
!    then its behavior is not specified.  That may be OK for you.  But
!    you may want to be able to force the same sequence to be generated
!    each time, or to force a different sequence.  
!
!    To control the sequence of random numbers, you need to set the seed.
!    In FORTRAN90, this is done by calling the RANDOM+SEED routine.
!    You can call it with no arguments, in fact.  But if you call
!    it with no arguments:
!
!      call random_seed ( )
!
!    then its behavior (or more particularly, the behavior of RANDOM_NUMBER)
!    is still not specified.  You might hope that the system will go out
!    and pick a nice random seed for you, but there's no guarantee. 
!
!
!    For example, on the DEC ALPHA, if you compile a program that calls
!    RANDOM_NUMBER, then every time you run it, you get the same sequence
!    of "random" values.  If you compile a program that calls RANDOM_SEED
!    with no arguments, and then calls RANDOM_NUMBER, you still get the
!    same sequence each time you run the program.
!
!    In order to actually try to scramble up the random number generator 
!    a bit, this routine goes through the tedious process of getting the 
!    size of the random number seed, making up values based on the current 
!    time, and setting the random number seed.
!
!    Unfortunately, the RANDOM_SEED routine has a very elastic definition.
!    It does not use a single scalar integer SEED.  Instead, it communicates
!    with the user through an integer vector whose size is not specified.
!    You actually have to "ask" the routine to tell you the size of this
!    vector.  Then you can fill up the vector with values to be used to
!    influence the seeding of the random number routine.  The details of
!    how the seed affects the sequence are also unspecified, but the only
!    thing we are reasonably confident about is that specifying the same
!    seed should result in the same sequence, and specifying different
!    seeds should result in different sequences!
!
!    I assume this is far more than you wanted to know.  (It's certainly
!    more than I wanted to know!)
!
!    The argument SEED is an input quantity only, so it is legal
!    to type
!
!      call random_initialize ( 0 )
!
!    or
!
!      call random_initialize ( 18867 )
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED_INPUT, the user's "suggestion" for a seed.
!    However, if the input value is 0, the routine will come up with
!    its own "suggestion", based on the system clock.
!
  implicit none

  integer count
  integer count_max
  integer count_rate
  logical, parameter :: debug = .false.
  integer i
  integer seed
  integer ( kind = 4 ) seed_input
  integer, allocatable :: seed_vector(:)
  integer              seed_size
  real t
  integer, parameter :: warm_up = 100
!
!  ALL EXECUTABLE CODE COMMENTED OUT!!!
!
!!seed = seed_input
!
!  Initialize the random seed routine.
!
!!call random_seed ( )
!
!  Determine the size of the random number seed vector.
!
!!call random_seed ( size = seed_size )
!
!  Allocate a vector of the right size to be used as a random seed.
!
!!allocate ( seed_vector(seed_size) )
!
!  If the user supplied a SEED value, use that.
!
!  Otherwise, use the system clock value to make up a value that is
!  likely to change based on when this routine is called.
!
!!if ( seed /= 0 ) then

!!  if ( debug ) then
!!    write ( *, '(a)' ) ' '
!!    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
!!    write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, user SEED = ', seed
!!  end if

!!else

!!  call system_clock ( count, count_rate, count_max )

!!  seed = count

!!  if ( debug ) then
!!    write ( *, '(a)' ) ' '
!!    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
!!    write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, arbitrary SEED = ', &
!!      seed
!!  end if

!!end if
!
!  Set the seed vector.  We don't know the significance of the
!  individual entries of the internal seed vector, so we'll just set 
!  all entries to SEED.
!
!!seed_vector(1:seed_size) = seed
!
!  Now call RANDOM_SEED, and tell it to use this seed vector.
!
!!call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
!!deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times just to "warm it up".
!
!!do i = 1, warm_up
!!  call random_number ( harvest = t )
!!end do

  return
end
subroutine rat_factor ( m, n, maxfactor, factor_num, factor, power, &
  mleft, nleft )

!*****************************************************************************80
!
!! RAT_FACTOR factors a rational value into a product of prime factors.
!
!  Formula:
!
!    ( M / N ) = ( MLEFT / NLEFT ) * Product ( 1 <= I <= FACTOR_NUM )
!      FACTOR(I)**POWER(I).
!
!  Modified:
!
!    21 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the top and bottom of a rational value.
!    The ratio of M and N must be positive.
!
!    Input, integer ( kind = 4 ) MAXFACTOR, the maximum number of factors for
!    which storage has been allocated.
!
!    Output, integer ( kind = 4 ) FACTOR_NUM, the number of prime factors 
!    of M/N.
!
!    Output, integer ( kind = 4 ) FACTOR(MAXFACTOR), the prime factors of M/N.
!
!    Output, integer ( kind = 4 ) POWER(MAXFACTOR).  POWER(I) is the power of
!    the FACTOR(I) in the representation of M/N.
!
!    Output, integer ( kind = 4 ) MLEFT, NLEFT, the top and bottom of 
!    the factor of M / N that remains.  If ABS ( MLEFT / NLEFT ) is not 1, then
!    the rational value was not completely factored.
!
  implicit none

  integer ( kind = 4 ) maxfactor

  integer ( kind = 4 ) factor(maxfactor)
  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mleft
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nleft
  integer ( kind = 4 ) p
  integer ( kind = 4 ) power(maxfactor)
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) prime_max

  factor_num = 0

  mleft = m
  nleft = n
!
!  NLEFT should be nonnegative.
!
  if ( nleft < 0 ) then
    mleft = -mleft
    nleft = -nleft
  end if

  if ( m == 0 .or. n == 0 ) then
    return
  end if

  if ( m == n ) then
    factor_num = 1
    factor(1) = 1
    power(1) = 1
    return
  end if
!
!  Find out how many primes we stored.
!
  prime_max = prime ( -1 )

  do i = 1, prime_max

    p = prime ( i )

    if ( mod ( nleft, p ) == 0 .or. &
         mod ( abs ( mleft ), p ) == 0 ) then

      if ( factor_num < maxfactor ) then

        factor_num = factor_num + 1
        factor(factor_num) = p
        power(factor_num) = 0
!
!  Divide MLEFT by PRIME(I) as often as you can.
!
        if ( mod ( abs ( mleft ), p ) == 0  ) then

          do

            power(factor_num) = power(factor_num) + 1
            mleft = mleft / p

            if ( mod ( abs ( mleft ), p ) /= 0 ) then
              exit
            end if

          end do

        end if
!
!  Divide NLEFT by PRIME(I) as often as you can.
!
        if ( mod ( nleft, p ) == 0  ) then

          do

            power(factor_num) = power(factor_num) - 1
            nleft = nleft / p

            if ( mod ( nleft, p ) /= 0 ) then
              exit
            end if

          end do

        end if

        if ( power(factor_num) == 0 ) then
          factor_num = factor_num - 1
        end if

      end if

    end if

  end do

  return
end
subroutine rickey ( ab, bb, er, f, h, hb, hp, r, so, tb, g )

!*****************************************************************************80
!
!! RICKEY evaluates Branch Rickey's baseball index.
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Schwarz,
!    Looking Beyond the Batting Average,
!    The New York Times, 
!    Sunday, 1 August 2004.
!    
!    Branch Rickey,
!    Goodby to Some Old Baseball Ideas,
!    Life Magazine,
!    2 August 1954.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) AB, number of at-bats.
!
!    Input, integer ( kind = 4 ) BB, base on balls.
!
!    Input, integer ( kind = 4 ) ER, earned runs.
!
!    Input, real ( kind = 8 ) F, the fielding rating.
!
!    Input, integer ( kind = 4 ) H, number of hits.
! 
!    Input, integer ( kind = 4 ) HB, hit batsmen.
!
!    Input, integer ( kind = 4 ) HP, hit by pitcher.
!
!    Input, integer ( kind = 4 ) R, runs.
!
!    Input, integer ( kind = 4 ) SO, strike outs.
!
!    Input, integer ( kind = 4 ) TB, total bases.
!
!    Output, real ( kind = 8 ) G, the Branch Rickey index, an estimate for the
!    expected winning percentage of a team with the given statistics.
!    (0.5 has already been subtracted from this value.)
!
  implicit none

  integer ( kind = 4 ) ab
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) er
  real    ( kind = 8 ) f
  real    ( kind = 8 ) g
  integer ( kind = 4 ) h
  integer ( kind = 4 ) hb
  real    ( kind = 8 ) hitting
  integer ( kind = 4 ) hp
  real    ( kind = 8 ) pitching
  integer ( kind = 4 ) r
  integer ( kind = 4 ) so
  integer ( kind = 4 ) tb

  hitting = &
      real (  h + bb + hp,   kind = 8 ) / real ( ab + bb + hp, kind = 8 ) &
    + real ( 3 * ( tb - h ), kind = 8 ) / real ( 4 * ab,       kind = 8 ) &
    + real (              r, kind = 8 ) / real ( h + bb + hp,  kind = 8 )

  pitching = &
      real ( h,       kind = 8 ) / real ( ab,                   kind = 8 ) &
    + real ( bb + hb, kind = 8 ) / real ( ab + bb + hb,         kind = 8 ) &
    + real ( er,      kind = 8 ) / real ( h + bb + hb,          kind = 8 ) &
    - real ( so,      kind = 8 ) / real ( 8 * ( ab + bb + hb ), kind = 8 )

  g = hitting - pitching - f

  return
end
subroutine roots_to_r8poly ( n, x, c )

!*****************************************************************************80
!
!! ROOTS_TO_R8POLY converts polynomial roots to polynomial coefficients.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of roots specified.
!
!    Input, real ( kind = 8 ) X(N), the roots.
!
!    Output, real ( kind = 8 ) C(0:N), the coefficients of the polynomial.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) c(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x(n)
!
!  Initialize C to (0, 0, ..., 0, 1).
!  Essentially, we are setting up a divided difference table.
!
  c(0:n-1) = 0.0D+00
  c(n) = 1.0D+00
!
!  Convert to standard polynomial form by shifting the abscissas
!  of the divided difference table to 0.
!
  do j = 1, n
    do i = 1, n + 1 - j
      c(n-i) = c(n-i) - x(n+1-i-j+1) * c(n-i+1)
    end do
  end do

  return
end
subroutine roots_to_i4poly ( n, x, c )

!*****************************************************************************80
!
!! ROOTS_TO_I4POLY converts polynomial roots to polynomial coefficients.
!
!  Modified:
!
!    22 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of roots specified.
!
!    Input, integer ( kind = 4 ) X(N), the roots.
!
!    Output, integer ( kind = 4 ) C(0:N), the coefficients of the polynomial.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) x(n)
!
!  Initialize C to (0, 0, ..., 0, 1).
!  Essentially, we are setting up a divided difference table.
!
  c(0:n-1) = 0
  c(n) = 1
!
!  Convert to standard polynomial form by shifting the abscissas
!  of the divided difference table to 0.
!
  do j = 1, n
    do i = 1, n + 1 - j
      c(n-i) = c(n-i) - x(n+1-i-j+1) * c(n-i+1)
    end do
  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
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
!    05 February 2004
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
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
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements 
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Modified:
!
!    06 August 2005
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

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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
end
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = * ) string
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

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

  write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tuple_next2 ( n, xmin, xmax, x, rank )

!*****************************************************************************80
!
!! TUPLE_NEXT2 computes the next element of an integer tuple space.
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
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, integer ( kind = 4 ) XMIN(N), XMAX(N), the "minimum" and "maximum"
!    entry values.  These values are minimum and maximum only in the sense of 
!    the lexicographic ordering.  In fact, XMIN(I) may be less than, equal to,
!    or greater than XMAX(I).
!
!    Input/output, integer ( kind = 4 ) X(N), on input the previous tuple.
!    On output, the next tuple.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank of the item.  On first
!    call, set RANK to 0 to start up the sequence.  On return, if RANK is zero,
!    there are no more items in the sequence.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_one = 1
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) xmin(n)
  integer ( kind = 4 ) xmax(n)

  if ( rank < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RANK = ', rank
    stop
  end if

  if ( product ( 1 + abs ( xmax(1:n) - xmin(1:n) ) ) < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RANK = ', rank
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
      x(i) = x(i) + sign ( i4_one, xmax(i) - xmin(i) )
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
subroutine tvec_even ( nt, t )

!*****************************************************************************80
!
!! TVEC_EVEN computes evenly spaced angles between 0 and 2*PI.
!
!  Discussion:
!
!    The computation realizes that 0 = 2 * PI.
!
!  Example:
!
!    NT = 4
!
!    T = ( 0, PI/2, PI, 3*PI/2 )
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of values to compute.
!
!    Output, real ( kind = 8 ) TVEC(NT), the evenly spaced angles, in radians.
!
  implicit none

  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) t(nt)

  do i = 1, nt
    t(i) = real ( 2 * ( i - 1 ), kind = 8 ) * pi &
         / real ( nt, kind = 8 )
  end do

  return
end
subroutine tvec_even2 ( nt, t )

!*****************************************************************************80
!
!! TVEC_EVEN2 computes evenly spaced angles between 0 and 2*PI.
!
!  Discussion:
!
!    The computation realizes that 0 = 2 * PI.  The values are equally
!    spaced in the circle, do not include 0, and are symmetric about 0.
!
!  Example:
!
!    NT = 4
!
!    T = ( PI/4, 3*PI/4, 5*PI/4, 7*PI/4 )
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of values to compute.
!
!    Output, real ( kind = 8 ) TVEC(NT), the evenly spaced angles, in radians.
!
  implicit none

  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) t(nt)

  do i = 1, nt
    t(i) = real ( 2 * i - 1, kind = 8 ) * pi &
         / real ( nt, kind = 8 )
  end do

  return
end
subroutine tvec_even3 ( nt, t )

!*****************************************************************************80
!
!! TVEC_EVEN3 computes evenly spaced angles between 0 and 2*PI.
!
!  Discussion:
!
!    The angles begin with 0 and end with 2*PI.
!
!  Example:
!
!    NT = 4
!
!    T = ( 0, 2*PI/3, 4*PI/3 2*PI )
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of values to compute.
!
!    Output, real ( kind = 8 ) TVEC(NT), the evenly spaced angles, in radians.
!
  implicit none

  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) t(nt)

  if ( nt == 1 ) then
    t(1) = pi
  else
    do i = 1, nt
      t(i) = real ( 2 * ( i - 1 ), kind = 8 ) * pi &
           / real ( nt - 1, kind = 8 )
    end do
  end if

  return
end
subroutine tvec_even_bracket ( nt, theta1, theta2, t )

!*****************************************************************************80
!
!! TVEC_EVEN_BRACKET computes evenly spaced angles between THETA1 and THETA2.
!
!  Example:
!
!    NT = 4
!    THETA1 = 30
!    THETA2 = 90
!
!    T = ( 30, 50, 70, 90 )
!
!  Discussion:
!
!    The interval between THETA1 and THETA2 is divided into NT-1 subintervals.
!
!    The angles returned are the breakpoints of these subintervals,
!    including THETA1 and THETA2.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of values to compute.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the limiting angles.
!
!    Output, real ( kind = 8 ) TVEC(NT), the evenly spaced angles.
!
  implicit none

  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  real    ( kind = 8 ) t(nt)
  real    ( kind = 8 ) theta1
  real    ( kind = 8 ) theta2

  if ( nt == 1 ) then

    t(1) = ( theta1 + theta2 ) / 2.0D+00

  else

    do i = 1, nt
      t(i) = ( real ( nt - i,     kind = 8 ) * theta1   &
             + real (      i - 1, kind = 8 ) * theta2 ) &
             / real ( nt     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine tvec_even_bracket2 ( nt, theta1, theta2, t )

!*****************************************************************************80
!
!! TVEC_EVEN_BRACKET2 computes evenly spaced angles from THETA1 to THETA2.
!
!  Discussion:
!
!    The interval between THETA1 and THETA2 is divided into NT+1 subintervals.
!
!    The angles returned are the internal NT breakpoints of the subintervals.
!
!  Example:
!
!    NT = 5
!    THETA1 = 30
!    THETA2 = 90
!
!    T = ( 40, 50, 60, 70, 80 )
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of values to compute.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the limiting angles.
!
!    Output, real ( kind = 8 ) TVEC(NT), the evenly spaced angles.
!
  implicit none

  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  real    ( kind = 8 ) t(nt)
  real    ( kind = 8 ) theta1
  real    ( kind = 8 ) theta2

  do i = 1, nt
    t(i) = ( real ( nt + 1 - i, kind = 8 ) * theta1   &
           + real (          i, kind = 8 ) * theta2 ) &
           / real ( nt + 1,     kind = 8 )
  end do

  return
end
subroutine tvec_even_bracket3 ( nt, theta1, theta2, t )

!*****************************************************************************80
!
!! TVEC_EVEN_BRACKET3 computes evenly spaced angles from THETA1 to THETA2.
!
!  Discussion:
!
!    The interval between THETA1 and THETA2 is divided into NT subintervals.
!
!    The angles returned are the midpoints of each subinterval.
!
!  Example:
!
!    NT = 3
!    THETA1 = 30
!    THETA2 = 90
!
!    T = ( 40, 60, 80 )
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the number of values to compute.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the limiting angles.
!
!    Output, real ( kind = 8 ) TVEC(NT), the evenly spaced angles.
!
  implicit none

  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  real    ( kind = 8 ) t(nt)
  real    ( kind = 8 ) theta1
  real    ( kind = 8 ) theta2

  do i = 1, nt
    t(i) = ( real ( 2 * nt - 2 * i + 1, kind = 8 ) * theta1   &
           + real (          2 * i - 1, kind = 8 ) * theta2 ) &
           / real ( 2 * nt,             kind = 8 )
  end do

  return
end
subroutine upc_check_digit ( p, l, r, c )

!*****************************************************************************80
!
!! UPC_CHECK_DIGIT returns the check digit of a UPC.
!
!  Discussion:
!
!    UPC stands for Universal Price Code.
!
!    A full UPC is a string of 12 digits, in groups of size 1, 5, 5, and 1,
!    of the form P-LLLLL-RRRRR-C, where:
!
!      P is the one-digit product type code.
!      L is the five-digit manufacturer code.
!      R is the five_digit product code
!      C is the check digit.
!
!  Example:
!
!    0-72890-00011-8
!    0-12345-67890-5
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
!    Input, integer ( kind = 4 ) P, the one-digit product type code.
!
!    Input, integer ( kind = 4 ) L, the five-digit manufacturer code.
!
!    Input, integer ( kind = 4 ) R, the five-digit product code.
!
!    Output, intege ( kind = 4 )r C, the check digit.
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ), parameter :: i4_ten = 10
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc(5)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) r
  integer ( kind = 4 ) rc(5)

  if ( p < 0 .or. 9 < p ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UPC_CHECK_DIGIT - Fatal error!'
    write ( *, '(a)' ) '  P < 0 or 9 < P!'
    stop
  end if

  if ( l < 0 .or. 99999 < l ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UPC_CHECK_DIGIT - Fatal error!'
    write ( *, '(a)' ) '  L < 0 or 99999 < L!'
    stop
  end if

  if ( r < 0 .or. 99999 < r ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UPC_CHECK_DIGIT - Fatal error!'
    write ( *, '(a)' ) '  R < 0 or 99999 < R!'
    stop
  end if

  call i4_to_digits_decimal ( l, 5, lc )
  call i4_to_digits_decimal ( r, 5, rc )

  c = ( p + lc(2) + lc(4) + rc(1) + rc(3) + rc(5) ) * 3 &
          + lc(1) + lc(3) + lc(5) + rc(2) + rc(4)

  c = mod ( c, i4_ten )

  c = mod ( i4_ten - c, i4_ten )

  return
end

