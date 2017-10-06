
/*
disnan.c  Fortran LOGICAL function disnan(x) implemented in C to return .TRUE. 
if DOUBLE PRECISION x is a NaN, and .FALSE. otherwise.  

isnan.c  Fortran LOGICAL function isnan(x) implemented in C to return .TRUE. 
if REAL x is a NaN, and .FALSE. otherwise.  

second.c  Fortran DOUBLE PRECISION function implemented in C i
to return the CPU time in seconds since job start.  

setxp.c  Fortran REAL function setxp(x,n) implemented in C to return 
(fraction of x) * base**n, that is, to replace the exponent field of a copy of x.  

dsetxp for double

eps2.c  Fortran REAL function implemented in C to compute the generalized machine 
epsilon. This implementation works by fast low-level bit manipulation.  

deps2.c  Fortran DOUBLE PRECISION function implemented in C to 
compute the generalized machine epsilon. This implementation works by fast low-level bit manipulation.  

isinf.c  Fortran LOGICAL function isinf(x) implemented in C to return .TRUE. 
if REAL x is infinite, and .FALSE. otherwise.  

disinf.c  Fortran LOGICAL function disinf(x) implemented in C to return 
.TRUE. if DOUBLE PRECISION x is infinite, and .FALSE. otherwise.  

*/


/***********************************************************************
This file defines typedefs and symbols for interfacing C IEEE 754
support functions to Fortran code.
The model representation of a t-digit floating-point number is
	x = (-1)**s * 0.d1 d2 d3 ... dt * beta**e
where the digits dk satisfy
	0 <= dk < beta
The fractional part, which we call the significand, is defined to lie
in the range
	1/beta <= significand < 1
For IEEE floating-point, with its hidden bit and denormalized numbers,
we adjust parameters to conform to our model.  Denormalized numbers
are normalized by expanding their exponent range.
IEEE floating point arithmetic has these formats, where s is the sign
bit, e is an exponent bit, and f is a fraction bit.
Single precision:
	seee eeee efff ffff ffff ffff ffff ffff
	significand = 1.fff ffff ffff ffff ffff ffff (1 + 23 bits)
	exponent bias = 127
Double precision:
	seee eeee eeee ffff ffff ffff ffff ffff
	ffff ffff ffff ffff ffff ffff ffff ffff
	significand = 1.ffff ffff ffff ffff ffff ffff ffff ffff
			ffff ffff ffff ffff ffff (1 + 52 bits)
	exponent bias = 1023
Here are some sample IEEE bit patterns:
========================================================================
			    LITTLE ENDIAN
Single precision
	               0	0x00000000
	               1	0x3f800000
	              -1	0xbf800000
	               2	0x40000000
	              -2	0xc0000000
	     1.19209e-07	0x34000000	eps(1.0)
	    -1.19209e-07	0xb4000000
	     1.17549e-38	0x00800000	smallest normal
	    -1.17549e-38	0x80800000
	      1.4013e-45	0x00000001	smallest subnormal
	     -1.4013e-45	0x80000001
	   3.4028235e+38	0x7f7fffff	largest normal
	  -3.4028235e+38	0xff7fffff
	        Infinity	0x7f800000
	       -Infinity	0xff800000
	             NaN	0x7f80ffff
Double precision
	               0	0x00000000 00000000
	               1	0x00000000 3ff00000
	              -1	0x00000000 bff00000
	               2	0x00000000 40000000
	              -2	0x00000000 c0000000
	     1.11022e-16	0x00000002 3ca00000	eps(1.0)
	    -1.11022e-16	0x00000002 bca00000
	    2.22507e-308	0x00000000 00100000	smallest normal
	   -2.22507e-308	0x00000000 80100000
	    4.94066e-324	0x00000001 00000000	smallest subnormal
	   -4.94066e-324	0x00000001 80000000
	    1.79769e+308	0xffffffff 7fefffff	largest normal
	   -1.79769e+308	0xffffffff ffefffff
	        Infinity	0x00000000 7ff00000
	       -Infinity	0x00000000 fff00000
	             NaN	0xffffffff 7ff7ffff
========================================================================
			     BIG ENDIAN
Single precision
	               0	0x00000000
	               1	0x3f800000
	              -1	0xbf800000
	               2	0x40000000
	              -2	0xc0000000
	     1.19209e-07	0x34000000	eps(1.0)
	    -1.19209e-07	0xb4000000
	     1.17549e-38	0x00800000	smallest normal
	    -1.17549e-38	0x80800000
	      1.4013e-45	0x00000001	smallest subnormal
	     -1.4013e-45	0x80000001
	   3.4028235e+38	0x7f7fffff	largest normal
	  -3.4028235e+38	0xff7fffff
	             Inf	0x7f800000
	            -Inf	0xff800000
	             NaN	0x7fffffff

Double precision
	               0	0x00000000 00000000
	               1	0x3ff00000 00000000
	              -1	0xbff00000 00000000
	               2	0x40000000 00000000
	              -2	0xc0000000 00000000
	     1.11022e-16	0x3ca00000 00000002	eps(1.0)
	    -1.11022e-16	0xbca00000 00000002
	    2.22507e-308	0x00100000 00000000	smallest normal
	   -2.22507e-308	0x80100000 00000000
	    4.94066e-324	0x00000000 00000001	smallest subnormal
	   -4.94066e-324	0x80000000 00000001
	    1.79769e+308	0x7fefffff ffffffff	largest normal
	   -1.79769e+308	0xffefffff ffffffff
	             Inf	0x7ff00000 00000000
	            -Inf	0xfff00000 00000000
	             NaN	0x7fffffff ffffffff
========================================================================
***********************************************************************/

/* type mappings from Fortran to C */
typedef double DOUBLEPRECISION;

typedef float REAL;

#if defined(__alpha)
typedef int INTEGER;	/* need 32-bit integers, not 64-bit ones */
#else
typedef long INTEGER;
#endif

typedef int LOGICAL;	/* use int, not long, to avoid conflicts on HP-UX with */
			/* system header file declarations of isinf(), isnan() */
typedef union
{
    REAL r;
    INTEGER i;
} REAL_PARTS;				/* decomposition of REAL */

typedef union
{
    DOUBLEPRECISION r;
    INTEGER i[2];
} DOUBLEPRECISION_PARTS;		/* decomposition of DOUBLEPRECISION */

/* Fortran LOGICAL values -- compiler dependent!  Most (all?) UNIX */
/* Fortran compilers use 0 for .FALSE. and non-zero for .TRUE., like C. */

#ifdef _FALSE_
#undef _FALSE_
#endif

#ifdef _TRUE_
#undef _TRUE_
#endif

#define _FALSE_				((LOGICAL)0)
#define _TRUE_				((LOGICAL)1)

#define BASE				2 /* number base */

/* stored significand bits in single and double precision */

#define	T_SP				23
#define	T_DP				52

#define BASE_TO_THE_T_SP		((REAL)8388608.0)
#define BASE_TO_THE_T_DP		((DOUBLEPRECISION)4503599627370496.0)

#define _SHIFTED_EXPONENT_MASK_SP	0xff
#define EXPONENT_MASK_SP		0x7f800000L
#define EXPONENT_MASK_DP		0x7ff00000L

#define _SHIFTED_EXPONENT_MASK_DP	0x7ff
#define _SIGNIFICAND_MASK_SP		0x007fffffL
#define _SIGNIFICAND_MASK_DP		0x000fffffL

/* Exponent biases such that significand lies in (1/beta) <= significand < 1.
These are 1 less than the IEEE biases, because its stored significand
lies in 1 <= significand < beta due to the hidden bit.  We define them
with a leading underscore because they are for internal use only. */

#define _BIAS_SP			126
#define _BIAS_DP			1022

/* Indexes into two-word INTEGER array to account for addressing order */

#if i386 || sun386 || __i386__ || __sun386__ || msdos || MIPSEL || __alpha
					/* Intel 80xxx or MIPS little endian */
#define DP_LOW				0
#define DP_HIGH				1
#else				/* big endian (MIPS, Motorola, SPARC, ...) */
#define DP_LOW				1
#define DP_HIGH				0
#endif

/* macros to extract (high-order) significand and exponent as integer values */

#define GET_EXPONENT_SP(x) ((((x) >> T_SP) & \
				_SHIFTED_EXPONENT_MASK_SP) - _BIAS_SP)
#define GET_EXPONENT_DP(x) ((((x) >> (T_DP - 32)) & \
				_SHIFTED_EXPONENT_MASK_DP) - _BIAS_DP)

#define SET_EXPONENT_SP(x)	(((x) + _BIAS_SP) << T_SP)
#define SET_EXPONENT_DP(x)	(((x) + _BIAS_DP) << (T_DP - 32))

#define EXPONENT_DENORM_SP		(-_BIAS_SP)
#define EXPONENT_DENORM_DP		(-_BIAS_DP)

#define EXPONENT_INFNAN_SP		(255 - _BIAS_SP)
#define EXPONENT_INFNAN_DP		(2047 - _BIAS_DP)

#define SIGNIFICAND_SP(x)		(((x) & _SIGNIFICAND_MASK_SP))
#define SIGNIFICAND_DP(x)		(((x) & _SIGNIFICAND_MASK_DP))

#define MAX_NORMAL_SP			0x7f7fffffL
#define MAX_NORMAL_DP			0x7fefffffL
#define MAX_NORMAL_Low_DP		0xffffffffL

#define MIN_NORMAL_SP			0x00800000L
#define MIN_DENORMAL_SP			0x00000001L

#define MIN_NORMAL_DP			0x00100000L
#define MIN_NORMAL_Low_DP		0x00000000L
#define MIN_DENORMAL_DP			0x00000000L
#define MIN_DENORMAL_Low_DP		0x00000001L

#define Inf_SP				0x7f800000L
#define NegInf_SP			0xff800000L
#define NaN_SP				0x7fffffffL /* significand is */
						    /* arbitrary non-zero */

/* High-order words for double-precision Infinity and NaN. */
#define Inf_DP				0x7ff00000L
#define Inf_Low_DP			0x00000000L
#define NegInf_DP			0xfff00000L
#define NegInf_Low_DP			0x00000000L
#define NaN_DP				0x7fffffffL /* significand is */
#define NaN_Low_DP			0xffffffffL /* arbitrary non-zero */

#define ISNEG_SP(x)			((x) & 0x80000000L)
#define ISNEG_DP(x)			((x) & 0x80000000L)

#ifndef ABS
#define ABS(x)				(((x) < 0.0) ? -(x) : (x))
#endif

/* Map external names onto the conventions of the
local Fortran compiler.  */

#define adx	adx_
#define dadx	dadx_
#define deps	deps_
#define deps2	deps2_
#define dintxp	dintxp_
#define disden	disden_
#define disinf	disinf_
#define disnan	disnan_
#define dsetxp	dsetxp_
#define eps	eps_
#define eps2	eps2_
#define intxp	intxp_
#define isden	isden_
#define isinf	isinf_
#define isnan	isnan_
#define setxp	setxp_

#if defined(__STDC__) || defined(__cplusplus)
#define ARGS(plist)	plist
#define STDC		1
#define VOID_ARG	void
#else
#define const
#define ARGS(plist)	()
#define STDC		0
#define VOID_ARG
#endif

#if defined(__cplusplus)
extern "C" {
#endif

DOUBLEPRECISION	dadx ARGS((DOUBLEPRECISION *x__, INTEGER *n__));
DOUBLEPRECISION	deps ARGS((DOUBLEPRECISION *x__));
INTEGER		dintxp ARGS((DOUBLEPRECISION *x__));
DOUBLEPRECISION	dsetxp ARGS((DOUBLEPRECISION *x__, INTEGER *n__));

LOGICAL		disden ARGS((DOUBLEPRECISION *x__));
LOGICAL		disinf ARGS((DOUBLEPRECISION *x__));
LOGICAL		disnan ARGS((DOUBLEPRECISION *x__));

REAL		adx ARGS((REAL *x__, INTEGER *n__));
REAL		eps ARGS((REAL *x__));
INTEGER		intxp ARGS((REAL *x__));
REAL		setxp ARGS((REAL *x__, INTEGER *n__));

LOGICAL		isden ARGS((REAL *x__));
LOGICAL		isinf ARGS((REAL *x__));
LOGICAL		isnan ARGS((REAL *x__));

#if defined(__cplusplus)
};
#endif

#if STDC
#include <stdlib.h>
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif



#if STDC
LOGICAL
disnan(DOUBLEPRECISION *x)
#else /* NOT STDC */
LOGICAL
disnan(x)
DOUBLEPRECISION *x;
#endif /* STDC */
{
#if defined(__alpha)
    if (disden(x))
	return (_FALSE_);
    else
	return ((*x != *x) ? _TRUE_ : _FALSE_);
#else
    return ((*x != *x) ? _TRUE_ : _FALSE_);
#endif
}


#if STDC
LOGICAL
isnan(REAL *x)
#else /* NOT STDC */
LOGICAL
isnan(x)
REAL *x;
#endif /* STDC */
{
#if defined(__alpha)
    if (isden(x))
	return (_FALSE_);
    else
	return ((*x != *x) ? _TRUE_ : _FALSE_);
#else
    return ((*x != *x) ? _TRUE_ : _FALSE_);
#endif
}

#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>

#if defined(__alpha)
#include <time.h>
#define HZ	CLK_TCK
#endif


#if __convex__
#include <time.h>
#if !defined(HZ) && defined(CLOCKS_PER_SEC)
#define HZ CLOCKS_PER_SEC
#endif
#if !defined(HZ) && defined(CLK_TCK)
#define HZ CLK_TCK
#endif
#endif /* __convex__ */

#if STDC
#include <unistd.h>
#endif

#ifndef HZ
#ifdef AHZ
#define HZ AHZ				/* NeXT (CMU MACH O/S) uses AHZ */
#else
#define HZ 60
#endif
#endif

double
second(VOID_ARG)
{
    struct tms buffer;

    (void)times(&buffer);
    return((double)(buffer.tms_utime + buffer.tms_stime)/(double)(HZ));
}



#if STDC
REAL
setxp(REAL *x, INTEGER *n)
#else /* NOT STDC */
REAL
setxp(x,n)			/* return (fraction of x) * base**n */
REAL *x;			/* i.e. replace the exponent field of x */
INTEGER *n;
#endif /* STDC */
{
    REAL_PARTS w;

    if (isden(x))			/* then x is denormalized */
        w.r = *x * BASE_TO_THE_T_SP;	/* so must first normalize fraction */
    else
	w.r = *x;			/* x is normal, NaN, or Inf */
    if (w.r == 0.0)			/* zeroes must be left intact */
        /* NO-OP */;
    else if (isnan(x))			/* NaNs must be left intact */
        /* NO-OP */;
    else if (isinf(x))			/* Infs must be left intact */
        /* NO-OP */;
    else if (*n <= (EXPONENT_DENORM_SP - T_SP)) /* underflow to zero */
	w.r = (REAL)0.0;
    else if (*n <= EXPONENT_DENORM_SP)	/* generate denormalized value */
    {
	w.i &= ~EXPONENT_MASK_SP;	/* clear exponent field */
	w.i |= SET_EXPONENT_SP(*n + T_SP);	/* set new exponent */
					/* of scaled normal value */
	w.r /= BASE_TO_THE_T_SP;	/* and denormalize by division */
    }
    else if (*n >= EXPONENT_INFNAN_SP)	/* generate infinity */
	w.i = ISNEG_SP(w.i) ? NegInf_SP : Inf_SP;
    else				/* result is normal number */
    {					/* and exponent in range */
	w.i &= ~EXPONENT_MASK_SP;	/* clear exponent field */
	w.i |= SET_EXPONENT_SP(*n);	/* set new exponent */
    }

    return (w.r);
}


#if STDC
REAL
eps(REAL *x)
#else /* NOT STDC */
REAL
eps(x)
REAL	*x;
#endif /* STDC */
{
    REAL	half = 0.5;
    INTEGER	n = intxp(x) - T_SP;
    INTEGER	zero = 0;
    REAL_PARTS	w;

    if ( (*x == 0.0) || isden(x) )
    {
	w.i = MIN_DENORMAL_SP;
	return (w.r);
    }
    else if ( isnan(x) || isinf(x) )
	return *x;
    else if ( (*x < 0.0) && (setxp(x,&zero) == -half) )
    {					/* special boundary case */
	w.r = setxp(&half,(n--, &n));
	if (w.r == 0.0)			/* then x = -1.17549e-38 80800000 */
	    w.i = MIN_DENORMAL_SP;	/* and setxp() -> 0, so fix up */
	return (w.r);
    }
    else
	return (setxp(&half,&n));
}


#if STDC
LOGICAL
isinf(REAL *x)
#else /* NOT STDC */
LOGICAL
isinf(x)
REAL *x;
#endif /* STDC */
{
    REAL_PARTS w;

    w.r = *x;
    return (((GET_EXPONENT_SP(w.i) == EXPONENT_INFNAN_SP) &&
	     (SIGNIFICAND_SP(w.i) == 0))
	    ? _TRUE_ : _FALSE_);
}


#if STDC
DOUBLEPRECISION
deps(DOUBLEPRECISION *x)
#else /* NOT STDC */
DOUBLEPRECISION
deps(x)
DOUBLEPRECISION	*x;
#endif /* STDC */
{
    DOUBLEPRECISION		half = 0.5;
    INTEGER			n = dintxp(x) - T_DP;
    INTEGER			zero = 0;
    DOUBLEPRECISION_PARTS	w;

    if ( (*x == 0.0) || disden(x) )
    {
	w.i[DP_HIGH] = MIN_DENORMAL_DP;
	w.i[DP_LOW] = MIN_DENORMAL_Low_DP;
	return (w.r);
    }
    else if ( disnan(x) || disinf(x) )
	return *x;
    else if ( (*x < 0.0) && (dsetxp(x,&zero) == -half) )
    {					/* special boundary case */
	w.r = dsetxp(&half,(n--, &n));
	if (w.r == 0.0)		/* then x = 2.22507e-308 0x00000000 00100000 */
	{
	    w.i[DP_HIGH] = MIN_DENORMAL_DP; /* and setxp() -> 0, so fix up */
	    w.i[DP_LOW] = MIN_DENORMAL_Low_DP;
	}
	return (w.r);
    }
    else
	return (dsetxp(&half,&n));
}


#if STDC
DOUBLEPRECISION
dsetxp(DOUBLEPRECISION *x, INTEGER *n)
#else /* NOT STDC */
DOUBLEPRECISION
dsetxp(x,n)			/* return (fraction of x) * base**n */
DOUBLEPRECISION *x;		/* i.e. replace the exponent field of x */
INTEGER *n;
#endif /* STDC */
{
    DOUBLEPRECISION_PARTS w;

    if (disden(x))			/* then x is denormalized */
        w.r = *x * BASE_TO_THE_T_DP;	/* so must first normalize fraction */
    else
	w.r = *x;			/* x is normal, NaN, or Inf */
    if (w.r == 0.0)			/* zeroes must be left intact */
        /* NO-OP */;
    else if (disnan(x))		/* NaNs must be left intact */
        /* NO-OP */;
    else if (disinf(x))		/* Infs must be left intact */
        /* NO-OP */;
    else if (*n <= (EXPONENT_DENORM_DP - T_DP)) /* underflow to zero */
	w.r = (DOUBLEPRECISION)0.0;
    else if (*n <= EXPONENT_DENORM_DP)	/* generate denormalized value */
    {
	w.i[DP_HIGH] &= ~EXPONENT_MASK_DP; /* clear exponent field */
	w.i[DP_HIGH] |= SET_EXPONENT_DP(*n + T_DP); /* set new */
					/* exponent of scaled normal value */
	w.r /= BASE_TO_THE_T_DP;	/* and denormalize by division */
    }
    else if (*n >= EXPONENT_INFNAN_DP)	/* generate infinity */
    {
	w.i[DP_HIGH] = ISNEG_DP(w.i[DP_HIGH]) ? NegInf_DP : Inf_DP;
	w.i[DP_LOW] = ISNEG_DP(w.i[DP_HIGH]) ? NegInf_Low_DP : Inf_Low_DP;
    }
    else				/* result is normal number */
    {					/* and exponent in range */
	w.i[DP_HIGH] &= ~EXPONENT_MASK_DP;	/* clear exponent field */
	w.i[DP_HIGH] |= SET_EXPONENT_DP(*n);	/* set new exponent */
    }

    return (w.r);
}

#if STDC
LOGICAL
disinf(DOUBLEPRECISION *x)
#else /* NOT STDC */
LOGICAL
disinf(x)
DOUBLEPRECISION *x;
#endif /* STDC */
{
    DOUBLEPRECISION_PARTS w;

    w.r = *x;
    return (((GET_EXPONENT_DP(w.i[DP_HIGH]) == EXPONENT_INFNAN_DP) &&
	     (SIGNIFICAND_DP(w.i[DP_HIGH]) == 0) && (w.i[DP_LOW] == 0))
	    ? _TRUE_ : _FALSE_);
}




#if STDC
LOGICAL
isden(REAL *x)
#else /* NOT STDC */
LOGICAL
isden(x)			/* return (x is denormalized) */
REAL	*x;
#endif /* STDC */
{
    return (intxp(x) <= EXPONENT_DENORM_SP) ? _TRUE_ : _FALSE_;
}


#if STDC
INTEGER
intxp(REAL *x)
#else /* NOT STDC */
INTEGER
intxp(x) /* return the exponent of the base in the representation of x */
REAL *x;
#endif /* STDC */
{
    REAL_PARTS w;
    register INTEGER e;

    w.r = *x;
    e = GET_EXPONENT_SP(w.i);

    if (*x == 0.0)			/* handle zero specially */
	e = 0;
    else if (e == EXPONENT_DENORM_SP)	/* have denormalized number */
    {
	w.r *= BASE_TO_THE_T_SP;	/* make normal number */
	e = GET_EXPONENT_SP(w.i) - T_SP;
    }
    return (e);
}



#if STDC
INTEGER
dintxp(DOUBLEPRECISION *x)
#else /* NOT STDC */
INTEGER
dintxp(x) /* return the exponent of the base in the representation of x */
DOUBLEPRECISION *x;
#endif /* STDC */
{
    DOUBLEPRECISION_PARTS w;
    register INTEGER e;

    w.r = *x;
    e = GET_EXPONENT_DP(w.i[DP_HIGH]);

    if (*x == 0.0)			/* handle zero specially */
	e = 0;
    else if (e == EXPONENT_DENORM_DP)	/* have denormalized number */
    {
	w.r *= BASE_TO_THE_T_DP;	/* make normal number */
	e = GET_EXPONENT_DP(w.i[DP_HIGH]) - T_DP;
    }
    return (e);
}



#if STDC
DOUBLEPRECISION
deps2(DOUBLEPRECISION *x)
#else /* NOT STDC */
DOUBLEPRECISION
deps2(x)
DOUBLEPRECISION	*x;
#endif /* STDC */
{
    DOUBLEPRECISION		half = 0.5;
    INTEGER			n = dintxp(x) - T_DP;
    INTEGER			zero = 0;
    DOUBLEPRECISION_PARTS	w;

    if ( (*x == 0.0) || disden(x) )
    {
	w.i[DP_HIGH] = MIN_DENORMAL_DP;
	w.i[DP_LOW] = MIN_DENORMAL_Low_DP;
	return (w.r);
    }
    else if ( disnan(x) || disinf(x) )
	return *x;
    else if ( (*x < 0.0) && (dsetxp(x,&zero) == -half) )
    {					/* special boundary case */
	w.r = dsetxp(&half,(n--, &n));
	if (w.r == 0.0)		/* then x = 2.22507e-308 0x00000000 00100000 */
	{
	    w.i[DP_HIGH] = MIN_DENORMAL_DP; /* and setxp() -> 0, so fix up */
	    w.i[DP_LOW] = MIN_DENORMAL_Low_DP;
	}
	return (w.r);
    }
    else
	return (dsetxp(&half,&n));
}



#if STDC
LOGICAL
disden(DOUBLEPRECISION *x)
#else /* NOT STDC */
LOGICAL
disden(x)			/* return (x is denormalized) */
DOUBLEPRECISION	*x;
#endif /* STDC */
{
    return (dintxp(x) <= EXPONENT_DENORM_DP) ? _TRUE_ : _FALSE_;
}

