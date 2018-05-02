module genvar

 use mpi
 use openmpmod

 implicit none

!include 'mpif.h'  --> FOR BG/L


!******************************************************************
!variables d utilite generale 

CHARACTER(LEN=100)                       :: LOGfile ! LOG FILE NAME 

integer                                  :: alternate_display_count=0,alternate_display_cycle=5
logical                                  :: flag_alternate_display=.false.,flag_introduce_only_noise_in_minimization=.false.
logical                                  :: force_no_spline_plots=.false.
integer(4),parameter                     :: int1=0
integer(8),parameter                     :: int2=0
real(8),parameter                        :: real1=0.d0
real(16),parameter                       :: real2=0._16
complex(8),parameter                     :: comp1=0.d0
complex(16),parameter                    :: comp2=0._16

integer,parameter                        :: wp      = kind (0d0)
integer,parameter                        :: int_k   = kind(12345678)
integer,parameter                        :: byte_k  = kind(.true.)
integer,parameter                        :: short_k = kind(12)
integer,parameter                        :: long_k  = kind(1234567891)

integer                                  :: indepLink
real(8),parameter                        :: errinout=0.0001d0

real(4),parameter                        :: pir = 3.14159265358979323846264338327950288419716939937510E+00
real(8),parameter                        :: pi  = 3.14159265358979323846264338327950288419716939937510d0
real(16),parameter                       :: piq = 3.14159265358979323846264338327950288419716939937510_16

real(8),parameter                        :: gamma_euler = 0.57721566490153286060d0  !Euler s constant
logical                                  :: flag_introduce_noise_in_minimization=.false.,no_mpi=.false.
logical                                  :: messages=.false.,messages2=.false.
logical                                  :: messages3=.false.,messages4=.false.
logical                                  :: flag_force_invmat_lapack=.false.
logical                                  :: fast_invmat=.true.
logical                                  :: use_cuda_routines=.false.
logical                                  :: use_cula_routines=.false.
logical                                  :: force_invmat_single_prec=.false.
logical                                  :: use_openmp_invmat=.false.
logical                                  :: plot_dislin_flag=.false.,plot_dislin_3d_movie=.false.
logical                                  :: diag_use_LU_instead_of_pivot=.false.
logical                                  :: strongstop=.false.
logical                                  :: enable_mpi_dot=.false.
logical                                  :: flag_use_invmat_jordan=.false.,flag_use_invmat_jordan_real=.false.

integer                                  :: MAX_INT  = huge(1) 
real(8),parameter                        :: MAX_REAL = huge(1.d0)
real(8),parameter                        :: MAX_EXP  =  700.d0
real(8),parameter                        :: MIN_EXP  = -700.d0
real(16),parameter                       :: MAX_EXP_QUAD = 11000.d0
real(16),parameter                       :: MIN_EXP_QUAD =-11000.d0
real(8),parameter                        :: MAX_EXP_R =  81.d0
real(8),parameter                        :: MIN_EXP_R = -81.d0
real(8),parameter                        :: error=epsilon(1.d0)*100.d0
real(8),parameter                        :: epsilonr=epsilon(1.d0),epsilonq=1.d-30
real(8),parameter                        :: euler_r= 2.71828182845904523536028747135266249&
                                                &77572470936999595749669676277240766303535d0
real(16),parameter                       :: euler_q= 2.71828182845904523536028747135266249&
                                                &77572470936999595749669676277240766303535_16
parameter(master=0)
real(8)                                  :: ran_tab(0:10000)
integer                                  :: iseed=41845213
logical                                  :: testing=.false.
logical                                  :: MPIseparate=.false.
complex(8),parameter                     :: imi  = CMPLX( 0.d0 , 1.0d0,8)
complex(8),parameter                     :: jimi = CMPLX(-0.5d0, 0.5d0*SQRT(3.d0),8) 


INTEGER,      PARAMETER                 :: iwprint=47         ! arbitrary frequency (to print debug)
CHARACTER(LEN=1),DIMENSION(2),PARAMETER ::     pm = (/"+","-"/)
CHARACTER(LEN=2),DIMENSION(2),PARAMETER ::  cspin = (/" 1","-1"/)
CHARACTER(LEN=2),DIMENSION(2),PARAMETER :: ccspin = (/"up","do"/)

CHARACTER(LEN=9), PARAMETER             :: BOSONIC    = 'BOSONIC  '
CHARACTER(LEN=9), PARAMETER             :: FERMIONIC  = 'FERMIONIC'
CHARACTER(LEN=9), PARAMETER             :: MATSUBARA  = 'MATSUBARA'
CHARACTER(LEN=9), PARAMETER             :: RETARDED   = 'RETARDED '
CHARACTER(LEN=9), PARAMETER             :: ADVANCED   = 'ADVANCED '

!******************************************************************

  INTEGER,      PARAMETER                 :: DBL=8,DP=8        ! "double" precision
  INTEGER,      PARAMETER                 :: DDP=16            ! "quad"   precision
  INTEGER,      PARAMETER                 :: SP = KIND(1.0)    ! "single" precision

  REAL(DBL),    PARAMETER                 :: tiny_   = 1.D-10
  REAL(DBL),    PARAMETER                 :: huge_   = 1.D+10
  REAL(DBL),    PARAMETER                 :: quarter = 0.25_DBL
  REAL(DBL),    PARAMETER                 :: third   = 0.3333333333333_DBL
  REAL(DBL),    PARAMETER                 :: half    = 0.5_DBL

  REAL(DBL)                               :: pi2,oneoverpi,s2,s3,s6
  REAL(DBL)                               :: one_over_clock_rate

  INTEGER                                 :: log_unit 
  INTEGER,    ALLOCATABLE                 :: ramp_proc(:),all_log_unit(:)

!******************************************************************

!MPI CONSTANTS

INTEGER                                  :: nproc
INTEGER                                  :: iproc
INTEGER                                  :: ierr,nname
CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME)    :: procname
integer                                  :: rank,size2,myid,numprocs
integer                                  :: status(MPI_STATUS_SIZE)

!******************************************************************


end module




