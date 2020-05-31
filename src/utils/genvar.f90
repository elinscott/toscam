module genvar

   use mpi !! External dependency
   use openmpmod

   implicit none

!include 'mpif.h'  --> FOR BG/L

!******************************************************************
!variables d utilite generale

   CHARACTER(LEN=100)                       :: LOGfile ! LOG FILE NAME

   integer, parameter                       :: SP = kind(1.0)    ! single precision type
   integer, parameter                       :: DP = kind(1.0d0)  ! double precision type
   integer, parameter                       :: DDP = 16          ! quad precision type
   integer, parameter                       :: LONG = selected_int_kind(18)

   integer                                  :: alternate_display_count = 0, alternate_display_cycle = 5
   logical                                  :: flag_alternate_display = .false., flag_introduce_only_noise_in_minimization = .false.
   logical                                   :: force_no_spline_plots = .false.
   integer(4), parameter                     :: int1 = 0
   integer(8), parameter                     :: int2 = 0
   real(kind=DP), parameter                  :: real1 = 0.d0
   real(16), parameter                       :: real2 = 0._16
   complex(kind=DP), parameter                     :: comp1 = 0.d0
   complex(16), parameter                    :: comp2 = 0._16

   integer                                  :: indepLink
   real(kind=DP), parameter                        :: errinout = 0.0001d0

   real(kind=SP), parameter                        :: pir = 3.14159265358979323846264338327950288419716939937510E+00
   real(kind=DP), parameter                  :: pi = 3.14159265358979323846264338327950288419716939937510d0
   real(16), parameter                       :: piq = 3.14159265358979323846264338327950288419716939937510_16

   real(kind=DP), parameter                        :: gamma_euler = 0.57721566490153286060d0  !Euler s constant
   logical                                  :: flag_introduce_noise_in_minimization = .false., no_mpi = .false.
   logical                                  :: messages = .false., messages2 = .false.
   logical                                  :: messages3 = .false., messages4 = .false.
   logical                                  :: flag_force_invmat_lapack = .false.
   logical                                  :: fast_invmat = .true.
   logical                                  :: use_cuda_routines = .false.
   logical                                  :: use_cula_routines = .false.
   logical                                  :: force_invmat_single_prec = .false.
   logical                                  :: use_openmp_invmat = .false.
   logical                                  :: plot_dislin_flag = .false., plot_dislin_3d_movie = .false.
   logical                                  :: diag_use_LU_instead_of_pivot = .false.
   logical                                  :: strongstop = .false.
   logical                                  :: enable_mpi_dot = .false.
   logical                                  :: flag_use_invmat_jordan = .false., flag_use_invmat_jordan_real = .false.

   integer                                  :: MAX_INT = huge(1)
   real(kind=DP), parameter                        :: MAX_REAL = huge(1.d0)
   real(kind=DP), parameter                        :: MAX_EXP = 700.d0
   real(kind=DP), parameter                        :: MIN_EXP = -700.d0
   real(16), parameter                       :: MAX_EXP_QUAD = 11000.d0
   real(16), parameter                       :: MIN_EXP_QUAD = -11000.d0
   real(kind=DP), parameter                        :: MAX_EXP_R = 81.d0
   real(kind=DP), parameter                        :: MIN_EXP_R = -81.d0
   real(kind=DP), parameter                        :: error = epsilon(1.d0)*100.d0
   real(kind=DP), parameter                        :: epsilonr = epsilon(1.d0), epsilonq = 1.d-30
   real(kind=DP), parameter                        :: euler_r = 2.71828182845904523536028747135266249&
                                                   &77572470936999595749669676277240766303535d0
   real(16), parameter                       :: euler_q = 2.71828182845904523536028747135266249&
                                                   &77572470936999595749669676277240766303535_16
   integer, parameter                        :: master = 0
   real(kind=DP)                                  :: ran_tab(0:10000)
   integer                                  :: iseed = 41845213
   logical                                  :: testing = .false.
   logical                                  :: MPIseparate = .false.
   complex(kind=DP), parameter                     :: imi = CMPLX(0.d0, 1.0d0, 8)
   complex(kind=DP), parameter                     :: jimi = CMPLX(-0.5d0, 0.5d0*SQRT(3.d0), 8)

   INTEGER, PARAMETER                 :: iwprint = 47         ! arbitrary frequency (to print debug)
   CHARACTER(LEN=1), DIMENSION(2), PARAMETER ::     pm = (/"+", "-"/)
   CHARACTER(LEN=2), DIMENSION(2), PARAMETER ::  cspin = (/" 1", "-1"/)
   CHARACTER(LEN=2), DIMENSION(2), PARAMETER :: ccspin = (/"up", "do"/)

   CHARACTER(LEN=9), PARAMETER             :: BOSONIC = 'BOSONIC  '
   CHARACTER(LEN=9), PARAMETER             :: FERMIONIC = 'FERMIONIC'
   CHARACTER(LEN=9), PARAMETER             :: MATSUBARA = 'MATSUBARA'
   CHARACTER(LEN=9), PARAMETER             :: RETARDED = 'RETARDED '
   CHARACTER(LEN=9), PARAMETER             :: ADVANCED = 'ADVANCED '

!******************************************************************


   REAL(DP), PARAMETER                 :: tiny_ = 1.D-10
   REAL(DP), PARAMETER                 :: huge_ = 1.D+10
   REAL(DP), PARAMETER                 :: quarter = 0.25_DP
   REAL(DP), PARAMETER                 :: third = 0.3333333333333_DP
   REAL(DP), PARAMETER                 :: half = 0.5_DP

   REAL(DP)                               :: pi2, oneoverpi, s2, s3, s6
   REAL(DP)                               :: one_over_clock_rate

   INTEGER                                 :: log_unit
   INTEGER, ALLOCATABLE                 :: ramp_proc(:), all_log_unit(:)

!******************************************************************

!MPI CONSTANTS

   INTEGER                                  :: nproc
   INTEGER                                  :: iproc
   INTEGER                                  :: ierr, nname
   CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME)    :: procname
   integer                                  :: rank, size2, myid, numprocs
   integer                                  :: status(MPI_STATUS_SIZE)

!******************************************************************

   ! Quality control testing
   logical :: running_qc_tests

end module

