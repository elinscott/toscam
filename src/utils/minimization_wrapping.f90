module minimization_wrapping

   use genvar
   use mpi_mod
   use conj_grad, only: minimize
   use StringManip, only: tostring

   IMPLICIT NONE
   PRIVATE

   PUBLIC                           :: minimize_func_wrapper
   REAL(8), PARAMETER, PRIVATE     :: zero = 0.0_8, one = 1.0_8, two = 2.0_8, three = 3.0_8, four = 4.0_8
   LOGICAL, PARAMETER, PRIVATE   :: F = .FALSE., T = .TRUE.

contains

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE minimize_func_wrapper(func_, test, nnn, FIT_METH, Niter_search_max_, dist_min, dist_max, search_step, use_mpi, pow_in_)

      use random, only: random_float_from_interval
      use linalg, only: minloci

      implicit none
      CHARACTER(LEN=*)               :: FIT_METH
      INTEGER                        :: iparam, ii, nnn, Niter_search_max, Niter_search_max_
      INTEGER                        :: iexit, ifault, i, j
      real(kind=DP)                        :: dist_min, dist_max, search_step
      INTEGER                        :: mode, iprint
      REAL(8)                        :: dfn
      real(kind=DP), allocatable            :: min_ar(:), g(:), xtemp(:), xprmt(:), hess(:), w(:), wcivelli(:), WORK(:)
      real(kind=DP), allocatable            :: var(:), stepvec(:)
      real(kind=DP)                        :: test(:)
      LOGICAL                        :: use_mpi
      REAL(8)                        :: pow_in
      REAL(8), OPTIONAL               :: pow_in_
      real(kind=DP)                        :: dist_init

      INTERFACE
         SUBROUTINE func_(dist, n, vec)
            use genvar
            REAL(8), INTENT(OUT)   :: dist
            INTEGER, INTENT(IN)  :: n
            REAL(8), INTENT(IN)    :: vec(n)
         END SUBROUTINE
      END INTERFACE

!BUG: 2015 JUNE, w(..) size is wrong, it should be 4*nnn
!  allocate(g(nnn),xtemp(nnn),xprmt(nnn),hess(nnn**2),w(nnn**2),WORK(max(5*nnn+2,nnn*(nnn+7)/2)))
      allocate (g(nnn), xtemp(nnn), xprmt(nnn), hess(nnn**2), w(nnn**2), wcivelli(4*nnn), WORK(max(5*nnn + 2, nnn*(nnn + 7)/2)))

      allocate (var(nnn), stepvec(nnn))
      allocate (min_ar(size2))

      if (size(test) /= nnn) then
         write (*, *) 'ERROR - minimize_func_wrapper, dimension of test do no match the argument'
         write (*, *) 'size of test : ', size(test)
         write (*, *) 'nnn          : ', nnn
         stop
      endif

      min_ar = 0.d0

      Niter_search_max = Niter_search_max_

      write (*, *) '============================================='
      write (*, *) '============================================='
      write (*, *) '============================================='
      write (*, *) '---- MINIMIZATION WRAPPER -----'
      write (*, *) '          RANK=', rank
      write (*, *) ' METHOD         : ', FIT_METH
      write (*, *) ' NPARAM         : ', nnn
      write (*, *) ' size test      : ', size(test)
      write (*, *) ' NITER          : ', Niter_search_max
      write (*, *) ' dist max       : ', dist_max
      write (*, *) ' dist min       : ', dist_min
      write (*, *) ' use mpi        : ', use_mpi
      write (*, *) ' present(pow_in): ', present(pow_in_)
      write (*, *) ' noise  ?       : ', flag_introduce_noise_in_minimization
      write (*, *) ' only noise ?   : ', flag_introduce_only_noise_in_minimization
      write (*, *) ' allocated minar: ', allocated(min_ar)
      write (*, *) '============================================='
      write (*, *) '============================================='
      write (*, *) '============================================='

      if (rank == 0) write (*, *) 'SIZE OF MPI BATCH : ', size2

      if (size2 > 1 .and. .not. no_mpi .and. use_mpi) then
         call mpibarrier
      endif

      stepvec = 0.05d0
      iprint = 1
      mode = 1
      dist_min = 0.d0
      hess = 0.d0
      g = 0.d0
      xprmt = 0.d0
      w = 0.d0
      var = 0.d0

      if (rank == 0) write (*, *) 'VARIABLE RESETTED '
      if (size2 > 1 .and. .not. no_mpi .and. use_mpi) then
         call mpibarrier
      endif

      if (flag_introduce_noise_in_minimization) test = test + [(random_float_from_interval(-1.d-1, 1.d-1), i=1, size(test))]
      if (flag_introduce_only_noise_in_minimization) test = [(random_float_from_interval(-1.d0, 1.d0), i=1, size(test))]

      if (rank == 0) write (*, *) 'STARTING MINIMIZATION '

      if (size2 > 1 .and. .not. no_mpi .and. use_mpi) then
         call mpibarrier
      endif

      !---------------------------------------------------------------------------------------------!
      SELECT CASE (TRIM(ADJUSTL(FIT_METH)))
      CASE ('MINIMIZE')
         call init_minimize
    CALL minimize(func_, nnn, xtemp, dist_min, g, hess, w, dfn, xprmt, search_step, dist_max, mode, Niter_search_max, iprint, iexit)
         test = xtemp(1:nnn)
         ! CASE('CIVELLI')
         !    call init_minimize
         !    CALL minimize_civelli(distance_func______,nnn,xtemp,dist_min,g,hess,wcivelli,dfn,xprmt,search_step,dist_max,mode,Niter_search_max,iprint,iexit)
         !    test = xtemp(1:nnn)
      CASE DEFAULT
         write (*, *) 'error minimize wrapper case not found'
         write (*, *) 'MY RANK : ', rank
         stop
      END SELECT
      !---------------------------------------------------------------------------------------------!

      if (use_mpi .and. size2 > 1 .and. .not. no_mpi) then
         min_ar = 0.d0; min_ar(rank + 1) = dist_min; call mpisum(min_ar)
         write (*, *) 'FIT DONE FOR RANK : ', rank
         write (*, *) ' now waiting ..................... '
         call MPI_BCAST(test, size(test), MPI_DOUBLE_PRECISION, minloci(min_ar) - 1, MPI_COMM_WORLD, ierr)
      endif

      call func_(dist_min, nnn, test); 
      write (*, *) ' ------- '
      write (*, *) ' value at minimum is : ', dist_min
      write (*, *) ' ...end....'
      write (*, *) '============================================='
      write (*, *) '============================================='
      write (*, *) '============================================='

      deallocate (g, xtemp, xprmt, hess, w, wcivelli, WORK)
      deallocate (var, stepvec, min_ar)

   contains

      !---------------------!

      subroutine dump_screen_param_
         write (6, *) nnn
         write (6, *) xtemp(1:nnn)
         write (6, *) dist_min, dfn, search_step, dist_max, mode, Niter_search_max, iprint, iexit
         write (6, *) g(1:nnn)
         write (6, *) hess(1:nnn)
         write (6, *) w(1:nnn)
         write (6, *) xprmt(1:nnn)
      end subroutine

      !---------------------!

      subroutine distance_func_(n, m, u, f, con)
         implicit none
         integer, intent(in)  :: n, m
         real(kind=DP), intent(in)  :: u(n)
         real(kind=DP), intent(out) :: con(m), f
         call func_(f, n, u)
         con = 0.d0
      end subroutine

      !---------------------!

      real(kind=DP) function distance_func__(u)
         implicit none
         real(kind=DP) :: u(:)
         call func_(distance_func__, nnn, u(1:nnn))
      end function

      !---------------------!

      subroutine distance_func___(u, d)
         implicit none
         real(kind=DP), intent(in) :: u(:)
         real(kind=DP), intent(out) :: d
         call func_(d, nnn, u(1:nnn))
      end subroutine

      !---------------------!

      subroutine distance_func____(n, u, d)
         implicit none
         integer, intent(in) :: n
         real(kind=DP), intent(in) :: u(:)
         real(kind=DP), intent(out) :: d
         call func_(d, nnn, u(1:nnn))
      end subroutine

      !---------------------!

      SUBROUTINE distance_func_____(N, X, F, G)
         real(kind=DP) :: X(N), G(N), F, Fder
         integer :: N, i
         if (size(X) /= nnn) then
            write (*, *) '  N    : ', N
            write (*, *) 'size x : ', size(X)
            write (*, *) 'nnn    : ', nnn
            write (*, *) 'error minimization wrapper : distance_func_____'
            stop
         endif
         call func_(F, nnn, X(1:nnn))
         G = 0.d0
         do i = 1, nnn
            X(i) = X(i) + search_step
            call func_(Fder, nnn, X(1:nnn))
            X(i) = X(i) - search_step
            G(i) = (Fder - F)/search_step
         enddo
      END SUBROUTINE

      !---------------------!

      SUBROUTINE distance_func______(n, x, f)
         use genvar
         REAL(8), INTENT(OUT) :: f
         INTEGER, INTENT(IN)  :: n
         REAL(8), INTENT(IN)  :: x(n)
         call func_(f, n, x)
      END SUBROUTINE

      !---------------------!

      subroutine init_minimize
         implicit none
         write (*, *) 'INIT MINIMIZE - RANK', rank
         if (size2 > 1 .and. .not. no_mpi .and. use_mpi) then
            call mpibarrier
         endif
         iexit = 0
         dfn = -half
         xtemp = zero
         xtemp(1:nnn) = test
         xprmt = zero
         xprmt(1:nnn) = ABS(test) + 1.d-12
         write (*, *) 'calling initial distance'
         call func_(dist_init, nnn, test)
         write (*, *) 'initial distance is : ', dist_init
         !call dump_screen_param_
         write (*, *) 'DONE - RANK ', rank
         if (size2 > 1 .and. .not. no_mpi .and. use_mpi) then
            call mpibarrier
         endif
      end subroutine

      !---------------------!

   END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

end module
