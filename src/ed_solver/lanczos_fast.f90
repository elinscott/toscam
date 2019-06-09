MODULE Lanczos_fast

   use genvar, only: DBL

   IMPLICIT NONE

   private

   REAL(DBL), PARAMETER, PRIVATE  ::  zero = 0.0_DBL, one = 1.0_DBL, rerror = 1.d-13

   ! DIAGONALIZATION USING PLAIN LANCZOS ALGORITHM TO GET LOWEST EIGENVALUE
   ! ONLY !

   public :: lanczos_fast_diagonalize
   public :: lanczos_get_gs_sector
   public :: one_step_lanczos_fast

contains

   subroutine Lanczos_fast_diagonalize(lowest)

      use common_def, only: c2s, dump_message, i2c, reset_timer, timer_fortran
      use eigen_class, only: add_eigen, delete_eigen, eigen_type, &
         eigenlist_type, is_eigen_in_window, new_eigen
      use genvar, only: huge_, iproc, log_unit, messages3, testing
      use globalvar_ed_solver, only: dEmax, neigen, nitermax, tolerance, &
         USE_TRANSPOSE_TRICK_MPI
      use h_class, only: dimen_H, sector_H, title_H_
      use rcvector_class, only: create_fix_initial_vector, delete_rcvector, &
         new_rcvector, norm_rcvector, rcvector_type
      use tridiag_class, only: delete_tridiag, diagonalize_tridiag, &
         new_tridiag, submatrix_tridiag, tridiag_type

      implicit none

      TYPE(eigenlist_type), INTENT(INOUT) :: lowest
      REAL(DBL), ALLOCATABLE, TARGET :: VECP(:, :), VALP(:)
      TYPE(eigen_type)               :: eigen(Neigen)
      TYPE(rcvector_type)            :: initvec, lastvec, tmp
      INTEGER                        :: iter, Niter, iter2, start_diagH
      REAL(DBL)                      :: rdist(Neigen), lasteigenval(Neigen), aa
      LOGICAL                        :: converged(Neigen), conv_one_step
      TYPE(tridiag_type)             :: Lmatrix, subLmatrix
      INTEGER                        :: istate, jstate, stateup, statedo, &
                                        Neigen_, jj
      REAL(DBL)                      :: coeff
      INTEGER                        :: dimenvec

      Niter = MIN(dimen_H(), Nitermax)
      Neigen_ = min(Neigen, dimen_H())
      Neigen_ = min(Neigen_, Niter)

      CALL reset_timer(start_diagH)

      ALLOCATE (VECP(Niter, Niter), VALP(Niter))
      lasteigenval = huge_

      CALL new_tridiag(Lmatrix, Niter) ! Lanczos matrix

      if (dimen_H() == 0) stop 'error Hilbert space has 0 dimension'

      if (.not. USE_TRANSPOSE_TRICK_MPI) then
         dimenvec = dimen_H()
      else
         dimenvec = (sector_h%updo%up%istatemax(iproc) - &
                     sector_h%updo%up%istatemin(iproc) + 1)* &
                    (sector_h%updo%down%istatemax(iproc) - &
                     sector_h%updo%down%istatemin(iproc) + 1)
         if (messages3) then
            write (log_unit, *) 'allocate Lanczos vector with size : ', dimenvec
            write (log_unit, *) 'up istatemin, istatemax : ', &
               sector_h%updo%up%istatemin(iproc), &
               sector_h%updo%up%istatemax(iproc)
            write (log_unit, *) 'dn istatemin, istatemax : ', &
               sector_h%updo%down%istatemin(iproc), &
               sector_h%updo%down%istatemax(iproc)
         endif
      endif

      CALL create_fix_initial_vector(initvec, dimenvec)
      CALL new_rcvector(lastvec, dimenvec)
      CALL new_rcvector(tmp, dimenvec)

      if (testing) then
         write (*, *) 'INITIAL VEC : ', initvec%rc
         write (*, *) 'dimenvec    : ', dimenvec
      endif

      DO iter = 1, Niter

         CALL one_step_Lanczos_fast(iter, initvec%rc, tmp%rc, lastvec%rc, &
                                    Lmatrix, conv_one_step)

         IF (norm_rcvector(lastvec) < tolerance) THEN
            CALL dump_message(TEXT="WARNING IN Lanczos_diagonalize: LANCZOS &
                 &LEAD TOWARDS NULL VECTOR!")
            WRITE (log_unit, '(2(a, E10.2))') "norm of the last Lanczos vector &
                 &= ", norm_rcvector(lastvec), " > tolerance = ", tolerance
            CALL flush(log_unit)
         ENDIF

         CALL submatrix_tridiag(subLmatrix, Lmatrix, (/1, iter/))
         CALL diagonalize_tridiag(subLmatrix, VALP(1:iter), VECP(1:iter, &
                                                                 1:iter), EIGENVAL_ONLY=.true.)
         CALL delete_tridiag(subLmatrix)
         rdist(1:min(Neigen_, iter)) = (VALP(1:min(Neigen_, iter)) - &
                                        lasteigenval(1:min(Neigen_, iter)))

         converged = ABS(rdist) < tolerance .and. iter > 2

         IF (ALL(converged(1:Neigen_)) .OR. iter == Niter .or. conv_one_step) &
            THEN
            write (log_unit, *) '.... fast lanczos converged ....'
            CALL submatrix_tridiag(subLmatrix, Lmatrix, (/1, iter/))
            CALL diagonalize_tridiag(subLmatrix, VALP(1:iter), VECP(1:iter, &
                                                                    1:iter), EIGENVAL_ONLY=.false.)
            CALL delete_tridiag(subLmatrix)

            do jj = 1, Neigen_
               if (jj == 1 .or. (converged(jj) .and. &
                                 is_eigen_in_window(VALP(jj), [VALP(1), VALP(1) + &
                                                               dEmax]))) then
                  CALL new_eigen(eigen(jj), VALP(jj), initvec, converged(jj), &
                                 RANK=jj, no_vector=.true.)
                  eigen(jj)%lanczos_iter = iter
                  eigen(jj)%lanczos_vecp(1:eigen(jj)%lanczos_iter) = &
                     VECP(1:iter, jj)
                  eigen(jj)%rdist = rdist(jj)
                  eigen(jj)%dim_space = dimen_H()
                  eigen(jj)%dim_sector = size(initvec%rc)
                  CALL add_eigen(eigen(jj), lowest)
               endif
            enddo
            EXIT
         ELSE
            lasteigenval(1:min(Neigen_, iter)) = VALP(1:min(Neigen_, iter))
         ENDIF

      ENDDO

      IF (ALLOCATED(VALP)) DEALLOCATE (VALP)
      IF (ALLOCATED(VECP)) DEALLOCATE (VECP)

      if (messages3) write (log_unit, *) ' Lanczos iterations .... delete &
           &vectors ....'

      do jj = 1, Neigen_
         CALL delete_eigen(eigen(jj))
      enddo

      CALL delete_rcvector(initvec)
      CALL delete_rcvector(lastvec)

      CALL delete_tridiag(Lmatrix)
      CALL delete_tridiag(subLmatrix)

      CALL timer_fortran(start_diagH, "# DIAGONALIZATION OF "// &
           TRIM(ADJUSTL(title_H_()))//" TOOK "//c2s(i2c(iter))//" &
           &ITERATIONS AND ")

   end subroutine

   subroutine Lanczos_get_GS_sector(lowest)

      use common_def, only: c2s, dump_message, i2c
      use eigen_class, only: eigen_allocate_vec, eigen_type, eigenlist_type
      use genvar, only: iproc, log_unit, testing
      use globalvar_ed_solver, only: nitermax, tolerance, USE_TRANSPOSE_TRICK_MPI
      use H_CLASS, only: dimen_H, sector_H
      use rcvector_class, only: create_fix_initial_vector, delete_rcvector, &
         new_rcvector, norm_rcvector, rcvector_type
      use tridiag_class, only: delete_tridiag, new_tridiag, tridiag_type
      use timer_mod, only: start_timer, stop_timer

      implicit none

      TYPE(eigenlist_type), intent(in) :: lowest
      TYPE(eigen_type), POINTER :: eigen => NULL()
      TYPE(rcvector_type)       :: initvec, lastvec, tmp
      TYPE(tridiag_type)        :: Lmatrix
      INTEGER                   :: iter, iter2, jj
      LOGICAL                   :: conv_one_step
      REAL(DBL)                 :: aa
      INTEGER                   :: dimenvec

      if (lowest%neigen == 0) return

      call start_timer("lanczos_get_gs_sector")

      write (log_unit, *) '====================================================='
      write (log_unit, *) '  GETTING THE GROUD STATE                            '
      write (log_unit, *) '  N eigenvalues    : ', lowest%neigen
      write (log_unit, *) '  iter convergence : ', lowest%eigen(:)%lanczos_iter
      write (log_unit, *) '  eigenvalues      : ', lowest%eigen(:)%val
      write (log_unit, *) '  SECTOR dimension : ', lowest%eigen(:)%dim_space
      write (log_unit, *) '  Lanczos vec. dim.: ', lowest%eigen(:)%dim_sector
      write (log_unit, *) '====================================================='

      if (.not. USE_TRANSPOSE_TRICK_MPI) then
         dimenvec = dimen_H()
      else
         dimenvec = (sector_h%updo%up%istatemax(iproc) - &
                     sector_h%updo%up%istatemin(iproc) + 1)* &
                    (sector_h%updo%down%istatemax(iproc) - &
                     sector_h%updo%down%istatemin(iproc) + 1)
      endif

      CALL new_rcvector(lastvec, dimenvec)

      CALL new_tridiag(Lmatrix, maxval(lowest%eigen(:)%lanczos_iter))
      CALL new_rcvector(tmp, dimenvec)

      do jj = 1, lowest%neigen

         eigen => lowest%eigen(jj)

         iter = eigen%lanczos_iter

         CALL create_fix_initial_vector(initvec, dimenvec)
         if (testing) then
            write (*, *) 'INITIAL VEC : ', initvec%rc
            write (*, *) 'dimenvec :', dimenvec
         endif
         CALL eigen_allocate_vec(eigen, initvec)

         aa = norm_rcvector(initvec)
         if (abs(aa) < rerror) aa = rerror
         aa = 1.d0/aa
         eigen%vec%rc = eigen%lanczos_vecp(1)*initvec%rc*aa
         CALL one_step_Lanczos_fast(1, initvec%rc, tmp%rc, lastvec%rc, &
                                    Lmatrix, conv_one_step)
         if (.not. conv_one_step) then
            DO iter2 = 2, iter - 1
               aa = norm_rcvector(lastvec)
               if (abs(aa) < rerror) aa = rerror
               aa = 1.d0/aa
               lastvec%rc = lastvec%rc*aa
               eigen%vec%rc = eigen%vec%rc + eigen%lanczos_vecp(iter2)* &
                              lastvec%rc
               CALL one_step_Lanczos_fast(iter2, initvec%rc, tmp%rc, &
                                          lastvec%rc, Lmatrix, conv_one_step)
            ENDDO
            aa = norm_rcvector(lastvec)
            if (abs(aa) < rerror) aa = rerror
            aa = 1.d0/aa
            lastvec%rc = lastvec%rc*aa
            eigen%vec%rc = eigen%vec%rc + eigen%lanczos_vecp(iter)*lastvec%rc
            aa = norm_rcvector(eigen%vec)
            if (abs(aa) < rerror) aa = rerror
            aa = 1.d0/aa
            eigen%vec%rc = eigen%vec%rc*aa
         endif

         IF (iter == Nitermax) THEN
            WRITE (log_unit, '(2(a, E10.2))') "Accuracy not met... relative &
                 &distance = ", eigen%rdist, " > tolerance = ", tolerance
            CALL flush(log_unit)
            CALL dump_message(TEXT="Increase Nitermax : Lanczos matrix size &
                 &= "//c2s(i2c(Nitermax)))
         ENDIF

      enddo

      CALL delete_rcvector(initvec)
      CALL delete_rcvector(lastvec)
      CALL delete_rcvector(tmp)
      CALL delete_tridiag(Lmatrix)

      call stop_timer("lanczos_get_gs_sector")

   end subroutine

   subroutine one_step_Lanczos_fast(iter, vec_in, tmp, vec_out, tri, &
                                    converged_in_one_step)

      use common_def, only: c2s, dump_message, i2c
      use genvar, only: log_unit
      use globalvar_ed_solver, only: USE_TRANSPOSE_TRICK_MPI, DO_NOT_USE_OPT_LANCZOS
      use mpi_mod, only: mpi_dot_product, split
      use tridiag_class, only: tridiag_type

      implicit none

#ifdef _complex
      COMPLEX(DBL), INTENT(INOUT)       :: vec_out(:), tmp(:), vec_in(:)
#else
      REAL(DBL), INTENT(INOUT)          :: vec_out(:), tmp(:), vec_in(:)
#endif
      TYPE(tridiag_type), INTENT(INOUT) :: tri
      REAL(DBL), SAVE :: normv_in = 0.d0, normv_out = 0.d0
      INTEGER         :: iter, i
      INTEGER         :: start_Hmult, ist
      LOGICAL         :: verbose
      LOGICAL         :: converged_in_one_step
      REAL(DBL)       :: aa, ddiag

      IF (iter == 1) THEN
         normv_in = SQRT(ABS(MPI_DOT_PRODUCT(vec_in, vec_in, split= &
                                             USE_TRANSPOSE_TRICK_MPI)))
         normv_out = 0.d0
         tri%diag = 0.0_DBL
         tri%subdiag = 0.0_DBL
         tri%subdiag(1) = normv_in
      ELSE
         IF (iter > tri%N) THEN
            CALL dump_message(TEXT="ERROR IN one_step_Lanczos_fast: TOO MANY &
                 &ITERATIONS!")
            CALL dump_message(TEXT="Lanczos iteration = "//c2s(i2c(iter)) &
                              //" > Lanczos matrix size = "//c2s(i2c(tri%N)))
            STOP 'CRITICAL'
         ENDIF
      ENDIF

      verbose = .false.
      if (verbose .or. mod(iter, 100) == 0) write (log_unit, *) 'iter Lanczos : &
           &', iter

      if (.not. USE_TRANSPOSE_TRICK_MPI .and. .not. DO_NOT_USE_OPT_LANCZOS) then
         call lanc_opt
      else
         call lanc_simp
      endif

      converged_in_one_step = abs(normv_out) < 1.d-24 .and. iter == 1

   contains

#include "Lanczos_plain.h"

   end subroutine

end module
