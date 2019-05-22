
!======================!
!======================!

program dmft_standalone
!---------------!
   use dmft_solver_ed
   use genvar
   use init_and_close_my_sim
!---------------!
   implicit none

   call initialize_my_simulation

   testing = .false.
   fast_invmat = .true.
   ! use_cuda_routines=.false.
   ! use_cula_routines=.false.
   force_invmat_single_prec = .false.
   use_openmp_invmat = .false.
   diag_use_LU_instead_of_pivot = .false.
   flag_use_invmat_jordan = .true.
   flag_use_invmat_jordan_real = .true.
   enable_mpi_dot = .false.

   call stand_alone_ed
   call finalize_my_simulation

   write (*, *) 'ED_SOLVER_STANDALONE_FINISHED'

end program

!======================!
!======================!

