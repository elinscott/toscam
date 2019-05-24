!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

include 'dmft_one_iteration_split_module.h'

!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

program dmftonetepsplit
   use dmft_split_variables
   use openmpmod, only: init_openmp
   use genvar,    only: iproc
   implicit none
   integer                    :: i, j, k, ww
   character(200)             :: files, filename_sigma_source, filename_sigma
   character(200), allocatable :: list_of_greens(:, :)
   logical                    :: check, flag_onetep_producing_only_up_spin
   type(string)               :: cc_
   integer, allocatable        :: list_uniform(:, :)

   !-------------------------------!
   call init_my_input_variables
   call initialize_my_simulation_
   call init_openmp(silent = (iproc /= 1))
   !-------------------------------!

   write (*, *) 'RUNNING THE DMFT CALCULATIONS WITH [x] cpus : ', size2
   if (uniform_sigma .and. size2 > 1) then
      write (*, *) 'USING MPI WITH UNIFORM_SIGMA FLAG IS A WASTE OF RESSOURCE'
      write (*, *) 'MPI PARALLELIZATION IS ONLY EFFICIENT WHEN SIGMA IS NON-UNIFORM'
      write (*, *) 'PLEASE CHANGE INPUTS ACCORDINGLY'
      stop
   endif

   call build_list_of_files

   call mpibarrier

   call loop_over_atoms

   call copy_files

   call mpibarrier

   !-------------------------------!
   if (reboot_mpi) call finalize_my_simulation
   !-------------------------------!

contains

   include 'dmft_one_iteration_split.h'

end program

!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

