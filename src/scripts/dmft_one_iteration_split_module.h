!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

module dmft_split_variables
   use namelistmod, only: namelist_set, namelist_init, putel_in_namelist, &
                                    & look_for_namelist_in_file, look_for_command_line_argument
   use StringManip, only: toString, StrInt2
   use genvar, only: size2, rank, ierr
   use init_and_close_my_sim, only: initialize_my_simulation, finalize_my_simulation
   use mpi_mod, only: mpibarrier
   use strings, only: replace_in_string, string, assignment(=)
   implicit none
   integer        ::  paramagnetic, iter_dmft, niter_dmft, rankin, size2in
   logical        ::  uniform_sigma, onetep_only_up, flag_donot_keep_all_files, reboot_mpi
   character(200) ::  input_temp_dir

contains

   !-------------------------!
   subroutine init_my_input_variables
      implicit none
      type(namelist_set) :: nm
      call namelist_init(nm, 200, name_of_namelist='dmftonetep_variables_split')
      call putel_in_namelist(nm, niter_dmft, 'niter_dmft', 5, 'Definition:Number of DMFT iteration')
      call putel_in_namelist(nm, paramagnetic, 'paramagnetic', 1, 'Definition:=1 is paramagnetic,=0 is non-paramagnetic')
   call putel_in_namelist(nm,uniform_sigma,'uniform_sigma',.false.,'Definition:=true sigma is imposed as spatially uniform, so only one dmft calculation is run,false every dmft is done')
   call putel_in_namelist(nm,flag_donot_keep_all_files,'flag_donot_keep_all_files',.false.,'Definition: if true will not keep all the details and every files of the calculations')
   call putel_in_namelist(nm,onetep_only_up,'onetep_only_up',.false.,'Definition:=true if onetep produces only up spin, keep false if you dont know') 
      call putel_in_namelist(nm, iter_dmft, 'iter_dmft', -1, 'Definition:=dmft iteration number')
   call putel_in_namelist(nm,input_temp_dir,'input_temp_dir',".",'Definition:=local directory to store the atom subdirectories, which contains the dmft calculations for each individual correlated atoms. Since there is a lot of writing-reading involved, when running on a cluster with shared nfs-file system, it might be a good idea to indicate a local directory on the nodes such as /tmp/dmft_XXX, take care to give an individual name for this dir related to the job number, such that 2 jobs dont collide')
   call putel_in_namelist(nm,rankin,'rankin',0,'Definition:=rank given as an input, in this case this module does not boot mpi, but uses an existing running mpi structure')
   call putel_in_namelist(nm,size2in,'size2in',0,'Definition:=size2 given as an input, in this case this module does not boot mpi, but uses an existing running mpi structure')
      call look_for_namelist_in_file(nm, './input_onetep_dmft.txt')
      call look_for_command_line_argument(nm)

      if (size2in > 0) then
         reboot_mpi = .false.
      else
         reboot_mpi = .true.
      endif
      if (.not. reboot_mpi) then
         write (*, *) '###### splitting Hubbard atoms : do not reboot MPI ######'
         rank = rankin
         size2 = size2in
      else
         write (*, *) '###### splitting Hubbard atoms :        reboot MPI ######'
      endif

   end subroutine
   !-------------------------!

end module

!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
