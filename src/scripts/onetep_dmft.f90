!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

module dmft_variables
   use namelistmod, only: namelist_set, namelist_init, putel_in_namelist, &
                                    & look_for_namelist_in_file, look_for_command_line_argument
   use strings, only: replace_in_string, string, assignment(=)
   use StringManip, only: toString
   use init_and_close_my_sim, only: initialize_my_simulation, finalize_my_simulation

   integer            ::  mpi_onetep_type, restart_from_iteration_number, ed_real_frequ_last, niter_dmft, &
                          nproc, nproc_onetep, nproc_mpi_solver, openmp_solver
   character(2000)    ::  CASE_ONETEP, dir_onetep, exec_onetep, dir_onetep_mpi
   character(60)      ::  mach_onetep
   real(8)            ::  ed_frequ_min, ed_frequ_max
logical            ::  split_onetep,start_from_an_old_sim,all_local_host,just_onetep,compute_dos,numa,dmft_split,dmft_splitkdmftall,nomachinefile
   integer            ::  dmft_splitk_batch, nproc_onetep_openmp_
   logical            ::  hide_errors
   type(namelist_set) ::  nm

contains

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

 character(2000) function build_mpi_command_line(prefix,np,p,omp,mach,mach_arg,args,exe,outputin,hide_errors,localhost,ampersand)
      implicit none
      integer         ::  np, omp, p
      character*(*)   ::  mach, mach_arg, args, exe, outputin, prefix, localhost
      character(2000) ::  output
      logical         ::  hide_errors
      logical, optional::  ampersand

      output = trim(adjustl(outputin))
      if (hide_errors) then
         output = trim(adjustl(output))//" 2> /dev/null "
      else
         output = trim(adjustl(output))//" 2>&1 "
      endif

      write (*, *) 'CALLING get_mpi_command_line with arguments : '
      write (*, *) trim(adjustl(exe))//' "   '&
                                         & //trim(adjustl(args))//' "   '&
                                         & //trim(adjustl(tostring(np)))//'     '&
                                         & //trim(adjustl(tostring(p)))//'     '&
                                         & //trim(adjustl(tostring(omp)))//'     '&
                                         & //trim(adjustl(mach))//' "   '&
                                         & //trim(adjustl(mach_arg))//' " " '&
                                         & //trim(adjustl(output))//' "   '&
                                         & //trim(adjustl(localhost))

      build_mpi_command_line = " "
      call system('  get_mpi_command_line '//trim(adjustl(exe))//' "     '&
                                         & //trim(adjustl(args))//' "     '&
                                         & //trim(adjustl(tostring(np)))//'       '&
                                         & //trim(adjustl(tostring(p)))//'       '&
                                         & //trim(adjustl(tostring(omp)))//'       '&
                                         & //trim(adjustl(mach))//' "     '&
                                         & //trim(adjustl(mach_arg))//' " "   '&
                                         & //trim(adjustl(output))//' "     '&
                                         & //trim(adjustl(localhost))//'  >  __mpi_cmd__  ')
      open (unit=9911, file='__mpi_cmd__')
      read (9911, '(a)') build_mpi_command_line
      close (9911)
      build_mpi_command_line = trim(adjustl(prefix))//" "//trim(adjustl(build_mpi_command_line))
      if (present(ampersand)) then
         build_mpi_command_line = trim(adjustl(build_mpi_command_line))//"  &   "
      endif

      call system(" rm __mpi_cmd__ ")

   end function

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

   subroutine replace_it(string1, string2)
      implicit none
      character*(*) :: string1, string2
      call system(" sed 's/"//string1//"/"//string2//"/' *.dat > scratch ")
      call system(" mv scratch "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
   end subroutine

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

   subroutine flip_flag(sstring, true_to_false)
      implicit none
      character*(*)  :: sstring
      character(200) :: match1, match2
      logical        :: true_to_false
      type(string)   :: cc_
      write (*, *) 'looking for strings with : ', sstring
      call system("grep_in_dat "//TRIM(ADJUSTL(sstring)))
      write (*, *) 'opening file scratch'
      open (unit=10002, file='scratch')
      do
         read (10002, '(a)', end=66) match1
         cc_ = trim(adjustl(match1))
         write (*, *) 'flipping flag : ', TRIM(ADJUSTL(match1))
         if (true_to_false) then
            call replace_in_string(cc_, ' T ', ' F ', 'last')
            call replace_in_string(cc_, ' TRUE', ' FALSE ', 'last')
         else
            call replace_in_string(cc_, ' F ', ' T ', 'last')
            call replace_in_string(cc_, ' FALSE', ' TRUE ', 'last')
         endif
         match2 = cc_
         write (*, *) 'modified string : ', TRIM(ADJUSTL(match2))
         call system(" sed 's/"//TRIM(ADJUSTL(match1))//"/"//TRIM(ADJUSTL(match2))//"/' *.dat > scratch2 ")
         call system(" mv scratch2 "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
      enddo
66    continue
      close (10002); 
      call system("rm scratch  > /dev/null 2>&1")
      call system("rm scratch2 > /dev/null 2>&1 ")
   end subroutine

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

   subroutine my_dmft_exe
      implicit none
      call system("get_onetep_env")
      open (unit=10001, file='dir_onetep')
      read (10001, '(a)') dir_onetep
      write (*, *) 'My TOSCAM base directory is ', TRIM(ADJUSTL(dir_onetep))
      close (10001)
      call system("rm dir_onetep > /dev/null 2>&1")

      call system("get_mpi_env")
      open (unit=10001, file='dir_onetep_mpi')
      read (10001, '(a)') dir_onetep_mpi
      write (*, *) 'My MPI executable is ', TRIM(ADJUSTL(dir_onetep_mpi))
      close (10001)
      call system("rm dir_onetep_mpi > /dev/null 2>&1")

      call system("get_onetep_exec")
      open (unit=10001, file='exec_onetep')
      read (10001, '(a)') exec_onetep
      write (*, *) 'My ONETEP executable is ', TRIM(ADJUSTL(exec_onetep))
      close (10001)
      call system("rm exec_onetep > /dev/null 2>&1")
   end subroutine

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

   subroutine check_if_input_file_changed
      implicit none
      call look_for_namelist_in_file(nm, './input_onetep_dmft.txt')
   end subroutine

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

   subroutine init_my_input_variables
      implicit none
      integer :: numb_dat

      call namelist_init(nm, 200, name_of_namelist='dmftonetep_variables')
      call putel_in_namelist(nm, niter_dmft, 'niter_dmft', 5, 'Definition:Number of DMFT iteration')
      call putel_in_namelist(nm, nproc, 'nproc', 1, 'Definition:=number of processors for dmft')
   call putel_in_namelist(nm,nproc_mpi_solver,'nproc_mpi_solver',1,'Definition:=number of processors (MPI type) used for DMFT solver')
   call putel_in_namelist(nm,nproc_onetep_openmp_,'nproc_onetep_openmp_',1,'Definition: will run mpi for the onetep dmft interface, but each mpi job will run with nproc_onetep_openmp openmp threads')
      call putel_in_namelist(nm, nproc_onetep, 'nproc_onetep', 1, 'Definition:=number of processors (MPI type) for onetep')
     call putel_in_namelist(nm, openmp_solver, 'openmp_solver', 8, 'Definition:Number of open-mp cores running for the dmft-solver')
      call putel_in_namelist(nm, split_onetep, 'split_onetep', .false., 'Definition:=asynchronous calls to onetep')
   call putel_in_namelist(nm,start_from_an_old_sim,'start_from_an_old_sim',.false.,'Definition:=if true start from an old calculation, sigma files must be present')
    call putel_in_namelist(nm, all_local_host, 'all_local_host', .false., 'Definition:=if true will send all the job on local host')
   call putel_in_namelist(nm,just_onetep,'just_onetep',.false.,'Definition:=if true it will only compute the onetep part, no dmft at all, usefull to compute the DOS after a converged DMFT with a matsubara solver')
   call putel_in_namelist(nm,compute_dos,'compute_dos',.false.,'Definition:=if true will run onetep only (1 iteration) to compute the density of states from a real frequency Solver, or from ED if the run reached the last iteration and if last_iter_real=.true. was set during the DMFT calculation. There is no need to adapt the number of dmft points or the dmft temperature ih the onetep file')
   call putel_in_namelist(nm,ed_real_frequ_last ,'ed_real_frequ_last', 100  , 'Definition:ed number of frequencies for last DMFT iter, going back to onetep for FULL_DOS')
      call putel_in_namelist(nm, ed_frequ_min, 'ed_frequ_min ', -10.d0, 'Definition: min ED real frequ')
      call putel_in_namelist(nm, ed_frequ_max, 'ed_frequ_max ', +10.d0, 'Definition: max ED real frequ')
   call putel_in_namelist(nm,restart_from_iteration_number,'restart_from_iteration_number',0,'Definition:restart from iteration number given here instead of starting from first iteration, note that you might need to feed in the right sigma_output files in the present directory, just copy files from dir_onetepX before_onetep/ to present directory')
   call putel_in_namelist(nm, mpi_onetep_type, 'mpi_onetep_type', 1, 'Definition: if 1 will run a home made mpi, if 2 will run a more regular mpi, hand made is to use in difficult situations when clusters are giving some connection troubles')
     call putel_in_namelist(nm, numa, 'numa', .false., 'Definition: numa avoids remote memory access in multi-socket architectures')
   call putel_in_namelist(nm,dmft_split,'dmft_split',.false.,'Definition: if true splits the mpi onetep dmft interface and each cpus computes a frequency, instead of running together and inverting matrices by SCALAPACK')
   call putel_in_namelist(nm,dmft_splitkdmftall,'dmft_splitkdmftall',.false.,'Definition: if true splits the mpi onetep dmft interface in several batches each of them running different K points')
   call putel_in_namelist(nm,dmft_splitk_batch,'dmft_splitk_batch',1,'Definition: number of cpus in each of the batch when splitting k points')
   call putel_in_namelist(nm,hide_errors,'hide_errors',.true.,'Definition: this iterface will not show the list of minor errors occuring when for instance the mpi aborts')
   call putel_in_namelist(nm,mach_onetep,'mach_onetep','machines_onetep','Definition: name of the machine file for the GF calculation')
 call putel_in_namelist(nm, nomachinefile, 'nomachinefile', .false., 'Definition: if true removes -machinefile from the mpi syntax')

      call look_for_namelist_in_file(nm, './input_onetep_dmft.txt')
      call look_for_command_line_argument(nm)

      call system("get_onetep_case")

      open (unit=10002, file='case_onetep')
      read (10002, *) CASE_ONETEP
      close (10002)
      write (*, *) 'my ONETEP CASE is : ', trim(adjustl(CASE_ONETEP))

      call system(" ls *.dat | wc -l >> other_dat_file_  ")
      open (unit=118229, file='other_dat_file_')
      read (118229, *) numb_dat
      close (118229)
      call system(" rm other_dat_file_ > /dev/null 2>&1 ")
      if (numb_dat > 1) then
         write (*, *) 'PLEASE ONLY KEEP ONE *.dat FILE IN THE RUNNING DIRECTORY'
         write (*, *) 'ERROR,CRITICAL'
         stop
      endif

      call check_flag_consistency

      if (.not. just_onetep) then
         if (.not. start_from_an_old_sim) then
            call system("rm i_am_dimer*         > /dev/null 2>&1")
            call system("rm sigma_output*       > /dev/null 2>&1")
            call system("rm mask_loc_rot_atom_* > /dev/null 2>&1 ")
            write (*, *) ' setting ONETEP in PROPERTIES mode : '
            call replace_it("SINGLEPOINT", "PROPERTIES")
            call replace_it("singlepoint", "PROPERTIES")
            call flip_flag("write", true_to_false=.true.)
            call flip_flag("read", true_to_false=.false.)
            call flip_flag("ngwf_analysis", true_to_false=.true.)
            call flip_flag("do_properties", true_to_false=.true.)
            call flip_flag("popn_calculate", true_to_false=.true.)
            call flip_flag("cond_calculate", true_to_false=.true.)
            call flip_flag("polarisation_calculate", true_to_false=.true.)
         endif
      endif

      call my_dmft_exe

      if (.not. just_onetep) call system(TRIM(ADJUSTL(dir_onetep))//"/disable_write ")

   end subroutine

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

   subroutine check_flag_consistency

      if (dmft_splitkdmftall) then
         if (dmft_splitk_batch < 1) dmft_splitk_batch = 1
         dmft_split = .true.
         split_onetep = .true.
         mpi_onetep_type = 1
         all_local_host = .false.
      else
         dmft_splitk_batch = 1
      endif

      if (nproc_onetep == 1 .and. .not. dmft_splitkdmftall) split_onetep = .false.

      if (compute_dos) then
         write (*, *) 'compute dos, split onetep set to true'
         call replace_it("SINGLEPOINT", "PROPERTIES")
         call replace_it("singlepoint", "PROPERTIES")
         split_onetep = .true. ! to receive the right command line arguments in onetep dmft module (to get 5 command line arguments)
         just_onetep = .true.
         restart_from_iteration_number = 0
      endif
      if (nproc > 1 .and. nproc_mpi_solver > 1) then
         write (*, *) 'cannot use MPI for CTQMC and for ONETEP+DMFT interface (where mpi sends one job for each atoms)'
         write (*, *) 'please choose (nproc,nproc_mpi_solver) with (>1,1) or (1,>1) choices'
         stop
      endif
      if (restart_from_iteration_number > 0) then
         write (*, *) 'RESTARTING FROM ITERATION NUMBER [x] = ', restart_from_iteration_number
         write (*, *) 'switching on flag start_from_an_old_sim'
         start_from_an_old_sim = .true.
      endif
      if (all_local_host) mpi_onetep_type = 1

#ifdef debug
      hide_errors = .false.
#endif

   end subroutine

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

end module

!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

program dmftonetep
   use dmft_variables
   use openmpmod, only: init_openmp, omp_num_threads
   implicit none
   integer                    :: i, j, k, ww, iter_dmft
   character(200)             :: files, filename_sigma_source, filename_sigma
   character(8000)            :: command_line
   logical                    :: flag_onetep_producing_only_up_spin
   type(string)               :: cc_

   call init_my_input_variables

   write (*, *) 'RUNNING THE DMFT   CALCULATIONS WITH [x] cpus             : ', nproc
   write (*, *) 'RUNNING THE ONETEP CALCULATIONS WITH [x] cpus             : ', nproc_onetep
   write (*, *) 'RUNNING THE ONETEP CALCULATIONS WITH [x] threads (openmp) : ', nproc_onetep_openmp_

   if (.not. just_onetep) call system("rm -r dir_onetep_iter* > /dev/null 2>&1 ")

   !=========================================================================!
   !=========================================================================!

   iter_dmft = 0

   if (restart_from_iteration_number > 0) then
      iter_dmft = restart_from_iteration_number
      write (*, *) 'RESTARTING RUNS FROM A FORMER CALCULATIONS, ITERATION : ', restart_from_iteration_number
   endif

   do

      iter_dmft = iter_dmft + 1

      !call check_if_input_file_changed

      if (iter_dmft > niter_dmft) then
         write (*, *) 'MAX ITER REACHED, iter_dmft, niter : ', iter_dmft, niter_dmft
         exit
      endif

      if (.not. just_onetep) then
         write (*, '(a)') 'DMFT ITERATION '//trim(adjustl(toString(iter_dmft))) &
            //' of '//toString(niter_dmft)
      endif

      write (*, *) 'running onetep for case : ', trim(adjustl(CASE_ONETEP))

      if (.not. just_onetep) then
         call system(" mkdir dir_onetep_iter"//TRIM(ADJUSTL(toString(iter_dmft))))
         call system(" mkdir            before_onetep ")
         call system(" mv green_output* before_onetep > /dev/null 2>&1 ")
         call system(" cp edc_output*   before_onetep > /dev/null 2>&1 ")
         call system(" cp sigma_output* before_onetep > /dev/null 2>&1 ")
         call system(" mv before_onetep dir_onetep_iter"//TRIM(ADJUSTL(toString(iter_dmft))))
      endif

      call system(" onetep.dmft.collect.edc ")

      call onetep_normal_mode

      if (just_onetep) stop

      call system(" dmft_collect_script_dimer.out > onetep_dimer_splitting"//TRIM(ADJUSTL(toString(iter_dmft)))//" 2>&1  ")

      call system(" mkdir after_onetep ")
      call system(" cp green_output* after_onetep ")
      write (*, *) 'did we get the green functions ?'
      call system(" ls green_output* && echo 'yes there are [x] green functions : ' && ls green_output* |wc -l  ")
      call system(" ls green_output* | wc -l > list_of_greens ")
      call system(" ls green_output* >> list_of_greens ")
      call system(" order_list_of_greens.out ")
      call onetepdmft
      call system(" cp sigma_output* after_onetep ")
      call system(" mv after_onetep dir_onetep_iter"//TRIM(ADJUSTL(toString(iter_dmft))))
   enddo
   !=========================================================================!
   !=========================================================================!

contains

   !-----------------!
   !-----------------!

   subroutine check_machines_onetep_file
      implicit none
      character(200) :: myhost

877   continue

      open (unit=123, file=trim(adjustl(mach_onetep)))
      do i = 1, nproc_onetep*dmft_splitk_batch
         read (123, *, end=467) myhost
         write (*, *) 'HOST IN MACHINES_ONETEP FILE : ', trim(adjustl(myhost))
         if (.false.) then
467         continue
            write (*, *) 'ERROR file machines containing hosts does not exist or does not contain all nodes, error'
            write (*, *) 'nproc for onetep : ', nproc_onetep*dmft_splitk_batch
            write (*, *) 'nodes in file machines : '
            call system(" cat "//trim(adjustl(mach_onetep)))
            if (all_local_host) then
               close (123)
               call system(" rm "//trim(adjustl(mach_onetep)))
               call fill_onetep_machine_files_rout
               goto 877
            else
               write (*, *) 'PLEASE PROVIDE ADEQUATE '//trim(adjustl(mach_onetep))//' files with list of nodes'
               write (*, *) 'THERE SHOULD BE [X] ENTRIES : ', nproc_onetep*dmft_splitk_batch
               stop
            endif
         endif
      enddo
      close (123)

   end subroutine

   !-----------------!
   !-----------------!

   subroutine fill_onetep_machine_files_rout
      implicit none

   call system( "ls "//trim(adjustl(mach_onetep))//" || fill_machine_file "//trim(adjustl(mach_onetep))//" "//TRIM(ADJUSTL(toString(nproc_onetep*dmft_splitk_batch))) )

      if (nproc_onetep_openmp_ > 1) then !generate machines_onetep
    call system("duplicate_machine_file "//trim(adjustl(mach_onetep))//" "//TRIM(ADJUSTL(toString(nproc_onetep*dmft_splitk_batch))))
         call system("mv "//trim(adjustl(mach_onetep))//".dup "//trim(adjustl(mach_onetep)))
      endif

   end subroutine

   !-----------------!
   !-----------------!

   subroutine onetep_normal_mode

      use common_def, only: utils_system_call

      implicit none

      integer         :: jj, k, nprocess
      character(2000) :: myhost, prefix, args, output, mach_arg, mach_file
      character(3)    :: lochost
      integer         :: ierr

      call system("rm ./onetep_confirmation*  > /dev/null 2>&1 ")
      call system("rm ./green_output*         > /dev/null 2>&1 ")

      if (nproc_onetep > 1 .or. compute_dos .or. dmft_splitkdmftall) then

         call fill_onetep_machine_files_rout

         call check_machines_onetep_file

         if (.not. split_onetep) then
            if (.not. nomachinefile) then; mach_arg = ' -machinefile '; else; mach_arg = '     '; endif; 
            command_line = build_mpi_command_line(prefix=" ", np=nproc_onetep, p=1, omp=nproc_onetep_openmp_, mach=mach_onetep,&
                                           & mach_arg=mach_arg, args=CASE_ONETEP, exe=exec_onetep,&
  & outputin=trim(adjustl(case_onetep))//"_"//TRIM(ADJUSTL(toString(iter_dmft)))//".onetep", hide_errors=hide_errors, localhost="F")
            call utils_system_call(trim(adjustl(command_line)), abort=.true., &
                                   report=.true.)
         else

            !------------------------------------------------!
            !------------------------------------------------!
            if (mpi_onetep_type == 1) then
               call system("rm ./onetep_confirmation* > /dev/null 2>&1 ")

               open (unit=123, file=trim(adjustl(mach_onetep)))
               if (dmft_splitkdmftall) then
                  do i = 1, nproc_onetep
                     open (unit=691, file='splitk'//TRIM(ADJUSTL(tostring(i))))
                     do j = 1, dmft_splitk_batch
                        read (123, '(a)', end=665) myhost
                        write (*, *) 'HOST IN MACHINES_ONETEP FILE / BATCH : ', trim(adjustl(myhost)), i
                        write (691, *) trim(adjustl(myhost))
                     enddo
                     close (691)
                  enddo
               endif

               if (.false.) then
665               continue
                  write (*, *) 'ERROR not enough entries in '//trim(adjustl(mach_onetep))
                  write (*, *) 'expecting in total with [x] cpus : ', dmft_splitk_batch*nproc_onetep
                  stop
               endif
               close (123)

               open (unit=123, file=trim(adjustl(mach_onetep)))

               do i = 1, nproc_onetep
                  read (123, *, end=65) myhost
                  write (*, *) 'HOST IN MACHINES_ONETEP FILE : ', trim(adjustl(myhost))
                  if (.false.) then
65                   continue
                     write (*, *) 'ERROR file machines containing hosts does not exist or does not contain all nodes, error'
                     write (*, *) 'nproc for onetep       : ', nproc_onetep
                     write (*, *) 'batch of splik nodes   : ', dmft_splitk_batch
                     write (*, *) 'using splik nodes ?    : ', dmft_splitkdmftall
                     write (*, *) 'using openmp      ?    : ', nproc_onetep_openmp_ > 1
                     write (*, *) 'nodes in file machines : '
                     call system(" cat "//trim(adjustl(mach_onetep)))
                     stop
                  endif

                  write (*, *) ' SENDING A JOB ON HOST : ', myhost
                  write (*, *) ' ASYNC CALLS TO ONETEP : ', i

                  if (all_local_host) then
                     prefix = " export OMP_NUM_THREADS="//trim(adjustl(tostring(nproc_onetep_openmp_)))//" ; "
                     if (numa) prefix = trim(adjustl(prefix))//"  numactl  --localalloc  "
                  else
                     prefix = " "
                  endif
                  if (.not. dmft_splitkdmftall) then
                     mach_arg = " -host "
                     mach_file = TRIM(ADJUSTL(myhost))
                     nprocess = 1
                  else
                     if (.not. nomachinefile) then; mach_arg = ' -machinefile '; else; mach_arg = '     '; endif; 
                     mach_file = "splitk"//TRIM(ADJUSTL(tostring(i)))
                     nprocess = dmft_splitk_batch
                  endif

                  args = trim(adjustl(CASE_ONETEP))//" "//TRIM(ADJUSTL(toString(i)))//" "//TRIM(ADJUSTL(toString(nproc_onetep)))
                  if (compute_dos) then
                     args = trim(adjustl(args))//" "//TRIM(ADJUSTL(toString(ed_real_frequ_last)))//" 1 "//&
                            & TRIM(ADJUSTL(toString(ed_frequ_min)))//"  "//TRIM(ADJUSTL(toString(ed_frequ_max)))//"  "
                  endif
                  output = " onetep_output_iter"//TRIM(ADJUSTL(toString(iter_dmft)))//"_rank__"//TRIM(ADJUSTL(toString(i)))
                  if (all_local_host) then
                     lochost = 'T'
                  else
                     lochost = 'F'
                  endif

                  command_line = build_mpi_command_line(prefix=prefix, np=nprocess, p=1, omp=nproc_onetep_openmp_, mach=mach_file,&
                                                   & mach_arg=mach_arg, args=args, exe=exec_onetep,&
                                                   & outputin=output, hide_errors=hide_errors, localhost=lochost, ampersand=.true.)

                  write (*, *) ' command line for proc [x], total : ', i, nproc_onetep
                  call utils_system_call(command_line, abort=.true.)

               enddo

               do
                  call system("sleep 3")
                  call system("ls -l onetep_confirmation_* 2>&1 | grep -v 'No such' | wc -l > onetep_confirmation ")
                  open (unit=10, file='onetep_confirmation')
                  jj = 0
                  do k = 1, 1
                     read (10, *, end=67) jj
                  enddo
67                continue
                  close (10)
                  write (*, *) 'N PROCESS ARE DONE = ', jj
                  if (jj == nproc_onetep) then
                     write (*, *) 'ALL JOBS ARE DONE'
                     exit
                  endif
               enddo
               close (123)
               !------------------------------------------------!
               !------------------------------------------------!
            else
               open (unit=55, file='temp_onetep_part', form='unformatted')
               write (55) nproc_onetep, iter_dmft, exec_onetep, &
                      & CASE_ONETEP, compute_dos, ed_frequ_min, ed_frequ_max, ed_real_frequ_last
               close (55)

               prefix = " "
               if (.not. nomachinefile) then; mach_arg = ' -machinefile '; else; mach_arg = '     '; endif; 
               mach_file = mach_onetep
               nprocess = nproc_onetep
               args = "  "
               output = " dmft_all_iterations_output "
               command_line = build_mpi_command_line(prefix=prefix, np=nprocess, p=1, omp=nproc_onetep_openmp_, mach=mach_file,&
                                                & mach_arg=mach_arg, args=args, exe="dmft_all_iterations.onetep.out",&
                                                & outputin=output, hide_errors=.false., localhost='F')

               write (*, *) 'running onetep in mpi-2 mode : ', TRIM(ADJUSTL(command_line))
               call utils_system_call(command_line, abort=.true.)
               write (*, *) 'onetep part done'
            endif
            !------------------------------------------------!
            !------------------------------------------------!
         endif

      else

         command_line = trim(adjustl(exec_onetep))//" "//trim(adjustl(case_onetep))// &
                        " > "//trim(adjustl(case_onetep))//"_"// &
                        trim(adjustl(toString(iter_dmft)))//".onetep"

         call utils_system_call(command_line, abort=.true.)

      endif

      call system(" dmft_collect_script.out 2>&1 ")

   end subroutine

   !-----------------!
   !-----------------!

   subroutine onetepdmft

      use common_def, only: utils_system_call

      implicit none
      integer         :: nprocess
      character(2000) :: prefix, args, output, mach_arg, mach_file
      integer         :: ierr
      integer         :: omp_threads


      if (nproc > 1) then
         call system("ls machines_dmft || fill_machine_file machines_dmft "//TRIM(ADJUSTL(toString(nproc))))
         prefix = " "
         mach_arg = ' -machinefile '
         mach_file = " machines_dmft "
         nprocess = nproc
         args = "  iter_dmft="//TRIM(ADJUSTL(toString(iter_dmft)))
         output = "  onetep_dmft_part_"//TRIM(ADJUSTL(toString(iter_dmft)))
         ! ebl: make sure OMP threads get inherited by onetep_split.out
         command_line = build_mpi_command_line(prefix=prefix, np=nprocess, p=1, omp=openmp_solver, mach=mach_file,&
                                          & mach_arg=mach_arg, args=args, exe="onetep_split.out",&
                                          & outputin=output, hide_errors=.false., localhost='F')
      else
         ! ebl: make sure OMP threads get inherited by onetep_split_serial.out
         call utils_system_call("export OMP_NUM_THREADS=" // trim(adjustl(tostring(openmp_solver))), abort=.true.)
         command_line=" onetep_split_serial.out iter_dmft=" // TRIM(ADJUSTL(toString(iter_dmft))) // &
               " > onetep_dmft_part_"//TRIM(ADJUSTL(toString(iter_dmft)))
      endif

      call utils_system_call(command_line, abort=.true.)

   end subroutine

   !-----------------!
   !-----------------!

end program

!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

