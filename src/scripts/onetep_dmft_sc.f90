!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

module dmft_variables_sc
! #ifdef _plot
!  use plotlib, only : plotarray
! #endif

   use genvar, only: DP
   use namelistmod, only: namelist_set, namelist_init, putel_in_namelist, &
                                    & look_for_namelist_in_file, look_for_command_line_argument
   use strings, only: replace_in_string, string, assignment(=)
   use StringManip, only: toString
   use init_and_close_my_sim, only: initialize_my_simulation, finalize_my_simulation

integer            ::  iter_restart_sc,niter_dft_dmft_sc,niter_sc_dmft,niter_sc_dft,nproc_mpi_solver,nproc_gpu,nproc_store,nproc,niter_sc_dft_first,niter_dmft_mu,niter_kernel_mu,kerneliter,kerneliterINIT,NGWF_CGiter,diis_max
   integer            ::  spin_breaking_iter, ncpt_two_step_iter
   real(kind=DP)            ::  spin_breaking_amp
   character(2000)    ::  CASE_ONETEP, dir_onetep, exec_onetep, dir_onetep_mpi
   type(namelist_set) ::  nm
   logical            :: sandwitch_embedding, erase_chem, monitor_gpu_temperature, copy_kernel, protect_projectors,onlygammakernel
   logical            :: sc_start_from_previous_run, flag_turn_off_dmft, flag_turn_off_store
   logical            :: use_previous_dmft_files, use_same_self_energy, debug_mode_erase_sigma, fully_sc_h
   real(kind=DP)            :: hydrogenic_projectors
   logical            :: restart_from_older_dft
   integer            :: dmft_kernel_process, dmft_spin, onetep_spin, nproc_onetep_openmp
logical            :: KS_shift,purify_sc,tough_converge,kernelfix_only_first_iter,nokernelupdate,fully_sc,numa,store_sig_in_scratch,real_axis_only_last_step
   real(kind=DP)            :: mixing_dft_dmft, kernel_cutoff, cutoff_energy, mu_diff
   real(kind=DP)            :: dpmin, dpmax
   integer            :: iter_lin, nkpoints
   logical            :: matsu_solver, improve_mu_conv, save_atoms_at_each_steps
   logical            :: lin_scaling, dmft_split, dmft_splitk, no_frequ_split
   real(kind=DP)            :: lin_window
   integer            :: lin_nval, nfrequencies_dmft_dft, dmft_splitk_batch
   integer            :: nprocess_
   character(2000)    :: prefix_, args_, output_, mach_arg_, mach_file_

contains

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

   subroutine header
      write (*, *)
      write (*, *)
      call system(" echo '=========================================================='")
      call system(" echo '=========================================================='")
      call system(" echo '=========================================================='")
      call system(" echo ")
      call system(" echo ' TTTTTTT   OOO    SSSS    CCC     A    M     M   ' ")
      call system(" echo '    T     O   O  S    S  C   C   A A   MM   MM   ' ")
      call system(" echo '    T    O     O S      C        A A   MM   MM   ' ")
      call system(" echo '    T    O     O  S     C        A A   MM   MM   ' ")
      call system(" echo '    T    O     O   SS   C       A   A  M M M M   ' ")
      call system(" echo '    T    O     O     S  C       A   A  M M M M   ' ")
      call system(" echo '    T    O     O      S C      AAAAAAA M M M M   ' ")
      call system(" echo '    T     O   O  S    S  C   C A     A M  M  M   ' ")
      call system(" echo '    T      OOO    SSSS    CCC  A     A M     M   ' ")
      call system(" echo")
      call system(" echo ' a -TO-lbox of -S-trongly -C-orrelated -A-pproaches to -M-aterials '")
      call system(" echo")
      call system(" echo ' A package written by Cedric Weber - KCL - 2013 '")
      call system(" echo")
      call system(" echo '=========================================================='")
      call system(" echo '=========================================================='")
      call system(" echo '=========================================================='")
      call system(" echo ' '")
      call system(" echo ' '")
      write (*, *)
      call system("echo 'YOU STARTED YOUR JOB ON : ' `date`")
      write (*, *)
      write (*, *) 'Your input file               : input_onetep_dmft.txt '
      write (*, *) ' '
      call system("echo '=========================================================='")
      write (*, *) 'YOU WILL NEED TO DEFINE IN YOUR BASHRC FILE THE FOLLOWING ENV VARIABLES : '
      write (*, *) ' '
      write (*, *) 'ONETEP                  : the path to your onetep executble (example : /bin/onetep.linux)'
      write (*, *) 'TOSCAM                  : the path to the TOSCAM directory (example : ~/MYLIBS/TOSCAM/)'
      write (*, *) 'DMFT_ONETEP_SCRATCH_DIR : (optional) a path to a scratch directory (example : /tmp/cedric/ )'
      write (*, *) 'DMFT_ONETEP_MPI_EXEC    : your MPI prefix (example mpirun )'
      write (*, *) 'DMFT_ONETEP_BG          : set it to 1 if you are using BlueGene/IBM, 0 otherwise'
      write (*, *)
      write (*, *) 'TO RUN THE TOSCAM PACKAGE ON A PARALLEL ARCHITECTURE, YOU WILL'
      write (*, *) 'NEED TO PROVIDE AT RUN TIME THE FILES : '
      write (*, *) ' '
     write(*,*)  '  1) machines_onetep_gf_ker  : list of threads hostnames (repeated when on same node) for the GF calculations and the Kernel - nproc_gpu entries '
     write(*,*)  '  2) machines_onetep_embed   : list of threads hostnames (repeated when on same node) for the Embedding self energy calculations - nproc_gpu entries '
     write(*,*)  '  3) machines_onetep_dft     : list of threads hostnames (repeated when on same node) for the GF calculation initialization and DFT calculations - nproc_store entries '
     write(*,*)  '  4) machines_dmft           : list of threads hostnames (repeated or not when on same node, openmp will fill partially the node, see below) - nproc entries '
     write(*,*)  '  5) machines_{Solver}       : {Solver=ed,ctqmc,oca,nca} The list of threads used by the DMFT solver. ED can also use openmp, this is given by the flag openmp_solver in input_onetep_dmft.txt'
      write (*, *) '                               for the DMFT solver (each thread solves one atom) '
      write (*, *) ' '
      write (*, *) '  NB : files 1-5 are not needed for localhost calculations, they will be generated automatically'
      write (*, *) ' '
      write (*, *) ' associated variables      :  '
      write (*, *) '       -------------------- '
      write (*, *) '       nproc_gpu           : number of mpi threads for the GF calculation and Kernel calculation'
     write(*,*)  '       nproc_onetep_openmp : openmp threads associated to nproc_gpu ( total threads=nproc_gpu x nproc_onetep_openmp ) '
      write (*, *) '       nproc_store         : number of mpi threads for the GF initialization and DFT calculations'
      write (*, *) '       nproc               : number of mpi threads for the DMFT solver '
     write(*,*)  '       openmp_solver       : number of openmp threads associated to the DMFT solver (total threads= nproc x openmp_solver)'
   write (*, *) '       nproc_mpi_solver    : if nproc=1, number of mpi threads for the DMFT solver (atoms are solved sequentially)'
     write(*,*)  '       dmft_splitk_batch   : number of threads in a MPI batch (paralellization over frequ), parallelization over kpoints is '
      write (*, *) '                             done among those batches (see splitk below) for both the GF and kernel'
      write (*, *) '       --------------------                     '
      write (*, *) '       no_frequ_split      : T - the GF calculations is not parrelized along frequencies nor kpoint'
     write(*,*)  '                                 it parallelizes the matrix-vector operations with SCALAPACK if onetep compiled against it'
    write (*, *) '       dmft_split          : T - the GF calculations is parallelized along frequencies in MPI (does not use GPUs)'
     write(*,*)  '       dmft_split          : F - the GF calculations is parallelized along frequencies with NFS communications (uses GPUs if present)'
     write(*,*)  '       dmft_splitk         : T - the GF calculations is parallelized along frequencies with MPI and along kpoints with NFS (does not use GPUs)'
     write(*,*)  '       uniform_sigma       : T - the DMFT solver assumes the self energy is uniform among Hubbard atoms, it will stop if nproc=1 (to avoid waste of cpu time) '
      write (*, *) '       -------------------- '
      write (*, *) ' '
      call system("echo '=========================================================='")

   end subroutine

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!
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
      call system("rm scratch2 > /dev/null 2>&1")
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
      write (*, *) '----------------------------------------------------------'
      write (*, *) 'My TOSCAM base directory is ', TRIM(ADJUSTL(dir_onetep))
      write (*, *) '----------------------------------------------------------'
      close (10001)
      call system("rm dir_onetep > /dev/null 2>&1")

      call system("get_mpi_env")
      open (unit=10001, file='dir_onetep_mpi')
      read (10001, '(a)') dir_onetep_mpi
      write (*, *) '----------------------------------------------------------'
      write (*, *) 'My MPI executable is ', TRIM(ADJUSTL(dir_onetep_mpi))
      write (*, *) '----------------------------------------------------------'
      close (10001)
      call system("rm dir_onetep_mpi > /dev/null 2>&1")

      call system("get_onetep_exec")
      open (unit=10001, file='exec_onetep')
      read (10001, '(a)') exec_onetep
      write (*, *) '----------------------------------------------------------'
      write (*, *) 'My ONETEP executable is ', TRIM(ADJUSTL(exec_onetep))
      write (*, *) '----------------------------------------------------------'
      close (10001)
      call system("rm exec_onetep > /dev/null 2>&1")
   end subroutine

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

   subroutine check_if_input_file_changed
      implicit none
      write (*, *) 'checking for command line argument ....................'
      call look_for_namelist_in_file(nm, './input_onetep_dmft.txt')

   end subroutine

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

   subroutine init_my_input_variables
      implicit none
      integer :: numb_dat

      write (*, *) 'Setting up your list of input variables ......'

      call namelist_init(nm, 200, name_of_namelist='dmft_sc_variables')
   call putel_in_namelist(nm,improve_mu_conv,'improve_mu_conv',.true.,'Definition: if true, it will used 4th order Newton method to find mu in the onetep_dmft interface')
   call putel_in_namelist(nm,nkpoints,'nkpoints',1,'Definition: although onetep is Gamma only calculation, the summation over k-points can be done in the DMFT interface to get a better estimation of Gloc, the averaged lattice GF')
      call putel_in_namelist(nm, onetep_spin, 'onetep_spin', 1, 'Definition: spin of the ONETEP calculation (1=pm or 2=af)')
   call putel_in_namelist(nm,dmft_spin,'dmft_spin',1,'Definition: spin of the DMFT self consistence. By default doing PM calculations, is =2 will do the AF calculations, irrespectively of the ONETEP parametrization (pm or af)') 
   call putel_in_namelist(nm,copy_kernel,'copy_kernel',.true.,'Definition: if true will copy the kernel at each iteration to sc_dmft_full_iter')
      call putel_in_namelist(nm, niter_dft_dmft_sc, 'niter_dft_dmft_sc', 5, 'Definition:Number of mixed DFT+DMFT iteration')
      call putel_in_namelist(nm, nproc, 'nproc', 1, 'Definition:=number of processors for dmft')
   call putel_in_namelist(nm,nproc_mpi_solver,'nproc_mpi_solver',1,'Definition:=number of processors (MPI type) used for DMFT solver')
      call putel_in_namelist(nm, niter_sc_dmft, 'niter_sc_dmft', 1, 'Definition:Number of iteration at the dmft level')
      call putel_in_namelist(nm, niter_sc_dft, 'niter_sc_dft', 1, 'Definition:Number of iteration at the DFT level')
   call putel_in_namelist(nm,sc_start_from_previous_run,'sc_start_from_previous_run',.false.,'Definition: if true starts from a previous run and reads tightbox and dkn files')
   call putel_in_namelist(nm,nproc_gpu,'nproc_gpu',1,'Definition: number of processor for the projection of the green function part in onetep, starting from the store files for the DMFT iterations')
   call putel_in_namelist(nm, nproc_store, 'nproc_store', 8, 'Definition: number of cpu used to obtain the store files with onetep')
   call putel_in_namelist(nm,flag_turn_off_dmft,'flag_turn_off_dmft',.false.,'Definition: if true will not copy the DMFT kernel for the next ONETEP iteration')
   call putel_in_namelist(nm,flag_turn_off_store,'flag_turn_off_store',.false.,'Definition: if true does not use the store file system to speed up calculations')
   call putel_in_namelist(nm,niter_sc_dft_first,'niter_sc_dft_first',4,'Definition: number of DFT iterations for the first DFT set')
   call putel_in_namelist(nm,niter_dmft_mu,'niter_dmft_mu',1,'Definition: number of iterations used to adjust the chemical potential in DFT when preparing the input files for the DMFT')   
   call putel_in_namelist(nm,niter_kernel_mu,'niter_kernel_mu',1,'Definition: number of iterations for the Newton method to adapt the chemical potential, when the new dmft kernel is computed')
   call putel_in_namelist(nm,dmft_kernel_process,'dmft_kernel_process',1,'Definition: 1-standard dmft kernel, -1-purified dmft kernel, 2:pullay mixing') 
   call putel_in_namelist(nm,use_previous_dmft_files,'use_previous_dmft_files',.true.,'Definition: if true the code will use the previous dmft files for the next iteration')
   call putel_in_namelist(nm,use_same_self_energy,'use_same_self_energy',.false.,'Definition: if true always copy the self energy from the 1st iteration back')
   call putel_in_namelist(nm,debug_mode_erase_sigma,'debug_mode_erase_sigma',.false.,'Definition: if true erases the self energy files, for debug mode')
   call putel_in_namelist(nm,tough_converge,'tough_converge',.false.,'Definition: if true will use some flags in onetep to help fixing the density kernel, can be useful if the charge gets out of control')
   call putel_in_namelist(nm,kernelfix_only_first_iter,'kernelfix_only_first_iter',.false.,'Definition: if true will use kernel_fix:2 only for first DFT iter')
 call putel_in_namelist(nm, kerneliter, 'kerneliter', 4, 'Definition: number of iterations for kernel updates, for dft_dmft iter>2')
   call putel_in_namelist(nm,kerneliterINIT,'kerneliterINIT',10,'Definition: number of iterations for kernel updates, for dft_dmft firs titer')
   call putel_in_namelist(nm,NGWF_CGiter,'NGWF_CGiter',3,'Definition: number of iterations for NGWF update after dmft iter, and before we update the DFT kernel')
   call putel_in_namelist(nm,nokernelupdate,'nokernelupdate',.false.,'Definition: if true will not update the density kernel in the dft part')
   call putel_in_namelist(nm,fully_sc,'fully_sc',.false.,'Definition: If true will use the self energy in the energy functional in the DFT kernel/NGWF optimization')
      call putel_in_namelist(nm, mixing_dft_dmft, 'mixing_dft_dmft', 0.2d0, 'Definition: mixing of DFT AND DMFT kernel')
   call putel_in_namelist(nm,purify_sc,'purify_sc',.false.,'Definition: if true uses the purifed kernel in the DFT module properties for dft+dmft')
   call putel_in_namelist(nm,KS_shift,'KS_shift',.true.,'Definition: if true adds the correction induced by shift of occupations in the Kernel DFT minimization where the DFT energy is introduced')
   call putel_in_namelist(nm,fully_sc_h,'fully_sc_h',.false.,'Definition: if true will use the one shot Kohn Sham hamiltonian obtained by DMFT in onetep')
      call putel_in_namelist(nm, diis_max, 'diis_max', 5, 'Definition: number of kernel iterations kept for pullay mixing')
   call putel_in_namelist(nm,sandwitch_embedding,'sandwitch_embedding',.false.,'Definition: if true it will generate both left and right embedding files, used for building sandwitched hetero-structures')
   call putel_in_namelist(nm,erase_chem,'erase_chem',.false.,'Definition: if true erases the chemical potential between different DMFT iteration')
     call putel_in_namelist(nm, numa, 'numa', .false., 'Definition: numa avoids remote memory access in multi-socket architectures')
   call putel_in_namelist(nm,monitor_gpu_temperature,'monitor_gpu_temperature',.false.,'Definition: if true will prepare a graphic of the GPU temperature during the ONETEP+DMFT run, critical temperature is about 75 Celsius degree')
      call putel_in_namelist(nm, kernel_cutoff, 'kernel_cutoff', 4000.0d0, 'Definition:cutoff for kernel truncation')
      call putel_in_namelist(nm, cutoff_energy, 'cutoff_energy', 850.0d0, 'Definition:energy plane wave cutoff')
      call putel_in_namelist(nm, iter_lin, 'iter_lin', 3, 'Definition: from 1-iter_lin the mixing is linear')
   call putel_in_namelist(nm,mu_diff,'mu_diff',0.04d0,'Definition:precision (max deviation) from target density in the DMFT when the chemical potential is obtained by the Newton-Parston algorithm')
   call putel_in_namelist(nm,matsu_solver,'matsu_solver',.false.,'Definition: if true does not generate real axis quantities (matsubara solver)')
      call putel_in_namelist(nm, dpmin, 'dpmin', 0.05d0, 'Definition:')
      call putel_in_namelist(nm, dpmax, 'dpmax', 0.60d0, 'Definition:')
   call putel_in_namelist(nm,iter_restart_sc,'iter_restart_sc',1,'Definition: restart from this iter the calculation, this should be the iter number which is NOT yet done, iter_restart_sc-1 should be available')
   call putel_in_namelist(nm,restart_from_older_dft,'restart_from_older_dft',.false.,'Definition: if true it will run a usual onetep+dmft calculations, and ALSO carry out DFT iterations at the beginning, but it will restart from a previous DFT run')
   call putel_in_namelist(nm,protect_projectors,'protect_projectors',.true.,'Definition: will make the file run.tightbox_hubbard read only, such that it cannot be rewritten during a self consistent simulation, this is for debug mode only and should not affect a simulation')
   call putel_in_namelist(nm,lin_scaling,'lin_scaling',.false.,'Definition: if true the code will use for the DMFT part (not the kernel) a linear scaling algorithm to invert the GF')
   call putel_in_namelist(nm,lin_window,'lin_window',0.4d0,'Definition: window of energy around fermi level for linear scaling inversion')
      call putel_in_namelist(nm, lin_nval, 'lin_nval', 40, 'Definition: number of eigenvalues for linear inversion of GF')
   call putel_in_namelist(nm,dmft_split,'dmft_split',.false.,'Definition: if true will split the onetep dmft interface over cpus with MPI rather than NFS')
   call putel_in_namelist(nm,dmft_splitk,'dmft_splitk',.false.,'Definition: if true splits the mpi onetep dmft interface in several batches each of them running different K points')
   call putel_in_namelist(nm,store_sig_in_scratch,'store_sig_in_scratch',.false.,'Definition: if true the sigma_output files will be copied to local scratch when running onetep, it will aleviate NFS communications')
   call putel_in_namelist(nm,save_atoms_at_each_steps,'save_atoms_at_each_steps',.true.,'Definition: if true TOSCAM will save the atoms directories at every self consistent iteration, this can generate a large amount of data')
   call putel_in_namelist(nm,nproc_onetep_openmp,'nproc_onetep_openmp',1,'Definition: will run mpi for the onetep dmft interface, but each mpi job will run with nproc_onetep_openmp openmp threads')
   call putel_in_namelist(nm,no_frequ_split,'no_frequ_split',.false.,'Definition: if true the onetep dmft interface will not split over frequencies, but use the sparse parallelization implemented within onetep (parallelization of matrix-vector operations, if onetep compiled with SCALAPACK flag, it will use it')
   call putel_in_namelist(nm,nfrequencies_dmft_dft,'nfrequencies_dmft_dft',160,'Definition: number of frequencies for the DFT+DMFT calculations')
   call putel_in_namelist(nm,dmft_splitk_batch,'dmft_splitk_batch',1,'Definition: number of cpus in each of the batch when splitting k points')
   call putel_in_namelist(nm,onlygammakernel,'onlygammakernel',.false.,'Definition: if true will not estimate the DFT kernel from the k-point average obtained from DMFT, but directly from G=0, however the chemical potential will still be fixed by using the k-point averaged GF')
   call putel_in_namelist(nm,real_axis_only_last_step,'real_axis_only_last_step',.false.,'Definition: if true the real axis quantities obtained by the ED solver are only obtained for the last two DMFT steps')
   call putel_in_namelist(nm,spin_breaking_amp,'spin_breaking_amp',0.05d0,'Definition: amplitude of the initial symmetry breaking for AF state')
   call putel_in_namelist(nm,spin_breaking_iter,'spin_breaking_iter',1,'Definition: number of iterations where the symmetry breaking is applied on the impurity energy levels') 
   call putel_in_namelist(nm,hydrogenic_projectors,'hydrogenic_projectors',-1.d0,'Definition: if >0 will use hydrogenic projectors with given effective charge instead of atomic orbitals used by the solver')
   call putel_in_namelist(nm,ncpt_two_step_iter,'ncpt_two_step_iter',100000,'Definition: if iter_dmft_sc < ncpt_two_step_iter then the CPT ED solver will use a two step procedure to fit the hybridization, where in the first step the hybrid. is fitted with the standard ED method, and in the second stage the CPT fit is turned on')

      !-------------------------------!

      call look_for_namelist_in_file(nm, './input_onetep_dmft.txt')
      call look_for_command_line_argument(nm)

      call system("get_onetep_case")
      call check_flag_consistency

      write (*, *) 'Variables defined, now proceed with consistency checks '

      open (unit=10002, file='case_onetep')
      read (10002, *) CASE_ONETEP
      close (10002)

      write (*, *) '----------------------------------------------------------'
      write (*, *) 'my ONETEP CASE is : ', TRIM(ADJUSTL(CASE_ONETEP))
      write (*, *) '----------------------------------------------------------'

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

      call my_dmft_exe

      if (monitor_gpu_temperature) call system("monitor_temp.out &")

   end subroutine

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

   subroutine check_flag_consistency
      if (nproc > 1 .and. nproc_mpi_solver > 1) then
         write (*, *) 'cannot use MPI for CTQMC and for ONETEP+DMFT interface (where mpi sends one job for each atoms)'
         write (*, *) 'please choose (nproc,nproc_mpi_solver) with (>1,1) or (1,>1) choices'
         stop
      endif
      if (iter_restart_sc > 1) then
         sc_start_from_previous_run = .true.
      endif
      if (onetep_spin == 2) dmft_spin = 2
      if (dmft_spin < 1 .or. dmft_spin > 2) then
         write (*, *) 'ERROR dmft_spin should be 1 or 2'
         stop
      endif
      if (onetep_spin < 1 .or. onetep_spin > 2) then
         write (*, *) 'ERROR onetep_spin should be 1 or 2'
         stop
      endif
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
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

program dmftonetep
   use dmft_variables_sc
   use common_def, only: utils_system_call
   implicit none
   integer                    :: i, j, k, ww, iter_dmft_sc, iji
   character(200)             :: files, filename_sigma_source, filename_sigma
   character(8000)            :: command_line
   logical                    :: flag_onetep_producing_only_up_spin, check_mu
   type(string)               :: cc_
   real(kind=DP), allocatable        :: chem_pots(:, :)
   logical                    :: spoil, diis

   call header

   spoil = .false.; diis = .false.

   call init_my_input_variables

   if (fully_sc .and. spoil) then
      write (*, *) 'FULLY_SC and SPOIL together, remove spoil'
      stop
   endif

   write (*, *) '----------------------------------------------------------'
   write (*, *) 'If necessary, create your scratch directory : '
   call system(" get_scratch_dir ")
   write (*, *) '----------------------------------------------------------'

   call system(" rm updating_variables > /dev/null 2>&1 ")

   write (*, *) '----------------------------------------------------------'
   write (*, *) 'add missing variables to your onetep case input file ....'
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'spin_polarized '     > /dev/null 2>&1 || echo 'spin_polarized  : F'        >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_paramagnetic '  > /dev/null 2>&1 || echo 'dmft_paramagnetic : T'      >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_cutoff_small '  > /dev/null 2>&1 || echo 'dmft_cutoff_small  : 0.01'  >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_cutoff_tail  '  > /dev/null 2>&1 || echo 'dmft_cutoff_tail  : 0.60'   >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'kerfix '             > /dev/null 2>&1 || echo 'kerfix : 2'                 >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'maxit_pen '          > /dev/null 2>&1 || echo 'maxitpen : 5'               >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'maxit_kernel_fix '   > /dev/null 2>&1 || echo 'maxit_kernel_fix : 5'       >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_spoil_kernel '  > /dev/null 2>&1 || echo 'dmft_spoil_kernel : F'      >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'kernel_diis_maxit '  > /dev/null 2>&1 || echo 'kernel_diis_maxit : 1'      >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'kernel_diis_max '    > /dev/null 2>&1 || echo 'kernel_diis_max : 5'        >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'kernel_diis '        > /dev/null 2>&1 || echo 'kernel_diis : F'            >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'maxit_ngwf_diis '    > /dev/null 2>&1 || echo 'maxit_ngwf_diis : 0'        >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_fully_sc '      > /dev/null 2>&1 || echo 'dmft_fully_sc : F'          >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_kernel_mix '    > /dev/null 2>&1 || echo 'dmft_kernel_mix : 0.1'      >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_purify_sc '     > /dev/null 2>&1 || echo 'dmft_purify_sc : F '        >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_KS_shift '      > /dev/null 2>&1 || echo 'dmft_KS_shift : F '         >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_fully_sc_h '    > /dev/null 2>&1 || echo 'dmft_fully_sc_h : F '       >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_mu_diff_max '   > /dev/null 2>&1 || echo 'dmft_mu_diff_max : 0.001 '  >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'kernel_cutoff '      > /dev/null 2>&1 || echo 'kernel_cutoff : 4000.0 '    >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
   if (cutoff_energy > 0.0) &
& call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'cutoff_energy '      > /dev/null 2>&1 || echo 'cutoff_energy : 850.0 eV'   >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'kernel_diis_liter '  > /dev/null 2>&1 || echo 'kernel_diis_liter : 3 '     >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_mu_order '      > /dev/null 2>&1 || echo 'dmft_mu_order : 4 '         >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_nkpoints '      > /dev/null 2>&1 || echo 'dmft_nkpoints : 1 '         >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_skip_energy '   > /dev/null 2>&1 || echo 'dmft_skip_energy : T '      >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_kpoints_sym '   > /dev/null 2>&1 || echo 'dmft_kpoints_sym : T '      >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_lin_scaling '   > /dev/null 2>&1 || echo 'dmft_lin_scaling : F '      >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_nval '          > /dev/null 2>&1 || echo 'dmft_nval : 40 '            >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_win '           > /dev/null 2>&1 || echo 'dmft_win : 0.1 '            >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_split '         > /dev/null 2>&1 || echo 'dmft_split : F '            >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_splitk '        > /dev/null 2>&1 || echo 'dmft_splitk : F '           >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_local_scratch ' > /dev/null 2>&1 || echo 'dmft_local_scratch : F '    >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_scaling_nmpi '  > /dev/null 2>&1 || echo 'dmft_scaling_nmpi : 16 '    >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_scaling_meth '  > /dev/null 2>&1 || echo 'dmft_scaling_meth : 5 '     >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_scaling_tail '  > /dev/null 2>&1 || echo 'dmft_scaling_tail : 0.3 '   >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_scaling_cutoff '> /dev/null 2>&1 || echo 'dmft_scaling_cutoff : 0.1 ' >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
 call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_kpoints_kernel_gamma '> /dev/null 2>&1 || echo 'dmft_kpoints_kernel_gamma : F ' >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")

   write (*, *) '----------------------------------------------------------'

   write (*, *) 'update your onetep case input file ....'

   if (onlygammakernel) Then
      call system(" update_case_file 'dmft_kpoints_kernel_gamma ' T ")
   else
      call system(" update_case_file 'dmft_kpoints_kernel_gamma ' F ")
   endif
   call system(" update_case_file 'dmft_points ' "//trim(adjustl(tostring(nfrequencies_dmft_dft))))

   call system(" update_case_file 'dmft_lin_scaling ' F ")
   call system(" update_case_file 'dmft_split '  F ")
   call system(" update_case_file 'dmft_splitk ' F ")
   call system("echo 'dmft_splitkdmftall=.false.' >> input_onetep_dmft.txt ")

   if (store_sig_in_scratch) then
      call system(" update_case_file 'dmft_local_scratch ' T ")
   else
      call system(" update_case_file 'dmft_local_scratch ' F ")
   endif

   call system(" update_case_file 'dmft_win ' "//trim(adjustl(tostring(lin_window))))
   call system(" update_case_file 'dmft_nval ' "//trim(adjustl(tostring(lin_nval))))
   call system(" update_case_file 'dmft_scaling_cutoff ' 0.001")
   call system(" update_case_file 'dmft_scaling_tail '    0.4")

   call system(" update_case_file 'dmft_skip_energy ' T")
   call system(" update_case_file 'kernel_cutoff ' "//trim(adjustl(tostring(kernel_cutoff))))
   if (cutoff_energy > 0.0) &
 & call system(" update_case_file 'cutoff_energy ' "//trim(adjustl(tostring(cutoff_energy)))//' eV')
   call system(" update_case_file 'kernel_diis_liter ' "//trim(adjustl(tostring(iter_lin))))
   call system(" update_case_file 'dmft_mu_diff_max ' "//trim(adjustl(tostring(mu_diff))))

   call system(" update_case_file 'dmft_cutoff_tail ' "//trim(adjustl(tostring(dpmax))))
   call system(" update_case_file 'dmft_cutoff_small ' "//trim(adjustl(tostring(dpmin))))

   call system(" update_case_file 'kernel_diis_max ' 0 ")
   if (KS_shift) then
      call system(" update_case_file 'dmft_KS_shift ' T ")
   else
      call system(" update_case_file 'dmft_KS_shift ' F ")
   endif
   if (purify_sc) then
      call system(" update_case_file 'dmft_purify_sc ' T ")
   else
      call system(" update_case_file 'dmft_purify_sc ' F ")
   endif
   call system(" update_case_file 'dmft_kernel_mix ' "//trim(adjustl(tostring(mixing_dft_dmft))))
   call system(" update_case_file 'kerfix ' 2 ")

   if (improve_mu_conv) then
      call system(" update_case_file 'dmft_mu_order ' 4 ")
   else
      call system(" update_case_file 'dmft_mu_order ' 2 ")
   endif

   call system(" update_case_file 'dmft_nkpoints ' "//trim(adjustl(tostring(nkpoints))))

   if (tough_converge) then
      call system(" update_case_file 'maxit_pen '        7 ")
      call system(" update_case_file 'maxit_kernel_fix ' 7 ")
   else
      call system(" update_case_file 'maxit_pen '        0 ")
      call system(" update_case_file 'maxit_kernel_fix ' 0 ")
   endif

   if (sandwitch_embedding) then
      write (*, *) '----------------------------------------------------------'
      write (*, *) 'preparing embedding scheme ....'
      call system("mkdir embedding_dir")
      call system(" mv connections_D_L* connections_D_R* embedding_dir ")
      write (*, *) '----------------------------------------------------------'
   endif

   write (*, *) '----------------------------------------------------------'
   write (*, *) 'checking your frequency profiles - if present (mask_frequ files) '
   write (*, *) 'Note: mask_frequ files are optional, but they are important      '
   write (*, *) '      to balance the load properly ona mixed GPU/CPU architecture'
   write (*, *) ' '
   call system(" check_mask_frequ ")
   write (*, *)
   write (*, *) '----------------------------------------------------------'
   write (*, *) ' about to start the DMFT loop .... '
   write (*, *) 'RUNNING THE DMFT   CALCULATIONS             WITH [x] cpus         : ', nproc
   write (*, *) 'RUNNING THE ONETEP CALCULATIONS STORE FILES WITH [x] cpus         : ', nproc_store
   write (*, *) 'RUNNING THE ONETEP CALCULATIONS GREEN FUNC  WITH [x] cpus         : ', nproc_gpu
   write (*, *) 'RUNNING THE DMFT   CALCULATIONS GREEN FUNC  WITH [x] openmp proc. : ', nproc_onetep_openmp
   write (*, *) 'RUNNING THE DMFT   CALCULATIONS GREEN FUNC  WITH [x] gpus         : '
   call system(" echo  `cat gpu_max  2>/dev/null || echo '0' ` ")
   if (mod(nfrequencies_dmft_dft, nproc_gpu) /= 0 .and. dmft_split) then
      write (*, *) 'number of frequencies not commensurate with the number of CPUs'
      write (*, *) 'number of CPUs associated with dmft_split (nproc_gpu) : ', nproc_gpu
      write (*, *) 'number of frequencies is (nfrequencies_dmft_dft)      : ', nfrequencies_dmft_dft
      stop
   endif
   if (mod(nfrequencies_dmft_dft, dmft_splitk_batch) /= 0 .and. dmft_splitk) then
      write (*, *) 'number of frequencies not commensurate with the number of CPUs in dmft_splitk'
      write (*, *) 'number of CPUs associated with dmft_splitk (dmft_splitk_batch) : ', nproc_gpu
      write (*, *) 'number of frequencies is (nfrequencies_dmft_dft)               : ', nfrequencies_dmft_dft
      stop
   endif
   write (*, *) '----------------------------------------------------------'
   write (*, *) '----------------------------------------------------------'
   write (*, *) 'purging AGR output graphic directory'
   call system(" rm -r AGR > /dev/null 2>&1 ")
   call system(" mkdir AGR > /dev/null 2>&1 ")
   write (*, *) '----------------------------------------------------------'

   !=========================================================================!
   !=========================================================================!
   !=========================================================================!
   !=========================================================================!

   write (*, *) '----------------------------------------------------------'
   write (*, *) 'making TOSCAM interface aware of the DFT calculations '
   call system("echo 'start_from_an_old_sim=.true.'  >> input_onetep_dmft.txt ")
   call system("echo 'last_iter_is_real=.false.'     >> input_onetep_dmft.txt ")
   if (.not. real_axis_only_last_step) then
      call system("echo 'ed_no_real_overide=.false.'   >> input_onetep_dmft.txt ")
   else
      call system("echo 'ed_no_real_overide=.true.'    >> input_onetep_dmft.txt ")
   endif
   write (*, *) '----------------------------------------------------------'

   write (*, *) '----------------------------------------------------------'
   write (*, *) 'did we restart from an older run ? '
   iter_dmft_sc = iter_restart_sc - 1

   if (iter_dmft_sc < 0) then
      write (*, *) 'ERROR : first iteration number should be positive'
      stop
   endif
   if (iter_dmft_sc > 0) then
      write (*, *) 'WARNING : starting back from iteration (will then go to the next iter)', iter_dmft_sc
   write(*,*) 'WARNING : you should copy run.dkn atom* sigma_output* and chem_potential* from the directory sc_dmft_full_iter'//trim(adjustl(toString(iter_dmft_sc)))
   endif
   write (*, *) '----------------------------------------------------------'
   write (*, *) 'Setting parameters for chemical potential convergence'
   call system(" echo niter_dmft="//trim(adjustl(toString(niter_sc_dmft)))//" >> input_onetep_dmft.txt  ")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_sc '                  > /dev/null 2>&1|| echo 'dmft_sc     : F' >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_kernel '              > /dev/null 2>&1|| echo 'dmft_kernel : 0' >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'dmft_nmu_loop '            > /dev/null 2>&1|| echo 'dmft_nmu_loop : 1' >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'write_converged_dkngwfs '  > /dev/null 2>&1|| echo 'write_converged_dkngwfs     : T' >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")
  call system(" cat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat | grep 'kernel_update '            > /dev/null 2>&1|| echo 'kernel_update : T' >> "//TRIM(ADJUSTL(CASE_ONETEP))//".dat")

   call system(" update_case_file 'write_converged_dkngwfs ' T ")
   call system(" update_case_file 'dmft_nmu_loop ' "//trim(adjustl(toString(niter_dmft_mu))))
   call system(" update_case_file 'dmft_kernel ' "//trim(adjustl(toString(dmft_kernel_process))))
   call system(" update_case_file 'kernel_update '          T ")

   allocate (chem_pots(niter_dft_dmft_sc, 2))

   write (*, *) '----------------------------------------------------------'
   write (*, *) 'cleaning files which could have been obtained from a previous run'
   call system(" rm dens_kern_tmp* > /dev/null 2>&1 ")
   call system(" rm store_*        > /dev/null 2>&1 ")
   call system(" rm fort.*         > /dev/null 2>&1 ")
   call system(" rm splitk*        > /dev/null 2>&1 ")
   write (*, *) '----------------------------------------------------------'

   if (fully_sc_h .and. .not. fully_sc) then
      write (*, *) 'fully_sc_h requires fully_sc to be true as well'
      stop
   endif

   write (*, *) '----------------------------------------------------------'
   write (*, *) '-------------------STARTING DMFT--------------------------'
   write (*, *) '----------------------------------------------------------'

   do

      !-------------------------!
      !-------------------------!
      !-------------------------!
      !-------------------------!

      iter_dmft_sc = iter_dmft_sc + 1

      write (*, *) '----------------------------------------------------------'
      write (*, *) 'DOING DFT+DMFT ITERATION NUMBER : ', iter_dmft_sc
      write (*, *) 'OUT OF                          : ', niter_dft_dmft_sc
      write (*, *) '----------------------------------------------------------'

      if (iter_dmft_sc > niter_dft_dmft_sc) then
         write (*, *) '----------------------------------------------------------'
         write (*, *) 'MAX ITER REACHED, iter_dmft_sc, niter : ', iter_dmft_sc, niter_dft_dmft_sc
         write (*, *) '----------------------------------------------------------'
         exit
      endif
      write (*, *) 'RUNNING DMFT ITERATION [x] : ', iter_dmft_sc, niter_dft_dmft_sc
      write (*, *) 'running onetep for case    : ', trim(adjustl(CASE_ONETEP))

      !-------------------------!
      !-------------------------!
      !-------------------------!
      !-------------------------!

      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      write (*, *) 'STEP1 : ONETEP ITER'
      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

      call system(" update_case_file 'task ' SINGLEPOINT ")
      call system(" update_case_file 'ngwf_analysis '    F ")
      call system(" update_case_file 'do_properties '    F ")

      if (iter_dmft_sc > 1 .or. sc_start_from_previous_run) then
         if (protect_projectors) then
            call system(" ls *.tightbox_hub_projs  > /dev/null 2>&1 && chmod au-w *.tightbox_hub_projs ")
         endif
         if (fully_sc) then
            call system(" update_case_file 'dmft_fully_sc ' T ")
         else
            call system(" update_case_file 'dmft_fully_sc ' F ")
         endif
         if (fully_sc_h) then
            call system(" update_case_file 'dmft_fully_sc_h ' T ")
         else
            call system(" update_case_file 'dmft_fully_sc_h ' F ")
         endif
         if (spoil) then
            call system(" update_case_file 'dmft_spoil_kernel ' T ")
         else
            call system(" update_case_file 'dmft_spoil_kernel ' F ")
         endif
         if (diis) then
            call system(" update_case_file 'kernel_diis ' T ")
            call system(" update_case_file 'kernel_diis_maxit ' 1 ")
            call system(" update_case_file 'maxit_ngwf_diis '  0 ")
         else
            call system(" update_case_file 'kernel_diis ' F ")
            call system(" update_case_file 'kernel_diis_maxit ' 0 ")
            call system(" update_case_file 'maxit_ngwf_diis '  0 ")
         endif
         call system(" update_case_file 'minit_lnv ' "//trim(adjustl(toString(kerneliter))))
         call system(" update_case_file 'maxit_lnv ' "//trim(adjustl(toString(kerneliter))))
         call system(" update_case_file 'read_denskern '         T ")
         call system(" update_case_file 'read_tightbox_ngwfs '   T ")
         call system(" update_case_file 'write_denskern '        T ")
         call system(" update_case_file 'write_tightbox_ngwfs '  T ")
         call system(" update_case_file 'write_forces '          T ")
         if (hydrogenic_projectors < 0.0) then
            call system(" update_case_file 'hubbard_proj_mixing ' -2.0 ")
            call system(" update_case_file_flip_hub +10             ")
         else
            call system(" update_case_file 'hubbard_proj_mixing ' 0.0 ")
            call system(" update_case_file_flip_hub "//trim(adjustl(tostring(hydrogenic_projectors))))
         endif
         call system(" update_case_file 'maxit_ngwf_cg ' "//trim(adjustl(toString(niter_sc_dft))))
         if (niter_sc_dft == 0) Then
            call system(" update_case_file 'write_tightbox_ngwfs '  F ")
         endif
         if (kernelfix_only_first_iter) call system(" update_case_file 'kerfix ' 0 ")
         if (nokernelupdate) call system(" update_case_file 'kernel_update ' F ")
         if (onetep_spin == 1) then
            if (dmft_spin == 1) Then
               call system(" update_case_file 'spin_polarized '    F ")
               call system(" update_case_file 'dmft_paramagnetic ' F ")
            else
               call system(" update_case_file 'spin_polarized '    T ")
               call system(" update_case_file 'dmft_paramagnetic ' F ")
            endif
         endif
         if (onetep_spin == 2) then
            if (dmft_spin == 1) Then
               call system(" update_case_file 'spin_polarized '    T ")
               call system(" update_case_file 'dmft_paramagnetic ' T ")
            else
               call system(" update_case_file 'spin_polarized '    T ")
               call system(" update_case_file 'dmft_paramagnetic ' F ")
            endif
         endif
      else
         call system(" update_case_file 'dmft_fully_sc_h ' F ")
         call system(" update_case_file 'dmft_fully_sc ' F ")
         call system(" update_case_file 'kernel_diis '    F ")
         call system(" update_case_file 'kernel_diis_maxit ' 0 ")
         call system(" update_case_file 'maxit_ngwf_diis '  0 ")
         call system(" update_case_file 'dmft_spoil_kernel ' F ")
         call system(" update_case_file 'minit_lnv ' "//trim(adjustl(toString(kerneliterINIT))))
         call system(" update_case_file 'maxit_lnv ' "//trim(adjustl(toString(kerneliterINIT))))
         if (.not. restart_from_older_dft) then
            call system(" update_case_file 'read_denskern '         F ")
            call system(" update_case_file 'read_tightbox_ngwfs '   F ")
         else
            call system(" update_case_file 'read_denskern '         T ")
            call system(" update_case_file 'read_tightbox_ngwfs '   T ")
         endif
         call system(" update_case_file 'write_denskern '        T ")
         call system(" update_case_file 'write_tightbox_ngwfs '   T ")
         call system(" update_case_file 'write_forces '           T ")

         if (hydrogenic_projectors < 0.0) then
            call system(" update_case_file 'hubbard_proj_mixing ' 0.0  ")
            call system(" update_case_file_flip_hub -10             ")
         else
            call system(" update_case_file 'hubbard_proj_mixing ' 0.0 ")
            call system(" update_case_file_flip_hub "//trim(adjustl(tostring(hydrogenic_projectors))))
         endif

         call system(" update_case_file 'maxit_ngwf_cg ' "//trim(adjustl(toString(niter_sc_dft_first))))
         if (onetep_spin == 1) then
            if (dmft_spin == 1) Then
               call system(" update_case_file 'spin_polarized '    F ")
               call system(" update_case_file 'dmft_paramagnetic ' F ")
            else
               call system(" update_case_file 'spin_polarized '    F ")
               call system(" update_case_file 'dmft_paramagnetic ' F ")
            endif
         endif
         if (onetep_spin == 2) then
            if (dmft_spin == 1) Then
               call system(" update_case_file 'spin_polarized '    T ")
               call system(" update_case_file 'dmft_paramagnetic ' T ")
            else
               call system(" update_case_file 'spin_polarized '    T ")
               call system(" update_case_file 'dmft_paramagnetic ' F ")
            endif
         endif
      endif

      if (erase_chem) call system(" rm chem.potential.nmu.iter* > /dev/null 2>&1 ")

      if ((iter_dmft_sc == 1 .and. niter_sc_dft_first > 0 .and. .not. sc_start_from_previous_run) &
         & .or. (iter_dmft_sc == 1 .and. niter_sc_dft > 0 .and. sc_start_from_previous_run) &
         & .or. (iter_dmft_sc /= 1 .and. niter_sc_dft > 0)) then

         if (numa) then
            prefix_ = " numactl  --localalloc "
         else
            prefix_ = " "
         endif
         mach_arg_ = ' -machinefile '
         mach_file_ = " machines_onetep_dft "
         nprocess_ = nproc_store
         args_ = trim(adjustl(CASE_ONETEP))
         output_ = "   sc_onetep_output_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc)))
         command_line = build_mpi_command_line(prefix=prefix_, np=nprocess_, p=1, omp=1, mach=mach_file_,&
                                          & mach_arg=mach_arg_, args=args_, exe=exec_onetep,&
                                          & outputin=output_, hide_errors=.false., localhost='F')

         call system("ls machines_onetep_dft || fill_machine_file machines_onetep_dft "//TRIM(ADJUSTL(toString(nproc_store))))
         call utils_system_call(command_line, abort=.true.)

         if (iter_dmft_sc == 1) then
            call system("mkdir backup_dft_run")
            call system("cp sc_onetep_output_iter* *.dat *.txt *.tightbox* *.dkn* backup_dft_run/ ")
         endif
      endif

      if (flag_turn_off_store) call system(" rm store_* > /dev/null 2>&1 ")

      !-------------------------!
      !-------------------------!
      !-------------------------!
      !-------------------------!

      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      write (*, *) 'STEP2 : cp run.dat run.dat.backup'
      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

      call system(" cp "//TRIM(ADJUSTL(CASE_ONETEP))//".dat "//TRIM(ADJUSTL(CASE_ONETEP))//".dat.backup ")
      call system(" cp input_onetep_dmft.txt input_onetep_dmft.txt.backup ")

      !-------------------------!
      !-------------------------!
      !-------------------------!
      !-------------------------!

      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      write (*, *) 'STEP3 : STORE FILES'
      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

      if (.false.) Then
         call system(" update_case_file 'dmft_fully_sc '    F ")
         call system(" update_case_file 'dmft_fully_sc_h '  F ")
      else
         !here if those are set to TRUE, it will use H_oneshotdmft instead of H for
         !the DMFT iteration
         !.....!
      endif
      !SET THE SPIN CONDITIONS
      if (dmft_spin == 1) then
         call system("echo 'paramagnetic=1'        >> input_onetep_dmft.txt ")
      else
         call system("echo 'paramagnetic=2'        >> input_onetep_dmft.txt ")
      endif
      if ((dmft_spin == 1 .and. onetep_spin == 1) .or. &
           & (dmft_spin == 1 .and. onetep_spin == 2)) then
         !in the second case, dmft_paramagnetic=T should be present in run.dat to only get out the UP files
         call system("echo 'onetep_only_up=.true.'  >> input_onetep_dmft.txt ")
      else
         call system("echo 'onetep_only_up=.false.' >> input_onetep_dmft.txt ")
      endif
      if (dmft_spin == 2 .and. iter_dmft_sc <= spin_breaking_iter) then
         call system("echo 'amp_slight_sym_breaking='"//TRIM(ADJUSTL(toString(spin_breaking_amp)))//" >> input_onetep_dmft.txt ")
      else
         call system("echo 'amp_slight_sym_breaking=0.00' >> input_onetep_dmft.txt ")
      endif
      if (abs(nkpoints) > 1) then
         call system("echo 'sum_over_k_dft=.true.'  >> input_onetep_dmft.txt ")
      else
         call system("echo 'sum_over_k_dft=.false.' >> input_onetep_dmft.txt ")
      endif
      !END SET

      call system("echo 'nproc_onetep_openmp_=1' >> input_onetep_dmft.txt ")

      call system(" update_case_file 'dmft_sc '                T ")

      call system(" update_case_file 'ngwf_analysis '          F ")
      call system(" update_case_file 'do_properties '          F ")
      call system(" update_case_file 'popn_calculate '         F ")
      call system(" update_case_file 'cond_calculate '         F ")
      call system(" update_case_file 'polarisation_calculate ' F ")

      call system(" update_case_file task PROPERTIES ")

      if (hydrogenic_projectors < 0.0) then
         call system(" update_case_file 'hubbard_proj_mixing ' -2.0")
         call system(" update_case_file_flip_hub +10")
      else
         call system(" update_case_file 'hubbard_proj_mixing ' 0.0 ")
         call system(" update_case_file_flip_hub "//trim(adjustl(tostring(hydrogenic_projectors))))
      endif

      call system(" update_case_file 'read_denskern '         T ")
      call system(" update_case_file 'read_tightbox_ngwfs '   T ")
      call system(" update_case_file 'write_denskern '        F ")
      call system(" update_case_file 'write_tightbox_ngwfs '  F ")
      call system(" update_case_file 'write_forces '          F ")

      call system(" update_case_file 'kernel_diis '     F ")
      call system(" update_case_file 'kernel_diis_maxit ' 0 ")
      call system(" update_case_file 'maxit_ngwf_diis '  0 ")
      call system(" update_case_file 'dmft_spoil_kernel ' F ")

      if (.not. flag_turn_off_store) then
         call system("echo 'just_onetep=.true.'   >> input_onetep_dmft.txt ")
         call system("echo 'split_onetep=.false.' >> input_onetep_dmft.txt ")
         call system("echo 'nproc_onetep='"//TRIM(ADJUSTL(toString(nproc_store)))//" >> input_onetep_dmft.txt ")
         call system("echo mach_onetep='machines_onetep_dft' >> input_onetep_dmft.txt ")

         call system("onetep.dmft > sc_onetep_store_files_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc)))//" 2>&1 ")
         call system("mkdir sc_store_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         call system("mv onetep_* sc_store_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         if (hydrogenic_projectors < 0.0) then
            call system(" update_case_file 'hubbard_proj_mixing ' 1.0")
            call system(" update_case_file_flip_hub +10")
         else
            call system(" update_case_file 'hubbard_proj_mixing ' 0.0 ")
            call system(" update_case_file_flip_hub "//trim(adjustl(tostring(hydrogenic_projectors))))
         endif
      else
         call system(" rm store_* > /dev/null 2>&1 ")
         if (hydrogenic_projectors < 0.0) then
            call system(" update_case_file 'hubbard_proj_mixing ' -2.0")
            call system(" update_case_file_flip_hub +10")
         else
            call system(" update_case_file 'hubbard_proj_mixing ' 0.0 ")
            call system(" update_case_file_flip_hub "//trim(adjustl(tostring(hydrogenic_projectors))))
         endif
      endif

      !-------------------------!
      !-------------------------!
      !-------------------------!
      !-------------------------!

      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      write (*, *) 'STEP4: DMFT ITERATION'
      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

      if (iter_dmft_sc <= ncpt_two_step_iter) then
         call system("echo 'ncpt_two_step=.true.'   >> input_onetep_dmft.txt ")
      else
         call system("echo 'ncpt_two_step=.false.'  >> input_onetep_dmft.txt ")
      endif

      call system("echo 'nproc_onetep_openmp_='"//TRIM(ADJUSTL(toString(nproc_onetep_openmp)))//" >> input_onetep_dmft.txt ")
      if (lin_scaling) then
         call system(" update_case_file 'dmft_lin_scaling ' T ")
         call system("echo 'use_eimp_from_onetep=.true.'   >> input_onetep_dmft.txt ")
         call system("echo 'use_simp_from_onetep=.true.'   >> input_onetep_dmft.txt ")
      endif
      if (dmft_split) call system(" update_case_file 'dmft_split '  T ")
      if (dmft_splitk) call system(" update_case_file 'dmft_splitk ' T ")
      if (dmft_splitk) call system("echo 'dmft_splitkdmftall=.true.' >> input_onetep_dmft.txt ")

      if (real_axis_only_last_step .and. iter_dmft_sc >= niter_dft_dmft_sc - 1) then
         call system("echo 'ed_no_real_overide=.false.' >> input_onetep_dmft.txt ")
      endif

      call system("update_case_file 'dmft_nmu_loop ' "//trim(adjustl(toString(niter_dmft_mu))))
      call system("update_case_file 'dmft_sc ' F ")
      call system("echo 'nproc_onetep='"//TRIM(ADJUSTL(toString(nproc_gpu)))//" >> input_onetep_dmft.txt ")
      call system("echo 'just_onetep=.false.' >> input_onetep_dmft.txt ")
      call system("echo mach_onetep='machines_onetep_gf_ker' >> input_onetep_dmft.txt ")

      if (.not. dmft_split .and. .not. dmft_splitk .and. .not. no_frequ_split) then
         call system("echo 'split_onetep=.true.' >> input_onetep_dmft.txt ")
      else
         call system("echo 'split_onetep=.false.' >> input_onetep_dmft.txt ")
      endif
      call system(" update_case_file 'kernel_diis_max ' 0 ")
      call system(" update_case_file 'dmft_kernel '     0 ")

      if (iter_dmft_sc == niter_dft_dmft_sc .and. .not. matsu_solver) then
         call system("echo 'last_iter_is_real=.true.' >> input_onetep_dmft.txt ")
      endif

      !preparing input files for later DOS calculations
      call system(" cp input_onetep_dmft.txt input_onetep_dmft.txt.dos "); call system(" cp *.dat case.for.dos ")
      call system(" echo 'just_onetep=.true.' >> input_onetep_dmft.txt.dos ")
      !done

      if (.not. use_same_self_energy .or. (iter_dmft_sc == 1 .and. .not. debug_mode_erase_sigma)) then
         call system("onetep.dmft> sc_dmft_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc)))//" 2>&1 ")
         call system("mkdir       sc_dmft_full_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         if (copy_kernel) then
            call system("cp *.dkn sc_dmft_full_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc)))//" 2>/dev/null ")
         endif
         call system("cp *.dat *.txt sc_dmft_full_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc)))//" 2>/dev/null ")
        call system("cp chem.potential.nmu.iter* embedding_potentials* sigma_output* sc_dmft_full_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc)))//" 2>/dev/null ")
        call system("mv fort.* onetep_output* dir_onetep_iter* onetep_* sc_dmft_full_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc)))//" 2>/dev/null ")
         if (use_previous_dmft_files) then
          if( save_atoms_at_each_steps .or. iter_dmft_sc>=niter_dft_dmft_sc-1 ) call system("cp -r atom* sc_dmft_full_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         else
            call system("mv atom* sc_dmft_full_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         endif
         if (use_same_self_energy) then
            call system(" mkdir SIG_BACKUP ")
            call system(" cp sigma_output* SIG_BACKUP/ ")
         endif
      else
         call system("cp SIG_BACKUP/* ./ ")
         if (debug_mode_erase_sigma) call system("rm ./sigma_output* > /dev/null 2>&1 ")
      endif

      if (lin_scaling) call system(" update_case_file 'dmft_lin_scaling ' F ")
      if (dmft_split) call system(" update_case_file 'dmft_split '  F ")
      if (dmft_splitk) call system(" update_case_file 'dmft_splitk ' F ")
      if (dmft_splitk) call system("echo 'dmft_splitkdmftall=.false.' >> input_onetep_dmft.txt ")

      !-------------------------!
      !-------------------------!
      !-------------------------!
      !-------------------------!

      if (sandwitch_embedding .and. iter_dmft_sc == niter_dft_dmft_sc - 1) then
         write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         write (*, *) 'STEP4-embedding : generate the embedding self energies'
         write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         ! NOTE 1 : NOT compatible with dmft_split and dmft_splitk
         ! NOte 2 : NOT DOING THE FULL DMFT - JUST ADJUSTING CHEM AND DUMPING OUT EMBEDDING TABLES
         if (dmft_split .or. dmft_splitk) then
            write (*, *) 'ERROR embedding not compatible with splitk or split'; stop
         endif
         if (.not. dmft_split .and. .not. dmft_splitk .and. .not. no_frequ_split) then
            call system("echo 'split_onetep=.true.' >> input_onetep_dmft.txt ")
         else
            call system("echo 'split_onetep=.false.' >> input_onetep_dmft.txt ")
         endif
         call system("echo 'nproc_onetep_openmp_='"//TRIM(ADJUSTL(toString(nproc_onetep_openmp)))//" >> input_onetep_dmft.txt ")

         call system(" echo mach_onetep='machines_onetep_embed' >> input_onetep_dmft.txt ")
         call system(" update_case_file 'kernel_diis_max ' 0 ")
         call system(" update_case_file 'dmft_kernel '     0 ")
         call system(" echo 'just_onetep=.true.' >> input_onetep_dmft.txt ")
         !LEFT
         call system("mv embedding_dir/connections_D_L* . ")
         call system("onetep.dmft> sc_dmft_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc)))//" 2>&1 ")
         call system("mkdir sc_dmft_full_LEFT_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         call system("cp sigma_output* sc_dmft_full_LEFT_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         call system("mv _H* embedding_potentials* sc_dmft_full_LEFT_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
     call system("mv fort.* onetep_output* dir_onetep_iter* onetep_* sc_dmft_full_LEFT_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         call system("mv connections_D_L* embedding_dir/ ")
         !RIGHT
         call system("mv embedding_dir/connections_D_R* . ")
         call system("onetep.dmft> sc_dmft_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc)))//" 2>&1 ")
         call system("mkdir sc_dmft_full_RIGHT_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         call system("cp sigma_output* sc_dmft_full_RIGHT_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         call system("mv _H* embedding_potentials* sc_dmft_full_RIGHT_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
    call system("mv fort.* onetep_output* dir_onetep_iter* onetep_* sc_dmft_full_RIGHT_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         call system("mv connections_D_R* embedding_dir/ ")
      endif

      !-------------------------!
      !-------------------------!
      !-------------------------!
      !-------------------------!

      if (iter_dmft_sc < niter_dft_dmft_sc .or. matsu_solver) then ! last iter SIGMA is in real frequency
         write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         write (*, *) 'STEP4-postprocess : generate the density kernel associated to sigma'
         write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

         !for kernel mixing with last iteration!
         if (dmft_split) call system(" update_case_file 'dmft_split '  T ")
         if (dmft_splitk) call system(" update_case_file 'dmft_splitk ' T ")
         if (dmft_splitk) call system("echo 'dmft_splitkdmftall=.true.' >> input_onetep_dmft.txt ")

         if (.not. dmft_split .and. .not. dmft_splitk .and. .not. no_frequ_split) then
            call system("echo 'split_onetep=.true.' >> input_onetep_dmft.txt ")
         else
            call system("echo 'split_onetep=.false.' >> input_onetep_dmft.txt ")
         endif
         call system("echo 'nproc_onetep_openmp_='"//TRIM(ADJUSTL(toString(nproc_onetep_openmp)))//" >> input_onetep_dmft.txt ")
         call system("echo mach_onetep='machines_onetep_gf_ker' >> input_onetep_dmft.txt ")

         call system(" update_case_file 'dmft_nmu_loop ' "//trim(adjustl(toString(niter_kernel_mu))))
         call system(" update_case_file 'dmft_kernel ' "//trim(adjustl(toString(dmft_kernel_process))))
         call system(" update_case_file 'kernel_diis_max ' "//TRIM(ADJUSTL(toString(diis_max))))
         call system(" echo 'just_onetep=.true.' >> input_onetep_dmft.txt ")
         call system(" onetep.dmft> sc_dmft_kernel_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc)))//" 2>&1 ")
         call system(" mkdir sc_dmft_kernel_full_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         call system(" mv fort.* onetep_output* sc_dmft_kernel_full_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         call system(" cp *.dat *.txt sc_dmft_kernel_full_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc))))
         call system(" for i in `ls | grep  _diis_ ` ; do j=`echo $i | sed 's/_diis_/diis_/g' ` ; mv $i $j ; done ; ")
         if (dmft_split) call system(" update_case_file 'dmft_split '  F ")
         if (dmft_splitk) call system(" update_case_file 'dmft_splitk ' F ")
         if (dmft_splitk) call system("echo 'dmft_splitkdmftall=.false.' >> input_onetep_dmft.txt ")

         inquire (file='chem.potential.nmu.iter1', exist=check_mu)
         if (check_mu) then
            open (unit=81119, file='chem.potential.nmu.iter1')
            read (81119, *) chem_pots(iter_dmft_sc, 1)
            close (81119)
#ifdef _plot
            call plotarray((/(dble(iji), iji=1, iter_dmft_sc)/), chem_pots(1:iter_dmft_sc, 1), 'Potential spin 1')
#endif
         endif
         inquire (file='chem.potential.nmu.iter2', exist=check_mu)
         if (check_mu) then
            open (unit=81119, file='chem.potential.nmu.iter2')
            read (81119, *) chem_pots(iter_dmft_sc, 2)
            close (81119)
#ifdef _plot
            call plotarray((/(dble(iji), iji=1, iter_dmft_sc)/), chem_pots(1:iter_dmft_sc, 2), 'Potential spin 2')
#endif
         endif

      endif

      !-------------------------!
      !-------------------------!
      !-------------------------!
      !-------------------------!

      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      write (*, *) 'STEP5: CLEAN AFTER DMFT ITERATION AND COPY KERNEL'
      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'


      if (iter_dmft_sc == niter_dft_dmft_sc) then
         call system("      mkdir _store_last_iter  ")
         call system(" cp store_* _store_last_iter/ ")
      endif

      call system(" rm store_*        > /dev/null 2>&1 ")
      call system(" rm dens_kern_tmp* > /dev/null 2>&1 ")

      !-------------------------!
      !-------------------------!
      !-------------------------!
      !-------------------------!

      if (iter_dmft_sc < niter_dft_dmft_sc .or. matsu_solver) then ! last iter SIGMA is in real frequency
         write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         write (*, *) 'STEP6: COPY DMFT KERNEL TO DFT KERNEL  '
         write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

         call system(" mv "//TRIM(ADJUSTL(CASE_ONETEP))//".dkn "//TRIM(ADJUSTL(CASE_ONETEP))//".dkn_onetep ")
         if (.not. flag_turn_off_dmft) then
            call system(" cp "//TRIM(ADJUSTL(CASE_ONETEP))//".dkn_dmft "//TRIM(ADJUSTL(CASE_ONETEP))//".dkn     ")
         else
            call system(" cp "//TRIM(ADJUSTL(CASE_ONETEP))//".dkn_onetep "//TRIM(ADJUSTL(CASE_ONETEP))//".dkn   ")
         endif
      endif

      !-------------------------!
      !-------------------------!
      !-------------------------!
      !-------------------------!

      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      write (*, *) 'STEP7: RESTORE CONFIG FILES            '
      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

      call system(" cp  "//TRIM(ADJUSTL(CASE_ONETEP))//".dat.backup "//TRIM(ADJUSTL(CASE_ONETEP))//".dat     ")
      call system(" cp                 input_onetep_dmft.txt.backup input_onetep_dmft.txt       ")
      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

      !-------------------------!
      !-------------------------!
      !-------------------------!
      !-------------------------!

      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      write (*, *) 'STEP8: UPDATE NGWFs with new DMFT KERNEL if requested'
      write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

      if (NGWF_CGiter > 0 .and. iter_dmft_sc < niter_dft_dmft_sc) then
         call system(" update_case_file 'maxit_ngwf_cg '"//TRIM(ADJUSTL(toString(NGWF_CGiter))))
         call system(" update_case_file 'maxit_lnv '"//trim(adjustl(toString(kerneliter))))
         call system(" update_case_file 'minit_lnv '"//trim(adjustl(toString(kerneliter))))
         if (spoil) then
            call system(" update_case_file 'dmft_spoil_kernel ' T ")
         else
            call system(" update_case_file 'dmft_spoil_kernel ' F ")
         endif
         if (diis) then
            call system(" update_case_file 'kernel_diis ' T     ")
            call system(" update_case_file 'kernel_diis_maxit ' 1 ")
            call system(" update_case_file 'maxit_ngwf_diis '   1 ")
         else
            call system(" update_case_file 'kernel_diis ' F     ")
            call system(" update_case_file 'kernel_diis_maxit ' 0 ")
            call system(" update_case_file 'maxit_ngwf_diis '  0 ")
         endif
         call system(" update_case_file 'kernel_update '         F ")
         call system(" update_case_file 'read_denskern '         T ")
         call system(" update_case_file 'read_tightbox_ngwfs '   T ")
         call system(" update_case_file 'write_denskern '        F ")
         call system(" update_case_file 'write_tightbox_ngwfs '  T ")
         call system(" update_case_file 'write_forces '          T ")
         if (hydrogenic_projectors < 0.0) then
            call system(" update_case_file 'hubbard_proj_mixing ' -2.0 ")
            call system(" update_case_file_flip_hub                +10 ")
         else
            call system(" update_case_file 'hubbard_proj_mixing '  0.0 ")
            call system(" update_case_file_flip_hub "//trim(adjustl(tostring(hydrogenic_projectors))))
         endif
         if (numa) then
            prefix_ = " numactl  --localalloc "
         else
            prefix_ = " "
         endif
         mach_arg_ = ' -machinefile '
         mach_file_ = " machines_onetep_dft "
         nprocess_ = nproc_store
         args_ = trim(adjustl(CASE_ONETEP))
         output_ = " sc_ngwf_output_iter"//TRIM(ADJUSTL(toString(iter_dmft_sc)))
         command_line = build_mpi_command_line(prefix=prefix_, np=nprocess_, p=1, omp=1, mach=mach_file_,&
                                          & mach_arg=mach_arg_, args=args_, exe=exec_onetep,&
                                          & outputin=output_, hide_errors=.false., localhost='F')

         call system("ls machines_onetep_dft || fill_machine_file machines_onetep_dft "//TRIM(ADJUSTL(toString(nproc_store))))
         call utils_system_call(command_line, abort=.true.)

         !-------------------------!
         !-------------------------!
         !-------------------------!
         !-------------------------!

         write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         write (*, *) 'STEP9: RESTORE CONFIG FILES'
         write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         call system(" cp  "//TRIM(ADJUSTL(CASE_ONETEP))//".dat.backup "//TRIM(ADJUSTL(CASE_ONETEP))//".dat     ")
         write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         write (*, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'

      endif

   enddo
   !=========================================================================!
   !=========================================================================!
   !=========================================================================!
   !=========================================================================!

   write (*, *) '----------------------------------------------------------'
   write (*, *) '-------------------ENDING DMFT----------------------------'
   write (*, *) '----------------------------------------------------------'

   write (*, *) '----------------------------------------------------------'
   write (*, *) 'terminating the temperature monitoring of GPUs'
   if (monitor_gpu_temperature) call system("killall -9 monitor_temp.out")
   write (*, *) '----------------------------------------------------------'
   write (*, *) 'DMFT+ONETEP done, The End.'
   write (*, *) '----------------------------------------------------------'

contains

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

   character(2000) function build_mpi_command_line(prefix, np, p, omp, mach, mach_arg, args, exe, outputin, hide_errors, localhost)
      implicit none
      integer         ::  np, omp, p
      character*(*)   ::  mach, mach_arg, args, exe, outputin, prefix, localhost
      character(2000) ::  output
      logical         ::  hide_errors
      output = trim(adjustl(outputin))
      if (hide_errors) then
         output = trim(adjustl(output))//" 2> /dev/null "
      else
         output = trim(adjustl(output))//" 2>&1 "
      endif
      build_mpi_command_line = " "
      call system('  get_mpi_command_line '//trim(adjustl(exe))//'  "    '&
                                         & //trim(adjustl(args))//'  "    '&
                                         & //trim(adjustl(tostring(np)))//'       '&
                                         & //trim(adjustl(tostring(p)))//'       '&
                                         & //trim(adjustl(tostring(omp)))//'       '&
                                         & //trim(adjustl(mach))//'  "    '&
                                         & //trim(adjustl(mach_arg))//'  "  " '&
                                         & //trim(adjustl(output))//'  "    '&
                                         & //trim(adjustl(localhost))//'  > __mpi_cmd__  ')
      open (unit=9911, file='__mpi_cmd__')
      read (9911, '(a)') build_mpi_command_line
      close (9911)
      build_mpi_command_line = trim(adjustl(prefix))//" "//trim(adjustl(build_mpi_command_line))
      call system(" rm __mpi_cmd__ ")
   end function

   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!
   !-------------------------!

end program

!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

