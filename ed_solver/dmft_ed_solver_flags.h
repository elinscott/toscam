
   if(flag_idelta_two_scales_ed/=0)then
     write(*,*) 'SETTING TWO ENERGY SCALE BROADENING'
     rdelta_frequ_eta2 = dble(flag_idelta_two_scales_ed)
    else
     write(*,*) 'SETTING SINGLE ENERGY SCALE BROADENING'
     rdelta_frequ_eta2 = 0.d0
   endif

   if(MAXT>1) then
     write(*,*) 'USING OPENMP WITH [x] THREADS : ', MAXT
     OPEN_MP=.true.
   else
     write(*,*) 'NOT USING OPENMP'
     OPEN_MP=.false.
   endif
   if(MAXT>1.and.size2>1) write(*,*) 'WARNING : TRY TO RUN MPI AND OPEN-mP AT THE SAME TIME'
   if(OPEN_MP) then
     if(.not.USE_CC) stop 'error cannot use OPEN_MP if USE_CC is off'
   endif
   if(use_cuda_lanczos.or.which_lanczos=='GPU') then
     if(.not.USE_CC) stop 'error cannot use GPU if USE_CC is off'
   endif
   if(USE_TRANSPOSE_TRICK_MPI.and..not.enable_mpi_dot)then
     stop 'ERROR USE_TRANSPOSE_TRICK_MPI needs enable_mpi_dot enabled, because the lanczos vector are split among the cpu'
   endif
   if(USE_TRANSPOSE_TRICK_MPI.and..not.USE_CC)then
     stop 'ERROR USE_TRANSPOSE_TRICK_MPI needs USE_CC'
   endif
   if(FLAG_ALL_GREEN_FUNC_COMPUTED)then
    FLAG_BUILD_CORREL_LOW_PART=.true.
   endif
   if(FLAG_MPI_GREENS>0) then
     FLAG_FORBID_SPLITTING=.true.
     enable_mpi_dot=.false.
   endif

   if(abs(beta_ED_)>1.d-7) beta_ED=beta_ED_

   call init_frequency_module(rdelta_frequ_eta1,rdelta_frequ_w0,rdelta_frequ_T,rdelta_frequ_eta2_=rdelta_frequ_eta2)

   if(Niter_search_max_0==0) Niter_search_max_0=Niter_search_max

   if(ncpt_approx/=0)then
      PAIRING_IMP_TO_BATH=.false.
      write(*,*) 'WARNING turning off pairing imp to bath, not compatible with CPT approximation yet'
   endif

   do_quench_=do_quench

   if(keldysh_pert_ground_sector) then
      write(logfile,*) 'WARNING : for keldysh, setting the flag keldysh_pert_ground_sector'
      write(logfile,*) '          this sets dEmax0 to a small value, in order to only get the ground state in each sector'
      write(logfile,*) '          and this also cancel the flag quench_cancel_statistics'
      dEmax0=0.00001
      quench_cancel_statistics=.true.
   else
      do_quench=0
   endif

   if(fmos)then
    ncpt_approx=0
    fast_fit=.false.
    diag_bath=.true.
    diag_V=.true.
    bath_nearest_hop=.false. 
    min_all_bath_param=min_all_bath_param-mod(min_all_bath_param,Nc)
   endif


