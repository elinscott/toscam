! Parameters for ED solver, listed as:
!   - variable name
!   - variable label
!   - default
!   - description

call putel_in_namelist(nm, &
                       which_lanczos, &
                       'which_lanczos', &
                       'NORMAL', &
                       'NORMAL, FULL_ED, or ARPACK')
!   'NORMAL, GPU, FULL_ED, or ARPACK')

call putel_in_namelist(nm, &
                       EDfile, &
                       'EDfile', &
                       './ED/ED.in', &
                       'ED-SOLVER PARAMETERS FILE')

call putel_in_namelist(nm, &
                       BATHfile, &
                       'BATHfile', &
                       TRIM(ADJUSTL('./ED/'//fileb)), &
                       'BATH PARAMETERS FILE')

call putel_in_namelist(nm, &
                       CORRELfile, &
                       'CORRELfile', &
                       TRIM(ADJUSTL('./ED/'//filec)), &
                       'CORRELATIONS PARAMETERS FILE')

call putel_in_namelist(nm, &
                       OLDGSFILE, &
                       'OLDGSFILE', &
                       './ED_out/GS.raw', &
                       'OLD GS FILE')

call putel_in_namelist(nm, &
                       fit_weight_power, &
                       'fit_weight_power', &
                       0.5d0, &
                       'Bath fit : give the exponent of the 1/w**a to weight the frequencies to fit')

call putel_in_namelist(nm, &
                       fit_shift, &
                       'fit_shift', &
                       0.01d0, &
                       'Bath fit : give the shift of the 1/(w**a+shift) fit')

call putel_in_namelist(nm, &
                       start_from_old_gs, &
                       'start_from_old_gs', &
                       .false., &
                       'START FROM OLD GS')

call putel_in_namelist(nm, &
                       Nitergreenmax, &
                       'Nitergreenmax', &
                       36, &
                       'max.numb.of Lanczos iterations in Green s fct. computation')

call putel_in_namelist(nm, &
                       Niter_search_max, &
                       'Niter_search_max', &
                       100, &
                       'Niter_search_max = max. numb. of function calls in conjugate gradient')

call putel_in_namelist(nm, &
                       search_step, &
                       'search_step', &
                       1.d-4, &
                       'small step in search direction')

call putel_in_namelist(nm, &
                       dist_max, &
                       'dist_max', &
                       1.d-10, &
                       'max. error on hybridization functions')

call putel_in_namelist(nm, &
                       nsec, &
                       'nsec  ', &
                       -1, &
                       'number of sectors to diagonalize, -1 to compute them all')

call putel_in_namelist(nm, &
                       nsec0, &
                       'nsec0 ', &
                       1, &
                       'If list_sectors not provided, start to scan at nsec0 and stop to scan at nsec0+nsec')

call putel_in_namelist(nm, &
                       dEmax0, &
                       'dEmax0', &
                       5.d0, &
                       'max. energy of excited state to consider, normalized by beta')

call putel_in_namelist(nm, &
                       Neigen, &
                       'Neigen', &
                       1, &
                       'number of eigenvalues to compute')

call putel_in_namelist(nm, &
                       Nitermax, &
                       'Nitermax  ', &
                       128, &
                       'max. number of Lanczos iterations')

call putel_in_namelist(nm, &
                       Block_size, &
                       'Block_size', &
                       0, &
                       'Block size (0 lets the routine decide), -1 use Cullum instead')

call putel_in_namelist(nm, &
                       tolerance, &
                       'tolerance ', &
                       1.d-12, &
                       'Lanczos convergence tolerance')

call putel_in_namelist(nm, &
   FIT_METH, &
   'FIT_METH', &
   'MINIMIZE', &
   'FITTING METHODS, KEYWORD CAN BE (default is MINIMIZE) : POWELL, CONJGGRAD, &
   &MINIMIZE, CIVELLI, BFGS, LIN_APPROX, POWELL_OPT, CONJG_GRAD, BFGS_, CG1, &
   &CG2, CG3')

call putel_in_namelist(nm, &
                       verbose_graph, &
                       'verbose_graph', &
                       .false., &
                       'plot hybridization fit at each step of the minimization')

call putel_in_namelist(nm, &
                       window_hybrid, &
                       'window_hybrid', &
                       0, &
                       'window (left)  of matsubara frequencies to fit the hybridization')

call putel_in_namelist(nm, &
                       window_hybrid2, &
                       'window_hybrid2', &
                       0, &
                       'window (right) of matsubara frequencies to fit the hybridization')

call putel_in_namelist(nm, &
                       window_weight, &
                       'window_weight', &
                       1, &
                       'ratio weight point inside window over point outside window')

call putel_in_namelist(nm, &
   min_all_bath_param, &
   'min_all_bath_param', &
   8, &
   'if /= 0 will overtake the ed_hybrid file and minimize all the bath &
   &parameters, excluding (including) superconductivity for >0 (<0). the &
   &value is the number of site in the bath. it will adjust for para/sc/not &
   &sc/not para depending on the chosen order type in the main library.')

call putel_in_namelist(nm, &
   force_sz_basis, &
   'force_sz_basis', &
   .false., &
   'Force the use of the Sz quantum number (mandatory for SC state, since &
   &number of part. not fixed) instead of using the nup/ndn quantum numbers &
   &(ok for paramagn. state)')

call putel_in_namelist(nm, &
                       force_nupdn_basis, &
                       'force_nupdn_basis', &
                       .false., &
                       'Force the use of the nup ndn quantum number')

call putel_in_namelist(nm, &
   first_iter_use_edinput, &
   'first_iter_use_edinput', &
   .false., &
   'start the dmft iterations with the initial bath parameters as a starting &
   &bath instead of the Delta(iw) input ')

call putel_in_namelist(nm, &
                       force_no_pairing, &
                       'force_no_pairing', &
                       .false., &
                       'if true kills the pairing')

call putel_in_namelist(nm, &
                       force_para_state, &
                       'force_para_state', &
                       .false., &
                       'force the paramagnetic symmetry (up/dn spin the same)')

call putel_in_namelist(nm, &
   SCAN_FULL_NUP_NDN, &
   'SCAN_FULL_NUP_NDN', &
   .true., &
   'if true, will scan all nup ndn sectors (also with nup/=ndn); if false, &
   &will only scan nup=ndn sectors')

call putel_in_namelist(nm, &
                       fit_nw, &
                       'fit_nw', &
                       650, &
                       'only consider the fit_nw matsubara frequencies for the fit')

call putel_in_namelist(nm, &
                       FLAG_BUILD_CORREL_LOW_PART, &
                       'FLAG_BUILD_CORREL_LOW_PART', &
                       .true., &
                       'build automatically the lower part of the correl mask')

call putel_in_namelist(nm, &
   cutoff_rvb, &
   'cutoff_rvb', &
   0.01d0, &
   'cutoff under which it is considered no more as rvb state but as normal &
   &state, under this cutoff the inversion of the Green Function matrix could &
   &lead to artificial strong anomalous weight in the self energt')

call putel_in_namelist(nm, &
                       dump_ground_state, &
                       'dump_ground_state', &
                       .false., &
                       'write the ground state to a file during the Lanczos process')

! ebl: removing GPU functionality, &
! call putel_in_namelist(nm, &
!    use_cuda_lanczos, &
!    'use_cuda_lanczos', &
!    .false., &
!    'run Lanczos on the GPU if present in the hardware')

call putel_in_namelist(nm, &
                       weight_expo, &
                       'weight_expo', &
                       2.d0, &
                       'exponent power of the fitting difference |a-b|^weight_exop')

! ebl: removing GPU functionality, &
! call putel_in_namelist(nm, &
!    cuda_blocksize, &
!    'cuda_blocksize', &
!    8, &
!    'the number of GPU threads working in cuda')

call putel_in_namelist(nm, &
   track_sectors, &
   'track_sectors', &
   .false., &
   'if true the code is tracking the ground state sector along the DMFT &
   &iterations, and only scanning the previous sector +- 1 sectors')

call putel_in_namelist(nm, &
                       ON_FLY, &
                       'ON_FLY', &
                       .false., &
                       ' do not store the sparse hamiltonian matrix for the Sz case')

call putel_in_namelist(nm, &
   diag_bath, &
   'diag_bath', &
   .false., &
   'takes into account only the diagonal elements of the bath, no bath to bath &
   &hoppings')

call putel_in_namelist(nm, &
   bath_nearest_hop, &
   'bath_nearest_hop', &
   .false., &
   ' whatever bath parametrization, includes the bath to bath nearest neighbor &
   &hoppings')

call putel_in_namelist(nm, &
                       USE_TRANSPOSE_TRICK_MPI, &
                       'USE_TRANSPOSE_TRICK_MPI', &
                       .false., &
                       'each cpu keeps only a chunk of the lanczos vector in memory')

call putel_in_namelist(nm, &
                       start_para, &
                       'start_para', &
                       .false., &
                       'force starting the first iteration fit from a paramagnetic state')

call putel_in_namelist(nm, &
                       cutoff_dynamic, &
                       'cutoff_dynamic', &
                       1.d-8, &
                       'cutoff for computing dynamical green function')

call putel_in_namelist(nm, &
                       cutoff_min_lanczos_vec, &
                       'cutoff_min_lanczos_vec', &
                       1.d-25, &
                       'cutoff for the minimal norm of a vector obtained in Lanczos')

call putel_in_namelist(nm, &
   cutoff_hamilt_param, &
   'cutoff_hamilt_param', &
   1.d-3, &
   'cutoff_hamilt_param: under this param the hamiltonian parameters are &
   &considered as 0')

call putel_in_namelist(nm, &
   FLAG_ALL_GREEN_FUNC_COMPUTED, &
   'FLAG_ALL_GREEN_FUNC_COMPUTED', &
   .false., &
   'will compute the full Green function correlations, and override the MASK &
   &given in ed_correl, the lower part will be symmetric though')

call putel_in_namelist(nm, &
   FLAG_FULL_ED_GREEN, &
   'FLAG_FULL_ED_GREEN', &
   .false., &
   'if set to true, the code with also use full diagonalization for getting &
   &the green function, and therefore no states will be filtered out of the &
   &sectors because of the temperature cutoff')

call putel_in_namelist(nm, &
                       PAIRING_IMP_TO_BATH, &
                       'PAIRING_IMP_TO_BATH', &
                       .false., &
                       'pairing between impurity and bath')

call putel_in_namelist(nm, &
   force_no_bcs_pairing, &
   'force_no_bcs_pairing', &
   .false., &
   'if true will run the superconducting code, but imposing the &
   &superconducting variables/parameters to be 0, the purpose is for testing')

call putel_in_namelist(nm, &
                       fast_fit, &
                       'fast_fit', &
                       .false., &
                       'use a faster way to obtain the fit through eigenvalue decomposition')

call putel_in_namelist(nm, &
                       freeze_pole_lambda, &
                       'freeze_pole_lambda', &
                       1.d0, &
                       'between 0 and 1, gives the shift of the frozen poles')

call putel_in_namelist(nm, &
                       force_singlet_state, &
                       'force_singlet_state', &
                       .true., &
                       'force the singlet state in the superconducting calculations')

call putel_in_namelist(nm, &
   force_pairing_from_mask, &
   'force_pairing_from_mask', &
   .false., &
   'force the use of the pairing defined in the MASK of ed_hybrid instead of &
   &relaxing all parameters')

call putel_in_namelist(nm, &
                       FLAG_MPI_GREENS, &
                       'FLAG_MPI_GREENS', &
                       0, &
                       'compute each green functions on a different CPU')

call putel_in_namelist(nm, &
                       beta_ED_, &
                       'beta_ED_', &
                       0.d0, &
                       'temperature for Boltzman weight. If =0. it will use the DMFT temperature')

call putel_in_namelist(nm, &
   lambda_sym_fit, &
   'lambda_sym_fit', &
   0.d0, &
   'coefficient that contributes to the fitting distance for finding the AIM &
   &parameters. Basically, it adds a term for how far the corresponding &
   &parameters leads to a assymetric solution')

call putel_in_namelist(nm, &
   flag_introduce_noise_in_minimization, &
   'flag_introduce_noise_in_minimization', &
   .false., &
   'if yes it will introduce some noise in the minimization between 2 dmft &
   &iterations')

call putel_in_namelist(nm, &
   flag_introduce_only_noise_in_minimization, &
   'flag_introduce_only_noise_in_minimization', &
   .false., &
   'if yes it will start from a fresh random guess for the bath minimization &
   &at every DMFT iteration')

call putel_in_namelist(nm, &
   iwindow, &
   'iwindow', &
   1, &
   'for track_sectors, the window +- iwindow around the last iterations kept &
   &states, where we should look at for the new low energy state at the &
   &present iteration')

call putel_in_namelist(nm, &
   FLAG_GUP_IS_GDN, &
   'FLAG_GUP_IS_GDN', &
   .false., &
   'take care, if true it will explicitely enforce Gup to be the same as Gdn, &
   &useful to speed up calculations for force_para_state')

call putel_in_namelist(nm, &
   always_compute_static_obs, &
   'always_compute_static_obs', &
   .true., &
   'if true will always compute static bosonic observables, whatever flags in &
   &ed_correl for N, Sz ...')

call putel_in_namelist(nm, &
   Niter_search_max_0, &
   'Niter_search_max_0', &
   0, &
   'if non zero, will use Niter_search_max_0 steps for the first iter for the &
   &minimization, if 0 Niter_search_max_0=Niter_search_max')

call putel_in_namelist(nm, &
                       FLAG_NCUP, &
                       'FLAG_NCUP', &
                       .false., &
                       'compute cor hopping green function n_dn C_iup, _iup^dagger')

call putel_in_namelist(nm, &
   FLAG_DUMP_INFO_FOR_GAMMA_VERTEX, &
   'FLAG_DUMP_INFO_FOR_GAMMA_VERTEX', &
   .false., &
   'if true will dump out the necessary files to compute Gamma (four leg &
   &vertex) for later postprocessing, only works with FULL_ED on')

call putel_in_namelist(nm, &
   fit_all_elements_show_graphs, &
   'fit_all_elements_show_graphs', &
   .false., &
   'if true will show the fit in agr files for all matrix elements, not only &
   &the diagonal ones')

call putel_in_namelist(nm, &
                       keldysh_t0, &
                       'keldysh_t0  ', &
                       0.d0, &
                       'keldysh : t0')

call putel_in_namelist(nm, &
                       keldysh_tmax, &
                       'keldysh_tmax', &
                       10.d0, &
                       'keldysh : tmax  ')

call putel_in_namelist(nm, &
                       keldysh_n, &
                       'keldysh_n   ', &
                       5, &
                       'keldysh : number of time frames ')

call putel_in_namelist(nm, &
   do_keldysh, &
   'do_keldysh', &
   .false., &
   'keldysh : if true will run the keldysh calculations, NB. the code should &
   &be compiled in complex mode with the -D_complex flag in the compilation line')

call putel_in_namelist(nm, &
                       DO_NOT_USE_OPT_LANCZOS, &
                       'DO_NOT_USE_OPT_LANCZOS', &
                       .false., &
                       'if true will not use the optimized version of the Lanczos algorithm')

call putel_in_namelist(nm, &
   keldysh_ortho, &
   'keldysh_ortho', &
   0, &
   'keldysh : if positive number will re-orthogonalize the Lanczos vector when &
   &building the local Kryslov space for the time evolution, it will use the &
   &[x] first Lanczos vector, if 0 it will not use any re-orthogonalization, if &
   &-1 it will use all Lanczos vector')

call putel_in_namelist(nm, &
   do_quench, &
   'do_quench', &
   0, &
   'keldysh : if 1 will do a quench in magnetic field, if 2 in interaction, &
   &for the keldysh evolution, if 0 no quench is used ')

call putel_in_namelist(nm, &
                       quench_mag, &
                       'quench_mag', &
                       1.d0, &
                       'keldysh : magnetic field along z axis used for the quench')

call putel_in_namelist(nm, &
                       quench_U, &
                       'quench_U', &
                       10.d0, &
                       'keldysh : Coulomb repulsion used for the quench')

call putel_in_namelist(nm, &
                       do_keldysh_gbigger, &
                       'do_keldysh_gbigger', &
                       .true., &
                       'keldysh : if true will also compute G^>, if false will only compute G^<')

call putel_in_namelist(nm, &
   quench_orb, &
   'quench_orb', &
   -1, &
   'keldysh : which orbital should be computed for the evolution, if -1 all &
   &orbitals are computed')

call putel_in_namelist(nm, &
   quench_cancel_statistics, &
   'quench_cancel_statistics', &
   .false., &
   'keldysh : if true will not use the initial statistics from the modified &
   &initial hamiltonian, but a constant statistics, the idea is to start with &
   &the sectors of the true ground state of the unperturbated hamiltonian')

call putel_in_namelist(nm, &
   donot_compute_holepart_spm, &
   'donot_compute_holepart_spm', &
   .true., &
   'if false will compute both hole and particle part of S+S-, it is redundant &
   &but useful for debugging purposes')

call putel_in_namelist(nm, &
   ncpt_approx, &
   'ncpt_approx', &
   0, &
   'this is an additional degree of approximation, if larger than 0 it will &
   &add ncpt sites next to each bath site, and this will be treated by cluster &
   &perturbation theory, it is still experimental, please check this level of &
   &approximation')

call putel_in_namelist(nm, &
   cpt_upper_bound, &
   'cpt_upper_bound', &
   10000.d0, &
   'this is the upper bound on the V cpt parameters, which connect the cluster &
   &to the cpt legs')

call putel_in_namelist(nm, &
   cpt_lagrange, &
   'cpt_lagrange', &
   0.00d0, &
   'this is the lagrangian parameter for the weight V connecting to the large &
   &cluster, they should be as small as possible')

call putel_in_namelist(nm, &
   cpt_correct_green_out, &
   'cpt_correct_green_out', &
   .true., &
   'if true will correct the output impurity green function with the CPT &
   &formula, when the cpt approx is used')

call putel_in_namelist(nm, &
   ncpt_flag_two_step_fit, &
   'ncpt_flag_two_step_fit', &
   .true., &
   'if true will do the ncpt fit in two steps : first ncpt is turned off, the &
   &fit is carried out to get the AIM parameter, then ncpt is switched on, the &
   &fit is done again by including this time the ncpt parameters')

call putel_in_namelist(nm, &
   keldysh_pert_ground_sector, &
   'keldysh_pert_ground_sector', &
   .false., &
   'if true, the perturbated state is in the sectors of the perturbated &
   &Hamiltonian, if false in the sectors of the unperturbated hamiltonian')

call putel_in_namelist(nm, &
   gen_cpt, &
   'gen_cpt', &
   .false., &
   'if true, a generalized version of cpt is used, where the V_cpt connects &
   &all the sites of the chains to the bath orbitals')

call putel_in_namelist(nm, &
   fmos, &
   'fmos', &
   .false., &
   'Definition: if true will run a mean field solver instead of ED, useful if &
   &too many bath sites are required to get a good fit ')

call putel_in_namelist(nm, &
                       fmos_iter, &
                       'fmos_iter', &
                       1, &
                       'Definition: number of mean field iterations')

call putel_in_namelist(nm, &
                       diag_V, &
                       'diag_V', &
                       .false., &
                       'Definition: if true bath V is diagonal, for fmos in particular')

call putel_in_namelist(nm, &
                       fmos_mix, &
                       'fmos_mix', &
                       0.3d0, &
                       'Definition: mixing for FMOS iterations')

call putel_in_namelist(nm, &
                       fmos_fluc, &
                       'fmos_fluc', &
                       .false., &
                       'Definition: if true takes also into account the charge fluctuations in the FMOS approximation')

call putel_in_namelist(nm, &
                       fmos_hub1, &
                       'fmos_hub1', &
                       .false., &
                       'Definition:  if true uses the Hubbard I solver as a fast multi-orbital solver')

call putel_in_namelist(nm, &
                       flag_idelta_two_scales_ed, &
                       'flag_idelta_two_scales_ed', &
                       0, &
                       'Definition: if true the ED solver will use a double energy scale for the idelta frequency dependence')

call putel_in_namelist(nm, ed_num_eigenstates_print, 'ed_num_eigenstates_print', &
                      16, 'Definition: the number of lowest-energy eigenstates of &
                      &the reduced spectral density to list in the logfile report')
                  
call putel_in_namelist(nm, print_qc, 'print_qc', &
                      .false., 'Definition: if true, print out various quantities for the purposes of quality control')
