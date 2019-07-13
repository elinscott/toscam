module globalvar_ed_solver

   use genvar, only: DBL

#ifdef ALLFIRSTCALL
   logical, parameter      :: ALL_FIRST_CALL = .true.
#else
   logical, parameter      :: ALL_FIRST_CALL = .false.
#endif
#ifdef debug
   logical, parameter      :: verboseall = .true.
#else
   logical, parameter      :: verboseall = .false.
#endif
   logical                 :: fmos, fmos_fluc, fmos_hub1
   integer                 :: fmos_iter, nvecout
   real(8)                 :: fmos_mix
   integer                 :: flag_idelta_two_scales_ed
   real(8), allocatable     :: UCCr(:, :, :)
   complex(8), allocatable :: UCCc(:, :, :)
   integer, allocatable    :: UCC(:, :, :, :)
   integer                 :: ncpt_chain_coup
   real(8)                 :: quench_mag, quench_U
   logical                 :: keldysh_pert_ground_sector, flag_slater_int, &
                              use_precomputed_slater_matrix, &
                              cpt_correct_green_out, ncpt_flag_two_step_fit, &
                              gen_cpt
   complex(8), allocatable :: Slater_Coulomb_c(:, :, :, :)
   real(8), allocatable    :: Slater_Coulomb_r(:, :, :, :)
   integer                 :: do_quench, do_quench_, quench_orb, ncpt_approx, &
                              ncpt_tot, ncpt_approx_tot, ncpt_para
   real(8)                 :: cpt_upper_bound, cpt_lagrange
   complex(8), allocatable :: epsilon_cpt(:, :), T_cpt(:, :)
   logical                 :: DO_NOT_USE_OPT_LANCZOS, do_keldysh_gbigger, &
                              quench_cancel_statistics
   logical                 :: donot_compute_holepart_spm
   integer                 :: keldysh_ortho
   real(8)                 :: keldysh_t0, keldysh_tmax, keldysh_delta
   integer                 :: keldysh_n
   integer                 :: freeze_poles_delta_iter, &
                              freeze_poles_delta_niter, Niter_search_max_0
   real(8)                 :: rdelta_frequ_eta1, rdelta_frequ_T, &
                              rdelta_frequ_w0, Jhund, rdelta_frequ_eta2
   real(8), allocatable    :: UUmatrix(:, :), JJmatrix(:, :)
   logical                 :: use_input_Eimp_matrices, &
                              fit_all_elements_show_graphs
   logical                 :: do_keldysh, FLAG_DUMP_INFO_FOR_GAMMA_VERTEX
   logical                 :: FLAG_NCUP, fit_green_func_and_not_hybrid
   integer                 :: iterdmft, niterdmft, iwindow
   integer                 :: cuda_blocksize
   integer, allocatable    :: MASK_AVERAGE(:, :), MASK_AVERAGE_SIGMA(:, :)
   logical                 :: freeze_poles_delta, FLAG_GUP_IS_GDN, &
                              always_compute_static_obs
   integer                 :: FLAG_MPI_GREENS
   real(DBL)               :: freeze_pole_lambda, beta_ED, beta_ED_, &
                              lambda_sym_fit
   real(DBL)               :: weight_expo, cutoff_dynamic, &
                              cutoff_min_lanczos_vec
   real(DBL), allocatable  :: FROZEN_POLES(:)
   character(20)           :: which_lanczos
   logical                 :: OPEN_MP, ON_FLY, USE_TRANSPOSE_TRICK_MPI, &
                              restarted, only_compute_density, skip_fit_
   logical                 :: start_para, PAIRING_IMP_TO_BATH, &
                              use_specific_set_parameters, force_no_bcs_pairing
   logical                 :: superconducting_state
   logical, parameter      :: Jhund_Slater_type = .true.
   real(DBL), allocatable  :: param_input(:), param_output(:)
   integer                 :: average_G
   logical                 :: diag_bath, diag_V, bath_nearest_hop, &
                              FLAG_ALL_GREEN_FUNC_COMPUTED
   logical                 :: FLAG_BUILD_CORREL_LOW_PART, dump_ground_state
   logical                 :: use_cuda_lanczos = .false.
   integer                 :: szmin, szmax, nupmin, nupmax, ndnmin, ndnmax
   logical                 :: track_sectors = .false., FLAG_FULL_ED_GREEN = .false.
   logical, parameter      :: USE_CC = .true.
   logical                 :: SCAN_FULL_NUP_NDN
   integer                 :: Niter_search_max, min_all_bath_param, fit_nw, istati
   logical                 :: supersc_state
   real(DBL)               :: tot_repulsion
   real(DBL)               :: dist_max, search_step, fit_weight_power, cutoff_rvb, &
                              cutoff_hamilt_param
   complex(DBL)            :: energy_global_shift, energy_global_shift2
   character(20)           :: FIT_METH
   integer                 :: Nitergreenmax
   integer                 :: Nitermax                    ! max # of Lanczos iterations
   logical                 :: force_pairing_from_mask
   logical                 :: para_state, donot_compute_holepart, &
                              force_sz_basis, force_nupdn_basis
   logical                 :: fast_fit, first_iter_use_edinput, &
                              force_no_pairing, force_para_state, &
                              force_singlet_state
   real(DBL)               :: tolerance               ! Lanczos tolerance
   real(DBL)               :: dEmax = 0.d0, dEmax0 = 0.d0 ! max.energy of excited states to consider
   integer                 :: Neigen = 0   ! max. # of eigenvalues computed in this window
   integer                 :: Block_size = 0   ! Block size (ignored if Neigen=1: Block_size=1)
   integer                 :: nsec = 1               ! Number of sectors to scan
   integer                 :: nsec0 = 1               ! Scan of the sectors start at sec=sec0
   integer                 :: window_hybrid = 0       ! window of matsubara frequencies to fit
   integer                 :: window_hybrid2 = 0      ! window of matsubara frequencies to fit
   integer                 :: window_weight = 1       ! ratio weight point inside window / point outside window
   real(DBL), allocatable  :: dens(:)                 ! onsite charge density
   integer, allocatable    :: list_sectors(:)
   logical                 :: verbose_graph           ! plot fits of the hybridization at each step of the minimization
   real(DBL)               :: fit_shift
   integer                 :: ed_num_eigenstates_print
end module
