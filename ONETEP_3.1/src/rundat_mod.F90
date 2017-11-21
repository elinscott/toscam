! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-

module rundat

  use constants, only: DP, HARTREE_IN_EVS

  implicit none

  ! ndmh: General properties
  character(len=80) :: pub_rootname
  character(len=80) :: task
  real(kind=DP)     :: cutoff_energy
  real(kind=DP)     :: kernel_cutoff
  character(len=80) :: xc_functional
  integer           :: pub_libxc_x_func_id
  integer           :: pub_libxc_c_func_id
  real(kind=DP)     :: pub_charge
  integer           :: pub_spin
  logical           :: pub_spin_polarised
  real(kind=DP)     :: pub_constant_efield(3)
  integer           :: fftbox_pref(3)
  integer           :: ppd_npoints(3)
  real(kind=DP)     :: psinc_spacing(3)
  logical           :: pub_old_input
  character(len=80) :: pub_devel_code
  logical           :: pub_nnho
  real(kind=DP)     :: pub_ngwf_halo
  logical           :: pub_nonsc_forces
  real(kind=DP)     :: pub_external_pressure
  real(kind=DP)     :: pub_smoothing_factor
  real(kind=DP)     :: pub_isosurface_cutoff


  ! ndmh: Kernel Optimisation parameters
  integer           :: pub_maxit_palser_mano
  integer           :: maxit_pen
  real(kind=DP)     :: pen_param
  integer           :: maxit_lnv
  integer           :: minit_lnv
  integer           :: maxit_lnv_coarse
  real(kind=DP)     :: lnv_threshold_orig
  character(len=80) :: lnv_cg_type
  real(kind=DP)     :: lnv_cg_max_step
  logical           :: exact_lnv
  logical           :: old_lnv
  logical           :: pub_lnv_check_trial_steps
  integer           :: pub_kerfix
  integer           :: pub_maxit_kernel_fix
  logical           :: pub_initial_dens_realspace
  logical           :: pub_kernel_diis          !use kernel DIIS
  character(len=1)  :: pub_kernel_diis_type     !type of DIIS (D=diag only, P=Pulay, L=linear)
  integer           :: pub_kernel_diis_max      !max # of kernels to store in memory
  integer           :: pub_kernel_diis_maxit    !max # of DIIS iterations
  real(kind=DP)     :: pub_kernel_diis_thres    !DIIS convergence threshold
  integer           :: pub_kernel_diis_liter    !#linear DIIS iter before Pulay DIIS
  real(kind=dp)     :: pub_kernel_diis_c_in     !coefficients for linear DIIS
  real(kind=dp)     :: pub_kernel_diis_c_out    !
  character(len=4)  :: pub_kernel_diis_criteria !convergence criteria

  ! ndmh: NGWF optimisation parameters
  integer           :: maxit_ngwf_cg
  character(len=80) :: ngwf_cg_type
  real(kind=DP)     :: ngwf_cg_max_step
  integer           :: pub_elec_cg_max
  character(len=80) :: precond_scheme
  real(kind=DP)     :: k_zero
  real(kind=DP)     :: r_precond
  logical           :: precond_recip
  logical           :: precond_real
  character(len=80) :: smooth_scheme
  real(kind=DP)     :: r_smooth
  real(kind=DP)     :: k_smooth
  real(kind=DP)     :: occ_mix
  logical           :: pub_kernel_update
  logical           :: pub_kernel_christoffel_update
  integer           :: maxit_hotelling
  real(kind=DP)     :: max_resid_hotelling
  integer           :: maxit_ngwf_diis
  integer           :: ngwf_diis_size
  integer           :: ngwf_diis_reset
  real(kind=DP)     :: ngwf_diis_max_step
  logical           :: pub_use_aux_ngwfs

  ! ndmh: Convergence criteria
  real(kind=DP)     :: ngwf_threshold_orig
  real(kind=DP)     :: pub_elec_energy_tol   ! Max energy/atom change per iteration
  real(kind=DP)     :: pub_elec_force_tol    ! Max force change per iteration
  real(kind=DP)     :: pub_ngwf_max_grad     ! Maximum permissable NGWF gradient
  logical           :: delta_e_conv

  ! ndmh: Calculation parameters not directly related to physics
  integer           :: locpot_int_batch_size
  integer           :: kinetic_int_batch_size
  integer           :: density_batch_size
  integer           :: ngwf_grad_batch_size
  real(kind=DP)     :: pub_dense_threshold
  logical           :: pub_ovlp_for_nonlocal
  logical           :: use_space_filling_curve
  logical           :: coreham_denskern_guess
  logical           :: pub_check_atoms
  character(len=80) :: pub_locpot_scheme
  real(kind=DP)     :: pub_smooth_projectors
  real(kind=DP)     :: pub_smooth_loc_pspot
  logical           :: pub_odd_psinc_grid
  logical           :: pub_even_psinc_grid
  logical           :: pub_realspace_projectors

  ! ndmh: Output/Visualisation parameters
  integer           :: pub_output_detail
  integer           :: pub_paw_output_detail
  integer           :: timings_level
  logical           :: pub_write_params
  logical           :: pub_write_forces
  logical           :: pub_write_positions
  logical           :: pub_write_xyz
  logical           :: print_qc
  logical           :: pub_cube_format
  logical           :: pub_dx_format
  logical           :: pub_dx_coarse
  integer           :: pub_dx_sig_digits
  logical           :: pub_grd_format
  logical           :: pub_write_density_plot
  logical           :: pub_write_ngwf_plot
  logical           :: pub_write_ngwf_grad_plot
  integer           :: pub_write_ngwf_radial
  integer           :: pub_write_ngwf_grad_radial
  real(kind=DP)     :: pub_write_radial_step
  real(kind=DP)     :: pub_write_radial_smear

  ! NGWF/kernel I/O parameters
  logical           :: read_denskern
  logical           :: write_denskern
  logical           :: write_tightbox_ngwfs
  logical           :: read_tightbox_ngwfs
  logical           :: pub_write_converged_dk_ngwfs
  logical           :: pub_write_sw_ngwfs
  logical           :: pub_read_sw_ngwfs
  integer           :: pub_write_max_l
  integer           :: pub_read_max_l
  integer           :: pub_extra_n_sw

  ! ndmh: System parameters relating to grids, augmentation and pseudopotentials
  logical           :: pub_paw
  real(kind=DP)     :: pub_fine_grid_scale
  real(kind=DP)     :: pub_dbl_grid_scale
  logical           :: pub_aug_funcs_recip
  logical           :: pub_usp         ! not input
  logical           :: pub_aug         ! not input
  logical           :: pub_any_nl_proj ! not input
  logical           :: pub_nlcc        ! not input
  logical           :: pub_fine_is_dbl ! not input
  logical           :: pub_dbl_is_std  ! not input
  integer           :: pub_aug_den_dim ! not input
  logical           :: pub_nhat_in_xc  ! not input

  ! qoh: for dispersion correction
  integer           :: pub_dispersion  ! damping funcion to use
  real(kind=DP)     :: pub_vdw_dcoeff  ! override of damping coefficient

  ! ddor: parameters for DFT+U with self-consistent projectors
  logical           :: pub_hubbard
  integer           :: pub_hub_max_iter      ! Maximum number of DFT+U projector optimisation steps
  real(kind=DP)     :: pub_hub_energy_tol    ! Energy tolerance when using DFT+U projector optimisation
  integer           :: pub_hub_conv_win      ! Energy convergence window when using DFT+U projector optimisation
  real(kind=DP)     :: pub_hub_proj_mixing   ! Proportion of old Hubbard projector to mix with new
  integer           :: pub_hub_functional    ! DFT+U correction functional to use
  integer           :: pub_hub_tensor_corr   ! DFT+U tensorial correction to use
  real(kind=DP)     :: pub_hub_ngwf_spin_thr ! NGWF RMS gradient at which to switch off DFT+U spin-splitting
  logical           :: pub_hubbard_restart   !  If restarting a HUBBARDSCF calculation (not input)
  logical           :: pub_hubbard_atomsolve ! If using atomic guesses as DFT+U projectors
  logical           :: pub_hub_calculating_u ! If calculating a bare static response (fixed Hamiltonian) (not input)

!CW
  real(8)           :: pub_dmft_doping,pub_dmft_scaling_cutoff,pub_dmft_scaling_cutoff_h,pub_dmft_scaling_tol,pub_dmft_scaling_tail
  integer           :: pub_dmft_scaling_maxspace, pub_dmft_scaling_iter,pub_dmft_scaling_nmpi,pub_dmft_scaling_meth
  logical           :: pub_dmft_impose_same_coeffs,pub_dmft_impose_chem_spin
  logical           :: pub_dmft_split,pub_dmft_splitk,pub_dmft_nbo,pub_dmft_kpoints_kernel_gamma,pub_dmft_lin_scaling
  logical           :: pub_dmft_invert_overlap,pub_dmft_local_scratch
  real(8)           :: pub_dmft_win
  integer           :: pub_dmft_nval
  logical           :: pub_dmft_kpoints_sym,pub_dmft_switch_off_proj_order,pub_dmft_fully_sc_H,pub_dmft_KS_shift,pub_dmft_purify_sc,pub_dmft_spoil_kernel,pub_dmft_fully_sc
  logical           :: pub_write_polarisation_plot
  real(kind=DP)     :: pub_dmft_mu_diff_max,pub_dmft_kernel_mix,pub_dmft_optics_window,pub_dmft_optics_x1,pub_dmft_optics_y1,pub_dmft_optics_z1,pub_dmft_order_proj
  integer           :: pub_dmft_nkpoints,pub_dmft_mu_order,pub_dmft_embed_iter,pub_dmft_optics_i1,pub_dmft_nmu_loop
  real(8)           :: pub_dmft_embed_mix
  integer           :: pub_dmft_optics_i2
  integer           :: pub_dmft_points
  logical           :: pub_dmft_skip_energy,pub_dmft_in_bohr,pub_dmft_plot_all_proj,pub_dmft_integrate_green,pub_dmft_ignore_type
  logical           :: pub_dmft_sc,pub_dmft_paramagnetic,pub_dmft_rotate_green,pub_dmft_plot_real_space,pub_dmft_optics,pub_dmft_plot_real_space_sigma
  integer           :: pub_dmft_norm_proj,pub_dmft_kernel
  real(kind=DP)     :: pub_dmft_smear,pub_dmft_smear_T,pub_dmft_smear_shift,pub_dmft_smear_w,pub_dmft_smear_eta
  real(kind=DP)     :: pub_dmft_temp
  real(kind=DP)     :: pub_dmft_chem_shift,pub_dmft_e_plot
  real(kind=DP)     :: pub_dmft_emin,pub_dmft_dos_min
  real(kind=DP)     :: pub_dmft_emax,pub_dmft_dos_max
  real(kind=DP)     :: pub_dmft_cutoff_tail,pub_dmft_cutoff_small,pub_dmft_free_green_frequ
!END CW

  ! gibo: parameters for Constrained_DFT (cDFT) with self-consistent projectors
  logical           :: pub_cdft                          ! Constrained DFT run (if .true.)
  logical           :: pub_ci_cdft                       ! Configuration-Interaction cDFT run (if .true.)
  logical           :: pub_cdft_hubbard                  ! Constrained DFT+U run (if .true.)
  logical           :: pub_cdft_atom_charge              ! ATOM-CHARGE-constrained CDFT simulation (if .true.)
  logical           :: pub_cdft_atom_spin                ! ATOM-SPIN-constrained CDFT simulation (if .true.)
  logical           :: pub_cdft_group_charge             ! GROUP-CHARGE-constrained CDFT simulation (if .true.)
  logical           :: pub_cdft_group_spin               ! GROUP-SPIN-constrained CDFT simulation (if .true.)
  logical           :: pub_cdft_group_charge_diff        ! GROUP-CHARGE-DIFFERENCE-constrained CDFT (if .true.)
  logical           :: pub_cdft_group_spin_diff          ! GROUP-SPIN-DIFFERENCE-constrained CDFT (if .true.)

  integer           :: pub_cdft_type                     ! integer to CASE-SELECT deal with different cdft-modes
  integer           :: maxit_cdft_u_cg                   ! maximum number of cDFT-U CG iterations
  real(kind=DP)     :: pub_cdft_max_grad     ! Maximum permissable cDFT U-gradient

  ! Targeted group-CHARGE for GROUP-CHARGE-constrained cDFT
  real(kind=DP)     :: pub_cdft_group_charge_target  
  ! Targeted group-SPIN for GROUP-SPIN-constrained cDFT
  real(kind=DP)     :: pub_cdft_group_spin_target  
  ! Targeted (acceptor-donor) CHARGE difference for  GROUP-CHARGE-DIFFERENCE-constrained cDFT
  real(kind=DP)     :: pub_cdft_group_charge_diff_target  
  ! Targeted (acceptor-donor) SPIN difference for  GROUP-SPIN-DIFFERENCE-constrained cDFT
  real(kind=DP)     :: pub_cdft_group_spin_diff_target  
  ! Constraining potential (eV) for GROUP-CHARGE/SPIN-DIFFERENCE-constrained CDFT simulation
  real(kind=DP)     :: pub_cdft_group_diff_u
  ! Constraining potential (eV) for GROUP-CHARGE/SPIN-constrained CDFT simulation
  real(kind=DP)     :: pub_cdft_group_u

  real(kind=DP)     :: cdft_cg_threshold      ! U-opt RMS gradient convergence threshold
  character(len=80) :: cdft_cg_type           ! type for CG CDFT U-optimisation
  integer           :: pub_cdft_cg_max        ! max number iterations before U-opt CG reset
  ! Max energy/atom change per cDFT iteration
  real(kind=DP)     :: pub_cdft_elec_energy_tol
  real(kind=DP)     :: cdft_cg_max_step       ! Max trial-step in CG U-opt
  logical           :: pub_cdft_guru          ! The user is a cDFT_guru, let her/him free to mess around...
  logical           :: pub_cdft_continuation  ! Continuate a previous cDFT-optimisation

  ! ndmh: parameters for Correction of PBCs
  character(len=80) :: pub_coulomb_cutoff_type
  real(kind=DP)     :: pub_coulomb_radius
  real(kind=DP)     :: pub_coulomb_length
  logical           :: pub_coulomb_cutoff
  logical           :: pub_coulomb_cutoff_write_int
  real(kind=DP)     :: pub_mt_cutoff

  ! ndmh: Properties calculation parameters
  logical           :: pub_do_properties
  integer           :: pub_num_eigenvalues
  integer           :: pub_homo_dens_plot
  integer           :: pub_lumo_dens_plot
  integer           :: pub_homo_plot
  integer           :: pub_lumo_plot
  real(kind=DP)     :: pub_dos_smear
  real(kind=DP)     :: pub_ldos_smear
  integer           :: pub_ldos_ngroups
  logical           :: pub_popn_calculate
  real(kind=DP)     :: pub_popn_bond_cutoff
  logical           :: pub_ngwf_analysis
  logical           :: pub_polarisation_calculate
  logical           :: pub_spread_calculate
  integer, allocatable :: pub_ldos_group_nsp(:) ! not input
  character(len=4), allocatable :: pub_ldos_groups(:,:)

  ! aam: geometry optimiser parameters
  real(kind=DP)     :: geom_modulus_est
  real(kind=DP)     :: geom_frequency_est
  real(kind=DP)     :: geom_energy_tol
  real(kind=DP)     :: geom_force_tol
  real(kind=DP)     :: geom_disp_tol
  integer           :: geom_max_iter
  integer           :: geom_convergence_win
  integer           :: geom_backup_iter
  integer           :: geom_reset_dk_ngwfs_iter
  integer           :: geom_lbfgs_max_updates
  integer           :: geom_lbfgs_block_length
  logical           :: geom_continuation
  logical           :: geom_print_inv_hessian
  logical           :: pub_geom_reuse_dk_ngwfs
  character(len=80) :: geom_method
  logical           :: geom_lbfgs


  ! vm: transition state search parameters
  character(len=8)  :: tssearch_method          ! transition state search method (LSTQST for now)
  character(len=20) :: tssearch_lstqst_protocol ! which protocol to use
  integer           :: tssearch_qst_max_iter    !
  integer           :: tssearch_cg_max_iter     !
  real (kind=dp)    :: tssearch_force_tol       ! force tolerance
  real (kind=dp)    :: tssearch_disp_tol        ! displacement tolerance

  ! aam: molecular dymanics
  real(kind=dp)     :: md_delta_t
  integer           :: md_num_iter
  logical           :: md_init_velocities
  logical           :: md_restart
  integer           :: md_reset_dkn_ngwfs
  logical           :: pub_md_properties
  logical           :: mts_xi
  integer           :: mts_nstep
  real(kind=dp)     :: mts_delta_t
  real(kind=dp)     :: mts_ngwf_threshold
  real(kind=dp)     :: mts_lnv_threshold
  integer           :: mix_dkn_num
  integer           :: mix_ngwfs_num
  real(kind=dp)     :: mix_ngwfs_coeff
  integer           :: mix_dkn_type
  integer           :: mix_ngwfs_type
  real(kind=dp)     :: mix_local_length
  real(kind=dp)     :: mix_local_smear

  ! lr408: Conduction parameters
  logical           :: pub_cond_calculate
  logical           :: cond_read_denskern
  logical           :: cond_read_tightbox_ngwfs
  logical           :: cond_fixed_shift
  logical           :: cond_calc_max_eigen
  integer           :: cond_num_states
  real(kind=DP)     :: cond_kernel_cutoff
  real(kind=DP)     :: cond_init_shift
  real(kind=DP)     :: cond_shift_buffer
  integer           :: cond_num_extra_states
  integer           :: cond_num_extra_its
  logical           :: cond_plot_joint_orbitals
  logical           :: cond_plot_vc_orbitals

  ! smmd: Transmission coefficients
  logical           :: pub_etrans_calculate
  logical           :: pub_etrans_same_leads
  logical           :: pub_etrans_bulk
  logical           :: pub_etrans_write_setup
  real(kind=dp)     :: pub_etrans_ecmplx

  real(kind=dp)     :: pub_etrans_emax
  real(kind=dp)     :: pub_etrans_emin
  integer           :: pub_etrans_enum
  integer           :: pub_etrans_source
  integer           :: pub_etrans_drain

  ! lr408: Conduction / optical spectra parameters
  logical :: pub_spectra_calculate
  ! lr408: if true calculate momentum matrix elements for spectra,
  ! lr408: otherwise calculate position matrix elements  (for molecules)
  logical :: pub_calc_mom_mat_els
  logical :: pub_calc_nonloc_comm
  logical :: pub_spec_cont_deriv
  real(kind=dp) :: pub_spec_nonloc_fin_diff
  logical :: pub_spectra_print_mat_els
  real(kind=dp) :: pub_scissor
  real(kind=dp) :: pub_opt_smear

  logical :: pub_eels_calculate

  ! fc: phonon parameters
  logical :: pub_have_phonon_disp_list
  integer :: pub_phonon_farming_task, pub_num_disp
  integer, allocatable :: pub_phonon_disp_list(:)
  real(kind=DP) :: pub_phonon_disp, pub_phonon_fmax
  real(kind=DP) :: pub_phonon_tmin, pub_phonon_tmax
  real(kind=DP) :: pub_phonon_deltat, pub_phonon_min_freq

  ! pdh: bandstructure parameters
  integer           :: pub_bs_kpoint_path_length
  integer           :: pub_bs_num_eigenvalues
  integer           :: pub_bs_unfold(3)
  logical           :: pub_do_bandstructure
  character(len=80) :: pub_bs_method
  real(kind=DP)     :: pub_bs_kpoint_path_spacing
  real(kind=DP), allocatable :: pub_bs_kpoint_path_start(:,:)
  real(kind=DP), allocatable :: pub_bs_kpoint_path_end(:,:)

  ! qoh: Hartree-Fock exchange parameters:
  ! qoh: Whether Hartree-Fock exchange is being used, set in xc module
  logical           :: pub_usehfx ! (not input)
  logical           :: pub_hfxsw  ! Are we using spherical waves for hf exchange
  ! jd: Variant of integration routine to use
  integer           :: pub_hfx_integration_variant
  ! qoh: Use non local pseudopotential sparsity for exchange matrix
  logical           :: pub_hfx_nlpp_for_exchange
  ! jd: Number of segments for radial integration
  integer           :: pub_hfx_radial_segments, pub_hfx_angular_segments
  ! jd: Details of SW expansion
  integer           :: pub_hfx_max_l, pub_hfx_max_zeros
  ! jd: Details of Chebyshev expansion of spherical waves
  integer           :: pub_hfx_cheb_order, pub_hfx_cheb_intervals
  ! jd: Batch sizes for products of Chebyshev expansions
  integer           :: pub_hfx_cheb_a_batchsize, pub_hfx_cheb_b_batchsize
  ! jd: V matrix read/write for restarts
  logical           :: pub_hfx_read_vmatrix, pub_hfx_write_vmatrix
  ! jd: X matrix read/write for restarts
  logical           :: pub_hfx_read_xmatrix, pub_hfx_write_xmatrix
  ! jd: HFx debug
  logical           :: pub_hfx_debug
  ! qoh: Whether we are using tightboxes to do FFTs
  logical           :: pub_tightbox_fft_coarse ! (not input)
  logical           :: pub_tightbox_fft_fine   ! (not input)

  ! ddor: TDDFT flags and parameters
  logical           :: pub_do_tddft
  real(kind=DP)     :: pub_tddft_maximum_energy
  real(kind=DP)     :: pub_tddft_resolution
  character(len=20) :: pub_tddft_propagation_method
  integer           :: pub_tddft_sparsity_level
  logical           :: pub_tddft_tammdancoff
  real(kind=DP)     :: pub_tddft_dipole_kick_strength(3)
  character(len=80) :: tddft_xc_functional
  real(kind=DP)     :: pub_tddft_hamiltonian_mixing
  real(kind=DP)     :: pub_tddft_damping
  logical           :: pub_tddft_enforced_idempotency
  integer           :: pub_tddft_maxit_hotelling
  real(kind=DP)     :: pub_tddft_max_resid_hotelling
  logical           :: pub_tddft_inv_overlap_exact

  ! jd: for implicit solvent
  real(kind=DP)     :: pub_is_density_threshold
  real(kind=DP)     :: pub_is_solvation_beta
  real(kind=DP)     :: pub_is_bulk_permittivity
  real(kind=DP)     :: pub_is_multigrid_defect_error_t
  real(kind=DP)     :: pub_is_multigrid_error_tol
  integer           :: pub_is_multigrid_max_iters
  integer           :: pub_is_multigrid_nlevels
  integer           :: pub_is_discretization_order
  integer           :: pub_is_bc_coarseness
  integer           :: pub_is_bc_surface_coarseness
  real(kind=DP)     :: pub_is_bc_threshold
  real(kind=DP)     :: pub_is_smeared_ion_width
  real(kind=DP)     :: pub_is_surface_thickness
  real(kind=DP)     :: pub_is_solvent_surface_tension
  logical           :: pub_is_implicit_solvent
  logical           :: pub_is_check_solv_energy_grad
  logical           :: pub_is_smeared_ion_rep
  logical           :: pub_is_include_cavitation
  character(len=80) :: pub_is_solvation_method
  character(len=80) :: pub_is_dielectric_model
  character(len=80) :: pub_is_dielectric_function
  character(len=80) :: pub_is_solvation_output_detail
  real(kind=DP)     :: is_default_bulk_permittivity ! jd: not used outside here

  ! jd: for open BCs in local pseudo and ion-ion energies
  !   - is ion-ion energy calculated with open BC (by direct summation)
  logical           :: pub_ii_energy_direct = .false.
  !   - is hartree calculated with MG? (currently this implies open BC)
  logical           :: pub_multigrid_hartree = .false.
  !   - is multgrid used? (affects cell_grid_distribute)
  logical           :: pub_multigrid_in_use = .false.
  !   - is local pseudo calculated with open BC?
  logical           :: pub_open_localpseudo = .false.
  !   - did the user force open BCs in ion-ion energy?
  logical           :: pub_openbc_ion_ion
  !   - did the user force open BCs in hartree calculation?
  logical           :: pub_openbc_hartree
  !   - did the user force open BCs in local pseudo calculation?
  logical           :: pub_openbc_pspot
  !   - parameter 'npts_x' in the open BC local pseudo calculation
  integer           :: pub_openbc_pspot_finetune_nptsx
  !   - parameter 'alpha' in the open BC local pseudo calculation
  real(kind=DP)     :: pub_openbc_pspot_finetune_alpha
  !   - parameter 'fineness' in the open BC local pseudo calculation
  integer           :: pub_openbc_pspot_finetune_f

  ! lpl: NBO stuff (part of properties_calculate)
  logical :: pub_write_nbo
  logical :: pub_nbo_init_lclowdin
  logical :: pub_nbo_write_lclowdin
  logical :: pub_nbo_write_npacomp
  logical :: pub_nbo_scale_dm

  logical, allocatable :: pub_nbo_write_species(:)      ! not input
  character(len=256), allocatable :: pub_nbo_ngwf_label(:) ! not input

  logical :: pub_nbo_pnao_analysis
  character(len=80) :: pub_nbo_aopnao_scheme

  logical :: pub_plot_nbo
  character(len=8) :: pub_nbo_plotorbtype
  integer, allocatable :: pub_nbo_list_plotnbo(:) ! not input
 ! lpl: NBO stuff (part of properties_calculate)


contains

  subroutine get_rundat

    use comms, only: comms_abort, comms_bcast, pub_on_root, pub_root_node_id, &
         pub_comms_group_size, pub_total_num_nodes
    use constants, only: DP, stdout, ANGSTROM, NORMAL, BRIEF, VERBOSE
    use esdf, only: esdf_double, esdf_string, esdf_integer, esdf_boolean,&
         esdf_block, block_data, esdf_physical, esdf_convfac
    use utils, only: utils_alloc_check, utils_assert

    character(len=80) :: buf, tdbuf ! ddor: tdbuf for TDDFT
    character(len=80) :: output_detail
    character(len=80) :: paw_output_detail
    character(len=40) :: psinc_spacing_unit
    integer :: ierr
    integer :: nlines,iline,ipath
    logical :: bs_start
    logical :: inconsistent_bcs
    real(kind=DP) :: dummy_real(3)

    ! Root node reads from input file...

    if (pub_on_root) then
       task                    = esdf_string('task', 'SINGLEPOINT')
       cutoff_energy           = esdf_physical('cutoff_energy', 20.0_dp, 'hartree')
       kernel_cutoff           = esdf_double('kernel_cutoff',1000.0_dp)
       xc_functional           = esdf_string('xc_functional','LDA')
       pub_libxc_x_func_id     = esdf_integer('libxc_x_func_id', 0)
       pub_libxc_c_func_id     = esdf_integer('libxc_c_func_id', 0)
       pub_charge              = esdf_double('charge',0.0_DP)
       pub_spin                = esdf_integer('spin', 0)
       pub_spin_polarised      = esdf_boolean('spin_polarized',.false.)
       pub_spin_polarised      = esdf_boolean('spin_polarised', &
            pub_spin_polarised)
       buf                     = esdf_string('constant_efield','0.0 0.0 0.0')
       read(buf,*) pub_constant_efield(:)
       buf                     = esdf_string('fftbox_pref','0 0 0')
       read(buf,*) fftbox_pref(:)
       buf                     = esdf_string('ppd_npoints', '0 0 1')
       read(buf,*) ppd_npoints(:)
       buf                     = esdf_string('psinc_spacing','0.0 0.0 0.0')
       ! ndmh: append "bohr" to psinc_spacing as default units
       buf = trim(adjustl(buf))//' bohr'
       read(buf,*) psinc_spacing(:),psinc_spacing_unit
       psinc_spacing(:) = psinc_spacing(:) * &
            esdf_convfac(psinc_spacing_unit,'bohr')
       pub_old_input           = esdf_boolean('old_input',.false.)
       pub_devel_code          = esdf_string('devel_code','')
       pub_nnho                = esdf_boolean('nnho',.false.)
       pub_ngwf_halo           = esdf_double('ngwf_halo', -1.0_DP)
       pub_nonsc_forces        = esdf_boolean('nonsc_forces',.false.)
       pub_external_pressure   = esdf_physical('external_pressure',0.0_DP, &
                   'ha/bohr**3')
       pub_smoothing_factor    = esdf_double('smoothing_factor',5.0_DP)            
       pub_isosurface_cutoff   = esdf_double('isosurface_cutoff',0.0005_DP) 
       pub_maxit_palser_mano   = esdf_integer('maxit_palser_mano', 50)
       maxit_pen               = esdf_integer('maxit_pen', 3)
       pen_param               = esdf_double('pen_param', 4.0_DP)
       minit_lnv               = esdf_integer('minit_lnv',3)
       maxit_lnv               = esdf_integer('maxit_lnv',8)
       maxit_lnv_coarse        = esdf_integer('maxit_lnv_coarse',0)
       lnv_threshold_orig      = esdf_double('lnv_threshold_orig',1.0e-9_dp)
       lnv_cg_type             = esdf_string('lnv_cg_type','LNV_FLETCHER')
       lnv_cg_max_step         = esdf_double('lnv_cg_max_step', 3.00_DP)
       exact_lnv               = esdf_boolean('exact_lnv',.true.)
       old_lnv                 = esdf_boolean('old_lnv',.false.)
       pub_lnv_check_trial_steps = esdf_boolean('lnv_check_trial_steps', &
            .false.)
       pub_kerfix              = esdf_integer('kerfix', 1)
       pub_maxit_kernel_fix    = esdf_integer('maxit_kernel_fix', 3)
       pub_initial_dens_realspace = esdf_boolean('initial_dens_realspace', &
            .false.)
       pub_kernel_diis         = esdf_boolean('kernel_diis',.false.)
       buf = esdf_string('kernel_diis_type', 'P')
       pub_kernel_diis_type    = buf(1:1)
       pub_kernel_diis_max     = esdf_integer('kernel_diis_max', 10)
       pub_kernel_diis_maxit   = esdf_integer('kernel_diis_maxit', 25)
       pub_kernel_diis_thres   = esdf_double('kernel_diis_thres', 1.0e-9_DP)
       pub_kernel_diis_liter   = esdf_integer('kernel_diis_liter', 5)
       pub_kernel_diis_c_in    = esdf_double('kernel_diis_c_in', 0.9_DP)
       pub_kernel_diis_c_out   = 1.0_DP - pub_kernel_diis_c_in
       buf = esdf_string('kernel_diis_criteria', '1000')
       pub_kernel_diis_criteria = buf(1:4)

       maxit_ngwf_cg           = esdf_integer('maxit_ngwf_cg',60)
       ngwf_cg_type            = esdf_string('ngwf_cg_type','NGWF_FLETCHER')
       ngwf_cg_max_step        = esdf_double('ngwf_cg_max_step',-8.00_DP)
       pub_elec_cg_max         = esdf_integer('elec_cg_max', 5)
       precond_scheme          = esdf_string('precond_scheme','TETER')
       k_zero                  = esdf_double('k_zero', 3.0_dp)
       r_precond               = esdf_double('r_precond', 2.0_dp)
       precond_recip           = esdf_boolean('precond_recip',.true.)
       precond_real            = esdf_boolean('precond_real',.false.)
       smooth_scheme           = esdf_string('smooth_scheme','NONE')
       r_smooth                = esdf_physical('r_smooth',1.5_dp,'bohr')
       k_smooth                = esdf_double('k_smooth',5.0_dp)
       occ_mix                 = esdf_double('occ_mix', 0.25_DP)
       pub_kernel_update       = esdf_boolean('kernel_update', .false.)
       pub_kernel_christoffel_update = &
            &esdf_boolean('kernel_christoffel_update',.false.)
       maxit_hotelling         = esdf_integer('maxit_hotelling', 50)
       max_resid_hotelling     = esdf_double('max_resid_hotelling', 1.0e-12_DP)
       maxit_ngwf_diis         = esdf_integer('maxit_ngwf_diis',0)
       ngwf_diis_size          = esdf_integer('ngwf_diis_size',4)
       ngwf_diis_reset         = esdf_integer('ngwf_diis_reset',8)
       ngwf_diis_max_step      = esdf_double('ngwf_diis_max_step', 2.00_DP)

       ngwf_threshold_orig     = esdf_double('ngwf_threshold_orig',2.0e-6_dp)
       pub_elec_energy_tol     = esdf_physical('elec_energy_tol',-0.001_DP, &
            'hartree')
       pub_elec_force_tol      = esdf_physical('elec_force_tol',-0.001_DP, &
            'ha/bohr')
       pub_ngwf_max_grad       = esdf_double('ngwf_max_grad',-2.0e-5_dp)
       delta_e_conv            = esdf_boolean('delta_e_conv', .true.)

       locpot_int_batch_size   = esdf_integer('locpot_int_batch_size', 10)
       kinetic_int_batch_size  = esdf_integer('kinetic_int_batch_size', 10)
       density_batch_size      = esdf_integer('density_batch_size', 10)
       ngwf_grad_batch_size    = esdf_integer('ngwf_grad_batch_size', 10)
       pub_dense_threshold     = esdf_double('dense_threshold',0.40_DP)
       pub_comms_group_size    = esdf_integer('comms_group_size',4)
       pub_ovlp_for_nonlocal   = esdf_boolean('ovlp_for_nonlocal',.false.)
       use_space_filling_curve = esdf_boolean('use_space_filling_curve', .true.)
       coreham_denskern_guess  = esdf_boolean('coreham_denskern_guess', .true.)
       pub_check_atoms         = esdf_boolean('check_atoms', .true.)
       pub_locpot_scheme       = esdf_string('locpot_scheme','FULL')
       pub_smooth_projectors   = esdf_double('smooth_projectors', -0.4_DP)
       pub_smooth_loc_pspot    = esdf_double('smooth_loc_pspot', -0.4_DP)
       pub_odd_psinc_grid      = esdf_boolean('odd_psinc_grid',.false.)
       pub_even_psinc_grid      = esdf_boolean('even_psinc_grid',.false.)
       pub_realspace_projectors= esdf_boolean('realspace_projectors',.false.)

       output_detail           = esdf_string('output_detail','NORMAL')
       paw_output_detail       = esdf_string('paw_output_detail','DEFAULT')
       timings_level           = esdf_integer('timings_level', 1)
       pub_write_params        = esdf_boolean('write_params', .true.)
       pub_write_forces        = esdf_boolean('write_forces',.false.)
       pub_write_positions     = esdf_boolean('write_positions',.true.)
       pub_write_xyz           = esdf_boolean('write_xyz', .false.)
       print_qc                = esdf_boolean('print_qc',.false.)
#ifdef ACCELRYS
       pub_cube_format         = esdf_boolean('cube_format', .false.)
#else
       pub_cube_format         = esdf_boolean('cube_format', .true.)
#endif
       pub_dx_format           = esdf_boolean('dx_format', .false.)
       pub_dx_coarse           = esdf_boolean('dx_format_coarse', .false.)
       pub_dx_sig_digits       = esdf_integer('dx_format_digits', 7)
#ifdef ACCELRYS
       pub_grd_format          = esdf_boolean('grd_format', .true.)
#else
       pub_grd_format          = esdf_boolean('grd_format', .false.)
#endif
       pub_write_density_plot  = esdf_boolean('write_density_plot', .true.)
       pub_write_ngwf_plot     = esdf_boolean('write_ngwf_plot',.false.)
       pub_write_ngwf_grad_plot = esdf_boolean('write_ngwf_grad_plot',.false.)
       pub_write_ngwf_radial    = esdf_integer('write_ngwf_radial',0)
       pub_write_ngwf_grad_radial = esdf_integer('write_ngwf_grad_radial',0)
       pub_write_radial_step    = esdf_physical('write_radial_step',0.005_dp,'bohr')
       pub_write_radial_smear   = esdf_physical('write_radial_smear',0.01_dp,'bohr')

!CW
       pub_dmft_mu_diff_max            = esdf_double('dmft_mu_diff_max',0.00d0)
       pub_write_polarisation_plot     = esdf_boolean('write_polarisation_plot',.false.)
       pub_dmft_ignore_type            = esdf_boolean('dmft_ignore_type',.false.)
       pub_dmft_KS_shift               = esdf_boolean('dmft_KS_shift',.true.)
       pub_dmft_purify_sc              = esdf_boolean('dmft_purify_sc',.false.)
       pub_dmft_kernel_mix             = esdf_double('dmft_kernel_mix',0.1d0)
       pub_dmft_switch_off_proj_order  = esdf_boolean('dmft_switch_off_proj_order',.false.)
       pub_dmft_fully_sc_h             = esdf_boolean('dmft_fully_sc_h',.false.)
       pub_dmft_fully_sc               = esdf_boolean('dmft_fully_sc',.false.)
       pub_dmft_spoil_kernel           = esdf_boolean('dmft_spoil_kernel',.false.)
       pub_dmft_kernel                 = esdf_integer('dmft_kernel',0)
       pub_dmft_sc                     = esdf_boolean('dmft_sc',.false.)
       pub_dmft_order_proj             = esdf_double('dmft_order_proj',0.00_dp)
       pub_dmft_in_bohr                = esdf_boolean('dmft_in_bohr',.false.)
       pub_dmft_temp                   = esdf_physical('dmft_temp',-0.01_dp,'hartree')
       pub_dmft_doping                 = esdf_double('dmft_doping',0.00d0)
       pub_dmft_smear                  = esdf_physical('dmft_smear',0.00018_dp,'hartree')
       pub_dmft_smear_T                = esdf_physical('dmft_smear_T',0.008_dp,'hartree')
       pub_dmft_smear_shift            = esdf_physical('dmft_smear_shift',0.000_dp,'hartree')
       pub_dmft_smear_eta              = esdf_physical('dmft_smear_eta',0.01_dp,'hartree')
       pub_dmft_smear_w                = esdf_physical('dmft_smear_w',0.035_dp,'hartree')
       pub_dmft_points                 = esdf_integer('dmft_points',0)
       pub_dmft_plot_real_space        = esdf_boolean('dmft_plot_real_space',.false.)
       pub_dmft_plot_real_space_sigma  = esdf_boolean('dmft_plot_real_space_sigma',.false.)
       pub_dmft_optics_window          = esdf_physical('dmft_optics_window',0.10_dp,'hartree')
       pub_dmft_optics_x1              = esdf_double('dmft_optics_x1',0.00_dp)
       pub_dmft_optics_y1              = esdf_double('dmft_optics_y1',0.00_dp)
       pub_dmft_optics_z1              = esdf_double('dmft_optics_z1',0.00_dp)
       pub_dmft_nmu_loop               = esdf_integer('dmft_nmu_loop',1)
       pub_dmft_optics_i1              = esdf_integer('dmft_optics_i1',1)
       pub_dmft_embed_iter             = esdf_integer('dmft_embed_iter',10)
       pub_dmft_nkpoints               = esdf_integer('dmft_nkpoints',1)
       pub_dmft_kpoints_sym            = esdf_boolean('dmft_kpoints_sym',.false.)
       pub_dmft_kpoints_kernel_gamma   = esdf_boolean('dmft_kpoints_kernel_gamma',.false.)
       pub_dmft_lin_scaling            = esdf_boolean('dmft_lin_scaling',.false.)
       pub_dmft_nval                   = esdf_integer('dmft_nval',40)
       pub_dmft_win                    = esdf_double('dmft_win',0.3d0)
       pub_dmft_scaling_cutoff         = esdf_double('dmft_scaling_cutoff',1.d-8)
       pub_dmft_scaling_cutoff_h       = esdf_double('dmft_scaling_cutoff_h',1.d-5)
       pub_dmft_scaling_tail           = esdf_double('dmft_scaling_tail',10.d0)
       pub_dmft_scaling_tol            = esdf_double('dmft_scaling_tol',1.d-8)
       pub_dmft_scaling_maxspace       = esdf_integer('dmft_scaling_maxspace',20)
       pub_dmft_scaling_iter           = esdf_integer('dmft_scaling_iter',2000)
       pub_dmft_scaling_nmpi           = esdf_integer('dmft_scaling_nmpi',1)
       pub_dmft_scaling_meth           = esdf_integer('dmft_scaling_meth',4)
       pub_dmft_split                  = esdf_boolean('dmft_split',.false.)
       pub_dmft_splitk                 = esdf_boolean('dmft_splitk',.false.)
       pub_dmft_nbo                    = esdf_boolean('dmft_nbo',.false.)
       pub_dmft_impose_same_coeffs     = esdf_boolean('dmft_impose_same_coeffs',.false.)
       pub_dmft_impose_chem_spin       = esdf_boolean('dmft_impose_chem_spin',.false.)
       pub_dmft_invert_overlap         = esdf_boolean('dmft_invert_overlap',.false.)
       pub_dmft_local_scratch          = esdf_boolean('dmft_local_scratch',.false.)
       pub_dmft_mu_order               = esdf_integer('dmft_mu_order',2)
       pub_dmft_embed_mix              = esdf_double('dmft_embed_mix',0.7d0)
       pub_dmft_optics_i2              = esdf_integer('dmft_optics_i2',1)
       pub_dmft_skip_energy            = esdf_boolean('dmft_skip_energy',.false.)
       pub_dmft_optics                 = esdf_boolean('dmft_optics',.false.)
       pub_dmft_plot_all_proj          = esdf_boolean('dmft_plot_all_proj',.false.)
       pub_dmft_integrate_green        = esdf_boolean('dmft_integrate_green',.false.)
       pub_dmft_paramagnetic           = esdf_boolean('dmft_paramagnetic',.false.)
       pub_dmft_rotate_green           = esdf_boolean('dmft_rotate_green',.false.)
       pub_dmft_norm_proj              = esdf_integer('dmft_norm_proj',0)
       pub_dmft_chem_shift             = esdf_physical('dmft_chem_shift',0.0_dp,'hartree')
       pub_dmft_e_plot                 = esdf_physical('dmft_e_plot',0.0_dp,'hartree')
       pub_dmft_emin                   = esdf_physical('dmft_emin',-1.0_dp,'hartree')
       pub_dmft_emax                   = esdf_physical('dmft_emax', 1.0_dp,'hartree')
       pub_dmft_dos_min                = esdf_physical('dmft_dos_min',-10.0_dp,'hartree')
       pub_dmft_dos_max                = esdf_physical('dmft_dos_max', 10.0_dp,'hartree')
       pub_dmft_cutoff_tail            = esdf_physical('dmft_cutoff_tail',10.0_dp,'hartree')
       pub_dmft_cutoff_small           = esdf_physical('dmft_cutoff_small',0.0_dp,'hartree')
       pub_dmft_free_green_frequ       = esdf_physical('dmft_free_green_frequ',200.0_dp,'hartree')
!END CW

       read_denskern           = esdf_boolean('read_denskern', .false.)
       write_denskern          = esdf_boolean('write_denskern', .true.)
       read_tightbox_ngwfs     = esdf_boolean('read_tightbox_ngwfs', .false.)
       write_tightbox_ngwfs    = esdf_boolean('write_tightbox_ngwfs', .true.)
       pub_write_converged_dk_ngwfs = esdf_boolean('write_converged_dk_ngwfs', .false.)
       pub_read_sw_ngwfs       = esdf_boolean('read_sw_ngwfs', .false.)
       pub_write_sw_ngwfs      = esdf_boolean('write_sw_ngwfs', .false.)
       pub_write_max_l         = esdf_integer('write_max_l', 3)
       pub_read_max_l          = esdf_integer('read_max_l', 3)
       pub_extra_n_sw          = esdf_integer('extra_n_sw', 0)

       pub_paw                 = esdf_boolean('paw', .false.)
       pub_fine_grid_scale     = esdf_double('fine_grid_scale',2.0_DP)
       pub_dbl_grid_scale      = esdf_double('dbl_grid_scale',2.0_DP)
       pub_aug_funcs_recip     = esdf_boolean('aug_funcs_recip',.true.)

       pub_dispersion          = esdf_integer('dispersion', 0)
       pub_vdw_dcoeff          = esdf_double('vdw_dcoeff',-1.0_DP)

       pub_hub_max_iter        = esdf_integer('hubbard_max_iter', 0)
       pub_hub_energy_tol      = esdf_physical('hubbard_energy_tol', 1.0e-8_dp,&
            'hartree')
       pub_hub_conv_win        = esdf_integer('hubbard_conv_win', 2)
       pub_hub_proj_mixing     = esdf_double('hubbard_proj_mixing',0.0_DP)
       pub_hub_functional      = esdf_integer('hubbard_functional',1)
       pub_hub_tensor_corr     = esdf_integer('hubbard_tensor_corr',1)
       pub_hub_ngwf_spin_thr   = esdf_double('hubbard_ngwf_spin_threshold', &
            2.0e-5_dp)

       pub_cdft_atom_charge       = esdf_boolean('cdft_atom_charge',.false.)
       pub_cdft_atom_spin         = esdf_boolean('cdft_atom_spin',.false.)
       pub_cdft_group_charge      = esdf_boolean('cdft_group_charge',.false.)
       pub_cdft_group_spin        = esdf_boolean('cdft_group_spin',.false.)
       pub_cdft_group_charge_diff = esdf_boolean('cdft_group_charge_diff', &
            .false.)
       pub_cdft_group_spin_diff   = esdf_boolean('cdft_group_spin_diff', &
            .false.)
       pub_cdft_hubbard           = esdf_boolean('cdft_hubbard',.false.)
       pub_ci_cdft                = esdf_boolean('ci_cdft',.false.)

       ! set cdft_type internally and broadcast
       pub_cdft_type = 0
       if (pub_cdft_atom_charge)       pub_cdft_type = 1
       if (pub_cdft_atom_spin)         pub_cdft_type = 2
       if (pub_cdft_group_charge)      pub_cdft_type = 3
       if (pub_cdft_group_spin)        pub_cdft_type = 4
       if (pub_cdft_group_charge_diff) pub_cdft_type = 5
       if (pub_cdft_group_spin_diff)   pub_cdft_type = 6

       pub_cdft_group_charge_target      = esdf_double( &
            'cdft_group_charge_target',0._DP)
       pub_cdft_group_spin_target        = esdf_double( &
            'cdft_group_spin_target',0._DP)
       pub_cdft_group_charge_diff_target = esdf_double( &
            'cdft_group_charge_diff_target',0._DP)
       pub_cdft_group_spin_diff_target   = esdf_double( &
            'cdft_group_spin_diff_target',0._DP)
       pub_cdft_group_u                  = esdf_physical( &
            'cdft_group_u',0._DP,'hartree')
       pub_cdft_group_diff_u             = esdf_physical( &
            'cdft_group_diff_u',0._DP,'hartree')
       maxit_cdft_u_cg         = esdf_integer('maxit_cdft_u_cg',60)
       cdft_cg_type           = esdf_string('cdft_cg_type','NGWF_FLETCHER')
       cdft_cg_threshold      = esdf_double('cdft_cg_threshold',2.0e-6_dp)
       pub_cdft_cg_max         = esdf_integer('cdft_cg_max', 5)
       pub_cdft_max_grad       = esdf_double('cdft_max_grad', 2.0e-5_dp)
       pub_cdft_elec_energy_tol= esdf_physical('cdft_elec_energy_tol', -0.001_DP, &
            'hartree')
       cdft_cg_max_step        = esdf_double('cdft_cg_max_step',8.00_DP)
       pub_cdft_guru           = esdf_boolean('cdft_guru', .false.)
       pub_cdft_continuation   = esdf_boolean('cdft_continuation', .false.)

       pub_coulomb_radius      = esdf_physical('coulomb_cutoff_radius', &
            0.0_DP,'bohr')
       pub_coulomb_length      = esdf_physical('coulomb_cutoff_length', &
            0.0_DP,'bohr')
       pub_coulomb_cutoff_type = esdf_string('coulomb_cutoff_type', 'NONE')
       pub_coulomb_cutoff_write_int = esdf_boolean('coulomb_cutoff_write_int',&
            .false.)
       pub_mt_cutoff           = esdf_double('pbc_correction_cutoff', 0.0_DP)

       pub_do_properties       = esdf_boolean('do_properties', .false.)
       pub_num_eigenvalues     = esdf_integer('num_eigenvalues', 10)
       pub_homo_dens_plot      = esdf_integer('homo_dens_plot', -1)
       pub_lumo_dens_plot      = esdf_integer('lumo_dens_plot', -1)
       pub_homo_plot           = esdf_integer('homo_plot', 5)
       pub_lumo_plot           = esdf_integer('lumo_plot', 5)
       pub_dos_smear           = esdf_physical('dos_smear',0.1_DP / &
            HARTREE_IN_EVS, 'hartree')
       pub_ldos_smear          = esdf_physical('ldos_smear',0.1_DP / &
            HARTREE_IN_EVS, 'hartree')
       pub_popn_calculate      = esdf_boolean('popn_calculate', &
            pub_do_properties)
       pub_popn_bond_cutoff    = esdf_physical('popn_bond_cutoff',&
            3.0_DP*ANGSTROM,'bohr')
       pub_ngwf_analysis       = esdf_boolean('ngwf_analysis', .false.)
       pub_polarisation_calculate = esdf_boolean('polarisation_calculate',&
            .false.)
       pub_spread_calculate    = esdf_boolean('spread_calculate',.false.)

       geom_modulus_est        = esdf_physical('geom_modulus_est', 0.017_dp, 'ha/bohr**3')
       geom_frequency_est      = esdf_physical('geom_frequency_est', 0.0076_dp, 'hartree')
       geom_energy_tol         = esdf_physical('geom_energy_tol', 1.0e-6_dp, 'hartree')
       geom_force_tol          = esdf_physical('geom_force_tol', 0.002_dp, 'ha/bohr')
       geom_disp_tol           = esdf_physical('geom_disp_tol', 0.005_dp, 'bohr')
       geom_max_iter           = esdf_integer('geom_max_iter', 50)
       geom_convergence_win    = esdf_integer('geom_convergence_win', 2)
       geom_continuation       = esdf_boolean('geom_continuation', .false.)
       geom_backup_iter        = esdf_integer('geom_backup_iter', 1)
       geom_reset_dk_ngwfs_iter= esdf_integer('geom_reset_dk_ngwfs_iter',6)
       geom_lbfgs_max_updates  = esdf_integer('geom_lbfgs_max_updates',30)
       geom_lbfgs_block_length = esdf_integer('geom_lbfgs_block_length',30)
       pub_geom_reuse_dk_ngwfs = esdf_boolean('geom_reuse_dk_ngwfs', .true.)
       geom_print_inv_hessian  = esdf_boolean('geom_print_inv_hessian',.false.)
       geom_method             = esdf_string('geom_method', 'CARTESIAN')
       geom_lbfgs              = esdf_boolean('geom_lbfgs',.false.)

       ! Set up LBFGS block allocation size...
       if(geom_lbfgs_max_updates.ne.0) then
          geom_lbfgs_block_length = geom_lbfgs_max_updates
       end if

       buf                     = esdf_string('tssearch_method', 'LSTQST')
       tssearch_method         = buf(1:8)
       buf                     = esdf_string('tssearch_lstqst_protocol', &
            'LSTMAXIMUM')
       tssearch_lstqst_protocol= buf(1:20)
       tssearch_qst_max_iter   = esdf_integer('tssearch_qst_max_iter', 5)
       tssearch_cg_max_iter    = esdf_integer('tssearch_cg_max_iter', 20)
       tssearch_force_tol      = esdf_physical('tssearch_force_tol', 0.005_DP,&
            'ha/bohr')
       tssearch_disp_tol       = esdf_physical('tssearch_disp_tol', 0.01_DP, &
            'bohr')

       pub_md_properties       = esdf_boolean('md_properties',.false.)
       md_restart              = esdf_boolean('md_restart',.false.)
       md_delta_t              = esdf_physical('md_delta_t',40.0_DP,'aut') !~1fs
       md_num_iter             = esdf_integer('md_num_iter',100)
       md_reset_dkn_ngwfs      = esdf_integer('md_reset_denskern_ngwfs',5)
       mts_xi                  = esdf_boolean('mts_xi',.false.)
       mts_nstep               = esdf_integer('mts_nstep',1)
       mts_ngwf_threshold      = esdf_double('mts_ngwf_threshold',5.0e-4_DP)
       mts_lnv_threshold       = esdf_double('mts_ngwf_threshold',5.0e-6_DP)

       mix_dkn_num             = esdf_integer('mix_denskern_num',0)
       mix_ngwfs_num           = esdf_integer('mix_ngwfs_num',0)
       mix_dkn_type            = esdf_integer('mix_denskern_type',0)
       mix_ngwfs_coeff           = esdf_double('mix_ngwfs_coeff',1.0_dp)
       mix_ngwfs_type          = esdf_integer('mix_ngwfs_type',0)
       mix_local_length        = esdf_physical('mix_local_length',0.0_DP,'bohr')
       mix_local_smear         = esdf_physical('mix_local_smear',0.0_DP,'bohr')


       cond_read_denskern      = esdf_boolean('cond_read_denskern',.false.)
       cond_read_tightbox_ngwfs= esdf_boolean('cond_read_tightbox_ngwfs', &
            .false.)
       cond_fixed_shift        = esdf_boolean('cond_fixed_shift',.false.)
       cond_calc_max_eigen     = esdf_boolean('cond_calc_max_eigen',.true.)
       cond_num_states         = esdf_integer('cond_num_states', 0)
       cond_kernel_cutoff      = esdf_double('cond_kernel_cutoff',1000.0_DP)
       cond_init_shift         = esdf_double('cond_init_shift',0.0_DP)
       cond_shift_buffer       = esdf_double('cond_shift_buffer',0.1_DP)
       cond_num_extra_states   = esdf_integer('cond_num_extra_states', 0)
       cond_num_extra_its      = esdf_integer('cond_num_extra_its', 0)
       cond_plot_joint_orbitals = esdf_boolean('cond_plot_joint_orbitals',.true.)
       cond_plot_vc_orbitals    = esdf_boolean('cond_plot_vc_orbitals',.true.)

       pub_spectra_calculate    = esdf_boolean('cond_calc_optical_spectra',.false.)
       pub_calc_mom_mat_els     = esdf_boolean('cond_spec_calc_mom_mat_els',.true.)
       pub_calc_nonloc_comm     = esdf_boolean('cond_spec_calc_nonloc_comm',.true.)
       pub_spec_cont_deriv      = esdf_boolean('cond_spec_cont_deriv',.true.)
       pub_spec_nonloc_fin_diff = esdf_double('cond_spec_nonloc_comm_shift',0.0001_dp)
       pub_spectra_print_mat_els= esdf_boolean('cond_spec_print_mat_els',.true.)
       pub_opt_smear            = esdf_physical('cond_spec_opt_smear',0.1_DP / &
            HARTREE_IN_EVS, 'hartree')
       pub_scissor              = esdf_physical('cond_spec_sciss_op',0.0_DP / &
            HARTREE_IN_EVS, 'hartree')

       pub_eels_calculate       = esdf_boolean('cond_calc_eels',.false.)


       ! smmd: Transmission coefficients
       pub_etrans_calculate     = esdf_boolean('etrans_calculate',.false.)
       pub_etrans_same_leads    = esdf_boolean('etrans_same_leads',.true.)
       pub_etrans_bulk          = esdf_boolean('etrans_bulk',.false.)
       pub_etrans_write_setup   = esdf_boolean('etrans_write_setup',.true.)
       pub_etrans_ecmplx        = esdf_physical('etrans_ecmplx',0.001_DP, 'hartree')
       pub_etrans_emax          = esdf_physical('etrans_emax',0.2_DP, 'hartree')
       pub_etrans_emin          = esdf_physical('etrans_emin',-0.2_DP, 'hartree')
       pub_etrans_enum          = esdf_integer('etrans_enum',50)
       pub_etrans_source        = esdf_integer('etrans_source',1)
       pub_etrans_drain         = esdf_integer('etrans_drain',2)

       ! fc: parameters for phonon task
       pub_phonon_farming_task = esdf_integer('phonon_farming_task',0)
       pub_phonon_disp = esdf_physical('phonon_finite_disp',0.1_DP,'bohr')
       pub_phonon_fmax = esdf_physical('phonon_fmax',0.005_DP,'ha/bohr')
       pub_have_phonon_disp_list = .false.
       pub_num_disp = -1
       pub_have_phonon_disp_list = esdf_block('phonon_disp_list',pub_num_disp)
       pub_phonon_tmin = esdf_physical('phonon_tmin',0.0_DP,'hartree')
       pub_phonon_tmax = esdf_physical('phonon_tmax',0.002_DP,'hartree')
       pub_phonon_deltat = esdf_physical('phonon_deltat',1.5e-5_DP,'hartree')
       pub_phonon_min_freq = esdf_physical('phonon_min_freq',3.6e-06_DP,'hartree')

       pub_bs_kpoint_path_spacing = esdf_physical('bs_kpoint_path_spacing', &
            0.1889727_DP,'1/bohr')
       pub_bs_num_eigenvalues  = esdf_integer('bs_num_eigenvalues',-1)
       buf                     = esdf_string('bs_unfold','0 0 0')
       read(buf,*) pub_bs_unfold(:)
       pub_bs_method           = esdf_string('bs_method', 'TB')

       pub_do_tddft                   = esdf_boolean('do_tddft', .false.)
       pub_tddft_maximum_energy       = &
            esdf_physical('tddft_maximum_energy',1.0_DP,'hartree')
       pub_tddft_resolution           = &
            esdf_physical('tddft_resolution',0.001_DP,'hartree')
       buf = esdf_string('tddft_propagation_method', 'CRANKNICHOLSON')
       pub_tddft_propagation_method   = buf(1:20)
       pub_tddft_sparsity_level       = esdf_integer('tddft_sparsity_level', 0)
       pub_tddft_tammdancoff          = &
            esdf_boolean('tddft_tammdancoff', .false.)
       tdbuf                          = &
            esdf_string('tddft_dipole_kick_strength','0.0 0.0 0.0')
       read(tdbuf,*) pub_tddft_dipole_kick_strength(:)
       tddft_xc_functional            = esdf_string('tddft_xc_functional','LDA')
       pub_tddft_hamiltonian_mixing   = &
            esdf_integer('tddft_hamiltonian_mixing', 0)
       pub_tddft_damping              = &
            esdf_physical('tddft_damping',0.0_DP,'hartree')
       pub_tddft_enforced_idempotency = &
            esdf_boolean('tddft_enforced_idempotency', .false.)
       pub_tddft_maxit_hotelling      = &
            esdf_integer('tddft_maxit_hotelling', 50)
       pub_tddft_max_resid_hotelling  = &
            esdf_double('tddft_max_resid_hotelling', 1.0e-18_DP)
       pub_tddft_inv_overlap_exact    = &
            esdf_boolean('tddft_inv_overlap_exact', .true.)

       pub_is_implicit_solvent        = esdf_boolean('is_implicit_solvent',.false.)
       pub_is_smeared_ion_rep         = esdf_boolean('is_smeared_ion_rep', .false.)
       pub_is_include_cavitation      = esdf_boolean('is_include_cavitation', &
            .false.)
       pub_is_density_threshold       = esdf_double('is_density_threshold',0.00078_DP)
       pub_is_solvation_beta          = esdf_double('is_solvation_beta',1.3_DP)

       ! jd: Default bulk permittivity is 80.0 if using implicit solvent,
       !     1.0 otherwise
       if(pub_is_implicit_solvent) then
          is_default_bulk_permittivity = 80.0_DP
       else
          is_default_bulk_permittivity = 1.0_DP
       end if
       pub_is_bulk_permittivity = esdf_double('is_bulk_permittivity', &
            is_default_bulk_permittivity)
       pub_is_multigrid_defect_error_t = &
            esdf_double('is_multigrid_defect_error_tol', 1D-2)
       pub_is_multigrid_error_tol = esdf_double('is_multigrid_error_tol', 1D-5)
       pub_is_smeared_ion_width = esdf_double('is_smeared_ion_width', 0.8_DP)
       pub_is_surface_thickness = esdf_double('is_surface_thickness', 0.0002_DP)
       pub_is_solvent_surface_tension = &
            esdf_physical('is_solvent_surface_tension', 4.7624D-5,'ha/bohr**2')
       ! jd: Default is an estimate for H2O from J.Chem.Phys.124 074103 (2006).
       pub_is_density_threshold = esdf_double('is_density_threshold',0.00078_DP)

       pub_is_multigrid_max_iters = esdf_integer('is_multigrid_max_iters', 100)
       pub_is_multigrid_nlevels = esdf_integer('is_multigrid_nlevels', 4)
       pub_is_discretization_order = esdf_integer('is_discretization_order', 8)
       pub_is_bc_coarseness = esdf_integer('is_bc_coarseness', 5)
       pub_is_bc_surface_coarseness = esdf_integer('is_bc_surface_coarseness', &
            1)
       pub_is_bc_threshold = esdf_double('is_bc_threshold', 1D-9)
       pub_is_check_solv_energy_grad = &
            esdf_boolean('is_check_solv_energy_grad', .false.)
       pub_is_solvation_method = esdf_string('is_solvation_method', 'direct')
       pub_is_dielectric_model = esdf_string('is_dielectric_model', &
            'fix_initial')
       pub_is_dielectric_function = esdf_string('is_dielectric_function', 'fgf')
       pub_is_solvation_output_detail = &
            esdf_string('is_solvation_output_detail', 'none')

       ! jd: Turn on smeared ions automatically if using implicit solvent
       if (pub_is_implicit_solvent .and. .not. pub_is_smeared_ion_rep) then
          write(stdout,'(/a)') 'WARNING: is_smeared_ion_rep was automatically &
               &set to true for you because you''ve used is_implicit_solvent .true.,&
               & and an implicit solvent calculation cannot be performed without the&
               & smeared-ion representation. Add ''is_smeared_ion_rep T'' to get rid&
               & of this warning.'
          pub_is_smeared_ion_rep = .true.
       end if
       ! jd: --------- End of implicit solvent stuff ---------

       ! jd: --------- Open BC keywords ---------
       pub_openbc_ion_ion = esdf_boolean('openbc_ion_ion', .false.)
       pub_openbc_hartree = esdf_boolean('openbc_hartree', .false.)
       pub_openbc_pspot = esdf_boolean('openbc_pspot', .false.)
       pub_openbc_pspot_finetune_nptsx = &
            esdf_integer('openbc_pspot_finetune_nptsx', 100000)
       pub_openbc_pspot_finetune_alpha = &
            esdf_double('openbc_pspot_finetune_alpha', 0.3_DP)
       pub_openbc_pspot_finetune_f = &
            esdf_integer('openbc_pspot_finetune_f', -1)
       ! jd: --------- End of open BC stuff ---------


       ! ars: Select between tightbox and spherical waves representation when
       ! ars: doing restart
       if (pub_read_sw_ngwfs.and.read_tightbox_ngwfs) then
          ! ars : the program uses tightbox_ngwfs unless we specify in the
          ! ars : input read_tightbox_ngwfs = FALSE pub_read_sw_ngwfs = TRUE
          write(stdout,'(/a)') 'WARNING: pub_read_sw_ngwfs set to FALSE. &
               &Restart will be done using .tightbox_ngwfs file'
          pub_read_sw_ngwfs=.false.
       endif

       ! jd: Select integration routine for HF exhange
       pub_hfx_integration_variant = esdf_integer('hfx_integration', 2)
       ! qoh: HF exchange parameters
       pub_hfx_nlpp_for_exchange = esdf_boolean('hfx_nlpp_for_exchange',.false.)
       pub_hfx_radial_segments = esdf_integer('hfx_radial_segments',50)
       pub_hfx_angular_segments = esdf_integer('hfx_angular_segments',50)       
       pub_hfx_cheb_order = esdf_integer('hfx_cheb_order', 16)
       pub_hfx_cheb_intervals = esdf_integer('hfx_cheb_intervals', 14)
       pub_hfx_cheb_a_batchsize = esdf_integer('hfx_cheb_a_batchsize', 10)
       pub_hfx_cheb_b_batchsize = esdf_integer('hfx_cheb_b_batchsize', 10)
       pub_hfx_max_l = esdf_integer('hfx_max_l', 4)
       pub_hfx_max_zeros = esdf_integer('hfx_max_zeros', 10)
       pub_hfx_read_vmatrix = esdf_boolean('hfx_read_vmatrix', .false.)
       pub_hfx_write_vmatrix = esdf_boolean('hfx_write_vmatrix', .false.)
       pub_hfx_read_xmatrix = esdf_boolean('hfx_read_xmatrix', .false.)
       pub_hfx_write_xmatrix = esdf_boolean('hfx_write_xmatrix', .false.)
       pub_hfx_debug = esdf_boolean('hfx_debug', .false.)


       pub_do_bandstructure = esdf_block('bs_kpoint_path',nlines)
       if (pub_do_bandstructure) then
          bs_start = .true.
          pub_bs_kpoint_path_length = 0
          do iline=1,nlines
             if (index(block_data(iline),'BREAK') > 0) then
                bs_start = .true.
             else
                if (.not. bs_start) pub_bs_kpoint_path_length = &
                     pub_bs_kpoint_path_length + 1
                bs_start = .false.
             end if
          end do
          if (pub_bs_kpoint_path_length < 1) then
             write(stdout,'(/a/23x,a)') 'WARNING in get_rundat: &
                  &invalid path in bs_kpoint_path block', &
                  'bandstructure calculation will be skipped'
             pub_do_bandstructure = .false.
          end if
          if (pub_bs_method /= 'TB' .and. pub_bs_method /= 'KP') then
             write(stdout,'(/2a/23x,a)') 'WARNING in get_rundat: &
                  &unknown bs_method: ',pub_bs_method, &
                  'bandstructure calculation will be skipped'
             pub_do_bandstructure = .false.
          end if
       end if

       ! cks: Further logic for setting maxit_pen to avoid spoiling denskern
       ! cks: in single-point energy restarts
       if (task == 'SINGLEPOINT') then
          if (read_denskern .and. read_tightbox_ngwfs) then
             maxit_pen = 0
          elseif(read_denskern .and. pub_read_sw_ngwfs) then ! ars
             maxit_pen = 0                                   ! ars
             ! ddor: Bare response calculation for DFT+U parameters
          elseif (pub_hub_max_iter .eq. -1) then
             maxit_pen = 0
          endif
       endif

       ! lpl: NBO stuff
       ! lpl: 26092011 - nbo_init_lclowdin now defaults to TRUE
       pub_write_nbo = esdf_boolean('write_nbo',.false.)
       pub_nbo_init_lclowdin  = esdf_boolean('nbo_init_lclowdin',.true.)
       pub_nbo_write_lclowdin = esdf_boolean('nbo_write_lclowdin',.false.)
       pub_nbo_write_npacomp  = esdf_boolean('nbo_write_npacomp',.false.)
       pub_nbo_scale_dm       = esdf_boolean('nbo_scale_dm',.true.)

       pub_nbo_pnao_analysis  = esdf_boolean('nbo_pnao_analysis',.false.)

       pub_nbo_aopnao_scheme = esdf_string('nbo_aopnao_scheme','ORIGINAL')
       if(pub_nbo_aopnao_scheme /= 'ORIGINAL' .and. &
          pub_nbo_aopnao_scheme /= 'DIAGONALIZATION' .and. &
          pub_nbo_aopnao_scheme /= 'NONE') then
             write(stdout,'(a,a,a)') &
                  'ERROR: Invalid nbo_aopnao_scheme specified (', &
                       pub_nbo_aopnao_scheme,')'
             call comms_abort
       end if

       pub_plot_nbo   = esdf_boolean('plot_nbo',.false.)
       if(pub_write_nbo .and. pub_plot_nbo) write(stdout,'(a)') 'ERROR: &
            &pub_write_nbo and pub_plot_nbo cannot be simulataneously T'
       pub_nbo_plotorbtype = esdf_string('nbo_plot_orbtype','')
       ! lpl: NBO stuff

    end if

    !...and broadcasts to all nodes
    call comms_bcast(pub_root_node_id, task)
    call comms_bcast(pub_root_node_id, cutoff_energy)
    call comms_bcast(pub_root_node_id, kernel_cutoff)
    call comms_bcast(pub_root_node_id, xc_functional)
    call comms_bcast(pub_root_node_id, pub_libxc_x_func_id)
    call comms_bcast(pub_root_node_id, pub_libxc_c_func_id)
    call comms_bcast(pub_root_node_id, pub_charge)
    call comms_bcast(pub_root_node_id, pub_spin)
    call comms_bcast(pub_root_node_id, pub_spin_polarised)
    call comms_bcast(pub_root_node_id, pub_spin_polarised)
    call comms_bcast(pub_root_node_id, pub_constant_efield)
    call comms_bcast(pub_root_node_id, fftbox_pref)
    call comms_bcast(pub_root_node_id, ppd_npoints)
    call comms_bcast(pub_root_node_id, psinc_spacing)
    call comms_bcast(pub_root_node_id, pub_old_input)
    call comms_bcast(pub_root_node_id, pub_devel_code)
    call comms_bcast(pub_root_node_id, pub_nnho)
    call comms_bcast(pub_root_node_id, pub_ngwf_halo)
    call comms_bcast(pub_root_node_id, pub_nonsc_forces)
    call comms_bcast(pub_root_node_id, pub_external_pressure)
    call comms_bcast(pub_root_node_id, pub_isosurface_cutoff) 
    call comms_bcast(pub_root_node_id, pub_smoothing_factor) 
    call comms_bcast(pub_root_node_id, pub_maxit_palser_mano)
    call comms_bcast(pub_root_node_id, maxit_pen)
    call comms_bcast(pub_root_node_id, pen_param)
    call comms_bcast(pub_root_node_id, minit_lnv)
    call comms_bcast(pub_root_node_id, maxit_lnv)
    call comms_bcast(pub_root_node_id, maxit_lnv_coarse)
    call comms_bcast(pub_root_node_id, lnv_threshold_orig)
    call comms_bcast(pub_root_node_id, lnv_cg_type)
    call comms_bcast(pub_root_node_id, lnv_cg_max_step)
    call comms_bcast(pub_root_node_id, exact_lnv)
    call comms_bcast(pub_root_node_id, old_lnv)
    call comms_bcast(pub_root_node_id, pub_lnv_check_trial_steps)
    call comms_bcast(pub_root_node_id, pub_kerfix)
    call comms_bcast(pub_root_node_id, pub_maxit_kernel_fix)
    call comms_bcast(pub_root_node_id, pub_initial_dens_realspace)
    call comms_bcast(pub_root_node_id, pub_kernel_diis)
    call comms_bcast(pub_root_node_id, pub_kernel_diis_type)
    call comms_bcast(pub_root_node_id, pub_kernel_diis_max)
    call comms_bcast(pub_root_node_id, pub_kernel_diis_maxit)
    call comms_bcast(pub_root_node_id, pub_kernel_diis_thres)
    call comms_bcast(pub_root_node_id, pub_kernel_diis_liter)
    call comms_bcast(pub_root_node_id, pub_kernel_diis_c_in)
    call comms_bcast(pub_root_node_id, pub_kernel_diis_c_out)
    call comms_bcast(pub_root_node_id, pub_kernel_diis_criteria)
    call comms_bcast(pub_root_node_id, maxit_ngwf_cg)
    call comms_bcast(pub_root_node_id, ngwf_cg_type)
    call comms_bcast(pub_root_node_id, ngwf_cg_max_step)
    call comms_bcast(pub_root_node_id, pub_elec_cg_max)
    call comms_bcast(pub_root_node_id, precond_scheme)
    call comms_bcast(pub_root_node_id, k_zero)
    call comms_bcast(pub_root_node_id, r_precond)
    call comms_bcast(pub_root_node_id, precond_recip)
    call comms_bcast(pub_root_node_id, precond_real)
    call comms_bcast(pub_root_node_id, smooth_scheme)
    call comms_bcast(pub_root_node_id, r_smooth)
    call comms_bcast(pub_root_node_id, k_smooth)
    call comms_bcast(pub_root_node_id, occ_mix)
    call comms_bcast(pub_root_node_id, pub_kernel_update)
    call comms_bcast(pub_root_node_id, pub_kernel_christoffel_update)
    call comms_bcast(pub_root_node_id, maxit_hotelling)
    call comms_bcast(pub_root_node_id, max_resid_hotelling)
    call comms_bcast(pub_root_node_id, maxit_ngwf_diis)
    call comms_bcast(pub_root_node_id, ngwf_diis_size)
    call comms_bcast(pub_root_node_id, ngwf_diis_reset)
    call comms_bcast(pub_root_node_id, ngwf_diis_max_step)
    call comms_bcast(pub_root_node_id, ngwf_threshold_orig)
    call comms_bcast(pub_root_node_id, pub_elec_energy_tol)
    call comms_bcast(pub_root_node_id, pub_elec_force_tol)
    call comms_bcast(pub_root_node_id, pub_ngwf_max_grad)
    call comms_bcast(pub_root_node_id, delta_e_conv)
    call comms_bcast(pub_root_node_id, locpot_int_batch_size)
    call comms_bcast(pub_root_node_id, kinetic_int_batch_size)
    call comms_bcast(pub_root_node_id, density_batch_size)
    call comms_bcast(pub_root_node_id, ngwf_grad_batch_size)
    call comms_bcast(pub_root_node_id, pub_dense_threshold)
    call comms_bcast(pub_root_node_id, pub_comms_group_size)
    call comms_bcast(pub_root_node_id, pub_ovlp_for_nonlocal)
    call comms_bcast(pub_root_node_id, use_space_filling_curve)
    call comms_bcast(pub_root_node_id, coreham_denskern_guess)
    call comms_bcast(pub_root_node_id, pub_check_atoms)
    call comms_bcast(pub_root_node_id, pub_locpot_scheme)
    call comms_bcast(pub_root_node_id, pub_smooth_projectors)
    call comms_bcast(pub_root_node_id, pub_smooth_loc_pspot)
    call comms_bcast(pub_root_node_id, pub_odd_psinc_grid)
    call comms_bcast(pub_root_node_id, pub_even_psinc_grid)
    call comms_bcast(pub_root_node_id, pub_realspace_projectors)
    call comms_bcast(pub_root_node_id, output_detail)
    call comms_bcast(pub_root_node_id, paw_output_detail)
    call comms_bcast(pub_root_node_id, timings_level)
    call comms_bcast(pub_root_node_id, pub_write_params)
    call comms_bcast(pub_root_node_id, pub_write_forces)
    call comms_bcast(pub_root_node_id, pub_write_positions)
    call comms_bcast(pub_root_node_id, pub_write_xyz)
    call comms_bcast(pub_root_node_id, print_qc)
    call comms_bcast(pub_root_node_id, pub_cube_format)
    call comms_bcast(pub_root_node_id, pub_dx_format)
    call comms_bcast(pub_root_node_id, pub_dx_coarse)
    call comms_bcast(pub_root_node_id, pub_dx_sig_digits)
    call comms_bcast(pub_root_node_id, pub_grd_format)
    call comms_bcast(pub_root_node_id, pub_write_density_plot)
    call comms_bcast(pub_root_node_id, pub_write_ngwf_plot)
    call comms_bcast(pub_root_node_id, pub_write_ngwf_grad_plot)
    call comms_bcast(pub_root_node_id, pub_write_ngwf_radial)
    call comms_bcast(pub_root_node_id, pub_write_ngwf_grad_radial)
    call comms_bcast(pub_root_node_id, pub_write_radial_step)
    call comms_bcast(pub_root_node_id, pub_write_radial_smear)
    call comms_bcast(pub_root_node_id, read_denskern)
    call comms_bcast(pub_root_node_id, write_denskern)
    call comms_bcast(pub_root_node_id, read_tightbox_ngwfs)
    call comms_bcast(pub_root_node_id, write_tightbox_ngwfs)
    call comms_bcast(pub_root_node_id, pub_write_converged_dk_ngwfs)
    call comms_bcast(pub_root_node_id, pub_read_sw_ngwfs)
    call comms_bcast(pub_root_node_id, pub_write_sw_ngwfs)
    call comms_bcast(pub_root_node_id, pub_write_max_l)
    call comms_bcast(pub_root_node_id, pub_read_max_l)
    call comms_bcast(pub_root_node_id, pub_extra_n_sw)
    call comms_bcast(pub_root_node_id, pub_paw)
    call comms_bcast(pub_root_node_id, pub_fine_grid_scale)
    call comms_bcast(pub_root_node_id, pub_dbl_grid_scale)
    call comms_bcast(pub_root_node_id, pub_aug_funcs_recip)
    call comms_bcast(pub_root_node_id, pub_dispersion)
    call comms_bcast(pub_root_node_id, pub_vdw_dcoeff)
    call comms_bcast(pub_root_node_id, pub_hub_max_iter)
    call comms_bcast(pub_root_node_id, pub_hub_energy_tol)
    call comms_bcast(pub_root_node_id, pub_hub_conv_win)
    call comms_bcast(pub_root_node_id, pub_hub_proj_mixing)
    call comms_bcast(pub_root_node_id, pub_hub_functional)
    call comms_bcast(pub_root_node_id, pub_hub_tensor_corr)
    call comms_bcast(pub_root_node_id, pub_hub_ngwf_spin_thr)

    call comms_bcast(pub_root_node_id, pub_cdft_atom_charge)
    call comms_bcast(pub_root_node_id, pub_cdft_atom_spin)
    call comms_bcast(pub_root_node_id, pub_cdft_group_charge)
    call comms_bcast(pub_root_node_id, pub_cdft_group_spin)
    call comms_bcast(pub_root_node_id, pub_cdft_group_charge_diff)
    call comms_bcast(pub_root_node_id, pub_cdft_group_spin_diff)
    call comms_bcast(pub_root_node_id, pub_cdft_hubbard)
    call comms_bcast(pub_root_node_id, pub_ci_cdft)
    call comms_bcast(pub_root_node_id, pub_cdft_type)
    call comms_bcast(pub_root_node_id, pub_cdft_group_charge_target)
    call comms_bcast(pub_root_node_id, pub_cdft_group_spin_target)
    call comms_bcast(pub_root_node_id, pub_cdft_group_charge_diff_target)
    call comms_bcast(pub_root_node_id, pub_cdft_group_spin_diff_target)
    call comms_bcast(pub_root_node_id, pub_cdft_group_u)
    call comms_bcast(pub_root_node_id, pub_cdft_group_diff_u)
    call comms_bcast(pub_root_node_id, maxit_cdft_u_cg)

!CW
    call comms_bcast(pub_root_node_id, pub_write_polarisation_plot)
    call comms_bcast(pub_root_node_id, pub_dmft_ignore_type)
    call comms_bcast(pub_root_node_id, pub_dmft_KS_shift)
    call comms_bcast(pub_root_node_id, pub_dmft_purify_sc)
    call comms_bcast(pub_root_node_id, pub_dmft_kernel_mix)
    call comms_bcast(pub_root_node_id, pub_dmft_mu_diff_max)
    call comms_bcast(pub_root_node_id, pub_dmft_switch_off_proj_order)
    call comms_bcast(pub_root_node_id, pub_dmft_fully_sc)
    call comms_bcast(pub_root_node_id, pub_dmft_fully_sc_h)
    call comms_bcast(pub_root_node_id, pub_dmft_spoil_kernel)
    call comms_bcast(pub_root_node_id, pub_dmft_kernel)
    call comms_bcast(pub_root_node_id, pub_dmft_sc)
    call comms_bcast(pub_root_node_id, pub_dmft_order_proj)
    call comms_bcast(pub_root_node_id, pub_dmft_in_bohr)
    call comms_bcast(pub_root_node_id, pub_dmft_temp)
    call comms_bcast(pub_root_node_id, pub_dmft_smear)
    call comms_bcast(pub_root_node_id, pub_dmft_smear_T)
    call comms_bcast(pub_root_node_id, pub_dmft_smear_shift)
    call comms_bcast(pub_root_node_id, pub_dmft_smear_eta)
    call comms_bcast(pub_root_node_id, pub_dmft_smear_w)
    call comms_bcast(pub_root_node_id, pub_dmft_smear_shift)
    call comms_bcast(pub_root_node_id, pub_dmft_points)
    call comms_bcast(pub_root_node_id, pub_dmft_optics)
    call comms_bcast(pub_root_node_id, pub_dmft_skip_energy)
    call comms_bcast(pub_root_node_id, pub_dmft_optics_window)
    call comms_bcast(pub_root_node_id, pub_dmft_optics_x1)
    call comms_bcast(pub_root_node_id, pub_dmft_optics_y1)
    call comms_bcast(pub_root_node_id, pub_dmft_optics_z1)
    call comms_bcast(pub_root_node_id, pub_dmft_embed_iter)
    call comms_bcast(pub_root_node_id, pub_dmft_nkpoints)
    call comms_bcast(pub_root_node_id, pub_dmft_kpoints_sym)
    call comms_bcast(pub_root_node_id, pub_dmft_kpoints_kernel_gamma)
    call comms_bcast(pub_root_node_id, pub_dmft_lin_scaling )
    call comms_bcast(pub_root_node_id, pub_dmft_nval)
    call comms_bcast(pub_root_node_id, pub_dmft_win)
    call comms_bcast(pub_root_node_id, pub_dmft_scaling_tol)
    call comms_bcast(pub_root_node_id, pub_dmft_split)
    call comms_bcast(pub_root_node_id, pub_dmft_splitk)
    call comms_bcast(pub_root_node_id, pub_dmft_nbo)
    call comms_bcast(pub_root_node_id, pub_dmft_impose_same_coeffs)
    call comms_bcast(pub_root_node_id, pub_dmft_impose_chem_spin)
    call comms_bcast(pub_root_node_id, pub_dmft_scaling_iter)
    call comms_bcast(pub_root_node_id, pub_dmft_scaling_nmpi)
    call comms_bcast(pub_root_node_id, pub_dmft_scaling_meth)
    call comms_bcast(pub_root_node_id, pub_dmft_invert_overlap)
    call comms_bcast(pub_root_node_id, pub_dmft_local_scratch)
    call comms_bcast(pub_root_node_id, pub_dmft_scaling_maxspace)
    call comms_bcast(pub_root_node_id, pub_dmft_scaling_cutoff)
    call comms_bcast(pub_root_node_id, pub_dmft_scaling_cutoff_h)
    call comms_bcast(pub_root_node_id, pub_dmft_scaling_tail)
    call comms_bcast(pub_root_node_id, pub_dmft_mu_order)
    call comms_bcast(pub_root_node_id, pub_dmft_embed_mix)
    call comms_bcast(pub_root_node_id, pub_dmft_optics_i1)
    call comms_bcast(pub_root_node_id, pub_dmft_optics_i2)
    call comms_bcast(pub_root_node_id, pub_dmft_nmu_loop)
    call comms_bcast(pub_root_node_id, pub_dmft_plot_real_space)
    call comms_bcast(pub_root_node_id, pub_dmft_plot_real_space_sigma)
    call comms_bcast(pub_root_node_id, pub_dmft_paramagnetic)
    call comms_bcast(pub_root_node_id, pub_dmft_integrate_green)
    call comms_bcast(pub_root_node_id, pub_dmft_plot_all_proj)
    call comms_bcast(pub_root_node_id, pub_dmft_rotate_green)
    call comms_bcast(pub_root_node_id, pub_dmft_chem_shift)
    call comms_bcast(pub_root_node_id, pub_dmft_e_plot)
    call comms_bcast(pub_root_node_id, pub_dmft_emin)
    call comms_bcast(pub_root_node_id, pub_dmft_emax)
    call comms_bcast(pub_root_node_id, pub_dmft_dos_min)
    call comms_bcast(pub_root_node_id, pub_dmft_dos_max)
    call comms_bcast(pub_root_node_id, pub_dmft_cutoff_tail)
    call comms_bcast(pub_root_node_id, pub_dmft_cutoff_small)
    call comms_bcast(pub_root_node_id, pub_dmft_free_green_frequ)
    call comms_bcast(pub_root_node_id, pub_dmft_norm_proj)
!END CW

    call comms_bcast(pub_root_node_id, cdft_cg_type)
    call comms_bcast(pub_root_node_id, cdft_cg_threshold)
    call comms_bcast(pub_root_node_id, pub_cdft_cg_max)
    call comms_bcast(pub_root_node_id, pub_cdft_max_grad)
    call comms_bcast(pub_root_node_id, pub_cdft_elec_energy_tol)
    call comms_bcast(pub_root_node_id, cdft_cg_max_step)
    call comms_bcast(pub_root_node_id, pub_cdft_guru)
    call comms_bcast(pub_root_node_id, pub_cdft_continuation)

    call comms_bcast(pub_root_node_id, pub_coulomb_radius)
    call comms_bcast(pub_root_node_id, pub_coulomb_length)
    call comms_bcast(pub_root_node_id, pub_coulomb_cutoff_type)
    call comms_bcast(pub_root_node_id, pub_coulomb_cutoff_write_int)
    call comms_bcast(pub_root_node_id, pub_mt_cutoff)
    call comms_bcast(pub_root_node_id, pub_do_properties)
    call comms_bcast(pub_root_node_id, pub_num_eigenvalues)
    call comms_bcast(pub_root_node_id, pub_homo_dens_plot)
    call comms_bcast(pub_root_node_id, pub_lumo_dens_plot)
    call comms_bcast(pub_root_node_id, pub_homo_plot)
    call comms_bcast(pub_root_node_id, pub_lumo_plot)
    call comms_bcast(pub_root_node_id, pub_dos_smear)
    call comms_bcast(pub_root_node_id, pub_ldos_smear)
    call comms_bcast(pub_root_node_id, pub_popn_calculate)
    call comms_bcast(pub_root_node_id, pub_popn_bond_cutoff)
    call comms_bcast(pub_root_node_id, pub_ngwf_analysis)
    call comms_bcast(pub_root_node_id, pub_polarisation_calculate)
    call comms_bcast(pub_root_node_id, pub_spread_calculate)
    call comms_bcast(pub_root_node_id, geom_modulus_est)
    call comms_bcast(pub_root_node_id, geom_frequency_est)
    call comms_bcast(pub_root_node_id, geom_energy_tol)
    call comms_bcast(pub_root_node_id, geom_force_tol)
    call comms_bcast(pub_root_node_id, geom_disp_tol)
    call comms_bcast(pub_root_node_id, geom_max_iter)
    call comms_bcast(pub_root_node_id, geom_convergence_win)
    call comms_bcast(pub_root_node_id, geom_continuation)
    call comms_bcast(pub_root_node_id, geom_backup_iter)
    call comms_bcast(pub_root_node_id, geom_reset_dk_ngwfs_iter)
    call comms_bcast(pub_root_node_id, geom_lbfgs_max_updates)
    call comms_bcast(pub_root_node_id, geom_lbfgs_block_length)
    call comms_bcast(pub_root_node_id, pub_geom_reuse_dk_ngwfs)
    call comms_bcast(pub_root_node_id, geom_print_inv_hessian)
    call comms_bcast(pub_root_node_id, geom_method)
    call comms_bcast(pub_root_node_id, geom_lbfgs)
    call comms_bcast(pub_root_node_id, tssearch_method)
    call comms_bcast(pub_root_node_id, tssearch_lstqst_protocol)
    call comms_bcast(pub_root_node_id, tssearch_qst_max_iter)
    call comms_bcast(pub_root_node_id, tssearch_cg_max_iter)
    call comms_bcast(pub_root_node_id, tssearch_force_tol)
    call comms_bcast(pub_root_node_id, tssearch_disp_tol)
    call comms_bcast(pub_root_node_id, md_restart)
    call comms_bcast(pub_root_node_id, pub_md_properties)
    call comms_bcast(pub_root_node_id, md_delta_t)
    call comms_bcast(pub_root_node_id, md_num_iter)
    call comms_bcast(pub_root_node_id, md_reset_dkn_ngwfs)
    call comms_bcast(pub_root_node_id, mts_xi)
    call comms_bcast(pub_root_node_id, mts_nstep)
    call comms_bcast(pub_root_node_id, mts_ngwf_threshold)
    call comms_bcast(pub_root_node_id, mts_lnv_threshold)
    call comms_bcast(pub_root_node_id, mix_dkn_num)
    call comms_bcast(pub_root_node_id, mix_ngwfs_num)
    call comms_bcast(pub_root_node_id, mix_ngwfs_coeff)
    call comms_bcast(pub_root_node_id, mix_dkn_type)
    call comms_bcast(pub_root_node_id, mix_ngwfs_type)
    call comms_bcast(pub_root_node_id, mix_local_length)
    call comms_bcast(pub_root_node_id, mix_local_smear)
    call comms_bcast(pub_root_node_id, cond_read_denskern)
    call comms_bcast(pub_root_node_id, cond_read_tightbox_ngwfs)
    call comms_bcast(pub_root_node_id, cond_fixed_shift)
    call comms_bcast(pub_root_node_id, cond_calc_max_eigen)
    call comms_bcast(pub_root_node_id, cond_num_states)
    call comms_bcast(pub_root_node_id, cond_num_extra_states)
    call comms_bcast(pub_root_node_id, cond_num_extra_its)
    call comms_bcast(pub_root_node_id, cond_kernel_cutoff)
    call comms_bcast(pub_root_node_id, cond_init_shift)
    call comms_bcast(pub_root_node_id, cond_shift_buffer)
    call comms_bcast(pub_root_node_id, cond_plot_joint_orbitals)
    call comms_bcast(pub_root_node_id, cond_plot_vc_orbitals)
    call comms_bcast(pub_root_node_id, pub_spectra_calculate)
    call comms_bcast(pub_root_node_id, pub_calc_mom_mat_els)
    call comms_bcast(pub_root_node_id, pub_calc_nonloc_comm)
    call comms_bcast(pub_root_node_id, pub_spec_cont_deriv)
    call comms_bcast(pub_root_node_id, pub_spec_nonloc_fin_diff)
    call comms_bcast(pub_root_node_id, pub_spectra_print_mat_els)
    call comms_bcast(pub_root_node_id, pub_scissor)
    call comms_bcast(pub_root_node_id, pub_opt_smear)
    call comms_bcast(pub_root_node_id, pub_bs_kpoint_path_spacing)
    call comms_bcast(pub_root_node_id, pub_do_bandstructure)
    call comms_bcast(pub_root_node_id, pub_bs_kpoint_path_length)
    call comms_bcast(pub_root_node_id, pub_bs_num_eigenvalues)
    call comms_bcast(pub_root_node_id, pub_bs_unfold)
    call comms_bcast(pub_root_node_id, pub_bs_method)
    call comms_bcast(pub_root_node_id, pub_do_tddft)
    call comms_bcast(pub_root_node_id, pub_tddft_maximum_energy)
    call comms_bcast(pub_root_node_id, pub_tddft_resolution)
    call comms_bcast(pub_root_node_id, pub_tddft_propagation_method)
    call comms_bcast(pub_root_node_id, pub_tddft_sparsity_level)
    call comms_bcast(pub_root_node_id, pub_tddft_tammdancoff)
    call comms_bcast(pub_root_node_id, pub_tddft_dipole_kick_strength)
    call comms_bcast(pub_root_node_id, tddft_xc_functional)
    call comms_bcast(pub_root_node_id, pub_tddft_hamiltonian_mixing)
    call comms_bcast(pub_root_node_id, pub_tddft_damping)
    call comms_bcast(pub_root_node_id, pub_tddft_enforced_idempotency)
    call comms_bcast(pub_root_node_id, pub_tddft_maxit_hotelling)
    call comms_bcast(pub_root_node_id, pub_tddft_max_resid_hotelling)
    call comms_bcast(pub_root_node_id, pub_tddft_inv_overlap_exact)
    call comms_bcast(pub_root_node_id, pub_hfx_integration_variant)
    call comms_bcast(pub_root_node_id, pub_hfx_nlpp_for_exchange)
    call comms_bcast(pub_root_node_id, pub_hfx_radial_segments)
    call comms_bcast(pub_root_node_id, pub_hfx_angular_segments)
    call comms_bcast(pub_root_node_id, pub_hfx_cheb_order)
    call comms_bcast(pub_root_node_id, pub_hfx_cheb_intervals)
    call comms_bcast(pub_root_node_id, pub_hfx_cheb_a_batchsize)
    call comms_bcast(pub_root_node_id, pub_hfx_cheb_b_batchsize)
    call comms_bcast(pub_root_node_id, pub_hfx_max_l)
    call comms_bcast(pub_root_node_id, pub_hfx_max_zeros)
    call comms_bcast(pub_root_node_id, pub_hfx_read_vmatrix)
    call comms_bcast(pub_root_node_id, pub_hfx_write_vmatrix)
    call comms_bcast(pub_root_node_id, pub_hfx_read_xmatrix)
    call comms_bcast(pub_root_node_id, pub_hfx_write_xmatrix)
    call comms_bcast(pub_root_node_id, pub_hfx_debug)
    call comms_bcast(pub_root_node_id, pub_is_implicit_solvent)
    call comms_bcast(pub_root_node_id, pub_is_smeared_ion_rep)
    call comms_bcast(pub_root_node_id, pub_is_include_cavitation)
    call comms_bcast(pub_root_node_id, pub_is_density_threshold)
    call comms_bcast(pub_root_node_id, pub_is_solvation_beta)
    call comms_bcast(pub_root_node_id, pub_is_bulk_permittivity)
    call comms_bcast(pub_root_node_id, pub_is_multigrid_defect_error_t)
    call comms_bcast(pub_root_node_id, pub_is_multigrid_error_tol)
    call comms_bcast(pub_root_node_id, pub_is_smeared_ion_width)
    call comms_bcast(pub_root_node_id, pub_is_surface_thickness)
    call comms_bcast(pub_root_node_id, pub_is_solvent_surface_tension)
    call comms_bcast(pub_root_node_id, pub_is_multigrid_max_iters)
    call comms_bcast(pub_root_node_id, pub_is_multigrid_nlevels)
    call comms_bcast(pub_root_node_id, pub_is_discretization_order)
    call comms_bcast(pub_root_node_id, pub_is_bc_coarseness)
    call comms_bcast(pub_root_node_id, pub_is_bc_surface_coarseness)
    call comms_bcast(pub_root_node_id, pub_is_bc_threshold)
    call comms_bcast(pub_root_node_id, pub_is_check_solv_energy_grad)
    call comms_bcast(pub_root_node_id, pub_is_solvation_method)
    call comms_bcast(pub_root_node_id, pub_is_dielectric_model)
    call comms_bcast(pub_root_node_id, pub_is_dielectric_function)
    call comms_bcast(pub_root_node_id, pub_is_solvation_output_detail)
    call comms_bcast(pub_root_node_id, pub_openbc_ion_ion)
    call comms_bcast(pub_root_node_id, pub_openbc_hartree)
    call comms_bcast(pub_root_node_id, pub_openbc_pspot)
    call comms_bcast(pub_root_node_id, pub_openbc_pspot_finetune_nptsx)
    call comms_bcast(pub_root_node_id, pub_openbc_pspot_finetune_alpha)
    call comms_bcast(pub_root_node_id, pub_openbc_pspot_finetune_f)

    ! lpl: NBO stuff
    call comms_bcast(pub_root_node_id, pub_write_nbo)
    call comms_bcast(pub_root_node_id, pub_nbo_init_lclowdin)
    call comms_bcast(pub_root_node_id, pub_nbo_write_lclowdin)
    call comms_bcast(pub_root_node_id, pub_nbo_write_npacomp)
    call comms_bcast(pub_root_node_id, pub_nbo_scale_dm)

    call comms_bcast(pub_root_node_id, pub_nbo_pnao_analysis)

    call comms_bcast(pub_root_node_id, pub_nbo_aopnao_scheme)

    call comms_bcast(pub_root_node_id, pub_plot_nbo)
    call comms_bcast(pub_root_node_id, pub_nbo_plotorbtype)
    ! lpl: NBO stuff



    ! smmd: Transmission coefficients
    call comms_bcast(pub_root_node_id, pub_etrans_calculate) 
    call comms_bcast(pub_root_node_id, pub_etrans_same_leads)
    call comms_bcast(pub_root_node_id, pub_etrans_bulk)
    call comms_bcast(pub_root_node_id, pub_etrans_write_setup)
    call comms_bcast(pub_root_node_id, pub_etrans_ecmplx)
    call comms_bcast(pub_root_node_id, pub_etrans_emax)
    call comms_bcast(pub_root_node_id, pub_etrans_emin)
    call comms_bcast(pub_root_node_id, pub_etrans_enum)
    call comms_bcast(pub_root_node_id, pub_etrans_source)
    call comms_bcast(pub_root_node_id, pub_etrans_drain)


    ! fc: deal with phonon parameters
    call comms_bcast(pub_root_node_id, pub_phonon_farming_task)
    call comms_bcast(pub_root_node_id, pub_phonon_disp)
    call comms_bcast(pub_root_node_id, pub_phonon_fmax)
    call comms_bcast(pub_root_node_id, pub_have_phonon_disp_list)
    if (pub_have_phonon_disp_list) then
       call comms_bcast(pub_root_node_id, pub_num_disp)
       allocate(pub_phonon_disp_list(1:pub_num_disp))
       if (pub_on_root) then
          do iline=1,pub_num_disp
             read(block_data(iline),*) pub_phonon_disp_list(iline)
          end do
       end if
       call comms_bcast(pub_root_node_id, pub_phonon_disp_list)
    end if
    call comms_bcast(pub_root_node_id, pub_phonon_tmin)
    call comms_bcast(pub_root_node_id, pub_phonon_tmax)
    call comms_bcast(pub_root_node_id, pub_phonon_deltat)
    call comms_bcast(pub_root_node_id, pub_phonon_min_freq)

    ! lr408: Set value of pub_cond_calculate according to task
    ! lr408: and check parameters are sensible
    if (task == 'COND' .or. task == 'PROPERTIES_COND') then
       pub_cond_calculate = .true.
       if (.not. cond_fixed_shift) then
          if (.not. cond_calc_max_eigen) then
             if (pub_on_root) write(stdout,*) 'Error in get_rundat: &
                  &for a conduction calculation with updating shift'
             if (pub_on_root) write(stdout,*) 'cond_calc_max_eigen must be set to true'
             cond_calc_max_eigen = .true.
          end if
       end if
    else
       pub_cond_calculate = .false.
    end if

    ! ndmh: check value of comms_group_size
    if (modulo(pub_total_num_nodes,pub_comms_group_size)/=0) then
       if (pub_on_root) write(stdout,*) 'WARNING in get_rundat: &
            &Comms group size is incompatible with number of nodes'
       if (pub_on_root) write(stdout,*) 'WARNING in get_rundat: &
            &Overriding to comms_group_size : 1'
       pub_comms_group_size = 1
    end if

    ! aam: checks
    if ( task=='MOLECULARDYNAMICS' ) then
       if (md_delta_t.lt.0.0_dp) then
          if (pub_on_root) write(stdout,*) 'Error in get_rundat: &
               &md_delta_t must be positive'
       endif
       if (md_num_iter.lt.0) then
          if (pub_on_root) write(stdout,*) 'Error in get_rundat: &
               &md_num_iter must be positive'
       endif
    endif

    ! pdh: sanity checks
    if (old_lnv .and. exact_lnv) then
       if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
            &setting both old_lnv and exact_lnv is forbidden'
       call comms_abort
    end if

    ! aam: more sanity checking
    if (geom_convergence_win < 2) then
       if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
            &geom_convergence_win < 2'
       call comms_abort
    end if

    if (task=='GEOMETRYOPTIMIZATION' .and. &
         (geom_max_iter < geom_convergence_win)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
            &geom_max_iter < geom_convergence_win'
       call comms_abort
    end if

    ! cks: set output detail integer parameter
    if (output_detail == 'BRIEF') then
       pub_output_detail = BRIEF
    else if (output_detail == 'NORMAL') then
       pub_output_detail = NORMAL
    else if (output_detail == 'VERBOSE') then
       pub_output_detail = VERBOSE
    end if

    ! ndmh: set PAW output detail integer parameter
    if (paw_output_detail == 'DEFAULT') then
       if (pub_output_detail==VERBOSE) then
          pub_paw_output_detail = NORMAL
       else if (pub_output_detail==NORMAL) then
          pub_paw_output_detail = BRIEF
       else if (pub_output_detail==BRIEF) then
          pub_paw_output_detail = BRIEF
       end if
    else if (paw_output_detail == 'BRIEF') then
       pub_paw_output_detail = BRIEF
    else if (paw_output_detail == 'NORMAL') then
       pub_paw_output_detail = NORMAL
    else if (paw_output_detail == 'VERBOSE') then
       pub_paw_output_detail = VERBOSE
    end if


    ! cks: automatic setup for properties task
    ! ddor: included TDDFT task
    ! lr408: Included conduction task
    ! ndmh: included properties_cond task
    if (task == 'PROPERTIES' .or. task == 'TDDFT' .or. &
         task == 'COND' .or. task == 'PROPERTIES_COND') then
       ! cks: force ngwf and denskern read from file
       ! jd: Warn the user that this is happening
       if(.not. read_denskern .and. pub_on_root) then
          write(stdout,'(/a,a,a/)') 'WARNING: read_denskern overridden to T&
               & because of task ',trim(task),'.'
       end if
       read_denskern =.true.

       ! ars: ONETEP uses tightbox_ngwfs unless we specify in the input
       ! ars: read_tightbox_ngwfs = FALSE;  pub_read_sw_ngwfs = TRUE
       if (read_tightbox_ngwfs) then
          pub_read_sw_ngwfs =.false.
       elseif (.not. pub_read_sw_ngwfs) then
          ! jd: Warn the user that this is happening
          if(pub_on_root) then
             write(stdout,'(/a,a,a/)') 'WARNING: read_tightbox_ngwfs overridden&
                  & to T because of task ',trim(task),'.'
          end if
          read_tightbox_ngwfs =.true.       ! ars
       endif                                ! ars
       maxit_pen =0
    endif

    if (task == 'PROPERTIES_COND') then
       ! ndmh: force cond ngwf and denskern read from file
       if(.not. read_denskern .and. pub_on_root) then
          write(stdout,'(/a,a,a/)') 'WARNING: cond_read_denskern overridden &
               &to T because of task ',trim(task),'.'
       end if
       cond_read_denskern =.true.

       if(pub_on_root) then
          write(stdout,'(/a,a,a/)') 'WARNING: cond_read_tightbox_ngwfs &
               &overridden to T because of task ',trim(task),'.'
       end if
       cond_read_tightbox_ngwfs = .true.
    endif

    ! pdh: spin polarisation checks
    if (pub_spin > 0) pub_spin_polarised = .true.

    ! ndmh: set boolean flag for cutoff coulomb
    if ((index(pub_coulomb_cutoff_type,'NONE')>0).or. &
         (index(pub_coulomb_cutoff_type,'none')>0))then
       pub_coulomb_cutoff=.false.
    else
       pub_coulomb_cutoff=.true.
    end if

    ! ndmh: sanity check for Cutoff Coulomb
    if (pub_coulomb_cutoff) then
      call utils_assert(pub_coulomb_radius >= 0.0_DP, &
           'coulomb_cutoff_radius must be positive')
      call utils_assert(pub_coulomb_length >= 0.0_DP, &
           'coulomb_cutoff_length must be positive')
    end if

    ! pdh: bandstructure path
    if (pub_do_bandstructure) then
       allocate(pub_bs_kpoint_path_start(3,pub_bs_kpoint_path_length), &
            stat=ierr)
       call utils_alloc_check('get_rundat','pub_bs_kpoint_path_start',ierr)
       allocate(pub_bs_kpoint_path_end(3,pub_bs_kpoint_path_length), &
            stat=ierr)
       call utils_alloc_check('get_rundat','pub_bs_kpoint_path_end',ierr)

       if (pub_on_root) then
          pub_do_bandstructure = esdf_block('bs_kpoint_path',nlines)
          bs_start = .true.
          ipath = 0
          do iline=1,nlines
             if (index(block_data(iline),'BREAK') > 0) then
                bs_start = .true.
             else
                read(block_data(iline),*) dummy_real
                if (.not. bs_start) then
                   ipath = ipath + 1
                   pub_bs_kpoint_path_end(:,ipath) = dummy_real
                else
                   bs_start = .false.
                end if
                if (ipath < pub_bs_kpoint_path_length) &
                     pub_bs_kpoint_path_start(:,ipath+1) = dummy_real
             end if
          end do
       end if

       call comms_bcast(pub_root_node_id, pub_bs_kpoint_path_start)
       call comms_bcast(pub_root_node_id, pub_bs_kpoint_path_end)

    end if

    ! qoh: Temporary assignment of tightbox FFT variables
    pub_tightbox_fft_coarse = ( pub_write_sw_ngwfs .or. pub_read_sw_ngwfs )
    pub_tightbox_fft_fine   = .false.

    ! ddor: Hubbard DFT+U input check
    pub_hubbard = .false. ! ddor: by default we do not carry out DFT+U
    pub_hub_calculating_u = .false. ! ddor: by default we renew the Hamiltonian
    if ( pub_hub_max_iter .gt. 1 ) then
       task = 'HUBBARDSCF'
    elseif (pub_hub_max_iter .eq. -1) then
       pub_hub_calculating_u = .true.
    endif
    if ( task == 'HUBBARDSCF' ) then
       pub_hubbard = .true.
       delta_e_conv = .false.
       if (pub_nnho) then
          if (pub_on_root) then
             write(stdout,'(a)') 'Error in get_rundat: pub_nnho==.true. is &
                  &incompatible with projector-self-consistent DFT+U.'
          end if
          call comms_abort
       endif
       if (ABS(pub_hub_proj_mixing) .gt. 1.0_DP) then
          if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
               &pub_hub_proj_mixing out of bounds, aborting'
          call comms_abort
       end if
    endif
    pub_hubbard_restart = .false.
    if (pub_hub_proj_mixing .lt. 0.0_DP) pub_hubbard_restart = .true.

    ! qoh: pub_hfx_nlpp_for_exchange isn't compatible with Hubbard
    pub_hfx_nlpp_for_exchange = ( pub_hfx_nlpp_for_exchange .and. &
         (.not. pub_hubbard) )
    if (pub_hfx_nlpp_for_exchange .and. pub_hubbard .and. pub_on_root) then
       write(stdout,'(/a)') 'WARNING: pub_hfx_nlpp_for_exchange set to FALSE, &
            &because it''s not compatible with HUBBARDSCF.'
    end if

    ! jd: Sanity check for HFX
    call utils_assert(pub_hfx_cheb_intervals > 0, &
         'hfx_cheb_intervals must be positive.')
    call utils_assert(pub_hfx_cheb_order > 1, &
         'hfx_cheb_order must be at least 2.')
    call utils_assert(pub_hfx_cheb_a_batchsize > 0, &
         'hfx_cheb_a_batchsize must be positive.')
    call utils_assert(pub_hfx_cheb_b_batchsize > 0, &
         'hfx_cheb_b_batchsize must be positive.')
    call utils_assert(pub_hfx_max_l >= 0, &
         'hfx_cheb_max_l must be non-negative.')
    call utils_assert(pub_hfx_max_zeros > 0, &
         'hfx_cheb_max_zeros must be positive.')

    ! gibo: CONSTRAINED_DFT checks
    pub_cdft = .false.              ! do not carry out cDFT by default

    ! gibo: avoid simultanous activation of conflicting cDFT-modes ==== START
    if (pub_cdft_atom_charge .AND. &
         (pub_cdft_atom_spin .OR. pub_cdft_group_charge .OR. &
         pub_cdft_group_spin .OR. pub_cdft_group_charge_diff .OR. &
         pub_cdft_group_spin_diff)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
            &pub_cdft_atom_charge is not compatible with other &
            &constrained-DFT modes, aborting'
       call comms_abort     
    end if

    if (pub_cdft_atom_spin .AND. &
         (pub_cdft_atom_charge .OR. pub_cdft_group_charge .OR. &
         pub_cdft_group_spin .OR. pub_cdft_group_charge_diff .OR. &
         pub_cdft_group_spin_diff)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
            &pub_cdft_atom_spin is not compatible with other &
            &constrained-DFT modes, aborting'
       call comms_abort
    end if

    if (pub_cdft_group_charge .AND. &
         (pub_cdft_atom_charge .OR. pub_cdft_atom_spin .OR. &
         pub_cdft_group_spin .OR. pub_cdft_group_charge_diff .OR. &
         pub_cdft_group_spin_diff ) ) then
        if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
             &pub_cdft_group_charge is not compatible with other &
             &constrained-DFT modes, aborting'
        call comms_abort
    end if

    if (pub_cdft_group_spin .AND. &
         (pub_cdft_atom_charge .OR. pub_cdft_atom_spin .OR. &
         pub_cdft_group_charge .OR. pub_cdft_group_charge_diff .OR. &
         pub_cdft_group_spin_diff)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
            &pub_cdft_group_spin is not compatible with other &
            &constrained-DFT modes, aborting'
       call comms_abort
    end if

    if (pub_cdft_group_charge_diff .AND. &
          (pub_cdft_atom_charge .OR. pub_cdft_atom_spin .OR. &
          pub_cdft_group_charge .OR. pub_cdft_group_spin .OR. &
          pub_cdft_group_spin_diff)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
            &pub_cdft_group_charge_diff is not compatible with &
            &other constrained-DFT modes, aborting'
       call comms_abort
    end if

    if (pub_cdft_group_spin_diff .AND. &
         (pub_cdft_atom_charge .OR. pub_cdft_atom_spin .OR. &
         pub_cdft_group_charge .OR. pub_cdft_group_spin .OR. &
         pub_cdft_group_charge_diff)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
            &pub_cdft_group_spin_diff is not compatible with other &
            &constrained-DFT modes , aborting'
       call comms_abort
    end if
    ! gibo: avoid simultanous activation of conflicting cDFT-modes ==== END


    ! gibo: force the user to set cdft_group_charge_diff=.T. if 
    ! gibo: cdft_group_charge_diff_target is found
    if ((ABS(pub_cdft_group_charge_diff_target) > 0._DP) .AND. &
         (.not.pub_cdft_group_charge_diff))  then
       if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
            &set pub_cdft_group_charge_diff=.T. for a group-charge-&
            &difference-constrained cDFT simulation, aborting'
       call comms_abort     
    end if

    ! gibo: force the user to set cdft_group_charge_diff_target > 0. if 
    ! gibo: cdft_group_charge_diff=.T. 
    if ((pub_cdft_group_charge_diff_target <= 0._DP) .AND. &
         (pub_cdft_group_charge_diff))  then
       if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
            &set cdft_group_charge_diff_target > 0. for a group-charge-&
            &difference-constrained cDFT simulation, aborting'
       call comms_abort     
    end if

    ! gibo: force the user to set cdft_group_spin_diff=.T. if
    ! gibo: cdft_group_spin_diff_target is found
    if ((ABS(pub_cdft_group_spin_diff_target) > 0._DP) .AND. &
         (.not.pub_cdft_group_spin_diff))  then
       if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
            & set pub_cdft_group_spin_diff=.T. for a group-spin-&
            &difference-constrained cDFT simulation, aborting'
       call comms_abort     
    end if

    ! gibo: force the user to set cdft_group_spin_diff_target > 0. if 
    ! gibo: cdft_group_spin_diff=.T. 
    if ((pub_cdft_group_spin_diff_target <= 0._DP) .AND. &
         (pub_cdft_group_spin_diff))  then
       if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
            &set cdft_group_spin_diff_target > 0. for a group-charge-&
            &difference-constrained  cDFT simulation, aborting'
       call comms_abort     
    end if

    ! jd: Sanity check for pbc_correction_cutoff
    if (pub_mt_cutoff /= 0.0_DP) then

       if (pub_mt_cutoff < 0.0_DP) then
          if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
               &pbc_correction_cutoff must be non-negative'
          call comms_abort
       end if

       if (pub_coulomb_cutoff) then
          if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
               &Cannot use cutoff Coulomb and pbc_correction_cutoff &
               &simultaneously'
          call comms_abort
       end if

       ! jd: When MT in use, cutoff Coulomb routine for Ion-ion energy is used
       pub_coulomb_cutoff_type = 'SPHERE'

    end if

    ! jd: Sanity check for dx_format_digits
    if (pub_dx_sig_digits < 1 ) then
       if (pub_on_root) write(stdout,'(a)') 'Error in get_rundat: &
            &number of significant digits specified with dx_format_digits &
            &must be positive.'
       call comms_abort
    end if

    ! jd: Ion-ion energy is calculated by direct summation iff
    !     (MT correction in use or cutoff Coulomb in use or smeared ions in use
    !      or implicit solvent in use or user overrode by pub_openbc_ion_ion.)
    pub_ii_energy_direct = &
         (pub_mt_cutoff /= 0.0_DP) .or. pub_coulomb_cutoff &
         .or. pub_is_smeared_ion_rep .or. pub_is_implicit_solvent &
         .or. pub_openbc_ion_ion

    ! jd: Hartree is calculated with multigrid iff
    !     (smeared ions in use or implicit solvent in use or user overrode by
    !     pub_openbc_hartree)
    pub_multigrid_hartree = &
         pub_is_smeared_ion_rep .or. pub_is_implicit_solvent &
         .or. pub_openbc_hartree

    ! jd: Local pseudo is calculated with open boundary conditions iff
    !     (smeared ions in use or implicit solvent in use or user overrode by
    !     pub_openbc_pspot)
    pub_open_localpseudo = &
         pub_is_smeared_ion_rep .or. pub_is_implicit_solvent &
         .or. pub_openbc_pspot

    ! jd: Force the cutoff type to be SPHERE, to get the usual Coulombic sum
    !     for ion-ion energy, except for vanilla cutoff Coulomb, where it's up
    !     to the user to decide
    if (pub_ii_energy_direct .and. .not. pub_coulomb_cutoff) then
       pub_coulomb_cutoff_type = 'SPHERE'
    end if

    ! jd: This flag tells cell_grid_distribute to ensure that the slab
    !     distribution is compatible with the multigrid solver. Set this flag
    !     to .true. anytime multigrid is used. Currently this is equivalent
    !     to using multigrid_hartree, but if there are other uses for the
    !     multigrid in the future (Simon?), this can be extended.
    pub_multigrid_in_use = pub_multigrid_hartree

    ! jd: Sanity check for implicit solvent keywords
    call utils_assert(pub_is_density_threshold > 0.0_DP, &
         'is_density_threshold must be positive')
    call utils_assert(pub_is_solvation_beta > 0.0_DP, &
         'is_solvation_beta must be positive')
    call utils_assert(pub_is_bulk_permittivity > 0.0_DP, &
         'is_bulk_permittivity must be positive')
    call utils_assert(pub_is_multigrid_error_tol > 0.0_DP, &
         'is_multigrid_error_tol must be positive')
    call utils_assert(pub_is_multigrid_defect_error_t > 0.0_DP, &
         'is_multigrid_defect_error_tol must be positive')
    call utils_assert(pub_is_smeared_ion_width > 0.0_DP, &
         'is_smeared_ion_width must be positive')
    call utils_assert(pub_is_solvation_method == 'DIRECT' .or. &
         pub_is_solvation_method == 'CORRECTIVE', &
         'is_solvation_method must be either ''direct'' or ''corrective''.')
    call utils_assert(pub_is_dielectric_model == 'FIX_INITIAL' .or. &
         pub_is_dielectric_model == 'GAUSSIAN_IONS' .or. &
         pub_is_dielectric_model == 'SELF_CONSISTENT', &
         'is_dielectric_model must be one of ''fix_initial'', &
         &''gaussian_ions &
         &'' or ''self_consistent''.')
    call utils_assert(pub_is_dielectric_function == 'FGF' .or. &
         pub_is_dielectric_function == 'GAUSSIAN', &
         'is_dielectric_function must be either ''fgf'' or ''gaussian''.')
    call utils_assert(pub_is_surface_thickness > 0.0_DP, &
         'is_surface_thickness must be positive')
    call utils_assert(pub_is_solvent_surface_tension > 0.0_DP, &
         'is_solvent_surface_tension must be positive')
    call utils_assert(pub_is_multigrid_max_iters > 0, &
         'is_multigrid_max_iters must be positive')
    call utils_assert(pub_is_multigrid_nlevels > 2, &
         'is_multigrid_nlevels must be at least 2')
    call utils_assert(pub_is_discretization_order > 0 .and. &
         mod(pub_is_discretization_order,2) == 0, &
         'is_discretization_order must be positive and even')
    call utils_assert(pub_is_bc_coarseness >= 0, &
         'is_bc_coarseness must be non-negative')
    call utils_assert(pub_is_bc_surface_coarseness > 0, &
         'is_bc_surface_coarseness must be positive')
    call utils_assert(pub_is_bc_threshold >= 0, &
         'is_bc_threshold must be non-negative')
    call utils_assert(.not. (pub_is_smeared_ion_rep .and. pub_coulomb_cutoff), &
         'is_smeared_ion_rep cannot be used simultaneously with &
         &coulomb_cutoff_type other than "NONE"')
    call utils_assert(.not. (pub_is_smeared_ion_rep .and. pub_mt_cutoff>0.0D0),&
         'is_smeared_ion_rep cannot be used simultaneously with &
         &pbc_correction_cutoff')
    call utils_assert(.not. (pub_is_include_cavitation .and. .not. &
         pub_is_implicit_solvent),'Cannot have is_include_cavitation without &
         &is_implicit_solvent')

    ! jd: Sanity check for open BC keywords
    call utils_assert(pub_openbc_pspot_finetune_nptsx > 1, &
         'openbc_pspot_finetune_nptsx must be >1')
    call utils_assert(pub_openbc_pspot_finetune_alpha >= 0.0_DP, &
         'openbc_pspot_finetune_alpha must be non-negative')
    call utils_assert(pub_openbc_pspot_finetune_f > 0 .or. &
         pub_openbc_pspot_finetune_f == -1, &
         'openbc_pspot_finetune_f must be positive or -1')

    ! jd: Warn the user that either PBC or OBC should be used, but not a mix
    ! ndmh: Warning was catching standard CC runs, so needed modified criteria
    inconsistent_bcs = (.not. pub_open_localpseudo .or. .not. &
         pub_multigrid_hartree) .and. (pub_open_localpseudo .or. &
         pub_multigrid_hartree)
    if (pub_open_localpseudo .and. pub_multigrid_hartree) then
       inconsistent_bcs = inconsistent_bcs .or. (.not. pub_ii_energy_direct)
    end if
    if (inconsistent_bcs .and. pub_on_root) then
       write(stdout,'(//a/)') 'WARNING: You are attempting a calculation with &
            &inconsistent boundary conditions (neither fully periodic, nor &
            &fully open). Make sure you know what you are doing.'
    end if

    ! jd: Advise agains using is_bc_surface_coarseness for charged systems
    if (pub_is_bc_surface_coarseness > 1 .and. pub_charge /= 0.0_DP &
         .and. pub_on_root) then
       write(stdout,'(//a)') 'WARNING: is_bc_surface_coarseness is larger than&
            & 1. Values larger than 1'
       write(stdout,'(a)') '         should only be used for charge-neutral&
            & systems -- for charged systems'
       write(stdout,'(a)') '         the accuracy of the energy&
            & calculation will likely be negatively'
       write(stdout,'(a)') '         impacted. Seriously&
            & consider setting is_bc_surface_coarseness to 1'
       write(stdout,'(a)') '         and if you find&
            & too much time is spent in multigrid_prepare_bound_cond,'
       write(stdout,'(a)') '         it''s&
            & better to (further) increase is_bc_coarseness rather than'
       write(stdout,'(a/)') '         is_bc_surface_coarseness.'
    end if

    ! ndmh: if ngwf_cg_max_step is left at its default, rescale for varying ecut
    ! ndmh: 2.0 is about right at 600eV - higher cutoff requires higher max step
    if (ngwf_cg_max_step < 0.0_DP) then
       ngwf_cg_max_step = -ngwf_cg_max_step * (cutoff_energy / 22.04959837_DP)
    end if


  end subroutine get_rundat

  subroutine rundat_exit()

    !=======================================================================!
    ! This subroutine deallocates those the allocatable module variables in !
    ! the rundat module (if required).                                      !
    !-----------------------------------------------------------------------!
    ! Written by Quintin Hill on 27/05/2009.                                !
    !=======================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    integer :: ierr

    if (allocated(pub_ldos_groups)) then
       deallocate(pub_ldos_groups,stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_ldos_groups',ierr)
    end if
    if (allocated(pub_ldos_group_nsp)) then
       deallocate(pub_ldos_group_nsp,stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_ldos_group_nsp',ierr)
    end if
    if (allocated(pub_bs_kpoint_path_start)) then
       deallocate(pub_bs_kpoint_path_start, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_bs_kpoint_path_start',ierr)
    end if
    if (allocated(pub_bs_kpoint_path_end)) then
       deallocate(pub_bs_kpoint_path_end, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_bs_kpoint_path_end',ierr)
    end if

    ! lpl: NBO stuff
    if (allocated(pub_nbo_write_species)) then
       deallocate(pub_nbo_write_species, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_nbo_write_species',ierr)
    end if
    if (allocated(pub_nbo_ngwf_label)) then
       deallocate(pub_nbo_ngwf_label, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_nbo_ngwf_label',ierr)
    end if
    if (allocated(pub_nbo_list_plotnbo)) then
       deallocate(pub_nbo_list_plotnbo, stat=ierr)
       call utils_dealloc_check('rundat_exit','pub_nbo_list_plotnbo',ierr)
    end if
    ! lpl: NBO stuff

  end subroutine rundat_exit


end module rundat
