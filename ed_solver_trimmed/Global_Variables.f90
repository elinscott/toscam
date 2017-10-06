module globalvar_ed_solver

  USE genvar
  USE masked_matrix_class_mod
  use mpirout
  use openmpmod 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
#ifdef ALLFIRSTCALL
 LOGICAL,parameter                                 :: ALL_FIRST_CALL=.true.
#else
 LOGICAL,parameter                                 :: ALL_FIRST_CALL=.false.
#endif
#ifdef debug
  LOGICAL, parameter                               :: verboseall=.true.
#else
  LOGICAL, parameter                               :: verboseall=.false.
#endif
  LOGICAL                                          :: fmos,fmos_fluc,fmos_hub1
  integer                                          :: fmos_iter,nvecout
  REAL(8)                                          :: fmos_mix
  INTEGER                                          :: flag_idelta_two_scales_ed
  REAL(8),allocatable                              :: UCCr(:,:,:)
  COMPLEX(8),allocatable                           :: UCCc(:,:,:)
  INTEGER,allocatable                              :: UCC(:,:,:,:)
  INTEGER                                          :: ncpt_chain_coup
  REAL(8)                                          :: quench_mag,quench_U
  LOGICAL                                          :: keldysh_pert_ground_sector,flag_slater_int,use_precomputed_slater_matrix,cpt_correct_green_out,ncpt_flag_two_step_fit,gen_cpt
  COMPLEX(8),allocatable                           :: Slater_Coulomb_c(:,:,:,:)
  REAL(8),allocatable                              :: Slater_Coulomb_r(:,:,:,:)
  INTEGER                                          :: do_quench,do_quench_,quench_orb,ncpt_approx,ncpt_tot,ncpt_approx_tot,ncpt_para
  REAL(8)                                          :: cpt_upper_bound,cpt_lagrange
  COMPLEX(8),allocatable                           :: epsilon_cpt(:,:),T_cpt(:,:)
  LOGICAL                                          :: DO_NOT_USE_OPT_LANCZOS,do_keldysh_gbigger,quench_cancel_statistics
  LOGICAL                                          :: donot_compute_holepart_spm
  INTEGER                                          :: keldysh_ortho
  REAL(8)                                          :: keldysh_t0,keldysh_tmax,keldysh_delta
  INTEGER                                          :: keldysh_n
  INTEGER                                          :: freeze_poles_delta_iter,freeze_poles_delta_niter,Niter_search_max_0
  REAL(8)                                          :: rdelta_frequ_eta1,rdelta_frequ_T,rdelta_frequ_w0,Jhund,rdelta_frequ_eta2
  REAL(8),allocatable                              :: UUmatrix(:,:),JJmatrix(:,:)
  logical                                          :: use_input_Eimp_matrices,fit_all_elements_show_graphs
  LOGICAL                                          :: do_keldysh,FLAG_DUMP_INFO_FOR_GAMMA_VERTEX
  LOGICAL                                          :: FLAG_NCUP,fit_green_func_and_not_hybrid
  INTEGER                                          :: iterdmft,niterdmft,iwindow
  INTEGER                                          :: cuda_blocksize 
  INTEGER,allocatable                              :: MASK_AVERAGE(:,:),MASK_AVERAGE_SIGMA(:,:)
  LOGICAL                                          :: freeze_poles_delta,FLAG_GUP_IS_GDN,always_compute_static_obs
  INTEGER                                          :: FLAG_MPI_GREENS
  REAL(DBL)                                        :: freeze_pole_lambda,beta_ED,beta_ED_,lambda_sym_fit
  REAL(DBL)                                        :: weight_expo,cutoff_dynamic,cutoff_min_lanczos_vec
  REAL(DBL),allocatable                            :: FROZEN_POLES(:)
  CHARACTER(20)                                    :: which_lanczos
  LOGICAL                                          :: OPEN_MP,ON_FLY,USE_TRANSPOSE_TRICK_MPI,restarted,only_compute_density,skip_fit_  
  LOGICAL                                          :: start_para,PAIRING_IMP_TO_BATH,use_specific_set_parameters,force_no_bcs_pairing
  LOGICAL                                          :: superconducting_state
  LOGICAL,parameter                                :: Jhund_Slater_type=.true.
  REAL(DBL),ALLOCATABLE                            :: param_input(:),param_output(:)
  INTEGER                                          :: average_G
  LOGICAL                                          :: diag_bath,diag_V,bath_nearest_hop,FLAG_ALL_GREEN_FUNC_COMPUTED
  LOGICAL                                          :: FLAG_BUILD_CORREL_LOW_PART,dump_ground_state
  LOGICAL                                          :: use_cuda_lanczos=.false.
  INTEGER                                          :: szmin,szmax,nupmin,nupmax,ndnmin,ndnmax
  LOGICAL                                          :: track_sectors=.false.,FLAG_FULL_ED_GREEN=.false.
  LOGICAL,PARAMETER                                :: USE_CC=.true.
  LOGICAL                                          :: SCAN_FULL_NUP_NDN
  INTEGER                                          :: Niter_search_max,min_all_bath_param,fit_nw,istati
  logical                                          :: supersc_state
  REAL(DBL)                                        :: tot_repulsion
  REAL(DBL)                                        :: dist_max,search_step,fit_weight_power,cutoff_rvb,cutoff_hamilt_param
  COMPLEX(DBL)                                     :: energy_global_shift,energy_global_shift2
  CHARACTER(20)                                    :: FIT_METH
  INTEGER                                          :: Nitergreenmax
  INTEGER                                          :: Nitermax                    ! max # of Lanczos iterations
  LOGICAL                                          :: force_pairing_from_mask
  LOGICAL                                          :: para_state,donot_compute_holepart,force_sz_basis,force_nupdn_basis
  LOGICAL                                          :: fast_fit,first_iter_use_edinput,force_no_pairing,force_para_state,force_singlet_state
  REAL(DBL)                                        :: tolerance                   ! Lanczos tolerance
  REAL(DBL)                                        :: dEmax=0.d0,dEmax0=0.d0      ! max.energy of excited states to consider
  INTEGER                                          :: Neigen            = 0       ! max. # of eigenvalues computed in this window
  INTEGER                                          :: Block_size        = 0       ! Block size (ignored if Neigen=1: Block_size=1)
  INTEGER                                          :: nsec  = 1                   ! Number of sectors to scan 
  INTEGER                                          :: nsec0 = 1                   ! Scan of the sectors start at sec=sec0 
  INTEGER                                          :: window_hybrid = 0           ! window of matsubara frequencies to fit
  INTEGER                                          :: window_hybrid2 = 0
! window of matsubara frequencies to fit
  INTEGER                                          :: window_weight = 1           ! ratio weight point inside window / point outside window
  REAL(DBL), ALLOCATABLE                           :: dens(:)                     ! onsite charge density
  INTEGER, ALLOCATABLE                             :: list_sectors(:)
  LOGICAL                                          :: verbose_graph               ! plot fits of the hybridization at each step of the minimization
  REAL(DBL)                                        :: fit_shift 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

end module
