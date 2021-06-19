!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

module fourier_transform_dmft_one_iter

   use linalg
   use specialFunction
   use genvar, only: DP

   PRIVATE
   PUBLIC :: densfourier

   logical :: green_
   real(kind=DP), parameter :: ddsign = -1.0

contains

!***********************************************!
!***********************************************!
!***********************************************!
!***********************************************!
!***********************************************!
!***********************************************!
!***********************************************!
!***********************************************!

   real(kind=DP) function densfourier(nw, beta, omi, green)
      implicit none
      integer           ::  i, j, k, nw
      integer, parameter ::  Ntau = 1
      real(kind=DP)           ::  beta, tau(Ntau), omi(nw)
      complex(kind=DP)        ::  Gtau(Ntau)
      complex(kind=DP)        ::  green(1:nw), greensave(1:nw)
      greensave = green
      tau = (/(ddsign*0.00001d0, i=1, 1)/)
      call Fourier_(Ntau, nw, greensave, omi, Gtau, tau, beta, green=.true.)
      densfourier = Gtau(1)
      return
   end function

!***********************************************!
!***********************************************!
!***********************************************!
!***********************************************!

   subroutine add_substract_tail(nw, Giom, iom, ahh, ddd, add)
      use genvar, only: imi
      implicit none
      integer          :: nw
      complex(kind=DP)       :: Giom(1:nw), temp(1:nw)
      real(kind=DP)          :: iom(1:nw), ah, ahh, drsign, ddd
      integer          :: t, n, i, j
      logical          :: add

      if (add) then
         drsign = 1.d0
      else
         drsign = -1.d0
      endif
      temp = 0.d0
      do n = 1, nw
         if (n /= 0) then
            temp(n) = Giom(n) + ddd*drsign/(imi*iom(n) + ahh)
         endif
      enddo

      Giom = temp

   end subroutine

!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/

   subroutine get_tail_parameters(nw, Giom, iom, ahh, ddd)
      implicit none
      integer          :: nw, kk
      complex(kind=DP)       :: Giom(1:nw)
      real(kind=DP)          :: iom(1:nw), ddd, ahh
      ddd = get_ddd(nw, iom, Giom)
      ahh = get_ah(nw, iom, Giom)
      ahh = ahh/ddd
      return
   end subroutine

!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/

   subroutine Fourier_(Ntau, nw, Giom, iom, Gtau_, tau, beta, green)
      implicit none
      integer    :: Ntau, nw, kk
      complex(kind=DP) :: Giom(1:nw), Giomback(1:nw)
      real(kind=DP)    :: iom(1:nw), tau(Ntau), ah, ahh, mm
      complex(kind=DP) :: Gtau_(Ntau)
      real(kind=DP)    :: beta, df, temp
      complex(kind=DP) :: csum
      real(kind=DP)    :: bb, ddd, ahh2, ddd2, g1, gn
      integer    :: t, n, i, j, itail
      real(16)   :: t1, t2, pp1, pp2
      logical    :: green
      integer    :: ifirst, ifac

      green_ = green
      ifirst = 1
      ifac = 2

      Giomback = Giom

      call get_tail_parameters(nw, Giom, iom, ahh, ddd)
      call add_substract_tail(nw, Giom, iom, ahh, ddd, add=.false.)

      do t = 1, Ntau
         csum = 0.d0
         do n = ifirst, nw
            if (n /= 0) csum = csum + MPLX(ddsign*iom(n)*tau(t))*Giom(n)
         enddo
         Gtau_(t) = csum/beta*dble(ifac)
      enddo

      do t = 1, Ntau
         if (ahh > 0.d0) then
            pp1 = DEXPc(ahh*(tau(t) - beta))
            pp1 = pp1/(DEXPc(-beta*ahh) + 1.d0)
            pp2 = DEXPc(ahh*(tau(t)))
            pp2 = pp2/(DEXPc(-beta*ahh) + 1.d0)
         else
            pp1 = DEXPc(ahh*tau(t))
            pp1 = pp1/(DEXPc(beta*ahh) + 1.d0)
            pp2 = DEXPc(ahh*(tau(t) + beta))
            pp2 = pp2/(DEXPc(beta*ahh) + 1.d0)
         endif
         Gtau_(t) = Gtau_(t) - ddd*(step_func_(tau(t))*pp1 - step_func_(-tau(t))*pp2)
      enddo

      Giom = Giomback

      return
   end subroutine

!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/

   real(kind=DP) function get_ddd(nw, iom, Giom)
      implicit none
      integer           :: nw
      real(kind=DP)           :: iom(1:nw), ddd
      complex(kind=DP)        :: Giom(1:nw)
      if (green_) then
         get_ddd = 1.d0
      else
         get_ddd = -aimag(Giom(nw))*iom(nw)
      endif
   end function

!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/

   real(kind=DP) function get_ah(nw, iom, Giom)
      implicit none
      integer    :: nw
      real(kind=DP)    :: iom(1:nw), ddd
      complex(kind=DP) :: Giom(1:nw)
      integer    :: i, kk, p, ii, k
      get_ah = Real(Giom(nw))*(iom(nw)**2)
   end function

!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/
!************************************************/

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
!==============================================================!

module onetep_variables

   use fourier_transform_dmft_one_iter ! local module

   use linalg, only: minloci, MPLX
   use namelistmod, only: namelist_set, namelist_init, putel_in_namelist, &
                              & look_for_namelist_in_file, look_for_command_line_argument
   use genvar, only: imi, pi, DP
   use matrix, only: diag, rearrange_columns_to_identity, QR_decomp, eigenvector_matrix_c_
   use StringManip, only: toString, StrInt2
   use strings, only: replace_in_string, string, assignment(=)
   use matrix, only: offdiag, Id, invmat, write_array, eigenvector_matrix
   use specialFunction, only: SlaterF, cmp_all_Gaunt

   implicit none
   logical           ::   real_sproj_from_onetep, real_eimp_from_tail, real_eimp_from_onetep, real_smat_from_onetep
   integer           ::   flag_use_broadening_two_scales_ed
   logical           ::   average_green_ed, verysilent
   integer           ::   fmos_iter, fmos_hub_iter_mu
logical           ::   fmos_use,fmos_fluc,fmos_hub1,sum_over_k_dft,keep_both_real_and_matsu_last,force_self_infty_real,double_counting_zero_self_from_matsu
   real(kind=DP)           ::   fmos_mix, fmos_hub_range_mu
   real(kind=DP)           ::   cutoff_simp_offdiag
   logical           ::   hf_average, hf_krauth, hf_hartree, ed_no_real_overide
   integer, parameter ::   nspec_frequ = 7   !The last nspec_frequ frequencies are for testing puroposes
   integer           ::   iwindow, nproc_mpi_solver, cubic, im_solver, real_solver, max_steps, followPeak, nsec, nsec0, ntauhf
real(kind=DP)           ::   slater_inter_cutoff,amp_slight_sym_breaking,QQ,spin_orbit,UU,Jhund,StartLambda,EndLambda,alpha,dist_max
   logical           ::   real_smat_coef
logical           ::   ncpt_two_step_all_iter,rotation_scheme_read_write,double_counting_zero_self,double_counting_zero_self_av,ncpt_two_step,assume_proj_overlap_is_diagonal,amp_slight_sym_breaking_all_iter,use_precomputed_slater_matrix
   real(kind=DP)           ::   fit_weight_power, double_counting_nf, mixing, ed_rdelta, ed_frequ_min, ed_frequ_max
   real(kind=DP)           ::   ed_rdelta_frequ_eta1, ed_rdelta_frequ_T, ed_rdelta_frequ_w0
   integer           ::   rotation_scheme, nmatsu_ed, nmatsu_ctqmc, nmatsu_long, mcs_ctqmc, openmp_solver
   integer           ::   n_frequ, ed_real_frequ, ed_real_frequ_last, paramagnetic_ed, sites_ed, niter_dmft, fit_nw
logical           ::   rotation_scheme_pm,dmft_for_dimer,rotation_green_function,cluster_dmft_green_for_self_consistence,last_iter_is_real,ed_solver_compute_all_green_functions
   logical           ::   ed_star_geom, use_custom_command_for_atomd, compute_ed_spin_correlation
logical           ::   ctqmc_erase_status,use_input_delta,ed_do_not_keep_previous_fit_param,read_cix_file,rotate_int_after_earlier_transfo,second_order_correction_to_eimp
   logical           ::   ed_compute_all, flag_symmetrize_green, flag_blank_out_green_offdiag_for_testing
   logical           ::   flag_blank_out_sigma_offdiag_for_testing, flag_use_jj_basis_for_f, flag_correct_eimp_spin_orbit
logical           ::   double_counting_with_no_average,dimer_average_occupation,nohole,fit_green,impose_no_cdw,ed_compute_retarded_every_step,use_eimp_from_onetep
 logical           ::   rotate_ortho_av_renorm_int, rotate_to_ortho_basis, use_simp_from_onetep, use_eimp_from_onetep_with_sigma_cor
   logical           ::   flag_get_t2c_real, flag_use_slater_in_ed
   logical           ::   use_simp_from_onetep_for_ortho
   character(200)    ::   atom_d_command
   character(30)     ::   which_lanczos
   integer           ::   ed_nsearch, Neigen, Block_size, lowest_occupancy, highest_occupancy
   integer           ::   ncpt_approx
   real(kind=DP)           ::   cpt_upper_bound, cpt_lagrange, tail_linear_scaling
   complex(kind=DP)        ::   ccctemp, ccctemp_av
   logical           :: print_qc

contains

   subroutine init_my_input_variables
      type(namelist_set) :: nm
      logical            :: checkrot
!#######################################!
!#######################################!
!#######################################!
      call namelist_init(nm, 200, name_of_namelist='dmftonetep_variables')
      call putel_in_namelist(nm, cubic, 'cubic', 1, 'Definition:=1 for cubic harmonics, =2 for spherical harmonics')
      call putel_in_namelist(nm, spin_orbit, 'spin_orbit', 0.d0, 'Definition:spin-orbit coupling')
      call putel_in_namelist(nm, UU, 'UU', 5.d0, 'Definition:Coulomb repulsion')
      call putel_in_namelist(nm, Jhund, 'Jhund', 0.d0, 'Definition:Hunds coupling')
 call putel_in_namelist(nm, im_solver     ,'im_solver'      ,4,        'Definition:choice of the solver for matsubara calculations, =1 for CTQMC, =4 for ED, =5 for HF ')
 call putel_in_namelist(nm, dmft_for_dimer,'dmft_for_dimer' ,.false.,  'Definition:if true the library treats the problem as a dimer')
 call putel_in_namelist(nm, real_solver   ,'real_solver'    ,2,        'Definition:choice of the solver for real-axis calculations, =2 for NCA, =3 for OCA')
      call putel_in_namelist(nm, StartLambda, 'StartLambda', -300.d0, 'Definition:NCA/OCA StartLambda')
      call putel_in_namelist(nm, EndLambda, 'EndLambda  ', 200.d0, 'Definition:NCA/OCA EndLambda')
      call putel_in_namelist(nm, alpha, 'alpha      ', 0.3d0, 'Definition:NCA/OCA mixing for internal self consistence ')
      call putel_in_namelist(nm, mixing, 'mixing     ', 0.5d0, 'Definition:DMFT mixing for Self-Energy ')
      call putel_in_namelist(nm, max_steps, 'max_steps  ', 30, 'Definition:NCA/OCA max numer of internal steps')
      call putel_in_namelist(nm, followPeak, 'followPeak ', -1, 'Definition:NCA/OCA follow peak')
      call putel_in_namelist(nm, QQ, 'QQ', 8.d0, 'Definition:NCA/OCA QQ parameter for seeking out lambda')
      call putel_in_namelist(nm, ed_real_frequ, 'ed_real_frequ', 1000, 'Definition:ed number of frequencies')
 call putel_in_namelist(nm, ed_real_frequ_last ,'ed_real_frequ_last', 100     , 'Definition:ed number of frequencies for last DMFT iter, going back to onetep for FULL_DOS')
      call putel_in_namelist(nm, ed_frequ_min, 'ed_frequ_min ', -10.d0, 'Definition: min ED real frequ')
      call putel_in_namelist(nm, ed_frequ_max, 'ed_frequ_max ', +10.d0, 'Definition: max ED real frequ')
     call putel_in_namelist(nm, ed_rdelta, 'ed_rdelta', +0.0001d0, 'Definition: small i*delta to shift off the real frequency axis')
      call putel_in_namelist(nm, ed_rdelta_frequ_eta1, 'ed_rdelta_frequ_eta1', +0.002d0, 'Definition: eta1, rdelta far from 0')
      call putel_in_namelist(nm, ed_rdelta_frequ_T, 'ed_rdelta_frequ_T', +0.0003d0, 'Definition: ramp to move from 0 to eta1')
      call putel_in_namelist(nm, ed_rdelta_frequ_w0, 'ed_rdelta_frequ_w0', +0.0006d0, 'Definition: frequency at which we have eta1')
      call putel_in_namelist(nm, sites_ed, 'sites_ed', 5, 'Definition: number of sites in the bath for ED solver')
      call putel_in_namelist(nm, niter_dmft, 'niter_dmft', 5, 'Definition:Number of DMFT iteration')
 call putel_in_namelist(nm, fit_nw,        'fit_nw',        0        , 'Definition:Number of matsubara frequency to fit for ED solver, 0=fits all frequencies')
 call putel_in_namelist(nm, double_counting_nf, 'double_counting_nf', -1.d0, 'Definition:if negative, double counting is computed from orbital occupation, if positive double counting is obtained from positive input double_counting_nf')
     call putel_in_namelist(nm, openmp_solver, 'openmp_solver', 8, 'Definition:Number of open-mp cores running for the dmft-solver')
 call putel_in_namelist(nm, rotation_green_function,  'rotation_green_function', .false. , 'Definition: if true the code will rotate the Green Function such that it is diagonal')
 call putel_in_namelist(nm, rotate_int_after_earlier_transfo, 'rotate_int_after_earlier_transfo', .false. ,  'Definition: if true will rotate the CTQMC interaction according to an earlier transformation, when the file Trans.loc.{pm,spin} are present')
 call putel_in_namelist(nm, cluster_dmft_green_for_self_consistence, 'cluster_dmft_green_for_self_consistence',.true. , 'Definition: if true keeps the off diagonal self energy terms when sending back to onetep')
 call putel_in_namelist(nm, fit_weight_power,  'fit_weight_power', 0.5d0, 'Definition: weight use for the fitting distance in ED solver')
 call putel_in_namelist(nm, last_iter_is_real, 'last_iter_is_real', .true. , 'Definition: if true and if solver is ED, last iteration dumps out real frequencies in sigma_output, such that FULL_DOS can be computed with onetep.System')
      call putel_in_namelist(nm, nmatsu_ed, 'nmatsu_ed', 4000, 'Definition: number of matsubara frequencies for ED solver')
      call putel_in_namelist(nm, nmatsu_ctqmc, 'nmatsu_ctqmc', 20, 'Definition: number of matsubara frequencies for CTQMC solver')
 call putel_in_namelist(nm, nmatsu_long, 'nmatsu_long',    0,  'Definition: if 0 does not have any effect, if >0 then it is the number of additional matsubara frequencies generated following the tail expansion for CTQMC solver')
      call putel_in_namelist(nm, mcs_ctqmc, 'mcs_ctqmc', 20000000, 'Definition: number of MonteCarlo steps for CTQMC solver')
 call putel_in_namelist(nm, rotation_scheme,   'rotation_scheme',     1,   'Definition: if 1 rotates the green function according to the occupation matrix, if =2 it rotates the green function according to its max off diagonal elements, 3:rotates to diagonalize H, 4:same as 2 but takes the real part of the rotation matrix only')
      call putel_in_namelist(nm, nproc_mpi_solver, 'nproc_mpi_solver', 1, 'Definition: number of MPI cpus for solver')
 call putel_in_namelist(nm, second_order_correction_to_eimp, 'second_order_correction_to_eimp', .true. , 'Definition: if true includes second order corrections to Eimp, especially large if nfrequ is small and T is small')
 call putel_in_namelist(nm, read_cix_file,     'read_cix_file', .false.,   'Definition: if true will read the cix file instead of generating it')
 call putel_in_namelist(nm, ed_do_not_keep_previous_fit_param, 'ed_do_not_keep_previous_fit_param', .false., 'Definition: if true the previous ED fit will not be used for the next iteration')
 call putel_in_namelist(nm, atom_d_command,'atom_d_command', 'atom_d.py',   'Definition: name of the command to run atom_d.py, note: can be useful to run atom_d.py on a different node if scipy and numpy not installed, in that case you can use a mpi command line , e.g. mpiexec -host XX -np 1 atom_d.py ')
 call putel_in_namelist(nm, which_lanczos ,'which_lanczos', 'NORMAL' ,      'Definition : FULL_ED,NORMAL,GPU (no hunds coupling implemented),ARPACK') 
 call putel_in_namelist(nm, Neigen        ,'Neigen',           10    ,      'Definition : for finite temperature Lanczos, number of max converged lowest eigenvalues')   
 call putel_in_namelist(nm, Block_size    ,'Block_size',        0    ,      'Definition : if =0 will not use Block Lanczos, -1 use Block-scheme Cullum, if +i will use Block Lanczos library with i Lanczos vectors, typical value is +5')  
 call putel_in_namelist(nm,use_custom_command_for_atomd,'use_custom_command_for_atomd', .false., 'Definition : if true uses custom command line for atom_d.py instead of the default, custom command is defined in variable atom_d_command')
 call putel_in_namelist(nm,use_input_delta, 'use_input_delta', .false.,     'Definition : if true will use the input hybridization to compute the self energy in the ED solver, instead of the fitted hybridization, this can beak causality when the fit is not good enough, but can help capturing more elements of the lattice physics')
 call putel_in_namelist(nm,impose_no_cdw, 'impose_no_cdw', .false.,         'Definition : if true will impose that the self energy in the dimer case is not a charge density wave (symmetrize sigma at the output of ED)')
 call putel_in_namelist(nm,fit_green,     'fit_green',     .false.,         'Definition : if true ED solver will fit the Weiss Field instead of the hybridization, this can have an impact especially when some of the orbitals are weakly correlated, then it is better to fit the Weiss Field (flag=.true.)')
 call putel_in_namelist(nm,compute_ed_spin_correlation,'compute_ed_spin_correlation', .false. , 'Definition: if true will compute spin-spin correlations within ED, useful for dimers especially')
 call putel_in_namelist(nm,ed_compute_retarded_every_step,'ed_compute_retarded_every_step', .false., 'Definition: if true will compute the retarded Green function at every dmft step, it slows down the code by a factor of two, if false the retarded GF is computed only at the last DMFT iteration')
 call putel_in_namelist(nm,nohole,'nohole',.false.,'Definition: if true will not compute the hole part of the green function, it is a safe assumption, but always good to make a check about that')
 call putel_in_namelist(nm,ed_star_geom,'ed_star_geom',.false.,'Definition: if true will use the so-called star geometry for the hybridization for ED solver')
 call putel_in_namelist(nm,ed_nsearch,'ed_nsearch',200000,'Definition: number of conjugate gradient iterations for fitting the hybridization for ED solver')
 call putel_in_namelist(nm,dimer_average_occupation,'dimer_average_occupation',.true.,'Definition: if true the occupation will be averaged on both sites of the dimer for the double counting correction')
 call putel_in_namelist(nm,lowest_occupancy , 'lowest_occupancy'  ,2 ,'Definition: lowest  occupancy of the atom, used for CTQMC for the L=3 case')
 call putel_in_namelist(nm,highest_occupancy, 'highest_occupancy' ,12,'Definition: highest occupancy of the atom, used for CTQMC for the L=3 case')
 call putel_in_namelist(nm,flag_correct_eimp_spin_orbit,'flag_correct_eimp_spin_orbit',.true.,'Definition: if true will add a correction to the impurity level due to the spin-orbit coupling')
 call putel_in_namelist(nm,flag_use_jj_basis_for_f,'flag_use_jj_basis_for_f',.false.,'Definition: if true will use the L+S j-j basis for the f orbital, needed for f orbitals for the CTQMC calculation')
 call putel_in_namelist(nm,flag_blank_out_sigma_offdiag_for_testing,'flag_blank_out_sigma_offdiag_for_testing',.false.,'Definition: if true will blank out the off-diagonal elements of the cluster self energy, so the diagonal single site dmft should be recovered')
 call putel_in_namelist(nm,flag_blank_out_green_offdiag_for_testing,'flag_blank_out_green_offdiag_for_testing',.false.,'Definition: if true will blank out the off-diagonal elements of the cluster hybridization, so the diagonal single site dmft should be recovered')
 call putel_in_namelist(nm,flag_symmetrize_green,'flag_symmetrize_green',.false.,'Definition: if true will symmetrize input green function')
 call putel_in_namelist(nm,amp_slight_sym_breaking,'amp_slight_sym_breaking',0.d0,'Definition: if finite it will induce a small symmetry breaking will be used for the polarized calculations in the impurity levels')
 call putel_in_namelist(nm,ctqmc_erase_status,'ctqmc_erase_status',.false.,'Definition: if true will erase status file after each CTQMC iteration, so there is no correlations between CTQMC iterations')
 call putel_in_namelist(nm, ed_compute_all, 'ed_compute_all', .false., 'Definition: if true will compute all observables within ED')
 call putel_in_namelist(nm,rotation_scheme_pm,'rotation_scheme_pm',.true.,'Definition: if true will use the same rotation for spin up and spin dn, it should be the case to respect the symmetry of the interacting Hamiltonian')
 call putel_in_namelist(nm,amp_slight_sym_breaking_all_iter,'amp_slight_sym_breaking_all_iter',.false.,'Definition: if true the small symmetry breaking amp_slight_sym_breaking will be used at all DMFT iteration and not only at the first iteration')
 call putel_in_namelist(nm,double_counting_with_no_average,'double_counting_with_no_average',.false.,'Definition: if true will not use the averaged U value which takes into accoung the hunds coupling rule, but instead just the screend U')
 call putel_in_namelist(nm,use_eimp_from_onetep,'use_eimp_from_onetep',.false.,'Definition:if true will use the analytical formula for eimp obtained from Onetep')
 call putel_in_namelist(nm,use_simp_from_onetep,'use_simp_from_onetep',.false.,'Definition:if true will use the analytical formula for simp obtained from Onetep')
 call putel_in_namelist(nm,use_eimp_from_onetep_with_sigma_cor,'use_eimp_from_onetep_with_sigma_cor',.false.,'Definition:if true, and if using use_eimp_from_onetep , then it wil use the formula himp + Proj Self(oo) Proj - Self_imp (oo), in principle the corrections induced by Sigma are by construction void, but it can be checked in the output file')
 call putel_in_namelist(nm,rotate_to_ortho_basis,'rotate_to_ortho_basis',.false.,'Definition: if true the code will perform a rotation to the orthogonal basis back and forth, please use the corresponding script to check the convergence')
 call putel_in_namelist(nm,rotate_ortho_av_renorm_int,'rotate_ortho_av_renorm_int',.false.,'Definition: if true will use averaged renormalization induced by the non orthogonal basis set, if false will use an orbital resolved renormalization')
 call putel_in_namelist(nm,flag_get_t2c_real,'flag_get_t2c_real',.false.,'Definition: if true will try to correct the T2C matrix for the rotation with the right phases such that it reduces the imaginary part as much as possible')
 call putel_in_namelist(nm,flag_use_slater_in_ed,'flag_use_slater_in_ed',.false.,'Definition: if true uses Slater interaction for ED')
 call putel_in_namelist(nm,slater_inter_cutoff,'slater_inter_cutoff',1.d-5,'Definition: the cutoff under which the Slater terms are not taken into account in the Hamiltonian matrix')
 call putel_in_namelist(nm,use_precomputed_slater_matrix, 'use_precomputed_slater_matrix', .false., 'Definition: if true will use a precomputed Hamiltonian matrix for the Slater interaction')
 call putel_in_namelist(nm,assume_proj_overlap_is_diagonal,'assume_proj_overlap_is_diagonal',.false.,'Definition: if true the code will simplify the calculations by assuming that the projected overlap matrix is diagonal')
      call putel_in_namelist(nm, ncpt_approx, 'ncpt_approx', 0, 'Definition : number of sites for CPT approximation in ED solver')
      call putel_in_namelist(nm, cpt_upper_bound, 'cpt_upper_bound', 0.d0, 'Definition : max value for CPT V parameter')
      call putel_in_namelist(nm, cpt_lagrange, 'cpt_lagrange', 0.d0, 'Definition : lagrange parameter to minimize V parameter')
 call putel_in_namelist(nm,  ncpt_two_step,   'ncpt_two_step',  .true., 'Definition : if true always do the fits in two steps procedures')
 call putel_in_namelist(nm,  double_counting_zero_self, 'double_counting_zero_self', .false., 'Definition: if true it will impose that the self energy is zero at infinite frequency, basically it is another way to substract the double counting correction, note: this should be used if the chemical potential is shifted within onetep with pub_dmft_shift_pot variable, since the occupancy matrix is not recomputed accordingly')
 call putel_in_namelist(nm, double_counting_zero_self_av,'double_counting_zero_self_av',.false., 'Definition: if true remove the orbital averaged Sigma(w=oo), and not the individual components')
      call putel_in_namelist(nm, iwindow, 'iwindow', 1, 'Definition: the window of quantum eigen sector that ED is tracking')
 call putel_in_namelist(nm,  rotation_scheme_read_write, 'rotation_scheme_read_write', .false., 'Definition: if true and rotation==4, it will after first dmft iteration write the rotation in file mask_dmft_rot, and use the same rotation for the next DMFT iteration')
call putel_in_namelist(nm, nsec0, 'nsec0', 0, 'Definition: first sector to scan for ED solver, then goes upward with nsec0 sectors')
      call putel_in_namelist(nm, dist_max, 'dist_max', 1.d-10, 'Definition: max distance for ED fit')
      call putel_in_namelist(nm, nsec, 'nsec', -1, 'Definition: list of sector to scan, if -1 scan all sectors')
 call putel_in_namelist(nm,  ncpt_two_step_all_iter, 'ncpt_two_step_all_iter', .false., 'Definition: if true will enforce the two step fit for CPT at every dmft iteration, only effective if ncpt_two_step is true')
      call putel_in_namelist(nm, ntauhf, 'ntauhf', 50, 'Definition: imaginary time discretization for HF solver')
      call putel_in_namelist(nm, fmos_use, 'fmos_use', .false., 'Definition: if true use FMOS solver instead of ED')
      call putel_in_namelist(nm, fmos_iter, 'fmos_iter', 10, 'Definition: number of FMOS internal iterations')
      call putel_in_namelist(nm, fmos_mix, 'fmos_mix', 0.4d0, 'Definition: mixing for FMOS internal iterations')
call putel_in_namelist(nm, fmos_fluc, 'fmos_fluc', .false., 'Definition: if true introduces charge flucutations in the FMOS solver')
      call putel_in_namelist(nm, hf_krauth, 'hf_krauth', .false., 'Definition: Hirsch Fye with Krauth')
      call putel_in_namelist(nm, hf_average, 'hf_average', .false., 'Definition: Hirsch Fye average over time slices')
 call putel_in_namelist(nm, hf_hartree,'hf_hartree',.true.,'Definition: if true will impose the hartree shift in the self energy obtained by Quantum Monte Carlo, the high energy component is indeed hard to obtain by brute force calculations')
 call putel_in_namelist(nm, fmos_hub1,'fmos_hub1',.false.,'Definition: if true will use the hubbard 1 solver as a fast multi orbital solver, only compatible with the Slater interaction vertex')
 call putel_in_namelist(nm, fmos_hub_iter_mu,'fmos_hub_iter_mu',1,'Definition: number of iteration to impose a target density in Hubbard I solver')
 call putel_in_namelist(nm, fmos_hub_range_mu,'fmos_hub_range_mu',0.001d0,'Definition: interval to adjust the chemical potential to reach a target density in the Hubbard I solver')
 call putel_in_namelist(nm, cutoff_simp_offdiag,'cutoff_simp_offdiag',0.01d0,'Definition: if off diag elements of simp smaller than cutoff, simp is considered as diagonal. The problem is that when simp-->Id, the unitary matrix become degenerate, and extremely sensitive to the small off diag elements.')
 call putel_in_namelist(nm, sum_over_k_dft,'sum_over_k_dft',.false.,'Definition: if true the code is expecting that a summation over k-points was done at the DFT level, this complicates the upfolding of the self energy')
 call putel_in_namelist(nm, keep_both_real_and_matsu_last,'keep_both_real_and_matsu_last',.false.,'Definition: if true keeps both the sigma_output in matsu and real representations')
 call putel_in_namelist(nm, real_smat_coef, 'real_smat_coef', .true., 'Definition: if true the S matrix (projected overlap matrix) will be enforced to be real')
 call putel_in_namelist(nm, force_self_infty_real,'force_self_infty_real',.true.,'Definition: if true it will estimate the self energy at infinite frequency as a real number, and discard the imaginary part')
 call putel_in_namelist(nm,double_counting_zero_self_from_matsu,'double_counting_zero_self_from_matsu',.true.,'Definition: if true the self energy estimated at w=oo for the double counting correction is taken from the matsubara self energy, and used also for the real axis self energy. The point being that the 1/w decay of Sigma(w) prevents a good estimate of Sigma(w=oo)')
 call putel_in_namelist(nm,ed_solver_compute_all_green_functions,'ed_solver_compute_all_green_functions',.true.,'Definition: if true it will tell the ED solver to compute all matrix elements of the Green function, if false only the diagonal will be computed, or if mask_sym_green_ed is present, it will use this mask instead')
 call putel_in_namelist(nm,average_green_ed,'average_green_ed',.false.,'definition: if true will average the green function with the mask_sym_green_ed')
 call putel_in_namelist(nm,tail_linear_scaling,'tail_linear_scaling',-1.d0,'Definition: tail to match the free GF, used for linear scaling where only the low energy part of GF is known')
 call putel_in_namelist(nm,verysilent,'verysilent',.false.,'Definition: if true the solver will be as silent and interact as little as possible with the local directory and outputs')
 call putel_in_namelist(nm,real_sproj_from_onetep,'real_sproj_from_onetep', .true. , 'Definition: condition on the quantities extracted from onetep') 
 call putel_in_namelist(nm,real_eimp_from_tail,   'real_eimp_from_tail',    .true. , 'Definition: condition on the quantities extracted from onetep') 
 call putel_in_namelist(nm,real_eimp_from_onetep, 'real_eimp_from_onetep',  .true. , 'Definition: condition on the quantities extracted from onetep')
 call putel_in_namelist(nm,real_smat_from_onetep, 'real_smat_from_onetep',  .true. , 'Definition:condition on the quantities extracted from onetep')
 call putel_in_namelist(nm,use_simp_from_onetep_for_ortho,'use_simp_from_onetep_for_ortho',.true., 'Definition: uses the Simp matrix obtained by onetep to compute the transformations to non-ortho basis')
 call putel_in_namelist(nm,flag_use_broadening_two_scales_ed,'flag_use_broadening_two_scales_ed',-1,'Definition: if >0 the ED solver will use a broadening factor which assumes two different energy scales')
 call putel_in_namelist(nm,ed_no_real_overide,'ed_no_real_overide',.false.,'Definition: if true, the ED solver will not provide the real axis quantities, this overrides any other combination of options')
      call putel_in_namelist(nm, print_qc, 'print_qc', .false., 'Definition: if true, print out various quantities for the purposes of quality control')

      call look_for_namelist_in_file(nm, './input_onetep_dmft.txt')

      if (fmos_hub1) then
         if (.not. flag_use_slater_in_ed) write (*, *) 'WARNING : FMOS HUBBARD SOLVER : SWITCHING SLATER INTERACTION TO TRUE'
         if (.not. fmos_use) write (*, *) 'WARNING : FMOS HUBBARD SOLVER : SWITCHING FMOS_USE TO TRUE'
         fmos_use = .true.
         flag_use_slater_in_ed = .true.
      endif

      if (nmatsu_long /= 0) then
         if (nmatsu_long < nmatsu_ctqmc) nmatsu_long = nmatsu_ctqmc + 200
      endif

      if (.not. cluster_dmft_green_for_self_consistence) assume_proj_overlap_is_diagonal = .true.

      if (flag_use_jj_basis_for_f) rotation_green_function = .true.

      if (flag_use_slater_in_ed) then
         write (*, *) 'WARNING, using slater interaction for ED solver'
         write (*, *) 'so we impose flag : double_counting_with_no_average=.true.'
         write (*, *) 'the latter flag was indeed introduced for Kanamori interaction'
         double_counting_with_no_average = .true.
      endif

      if (rotation_scheme == 5) then
         INQUIRE (file='mask_user_rot', exist=checkrot)
         if (.not. checkrot) rotation_scheme = 4
      endif

      if (sum_over_k_dft) then
         real_smat_coef = .true.
         real_sproj_from_onetep = .true.
         real_eimp_from_tail = .true.
         real_eimp_from_onetep = .true.
         real_smat_from_onetep = .true.
      endif

!#######################################!
!#######################################!
!#######################################!
   end subroutine

end module

!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

program onestep_dmft_iteration
   !---------------!
   use common_def, only: utils_system_call
   use onetep_variables
   use openmpmod, only: init_openmp
   use genvar, only: iproc
   use timer_mod, only: initialize_timing, finalize_timing
   !---------------!
   implicit none
   logical                :: checkujmat
   complex(kind=DP), allocatable :: green_(:, :, :, :), green_diag(:, :, :), green_temp(:, :), frequ_(:, :), sigma_mat(:, :, :, :)
   real(kind=DP)                :: UU0, Jhund0
   real(kind=DP), allocatable    :: UU_ren0(:, :), JJ_ren0(:, :), UU_ren(:, :), JJ_ren(:, :), double_counting(:, :), edcmatsu(:, :, :)
   real(kind=DP), allocatable    :: eimp(:, :), eimp_nca(:, :)
   real(kind=DP), allocatable    :: occupation_matrix(:, :), occupation_matrixup(:, :), occupation_matrixdn(:, :)
real(kind=DP),allocatable    :: eimp_ed(:,:,:),Zimp_ren_p(:,:,:),Zimp_ren_m(:,:,:),Simp_back(:,:,:),Simp_m1(:,:,:),Simp_rot(:,:,:)
   real(kind=DP), allocatable    :: Simp_root_m(:, :, :), Simp_root_p(:, :, :), Simp_root_gm(:, :, :), Simp_root_gp(:, :, :)
   complex(kind=DP), allocatable :: ddiag(:, :), rotation(:, :, :), diagdens(:, :), mmat(:, :)
   real(kind=DP)                :: denstarget(2)
   integer                :: i, j, k, kk, channels, channelsb, atom, atomb
   integer                :: num, len__, status__, LL, paramagnetic
   character(200)         :: value(100), filename, filename_sigma(2), filename_edc
   complex(kind=DP)             :: frequ
   character(3)           :: label
   character(400)         :: stringdat
   real(kind=DP)                :: bbeta, mmu(2), occupation_orbital
   integer                :: iter_dmft, solver, real_or_imaginary_solver
   character(200)         :: dir_onetep
   character(200)         :: chartemp
   type(string)           :: cc_
   real(kind=DP)                :: densupdn(2), Uaverage
   logical, allocatable    :: list_orb(:)
   complex(kind=DP), allocatable :: UC(:, :, :, :, :)
   real(kind=DP)                :: gck(0:3, -3:3, -3:3, 0:3)
   real(kind=DP)                :: Fn(0:3, 0:3)
   complex(kind=DP), allocatable :: T2C(:, :), T2Cp(:, :)
   integer, allocatable    :: UCC(:, :, :, :)
   complex(kind=DP), allocatable :: UCCr(:, :, :)

   call initialize_timing()
   
   call init_my_input_variables
   write (*, *) 'init data'
   call init_run

   call init_openmp(silent = (iproc /= 1))

   if (iter_dmft == niter_dmft) ed_real_frequ = ed_real_frequ_last

   write (*, *) 'collect greens function and Sigma'
   call collect_greens(1)
   call collect_sigmas(1)
   if (paramagnetic /= 1) then
      call collect_greens(2)
      call collect_sigmas(2)
   endif

   call rotations_and_interactions_and_double_counting

   write (*, *) 'compute hybridizations'
   call get_hybridizations(1)
   if (paramagnetic /= 1) then
      call get_hybridizations(2)
      if (solver == 2 .or. solver == 3) then
         call utils_system_call("merge_ac Ac.inp "//TRIM(ADJUSTL(toString(channels))))
      endif
      call utils_system_call("merge_delta delta_input "//TRIM(ADJUSTL(toString(channels))))
   endif

   write (*, *) 'start DMFT run'
   call dmft_run

   if (.not. verysilent) call utils_system_call(" echo 'sigma files after dmft run' `ls sigma_output* 2> /dev/null` ")

   write (*, *) 'plug back double counting in self energy'

   if (solver /= 4 .or. iter_dmft /= niter_dmft .or. .not. last_iter_is_real) then
      call plug_back_double_counting_and_fermi_e_in_self
   else
      if (keep_both_real_and_matsu_last) then
         call plug_back_double_counting_and_fermi_e_in_self(keepboth=.true.)
         call plug_back_double_counting_and_fermi_e_in_self_real(keepboth=.true.)
      else
         call plug_back_double_counting_and_fermi_e_in_self_real
      endif
   endif

   write (*, *) 'now cleanup'

   if (.not. verysilent) call utils_system_call(" echo 'sigma files after double counting correction' `ls sigma_output* 2> /dev/null` ")

   call cleanup

   if (paramagnetic == 1) then
      write (*, *) 'paramagnetic case'
      write (*, *) 'copying sigma file'
      cc_ = trim(adjustl(filename_sigma(1)))
      call replace_in_string(cc_, '_1', '_2', 'last')
      chartemp = cc_
      write (*, *) 'copy : ', TRIM(ADJUSTL(filename_sigma(1)))
      write (*, *) 'to   : ', TRIM(ADJUSTL(chartemp))
      call utils_system_call("cp "//TRIM(ADJUSTL(filename_sigma(1)))//" "//TRIM(ADJUSTL(chartemp)))
      write (*, *) 'done'
   endif

   if (.not. verysilent) call utils_system_call(" echo 'sigma files after paramagnetic copy' `ls sigma_output* 2> /dev/null` ")

   call finalize_timing()

contains

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine plug_back_double_counting_and_fermi_e_in_self_real(keepboth)

      use common_def, only: utils_assert, utils_abort

      implicit none
      complex(kind=DP)       ::  mat(channels, channels), sigout(2, ed_real_frequ, channels, channels)
      integer          ::  k, kk_, tt, j_
      logical          ::  test_
      logical, optional ::  keepboth

      call utils_assert(solver == 4, 'Plugging back DC for real self energy is only valid for ed solver')

      if (paramagnetic /= 1) then
         j_ = 2
      else
         j_ = 1
      endif

      !===================================================================================!
      do kk_ = 1, j_

         if (kk_ == 1) INQUIRE (file='_sigma_output_full_real_1', EXIST=test_)
         if (kk_ == 2) INQUIRE (file='_sigma_output_full_real_2', EXIST=test_)
         if (kk_ == 1 .and. test_) write (*, *) 'FILE _sigma_output_full_real_1 exists'
         if (kk_ == 2 .and. test_) write (*, *) 'FILE _sigma_output_full_real_2 exists'
         if (kk_ == 1 .and. test_) open (unit=10004, file='_sigma_output_full_real_1', form='unformatted')
         if (kk_ == 2 .and. test_) open (unit=10004, file='_sigma_output_full_real_2', form='unformatted')
         call utils_assert(test_, 'file _sigma_output_full_real_(1,2) is missing')
         !----------------------------------------------------------------------------------!
         do i = 1, ed_real_frequ
            read (10004) mat
            call utils_assert(size(sigout, 2) >= i, 'Mismatched arrays when dumping back double counting to Sigma')
            sigout(kk_, i, :, :) = mat
            if (force_self_infty_real .and. i == ed_real_frequ - 1) sigout(kk_, i, :, :) = real(sigout(kk_, i, :, :))
         enddo
         !----------------------------------------------------------------------------------!
         close (10004); 
         write (*, *) 'number of frequencies, frequmin, frequmax : ', ed_real_frequ, ed_frequ_min, ed_frequ_max

         if (kk_ == 1) call utils_system_call("rm _sigma_output_full_real_1")
         if (kk_ == 2) call utils_system_call("rm _sigma_output_full_real_2")
         if (.not. present(keepboth)) then
            if (kk_ == 1) call utils_system_call("rm _sigma_output_full_1")
            if (kk_ == 2) call utils_system_call("rm _sigma_output_full_2")
         endif
      enddo
      !===================================================================================!

      if (rotate_to_ortho_basis) then
         if (.not. flag_use_jj_basis_for_f) then
            do kk_ = 1, j_
               do i = 1, ed_real_frequ
                  call rotate_back_sigma_ortho(sigout(kk_, i, :, :), kk_)
               enddo
            enddo
         else
            call utils_abort('Using F with rotation ortho is invalid')
         endif
      endif

      do kk_ = 1, j_
         if (double_counting_zero_self_av) then
            ccctemp_av = sum(diag(sigout(kk_, ed_real_frequ, :, :)))/dble(size(sigout, 3))
            if (double_counting_zero_self_from_matsu) &
          & ccctemp_av = sum(diag(edcmatsu(kk_, :, :)))/dble(size(sigout, 3))

         endif
         do tt = 1, channels
            if (double_counting_zero_self) then
               ccctemp = real(sigout(kk_, ed_real_frequ, tt, tt))
               if (double_counting_zero_self_from_matsu) ccctemp = edcmatsu(kk_, tt, tt)
               if (double_counting_zero_self_av) ccctemp = real(ccctemp_av)
            else
               ccctemp = double_counting(kk_, tt)
            endif
            write (*, *) 'REAL FREQU CHANNEL - CORRECTION : ', tt, ccctemp
            do i = 1, ed_real_frequ
               sigout(kk_, i, tt, tt) = sigout(kk_, i, tt, tt) - ccctemp
            enddo
            sigout(kk_, ed_real_frequ, tt, tt) = ccctemp
         enddo
      enddo

      if (rotation_green_function) then
         if (.not. flag_use_jj_basis_for_f) then
            do kk_ = 1, j_
               do i = 1, ed_real_frequ
                  call rotate_back_sigma(sigout(kk_, i, :, :), kk_)
               enddo
            enddo
         else
            ! back from J-J to cubic harmonics !
            do i = 1, ed_real_frequ
               call rotate_sigma_green_to_f_back(sigout(1, i, :, :), sigout(2, i, :, :))
            enddo
         endif
      endif

      if (sum_over_k_dft) Then
         do kk_ = 1, j_
            sigout(kk_, ed_real_frequ - 6, :, :) = Simp_m1(kk_, :, :)
         enddo
      endif

      do kk_ = 1, j_
         open (unit=10003, file='sigma_real_scratch', form='unformatted')
         do i = 1, ed_real_frequ
            write (10003) sigout(kk_, i, :, :)
         enddo
         close (10003)
         if (.not. present(keepboth)) then
            call utils_system_call("mv sigma_real_scratch "//trim(adjustl(filename_sigma(kk_))))
         else
            call utils_system_call("mv sigma_real_scratch real_"//trim(adjustl(filename_sigma(kk_))))
         endif
      enddo

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine plug_back_double_counting_and_fermi_e_in_self(keepboth)

      use common_def, only: utils_assert, utils_abort

      implicit none
      real(kind=DP)          ::  temp(1 + 2*channels)
      complex(kind=DP)       ::  mat(channels, channels), sigout(2, n_frequ, channels, channels)
      integer          ::  k, kk_, tt, ttt, j_
      logical          ::  full_, test_1, test_2
      logical, optional ::  keepboth

      if (paramagnetic /= 1) then
         j_ = 2
      else
         j_ = 1
      endif

      write (*, *) 'writing sigma in file'

      !==================================================================================!

      do kk_ = 1, j_

         write (*, *) 'subtracting double counting from Sigma : ', double_counting(kk_, :)

         if ((solver == 1 .or. solver == 2 .or. solver == 3) .and. nmatsu_long > 0) then
            call copy_sigma_position_nlong_to_position_nw_minus_1(filename_sigma(kk_))
         endif

         open (unit=10002, file=trim(adjustl(filename_sigma(kk_))))

         full_ = .false.
         if (cluster_dmft_green_for_self_consistence) then
            full_ = .true.
            if (kk_ == 1) INQUIRE (file='_sigma_output_full_1', EXIST=test_1)
            if (kk_ == 2) INQUIRE (file='_sigma_output_full_2', EXIST=test_2)
            if (kk_ == 1) full_ = test_1
            if (kk_ == 2) full_ = test_2
            if (kk_ == 1 .and. full_) write (*, *) 'FILE _sigma_output_full_1 exists'
            if (kk_ == 2 .and. full_) write (*, *) 'FILE _sigma_output_full_2 exists'
         endif

         if (full_ .and. kk_ == 1) open (unit=10004, file='_sigma_output_full_1', form='unformatted')
         if (full_ .and. kk_ == 2) open (unit=10004, file='_sigma_output_full_2', form='unformatted')

         if (solver == 1) then
            read (10002, *)
         elseif (solver == 2 .or. solver == 3 .or. solver == 4 .or. solver == 5) then
         else
            call utils_abort("Solver not yet implemented")
         endif

         write (*, *) 'writing down new SELF ENERGY Sigma, mixing is : ', mixing
         !----------------------------------------------------------------------------------!
         do i = 1, n_frequ

            read (10002, *) (temp(k), k=1, 1 + 2*channels)
            mat = 0.d0
            if (full_) read (10004) mat
            if (i == n_frequ - 1) then
               write (*, '(a,200f12.7)') 'uncorrected self energy is at oo RE     : ', (temp(k), k=2, size(temp), 2)
               write (*, '(a,200f12.7)') 'uncorrected self energy is at oo IM     : ', (temp(k), k=3, size(temp), 2)
               write (*, *) 'double counting is                      : ', double_counting(kk_, :)
               write (*, *) 'here we send Sigma-Edc back to onetep   '
               write (*, *) 'This is done to get the double counting correction in onetep'
            endif
            if (i == n_frequ) then
               write (*, '(a,200f12.3)') 'uncorrected self energy is at n_frequ : ', (temp(k), k=2, size(temp), 2)
               write (*, *) 'double counting is                    : ', double_counting(kk_, :)
            endif
            if (.not. full_) then
               do tt = 1, channels
                  mat(tt, tt) = temp(1 + 2*tt - 1) + imi*temp(1 + 2*tt)
               enddo
            endif

            if (i <= n_frequ - nspec_frequ .or. i == n_frequ - 1 .or. i == n_frequ - 5) then
               mat = mixing*mat + (1.d0 - mixing)*sigma_mat(kk_, i, :, :)  ! mixing with former iteration
            endif

        if(i==n_frequ-nspec_frequ) write(*,*) 'CHECK MIXING : ', maxval(abs(real(mat-sigma_mat(kk_,i,:,:)))),maxval(abs(aimag(mat-sigma_mat(kk_,i,:,:))))

            if (full_ .and. flag_blank_out_sigma_offdiag_for_testing) then
               do tt = 1, channels
                  do ttt = 1, channels
                     if (tt /= ttt) then
                        mat(tt, ttt) = 0.d0
                     endif
                  enddo
               enddo
            endif

            call utils_assert(size(sigout, 2) >= i, 'Array size mismatch when dmping back double counting to Sigma')

            sigout(kk_, i, :, :) = mat

            if (force_self_infty_real .and. i == n_frequ - 1) sigout(kk_, i, :, :) = real(sigout(kk_, i, :, :))

         enddo
         !----------------------------------------------------------------------------------!

         close (10002)

         if (full_) then
            close (10004)
            if (kk_ == 1) call utils_system_call("rm _sigma_output_full_1")
            if (kk_ == 2) call utils_system_call("rm _sigma_output_full_2")
            if (.not. present(keepboth)) then
               if (kk_ == 1) call utils_system_call("rm _sigma_output_full_real_1")
               if (kk_ == 2) call utils_system_call("rm _sigma_output_full_real_2")
            endif
         endif

      enddo
      !==================================================================================!

      if (rotate_to_ortho_basis) then
         if (.not. flag_use_jj_basis_for_f) then
            do kk_ = 1, j_
               do i = 1, n_frequ
                  call rotate_back_sigma_ortho(sigout(kk_, i, :, :), kk_)
               enddo
            enddo
         else
            call utils_abort("Using f electron and ortho rotation")
         endif
      endif

      do kk_ = 1, j_
         edcmatsu(kk_, :, :) = sigout(kk_, n_frequ - 1, :, :)
         if (double_counting_zero_self_av) then
            ccctemp_av = sum(diag(sigout(kk_, n_frequ - 1, :, :)))/dble(size(sigout, 3))
         endif
         do tt = 1, channels
            if (double_counting_zero_self) then
               ccctemp = real(sigout(kk_, n_frequ - 1, tt, tt))
               if (double_counting_zero_self_av) ccctemp = real(ccctemp_av)
            else
               ccctemp = double_counting(kk_, tt)
            endif
            write (*, *) 'CHANNEL - CORRECTION : ', tt, ccctemp
            do i = 1, n_frequ
               sigout(kk_, i, tt, tt) = sigout(kk_, i, tt, tt) - ccctemp
            enddo
            sigout(kk_, n_frequ, tt, tt) = ccctemp
         enddo
      enddo

      if (rotation_green_function) then
         if (.not. flag_use_jj_basis_for_f) then
            do kk_ = 1, j_
               do i = 1, n_frequ
                  call rotate_back_sigma(sigout(kk_, i, :, :), kk_)
               enddo
            enddo
         else
            ! back from J-J to cubic harmonics !
            do i = 1, n_frequ
               call rotate_sigma_green_to_f_back(sigout(1, i, :, :), sigout(2, i, :, :))
            enddo
         endif
      endif

      if (sum_over_k_dft) Then
         do kk_ = 1, j_
            sigout(kk_, n_frequ - 6, :, :) = Simp_m1(kk_, :, :)
         enddo
      endif

      do kk_ = 1, j_
         open (unit=10003, file='sigma_scratch', form='unformatted')
         do i = 1, n_frequ
            write (10003) sigout(kk_, i, :, :)
         enddo
         close (10003)
         call utils_system_call("mv sigma_scratch "//trim(adjustl(filename_sigma(kk_))))
      enddo

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!

   subroutine cleanup
      deallocate (frequ_, green_, green_diag, sigma_mat, eimp, eimp_ed, eimp_nca, Simp_back, Simp_m1, Simp_rot)
      deallocate (Zimp_ren_p, Zimp_ren_m, Simp_root_m, Simp_root_p, UU_ren, JJ_ren, UU_ren0, JJ_ren0, edcmatsu)
      deallocate (Simp_root_gm, Simp_root_gp)
   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!

   subroutine dmft_run

      use common_def, only: utils_assert, utils_abort, utils_system_call
      use timer_mod,  only: start_timer, stop_timer

      implicit none
      complex(kind=DP)   ::  mat_tmp0(2*LL + 1, 2*LL + 1), mat_tmp20(2*(2*LL + 1), 2*(2*LL + 1))
      complex(kind=DP)   ::   mat_tmp(2*LL + 1, 2*LL + 1), mat_tmp2(2*(2*LL + 1), 2*(2*LL + 1))
      real(kind=DP)      ::  ttmp(4*(2*LL + 1))
      real(kind=DP)      ::  U_(channels, channels), rot_trans_pm(2*LL + 1, 2*LL + 1), rot_trans_af((2*LL + 1)*2, (2*LL + 1)*2)
      integer      ::  i_, i, ii, jj, kk, nheaders, ngwfs_indexes(2*LL + 1)
      logical      ::  check, check_ngwfs, extend

      inquire (file='mask_ngwfs', exist=check_ngwfs)
      if (check_ngwfs) then
         open (unit=9812, file='mask_ngwfs')
         read (9812, *) (ngwfs_indexes(ii), ii=1, 2*LL + 1)
         close (9812)
      else
         ngwfs_indexes = (/(ii, ii=1, 2*LL + 1)/)
      endif

      write (*, *) 'NGWFS indices : ', ngwfs_indexes

!###################COPYING FILES FOR ATOMIC SOLVER##################!
!###################COPYING FILES FOR ATOMIC SOLVER##################!
!###################COPYING FILES FOR ATOMIC SOLVER##################!

      extend = .not. ((channels == 3 .and. LL == 1) .or. (channels == 5 .and. LL == 2) .or. (channels == 7 .and. LL == 3))

      if (LL == 0 .or. LL > 3) goto 1654

      ! Defines rotated/transformed/renormalized cubic harmonics in terms of spherical harmonics
      ! This transformation is used to compute the Slater interaction vertex

      write (*, *) 'Cubic Harmonics in terms of spherical harmonics'

      !===============================================================================================!
      !===============================================================================================!
      !===============================================================================================!

      if (paramagnetic == 1) then

         if (allocated(T2C)) deallocate (T2C, T2Cp)
         allocate (T2C(2*LL + 1, channels), T2Cp(channels, 2*LL + 1))
         nheaders = 3 + 2*LL + 1
         if (cubic == 1) then
        call utils_system_call("cp "//TRIM(ADJUSTL(dir_onetep))//"/utils/INPUTS/Trans_cubicL"//TRIM(ADJUSTL(toString(LL)))//".dat ./Trans.dat")
         else
    call utils_system_call("cp "//TRIM(ADJUSTL(dir_onetep))//"/utils/INPUTS/Trans_sphericalL"//TRIM(ADJUSTL(toString(LL)))//".dat ./Trans.dat")
         endif

         inquire (file='Trans.loc.pm.dat', exist=check)
         if (check) then
            open (unit=1745, file='Trans.loc.pm.dat', form='unformatted')
            read (1745)
            read (1745) rot_trans_pm
            close (1745)
         else
            rot_trans_pm = Id(size(rot_trans_pm, 1))
         endif

         if (rotation_green_function) rot_trans_pm = matmul_ex(rot_trans_pm, real(rotation(:, :, 1), kind=8), extend)
         if (rotate_to_ortho_basis) rot_trans_pm = matmul_ex(rot_trans_pm, real(Simp_rot(1, :, :), kind=8), extend)

         open (unit=1745, file='Trans.dat')
         do i = 1, nheaders; read (1745, *); enddo
         write (*, *) 'READING TRANS.DAT'
         do i = 1, 2*LL + 1
            write (*, *) 'reading Trans.dat, line : ', i, 2*LL + 1
            read (1745, *) (ttmp(j), j=1, 2*(2*LL + 1))
            do j = 1, (2*LL + 1)
               mat_tmp0(i, j) = ttmp(2*j - 1) + imi*ttmp(2*j)
            enddo
            write (*, '(200f8.4)') (mat_tmp0(i, j), j=1, 2*LL + 1)
         enddo

         do i = 1, 2*LL + 1
            i_ = ngwfs_indexes(i)
            mat_tmp(i, :) = mat_tmp0(i_, :)
            write (*, *) 'state i (re-ordered) : ', mat_tmp(i, :)
         enddo

         if (flag_get_t2c_real) call change_phase_T2C(mat_tmp(1:channels, 1:2*LL + 1))

         do i = 1, 2*LL + 1
            write (*, *) 'channel i, norm of state : ', sum(abs(mat_tmp(i, :))**2)
         enddo
         if (maxval(abs(matmul(rot_trans_pm, transpose(rot_trans_pm)) - Id(size(rot_trans_pm, 1)))) > 1.d-5) then
            call utils_abort('Rotation matrix not orthogonal')
         endif

         if (rotate_int_after_earlier_transfo) mat_tmp = matmul(transpose(rot_trans_pm), mat_tmp)

         close (1745)
         call utils_system_call(" head -n "//trim(adjustl(toString(nheaders)))//" Trans.dat > tmp_trans ")
         call utils_system_call(" mv tmp_trans Trans.dat ")

         open (unit=1745, file='Trans.dat', ACCESS='append')
         i_ = 0
         do i = 1, 2*LL + 1
            if (list_orb(i)) then
               i_ = i_ + 1
               T2C(:, i_) = mat_tmp(i, :)
               T2Cp(i_, :) = conjg(mat_tmp(i, :))
            endif
         enddo

         do i = 1, channels
            if (.not. rotate_ortho_av_renorm_int .and. rotate_to_ortho_basis) then
               T2C(:, i) = T2C(:, i)*Zimp_ren_m(1, i, i)
               T2Cp(i, :) = conjg(T2C(:, i))
            endif
            write (1745, '(200f15.8)') (T2C(j, i), j=1, 2*LL + 1)
         enddo
         close (1745)

         !===============================================================================================!
         !===============================================================================================!
         !===============================================================================================!
      else
         !===============================================================================================!
         !===============================================================================================!
         !===============================================================================================!

         if (allocated(T2C)) deallocate (T2C, T2Cp)
         allocate (T2C(2*(2*LL + 1), 2*channels), T2Cp(2*channels, 2*(2*LL + 1)))

         nheaders = 3 + 2*(2*LL + 1)
         if (cubic == 1) then
            call utils_system_call("cp " // TRIM(ADJUSTL(dir_onetep)) // "/utils/INPUTS/Trans_cubic_spinorbL" &
                  // TRIM(ADJUSTL(toString(LL)))//".dat ./Trans.dat")
         else
            call utils_system_call("cp " // TRIM(ADJUSTL(dir_onetep)) // "/utils/INPUTS/Trans_spherical_spinorbL" &
                  // TRIM(ADJUSTL(toString(LL)))//".dat ./Trans.dat")
         endif
         inquire (file='Trans.loc.spin.dat', exist=check)

         if (check) then
            open (unit=1745, file='Trans.loc.spin.dat', form='unformatted')
            read (1745)
            read (1745) rot_trans_af
            close (1745)
         else
            rot_trans_af = Id(size(rot_trans_af, 1))
         endif
         if (rotation_green_function) then
     rot_trans_af(         1:  (2*LL+1),         1:  (2*LL+1))=matmul_ex(rot_trans_af(         1:  (2*LL+1),         1:  (2*LL+1)),real(rotation(:,:,1),kind=8),extend)
     rot_trans_af((2*LL+1)+1:2*(2*LL+1),(2*LL+1)+1:2*(2*LL+1))=matmul_ex(rot_trans_af((2*LL+1)+1:2*(2*LL+1),(2*LL+1)+1:2*(2*LL+1)),real(rotation(:,:,2),kind=8),extend)
         endif
         if (rotate_to_ortho_basis) then
     rot_trans_af(         1:  (2*LL+1),         1:  (2*LL+1))=matmul_ex(rot_trans_af(         1:  (2*LL+1),         1:  (2*LL+1)),real(Simp_rot(1,:,:),kind=8),extend)
     rot_trans_af((2*LL+1)+1:2*(2*LL+1),(2*LL+1)+1:2*(2*LL+1))=matmul_ex(rot_trans_af((2*LL+1)+1:2*(2*LL+1),(2*LL+1)+1:2*(2*LL+1)),real(Simp_rot(2,:,:),kind=8),extend)
         endif

         open (unit=1745, file='Trans.dat')
         do i = 1, nheaders; read (1745, *); enddo
         write (*, *) 'READING Trans.dat'

         do i = 1, 2*(2*LL + 1)
            read (1745, *) (ttmp(j), j=1, 4*(2*LL + 1))
            do j = 1, 2*(2*LL + 1)
               mat_tmp20(i, j) = ttmp(2*j - 1) + imi*ttmp(2*j)
            enddo
            write (*, '(200f8.3)') (mat_tmp20(i, j), j=1, 2*(2*LL + 1))
         enddo

         do i = 1, 2*(2*LL + 1)
            if (i <= (2*LL + 1)) then
               i_ = ngwfs_indexes(i)
            else
               i_ = ngwfs_indexes(i - (2*LL + 1)) + (2*LL + 1)
            endif
            mat_tmp2(i, :) = mat_tmp20(i_, :)
         enddo

         if (flag_get_t2c_real) call change_phase_T2C(mat_tmp2(1:channels, 1:2*LL + 1))
         if (flag_get_t2c_real) call change_phase_T2C(mat_tmp2(channels + 1:2*channels, (2*LL + 1) + 1:2*(2*LL + 1)))

         if (rotate_int_after_earlier_transfo) mat_tmp2 = matmul(transpose(rot_trans_af), mat_tmp2)

         close (1745)

         call utils_system_call(" head -n "//trim(adjustl(toString(nheaders)))//" Trans.dat > tmp_trans ")
         call utils_system_call(" mv tmp_trans Trans.dat ")

         open (unit=1745, file='Trans.dat', ACCESS='append')

         i_ = 0
         do i = 1, 2*(2*LL + 1)
            if (i <= (2*LL + 1)) then
               if (list_orb(i)) then
                  i_ = i_ + 1
                  T2C(:, i_) = mat_tmp2(i, :)
                  T2Cp(i_, :) = conjg(mat_tmp2(i, :))
               endif
            else
               if (list_orb(i - (2*LL + 1))) then
                  i_ = i_ + 1
                  T2C(:, i_) = mat_tmp2(i, :)
                  T2Cp(i_, :) = conjg(mat_tmp2(i, :))
               endif
            endif
         enddo

         do i = 1, 2*channels
            if (.not. rotate_ortho_av_renorm_int .and. rotate_to_ortho_basis) then
               if (i <= channels) T2C(:, i) = T2C(:, i)*Zimp_ren_m(1, i, i)
               if (i > channels) T2C(:, i) = T2C(:, i)*Zimp_ren_m(2, i - channels, i - channels)
               T2Cp(i, :) = conjg(T2C(:, i))
            endif
            write (1745, '(200f15.8)') (T2C(j, i), j=1, 2*(2*LL + 1))
         enddo

         close (1745)

      endif
      !===============================================================================================!
      !===============================================================================================!
      !===============================================================================================!

      if (extend) then
         if (solver == 1 .or. solver == 2 .or. solver == 3) then
            write (*, *) 'ERROR L=1,2,3 only were implemented for CTQMC/NCA/OCA, please update the code, critical'
            write (*, *) 'Note that projections were not yet implemented for these solvers'
            call utils_abort('Projectors not yet implemented')
         endif
      endif

      ! now come the calculation of the vertex

      if ((solver == 4 .and. flag_use_slater_in_ed) .or. &
           & (solver < 4 .and. (  &
                               & (.not. rotate_ortho_av_renorm_int .and. rotate_to_ortho_basis) &
                         & .or. checkujmat &
                         &  )  &
           &  )) then
         call generate_UC_dat_files(UU, Jhund) !--> this goes to ED-Slater or CTQMC-renormalized-interaction
         UU_ren = 0.d0
         JJ_ren = 0.d0
         UU = 0.d0
         if (solver == 4) Jhund = 0.d0  !Jhund is needed for atom_d.py script for CTQMC/OCA/NCA solvers
      endif

1654  continue

!###################PREPARING ATOMIC SOLVER##########################!
!###################PREPARING ATOMIC SOLVER##########################!
!###################PREPARING ATOMIC SOLVER##########################!

!-----------------------------------------------------------------!
!-----------------------------------------------------------------!
!-----------------------------------------------------------------!
      if ((solver == 1 .or. solver == 2 .or. solver == 3) .and. .not. read_cix_file) then
         write (*, *) '##################### ATOMIC SOLVER #########################'
         write (*, *) 'atomic number J U and SO coupling : ', Jhund, UU, spin_orbit

         if (LL == 1 .or. LL == 2) then
            stringdat = trim(adjustl(atom_d_command))//' Eimp=['
         elseif (LL == 3) then
            if (solver /= 1) then
               call utils_abort('L = 3 not implemented for other solver than CTQMC')
            endif
            stringdat = 'atomh -Eimp ['
         endif

         do i = 1, channels
            if (solver /= 2 .and. solver /= 3) then
               stringdat = TRIM(ADJUSTL(stringdat))//TRIM(ADJUSTL(toString(eimp(1, i))))
            else
               stringdat = TRIM(ADJUSTL(stringdat))//TRIM(ADJUSTL(toString(eimp_nca(1, i))))
            endif
            if (i < channels) stringdat = TRIM(ADJUSTL(stringdat))//","
         enddo
         if (paramagnetic /= 1) then
            stringdat = TRIM(ADJUSTL(stringdat))//','
            do i = 1, channels
               if (solver /= 2 .and. solver /= 3) then
                  stringdat = TRIM(ADJUSTL(stringdat))//TRIM(ADJUSTL(toString(eimp(2, i))))
               else
                  stringdat = TRIM(ADJUSTL(stringdat))//TRIM(ADJUSTL(toString(eimp_nca(2, i))))
               endif
               if (i < channels) stringdat = TRIM(ADJUSTL(stringdat))//","
            enddo
         endif

         if (LL == 1 .or. LL == 2) then
        stringdat=TRIM(ADJUSTL(stringdat))//'] qOCA=1 Nmax=2625 J='//trim(adjustl(toString(Jhund)))//' Eoca=1.0 l='//trim(adjustl(toString(LL)))//' n=[6,7,8,9,10] cx='//trim(ADJUSTL(toString(spin_orbit)))//' mOCA=0.001 Ex=[0.5,0.5,0.5] Ncentral=[6,7,8,9,10] qatom=0 Ep=[3.0,3.0,3.0]' 
            if (solver == 2 .or. solver == 3) then
               stringdat = TRIM(ADJUSTL(stringdat))//' OCA_G=True '
            else
               stringdat = TRIM(ADJUSTL(stringdat))//' OCA_G=False '
            endif
         elseif (LL == 3) then
        stringdat=TRIM(ADJUSTL(stringdat))//'] -Nmax 2625 -ImpTot -l 3 -J '//trim(adjustl(toString(Jhund)))//' -cx '//trim(ADJUSTL(toString(spin_orbit)))
            if (paramagnetic == 1) then
               stringdat = TRIM(ADJUSTL(stringdat))//' -Impq [0,1,2,3,4,5,6,0,1,2,3,4,5,6]     '
            else
               stringdat = TRIM(ADJUSTL(stringdat))//' -Impq [0,1,2,3,4,5,6,7,8,9,10,11,12,13] '
            endif
        stringdat=TRIM(ADJUSTL(stringdat))//' -ns '//trim(adjustl(toString(lowest_occupancy)))//' -ne '//trim(adjustl(toString(highest_occupancy)))
         endif

         stringdat = TRIM(ADJUSTL(stringdat))//' > gen_atom.log '

         write (*, *) TRIM(ADJUSTL(stringdat))

         if ((LL == 1 .or. LL == 2) .and. use_custom_command_for_atomd) then
            call utils_system_call(" rm atomic_diag_script.out")
            call utils_system_call(" echo '#!/bin/bash' >> atomic_diag_script.out")
            call utils_system_call(" echo 'j=`pwd`' >> atomic_diag_script.out")
            call utils_system_call(" echo 'k=`basename $j`' >> atomic_diag_script.out")
            call utils_system_call(" ls ~/tmp || mkdir ~/tmp")
            call utils_system_call(" echo 'cd ~/tmp'  >> atomic_diag_script.out")
            call utils_system_call(" echo 'mkdir tmp_atomd$k' >> atomic_diag_script.out")
            call utils_system_call(" echo 'cd tmp_atomd$k' >> atomic_diag_script.out")
            call utils_system_call(" echo 'cp $j/Trans* ./ ' >> atomic_diag_script.out")
            call utils_system_call(" echo 'cp $j/Uc.* ./ ' >> atomic_diag_script.out")
            call utils_system_call(" echo ' "//TRIM(ADJUSTL(stringdat))//" ' >> atomic_diag_script.out")
            call utils_system_call(" echo 'cd ..' >> atomic_diag_script.out")
            call utils_system_call(" echo 'cp -r tmp_atomd$k/* $j' >> atomic_diag_script.out")
            call utils_system_call(" echo 'rm -r tmp_atomd$k' >> atomic_diag_script.out")
            call utils_system_call(" echo 'cd $j' >> atomic_diag_script.out")
            call utils_system_call(" chmod +x atomic_diag_script.out")
            call utils_system_call(" ./atomic_diag_script.out")
            call utils_system_call(" rm ./atomic_diag_script.out")
         else
            call utils_system_call(TRIM(ADJUSTL(stringdat)))
         endif

      endif

!######################CALLING ATOMIC SOLVER##########################!
!######################CALLING ATOMIC SOLVER##########################!
!######################CALLING ATOMIC SOLVER##########################!
!######################CALLING ATOMIC SOLVER##########################!

      call utils_system_call("mv LOCAL_DOS* FULL_DOS* "//TRIM(ADJUSTL(filename))//"  > /dev/null 2>&1 ")

      if (solver == 1) then
         call utils_system_call("cp "//TRIM(ADJUSTL(dir_onetep))//"/utils/INPUTS/PARAMS.ctqmc ./PARAMS")
      elseif (solver == 4 .or. solver == 5) then
         write (*, *) 'DO NOT COPY ED PARAM, GENERATE THEM ON THE FLY'
      elseif (solver == 2 .or. solver == 3) then
         call utils_system_call("cp "//TRIM(ADJUSTL(dir_onetep))//"/utils/INPUTS/PARAMS.oca ./PARAMS")
      else
         call utils_abort('solver not implemented')
      endif

      call utils_system_call(" echo 'sweep_exact='"//trim(adjustl(toString(1000)))//"     > hf.in ")
      call utils_system_call(" echo 'spline_numb=10000'                                  >> hf.in ")
      call utils_system_call(" echo 'fourier_smoothing=0.00001'                          >> hf.in ")
      if (hf_hartree) then
         call utils_system_call(" echo 'impose_hartree=.true.'                              >> hf.in ")
      else
         call utils_system_call(" echo 'impose_hartree=.false.'                             >> hf.in ")
      endif
      if (hf_average) then
         call utils_system_call(" echo 'average_over_time_slices=.true.'                    >> hf.in ")
      else
         call utils_system_call(" echo 'average_over_time_slices=.false.'                   >> hf.in ")
      endif
      if (hf_krauth) then
         call utils_system_call(" echo 'use_krauth_splining=.true.'                         >> hf.in ")
      else
         call utils_system_call(" echo 'use_krauth_splining=.false.'                        >> hf.in ")
      endif
      call utils_system_call(" echo 'shift_tau_of_boundary=0.0'                          >> hf.in ")
      call utils_system_call(" echo 'U_coul_param=1.'                                    >> hf.in ")
      call utils_system_call(" echo 'enforce_boundary=.true.'                            >> hf.in ")
      call utils_system_call(" echo 'substract_hf=.true.'                                >> hf.in ")
      call utils_system_call(" echo 'substract_hf_im=.true.'                             >> hf.in ")
      call utils_system_call(" echo 'pole_self=0.1'                                      >> hf.in ")
      call utils_system_call(" echo 'impose_causality=.true.'                            >> hf.in ")
      call utils_system_call(" echo 'nsweep='"//trim(adjustl(toString(mcs_ctqmc)))//"    >> hf.in ")
      call utils_system_call(" echo 'ntau='"//trim(adjustl(toString(ntauhf)))//"         >> hf.in ")
      call utils_system_call(" echo 'tail_energy='"//trim(adjustl(toString(3.5d0)))//"   >> hf.in ")
      call utils_system_call(" echo 'warmup='"//trim(adjustl(toString(mcs_ctqmc/20)))//" >> hf.in ")
      call utils_system_call(" echo 'measure_cycle='"//trim(adjustl(toString(100)))//"    >> hf.in ")
      if (paramagnetic_ed == 1) then
         call utils_system_call("echo paramagnetic=.true. >> hf.in")
      else
         call utils_system_call("echo paramagnetic=.false. >> hf.in")
      endif

      call utils_system_call("cp "//TRIM(ADJUSTL(dir_onetep))//"/utils/INPUTS/ed.* .")

      if (average_green_ed) then
         call utils_system_call("echo 1 > ed.average ")
      endif

      if (flag_use_broadening_two_scales_ed <= 0) then
         call utils_system_call("echo flag_idelta_two_scales_ed=0 >> ed.in ")
      else
         call utils_system_call("echo flag_idelta_two_scales_ed="//trim(adjustl(toString(flag_use_broadening_two_scales_ed)))//"  >> ed.in ")
      endif
      call utils_system_call("echo dist_max="//trim(adjustl(toString(dist_max)))//" >> ed.in")
      call utils_system_call("echo nsec="//trim(adjustl(toString(nsec)))//" >> ed.in")
      call utils_system_call("echo nsec0="//trim(adjustl(toString(nsec0)))//" >> ed.in")
      call utils_system_call("echo which_lanczos="//trim(adjustl(which_lanczos))//" >> ed.in")
      call utils_system_call("echo Neigen="//trim(adjustl(toString(Neigen)))//" >> ed.in")
      call utils_system_call("echo Block_size="//trim(adjustl(toString(Block_size)))//" >> ed.in")
      call utils_system_call("echo ncpt_approx="//trim(adjustl(toString(ncpt_approx)))//" >> ed.in")
      call utils_system_call("echo cpt_upper_bound="//trim(adjustl(toString(cpt_upper_bound)))//" >> ed.in")
      call utils_system_call("echo cpt_lagrange="//trim(adjustl(toString(cpt_lagrange)))//" >> ed.in")
      call utils_system_call("echo iwindow="//trim(adjustl(toString(iwindow)))//" >> ed.in")
      call utils_system_call("echo fmos_iter="//trim(adjustl(toString(fmos_iter)))//" >> ed.in")
      call utils_system_call("echo  fmos_mix="//trim(adjustl(toString(fmos_mix)))//" >> ed.in")
      if (fmos_use) then
         call utils_system_call("echo fmos=.true.       >> ed.in")
      else
         call utils_system_call("echo fmos=.false.      >> ed.in")
      endif
      if (fmos_fluc) then
         call utils_system_call("echo fmos_fluc=.true.  >> ed.in")
      else
         call utils_system_call("echo fmos_fluc=.false. >> ed.in")
      endif
      if (fmos_hub1) then
         call utils_system_call("echo fmos_hub1=.true.  >> ed.in")
      else
         call utils_system_call("echo fmos_hub1=.false. >> ed.in")
      endif
      if (((ed_do_not_keep_previous_fit_param .or. iter_dmft == 1) .and. ncpt_two_step) .or. ncpt_two_step_all_iter) then
         call utils_system_call("echo ncpt_flag_two_step_fit=.true. >> ed.in")
      else
         call utils_system_call("echo ncpt_flag_two_step_fit=.false. >> ed.in")
      endif

      if (nproc_mpi_solver > 1) then
         call utils_system_call("echo FLAG_MPI_GREENS=1 >> ed.in")
      endif
      if (ed_star_geom) then
         call utils_system_call("echo diag_bath=.true. >> ed.in")
      endif
      call utils_system_call("echo Niter_search_max_0="//trim(adjustl(tostring(ed_nsearch)))//" >> ed.in ")
      if (cluster_dmft_green_for_self_consistence .and. ed_solver_compute_all_green_functions) then
         call utils_system_call("echo FLAG_ALL_GREEN_FUNC_COMPUTED=.true. >> ed.in")
      endif
      if (compute_ed_spin_correlation) then
         call utils_system_call(" dmft_gen_correl.out  "//TRIM(ADJUSTL(toString(channels)))//" 1", abort=.true.)
      else
         call utils_system_call(" dmft_gen_correl.out  "//TRIM(ADJUSTL(toString(channels)))//" 0", abort=.true.)
      endif
      call utils_system_call("echo fit_weight_power="//trim(adjustl(tostring(fit_weight_power)))//" >> ed.in")
      if (paramagnetic_ed == 1) then
         call utils_system_call("echo force_para_state=.true. >> ed.in")
         call utils_system_call("echo FLAG_GUP_IS_GDN=.true. >> ed.in")
      endif
      if (fit_nw == 0) then
         call utils_assert(n_frequ >= 4, 'Not enough matsubara frequencies')
         call utils_system_call("echo fit_nw="//TRIM(ADJUSTL(toString(n_frequ - nspec_frequ)))//" >> ed.in")
      else
         call utils_system_call("echo fit_nw="//TRIM(ADJUSTL(toString(fit_nw)))//" >> ed.in")
      endif

      call utils_system_call("echo min_all_bath_param="//TRIM(ADJUSTL(toString(sites_ed)))//" >> ed.in")
      if (ed_do_not_keep_previous_fit_param) then
         call utils_system_call(" rm ed.fit.param ed.fmos ")
      endif

      if (print_qc) call utils_system_call('echo print_qc=.true. >> ed.in')

      call utils_system_call("mv hf.in "//TRIM(ADJUSTL(filename)))
      call utils_system_call("mv ed.* "//TRIM(ADJUSTL(filename)))
      call utils_system_call("ls ./Sigma.000 || cp "//TRIM(ADJUSTL(dir_onetep))//"/utils/INPUTS/Sigma.000 .")
      call utils_system_call("mv Sigma.000 Trans.dat Uc.dat.* gen_atom.log info_atom_d.dat "//TRIM(ADJUSTL(filename)))

!-----------------------------------------------------------------!
      open (unit=10001, file="PARAMS", ACCESS='append')
      write (*, *) '############# CHECK ED in PARAMS #################'
      if (solver == 1) then
         write (10001, *) 'nom ', nmatsu_ctqmc
         write (10001, *) 'M   ', mcs_ctqmc
         write (10001, *) 'nf0 ', occupation_orbital
         write (10001, *) 'mu  ', mmu(1)
         write (10001, *) 'beta ', bbeta
         write (10001, *) 'U ', UU
         write (10001, *) 'J ', Jhund
         write (10001, *) 'cx ', spin_orbit
         write (10001, *) 'T ', 1.d0/bbeta
         write (10001, *) 'exe ctqmc'
         write (10001, *) 'Sig sigma_output'
         write (10001, *) 'Gf green_output'
         write (10001, *) 'Delta delta_input'
         write (10001, *) 'cix actqmc.cix'
         write (10001, *) 'nc [6, 7, 8, 9, 10]      # Impurity occupancies'
!-----------------------------------------------------------------!
      elseif (solver == 4 .or. solver == 5) then
         write (10001, *) channels
         do ii = 1, channels
            write (10001, '(200f13.7)') (eimp_ed(1, ii, jj), jj=1, channels)
         enddo
         do ii = 1, channels
            write (10001, '(200f13.7)') (eimp_ed(2, ii, jj), jj=1, channels)
         enddo

         if (.not. dmft_for_dimer) then
            U_ = UU_ren
         else
            ii = size(U_, 1)
            ii = ii/2
            U_ = 0.d0
            U_(1:ii, 1:ii) = UU_ren(1:ii, 1:ii)
            U_(ii + 1:2*ii, ii + 1:2*ii) = UU_ren(ii + 1:2*ii, ii + 1:2*ii)
            JJ_ren(ii + 1:2*ii, 1:ii) = 0.d0
            JJ_ren(1:ii, ii + 1:2*ii) = 0.d0
         endif
         do ii = 1, channels
            write (10001, '(200 f13.7)') (U_(ii, jj), jj=1, channels)
         enddo

         write (10001, *) n_frequ
         write (10001, *) nmatsu_ed
         write (10001, *) ed_real_frequ
         if (iter_dmft /= niter_dmft) ed_compute_all = .false.
         write (10001, *) ed_compute_all
         write (10001, *) paramagnetic_ed ! FLAG_ORDER_TYPE
         if ((iter_dmft == niter_dmft .or. ed_compute_retarded_every_step) .and. .not. ed_no_real_overide) then
            write (10001, *) .true.          ! compute retarded GF at the last DMFT iteration
         else
            write (10001, *) .false.         ! no retarded GF
         endif
         write (10001, *) ed_frequ_min
         write (10001, *) ed_frequ_max
         write (10001, *) mmu(1)
         write (10001, *) bbeta
         write (10001, *) UU
         write (10001, *) ed_rdelta
         write (10001, *) ed_rdelta_frequ_eta1
         write (10001, *) ed_rdelta_frequ_T
         write (10001, *) ed_rdelta_frequ_w0
         write (10001, *) iter_dmft
         write (10001, *) Jhund
         write (10001, *) use_input_delta
         write (10001, *) (impose_no_cdw .and. dmft_for_dimer)
         write (10001, *) fit_green
         if (iter_dmft == niter_dmft .or. ed_compute_retarded_every_step) nohole = .false.
         write (10001, *) nohole
         do ii = 1, channels
            write (10001, '(200 f13.7)') (JJ_ren(ii, jj), jj=1, channels)
         enddo
!-----------------------------------------------------------------!
      elseif (solver == 2 .or. solver == 3) then
         if (solver == 2) then
            write (10001, *) 'exe=nca'
         else
            write (10001, *) 'exe=oca'
         endif
         stringdat = ' Ed=['
         do i = 1, channels
            stringdat = TRIM(ADJUSTL(stringdat))//TRIM(ADJUSTL(toString(eimp_nca(1, i))))
            if (i < channels) stringdat = TRIM(ADJUSTL(stringdat))//","
         enddo
         if (paramagnetic /= 1) then
            stringdat = TRIM(ADJUSTL(stringdat))//','
            do i = 1, channels
               stringdat = TRIM(ADJUSTL(stringdat))//TRIM(ADJUSTL(toString(eimp_nca(2, i))))
               if (i < channels) stringdat = TRIM(ADJUSTL(stringdat))//","
            enddo
         endif
         stringdat = TRIM(ADJUSTL(stringdat))//']'
         write (10001, '(a112)') TRIM(ADJUSTL(stringdat))
         write (10001, *) 'nf0='//TRIM(ADJUSTL(toString(occupation_orbital)))
         write (10001, *) 'mu=0.' !and NOT mmu, we included mmu in eimp_nca, the reason is that NCA cannot read any mmu
         write (10001, *) 'beta='//TRIM(ADJUSTL(toString(bbeta)))
         write (10001, *) 'U='//TRIM(ADJUSTL(toString(UU)))
         write (10001, *) 'J='//TRIM(ADJUSTL(toString(Jhund)))
         write (10001, *) 'cx='//TRIM(ADJUSTL(toString(spin_orbit)))
         write (10001, *) 'T='//TRIM(ADJUSTL(toString(1.d0/bbeta)))
         write (10001, *) 'Delta=delta_input'
         write (10001, *) 'Ac=Ac.inp'
         write (10001, *) 'SigOut=Sigma.000'
         write (10001, *) 'Sig=Sigma.000'
         write (10001, *) 'sig=sigma_output'
         write (10001, *) 'AlocOut=Aloc.imp'
         write (10001, *) 'gloc=gloc.out'
         write (10001, *) 'StartLambda='//TRIM(ADJUSTL(toString(StartLambda)))
         write (10001, *) 'EndLambda='//TRIM(ADJUSTL(toString(EndLambda)))
         write (10001, *) 'alpha='//TRIM(ADJUSTL(toString(alpha)))
         write (10001, *) 'max_steps='//TRIM(ADJUSTL(toString(max_steps)))
         write (10001, *) 'followPeak='//TRIM(ADJUSTL(toString(followPeak)))
         write (10001, *) 'Q='//TRIM(ADJUSTL(toString(QQ)))
         write (10001, *) 'cix=out.cix'
         write (10001, *) 'nc=[4, 5, 6, 7]         #Impurity occupancies'
      else
         call utils_abort('ERROR not implemented yet')
      endif
      close (10001)
      write (*, *) '############# END PARAMS #################'
!-----------------------------------------------------------------!

      if (.not. read_cix_file) then
         call utils_system_call("cp actqmc.cix _actqmc.cix_")
      endif

      call utils_system_call("mv green_onetep_spin* Ac.inp* "//TRIM(ADJUSTL(filename))//" || echo error_no_green_onetep_spin_files ")
      call utils_system_call("mv delta_input* out.cix actqmc.cix "//TRIM(ADJUSTL(filename)))
      call utils_system_call("cp PARAMS PARAMS.oca ")
      call utils_system_call("mv PARAMS* "//TRIM(ADJUSTL(filename)))

      if (.not. ctqmc_erase_status) call utils_system_call("cp status* "//TRIM(ADJUSTL(filename)))

      if (solver == 1) then
         call utils_system_call("run_iter_ctqmc "//TRIM(ADJUSTL(filename))//" "//TRIM(ADJUSTL(toString(openmp_solver)))//" "//TRIM(ADJUSTL(toString(nproc_mpi_solver))), abort=.true.) 
      elseif (solver == 2) then
         call utils_system_call("run_iter_nca "//TRIM(ADJUSTL(filename))//" "//TRIM(ADJUSTL(toString(openmp_solver)))//" "//TRIM(ADJUSTL(toString(nproc_mpi_solver))), abort=.true.)
      elseif (solver == 3) then
         call utils_system_call("run_iter_oca "//TRIM(ADJUSTL(filename))//" "//TRIM(ADJUSTL(toString(openmp_solver)))//" "//TRIM(ADJUSTL(toString(nproc_mpi_solver))), abort=.true.)
      elseif (solver == 4) then
         call start_timer("run_iter_ed")
         call utils_system_call("run_iter_ed "//TRIM(ADJUSTL(filename))//" "//TRIM(ADJUSTL(toString(openmp_solver)))//" "//TRIM(ADJUSTL(toString(nproc_mpi_solver))), abort=.true.)
         call stop_timer("run_iter_ed")
      elseif (solver == 5) then
         call utils_system_call("run_iter_hf "//TRIM(ADJUSTL(filename))//" "//TRIM(ADJUSTL(toString(openmp_solver)))//" "//TRIM(ADJUSTL(toString(nproc_mpi_solver))), abort=.true.)
      else
         call utils_abort('Solver not defined')
      endif

      write (*, *) 'done ... now formatting data, paramagnetic = ', paramagnetic

      if (paramagnetic /= 1) then
         write (*, *) 'need to split sigma tot in sigma up and dn'
         write (*, *) 'in onetep need to read sigma up and down separately'
         call utils_system_call("split_sigma_1 " // trim(adjustl(filename_sigma(1))) // " " &
               // TRIM(ADJUSTL(toString(channels))), abort=.true.)
         call utils_system_call("split_sigma_2 " // trim(adjustl(filename_sigma(2))) // " " &
               // TRIM(ADJUSTL(toString(channels))), abort=.true.)
         call utils_system_call(" rm sigma_output ")
      else
         call utils_system_call("mv sigma_output "//trim(adjustl(filename_sigma(1))))
      endif

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine getting_S_matrix_from_Green(Smat, COEF, channels, kk_, silent, inverse)
      implicit none
      integer          :: kk_, channels, kkk
      complex(kind=DP)       :: COEF(channels), dummy(channels, channels), Smat(channels, channels)
      logical, optional :: silent, inverse

      !-------------------------------!
      !-------------------------------!
      if (cluster_dmft_green_for_self_consistence) then
         if (real_or_imaginary_solver == 2) then
            !---------------------!
            ! S    = G^-1 * 1/iw  !
            !---------------------!
            Smat = green_(kk_, n_frequ, :, :)*imi*aimag(frequ_(kk_, n_frequ))
            if (.not. present(inverse)) call invmat(channels, Smat)
            if (.not. present(silent)) call write_array(Smat, 'S matrix', short=.true., unit=6)
            !-------------------------------!
         else
            !-------------------------------!
            Smat = green_(kk_, n_frequ, :, :)*real(frequ_(kk_, n_frequ))
            if (.not. present(inverse)) call invmat(channels, Smat)
            dummy = Smat
            if (.not. present(silent)) call write_array(Smat, 'S matrix large frequ', short=.true., unit=6)
            Smat = green_(kk_, 1, :, :)*real(frequ_(kk_, 1))
            if (.not. present(inverse)) call invmat(channels, Smat)
            Smat = (Smat + dummy)/2.d0
            if (.not. present(silent)) call write_array(Smat, 'S matrix average', short=.true., unit=6)
         endif
         COEF = DIAG(Smat)
         !-------------------------------!
         !-------------------------------!
      else
         !-------------------------------!
         !-------------------------------!
         if (real_or_imaginary_solver == 2) then
            COEF(:) = aimag((1.d0/green_diag(kk_, n_frequ, :)))/aimag(frequ_(kk_, n_frequ))
         else
            !-------------------------------!
            COEF(:) = (1.d0/green_diag(kk_, n_frequ, :))/frequ_(kk_, n_frequ)
            COEF = COEF + (1.d0/green_diag(kk_, 1, :))/frequ_(kk_, 1)
            COEF = COEF/2.d0
           if (.not. present(silent)) write (*, *) 'COEF POSITIVE FREQU : ', (1.d0/green_diag(kk_, n_frequ, :))/frequ_(kk_, n_frequ)
            if (.not. present(silent)) write (*, *) 'COEF NEGATIVE FREQU : ', (1.d0/green_diag(kk_, 1, :))/frequ_(kk_, 1)
         endif
         if (present(inverse)) COEF = 1.d0/COEF
         Smat = 0.
         do kkk = 1, channels
            Smat(kkk, kkk) = COEF(kkk)
         enddo
         if (.not. present(silent)) write (*, *) 'HYBRIDIZATION COEFFICIENTS DUE TO NON ORTHOGONALITY ARE : ', COEF
      endif
      !-------------------------------!
      !-------------------------------!

      if (real_smat_coef) then
         Smat = real(Smat)
         COEF = real(COEF)
      endif

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!

   function s_mat_from_onetep(kk, channels)
      implicit none
      integer          :: channels
      complex(kind=DP)       :: s_mat_from_onetep(channels, channels)
      integer          :: kk
      s_mat_from_onetep = green_(kk, n_frequ - 3, :, :)
      if (real_smat_from_onetep) s_mat_from_onetep = real(s_mat_from_onetep)
   end function

   !====================!
   !====================!

   function eimp_mat_from_onetep(kk, channels)
      implicit none
      integer          :: channels
      complex(kind=DP)       :: eimp_mat_from_onetep(channels, channels)
      integer          :: kk
      eimp_mat_from_onetep = green_(kk, n_frequ - 4, :, :)
      if (real_eimp_from_onetep) eimp_mat_from_onetep = real(eimp_mat_from_onetep)
   end function

   !====================!
   !====================!

   function eimp_mat_from_tail(kk, channels)
      implicit none
      integer          ::  channels
      complex(kind=DP)       ::  eimp_mat_from_tail(channels, channels)
      integer          ::  kk
      eimp_mat_from_tail = green_(kk, n_frequ - 1, :, :)
      if (real_eimp_from_tail) eimp_mat_from_tail = real(eimp_mat_from_tail)
   end function

   !====================!
   !====================!

   function eimp_diag_from_tail(kk, channels)
      implicit none
      integer          ::  channels
      complex(kind=DP)       ::  eimp_diag_from_tail(channels)
      integer          ::  kk
      eimp_diag_from_tail = diag(green_(kk, n_frequ - 1, :, :))
      if (real_eimp_from_tail) eimp_diag_from_tail = real(eimp_diag_from_tail)
   end function

   !====================!
   !====================!

   function sproj_mat_from_onetep(kk, channels)
      implicit none
      integer          :: channels
      complex(kind=DP)       :: sproj_mat_from_onetep(channels, channels)
      integer          :: kk
      sproj_mat_from_onetep = green_(kk, n_frequ - 5, :, :)
      if (real_sproj_from_onetep) sproj_mat_from_onetep = real(sproj_mat_from_onetep)
   end function

   !====================!
   !====================!
   !====================!
   !====================!

   subroutine get_hybridizations(kk_)

      use common_def, only: utils_abort

      implicit none
      logical       :: check
      integer       :: kk_, kkk
      real(kind=DP)       :: tt, tmp_read(1 + 2*channels), tmp1, tmp2, ah(channels), bh(channels), vvv(channels), temp(1 + 2*channels)
      real(kind=DP)       :: ahl(channels), bhl(channels), ahlmat(channels, channels), bhlmat(channels, channels)
      real(kind=DP)       :: ahmat(channels, channels), bhmat(channels, channels)
      complex(kind=DP)    :: COEF(channels), dummy(channels, channels), Smat(channels, channels), iwplusmu(channels, channels)
      complex(kind=DP)    :: himp(channels, channels), Simp(channels, channels), sinfty(channels, channels), tmp__(channels, channels)
      character(30) :: delta_filename
      integer       :: i1__, j1__

      !###############################!
      !###############################!
      !###############################!

      write (*, *) '####################################'
      write (*, *) '####################################'
      write (*, *) '####################################'
      write (*, *) 'computing high frequency corrections'

      dummy = s_mat_from_onetep(kk_, channels)
      call write_array(dummy, ' Simp =   U^+ S^-1 U ', short=.true., unit=6)
      if (maxval(abs(dummy)) > 1.d-14) then
         call invmat(channels, dummy)
         call write_array(dummy, ' Simp = ( U^+ S^-1 U )^-1 ', short=.true., unit=6)
         Simp = dummy
      else
         Simp = 0.d0
      endif
      dummy = eimp_mat_from_onetep(kk_, channels)
      call write_array(-dummy, '  U^+ [S^-1 H S^-1] U      ', short=.true., unit=6)
      dummy = matmul(matmul(Simp, dummy), Simp)
      call write_array(-dummy, ' himp = Simp ( U^+ [S^-1 H S^-1] U ) Simp  ', short=.true., unit=6)
      himp = -dummy

      dummy = sproj_mat_from_onetep(kk_, channels)
      call write_array(-dummy, '  U^+ [S^-1 Soo S^-1] U      ', short=.true., unit=6)
      dummy = matmul(matmul(Simp, dummy), Simp)
      call write_array(-dummy, ' Sproj(oo) = Simp ( U^+ [S^-1 Soo S^-1] U ) Simp ', short=.true., unit=6)
      sinfty = -dummy

      write (*, *) '####################################'
      write (*, *) '####################################'
      write (*, *) '####################################'

      !###############################!
      !###############################!
      !###############################!

      write (*, *) 'some info for checking/debugging purposes'
      do i = n_frequ - nspec_frequ, n_frequ
         if (i == n_frequ - 5) then
            write (*, *) '(for next line, if you want to use a finite frequency for comparison, turn off nmatsu_long)'
     write(*,'(a,200f13.3)') 'cor.space Sigma-Edc(w=oo) is at n_frequ-5 (sigma_mat) RE : ',(real(sigma_mat(kk_,i,k,k)),k=1,channels)
        write (*, '(a,200f13.3)') 'projected Sigma-Edc(w=oo) is at n_frequ-5 (sinfty)    RE : ', (real(sinfty(k, k)), k=1, channels)
     write(*,'(a,200f13.3)') 'cor.space Sigma-Edc(w=oo) is at n_frequ-5 (sigma_mat) IM : ',(aimag(sigma_mat(kk_,i,k,k)),k=1,channels)
       write (*, '(a,200f13.3)') 'projected Sigma-Edc(w=oo) is at n_frequ-5 (sinfty)    IM : ', (aimag(sinfty(k, k)), k=1, channels)
         endif
  if(i==n_frequ-1) write(*,'(a,200f13.3)') 'self energy is at oo (nfreq-1 - check ) RE : ',(real(sigma_mat(kk_,i,k,k)),k=1,channels)
  if(i==n_frequ-5) write(*,'(a,200f13.3)') 'self energy is at oo (nfreq-5 - check ) RE : ',(real(sigma_mat(kk_,i,k,k)),k=1,channels)
  if(i==n_frequ)   write(*,'(a,200f13.3)') 'self energy is at n_frequ (E_DC)        RE : ',(real(sigma_mat(kk_,i,k,k)),k=1,channels)
  if(i==n_frequ-1) write(*,'(a,200f13.3)') 'self energy is at oo (nfreq-1 - check ) IM : ',(aimag(sigma_mat(kk_,i,k,k)),k=1,channels)
  if(i==n_frequ-5) write(*,'(a,200f13.3)') 'self energy is at oo (nfreq-5 - check ) IM : ',(aimag(sigma_mat(kk_,i,k,k)),k=1,channels)
  if(i==n_frequ)   write(*,'(a,200f13.3)') 'self energy is at n_frequ (E_DC)        IM : ',(aimag(sigma_mat(kk_,i,k,k)),k=1,channels)
      enddo

      write (*, *) '------ done -----'
      write (*, *) 'dumping diagonal part of the green function'
      open (unit=10102, file="green_onetep_spin"//TRIM(ADJUSTL(toString(kk_))))
      do i = 1, n_frequ - nspec_frequ
       write (10102, '(200f13.5)') myfrequ(kk_, i), (real(green_diag(kk_, i, kkk)), aimag(green_diag(kk_, i, kkk)), kkk=1, channels)
      enddo
      close (10102)
      write (*, *) '------ done -----'

      !###############################!
      !###############################!
      !###############################!

      call getting_S_matrix_from_Green(Smat, COEF, channels, kk_)

      if (cluster_dmft_green_for_self_consistence .and. use_simp_from_onetep) then
         write (*, *) 'USING ONETEP Simp instead of extracted matrix from projected GF'
         Smat = Simp
         COEF = DIAG(Smat)
      endif

      !###############################!
      !###############################!
      !###############################!

      denstarget(kk_) = 2.0*sum(real(green_diag(kk_, 1:n_frequ - nspec_frequ - 1, :)))/bbeta + 0.5d0*dble(channels)
      write (*, *) 'WARNING : HUB1 SOLVER WILL NOT WORK PROPERLY IN NON ORTHOGONAL BASIS - please use rotate ortho flag'
      write (*, *) 'DENSTARGET OBTAINED BY SUMMING GF              : ', denstarget(kk_)
      write (*, *) 'COMPARISON - DENSTARGET OBTAINED WITH TAIL     : ', &
             & sum((/(densfourier(n_frequ - nspec_frequ - 1, bbeta, aimag(frequ_(kk_, 1:n_frequ - nspec_frequ - 1)),&
                                   & green_diag(kk_, 1:n_frequ - nspec_frequ - 1, j)), j=1, channels)/))

      if (kk_ == 2 .or. paramagnetic == 1) then
         open (unit=10015, file='ed.hub1.target')
         if (paramagnetic == 1) then
            write (10015, *) denstarget(1)
            write (10015, *) denstarget(1)
         else
            write (10015, *) denstarget(1)
            write (10015, *) denstarget(2)
         endif
         write (10015, *) fmos_hub_iter_mu
         write (10015, *) fmos_hub_range_mu
         close (10015)
      endif

      if (.not. rotate_to_ortho_basis) then
         open (unit=10012, file="green_onetep_spin_check_convergence"//TRIM(ADJUSTL(toString(kk_))), form='unformatted')
      else
         open (unit=10012, file="green_onetep_spin_check_convergence_ren"//TRIM(ADJUSTL(toString(kk_))), form='unformatted')
      endif
  vvv(:) = real(REAL(COEF(:))*(frequ_(kk_,n_frequ-nspec_frequ)+mmu(kk_)) - diag(sigma_mat(kk_,n_frequ-nspec_frequ,:,:)) - 1.d0/green_diag(kk_,n_frequ-nspec_frequ,:))
      write (10012) n_frequ - nspec_frequ, channels
      if (.not. rotate_to_ortho_basis) then
         write (10012) vvv, mmu(kk_), COEF
      else
         write (*, *) 'Writing in file nfrequ : ', n_frequ - nspec_frequ
         write (*, *) 'Writing in file (chan) : ', channels
         write (*, *) 'Writing in file (vvv)  : ', vvv
         write (*, *) 'Writing in file (mmu)  : ', mmu(kk_)
         write (*, *) 'Writing in file (Simp) : ', diag(Simp_back(kk_, :, :))
         write (*, *) 'Writing in file (gp)   : ', diag(Simp_root_gp(kk_, :, :))
         write (10012) vvv
         write (10012) mmu(kk_)
         write (10012) diag(Simp_back(kk_, :, :))
         write (10012) diag(Simp_root_gp(kk_, :, :))
      endif
      do i = 1, n_frequ - nspec_frequ
         if (.not. rotate_to_ortho_basis) then
            write (10012) frequ_(kk_, i), green_diag(kk_, i, :), diag(sigma_mat(kk_, i, :, :))
         else
            write (10012) frequ_(kk_, i), green_diag(kk_, i, :)
         endif
      enddo
      close (10012)

      if (.not. rotate_to_ortho_basis) then
         open (unit=10113, file="green_onetep_spin_full_check_convergence"//TRIM(ADJUSTL(toString(kk_))), form='unformatted')
      else
         open (unit=10113, file="green_onetep_spin_full_check_convergence_ren"//TRIM(ADJUSTL(toString(kk_))), form='unformatted')
      endif
      write (10113) n_frequ - nspec_frequ, channels
      if (.not. rotate_to_ortho_basis) then
         write (10113) mmu(kk_), Smat
      else
         write (10113) mmu(kk_), Simp_back(kk_, :, :), Simp_root_gp(kk_, :, :)
      endif
      do i = 1, n_frequ - nspec_frequ
         if (.not. rotate_to_ortho_basis) then
            write (10113) frequ_(kk_, i), green_(kk_, i, :, :), sigma_mat(kk_, i, :, :)
         else
            write (10113) frequ_(kk_, i), green_(kk_, i, :, :)
         endif
      enddo
      close (10113)

      write (*, *) 'DENSITY MATRIX AT FERMI LEVEL'
      write (*, *) 'obtained from G_proj from onetep, with the overlap matrix Smat'
      i = my_zero_frequ(kk_)
      call write_array(matmul(green_(kk_, i, :, :), Smat), 'Gproj(w=0)*Smat', unit=6)
      call write_array(green_(kk_, i, :, :), 'Gproj(w=0)', unit=6)

      !###############################!
      !###############################!
      !###############################!

      !---------------------------------------!
      ! compute Hybridization of the impurity !
      !---------------------------------------!

      if (cluster_dmft_green_for_self_consistence) then
         do i = 1, n_frequ
            dummy = green_(kk_, i, :, :)

            if (flag_blank_out_green_offdiag_for_testing) then
               dummy = 0.d0
               do k = 1, channels
                  dummy(k, k) = green_(kk_, i, k, k)
               enddo
            endif
            call invmat(channels, dummy)
            iwplusmu = 0.
            do k = 1, channels
               iwplusmu(k, k) = frequ_(kk_, i) + mmu(kk_)
            enddo

            green_(kk_, i, :, :) = MATMUL(iwplusmu, Smat) - dummy - sigma_mat(kk_, i, :, :)

            if (flag_blank_out_green_offdiag_for_testing) then
               do k = 1, channels
                  do kkk = 1, channels
                     if (k /= kkk) green_(kk_, i, k, kkk) = 0.d0
                  enddo
               enddo
            endif

            green_diag(kk_, i, :) = diag(green_(kk_, i, :, :))
         enddo
      else
         do i = 1, n_frequ
      green_diag(kk_, i, :) = REAL(COEF(:))*(frequ_(kk_, i) + mmu(kk_)) - diag(sigma_mat(kk_, i, :, :)) - 1.d0/green_diag(kk_, i, :)
            green_(kk_, i, :, :) = 0.d0
            do k = 1, channels
               green_(kk_, i, k, k) = green_diag(kk_, i, k)
            enddo
         enddo
      endif

      write (*, *) 'HYBRID AT FERMI LEVEL'
      write (*, *) 'obtained from G_proj from onetep, with the overlap matrix Smat'
      call write_array(Smat, ' Smat ', unit=6)
      i = my_zero_frequ(kk_)
      call write_array(green_(kk_, i, :, :), 'Hybrid(w=0)', unit=6)
      call write_array(sigma_mat(kk_, i, :, :), ' Sigma(w=0)', unit=6)
      i = min(10, n_frequ)
      call write_array(sigma_mat(kk_, i, :, :), ' Sigma(iw_{n=10})', unit=6)
      call write_array(green_(kk_, i, :, :), 'Hybrid(iw_{n=10})', unit=6)
      call write_array(green_(kk_, i, :, :) - green_(kk_, n_frequ, :, :), 'Hybrid(iw_{n=10})-Hybrid(w=oo)', unit=6)

      !###############################!
      !###############################!
      !###############################!

      !-----------------------------------------------------------------------------!
      ! Compute Impurity Green Function (should converge to Lattice green function) !
      !-----------------------------------------------------------------------------!

      open (unit=340, file="green_onetep_spin_imp"//TRIM(ADJUSTL(toString(kk_))))
      do i = 1, n_frequ - nspec_frequ
         dummy = green_(kk_, i, :, :)
         iwplusmu = 0.
         do k = 1, channels
            iwplusmu(k, k) = frequ_(kk_, i) + mmu(kk_)
         enddo
         dummy = MATMUL(iwplusmu, Smat) - dummy - sigma_mat(kk_, i, :, :)
         if (.not. cluster_dmft_green_for_self_consistence) then
            do kkk = 1, channels
               dummy(kkk, kkk) = 1.d0/dummy(kkk, kkk)
            enddo
         else
            call invmat(channels, dummy)
         endif
         temp(1) = myfrequ(kk_, i)
         do kkk = 1, channels
            temp(1 + 2*kkk - 1) = real(dummy(kkk, kkk))
            temp(1 + 2*kkk) = aimag(dummy(kkk, kkk))
         enddo
         write (340, '(200f13.5)') (temp(kkk), kkk=1, 1 + 2*channels)
      enddo
      close (340)

      !###############################!
      !###############################!
      !###############################!

      eimp_ed(kk_, :, :) = eimp_mat_from_tail(kk_, channels)
      eimp(kk_, :) = eimp_diag_from_tail(kk_, channels)
      eimp_nca(kk_, :) = eimp(kk_, :) - mmu(kk_)

      if (second_order_correction_to_eimp) then
         do i = 1, channels
            do j = 1, channels
               tt = eimp_ed(kk_, i, j)
               call extract_Eimp_from_tail_of_hybridization(frequ_(kk_, :), green_(kk_, :, i, j), tt, tmp1, tmp2)
               write (*, '(a,2i4,2f12.4)') 'CORRECTING eimp with second order terms : ', i, j, eimp_ed(kk_, i, j) - tt
               eimp_ed(kk_, i, j) = tt
               write (*, *) 'a,b coefs : ', tmp1, tmp2
            enddo
            call extract_Eimp_from_tail_of_hybridization(frequ_(kk_, :), green_diag(kk_, :, i), eimp(kk_, i), tmp1, tmp2)
            eimp_nca(kk_, i) = eimp(kk_, i) - mmu(kk_)
         enddo
      call write_array( real(eimp_ed(kk_,:,:) - green_(kk_,n_frequ-1,:,:)  ), ' Real correction to eimp'//trim(adjustl(tostring(kk_))),short=.true.,unit=6) 
         call write_array(real(eimp(kk_, :) - green_diag(kk_, n_frequ - 1, :)), ' Real correction to eimp_ ', short=.true., unit=6)
      call write_array( aimag(eimp_ed(kk_,:,:) - green_(kk_,n_frequ-1,:,:)  ), ' Im correction to eimp'//trim(adjustl(tostring(kk_))),short=.true.,unit=6)
         call write_array(aimag(eimp(kk_, :) - green_diag(kk_, n_frequ - 1, :)), ' Im correction to eimp_ ', short=.true., unit=6)

      endif

      dummy = 0.0d0
      do i = 1, channels
         dummy(i, i) = double_counting(kk_, i)
      enddo
      if (double_counting_zero_self) Then
         do i = 1, channels
            dummy(i, i) = sigma_mat(kk_, n_frequ, i, i)
         enddo
      endif

      if (rotate_to_ortho_basis) then
         dummy = MATMUL(MATMUL(transpose(Simp_root_m(kk_, :, :)), dummy), Simp_root_m(kk_, :, :))
      endif

      if (use_eimp_from_onetep) then
         if (use_eimp_from_onetep_with_sigma_cor) then
            eimp_ed(kk_, :, :) = himp + sinfty - real(sigma_mat(kk_, n_frequ - 1, :, :))
         else
            eimp_ed(kk_, :, :) = himp - real(dummy)
         endif
         eimp(kk_, :) = diag(eimp_ed(kk_, :, :))
         eimp_nca(kk_, :) = eimp(kk_, :) - mmu(kk_)
      endif

      write (*, *) '#####################################################################################'
      write (*, *) 'The following gives an idea how Sigma(oo) upfolded contributes to himp when projected'
      write (*, *) 'it should not contribute much, as eimp should not be affected by Sigma.              '
      call write_array(eimp_ed(kk_, :, :) - (real(himp) - real(dummy)),&
                                 & ' eimp - (himp - Edc)        for spin '//trim(adjustl(tostring(kk_))), short=.true., unit=6)
      call write_array(eimp_ed(kk_, :, :) - (real(himp + sinfty) - real(sigma_mat(kk_, n_frequ - 1, :, :))),&
                                 & ' eimp - (himp+sinfty-sigma) for spin '//trim(adjustl(tostring(kk_))), short=.true., unit=6)
      call write_array((real(sinfty) + real(dummy)) - real(sigma_mat(kk_, n_frequ - 1, :, :)),  &
                                 & '(s_infty+Edc) - Sig_imp(oo) for spin '//trim(adjustl(tostring(kk_))), short=.true., unit=6)
      write (*, *) '#####################################################################################'

      !###############################!
      !###############################!
      !###############################!

      if (real_or_imaginary_solver == 2) then
         bbeta = 1.d0/(aimag(frequ_(kk_, 1))/pi)
         write (*, *) 'temperature is (obtained from matsubara frequ) [Hartree] : ', 1.d0/bbeta
         write (*, *) 'beta is                                                  : ', bbeta
      endif

      do i = 1, n_frequ
         green_diag(kk_, i, :) = green_diag(kk_, i, :) - eimp(kk_, :)
         green_(kk_, i, :, :) = green_(kk_, i, :, :) - eimp_ed(kk_, :, :)
      enddo

      !###############################!
      !###############################!
      !###############################!

      do i = n_frequ - nspec_frequ + 1, n_frequ
         frequ_(kk_, i) = frequ_(kk_, i - 1) + (frequ_(kk_, i - 1) - frequ_(kk_, i - 2))
         green_diag(kk_, i, :) = green_diag(kk_, i - 1, :)
         green_(kk_, i, :, :) = green_(kk_, i - 1, :, :)
      enddo

      if (real_or_imaginary_solver == 1) then
         green_diag(kk_, 1, :) = green_diag(kk_, 2, :)
         green_(kk_, 1, :, :) = green_(kk_, 2, :, :)
         frequ_(kk_, 1) = frequ_(kk_, 2) + (frequ_(kk_, 2) - frequ_(kk_, 3))
      endif

      !###############################!
      !###############################!
      !###############################!

      if (flag_correct_eimp_spin_orbit) then
         eimp(kk_, :) = eimp(kk_, :) + 2.d0*spin_orbit
         eimp_nca(kk_, :) = eimp_nca(kk_, :) + 2.d0*spin_orbit
         do i = 1, channels
            eimp_ed(kk_, i, i) = eimp_ed(kk_, i, i) + 2.d0*spin_orbit
         enddo
      endif
      if (kk_ == 2 .and. (iter_dmft == 1 .or. amp_slight_sym_breaking_all_iter)) then
         do i = 1, channels
            eimp_ed(kk_, i, i) = eimp_ed(kk_, i, i) - amp_slight_sym_breaking
         enddo
         eimp(kk_, :) = eimp(kk_, :) - amp_slight_sym_breaking
         eimp_nca(kk_, :) = eimp_nca(kk_, :) - amp_slight_sym_breaking
      endif

      if (abs(spin_orbit) > 1.d-5 .and. (solver == 4 .or. solver == 5)) then
         call utils_abort('spin orbit not implemented in ED solver or in HF')
      endif

      write (*, *) '---------------------------------------------------------'
      write (*, *) ' paramagnetic          : ', paramagnetic
      write (*, *) ' EIMP IS               : ', eimp(kk_, :)
      write (*, *) ' EIMP NCA IS           : ', eimp_nca(kk_, :)
      write (*, *) ' COEF                  : ', real(COEF)
      write (*, *) ' frequ                 : ', frequ_(kk_, n_frequ - 1)
      write (*, *) ' sigma_diag last frequ : ', diag(sigma_mat(kk_, n_frequ - 1, :, :))
      write (*, *) ' green_diag last frequ : ', green_diag(kk_, n_frequ - 1, :)
      write (*, *) ' green_ last frequ     : ', diag(green_(kk_, n_frequ - 1, :, :))
      write (*, *) '---------------------------------------------------------'

      !###############################!
      !###############################!
      !###############################!

      write (*, *) 'starting test on diagonal part of hybridization'
      do i = 1, n_frequ
         if (maxval(abs(green_diag(kk_, i, :) - diag(green_(kk_, i, :, :)))) > 1.d-4) then
            write (*, *) 'ERROR DIAGONAL HYBRIDIZATION NOT MATCHING DIAG_ARRAY'
            write (*, *) 'diag green matrix (Re) : ', diag(real(green_(kk_, i, :, :), kind=8))
            write (*, *) 'diag array        (Re) : ', real(green_diag(kk_, i, :))
            call utils_abort('Diagonal hybridization does not match diagonal array')
         endif
      enddo

      write (*, *) ' hybridization function '
      if (paramagnetic == 1) then
         delta_filename = 'delta_input'
      else
         delta_filename = 'delta_input'//TRIM(ADJUSTL(toString(kk_)))
      endif
      open (unit=10001, file=trim(adjustl(delta_filename)))
      do i = 1, n_frequ
         write (10001, '(200f19.9)') myfrequ(kk_, i), (real(green_diag(kk_, i, k)), aimag(green_diag(kk_, i, k)), k=1, channels)
      enddo
      close (10001)

      if (kk_ == 1) open (unit=10001, file='delta_input_full_1', form='unformatted')
      if (kk_ == 2) open (unit=10001, file='delta_input_full_2', form='unformatted')
      if (tail_linear_scaling < 0.) then
         do i = 1, n_frequ
            write (10001) myfrequ(kk_, i), green_(kk_, i, :, :)
         enddo
      else
         write (*, *) 'SETTING UP THE TAIL PART OF THE GF FOR LINEAR SCALING'
         do i = 1, channels
            do j = 1, channels
               call extract_Eimp_from_tail_of_hybridization(frequ_(kk_, :), green_(kk_, :, i, j), tmp1, ahmat(i, j), bhmat(i, j))
               if (real_or_imaginary_solver == 1) then
   call extract_Eimp_from_tail_of_hybridization(frequ_(kk_, :), green_(kk_, :, i, j), tmp1, ahlmat(i, j), bhlmat(i, j), left=.true.)
               endif
            enddo
         enddo
         do i = 1, n_frequ
            tt = myfrequ(kk_, i)
            if (abs(tt) < abs(tail_linear_scaling)) then
               write (10001) tt, green_(kk_, i, :, :)
            else
               do i1__ = 1, channels
                  do j1__ = 1, channels
                     tmp__(i1__, j1__) = extended_delta_el(kk_, i, ahmat(i1__, j1__), bhmat(i1__, j1__))
                  enddo
               enddo
               write (10001) tt, tmp__
            endif
         enddo
      endif
      close (10001)

      if ((solver == 1 .or. solver == 2 .or. solver == 3) .and. nmatsu_long > 0) then
         call utils_system_call("cp "//trim(adjustl(delta_filename))//" "//trim(adjustl(delta_filename))//"_backup")
         open (unit=10001, file=trim(adjustl(delta_filename)))
         do i = 1, channels
            call extract_Eimp_from_tail_of_hybridization(frequ_(kk_, :), green_diag(kk_, :, i), tmp1, ah(i), bh(i))
            if (real_or_imaginary_solver == 1) then
              call extract_Eimp_from_tail_of_hybridization(frequ_(kk_, :), green_diag(kk_, :, i), tmp1, ahl(i), bhl(i), left=.true.)
            endif
         enddo
         write (*, *) 'extending hybridization with coefficients : '
         write (*, '(a,200f10.3)') 'Ah = ', ah
         write (*, '(a,200f10.3)') 'Bh = ', bh

         if (real_or_imaginary_solver == 1) then
            j = 1
            do i = 1, nmatsu_long
               tt = extended_delta_real_frequ(kk_, i)
               if (tt >= myfrequ(kk_, 2) .and. tt <= myfrequ(kk_, n_frequ - nspec_frequ)) then
                  j = j + 1
             write (10001, '(200f19.9)') myfrequ(kk_, j), (real(green_diag(kk_, j, k)), aimag(green_diag(kk_, j, k)), k=1, channels)
               else
                  if (tt < myfrequ(kk_, 2)) then
                     write (10001, '(200f19.9)') extended_delta(kk_, i, ahl, bhl)
                  else
                     write (10001, '(200f19.9)') extended_delta(kk_, i, ah, bh)
                  endif
               endif
            enddo
         else
            do i = 1, nmatsu_long
               if (i <= n_frequ - nspec_frequ) then
             write (10001, '(200f19.9)') myfrequ(kk_, i), (real(green_diag(kk_, i, k)), aimag(green_diag(kk_, i, k)), k=1, channels)
               else
                  write (10001, '(200f19.9)') extended_delta(kk_, i, ah, bh)
               endif
            enddo
         endif
         close (10001)
      endif

      if (paramagnetic == 1) then
         open (unit=10001, file='Ac.inp')
      else
         open (unit=10001, file='Ac.inp'//TRIM(ADJUSTL(toString(kk_))))
      endif
      open (unit=10002, file=trim(adjustl(delta_filename)))
      do
         read (10002, *, end=5678) (tmp_read(i), i=1, 2*channels + 1)
         write (10001, '(200f13.6)') tmp_read(1), (-tmp_read(2*i + 1)/pi, i=1, channels)
      enddo
5678  continue
      close (10001)
      close (10002)

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine collect_sigmas(kk_)
      implicit none
      logical :: check
      integer :: kk_, i

      write (*, *) 'reading sigma output of the previous dmft step'

      check = .false.
      INQUIRE (file=trim(Adjustl(filename_sigma(kk_))), EXIST=check)
      if (check) then
         write (*, *) '======== SIGMA FILE EXISTS ========='
         open (unit=10002, file=trim(Adjustl(filename_sigma(kk_))), form='unformatted')
         do i = 1, n_frequ
            read (10002) sigma_mat(kk_, i, :, :)
         enddo
         close (10002)
         write (*, *) 'moving sigma file to : ', TRIM(ADJUSTL(filename))//"/"//trim(adjustl(filename_sigma(kk_)))//".backup"
         call utils_system_call("mv " // trim(adjustl(filename_sigma(kk_))) // " " // TRIM(ADJUSTL(filename)) &
               // "/" // trim(adjustl(filename_sigma(kk_))) // ".backup")
      else
         write (*, *) '======== SIGMA FILE DOES NOT EXISTS ========='
         sigma_mat(kk_, :, :, :) = 0.d0
      endif

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine init_run

      use common_def, only: utils_assert, utils_abort

      type(string) :: cc_
      integer      :: kk_, kkk

      call utils_system_call("get_onetep_env", abort=.true.)
      open (unit=10001, file='dir_onetep')
      read (10001, '(a)') dir_onetep
      write (*, *) 'My TOSCAM base directory is ', TRIM(ADJUSTL(dir_onetep))
      close (10001, status='delete')

      num = COMMAND_ARGUMENT_COUNT()
      if (num == 2) then
         paramagnetic = 1
         write (*, *) '###### RUNNING PARAMAGNETIC CASE ######'
         call utils_assert(abs(spin_orbit) < 1.d-3, 'Spin orbit coupling cannot be done with paramagnetic case')
      elseif (num == 3) then
         write (*, *) '###### RUNNING NON-PARAMAGNETIC CASE ######'
         paramagnetic = 2
      else
         write (*, *) 'wrong number of arguments'
         write (*, *) 'NUM : ', num
         write (*, *) 'ARGUMENTS : ', value(1:num)
         call utils_abort('Incorrect number of arguments')
      endif

      do i = 1, num
         call GET_COMMAND_ARGUMENT(i, value(i), len__, status__)
      enddo
      write (*, *) 'first  argument is the filenames of the green functions : ', value(1:num - 1)
      write (*, *) 'second argument is the dmft iteration number            : ', value(num)
      if (paramagnetic == 1) then
         call utils_assert(num == 2, 'Expecting 2 arguments for paramagnetic case')
      else
         call utils_assert(num == 3, 'Expecting 3 arguments for non-paramagnetic case')
      endif

      write (*, *) 'SOLVER CHOICE REAL : ', real_solver
      write (*, *) 'SOLVER CHOICE IM   : ', im_solver

      paramagnetic_ed = 2
      if ((im_solver == 5 .or. im_solver == 4 .or. im_solver == 5) .and. paramagnetic == 1) then
         paramagnetic_ed = 1
         paramagnetic = 2
         !copy GREENS
         cc_ = TRIM(ADJUSTL(value(1)))
         call replace_in_string(cc_, '_1', '_2', 'last')
         value(100) = cc_
         write (*, *) 'paramagnetic ED, copy input spin up to spin dn'
         write (*, *) 'filename : ', value(100)
         !UPDATE COMMAND LINE INPUTS
         value(3) = value(2)
         value(2) = TRIM(ADJUSTL(value(100)))
         num = 3
         write (*, *) 'ARGUMENTS : ', TRIM(ADJUSTL(value(1))), " ", TRIM(ADJUSTL(value(2))), " ", TRIM(ADJUSTL(value(3)))
         call utils_system_call(" cp  "//TRIM(ADJUSTL(value(1)))//" "//TRIM(ADJUSTL(value(100))))
         !copy SIGMAS
         write (*, *) 'copy sigmas ...'
         cc_ = TRIM(ADJUSTL(value(1)))
         call replace_in_string(cc_, 'green', 'sigma', 'first')
         value(100) = cc_
         call replace_in_string(cc_, '_1', '_2', 'last')
         value(99) = cc_
         write (*, *) 'copying sigmas, line : ', " cp  "//TRIM(ADJUSTL(value(100)))//" "//TRIM(ADJUSTL(value(99)))
         call utils_system_call(" cp  "//TRIM(ADJUSTL(value(100)))//" "//TRIM(ADJUSTL(value(99))))
      endif

      kkk = 1; if (paramagnetic /= 1) kkk = 2
      do kk_ = 1, kkk
         filename_sigma(kk_) = trim(adjustl(value(kk_)))
         cc_ = trim(adjustl(filename_sigma(kk_)))
         call replace_in_string(cc_, 'green', 'sigma', 'first')
         filename_sigma(kk_) = cc_
         write (*, *) 'name of the self energy file : ', filename_sigma(kk_)
      enddo

      cc_ = trim(Adjustl(filename_sigma(1)))
      call replace_in_string(cc_, 'sigma', 'edc', 'first')
      filename_edc = cc_

      write (*, *) 'output double counting filename is : ', TRIM(ADJUSTL(filename_edc))

      write (*, *) '==========================================='
      write (*, *) 'getting rid of other sigmas present in dir'
      if (.not. verysilent) then
      if (kkk == 1) then
         call utils_system_call(" get_rid_other_sigmas.out "//TRIM(ADJUSTL(filename_sigma(1))), abort=.true.)
      elseif (kkk == 2) then
         call utils_system_call(" get_rid_other_sigmas.out "//TRIM(ADJUSTL(filename_sigma(1)))//"  "//TRIM(ADJUSTL(filename_sigma(2))), abort=.true.)
      else
         call utils_abort('Error when getting rid of sigmas; this case should not arise')
      endif
      endif
      write (*, *) '==========================================='

      iter_dmft = StrInt2(value(num))
      filename = "dir_"//TRIM(ADJUSTL(value(1)))//'_iter'//TRIM(ADJUSTL(toString(iter_dmft)))
      call utils_system_call("ls -F "//TRIM(ADJUSTL(filename))//" && rm -r "//TRIM(ADJUSTL(filename)))
      call utils_system_call("mkdir "//TRIM(ADJUSTL(filename)))
      call utils_system_call("mkdir "//TRIM(ADJUSTL(filename))//"/ED_out    ")
      call utils_system_call("mkdir "//TRIM(ADJUSTL(filename))//"/CTQMC_out ")
      call utils_system_call("mkdir "//TRIM(ADJUSTL(filename))//"/NCA_out   ")
      call utils_system_call("mkdir "//TRIM(ADJUSTL(filename))//"/OCA_out   ")
      call utils_system_call("mkdir "//TRIM(ADJUSTL(filename))//"/AGR       ")
      call utils_system_call("mkdir "//TRIM(ADJUSTL(filename))//"/PNG       ")

      write (*, *) 'converted to char-to-integer dmft iteration is : ', iter_dmft
      write (*, *) 'PARAMAGNETIC                                   : ', paramagnetic
      write (*, *) 'DMFT ITERATION                                 : ', iter_dmft
      write (*, *) 'PREPARING RUN FOR DIRECTORY                    : ', trim(adjustl(filename))

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine collect_greens(kk_)

      use common_def, only: utils_unit, utils_assert, utils_abort

      implicit none
      real(kind=DP)                :: ttemp
      integer                :: klm, ien, pub_dmft_points, uuu(1), kk_
      integer                :: i_, iii, jjj
      complex(kind=DP), allocatable :: Qmat(:, :), Rmat(:, :)
      logical                :: check
      integer                :: funit

      inquire(file=trim(adjustl(value(kk_))), exist=check)

      call utils_assert(check, 'Error in collect_greens: ' &
            // trim(adjustl(value(kk_))) // ' file not found')

      funit = utils_unit()

      open (unit=funit, file=trim(adjustl(value(kk_))), form='unformatted')
      j = 0
      do
         read (funit, end=96) ien, pub_dmft_points, ttemp, channels, channelsb, atom, atomb
         if (.not. allocated(green_temp)) then
            allocate (green_temp(channels, channels),                     &
                   &  mmat(channels, channels),                           &
                   &  occupation_matrix(channels, channels),              &
                   &  occupation_matrixup(channels, channels),            &
                   &  occupation_matrixdn(channels, channels),            &
                   &  ddiag(channels, channels),                          &
                   &  diagdens(2, channels), rotation(channels, channels, 2),&
                   &  double_counting(2, channels))
            double_counting = 0.d0
         endif
         read (funit, end=96) occupation_matrixup, occupation_matrixdn
         read (funit, end=96) mmu(kk_), frequ, green_temp
         j = j + 1
      enddo
96    continue
      rewind (funit)
      n_frequ = j
      write (*, *) '---------------> N_FREQU = ', n_frequ
      if (.not. allocated(green_)) then
         allocate (green_(2, n_frequ, channels, channels), frequ_(2, n_frequ), green_diag(2, n_frequ, channels), &
              & sigma_mat(2, n_frequ, channels, channels), eimp(2, channels), eimp_nca(2, channels), eimp_ed(2, channels, channels))
         allocate (JJ_ren(channels, channels), UU_ren(channels, channels), edcmatsu(2, channels, channels))
         allocate (JJ_ren0(channels, channels), UU_ren0(channels, channels))
         allocate (Zimp_ren_p(2, channels, channels), Zimp_ren_m(2, channels, channels))
     allocate(Simp_rot(2,channels,channels),Simp_back(2,channels,channels),Simp_m1(2,channels,channels),Simp_root_m(2,channels,channels),Simp_root_p(2,channels,channels))
         allocate (Simp_root_gm(2, channels, channels), Simp_root_gp(2, channels, channels))
         eimp = 0.; eimp_ed = 0.; 
      endif

      !=====================================================================!
      !=====================================================================!
      !=====================================================================!
      !=====================================================================!
      !=====================================================================!
      !=====================================================================!
      !=====================================================================!
      !=====================================================================!

      j = 0
      do
         read (funit, end=66) ien, pub_dmft_points, ttemp, channels, channelsb, atom, atomb
         call utils_assert(channelsb == channels, 'channels/=channelsb not yet implemented')
         if ((channels /= 3 .and. channels /= 5 .and. channels /= 7) .and. (solver == 1 .or. solver == 2 .or. solver == 3)) then
            write (*, *) 'CHANNELS =  ', channels
            write (*, *) 'ERROR not implemented for solver = ', solver
            call utils_abort('Mismatch between channels and solver')
         endif
         bbeta = 1.d0/abs(ttemp)
         if (ttemp < 0.) then
            if (j == 0) write (*, *) 'REAL AXIS CALCULATIONS ....'
            solver = real_solver
            real_or_imaginary_solver = 1
         else
            if (j == 0) write (*, *) 'MATSUBARA AXIS CALCULATIONS ....'
            solver = im_solver
            real_or_imaginary_solver = 2
         endif

         if (j == 0) write (*, *) 'BETA OBTAINED FROM ONETEP IS [Hartree]: ', bbeta
         if (j == 0) write (*, *) 'REAL AXIS OR MATSUBARA CALCULATIONS   : ', ttemp
         if (j == 0) write (*, *) 'channels                              : ', channels

         if (allocated(green_temp) .and. size(green_temp, 1) /= channels) then
            write (*, *) 'ERROR something wrong in dmft_one_iteration'
            call utils_abort('Mismatched array sizes')
         endif

         read (funit, end=66) occupation_matrixup, occupation_matrixdn
         read (funit, end=66) mmu(kk_), frequ, green_temp

         if (flag_symmetrize_green) then
            occupation_matrixup = (occupation_matrixup + transpose(occupation_matrixup))/2.d0
            occupation_matrixdn = (occupation_matrixdn + transpose(occupation_matrixdn))/2.d0
         endif

         occupation_matrix = occupation_matrixup + occupation_matrixdn

         if (j == 0) then
            if (dmft_for_dimer) then
               klm = channels/2
               densupdn(1) = sum(diag(occupation_matrixup(1:klm, 1:klm)))
               densupdn(2) = sum(diag(occupation_matrixdn(1:klm, 1:klm)))
            else
               densupdn(1) = sum(diag(occupation_matrixup))
               densupdn(2) = sum(diag(occupation_matrixdn))
            endif

            occupation_orbital = sum(densupdn)

            if (kk_ == 1) call define_double_counting

            write (*, *) 'ien,pub_dmft_points   : ', ien, pub_dmft_points
            if (ien == 1) then
               write (*, *) 'bare U                : ', UU
               write (*, *) 'averaged U            : ', Uaverage
               write (*, *) 'atom number           : ', atom
               write (*, *) 'chem.pot(from onetep) : ', mmu(kk_)
               write (*, *) 'double counting up    : ', double_counting(1, :)
               write (*, *) 'double counting dn    : ', double_counting(2, :)
               write (*, *) 'size of matrices      : ', channels
               write (*, *) 'frequency             : ', frequ
               write (*, *) 'orbital densities tot : ', diag(occupation_matrix)
               write (*, *) 'orbital densities up  : ', diag(occupation_matrixup)
               write (*, *) 'orbital densities dn  : ', diag(occupation_matrixdn)
               write (*, *) 'total occupation      : ', occupation_orbital
            endif
         endif

         j = j + 1
         green_(kk_, j, :, :) = green_temp
         if (flag_symmetrize_green) then
            green_(kk_, j, :, :) = (green_(kk_, j, :, :) + transpose(green_(kk_, j, :, :)))/2.d0
         endif
         frequ_(kk_, j) = frequ - mmu(kk_)
         green_diag(kk_, j, :) = diag(green_(kk_, j, :, :))

      enddo

66    continue

      write (*, *) 'number of frequencies : ', j

      if (kk_ == 2) then
         if (abs(mmu(2) - mmu(1)) > 1.d-5) then
            write (*, *) 'WARNING : different chemical potentials for spin up and dn'
            write (*, *) '          here we put the difference in Eimp, the impurity level of the DMFT impurity'
            write (*, *) '          we do not use in this case the analytic formulas to extract Eimp, since they...'
            write (*, *) '          assume the same chemical potential for up and down spins'
            use_eimp_from_onetep = .false.
         endif
         mmu(2) = mmu(1)
      endif

      !=====================================================================!
      !=====================================================================!
      !=====================================================================!
      !=====================================================================!
      !=====================================================================!
      !=====================================================================!
      !=====================================================================!
      !=====================================================================!

      ! some info below on the green function matrix is printed to the screen

      if (rotation_scheme == 1) then
         write (*, *) 'rotation matrix to diagonalize occupation matrix is : '
         mmat = occupation_matrix
      elseif (rotation_scheme == 2 .or. rotation_scheme == 4 .or. rotation_scheme == 5) then
         uuu = maxloc((/(maxval(abs(offdiag(green_(kk_, i, :, :)))), i=1, n_frequ - nspec_frequ)/))
         mmat = green_(kk_, uuu(1), :, :)
     write(*,*) '============> max value       : ', maxval( (/( maxval(abs(offdiag(green_(kk_,i,:,:)))), i=1,n_frequ-nspec_frequ )/) )
         write (*, *) '============> chose frequency : ', uuu(1)
      elseif (rotation_scheme == 3) then
         mmat = green_(kk_, n_frequ - 2, :, :)
         call invmat(channels, mmat)
         mmat = -mmat ! this is now the projected hamiltonian
      else
         call utils_abort('Rotation scheme not implemented')
      endif

      if (channels > 1) then
         call write_array(real(mmat), 'Re Matrix to Diag', short=.true., unit=6)
         call write_array(aimag(mmat), 'Im Matrix to Diag', short=.true., unit=6)
         call eigenvector_matrix_c_(channels, mmat, diagdens(kk_, :), rotation(:, :, kk_), symmetric=.true.)
      else
         rotation(:, :, kk_) = 1; diagdens(kk_, 1) = mmat(1, 1)
      endif
      call write_array(rotation(:, :, kk_), ' COMPLEX ORTHOGONAL TRANSFORM ', short=.true., unit=6)

      if (rotation_scheme == 4) then
         rotation(:, :, kk_) = real(rotation(:, :, kk_))
         if (allocated(Qmat)) deallocate (Qmat)
         if (allocated(Rmat)) deallocate (Rmat)
         allocate (Qmat(channels, channels), Rmat(channels, channels))
         call QR_decomp(channels, channels, rotation(:, :, kk_), Qmat, Rmat)
         rotation(:, :, kk_) = Qmat
      endif

      if (channels > 1) call rearrange_columns_to_identity(channels, rotation(:, :, kk_), diagdens(kk_, :))

      write (*, *) ' ROTATION SCHEME : ', rotation_scheme
      write (*, *) ' PARAMAGNETIC ?  : ', rotation_scheme_pm
      write (*, *) ' FOR SPIN        : ', kk_

      if (rotation_scheme_pm .and. kk_ == 2) then
         rotation(:, :, 2) = rotation(:, :, 1)
         diagdens(2, :) = diagdens(1, :)
      endif

      if (channels > 1) then
   call write_array(MATMUL(MATMUL(transpose(rotation(:,:,kk_)),mmat),rotation(:,:,kk_)),'DIAG MATRIX RE-OBTAINED, AFTER COLUMN REORDERING ', short=.true., unit=6)
         call write_array(rotation(:, :, kk_), ' Q DECOMP ORTHOGONAL TRANSFORM ', short=.true., unit=6)
      endif

      write (*, *) ' DIAGONAL DENSITIES ARE : ', diagdens(kk_, :)
      ddiag = 0.
      do i = 1, channels
         ddiag(i, i) = diagdens(kk_, i)
      enddo

      if (rotation_scheme == 4 .and. rotation_scheme_read_write) then
         inquire (file='mask_dmft_rot'//trim(adjustl(toString(kk_))), exist=check)
         if (check) then
            open (unit=5000, file='mask_dmft_rot'//trim(adjustl(toString(kk_))), form='unformatted')
            read (5000) rotation(:, :, kk_)
            close (5000)
            call write_array(rotation(:, :, kk_), ' READ ORTHOGONAL TRANSFORM ', short=.true., unit=6)
         endif
         if (.not. check) then
            open (unit=5000, file='mask_dmft_rot'//trim(adjustl(toString(kk_))), form='unformatted')
            write (5000) rotation(:, :, kk_)
            close (5000)
         endif
      endif

      if (rotation_scheme == 5) then
         open (unit=5000, file='mask_user_rot')
         do i = 1, channels
            read (5000, *) (rotation(i, jjj, kk_), jjj=1, channels)
         enddo
         close (5000)
      endif

      write (*, *) '------------------------------------------------'
      call write_array(mmat, ' ONETEP DENSITY MATRIX    ', short=.true., unit=6)
      call write_array((rotation(:, :, kk_) + transpose(rotation(:, :, kk_)))/2., ' O + O^T ', ultrashort=.true., unit=6)
      call write_array((rotation(:, :, kk_) - transpose(rotation(:, :, kk_)))/2., ' O - O^T ', ultrashort=.true., unit=6)
      call write_array((rotation(:, :, kk_) + conjg(transpose(rotation(:, :, kk_))))/2., ' O + O^dag ', ultrashort=.true., unit=6)
      call write_array((rotation(:, :, kk_) - conjg(transpose(rotation(:, :, kk_))))/2., ' O - O^dag ', ultrashort=.true., unit=6)
  call write_array(matmul(rotation(:,:,kk_)-Id(channels),rotation(:,:,kk_)-Id(channels))  ,' ,Lambda^2', ultrashort=.true., unit=6)  
  call write_array(matmul(mmat,rotation(:,:,kk_)-Id(channels))-matmul(rotation(:,:,kk_)-Id(channels),mmat)  ,'[G,Lambda], O=1+Lambda',ultrashort=.true., unit=6)
      write (*, *) '------------------------------------------------'
      call write_array(matmul(rotation(:, :, kk_), transpose(rotation(:, :, kk_))), ' U * U^T       ', short=.true., unit=6)
      call write_array(matmul(rotation(:, :, kk_), transpose(conjg(rotation(:, :, kk_)))), ' U * U\dagger  ', short=.true., unit=6)
  call write_array(MATMUL(MATMUL(rotation(:,:,kk_),ddiag),transpose(rotation(:,:,kk_)))-mmat,' CHECK INVERSE (should be 0) ', short=.true., unit=6)
  call write_array(MATMUL(MATMUL(transpose(rotation(:,:,kk_)),mmat),rotation(:,:,kk_)),' should give diagonal densities ', short=.true., unit=6)

      uuu = maxloc((/(maxval(abs(    &
                     &  offdiag(matmul(matmul(transpose(rotation(:, :, kk_)), green_(kk_, i, :, :)), rotation(:, :, kk_)))  &
                     &  )), i=1, n_frequ - nspec_frequ)/))
      mmat = matmul(matmul(transpose(rotation(:, :, kk_)), green_(kk_, uuu(1), :, :)), rotation(:, :, kk_))
      call write_array(mmat, 'green of frequ for max offdiagonal elements after the rotation', short=.true., unit=6)

      uuu = minloc((/(maxval(abs(    &
                     &  offdiag(matmul(matmul(transpose(rotation(:, :, kk_)), green_(kk_, i, :, :)), rotation(:, :, kk_)))  &
                     &  )), i=1, n_frequ - nspec_frequ)/))
      mmat = matmul(matmul(transpose(rotation(:, :, kk_)), green_(kk_, uuu(1), :, :)), rotation(:, :, kk_))
      call write_array(mmat, 'green of frequ for min offdiagonal elements after the rotation', short=.true., unit=6)

      close (funit)

      if (kk_ == 1) then
         inquire (file='atom_number_of_orbitals', exist=check)
         if (check) then
            funit = utils_unit()
            open (unit=funit, file='atom_number_of_orbitals')
            read (funit, *) LL
            close (funit)
            LL = (LL - 1)/2
            if (LL /= (channels - 1)/2) then
               write (*, *) 'A PROJECTION WAS USED, TRUE ANGULAR MOMENTUM IS : ', LL
               write (*, *) 'AND [X] ORBITALS WERE KEPT                      : ', channels
               call utils_assert(solver == 4, 'Projections not implemented for solvers other than ED')
            endif
         else
            call utils_abort('./atom?/atom_number_of_orbitals file is missing')
         endif
      endif

      if (kk_ == 1) then
         allocate (list_orb(2*LL + 1))
         inquire (file='mask_projections', exist=check)
         call utils_assert(check, './atom?/mask_projection file is missing')
         funit = utils_unit()
         open (unit=funit, file='mask_projections')
         read (funit, *) (list_orb(i), i=1, 2*LL + 1)
         close (funit)
      endif

      return
   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine rotate_green_sigma(mat, kk_)
      implicit none
      complex(kind=DP) :: mat(:, :)
      integer    :: kk_
      if (.not. rotation_green_function) return
      mat = MATMUL(MATMUL(transpose(rotation(:, :, kk_)), mat), rotation(:, :, kk_))
   end subroutine

   !-----------------!

   subroutine rotate_back_sigma(mat, kk_)
      implicit none
      complex(kind=DP) :: mat(:, :)
      integer    :: kk_
      if (.not. rotation_green_function) return
      mat = MATMUL(MATMUL(rotation(:, :, kk_), mat), transpose(rotation(:, :, kk_)))
   end subroutine

   !-----------------!

   subroutine rotate_sigma_ortho(mat, kk_)
      implicit none
      complex(kind=DP) :: mat(:, :)
      integer    :: kk_
      mat = MATMUL(MATMUL(transpose(Simp_root_m(kk_, :, :)), mat), Simp_root_m(kk_, :, :))
   end subroutine

   subroutine rotate_green_ortho(mat, kk_)
      implicit none
      complex(kind=DP) :: mat(:, :)
      integer    :: kk_
      mat = MATMUL(MATMUL(transpose(Simp_root_gm(kk_, :, :)), mat), Simp_root_gm(kk_, :, :))
   end subroutine

   !-----------------!

   subroutine rotate_back_sigma_ortho(mat, kk_)
      implicit none
      complex(kind=DP) :: mat(:, :)
      integer    :: kk_
      mat = MATMUL(MATMUL(transpose(Simp_root_p(kk_, :, :)), mat), Simp_root_p(kk_, :, :))
   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine extract_Eimp_from_tail_of_hybridization(frequ, delta, my_eimp, a, b, left)

      use common_def, only: utils_abort

      implicit none
      complex(kind=DP)       :: delta(:), frequ(:)
      integer          :: ii, i, j, k, nn, i1, i2, i3, countit
      real(kind=DP)          :: my_eimp_cor, my_eimp, w1, w2, w3, a, b, c, d1i, d2i, d3i, test, d1r, d2r, d3r
      real(kind=DP)          :: distance(2), temp, a_(2), b_(2), e_(2), t1(2), t2(2), d
      logical, optional :: left

      countit = 0
      nn = size(frequ)

      !====================!
      !====================!
      !====================!
      !====================!

      if (real_or_imaginary_solver == 2) then

         i1 = nn - (nspec_frequ + 1)
         i2 = nn - nspec_frequ

82       continue
         countit = countit + 1
         i1 = i1 - 1
         i2 = i2 - 1

         if (i1 < 2 .or. i2 < 2) then
            a = 0; b = 0
            return
         endif
         w1 = aimag(frequ(i1))
         w2 = aimag(frequ(i2))

         if (tail_linear_scaling > 0.) then
            if (abs(w1) > abs(tail_linear_scaling)) then
#ifdef debug
               write (*, *) 'FREQUENCY TO FIND TAIL PARAMETERS : ', w1
#endif
               goto 82
            endif
         endif

         d1i = aimag(delta(i1))
         d2i = aimag(delta(i2))
         d1r = real(delta(i1))
         d2r = real(delta(i2))

         if (abs(w1) < 1.d-7 .or. abs(w2) < 1.d-7) then
            write (*, *) 'w(i) and w(j), for i,j = :', i1, i2
            write (*, *) 'have some trouble, total number of frequency is : ', n_frequ
            call utils_abort('w1 or w2 is zero, this should not happen for matsubara frequencies')
         endif
         if (abs(d1i) < 1.d-5 .or. abs(d2i) < 1.d-5) then
            if (maxval(abs(aimag(delta))) < 1.d-5) then
               a = 0.
               b = 0.
               my_eimp = d2r
               return
            endif
            goto 82
         endif
         if (abs(d1r) < 1.d-5 .or. abs(d2r) < 1.d-5) then
            write (*, *) 'd1r or d2r zero'
            write (*, *) 'd1r, d2r           : ', d1r, d2r
            write (*, *) 'd1i, d2i           : ', d1i, d2i
            write (*, *) 'w1,w2              : ', w1, w2
            a = -w1*d1i
            my_eimp = d1r
            b = 0.d0
            return
         endif
         if (abs(w1*d2i - w2*d1i) < 1.d-5) then
            b = 0.
            a = -w1*d1i
            my_eimp = d1r
            return
         endif
         d = w2/w1*d1i/d2i
         temp = d*w1**2 - w2**2
         temp = temp/(1.d0 - d)

         if (temp < -1.d-3) then
            if (mod(countit, 100) == 0) then
               write (*, *) 'error matching tail, negative sqrt'
               write (*, *) 'w1*d1i-w2*d2i  : ', w1*d1i - w2*d2i
               write (*, *) 'd1r,d2r        : ', d1r, d2r
               write (*, *) 'w1,w2          : ', w1, w2
               write (*, *) 'd1i,d2i        : ', d1i, d2i
               write (*, *) 'w1**2,w2**2    : ', w1**2, w2**2
               write (*, *) 'd              : ', d
               write (*, *) '1-d            : ', 1.d0 - d
               write (*, *) 'temp           : ', temp
            endif
            goto 82
         elseif (abs(temp) < 1.d-3) then
            temp = 0.d0
            my_eimp_cor = d1r
            a = -w1*d1i
            b = 0.d0
            b_(1) = 0.
            a_(1) = a
            e_(1) = my_eimp_cor
            distance(2) = 1.d10
            distance(1) = abs(delta(i1) - (a_(1)/(frequ(i1) + b_(1)) + e_(1)))
            goto 102
         endif

         do ii = 1, 2
         if (ii == 1) then
            b_(ii) = sqrt(temp)
         else
            b_(ii) = -sqrt(temp)
         endif
         t1(ii) = w1**2 + b_(ii)**2
         t2(ii) = w2**2 + b_(ii)**2
         a_(ii) = (d1r - d2r)/b_(ii)*t1(ii)*t2(ii)/(t2(ii) - t1(ii))
         a_(ii) = -d1i*t1(ii)/w1
         e_(ii) = d1r - a_(ii)*b_(ii)/t1(ii)
         distance(ii) = abs(delta(i1) - (a_(ii)/(frequ(i1) + b_(ii)) + e_(ii)))
         enddo

102      continue

         if (distance(1) < 1.d-7) then
            a = a_(1); b = b_(1); my_eimp = e_(1)
         elseif (distance(2) < 1.d-7) then
            a = a_(2); b = b_(2); my_eimp = e_(2)
         else
            write (*, *) 'ERROR : error couldnt match tail, distances : ', distance
            write (*, *) 'frequ 1 : ', frequ(1)
            do ii = 1, 2
               write (*, *) 'ii      : ', ii
               write (*, *) ' a,b,e  : ', a_(ii), b_(ii), e_(ii)
               write (*, *) 'im      : ', aimag(delta(i1) - (a_(ii)/(frequ(i1) + b_(ii)) + e_(ii)))
               write (*, *) 're      : ', real(delta(i1) - (a_(ii)/(frequ(i1) + b_(ii)) + e_(ii)))
            enddo
         endif

         !====================!
         !====================!
         !====================!
         !====================!
      else
         !====================!
         !====================!
         !====================!
         !====================!

         if (.not. present(left)) then
            i1 = nn - (nspec_frequ + 2)
            i2 = nn - (nspec_frequ + 1)
            i3 = nn - (nspec_frequ)
         else
            i1 = 2
            i2 = 3
            i3 = 4
         endif

91       continue

         if (.not. present(left)) then
            i1 = i1 - 1
            i2 = i2 - 1
            i3 = i3 - 1
         else
            i1 = i1 + 1
            i2 = i2 + 1
            i3 = i3 + 1
         endif

         if (i1 < 2) then
            a = 0; b = 0
            return
         endif

         w1 = real(frequ(i1))
         w2 = real(frequ(i2))
         w3 = real(frequ(i3))
         d1r = real(delta(i1))
         d2r = real(delta(i2))
         d3r = real(delta(i3))

         if (tail_linear_scaling > 0.) then
            if (abs(w1) > abs(tail_linear_scaling)) then
#ifdef debug
               write (*, *) 'FREQUENCY TO FIND TAIL PARAMETERS : ', w1
#endif
               goto 91
            endif
         endif

         if (abs(w1) < 1.d-5 .or. abs(w2) < 1.d-5 .or. abs(w3) < 1.d-5) then
            write (*, *) 'w(i) and w(j), for i,j = :', i1, i2
            write (*, *) 'have some trouble, total number of frequency is : ', n_frequ
            write (*, *) 'ERROR : w1 or w2 or w3 zero'
            goto 91
         endif

         if (abs(d1r) < 1.d-5 .or. abs(d2r) < 1.d-5 .or. abs(d3r) < 1.d-5) then
            if (maxval(abs(real(delta))) < 1.d-5) then
               a = 0.
               b = 0.
               return
            endif
            goto 91
         endif

         c = (d2r - d3r)/(d1r - d2r)
         b = (w3*c - w1)/(1.d0 - c)
         a = (w1 + b)*(w2 + b)*(d1r - d2r)/(w2 - w1)
         my_eimp_cor = d1r - a/(w1 + b)
         test = d2r - a/(w2 + b)

         if (abs(my_eimp_cor - test) > 1.d-5) then
            write (*, *) 'ERROR : extract Eimp : test failed'
            goto 91
         endif

         distance(1) = abs(delta(i1) - (a/(frequ(i1) + b) + my_eimp_cor))
         if (distance(1) > 1.d-3) then
            write (*, *) 'ERROR : real axis calculation, error couldnt match tail, distances : ', distance(1)
            write (*, *) 'frequ 1 : ', frequ(1)
            write (*, *) 'Left?   : ', present(left)
            write (*, *) ' a,b,e  : ', a, b, my_eimp_cor
            write (*, *) 'im      : ', aimag(delta(i1) - (a/(frequ(i1) + b) + my_eimp_cor))
            write (*, *) 're      : ', real(delta(i1) - (a/(frequ(i1) + b) + my_eimp_cor))
            goto 91
         endif

         my_eimp = my_eimp_cor

         !====================!
         !====================!
         !====================!
         !====================!
      endif

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   function extended_delta_real_frequ(kk_, ii)
      implicit none
      real(kind=DP)    :: extended_delta_real_frequ
      integer    :: i, j, k, l, ii, kk_
      real(kind=DP)    :: w0, wf, wstep
      w0 = myfrequ(kk_, 1)
      wf = myfrequ(kk_, n_frequ)
      wstep = (wf - w0)/dble(n_frequ - 1)
      extended_delta_real_frequ = w0 - (nmatsu_long - n_frequ)/2*wstep + dble(ii - 1)*wstep
   end function

   !====================!
   !====================!

   function extended_delta(kk_, ii, ah, bh)
      implicit none
      real(kind=DP)    :: extended_delta(2*channels + 1), ah(channels), bh(channels)
      integer    :: i, j, k, l, ii, kk_
      real(kind=DP)    :: w0, wf, wstep
      complex(kind=DP) :: ww, tmp(channels)

      if (real_or_imaginary_solver == 1) then
         w0 = myfrequ(kk_, 1)
         wf = myfrequ(kk_, n_frequ)
         wstep = (wf - w0)/dble(n_frequ - 1)
         extended_delta(1) = w0 - (nmatsu_long - n_frequ)/2*wstep + dble(ii - 1)*wstep
         ww = extended_delta(1)
      else
         extended_delta(1) = dacos(-1.d0)/bbeta*dble(2*ii - 1)   ! pi/beta * (2n-1)
         ww = imi*extended_delta(1)
      endif

      do i = 1, channels
         if (real_or_imaginary_solver == 1) then
            tmp(i) = ah(i)/(ww + bh(i))
         else
            tmp(i) = ah(i)/((ww + aimag(frequ_(kk_, 1))) + bh(i))
         endif
      enddo

      do i = 1, channels
         extended_delta(1 + 2*i - 1) = real(tmp(i))
         extended_delta(1 + 2*i) = aimag(tmp(i))
      enddo

   end function

   !====================!
   !====================!

   function extended_delta_el(kk_, ii, ah, bh)
      implicit none
      complex(kind=DP) :: extended_delta_el
      real(kind=DP)    :: ah, bh
      integer    :: i, j, k, l, ii, kk_
      real(kind=DP)    :: w0, wf, wstep, wwext
      complex(kind=DP) :: ww, tmp

      if (real_or_imaginary_solver == 1) then
         w0 = myfrequ(kk_, 1)
         wf = myfrequ(kk_, n_frequ)
         wstep = (wf - w0)/dble(n_frequ - 1)
         wwext = w0 - (nmatsu_long - n_frequ)/2*wstep + dble(ii - 1)*wstep
         ww = wwext
      else
         wwext = dacos(-1.d0)/bbeta*dble(2*ii - 1)   ! pi/beta * (2n-1)
         ww = imi*wwext
      endif

      if (real_or_imaginary_solver == 1) then
         tmp = ah/(ww + bh)
      else
         tmp = ah/((ww + aimag(frequ_(kk_, 1))) + bh)
      endif

      extended_delta_el = tmp

   end function

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine copy_sigma_position_nlong_to_position_nw_minus_1(filenamesig)

      use common_def, only: utils_unit, utils_abort

      implicit none
      character*(*)  ::  filenamesig
      real(kind=DP)        ::  temp(nmatsu_long, 2*channels + 1)
      integer        :: funit

      funit = utils_unit()
      open (unit=funit, file=trim(adjustl(filenamesig)))
      if (solver == 1) then
         !header in self energy of CTQMC solver
         read (funit, *)
      endif
      do i = 1, nmatsu_long
         read (funit, *, end=88) (temp(i, k), k=1, 1 + 2*channels)
      enddo
      if (.false.) then
88       write (*, *) 'ERROR : end of sigma file, it should have been larger due to nmatsu_long extension of bath'
         call utils_abort('sigma file ended prematurely')
      endif
      close (funit)

      funit = utils_unit()
      call utils_system_call(" rm "//trim(adjustl(filenamesig)))
      open (unit=funit, file=trim(adjustl(filenamesig)))
      write (funit, *)
      do i = 1, n_frequ
         if (i /= n_frequ - 1 .and. i /= n_frequ - 5) then
            write (funit, *) (temp(i, k), k=1, 1 + 2*channels)
         else
            write (funit, *) (temp(nmatsu_long, k), k=1, 1 + 2*channels)
            write (*, *) 'CORRECTION IMPROVED BY MATSU_LONG : '
            write (*, '(a,200f10.3)') 'sigma at long   : ', temp(nmatsu_long, 2:)
            write (*, '(a,200f10.3)') 'sigma at nfrequ : ', temp(i, 2:)
         endif
      enddo
      close (funit)

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   integer function my_zero_frequ(kk_)
      implicit none
      integer :: ii, kk_
      if (real_or_imaginary_solver == 2) then
         my_zero_frequ = 1
      else
         my_zero_frequ = minloci(abs((/(real(frequ_(kk_, ii)) - 0.d0, ii=1, n_frequ - nspec_frequ)/)))
      endif
   end function

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   real(kind=DP) function myfrequ(kk_, ii)
      implicit none
      integer :: ii, kk_
      if (real_or_imaginary_solver == 2) then
         myfrequ = aimag(frequ_(kk_, ii))
      else
         myfrequ = real(frequ_(kk_, ii))
      endif
   end function

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   function Rotate_jj_spin_orbit_for_f()
      implicit none
      real(kind=DP)    :: tt(14, 28), ttt(14, 28)
      complex(kind=DP) :: tj(14, 14), tc(14, 14), Rc(14, 14), Rotate_jj_spin_orbit_for_f(14, 14)
      integer    :: i, k, l, m, ispin

       tt(1, 1:28) =  [ 0.92582010d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,  -0.37796447d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0]
       tt(2, 1:28) =  [ 0d0, 0d0,   0.84515425d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,  -0.53452248d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0]
       tt(3, 1:28) =  [ 0d0, 0d0,   0d0, 0d0,   0.75592895d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,  -0.65465367d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0]
       tt(4, 1:28) =  [ 0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0.65465367d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,  -0.75592895d0, 0d0,   0d0, 0d0,   0d0, 0d0]
       tt(5, 1:28) =  [ 0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0.53452248d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,  -0.84515425d0, 0d0,   0d0, 0d0]
       tt(6, 1:28) =  [ 0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0.37796447d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,  -0.92582010d0, 0d0]
       tt(7, 1:28) =  [ 0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   1.00000000d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0.d0        ]
       tt(8, 1:28) =  [ 0.37796447d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0.92582010d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0]
       tt(9, 1:28) =  [ 0d0, 0d0,   0.53452248d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0.84515425d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0]
       tt(10,1:28) =  [ 0d0, 0d0,   0d0, 0d0,   0.65465367d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0.75592895d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0]
       tt(11,1:28) =  [ 0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0.75592895d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0.65465367d0, 0d0,   0d0, 0d0,   0d0, 0d0]
       tt(12,1:28) =  [ 0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0.84515425d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0.53452248d0, 0d0,   0d0, 0d0]
       tt(13,1:28) =  [ 0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0.92582010d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0.37796447d0, 0d0]
       tt(14,1:28) =  [ 0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   1.00000000d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0,   0d0, 0d0]

      do i = 1, 14
         do j = 1, 14
            tj(i, j) = CMPLX(tt(i, 2*j - 1), tt(i, 2*j))
         enddo
      enddo
      tj = transpose(tj)

      !tj : spherical to j-j, spin up+dn

       ttt(1 ,1:28) = [0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.70710678d0,  0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 ]
       ttt(2 ,1:28) = [0.00d0, 0.00d0, 0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, -0.70710678d0, 0.00d0, 0.00d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 ]
       ttt(3 ,1:28) = [0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,  0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 ]
       ttt(4 ,1:28) = [0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 1.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,              0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 ]
       ttt(5 ,1:28) = [0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0, -0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 ]
       ttt(6 ,1:28) = [0.00d0, 0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0,  0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 ]
       ttt(7 ,1:28) = [0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, -0.70710678d0, 0.00d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 ]
       ttt(8 ,1:28) = [0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.70710678d0  ]
       ttt(9 ,1:28) = [0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0.00d0, 0.00d0, 0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, -0.70710678d0, 0.00d0, 0.00d0 ]
       ttt(10,1:28) = [0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0  ]
       ttt(11,1:28) = [0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 1.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0              ]
       ttt(12,1:28) = [0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0, -0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0 ]
       ttt(13,1:28) = [0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0.00d0, 0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0  ]
       ttt(14,1:28) = [0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0.70710678d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, -0.70710678d0, 0.00d0 ]

      do i = 1, 14
         do j = 1, 14
            tc(i, j) = CMPLX(ttt(i, 2*j - 1), ttt(i, 2*j))
         enddo
      enddo
      tc = transpose(tc)
      call invmat(14, tc)

      !tc : f orbitals to spherical harmonics, spin up+dn

      Rc = MATMUL(tc, tj)

      Rotate_jj_spin_orbit_for_f = Rc

      return
   end function

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine rotate_sigma_green_to_f(matup, matdn)

      use common_def, only: utils_abort

      implicit none
      complex(kind=DP) :: matup(:, :), matdn(:, :), mat(size(matup, 1) + size(matdn, 1), size(matup, 2) + size(matdn, 2))
      complex(kind=DP) :: rotmat(size(matup, 1) + size(matdn, 1), size(matup, 2) + size(matdn, 2))
      integer    :: kk_, i
      if (.not. rotation_green_function) return
      i = size(matup, 1)
      mat = 0.d0
      mat(1:i, 1:i) = matup
      mat(i + 1:2*i, i + 1:2*i) = matdn
      rotmat = Rotate_jj_spin_orbit_for_f()
      mat = MATMUL(MATMUL(transpose(conjg(rotmat)), mat), rotmat)
      matup = mat(1:i, 1:i)
      matdn = mat(i + 1:2*i, i + 1:2*i)
      if (maxval(abs(rotmat - transpose(conjg(rotmat)))) > 1.d-4) then
         write (*, *) 'TRANSFORM FROM F ORBITALS TO j-j basis not unitary'
         write (*, *) 'U U^dagger - 1 = ', maxval(abs(rotmat - transpose(conjg(rotmat))))
         call utils_abort('Transform from f orbitals to j-j basis not unitary')
      endif
   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine rotate_sigma_green_to_f_back(matup, matdn)

      use common_def, only: utils_abort

      implicit none
      complex(kind=DP) :: matup(:, :), matdn(:, :), mat(size(matup, 1) + size(matdn, 1), size(matup, 2) + size(matdn, 2))
      complex(kind=DP) :: rotmat(size(matup, 1) + size(matdn, 1), size(matup, 2) + size(matdn, 2))
      integer    :: kk_, i
      if (.not. rotation_green_function) return
      i = size(matup, 1)
      mat = 0.d0
      mat(1:i, 1:i) = matup
      mat(i + 1:2*i, i + 1:2*i) = matdn
      rotmat = Rotate_jj_spin_orbit_for_f()
      mat = MATMUL(MATMUL(rotmat, mat), transpose(conjg(rotmat)))
      matup = mat(1:i, 1:i)
      matdn = mat(i + 1:2*i, i + 1:2*i)
      if (maxval(abs(rotmat - transpose(conjg(rotmat)))) > 1.d-4) then
         write (*, *) 'TRANSFORM FROM F ORBITALS TO j-j basis not unitary'
         write (*, *) 'U U^dagger - 1 = ', maxval(abs(rotmat - transpose(conjg(rotmat))))
         call utils_abort('Transform from f orbitals to j-j basis not unitary')
      endif
   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine rotations_and_interactions_and_double_counting

      use common_def, only: utils_assert, utils_abort

      implicit none
      integer    :: j_, i, j, k, kk_, iii, jjj
      real(kind=DP)    :: zav, dummy(channels, channels)
      complex(kind=DP) :: temp_mat(channels, channels)
      complex(kind=DP) :: Smatloc(2, channels, channels), COEFloc(2, channels)
      logical    :: test

      if (paramagnetic == 1) then
         j_ = 1
      else
         j_ = 2
      endif

      do kk_ = 1, j_
         !below dummy is simp^-1, not simp
         call getting_S_matrix_from_Green(Smatloc(kk_, :, :), COEFloc(kk_, :), channels, kk_, silent=.true., inverse=.true.)
      enddo

      do kk_ = 1, j_
         !WE NEED THOSE BEFORE THE ROTATIONS - SINCE THIS TRANSFORMATION IS APPLIED
         !AT THE VERY END, TO THE FINAL SIGMA (IT IS CONNECTED WITH THE UPFOLDING OF SIGMA)
         if (use_simp_from_onetep_for_ortho) then
            dummy = s_mat_from_onetep(kk_, channels)
         else
            dummy = Smatloc(kk_, :, :)
         endif
         call write_array(dummy, 'Simp (unrotated) = U^+ S^-1 U ', short=.true., unit=6)
         Simp_m1(kk_, :, :) = dummy
      enddo

      if (rotation_green_function) then
         if (.not. flag_use_jj_basis_for_f) then
            do kk_ = 1, j_
               do i = 1, n_frequ
                  call rotate_green_sigma(sigma_mat(kk_, i, :, :), kk_)
                  call rotate_green_sigma(green_(kk_, i, :, :), kk_)
                  green_diag(kk_, i, :) = diag(green_(kk_, i, :, :))
               enddo
            enddo
         else
            if (LL /= 3) then
               write (*, *) 'ERROR, flag flag_use_jj_basis_for_f only for F orbitals (L=3)'
               write (*, *) 'here L=', LL
               call utils_abort("flag_use_jj_basis_for_f is only valid for f-shell orbitals")
            endif
            call utils_assert(paramagnetic /= 1, 'Spin orbit is not valid for paramagnetic case')
            do i = 1, n_frequ
               call rotate_sigma_green_to_f(sigma_mat(1, i, :, :), sigma_mat(2, i, :, :))
               call rotate_sigma_green_to_f(green_(1, i, :, :), green_(2, i, :, :))
               do kk_ = 1, 2
                  green_diag(kk_, i, :) = diag(green_(kk_, i, :, :))
               enddo
            enddo
         endif
      endif

! What follows is a simple correction : "Sigma" from file below ...
!  ... is actually Sigma-Edc

      write (*, *) 'double counting correction'

      do kk_ = 1, j_
         do i = 1, n_frequ
            if (.not. double_counting_zero_self .and. maxval(abs(double_counting(kk_, :))) < 1.d-5 .and. abs(UU) > 1.d-5) then
               write (*, '(a)') ' WARNING: double counting is undefined/zero'
               ! write (*, *) ' UU,double counting                 : ', UU, double_counting(kk_, :)
               ! call utils_abort('Double counting is undefined')
            endif
            if (maxval(abs(sigma_mat(kk_, n_frequ, :, :))) > 1.d-5) then
               do k = 1, channels
                  ccctemp = sigma_mat(kk_, n_frequ, k, k)
                  if (i /= n_frequ) sigma_mat(kk_, i, k, k) = sigma_mat(kk_, i, k, k) + ccctemp
               enddo
            else
               do k = 1, channels
                  sigma_mat(kk_, i, k, k) = sigma_mat(kk_, i, k, k) + double_counting(kk_, k)
               enddo
            endif
         enddo
      enddo

      do kk_ = 1, j_
         !below dummy is simp^-1, not simp
         call getting_S_matrix_from_Green(Smatloc(kk_, :, :), COEFloc(kk_, :), channels, kk_, silent=.true., inverse=.true.) !below simp^-1 is used, not simp
      enddo

      UU_ren = UU_ren0
      JJ_ren = JJ_ren0

      Zimp_ren_m = 1.d0

      if (rotate_to_ortho_basis) then
         write (*, *) 'ROTATING BACK TO ORTHONORMAL AXIS'
         if (.not. flag_use_jj_basis_for_f) then
            do kk_ = 1, j_
               if (use_simp_from_onetep_for_ortho) then
                  Simp_back(kk_, :, :) = s_mat_from_onetep(kk_, channels)
               else
                  Simp_back(kk_, :, :) = Smatloc(kk_, :, :)
               endif
               call invmat(channels, Simp_back(kk_, :, :))
               call matrix_square_root(Simp_back(kk_, :, :), Simp_root_gp(kk_, :, :), Simp_root_gm(kk_, :, :), &
                                                        & Simp_root_p(kk_, :, :), Simp_root_m(kk_, :, :),  &
                                                        & Zimp_ren_p(kk_, :, :), Zimp_ren_m(kk_, :, :), Simp_rot(kk_, :, :))
    call write_array(matmul(matmul(transpose(Simp_root_m(kk_,:,:)),Simp_back(kk_,:,:)),Simp_root_m(kk_,:,:)),'S_m^T Simp S_m = Id',unit=6,short=.true.)

               if (use_simp_from_onetep_for_ortho) then
            temp_mat = matmul(matmul(transpose(Simp_root_gm(kk_, :, :)), s_mat_from_onetep(kk_, channels)), Simp_root_gm(kk_, :, :))
               else
                  temp_mat = matmul(matmul(transpose(Simp_root_gm(kk_, :, :)), Smatloc(kk_, :, :)), Simp_root_gm(kk_, :, :))
               endif

               call invmat(channels, temp_mat)
               call write_array(temp_mat, '( (Sp^+ G Sp)^-1', unit=6, short=.true.)
               !---------------!
               if (kk_ == 1) then
                  do iii = 1, channels
                     do jjj = 1, channels
                        UU_ren(iii, jjj) = UU_ren0(iii, jjj)*(Zimp_ren_m(kk_, iii, iii)**2)*(Zimp_ren_m(kk_, jjj, jjj)**2) !--> this goes to ED-Kanamori
                        JJ_ren(iii, jjj) = JJ_ren0(iii, jjj)*(Zimp_ren_m(kk_, iii, iii)**2)*(Zimp_ren_m(kk_, jjj, jjj)**2) !--> this goes to ED-Kanamori
                     enddo
                  enddo
                  zav = sum(diag(Zimp_ren_m(kk_, :, :))**2)/dble(channels)
                  if (rotate_ortho_av_renorm_int) then
                     if (LL == 3) then
                        write (*, *) 'SORRY RENORMALIZATION NOT YET CODED FOR L=3'
                        write (*, *) 'the problem is the atomic solver which needs to be modified'
                        write (*, *) 'in particular for L=1,2 atom_d.py was modified to accept home made Uc.dat file'
                        call utils_abort('Renormalization not yet implemented for L = 3')
                     endif
                     UU = UU0/zav/zav !--> this goes to CTQMC,NCA,OCA
                     Jhund = Jhund0/zav/zav !--> this goes to CTQMC,NCA,OCA
                     UU_ren = UU_ren0/zav/zav !--> this goes to ED-Kanamori
                     JJ_ren = JJ_ren0/zav/zav !--> this goes to ED-Kanamori
                  endif
               endif
               !---------------!

               !h_imp and s_proj do not transform like G, this is not a unitary
               !transformation
         green_(kk_, n_frequ - 4, :, :) = matmul(matmul(Simp_back(kk_, :, :), green_(kk_, n_frequ - 4, :, :)), Simp_back(kk_, :, :))
         green_(kk_, n_frequ - 5, :, :) = matmul(matmul(Simp_back(kk_, :, :), green_(kk_, n_frequ - 5, :, :)), Simp_back(kk_, :, :))
               call rotate_sigma_ortho(green_(kk_, n_frequ - 4, :, :), kk_)
               call rotate_sigma_ortho(green_(kk_, n_frequ - 5, :, :), kk_)

               do i = 1, n_frequ
                  if (i /= n_frequ - 4 .and. i /= n_frequ - 5) then
                     call rotate_sigma_ortho(sigma_mat(kk_, i, :, :), kk_)
                     call rotate_green_ortho(green_(kk_, i, :, :), kk_)
                     green_diag(kk_, i, :) = diag(green_(kk_, i, :, :))
                  endif
               enddo
            enddo
         else
            call utils_abort('Rotation back to orthogonal basis is not yet implemented for f-shell orbitals')
         endif
      endif

      write (*, *) 'RENORMALIZATION ZIMP SPIN 1 : ', diag(Zimp_ren_m(1, :, :))
      if (.not. paramagnetic == 1) then
         write (*, *) 'RENORMALIZATION ZIMP SPIN 2 : ', diag(Zimp_ren_m(2, :, :))
         if (maxval(abs(Zimp_ren_m(1, :, :) - Zimp_ren_m(2, :, :))) > 3.d-3) then
            write (*, *) 'ERROR, renormalization different for up and dn spin, not sure'
            write (*, *) '       how the code fares with that....critical...'
            write (*, *) '       -----------------> '
            write (*, *) '       we symmetrize the renormalization          '
            Zimp_ren_m(1, :, :) = (Zimp_ren_m(1, :, :) + Zimp_ren_m(2, :, :))/2.d0
            Zimp_ren_m(2, :, :) = Zimp_ren_m(1, :, :)
            !stop
         endif
      endif

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine matrix_square_root(matin, root_gp, root_gm, root_p, root_m, ZZ_p, ZZ_m, rot)
      real(kind=DP)  ::  matin(:, :)
      real(kind=DP)  ::  root_p(size(matin, 1), size(matin, 2)), vaps(size(matin, 1))
      real(kind=DP)  ::  root_m(size(matin, 1), size(matin, 2))
      real(kind=DP)  ::  root_gp(size(matin, 1), size(matin, 2))
      real(kind=DP)  ::  root_gm(size(matin, 1), size(matin, 2))
      real(kind=DP)  ::  tmp_mat_offdiag(size(matin, 1), size(matin, 2))
      real(kind=DP)  ::  vec(size(matin, 1), size(matin, 2)), temp(size(matin, 1), size(matin, 2))
      real(kind=DP)  ::  WORK(3*size(matin, 1)), ZZ_p(size(matin, 1), size(matin, 2)), ZZ_m(size(matin, 1), size(matin, 2))
      real(kind=DP)  ::  rot(size(matin, 1), size(matin, 2))
      integer  ::  i, j, k, l, nn, ierr

      nn = size(matin, 1)
      vec = matin

      if (maxval(abs(vec - Id(nn))) < 1.d-3) then
         vaps = 1.d0
         vec = Id(nn)
         goto 30
      endif

      tmp_mat_offdiag = matin
      do i = 1, size(matin, 1)
         tmp_mat_offdiag(i, i) = 0.d0
      enddo
      if (maxval(abs(tmp_mat_offdiag)) < cutoff_simp_offdiag) then
         vaps = diag(matin)
         vec = Id(nn)
         goto 30
      endif

      if (.not. assume_proj_overlap_is_diagonal) then
         call DSYEV('V', 'U', nn, vec, nn, vaps, WORK, 3*nn, ierr)
         call rearrange_columns_to_identity(nn, vec, vaps)
         do i = 1, nn
            if (vec(i, i) < -1.d-2) vec(:, i) = -vec(:, i)
         enddo
      else
         vaps = diag(matin)
         vec = Id(nn)
      endif

30    continue

      call write_array(matin, 'Simp : original matrix', unit=6, short=.true.)
      call write_array(vec, 'Simp : unitary  matrix', unit=6, short=.true.)

      rot = vec

      !sigma
      temp = 0.d0; do i = 1, nn; temp(i, i) = 1.d0/sqrt(abs(vaps(i))); enddo; 
      root_m = matmul(vec, temp)
      ZZ_m = temp
      !green
      temp = 0.d0; do i = 1, nn; temp(i, i) = sqrt(abs(vaps(i))); enddo; 
      root_gm = matmul(vec, temp)
      !sigma back
      temp = 0.d0; do i = 1, nn; temp(i, i) = sqrt(abs(vaps(i))); enddo; 
      root_p = matmul(temp, transpose(vec))
      ZZ_p = temp
      !green back
      temp = 0.d0; do i = 1, nn; temp(i, i) = 1.d0/sqrt(abs(vaps(i))); enddo; 
      root_gp = matmul(temp, transpose(vec))

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine define_double_counting

      use common_def, only: utils_abort

      implicit none
      integer                :: klm
      integer                :: i_, ii1, ii2
      real(kind=DP)                :: tt1_, tt2_(2), occup, Utmp, Jtmp, occup1, occup2
      logical                :: test

      if (dmft_for_dimer) then
         if (channels == 2) then
            Jhund = 0.d0
            write (*, *) 'WARNING : SWITCHING OFF HUNDs coupling, only 2 orbitals and dimers'
         endif
      else
         if (channels == 1) then
            Jhund = 0.d0
            write (*, *) 'WARNING : SWITCHING OFF HUNDs coupling, only 1 orbital'
         endif
      endif

      inquire (file='mask_u_matrix', exist=checkujmat)
      if (checkujmat) then
         open (unit=19881, file='mask_u_matrix')
         do k = 1, channels
            read (19881, *) (UU_ren0(k, i), i=1, channels)
         enddo
         close (19881)
      endif
      inquire (file='mask_j_matrix', exist=test)
      if (checkujmat .and. .not. test) then
         call utils_abort('mask_u_matrix present but not mask_j_matrix. Please provide both files.')
      endif
      if (checkujmat) then
         open (unit=19881, file='mask_j_matrix')
         do k = 1, channels
            read (19881, *) (JJ_ren0(k, i), i=1, channels)
         enddo
         close (19881)
      endif
      if (.not. checkujmat) then
         JJ_ren0 = Jhund
         UU_ren0 = UU
      endif
      Jhund0 = Jhund
      UU0 = UU

      occup = 0.d0

      if ((solver == 4 .or. solver == 5) .and. .not. double_counting_with_no_average) then
         if (dmft_for_dimer) then
            klm = channels/2
            Uaverage = (UU + 2.d0*dble(klm - 1)*(UU - 2.d0*Jhund))/dble(2*klm - 1)
         else
            Uaverage = (UU + 2.d0*dble(channels - 1)*(UU - 2.d0*Jhund))/dble(2*channels - 1)
         endif
      else
         Uaverage = UU
      endif

      if (double_counting_nf < -0.0001d0) then
         if (dmft_for_dimer) then
            klm = channels/2
            tt2_(1) = sum(diag(occupation_matrixup(1:klm, 1:klm)))
            tt2_(2) = sum(diag(occupation_matrixdn(1:klm, 1:klm)))
            if (dimer_average_occupation) then
               tt2_(1) = (tt2_(1) + sum(diag(occupation_matrixup(klm + 1:2*klm, klm + 1:2*klm))))/2.d0
               tt2_(2) = (tt2_(2) + sum(diag(occupation_matrixdn(klm + 1:2*klm, klm + 1:2*klm))))/2.d0
            endif
            tt1_ = sum(tt2_)
            occup1 = tt1_
            do i_ = 1, 2
               do ii1 = 1, klm
                  if (checkujmat) then; Utmp = UU_ren0(ii1, ii1); Jtmp = JJ_ren0(ii1, ii1); else; Utmp = Uaverage; Jtmp = Jhund; endif; 
                  double_counting(i_, ii1) = Utmp*(tt1_ - 0.5d0) - Jtmp*(tt2_(i_) - 0.5d0)
               enddo
            enddo
            tt2_(1) = sum(diag(occupation_matrixup(klm + 1:2*klm, klm + 1:2*klm)))
            tt2_(2) = sum(diag(occupation_matrixdn(klm + 1:2*klm, klm + 1:2*klm)))
            if (dimer_average_occupation) then
               tt2_(1) = (tt2_(1) + sum(diag(occupation_matrixup(1:klm, 1:klm))))/2.d0
               tt2_(2) = (tt2_(2) + sum(diag(occupation_matrixdn(1:klm, 1:klm))))/2.d0
            endif
            tt1_ = sum(tt2_)
            occup2 = tt1_
            do i_ = 1, 2
               do ii1 = 1, klm
                  if (checkujmat) then; Utmp = UU_ren0(klm + ii1, klm + ii1); Jtmp = JJ_ren0(klm + ii1, klm + ii1); else; Utmp = Uaverage; Jtmp = Jhund; endif; 
                  double_counting(i_, klm + ii1) = Utmp*(tt1_ - 0.5d0) - Jtmp*(tt2_(i_) - 0.5d0)
               enddo
            enddo
         else
            occup = occup + occupation_orbital
            do i_ = 1, channels
               if (checkujmat) then
                  Utmp = UU_ren0(i_, i_)
                  Jtmp = JJ_ren0(i_, i_)
               else
                  Utmp = Uaverage
                  Jtmp = Jhund
               endif 
               double_counting(:, i_) = Utmp*(occupation_orbital - 0.5d0) - Jtmp*(densupdn(:) - 0.5d0)
            enddo
         endif
      elseif (double_counting_nf > 0.0001) then
         occup = occup + double_counting_nf
         double_counting = Uaverage*(double_counting_nf - 0.5d0) - 0.5*Jhund*(double_counting_nf - 1.d0)
      elseif (abs(double_counting_nf) < 0.0001) then
         write (*, *) 'WARNING : NO DOUBLE COUNTING SCHEME WAS USED (double_counting_nf=0)'
         double_counting = 0.d0
         occup = 0.d0
      endif
      !correction to E_DC
      open(unit=1981,file=trim(adjustl(filename_edc)))

      if(dmft_for_dimer)then
         if(double_counting_zero_self)then
            write(1981,*)   0.d0, 0.d0
         else
            write(1981,*) (-Uaverage/4.d0 + Jhund/8.d0)*occup1, (-Uaverage/4.d0 + Jhund/8.d0)*occup2
         end if
      else
         if (double_counting_zero_self) then
            write (1981, *) 0.d0
         else
            write (1981, *) (-Uaverage/4.d0 + Jhund/8.d0)*occup
         end if
      end if

      write (1981, *) (-Uaverage/4.d0 + Jhund/8.d0)
      write (1981, *) '# Uav J n_LDA '
      write (1981, *) Uaverage
      write (1981, *) Jhund
      write (1981, *) occup

      close (1981)

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   function matmul_ex(A_trans, rot_dmft, extend)

      use common_def, only: utils_assert

      real(kind=DP)    :: A_trans(:, :), rot_dmft(:, :)
      logical    :: extend, check
      real(kind=DP)    :: matmul_ex(size(A_trans, 1), size(A_trans, 2))
      real(kind=DP)    :: temp(size(A_trans, 1), size(A_trans, 2))
      integer    :: i, j, k, l, i1, i2, j1, j2

      i1 = size(A_trans, 1); i2 = size(A_trans, 2)
      j1 = size(rot_dmft, 1); j2 = size(rot_dmft, 2)
      if (.not. extend) then
         call utils_assert(i1 == j1 .and. i2 == j2, 'Array size mismatch in matmul_ex')
         matmul_ex = matmul(A_trans, rot_dmft)
         return
      endif

      !  Trans.dat in onetep ngwfs basis.....
      !  rot_dmft in subspace of selected orbitals (mask_projections)
      !  we extend the dmft subspace by adding 1 in the matrix

      temp = Id(2*LL + 1)

      k = 0
      do i = 1, 2*LL + 1
         if (list_orb(i)) then
            k = k + 1
            l = 0
            do j = 1, 2*LL + 1
               if (list_orb(j)) then
                  l = l + 1
                  temp(i, j) = rot_dmft(k, l)
               endif
            enddo
         endif
      enddo
      matmul_ex = matmul(A_trans, temp)

   end function

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine change_phase_T2C(T2Cp)

      use common_def, only: utils_assert

      implicit none
      real(kind=DP)    :: ctg, Qp, Rm, Rp, Qm, xb, xa, phi
      integer    :: i, m
      complex(kind=DP) :: T2Cp(:, :), T2C(size(T2Cp, 2), size(T2Cp, 1))

      write (*, *) 'Updating phases of T2C'

      T2C = transpose(T2Cp)

      !   Imag( T2C[m,i] + (-1)**m * T2C[-m,i] ) == 0
      !   Real( T2C[m,i] - (-1)**m * T2C[-m,i] ) ==0
      !  \sum_m T2C[m,i]*exp(i*m*phi) is real for any phi
      !   T2C[m,i] -> T2C[m,i]*exp(i*phi_i)
      !   This leads to the following 2x2 system of equations:
      !    ( Rp[m,i], Qp[m,i] ) ( sin(phi_i) )
      !    ( Qm[m,i], Rm[m,i] ) ( cos(phi_i) ) = 0
      !   where
      !   Qp[m,i] = Imag( T2C[m,i] + (-1)**m * T2C[-m,i] )
      !   Rm[m,i] = Real( T2C[m,i] - (-1)**m * T2C[-m,i] )
      !   Rp[m,i] = Real( T2C[m,i] + (-1)**m * T2C[-m,i] )
      !   Qm[m,i] = Imag(-T2C[m,i] + (-1)**m * T2C[-m,i] )

      call utils_assert((size(T2C, 1) == 2*LL + 1 .and. size(T2C, 2) == channels), "Array shape mismatch in phase T2C routine")

      do i = 1, channels
         ctg = 0.
         write (*, '(a,100f8.3)') 'scanning line : ', T2C(:, i)
         do m = 0, LL
            Qp = aimag(T2C(m + LL + 1, i)) + (-1)**m*aimag(T2C(-m + LL + 1, i))
            Rm = real(T2C(m + LL + 1, i)) - (-1)**m*real(T2C(-m + LL + 1, i))
            Rp = real(T2C(m + LL + 1, i)) + (-1)**m*real(T2C(-m + LL + 1, i))
            Qm = -aimag(T2C(m + LL + 1, i)) + (-1)**m*aimag(T2C(-m + LL + 1, i))
            write (*, '(a,4f13.4)') 'Qp Rm Rp Qm : ', Qp, Rm, Rp, Qm
            if (abs(Qp) > 1.d-5 .or. abs(Rm) > 1.d-5) then
               if (abs(Qp) > 1.d-5) then
                  ctg = -Rp/Qp
                  xb = -Rp
                  xa = Qp
               endif
               if (abs(Rm) > 1.d-5) then
                  ctg = -Qm/Rm
                  xb = -Qm
                  xa = Rm
               endif
            endif
         enddo

         write (*, *) 'i,ctg : ', i, ctg

         if (abs(ctg) > 1.d-5) then
            do m = 0, LL
               Qp = aimag(T2C(m + LL + 1, i)) + (-1)**m*aimag(T2C(-m + LL + 1, i))
               Rm = real(T2C(m + LL + 1, i)) - (-1)**m*real(T2C(-m + LL + 1, i))
               Rp = real(T2C(m + LL + 1, i)) + (-1)**m*real(T2C(-m + LL + 1, i))
               Qm = -aimag(T2C(m + LL + 1, i)) + (-1)**m*aimag(T2C(-m + LL + 1, i))
               if (abs(Rp + Qp*ctg) > 1.d-5 .or. abs(Qm + Rm*ctg) > 1.d-5) then
                  write (*, *) 'ERROR: Could not find an angle to make all cubic harmonics real'
               endif
            enddo
            phi = atan2(xa, xb)
            write (*, *) 'found a coefficient : ', phi
            write (*, *) 'T2C before : ', T2C(:, i)
            T2C(:, i) = T2C(:, i)*MPLX(phi)
            write (*, *) 'T2C after  : ', T2C(:, i)
         endif

      enddo

      T2Cp = transpose(T2C)

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine compute_Slater_interaction(U, J)
      implicit none
      INTEGER                 :: i1, i2, i3, i4, m1, m2, m3, m4, shft
      INTEGER                 :: nw
      complex(kind=DP)              :: dsum
      real(kind=DP)                 :: U, J
      integer                 :: imin, imax, kkk, kkk_max

      call cmp_all_Gaunt(gck)
      nw = 2*LL + 1
      if (allocated(UC)) deallocate (UC)
      allocate (UC(0:LL, channels, channels, channels, channels))
      UC = 0.d0

      if (dmft_for_dimer) then
         kkk_max = 2
      else
         kkk_max = 1
      endif

      do kkk = 1, kkk_max

         if (dmft_for_dimer) then
            if (kkk == 1) then
               imin = 1
               imax = channels/2
            else
               imin = channels/2 + 1
               imax = channels
            endif
         else
            imin = 1
            imax = channels
         endif

         do i4 = imin, imax
         do i3 = imin, imax
         do i2 = imin, imax
         do i1 = imin, imax
         do k = 0, LL
            dsum = 0
            !gck -> m4 -L:L <=  m4 0:2*L
            do m4 = -LL, LL
            do m3 = -LL, LL
            do m2 = -LL, LL
            do m1 = -LL, LL
               if (m1 + m2 /= m3 + m4) goto 25
        dsum = dsum + T2Cp(i4,m4+LL+1)*gck(LL,m4,m1,k)*T2C(m1+LL+1,i1) * T2Cp(i3,m3+LL+1)*gck(LL,m2,m3,k)*T2C(m2+LL+1,i2)
25          enddo
            enddo
            enddo
            enddo
            UC(k, i4, i3, i2, i1) = dsum
         enddo
         enddo
         enddo
         enddo
         enddo

      enddo

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine generate_UC_dat_files(U, Jh)

      use common_def, only: utils_assert

      implicit none
      real(kind=DP) :: U, Jh
      integer :: i, j, k, checkit

      if (solver < 4) then
         call utils_system_call("cat `which atom_d.py` |grep onetep |head -n 1  | wc -l >> check_script ")
         open (unit=30112, file='check_script')
         read (30112, *) checkit
         close (30112)
         call utils_system_call("rm check_script")
         call utils_assert(checkit /= 0, 'Your CTQMC code (in particular the atom_d.py script) was not updated for onetep')
         open (unit=30112, file='Uc.dat.onetep')
         write (30112, *) U
         close (30112)
      endif

      if (solver == 4) then
         call SlaterF(Fn, U, Jh)
         call compute_Slater_interaction(U, Jh)
         write (*, *) 'BEFORE SUMMATION ON K'
         write (*, *) 'Coulomb interaction is complex? max imaginary part?', maxval(abs(aimag(UC))), maxloc(abs(aimag(UC)))
         write (*, *) '                                     max real part?', maxval(abs(real(UC))), maxloc(abs(real(UC)))
         UC(0, :, :, :, :) = UC(0, :, :, :, :)*Fn(0, LL)*0.5d0
         do k = 1, LL
            UC(0, :, :, :, :) = UC(0, :, :, :, :) + UC(k, :, :, :, :)*Fn(k, LL)*0.5d0
         enddo
         write (*, *) 'Slater coefficients F0,F2,F4,F6 : ', Fn(0:LL, LL)
         write (*, *) 'AFTER SUMMATION ON K'
         write (*, *) 'Coulomb interaction is complex? max imaginary part?', maxval(abs(aimag(UC(0, :, :, :, :))))
         write (*, *) '                                     max real part?', maxval(abs(real(UC(0, :, :, :, :))))
         open (unit=30112, file='Uc.dat.ED', form='unformatted')
         write (30112) channels
         write (30112) UC(0, :, :, :, :)
         close (30112)
         write (*, *) 'Now computing impurity many-body Hamiltonian matrix'
         if (use_precomputed_slater_matrix) call generate_slater_many_body_matrix
         write (*, *) 'done....'
      endif

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

   subroutine generate_slater_many_body_matrix
      implicit none
      LOGICAL                     :: bn(2), bn1(2), bn2(2)
      INTEGER                     :: jj, k_, kmax
      INTEGER                     :: istate_, jstate_, iup, idn
      INTEGER                     :: up_in_, do_in_, i1, i2, do_out_, up_out_, fermion_sign_
      integer                     :: fermion_signb_, site1, site2, site3, site4, n3, n4, is1, is2, up_outb_, do_outb_
      logical                     :: bn3(2), bn4(2)
      integer                     :: up_outbb_, do_outbb_, up_outbbb_, do_outbbb_, fermion_signbb_, fermion_signbbb_
      integer                     :: nstate, Nc

      if (allocated(UCCr)) deallocate (UCC, UCCr)

      kmax = -10000; nstate = (2**(channels)); Nc = channels

      if (.false.) then
180      continue
         allocate (UCC(0:nstate - 1, 0:nstate - 1, kmax, 3), UCCr(0:nstate - 1, 0:nstate - 1, kmax))
         UCC(0:, 0:, 1:, 1:3) = 0
      endif

      DO iup = 0, nstate - 1
      DO idn = 0, nstate - 1

         k_ = 0

         DO site4 = 1, Nc
            bn4 = (/BTEST(iup, site4 - 1), BTEST(idn, site4 - 1)/)
            up_in_ = iup
            do_in_ = idn
            do is2 = 1, 2
               if (.not. bn4(is2)) cycle
               if (is2 == 1) then
                  up_out_ = IBCLR(up_in_, site4 - 1)
                  do_out_ = do_in_
                  fermion_sign_ = 1
                  i2 = site4
                  DO jj = 1, i2 - 1; IF (BTEST(up_in_, jj - 1)) fermion_sign_ = -fermion_sign_; ENDDO
               else
                  do_out_ = IBCLR(do_in_, site4 - 1)
                  up_out_ = up_in_
                  fermion_sign_ = 1
                  i2 = site4
                  DO jj = 1, i2 - 1; IF (BTEST(do_in_, jj - 1)) fermion_sign_ = -fermion_sign_; ENDDO
               endif

               DO site3 = 1, Nc
                  bn3 = (/BTEST(up_out_, site3 - 1), BTEST(do_out_, site3 - 1)/)
                  do is1 = 1, 2
                     if (.not. bn3(is1)) cycle
                     if (is1 == 1) then
                        up_outb_ = IBCLR(up_out_, site3 - 1)
                        do_outb_ = do_out_
                     else
                        do_outb_ = IBCLR(do_out_, site3 - 1)
                        up_outb_ = up_out_
                     endif

                     DO site2 = 1, Nc
                        bn2 = (/BTEST(up_outb_, site2 - 1), BTEST(do_outb_, site2 - 1)/)
                        if (bn2(is1)) cycle
                        if (is1 == 1) then
                           up_outbb_ = IBSET(up_outb_, site2 - 1)
                           do_outbb_ = do_outb_
                           fermion_signbb_ = fermion_sign_
                           i1 = min(site2, site3); i2 = max(site2, site3)
                           DO jj = i1 + 1, i2 - 1; IF (BTEST(up_outb_, jj - 1)) fermion_signbb_ = -fermion_signbb_; ENDDO
                        else
                           up_outbb_ = up_outb_
                           do_outbb_ = IBSET(do_outb_, site2 - 1)
                           fermion_signbb_ = fermion_sign_
                           i1 = min(site2, site3); i2 = max(site2, site3)
                           DO jj = i1 + 1, i2 - 1; IF (BTEST(do_outb_, jj - 1)) fermion_signbb_ = -fermion_signbb_; ENDDO
                        endif

                        DO site1 = 1, Nc
                           bn1 = (/BTEST(up_outbb_, site1 - 1), BTEST(do_outbb_, site1 - 1)/)
                           if (bn1(is2)) cycle
                           if (is2 == 1) then
                              up_outbbb_ = IBSET(up_outbb_, site1 - 1)
                              do_outbbb_ = do_outbb_
                              fermion_signbbb_ = fermion_signbb_
                              i2 = site1
                              DO jj = 1, i2 - 1; IF (BTEST(up_outbb_, jj - 1)) fermion_signbbb_ = -fermion_signbbb_; ENDDO
                           else
                              up_outbbb_ = up_outbb_
                              do_outbbb_ = IBSET(do_outbb_, site1 - 1)
                              fermion_signbbb_ = fermion_signbb_
                              i2 = site1
                              DO jj = 1, i2 - 1; IF (BTEST(do_outbb_, jj - 1)) fermion_signbbb_ = -fermion_signbbb_; ENDDO
                           endif

                           if (abs(UC(0, site1, site2, site3, site4)) > slater_inter_cutoff) then
                              k_ = k_ + 1
                              if (k_ > kmax) kmax = k_
                              if (allocated(UCC)) then
                                 UCC(iup, idn, k_, 1) = up_outbbb_
                                 UCC(iup, idn, k_, 2) = do_outbbb_
                                 UCC(iup, idn, 1, 3) = k_
                                 UCCr(iup, idn, k_) = dble(fermion_signbbb_)*conjg(UC(0, site1, site2, site3, site4)) !for Lanczos, which will use its conjugate
                              endif
                           endif
                        enddo !site1

                     enddo !site2

                  enddo !is1
               enddo !site3

            enddo !is2
         ENDDO !site4

      ENDDO
      ENDDO

      if (.not. allocated(UCC)) then
         write (*, *) 'THERE ARE [x] MAX CONNECTIONS : ', kmax
         goto 180
      endif

      do i = 0, nstate - 1
      do j = 1, UCC(i, 1, 1, 3)
         write (*, '(a,2i6,4f10.4)') 'istate up connected to : ', UCC(i, 1, j, 1:2), UCCr(i, 2, j)
      enddo
      enddo

      open (unit=737, file='Uc.dat.ED.matrix', form='unformatted')
      kmax = size(UCC, 3)
      write (737) nstate, nstate, kmax, 3
      write (737) UCC(0:nstate - 1, 0:nstate - 1, 1:kmax, 1:3)
      write (737) Nc, nstate, nstate, kmax
      write (737) UCCr(0:nstate - 1, 0:nstate - 1, 1:kmax)
      close (737)

   end subroutine

   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!
   !====================!

end program

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
!==============================================================!
