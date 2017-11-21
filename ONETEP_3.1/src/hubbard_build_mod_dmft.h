!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hubbard_dmft_interface(eigen_en_input, n_occ_input, &
       ham,overlap_input,inv_overlap_input,ngwf_basis,hub_proj_basis,hub,rep,elements,denskern,proj_basis,dmft_energy_cor,dmft_self,dmft_kernel,dmft_z)

    !=========================================================================!
    ! This subroutine calculates the projection of the Kohn-Sham Greens       !
    ! G(i omega)=([mu + i omega]S - H )^-1 onto each DFT+U correlate subspace.!
    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!
!    use IFPORT, only: system
    use comms, only: comms_abort, comms_barrier, comms_reduce, &
         pub_my_node_id, pub_on_root,pub_total_num_nodes, comms_bcast
    use simulation_cell, only: pub_cell
    use integrals, only: integrals_grad
    use constants, only: DP, UP, DN, PI, stdout, max_spins, HARTREE_IN_EVS
    use wrappers, only : wrappers_invert_sym_cmatrix
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use dense, only: dense_eigensolve
    use ion, only: ELEMENT
    use optics, only : optics_grad_mat_els
    use parallel_strategy, only: pub_first_atom_on_node, &
         pub_num_atoms_on_node, pub_distr_atom, &
         pub_num_hub_atoms_on_node, pub_hub_atoms_on_node,pub_elements_on_node,pub_orig_atom,pub_node_of_atom
    use rundat, only: max_resid_hotelling, maxit_hotelling, pub_rootname, &
         pub_dmft_chem_shift,pub_dmft_smear,pub_dmft_smear_T,pub_dmft_smear_shift,pub_dmft_smear_eta,pub_dmft_smear_w, pub_dmft_points,  &
         pub_dmft_emin,pub_dmft_emax,pub_dmft_temp,pub_dmft_paramagnetic, &
         pub_dmft_cutoff_small,pub_dmft_cutoff_tail,pub_dmft_rotate_green,pub_dmft_free_green_frequ, &
         pub_dmft_plot_real_space,pub_dmft_plot_real_space_sigma,pub_dmft_optics,pub_dmft_optics_window,pub_dmft_optics_i1,pub_dmft_optics_i2, &
         pub_dmft_dos_min,pub_dmft_dos_max,pub_dmft_norm_proj,pub_dmft_plot_all_proj,pub_dmft_integrate_green,pub_dmft_optics_x1,pub_dmft_optics_y1, &
         pub_dmft_optics_z1,pub_dmft_e_plot,pub_dmft_in_bohr,pub_dmft_nmu_loop,pub_dmft_sc,pub_cond_calculate,pub_dmft_kernel,pub_rootname, &
         pub_dmft_kernel_mix,pub_dmft_KS_shift,pub_dmft_fully_sc_h,pub_kernel_diis_max,pub_kernel_diis_c_in,pub_kernel_diis_c_out,&
         pub_dmft_doping,pub_dmft_mu_diff_max,pub_dmft_embed_iter,pub_dmft_embed_mix,pub_dmft_nkpoints,pub_dmft_mu_order,pub_dmft_skip_energy,pub_dmft_kpoints_sym,&
         pub_dmft_kpoints_kernel_gamma,pub_dmft_lin_scaling,pub_dmft_nval,pub_dmft_win,pub_dmft_scaling_maxspace, pub_dmft_scaling_iter,pub_dmft_scaling_cutoff,&
         pub_dmft_scaling_tol,pub_dmft_scaling_nmpi,pub_dmft_scaling_meth,pub_dmft_scaling_tail,pub_dmft_split,pub_dmft_invert_overlap,pub_dmft_local_scratch, &
         pub_dmft_splitk,pub_dmft_scaling_cutoff_h,pub_dmft_impose_same_coeffs,pub_dmft_impose_chem_spin,pub_dmft_ignore_type
    use restart, only: restart_kernel_write,restart_kernel_read
#ifdef GPU_SPEEDUP
    use fortran_cuda, only : diago_cuda_it_c
#endif
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_axpy, sparse_copy, sparse_create, &
         sparse_destroy, sparse_get_element, sparse_put_element, &
         sparse_hotelling_init, sparse_hotelling_invert, sparse_product, &
         sparse_scale, sparse_trace, sparse_show_matrix, sparse_read, sparse_write,sparse_convert,sparse_num_rows
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit
    use ngwf_representation, only: NGWF_REP
    use ion, only: ELEMENT
    use kernel, only: kernel_purify,kernel_init_core_ham
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!

    implicit none

#define hub_overlap_matrix   hub_overlap
#define hub_overlap_tmatrix  hub_overlap_t
#define pub_hub_nat          pub_cell%nat_hub

    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!
    type(NGWF_REP), intent(in)        :: rep
    type(ELEMENT), intent(in)         :: elements(:)
    real(kind=DP), intent(in)         :: eigen_en_input(:,:)       ! molecular orbital energies
    integer, intent(in)               :: n_occ_input(:)            ! number of occupied orbitals
    real(kind=DP),allocatable         :: eigen_en(:,:)             ! molecular orbital energies
    integer,allocatable               :: n_occ(:)                  ! number of occupied orbitals
    integer                           :: ii,num_spins,i,n1_eigen,n2_eigen,n3_eigen
    integer,allocatable               :: dimer_table(:,:),merge_table(:,:),merge_neigh(:)
    real(8), allocatable              :: nabla_square(:,:,:),temp1_real(:,:),temp2_real(:,:)
    complex(8),allocatable            :: temp2_complex(:,:)
    type(SPAM3)                       :: nabla(3)
    type(SPAM3), intent(inout)        ::       ham(pub_cell%num_spins)
    type(SPAM3)                       :: ham_embed(pub_cell%num_spins)
    type(SPAM3)                       :: hamK
    type(SPAM3)                       :: ovK,invK,selfK,selfK1,selfK2,selfK3,O_bar,self_energy_vK
    type(SPAM3), intent(in)           :: overlap_input
    type(SPAM3), intent(in)           :: inv_overlap_input
    type(SPAM3)                       :: overlap
    type(SPAM3)                       :: inv_overlap
    type(FUNC_BASIS), intent(in)      :: ngwf_basis
    type(FUNC_BASIS), optional        :: proj_basis
    type(FUNC_BASIS), intent(in)      :: hub_proj_basis
    type(HUBBARD_MODEL),intent(inout) :: hub
    complex(kind=DP),allocatable      :: site_greenf_buffer(:,:,:,:)
    complex(kind=DP)                  :: greenf_element
    type(SPAM3)                       :: greenf,greenf_backup
    type(SPAM3)                       :: inv_greenf,inv_greenf_backup
    type(SPAM3)                       :: overlap_cplx
    type(SPAM3)                       :: hub_overlap_matrix
    type(SPAM3)                       :: hub_overlap_tmatrix
    type(SPAM3)                       :: hub_overlap_matrix_cplx
    type(SPAM3)                       :: hub_overlap_tmatrix_cplx
    type(SPAM3)                       :: w_greenf
    type(SPAM3)                       :: w_greenf_v,w_S_v,invOk,Ok
    type(SPAM3)                       :: hub_inv_overlap_matrix
    type(SPAM3)                       :: hub_inv_overlap_tmatrix
    type(SPAM3)                       :: hub_inv_overlap_matrix_cplx
    type(SPAM3)                       :: hub_inv_overlap_tmatrix_cplx
    type(SPAM3)                       :: w_inv_greenf
    type(SPAM3)                       :: w_inv_greenf_v
    complex(kind=DP), allocatable     :: site_self_energy_buffer(:,:,:,:)
    complex(kind=DP)                  :: self_energy_element
    type(SPAM3)                       :: self_energy
    type(SPAM3)                       :: self_energy_v,greenf_v
    type(SPAM3)                       :: w_self_energy_v 
    integer                           :: is, ierr
    integer                           :: ien
    real(kind=DP)                     :: en_start, en_range, en_step
    complex(kind=DP)                  :: energy
    logical                           :: check,real_space_exit
    integer                           :: mpi_size,mpi_rank
    character(200)                    :: filename
    integer                           :: kk_,error
    integer                           :: num,len__,status__
    real(kind=DP)                     :: dist
    character(200)                    :: value(100)
    logical                           :: cluster,show_matrices
    integer                           :: cluster_scheme,fopen,j,totn
    integer                           :: hub_atomb,hatb,theatomb,channelsb,spb
    real(kind=DP)                     :: max_resid
    real(kind=DP)                     :: fermi_e(2)
    real(kind=DP)                     :: temp_times_pi, en_step_over_pi
    real(kind=DP)                     :: mytrace,profile_ratio
    character(len=32)                 :: string 
    integer,parameter                 :: nmerge=3
    complex(8),allocatable            :: self_temp(:,:)
    integer                           :: k,l,k2,l2,i2,jj,ressys
#ifdef debug
    logical,parameter                 :: verboseall=.true.
#else
#ifndef debugMIXv2
    logical,parameter                 :: verboseall=.false.
#else
    logical,parameter                 :: verboseall=.true.
#endif
#endif
    logical,allocatable               :: connection_table(:,:)
    integer,allocatable               :: the_chosen(:)
    real(8),allocatable               :: rot_vec_angles(:,:,:),rot_vec_angles_split(:,:,:)
    integer                           :: hub_atom, hat, theatom, channels, sp, kk
    integer                           :: row_proj, col_proj
    integer                           :: row_counter, col_counter
    integer,allocatable               :: myfrequencies(:)
    integer                           :: ffirst,llast,iien,myfrequencies_numb,energyzero
    logical                           :: check_profile
    complex(8),allocatable            :: overlap_hubbard(:,:,:)
    real(kind=DP), allocatable        :: sigma_density_fine(:,:,:,:)
    real(kind=DP), allocatable        :: green_density_fine(:,:,:,:),green_density_fine_proj(:,:,:,:),green_tot_fine(:,:,:,:)
    real(kind=DP), allocatable        :: Z_density_fine(:,:,:,:)
    type(SPAM3)                       :: Z_density_spam(2),green_density_spam(2),green_density_spam_proj(2),sigma_density_spam(2),green_tot_spam(2)
    integer                           :: gpu_max
    integer                           :: sign_temp
    logical                           :: use_gpu,use_gpu_onlyme
    logical                           :: reopen_sigma_file_each_step
    logical                           :: restrict_window
    integer,parameter                 :: nspecial_frequ=7
    integer                           :: pub_hub_nproj
    TYPE(array_of_matrices)           :: h_atoms(pub_cell%nat_hub) 
    real(kind=DP)                     :: E_sum_eigenvalues_DFT(2), total_energy, lhxc_energy, hubbard_energy, dmft_Energy(2)
    integer                           :: kkk1,kkk2,nmu_step
    real(8)                           :: target_N
    integer                           :: nmu
    real(8),parameter                 :: chem_left_shift=0.1d0,chem_right_shift=0.1d0
    real(8)                           :: chem_left,chem_right,chem_shift
    real(8),allocatable               :: vectorAA(:,:)
#ifdef debug
    real(8),allocatable               :: vectorAAback(:,:)
#endif
    real(8)                           :: vectorMu(1000),vectorNmu(1000)
    complex(8),allocatable            :: matrixAAc(:,:,:),matrixHc(:,:,:),matrixOverlapc(:,:),matrixContributionTailc(:,:)
    real(8),allocatable               :: matrixAA(:,:,:),matrixOverlap(:,:),matrixH(:,:,:),matrixContributionTail(:,:)
    complex(8),allocatable            :: matrixAAt(:,:,:)
    real(8),allocatable               :: matrixDensKernel(:,:),tmpr(:,:),tmpr2(:,:)
    complex(8),allocatable            :: matrixDensKernelc(:,:)
    complex(8),allocatable            :: tmp2(:,:),matrixAsig(:,:,:),sigmabackup(:,:),greenfbackup(:,:),tmp1(:,:),matrixSigStat(:,:,:),matrixGan(:,:)
    logical                           :: force_reopen_file,greenfbackupflag
    type(SPAM3), intent(inout)        :: denskern(pub_cell%num_spins)
    type(SPAM3)                       :: denskern_tmp(pub_cell%num_spins),pur_denskern_tmp(pub_cell%num_spins),denskern_tmpc
    logical                           :: partial_dens_present
    type(SPAM3)                       :: denskern_partial(pub_cell%num_spins)
    real(8),optional                  :: dmft_energy_cor
    real(8)                           :: dmft_energy_cor_(2),contrib_energy,contrib_N
    logical,parameter                 :: full_prec_dmft_energy=.false.  !BUG should be false
    logical,parameter                 :: use_gpu_eigenvectors=.true.
    type(SPAM3)                       :: self_infinity(pub_cell%num_spins)
    logical                           :: silent=.false.,split
    integer                           :: kstart,kstep,kien,kkien
    type(SPAM3),optional              :: dmft_kernel(pub_cell%num_spins),dmft_self(pub_cell%num_spins),dmft_z(pub_cell%num_spins)
    logical                           :: checkfirstiter
    type(SPAM3)                       :: dkn_in(1:pub_cell%num_spins,  1:pub_kernel_diis_max)
    type(SPAM3)                       :: dkn_out(1:pub_cell%num_spins, 1:pub_kernel_diis_max)
    type(SPAM3)                       :: residues(1:pub_cell%num_spins,1:pub_kernel_diis_max)
    type(SPAM3)                       :: next_dkn_in(1:pub_cell%num_spins)
    integer                           :: ientry,iterdiis,fopenE,fopenEL,fopenER
    logical                           :: shifted,checkdiis,check_embed,check_embed_h
    integer,allocatable               :: maskH1_embed(:,:,:),maskH1_embed_atoms(:,:,:),maskH0_embed(:,:,:),maskH0_embed_atoms(:,:,:),embed_index(:),embed_size(:),embed_pos(:),embed_pos_all(:)
    complex(8),allocatable            :: ttt0L(:,:),ttt0R(:,:)
#ifdef savememory
    complex(8),allocatable            :: SSL_(:,:),SSR_(:,:)
#else
    complex(8),allocatable            :: SSL_(:,:,:),SSR_(:,:,:)
#endif
    real(8),allocatable               :: mmuRar(:),mmuLar(:)
    complex(8),allocatable            :: wwRar(:),wwLar(:)
    real(8),allocatable               :: S0L_(:,:),S1L_(:,:),H0L_(:,:),H1L_(:,:)
    real(8),allocatable               :: S0R_(:,:),S1R_(:,:),H0R_(:,:),H1R_(:,:)
    integer,allocatable               :: maskTL(:,:,:),maskTLatom(:,:,:),mask_vectorL(:),maskTR(:,:,:),maskTRatom(:,:,:),mask_vectorR(:)
    real(8),allocatable               :: ttL_(:,:),tttL_(:,:),ttR_(:,:),tttR_(:,:)
    complex(8),allocatable            :: ttRc_(:,:),tttRc_(:,:),ttLc_(:,:),tttLc_(:,:)
    complex(8),allocatable            :: GGL_(:,:),GGR_(:,:),GGLback_(:,:),GGRback_(:,:)
    integer                           :: size_indexL,size_indexR,size_mat0L,size_mat0R
    integer                           :: nplanes
    logical                           :: planes_check
    logical,allocatable               :: orb_in_plane(:,:),atom_in_plane(:,:)
    complex(8),allocatable            :: tmpc(:,:)
    real(8),parameter                 :: slope_mu_mixing=1.0,slope_mu_max_step=0.10
   !real(8),parameter                 :: slope_mu_mixing=0.3,slope_mu_max_step=0.01
    logical                           :: wont_compute_residues,use_gpu_some_only
    integer                           :: sizdiis,ikpoint,nkpoints,lastkpoint,kp_file_shift,ii1,ii2,ii3,jp
    real(8)                           :: tempvar,totkpoint
    real(8),allocatable               :: norm_kpoints(:),kpt(:,:)
    real(8)                           :: norm,t1_,t2_,t3_
    complex(8),allocatable            :: tmpc2(:,:)
    logical                           :: same_self_for_all
#ifdef debug
    complex(8),allocatable            :: tmpc3(:,:)
#endif
#ifdef LINEAR_SCALING_INV
    complex(kind=DP), allocatable     :: lin_ovSq(:,:),lin_HSq(:,:),lin_sigSq(:,:),lin_GSq(:,:),lin_vecs(:,:) 
    real(8),allocatable               :: lin_HrSq(:,:)
    real(8),allocatable               :: lin_vec(:)
    integer                           :: lin_nconv_back
    real(8)                           :: lin_rconv_back
#endif
    complex(8),allocatable            :: split_hubt(:,:),split_hub(:,:),split_wgv(:,:),split_Sm(:,:),split_S(:,:),split_H(:,:)
    real(8),allocatable               :: split_Hr(:,:),split_Sr(:,:)
    character(2000)                   :: scratchdir
    logical                           :: breaknfscom
    integer                           :: splitk_start,splitk_step,splitk_end,splitk_last_batch,dmft_splitk_iter,splitk_first
    logical,parameter                 :: one_k_one_file=.false.,optimisation_parallel=.true.
    logical                           :: inv_projectors
    logical,parameter                 :: plot_all_orbitals=.true. !to plot the DOS of all orbitals 
    logical,parameter                 :: gbooth=.false.
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!

    call set_parallel_scheme

    call set_kpoints

    call catch_input_arguments_and_set_module_variables

    call define_cluster

    call init_routine_variables

    call init_scratch_dir

    call init_sparse_and_dimers

    call init_optics_calculation
    
    call check_if_we_can_catch_hamiltonian_in_a_file

    call fix_fermi_and_total_energies

    call embedding_prepare_tables

    call define_my_frequencies

    call embedding_setup_planes

    if(pub_dmft_sc) then
      write(*,*) 'WE LEAVE HERE THE DMFT' 
      write(*,*) 'THE PURPOSE IS TO GENERATE THE STORE FILES'
      write(*,*) 'FOR FURTHER DMFT CALCULATIONS'
      goto 9801 !aborting
    endif


    !---------!
    !---------!
    !---------!
    !---------!
    !---------!
    !---------!
    !---------!
    !---------!


    nmu_step=0 ; vectorMu=0.0 ; vectorNmu=0.0 ; chem_shift=0.d0
    2027 continue

    if(pub_dmft_temp>0.0_DP.and..not.pub_dmft_impose_chem_spin) call fix_chemical_potential_loop

    do is=1,num_spins

    if(pub_dmft_temp>0.0_DP.and.pub_dmft_impose_chem_spin) then; nmu_step=0 ; vectorMu=0.0 ; vectorNmu=0.0 ; chem_shift=0.d0 ; endif;
    2029 continue
    if(pub_dmft_temp>0.0_DP.and.pub_dmft_impose_chem_spin) call fix_chemical_potential_loop_spin

    if(pub_dmft_splitk)then
      splitk_start=mpi_rank-1; splitk_step=mpi_size
      splitk_end=nkpoints/mpi_size+1 !number of batch
      splitk_end=splitk_end*mpi_size
      splitk_last_batch=splitk_end-mpi_size+1 !>= splitk_last_batch is in last batch
      lastkpoint=splitk_last_batch
      dmft_splitk_iter=0
      if(pub_my_node_id==0)then
         write(*,*) 'SPLITTING K-POINTS - START/END - my rank : ', splitk_start,splitk_end,mpi_rank
         write(*,*) 'SPLITTING K-POINTS - STEP                : ', splitk_step
         write(*,*) 'SPLITTING K-POINTS - LAST BATCH START AT : ', splitk_last_batch
         write(*,*) 'BREAK NFS COMM                           : ', breaknfscom
      endif
      splitk_first=nkpoints+1
      do ikpoint=splitk_start+1,splitk_end,splitk_step
         if(abs(norm_kpoints(ikpoint))>1.d-4) then
           splitk_first=ikpoint
           exit
         endif
      enddo
      write(*,*) 'SPLITTING K-POINTS - FIRST K POINT : ', splitk_first
    else
     splitk_start=0; splitk_step=1;splitk_end=nkpoints
    endif

    do ikpoint=splitk_start+1,splitk_end,splitk_step

       call reset_kernel_and_energy

       if(pub_dmft_splitk)then
         if(abs(norm_kpoints(ikpoint))<1.d-4) then    
           write(*,*) 'DMFT SPLIT K : skipping k-point, due to symmetry'
           goto 1030
         endif
         if(ikpoint>nkpoints)then
           write(*,*) 'DMFT SPLIT K : skipping k-point - outbound index'
           goto 1030
         endif
       endif

       if(abs(norm_kpoints(ikpoint))<1.d-4) then
         write(*,*) 'SKIPPING K POINT DUE TO SYMMETRY'; cycle 
       endif

       write(*,*)            'DOING KPOINT NUMBER : ', ikpoint
       if(ikpoint<=nkpoints)then
        write(*,'(a,3f14.4)') ' KP                 : ', kpt(ikpoint,:)
        write(*,*)            'NORM                : ', norm_kpoints(ikpoint)
       else
        write(*,'(a,3f14.4)') ' KP-out-of-range    : ', ikpoint,lastkpoint
       endif

       call load_k_ham_and_overlap

       if(pub_dmft_temp>0.0_DP .and. nmu>0) call prepare_for_dmft_density_kernel
    
       real_space_exit=.false.
 
#ifdef debug 
       call compute_projected_overlap_matrix_Sproj
#endif

       if(split)then  
          kstart=pub_my_node_id; kstep=pub_total_num_nodes
          if(mod(pub_dmft_points,kstep)/=0) then
             write(*,*) 'total number of nodes : ',  kstep
             write(*,*) 'total number of frequ : ',  pub_dmft_points
             write(*,*) 'the two should be commensurate, abort'
             goto 9801 !aborting 
          endif 
       else
          kstart=0; kstep=1 
       endif

       do iien=kstart+1,myfrequencies_numb,kstep

          ien=myfrequencies(iien) 

          if(present(dmft_energy_cor).and.split)then
             if(iien>pub_dmft_points-nspecial_frequ) then 
               write(*,*) '...processor done with dmft_energy part so cycle now...'
               exit
             endif
          endif

                          if(pub_my_node_id==0) write(*,*) 'doing frequency       : ', ien
          if(.not.silent.and.pub_my_node_id==0) write(*,*) 'number [x] out of [y] : ', iien,myfrequencies_numb

          call define_energy

          !##################################################################!
          !##################################################################!
          !##################################################################!
          if( .not. restrict_window .or. split .or. check_embed .or.                         &
            &        ( real(energy) - (fermi_e(is)+pub_dmft_chem_shift) > pub_dmft_dos_min   &
            &  .and.   real(energy) - (fermi_e(is)+pub_dmft_chem_shift) < pub_dmft_dos_max ) ) then

          call embedding_read_potential

          call build_inv_green_matrix

          call wrapper_build_self_energy

          if(verboseall) write(*,*) 'Subtract it before computation of Greens function'

          if(.not.split) call show_the_self_matrix

          call wrapper_add_self_and_invert_green
 
          if(ien<pub_dmft_points-nspecial_frequ+1.and.pub_dmft_temp>0.0_DP.and.nmu>0)then 
             call compute_dmft_energy
             if(present(dmft_energy_cor)) cycle
          endif

          call dump_data_optics(ien,is,greenf,energy,fermi_e(is),pub_dmft_chem_shift)
 
          call project_green_function(w_greenf_v,greenf)

#ifdef debug
          if(.not.split) call project_hamiltonian
#endif

          if(.not.split) call comms_barrier
          if(.not.split) call show_the_matrices

          call build_projected_green

          call embedding_write_sigma

          call plotting_real_space_quantities(1)

          if(pub_dmft_plot_real_space.and.real_space_exit) exit

          endif
          !##################################################################!
          !##################################################################!
          !##################################################################!

       end do

       write(*,*) 'LOOP OVER FREQUENCIES DONE'

       if(pub_dmft_temp>0.0_DP.and.nmu>0) call wrapper_comm_reduce_energy_and_kernel

       1030 continue
       if(pub_dmft_temp>0.0_DP.and.nmu>0) call clean_dmft_energy

    end do

    if(pub_dmft_temp>0.0_DP.and.pub_dmft_impose_chem_spin)then
       if(nmu_step < nmu) goto 2029
    endif

    end do

    if(pub_dmft_temp>0.0_DP.and..not.pub_dmft_impose_chem_spin)then
       if(nmu_step < nmu) goto 2027
    endif


    call plotting_real_space_quantities(2)

    !---------!
    !---------!
    !---------!
    !---------!
    !---------!
    !---------!
    !---------!
    !---------!

    wont_compute_residues =  use_gpu_some_only .and. &
                    & .not. (present(dmft_kernel).or.present(dmft_energy_cor))
                    
    call finalize_kernel

    if(use_gpu_some_only.and.pub_dmft_kernel/=0.and.pub_dmft_temp>0.0_DP.and.mpi_size>1)then
      tempvar=0.0 ; 
      if(.not.pub_dmft_splitk)then
       call sync_nodes_by_shared_nfs(tempvar,'sum_tempvar_bef_exit',mpi_rank,mpi_size,longdelay=.true.)
      else
       if(pub_my_node_id==0) call sync_nodes_by_shared_nfs(tempvar,'sum_tempvar_bef_exit',mpi_rank,mpi_size,longdelay=.true.)
       call comms_barrier
      endif
    endif

    if(present(dmft_energy_cor)) dmft_energy_cor=sum(dmft_energy_cor_)

    !---------!
    !---------!
    !---------!
    !---------!
    !---------!
 
     if(pub_my_node_id==0.and.mpi_size>1.and.mpi_rank==1) &
   & ressys=system("rm ./onetep_dmft_confirmation_dmft_density_kernel* onetep_dmft_confirmation_dmft_density_count > /dev/null 2>&1 ")

     call correct_dimer_file

     if(.false.)then
       9801 continue
       write(*,*) 'ABORTING DMFT INTERFACE IN ONETEP'
     endif

     call read_gbooth_write_kernel

     call cleanitall

    return

  contains

#include "hubbard_build_mod_dmft_routines.h"
#include "hubbard_build_mod_dmft_diis.h"
    
  end subroutine hubbard_dmft_interface


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
