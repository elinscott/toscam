!@TODO:
!@@dealloc cheb nodes, cheb coeffs
!@ Add MIC to Chebyshevs

! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                   Hartree-Fock exchange module                 !
!                                                                !
! This module contains subroutines for calculating a Hartree-Fock!
! exchange contribution to the energy using FFT boxes which      !
! coincide with the simulation cell.                             !
!----------------------------------------------------------------!
! Written by Quintin Hill in 2008/9 with supervision by          !
! Chris-Kriton Skylaris.                                         !
!================================================================!

module hf_exchange
#define CUBECC
  use constants, only: DP, stdout
  use geometry, only: point
  use utils, only: utils_trace_in, utils_trace_out, utils_flush
  use rundat, only: pub_hfx_debug

  implicit none

  private

  public :: hf_exchange_calculate
  public :: hf_exchange_set_pub_hfxsw
  public :: hf_exchange_fill_vmatrix
  public :: hf_exchange_init_vmatrix

  !qoh: Parameters of a Bessel function:

  type bessel
     integer        :: lval   ! Quantum number l
     real(kind=DP)  :: qval   ! Set so that j_l(qa) = 0
     real(kind=DP)  :: aval   ! Cut off radius of spherical Bessel function
     real(kind=DP)  :: farpotint ! Radial potential integral when r > a
     real(kind=DP)  :: nearpotint !Constant part of radial potential when r < a
  end type bessel

  type atom_centre
     integer     :: species_number
     integer     :: tbstart(3) ! Start point of tightbox
     integer     :: tbend(3)   ! End point of tightbox
     integer     :: tbcellstart(3) ! Start point of tightbox in cell
     real(kind=DP) :: radius ! jd: Localization radius of the NGWFs
     type(point) :: incell  ! Centre of spherical wave in cell
     type(point) :: intb    ! Centre of spherical wave in tightbox
  end type atom_centre

  integer :: max_num_bessels ! The number of Bessels to be used
  integer :: max_sw_set_size ! The maximum possible SW set size for allocation

  logical :: vmatrix_filled = .false.

  ! qoh: Running options
  logical,parameter :: allowfftonly = .true. ! Use FFT if cell = FFT box
  logical,parameter :: usefftforsc = .true. ! Use FFT method when atom A=atom B
  logical,parameter :: useoverlapmetric = .false. ! Use overlap metric
  logical,parameter :: usenpa = .false.  ! Use numerical pointwise approach
  integer,parameter :: lmin = 0 ! Minimum value of l

  ! qoh: Metric matrix calculation options
  ! qoh: sph_grid_metric true trumps recip_grid_metric true
  ! qoh: If both are false than use the real Cartesian grid
  logical,parameter :: chebyshev_grid_metric = .true.    ! Use chebyshev grid
  logical,parameter :: sph_grid_metric = .false.   ! Use spherical grid
  logical,parameter :: recip_grid_metric = .false.  ! Use recip Cartesian grid
  integer,parameter :: finer = 1 ! Use an x times finer real grid for metric

  ! qoh: Experimental options
  logical,parameter :: shavesw = .false. ! Shave generated spherical waves
  logical,parameter :: singlecentre = .false. ! Put SW only on B in ES metric
  logical,parameter :: force_con = .true. ! Force metric matrix consistency
  logical,parameter :: overlappinginfft = .false. ! If atoms overlap use FFT box

  ! qoh: Debugging output options
  logical,parameter, public :: printdkgrad = .false. ! Print total DK gradient
  logical,parameter :: printcoeffs = .false. ! Print expansion coefficients
  logical,parameter :: printxmat   = .false. ! Print HF matrix elements
  logical,parameter :: printdk     = .false. ! Print density kernel
  logical,parameter :: printvmat   = .false. ! Print metric matrix

  ! qoh: Debugging options
  logical,parameter :: justoneblock = .false.  ! only compute a single atomblock

  ! jd: Internal counters of integration routine, for verbose output
  integer :: hfx_int_total_calls
  integer :: hfx_int_total_cg_instances
  integer :: hfx_int_total_no_cg_instances
  real(kind=DP) :: hfx_int_total_packed
  real(kind=DP) :: hfx_int_total_points_tb1
  real(kind=DP) :: hfx_int_total_large
  real(kind=DP) :: hfx_int_total_margin
  real(kind=DP) :: hfx_int_total_coarse_from
  real(kind=DP) :: hfx_int_total_coarse_to
  real(kind=DP) :: hfx_int_total_zero

  ! jd: Tags for comms
  integer, parameter :: SW_REQUEST_TAG = 10000
  integer, parameter :: SW_DATA_TAG = 100000


contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine hf_exchange_calculate(hfexchange,& !output
       denskern,overlap, ngwfs_on_grid, ngwf_basis, & !input
       elements, full_vmatrix, calcgradient,& !input - last argument if false
       tc_denskern, cov_grad, contra_grad, & !Opt arguments
       precond_func_recip) !Opt arguments

    !==========================================================================!
    ! This subroutine calculates the Hartree-Fock exchange matrix and          !
    ! optionally the exchange contribution to the NGWF gradient.  It loops over!
    ! atom As that are local to this node and their corresponding atom Bs      !
    ! (those with non-zero density kernel atom blocks in the column of A.      !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2008 and January 2009.                        !
    ! Modified by Quintin Hill on 12/02/2009 to make hfexchange symmetric.     !
    !==========================================================================!

    use comms, only: pub_total_num_nodes, pub_my_node_id, pub_on_root, comms_abort !@remove comms_abort
    use constants, only: DP, NORMAL !@remove NORMALS
    use function_basis, only: FUNC_BASIS
    use geometry, only: point
    use ion, only: element
    use parallel_strategy, only: pub_max_atoms_on_node, pub_num_atoms_on_node,&
         pub_first_atom_on_node
    use rundat, only: density_batch_size, pub_hfxsw, &
         pub_hfx_read_xmatrix, pub_hfx_write_xmatrix, pub_rootname, pub_output_detail !@jd remove last 4
    use simulation_cell, only: pub_fftbox, pub_cell, pub_maxtight_pts1, &
         pub_maxtight_pts2, pub_maxtight_pts3
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index, &
         sparse_put_block, sparse_create, sparse_destroy,&
         sparse_scale, sparse_transpose, sparse_axpy, sparse_show_matrix, &
         sparse_read, sparse_write !@jd remove last 2
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert !@removelast
    use xc, only: pub_hfxfraction

    implicit none

    !qoh: Arguments
    type(SPAM3),intent(inout) :: hfexchange(pub_cell%num_spins)! Exchange matrix
    type(SPAM3),   intent(in) :: denskern(pub_cell%num_spins)  ! Density kernel
    type(SPAM3),   intent(in) :: overlap                       ! Overlap matrix
    type(SPAM3),   intent(in) :: full_vmatrix ! Precalculated vmatrix atomblocks
    type(FUNC_BASIS), intent(in)  :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(element), intent(in) :: elements(pub_cell%nat)
    logical,       intent(in) :: calcgradient!Should NGWF gradient be calculated
    !qoh: Arguments for calculating NGWF gradient
    !qoh: (required if calcgradient is true)
    type(SPAM3),  optional,intent(in)    :: tc_denskern(pub_cell%num_spins)! KS
    real(kind=DP),optional,intent(inout) :: cov_grad(ngwf_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP),optional,intent(inout) :: contra_grad&
         (ngwf_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP),optional,intent(in)    :: precond_func_recip(:,:,:)

    !qoh: Module Workspaces
    real(kind=DP), allocatable :: bb_fftbox(:,:,:)
    real(kind=DP), allocatable :: cc_fftbox_sum_batch(:,:,:,:,:)
    real(kind=DP), allocatable :: work1_fftbox(:,:,:)
    real(kind=DP), allocatable :: work2_fftbox(:,:,:)
    real(kind=DP), allocatable :: prod_in_fftbox(:,:,:,:)
    real(kind=DP), allocatable :: ket_in_fftbox(:,:,:,:)
    real(kind=DP), allocatable :: ngwf_on_grid_buffer(:)
    real(kind=DP), allocatable :: k_blk(:,:,:) ! K atomblock
    real(kind=DP), allocatable :: tc_k_blk(:,:,:) ! KS atomblock
    complex(kind=DP), allocatable :: comp_fftbox(:,:,:)
    real(kind=DP), allocatable :: xatomblock(:,:,:) ! Atomblock of hfx matrix
    real(kind=DP), allocatable :: dkrecvbuf(:,:,:) ! Density kernel buffer
    real(kind=DP), allocatable :: dklinrecvbuf(:) ! Density kernel linear buffer

    type(bessel),  allocatable :: sphbessels(:,:) !Spherical bessel data
    real(kind=DP), allocatable :: radtable(:) ! Radius for each species
    ! Set of spherical waves that product is expanded in
    real(kind=DP), allocatable :: swswoverlap(:) ! Overlap vector of
    !set of spherical waves
    real(kind=DP), allocatable :: swcoeff(:,:,:)
    real(kind=DP), allocatable :: swa_tightbox(:,:,:)
    real(kind=DP), allocatable :: tb_batch(:,:,:,:,:)
    real(kind=DP), allocatable :: tb_batch2(:,:,:,:,:)
    complex(kind=DP), allocatable :: tb_zwork(:,:,:)
    !qoh: Plan arrays
    logical, allocatable :: dd_send_plan(:,:)
    logical, allocatable :: dd_recv_plan(:,:)
    logical, allocatable :: cc_send_plan(:,:)
    logical, allocatable :: cc_recv_plan(:,:)
    integer, allocatable :: atomb_on_node(:)
    integer, allocatable :: dds_on_node(:,:)
    integer, allocatable :: atomds_on_node(:,:)
    integer, allocatable :: dds_in_batch(:)

    !qoh: Local Variables
    integer :: k_idxlen ! length of index of denskern
    integer :: o_idxlen ! length of index of overlap
    integer, allocatable, dimension(:) :: denskern_idx  ! Denskern index
    integer, allocatable, dimension(:) :: overlap_idx ! Overlap index
    integer :: local_a ! Local index of atom A
    integer :: atoma ! Global index of atom A
    integer :: atomb ! Global index of atom B
    !qoh: Column of atom A in density kernel
    integer :: k_fstnzbridx_cola ! first kernel non-zero block row index
    integer :: k_lstnzbridx_cola !last kernel non-zero block row index
    integer :: k_nzbridx_cola ! kernel non-zero block row index
    integer :: batch_size ! Size of the batch
    integer :: is ! Spin counter
    integer :: ierr !error code
    real(kind=DP) :: spin_fac !Spin factor
    type(SPAM3) :: transhfx  !Transposed HF exchange matrix

    character(len=512) :: filename !@removeme
    logical :: fileexists !@removeme
    integer, save :: niter = 1 !@removeme


    call utils_trace_in('hf_exchange_calculate')

    batch_size = density_batch_size

    call internal_init_stats()

    !qoh: Get index arrays for density kernel and overlap matrix

    !qoh: Allocate arrays to store indices
    k_idxlen = sparse_index_length(denskern(1))
    o_idxlen = sparse_index_length(overlap)
    allocate(denskern_idx(k_idxlen),stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','denskern_idx',ierr)
    allocate(overlap_idx(o_idxlen),stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','overlap_idx',ierr)

    ! qoh: Populate these arrays
    call sparse_generate_index(denskern_idx,denskern(1))
    call sparse_generate_index(overlap_idx,overlap)

    !qoh: Allocate workspace

    allocate(bb_fftbox(pub_fftbox%total_ld1,pub_fftbox%total_ld2,&
         pub_fftbox%total_pt3), stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','bb_fftbox',ierr)
    allocate(cc_fftbox_sum_batch(pub_fftbox%total_ld1, pub_fftbox%total_ld2,&
         pub_fftbox%total_pt3, pub_cell%num_spins, batch_size), stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','cc_fftbox_sum_batch',ierr)
    allocate(work1_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2, pub_fftbox%total_pt3), stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','work1_fftbox',ierr)
    allocate(work2_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2, pub_fftbox%total_pt3), stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','work2_fftbox',ierr)
    allocate(prod_in_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2, pub_fftbox%total_pt3, 2), stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','prod_in_fftbox',ierr)
    allocate(ket_in_fftbox(pub_fftbox%total_ld1,pub_fftbox%total_ld2,&
         pub_fftbox%total_pt3,pub_cell%num_spins), stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','ket_in_fftbox',ierr)
    allocate(ngwf_on_grid_buffer(ngwf_basis%max_n_ppds_sphere*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','ngwf_on_grid_buffer',ierr)
    allocate(xatomblock(ngwf_basis%max_on_atom, ngwf_basis%max_on_atom,&
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','xatomblock',ierr)
    allocate(k_blk(ngwf_basis%max_on_atom, ngwf_basis%max_on_atom,&
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('hf_exchange_calculate', 'k_blk', ierr)
    allocate(tc_k_blk(ngwf_basis%max_on_atom, ngwf_basis%max_on_atom,&
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('hf_exchange_calculate', 'tc_k_blk', ierr)
    allocate(comp_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2, pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','comp_fftbox',ierr)
    allocate(dkrecvbuf(batch_size, ngwf_basis%max_on_atom, pub_cell%num_spins),&
         stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','dkrecvbuf',ierr)
    allocate(dklinrecvbuf(batch_size*ngwf_basis%max_on_atom*pub_cell%num_spins),&
         stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','dklinrecvbuf',ierr)

    if (pub_hfxsw .or. usenpa) then
       allocate(tb_batch(pub_maxtight_pts1,pub_maxtight_pts2,&
            pub_maxtight_pts3,batch_size,pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('hf_exchange_calculate','tb_batch',ierr)
       allocate(tb_batch2(pub_maxtight_pts1,pub_maxtight_pts2,&
            pub_maxtight_pts3,batch_size,pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('hf_exchange_calculate','tb_batch2',ierr)
    end if

    !qoh: Allocate arrays for plan
    allocate(dd_send_plan(ngwf_basis%max_on_node,0:pub_total_num_nodes-1),&
         stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','dd_send_plan',ierr)
    allocate(dd_recv_plan(ngwf_basis%max_on_node,0:pub_total_num_nodes-1),&
         stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','dd_recv_plan',ierr)
    allocate(cc_send_plan(pub_max_atoms_on_node,0:pub_total_num_nodes-1), &
         stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','cc_send_plan',ierr)
    allocate(cc_recv_plan(pub_max_atoms_on_node,0:pub_total_num_nodes-1), &
         stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','cc_recv_plan',ierr)
    allocate(atomb_on_node(0:pub_total_num_nodes-1), stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','atomb_on_node',ierr)
    allocate(dds_on_node(batch_size,0:pub_total_num_nodes-1), stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','dds_on_node',ierr)
    allocate(atomds_on_node(batch_size,0:pub_total_num_nodes-1), stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','atomds_on_node',ierr)
    allocate(dds_in_batch(batch_size), stat=ierr)
    call utils_alloc_check('hf_exchange_calculate','dds_in_batch',ierr)

    if (pub_hfxsw) then
       ! qoh: Allocate arrays for Bessel and Spherical Wave data
       allocate(radtable(pub_cell%num_species),stat=ierr)
       call utils_alloc_check('hf_exchange_calculate','radtable',ierr)
       ! qoh: Initialise max_num_bessels and max_num_sphwaves
       call hf_exchange_num_sph_functions(radtable)
       allocate(sphbessels(pub_cell%num_species,max_num_bessels),stat=ierr)
       call utils_alloc_check('hf_exchange_calculate','sphbessels',ierr)
       call hf_exchange_sph_bessels_init(sphbessels,radtable)
       allocate(swcoeff(max_sw_set_size,batch_size,pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('hf_exchange_calculate','swcoeff',ierr)
       allocate(swswoverlap(max_sw_set_size),stat=ierr)
       call utils_alloc_check('hf_exchange_calculate','swswoverlap',ierr)
       allocate(swa_tightbox(pub_maxtight_pts1,pub_maxtight_pts2,&
            pub_maxtight_pts3),stat=ierr)
       call utils_alloc_check('hf_exchange_calculate','swa_tightbox',ierr)
       allocate(tb_zwork(pub_maxtight_pts1,pub_maxtight_pts2,&
            pub_maxtight_pts3),stat=ierr)
       call utils_alloc_check('hf_exchange_calculate','tb_zwork',ierr)
    end if

    !qoh: Loop over all atoms on this core
    loop_A: do local_a=1,pub_num_atoms_on_node(pub_my_node_id)
       atoma = pub_first_atom_on_node(pub_my_node_id) + local_a -1

       !qoh: Get kernel and overlap blocks in column of atom A
       k_fstnzbridx_cola = denskern_idx(local_a)
       k_lstnzbridx_cola = denskern_idx(local_a+1) - 1

       !qoh: Loop over atoms (atom B) that have non zero blocks K_{AB}
       loop_B: do k_nzbridx_cola = k_fstnzbridx_cola, k_lstnzbridx_cola
          atomb = denskern_idx(k_nzbridx_cola)
          xatomblock = 0.0_DP

          !qoh: The uglyness that follows is required to deal with the varying
          !qoh: number of arguements required for hf_exchange_atomblock
          if (pub_hfxsw) then
             call hf_exchange_atomblock(xatomblock,atoma, atomb, & !output
                  denskern,denskern_idx,overlap_idx,batch_size,&!in
                  ngwfs_on_grid, ngwf_basis, elements, calcgradient, &!input
                  bb_fftbox, cc_fftbox_sum_batch, work1_fftbox, work2_fftbox, &
                  prod_in_fftbox, ket_in_fftbox, ngwf_on_grid_buffer, & !ws
                  comp_fftbox, k_blk, tc_k_blk, &
                  dkrecvbuf, dklinrecvbuf, & !ws
                  dd_send_plan, dd_recv_plan, cc_send_plan, cc_recv_plan,& !plan
                  atomb_on_node, dds_on_node, & !Used by plan
                  atomds_on_node, dds_in_batch, & !Used by plan.
                  tb_batch, tb_batch2, tb_zwork, swa_tightbox, &!workspace
                  sphbessels, full_vmatrix, swcoeff, swswoverlap,&
                  tc_denskern, cov_grad, contra_grad, & !Grad args
                  precond_func_recip)                   !Grad args
          else if (usenpa) then
             call hf_exchange_atomblock(xatomblock,atoma, atomb, & !output
                  denskern,denskern_idx,overlap_idx,batch_size,&!in
                  ngwfs_on_grid, ngwf_basis, elements, calcgradient, &!input
                  bb_fftbox, cc_fftbox_sum_batch, work1_fftbox, work2_fftbox, &
                  prod_in_fftbox, ket_in_fftbox, ngwf_on_grid_buffer, & !ws
                  comp_fftbox, k_blk, tc_k_blk, &
                  dkrecvbuf, dklinrecvbuf, & !ws
                  dd_send_plan, dd_recv_plan, cc_send_plan, cc_recv_plan,& !plan
                  atomb_on_node, dds_on_node, & !Used by plan
                  atomds_on_node, dds_in_batch, & !Used by plan.
                  tb_batch, tb_batch2, & ! workspace
                  tc_denskern=tc_denskern,&!Grad args
                  cov_grad=cov_grad, contra_grad=contra_grad, & !Grad args
                  precond_func_recip=precond_func_recip)        !Grad args
          else
             call hf_exchange_atomblock(xatomblock,atoma, atomb, & !output
                  denskern,denskern_idx,overlap_idx,batch_size,&!in
                  ngwfs_on_grid, ngwf_basis, elements, calcgradient, &!input
                  bb_fftbox, cc_fftbox_sum_batch, work1_fftbox, work2_fftbox, &
                  prod_in_fftbox, ket_in_fftbox, ngwf_on_grid_buffer, & !ws
                  comp_fftbox, k_blk, tc_k_blk, &
                  dkrecvbuf, dklinrecvbuf, & !ws
                  dd_send_plan, dd_recv_plan, cc_send_plan, cc_recv_plan,& !plan
                  atomb_on_node, dds_on_node, & !Used by plan
                  atomds_on_node, dds_in_batch, & !Used by plan.
                  tc_denskern=tc_denskern,&!Grad args
                  cov_grad=cov_grad, contra_grad=contra_grad, & !Grad args
                  precond_func_recip=precond_func_recip)        !Grad args
          end if

          !qoh: Scale by spin factor
          if (calcgradient) then
             spin_fac = 2.0_DP
          else
             spin_fac = real(pub_cell%num_spins,kind=DP)
          end if

          !qoh: Scale by 0.5 (for Hartree to exchange ratio) and pub_hfxfraction
          xatomblock = xatomblock * 0.5_DP * pub_hfxfraction * spin_fac

          !qoh: Insert atom block in to sparse matrix structure.
          do is=1, pub_cell%num_spins

             !if(xatomblock(1,1,1) /= 0D0) xatomblock(:,:,:) = 1.0_DP !@@@@@@@@@@@@@@@@@@@@
          
             call sparse_put_block(xatomblock(:,:,is), &
                  hfexchange(is), atomb, atoma)
          end do

       end do loop_B

    end do loop_A

    !qoh: This node has finished calculating atomblocks, but NGWFs may need to
    !qoh: be sent to those that have not finished.

    call hf_exchange_sendngwfs(denskern,denskern_idx,overlap_idx,&
         batch_size,ngwfs_on_grid,ngwf_basis,&!in
         bb_fftbox, cc_fftbox_sum_batch, & !ws
         work1_fftbox,work2_fftbox, & !ws
         prod_in_fftbox, ket_in_fftbox,  ngwf_on_grid_buffer, & !ws
         comp_fftbox,k_blk,dkrecvbuf,dklinrecvbuf, & !ws
         dd_send_plan, dd_recv_plan, cc_send_plan, cc_recv_plan, & !plan
         atomb_on_node, dds_on_node, & !Used by plan
         atomds_on_node, dds_in_batch) !Used by plan.

    if (printxmat .and. pub_on_root) then
       print *, ""
       print *, "-XMAT:-------------------------------"
    end if
    if (printxmat) call sparse_show_matrix(hfexchange(1))
    if (printxmat .and. pub_on_root) then
       print *, "-------------------------------------"
       print *, "-DKN:--------------------------------"
    end if
    if (printxmat) call sparse_show_matrix(denskern(1))
    if (printxmat .and. pub_on_root) then
       print *, "-------------------------------------"
    end if

    ! jd: Read xmatrix from a file rather than computing it, if asked to
    if(pub_hfx_read_xmatrix) then
       write(filename,'(2a,i0)') trim(pub_rootname),trim('.xmatrix.'),niter
       if (pub_on_root) then
          write(stdout,'(/3a)',advance='no') &
               'Reading xmatrix (@just 1 spin@, also not using this in &
               &NGWF gradient!@) from file "', trim(filename),'" ...'
       end if

       ! Check that the file exists
       ! ndmh: only root node needs to be able to see the file
       if (pub_on_root) then
          inquire(file=filename,exist=fileexists)
       else
          fileexists = .true.
       end if

       if (fileexists) then
          ! Read density kernel from this file
          call sparse_read(hfexchange(1),trim(filename))
       else
          if (pub_on_root) write(stdout,'(/a/)') ' File not found, quitting.'
          call comms_abort
       end if

       if (pub_on_root) write(stdout,'(a)') ' done'

    end if


    ! jd: Write xmatrix to a file, if asked to
    if(pub_hfx_write_xmatrix) then
       write(filename,'(2a,i0)') trim(pub_rootname),trim('.xmatrix.'),niter
       if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
          write(stdout,'(/3a)',advance='no') &
               'Writing xmatrix to file "', trim(filename),'" ...'
       end if

       call sparse_write(hfexchange(1),trim(filename))

       if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
          write(stdout,'(a)') ' done'
       end if
    end if

    niter = niter + 1 !@@@@@

    !qoh: Make hfexchange matrix symmetric
    call sparse_create(transhfx,hfexchange(1))
    do is = 1,pub_cell%num_spins
       call sparse_scale(hfexchange(is),0.5_DP)
       call sparse_transpose(transhfx,hfexchange(is))
       call sparse_axpy(hfexchange(is),transhfx,1.0_DP)
    end do
    call sparse_destroy(transhfx)

    ! qoh: Deallocate plan arrays
    deallocate(dd_send_plan, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','dd_send_plan',ierr)
    deallocate(dd_recv_plan, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','dd_recv_plan',ierr)
    deallocate(cc_send_plan, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','cc_send_plan',ierr)
    deallocate(cc_recv_plan, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','cc_recv_plan',ierr)
    deallocate(atomb_on_node, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','atomb_on_node',ierr)
    deallocate(dds_on_node, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','dds_on_node',ierr)
    deallocate(atomds_on_node, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','atomds_on_node',ierr)
    deallocate(dds_in_batch, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','dds_in_batch',ierr)

    ! qoh: Deallocate workspace
    deallocate(dklinrecvbuf, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','dklinrecvbuf',ierr)
    deallocate(dkrecvbuf,  stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','dkrecvbuf',ierr)
    deallocate(bb_fftbox, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','bb_fftbox',ierr)
    deallocate(cc_fftbox_sum_batch, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','cc_fftbox_sum_batch',ierr)
    deallocate(work1_fftbox, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','work1_fftbox',ierr)
    deallocate(work2_fftbox, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','work2_fftbox',ierr)
    deallocate(prod_in_fftbox, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','prod_in_fftbox',ierr)
    deallocate(ket_in_fftbox, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','ket_in_fftbox',ierr)
    deallocate(ngwf_on_grid_buffer, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','ngwf_on_grid_buffer',&
         ierr)
    deallocate(xatomblock, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','xatomblock',ierr)
    deallocate(k_blk, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate', 'k_blk', ierr)
    deallocate(tc_k_blk, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate', 'tc_k_blk', ierr)
    deallocate(comp_fftbox,stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','comp_fftbox',ierr)

    if (pub_hfxsw) then
       deallocate(radtable,stat=ierr)
       call utils_dealloc_check('hf_exchange_calculate','radtable',ierr)
       deallocate(sphbessels,stat=ierr)
       call utils_dealloc_check('hf_exchange_calculate','sphbessels',ierr)
       deallocate(swcoeff,stat=ierr)
       call utils_dealloc_check('hf_exchange_calculate','swcoeff',ierr)
       deallocate(swswoverlap,stat=ierr)
       call utils_dealloc_check('hf_exchange_calculate','swswoverlap',ierr)
       deallocate(swa_tightbox,stat=ierr)
       call utils_dealloc_check('hf_exchange_calculate','swa_tightbox',ierr)
       deallocate(tb_zwork,stat=ierr)
       call utils_dealloc_check('hf_exchange_calculate','tb_zwork',ierr)
    end if

    if (pub_hfxsw .or. usenpa) then
       deallocate(tb_batch,stat=ierr)
       call utils_dealloc_check('hf_exchange_calculate','tb_batch',ierr)
       deallocate(tb_batch2,stat=ierr)
       call utils_dealloc_check('hf_exchange_calculate','tb_batch2',ierr)
    end if

    !qoh: Deallocate index arrays
    deallocate(overlap_idx, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','overlap_idx',ierr)
    deallocate(denskern_idx, stat=ierr)
    call utils_dealloc_check('hf_exchange_calculate','denskern_idx',ierr)

    call internal_print_stats

    call utils_trace_out('hf_exchange_calculate')

  contains

    subroutine internal_init_stats()

      ! jd: Initialize internal counters for verbose output
      hfx_int_total_calls = 0
      hfx_int_total_cg_instances = 0
      hfx_int_total_no_cg_instances = 0
      hfx_int_total_packed = 0.0_DP
      hfx_int_total_packed = 0.0_DP
      hfx_int_total_points_tb1 = 0.0_DP
      hfx_int_total_large = 0.0_DP
      hfx_int_total_margin = 0.0_DP
      hfx_int_total_coarse_from = 0.0_DP
      hfx_int_total_coarse_to = 0.0_DP
      hfx_int_total_zero = 0.0_DP

    end subroutine internal_init_stats

    subroutine internal_print_stats()

      ! jd: Gather and show some statistics in verbose mode

      use comms,  only: comms_reduce, pub_on_root
      use constants, only: DP, VERBOSE, stdout
      use rundat, only: pub_output_detail, pub_hfx_integration_variant

      if (pub_output_detail /= VERBOSE) return

      call comms_reduce('SUM',hfx_int_total_calls)
      call comms_reduce('SUM',hfx_int_total_cg_instances)
      call comms_reduce('SUM',hfx_int_total_no_cg_instances)
      call comms_reduce('SUM',hfx_int_total_packed)
      call comms_reduce('SUM',hfx_int_total_points_tb1)
      call comms_reduce('SUM',hfx_int_total_large)
      call comms_reduce('SUM',hfx_int_total_margin)
      call comms_reduce('SUM',hfx_int_total_coarse_from)
      call comms_reduce('SUM',hfx_int_total_coarse_to)
      call comms_reduce('SUM',hfx_int_total_zero)

      if (pub_on_root .and. hfx_int_total_calls /= 0) then
         write(stdout,'(a,i1,a)') 'hfx: Variant ', &
              pub_hfx_integration_variant,' was used for integration.'
         if (pub_hfx_integration_variant > 1) then
            write(stdout,'(a,i7)') 'hfx: Total calls: ', hfx_int_total_calls
            ! avoid /0, if there were no calls to the routine
            if(hfx_int_total_calls == 0) hfx_int_total_calls = 1
            write(stdout,'(a,f5.1,a)') 'hfx: Average filling of tightbox 1: ',&
                 100.0_DP * hfx_int_total_points_tb1 / &
                 real(hfx_int_total_calls,kind=DP),' %'
         end if
         if (pub_hfx_integration_variant == 2) then
            write(stdout,'(a,f5.1,a)') 'hfx: Average filling of tightbox 2: ',&
                 100.0_DP * hfx_int_total_packed / &
                 real(hfx_int_total_calls,kind=DP),' %'
         end if
         if (pub_hfx_integration_variant == 3) then
            write(stdout,'(a,f5.1,a)') 'hfx: Average filling of tightbox 2 &
                 &after coarse-graining: ', 100.0_DP * hfx_int_total_packed / &
                 real(hfx_int_total_calls,kind=DP),' %'
            write(stdout,'(a,f5.1,a)') 'hfx: Incidence of coarse-graining &
                 &attempts: ', 100.0_DP* &
                 real(hfx_int_total_cg_instances,kind=DP) / &
                 max(1.0_DP,real(hfx_int_total_cg_instances + &
                 hfx_int_total_no_cg_instances,kind=DP)),' %'
            write(stdout,'(a,f5.1,a)') 'hfx: Values in tb2 that were zero: ', &
                 100.0_DP* &
                 hfx_int_total_zero/real(hfx_int_total_calls,kind=DP),' %'
            write(stdout,'(a,f5.1,a,f5.1,a)') 'hfx: Values in tb2 that were &
                 &coarse grained: ', 100.0_DP* &
                 hfx_int_total_coarse_from/real(hfx_int_total_calls,kind=DP), &
                 ' % -> ', 100.0_DP* &
                 hfx_int_total_coarse_to/real(hfx_int_total_calls,kind=DP),' %'
            write(stdout,'(a,f5.1,a,f5.1,a)') 'hfx: Values in tb2 that &
                 &were not coarse-grained: ', 100.0_DP* &
                 hfx_int_total_large/real(hfx_int_total_calls,kind=DP), &
                 ' % (on purpose), ', 100.0_DP* &
                 hfx_int_total_margin/real(hfx_int_total_calls,kind=DP), &
                 ' % (due to being on the margin)'

         end if
      end if

    end subroutine internal_print_stats

  end subroutine hf_exchange_calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_atomblock(xatomblock, atoma, atomb, & !output,2*input
       denskern, denskern_idx, overlap_idx, batch_size, & !input
       ngwfs_on_grid, ngwf_basis, elements, calcgradient,& !input
       bb_fftbox, cc_fftbox_sum_batch, & !ws
       work1_fftbox, work2_fftbox, & !ws
       prod_in_fftbox, ket_in_fftbox,& !ws
       ngwf_on_grid_buffer,comp_fftbox,&!workspace
       k_blk, tc_k_blk, dkrecvbuf, dklinrecvbuf, & !workspace
       dd_send_plan, dd_recv_plan, cc_send_plan, cc_recv_plan,& !plan
       atomb_on_node, dds_on_node, atomds_on_node, dds_in_batch,&
       tb_batch, tb_batch2, tb_zwork, swa_tightbox, & !opt sw arg
       sphbessels, full_vmatrix, swcoeff, swswoverlap, & !opt sw args
       tc_denskern, cov_grad, contra_grad, & !Opt grad arguments
       precond_func_recip) !Opt grad arguments

    !==========================================================================!
    ! This subroutine calculates the BA atom block of the Hartree-Fock         !
    ! exchange matrix.  The expression for this is:                            !
    ! V_{B,b,A,a} = \int \sum_{\substack{D \\ S_{AD} \neq 0}}                  !
    ! \sum_d \phi_{A,a}^*(1)\phi_{D,d}(1)                                      !
    ! \left[\sum_{\substack{C \\ S_{BC} \neq 0 \\K^{CD} \neq 0}}               !
    ! \int\frac{\phi_{B,b}(2)\sum_c\phi_{C,c}^*(2)K^{C,cD,d}}                  !
    ! {\vert \vect{r}_1 - \vect{r}_2 \vert} d\vect{r}_2 \right] d\vect{r}_1    !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in December 2008 and January 2009.               !
    !==========================================================================!

    use comms, only: comms_abort, pub_total_num_nodes, comms_reduce,&
         pub_my_node_id !, pub_on_root
    use constants, only: DP, stderr
    use function_basis, only: FUNC_BASIS
    use geometry, only: point
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_max_atoms_on_node, &
         pub_first_atom_on_node
    use simulation_cell, only: pub_fftbox, pub_cell, pub_maxtight_pts1, &
         pub_maxtight_pts2, pub_maxtight_pts3
    use sparse, only: SPAM3
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    integer,        intent(in) :: atoma ! Global index of atom A
    integer,        intent(in) :: atomb ! Global index of atom B
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(out) :: xatomblock(ngwf_basis%max_on_atom,&
         ngwf_basis%max_on_atom,pub_cell%num_spins)
    type(SPAM3),    intent(in) :: denskern(pub_cell%num_spins)
    integer,        intent(in) :: denskern_idx(:)
    integer,        intent(in) :: overlap_idx(:)
    integer,        intent(in) :: batch_size
    real(kind=DP),  intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(ELEMENT),  intent(in) :: elements(pub_cell%nat)
    logical,        intent(in) :: calcgradient

    ! qoh: Spherical Bessel functions
    type(bessel),optional,intent(in) :: sphbessels(pub_cell%num_species,&
         max_num_bessels)
    type(SPAM3), optional,intent(in) :: full_vmatrix

    ! Arguments for calculating NGWF gradient (required if calcgradient is true)
    type(SPAM3),  optional,intent(in)    :: tc_denskern(pub_cell%num_spins)! KS
    real(kind=DP),optional,intent(inout) :: cov_grad(ngwf_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP),optional,intent(inout) :: contra_grad&
         (ngwf_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP),optional,intent(in)    :: precond_func_recip(:,:,:)

    ! qoh: Workspace

    real(kind=DP), intent(out) :: bb_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2, pub_fftbox%total_pt3)
    real(kind=DP),intent(out) :: cc_fftbox_sum_batch(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3,pub_cell%num_spins,batch_size)
    real(kind=DP), intent(out)  :: work1_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(out)  :: work2_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(out) :: prod_in_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2, pub_fftbox%total_pt3,2)
    real(kind=DP), intent(out):: ket_in_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3, pub_cell%num_spins)
    real(kind=DP), intent(out) :: ngwf_on_grid_buffer&
         (ngwf_basis%max_n_ppds_sphere*pub_cell%n_pts)
    real(kind=DP), intent(out) :: k_blk(ngwf_basis%max_on_atom,&
         ngwf_basis%max_on_atom,pub_cell%num_spins)
    real(kind=DP), intent(out) :: tc_k_blk(ngwf_basis%max_on_atom,&
         ngwf_basis%max_on_atom,pub_cell%num_spins)
    real(kind=DP), intent(out) :: dkrecvbuf(batch_size, ngwf_basis%max_on_atom, &
         pub_cell%num_spins)
    real(kind=DP), intent(out) :: dklinrecvbuf&
         (batch_size*ngwf_basis%max_on_atom*pub_cell%num_spins)
    complex(kind=DP), intent(out) :: comp_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2, pub_fftbox%total_pt3)

    !qoh: Spherical wave arrays only required if using spherical waves
    real(kind=DP), optional, intent(out) :: swcoeff&
         (max_sw_set_size,batch_size,pub_cell%num_spins)
    real(kind=DP), optional, intent(out) :: swswoverlap(max_sw_set_size)
    real(kind=DP), optional, intent(out) :: swa_tightbox&
         (pub_maxtight_pts1,pub_maxtight_pts2,pub_maxtight_pts3)
    real(kind=DP), optional, intent(out) :: tb_batch&
         (pub_maxtight_pts1,pub_maxtight_pts2,pub_maxtight_pts3, &
         batch_size, pub_cell%num_spins)
    real(kind=DP), optional, intent(out) :: tb_batch2&
         (pub_maxtight_pts1,pub_maxtight_pts2,pub_maxtight_pts3, &
         batch_size, pub_cell%num_spins)
    complex(kind=DP), optional, intent(out) :: tb_zwork&
         (pub_maxtight_pts1,pub_maxtight_pts2,pub_maxtight_pts3)

    !qoh: Plan arrays
    logical, intent(out) :: dd_send_plan(ngwf_basis%max_on_node,&
         0:pub_total_num_nodes-1)
    logical, intent(out) :: dd_recv_plan(ngwf_basis%max_on_node,&
         0:pub_total_num_nodes-1)
    logical, intent(out) :: cc_send_plan(pub_max_atoms_on_node,&
         0:pub_total_num_nodes-1)
    logical, intent(out) :: cc_recv_plan(pub_max_atoms_on_node,&
         0:pub_total_num_nodes-1)
    integer, intent(out) :: atomb_on_node(0:pub_total_num_nodes-1)
    integer, intent(out) :: dds_on_node(batch_size,0:pub_total_num_nodes-1)
    integer, intent(out) :: atomds_on_node(batch_size,0:pub_total_num_nodes-1)
    integer, intent(out) :: dds_in_batch(batch_size)

    !qoh: Local Variables
    integer :: atomd ! Global index of atom D
    integer :: batch_count !Number of items in batch
    integer :: bbngwfidx  ! Global index of current B,b NGWF
    integer :: ngwf_atomb, ngwf_atomd ! Index of ngwf on atomx
    logical :: bb_interpolated ! Has the current B,b NGWF been interpolated
    integer :: current_batch_size
    integer :: num_dds ! Number of D,d NGWFs overlapping with A
    integer :: n_batches ! Number of batches
    integer :: current_dd ! Current D,d in the batch [1,num_dds]
    integer :: batch_next_dd ! Next D,d in the batch [1,num_dds+1]
    integer :: batch ! Batch index
    integer :: local_a ! Local index of atom A
    logical :: all_finished ! Are all nodes finished?
    logical :: swswovlpdone ! Has SW-SW overlap been calculated
    logical :: boxiscell ! FFT box == cell
    logical :: use_fft ! Use FFT method

    !qoh: Column of atoma in overlap matrix
    integer :: o_fstnzbridx_cola ! first overlap non-zero block row index
    integer :: o_lstnzbridx_cola !last overlap non-zero block row index
    integer :: o_nzbridx_cola !overlap non-zero block row index
    type(atom_centre) :: centrea
    type(atom_centre) :: centreb
    integer :: current_sw_set_size
    real(kind=DP), allocatable :: vmatrix(:,:)
    integer :: ierr
    logical :: ab_overlap

    integer :: only_c, cmax !@@
    integer :: only_cf, cfmax !@@

    logical, parameter :: hack = .false. !@@ also sendngwfs turn off, also 423

    call utils_trace_in('hf_exchange_atomblock')
    call timer_clock('hf_exchange_atomblock',1)

    local_a = atoma - pub_first_atom_on_node(pub_my_node_id) + 1

    o_fstnzbridx_cola = overlap_idx(local_a)
    o_lstnzbridx_cola = overlap_idx(local_a+1) - 1

    !qoh: Does B overlap with atom A?
    ab_overlap = any(overlap_idx(o_fstnzbridx_cola:o_lstnzbridx_cola) == &
         atomb)
    boxiscell = (pub_fftbox%coin1 .and. pub_fftbox%coin2 .and. pub_fftbox%coin3)
    use_fft = ( (boxiscell .and. allowfftonly) &
         .or. (atoma == atomb .and. usefftforsc) &
         .or. (ab_overlap .and. overlappinginfft) )

    if (.not. use_fft) then
       call hf_exchange_init_centre(centreb, ngwf_basis, elements, atomb)
       call hf_exchange_init_centre(centrea, ngwf_basis, elements, atoma)
       if (.not. useoverlapmetric .and. .not. usenpa) then
          current_sw_set_size = hf_exchange_sw_set_size(sphbessels,centrea,&
               centreb)
          allocate(vmatrix(current_sw_set_size,current_sw_set_size),stat=ierr)
          call utils_alloc_check('hf_exchange_atomblock','vmatrix',ierr)
          call hf_exchange_make_vmatrix(vmatrix, full_vmatrix, &
               current_sw_set_size, sphbessels, centrea, centreb, &
               atoma, atomb)
       else
          current_sw_set_size = 0
       end if
       swcoeff = 0.0_DP
       swswoverlap = 0.0_DP
    end if

    swswovlpdone = .false.
    xatomblock = 0.0_DP

    !qoh: Count the number of NGWFs that overlap with atom A
    !qoh: Loop over atoms (atom D) that overlap with atom A
    num_dds = 0
    countdds: do o_nzbridx_cola = o_fstnzbridx_cola, o_lstnzbridx_cola
       atomd = overlap_idx(o_nzbridx_cola)
       num_dds = num_dds + ngwf_basis%num_on_atom(atomd)
    end do countdds

    !qoh: Calculate number of batches
    n_batches = num_dds / batch_size

    if (mod(num_dds, batch_size) > 0) n_batches = n_batches + 1

    loop_bf: do ngwf_atomb =1,ngwf_basis%num_on_atom(atomb)
       bbngwfidx =  ngwf_basis%first_on_atom(atomb) + ngwf_atomb - 1

       !qoh: initialise variables.
       ket_in_fftbox= 0.0_DP
       bb_interpolated=.false.

       batch_next_dd = 1
       !qoh: Loop over batches
       batches: do batch=1,n_batches

          !qoh: Tell other nodes that we are still not finished.
          all_finished = .false.
          call timer_clock('hf_exchange_comms_wait',1)
          call comms_reduce("AND",all_finished)
          call timer_clock('hf_exchange_comms_wait',2)

          !qoh: Decide which Dd NGWFs are in this batch.
          current_dd = 0
          batch_count = 0
          dds_in_batch = 0
          loop_D: do o_nzbridx_cola = o_fstnzbridx_cola, o_lstnzbridx_cola
             atomd = overlap_idx(o_nzbridx_cola)
             loop_df: do ngwf_atomd = 1,ngwf_basis%num_on_atom(atomd)
                current_dd = current_dd + 1
                ! qoh: Only if include if the Dd NGWF is not already in a batch
                if (current_dd == batch_next_dd) then
                   batch_count = batch_count + 1
                   dds_in_batch(batch_count) =  &
                        ngwf_basis%first_on_atom(atomd) + ngwf_atomd - 1
                   batch_next_dd = current_dd + 1
                   if (batch_count == batch_size) exit
                end if
             end do loop_df
             if (batch_count == batch_size) exit
          end do loop_D
          current_batch_size = batch_count

          call timer_clock('hf_exchange_ccsum_batch',1)
          call hf_exchange_ccsum_batch(cc_fftbox_sum_batch, bb_fftbox, &
               k_blk, denskern, denskern_idx, overlap_idx, ngwfs_on_grid, &
               ngwf_basis, bbngwfidx, dds_in_batch, &
               batch_size, current_batch_size, bb_interpolated, &
               dkrecvbuf,dklinrecvbuf, & !comms buffers
               ngwf_on_grid_buffer, &!comms buffers
               dd_send_plan, dd_recv_plan, cc_send_plan, cc_recv_plan,& !plan
               atomb_on_node, dds_on_node, atomds_on_node) !for plan
              
          call timer_clock('hf_exchange_ccsum_batch',2)

          if (use_fft) then
             call timer_clock('hf_exchange_fftpot_ngwf_ket',1)
             call hf_exchange_fftpot_ngwf_ket(ket_in_fftbox,&
                  prod_in_fftbox, bb_fftbox, cc_fftbox_sum_batch,&
                  work1_fftbox, work2_fftbox, &
                  ngwf_basis, ngwfs_on_grid, &
                  batch_size,  current_batch_size, &
                  ngwf_on_grid_buffer, &!commsbufs
                  dds_in_batch, dd_send_plan, dd_recv_plan, &
                  comp_fftbox,atoma,atomb)
             call timer_clock('hf_exchange_fftpot_ngwf_ket',2)
          else
             call hf_exchange_ngwf_product(tb_batch2, prod_in_fftbox, &
                  bb_fftbox, cc_fftbox_sum_batch, centreb, batch_size, &
                  current_batch_size)
             if (.not. usenpa) &
                  call hf_exchange_expansion(sphbessels, centrea, centreb, &
                  swswoverlap, swcoeff, batch_size, current_batch_size, &
                  swswovlpdone, swa_tightbox, tb_batch2, tb_zwork, &
                  current_sw_set_size, vmatrix)

             call timer_clock('hf_exchange_tbpot_ngwf_ket',1)
             call hf_exchange_tbpot_ngwf_ket(ket_in_fftbox, &
                  prod_in_fftbox, tb_batch, tb_batch2,&
                  work1_fftbox, work2_fftbox, &
                  ngwf_basis, ngwfs_on_grid, batch_size, current_batch_size,&
                  ngwf_on_grid_buffer, &!commsbufs
                  dds_in_batch, dd_send_plan, dd_recv_plan, atoma, atomb,&
                  centrea, centreb, sphbessels, swcoeff)
             call timer_clock('hf_exchange_tbpot_ngwf_ket',2)

           end if

       end do batches

       !qoh: Calculate matrix elements in row b of atom block
       call timer_clock('hf_exchange_atomblock_elements',1)
       call hf_exchange_atomblock_elements(xatomblock(ngwf_atomb,:,:),&
            ket_in_fftbox, ngwfs_on_grid, atoma, ngwf_basis, &
            ngwf_on_grid_buffer)
       call timer_clock('hf_exchange_atomblock_elements',2)

       !qoh: Calculate gradient if required
       if(calcgradient) then
          call timer_clock('hf_exchange_gradient',1)
          call hf_exchange_gradient(contra_grad, cov_grad, & !output
               ket_in_fftbox, atoma, atomb, ngwf_atomb,ngwf_basis,&!input
               denskern, tc_denskern, & !input
               precond_func_recip,&!input
               ngwf_on_grid_buffer, k_blk, tc_k_blk, work1_fftbox) !workspace
          call timer_clock('hf_exchange_gradient',2)
       end if

    end do loop_bf

    if (printxmat .and. justoneblock) print *, xatomblock

    call timer_clock('hf_exchange_atomblock',2)
    if (justoneblock) call timer_clock('total_time',3)
    if (justoneblock) call comms_abort

    if (allocated(vmatrix)) then
       deallocate(vmatrix,stat=ierr)
       call utils_dealloc_check('hf_exchange_atomblock','vmatrix',ierr)
    end if

    call utils_trace_out('hf_exchange_atomblock')

  end subroutine hf_exchange_atomblock

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_sendngwfs(denskern,denskern_idx,overlap_idx,&
       batch_size,ngwfs_on_grid,ngwf_basis, &!in
       bb_fftbox, cc_fftbox_sum_batch, & !ws
       work1_fftbox, work2_fftbox, & !ws
       prod_in_fftbox, ket_in_fftbox,  ngwf_on_grid_buffer, & !ws
       comp_fftbox,k_blk,dkrecvbuf,dklinrecvbuf, & !ws
       dd_send_plan, dd_recv_plan, cc_send_plan, cc_recv_plan,& !plan
       atomb_on_node, dds_on_node, & !Used by plan
       atomds_on_node, dds_in_batch) !Used by plan.

    !==========================================================================!
    ! This subroutine exists to send the required NGWFs and density kernel     !
    ! elements to other nodes that are still in hf_exchange_atomblock after    !
    ! this node has finished.                                                  !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in January 2009.                                 !
    !==========================================================================!

    use comms, only: pub_total_num_nodes, comms_reduce
    use constants, only: DP
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: pub_max_atoms_on_node
    use simulation_cell, only: pub_fftbox, pub_cell
    use sparse, only: SPAM3
    use timer, only: timer_clock
    use utils, only: utils_abort ! @removeme

    implicit none

    !qoh: Arguments
    type(SPAM3),   intent(in) :: denskern(pub_cell%num_spins)
    integer,       intent(in) :: denskern_idx(:)
    integer,       intent(in) :: overlap_idx(:)
    integer,       intent(in) :: batch_size
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)

    !qoh: Workspace
    real(kind=DP), intent(out) :: bb_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2, pub_fftbox%total_pt3)
    real(kind=DP),intent(out) :: cc_fftbox_sum_batch(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3,pub_cell%num_spins,batch_size)
    real(kind=DP), intent(out)  :: work1_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(out)  :: work2_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(out) :: prod_in_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2, pub_fftbox%total_pt3,2)
    real(kind=DP), intent(out) :: ket_in_fftbox(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2, pub_fftbox%total_pt3,pub_cell%num_spins)
    real(kind=DP), intent(out) :: ngwf_on_grid_buffer&
         (ngwf_basis%max_n_ppds_sphere*pub_cell%n_pts)
    real(kind=DP), intent(out) :: k_blk(ngwf_basis%max_on_atom,&
         ngwf_basis%max_on_atom,pub_cell%num_spins)
    real(kind=DP), intent(out) :: dkrecvbuf(batch_size, ngwf_basis%max_on_atom, &
         pub_cell%num_spins)
    real(kind=DP), intent(out) :: dklinrecvbuf&
         (batch_size*ngwf_basis%max_on_atom*pub_cell%num_spins)
    complex(kind=DP), intent(out) :: comp_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2, pub_fftbox%total_pt3)

    !qoh: Plan arrays
    logical, intent(out) :: dd_send_plan(ngwf_basis%max_on_node,&
         0:pub_total_num_nodes-1)
    logical, intent(out) :: dd_recv_plan(ngwf_basis%max_on_node,&
         0:pub_total_num_nodes-1)
    logical, intent(out) :: cc_send_plan(pub_max_atoms_on_node,&
         0:pub_total_num_nodes-1)
    logical, intent(out) :: cc_recv_plan(pub_max_atoms_on_node,&
         0:pub_total_num_nodes-1)
    integer, intent(out) :: atomb_on_node(0:pub_total_num_nodes-1)
    integer, intent(out) :: dds_on_node(batch_size,0:pub_total_num_nodes-1)
    integer, intent(out) :: atomds_on_node(batch_size,0:pub_total_num_nodes-1)
    integer, intent(out) :: dds_in_batch(batch_size)

    !qoh: Local variables and parameters
    logical :: all_finished
    logical :: bb_interpolated

    integer, parameter :: current_batch_size  = 0 ! Null batch
    integer, parameter :: bbngwfidx = 0           ! Null NGWF B,b
    integer, parameter :: atoma = 0               ! Null atom A

    call utils_trace_in('hf_exchange_sendngwfs')

    do
       all_finished = .true. ! We are finished
       !qoh: Are any other nodes still busy calculating atomblocks?
       call timer_clock('hf_exchange_comms_wait',1)
       call comms_reduce("AND",all_finished)
       call timer_clock('hf_exchange_comms_wait',2)
       if (all_finished) exit
       bb_interpolated = .true. ! It's not but we don't want to interpolate it.
       dds_in_batch = 0 ! Nothing in our batch!

       !qoh: Send B,b and C,c NGWFs and prepare plan
       call timer_clock('hf_exchange_ccsum_batch',1)

       call hf_exchange_ccsum_batch(cc_fftbox_sum_batch, bb_fftbox, &
            k_blk, denskern, denskern_idx, overlap_idx, ngwfs_on_grid, &
            ngwf_basis, bbngwfidx, dds_in_batch, &
            batch_size, current_batch_size, bb_interpolated, &
            dkrecvbuf, dklinrecvbuf, & ! comms buffers
            ngwf_on_grid_buffer, & ! comms buffers
            dd_send_plan, dd_recv_plan, cc_send_plan, cc_recv_plan,& !plan
            atomb_on_node, dds_on_node, atomds_on_node)!for plan
       call timer_clock('hf_exchange_ccsum_batch',2)

       !qoh: Send D,d NGWFs
       call timer_clock('hf_exchange_fftpot_ngwf_ket',1)
       call hf_exchange_fftpot_ngwf_ket(ket_in_fftbox,&
            prod_in_fftbox, bb_fftbox, cc_fftbox_sum_batch,&
            work1_fftbox, work2_fftbox, &
            ngwf_basis, ngwfs_on_grid, &
            batch_size,  current_batch_size, &
            ngwf_on_grid_buffer, &!commsbufs
            dds_in_batch, dd_send_plan, dd_recv_plan, &
            comp_fftbox,atoma,0)
       call timer_clock('hf_exchange_fftpot_ngwf_ket',2)
    end do

    call utils_trace_out('hf_exchange_sendngwfs')

  end subroutine hf_exchange_sendngwfs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_ccsum_batch(cc_fftbox_sum_batch, bb_fftbox, &
       dck_blk, denskern, denskern_idx, overlap_idx, ngwfs_on_grid, &
       ngwf_basis, bbngwfidx, dds_in_batch, &
       batch_size, current_batch_size, bb_interpolated, &
       dkrecvbuf,dklinrecvbuf,cc_on_grid_buffer,&!comms buffers
       dd_send_plan, dd_recv_plan, cc_send_plan, cc_recv_plan,& !plan
       atomb_on_node, dds_on_node, atomds_on_node) !for plan

    !==========================================================================!
    ! This subroutine calculates:                                              !
    ! \sum_{\substack{C \\ S_{BC} \neq 0 \\K^{CD} \neq 0}}                     !
    ! \sum_c\phi_{C,c}^*(2)K^{C,cD,d}}                                         !
    ! A communications plan is produced for the C,c and D,d NGWFs by calling   !
    ! hf_exchange_construct_plan.                                              !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2008 and January 2009.                        !
    ! Modified by Quintin Hill on 11/02/2009 to use linear DK recv buffer.     !
    !==========================================================================!

    use basis, only: basis_add_function_to_box,basis_ket_start_wrt_fftbox,&
         basis_copy_function_to_box, basis_location_func_wrt_cell, &
         basis_location_fb_wrt_box
    use comms, only: pub_total_num_nodes,comms_send, comms_recv, &
         pub_my_node_id, comms_free
    use constants, only: DP
    use function_basis, only: FUNC_BASIS, pub_buffer_sphere
    use parallel_strategy, only: pub_max_atoms_on_node, &
         pub_first_atom_on_node, pub_num_atoms_on_node
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: SPAM3, sparse_get_block, sparse_get_element
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    type buffptr
       real(kind=DP), pointer :: dk(:,:,:) ! Denskern send buffer
    end type buffptr

    integer,        intent(in) :: batch_size
    type(FUNC_BASIS),   intent(in) :: ngwf_basis
    real(kind=DP), intent(out) :: cc_fftbox_sum_batch(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3,pub_cell%num_spins,batch_size)
    real(kind=DP),intent(inout):: bb_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(out) :: dck_blk(ngwf_basis%max_on_atom,&
         ngwf_basis%max_on_atom,pub_cell%num_spins)
    real(kind=DP), intent(out) :: dkrecvbuf(batch_size,ngwf_basis%max_on_atom, &
         pub_cell%num_spins)
    real(kind=DP), intent(out) :: dklinrecvbuf&
         (batch_size*ngwf_basis%max_on_atom*pub_cell%num_spins)
    type(SPAM3),    intent(in) :: denskern(pub_cell%num_spins)
    integer,        intent(in) :: denskern_idx(:)
    integer,        intent(in) :: overlap_idx(:)
    real(kind=DP),  intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    integer,        intent(in) :: bbngwfidx ! Global index of NGWF B,b
    integer,        intent(in) :: dds_in_batch(batch_size)
    integer,        intent(in) :: current_batch_size
    logical,     intent(inout) :: bb_interpolated !Has NGWF B,b been interpolated
    real(kind=DP), intent(out) :: cc_on_grid_buffer&
         (ngwf_basis%max_n_ppds_sphere*pub_cell%n_pts)

    ! Plan arrays
    logical, intent(out) :: dd_send_plan(ngwf_basis%max_on_node,&
         0:pub_total_num_nodes-1)
    logical, intent(out) :: dd_recv_plan(ngwf_basis%max_on_node,&
         0:pub_total_num_nodes-1)
    logical, intent(out) :: cc_send_plan(pub_max_atoms_on_node,&
         0:pub_total_num_nodes-1)
    logical, intent(out) :: cc_recv_plan(pub_max_atoms_on_node,&
         0:pub_total_num_nodes-1)
    integer, intent(out) :: atomb_on_node(0:pub_total_num_nodes-1)
    integer, intent(out) :: dds_on_node(batch_size,0:pub_total_num_nodes-1)
    integer, intent(out) :: atomds_on_node(batch_size,0:pub_total_num_nodes-1)

    ! Local variables:
    type(buffptr), allocatable :: dkbuffer(:,:) ! Send buffer of density kernel
    real(kind=DP) :: dk_ccddel  ! K_{C,cD,d}
    integer :: first_ngwf_idx_of_c !Global NGWF index of the first NGWF on atom C
    integer :: atomc ! Global index of atom C
    integer :: ngwf_atomc ! NGWF on atom C
    integer :: bb_start1, bb_start2, bb_start3
    integer :: bb_cell_start1, bb_cell_start2, bb_cell_start3
    integer :: cc_cell_start1, cc_cell_start2, cc_cell_start3
    integer :: cc_start1, cc_start2, cc_start3
    integer :: is  ! Spin counter
    integer :: batch_count
    integer :: atomd ! Global index of atom D
    integer :: ngwf_atomd ! NGWF on atom D
    integer :: local_c ! Local index of atom C
    integer :: local_cc ! Local index of Cc NGWF
    integer :: first_local_cc ! First NGWF on atom C
    integer :: last_local_cc ! Last NGWF on atom C
    integer :: local_cc_start ! Starting index in ngwfs_on_grid of Cc
    integer :: send_cc ! MPI send tag
    integer :: remote_c ! Remote index of atom C
    integer :: recv_cc ! MPI recv tag
    integer :: n_cc_ppds ! Number of Cc ppds
    integer :: cc_npts ! Number of points in Cc
    integer :: previous_d ! Previous atom D
    integer :: node
    integer :: nodeshift
    integer :: ierr ! Error flag
    integer :: nels ! number of elements in array

    call utils_trace_in('hf_exchange_ccsum_batch')

    local_cc = 0

    !qoh: Initialise dkbuffer
    allocate(dkbuffer(pub_num_atoms_on_node(pub_my_node_id),&
         0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('hf_exchange_ccsum_batch','dkbuffer',ierr)
    do node=0,pub_total_num_nodes-1
       do local_c=1,pub_num_atoms_on_node(pub_my_node_id)
          nullify(dkbuffer(local_c,node)%dk)
       end do
    end do

    cc_fftbox_sum_batch = 0.0_DP

    call hf_exchange_construct_plan(dd_send_plan, dd_recv_plan, &
         cc_send_plan, cc_recv_plan, bbngwfidx, batch_size, dds_in_batch,&
         denskern_idx, overlap_idx, atomb_on_node, dds_on_node,&
         atomds_on_node, ngwf_basis)

    ! cks: the col function will always be at the centre of the fftbox
    ! cks:or  left where it is in the simulation cell
    ! cks: (when the simulation cell and FFT-box coincide)
    call basis_ket_start_wrt_fftbox(bb_start1,bb_start2,bb_start3,&
         pub_fftbox%total_pt1,pub_fftbox%total_pt2,pub_fftbox%total_pt3)

    if (bbngwfidx /=0) call basis_location_func_wrt_cell(bb_cell_start1, &
         bb_cell_start2, bb_cell_start3, ngwf_basis%all_tbs(bbngwfidx))

    !qoh: Send relevant NGWFs and density kernel elements
    send_node: do nodeshift=1, pub_total_num_nodes-1
       node = modulo(pub_my_node_id+nodeshift,pub_total_num_nodes)
       loc_c: do local_c=1,pub_num_atoms_on_node(pub_my_node_id)
          insplan: if (cc_send_plan(local_c,node)) then
             atomc = pub_first_atom_on_node(pub_my_node_id) + local_c - 1
             allocate(dkbuffer(local_c,node)%dk(batch_size,&
                  ngwf_basis%num_on_atom(atomc),pub_cell%num_spins),stat=ierr)
             call utils_alloc_check('hf_exchange_ccsum_batch','dkbuffer%dk',ierr)
             dkbuffer(local_c,node)%dk = 0.0_DP
             dck_blk = 0.0_DP
             previous_d = 0
             rem_dd: do batch_count = 1, batch_size
                atomd = atomds_on_node(batch_count,node)
                if (atomd == 0) exit
                if (atomd /= previous_d) then
                   do is=1,pub_cell%num_spins
                      call sparse_get_block(dck_blk(:,:,is),denskern(is),&
                           atomd,atomc)
                   end do
                   previous_d = atomd
                end if
                ngwf_atomd = dds_on_node(batch_count,node) - &
                     ngwf_basis%first_on_atom(atomd) +1
                dkbuffer(local_c,node)%dk(batch_count,:,:) = &
                     dck_blk(ngwf_atomd,1:ngwf_basis%num_on_atom(atomc),:)
             end do rem_dd

             call comms_send(node,dkbuffer(local_c,node)%dk, tag=atomc)

             first_local_cc = ngwf_basis%first_on_atom(atomc)&
                  - ngwf_basis%first_on_node(pub_my_node_id)+1
             last_local_cc = first_local_cc + ngwf_basis%num_on_atom(atomc)-1

             !qoh: Pack into buffer
             n_cc_ppds = ngwf_basis%spheres(first_local_cc)%n_ppds_sphere
             cc_npts = n_cc_ppds * pub_cell%n_pts
             sendcc: do local_cc=first_local_cc, last_local_cc
                send_cc = local_cc + ngwf_basis%first_on_node(pub_my_node_id) - 1
                local_cc_start = ngwf_basis%spheres(local_cc)%offset
                ! ndmh: send list of ppds straight from sphere
                call comms_send(node, ngwf_basis%spheres(local_cc)%ppd_list, &
                     2*ngwf_basis%n_ppds_sphere(send_cc), tag=send_cc)
                !qoh: send list of values in points in ppds
                call comms_send(node, ngwfs_on_grid(local_cc_start), &
                     cc_npts, tag=send_cc)
             end do sendcc

          end if insplan
       end do loc_c
    end do send_node

    !qoh: Receive relevant NGWFs
    recv_node: do nodeshift=1, pub_total_num_nodes
       node = modulo(pub_my_node_id+nodeshift,pub_total_num_nodes)
       rem_c: do remote_c=1,pub_num_atoms_on_node(node)
          inrplan: if (cc_recv_plan(remote_c,node)) then

             atomc = remote_c + pub_first_atom_on_node(node) - 1
             first_ngwf_idx_of_c = ngwf_basis%first_on_atom(atomc)

             !qoh: Find position of \phi_{C,c} wrt simulation cell
             !qoh: This is the same for all NGWFs on atom C
             call basis_location_func_wrt_cell(cc_cell_start1, &
                  cc_cell_start2,cc_cell_start3, &
                  ngwf_basis%all_tbs(first_ngwf_idx_of_c))

             !qoh: Find where to deposit \phi_{C,c} function in fftbox
             !qoh: This is the same for all NGWFs on atom C
             call basis_location_fb_wrt_box( &
                  cc_start1, cc_start2, cc_start3, &
                  bb_start1, bb_start2, bb_start3, &
                  bb_cell_start1, bb_cell_start2, bb_cell_start3, &
                  cc_cell_start1, cc_cell_start2, cc_cell_start3, &
                  pub_cell%total_pt1, pub_cell%total_pt2, pub_cell%total_pt3)

             dklinrecvbuf = 0.0_DP
             dkrecvbuf = 0.0_DP
             if ( node /= pub_my_node_id ) then
                call comms_recv(node,dklinrecvbuf,tag=atomc)
                nels = batch_size*ngwf_basis%num_on_atom(atomc)*pub_cell%num_spins
                dkrecvbuf(1:batch_size,1:ngwf_basis%num_on_atom(atomc),&
                     1:pub_cell%num_spins)=reshape(dklinrecvbuf(1:nels),&
                     (/batch_size,ngwf_basis%num_on_atom(atomc),&
                     pub_cell%num_spins/))
             end if

             loop_cf: do ngwf_atomc = 1, ngwf_basis%num_on_atom(atomc)
                recv_cc = first_ngwf_idx_of_c + ngwf_atomc - 1
                cc_on_grid_buffer = 0.0_DP

                if ( node == pub_my_node_id )  then
                   local_cc = recv_cc - ngwf_basis%first_on_node(node) + 1
                else
                   ! ndmh: receive straight into pub_buffer_sphere
                   call comms_recv(node,pub_buffer_sphere%ppd_list, &
                        2*ngwf_basis%n_ppds_sphere(recv_cc), tag=recv_cc)
                   n_cc_ppds = ngwf_basis%n_ppds_sphere(recv_cc)
                   pub_buffer_sphere%n_ppds_sphere = n_cc_ppds
                   pub_buffer_sphere%offset = 1
                   !qoh: Receive NGWFs from node
                   cc_npts = n_cc_ppds * pub_cell%n_pts
                   call comms_recv(node, cc_on_grid_buffer, &
                        cc_npts, tag=recv_cc)
                end if

                batch: do batch_count = 1, current_batch_size

                   !qoh: Put \phi_{C,c}* K^{C,cD,d} in batch of FFT boxes.

                   ! jd: - cc_fftbox_sum_batch: box where data will be stored
                   ! jd: - ngwf_basis%all_tbs(recv_cc): tightbox of Cc
                   ! jd: - ngwfs_on_grid or cc_on_grid_buffer: \phi_{C,c}
                   ! jd: - spheres...
                   ! jd: - dk_ccddel - K^{C,cD,d}
                   spins: do is=1,pub_cell%num_spins
                      if ( node == pub_my_node_id )  then
                         !qoh: Denskern is symmetric K^{C,cD,d} = K^{D,dC,c}
                         call sparse_get_element(dk_ccddel,denskern(is),&
                              dds_in_batch(batch_count),recv_cc)

                         call basis_add_function_to_box(&
                              cc_fftbox_sum_batch(:,:,:,is,batch_count), &
                              pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                              pub_fftbox%total_pt3, &
                              cc_start1, cc_start2, cc_start3, &
                              ngwf_basis%all_tbs(recv_cc), ngwfs_on_grid,&
                              ngwf_basis%spheres(local_cc), dk_ccddel)
                      else
                         dk_ccddel = dkrecvbuf(batch_count,ngwf_atomc,is)

                         call basis_add_function_to_box(&
                              cc_fftbox_sum_batch(:,:,:,is,batch_count), &
                              pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                              pub_fftbox%total_pt3, &
                              cc_start1, cc_start2, cc_start3, &
                              ngwf_basis%all_tbs(recv_cc), cc_on_grid_buffer,&
                              pub_buffer_sphere, dk_ccddel)
                      end if
                   end do spins
                end do batch

                if (.not. bb_interpolated .and. recv_cc == bbngwfidx) then
                   !qoh: use new basis_copy routine to transfer \phi_{B,b}
                   !qoh: to FFT box if we haven't done so yet
                   bb_fftbox = 0.0_DP
                   if ( node == pub_my_node_id )  then
                      call basis_copy_function_to_box(bb_fftbox, &
                           pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                           pub_fftbox%total_pt3, &
                           bb_start1, bb_start2, bb_start3, &
                           ngwf_basis%all_tbs(bbngwfidx), ngwfs_on_grid, &
                           ngwf_basis%spheres(local_cc))
                   else
                      call basis_copy_function_to_box(bb_fftbox, &
                           pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                           pub_fftbox%total_pt3, &
                           bb_start1, bb_start2, bb_start3, &
                           ngwf_basis%all_tbs(bbngwfidx), cc_on_grid_buffer, &
                           pub_buffer_sphere)
                   end if

                end if

             end do loop_cf
          end if inrplan
       end do rem_c
    end do recv_node

    !qoh: Interpolate \phi_{B,b} if this hasn't been done yet.
    bb_interpolated = .true.
    call timer_clock('hf_exchange_comms_wait',1)
    call comms_free
    call timer_clock('hf_exchange_comms_wait',2)

    !qoh: Destroy the density kernel buffer.
    do local_c=1,pub_num_atoms_on_node(pub_my_node_id)
       do node=0,pub_total_num_nodes-1
          if(associated(dkbuffer(local_c,node)%dk)) then
             deallocate(dkbuffer(local_c,node)%dk, stat=ierr)
             call utils_dealloc_check('hf_exchange_ccsum_batch',&
                  'dkbuffer%dk',ierr)
          end if
       end do
    end do
    deallocate(dkbuffer,stat=ierr)
    call utils_dealloc_check('hf_exchange_ccsum_batch','dkbuffer',ierr)

    call utils_trace_out('hf_exchange_ccsum_batch')

  end subroutine hf_exchange_ccsum_batch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_construct_plan(dd_send_plan, dd_recv_plan, &
       cc_send_plan, cc_recv_plan, bbngwfidx, batch_size, dds_in_batch,&
       denskern_idx, overlap_idx, atomb_on_node, dds_on_node,&
       atomds_on_node,ngwf_basis)

    !==========================================================================!
    ! This subroutine constructs the plan for communication of the B,b NGWF,   !
    ! the C,c NGWFs and the D,d NGWFs. The currently implemented plan has a    !
    ! communications suitable for calculating the exchange potential using     !
    ! FFTs.                                                                    !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in December 2008 and January 2009.               !
    !==========================================================================!

    use comms, only: pub_my_node_id, comms_bcast, comms_send, comms_recv, &
         pub_total_num_nodes, comms_barrier, comms_free
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: pub_max_atoms_on_node, pub_num_atoms_on_node, &
         pub_first_atom_on_node

    implicit none

    type(FUNC_BASIS), intent(in) :: ngwf_basis
    logical, intent(out) :: dd_send_plan(ngwf_basis%max_on_node,&
         0:pub_total_num_nodes-1)
    logical, intent(out) :: dd_recv_plan(ngwf_basis%max_on_node,&
         0:pub_total_num_nodes-1)
    logical, intent(out) :: cc_send_plan(pub_max_atoms_on_node,&
         0:pub_total_num_nodes-1)
    logical, intent(out) :: cc_recv_plan(pub_max_atoms_on_node,&
         0:pub_total_num_nodes-1)
    integer, intent(in)  :: bbngwfidx ! Global index of current B,b
    integer, intent(in)  :: batch_size
    integer, intent(in)  :: dds_in_batch(batch_size)
    integer, intent(in)  :: denskern_idx(:)
    integer, intent(in)  :: overlap_idx(:)
    integer, intent(out) :: atomb_on_node(0:pub_total_num_nodes-1)
    integer, intent(out) :: dds_on_node(batch_size,0:pub_total_num_nodes-1)
    integer, intent(out) :: atomds_on_node(batch_size,0:pub_total_num_nodes-1)

    integer :: o_fstnzbridx_col ! First non-zero block index in overlap matrix
    integer :: o_lstnzbridx_col ! Last non-zero block index in overlap matrix
    integer :: k_fstnzbridx_col ! First non-zero block index in kernel matrix
    integer :: k_lstnzbridx_col ! Last non-zero block index in kernel matrix
    integer :: node
    integer :: nodeshift ! So to spread out accesses to each node
    integer :: local_atom ! Local index of atom
    integer :: local_atomidx ! Global index of atom
    integer :: local_ngwf ! Local index of NGWF
    integer :: ddidx ! Global index of NGWF D,d

    call utils_trace_in('hf_exchange_construct_plan')

    !qoh: Distribute information about batches to all processors
    if (bbngwfidx /=0) &
         atomb_on_node(pub_my_node_id) = ngwf_basis%atom_of_func(bbngwfidx)
    dds_on_node(:,pub_my_node_id) = dds_in_batch
    atomds_on_node = 0
    do node=0,pub_total_num_nodes-1
       call comms_bcast(node,atomb_on_node(node))
       call comms_bcast(node,dds_on_node(:,node))
       where (dds_on_node(:,node) /= 0) &
            atomds_on_node(:,node) = ngwf_basis%atom_of_func(dds_on_node(:,node))
    end do

    cc_send_plan = .false.
    dd_send_plan = .false.

    !qoh: Construct plan
    do local_atom=1,pub_num_atoms_on_node(pub_my_node_id)
       !qoh: First and last non-zero blocks in sparse matrices
       o_fstnzbridx_col = overlap_idx(local_atom)
       o_lstnzbridx_col = overlap_idx(local_atom+1) - 1
       k_fstnzbridx_col = denskern_idx(local_atom)
       k_lstnzbridx_col = denskern_idx(local_atom+1) - 1
       local_atomidx = local_atom+pub_first_atom_on_node(pub_my_node_id)-1
       nodes: do node=0, pub_total_num_nodes-1
          !qoh: Does B overlap with local atom?
          if (any(overlap_idx(o_fstnzbridx_col:o_lstnzbridx_col) == &
               atomb_on_node(node))) then
             !qoh: Is K^{DC} non zero?
             ccdds: do ddidx=1,batch_size
                if (any(denskern_idx(k_fstnzbridx_col:k_lstnzbridx_col) == &
                     atomds_on_node(ddidx,node))) &
                     cc_send_plan(local_atom,node) =  .true.
             end do ccdds
             !qoh: Make sure NGWF B,b is sent!
             if (local_atomidx == atomb_on_node(node)) &
                  cc_send_plan(local_atom,node) =  .true.
          end if
          !qoh: Do we have to a dd NGWF to send?
          dds: do ddidx=1,batch_size
             if (atomds_on_node(ddidx,node) == local_atomidx) then
                local_ngwf = dds_on_node(ddidx,node) &
                     - ngwf_basis%first_on_node(pub_my_node_id) + 1
                dd_send_plan(local_ngwf,node) = .true.
             end if
          end do dds
       end do nodes
    end do

    !qoh: Send relevant parts of plans to other nodes
    do nodeshift=1, pub_total_num_nodes-1
       node = modulo(pub_my_node_id+nodeshift,pub_total_num_nodes)
       call comms_send(node,cc_send_plan(:,node))
       call comms_send(node,dd_send_plan(:,node))
    end do

    !qoh: Receive relevant parts of plans from other nodes
    do nodeshift=1, pub_total_num_nodes-1
       node = modulo(pub_my_node_id+nodeshift,pub_total_num_nodes)
       call comms_recv(node,cc_recv_plan(:,node))
       call comms_recv(node,dd_recv_plan(:,node))
    end do

    !qoh: Local transfer of plans
    cc_recv_plan(:,pub_my_node_id) = cc_send_plan(:,pub_my_node_id)
    dd_recv_plan(:,pub_my_node_id) = dd_send_plan(:,pub_my_node_id)

    call comms_free
    call comms_barrier

    call utils_trace_out('hf_exchange_construct_plan')

  end subroutine hf_exchange_construct_plan

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_fftpot_ngwf_ket(ket_in_fftbox, &
       prod_in_fftbox, bb_fftbox, cc_fftbox_sum_batch, &
       dd1_fftbox, dd2_fftbox, & !ws
       ngwf_basis, ngwfs_on_grid, batch_size, current_batch_size, &
       dd_on_grid_buffer, &!commsbufs
       dds_in_batch, dd_send_plan, dd_recv_plan, zwork, atoma, atomb)

    !==========================================================================!
    ! This subroutine calculates and accumulates (for all \phi_{D,d}s in the   !
    ! current batch) the ket for the current \phi_{B,b} in the AB exchange     !
    ! matrix atom block using FFTs for the integral.                           !
    ! We calculate: \sum_{\substack{D \\ S_{AD} \neq 0}}\sum_d\phi_{D,d}(1)    !
    ! \left[\int\frac{\phi^*_{B,b}(2)                                          !
    ! \sum_{\substack{C \\ S_{BC} \neq 0 \\K^{CD} \neq 0}}\sum_c\phi_{C,c}^*(2)!
    ! K^{C,cD,d}}{\vert \vect{r}_1 - \vect{r}_2 \vert} d\vect{r}_2 \right]     !
    ! The sum of C,c NGWFs (cc_fftbox_sum_batch) is interpolated and multiplied!
    ! with the interpolated B,b NGWF (bb_fftbox).  The integral with the       !
    ! (cutoff) coulomb operator is then calculated in reciprocal space using   !
    ! FFTs.  This potential is then multiplied by the corresponding D,d NGWF   !
    ! and accumulated in ket_in_fftbox.                                        !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2008                                          !
    ! Modified to run in parallel by Quintin Hill in January 2009              !
    !==========================================================================!

    use basis, only: basis_copy_function_to_box, basis_ket_start_wrt_fftbox,&
         basis_location_func_wrt_cell, basis_location_fb_wrt_box
    use comms, only: pub_total_num_nodes, pub_my_node_id,&
         comms_recv, comms_free, comms_send
    use constants, only: DP, PI
    use fourier, only: fourier_apply_box_pair
    use function_basis, only: FUNC_BASIS, pub_buffer_sphere
    use simulation_cell, only: pub_fftbox, pub_cell
    use timer, only: timer_clock

    implicit none

    !qoh: Arguments
    integer,       intent(in)  :: batch_size
    real(kind=DP), intent(in)  :: cc_fftbox_sum_batch(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3,pub_cell%num_spins,batch_size)
    real(kind=DP),intent(inout):: ket_in_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3,pub_cell%num_spins)
    real(kind=DP),intent(out)  :: prod_in_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3,2)
    real(kind=DP), intent(in)  :: bb_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(out)  :: dd1_fftbox(pub_fftbox%total_ld1,& ! jd: workspace
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(out)  :: dd2_fftbox(pub_fftbox%total_ld1,& ! jd: workspace
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    integer,       intent(in)  :: dds_in_batch(batch_size) ! index of \phi_{D,d}
    type(FUNC_BASIS),  intent(in)  :: ngwf_basis
    real(kind=DP), intent(in)  :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    integer,       intent(in)  :: current_batch_size
    integer,       intent(in)  :: atoma, atomb
    logical,       intent(in)  :: dd_send_plan(ngwf_basis%max_on_node,&
         0:pub_total_num_nodes-1)
    logical,       intent(in)  :: dd_recv_plan(ngwf_basis%max_on_node,&
         0:pub_total_num_nodes-1)

    !qoh: Workspace
    real(kind=DP), intent(out) :: dd_on_grid_buffer&
         (ngwf_basis%max_n_ppds_sphere*pub_cell%n_pts)
    complex(kind=DP),intent(out) :: zwork(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    !Allocated as comp_fftbox

    ! Local variables
    integer :: aa_start1, aa_start2, aa_start3
    integer :: aa_cell_start1, aa_cell_start2, aa_cell_start3
    integer :: dd_start1, dd_start2, dd_start3
    integer :: dd_cell_start1, dd_cell_start2, dd_cell_start3
    integer :: is ! Spin counter
    real(kind=DP), parameter :: fourpi = 4.0_DP*PI
    integer :: batch_count ! Number of batch items processed so far
    integer :: ddngwfidx ! Global index of NGWF D,d
    integer :: local_dd ! Local index of NGWF D,d
    integer :: local_dd_start ! Starting index of dd in ngwfs_on_grid
    integer :: remote_dd ! Remote index of NGWF D,d
    integer :: n_dd_ppds ! Number of ppds in sphere of D,d
    integer :: dd_npts ! Number of points in D,d
    integer :: send_dd ! Global NGWF index of sent D,d NGWF (MPI tag)
    integer :: node
    integer :: nodeshift
    integer :: dd1_batchidx(1) ! Batch index of D,d in fftbox 1
    integer :: dd2_batchidx(1) ! Batch index of D,d in fftbox 2

    call utils_trace_in('hf_exchange_fftpot_ngwf_ket')

    local_dd = 0

    ! qoh: Send relevant NGWFs
    send_node: do nodeshift=1, pub_total_num_nodes-1
       node = modulo(pub_my_node_id+nodeshift,pub_total_num_nodes)
       loc_dd: do local_dd=1,ngwf_basis%num_on_node(pub_my_node_id)
          insplan: if (dd_send_plan(local_dd,node)) then
             n_dd_ppds = ngwf_basis%spheres(local_dd)%n_ppds_sphere
             !qoh: Send to node
             send_dd = local_dd + ngwf_basis%first_on_node(pub_my_node_id) - 1
             local_dd_start = ngwf_basis%spheres(local_dd)%offset
             ! ndmh: send list of ppds straight from sphere
             call comms_send(node, ngwf_basis%spheres(local_dd)%ppd_list, &
                  2*ngwf_basis%n_ppds_sphere(send_dd), tag=send_dd)
             ! qoh: send list of values in points in ppds
             dd_npts = ngwf_basis%n_ppds_sphere(send_dd) * pub_cell%n_pts
             call comms_send(node, ngwfs_on_grid(local_dd_start), &
                  dd_npts, tag=send_dd)

          end if insplan
       end do loc_dd
    end do send_node

    !qoh: Find position of \phi_{A,a} function wrt simulation cell
    if (atoma /=0) &
         call basis_location_func_wrt_cell(aa_cell_start1, aa_cell_start2, &
         aa_cell_start3, ngwf_basis%all_tbs(ngwf_basis%first_on_atom(atoma)))

    !qoh: Find position of \phi_{A,a} function wrt fftbox
    call basis_ket_start_wrt_fftbox(aa_start1,aa_start2,aa_start3, &
         pub_fftbox%total_pt1, pub_fftbox%total_pt2, pub_fftbox%total_pt3)

    batch_count = 0
    ! qoh: Receive relevant NGWFs
    recv_node: do nodeshift=1, pub_total_num_nodes
       node = modulo(pub_my_node_id+nodeshift,pub_total_num_nodes)
       rem_dd: do remote_dd=1,ngwf_basis%num_on_node(node)
          inrplan: if (dd_recv_plan(remote_dd,node)) then
             ddngwfidx = ngwf_basis%first_on_node(node) + remote_dd -1
             dd_on_grid_buffer=0.0_DP

             if ( node == pub_my_node_id )  then
                local_dd = ddngwfidx - ngwf_basis%first_on_node(node) + 1
             else
                ! ndmh: receive straight into pub_buffer_sphere
                call comms_recv(node, pub_buffer_sphere%ppd_list, &
                     2*ngwf_basis%n_ppds_sphere(ddngwfidx), tag=ddngwfidx)
                n_dd_ppds = ngwf_basis%n_ppds_sphere(ddngwfidx)
                pub_buffer_sphere%n_ppds_sphere = n_dd_ppds
                pub_buffer_sphere%offset = 1
                !qoh: Receive NGWFs from node
                dd_npts = n_dd_ppds * pub_cell%n_pts
                call comms_recv(node, dd_on_grid_buffer, &
                     dd_npts, tag=ddngwfidx)
             end if

             !qoh: Find position of \phi_{D,d} wrt simulation cell
             call basis_location_func_wrt_cell( &
                  dd_cell_start1, dd_cell_start2, dd_cell_start3, &
                  ngwf_basis%all_tbs(ddngwfidx))

             !qoh: Find where to deposit \phi_{D,d} in fftbox
             call basis_location_fb_wrt_box( &
                  dd_start1, dd_start2, dd_start3, &
                  aa_start1, aa_start2, aa_start3, &
                  aa_cell_start1, aa_cell_start2, &
                  aa_cell_start3, dd_cell_start1, &
                  dd_cell_start2, dd_cell_start3, &
                  pub_cell%total_pt1, pub_cell%total_pt2, pub_cell%total_pt3)

             batch_count = batch_count + 1

             if (modulo(batch_count,2) == 1) then

                dd1_batchidx = minloc(dds_in_batch, &
                     mask = dds_in_batch==ddngwfidx)
                dd1_fftbox =0.0_DP
                dd1_fftbox = 0.0_DP

                !qoh: Use new basis_copy function to transfer \phi_{D,d}
                !qoh: to an FFT box
                if (node == pub_my_node_id) then
                   call basis_copy_function_to_box(dd1_fftbox, &
                        pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                        pub_fftbox%total_pt3, &
                        dd_start1, dd_start2, dd_start3, &
                        ngwf_basis%all_tbs(ddngwfidx),&
                        ngwfs_on_grid, ngwf_basis%spheres(local_dd))
                else
                   call basis_copy_function_to_box(dd1_fftbox, &
                        pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                        pub_fftbox%total_pt3, &
                        dd_start1, dd_start2, dd_start3, &
                        ngwf_basis%all_tbs(ddngwfidx),&
                        dd_on_grid_buffer, pub_buffer_sphere)
                end if

                if (batch_count /= current_batch_size) cycle rem_dd
             else
                dd2_batchidx = minloc(dds_in_batch, &
                     mask = dds_in_batch==ddngwfidx)
                dd2_fftbox = 0.0_DP
                dd2_fftbox = 0.0_DP
                !qoh: use new basis_copy function to transfer \phi_{D,d} to
                !qoh: an FFT box
                if (node == pub_my_node_id) then
                   call basis_copy_function_to_box(dd2_fftbox, &
                        pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                        pub_fftbox%total_pt3, &
                        dd_start1, dd_start2, dd_start3, &
                        ngwf_basis%all_tbs(ddngwfidx),&
                        ngwfs_on_grid, ngwf_basis%spheres(local_dd))
                else
                   call basis_copy_function_to_box(dd2_fftbox, &
                        pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                        pub_fftbox%total_pt3, &
                        dd_start1, dd_start2, dd_start3, &
                        ngwf_basis%all_tbs(ddngwfidx),&
                        dd_on_grid_buffer, pub_buffer_sphere)
                end if
             end if

             !qoh:  Interpolate \phi_{D,d}

             ! jd: I think the deal with the modulo is to stuff things in
             !     pairs to use the fact that apply_box_pair can process
             !     two sets of data at once. 

             spins: do is=1,pub_cell%num_spins
                prod_in_fftbox = 0.0_DP
                if (modulo(batch_count,2) == 0) then

                   !qoh: interpolate sums of \phi_{C,c}
                   !jd: Rather looks like phi_Bb * cc_sum
                   prod_in_fftbox(:,:,:,2) = &
                        cc_fftbox_sum_batch(:,:,:,is,dd2_batchidx(1)) &
                        * bb_fftbox
                else
                   !qoh: interpolate \phi_{D,d} and sum of \phi_{C,c}
                   prod_in_fftbox(:,:,:,2) = 0.0_DP
                end if

                prod_in_fftbox(:,:,:,1) = cc_fftbox_sum_batch&
                     (:,:,:,is,dd1_batchidx(1)) &
                     * bb_fftbox

                if ((atoma /= atomb) .and. .not. (pub_fftbox%coin1 .and. &
                     pub_fftbox%coin2 .and. pub_fftbox%coin3)) &
                     call hf_exchange_shift_product(prod_in_fftbox, &
                     ngwf_basis, atoma, atomb)

                zwork = (0.0_DP,0.0_DP)

                ! jd: 'C' is coarse, not 'complex'

                !qoh: Fourier transform pseudodensity to reciprocal space
                call fourier_apply_box_pair('C','Forward',&
                     prod_in_fftbox(:,:,:,1),&
                     prod_in_fftbox(:,:,:,2),zwork)

                !qoh: Multiply by Coulomb potential in reciprocal space
                zwork = zwork* pub_fftbox%recip_grid(6,:,:,:) * fourpi
                prod_in_fftbox = 0.0_DP

                !qoh: Fourier transform pseudodensity to reciprocal space
                !jd:  Probably meant 'back to real space'
                call fourier_apply_box_pair('C','Backward',&
                     prod_in_fftbox(:,:,:,1), prod_in_fftbox(:,:,:,2),&
                     zwork)

                ket_in_fftbox(:,:,:,is) = ket_in_fftbox(:,:,:,is) + &
                     (dd1_fftbox * prod_in_fftbox(:,:,:,1))
                if (modulo(batch_count,2) == 0 ) &
                     ket_in_fftbox(:,:,:,is) = ket_in_fftbox(:,:,:,is)&
                     + (dd2_fftbox * prod_in_fftbox(:,:,:,2))
             end do spins

          end if inrplan
       end do rem_dd
    end do recv_node

    !qoh: Wait for the other nodes to catch up!
    call timer_clock('hf_exchange_comms_wait',1)
    call comms_free
    call timer_clock('hf_exchange_comms_wait',2)

    call utils_trace_out('hf_exchange_fftpot_ngwf_ket')

  end subroutine hf_exchange_fftpot_ngwf_ket

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_atomblock_elements(bracket_els, &       ! output
       ket_in_fftbox, ngwfs_on_grid, atoma, ngwf_basis, &   ! input
       ket_on_bragrid_buffer)   ! Workspace

    !==========================================================================!
    ! This subroutine calculates the elements in the a column of the BA        !
    ! atomblock. The filrered ket_in_fftbox is outputted for potential use in  !
    ! the calculation of the NGWF gradient.  Inspired by                       !
    !  integrals_brappd_ketfftbox.                                             !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2008.                                         !
    !==========================================================================!

    use basis, only: basis_ket_start_wrt_fftbox, basis_extract_function_from_box
    use comms, only: pub_my_node_id
    use constants, only: DP
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: pub_fftbox, pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    type(FUNC_BASIS), intent(in)  :: ngwf_basis
    real(kind=DP), intent(out)    :: bracket_els(ngwf_basis%max_on_atom, &
         pub_cell%num_spins)
    real(kind=DP), intent(in)     :: ket_in_fftbox(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2, pub_fftbox%total_pt3,pub_cell%num_spins)
    real(kind=DP), intent(in)     :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    integer,       intent(in)     :: atoma

    ! Workspace
    real(kind=DP), intent(out) :: ket_on_bragrid_buffer&
         (ngwf_basis%max_n_ppds_sphere*pub_cell%n_pts)

    real(kind=DP) :: bracket_el ! matrix element
    integer :: is ! spin counter
    integer :: bra_start1, bra_start2, bra_start3
    integer :: dot_npts ! number of points to dot
    integer :: ipt  ! point loop counter
    integer :: ngwf_atoma ! NGWF on atom A
    integer :: aangwfidx ! Global index of NGWF A,a
    integer :: local_aa ! Local index of NGWF A,a
    integer :: bra_start ! Index of start of bra points in list of NGWF points

    call utils_trace_in('hf_exchange_atomblock_elements')

    !qoh: Find the position of the bra in the FFT box
    call basis_ket_start_wrt_fftbox(bra_start1,bra_start2,bra_start3,&
         pub_fftbox%total_pt1,pub_fftbox%total_pt2,pub_fftbox%total_pt3)

    spins: do is=1,pub_cell%num_spins

       aangwfidx = ngwf_basis%first_on_atom(atoma)
       local_aa = aangwfidx - ngwf_basis%first_on_node(pub_my_node_id) + 1
       ket_on_bragrid_buffer=0.0_DP

       !cks: extract ppds belonging to bra function from ket fftbox
       call basis_extract_function_from_box(ket_on_bragrid_buffer, &
            pub_fftbox%total_ld1, pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
            ket_in_fftbox(:,:,:,is), ngwf_basis%spheres(local_aa), &
            ngwf_basis%all_tbs(aangwfidx), bra_start1, bra_start2, bra_start3, 1)

       !qoh: Loop over ngwfs of atom A
       loop_af: do ngwf_atoma = 1,ngwf_basis%num_on_atom(atoma)

          !cks: ddot ppds - these bra and ket representations have the same ppds!
          dot_npts = ngwf_basis%spheres(local_aa)%n_ppds_sphere * pub_cell%n_pts

          bracket_el = 0.0_DP
          bra_start = ngwf_basis%spheres(local_aa)%offset
          do ipt=0,dot_npts-1
             bracket_el = bracket_el + &
                  ngwfs_on_grid(bra_start+ipt) * ket_on_bragrid_buffer(ipt+1)
          end do

          !cks: scale with grid point weight
          bracket_els(ngwf_atoma,is) = pub_cell%weight * bracket_el
          local_aa = local_aa + 1
          aangwfidx = aangwfidx + 1

       end do loop_af

    end do spins

    call utils_trace_out('hf_exchange_atomblock_elements')

  end subroutine hf_exchange_atomblock_elements

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_gradient(contra_grad, cov_grad, & !output
       ket_in_fftbox, atoma, atomb, ngwf_atomb, ngwf_basis, & !input
       denskern, tc_denskern, precond_func_recip, & ! input
       ket_on_bragrid_buffer, k_bablk, tc_k_bablk, work_fftbox) !workspace

    !==========================================================================!
    ! This subroutine calculates the NGWF gradinet for the a NGWFs in the BA   !
    ! atomblock. The filtered ket_in_fftbox from the calculation of the matrix !
    ! elements.  Reciprocal space preconditioning is applied if necessary.     !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in late 2008 and January 2009.                   !
    !==========================================================================!

    use basis, only: basis_ket_start_wrt_fftbox, sphere, &
         basis_extract_function_from_box, basis_copy_sphere, &
         basis_sphere_deallocate, basis_clean_function
    use comms, only: pub_my_node_id
    use constants, only: DP
    use fourier, only: fourier_apply_box_pair
    use function_basis, only: FUNC_BASIS
    use geometry, only: point
    use rundat, only: precond_recip
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: sparse_get_block, SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check
    use xc, only: pub_hfxfraction

    implicit none

    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP),intent(inout) :: contra_grad(ngwf_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP),intent(inout) :: cov_grad(ngwf_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP),intent(inout) :: ket_in_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3, pub_cell%num_spins)
    integer,      intent(in)    :: atoma
    integer,      intent(in)    :: atomb
    integer,      intent(in)    :: ngwf_atomb
    type(SPAM3),  intent(in)    :: denskern(pub_cell%num_spins)
    type(SPAM3),  intent(in)    :: tc_denskern(pub_cell%num_spins)
    real(kind=DP),intent(in)    :: precond_func_recip(:,:,:)

    ! Workspace
    real(kind=DP), intent(out)  :: ket_on_bragrid_buffer&
         (ngwf_basis%max_n_ppds_sphere*pub_cell%n_pts)
    real(kind=DP), intent(out)  :: k_bablk(ngwf_basis%max_on_atom,&
         ngwf_basis%max_on_atom,pub_cell%num_spins) !K^{A,B}
    real(kind=DP), intent(out)  :: tc_k_bablk(ngwf_basis%max_on_atom,&
         ngwf_basis%max_on_atom,pub_cell%num_spins) !KS_A^B
    real(kind=DP), intent(out)  :: work_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)

    complex(kind=DP), allocatable  :: zwork_fftbox(:,:,:) ! complex workspace

    type(SPHERE)  :: asphere ! Sphere for atom A
    real(kind=DP) :: dk_aabbel  ! K_{A,aB,b}
    real(kind=DP) :: tc_dk_aabbel ! SK_{A,aB,b}
    integer :: is ! spin counter
    integer :: bra_start1, bra_start2, bra_start3
    integer :: dot_npts ! Number of points to be dotted
    integer :: ipt  ! point loop counter
    integer :: ngwf_atoma ! NGWF on atom A
    integer :: aangwfidx ! Global index of NGWF A,a
    integer :: local_aa ! Local index of NGWF A,a
    integer :: bra_start ! Start of bra NGWF
    integer :: ierr ! error flag

    call utils_trace_in('hf_exchange_gradient')

    !qoh: initialise k block variable
    k_bablk = 0.0_DP
    tc_k_bablk = 0.0_DP
    aangwfidx = ngwf_basis%first_on_atom(atoma)
    local_aa = aangwfidx - ngwf_basis%first_on_node(pub_my_node_id) + 1

    !qoh: Find the position of the bra in the FFT box.
    call basis_ket_start_wrt_fftbox(bra_start1,bra_start2,bra_start3,&
         pub_fftbox%total_pt1,pub_fftbox%total_pt2,pub_fftbox%total_pt3)

    call basis_copy_sphere(asphere,ngwf_basis%spheres(local_aa),1)

    do is=1,pub_cell%num_spins
       call sparse_get_block(k_bablk(:,:,is),denskern(is),atomb,atoma)
       call sparse_get_block(tc_k_bablk(:,:,is),tc_denskern(is),atomb,atoma)
    end do

    !qoh: Apply scaling for exchange to all density kernel elements
    !qoh: Scaling is 4 (for the 4 NGWFs) * 0.5 (Hartree:Exchange ratio)
    !qoh: * 0.5 (double counting of electron pairs) * (2/num_spins)^2
    !qoh: (scale both density kernels by spin factor) * hfxfraction *
    !qoh: grid point weight * num_spins
    k_bablk = k_bablk * pub_hfxfraction * pub_cell%weight * 4.0_DP /&
         real((pub_cell%num_spins),kind=DP)
    tc_k_bablk = tc_k_bablk * pub_hfxfraction * pub_cell%weight * 4.0_DP /&
         real((pub_cell%num_spins),kind=DP)

    spins: do is=1,pub_cell%num_spins

       aangwfidx = ngwf_basis%first_on_atom(atoma)
       local_aa = aangwfidx - ngwf_basis%first_on_node(pub_my_node_id) + 1
       ket_on_bragrid_buffer=0.0_DP

       !cks: extract ppds belonging to bra function from ket fftbox
       call basis_extract_function_from_box(ket_on_bragrid_buffer, &
            pub_fftbox%total_ld1, pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
            ket_in_fftbox(:,:,:,is), ngwf_basis%spheres(local_aa), &
            ngwf_basis%all_tbs(aangwfidx), bra_start1, bra_start2, bra_start3, 1)

       !cks: shaving - stage 2 / zero points outside NGWF sphere in PPD rep.
       if (pub_cell%n_pts > 1) call basis_clean_function(ket_on_bragrid_buffer, &
            asphere, ngwf_basis%max_n_ppds_sphere)

       !qoh: Loop over NGWFs of atom A
       loop_af: do ngwf_atoma = 1,ngwf_basis%num_on_atom(atoma)
          !qoh: Denskern is symmetric K^{A,aB,b} = K^{B,bA,a}
          dk_aabbel = k_bablk(ngwf_atomb,ngwf_atoma,is)
          !qoh: TC Denskern is anti-symmetric K^{A,aB,b} = -K^{B,bA,a}
          tc_dk_aabbel = -tc_k_bablk(ngwf_atomb,ngwf_atoma,is)

          dot_npts = ngwf_basis%spheres(local_aa)%n_ppds_sphere * pub_cell%n_pts

          bra_start = ngwf_basis%spheres(local_aa)%offset
          do ipt=0,dot_npts-1
             contra_grad(bra_start+ipt) = contra_grad(bra_start+ipt)- &
                  dk_aabbel * ket_on_bragrid_buffer(ipt+1)
             if (.not. precond_recip) &
                  cov_grad(bra_start+ipt) = cov_grad(bra_start+ipt)&
                  - tc_dk_aabbel * ket_on_bragrid_buffer(ipt+1)
          end do

          aangwfidx = aangwfidx + 1
          local_aa = local_aa + 1

       end do loop_af

    end do spins

    if (precond_recip) then

       work_fftbox = 0.0_DP
       allocate(zwork_fftbox(pub_fftbox%total_ld1,&
            pub_fftbox%total_ld2, pub_fftbox%total_pt3), stat=ierr)
       call utils_alloc_check('hf_exchange_gradient','zwork_fftbox',ierr)
       zwork_fftbox = (0.0_DP,0.0_DP)

       if (pub_cell%num_spins == 2) then
          !qoh: Forward FFT the covariant gradient to reciprocal space
          call fourier_apply_box_pair('C','Forward',ket_in_fftbox(:,:,:,1),&
               ket_in_fftbox(:,:,:,2),zwork_fftbox)
          !qoh: Apply kinetic energy preconditioning to covariant gradient
          zwork_fftbox = zwork_fftbox * precond_func_recip
          !qoh:  Backward FFT the covariant gradient to real space
          call fourier_apply_box_pair('C','Backward',ket_in_fftbox(:,:,:,1),&
               ket_in_fftbox(:,:,:,2),zwork_fftbox)
       else
          !qoh: Forward FFT the covariant gradient to reciprocal space
          call fourier_apply_box_pair('C','Forward',ket_in_fftbox(:,:,:,1),&
               work_fftbox,zwork_fftbox)
          !qoh: Apply kinetic energy preconditioning to covariant gradient
          zwork_fftbox = zwork_fftbox * precond_func_recip
          !qoh: Backward FFT the covariant gradient to real space
          call fourier_apply_box_pair('C','Backward',ket_in_fftbox(:,:,:,1),&
               work_fftbox,zwork_fftbox)
          deallocate(zwork_fftbox, stat=ierr)
          call utils_dealloc_check('hf_exchange_gradient','zwork_fftbox',ierr)
       end if

       spins2: do is=1,pub_cell%num_spins

          aangwfidx = ngwf_basis%first_on_atom(atoma)
          local_aa = aangwfidx - ngwf_basis%first_on_node(pub_my_node_id) + 1

          ket_on_bragrid_buffer=0.0_DP
          !cks: extract ppds belonging to bra function from ket fftbox
          call basis_extract_function_from_box(ket_on_bragrid_buffer, &
               pub_fftbox%total_ld1,pub_fftbox%total_ld2,pub_fftbox%total_pt3, &
               ket_in_fftbox(:,:,:,is), ngwf_basis%spheres(local_aa), &
               ngwf_basis%all_tbs(aangwfidx), bra_start1, bra_start2, bra_start3, 1)

          !cks: shaving - stage 2 / zero points outside NGWF sphere in PPD rep.
          if (pub_cell%n_pts > 1) call basis_clean_function(ket_on_bragrid_buffer,&
               asphere, ngwf_basis%max_n_ppds_sphere)

          !qoh: Loop over NGWFs of atom A
          loop_af2: do ngwf_atoma = 1,ngwf_basis%num_on_atom(atoma)
             !qoh: Denskern is symmetric K^{A,aB,b} = K^{B,bA,a}
             dk_aabbel = k_bablk(ngwf_atomb,ngwf_atoma,is)
             tc_dk_aabbel = tc_k_bablk(ngwf_atomb,ngwf_atoma,is)

             dot_npts = ngwf_basis%spheres(local_aa)%n_ppds_sphere * pub_cell%n_pts

             bra_start = ngwf_basis%spheres(local_aa)%offset
             do ipt=0,dot_npts-1
                cov_grad(bra_start+ipt) = cov_grad(bra_start+ipt)- &
                     tc_dk_aabbel * ket_on_bragrid_buffer(ipt+1)
             end do

             aangwfidx = aangwfidx + 1
             local_aa = local_aa + 1

          end do loop_af2

       end do spins2

    end if

    !qoh: Destroy the sphere we created above
    call basis_sphere_deallocate(asphere)

    call utils_trace_out('hf_exchange_gradient')

  end subroutine hf_exchange_gradient

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_num_sph_functions(radtable, num_sw)

    !=========================================================================!
    ! This subroutine sets max_num_bessels and max_sw_set_size.               !
    !-------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2009.                                        !
    !=========================================================================!

    use comms, only: comms_reduce, pub_my_node_id, pub_on_root
    use constants, only: stdout
    use parallel_strategy, only: pub_elements_on_node, pub_num_atoms_on_node
    use rundat, only: cutoff_energy, pub_hfx_max_l, pub_hfx_max_zeros
    use simulation_cell, only: pub_cell
    use spherical_wave, only: sw_bessel_zeros_init, sw_bessel_zeros

    implicit none

    real(kind=DP), intent(out) :: radtable(pub_cell%num_species)
    ! Table giving radius of each species
    integer, intent(out), optional :: num_sw(pub_cell%num_species)

    integer :: lval ! Quantum number l
    integer :: zidx ! Zero index
    integer :: ridx ! Radius index
    integer :: eidx ! Element index
    integer :: nbessels  ! Number of spherical bessels (per species)
    integer :: nsphwaves ! Number of spherical waves (per species)
    real(kind=DP) :: maxzero ! Maximum acceptable Bessel zero

    call utils_trace_in('hf_exchange_num_sph_functions')

    max_num_bessels = 0
    max_sw_set_size = 0

    call sw_bessel_zeros_init(pub_hfx_max_zeros,pub_hfx_max_l)

    radtable = -1.0_DP
    do ridx = 1, pub_cell%num_species
       do eidx=1,pub_num_atoms_on_node(pub_my_node_id)
          if (pub_elements_on_node(eidx)%species_number == ridx) then
             radtable(ridx)=pub_elements_on_node(eidx)%radius
             exit
          end if
       end do
    end do
    !qoh: Get radii of all species
    call comms_reduce('MAX',radtable,pub_cell%num_species)

    do ridx = 1,pub_cell%num_species
       nbessels = 0
       nsphwaves = 0

       ! qoh: qmax = sqrt(2.0_DP*cutoff_energy)
       maxzero = sqrt(2.0_DP*cutoff_energy)*radtable(ridx)

       do lval = lmin, pub_hfx_max_l
          do zidx=1, pub_hfx_max_zeros
             if (sw_bessel_zeros(zidx,lval) .lt. maxzero) then
                nbessels = nbessels + 1
                nsphwaves=nsphwaves + 1 + 2*lval
             else
                if (pub_on_root) then
                   write(stdout,'(a,i0,a,i0,a)') &
                        'WARNING: Spherical waves with l=',lval,' and qi=', &
                        zidx,' will be ignored.'
                end if
             end if
          end do
       end do
       max_num_bessels = max(nbessels,max_num_bessels)
       max_sw_set_size = max(nsphwaves,max_sw_set_size)
       if (present(num_sw)) num_sw(ridx) = nsphwaves
    end do

    if (.not. useoverlapmetric .and. .not. singlecentre) &
         max_sw_set_size = 2*max_sw_set_size

    call utils_trace_out('hf_exchange_num_sph_functions')

  end subroutine hf_exchange_num_sph_functions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_sph_bessels_init(sphbessels,radtable)

    !==========================================================================!
    ! This subroutine initialises sphbessels and the spherical wave module.    !
    ! It defines the spherical Bessel functions and calculates parts of the    !
    ! potential integrals.                                                     !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2009.                                         !
    !==========================================================================!

    use rundat, only: cutoff_energy, pub_hfx_max_l, pub_hfx_max_zeros
    use simulation_cell, only: pub_cell, pub_tb_recip_grid, &
         pub_maxtight_pts1, pub_maxtight_pts2, pub_maxtight_pts3
    use spherical_wave, only: sw_bessel_zeros, sw_init

    implicit none

    type(bessel),intent(out) :: sphbessels(pub_cell%num_species,max_num_bessels)
    real(kind=DP),intent(in) :: radtable(pub_cell%num_species)

    integer :: lval  ! Quantum number l
    integer :: zidx  ! zero index
    integer :: ridx  ! Radius index
    integer :: sbidx ! spherical bessel idx
    real(kind=DP) :: maxzero ! maximum acceptable Bessel zero
    real(kind=DP) :: rmaxsq  ! Maximum arguement of Bessel function squared
    real(kind=DP) :: gmaxsq  ! Maximum square magnitude of G vector
    integer :: n1half, n2half, n3half

    call utils_trace_in('hf_exchange_sph_bessels_init')

    sphbessels%lval = -1

    do ridx=1,pub_cell%num_species

       ! qoh: qmax = sqrt(2.0_DP*cutoff_energy)
       maxzero = sqrt(2.0_DP*cutoff_energy)*radtable(ridx)
       sbidx = 0
       loopl: do lval = lmin, pub_hfx_max_l
          loopz: do zidx = 1, pub_hfx_max_zeros
             if ( sw_bessel_zeros(zidx,lval) .lt. maxzero) then
                sbidx = sbidx + 1
                sphbessels(ridx,sbidx)%lval = lval
                sphbessels(ridx,sbidx)%qval = sw_bessel_zeros(zidx,lval)&
                     /radtable(ridx)
                sphbessels(ridx,sbidx)%aval = radtable(ridx)
                call internal_const_pot
             end if
          end do loopz
       end do loopl
    end do

    ! qoh: Need to set rmax as default is too small

    n1half = pub_maxtight_pts1/2+1
    n2half = pub_maxtight_pts2/2+1
    n3half = pub_maxtight_pts3/2+1

    gmaxsq = 2.0_DP*pub_tb_recip_grid(5,n1half,&
         n2half,n3half)
    rmaxsq = gmaxsq*(maxval(radtable)**2)

    ! qoh: Initialise SW module
    call sw_init(pub_hfx_max_l, pub_hfx_max_zeros, rmaxsq)

    call utils_trace_out('hf_exchange_sph_bessels_init')

  contains

    subroutine internal_const_pot

      ! qoh: Calculate the constant part of the spherical Bessel potential
      ! qoh: integrals.

      ! jd: Formulas verified with Mathematica

      implicit none
      real(kind=DP) :: qa
      real(kind=DP) :: qval

      qa = sphbessels(ridx,sbidx)%qval * sphbessels(ridx,sbidx)%aval
      qval = sphbessels(ridx,sbidx)%qval

      select case (lval)
      case(0)
         sphbessels(ridx,sbidx)%farpotint = (sin(qa) -qa*cos(qa))/ (qval**3) ! jd: OK
         sphbessels(ridx,sbidx)%nearpotint = - cos(qa) / (qval*qval) ! jd: OK
      case(1)
         sphbessels(ridx,sbidx)%farpotint = ((3.0_DP- qa*qa )*sin(qa) &
              - 3.0_DP*qa*cos(qa))/(qval**4) ! jd: OK
         sphbessels(ridx,sbidx)%nearpotint = - sin(qa) /(qval*qa) ! jd: OK
      case(2)
         sphbessels(ridx,sbidx)%farpotint = (qa*(qa*qa-15.0_DP)*cos(qa) - &
              (6.0_DP*qa*qa-15.0_DP)*sin(qa))/(qval**5) ! jd: OK
         sphbessels(ridx,sbidx)%nearpotint = (qa*cos(qa) - sin(qa))/(qa**3) ! jd: OK
      case(3)
         sphbessels(ridx,sbidx)%farpotint = ((qa**4 -45.0_DP*qa*qa + 105.0_DP) &
              * sin(qa) + qa*(10.0_DP*qa*qa - 105.0_DP)*cos(qa))/(qval**6) ! jd: OK
         sphbessels(ridx,sbidx)%nearpotint = qval*(3.0_DP * qa * cos (qa) +&
              (qa*qa-3.0_DP)*sin(qa))/(qa**5) ! jd: OK
      case(4)
         sphbessels(ridx,sbidx)%farpotint = ((15.0_DP * qa**4 - 420.0_DP*qa*qa&
              + 945.0_DP) *sin(qa) &
              -qa*(qa**4 - 105.0_DP*qa*qa + 945.0_DP)*cos(qa))/ (qval**7) ! jd: OK
         sphbessels(ridx,sbidx)%nearpotint = qval*qval*(qa*(15.0_DP - qa*qa)&
              *cos(qa) + (6.0_DP * qa*qa - 15.0_DP) * sin(qa)) / (qa**7) ! jd: OK

      end select

    end subroutine internal_const_pot

  end subroutine hf_exchange_sph_bessels_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function hf_exchange_sw_set_size(sphbessels,centrea,centreb)

    !==========================================================================!
    ! This subroutine calculates the size of the set of spherical waves        !
    ! centred on centrea and (optionally b).                                   !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in September 2009                                !
    !==========================================================================!

    use geometry, only: geometry_distance
    use simulation_cell, only: pub_cell

    implicit none

    integer :: hf_exchange_sw_set_size
    type(bessel), intent(in) :: sphbessels(pub_cell%num_species,&
         max_num_bessels)
    type(atom_centre), intent(in) :: centrea
    type(atom_centre), intent(in), optional ::centreb

    integer :: sbidx
    integer :: lval

    call utils_trace_in('hf_exchange_sw_set_size')

    hf_exchange_sw_set_size = 0

    do sbidx=1, max_num_bessels
       lval = sphbessels(centrea%species_number,sbidx)%lval
       if (lval == -1) exit
       hf_exchange_sw_set_size = hf_exchange_sw_set_size + 2*lval + 1
    end do

    ! qoh: Also include centre B if B is given and A neq B
    if (present(centreb)) then
       if (geometry_distance(centrea%incell,centreb%incell) &
            > epsilon(1.0_DP)) then
          if (singlecentre) hf_exchange_sw_set_size = 0
          do sbidx=1, max_num_bessels
             lval = sphbessels(centreb%species_number,sbidx)%lval
             if (lval == -1) exit
             hf_exchange_sw_set_size = hf_exchange_sw_set_size + 2*lval + 1
          end do
       end if
    end if

    call utils_trace_out('hf_exchange_sw_set_size')

  end function hf_exchange_sw_set_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_init_centre(atomcentre, ngwf_basis, elements, atom, &
       shifttopoint)

    !==========================================================================!
    ! This subroutine initialises centre for the current atom.                 !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in September 2009                                !
    ! Modified to run on multiple cores on 13 October 2009 by Quintin Hill.    !
    !==========================================================================!

    use basis, only: basis_function_origin_wrt_tb, &
         basis_func_centre_wrt_fftbox, basis_ket_start_wrt_fftbox, &
         basis_location_func_wrt_cell
    use function_basis, only: FUNC_BASIS
    use geometry, only: point, operator(+), operator(-), operator(*), &
         operator(.DOT.)
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_orig_atom
    use simulation_cell, only: pub_fftbox, pub_cell, pub_maxtight_pts1, &
         pub_maxtight_pts2, pub_maxtight_pts3

    implicit none

    type(atom_centre), intent(out) :: atomcentre
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(ELEMENT),    intent(in) :: elements(pub_cell%nat)
    integer,          intent(in) :: atom
    type(POINT), optional, intent(out) :: shifttopoint

    integer :: first_ngwf_idx ! Index of first NGWF on atom B
    integer :: atom_species_number
    type(POINT) :: atomcentreinfftbox !Coordinates of NGWF Centre in FFT box
    type(POINT) :: tborigin ! Origin of tightbox
    real(kind=DP) :: tb_orig1, tb_orig2, tb_orig3 ! Origin of SW in TB in pts
    integer :: cell_start1, cell_start2, cell_start3 ! Point in cell of TB start
    integer :: start1, start2, start3 ! Point in FFT box of TB start
    integer :: orig_atom

    call utils_trace_in('hf_exchange_init_centre')

    orig_atom = pub_orig_atom(atom)
    atom_species_number = elements(orig_atom)%species_number
    atomcentre%species_number = atom_species_number
    atomcentre%incell = elements(orig_atom)%centre
    first_ngwf_idx = ngwf_basis%first_on_atom(atom)

    call basis_ket_start_wrt_fftbox(start1,start2,start3,&
         pub_fftbox%total_pt1,pub_fftbox%total_pt2,pub_fftbox%total_pt3)

    call basis_location_func_wrt_cell(cell_start1, &
         cell_start2, cell_start3, ngwf_basis%all_tbs(first_ngwf_idx))

    atomcentre%tbcellstart(1) = cell_start1
    atomcentre%tbcellstart(2) = cell_start2
    atomcentre%tbcellstart(3) = cell_start3

    atomcentreinfftbox = basis_func_centre_wrt_fftbox&
         (atomcentre%incell, start1, start2, start3, &
         cell_start1, cell_start2, cell_start3)

    call basis_function_origin_wrt_tb(tb_orig1, tb_orig2, tb_orig3,  &
         atomcentre%incell, ngwf_basis%all_tbs(first_ngwf_idx))

    atomcentre%intb = (tb_orig1 * pub_fftbox%d1)*pub_fftbox%a1_unit &
         + (tb_orig2 * pub_fftbox%d2)*pub_fftbox%a2_unit &
         + (tb_orig3 * pub_fftbox%d3)*pub_fftbox%a3_unit

    if (present(shifttopoint)) shifttopoint = &
         ((1.0_DP - tb_orig1 + aint(tb_orig1,kind=DP))*pub_fftbox%d1) &
         *pub_fftbox%a1_unit + ((1.0_DP - tb_orig2 + aint(tb_orig2,kind=DP))&
         *pub_fftbox%d2)*pub_fftbox%a2_unit + ((1.0_DP - tb_orig3 &
         + aint(tb_orig3,kind=DP))*pub_fftbox%d3)*pub_fftbox%a3_unit

    tborigin = atomcentreinfftbox - atomcentre%intb
    atomcentre%tbstart(1) = &
         1 + nint( (tborigin .DOT. pub_fftbox%a1_unit) / pub_fftbox%d1)
    atomcentre%tbend(1) = atomcentre%tbstart(1) + pub_maxtight_pts1 - 1
    atomcentre%tbstart(2) = &
         1 + nint( (tborigin .DOT. pub_fftbox%a2_unit) / pub_fftbox%d2)
    atomcentre%tbend(2) = atomcentre%tbstart(2) + pub_maxtight_pts2 - 1
    atomcentre%tbstart(3) = &
         1 + nint( (tborigin .DOT. pub_fftbox%a3_unit) / pub_fftbox%d3)
    atomcentre%tbend(3) = atomcentre%tbstart(3) + pub_maxtight_pts3 - 1

    atomcentre%radius = elements(orig_atom)%radius

    call utils_trace_out('hf_exchange_init_centre')

  end subroutine hf_exchange_init_centre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_ngwf_product(tb_batch, prod_in_fftbox, bb_fftbox,&
       cc_fftbox_sum_batch, centreb, batch_size, &
       current_batch_size)

    !==========================================================================!
    ! This subroutine calculates a batch of NGWF products depositing them into !
    ! a batch of tightboxes.                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2009.                                         !
    !==========================================================================!

    use constants, only: DP
    use simulation_cell, only: pub_fftbox, pub_cell, &
         pub_maxtight_pts1, pub_maxtight_pts2, pub_maxtight_pts3
    use timer, only: timer_clock

    implicit none

    integer, intent(in)          :: batch_size
    real(kind=DP), intent(in)    :: cc_fftbox_sum_batch(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3,pub_cell%num_spins,&
         batch_size)
    real(kind=DP), intent(in)    :: bb_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(out)   :: tb_batch(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3,batch_size, pub_cell%num_spins)
    type(atom_centre),intent(in) :: centreb
    integer, intent(in)          :: current_batch_size

    ! qoh: Workspace
    real(kind=DP),intent(out)    :: prod_in_fftbox &
         (pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3,2)


    integer :: batch_count
    integer :: is  ! spin counter

    call utils_trace_in('hf_exchange_ngwf_product')

    spins: do is=1,pub_cell%num_spins
       do batch_count = 1,current_batch_size,2
          prod_in_fftbox = 0.0_DP
          if (batch_count /= current_batch_size) then

             prod_in_fftbox(:,:,:,1) = cc_fftbox_sum_batch&
                  (:,:,:,is,batch_count) * bb_fftbox
             prod_in_fftbox(:,:,:,2) = &
                  cc_fftbox_sum_batch(:,:,:,is,batch_count+1) &
                  * bb_fftbox

             tb_batch(:,:,:,batch_count:batch_count+1,is) = &
                  prod_in_fftbox(centreb%tbstart(1):centreb%tbend(1),&
                  centreb%tbstart(2):centreb%tbend(2),&
                  centreb%tbstart(3):centreb%tbend(3),:)

          else
             ! qoh: One sum of \phi_{C,c}
             prod_in_fftbox(:,:,:,2) = 0.0_DP

             prod_in_fftbox(:,:,:,1) = cc_fftbox_sum_batch&
                  (:,:,:,is,batch_count) * bb_fftbox

             tb_batch(:,:,:,batch_count,is) = &
                  prod_in_fftbox(centreb%tbstart(1):centreb%tbend(1),&
                  centreb%tbstart(2):centreb%tbend(2),&
                  centreb%tbstart(3):centreb%tbend(3),1)

          end if
       end do

    end do spins

    call utils_trace_out('hf_exchange_ngwf_product')

  end subroutine hf_exchange_ngwf_product

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_expansion(sphbessels, centrea, centreb, &
       swswoverlap, swcoeff, batch_size, current_batch_size, swswovlpdone, &
       swa_tightbox, tb_batch, tb_zwork, current_sw_set_size, vmatrix)

    !==========================================================================!
    ! This subroutine expands the NGWF product in terms of a set of spherical  !
    ! waves using the appropriate metric.                                      !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2010.                                         !
    !==========================================================================!

    use constants, only: DP
    use simulation_cell, only: pub_cell, pub_maxtight_pts1, pub_maxtight_pts2, &
         pub_maxtight_pts3
    use timer, only: timer_clock

    implicit none

    integer, intent(in) :: batch_size
    real(kind=DP), intent(in) :: tb_batch(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3,batch_size, pub_cell%num_spins)

    ! qoh: Spherical Bessel functions
    type(bessel),intent(in) :: sphbessels(pub_cell%num_species,&
         max_num_bessels)
    type(atom_centre),intent(in) :: centrea, centreb
    real(kind=DP), intent(inout) :: swswoverlap(max_sw_set_size)
    real(kind=DP), intent(out) :: swcoeff(max_sw_set_size,batch_size,&
         pub_cell%num_spins)
    logical, intent(inout) :: swswovlpdone
    integer, intent(in) :: current_batch_size
    integer, intent(in) :: current_sw_set_size
    real(kind=DP), intent(in) :: vmatrix(:,:)

    ! qoh: Workspace
    real(kind=DP), intent(out) :: swa_tightbox(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3)
    complex(kind=DP), intent(out) :: tb_zwork(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3)

    call utils_trace_in('hf_exchange_expansion')

    ! qoh: Find spherical wave coefficients
    if (useoverlapmetric) then
       call timer_clock('hf_exchange_sw_ovlp_metric',1)
       call hf_exchange_sw_ovlp_metric(swcoeff,swswoverlap, &
            swa_tightbox, tb_batch, centreb, sphbessels, &
            batch_size, current_batch_size, tb_zwork, &
            swswovlpdone)
       call timer_clock('hf_exchange_sw_ovlp_metric',2)
    else
       call timer_clock('hf_exchange_sw_es_metric',1)
       call hf_exchange_sw_es_metric(swcoeff,vmatrix,&
            centrea,centreb,sphbessels,tb_batch,batch_size, &
            current_batch_size, current_sw_set_size)
       call timer_clock('hf_exchange_sw_es_metric',2)
    end if

    call utils_trace_out('hf_exchange_expansion')

  end subroutine hf_exchange_expansion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_sw_ovlp_metric(swcoeff,swswoverlap, &
       swa_tightbox, tb_batch, centreb, &
       sphbessels, batch_size, current_batch_size, tb_zwork, &
       swswovlpdone)

    !==========================================================================!
    ! This subroutine expands a batch of NGWF products in tightboxes in terms  !
    ! of a set of spherical waves, returning a batch of expansion coefficients.!
    ! It uses the overlap metric to do this.                                   !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2009.                                         !
    !==========================================================================!

    use rundat, only: pub_hfx_max_zeros
    use simulation_cell, only: pub_cell, pub_maxtight_pts1, pub_maxtight_pts2, &
         pub_maxtight_pts3
    use spherical_wave, only: sw_recp_generate_in_tb

    implicit none

    integer, intent(in) :: batch_size
    real(kind=DP),intent(inout) :: &
         swcoeff(max_sw_set_size,batch_size,pub_cell%num_spins)
    real(kind=DP), intent(inout) :: swswoverlap(max_sw_set_size)
    real(kind=DP), intent(out) :: swa_tightbox&
         (pub_maxtight_pts1,pub_maxtight_pts2,pub_maxtight_pts3)
    real(kind=DP), intent(in) :: tb_batch&
         (pub_maxtight_pts1,pub_maxtight_pts2,pub_maxtight_pts3,&
         batch_size, pub_cell%num_spins)
    complex(kind=DP), intent(out) :: tb_zwork&
         (pub_maxtight_pts1,pub_maxtight_pts2,pub_maxtight_pts3)
    type(atom_centre), intent(in) :: centreb
    type(bessel), intent(in)    :: sphbessels(pub_cell%num_species,&
         max_num_bessels)
    integer, intent(in) :: current_batch_size
    logical, intent(inout) :: swswovlpdone

    ! Local variables
    integer :: swsidx  ! spherical wave set index
    integer :: sbidx   ! spherical Bessel index
    integer :: batch_count
    integer :: lval
    integer :: mval
    integer :: is      ! spin counter

    call utils_trace_in('hf_exchange_sw_ovlp_metric')

    swcoeff = 0.0_DP

    swsidx = 0
    sb: do sbidx=1,max_num_bessels
       lval = sphbessels(centreb%species_number,sbidx)%lval
       if (lval == -1) exit
       mloop: do mval=-lval,lval
          swsidx = swsidx + 1
          call sw_recp_generate_in_tb(swa_tightbox, tb_zwork, &
               lval, mval, sphbessels(centreb%species_number,sbidx)%qval, &
               sphbessels(centreb%species_number,sbidx)%aval,centreb%intb, &
               pub_maxtight_pts1, pub_maxtight_pts2, pub_maxtight_pts3)
          if (shavesw) call hf_exchange_shave_tb(swa_tightbox, centreb, &
               sphbessels(centreb%species_number,1)%aval)
          if (.not. swswovlpdone) &
               swswoverlap(swsidx) = sum(swa_tightbox * swa_tightbox)&
               /pub_cell%weight
          do batch_count = 1,current_batch_size
             do is=1,pub_cell%num_spins
                swcoeff(swsidx,batch_count,is) = &
                     sum(swa_tightbox*tb_batch(:,:,:,batch_count,is)) &
                     /swswoverlap(swsidx)
             end do
          end do
       end do mloop
    end do sb

    if (.not. swswovlpdone .and. pub_hfx_max_zeros == 1) then
       ! @@@@jd: WTF is this?
       print *, "#############################"
       print *, swcoeff
       print *, swswoverlap
       print *, "#############################"
    end if
    swswovlpdone = .true.

    call utils_trace_out('hf_exchange_sw_ovlp_metric')

  end subroutine hf_exchange_sw_ovlp_metric

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_sw_es_metric(swcoeff,vmatrix, &
       centrea,centreb,sphbessels,tb_batch,batch_size, &
       current_batch_size, current_sw_set_size)

    !==========================================================================!
    ! This subroutine expands a batch of NGWF products in tightboxes in terms  !
    ! of a set of spherical waves, returning a batch of expansion coefficients.!
    ! It uses the electrostatic metric to do this.                             !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2009.                                         !
    !==========================================================================!

    use simulation_cell, only: pub_cell, pub_maxtight_pts1, pub_maxtight_pts2, &
         pub_maxtight_pts3
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_dump_array3D_to_file !@@
    use wrappers, only: wrappers_dgesv !, wrappers_dgelss

    implicit none

    integer, intent(in) :: batch_size
    real(kind=DP),intent(inout) :: &
         swcoeff(max_sw_set_size,batch_size,pub_cell%num_spins)
    integer, intent(in) :: current_sw_set_size
    real(kind=DP), intent(in) :: vmatrix(current_sw_set_size,&
         current_sw_set_size)
    real(kind=DP), intent(in) :: tb_batch&
         (pub_maxtight_pts1,pub_maxtight_pts2,pub_maxtight_pts3,&
         batch_size, pub_cell%num_spins)
    type(atom_centre), intent(in) :: centrea, centreb
    type(bessel), intent(in)    :: sphbessels(pub_cell%num_species,&
         max_num_bessels)
    integer, intent(in) :: current_batch_size

    ! Local variables
    integer :: batch_count
    integer :: ierr ! error flag
    integer :: is   ! spin counter
    real(kind=DP), allocatable :: es_metric(:,:,:)
    real(kind=DP), allocatable :: tb_batch_fine(:,:,:,:,:)

    integer :: i,j, in_batch !@removeme

    call utils_trace_in('hf_exchange_sw_es_metric')

    allocate(es_metric(current_sw_set_size,current_batch_size, &
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('hf_exchange_sw_es_metric','es_metric',ierr)

    swcoeff = 0.0_DP

!    allocate(tb_batch_fine(finer*pub_maxtight_pts1, finer*pub_maxtight_pts2, &
!         finer*pub_maxtight_pts3,batch_size, pub_cell%num_spins),stat=ierr)
!    call utils_alloc_check('hf_exchange_sw_es_metric','tb_batch_fine',ierr)

!    write(*,*) '@Interpolating x 3... (spin 1 only)'
!    do in_batch = 1, batch_size
!       call fourier_interpolate_tb_3( &
!           tb_batch(:,:,:,in_batch,1), tb_batch(:,:,:,in_batch,1), &
!           tb_batch_fine(:,:,:,in_batch,1), tb_batch_fine(:,:,:,in_batch,1)) 
!       if(in_batch==1) then
!          call utils_dump_array3D_to_file(tb_batch(:,:,:,in_batch,1),'tb_batch_orig') !@@
!          call utils_dump_array3D_to_file(tb_batch_fine(:,:,:,in_batch,1),'tb_batch_fine') !@@
!       end if
!    end do
!    write(*,*) '@Interpolated...'

    ! jd: Probably calculates a batch of "b"'s (eq. 5.3.18),
    !     as (f_i | P_{beta,delta} )
    call hf_exchange_swpot_batch(tb_batch, es_metric, centrea,&
         centreb, sphbessels, batch_size, current_batch_size,  &
         current_sw_set_size, pub_cell%num_spins, .false.)
!@@@@@@@@@@
!    call hf_exchange_swpot_batch_finer(tb_batch_fine, es_metric, centrea,&
!         centreb, sphbessels, batch_size, current_batch_size,  &
!         current_sw_set_size, pub_cell%num_spins, .false.)

!    es_metric = es_metric * 0.015625_DP ! @@ take fineness into account
!    es_metric = es_metric * 0.125_DP ! @@ take fineness into account
!    es_metric = es_metric * 0.03703703703703703704

!    deallocate(tb_batch_fine,stat=ierr)
!    call utils_dealloc_check('hf_exchange_sw_es_metric','tb_batch_fine',ierr)

    if (printcoeffs) then
       write(*,*) current_sw_set_size,' ',current_batch_size
       print *, "-VVMAT:-============================@"
       do j = 1, current_sw_set_size
          do i = 1, current_sw_set_size
             print *, vmatrix(i,j)
          end do
       end do
       print *, "-BVEC:-=============================@", current_batch_size
       do batch_count=1,current_batch_size
          print *, es_metric(:,batch_count,1)
       end do
       print *, "---------------------------------------"
       print *, "-CVEC:-=============================@"
    end if

    ! jd: solving  Vc = b for c
    ! V = vmatrix
    ! c = swcoeff
    ! b = es_metric

    do is=1,pub_cell%num_spins
       call wrappers_dgesv(&
            swcoeff(1:current_sw_set_size,1:current_batch_size,is), &
            vmatrix, es_metric(:,:,is), current_sw_set_size, current_batch_size)
    end do

    if (printcoeffs) then
       do batch_count=1,current_batch_size
          print *, swcoeff(1:current_sw_set_size,batch_count,1)
       end do
       print *, "---------------------------------------"
    end if

    deallocate(es_metric,stat=ierr)
    call utils_dealloc_check('hf_exchange_sw_es_metric','es_metric',ierr)

    call utils_trace_out('hf_exchange_sw_es_metric')

  end subroutine hf_exchange_sw_es_metric

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_tbpot_ngwf_ket(ket_in_fftbox, &
       pot_in_fftbox, tb_batch, tb_batch2, dd1_fftbox, dd2_fftbox, & !ws
       ngwf_basis, ngwfs_on_grid, batch_size, current_batch_size, &
       dd_on_grid_buffer, dds_in_batch, dd_send_plan, dd_recv_plan, atoma, &
       atomb, centrea, centreb, sphbessels, swcoeff)

    !==========================================================================!
    ! This subroutine calculates the exchange potential from an NGWF product   !
    ! and applies it to batch of NGWFs (\phi_{D,d}).  The potential is         !
    ! calculated from the spherical wave expansion or using the NPA as         !
    ! appropriate.                                                             !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in 2009 and 2010.                                !
    !==========================================================================!

    use basis, only: basis_copy_function_to_box, &
         basis_location_func_wrt_cell, basis_location_fb_wrt_box, &
         basis_ket_start_wrt_fftbox
    use comms, only: pub_total_num_nodes, pub_my_node_id,&
         comms_recv, comms_free, comms_send, comms_abort
    use constants, only: DP, stdout
    use function_basis, only: pub_buffer_sphere, FUNC_BASIS
    use rundat, only: pub_hfx_integration_variant
    use simulation_cell, only: pub_fftbox, pub_cell , pub_maxtight_pts1, &
         pub_maxtight_pts2, pub_maxtight_pts3
    use timer, only: timer_clock

    implicit none

    !qoh: Arguments
    integer,       intent(in)  :: batch_size
    real(kind=DP),intent(inout):: ket_in_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3,pub_cell%num_spins)
    real(kind=DP), intent(out) :: dd1_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    real(kind=DP), intent(out) :: dd2_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3)
    integer,       intent(in)  :: dds_in_batch(batch_size) ! index of \phi_{D,d}
    real(kind=DP), intent(out) :: tb_batch(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3,batch_size, pub_cell%num_spins)
    real(kind=DP), intent(in) :: tb_batch2(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3,batch_size, pub_cell%num_spins)
    type(FUNC_BASIS), intent(in)  :: ngwf_basis
    real(kind=DP), intent(in)  :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    integer,       intent(in)  :: current_batch_size
    integer,       intent(in)  :: atoma
    integer,       intent(in)  :: atomb
    logical,       intent(in)  :: dd_send_plan(ngwf_basis%max_on_node,&
         0:pub_total_num_nodes-1)
    logical,       intent(in)  :: dd_recv_plan(ngwf_basis%max_on_node,&
         0:pub_total_num_nodes-1)

    type(atom_centre), intent(in) :: centrea
    type(atom_centre), intent(in) :: centreb
    ! qoh: Spherical Bessel functions
    type(bessel),  optional, intent(in) :: sphbessels(pub_cell%num_species,&
         max_num_bessels)
    real(kind=DP), optional, intent(in) :: swcoeff(max_sw_set_size,&
         batch_size, pub_cell%num_spins)

    !qoh: Workspace
    real(kind=DP),intent(out)  :: pot_in_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3,2)
    real(kind=DP), intent(out) :: dd_on_grid_buffer&
         (ngwf_basis%max_n_ppds_sphere*pub_cell%n_pts)


    ! Local variables
    integer :: first_a_ngwf_idx
    integer :: aa_start1, aa_start2, aa_start3
    integer :: aa_cell_start1, aa_cell_start2, aa_cell_start3
    integer :: dd_start1, dd_start2, dd_start3
    integer :: dd_cell_start1, dd_cell_start2, dd_cell_start3
    integer :: is ! Spin counter
    integer :: batch_count ! Number of batch items processed so far
    integer :: recv_dd ! Global index of NGWF D,d received
    integer :: local_dd ! Local index of NGWF D,d
    integer :: local_dd_start ! Starting index of dd in ngwfs_on_grid
    integer :: remote_dd ! Remote index of NGWF D,d
    integer :: n_dd_ppds ! Number of ppds in sphere of D,d
    integer :: dd_npts ! Number of points in D,d
    integer :: send_dd ! Global NGWF index of sent D,d NGWF (MPI tag)
    integer :: node
    integer :: nodeshift
    integer :: dd1_batchidx(1) ! Batch index of D,d in fftbox 1
    integer :: dd2_batchidx(1) ! Batch index of D,d in fftbox 2

    call utils_trace_in('hf_exchange_tbpot_ngwf_ket')

    local_dd = 0

    ! qoh: Send relevant NGWFs
    send_node: do nodeshift=1, pub_total_num_nodes-1
       node = modulo(pub_my_node_id+nodeshift,pub_total_num_nodes)
       loc_dd: do local_dd=1,ngwf_basis%num_on_node(pub_my_node_id)
          insplan: if (dd_send_plan(local_dd,node)) then
             n_dd_ppds = ngwf_basis%spheres(local_dd)%n_ppds_sphere
             dd_npts = n_dd_ppds * pub_cell%n_pts
             !qoh: Send to node
             send_dd = local_dd + ngwf_basis%first_on_node(pub_my_node_id) - 1
             local_dd_start = ngwf_basis%spheres(local_dd)%offset
             call comms_send(node,ngwf_basis%spheres(local_dd)%ppd_list, &
                  2*ngwf_basis%n_ppds_sphere(send_dd), tag=send_dd)

             ! qoh: send list of values in points in ppds
             call comms_send(node, ngwfs_on_grid(local_dd_start), &
                  dd_npts, tag=send_dd)

          end if insplan
       end do loc_dd
    end do send_node

    null_atoma: if (atoma /=0) then

       if (usenpa) then
          ! jd: Select the right integration routine
          select case (pub_hfx_integration_variant)
          case(1)
             call timer_clock('hf_exchange_npa_batch_var1',1)
             call hf_exchange_npa_batch_var1(tb_batch, tb_batch2, centrea,&
                  centreb, batch_size)
             call timer_clock('hf_exchange_npa_batch_var1',2)
          case(2)
             call timer_clock('hf_exchange_npa_batch_var2',1)
             call hf_exchange_npa_batch_var2(tb_batch, tb_batch2, centrea,&
                  centreb, batch_size)
             call timer_clock('hf_exchange_npa_batch_var2',2)
          case(3)
             call timer_clock('hf_exchange_npa_batch_var3',1)
             call hf_exchange_npa_batch_var3(tb_batch, tb_batch2, centrea,&
                  centreb, batch_size)
             call timer_clock('hf_exchange_npa_batch_var3',2)
          case default
             write (stdout,'(2a)')&
                  'Unrecognized variant of HF exchange integration.',&
                  'Program execution stops.'
             call comms_abort
          end select
       else
          call timer_clock('hf_exchange_swsetpot_batch',1)
          call hf_exchange_swsetpot_batch(tb_batch, swcoeff, centrea,&
               centreb, sphbessels, batch_size, current_batch_size,atoma,atomb)
          call timer_clock('hf_exchange_swsetpot_batch',2)

       end if

       !qoh: Find position of \phi_{A,a} function wrt simulation cell and fftbox
       first_a_ngwf_idx = ngwf_basis%first_on_atom(atoma)

       call basis_location_func_wrt_cell(aa_cell_start1, aa_cell_start2, &
            aa_cell_start3, ngwf_basis%all_tbs(first_a_ngwf_idx))

       call basis_ket_start_wrt_fftbox(aa_start1,aa_start2,aa_start3, &
            pub_fftbox%total_pt1, pub_fftbox%total_pt2, pub_fftbox%total_pt3)

       batch_count = 0
       ! qoh: Receive relevant NGWFs
       recv_node: do nodeshift=1, pub_total_num_nodes
          node = modulo(pub_my_node_id+nodeshift,pub_total_num_nodes)
          rem_dd: do remote_dd=1,ngwf_basis%num_on_node(node)
             inrplan: if (dd_recv_plan(remote_dd,node)) then
                recv_dd = ngwf_basis%first_on_node(node) + remote_dd -1
                dd_on_grid_buffer=0.0_DP

                if ( node == pub_my_node_id )  then
                   local_dd = recv_dd - ngwf_basis%first_on_node(node) + 1
                else
                   !qoh: Receive and unpack from buffer
                   call comms_recv(node, pub_buffer_sphere%ppd_list, &
                        2*ngwf_basis%n_ppds_sphere(recv_dd), tag=recv_dd)
                   n_dd_ppds = ngwf_basis%n_ppds_sphere(recv_dd)
                   pub_buffer_sphere%n_ppds_sphere = n_dd_ppds
                   dd_npts = n_dd_ppds * pub_cell%n_pts
                   pub_buffer_sphere%offset = 1
                   !qoh: Receive NGWFs from node
                   call comms_recv(node, dd_on_grid_buffer, &
                        dd_npts, tag=recv_dd)
                end if

                !qoh: Find position of \phi_{D,d} wrt simulation cell
                call basis_location_func_wrt_cell( &
                     dd_cell_start1, dd_cell_start2, dd_cell_start3, &
                     ngwf_basis%all_tbs(recv_dd))

                !qoh: Find where to deposit \phi_{D,d} in fftbox
                call basis_location_fb_wrt_box( &
                     dd_start1, dd_start2, dd_start3, &
                     aa_start1, aa_start2, aa_start3, &
                     aa_cell_start1, aa_cell_start2, &
                     aa_cell_start3, dd_cell_start1, &
                     dd_cell_start2, dd_cell_start3, &
                     pub_cell%total_pt1, pub_cell%total_pt2, pub_cell%total_pt3)

                batch_count = batch_count + 1

                if (modulo(batch_count,2) == 1) then

                   dd1_batchidx = minloc(dds_in_batch, &
                        mask = dds_in_batch==recv_dd)
                   dd1_fftbox =0.0_DP
                   dd1_fftbox = 0.0_DP
                   dd2_batchidx(1) = -1

                   !qoh: Use new basis_copy function to transfer \phi_{D,d}
                   !qoh: to an FFT box
                   if (node == pub_my_node_id) then
                      call basis_copy_function_to_box(dd1_fftbox, &
                           pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                           pub_fftbox%total_pt3, &
                           dd_start1, dd_start2, dd_start3, &
                           ngwf_basis%all_tbs(recv_dd),&
                           ngwfs_on_grid, ngwf_basis%spheres(local_dd))
                   else
                      call basis_copy_function_to_box(dd1_fftbox, &
                           pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                           pub_fftbox%total_pt3, &
                           dd_start1, dd_start2, dd_start3,&
                           ngwf_basis%all_tbs(recv_dd),&
                           dd_on_grid_buffer, pub_buffer_sphere)
                   end if

                   if (batch_count /= current_batch_size) cycle rem_dd
                else
                   dd2_batchidx = minloc(dds_in_batch, &
                        mask = dds_in_batch==recv_dd)
                   dd2_fftbox = 0.0_DP
                   dd2_fftbox = 0.0_DP
                   !qoh: use new basis_copy function to transfer \phi_{D,d} to
                   !qoh: an FFT box
                   if (node == pub_my_node_id) then
                      call basis_copy_function_to_box(dd2_fftbox, &
                           pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                           pub_fftbox%total_pt3, &
                           dd_start1,dd_start2,dd_start3,&
                           ngwf_basis%all_tbs(recv_dd),&
                           ngwfs_on_grid, ngwf_basis%spheres(local_dd))
                   else
                      call basis_copy_function_to_box(dd2_fftbox, &
                           pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                           pub_fftbox%total_pt3, &
                           dd_start1,dd_start2,dd_start3,&
                           ngwf_basis%all_tbs(recv_dd),&
                           dd_on_grid_buffer, pub_buffer_sphere)
                   end if
                end if

                do is =1 ,pub_cell%num_spins
                   pot_in_fftbox(centrea%tbstart(1):centrea%tbend(1),&
                        centrea%tbstart(2):centrea%tbend(2),&
                        centrea%tbstart(3):centrea%tbend(3),1) &
                        = tb_batch(:,:,:,dd1_batchidx(1),is)

                   if (modulo(batch_count,2) == 0) &
                        pot_in_fftbox(&
                        centrea%tbstart(1):centrea%tbend(1),&
                        centrea%tbstart(2):centrea%tbend(2),&
                        centrea%tbstart(3):centrea%tbend(3),2) &
                        = tb_batch(:,:,:,dd2_batchidx(1),is)

                   ket_in_fftbox(:,:,:,is) = ket_in_fftbox(:,:,:,is) &
                        + (dd1_fftbox * pot_in_fftbox(:,:,:,1))
                   if (modulo(batch_count,2) == 0 ) &
                        ket_in_fftbox(:,:,:,is) = &
                        ket_in_fftbox(:,:,:,is)&
                        + (dd2_fftbox * pot_in_fftbox(:,:,:,2))
                end do

             end if inrplan
          end do rem_dd
       end do recv_node

    end if null_atoma
    !qoh: Wait for the other nodes to catch up!
    call timer_clock('hf_exchange_comms_wait',1)
    call comms_free
    call timer_clock('hf_exchange_comms_wait',2)

    call utils_trace_out('hf_exchange_tbpot_ngwf_ket')

  end subroutine hf_exchange_tbpot_ngwf_ket

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_swsetpot_batch(tb_batch, swcoeff, centrea,&
       centreb, sphbessels, batch_size, current_batch_size, atoma, atomb)

    !==========================================================================!
    ! This subroutine computes the (exchange) potential in tightbox that       !
    ! results from a fitted set of spherical waves.                            !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in late 2009 and early 2010.                     !
    !==========================================================================!

    use constants, only: DP
    use geometry, only: POINT, operator(*), operator(-), operator(+)
    use simulation_cell, only: pub_cell, pub_fftbox, pub_maxtight_pts1, &
         pub_maxtight_pts2, pub_maxtight_pts3
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    integer,       intent(in)  :: batch_size
    real(kind=DP), intent(out) :: tb_batch(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3,batch_size, pub_cell%num_spins)
    integer,       intent(in)  :: current_batch_size
    integer,       intent(in)  :: atoma
    integer,       intent(in)  :: atomb
    ! qoh: Spherical Bessel functions
    type(bessel),    intent(in) :: sphbessels(pub_cell%num_species,&
         max_num_bessels)
    type(atom_centre), intent(in) :: centrea
    type(atom_centre), intent(in) :: centreb
    real(kind=DP),   intent(in) :: swcoeff(max_sw_set_size,batch_size,&
         pub_cell%num_spins)

    type(point) :: curpoint
    type(point) :: a1, a2, a3
    real(kind=DP) :: coulomb_cutoff
    integer :: ierr
    integer :: d1idx, d2idx, d3idx
    integer :: theoffset

    real(kind=DP), allocatable :: sw_pot(:,:)

    call utils_trace_in('hf_exchange_swsetpot_batch')

    allocate(sw_pot(current_batch_size,pub_cell%num_spins), stat=ierr)
    call utils_alloc_check('hf_exchange_swsetpot_batch','swpot',ierr)
    sw_pot = 0.0_DP

    coulomb_cutoff = 0.499_DP*min(&
         (real(pub_cell%total_pt1,kind=DP)*pub_cell%d1),&
         (real(pub_cell%total_pt2,kind=DP)*pub_cell%d2),&
         (real(pub_cell%total_pt3,kind=DP)*pub_cell%d3))

    a1 = pub_fftbox%d1 * pub_fftbox%a1_unit
    a2 = pub_fftbox%d2 * pub_fftbox%a2_unit
    a3 = pub_fftbox%d3 * pub_fftbox%a3_unit
    curpoint = centrea%incell - centrea%intb - a1 &
         - a2 - a3

    if (atoma == atomb .or. useoverlapmetric .or. singlecentre) then
       do d1idx=1,pub_maxtight_pts1
          curpoint = curpoint + a1
          do d2idx=1,pub_maxtight_pts2
             curpoint = curpoint + a2
             do d3idx=1,pub_maxtight_pts3
                curpoint = curpoint + a3
                call internal_swpot_point(sw_pot,curpoint,centreb,0)
                tb_batch(d1idx,d2idx,d3idx,1:current_batch_size,:) = sw_pot
             end do
             curpoint = curpoint &
                  - real(pub_maxtight_pts3, kind=DP)*a3
          end do
          curpoint = curpoint &
               - real(pub_maxtight_pts2,kind=DP) * a2
       end do

    else
       theoffset = hf_exchange_sw_set_size(sphbessels,centrea)
       do d1idx=1,pub_maxtight_pts1
          curpoint = curpoint + a1
          do d2idx=1,pub_maxtight_pts2
             curpoint = curpoint + a2
             do d3idx=1,pub_maxtight_pts3
                curpoint = curpoint + a3
                call internal_swpot_point(sw_pot,curpoint,centrea,0)
                tb_batch(d1idx,d2idx,d3idx,1:current_batch_size,:) = sw_pot
                call internal_swpot_point(sw_pot,curpoint,centreb,theoffset)
                tb_batch(d1idx,d2idx,d3idx,1:current_batch_size,:) = &
                     tb_batch(d1idx,d2idx,d3idx,1:current_batch_size,:)&
                     + sw_pot
             end do
             curpoint = curpoint &
                  - real(pub_maxtight_pts3, kind=DP)*a3
          end do
          curpoint = curpoint &
               - real(pub_maxtight_pts2,kind=DP) * a2
       end do

    end if

    deallocate(sw_pot, stat=ierr)
    call utils_dealloc_check('hf_exchange_swsetpot_batch','swpot',ierr)

    call utils_trace_out('hf_exchange_swsetpot_batch')

  contains

    subroutine internal_swpot_point(sw_pot,curpoint,centrex, swioffset)

      use constants, only: DP, PI
      use geometry, only: point, operator(*)
      use spherical_wave, only: sw_real_sph_harm_unit
      implicit none
      real(kind=DP),  intent(out) :: sw_pot(:,:)
      type(point),     intent(in) :: curpoint
      type(atom_centre), intent(in) :: centrex
      integer,         intent(in) :: swioffset ! SW index offset
      integer :: swsidx    ! SW index
      real(kind=DP) :: invrad, rad
      integer :: lval
      integer :: species_number
      integer :: sbidx, mval
      real(kind=DP) :: factor, bessint
      type(POINT) :: disp, unit_disp
      real(kind=DP) :: common_factor

      sw_pot = 0.0_DP
      swsidx = swioffset
      species_number = centrex%species_number
      common_factor = 4.0_DP * PI

      call hf_exchange_find_disp_mic(rad, disp, curpoint, centrex%incell, &
           coulomb_cutoff)

      !qoh: Find pointwise contribution to integral
      badr: if (rad > epsilon(1.0_DP)) then
         invrad = 1.0_DP / rad
         unit_disp = invrad * disp

         sb: do sbidx=1,  max_num_bessels
            lval = sphbessels(species_number,sbidx)%lval
            if (lval == -1) exit

            bessint = hf_exchange_sph_bess_pot_int(rad,invrad,&
                 sphbessels(species_number,sbidx),.false.) ! @@
            factor = common_factor / real(((2*lval+1) * (-1)**lval),kind=DP)
            mloop: do mval=-lval,lval
               swsidx = swsidx + 1
               sw_pot = sw_pot + bessint * factor * &
                    swcoeff(swsidx,1:current_batch_size,:) * &
                    sw_real_sph_harm_unit(unit_disp%x,unit_disp%y,unit_disp%z,&
                    lval,mval)
            end do mloop
         end do sb
      else if (rad > -epsilon(1.0_DP)) then badr
         !qoh: At rad = 0 only consider l=0 (we have 0/0 for other values of l)
         sb2: do sbidx=1,  max_num_bessels
            if (sphbessels(species_number,sbidx)%lval /= 0) exit
            swsidx = swsidx + 1
            bessint = 1.0_DP / (sphbessels(species_number,sbidx)%qval**2) &
                 + sphbessels(species_number,sbidx)%nearpotint
            ! z_00 = 0.5_DP/sqrt(pi) = 0.282094791773878_DP
            sw_pot = sw_pot + bessint * common_factor * 0.282094791773878_DP *&
                 swcoeff(swsidx,1:current_batch_size,:)
         end do sb2
      end if badr

    end subroutine internal_swpot_point

  end subroutine hf_exchange_swsetpot_batch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_swpot_batch(tb_batch, es_metric, centrea,&
       centreb, sphbessels, batch_size, current_batch_size, &
       current_sw_set_size, spins, vmat, shifttopoint)

    !==========================================================================!
    ! This subroutine evaluates the electrostatic integrals that result from   !
    ! an NGWF product or spherical wave interacting with a set of spherical    !
    ! waves.                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in late 2009 and early 2010.                     !
    !==========================================================================!

    use constants, only: DP, PI
    use geometry, only: POINT, operator(*), operator(-), operator(+), &
         geometry_distance
    use simulation_cell, only: pub_cell, pub_fftbox, pub_maxtight_pts1, &
         pub_maxtight_pts2, pub_maxtight_pts3
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    integer,        intent(in) :: batch_size
    integer,        intent(in) :: spins ! Set to one for the vmatrix
    real(kind=DP),  intent(in) :: tb_batch(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3,batch_size, spins)
    integer,        intent(in) :: current_batch_size
    integer,        intent(in) :: current_sw_set_size
    ! qoh: Spherical Bessel functions
    type(bessel),   intent(in) :: sphbessels(pub_cell%num_species,&
         max_num_bessels)
    type(atom_centre),intent(in) :: centrea
    type(atom_centre),intent(in) :: centreb
    real(kind=DP), intent(out) :: es_metric(current_sw_set_size,&
         current_batch_size, spins)
    logical,        intent(in) :: vmat
    type(POINT), optional, intent(in) :: shifttopoint

    type(point) :: curpoint
    type(point) :: a1, a2, a3
    real(kind=DP) :: coulomb_cutoff
    real(kind=DP) :: common_factor
    integer :: d1idx, d2idx, d3idx
    integer :: theoffset

    call utils_trace_in('hf_exchange_swpot_batch')

    es_metric = 0.0_DP

    if (vmat) then
       common_factor = 4.0_DP * PI
    else
       common_factor = 4.0_DP * PI * pub_fftbox%weight
    end if

    coulomb_cutoff = 0.499_DP*min(&
         (real(pub_cell%total_pt1,kind=DP)*pub_cell%d1),&
         (real(pub_cell%total_pt2,kind=DP)*pub_cell%d2),&
         (real(pub_cell%total_pt3,kind=DP)*pub_cell%d3))

    a1 = pub_fftbox%d1 * pub_fftbox%a1_unit
    a2 = pub_fftbox%d2 * pub_fftbox%a2_unit
    a3 = pub_fftbox%d3 * pub_fftbox%a3_unit
    curpoint = centreb%incell - centreb%intb - a1 &
         - a2 - a3

    if (present(shifttopoint)) curpoint = curpoint - shifttopoint

    if (geometry_distance(centrea%incell,centreb%incell) &
         < epsilon(1.0_DP) .or. vmat) then
       do d1idx=1,pub_maxtight_pts1
          curpoint = curpoint + a1
          do d2idx=1,pub_maxtight_pts2
             curpoint = curpoint + a2
             do d3idx=1,pub_maxtight_pts3
                curpoint = curpoint + a3
                call internal_swpot_point(curpoint,centrea,0)
             end do
             curpoint = curpoint &
                  - real(pub_maxtight_pts3, kind=DP)*a3
          end do
          curpoint = curpoint &
               - real(pub_maxtight_pts2,kind=DP) * a2
       end do

    else if (singlecentre) then
       do d1idx=1,pub_maxtight_pts1
          curpoint = curpoint + a1
          do d2idx=1,pub_maxtight_pts2
             curpoint = curpoint + a2
             do d3idx=1,pub_maxtight_pts3
                curpoint = curpoint + a3
                call internal_swpot_point(curpoint,centreb,0)
             end do
             curpoint = curpoint &
                  - real(pub_maxtight_pts3, kind=DP)*a3
          end do
          curpoint = curpoint &
               - real(pub_maxtight_pts2,kind=DP) * a2
       end do

    else
       theoffset = hf_exchange_sw_set_size(sphbessels,centrea)
       do d1idx=1,pub_maxtight_pts1
          curpoint = curpoint + a1
          do d2idx=1,pub_maxtight_pts2
             curpoint = curpoint + a2
             do d3idx=1,pub_maxtight_pts3
                curpoint = curpoint + a3
                call internal_swpot_point(curpoint,centrea,0)
                call internal_swpot_point(curpoint,centreb,theoffset)
             end do
             curpoint = curpoint &
                  - real(pub_maxtight_pts3, kind=DP)*a3
          end do
          curpoint = curpoint &
               - real(pub_maxtight_pts2,kind=DP) * a2
       end do

    end if

    call utils_trace_out('hf_exchange_swpot_batch')

  contains

    subroutine internal_swpot_point(curpoint,centrex, swioffset)

      use geometry, only: point, operator(*)
      use spherical_wave, only: sw_real_sph_harm_unit
      implicit none

      type(point),     intent(in) :: curpoint
      type(atom_centre), intent(in) :: centrex
      integer,         intent(in) :: swioffset ! SW index offset
      integer :: swsidx    ! SW index
      real(kind=DP) :: invrad, rad
      integer :: lval
      integer :: species_number
      integer :: sbidx, mval
      real(kind=DP) :: factor, bessint
      type(POINT) :: disp, unit_disp

      swsidx = swioffset
      species_number = centrex%species_number

      call hf_exchange_find_disp_mic(rad, disp, curpoint, centrex%incell, &
           coulomb_cutoff)

      !qoh: Find pointwise contribution to integral
      badr: if (rad > epsilon(1.0_DP)) then
         invrad = 1.0_DP / rad
         unit_disp = invrad * disp

         sb: do sbidx=1,  max_num_bessels
            lval = sphbessels(species_number,sbidx)%lval
            if (lval == -1 .or. swsidx > current_sw_set_size) exit

            bessint = hf_exchange_sph_bess_pot_int(rad,invrad,&
                 sphbessels(species_number,sbidx),.false.) ! @@
            factor = common_factor / real(((2*lval+1) * (-1)**lval),kind=DP)
            mloop: do mval=-lval,lval
               swsidx = swsidx + 1
               if (swsidx > current_sw_set_size) exit
               es_metric(swsidx,:,:) = es_metric(swsidx,:,:) + &
                    tb_batch(d1idx,d2idx,d3idx,1:current_batch_size,:) * &
                    bessint * factor * sw_real_sph_harm_unit&
                    (unit_disp%x,unit_disp%y,unit_disp%z,lval,mval)
            end do mloop
         end do sb
      else if (rad > -epsilon(1.0_DP)) then badr
         !qoh: At rad = 0 only consider l=0 (we have 0/0 for other values of l)
         sb2: do sbidx=1,  max_num_bessels
            if (sphbessels(species_number,sbidx)%lval /= 0) exit
            swsidx = swsidx + 1
            if (swsidx > current_sw_set_size) exit
            bessint = 1.0_DP / (sphbessels(species_number,sbidx)%qval**2) &
                 + sphbessels(species_number,sbidx)%nearpotint
            ! z_00 = 0.5_DP/sqrt(pi) = 0.282094791773878_DP
            es_metric(swsidx,:,:) = es_metric(swsidx,:,:) + &
                 tb_batch(d1idx,d2idx,d3idx,1:current_batch_size,:) * &
                 common_factor * 0.282094791773878_DP * bessint
         end do sb2
      end if badr

    end subroutine internal_swpot_point

  end subroutine hf_exchange_swpot_batch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_swpot_batch_finer(tb_batch, es_metric, centrea,&
       centreb, sphbessels, batch_size, current_batch_size, &
       current_sw_set_size, spins, vmat, shifttopoint)

    !==========================================================================!
    ! This subroutine evaluates the electrostatic integrals that result from   !
    ! an NGWF product or spherical wave interacting with a set of spherical    !
    ! waves. It supports the use of a finer than normal grid.                  !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in  early 2010.                                  !
    !==========================================================================!

    use constants, only: DP, PI
    use geometry, only: POINT, operator(*), operator(-), operator(+), &
         geometry_distance
    use simulation_cell, only: pub_cell, pub_fftbox, pub_maxtight_pts1, &
         pub_maxtight_pts2, pub_maxtight_pts3
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    integer,         intent(in)  :: batch_size
    integer,         intent(in)  :: spins ! Set to one for the vmatrix
    real(kind=DP),   intent(in)  :: tb_batch(pub_maxtight_pts1*finer,&
         pub_maxtight_pts2*finer,pub_maxtight_pts3*finer,batch_size, spins)
    integer,         intent(in)  :: current_batch_size
    integer,         intent(in)  :: current_sw_set_size
    ! qoh: Spherical Bessel functions
    type(bessel),    intent(in)  :: sphbessels(pub_cell%num_species,&
         max_num_bessels)
    type(atom_centre), intent(in)  :: centrea
    type(atom_centre), intent(in)  :: centreb
    real(kind=DP),   intent(out) :: es_metric(current_sw_set_size,&
         current_batch_size, spins)
    logical,         intent(in)  :: vmat
    type(POINT),     intent(in), optional  :: shifttopoint


    integer, save :: ncall=0

    type(point) :: curpoint
    type(point) :: a1, a2, a3
    real(kind=DP) :: coulomb_cutoff
    real(kind=DP) :: common_factor
    integer :: d1idx, d2idx, d3idx
    integer :: theoffset

    call utils_trace_in('hf_exchange_swpot_batch_finer')

    es_metric = 0.0_DP

    if (vmat) then
       common_factor = 4.0_DP * PI
    else
       common_factor = 4.0_DP * PI * pub_fftbox%weight
    end if

    coulomb_cutoff = 0.499_DP*min(&
         (real(pub_cell%total_pt1,kind=DP)*pub_cell%d1),&
         (real(pub_cell%total_pt2,kind=DP)*pub_cell%d2),&
         (real(pub_cell%total_pt3,kind=DP)*pub_cell%d3))

    a1 = pub_fftbox%d1/real(finer,kind=DP) * pub_fftbox%a1_unit
    a2 = pub_fftbox%d2/real(finer,kind=DP) * pub_fftbox%a2_unit
    a3 = pub_fftbox%d3/real(finer,kind=DP) * pub_fftbox%a3_unit

    curpoint = centreb%incell - centreb%intb - a1 &
         - a2 - a3

    ! @jd: made optional, like in hf_exchange_swpot_batch
    if (present(shifttopoint)) curpoint = curpoint - shifttopoint

!    write(*,*) '@@gd: ',geometry_distance(centrea%incell,centreb%incell)

!    write(*,*) '@ncall: ',ncall

    if (geometry_distance(centrea%incell,centreb%incell) &
         < epsilon(1.0_DP) .or. vmat) then
!       write(*,*) '@@if1'
       ! jd: Main branch
       do d1idx=1,pub_maxtight_pts1*finer
          curpoint = curpoint + a1
          do d2idx=1,pub_maxtight_pts2*finer
             curpoint = curpoint + a2
             do d3idx=1,pub_maxtight_pts3*finer
                curpoint = curpoint + a3
!                if(d1idx == 29 .and. d2idx == 37 .and. d3idx == 39) then
!                   write(*,*) '@curpoint is ',curpoint%X,' ',curpoint%Y,' ',curpoint%Z
!                   stop
!                end if

!                if(ncall==4) then
!                  write(*,*) "@cp: ", curpoint%X," ",curpoint%Y," ",curpoint%Z," "
!                  write(*,*) "@ca: ", centrea%incell%X," ",centrea%incell%Y," ",centrea%incell%Z," "
!                  write(*,*) "@cb: ", centreb%incell%X," ",centreb%incell%Y," ",centreb%incell%Z," "
!                end if
                call internal_swpot_point(curpoint,centrea,0,ncall)
             end do
             curpoint = curpoint &
                  - real(pub_maxtight_pts3*finer, kind=DP)*a3
          end do
          curpoint = curpoint &
               - real(pub_maxtight_pts2*finer,kind=DP) * a2
       end do

    else
!       write(*,*) '@@else'
       theoffset = hf_exchange_sw_set_size(sphbessels,centrea)
       do d1idx=1,pub_maxtight_pts1*finer
          curpoint = curpoint + a1
          do d2idx=1,pub_maxtight_pts2*finer
             curpoint = curpoint + a2
             do d3idx=1,pub_maxtight_pts3*finer
                curpoint = curpoint + a3
                call internal_swpot_point(curpoint,centrea,0,ncall)
                call internal_swpot_point(curpoint,centreb,theoffset,ncall)
             end do
             curpoint = curpoint &
                  - real(pub_maxtight_pts3*finer, kind=DP)*a3
          end do
          curpoint = curpoint &
               - real(pub_maxtight_pts2*finer,kind=DP) * a2
       end do

    end if

    ncall = ncall+1

    call utils_trace_out('hf_exchange_swpot_batch_finer')

  contains

    subroutine internal_swpot_point(curpoint,centrex, swioffset,ncall)

      use geometry, only: point, operator(*)
      use spherical_wave, only: sw_real_sph_harm_unit
      implicit none

      type(point),     intent(in) :: curpoint
      type(atom_centre), intent(in) :: centrex
      integer,         intent(in) :: swioffset ! SW index offset
      integer, intent(in) :: ncall !@@

      ! jd: ncall is number of batches already processed
      !     thus ncall=4 means we're now processing indices (pos[]) 41-50


      integer :: swsidx    ! SW index
      real(kind=DP) :: invrad, rad
      integer :: lval
      integer :: species_number
      integer :: sbidx, mval
      real(kind=DP) :: factor, bessint
      type(POINT) :: disp, unit_disp


      swsidx = swioffset
      species_number = centrex%species_number

      call hf_exchange_find_disp_mic(rad, disp, curpoint, centrex%incell, &
           coulomb_cutoff)

      !write(*,*) '@@lmq: ',l_l,' ',m_m,' ',q_q

!      write(*,*) '@tbb: ',tb_batch(29,37,39,1,1)
!      STOP

      !qoh: Find pointwise contribution to integral
      badr: if (rad > epsilon(1.0_DP)) then
         invrad = 1.0_DP / rad
         unit_disp = invrad * disp

         sb: do sbidx=1,  max_num_bessels
            lval = sphbessels(species_number,sbidx)%lval
            if (lval == -1 .or. swsidx > current_sw_set_size) exit

            bessint = hf_exchange_sph_bess_pot_int(rad,invrad,&
                 sphbessels(species_number,sbidx), rad<4.0_DP) ! @@
            factor = common_factor / real(((2*lval+1) * (-1)**lval),kind=DP)

! @@@@@@@@@@@@@@@@@@@@@@@@@@2
!            factor = common_factor / real(((2*lval+1)),kind=DP)

            mloop: do mval=-lval,lval
               swsidx = swsidx + 1
               if (swsidx > current_sw_set_size) exit
! bylo: ncall 4
#if 0
               if(lval==3 .and. mval==0 .and. ncall==4. .and.&
                    abs(sphbessels(species_number,sbidx)%qval-2.417660)<0.0001_DP ) then
                  write(*,*) '@@tb_batch: ',tb_batch(d1idx,d2idx,d3idx,1,1)
                  write(*,*) '@@bessint: ', bessint
                  write(*,*) '@@factor: ',factor
                  write(*,*) '@@harmonic: ', sw_real_sph_harm_unit&
                       (unit_disp%x,unit_disp%y,unit_disp%z,lval,mval)
                  write(*,*) '@@l: ',lval
                  write(*,*) '@@m: ',mval
                  write(*,*) '@@disp: ',disp%x,' ',disp%y,' ',disp%z
                  write(*,*) '@@mult2: ',bessint * factor * sw_real_sph_harm_unit&
                    (unit_disp%x,unit_disp%y,unit_disp%z,lval,mval)
                  write(*,'(a,i0,a,i0,a,f15.10,a,f15.10,a,f15.10,a,f20.15,a,f20.15,a,f20.15,a,f20.15,a,f20.15)') &
                    '@@# ',lval,' ',mval,' ',disp%x,' ',disp%y,' ',disp%z,' ',bessint * factor * sw_real_sph_harm_unit&
                    (unit_disp%x,unit_disp%y,unit_disp%z,lval,mval),' ',&
                    sphbessels(species_number,sbidx)%qval,' ',bessint,' ',sw_real_sph_harm_unit&
                    (unit_disp%x,unit_disp%y,unit_disp%z,lval,mval),' ',tb_batch(d1idx,d2idx,d3idx,1,1)
                  write(*,*) '@@contrib : ',tb_batch(d1idx,d2idx,d3idx,1,1) * &
                    bessint * factor * sw_real_sph_harm_unit&
                    (unit_disp%x,unit_disp%y,unit_disp%z,lval,mval)
               end if
#endif

               es_metric(swsidx,:,:) = es_metric(swsidx,:,:) + &
                    tb_batch(d1idx,d2idx,d3idx,1:current_batch_size,:) * &
                    bessint * factor * sw_real_sph_harm_unit&
                    (unit_disp%x,unit_disp%y,unit_disp%z,lval,mval)

            end do mloop
         end do sb
      else if (rad > -epsilon(1.0_DP)) then badr
         !qoh: At rad = 0 only consider l=0 (we have 0/0 for other values of l)
         sb2: do sbidx=1,  max_num_bessels
            if (sphbessels(species_number,sbidx)%lval /= 0) exit
            swsidx = swsidx + 1
            if (swsidx > current_sw_set_size) exit
            bessint = 1.0_DP / (sphbessels(species_number,sbidx)%qval**2) &
                 + sphbessels(species_number,sbidx)%nearpotint
            ! z_00 = 0.5_DP/sqrt(pi) = 0.282094791773878_DP
            es_metric(swsidx,:,:) = es_metric(swsidx,:,:) + &
                 tb_batch(d1idx,d2idx,d3idx,1:current_batch_size,:) * &
                 common_factor * 0.282094791773878_DP * bessint
         end do sb2
      end if badr

    end subroutine internal_swpot_point

  end subroutine hf_exchange_swpot_batch_finer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_swpot_batch_sg(es_metric, centrea,&
       centreb, sphbessels, current_batch_size, &
       current_sw_set_size)

    !==========================================================================!
    ! This subroutine evaluates the electrostatic integrals that result from   !
    ! a spherical wave interacting with a set of spherical waves.  This is     !
    ! done on a spherical grid (in spherical polar co-ordinates) using an      !
    ! Euler-Maclaurin quadrature for the angular integration and a uniform     !
    ! radial grid with an 11 point Newton-Cotes integration scheme for         !
    ! the radial integrals.                                                    !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in early 2010.                                   !
    !==========================================================================!

    use constants, only: DP, PI
    use geometry, only: POINT, operator(*), operator(-), operator(+), &
         unit_vector
    use rundat, only: pub_hfx_radial_segments, pub_hfx_angular_segments
    use spherical_wave, only: sw_bessel_fast, sw_real_sph_harm_unit
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! qoh: arguments
    integer,         intent(in)  :: current_batch_size
    integer,         intent(in)  :: current_sw_set_size
    ! qoh: Spherical Bessel functions
    type(bessel),    intent(in)  :: sphbessels(pub_cell%num_species,&
         max_num_bessels)
    type(atom_centre), intent(in)  :: centrea
    type(atom_centre), intent(in)  :: centreb
    real(kind=DP),   intent(out) :: es_metric(current_sw_set_size,&
         current_batch_size, 1)

    ! qoh: local variables
    real(kind=DP), allocatable :: sw_batch(:)
    type(point) :: curpoint
    type(point) :: unit_rad
    real(kind=DP) :: coulomb_cutoff
    real(kind=DP) :: common_factor
    real(kind=DP) :: theta
    real(kind=DP) :: phi
    integer :: bspecies
    integer :: swb
    integer :: ierr
    integer :: rpt
    integer :: sbidxb
    integer :: lvalb
    integer :: mvalb
    integer :: phipt
    integer :: thetapt
    integer :: max_rpt
    real(kind=DP) :: avalb
    real(kind=DP), allocatable :: r_points(:,:)
    real(kind=DP) :: sphbess
    real(kind=DP) :: qrsq
    real(kind=DP) :: weight_int
    real(kind=DP) :: theta_der
    real(kind=DP) :: q_theta
    real(kind=DP) :: quad_pts

    call utils_trace_in('hf_exchange_swpot_batch_sg')

    allocate(sw_batch(current_batch_size), stat=ierr)
    call utils_alloc_check('hf_exchange_swpot_batch_sg','sw_batch',ierr)

    sw_batch = 0.0_DP
    es_metric = 0.0_DP
    common_factor = 4.0_DP * PI

    coulomb_cutoff = 0.499_DP*min(&
         (real(pub_cell%total_pt1,kind=DP)*pub_cell%d1),&
         (real(pub_cell%total_pt2,kind=DP)*pub_cell%d2),&
         (real(pub_cell%total_pt3,kind=DP)*pub_cell%d3))

    quad_pts = real(pub_hfx_radial_segments,kind=DP)*2.0_DP*&
         real(pub_hfx_angular_segments,kind=DP)**2
    bspecies=centreb%species_number
    avalb = sphbessels(bspecies,1)%aval
    max_rpt = pub_hfx_radial_segments-1
    allocate(r_points(max_rpt,2), stat=ierr)
    call utils_alloc_check('hf_exchange_swpot_batch_sg','r_points',ierr)

    do rpt= 1, max_rpt
       r_points(rpt,1) = (real(rpt,kind=DP))*avalb/&
            real(pub_hfx_radial_segments,kind=DP)
       ! qoh: 11 point Newton-Cotes integration scheme
       select case (mod(rpt,10))
       case (0)
          r_points(rpt,2) = avalb*16067.0_DP*10.0_DP/299376.0_DP
       case (1,9)
          r_points(rpt,2) = avalb*106300.0_DP*5.0_DP/299376.0_DP
       case (2,8)
          r_points(rpt,2) = avalb*(-48525.0_DP)*5.0_DP/299376.0_DP
       case (3,7)
          r_points(rpt,2) = avalb*272400.0_DP*5.0_DP/299376.0_DP
       case (4,6)
          r_points(rpt,2) = avalb*(-260550.0_DP)*5.0_DP/299376.0_DP
       case (5)
          r_points(rpt,2) = avalb*427368.0_DP*5.0_DP/299376.0_DP
       end select
    end do

    phis: do phipt=1,2*pub_hfx_angular_segments
       phi = PI*phipt/real(pub_hfx_angular_segments,kind=DP)

       thetas: do thetapt=1, pub_hfx_angular_segments
          q_theta =  PI * real(thetapt,kind=DP) / &
               real(pub_hfx_angular_segments,kind=DP)
          theta = q_theta**2 * (3.0_DP*PI - 2.0_DP*q_theta) / (PI**2)
          theta_der = 6.0_DP*q_theta*(PI - q_theta) / (PI**2)
          unit_rad%x = sin(theta)*cos(phi)
          unit_rad%y = sin(theta)*sin(phi)
          unit_rad%z = cos(theta)

          radpts: do rpt=1,max_rpt
             curpoint = centreb%incell + r_points(rpt,1)*unit_rad
             weight_int = sin(theta)*r_points(rpt,1)**2 * r_points(rpt,2) * &
                  theta_der * 2.0_DP* PI*PI/ quad_pts
             swb = 0
             sb: do sbidxb=1,max_num_bessels
                lvalb = sphbessels(bspecies,sbidxb)%lval
                if (lvalb == -1) exit
                qrsq = (r_points(rpt,1)*sphbessels(bspecies,sbidxb)%qval)**2
                sphbess = sw_bessel_fast(lvalb,qrsq) * (-1.0_DP)**lvalb
                m: do mvalb=-lvalb,lvalb
                   swb = swb + 1
                   sw_batch(swb) = weight_int*sphbess*sw_real_sph_harm_unit&
                        (unit_rad%x,unit_rad%y,unit_rad%z,lvalb,mvalb)
                end do m
             end do sb
             call internal_swpot_point(curpoint,centrea,0)

          end do radpts

       end do thetas

    end do phis

    deallocate(r_points, stat=ierr)
    call utils_dealloc_check('hf_exchange_swpot_batch_sg','r_points',ierr)
    deallocate(sw_batch, stat=ierr)
    call utils_dealloc_check('hf_exchange_swpot_batch_sg','sw_batch',ierr)

    call utils_trace_out('hf_exchange_swpot_batch_sg')

  contains

    subroutine internal_swpot_point(curpoint,centrex, swioffset)

      use geometry, only: point, operator(*)
      implicit none

      type(point),     intent(in) :: curpoint
      type(atom_centre), intent(in) :: centrex
      integer,         intent(in) :: swioffset ! SW index offset
      integer :: swsidx    ! SW index
      real(kind=DP) :: invrad, rad
      integer :: lval
      integer :: species_number
      integer :: sbidx, mval
      real(kind=DP) :: factor, bessint
      type(POINT) :: disp, unit_disp

      swsidx = swioffset
      species_number = centrex%species_number

      call hf_exchange_find_disp_mic(rad, disp, curpoint, centrex%incell, &
           coulomb_cutoff)

      !qoh: Find pointwise contribution to integral
      badr: if (rad > epsilon(1.0_DP)) then
         invrad = 1.0_DP / rad
         unit_disp = invrad*disp

         sb: do sbidx=1,  max_num_bessels
            lval = sphbessels(species_number,sbidx)%lval
            if (lval == -1 .or. swsidx > current_sw_set_size) exit

            bessint = hf_exchange_sph_bess_pot_int(rad,invrad,&
                 sphbessels(species_number,sbidx),.false.)
            factor = common_factor / real(((2*lval+1) * (-1)**lval),kind=DP)

! @@@@@@@@@@@@@@@@@@@@@@@@@@2
!            factor = common_factor / real(((2*lval+1)),kind=DP)


            mloop: do mval=-lval,lval
               swsidx = swsidx + 1
               if (swsidx > current_sw_set_size) exit
               es_metric(swsidx,:,1) = es_metric(swsidx,:,1) + &
                    sw_batch * &
                    bessint * factor * sw_real_sph_harm_unit&
                    (unit_disp%x,unit_disp%y,unit_disp%z,lval,mval)
            end do mloop
         end do sb
      else if (rad > -epsilon(1.0_DP)) then badr
         !qoh: At rad = 0 only consider l=0 (we have 0/0 for other values of l)
         sb2: do sbidx=1,  max_num_bessels
            if (sphbessels(species_number,sbidx)%lval /= 0) exit
            swsidx = swsidx + 1
            if (swsidx > current_sw_set_size) exit
            bessint = 1.0_DP / (sphbessels(species_number,sbidx)%qval**2) &
                 + sphbessels(species_number,sbidx)%nearpotint
            ! z_00 = 0.5_DP/sqrt(pi) = 0.282094791773878_DP
            es_metric(swsidx,:,1) = es_metric(swsidx,:,1) + &
                 sw_batch * &
                 common_factor * 0.282094791773878_DP * bessint
         end do sb2
      end if badr

    end subroutine internal_swpot_point

  end subroutine hf_exchange_swpot_batch_sg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_eval_sw_at_nodes(values, nodes, l, m, q, a)

    !==========================================================================!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2012.                               !
    !==========================================================================!

    use chebyshev_rep, only: SPHERE_NODES
    use geometry, only: geometry_magnitude, unit_vector
    use spherical_wave, only: sw_bessel_fast, sw_real_sph_harm_unit

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: values(:) ! Values for all xnodes of this SW
    type(SPHERE_NODES), intent(in) :: nodes
    integer, intent(in) :: l, m
    real(kind=DP), intent(in) :: q, a

    ! jd: Local variables
    character(len=*), parameter :: myself = 'hf_exchange_eval_sw_at_nodes'
    integer :: xi, yi, zi
    integer :: n_stripes
    integer :: loc
    real(kind=DP) :: x, y, z, r
    real(kind=DP) :: sphbess
    type(POINT) :: pos, unit_pos

    ! --------------------------------------------------------------------------
    call utils_trace_in(myself)

    n_stripes = nodes%sph_ranges%n_stripes

    do zi = 1, n_stripes
       pos%z = nodes%znodes(zi)

       do yi = 1, n_stripes
          pos%y = nodes%ynodes(yi,zi)

          do xi = 1, n_stripes
             pos%x = nodes%xnodes(xi,yi,zi)

             unit_pos = unit_vector(pos) ! @optimize: unit_vector calculates magnitude already
             r = geometry_magnitude(pos)

             loc = xi + (yi-1) * n_stripes + (zi-1) * n_stripes*n_stripes

             sphbess = sw_bessel_fast(l,r*r*q*q) * (-1.0_DP)**l

             values(loc) = sphbess * &
                  sw_real_sph_harm_unit(unit_pos%x, unit_pos%y, unit_pos%z, l,m)

          end do
       end do
    end do

    call utils_trace_out(myself)

  end subroutine hf_exchange_eval_sw_at_nodes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_eval_swpot_at_nodes(values, nodes, l, m, q, a, disp, &
       species_number, bessel_idx, sphbessels)

    !==========================================================================!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2012.                               !
    !==========================================================================!

    use chebyshev_rep, only: SPHERE_NODES
    use constants, only: PI
    use geometry, only: geometry_magnitude, unit_vector, operator(-)
    use spherical_wave, only: sw_bessel_fast, sw_real_sph_harm_unit
    use utils, only: utils_assert

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)     :: values(:) ! Values for all xnodes of this SWpot
    type(SPHERE_NODES), intent(in) :: nodes
    integer, intent(in)            :: l, m
    real(kind=DP), intent(in)      :: q, a
    type(POINT), intent(in)        :: disp
    integer, intent(in)            :: species_number
    integer, intent(in)            :: bessel_idx
    type(bessel), intent(in)       :: sphbessels(:,:) !Spherical bessel data

    ! jd: Local variables
    character(len=*), parameter :: myself = 'hf_exchange_eval_swpot_at_nodes'
    integer :: xi, yi, zi
    integer :: n_stripes
    integer :: loc
    real(kind=DP) :: x, y, z, r
    real(kind=DP) :: bessint
    type(POINT) :: pos, vec, unit_vec
    real(kind=DP) :: factor

    ! --------------------------------------------------------------------------
    call utils_trace_in(myself)

    n_stripes = nodes%sph_ranges%n_stripes

    factor = 4.0_DP * PI / real(((2*l+1) * (-1)**l),kind=DP)

    ! write(*,'(a,f14.10,a,f14.10,a,f14.10)') 'disp: ',disp%X,' ',disp%Y,' ',disp%Z

    do zi = 1, n_stripes
       pos%z = nodes%znodes(zi)

       do yi = 1, n_stripes
          pos%y = nodes%ynodes(yi,zi)

          do xi = 1, n_stripes
             pos%x = nodes%xnodes(xi,yi,zi)

             vec = pos - disp

             ! write(*,'(a,f14.10,a,f14.10,a,f14.10)') 'pos: ',pos%X,' ',pos%Y,' ',pos%Z
             ! write(*,'(a,f14.10,a,f14.10,a,f14.10)') 'vec: ',vec%X,' ',vec%Y,' ',vec%Z

             unit_vec = unit_vector(vec) ! @optimize: unit_vector calculates magnitude already
             r = geometry_magnitude(vec)

             call utils_assert(r /= 0.0_DP,'r is unexpectedly zero')

             loc = xi + (yi-1) * n_stripes + (zi-1) * n_stripes*n_stripes

             bessint = hf_exchange_sph_bess_pot_int(r,1.0_DP/r, &
                 sphbessels(species_number, bessel_idx), .false.)

             values(loc) = bessint * factor * &
                  sw_real_sph_harm_unit(unit_vec%x, unit_vec%y, unit_vec%z, l,m)

             ! write(*,*) 'val: ',values(loc)

          end do
       end do
    end do

    call utils_trace_out(myself)

  end subroutine hf_exchange_eval_swpot_at_nodes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_fill_vmatrix(elements, full_vmatrix, ngwf_basis)

    !==========================================================================!
    ! This subroutine fills the V matrix (the electrostatic metric matrix) at  !
    ! the start of the calculation.                                            !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in late 2009 and early 2010.                     !
    !==========================================================================!

    use chebyshev_rep, only: SPHERE_NODES, SPHERE_COEFFS, &
         cheb_alloc_coeffs, cheb_dealloc_coeffs, &
         cheb_gen_nodes_for_sphere, cheb_expansion_for_sphere, &
         cheb_int_product_sphere
    use comms, only: pub_my_node_id, comms_barrier, pub_on_root, &
         pub_total_num_nodes, comms_send, comms_wait, comms_abort
    use constants, only: DP, NORMAL, stdout
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(+)
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_first_atom_on_node, pub_num_atoms_on_node
    use rundat, only: density_batch_size, &
         pub_hfx_cheb_intervals, pub_hfx_cheb_order, &
         pub_hfx_cheb_a_batchsize, pub_hfx_cheb_b_batchsize, pub_rootname, &
         pub_output_detail, pub_hfx_write_vmatrix, pub_hfx_read_vmatrix
    use simulation_cell, only: pub_cell, pub_maxtight_pts1, pub_maxtight_pts2, &
         pub_maxtight_pts3
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index, &
         sparse_create, sparse_put_block, sparse_axpy, sparse_destroy, &
         sparse_transpose, sparse_scale, sparse_show_matrix, &
         sparse_write, sparse_read
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
        utils_abort

    implicit none

    !qoh: Arguments
    type(SPAM3),intent(inout) :: full_vmatrix
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    type(FUNC_BASIS), intent(in)  :: ngwf_basis

    real(kind=DP), allocatable :: vatomblock(:,:) ! Atomblock of v matrix
    type(bessel),  allocatable :: sphbessels(:,:) !Spherical bessel data
    real(kind=DP), allocatable :: radtable(:) ! Radius for each species
    real(kind=DP), allocatable :: es_metric(:,:) ! Electrostatic metric vector
    real(kind=DP), allocatable :: tb_batch(:,:,:,:,:) ! Tightboxes
    complex(kind=DP), allocatable :: tb_zwork(:,:,:)  ! Complex workspace

    !qoh: Local Variables
    integer :: v_idxlen ! length of index of denskern
    integer, allocatable, dimension(:) :: v_idx  ! Denskern index
    integer :: local_a ! Local index of atom A
    integer :: atoma ! Global index of atom A
    integer :: atomb ! Global index of atom B
    !qoh: Column of atom A in density kernel
    integer :: v_fstnzbridx_cola ! first kernel non-zero block row index
    integer :: v_lstnzbridx_cola !last kernel non-zero block row index
    integer :: v_nzbridx_cola ! kernel non-zero block row index
    integer :: batch_size ! Size of the batch
    integer :: ierr !error code

    type(atom_centre) :: centrea, centreb ! Centre data
    type(SPAM3) :: tfull_vmatrix            ! Transposed V matrix
    type(POINT) :: shifttopoint             ! Shift SW centre to point

    ! jd:
    character(len=*), parameter :: myself = 'hf_exchange_fill_vmatrix'
    integer :: species
    real(kind=DP) :: rad
    type(SPHERE_NODES) :: sph_nodes_template
    integer :: sw_idx, sw_idx_loc
    integer :: bessel_idx
    integer :: n_bessels, n_sws
    integer :: l, m
    real(kind=DP) :: q, a
    ! coeffs of SWs belonging to this node
    type(SPHERE_COEFFS), allocatable :: my_sws_sph_coeffs(:) ! sw_idx
    ! coeffs of a batch of SWpots (b,x)
    type(SPHERE_COEFFS), allocatable :: swpot_batch_sph_coeffs(:) ! sw_idx
    ! coeffs of a batch of SWs (a,y)
    type(SPHERE_COEFFS), allocatable :: sw_batch_sph_coeffs(:) ! sw_idx
    real(kind=DP), allocatable :: sph_values(:) ! indexed by xnode_idx
    integer :: n_stripes
    real(kind=DP) :: t1, t2 !@@remove
    integer :: my_first_sw, n_sws_of_mine
    logical :: who_is_done(0:pub_total_num_nodes-1)
    integer :: sw_first(0:pub_total_num_nodes-1)
    integer :: sw_count(0:pub_total_num_nodes-1)
    integer :: sw_where(1:max_sw_set_size)
    integer :: sw_idx_to_bessel_idx(1:max_sw_set_size)
    integer :: sw_idx_to_m(1:max_sw_set_size)
    integer :: i
    logical :: done
    integer :: node
    integer :: send_handle
    character(len=512) :: filename
    logical :: fileexists

    ! -----------------------------------------------------------------------
    if (useoverlapmetric .or. usenpa .or. allowfftonly .or. vmatrix_filled) &
         return

    call utils_trace_in(myself)

    call timer_clock(myself,1)

    full_vmatrix%structure = 'M'
    call sparse_create(full_vmatrix)

    ! jd: Read vmatrix from a file rather than computing it, if asked to
    if(pub_hfx_read_vmatrix) then
       write(filename,'(2a)') trim(pub_rootname),'.vmatrix'
       if (pub_on_root) then
          write(stdout,'(/3a)',advance='no') &
               'Reading vmatrix from file "', trim(filename),'" ...'
       end if

       ! Check that the file exists
       ! ndmh: only root node needs to be able to see the file
       if (pub_on_root) then
          inquire(file=filename,exist=fileexists)
       else
          fileexists = .true.
       end if

       if (fileexists) then
          ! Read density kernel from this file
          call sparse_read(full_vmatrix,trim(filename))
       else
          if (pub_on_root) write(stdout,'(/a/)') ' File not found, quitting.'
          call comms_abort
       end if

       if (pub_on_root) write(stdout,'(a)') ' done'

       return

    end if

    who_is_done(:) = .false.

    batch_size = density_batch_size

    !qoh: Get index array for v matrix
    v_idxlen = sparse_index_length(full_vmatrix)
    allocate(v_idx(v_idxlen),stat=ierr)
    call utils_alloc_check(myself,'v_idx',ierr)
    call sparse_generate_index(v_idx,full_vmatrix)
    if (.not. sph_grid_metric .and. .not. singlecentre) then
       if (recip_grid_metric) then
          allocate(tb_batch(pub_maxtight_pts1,pub_maxtight_pts2,&
               pub_maxtight_pts3,batch_size,1),stat=ierr)
          call utils_alloc_check(myself,'tb_batch',ierr)

          allocate(tb_zwork(pub_maxtight_pts1,pub_maxtight_pts2,&
               pub_maxtight_pts3),stat=ierr)
          call utils_alloc_check(myself,'tb_zwork',ierr)

       else
          allocate(tb_batch(pub_maxtight_pts1*finer,pub_maxtight_pts2*finer,&
               pub_maxtight_pts3*finer,batch_size,1),stat=ierr)
          call utils_alloc_check(myself,'tb_batch',ierr)
       end if
    end if

    allocate(radtable(pub_cell%num_species),stat=ierr)
    call utils_alloc_check(myself,'radtable',ierr)
    ! qoh: Initialise max_num_bessels and max_num_sphwaves
    call hf_exchange_num_sph_functions(radtable)
    allocate(sphbessels(pub_cell%num_species,max_num_bessels),stat=ierr)
    call utils_alloc_check(myself,'sphbessels',ierr)
    call hf_exchange_sph_bessels_init(sphbessels,radtable)

    allocate(vatomblock(max_sw_set_size,max_sw_set_size), stat=ierr)
    call utils_dealloc_check(myself,'vatomblock',ierr)
    allocate(es_metric(max_sw_set_size,batch_size),stat=ierr)
    call utils_alloc_check('hf_exchange_fillvmatrix','es_metric',ierr)

    ! @jd: Generalize later
    call utils_assert(minval(radtable) == maxval(radtable), &
         "Current implementation of Hartree-Fock exchange only works if all &
         &NGWF radii are identical, sorry.")
    rad = radtable(1)

    species = 1 ! @ generalize this later

    ! jd: Count bessels, spherical waves
    n_sws = 0
    n_bessels = 0
    do bessel_idx = 1, max_num_bessels
       l = sphbessels(species,bessel_idx)%lval
       if(l == -1) exit
       n_bessels = n_bessels + 1
       n_sws = n_sws + (2*l+1) ! jd: account for m = -l..l
    end do
    call utils_assert(n_bessels > 0,'n_bessels must be positive')
    call utils_assert(n_sws > 0,'n_sws must be positive')

    call hf_exchange_memory_estimate(pub_hfx_cheb_intervals, &
         pub_hfx_cheb_order, pub_hfx_cheb_a_batchsize, &
         pub_hfx_cheb_b_batchsize, n_sws)

    ! jd: Generate Chebyshev nodes
    call cheb_gen_nodes_for_sphere(sph_nodes_template, rad, &
         pub_hfx_cheb_intervals, pub_hfx_cheb_order)

    n_stripes = sph_nodes_template%sph_ranges%n_stripes

    ! jd: Distribute SW numbers across nodes
    call hf_exchange_distribute_sws(sw_first,sw_count,sw_where,n_sws)
    my_first_sw = sw_first(pub_my_node_id)
    n_sws_of_mine = sw_count(pub_my_node_id)

    allocate(my_sws_sph_coeffs(n_sws_of_mine), stat=ierr)
    call utils_alloc_check(myself,'my_sws_sph_coeffs',ierr)
    allocate(sph_values(sph_nodes_template%n_points_total), stat=ierr)
    call utils_alloc_check(myself,'sph_values',ierr)
    allocate(swpot_batch_sph_coeffs(pub_hfx_cheb_b_batchsize), stat=ierr)
    call utils_alloc_check(myself,'swpot_batch_sph_coeffs',ierr)
    allocate(sw_batch_sph_coeffs(pub_hfx_cheb_a_batchsize), stat=ierr)
    call utils_alloc_check(myself,'sw_batch_sph_coeffs',ierr)

    ! jd: Allocate memory for coefficients (members of swpot_batch_sph_coeffs)
    do i = 1, pub_hfx_cheb_a_batchsize
       call cheb_alloc_coeffs(sw_batch_sph_coeffs(i), &
            pub_hfx_cheb_intervals, pub_hfx_cheb_order)
    end do
    do i = 1, pub_hfx_cheb_b_batchsize
       call cheb_alloc_coeffs(swpot_batch_sph_coeffs(i), &
            pub_hfx_cheb_intervals, pub_hfx_cheb_order)
    end do

    ! jd: Fill sw_idx_to_bessel_idx and sw_idx_to_m
    sw_idx = 1
    do bessel_idx = 1, n_bessels
       l = sphbessels(species,bessel_idx)%lval
       do m = -l, l ! @ optimize: most of evaluation is independent of m
          sw_idx_to_bessel_idx(sw_idx) = bessel_idx
          sw_idx_to_m(sw_idx) = m
          sw_idx = sw_idx + 1
       end do
    end do

    ! jd: Go over all local SWs and evaluate them at Chebyshev nodes, &
    !     find out the expansion coefficients
    sw_idx = 1
    do bessel_idx = 1, n_bessels
       l = sphbessels(species,bessel_idx)%lval

       q = sphbessels(species,bessel_idx)%qval
       a = sphbessels(species,bessel_idx)%aval

       do m = -l, l ! @ optimize: most of evaluation is independent of m

          if(sw_idx >= my_first_sw .and. &
              sw_idx < my_first_sw + n_sws_of_mine) then

             sw_idx_loc = sw_idx - my_first_sw + 1

             write(*,*) 'Node ',pub_my_node_id,' expanding SW #',sw_idx,' (',sw_idx_loc,')'

             call hf_exchange_eval_sw_at_nodes(sph_values(:), &
                  sph_nodes_template, l, m, q, a)

             call cheb_alloc_coeffs(my_sws_sph_coeffs(sw_idx_loc), &
                  pub_hfx_cheb_intervals, pub_hfx_cheb_order)

             call cheb_expansion_for_sphere(my_sws_sph_coeffs(sw_idx_loc), &
                  sph_nodes_template, sph_values(:))
          end if ! of if SW is mine

          sw_idx = sw_idx + 1

       end do ! jd: over m
    end do ! jd: over bessels

    ! jd: Fill vmatrix

    ! Over all atoms local to this node
    loop_A: do local_a=1,pub_num_atoms_on_node(pub_my_node_id)
       atoma = pub_first_atom_on_node(pub_my_node_id) + local_a -1

       call hf_exchange_init_centre(centrea, ngwf_basis, elements, atoma, &
            shifttopoint)

       !qoh: Get kernel and overlap blocks in column of atom A
       v_fstnzbridx_cola = v_idx(local_a)        ! first kernel non-zero block row index
       v_lstnzbridx_cola = v_idx(local_a+1) - 1  ! last  kernel non-zero block row index

       !qoh: Loop over atoms (atom B) that have non zero blocks K_{AB}
       loop_B: do v_nzbridx_cola = v_fstnzbridx_cola, v_lstnzbridx_cola
          atomb = v_idx(v_nzbridx_cola)
          vatomblock = 0.0_DP

          if ((atoma .ge. atomb .and. mod(atoma,2) == mod(atomb,2)) .or. &
               ( atoma .lt. atomb .and. mod(atoma,2) /= mod(atomb,2))) then
             call hf_exchange_init_centre(centreb, ngwf_basis, elements, atomb)

             if (atoma == atomb) then
                call hf_exchange_vmatrixblock_sc(vatomblock,centrea,&
                     sphbessels)
             else if (singlecentre) then
                vatomblock = 0.0_DP
             else if (sph_grid_metric) then
                call internal_vmatrixblock_sg(centreb,centrea)
             else if (chebyshev_grid_metric) then
                call internal_vmatrixblock_cheb(centreb,centrea,atomb,atoma)
             else
                call internal_vmatrixblock(centreb,centrea,atomb,atoma)
             end if

             !qoh: Insert atom block in to sparse matrix structure.
             ! jd: NB 'atoma' must belong to this node
             call sparse_put_block(vatomblock, full_vmatrix, atomb, atoma)

          end if

       end do loop_B

    end do loop_A

    ! jd: After we're done, first notify everyone, then keep serving
    !     other nodes with SWs, until everyone is done
    do node = 0, pub_total_num_nodes-1
       call comms_send(node, -1, 1, tag = SW_REQUEST_TAG, &
            return_handle = send_handle, add_to_stack = .false.)
       call comms_wait(send_handle) !@podejrzane - blokuje az wszyscy odbiora
    end do

    write(*,*) '### Node ',pub_my_node_id,' done computing and will now serve'
    done = .false.
    do while(.not. done)
       call hf_exchange_serve_SW_coeffs(done, my_sws_sph_coeffs, &
            sw_first, sw_count, sw_where, who_is_done)
    end do

    call comms_barrier

    deallocate(sph_values,stat=ierr)
    call utils_dealloc_check(myself,'sph_values',ierr)

    ! jd: Deallocate memory for batches of coefficients
    do i = 1, pub_hfx_cheb_a_batchsize
       call cheb_dealloc_coeffs(sw_batch_sph_coeffs(i))
    end do
    do i = 1, pub_hfx_cheb_b_batchsize
       call cheb_dealloc_coeffs(swpot_batch_sph_coeffs(i))
    end do

    deallocate(my_sws_sph_coeffs, stat=ierr)
    call utils_dealloc_check(myself,'my_sws_sph_coeffs',ierr)
    deallocate(swpot_batch_sph_coeffs, stat=ierr)
    call utils_dealloc_check(myself,'swpot_batch_sph_coeffs',ierr)
    deallocate(sw_batch_sph_coeffs, stat=ierr)
    call utils_dealloc_check(myself,'sw_batch_sph_coeffs',ierr)

    deallocate(es_metric,stat=ierr)
    call utils_dealloc_check('hf_exchange_fillvmatrix','es_metric',ierr)
    deallocate(radtable,stat=ierr)
    call utils_dealloc_check(myself,'radtable',ierr)
    deallocate(sphbessels,stat=ierr)
    call utils_dealloc_check(myself,'sphbessels',ierr)
    if (.not. sph_grid_metric .and. .not. singlecentre) then
       deallocate(tb_batch,stat=ierr)
       call utils_dealloc_check(myself,'tb_batch',ierr)
       if(recip_grid_metric) then
          deallocate(tb_zwork,stat=ierr)
          call utils_dealloc_check(myself,'tb_zwork',ierr)
       end if
    end if
    deallocate(v_idx, stat=ierr)
    call utils_dealloc_check(myself,'v_idx',ierr)

    !qoh: Fill in missing blocks
    call comms_barrier
    vatomblock = 0.0_DP
    call sparse_create(tfull_vmatrix,full_vmatrix)
    !call sparse_scale(full_vmatrix,0.5_DP)
    call sparse_transpose(tfull_vmatrix,full_vmatrix)
    !qoh: Zero diagonal atomblocks in tfull_vmatrix
    do local_a=1,pub_num_atoms_on_node(pub_my_node_id)
       atoma = pub_first_atom_on_node(pub_my_node_id) + local_a -1
       call sparse_put_block(vatomblock,tfull_vmatrix,atoma,atoma)
    end do
    call sparse_axpy(full_vmatrix,tfull_vmatrix,1.0_DP)
    if (printvmat .and. pub_on_root) then
       print *, ""
       print *, "-VMAT:-=============================@"
    end if
    if (printvmat) call sparse_show_matrix(full_vmatrix)
    if (printvmat .and. pub_on_root) then
       print *, ""
       print *, "-------------------------------------"
    end if

    call sparse_destroy(tfull_vmatrix)

    deallocate(vatomblock, stat=ierr)
    call utils_dealloc_check(myself,'vatomblock',ierr)

    ! jd: Write vmatrix to a file, if asked to
    if(pub_hfx_write_vmatrix) then
       write(filename,'(2a)') trim(pub_rootname),'.vmatrix'
       if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
          write(stdout,'(/3a)',advance='no') &
               'Writing vmatrix to file "', trim(filename),'" ...'
       end if

       call sparse_write(full_vmatrix,trim(filename))

       if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
          write(stdout,'(a)') ' done'
       end if
    end if

    vmatrix_filled = .true.

    call utils_trace_out(myself)
    call timer_clock(myself,2)

  contains

    subroutine internal_vmatrixblock_cheb(centrex,centrey,atomx,atomy)

      !==========================================================================!
      !--------------------------------------------------------------------------!
      !==========================================================================!

      use geometry, only: operator(-)
      use rundat, only: pub_hfx_cheb_a_batchsize, pub_hfx_cheb_b_batchsize

      implicit none

      ! Arguments
      type(atom_centre), intent(in) :: centrex
      type(atom_centre), intent(in) :: centrey
      integer, intent(in)           :: atomx
      integer, intent(in)           :: atomy

      ! Local variables
      character(len=*), parameter :: myself = 'internal_vmatrixblock_cheb'

      integer :: current_swy
      integer :: sbidx
      integer :: batch_count
      integer :: lval
      integer :: mval
      integer :: n_batches
      integer :: batch
      integer :: current_swy_set_size, current_swx_set_size
      integer :: batch_next_swy
      integer :: first_swy
      integer :: current_batch_size
      integer :: swx, swy
      integer :: mvalx, lvalx, xsbidx
      type(POINT) :: disp

      integer :: n_xbatches, n_ybatches
      integer :: xbatch, ybatch
      integer :: swx_first, swx_last, swx_in_batch
      integer :: swy_first, swy_last, swy_in_batch
      integer :: besselx, bessely
      integer :: lx, mx, ly, my
      real(kind=DP) :: qx, ax, qy, ay

      ! --------------------------------------------------------------------
      call utils_trace_in(myself)
      call utils_assert(atomx /= atomy,'Internal error [1] in '//myself)

      disp = centrex%incell - centrey%incell

      ! jd: Determine the number of batches over x (SWpot)
      n_xbatches = n_sws / pub_hfx_cheb_b_batchsize
      if(mod(n_sws,pub_hfx_cheb_b_batchsize) /= 0) n_xbatches  = n_xbatches + 1

      ! jd: Determine the number of batches over y (SW)
      n_ybatches = n_sws / pub_hfx_cheb_a_batchsize
      if(mod(n_sws,pub_hfx_cheb_a_batchsize) /= 0) n_ybatches  = n_ybatches + 1

      ! ----------------------------------
      ! jd: Loop over batches of SW pots
      ! ----------------------------------
      do xbatch = 1, n_xbatches

         swx_first = (xbatch-1) * pub_hfx_cheb_b_batchsize + 1
         swx_last = swx_first + pub_hfx_cheb_b_batchsize - 1
         if(swx_last > n_sws) swx_last = n_sws

         write(*,'(a,i3,a,i3,a,i3,a,i3)') &
              'b-batch ',xbatch,'/',n_xbatches,': ',swx_first,':',swx_last

         ! ------------------------------
         ! 1. Prepare coeffs for SW pots
         ! ------------------------------

         ! jd: Loop over SW pots in this batch
         do swx = swx_first, swx_last

            swx_in_batch = swx - swx_first + 1

            besselx = sw_idx_to_bessel_idx(swx)
            mx = sw_idx_to_m(swx)

            lx = sphbessels(centrex%species_number,besselx)%lval
            qx = sphbessels(centrex%species_number,besselx)%qval
            ax = sphbessels(centrex%species_number,besselx)%aval

!            write(*,'(a,i0)') 'Evaluating SWpot ', swx
            ! jd: Evaluate this SWpot on x acting on y
            call hf_exchange_eval_swpot_at_nodes(sph_values(:), &
                 sph_nodes_template, lx, mx, qx, ax, disp, &
                 centrex%species_number, besselx, sphbessels)

!            write(*,'(a,i0)') 'Expanding SWpot ', swx
            ! jd: Perform Chebyshev expansion for SWpot on x acting on Y
            call cheb_expansion_for_sphere( &
                 swpot_batch_sph_coeffs(swx_in_batch), &
                 sph_nodes_template, sph_values(:))

         end do

         ! --------------------------
         ! 2. Prepare coeffs for SWs
         ! --------------------------

         ! -----------------------------
         ! jd: Loop over batches of SWs
         ! -----------------------------
         do ybatch = 1, n_ybatches

            swy_first = (ybatch-1) * pub_hfx_cheb_a_batchsize + 1
            swy_last = swy_first + pub_hfx_cheb_a_batchsize - 1
            if(swy_last > n_sws) swy_last = n_sws

!            write(*,'(a,i3,a,i3,a,i3,a,i3)') &
!                 'a-batch ',ybatch,'/',n_ybatches,': ',swy_first,':',swy_last

            ! jd: Loop over SWs in this batch
            do swy = swy_first, swy_last

               swy_in_batch = swy - swy_first + 1

               bessely = sw_idx_to_bessel_idx(swy)
               my = sw_idx_to_m(swy)

               ly = sphbessels(centrey%species_number,besselx)%lval
               qy = sphbessels(centrey%species_number,besselx)%qval
               ay = sphbessels(centrey%species_number,besselx)%aval

               call hf_exchange_obtain_SW_coeffs( &
                    sw_batch_sph_coeffs(swy_in_batch), my_sws_sph_coeffs, &
                    swy, sw_first, sw_count, sw_where, who_is_done)

            end do ! SWs in batch

            ! --------------------------
            ! Deal with this batchblock
            ! --------------------------

            do swx = swx_first, swx_last
               swx_in_batch = swx - swx_first + 1

               do swy = swy_first, swy_last
                  swy_in_batch = swy - swy_first + 1

                  !write(*,*) 'BATCHBLOCK: ', swx,' ',swy

                  vatomblock(swx,swy) = cheb_int_product_sphere( &
                       sw_batch_sph_coeffs(swy_in_batch), &
                       swpot_batch_sph_coeffs(swx_in_batch), &
                       sph_nodes_template)

               end do ! SWs in batch

            end do ! SWpots in batch
 
         end do ! batches of SWs

      end do ! batches of SWpots

      ! jd: Forced consistency of V matrix
      if (force_con) then
         swy = 0
         do sbidx=1, max_num_bessels
            lval = sphbessels(centrey%species_number,sbidx)%lval
            if (lval == -1) exit
            do mval=-lval,lval
               swy = swy + 1
               swx = 0
               do xsbidx=1,max_num_bessels
                  lvalx = sphbessels(centrex%species_number,xsbidx)%lval
                  if (lval == -1) exit
                  do mvalx=-lvalx,lvalx
                     swx = swx + 1
                     if (swy > swx) &
                          vatomblock(swy,swx) = &
                               (-1)**(lval+lvalx)*vatomblock(swx,swy)
                  end do
               end do
            end do
         end do
      end if

      call utils_trace_out(myself)

    end subroutine internal_vmatrixblock_cheb


    subroutine internal_vmatrixblock(centrex,centrey,atomx,atomy)

      !==========================================================================!
      ! This subroutine calculates a V matrix atomblock using the real Cartesian !
      ! grid or the reciprocal Cartesian grid.                                   !
      !--------------------------------------------------------------------------!
      ! Written by Quintin Hill in October 2009.                                 !
      !==========================================================================!

      use spherical_wave, only: sw_recp_generate_in_tb

      type(atom_centre), intent(in) :: centrex
      type(atom_centre), intent(in) :: centrey
      integer,intent(in) :: atomx
      integer,intent(in) :: atomy

      integer :: current_swy
      integer :: sbidx
      integer :: batch_count
      integer :: lval
      integer :: mval
      integer :: n_batches
      integer :: batch
      integer :: current_swy_set_size, current_swx_set_size
      integer :: batch_next_swy
      integer :: first_swy
      integer :: current_batch_size
      integer :: swx, swy
      integer :: mvalx, lvalx, xsbidx
      type(POINT) :: disp

      call utils_trace_in('internal_vmatrixblock')

      current_swx_set_size = hf_exchange_sw_set_size(sphbessels,centrex)
      current_swy_set_size = hf_exchange_sw_set_size(sphbessels,centrey)
      disp = centrey%intb + shifttopoint

      n_batches = current_swy_set_size / batch_size
      if (mod(current_swy_set_size, batch_size) > 0) n_batches = n_batches + 1

      batch_next_swy = 1
      ! Do XY block
      do batch=1,n_batches
         current_swy = 0
         first_swy = batch_next_swy
         batch_count = 0
         do sbidx=1, max_num_bessels
            lval = sphbessels(centrey%species_number,sbidx)%lval
            if (lval == -1) exit
            mloopy: do mval=-lval,lval
               current_swy = current_swy + 1
               if (current_swy /= batch_next_swy) cycle mloopy
               batch_count = batch_count + 1
               if (recip_grid_metric) then
                  call sw_recp_generate_in_tb(&
                       tb_batch(:,:,:,batch_count,1), tb_zwork, lval,mval, &
                       sphbessels(centrey%species_number,sbidx)%qval,&
                       sphbessels(centrey%species_number,sbidx)%aval,disp,&
                       pub_maxtight_pts1,pub_maxtight_pts2,pub_maxtight_pts3)
               else
                  call hf_exchange_sw_tb(&
                       tb_batch(:,:,:,batch_count,1), lval,mval, &
                       sphbessels(centrey%species_number,sbidx)%qval,&
                       sphbessels(centrey%species_number,sbidx)%aval,disp)
!                  write(*,*) '@@b: ',batch_count,' ',lval,' ',mval,' ',sphbessels(centrey%species_number,sbidx)%qval
               end if
               batch_next_swy = current_swy + 1
               if (batch_count == batch_size) exit
            end do mloopy
            if (batch_count == batch_size) exit
         end do
         current_batch_size = batch_count
         if (shavesw) call hf_exchange_shave_tb_batch(tb_batch, centrey, &
              batch_size, sphbessels(centrey%species_number,1)%aval)
         es_metric = 0.0_DP
         if (atomx == atomy .or. force_con) current_swx_set_size = current_swy

         if (recip_grid_metric) then
            call hf_exchange_swpot_batch(tb_batch, &
                 es_metric(1:current_swx_set_size,1:current_batch_size), &
                 centrex, centrey, sphbessels, batch_size, current_batch_size, &
                 current_swx_set_size,1,.true.,shifttopoint)
         else
            call hf_exchange_swpot_batch_finer(tb_batch, &
                 es_metric(1:current_swx_set_size,1:current_batch_size), &
                 centrex, centrey, sphbessels, batch_size, current_batch_size, &
                 current_swx_set_size,1,.true.,shifttopoint)
         end if

         vatomblock(1:current_swx_set_size,first_swy:current_swy) = &
              es_metric(1:current_swx_set_size,&
              1:current_batch_size)

      end do

      if (atomx == atomy) then
         do swy= 1, current_swy_set_size
            do swx = 1,swy-1
               vatomblock(swy,swx) = vatomblock(swx,swy)
            end do
         end do
         swy = 0
         do sbidx=1, max_num_bessels
            lval = sphbessels(centrey%species_number,sbidx)%lval
            if (lval == -1) exit
            do mval=-lval,lval
               swy = swy + 1
               swx = 0
               do xsbidx=1,max_num_bessels
                  lvalx = sphbessels(centrex%species_number,xsbidx)%lval
                  do mvalx=-lvalx,lvalx
                     swx = swx + 1
                     if (mval /= mvalx .or. lval /= lvalx) &
                          vatomblock(swx,swy)= 0.0_DP
                  end do
               end do
            end do
         end do
      end if

      ! jd: Forced consistency of V matrix
      if (force_con) then
         swy = 0
         do sbidx=1, max_num_bessels
            lval = sphbessels(centrey%species_number,sbidx)%lval
            if (lval == -1) exit
            do mval=-lval,lval
               swy = swy + 1
               swx = 0
               do xsbidx=1,max_num_bessels
                  lvalx = sphbessels(centrex%species_number,xsbidx)%lval
                  if (lval == -1) exit
                  do mvalx=-lvalx,lvalx
                     swx = swx + 1
                     if (swy > swx) &
                          vatomblock(swy,swx) = (-1)**(lval+lvalx)*vatomblock(swx,swy)
                  end do
               end do
            end do
         end do
      end if

      call utils_trace_out('internal_vmatrixblock')

    end subroutine internal_vmatrixblock

    subroutine internal_vmatrixblock_sg(centrex,centrey)

      !==========================================================================!
      ! This subroutine calculates a V matrix atomblock on the spherical grid.   !
      !--------------------------------------------------------------------------!
      ! Written by Quintin Hill in March 2010.                                   !
      !==========================================================================!

      ! qoh: Calculate atom block of V matrix on radial grid

      type(atom_centre), intent(in) :: centrex
      type(atom_centre), intent(in) :: centrey
      integer :: current_swy_set_size, current_swx_set_size

      call utils_trace_in('internal_vmatrixblock_sg')

      current_swx_set_size = hf_exchange_sw_set_size(sphbessels,centrex)
      current_swy_set_size = hf_exchange_sw_set_size(sphbessels,centrey)

      call hf_exchange_swpot_batch_sg(&
           vatomblock(1:current_swx_set_size,1:current_swy_set_size), &
           centrex, centrey, sphbessels, current_swy_set_size, &
           current_swx_set_size)

      call utils_trace_out('internal_vmatrixblock_sg')

    end subroutine internal_vmatrixblock_sg

  end subroutine hf_exchange_fill_vmatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_vmatrixblock_sc(vatomblock,centrex,sphbessels)

    !==========================================================================!
    ! This subroutine calculates a same centre V matrix atomblock analytically.!
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in March 2010.                                   !
    !==========================================================================!

    use constants, only: DP, PI
    use simulation_cell, only: pub_cell

    real(kind=DP), intent(out) :: vatomblock(:,:)
    type(atom_centre), intent(in) :: centrex
    type(bessel), intent(in)    :: sphbessels(pub_cell%num_species,&
         max_num_bessels)

    integer :: swx, swy       ! x and y spherical wave index
    integer :: species_number
    integer :: xlval, ylval   ! x and y l value
    integer :: xsbidx, ysbidx ! x and y spherical bessel index
    integer :: xmval, ymval   ! x and y m value
    real(kind=DP) :: aval
    real(kind=DP) :: qval
    real(kind=DP) :: qa
    real(kind=DP) :: swxpart
    real(kind=DP) :: diagel

    call utils_trace_in('hf_exchange_vmatrixblock_sc')

    species_number = centrex%species_number

    vatomblock = 0.0_DP
    swx = 0
    aval = sphbessels(species_number,1)%aval
    xsb: do xsbidx=1,max_num_bessels
       qval = sphbessels(species_number,xsbidx)%qval
       xlval = sphbessels(species_number,xsbidx)%lval
       qa = qval * aval

       ! jd: This likely realizes (5.7.3)

       select case (xlval)
       case(0)
          swxpart = -4.0_DP * PI * aval * cos(qa) / (qval**2)
          diagel = 4.0_DP*PI * aval/(2.0_DP*(qval**4)) &
               +swxpart* sphbessels(species_number,xsbidx)%nearpotint
       case(1)
          swxpart = -4.0_DP*PI/(3.0_DP*qval**4) * sin(qa)*(qa*qa)
          diagel = 4.0_DP*PI/(qval**6) *(qa*qa - sin(qa)**2)/&
               (2.0_DP*aval)&
               + swxpart * sphbessels(species_number,xsbidx)%nearpotint
       case(2)
          swxpart = 4.0_DP*PI/5.0_DP*&
               sphbessels(species_number,xsbidx)%farpotint
          diagel = PI /(qa**3 * qval**5) * (qa*qa*(2.0_DP*qa*qa - 12.0_DP) + &
               (sin(qa))**2 *(qa*qa*(qa*qa*2.0_DP/3.0_DP + 2.0_DP) +12.0_DP))&
               + swxpart * sphbessels(species_number,xsbidx)%nearpotint
       case(3)
          swxpart = 4.0_DP*PI/7.0_DP*&
               sphbessels(species_number,xsbidx)%farpotint
          diagel = 4.0_DP*PI*((2.0_DP*(-45.0_DP - 15.0_DP*qa**2 - &
               6.0_DP*qa**4 + qa**6) + 6.0_DP*(15.0_DP - 25.0_DP*qa**2 + &
               2.0_DP*qa**4)*cos(2.0_DP*qa) + qa*(180.0_DP - 60.0_DP*qa**2 + &
               qa**4)*sin(2.0_DP*qa)))/ (4.0_DP*qa**5*qval**5) + &
               swxpart * sphbessels(species_number,xsbidx)%nearpotint
       case(4)
          swxpart = 4.0_DP*PI/9.0_DP*&
               sphbessels(species_number,xsbidx)%farpotint
          diagel = 4.0_DP*PI*((-3150.0_DP - 630.0_DP*qa**2 &
               - 90.0_DP*qa**4 - 20.0_DP*qa**6 + 2.0_DP*qa**8) &
               + 10.0_DP*(315.0_DP - 567.0_DP*qa**2 + 93.0_DP*qa**4 &
               - 2.0_DP*qa**6)*cos(2.0_DP*qa) &
               + qa*(6300.0_DP - 2940.0_DP*qa**2 + &
               180.0_DP*qa**4 - qa**6)*sin(2.0_DP*qa))/ (4.0_DP*qa**7*qval**5)&
               + swxpart * sphbessels(species_number,xsbidx)%nearpotint
       case default
          swxpart = 0.0_DP
          diagel = 0.0_DP
          exit
       end select
       xm: do xmval=-xlval,xlval
          swx = swx + 1
          vatomblock(swx,swx) = diagel
          swy = 0
          ysb: do ysbidx= 1, max_num_bessels
             ylval = sphbessels(species_number,ysbidx)%lval
             ym: do ymval= -ylval,ylval
                swy = swy + 1
                if (xlval /= ylval .or. xmval /= ymval .or. swx == swy) cycle
                vatomblock(swx,swy) =  swxpart * &
                     sphbessels(species_number,ysbidx)%nearpotint
             end do ym ! jd: m's on y
          end do ysb ! jd: Bessels on y
       end do xm ! jd: m's on x
    end do xsb ! jd: Bessels on x

    call utils_trace_out('hf_exchange_vmatrixblock_sc')

  end subroutine hf_exchange_vmatrixblock_sc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_sw_tb(sw_out, lval, mval, qval, aval, disp)

    !==========================================================================!
    ! This subroutine calculates in real space a spherical wave in a tightbox. !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill on 20/01/2010.                                   !
    !==========================================================================!

    use geometry, only: point, operator(.dot.), operator(-), operator(+), &
         operator(*)
    use simulation_cell, only: pub_fftbox, pub_maxtight_pts1, &
         pub_maxtight_pts2, pub_maxtight_pts3
    use spherical_wave, only: sw_bessel_fast, sw_real_sph_harm
    use timer, only: timer_clock

    implicit none

    real(kind=DP), intent(out) :: sw_out(pub_maxtight_pts1*finer,&
         pub_maxtight_pts2*finer,pub_maxtight_pts3*finer)
    integer,       intent(in) :: lval
    integer,       intent(in) :: mval
    real(kind=DP), intent(in) :: qval ! j_l(qr)
    real(kind=DP), intent(in) :: aval ! Cutoff radius of spherical wave
    type(POINT),   intent(in) :: disp ! Tightbox origin relative to origin of sw

    ! Internal variables

    integer :: d1idx, d2idx, d3idx ! Tightbox grid points
    type(POINT)      :: a1, a2, a3 ! Vectors between grid points
    type(POINT)      :: curpoint ! current point
    real(kind=DP)    :: asq ! a squared
    real(kind=DP)    :: qsq ! q squared
    real(kind=DP)    :: sph_bess
    real(kind=DP)    :: sph_harm
    real(kind=DP)    :: rad, radsq

    call utils_trace_in('hf_exchange_sw_tb')
    call timer_clock('hf_exchange_sw_tb',1)

    sw_out = 0.0_DP
    ! calculate factor common to all points
    qsq = qval*qval
    asq = aval*aval

    a1 = pub_fftbox%d1/real(finer,kind=DP) * pub_fftbox%a1_unit
    a2 = pub_fftbox%d2/real(finer,kind=DP) * pub_fftbox%a2_unit
    a3 = pub_fftbox%d3/real(finer,kind=DP) * pub_fftbox%a3_unit
    curpoint = POINT(0.0_DP,0.0_DP,0.0_DP) - disp - a1 &
         - a2 - a3

    do d1idx=1,pub_maxtight_pts1*finer
       curpoint = curpoint + a1
       do d2idx=1,pub_maxtight_pts2*finer
          curpoint = curpoint + a2
          do d3idx=1,pub_maxtight_pts3*finer
             curpoint = curpoint + a3
             radsq = curpoint .DOT. curpoint

             if (radsq < asq) then
                rad = sqrt(radsq)
                ! calculate spherical harmonic
                sph_harm = sw_real_sph_harm(curpoint%x,curpoint%y,&
                     curpoint%z,rad, lval,mval)
                ! calculate spherical bessel function
                sph_bess = sw_bessel_fast(lval,radsq*qsq)
                ! calculate spherical wave
                sw_out(d1idx,d2idx,d3idx) = sph_harm * sph_bess &
                     * pub_fftbox%weight/real(finer**3,kind=DP) * (-1.0_DP)**lval

!                write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,f15.10,a,f7.3,a,f10.7,a,f10.7,a,f10.7,a,f15.10,a,f15.10,a,f15.10,a,f15.10)') &
!                    '@sw ',d1idx,' ',d2idx,' ',d3idx,' ',lval,' ',mval,' ',qval,' ',aval,' ', &
!                     curpoint%x,' ',curpoint%y,' ',curpoint%z,' ',sph_harm,' ',sph_bess,' ',&
!                     pub_fftbox%weight/real(finer**3,kind=DP),' ',sw_out(d1idx,d2idx,d3idx)
             end if
          end do
          curpoint = curpoint &
               - real(pub_maxtight_pts3*finer, kind=DP)*a3
       end do
       curpoint = curpoint &
            - real(pub_maxtight_pts2*finer,kind=DP) * a2
    end do

    call utils_trace_out('hf_exchange_sw_tb')
    call timer_clock('hf_exchange_sw_tb',2)

  end subroutine hf_exchange_sw_tb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_init_vmatrix(elements)

    !==========================================================================!
    ! This subroutine initialises the structure used to hold all atomblocks of !
    ! the V matrix.                                                            !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill on 29/10/2010.                                   !
    !==========================================================================!

    use comms, only: pub_total_num_nodes
    use constants, only: DP
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_orig_atom, parallel_strategy_distr_funcs
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_init_blocking_scheme, BLKS_SW
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    real(kind=DP), allocatable :: radtable(:)
    integer, allocatable :: first_func_on_node(:)
    integer, allocatable :: num_funcs_on_node(:)
    integer, allocatable ::first_func_on_atom(:)
    integer, allocatable ::num_funcs_on_atom(:)
    integer, allocatable ::node_of_func(:)
    integer, allocatable ::atom_of_func(:)
    integer, allocatable :: num_sw(:)
    integer, allocatable :: nfuncs_orig(:)
    integer :: max_funcs_on_node
    integer :: max_funcs_on_atom
    integer :: ierr      ! Error flag
    integer :: eidx      ! Element idx
    integer :: total_sw  ! Total number of spherical waves

    call utils_trace_in('hf_exchange_init_vmatrix')

    allocate(radtable(pub_cell%num_species), stat=ierr)
    call utils_alloc_check('hf_exchange_init_vmatrix', 'radtable',ierr)
    allocate(num_sw(pub_cell%num_species), stat=ierr)
    call utils_alloc_check('hf_exchange_init_vmatrix', 'num_sw',ierr)
    allocate(nfuncs_orig(pub_cell%nat), stat=ierr)
    call utils_alloc_check('hf_exchange_init_vmatrix', 'nfuncs_orig',ierr)

    call hf_exchange_num_sph_functions(radtable, num_sw)

    do eidx = 1,pub_cell%nat
       nfuncs_orig(eidx) = num_sw(elements(eidx)%species_number)
    end do
    total_sw = sum(nfuncs_orig)
    pub_cell%num_sw = total_sw

    allocate(num_funcs_on_node(0:pub_total_num_nodes-1), stat=ierr)
    call utils_alloc_check('hf_exchange_init_vmatrix', 'num_funcs_on_node',ierr)
    allocate(first_func_on_node(0:pub_total_num_nodes),stat=ierr)
    call utils_alloc_check('hf_exchange_init_vmatrix', 'first_func_on_node',ierr)
    allocate(num_funcs_on_atom(1:pub_cell%nat),stat=ierr)
    call utils_alloc_check('hf_exchange_init_vmatrix', 'num_funcs_on_atom',ierr)
    allocate(first_func_on_atom(1:pub_cell%nat),stat=ierr)
    call utils_alloc_check('hf_exchange_init_vmatrix', 'first_func_on_atom',ierr)
    allocate(node_of_func(1:total_sw),stat=ierr)
    call utils_alloc_check('hf_exchange_init_vmatrix', 'node_of_func',ierr)
    allocate(atom_of_func(1:total_sw),stat=ierr)
    call utils_alloc_check('hf_exchange_init_vmatrix', 'atom_of_func',ierr)

    call parallel_strategy_distr_funcs(total_sw, nfuncs_orig, pub_orig_atom, &
         first_func_on_node, num_funcs_on_node, first_func_on_atom, &
         num_funcs_on_atom, node_of_func, atom_of_func, max_funcs_on_node, &
         max_funcs_on_atom)

    call sparse_init_blocking_scheme(BLKS_SW,total_sw,num_funcs_on_node, &
         num_funcs_on_atom, first_func_on_node, first_func_on_atom, &
         atom_of_func, node_of_func)

    deallocate(num_funcs_on_node, stat=ierr)
    call utils_dealloc_check('hf_exchange_init_vmatrix', 'num_funcs_on_node',&
         ierr)
    deallocate(first_func_on_node,stat=ierr)
    call utils_dealloc_check('hf_exchange_init_vmatrix', 'first_func_on_node',&
         ierr)
    deallocate(num_funcs_on_atom,stat=ierr)
    call utils_dealloc_check('hf_exchange_init_vmatrix', 'num_funcs_on_atom',&
         ierr)
    deallocate(first_func_on_atom,stat=ierr)
    call utils_dealloc_check('hf_exchange_init_vmatrix', 'first_func_on_atom',&
         ierr)
    deallocate(node_of_func,stat=ierr)
    call utils_dealloc_check('hf_exchange_init_vmatrix', 'node_of_func',ierr)
    deallocate(atom_of_func,stat=ierr)
    call utils_dealloc_check('hf_exchange_init_vmatrix', 'atom_of_func',ierr)
    deallocate(radtable, stat=ierr)
    call utils_dealloc_check('hf_exchange_init_vmatrix', 'radtable',ierr)
    deallocate(num_sw, stat=ierr)
    call utils_dealloc_check('hf_exchange_init_vmatrix', 'num_sw',ierr)
    deallocate(nfuncs_orig, stat=ierr)
    call utils_dealloc_check('hf_exchange_init_vmatrix', 'nfuncs_orig',ierr)

    call utils_trace_out('hf_exchange_init_vmatrix')

  end subroutine hf_exchange_init_vmatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_set_pub_hfxsw()

    !========================================!
    ! This subroutine sets pub_hfxsw.        !
    !----------------------------------------!
    ! Written by Quintin Hill on 29/10/2010. !
    !========================================!

    use rundat, only: pub_hfxsw, pub_usehfx
    use simulation_cell, only: pub_fftbox

    implicit none

    logical :: boxiscell

    call utils_trace_in('hf_exchange_set_pub_hfxsw')

    boxiscell = (pub_fftbox%coin1 .and. pub_fftbox%coin2 .and. pub_fftbox%coin3)
    pub_hfxsw = ( (.not. allowfftonly .or. .not. boxiscell) &
         .and. pub_usehfx .and. .not. usenpa )

    call utils_trace_out('hf_exchange_set_pub_hfxsw')

  end subroutine hf_exchange_set_pub_hfxsw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_make_vmatrix(vmatrix, full_vmatrix, &
       current_sw_set_size, sphbessels, centrea, centreb, atoma, atomb)

    !==========================================================================!
    ! This subroutine constructs the V matrix for the current pair of atoms.   !
    ! To avoid communication if atom B is on a different core then the BB      !
    ! atom block is calculated.  Otherwise the atomblocks are retrieved from   !
    ! full_vmatrix.                                                            !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in late 2009 and early 2010.                     !
    !==========================================================================!

    use comms, only: pub_my_node_id
    use constants, only: DP
    use parallel_strategy, only: pub_node_of_atom, pub_orig_atom
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_get_block
    use timer, only: timer_clock

    implicit none

    integer, intent(in)           :: current_sw_set_size
    real(kind=DP), intent(out)    :: vmatrix(current_sw_set_size, &
         current_sw_set_size)
    type(SPAM3), intent(in)       :: full_vmatrix
    type(atom_centre), intent(in) :: centrea
    type(atom_centre), intent(in) :: centreb
    type(bessel), intent(in)      :: sphbessels(pub_cell%num_species,&
         max_num_bessels)
    integer, intent(in)           :: atoma, atomb

    ! Local variables
    integer :: theoffset

    call utils_trace_in('hf_exchange_make_vmatrix')
    call timer_clock('hf_exchange_make_vmatrix',1)

    if (singlecentre) then
       vmatrix = 0.0_DP
       call hf_exchange_vmatrixblock_sc(vmatrix,centreb,sphbessels)
       call utils_trace_out('hf_exchange_make_vmatrix')
       call timer_clock('hf_exchange_make_vmatrix',2)
       return
    end if

    vmatrix = 0.0_DP
    theoffset = hf_exchange_sw_set_size(sphbessels,centrea)
    call sparse_get_block(vmatrix(1:theoffset,1:theoffset),full_vmatrix, &
         atoma, atoma)
    if (atoma /= atomb) then
       call sparse_get_block(vmatrix(theoffset+1:current_sw_set_size,&
            1:theoffset),full_vmatrix, atomb, atoma)
       if (pub_node_of_atom(pub_orig_atom(atomb)) == pub_my_node_id) then
          call sparse_get_block(vmatrix(1:theoffset,&
               theoffset+1:current_sw_set_size),full_vmatrix, atoma, atomb)
          call sparse_get_block(vmatrix(theoffset+1:current_sw_set_size,&
               theoffset+1:current_sw_set_size),full_vmatrix, atomb, atomb)

       else
          vmatrix(1:theoffset,theoffset+1:current_sw_set_size) = &
               transpose(vmatrix(theoffset+1:current_sw_set_size,1:theoffset))
          call hf_exchange_vmatrixblock_sc&
               (vmatrix(theoffset+1:current_sw_set_size,&
               theoffset+1:current_sw_set_size),centreb,sphbessels)
       end if
    end if

    call utils_trace_out('hf_exchange_make_vmatrix')
    call timer_clock('hf_exchange_make_vmatrix',2)

  end subroutine hf_exchange_make_vmatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_shave_tb_batch(tb_batch, atomcentre, batch_size,radius)

    !==========================================================================!
    ! This subroutine sets all points in a batch of tightboxes outside the     !
    ! localisation sphere of a centre to zero.                                 !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill on 20/01/2010.                                   !
    !==========================================================================!

    use geometry, only: POINT, operator(*), operator(-), operator(+), &
         geometry_magnitude
    use simulation_cell, only: pub_fftbox, pub_maxtight_pts1, &
         pub_maxtight_pts2, pub_maxtight_pts3

    implicit none

    integer, intent(in) :: batch_size
    real(kind=DP), intent(inout) ::  tb_batch&
         (pub_maxtight_pts1,pub_maxtight_pts2,pub_maxtight_pts3,&
         batch_size, 1)
    type(atom_centre), intent(in) :: atomcentre
    real(kind=DP), intent(in) :: radius ! cutoff radius
    integer :: d1idx, d2idx, d3idx
    type(point) :: curpoint
    type(point) :: a1, a2, a3

    call utils_trace_in('hf_exchange_shave_tb_batch')

    a1 = pub_fftbox%d1 * pub_fftbox%a1_unit
    a2 = pub_fftbox%d2 * pub_fftbox%a2_unit
    a3 = pub_fftbox%d3 * pub_fftbox%a3_unit
    curpoint = POINT(0.0_DP,0.0_DP,0.0_DP) - atomcentre%intb - a1 &
         - a2 - a3

    do d1idx=1,pub_maxtight_pts1
       curpoint = curpoint + a1
       do d2idx=1,pub_maxtight_pts2
          curpoint = curpoint + a2
          do d3idx=1,pub_maxtight_pts3
             curpoint = curpoint + a3
             if(geometry_magnitude(curpoint) > radius) &
                  tb_batch(d1idx,d2idx,d3idx,:,1) = 0.0_DP
          end do
          curpoint = curpoint &
               - real(pub_maxtight_pts3, kind=DP)*a3
       end do
       curpoint = curpoint &
            - real(pub_maxtight_pts2,kind=DP) * a2
    end do

    call utils_trace_out('hf_exchange_shave_tb_batch')

  end subroutine hf_exchange_shave_tb_batch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_shave_tb(tb, atomcentre, radius)

    !==========================================================================!
    ! This subroutine sets all points in a tightbox outside the localisation   !
    ! sphere of a centre to zero.                                              !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill on 20/01/2010.                                   !
    !==========================================================================!

    use geometry, only: POINT, operator(*), operator(-), operator(+), &
         geometry_magnitude
    use simulation_cell, only: pub_fftbox, pub_maxtight_pts1, &
         pub_maxtight_pts2, pub_maxtight_pts3

    implicit none

    real(kind=DP), intent(inout) ::  tb&
         (pub_maxtight_pts1,pub_maxtight_pts2,pub_maxtight_pts3)
    type(atom_centre), intent(in) :: atomcentre
    real(kind=DP), intent(in) :: radius
    integer :: d1idx, d2idx, d3idx
    type(point) :: curpoint
    type(point) :: a1, a2, a3

    call utils_trace_in('hf_exchange_shave_tb')

    a1 = pub_fftbox%d1 * pub_fftbox%a1_unit
    a2 = pub_fftbox%d2 * pub_fftbox%a2_unit
    a3 = pub_fftbox%d3 * pub_fftbox%a3_unit
    curpoint = POINT(0.0_DP,0.0_DP,0.0_DP) - atomcentre%intb - a1 &
         - a2 - a3

    do d1idx=1,pub_maxtight_pts1
       curpoint = curpoint + a1
       do d2idx=1,pub_maxtight_pts2
          curpoint = curpoint + a2
          do d3idx=1,pub_maxtight_pts3
             curpoint = curpoint + a3
             if(geometry_magnitude(curpoint) > radius) &
                  tb(d1idx,d2idx,d3idx) = 0.0_DP
          end do
          curpoint = curpoint &
               - real(pub_maxtight_pts3, kind=DP)*a3
       end do
       curpoint = curpoint &
            - real(pub_maxtight_pts2,kind=DP) * a2
    end do

    call utils_trace_out('hf_exchange_shave_tb')

  end subroutine hf_exchange_shave_tb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_shift_product(prod_in_fftbox, ngwf_basis, atoma, atomb)

    !==========================================================================!
    ! This subroutine shifts an NGWF product in an FFT box.                    !
    ! Input : Product centred on B in FFT box centred on B.                    !
    ! Output: Product centred on B in FFT box centred on A.                    !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill on 16/02/2010.                                   !
    !==========================================================================!

    use basis, only: basis_ket_start_wrt_fftbox, basis_location_func_wrt_cell, &
         basis_location_fb_wrt_box, basis_extract_function_from_box, &
         basis_copy_function_to_box, basis_copy_sphere, &
         basis_sphere_deallocate, SPHERE
    use constants, only: DP
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: pub_cell, pub_fftbox
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    real(kind=DP), intent(inout) :: prod_in_fftbox(pub_fftbox%total_ld1,&
         pub_fftbox%total_ld2,pub_fftbox%total_pt3,2)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    integer,          intent(in) :: atoma, atomb

    real(kind=DP), allocatable :: prod_on_grid_buffer(:)
    type(SPHERE) :: bsphere
    integer :: aa_start1, aa_start2, aa_start3
    integer :: aa_cell_start1, aa_cell_start2, aa_cell_start3
    integer :: bb_start1, bb_start2, bb_start3
    integer :: bb_cell_start1, bb_cell_start2, bb_cell_start3
    integer :: ierr, idx

    call utils_trace_in('hf_exchange_shift_product')

    allocate(prod_on_grid_buffer&
         (ngwf_basis%max_n_ppds_sphere*pub_cell%n_pts), stat=ierr)
    call utils_alloc_check('hf_exchange_shift_product','prod_on_grid_buffer',&
         ierr)

    call basis_copy_sphere(bsphere,&
         ngwf_basis%spheres(ngwf_basis%first_on_atom(atomb)),1)
    call basis_location_func_wrt_cell(aa_cell_start1, aa_cell_start2, &
         aa_cell_start3, ngwf_basis%all_tbs(ngwf_basis%first_on_atom(atoma)))

    !qoh: Find position of \phi_{A,a} function wrt fftbox
    call basis_ket_start_wrt_fftbox(aa_start1,aa_start2,aa_start3, &
         pub_fftbox%total_pt1, pub_fftbox%total_pt2, pub_fftbox%total_pt3)

    call basis_location_func_wrt_cell( &
         bb_cell_start1, bb_cell_start2, bb_cell_start3, &
         ngwf_basis%all_tbs(ngwf_basis%first_on_atom(atomb)))

    call basis_location_fb_wrt_box( &
         bb_start1, bb_start2, bb_start3, &
         aa_start1, aa_start2, aa_start3, &
         aa_cell_start1, aa_cell_start2, &
         aa_cell_start3, bb_cell_start1, &
         bb_cell_start2, bb_cell_start3, &
         pub_cell%total_pt1, pub_cell%total_pt2, pub_cell%total_pt3)
    do idx=1,2
       prod_on_grid_buffer=0.0_DP

       ! qoh: extract ppds belonging to product function from product fftbox
       call basis_extract_function_from_box(prod_on_grid_buffer, &
            pub_fftbox%total_ld1, pub_fftbox%total_ld2, pub_fftbox%total_pt3, &
            prod_in_fftbox(:,:,:,idx), &
            ngwf_basis%spheres(ngwf_basis%first_on_atom(atomb)), &
            ngwf_basis%all_tbs(ngwf_basis%first_on_atom(atomb)), &
            aa_start1, aa_start2, aa_start3, 1)
       prod_in_fftbox(:,:,:,idx) = 0.0_DP

       call basis_copy_function_to_box(prod_in_fftbox(:,:,:,idx), &
            pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
            pub_fftbox%total_pt3, &
            bb_start1, bb_start2, bb_start3, &
            ngwf_basis%all_tbs(ngwf_basis%first_on_atom(atomb)),&
            prod_on_grid_buffer, bsphere)
    end do

    call basis_sphere_deallocate(bsphere)

    deallocate(prod_on_grid_buffer, stat=ierr)
    call utils_dealloc_check('hf_exchange_shift_product',&
         'prod_on_grid_buffer',ierr)

    call utils_trace_out('hf_exchange_shift_product')

  end subroutine hf_exchange_shift_product

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_find_disp_mic(rad, disp, curpoint, centre_in_cell, &
       coulomb_cutoff)

    !==========================================================================!
    ! This subroutine finds the displacement vector using the minimum image    !
    ! convention under the constraint of a cutoff Coulomb operator.            !
    !--------------------------------------------------------------------------!
    ! Split from internal_swpot_point by Quintin Hill on 29/04/2010.           !
    !==========================================================================!

    use constants, only: DP
    use geometry, only: POINT, operator(+), operator(*), operator(-), &
         geometry_MAGNITUDE
    use simulation_cell, only: pub_cell

    implicit none

    real(kind=DP), intent(out) :: rad  ! Magnitude of displacement
    type(POINT), intent(out) :: disp   ! Displacement vector of chosen point
    type(POINT), intent(in) :: curpoint ! Curent point wrt cell
    type(POINT), intent(in) :: centre_in_cell ! Centre of reference wrt cell
    real(kind=DP), intent(in) :: coulomb_cutoff ! Coulomb cutoff distance

    integer :: a1_neighbour,a2_neighbour,a3_neighbour
    type(POINT) :: periodic_centre
#ifndef CUBECC
    real(kind=DP) :: trial_rad
#endif

    rad = -1.0_DP
    !qoh: loop over periodic images to find appropriate one
    a1: do a1_neighbour = -1,1
       a2: do a2_neighbour = -1,1
          a3: do a3_neighbour = -1,1

             periodic_centre = centre_in_cell &
                  + real(a1_neighbour,kind=DP)*pub_cell%a1 &
                  + real(a2_neighbour,kind=DP)*pub_cell%a2 &
                  + real(a3_neighbour,kind=DP)*pub_cell%a3
             disp = curpoint - periodic_centre
#ifdef CUBECC
             if (max(abs(disp%x),abs(disp%y),abs(disp%z)) < coulomb_cutoff) &
                  then
                rad = geometry_MAGNITUDE(disp)
                exit a1
             end if
#else
             trial_rad = geometry_MAGNITUDE(disp)
             if (trial_rad < coulomb_cutoff) then
                rad = trial_rad
                exit a1
             end if
#endif
          end do a3
       end do a2
    end do a1

  end subroutine hf_exchange_find_disp_mic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function hf_exchange_sph_bess_pot_int(rad,invrad,sphbessel,dbg)

    !==========================================================================!
    ! This function calculates a spherical Bessel potential integral at a      !
    ! point.                                                                   !
    !--------------------------------------------------------------------------!
    ! Split from internal_swpot_point by Quintin Hill on 29/04/2010.           !
    !==========================================================================!

    implicit none

    real(kind=DP)             :: hf_exchange_sph_bess_pot_int
    type(bessel),  intent(in) :: sphbessel ! Spherical Bessel  parameters
    real(kind=DP), intent(in) :: rad       ! Radius from centre
    real(kind=DP), intent(in) :: invrad    ! 1/rad
    logical, intent(in) :: dbg

    real(kind=DP) :: invq
    real(kind=DP) :: q, qr

    if(dbg) then
!       write(*,*) '@@@ rad: ',rad
    end if

    if (rad > sphbessel%aval) then
       hf_exchange_sph_bess_pot_int = invrad**(sphbessel%lval+1) * &
            sphbessel%farpotint ! jd: (5.6.4, downstairs)

!       if(dbg) then
!          write(*,*) '@@@ FAR: ', hf_exchange_sph_bess_pot_int
!          write(*,*) '@@@ invrad: ',invrad
!          write(*,*) '@@@ ^ ',(sphbessel%lval+1)
!          write(*,*) '@@@ prefactor: ',invrad**(sphbessel%lval+1)
!          write(*,*) '@@@ int: ',sphbessel%farpotint
!       end if

    else

!       if(dbg) then
!          write(*,*) '@@@ NEAR'
!       end if

       ! jd: Formulas below verified with Mathematica

       invq = 1.0_DP / sphbessel%qval
       q = sphbessel%qval
       qr = q*rad
       select case(sphbessel%lval)
       case (0)
          hf_exchange_sph_bess_pot_int = invrad * invq**3 * sin(qr) &
               + sphbessel%nearpotint
!          if(dbg) then
!             write(*,*) '@@@ q: ',q
!             write(*,*) '@@@ rad: ',rad
!             write(*,*) '@@@ qr: ',qr
!             write(*,*) '@@@ first term: ',invrad * invq**3 * sin(qr)
!             write(*,*) '@@@ nearpotint: ',sphbessel%nearpotint
!          end if

       case(1)
          hf_exchange_sph_bess_pot_int = 3.0_DP * invrad * invq**3 * &
               (sin(qr)*invq*invrad - cos(qr)) &
               + rad*sphbessel%nearpotint
       case(2)
          hf_exchange_sph_bess_pot_int = invq**5 * invrad**3 *&
               ((15.0_DP - 5.0_DP*qr*qr) * sin(qr) - 15.0_DP*qr*cos(qr))  &
               +rad*rad*sphbessel%nearpotint

       case(3)
          hf_exchange_sph_bess_pot_int = invq**6 * invrad**4 * &
               ((105.0_DP-42.0_DP*qr*qr) * sin(qr) &
               - (105.0_DP - 7.0_DP*qr*qr)*qr*cos(qr)) &
               +rad**3 * sphbessel%nearpotint
       case(4)
          hf_exchange_sph_bess_pot_int = invq**7 * invrad**5 * 9.0_DP *&
               ((10.0_DP*qr*qr-105.0_DP)*qr*cos(qr) + &
               (qr*qr*(qr*qr-45.0_DP) + 105.0_DP)*sin(qr)) &
               + rad**4 * sphbessel%nearpotint
       case default
          hf_exchange_sph_bess_pot_int = 0.0_DP
       end select
    end if

  end function hf_exchange_sph_bess_pot_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_npa_batch_var1(tb_batch, tb_batch_prod, centrea,&
       centreb, batch_size)
    !========================================================================!
    ! Variant 1 of hf_exchange_npa_batch.                                    !
    ! Calculates the potential U_{\beta,\delta} in a tightbox (tb) on atom A,!
    ! by integrating P(r2)/|r1-r2| over tb 2, for every point in tb 1, with  !
    ! batches over \delta.                                                   !
    ! This variant uses a somewhat optimized version of the integration,     !
    ! which changes the ordering of sums, so that the innermost sum is over  !
    ! |r1-r2|, which can then be reused.                                     !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    ! @docme                                                                 !
    !------------------------------------------------------------------------!
    ! Written by Quintin Hill in January 2010                                !
    !========================================================================!

    use constants, only: DP, PI
    use simulation_cell, only: pub_cell, pub_maxtight_pts1, pub_maxtight_pts2, &
         pub_maxtight_pts3

    implicit none

    ! qoh: Arguments

    integer,       intent(in)  :: batch_size
    real(kind=DP), intent(out) :: tb_batch(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3,batch_size, pub_cell%num_spins)
    real(kind=DP), intent(in) :: tb_batch_prod(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3,batch_size, pub_cell%num_spins)
    type(atom_centre), intent(in) :: centrea
    type(atom_centre), intent(in) :: centreb

    ! qoh: Local variables
    integer :: startpnt(3) ! starting point
    integer :: pnt(3)      ! current point
    real(kind=DP) :: d1sq, d2sq, d3sq ! grid spacing squared
    integer :: tbp1, tbp2, tbp3 ! number of tightbox grid points in direction
    integer :: d1idx, d2idx, d3idx
    integer :: invd1idx, invd2idx, invd3idx
    real(kind=DP) :: invrad
    real(kind=DP) :: sincint ! sinc integral

    call utils_trace_in('hf_exchange_npa_batch_var1')

    ! Get useful constants

    tbp1 = pub_maxtight_pts1
    tbp2 = pub_maxtight_pts2
    tbp3 = pub_maxtight_pts3
    d1sq = pub_cell%d1**2
    d2sq = pub_cell%d2**2
    d3sq = pub_cell%d3**2

    ! qoh: Potential at centre of a psinc function due to the psinc function.
    ! qoh: This is calculated by treating the psinc function as an l=0 spherical
    ! qoh: wave with a cutoff radius equal to the grid spacing. (The spherical
    ! qoh: harmonic factor is removed.)
    sincint = 8.0_DP * min(pub_cell%d1,pub_cell%d2,pub_cell%d3)**2 / &
         (PI*pub_cell%weight)

    tb_batch = 0.0_DP
    startpnt = centrea%tbcellstart - centreb%tbcellstart

    d1: do d1idx=1,tbp1

       pnt(1) = d1idx + startpnt(1) - 1
       invd1idx = tbp1 - d1idx + 1

       d2: do d2idx=1,tbp2

          pnt(2) = d2idx + startpnt(2) - 1
          invd2idx = tbp2 - d2idx + 1
          d3: do d3idx=1,tbp3

             pnt(3) = d3idx + startpnt(3) -1
             invd3idx = tbp3 - d3idx + 1

             !qoh: 8 areas to deal with.

             ! Area 1
             if (abs(pnt(1))+abs(pnt(2))+abs(pnt(3)) /= 0) then
                invrad = 1.0_DP/sqrt(d1sq*real((pnt(1))**2,kind=DP) &
                     + d2sq*real((pnt(2))**2,kind=DP) &
                     + d3sq*real((pnt(3))**2,kind=DP))
             else
                invrad = sincint
             end if

             tb_batch(d1idx:tbp1,d2idx:tbp2,d3idx:tbp3,:,:) = &
                  tb_batch(d1idx:tbp1,d2idx:tbp2,d3idx:tbp3,:,:)  &
                  + invrad*tb_batch_prod(1:invd1idx,1:invd2idx,1:invd3idx,:,:)   ! (3) 9.5%

             if (d1idx /= 1) then

                ! Area 2

                if (abs(pnt(1)-tbp1)+abs(pnt(2))+abs(pnt(3)) /= 0) then
                   invrad = 1.0_DP/sqrt(d1sq*real((pnt(1)-tbp1)**2,kind=DP) &
                        + d2sq*real((pnt(2))**2,kind=DP) &
                        + d3sq*real((pnt(3))**2,kind=DP))
                else
                   invrad = sincint
                end if

                tb_batch(1:(d1idx-1),d2idx:tbp2,d3idx:tbp3,:,:) = &
                     tb_batch(1:(d1idx-1),d2idx:tbp2,d3idx:tbp3,:,:)  &
                     + invrad*tb_batch_prod((invd1idx+1):tbp1,&
                     1:invd2idx,1:invd3idx,:,:)

                if (d2idx /= 1) then

                   ! Area 3

                   if (abs(pnt(1)-tbp1)+abs(pnt(2)-tbp1)+abs(pnt(3)) &
                        /= 0)  then
                      invrad = 1.0_DP/sqrt(d1sq*real((pnt(1)-tbp1)**2,kind=DP)&
                           + d2sq*real((pnt(2)-tbp2)**2,kind=DP) &
                           + d3sq*real((pnt(3))**2,kind=DP))
                   else
                      invrad = sincint
                   end if

                   tb_batch(1:(d1idx-1),1:(d2idx-1),d3idx:tbp3,:,:) =&
                        tb_batch(1:(d1idx-1),1:(d2idx-1), d3idx:tbp3,:,:) &
                        + invrad*tb_batch_prod((invd1idx+1):tbp1,&
                        (invd2idx+1):tbp2,1:invd3idx,:,:)

                   if (d3idx /= 1) then

                      ! Area 4

                      if (abs(pnt(1)-tbp1)+abs(pnt(2)-tbp2)+ &
                           abs(pnt(3)-tbp3) /= 0) then
                         invrad=1.0_DP/sqrt(d1sq*real((pnt(1)-tbp1)**2,kind=DP)&
                              + d2sq*real((pnt(2)-tbp2)**2,kind=DP) &
                              + d3sq*real((pnt(3)-tbp3)**2,kind=DP))
                      else
                         invrad = sincint
                      end if

                      tb_batch(1:(d1idx-1),1:(d2idx-1),1:(d3idx-1),:,:)=&
                           tb_batch(1:(d1idx-1),1:(d2idx-1),1:(d3idx-1),:,:) &
                           + invrad*tb_batch_prod((invd1idx+1):tbp1,&
                           (invd2idx+1):tbp2,(invd3idx+1):tbp3,:,:)

                   end if

                end if

                if (d3idx /= 1) then

                   ! Area 5

                   if (abs(pnt(1)-tbp1)+abs(pnt(2))+abs(pnt(3)-tbp3)&
                        /= 0) then
                      invrad = 1.0_DP/sqrt(d1sq*real((pnt(1)-tbp1)**2,kind=DP)&
                           + d2sq*real((pnt(2))**2,kind=DP) &
                           + d3sq*real((pnt(3)-tbp3)**2,kind=DP))
                   else
                      invrad = sincint
                   end if

                   tb_batch(1:(d1idx-1),d2idx:tbp2,1:(d3idx-1),:,:) =&
                        tb_batch(1:(d1idx-1),d2idx:tbp2,1:(d3idx-1),:,:)&
                        + invrad*tb_batch_prod((invd1idx+1):tbp1,&
                        1:invd2idx,(invd3idx+1):tbp3,:,:)                        ! (2) 11%

                end if
             endif

             if (d2idx /= 1) then

                ! Area 6
                if (abs(pnt(1))+abs(pnt(2)-tbp2)+abs(pnt(3)) /= 0) then
                   invrad = 1.0_DP/sqrt(d1sq*real((pnt(1))**2,kind=DP) &
                        + d2sq*real((pnt(2)-tbp2)**2,kind=DP) &
                        + d3sq*real((pnt(3))**2,kind=DP))
                else
                   invrad = sincint
                end if

                tb_batch(d1idx:tbp1, 1:(d2idx-1),d3idx:tbp3,:,:) =&
                     tb_batch(d1idx:tbp1,1:(d2idx-1),d3idx:tbp3,:,:) &
                     + invrad*tb_batch_prod(1:invd1idx,&
                     (invd2idx+1):tbp2,1:invd3idx,:,:)                          ! (4) 8.7%
                if (d3idx /=1) then

                   !Area 7
                   if (abs(pnt(1))+abs(pnt(2)-tbp2)+&
                        abs(pnt(3)-tbp3) /= 0) then
                      invrad = 1.0_DP/sqrt(d1sq*real((pnt(1))**2,kind=DP) &
                           + d2sq*real((pnt(2)-tbp2)**2,kind=DP) &
                           + d3sq*real((pnt(3)-tbp3)**2,kind=DP))
                   else
                      invrad = sincint
                   end if

                   tb_batch(d1idx:tbp1, 1:(d2idx-1),1:(d3idx-1),:,:)=&
                        tb_batch(d1idx:tbp1,1:(d2idx-1), 1:(d3idx-1),:,:) &
                        + invrad*tb_batch_prod(1:invd1idx,&
                        (invd2idx+1):tbp2,(invd3idx+1):tbp3,:,:)

                end if
             end if

             if (d3idx /= 1) then

                !Area 8
                if (abs(pnt(1))+abs(pnt(2))+ abs(pnt(3)-tbp3) /= 0) then
                   invrad = 1.0_DP/sqrt(d1sq*real((pnt(1))**2,kind=DP) &
                        + d2sq*real((pnt(2))**2,kind=DP) &
                        + d3sq*real((pnt(3)-tbp3)**2,kind=DP))
                else
                   invrad = sincint
                end if

                tb_batch(d1idx:tbp1, d2idx:tbp2, 1:(d3idx-1),:,:)=&
                     tb_batch(d1idx:tbp1,d2idx:tbp2, 1:(d3idx-1),:,:) &
                     + invrad*tb_batch_prod(1:invd1idx,&
                     1:invd2idx,(invd3idx+1):tbp3,:,:)                                ! (1) 12%

             end if

          end do d3
       end do d2
    end do d1

    tb_batch = tb_batch * pub_cell%weight

    call utils_trace_out('hf_exchange_npa_batch_var1')

  end subroutine hf_exchange_npa_batch_var1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_npa_batch_var2(tb_batch, tb_batch_prod, centrea,&
       centreb, batch_size)
    !========================================================================!
    ! Variant 2 of hf_exchange_npa_batch.                                  !
    ! Calculates the potential U_{\beta,\delta} in a tightbox (tb) on atom A,!
    ! by integrating P(r2)/|r1-r2| over tb 2, for every point in tb 1, with  !
    ! batches over \delta.                                                   !
    ! This variant is an optimized version of variant 1, and does NOT intro- !
    ! duce any further approximations. The optimization consists in three    !
    ! things: 1) only the points with non-zero P(r2) in tb2 are taken into   !
    ! account, these are copied over to a packed array to improve cache      !
    ! efficiency. 2) only the points within the localization sphere on A     !
    ! (in tb1) are taken into account, the potential is left at zero for     !
    ! points in tb 1 outside the localization sphere. 3) a lookup table is   !
    ! employed for the inverse distance between points in tb 1 and tb 2,     !
    ! exploiting the fact that 1/r depends only on the difference between    !
    ! the point indices on the grid.                                         !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    ! Refer to hf_exchange_npa_batch_var1, the original routine.           !
    !------------------------------------------------------------------------!
    ! Brief explanation of approach:                                         !
    ! 1) By examining the corners of tb1 and tb2, find the range of values   !
    !    which the distance components (along a1,a2,a3), measured in         !
    !    gridpoints, between points in tb1 and tb2 may assume.               !
    ! 2) Populate a lookup array of 1/r, indexed by the distance components. !
    ! 3) Examine tb2 and copy all non-zero points to a packed array. To a    !
    !    separate packed array, d_packed, copy the offset of these points    !
    !    wrt the origin of tb1 (not tb2).                                    !
    ! 4) Go over all points in tb1, immediately excluding the ones outside   !
    !    the localization sphere of A. For the remaining points, add the     !
    !    contribution from the points in packed version of tb2, looking up   !
    !    1/r in the lookup table on the fly.                                 !
    !------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2010, using the original routine !
    ! by Quintin Hill as a template.                                         !
    !========================================================================!

    use constants, only: DP, PI
    use geometry, only: magnitude, point, operator(+), operator(*)
    use simulation_cell, only: pub_cell, pub_maxtight_pts1, pub_maxtight_pts2, &
         pub_maxtight_pts3
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! qoh: Arguments
    integer,       intent(in)  :: batch_size
    real(kind=DP), intent(out) :: tb_batch(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3,batch_size, pub_cell%num_spins)
    real(kind=DP), intent(in) :: tb_batch_prod(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3,batch_size, pub_cell%num_spins)
    type(atom_centre), intent(in) :: centrea
    type(atom_centre), intent(in) :: centreb

    ! jd: Local variables
    real(kind=DP) :: d1sq, d2sq, d3sq ! grid spacing squared
    integer :: tbp1, tbp2, tbp3 ! number of tightbox grid points in direction
    integer :: startpnt1(3), startpnt2(3) ! bottom-lhs corners of tb1 and tb2
    real(kind=DP) :: sincint ! sinc integral

    integer :: tbdisplacement(3) ! displacement vector between the tb's
    integer :: tb1corner(3), tb2corner(3) ! temporaries, corners of tb's
    integer :: d1idx, d2idx, d3idx, e1idx, e2idx, e3idx ! indices in tb's
    integer :: idx, nidx ! indices in packed array for tb2
    integer, allocatable :: d_packed(:,:) ! distance components in packed array
    real(kind=DP), allocatable :: d1d2d32invr(:,:,:) ! lookup array for 1/r
    real(kind=DP), allocatable :: tb_batch_prod_packed(:,:,:) ! packed tb2
    real(kind=DP) :: invr ! temporary, 1/r
    integer :: ierr ! Error flag
    integer, parameter :: anyinteger = 1 ! any integer, for huge()
    ! jd: minimum, maximum and current (temporary) distance components between
    !     points in tb1 and tb2
    integer :: mind1, mind2, mind3
    integer :: maxd1, maxd2, maxd3
    integer :: curd1, curd2, curd3
    real(kind=DP) :: tb1_rad ! NGWF radius of atom A (the one in tb1)
    real(kind=DP) :: r ! temporary, distance
    type(POINT) :: minus_rvec ! temporary

    integer :: n_points_tb1 ! counter for verbose output

    ! -----------------------------------------------------------------------

    call utils_trace_in('hf_exchange_npa_batch_var2')

    ! qoh: Get useful constants
    tbp1 = pub_maxtight_pts1
    tbp2 = pub_maxtight_pts2
    tbp3 = pub_maxtight_pts3
    d1sq = pub_cell%d1**2
    d2sq = pub_cell%d2**2
    d3sq = pub_cell%d3**2

    ! qoh: Potential at centre of a psinc function due to the psinc function.
    ! qoh: This is calculated by treating the psinc function as an l=0 spherical
    ! qoh: wave with a cutoff radius equal to the grid spacing. (The spherical
    ! qoh: harmonic factor is removed.)
    sincint = 8.0_DP * min(pub_cell%d1,pub_cell%d2,pub_cell%d3)**2 / &
         (PI*pub_cell%weight)

    startpnt1 = centrea%tbcellstart
    startpnt2 = centreb%tbcellstart
    tb1_rad = centrea%radius

    ! jd: Clear the output array
    tb_batch = 0.0_DP

    ! jd: Find range of d1, d2, d3 by looking at every 8x8 corner
    !     combination between tb1 and tb2
    tbdisplacement = startpnt2 - startpnt1
    mind1=huge(anyinteger)
    mind2=huge(anyinteger)
    mind3=huge(anyinteger)
    maxd1=-huge(anyinteger)
    maxd2=-huge(anyinteger)
    maxd3=-huge(anyinteger)
    do d1idx=0,1
       tb1corner(1) = startpnt1(1) + d1idx * tbp1
       do d2idx=0,1
          tb1corner(2) = startpnt1(2) + d2idx * tbp2
          do d3idx=0,1
             tb1corner(3) = startpnt1(3) + d3idx * tbp3
             do e1idx=0,1
                tb2corner(1) = startpnt2(1) + e1idx * tbp1
                do e2idx=0,1
                   tb2corner(2) = startpnt2(2) + e2idx * tbp2
                   do e3idx=0,1
                      tb2corner(3) = startpnt2(3) + e3idx * tbp3
                      curd1 = tb2corner(1) - tb1corner(1)
                      curd2 = tb2corner(2) - tb1corner(2)
                      curd3 = tb2corner(3) - tb1corner(3)
                      if(curd1 < mind1) mind1=curd1
                      if(curd2 < mind2) mind2=curd2
                      if(curd3 < mind3) mind3=curd3
                      if(curd1 > maxd1) maxd1=curd1
                      if(curd2 > maxd2) maxd2=curd2
                      if(curd3 > maxd3) maxd3=curd3
                   end do
                end do
             end do
          end do
       end do
    end do

    ! jd: Allocate a lookup array to store 1/r for all possible d1, d2, d3
    allocate(d1d2d32invr(mind1:maxd1,mind2:maxd2,mind3:maxd3),stat=ierr)
    call utils_alloc_check('hf_exchange_npa_batch_var2','d1d2d32invr',ierr)

    ! jd: Populate the lookup array
    !     with all possible distances between points in tb1 and tb2
    do d3idx=mind3,maxd3
       do d2idx=mind2,maxd2
          do d1idx=mind1,maxd1
             if(d1idx == 0 .and. d2idx == 0 .and. d3idx == 0) then
                d1d2d32invr(d1idx,d2idx,d3idx) = sincint ! jd: Avoid singularity
             else
                d1d2d32invr(d1idx,d2idx,d3idx) = 1.0/sqrt( &
                     d1sq*d1idx*d1idx + d2sq*d2idx*d2idx + d3sq*d3idx*d3idx)
             end if
          end do
       end do
    end do

    ! jd: Allocate temporary arrays for packed versions
    !     - of partial distance components for tb2
    allocate(d_packed(3,tbp1*tbp2*tbp3),stat=ierr)
    call utils_alloc_check('hf_exchange_npa_batch_var2','d_packed',ierr)
    !     - of data in tb_batch_prod
    allocate(tb_batch_prod_packed(batch_size, &
         pub_cell%num_spins,tbp1*tbp2*tbp3),stat=ierr)
    call utils_alloc_check('hf_exchange_npa_batch_var2', &
         'tb_batch_prod_packed',ierr)

    ! jd: Pack the non-zero P's from tb2 into the packed array, at idx
    idx=1
    do d3idx=1,tbp3
       do d2idx=1,tbp2
          do d1idx=1,tbp1

             ! jd: For non-zero P's store their location and value
             !     It suffices to check the first-of-batch-first-spin
             !     element, but all must be stored
             if(tb_batch_prod(d1idx,d2idx,d3idx,1,1) /= 0.0_DP) then
                tb_batch_prod_packed(:,:,idx) = &
                     tb_batch_prod(d1idx,d2idx,d3idx,:,:)
                d_packed(1,idx) = (d1idx-1)+tbdisplacement(1)
                d_packed(2,idx) = (d2idx-1)+tbdisplacement(2)
                d_packed(3,idx) = (d3idx-1)+tbdisplacement(3)
                idx = idx + 1
             end if

          end do
       end do
    end do

    ! jd: Remember the number of elements in the packed array
    nidx = idx - 1

    ! jd: Perform the actual computation, using the packed array
    !     and ignoring the elements outside the localization sphere on A
    n_points_tb1 = 0
    do d3idx=1,tbp3
       do d2idx=1,tbp2
          do d1idx=1,tbp1

             ! jd: Find the vector 'r', from the corner of tb1 to current point
             minus_rvec = &
                  real(-(d1idx-1),kind=DP) * pub_cell%d1 * pub_cell%a1_unit + &
                  real(-(d2idx-1),kind=DP) * pub_cell%d2 * pub_cell%a2_unit + &
                  real(-(d3idx-1),kind=DP) * pub_cell%d3 * pub_cell%a3_unit

             ! jd: Subtract 'r' from the vector to the centre, the magnitude of
             !     the result is the distance of current point from tb1 centre
             r = magnitude(centrea%intb + minus_rvec)

             ! jd: Ignore all points outside the localization sphere on A
             if(r > tb1_rad) cycle

             n_points_tb1 = n_points_tb1 + 1

             ! jd: Go over contrib. from all points in the packed array for tb2
             do idx=1,nidx
                invr = d1d2d32invr( &
                     d_packed(1,idx)-d1idx+1, &
                     d_packed(2,idx)-d2idx+1, &
                     d_packed(3,idx)-d3idx+1 )

                tb_batch(d1idx,d2idx,d3idx,:,:) = &
                     tb_batch(d1idx,d2idx,d3idx,:,:) + &
                     invr * tb_batch_prod_packed(:,:,idx)
             end do

          end do
       end do
    end do

    tb_batch = tb_batch * pub_cell%weight

    ! jd: Save some statistics on the optimization for verbose output
    hfx_int_total_calls = hfx_int_total_calls + 1
    hfx_int_total_packed = hfx_int_total_packed + &
         real(nidx,kind=DP)/real(tbp1*tbp2*tbp3,kind=DP)
    hfx_int_total_points_tb1 = hfx_int_total_points_tb1 + &
         real(n_points_tb1,kind=DP)/real(tbp1*tbp2*tbp3,kind=DP)

    ! jd: Clean up the temporary arrays
    deallocate(tb_batch_prod_packed,stat=ierr)
    call utils_dealloc_check('hf_exchange_npa_batch_var2',&
         'tb_batch_prod_packed',ierr)

    deallocate(d_packed,stat=ierr)
    call utils_dealloc_check('hf_exchange_npa_batch_var2','d_packed',ierr)

    deallocate(d1d2d32invr,stat=ierr)
    call utils_dealloc_check('hf_exchange_npa_batch_var2','d1d2d32invr',ierr)

    call utils_trace_out('hf_exchange_npa_batch_var2')

  end subroutine hf_exchange_npa_batch_var2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_npa_batch_var3(tb_batch, tb_batch_prod, centrea,&
       centreb, batch_size)
    !========================================================================!
    ! Variant 3 of hf_exchange_npa_batch.                                  !
    ! Calculates the potential U_{\beta,\delta} in a tightbox (tb) on atom A,!
    ! by integrating P(r2)/|r1-r2| over tb 2, for every point in tb 1, with  !
    ! batches over \delta.                                                   !
    ! This variant is an optimized version of variant 2 and it DOES          !
    ! introduce additional approximations, by peforming coarse-graining on   !
    ! the values of P(r2) that are small in magnitude and far enough.        !
    ! With the default parameters (default_coarse_grain_threshold = 1E-6 and !
    ! padmargin=5) the exchange matrix was found to be accurate to no worse  !
    ! than 0.0000005 (5E-7) in any element and to 2E-9 on average.           !
    ! Less accuracy than this will likely impact convergence, thus increasing!
    ! default_coarse_grain_threshold or lowering the padmargin is not        !
    ! recommended. Decreasing default_coarse_grain_threshold will give better!
    ! accuracy, but fewer points will be coarse-grained, lowering perfor-    !
    ! mance. Increasing padmargin will give better accuracy, but only more   !
    ! distant points will be coarse-grained, lowering performance.           !
    ! The coarse-graining in this variant is employed on top of the optimiza-!
    ! tions described in variant 2. Since coarse-graining is performed only  !
    ! for points that are sufficiently distant, it is pointless to use this  !
    ! variant for small (Natoms < roughly 80) systems. In smaller systems it !
    ! is better to use variant 2, which does not introduce any approximations!
    ! and will likely be faster.                                             !
    !------------------------------------------------------------------------!
    ! Arguments:                                                             !
    ! Refer to hf_exchange_npa_batch_var1, the original routine.           !
    !------------------------------------------------------------------------!
    ! Brief explanation of approach:                                         !
    ! 1) By examining the corners of tb1 and tb2, find the range of values   !
    !    which the distance components (along a1,a2,a3), measured in         !
    !    gridpoints, between points in tb1 and tb2 may assume.               !
    ! 2) Populate a lookup array of 1/r, indexed by the distance components. !
    ! 3) Divide tb1 into a 3x3x3 array of equally sized subtightboxes        !
    !    (last subtb may be slightly larger).                                !
    ! 4) For every subtb of tb1...                                           !
    !    4a) Imagine a padded version, with padmargin extra points in every  !
    !        direction.                                                      !
    !    4b) Determine if the padded subtb overlaps with tb2, taking PBC     !
    !        into account. If yes, no coarse-graining will be performed,     !
    !        skip to 4x).                                                    !
    !    4c) Examine tb2 and copy all 'large' points to a packed array. To a !
    !        separate packed array, d_packed, copy the offset of these points!
    !        wrt the origin of tb1 (not tb2). Use different packed arrays for!
    !        each subtb of tb1. 'Large' points are those whose P(r2) is above!
    !        default_coarse_grain_threshold.                                 !
    !    4d) Examine the remaining points in tb2, by looking at 3x3x3-point  !
    !        blocks. Coarse grain the contents of each 27-point block by     !
    !        putting the 27 points into the packed array as in 4c), but as   !
    !        one point. Ignore blocks of all zeroes. Use different packed    !
    !        arrays for each subtb of tb1.                                   !
    !    4e) If tb2 did not divide without remainder into 3x3x3 blocks, take !
    !        care of the remaining margins, _not_ performing coarse-graining !
    !        on them, copying them to the packed array as in 4c). Use        !
    !        different packed arrays for each subtb of tb1.                  !
    !    4f) Loop to 4a until all subtbs are exhausted, then go to 5.        !
    !    4x) Examine tb2 and copy all non-zero points to a packed array. To a!
    !        separate packed array, d_packed, copy the offset of these points!
    !        wrt the origin of tb1 (not tb2). Use different packed arrays    !
    !        for each subtb of tb1. Loop to 4a until all subtbs are          !
    !        exhausted, then go to 5.                                        !
    ! 5) Go over all subtbs of tb1 and over all points in every subtb1,      !
    !    immediately excluding the points outside the localization sphere    !
    !    of A. For the remaining points, add the contribution from the       !
    !    points in appropriate packed version of tb2, looking up 1/r in the  !
    !    lookup table on the fly.                                            !
    !------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in February 2010, using the original routine !
    ! by Quintin Hill as a template.                                         !
    !========================================================================!

    use comms, only: comms_abort
    use constants, only: DP, PI, stdout
    use geometry, only: magnitude, point, operator(+), operator(*)
    use simulation_cell, only: pub_cell, pub_maxtight_pts1, pub_maxtight_pts2, &
         pub_maxtight_pts3
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! qoh: Arguments
    integer,       intent(in)  :: batch_size
    real(kind=DP), intent(out) :: tb_batch(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3,batch_size, pub_cell%num_spins)
    real(kind=DP), intent(in) :: tb_batch_prod(pub_maxtight_pts1,&
         pub_maxtight_pts2,pub_maxtight_pts3,batch_size, pub_cell%num_spins)
    type(atom_centre), intent(in) :: centrea
    type(atom_centre), intent(in) :: centreb

    ! jd: Local variables
    real(kind=DP) :: d1sq, d2sq, d3sq ! grid spacing squared
    integer :: tbp1, tbp2, tbp3 ! number of tightbox grid points in direction
    integer :: startpnt1(3), startpnt2(3) ! bottom-lhs corners of tb1 and tb2
    real(kind=DP) :: sincint ! sinc integral

    integer :: tbdisplacement(3) ! displacement vector between the tb's
    integer :: tb1corner(3), tb2corner(3) ! temporaries, corners of tb's
    integer :: d1idx, d2idx, d3idx, e1idx, e2idx, e3idx ! indices in tb's
    integer :: d1j, d2j, d3j, d1ij, d2ij, d3ij ! indices in tb's

    integer :: idx ! index in packed arrays for tb2
    integer :: nidx(0:2,0:2,0:2) ! lenghts of packed arrays for tb2
    ! jd: distance components in packed arrays
    integer, allocatable :: d_packed(:,:,:,:,:)
    real(kind=DP), allocatable :: d1d2d32invr(:,:,:) ! lookup array for 1/r
    real(kind=DP), allocatable :: tb_batch_prod_packed(:,:,:,:,:,:) ! packed tb2
    ! jd: Working copy of tb_batch_prod
    real(kind=DP), allocatable :: local_tb_batch_prod(:,:,:,:,:)
    real(kind=DP), allocatable :: accum(:,:) ! temporary
    real(kind=DP) :: invr ! temporary, 1/r
    integer :: ierr ! error flag
    integer, parameter :: anyinteger = 1 ! any integer, for huge()
    ! jd: minimum, maximum and current (temporary) distance components between
    !     points in tb1 and tb2
    integer :: mind1, mind2, mind3
    integer :: maxd1, maxd2, maxd3
    integer :: curd1, curd2, curd3
    real(kind=DP) :: tb1_rad ! NGWF radius of atom A (the one in tb1)
    real(kind=DP) :: r ! temporary, distance
    type(POINT) :: minus_rvec ! temporary
    ! jd: Indices for dividing tb1 into subtightboxes
    integer :: d1subidx, d2subidx, d3subidx, e1subidx, e2subidx, e3subidx
    ! jd: Indices and numbers of points in subtightboxes
    integer :: subtbp1, subtbp2, subtbp3
    integer :: subtbp1_last, subtbp2_last, subtbp3_last
    integer :: subtbp1_now, subtbp2_now, subtbp3_now
    integer :: subbox1corner1(3), subbox1corner2(3) ! corners of subtb's and
    integer :: subbox1paddedcorner1(3), subbox1paddedcorner2(3) ! padded subtb's

    ! jd: Counters for coarse-graining
    integer :: n_nonzero, n_coarse_from, n_coarse_to, n_fine_large, &
         n_fine_margin, n_zero, n_coarse_grain_instances, &
         n_no_coarse_grain_instances
    logical :: coarse_grain ! True, if coarse-graining tb2 in this subtb1
    real(kind=DP) :: coarse_grain_threshold ! P's below this will be CG'ed
    integer :: n_points_tb1 ! counter for verbose output

    ! jd: Parameters
    real(kind=DP) :: default_coarse_grain_threshold = 1D-6
    integer, parameter :: padmargin = 5 ! Minimum distance for coarse-graining

    ! -----------------------------------------------------------------------

    call utils_trace_in('hf_exchange_npa_batch_var3')

    ! qoh: Get useful constants
    tbp1 = pub_maxtight_pts1
    tbp2 = pub_maxtight_pts2
    tbp3 = pub_maxtight_pts3
    d1sq = pub_cell%d1**2
    d2sq = pub_cell%d2**2
    d3sq = pub_cell%d3**2

    ! qoh: Potential at centre of a psinc function due to the psinc function.
    ! qoh: This is calculated by treating the psinc function as an l=0 spherical
    ! qoh: wave with a cutoff radius equal to the grid spacing. (The spherical
    ! qoh: harmonic factor is removed.)
    sincint = 8.0_DP * min(pub_cell%d1,pub_cell%d2,pub_cell%d3)**2 / &
         (PI*pub_cell%weight)

    startpnt1 = centrea%tbcellstart
    startpnt2 = centreb%tbcellstart
    tb1_rad = centrea%radius

    ! jd: Clear the output array, counters
    tb_batch = 0.0_DP
    n_coarse_grain_instances = 0
    n_no_coarse_grain_instances = 0

    ! jd: Make a working copy of the array that contains P(r2) in tb2
    allocate(local_tb_batch_prod(tbp1,tbp2,tbp3,batch_size,pub_cell%num_spins),&
         stat=ierr)
    call utils_alloc_check('hf_exchange_npa_batch_var3',&
         'local_tb_batch_prod',ierr)
    local_tb_batch_prod = tb_batch_prod

    ! jd: Find range of d1, d2, d3 by looking at every 8x8 corner
    !     combination between tb1 and tb2
    tbdisplacement = startpnt2 - startpnt1
    mind1=huge(anyinteger)
    mind2=huge(anyinteger)
    mind3=huge(anyinteger)
    maxd1=-huge(anyinteger)
    maxd2=-huge(anyinteger)
    maxd3=-huge(anyinteger)
    do d1idx=0,1
       tb1corner(1) = startpnt1(1) + d1idx * tbp1
       do d2idx=0,1
          tb1corner(2) = startpnt1(2) + d2idx * tbp2
          do d3idx=0,1
             tb1corner(3) = startpnt1(3) + d3idx * tbp3
             do e1idx=0,1
                tb2corner(1) = startpnt2(1) + e1idx * tbp1
                do e2idx=0,1
                   tb2corner(2) = startpnt2(2) + e2idx * tbp2
                   do e3idx=0,1
                      tb2corner(3) = startpnt2(3) + e3idx * tbp3
                      curd1 = tb2corner(1) - tb1corner(1)
                      curd2 = tb2corner(2) - tb1corner(2)
                      curd3 = tb2corner(3) - tb1corner(3)
                      if(curd1 < mind1) mind1=curd1
                      if(curd2 < mind2) mind2=curd2
                      if(curd3 < mind3) mind3=curd3
                      if(curd1 > maxd1) maxd1=curd1
                      if(curd2 > maxd2) maxd2=curd2
                      if(curd3 > maxd3) maxd3=curd3
                   end do
                end do
             end do
          end do
       end do
    end do

    ! jd: Allocate a lookup array to store 1/r for all possible d1, d2, d3
    allocate(d1d2d32invr(mind1:maxd1,mind2:maxd2,mind3:maxd3),stat=ierr)
    call utils_alloc_check('hf_exchange_npa_batch_var3','d1d2d32invr',ierr)

    ! jd: Populate the lookup array
    !     with all possible distances between points in tb1 and tb2
    do d3idx=mind3,maxd3
       do d2idx=mind2,maxd2
          do d1idx=mind1,maxd1
             if(d1idx == 0 .and. d2idx == 0 .and. d3idx == 0) then
                d1d2d32invr(d1idx,d2idx,d3idx) = sincint ! jd: Avoid singularity
             else
                d1d2d32invr(d1idx,d2idx,d3idx) = 1.0/sqrt( &
                     d1sq*d1idx*d1idx + d2sq*d2idx*d2idx + d3sq*d3idx*d3idx)
             end if
          end do
       end do
    end do

    ! jd: Allocate temporary arrays for packed versions
    !     - of partial distance components for tb2, as seen
    allocate(d_packed(3,tbp1*tbp2*tbp3,0:2,0:2,0:2),stat=ierr)
    call utils_alloc_check('hf_exchange_npa_batch_var3','d_packed',ierr)
    !     - of data in local_tb_batch_prod
    allocate(tb_batch_prod_packed(batch_size, &
         pub_cell%num_spins,tbp1*tbp2*tbp3,0:2,0:2,0:2),stat=ierr)
    call utils_alloc_check('hf_exchange_npa_batch_var3', &
         'tb_batch_prod_packed',ierr)
    ! jd: Allocate a temporary array for coarse-graining P(r2)
    allocate(accum(batch_size,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('hf_exchange_npa_batch.','accum',ierr)

    ! jd: Tightbox 1 will now be split into 3x3x3 blocks (subtb's).
    !     Calculate the lengths of the blocks, last one may be larger.
    subtbp1=tbp1/3
    subtbp2=tbp2/3
    subtbp3=tbp3/3
    subtbp1_last=tbp1-2*subtbp1
    subtbp2_last=tbp2-2*subtbp2
    subtbp3_last=tbp3-2*subtbp3

    ! ---------------------------------------------
    ! jd: Loop over the 27 subtightboxes of tb 1
    ! ---------------------------------------------
    do d3subidx=0,2 ! over d3-components

       ! jd: Find the d3-component of the corners of the current subtb.
       subbox1corner1(3) = startpnt1(3) + d3subidx*subtbp3
       if(d3subidx .lt. 2) then
          subbox1corner2(3) = startpnt1(3) + (d3subidx+1)*subtbp3 - 1
       else
          subbox1corner2(3) = tb1corner(3)
       end if

       do d2subidx=0,2 ! over d2-components

          ! jd: Find the d2-component of the corners of the current subtb.
          subbox1corner1(2) = startpnt1(2) + d2subidx*subtbp2
          if(d2subidx .lt. 2) then
             subbox1corner2(2) = startpnt1(2) + (d2subidx+1)*subtbp2 - 1
          else
             subbox1corner2(2) = tb1corner(2)
          end if

          do d1subidx=0,2 ! over d1-component

             ! jd: Find the d1-component of the corners of the current subtb.
             subbox1corner1(1) = startpnt1(1) + d1subidx*subtbp1
             if(d1subidx .lt. 2) then
                subbox1corner2(1) = startpnt1(1) + (d1subidx+1)*subtbp1 - 1
             else
                subbox1corner2(1) = tb1corner(1)
             end if

             ! jd: Work on a local copy of P(r2)
             local_tb_batch_prod = tb_batch_prod

             ! jd: Determine the corners of a subtightbox, padded by a margin
             subbox1paddedcorner1(:) = subbox1corner1(:) - padmargin
             subbox1paddedcorner2(:) = subbox1corner2(:) + padmargin

             ! jd: Perform coarse-graining only sufficiently far from the
             !     subtightbox, that is, only when the padded subtightbox
             !     and tb2 do not overlap.
             if(.not. internal_do_boxes_overlap(&
                  subbox1paddedcorner1,subbox1paddedcorner2,&
                  startpnt2,tb2corner)) then
                coarse_grain = .true.
                coarse_grain_threshold = default_coarse_grain_threshold
                n_coarse_grain_instances = n_coarse_grain_instances + 1

             else
                coarse_grain = .false.
                coarse_grain_threshold = tiny(1.0_DP)
                n_no_coarse_grain_instances = n_no_coarse_grain_instances + 1
             end if

             ! jd: Pack the P's from tb2 into the packed array, at idx
             !     ... but only the _large_ P's, that is, those above the
             !         coarse-graining threshold
             n_fine_large = 0
             n_zero = 0

             idx=1
             do d3idx=1,tbp3
                do d2idx=1,tbp2
                   do d1idx=1,tbp1

                      ! jd: For large P's store their location and value.
                      !     It suffices to check the first-of-batch-first-spin
                      !     element, but all must be stored.
                      !     Then, remove them from (local copy of) tb_batch_prod
                      if(abs(local_tb_batch_prod(d1idx,d2idx,d3idx,1,1)) &
                           >= coarse_grain_threshold) then
                         tb_batch_prod_packed(:,:,idx,d1subidx,d2subidx,&
                              d3subidx) = local_tb_batch_prod(d1idx,d2idx,&
                              d3idx,:,:)
                         local_tb_batch_prod(d1idx,d2idx,d3idx,:,:) = 0.0_DP
                         d_packed(1,idx,d1subidx,d2subidx,d3subidx) = &
                              (d1idx-1)+tbdisplacement(1)
                         d_packed(2,idx,d1subidx,d2subidx,d3subidx) = &
                              (d2idx-1)+tbdisplacement(2)
                         d_packed(3,idx,d1subidx,d2subidx,d3subidx) = &
                              (d3idx-1)+tbdisplacement(3)
                         idx = idx + 1
                      else
                         if(local_tb_batch_prod(d1idx,d2idx,d3idx,1,1) &
                              == 0.0_DP) then
                            n_zero = n_zero + 1 ! Store for verbose output
                         end if
                      end if
                   end do
                end do
             end do
             n_fine_large = idx-1 ! Useful to store, for verbose output

             ! jd: If not coarse-graining, than we're done
             !     If coarse-graining, coarse-grain the _small_ P(r2)

             n_fine_margin = 0
             n_coarse_from = 0
             n_coarse_to = 0

             if(coarse_grain) then

                ! jd: Go over 3x3x3 blocks in tb2. For now, ignore any margin
                !     that is left if tbp{1,2,3} modulo 3 is not zero.
                do d3idx=2,tbp3-1,3
                   do d2idx=2,tbp2-1,3
                      do d1idx=2,tbp1-1,3

                         n_nonzero=0

                         ! jd: Coarse grain P's in this 3x3x3 block into one
                         !     value, accumulating them in accum() and removing
                         !     them from (local copy of) tb_batch_prod.
                         accum = 0.0_DP
                         do d3j=-1,1
                            d3ij = d3idx+d3j
                            do d2j=-1,1
                               d2ij = d2idx+d2j
                               do d1j=-1,1
                                  d1ij = d1idx+d1j

                                  if(local_tb_batch_prod(d1ij,d2ij,d3ij,1,1) &
                                       /= 0.0_DP) then
                                     n_nonzero = n_nonzero + 1
                                     accum(:,:) = accum(:,:) + &
                                          local_tb_batch_prod(d1ij,d2ij,d3ij,&
                                          :,:)
                                     local_tb_batch_prod(d1ij,d2ij,d3ij,:,:) = &
                                          0.0_DP
                                  end if

                               end do
                            end do
                         end do

                         ! jd: Store the coarse-grained point into the packed
                         !     array, except when it's made up of zeroes
                         if(n_nonzero /= 0) then
                            tb_batch_prod_packed(:,:,idx,d1subidx,d2subidx,&
                                 d3subidx) = accum(:,:)
                            d_packed(1,idx,d1subidx,d2subidx,d3subidx) = &
                                 (d1idx-1)+tbdisplacement(1)
                            d_packed(2,idx,d1subidx,d2subidx,d3subidx) = &
                                 (d2idx-1)+tbdisplacement(2)
                            d_packed(3,idx,d1subidx,d2subidx,d3subidx) = &
                                 (d3idx-1)+tbdisplacement(3)
                            n_coarse_from = n_coarse_from + n_nonzero
                            idx = idx + 1
                         end if
                      end do
                   end do
                end do
                n_coarse_to = idx-1 - n_fine_large

                ! jd: We have coarse-grained the bulk of the points in tb2,
                !     but non-zero, small P(r2) values are still left in
                !     the margins. Go over tb2, pick these values and add
                !     them to the packed array, without coarse-graining.

                do d3idx=1,tbp3
                   do d2idx=1,tbp2
                      do d1idx=1,tbp1

                         if(local_tb_batch_prod(d1idx,d2idx,d3idx,1,1) &
                              /= 0.0_DP) then
                            ! jd: Is this non-zero P on the margin as expected?
                            if(d3idx>=tbp3-1 .or. d2idx>=tbp2-1 &
                                 .or. d1idx>=tbp1-1) then
                               tb_batch_prod_packed(:,:,idx,d1subidx,d2subidx,&
                                    d3subidx) = local_tb_batch_prod(d1idx,&
                                    d2idx,d3idx,:,:)
                               local_tb_batch_prod(d1idx,d2idx,d3idx,:,:) = &
                                    0.0_DP
                               d_packed(1,idx,d1subidx,d2subidx,d3subidx) = &
                                    (d1idx-1)+tbdisplacement(1)
                               d_packed(2,idx,d1subidx,d2subidx,d3subidx) = &
                                    (d2idx-1)+tbdisplacement(2)
                               d_packed(3,idx,d1subidx,d2subidx,d3subidx) = &
                                    (d3idx-1)+tbdisplacement(3)
                               idx = idx + 1
                               ! jd: If not, something is very wrong
                            else
                               write (stdout,'(a)') 'Internal error in &
                                    &hf_exchange_npa_batch_var3'
                               call comms_abort
                            end if
                         end if

                      end do
                   end do
                end do
                n_fine_margin = idx-1 - n_coarse_to - n_fine_large

             end if ! if(coarse-grain)

             ! jd: Remember the number of elements in the packed array
             !     for this subtightbox of tb1
             nidx(d1subidx,d2subidx,d3subidx) = idx - 1

          end do
       end do
    end do
    ! --------------------------------------------------
    ! jd: End of loop over the 27 subtightboxes of tb 1
    ! --------------------------------------------------

    ! -------------------------------------------------------------------------

    ! jd: Perform the actual computation, using the packed array
    !     and ignoring the elements outside the localization sphere on A.
    !     Again, must loop over 27 subtightboxes of tb 1

    n_points_tb1 = 0
    do d3subidx=0,2 ! jd: go over subtightboxes along d3
       ! jd: How long is this subtightbox?
       if(d3subidx<2) then
          subtbp3_now = subtbp3
       else
          subtbp3_now = subtbp3_last
       end if

       do d2subidx=0,2 ! jd: go over subtightboxes along d2
          ! jd: How long is this subtightbox?
          if(d2subidx<2) then
             subtbp2_now = subtbp2
          else
             subtbp2_now = subtbp2_last
          end if

          do d1subidx=0,2 ! jd: go over subtightboxes along d1
             ! jd: How long is this subtightbox?
             if(d1subidx<2) then
                subtbp1_now = subtbp1
             else
                subtbp1_now = subtbp1_last
             end if

             ! jd: Go over the points in this subtightbox
             do e3subidx=1,subtbp3_now
                d3idx = d3subidx*subtbp3 + e3subidx

                do e2subidx=1,subtbp2_now
                   d2idx = d2subidx*subtbp2 + e2subidx

                   do e1subidx=1,subtbp1_now
                      d1idx = d1subidx*subtbp1 + e1subidx

                      ! jd: Find the vector 'r', from the corner of tb1
                      !     to the current point
                      minus_rvec = &
                           real(-(d1idx-1),kind=DP) * pub_cell%d1 * &
                           pub_cell%a1_unit + &
                           real(-(d2idx-1),kind=DP) * pub_cell%d2 * &
                           pub_cell%a2_unit + &
                           real(-(d3idx-1),kind=DP) * pub_cell%d3 * &
                           pub_cell%a3_unit

                      ! jd: Subtract 'r' from the vector to the centre,
                      !     the magnitude of the result is the distance
                      !     of the current point from tb1 centre
                      r=magnitude(centrea%intb + minus_rvec)

                      ! jd: Ignore all points outside the localization
                      !     sphere on A
                      if(r > tb1_rad) cycle

                      n_points_tb1 = n_points_tb1 + 1

                      ! jd: Go over contrib. from all points in the packed array
                      !     representing tb2. Pick the packed array suitable
                      !     for this subtightbox of tb1.
                      do idx=1,nidx(d1subidx,d2subidx,d3subidx)
                         invr = d1d2d32invr( &
                              d_packed(1,idx,d1subidx,d2subidx,d3subidx)-&
                              d1idx+1, &
                              d_packed(2,idx,d1subidx,d2subidx,d3subidx)-&
                              d2idx+1, &
                              d_packed(3,idx,d1subidx,d2subidx,d3subidx)-&
                              d3idx+1 )
                         tb_batch(d1idx,d2idx,d3idx,:,:) = &
                              tb_batch(d1idx,d2idx,d3idx,:,:) + &
                              invr * tb_batch_prod_packed(:,:,idx,&
                              d1subidx,d2subidx,d3subidx)
                      end do

                   end do  ! Loops over the points
                end do     ! in the current
             end do        ! subtightbox of tb1

          end do  ! Loops over the
       end do     ! subtightboxes
    end do        ! of tb1

    tb_batch = tb_batch * pub_cell%weight

    ! jd: Save some statistics on the optimization for verbose output
    hfx_int_total_calls = hfx_int_total_calls + 1
    hfx_int_total_packed = hfx_int_total_packed + &
         real(sum(nidx)/max(1,size(nidx)),kind=DP)/real(tbp1*tbp2*tbp3,kind=DP)
    hfx_int_total_points_tb1 = hfx_int_total_points_tb1 + &
         real(n_points_tb1,kind=DP)/real(tbp1*tbp2*tbp3,kind=DP)
    hfx_int_total_large = hfx_int_total_large + &
         real(n_fine_large,kind=DP)/real(tbp1*tbp2*tbp3,kind=DP)
    hfx_int_total_margin = hfx_int_total_margin + &
         real(n_fine_margin,kind=DP)/real(tbp1*tbp2*tbp3,kind=DP)
    hfx_int_total_coarse_from = hfx_int_total_coarse_from + &
         real(n_coarse_from,kind=DP)/real(tbp1*tbp2*tbp3,kind=DP)
    hfx_int_total_coarse_to = hfx_int_total_coarse_to + &
         real(n_coarse_to,kind=DP)/real(tbp1*tbp2*tbp3,kind=DP)
    hfx_int_total_zero = hfx_int_total_zero + &
         real(n_zero,kind=DP)/real(tbp1*tbp2*tbp3,kind=DP)
    hfx_int_total_cg_instances = hfx_int_total_cg_instances + &
         n_coarse_grain_instances
    hfx_int_total_no_cg_instances = hfx_int_total_no_cg_instances + &
         n_no_coarse_grain_instances

    ! jd: Clean up the temporary arrays
    deallocate(accum,stat=ierr)
    call utils_dealloc_check('hf_exchange_npa_batch_var3','accum',ierr)

    deallocate(tb_batch_prod_packed,stat=ierr)
    call utils_dealloc_check('hf_exchange_npa_batch_var3',&
         'tb_batch_prod_packed',ierr)

    deallocate(d_packed,stat=ierr)
    call utils_dealloc_check('hf_exchange_npa_batch_var3','d_packed',ierr)

    deallocate(d1d2d32invr,stat=ierr)
    call utils_dealloc_check('hf_exchange_npa_batch_var3','d1d2d32invr',ierr)

    deallocate(local_tb_batch_prod,stat=ierr)
    call utils_dealloc_check('hf_exchange_npa_batch_var3', &
         'local_tb_batch_prod',ierr)

    call utils_trace_out('hf_exchange_npa_batch_var3')

  contains
    logical function internal_do_boxes_overlap(box1corner1, box1corner2, &
         box2corner1, box2corner2)
      !========================================================================!
      ! Returns true if two boxes (cuboids), specified by corners, overlap,    !
      ! and false otherwise. Corners are specified as triplets of integers,    !
      ! which are assumed to be indices on the coarse grid [*], this lets PBCs !
      ! be taken into account. This routine works as well in non-cuboid simu-  !
      ! lation cells, where 'box' should be interpreted as 'parallelepiped,    !
      ! not 'cuboid'.                                                          !
      ! [*] Boxes extending _outside_ the coarse grid are handled correctly,   !
      !     this may happen when padded boxes are inspected. This is when PBCs !
      !     come into play -- e.g. a box that starts at (1,1,1) and is padded  !
      !     with a margin of 5, actually wraps around to the other side of the !
      !     cell. Such cases are supported, unless the box wraps around twice  !
      !     or more, which produces an error.                                  !
      !------------------------------------------------------------------------!
      ! Arguments (all input):                                                 !
      !   box1corner1: lhs-bottom corner of box 1                              !
      !   box1corner2: rhs-top corner of box 1                                 !
      !   box2corner1: lhs-bottom corner of box 2                              !
      !   box2corner2: rhs-top corner of box 2                                 !
      !------------------------------------------------------------------------!
      ! Caveat:                                                                !
      !   The routine WILL NOT work when corners are switched, i.e. it is      !
      !   wrong to supply the rhs-top corner in box{1,2}corner1 and lhs-bottom !
      !   corner in box{1,2}corner2. No check is made, to keep things simple.  !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in February 2010.                            !
      !========================================================================!

      use constants, only: stdout
      use comms, only: comms_abort
      use simulation_cell, only: pub_cell

      implicit none

      ! Arguments
      integer, intent(in) :: box1corner1(3), box1corner2(3)
      integer, intent(in) :: box2corner1(3), box2corner2(3)

      ! Local variables
      integer :: total_pt(3) ! Extent of the coarse grid
      ! jd: Describe where the projections begin and end
      !     b1, b2: box number
      !     (1) refers to the part inside cell, (2) is the periodic image
      integer :: b1min(2), b2min(2), b1max(2), b2max(2)
      logical :: box_has_periodic_image(2)
      logical :: overlap(3) ! Do we have overlap of projections onto axes?
      integer :: axis, n, n1, n2 ! Indices

      ! ---------------------------------------------------------------------

      internal_do_boxes_overlap = .false.

      ! jd: Useful constants
      total_pt(1) = pub_cell%total_pt1
      total_pt(2) = pub_cell%total_pt2
      total_pt(3) = pub_cell%total_pt3

      overlap(:) = .false.
      ! jd: Examine the projections of the box onto all three axes
      do axis=1,3

         ! jd: Project the box
         b1min(1)=box1corner1(axis)
         b1max(1)=box1corner2(axis)
         b2min(1)=box2corner1(axis)
         b2max(1)=box2corner2(axis)

         ! jd: Assume no periodicity issues, until we detect them
         box_has_periodic_image(:) = .false.

         ! jd: Then deal with another unlikely scenario, where the
         !     box is completely outside the cell (both corners should
         !     be periodized)
         !     Example: +--+ |            |   <- this
         !                   |       +--+ |   <- gets mapped to this
         if(b1min(1)<1 .and. b1max(1)<1) then
            b1min(1) = b1min(1) + total_pt(axis)
            b1max(1) = b1max(1) + total_pt(axis)
         end if
         if(b2min(1)<1 .and. b2max(1)<1) then
            b2min(1) = b2min(1) + total_pt(axis)
            b2max(1) = b2max(1) + total_pt(axis)
         end if

         ! jd: Reverse of the previous scenario, where the
         !     box is completely outside the cell (both corners should
         !     be periodized)
         !     Example: |            | +--+  <- this
         !              | +--+       |       <- gets mapped to this
         if(b1min(1)>total_pt(axis) .and. b1max(1)>total_pt(axis)) then
            b1min(1) = b1min(1) - total_pt(axis)
            b1max(1) = b1max(1) - total_pt(axis)
         end if
         if(b2min(1)>total_pt(axis) .and. b2max(1)>total_pt(axis)) then
            b2min(1) = b2min(1) - total_pt(axis)
            b2max(1) = b2max(1) - total_pt(axis)
         end if

         ! jd: Now, apply PBCs to create a projection of the _image_, if any
         !     Example: +--|-------+     |   <- this
         !                 +-------+   +-+   <- gets mapped to this
         !     Legend: || denotes the cell, +----+ denotes the box
         if(b1min(1)<1) then
            b1min(2) = b1min(1) + total_pt(axis)
            b1max(2) = total_pt(axis)
            box_has_periodic_image(1) = .true.
            b1min(1) = 1
         end if
         if(b2min(1)<1) then
            b2min(2) = b2min(1) + total_pt(axis)
            b2max(2) = total_pt(axis)
            box_has_periodic_image(2) = .true.
            b2min(1) = 1
         end if
         if(b1max(1)>total_pt(axis)) then
            b1min(2) = 1
            b1max(2) = b1max(1) - total_pt(axis)
            box_has_periodic_image(1) = .true.
            b1max(1) = total_pt(axis)
         end if
         if(b2max(1)>total_pt(axis)) then
            b2min(2) = 1
            b2max(2) = b2max(1) - total_pt(axis)
            box_has_periodic_image(2) = .true.
            b2max(1) = total_pt(axis)
         end if

         ! jd: Verify if PBCs were successfully applied. If not, the boxes
         !     must have wrapped around more than once, which is not allowed
         do n=1,2
            if( (.not. box_has_periodic_image(1)) .and. n == 2) cycle
            if(  b1min(n)<1 .or. b1max(n)<1 .or. &
                 b1min(n)>total_pt(axis) .or. b1max(n)>total_pt(axis) ) then
               write (stdout,'(a,i4,a,i4,2a)') 'Illegal box in &
                    &internal_do_boxes_overlap:', b1min(n),' - ',b1max(n), &
                    ' This indicates that padmargin is too large or the cell &
                    &is ridiculously small or the coarse grid is ridiculously &
                    &coarse or an internal error.','Program execution stops.'
               call comms_abort
            end if
         end do

         do n=1,2
            if( (.not. box_has_periodic_image(2)) .and. n == 2) cycle
            if(  b2min(n)<1 .or. b2max(n)<1 .or. &
                 b2min(n)>total_pt(axis) .or. b2max(n)>total_pt(axis) ) then
               write (stdout,'(a,i4,a,i4,a,a)') 'Illegal box in &
                    &internal_do_boxes_overlap:', b2min(n),' - ',b2max(n), &
                    ' This indicates that padmargin is too large or the cell &
                    &is ridiculously small or the coarse grid is ridiculously &
                    &coarse or an internal error.','Program execution stops.'
               call comms_abort
            end if
         end do

         ! jd: Finally, see if either the projection (or its image) of box 1
         !     overlaps the projection (or its image) of box 2
         do n1=1,2
            if( (.not. box_has_periodic_image(1)) .and. n1 == 2) cycle
            do n2=1,2
               if( (.not. box_has_periodic_image(2)) .and. n2 == 2) cycle
               if( ((b2min(n2) .ge. b1min(n1)) .and. &
                    (b2min(n2) .le. b1max(n1))) .or. &
                    ((b2max(n2) .ge. b1min(n1)) .and. &
                    (b2max(n2) .le. b1max(n1))) .or. &
                    ((b1min(n1) .ge. b2min(n2)) .and. &
                    (b1min(n1) .le. b2max(n2))) .or. &
                    ((b1max(n1) .ge. b2min(n2)) .and. &
                    (b1max(n1) .le. b2max(n2))) &
                    ) overlap(axis) = .true.
            end do
         end do
      end do

      ! jd: Only if all projections overlap, do the boxes themselves overlap
      internal_do_boxes_overlap = (overlap(1) .and. overlap(2) .and. overlap(3))

    end function internal_do_boxes_overlap


  end subroutine hf_exchange_npa_batch_var3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_distribute_sws(sw_first,sw_count,sw_where,n_sws)
    !==========================================================================!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2012.                               !
    !==========================================================================!

    use comms, only: pub_total_num_nodes, pub_on_root
    use constants, only: DP
    use utils, only: utils_assert

    implicit none

    ! Arguments
    integer, intent(in)  :: n_sws
    integer, intent(out) :: sw_first(0:pub_total_num_nodes-1)
    integer, intent(out) :: sw_count(0:pub_total_num_nodes-1)
    integer, intent(out) :: sw_where(1:max_sw_set_size)

    ! Local variables
    character(len=*), parameter :: myself = 'hf_exchange_distribute_sws'
    integer :: portion, remainder
    integer :: p
    integer :: k
    real(kind=DP) :: delta

    ! -------------------------------------------------------------------------
    sw_where(:) = -1
    portion = n_sws / pub_total_num_nodes
    remainder = n_sws - portion * pub_total_num_nodes

    ! Deal out 'portion' to every node
    sw_count(:) = portion

    ! Deal out the remainder to every delta'th processor, if necessary
    if(remainder /= 0) then
       delta = real(pub_total_num_nodes,kind=DP) / real(remainder,kind=DP)
       do k = 1, remainder
          p = nint(real(k-1,kind=DP) * delta)
          sw_count(p) = sw_count(p) + 1
       end do
    end if

    call utils_assert(sum(sw_count) == n_sws,'Sanity check failed in '//myself)

    ! Construct sw_first, sw_here from sw_count
    do p = 0, pub_total_num_nodes - 1
       if(p == 0) then
          sw_first(p) = 1
       else
          sw_first(p) = sw_first(p-1) + sw_count(p-1)
       end if

       if(sw_count(p) == 0) sw_first(p) = -1

       sw_where(sw_first(p):sw_first(p)+sw_count(p)-1) = p

       if(pub_on_root) then
          if(sw_count(p) > 0) then
             write(*,'(a,i4,a,i4,a,i4,a)') 'Node ',p,': SWs ',sw_first(p), &
                  ' to ',sw_first(p)+sw_count(p)-1,'.'
          else
             write(*,'(a,i4,a)') 'Node ',p,': no SWs.'
          end if
       end if

    end do

  end subroutine hf_exchange_distribute_sws

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_obtain_sw_coeffs(out_sph_coeffs, my_sph_coeffs, sw, &
       sw_first, sw_count, sw_where, who_is_done)
    !==========================================================================!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2012.                               !
    !==========================================================================!

    use chebyshev_rep, only: SPHERE_COEFFS, HANDLE_SET, &
        cheb_recv_coeffs_initiate, cheb_test_completion
    use comms, only: pub_my_node_id, pub_total_num_nodes, comms_send, &
        comms_irecv, comms_recv, comms_probe, comms_wait
    use utils, only: utils_assert, utils_abort

    implicit none

    ! Arguments
    type(SPHERE_COEFFS), intent(inout) :: out_sph_coeffs
    type(SPHERE_COEFFS), intent(in)  :: my_sph_coeffs(:)
    integer, intent(in) :: sw
    integer, intent(in) :: sw_first(0:pub_total_num_nodes-1)
    integer, intent(in) :: sw_count(0:pub_total_num_nodes-1)
    integer, intent(in) :: sw_where(1:max_sw_set_size)
    logical, intent(inout) :: who_is_done(0:pub_total_num_nodes-1)

    ! Local variables
    integer :: sw_loc   ! Node-local index corresponding to sw
    integer :: sw_owner ! Node that owns this sw
    integer :: node
    logical :: we_requested_a_sw
    integer :: send_handle
    logical :: done ! ignored here
    logical :: ready
    type(HANDLE_SET) :: handles

    ! -------------------------------------------------------------------------
    call utils_assert(sw > 0 .and. sw <= max_sw_set_size, &
         'Internal error in hf_exchange_obtain_sw_coeffs')

    if(.false. .and. pub_hfx_debug) then
       write(stdout,*) 'Node ', pub_my_node_id,' enters obtain_sw_coeffs'
       call utils_flush(stdout,.true.)
    end if

    ! jd: If sw is local to this node, then that's very convenient
    if(sw_where(sw) == pub_my_node_id) then
       sw_loc = sw - sw_first(pub_my_node_id) + 1
       out_sph_coeffs = my_sph_coeffs(sw_loc)
       we_requested_a_sw = .false.

       if(.false. .and. pub_hfx_debug) then
          write(stdout,*) 'Node ', pub_my_node_id,' used LOCAL in obtain_sw_coeffs'
          call utils_flush(stdout,.true.)
       end if

    else
       we_requested_a_sw = .true.

       ! jd: If not, find out who owns it and ask them for it
       sw_owner = sw_where(sw)

       if(.false. .and. pub_hfx_debug) then
          write(stdout,*) 'Node ', pub_my_node_id,' will request SW ', sw, ' from ',sw_owner
          call utils_flush(stdout,.true.)
       end if

       call comms_send(sw_owner, sw, 1, tag = SW_REQUEST_TAG, &
            return_handle = send_handle, add_to_stack = .false.)
       call comms_wait(send_handle)

       if(.false. .and. pub_hfx_debug) then
          write(stdout,*) 'Node ', pub_my_node_id,' will successfully requested SW ', sw, ' from ',sw_owner
          call utils_flush(stdout,.true.)
       end if

    end if

    if(.false. .and. pub_hfx_debug) then
       write(stdout,*) 'Node ', pub_my_node_id,' will now serve requests, if any.'
       call utils_flush(stdout,.true.)
    end if

    ! jd: Satisfy the requests of others
    call hf_exchange_serve_sw_coeffs(done, my_sph_coeffs, &
         sw_first, sw_count, sw_where, who_is_done)

    ! jd: Get our requested sw, if any
    if(we_requested_a_sw) then

       if(.false. .and. pub_hfx_debug) then
          write(stdout,*) 'Node ', pub_my_node_id,' will initiate recv of data from ', sw_owner
          call utils_flush(stdout,.true.)
       end if

       call cheb_recv_coeffs_initiate(handles, sw_owner, out_sph_coeffs, SW_DATA_TAG)

       if(.false. .and. pub_hfx_debug) then
          write(stdout,*) 'Node ', pub_my_node_id,' did  initiate recv of data from ', sw_owner
          call utils_flush(stdout,.true.)
       end if

       ready = .false.
       do while(.not. ready)

          if(.false. .and. pub_hfx_debug) then
             write(stdout,*) 'Node ', pub_my_node_id,' waits for completion'
             call utils_flush(stdout,.true.)
          end if

          call cheb_test_completion(ready, handles)

          if(.false. .and. pub_hfx_debug) then
             if(ready) then
                write(stdout,*) 'Node ', pub_my_node_id,' GOT completion'
             else
                write(stdout,*) 'Node ', pub_my_node_id,' DID NOT GET completion'
             end if
             call utils_flush(stdout,.true.)
          end if

          ! jd: Satisfy the requests of others
          call hf_exchange_serve_sw_coeffs(done, my_sph_coeffs, &
               sw_first, sw_count, sw_where, who_is_done)

       end do
   
    end if

    if(.false. .and. pub_hfx_debug) then
       write(stdout,*) 'Node ', pub_my_node_id,' obtained what it wanted: ',sw
       call utils_flush(stdout,.true.)
    end if

  end subroutine hf_exchange_obtain_sw_coeffs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hf_exchange_serve_sw_coeffs(done, my_sph_coeffs, &
       sw_first, sw_count, sw_where, who_is_done)
    !==========================================================================!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2012.                               !
    !==========================================================================!

    use chebyshev_rep, only: SPHERE_COEFFS, HANDLE_SET, &
         cheb_send_coeffs_initiate
    use comms, only: pub_my_node_id, pub_total_num_nodes, comms_send, &
        comms_irecv, comms_recv, comms_probe, comms_wait
    use utils, only: utils_assert, utils_abort

    implicit none

    ! Arguments
    logical, intent(out) :: done
    type(SPHERE_COEFFS), intent(in)  :: my_sph_coeffs(:)
    integer, intent(in) :: sw_first(0:pub_total_num_nodes-1)
    integer, intent(in) :: sw_count(0:pub_total_num_nodes-1)
    integer, intent(in) :: sw_where(1:max_sw_set_size)
    logical, intent(inout) :: who_is_done(0:pub_total_num_nodes-1)

    ! Local variables
    integer :: sw_loc   ! Node-local index corresponding to sw
    integer :: sw_owner ! Node that owns this sw
    logical :: who_wants_sws(0:pub_total_num_nodes-1)
    integer :: who_wants_which_sws(0:pub_total_num_nodes-1)
    integer :: sw_to_send, sw_to_send_loc
    integer :: requested_sw
    integer :: send_handle
    integer :: node
    type(HANDLE_SET) :: handles

    ! -------------------------------------------------------------------------

    if(.false. .and. pub_hfx_debug) then
       write(stdout,*) 'Node ', pub_my_node_id,' enters serve_sw_coeffs.'
       call utils_flush(stdout,.true.)
    end if

    ! jd: Probe for requests from others
    who_wants_sws(:) = .false.
    do node = 0, pub_total_num_nodes - 1
       call comms_probe(who_wants_sws(node), node, SW_REQUEST_TAG)

       if(.false. .and. pub_hfx_debug) then
          write(stdout,*) 'Node ', pub_my_node_id,' probed out a request from ', node
          call utils_flush(stdout,.true.)
       end if

    end do

    if(.false. .and. pub_hfx_debug .and. .not. any(who_wants_sws)) then
       write(stdout,*) 'Node ', pub_my_node_id,' sees NO requests.'
       call utils_flush(stdout,.true.)
    end if

    ! jd: Find out what SWs do the requests refer to, by actually receiving
    !     the requests
    who_wants_which_sws(:) = -1
    do node = 0, pub_total_num_nodes - 1
       if(who_wants_sws(node)) then

          if(.false. .and. pub_hfx_debug .and. .not. any(who_wants_sws)) then
             write(stdout,*) 'Node ', pub_my_node_id,' will blocking-recv the request from ',node
             call utils_flush(stdout,.true.)
          end if

          call comms_recv(node, requested_sw, 1, SW_REQUEST_TAG)

          if(.false. .and. pub_hfx_debug .and. .not. any(who_wants_sws)) then
             write(stdout,*) 'Node ', pub_my_node_id,' DID  blocking-recv the request from ',node
             call utils_flush(stdout,.true.)
          end if

          who_wants_which_sws(node) = requested_sw

          ! jd: Is it a 'done' notification?
          if(requested_sw == -1) then
             who_is_done(node) = .true.
             who_wants_sws(node) = .false.

             if(.false. .and. pub_hfx_debug .and. .not. any(who_wants_sws)) then
                write(stdout,*) 'Node ', pub_my_node_id,' got *DONE* from ',node
                call utils_flush(stdout,.true.)
             end if

          else

             if(.false. .and. pub_hfx_debug .and. .not. any(who_wants_sws)) then
                write(stdout,*) 'Node ', pub_my_node_id,' got a SW request from ',node, ' for ',requested_sw
                call utils_flush(stdout,.true.)
             end if
          
          end if
       end if
    end do

    done = all(who_is_done)

    ! jd: Satisfy the requests from other nodes
    do node = 0, pub_total_num_nodes - 1
       if(who_wants_sws(node)) then
          sw_to_send = who_wants_which_sws(node)
          sw_to_send_loc = sw_to_send - sw_first(pub_my_node_id) + 1
          if(.false. .and. pub_hfx_debug) then
             write(stdout,*) 'Node ', pub_my_node_id,' serving ', node,' for ', sw_to_send
             call utils_flush(stdout,.true.)
          end if
          call cheb_send_coeffs_initiate(handles, node, my_sph_coeffs(sw_to_send_loc), SW_DATA_TAG)
          if(.false. .and. pub_hfx_debug) then
             write(*,*) 'Node ', pub_my_node_id,' initiated send ', node,' for ', sw_to_send
             call utils_flush(stdout,.true.)
          end if

       end if
    end do

    if(.false. .and. pub_hfx_debug) then
       write(*,*) 'Node ', pub_my_node_id,' exiting serve_sw_coeffs.'
       call utils_flush(stdout,.true.)
    end if

  end subroutine hf_exchange_serve_sw_coeffs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine hf_exchange_memory_estimate(intervals, order, a_batchsize, &
       b_batchsize, n_sws)
    !==========================================================================!
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in January 2012.                               !
    !==========================================================================!

    use comms, only: pub_on_root, pub_total_num_nodes

    implicit none

    ! Arguments
    integer, intent(in) :: intervals, order
    integer, intent(in) :: a_batchsize, b_batchsize
    integer, intent(in) :: n_sws

    ! Local variables
    integer       :: n_stripes
    integer       :: n_points
    real(kind=DP) :: size_of_rep

    ! ----------------------------------------------------------------------
    if(pub_on_root) then 
       n_stripes = intervals * order
       n_points = n_stripes**3
       size_of_rep = real(n_points,kind=DP) * 8.0_DP / 1024.0_DP / 1024.0_DP

       write(*,'(a)') ' --- Memory requirements of HFx Chebyshev engine --- '
       write(*,'(a,f7.1,a)') ' Representation of SWs: ', &
            size_of_rep * real(n_sws,kind=DP) / pub_total_num_nodes, &
            ' MiB per core (estimate).'
       write(*,'(a,f7.1,a)') ' A-Batch of SWs:        ', &
            size_of_rep * real(a_batchsize,kind=DP), ' MiB per core.'
       write(*,'(a,f7.1,a)') ' B-Batch of SWpots:     ', &
            size_of_rep * real(b_batchsize,kind=DP), ' MiB per core.'
       write(*,'(a)') ' --------------------------------------------------- '
       write(*,'(a,f7.1,a)') ' TOTAL                 ', &
            size_of_rep * (real(n_sws,kind=DP) / pub_total_num_nodes + &
            real(a_batchsize,kind=DP) + real(b_batchsize,kind=DP)), &
            ' MiB per core (estimate).'
       write(*,'(a)') ' --------------------------------------------------- '

    end if

  end subroutine hf_exchange_memory_estimate


end module hf_exchange
