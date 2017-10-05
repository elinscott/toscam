! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!          Natural Atomic Orbital analysis module                !
!                                                                !
! This module performs NPA analysis using the converged          !
! density kernel from general ONETEP calculation. Called from    !
! subroutine properties_calculate in properties_mod              !
!----------------------------------------------------------------!
! Written by Louis Lee                                           !
! lpl: 23/09/2011 - Set W range in WWSW to be scaled to [0,1]    !
!                   with W_i < w_min = w_min = 1e-4.             !
!                 - Fixed bug in determining # non-pseudised     !
!                   NMBs (num shells referred to n, not nlm)     !
!                 - Added various AO to PNAO schemes (not to be  !
!                   used under normal circumstances              !
! lpl: 02/09/2011 - Implemented NRB                              !
!                 - inittr format f16.12 to avoid *****          !
! lpl: 26/08/2011 - Converted to dense                           !
! lpl: 19/08/2011 - Attempting to add NRB & revised plotting     !
!                 - Memeory efficiency mod                       !
! lpl: 17/06/2011 with kernel_purify                             !
! lpl: 14/07/2011 - Minor revisions to comment syntax            !
!                 - Corrected 'CMO$END' error in GENNBO output   !
!                 - Density matrix charge is now rounded         !
!                   downwards instead of to the nearest integer  !
!================================================================!

module npa

  implicit none

  private

  public :: npa_main
  public :: npa_read_nbodata

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine npa_main(denskern,rep,ham,elements,ngwf_basis,pnao_tr_output)

    !======================================================================!
    ! This subroutine performs Natural population analysis using the       !
    ! NGWFs as the set of local orbitals.                                  !
    ! Follows method described in J. Chem. Phys. vol 83,735 (1985)         !
    !                                                                      !
    ! http://www.chem.wisc.edu/~nbo5/                                      !
    !----------------------------------------------------------------------!
    ! Written by Louis Lee 26-May-2011                                     !
    ! Currently only works for NMB (no NMB/NRB splitting yet)              !
    ! SPAM3 all the way up till before W(WSW)^(-1/2), DEM thereafter       !
    !======================================================================!

    use comms, only: comms_abort, comms_barrier, comms_bcast, comms_reduce, &
         pub_on_root, pub_my_node_id, pub_root_node_id, pub_total_num_nodes
    use constants, only: DP, stdout
    use dense, only: DEM, dense_convert, dense_copy, dense_create, &
         dense_destroy, dense_eigensolve, dense_get_element, dense_product
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use ngwf_data, only: GTO_SET
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use parallel_strategy, only: pub_orig_atom, pub_distr_atom, &
         pub_first_atom_on_node, pub_num_atoms_on_node
    use rundat, only: print_qc, pub_rootname, pub_nbo_aopnao_scheme, &
         pub_nbo_write_species, pub_nbo_ngwf_label, pub_nbo_init_lclowdin, &
         pub_nbo_write_lclowdin, pub_nbo_write_npacomp, pub_nbo_scale_dm
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_axpy, sparse_copy, sparse_create, &
         sparse_destroy, sparse_product, sparse_scale, sparse_transpose
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_open_unit_check, utils_close_unit_check
    use wrappers, only: wrappers_dsygv_lt
!CW
    use rundat, only : pub_dmft_nbo
!END CW

    implicit none

    ! lpl: Inputs
    ! lpl: Densiy kernel input
    type(SPAM3),    intent(inout) :: denskern(pub_cell%num_spins)
    type(NGWF_REP), intent(in) :: rep
    type(NGWF_HAM), intent(in) :: ham
    ! lpl: Cell elements info
    type(ELEMENT),    intent(in) :: elements(pub_cell%nat)
    ! lpl: NGWF Basis info
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    ! lpl: Optional pnao_tr_output to be passed to properties_mod
    !      for NGWF character analysis
    type(SPAM3), optional, intent(inout) :: pnao_tr_output

    ! lpl: Internal variables
    ! lpl: NGWF labelling
    integer, allocatable :: ngwf_lm_label(:) ! NGWF false lm labels
    integer, allocatable :: atom_order(:)    ! NGWF atomic centers

    ! lpl: Sparse matrices (5 + is)
    type(SPAM3), allocatable :: pure_dm(:) ! Purified density matrix

    type(SPAM3) :: init_tr, aopnao_tr ! tr's we need to store
    type(SPAM3) :: d_tr  ! Temporary storage for intermediate tr

    ! lpl: Sparse buffers
    ! lpl: Used only at the beginning
    type(SPAM3) :: sbuffer_sks
    type(SPAM3) :: sbuffer_s

    ! lpl: Peak usage = 3*DEM + is*DEM + 3*D + 1*SKS + 1*S
    ! lpl: Used throughout the code
    type(SPAM3) :: sbuffer_d1, sbuffer_d2
    type(SPAM3) :: init_dm, init_ovlp

    ! lpl: Full dense buffers (4 + is)
    type(DEM), allocatable :: dbuffer_1(:)
    type(DEM) :: dbuffer_2, dbuffer_3
    type(DEM) :: full_tr

    ! lpl: Column vector buffers
    real(kind=DP), allocatable :: col_buffer(:)
    ! lpl: PNAO weights (tracked throughout code)
    real(kind=DP), allocatable :: pnao_w(:)

    ! lpl: lm-label & NMB/NRB list variables
    logical, allocatable :: is_nmb(:,:) ! Logical of NGWFs that are NMBs
    integer, allocatable :: nmb_list(:), nrb_list(:)  ! List of NMB/NRBs
    integer, allocatable :: high_occ_list(:), low_occ_list(:)

    character(len=256), allocatable :: ngwf_label_option(:)
    integer, allocatable   :: species_ngwf_label(:,:)
    logical, allocatable   :: species_is_nmb(:,:)
    type(GTO_SET), allocatable :: gbasis(:)

    ! lpl: GENNBO labels sorted according to NGWF sequence
    integer :: nbo_lm_label(25)

    ! lpl: Misc scalars
    integer :: is
    integer :: iat, iat_ngwf, orig_iat
    integer :: it, lc_it, num_ngwfs
    integer :: irow, icol
    integer :: num_nmb, num_nrb
    integer :: num_high_occ, num_low_occ

    integer :: list_nat, list_ngwfnum

    real(kind=DP) :: mtx_el, q_total, q_dev

    ! NPA matrix & partial matrix output list
    real(kind=DP), allocatable :: q_atom(:,:)
    integer, allocatable :: atom_list(:), ngwf_list(:), ngwf_orig_list(:)

    integer :: ierr = 0
    integer, parameter :: max_ngwf_label = 95

    ! lpl: Threshold for WSW eigenvalue before program gives up
    !      Since 'zero' in S is constantly observed to be ~<1e-14
    !      set it 100x larger
    ! lpl: 23092011 - W_i < w_min = w_min in wwsw_tr
    real(kind=DP), parameter :: s_threshold = 1.0e-12_DP
    real(kind=DP), parameter :: w_min = 1.0e-4_DP
    ! lpl: Threshold for scaling density matrix to integral charge
    !      J.Chem.Phys.83,735 mentioned tolerence check of 10e-6
    real(kind=DP), parameter :: qint_thresh = 1.0e-6_DP
    ! lpl: Occ thresh for defining low-occ NRBs
    real(kind=DP), parameter :: occ_thresh = 1.0e-4_DP

    ! lpl: Debug stuff
    integer :: debug_output_unit
    character(len=256) :: debug_filename

    ! lpl: Start timer
    call timer_clock('npa_main',1)

#ifdef DEBUG
    if(pub_on_root) then
       write(stdout,'(a)') 'DEBUG: Entering properties_nat_popn_analysis'

       debug_filename = trim(pub_rootname)//'_nbo_DEBUG.dat'
       debug_filename = adjustl(debug_filename)
       debug_output_unit = utils_unit()
       open(unit=debug_output_unit,form="formatted", &
            file=trim(debug_filename),action="write",iostat=ierr)
       call utils_open_unit_check('properties_nat_popn_analysis', &
            'debug_output_unit',ierr)
    end if
#endif

  if(pub_on_root) write(stdout,'(a)') &
       '======= NBO Output (Version 23-09-2011) ========'

 !--------------------------------------------------------------------------!
 ! Matrices initialization                                                  !
 !--------------------------------------------------------------------------!
    ! lpl: Transformation matrices & buffers
    init_tr%structure='D'        ! lpl: init_tr stores initial NGWF --> AO
    call sparse_create(init_tr)  !      transformation (lclowdin or norm)
    aopnao_tr%structure='D'        ! lpl: aopnao_tr stores PNAO transform
    call sparse_create(aopnao_tr)  !      needed for plotting
    d_tr%structure='D'           ! lpl: d_tr acts as buffer to store all
    call sparse_create(d_tr)     !      atom-local transformations

    ! lpl: Dense matrices (used throughout code)
    allocate(dbuffer_1(pub_cell%num_spins))
    call utils_alloc_check('npa_main','dbuffer_1',ierr)
    do is=1,pub_cell%num_spins
       call dense_create(dbuffer_1(is),ngwf_basis%num,ngwf_basis%num)
    end do
    call dense_create(dbuffer_2,ngwf_basis%num,ngwf_basis%num)
    call dense_create(dbuffer_3,ngwf_basis%num,ngwf_basis%num)
    call dense_create(full_tr,ngwf_basis%num,ngwf_basis%num)

    ! lpl: Column vectors/buffers
    allocate(col_buffer(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('npa_main','col_buffer',ierr)
    allocate(pnao_w(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('npa_main','pnao_w',ierr)

    ! lpl: NGWF lm-labelling and NMB/NRB indexing arrays
    allocate(is_nmb(pub_cell%nat,ngwf_basis%max_on_atom),stat=ierr)
    call utils_alloc_check('npa_main','is_nmb',ierr)
    is_nmb = .false.

    allocate(nmb_list(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('npa_main','nmb_list',ierr)
    allocate(nrb_list(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('npa_main','nrb_list',ierr)

    allocate(high_occ_list(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('npa_main','high_occ_list',ierr)
    allocate(low_occ_list(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('npa_main','low_occ_list',ierr)

    allocate(ngwf_label_option(pub_cell%num_species),stat=ierr)
    call utils_alloc_check('npa_main','ngwf_label_option',ierr)

    allocate(species_ngwf_label(pub_cell%num_species, &
         0:ngwf_basis%max_on_atom),stat=ierr)
    call utils_alloc_check('npa_main','species_ngwf_label',ierr)

    allocate(species_is_nmb(pub_cell%num_species,ngwf_basis%max_on_atom), &
         stat=ierr)
    call utils_alloc_check('npa_main','species_is_nmb',ierr)

    allocate(atom_order(ngwf_basis%num))
    call utils_alloc_check('npa_main','atom_order',ierr)
    atom_order=0

    allocate(ngwf_lm_label(ngwf_basis%num))
    call utils_alloc_check('npa_main','ngwf_lm_label',ierr)
    ngwf_lm_label=0

    data nbo_lm_label(:)  /   1,152,153,151, &
         251,253,255,252,254,357,355,353,351,352,354,356, &
         459,457,455,453,451,452,454,456,458 /

    ! Partial matrix selection array & NPA analysis
    allocate(q_atom(pub_cell%nat,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('npa_main','q_atom',ierr)
    q_atom = 0.0_DP

    allocate(atom_list(pub_cell%nat),stat=ierr)
    call utils_alloc_check('npa_main','atom_list',ierr)
    atom_list = -1   ! Set for error checking

    allocate(ngwf_list(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('npa_main','ngwf_list',ierr)
    ngwf_list = -1   ! Set for error checking

    allocate(ngwf_orig_list(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('npa_main','ngwf_orig_list',ierr)
    ngwf_orig_list = -1   ! Set for error checking

    ! lpl: Matrices neede in only part of the code
    ! lpl: Sparse buffers for atom-local transformation inputs
    sbuffer_d1%structure='D'       ! lpl: sbuffer_d1 & sbuffer_d2 stores
    call sparse_create(sbuffer_d1) !      only block-diagonal elements from
    sbuffer_d2%structure='D'       !      dense P and S required for all
    call sparse_create(sbuffer_d2) !      atom-local subroutines

    ! lpl: P & S after init_tr (used throughout code)
    init_dm%structure='SKS'
    call sparse_create(init_dm)
    init_ovlp%structure='S'
    call sparse_create(init_ovlp)

    call comms_barrier

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

 !--------------------------------------------------------------------------!
 ! NGWF lm-labelling and NMB/NRB indexing                                   !
 !--------------------------------------------------------------------------!

    ! lpl: Find NGWF labelling option for all species
    if(pub_on_root) then

       call internal_ngwf_sto3g_create

       ngwf_label_option  = 'ERR'     ! Default if nothing is read
       species_ngwf_label = 0
       species_is_nmb     = .false.

       ! lpl: Find NGWF labelling option for this species
       do it=1,pub_cell%num_species
          ! lpl: Find e.g. of this species
          do orig_iat=1,pub_cell%nat
             if(pub_nbo_ngwf_label(it)(1:4) &
                  == elements(orig_iat)%species_id) then

                ngwf_label_option(elements(orig_iat)%species_number)(:) &
                     = pub_nbo_ngwf_label(it)(5:)
                ! lpl: species_ngwf_label(0) stores e.g. orig_iat of species
                species_ngwf_label(elements(orig_iat)%species_number,0) &
                     = orig_iat
                exit

             end if
          end do
       end do

       ! lpl: Get NGWF labels for each species
       do it=1,pub_cell%num_species

          if(ngwf_label_option(it) == 'ERR') then ! lpl:  Error checking
             write(stdout,'(a)') 'ERROR: Labelling error in ngwf_label_option'
             ierr = 1
          else if(ngwf_label_option(it) == 'AUTO') then ! lpl: If 'AUTO'
             call internal_generatelabel_auto(it,species_ngwf_label(it,0:), &
                  species_is_nmb(it,:),gbasis)
          else ! If user-specified label
             call internal_generatelabel_manual(it,species_ngwf_label(it,0:), &
                  species_is_nmb(it,:),ngwf_label_option(it))
          end if

       end do

       call internal_ngwf_sto3g_destroy

       ! lpl: Assign lm-labels to each atom based on species lm-label
       ngwf_lm_label = 0
       atom_order    = 0
       is_nmb = .false.
       do iat=1,pub_cell%nat
          orig_iat  = pub_orig_atom(iat)

          do lc_it=1,ngwf_basis%num_on_atom(iat)
             ! lpl: Assign NGWF lm-label
             it = lc_it + ngwf_basis%first_on_atom(iat) - 1
             ngwf_lm_label(it) &
                  = species_ngwf_label(elements(orig_iat)%species_number,lc_it)
             ! lpl: Designate NMBs
             is_nmb(iat,lc_it) &
                  = species_is_nmb(elements(orig_iat)%species_number,lc_it)
             atom_order(it) = orig_iat   ! lpl: Assign NGWF atomic centers

          end do

       end do

       ! lpl: Obtain global NMB/NRB lists
       num_nrb = 0
       num_nmb = 0
       nmb_list = -1   ! lpl: Set mapping to -1 for error checking
       nrb_list = -1

       ! Get list of NMB and NRB NGWF indices
       do it=1,ngwf_basis%num
          iat      = ngwf_basis%atom_of_func(it)
          iat_ngwf = it - ngwf_basis%first_on_atom(iat) + 1

          ! Get NGWF to NMB/NRB sub-matrix indices
          if(is_nmb(iat,iat_ngwf)) then ! NMB
             num_nmb = num_nmb + 1
             nmb_list(num_nmb) = it
          else                          ! NRB
             num_nrb = num_nrb + 1
             nrb_list(num_nrb) = it
          end if
       end do

    end if ! END if(pub_on_root)

    if(ierr == 1) then
       write(stdout,'(a)') &
            'ERROR: Something broke in the NGWF labelling routine.'
       call comms_abort
    end if

    ! lpl: Broadcast relevant lists/mappings
    call comms_bcast(pub_root_node_id,ngwf_lm_label)
    call comms_bcast(pub_root_node_id,atom_order)

    ! lpl: For some reason bcast of logical vectors didn't work
    !      so do it explicitly for each element
    do iat=1,pub_cell%nat
       do it=1,ngwf_basis%max_on_atom
          call comms_bcast(pub_root_node_id,is_nmb(iat,it))
       end do
    end do

    call comms_bcast(pub_root_node_id,nmb_list)
    call comms_bcast(pub_root_node_id,nrb_list)
    call comms_bcast(pub_root_node_id,num_nmb)
    call comms_bcast(pub_root_node_id,num_nrb)

    call comms_barrier

    ! lpl: Check that all NGWFs are being accounted for
    if(num_nrb + num_nmb /= ngwf_basis%num) then
       write(stdout,'(a,I5,a,I5,a,I5,a)') &
           'ERROR: # NRB (',num_nrb,') + # NMB (',num_nmb,') &
                 &/= # NGWFs (',ngwf_basis%num,')'
       call comms_abort
    end if

    ! lpl: Get NGWF labelling in 'original' input order (inter-atomically
    !      sorted) - makes NBO printing easier when reading FILE.47
    icol = 0
    do orig_iat=1,pub_cell%nat
       iat = pub_distr_atom(orig_iat)

       num_ngwfs = ngwf_basis%num_on_atom(iat)
       do lc_it=1,num_ngwfs
          it = lc_it + ngwf_basis%first_on_atom(iat) - 1
          icol = icol + 1
          ngwf_orig_list(icol) = it
       end do
    end do
    if(icol /= ngwf_basis%num) then
       write(stdout,'(a)') 'ERROR: Iterator /= ngwf_basis%num'
       call comms_abort
    end if

    call comms_barrier

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

 !--------------------------------------------------------------------------!
 ! K & S initialization                                                     !
 !--------------------------------------------------------------------------!

    call internal_pure_dm_init

!CW
    if(.not.pub_dmft_nbo)then
     call internal_get_puredm(pure_dm(:),denskern(:),rep)
    else
     do is=1,pub_cell%num_spins
       call sparse_copy(pure_dm(is),denskern(is))
     end do
    endif
!END CW

    ! Reduce to spinless into init_dm
    do is=1,pub_cell%num_spins
       call sparse_axpy(init_dm,pure_dm(is),1.0_DP)
    end do
    ! Get overlap matrix init_ovlp
    call sparse_copy(init_ovlp,rep%overlap)

    call internal_pure_dm_destroy ! Save memory especially in large systems

    call comms_barrier

#ifdef DEBUG
    call dense_convert(dbuffer_3,init_dm)
    call internal_testprint_dem(debug_output_unit,dbuffer_3, &
         'SPINLESS DM')
#endif

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

 !--------------------------------------------------------------------------!
 ! Atom-locacl Lowdin transformation                                        !
 !--------------------------------------------------------------------------!

    if(pub_nbo_init_lclowdin) then
       call internal_lclowdin_tr(init_tr,rep%overlap)
    else
       call internal_norm_tr(init_tr,rep%overlap)
    end if

    ! lpl: Update init_dm and init_ovlp - seperate init_tr into dm
    !      and ovlp
    ! lpl: sbuffers only needed here
    sbuffer_sks%structure = 'SKS'
    call sparse_create(sbuffer_sks)
    sbuffer_s%structure = 'S'
    call sparse_create(sbuffer_s)

    call sparse_transpose(sbuffer_d1,init_tr)

    call sparse_product(sbuffer_sks,init_dm,init_tr)
    call sparse_product(init_dm,sbuffer_d1,sbuffer_sks)

    call sparse_product(sbuffer_s,init_ovlp,init_tr)
    call sparse_product(init_ovlp,sbuffer_d1,sbuffer_s)

    ! lpl: Not needed anymore
    call sparse_destroy(sbuffer_s)
    call sparse_destroy(sbuffer_sks)

#ifdef DEBUG
    call dense_convert(dbuffer_1(1),init_dm)
    call dense_convert(dbuffer_2,init_ovlp)
    call internal_testprint_dem(debug_output_unit,dbuffer_1(1),'INIT DM')
    call internal_testprint_dem(debug_output_unit,dbuffer_2,'INIT OVLP')
#endif

    call comms_barrier

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

 !--------------------------------------------------------------------------!
 ! AO --> PNAO transformation                                               !
 !--------------------------------------------------------------------------!

    ! lpl: pnao_w stores dioganal PNAO occupancies for the next step
    pnao_w = 0.0_DP

    ! lpl: Atom-local transformation needs only block-diagonal info
    call sparse_copy(sbuffer_d1,init_dm)
    call sparse_copy(sbuffer_d2,init_ovlp)

!    ! lpl: If sparse_copy gives error then use this (commented out)
!       call dense_convert(dbuffer_1(1),init_dm)
!       call dense_convert(dbuffer_2,init_ovlp)
!
!       call dense_convert(sbuffer_d1,dbuffer_1(1))
!       call dense_convert(sbuffer_d2,dbuffer_2)
!    ! lpl: If sparse_copy gives error then use this (commented out)


    ! lpl: Select initial AO to PNAO transformation scheme
    !      NOTE: 'DIAGONALIZATION' and 'NONE' for testing purposes only
    !      but might become useful later. Should never be used
    select case(pub_nbo_aopnao_scheme)
       case ('ORIGINAL')
          if(pub_on_root) write(stdout,'(a)') &
             ' AOPNAO Scheme: ORIGINAL'
          call internal_aopnao_tr(aopnao_tr,pnao_w,sbuffer_d1, &
               sbuffer_d2,ngwf_lm_label,is_nmb,'Full')
       case ('DIAGONALIZATION')
          if(pub_on_root) write(stdout,'(a)') &
             ' AOPNAO Scheme: DIAGONALIZATION'
          call internal_diag_tr(aopnao_tr,pnao_w,sbuffer_d1, &
               sbuffer_d2,is_nmb,'Full')
       case ('NONE')
          if(pub_on_root) write(stdout,'(a)') &
             ' AOPNAO Scheme: NONE'
          call internal_sparse_id(aopnao_tr)
          ! lpl: Get PNAO weights directly into col_buffer
          do icol=1,ngwf_basis%num
             call dense_get_element(pnao_w(icol), &
                  dbuffer_1(1),icol,icol)
          end do
       case default
          if(pub_on_root) write(stdout,'(a)') 'ERROR: Invalid &
               &nbo_aopnao_scheme specified'
          call comms_abort
    end select

    ! lpl: Update transformations for next step
    call dense_convert(dbuffer_3,aopnao_tr)

#ifdef DEBUG
    call internal_testprint_dem(debug_output_unit,dbuffer_3,'PNAO TR')
#endif

    call internal_update_tr(init_dm,init_ovlp,dbuffer_3, &
         dbuffer_1(1),dbuffer_2,full_tr,first_tr=.true.)

#ifdef DEBUG
    call internal_testprint_dem(debug_output_unit,dbuffer_1(1), &
         'PNAO DM')
    call internal_testprint_dem(debug_output_unit,dbuffer_2, &
         'PNAO OVLP')

    if(pub_on_root) then
       write(debug_output_unit,'(a)') 'PNAO Occ 1'
       write(debug_output_unit,'(a)') '----------------------------'
       write(debug_output_unit,'(a)') '  NGWF #          Weight    '
       write(debug_output_unit,'(a)') '----------------------------'
       do it=1,ngwf_basis%num
          write(debug_output_unit,'(1x,I5,1x,f14.7)') it,pnao_w(it)
       end do
       write(debug_output_unit,'(a)') '----------------------------'
    end if
#endif

    call comms_barrier

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

 !--------------------------------------------------------------------------!
 ! NMB Weighted Orthogonalisation                                           !
 !--------------------------------------------------------------------------!
    ! lpl: Transformation in dbuffer_3
    !      Overlap mat in dbuffer_2 gets used as a buffer
    call internal_wwsw_tr(dbuffer_3,dbuffer_2, &
         pnao_w,nmb_list,num_nmb,dbuffer_1(1),col_buffer)

#ifdef DEBUG
    call internal_testprint_dem(debug_output_unit,dbuffer_3,'NMB WWSW TR')
#endif

    ! lpl: Update mat, dbuffer_2 contains updated ovlp from init_ovlp*full_tr
    !      dbuffer_1 contains updated dm from init_dm*full_tr etc.
    call internal_update_tr(init_dm,init_ovlp,dbuffer_3, &
         dbuffer_1(1),dbuffer_2,full_tr)

#ifdef DEBUG
    call internal_testprint_dem(debug_output_unit,dbuffer_2,'NMB WWSW OVLP')
#endif

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

 !--------------------------------------------------------------------------!
 ! Skip all NRB-related transformations if minimal basis is given           !
 !--------------------------------------------------------------------------!

    if(num_nrb > 0) then

    !-----------------------------------------------------------------------!
    ! NMB/NRB Schmidt Orthogonalisation                                     !
    !-----------------------------------------------------------------------!

       ! lpl: NOTE - For some weird reason GENNBO does not use the normalised
       !      basis for lm-averaging
       call internal_schmidt_tr(dbuffer_3,dbuffer_2, &
            nrb_list(1:num_nrb),num_nrb,nmb_list(1:num_nmb),num_nmb, &
            col_buffer,norm=.false.)

#ifdef DEBUG
       call internal_testprint_dem(debug_output_unit,dbuffer_3, &
            'NMB/NRB SCH TR')
#endif

       call internal_update_tr(init_dm,init_ovlp,dbuffer_3, &
            dbuffer_1(1),dbuffer_2,full_tr)

#ifdef DEBUG
       call internal_testprint_dem(debug_output_unit,dbuffer_1(1), &
            'NMB/NRB SCH DM')
       call internal_testprint_dem(debug_output_unit,dbuffer_2, &
            'NMB/NRB SCH OVLP')
#endif

    !-----------------------------------------------------------------------!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !-----------------------------------------------------------------------!

    !-----------------------------------------------------------------------!
    ! NRB lm-averaging                                                      !
    !-----------------------------------------------------------------------!

       ! lpl: pnao_w containes symmetry-averaged occ for NMBs and Schmidt-
       !      ortho occs for NRBs that will soon be overwritten by lm-avg
       ! lpl: Need only atom-block info
       call dense_convert(sbuffer_d1,dbuffer_1(1))
       call dense_convert(sbuffer_d2,dbuffer_2)

       ! lpl: Select NRB 'PNAO' transformation scheme
       !      NOTE: 'DIAGONALIZATION' and 'NONE' for testing purposes only
       !      but might become useful later. Should never be used
       select case(pub_nbo_aopnao_scheme)
          case ('ORIGINAL')
             call internal_aopnao_tr(d_tr,pnao_w,sbuffer_d1, &
                  sbuffer_d2,ngwf_lm_label,is_nmb,'NRB')
          case ('DIAGONALIZATION')
             call internal_diag_tr(d_tr,pnao_w,sbuffer_d1, &
                  sbuffer_d2,is_nmb,'NRB')
          case ('NONE')
             call internal_sparse_id(d_tr)
             ! lpl: Re-obtain pnao_w that has changed due to schmidt_tr
             do icol=1,ngwf_basis%num
                call dense_get_element(pnao_w(icol), &
                     dbuffer_1(1),icol,icol)
             end do
          case default
             if(pub_on_root) write(stdout,'(a)') 'ERROR: Invalid &
                  &nbo_aopnao_scheme specified'
             call comms_abort
       end select

       ! lpl: Transfer transformation to dbuffer_3
       call dense_convert(dbuffer_3,d_tr)

#ifdef DEBUG
       call internal_testprint_dem(debug_output_unit,dbuffer_3,'NRB LM TR')
#endif

       call internal_update_tr(init_dm,init_ovlp,dbuffer_3, &
            dbuffer_1(1),dbuffer_2,full_tr)

#ifdef DEBUG
       call internal_testprint_dem(debug_output_unit,dbuffer_1(1), &
            'NRB LM DM')
       call internal_testprint_dem(debug_output_unit,dbuffer_2, &
            'NRB LM OVLP')

       if(pub_on_root) then
          write(debug_output_unit,'(a)') 'PNAO Occ 2'
          write(debug_output_unit,'(a)') '----------------------------'
          write(debug_output_unit,'(a)') '  NGWF #          Weight    '
          write(debug_output_unit,'(a)') '----------------------------'
          do it=1,ngwf_basis%num
             write(debug_output_unit,'(1x,I5,1x,f14.7)') it,pnao_w(it)
          end do
          write(debug_output_unit,'(a)') '----------------------------'
       end if
#endif

    !-----------------------------------------------------------------------!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !-----------------------------------------------------------------------!

    !-----------------------------------------------------------------------!
    ! High occupancy NRB Weighted Orthogonalisation                         !
    !-----------------------------------------------------------------------!

       ! Identify low-occupancy NRBs
       high_occ_list = 0
       low_occ_list  = 0

       num_high_occ = 0
       num_low_occ  = 0

       do icol=1,num_nrb
          if(pnao_w(nrb_list(icol)) >= occ_thresh) then
             num_high_occ = num_high_occ + 1
             high_occ_list(num_high_occ) = nrb_list(icol)
          else
             num_low_occ = num_low_occ + 1
             low_occ_list(num_low_occ) = nrb_list(icol)

             ! lpl: Set equal weightings for Lowdin orthogonlisation
             !      Might do away with this in a seperate Lowdin subroutine
             pnao_w(nrb_list(icol)) = 1.0_DP
          end if
       end do
       ! lpl: Sanity check
       if(pub_on_root) then
          if( num_high_occ + num_low_occ /= num_nrb ) then
             write(stdout,'(a)') 'ERROR: num_high_occ + num_low_occ /= &
                  &num_nrb in npa_main'
             call comms_abort
          end if
       end if

       call internal_wwsw_tr(dbuffer_3,dbuffer_2, &
            pnao_w,high_occ_list,num_high_occ,dbuffer_1(1),col_buffer)

#ifdef DEBUG
       call internal_testprint_dem(debug_output_unit,dbuffer_3,'NRB WWSW TR')
#endif

       call internal_update_tr(init_dm,init_ovlp,dbuffer_3, &
            dbuffer_1(1),dbuffer_2,full_tr)

#ifdef DEBUG
       call internal_testprint_dem(debug_output_unit,dbuffer_2, &
            'NRB WWSW OVLP')
#endif

    !-----------------------------------------------------------------------!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !-----------------------------------------------------------------------!

    !-----------------------------------------------------------------------!
    ! Schmidt-orthogonalise low-occupancy NRBs                              !
    !-----------------------------------------------------------------------!

       if(num_low_occ > 0) then
          call internal_schmidt_tr(dbuffer_3,dbuffer_2, &
               low_occ_list(1:num_low_occ),num_low_occ, &
               high_occ_list(1:num_high_occ),num_high_occ, &
               col_buffer,norm=.false.)

#ifdef DEBUG
          call internal_testprint_dem(debug_output_unit,dbuffer_3, &
               'NRB SCH TR')
#endif

          call internal_update_tr(init_dm,init_ovlp,dbuffer_3, &
               dbuffer_1(1),dbuffer_2,full_tr)

#ifdef DEBUG
          if(pub_on_root .and. num_low_occ > 0) &
               write(debug_output_unit,'(a,I5)') &
                    '# low-occ NRB: ',low_occ_list

          call internal_testprint_dem(debug_output_unit,dbuffer_2, &
               'LOW OCC NRB OVLP')

          if(pub_on_root) then
             write(debug_output_unit,'(a)') 'PNAO Occ 3'
             write(debug_output_unit,'(a)') '----------------------------'
             write(debug_output_unit,'(a)') '  NGWF #          Weight    '
             write(debug_output_unit,'(a)') '----------------------------'
             do it=1,ngwf_basis%num
                write(debug_output_unit,'(1x,I5,1x,f14.7)') it,pnao_w(it)
             end do
             write(debug_output_unit,'(a)') '----------------------------'
          end if
#endif
       end if

    !-----------------------------------------------------------------------!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !-----------------------------------------------------------------------!

    !-----------------------------------------------------------------------!
    ! Lowdin-orthogonalise low-occupancy NRBs                               !
    !-----------------------------------------------------------------------!

       if(num_low_occ > 0) then
          ! lpl: Transformation in dbuffer_3
          call internal_wwsw_tr(dbuffer_3,dbuffer_2, &
               pnao_w,low_occ_list,num_low_occ,dbuffer_1(1),col_buffer)

#ifdef DEBUG
          call internal_testprint_dem(debug_output_unit,dbuffer_3, &
               'NRB LOWDIN TR')
#endif

          call internal_update_tr(init_dm,init_ovlp,dbuffer_3, &
               dbuffer_1(1),dbuffer_2,full_tr)

#ifdef DEBUG
          call internal_testprint_dem(debug_output_unit,dbuffer_1(1), &
               'NRB LOWDIN DM')
          if(pub_on_root) then
             q_total = 0.0_DP
             do icol=1,ngwf_basis%num
                call dense_get_element(mtx_el,dbuffer_1(1),icol,icol)
                q_total = q_total + mtx_el
             end do
             write(stdout,'(a,F14.7)') 'Tr[P] = ',q_total
             q_total = 0.0_DP
          end if
          call comms_barrier
          call internal_testprint_dem(debug_output_unit,dbuffer_2, &
               'NRB LOWDIN OVLP')
#endif
       end if

    !-----------------------------------------------------------------------!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !-----------------------------------------------------------------------!

    end if ! END if(num_nrb > 0)

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

 !--------------------------------------------------------------------------!
 ! Rediagonalisation (Since O_W destroys 'natural' character)               !
 !--------------------------------------------------------------------------!

   ! lpl: 12092011 - Rediagonalisation to be performed even on minimal
   !      basis since O_W apparently destroys 'natural' character
   ! lpl: Atom-local transformation needs only block-diagonal info
    call dense_convert(sbuffer_d1,dbuffer_1(1))
    call dense_convert(sbuffer_d2,dbuffer_2)

    ! lpl: Select initial AO to PNAO transformation scheme
    ! lpl: Final step doesn't require pnao_w anymore so ignore
    ! lpl: Commented out re-printing of AOPNAO scheme
    select case(pub_nbo_aopnao_scheme)
       case ('ORIGINAL')
!          if(pub_on_root) write(stdout,'(a)') &
!             ' AOPNAO Scheme: ORIGINAL'
          call internal_aopnao_tr(d_tr,col_buffer,sbuffer_d1, &
               sbuffer_d2,ngwf_lm_label,is_nmb,'Full')
       case ('DIAGONALIZATION')
!          if(pub_on_root) write(stdout,'(a)') &
!             ' AOPNAO Scheme: DIAGONALIZATION'
          call internal_diag_tr(d_tr,col_buffer,sbuffer_d1, &
               sbuffer_d2,is_nmb,'Full')
       case ('NONE')
!          if(pub_on_root) write(stdout,'(a)') &
!             ' AOPNAO Scheme: NONE'
          call internal_sparse_id(d_tr)
       case default
          if(pub_on_root) write(stdout,'(a)') 'ERROR: Invalid &
               &nbo_aopnao_scheme specified'
          call comms_abort
    end select

    ! lpl: Update transformations for next step
    call dense_convert(dbuffer_3,d_tr)
    call internal_update_tr(init_dm,init_ovlp,dbuffer_3, &
         dbuffer_1(1),dbuffer_2,full_tr)

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

    ! lpl: Not needed anymore
    call sparse_destroy(init_ovlp)
    call sparse_destroy(init_dm)

    call sparse_destroy(sbuffer_d2)
    call sparse_destroy(sbuffer_d1)

#ifdef DEBUG
   ! Transformation w/o init_tr (i.e. full_tr matrix)
    if(pub_on_root) then
       write(debug_output_unit,'(a)') ' NAO in run-time AO basis'
       write(debug_output_unit,'(a)') '---------------------------------'
       write(debug_output_unit,'(a)') ' orig_iat NGWF#'
       write(debug_output_unit,'(a)') '---------------------------------'
       do orig_iat=1,pub_cell%nat
          iat = pub_distr_atom(orig_iat)
          num_ngwfs = ngwf_basis%num_on_atom(iat)
          it = ngwf_basis%first_on_atom(iat)

          write(debug_output_unit,'(1x,I5)',advance='no') orig_iat

          do lc_it=1,num_ngwfs
             if(lc_it > 1) write(debug_output_unit,'(6x,a)',advance='no') ''
             write(debug_output_unit,'(I3)',advance='no') lc_it
             icol = lc_it + it - 1
             do irow=1,ngwf_basis%num
                call dense_get_element(mtx_el,full_tr,irow,icol)
                write(debug_output_unit,'(1x,f16.9)',advance='no') mtx_el
             end do
             write(debug_output_unit,'(a)') ''
          end do
       end do
       write(debug_output_unit,'(a)') '--------------------------------'

       q_total = 0.0_DP
       do icol=1,num_nrb
          call dense_get_element(mtx_el,dbuffer_1(1), &
               nrb_list(icol),nrb_list(icol))
          q_total = q_total + mtx_el
       end do
       write(debug_output_unit,'(a,f11.9)') 'NRB Population = ',q_total
       q_total = 0.0_DP
    end if
#endif

 !--------------------------------------------------------------------------!
 ! Full AO --> NAO transformation on original (spin-polarized)              !
 ! density matrix                                                           !
 !--------------------------------------------------------------------------!

    ! lpl: Get full spinless AO --> NAO transformation
    ! lpl: Now combine init_tr to get full transformation
    call dense_copy(dbuffer_1(1),full_tr)
    call dense_convert(dbuffer_2,init_tr)
    call dense_product(full_tr,dbuffer_2,dbuffer_1(1))   ! init_tr*full_tr

    ! lpl: Get pure_dm again
    call internal_pure_dm_init

!CW
    if(.not.pub_dmft_nbo)then
     call internal_get_puredm(pure_dm(:),denskern(:),rep)
    else
     do is=1,pub_cell%num_spins
       call sparse_copy(pure_dm(is),denskern(is))
     end do
    endif
!ENDCW

    ! lpl: Transform density matrix and store in dbuffer_1(is)
    do is=1,pub_cell%num_spins
       call dense_convert(dbuffer_1(is),pure_dm(is))
       call internal_dense_transform(dbuffer_1(is),full_tr,dbuffer_2)
    end do

    ! lpl: Write NPA information
    if(pub_on_root) then
       if(pub_nbo_write_npacomp) then
          write(stdout,'(a)') '================================================'
          write(stdout,'(a)') '               Natural Population               '
          write(stdout,'(a)') '------------------------------------------------'
          write(stdout,'(a)') ' Component charges'
          write(stdout,'(a)') '------------------------------------------------'
          write(stdout,'(a)') '    Atom     NGWF  lm-Label    Population (e)   '
          write(stdout,'(a)') '------------------------------------------------'
       end if
    end if

    ! lpl: Get NPA from density matrix [dbuffer_1(is)]
    do is=1,pub_cell%num_spins
       if(pub_on_root) then

          if(pub_nbo_write_npacomp) then
             write(stdout,'(1x,a,I1)') 'Spin ',is
             write(stdout,'(a)') '----------'
          end if

          do orig_iat=1,pub_cell%nat
            if(pub_nbo_write_npacomp) &
                 write(stdout,'(1x,a4,1x,I5)',advance='no') &
                      elements(orig_iat)%species_id,orig_iat

             iat = pub_distr_atom(orig_iat)
             num_ngwfs = ngwf_basis%num_on_atom(iat)

             q_atom(orig_iat,is) = 0.0_DP
             do lc_it=1,num_ngwfs
                it = lc_it + ngwf_basis%first_on_atom(iat) - 1
                call dense_get_element(mtx_el,dbuffer_1(is),it,it)

                if(pub_nbo_write_npacomp) then
                   if(lc_it > 1) write(stdout,'(11x,a)',advance='no') ''
                   write(stdout,'(4x,I3,4x,I3,9x,f11.7)') &
                        lc_it,ngwf_lm_label(it),mtx_el
                end if

                q_atom(orig_iat,is) = q_atom(orig_iat,is) + mtx_el
             end do

          end do ! END do orig_iat=1,pub_cell%nat
          if(pub_nbo_write_npacomp) write(stdout,'(a)') &
               '------------------------------------------------'

       end if ! END if(pub_on_root)
    end do ! END do is=1,pub_cell%num_spins

    if(pub_on_root) then
       ! lpl: QC test for 1st 3 atoms
       if(print_qc) then
          write(stdout,'(a)') ''
          do orig_iat=1,min(3,pub_cell%nat)
             mtx_el = 0.0_DP
             do is=1,pub_cell%num_spins
                mtx_el = mtx_el + q_atom(orig_iat,is)
             end do

             write(stdout,'(a,I1,a,f14.7)') &
                '<QC>    [atom_',orig_iat,'_natural_population]: ',mtx_el
          end do
       end if

       q_total = 0.0_DP
       write(stdout,'(a)') ' Summary                                        '
       write(stdout,'(a)') '------------------------------------------------'
       write(stdout,'(a)') '   Atom        Population (e)      Charge (e)   '
       write(stdout,'(a)') '------------------------------------------------'
       do orig_iat=1,pub_cell%nat
          mtx_el = 0.0_DP
          do is=1,pub_cell%num_spins
             mtx_el = mtx_el + q_atom(orig_iat,is)
          end do
          write(stdout,'(1x,a4,1x,I5,7x,f11.7,5x,f11.7)') &
               elements(orig_iat)%species_id,orig_iat,mtx_el, &
                    elements(orig_iat)%ion_charge - mtx_el
          q_total = q_total + mtx_el
       end do
       write(stdout,'(a)') '------------------------------------------------'
       write(stdout,'(a,f14.7)') ' Total charge (e): ',q_total
       write(stdout,'(a)') '================================================'

    end if

    call comms_barrier

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

 !--------------------------------------------------------------------------!
 ! Output matrices in the NAO basis into FILE.47                            !
 !--------------------------------------------------------------------------!

    if(pub_on_root) then
       list_nat     = 0
       list_ngwfnum = 0
       q_total = 0.0_DP
       do orig_iat=1,pub_cell%nat
          if(pub_nbo_write_species(orig_iat)) then
             list_nat = list_nat + 1
             atom_list(list_nat) = orig_iat

             do is=1,pub_cell%num_spins
                q_total = q_total + q_atom(orig_iat,is)
             end do

             iat = pub_distr_atom(orig_iat)

             num_ngwfs = ngwf_basis%num_on_atom(iat)
             do lc_it=1,num_ngwfs
                iat_ngwf = lc_it + ngwf_basis%first_on_atom(iat) - 1
                list_ngwfnum = list_ngwfnum + 1
                ngwf_list(list_ngwfnum) = iat_ngwf
             end do
          end if
       end do

       ! lpl: Modified nint to int to round DOWN charges - else we migh have
       !      orbital populations > 2.0 e, which causes GENNBO to terminate
       mtx_el = 1.0_DP*int(q_total)/q_total   ! Scaling parameter
       q_dev  = q_total - 1.0_DP*int(q_total) ! Deviation from integral charge
       if(abs(q_dev) > qint_thresh .and. (list_ngwfnum /= ngwf_basis%num)) then
          write(stdout,'(a,/a)') &
               ' WARNING: Non-integer total charge in partial &
                &density kernel. This will','not be accepted by GENNBO.'

          if(pub_nbo_scale_dm) then ! lpl: Default behaviour (scale DM)
             write(stdout,'(a)') 'Density matrix will be re-scaled by:'
             write(stdout,'(1x,a,f14.7)') 'Scaling pamater : ',mtx_el
             write(stdout,'(1x,a,f14.7)') 'Real charge     : ',q_total
             write(stdout,'(1x,a,f14.7)') 'Scaled charge   : ',mtx_el*q_total
             write(stdout,'(a)') 'Please scale your NBO occupancies &
                  &accordingly'
          else ! lpl: User requested to skip rescaling
             write(stdout,'(a)') 'WARNING: Density matrix re-scaling has been &
                  &manually skipped.'
             mtx_el = 1.0_DP
          end if
       else
          mtx_el = 1.0_DP
       end if

    end if ! END if(pub_on_root)
    call comms_bcast(pub_root_node_id,mtx_el)
    call comms_bcast(pub_root_node_id,list_nat)
    call comms_bcast(pub_root_node_id,list_ngwfnum)
    call comms_bcast(pub_root_node_id,atom_list)
    call comms_bcast(pub_root_node_id,ngwf_list)

    ! lpl: Write FILE.47 in the NAO basis for the potentially partial
    !      system, dictated by ngwf_list(:) and atom_list(:)
    ! lpl: 30082011 - set ortho=.true. (3rd argument)
    call internal_write_nbo(full_tr,'nao',.true.,atom_list,list_nat, &
                            ngwf_list,list_ngwfnum,mtx_el, &
                            dbuffer_1(:),dbuffer_2,dbuffer_3, &
                            tr_dm_input=.true.,tr_ovlp_input=.false., &
                            print_tr_input=.true.,print_47_input=.true.)

   !--------------------------------------------------------------------------!

    ! lpl: Print PNAO matrix (for plotting)
    call sparse_product(d_tr,init_tr,aopnao_tr)
    call dense_convert(full_tr,d_tr)

    ! lpl: Pass to pnao_tr_output if present
    if(present(pnao_tr_output)) call sparse_copy(pnao_tr_output,d_tr)

    call internal_write_nbo(full_tr,'pnao',.false.,atom_list,list_nat, &
                            ngwf_list,list_ngwfnum,mtx_el, &
                            dbuffer_1(:),dbuffer_2,dbuffer_3, &
                            tr_dm_input=.false.,tr_ovlp_input=.false., &
                            print_tr_input=.true.,print_47_input=.false.)

   !--------------------------------------------------------------------------!


    ! lpl: Get lclowdin transformation if this has not already been done
    !      only if nbo_write_lclowdin = T
    if(pub_nbo_write_lclowdin) then
       if(.not. pub_nbo_init_lclowdin) &
            call internal_lclowdin_tr(init_tr,rep%overlap)
       call dense_convert(full_tr,init_tr)

       ! lpl: Reset ngwf_list and atom_list
       !      23/06/11 - now follows orig_iat
       atom_list = -1
       ngwf_list = -1
       it = 0
       icol = 0
       do orig_iat=1,pub_cell%nat
          it = it + 1
          atom_list(it) = orig_iat

          iat = pub_distr_atom(orig_iat)

          num_ngwfs = ngwf_basis%num_on_atom(iat)
          do lc_it=1,num_ngwfs
             iat_ngwf = lc_it + ngwf_basis%first_on_atom(iat) - 1
             icol = icol + 1
             ngwf_list(icol) = iat_ngwf
          end do
       end do

       call comms_barrier

       ! lpl: Print full matrix into FILE.47 for reference
       !      In order to maintain compatibility with GENNBO matrices are
       !      expressed in the lclowdin basis
       call internal_write_nbo(full_tr,'lclowdin',.false., &
            atom_list,pub_cell%nat,ngwf_list,ngwf_basis%num,1.0_DP, &
            dbuffer_1(:),dbuffer_2,dbuffer_3, &
            tr_dm_input=.false.,tr_ovlp_input=.false., &
            print_tr_input=.false.,print_47_input=.true.)
    end if

    ! lpl: Deaellocated pure_dm
    call internal_pure_dm_destroy

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

 !--------------------------------------------------------------------------!
 ! Deallocation & exit subroutines                                          !
 !--------------------------------------------------------------------------!

    deallocate(ngwf_orig_list,stat=ierr)
    call utils_dealloc_check('npa_main','ngwf_orig_list',ierr)

    deallocate(ngwf_list,stat=ierr)
    call utils_dealloc_check('npa_main','ngwf_list',ierr)
    deallocate(atom_list,stat=ierr)
    call utils_dealloc_check('npa_main','atom_list',ierr)
    deallocate(q_atom,stat=ierr)
    call utils_dealloc_check('npa_main','q_atom',ierr)

    ! lpl: 19082011 - Forgot to deallocate in previous revision
    deallocate(atom_order,stat=ierr)
    call utils_dealloc_check('npa_main','atom_order',ierr)
    deallocate(ngwf_lm_label,stat=ierr)
    call utils_dealloc_check('npa_main','ngwf_lm_label',ierr)

    deallocate(ngwf_label_option,stat=ierr)
    call utils_dealloc_check('npa_main','ngwf_label_option',ierr)
    deallocate(species_ngwf_label,stat=ierr)
    call utils_dealloc_check('npa_main','species_ngwf_label',ierr)
    deallocate(species_is_nmb,stat=ierr)
    call utils_dealloc_check('npa_main','species_is_nmb',ierr)

    deallocate(low_occ_list,stat=ierr)
    call utils_dealloc_check('npa_main','low_occ_list',ierr)
    deallocate(high_occ_list,stat=ierr)
    call utils_dealloc_check('npa_main','high_occ_list',ierr)

    deallocate(nrb_list,stat=ierr)
    call utils_dealloc_check('npa_main','nrb_list',ierr)
    deallocate(nmb_list,stat=ierr)
    call utils_dealloc_check('npa_main','nmb_list',ierr)
    deallocate(is_nmb,stat=ierr)
    call utils_dealloc_check('npa_main','is_nmb',ierr)

    deallocate(col_buffer,stat=ierr)
    call utils_dealloc_check('npa_main','col_buffer',ierr)
    deallocate(pnao_w,stat=ierr)
    call utils_dealloc_check('npa_main','pnao_w',ierr)

    ! lpl: Deallocate dense matrices
    call dense_destroy(full_tr)
    call dense_destroy(dbuffer_3)
    call dense_destroy(dbuffer_2)
    do is=pub_cell%num_spins,1,-1
       call dense_destroy(dbuffer_1(is))
    end do
    deallocate(dbuffer_1,stat=ierr)
    call utils_dealloc_check('npa_main','dbuffer_1',ierr)

    ! lpl: Deallocate sparse matrices
    call sparse_destroy(d_tr)
    call sparse_destroy(aopnao_tr)
    call sparse_destroy(init_tr)

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

    if(pub_on_root) write(stdout,'(a)') &
         '================================================'

#ifdef DEBUG
    if(pub_on_root) then
       close(unit=debug_output_unit,iostat=ierr)
       call utils_close_unit_check('properties_nat_popn_analysis', &
            'debug_output_unit',ierr)
       write(stdout,'(a)') 'DEBUG: Leaving properties_nat_popn_analysis'
    end if
#endif

    ! lpl: Stop timer
    call timer_clock('npa_main',2)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains

 !--------------------------------------------------------------------------!
 ! Misc. matrix orperation                                                  !
 ! Mostly to minimize clutter in main code                                  !
 !    1. internal_testprint_dem   - Prints out various matrices. Used only  !
 !                                  in debug statements                     !
 !    2. internal_dense_transform - Transform dense matrices                !
 !                                  B' = (A^dagger)*B*A                     !
 !    3. internal_dense_coldiag   - Column vector --> Diagonal matrix       !
 !    4. internal_dense_id        - Generate identity matrix                !
 !--------------------------------------------------------------------------!

    subroutine internal_testprint_dem(output_unit,mat,mat_name)
      ! lpl: Prints out various matrices. Used in debug statements
      use dense, only: dense_get_element
      implicit none

      type(DEM),    intent(in) :: mat
      character(*), intent(in) :: mat_name
      integer, intent(in) :: output_unit

      integer :: irow, icol

      call comms_barrier
      if(pub_on_root) then
         write(output_unit,'(a,a,a,I5,a,I5,a)') 'MATRIX: ',mat_name, &
                                           ' (',mat%nrows,'x',mat%mcols,')'
         write(output_unit,'(a)') '--------------------'
         do irow=1,mat%nrows
            do icol=1,mat%mcols
               call dense_get_element(mtx_el,mat,irow,icol)
               write(output_unit,'(1x,e14.7)',advance='no') mtx_el
            end do
            write(output_unit,'(a)') ''
         end do
         write(output_unit,'(a)') '--------------------'
      end if
      call comms_barrier

    end subroutine internal_testprint_dem

    subroutine internal_dense_transform(bmat,amat,buffmat)
      ! lpl: Performs B' = A^(dagger)*B*A
      implicit none

      ! Input/Output
      type(DEM), intent(inout) :: bmat    ! Matrix to be transformed
      type(DEM), intent(inout) :: buffmat ! Buffer
      type(DEM), intent(in) :: amat       ! Transformation matrix

      ! lpl:  MT then (T^dagger)MT
      call dense_product(buffmat,bmat,amat)
      call dense_product(bmat,amat,buffmat, &
                         transpose_amat=.true.,transpose_bmat=.false.)

    end subroutine internal_dense_transform

    subroutine internal_dense_coldiag(amat,col_diag,num)
      ! lpl: Gets column list into diagonal matrix form to minimize clutter
      use dense, only: dense_put_element, dense_scale

      implicit none

      type(DEM),  intent(inout) :: amat
      real(kind=DP), intent(in) :: col_diag(:)
      integer,       intent(in) :: num

      integer :: irow, icol

      do irow=1,amat%nrows
         do icol=1,amat%mcols
            call dense_put_element(0.0_DP,amat,irow,icol)
         end do
      end do
      call comms_barrier
      do icol=1,num
         call dense_put_element(col_diag(icol),amat,icol,icol)
      end do
      call comms_barrier

    end subroutine internal_dense_coldiag

    subroutine internal_dense_id(amat,num)
      ! lpl: Constructs identity matrix
      use dense, only: dense_put_element, dense_scale
      implicit none

      type(DEM),  intent(inout) :: amat
      integer,       intent(in) :: num

      integer :: irow, icol

      do irow=1,amat%nrows
         do icol=1,amat%mcols
            call dense_put_element(0.0_DP,amat,irow,icol)
         end do
      end do
      call comms_barrier
      do icol=1,num
         call dense_put_element(1.0_DP,amat,icol,icol)
      end do
      call comms_barrier

    end subroutine internal_dense_id

    subroutine internal_sparse_id(mat)
      ! lpl: Generate identity SPAM3
      use parallel_strategy, only: pub_first_atom_on_node, &
          pub_num_atoms_on_node
      use sparse, only: sparse_node_of_elem, sparse_put_element

      implicit none

      type(SPAM3), intent(inout) :: mat
      integer :: irow, icol

      do icol=1,ngwf_basis%num
         if(sparse_node_of_elem(icol,mat,'C') == pub_my_node_id) then
            do irow=1,ngwf_basis%num ! lpl: Explicitly zero elements
               call sparse_put_element(0.0_DP,mat,irow,icol)
            end do
            call sparse_put_element(1.0_DP,mat,icol,icol)
         end if
      end do

    end subroutine internal_sparse_id

    subroutine internal_update_tr(init_dm,init_ovlp,tr_buffer, &
         dbuffer_dm,dbuffer_ovlp,full_tr,first_tr)
      ! lpl: Updates dm, ovlp, and current accumulated tr 'full_tr'
      implicit none

      type(SPAM3), intent(in) :: init_dm, init_ovlp ! Initial mat
      ! lpl: tr_buffer contains tr that is used as a buffer later on
      type(DEM), intent(inout) :: tr_buffer

      ! lpl: Outputs updated dm, ovlp, and full_tr for next step
      type(DEM), intent(inout) :: dbuffer_dm, dbuffer_ovlp
      type(DEM), intent(inout) :: full_tr ! Full tr to be updated

      ! lpl: If this is the 1st tr just copy tr to full_tr
      logical, optional, intent(in) :: first_tr

      logical :: is_first

      is_first = .false. ! Default behaviour
      if(present(first_tr)) is_first = first_tr

      ! lpl: Update full_tr
      if(is_first) then
         call dense_copy(full_tr,tr_buffer)
      else
         call dense_copy(dbuffer_dm,full_tr)
         call dense_product(full_tr,dbuffer_dm,tr_buffer)
      end if

      ! Update DM and S
      call dense_convert(dbuffer_dm,init_dm)
      call dense_convert(dbuffer_ovlp,init_ovlp)

      call internal_dense_transform(dbuffer_dm,full_tr,tr_buffer)
      call internal_dense_transform(dbuffer_ovlp,full_tr,tr_buffer)

      call comms_barrier

    end subroutine internal_update_tr

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

 !--------------------------------------------------------------------------!
 ! Major subroutines                                                        !
 ! Contains transformation matrix generation procedures                     !
 !   1. internal_aopnao_tr     - AO to PNAO transformation with             !
 !                               lm-averaging                               !
 !   2. internal_norm_tr       - AO orbital normalization transformation    !
 !   3. internal_lclowdin_tr   - Atom-localised Lowdin transformation       !
 !   4. internal_wwsw_tr       - Weighted orthogoanlization (output in DEM) !
 !   5. internal_schmidt_tr    - Schmidt orthogonaliastion between sets     !
 !--------------------------------------------------------------------------!

    subroutine internal_aopnao_tr(aopnao_tr,pnao_w, &
         spinless_denskern,overlap,ngwf_lm_label,is_nmb,tr_type)
      ! lpl: Performs lm-averaging of density and overlap matrices to
      !      transform Alm blocks into Al blocks before diagonalising
      !      to obtain symmetry-averaged spinless AO to PNAO transformation.
      !      Returns both AO to PNAO transformation and PNAO diagonal weights
      !      obtained from diagonalising the Al blocks sorted from highest to
      !      lowest occupancies. NGWFs labelled with the same l will have the
      !      same transformation and diagonal PNAO weights
      ! lpl: 19082011 - Extension to NRB under construction
      use comms, only: pub_my_node_id
      use parallel_strategy, only: pub_first_atom_on_node, &
           pub_num_atoms_on_node
      use sparse, only: sparse_get_block, sparse_put_block
      use wrappers, only: wrappers_dsygv_lt

      implicit none

      ! lpl: Outputs
      ! AO to PNAO transformation matrix
      type(SPAM3), intent(inout)   :: aopnao_tr
      ! PNAO diagonal weights
      real(kind=DP), intent(inout) :: pnao_w(:)

      ! lpl: Inputs
      type(SPAM3), intent(in) :: spinless_denskern ! Spinless AO density matrix
      type(SPAM3), intent(in) :: overlap           ! AO overlap matrix

      integer, intent(in) :: ngwf_lm_label(:)      ! NGWF lm-symmetry labels
      logical, intent(in) :: is_nmb(:,:)           ! NMB/NRB flag
      character(len=*), intent(in) :: tr_type      ! Full or NRB only

      ! lpl: Local variables
      integer :: it, lc_it, iat, lc_iat
      integer :: irow, icol, nl_irow, nl_icol
      integer :: num_ngwfs, max_ngwfs

      integer :: l_max = 4
      integer :: l, n, m
      integer :: n_check, current_label

      ! lpl: Local arrays and matrices
      integer, allocatable :: lnm_list(:,:,:)
      integer, allocatable :: lm_n_num(:)

      real(kind=DP), allocatable :: lc_aopnao_tr(:,:)

      real(kind=DP), allocatable :: lc_dm(:,:)
      real(kind=DP), allocatable :: lc_ovlp(:,:)

      real(kind=DP), allocatable :: lc_buffer1(:,:)
      real(kind=DP), allocatable :: lc_buffer2(:,:)

      real(kind=DP), allocatable :: lc_eigval(:)

      logical, allocatable :: selected_orb(:,:)

      max_ngwfs = ngwf_basis%max_on_atom

      ! lpl: selected_orb contains index of orbital selected for lm-averaging
      allocate(selected_orb(pub_num_atoms_on_node(pub_my_node_id), &
           max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_aopnao_tr','selected_orb',ierr)

      ! lpl: For some reason having a seperate array for the selection
      !      criterion seems less complicated
      selected_orb = .false.
      select case(tr_type)
         case ('Full') ! lpl: Select all orbitals for aopnao_tr
            selected_orb = .true.
         case ('NRB')  ! lpl: Select only NRBs
            do lc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
               iat = lc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
               selected_orb(lc_iat,:) = .not.is_nmb(iat,:)
            end do
         case default
            write(stdout,'(a)') 'ERROR: Invalid internal_aopnao_tr tr_type'
            call comms_abort
      end select

      ! lpl: Orbital labelling lists
      ! lpl: (l,nm) mapping
      allocate(lnm_list(0:l_max,max_ngwfs,2*(l_max-1)+1),stat=ierr)
      call utils_alloc_check('internal_aopnao_tr','lnm_list',ierr)

      ! # n in each l
      allocate(lm_n_num(0:l_max),stat=ierr)
      call utils_alloc_check('internal_aopnao_tr','lm_n_num',ierr)

      ! lpl: Atom-local AO to PNAO transformation block
      allocate(lc_aopnao_tr(max_ngwfs,max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_aopnao_tr','lc_aopnao_tr',ierr)

      ! lpl: Atom-local symmetry-averaged density matrix & overlap matrices
      allocate(lc_dm(max_ngwfs,max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_aopnao_tr','lc_dm',ierr)
      allocate(lc_ovlp(max_ngwfs,max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_aopnao_tr','lc_ovlp',ierr)

      ! lpl: Atom-local buffer matrices
      allocate(lc_buffer1(max_ngwfs,max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_aopnao_tr','lc_buffer1',ierr)
      allocate(lc_buffer2(max_ngwfs,max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_aopnao_tr','lc_buffer2',ierr)

      ! lpl: Eigenvalue buffer
      allocate(lc_eigval(max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_aopnao_tr','lc_eigval',ierr)

      ! lpl: Iterate over all atoms on this node
      do lc_iat=1,pub_num_atoms_on_node(pub_my_node_id)

         lc_aopnao_tr = 0.0_DP

         iat = lc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
         num_ngwfs = ngwf_basis%num_on_atom(iat)

         lnm_list   = -1 ! Set all mapping indices to -1 for error-checking
         lm_n_num   =  0

         ! lpl: 19082011 - Initialisse ao_pnao_tr as identity so that elements
         !      that are not overwritten (e.g. because this is a NRB-only
         !      diagonalisation that does not span the whole atomic space)
         !      does nothing in terms of transformation
         do icol=1,ngwf_basis%num_on_atom(iat)
            lc_aopnao_tr(icol,icol) = 1.0_DP
         end do

         ! Get (Alm) to (Al) block mapping
         ! Iterate over all possible l channels (currently 0-4 (s to g))
         do l=0,l_max
            n = 0
            do m=1,2*l+1 ! Iterate over all 2l+1 m

               ! n_check to ensure equal of n for each m in current l
               n_check = n
               ! n = number of n for each l
               n = 0

               ! lpl: Get GENNBO labels
               current_label = l*100 + 50 + m
               ! lpl: If s then label is '1'
               if(current_label == 51) current_label = 1

               ! lpl: Iterate over all NGWFs on this atom
               do lc_it=1,ngwf_basis%num_on_atom(iat)

                  ! lpl: 19082011 - if this is an NRB-only diagonaliation
                  !      ignore all NMBs i.e. don't add to n or lnm_list
                  !      Decided by selected_orb = .true.
                  if(selected_orb(lc_iat,lc_it)) then
                     ! Get global NGWF index
                     it = ngwf_basis%first_on_atom(iat) + lc_it - 1
                     if(current_label == ngwf_lm_label(it)) then
                        n = n + 1
                        ! lnm_list stores local NGWF index mapping
                        lnm_list(l,n,m) = lc_it
                     end if
                  end if ! END if(selected_orb(lc_it))

               end do ! END do lc_it=1,ngwf_basis%num_on_atom(iat)

               ! lpl: Check that each m has the same number of nl
               if(m /= 1) then
                  if(n_check /= n) then
                     write(stdout,'(a,/a)') &
                          'ERROR: Inconsistent nl for different m in', &
                               'ao_pnao_tr lm-averaging'
                     call comms_abort
                  end if
               end if ! END if(m /= 1)

            end do ! END do m=1,2*l+1

            lm_n_num(l) =  n ! Number of n in this l

         end do ! END do l=0,l_max

         ! lpl: Perform lm-averaging and diagonalize Al block to get
         !      AO to PNAO transformation
         lc_dm   = 0.0_DP
         lc_ovlp = 0.0_DP
         call sparse_get_block(lc_dm(1:num_ngwfs,1:num_ngwfs), &
              spinless_denskern,iat,iat)
         call sparse_get_block(lc_ovlp(1:num_ngwfs,1:num_ngwfs), &
              overlap,iat,iat)

         do l=0,l_max

            lc_buffer1 = 0.0_DP
            lc_buffer2 = 0.0_DP

            n = lm_n_num(l) ! lpl: n = number of orbitals in this l
            ! lpl: Average over 2l+1 channels if orbitals in this l exist
            if(n > 0) then

               ! lpl: Symmetry-average Alm componenets into Al blocks
               do nl_irow=1,n
                  do nl_icol=1,n
                     do m=1,2*l+1
                        irow = abs(lnm_list(l,nl_irow,m))
                        icol = abs(lnm_list(l,nl_icol,m))

                        ! lpl: Sum density matrix
                        lc_buffer1(nl_irow,nl_icol) &
                             = lc_buffer1(nl_irow,nl_icol) + lc_dm(irow,icol)
                        ! lpl: Sum overlap
                        lc_buffer2(nl_irow,nl_icol) &
                             = lc_buffer2(nl_irow,nl_icol) + lc_ovlp(irow,icol)
                     end do
                  end do
               end do

               ! lpl: Average density matrix
               lc_buffer1 = lc_buffer1/(1.0_DP*(2*l+1))
               ! lpl: Average overlap matrix
               lc_buffer2 = lc_buffer2/(1.0_DP*(2*l+1))

               ! lpl: Diagonalize symmetry-avergaed Al blocks
               call wrappers_dsygv_lt(lc_buffer1(1:n,1:n),lc_eigval(1:n), &
                    lc_buffer2(1:n,1:n),n)

               ! lpl: Sort eigenvectors from highest to lowest occupancies
               !      and store into lc_buffer2 overriding averaged overlap
               lc_buffer2 = 0.0_DP
               do nl_icol=1,n
                  lc_buffer2(1:n,nl_icol) = lc_buffer1(1:n,n-nl_icol+1)
               end do

               ! lpl: Sort final transformation back in original NGWF order
               do nl_irow=1,n     ! Loop over dimension
                  do nl_icol=1,n  ! of this 'reduced' Al block

                     ! lpl: Loop over m to expand block back to original Alm
                     !      dimension
                     do m=1,2*l+1
                        irow = lnm_list(l,nl_irow,m)
                        icol = lnm_list(l,nl_icol,m)
                        lc_aopnao_tr(irow,icol) = lc_buffer2(nl_irow,nl_icol)

                        ! PNAO weights in highest to lowest order
                        it = irow + ngwf_basis%first_on_atom(iat) - 1
                        pnao_w(it) = lc_eigval(n-nl_irow+1)
                     end do

                  end do ! END do nl_irow=1,n
               end do ! END do nl_icol=1,n

            end if ! END if(n > 0)
         end do ! END do l=0,l_max

         call sparse_put_block(lc_aopnao_tr(1:num_ngwfs,1:num_ngwfs), &
              aopnao_tr,iat,iat)

      end do ! END do lc_iat=1,pub_num_atoms_on_node(pub_my_node_id)

      call comms_barrier

      ! lpl: Discard all elements not local to this node to avoid
      !      recombination of pnao_w that was redistributed from
      !      any previous routines
      do it=1,ngwf_basis%num
         if(ngwf_basis%node_of_func(it) /= pub_my_node_id) &
              pnao_w(it) = 0.0_DP
      end do
      call comms_reduce('SUM',pnao_w)

      ! Deallocate all matrices
      deallocate(lc_eigval,stat=ierr)
      call utils_dealloc_check('internal_aopnao_tr','lc_eigval',ierr)

      deallocate(lc_buffer2,stat=ierr)
      call utils_dealloc_check('internal_aopnao_tr','lc_buffer2',ierr)
      deallocate(lc_buffer1,stat=ierr)
      call utils_dealloc_check('internal_aopnao_tr','lc_buffer1',ierr)
      deallocate(lc_ovlp,stat=ierr)
      call utils_dealloc_check('internal_aopnao_tr','lc_ovlp',ierr)
      deallocate(lc_dm,stat=ierr)
      call utils_dealloc_check('internal_aopnao_tr','lc_dm',ierr)

      deallocate(lc_aopnao_tr,stat=ierr)
      call utils_dealloc_check('internal_aopnao_tr','lc_aopnao_tr',ierr)

      deallocate(lm_n_num,stat=ierr)
      call utils_dealloc_check('internal_aopnao_tr','lm_n_num',ierr)
      deallocate(lnm_list,stat=ierr)
      call utils_dealloc_check('internal_aopnao_tr','lnm_list',ierr)

      deallocate(selected_orb,stat=ierr)
      call utils_dealloc_check('internal_aopnao_tr','selected_orb',ierr)

      call comms_barrier

    end subroutine internal_aopnao_tr

    subroutine internal_norm_tr(norm_tr,overlap)
      ! lpl: Trivial normalization transformation
      use comms, only: pub_my_node_id
      use sparse, only: sparse_get_element, sparse_put_element
      implicit none

      type(SPAM3), intent(inout) :: norm_tr
      type(SPAM3), intent(in)    :: overlap

      integer :: it, lc_it
      real(kind=DP) mtx_el

      do lc_it=1,ngwf_basis%num_on_node(pub_my_node_id)
         it = lc_it + ngwf_basis%first_on_node(pub_my_node_id) - 1

         call sparse_get_element(mtx_el,overlap,it,it)
         mtx_el = 1.0_DP/sqrt(mtx_el)
         call sparse_put_element(mtx_el,norm_tr,it,it)
      end do

      call comms_barrier

    end subroutine internal_norm_tr

    subroutine internal_lclowdin_tr(lowdin_tr,overlap)
      ! lpl: Lowdin-orthogonalization transformation between orbitals
      !      sharing the same centre
      use comms, only: pub_my_node_id
      use parallel_strategy, only: pub_first_atom_on_node, &
           pub_num_atoms_on_node
      use sparse, only: sparse_get_block, sparse_put_block
      use wrappers, only: wrappers_dsyev_lt

      implicit none

      type(SPAM3), intent(inout) :: lowdin_tr  ! Output
      type(SPAM3), intent(in)    :: overlap    ! Input

      ! Local arrays
      real(kind=DP), allocatable :: lc_buffer(:,:)
      real(kind=DP), allocatable :: lc_eigvec(:,:)
      real(kind=DP), allocatable :: lc_eigval(:)
      real(kind=DP), allocatable :: lc_lowdin(:,:)
      real(kind=DP), allocatable :: lc_diag(:,:)

      ! Local variables
      integer :: num_ngwfs, max_ngwfs
      integer :: it, iat, lc_iat

      ! lpl: Allocate arrays
      max_ngwfs=ngwf_basis%max_on_atom

      allocate(lc_eigvec(max_ngwfs,max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_lclowdin_tr','lc_eigvec',ierr)
      allocate(lc_buffer(max_ngwfs,max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_lclowdin_tr','lc_buffer',ierr)
      allocate(lc_eigval(max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_lclowdin_tr','lc_eigval',ierr)
      allocate(lc_lowdin(max_ngwfs,max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_lclowdin_tr','lc_lowdin',ierr)
      allocate(lc_diag(max_ngwfs,max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_lclowdin_tr','lc_diag',ierr)

      ! lpl: Obtain atom-local Lowdin transformation matrix
      do lc_iat=1,pub_num_atoms_on_node(pub_my_node_id)

         iat = lc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
         num_ngwfs = ngwf_basis%num_on_atom(iat)

         lc_eigval = 0.0_DP
         lc_eigvec = 0.0_DP
         lc_diag   = 0.0_DP
         lc_buffer = 0.0_DP
         lc_lowdin = 0.0_DP

         ! lpl: Obtain atom-local S^(-1/2)
         if(num_ngwfs > 1) then
            call sparse_get_block(lc_eigvec(1:num_ngwfs,1:num_ngwfs), &
                 overlap,iat,iat)
            call wrappers_dsyev_lt(lc_eigvec(1:num_ngwfs,1:num_ngwfs), &
                 lc_eigval(1:num_ngwfs),num_ngwfs)

            do it=1,num_ngwfs ! lpl: D --> D^(-1/2)
               lc_diag(it,it) = 1.0_DP/sqrt(lc_eigval(it))
            end do

            ! V*D
            lc_buffer(1:num_ngwfs,1:num_ngwfs) = &
                 matmul( lc_eigvec(1:num_ngwfs,1:num_ngwfs), &
                      lc_diag(1:num_ngwfs,1:num_ngwfs) )
            ! (V*D)*V^(dagger)
            lc_lowdin(1:num_ngwfs,1:num_ngwfs) = &
                 matmul( lc_buffer(1:num_ngwfs,1:num_ngwfs), &
                      transpose(lc_eigvec(1:num_ngwfs,1:num_ngwfs)) )

            call sparse_put_block(lc_lowdin(1:num_ngwfs,1:num_ngwfs), &
                 lowdin_tr,iat,iat)
         else
            call sparse_get_block(lc_eigvec(1:num_ngwfs,1:num_ngwfs), &
                 overlap,iat,iat)
            lc_lowdin(1,1)=1/sqrt(lc_eigvec(1,1))
            call sparse_put_block(lc_lowdin(1:num_ngwfs,1:num_ngwfs), &
                 lowdin_tr,iat,iat)
         end if
      end do

      ! lpl: Deallocate arrays
      deallocate(lc_diag,stat=ierr)
      call utils_dealloc_check('internal_lclowdin_tr','lc_diag',ierr)
      deallocate(lc_lowdin,stat=ierr)
      call utils_dealloc_check('internal_lclowdin_tr','lc_lowdin',ierr)
      deallocate(lc_buffer,stat=ierr)
      call utils_dealloc_check('internal_lclowdin_tr','lc_buffer',ierr)
      deallocate(lc_eigval,stat=ierr)
      call utils_dealloc_check('internal_lclowdin_tr','lc_eigval',ierr)
      deallocate(lc_eigvec,stat=ierr)
      call utils_dealloc_check('internal_lclowdin_tr','lc_eigvec',ierr)

      call comms_barrier

    end subroutine internal_lclowdin_tr

    subroutine internal_diag_tr(diag_tr,pnao_w, &  ! Outputs
         spinless_denskern,overlap,is_nmb,tr_type) ! Inputs
      ! lpl: Diagonalisation transformation w/o lm-averaging
      use comms, only: pub_my_node_id
      use parallel_strategy, only: pub_first_atom_on_node, &
           pub_num_atoms_on_node
      use sparse, only: sparse_get_element, sparse_put_block, &
           sparse_put_element
      use wrappers, only: wrappers_dsygv_lt

      implicit none

      ! lpl: Outputs
      type(SPAM3), intent(inout)   :: diag_tr ! AO to PNAO transformation
      real(kind=DP), intent(inout) :: pnao_w(:) ! PNAO diagonal weights

      ! lpl: Inputs
      type(SPAM3), intent(in) :: spinless_denskern ! Spinless AO DM
      type(SPAM3), intent(in) :: overlap           ! AO overlap matrix

      logical, intent(in) :: is_nmb(:,:)           ! NMB/NRB flag
      character(len=*), intent(in) :: tr_type      ! Full or NRB only

      ! Local arrays
      real(kind=DP), allocatable :: lc_buffer(:,:)
      real(kind=DP), allocatable :: lc_eigvec(:,:)
      real(kind=DP), allocatable :: lc_eigval(:)

      logical, allocatable :: selected_orb(:,:)

      ! Local variables
      integer :: num_ngwfs, max_ngwfs, n_orbs
      integer, allocatable :: n_list(:)
      integer :: it, iat, lc_iat, lc_irow, lc_icol

      ! lpl: Indicate in stdout that non-standard aopnao_tr is being used
      if(pub_on_root) write(stdout,'(a)') 'WARNING: Non-standard AO-PNAO &
         &transformation (DIAG) is being used'

      ! lpl: Allocate arrays
      max_ngwfs=ngwf_basis%max_on_atom

      ! lpl: selected_orb contains index of orbital selected for lm-averaging
      allocate(selected_orb(pub_num_atoms_on_node(pub_my_node_id), &
           max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_diag_tr','selected_orb',ierr)

      allocate(n_list(max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_diag_tr','n_list',ierr)

      ! lpl: For some reason having a seperate array for the selection
      !      criterion seems less complicated
      selected_orb = .false.
      select case(tr_type)
         case ('Full') ! lpl: Select all orbitals for aopnao_tr
            selected_orb = .true.
         case ('NRB')  ! lpl: Select only NRBs
            do lc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
               iat = lc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
               selected_orb(lc_iat,:) = .not.is_nmb(iat,:)
            end do
         case default
            write(stdout,'(a)') 'ERROR: Invalid internal_aopnao_tr tr_type'
            call comms_abort
      end select

      allocate(lc_eigvec(max_ngwfs,max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_diag_tr','lc_eigvec',ierr)
      allocate(lc_buffer(max_ngwfs,max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_diag_tr','lc_buffer',ierr)
      allocate(lc_eigval(max_ngwfs),stat=ierr)
      call utils_alloc_check('internal_diag_tr','lc_eigval',ierr)

      ! lpl: Obtain atom-local transformation matrix
      do lc_iat=1,pub_num_atoms_on_node(pub_my_node_id)

         iat = lc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
         num_ngwfs = ngwf_basis%num_on_atom(iat)

         ! lpl: Obtain list of selected orbitals
         n_orbs = 0
         n_list = -1 ! lpl: Set to -1 for error-catching
         do lc_it=1,num_ngwfs
            if(selected_orb(lc_iat,lc_it)) then
               n_orbs = n_orbs + 1
               n_list(n_orbs) = lc_it + ngwf_basis%first_on_atom(iat) - 1
            end if
         end do

         ! lpl: Explicitly place identity into initial transformation
         lc_buffer = 0.0_DP
         do lc_icol=1,num_ngwfs
            lc_buffer(lc_icol,lc_icol) = 1.0_DP
         end do
         call sparse_put_block(lc_buffer(1:num_ngwfs,1:num_ngwfs), &
              diag_tr,iat,iat)

         lc_eigval = 0.0_DP
         lc_eigvec = 0.0_DP
         lc_buffer = 0.0_DP

         ! lpl: Obtain atom-local S^(-1/2)
         if(n_orbs > 1) then

            do lc_irow=1,n_orbs
               do lc_icol=1,n_orbs
                  call sparse_get_element(lc_eigvec(lc_irow,lc_icol), &
                       spinless_denskern,n_list(lc_irow),n_list(lc_icol))
                  call sparse_get_element(lc_buffer(lc_irow,lc_icol), &
                       overlap,n_list(lc_irow),n_list(lc_icol))
               end do
            end do

            call wrappers_dsygv_lt(lc_eigvec(1:n_orbs,1:n_orbs), &
                                   lc_eigval(1:n_orbs), &
                                   lc_buffer(1:n_orbs,1:n_orbs), &
                                   n_orbs)

         else if(n_orbs == 1) then ! Simply normalize if only 1 NGWF on atom

            call sparse_get_element(lc_eigvec(1,1),spinless_denskern, &
                 n_list(1),n_list(1))
            call sparse_get_element(lc_buffer(1,1),overlap, &
                 n_list(1),n_list(1))
            ! W = P/S
            lc_eigval(1) = 1.0_DP*lc_eigvec(1,1)/lc_buffer(1,1)
            ! M = S^(-1/2)
            lc_eigvec(1,1)=1/sqrt(lc_buffer(1,1))

         end if

         ! lpl: Sort eigenvectors & eigenvalues from highest to lowest
         !      occupancy - correct transformation now in lc_buffer
         do lc_it=1,n_orbs
            lc_buffer(1:n_orbs,lc_it) &
               = lc_eigvec(1:n_orbs,n_orbs-lc_it+1)

            ! lpl: Replace relevant occupancies in global pnao_w
            pnao_w(n_list(lc_it)) = lc_eigval(n_orbs-lc_it+1)
         end do

         ! lpl: Replace relevant elements in diag_tr if required
         do lc_irow=1,n_orbs
            do lc_icol=1,n_orbs
               call sparse_put_element(lc_buffer(lc_irow,lc_icol), &
                    diag_tr,n_list(lc_irow),n_list(lc_icol))
            end do
         end do

      end do ! END do lc_iat=1,pub_num_atoms_on_node(pub_my_node_id)

      call comms_barrier

      ! lpl: Discard all elements not local to this node to avoid
      !      recombination of pnao_w that was redistributed from
      !      any previous routines
      do it=1,ngwf_basis%num
         if(ngwf_basis%node_of_func(it) /= pub_my_node_id) &
              pnao_w(it) = 0.0_DP
      end do
      call comms_reduce('SUM',pnao_w)

      ! lpl: Deallocate arrays
      deallocate(lc_buffer,stat=ierr)
      call utils_dealloc_check('internal_diag_tr','lc_buffer',ierr)
      deallocate(lc_eigval,stat=ierr)
      call utils_dealloc_check('internal_diag_tr','lc_eigval',ierr)
      deallocate(lc_eigvec,stat=ierr)
      call utils_dealloc_check('internal_diag_tr','lc_eigvec',ierr)

      deallocate(n_list,stat=ierr)
      call utils_dealloc_check('internal_diag_tr','n_list',ierr)
      deallocate(selected_orb,stat=ierr)
      call utils_dealloc_check('internal_diag_tr','selected_orb',ierr)

      call comms_barrier

    end subroutine internal_diag_tr

    subroutine internal_schmidt_tr(schmidt_tr,ovlp,list,num, &
         orth_list,orth_num,col_buffer,norm)
      ! lpl: Schmidt orthogonalisation between orbitals in list(num) w.r.t.
      !      all other orbitals, transforming only orbitals in list(num).
      !      Orbitals outside list assumed to be orthogonal
      use dense, only: dense_get_col, dense_put_col, dense_put_element
      implicit none

      type(DEM), intent(inout) :: schmidt_tr ! Transformation
      type(DEM), intent(in) :: ovlp ! Input DEM overlap
      real(kind=DP), intent(inout) :: col_buffer(:)

      integer, intent(in) :: num, orth_num
      integer, intent(in) :: list(num), orth_list(orth_num)

      ! Normalization flag (assumes ortho outside list)
      logical, optional, intent(in) :: norm

      integer :: irow, icol, list_irow, list_icol
      logical :: normalize
      real(kind=DP) :: norm_fac, mtx_el

      ! lpl: Check that this matrix is appropriate
      if(schmidt_tr%nrows /= schmidt_tr%mcols) then
         if(pub_on_root) write(stdout,'(a)') 'ERROR: nrows /= mcols in &
              &schmidt_tr in internal_schmidt_tr in npa_main'
         call comms_abort
      else if(schmidt_tr%nrows /= ngwf_basis%num) then
         if(pub_on_root) write(stdout,'(a)') 'ERROR: nrows /= ngwf_basis%num&
              & in schmidt_tr in internal_schmidt_tr in npa_main'
         call comms_abort
      end if

      normalize = .false. ! Default behaviour
      if(present(norm)) normalize = norm

      ! lpl: Initialise as identity, then replace with appropriate elements
      call internal_dense_id(schmidt_tr,ngwf_basis%num)

      if(num > 0) then
         ! lpl: Loop over all overlap of orbitals in list(num) only
         do list_icol=1,num

            col_buffer = 0.0_DP
            norm_fac = 1.0_DP

            icol = list(list_icol) ! Get correct column index for this orbital

            ! Subtract (orthogonalize) S_ij w.r.t. orth_list
            do list_irow=1,orth_num
               irow = orth_list(list_irow)
               if(irow == icol) then ! Sanity check
                  write(stdout,'(a)') 'ERROR: orth_list(:) contains &
                       &orbital in list(:) in internal_schmidt_tr in npa_main'
                  call comms_abort
               end if
               call dense_get_element(mtx_el,ovlp,irow,icol)

               col_buffer(irow) = -1.0_DP*mtx_el
               norm_fac = norm_fac - mtx_el**2
            end do

            col_buffer(icol) = 1.0_DP ! Do nothing with itself

            if(.not. normalize) norm_fac = 1.0_DP
            call dense_put_col(col_buffer/sqrt(norm_fac),schmidt_tr,icol)
         end do
      end if

      call comms_barrier

    end subroutine internal_schmidt_tr

    subroutine internal_wwsw_tr(wwsw_tr,ovlp_buff, &
         ao_weights,ao_list,num,buff1,col_buffer)
      ! lpl: Weighted orthogonalization transformation W(WSW)^(-1/2)
      !      Transformation matrix is of type(DEM)
      !      ovlp_buff contains input overlap mat and doubles as buffer
      !      wwsw_tr start as buffer and ends up as wwsw_tr
      ! lpl: 23092011 - pnao_w mat weights (W) < scaled to [0,1], and
      !      all W_i < w_min set to w_min. This avoids possible errors
      !      in 1/WSW esp. when W_i ~0, leaving such problems to only
      !      manifest in linear-dependence of the AO basis
      implicit none

      ! lpl: Input
      type(DEM), intent(inout)   :: ovlp_buff
      real(kind=DP), intent(in) :: ao_weights(:)

      integer, intent(in)       :: num
      integer, intent(in)       :: ao_list(num)

      ! lpl: Output
      type(DEM), intent(inout) :: wwsw_tr

      ! lpl: Extra buffers from main routine
      type(DEM), intent(inout) :: buff1
      real(kind=DP), intent(inout) :: col_buffer(:)

      integer :: it, non_wwsw_num
      real(kind=DP) :: w_max

      if(num > 0) then
         ! lpl: Check that sensible AO list is being given
         non_wwsw_num = ngwf_basis%num - num
         if(non_wwsw_num < 0) then
            write(stdout,'(a)') 'ERROR: # orbitals excluded from &
                 &internal_wwsw_tr < 0'
            call comms_abort
         end if

         ! lpl: Read W diagonal into buff1
         col_buffer = 0.0_DP
         do it=1,num
            col_buffer(it) = ao_weights(ao_list(it))
         end do
         ! lpl: 23092011 - Scale max weight in W to 1 and set all
         !      W_i < w_min to w_min
         w_max = maxval(col_buffer(1:num))
         col_buffer = col_buffer/w_max
         do it=1,num
            if(col_buffer(it) < w_min) col_buffer(it) = w_min
         end do
         call internal_dense_coldiag(buff1,col_buffer(1:num),num)

         ! lpl: Read S into local block
         call internal_dense_reblock(wwsw_tr,ovlp_buff,ao_list,num)

         ! lpl: Get WSW into wwsw_tr
         !      olvp_buff now used as a regular buffer
         call dense_product(ovlp_buff,wwsw_tr,buff1, &
              transpose_amat=.false.,transpose_bmat=.false., &
              first_k=1,last_k=num) ! SW
         call dense_product(wwsw_tr,buff1,ovlp_buff, &
              transpose_amat=.false.,transpose_bmat=.false., &
              first_k=1,last_k=num) ! W(SW)

#ifdef DEBUG
   call internal_testprint_dem(debug_output_unit,buff1, &
            'W MAT')
   call internal_testprint_dem(debug_output_unit,wwsw_tr, &
            '(WSW) MAT')
#endif

         ! lpl: Solve eigenvalue problem for WSW to find (WSW)^(-1/2)
         col_buffer = 0.0_DP
         ! lpl: wwsw_tr now contains eigenvectors
         call internal_dense_block_dsyev(wwsw_tr,col_buffer,1,num)

#ifdef DEBUG
   if(pub_on_root) then
      write(stdout,'(a)') 'WSW eigenvalues'
      do it=1,ngwf_basis%num
         write(debug_output_unit,'(f16.10)') col_buffer(it)
      end do
   end if
#endif

         ! lpl: Calculate SQRT of eigenvalues diagonal to obtain (WSW)^(-1/2)
         do it=1,num
            if(col_buffer(it) < s_threshold) write(stdout,'(a,/a)') &
                 'WARNING: WSW eigenvalue < numerical threshold', &
                 ' - Final NAOs might not be accurate'
            col_buffer(it) = 1.0_DP/sqrt(col_buffer(it))
         end do
         call comms_barrier
         ! lpl: buff1 contains eigenvalue diagonals
         call internal_dense_coldiag(buff1,col_buffer(1:num),num)

#ifdef DEBUG
   call internal_testprint_dem(debug_output_unit,wwsw_tr, &
            'NMB WSW EIGVECS')
#endif

         ! lpl: Get (WSW)^-(1/2) (buff1)
         ! D*V^(dagger)
         call dense_product(ovlp_buff,buff1,wwsw_tr, &
              transpose_amat=.false.,transpose_bmat=.true., &
              first_k=1,last_k=num)
         ! V*[D*V^(dagger)] = (WSW)^(-1/2)
         call dense_product(buff1,wwsw_tr,ovlp_buff, &
              transpose_amat=.false.,transpose_bmat=.false., &
              first_k=1,last_k=num)
         call comms_barrier

         ! lpl: Read diagonal AO occupancies into W again (wwsw_tr)
         col_buffer = 0.0_DP
         do it=1,num
            col_buffer(it) = ao_weights(ao_list(it))
         end do
         ! lpl: 23092011 - Scale max weight in W to 1 and set all
         !      W_i < w_min to w_min
         col_buffer = col_buffer/w_max
         do it=1,num
            if(col_buffer(it) < w_min) col_buffer(it) = w_min
         end do
         call internal_dense_coldiag(wwsw_tr,col_buffer(1:num),num)

         ! lpl: W(WSW)^(-1/2)
         call dense_product(ovlp_buff,wwsw_tr,buff1, &
              transpose_amat=.false.,transpose_bmat=.false., &
              first_k=1,last_k=num)

         ! lpl: Nominal identity transform for orbitals not in AO list
         call internal_dense_id(wwsw_tr,ngwf_basis%num)
         ! lpl: Replace transformation for orbitals in AO list with wwsw_tr
         call internal_dense_deblock(ovlp_buff,wwsw_tr,ao_list,num)

      else
         call internal_dense_id(wwsw_tr,ngwf_basis%num)
      end if

      call comms_barrier

    end subroutine internal_wwsw_tr

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

 !--------------------------------------------------------------------------!
 ! Dense matrices block subroutines                                         !
 ! Extracts and inserts dense blocks based on user-specified list of        !
 ! indices                                                                  !
 !    1. internal_dense_get_block                                           !
 !    2. internal_dense_put_block                                           !
 !--------------------------------------------------------------------------!

    subroutine internal_dense_get_block(sub_mat,full_mat,list,num)
      ! lpl: Extracts num x num block from full square dense matrix based on
      !      set of indices given [list(:)]
      use dense, only: dense_scale, dense_put_element

      implicit none

      ! Output
      type(DEM), intent(inout) :: sub_mat ! Output dense sub-matrix

      ! Input
      type(DEM), intent(in) :: full_mat ! Full dense matrix input
      integer, intent(in)   :: list(:)  ! List of elements spanning sub-matrix
      integer, intent(in)   :: num      ! Dimension of sub-matrix

      ! Local variables
      integer :: irow, icol
      real(kind=DP) :: mtx_el

      ! lpl: Check that inputs are sensible
      if(full_mat%nrows /= full_mat%mcols) then
         write(stdout,'(a)') &
              'ERROR: Non-square full_mat in internal_dense_get_block'
         call comms_abort
      end if

      if(num < 0) then
         write(stdout,'(a)') &
              'ERROR: Sub-matrix size < 0 in internal_dense_get_block'
         call comms_abort
      elseif(num > full_mat%nrows) then
         write(stdout,'(a,/a)') 'ERROR: Sub-matrix size > full matrix &
              &dimension','in internal_dense_get_block'
         call comms_abort
      end if

      if(sub_mat%nrows /= sub_mat%mcols) then
         write(stdout,'(a)') &
              'ERROR: Non-square sub_mat in internal_dense_get_block'
         call comms_abort
      end if

      if(sub_mat%nrows /= num) then
         write(stdout,'(a,/a)') 'ERROR: Sub-mat size inconsistent with num', &
              'in internal_dense_get_block'
         call comms_abort
      end if

      call dense_scale(sub_mat,0.0_DP)
      call comms_barrier
      do irow=1,num
         do icol=1,num
            call dense_get_element(mtx_el,full_mat,list(irow),list(icol))
            call dense_put_element(mtx_el,sub_mat,irow,icol)
         end do
      end do
      call comms_barrier

    end subroutine internal_dense_get_block

    subroutine internal_dense_put_block(sub_mat,full_mat,list,num)
      ! lpl: Inserts num x num block from dense sub-matrix into a full dense
      !      matrix based on set of indices given [list(:)], overriding the
      !      full matrix elements
      use dense, only: dense_put_element, dense_scale

      implicit none

      ! Output
      type(DEM), intent(inout) :: full_mat ! Full dense matrix output

      ! Input
      type(DEM), intent(in) :: sub_mat  ! Input dense sub-matrix
      integer, intent(in)   :: list(:)  ! List of elements spanning sub-matrix
      integer, intent(in)   :: num      ! Dimension of sub-matrix

      ! Local variables
      integer :: irow, icol
      real(kind=DP) :: mtx_el

      ! lpl: Check that inputs are sensible
      if(full_mat%nrows /= full_mat%mcols) then
         write(stdout,'(a)') &
              'ERROR: Non-square full_mat in internal_dense_get_block'
         call comms_abort
      end if

      if(num < 0) then
         write(stdout,'(a)') &
              'ERROR: Sub-matrix size < 0 in internal_dense_get_block'
         call comms_abort
      elseif(num > full_mat%nrows) then
         write(stdout,'(a,/a)') 'ERROR: Sub-matrix size > full matrix &
              &dimension','in internal_dense_get_block'
         call comms_abort
      end if

      if(sub_mat%nrows /= sub_mat%mcols) then
         write(stdout,'(a)') &
              'ERROR: Non-square sub_mat in internal_dense_get_block'
         call comms_abort
      end if

      if(sub_mat%nrows /= num) then
         write(stdout,'(a,/a)') 'ERROR: Sub-mat size inconsistent with num', &
              'in internal_dense_get_block'
       call comms_abort
      end if

      do irow=1,num
         do icol=1,num
            call dense_get_element(mtx_el,sub_mat,irow,icol)
            call dense_put_element(mtx_el,full_mat,list(irow),list(icol))
         end do
      end do
      call comms_barrier

    end subroutine internal_dense_put_block

    subroutine internal_dense_reblock(buff_mat,full_mat,list,num)
      ! lpl: Extracts num x num block from full square dense matrix based on
      !      set of indices given [list(:)]
      use dense, only: dense_put_element

      implicit none

      ! Output
      type(DEM), intent(inout) :: buff_mat ! Reblocked mat

      ! Input
      type(DEM), intent(in) :: full_mat ! Full dense matrix input
      integer, intent(in)   :: list(:)  ! List of elements spanning sub-matrix
      integer, intent(in)   :: num      ! Dimension of sub-matrix

      ! Local variables
      integer :: irow, icol
      real(kind=DP) :: mtx_el

      ! lpl: Check that inputs are sensible
      if(full_mat%nrows /= full_mat%mcols) then
         write(stdout,'(a)') &
              'ERROR: Non-square full_mat in internal_dense_get_block'
         call comms_abort
      end if

      if(num < 0) then
         write(stdout,'(a)') &
              'ERROR: Sub-matrix size < 0 in internal_dense_get_block'
         call comms_abort
      elseif(num > full_mat%nrows) then
         write(stdout,'(a,/a)') 'ERROR: Sub-matrix size > full matrix &
              &dimension','in internal_dense_get_block'
         call comms_abort
      end if

      if(buff_mat%nrows /= buff_mat%mcols) then
         write(stdout,'(a)') &
              'ERROR: Non-square sub_mat in internal_dense_get_block'
         call comms_abort
      end if

      do irow=1,buff_mat%nrows ! lpl: Explicitly zero mat
         do icol=1,buff_mat%mcols
            call dense_put_element(0.0_DP,buff_mat,irow,icol)
         end do
      end do
      call comms_barrier

      do irow=1,num ! lpl: Re-block full_mat using list into buff_mat
         do icol=1,num
            call dense_get_element(mtx_el,full_mat,list(irow),list(icol))
            call dense_put_element(mtx_el,buff_mat,irow,icol)
         end do
      end do
      call comms_barrier

    end subroutine internal_dense_reblock

    subroutine internal_dense_deblock(buff_mat,full_mat,list,num)
      ! lpl: Inserts num x num block from dense sub-matrix into a full dense
      !      matrix based on set of indices given [list(:)], overriding the
      !      targeted full matrix elements
      use dense, only: dense_put_element

      implicit none

      ! Output
      type(DEM), intent(inout) :: full_mat ! Full dense matrix output

      ! Input
      type(DEM), intent(in) :: buff_mat  ! Input dense blocked mat
      integer, intent(in)   :: list(:)  ! List of elements spanning sub-matrix
      integer, intent(in)   :: num      ! Dimension of sub-matrix

      ! Local variables
      integer :: irow, icol
      real(kind=DP) :: mtx_el

      ! lpl: Check that inputs are sensible
      if(full_mat%nrows /= full_mat%mcols) then
         write(stdout,'(a)') &
              'ERROR: Non-square full_mat in internal_dense_get_block'
         call comms_abort
      end if

      if(num < 0) then
         write(stdout,'(a)') &
              'ERROR: Sub-matrix size < 0 in internal_dense_get_block'
         call comms_abort
      elseif(num > full_mat%nrows) then
         write(stdout,'(a,/a)') 'ERROR: Sub-matrix size > full matrix &
              &dimension','in internal_dense_get_block'
         call comms_abort
      end if

      if(buff_mat%nrows /= buff_mat%mcols) then
         write(stdout,'(a)') &
              'ERROR: Non-square sub_mat in internal_dense_get_block'
         call comms_abort
      end if

      do irow=1,num
         do icol=1,num
            call dense_get_element(mtx_el,buff_mat,irow,icol)
            call dense_put_element(mtx_el,full_mat,list(irow),list(icol))
         end do
      end do
      call comms_barrier

    end subroutine internal_dense_deblock


 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

 !--------------------------------------------------------------------------!
 ! NGWF labelling subroutines                                               !
 ! Labels NGWFs with false lm-symmetry based on ONETEP's default            !
 ! initialization order                                                     !
 !    1. internal_generatelabel_manual                                      !
 !    2. internal_generatelabel_auto                                        !
 !    3. internal_ngwf_sto3g_create  } Copied from similarly-named internal !
 !    4. internal_ngwf_sto3g_destroy } subroutines in ngwfs_mod             !
 !--------------------------------------------------------------------------!

    subroutine internal_generatelabel_manual(species_number, &
         species_ngwf_label,species_is_nmb,ngwf_label_option)
      ! lpl: Label NGWFs with lm-labels and as NMB/NRBs based on user-
      !      specified list in input file
      implicit none

      ! lpl: Inputs/Outputs
      ! lpl: Label string - discarded at the end
      character(len=*), intent(inout) :: ngwf_label_option
      ! lpl: NGWF lm-labels and NMB logical flag for this species
      integer, intent(inout) :: species_ngwf_label(0:)
      logical, intent(inout) :: species_is_nmb(:)
      ! lpl: Current species number
      integer, intent(in) :: species_number

      ! lpl: Local variables
      integer :: ipos, jpos, end_pos, sub_key
      integer :: lm, num, sub_it
      character(len=8) :: label
      character :: occ_type
      logical   :: valid_lm

      ngwf_label_option = trim(adjustl(ngwf_label_option))
      end_pos = len_trim(ngwf_label_option)

      num  = 0  ! lpl: Local NGWF index of lm labels in this string
      jpos = 1
      do while(jpos < end_pos)
         ipos = jpos
         ! lpl: Skip blank spaces
         do while(ngwf_label_option(ipos:ipos) == ' ')
            ipos = ipos + 1
         end do

         jpos = index(ngwf_label_option(ipos:),' ')
         jpos = ipos + jpos - 1

         label   = trim(adjustl(ngwf_label_option(ipos:jpos)))
         sub_key = index(label,'N')
         if(sub_key == 0) then
            occ_type = 'R'
         else
            occ_type = label(sub_key:sub_key)
            label    = label(1:sub_key-1)
         end if

         occ_type = trim(adjustl(occ_type))
         label    = trim(adjustl(label))

         ! lpl: Get integer format of this label
         read(label,*) lm
         num = num + 1

         valid_lm = .false.
         ! lpl: Check that this label is valid
         do sub_it=1,size(nbo_lm_label)
            if(lm == nbo_lm_label(sub_it)) then
               valid_lm = .true.
               exit
            end if
         end do
         if(.not. valid_lm) then
            write(stdout,'(a,/a)') 'ERROR: Invalid user-specified lm label &
                 &detected','in internal_read_ngwf_labels'
            call comms_abort
         end if

         ! lpl: Put info into local arrays lc_ngwf_lm_label and lc_is_nmb
         species_ngwf_label(num) = lm
         if(occ_type == 'N') species_is_nmb(num) = .true.
      end do ! END do while(jpos < end_pos)

      ! lpl: If # NGWFs /= # labels in list
      if(num /= &
           ngwf_basis%num_on_atom(pub_distr_atom(species_ngwf_label(0)))) then
         write(stdout,'(a,/a)') 'ERROR: Inconsistent # NGWFs in &
              &user-specified lm label detected','in internal_read_ngwf_labels'
         call comms_abort
      end if

    end subroutine internal_generatelabel_manual

    subroutine internal_generatelabel_auto(species_number,species_ngwf_label, &
        species_is_nmb,gbasis)
      ! lpl: Generate lm-labels based on ONETEP's default STO-3G initial basis
      implicit none

      ! lpl: Outputs
      ! lpl: NGWF lm-labels and NMBs for this species
      integer, intent(inout) :: species_ngwf_label(0:)
      logical, intent(inout) :: species_is_nmb(:)
      ! lpl: Current species number
      integer, intent(in) :: species_number

      ! lpl: Input
      ! lpl: GTO set generated once from main subroutine
      type(GTO_SET), intent(in) :: gbasis(:)

      ! lpl: Local variables
      integer :: shell, start_shell, fbl_shells, l, m, sub_it, pseud_shells
      integer :: atnum, shell_index
      real(kind=DP) :: shell_charge, ion_charge

      ! lpl: Local arrays
      ! lpl: Number of g.s. shells
      integer :: num_gs_shells(109)
      ! lpl: Full sequence of NGWF lm-labels from 1:n, where n >= num_ngwfs
      integer :: species_label_full(0:max_ngwf_label)

      ! lpl: Number of g.s. shells - derived from gs_occ in atoms_mod
      !      Assumes correct energetic ordering of initial NGWFs
      !      (e.g. 3p,4s,3d)
      data num_gs_shells(:) /   1,1,2,2,5,5,5,5,5,5,6,6,9,9,9,9,9,9,10,10, &
           15,15,15,15,15,15,15,15,15,15,18,18,18,18,18,18,19,19,24,24,24, &
           24,24,24,24,23,24,24,27,27,27,27,27,27,28,28,33,40,35,35,35,35, &
           35,40,35,35,35,35,35,35,40,40,40,40,40,40,40,40,40,40,43,43,43, &
           43,43,43,44,44,49,49,56,56,56,51,51,56,51,51,51,51,51,51,56,56, &
           56,56,56,56,56 /

      ! lpl: E.g. of species
      orig_iat = species_ngwf_label(0)
      num_ngwfs = ngwf_basis%num_on_atom(pub_distr_atom(orig_iat))

      ! lpl: Generate NGWF lm-label for each species based on ONETEP's STO-3G
      !      initialization pattern taking into account pseudised core charges
      !      Copied & modified from ngwfs_generate_sto3g in ngwfs_mod. Needs
      !      gbasis(:). Outputs species_label_full (contains lm label from 1st
      !      unpseudized shell to max available from ONETEP ngwf_data_mod
      !
      !      Original comments:
      !      cks: find starting shell of STO-3G basis from which to start
      !           fireball set
      !      cks: (depends on ion charge of pseudopotential)

      ! lpl: Get effective Z and A for current species
      atnum      = elements(orig_iat)%atomic_number
      ion_charge = elements(orig_iat)%ion_charge

      start_shell  = -1
      shell_charge =  0
      find_start_shell_loop: do shell=1,gbasis(atnum)%nshells
         shell_charge = shell_charge + 2*(2*gbasis(atnum)%angmom(shell) + 1)
         if (shell_charge >= ( atnum-ion_charge ) ) then
            if (shell_charge == ( atnum-ion_charge ) ) start_shell = shell + 1
            if (shell_charge > ( atnum-ion_charge ) ) start_shell = shell
            exit find_start_shell_loop
         endif
      end do find_start_shell_loop

      ! lpl: 23092011 - Corrected bug that gives wrong #NMBs (pseud_shells)
      species_label_full = 0
      sub_it = 0
      pseud_shells = 0 ! lpl: Contains # pseudised NGWFs (nlm)
      do shell=1,gbasis(atnum)%nshells
         l = gbasis(atnum)%angmom(shell)
         ! cks: initialise angular momentum for current shell
         do m=-l,l
            if(shell >= start_shell) then ! lpl: If shell is not pseudised
               ! lpl: sub_it tallies num generated shells starting from
               !      start_shell
               sub_it = sub_it + 1
               if(sub_it > max_ngwf_label) then
                  if(pub_on_root) write(stdout,'(a,I5,a)') &
                       'Error: # sub-shells > max_ngwf_label (', &
                            max_ngwf_label,').'
                  call comms_abort
               end if
               ! lpl: Store NBO lm-label
               species_label_full(sub_it) = nbo_lm_label((l*l)+l+m+1)
            else
               pseud_shells = pseud_shells + 1
            end if ! END if(shell >= start_shell)
         end do ! END do m=-l,l
      end do

      ! lpl: Store number of generated shells into sspecies_label_full(0)
      species_label_full(0) = sub_it

#ifdef DEBUG
      write(debug_output_unit,'(1x,a)') '--------------------'
      write(debug_output_unit,'(1x,a,I4,1x,a)') 'Atom: ', &
           orig_iat,elements(orig_iat)%species_id
      write(debug_output_unit,'(1x,a)') '--------------------'
      write(debug_output_unit,'(1x,a)') 'NGWF#       lm-label'
      write(debug_output_unit,'(1x,a)') '--------------------'
      do sub_it=1,species_label_full(0)
         write(debug_output_unit,'(1x,I3,1x,I5)') &
              sub_it,species_label_full(sub_it)
      end do
      write(debug_output_unit,'(1x,a)') '--------------------'
# endif

      ! lpl: Label NGWFs on current species from 1:num_ngwfs using sequence
      !      from species_label_full. If num_ngwfs > max available from
      !      ngwf_data_mod excess orbitals are labelled as s

      ! lpl: Warn if # NGWFs on atom > max lm available for current species
      if(num_ngwfs >  species_label_full(0)) then
         if(pub_on_root) write(stdout,'(a,I5,a,/a)') &
              'Warning: # NGWFs on species ', &
                   species_number,' > max. orbital labels available', &
                        'Labelling excess orbitals as s.'
      end if

      ! lpl: Iterate over all NGWFs for current species
      do sub_it=1,num_ngwfs
         ! lpl: If # NGWFs > max supported label extra as
         if(sub_it > max_ngwf_label) then
            species_ngwf_label(sub_it) = nbo_lm_label(1)
         else
            species_ngwf_label(sub_it) = species_label_full(sub_it)
         end if
      end do

      ! lpl: Label 1st n unpseudized g.s. NGWFs as NMBs
      ! lpl: 23092011 - corrected bug (start_shell labels n, not nlm)
      !      so we were subtracting wrong # pseudised NGWFs
      shell = num_gs_shells(atnum) - pseud_shells

#ifdef DEBUG
      write(debug_output_unit,'(a)') '--------------------------------'
      write(debug_output_unit,'(a)') 'Species ID: ', &
           elements(orig_iat)%species_id
      write(debug_output_unit,'(a,I5)') 'Atomic number      : ', atnum
      write(debug_output_unit,'(a,I5)') '# full g.s. shells : ', &
           num_gs_shells(atnum)
      write(debug_output_unit,'(a,I5)') 'Start shell        : ', start_shell
      write(debug_output_unit,'(a,I5)') 'Psuedised # nlm    : ', pseud_shells
      write(debug_output_unit,'(a,I5)') '# Actual max # NMB : ', shell
      write(debug_output_unit,'(a,I5)') '# NGWFs            : ', num_ngwfs
      write(debug_output_unit,'(a)') '--------------------------------'
# endif

      ! lpl: Label 1st n unpseudized g.s. NGWFs as NMBs
      do sub_it=1,num_ngwfs
         if(sub_it <= shell) species_is_nmb(sub_it) = .true.
      end do

#ifdef DEBUG
      write(debug_output_unit,'(a)') '----------------------------'
      write(debug_output_unit,'(1x,a,I5,1x,a,a)') 'Atom ', &
           orig_iat,elements(orig_iat)%species_id,' NGWF labels:'
      write(debug_output_unit,'(a)') '----------------------------'
      write(debug_output_unit,'(a)') 'NGWF#  lm-label    is_nmb   '
      write(debug_output_unit,'(a)') '----------------------------'
      do sub_it=1,num_ngwfs
         if(species_is_nmb(sub_it)) then
            write(debug_output_unit,'(1x,I3,4x,I5,5x,a)') &
                 sub_it,species_ngwf_label(sub_it),'Y'
         else
            write(debug_output_unit,'(1x,I3,4x,I5,5x,a)') &
                 sub_it,species_ngwf_label(sub_it),'N'
         end if
      end do
      write(debug_output_unit,'(a)') '----------------------------'
# endif

    end subroutine internal_generatelabel_auto

    subroutine internal_ngwf_sto3g_create

      !===================================================================!
      ! Allocate appropriate memory and initialise it to hold Gaussian    !
      ! basis set parameters (all STO-3G and only the polarisation        !
      ! functions of 6-31G*) for each element up to atomic number         !
      ! pub_num_gtatoms                                                   !
      !-------------------------------------------------------------------!
      ! lpl: Copied from ngwfs_mod (23-May-2011)                          !
      !===================================================================!

      use constants, only: DP
      use ngwf_data, only: GTO_SET, pub_num_gtatoms, &
                           ngwf_data_cocos_and_expos01_20, &
                           ngwf_data_cocos_and_expos21_40, &
                           ngwf_data_cocos_and_expos41_60, &
                           ngwf_data_cocos_and_expos61_80, &
                           ngwf_data_cocos_and_expos81_103
      use utils, only: utils_alloc_check
      implicit none

      ! cks: <<local variables>>
      integer :: atom     ! atom counter
      integer :: n_shells ! number of shells
      integer :: n_prim   ! number of primitives in given shell
      integer :: ierr     ! allocation error flag

      ! lpl: Clear pointer before allocating if necessary
      if(allocated(gbasis)) then
         deallocate(gbasis,stat=ierr)
         call utils_dealloc_check('internal_ngwf_sto3g_create', &
              'gbasis',ierr)
      end if

      ! cks: allocate gbasis array for all supported elements
      allocate(gbasis(pub_num_gtatoms), stat=ierr)
      call utils_alloc_check('ngwf_data_sto3g_create','gbasis',ierr)

      ! cks: Number of shells per element, including one polarisation shell
      ! cks: This data was typed by Shyong Chen

      ! cks: Number of shells per element, including one polarisation shell
      ! cks: This data was typed by Shyong Chen
      gbasis(1)%nshells=2
      gbasis(2)%nshells=2
      gbasis(3)%nshells=4
      gbasis(4)%nshells=4
      gbasis(5)%nshells=4
      gbasis(6)%nshells=4
      gbasis(7)%nshells=4
      gbasis(8)%nshells=4
      gbasis(9)%nshells=4
      gbasis(10)%nshells=4
      gbasis(11)%nshells=6
      gbasis(12)%nshells=6
      gbasis(13)%nshells=6
      gbasis(14)%nshells=6
      gbasis(15)%nshells=6
      gbasis(16)%nshells=6
      gbasis(17)%nshells=6
      gbasis(18)%nshells=6
      gbasis(19)%nshells=8
      gbasis(20)%nshells=8
      gbasis(21)%nshells=9
      gbasis(22)%nshells=9
      gbasis(23)%nshells=9
      gbasis(24)%nshells=9
      gbasis(25)%nshells=9
      gbasis(26)%nshells=9
      gbasis(27)%nshells=9
      gbasis(28)%nshells=9
      gbasis(29)%nshells=9
      gbasis(30)%nshells=9
      gbasis(31)%nshells=9
      gbasis(32)%nshells=9
      gbasis(33)%nshells=9
      gbasis(34)%nshells=9
      gbasis(35)%nshells=9
      gbasis(36)%nshells=9
      gbasis(37)%nshells=11
      gbasis(38)%nshells=11
      gbasis(39)%nshells=12
      gbasis(40)%nshells=12
      gbasis(41)%nshells=12
      gbasis(42)%nshells=12
      gbasis(43)%nshells=12
      gbasis(44)%nshells=12
      gbasis(45)%nshells=12
      gbasis(46)%nshells=12
      gbasis(47)%nshells=12
      gbasis(48)%nshells=12
      gbasis(49)%nshells=12
      gbasis(50)%nshells=12
      gbasis(51)%nshells=12
      gbasis(52)%nshells=12
      gbasis(53)%nshells=12
      gbasis(54)%nshells=13
      gbasis(55)%nshells=14
      gbasis(56)%nshells=14
      gbasis(57)%nshells=15
      gbasis(58)%nshells=15
      gbasis(59)%nshells=15
      gbasis(60)%nshells=15
      gbasis(61)%nshells=15
      gbasis(62)%nshells=15
      gbasis(63)%nshells=15
      gbasis(64)%nshells=15
      gbasis(65)%nshells=15
      gbasis(66)%nshells=15
      gbasis(67)%nshells=15
      gbasis(68)%nshells=15
      gbasis(69)%nshells=15
      gbasis(70)%nshells=15
      gbasis(71)%nshells=15
      gbasis(72)%nshells=15
      gbasis(73)%nshells=15
      gbasis(74)%nshells=15
      gbasis(75)%nshells=15
      gbasis(76)%nshells=15
      gbasis(77)%nshells=15
      gbasis(78)%nshells=15
      gbasis(79)%nshells=15
      gbasis(80)%nshells=15
      gbasis(81)%nshells=15
      gbasis(82)%nshells=15
      gbasis(83)%nshells=15
      gbasis(84)%nshells=15
      gbasis(85)%nshells=15
      gbasis(86)%nshells=16
      gbasis(87)%nshells=16
      gbasis(88)%nshells=16
      gbasis(89)%nshells=16
      gbasis(90)%nshells=17
      gbasis(91)%nshells=18
      gbasis(92)%nshells=18
      gbasis(93)%nshells=18
      gbasis(94)%nshells=18
      gbasis(95)%nshells=17
      gbasis(96)%nshells=17
      gbasis(97)%nshells=17
      gbasis(98)%nshells=17
      gbasis(99)%nshells=17
      gbasis(100)%nshells=17
      gbasis(101)%nshells=17
      gbasis(102)%nshells=17
      gbasis(103)%nshells=17

      ! cks: allocate memory for the components of each shell
      do atom=1, pub_num_gtatoms

         n_shells = gbasis(atom)%nshells
         n_prim =3

         ! cks: allocate angular momentum memory
         allocate(gbasis(atom)%angmom(n_shells), stat=ierr)
         call utils_alloc_check('ngwfs_sto3g_create','gbasis%angmom',ierr)
         ! cks: initialise
         gbasis(atom)%angmom =-1


         ! cks: allocate exponents memory
         allocate(gbasis(atom)%expo(n_shells, n_prim), stat=ierr)
         call utils_alloc_check('ngwfs_sto3g_create','gbasis%expo',ierr)
         ! cks: initialise
         gbasis(atom)%expo =0.0_DP

         ! cks: allocate contraction coefficients memory
         allocate(gbasis(atom)%coco(n_shells, n_prim), stat=ierr)
         call utils_alloc_check('ngwfs_sto3g_create','gbasis%coco',ierr)
         ! cks: initialise
         gbasis(atom)%coco =0.0_DP

      enddo

      ! cks: load cocos and expos in the allocated memory
      call ngwf_data_cocos_and_expos01_20(gbasis(:))
      call ngwf_data_cocos_and_expos21_40(gbasis(:))
      call ngwf_data_cocos_and_expos41_60(gbasis(:))
      call ngwf_data_cocos_and_expos61_80(gbasis(:))
      call ngwf_data_cocos_and_expos81_103(gbasis(:))

    end subroutine internal_ngwf_sto3g_create

    subroutine internal_ngwf_sto3g_destroy

      !=======================================================!
      ! Deallocate all memory of gbasis array.                !
      !-------------------------------------------------------!
      ! lpl: Copied from ngwfs_mod (23-May-2011)              !
      !=======================================================!

      use ngwf_data, only: GTO_SET, pub_num_gtatoms
      use utils, only: utils_dealloc_check
      implicit none

      ! cks: <<local variables>>
      integer :: atom     ! atom counter
      integer :: ierr     ! allocation error flag

      ! cks: deallocate memory for the components of each shell
      do atom=1, pub_num_gtatoms

         ! cks: deallocate angular momentum memory
         deallocate(gbasis(atom)%angmom, stat=ierr)
         call utils_dealloc_check('ngwfs_sto3g_destroy','gbasis%angmom',ierr)

         ! cks: deallocate exponents memory
         deallocate(gbasis(atom)%expo, stat=ierr)
         call utils_dealloc_check('ngwfs_sto3g_destroy','gbasis%expo',ierr)

         ! cks: deallocate contraction coefficients memory
         deallocate(gbasis(atom)%coco, stat=ierr)
         call utils_dealloc_check('ngwfs_sto3g_destroy','gbasis%coco',ierr)

      enddo

      ! cks: deallocate the gbasis array
      deallocate(gbasis, stat=ierr)
      call utils_dealloc_check('ngwfs_sto3g_destroy','gbasis',ierr)

    end subroutine internal_ngwf_sto3g_destroy

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

 !--------------------------------------------------------------------------!
 ! FILE.47 output subroutines                                               !
 ! Prints matrices and other information in GENNBO's FILE.47 input format   !
 !    1. internal_write_nbo                                                 !
 !    2. internal_print_nbo_matrix                                          !
 !--------------------------------------------------------------------------!

    subroutine internal_write_nbo(tr_mat,file_suffix,ortho, &
         orig_atom_list,list_nat,ngwf_list,list_ngwfnum,scaling, &
         dbuffer_1,dbuffer_2,dbuffer_3, &
         tr_dm_input,tr_ovlp_input,print_tr_input,print_47_input)
      ! lpl: Transforms and prints matrices & other info to GENNBO FILE.47
      ! lpl: 24/06/2011: Added print_tr and routine to print actual
      !      transformation used to generate current matrices (prints out
      !      partial basis in the full NGWF basis to a seperate file)

      use constants, only: ANGSTROM, stderr
      use dense, only: dense_copy, dense_get_element, dense_scale
      use utils, only: utils_unit, utils_open_unit_check, &
           utils_close_unit_check

      implicit none

      ! lpl: Input
      type(DEM), intent(in) :: tr_mat ! Full transformation matrix
      character(len=*), intent(in) :: file_suffix ! File suffix

      logical, intent(in) :: ortho ! Flag to suppress printing of overlap
      integer, intent(in) :: orig_atom_list(:) ! List of atoms to be printed
      integer, intent(in) :: list_nat          ! (in the orig_iat index)

      integer, intent(in) :: ngwf_list(:) ! List of NGWFs to be printed
      integer, intent(in) :: list_ngwfnum

      real(kind=DP), intent(in) :: scaling ! Density matrix scaling

      ! lpl: Buffers
      ! lpl: dbuffer_1(:) doubles as transformed density matrix input
      !      dbuffer_3 doubles as transformed overlap matrix input
      type(DEM), intent(inout) :: dbuffer_1(:), dbuffer_3
      type(DEM), intent(inout) :: dbuffer_2
      logical, optional, intent(in) :: tr_dm_input, tr_ovlp_input
      ! lpl: Print transformation matrix
      logical, optional, intent(in) :: print_tr_input, print_47_input

      ! lpl: Local variables
      character(len=64)  :: nbo_extraflags, nbo_str1, nbo_str2
      character(len=256) :: nbo_output
      integer :: nbo_is, nbo_orig_iat, ingwf, list_it
      integer :: nbo_output_unit
      logical :: nbo_warn, tr_dm, tr_ovlp, print_tr, print_47
      real(kind=DP) :: nbo_rderr

      integer, allocatable :: rev_orig_atom_list(:) ! Global iat to local iat

      ! lpl : Input flags option
      tr_dm = .false.
      tr_ovlp = .false.
      if(present(tr_dm_input)) tr_dm = tr_dm_input
      if(present(tr_ovlp_input)) tr_ovlp = tr_ovlp_input

      print_tr = .false.
      if(present(print_tr_input)) print_tr = print_tr_input
      print_47 = .true.
      if(present(print_47_input)) print_47 = print_47_input

      ! lpl: Write everything apart from matrices to FILE.47
      if(print_47) then

         ! lpl: temporary - reverse atomic list (global to local)
         allocate(rev_orig_atom_list(pub_cell%nat),stat=ierr)
         call utils_alloc_check('npa_main','rev_orig_atom_list',ierr)

         rev_orig_atom_list = -1
         do list_it=1,list_nat
            nbo_orig_iat = orig_atom_list(list_it)
            rev_orig_atom_list(nbo_orig_iat) = list_it
         end do
         ! lpl: temporary - reverse atomic list (global to local)

         if(pub_on_root) then
            nbo_output = trim(pub_rootname)//'_'//file_suffix//'_nbo.47'
            nbo_output = adjustl(nbo_output)
            nbo_output_unit = utils_unit()
            open(unit=nbo_output_unit,form="formatted",file=trim(nbo_output), &
                 action="write",iostat=ierr)
            call utils_open_unit_check('npa_main','nbo_output',ierr)

            ! lpl: Announce writing to file in stdout
            write(stdout,'(3a)',advance='no') &
                 'Writing "',trim(adjustl(nbo_output)),'" ...'

            ! lpl: Trim blanks for legible output
            write(nbo_str1,*) list_nat
            nbo_str1 = adjustl(nbo_str1)
            write(nbo_str2,*) list_ngwfnum
            nbo_str2 = adjustl(nbo_str2)

            nbo_extraflags = ''
            if(pub_cell%num_spins > 1) &
                 nbo_extraflags = trim(adjustl(nbo_extraflags))//'  OPEN'
            if(ortho) nbo_extraflags = trim(adjustl(nbo_extraflags))//' ORTHO'
            nbo_extraflags = adjustl(nbo_extraflags)

            ! lpl: Write FILE.47 header
            write(nbo_output_unit,'(a)') ' $GENNBO  NATOMS='//trim(nbo_str1)// &
                  ' NBAS='//trim(nbo_str2)//'  '//trim(nbo_extraflags)//'  $END'

            ! lpl: Include CMO flag if full system is being printed
            nbo_extraflags = ''
            if(list_ngwfnum == ngwf_basis%num) nbo_extraflags = 'CMO'
            write(nbo_output_unit,'(a)') ' $NBO  FILE='//trim(pub_rootname)// &
                  '_'//trim(adjustl(file_suffix))//'  '// &
                       trim(adjustl(nbo_extraflags))//'  $END'

            ! lpl: Write atomic coords in the original input order
            nbo_warn = .false.
            write(nbo_output_unit,'(a)') ' $COORD'
            write(nbo_output_unit,'(a)') pub_rootname

            do list_it=1,list_nat

               nbo_orig_iat = orig_atom_list(list_it)

               ! lpl: NBO program requires integer core charges
               nbo_rderr = abs(elements(nbo_orig_iat)%ion_charge &
                    - nint(elements(nbo_orig_iat)%ion_charge))
               if(nbo_rderr /= 0) then ! Round Q to to nearest integer and warn
                  if(nbo_warn) then
                     write(stdout,*) 'Warning: Non-integer ionic charge &
                          &detected in the following atomic species:'
                     write(stdout,*) ''
                     write(stdout,*) ' Species   Label     Ionic Charge &
                          &  Rounded Val  |Rounding Err|'
                     write(stdout,*) &
                          ' ------------------------------------------------&
                               &----------------'
                     nbo_warn = .true.
                  end if
                  write(stdout,'(4x,a,5x,I2,7x,f10.5,6x,I3,9x,f11.6)') &
                        elements(nbo_orig_iat)%species_id, &
                        elements(nbo_orig_iat)%species_number, &
                        elements(nbo_orig_iat)%ion_charge, &
                        nint(elements(nbo_orig_iat)%ion_charge), &
                        nbo_rderr
               end if

               ! lpl: Write atomic coordinates (default ANGSTROMS)
               write(nbo_str1,*) elements(nbo_orig_iat)%atomic_number
               nbo_str1 = adjustl(nbo_str1)
               write(nbo_output_unit,'(a,I3,3(1x,f12.6))',iostat=ierr) &
                     '   '//trim(nbo_str1)//' ', &
                     nint(elements(nbo_orig_iat)%ion_charge),  &
                     elements(nbo_orig_iat)%centre%x/ANGSTROM, &
                     elements(nbo_orig_iat)%centre%y/ANGSTROM, &
                     elements(nbo_orig_iat)%centre%z/ANGSTROM
               if(ierr /= 0) then
                  write(stderr,'(a)') 'Error in internal_write_nbo : &
                        &writing elements(row)%centre%x/ANGSTROM failed'
                  call comms_abort
               end if
            end do

            if(nbo_warn) then
               write(stdout,'(a)') ' --------------------------------------&
                    &-----------------------------'
               write(stdout,'(a)') 'WARNING: Ionic charges were rounded up &
                    &to the nearest integer in FILE.47.'
               write(stdout,'(a)') 'WARNING: Proceed with caution in your &
                    &NBO analysis.'
            end if

            write(nbo_output_unit,'(a)') ' $END'
            ! lpl: END Write atomic coords in the original input order

            ! lpl: Write atomic centres w.r.t. NGWF rep based on original
            !      input order
            write(nbo_output_unit,'(a)') ' $BASIS'
            do list_it=1,list_ngwfnum

               ingwf = ngwf_list(list_it)

               if(list_it == 1) then
                  write(nbo_output_unit,'(a)',advance="no") '  CENTER = '
               else if(MOD(list_it,16) == 1) then
                  write(nbo_output_unit,'(/a)',advance="no") '           '
               end if

               write(nbo_str1,*) rev_orig_atom_list(atom_order(ingwf))
               nbo_str1=adjustl(nbo_str1)
               write(nbo_output_unit,'(3x,a)',advance="no") trim(nbo_str1)

            end do
            write(nbo_output_unit,'(a)') ''

            ! lpl: Write (false) Y(l,m) labels for NGWFs
            do list_it=1,list_ngwfnum

               ingwf = ngwf_list(list_it)

               if(list_it == 1) then
                  write(nbo_output_unit,'(a)',advance="no") '  LABEL  =  '
               else if(MOD(list_it,16) == 1) then
                  write(nbo_output_unit,'(/a)',advance="no") '            '
               end if

               write(nbo_output_unit,'(I3,1x)',advance="no") &
                    ngwf_lm_label(ingwf)

            end do

            write(nbo_output_unit,'(/a)') ' $END'
            ! lpl: END Write atomic centres w.r.t. NGWF rep based on original
            !      input order
         end if ! END if(pub_on_root)
         ! lpl: END Write everything apart from matrices to FILE.47

         call comms_barrier

         ! lpl: Transform and write matrices

         ! lpl: Write covariant density matrix [dbuffer_1(:)]
         ! lpl: If tr_dm = .false., perform density matrix transformation and
         !      store in dbuffer_1(:)
         if(.not. tr_dm) then
            do nbo_is=1,pub_cell%num_spins
               call dense_convert(dbuffer_1(nbo_is),pure_dm(nbo_is))
               call internal_dense_transform(dbuffer_1(nbo_is), &
                    full_tr,dbuffer_2)
            end do
         end if

         ! lpl: QC for 1st 5 column element of row 1 of partial scaled DM
         if(print_qc .and. file_suffix == 'nao') then
            if(pub_on_root) then
               write(stdout,'(a)') ''
               do list_it=1,min(5,list_ngwfnum)
                  call dense_get_element(nbo_rderr,dbuffer_1(1),1, &
                       ngwf_list(list_it))
                  write(stdout,'(a,I1,a,f14.9)') &
                     '<QC> [scaled_nao_dm_el_',list_it,']: ',nbo_rderr
               end do
            end if
         end if

         ! lpl: Scale and print density matrix
         if(pub_on_root) write(nbo_output_unit,'(a)',advance="no") ' $DENSITY'
         do nbo_is=1,pub_cell%num_spins
            ! lpl: Scale density matrix if required
            call dense_scale(dbuffer_1(nbo_is),scaling)
            call internal_print_nbo_matrix(dbuffer_1(nbo_is),ngwf_list, &
                 list_ngwfnum,nbo_output_unit,'DENSITY')
         end do
         if(pub_on_root) write(nbo_output_unit,'(/a)') ' $END'

         ! lpl: Write overlap matrix
         ! lpl: If tr_ovlp = .false., perform overlap matrix transformation and
         !      store in dbuffer_3
         if(.not. tr_ovlp) then
            call dense_convert(dbuffer_3,rep%overlap)
            call internal_dense_transform(dbuffer_3,full_tr,dbuffer_2)
         end if

         ! lpl: If ortho = .false., print overlap matrix in FILE.47
         if(.not. ortho) then
            if(pub_on_root) write(nbo_output_unit,'(a)',advance="no") &
                 ' $OVERLAP'
            call internal_print_nbo_matrix(dbuffer_3,ngwf_list,list_ngwfnum, &
                 nbo_output_unit,'OVERLAP')
            if(pub_on_root) write(nbo_output_unit,'(/a)') ' $END'
         end if

         ! lpl: Transform and write kinetic integral matrix
         call dense_convert(dbuffer_2,rep%kinet)
         call internal_dense_transform(dbuffer_2,full_tr,dbuffer_1(1))

         if(pub_on_root) write(nbo_output_unit,'(a)',advance="no") ' $KINETIC'
         call internal_print_nbo_matrix(dbuffer_2,ngwf_list,list_ngwfnum, &
              nbo_output_unit,'KINETIC')
         if(pub_on_root) write(nbo_output_unit,'(/a)') ' $END'

         ! lpl: Transform and write Hamiltonian, and calculate CMO
         if(pub_on_root) write(nbo_output_unit,'(a)',advance="no") ' $FOCK'
         do nbo_is=1,pub_cell%num_spins
            call dense_convert(dbuffer_2,ham%ham(nbo_is))
            call internal_dense_transform(dbuffer_2,full_tr,dbuffer_1(1))

            ! lpl: Transformed Hamiltonian for this spin in dbuffer_2
            call internal_print_nbo_matrix(dbuffer_2,ngwf_list,list_ngwfnum, &
                 nbo_output_unit,'FOCK')

            ! lpl: Diagonalize H to get CMO [ HM = SME ] if full system is
            !      being printed
            !      dbuffer_2 = H; dbuffer_1(nbo_is) = eigenvectors for this spin
            !      dbuffer_3 = S; col_buffer = eigenvalues (not needed) buffer
            ! lpl: Diagonalize H if full system is to be printed and store
            !      MO in dbuffer_1(nbo_is)
            if(list_ngwfnum == ngwf_basis%num) then
               ! lpl: Copy S in 'nbo_is+1''th buffer unless this is the
               !      last spin index, then use S directly. For some reason
               !      S gets disturbed after dense_eigensolve
               if(nbo_is < pub_cell%num_spins) then
                  call dense_copy(dbuffer_1(nbo_is+1),dbuffer_3)
                  call dense_eigensolve(list_ngwfnum,col_buffer(:),dbuffer_2, &
                       dbuffer_1(nbo_is+1),1,dbuffer_1(nbo_is))
               else
                  call dense_eigensolve(list_ngwfnum,col_buffer(:),dbuffer_2, &
                       dbuffer_3,1,dbuffer_1(nbo_is))
               end if

#ifdef DEBUG
               if(pub_on_root) then
                  write(debug_output_unit,'(a)') &
                       'CMO Energies for '&
                       //trim(adjustl(file_suffix))//' (a.u.)'
                  write(debug_output_unit,'(a,I3)') 'Spin ',nbo_is
                  write(debug_output_unit,'(a)') '-----------------------------'
                  do list_it=1,list_ngwfnum
                     write(debug_output_unit,'(1x,I3,1x,f14.7)') &
                          list_it,col_buffer(list_it)
                  end do
                  write(debug_output_unit,'(a)') '-----------------------------'
               end if
               call comms_barrier
#endif

            end if ! END if(list_ngwfnum == ngwf_basis%num)
         end do ! END nbo_is=1,pub_cell%num_spins
         if(pub_on_root) write(nbo_output_unit,'(/a)') ' $END'

         ! lpl: Write canonical MO if full system is being printed
         if(list_ngwfnum == ngwf_basis%num) then
            if(pub_on_root) write(nbo_output_unit,'(a)',advance="no") ' $LCAOMO'
            do nbo_is=1,pub_cell%num_spins
                  call internal_print_nbo_matrix(dbuffer_1(nbo_is),ngwf_list, &
                       list_ngwfnum,nbo_output_unit,'LCAOMO',cmo=.true.)
            end do
            if(pub_on_root) write(nbo_output_unit,'(/a)') ' $END'
         end if

         ! lpl: Close FILE.47
         if(pub_on_root) then
            close(unit=nbo_output_unit,iostat=ierr)
            call utils_close_unit_check('npa_main','nbo_output_unit',ierr)
            write(stdout,'(a)') ' done'
         end if

         ! lpl: Write partial atomic index to orig_iat mapping only if
         !      partial system is being printed
         if(list_ngwfnum < ngwf_basis%num) then
            if(pub_on_root) then
               nbo_output = trim(pub_rootname)//'_'//file_suffix// &
                    '_atomindex.dat'
               nbo_output = adjustl(nbo_output)
               nbo_output_unit = utils_unit()

               open(unit=nbo_output_unit,form="formatted", &
                    file=trim(nbo_output),action="write",iostat=ierr)
               call utils_open_unit_check('npa_main','nbo_output',ierr)
               write(stdout,'(3a)',advance='no') &
                    'Writing "',trim(adjustl(nbo_output)),'" ...'

               write(nbo_output_unit,'(a,f14.7)') &
                    'Density matrix scaling: ',scaling
               write(nbo_output_unit,'(a)') '---------------------------------&
                    &---------------------------'
               write(nbo_output_unit,'(a)') 'Original atomic index to partial &
                    &matrix atomic index mapping'
               write(nbo_output_unit,'(a)') '---------------------------------&
                    &---------------------------'
               write(nbo_output_unit,'(a)') ' Partial index     Original index '
               write(nbo_output_unit,'(a)') '----------------------------------'

               do list_it=1,list_nat
                  write(nbo_output_unit,'(5x,I5,13x,I5)') &
                       list_it,orig_atom_list(list_it)
               end do
               write(nbo_output_unit,'(a)') '---------------------------------&
                    &---------------------------'
               close(unit=nbo_output_unit,iostat=ierr)
               call utils_close_unit_check('npa_main','nbo_output_unit',ierr)
               write(stdout,'(a)') ' done'
            end if
         end if

         deallocate(rev_orig_atom_list,stat=ierr)
         call utils_dealloc_check('npa_main','rev_orig_atom_list',ierr)

      end if ! END if(print_47)

      ! lpl: 23/06/2011 - Write NGWF to (partial) AO transformation
      !      ngwf_list = sorted AO index
      !      ngwf_orig_list = sorted NGWF index
      if(print_tr) then
         if(pub_on_root) then
            nbo_output = trim(pub_rootname)//'_inittr_'// &
                 file_suffix//'_nbo.dat'
            nbo_output = adjustl(nbo_output)
            nbo_output_unit = utils_unit()

            open(unit=nbo_output_unit,form="formatted", &
                 file=trim(nbo_output),action="write",iostat=ierr)
            call utils_open_unit_check('npa_main','nbo_output',ierr)

            ! lpl: Write header
            write(nbo_output_unit,'(a,/a)') trim(pub_rootname)//'_'// &
                 file_suffix//'_nbo.47','NGWF to AO transformation'
            ! lpl: Write info used in npa_plotnbo
            write(nbo_output_unit,'(1x,I5,1x,I5)',advance='no') &
                 ngwf_basis%num,list_ngwfnum

            ! lpl: Write (partial) NGWF to AO transformation
            do icol=1,list_ngwfnum ! col = AO index (orig order)
               it = 0
               do irow=1,ngwf_basis%num ! row = NGWF index (orig order)
                  call dense_get_element(mtx_el,full_tr, &
                       ngwf_orig_list(irow),ngwf_list(icol))
                  if(MOD(it,5) == 0) write(nbo_output_unit,'(a)') ''
                  write(nbo_output_unit,'(1x,E21.14)',advance='no') mtx_el
                  it = it + 1
               end do
            end do

            close(unit=nbo_output_unit,iostat=ierr)
            call utils_close_unit_check('npa_main','nbo_output',ierr)
         end if ! END if(pub_on_root)
      end if

      call comms_barrier

    end subroutine internal_write_nbo

    subroutine internal_print_nbo_matrix(mat,list,num,out_unit,errname,cmo)
      ! lpl: Subroutine to print matrices in FILE.47 format
      implicit none

      type(DEM),  intent(in) :: mat
      integer,    intent(in) :: list(:)
      integer,    intent(in) :: num
      integer,    intent(in) :: out_unit
      character(len=*), intent(in) :: errname

      ! lpl: Flag to Write CMO eigenvectors (skip column sorting)
      logical, optional, intent(in) :: cmo

      logical :: write_cmo
      integer :: it, irow, icol, ierr
      real(kind=DP) :: mtx_el

      ierr = 0

      it = 0
      write_cmo = .false.

      ! lpl: Deal with optional cmo flag
      if(present(cmo)) write_cmo = cmo

      if(pub_on_root) then
         if(write_cmo) then ! lpl: Write transpose
            do icol=1,num
               do irow=1,num
                  call dense_get_element(mtx_el,mat,list(irow),icol)
                  if(MOD(it,5) == 0) write(out_unit,'(a)') ' '
                  write(out_unit,'(a,e14.7)',advance="no",iostat=ierr) &
                       ' ',mtx_el
                  it = it + 1
                  if(ierr /= 0) then
                     write(stdout,'(a)') 'ERROR: Writing '//errname// &
                          ' block failed'
                     call comms_abort
                  end if
               end do
            end do
         else ! lpl: Write normal matrix
            do irow=1,num
               do icol=1,num
                  call dense_get_element(mtx_el,mat,list(irow),list(icol))
                  if(MOD(it,5) == 0) write(out_unit,'(a)') ' '
                  write(out_unit,'(a,e14.7)',advance="no",iostat=ierr) &
                       ' ',mtx_el
                  it = it + 1
                  if(ierr /= 0) then
                     write(stdout,'(a)') 'ERROR: Writing '//errname// &
                          ' block failed'
                     call comms_abort
                  end if
               end do
            end do
         end if
      end if

    end subroutine internal_print_nbo_matrix

    subroutine internal_dense_block_dsyev(eigenvecs,eigenvals, &
         start_idx,end_idx)
       ! lpl: Wrapper for LAPACK routine to solve normal eigenvalue equation
       !      AV=VL for potentially sub-matrix block of full DEM
       use wrappers, only: wrappers_dsyev_lt

       implicit none
       type(DEM), intent(inout) :: eigenvecs
       real(kind=DP), intent(inout) :: eigenvals(:)
       integer, intent(in) :: start_idx, end_idx

       if(eigenvecs%iscmplx) then
          if(pub_on_root) write(stdout,'(a)') 'ERROR: Complex matrix given in &
               &internal_dense_block_dsyev in npa_mod'
          call comms_abort
       else
          call wrappers_dsyev_lt&
               (eigenvecs%dmtx(start_idx:end_idx,start_idx:end_idx), &
               eigenvals(start_idx:end_idx),end_idx-start_idx+1)
       end if
       ! lpl: Synchronize degenerate eigenvectos (from dense_eigensolve)
       call comms_bcast(pub_root_node_id,eigenvecs%dmtx, &
            eigenvecs%nrows*eigenvecs%mcols)

       call comms_barrier

    end subroutine internal_dense_block_dsyev

    subroutine internal_pure_dm_init
      ! lpl: Minimize clutter in main code
      implicit none

      integer :: is

      if(allocated(pure_dm)) then
         if(pub_on_root) write(stdout,'(a)') 'ERROR: pure_dm is already &
              &allocated in internal_pure_dm_init in npa_main'
         call comms_abort
      end if

      allocate(pure_dm(pub_cell%num_spins),stat=ierr)
      call utils_alloc_check('internal_pure_dm_init','pure_dm',ierr)
      do is=1,pub_cell%num_spins
         call sparse_create(pure_dm(is),denskern(is))
      end do

    end subroutine internal_pure_dm_init

    subroutine internal_pure_dm_destroy
      ! lpl: Minimize clutter in main code
      implicit none

      integer :: is

      if(.not.allocated(pure_dm)) then
         if(pub_on_root) write(stdout,'(a)') 'ERROR: pure_dm is already &
              &unallocated in internal_pure_dm_init in npa_main'
         call comms_abort
      end if

      do is=pub_cell%num_spins,1,-1
         call sparse_destroy(pure_dm(is))
      end do
      deallocate(pure_dm,stat=ierr)
      call utils_dealloc_check('internal_pure_dm_destroy','pure_dm',ierr)

    end subroutine internal_pure_dm_destroy

    subroutine internal_get_puredm(pure_dm,denskern,rep)
      ! lpl: Purify K and transform to covariant DM 'pure_dm'
      use kernel, only: kernel_purify

      implicit none

      type(SPAM3), intent(inout) :: pure_dm(pub_cell%num_spins)
      type(SPAM3), intent(inout) :: denskern(pub_cell%num_spins)
      type(NGWF_REP), intent(in) :: rep

      type(SPAM3), allocatable :: purekern(:)
      type(SPAM3) :: sbuffer_ks

      integer :: is

      sbuffer_ks%structure='KS'
      call sparse_create(sbuffer_ks)

      allocate(purekern(pub_cell%num_spins),stat=ierr)
      call utils_alloc_check('npa_main','purekern',ierr)
      do is=1,pub_cell%num_spins
         call sparse_create(purekern(is),denskern(is))
      end do

      ! Purify density kernel
      call kernel_purify(purekern,denskern,rep%overlap, &
           rep%inv_overlap,rep%n_occ)

      ! Get covariant pure_dm(is)
      do is=1,pub_cell%num_spins
         call sparse_product(sbuffer_ks,purekern(is),rep%overlap)
         call sparse_product(pure_dm(is),rep%overlap,sbuffer_ks)
      end do
      if(pub_cell%num_spins == 1) call sparse_scale(pure_dm(1),2.0_DP)

      do is=pub_cell%num_spins,1,-1
         call sparse_destroy(purekern(is))
      end do
      deallocate(purekern,stat=ierr)
      call utils_dealloc_check('npa_main','purekern',ierr)

      call sparse_destroy(sbuffer_ks)
      call comms_barrier

    end subroutine internal_get_puredm

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!

  end subroutine npa_main

  subroutine npa_read_nbodata(nbongwf_tr,ngwf_basis,nbo_trfile,pre_ao,num_ao)
    ! lpl: Calculates NGWF to NBO transformation. Called in properties_mod
    !      to plot NBOs
    use comms, only: comms_abort, comms_barrier, comms_bcast, pub_on_root, &
         pub_root_node_id
    use constants, only: DP, stdout
    use dense, only: DEM, dense_create, dense_destroy, dense_product, &
         dense_put_element, dense_get_element, dense_scale
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: pub_orig_atom, pub_distr_atom
    use simulation_cell, only: pub_cell
    use rundat, only: pub_rootname
    use utils, only: utils_unit, utils_alloc_check, utils_dealloc_check, &
         utils_open_unit_check, utils_close_unit_check

    implicit none

    ! lpl: Input/Output
    type(DEM), intent(inout) :: nbongwf_tr ! NGWF to NBO transformation
    type(FUNC_BASIS), intent(in) :: ngwf_basis ! Run-time NGWF basis
    ! lpl: Orbital type in AO basis file generated by GENNBO
    character(len=*), intent(in) :: nbo_trfile
    logical, intent(in) :: pre_ao ! If true, inittr = PNAO, not NAO
    integer, intent(out) :: num_ao ! For index checking in properties_mod

    ! lpl: Local variables
    character(len=256) :: c_buff
    character(len=64)  :: ngwfao_input_file

    integer :: it, lc_it, irow, icol
    integer :: iat, orig_iat
    integer :: nbo_input_unit, ierr
    integer :: num_ngwfs

    ! lpl: NGWF index mapping
    integer, allocatable :: orig_to_runtime(:)

    real(kind=DP), allocatable :: col_buffer(:)

    ! lpl: Transformation matrices
    type(DEM) :: aonbo_tr
    type(DEM) :: init_tr

#ifdef DEBUG
    if(pub_on_root) write(stdout,'(a)') 'Entering npa_read_nbodata'
#endif

    ierr = 0
    ! lpl: Read in appropraite NGWF to AO transformation
    if(pre_ao) then
       ngwfao_input_file = trim(pub_rootname)//'_inittr_pnao_nbo.dat'
    else
       ngwfao_input_file = trim(pub_rootname)//'_inittr_nao_nbo.dat'
    end if

    ! lpl: Generate orig_to_runtime NGWF mapping
    allocate(orig_to_runtime(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('npa_read_nbodata','orig_to_runtime',ierr)
    orig_to_runtime = 0

    icol = 0
    do orig_iat=1,pub_cell%nat
       iat = pub_distr_atom(orig_iat)
       do lc_it=1,ngwf_basis%num_on_atom(iat)
          it = ngwf_basis%first_on_atom(iat) + lc_it - 1
          icol = icol + 1
          orig_to_runtime(icol) = it
       end do
    end do

    allocate(col_buffer(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('npa_read_nbodata','col_buffer',ierr)

    ! lpl: Explicitly initialize nbongwf_tr with zeroes
    do irow=1,ngwf_basis%num
       do icol=1,ngwf_basis%num
          call dense_put_element(0.0_DP,nbongwf_tr,irow,icol)
       end do
    end do

    if(nbo_trfile == 'NONE') then ! If we don't need to read GENNBO tr files
       if(pub_on_root) then
          nbo_input_unit = utils_unit()

          open(unit=nbo_input_unit,file=trim(ngwfao_input_file), &
               form='formatted',action='read',iostat=ierr)
          call utils_open_unit_check('npa_read_nbodata','inittr',ierr)

          ! lpl: Read basis number info & perform sanity check
          do it=1,2
             read(nbo_input_unit,*) c_buff
          end do
          read(nbo_input_unit,*) num_ngwfs,num_ao
          if(num_ngwfs /= ngwf_basis%num) then
             write(stdout,'(a)') 'ERROR: num_ngwfs /= ngwf_basis%num'
             call comms_abort
          end if
          if(num_ao > num_ngwfs) then
             write(stdout,'(a)') 'ERROR: num_ao > num_ngwfs'
             call comms_abort
          end if
       end if ! END if(pub_on_root)
       call comms_bcast(pub_root_node_id,num_ao)

       ! lpl: Read and sort nbongwf_tr matrix from NGWF to AO
       do icol=1,num_ao ! irow = AO index
          if(pub_on_root) read(nbo_input_unit,*) col_buffer(:)
          call comms_bcast(pub_root_node_id,col_buffer)
          call comms_barrier
          do irow=1,ngwf_basis%num ! irow = NGWF index
             call dense_put_element(col_buffer(irow), &
                  nbongwf_tr,orig_to_runtime(irow),icol)
          end do
       end do

       ! lpl: Close inittr
       if(pub_on_root) then
          close(unit=nbo_input_unit,iostat=ierr)
          call utils_close_unit_check('npa_read_nbodata','nbo_input_unit',ierr)
       end if

    else ! else read nbo_trfile

       ! lpl: Create & zero dense matrices
       call dense_create(init_tr,ngwf_basis%num,ngwf_basis%num)
       call dense_create(aonbo_tr,ngwf_basis%num,ngwf_basis%num)

       do irow=1,ngwf_basis%num
          do icol=1,ngwf_basis%num
             call dense_put_element(0.0_DP,init_tr,irow,icol)
             call dense_put_element(0.0_DP,aonbo_tr,irow,icol)
          end do
       end do

       call comms_barrier

       if(pub_on_root) then
          nbo_input_unit = utils_unit()

          open(unit=nbo_input_unit,file=trim(ngwfao_input_file), &
               form='formatted',action='read',iostat=ierr)
          call utils_open_unit_check('npa_read_nbodata','inittr',ierr)

          ! lpl: Read basis number info & perform sanity check
          do it=1,2
             read(nbo_input_unit,*) c_buff
          end do
          read(nbo_input_unit,*) num_ngwfs,num_ao
          if(num_ngwfs /= ngwf_basis%num) then
             write(stdout,'(a)') 'ERROR: num_ngwfs /= ngwf_basis%num'
             call comms_abort
          end if
          if(num_ao > num_ngwfs) then
             write(stdout,'(a)') 'ERROR: num_ao > num_ngwfs'
             call comms_abort
          end if
       end if ! END if(pub_on_root)
       call comms_bcast(pub_root_node_id,num_ao)

       ! lpl: Read and sort init_tr matrix from NGWF to AO
       do irow=1,num_ao ! irow = AO index
          if(pub_on_root) read(nbo_input_unit,*) col_buffer(:)
          call comms_bcast(pub_root_node_id,col_buffer)
          call comms_barrier
          do icol=1,ngwf_basis%num ! icol = NGWF index
             call dense_put_element(col_buffer(icol), &
                  init_tr,irow,orig_to_runtime(icol))
          end do
       end do

       ! lpl: Close inittr and open GENNBO FILE.xx
       if(pub_on_root) then
          close(unit=nbo_input_unit,iostat=ierr)
          call utils_close_unit_check('npa_read_nbodata','nbo_input_unit',ierr)

          ! lpl: Open FILE.xx
          open(unit=nbo_input_unit,file=trim(nbo_trfile), &
               form='formatted',action='read',iostat=ierr)
          call utils_open_unit_check('npa_read_nbodata','nbo_filexx',ierr)

          ! lpl: Skip headers in FILE.xx
          do it=1,3
             read(nbo_input_unit,*) c_buff
          end do
       end if

       ! lpl: Read orbital vector comp in AO basis from FILE.xx
       do irow=1,num_ao ! irow = NBO's orbital index
          col_buffer = 0.0_DP
          if(pub_on_root) read(nbo_input_unit,*,iostat=ierr) &
               col_buffer(1:num_ao)
          call comms_bcast(pub_root_node_id,col_buffer)
          call comms_barrier
          do icol=1,num_ao ! icol = AO index
             call dense_put_element(col_buffer(icol),aonbo_tr,irow,icol)
          end do

          if(pub_on_root) then
             if(ierr /= 0) then
                write(stdout,'(a,I5)') 'WARNING: EOF reached in FILE.xx when &
                     &attempting to read orbital # ', it
                write(stdout,'(a)') 'This means # orbitals < nbas. Make sure &
                     &non-existent orbitals','are not selected for plotting.'
                exit
             end if
          end if ! END if(pub_on_root)
       end do

       ! lpl: Close FILE.xx
       if(pub_on_root) then
          close(unit=nbo_input_unit,iostat=ierr)
          call utils_close_unit_check('npa_read_nbodata','nbo_input_unit',ierr)
       end if

       ! lpl: Get NGWF to NBO's orbital transformation
       call dense_scale(nbongwf_tr,0.0_DP)
       call dense_product(nbongwf_tr,init_tr,aonbo_tr,transpose_amat=.true., &
            transpose_bmat=.true.)

       call dense_destroy(aonbo_tr)
       call dense_destroy(init_tr)

    end if ! END if(nbo_trfile == 'NONE')

    deallocate(orig_to_runtime,stat=ierr)
    call utils_dealloc_check('npa_read_nbodata','orig_to_runtime',ierr)
    deallocate(col_buffer,stat=ierr)
    call utils_dealloc_check('npa_read_nbodata','col_buffer',ierr)

    call comms_barrier

#ifdef DEBUG
    if(pub_on_root) write(stdout,'(a)') 'Exiting npa_read_nbodata'
#endif

  end subroutine npa_read_nbodata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module npa
