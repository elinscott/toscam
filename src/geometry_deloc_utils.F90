! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!              G E O M E T R Y _ D E L O C  _ U T I L S                       !
!=============================================================================!
!                                                                             !
! $Id: geometry_deloc_utils.F90,v 1.8 2010/05/21 13:30:32 vmilman Exp $         !
!                                                                             !
!-----------------------------------------------------------------------------!
! This module contains utilities for geometry optimization in delocalized     !
! coordinates. A similar module will be required in any other application     !
! that wishes to use geometry_deloc_algor module: i.e., the public routines   !
! provided in this module would have to be implemented in other codes.        !
!                                                                             !
! Note that the data structures in this module should be copied by other apps !
!                                                                             !
!-----------------------------------------------------------------------------!
! Written from  "Module Specification for New Plane Wave Code",  M. Segall,   !
!    P. Lindan, M. Probert, C. Pickard, P. Hasnip, S. Clark and M. Payne      !
!                           Copyright 1999/2000                               !
!-----------------------------------------------------------------------------!
! Written by Victor Milman, v0.1, 27/03/2003                                  !
!-----------------------------------------------------------------------------!
! modification information                                                    !
!=============================================================================!

module geometry_deloc_utils

  use constants, only : dp,pi,periodic_table_name,ANGSTROM
  use comms, only : comms_bcast,pub_root_node_id,pub_on_root,comms_abort
  use simulation_cell, only: castep_cell_cart_to_frac, castep_model, &
       castep_model_cell_changed

  implicit none                                 !Impose strong typing

  private                                       !Everything is private ...

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: deloc_utils_initialize
  public :: deloc_utils_deallocate
  public :: deloc_utils_nullify
  public :: deloc_utils_energy_gradient
  public :: deloc_utils_rationalize_coords
  public :: deloc_utils_mdl_to_internals
  public :: deloc_utils_internals_to_mdl
  public :: deloc_utils_io_abort
  public :: deloc_utils_output_converged
  public :: deloc_utils_geom_converged
  public :: deloc_utils_update_params
  public :: deloc_utils_read_hessian        ! rdhessian
  public :: deloc_utils_write_hessian       ! wrhessian
  public :: deloc_utils_read_di_data        ! rdpchk
  public :: deloc_utils_read_np_p           ! short version of rdpchk
  public :: deloc_utils_write_di_data       ! wrpchk
  public :: deloc_utils_write_trajectory    ! opttbl_iter (dummy for CASTEP purposes)
  public :: deloc_utils_get_alt_bond_list   ! getMDF
  public :: deloc_utils_pack_atom_indices   ! pacel
  public :: deloc_utils_unpack_atom_indices ! unpacel
  public :: deloc_utils_save_structure      ! dummy for wrcoord
  public :: deloc_utils_read_opt_mode       ! dummy for rdoptmode
  public :: deloc_utils_read_DI_constraints ! dummy for rdcon_p
  public :: deloc_utils_read_fix_atoms      ! dummy for rdfixat
  public :: deloc_utils_read_fix_cartesians ! dummy for rdfixp
  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  type(castep_model),pointer,public,save  :: mdl_ptr
  logical,public,save                     :: di_on_root
  integer,public,save                     :: ip_error
  integer,public,save                     :: di_stdout
  integer,public,save                     :: di_iprint
  real(kind=dp),public,save               :: length_conv
  real(kind=dp),public,save               :: bohr_di
  real(kind=dp),public,save               :: energy_conv
  real(kind=dp),public,save               :: force_conv
  character (len=30),public,save  :: di_energy_label
  character (len=30),public,save  :: di_force_label
  character (len=30),public,save  :: di_pressure_label
  character (len=30),public,save  :: di_length_label

  logical,public,save :: use_deloc_int_output ! if TRUE, then printout from DI module
                                              ! is used. if FALSE, then all information
                                              ! about convergence and structure is printed out
                                              ! by the driver-application (e.g., CASTEP)

  integer,parameter,public           :: di_dp = dp
  real(kind=dp),parameter,public     :: di_pi = pi

  integer,public,save                :: ncycle  ! number of cycles in DI minimizer
                                                ! 0 just to verify that delocalized can be done
                                                !   make internals and store (memory or PCHK)
                                                ! 1 finish delocalized with accurate gradients
                                                !   make final delocalized, deal with starting Hessian
                                                !>1 continue, just update Hessian  
  integer,public,save              :: number_atoms           ! # of real atoms
  integer,public,save              :: lattice_dimension = 3  ! only 3D systems
  integer,public,save              :: maxdiis = 0            ! don't use DIIS for now
  integer,public,save              :: ncons_uns              ! # of unsatisfied internal constraints
  integer,public,save              :: ncons                  ! # of internal constraints
  integer,public,save              :: ioptc = 2              ! type of optimizer (2=deloc)
  integer,public,save              :: iupdat                 ! type of Hessian update
                                                             ! 0 - do not update the Hessian (!?)
                                                             ! 1 - Powell update (default for TS search)
                                                             ! 2 - BFGS update (default for minimization)
                                                             ! 3 - BFGS with safeguards to ensure retention of
                                                             !     positive definiteness (default for GDIIS)
                                                             ! 4 - Murtagh-Sargent update
                                                             ! 5 - Powell/Murtagh-Sargent update
  real(kind=dp),public,save               :: tolg = 0.001_dp        ! gradient tolerance
  real(kind=dp),public,save               :: told = 0.001_dp        ! displacement tolerance
  real(kind=dp),public,save               :: tole = 0.0001_dp       ! energy tolerance

  integer,public,save                     :: num_fixed_coords       ! number of fixed Cartesian coordinates
  integer,public,save                     :: num_fixed_constr_rigid_body ! number of internals from num_rigid_body_atoms atoms
  integer,public,save                     :: num_rigid_body_atoms   ! number of atoms with all the internals frozen
  integer,public,save                     :: num_dummies            ! number of dummy atoms
  integer,public,save                     :: ntrans                 ! number of symmetry operations
  logical,public,save                     :: xc_redressed = .false. ! did we rationalize coords?
  real (kind=dp),dimension(1:3,1:3),public,save :: elat             ! real space lattice
  real (kind=dp),dimension(1:3,1:3),public,save :: elatb            ! reciprocal space lattice
  real (kind=dp),dimension(1:3,1:3),public,save :: elatp            ! diagonalized real space basis
  integer,dimension(1:3),public,save      :: ipvlt                  ! pivoting for real-space basis diagonalization

  !------------------------------ logical switches
  logical,public,save                :: lp_mdf = .false.       ! only mdf read
  logical,public,save                :: lp_mdf_del = .false.   ! mdf read follow by automatic internal generation
  logical,public,save                :: lp_del = .true.        ! only automatic internal generation
  logical,public,save                :: lp_discf = .true.      ! automatic disconnected fragments
  logical,public,save                :: lp_idb = .false.       ! derivative terms in Hessian transform
  logical,public,save                :: lp_del_test = .false.  ! only test delocalized
  logical,public,save                :: ldokpt
  logical,public,save                :: period
  logical,public,save                :: tsflag                 ! true.  - TS search, .false. - minimization
  logical,public,save                :: symflag
  logical,public,save                :: l_scan_pes
  logical,public,save                :: ldoqst
  logical,public,save                :: ldoneb
  logical,public,save                :: global_data_exists = .false.  ! do we have a memory backup?

  real(kind=dp),public,save               :: ec                     ! current total energy
  real(kind=dp),public,save               :: eold                   ! previous total energy

  ! arrays
  integer,dimension(:,:),allocatable,public,save        :: np_p   ! for rationalizing coordinates
  real(kind=dp),dimension(:,:),allocatable,public,save  :: coords_cart     ! cartesian coordinates
  real(kind=dp),dimension(:,:),allocatable,public,save  :: gradient_cart     ! cartesian gradients
  real(kind=dp),dimension(:,:),allocatable,public,save  :: hess_cartesian ! cartesian hessian (the full 3*NATOMS by 3*NATOMS matrix)
  real(kind=dp),dimension(:,:,:),allocatable,public,save  :: trans  ! symmetry operations
  integer,dimension(:,:),allocatable,public,save        :: neqatm ! list of symmetry equivalent atoms
  integer,dimension(:),allocatable,public,save          :: map_fixed_coords
  integer,dimension(:,:),allocatable,public,save        :: list_fixed_coords ! 0 - coordinate active, 1 - fixed
  integer,dimension(:),allocatable,public,save          :: list_rigid_body_atoms
  character(len=8),dimension(:),allocatable,public,save :: atsymb ! atomic symbols
  integer,dimension(:),allocatable,public,save          :: atomic_numbers
  real (kind=dp),dimension(:,:),allocatable,public,save :: bmat_disconn_fragm ! B matrix for disconnected fragments

  !------------------------------ data storage (instead of PCHK etc. files)
  integer,dimension(:,:),allocatable,save              :: np_p_global ! for rationalizing coordinates
  real(kind=di_dp),dimension(:,:),allocatable,save     :: ut_global ! U matrix
  integer,dimension(:,:),allocatable,save              :: klist_global
  real(kind=di_dp),dimension(:),allocatable,save       :: xprim_global
  real(kind=di_dp),dimension(:),allocatable,save       :: savtor_global
  integer,dimension(:),allocatable,save                :: ktyp_global
  integer,dimension(:),allocatable,save                :: map_fixed_coords_global ! maps fixed internal to Cartesian coordinate
  real(kind=di_dp),dimension(:,:),allocatable,save     :: bmat_disconn_fragm_global
  integer,dimension(:,:),allocatable,save              :: ictyp_global
  integer,dimension(:,:),allocatable,save              :: icc_global
  real(kind=di_dp),dimension(:,:),allocatable,save     :: rcon_global
  real(kind=di_dp),dimension(:),allocatable,save       :: disp_DI_global
  real(kind=di_dp),dimension(:),allocatable,save       :: gradient_DI_global
  real(kind=di_dp),dimension(:,:),allocatable,save     :: hess_DI_global
  real(kind=di_dp),dimension(:,:),allocatable,save     :: hpad_global
  real(kind=di_dp),dimension(:,:),allocatable,save     :: xstr_global
  real(kind=di_dp),dimension(:,:),allocatable,save     :: gstr_global
  integer,dimension(:,:),allocatable,save              :: list_fixed_coords_global ! 0 - coordinate active, 1 - fixed
  integer,dimension(:),allocatable,save                :: list_rigid_body_atoms_global
  real(kind=di_dp),dimension(:),allocatable,save       :: vmdel_global
  real(kind=dp),dimension(:,:),allocatable,save        :: hess_cartesian_global ! backup storage

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  real(kind=dp),parameter :: one=1.0_dp
  real(kind=dp),parameter :: zero=0.0_dp

  !VM: temporary location
  logical :: fix_all_ions=.false.

contains

  subroutine deloc_utils_initialize(mdl,iprint)
    !=========================================================================!
    ! Initialize the module                                                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) mdl (in) - current model to be optimized                             !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   model                                                                 !
    !   parameters                                                            !
    !   cell                                                                  !
    !   io                                                                    !
    !   comms                                                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 27/03/2003                              !
    !=========================================================================!
    use rundat, only    : geom_force_tol, geom_disp_tol, geom_energy_tol
    use constants, only : stdout

    implicit none
    !-------------------------------------------------------------------------!
    ! Arguments
    type(castep_model), intent(inout), target :: mdl
    integer, intent(in)                       :: iprint
    !-------------------------------------------------------------------------!
    integer :: ierr,ntrans1
    !-------------------------------------------------------------------------!

    !point mdl_ptr to structure
    mdl_ptr => mdl

    !Set parallel output flag as appropriate
    di_on_root=pub_on_root

    !define output unit and verbosity level
    di_stdout = stdout 
    di_iprint = iprint
    use_deloc_int_output = .false.

    ! Store labels
    di_energy_label     = 'Ha'         ! call io_unit_label(energy_unit,di_energy_label)
    di_force_label      = 'Ha/Bohr'    ! call io_unit_label(force_unit,di_force_label)
    di_pressure_label   = 'Ha/Bohr**3' ! call io_unit_label(pressure_unit,di_pressure_label)
    di_length_label     = 'Bohr'       ! call io_unit_label(length_unit,di_length_label)
    !di_energy_label = trim(di_energy_label)
    !di_force_label = trim(di_force_label)
    !di_pressure_label = trim(di_pressure_label)
    !di_length_label = trim(di_length_label)

    ! no need for unit conversion, since ONETEP operates only in atomic units
    if (di_on_root) then
       length_conv = 1.0_dp ! io_atomic_to_unit(1.0_dp,length_unit)
       force_conv  = 1.0_dp ! io_atomic_to_unit(1.0_dp,force_unit)
       energy_conv = 1.0_dp ! io_atomic_to_unit(1.0_dp,energy_unit)
       bohr_di     = ANGSTROM ! io_unit_to_atomic(1.0_dp,"ang")
    endif
    call comms_bcast(pub_root_node_id,length_conv)
    call comms_bcast(pub_root_node_id,force_conv)
    call comms_bcast(pub_root_node_id,energy_conv)
    call comms_bcast(pub_root_node_id,bohr_di)

    number_atoms = mdl%cell%num_ions

    ! Set parameters for delocalized module
    period = .true.
    ldokpt = .true.
    tsflag = .false.
    l_scan_pes = .false.
    ldoqst = .false.
    ldoneb = .false.
    ip_error = 0
    lattice_dimension = 3
    maxdiis = 0
    ncons = 0
    ncons_uns = 0
    ioptc = 2
    tolg = geom_force_tol
    told = geom_disp_tol
    tole = geom_energy_tol
    iupdat = 2
    global_data_exists = .false.

    ! allocate arrays that are not going to change size during the optimization
    allocate (np_p(1:3,1:number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating np_p in deloc_utils_initialize')
    np_p = 0

    allocate(gradient_cart(1:3,1:number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating gradient_cart in deloc_utils_initialize')
    gradient_cart = zero

    allocate(coords_cart(1:3,1:number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating coords_cart in deloc_utils_initialize')
    coords_cart = zero

    allocate(hess_cartesian(1:3*number_atoms,1:3*number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating hess_cartesian in deloc_utils_initialize')
    hess_cartesian = zero

    ! vm: no symmetry in ONETEP
    !ntrans1 = max(1,num_symmetry_operations)
    ntrans1 = 1

    allocate(trans(1:3,1:3,1:ntrans1),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating trans in deloc_utils_initialize')
    trans = zero

    allocate(neqatm(1:number_atoms,1:ntrans1),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating neqatm in deloc_utils_initialize')
    neqatm = 0

    num_fixed_coords = 0
    num_fixed_constr_rigid_body = 0
    allocate(list_fixed_coords(1:3,1:number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating list_fixed_coords in deloc_utils_initialize')
    list_fixed_coords = 0

    num_rigid_body_atoms = 0
    allocate(list_rigid_body_atoms(1:number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating list_rigid_body_atoms in deloc_utils_initialize')
    list_rigid_body_atoms = 0

    ! atomic symbols
    allocate(atsymb(1:number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating atsymb in deloc_utils_initialize')

    allocate(atomic_numbers(1:number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating atomic_numbers in deloc_utils_initialize')

    return
  end subroutine deloc_utils_initialize

  subroutine deloc_utils_update_params
    !=========================================================================!
    ! Reinitialize the module if parameters file has changed                  !
    !-------------------------------------------------------------------------!
    ! Arguments:  none                                                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   io                                                                    !
    !   parameters                                                            !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 27/03/2003                              !
    !=========================================================================!
    implicit none
    !-------------------------------------------------------------------------!

    ! Dummy routine for now

    return
  end subroutine deloc_utils_update_params

  subroutine deloc_utils_deallocate
    !=========================================================================!
    ! Deallocate the module's memory                                          !
    !-------------------------------------------------------------------------!
    ! Arguments: none                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !  io                                                                     !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 27/03/2003                              !
    !=========================================================================!
    implicit none
    !-------------------------------------------------------------------------!
    integer :: ierr
    !-------------------------------------------------------------------------!

    if (allocated(np_p)) then
       deallocate (np_p,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating np_p in deloc_utils_deallocate')
    endif

    if (allocated(gradient_cart)) then
       deallocate (gradient_cart,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating gradient_cart in deloc_utils_deallocate')
    endif

    if (allocated(coords_cart)) then
       deallocate (coords_cart,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating coords_cart in deloc_utils_deallocate')
    endif

    if (allocated(hess_cartesian)) then
       deallocate (hess_cartesian,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hess_cartesian in deloc_utils_deallocate')
    endif

    if (allocated(hess_cartesian_global)) then
       deallocate (hess_cartesian_global,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hess_cartesian_global in deloc_utils_deallocate')
    endif

    if (allocated(trans)) then
       deallocate (trans,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating trans in deloc_utils_deallocate')
    endif

    if (allocated(neqatm)) then
       deallocate (neqatm,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating neqatm in deloc_utils_deallocate')
    endif

    if (allocated(atsymb)) then
       deallocate (atsymb,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating atsymb in deloc_utils_deallocate')
    endif

    if (allocated(atomic_numbers)) then
       deallocate (atomic_numbers,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating atomic_numbers in deloc_utils_deallocate')
    endif

    if (allocated(list_fixed_coords)) then
       deallocate(list_fixed_coords,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating list_fixed_coords in deloc_utils_deallocate')
    endif

    if (allocated(list_rigid_body_atoms)) then
       deallocate(list_rigid_body_atoms,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating list_rigid_body_atoms in deloc_utils_deallocate')
    endif

    if (allocated(bmat_disconn_fragm)) then
        deallocate(bmat_disconn_fragm,stat=ierr)
        if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bmat_disconn_fragm in deloc_utils_deallocate')
    endif

    return
  end subroutine deloc_utils_deallocate

  subroutine deloc_utils_nullify
    !=========================================================================!
    ! Nullify the module's pointers                                           !
    !-------------------------------------------------------------------------!
    ! Arguments: none                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 27/03/2003                              !
    !=========================================================================!
    implicit none
    !-------------------------------------------------------------------------!

    ! nullify the pointer
    nullify(mdl_ptr)

    return
  end subroutine deloc_utils_nullify



  subroutine deloc_utils_energy_gradient
    !=========================================================================!
    ! Returns the energy and forces on atoms for the current structure        !
    ! and stores them in EC and gradient_cart variables                       !
    !                                                                         !
    ! DI module does not need internal re-evaluation of energy/gradients,     !
    ! so this routine is needed at the beginning of each DI step.             !
    !-------------------------------------------------------------------------!
    ! Arguments: none                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: Model should be at ground state and have forces   !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 27/03/2003                              !
    !=========================================================================!

    implicit none
    !-------------------------------------------------------------------------!

    ! Local variables

    integer :: i0
    integer :: ispec
    integer :: iatom
    !-------------------------------------------------------------------------!

    ! Total energy from the model
    ec = mdl_ptr%total_energy

    ! get the gradients
    if (associated(mdl_ptr%forces)) then
       i0 = 0
       do ispec=1,mdl_ptr%cell%num_species
          do iatom=1,mdl_ptr%cell%num_ions_in_species(ispec)
             i0 = i0 + 1
             gradient_cart(1:3,i0) = - mdl_ptr%forces(1:3,iatom,ispec) ! note sign flip
          end do
       end do
    end if

    return
  end subroutine deloc_utils_energy_gradient



!------------------------------------------------------------------------------
  subroutine deloc_utils_mdl_to_internals(mdl,elements)
    !=========================================================================!
    ! Fill data structures for delocalized internals optimization based on    !
    ! the model data.                                                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mdl (inout) : The model to be used.                                   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   model                                                                 !
    !   cell                                                                  !
    !   io                                                                    !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) elat                                                                 !
    ! 2) elatb                                                                !
    ! 3) elatp                                                                !
    ! 4) num_dummies                                                          !
    ! 5) ipvlt                                                                !
    ! 6) ntrans                                                               !
    ! 7) symflag                                                              !
    ! 8) trans                                                                !
    ! 9) neqatm                                                               !
    !10) atsymb                                                               !
    !11) atomic_numbers                                                       !
    !12) coords_cart                                                          !
    !13) num_fixed_coords                                                     !
    !14) list_fixed_coords                                                    !
    !15) xc_redressed                                                         !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    use simulation_cell, only: castep_model, pub_cell, castep_cell_frac_to_cart
    use ion, only: element
    implicit none
    !-------------------------------------------------------------------------!
    type(castep_model), intent(in)        :: mdl
    type(ELEMENT), intent(in)             :: elements(pub_cell%nat)
    !-------------------------------------------------------------------------!
    integer :: ispec,iatom,indx,ierr,i,j
    real(kind=dp) :: constr(1:3)
    !-------------------------------------------------------------------------!

    ! dummy atoms
    num_dummies = 0

    ! cell data
    elat  = transpose(mdl%cell%real_lattice)
    elatb = transpose(mdl%cell%recip_lattice)
    elatp = elat

    ! factor elatp matrix using gaussian elimination
    call dgetrf(3,3,elatp,3,ipvlt,ierr)
    do i=1,3
       do j=1,3
          elatb(j,i) = zero
       enddo
       elatb(i,i) = one
       call dgetrs('N',3,1,elatp,3,ipvlt,elatb(1,i),3,ierr)
    enddo
    
    ! symmetry data 
    ! vm - no symmetry in ONETEP
    !ntrans = num_symmetry_operations
    ntrans = 1

    !symflag = num_symmetry_operations>1
    symflag = ntrans>1

    indx = 1
    do ispec=1,mdl%cell%num_species
       do iatom=1,mdl%cell%num_ions_in_species(ispec)
          atsymb(indx) = mdl%cell%species_symbol(ispec)

          do j=1,size(periodic_table_name)
             ! Get the atomic_number for the symbol
             if (atsymb(indx)==periodic_table_name(j)) atomic_numbers(indx)=j
          enddo

          indx = indx + 1
       enddo
    enddo

    ! coordinates (Cartesian)

    indx = 1
    do ispec = 1,mdl%cell%num_species
       do iatom = 1,mdl%cell%num_ions_in_species(ispec)
          call castep_cell_frac_to_cart(mdl%cell,mdl%cell%ionic_positions(:,iatom,ispec),coords_cart(:,indx))
          indx = indx + 1
       enddo
    enddo

    ! find number of fixed cartesian coordinates and make a list of them
    num_fixed_coords = 0
    do iatom=1,pub_cell%nat
       constr(1:3) = elements(iatom)%ion_constraint(1:3)
       select case (elements(iatom)%ion_constraint_type)
       case ('NONE')  ; continue
       ! Temporary - the formalism for now can handle only fixed cartesians, so look for the special cases
       case ('LINE')  
          if (constr(1) == zero .and. constr(2) == zero .and. constr(3) == one) then ! along 001: X and Y are fixed
              num_fixed_coords = num_fixed_coords + 2
              list_fixed_coords(1,iatom) = 1             
              list_fixed_coords(2,iatom) = 1             
          else if (constr(1) == zero .and. constr(2) == one .and. constr(3) == zero) then ! along 010: X and Z are fixed
              num_fixed_coords = num_fixed_coords + 2
              list_fixed_coords(1,iatom) = 1             
              list_fixed_coords(3,iatom) = 1             
          else if (constr(1) == one .and. constr(2) == zero .and. constr(3) == zero) then ! along 100: Y and Z are fixed
              num_fixed_coords = num_fixed_coords + 2
              list_fixed_coords(2,iatom) = 1             
              list_fixed_coords(3,iatom) = 1             
          end if
       case ('PLANE') 
          if (constr(1) == one .and. constr(2) == zero .and. constr(3) == zero) then ! normal 100: X is fixed
              num_fixed_coords = num_fixed_coords + 1
              list_fixed_coords(1,iatom) = 1             
          else if (constr(1) == zero .and. constr(2) == one .and. constr(3) == zero) then ! normal 010: Y is fixed
              num_fixed_coords = num_fixed_coords + 1
              list_fixed_coords(2,iatom) = 1             
          else if (constr(1) == zero .and. constr(2) == zero .and. constr(3) == one) then ! normal 001: Z is fixed
              num_fixed_coords = num_fixed_coords + 1
              list_fixed_coords(3,iatom) = 1             
          end if
       case ('FIXED') 
          num_fixed_coords = num_fixed_coords + 3
          do i = 1,3
              list_fixed_coords(i,iatom) = 1
          enddo
       case default
          if (di_on_root) then
             write(di_stdout,'(a,i6,a)') 'Error in deloc_utils_mdl_to_internals: &
                  &illegal value for elements(',iatom,')%ion_constraint_type'
             call comms_abort
          endif
       end select
    enddo

 
    ! Set the state of deloc_utils_rationalize_coords
    xc_redressed = .false.

    return
  end subroutine deloc_utils_mdl_to_internals

!------------------------------------------------------------------------------
  subroutine deloc_utils_internals_to_mdl(mdl)
    !=========================================================================!
    ! Fill the model data based on the data structures for delocalized        !
    ! internals optimization.                                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mdl (inout) : The model to be updated.                                !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   model                                                                 !
    !   cell                                                                  !
    !   io                                                                    !
    !   basis                                                                 !
    !   parameters                                                            !
    !   density                                                               !
    !   dm                                                                    !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) xc_redressed                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 27/03/2003                              !
    !=========================================================================!
    implicit none
    !-------------------------------------------------------------------------!
    type(castep_model), intent(inout) :: mdl
    !-------------------------------------------------------------------------!
    integer :: ispec,iatom,indx,ierr,ndim,i,j
    real(kind=dp),dimension(:,:), allocatable :: hessian
    !-------------------------------------------------------------------------!

    ! Return atoms to their original cells, too (and restore afterwards)
    if(di_on_root)then
        if (xc_redressed) call deloc_utils_rationalize_coords(-1)

        indx = 1
        do ispec = 1,mdl%cell%num_species
           do iatom = 1,mdl%cell%num_ions_in_species(ispec)
              call castep_cell_cart_to_frac(mdl%cell,coords_cart(:,indx),mdl%cell%ionic_positions(:,iatom,ispec))
              indx = indx + 1
           enddo
        enddo

        ! Restore coords_cart after the call to redres
        if (.not.xc_redressed) call deloc_utils_rationalize_coords(2)
    end if
    call comms_bcast(pub_root_node_id,mdl%cell%ionic_positions,3*mdl%cell%num_species)
    
    !... save hessian on the model
    ndim = 3*number_atoms
    allocate(hessian(1:ndim,1:ndim),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating hessian in deloc_utils_internals_to_mdl')
    hessian = zero

    if(di_on_root)then
       hessian = hess_cartesian

       do i = 1,ndim
          do j = 1,ndim
             mdl%bfgs_inv_Hessian(i,j) = hessian(i,j)
          enddo
       enddo
    end if
       
    !... and set the appropriate flags in model
    call castep_model_cell_changed(mdl)
 
    deallocate(hessian,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hessian in deloc_utils_internals_to_mdl')

    return
  end subroutine deloc_utils_internals_to_mdl

!------------------------------------------------------------------------------

  subroutine deloc_utils_rationalize_coords(nway)
    !=========================================================================!
    ! Redress atoms into unit cell [0,1]                                      !
    ! The actions of this routine is as follows:                              !
    !                                                                         !
    ! nway==1, xc_redressed==.false.                                          !
    !   - assumes that atoms are in arbitrary cells, calculates np_p set of   !
    !     translations, brings atoms into [0,1], sets xc_redressed=.true.     !
    !                                                                         !
    ! nway==2, xc_redressed==.false.                                          !
    !   - assumes that atoms are in arbitrary cells, uses stored np_p set of  !
    !     translations, brings atoms into [0,1], sets xc_redressed=.true.     !
    !                                                                         !
    ! nway==-1, xc_redressed==.true.                                          !
    !   - assumes that atoms are in [0,1] cell, uses stored np_p set of       !
    !     translations, brings atoms back to the original cells,              !
    !     sets xc_redressed=.false.                                           !
    !                                                                         !
    ! There is no action for any other combination of nway and xc_redressed   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) nway (in) : defines the action of the routine                        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) lattice_dimension                                                    !
    ! 3) number_atoms                                                         !
    ! 4) elat                                                                 !
    ! 5) elatp                                                                !
    ! 6) ipvlt                                                                !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) np_p                                                                 !
    ! 2) coords_cart                                                          !
    ! 3) xc_redressed                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 27/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)        :: nway
    !-------------------------------------------------------------------------!
    real(kind=dp),dimension(1:3) :: prj,prma_p,prmi_p
    logical :: amov
    integer :: i,j,np,ierr
    real(kind=dp),parameter :: eps_edg=1.0e-14_dp
    !-------------------------------------------------------------------------!

    amov = .false.

    ! redres using new np_p just developed
    if (nway==1 .and. .not.xc_redressed) then

       ! start

       do i=1,lattice_dimension
          ! set borders of cell
          prma_p(i) = one
          prmi_p(i) = zero
       enddo

       do i=1,number_atoms
          prj(1) = coords_cart(1,i)
          prj(2) = coords_cart(2,i)
          prj(3) = coords_cart(3,i)
          call dgetrs('N',3,1,elatp,3,ipvlt,prj,3,ierr)

          do j=1,lattice_dimension
             np_p(j,i) = 0  
             if ((prj(j)-prma_p(j))>eps_edg) then
                amov=.true.
                np = prj(j) - prma_p(j) + one
                np_p(j,i) = np 
                coords_cart(1,i) = coords_cart(1,i) - np*elat(1,j)
                coords_cart(2,i) = coords_cart(2,i) - np*elat(2,j)
                coords_cart(3,i) = coords_cart(3,i) - np*elat(3,j)
                if (di_on_root .and. di_iprint>5) write(di_stdout,'(a,3i5,3f10.5,1p,e10.2)')' transl -', &
     &             i,j,np,coords_cart(1,i),coords_cart(2,i),coords_cart(3,i),prj(j)
                elseif ((prmi_p(j)-prj(j))>eps_edg) then
                   amov=.true.
                   np = - prj(j) + prmi_p(j) + one
                   np_p(j,i) = -np 
                   coords_cart(1,i) = coords_cart(1,i) + np*elat(1,j)
                   coords_cart(2,i) = coords_cart(2,i) + np*elat(2,j)
                   coords_cart(3,i) = coords_cart(3,i) + np*elat(3,j)
                   if (di_on_root .and. di_iprint>5) write(di_stdout,'(a,3i5,3f10.5,1p,e10.2)')' transl +', &
     &                i,j,np,coords_cart(1,i),coords_cart(2,i),coords_cart(3,i),prj(j)
                end if
             enddo ! j=1,lattice_dimension
          enddo ! i=1,number_atoms
 
          ! set the flag for the state of deloc_utils_rationalize_coords
          xc_redressed = .true.

          if (di_on_root .and. amov .and. di_iprint>1) write(di_stdout,'(/a/)') &
     &       'Note: One or more atoms were translated into the central cell'

       elseif (nway==2 .and. .not.xc_redressed) then
          ! redres using old np_p
          do i=1,number_atoms
             do j=1,lattice_dimension
                np = np_p(j,i) 
                coords_cart(1,i) = coords_cart(1,i) - np*elat(1,j)
                coords_cart(2,i) = coords_cart(2,i) - np*elat(2,j)
                coords_cart(3,i) = coords_cart(3,i) - np*elat(3,j)
                if (di_on_root .and. di_iprint>5) write(di_stdout,'(a,3i5,3f10.5)')' translat', &
     &             i,j,np,coords_cart(1,i),coords_cart(2,i),coords_cart(3,i)
             enddo ! j=1,lattice_dimension
          enddo ! i=1,number_atoms

          ! set the flag for the state of deloc_utils_rationalize_coords
          xc_redressed = .true.
 
       elseif (nway==-1 .and. xc_redressed) then

          ! redres back using old np_p
          do i=1,number_atoms
             do j=1,lattice_dimension
                np = -np_p(j,i) 
                coords_cart(1,i) = coords_cart(1,i) - np*elat(1,j)
                coords_cart(2,i) = coords_cart(2,i) - np*elat(2,j)
                coords_cart(3,i) = coords_cart(3,i) - np*elat(3,j)
                if (di_on_root .and. di_iprint>5) write(di_stdout,'(a,3i5,3f10.5)')' trans back', &
     &             i,j,np,coords_cart(1,i),coords_cart(2,i),coords_cart(3,i)
             enddo ! j=1,lattice_dimension
          enddo ! i=1,number_atoms
 
          ! set the flag for the state of deloc_utils_rationalize_coords
          xc_redressed = .false.

       endif

    return
  end subroutine deloc_utils_rationalize_coords



  subroutine deloc_utils_io_abort(message)
    !=========================================================================!
    ! Called on abnormal termination of the program, ensures completion of    !
    ! all outstanding I/O operations (wrapper to io_abort)                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) message : input : error message written on abort                     !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    ! 1) io : io_abort : aborts execution                                     !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  none                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 27/03/2003                              !
    !=========================================================================!

    implicit none

    character(len=*), intent(in) :: message
    !-------------------------------------------------------------------------!

    if (di_on_root) then
       write (di_stdout,'(a)') message
    end if
    call comms_abort

    return
  end subroutine deloc_utils_io_abort

  subroutine deloc_utils_write_trajectory(energy,denergy,gmax,step,iteration)
    !=========================================================================!
    ! Writes latest SCF iteration data to trajectory type file.               !
    ! Dummy routine in CASTEP since the writer used is in geometry module.    !
    ! We cannot put a wrapper here since it would have introduced a           !
    ! circular module dependency.                                             !
    !                                                                         !
    ! Other applications that want to save trajectory-style information       !
    ! have to provide a real routine here.                                    !
    !                                                                         !
    !                             [opttbl_iter in DMol]                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) energy (input)    = total energy (Ha)                                !
    ! 2) denergy (input)   = energy change (Ha)                               !
    ! 3) gmax (input)      = max grad (Ha/Bohr)                               !
    ! 4) step (input)      = max step (Bohr)                                  !
    ! 5) iteration (input) = iteration #                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  none                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 27/03/2003                              !
    !=========================================================================!
    use constants, only: stdout
    implicit none
    real(kind=dp),intent(in)     :: energy
    real(kind=dp),intent(in)     :: denergy
    real(kind=dp),intent(in)     :: gmax
    real(kind=dp),intent(in)     :: step
    integer,intent(in)           :: iteration
    !-------------------------------------------------------------------------!

    if (di_on_root) then
       ! some kind of i/o
       ! qoh: Do some stuff with inputs to avoid compiler warnings
       if (max(energy,denergy,gmax,step,real(iteration,kind=DP)) == 0.5_DP) &
          write (stdout,*) "."
    endif
    return

  end subroutine deloc_utils_write_trajectory

  subroutine deloc_utils_output_converged(iteration,enthalpy,dE,Fmax,dRmax, &
       & converged_dE,converged_Fmax,converged_dRmax)
    !=========================================================================!
    ! Output relevant convergence data to stdout                              !
    ! (based on geom_output_converged)                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v1.0, 11/03/2003                              !
    !=========================================================================!
    use constants, only: DP, stdout
    use services, only : services_flush
    implicit none

    integer       , intent(in) :: iteration
    real (kind=dp), intent(in) :: enthalpy

    real (kind=dp), intent(in) :: dE
    real (kind=dp), intent(in) :: Fmax   !max Force
    real (kind=dp), intent(in) :: dRmax  !max displacement
    !real (kind=dp), intent(in) :: Smax   !max Stress

    logical,        intent(in) :: converged_dE
    logical,        intent(in) :: converged_Fmax
    logical,        intent(in) :: converged_dRmax
    !logical,        intent(in) :: converged_Smax

    !local vars
    character(len=80) :: data_string
    character(len=80) :: divider_string
    character(len=80) :: label_string
    integer           :: string_index
    integer           :: len_label, len_value, len_tol, len_unit, len_flag

#ifdef debug
    write(stdout,*) 'starting deloc_utils_output_converged'
#endif

    !initialise strings
    data_string     = repeat(' ',len(data_string))
    divider_string  = repeat(' ',len(divider_string))
    label_string    = repeat(' ',len(label_string))
    !                  1234567890123 23456789012345678 23456789012345678 234567890123 23456
    divider_string  = '+-----------+-----------------+-----------------+------------+-----+'
    label_string    = '| Parameter |      value      |    tolerance    |    units   | OK? |'

    !write out all the relevant data (root node only)
    if (di_on_root) then

       !explain what we are testing ...
       write (di_stdout,1) 'DI: finished iteration',iteration,' with enthalpy=', &
            & enthalpy,trim(di_energy_label)

       !header
       write(di_stdout,*) ' '
       write(di_stdout,7) divider_string
       write(di_stdout,7) label_string
       write(di_stdout,7) divider_string

       !dE                                                       '1234567890'
       string_index=1
       len_label   =13
       len_value   =18
       len_tol     =18
       len_unit    =13
       len_flag    = 6
       write(data_string(string_index:string_index+len_label),2) '  dE/ion  '
       string_index=string_index+len_label
       write(data_string(string_index:string_index+len_value),3) dE
       string_index=string_index+len_value
       write(data_string(string_index:string_index+len_tol),  4) tole
       string_index=string_index+len_tol
       write(data_string(string_index:string_index+len_unit), 5) trim(di_energy_label)
       string_index=string_index+len_unit
       if (converged_dE) then
          write(data_string(string_index:string_index+len_flag), 6) 'Yes'
       else
          write(data_string(string_index:string_index+len_flag), 6) 'No '
       end if
       string_index=string_index+len_flag
       !cross check all is well (NB strings start from 1 not 0 ...)
       if (string_index>71) write (di_stdout,*) 'geom_output_converged: format problem?'
       write(di_stdout,7) data_string

       !|F|max                                                      '1234567890'
       if (.not.fix_all_ions) then
          string_index=1
          write(data_string(string_index:string_index+len_label),2) '  |F|max  '
          string_index=string_index+len_label
          write(data_string(string_index:string_index+len_value),3) Fmax
          string_index=string_index+len_value
          write(data_string(string_index:string_index+len_tol),  4) tolg
          string_index=string_index+len_tol
          write(data_string(string_index:string_index+len_unit), 5) trim(di_force_label)
          string_index=string_index+len_unit
          if (converged_Fmax) then
             write(data_string(string_index:string_index+len_flag), 6) 'Yes'
          else
             write(data_string(string_index:string_index+len_flag), 6) 'No '
          end if
          write(di_stdout,7) data_string
       end if

       !|dR|max                                                     '1234567890'
       if (.not.fix_all_ions) then
          string_index=1
          write(data_string(string_index:string_index+len_label),2) '  |dR|max '
          string_index=string_index+len_label
          write(data_string(string_index:string_index+len_value),3) dRmax
          string_index=string_index+len_value
          write(data_string(string_index:string_index+len_tol),  4) told
          string_index=string_index+len_tol
          write(data_string(string_index:string_index+len_unit), 5) trim(di_length_label)
          string_index=string_index+len_unit
          if (converged_dRmax) then
             write(data_string(string_index:string_index+len_flag), 6) 'Yes'
          else
             write(data_string(string_index:string_index+len_flag), 6) 'No '
          end if
          write(di_stdout,7) data_string
       end if
       
       !footer
       write(di_stdout,7) divider_string
       write(di_stdout,*) ' '

       !back to all nodes
    end if

    call services_flush

1   format(1x,a,i4,1x,a,es17.8e3,1x,a)

2   format('|',a10,    1x,'|')       !1+10+1+1=13
3   format(1x,es15.6e3,1x,'|')       !1+15+1+1=18
4   format(1x,es15.6e3,1x,'|')       !1+15+1+1=18
!NB In IO the phys_unit is defined as 30 characters long but max. is 10 in practice
5   format(1x,a10,     1x,'|')       !1+10+1+1=13
6   format(1x,a3,      1x,'|')       !1+3+1+1=6
                                     !13+18+18+13+6=68
7   format(1x,a68,' <-- DI')

    return
  end subroutine deloc_utils_output_converged

!--------------------------------------------------------
  subroutine deloc_utils_geom_converged(coords_cart0,ndeg,neg,negreq,gradient_DI,cnvgd)
    !=========================================================================!
    ! check convergence of cartesian gradients and displacement               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) coords_cart0 (input)    = old cartesians                             !
    ! 2) ndeg (input)   = number of degrees of freedom                        !
    ! 3) neg (input)    = # of negative hessian eigenvalues                   !
    ! 4) negreq (input) = required # of negative hessian eigenvalues          !
    ! 5) gradient_DI (input) = derivatives in primitive space                 !
    ! 6) cnvgd (output) = are we converged?                                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) use_deloc_int_output                                                 !
    ! 3) ec                                                                   !
    ! 4) ncycle                                                               !
    ! 5) number_atoms                                                         !
    ! 6) num_fixed_coords                                                     !
    ! 7) map_fixed_coords                                                     !
    ! 8) coords_cart                                                          !
    ! 9) gradient_cart                                                        !
    !10) atsymb                                                               !
    !11) tole                                                                 !
    !12) tolg                                                                 !
    !13) told                                                                 !
    !14) eold                                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  none                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 27/03/2003                              !
    !=========================================================================!
    implicit none

    real(kind=dp),dimension(:,:),intent(in)       :: coords_cart0
    integer,intent(in)                            :: ndeg
    integer,intent(in)                            :: neg
    integer,intent(in)                            :: negreq
    real(kind=dp),dimension(:),intent(in)         :: gradient_DI
    logical,intent(out)                           :: cnvgd
    !-------------------------------------------------------------------------!
    character(len=3) :: dcnvgd,gcnvgd,ecnvgd
    character(len=1),dimension(1:3) :: c
    data c(1)/'x'/, c(2)/'y'/, c(3)/'z'/
    real(kind=dp) :: gx,dx,gmax,disp_max,cx,cn,edif,tmin
    real(kind=dp),dimension(1:3) :: disp_centre_mass
    integer :: iatm,ii,j,it,i,ix
    logical :: skip
    !-------------------------------------------------------------------------!
    !  Checks convergence during geometry optimization
    !  in Cartesian coordinates
    !  The flag Cnvgd is returned .true. if the maximum component
    !  of the gradient vector, gradient_cart, is below TolG AND either the
    !  maximum component of the displacement vector, D, is below
    !  TolD OR the energy change is below TolE
    !  (TolG, TolD and TolE are currently in atomic units)
    !  Additionally, of course, the Hessian must have the
    !  correct eigenvalue structure
    !
    ! if ncons substitute cartesian forces by internals
    ! for the final test of convergence

    if (di_iprint>2 .and. di_on_root .and. use_deloc_int_output) write(di_stdout,1000)

    gmax = zero
    disp_max = zero 

    disp_centre_mass = zero

    do iatm=1,number_atoms
       do j=1,3
          disp_centre_mass(j) = disp_centre_mass(j) + coords_cart(j,iatm) - coords_cart0(j,iatm)
       enddo
    enddo
    do j=1,3
       disp_centre_mass(j) = disp_centre_mass(j) / number_atoms
    enddo

    do iatm=1,number_atoms
       ii = 3*(iatm-1)
       gx = 0.0_dp
       dx = 0.0_dp
       do j=1,3
          it = ii+j
 
          ! check if it is fixed do not consider
          skip = .false.
          if (num_fixed_coords/=0) then
             do ix=1,num_fixed_coords
                if (map_fixed_coords(ix)==it) skip=.true.
             enddo
          endif

          if (.not.skip) then
             gx = gx + gradient_cart(j,iatm)**2
             dx = dx + (coords_cart(j,iatm)-coords_cart0(j,iatm)-disp_centre_mass(j))**2
          endif
   
          if (di_iprint>2 .and. di_on_root .and. use_deloc_int_output) then
             cx = coords_cart(j,iatm)
             cn = cx + dx            
             if (j==1) then
                write(di_stdout,1010) iatm,atsymb(iatm),c(j),cn,gx,dx,cx
             else 
                write(di_stdout,1020)                   c(j),cn,gx,dx,cx
             endif
          endif
       enddo ! j=1,3

       gmax = max(gmax,gx)
       disp_max = max(disp_max,dx)
    enddo ! iatm=1,number_atoms

    gmax = sqrt(gmax)
    disp_max = sqrt(disp_max)

    ! ncons part
    gx = zero
    if (ncons>0) then
       ! find maximum of internal force
       gx = zero
       do i=1,ndeg
          gx = max(gx,abs(gradient_DI(i)))
       enddo
       ! substitute!!
       gmax = gx
    endif

    edif = ec - eold

    if (di_iprint>2 .and. di_on_root .and. use_deloc_int_output) then
       gcnvgd = ' NO'
       dcnvgd = ' NO'
       ecnvgd = ' NO'
       if (disp_max<told) dcnvgd = 'YES'
       if (abs(edif)<tole) ecnvgd = 'YES'
       write(di_stdout,1030)

       ! ncons
       if (ncons/=0) then
          if (gx<tolg) gcnvgd = 'YES'
          write(di_stdout,1045) gx,tolg,gcnvgd
       else
          if (gmax<tolg) gcnvgd = 'YES'
          write(di_stdout,1040) gmax,tolg,gcnvgd
       endif

       write(di_stdout,1050) disp_max,told,dcnvgd
       if (eold==zero) then
          write(di_stdout,1059) tole,ecnvgd
       else
          write(di_stdout,1060) edif,tole,ecnvgd
       endif
    endif

    ! no print for NEB

    if (di_on_root .and. use_deloc_int_output) then
       if (ncycle==1) then
          write(di_stdout,1200) tole,tolg,told,ncycle,ec,zero,gmax,disp_max
          if (.not.ldoneb)  call deloc_utils_write_trajectory(ec,ec,gmax,disp_max,ncycle)
       else
          write(di_stdout,1100) ncycle,ec,edif,gmax,disp_max
          if (.not.ldoneb)  call deloc_utils_write_trajectory(ec,edif,gmax,disp_max,ncycle)
       endif
    endif


    !  Converged?

    cnvgd = abs(edif)<tole .and. (disp_max<told .or. gmax<tolg) .and. neg==negreq

    if (ncycle==1) then
       tmin = 0.1_dp
       cnvgd = gmax<tolg*tmin .and. disp_max<told*tmin .and. neg==negreq
    endif

    return

 1000 format(//,21X,' Coordinates and Displacements in Atomic Units',/, &
     &       '   ATOM            Current Value    Gradient  Displacement   New Value')
 1010 format(I3,2X,A8,2X,A1,4(2x,f12.6))
 1020 format(15X,A1,4(2x,f12.6))
 1030 format(/,29X,'Maximum     Tolerance    Cnvgd?')
 1040 format(9X,'Gradient           ',F8.5,6X,F8.5,5X,A3)
 1045 format(9X,'Internal Gradient  ',F8.5,6X,F8.5,5X,A3)
 1050 format(9X,'Displacement       ',F8.5,6X,F8.5,5X,A3)
 1059 format(9X,'Energy change     ',(3X,'---',3X),6X,F8.5,5X,A3,/)
 1060 format(9X,'Energy change     ',F9.6,6X,F8.5,5X,A3,/)

 1100 format(/,7x,'Cycle',4x,'Total Energy',3x, &
     & 'Energy change   Max Gradient   Max Displacement', &
     &  /,     'opt==  ',i3,3x,f15.7,4x,f11.7,7x,f9.6,6x,f9.6,/)
 1200 format(/,'opt==',2x,'Cycle',4x,'Total Energy',3x, &
     & 'Energy change   Max Gradient   Max Displacement', &
     &  /,     'opt==  ',4x,'tolerance:.......',4x,f11.7,7x,f9.6,6x,f9.6, &
     &  /,     'opt==  ',i3,3x,f15.7,4x,f11.7,7x,f9.6,6x,f9.6,/)

  end subroutine deloc_utils_geom_converged

  subroutine deloc_utils_get_alt_bond_list(bonds,nbonds,ierr,connectivity)
    !=========================================================================!
    ! Read bond list from an alternative file (could be MDF)                  !
    ! bonds is really all we need; the actual getmdf can return lots of       !
    ! other stuff. all of that is classified as optional here.                !
    !                                                                         !
    !                       [getmdf in dmol]                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) bonds(:,2) (out) : List of bonds found.                              !
    !    bonds(i,1) - first atom in the bond                                  !
    !    bonds(i,2) - packed index of the 2nd atom (includes cell translation)!
    ! 2) nbonds (out) : Number of bonds found                                 !
    ! 3) ierr (out)                                                           !
    ! 4) connectivity(number_atoms,number_atoms) (out,optional)               !
    !    Matrix of atomic connectivity                                        !
    !    connectivity(i,j) = 1 : i and j are connected, 0 - not connected     !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:  none                                     !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    implicit none
    integer,dimension(:,:),intent(out)                    :: bonds
    integer,intent(out)                                   :: nbonds
    integer,intent(out)                                   :: ierr
    integer,dimension(:,:),intent(out),optional           :: connectivity
    !-------------------------------------------------------------------------!
    ! Local variables
    !-------------------------------------------------------------------------!
    ierr = 0
    nbonds = 0
    bonds = 0
    if (present(connectivity)) connectivity = 0

    return
  end subroutine deloc_utils_get_alt_bond_list

  subroutine deloc_utils_save_structure
    !=========================================================================!
    ! Save structural data if necessary - CASTEP does not do anything here,   !
    ! DMol wants to save INCOOR file for some reason.                         !
    !-------------------------------------------------------------------------!
    ! Arguments: none                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    implicit none
    !-------------------------------------------------------------------------!
    ! Local variables
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!

    !  We already have atsymb, coords_cart and elatb (for 'cell') in the module,
    !  all that is needed is to write them to file
    !          call getatno(natoms,atsymb,zz(1))
    !          call wrcoord(natoms,atsymb,zz(1),coords_cart,period,cell)

    return
  end subroutine deloc_utils_save_structure


  subroutine deloc_utils_read_opt_mode(vmode,modeok)
    !=========================================================================!
    ! Read from file (or from memory?) the mode to be followed.               !
    ! Dummy in CASTEP, rdoptmode wrapper in DMol                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) vmode(out) - the mode (cartesian coordinates)                        !
    ! 2) modeok(out) - logical switch, did we get the mode?                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    implicit none
    real(kind=dp),dimension(:),intent(out)       :: vmode
    logical,intent(out)                          :: modeok
    !-------------------------------------------------------------------------!
    ! Local variables
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!

    modeok = .false.
    vmode = zero

    return
  end subroutine deloc_utils_read_opt_mode
  
  subroutine deloc_utils_read_DI_constraints(ncons,icc,rcon)
    !=========================================================================!
    ! Read from file (or from memory?) internal constraints                   !
    !                                                                         !
    !                              [rdcon_p in DMol]                          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ncons(in) - number of constraints                                    !
    ! 2) icc(out) - constraint description                                    !
    !    icc(19,   -  constraint type                                         !
    !               1 - fixed distance                                        !
    !               2 - fixed bond angle                                      !
    !               3 - fixed dihedral angle                                  !
    !  icc      -  atoms involved in constraint                               !
    !                ic1               atom constraint                        !
    !                ic1-ic2           distance constraint                    !
    !                ic1-ic2-ic3       angle constraint                       !
    !                ic1-ic2-ic3-ic4   dihedral constraint                    !
    ! 3) rcon(out)  -  constraint values (a.u.)                               !
    !   Rcon(.,1) is a coordinate                                             ! 
    !   Rcon(.,2) is a gradient                                               !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    implicit none
    integer,intent(in)                           :: ncons
    real(kind=dp),dimension(:,:),intent(out)     :: rcon
    integer,dimension(:,:),intent(out)           :: icc
    !-------------------------------------------------------------------------!
    ! Local variables
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!
    
    icc = ncons  ! qoh: use ncons to prevent compiler warning
    rcon = zero
    icc = 0

    return
  end subroutine deloc_utils_read_DI_constraints

  subroutine deloc_utils_read_fix_atoms
    !=========================================================================!
    ! Read from file (or from memory?) atoms with all internals               !
    ! completely fixed. This is obviously more strict than simply fixing      !
    ! the atom, and it is not currently supported in CASTEP input             !
    ! Dummy in CASTEP (this information cannot be deduced from CASTEP input)  !
    !                                                                         !
    !                         [rdfixat in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments: none                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:  none                                     !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) num_rigid_body_atoms                                                 !
    ! 2) list_rigid_body_atoms                                                !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   cell                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    implicit none
    !-------------------------------------------------------------------------!
    ! Local variables
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!
    ! Atoms with all fixed internals (we do not have this functionality)
    ! It's not enough to know that all three coordinates of the atom are
    ! fixed, it has to be tagged as a special constraint.

    num_rigid_body_atoms = 0
    list_rigid_body_atoms = 0
  !  iat = 1
  !  do nsp = 1,current_cell%num_species
  !     do ni = 1,current_cell%num_ions_in_species(nsp)
  !        if (list_fixed_coords(1,iat)==1 .and. list_fixed_coords(2,iat)==1 .and. list_fixed_coords(3,iat)==1) then
  !           num_rigid_body_atoms = num_rigid_body_atoms + 1
  !           list_rigid_body_atoms(num_rigid_body_atoms) = iat
  !        endif
  !        iat = iat + 1
  !     end do
  !  end do


    return
  end subroutine deloc_utils_read_fix_atoms

  subroutine deloc_utils_read_fix_cartesians
    !=========================================================================!
    ! Read from file (or from memory?) fixed cartesian coordinates            !
    ! Dummy in ONETEP (this information is already deduced from               !
    ! the deloc_utils_mdl_to_internals call in GEOMETRY module)               !
    !                          [rdfixp in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments: none                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) list_fixed_coords                                                    !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    implicit none

    !-------------------------------------------------------------------------!
    ! Local variables
    !-------------------------------------------------------------------------!
    !integer :: i,iat
    !-------------------------------------------------------------------------!

    return
  end subroutine deloc_utils_read_fix_cartesians

!------------------------------------------------------------------------------
  subroutine deloc_utils_read_hessian(ndim,hess,status)
    !=========================================================================!
    ! Read hessian matrix from memory (instead of HESSIAN file)               !
    !                                                                         !
    !                       [rdhessian in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndim (in) : Dimension of the hessian matrix                          !
    ! 2) hess (in) : Hessian matrix to be restored from memory                !
    ! 3) status (out) : Flag to tell whether we found hessian in memory       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) hess_cartesian_global                                                !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    implicit none
    integer,intent(in)                                    :: ndim
    real(kind=dp),dimension(1:ndim,1:ndim),intent(out)    :: hess
    integer,intent(out)                                   :: status
    !-------------------------------------------------------------------------!
    integer :: ierr
    !-------------------------------------------------------------------------!
    status = -1
    if (di_on_root) then
       status = -1 ! assume no hessian
       if (allocated(hess_cartesian_global)) then
          if (size(hess_cartesian_global,1)==ndim .and. size(hess_cartesian_global,2)==ndim) then
             hess = hess_cartesian_global
             status = 0 ! found hessian
          else 
             deallocate(hess_cartesian_global,stat=ierr) 
             if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hess_cartesian_global in deloc_utils_read_hessian')
          endif
       endif
    endif ! di_on_root

    return
  end subroutine deloc_utils_read_hessian


!------------------------------------------------------------------------------
  subroutine deloc_utils_write_hessian(ndim,hess)
    !=========================================================================!
    ! Saves hessian matrix in memory (instead of HESSIAN file)                !
    !                                                                         !
    !                       [wrhessian in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndim (in) : Dimension of the hessian matrix                          !
    ! 2) hess (in) : Hessian matrix to be stored                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) hess_cartesian_global                                                !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    implicit none
    integer,intent(in)                                   :: ndim
    real(kind=dp),dimension(1:ndim,1:ndim),intent(in)    :: hess
    !-------------------------------------------------------------------------!
    integer :: ierr
    !-------------------------------------------------------------------------!
    if (.not. di_on_root) return

    if (allocated(hess_cartesian_global)) then
       deallocate(hess_cartesian_global,stat=ierr) 
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hess_cartesian_global in deloc_utils_write_hessian')
    endif

    allocate(hess_cartesian_global(1:ndim,1:ndim),stat=ierr) 
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating hess_cartesian_global in deloc_utils_write_hessian')

    hess_cartesian_global = hess
    return
  end subroutine deloc_utils_write_hessian


!------------------------------------------------------------------------------

  subroutine deloc_utils_clean_backup
    !=========================================================================!
    ! Deallocate arrays that are used as a backup (used to be PCHK file)      !
    ! Conceivably some other action will be needed if other means of storage  !
    ! are used (e.g., close files)                                            !
    !-------------------------------------------------------------------------!
    ! Arguments:  None                                                        !
    !-------------------------------------------------------------------------!
    ! Modules used: None                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:   All _global arrays                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: Global arrays have to be allocated                !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 27/03/2003                              !
    !=========================================================================!
    implicit none

    !-------------------------------------------------------------------------!
    integer :: ierr
    !-------------------------------------------------------------------------!

    deallocate(np_p_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating np_p_global in deloc_deallocate_global')
    deallocate(ut_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ut_global in deloc_deallocate_global')
    deallocate(klist_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating klist_global in deloc_deallocate_global')
    deallocate(xprim_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating xprim_global in deloc_deallocate_global')
    deallocate(savtor_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating savtor_global in deloc_deallocate_global')
    deallocate(ktyp_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ktyp_global in deloc_deallocate_global')
    deallocate(map_fixed_coords_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating map_fixed_coords_global in deloc_deallocate_global')
    deallocate(bmat_disconn_fragm_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bmat_disconn_fragm_global in deloc_deallocate_global')
    deallocate(ictyp_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ictyp_global in deloc_deallocate_global')
    deallocate(icc_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icc_global in deloc_deallocate_global')
    deallocate(rcon_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating rcon_global in deloc_deallocate_global')
    deallocate(gradient_DI_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating gradient_DI_global in deloc_deallocate_global')
    deallocate(disp_DI_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating disp_DI_global in deloc_deallocate_global')
    deallocate(hess_DI_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hess_DI_global in deloc_deallocate_global')
    deallocate(xstr_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating xstr_global in deloc_deallocate_global')
    deallocate(gstr_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating gstr_global in deloc_deallocate_global')
    deallocate(list_fixed_coords_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating list_fixed_coords_global in deloc_deallocate_global')
    deallocate(list_rigid_body_atoms_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating list_rigid_body_atoms_global in deloc_deallocate_global')
    deallocate(vmdel_global,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating vmdel_global in deloc_deallocate_global')

    global_data_exists = .false.
    return

  end subroutine deloc_utils_clean_backup

!------------------------------------------------------------------------------
  
  subroutine deloc_utils_read_di_data (np_p_loc,ut_loc,klist_loc,xprim_loc,savtor_loc,ktyp_loc, &
     &     map_fixed_coords_loc,bmat_disconn_fragm_loc,ictyp_loc,icc_loc,rcon_loc, &
     &     disp_DI_loc,gradient_DI_loc,hess_DI_loc,hpad_loc,vmode_loc,xstr_loc,gstr_loc,list_fixed_coords_loc, &
     &     list_rigid_body_atoms_loc,num_disconn_fragm,nile,ierr)
    !=========================================================================!
    ! Restore a lot of saved data (used to be from PCHK file). If the         !
    ! corresponding global array exists and has nonzero dimensions, its       !
    ! contents is copied into the local (<>_loc) array and returned.          !
    ! Only root node does the reading                                         !
    !                                                                         !
    ! NOTE: local arrays are re-allocated if they are of the wrong size, or   !
    !       allocated if they are no initialize. Any wrapper that uses disk   !
    !       i/o here should be beware of this.                                !
    !                                                                         !
    !                           [rdpchk in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) np_p_loc (out)                                                       !
    ! 2) ut_loc (out)                                                         !
    ! 3) klist_loc (out)                                                      !
    ! 4) xprim_loc (out)                                                      !
    ! 5) savtor_loc (out)                                                     !
    ! 6) ktyp_loc (out)                                                       !
    ! 7) map_fixed_coords_loc (out)                                           !
    ! 8) bmat_disconn_fragm_loc (out)                                         !
    ! 9) ictyp_loc (out)                                                      !
    !10) icc_loc (out)                                                        !
    !11) rcon_loc (out)                                                       !
    !12) disp_DI_loc (out)                                                    !
    !13) gradient_DI_loc (out)                                                ! 
    !14) hess_DI_loc (out)                                                    !
    !15) hpad_loc (out)                                                       !
    !16) vmode_loc (out)                                                      !
    !17) xstr_loc (out)                                                       !
    !18) gstr_loc (out)                                                       !
    !19) list_fixed_coords_loc (out)                                          !
    !20) list_rigid_body_atoms_loc (out)                                      !
    !21) nile (in) : only read gradients/hessian/disps when nile>1            !
    !22) num_disconn_fragm (in)                                               !
    !23) ierr (out)                                                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) num_fixed_coords                                                     !
    ! 2) num_fixed_constr_rigid_body                                          !
    ! 3) ncons                                                                !
    ! 4) tsflag                                                               !
    ! 5) maxdiis                                                              !
    ! 6) np_p_global                                                          !
    ! 7) klist_global                                                         !
    ! 8) ktyp_global                                                          !
    ! 9) ut_global                                                            !
    !10) savtor_global                                                        !
    !11) xprim_global                                                         !
    !12) bmat_disconn_fragm_global                                            !
    !13) map_fixed_coords_global                                              !
    !14) list_fixed_coords_global                                             !
    !15) list_rigid_body_atoms_global                                         !
    !16) ictyp_global                                                         !
    !17) icc_global                                                           !
    !18) rcon_global                                                          !
    !19) disp_DI_global                                                       !
    !20) gradient_DI_global                                                   !
    !21) hess_DI_global                                                       !
    !22) hpad_global                                                          !
    !23) vmdel_global                                                         !
    !24) xstr_global                                                          !
    !25) gstr_global                                                          !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    implicit none

    integer,dimension(:,:),intent(out)               :: np_p_loc
    real(kind=di_dp),dimension(:,:),intent(out)      :: ut_loc
    integer,dimension(:,:),intent(out)               :: klist_loc
    real(kind=di_dp),dimension(:),intent(out)        :: xprim_loc
    real(kind=di_dp),dimension(:),intent(out)        :: savtor_loc
    integer,dimension(:),intent(out)                 :: ktyp_loc
    integer,dimension(:),intent(out)                 :: map_fixed_coords_loc
    real(kind=di_dp),dimension(:,:),intent(out)      :: bmat_disconn_fragm_loc
    integer,dimension(:,:),intent(out)               :: ictyp_loc
    integer,dimension(:,:),intent(out)               :: icc_loc
    real(kind=di_dp),dimension(:,:),intent(out)      :: rcon_loc
    real(kind=di_dp),dimension(:),intent(out)        :: disp_DI_loc
    real(kind=di_dp),dimension(:),intent(out)        :: gradient_DI_loc
    real(kind=di_dp),dimension(:,:),intent(out)      :: hess_DI_loc
    real(kind=di_dp),dimension(:,:),intent(out)      :: hpad_loc
    real(kind=di_dp),dimension(:),intent(out)        :: vmode_loc
    real(kind=di_dp),dimension(:,:),intent(out)      :: xstr_loc
    real(kind=di_dp),dimension(:,:),intent(out)      :: gstr_loc
    integer,dimension(:,:),intent(out)               :: list_fixed_coords_loc
    integer,dimension(:),intent(out)                 :: list_rigid_body_atoms_loc
    integer,intent(in)                                           :: num_disconn_fragm
    integer,intent(in)                                           :: nile
    integer,intent(out)                                          :: ierr
    !-------------------------------------------------------------------------!
    integer :: n1,n2,j
    logical :: data_found,failed
    !-------------------------------------------------------------------------!
    ierr = 0
    failed = .false.
    data_found = .false.
    n2 = -1 ! qoh: Initialise to remove compiler warnings

    ! NP_P
    if (di_on_root) then
       n1 = size(np_p_global,1)
       n2 = size(np_p_global,2)
       failed = (n1 /= size(np_p_loc,1) .or. n2 /= size(np_p_loc,2))
       data_found = n1>0 .and. n2>0
    endif
    if (failed) call deloc_utils_io_abort('Error for np_p_loc in deloc_utils_read_di_data')

    if (data_found) then
       if (di_on_root) np_p_loc = np_p_global
    endif ! data_found

    ! KLIST (known also as LIDELP)
    if (di_on_root) then
       n1 = size(klist_global,1)
       n2 = size(klist_global,2)
       failed = (n1 /= size(klist_loc,1) .or. n2 /= size(klist_loc,2)) 
       data_found = n1>0 .and. n2>0
    endif ! di_on_root
    if (failed) call deloc_utils_io_abort('Error for klist_loc in deloc_utils_read_di_data')

    if (data_found) then
       if (di_on_root) klist_loc  = klist_global 
    endif ! data_found

    ! KTYP
    if (di_on_root) then
       n1 = size(ktyp_global,1)
       failed = (n1 /= size(ktyp_loc,1)) 
       data_found = n1>0
    endif ! di_on_root
    if (failed) call deloc_utils_io_abort('Error for ktyp_loc in deloc_utils_read_di_data')

    if (data_found) then
       if (di_on_root) ktyp_loc   = ktyp_global
    endif ! data_found

    ! UT 2nd dimension (ndeg) could have changed
    if (di_on_root) then
       n1 = size(ut_global,1)
       failed = (n1 /= size(ut_loc,1))
       n2 = min(size(ut_global,2),size(ut_loc,2))
       data_found = n1>0 .and. n2>0
    endif ! di_on_root
    if (failed) call deloc_utils_io_abort('Error for ut_loc in deloc_utils_read_di_data')

    if (data_found) then
       if (di_on_root) then
          do j = 1,n2
             ut_loc(:,j)  = ut_global(:,j)
          enddo
       endif
    endif ! data_found

    ! SAVTOR
    if (di_on_root) then
       n1 = size(savtor_global,1)
       data_found = n1>0
       failed = (n1 /= size(savtor_loc,1)) 
    endif ! di_on_root
    if (failed) call deloc_utils_io_abort('Error for savtor_loc in deloc_utils_read_di_data')

    if (data_found) then
       if (di_on_root) savtor_loc = savtor_global
    endif ! data_found

    ! XPRIM
    if (di_on_root) then
       n1 = size(xprim_global,1)
       data_found = n1>0
       failed = (n1 /= size(xprim_loc,1)) 
    endif ! di_on_root
    if (failed) call deloc_utils_io_abort('Error for xprim_loc in deloc_utils_read_di_data')

    if (data_found) then
       if (di_on_root) xprim_loc  = xprim_global
    endif ! data_found

    if (num_disconn_fragm/=0) then

       ! BMAT_DISCONN_FRAGM
       if (di_on_root) then
          n1 = size(bmat_disconn_fragm_global,1)
          n2 = size(bmat_disconn_fragm_global,2)
          failed = (n1 /= size(bmat_disconn_fragm_loc,1) .or. n2 /= size(bmat_disconn_fragm_loc,2))
          data_found = n1>0 .and. n2>0
       endif ! di_on_root
       if (failed) call deloc_utils_io_abort('Error for bmat_disconn_fragm_loc in deloc_utils_read_di_data')

       if (data_found) then
          if (di_on_root) bmat_disconn_fragm_loc  = bmat_disconn_fragm_global
       endif ! data_found

    endif ! num_disconn_fragm

    if (num_fixed_coords/=0) then

       ! MAP_FIXED_COORDS
       if (di_on_root) then
          n1 = size(map_fixed_coords_global,1)
          failed = (n1 /= size(map_fixed_coords_loc,1))
          data_found = n1>0
       endif ! di_on_root
       if (failed) call deloc_utils_io_abort('Error for map_fixed_coords_loc in deloc_utils_read_di_data')

       if (data_found) then
          if (di_on_root) map_fixed_coords_loc = map_fixed_coords_global
       endif ! data_found

       ! LIST_FIXED_COORDS
       if (di_on_root) then
          n1 = size(list_fixed_coords_global,1)
          n2 = size(list_fixed_coords_global,2)
          failed = (n1 /= size(list_fixed_coords_loc,1) .or. n2 /= size(list_fixed_coords_loc,2))
          data_found = n1>0 .and. n2>0
       endif ! di_on_root
       if (failed) call deloc_utils_io_abort('Error for list_fixed_coords_loc in deloc_utils_read_di_data')

       if (data_found) then
          if (di_on_root) list_fixed_coords_loc   = list_fixed_coords_global
       endif ! data_found

    endif ! num_fixed_coords

    if (num_fixed_constr_rigid_body/=0) then

       ! LIST_RIGID_BODY_ATOMS
       if (di_on_root) then
          n1 = size(list_rigid_body_atoms_global,1)
          failed = (n1 /= size(list_rigid_body_atoms_loc,1)) 
          data_found = n1>0
       endif ! di_on_root
       if (failed) call deloc_utils_io_abort('Error for list_rigid_body_atoms_loc in deloc_utils_read_di_data')

       if (data_found) then
          if (di_on_root) list_rigid_body_atoms_loc  = list_rigid_body_atoms_global
       endif ! data_found

    endif ! num_fixed_constr_rigid_body

    if ((num_fixed_coords+ncons)/=0) then

       ! ICTYP
       if (di_on_root) then
          n1 = size(ictyp_global,1)
          n2 = size(ictyp_global,2)
          failed = (n1 /= size(ictyp_loc,1) .or. n2 /= size(ictyp_loc,2)) 
          data_found = n1>0 .and. n2>0
       endif ! di_on_root
       if (failed) call deloc_utils_io_abort('Error for ictyp_loc in deloc_utils_read_di_data')

       if (data_found) then
          if (di_on_root) ictyp_loc  = ictyp_global
       endif ! data_found

   endif ! num_fixed_coords+ncons

   if (ncons/=0) then 

       ! ICC
       if (di_on_root) then
          n1 = size(icc_global,1)
          n2 = size(icc_global,2)
          failed = (n1 /= size(icc_loc,1) .or. n2 /= size(icc_loc,2))
          data_found = n1>0 .and. n2>0
       endif ! di_on_root
       if (failed) call deloc_utils_io_abort('Error for icc_loc in deloc_utils_read_di_data')

       if (data_found) then
          if (di_on_root) icc_loc    = icc_global
       endif ! data_found

       ! RCON
       if (di_on_root) then
          n1 = size(rcon_global,1)
          n2 = size(rcon_global,2)
          failed = (n1 /= size(rcon_loc,1) .or. n2 /= size(rcon_loc,2))
          data_found = n1>0 .and. n2>0
       endif ! di_on_root
       if (failed) call deloc_utils_io_abort('Error for rcon_loc in deloc_utils_read_di_data')

       if (data_found) then
          if (di_on_root) rcon_loc   = rcon_global
       endif ! data_found

    endif ! ncons

    ! read if you have accurate gradients 
    if (nile>1) then

       ! DISP_DI
       if (di_on_root) then
          n1 = size(disp_DI_global,1)
          failed = (n1 /= size(disp_DI_loc,1))
          data_found = n1>0
       endif ! di_on_root
       if (failed) call deloc_utils_io_abort('Error for disp_DI_loc in deloc_utils_read_di_data')

       if (data_found) then
          if (di_on_root) disp_DI_loc   = disp_DI_global
       endif ! data_found

       ! GRADIENT_DI
       if (di_on_root) then
          n1 = size(gradient_DI_global,1)
          failed = (n1 /= size(gradient_DI_loc,1)) 
          data_found = n1>0
       endif ! di_on_root
       if (failed) call deloc_utils_io_abort('Error for gradient_DI_loc in deloc_utils_read_di_data')

       if (data_found) then
          if (di_on_root) gradient_DI_loc   = gradient_DI_global
       endif ! data_found

       ! HESS_DI
       if (di_on_root) then
          n1 = size(hess_DI_global,1)
          n2 = size(hess_DI_global,2)
          failed = (n1 /= size(hess_DI_loc,1) .or. n2 /= size(hess_DI_loc,2))
          data_found = n1>0 .and. n2>0
       endif ! di_on_root
       if (failed) call deloc_utils_io_abort('Error for hess_DI_loc in deloc_utils_read_di_data')

       if (data_found) then
          if (di_on_root) hess_DI_loc   = hess_DI_global
       endif ! data_found


       if (ncons>0) then

          ! HPAD
          if (di_on_root) then
             n1 = size(hpad_global,1)
             n2 = size(hpad_global,2)
             failed = (n1 /= size(hpad_loc,1) .or. n2 /= size(hpad_loc,2))
             data_found = n1>0 .and. n2>0
          endif ! di_on_root
          if (failed) call deloc_utils_io_abort('Error for hpad_loc in deloc_utils_read_di_data')

          if (data_found) then
             if (di_on_root) hpad_loc  = hpad_global
          endif ! data_found

       endif ! ncons

       if (tsflag) then 

          ! VMDEL
          if (di_on_root) then
             n1 = size(vmdel_global,1)
             failed = (n1 /= size(vmode_loc,1))
             data_found = n1>0
          endif ! di_on_root
          if (failed) call deloc_utils_io_abort('Error for vmode_loc in deloc_utils_read_di_data')

          if (data_found) then
             if (di_on_root) vmode_loc = vmdel_global
          endif ! data_found

       endif ! tsflag

       if (maxdiis>1) then

          ! XSTR
          if (di_on_root) then
             n1 = size(xstr_global,1)
             n2 = size(xstr_global,2)
             failed =  (n1 /= size(xstr_loc,1) .or. n2 /= size(xstr_loc,2))
             data_found = n1>0 .and. n2>0
          endif ! di_on_root
          if (failed) call deloc_utils_io_abort('Error for xstr_loc in deloc_utils_read_di_data')

          if (data_found) then
             if (di_on_root) xstr_loc = xstr_global
          endif ! data_found

          ! GSTR
          if (di_on_root) then
             n1 = size(gstr_global,1)
             n2 = size(gstr_global,2)
             failed = (n1/= size(gstr_loc,1) .or. n2 /= size(gstr_loc,2))
             data_found = n1>0 .and. n2>0
          endif ! di_on_root
          if (failed) call deloc_utils_io_abort('Error for gstr_loc in deloc_utils_read_di_data')

          if (data_found) then
             if (di_on_root) gstr_loc = gstr_global
          endif ! data_found

       endif ! maxdiis

    endif ! nile>1


    return
  end subroutine deloc_utils_read_di_data

!------------------------------------------------------------------------------
  subroutine deloc_utils_write_di_data(np_p_loc,ut_loc,klist_loc,xprim_loc,savtor_loc,ktyp_loc, &
     &     map_fixed_coords_loc,bmat_disconn_fragm_loc,ictyp_loc,icc_loc,rcon_loc, &
     &     disp_DI_loc,gradient_DI_loc,hess_DI_loc,list_fixed_coords_loc,list_rigid_body_atoms_loc, &
     &     num_disconn_fragm,nile,status,hpad_loc,vmode_loc,xstr_loc,gstr_loc)
    !=========================================================================!
    ! Save to memory a lot of data (used to be from PCHK file). If the        !
    ! corresponding global array exists and has nonzero dimensions, its       !
    ! contents is copied into the local (<>_loc) array and returned.          !
    ! _global arrays that are used for storage are initialized if necessary   !
    ! Only root node does the saving                                          !
    !                                                                         !
    !                           [wrpchk in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) np_p (in)                                                            !
    ! 2) ut_loc (in)                                                          !
    ! 3) klist_loc (in)                                                       !
    ! 4) xprim_loc (in)                                                       !
    ! 5) savtor_loc (in)                                                      !
    ! 6) ktyp_loc (in)                                                        !
    ! 7) map_fixed_coords_loc (in)                                            !
    ! 8) bmat_disconn_fragm_loc (in)                                          !
    ! 9) ictyp_loc (in)                                                       !
    !10) icc_loc (in)                                                         !
    !11) rcon_loc (in)                                                        !
    !12) disp_DI_loc (in)                                                     !
    !13) gradient_DI_loc (in)                                                 ! 
    !14) hess_DI_loc (in)                                                     !
    !15) list_fixed_coords_loc (in)                                           !
    !16) list_rigid_body_atoms_loc (in)                                       !
    !17) num_disconn_fragm (in)                                               !
    !18) nile (in) : save gradients/hessian/disps and update ncycle if nile>1 !
    !19) ierr (out)                                                           !
    !20) hpad_loc (optional in)                                               !
    !21) vmode_loc (optional in)                                              !
    !22) xstr_loc (optional in)                                               !
    !23) gstr_loc (optional in)                                               !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) num_disconn_fragm                                                    !
    ! 2) num_fixed_coords                                                     !
    ! 3) num_fixed_constr_rigid_body                                          !
    ! 4) ncons                                                                !
    ! 5) nile                                                                 !
    ! 6) tsflag                                                               !
    ! 7) maxdiis                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) np_p_global                                                          !
    ! 2) klist_global                                                         !
    ! 3) ktyp_global                                                          !
    ! 4) ut_global                                                            !
    ! 5) savtor_global                                                        !
    ! 6) xprim_global                                                         !
    ! 7) bmat_disconn_fragm_global                                            !
    ! 8) map_fixed_coords_global                                              !
    ! 9) list_fixed_coords_global                                             !
    !10) list_rigid_body_atoms_global                                         !
    !11) ictyp_global                                                         !
    !12) icc_global                                                           !
    !13) rcon_global                                                          !
    !14) disp_DI_global                                                       !
    !15) gradient_DI_global                                                   !
    !16) hess_DI_global                                                       !
    !17) hpad_global                                                          !
    !18) vmdel_global                                                         !
    !19) xstr_global                                                          !
    !20) gstr_global                                                          !
    !21) ncycle                                                               !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    implicit none

    integer,dimension(:,:),intent(in)                     :: np_p_loc
    real(kind=di_dp),dimension(:,:),intent(in)            :: ut_loc
    integer,dimension(:,:),intent(in)                     :: klist_loc
    real(kind=di_dp),dimension(:),intent(in)              :: xprim_loc
    real(kind=di_dp),dimension(:),intent(in)              :: savtor_loc
    integer,dimension(:),intent(in)                       :: ktyp_loc
    integer,dimension(:),intent(in)                       :: map_fixed_coords_loc
    real(kind=di_dp),dimension(:,:),intent(in)            :: bmat_disconn_fragm_loc
    integer,dimension(:,:),intent(in)                     :: ictyp_loc
    integer,dimension(:,:),intent(in)                     :: icc_loc
    real(kind=di_dp),dimension(:,:),intent(in)            :: rcon_loc
    real(kind=di_dp),dimension(:),intent(in)              :: disp_DI_loc
    real(kind=di_dp),dimension(:),intent(in)              :: gradient_DI_loc
    real(kind=di_dp),dimension(:,:),intent(in)            :: hess_DI_loc
    integer,dimension(:,:),intent(in)                     :: list_fixed_coords_loc
    integer,dimension(:),intent(in)                       :: list_rigid_body_atoms_loc
    integer,intent(in)                                    :: num_disconn_fragm
    integer,intent(in)                                    :: nile
    integer,intent(out)                                   :: status 
    real(kind=di_dp),dimension(:,:),intent(in),optional   :: hpad_loc
    real(kind=di_dp),dimension(:),intent(in),optional     :: vmode_loc
    real(kind=di_dp),dimension(:,:),intent(in),optional   :: xstr_loc
    real(kind=di_dp),dimension(:,:),intent(in),optional   :: gstr_loc
    !-------------------------------------------------------------------------!
    integer :: n1,n2,ierr
    logical :: bad_array
    !-------------------------------------------------------------------------!
    status = 0
    if (.not. di_on_root) return
       
    ! augment ncycle only if there is no failure in di_algor_back_trans or get_step
    if (nile>=0) ncycle = ncycle + 1

    n1 = size(np_p_loc,1)
    n2 = size(np_p_loc,2)
    bad_array = .true.
    if (allocated(np_p_global)) then
       if (size(np_p_global,1)/=n1 .or. size(np_p_global,2)/=n2) then
          deallocate(np_p_global,stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in deallocating np_p_global in deloc_utils_write_di_data')
       else
          bad_array = .false.
       endif
    endif
    if (bad_array) then
       allocate(np_p_global(1:n1,1:n2),stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating np_p_global in deloc_utils_write_di_data')
    endif
    np_p_global   = np_p_loc


    n1 = size(klist_loc,1)
    n2 = size(klist_loc,2)
    bad_array = .true.
    if (allocated(klist_global)) then
       if (size(klist_global,1)/=n1 .or. size(klist_global,2)/=n2) then
          deallocate(klist_global,stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in deallocating klist_global in deloc_utils_write_di_data')
       else
          bad_array = .false.
       endif
    endif
    if (bad_array) then
       allocate(klist_global(1:n1,1:n2),stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating klist_global in deloc_utils_write_di_data')
    endif
    klist_global  = klist_loc 


    n1 = size(ktyp_loc,1)
    bad_array = .true.
    if (allocated(ktyp_global)) then
       if (size(ktyp_global,1)/=n1) then
          deallocate(ktyp_global,stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ktyp_global in deloc_utils_write_di_data')
       else
          bad_array = .false.
       endif
    endif
    if (bad_array) then
       allocate(ktyp_global(1:n1),stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating ktyp_global in deloc_utils_write_di_data')
    endif
    ktyp_global   = ktyp_loc


    n1 = size(ut_loc,1)
    n2 = size(ut_loc,2)
    bad_array = .true.
    if (allocated(ut_global)) then
       if (size(ut_global,1)/=n1 .or. size(ut_global,2)/=n2) then
          deallocate(ut_global,stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ut_global in deloc_utils_write_di_data')
       else
          bad_array = .false.
       endif
    endif
    if (bad_array) then
       allocate(ut_global(1:n1,1:n2),stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating ut_global in deloc_utils_write_di_data')
    endif
    ut_global     = ut_loc


    n1 = size(savtor_loc,1)
    bad_array = .true.
    if (allocated(savtor_global)) then
       if (size(savtor_global,1)/=n1) then
          deallocate(savtor_global,stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in deallocating savtor_global in deloc_utils_write_di_data')
       else
          bad_array = .false.
       endif
    endif
    if (bad_array) then
       allocate(savtor_global(1:n1),stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating savtor_global in deloc_utils_write_di_data')
    endif
    savtor_global = savtor_loc

    n1 = size(xprim_loc,1)
    bad_array = .true.
    if (allocated(xprim_global)) then
       if (size(xprim_global,1)/=n1) then
          deallocate(xprim_global,stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in deallocating xprim_global in deloc_utils_write_di_data')
       else
          bad_array = .false.
       endif
    endif
    if (bad_array) then
       allocate(xprim_global(1:n1),stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating xprim_global in deloc_utils_write_di_data')
    endif
    xprim_global  = xprim_loc

    if (num_disconn_fragm/=0) then
       n1 = size(bmat_disconn_fragm_loc,1)
       n2 = size(bmat_disconn_fragm_loc,2)
       bad_array = .true.
       if (allocated(bmat_disconn_fragm_global)) then
          if (size(bmat_disconn_fragm_global,1)/=n1 .or. size(bmat_disconn_fragm_global,2)/=n2) then
             deallocate(bmat_disconn_fragm_global,stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bmat_disconn_fragm_global in deloc_utils_write_di_data')
          else
             bad_array = .false.
          endif
       endif
       if (bad_array) then
          allocate(bmat_disconn_fragm_global(1:n1,1:n2),stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in allocating bmat_disconn_fragm_global in deloc_utils_write_di_data')
       endif
       bmat_disconn_fragm_global  = bmat_disconn_fragm_loc
    endif

    if (num_fixed_coords/=0) then
       n1 = size(map_fixed_coords_loc,1)
       bad_array = .true.
       if (allocated(map_fixed_coords_global)) then
          if (size(map_fixed_coords_global,1)/=n1) then
             deallocate(map_fixed_coords_global,stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort('Error in deallocating map_fixed_coords_global in deloc_utils_write_di_data')
          else
             bad_array = .false.
          endif
       endif
       if (bad_array) then
          allocate(map_fixed_coords_global(1:n1),stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in allocating map_fixed_coords_global in deloc_utils_write_di_data')
       endif
       map_fixed_coords_global = map_fixed_coords_loc

       n1 = size(list_fixed_coords_loc,1)
       n2 = size(list_fixed_coords_loc,2)
       bad_array = .true.
       if (allocated(list_fixed_coords_global)) then
          if (size(list_fixed_coords_global,1)/=n1 .or. size(list_fixed_coords_global,2)/=n2) then
             deallocate(list_fixed_coords_global,stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort('Error in deallocating list_fixed_coords_global in deloc_utils_write_di_data')
          else
             bad_array = .false.
          endif
       endif
       if (bad_array) then
          allocate(list_fixed_coords_global(1:n1,1:n2),stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in allocating list_fixed_coords_global in deloc_utils_write_di_data')
       endif
       list_fixed_coords_global   = list_fixed_coords_loc
    endif

    if (num_fixed_constr_rigid_body/=0) then
       n1 = size(list_rigid_body_atoms_loc,1)
       bad_array = .true.
       if (allocated(list_rigid_body_atoms_global)) then
          if (size(list_rigid_body_atoms_global,1)/=n1) then
             deallocate(list_rigid_body_atoms_global,stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort &
                  & ('Error in deallocating list_rigid_body_atoms_global in deloc_utils_write_di_data')
          else
             bad_array = .false.
          endif
       endif
       if (bad_array) then
          allocate(list_rigid_body_atoms_global(1:n1),stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in allocating list_rigid_body_atoms_global in deloc_utils_write_di_data')
       endif
       list_rigid_body_atoms_global  = list_rigid_body_atoms_loc
    endif

    if ((num_fixed_coords+ncons)/=0) then
       n1 = size(ictyp_loc,1)
       n2 = size(ictyp_loc,2)
       bad_array = .true.
       if (allocated(ictyp_global)) then
          if (size(ictyp_global,1)/=n1 .or. size(ictyp_global,2)/=n2) then
             deallocate(ictyp_global,stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ictyp_global in deloc_utils_write_di_data')
          else
             bad_array = .false.
          endif
       endif
       if (bad_array) then
          allocate(ictyp_global(1:n1,1:n2),stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in allocating ictyp_global in deloc_utils_write_di_data')
       endif
       ictyp_global  = ictyp_loc
    endif

    if (ncons/=0) then 
       n1 = size(icc_loc,1)
       n2 = size(icc_loc,2)
       bad_array = .true.
       if (allocated(icc_global)) then
          if (size(icc_global,1)/=n1 .or. size(icc_global,2)/=n2) then
             deallocate(icc_global,stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icc_global in deloc_utils_write_di_data')
          else
             bad_array = .false.
          endif
       endif
       if (bad_array) then
          allocate(icc_global(1:n1,1:n2),stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in allocating icc_global in deloc_utils_write_di_data')
       endif
       icc_global    = icc_loc

       n1 = size(rcon_loc,1)
       n2 = size(rcon_loc,2)
       bad_array = .true.
       if (allocated(rcon_global)) then
          if (size(rcon_global,1)/=n1 .or. size(rcon_global,2)/=n2) then
             deallocate(rcon_global,stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort('Error in deallocating rcon_global in deloc_utils_write_di_data')
          else
             bad_array = .false.
          endif
       endif
       if (bad_array) then
          allocate(rcon_global(1:n1,1:n2),stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in allocating rcon_global in deloc_utils_write_di_data')
       endif
       rcon_global   = rcon_loc
    endif

    ! write if you have accurate gradients 
    if (nile>=1) then
       n1 = size(disp_DI_loc,1)
       bad_array = .true.
       if (allocated(disp_DI_global)) then
          if (size(disp_DI_global,1)/=n1) then
             deallocate(disp_DI_global,stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort('Error in deallocating disp_DI_global in deloc_utils_write_di_data')
          else
             bad_array = .false.
          endif
       endif
       if (bad_array) then
          allocate(disp_DI_global(1:n1),stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in allocating disp_DI_global in deloc_utils_write_di_data')
       endif
       disp_DI_global   = disp_DI_loc

       n1 = size(gradient_DI_loc,1)
       bad_array = .true.
       if (allocated(gradient_DI_global)) then
          if (size(gradient_DI_global,1)/=n1) then
             deallocate(gradient_DI_global,stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort('Error in deallocating gradient_DI_global in deloc_utils_write_di_data')
          else
             bad_array = .false.
          endif
       endif
       if (bad_array) then
          allocate(gradient_DI_global(1:n1),stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in allocating gradient_DI_global in deloc_utils_write_di_data')
       endif
       gradient_DI_global   = gradient_DI_loc

       n1 = size(hess_DI_loc,1)
       n2 = size(hess_DI_loc,2)
       bad_array = .true.
       if (allocated(hess_DI_global)) then
          if (size(hess_DI_global,1)/=n1 .or. size(hess_DI_global,2)/=n2) then
             deallocate(hess_DI_global,stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hess_DI_global in deloc_utils_write_di_data')
          else
             bad_array = .false.
          endif
       endif
       if (bad_array) then
          allocate(hess_DI_global(1:n1,1:n2),stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in allocating hess_DI_global in deloc_utils_write_di_data')
       endif
       hess_DI_global   = hess_DI_loc

       if (ncons>0 .and. present(hpad_loc)) then
          n1 = size(hpad_loc,1)
          n2 = size(hpad_loc,2)
          bad_array = .true.
          if (allocated(hpad_global)) then
             if (size(hpad_global,1)/=n1 .or. size(hpad_global,2)/=n2) then
                deallocate(hpad_global,stat=ierr)
                if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hpad_global in deloc_utils_write_di_data')
             else
                bad_array = .false.
             endif
          endif
          if (bad_array) then
             allocate(hpad_global(1:n1,1:n2),stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort('Error in allocating hpad_global in deloc_utils_write_di_data')
          endif
          hpad_global  = hpad_loc
       endif

       if (tsflag .and. present(vmode_loc)) then 
          n1 = size(vmode_loc,1)
          bad_array = .true.
          if (allocated(vmdel_global)) then
             if (size(vmdel_global,1)/=n1) then
                deallocate(vmdel_global,stat=ierr)
                if (ierr/=0) call deloc_utils_io_abort('Error in deallocating vmdel_global in deloc_utils_write_di_data')
             else
                bad_array = .false.
             endif
          endif
          if (bad_array) then
             allocate(vmdel_global(1:n1),stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort('Error in allocating vmdel_global in deloc_utils_write_di_data')
          endif
          vmdel_global = vmode_loc
       endif

       if (maxdiis>1) then
          if (present(xstr_loc)) then
             n1 = size(xstr_loc,1)
             n2 = size(xstr_loc,2)
             bad_array = .true.
             if (allocated(xstr_global)) then
                if (size(xstr_global,1)/=n1 .or. size(xstr_global,2)/=n2) then
                   deallocate(xstr_global,stat=ierr)
                   if (ierr/=0) call deloc_utils_io_abort('Error in deallocating xstr_global in deloc_utils_write_di_data')
                else
                   bad_array = .false.
                endif
             endif
             if (bad_array) then
                allocate(xstr_global(1:n1,1:n2),stat=ierr)
                if (ierr/=0) call deloc_utils_io_abort('Error in allocating xstr_global in deloc_utils_write_di_data')
             endif
             xstr_global = xstr_loc
          endif

          if (present(gstr_loc)) then
             n1 = size(gstr_loc,1)
             n2 = size(gstr_loc,2)
             bad_array = .true.
             if (allocated(gstr_global)) then
                if (size(gstr_global,1)/=n1 .or. size(gstr_global,2)/=n2) then
                   deallocate(gstr_global,stat=ierr)
                   if (ierr/=0) call deloc_utils_io_abort('Error in deallocating gstr_global in deloc_utils_write_di_data')
                else
                   bad_array = .false.
                endif
             endif
             if (bad_array) then
                allocate(gstr_global(1:n1,1:n2),stat=ierr)
                if (ierr/=0) call deloc_utils_io_abort('Error in allocating gstr_global in deloc_utils_write_di_data')
             endif
             gstr_global = gstr_loc
          endif
       endif
    endif

    return
  end subroutine deloc_utils_write_di_data

!------------------------------------------------------------------------------
  
  subroutine deloc_utils_read_np_p (np_p_loc,ierr)
    !=========================================================================!
    ! Restore indices for coordinate rationalization                          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) np_p_loc (out)                                                       !
    ! 2) ierr (out)                                                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) np_p_global                                                          !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    implicit none
    integer,dimension(:,:),intent(out)           :: np_p_loc
    integer,intent(out)                          :: ierr

    !-------------------------------------------------------------------------!
    integer :: n1,n2
    logical :: data_found
    !-------------------------------------------------------------------------!
    ierr = 0
    data_found = .false.
    if (di_on_root) then
       n1 = size(np_p_global,1)
       n2 = size(np_p_global,2)
       if (n1 /= size(np_p_loc,1) .or. n2 /= size(np_p_loc,2)) &
          call deloc_utils_io_abort('Error for np_p_loc in deloc_utils_read_np_p')
       data_found = n1>0 .and. n2>0
    endif

    if (data_found) then
       if (di_on_root) np_p_loc = np_p_global
    endif

    return
  end subroutine deloc_utils_read_np_p

  subroutine di_utils_svd_matrix_invert(a,n,ierr)
    !=========================================================================!
    ! Invert A matrix using Singular value decomposition algorithm            !
    !                                                                         !
    !                        [invmatd in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) A (inout)  -  input matrix (returns BB(T)-1  ==A-1)                  !
    ! 2) N (in)     -  dimension of A                                         !
    ! 3) IErr (out) -  error flag                                             !
    !             0 - matrix successfully inverted                            !
    !            -1 - SVD decomposition failed                                !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                 :: n
    real (kind=di_dp),dimension(1:n,1:n),intent(inout) :: a
    integer,intent(out)                                :: ierr

    integer :: lw,i,k,j,istat
    real (kind=di_dp),dimension(:),allocatable         :: w,s
    real (kind=di_dp),dimension(:,:),allocatable       :: u,v
    !-------------------------------------------------------------------------!
    ! Relative Threshold of Zero
    real(kind=di_dp), parameter :: svdthr = 1.0e-14_di_dp
    real(kind=di_dp) :: ss
    !-------------------------------------------------------------------------!

    !  Singular value decomposition:

    !       A = U * S * transpose(V)
    !

    lw = 6*n
    allocate(w(1:lw),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating W in di_utils_svd_matrix_invert')
    w = zero
    allocate(s(1:n),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating S in di_utils_svd_matrix_invert')
    s = zero
    allocate(u(1:n,1:n),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating V in di_utils_svd_matrix_invert')
    u = zero
    allocate(v(1:n,1:n),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating V in di_utils_svd_matrix_invert')
    v = zero

    call dgesvd('A','A', n,n, a,n, s,u,n, v,n, w,lw, ierr)
    if (ierr/=0) go to 999

    !  u*u(T) = 1;   v*v(T) = 1

    !   calculate A-1
    !      A^(-1) = V*S^(-1)*U(T)

    do i=1,n
       if (s(i)<s(1)*svdthr ) then
          if (di_iprint>4 .and. di_on_root) write(di_stdout,*) 'small diagonal element ',i,s(i)
          s(i)=zero
       endif
    enddo

    !  find inverse of A  by columns
    do k = 1,n
       do j = 1,n

          ! 1/s * U(T)
          if (s(j)/=zero) then
              w(j) = u(k,j)/s(j)
          else
              w(j) = zero
          endif
       end do

       ! V* A(J,K)
       do j=1,n
          ss = zero
          do i = 1,n
             ss = ss + v(i,j)*w(i)
          end do
          a(j,k) = ss
       end do

    end do  ! k column

 999 continue

    deallocate (w,stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in deallocating W in di_utils_svd_matrix_invert')
    deallocate (s,stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in deallocating S in di_utils_svd_matrix_invert')
    deallocate (u,stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in deallocating U in di_utils_svd_matrix_invert')
    deallocate (v,stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in deallocating V in di_utils_svd_matrix_invert')

    return
  end subroutine di_utils_svd_matrix_invert

!------------------------------------------------------------------------------

  subroutine deloc_utils_pack_atom_indices (compound_index,i1,i2,i3,i4)
    !=========================================================================!
    ! Pack indices into integer compound_index                                !
    ! Belongs in urils, rather than algor, because it can be used by such     !
    ! tools as MDF reader                                                     !
    !                                                                         !
    !                           [pacel in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) compound_index (out) : compound index                                !
    ! 2) i1 (in) : 1st integer (cell index along 'a')                         !
    ! 3) i2 (in) : 2nd integer (cell index along 'b')                         !
    ! 4) i3 (in) : 3rd integer (cell index along 'c')                         !
    ! 5) i4 (in) : 4th integer (atom number)                                  !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: i1,i2,i3 < 31, i4 < 4095                          !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    integer,intent(out)  :: compound_index
    integer,intent(in)   :: i1
    integer,intent(in)   :: i2
    integer,intent(in)   :: i3
    integer,intent(in)   :: i4
    !-------------------------------------------------------------------------!

    if (i1>=31) then
       call deloc_utils_io_abort('Error in deloc_utils_pack_atom_indices - 1st input argument out of bounds')
    elseif (i2>=31) then
       call deloc_utils_io_abort('Error in deloc_utils_pack_atom_indices - 2nd input argument out of bounds')
    elseif (i3>=31) then
       call deloc_utils_io_abort('Error in deloc_utils_pack_atom_indices - 3rd input argument out of bounds')
    elseif (i4>=4095) then
       call deloc_utils_io_abort('Error in deloc_utils_pack_atom_indices - 4th input argument out of bounds')
    endif

    compound_index = 32 + i1 + ishft(32+i2,6) + ishft(32+i3,12) +  ishft(32+i4,18)

    return
  end subroutine deloc_utils_pack_atom_indices
!------------------------------------------------------------------------------

  subroutine deloc_utils_unpack_atom_indices (compound_index,i1,i2,i3,i4)
    !=========================================================================!
    ! unpack from integer compound_index  i1,i2,i3 < 31                       !
    !                                     i4       < 4095                     !
    ! Belongs in urils, rather than algor, because it can be used by such     !
    ! tools as MDF reader                                                     !
    !                                                                         !
    !                          [unpacel in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) compound_index (in) : compound index                                 !
    ! 2) i1 (out) : 1st integer (cell index along 'a')                        !
    ! 3) i2 (out) : 2nd integer (cell index along 'b')                        !
    ! 4) i3 (out) : 3rd integer (cell index along 'c')                        !
    ! 5) i4 (out) : 4th integer (atom number)                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  none                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    integer,intent(in)  :: compound_index
    integer,intent(out) :: i1
    integer,intent(out) :: i2
    integer,intent(out) :: i3
    integer,intent(out) :: i4
    !-------------------------------------------------------------------------!

    i1 = iand(        compound_index,       63) -32
    i2 = iand( ishft( compound_index,  -6), 63) -32
    i3 = iand( ishft( compound_index, -12), 63) -32
    i4 = iand( ishft( compound_index, -18), 4095) -32

    if (i1>=31) then
       call deloc_utils_io_abort('Error in deloc_utils_unpack_atom_indices - 1st output argument out of bounds')
    elseif (i2>=31) then
       call deloc_utils_io_abort('Error in deloc_utils_unpack_atom_indices - 2nd output argument out of bounds')
    elseif (i3>=31) then
       call deloc_utils_io_abort('Error in deloc_utils_unpack_atom_indices - 3rd output argument out of bounds')
    elseif (i4>=4095) then
       call deloc_utils_io_abort('Error in deloc_utils_unpack_atom_indices - 4th output argument out of bounds')
    endif

    return
  end subroutine deloc_utils_unpack_atom_indices

end module geometry_deloc_utils
