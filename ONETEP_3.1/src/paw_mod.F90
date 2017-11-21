! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!             Projector-Augmented Wave module                    !
!                                                                !
! This module contains routines associated with the use of the   !
! Projector-Augmented Wave method, including initialisation,     !
! calculation of the augmentation density, calculation of the    !
! nonlocal energy terms, and calculation of the forces.          !
!----------------------------------------------------------------!
! This module was created by Nicholas Hine in May 2010.          !
!================================================================!

module paw

  use constants, only: DP, PI, SQRT_PI, NORMAL, VERBOSE, stdout
  use paw_shape, only: PAW_SHAPE_INFO
  use projectors, only: PROJECTOR_SET

  implicit none

  private

  ! Type describing radial grids
  type RADIAL_GRID

     ! ndmh: number of points in this radial grid
     integer :: npt

     ! ndmh: step sizes for radial grid
     real(kind=DP) :: rad_step, log_step

     ! ndmh: positions on radial grid
     real(kind=DP), pointer :: r(:)

     ! ndmh: "spacings" of radial grid
     real(kind=DP), pointer :: rab(:)

  end type RADIAL_GRID

  ! Type describing PAW species
  type PAW_SPECIES

     ! ndmh: dataset name
     character(len=64) :: ds_name

     ! ndmh: core wavefunctions filename
     character(len=64) :: core_wf_name

     ! ndmh: charge of core of atom
     real(kind=DP) :: ion_charge

     ! ndmh: atomic number of atom
     integer :: atomic_number

     ! ndmh: PAW cutoff radius
     real(kind=DP) :: rcut

     ! ndmh: PAW projector max radius
     real(kind=DP) :: proj_rcut

     ! ndmh: Number of projectors / partial waves (counting n, l only)
     integer :: npw

     ! ndmh: Number of projectors / partial waves (counting n, l, m)
     integer :: npw_tot

     ! ndmh: Angular momentum of each partial wave / projector pair
     integer, pointer :: l_pw(:)
     integer, pointer :: l_pw_tot(:),m_pw_tot(:)
     integer, pointer :: ipw_tot(:)

     ! ndmh: Maximum angular momentum of any partial wave / projector pair
     integer :: lmax

     ! ndmh: Number of meshes defined
     integer :: ngrid

     ! ndmh: Which mesh the partial waves use
     integer :: phi_grid

     ! ndmh: Which mesh the projectors use
     integer :: proj_grid

     ! ndmh: Which mesh the core densities use
     integer :: core_den_grid

     ! ndmh: Which mesh the local potential uses
     integer :: vhntzc_grid

     ! ndmh: The grid to use for the shape function
     integer :: shape_grid

     ! ndmh: The format of the vhntzc potential
     integer :: vhntzc_format

     ! ndmh: Array to hold information about grids
     type(RADIAL_GRID), pointer :: grid(:)

     ! ndmh: Type to hold information about shape function
     type(PAW_SHAPE_INFO) :: shape

     integer :: n_recip_pts
     real(kind=DP) :: g_max
     real(kind=DP) :: inv_g_spacing

     ! ndmh: AE and PS partial waves on radial grid
     real(kind=DP), pointer :: phi_rad(:,:)
     real(kind=DP), pointer :: tphi_rad(:,:)

     ! ndmh: projectors on radial grid
     real(kind=DP), pointer :: tproj_rad(:,:)
     real(kind=DP), pointer :: tproj_recip(:,:)

     ! ndmh: core densities on radial grids (real/recip)
     real(kind=DP), pointer :: core_den_rad(:)
     real(kind=DP), pointer :: tcore_den_rad(:)
     real(kind=DP), pointer :: tcore_den_recip(:)

     ! ndmh: Frozen part of nonlocal energies Dij
     real(kind=DP), pointer :: dij0(:,:)

     ! ndmh: Initial guess for projector density kernel
     real(kind=DP), pointer :: rhoij0(:,:)

     ! ndmh: Hartree potential of core, real and recip radial grids
     real(kind=DP), pointer :: vhntzc_rad(:)
     real(kind=DP), pointer :: vhntzc_recip(:)

     ! ndmh: for augmentation charge:
     ! ndmh: \int (phi_i(r)\phi_j(r)-\tphi_i(r)\tphi_j(r))r^L dr
     real(kind=DP), pointer :: aug_nLij(:,:,:)

     ! ndmh: for onsite Hartree energy
     real(kind=DP), pointer :: e_ijkl(:,:,:,:)

     ! ndmh: whether AE core density is nonzero
     logical :: core_charge
     logical :: core_charge_calculated

     ! ndmh: whether PS core density is nonzero
     logical :: tcore_charge

     ! ndmh: XC energy of AE core density
     real(kind=DP) :: exc_core

     ! ---- FOR CORE LEVEL RECONSTRUCTION ----
     logical :: core_wvfns_exist

     ! ndmh: number core orbitals on radial grid (and number including m)
     integer :: n_core_wfs, n_core_wfs_tot

     ! ndmh: angular momentum values of each core wf
     integer, pointer :: l_core_wf(:)
     integer, pointer :: l_core_wf_tot(:)
     integer, pointer :: m_core_wf_tot(:)
     integer, pointer :: icore_wf_tot(:)

     ! ndmh: principle quantum number of each state (for labels)
     integer, pointer :: n_core_wf(:)

     ! ndmh: Array to hold information about core wf grids
     integer :: ngrid_core
     integer :: core_wf_grid
     type(RADIAL_GRID), pointer :: grid_core(:)
     real(kind=DP) :: rcut_core

     ! ndmh: charge of core of atom
     real(kind=DP) :: core_wf_charge

     ! ndmh: Core orbitals on radial grid
     real(kind=DP), pointer :: core_wf_rad(:,:)
     real(kind=DP), pointer :: core_wf_recip(:,:)
     real(kind=DP), pointer :: core_wf_eig(:)
     real(kind=DP), pointer :: core_wf_occ(:)

     ! lr408: \int (psi_c(r)\phi_j(r)-\psi_c(r)\tphi_j(r))r^L dr
     real(kind=DP), pointer :: core_aug_nLij(:,:,:)

  end type PAW_SPECIES

  ! Array to hold all the PAW species in the system
  type(PAW_SPECIES),allocatable :: paw_sp(:)

  ! Maximum number of partial waves on any atom (including m_i)
  integer :: max_paw_proj_tot
  integer :: max_core_wf_tot

  real(kind=DP), parameter :: sqrt_4pi = 2.0_DP*SQRT_PI
  real(kind=DP), parameter :: inv_sqrt_4pi = 1.0_DP/sqrt_4pi

  ! Projector type to describe PAW projectors
  type(PROJECTOR_SET), public :: paw_projectors
  type(PROJECTOR_SET), public :: paw_core_wvfns

  ! Public subroutines
  public :: paw_read_species
  public :: paw_species_exit
  public :: paw_tcore_hartree_on_grid
  public :: paw_tcore_density
  public :: paw_species_init_proj
  public :: paw_projector_overlap
  public :: paw_position_operator
  public :: paw_grad_operator
  public :: paw_projector_denskern_init
  public :: paw_nonlocal_energies
  public :: paw_tcore_hartree_calc_forces
  public :: paw_nlcc_calculate_forces
  public :: paw_sphere_density_on_grid
  public :: paw_atom_aug_den
  public :: paw_atom_aug_integrals
  public :: paw_atom_aug_force

  public :: paw_get_projector_info
  public :: paw_get_locpot_rad
  public :: paw_get_core_den_rad
  public :: paw_get_projectors_q
  public :: paw_get_aug_funcs
  public :: paw_exc_core_atom
  public :: paw_dij_hartree_atom
  public :: paw_dij_xc_atom

  public :: paw_species_init_core_wvfns !lr408
  public :: paw_core_position_operator !lr408

contains

  subroutine paw_read_species(elements)

    !==================================================================!
    ! This subroutine reads the elements array and works out what PAW  !
    ! species are present, allocates storage for them in the paw_sp    !
    ! array and then read the PAW dataset file for each species into   !
    ! the paw_sp array                                                 !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  elements (inout) : list of atoms present                        !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 17/05/2010.                          !
    !==================================================================!

    use comms, only: pub_on_root
    use gaunt_coeff, only: gaunt_init
    use ion, only: ELEMENT
    use rundat, only: pub_nlcc, pub_aug, pub_usp, pub_nhat_in_xc
    use simulation_cell, only: pub_cell
    use utils, only: utils_abort, utils_unit, utils_alloc_check, &
         utils_dealloc_check

    implicit none

    ! Arguments
    type(ELEMENT), intent(inout)   :: elements(pub_cell%nat)

    ! Local Variables
    logical :: found
    integer :: ierr
    integer :: iat, jat
    integer :: isp,nsp,msp
    integer :: lmax

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering paw_read_species'
#endif

    ! count unique species
    nsp = 0
    do iat=1,pub_cell%nat
       found = .false.
       do jat=1,iat-1
          if (elements(iat)%pseudo_name==elements(jat)%pseudo_name) then
             found = .true.
             exit
          end if
       end do
       if (.not.found) then
          nsp = nsp + 1
       end if
    end do

    allocate(paw_sp(nsp),stat=ierr)
    call utils_alloc_check('paw_read_species','paw_sp',ierr)

    ! set dataset file names
    paw_sp(:)%ds_name = ''
    msp = 0
    do iat=1,pub_cell%nat
       found = .false.
       do isp=1,nsp
          if (elements(iat)%pseudo_name==paw_sp(isp)%ds_name) then
             found = .true.
             exit
          end if
       end do
       if (.not.found) then
          msp = msp + 1
          paw_sp(msp)%ds_name = elements(iat)%pseudo_name
          paw_sp(msp)%core_wf_name = elements(iat)%core_wf_name
       end if
    end do

    ! Consistency check
    if (msp/=nsp) then
       call utils_abort('Error: Inconsistent dataset names in paw_read_species')
    end if

    ! Read Datasets
    if (pub_on_root) then
       do isp=1,nsp
          call paw_read_dataset(paw_sp(isp))
       end do
    end if

    ! Read Core wavefunctions
    if (pub_on_root) then
       paw_sp(:)%core_wvfns_exist = .false.
       do isp=1,nsp
          if (paw_sp(isp)%core_wf_name/='NONE') then

             call paw_read_core_wfs(paw_sp(isp))

             if (paw_sp(isp)%n_core_wfs > 0) then
                paw_sp(isp)%core_wvfns_exist = .true.
             else
                paw_sp(isp)%n_core_wfs = 0
                paw_sp(isp)%n_core_wfs_tot = 0
                paw_sp(isp)%ngrid_core = 0
                paw_sp(isp)%core_wf_grid = 0
                paw_sp(isp)%rcut_core = 0.0_DP
             end if
          else
             paw_sp(isp)%n_core_wfs = 0
             paw_sp(isp)%n_core_wfs_tot = 0
             paw_sp(isp)%ngrid_core = 0
             paw_sp(isp)%core_wf_grid = 0
             paw_sp(isp)%rcut_core = 0.0_DP
          end if
       end do
    end if

    ! Broadcast data read from dataset from root node to all other nodes
    do isp=1,nsp
       call paw_bcast_dataset(paw_sp(isp))
    end do

    ! Copy relevant information back into the species array on all nodes
    elements(:)%pspecies_number = -1
    do iat=1,pub_cell%nat

       ! Loop over species found and check if name matches
       do isp=1,nsp
          if (paw_sp(isp)%ds_name == elements(iat)%pseudo_name) then
             elements(iat)%pspecies_number = isp
             elements(iat)%atomic_number = paw_sp(isp)%atomic_number
             elements(iat)%ion_charge = paw_sp(isp)%ion_charge
             elements(iat)%npawpws = paw_sp(isp)%npw_tot
             elements(iat)%nprojectors = 0
             elements(iat)%max_core_radius = paw_sp(isp)%proj_rcut
             elements(iat)%ncorewfs = paw_sp(isp)%n_core_wfs_tot
             elements(iat)%max_core_wf_radius = paw_sp(isp)%rcut_core
             exit
          end if
       end do

       ! Check we have found a species number for every species
       if (elements(iat)%pspecies_number == -1) then
          call utils_abort('Error in paw_read_species: No species found to &
               &match'//elements(iat)%pseudo_name)
       end if

    end do

    ! Count total number of PAW partial waves in system
    pub_cell%num_projectors = 0
    pub_cell%num_pawpws = 0
    pub_cell%num_corewfs = 0
    do iat=1,pub_cell%nat
       pub_cell%num_pawpws = pub_cell%num_pawpws + &
            paw_sp(elements(iat)%pspecies_number)%npw_tot

       pub_cell%num_corewfs = pub_cell%num_corewfs + &
            paw_sp(elements(iat)%pspecies_number)%n_core_wfs_tot
    end do

    ! Display species information in table
    if (pub_on_root) call internal_print_datasets

    ! Count highest angular momentum over all species
    lmax = 0
    do isp=1,pub_cell%num_pspecies
       lmax=max(paw_sp(isp)%lmax,lmax)
    end do

    ! Initialise Gaunt Coefficients
    call gaunt_init(lmax+1)

    ! Set module-level variables
    max_core_wf_tot = maxval(paw_sp(:)%n_core_wfs_tot) ! lr408
    max_paw_proj_tot = maxval(paw_sp(:)%npw_tot)
    if (all(paw_sp(:)%vhntzc_format==2)) then
       pub_nhat_in_xc = .false.
    else if (all(paw_sp(:)%vhntzc_format==1)) then
       pub_nhat_in_xc = .true.
    else
       call utils_abort('Error in paw_read_species: vhntzc_format must be &
            &same for all datasets')
    end if

    ! No NLCC unless we find a species with a nonzero PS core charge
    pub_nlcc = .false.

    ! Calculate all the pre-calculated information required to perform
    ! a PAW calculation, for each species
    do isp=1,nsp
       call paw_dataset_init(paw_sp(isp),paw_sp(isp)%core_wvfns_exist)
    end do

    ! Charge augmentation will always be active in PAW calculations
    pub_aug = .true.

    ! USP and PAW cannot be simultaneously active
    pub_usp = .false.

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving paw_read_species'
#endif

contains

    subroutine internal_print_datasets

      !=========================================================-!
      ! This subroutine prints out the number of atoms, NGWFs    !
      ! and partial waves for each species.                      !
      !----------------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 25/03/2007.          !
      ! Adapted for PAW by Nicholas Hine on 20/05/2010.          !
      !==========================================================!

      implicit none

      ! ndmh: Local Variables
      integer :: n_atoms
      integer :: n_ngwfs
      integer :: n_pws
      integer :: species_counter
      integer :: row
      integer :: ipw
      character(len=2) :: el_symbol

      character(len=64)    :: current_file


      write(stdout,'(a)') '<<<<<<<<<<<<<<<<<<<<<<<<<<< &
           &PAW Dataset information >>>>>>>>>>>>>>>>>>>>>>>>>>>'

      do isp=1,pub_cell%num_pspecies

         ! ndmh: get file name for this species
         current_file = paw_sp(isp)%ds_name

         ! ndmh: print basic information about this pseudopotential
         write(stdout,'(3a,f5.2,a,i3,a,f5.2,a)') 'File: ',trim(current_file), &
                 '; rc =',paw_sp(isp)%rcut,' bohr; shape type =', &
                 paw_sp(isp)%shape%shape_type, '; rshape =', &
                 paw_sp(isp)%shape%rshape,' bohr;'
         write(stdout,'(a,i3,a,f10.6)') '  Atomic number:', &
              paw_sp(isp)%atomic_number,';  ionic charge:', &
              paw_sp(isp)%ion_charge
         do ipw=1,paw_sp(isp)%npw
            write(stdout,'(2(a,i2))') '    Partial Wave',ipw,': l =', &
                 paw_sp(isp)%l_pw(ipw)
         end do
         if (paw_sp(isp)%tcore_charge) write(stdout,'(a)') &
              '  Core charge supplied for Nonlinear Core Corrections'

      end do

      write(stdout,'(a/)') '<<<<<<<<<<<<<<<<<<<<<<<<<<&
           &<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

    end subroutine internal_print_datasets

  end subroutine paw_read_species


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_read_dataset(species)

    !==================================================================!
    ! This subroutine reads one PAW dataset file into the PAW_SPECIES  !
    ! type 'species'.                                                  !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  species (inout)  : Dataset for the PAW species to load.         !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 17/05/2010.                          !
    !==================================================================!

    use utils, only: utils_abort, utils_unit, utils_alloc_check, &
         utils_dealloc_check, utils_open_unit_check, utils_close_unit_check

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(inout) :: species

    ! Local Variables
    character(len=100) :: dummy,dummy2
    integer :: ierr
    integer :: iunit
    integer :: i

    integer :: ipw, ipw_tot
    integer :: mi,li

    ! data to read in
    character(len=32) :: pspfmt
    real(kind=DP) :: r2well
    real(kind=DP) :: atnum_real
    integer :: pspdat, pspcod, pspxc, lmax, lloc, mmax
    integer :: creatorID
    integer :: igrid, jgrid

    ! ndmh: open the file
    iunit = utils_unit()
    open(iunit,file=trim(species%ds_name),status='old',position='rewind', &
         iostat=ierr)
    call utils_open_unit_check('paw_read_dataset',trim(species%ds_name),ierr)

#ifdef DEBUG
    write(stdout,'(a)') 'DEBUG: Entering paw_read_dataset'
#endif

    ! Read comment and header
    read(iunit,*) dummy
    read(iunit,*) atnum_real,species%ion_charge,pspdat
    species%atomic_number = int(atnum_real)
    read(iunit,*) pspcod,pspxc,lmax,lloc,mmax,r2well
    read(iunit,*) pspfmt,creatorID
    read(iunit,*) species%npw,species%npw_tot

    ! Allocate arrays to store angular momentum values of each partial wave
    ! numbered either by ni,li or by ni,li,mi. Also store the ni,li value
    ! corresponding to each ni,li,mi value (in ipwtot)
    allocate(species%l_pw(species%npw),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%l_pw',ierr)
    allocate(species%ipw_tot(species%npw_tot),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%ipw_tot',ierr)
    allocate(species%l_pw_tot(species%npw_tot),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%l_pw_tot',ierr)
    allocate(species%m_pw_tot(species%npw_tot),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%m_pw_tot',ierr)

    read(iunit,*,iostat=ierr) (species%l_pw(i),i=1,species%npw)
    species%lmax = maxval(species%l_pw(:))
    species%lmax = max(species%lmax*2,species%lmax+1)

    ! Fill in values of ipw_tot, l_pw_tot, m_pw_tot
    ipw_tot = 1
    do ipw=1,species%npw
       li = species%l_pw(ipw)
       do mi=-li,li
          species%ipw_tot(ipw_tot) = ipw
          species%l_pw_tot(ipw_tot) = li
          species%m_pw_tot(ipw_tot) = mi
          ipw_tot = ipw_tot + 1
       end do
    end do

    ! Read number of grids
    read(iunit,*) species%ngrid
    allocate(species%grid(species%ngrid),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%grid',ierr)

    ! Read in and initialise the grids
    call paw_read_grids(species%ngrid,species%grid,iunit)

    ! Read augmentation region radius
    read(iunit,*) species%rcut

    ! Read compensation density shape type and radius
    read(iunit,*) species%shape%shape_type,species%shape%rshape
    if (species%shape%rshape<1e-8_DP) species%shape%rshape = species%rcut

    ! Read the AE partial Waves
    do ipw=1,species%npw

       read(iunit,*) dummy, dummy2
       if (index(dummy2,'PHI')==0) then
          call utils_abort('Error reading dataset '//trim(species%ds_name) &
               //': Expected "PHI"')
       end if

       read(iunit,*) igrid
       if (ipw==1) then
          species%phi_grid = igrid
          allocate(species%phi_rad(species%grid(igrid)%npt,species%npw), &
               stat=ierr)
          call utils_alloc_check('paw_read_dataset','species%phi_rad',ierr)
          allocate(species%tphi_rad(species%grid(igrid)%npt,species%npw), &
               stat=ierr)
          call utils_alloc_check('paw_read_dataset','species%tphi_rad',ierr)
       else
          if (igrid/=species%phi_grid) then
             call utils_abort('Error in paw_read_dataset: Non-matching grids &
                  &for partial waves in '//trim(species%ds_name))
          end if
       end if

       read(iunit,*) (species%phi_rad(i,ipw),i=1,species%grid(igrid)%npt)

    end do

    ! Read the PS partial waves
    do ipw=1,species%npw

       read(iunit,*) dummy, dummy2
       if (index(dummy2,'TPHI')==0) then
          call utils_abort('Error reading dataset '//trim(species%ds_name) &
               //': Expected "TPHI"')
       end if

       read(iunit,*) igrid
       if (igrid/=species%phi_grid) then
          call utils_abort('Error in paw_read_dataset: Non-matching grids &
               &for partial waves in '//trim(species%ds_name))
       end if

       read(iunit,*) (species%tphi_rad(i,ipw),i=1,species%grid(igrid)%npt)

    end do

    ! Read the projectors
    species%proj_rcut = -1.0_DP
    do ipw=1,species%npw

       read(iunit,*) dummy, dummy2
       if (index(dummy2,'TPROJ')==0) then
          call utils_abort('Error reading dataset '//trim(species%ds_name) &
               //': Expected "TPROJ"')
       end if
       read(iunit,*) igrid

       if (ipw==1) then
          species%proj_grid = igrid
          allocate(species%tproj_rad(species%grid(igrid)%npt,species%npw), &
               stat=ierr)
          call utils_alloc_check('paw_read_dataset','species%tproj_rad',ierr)
       else
          if (igrid/=species%proj_grid) then
             call utils_abort('Error in paw_read_dataset: Non-matching grids &
                  &for projectors in '//trim(species%ds_name))
          end if
       end if

       read(iunit,*) (species%tproj_rad(i,ipw),i=1,species%grid(igrid)%npt)

    end do
    species%proj_rcut = species%rcut

    ! Find the AE Core density
    read(iunit,*) dummy, dummy2
    if (index(dummy2,'CORE_DENSITY')==0) then
       call utils_abort('Error reading dataset '//trim(species%ds_name) &
            //': Expected "CORE_DENSITY"')
    end if
    read(iunit,*) igrid
    species%core_den_grid = igrid

    ! Check that the Core density grid has more or equal points to the phi grid
    ! and that it has the same spacings
    jgrid = species%phi_grid
    if (species%grid(igrid)%npt<species%grid(jgrid)%npt) then
       call utils_abort('Error in paw_read_dataset: Core density grid has &
            &fewer points than partial wave grid')
    end if
    if (species%grid(igrid)%rad_step/=species%grid(jgrid)%rad_step) then
       call utils_abort('Error in paw_read_dataset: Core density grid has &
            &different spacings from partial wave grid')
    end if
    if (species%grid(igrid)%log_step/=species%grid(jgrid)%log_step) then
       call utils_abort('Error in paw_read_dataset: Core density grid has &
            &different spacings from partial wave grid')
    end if

    ! Allocate storage for core densities
    allocate(species%core_den_rad(species%grid(igrid)%npt),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%core_den_rad',ierr)
    allocate(species%tcore_den_rad(species%grid(igrid)%npt),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%tcore_den_rad',ierr)
    read(iunit,*) (species%core_den_rad(i),i=1,species%grid(igrid)%npt)

    ! Read the Pseudo Core density
    read(iunit,*) dummy, dummy2
    if (index(dummy2,'CORE_DENSITY')==0) then
       call utils_abort('Error reading dataset '//trim(species%ds_name) &
            //': Expected "PSEUDO_CORE_DENSITY" or "TCORE_DENSITY"')
    end if
    read(iunit,*) igrid
    if (igrid/=species%core_den_grid) then
       call utils_abort('Error in paw_read_dataset: Non-matching grids &
            &for core densities in '//trim(species%ds_name))
    end if
    read(iunit,*) (species%tcore_den_rad(i),i=1,species%grid(igrid)%npt)
    if (any(species%tcore_den_rad/=0.0_DP)) then
       species%tcore_charge = .true.
    else
       species%tcore_charge = .false.
    end if

    ! Read Dij0
    allocate(species%dij0(species%npw_tot,species%npw_tot),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%dij0',ierr)
    species%dij0 = 0.0_DP
    read(iunit,*) dummy, dummy2
    if (index(dummy2,'Dij0')==0) then
       call utils_abort('Error reading dataset '//trim(species%ds_name) &
            //': Expected "Dij0"')
    end if
    do ipw=1,species%npw_tot
       read(iunit,*) (species%dij0(i,ipw),i=1,ipw)
    end do
    do ipw=1,species%npw_tot
       do i=ipw+1,species%npw_tot
          species%dij0(i,ipw) = species%dij0(ipw,i)
       end do
    end do

    ! Read rhoij0
    allocate(species%rhoij0(species%npw_tot,species%npw_tot),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%rhoij0',ierr)
    species%rhoij0 = 0.0_DP
    read(iunit,*) dummy, dummy2
    if (index(dummy2,'Rhoij0')==0) then
       call utils_abort('Error reading dataset '//trim(species%ds_name) &
            //': Expected "Rhoij0"')
    end if
    do ipw=1,species%npw_tot
       read(iunit,*) (species%rhoij0(i,ipw),i=1,ipw)
    end do
    do ipw=1,species%npw_tot
       do i=ipw+1,species%npw_tot
          species%rhoij0(i,ipw) = species%rhoij0(ipw,i)
       end do
    end do

    ! Read the vhntzc potential
    read(iunit,*) dummy, dummy2
    if (index(dummy2,'VHntZC')==0) then
       call utils_abort('Error reading dataset '//trim(species%ds_name) &
            //': Expected "VHntZC"')
    end if
    read(iunit,'(a)') dummy
    if (index(pspfmt,'paw5')>0) then
       read(dummy,*) igrid, species%vhntzc_format
       if ((species%vhntzc_format<1).or.(species%vhntzc_format>2)) then
          call utils_abort('Error reading dataset '//trim(species%ds_name) &
            //': Expected "vhntzc format" to be 1 or 2')
       end if
    else
       species%vhntzc_format = 1
       read(dummy,*) igrid
    end if
    species%vhntzc_grid = igrid
    allocate(species%vhntzc_rad(species%grid(igrid)%npt),stat=ierr)
    call utils_alloc_check('paw_read_dataset','species%vhntzc_rad',ierr)
    read(iunit,*) (species%vhntzc_rad(i),i=1,species%grid(igrid)%npt)

    ! ndmh: close the file
    close(iunit,iostat=ierr)
    call utils_close_unit_check('paw_read_dataset',trim(species%ds_name),ierr)

#ifdef DEBUG
    write(stdout,'(a)') 'DEBUG: Leaving paw_read_dataset'
#endif

  end subroutine paw_read_dataset


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_read_core_wfs(species)

    !==================================================================!
    ! This subroutine reads one PAW dataset file into the PAW_SPECIES  !
    ! type 'species'.                                                  !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  species (inout)  : Dataset for the PAW species to load.         !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 17/05/2010.                          !
    !==================================================================!

    use utils, only: utils_abort, utils_unit, utils_alloc_check, &
         utils_dealloc_check, utils_open_unit_check, utils_close_unit_check

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(inout) :: species

    ! Local Variables
    character(len=100) :: dummy,dummy2
    integer :: ierr
    integer :: iunit
    integer :: i

    integer :: icore_wf, icore_wf_tot
    integer :: mi,li

    ! data to read in
    character(len=32) :: pspfmt
    integer :: method,nspinor,nsppol
    integer :: pspdat, pspcod, pspxc, lmax, lloc
    integer :: creatorID
    integer :: igrid
    real(kind=DP) :: atnum_real

    ! ndmh: open the file
    iunit = utils_unit()
    open(iunit,file=trim(species%core_wf_name),status='old',position='rewind', &
         iostat=ierr)
    call utils_open_unit_check('paw_read_core_wfs', &
         trim(species%core_wf_name),ierr)

#ifdef DEBUG
    write(stdout,'(a)') 'DEBUG: Entering paw_read_core_wfs'
#endif

    ! Read comment and header
    read(iunit,*) dummy
    read(iunit,*) method,nspinor,nsppol
    read(iunit,*) atnum_real,species%core_wf_charge,pspdat
    if (species%atomic_number /= int(atnum_real)) then
       call utils_abort('Error in paw_read_core_wfs: Atomic number of core &
            &wvfn file does not match PAW dataset')
    end if
    read(iunit,*) pspcod,pspxc,lmax
    read(iunit,*) pspfmt,creatorID
    read(iunit,*) species%n_core_wfs,species%n_core_wfs_tot

    allocate(species%l_core_wf(species%n_core_wfs),stat=ierr)
    call utils_alloc_check('paw_read_core_wfs','species%l_core_wf',ierr)
    allocate(species%l_core_wf_tot(species%n_core_wfs_tot),stat=ierr)
    call utils_alloc_check('paw_read_core_wfs','species%l_core_wf_tot',ierr)
    allocate(species%m_core_wf_tot(species%n_core_wfs_tot),stat=ierr)
    call utils_alloc_check('paw_read_core_wfs','species%m_core_wf_tot',ierr)
    allocate(species%icore_wf_tot(species%n_core_wfs_tot),stat=ierr)
    call utils_alloc_check('paw_read_core_wfs','species%icore_wf_tot',ierr)
    allocate(species%n_core_wf(species%n_core_wfs),stat=ierr)
    call utils_alloc_check('paw_read_core_wfs','species%n_core_wf',ierr)

    read(iunit,*,iostat=ierr) (species%l_core_wf(i),i=1,species%n_core_wfs)

    ! Fill in values of icore_wf_tot, l_core_wf_tot, m_core_wf_tot
    icore_wf_tot = 1
    do icore_wf=1,species%n_core_wfs
       li = species%l_core_wf(icore_wf)
       do mi=-li,li
          species%icore_wf_tot(icore_wf_tot) = icore_wf
          species%l_core_wf_tot(icore_wf_tot) = li
          species%m_core_wf_tot(icore_wf_tot) = mi
          icore_wf_tot = icore_wf_tot + 1
       end do
    end do

    if (icore_wf_tot-1/=species%n_core_wfs_tot) then
       call utils_abort('Error reading dataset '//trim(species%core_wf_name) &
            //': Core wavefunction total counts do not match')
    end if

    ! Read number of grids
    read(iunit,*) species%ngrid_core
    allocate(species%grid_core(species%ngrid_core),stat=ierr)
    call utils_alloc_check('paw_read_core_wfs','species%grid_core',ierr)

    ! Read and initialise the grids
    call paw_read_grids(species%ngrid_core,species%grid_core,iunit)

    ! Read augmentation region radius
    read(iunit,*) species%rcut_core

    ! Read the Core wavefunctions
    species%core_wf_grid = 0
    do icore_wf=1,species%n_core_wfs

       read(iunit,*) dummy, dummy2
       if (index(dummy2,'Core')==0) then
          call utils_abort('Error reading dataset '//trim(species%core_wf_name) &
               //': Expected "Core"')
       end if

       ! Read in grid and check it matches previous grids
       read(iunit,*) igrid
       if (icore_wf==1) then
          species%core_wf_grid = igrid
          allocate(species%core_wf_rad(species%grid_core(igrid)%npt, &
               species%n_core_wfs),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%core_wf_rad',ierr)
          allocate(species%core_wf_eig(species%n_core_wfs),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%core_wf_eig',ierr)
          allocate(species%core_wf_occ(species%n_core_wfs),stat=ierr)
          call utils_alloc_check('paw_read_core_wfs','species%core_wf_occ',ierr)
       else
          if (igrid/=species%core_wf_grid) then
             call utils_abort('Error in paw_read_core_wfs: Non-matching grids &
                  &for core wavefunctions in '//trim(species%core_wf_name))
          end if
       end if

       read(iunit,*) species%n_core_wf(icore_wf), lloc, nsppol
       if (lloc/=species%l_core_wf(icore_wf)) then
          call utils_abort('Error reading dataset '//trim(species%ds_name) &
               //': Core wavefunction angular momenta do not match')
       end if
       if (nsppol/=1) then
          call utils_abort('Error reading dataset '//trim(species%ds_name) &
               //': Unexpected value of spin')
       end if

       read(iunit,*) species%core_wf_eig(icore_wf),species%core_wf_occ(icore_wf)

       read(iunit,*) (species%core_wf_rad(i,icore_wf),i=1, &
            species%grid_core(igrid)%npt)

    end do

    ! ndmh: close the file
    close(iunit,iostat=ierr)
    call utils_close_unit_check('paw_read_core_wfs', &
         trim(species%core_wf_name),ierr)

#ifdef DEBUG
    write(stdout,'(a)') 'DEBUG: Leaving paw_read_core_wfs'
#endif

  end subroutine paw_read_core_wfs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_read_grids(ngrid,grids,iunit)

    !==================================================================!
    ! This subroutine reads and initialises grids from PAW datasets    !
    ! and core wavefunction files.                                     !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  ngrid (in)     : Number of grids to be read                     !
    !  grids (inout)  : RADIAL_GRID type (already allocated) to create !
    !  iunit (in)     : Unit number to read from                       !
    !------------------------------------------------------------------!
    ! Originally ritten by Nicholas Hine on 17/05/2010 as part of      !
    ! paw_read_species.                                                !
    ! Moved to its own routine on 25/10/2011 so that it can also be    !
    ! used by paw_read_core_wfs                                        !
    !==================================================================!

    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: ngrid
    type(RADIAL_GRID), intent(inout) :: grids(ngrid)
    integer, intent(in) :: iunit

    ! Local Variables
    integer :: grid_type
    integer :: igrid,jgrid
    integer :: ierr
    integer :: i
    character(len=100) :: dummy

    ! Read in and initialise the grids
    do igrid=1,ngrid

       read(iunit,'(a)') dummy
       read(dummy,*) jgrid, grid_type

       select case (grid_type)
       case (1)
          grids(igrid)%log_step = 0.0_DP
          read(dummy,*) jgrid, grid_type, grids(igrid)%npt, &
               grids(igrid)%rad_step
       case (2)
          read(dummy,*) jgrid, grid_type, grids(igrid)%npt, &
               grids(igrid)%rad_step, grids(igrid)%log_step
       case default
          call utils_abort('Error in paw_read_grids: Unsupported grid type')
       end select

       allocate(grids(igrid)%r(grids(igrid)%npt),stat=ierr)
       call utils_alloc_check('paw_read_grids','grids(igrid)%r',ierr)
       allocate(grids(igrid)%rab(grids(igrid)%npt),stat=ierr)
       call utils_alloc_check('paw_read_grids','grids(igrid)%rab', &
            ierr)

       select case (grid_type)
       case (1)
          ! Regular grid
          do i=1,grids(igrid)%npt
             grids(igrid)%r(i) = real(i-1,kind=DP)*grids(igrid)%rad_step
             grids(igrid)%rab(i) = grids(igrid)%rad_step
          end do
       case (2)
          ! Logarithmic grid
          do i=1,grids(igrid)%npt
             grids(igrid)%r(i) = (exp(grids(igrid)%log_step &
                  *real(i-1,kind=DP))-1.0_DP)*grids(igrid)%rad_step
             grids(igrid)%rab(i) = grids(igrid)%log_step* &
                  (grids(igrid)%r(i)+grids(igrid)%rad_step)
          end do
       case default
          call utils_abort('Error in paw_read_grids: Unsupported grid type')
       end select

    end do

  end subroutine paw_read_grids


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_bcast_dataset(species)

    !==================================================================!
    ! This subroutine broadcasts the contents of one PAW_SPECIES type  !
    ! from the root node (which read it in) to all the other nodes.    !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  species (inout) : Dataset for the PAW species to share between  !
    !                    nodes.                                        !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 18/05/2010.                          !
    !==================================================================!

    use comms, only: comms_bcast, pub_root_node_id, pub_on_root, pub_my_node_id
    use utils, only: utils_abort, utils_unit, utils_alloc_check, &
         utils_dealloc_check, utils_open_unit_check, utils_close_unit_check

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(inout) :: species

    ! Local Variables
    integer :: ierr
    integer :: npt

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering paw_bcast_dataset'
#endif

    call comms_bcast(pub_root_node_id,species%atomic_number)
    call comms_bcast(pub_root_node_id,species%ion_charge)
    call comms_bcast(pub_root_node_id,species%npw)
    call comms_bcast(pub_root_node_id,species%npw_tot)

    if (.not.pub_on_root) then
       allocate(species%l_pw(species%npw),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%l_pw',ierr)
       allocate(species%ipw_tot(species%npw_tot),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%ipw_tot',ierr)
       allocate(species%l_pw_tot(species%npw_tot),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%l_pw_tot',ierr)
       allocate(species%m_pw_tot(species%npw_tot),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%m_pw_tot',ierr)
    end if
    call comms_bcast(pub_root_node_id,species%l_pw,species%npw)
    call comms_bcast(pub_root_node_id,species%ipw_tot,species%npw_tot)
    call comms_bcast(pub_root_node_id,species%l_pw_tot,species%npw_tot)
    call comms_bcast(pub_root_node_id,species%m_pw_tot,species%npw_tot)
    call comms_bcast(pub_root_node_id,species%lmax)

    call comms_bcast(pub_root_node_id,species%ngrid)
    if (.not.pub_on_root) then
       allocate(species%grid(species%ngrid),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%grid',ierr)
    end if
    call paw_bcast_grids(species%ngrid,species%grid)

    call comms_bcast(pub_root_node_id,species%rcut)
    call comms_bcast(pub_root_node_id,species%proj_rcut)
    call comms_bcast(pub_root_node_id,species%shape%shape_type)
    call comms_bcast(pub_root_node_id,species%shape%rshape)

    call comms_bcast(pub_root_node_id,species%phi_grid)
    call comms_bcast(pub_root_node_id,species%proj_grid)
    call comms_bcast(pub_root_node_id,species%core_den_grid)
    call comms_bcast(pub_root_node_id,species%vhntzc_grid)

    if (.not.pub_on_root) then
       allocate(species%phi_rad(species%grid(species%phi_grid)%npt, &
            species%npw),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%phi_rad',ierr)
       allocate(species%tphi_rad(species%grid(species%phi_grid)%npt, &
            species%npw),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%tphi_rad',ierr)
       allocate(species%tproj_rad(species%grid(species%proj_grid)%npt, &
            species%npw),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%tproj_rad',ierr)
       allocate(species%core_den_rad(species%grid(species%core_den_grid)%npt), &
            stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%core_den_rad',ierr)
       allocate(species%tcore_den_rad(species%grid(species%core_den_grid)%npt),&
            stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%tcore_den_rad',ierr)
    end if
    call comms_bcast(pub_root_node_id,species%phi_rad)
    call comms_bcast(pub_root_node_id,species%tphi_rad)
    call comms_bcast(pub_root_node_id,species%tproj_rad)
    call comms_bcast(pub_root_node_id,species%core_den_rad)
    call comms_bcast(pub_root_node_id,species%tcore_den_rad)
    call comms_bcast(pub_root_node_id,species%tcore_charge)

    if (.not.pub_on_root) then
       allocate(species%dij0(species%npw_tot,species%npw_tot),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%dij0',ierr)
       allocate(species%rhoij0(species%npw_tot,species%npw_tot),stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%rhoij0',ierr)
    end if
    call comms_bcast(pub_root_node_id,species%dij0)
    call comms_bcast(pub_root_node_id,species%rhoij0)

    call comms_bcast(pub_root_node_id,species%vhntzc_format)
    if (.not.pub_on_root) then
       allocate(species%vhntzc_rad(species%grid(species%vhntzc_grid)%npt), &
            stat=ierr)
       call utils_alloc_check('paw_read_dataset','species%vhntzc_rad',ierr)
    end if
    call comms_bcast(pub_root_node_id,species%vhntzc_rad)

    ! Broadcast core WF info
    call comms_bcast(pub_root_node_id,species%core_wvfns_exist) ! lr408
    call comms_bcast(pub_root_node_id,species%n_core_wfs)
    call comms_bcast(pub_root_node_id,species%n_core_wfs_tot)

    if ((.not.pub_on_root).or.(species%n_core_wfs==0)) then
       allocate(species%l_core_wf(species%n_core_wfs),stat=ierr)
       call utils_alloc_check('paw_read_core_wfs','species%l_core_wf',ierr)
       allocate(species%l_core_wf_tot(species%n_core_wfs_tot),stat=ierr)
       call utils_alloc_check('paw_read_core_wfs','species%l_core_wf_tot', ierr)
       allocate(species%m_core_wf_tot(species%n_core_wfs_tot),stat=ierr)
       call utils_alloc_check('paw_read_core_wfs','species%m_core_wf_tot',ierr)
       allocate(species%icore_wf_tot(species%n_core_wfs_tot),stat=ierr)
       call utils_alloc_check('paw_read_core_wfs','species%icore_wf_tot',ierr)
       allocate(species%n_core_wf(species%n_core_wfs),stat=ierr)
       call utils_alloc_check('paw_read_core_wfs','species%n_core_wf',ierr)
    end if

    call comms_bcast(pub_root_node_id,species%l_core_wf)
    call comms_bcast(pub_root_node_id,species%l_core_wf_tot)
    call comms_bcast(pub_root_node_id,species%m_core_wf_tot)
    call comms_bcast(pub_root_node_id,species%icore_wf_tot)
    call comms_bcast(pub_root_node_id,species%n_core_wf)
    call comms_bcast(pub_root_node_id,species%ngrid_core)
    if ((.not.pub_on_root).or.(species%ngrid_core==0)) then
       allocate(species%grid_core(species%ngrid_core),stat=ierr)
       call utils_alloc_check('paw_read_core_wfs','species%grid_core',ierr)
    end if
    call paw_bcast_grids(species%ngrid_core,species%grid_core)
    call comms_bcast(pub_root_node_id,species%core_wf_grid)
    call comms_bcast(pub_root_node_id,species%rcut_core)

    if ((.not.pub_on_root).or.(species%n_core_wfs==0)) then
       if (species%core_wf_grid>0) then
          npt = species%grid_core(species%core_wf_grid)%npt
       else
          npt = 0
       end if
       allocate(species%core_wf_rad(npt,species%n_core_wfs),stat=ierr)
       call utils_alloc_check('paw_read_core_wfs','species%core_wf_rad',ierr)
       allocate(species%core_wf_eig(species%n_core_wfs),stat=ierr)
       call utils_alloc_check('paw_read_core_wfs','species%core_wf_eig', ierr)
       allocate(species%core_wf_occ(species%n_core_wfs),stat=ierr)
       call utils_alloc_check('paw_read_core_wfs','species%core_wf_occ',ierr)
    end if
    call comms_bcast(pub_root_node_id,species%core_wf_rad)
    call comms_bcast(pub_root_node_id,species%core_wf_eig)
    call comms_bcast(pub_root_node_id,species%core_wf_occ)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving paw_bcast_dataset'
#endif

  end subroutine paw_bcast_dataset


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_bcast_grids(ngrid,grids)

    !==================================================================!
    ! This subroutine broadcasts grids from PAW datasets to all nodes. !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  ngrid (in)     : Number of grids to be bcast                    !
    !  grids (inout)  : RADIAL_GRID type (already allocated) to bcast. !
    !------------------------------------------------------------------!
    ! Originally ritten by Nicholas Hine on 17/05/2010 as part of      !
    ! paw_bcast_dataset.                                               !
    ! Moved to its own routine on 25/10/2011 so that it can also be    !
    ! used for core wf grids .                                         !
    !==================================================================!

    use comms, only: comms_bcast, pub_on_root, pub_root_node_id
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: ngrid
    type(RADIAL_GRID), intent(inout) :: grids(ngrid)

    ! Local Variables
    integer :: igrid
    integer :: ierr

    do igrid=1,ngrid
       call comms_bcast(pub_root_node_id,grids(igrid)%npt)
       call comms_bcast(pub_root_node_id,grids(igrid)%rad_step)
       call comms_bcast(pub_root_node_id,grids(igrid)%log_step)

       if (.not.pub_on_root) then
          allocate(grids(igrid)%r(grids(igrid)%npt),stat=ierr)
          call utils_alloc_check('paw_read_dataset','grids(igrid)%r', &
               ierr)
          allocate(grids(igrid)%rab(grids(igrid)%npt),stat=ierr)
          call utils_alloc_check('paw_read_dataset','grids(igrid)%rab', &
               ierr)
       end if
       call comms_bcast(pub_root_node_id,grids(igrid)%r)
       call comms_bcast(pub_root_node_id,grids(igrid)%rab)
    end do

  end subroutine paw_bcast_grids


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dataset_init(species,core_wvfns_exist)

    !==================================================================!
    ! This subroutine allocates and initialises those arrays in the    !
    ! PAW_SPECIES type which are not loaded in directly but must be    !
    ! calculated, but which are nevertheless independent of any        !
    ! quantities depending on the system.                              !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  species (inout) : Dataset for the PAW species to calculate all  !
    !                    quantities required which are not stored in   !
    !                    the dataset file.                             !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 19/05/2010.                          !
    !==================================================================!

    use comms, only: comms_abort, comms_barrier, pub_on_root
    use constants, only: ANGSTROM
    use gaunt_coeff, only: realgaunt
    use paw_shape, only: paw_shape_init
    use services, only: services_radial_transform, services_locate_interp, &
         services_linear_interpolation, services_radial_integral, &
         services_radial_integral_rmax, services_radial_derivative
    use rundat, only: pub_nlcc
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_erf

    implicit none

    ! Arguments
    type(PAW_SPECIES), intent(inout) :: species
    logical, intent(in) :: core_wvfns_exist ! lr408

    ! Local Variables
    integer :: ierr
    integer :: igrid
    integer :: iq,ir,ir_c
    integer :: ipw,jpw,kpw,lpw,ipwtot,jpwtot,kpwtot,lpwtot
    integer :: li,mi,lj,mj,lk,mk,ll,ml
    integer :: lup, mup
    integer :: npts
    integer :: nptsc, igridc ! lr408
    real(kind=DP) :: q, r
    real(kind=DP) :: rcmax, r2new
    real(kind=DP) :: rgij, rgkl
    real(kind=DP) :: int1, int2
    real(kind=DP) :: lfac
    real(kind=DP), allocatable :: work(:),work2(:),work3(:)
    real(kind=DP), allocatable :: inter(:),inter2(:)
    real(kind=DP), allocatable :: phir_phjr(:),tphir_tphjr(:)
    real(kind=DP), allocatable :: phkrp_phlrp(:),tphkrp_tphlrp(:)
    real(kind=DP), allocatable :: rwork(:,:)

    ! ndmh: temporary arrays
    real(kind=DP), allocatable :: v_shape_L(:,:)
    real(kind=DP), allocatable :: int_v_L_g_L(:)
    real(kind=DP), allocatable :: Vhat_L_ij(:,:,:)
    real(kind=DP), allocatable :: V_L_ijkl(:,:,:,:,:)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering paw_dataset_init'
#endif

    ! Start timer
    call timer_clock('paw_dataset_init',1)

    ! Prepare reciprocal space grid
    species%n_recip_pts = 3001
    species%g_max = 50.0_DP !13.145341380124_DP
    species%inv_g_spacing = real(species%n_recip_pts-1,kind=DP)/species%g_max

    ! Fix all grid sizes to be odd
    species%grid(:)%npt = species%grid(:)%npt + modulo(species%grid(:)%npt,2) -1

    ! Find interpolation points on phi_grid for integrations up to r_c
    igrid = species%phi_grid
    npts = species%grid(igrid)%npt
    ir_c = services_locate_interp(species%rcut,species%grid(igrid)%r,npts)
    if (ir_c > npts-3) then
       if (pub_on_root) then
          write(stdout,'(a)') 'Error in paw_dataset_init: Not enough grid &
               &points beyond PAW radius on partial'
          write(stdout,'(a)') 'wave grid for accurate integration of partial &
               &waves'
       end if
       call comms_barrier
       call comms_abort
    end if

    !=========================================================================!
    ! PREPARE SHAPE FUNCTION IN REAL SPACE                                    !
    !=========================================================================!

    ! Initialise the shape functions
    igrid = species%phi_grid
    species%shape_grid = igrid
    call paw_shape_init(species%shape,species%grid(igrid)%r, &
         species%grid(igrid)%rab,species%grid(igrid)%npt,species%lmax)

    !=========================================================================!
    ! PREPARE PROJECTORS IN RECIPROCAL SPACE                                  !
    !=========================================================================!

    allocate(species%tproj_recip(species%n_recip_pts,species%npw),stat=ierr)
    call utils_alloc_check('paw_dataset_init','species%tproj_recip',ierr)

    igrid = species%proj_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
    do ipw=1,species%npw
       call services_radial_transform(species%l_pw(ipw), 1, npts,  &
            species%grid(igrid)%r,species%grid(igrid)%rab, &
            species%n_recip_pts,species%g_max,species%tproj_rad(:,ipw), &
            species%tproj_recip(:,ipw))
    end do

    !=========================================================================!
    ! PREPARE HARTREE POTENTIAL OF CORE DENSITY n_Zc                          !
    !=========================================================================!

    allocate(species%vhntzc_recip(species%n_recip_pts),stat=ierr)
    call utils_alloc_check('paw_dataset_init','species%vhntzc_recip',ierr)

    igrid = species%vhntzc_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
    rcmax = 1.0_dp/sqrt(min(-100.0_dp/(4.0_dp*log(1.d-12)),0.7_dp))

    allocate(work(npts),stat=ierr)
    allocate(work2(npts),stat=ierr)
    do ir=2,npts
       r = species%grid(igrid)%r(ir)
       work(ir) = species%vhntzc_rad(ir)*r + &
            species%ion_charge*utils_erf(r/rcmax)
    end do
    r2new = 0.25_dp*(rcmax)**2
    species%vhntzc_recip(1) = services_radial_integral(npts, &
         species%grid(igrid)%rab,work*species%grid(igrid)%r(1:npts)) &
         + species%ion_charge*r2new
    do iq=2,species%n_recip_pts
       q = real(iq-1,dp)/real(species%n_recip_pts-1,dp)*species%g_max
       do ir=1,npts
          work2(ir)=work(ir)*sin(species%grid(igrid)%r(ir)*q)
       enddo
       species%vhntzc_recip(iq) = &
            services_radial_integral(npts,species%grid(igrid)%rab,work2)/q &
            + species%ion_charge*(1.0_dp-exp(-r2new*q*q))/(q**2)
    end do
    deallocate(work,stat=ierr)
    deallocate(work2,stat=ierr)

    !=========================================================================!
    ! PREPARE TRANSFORM OF PSEUDO CORE DENSITY \tilde{n}_c                    !
    !=========================================================================!

    if (species%tcore_charge) then

       ! Set global flag if not yet already true
       pub_nlcc = .true.

       ! Storage for reciprocal space core density
       allocate(species%tcore_den_recip(species%n_recip_pts),stat=ierr)
       call utils_alloc_check('paw_dataset_init','species%tcore_den_recip',ierr)

       igrid = species%core_den_grid
       npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1

       ! Transform the real-space core charge to reciprocal space
       call services_radial_transform(0,2,npts,species%grid(igrid)%r, &
            species%grid(igrid)%rab,species%n_recip_pts,&
            species%g_max,species%tcore_den_rad,species%tcore_den_recip)

       species%tcore_den_recip = species%tcore_den_recip * 4.0_DP * PI

    else
       species%tcore_charge = .false.
       nullify(species%tcore_den_recip)
    end if

    if (any(species%core_den_rad/=0.0_DP)) then
       species%core_charge = .true.
       species%core_charge_calculated = .false.
    else
       species%core_charge = .false.
       species%core_charge_calculated = .false.
    end if

    !=========================================================================!
    ! PREPARE n^L_{n_i l_i n_j l_j} ARRAY                                     !
    !=========================================================================!

    ! This needs to be stored in the species array, for calculating nhat later
    allocate(species%aug_nLij(species%npw,species%npw,-2:species%lmax), &
         stat=ierr)
    call utils_alloc_check('paw_dataset_init','species%aug_nLij',ierr)

    igrid = species%phi_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
    allocate(work(npts),stat=ierr)
    allocate(inter(npts),stat=ierr)
    allocate(rwork(npts,2),stat=ierr)
    ! sum runs from -1 so that we get the (phi phj-tphi tphj)/r term for
    ! use in the grad operator
    do lup=-1,species%lmax
       rwork(:,1) = species%grid(igrid)%r(:)**lup
       do ipw=1,species%npw
          do jpw=1,ipw
             do ir=1,npts
                work(ir) = (species%phi_rad(ir,ipw)*species%phi_rad(ir,jpw) - &
                     species%tphi_rad(ir,ipw)*species%tphi_rad(ir,jpw)) * &
                     rwork(ir,1)
             end do
             if (modulo(ir_c,2)==1) then
                species%aug_nLij(ipw,jpw,lup) = services_radial_integral( &
                     ir_c,species%grid(igrid)%rab(1:),work(1:))
             else
                species%aug_nLij(ipw,jpw,lup) = services_radial_integral( &
                     ir_c-1,species%grid(igrid)%rab(2:),work(2:))
             end if

             species%aug_nLij(jpw,ipw,lup) = species%aug_nLij(ipw,jpw,lup)
          end do
       end do
    end do
    ! extra set of values, stored in L=-2, for use in grad operator
    do ipw=1,species%npw
       do jpw=1,species%npw
          call services_radial_derivative(rwork(:,1),species%phi_rad(:,jpw), &
               npts,real(npts,kind=DP))
          rwork(:,1) = rwork(:,1) / species%grid(igrid)%rab(:)
          call services_radial_derivative(rwork(:,2),species%tphi_rad(:,jpw), &
               npts,real(npts,kind=DP))
          rwork(:,2) = rwork(:,2) / species%grid(igrid)%rab(:)
          do ir=1,npts
             work(ir) = (species%phi_rad(ir,ipw)*rwork(ir,1) - &
                  species%tphi_rad(ir,ipw)*rwork(ir,2))
          end do
          if (modulo(ir_c,2)==1) then
             species%aug_nLij(ipw,jpw,-2) = services_radial_integral( &
                  ir_c,species%grid(igrid)%rab(1:),work(1:))
          else
             species%aug_nLij(ipw,jpw,-2) = services_radial_integral( &
                  ir_c-1,species%grid(igrid)%rab(2:),work(2:))
          end if
       end do
    end do

    deallocate(rwork,stat=ierr)
    deallocate(inter,stat=ierr)
    deallocate(work,stat=ierr)

    ! Allocate e_ijkl now, since subsequent temporary arrays will be
    ! deallocated at the end of the routine
    allocate(species%e_ijkl(species%npw_tot,species%npw_tot,species%npw_tot, &
         species%npw_tot),stat=ierr)
    call utils_alloc_check('paw_dataset_init','species%e_ijkl',ierr)

    !=========================================================================!
    ! PREPARE v^L(r) FUNCTION                                                 !
    !=========================================================================!

    igrid = species%shape_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
    allocate(v_shape_L(npts,0:species%lmax),stat=ierr)
    call utils_alloc_check('paw_dataset_init','v_shape_L',ierr)
    allocate(work(npts),stat=ierr)
    allocate(work2(npts),stat=ierr)
    allocate(inter(npts),stat=ierr)
    allocate(inter2(npts),stat=ierr)
    allocate(rwork(npts,4),stat=ierr)
    do lup=0,species%lmax
       rwork(:,1) = species%grid(igrid)%r(:)**lup
       rwork(:,2) = 1.0_DP/(species%grid(igrid)%r(:)**(lup+1))
       rwork(:,3) = species%grid(igrid)%r(:)**(lup+2)
       rwork(:,4) = species%grid(igrid)%r(:)**(1-lup)
       if (species%grid(igrid)%r(1)==0.0_DP) then
          rwork(1,2) = 0.0_DP
          if (lup>1) rwork(1,4) = 0.0_DP
       end if
       do ir=1,npts
          r = species%grid(igrid)%r(ir)
          ! work(:)  = g_L(r) * rp^L / r ^(L+1) * rp^2
          ! work2(:) = g_L(r) * r ^L / rp^(L+1) * rp^2
          work(:)  = rwork(:,3) * species%shape%shape_rad(:,lup) * rwork(ir,2)
          work2(:) = rwork(:,4) * species%shape%shape_rad(:,lup) * rwork(ir,1)
          int1 = services_radial_integral(npts,species%grid(igrid)%rab, &
               work,inter)
          int2 = services_radial_integral(npts,species%grid(igrid)%rab, &
               work2,inter2)
          v_shape_L(ir,lup) = inter(ir) + (int2 - inter2(ir))
          !write(stdout,'(i5,7f20.12)') ir, r, species%shape%shape_rad(ir,lup), &
          !     inter(ir),int2,inter2(ir),work2(ir),v_shape_L(ir,lup)
       end do
       !write(stdout,*)
       v_shape_L(:,lup) = v_shape_L(:,lup) * 4.0_DP*PI / real(2*lup+1,kind=DP)
    end do
    deallocate(rwork,stat=ierr)
    deallocate(inter2,stat=ierr)
    deallocate(inter,stat=ierr)
    deallocate(work2,stat=ierr)
    deallocate(work,stat=ierr)

    !=========================================================================!
    ! PREPARE \int v^L(r) g_L(r) dr INTEGRALS                                 !
    !=========================================================================!

    igrid = species%shape_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
    allocate(int_v_L_g_L(0:species%lmax),stat=ierr)
    call utils_alloc_check('paw_dataset_init','int_v_L_g_L',ierr)

    allocate(work(npts),stat=ierr)
    do lup=0,species%lmax
       do ir=1,npts
          r = species%grid(igrid)%r(ir)
          work(ir) = v_shape_L(ir,lup)*species%shape%shape_rad(ir,lup)*r**2
       end do
       int_v_L_g_L(lup) = services_radial_integral(npts, &
            species%grid(igrid)%rab,work)
    end do
    deallocate(work,stat=ierr)

    !=========================================================================!
    ! PREPARE \hat{V}^L_ij ARRAY                                              !
    !=========================================================================!

    igrid = species%phi_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
    allocate(Vhat_L_ij(species%npw,species%npw,0:species%lmax),stat=ierr)
    call utils_alloc_check('paw_dataset_init','Vhat_L_ij',ierr)

    allocate(work(npts),stat=ierr)
    allocate(inter(npts),stat=ierr)
    do lup=0,species%lmax
       do ipw=1,species%npw
          do jpw=1,species%npw
             work(:) = 0.0_DP
             do ir=1,npts
                work(ir) = v_shape_L(ir,lup)*species%tphi_rad(ir,ipw)* &
                     species%tphi_rad(ir,jpw)
             end do
             Vhat_L_ij(ipw,jpw,lup) = &
                  services_radial_integral_rmax(npts,species%grid(igrid)%rab, &
                  species%grid(igrid)%r,species%rcut,work,inter)
          end do
       end do
    end do
    deallocate(inter,stat=ierr)
    deallocate(work,stat=ierr)

    !=========================================================================!
    ! PREPARE V^L_ijkl ARRAY                                                  !
    !=========================================================================!

    igrid = species%phi_grid
    npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
    allocate(V_L_ijkl(species%npw,species%npw,species%npw,species%npw, &
         0:species%lmax),stat=ierr)
    call utils_alloc_check('paw_dataset_init','V_L_ijkl',ierr)

    allocate(work(npts),stat=ierr)
    allocate(work2(npts),stat=ierr)
    allocate(work3(npts),stat=ierr)
    allocate(inter(npts),stat=ierr)
    allocate(inter2(npts),stat=ierr)
    allocate(phir_phjr(npts),stat=ierr)
    allocate(tphir_tphjr(npts),stat=ierr)
    allocate(phkrp_phlrp(npts),stat=ierr)
    allocate(tphkrp_tphlrp(npts),stat=ierr)
    allocate(rwork(npts,2),stat=ierr)

    ! Loop over L and quadruple loop over partial waves
    do lup=0,species%lmax
       lfac = 4.0_DP * PI / real(2*lup+1,kind=DP)
       rwork(:,1) = species%grid(igrid)%r(:)**lup
       rwork(:,2) = 1.0_DP / (species%grid(igrid)%r(:)**(lup+1))
       if (species%grid(igrid)%r(1)==0.0_DP) then
          rwork(1,2) = 0.0_DP
       end if
       do ipw=1,species%npw
          do jpw=1,ipw
             phir_phjr(:) = species%phi_rad(:,ipw)*species%phi_rad(:,jpw)
             tphir_tphjr(:) = species%tphi_rad(:,ipw)*species%tphi_rad(:,jpw)

             do kpw=1,species%npw
                do lpw=1,kpw
                   phkrp_phlrp(:) = species%phi_rad(:,kpw) * &
                        species%phi_rad(:,lpw)
                   tphkrp_tphlrp(:) = species%tphi_rad(:,kpw) * &
                        species%tphi_rad(:,lpw)
                   ! Loop over grid points
                   do ir=1,npts

                      inter2(:) = phir_phjr(ir) * phkrp_phlrp(:) - &
                                  tphir_tphjr(ir) * tphkrp_tphlrp(:)

                      ! work(r') = (phi_i(r).phi_j(r).phi_k(r').phi_l(r')
                      !          - tphi_i(r).tphi_j(r).tphi_k(r').tphi_l(r')) *
                      !            r'^L / r^(L+1)
                      work(:)  = inter2(:) * rwork(:,1) * rwork(ir,2)

                      ! work2(r') = (phi_i(r).phi_j(r).phi_k(r').phi_l(r')
                      !          - tphi_i(r).tphi_j(r).tphi_k(r').tphi_l(r')) *
                      !            r^L / r'^(L+1)
                      work2(:) = inter2(:) * rwork(ir,1) * rwork(:,2)

                      ! Integrate work(r) from 0 to r
                      int1 = services_radial_integral_rmax(npts, &
                           species%grid(igrid)%rab,species%grid(igrid)%r, &
                           species%grid(igrid)%r(ir),work,inter)

                      ! Integrate work2(r) from 0 to rc, storing intermediates
                      int2 = services_radial_integral_rmax(npts, &
                           species%grid(igrid)%rab,species%grid(igrid)%r, &
                           species%rcut,work2,inter2)

                      ! int1            is integral of work(r) from 0 to r
                      ! int2-inter2(ir) is integral of work2(r) from r to rc
                      work3(ir) = int1 + (int2 - inter2(ir))

                   end do

                   ! integrate work3(r) from 0 to rc
                   V_L_ijkl(ipw,jpw,kpw,lpw,lup) = lfac * &
                        services_radial_integral_rmax(npts, &
                        species%grid(igrid)%rab,species%grid(igrid)%r, &
                        species%rcut,work3,inter2)

                end do  ! lpw
             end do  ! kpw

             ! Fill in remaining values from symmetric equivalents for this ij
             do kpw=1,species%npw
                do lpw=kpw+1,species%npw
                   V_L_ijkl(ipw,jpw,kpw,lpw,lup) = V_L_ijkl(ipw,jpw,lpw,kpw,lup)
                end do  ! lpw
             end do  ! kpw

          end do  ! jpw
       end do  ! ipw
    end do  ! lup

    ! Fill in remaining values from symmetric equivalents
    do ipw=1,species%npw
       do jpw=ipw+1,species%npw
          V_L_ijkl(ipw,jpw,:,:,:) = V_L_ijkl(jpw,ipw,:,:,:)
       end do
    end do

    deallocate(rwork,stat=ierr)
    deallocate(tphkrp_tphlrp,stat=ierr)
    deallocate(phkrp_phlrp,stat=ierr)
    deallocate(tphir_tphjr,stat=ierr)
    deallocate(phir_phjr,stat=ierr)
    deallocate(inter2,stat=ierr)
    deallocate(inter,stat=ierr)
    deallocate(work3,stat=ierr)
    deallocate(work2,stat=ierr)
    deallocate(work,stat=ierr)

    !=========================================================================!
    ! PREPARE e_ijkl TENSOR                                                   !
    !=========================================================================!

    species%e_ijkl(:,:,:,:) = 0.0_DP
    do lup=0,species%lmax
       do mup=-lup,lup
          do ipwtot=1,species%npw_tot
             ipw = species%ipw_tot(ipwtot)
             li = species%l_pw_tot(ipwtot)
             mi = species%m_pw_tot(ipwtot)
             do jpwtot=1,species%npw_tot
                jpw = species%ipw_tot(jpwtot)
                lj = species%l_pw_tot(jpwtot)
                mj = species%m_pw_tot(jpwtot)
                rgij = realgaunt(lup,mup,li,mi,lj,mj)
                if (abs(rgij)<1e-16) cycle
                do kpwtot=1,species%npw_tot
                   kpw = species%ipw_tot(kpwtot)
                   lk = species%l_pw_tot(kpwtot)
                   mk = species%m_pw_tot(kpwtot)
                   do lpwtot=1,species%npw_tot
                      lpw = species%ipw_tot(lpwtot)
                      ll = species%l_pw_tot(lpwtot)
                      ml = species%m_pw_tot(lpwtot)
                      rgkl = realgaunt(lup,mup,lk,mk,ll,ml)
                      if (abs(rgkl)<1e-16) cycle
                      species%e_ijkl(ipwtot,jpwtot,kpwtot,lpwtot) = &
                           species%e_ijkl(ipwtot,jpwtot,kpwtot,lpwtot) &
                           + rgij*rgkl*(V_L_ijkl(ipw,jpw,kpw,lpw,lup) &
                           - species%aug_nLij(ipw,jpw,lup)*Vhat_L_ij(kpw,lpw,lup) &
                           - species%aug_nLij(kpw,lpw,lup)*Vhat_L_ij(ipw,jpw,lup) &
                           - species%aug_nLij(ipw,jpw,lup)*species%aug_nLij(kpw, &
                           lpw,lup)*int_v_L_g_L(lup))
#if 0
                      print '(2i3,4(i3,2i2),2f8.4,6f14.10)', &
                           lup,mup,ipw,li,mi,jpw,lj,mj,kpw,lk,mk,lpw,ll,ml, &
                           rgij,rgkl,V_L_ijkl(ipw,jpw,kpw,lpw,lup), &
                           Vhat_L_ij(kpw,lpw,lup), &
                           species%aug_nLij(ipw,jpw,lup), &
                           Vhat_L_ij(ipw,jpw,lup), &
                           species%aug_nLij(kpw,lpw,lup), &
                           int_v_L_g_L(lup)
                      print '(2i3,4(i3,2i2),2f8.4,3f14.9,2f20.12)', &
                           lup,mup,ipw,li,mi,jpw,lj,mj,kpw,lk,mk,lpw,ll,ml, &
                           rgij,rgkl,V_L_ijkl(ipw,jpw,kpw,lpw,lup), &
                           - species%aug_nLij(ipw,jpw,lup)*Vhat_L_ij(kpw,lpw,lup), &
                           - species%aug_nLij(kpw,lpw,lup)*Vhat_L_ij(ipw,jpw,lup), &
                           - species%aug_nLij(ipw,jpw,lup)*species%aug_nLij(kpw,lpw,lup) * int_v_L_g_L(lup), &
                           species%e_ijkl(ipwtot,jpwtot,kpwtot,lpwtot)
#endif
                   end do
                end do
             end do
          end do
       end do
    end do

    do ipwtot=1,species%npw_tot
       do jpwtot=1,species%npw_tot
          do kpwtot=1,species%npw_tot
             do lpwtot=1,species%npw_tot
                !print '(4i4,f20.12)',ipwtot,jpwtot,kpwtot,lpwtot,species%e_ijkl(ipwtot,jpwtot,kpwtot,lpwtot)
             end do
          end do
       end do
    end do


    ! lr408: assuming same grid here so basically not changing that much
    if (core_wvfns_exist) then

!    ! Prepare reciprocal space grid
!    species%n_recip_pts = 3001
!    species%g_max = 50.0_DP !13.145341380124_DP
!    species%inv_g_spacing = real(species%n_recip_pts-1,kind=DP)/species%g_max

       ! Fix all grid sizes to be odd
       species%grid_core(:)%npt = species%grid_core(:)%npt + modulo(species%grid_core(:)%npt,2) -1


       ! Find interpolation points on phi_grid for integrations up to r_c
       igrid = species%core_wf_grid
       npts = species%grid_core(igrid)%npt

       ir_c = services_locate_interp(species%rcut,species%grid_core(igrid)%r,npts)
       if (ir_c > npts-3) then
          if (pub_on_root) then
             write(stdout,'(a)') 'Error in paw_dataset_init: Not enough grid &
                  &points beyond PAW radius on partial'
             write(stdout,'(a)') 'wave grid for accurate integration of partial &
                  &waves'
          end if
          call comms_barrier
          call comms_abort
       end if


       !=========================================================================!
       ! PREPARE CORE WAVEFUNCTIONS IN RECIPROCAL SPACE                          !
       !=========================================================================!

       allocate(species%core_wf_recip(species%n_recip_pts,species%n_core_wfs),stat=ierr)
       call utils_alloc_check('paw_dataset_init','species%core_wf_recip',ierr)
       igrid = species%core_wf_grid
       npts = species%grid_core(igrid)%npt + modulo(species%grid_core(igrid)%npt,2) - 1
       do ipw=1,species%n_core_wfs
          call services_radial_transform(species%l_core_wf(ipw), 1, npts,  &
               species%grid_core(igrid)%r,species%grid_core(igrid)%rab, &
               species%n_recip_pts,species%g_max,species%core_wf_rad(:,ipw), &
               species%core_wf_recip(:,ipw))
       end do
    end if

    if (core_wvfns_exist) then

       !=========================================================================!
       ! PREPARE n^L_{n_i l_i n_j l_j} ARRAY FOR CORE WAVEFUNCTIONS              !
       !=========================================================================!

       ! RELIES ON BOTH GRIDS HAVING SAME SPACING - STILL NEED TO ADD ERROR MESSAGE FOR THAT!

       ! This needs to be stored in the species array, for calculating nhat later

       ! unclear what to do with lmax - hopefully the same for core wvfns but may
       ! need to return to
       allocate(species%core_aug_nLij(species%n_core_wfs, &
            species%npw,0:species%lmax),stat=ierr)
       call utils_alloc_check('paw_dataset_init','species%core_aug_nLij',ierr)

       igrid = species%phi_grid
       igridc = species%core_wf_grid
       npts = species%grid(igrid)%npt + modulo(species%grid(igrid)%npt,2) - 1
       nptsc = species%grid_core(igridc)%npt + &
            modulo(species%grid_core(igridc)%npt,2) - 1
       npts = min(npts,nptsc)

       allocate(work(npts),stat=ierr)
       allocate(inter(npts),stat=ierr)
       allocate(rwork(npts,1),stat=ierr)
       do lup=0,species%lmax
          rwork(1:npts,1) = species%grid(igrid)%r(1:npts)**lup
          do ipw=1,species%n_core_wfs
             do jpw=1,species%npw
                do ir=1,npts
                   work(ir) = (species%core_wf_rad(ir,ipw)*species%phi_rad(ir,jpw) - &
                        species%core_wf_rad(ir,ipw)*species%tphi_rad(ir,jpw)) * &
                        rwork(ir,1)
                end do
                if (modulo(ir_c,2)==1) then
                   species%core_aug_nLij(ipw,jpw,lup) = services_radial_integral( &
                        ir_c,species%grid(igrid)%rab(1:npts),work(1:npts))
                else
                   species%core_aug_nLij(ipw,jpw,lup) = services_radial_integral( &
                        ir_c-1,species%grid(igrid)%rab(2:npts),work(2:npts))
                end if

             end do
          end do
       end do
       deallocate(rwork,stat=ierr)
       deallocate(inter,stat=ierr)
       deallocate(work,stat=ierr)

    end if

    deallocate(V_L_ijkl,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','V_L_ijkl',ierr)
    deallocate(Vhat_L_ij,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','Vhat_L_ij',ierr)
    deallocate(int_v_L_g_L,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','int_v_L_g_L',ierr)
    deallocate(v_shape_L,stat=ierr)
    call utils_dealloc_check('paw_dataset_init','v_shape_L',ierr)

    ! Stop timer
    call timer_clock('paw_dataset_init',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving paw_dataset_init'
#endif

  end subroutine paw_dataset_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_species_exit

    !====================================================================!
    ! This subroutine deallocates the contents of the PAW species array, !
    ! deallocates the PAW projectors and cleans up storage for the Gaunt !
    ! Coefficients.                                                      !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  species (inout) : Dataset for the PAW species to share between    !
    !                    nodes.                                          !
    !--------------------------------------------------------------------!
    ! Written by Nicholas Hine on 18/05/2010.                            !
    !====================================================================!

    use comms, only: pub_on_root
    use gaunt_coeff, only: gaunt_exit
    use paw_shape, only: paw_shape_exit
    use projectors, only: projectors_deallocate_set
    use simulation_cell, only: pub_cell
    use utils, only: utils_dealloc_check

    implicit none

    ! Local Variables
    integer :: ierr
    integer :: isp
    integer :: igrid

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering paw_species_exit'
#endif

    ! Deallocate projector set
    call projectors_deallocate_set(paw_projectors)

    if (any(paw_sp(:)%core_wvfns_exist) .and. &
         allocated(paw_core_wvfns%fftbox_proj_recip)) then
       call projectors_deallocate_set(paw_core_wvfns)
    end if

    call gaunt_exit

    ! Deallocate reciprocal space arrays in paw_sp type
    do isp=pub_cell%num_pspecies,1,-1
       deallocate(paw_sp(isp)%e_ijkl,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%e_ijkl',ierr)
       deallocate(paw_sp(isp)%aug_nLij,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%aug_nLij',ierr)
       if (associated(paw_sp(isp)%tcore_den_recip)) then
          deallocate(paw_sp(isp)%tcore_den_recip,stat=ierr)
          call utils_dealloc_check('paw_species_exit', &
               'species%tcore_den_recip',ierr)
       end if
       deallocate(paw_sp(isp)%vhntzc_recip,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%vhntzc_recip',ierr)

       ! lr408
       if (paw_sp(isp)%core_wvfns_exist) then
          deallocate(paw_sp(isp)%core_wf_recip,stat=ierr)
          call utils_dealloc_check('paw_species_exit','species%core_wf_recip',ierr)
       end if

       deallocate(paw_sp(isp)%tproj_recip,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%tproj_recip',ierr)
    end do

    ! Initialise the shape functions
    do isp=pub_cell%num_pspecies,1,-1
       call paw_shape_exit(paw_sp(isp)%shape)
    end do

    ! Deallocate remaining arrays in paw_sp type
    do isp=pub_cell%num_pspecies,1,-1
       deallocate(paw_sp(isp)%vhntzc_rad,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%vhntzc_rad',ierr)
       deallocate(paw_sp(isp)%rhoij0,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%rhoij0',ierr)
       deallocate(paw_sp(isp)%dij0,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%dij0',ierr)
       deallocate(paw_sp(isp)%tcore_den_rad,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%tcore_den_rad',ierr)
       deallocate(paw_sp(isp)%core_den_rad,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%core_den_rad',ierr)
       deallocate(paw_sp(isp)%tproj_rad,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%tproj_rad',ierr)
       deallocate(paw_sp(isp)%tphi_rad,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%tphi_rad',ierr)
       deallocate(paw_sp(isp)%phi_rad,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%phi_rad',ierr)
       do igrid=paw_sp(isp)%ngrid,1,-1
          deallocate(paw_sp(isp)%grid(igrid)%rab,stat=ierr)
          call utils_dealloc_check('paw_species_exit', &
               'species%grid(igrid)%rab',ierr)
          deallocate(paw_sp(isp)%grid(igrid)%r,stat=ierr)
          call utils_dealloc_check('paw_species_exit', &
               'species%grid(igrid)%r',ierr)
       end do
       deallocate(paw_sp(isp)%grid,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%grid',ierr)
       deallocate(paw_sp(isp)%m_pw_tot,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%m_pw_tot',ierr)
       deallocate(paw_sp(isp)%l_pw_tot,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%l_pw_tot',ierr)
       deallocate(paw_sp(isp)%ipw_tot,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%ipw_tot',ierr)
       deallocate(paw_sp(isp)%l_pw,stat=ierr)
       call utils_dealloc_check('paw_species_exit','species%l_pw',ierr)

       ! lr408
       if (paw_sp(isp)%core_wvfns_exist) then

          deallocate(paw_sp(isp)%core_wf_occ,stat=ierr)
          call utils_dealloc_check('paw_species_exit','species%core_wf_occ',ierr)
          deallocate(paw_sp(isp)%core_wf_eig,stat=ierr)
          call utils_dealloc_check('paw_species_exit','species%core_wf_eig',ierr)
          deallocate(paw_sp(isp)%core_wf_rad,stat=ierr)
          call utils_dealloc_check('paw_species_exit','species%core_wf_rad',ierr)
          deallocate(paw_sp(isp)%grid_core,stat=ierr)
          call utils_dealloc_check('paw_species_exit','species%grid_core',ierr)

          deallocate(paw_sp(isp)%n_core_wf,stat=ierr)
          call utils_dealloc_check('paw_species_exit','species%n_core_wf',ierr)
          deallocate(paw_sp(isp)%icore_wf_tot,stat=ierr)
          call utils_dealloc_check('paw_species_exit','species%icore_wf_tot',ierr)
          deallocate(paw_sp(isp)%m_core_wf_tot,stat=ierr)
          call utils_dealloc_check('paw_species_exit','species%m_core_wf_tot',ierr)
          deallocate(paw_sp(isp)%l_core_wf_tot,stat=ierr)
          call utils_dealloc_check('paw_species_exit','species%l_core_wf_tot',ierr)
          deallocate(paw_sp(isp)%l_core_wf,stat=ierr)
          call utils_dealloc_check('paw_species_exit','species%l_core_wf',ierr)

       end if
    end do

    ! Deallocate paw_sp itself
    deallocate(paw_sp,stat=ierr)
    call utils_dealloc_check('paw_species_exit','paw_sp',ierr)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving paw_species_exit'
#endif

  end subroutine paw_species_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! SUBROUTINES FOR INTERFACING WITH THE ATOM SOLVER !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_get_locpot_rad(locpot,npts,rad,isp)

    !=====================================================================!
    ! This subroutine fetches the local pseudopotential on a regular grid !
    !---------------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                         !
    !=====================================================================!

    use comms, only: comms_abort
    use services, only: services_sbessj
    use utils, only: utils_alloc_check,utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP), intent(out) :: locpot(npts)
    real(kind=DP), intent(in) :: rad(npts)
    integer, intent(in) :: isp

    ! Local variables
    real(kind=DP), allocatable :: work(:)
    real(kind=DP) :: q,dq,Z
    integer :: npts_q
    integer :: ir,iq
    integer :: ierr

    npts_q = paw_sp(isp)%n_recip_pts
    allocate(work(npts_q),stat=ierr)
    call utils_alloc_check('paw_get_locpot_rad','work',ierr)

    dq = paw_sp(isp)%g_max / real(npts_q-1,kind=DP)
    Z = paw_sp(isp)%ion_charge

    do ir=1,npts
       do iq=1,npts_q
          q = real(iq-1,kind=DP)*dq
          work(iq) = (paw_sp(isp)%vhntzc_recip(iq)*q**2 - Z) * &
               services_sbessj(0,q*rad(ir))
       end do

       locpot(ir) = work(1)+work(npts_q)
       do iq = 2,npts_q-1,2
          locpot(ir) = locpot(ir) + 4.0_dp*work(iq)+2.0_dp*work(iq+1)
       enddo
       locpot(ir) = locpot(ir)*dq/3.0_dp
    end do

    locpot(:) = locpot(:)*(2.0_DP/PI)

    deallocate(work,stat=ierr)
    call utils_dealloc_check('paw_get_locpot_rad','work',ierr)

  end subroutine paw_get_locpot_rad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_get_core_den_rad(core_den,npts,rad,isp)

    !============================================================!
    ! This subroutine fetches the core density on a regular grid !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                !
    !============================================================!

    use comms, only: comms_abort
    use services, only: services_regular_integral, services_sbessj
    use utils, only: utils_alloc_check,utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP), intent(out) :: core_den(npts)
    real(kind=DP), intent(in) :: rad(npts)
    integer, intent(in) :: isp

    ! Local variables
    real(kind=DP), allocatable :: work(:)
    real(kind=DP) :: q,dq
    integer :: npts_q
    integer :: ir,iq
    integer :: ierr

    if (.not.paw_sp(isp)%tcore_charge) then
       core_den(:) = 0.0_DP
       return
    end if

    npts_q = paw_sp(isp)%n_recip_pts
    allocate(work(npts_q),stat=ierr)
    call utils_alloc_check('paw_get_core_den_rad','work',ierr)

    dq = paw_sp(isp)%g_max / real(npts_q-1,kind=DP)

    do ir=1,npts
       do iq=1,npts_q
          q = real(iq-1,kind=DP)*dq
          work(iq) = paw_sp(isp)%tcore_den_recip(iq)*q**2 * &
                services_sbessj(0,q*rad(ir))
       end do

       core_den(ir) = work(1) + work(npts_q)
       do iq=2,npts_q-1,2
          core_den(ir) = core_den(ir) + 4.0_DP*work(iq) + 2.0_DP*work(iq+1)
       end do
       core_den(ir) = core_den(ir)*dq/3.0_DP
    end do

    core_den(:) = core_den(:)*(2.0_DP/PI)

    deallocate(work,stat=ierr)
    call utils_dealloc_check('paw_get_core_den_rad','work',ierr)

  end subroutine paw_get_core_den_rad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_get_projector_info(isp,npw,npwtot,lmax,log_npts_max)

    !============================================================!
    ! This function fetches the number of shells of projectors.  !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                !
    !============================================================!

    ! Arguments
    integer,intent(in) :: isp
    integer,intent(out) :: npw
    integer,intent(out) :: npwtot
    integer,intent(out) :: lmax
    integer,intent(out) :: log_npts_max

    ! Local Variables
    integer :: igrid

    npw = paw_sp(isp)%npw
    npwtot = paw_sp(isp)%npw_tot
    lmax = paw_sp(isp)%lmax
    
    log_npts_max = 0
    igrid = paw_sp(isp)%phi_grid
    log_npts_max = max(log_npts_max,paw_sp(isp)%grid(igrid)%npt)
    igrid = paw_sp(isp)%shape_grid
    log_npts_max = max(log_npts_max,paw_sp(isp)%grid(igrid)%npt)
    igrid = paw_sp(isp)%core_den_grid
    log_npts_max = max(log_npts_max,paw_sp(isp)%grid(igrid)%npt)

  end subroutine paw_get_projector_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_get_projectors_q(proj_q,dij0,ang_mom,npw,nsws, &
       nsws_max,lmax,qb,isp)

    !============================================================!
    ! This subroutine fetches the projectors at chosen q-points  !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                !
    !============================================================!

    use comms, only: comms_abort
    use services, only: services_1d_interpolation
    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: npw
    integer, intent(in) :: lmax
    integer, intent(in) :: nsws(0:lmax)
    integer, intent(in) :: nsws_max
    integer, intent(in) :: isp
    integer, intent(out) :: ang_mom(npw)
    real(kind=DP), intent(in) :: qb(nsws_max,0:lmax)
    real(kind=DP), intent(out) :: proj_q(nsws_max,npw)
    real(kind=DP), intent(out) :: dij0(npw,npw)

    ! Local variables
    integer :: ipw,jpw,isw
    integer :: ipwtot,jpwtot
    integer :: li,lj,di
    real(kind=DP) :: qq

    ! Check matrices are right size
    if (npw/=paw_sp(isp)%npw) then
       call utils_abort('Error in paw_get_projectors_q: Wrong number of shells &
            &npw')
    end if

    ! Set the dij0 terms
    dij0(:,:) = 0.0_DP
    ipwtot = 1
    do ipw=1,npw
       li = paw_sp(isp)%l_pw(ipw)
       jpwtot = 1
       do jpw=1,npw
          lj = paw_sp(isp)%l_pw(jpw)
          if (li==lj) then
             do di=0,0
                dij0(ipw,jpw) = dij0(ipw,jpw) + &
                     paw_sp(isp)%dij0(ipwtot+di,jpwtot+di)
             end do
          end if
          jpwtot = jpwtot + 2*lj + 1
       end do
       ipwtot = ipwtot + 2*li + 1
    end do

    ! Interpolate the reciprocal-space projectors at the qb values
    do ipw=1,npw
       proj_q(:,ipw) = 0.0_DP
       ang_mom(ipw) = paw_sp(isp)%l_pw(ipw)
       do isw=1,nsws(ang_mom(ipw))
          qq = qb(isw,paw_sp(isp)%l_pw(ipw))
          proj_q(isw,ipw) = services_1d_interpolation( &
               paw_sp(isp)%tproj_recip(:,ipw), &
               paw_sp(isp)%n_recip_pts, &
               qq*paw_sp(isp)%inv_g_spacing, &
               paw_sp(isp)%l_pw(ipw))
       end do
    end do

  end subroutine paw_get_projectors_q


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_get_aug_funcs(qijL,shape_l,grad_shape_l,rhoij0,lmax, &
       npw,npts,rad,isp)

    !============================================================!
    ! This subroutine fetches the shape functions and the qijL   !
    ! terms for each angular momentum channel.                   !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                !
    !============================================================!

    use comms, only: comms_abort
    use gaunt_coeff, only: realgaunt
    use paw_shape, only: paw_shape_calculate
    use services, only: services_1d_interpolation
    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    integer, intent(in) :: npw
    real(kind=DP), intent(in) :: rad(npts)
    integer, intent(in) :: isp
    integer, intent(inout) :: lmax
    real(kind=DP), intent(out) :: qijL(npw,npw,0:lmax)
    real(kind=DP), intent(out) :: shape_l(npts,0:lmax)
    real(kind=DP), intent(out) :: grad_shape_l(npts,0:lmax)
    real(kind=DP), intent(out) :: rhoij0(npw,npw)

    ! Local variables
    integer :: lup, mup
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: li,lj,mi,mj

    ! Ensure the guessed lmax value is the same PAW species' version
    if (lmax < paw_sp(isp)%lmax) then
       call utils_abort('Error paw_get_aug_funcs: value of lmax provided is &
            &too small')
    else if (lmax > paw_sp(isp)%lmax) then
       lmax = paw_sp(isp)%lmax
    end if

    ! Add up qijL factors to reduce npwtot x npwtot down to npw x npw
    qijL(:,:,:) = 0.0_DP
    do lup=0,lmax
       do mup=0,0
          do ipwtot=1,paw_sp(isp)%npw_tot
             li = paw_sp(isp)%l_pw_tot(ipwtot)
             mi = paw_sp(isp)%m_pw_tot(ipwtot)
             if (mi/=0) cycle
             ipw = paw_sp(isp)%ipw_tot(ipwtot)
             do jpwtot=1,paw_sp(isp)%npw_tot
                lj = paw_sp(isp)%l_pw_tot(jpwtot)
                mj = paw_sp(isp)%m_pw_tot(jpwtot)
                if (mj/=0) cycle
                jpw = paw_sp(isp)%ipw_tot(jpwtot)
                qijL(ipw,jpw,lup) = qijL(ipw,jpw,lup) + &
                     paw_sp(isp)%aug_nLij(ipw,jpw,lup) * &
                     realgaunt(lup,mup,li,mi,lj,mj) * sqrt_4pi
             end do
          end do
       end do
    end do

    ! Evaluate the shape function on the radial grid
    do lup=0,lmax
       shape_l(:,lup) = 0.0_DP
       grad_shape_l(:,lup) = 0.0_DP
       call paw_shape_calculate(shape_l(:,lup),grad_shape_l(:,lup), &
            rad,npts,lup,paw_sp(isp)%shape%shape_type, &
            paw_sp(isp)%shape%rshape, &
            paw_sp(isp)%shape%shape_alpha(:,lup), &
            paw_sp(isp)%shape%shape_q(:,lup), &
            paw_sp(isp)%shape%shape_lambda(lup), &
            paw_sp(isp)%shape%shape_sigma(lup))
       shape_l(:,lup) = shape_l(:,lup) / paw_sp(isp)%shape%norm(lup)
       grad_shape_l(:,lup) = grad_shape_l(:,lup) / paw_sp(isp)%shape%norm(lup)
    end do

    ! Add up rhoij for ipwtot,jpwtot contributing to each ipw,jpw to get
    ! initial guess projector density matrix
    rhoij0(:,:) = 0.0_DP
    do ipwtot=1,paw_sp(isp)%npw_tot
       ipw = paw_sp(isp)%ipw_tot(ipwtot)
       do jpwtot=1,paw_sp(isp)%npw_tot
          jpw = paw_sp(isp)%ipw_tot(jpwtot)

          rhoij0(ipw,jpw) = rhoij0(ipw,jpw) + paw_sp(isp)%rhoij0(ipwtot,jpwtot)
       end do
    end do

  end subroutine paw_get_aug_funcs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_tcore_hartree_on_grid(tcore_hartree_potential,&
        struct_fac,struct_fac_classical,grid)

    !===================================================================!
    ! This subroutine generates the Hartree potential of the pseudized  !
    ! core charge \tilde{n_Zc} on the simulation cell fine grid. This   !
    ! is the PAW equivalent of the local potential.                     !
    !-------------------------------------------------------------------!
    ! Arguments:                                                        !
    !  tcore_hartree_potential (out) : V_H [n_Zc] (r) in real space     !
    !  struct_fac (in) : The structure factor for each species in       !
    !  reciprocal space.                                                !
    !  struct_fac_classical (in) : The structure factor for classical   !
    !  atoms in reciprocal space.                                       !
    !  grid (in) : The grid definition.                                 !
    !-------------------------------------------------------------------!
    ! Adapted on 20/05/2010.by Nicholas Hine from the routine           !
    ! pseudopotentials_local_on_fine, originally written by             !
    ! Chris-Kriton Skylaris on 21/2/2004 with subsequent modifications  !
    ! by Peter Haynes, Arash Mostofi and Nicholas Hine.                 !
    !===================================================================!

    use cell_grid, only: GRID_INFO
    use classical_pot, only: classical_pot_recip
    use fourier, only: fourier_apply_cell_backward
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)  :: grid
    complex(kind=DP), intent(in) :: struct_fac(pub_cell%num_pspecies, &
         grid%ld3, grid%ld2, grid%max_slabs23)
    complex(kind=DP), intent(in) :: struct_fac_classical( &
         grid%ld3, grid%ld2, grid%max_slabs23)
    real(kind=DP), intent(out) :: tcore_hartree_potential(grid%ld1, &
         grid%ld2, grid%max_slabs12)

    ! Local Variables
    integer :: ierr              ! Error flag
    complex(kind=DP), allocatable :: tcore_hartree_recip(:,:,:)

    call timer_clock('paw_tcore_hartree_on_grid',1)

    ! ndmh: expand to 3D in reciprocal fine grid and sum together
    !       (along with the structure factor) the pseudo-core Hartree
    !       potential for each species to obtain the total pseudo-core Hartree
    !       potential in reciprocal representation on the fine grid
    allocate(tcore_hartree_recip(grid%ld3,grid%ld2, &
         grid%max_slabs23),stat=ierr)
    call utils_alloc_check('paw_tcore_hartree_on_grid','tcore_hartree_recip',ierr)

    call paw_tcore_hartree_rec(tcore_hartree_recip, &      ! output
         struct_fac,grid)                                  ! input

    ! ndmh: include external potential from "classical" atoms
    if (pub_cell%nat_classical > 0) then
       call classical_pot_recip(tcore_hartree_recip, &       ! output
            struct_fac_classical,grid)                       ! input
    endif

    ! FFT the local ionic potential from reciprocal to real space
    call fourier_apply_cell_backward(tcore_hartree_potential, &
         tcore_hartree_recip,grid)

    deallocate(tcore_hartree_recip,stat=ierr)
    call utils_dealloc_check('paw_tcore_hartree_on_grid', &
         'tcore_hartree_recip',ierr)

    call timer_clock('paw_tcore_hartree_on_grid',2)

  end subroutine paw_tcore_hartree_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_tcore_hartree_rec(fine_complex,struct_fac,grid)

    !=================================================================!
    ! This subroutine generates in reciprocal space the Hartree       !
    ! potential in the simulation cell due to the pseudized core      !
    ! density \tilde{n_Zc} of all ions.                               !
    !-----------------------------------------------------------------!
    ! Arguments:                                                      !
    !  fine_complex (out) : the Hartree potential in reciprocal space !
    !  struct_fac (in) : The structure factor for each species in     !
    !  reciprocal space.                                              !
    !-----------------------------------------------------------------!
    ! Adapted on 20/05/2010.by Nicholas Hine from the routine         !
    ! pseudopotentials_local_on_fine, originally written by           !
    ! Chris-Kriton Skylaris in 2000 with subsequent modifications     !
    ! by Peter Haynes, Arash Mostofi and Nicholas Hine.               !
    !=================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: comms_abort, pub_my_node_id, pub_total_num_nodes
    use rundat, only: pub_mt_cutoff, pub_smooth_loc_pspot
    use services, only: services_1d_interpolation
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    complex(kind=DP), intent(in) :: struct_fac(pub_cell%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)
    complex(kind=DP), intent(out) :: fine_complex(grid%ld3,&
         grid%ld2, grid%max_slabs23)

    ! Local variables
    integer :: species              ! Atomic species counter
    integer :: i3,i2,islab23        ! Reciprocal grid loop counters
    real(kind=DP) :: gvec(3)        ! G vector
    real(kind=DP) :: g_length       ! Length of this G vector
    real(kind=DP) :: v_loc_value    ! Local potential at this G
    !real(kind=DP) :: g_cut, alpha   ! For filtering locps
    real(kind=DP),parameter :: fourpi = 4.0_DP * PI ! Constant multiplier

    ! Loop over reciprocal space grid on this node
    do islab23=1,grid%num_slabs23       ! along b1
       do i2=1,grid%n2                      ! along b2
          do i3=1,grid%n3                   ! along b3

             ! Get magnitude of this G-vector
             call cell_grid_recip_pt(gvec,islab23 + &
                  grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)
             g_length = sqrt(sum(gvec(1:3)**2))

             ! ndmh: initialise (cache-efficiently)
             fine_complex(i3,i2,islab23) = cmplx(0.0_DP,0.0_DP,kind=DP)

             ! Loop over atomic species
             do species=1,pub_cell%num_pspecies

                ! Get potential at this G-vector
                v_loc_value = services_1d_interpolation( &
                     paw_sp(species)%vhntzc_recip, &
                     paw_sp(species)%n_recip_pts,&
                     g_length*paw_sp(species)%inv_g_spacing,0)

                ! Add back Coulomb potential
                v_loc_value = v_loc_value - paw_sp(species)%ion_charge * &
                      grid%coulomb_recip(i3,i2,islab23)

                ! ndmh: scale by 4 pi/weight
                v_loc_value = v_loc_value * fourpi / grid%weight

                fine_complex(i3,i2,islab23) = fine_complex(i3,i2,islab23) + &
                     struct_fac(species,i3,i2,islab23) * v_loc_value

             end do    ! loop over species

          end do   ! b3
       end do      ! b2
    end do         ! b1


    ! G=0 element must be real
    if (pub_my_node_id==grid%node_slab23(1)) then
       if (aimag(fine_complex(1,1,1)) /= 0.0_DP) then
          write(stdout,'(a)') 'Error in paw_tcore_hartree_rec: &
               &potential not real'
          call comms_abort
       end if
    end if

    if (grid%num_slabs23 > 0) then
       ! Nyquist filter (fine grid is always going to be even)
       fine_complex(grid%n3/2+1,:,:) = (0.0_DP,0.0_DP)
       fine_complex(:,grid%n2/2+1,:) = (0.0_DP,0.0_DP)
    end if
    ! Nyquist filter for last slab
    if (pub_my_node_id==grid%node_slab23(grid%n1/2+1))&
         fine_complex(:,:,grid%num_slabs23) = (0.0_DP,0.0_DP)

  end subroutine paw_tcore_hartree_rec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_tcore_density(tcore_density,struct_fac,grid)

    !==================================================================!
    ! This subroutine reconstructs the core density for the whole      !
    ! supercell on the fine real space grid using the information      !
    ! stored for each dataset in the array tcore_den_recip.            !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  tcore_density (out) : data-parallelised core density            !
    !  struct_fac  (in)    : data-parallelised structure factor        !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 28/05/10.           !
    !==================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_node_id
    use fourier, only: fourier_apply_cell_backward
    use services, only: services_1d_interpolation
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(out) :: tcore_density(grid%ld1, &
         grid%ld2, grid%max_slabs12)
    complex(kind=DP), intent(in) :: struct_fac(pub_cell%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)

    ! Local variables
    integer :: ierr              ! Error flag
    complex(kind=DP), allocatable :: tcore_density_recip(:,:,:)
    real(kind=DP) :: gvec(3),g_length, tcore_den_value, factor
    integer :: species              ! Atomic species counter
    integer :: i3,i2,islab23        ! Reciprocal grid loop counters

    call timer_clock('paw_tcore_density',1)

    ! ndmh: allocate storage for core density in reciprocal space
    allocate(tcore_density_recip(grid%ld3, &
         grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('paw_tcore_density', &
         'tcore_density_recip',ierr)

    ! ndmh: loop over reciprocal space grid on this node
    do islab23=1,grid%num_slabs23           ! along b1
       do i2=1,grid%n2                      ! along b2
          do i3=1,grid%n3                   ! along b3

             ! Initialise
             tcore_density_recip(i3,i2,islab23) = (0.0_DP,0.0_DP)

             ! Loop over atomic species
             do species=1,pub_cell%num_pspecies

                ! Check if we have a core charge for this species
                if (.not.paw_sp(species)%tcore_charge) cycle

                ! Get magnitude of this G-vector
                call cell_grid_recip_pt(gvec,islab23 + &
                     grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)
                g_length = sqrt(sum(gvec(1:3)**2))

                ! Get core density at this G-vector
                tcore_den_value = services_1d_interpolation( &
                     paw_sp(species)%tcore_den_recip, &
                     paw_sp(species)%n_recip_pts,&
                     g_length*paw_sp(species)%inv_g_spacing,0)

                tcore_density_recip(i3,i2,islab23) = &
                     tcore_density_recip(i3,i2,islab23) + &
                     struct_fac(species,i3,i2,islab23) * tcore_den_value

             end do    ! loop over species

          end do   ! b3
       end do      ! b2
    end do         ! b1

    ! FFT the core density from reciprocal to real space
    call fourier_apply_cell_backward(tcore_density,tcore_density_recip,grid)

    ! ndmh: scale with 1.0/weight
    factor = 1.0_DP / grid%weight
    tcore_density = factor * tcore_density

    ! ndmh: deallocate storage for core density in reciprocal space
    deallocate(tcore_density_recip,stat=ierr)
    call utils_dealloc_check('paw_tcore_density', &
         'tcore_density_recip',ierr)

    call timer_clock('paw_tcore_density',2)

  end subroutine paw_tcore_density


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_species_init_proj(proj_set,elements)

    !=================================================================!
    ! This subroutine allocates and initialises all fftbox_proj_recip !
    ! for each paw_sp element. Each such fftbox_proj_recip is         !
    ! a full non-local projector in reciprocal space in the FFTbox.   !
    !-----------------------------------------------------------------!
    ! Adapted on 21/05/2010.by Nicholas Hine from the routine         !
    ! pseudopot_species_init_proj, originally written by              !
    ! Chris-Kriton Skylaris in 2004 with subsequent modifications     !
    ! by Peter Haynes, Arash Mostofi and Nicholas Hine.               !
    !=================================================================!

    use comms, only: pub_on_root
    use geometry, only: POINT, OPERATOR(.dot.)
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_orig_atom
    use projectors, only: PROJECTOR_SET, projectors_allocate_set, &
         projectors_init_fftbox_recip
    use rundat, only: pub_output_detail
    use simulation_cell, only: pub_cell, pub_fftbox
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    type(PROJECTOR_SET), intent(inout) :: proj_set

    ! Local Variables
    integer :: iat, orig_iat
    integer :: isp
    integer :: shell
    integer :: proj_count

    if (pub_on_root  .and. pub_output_detail == VERBOSE) &
         write(stdout,'(a)') '... PAW projector initialisation'

    ! ndmh: find number of unique projectors
    proj_count = 0
    do isp=1,pub_cell%num_pspecies
       proj_count = proj_count + paw_sp(isp)%npw_tot
    end do

    ! Set up the entries of proj_set
    proj_set%n_proj_species = pub_cell%num_pspecies
    call projectors_allocate_set(proj_set, &
         maxval(paw_sp(:)%npw),maxval(paw_sp(:)%n_recip_pts))

    ! ndmh: set species_num_proj and species_first_proj values
    ! ndmh: also set gmax, n_rad_pts, n_shells, ang_mom and rad_proj_recip
    proj_count = 1
    do isp=1,proj_set%n_proj_species
       proj_set%species_num_proj(isp) = paw_sp(isp)%npw_tot
       proj_set%species_first_proj(isp) = proj_count
       proj_set%gmax(isp) = paw_sp(isp)%g_max
       proj_set%n_rad_pts(isp) = paw_sp(isp)%n_recip_pts
       proj_set%num_shells(isp) = paw_sp(isp)%npw
       proj_set%ang_mom(:,isp) = 0
       proj_set%rad_proj_recip(:,:,isp) = 0.0_DP
       do shell=1,paw_sp(isp)%npw
          proj_set%ang_mom(shell,isp) = paw_sp(isp)%l_pw(shell)
          proj_set%rad_proj_recip(1:paw_sp(isp)%n_recip_pts,shell,isp) = &
               paw_sp(isp)%tproj_recip(1:paw_sp(isp)%n_recip_pts,shell)
       end do
       proj_count = proj_count + proj_set%species_num_proj(isp)
    end do

    ! ndmh: copy projector centre and radius from elements array
    do iat=1,pub_cell%nat
       orig_iat = pub_orig_atom(iat)
       proj_set%proj_centre(iat) = elements(orig_iat)%centre
       proj_set%proj_max_radius(iat) = elements(orig_iat)%max_core_radius
       proj_set%proj_species(iat) = elements(orig_iat)%pspecies_number
    end do

  end subroutine paw_species_init_proj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine paw_species_init_core_wvfns(elements,kpt,swap_rc)

    !=================================================================!
    ! This subroutine allocates and initialises all fftbox_proj_recip !
    ! for each paw_sp element. Each such fftbox_proj_recip is         !
    ! a full non-local projector in reciprocal space in the FFTbox.   !
    !-----------------------------------------------------------------!
    ! Adapted on 21/05/2010.by Nicholas Hine from the routine         !
    ! pseudopot_species_init_proj, originally written by              !
    ! Chris-Kriton Skylaris in 2004 with subsequent modifications     !
    ! by Peter Haynes, Arash Mostofi and Nicholas Hine.               !
    !=================================================================!

    use comms, only: pub_on_root
    use geometry, only: POINT, OPERATOR(.dot.)
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_orig_atom
    use projectors, only: PROJECTOR_SET, projectors_allocate_set, &
         projectors_init_fftbox_recip
    use simulation_cell, only: pub_cell, pub_fftbox
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    type(POINT), intent(in), optional :: kpt
    logical, intent(in), optional :: swap_rc

    ! Local Variables
    integer :: iat, orig_iat
    integer :: isp
    integer :: shell
    integer :: proj_count
    type(POINT) :: kpt_loc
    logical :: loc_swap_rc

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &paw_init_core_wvfns'
#endif

    ! Check for optional arguments
    if (present(kpt)) then
       kpt_loc = kpt
    else
       kpt_loc%x = 0.0_DP ; kpt_loc%y = 0.0_DP ; kpt_loc%z = 0.0_DP
    end if
    if (present(swap_rc)) then
       loc_swap_rc = swap_rc
    else
       loc_swap_rc = .false.
    end if

    ! ndmh: find number of unique projectors
    proj_count = 0
    do isp=1,pub_cell%num_pspecies
       proj_count = proj_count + paw_sp(isp)%n_core_wfs_tot
    end do

    ! Set up the entries of paw_projectors
    paw_core_wvfns%n_proj_species = pub_cell%num_pspecies
    call projectors_allocate_set(paw_core_wvfns, &
         maxval(paw_sp(:)%n_core_wfs),maxval(paw_sp(:)%n_recip_pts))

    ! ndmh: set species_num_proj and species_first_proj values
    proj_count = 1
    do isp=1,paw_core_wvfns%n_proj_species
       paw_core_wvfns%species_num_proj(isp) = paw_sp(isp)%n_core_wfs_tot
       paw_core_wvfns%gmax(isp) = paw_sp(isp)%g_max
       paw_core_wvfns%n_rad_pts(isp) = paw_sp(isp)%n_recip_pts
       paw_core_wvfns%num_shells(isp) = paw_sp(isp)%n_core_wfs
       paw_core_wvfns%ang_mom(:,isp) = 0
       paw_core_wvfns%rad_proj_recip(:,:,isp) = 0.0_DP
       do shell=1,paw_sp(isp)%n_core_wfs
          paw_core_wvfns%ang_mom(shell,isp) = paw_sp(isp)%l_core_wf(shell)
          paw_core_wvfns%rad_proj_recip(1:paw_sp(isp)%n_recip_pts,shell,isp) = &
               paw_sp(isp)%core_wf_recip(1:paw_sp(isp)%n_recip_pts,shell)
       end do
       paw_core_wvfns%species_first_proj(isp) = proj_count
       proj_count = proj_count + paw_core_wvfns%species_num_proj(isp)
    end do

    ! ndmh: copy projector centre and radius from elements array
    do iat=1,pub_cell%nat
       orig_iat = pub_orig_atom(iat)
       paw_core_wvfns%proj_centre(iat) = elements(orig_iat)%centre
       paw_core_wvfns%proj_max_radius(iat) = elements(orig_iat)%max_core_wf_radius
       paw_core_wvfns%proj_species(iat) = elements(orig_iat)%pspecies_number
    end do

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &paw_init_core_wvfns'
#endif

  end subroutine paw_species_init_core_wvfns


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_projector_overlap(proj_overlap)

    !==================================================================!
    ! This subroutine returns the PAW-projector overlap matrix in      !
    ! SPAM3 format.                                                    !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  proj_overlap (inout) : The SPAM3 block-diagonal overlap matrix  !
    !  of the partial waves of each atom                               !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 06/06/10.           !
    !==================================================================!

    use comms, only: pub_my_node_id
    use gaunt_coeff, only: realgaunt
    use parallel_strategy, only: pub_elements_on_node, &
         pub_num_atoms_on_node, pub_first_atom_on_node
    use sparse, only: SPAM3, sparse_put_block
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: proj_overlap

    ! Local Variables
    integer :: iat
    integer :: loc_iat
    integer :: isp
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: li,mi,lj,mj
    integer :: ierr
    real(kind=DP),allocatable :: ovlp_block(:,:)

    allocate(ovlp_block(max_paw_proj_tot,max_paw_proj_tot),stat=ierr)
    call utils_alloc_check('paw_projector_overlap','ovlp_block',ierr)

    ! Loop over atoms on this node
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1
       isp = pub_elements_on_node(loc_iat)%pspecies_number

       ! Double loop over partial waves i,j
       do ipwtot=1,paw_sp(isp)%npw_tot
          ipw = paw_sp(isp)%ipw_tot(ipwtot)
          li = paw_sp(isp)%l_pw_tot(ipwtot)
          mi = paw_sp(isp)%m_pw_tot(ipwtot)
          do jpwtot=1,paw_sp(isp)%npw_tot
             jpw = paw_sp(isp)%ipw_tot(jpwtot)
             lj = paw_sp(isp)%l_pw_tot(jpwtot)
             mj = paw_sp(isp)%m_pw_tot(jpwtot)
             ! Find O_ij = sqrt(4pi)*n_nilinjlj^0*G_limiljmj^00
             ovlp_block(ipwtot,jpwtot) = &
                  paw_sp(isp)%aug_nLij(ipw,jpw,0) * &
                  realgaunt(0,0,li,mi,lj,mj) * sqrt_4pi
          end do
       end do

       ! Put this atom's block into SPAM3
       call sparse_put_block(ovlp_block,proj_overlap,iat,iat)
    end do

    deallocate(ovlp_block,stat=ierr)
    call utils_dealloc_check('paw_projector_overlap','ovlp_block',ierr)

  end subroutine paw_projector_overlap


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_position_operator(pos_ij,axis)

    !==================================================================!
    ! This subroutine returns the PAW sphere part of the position      !
    ! operator, defined as <phi_i|r|phi_j> - <tphi_i|r|tphi_j> for     !
    ! each pair of partial waves i,j                                   !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  pos_ij(3) (inout) : The SPAM3 block-diagonal position matrix    !
    !  between the partial waves of each atom                          !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 06/10/11.           !
    !==================================================================!

    use comms, only: pub_my_node_id
    use gaunt_coeff, only: realgaunt
    use parallel_strategy, only: pub_elements_on_node, &
         pub_num_atoms_on_node, pub_first_atom_on_node
    use sparse, only: SPAM3, sparse_put_block
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: pos_ij(3)
    integer, intent(in), optional :: axis

    ! Local Variables
    integer :: axmin,axmax
    integer :: iat
    integer :: loc_iat
    integer :: isp
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: li,mi,lj,mj,lr,mr
    integer :: icart
    integer :: ierr
    real(kind=DP),allocatable :: pos_block(:,:,:)

    ! Optional variable
    if (present(axis)) then
       axmin = axis
       axmax = axis
    else
       axmin = 1
       axmax = 3
    end if

    allocate(pos_block(max_paw_proj_tot,max_paw_proj_tot,3),stat=ierr)
    call utils_alloc_check('paw_position_operator','pos_block',ierr)

    ! Expansion of r in spherical harmonics has only L=1 components
    lr = 1

    ! Loop over atoms on this node
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1
       isp = pub_elements_on_node(loc_iat)%pspecies_number

       ! Double loop over partial waves i,j
       do ipwtot=1,paw_sp(isp)%npw_tot
          ipw = paw_sp(isp)%ipw_tot(ipwtot)
          li = paw_sp(isp)%l_pw_tot(ipwtot)
          mi = paw_sp(isp)%m_pw_tot(ipwtot)
          do jpwtot=1,paw_sp(isp)%npw_tot
             jpw = paw_sp(isp)%ipw_tot(jpwtot)
             lj = paw_sp(isp)%l_pw_tot(jpwtot)
             mj = paw_sp(isp)%m_pw_tot(jpwtot)
             ! Find r_ij = sqrt(4pi)*n_nilinjlj^0*G_limiljmj^1M
             do icart=axmin,axmax
                if (icart==1) mr = 1
                if (icart==2) mr = -1
                if (icart==3) mr = 0
                pos_block(ipwtot,jpwtot,icart) = &
                     paw_sp(isp)%aug_nLij(ipw,jpw,lr) * &
                     realgaunt(lr,mr,li,mi,lj,mj) * sqrt_4pi / sqrt(3.0_DP)
             end do
          end do
       end do

       ! Put this atom's blocks into SPAM3 matrices
       do icart=axmin,axmax
          call sparse_put_block(pos_block(:,:,icart),pos_ij(icart),iat,iat)
       end do
    end do

    deallocate(pos_block,stat=ierr)
    call utils_dealloc_check('paw_position_operator','pos_block',ierr)

  end subroutine paw_position_operator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_grad_operator(grad_ij,axis)

    !==================================================================!
    ! This subroutine returns the PAW sphere part of the grad          !
    ! operator, defined as <phi_i|nabla|phi_j> - <tphi_i|nabla|tphi_j> !
    ! for each pair of partial waves i,j                               !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  grad_ij(3) (inout) : The SPAM3 block-diagonal grad operator     !
    !  matrix between the partial waves of each atom                   !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 02/12/11.           !
    !==================================================================!

    use comms, only: pub_my_node_id
    use gaunt_coeff, only: realgaunt, gaunt_grad_integrals
    use parallel_strategy, only: pub_elements_on_node, &
         pub_num_atoms_on_node, pub_first_atom_on_node
    use sparse, only: SPAM3, sparse_put_block
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: grad_ij(3)
    integer, intent(in), optional :: axis

    ! Local Variables
    integer :: axmin,axmax
    integer :: iat
    integer :: loc_iat
    integer :: isp
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: li,mi,lj,mj
    integer :: icart
    integer :: ierr
    real(kind=DP),allocatable :: grad_block(:,:,:)
    real(kind=DP) :: angint(3,3)

    ! Optional variable
    if (present(axis)) then
       axmin = axis
       axmax = axis
    else
       axmin = 1
       axmax = 3
    end if

    allocate(grad_block(max_paw_proj_tot,max_paw_proj_tot,3),stat=ierr)
    call utils_alloc_check('paw_grad_operator','grad_block',ierr)

    ! Loop over atoms on this node
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1
       isp = pub_elements_on_node(loc_iat)%pspecies_number

       ! Double loop over partial waves i,j
       do ipwtot=1,paw_sp(isp)%npw_tot
          ipw = paw_sp(isp)%ipw_tot(ipwtot)
          li = paw_sp(isp)%l_pw_tot(ipwtot)
          mi = paw_sp(isp)%m_pw_tot(ipwtot)
          do jpwtot=1,paw_sp(isp)%npw_tot
             jpw = paw_sp(isp)%ipw_tot(jpwtot)
             lj = paw_sp(isp)%l_pw_tot(jpwtot)
             mj = paw_sp(isp)%m_pw_tot(jpwtot)
             ! Find grad_ij
             call gaunt_grad_integrals(angint,li,mi,lj,mj)
             do icart=axmin,axmax
                grad_block(ipwtot,jpwtot,icart) = &
                     paw_sp(isp)%aug_nLij(ipw,jpw,-2) * angint(icart,1) &
                     + paw_sp(isp)%aug_nLij(ipw,jpw,-1) * (-angint(icart,1) &
                     + angint(icart,2) + angint(icart,3))
             end do
          end do
       end do

       ! Put this atom's blocks into SPAM3 matrices
       do icart=axmin,axmax
          call sparse_put_block(grad_block(:,:,icart),grad_ij(icart),iat,iat)
       end do
    end do

    deallocate(grad_block,stat=ierr)
    call utils_dealloc_check('paw_grad_operator','grad_block',ierr)

  end subroutine paw_grad_operator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_core_position_operator(pos_ij)

    !==================================================================!
    ! This subroutine returns the PAW sphere part of the position      !
    ! operator between partial waves and core wavefunctions, defined   !
    ! as <phi_c|r|phi_j> - <tphi_i|r|tphi_j> for each combination      !
    ! of partial wave i and core wavefunction c                        !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  pos_ij(3) (inout) : The SPAM3 block-diagonal position matrix    !
    !  between the partial waves and core wavefunctions of each atom   !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 06/10/11.           !
    ! Adapted for Core Wavefunctions by Laura Ratcliff.                !
    !==================================================================!

    use comms, only: pub_my_node_id
    use gaunt_coeff, only: realgaunt
    use parallel_strategy, only: pub_elements_on_node, &
         pub_num_atoms_on_node, pub_first_atom_on_node
    use sparse, only: SPAM3, sparse_put_block
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: pos_ij(3)

    ! Local Variables
    integer :: iat
    integer :: loc_iat
    integer :: isp
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: li,mi,lj,mj,lr,mr
    integer :: icart
    integer :: ierr
    real(kind=DP),allocatable :: pos_block(:,:,:)

    ! lr408: might need to change this - presumably max_paw_proj_tot will always be smaller
    ! than number of core wfs but really need to add some kind of check
    allocate(pos_block(max_core_wf_tot,max_paw_proj_tot,3),stat=ierr)
    call utils_alloc_check('paw_core_position_operator','pos_block',ierr)

    ! Expansion of r in spherical harmonics has only L=1 components
    lr = 1

    ! Loop over atoms on this node
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1
       isp = pub_elements_on_node(loc_iat)%pspecies_number

       ! but core wvfns don't exist on each atom so need to check that beforehand!
       if (.not. paw_sp(isp)%core_wvfns_exist) cycle

       ! Double loop over partial waves i,j
       do ipwtot=1,paw_sp(isp)%n_core_wfs_tot!npw_tot
          ipw = paw_sp(isp)%icore_wf_tot(ipwtot)
          li = paw_sp(isp)%l_core_wf_tot(ipwtot)
          mi = paw_sp(isp)%m_core_wf_tot(ipwtot)
          do jpwtot=1,paw_sp(isp)%npw_tot
             jpw = paw_sp(isp)%ipw_tot(jpwtot)
             lj = paw_sp(isp)%l_pw_tot(jpwtot)
             mj = paw_sp(isp)%m_pw_tot(jpwtot)
             ! Find r_ij = sqrt(4pi)*n_nilinjlj^0*G_limiljmj^1M
             do icart=1,3
                if (icart==1) mr = 1
                if (icart==2) mr = -1
                if (icart==3) mr = 0
                pos_block(ipwtot,jpwtot,icart) = &
                     paw_sp(isp)%core_aug_nLij(ipw,jpw,lr) * &
                     realgaunt(lr,mr,li,mi,lj,mj) * sqrt_4pi
             end do
          end do
       end do

       ! Put this atom's blocks into SPAM3 matrices
       do icart=1,3
          call sparse_put_block(pos_block(:,:,icart),pos_ij(icart),iat,iat)
       end do
    end do

    deallocate(pos_block,stat=ierr)
    call utils_dealloc_check('paw_core_position_operator','pos_block',ierr)

  end subroutine paw_core_position_operator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_projector_denskern_init(rhoij0)

    !=================================================================!
    ! This subroutine calculates an initial guess for the projector   !
    ! density kernel, based on the atomic rhoij0 values found in the  !
    ! PAW dataset. This is used to calculate the nonlocal energies    !
    ! for the Hamiltonian used in the Palser-Manolopoulos canonical   !
    ! purification.                                                   !
    !-----------------------------------------------------------------!
    ! Arguments:                                                      !
    !  rhoij0 (out) : The initial projector density kernel.           !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 04/06/10.                           !
    !=================================================================!

    use comms, only: pub_my_node_id
    use parallel_strategy, only: pub_num_atoms_on_node, &
         pub_first_atom_on_node, pub_elements_on_node
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_put_block
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: rhoij0(pub_cell%num_spins)

    ! Local Variables
    integer :: loc_iat, iat
    integer :: isp
    integer :: is
    integer :: ipwtot,jpwtot
    integer :: ierr
    real(kind=DP), allocatable :: rhoij0_sp(:,:,:)

    allocate(rhoij0_sp(max_paw_proj_tot,max_paw_proj_tot, &
         pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('paw_projector_denskern_init','rhoij0_sp',ierr)

    ! Loop over species, setting up array of rhoij0 terms
    do isp=1,pub_cell%num_pspecies
       if (any(pub_elements_on_node(:)%pspecies_number==isp)) then
          do ipwtot=1,paw_sp(isp)%npw_tot
             do jpwtot=1,paw_sp(isp)%npw_tot
                rhoij0_sp(ipwtot,jpwtot,isp) = &
                     paw_sp(isp)%rhoij0(ipwtot,jpwtot) &
                     / real(pub_cell%num_spins,kind=DP)
             end do
          end do
       end if
    end do

    ! Loop over atoms
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = loc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
       isp = pub_elements_on_node(loc_iat)%pspecies_number
       ! Put block of rhoij0 into SPAM3 matrix
       do is=1,pub_cell%num_spins
          call sparse_put_block(rhoij0_sp(:,:,isp),rhoij0(is),iat,iat)
       end do
    end do

    deallocate(rhoij0_sp,stat=ierr)
    call utils_dealloc_check('paw_projector_denskern_init','rhoij0_sp',ierr)

  end subroutine paw_projector_denskern_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_atom_aug_den(atom_nhat,atom_grad_nhat,total_nhat, &
       total_nhat_targ,rho_ij_block, &
       isp,atom_centre,grid,box_n1,box_n2,box_n3, &
       box_start1,box_start2,box_start3,atom_aug_func_real, &
       atom_aug_func_grad_real,atom_aug_func_recip,atom_aug_func_grad_recip, &
       atom_nhat_recip,atom_grad_nhat_recip)

    use cell_grid, only: GRID_INFO
    use constants, only: max_spins
    use fourier, only: fourier_apply_box
    use gaunt_coeff, only: realgaunt
    use geometry, only: POINT
    use paw_shape, only: paw_shape_gLSLM_recip, paw_shape_gLSLM_real, &
         paw_shape_grad_gLSLM_real, paw_shape_grad_gLSLM_recip
    use rundat, only: pub_aug_funcs_recip
    use simulation_cell, only: pub_cell
    use xc, only: pub_xc_gradient_corrected

    ! Arguments
    integer, intent(in) :: box_n1,box_n2,box_n3
    real(kind=DP), intent(out) :: atom_nhat(box_n1,box_n2,box_n3, &
         pub_cell%num_spins)
    real(kind=DP), intent(out) :: atom_grad_nhat(box_n1,box_n2,box_n3, &
         pub_cell%num_spins,3)
    real(kind=DP), intent(out) :: total_nhat(max_spins)
    real(kind=DP), intent(out) :: total_nhat_targ(max_spins)
    type(GRID_INFO), intent(in) :: grid
    type(POINT), intent(in) :: atom_centre
    real(kind=DP), intent(in) :: rho_ij_block(:,:,:)
    integer, intent(in) :: isp
    integer, intent(in) :: box_start1,box_start2,box_start3
    real(kind=DP), intent(out), optional :: atom_aug_func_real(box_n1,box_n2, &
         box_n3)
    real(kind=DP), intent(out), optional :: atom_aug_func_grad_real(box_n1, &
         box_n2,box_n3,3)
    complex(kind=DP), intent(out), optional :: atom_aug_func_recip(box_n1, &
         box_n2,box_n3)
    complex(kind=DP), intent(out), optional :: atom_aug_func_grad_recip(box_n1, &
         box_n2,box_n3,3)
    complex(kind=DP), intent(out), optional :: atom_nhat_recip(box_n1, &
         box_n2,box_n3,pub_cell%num_spins)
    complex(kind=DP), intent(out), optional :: atom_grad_nhat_recip(box_n1, &
         box_n2,box_n3,pub_cell%num_spins,3)


    ! Local Variables
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: lup,mup
    integer :: li,mi,lj,mj
    integer :: is
    integer :: cart
    real(kind=DP) :: rgij
    real(kind=DP) :: qijLM
    real(kind=DP) :: qLM(max_spins)
    real(kind=DP) :: total_gLSLM

    if (pub_aug_funcs_recip) then
       atom_nhat_recip = 0.0_DP
       if (pub_xc_gradient_corrected) atom_grad_nhat_recip = 0.0_DP
    else
       atom_nhat = 0.0_DP
       if (pub_xc_gradient_corrected) atom_grad_nhat = 0.0_DP
    end if

    do lup=0,paw_sp(isp)%lmax
       do mup=-lup,lup

          qLM(:) = 0.0_DP

          ! Loop over partial waves to find total target comp. charge
          ! and calculate q_LM^\sigma for this L,M,sigma channel
          do ipwtot=1,paw_sp(isp)%npw_tot
             ipw = paw_sp(isp)%ipw_tot(ipwtot)
             li = paw_sp(isp)%l_pw_tot(ipwtot)
             mi = paw_sp(isp)%m_pw_tot(ipwtot)
             do jpwtot=1,paw_sp(isp)%npw_tot
                jpw = paw_sp(isp)%ipw_tot(jpwtot)
                lj = paw_sp(isp)%l_pw_tot(jpwtot)
                mj = paw_sp(isp)%m_pw_tot(jpwtot)
                rgij = realgaunt(lup,mup,lj,mj,li,mi)
                qijLM = paw_sp(isp)%aug_nLij(ipw,jpw,lup)*rgij
                do is=1,pub_cell%num_spins
                   qLM(is) = qLM(is) + &
                        rho_ij_block(ipwtot,jpwtot,is) * qijLM
                   if (lup==0) &
                        total_nhat_targ(is) = total_nhat_targ(is) + &
                        qijLM*rho_ij_block(ipwtot,jpwtot,is) * sqrt_4pi
                end do

             end do
          end do

          if (any(abs(qLM(:))>1e-20_DP)) then

             ! Get shape function g_L(r) multiplied by spherical harmonic
             ! S_LM(r) for this L,M pair
             if (pub_aug_funcs_recip) then
                ! Get augmentation function in reciprocal space
                call paw_shape_gLSLM_recip(atom_aug_func_recip, &
                     paw_sp(isp)%shape, &
                     lup,mup,box_start1,box_start2,box_start3,grid, &
                     box_n1,box_n2,box_n3,atom_centre)
                total_gLSLM = real(atom_aug_func_recip(1,1,1),kind=DP)

                do is=1,pub_cell%num_spins
                   atom_nhat_recip(:,:,:,is) = atom_nhat_recip(:,:,:,is) &
                        + qLM(is)*atom_aug_func_recip(:,:,:)
                   total_nhat(is) = total_nhat(is) + total_gLSLM*qLM(is)
                end do

                ! Find gradient of \hat{n}(r) if required
                if (pub_xc_gradient_corrected) then
                   ! Get gradient of (shape function g_L(r) multiplied by
                   ! spherical harmonic S_LM(r) for this L,M pair)
                   call paw_shape_grad_gLSLM_recip(atom_aug_func_grad_recip, &
                        paw_sp(isp)%shape,lup,mup,box_start1,box_start2, &
                        box_start3,grid,box_n1,box_n2,box_n3,atom_centre)
                   do cart=1,3
                      do is=1,pub_cell%num_spins
                         atom_grad_nhat_recip(:,:,:,is,cart) = &
                              atom_grad_nhat_recip(:,:,:,is,cart) + &
                              qLM(is)*atom_aug_func_grad_recip(:,:,:,cart)
                      end do
                   end do
                end if


             else

                ! Get augmentation function in real space
                call paw_shape_gLSLM_real(atom_aug_func_real,paw_sp(isp)%shape, &
                     lup,mup,box_start1,box_start2,box_start3,grid, &
                     box_n1,box_n2,box_n3,atom_centre)
                total_gLSLM = sum(atom_aug_func_real(:,:,:))*grid%weight

                do is=1,pub_cell%num_spins
                   atom_nhat(:,:,:,is) = atom_nhat(:,:,:,is) &
                        + qLM(is)*atom_aug_func_real(:,:,:)
                   total_nhat(is) = total_nhat(is) + total_gLSLM*qLM(is)
                end do

                ! Find gradient of \hat{n}(r) if required
                if (pub_xc_gradient_corrected) then
                   ! Get gradient of (shape function g_L(r) multiplied by
                   ! spherical harmonic S_LM(r) for this L,M pair)
                   call paw_shape_grad_gLSLM_real(atom_aug_func_grad_real, &
                        paw_sp(isp)%shape,lup,mup,box_start1,box_start2, &
                        box_start3,grid,box_n1,box_n2,box_n3,atom_centre)
                   do cart=1,3
                      do is=1,pub_cell%num_spins
                         atom_grad_nhat(:,:,:,is,cart) = &
                              atom_grad_nhat(:,:,:,is,cart) + &
                              qLM(is)*atom_aug_func_grad_real(:,:,:,cart)
                      end do
                   end do
                end if

             end if
          end if
       end do
    end do

    if (pub_aug_funcs_recip) then

       ! Fourier transform to real space in augmentation box
       do is=1,pub_cell%num_spins
          call fourier_apply_box('F','B',atom_nhat_recip(:,:,:,is),aug=.true.)
          atom_nhat(:,:,:,is) = real(atom_nhat_recip(:,:,:,is),kind=DP) &
               / grid%weight
       end do

       ! Fourier transform gradient of \hat{n}(r) if required
       if (pub_xc_gradient_corrected) then
          do cart=1,3
             do is=1,pub_cell%num_spins
                call fourier_apply_box('F','B', &
                     atom_grad_nhat_recip(:,:,:,is,cart),aug=.true.)
                atom_grad_nhat(:,:,:,is,cart) = &
                     real(atom_grad_nhat_recip(:,:,:,is,cart),kind=DP)/grid%weight
             end do
          end do
       end if
    end if

  end subroutine paw_atom_aug_den


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_atom_aug_integrals(dijhat_at, &
       num_spins,isp,atom_centre,grid,box_n1,box_n2,box_n3, &
       box_start1,box_start2,box_start3,locpot_box_real,atom_aug_func_real, &
       locpot_box_recip,atom_aug_func_recip)

    use cell_grid, only: GRID_INFO
    use constants, only: max_spins
    use gaunt_coeff, only: realgaunt
    use geometry, only: POINT
    use paw_shape, only: paw_shape_gLSLM_recip, paw_shape_gLSLM_real
    use rundat, only: pub_aug_funcs_recip
    use utils, only: utils_abort

    ! Arguments
    integer, intent(in) :: box_n1,box_n2,box_n3
    integer, intent(in) :: num_spins
    type(GRID_INFO), intent(in) :: grid
    type(POINT), intent(in) :: atom_centre
    real(kind=DP), intent(inout) :: dijhat_at(:,:,:)
    integer, intent(in) :: isp
    integer, intent(in) :: box_start1,box_start2,box_start3
    real(kind=DP), intent(in), optional :: locpot_box_real(box_n1,box_n2, &
         box_n3,num_spins)
    real(kind=DP), intent(out), optional :: atom_aug_func_real(box_n1,box_n2, &
         box_n3)
    complex(kind=DP), intent(in), optional :: locpot_box_recip(box_n1,box_n2, &
         box_n3,num_spins)
    complex(kind=DP), intent(out), optional :: atom_aug_func_recip(box_n1, &
         box_n2,box_n3)

    ! Local Variables
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: lup,mup
    integer :: li,mi,lj,mj
    integer :: is
    real(kind=DP) :: rgij
    real(kind=DP) :: qijLM
    real(kind=DP) :: locpot_gLSLM_product(max_spins)

    ! Check arguments
    if (pub_aug_funcs_recip) then
       if ((.not.present(locpot_box_recip)) .or. &
            (.not.present(atom_aug_func_recip))) then
          call utils_abort('Error in paw_atom_aug_integrals: Incorrect &
               &arguments supplied')
       end if
    else
       if ((.not.present(locpot_box_real)) .or. &
            (.not.present(atom_aug_func_real))) then
          call utils_abort('Error in paw_atom_aug_integrals: Incorrect &
               &arguments supplied')
       end if
    end if

    ! Loop over angular momentum channels L,M
    do lup=0,paw_sp(isp)%lmax
       do mup=-lup,lup

          ! Get shape function g_L(r) multiplied by spherical harmonic
          ! S_LM(r) for this L,M pair and integrate with locpot

          if (pub_aug_funcs_recip) then

             ! Get augmentation function in reciprocal space
             call paw_shape_gLSLM_recip(atom_aug_func_recip,paw_sp(isp)%shape, &
                  lup,mup,box_start1,box_start2,box_start3,grid, &
                  box_n1,box_n2,box_n3,atom_centre)

             ! Calculate integral \int veff(r) g_L(r) S_LM(r) dr
             ! as sum in reciprocal space: \sum_G veff(G) (g_L S_LM)(G)
             do is=1,num_spins
                locpot_gLSLM_product(is) = sum(locpot_box_recip(:,:,:,is) &
                     * atom_aug_func_recip(:,:,:))
             end do

          else

             ! Get augmentation function in real space
             call paw_shape_gLSLM_real(atom_aug_func_real,paw_sp(isp)%shape, &
                  lup,mup,box_start1,box_start2,box_start3,grid, &
                  box_n1,box_n2,box_n3,atom_centre)

             ! Calculate integral \int veff(r) g_L(r) S_LM(r) dr
             do is=1,num_spins
                locpot_gLSLM_product(is) = sum(locpot_box_real(:,:,:,is) &
                     * atom_aug_func_real(:,:,:)) * grid%weight
             end do

          end if

          ! Loop over spins and double loop over partial waves
          do is=1,num_spins
             do ipwtot=1,paw_sp(isp)%npw_tot
                ipw = paw_sp(isp)%ipw_tot(ipwtot)
                li = paw_sp(isp)%l_pw_tot(ipwtot)
                mi = paw_sp(isp)%m_pw_tot(ipwtot)
                do jpwtot=1,paw_sp(isp)%npw_tot
                   jpw = paw_sp(isp)%ipw_tot(jpwtot)
                   lj = paw_sp(isp)%l_pw_tot(jpwtot)
                   mj = paw_sp(isp)%m_pw_tot(jpwtot)
                   rgij = realgaunt(lup,mup,li,mi,lj,mj)
                   qijLM = paw_sp(isp)%aug_nLij(ipw,jpw,lup)*rgij

                   dijhat_at(ipwtot,jpwtot,is) = &
                        dijhat_at(ipwtot,jpwtot,is) + &
                        locpot_gLSLM_product(is) * qijLM

                end do  ! jpwtot
             end do  ! ipwtot
          end do  ! is

       end do  ! mup
    end do  ! lup

  end subroutine paw_atom_aug_integrals


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_atom_aug_force(nhat_force, &
       rhoij_at,num_spins,isp,atom_centre,grid,box_n1,box_n2,box_n3, &
       box_start1,box_start2,box_start3,locpot_box_real, &
       atom_grad_aug_func_real,locpot_box_recip,atom_grad_aug_func_recip)

    use cell_grid, only: GRID_INFO
    use constants, only: max_spins
    use gaunt_coeff, only: realgaunt
    use geometry, only: POINT
    use paw_shape, only: paw_shape_grad_gLSLM_real, paw_shape_grad_gLSLM_recip
    use rundat, only: pub_aug_funcs_recip
    use simulation_cell, only: pub_cell

    ! Arguments
    real(kind=DP),intent(inout) :: nhat_force(3)
    integer,intent(in) :: box_n1,box_n2,box_n3
    integer,intent(in) :: num_spins
    type(GRID_INFO),intent(in) :: grid
    type(POINT),intent(in) :: atom_centre
    real(kind=DP),intent(in) :: rhoij_at(:,:,:)
    integer,intent(in) :: isp
    integer,intent(in) :: box_start1,box_start2,box_start3
    real(kind=DP),intent(out),optional :: atom_grad_aug_func_real(box_n1, &
         box_n2,box_n3,3)
    real(kind=DP),intent(in),optional :: locpot_box_real(box_n1,box_n2, &
         box_n3,num_spins)
    complex(kind=DP),intent(out),optional :: atom_grad_aug_func_recip(box_n1, &
         box_n2,box_n3,3)
    complex(kind=DP),intent(in),optional :: locpot_box_recip(box_n1,box_n2, &
         box_n3,num_spins)

    ! Local Variables
    integer :: ipw,jpw
    integer :: ipwtot,jpwtot
    integer :: lup,mup
    integer :: li,mi,lj,mj
    integer :: is
    integer :: cart
    real(kind=DP) :: rgij
    real(kind=DP) :: qijLM
    real(kind=DP) :: locpot_grad_gLSLM_product(3,max_spins)

    ! Loop over angular momentum channels L,M
    do lup=0,paw_sp(isp)%lmax
       do mup=-lup,lup

          ! Construct gradient d/dR_I (g_L(r)*S_LM(\hat{r})) of shape
          ! function g_L(r) multiplied by spherical harmonic S_LM(r) for
          ! this L,M pair, then integrate it with part of locpot in aug box

          if (pub_aug_funcs_recip) then

             call paw_shape_grad_gLSLM_recip(atom_grad_aug_func_recip, &
                  paw_sp(isp)%shape,lup,mup,box_start1,box_start2, &
                  box_start3,grid,box_n1,box_n2,box_n3, &
                  atom_centre)

             ! Loop over spins and loop over cartesian directions
             do is=1,num_spins
                do cart=1,3
                   ! Calculate integral \int veff(r) d/dR(g_L(r) S_LM(r)) dr
                   locpot_grad_gLSLM_product(cart,is) = &
                        sum(locpot_box_recip(:,:,:,is) &
                        * atom_grad_aug_func_recip(:,:,:,cart))
                end do
             end do

          else

             call paw_shape_grad_gLSLM_real(atom_grad_aug_func_real, &
                  paw_sp(isp)%shape,lup,mup,box_start1,box_start2, &
                  box_start3,grid,box_n1,box_n2,box_n3, &
                  atom_centre)

             ! Loop over spins and loop over cartesian directions
             do is=1,num_spins
                do cart=1,3
                   ! Calculate integral \int veff(r) d/dR(g_L(r) S_LM(r)) dr
                   locpot_grad_gLSLM_product(cart,is) = &
                        sum(locpot_box_real(:,:,:,is) &
                        * atom_grad_aug_func_real(:,:,:,cart)) * grid%weight
                end do
             end do

          end if

          ! Double loop over partial waves
          do ipwtot=1,paw_sp(isp)%npw_tot
             ipw = paw_sp(isp)%ipw_tot(ipwtot)
             li = paw_sp(isp)%l_pw_tot(ipwtot)
             mi = paw_sp(isp)%m_pw_tot(ipwtot)
             do jpwtot=1,paw_sp(isp)%npw_tot
                jpw = paw_sp(isp)%ipw_tot(jpwtot)
                lj = paw_sp(isp)%l_pw_tot(jpwtot)
                mj = paw_sp(isp)%m_pw_tot(jpwtot)
                rgij = realgaunt(lup,mup,li,mi,lj,mj)
                qijLM = paw_sp(isp)%aug_nLij(ipw,jpw,lup)*rgij

                do is=1,num_spins
                   nhat_force(:) = nhat_force(:) &
                        + (locpot_grad_gLSLM_product(:,is) * qijLM &
                        * rhoij_at(ipwtot,jpwtot,is))
                end do  ! is

             end do  ! jpwtot
          end do  ! ipwtot

       end do  ! mup
    end do  ! lup

  end subroutine paw_atom_aug_force


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_nonlocal_energies(dij,rho_ij,paw_sphere_energies,show_matrices)

    !==================================================================!
    ! This subroutine creates the PAW nonlocal energies Dij, given by  !
    !  D_ij = \hat{D}_ij + D^1_ij - \tilde{D}^1_ij                     !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  dij (out) : nonlocal matrix in Nproj x Nproj sized matrix       !
    !  rhoij (in) : projector-reduced density kernel (Nproj-sized)     !
    !  grid (in) : Grid definition for fine grid                       !
    !  paw_sphere_energies (inout) : array of energy contributions     !
    !    from sphere terms                                             !
    !  show_matrices (in,opt) : Whether it is a good time to display   !
    !    the intermediate matrices (eg during energy components report)!
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.           !
    !==================================================================!

    use comms, only: pub_my_node_id, pub_on_root
    use constants, only: paw_en_size, paw_en_dij0, &
         paw_en_ehart, paw_en_exc, paw_en_exc_dc, paw_en_etxc, paw_en_etxc_dc, &
         paw_en_dijxc, paw_en_exc_core
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: pub_num_atoms_on_node
    use rundat, only: pub_paw_output_detail, xc_functional
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, &
         sparse_trace, sparse_axpy, sparse_show_matrix
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: dij(pub_cell%num_spins)
    type(SPAM3), intent(in) :: rho_ij(pub_cell%num_spins)
    real(kind=DP), intent(inout), optional :: paw_sphere_energies(paw_en_size)
    logical, intent(in), optional :: show_matrices

    ! Local Variables
    type(SPAM3), allocatable :: dij0(:)
    type(SPAM3), allocatable :: dij_hartree(:)
    type(SPAM3), allocatable :: dij_xc(:)
    real(kind=DP) :: exc,etxc,exc_dc,etxc_dc,edijxc,exc_core
    integer :: ierr
    integer :: is
    logical :: loc_show_matrices

    ! Handle optional arguments
    if (present(paw_sphere_energies)) paw_sphere_energies(:) = 0.0_DP
    loc_show_matrices = .false.
    if (present(show_matrices)) loc_show_matrices = show_matrices

    ! Allocate temporary matrices
    allocate(dij0(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('paw_nonlocal_energies','dij0',ierr)
    allocate(dij_hartree(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('paw_nonlocal_energies','dij_hartree',ierr)
    allocate(dij_xc(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('paw_nonlocal_energies','dij_xc',ierr)

    ! Allocate temporary matrices
    do is=1,pub_cell%num_spins
       call sparse_create(dij0(is),rho_ij(is))
       call sparse_create(dij_hartree(is),rho_ij(is))
       call sparse_create(dij_xc(is),rho_ij(is))
    end do

    if ((loc_show_matrices).and.(pub_paw_output_detail>=VERBOSE)) then
       do is=1,pub_cell%num_spins
          if (pub_on_root) write(stdout,'(a,i4)') 'rho_ij', is
          call sparse_show_matrix(rho_ij(is),show_elems=.true.)
       end do
    end if

    ! Calculate the fixed kinetic and Hartree contributions to the nonlocal
    ! energies
    call paw_dij0(dij0)
    do is=1,pub_cell%num_spins
       if ((loc_show_matrices).and.(pub_paw_output_detail>=VERBOSE)) then
          if (pub_on_root) write(stdout,'(a,i4)') 'dij0', is
          call sparse_show_matrix(dij0(is),show_elems=.true.)
       end if
       call sparse_axpy(dij(is),dij0(is),1.0_DP)
       if (present(paw_sphere_energies)) then
          paw_sphere_energies(paw_en_dij0) = paw_sphere_energies(paw_en_dij0) &
               + sparse_trace(rho_ij(is),dij0(is))
       end if
    end do

    ! Calculate the remaining Hartree contribution to the nonlocal energies
    call paw_dij_hartree(dij_hartree,rho_ij)
    do is=1,pub_cell%num_spins
       if ((loc_show_matrices).and.(pub_paw_output_detail>=VERBOSE)) then
          if (pub_on_root) write(stdout,'(a,i4)') 'dij_hartree', is
          call sparse_show_matrix(dij_hartree(is),show_elems=.true.)
       end if
       call sparse_axpy(dij(is),dij_hartree(is),1.0_DP)
       if (present(paw_sphere_energies)) then
          paw_sphere_energies(paw_en_ehart) = &
               paw_sphere_energies(paw_en_ehart) &
               + 0.5_DP*sparse_trace(rho_ij(is),dij_hartree(is))
       end if
    end do

    ! Calculate the exchange-correlation part of the nonlocal energies
    if (xc_functional/='NONE') then
       call paw_exc_core_tot(exc_core)
       call paw_dij_xc(dij_xc,rho_ij,exc,exc_dc,etxc,etxc_dc)
       edijxc = 0.0_DP
       do is=1,pub_cell%num_spins
          if ((loc_show_matrices).and.(pub_paw_output_detail>=VERBOSE)) then
             if (pub_on_root) write(stdout,'(a,i4)') 'dijxc', is
             call sparse_show_matrix(dij_xc(is),show_elems=.true.)
          end if
          call sparse_axpy(dij(is),dij_xc(is),1.0_DP)
          edijxc = edijxc + sparse_trace(rho_ij(is),dij_xc(is))
       end do
    else
       exc = 0.0_DP; exc_dc = 0.0_DP
       etxc = 0.0_DP; etxc_dc = 0.0_DP
       exc_core = 0.0_DP; edijxc = 0.0_DP
    end if
    if (present(paw_sphere_energies)) then
       paw_sphere_energies(paw_en_exc) = exc
       paw_sphere_energies(paw_en_exc_dc) = exc_dc
       paw_sphere_energies(paw_en_etxc) = etxc
       paw_sphere_energies(paw_en_etxc_dc) = etxc_dc
       paw_sphere_energies(paw_en_dijxc) = edijxc
       paw_sphere_energies(paw_en_exc_core) = exc_core
    end if

    ! Clean up temporary matrices
    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(dij_xc(is))
       call sparse_destroy(dij_hartree(is))
       call sparse_destroy(dij0(is))
    end do
    deallocate(dij_xc,stat=ierr)
    call utils_dealloc_check('paw_nonlocal_energies','dij_xc',ierr)
    deallocate(dij_hartree,stat=ierr)
    call utils_dealloc_check('paw_nonlocal_energies','dij_hartree',ierr)
    deallocate(dij0,stat=ierr)
    call utils_dealloc_check('paw_nonlocal_energies','dij0',ierr)

  end subroutine paw_nonlocal_energies


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dij0(dij0)

    !=================================================================!
    ! This subroutine calculates the constant terms in the nonlocal   !
    ! energies Dij, which are given in the PAW dataset as Dij0. These !
    ! result from the kinetic energy of the sphere parts of the all-  !
    ! electron wavefunctions, and the interaction of these with the   !
    ! core Hartree potential.                                         !
    !-----------------------------------------------------------------!
    ! Arguments:                                                      !
    !   dij0(inout) : Constant contribution D_ij^0 to nonlocal energy !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 28/05/10.                           !
    !=================================================================!

    use comms, only: pub_my_node_id
    use parallel_strategy, only: pub_num_atoms_on_node, &
         pub_first_atom_on_node, pub_elements_on_node
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_put_block
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: dij0(pub_cell%num_spins)

    ! Local Variables
    integer :: loc_iat, iat
    integer :: isp
    integer :: is
    integer :: ipwtot,jpwtot
    integer :: ierr
    real(kind=DP), allocatable :: dij0_sp(:,:,:)
    complex(kind=DP), allocatable :: dij0_sp_cmplx(:,:,:)

    allocate(dij0_sp(max_paw_proj_tot,max_paw_proj_tot, &
         pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('paw_dij0','dij0_sp',ierr)

    ! Loop over species, setting up array of dij0 terms
    do isp=1,pub_cell%num_pspecies
       if (any(pub_elements_on_node(:)%pspecies_number==isp)) then
          do ipwtot=1,paw_sp(isp)%npw_tot
             do jpwtot=1,paw_sp(isp)%npw_tot
                dij0_sp(ipwtot,jpwtot,isp) = paw_sp(isp)%dij0(ipwtot,jpwtot)
             end do
          end do
       end if
    end do

    if (dij0(1)%iscmplx) then
       allocate(dij0_sp_cmplx(max_paw_proj_tot,max_paw_proj_tot, &
            pub_cell%num_pspecies),stat=ierr)
       call utils_alloc_check('paw_dij0','dij0_sp',ierr)
       dij0_sp_cmplx(:,:,:) = cmplx(dij0_sp(:,:,:),kind=DP)
    end if

    ! Loop over atoms
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = loc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
       isp = pub_elements_on_node(loc_iat)%pspecies_number
       ! Put block of dij0 into SPAM3 matrix
       do is=1,pub_cell%num_spins
          if (dij0(1)%iscmplx) then
             call sparse_put_block(dij0_sp_cmplx(:,:,isp),dij0(is),iat,iat)
          else
             call sparse_put_block(dij0_sp(:,:,isp),dij0(is),iat,iat)
          end if
       end do
    end do

    if (dij0(1)%iscmplx) then
       deallocate(dij0_sp_cmplx,stat=ierr)
       call utils_dealloc_check('paw_dij0','dij0_sp_cmplx',ierr)
    end if

    deallocate(dij0_sp,stat=ierr)
    call utils_dealloc_check('paw_dij0','dij0_sp',ierr)

  end subroutine paw_dij0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dij_hartree(dijhartree,rho_ij)

    !=================================================================!
    ! This subroutine calculates the Hartree terms in the nonlocal    !
    ! energies Dij, from the tensor e_ijkl calculated from the        !
    ! information in the PAW dataset, in paw_init_species.            !
    !-----------------------------------------------------------------!
    ! Arguments:                                                      !
    !   dijhartree(inout) : Hartree contribution to nonlocal energy   !
    !                       term hat{D}_ij.                           !
    !   rho_ij(in)        : Projector density kernel rho^ij           !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 28/05/10.                           !
    !=================================================================!

    use comms, only: pub_my_node_id
    use parallel_strategy, only: pub_num_atoms_on_node, pub_elements_on_node, &
         pub_first_atom_on_node
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_get_block, sparse_put_block
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: dijhartree(pub_cell%num_spins)
    type(SPAM3),intent(in) :: rho_ij(pub_cell%num_spins)

    ! Local Variables
    integer :: loc_iat, iat, ierr
    integer :: isp
    integer :: is
    real(kind=DP), allocatable :: rhoij_at(:,:,:), dijh_at(:,:)

    ! Start Timer
    call timer_clock('paw_dij_hartree',1)

    allocate(rhoij_at(max_paw_proj_tot,max_paw_proj_tot, &
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_hartree','rhoij_at',ierr)
    allocate(dijh_at(max_paw_proj_tot,max_paw_proj_tot),stat=ierr)
    call utils_alloc_check('paw_dij_hartree','dijh_at',ierr)
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = loc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
       isp = pub_elements_on_node(loc_iat)%pspecies_number

       ! Get rhoij_at for this atom (summed over spins if necessary
       rhoij_at(:,:,:) = 0.0_DP
       do is=1,pub_cell%num_spins
          call sparse_get_block(rhoij_at(:,:,is),rho_ij(is),iat,iat)
       end do
       if (pub_cell%num_spins==2) then
          rhoij_at(:,:,1) = rhoij_at(:,:,1) + rhoij_at(:,:,2)
       end if

       ! Calculate Hartree terms in nonlocal energies for this atom
       dijh_at(:,:) = 0.0_DP
       call paw_dij_hartree_atom(isp,max_paw_proj_tot,rhoij_at,dijh_at)
       do is=1,pub_cell%num_spins
          call sparse_put_block(dijh_at,dijhartree(is),iat,iat)
       end do

    end do

    deallocate(dijh_at,stat=ierr)
    call utils_dealloc_check('paw_dij_hartree','dijh_at',ierr)
    deallocate(rhoij_at,stat=ierr)
    call utils_dealloc_check('paw_dij_hartree','rhoij_at',ierr)

    ! Stop Timer
    call timer_clock('paw_dij_hartree',2)

  end subroutine paw_dij_hartree


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dij_hartree_atom(isp,npwdim,rhoij,dijh)

    !=================================================================!
    ! This subroutine calculates the Hartree terms in the nonlocal    !
    ! energies Dij for a given atom, from the tensor e_ijkl and the   !
    ! projector density matrix rhoij.                                 !
    !-----------------------------------------------------------------!
    ! Arguments:                                                      !
    !   dijh(inout)   : Compensation density contribution to nonlocal !
    !                   energy term \hat{D}_ij.                       !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 28/05/10.                           !
    !=================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: isp
    integer, intent(in) :: npwdim
    real(kind=DP), intent(in) :: rhoij(npwdim,npwdim)
    real(kind=DP), intent(out) :: dijh(npwdim,npwdim)

    ! Local Variables
    integer :: ipwtot, jpwtot, kpwtot, lpwtot

    dijh(:,:) = 0.0_DP
    do ipwtot=1,paw_sp(isp)%npw_tot
       do jpwtot=1,paw_sp(isp)%npw_tot
          do kpwtot=1,paw_sp(isp)%npw_tot
             do lpwtot=1,paw_sp(isp)%npw_tot
                dijh(ipwtot,jpwtot) = dijh(ipwtot,jpwtot) + &
                     paw_sp(isp)%e_ijkl(ipwtot,jpwtot,kpwtot,lpwtot) * &
                     rhoij(kpwtot,lpwtot)
             end do
          end do
       end do
    end do

  end subroutine paw_dij_hartree_atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dij_xc(dijxc,rho_ij,exc,exc_dc,etxc,etxc_dc)

    !=================================================================!
    ! This subroutine calculates the XC terms in the nonlocal         !
    ! energies Dij, and also the XC energies and double-counting      !
    ! terms.                                                          !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 28/05/10.                           !
    !=================================================================!

    use comms, only: comms_reduce, pub_my_node_id, pub_on_root
    use constants, only: max_spins, PI
    use parallel_strategy, only: pub_num_atoms_on_node, pub_elements_on_node, &
       pub_first_atom_on_node
    use rundat, only: pub_paw_output_detail, pub_aug_funcs_recip
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_get_block, sparse_put_block
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: dijxc(pub_cell%num_spins)
    type(SPAM3),intent(in) :: rho_ij(pub_cell%num_spins)
    real(kind=DP),intent(out) :: exc
    real(kind=DP),intent(out) :: exc_dc
    real(kind=DP),intent(out) :: etxc
    real(kind=DP),intent(out) :: etxc_dc

    ! Local Variables
    character(20) :: fmt,tmp
    real(kind=DP), allocatable :: vxc_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: density_rad_LM(:,:,:,:)
    real(kind=DP), allocatable :: rhoij_at(:,:,:)
    real(kind=DP), allocatable :: dijxc_at(:,:,:)
    real(kind=DP), allocatable :: dijtxc_at(:,:,:)
    real(kind=DP), allocatable :: pot_work(:,:,:,:)
    real(kind=DP), allocatable :: den_work(:,:,:,:)
    real(kind=DP), allocatable :: inter(:)
    integer :: loc_iat
    integer :: iat
    integer :: ierr
    integer :: is
    integer :: isp
    integer :: igrid
    integer :: ipwtot,jpwtot
    integer :: npts_max
    integer :: nLM_max
    real(kind=DP) :: exc_at
    real(kind=DP) :: exc_dc_at
    real(kind=DP) :: etxc_at
    real(kind=DP) :: etxc_dc_at
    real(kind=DP) :: edijxc_at
    real(kind=DP) :: edijtxc_at
    real(kind=DP) :: total_nhat(max_spins)
    real(kind=DP) :: total_nhat_at(max_spins)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering paw_dij_xc'
#endif

    ! Start Timer
    call timer_clock('paw_dij_xc',1)

    ! Find array sizes
    npts_max = 0
    nLM_max = 0
    do isp=1,pub_cell%num_pspecies
       if (.not.any(pub_elements_on_node(:)%pspecies_number==isp)) cycle
       igrid = paw_sp(isp)%phi_grid
       npts_max = max(npts_max,paw_sp(isp)%grid(igrid)%npt)
       igrid = paw_sp(isp)%shape_grid
       npts_max = max(npts_max,paw_sp(isp)%grid(igrid)%npt)
       igrid = paw_sp(isp)%core_den_grid
       npts_max = max(npts_max,paw_sp(isp)%grid(igrid)%npt)
       nLM_max = max(nLM_max,(paw_sp(isp)%lmax+1)**2)
    end do

    ! Allocate temporary arrays
    allocate(vxc_rad_LM(npts_max,pub_cell%num_spins,nLM_max,2),stat=ierr)
    call utils_alloc_check('paw_dij_xc','vxc_rad_LM',ierr)
    allocate(density_rad_LM(npts_max,pub_cell%num_spins,nLM_max,3),stat=ierr)
    call utils_alloc_check('paw_dij_xc','density_rad_LM',ierr)
    allocate(rhoij_at(max_paw_proj_tot,max_paw_proj_tot, &
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_xc','rhoij_at',ierr)
    allocate(dijxc_at(max_paw_proj_tot,max_paw_proj_tot, &
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_xc','dijxc_at',ierr)
    allocate(dijtxc_at(max_paw_proj_tot,max_paw_proj_tot, &
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('paw_dij_xc','dijtxc_at',ierr)
    allocate(pot_work(npts_max,pub_cell%num_spins,pub_cell%num_spins,6), &
         stat=ierr)
    call utils_alloc_check('paw_dij_xc','pot_work',ierr)
    allocate(den_work(npts_max,pub_cell%num_spins,pub_cell%num_spins,4), &
         stat=ierr)
    call utils_alloc_check('paw_dij_xc','den_work',ierr)
    allocate(inter(npts_max),stat=ierr)
    call utils_alloc_check('paw_dij_xc','inter',ierr)

    ! Initialisations
    exc = 0.0_DP
    exc_dc = 0.0_DP
    etxc = 0.0_DP
    etxc_dc = 0.0_DP
    total_nhat(:) = 0.0_DP

    ! Loop over atoms
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = loc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
       isp = pub_elements_on_node(loc_iat)%pspecies_number

       ! Get block of projector density matrix for this atom
       do is=1,pub_cell%num_spins
          call sparse_get_block(rhoij_at(:,:,is),rho_ij(is),iat,iat)
       end do  ! is

       ! Calculate exchange-correlation energy and potential for this atom
       call paw_dij_xc_atom(isp,pub_cell%num_spins,npts_max,nLM_max, &
            max_paw_proj_tot,rhoij_at,dijxc_at,dijtxc_at, &
            exc_at,etxc_at,exc_dc_at,etxc_dc_at,total_nhat_at, &
            density_rad_LM,vxc_rad_LM,den_work,pot_work,inter)

       ! Add up contributions to energy from this atom
       exc = exc + exc_at
       etxc = etxc + etxc_at
       exc_dc = exc_dc + exc_dc_at
       etxc_dc = etxc_dc + etxc_dc_at
       total_nhat(:) = total_nhat(:) + total_nhat_at(:)

       ! Create total dijxc and add it to matrix (then undo total)
       do is=1,pub_cell%num_spins
          dijxc_at(:,:,is) = dijxc_at(:,:,is) - dijtxc_at(:,:,is)
          call sparse_put_block(dijxc_at(:,:,is),dijxc(is),iat,iat)
          dijxc_at(:,:,is) = dijxc_at(:,:,is) + dijtxc_at(:,:,is)
       end do  ! is

       ! Add up AE dijxc energy \sum_ij \rho_ij D^xc_ij and
       ! and PS dijtxc energy \sum_ij \rho_ij \tilde{D}^xc_ij
       edijxc_at = 0.0_DP
       edijtxc_at = 0.0_DP
       do is=1,pub_cell%num_spins
          do ipwtot=1,paw_sp(isp)%npw_tot
             do jpwtot=1,paw_sp(isp)%npw_tot
                edijxc_at = edijxc_at + rhoij_at(ipwtot,jpwtot,is) * &
                     dijxc_at(ipwtot,jpwtot,is)
                edijtxc_at = edijtxc_at + rhoij_at(ipwtot,jpwtot,is) * &
                     dijtxc_at(ipwtot,jpwtot,is)
             end do
          end do
       end do

       ! Consistency Check
       if (abs(edijxc_at-exc_dc_at)>0.0000001_DP*abs(exc_dc_at)) then
          if (pub_on_root) then
             write(stdout,'(a,i6,a)') 'For atom ',iat,':'
             write(stdout,'(a,f20.12)') 'edijxc_at: ',edijxc_at
             write(stdout,'(a,f20.12)') 'exc_dc_at: ',exc_dc_at
          end if
          call utils_abort('Error in paw_dij_xc: Consistency check failed: &
               &edijxc_at /= exc_dc_at')
       end if
       if (abs(edijtxc_at-etxc_dc_at)>0.0000001_DP*abs(etxc_dc_at)) then
          if (pub_on_root) then
             write(stdout,'(a,i6,a)') 'For atom ',iat,':'
             write(stdout,'(a,f20.12)') 'edijtxc_at: ',edijtxc_at
             write(stdout,'(a,f20.12)') 'etxc_dc_at: ',etxc_dc_at
          end if
          call utils_abort('Error in paw_dij_xc: Consistency check failed: &
               &edijtxc_at /= etxc_dc_at')
       end if

    end do  ! loc_iat

    ! Sum energies across all nodes
    call comms_reduce('SUM',exc)
    call comms_reduce('SUM',exc_dc)
    call comms_reduce('SUM',etxc)
    call comms_reduce('SUM',etxc_dc)

    ! Check and show total compensation charge on radial grid
    call comms_reduce('SUM',total_nhat(:))
    if (pub_on_root.and.(pub_paw_output_detail>=VERBOSE).and. &
         (.not.pub_aug_funcs_recip)) then
       write(tmp,'(i5)') pub_cell%num_spins
       write(fmt,'(3a)') '(a,',trim(adjustl(tmp)),'f14.8)'
       write(stdout,fmt) 'Total Compensation Charge on Radial Grid : ',&
            total_nhat(1:pub_cell%num_spins)
    end if

    ! Deallocate temporary arrays
    deallocate(inter,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','inter',ierr)
    deallocate(den_work,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','den_work',ierr)
    deallocate(pot_work,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','pot_work',ierr)
    deallocate(dijtxc_at,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','dijtxc_at',ierr)
    deallocate(dijxc_at,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','dijxc_at',ierr)
    deallocate(rhoij_at,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','rhoij_at',ierr)
    deallocate(density_rad_LM,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','density_rad_LM',ierr)
    deallocate(vxc_rad_LM,stat=ierr)
    call utils_dealloc_check('paw_dij_xc','vxc_rad_LM',ierr)

    ! Stop Timer
    call timer_clock('paw_dij_xc',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving paw_dij_xc'
#endif

  end subroutine paw_dij_xc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_dij_xc_atom(isp,nspins,npts_max,nLM_max,npw, &    ! in
       rhoij,dijxc,dijtxc,exc,etxc,exc_dc,etxc_dc,total_nhat_at, & ! out
       density_rad_LM,vxc_rad_LM,den_work,pot_work,inter)          ! workspace

    !=================================================================!
    ! This subroutine calculates the xc terms in the nonlocal         !
    ! energies Dij for a given atom, using the projector density      !
    ! matrix rhoij.                                                   !
    !-----------------------------------------------------------------!
    ! Written by Nicholas Hine on 28/05/10.                           !
    !=================================================================!

    use comms, only: pub_on_root
    use constants, only: max_spins, PI
    use gaunt_coeff, only: realgaunt
    use paw_xc, only: paw_xc_pot_rad_LM, paw_xc_exc_dc
    use rundat, only: pub_paw_output_detail, pub_nhat_in_xc
    use services, only: services_radial_integral_rmax
    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: isp
    integer, intent(in) :: npw
    integer, intent(in) :: nspins
    integer, intent(in) :: npts_max, nLM_max
    real(kind=DP), intent(in) :: rhoij(npw,npw,nspins)
    real(kind=DP), intent(out) :: dijxc(npw,npw,nspins)
    real(kind=DP), intent(out) :: dijtxc(npw,npw,nspins)
    real(kind=DP), intent(out) :: exc, etxc
    real(kind=DP), intent(out) :: exc_dc, etxc_dc
    real(kind=DP), intent(out) :: density_rad_LM(npts_max,nspins,nLM_max,3)
    real(kind=DP), intent(out) :: vxc_rad_LM(npts_max,nspins,nLM_max,2)
    real(kind=DP), intent(out) :: den_work(npts_max,nspins,nspins,4)
    real(kind=DP), intent(out) :: pot_work(npts_max,nspins,nspins,6)
    real(kind=DP), intent(out) :: inter(npts_max)
    real(kind=DP), intent(out) :: total_nhat_at(max_spins)

    ! Local Variables
    integer :: ipwtot, jpwtot
    integer :: igrid, is
    integer :: lup, mup
    integer :: ipw,jpw
    integer :: li,mi,lj,mj
    integer :: npts_shp
    integer :: npts_phi
    integer :: npts_den
    integer :: nLM
    integer :: iLM
    integer :: lmax
    real(kind=DP) :: rgij
    real(kind=DP) :: total_nhat_at_targ

    ! Find numbers of grid points and angular momenta for this atom
    igrid = paw_sp(isp)%phi_grid
    npts_phi = paw_sp(isp)%grid(igrid)%npt
    igrid = paw_sp(isp)%shape_grid
    npts_shp = paw_sp(isp)%grid(igrid)%npt
    igrid = paw_sp(isp)%core_den_grid
    npts_den = paw_sp(isp)%grid(igrid)%npt
    npts_den = max(npts_den,npts_shp,npts_phi)
    nLM = (paw_sp(isp)%lmax + 1)**2
    lmax = paw_sp(isp)%lmax

    ! Reset dijxc blocks for this atom
    dijxc(:,:,:) = 0.0_DP
    dijtxc(:,:,:) = 0.0_DP
    exc_dc = 0.0_DP
    etxc_dc = 0.0_DP

    ! Create (n^1 + n_c)(r) in density(1)
    ! Create (\tilde{n}^1 + \tilde{n}_c)(r) in density(2)
    ! Create \hat{n}  in density(3)
    call paw_ae_density_rad_LM(density_rad_LM(:,:,:,1),isp,rhoij, &
         npw,npts_max,nspins,nLM_max,den_work)
    call paw_ps_density_rad_LM(density_rad_LM(:,:,:,2),isp,rhoij, &
         npw,npts_max,nspins,nLM_max,den_work)
    call paw_aug_density_rad_LM(density_rad_LM(:,:,:,3),isp,rhoij, &
         npw,npts_max,nspins,nLM_max,den_work)

    ! Calculate total compensation density and compare to target
    do is=1,nspins
       den_work(1:npts_den,1,1,1) = density_rad_LM(1:npts_den,is,1,3) * &
            paw_sp(isp)%grid(igrid)%r(1:npts_den)**2 * inv_sqrt_4pi
       total_nhat_at(is) = services_radial_integral_rmax(npts_den, &
            paw_sp(isp)%grid(igrid)%rab,paw_sp(isp)%grid(igrid)%r, &
            paw_sp(isp)%rcut,den_work(:,1,1,1),inter) * 4.0_DP*PI
       total_nhat_at_targ = 0.0_DP
       do ipwtot=1,paw_sp(isp)%npw_tot
          ipw = paw_sp(isp)%ipw_tot(ipwtot)
          li = paw_sp(isp)%l_pw_tot(ipwtot)
          mi = paw_sp(isp)%m_pw_tot(ipwtot)
          do jpwtot=1,paw_sp(isp)%npw_tot
             jpw = paw_sp(isp)%ipw_tot(jpwtot)
             lj = paw_sp(isp)%l_pw_tot(jpwtot)
             mj = paw_sp(isp)%m_pw_tot(jpwtot)

             total_nhat_at_targ = total_nhat_at_targ + sqrt(4.0_DP*PI) * &
                  realgaunt(0,0,li,mi,lj,mj) * rhoij(ipwtot,jpwtot,is) * &
                  paw_sp(isp)%aug_nLij(ipw,jpw,0)
          end do  ! jpwtot
       end do  ! ipwtot

       ! Check total compensation density does not differ by more than 0.1%
       ! from its intended value
       if (abs(total_nhat_at_targ-total_nhat_at(is)) > &
            0.001_DP*abs(total_nhat_at_targ)) then
            call utils_abort('Error in paw_dij_xc_atom: Total Compensation &
                 &charge on radial grid does not match target')
       end if
    end do  ! is

    ! Print AE density
    !call print_rad(density_rad_LM(:,:,:,1),paw_sp(isp)%grid(igrid)%r, &
    !     npts_den,nspins,nLM,'n^1+n_c')

    ! Calculate v_xc[density(1)]
    call paw_xc_pot_rad_LM(vxc_rad_LM(:,:,:,1),exc, &
         density_rad_LM(:,:,:,1), &
         paw_sp(isp)%grid(igrid)%r(:),paw_sp(isp)%grid(igrid)%rab(:), &
         paw_sp(isp)%rcut,npts_phi,npts_max,nspins,nLM_max,lmax, &
         den_work,pot_work,inter(:))

    ! Subtract off AE core density and calculate double-counting term
    do is=1,nspins
       density_rad_LM(1:npts_den,is,1,1) = density_rad_LM(1:npts_den,is,1,1) &
            - paw_sp(isp)%core_den_rad(1:npts_den) * sqrt_4pi &
            / real(nspins,kind=DP)
    end do

    ! Calculate Double Counting correction to AE XC energy
    call paw_xc_exc_dc(exc_dc,vxc_rad_LM(:,:,:,1), &
         density_rad_LM(:,:,:,1),paw_sp(isp)%grid(igrid)%r(:), &
         paw_sp(isp)%grid(igrid)%rab(:),paw_sp(isp)%rcut, &
         npts_phi,npts_max,nspins,nLM,nLM_max,den_work(:,1,1,1),inter(:))

    ! Print AE potential
    !call print_rad(vxc_rad_LM(:,:,:,1),paw_sp(isp)%grid(igrid)%r, &
    !     npts_den,nspins,nLM,'vxc(n^1+n_c)')

    ! Add compensation density nhat to valence density \tilde{n}^1
    ! (but only if we are using the nhat density in the xc terms)
    if (pub_nhat_in_xc) then
       density_rad_LM(1:npts_den,1:nspins,1:nLM,2) = &
            density_rad_LM(1:npts_den,1:nspins,1:nLM,2) + &
            density_rad_LM(1:npts_den,1:nspins,1:nLM,3)
    end if

    ! Print soft PS density
    !call print_rad(density_rad_LM(:,:,:,2),paw_sp(isp)%grid(igrid)%r, &
    !     npts_den,nspins,nLM,'tn^1+nhat+tn_c')

    ! Calculate v_xc[density(2)]
    call paw_xc_pot_rad_LM(vxc_rad_LM(:,:,:,2),etxc, &
         density_rad_LM(:,:,:,2), &
         paw_sp(isp)%grid(igrid)%r(:),paw_sp(isp)%grid(igrid)%rab(:), &
         paw_sp(isp)%rcut,npts_phi,npts_max,nspins,nLM_max,lmax, &
         den_work,pot_work,inter(:))

    ! Subtract off PS core density and calculate double-counting term
    do is=1,nspins
       density_rad_LM(1:npts_den,is,1,2) = density_rad_LM(1:npts_den,is,1,2) &
            - paw_sp(isp)%tcore_den_rad(1:npts_den) * sqrt_4pi &
            / real(nspins,kind=DP)
    end do

    ! Calculate Double Counting correction to PS XC energy
    call paw_xc_exc_dc(etxc_dc,vxc_rad_LM(:,:,:,2), &
         density_rad_LM(:,:,:,2),paw_sp(isp)%grid(igrid)%r(:), &
         paw_sp(isp)%grid(igrid)%rab(:),paw_sp(isp)%rcut, &
         npts_phi,npts_max,nspins,nLM,nLM_max,den_work(:,1,1,1),inter(:))

    ! Print PS potential
    !call print_rad(vxc_rad_LM(:,:,:,2),paw_sp(isp)%grid(igrid)%r, &
    !     npts_den,nspins,nLM,'vxc(tn^1+nhat+tn_c)')

    ! Calculate matrix elements of dij_xc
    do is=1,nspins
       do lup=0,lmax
          do mup=-lup,lup
             iLM = lup**2 + lup + 1 + mup

             ! Double loop over projectors ipwtot,jpwtot
             do ipwtot=1,paw_sp(isp)%npw_tot
                ipw = paw_sp(isp)%ipw_tot(ipwtot)
                li = paw_sp(isp)%l_pw_tot(ipwtot)
                mi = paw_sp(isp)%m_pw_tot(ipwtot)
                do jpwtot=1,paw_sp(isp)%npw_tot
                   jpw = paw_sp(isp)%ipw_tot(jpwtot)
                   lj = paw_sp(isp)%l_pw_tot(jpwtot)
                   mj = paw_sp(isp)%m_pw_tot(jpwtot)

                   ! Find Gaunt coefficient for this combination of indices
                   ! G = G^LM_{l_i m_i l_j m_j}
                   rgij = realgaunt(lup,mup,li,mi,lj,mj)

                   ! Term 1: \sum_LM G \int_0^{r_c} dr
                   !           v_xc[n^1 + n_c](r) \phi_i(r) \phi_j(r)
                   den_work(1:npts_phi,1,1,1) = &
                        paw_sp(isp)%phi_rad(1:npts_phi,ipw) * &
                        paw_sp(isp)%phi_rad(1:npts_phi,jpw) * &
                        vxc_rad_LM(1:npts_phi,is,iLM,1)
                   dijxc(ipwtot,jpwtot,is) = &
                        dijxc(ipwtot,jpwtot,is) + rgij * &
                        services_radial_integral_rmax(npts_phi, &
                        paw_sp(isp)%grid(igrid)%rab, &
                        paw_sp(isp)%grid(igrid)%r,paw_sp(isp)%rcut, &
                        den_work(1:npts_phi,1,1,1),inter)

                   ! Term 2: -\sum_LM G \int_0^{r_c} dr
                   !           v_xc[\tilde{n}^1 + \hat{n}* + \tilde{n}_c](r)
                   !          \tilde{\phi}_i(r) \tilde{\phi}_j(r)
                   ! * only if pub_nhat_in_xc is true
                   den_work(1:npts_phi,1,1,2) = &
                        paw_sp(isp)%tphi_rad(1:npts_phi,ipw) * &
                        paw_sp(isp)%tphi_rad(1:npts_phi,jpw) * &
                        vxc_rad_LM(1:npts_phi,is,iLM,2)
                   dijtxc(ipwtot,jpwtot,is) = &
                        dijtxc(ipwtot,jpwtot,is) + rgij * &
                        services_radial_integral_rmax(npts_phi, &
                        paw_sp(isp)%grid(igrid)%rab, &
                        paw_sp(isp)%grid(igrid)%r,paw_sp(isp)%rcut, &
                        den_work(1:npts_phi,1,1,2),inter)

                   ! Last part not present if nhat is not included in XC
                   if (pub_nhat_in_xc) then

                      ! Term 3: -\sum_LM G \int_0^{r_c} dr
                      !           v_xc[\tilde{n}^1 + \hat{n} + \tilde{n}_c](r)
                      !           n_{n_i l_i n_j l_j}^L g_L(r) * r^2
                      den_work(1:npts_phi,1,1,3) = &
                           paw_sp(isp)%aug_nLij(ipw,jpw,lup) * &
                           paw_sp(isp)%shape%shape_rad(1:npts_phi,lup) * &
                           vxc_rad_LM(1:npts_phi,is,iLM,2) * &
                           paw_sp(isp)%grid(igrid)%r(1:npts_phi)**2
                      dijtxc(ipwtot,jpwtot,is) = &
                           dijtxc(ipwtot,jpwtot,is) + rgij * &
                           services_radial_integral_rmax(npts_phi, &
                           paw_sp(isp)%grid(igrid)%rab, &
                           paw_sp(isp)%grid(igrid)%r,paw_sp(isp)%rcut, &
                           den_work(1:npts_phi,1,1,3),inter)

                   end if

                end do  ! jpwtot
             end do  ! ipwtot
          end do  ! mup
       end do  ! lup

    end do  ! is

    return

contains

    subroutine print_rad(rad_fun,rad,npt,ns,ncomp,name)

       ! Arguments
       integer, intent(in) :: npt
       integer, intent(in) :: ns
       integer, intent(in) :: ncomp
       real(kind=DP), intent(in) :: rad_fun(npt,ns,ncomp)
       real(kind=DP), intent(in) :: rad(npt)
       character(*), intent(in) :: name

       ! Local Variables
       integer :: ipt
       integer :: is
       integer :: icomp

       do is=1,ns
          write(stdout,*)
          write(stdout,*) name
          do icomp=1,ncomp
             do ipt=1,npts_den
                write(stdout,'(2i6,2f22.12)') icomp,ipt,rad(ipt),rad_fun(ipt,is,icomp)
             end do
             write(stdout,*)
          end do
       end do

    end subroutine print_rad

  end subroutine paw_dij_xc_atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_exc_core_tot(exc_core)

    !==================================================================!
    ! This subroutine adds up the core density XC energies for all the !
    ! atoms in the system                                              !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/05/10.                            !
    !==================================================================!

    use comms, only: comms_reduce, pub_my_node_id
    use parallel_strategy, only: pub_elements_on_node, pub_first_atom_on_node, &
         pub_num_atoms_on_node
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: exc_core

    ! Local Variables
    integer :: isp
    integer :: loc_iat
    integer :: iat

    ! Calculate exc core for each species if not already done
    do isp=1,pub_cell%num_pspecies
       if (.not.paw_sp(isp)%core_charge_calculated) then
          call paw_exc_core_atom(paw_sp(isp)%exc_core,isp)
          paw_sp(isp)%core_charge_calculated = .true.
       end if
    end do

    ! Add up exc core over all atoms in system
    exc_core = 0.0_DP
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = loc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
       isp = pub_elements_on_node(loc_iat)%pspecies_number
       exc_core = exc_core + paw_sp(isp)%exc_core
    end do
    call comms_reduce('SUM',exc_core)

  end subroutine paw_exc_core_tot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_exc_core_atom(exc_core_atom,isp)

    !==================================================================!
    ! This subroutine calculates the core density XC energies for a    !
    ! Single species of atom.                                          !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/05/10.                            !
    !==================================================================!

    use comms, only: comms_reduce, pub_my_node_id
    use parallel_strategy, only: pub_elements_on_node, pub_first_atom_on_node, &
         pub_num_atoms_on_node
    use paw_xc, only: paw_xc_pot_rad_LM
    use rundat, only: xc_functional
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: exc_core_atom
    integer, intent(in) :: isp

    ! Local Variables
    integer :: ierr
    integer :: is
    integer :: npts
    integer :: nspins
    integer :: igrid
    real(kind=DP), allocatable :: core_den(:,:,:)
    real(kind=DP), allocatable :: core_den_vxc(:,:,:)
    real(kind=DP), allocatable :: xc_pot_work(:,:,:,:)
    real(kind=DP), allocatable :: xc_den_work(:,:,:,:)
    real(kind=DP), allocatable :: inter(:)

    if (paw_sp(isp)%core_charge) then
       igrid = paw_sp(isp)%core_den_grid
       npts = paw_sp(isp)%grid(igrid)%npt
       nspins = 1
       allocate(core_den(npts,nspins,1),stat=ierr)
       call utils_alloc_check('paw_exc_core_atom','core_den',ierr)
       allocate(core_den_vxc(npts,nspins,1),stat=ierr)
       call utils_alloc_check('paw_exc_core_atom','core_den_vxc',ierr)
       allocate(xc_pot_work(npts,nspins,nspins,6),stat=ierr)
       call utils_alloc_check('paw_exc_core_atom','xc_pot_work',ierr)
       allocate(xc_den_work(npts,nspins,nspins,4),stat=ierr)
       call utils_alloc_check('paw_exc_core_atom','xc_den_work',ierr)
       allocate(inter(npts),stat=ierr)
       call utils_alloc_check('paw_exc_core_atom','inter',ierr)

       do is=1,nspins
          core_den(1:npts,is,1) = paw_sp(isp)%core_den_rad(1:npts) &
               * sqrt_4pi / real(nspins,kind=DP)
       end do

       ! Call radial XC function with just nLM_max = 1 (no moment expansion)
       call paw_xc_pot_rad_LM(core_den_vxc(:,:,:),exc_core_atom, &
            core_den(:,:,:),paw_sp(isp)%grid(igrid)%r, &
            paw_sp(isp)%grid(igrid)%rab,paw_sp(isp)%rcut,npts,npts, &
            nspins,1,0,xc_den_work,xc_pot_work,inter)

       deallocate(inter,stat=ierr)
       call utils_dealloc_check('paw_exc_core_atom','inter',ierr)
       deallocate(xc_den_work,stat=ierr)
       call utils_dealloc_check('paw_exc_core_atom','xc_den_work',ierr)
       deallocate(xc_pot_work,stat=ierr)
       call utils_dealloc_check('paw_exc_core_atom','xc_pot_work',ierr)
       deallocate(core_den_vxc,stat=ierr)
       call utils_dealloc_check('paw_exc_core_atom','core_den_vxc',ierr)
       deallocate(core_den,stat=ierr)
       call utils_dealloc_check('paw_exc_core_atom','core_den',ierr)

    else

       paw_sp(isp)%exc_core = 0.0_DP
       paw_sp(isp)%core_charge_calculated = .true.

    end if

  end subroutine paw_exc_core_atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_ae_density_rad_LM(ae_den,isp,rhoij_at,npw,npts_max,nspins, &
       nLM_max,den_work)

    !==================================================================!
    ! This subroutine calculates the all-electron sphere density in    !
    ! the PAW method, as a sum over the moments LM.                    !
    !  ae_den(:,is,iLM) = n^1 density + core density n_c               !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/05/10.                            !
    ! Split off to a separate routine by Nicholas Hine on 09/05/11.    !
    !==================================================================!

    use gaunt_coeff, only: realgaunt
    use services, only: services_analytic_limit

    implicit none

    ! Arguments
    integer, intent(in) :: npts_max
    integer, intent(in) :: nspins
    integer, intent(in) :: npw
    integer, intent(in) :: nLM_max
    real(kind=DP), intent(out) :: ae_den(npts_max,nspins,nLM_max)
    real(kind=DP), intent(out) :: den_work(npts_max)
    real(kind=DP), intent(in) :: rhoij_at(npw,npw,nspins)
    integer, intent(in) :: isp

    ! Local Variables
    integer :: iLM
    integer :: npts_den
    integer :: npts_phi
    integer :: kpwtot,lpwtot
    integer :: kpw,lpw
    integer :: lk,mk,ll,ml
    integer :: lup,mup
    integer :: is
    integer :: igrid
    real(kind=DP) :: rgkl
    real(kind=DP),parameter :: gaunt_tol = 1.0e-14_DP
    real(kind=DP),parameter :: rho_tol = 1.0e-14_DP

    ! Initialise, and find information about this species and its phi_grid
    ae_den(:,:,:) = 0.0_DP
    igrid = paw_sp(isp)%phi_grid
    npts_phi = paw_sp(isp)%grid(igrid)%npt
    igrid = paw_sp(isp)%core_den_grid
    npts_den = paw_sp(isp)%grid(igrid)%npt

    ! Double loop over partial waves on this atom
    do kpwtot=1,paw_sp(isp)%npw_tot
       kpw = paw_sp(isp)%ipw_tot(kpwtot)
       lk = paw_sp(isp)%l_pw_tot(kpwtot)
       mk = paw_sp(isp)%m_pw_tot(kpwtot)
       do lpwtot=1,paw_sp(isp)%npw_tot
          lpw = paw_sp(isp)%ipw_tot(lpwtot)
          ll = paw_sp(isp)%l_pw_tot(lpwtot)
          ml = paw_sp(isp)%m_pw_tot(lpwtot)

          if (all(abs(rhoij_at(kpwtot,lpwtot,:))<rho_tol)) cycle

          ! den_work(:) = tphi_i(r)*tphi_j(r) / r^2
          den_work(2:npts_phi) = paw_sp(isp)%phi_rad(2:npts_phi,kpw) &
               * paw_sp(isp)%phi_rad(2:npts_phi,lpw) &
               / paw_sp(isp)%grid(igrid)%r(2:npts_phi)**2
          ! Find analytic limit
          call services_analytic_limit(npts_phi,paw_sp(isp)%grid(igrid)%r(:), &
               den_work(:))

          ! Loop over angular momentum channels
          do lup=0,paw_sp(isp)%lmax
             do mup=-lup,lup
                rgkl = realgaunt(lup,mup,lk,mk,ll,ml)
                if (abs(rgkl)<gaunt_tol) cycle
                iLM = lup*lup + lup + 1 + mup
                do is=1,nspins
                   ae_den(1:npts_phi,is,iLM) = ae_den(1:npts_phi,is,iLM) + &
                        rhoij_at(kpwtot,lpwtot,is)*den_work(1:npts_phi)*rgkl
                end do
             end do
          end do
       end do
    end do

    ! Add tcore density to tn1_{LM=00}(r)
    do is=1,nspins
       ae_den(1:npts_den,is,1) = ae_den(1:npts_den,is,1) + &
            paw_sp(isp)%core_den_rad(1:npts_den) &
            * sqrt_4pi / real(nspins,kind=DP)
    end do

  end subroutine paw_ae_density_rad_LM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_ps_density_rad_LM(ps_den,isp,rhoij_at,npw,npts_max,nspins, &
       nLM_max,den_work)

    !==================================================================!
    ! This subroutine calculates the soft pseudo sphere density in     !
    ! the PAW method, as a sum over the moments LM.                    !
    !  ps_den(:,is,iLM) = \tilde{n^1} density + psc density \tilde{n_c}!
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/05/10.                            !
    ! Split off to a separate routine by Nicholas Hine on 09/05/11.    !
    !==================================================================!

    use gaunt_coeff, only: realgaunt
    use services, only: services_analytic_limit

    implicit none

    ! Arguments
    integer, intent(in) :: npts_max
    integer, intent(in) :: nspins
    integer, intent(in) :: npw
    integer, intent(in) :: nLM_max
    real(kind=DP), intent(out) :: ps_den(npts_max,nspins,nLM_max)
    real(kind=DP), intent(out) :: den_work(npts_max)
    real(kind=DP), intent(in) :: rhoij_at(npw,npw,nspins)
    integer, intent(in) :: isp

    ! Local Variables
    integer :: iLM
    integer :: npts_den
    integer :: npts_phi
    integer :: kpwtot,lpwtot
    integer :: kpw,lpw
    integer :: lk,mk,ll,ml
    integer :: lup,mup
    integer :: is
    integer :: igrid
    real(kind=DP) :: rgkl
    real(kind=DP),parameter :: gaunt_tol = 1.0e-14_DP
    real(kind=DP),parameter :: rho_tol = 1.0e-14_DP

    ! Initialise, and find information about this species and its phi_grid
    ps_den(:,:,:) = 0.0_DP
    igrid = paw_sp(isp)%phi_grid
    npts_phi = paw_sp(isp)%grid(igrid)%npt
    igrid = paw_sp(isp)%core_den_grid
    npts_den = paw_sp(isp)%grid(igrid)%npt

    ! Double loop over partial waves on this atom
    do kpwtot=1,paw_sp(isp)%npw_tot
       kpw = paw_sp(isp)%ipw_tot(kpwtot)
       lk = paw_sp(isp)%l_pw_tot(kpwtot)
       mk = paw_sp(isp)%m_pw_tot(kpwtot)
       do lpwtot=1,paw_sp(isp)%npw_tot
          lpw = paw_sp(isp)%ipw_tot(lpwtot)
          ll = paw_sp(isp)%l_pw_tot(lpwtot)
          ml = paw_sp(isp)%m_pw_tot(lpwtot)

          if (all(abs(rhoij_at(kpwtot,lpwtot,:))<rho_tol)) cycle

          ! den_work(:) = tphi_i(r)*tphi_j(r) / r^2
          den_work(2:npts_phi) = paw_sp(isp)%tphi_rad(2:npts_phi,kpw) &
               * paw_sp(isp)%tphi_rad(2:npts_phi,lpw) &
               / paw_sp(isp)%grid(igrid)%r(2:npts_phi)**2
          ! Find analytic limit
          call services_analytic_limit(npts_phi,paw_sp(isp)%grid(igrid)%r(:), &
               den_work(:))

          ! Loop over angular momentum channels
          do lup=0,paw_sp(isp)%lmax
             do mup=-lup,lup
                rgkl = realgaunt(lup,mup,lk,mk,ll,ml)
                if (abs(rgkl)<gaunt_tol) cycle
                iLM = lup*lup + lup + 1 + mup
                do is=1,nspins
                   ps_den(1:npts_phi,is,iLM) = ps_den(1:npts_phi,is,iLM) + &
                        rhoij_at(kpwtot,lpwtot,is)*den_work(1:npts_phi)*rgkl
                end do
             end do
          end do
       end do
    end do

    ! Add tcore density to tn1_{LM=00}(r)
    do is=1,nspins
       ps_den(1:npts_den,is,1) = ps_den(1:npts_den,is,1) + &
            paw_sp(isp)%tcore_den_rad(1:npts_den) &
            * sqrt_4pi / real(nspins,kind=DP)
    end do

  end subroutine paw_ps_density_rad_LM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine paw_aug_density_rad_LM(aug_den,isp,rhoij_at,npw,npts_max, &
       nspins,nLM_max,den_work)

    !==================================================================!
    ! This subroutine calculates the various sphere-only densities in  !
    ! the PAW as a sums over their moments LM.                         !
    !  den(:,is,iLM,1) = n^1 density + core density n_c                !
    !  den(:,is,iLM,2) = \tilde{n^1} density + psc density \tilde{n_c} !
    !  den(:,is,iLM,3) = \hat{n} compensation density                  !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/05/10.                            !
    !==================================================================!

    use gaunt_coeff, only: realgaunt

    implicit none

    ! Arguments
    integer, intent(in) :: npts_max
    integer, intent(in) :: nspins
    integer, intent(in) :: npw
    integer, intent(in) :: nLM_max
    real(kind=DP), intent(out) :: aug_den(npts_max,nspins,nLM_max)
    real(kind=DP), intent(out) :: den_work(npts_max)
    real(kind=DP), intent(in) :: rhoij_at(npw,npw,nspins)
    integer, intent(in) :: isp

    ! Local Variables
    integer :: iLM
    integer :: npts_shp
    integer :: kpwtot,lpwtot
    integer :: kpw,lpw
    integer :: lk,mk,ll,ml
    integer :: lup,mup
    integer :: is
    integer :: igrid
    real(kind=DP) :: rgkl
    real(kind=DP),parameter :: gaunt_tol = 1.0e-14_DP
    real(kind=DP),parameter :: rho_tol = 1.0e-14_DP

    ! Initialise, and find information about this species and its phi_grid
    aug_den(:,:,:) = 0.0_DP
    igrid = paw_sp(isp)%shape_grid
    npts_shp = paw_sp(isp)%grid(igrid)%npt

    ! Double loop over partial waves on this atom
    do kpwtot=1,paw_sp(isp)%npw_tot
       kpw = paw_sp(isp)%ipw_tot(kpwtot)
       lk = paw_sp(isp)%l_pw_tot(kpwtot)
       mk = paw_sp(isp)%m_pw_tot(kpwtot)
       do lpwtot=1,paw_sp(isp)%npw_tot
          lpw = paw_sp(isp)%ipw_tot(lpwtot)
          ll = paw_sp(isp)%l_pw_tot(lpwtot)
          ml = paw_sp(isp)%m_pw_tot(lpwtot)

          if (all(abs(rhoij_at(kpwtot,lpwtot,:))<rho_tol)) cycle

          ! Loop over angular momentum channels
          do lup=0,paw_sp(isp)%lmax
             ! work(3) = g_L(r)*q_ij
             den_work(1:npts_shp) = paw_sp(isp)%aug_nLij(kpw,lpw,lup) * &
                  paw_sp(isp)%shape%shape_rad(1:npts_shp,lup)
             do mup=-lup,lup
                rgkl = realgaunt(lup,mup,lk,mk,ll,ml)
                if (abs(rgkl)<gaunt_tol) cycle
                iLM = lup*lup + lup + 1 + mup
                do is=1,nspins
                   aug_den(1:npts_shp,is,iLM) = aug_den(1:npts_shp,is,iLM) + &
                        rhoij_at(kpwtot,lpwtot,is)*den_work(1:npts_shp)*rgkl
                end do
             end do
          end do
       end do
    end do

  end subroutine paw_aug_density_rad_LM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_tcore_hartree_calc_forces(den_slabs12,grid,elements,loc_forces)

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the Hartree potential of the pseudized core charge.                !
    !  F_I = d/dR_I (\int (\tilde{n} + \hat{n})(r) v_H[\tilde{n}_{Zc}](r)) dr !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  den_slabs12  : input  : ground state data-parallelised charge density  !
    !  elements     : input  : list of elements and corresponding info        !
    !  loc_forces   : output : ionic forces due to local part                 !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 09/06/2010.                                 !
    !=========================================================================!

    use cell_grid, only: cell_grid_recip_pt, GRID_INFO
    use comms, only: comms_reduce, pub_my_node_id, pub_on_root
    use constants, only: PI
    use fourier, only: fourier_apply_cell_forward
    use geometry, only: POINT, OPERATOR(.dot.)
    use ion, only: ELEMENT
    use rundat, only: pub_smooth_loc_pspot
    use services, only: services_1d_interpolation
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=dp), intent(inout) :: den_slabs12(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_cell%num_spins)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    real(kind=dp), intent(out) :: loc_forces(1:3,pub_cell%nat)

    ! Local Variables
    integer :: i2,i3,islab23
    integer :: species,atom
    integer :: ierr
    real(kind=DP) :: gvec(3)
    real(kind=DP) :: g_length,inv_gsq,gdotR
    real(kind=DP) :: factor,v_loc_value
    complex(kind=DP) :: eiGR
    complex(kind=DP), allocatable :: iG_vden(:,:)
    complex(kind=DP), allocatable :: recip(:,:,:)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering paw_tcore_hartree_calc_forces'
#endif

    ! Start timer
    call timer_clock('paw_tcore_hartree_calc_forces',1)

    ! Allocate
    allocate(iG_vden(1:3,pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('paw_tcore_hartree_calc_forces','iG_vden',ierr)
    allocate(recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('paw_tcore_hartree_calc_forces','recip',ierr)

    ! Initialise
    loc_forces = 0.0_DP

    ! If spin polarised, put total density in up spin
    if (pub_cell%num_spins == 2) den_slabs12(:,:,:,1) = &
         den_slabs12(:,:,:,1) + den_slabs12(:,:,:,2)

    ! Fourier transform the charge density to reciprocal space
    call fourier_apply_cell_forward(den_slabs12(:,:,:,1),recip,grid)

    ! For components g with symmetry points at -g
    factor = 2.0_DP

    ! For g1=0 slab, which has no symmetry points
    if (grid%first_slab23(pub_my_node_id)==1) factor=1.0_DP

    ! Loop along b1
    do islab23=1,grid%num_slabs23

       ! Loop along b2
       do i2=1,grid%n2

          ! Loop along b3
          do i3=1,grid%n3

             ! g-vector, |g| and 1/|g|^2
             call cell_grid_recip_pt(gvec,islab23 + &
                  grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)
             g_length = sqrt(sum(gvec(1:3)**2))
             inv_gsq = grid%coulomb_recip(i3,i2,islab23)

             ! ndmh: loop over species, calculating v(G) for each.
             ! ndmh: this bit is O(N) so as much as possible should be
             ! ndmh: pre-calculated here.
             do species=1,pub_cell%num_pspecies

                ! Interpolate value of local potential at current g
                v_loc_value = services_1d_interpolation( &
                     paw_sp(species)%vhntzc_recip, &
                     paw_sp(species)%n_recip_pts,g_length * &
                     paw_sp(species)%inv_g_spacing,0)

                ! Add back the Coulomb potential; set g=0 term to zero
                if (g_length .gt. 0.0_DP) then
                   ! pa: changed to allow fractional ionic charge
                   v_loc_value = v_loc_value - &
                        paw_sp(species)%ion_charge*inv_gsq
                else
                   v_loc_value = 0.0_DP
                endif

                ! ndmh: calculate iG.v(G).n*(G) for each species and
                ! ndmh: each Cartesian direction
                iG_vden(1:3,species) = factor * cmplx(0.0_DP,1.0_DP) &
                     * gvec(1:3) * v_loc_value * conjg(recip(i3,i2,islab23))

             end do

             ! ndmh: loop over atoms: this bit is O(N^2) so should be made as
             ! ndmh: short and simple as possible
             do atom=1,pub_cell%nat

                ! e^{-ig.R}
                gdotR =  -(gvec(1)*elements(atom)%centre%x + &
                     gvec(2)*elements(atom)%centre%y + &
                     gvec(3)*elements(atom)%centre%z)
                eiGR = cmplx(COS(gdotR),SIN(gdotR),kind=DP)

                ! Sum force over g in each Cartesian direction i
                ! ==>  f = sum_{g} i.g.e^{-ig.R}.V_H[n_Zc](g).den^{*}(g)
                loc_forces(:,atom) = loc_forces(:,atom) + &
                     real(iG_vden(:,elements(atom)%pspecies_number)*eiGR,kind=DP)

             enddo     ! End loop over atoms

          enddo     ! End loop along b3

       enddo     ! End loop along b2

       ! For g1/=0 slabs
       factor=2.0_DP

    enddo     ! End loop along b1

    ! Sum the result over all nodes
    call comms_reduce('SUM',loc_forces)

    ! Deallocate
    deallocate(recip,stat=ierr)
    call utils_dealloc_check('paw_tcore_hartree_calc_forces','recip',ierr)
    deallocate(iG_vden,stat=ierr)
    call utils_dealloc_check('paw_tcore_hartree_calc_forces','iG_vden',ierr)

    ! Scale
    loc_forces = loc_forces * 4.0_dp * PI / (grid%n1*grid%n2*grid%n3)

    ! If spin polarised, restore up spin density
    if (pub_cell%num_spins == 2) den_slabs12(:,:,:,1) = den_slabs12(:,:,:,1) - &
         den_slabs12(:,:,:,2)

    ! Stop timer
    call timer_clock('paw_tcore_hartree_calc_forces',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving paw_tcore_hartree_calc_forces'
#endif

    return
  end subroutine paw_tcore_hartree_calc_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_nlcc_calculate_forces(density,tcore_density,nhat_den_grad, &
       grid,elements,nlcc_forces)

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the nonlinear core correction core charge in PAW.                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  density : input : ground state charge density on grid                  !
    !  tcore_density   : input  : NLCC core charge on grid                    !
    !  elements        : input  : list of elements and corresponding info     !
    !  grid            : input  : GRID_INFO type describing the grid          !
    !  nlcc_forces     : output : ionic forces due to NLCC corrections        !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine in June 2010.                                  !
    !=========================================================================!

    use cell_grid, only: cell_grid_recip_pt, GRID_INFO
    use comms, only: comms_reduce, pub_my_node_id, pub_on_root
    use fourier, only: fourier_apply_cell_forward
    use geometry, only: POINT, OPERATOR(.dot.)
    use ion, only: ELEMENT
    use rundat, only: pub_aug_den_dim, pub_nhat_in_xc
    use services, only: services_1d_interpolation
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use xc, only: xc_energy_potential

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=dp), intent(inout) :: density(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_cell%num_spins)
    real(kind=DP), intent(in) :: tcore_density(grid%ld1, &
         grid%ld2,grid%max_slabs12)
    real(kind=DP), intent(in) :: nhat_den_grad(grid%ld1, &
         grid%ld2,grid%max_slabs12,pub_cell%num_spins,0:pub_aug_den_dim)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    real(kind=dp), intent(out) :: nlcc_forces(1:3,pub_cell%nat)

    ! Local Variables
    integer :: i2,i3,islab23,is
    integer :: species,atom
    integer :: ierr
    real(kind=DP) :: gvec(3)
    real(kind=DP) :: g_length
    real(kind=DP) :: factor
    real(kind=DP) :: xc_energy
    real(kind=DP) :: coreden, GdotR
    real(kind=DP), allocatable :: total_density(:,:,:,:)
    real(kind=DP), allocatable :: xc_pot(:,:,:,:)
    complex(kind=DP), allocatable :: iG_coreden_vxc(:,:)
    complex(kind=DP) :: eiGR
    complex(kind=DP), allocatable :: recip(:,:,:)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering paw_nlcc_calculate_forces'
#endif

    ! Start timer
    call timer_clock('paw_nlcc_calculate_forces',1)

    ! Allocate
    allocate(total_density(grid%ld1,grid%ld2,grid%max_slabs12, &
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces','total_density',ierr)
    allocate(xc_pot(grid%ld1,grid%ld2,grid%max_slabs12,pub_cell%num_spins), &
         stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces','xc_pot',ierr)
    allocate(recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces','recip',ierr)
    allocate(iG_coreden_vxc(1:3,pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('paw_nlcc_calculate_forces','iG_coreden_vxc',ierr)

    ! Calculate total density
    factor = 1.0_DP / pub_cell%num_spins
    do is=1,pub_cell%num_spins
       total_density(:,:,:,is) = density(:,:,:,is) + tcore_density*factor
    end do

    if (.not.pub_nhat_in_xc) then
       total_density = total_density - nhat_den_grad(:,:,:,:,0)
    end if

    ! Calculate the exchange correlation potential
    call xc_energy_potential(total_density,xc_energy,xc_pot,grid, &
         pub_aug_den_dim,nhat_den_grad)

    ! Average the spin channels in (:,:,:,1)
    if (pub_cell%num_spins == 2) then
       xc_pot(:,:,:,1) = factor*(xc_pot(:,:,:,1) &
            + xc_pot(:,:,:,2))
    end if

    ! Initialise
    nlcc_forces = 0.0_DP

    ! Fourier transform the xc potential to reciprocal space
    call fourier_apply_cell_forward(xc_pot(:,:,:,1),recip,grid)

    ! For components g with symmetry points at -g
    factor=2.0_DP
    ! For g1=0 slab which has no symmetry points
    if (pub_my_node_id==grid%node_slab23(1)) factor=1.0_DP

    ! Loop along b1
    do islab23=1,grid%num_slabs23

       ! Loop along b2
       do i2=1,grid%n2

          ! Loop along b3
          do i3=1,grid%n3

             ! Get g-vector and |g|
             call cell_grid_recip_pt(gvec,islab23 + &
                  grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)
             g_length = sqrt(sum(gvec(1:3)**2))

             ! Loop over atoms to find n(G) for each
             do species=1,pub_cell%num_pspecies

                ! Check if we actually have a core charge for this species
                if (.not.paw_sp(species)%tcore_charge) cycle

                ! Interpolate value of core density at current g
                coreden = services_1d_interpolation( &
                     paw_sp(species)%tcore_den_recip, &
                     paw_sp(species)%n_recip_pts, &
                     g_length*paw_sp(species)%inv_g_spacing,0)

                ! ndmh: calculate iG.n(G).vxc*(G) for each species and
                ! ndmh: each Cartesian direction
                iG_coreden_vxc(1:3,species) = factor * cmplx(0.0_DP,1.0_DP) &
                     * gvec(1:3) * coreden * conjg(recip(i3,i2,islab23))

             end do

             ! ndmh: loop over atoms: this bit is O(N^2) so should be made as
             ! ndmh: short and simple as possible
             do atom=1,pub_cell%nat

                species = elements(atom)%pspecies_number

                ! Check if we actually have a core charge for this atom
                if (.not.paw_sp(species)%tcore_charge) cycle

                ! e^{-ig.R}
                gdotR = -(gvec(1)*elements(atom)%centre%x + &
                     gvec(2)*elements(atom)%centre%y + &
                     gvec(3)*elements(atom)%centre%z)
                eiGR = cmplx(COS(gdotR),SIN(gdotR),kind=DP)

                ! Sum force over g in each Cartesian direction i
                ! ==>  f = sum_{g} i.g.e^{-ig.R}.den_core(g).Vxc^{*}(g)
                nlcc_forces(:,atom) = nlcc_forces(:,atom) + &
                     real(iG_coreden_vxc(:,species)*eiGR,kind=DP)

             enddo     ! End loop over atoms

          enddo     ! End loop along b3

       enddo     ! End loop along b2

       ! For g1/=0 slabs
       factor=2.0_DP

    enddo     ! End loop along b1

    ! Sum the result over all nodes
    call comms_reduce('SUM',nlcc_forces)

    ! Deallocate
    deallocate(iG_coreden_vxc,stat=ierr)
    call utils_dealloc_check('paw_nlcc_calculate_forces','iG_coreden_vxc',ierr)
    deallocate(recip,stat=ierr)
    call utils_dealloc_check('paw_nlcc_calculate_forces','recip',ierr)
    deallocate(xc_pot,stat=ierr)
    call utils_dealloc_check('paw_nlcc_calculate_forces','xc_pot',ierr)
    deallocate(total_density,stat=ierr)
    call utils_dealloc_check('paw_nlcc_calculate_forces', &
         'total_density',ierr)

    ! Scale factor
    nlcc_forces = nlcc_forces / real(grid%n1*grid%n2*grid%n3,dp)

    ! Stop timer
    call timer_clock('paw_nlcc_calculate_forces',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving paw_nlcc_calculate_forces'
#endif

    return
  end subroutine paw_nlcc_calculate_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_sphere_density_on_grid(density_on_grid,grid,rhoij, &
       ae_coeff,ps_coeff)

    !=========================================================================!
    ! This subroutine calculates the PAW sphere terms of the density on a     !
    ! regular real space grid (not used in the main calculation, but useful   !
    ! for analysing results afterwards).                                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  density_on_grid (inout) : on input, the soft PS density (or zero)      !
    !                            on output, the contents on input plus the    !
    !                            sphere ae density times ae_coeff plus the    !
    !                            sphere ps density times ps_coeff             !
    !  grid (input) : GRID_INFO type describing grid                          !
    !  rhoij (input) : projector density kernel                               !
    !  ae_coeff (input) : coefficient of ae density                           !
    !  ps_coeff (input) : coefficient of ps density                           !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 09/06/2011                                    !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_deposit_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: pub_on_root, pub_my_node_id
    use constants, only: VERBOSE
    use fourier, only: fourier_apply_box
    use parallel_strategy, only: pub_elements_on_node, pub_max_atoms_on_node, &
         pub_first_atom_on_node, pub_num_atoms_on_node
    use rundat, only: pub_paw, pub_usp, pub_output_detail, pub_aug_den_dim, &
         pub_aug_funcs_recip
    use services, only: services_radial_transform, services_radial_integral_rmax
    use simulation_cell, only: pub_cell, pub_aug_box_n1, pub_aug_box_n2, &
         pub_aug_box_n3
    use sparse, only: SPAM3, sparse_get_block
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: density_on_grid(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_cell%num_spins)
    type(SPAM3), intent(in) :: rhoij(pub_cell%num_spins)
    real(kind=DP), intent(in) :: ae_coeff, ps_coeff

    ! Local Variables
    real(kind=DP),allocatable :: density_box(:,:,:,:)
    real(kind=DP),allocatable :: buffer(:,:,:)
    real(kind=DP),allocatable :: rho_ij_block(:,:,:)
    real(kind=DP),allocatable :: ae_den_rad_LM(:,:,:)
    real(kind=DP),allocatable :: ps_den_rad_LM(:,:,:)
    real(kind=DP),allocatable :: den_work(:)
    integer :: iat, loc_iat
    integer :: ierr
    integer :: is
    integer :: isp
    integer :: npw, nspins
    integer :: box_n1,box_n2,box_n3
    integer :: box_start1,box_start2,box_start3
    integer :: nLM_max, npts_max, npts_recip, npts_phi
    integer :: igrid
    integer :: lup, mup, iLM
    logical :: i_have_box

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering paw_sphere_density_on_grid'
#endif

    ! Start timer
    call timer_clock('paw_sphere_density_on_grid',1)

    ! Find sizes of box, shorthands, grid sizes tec
    box_n1 = pub_aug_box_n1
    box_n2 = pub_aug_box_n2
    box_n3 = pub_aug_box_n3
    npw = maxval(paw_projectors%species_num_proj(:))
    nspins = pub_cell%num_spins
    npts_max = 0
    nLM_max = 0
    npts_recip = 0
    do isp=1,pub_cell%num_pspecies
       if (.not.any(pub_elements_on_node(:)%pspecies_number==isp)) cycle
       npts_max = max(npts_max,paw_sp(isp)%grid(paw_sp(isp)%phi_grid)%npt)
       npts_max = max(npts_max,paw_sp(isp)%grid(paw_sp(isp)%core_den_grid)%npt)
       nLM_max = max(nLM_max,(paw_sp(isp)%lmax+1)**2)
       npts_recip = max(npts_recip,paw_sp(isp)%n_recip_pts)
    end do

    ! Allocate workspace
    allocate(rho_ij_block(npw,npw,nspins),stat=ierr)
    call utils_alloc_check('paw_sphere_density_on_grid','rho_ij_block',ierr)
    allocate(density_box(box_n1,box_n2,box_n3,nspins),stat=ierr)
    call utils_alloc_check('paw_sphere_density_on_grid','density_box',ierr)
    allocate(buffer(box_n1,box_n2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('paw_sphere_density_on_grid','buffer',ierr)
    allocate(ae_den_rad_LM(npts_max,nspins,nLM_max),stat=ierr)
    call utils_alloc_check('paw_sphere_density_on_grid','ae_den_rad_LM',ierr)
    allocate(ps_den_rad_LM(npts_max,nspins,nLM_max),stat=ierr)
    call utils_alloc_check('paw_sphere_density_on_grid','ps_den_rad_LM',ierr)
    allocate(den_work(npts_max),stat=ierr)
    call utils_alloc_check('paw_sphere_density_on_grid','den_work',ierr)

    do loc_iat=1,pub_max_atoms_on_node

       if (loc_iat<=pub_num_atoms_on_node(pub_my_node_id)) then
          iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1
          isp = pub_elements_on_node(loc_iat)%pspecies_number
          igrid = paw_sp(isp)%phi_grid
          npts_phi = paw_sp(isp)%grid(igrid)%npt

          ! Find where box for this atom is located in simulation cell
          call cell_grid_box_start_wrt_atom( &
               box_start1, box_start2, box_start3, &
               pub_elements_on_node(loc_iat)%centre, box_n1, box_n2, box_n3, &
               grid)

          ! Get block of \rho_ij for this atom
          do is=1,nspins
             call sparse_get_block(rho_ij_block(:,:,is),rhoij(is),iat,iat)
          end do

          ! Generate the ae and ps densities for this atom in the augmentation
          ! box using rhoij
          call paw_ae_density_rad_LM(ae_den_rad_LM,isp,rho_ij_block,npw, &
               npts_max,nspins,nLM_max,den_work)
          call paw_ps_density_rad_LM(ps_den_rad_LM,isp,rho_ij_block,npw, &
               npts_max,nspins,nLM_max,den_work)

          density_box(:,:,:,:) = 0.0_DP
          do lup=0,paw_sp(isp)%lmax
             do mup=-lup,lup
                iLM = lup*lup + lup + 1 + mup

                do is=1,nspins

                   ! Create density to be created, from ae and ps coefficients
                   den_work(:) = ae_coeff * ae_den_rad_LM(:,is,iLM) + &
                        ps_coeff * ps_den_rad_LM(:,is,iLM)

                   call paw_den_to_box_real(density_box(:,:,:,is),den_work, &
                        paw_sp(isp)%grid(igrid)%r,paw_sp(isp)%rcut,npts_phi,lup,mup, &
                        box_start1,box_start2,box_start3,grid,box_n1, &
                        box_n2,box_n3,pub_elements_on_node(loc_iat)%centre)
                end do  ! is
             end do  ! mup
          end do  ! lup

          i_have_box = .true.
       else
          ! Nothing to deposit on this node
          i_have_box = .false.
       end if

       ! Deposit this box to the simulation cell if present, or just wait for
       ! data from other nodes if no box
       do is=1,pub_cell%num_spins
          call cell_grid_deposit_box(density_on_grid(:,:,:,is), &
               density_box(:,:,:,is), buffer, grid, &
               box_n1, box_n2, box_n3, box_n1, box_n2, &
               box_start1, box_start2, box_start3, i_have_box, .false.)
       end do
    end do

    ! Deallocate
    deallocate(den_work,stat=ierr)
    call utils_dealloc_check('paw_sphere_density_on_grid','den_work',ierr)
    deallocate(ps_den_rad_LM,stat=ierr)
    call utils_dealloc_check('paw_sphere_density_on_grid','ps_den_rad_LM',ierr)
    deallocate(ae_den_rad_LM,stat=ierr)
    call utils_dealloc_check('paw_sphere_density_on_grid','ae_den_rad_LM',ierr)
    deallocate(density_box,stat=ierr)
    call utils_dealloc_check('paw_sphere_density_on_grid','density_box',ierr)
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('paw_sphere_density_on_grid','buffer',ierr)
    deallocate(rho_ij_block,stat=ierr)
    call utils_dealloc_check('paw_sphere_density_on_grid','rho_ij_block',ierr)

    ! Stop timer
    call timer_clock('paw_sphere_density_on_grid',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving paw_sphere_density_on_grid'
#endif

    return
  end subroutine paw_sphere_density_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine paw_den_to_box_real(box,den_lm,r,rcut,npts,lup,mup, &
      cell_start1,cell_start2,cell_start3,grid,box_n1,box_n2,box_n3, &
      atom_origin)

    !==================================================================!
    ! This subroutine calculates a l,m component of the density in     !
    ! real space within a sphere of radius rcut, provided on a grid    !
    ! of points r(1:npts), where the density is given by               !
    !   n_LM(\vec r) = n_LM(r) S_LM(\vec r)                            !
    ! It does so on the simulation cell fine grid in a box centered on !
    ! the atom.                                                        !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  box_(out) :                                                     !
    !  den_lm (in) :                                                   !
    !  r (in) :                                                        !
    !  rcut (in) :                                                     !
    !  npts (in) :                                                     !
    !  lup (in) :                                                      !
    !  mup (in) :                                                      !
    !  cell_start1 (in) :                                              !
    !  cell_start2 (in) :                                              !
    !  cell_start3 (in) :                                              !
    !  grid (in) :                                                     !
    !  box_n1 (in) :                                                   !
    !  box_n2 (in) :                                                   !
    !  box_n3 (in) :                                                   !
    !  atom_origin (in) :                                              !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 09/05/11.           !
    !==================================================================!

    use basis, only: basis_box_origin_to_atom
    use cell_grid, only: GRID_INFO
    use constants, only: PI, cmplx_i
    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use services, only: services_locate_interp,services_linear_interpolation
    use simulation_cell, only: pub_cell
    use spherical_wave, only: sw_real_sph_harm, sw_init

    implicit none

    ! Arguments
    integer,intent(in) :: box_n1
    integer,intent(in) :: box_n2
    integer,intent(in) :: box_n3
    real(kind=DP),intent(inout) :: box(box_n1,box_n2,box_n3)
    integer,intent(in) :: npts
    real(kind=DP), intent(in) :: den_lm(npts)
    real(kind=DP), intent(in) :: r(npts), rcut
    integer,intent(in) :: lup
    integer,intent(in) :: mup
    integer,intent(in) :: cell_start1
    integer,intent(in) :: cell_start2
    integer,intent(in) :: cell_start3
    type(GRID_INFO),intent(in) :: grid
    type(POINT),intent(in) :: atom_origin

    ! Local Variables
    integer :: i1,i2,i3
    integer :: ipt
    type(POINT) :: box_origin
    type(POINT) :: box_to_atom
    real(kind=DP) :: denval
    real(kind=DP) :: slmval
    type(POINT) :: r_cell, r_sphere
    real(kind=DP) :: rmag

    if (lup>4) call sw_init(lup+1,1)

    call basis_box_origin_to_atom(box_to_atom,box_origin,atom_origin, &
         cell_start1,cell_start2,cell_start3,grid%da1,grid%da2,grid%da3)

    do i3=1,box_n3
       do i2=1,box_n2
          do i1=1,box_n1

             r_cell = box_origin + (i1-1)*grid%da1 &
                                 + (i2-1)*grid%da2 &
                                 + (i3-1)*grid%da3

             r_sphere = r_cell - atom_origin
             rmag = geometry_magnitude(r_sphere)

             if (rmag>rcut) then
                box(i1,i2,i3) = 0.0_DP
                cycle
             end if

             ipt = services_locate_interp(rmag,r,npts)
             denval = services_linear_interpolation(rmag,den_lm(ipt), &
                  den_lm(ipt+1),r(ipt),r(ipt+1))

             ! multiply value of shapefunc at each g-point by the appropriate
             ! phase factor (n.b. real part and imaginary part stored as
             ! separate consecutive elements in the x-direction of the array),
             ! and the appropriate real spherical harmonic factor.
             slmval = sw_real_sph_harm(r_sphere%x, &
                  r_sphere%y,r_sphere%z,rmag,lup,mup)

             box(i1,i2,i3) = box(i1,i2,i3) + &
                  denval * slmval

          end do
       end do
    end do

  end subroutine paw_den_to_box_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module paw
