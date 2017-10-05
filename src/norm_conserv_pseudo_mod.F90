! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris, Arash A. Mostofi and Nicholas D.M. Hine
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module pseudopotentials

  use constants, only : DP, PI, stdout
  use geometry, only : POINT
  use projectors, only : PROJECTOR_SET

  implicit none

  private

  ! cks, 22/1/2004: type for containing pseudopotential information
  type PSEUDO_SPECIES

     ! cks: name of pseudopotential file
     character(len=64) :: pseudo_name

     ! cks: atomic number of atom
     integer :: atomic_number

     ! cks: charge of "pseudised" core of atom
     real(kind=DP) :: ion_charge   ! pa: changed from integer to real

     ! ndmh: whether to subtract the coulomb potential -Z/r
     logical :: subtract_coul

     ! cks: number of points in the radial grid
     integer :: n_rad_pts

     ! ndmh: spacing of points in the radial grid
     real(kind=DP) :: inv_g_spacing

     ! cks: number of angular momentum shells of projectors
     integer :: n_shells

     ! cks: angular momentum of each shell, this
     ! cks: array is allocated to be ang_mom(n_shells)
     integer, allocatable, dimension(:) :: ang_mom

     ! cks: maximum radial reciprocal k-vector up to which
     ! cks: the potential and projectors are defined
     real(kind=DP) :: ps_gmax

     ! cks: core radius for each projector shell,
     ! cks: core_radius(n_shells)
     real(kind=DP), allocatable, dimension(:) :: core_radius

     ! ndmh: Kleinman-Bylander energies for each pair of projector
     ! ndmh: shells: D0(n_shells,n_shells)
     real(kind=DP), allocatable, dimension(:,:) :: D0

     ! cks: the radial recip part of the local potential
     ! cks: rad_locpot_recip(n_rad_pts)
     real(kind=DP), allocatable, dimension(:) :: rad_locpot_recip

     ! cks: the values of each projector (shell) on the radial
     ! cks: grid. The last column of this array contains the
     ! cks: grid points, so the size of the array is
     ! cks: rad_proj_recip(n_rad_pts, n_shells +1 )
     real(kind=DP), allocatable, dimension(:,:) :: rad_proj_recip

     ! ndmh: whether a core charge is present for this species
     logical :: core_charge

     ! ndmh: if core charge is present, array to store core charge in
     ! ndmh: reciprocal space
     real(kind=DP), allocatable, dimension(:) :: core_charge_recip

     ! ndmh: the augmentation charges for each shell
     real(kind=DP), allocatable, dimension(:,:) :: aug_q

     ! ndmh: the L-independent function
     integer :: kkbeta, mesh
     real(kind=DP), allocatable, dimension(:) :: rlog
     real(kind=DP), allocatable, dimension(:) :: rab
     real(kind=DP), allocatable, dimension(:,:,:) :: qfunc
     real(kind=DP), allocatable, dimension(:) :: rinner

     ! ndmh: the L-dependent coefficients
     integer :: nqfcoef
     integer :: qf_lmax
     real(kind=DP), allocatable, dimension(:,:,:,:) :: qfcoef

  end type PSEUDO_SPECIES

  ! cks: this array will be num_pspecies long and will be allocated in
  ! cks: this module
  type(PSEUDO_SPECIES), allocatable, dimension(:) :: p_species

  ! Projector type to describe nonlocal projectors
  type(PROJECTOR_SET), public :: nlps_projectors

  ! ndmh: subroutines relating to pseudopotential initialisation/exit
  public :: pseudopotentials_read_species
  public :: pseudo_species_init_proj
  public :: pseudopotentials_species_exit

  ! ndmh: routines for retrieving the pseudopotential for use by atom_mod
  public :: pseudo_get_locpot_rad
  public :: pseudo_get_core_den_rad
  public :: pseudo_get_projectors_q
  public :: pseudo_get_aug_funcs
  public :: pseudo_get_projector_info

  ! ndmh: subroutines relating to construction of local potential and forces
  public :: pseudo_make_structure_factor
  public :: pseudopotentials_local_on_grid
  public :: pseudopotentials_core_density
  public :: pseudo_local_calculate_forces
  public :: pseudo_nlcc_calculate_forces
  public :: pseudo_local_on_grid_openbc

  ! ndmh: subroutines relating to nonlocal potentials
  public :: pseudo_get_dij
  public :: pseudopotentials_nonlocal_mat
  public :: pseudo_nl_calculate_forces
  public :: pseudo_aug_Q_matrix
  public :: pseudo_atom_aug_den
  public :: pseudo_atom_aug_integrals
  public :: pseudo_atom_aug_force
  public :: pseudo_nonlocal_commutator_mat
  public :: pseudo_nonlocal_com_mat_fd
  ! lr408 - for debugging purposes only
  public :: pseudo_nonlocal_com_mat_direct           


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SUBROUTINES FOR PSEUDOPOTENTIAL INITIALISATION !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudopotentials_read_species(elements)

    !=================================================================!
    ! This subroutine allocates and initialises the p_species array   !
    ! whose elements are PSEUDO_SPECIES types - one for each distinct !
    ! ionic species. It also sets the pub_cell%num_projectors value.  !
    ! The norm-conserving pseudopotential for each ionic species is   !
    ! read from a file in CASTEP .recpot format.                      !
    !-----------------------------------------------------------------!
    ! WARNING: The fftbox_proj_recip corresponding to each p_species  !
    !          are not initialised by this subroutine. They are       !
    !          initialised by a separate call to the                  !
    !          pseudo_species_init_proj subroutine.                   !
    !-----------------------------------------------------------------!
    ! Written  by Chris-Kriton Skylaris on 23/1/2004.                 !
    ! Tidied up by Peter D. Haynes on 1/7/2004.                       !
    !=================================================================!

#ifdef ACCELRYS
    use constants, only: stderr
    use comms, only: pub_my_node_id
#endif
    use comms, only: comms_abort, comms_bcast, pub_on_root
    use ion, only: ELEMENT
    use rundat, only: pub_any_nl_proj, pub_nhat_in_xc
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(ELEMENT), intent(inout) :: elements(pub_cell%nat)

    ! Local Variables
    integer :: ierr                ! pdh: Error flag
    integer :: species_counter
    integer :: sp
    integer :: atom
    integer :: shell
    integer :: mom
    integer :: p_number
    logical :: we_have_it
    type(ELEMENT), allocatable :: element_buffer(:)

    ! Allocate workspace
    allocate(element_buffer(pub_cell%nat),stat=ierr)
    call utils_alloc_check('pseudopotentials_read_species', &
         'element_buffer',ierr)

    species_counter = 1
    element_buffer(1)%pseudo_name = elements(1)%pseudo_name
    element_buffer(1)%atomic_number = elements(1)%atomic_number

    ! cks: loop over all atoms as they appear in the input file
    ! cks: and fill element_buffer with information about each
    ! cks: distinct species.
    do atom=2,pub_cell%nat

       we_have_it = .false.

       ! cks: loop over all species found so far and check if we have this
       !      species
       species_loop: do sp=1,species_counter
          if (element_buffer(sp)%pseudo_name == elements(atom)%pseudo_name) then
             we_have_it = .true.
             exit species_loop
          end if
       end do species_loop

       ! cks: if we don't have it, let's get it!
       if (.not. we_have_it) then
          species_counter = species_counter + 1
          element_buffer(species_counter)%pseudo_name = elements(atom)%pseudo_name
          element_buffer(species_counter)%atomic_number = &
               elements(atom)%atomic_number
       end if

    end do

    ! cks: allocate p_species type(pseudo_species) array
    if (.not. allocated(p_species)) then
       allocate(p_species(species_counter),stat=ierr)
       call utils_alloc_check('pseudopotentials_read_species','p_species',ierr)
    else
       if (pub_on_root) write(stdout,'(a)') &
            'Error in pseudopotentials_read_species: &
            &p_species already allocated'
       call comms_abort
    end if

    ! cks: do some p_species initialisation from the element_buffer
    do sp=1,species_counter
       p_species(sp)%pseudo_name = element_buffer(sp)%pseudo_name
       p_species(sp)%atomic_number = element_buffer(sp)%atomic_number
    end do

    deallocate(element_buffer,stat=ierr)
    call utils_dealloc_check('pseudopotentials_read_species', &
         'element_buffer',ierr)

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! cks: Now connect the elements array to the species array
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    elements(:)%pspecies_number = -1
    do atom=1,pub_cell%nat

       ! cks: loop over all species found so far and check if we have this
       !      species
       number_set_loop: do sp=1,species_counter
          if (p_species(sp)%pseudo_name == elements(atom)%pseudo_name) then
             elements(atom)%pspecies_number = sp
             exit number_set_loop
          end if
       end do number_set_loop

       ! cks: sanity check
       if (elements(atom)%pspecies_number == -1) then
          if (pub_on_root) write(stdout,'(a,i6,a)') &
               'Error in pseudopotentials_read_species: elements(', &
               atom,')%pspecies_number = -1'
          call comms_abort
       end if

    end do
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    ! cks: end connect the elements array to the species array
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    ! cks: allocate memory for the allocatable components of p_species array
    call pseudo_species_alloc

    ! cks: complete the initialisation of the p_species array
    call pseudo_species_read_radial

    ! cks: find the total number of projectors
    pub_cell%num_projectors = 0
    pub_cell%num_pawpws = 0
    do atom=1,pub_cell%nat

#ifdef ACCELRYS
       ! ndmh: protection against dodgy compilers
       if(atom>size(elements))then
          write(stderr,"('ERROR in pseudopotentials_read_species error on &
               &node:',i4, ' atom ',i5,' does not fit element array(size:', &
               &i5,')')") pub_my_node_id,atom,size(elements)
          call comms_abort
       endif
#endif

       p_number = elements(atom)%pspecies_number
       do shell=1,p_species(p_number)%n_shells

          mom = p_species(p_number)%ang_mom(shell)
          pub_cell%num_projectors = pub_cell%num_projectors + 2*mom+1

       end do
    end do

    ! ndmh: set flag to denote presence of nl projectors
    ! ndmh 04/01/10: only if core radius of any species is > 0
    pub_any_nl_proj = .false.
    if (pub_cell%num_projectors > 0) then
       do sp=1,pub_cell%num_pspecies
          if (any(p_species(sp)%core_radius(:)>0.0_DP)) then
             pub_any_nl_proj = .true.
          end if
       end do
    end if
  !CW
    write(*,*) 'ARE THEY ANY NLPP ? : ', pub_any_nl_proj
  !END CW


    ! ndmh: set flag to include augmentation charge in XC calculation
    pub_nhat_in_xc = .true.

    ! cks: now set some redundant information in the elements
    ! cks: structure which we need for the time being but should
    ! cks: be removed in the future and be only available in elements
    do atom=1,pub_cell%nat

       elements(atom)%nprojectors     = 0
       elements(atom)%npawpws         = 0
       elements(atom)%ncorewfs        = 0
       elements(atom)%max_core_radius = 0.0_DP
       elements(atom)%max_core_wf_radius = 0.0_DP
       p_number = elements(atom)%pspecies_number

       elements(atom)%ion_charge = p_species(p_number)%ion_charge
       do shell=1,p_species(p_number)%n_shells

          mom = p_species(p_number)%ang_mom(shell)

          elements(atom)%nprojectors = elements(atom)%nprojectors + 2*mom+1

          elements(atom)%max_core_radius = &
               max(elements(atom)%max_core_radius, &
               p_species(p_number)%core_radius(shell))

       end do

    end do


  end subroutine pseudopotentials_read_species


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_species_alloc

    !============================================================!
    ! This subroutine reads the pseudopotentials for each atomic !
    ! species and allocates the necessary memory for them as     !
    ! elements of the p_species array.                           !
    !------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 23/1/2004.             !
    ! Tidied up by Peter Haynes on 1/7/2004.                     !
    ! Modified to read multiple file formats (recpots and usps)  !
    ! by Nicholas Hine on 29/01/2009.                            !
    !============================================================!

    use comms, only: comms_abort, comms_bcast, pub_on_root, pub_root_node_id
    use constants, only: ANGSTROM, HARTREE_IN_EVS
    use rundat, only: pub_nlcc, xc_functional, pub_aug, pub_usp
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check

    implicit none

    ! Local Variables
    integer, parameter  :: cfac=4
    integer             :: ierr
    integer             :: ps
    integer             :: tot_num_points
    integer             :: tot_num_projectors
    integer             :: num_projectors
    integer             :: i
    character(len=32)   :: string
    character(len=80)   :: line
    character(len=64)   :: current_file
    character(len=80)   :: current_job
    character(len=10)   :: psp_xc_functional
    logical             :: xc_mismatch
    logical             :: count_points
    logical             :: all_abort ! prevent hung procs if read of pseudo fails
    real(kind=DP)       :: gmax
    real(kind=DP)       :: temp
    real(kind=DP)       :: fact1
    real(kind=DP)       :: ionic_charge

    ! ndmh: determine if any species present contain core charges for NLCC
    pub_nlcc = .false.

    ! ndmh: Charge augmentation will not be active unless any species with
    ! ndmh: nonzero augmentation charges are found.
    pub_usp = .false.
    pub_aug = .false.

    ! cks: all processes will abort so none is left hung if read on root fails
    all_abort =.false.

    if (pub_on_root) then

       do ps=1,pub_cell%num_pspecies

          ! pdh: get file name for this species
          current_file = p_species(ps)%pseudo_name

          ! ndmh: determine what file format we are reading and scan it
          ! ndmh: to find sizes of various arrays

          ! ndmh: Old-style CASTEP .recpot
          if (index(current_file,'recpot') > 0) then

             call internal_scan_recpot

          ! ndmh: New-style CASTEP OTFG .usp (not necessarily ultrasoft!)
          else if (index(current_file,'usp') > 0) then

             call internal_scan_usp

          else
              write(stdout,"(a)")"Error in pseudopotentials_read_species: &
                   & Unrecognised file format: ", current_file
          end if

          ! ndmh: check if pseudopotential file contains a reference to a
          ! ndmh: different XC functional from the one we are using
          xc_mismatch = .false.
          select case (xc_functional)
          case ('LDA','CAPZ')
             if (psp_xc_functional=='GGA') xc_mismatch = .true.
          case ('GGA','PBE','PW91','RPBE','WC')
             if (psp_xc_functional=='LDA') xc_mismatch = .true.
          end select

          ! ndmh: report if there is a mismatch
          if (xc_mismatch.and.pub_on_root) then
             write(stdout,'(3a)')'WARNING in pseudopotentials_read_species: &
                  &string "',trim(psp_xc_functional),'" found in pseudopotential'
             write(stdout,'(5a)')'file "',trim(current_file), &
                  '", yet xc_functional = "',trim(xc_functional),'".'
             write(stdout,*)
          end if

          ! ndmh: check for NLCC core charge
          if (p_species(ps)%core_charge) pub_nlcc = .true.

       end do

    endif

    ! cks: abort all procs if some pseudo file was missing
    call comms_bcast(pub_root_node_id, all_abort)
    if (all_abort) call comms_abort


    ! ndmh: let all nodes know whether nonlinear core corrections are present
    call comms_bcast(pub_root_node_id,pub_nlcc)

    ! cks: let all nodes know what pub_root_node_id just got
    do ps=1,pub_cell%num_pspecies
       call comms_bcast(pub_root_node_id,p_species(ps)%n_shells)
       call comms_bcast(pub_root_node_id,p_species(ps)%n_rad_pts)
       call comms_bcast(pub_root_node_id,p_species(ps)%ion_charge)
       call comms_bcast(pub_root_node_id,p_species(ps)%ps_gmax)
       call comms_bcast(pub_root_node_id,p_species(ps)%inv_g_spacing)
       call comms_bcast(pub_root_node_id,p_species(ps)%core_charge)
       call comms_bcast(pub_root_node_id,p_species(ps)%subtract_coul)
       call comms_bcast(pub_root_node_id,p_species(ps)%kkbeta)
       call comms_bcast(pub_root_node_id,p_species(ps)%mesh)
       call comms_bcast(pub_root_node_id,p_species(ps)%nqfcoef)
       call comms_bcast(pub_root_node_id,p_species(ps)%qf_lmax)
    end do


    ! cks: now allocate memory for all allocatable components
    do ps=1,pub_cell%num_pspecies

       allocate(p_species(ps)%ang_mom(p_species(ps)%n_shells),stat=ierr)
       write(string,'(a,i2,a)') 'p_species(',ps,')%ang_mom'
       call utils_alloc_check('pseudo_species_alloc',trim(string),ierr)

       allocate(p_species(ps)%core_radius(p_species(ps)%n_shells),stat=ierr)
       write(string,'(a,i2,a)') 'p_species(',ps,')%core_radius'
       call utils_alloc_check('pseudo_species_alloc',trim(string),ierr)

       allocate(p_species(ps)%D0(p_species(ps)%n_shells, &
            p_species(ps)%n_shells),stat=ierr)
       write(string,'(a,i2,a)') 'p_species(',ps,')%D0'
       call utils_alloc_check('pseudo_species_alloc',trim(string),ierr)

       allocate(p_species(ps)%rad_locpot_recip(p_species(ps)%n_rad_pts), &
            stat=ierr)
       write(string,'(a,i2,a)') 'p_species(',ps,')%rad_locpot_recip'
       call utils_alloc_check('pseudo_species_alloc',trim(string),ierr)

       allocate(p_species(ps)%rad_proj_recip(p_species(ps)%n_rad_pts, &
            p_species(ps)%n_shells+1),stat=ierr)
       write(string,'(a,i2,a)') 'p_species(',ps,')%rad_proj_recip'
       call utils_alloc_check('pseudo_species_alloc',trim(string),ierr)

       ! ndmh: allocate memory for core charge if present
       if (p_species(ps)%core_charge) then
          allocate(p_species(ps)%core_charge_recip(p_species(ps)%n_rad_pts), &
               stat=ierr)
          write(string,'(a,i2,a)') 'p_species(',ps,')%core_charge_recip'
          call utils_alloc_check('pseudo_species_alloc',trim(string),ierr)
       end if

       ! ndmh: allocate memory for augmentation charges
       allocate(p_species(ps)%aug_q(p_species(ps)%n_shells, &
            p_species(ps)%n_shells), stat=ierr)
       write(string,'(a,i2,a)') 'p_species(',ps,')%aug_q'
       call utils_alloc_check('pseudo_species_alloc',trim(string),ierr)

       ! ndmh: allocate memory for L-independent function
       allocate(p_species(ps)%qfunc(p_species(ps)%kkbeta, &
            p_species(ps)%n_shells,p_species(ps)%n_shells), stat=ierr)
       write(string,'(a,i2,a)') 'p_species(',ps,')%qfunc'
       call utils_alloc_check('pseudo_species_alloc',trim(string),ierr)
       allocate(p_species(ps)%qfcoef(p_species(ps)%nqfcoef, &
            0:p_species(ps)%qf_lmax,p_species(ps)%n_shells, &
            p_species(ps)%n_shells),stat=ierr)
       write(string,'(a,i2,a)') 'p_species(',ps,')%qfcoef'
       call utils_alloc_check('pseudo_species_alloc',trim(string),ierr)
       allocate(p_species(ps)%rinner(0:p_species(ps)%qf_lmax),stat=ierr)
       write(string,'(a,i2,a)') 'p_species(',ps,')%rinner'
       call utils_alloc_check('pseudo_species_alloc',trim(string),ierr)
       allocate(p_species(ps)%rlog(p_species(ps)%mesh), stat=ierr)
       write(string,'(a,i2,a)') 'p_species(',ps,')%rlog'
       call utils_alloc_check('pseudo_species_alloc',trim(string),ierr)
       allocate(p_species(ps)%rab(p_species(ps)%mesh), stat=ierr)
       write(string,'(a,i2,a)') 'p_species(',ps,')%rab'
       call utils_alloc_check('pseudo_species_alloc',trim(string),ierr)
       p_species(ps)%qfcoef(:,:,:,:) = 0.0_DP
       p_species(ps)%qfunc(:,:,:) = 0.0_DP
       p_species(ps)%rlog(:) = 0.0_DP
       p_species(ps)%rab(:) = 0.0_DP

    end do

    return

  contains

    subroutine internal_scan_recpot

      !============================================================!
      ! This subroutine reads the pseudopotential for a .recpot    !
      ! and stores the sizes required to allocate the necessary    !
      ! memory for this element in the p_species array.            !
      !------------------------------------------------------------!
      ! Assembled by Nicholas Hine on 29/01/09 based on code by    !
      ! Chris-Kriton Skylaris from 23/1/2004 and tidied up by      !
      ! Peter Haynes on 1/7/2004.                                  !
      !============================================================!

      call pseudo_open_pspot_file(2,current_file,ierr)
      if (ierr /= 0) then
         write(stdout,'(3a,i6)') &
              'Error in pseudo_species_alloc: opening file "', &
              trim(current_file),'" failed with code ',ierr
!         call comms_abort
         all_abort =.true.
         return
      end if

      current_job = 'seeking marker "END COMMENT"'
      line = repeat(' ',80)
      psp_xc_functional = ''
      do while(index(line,'END COMMENT') == 0)
         read(2,'(a80)',err=100,end=200) line
         if (index(line,'LDA')>0) psp_xc_functional = 'LDA'
         if (index(line,'lda')>0) psp_xc_functional = 'LDA'
         if (index(line,'GGA')>0) psp_xc_functional = 'GGA'
         if (index(line,'gga')>0) psp_xc_functional = 'GGA'
         if (index(line,'PBE')>0) psp_xc_functional = 'GGA'
         if (index(line,'PW91')>0) psp_xc_functional = 'GGA'
      end do

      ! * Find out the number of projectors, and g-grid points
      !    --- stop at 1000 (or a flag of 4 characters) if the
      !        number of projectors is less than 32 (2 * s,p,d,f)

      count_points       = .true.
      tot_num_points     = 0
      tot_num_projectors = 0
      num_projectors     = 0

      current_job = 'counting projectors and points'
      line = repeat(' ',80)
      do while (len_trim(adjustl(line)) /= 4 .and. tot_num_projectors < 32)
         read(2,'(a80)',err=100,end=200) line

         if (len_trim(adjustl(line)) == 1) then
            read(line,*) i
            tot_num_projectors = tot_num_projectors + (2*i+1) ! 2l+1
            num_projectors = num_projectors + 1
            count_points = .false.
         end if

         if (count_points) then
            do i=1,len(line)
               if (line(i:i) == '.') tot_num_points = tot_num_points + 1
            end do
         end if
      end do

      ! cks: initialise
      p_species(ps)%n_shells  = num_projectors
      p_species(ps)%n_rad_pts = tot_num_points - 1

      ! ndmh: skip to the end of the projectors
      current_job = 'seeking marker "1000"'
      do while ((index(line,'1000') == 0).and.(len_trim(adjustl(line)) /= 4))
         read(2,'(a80)',err=100,end=200) line
      end do

      ! ndmh: check if there are numbers here - if so there is non-linear
      ! ndmh: core correction data to read in
      read(2,'(a80)',iostat=ierr) line
      if (ierr /= 0) then
         p_species(ps)%core_charge = .false.
      else
         p_species(ps)%core_charge = .true.
      end if

      ! ** Rewind and find out the ionic charge
      rewind(2,iostat=ierr,err=400)
400   if (ierr /= 0) then
         write(stdout,'(3a,i6)') &
              'Error in pseudo_species_alloc: rewinding file "', &
              trim(current_file),'" failed with code ',ierr
         call comms_abort
      end if

      current_job = 'seeking marker "END COMMENT"'
      line = repeat(' ',80)
      do while (index(line,'END COMMENT') == 0)
         read(2,'(a80)',err=100,end=200) line
      end do

      current_job = 'skipping version string'
      read(2,*,err=100,end=200) line
      current_job = 'reading Gmax'
      read(2,*,err=100,end=200) gmax
      current_job = 'reading first two points of local potential'
      read(2,*,err=100,end=200) temp,temp

      ! * Calculate the ionic charge by looking at the g -> 0 part
      !   of the local potential
      fact1 = 4.0_DP*PI*HARTREE_IN_EVS/ANGSTROM

      ! pa: changed from
      ! ionic_charge = real(nint(-temp*(gmax/real(tot_num_points-1,DP))**2/fact1),kind=DP)
      ! pa: to allow fractional ionic charges
      ionic_charge = -temp*(gmax/real(tot_num_points-1,DP))**2/fact1

      !write(stdout,'(a,f9.6)') &
      !     'ionic_charge from psp before rounding: ', ionic_charge

      ! pa: round ionic charge to the nearest 1/cfac
      ionic_charge = (real(nint(cfac*ionic_charge), kind=DP))/real(cfac,DP)

      !write(stdout,'(a,f9.6)') &
      !            'ionic_charge from psp after rounding is: ', ionic_charge

      ! cks: convert gmax to atomic units
      p_species(ps)%ps_gmax = gmax/ANGSTROM

      ! ndmh: calculate spacing of recip space points
      p_species(ps)%inv_g_spacing = (p_species(ps)%n_rad_pts-1) / &
           p_species(ps)%ps_gmax

      p_species(ps)%ion_charge = ionic_charge

      ! ndmh: recpots do include the -Z/r, so flag to subtract it
      p_species(ps)%subtract_coul = .true.

      ! ndmh: recpots do not use kkbeta, mesh, qf_lmax or nqfcoef
      p_species(ps)%kkbeta = 0
      p_species(ps)%mesh = 0
      p_species(ps)%qf_lmax = 0
      p_species(ps)%nqfcoef = 0

      ! ** Close the file
      close(2,iostat=ierr,err=500)
500   if (ierr /= 0) then
         write(stdout,'(3a,i6)') &
              'Error in pseudo_species_alloc: closing file "', &
              trim(current_file),'" failed with code ',ierr
         call comms_abort
      end if

      return

100   write(stdout,'(4a)') 'Error in pseudo_species_alloc: reading file "', &
           trim(current_file),'" failed while ',trim(current_job)
      call comms_abort

200   write(stdout,'(4a)') 'Error in pseudo_species_alloc: file "', &
           trim(current_file),'" ended unexpectedly while ',trim(current_job)
      call comms_abort

    end subroutine internal_scan_recpot

    subroutine internal_scan_usp

      !============================================================!
      ! This subroutine reads the pseudopotential for a .usp       !
      ! and stores the sizes required to allocate the necessary    !
      ! memory for this element in the p_species array.            !
      !------------------------------------------------------------!
      ! Written by Nicholas Hine on 29/01/09.                      !
      !============================================================!

      ! Local Variables
      integer, parameter :: max_proj=huge(1)
      integer :: version(2), usp_nlcc
      integer :: shell
      integer :: j,k
      integer :: lshell, lmax
      real(kind=DP) :: tmp1,tmp2

      call pseudo_open_pspot_file(2,current_file,ierr)
      if (ierr /= 0) then
         write(stdout,'(3a,i6)') &
              'Error in pseudo_species_alloc: opening file "', &
              trim(current_file),'" failed with code ',ierr
!         call comms_abort
          all_abort =.true.
          return
      end if

      current_job = 'seeking marker "END COMMENT"'
      line = repeat(' ',80)
      psp_xc_functional = ''
      do while(index(line,'END COMMENT') == 0)
         read(2,'(a80)',err=100,end=200) line
         if (index(line,'LDA')>0) psp_xc_functional = 'LDA'
         if (index(line,'lda')>0) psp_xc_functional = 'LDA'
         if (index(line,'GGA')>0) psp_xc_functional = 'GGA'
         if (index(line,'gga')>0) psp_xc_functional = 'GGA'
         if (index(line,'PBE')>0) psp_xc_functional = 'GGA'
         if (index(line,'PW91')>0) psp_xc_functional = 'GGA'
      end do

      ! ** Read the version
      read(2,*,err=100,end=200) version(1),version(2)

      ! ** Read the Ionic Charge
      read(2,*,err=100,end=200) ionic_charge

      p_species(ps)%ion_charge = ionic_charge

      ! ** Read the gmax, and nlcc flag
      read(2,*,err=100,end=200) gmax, usp_nlcc
      p_species(ps)%ps_gmax = gmax/ANGSTROM

      ! ** Determine whether to include nonlinear core corrections
      p_species(ps)%core_charge = (usp_nlcc==2)

      ! ** Find out the number of projectors, and g-grid points
      count_points = .true.
      tot_num_points = 0
      num_projectors = 0
      line = ''
      lmax = -100
      do while((len_trim(adjustl(line))/=4).and.(num_projectors<max_proj))
         read(2,'(a)',err=100,end=200) line
         if(len_trim(adjustl(line))==1) then
            read(line,*,err=100,end=200) lshell
            if (lshell>lmax) lmax = lshell
            num_projectors = num_projectors+1
            count_points = .false.
         end if
         if(count_points) then
            do i=1,len(line)
               if(line(i:i)=='.') tot_num_points = tot_num_points + 1
            end do
         end if
      end do

      p_species(ps)%n_rad_pts = tot_num_points
      p_species(ps)%n_shells  = num_projectors

      ! Skip over the augmentation charges
      write(current_job,'(a80)') 'skipping augmentation charges'
      do shell=1,p_species(ps)%n_shells
         read(2,*,err=100,end=200) line
      end do

       ! ndmh: read in grid sizes for kkbeta and mesh
      write(current_job,'(a80)') 'reading grid sizes'
      if(p_species(ps)%core_charge) then
         read(2,*,err=100,end=200) p_species(ps)%kkbeta,p_species(ps)%mesh
      else
         read(2,*,err=100,end=200) p_species(ps)%kkbeta
         p_species(ps)%mesh = p_species(ps)%kkbeta
      end if

      ! Skip over the augmentation charges
      write(current_job,'(a80)') 'skipping grid points'
      read(2,*,err=100,end=200) (tmp1,i=1,p_species(ps)%mesh)
      read(2,*,err=100,end=200) (tmp2,i=1,p_species(ps)%mesh)

      ! ndmh: skip the L independent function
      write(current_job,'(a80)') 'skipping L independent function'
      do i=1,p_species(ps)%n_shells
         do j=1,i
            read(2,*,err=100,end=200) (tmp1,k=1,p_species(ps)%kkbeta)
         end do
      end do

      ! ndmh: read the number of qfunc coefficients
      read(2,*,err=100,end=200) p_species(ps)%nqfcoef

      ! ndmh: calculate maximum L value of qfcoefs
      p_species(ps)%qf_lmax = 2*lmax

      ! ndmh: calculate inverse of spacing of recip space points
      p_species(ps)%inv_g_spacing = (p_species(ps)%n_rad_pts-1) / &
           p_species(ps)%ps_gmax

      ! ndmh: usps do not include the -Z/r, so flag to not subtract it
      p_species(ps)%subtract_coul = .false.

      ! ** Close the file
      close(2,iostat=ierr,err=500)
500   if (ierr /= 0) then
         write(stdout,'(3a,i6)') &
              'Error in pseudo_species_alloc: closing file "', &
              trim(current_file),'" failed with code ',ierr
         call comms_abort
      end if

      return

100   write(stdout,'(4a)') 'Error in pseudo_species_alloc: reading file "', &
           trim(current_file),'" failed while ',trim(current_job)
      call comms_abort

200   write(stdout,'(4a)') 'Error in pseudo_species_alloc: file "', &
           trim(current_file),'" ended unexpectedly while ',trim(current_job)
      call comms_abort

    end subroutine internal_scan_usp

  end subroutine pseudo_species_alloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_species_read_radial

    !=============================================================!
    ! This subroutine reads the radial parts (local and nonlocal) !
    ! of the pseudopotential for a particular ionic species.      !
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 23/1/2004.              !
    ! Tidied up by Peter D. Haynes on 1/7/2004.                   !
    ! Modified to read multiple file formats (recpots and usps)   !
    ! by Nicholas Hine on 29/01/2009.                             !
    !=============================================================!

    use comms, only: comms_abort, comms_bcast, pub_on_root, pub_root_node_id
    use constants, only: ANGSTROM, HARTREE_IN_EVS
    use rundat, only: pub_smooth_projectors
    use services, only: services_radial_transform
    use simulation_cell, only: pub_cell
    use rundat, only: pub_aug, pub_usp
    use utils, only: utils_alloc_check,utils_dealloc_check, &
         utils_open_unit_check, utils_close_unit_check

    implicit none

    ! Local Variables
    integer :: i
    integer :: ps
    integer :: shell
    integer :: ierr
    real(kind=DP) :: delta
    character(len=64)    :: current_file
    character(len=80)    :: line
    character(len=3)     :: shell_string
    character(len=80)    :: current_job

    if (pub_on_root) then

       write(stdout,'(a)') '<<<<<<<<<<<<<<<<<<<<<<<<< &
            &Pseudopotential information >>>>>>>>>>>>>>>>>>>>>>>>>'

       do ps=1,pub_cell%num_pspecies

          ! pdh: get file name for this species
          current_file = p_species(ps)%pseudo_name

          ! ndmh: determine what file format we are reading and scan it
          ! ndmh: to find sizes of various arrays

          ! ndmh: Old-style CASTEP .recpot
          if (index(current_file,'recpot') > 0) then

             call internal_read_recpot

          ! ndmh: New-style CASTEP OTFG .usp (not necessarily ultrasoft!)
          else if (index(current_file,'usp') > 0) then

             call internal_read_usp

          else
              write(stdout,"(a)")"Error in pseudopotentials_read_species: &
                   & Unrecognised file format: ", current_file
          end if

          ! pdh: print basic information about this pseudopotential
          write(stdout,'(3a,i5,a,f6.1,a)') 'File: ',trim(current_file),' [', &
               p_species(ps)%n_rad_pts,' points up to Gmax=', &
               p_species(ps)%ps_gmax,' (1/bohr)]'
          ! pa: changed from i3 to f10.6
          write(stdout,'(a,i3,a,f10.6)') '  Atomic number:', &
               p_species(ps)%atomic_number,';  ionic charge:', &
               p_species(ps)%ion_charge
          do shell=1,p_species(ps)%n_shells
             write(stdout,'(2(a,i2),a,f5.2,a)') '    Shell',shell,': l =', &
                  p_species(ps)%ang_mom(shell),'; rc =', &
                  p_species(ps)%core_radius(shell),' bohr'
          end do
          if (p_species(ps)%n_shells < 1) write(stdout,'(a)') '    Local potential'
          if (p_species(ps)%core_charge) write(stdout,'(a)') &
               '  Core charge supplied for Nonlinear Core Corrections'

       end do

       write(stdout,'(a/)') '<<<<<<<<<<<<<<<<<<<<<<<<<<&
            &<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

       ! Subtract Coloumb potential from any recpots
       call subtract_or_add_coul

    end if


    ! cks: broadcast results from root to the rest of the processors
    do ps=1,pub_cell%num_pspecies
       call comms_bcast(pub_root_node_id,p_species(ps)%core_radius)
       call comms_bcast(pub_root_node_id,p_species(ps)%ang_mom)
       call comms_bcast(pub_root_node_id,p_species(ps)%D0)
       call comms_bcast(pub_root_node_id,p_species(ps)%rad_locpot_recip)
       call comms_bcast(pub_root_node_id,p_species(ps)%rad_proj_recip)
       if (p_species(ps)%core_charge) then
          call comms_bcast(pub_root_node_id,p_species(ps)%core_charge_recip)
       end if
       ! ndmh: Broadcast augmentation charge info
       call comms_bcast(pub_root_node_id,p_species(ps)%aug_q)
       call comms_bcast(pub_root_node_id,p_species(ps)%qfunc)
       call comms_bcast(pub_root_node_id,p_species(ps)%rinner)
       if (p_species(ps)%nqfcoef > 0) then
          do i=1,p_species(ps)%n_shells
             call comms_bcast(pub_root_node_id,p_species(ps)%qfcoef(:,:,:,i))
          end do
       end if
       call comms_bcast(pub_root_node_id,p_species(ps)%rlog)
       call comms_bcast(pub_root_node_id,p_species(ps)%rab)

    end do

    return

  contains

    subroutine internal_read_recpot

      !==============================================================!
      ! This subroutine reads the pseudopotential for a .recpot      !
      ! and stores the local and nonlocal components of the          !
      ! potential in the arrays allocated in pseudo_alloc_species !
      !--------------------------------------------------------------!
      ! Assembled by Nicholas Hine on 29/01/09 based on code by      !
      ! Chris-Kriton Skylaris from 23/1/2004 and tidied up by        !
      ! Peter Haynes on 1/7/2004.                                    !
      !==============================================================!

      ! Local Variables
      ! cks: for scaling magic
      real(kind=DP) :: g, factor, denom

      call pseudo_open_pspot_file(2,current_file,ierr)
      call utils_open_unit_check('internal_read_recpot &
           &(pseudo_species_read_radial)',trim(current_file),ierr)

      ! * Get to the top of the pseudopotential proper
      current_job = 'seeking marker "END COMMENT"'
      line = repeat(' ',80)
      do while (index(line,'END COMMENT') == 0)
         read(2,'(a80)',err=100,end=200) line  ! To the end of the comments
      end do
      current_job = 'moving to end of header'
      do i=1,2
         read(2,'(a80)',err=100,end=200) line  ! To the end of the header
      end do

      ! * Read the local component of the pseudopotential
      current_job = 'reading local component'
      read(2,*,err=100,end=200) (p_species(ps)%rad_locpot_recip(i), &
           i=1,p_species(ps)%n_rad_pts)


      ! cks: scale radial locpot with suitable factors to make it
      !      consistent with atomic units
      ! cks, 28/3/2004: numerical treatment consistent with CASTEP
      factor = ((ANGSTROM**3)/HARTREE_IN_EVS) / (4.0_DP*PI)
      p_species(ps)%rad_locpot_recip(:) = &
           p_species(ps)%rad_locpot_recip(:) * factor

      ! cks, 28/3/2004: scaling near low g consistent with CASTEP
      g = p_species(ps)%ps_gmax / real(p_species(ps)%n_rad_pts-1,kind=DP)

      ! pa: removed real(...,DP) for fractional ion charges
      factor = 3.75_DP*p_species(ps)%ion_charge / (g*g)

      denom = p_species(ps)%rad_locpot_recip(3) - &
           4.0_DP * p_species(ps)%rad_locpot_recip(2) + &
           3.0_DP * p_species(ps)%rad_locpot_recip(1)

      if (abs(denom) > epsilon(1.0_DP)) factor = factor / denom

      p_species(ps)%rad_locpot_recip(:) = &
           p_species(ps)%rad_locpot_recip(:) * factor

      ! * Read the non-local component of the pseudopotential
      p_species(ps)%D0(:,:) = 0.0_DP
      do shell=1,p_species(ps)%n_shells
         write(shell_string,'(i3)') shell

         write(current_job,'(a80)') 'reading angular momentum of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(2,*,err=100,end=200) p_species(ps)%ang_mom(shell)

         write(current_job,'(a80)') 'reading KB denominator of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(2,*,err=100,end=200) p_species(ps)%D0(shell,shell)

         ! cks: convert to atomic units
         ! cks: scaling consistent with CASTEP
         p_species(ps)%D0(shell,shell) = p_species(ps)%D0(shell,shell) * &
              HARTREE_IN_EVS
         p_species(ps)%D0(shell,shell) = 1.0_DP / p_species(ps)%D0(shell,shell)

         write(current_job,'(a80)') 'reading projector of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(2,*,err=100,end=200) (p_species(ps)%rad_proj_recip(i,shell), &
              i=1,p_species(ps)%n_rad_pts)

         ! cks: convert projector units from A^{3/2}*eV to a0^{3/2}*Eh
         ! cks: scaling consistent with CASTEP
         p_species(ps)%rad_proj_recip(:,shell) = &
              p_species(ps)%rad_proj_recip(:,shell) * sqrt(ANGSTROM**3)

      end do

      ! cks: initialise the radial grid for the interpolation
      ! pdh: inline this call to avoid awkward copy of arguments
      delta = p_species(ps)%ps_gmax / real(p_species(ps)%n_rad_pts-1,DP)
      do i=1,p_species(ps)%n_rad_pts
         p_species(ps)%rad_proj_recip(i,p_species(ps)%n_shells+1) = &
              real(i-1,DP) * delta
      end do

      ! pdh: get core radii from projectors
      call internal_get_core_radii(p_species(ps))

      ! cks: Gauss filter projectors in reciprocal and real space
      if (pub_smooth_projectors >= 0.0_DP ) then
         call pseudo_gauss_filter_proj(p_species(ps))
      endif

      ! ndmh: find end of projectors
      current_job = 'seeking marker "1000"'
      do while ((index(line,'1000') == 0).and.(len_trim(adjustl(line)) /= 4))
         read(2,'(a80)',err=100,end=200) line
      end do

      ! ndmh: get NLCC core charge if present
      if (p_species(ps)%core_charge) then
         read(2,*,err=100,end=200) (p_species(ps)%core_charge_recip(i), &
              i=1,p_species(ps)%n_rad_pts)
      end if

      ! ndmh: recpots do not use kkbeta or mesh
      p_species(ps)%rlog = 0.0_DP
      p_species(ps)%rab = 0.0_DP
      p_species(ps)%qfunc = 0.0_DP
      p_species(ps)%aug_q = 0.0_DP

      close(2,iostat=ierr,err=400)
400   call utils_close_unit_check('internal_read_recpot &
           &(pseudo_species_read_radial)',trim(current_file),ierr)

      return

100   write(stdout,'(4a)') 'Error in internal_read_recpot (pseudo_species_read_radial):',&
           'reading file "',trim(current_file),'" failed while ',trim(current_job)
      call comms_abort

200   write(stdout,'(4a)') 'Error in internal_read_recpot (pseudo_species_read_radial):',&
           'file "',trim(current_file),'" ended unexpectedly while ',trim(current_job)
      call comms_abort

    end subroutine internal_read_recpot

    subroutine internal_read_usp

      !==============================================================!
      ! This subroutine reads the pseudopotential for a .usp         !
      ! and stores the local and nonlocal components of the          !
      ! potential in the arrays allocated in pseudo_alloc_species !
      !--------------------------------------------------------------!
      ! Written by Nicholas Hine on 29/01/09                         !
      !==============================================================!

      ! Local Variables
      integer,parameter :: max_proj=huge(1)
      integer :: j, k, l
      integer :: lshell
      real(kind=DP) :: factor
      real(kind=DP), allocatable :: core_charge_real(:)

      call pseudo_open_pspot_file(2,current_file,ierr)
      call utils_open_unit_check('internal_read_usp &
           &(pseudo_species_read_radial)',trim(current_file),ierr)

      ! ndmh: get to the top of the pseudopotential proper
      current_job = 'seeking marker "END COMMENT"'
      line = repeat(' ',80)
      do while (index(line,'END COMMENT') == 0)
         read(2,'(a80)',err=100,end=200) line  ! To the end of the comments
      end do
      current_job = 'moving to end of header'
      do i=1,3
         read(2,'(a80)',err=100,end=200) line  ! To the end of the header
      end do

      ! ndmh: read the local component of the pseudopotential
      read(2,*,err=100,end=200) (p_species(ps)%rad_locpot_recip(i), &
           i=1,p_species(ps)%n_rad_pts)

      ! ndmh: Scale radial locpot with suitable factors to convert it
      ! ndmh: to atomic units
      factor = ((ANGSTROM**3)/HARTREE_IN_EVS) / (4.0_DP*PI)
      p_species(ps)%rad_locpot_recip(:) = &
           p_species(ps)%rad_locpot_recip(:) * factor

      ! ndmh: Read the non-local component of the pseudopotential
      p_species(ps)%D0(:,:) = 0.0_DP
      do shell=1,p_species(ps)%n_shells

         write(shell_string,'(i3)') shell

         write(current_job,'(a80)') 'reading angular momentum of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(2,*,err=100,end=200) p_species(ps)%ang_mom(shell)

         ! Determine lower limit of shells of this l
         l = p_species(ps)%ang_mom(shell)
         do lshell=1,p_species(ps)%n_shells
            if (p_species(ps)%ang_mom(lshell)==l) exit
         end do

         ! Read KB energies
         write(current_job,'(a80)') 'reading KB energy of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(2,*,err=100,end=200) (p_species(ps)%D0(i,shell), &
              i=lshell,shell)

         ! Fill in other half of D0 matrix
         do i=lshell,shell
            p_species(ps)%D0(shell,i) = p_species(ps)%D0(i,shell)
         end do

         write(current_job,'(a80)') 'reading projector of shell '// &
              adjustl(shell_string)
         current_job = adjustl(current_job)
         read(2,*,err=100,end=200) (p_species(ps)%rad_proj_recip(i,shell), &
              i=1,p_species(ps)%n_rad_pts)

         ! cks: convert projector units from A^{3/2}*eV to a0^{3/2}*Eh
         ! cks: scaling consistent with CASTEP
         p_species(ps)%rad_proj_recip(:,shell) = &
              p_species(ps)%rad_proj_recip(:,shell) * sqrt(ANGSTROM**3)

      end do

      ! Convert units of D0
      p_species(ps)%D0(:,:) = p_species(ps)%D0(:,:) / HARTREE_IN_EVS

      ! cks: initialise the radial grid for the interpolation
      ! pdh: inline this call to avoid awkward copy of arguments
      delta = p_species(ps)%ps_gmax / real(p_species(ps)%n_rad_pts-1,DP)
      do i=1,p_species(ps)%n_rad_pts
         p_species(ps)%rad_proj_recip(i,p_species(ps)%n_shells+1) = &
              real(i-1,DP) * delta
      end do

      ! ndmh: check we are at the right point of the file
      !    -- There is no '1000' flag if there are 2 projectors per
      !       all channels up to l=3 (the assumed maximum)
      write(current_job,'(a80)') 'reading "1000" '
      if(p_species(ps)%n_shells < max_proj) then
         read(2,'(a)',err=100,end=200) line
         if(trim(adjustl(line))/='1000') then
            write(stdout,'(4a)') 'Error in internal_read_usp '//&
                '(pseudo_species_read_radial):','file "',trim(current_file), &
                '" ended unexpectedly while ',trim(current_job)
            call comms_abort
         end if
      end if

      ! ndmh: Read in augmentation charges
      write(current_job,'(a80)') 'reading augmentation charges'
      p_species(ps)%aug_q(:,:) = 0.0_DP
      do shell=1,p_species(ps)%n_shells
         l = p_species(ps)%ang_mom(shell)
         ! Determine lower limit of shells of this l
         do lshell=1,p_species(ps)%n_shells
            if (p_species(ps)%ang_mom(lshell)==l) exit
         end do
         ! Read in the q-values
         read(2,*,err=100,end=200) (p_species(ps)%aug_q(i,shell), &
              i=lshell,shell)
         ! Fill in other half of matrix
         do i=lshell,shell
            p_species(ps)%aug_q(shell,i) = p_species(ps)%aug_q(i,shell)
         end do
      end do
      if (any(abs(p_species(ps)%aug_q(:,:))>tiny(1.0_DP))) then
         write(stdout,'(a)') 'Error in internal_read_usp &
              &(pseudo_species_read_radial):'
         write(stdout,'(a,i4)') 'Nonzero augmentation charge found'
         write(stdout,'(a)') 'Ultrasoft Pseudopotentials with augmentation &
              &charges are not supported.'
         call comms_abort
         pub_usp = .true.
         pub_aug = .true.
      end if

      ! ndmh: read in grid sizes for kkbeta and mesh
      write(current_job,'(a80)') 'reading grid sizes'
      if(p_species(ps)%core_charge) then
         read(2,*,err=100,end=200) p_species(ps)%kkbeta,p_species(ps)%mesh
      else
         read(2,*,err=100,end=200) p_species(ps)%kkbeta
         p_species(ps)%mesh = p_species(ps)%kkbeta
      end if

      ! * Allocate temporary storage for core_charge_real if required
      if (p_species(ps)%core_charge) then
         allocate(core_charge_real(p_species(ps)%mesh),stat=ierr)
         call utils_alloc_check('internal_read_usp '//&
              '(pseudo_species_read_radial)','core_charge_real',ierr)
      end if

      ! ndmh: read in the grid points
      write(current_job,'(a80)') 'reading grid points'
      read(2,*,err=100,end=200) (p_species(ps)%rlog(i),i=1,p_species(ps)%mesh)
      read(2,*,err=100,end=200) (p_species(ps)%rab(i),i=1,p_species(ps)%mesh)

      ! ndmh: read the L independent function
      ! ndmh: this is zero for ncpp's
      write(current_job,'(a80)') 'reading L independent function'
      do i=1,p_species(ps)%n_shells
         do j=1,i
            read(2,*,err=100,end=200) (p_species(ps)%qfunc(k,i,j),k=1, &
                 p_species(ps)%kkbeta)
         end do
      end do
      do i=1,p_species(ps)%n_shells
         do j=i+1,p_species(ps)%n_shells
            p_species(ps)%qfunc(:,i,j) = p_species(ps)%qfunc(:,j,i)
         end do
      end do

      ! ndmh: the L dependent coeffs
      ! ndmh: these are zero for ncpp's
      p_species(ps)%qfcoef(:,:,:,:) = 0.0_DP
      write(current_job,'(a80)') 'reading L dependent coeffs'
      p_species(ps)%qf_lmax = 2*maxval(p_species(ps)%ang_mom( &
           1:p_species(ps)%n_shells))
      read(2,*,err=100,end=200) p_species(ps)%nqfcoef
      do l=0,p_species(ps)%qf_lmax
         read(2,*,err=100,end=200) p_species(ps)%rinner(l)
         do i=1,p_species(ps)%n_shells
            do j=1,i
               read(2,*,err=100,end=200) (p_species(ps)%qfcoef(k,l,i,j), &
                    k=1,p_species(ps)%nqfcoef)
            end do
         end do
      end do
      do l=0,p_species(ps)%qf_lmax
         do i=1,p_species(ps)%n_shells
            do j=i+1,p_species(ps)%n_shells
               p_species(ps)%qfcoef(:,l,j,i) = p_species(ps)%qfcoef(:,l,i,j)
            end do
         end do
      end do


      ! ndmh: finally read in the core charge on the log grid
      if (p_species(ps)%core_charge) then
         read(2,*,err=100,end=200) &
              (core_charge_real(i),i=1,p_species(ps)%mesh)

         ! ndmh: transform the real-space core charge to reciprocal space
         call services_radial_transform(0,0,p_species(ps)%mesh,p_species(ps)%rlog, &
              p_species(ps)%rab,p_species(ps)%n_rad_pts,&
              p_species(ps)%ps_gmax,core_charge_real,&
              p_species(ps)%core_charge_recip)

      end if

      ! pdh: get core radii from projectors
      call internal_get_core_radii(p_species(ps))

      ! cks: Gauss filter projectors in reciprocal and real space
      if (pub_smooth_projectors >= 0.0_DP ) then
         call pseudo_gauss_filter_proj(p_species(ps))
      endif

      ! ndmh: deallocate temporary storage
      if (p_species(ps)%core_charge) then
         deallocate(core_charge_real,stat=ierr)
         call utils_dealloc_check('internal_read_usp '//&
              '(pseudo_species_read_radial)','core_charge_real',ierr)
      end if

      ! ndmh: close the file
      close(2,iostat=ierr,err=400)
400   call utils_close_unit_check('internal_read_usp &
           &(pseudo_species_read_radial)',trim(current_file),ierr)

      return

100   write(stdout,'(4a)') 'Error in internal_read_usp '//&
           '(pseudo_species_read_radial):','reading file "',&
           trim(current_file),'" failed while ',trim(current_job)
      call comms_abort

200   write(stdout,'(4a)') 'Error in internal_read_usp '//&
           '(pseudo_species_read_radial):','file "',trim(current_file), &
           '" ended unexpectedly while ',trim(current_job)
      call comms_abort

    end subroutine internal_read_usp

    subroutine internal_get_core_radii(p_species)

      !==================================================================!
      ! This subroutine calculates the core radii of the pseudopotential !
      ! species by examining the projectors in real space.               !
      !------------------------------------------------------------------!
      ! Written by Peter D. Haynes in early 2004.                        !
      ! Tidied up by Peter D. Haynes on 1/7/2004.                        !
      !==================================================================!

      use services, only: services_sbessj
      use utils, only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Argument
      type(PSEUDO_SPECIES), intent(inout) :: p_species

      ! Local variables
      integer :: ierr                           ! Error flag
      integer :: i_shell                        ! Shell counter
      integer :: l                              ! Angular momentum of shell
      integer :: g_rad_pt,r_rad_pt              ! Grid point counters
      integer :: n_r_rad_pts                    ! Num points on real radial grid
      real(kind=DP), parameter :: rmax = 5.0_DP ! Maximum core radius
      real(kind=DP), parameter :: halfpi = 0.5_DP * PI
      real(kind=DP) :: gnorm,rnorm              ! Norms in recip/real space
      real(kind=DP) :: proj,proj_max            ! Projector (max value) in rspc
      real(kind=DP) :: val,fac                  ! Variables for integration
      real(kind=DP) :: delr,r,delg,g            ! Grid variables
      real(kind=DP), allocatable :: integrand(:)! gspc integrand for FT

      ! Set parameters for radial grid
      delg = p_species%ps_gmax / (p_species%n_rad_pts - 1)
      delr = pi / p_species%ps_gmax
      n_r_rad_pts = nint(rmax / delr) + 1
      fac = (7.0_DP/17280.0_DP) * delg / sqrt(halfpi)

      ! Allocate workspace for integration
      allocate(integrand(p_species%n_rad_pts),stat=ierr)
      call utils_alloc_check('internal_get_core_radii &
           &(pseudo_species_read_radial)','integrand',ierr)

      ! Loop over shells: angular momentum is l
      do i_shell=1,p_species%n_shells
         l = p_species%ang_mom(i_shell)

         ! Calculate norm from reciprocal space representation
         gnorm = 0.0_DP
         do g_rad_pt=1,p_species%n_rad_pts
            g = (g_rad_pt-1) * delg
            val = p_species%rad_proj_recip(g_rad_pt,i_shell) * g
            gnorm = gnorm + val*val
         end do
         gnorm = gnorm * delg

         if (gnorm < epsilon(1.0_DP)) then
            p_species%core_radius(i_shell) = 0.0_DP
            cycle
         end if

         ! Calculate projector on real-space radial grid, and cumulative
         ! norm of projector in real space
         rnorm = 0.0_DP
         proj_max = 0.0_DP

         do r_rad_pt=1,n_r_rad_pts
            r = (r_rad_pt-1) * delr
            p_species%core_radius(i_shell) = 1.8_DP
            integrand(1) = 0.0_DP
            do g_rad_pt=2,p_species%n_rad_pts
               g = (g_rad_pt-1) * delg
               integrand(g_rad_pt) = services_sbessj(l,g*r) * g*g * &
                    p_species%rad_proj_recip(g_rad_pt,i_shell)
            end do

            ! Due to the oscillatory nature of the integrand,
            ! an eight point Newton-Cotes integration method is used.
            ! See  Abramowitz and Stegun Eq. 25.4.17

            proj = 0.0_DP
            do g_rad_pt=1,p_species%n_rad_pts-7,7
               proj = proj + &
                    751.0_DP*(integrand(g_rad_pt)+integrand(g_rad_pt+7)) + &
                    3577.0_DP*(integrand(g_rad_pt+1)+integrand(g_rad_pt+6)) + &
                    1323.0_DP*(integrand(g_rad_pt+2)+integrand(g_rad_pt+5)) + &
                    2989.0_DP*(integrand(g_rad_pt+3)+integrand(g_rad_pt+4))
            end do
            proj = proj * fac
            val = r * proj
            rnorm = rnorm + val*val * delr

            ! Set core radius if magnitude of projector is less then 1% of
            ! maximum value and cumulative norm is at least 99.9% of
            ! reciprocal space norm

            proj_max = max(proj_max,abs(proj))
            if (r_rad_pt > 1) then
               if (abs(proj/proj_max) < 1.0e-2_DP .and. &
                    rnorm/gnorm > 0.999_DP) then
                  p_species%core_radius(i_shell) = r
                  exit
               end if
            end if

         end do

      end do

      ! Free workspace
      deallocate(integrand,stat=ierr)
      call utils_dealloc_check('internal_get_core_radii &
           &(pseudo_species_read_radial)','integrand',ierr)

    end subroutine internal_get_core_radii

  end subroutine pseudo_species_read_radial


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_open_pspot_file(iunit,filename,ierr)

    !========================================================================!
    ! This subroutine tries to open a pseudopotential file: first it tries   !
    ! the current directory. If the file does not exist, it then checks the  !
    ! environment variable PSPOT_DIR and looks for tbe file there (ACCELRYS  !
    ! version only).                                                         !
    !------------------------------------------------------------------------!
    ! Arguments                                                              !
    !   iunit (in): unit number to open file with.                           !
    !   filename (in): name of the pseudopotential file to try to open       !
    !   ierr (out): error flag                                               !
    !------------------------------------------------------------------------!
    ! Written by Victor Milman in May 2010.                                  !
    ! Cleaned up and formatted for ONETEP by Nicholas Hine in June 2010.     !
    !========================================================================!

    implicit none

    ! Arguments
    integer,intent(in) :: iunit
    character(len=*) :: filename
    integer,intent(out) :: ierr

    ! Local Variables
    logical :: file_exists
#ifdef ACCELRYS
    character(len=500) :: pspot_dir
#endif

    inquire(file=trim(filename),exist=file_exists)

    if (file_exists) then
       open(iunit,file=trim(filename),status='old',position='rewind', &
            iostat=ierr)
#ifdef ACCELRYS
    else
       ! Check if there is a defined pseudopotentials directory in the
       ! environment variable PSPOT_DIR
       call getenv('PSPOT_DIR',pspot_dir)
       open(iunit,file=trim(pspot_dir)//'/'//trim(filename),status='old', &
            position='rewind',iostat=ierr)
#endif
    end if

  end subroutine pseudo_open_pspot_file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_gauss_filter_proj(current_species)

    !=========================================================+=========!
    ! This subroutine applies a Gaussian smoothing filter to a          !
    ! non-local projector while it is still in radial form.             !
    !-------------------------------------------------------------------!
    !   *Description of the algorithm*                                  !
    ! First the smoothing is applied in reciprocal space.               !
    ! The part of the projector from 0.667g_grid to g_max is multiplied !
    ! by a Gaussian which begins with the value 1.0 at 0.667*g_grid     !
    ! and decays with 0.1 of an exponent whose half-width which is      !
    ! equal to pub_smooth_projectors*g_grid.                            !
    ! Then the projector is Fourier transformed to radial real space    !
    ! and again filtered with a Gaussian in the region r_core -> infty. !
    ! The halfwidth of the Gaussian is now equal to                     !
    ! pub_smooth_projectors*r_core.                                     !
    ! Finally the projector is put back to reciprocal space.            !
    !-------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 20/01/2005 using parts of     !
    ! Peter Haynes' subroutine "internal_get_core_radii".               !
    !===================================================================!

    use rundat, only: pub_smooth_projectors
    use services, only: services_sbessj
    use simulation_cell, only : pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Argument
    type(PSEUDO_SPECIES), intent(inout) :: current_species

    ! Local variables
    integer :: ierr                           ! Error flag
    integer :: i_shell                        ! Shell counter
    integer :: l                              ! Angular momentum of shell
    integer :: g_rad_pt,r_rad_pt              ! Grid point counters
    integer :: n_r_rad_pts                    ! Num points on real radial grid
    integer :: max_npts         ! maximum number of real and recip gridpoints
    real(kind=DP), parameter :: rmax = 10.0_DP ! Maximum real space distance
    real(kind=DP), parameter :: halfpi = 0.5_DP * PI
    real(kind=DP) :: proj                     ! Projector in rspc
    real(kind=DP) :: fac                      ! Variables for integration
    real(kind=DP) :: delr, delg               ! Grid variables
    real(kind=DP), allocatable :: integrand(:)! gspc integrand for FT
    real(kind=DP), allocatable :: rspace_projector (:)
    real(kind =DP) :: rval      ! grid distance in real space
    real(kind =DP) :: gval      ! grid distance in reciprocal space
    real(kind =DP) :: g_grid    ! 1-D maximum g-vector of psinc-grid
    real(kind =DP) :: fac_real  ! factor to multiply real space integral
    real(kind =DP) :: alpha     ! exponent of Gaussian
    real(kind =DP) :: core_rad  ! core radius of projector

    ! Set parameters for radial grid
    delg = current_species%ps_gmax / (current_species%n_rad_pts - 1)
    delr = PI / current_species%ps_gmax

    ! cks: Overkill to retain all precission in real-> recip FT
    delr =0.1_DP*delr

    n_r_rad_pts = nint(rmax / delr) + 1

    max_npts =max(n_r_rad_pts, current_species%n_rad_pts)

    g_grid =3.0_DP*PI/(pub_cell%d1 + pub_cell%d2 +pub_cell%d3)

    fac = (7.0_DP/17280.0_DP) * delg  / halfpi
    fac_real =(7.0_DP/17280.0_DP) * delr

    if (pub_smooth_projectors > tiny(1.0_DP) ) then
       alpha = 0.1_DP*log(2.0_DP)/( (pub_smooth_projectors*g_grid)**2 )
    else
       alpha = 0.0_DP
    end if

    ! Allocate workspace for integration
    allocate(integrand(max_npts),stat=ierr)
    call utils_alloc_check('pseudo_gauss_filter_proj','integrand',ierr)

    ! cks: Allocate rspace_projector
    allocate(rspace_projector(n_r_rad_pts),stat=ierr)
    call utils_alloc_check('pseudo_gauss_filter_proj', &
         'rspace_projector',ierr)

    ! Loop over shells: angular momentum is l
    do i_shell=1, current_species%n_shells
       l = current_species%ang_mom(i_shell)


       ! cks: ++++++++ Gauss filter projector in reciprocal space +++++
       do g_rad_pt =1, current_species%n_rad_pts

          gval =current_species%rad_proj_recip(g_rad_pt, current_species%n_shells+1)

          if ( gval .gt. (0.667_DP*g_grid ) ) then
             if (pub_smooth_projectors > tiny(1.0_DP) ) then
                current_species%rad_proj_recip(g_rad_pt, i_shell) = &
                     exp( -alpha*(( gval -0.667_DP*g_grid )**2) ) &
                     * current_species%rad_proj_recip(g_rad_pt, i_shell)
             else
                current_species%rad_proj_recip(g_rad_pt, i_shell) =0.0_DP
             endif
          endif

       enddo
       ! cks: +++++ END Gauss filter projector in reciprocal space +++++




       ! cks: ====== RECIP to REAL transform ===========================

       ! cks: Calculate projector on real-space radial grid
       do r_rad_pt =1, n_r_rad_pts

          rval = (r_rad_pt-1) * delr
          integrand(1) = 0.0_DP
          do g_rad_pt=2, current_species%n_rad_pts

             gval =current_species%rad_proj_recip(g_rad_pt, current_species%n_shells+1)

             integrand(g_rad_pt) = services_sbessj(l,gval*rval) * gval*gval * &
                  current_species%rad_proj_recip(g_rad_pt,i_shell)

          end do

          ! Due to the oscillatory nature of the integrand,
          ! an eight point Newton-Cotes integration method is used.
          ! See  Abramowitz and Stegun Eq. 25.4.17
          proj = 0.0_DP
          do g_rad_pt=1, current_species%n_rad_pts-7,7
             proj = proj + &
                  751.0_DP*(integrand(g_rad_pt)+integrand(g_rad_pt+7)) + &
                  3577.0_DP*(integrand(g_rad_pt+1)+integrand(g_rad_pt+6)) + &
                  1323.0_DP*(integrand(g_rad_pt+2)+integrand(g_rad_pt+5)) + &
                  2989.0_DP*(integrand(g_rad_pt+3)+integrand(g_rad_pt+4))
          end do
          proj = proj * fac

          rspace_projector(r_rad_pt) =proj

       end do

       ! cks: == END RECIP to REAL transform ===========================



       ! cks: ++++++++ Gauss filter projector in real space +++++
       core_rad =current_species%core_radius(i_shell)

       if (pub_smooth_projectors > tiny (1.0_DP) ) then
          alpha =log(2.0_DP)/( (pub_smooth_projectors*core_rad)**2 )
       endif


       do r_rad_pt =1, n_r_rad_pts

          rval = (r_rad_pt-1) * delr

          if ( rval .gt. core_rad ) then
             if (pub_smooth_projectors > tiny (1.0_DP) ) then
                rspace_projector(r_rad_pt) = &
                     exp( -alpha*(( rval -core_rad )**2) ) &
                     * rspace_projector(r_rad_pt)
             else
                rspace_projector(r_rad_pt) =0.0_DP
             endif
          endif

       enddo
       ! cks: ++++ END Gauss filter projector in real space +++++



       ! cks: ====== REAL to RECIP transform ===========================

       ! cks: put projector back to reciprocal space
       do g_rad_pt =1, current_species%n_rad_pts

          gval =current_species%rad_proj_recip(g_rad_pt, current_species%n_shells+1)
          integrand(1) = 0.0_DP
          do r_rad_pt= 2, n_r_rad_pts
             rval = (r_rad_pt-1) * delr
             integrand(r_rad_pt) = services_sbessj(l,gval*rval) *rval *rval * &
                  rspace_projector(r_rad_pt)
          end do


          ! Due to the oscillatory nature of the integrand,
          ! an eight point Newton-Cotes integration method is used.
          ! See  Abramowitz and Stegun Eq. 25.4.17
          proj = 0.0_DP
          do r_rad_pt= 1, n_r_rad_pts -7,7
             proj = proj + &
                  751.0_DP*(integrand(r_rad_pt)+integrand(r_rad_pt+7)) + &
                  3577.0_DP*(integrand(r_rad_pt+1)+integrand(r_rad_pt+6)) + &
                  1323.0_DP*(integrand(r_rad_pt+2)+integrand(r_rad_pt+5)) + &
                  2989.0_DP*(integrand(r_rad_pt+3)+integrand(r_rad_pt+4))
          end do

          current_species%rad_proj_recip(g_rad_pt,i_shell) =proj *fac_real
       end do
       ! cks: == END REAL to RECIP transform ===========================

    end do

    ! Free workspace
    deallocate(rspace_projector, stat=ierr)
    call utils_dealloc_check('pseudo_gauss_filter_proj', &
         'rspace_projector',ierr)

    ! Free workspace
    deallocate(integrand,stat=ierr)
    call utils_dealloc_check('pseudo_gauss_filter_proj','integrand',ierr)

  end subroutine pseudo_gauss_filter_proj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_species_init_proj(proj_set,elements)

    !=================================================================!
    ! This subroutine allocates and initialises all fftbox_proj_recip !
    ! for each p_species element. Each such fftbox_proj_recip is      !
    ! a full non-local projector in reciprocal space in the FFTbox.   !
    ! Optionally, this subroutine can also calculate instead Del_G of !
    ! the projectors along a given cartesian direction given by cart  !
    !-----------------------------------------------------------------!
    ! Original version by Chris-Kriton Skylaris, January 2004.        !
    ! Tidied up by Peter D. Haynes, July 2004.                        !
    ! Re-write by Nicholas Hine, creating PROJECTOR_SET type in       !
    ! July 2009.                                                      !
    ! Del_G code added by Laura Ratcliff in March 2011 as a separate  !
    ! routine (pseudo_species_calc_grad_proj)                         !
    ! Merged back into this routine to avoid code duplication by      !
    ! Nicholas Hine in April 2011.                                    !
    !=================================================================!

    use comms, only: pub_on_root
    use constants, only: VERBOSE
    use geometry, only: POINT, OPERATOR(.dot.)
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_orig_atom
    use projectors, only: PROJECTOR_SET, projectors_allocate_set, &
         projectors_init_fftbox_recip
    use rundat, only: pub_output_detail
    use simulation_cell, only : pub_cell, pub_fftbox
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    type(PROJECTOR_SET), intent(inout) :: proj_set

    ! Local Variables
    integer :: iat, orig_iat
    integer :: ps
    integer :: shell
    integer :: proj_count

    if (pub_on_root  .and. pub_output_detail == VERBOSE) &
         write(stdout,'(a)') &
         '... Nonlocal pseudopotential projector initialisation'

    ! ndmh: find number of unique projectors
    proj_count = 0
    do ps=1,pub_cell%num_pspecies
       do shell=1,p_species(ps)%n_shells
          proj_count = proj_count + 2*p_species(ps)%ang_mom(shell) + 1
       end do
    end do

    ! Set up the entries of proj_set
    proj_set%n_proj_species = pub_cell%num_pspecies

    if (.not. allocated(proj_set%fftbox_proj_recip)) then
       call projectors_allocate_set(proj_set, &
            maxval(p_species(:)%n_shells),maxval(p_species(:)%n_rad_pts))
    else
       call utils_abort('Error in pseudo_species_init_proj: &
            &proj_set already allocated')
    end if

    ! ndmh: set species_num_proj and species_first_proj values
    proj_count = 1
    do ps=1,proj_set%n_proj_species
       proj_set%species_num_proj(ps) = 0
       proj_set%species_first_proj(ps) = proj_count
       proj_set%gmax(ps) = p_species(ps)%ps_gmax
       proj_set%n_rad_pts(ps) = p_species(ps)%n_rad_pts
       proj_set%num_shells(ps) = p_species(ps)%n_shells
       proj_set%ang_mom(:,ps) = 0
       proj_set%rad_proj_recip(:,:,ps) = 0.0_DP
       do shell=1,p_species(ps)%n_shells
          proj_set%ang_mom(shell,ps) = p_species(ps)%ang_mom(shell)
          proj_set%species_num_proj(ps) = &
               proj_set%species_num_proj(ps) + &
               2*p_species(ps)%ang_mom(shell) + 1
          proj_set%rad_proj_recip(1:p_species(ps)%n_rad_pts,shell,ps) = &
               p_species(ps)%rad_proj_recip(1:p_species(ps)%n_rad_pts,shell)
       end do
       proj_count = proj_count + proj_set%species_num_proj(ps)
    end do

    ! ndmh: copy projector centre and radius from elements array
    do iat=1,pub_cell%nat
       orig_iat = pub_orig_atom(iat)
       proj_set%proj_centre(iat) = elements(orig_iat)%centre
       proj_set%proj_max_radius(iat) = elements(orig_iat)%max_core_radius
       proj_set%proj_species(iat) = elements(orig_iat)%pspecies_number
    end do

    ! ndmh: Initialise projectors in reciprocal-space fftbox representation
    !call projectors_init_fftbox_recip(proj_set,kpt,delta,cart,swap_rc)

  end subroutine pseudo_species_init_proj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudopotentials_species_exit

    !==========================================================!
    ! This subroutine deallocates all the memory allocated for !
    ! the p_species array.                                     !
    !----------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/1/2004.           !
    ! Tidied up by Peter D. Haynes on 1/7/2004.                !
    ! Modified to also deallocate nlps_projectors by Nicholas  !
    ! Hine on 02/11/2009.                                      !
    !==========================================================!

    use comms, only: comms_abort
    use projectors, only: projectors_deallocate_set
    use rundat, only: pub_any_nl_proj
    use utils, only: utils_dealloc_check

    implicit none

    ! Local Variables
    integer :: ierr
    integer :: sp
    character(len=32) :: string

    if (pub_any_nl_proj) call projectors_deallocate_set(nlps_projectors)

    if (allocated(p_species)) then

       do sp=size(p_species),1,-1

          if (allocated(p_species(sp)%aug_q)) then
             deallocate(p_species(sp)%aug_q,stat=ierr)
             write(string,'(a,i2,a)') 'p_species(',sp,')%aug_q'
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  trim(string),ierr)
          end if

          if (allocated(p_species(sp)%qfunc)) then
             deallocate(p_species(sp)%qfunc,stat=ierr)
             write(string,'(a,i2,a)') 'p_species(',sp,')%qfunc'
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  trim(string),ierr)
          end if

          if (allocated(p_species(sp)%qfcoef)) then
             deallocate(p_species(sp)%qfcoef,stat=ierr)
             write(string,'(a,i2,a)') 'p_species(',sp,')%qfcoef'
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  trim(string),ierr)
          end if

          if (allocated(p_species(sp)%rinner)) then
             deallocate(p_species(sp)%rinner,stat=ierr)
             write(string,'(a,i2,a)') 'p_species(',sp,')%rinner'
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  trim(string),ierr)
          end if

          if (allocated(p_species(sp)%rlog)) then
             deallocate(p_species(sp)%rlog,stat=ierr)
             write(string,'(a,i2,a)') 'p_species(',sp,')%rlog'
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  trim(string),ierr)
          end if

          if (allocated(p_species(sp)%rab)) then
             deallocate(p_species(sp)%rab,stat=ierr)
             write(string,'(a,i2,a)') 'p_species(',sp,')%rab'
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  trim(string),ierr)
          end if

          if (p_species(sp)%core_charge) then
             deallocate(p_species(sp)%core_charge_recip,stat=ierr)
             write(string,'(a,i2,a)') 'p_species(',sp,')%core_charge_recip'
             call utils_dealloc_check('pseudopotentials_species_exit', &
                  trim(string),ierr)
          end if

          deallocate(p_species(sp)%rad_proj_recip,stat=ierr)
          write(string,'(a,i2,a)') 'p_species(',sp,')%rad_proj_recip'
          call utils_dealloc_check('pseudopotentials_species_exit', &
               trim(string),ierr)

          deallocate(p_species(sp)%rad_locpot_recip,stat=ierr)
          write(string,'(a,i2,a)') 'p_species(',sp,')%rad_locpot_recip'
          call utils_dealloc_check('pseudopotentials_species_exit', &
               trim(string),ierr)

          deallocate(p_species(sp)%D0,stat=ierr)
          write(string,'(a,i2,a)') 'p_species(',sp,')%D0'
          call utils_dealloc_check('pseudopotentials_species_exit', &
               trim(string),ierr)

          deallocate(p_species(sp)%core_radius,stat=ierr)
          write(string,'(a,i2,a)') 'p_species(',sp,')%core_radius'
          call utils_dealloc_check('pseudopotentials_species_exit', &
               trim(string),ierr)

          deallocate(p_species(sp)%ang_mom,stat=ierr)
          write(string,'(a,i2,a)') 'p_species(',sp,')%ang_mom'
          call utils_dealloc_check('pseudopotentials_species_exit', &
               trim(string),ierr)

       end do

       deallocate(p_species,stat=ierr)
       call utils_dealloc_check('pseudopotentials_species_exit', &
            'p_species',ierr)

    end if

  end subroutine pseudopotentials_species_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! SUBROUTINES FOR INTERFACING WITH THE ATOM SOLVER !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_locpot_rad(locpot,npts,rad,isp)

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

    ! Local Variables
    real(kind=DP), allocatable :: work(:)
    real(kind=DP) :: q,dq,Z
    integer :: npts_q
    integer :: ir,iq
    integer :: ierr

    npts_q = p_species(isp)%n_rad_pts
    allocate(work(npts_q),stat=ierr)
    call utils_alloc_check('pseudo_get_locpot_rad','work',ierr)

    dq = p_species(isp)%ps_gmax / real(npts_q-1,kind=DP)
    Z = p_species(isp)%ion_charge

    do ir=1,npts
       do iq=1,npts_q
          q = real(iq-1,kind=DP)*dq
          work(iq) = (p_species(isp)%rad_locpot_recip(iq)*q**2 - Z) * &
               services_sbessj(0,q*rad(ir))
       end do

       locpot(ir) = work(1)+work(npts_q)
       do iq = 2,npts_q-1,2
           locpot(ir) = locpot(ir) + 4.0_DP*work(iq)+2.0_DP*work(iq+1)
       enddo
       locpot(ir) = locpot(ir)*dq/3.0_DP
    end do

    locpot(:) = locpot(:)*(2.0_DP/PI)

    deallocate(work,stat=ierr)
    call utils_dealloc_check('pseudo_get_locpot_rad','work',ierr)

  end subroutine pseudo_get_locpot_rad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_core_den_rad(core_den,npts,rad,isp)

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

    if (.not.p_species(isp)%core_charge) then
       core_den(:) = 0.0_DP
       return
    end if

    npts_q = p_species(isp)%n_rad_pts
    allocate(work(npts_q),stat=ierr)
    call utils_alloc_check('pseudo_get_core_den_rad','work',ierr)

    dq = p_species(isp)%ps_gmax / real(npts_q-1,kind=DP)

    do ir=1,npts
       do iq=1,npts_q
          q = real(iq-1,kind=DP)*dq
          work(iq) = p_species(isp)%core_charge_recip(iq)*q**2 * &
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
    call utils_dealloc_check('pseudo_get_core_den_rad','work',ierr)

  end subroutine pseudo_get_core_den_rad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_projector_info(isp,nshells,nproj_tot,lmax)

    !============================================================!
    ! This function fetches the number of shells of projectors.  !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                !
    !============================================================!

    ! Arguments
    integer,intent(in) :: isp
    integer,intent(out) :: nshells
    integer,intent(out) :: nproj_tot
    integer,intent(out) :: lmax

    ! Local Variables
    integer :: ishell

    nshells = p_species(isp)%n_shells
    nproj_tot = 0
    do ishell=1,nshells
       nproj_tot = nproj_tot + 2*p_species(isp)%ang_mom(ishell) + 1
    end do
    lmax = maxval(p_species(isp)%ang_mom)

  end subroutine pseudo_get_projector_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_aug_funcs(qijl,qfunc,nshells,npts,rad,isp)

    !============================================================!
    ! This subroutine fetches the q function and the qijL        !
    ! terms for each angular momentum channel.                   !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in February 2011.                 !
    !============================================================!

    use services, only: services_locate_interp, services_linear_interpolation

    ! Arguments
    integer, intent(in) :: npts
    integer, intent(in) :: nshells
    real(kind=DP), intent(in) :: rad(npts)
    integer, intent(in) :: isp
    real(kind=DP), intent(out) :: qijL(nshells,nshells)
    real(kind=DP), intent(out) :: qfunc(npts,nshells,nshells)

    ! Local variables
    integer :: ipt,jpt
    integer :: ishell,jshell

    ! Set the augmentation charges
    qijL(:,:) = p_species(isp)%aug_q(:,:)

    if (p_species(isp)%kkbeta==0) then
       qfunc(:,:,:) = 0.0_DP
       return
    end if

    ! Set the L-independent function
    do ishell=1,nshells
       do jshell=1,nshells
          do ipt=1,npts
             if (rad(ipt)>p_species(isp)%rlog(p_species(isp)%kkbeta-1)) then
                qfunc(ipt,ishell,jshell) = 0.0_DP
                cycle
             end if
             jpt = services_locate_interp(rad(ipt), &
                  p_species(isp)%rlog,p_species(isp)%mesh)
             qfunc(ipt,ishell,jshell) = services_linear_interpolation(rad(ipt), &
                  p_species(isp)%qfunc(jpt,ishell,jshell), &
                  p_species(isp)%qfunc(jpt+1,ishell,jshell), &
                  p_species(isp)%rlog(jpt),p_species(isp)%rlog(jpt+1))
          end do
       end do
    end do

  end subroutine pseudo_get_aug_funcs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_projectors_q(proj_q,dij,ang_mom,nshells,nsws, &
       nsws_max,lmax,qb,isp)

    !============================================================!
    ! This subroutine fetches the projectors at chosen q-points  !
    !------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                !
    !============================================================!

    use comms, only: comms_abort
    use services, only: services_1d_interpolation

    implicit none

    ! Arguments
    integer, intent(in) :: nshells
    integer, intent(in) :: lmax
    integer, intent(in) :: nsws(0:lmax)
    integer, intent(in) :: nsws_max
    integer, intent(in) :: isp
    integer, intent(out) :: ang_mom(nshells)
    real(kind=DP), intent(in) :: qb(nsws_max,0:lmax)
    real(kind=DP), intent(out) :: proj_q(nsws_max,nshells)
    real(kind=DP), intent(out) :: dij(nshells,nshells)

    ! Local variables
    integer :: ishell, isw
    real(kind=DP) :: qq

    ! Set the KB denominators
    dij(:,:) = p_species(isp)%D0(:,:)

    ! Interpolate the reciprocal-space projectors at the qb values
    do ishell=1,nshells
       proj_q(:,ishell) = 0.0_DP
       ang_mom(ishell) = p_species(isp)%ang_mom(ishell)
       do isw=1,nsws(ang_mom(ishell))
          qq = qb(isw,p_species(isp)%ang_mom(ishell))
          proj_q(isw,ishell) = services_1d_interpolation( &
               p_species(isp)%rad_proj_recip(:,ishell), &
               p_species(isp)%n_rad_pts, &
               qq*p_species(isp)%inv_g_spacing, &
               p_species(isp)%ang_mom(ishell))
       end do
    end do

  end subroutine pseudo_get_projectors_q



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! SUBROUTINES FOR THE LOCAL PSEUDOPOTENTIAL !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_make_structure_fac_old(struct_fac,    &  ! output
       elements, grid)

    !==============================================================!
    ! This subroutine generates the structure factor each distinct !
    ! ionic pseudopotential species.                               !
    !--------------------------------------------------------------!
    ! Written by Arash A. Mostofi 19/2/2004.                       !
    ! Modified by Chris-Kriton Skylaris on 22/2/2004 so that       !
    ! struct_fac is parallelised in memory.                        !
    ! Re-written by Peter D. Haynes on 30/6/2004 to use data       !
    ! parallelisation according to the new Fourier routines.       !
    ! Modified by Nicholas Hine on 19/09/2008 to optimise cache    !
    ! performance by re-ordering array                             !
    !==============================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: comms_barrier, pub_my_node_id
    use ion, only: ELEMENT
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)   :: grid
    type(ELEMENT), intent(in)     :: elements(pub_cell%nat)
    complex(kind=DP), intent(out) :: struct_fac(pub_cell%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)

    ! Local variables
    integer :: i3,i2,islab23
    integer :: species,atom
    integer :: ierr
    real(kind=DP) :: gvec(3),gdotr

    real(kind=DP),allocatable :: species_pos(:,:,:)
    integer, allocatable :: species_nat(:)

    ! Start timer
    call timer_clock('pseudo_make_structure_fac_old',1)

    ! ndmh: allocate storage for number of atoms in each species
    allocate(species_nat(pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('pseudo_make_structure_fac_old','species_nat',ierr)

    ! ndmh: count number of atoms in each species
    species_nat = 0
    do atom=1,pub_cell%nat
       species = elements(atom)%pspecies_number
       species_nat(species) = species_nat(species) + 1
    end do

    ! ndmh: allocate temporary storage for species atom positions
    allocate(species_pos(3,maxval(species_nat),pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('pseudo_make_structure_fac_old','species_pos',ierr)

    ! ndmh: copy atom positions to temporary array for cache-efficiency
    species_nat = 0
    do atom=1,pub_cell%nat
       species = elements(atom)%pspecies_number
       species_nat(species) = species_nat(species) + 1
       species_pos(1,species_nat(species),species) = elements(atom)%centre%x
       species_pos(2,species_nat(species),species) = elements(atom)%centre%y
       species_pos(3,species_nat(species),species) = elements(atom)%centre%z
    end do

    ! Loop over reciprocal space on this node
    do islab23=1,grid%num_slabs23                ! along b1
       do i2=1,grid%n2                           ! along b2
          ! ndmh: initialise this slab of the structure factor array
          struct_fac(:,:,i2,islab23) = (0.0_DP,0.0_DP)
          do i3=1,grid%n3                        ! along b3
          
             call cell_grid_recip_pt(gvec,islab23 + &
                  grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)

             ! ndmh: loop over species
             do species=1,pub_cell%num_pspecies
                ! ndmh: loop over atoms of this species
                do atom=1,species_nat(species)
                   gdotr = gvec(1)*species_pos(1,atom,species) + &
                        gvec(2)*species_pos(2,atom,species) + &
                        gvec(3)*species_pos(3,atom,species)
                   struct_fac(species,i3,i2,islab23) = &
                        struct_fac(species,i3,i2,islab23) + &
                        cmplx(cos(gdotr),-sin(gdotr),kind=DP)
                end do
             end do

          end do      ! loop along b3
       end do         ! loop along b2
    end do            ! loop along b1

    ! ndmh: this routine has no comms, so is not necessarily synchronised.
    ! ndmh: synchronise here before returning, otherwise timers report
    ! ndmh: inaccurate timings for next function called, due to load imbalance
    ! ndmh: in this routine
    call comms_barrier

    ! ndmh: deallocate temporary storage
    deallocate(species_pos,stat=ierr)
    call utils_dealloc_check('pseudo_make_structure_fac_old','species_pos',ierr)
    deallocate(species_nat,stat=ierr)
    call utils_dealloc_check('pseudo_make_structure_fac_old','species_nat',ierr)

    ! Stop timer
    call timer_clock('pseudo_make_structure_fac_old',2)

  end subroutine pseudo_make_structure_fac_old



  subroutine pseudo_make_structure_factor(struct_fac,    &  ! output
       elements, grid)

    !==============================================================!
    ! This subroutine generates the structure factor each distinct !
    ! ionic pseudopotential species.                               !
    !--------------------------------------------------------------!
    ! Written by Arash A. Mostofi 19/2/2004.                       !
    ! Modified by Chris-Kriton Skylaris on 22/2/2004 so that       !
    ! struct_fac is parallelised in memory.                        !
    ! Re-written by Peter D. Haynes on 30/6/2004 to use data       !
    ! parallelisation according to the new Fourier routines.       !
    ! Modified by Nicholas Hine on 19/09/2008 to optimise cache    !
    ! performance by re-ordering array                             !
    !==============================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: comms_barrier, pub_my_node_id
    use constants, only: cmplx_i, cmplx_0
    use ion, only: ELEMENT
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)   :: grid
    type(ELEMENT), intent(in)     :: elements(pub_cell%nat)
    complex(kind=DP), intent(out) :: struct_fac(pub_cell%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)

    ! Local variables
    integer :: i3,i2,islab23
    integer :: species,atom
    integer :: ierr
    real(kind=DP) :: gpt(3)
    real(kind=DP),allocatable :: species_pos(:,:,:)
    complex(kind=DP),allocatable :: phase_components(:,:,:,:)
    integer, allocatable :: species_nat(:)
    real(kind=DP) :: g1r1, g2r2, g3r3
    complex(kind=DP) :: accum

    ! Start timer
    call timer_clock('pseudo_make_structure_factor',1)

    ! jd: Take care of padding, the loop ignores it
    struct_fac = cmplx_0

    ! ndmh: allocate storage for number of atoms in each species
    allocate(species_nat(pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('pseudo_make_structure_factor','species_nat',ierr)

    ! ndmh: count number of atoms in each species
    species_nat = 0
    do atom=1,pub_cell%nat
       species = elements(atom)%pspecies_number
       species_nat(species) = species_nat(species) + 1
    end do

    ! ndmh: allocate temporary storage for species atom positions
    allocate(species_pos(3,maxval(species_nat),pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('pseudo_make_structure_factor','species_pos',ierr)

    ! jd: allocate temporary for components of the phase factor
    allocate(phase_components(3,maxval(species_nat),pub_cell%num_pspecies, &
         max(grid%num_slabs23,grid%n2,grid%n3)),stat=ierr)
    call utils_alloc_check('pseudo_make_structure_factor','phase_components',ierr)

    ! ndmh: copy atom positions to temporary array for cache-efficiency
    species_nat = 0
    do atom=1,pub_cell%nat
       species = elements(atom)%pspecies_number
       species_nat(species) = species_nat(species) + 1
       species_pos(1,species_nat(species),species) = elements(atom)%centre%x
       species_pos(2,species_nat(species),species) = elements(atom)%centre%y
       species_pos(3,species_nat(species),species) = elements(atom)%centre%z
    end do

    ! jd: Determine exp(-i g1*r1) for all atoms
    do islab23=1,grid%num_slabs23                ! along b1
       do species=1,pub_cell%num_pspecies
          do atom=1,species_nat(species)
             call cell_grid_recip_pt(gpt(1:3),islab23 + &
                  grid%first_slab23(pub_my_node_id) - 1,1,1,grid)
             g1r1 = sum(gpt(1:3)*species_pos(1:3,atom,species))
             phase_components(1,atom,species,islab23) = &
                  cmplx(cos(g1r1),-sin(g1r1),kind=DP)
          end do
       end do
    end do

    ! jd: Determine exp(-i g2*r2) for all atoms
    do i2=1,grid%n2                              ! along b2
       do species=1,pub_cell%num_pspecies
          do atom=1,species_nat(species)
             call cell_grid_recip_pt(gpt(1:3),1,i2,1,grid)
             g2r2 = sum(gpt(1:3)*species_pos(1:3,atom,species))
             phase_components(2,atom,species,i2) = &
                  cmplx(cos(g2r2),-sin(g2r2),kind=DP)
          end do
       end do
    end do

    ! jd: Determine exp(-i g3*r3) for all atoms
    do i3=1,grid%n3                             ! along b3
       do species=1,pub_cell%num_pspecies
          do atom=1,species_nat(species)
             call cell_grid_recip_pt(gpt(1:3),1,1,i3,grid)
             g3r3 = sum(gpt(1:3)*species_pos(1:3,atom,species))
             phase_components(3,atom,species,i3) = &
                  cmplx(cos(g3r3),-sin(g3r3),kind=DP)
          end do
       end do
    end do

    ! Loop over reciprocal space on this node
    ! jd: Use exp(-i g.r) = exp(-i g1*r1) * exp(-i g2*r2) * exp(-i g3*r3)
    !     to avoid calculating trigs every time.
    do islab23=1,grid%num_slabs23                ! along b1
       do i2=1,grid%n2                           ! along b2
          do i3=1,grid%n3                        ! along b3

             do species=1,pub_cell%num_pspecies
                accum = cmplx_0
                do atom=1,species_nat(species)
                   accum = accum + &
                        phase_components(1,atom,species,islab23) * &
                        phase_components(2,atom,species,i2) * &
                        phase_components(3,atom,species,i3)
                end do
                struct_fac(species,i3,i2,islab23) = accum
             end do

          end do      ! loop along b3
       end do         ! loop along b2
    end do            ! loop along b1

    ! ndmh: this routine has no comms, so is not necessarily synchronised.
    ! ndmh: synchronise here before returning, otherwise timers report
    ! ndmh: inaccurate timings for next function called, due to load imbalance
    ! ndmh: in this routine
    call comms_barrier

    ! ndmh: deallocate temporary storage
    deallocate(species_pos,stat=ierr)
    call utils_dealloc_check('pseudo_make_structure_factor','species_pos',ierr)
    deallocate(species_nat,stat=ierr)
    call utils_dealloc_check('pseudo_make_structure_factor','species_nat',ierr)
    deallocate(phase_components,stat=ierr)
    call utils_dealloc_check('pseudo_make_structure_factor','phase_components',ierr)

    ! Stop timer
    call timer_clock('pseudo_make_structure_factor',2)

  end subroutine pseudo_make_structure_factor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudopotentials_local_on_grid(v_local_fine,&
        struct_fac,struct_fac_classical,grid,elements)

    !==================================================================!
    ! This subroutine generates the local part of the total ionic      !
    ! pseudopotential on the fine grid simulation cell.                !
    !------------------------------------------------------------------!
    !                                                                  !
    !------------------------------------------------------------------!
    ! This subroutine was originally written by Chris-Kriton Skylaris  !
    ! in 2000.                                                         !
    ! Modified by A. A. Mostofi in 19/02/2004 so that it uses          !
    ! complex-complex FFTs.                                            !
    ! Modified by Chris-Kriton Skylaris on 22/2/2004 so that           !
    ! struct_fac is memory-parallelised                                !
    ! Modified for parallelisation with new fourier routines by        !
    ! Peter D. Haynes on 30/6/2004.                                    !
    ! Modified to include smeared-ion contribution if required, which  !
    ! necessitated passing elements, by Jacek Dziedzic on 13/05/2010.  !
    !==================================================================!

    use cell_grid, only: GRID_INFO
    use classical_pot, only: classical_pot_recip
    use fourier, only: fourier_apply_cell_backward
    use ion, only: ELEMENT
    use is_smeared_ions, only: smeared_ion_apply_vloc_corr, &
         smeared_ion_initialise
    use rundat, only: pub_is_smeared_ion_rep, pub_open_localpseudo
    use simulation_cell, only : pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(inout) :: grid
    complex(kind=DP), intent(in) :: struct_fac(pub_cell%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)
    complex(kind=DP), intent(in) :: struct_fac_classical( &
         grid%ld3,grid%ld2,grid%max_slabs23)
    real(kind=DP), intent(out) :: v_local_fine(grid%ld1, &
         grid%ld2, grid%max_slabs12)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! Local variables
    integer :: ierr              ! Error flag
    complex(kind=DP), allocatable :: v_local_recip(:,:,:)

    call timer_clock('pseudopotentials_local_on_grid',1)

    ! jd: This is a good place to initialise smeared ions, if in use
    if (pub_is_smeared_ion_rep) call smeared_ion_initialise(elements)

    ! jd: If using open boundary conditions, use a different method to get
    !     the vloc in real space
    if(pub_open_localpseudo) then
       call pseudo_local_on_grid_openbc(v_local_fine,grid,elements)

    else ! jd: Usual PBC method

       ! cks: expand to 3D in reciprocal fine grid and sum together
       !      (along with the structure factor) the local pseudopotentials
       !      for each species to obtain the total local ionic potential in
       !      reciprocal representation on the fine grid
       allocate(v_local_recip(grid%ld3,grid%ld2, &
            grid%max_slabs23),stat=ierr)
       call utils_alloc_check('pseudopotentials_local_on_grid', &
            'v_local_recip',ierr)

       call pseudopotentials_sum_local_rec(v_local_recip, & ! output
            struct_fac,grid)    ! input

       ! cks: include external potential from "classical" atoms
       if (pub_cell%nat_classical > 0) then
          call classical_pot_recip(v_local_recip, & ! output
               struct_fac_classical,grid)    ! input
       endif

       ! FFT the local ionic potential from reciprocal to real space
       call fourier_apply_cell_backward(v_local_fine,v_local_recip,grid)

       deallocate(v_local_recip,stat=ierr)
       call utils_dealloc_check('pseudopotentials_local_on_grid', &
            'v_local_recip',ierr)

    end if

    ! jd: Include the correction due to smeared ions, if any, in v_local_fine
    if (pub_is_smeared_ion_rep) then
       call smeared_ion_apply_vloc_corr(elements, v_local_fine)
    end if

    call timer_clock('pseudopotentials_local_on_grid',2)

  end subroutine pseudopotentials_local_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudopotentials_sum_local_rec(fine_complex,struct_fac,grid)

    !===============================================================!
    ! This subroutine generates in reciprocal space the total local !
    ! potential in the simulation cell due to the ionic             !
    ! pseudopotentials of all ions.                                 !
    !---------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris in 2000.          !
    ! Improvements by Chris-Kriton Skylaris on 30/5/2001.           !
    ! Modified by Chris-Kriton Skylaris on 24/1/2004 to work with   !
    ! the pseudo_species type.                                      !
    ! Modified by Arash A. Mostofi in Aprl 2003 to be compatible    !
    ! with complex-to-complex FFTs (rather than complex-to-real).   !
    ! Modified by Arash A. Mostofi on 19/2/2004 to work with        !
    ! structure factor as input.                                    !
    ! Modified by Chris-Kriton Skylaris on 22/2/2004 to work with   !
    ! data-parallel structure factor.                               !
    ! Modified for parallelisation with new fourier routines by     !
    ! Peter D. Haynes on 30/6/2004.                                 !
    ! Modified by Nicholas Hine in December 2009 to allow less than !
    ! one slab23 per node.                                          !
    !===============================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: comms_abort, pub_my_node_id, pub_total_num_nodes
    use pbc_corrections, only: pbc_corr_vloc
    use rundat, only: pub_mt_cutoff, pub_smooth_loc_pspot
    use services, only: services_1d_interpolation
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(GRID_INFO), intent(inout)   :: grid
    complex(kind=DP), intent(in)  :: struct_fac(pub_cell%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)
    complex(kind=DP), intent(out) :: fine_complex(grid%ld3,&
         grid%ld2, grid%max_slabs23)

    ! Local variables
    integer :: species              ! Atomic species counter
    integer :: i3,i2,islab23        ! Reciprocal grid loop counters
    real(kind=DP) :: gvec(3)        ! G-vector
    real(kind=DP) :: g_length       ! Length of this G vector
    real(kind=DP) :: v_loc_value    ! Local potential at this G
    real(kind=DP) :: fourpi         ! Constant multiplier
    real(kind=DP) :: g_cut, alpha   ! For filtering locps

    fourpi = 4.0_DP * PI

    if (pub_smooth_loc_pspot > 0.0_DP) then
       g_cut = 3.0_DP*PI/(pub_cell%d1+pub_cell%d2+pub_cell%d3)*(4.0_DP/7.0_DP)
       alpha = 0.1_DP*log(2.0_DP)/((pub_smooth_loc_pspot*g_cut)**2)
    else
       g_cut = huge(1.0_DP)
       alpha = 0.0_DP
    end if

    ! Loop over reciprocal space grid on this node
    do islab23=1,grid%num_slabs23       ! along b1
       do i2=1,grid%n2                  ! along b2
          do i3=1,grid%n3               ! along b3

             ! Get G-vector
             call cell_grid_recip_pt(gvec,islab23 + &
                  grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)

             ! Get magnitude of this G-vector
             g_length = sqrt(sum(gvec(1:3)**2))

             ! ndmh: initialise (cache-efficiently)
             fine_complex(i3,i2,islab23) = cmplx(0.0_DP,0.0_DP,kind=DP)

             ! Loop over atomic species
             do species=1,pub_cell%num_pspecies

                ! Get pseudopotential at this G-vector
                v_loc_value = services_1d_interpolation( &
                     p_species(species)%rad_locpot_recip, &
                     p_species(species)%n_rad_pts,&
                     g_length*p_species(species)%inv_g_spacing,0)

                ! Filter high G components if required
                if (g_length > g_cut) then
                   v_loc_value = v_loc_value * exp(-alpha*(g_length-g_cut)**2)
                end if

                ! Add back Coulomb potential
                v_loc_value = v_loc_value - p_species(species)%ion_charge * &
                      grid%coulomb_recip(i3,i2,islab23)

                ! ndmh: scale by 4 pi/fine_weight
                ! jd: split the scaling into two steps
                v_loc_value = v_loc_value * fourpi

                ! jd: Apply the Martyna-Tuckerman correction, if requested
                if (pub_mt_cutoff /= 0.0_DP) then
                   call pbc_corr_vloc(i3,i2,islab23, &
                        p_species(species)%ion_charge,v_loc_value,grid)
                end if

                ! jd: remaining part of scaling
                v_loc_value = v_loc_value / grid%weight

                fine_complex(i3,i2,islab23) = fine_complex(i3,i2,islab23) + &
                     struct_fac(species,i3,i2,islab23) * v_loc_value

             end do    ! loop over species

          end do   ! b3
       end do      ! b2
    end do         ! b1


    ! G=0 element must be real
    if (pub_my_node_id==grid%node_slab23(1)) then
       if (aimag(fine_complex(1,1,1)) /= 0.0_DP) then
          write(stdout,'(a)') 'Error in pseudopotentials_sum_local_rec: &
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
    if (pub_my_node_id==grid%node_slab23(grid%n1/2+1)) &
         fine_complex(:,:,grid%num_slabs23) = (0.0_DP,0.0_DP)

  end subroutine pseudopotentials_sum_local_rec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudopotentials_core_density(core_density_fine,struct_fac,grid)

    !==================================================================!
    ! This subroutine reconstructs the core density for the whole      !
    ! supercell on the fine real space grid using the information      !
    ! stored for each pseudopotential species in the array             !
    ! core_charge_recip.                                               !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !                                                                  !
    ! 1) core_density_fine  : output  : data-parallelised core density !
    ! 2) struct_fac : input : data-parallelised structure factor       !
    ! 3) grid : input : description of whole-cell grid                 !
    !                                                                  !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 29/01/09.           !
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
    real(kind=DP), intent(out) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    complex(kind=DP), intent(in) :: struct_fac(pub_cell%num_pspecies, &
         grid%ld3,grid%ld2,grid%max_slabs23)

    ! Local variables
    integer :: ierr              ! Error flag
    complex(kind=DP), allocatable :: core_density_recip(:,:,:)
    real(kind=DP) :: gvec(3), g_length, core_den_value, factor
    integer :: species              ! Atomic species counter
    integer :: i3,i2,islab23        ! Reciprocal grid loop counters

    call timer_clock('pseudopotentials_core_density',1)

    ! ndmh: allocate storage for core density in reciprocal space
    allocate(core_density_recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('pseudopotentials_core_density', &
         'core_density_recip',ierr)

    ! ndmh: loop over reciprocal space grid on this node
    do islab23=1,grid%num_slabs23       ! along b1
       do i2=1,grid%n2                  ! along b2
          do i3=1,grid%n3               ! along b3

             ! Initialise
             core_density_recip(i3,i2,islab23) = (0.0_DP,0.0_DP)

             ! Loop over atomic species
             do species=1,pub_cell%num_pspecies

                ! Check if we have a core charge for this species
                if (.not.p_species(species)%core_charge) cycle

                ! Get length of this G-vector
                call cell_grid_recip_pt(gvec,islab23 + &
                     grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)
                g_length = sqrt(sum(gvec(1:3)**2))

                ! Get core density at this G-vector
                core_den_value = services_1d_interpolation( &
                     p_species(species)%core_charge_recip, &
                     p_species(species)%n_rad_pts,&
                     g_length*p_species(species)%inv_g_spacing,0)

                core_density_recip(i3,i2,islab23) = &
                     core_density_recip(i3,i2,islab23) + &
                     struct_fac(species,i3,i2,islab23) * core_den_value

             end do    ! loop over species

          end do   ! b3
       end do      ! b2
    end do         ! b1

    ! FFT the core density from reciprocal to real space
    call fourier_apply_cell_backward(core_density_fine,core_density_recip,grid)

    ! ndmh: scale with 1.0/weight
    factor = 1.0_DP / grid%weight
    core_density_fine = factor * core_density_fine

    ! ndmh: deallocate storage for core density in reciprocal space
    deallocate(core_density_recip,stat=ierr)
    call utils_dealloc_check('pseudopotentials_core_density', &
         'core_density_recip',ierr)

    call timer_clock('pseudopotentials_core_density',2)

  end subroutine pseudopotentials_core_density


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_local_calculate_forces(den_slabs12,grid,elements, &
       locps_forces)

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the local part of the ionic pseudopotential.                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) den_slabs12  : input  : ground state data-parallelised charge density!
    ! 2) elements     : input  : list of elements and corresponding info      !
    ! 3) locps_forces : output : ionic forces due to local part of pseudopot  !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) p_species    : pseudopotential information of all ionic species      !
    !-------------------------------------------------------------------------!
    ! Written by Arash A. Mostofi, v0.0, 28th June 2004                       !
    ! Modified by Nicholas Hine on 04/09/2009 to re-order the loops for cache !
    ! efficiency and to pre-calculate v(G) for each species for each G vector.!
    ! Modified by Nicholas Hine in December 2009 to allow for less than one   !
    ! 23-slab per processor.                                                  !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: comms_reduce, pub_my_node_id, pub_on_root
    use fourier, only: fourier_apply_cell_forward
    use geometry, only: POINT, OPERATOR(.dot.)
    use ion, only: ELEMENT
    use rundat, only : pub_smooth_loc_pspot
    use services, only: services_1d_interpolation
    use simulation_cell, only : pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: den_slabs12(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_cell%num_spins)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    real(kind=DP), intent(out) :: locps_forces(1:3,pub_cell%nat)

    ! Local Variables
    integer :: i2,i3,islab23
    integer :: species,atom
    integer :: ierr
    real(kind=DP) :: gvec(3)
    real(kind=DP) :: g_length,gdotR
    real(kind=DP) :: factor,v_loc_value
    real(kind=DP) :: g_cut
    complex(kind=DP) :: eiGR
    complex(kind=DP), allocatable :: iG_vden(:,:)
    complex(kind=DP), allocatable :: recip(:,:,:)
    logical :: filter_locps

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering pseudo_local_calculate_forces'
#endif

    ! Start timer
    call timer_clock('pseudo_local_calculate_forces',1)

    ! Allocate
    allocate(iG_vden(1:3,pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('pseudo_local_calculate_forces','iG_vden',ierr)
    allocate(recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('pseudo_local_calculate_forces','recip',ierr)

    ! Initialise
    locps_forces = 0.0_DP

    if (pub_smooth_loc_pspot > 0.0_DP) then
       filter_locps = .true.
       g_cut = 3.0_DP*PI/(pub_cell%d1+pub_cell%d2+pub_cell%d3)*(4.0_DP/7.0_DP)
    else
       filter_locps = .false.
       g_cut = huge(1.0_DP)
    end if

    ! If spin polarised, put total density in up spin
    if (pub_cell%num_spins == 2) den_slabs12(:,:,:,1) = den_slabs12(:,:,:,1) + &
         den_slabs12(:,:,:,2)

    ! Fourier transform the charge density to reciprocal space
    call fourier_apply_cell_forward(den_slabs12(:,:,:,1),recip,grid)

    ! For components g with symmetry points at -g
    factor=2.0_DP

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

             ! ndmh: loop over species, calculating v(G) for each.
             ! ndmh: this bit is O(N) so as much as possible should be
             ! ndmh: pre-calculated here.
             do species=1,pub_cell%num_pspecies

                ! Interpolate value of local potential at current g
                v_loc_value = services_1d_interpolation( &
                     p_species(species)%rad_locpot_recip, &
                     p_species(species)%n_rad_pts,g_length * &
                     p_species(species)%inv_g_spacing,0)

                ! Filter high G components if required
                if (filter_locps.and.(g_length > g_cut)) then
                   v_loc_value = v_loc_value * &
                        exp(-pub_smooth_loc_pspot*(g_length/g_cut-1.0_DP)**2)
                end if

                ! Add back the Coulomb potential; set g=0 term to zero
                if (g_length .gt. 0.0_DP) then
                   ! pa: changed to allow fractional ionic charge
                   v_loc_value = v_loc_value - &
                        p_species(species)%ion_charge * &
                        grid%coulomb_recip(i3,i2,islab23)
                else
                   v_loc_value = 0.0_DP
                endif

                ! ndmh: calculate iG.v(G).n*(G) for each species and
                ! ndmh: each Cartesian direction
                iG_vden(1:3,species) = factor * cmplx(0.0_DP,1.0_DP) &
                     * gvec * v_loc_value * conjg(recip(i3,i2,islab23))

             end do

             ! ndmh: loop over atoms: this bit is O(N^2) so should be made as
             ! ndmh: short and simple as possible
             do atom=1,pub_cell%nat

                ! e^{-ig.R}
                gdotR = -(gvec(1)*elements(atom)%centre%x + &
                     gvec(2)*elements(atom)%centre%y + &
                     gvec(3)*elements(atom)%centre%z)
                eiGR = cmplx(COS(gdotR),SIN(gdotR),kind=DP)

                ! Sum force over g in each Cartesian direction i
                ! ==>  f = sum_{g} i.g.e^{-ig.R}.Vlocps(g).den^{*}(g)
                locps_forces(:,atom) = locps_forces(:,atom) + &
                     real(iG_vden(:,elements(atom)%pspecies_number)*eiGR,kind=DP)

             enddo     ! End loop over atoms

          enddo     ! End loop along b3

       enddo     ! End loop along b2

       ! For g1/=0 slabs
       factor=2.0_DP

    enddo     ! End loop along b1

    ! Sum the result over all nodes
    call comms_reduce('SUM',locps_forces)

    ! Deallocate
    deallocate(recip,stat=ierr)
    call utils_dealloc_check('pseudo_local_calculate_forces','recip',ierr)
    deallocate(iG_vden,stat=ierr)
    call utils_dealloc_check('pseudo_local_calculate_forces','iG_vden',ierr)

    ! Scale
    locps_forces = locps_forces * 4.0_DP * PI / (grid%n1*grid%n2*grid%n3)

    ! If spin polarised, restore up spin density
    if (pub_cell%num_spins == 2) den_slabs12(:,:,:,1) = den_slabs12(:,:,:,1) - &
         den_slabs12(:,:,:,2)

    ! Stop timer
    call timer_clock('pseudo_local_calculate_forces',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving pseudo_local_calculate_forces'
#endif

    return
  end subroutine pseudo_local_calculate_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_nlcc_calculate_forces(density_fine,core_density_fine, &
       grid,elements,nlcc_forces)

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the nonlinear core correction core charge.                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) density_fine : input  : ground state charge density                  !
    ! 2) core_density_fine  : input  : NLCC core charge                       !
    ! 3) elements     : input  : list of elements and corresponding info      !
    ! 4) nlcc_forces : output : ionic forces due to NLCC corrections          !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) p_species    : pseudopotential information of all ionic species      !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 6th February 2009                             !
    ! based on pseudo_local_calculate_forces by Arash Mostofi                 !
    ! Modified by Nicholas Hine on 04/09/2009 to re-order the loops for cache !
    ! efficiency and to pre-calculate n(G) for each species for each G vector.!
    ! Modified by Nicholas Hine in December 2009 to allow for less than one   !
    ! 23-slab per processor.                                                  !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: comms_reduce, pub_my_node_id, pub_on_root
    use fourier, only: fourier_apply_cell_forward
    use geometry, only: POINT, OPERATOR(.dot.)
    use ion, only: ELEMENT
    use services, only: services_1d_interpolation
    use simulation_cell, only : pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use xc, only : xc_energy_potential

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: density_fine(grid%ld1,grid%ld2,&
         grid%max_slabs12,pub_cell%num_spins)
    real(kind=DP), intent(in) :: core_density_fine(:,:,:)
    ! jd: NB: usually ld1,ld2,max_slabs12, but might be 1x1x1 if no NLCC
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    real(kind=DP), intent(out) :: nlcc_forces(1:3,pub_cell%nat)

    ! Local Variables
    integer :: i2,i3,islab23,is
    integer :: species,atom
    integer :: ierr
    real(kind=DP) :: gvec(3)
    real(kind=DP) :: g_length
    real(kind=DP) :: factor
    real(kind=DP) :: xc_energy
    real(kind=DP) :: coreden, GdotR
    real(kind=DP), allocatable :: total_density_fine(:,:,:,:)
    real(kind=DP), allocatable :: xc_pot_fine(:,:,:,:)
    complex(kind=DP), allocatable :: iG_coreden_vxc(:,:)
    complex(kind=DP) :: eiGR
    complex(kind=DP), allocatable :: recip(:,:,:)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering pseudo_nlcc_calculate_forces'
#endif

    ! Start timer
    call timer_clock('pseudo_nlcc_calculate_forces',1)

    ! Allocate
    allocate(total_density_fine(grid%ld1,grid%ld2,grid%max_slabs12, &
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('pseudo_nlcc_calculate_forces','total_density_fine',&
         ierr)
    allocate(xc_pot_fine(grid%ld1,grid%ld2,grid%max_slabs12, &
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('pseudo_nlcc_calculate_forces','xc_pot_fine',ierr)
    allocate(recip(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('pseudo_nlcc_calculate_forces','recip',ierr)
    allocate(iG_coreden_vxc(1:3,pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('pseudo_nlcc_calculate_forces','iG_coreden_vxc',ierr)

    ! Calculate total density
    factor = 1.0_DP / pub_cell%num_spins
    do is=1,pub_cell%num_spins
       total_density_fine(:,:,:,is) = density_fine(:,:,:,is) &
            + core_density_fine*factor
    end do

    ! Calculate the exchange correlation potential
    call xc_energy_potential(total_density_fine,xc_energy,xc_pot_fine,grid,0)

    ! Average the spin channels in (:,:,:,1)
    if (pub_cell%num_spins == 2) then
       xc_pot_fine(:,:,:,1) = factor*(xc_pot_fine(:,:,:,1) &
            + xc_pot_fine(:,:,:,2))
    end if

    ! Initialise
    nlcc_forces = 0.0_DP

    ! Fourier transform the xc potential to reciprocal space
    call fourier_apply_cell_forward(xc_pot_fine(:,:,:,1),recip,grid)

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

             ! g-vector, |g| and 1/|g|^2
             call cell_grid_recip_pt(gvec,islab23 + &
                  grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)

             g_length = sqrt(sum(gvec(1:3)**2))

             ! Loop over atoms to find n(G) for each
             do species=1,pub_cell%num_pspecies

                ! Check if we actually have a core charge for this species
                if (.not.p_species(species)%core_charge) cycle

                ! Interpolate value of core density at current g
                coreden = services_1d_interpolation( &
                     p_species(species)%core_charge_recip, &
                     p_species(species)%n_rad_pts, &
                     g_length*p_species(species)%inv_g_spacing,0)

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
                if (.not.p_species(species)%core_charge) cycle

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
    call utils_dealloc_check('pseudo_nlcc_calculate_forces','iG_coreden_vxc',ierr)
    deallocate(recip,stat=ierr)
    call utils_dealloc_check('pseudo_nlcc_calculate_forces','recip',ierr)
    deallocate(xc_pot_fine,stat=ierr)
    call utils_dealloc_check('pseudo_nlcc_calculate_forces','xc_pot_fine',ierr)
    deallocate(total_density_fine,stat=ierr)
    call utils_dealloc_check('pseudo_nlcc_calculate_forces', &
         'total_density_fine',ierr)

    ! Scale factor
    nlcc_forces = nlcc_forces / real(grid%n1*grid%n2*grid%n3,DP)

    ! Stop timer
    call timer_clock('pseudo_nlcc_calculate_forces',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving pseudo_nlcc_calculate_forces'
#endif

    return
  end subroutine pseudo_nlcc_calculate_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! SUBROUTINES FOR THE NON-LOCAL PSEUDOPOTENTIAL !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_dij(dij)

    !=======================================================================!
    ! This subroutine finds the Kleinman-Bylander denominators for the      !
    ! nonlocal pseudopotential projectors on this node and puts 1/kb as the !
    ! diagonal elements of the matrix kb_denominators.                      !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine on 15 June 2009.                             !
    !=======================================================================!

    use comms, only: pub_my_node_id
    use parallel_strategy, only: pub_first_atom_on_node, &
         pub_num_atoms_on_node, pub_elements_on_node
    use simulation_cell, only : pub_cell
    use sparse, only: SPAM3, sparse_put_block
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: dij

    ! Local Variables
    integer :: iat, loc_iat
    integer :: isp
    integer :: ishell, jshell
    integer :: iproj, jproj
    integer :: li, lj
    integer :: mi, mj
    integer :: ierr
    integer :: max_atom_proj
    real(kind=DP), allocatable :: dij0_sp(:,:,:)
    complex(kind=DP), allocatable :: dij0_sp_cmplx(:,:,:)
    logical :: iscmplx

    ! ndmh: initialisations
    iscmplx = dij%iscmplx
    max_atom_proj = maxval(nlps_projectors%species_num_proj(:))

    allocate(dij0_sp(max_atom_proj,max_atom_proj, &
         pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('pseudo_get_dij','dij0_sp',ierr)
    dij0_sp(:,:,:) = 0.0_DP

    ! Loop over species, setting up array of dij0 terms
    do isp=1,pub_cell%num_pspecies
       if (any(pub_elements_on_node(:)%pspecies_number==isp)) then
          iproj = 0
          do ishell=1,p_species(isp)%n_shells
             li = p_species(isp)%ang_mom(ishell)
             do mi=-li,li
                iproj = iproj + 1
                jproj = 0
                do jshell=1,p_species(isp)%n_shells
                   lj = p_species(isp)%ang_mom(jshell)
                   do mj=-lj,lj
                      jproj = jproj + 1
                      if (mj==mi) dij0_sp(iproj,jproj,isp) = &
                           p_species(isp)%D0(ishell,jshell)
                   end do
                end do
             end do
          end do
       end if
    end do

    if (iscmplx) then
       allocate(dij0_sp_cmplx(max_atom_proj,max_atom_proj, &
            pub_cell%num_pspecies),stat=ierr)
       call utils_alloc_check('pseudo_get_dij','dij0_sp_cmplx',ierr)
       dij0_sp_cmplx(:,:,:) = cmplx(dij0_sp(:,:,:),kind=DP)
    end if

    ! Loop over atoms
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = loc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
       isp = pub_elements_on_node(loc_iat)%pspecies_number
       ! Put block of dij0 into SPAM3 matrix if there is one
       if (p_species(isp)%n_shells>0) then
          if (iscmplx) then
             call sparse_put_block(dij0_sp_cmplx(:,:,isp),dij,iat,iat)
          else
             call sparse_put_block(dij0_sp(:,:,isp),dij,iat,iat)
          end if
       end if
    end do

    if (iscmplx) then
       deallocate(dij0_sp_cmplx,stat=ierr)
       call utils_dealloc_check('pseudo_get_dij','dij0_sp_cmplx',ierr)
    end if

    deallocate(dij0_sp,stat=ierr)
    call utils_dealloc_check('pseudo_get_dij','dij0_sp',ierr)

  end subroutine pseudo_get_dij


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudopotentials_nonlocal_mat(nonlocpot, sp_overlap)

    !==================================================================!
    ! This subroutine calculates the sparse nonlocal potential matrix  !
    ! from the ngwf-projector overlap matrix.                          !
    !------------------------------------------------------------------!
    ! Original version written by Arash A. Mostofi on 12/01/2004.      !
    ! Modified by Chris-Kriton Skylaris on 25/1/2004 to work with      !
    ! the pseudo_species type array.                                   !
    ! Modified by Chris-Kriton Skylaris on 30/1/2004 to use the        !
    ! sparsity pattern of nonlocpot for determining which elements     !
    ! to calculate and store.                                          !
    ! Modified by Peter Haynes to use parallel SPAM 2, July 2006       !
    ! Modified by Nick Hine to use sparse_symmetric_expand to even out !
    ! load between processes                                           !
    ! Rewritten using SPAM3 matrices by Nicholas Hine, June 2009       !
    ! Modified by Nicholas Hine on 03/11/2009 to no longer need        !
    ! elements array, since pseudo_get_dij no longer does              !
    ! Modified by Laura Ratcliff to allow for calculation of a         !
    ! non-local matrix between 2 different sets of NGWFs, October 2010 !
    !==================================================================!

    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create,sparse_destroy,sparse_product, &
         sparse_transpose, sparse_transpose_structure
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: nonlocpot
    type(SPAM3), intent(in) :: sp_overlap

    ! Local Variables
    type(SPAM3) :: dij, sp_overlap_dij, ps_overlap

    call timer_clock('pseudopotentials_nonlocal_mat',1)

    ! Create the matrix structures (set as real or complex depending on
    ! whether sp_overlap is real or complex)
    call sparse_create(sp_overlap_dij,sp_overlap,iscmplx=sp_overlap%iscmplx)
    dij%structure = 'E'
    call sparse_create(dij,iscmplx=sp_overlap%iscmplx)

    ! ndmh: Create ps_overlap (modifying structure code if required)
    call sparse_transpose_structure(ps_overlap%structure,sp_overlap)
    call sparse_create(ps_overlap,iscmplx=sp_overlap%iscmplx)

    ! Create the matrix of nonlocal energies
    call pseudo_get_dij(dij)

    ! Calculate the matrix <NGWF_a|Proj_i> * D_ij
    call sparse_product(sp_overlap_dij,sp_overlap,dij)

    ! Create the transpose of the NGWF-projector overlap matrix
    call sparse_transpose(ps_overlap,sp_overlap)

    ! Calculate the matrix \sum_i (<NGWF_a|Proj_i> D_ij <Proj_j|NGWF_b> )
    call sparse_product(nonlocpot,sp_overlap_dij,ps_overlap)

    ! Clean up temporary matrices
    call sparse_destroy(ps_overlap)
    call sparse_destroy(dij)
    call sparse_destroy(sp_overlap_dij)

    call timer_clock('pseudopotentials_nonlocal_mat',2)

  end subroutine pseudopotentials_nonlocal_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pseudo_nonlocal_commutator_mat(nonlocpot_com, proj_basis, &
       ngwf_basis, ngwfs_on_grid, sp_overlap, delta_in)

    !==================================================================!
    ! This subroutine calculates the commutator between the nonlocal   !
    ! potential and the position operator for the 3 Cartesian          !
    ! directions.                                                      !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff February 2011.                         !
    ! Modified by Nicholas Hine in April 2011 to not use NGWF_REP and  !
    ! to simplify and cleanup existing code, and for changes to        !
    ! projector initialisation routines.                               !
    ! Modified by Laura Ratcliff Dec 2011 to fix bug involving complex !
    ! matrices and phase factors.                                      !
    !==================================================================!

    use comms, only: pub_on_root
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(*)
    use ion, only: ELEMENT
    use projectors, only: PROJECTOR_SET, projectors_func_ovlp_box, &
         projectors_deallocate_set, projectors_commutator
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product, &
         sparse_transpose, sparse_axpy, sparse_scale, sparse_copy, &
         sparse_transpose_structure
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: nonlocpot_com(3)
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(SPAM3), intent(in) :: sp_overlap
    real(kind=DP), intent(in) :: delta_in ! finite difference shift

    ! Local Variables
    type(SPAM3) :: dij

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering pseudo_nonlocal_commutator_mat'
#endif

    call timer_clock('pseudo_nonlocal_commutator_mat',1)

    dij%structure = 'E'
    call sparse_create(dij,iscmplx=.true.)

    ! Create the matrix of Kleinman-Bylander denominators
    call pseudo_get_dij(dij)
    call projectors_commutator(nonlocpot_com, proj_basis, &
         ngwf_basis, ngwfs_on_grid, sp_overlap, nlps_projectors, &
         delta_in, dij)
    call sparse_destroy(dij)

    call timer_clock('pseudo_nonlocal_commutator_mat',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving pseudo_nonlocal_commutator_mat'
#endif


  end subroutine pseudo_nonlocal_commutator_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pseudo_nonlocal_com_mat_fd(nonlocpot_com, &
       proj_basis, ngwf_basis, ngwfs_on_grid, sp_overlap, &
       nonlocpot, delta)

    !==================================================================!
    ! This subroutine calculates the commutator between the nonlocal   !
    ! potential and the position operator for the 3 Cartesian          !
    ! directions.                                                      !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff February 2011.                         !
    ! Modified by Nicholas Hine in April 2011 to not use NGWF_REP      !
    ! Modified by Laura Ratcliff Dec 2011 to fix bug involving complex !
    ! matrices and phase factors and changed from centred difference   !
    ! to forward difference.                                           !
    !==================================================================!

    use comms, only: pub_on_root
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(*)
    use ion, only: ELEMENT
    use projectors, only: projectors_func_ovlp_box, projectors_deallocate_set
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product, &
         sparse_transpose, sparse_axpy, sparse_scale, sparse_copy, &
         sparse_transpose_structure
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: nonlocpot_com(3)
    type(SPAM3), intent(in) :: nonlocpot
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(SPAM3), intent(in) :: sp_overlap
    real(kind=DP), intent(in) :: delta ! finite difference shift

    ! Local Variables
    type(SPAM3) :: dij, sp_overlap_dij, ps_overlap, dij_ps_overlap
    type(SPAM3) :: sDGp_overlap, DGps_overlap, sp_overlap_cmplx
    type(SPAM3) :: ps_overlap_rc, sDGp_overlap_rc
    type(SPAM3) :: nonlocpot_com2
    type(POINT) :: q_vec, q_cart(3)
    real(kind=DP) :: inv_delta
    integer :: cart

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering pseudo_nonlocal_com_mat_fd'
#endif

    call timer_clock('pseudo_nonlocal_com_mat_fd',1)

    ! NB not allowing for sp_overlap or ps_overlap to be at diff k-points
    inv_delta = 1.0_dp / delta

    q_vec%x = delta
    q_vec%y = delta
    q_vec%z = delta

    do cart=1,3
      q_cart(cart)%x = 0.0_dp
      q_cart(cart)%y = 0.0_dp
      q_cart(cart)%z = 0.0_dp
    end do

    q_cart(1)%x = q_vec%x
    q_cart(2)%y = q_vec%y
    q_cart(3)%z = q_vec%z

    call sparse_create(nonlocpot_com2,nonlocpot_com(1),iscmplx=.true.)

    call sparse_create(sp_overlap_dij,sp_overlap,iscmplx=.true.)
    call sparse_create(sDGp_overlap,sp_overlap,iscmplx=.true.)
    call sparse_create(sDGp_overlap_rc,sp_overlap,iscmplx=.false.)

    ! ndmh: Create ps_overlap
    call sparse_transpose_structure(ps_overlap%structure,sp_overlap)
    ps_overlap_rc%structure = ps_overlap%structure
    dij_ps_overlap%structure = ps_overlap%structure
    DGps_overlap%structure = ps_overlap%structure

    call sparse_create(dij_ps_overlap,iscmplx=.true.)
    call sparse_create(DGps_overlap,iscmplx=.true.)

    call sparse_create(sp_overlap_cmplx,sp_overlap,iscmplx=.true.)
    call sparse_create(ps_overlap,iscmplx=.true.)
    call sparse_create(ps_overlap_rc,iscmplx=.false.)

    ! Create the transpose of the NGWF-projector overlap matrix
    call sparse_transpose(ps_overlap_rc,sp_overlap)

    call sparse_copy(ps_overlap,ps_overlap_rc)
    call sparse_copy(sp_overlap_cmplx,sp_overlap)

    call sparse_destroy(ps_overlap_rc)

    dij%structure = 'E'
    call sparse_create(dij,iscmplx=.true.)

    ! Create the matrix of Kleinman-Bylander denominators
    call pseudo_get_dij(dij)

    call sparse_product(sp_overlap_dij,sp_overlap_cmplx,dij)
    call sparse_product(dij_ps_overlap,dij,ps_overlap)

    call sparse_destroy(dij)
    call sparse_destroy(ps_overlap)
    call sparse_destroy(sp_overlap_cmplx)

    ! Loop over Cartesian directions
    do cart=1,3

       q_vec = q_cart(cart)

       ! Calculate <phi|Del_G(proj(G))> overlap matrix for ngwf_basis
       ! real part
       call projectors_func_ovlp_box(sDGp_overlap_rc, &
            ngwfs_on_grid, ngwf_basis, proj_basis, nlps_projectors, q_vec, &
            swap_rc=.false.)

       ! ndmh: add real part to complex sp_overlap matrix
       call sparse_copy(sDGp_overlap,sDGp_overlap_rc)

       ! imaginary part
       call projectors_func_ovlp_box(sDGp_overlap_rc, &
            ngwfs_on_grid, ngwf_basis, proj_basis, nlps_projectors, q_vec, &
            swap_rc=.true.)

       ! ndmh: add imaginary part to complex sp_overlap matrix
       call sparse_axpy(sDGp_overlap,sDGp_overlap_rc,cmplx(0.0,1.0,dp))

       ! Transpose <phi|Del_G(proj(G))> to get <Del_G(proj(G))|phi>
       call sparse_transpose(DGps_overlap,sDGp_overlap)

       ! Calculate the matrix \sum_ij <phi_a|Del_G(p_i(G))>D_ij<p_j|phi_b>
       call sparse_product(nonlocpot_com2,sDGp_overlap,dij_ps_overlap)

       ! Calculate the matrix \sum_ij <phi_a|p_i>D_ij<Del_G(p_j(G))|phi_b>
       call sparse_product(nonlocpot_com(cart),sp_overlap_dij,DGps_overlap)

       ! Sum the two
       call sparse_axpy(nonlocpot_com(cart),nonlocpot_com2,cmplx(1.0,0.0,dp))

       ! Subtract 2 * nonlocpot
       call sparse_axpy(nonlocpot_com(cart),nonlocpot,cmplx(-2.0,0.0,dp))

       ! Multiply by -i/|q|
       call sparse_scale(nonlocpot_com(cart),cmplx(0.0,-inv_delta,dp))

    end do

    ! Clean up temporary matrices
    call sparse_destroy(DGps_overlap)
    call sparse_destroy(dij_ps_overlap)
    call sparse_destroy(sDGp_overlap_rc)
    call sparse_destroy(sDGp_overlap)
    call sparse_destroy(sp_overlap_dij)
    call sparse_destroy(nonlocpot_com2)

    call timer_clock('pseudo_nonlocal_com_mat_fd',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving pseudo_nonlocal_com_mat_fd'
#endif


  end subroutine pseudo_nonlocal_com_mat_fd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pseudo_nonlocal_com_mat_direct(nonlocpot_com, proj_basis, &
       ngwf_basis, ngwfs_on_grid, sp_overlap)

    !==================================================================!
    ! This subroutine calculates the commutator between the nonlocal   !
    ! potential and the position operator for the 3 Cartesian          !
    ! directions.                                                      !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff February 2011.                         !
    !==================================================================!

    use comms, only: pub_on_root
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(*)
    use ion, only: ELEMENT
    use projectors, only: projectors_func_ovlp_box, projectors_deallocate_set,&
         projectors_func_pos_ovlp_box
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product, &
         sparse_transpose, sparse_axpy, sparse_scale, sparse_copy, &
         sparse_transpose_structure
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: nonlocpot_com(3)
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(SPAM3), intent(in) :: sp_overlap

    ! Local Variables
    type(SPAM3) :: dij, sp_overlap_KB, ps_overlap, KB_ps_overlap
    type(SPAM3) :: sDGp_overlap(3), DGps_overlap
    type(SPAM3) :: nonlocpot_com2

    integer :: cart

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering pseudo_nonlocal_com_mat_direct'
#endif

    call timer_clock('pseudo_nonlocal_com_mat_direct',1)

    ! Assuming all matrices are real

    call sparse_create(nonlocpot_com2,nonlocpot_com(1),iscmplx=.false.)
    call sparse_create(sp_overlap_KB,sp_overlap)
    do cart=1,3
       call sparse_create(sDGp_overlap(cart),sp_overlap)
    end do

    ! ndmh: Create ps_overlap (modifying structure code if required)
    call sparse_transpose_structure(ps_overlap%structure,sp_overlap)
    KB_ps_overlap%structure = ps_overlap%structure
    DGps_overlap%structure = ps_overlap%structure

    call sparse_create(KB_ps_overlap)
    call sparse_create(DGps_overlap)
    call sparse_create(ps_overlap)

    ! Create the transpose of the NGWF-projector overlap matrix
    call sparse_transpose(ps_overlap,sp_overlap)

    dij%structure = 'E'
    call sparse_create(dij,iscmplx=.false.)

    ! Create the matrix of Kleinman-Bylander denominators
    call pseudo_get_dij(dij)

    call sparse_product(sp_overlap_KB,sp_overlap,dij)
    call sparse_product(KB_ps_overlap,dij,ps_overlap)

    call sparse_destroy(dij)
    call sparse_destroy(ps_overlap)

    ! Calculate <phi|r|proj> overlap matrix
    call projectors_func_pos_ovlp_box(sDGp_overlap, &
         ngwfs_on_grid,ngwf_basis,proj_basis,nlps_projectors)

    ! Loop over Cartesian directions
    do cart=1,3 

       ! Transpose <phi|r|proj> to get <proj|r|phi> for ngwf_basis
       call sparse_transpose(DGps_overlap,sDGp_overlap(cart))

       ! Calculate the matrix \sum_i (<NGWF_a|r|Proj_i><Proj_i|NGWF_b>/D_i)
       call sparse_product(nonlocpot_com2,sDGp_overlap(cart),KB_ps_overlap) 

       ! Calculate the matrix \sum_i (<NGWF_a|Proj_i><Proj_i|r|NGWF_b>/D_i)
       call sparse_product(nonlocpot_com(cart),sp_overlap_KB,DGps_overlap) 

       ! Sum the two
       call sparse_axpy(nonlocpot_com(cart),nonlocpot_com2,-1.0_dp)

    end do

    ! Clean up temporary matrices
    call sparse_destroy(DGps_overlap)
    call sparse_destroy(KB_ps_overlap)
    do cart=3,1,-1
       call sparse_destroy(sDGp_overlap(cart))
    end do
    call sparse_destroy(sp_overlap_KB)
    call sparse_destroy(nonlocpot_com2)

    call timer_clock('pseudo_nonlocal_com_mat_direct',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving pseudo_nonlocal_com_mat_direct'
#endif

  end subroutine pseudo_nonlocal_com_mat_direct


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_aug_Q_matrix(aug_Q)

    !==================================================================!
    ! This subroutine returns the Q matrix for all the USP atoms in    !
    ! the simulation cell.                                             !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  proj_overlap (inout) : The SPAM3 block-diagonal overlap matrix  !
    !  of the partial waves of each atom                               !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 11/02/11.           !
    !==================================================================!

    use comms, only: pub_my_node_id
    use parallel_strategy, only: pub_elements_on_node, &
         pub_num_atoms_on_node, pub_first_atom_on_node
    use sparse, only: SPAM3, sparse_put_block
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: aug_Q

    ! Local Variables
    integer :: iat
    integer :: loc_iat
    integer :: isp
    integer :: ishell,jshell
    integer :: iproj_at,jproj_at
    integer :: li,lj,mi,mj
    integer :: ierr
    integer :: max_proj
    real(kind=DP),allocatable :: aug_Q_block(:,:)

    max_proj = maxval(nlps_projectors%species_num_proj(:))
    allocate(aug_Q_block(max_proj,max_proj),stat=ierr)
    call utils_alloc_check('pseudo_aug_Q_matrix','ovlp_block',ierr)

    ! Loop over atoms on this node
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1
       isp = pub_elements_on_node(loc_iat)%pspecies_number

       ! Double loop over partial waves i,j
       iproj_at = 0
       do ishell=1,p_species(isp)%n_shells
          li = p_species(isp)%ang_mom(ishell)
          do mi=-li,li
             iproj_at = iproj_at + 1
             jproj_at = 0
             do jshell=1,p_species(isp)%n_shells
                lj = p_species(isp)%ang_mom(jshell)
                do mj=-lj,lj
                   jproj_at = jproj_at + 1

                   ! Find Q_ij for this projector
                   aug_Q_block(iproj_at,jproj_at) = &
                        p_species(isp)%aug_q(ishell,jshell)
                enddo
             end do
          end do
       end do

       ! Put this atom's block into SPAM3 matrix
       call sparse_put_block(aug_Q_block,aug_Q,iat,iat)
    end do

    deallocate(aug_Q_block,stat=ierr)
    call utils_dealloc_check('pseudo_aug_Q_matrix','ovlp_block',ierr)

  end subroutine pseudo_aug_Q_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_atom_aug_den(atom_nhat,atom_grad_nhat,atom_aug_func, &
       atom_aug_func_grad,total_nhat,total_nhat_targ,rho_ij_block, &
       isp,atom_centre,grid,box_n1,box_n2,box_n3, &
       box_start1,box_start2,box_start3)

    use cell_grid, only: GRID_INFO
    use constants, only: max_spins
    use geometry, only: POINT
    use simulation_cell, only: pub_cell
    use xc, only: pub_xc_gradient_corrected

    ! Arguments
    integer, intent(in) :: box_n1,box_n2,box_n3
    real(kind=DP), intent(out) :: atom_aug_func(box_n1,box_n2,box_n3)
    real(kind=DP), intent(out) :: atom_aug_func_grad(box_n1,box_n2,box_n3,3)
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

    ! Local Variables
    integer :: ishell,jshell
    integer :: iproj,jproj
    integer :: li,mi,lj,mj
    integer :: is
    integer :: cart
    real(kind=DP) :: Qij(max_spins)
    real(kind=DP) :: total_aug_func

    ! Loop over projectors i,j
    iproj = 0
    total_nhat_targ = 0.0_DP
    total_nhat = 0.0_DP
    do ishell=1,p_species(isp)%n_shells
       li = p_species(isp)%ang_mom(ishell)
       do mi=-li,li
          iproj = iproj + 1
          jproj = 0
          do jshell=1,p_species(isp)%n_shells
             lj = p_species(isp)%ang_mom(jshell)
             do mj=-lj,lj
                jproj = jproj + 1

                Qij(:) = 0.0_DP

                ! Find contribution to total target augmentation charge
                ! for this i,j pair
                do is=1,pub_cell%num_spins
                   Qij(is) = Qij(is) + rho_ij_block(iproj,jproj,is)
                end do

                if (any(abs(Qij(:))>1e-20_DP)) then

                   ! Get shape function g_L(r) multiplied by spherical harmonic
                   ! S_LM(r) for this L,M pair
                   atom_aug_func = 0.0_DP
                   call pseudo_get_aug_func(atom_aug_func, &
                        isp,ishell,jshell,mi,mj, &
                        box_start1,box_start2,box_start3, &
                        grid,box_n1,box_n2,box_n3,atom_centre)

                   total_aug_func = sum(atom_aug_func(:,:,:))*grid%weight

                   do is=1,pub_cell%num_spins
                      atom_nhat(:,:,:,is) = atom_nhat(:,:,:,is) &
                           + Qij(is)*atom_aug_func(:,:,:)
                      total_nhat(is) = total_nhat(is) + total_aug_func*Qij(is)
                   end do

                   ! Find gradient of \hat{n}(r) if required
                   if (pub_xc_gradient_corrected) then
                      ! Get gradient of (shape function g_L(r) multiplied by
                      ! spherical harmonic S_LM(r) for this L,M pair)
                      atom_aug_func_grad = 0.0_DP
                      call pseudo_get_aug_func(atom_aug_func_grad, &
                           isp,ishell,jshell,mi,mj, &
                           box_start1,box_start2,box_start3, &
                           grid,box_n1,box_n2,box_n3,atom_centre)
                      do cart=1,3
                         do is=1,pub_cell%num_spins
                            atom_grad_nhat(:,:,:,is,cart) = &
                                 atom_grad_nhat(:,:,:,is,cart) + &
                                 Qij(is)*atom_aug_func_grad(:,:,:,cart)
                         end do
                      end do
                   end if
                end if

             end do
          end do
       end do
    end do

  end subroutine pseudo_atom_aug_den


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_atom_aug_integrals(locpot_box,atom_aug_func, &
       dij_at,num_spins,isp,atom_centre,grid,box_n1,box_n2,box_n3, &
       box_start1,box_start2,box_start3)

    use cell_grid, only: GRID_INFO
    use constants, only: max_spins
    use geometry, only: POINT

    implicit none

    ! Arguments
    integer, intent(in) :: box_n1,box_n2,box_n3
    integer, intent(in) :: num_spins
    real(kind=DP), intent(out) :: atom_aug_func(box_n1,box_n2,box_n3)
    real(kind=DP), intent(in) :: locpot_box(box_n1,box_n2,box_n3,num_spins)
    type(GRID_INFO), intent(in) :: grid
    type(POINT), intent(in) :: atom_centre
    real(kind=DP), intent(inout) :: dij_at(:,:,:)
    integer, intent(in) :: isp
    integer, intent(in) :: box_start1,box_start2,box_start3

    ! Local Variables
    integer :: iproj, jproj
    integer :: ishell, jshell
    integer :: li, lj
    integer :: mi, mj
    integer :: is
    integer :: locpot_Qij_product(max_spins)

    ! Loop over projectors i,j
    iproj = 0
    do ishell=1,p_species(isp)%n_shells
       li = p_species(isp)%ang_mom(ishell)
       do mi=-li,li
          iproj = iproj + 1
          jproj = 0
          do jshell=1,p_species(isp)%n_shells
             lj = p_species(isp)%ang_mom(jshell)
             do mj=-lj,lj
                jproj = jproj + 1

                ! Get the augmentation function Qij(r) corresponding
                ! to this pair of shells i,j and this mi,mj
                atom_aug_func = 0.0_DP
                call pseudo_get_aug_func(atom_aug_func, &
                     isp,ishell,jshell,mi,mj, &
                     box_start1,box_start2,box_start3, &
                     grid,box_n1,box_n2,box_n3,atom_centre)

                ! Calculate integral \int veff(r) Qij(r) dr
                do is=1,num_spins
                   locpot_Qij_product(is) = sum(locpot_box(:,:,:,is) &
                        * atom_aug_func(:,:,:)) * grid%weight
                end do

                dij_at(iproj,jproj,is) = dij_at(iproj,jproj,is) + &
                     locpot_Qij_product(is)

             end do
          end do
       end do
    end do

  end subroutine pseudo_atom_aug_integrals


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_atom_aug_force(nhat_force,locpot_box,atom_grad_aug_func, &
       rhoij_at,num_spins,isp,atom_centre,grid,box_n1,box_n2,box_n3, &
       box_start1,box_start2,box_start3)

    use cell_grid, only: GRID_INFO
    use constants, only: max_spins
    use geometry, only: POINT
    use simulation_cell, only: pub_cell

    ! Arguments
    integer, intent(in) :: box_n1,box_n2,box_n3
    integer, intent(in) :: num_spins
    real(kind=DP), intent(out) :: atom_grad_aug_func(box_n1,box_n2,box_n3,3)
    real(kind=DP),intent(inout) :: nhat_force(3)
    real(kind=DP), intent(in) :: locpot_box(box_n1,box_n2,box_n3,num_spins)
    type(GRID_INFO), intent(in) :: grid
    type(POINT), intent(in) :: atom_centre
    real(kind=DP), intent(in) :: rhoij_at(:,:,:)
    integer, intent(in) :: isp
    integer, intent(in) :: box_start1,box_start2,box_start3

    ! Local Variables
    integer :: ishell,jshell
    integer :: iproj,jproj
    integer :: li,mi,lj,mj
    integer :: is
    integer :: cart
    real(kind=DP) :: locpot_Qij_grad_product(3,max_spins)

    ! Loop over projectors i,j
    iproj = 0
    do ishell=1,p_species(isp)%n_shells
       li = p_species(isp)%ang_mom(ishell)
       do mi=-li,li
          iproj = iproj + 1
          jproj = 0
          do jshell=1,p_species(isp)%n_shells
             lj = p_species(isp)%ang_mom(jshell)
             do mj=-lj,lj
                jproj = jproj + 1

                ! Construct gradient d/dR_I (Qij(r)) of augmentation
                ! function for this i,j pair
                atom_grad_aug_func = 0.0_DP
                call pseudo_get_aug_func_grad(atom_grad_aug_func, &
                     isp,ishell,jshell,mi,mj,box_start1,box_start2, &
                     box_start3,grid,box_n1,box_n2,box_n3, &
                     atom_centre)

                ! Loop over spins and loop over cartesian directions
                do is=1,num_spins
                   do cart=1,3
                      ! Calculate integral \int veff(r) d/dR(Qij(r)) dr
                      locpot_Qij_grad_product(cart,is) = &
                           sum(locpot_box(:,:,:,is) &
                           * atom_grad_aug_func(:,:,:,cart)) * grid%weight

                   end do
                end do

                do is=1,num_spins
                   nhat_force(:) = nhat_force(:) &
                        + (locpot_Qij_grad_product(:,is) &
                        * rhoij_at(iproj,jproj,is))
                end do  ! is

             end do
          end do
       end do
    end do

  end subroutine pseudo_atom_aug_force


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_aug_func(aug_func,isp,ishell,jshell,mi,mj, &
       box_start1,box_start2,box_start3,grid,box_n1,box_n2,box_n3,atom_origin)

    !==================================================================!
    ! This subroutine calculates the function g_L(r)*S_LM(\hat{r})     !
    ! which is required to form the compensation charge \hat{n}.       !
    ! It does so on the simulation cell fine grid in a box centered on !
    ! the atom.                                                        !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  box_start1 (in) :                                               !
    !  box_start2 (in) :                                               !
    !  box_start3 (in) :                                               !
    !  grid (in) :                                                     !
    !  box_n1 (in) :                                                   !
    !  box_n2 (in) :                                                   !
    !  box_n3 (in) :                                                   !
    !  atom_origin (in) :                                              !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 15/02/11.           !
    ! NOT FINISHED YET - NON-FUNCTIONAL!                               !
    !==================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: PI
    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use simulation_cell, only: pub_cell
    use spherical_wave, only: sw_real_sph_harm

    implicit none

    ! Arguments
    integer,intent(in) :: box_n1
    integer,intent(in) :: box_n2
    integer,intent(in) :: box_n3
    real(kind=DP),intent(out) :: aug_func(box_n1,box_n2,box_n3)
    integer, intent(in) :: isp
    integer,intent(in) :: ishell
    integer,intent(in) :: jshell
    integer,intent(in) :: mi,mj
    integer,intent(in) :: box_start1
    integer,intent(in) :: box_start2
    integer,intent(in) :: box_start3
    type(GRID_INFO),intent(in) :: grid
    type(POINT),intent(in) :: atom_origin

    ! Local Variables
    integer :: li, lj
    integer :: i1,i2,i3
    type(POINT) :: box_origin
    type(POINT) :: box_to_atom
    type(POINT) :: r_cell
    type(POINT) :: r_sphere
    real(kind=DP) :: rmag
    real(kind=DP) :: slmval

    ! Find vector to origin of box
    box_origin = (box_start1-1)*grid%da1 + &
                 (box_start2-1)*grid%da2 + &
                 (box_start3-1)*grid%da3

    ! ndmh: vector from origin of box to centre of atom
    box_to_atom = atom_origin - box_origin

    li = p_species(isp)%ang_mom(ishell)
    lj = p_species(isp)%ang_mom(jshell)
    if (mj==mi) i3=1

    ! ndmh: check components of this vector in terms of lattice vectors
    ! ndmh: are all positive. If not, the box has been looped back into the
    ! ndmh: cell, so shift the origin back accordingly
    if ((box_to_atom.DOT.pub_cell%b1) < 0.0_DP) then
       box_origin = box_origin - pub_cell%a1
    end if
    if ((box_to_atom.DOT.pub_cell%b2) < 0.0_DP) then
       box_origin = box_origin - pub_cell%a2
    end if
    if ((box_to_atom.DOT.pub_cell%b3) < 0.0_DP) then
       box_origin = box_origin - pub_cell%a3
    end if

    do i3=1,box_n3
       do i2=1,box_n2
          do i1=1,box_n1

             r_cell = box_origin + (i1-1)*grid%da1 &
                                 + (i2-1)*grid%da2 &
                                 + (i3-1)*grid%da3

             r_sphere = r_cell - atom_origin
             rmag = geometry_magnitude(r_sphere)

             if (rmag>p_species(isp)%core_radius(ishell)) then
                aug_func(i1,i2,i3) = 0.0_DP
                cycle
             end if

             ! Not functional yet!

             slmval = sw_real_sph_harm(r_sphere%x, &
                  r_sphere%y,r_sphere%z,rmag,li,mi)

             aug_func(i1,i2,i3) = slmval

          end do
       end do
    end do

  end subroutine pseudo_get_aug_func


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_get_aug_func_grad(aug_func_grad,isp,ishell,jshell,mi,mj, &
       box_start1,box_start2,box_start3,grid,box_n1,box_n2,box_n3,atom_origin)

    !==================================================================!
    ! This subroutine calculates the function g_L(r)*S_LM(\hat{r})     !
    ! which is required to form the compensation charge \hat{n}.       !
    ! It does so on the simulation cell fine grid in a box centered on !
    ! the atom.                                                        !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  box_start1 (in) :                                               !
    !  box_start2 (in) :                                               !
    !  box_start3 (in) :                                               !
    !  grid (in) :                                                     !
    !  box_n1 (in) :                                                   !
    !  box_n2 (in) :                                                   !
    !  box_n3 (in) :                                                   !
    !  atom_origin (in) :                                              !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 15/02/11.           !
    ! NOT FINISHED YET - NON-FUNCTIONAL!                               !
    !==================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: PI
    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use simulation_cell, only: pub_cell
    use spherical_wave, only: sw_real_sph_harm, sw_grad_real_sph_harm, &
         sw_init

    implicit none

    ! Arguments
    integer,intent(in) :: box_n1
    integer,intent(in) :: box_n2
    integer,intent(in) :: box_n3
    real(kind=DP),intent(out) :: aug_func_grad(box_n1,box_n2,box_n3,3)
    integer, intent(in) :: isp
    integer,intent(in) :: ishell
    integer,intent(in) :: jshell
    integer,intent(in) :: mi,mj
    integer,intent(in) :: box_start1
    integer,intent(in) :: box_start2
    integer,intent(in) :: box_start3
    type(GRID_INFO),intent(in) :: grid
    type(POINT),intent(in) :: atom_origin

    ! Local Variables
    integer :: li, lj
    integer :: lup, mup
    integer :: i1,i2,i3
    type(POINT) :: box_origin
    type(POINT) :: box_to_atom
    type(POINT) :: r_cell
    type(POINT) :: r_sphere
    type(POINT) :: r_sphere_unit
    real(kind=DP) :: rmag
    real(kind=DP) :: slmval
    real(kind=DP) :: grad_SLM_R(3)

    ! Find vector to origin of box
    box_origin = (box_start1-1)*grid%da1 + &
                 (box_start2-1)*grid%da2 + &
                 (box_start3-1)*grid%da3

    ! ndmh: vector from origin of box to centre of atom
    box_to_atom = atom_origin - box_origin

    li = p_species(isp)%ang_mom(ishell)
    lj = p_species(isp)%ang_mom(jshell)
    if (mj==mi) i3=1

    ! ndmh: check components of this vector in terms of lattice vectors
    ! ndmh: are all positive. If not, the box has been looped back into the
    ! ndmh: cell, so shift the origin back accordingly
    if ((box_to_atom.DOT.pub_cell%b1) < 0.0_DP) then
       box_origin = box_origin - pub_cell%a1
    end if
    if ((box_to_atom.DOT.pub_cell%b2) < 0.0_DP) then
       box_origin = box_origin - pub_cell%a2
    end if
    if ((box_to_atom.DOT.pub_cell%b3) < 0.0_DP) then
       box_origin = box_origin - pub_cell%a3
    end if

    do i3=1,box_n3
       do i2=1,box_n2
          do i1=1,box_n1

             r_cell = box_origin + (i1-1)*grid%da1 &
                                 + (i2-1)*grid%da2 &
                                 + (i3-1)*grid%da3

             r_sphere = r_cell - atom_origin
             rmag = geometry_magnitude(r_sphere)

             if (rmag>1d-10) then
                r_sphere_unit = (1.0_DP/rmag)*r_sphere
             else
                r_sphere_unit = 0.0_DP*r_sphere
             end if

             if (rmag>p_species(isp)%core_radius(ishell)) then
                aug_func_grad(i1,i2,i3,:) = 0.0_DP
                cycle
             end if

             ! Not functional yet!

             ! Get spherical harmonic at this \hat{r}
             slmval = sw_real_sph_harm(r_sphere%x,r_sphere%y,r_sphere%z, &
                  rmag,lup,mup)

             ! Get gradients of spherical harmonic at this \hat{r}
             call sw_grad_real_sph_harm(grad_SLM_R(:),r_sphere%x,r_sphere%y, &
                  r_sphere%z,rmag,lup,mup)

             ! S_LM(\hat{r}) * d/dR(g(|r-R|) part
             aug_func_grad(i1,i2,i3,1) = 1.0_DP * slmval * r_sphere_unit%x
             aug_func_grad(i1,i2,i3,2) = 1.0_DP * slmval * r_sphere_unit%y
             aug_func_grad(i1,i2,i3,3) = 1.0_DP * slmval * r_sphere_unit%z

             ! g(|r-R|) * d S_LM(\hat{r}) / dR(g(|r-R|) part
             aug_func_grad(i1,i2,i3,:) = aug_func_grad(i1,i2,i3,:) + &
                  grad_SLM_R(:) * 1.0_DP

          end do
       end do
    end do

  end subroutine pseudo_get_aug_func_grad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine pseudo_nl_calculate_forces(nlps_forces, &
       sp_overlap,ngwfs_on_grid,ngwf_basis,proj_basis, &
       pur_denskern)

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the nonlocal part of the ionic pseudopotential.                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1)  nlps_forces     : output : nonlocal forces                          !
    ! 2)  sp_overlap      : input  : ngwf-projector overlap matrix            !
    ! 3)  ngwfs_on_grid   : input  : NGWF data in ppd format                  !
    ! 4)  ngwf_basis      : input  : Function basis description for NGWFs     !
    ! 5)  proj_basis      : input  : Function basis description for projectors!
    ! 6)  pur_denskern    : input  : purified density kernel SPAM3            !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !                                                                         !
    ! 1) basis                                                                !
    ! 2) comms                                                                !
    ! 3) constants                                                            !
    ! 4) ion                                                                  !
    ! 5) simulation_cell                                                      !
    ! 6) parallel_strategy                                                    !
    ! 7) timer                                                                !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Written by Arash A. Mostofi, V0.0, 28th June 2004                       !
    !-------------------------------------------------------------------------!
    ! Modified by Arash A. Mostofi, 15th September 2004                       !
    !   - Distribution of the derivative with respect to atomic positions of  !
    !     the non-local potential matrix (nlfmtx_node) amongst nodes          !
    ! Modified by Peter Haynes for parallel SPAM3, July 2006                  !
    !   - Now all nodes calculate contributions to the forces on all atoms    !
    !     and nlfmtx_node is eliminated                                       !
    ! Modified by Nicholas Hine to improve efficiency of inner loop,          !
    !     November 2008.                                                      !
    ! Modified by Nicholas Hine in July 2009 to work with SPAM3.              !
    !=========================================================================!

    use comms, only: comms_barrier, comms_reduce, pub_on_root, pub_my_node_id
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: pub_first_atom_on_node, pub_orig_atom
    use projectors, only: projectors_func_grad_ovlp_box
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_axpy, sparse_create, sparse_destroy, &
         sparse_get_element, sparse_product, sparse_transpose, sparse_scale, &
         sparse_transpose_structure
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: nlps_forces(1:3,pub_cell%nat)
    type(SPAM3), intent(in) :: sp_overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(SPAM3), intent(in) :: pur_denskern(pub_cell%num_spins)

    ! Local Variables
    type(SPAM3) :: siGp_overlap,iGps_overlap
    type(SPAM3) :: sp_overlap_dij,kb_ps_overlap
    type(SPAM3) :: kb_denom
    type(SPAM3) :: nl_force_mat(3)
    type(SPAM3) :: kq,rkq
    integer :: cart, is
    integer :: iat, orig_iat
    integer :: atom_proj, global_proj
    real(kind=DP) :: proj_force

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering pseudo_nl_calculate_forces'
#endif

    ! Start timer
    call timer_clock('pseudo_nl_calculate_forces',1)

    ! Initialise
    nlps_forces  = 0.0_DP

    ! Create result matrices to hold nonlocal forces, projector-by-projector
    do cart=1,3
       nl_force_mat(cart)%structure = 'E'
       call sparse_create(nl_force_mat(cart))
    end do

    ! Create matrices to hold <NGWF_a|Proj_i> / D_i and <Proj_i|NGWF_a> / D_i
    call sparse_transpose_structure(kb_ps_overlap%structure,sp_overlap)
    call sparse_create(kb_ps_overlap)
    call sparse_create(sp_overlap_dij,sp_overlap)

    ! Get 1/KB denominators in diagonal matrix
    kb_denom%structure = 'E'
    call sparse_create(kb_denom)
    call pseudo_get_dij(kb_denom)

    ! Calculate the matrices <NGWF_a|Proj_i> / D_i and <Proj_i|NGWF_a> / D_i
    call sparse_product(sp_overlap_dij,sp_overlap,kb_denom)
    call sparse_transpose(kb_ps_overlap,sp_overlap_dij)

    ! Destroy kb_denom matrix to free up memory
    call sparse_destroy(kb_denom)

    ! Create matrices to hold <phi|iG*proj> overlap matrix and transpose
    call sparse_create(siGp_overlap,sp_overlap)
    call sparse_create(iGps_overlap,kb_ps_overlap)

    ! Create temporary matrices kq and rkq
    call sparse_create(kq,pur_denskern(1),sp_overlap_dij)
    call sparse_create(rkq,nl_force_mat(1))

    ! Loop over Cartesian directions
    do cart=1,3

       ! Calculate <phi|iG*proj> overlap matrix
       call projectors_func_grad_ovlp_box(siGp_overlap, &
            ngwfs_on_grid,ngwf_basis,proj_basis,nlps_projectors,cart)

       ! Transpose it to get <iG*proj|phi> overlap matrix
       call sparse_transpose(iGps_overlap,siGp_overlap)

       do is=1,pub_cell%num_spins

          ! Calculate <proj_j|phi_b> K^ab <phi_a|iG.proj_i>
          call sparse_product(kq,pur_denskern(is),siGp_overlap)
          call sparse_product(rkq,kb_ps_overlap,kq)
          call sparse_axpy(nl_force_mat(cart),rkq,1.0_DP)

          ! Calculate <iG.proj_j|phi_b> K^ab <phi_a|proj_i>
          call sparse_product(kq,pur_denskern(is),sp_overlap_dij)
          call sparse_product(rkq,iGps_overlap,kq)
          call sparse_axpy(nl_force_mat(cart),rkq,1.0_DP)

       end do
    end do

    ! Destroy temporary matrices
    call sparse_destroy(rkq)
    call sparse_destroy(kq)
    call sparse_destroy(iGps_overlap)
    call sparse_destroy(siGp_overlap)
    call sparse_destroy(kb_ps_overlap)
    call sparse_destroy(sp_overlap_dij)

    ! Loop over atoms
    do iat=pub_first_atom_on_node(pub_my_node_id), &
         pub_first_atom_on_node(pub_my_node_id + 1) - 1

       ! Find atom number in input file order
       orig_iat = pub_orig_atom(iat)

       ! Loop over projectors on this atom
       do atom_proj=1,proj_basis%num_on_atom(iat)
          global_proj = proj_basis%first_on_atom(iat) + atom_proj - 1

          ! Loop over Cartesian co-ordinates
          do cart=1,3

             ! Find contribution of this projector to force on this atom
             ! from diagonal elements of nl_force_mat for this coordinate
             call sparse_get_element(proj_force,nl_force_mat(cart), &
                  global_proj, global_proj)
             nlps_forces(cart,orig_iat) = nlps_forces(cart,orig_iat) + proj_force

          end do  ! cart

       end do  ! atom_proj

    end do   ! loc_iat

    ! Reduce result across nodes
    call comms_barrier
    call comms_reduce('SUM',nlps_forces,3*pub_cell%nat)

    ! Destroy temporary matrices
    do cart=1,3
       call sparse_destroy(nl_force_mat(cart))
    end do

    ! Stop timer
    call timer_clock('pseudo_nl_calculate_forces',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving pseudo_nl_calculate_forces'
#endif

  end subroutine pseudo_nl_calculate_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pseudo_local_on_grid_openbc(v_local_fine, grid, elements)
    !========================================================================!
    ! Calculates the local part of the pseudopotential on the fine grid,     !
    ! without resorting to discrete Fourier transforms, and thus realizing   !
    ! OPEN boundary conditions.                                              !
    ! This is achieved by calculating a continuous Fourier transform, to     !
    ! obtain Vloc(x) from Vloc(g), both on radial grids. Interpolation is    !
    ! then used to obtain Vloc(r-R_i) for all fine grid points r, coming from!
    ! ions at R_i.                                                           !
    ! The integral is tricky, because:                                       !
    ! 1) for large values of x, the oscillations in sin(gx) become very      !
    !    frequent, requiring finer resolution than that of the recpot file,  !
    ! 2) the low-g part of the integral is singular,                         !
    ! 3) it must be computed to very, very high accuracy to get energies with!
    !    an accuracy of 0.001%.                                              !
    !------------------------------------------------------------------------!
    ! Arguments                                                              !
    ! v_local_fine (output): Contains v_loc on the fine grid on output.      !
    ! elements (input): The elements array req'd to perform the calculation. !
    !------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in May 2010, based on notes by Chris-Kriton  !
    ! Skylaris. Updated by Jacek Dziedzic in Sep 2010 to use the trick to    !
    ! split the potential into lr and sr parts in order to improve accuracy. !
    !========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_my_node_id, pub_on_root, pub_total_num_nodes, &
         comms_allgather
    use geometry, only: magnitude, POINT, operator(+), operator(-), operator(*)
    use ion, only: ELEMENT
    use rundat, only: pub_open_localpseudo, &
         pub_openbc_pspot_finetune_f, pub_openbc_pspot_finetune_nptsx, &
         pub_openbc_pspot_finetune_alpha
    use services, only: services_1d_interpolation
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_assert, utils_abort, &
         utils_alloc_check, utils_dealloc_check, utils_flush, utils_unit, &
         utils_erf

    implicit none

    ! jd: Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(out)  :: v_local_fine(grid%ld1, &
         grid%ld2, grid%max_slabs12)
    type(ELEMENT), intent(in)   :: elements(pub_cell%nat)

    ! jd: Local variables, parameters
    real(kind=DP), allocatable :: vloc_lookup(:,:)
    ! jd: Vloc(x) on radial grid, for all species

    real(kind=DP), allocatable :: vloc_lookup_local(:)
    ! jd: Portion of the above local to this core, for current species

    integer :: npts_x                    ! jd: Size of the radial realspace grid
    integer :: npts_x_per_core           ! jd: One core's share of the abovee
    integer :: species, x_rad_pt         ! jd: Indices
    integer :: islab12, i1, i2, i3, I    ! jd: Indices
    real(kind=DP) :: hx                  ! jd: Gridpoint separation
    real(kind=DP) :: xmax                ! jd: Last x stored in lookups
    real(kind=DP) :: x                   ! jd: Current x
    type(POINT)   :: R_I                 ! jd: Position of ion I
    type(POINT)   :: r                   ! jd: Position of curr. fine gridpoint
    real(kind=DP) :: accum               ! jd: Accumulator
    real(kind=DP) :: delta_g             ! jd: g spacing in the recpot file
    real(kind=DP) :: Z                   ! jd: Charge of current ion
    real(kind=DP) :: Vloc_here           ! jd: Interpolated vloc at curr. gridpt
    integer       :: x1, x2              ! jd: min and max x for this core
    real(kind=DP) :: coulombic           ! jd: Needed for self-test
    real(kind=DP) :: reldiff, absreldiff ! jd: Needed for self-test
    real(kind=DP) :: avgabsreldiff       ! jd: Needed for self-test
    real(kind=DP) :: avgreldiff          ! jd: Needed for self-test
    integer       :: avgpts              ! jd: Needed for self-test
    real(kind=DP) :: maxabsreldiff       ! jd: Needed for self-test
    integer       :: fineness            ! jd: fineness of the radial g-grid
    integer       :: ierr                ! jd: Error flag
    real(kind=DP) :: alpha               ! jd: Crossover between sr and lr term
    real(kind=DP) :: erf_alpha_x_over_x  ! jd: erf(alpha*x)/x
    real(kind=DP), parameter :: factor = 2.0_DP / PI
    ! jd: 4pi (from ang. integration) * 4pi (onetep real/recip space convention)
    !     * 1/(2pi)**3 (Fourier transform factor) = 2/pi
#ifdef DEBUG
    integer :: output_unit = 0           ! jd: Dummy init. to kill compiler warn
#endif

    ! -----------------------------------------------------------------------

    call timer_clock('pseudo_local_on_fine_openbc',1)

    if(pub_on_root) write(stdout,'(a)') &
         'Calculating local pseudopotential in real space ...'
    call utils_flush

    ! jd: Set up parameters
    npts_x = pub_openbc_pspot_finetune_nptsx
    alpha = pub_openbc_pspot_finetune_alpha / min( &
         magnitude(pub_cell%a1),magnitude(pub_cell%a2),magnitude(pub_cell%a3))

    ! jd: Sanity check
    call utils_assert(pub_open_localpseudo,'Cannot use &
         &pseudo_local_on_fine_openbc without having open boundary conditions')

    ! jd: Find out the maximum value of x which we might need, determine the
    !     fineness of the radial grid  depending on this
    xmax = magnitude(pub_cell%a1 + pub_cell%a2 + pub_cell%a3)
    hx = xmax/(npts_x-1)

    ! jd: Split the work across cores: every core will work with points x1..x2
    npts_x_per_core = npts_x / pub_total_num_nodes
    if(mod(npts_x,pub_total_num_nodes) /= 0) then
       npts_x_per_core = npts_x_per_core + 1
    end if
    x1 = pub_my_node_id*npts_x_per_core + 1
    x2 = (pub_my_node_id+1)*npts_x_per_core

    ! jd: Tail on last node may be shorter, adjust x2
    !     Leave npts_x_per_core alone, because it must be the same on all cores
    !     for Allgather to work
    if(x2 > npts_x .and. pub_my_node_id == pub_total_num_nodes-1) x2 = npts_x

    ! jd: Sanity check against a situation where the last core(s) is/are left
    !     with no work
    if(x1 > npts_x .or. x2 > npts_x) then
       call utils_abort('Too many cores to sensibly divide &
            &openbc_pspot_finetune_nptsx radial points in the open BC &
            &pseudopotential integration. Increase the number of points or &
            &reduce the number of cores.')
    end if

    ! jd: Allocate and zero lookup array. It must be large enough to hold
    !     not only npts_x, but (ncores * npts_x_per_core), which will be
    !     somewhat larger than npts_x.
    allocate(vloc_lookup(npts_x_per_core * pub_total_num_nodes, &
         pub_cell%num_pspecies),stat=ierr)
    call utils_alloc_check('pseudo_local_on_fine_openbc','vloc_lookup',ierr)
    vloc_lookup = 0.0_DP

    ! jd: Allocate the part of the lookup local to this core
    allocate(vloc_lookup_local(npts_x_per_core),stat=ierr)
    call utils_alloc_check('pseudo_local_on_fine_openbc','vloc_lookup_local', &
         ierr)
    vloc_lookup_local = 0.0_DP

    if (pub_on_root) then
       write(stdout,'(a,f0.3,a)') 'Discrepancy between calculated realspace &
            &psloc and -Z/x on [5a0,', (npts_x-1)*hx, 'a0]:'
       write(stdout,'(a)') 'max |rel. diff| # avg |rel. diff| # &
            & avg rel. diff #  f # pseudopotential file'
    end if

    ! jd: Fill up the lookup arrays, each core deals with its own portion of x's
    do species = 1, pub_cell%num_pspecies

       Z = p_species(species)%ion_charge
       delta_g = 1.0_DP / p_species(species)%inv_g_spacing

       ! jd: If the user did not choose an explicit value, select a sensible
       !     default for the fineness parameter. It is chosen such, that we
       !     still have at least 50 points to sample a full period of sin(gx)
       !     for the largest possible x, on the g-grid.
       if(pub_openbc_pspot_finetune_f == -1) then
          fineness = 50.0_DP * delta_g / (2.0_DP*PI/xmax) + 1
       else
          fineness = pub_openbc_pspot_finetune_f
       end if

       ! jd: Go over all of the x's that belong to this core
       do x_rad_pt = x1, x2
          x = (x_rad_pt-1) * hx

          if(x>1D-10) then
            erf_alpha_x_over_x = utils_erf(alpha*x)/x
          else
            ! jd: lim x->0 (erf(alpha*x)/x)
            erf_alpha_x_over_x = 2.0_DP*alpha/sqrt(PI)
          end if

          ! jd: Vloc(x) is separated into two parts: short-range and long-range,
          !     the short-range part is evaluated in internal_Is_of_x, the
          !     long-range part is -Z*erf(alpha*x)/x
          vloc_lookup_local(x_rad_pt-x1+1) = factor * &
               internal_Is_of_x(species,x,fineness,alpha) - &
               (Z*erf_alpha_x_over_x)

       end do

       ! jd: Allreduce the partial results vloc_lookup_local -> vloc_lookup
       call comms_allgather(vloc_lookup(:,species),vloc_lookup_local)

       ! jd: Sanity-check on the result and gather some statistics on it
       maxabsreldiff = 0.0_DP
       avgabsreldiff = 0.0_DP
       avgreldiff = 0.0_DP
       avgpts = 0
       do x_rad_pt = 1, npts_x

          x = (x_rad_pt-1) * hx

#ifdef DEBUG
          ! jd: Output the calculated pseudopotential to a file
          if(pub_on_root) then
             if(x_rad_pt == 1) then
                output_unit = utils_unit()
                open(unit = output_unit, file = trim(p_species(species)&
                     %pseudo_name)//'.calculated_realpot', &
                     action = "write", err=10)
             end if
             write(output_unit,'(f20.15,f20.15)') &
                  x, vloc_lookup(x_rad_pt,species)
             if(x_rad_pt == npts_x) close(output_unit,err=20)
          end if
#endif

          ! jd: See if the tail of the potential closely matches -Z/x, it should
          if(x > 5.0_DP) then
             coulombic = -Z/x
             reldiff=(vloc_lookup(x_rad_pt,species) - coulombic)/coulombic
             absreldiff=abs(reldiff)

             if(absreldiff > maxabsreldiff) maxabsreldiff = absreldiff
             avgreldiff = avgreldiff + reldiff
             avgabsreldiff = avgabsreldiff + absreldiff
             avgpts = avgpts + 1

             if(absreldiff>0.03) then
                if(pub_on_root) then
                   write(*,*) 'x=',x
                   write(*,*) 'V_Coul=',coulombic
                   write(*,*) 'V_obtained=',vloc_lookup(x_rad_pt,species)
                end if
                call utils_abort('Discrepancy of more than 3% between the &
                     &calculated pseudopotential and Coulombic potential, in &
                     &the tail of the potential, detected in &
                     &pseudo_local_on_fine_openbc. Try increasing &
                     &openbc_pspot_finetune_f or the KE cutoff.')
             end if
          end if

       end do ! radial points

       ! jd: Make sure we do not divide by zero
       call utils_assert(avgpts > 0, &
            'Internal error in pseudo_local_on_fine_openbc, check sim. cell.')

       ! jd: Calculate and print out maximum and average relative difference
       !     from Coulombic in potential tail
       avgreldiff = avgreldiff / avgpts
       avgabsreldiff = avgabsreldiff / avgpts

       if (pub_on_root) then
          write(stdout,'(a,f14.12,a,f15.12,a,f15.12,a,i3,a,a,a)') &
               ' ',maxabsreldiff, '   ',avgabsreldiff,'  ',avgreldiff, &
               '  ',fineness,'  ''',trim(p_species(species)%pseudo_name), ''''
       end if

    end do ! species

    ! jd: Get rid of vloc_lookup_local, as it is no longer needed
    deallocate(vloc_lookup_local,stat=ierr)
    call utils_dealloc_check('pseudo_local_on_fine_openbc','vloc_lookup_local',&
         ierr)

    ! jd: Fill, in a distributed fashion, the v_local_fine using the lookups
    ! jd: Loop over all my points on the fine grid
    v_local_fine = 0.0_DP ! jd: Takes care of padding between pt and ld
    do islab12 = 1, grid%num_my_slabs12
       i3 = grid%first_slab12(pub_my_node_id) + islab12 - 1

       do i2 = 1, grid%n2
          do i1 = 1, grid%n1

             r = &
                  real((i1-1),kind=DP) * grid%da1 + &
                  real((i2-1),kind=DP) * grid%da2 + &
                  real((i3-1),kind=DP) * grid%da3

             ! jd: Loop over all the ions
             accum = 0.0_DP
             do I = 1, pub_cell%nat

                R_I = elements(I)%centre
                x = magnitude(r-R_I)
                species = elements(I)%pspecies_number

                ! jd: Interpolate Vloc(x) from radial lookups
                vloc_here = services_1d_interpolation( &
                     vloc_lookup(:,species), &
                     npts_x, x/xmax*real(npts_x-1,kind=DP),0) ! -1, sic!

                ! jd: It's faster to accumulate in a variable first
                accum = accum + vloc_here

              end do ! over I

             ! jd: Add the contribution from all ions to the current point
             v_local_fine(i1,i2,islab12) = accum

          end do ! over i1
       end do ! over i2
    end do ! over 12-slabs

    ! jd: Clean up
    deallocate(vloc_lookup,stat=ierr)
    call utils_dealloc_check('pseudo_local_on_fine_openbc','vloc_lookup',ierr)

    if (pub_on_root) then
       write(stdout,'(a/)') '... done'
    end if

    call utils_flush

    call timer_clock('pseudo_local_on_fine_openbc',2)

    return

#ifdef DEBUG
10  call utils_abort('Could not open a debug file to write the&
         & real space pseudopotential to.')
20  call utils_abort('Could not close a debug file with the&
         & real space pseudopotential.')
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real(kind=DP) function internal_Is_of_x(s,x,fineness,alpha)
      !========================================================================!
      ! Calculates the integral of                                             !
      ! ( Vloc(g) * sin(gx)/x * g * (1-exp(-g^2/(4*alpha^2))dg from 0 to g_cut,!
      ! where g_cut is the highest g supported by the current reciprocal space !
      ! grid, which will be somewhat smaller than g_max, up to which the recpot!
      ! is defined in the recpot file.                                         !
      ! The integral is tricky in that for large x it starts to oscillate very !
      ! heavily and the g resolution in the recpot file is too low to just     !
      ! integrate over the points in rad_locpot_recip, even with the Newton-   !
      ! Cotes integrator. Thus, integration is performed on a finer grid, cf.  !
      ! the parameter 'fineness'. Unless the user chose to override it, this   !
      ! will have been set to a sensible default by the parent routine.        !
      ! A check in the parent routine checks the high-x value for adequate     !
      ! closeness to the Coulombic potential.                                  !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      ! s (input)       : number of the species whose Vloc we want.            !
      ! x (input)       : the parameter x of the integration.                  !
      ! fineness (input): the fine radial g-grid will be this much finer.      !
      ! alpha (input)   : a parameter controlling the crossover between the    !
      !                   short-range and long-range parts.                    !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic in May 2010, based on notes by Chris-Kriton  !
      ! Skylaris. Updated by Jacek Dziedzic in Sep 2010 to use the trick to    !
      ! split the potential into lr and sr parts in order to improve accuracy. !
      !========================================================================!

      use rundat, only: pub_openbc_pspot_finetune_f
      use services, only: services_1d_interpolation
      use simulation_cell, only: pub_cell
      use utils, only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! jd: Arguments
      integer, intent(in) :: s           ! jd: Species in question
      real(kind=DP), intent(in) :: x     ! jd: Parameter in the integration
      integer, intent(in) :: fineness    ! jd: Fineness of the radial g-grid
      real(kind=DP), intent(in) :: alpha ! jd: Parameter

      ! jd: Local variables
      real(kind=DP), allocatable :: integrand(:) ! jd: Integrated quantity
      real(kind=DP) :: fac               ! jd: Newton-Cotes integration factor
      real(kind=DP) :: dg                ! jd: g spacing in the recpot file
      real(kind=DP) :: dg_fine           ! jd: g spacing in the fine radial grid
      real(kind=DP) :: g                 ! jd: Current g
      real(kind=DP) :: g_cut             ! jd: Maximum g in the recip grid
      real(kind=DP) :: sin_gx_over_x     ! jd: sin(gx)/x or limit thereof
      real(kind=DP) :: integral          ! jd: Accumulator for the integral
      real(kind=DP) :: v_loc_value       ! jd: Vloc(g) for current g
      integer :: ngpts, ngpts_fine       ! jd: Number of g points for both grids
      integer :: g_rad_pt_fine           ! jd: Index
      integer :: k                       ! jd: Counter
      integer :: ierr                    ! jd: Error flag

      ! -----------------------------------------------------------------------

      ! jd: Basic variables
      Z = p_species(s)%ion_charge

      ! jd: Set up the integration domain
      ngpts = p_species(s)%n_rad_pts
      dg = 1.0_DP / p_species(s)%inv_g_spacing
      dg_fine = dg/fineness

      ! jd: Ignore all g's beyond the maximum that is representable in our
      !     reciprocal space grid
      g_cut = 2.0_DP*PI/max(pub_cell%d1,pub_cell%d2,pub_cell%d3)

      ! jd: Integration weight for Newton-Cotes (25.4.17, Abramowitz, Stegun)
      fac = (7.0_DP/17280.0_DP) * dg_fine

      ! jd: Allocate array holding the integrand
      allocate(integrand(ngpts*fineness),stat=ierr)
      call utils_alloc_check('internal_Is_of_x','integrand',ierr)
      integrand(:) = 0.0_DP

      ! jd: Prepare integrand by populating the array
      k=1
      do g_rad_pt_fine = 1,ngpts*fineness-1

         g = (g_rad_pt_fine-1) * dg_fine

         if(g>g_cut) exit

         v_loc_value = services_1d_interpolation( &
              p_species(s)%rad_locpot_recip, &
              ngpts, &
              real(g_rad_pt_fine-1,kind=DP)/real(fineness,kind=DP),0)! jd:-1,sic

         if(g < 1D-10) then
            v_loc_value = 0.0_DP  ! integrand(g=0) is zero
         else

            ! jd: Add back Coulombic part, subtracted earlier
            if(alpha /= 0.0_DP) then
               v_loc_value = v_loc_value -Z/(g*g)
            end if

         end if

         if(x < 1D-10) then ! jd: Avoid singularity at 0
            sin_gx_over_x = g
         else
            sin_gx_over_x = sin(g*x)/x
         end if

         if(alpha /= 0.0_DP) then
            integrand(k) = v_loc_value * (1-exp(-g*g/(4.0_DP*alpha*alpha))) * &
                 sin_gx_over_x * g
         else
            integrand(k) = v_loc_value * sin_gx_over_x * g
         end if

         k = k+1
      end do
      ngpts_fine = k-1

      ! jd: Integrate with 7-point Newton-Cotes
      integral = 0.0_DP
      do g_rad_pt_fine=1, ngpts_fine-7,7
         integral = integral + &
              751.0_DP * &
              (integrand(g_rad_pt_fine)+integrand(g_rad_pt_fine+7)) + &
              3577.0_DP* &
              (integrand(g_rad_pt_fine+1)+integrand(g_rad_pt_fine+6)) + &
              1323.0_DP* &
              (integrand(g_rad_pt_fine+2)+integrand(g_rad_pt_fine+5)) + &
              2989.0_DP* &
              (integrand(g_rad_pt_fine+3)+integrand(g_rad_pt_fine+4))
      end do

      internal_Is_of_x = integral * fac

      deallocate(integrand,stat=ierr)
      call utils_dealloc_check('internal_Is_of_x','integrand',ierr)

    end function internal_Is_of_x

end subroutine pseudo_local_on_grid_openbc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine subtract_or_add_coul(factor)

    !=================================================================!
    ! This subroutine is used to subtract the Coulomb (Z/r) potential !
    ! from the local part of each pseudopotential species before      !
    ! interpolation.                                                  !
    !-----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 22/8/2001.                  !
    ! Modified on 24/1/2004 by Chris-Kriton Skylaris to work with the !
    ! PSEUDO_SPECIES type.                                            !
    ! Generalized on 27/5/2010 by Jacek Dziedzic to allow for adding  !
    ! 1/r back. Default behaviour is unchanged.                       !
    !=================================================================!

    use simulation_cell, only: pub_cell

    implicit none

    ! jd: Arguments
    real(kind=DP), optional, intent(in) :: factor

    ! Local Variables
    integer :: row, species
    integer :: nsp1
    ! pa: changed ion_charge from int to real
    real(kind=DP) :: g_value, ion_charge
    real(kind=DP) :: fac

    ! ----------------------------------------------------------------------

    ! jd: Default behaviour is to subtract the Coulombic potential but alllow
    !     the caller to change it by passing a factor in an optional argument
    if (present(factor)) then
       fac = factor
    else
       fac = -1.0_DP
    end if

    ! cks: loop over different atomic species
    do species=1,pub_cell%num_pspecies

       ! ndmh: no need to subtract Coulomb potential from usps
       if (.not.p_species(species)%subtract_coul) cycle

       ion_charge = p_species(species)%ion_charge
       nsp1       = p_species(species)%n_shells + 1

       ! cks: loop over radial points of current species
       do row=2,p_species(species)%n_rad_pts

          g_value = p_species(species)%rad_proj_recip(row,nsp1)

          ! cks: subtract the attractive coulomb potential
          ! cks, 28/3/2004: modification for numerical consistency with CASTEP
          ! pa: changed to allow fractional ionic charge
          ! jd: changed to allow adding the potential back
          !     fac = -1.0 -> subtraction (default), fac = 1.0 -> addition
          p_species(species)%rad_locpot_recip(row) = &
               p_species(species)%rad_locpot_recip(row) - &
               fac * ion_charge / (g_value*g_value)

       end do

    end do

  end subroutine subtract_or_add_coul

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pseudopotentials

