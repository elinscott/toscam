!=============================================================================!
!                              P O I S S O N                                  !
!=============================================================================!
! This module defines a general data class for collecting and organizing      !
! arrays related to solving the generalized Poisson equation.  Initially      !
! developed for use with implicit solvent but can be used in other modules    !
! for straightforward data output.                                            !
!-----------------------------------------------------------------------------!
! Written by Hatem H Helal, starting in 11/2009                               !
! please report bugs to hhh23(at)cam.ac.uk                                    !
!                                                                             !
! Adapted for onetep by Jacek Dziedzic, 04-05/2010                            !
! Made parallel-ready by Jacek Dziedzic, 04-05/2010                           !
!-----------------------------------------------------------------------------!
module is_poisson

  use constants, only: DP, ANGSTROM
  use utils, only: utils_trace_in, utils_trace_out

  implicit none

  ! Everything is private ...
  private

  !... unless exposed here.
  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!

  public :: allocate_poisson_problem
  public :: deallocate_poisson_problem
  public :: write_poisson_problem

  !---------------------------------------------------------------------------!
  !                       P u b l i c   V a r i a b l e s                     !
  !---------------------------------------------------------------------------!
  ! jd: Holds all data relevant to the solution of the Poisson-Boltzmann      !
  !     equation. All arrays except eps_half represent quantities on the fine !
  !     grid and are stored in usual onetep 12slab-distributed fashion.       !
  !     eps_half holds three such arrays, with the last index indexing the    !
  !     direction in which the grid is shifted by half a spacing.             !
  !---------------------------------------------------------------------------!
  ! jd: Some of the fields in this data structure, denoted with (*), were     !
  !     retained for compatibility with the CASTEP version. They are not used !
  !     currently, but were left over in case the CORRECTIVE and/or FFT       !
  !     approach is implemented, or if write_poisson_problem is ever used.    !
  !---------------------------------------------------------------------------!
  type, public :: POISSON_PROBLEM

     ! jd: electronic density
     real(kind=DP), allocatable, dimension(:,:,:)       :: rho_elec
     ! jd: (smeared) ion density
     real(kind=DP), allocatable, dimension(:,:,:)       :: rho_ion
     ! jd: dielectric permittivity
     real(kind=DP), allocatable, dimension(:,:,:)       :: eps_full
     ! jd: dielectric permittivity, shifted by half the fine grid spacing
     !     last index picks the direction of the shift
     real(kind=DP), allocatable, dimension(:,:,:,:)     :: eps_half
     ! jd: Potential in PBCs (*)
     real(kind=DP), allocatable, dimension(:,:,:)       :: phi_pbc
     ! jd: Potential due to rho_elec+rho_ion
     real(kind=DP), allocatable, dimension(:,:,:)       :: phi_mol
     ! jd: Corrective potential (*)
     real(kind=DP), allocatable, dimension(:,:,:)       :: phi_corr
     ! jd: Correction to the energy gradient due to cavity changing shape
     real(kind=DP), allocatable, dimension(:,:,:)       :: V_eps
     ! jd: As above, the cavitational term
     real(kind=DP), allocatable, dimension(:,:,:)       :: V_cav
     ! jd: Temporary for |grad phi|^2
     real(kind=DP), allocatable, dimension(:,:,:)       :: grad_phi_sqd

     ! jd: Hartree energy for the above problem
     real(kind=DP)                                  :: E_Hartree
     ! jd: Hartree energy calculated from \int rho phi dr (*)
     real(kind=DP)                                  :: E_from_rho
     ! jd: Hartree energy calculated from \int grad phi^2 eps dr (*)
     real(kind=DP)                                  :: E_from_eps

     ! jd: Gradient correctness ratio
     real(kind=DP)                                  :: energy_gradient_ratio_rho
     ! jd: (*)
     real(kind=DP)                                  :: energy_gradient_ratio_eps

     ! jd: Surface area of the dielectric cavity
     real(kind=DP)                                  :: cavity_surface_area
     ! jd: Cavitation energy, possibly including dispersion-repulsion
     real(kind=DP)                                  :: cavitation_energy
     ! jd: Volume of the cavity
     real(kind=DP)                                  :: cavity_volume

     ! jd: .true. if datastructures allocated
     logical                                    :: is_allocated = .false.
     ! jd: .true. if grad_phi_sqd up to date
     logical                                    :: have_grad_phi_sqd = .false.
     ! jd: .true. if energy_gradient_ratio_rho up to date
     logical                                    :: have_energy_gradient= .false.

  end type POISSON_PROBLEM

contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine allocate_poisson_problem(problem,method)
    !=========================================================================!
    ! This subroutine will allocate the arrays associated with the            !
    ! poisson_problem data type.  Method argument is used to allocate extra   !
    ! arrays if the corrective method is chosen.                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   problem,        intent=inout,  poisson_problem data type              !
    !   method,         intent=in,     method for solving poisson problem     !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal starting 07/2009                               !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use rundat, only: pub_is_dielectric_model
    use utils, only: utils_alloc_check

    implicit none

    ! jd: Arguments
    type(POISSON_PROBLEM), intent(inout) :: problem
    character(len=*), intent(in) :: method

    ! jd: Internal variables
    integer :: ierr ! jd: Error flag
    integer :: ld1f, ld2f, max_slabs12

    !------------------------------------------------------------------------

    call utils_trace_in('allocate_poisson_problem')

    ld1f = pub_fine_grid%ld1
    ld2f = pub_fine_grid%ld2
    max_slabs12 = pub_fine_grid%max_slabs12

    allocate(problem%rho_elec(ld1f,ld2f,max_slabs12),stat=ierr)
    call utils_alloc_check('allocate_poisson_problem','rho_elec',ierr)

    ! Allow for standard Hartree evaluation
    if(method=='FFT') then
       allocate(problem%phi_pbc(ld1f,ld2f,max_slabs12),stat=ierr)
       call utils_alloc_check('allocate_poisson_problem','phi_pbc',ierr)
    else
       allocate(problem%rho_ion(ld1f,ld2f,max_slabs12),stat=ierr)
       call utils_alloc_check('allocate_poisson_problem','rho_ion',ierr)
       allocate(problem%eps_full(ld1f,ld2f,max_slabs12),stat=ierr)
       call utils_alloc_check('allocate_poisson_problem','eps_full',ierr)
       allocate(problem%eps_half(ld1f,ld2f,max_slabs12,3),stat=ierr)
       call utils_alloc_check('allocate_poisson_problem','eps_half',ierr)
       if(method /='INITIAL') then
          allocate(problem%phi_mol(ld1f,ld2f,max_slabs12),stat=ierr)
          call utils_alloc_check('allocate_poisson_problem','phi_mol',ierr)
          if(pub_is_dielectric_model == 'SELF_CONSISTENT') then
             allocate(problem%V_eps(ld1f,ld2f,max_slabs12),stat=ierr)
             call utils_alloc_check('allocate_poisson_problem','V_eps',ierr)
             allocate(problem%V_cav(ld1f,ld2f,max_slabs12),stat=ierr)
             call utils_alloc_check('allocate_poisson_problem','V_cav',ierr)
          end if
          allocate(problem%grad_phi_sqd(ld1f,ld2f,max_slabs12),stat=ierr)
          call utils_alloc_check('allocate_poisson_problem','grad_phi_sqd',ierr)
          ! Need more arrays if we are using the corrective potential method
          if(method=='CORRECTIVE') then
             allocate(problem%phi_pbc(ld1f,ld2f,max_slabs12),stat=ierr)
             call utils_alloc_check('allocate_poisson_problem','phi_pbc',ierr)
             allocate(problem%phi_corr(ld1f,ld2f,max_slabs12),stat=ierr)
             call utils_alloc_check('allocate_poisson_problem','phi_corr',ierr)
          end if
       end if
    end if

    problem%is_allocated=.true.

    call utils_trace_out('allocate_poisson_problem')

  end subroutine allocate_poisson_problem
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine deallocate_poisson_problem(problem,method)
    !=========================================================================!
    ! This subroutine will deallocate the arrays associated with the          !
    ! poisson_problem data type.  Method argument is used to deallocate       !
    ! extra arrays if the corrective method is chosen.                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   problem,        intent=inout,  poisson_problem data type              !
    !   method,         intent=in,     method for solving poisson problem     !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal starting 07/2009                               !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use rundat, only: pub_is_dielectric_model
    use utils, only: utils_dealloc_check

    implicit none

    ! jd: Arguments
    type(POISSON_PROBLEM), intent(inout) :: problem
    character(len=*), intent(in) :: method

    ! jd: Local variables
    integer :: ierr ! jd: Error flag

    !------------------------------------------------------------------------

    call utils_trace_in('deallocate_poisson_problem')

    !@improveme: should make an attempt to deallocate in reverse order
    if(method=='FFT') then
       deallocate(problem%phi_pbc,stat=ierr)
       call utils_dealloc_check('deallocate_poisson_problem','phi_pbc',ierr)
    else
       deallocate(problem%rho_ion,stat=ierr)
       call utils_dealloc_check('deallocate_poisson_problem','rho_ion',ierr)
       deallocate(problem%eps_full,stat=ierr)
       call utils_dealloc_check('deallocate_poisson_problem','eps_full',ierr)
       deallocate(problem%eps_half,stat=ierr)
       call utils_dealloc_check('deallocate_poisson_problem','eps_half',ierr)
       if(method /='INITIAL') then
          deallocate(problem%phi_mol,stat=ierr)
          call utils_dealloc_check('deallocate_poisson_problem','phi_mol',ierr)
          if(pub_is_dielectric_model == 'SELF_CONSISTENT') then
             deallocate(problem%V_cav,stat=ierr)
             call utils_dealloc_check('deallocate_poisson_problem','V_cav',ierr)
             deallocate(problem%V_eps,stat=ierr)
             call utils_dealloc_check('deallocate_poisson_problem','V_eps',ierr)
          end if
          deallocate(problem%grad_phi_sqd,stat=ierr)
          call utils_dealloc_check('deallocate_poisson_problem','grad_phi_sqd',&
               ierr)
          ! Need more arrays if we are using the corrective potential method
          if(method=='CORRECTIVE') then
             deallocate(problem%phi_pbc,stat=ierr)
             call utils_dealloc_check('deallocate_poisson_problem','phi_pbc',&
                  ierr)
             deallocate(problem%phi_corr,stat=ierr)
             call utils_dealloc_check('deallocate_poisson_problem','phi_corr',&
                  ierr)
          end if
       end if
    end if

    deallocate(problem%rho_elec,stat=ierr)
    call utils_dealloc_check('deallocate_poisson_problem','rho_elec',ierr)

    problem%is_allocated=.false.

    call utils_trace_out('deallocate_poisson_problem')

  end subroutine deallocate_poisson_problem
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine write_poisson_problem(problem,root,method,model)
    !=========================================================================!
    ! This subroutine will write out a bunch of data files associated with    !
    ! the given poisson problem.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   problem,        intent=in,     poisson_problem data type              !
    !   root,           intent=in,     'initial', 'current', etc.             !
    !   method,         intent=in,     method for solving poisson problem     !
    !   model,          intent=in,     the dielectric model being used        !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal starting 07/2009                               !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use visual, only: visual_scalarfield

    implicit none

    ! jd: Arguments
    type(POISSON_PROBLEM), intent(in) :: problem
    character(len=*), intent(in) :: method
    character(len=*), intent(in) :: root
    character(len=*), intent(in) :: model

    ! jd: Local variables
    character(len=200) :: suffix
    integer, save :: counter = 0

    !------------------------------------------------------------------------

    call utils_trace_in('write_poisson_problem')

    counter = counter+1
    write(suffix,'(i0)') counter

    if(method=='FFT') then
       ! Write out the charge density
       call visual_scalarfield(problem%rho_elec, &
            pub_fine_grid, 'Charge density',root//'_rho_'//suffix, &
            conversion_factor = ANGSTROM**3)
       ! Write out the periodic potential
       call visual_scalarfield(problem%phi_pbc, pub_fine_grid, &
            'Periodic potential found by FFTs',root//'_phi_pbc_'//suffix)
    end if

    if(method=='INITIAL') then
       ! Only want to write out the charge density used to initialize the
       ! dielectric, this could either be the initial electronic density
       ! or else a sum of Gaussians on ions
       call visual_scalarfield(problem%rho_elec, &
            pub_fine_grid, 'Initial charge density',root//'_rho_'//suffix, &
            conversion_factor = ANGSTROM**3)
    end if

    if(method=='DIRECT' .or. method=='CORRECTIVE') then
       ! Write out the solvation potential
       call visual_scalarfield(problem%phi_mol, pub_fine_grid, &
            'Solvation potential',root//'_phi_mol_'//suffix)

       ! Write out the dielectric on the full grid
       call visual_scalarfield(problem%eps_full, pub_fine_grid, &
            'Dielectric on full grid',root//'_eps_full_'//suffix)

       ! Write out the charge density - scaled with 1/V to get true density
       call visual_scalarfield(problem%rho_elec, pub_fine_grid, &
            'Electronic charge density',root//'_rho_'//suffix, &
            conversion_factor = ANGSTROM**3)

       ! Write out the ion charge density - scaled with 1/V to get true density
       call visual_scalarfield(problem%rho_ion, pub_fine_grid, &
            'Smeared ion charge density', root//'_rho_ion_'//suffix, &
            conversion_factor = ANGSTROM**3)

       ! Write out the molecular charge density - scaled with 1/V as above
       call visual_scalarfield(problem%rho_elec+problem%rho_ion, &
            pub_fine_grid, 'Molecular charge density = elec + ion', &
            root//'_rho_mol_'//suffix, conversion_factor = ANGSTROM**3)

       if(model=='SELF_CONSISTENT') then
          ! Write out the dielectric potential
          call visual_scalarfield(problem%V_eps, pub_fine_grid, &
               'Dielectric potential term V_eps',root//'_V_eps_'//suffix)
          call visual_scalarfield(problem%V_cav, pub_fine_grid, &
               'Dielectric potential term V_cav',root//'_V_cav_'//suffix)
       end if

       if(method=='CORRECTIVE') then
          ! Write out the periodic potential
          call visual_scalarfield(problem%phi_pbc, pub_fine_grid, &
               'Periodic potential found by FFTs',root//'_phi_pbc_'//suffix, &
               conversion_factor = ANGSTROM**3)

          ! Write out the corrective potential
          call visual_scalarfield(problem%phi_corr, pub_fine_grid, &
               'Corrective potential',root//'_phi_corr_'//suffix)
       end if
    end if

    call utils_trace_out('write_poisson_problem')

  end subroutine write_poisson_problem
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module is_poisson
