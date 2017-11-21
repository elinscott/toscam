!=============================================================================!
!                            S O L V A T I O N                                !
!=============================================================================!
! This module defines the implicit solvation model outlined in the following  !
! paper and references therein:                                               !
!                                                                             !
! [1] DA Scherlis, J-L Fattebert, F Gygi, M Cococcioni, and N Marzari         !
!     Journal of Chemical Physics 124, 074103  (2006)                         !
!                                                                             !
!-----------------------------------------------------------------------------!
! Written by Hatem H Helal, starting in 8/2007                                !
! please report bugs to hhh23@cam.ac.uk                                       !
!                                                                             !
! Adapted for onetep by Jacek Dziedzic, 04-05/2010                            !
! Made parallel-ready by Jacek Dziedzic, 04-05/2010                           !
!-----------------------------------------------------------------------------!

module is_solvation

  use constants, only: DP
  use is_poisson, only: POISSON_PROBLEM
  use utils, only: utils_trace_in, utils_trace_out

  implicit none

  ! Everything is private ...
  private

  !... unless exposed here.
  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!

  public :: implicit_solvent_exit
  public :: implicit_solvent_hartree

  !---------------------------------------------------------------------------!
  !                     P r i v a t e   V a r i a b l e s                     !
  !---------------------------------------------------------------------------!
  logical,                save                   :: have_initial_eps = .false.
  logical,                save                   :: have_rho_ion     = .false.
  logical,                save                   :: have_phi_ions    = .false.
  type(POISSON_PROBLEM),  save                   :: initial_problem

contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_hartree(phi_tot, rho_elec, E_Hartree, E_cav)
    !=========================================================================!
    ! This is the main driver for the implicit solvation method which will    !
    ! use the parameters 'solvation_method' and 'dielectric_model' to setup   !
    ! the solvation problem with the correct dielectric cavity and select     !
    ! either the direct or corrective approaches to solving the Poisson eq in !
    ! dielectric.  The dielectric is either self-consistently updated along   !
    ! along with the potential or else it is held fixed.  For the latter it   !
    ! is essential that the starting charge density is something sensible     !
    ! since this is used to initialize the dielectric.                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   phi_tot,          intent=out, total potential: hartree + dielectric   !
    !   rho_elec,         intent=in,  the input electronic density            !
    !   E_Hartree, (opt)  intent=out, calculated Hartree energy               !
    !   E_cav, (opt)      intent=out, cavitation energy, if requested         !
    ! Returns:                                                                !
    !   The (molecular) Hartree energy in the presence of dielectric.         !
    !   Note that in the case of self-consistently changing dielectric, this  !
    !   energy is different than the one obtained by integrating the returned !
    !   phi_tot with the molecular density. This is because the returned      !
    !   phi_tot contains additional contributions (V_eps, V_cav) from the     !
    !   derivative of the permittivity wrt changing density.                  !
    ! NB: This subroutine is spin-agnostic.                                   !
    ! Q. Why is E_cav calculated here?                                        !
    ! A. Because current_problem lives only for the duration of this routine. !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal starting 07/2009                               !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use constants, only: VERBOSE, stdout, DP
    use comms, only: pub_on_root
    use is_poisson, only: deallocate_poisson_problem
    use rundat, only: pub_is_check_solv_energy_grad, pub_is_dielectric_model, &
         pub_output_detail, pub_is_solvation_method, pub_is_include_cavitation
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out)  :: phi_tot(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP), intent(in)   :: rho_elec(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    real(kind=DP), intent(out), optional :: E_Hartree
    real(kind=DP), intent(out), optional :: E_cav

    ! Internal variables
    type(POISSON_PROBLEM) :: current_problem

    !------------------------------------------------------------------------

    call utils_trace_in('implicit_solvent_hartree')

    ! jd: Allocate arrays and evaluate dielectric.
    call initialize_solvation_problem(rho_elec, current_problem)

    select case(pub_is_solvation_method)
    case ('DIRECT')
       ! Directly solve the Poisson equation in dielectric to find the
       ! solvation potential
       call direct_solvation_potential(current_problem)

    case ('CORRECTIVE')
       ! Use the corrective potential method to get the solvation potential
       call utils_abort('Corrective approach not (yet) implemented in ONETEP.')

    case default
       call utils_abort('Unrecognized is_solvation_method in &
            &implicit_solvent_hartree.')
    end select

    if(present(E_cav)) then
       E_cav = 0.0_DP
       if(pub_is_include_cavitation) E_cav = current_problem%cavitation_energy
    end if

    ! Now we have the solvation potential we need to calculate the
    ! dielectric potential and energy correction
    select case(pub_is_dielectric_model)
    case ('FIX_INITIAL','GAUSSIAN_IONS') ! The simple case...
       if(present(E_cav)) E_cav = initial_problem%cavitation_energy

    case ('SELF_CONSISTENT') ! The more involved case...
       ! jd: Include V_eps correction (eq. (6) in [1]) in phi_mol.
       !     Energy uses uncorrected phi and has been calculated earlier, so OK.
       call calc_diel_correction(current_problem)

       current_problem%phi_mol = current_problem%phi_mol + current_problem%V_eps

       ! jd: Include V_cav correction (eq. (14) in [1]) in phi_mol.
       !     Energy uses uncorrected phi and has been calculated earlier, so OK.
       if(pub_is_include_cavitation) then
          current_problem%phi_mol = &
               current_problem%phi_mol + current_problem%V_cav
       end if

    case default
       call utils_abort('Unrecognized is_dielectric_model in &
            &implicit_solvent_hartree.')
    end select

    ! Finally we set up pot variable for returning
    phi_tot = current_problem%phi_mol

    ! Check electrostatic energy gradient if necessary
    if(pub_is_check_solv_energy_grad) then
       call check_elec_energy_gradient(current_problem)
    end if

    ! Write out data files for later analysis
    call write_solvation_data(current_problem)

    ! All done so deallocate solvation problem that is no longer needed
    call deallocate_poisson_problem(current_problem,pub_is_solvation_method)

    if(present(E_Hartree)) E_Hartree = current_problem%E_Hartree

    call utils_trace_out('implicit_solvent_hartree')

  end subroutine implicit_solvent_hartree
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine implicit_solvent_exit
    !=========================================================================!
    ! Cleans up after implicit solvent. Currently a no-op.                    !
    !=========================================================================!

    implicit none

    !------------------------------------------------------------------------

  end subroutine implicit_solvent_exit
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine initialize_solvation_problem(den,problem)
    !=========================================================================!
    ! This subroutine will allocate the necessary globally available arrays   !
    ! and initialize the dielectric depending on the chosen method.           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   den,              intent=in,     the input electron density           !
    !   problem,          intent(inout), the POISSON_PROBLEM to initialize,   !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal starting 07/2009                               !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use constants, only: DP
    use is_poisson, only: POISSON_PROBLEM, allocate_poisson_problem
    use is_smeared_ions, only: rho_ion
    use rundat, only: pub_is_dielectric_model, pub_is_solvation_method
    use utils, only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(in)            :: den(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)
    type(POISSON_PROBLEM), intent(inout) :: problem

    !------------------------------------------------------------------------

    call utils_trace_in('initialize_solvation_problem')

    ! Start by allocating solvation problem arrays
    call allocate_poisson_problem(problem,pub_is_solvation_method)

    ! Check if initial problem has been allocated and do so if needed
    if(.not. initial_problem%is_allocated) then
       call allocate_poisson_problem(initial_problem,'INITIAL')
    end if

    ! jd: Set up the rho_elec array
    problem%rho_elec = den

    ! Retrieve ion density if its already been stored in initial data
    if (have_rho_ion) then
       problem%rho_ion = initial_problem%rho_ion
    ! Otherwise this is our first time through so calculate and store it
    else
       problem%rho_ion = rho_ion
       initial_problem%rho_ion = problem%rho_ion
       have_rho_ion = .true.
    end if

    ! Find out if we are using a previously stored dielectric
    select case(pub_is_dielectric_model)
    case('FIX_INITIAL','GAUSSIAN_IONS')

      if (have_initial_eps) then
        ! The dielectric is fixed and we already have the initial data
        ! so we simply retrieve the saved dielectric and cavitation data
        problem%eps_full = initial_problem%eps_full
        problem%eps_half = initial_problem%eps_half

        problem%cavity_surface_area = initial_problem%cavity_surface_area
        problem%cavitation_energy   = initial_problem%cavitation_energy

      else

        if(pub_is_dielectric_model=='GAUSSIAN_IONS') then
           ! Use smeared ion density we calculated earlier
           initial_problem%rho_elec = -problem%rho_ion ! Note the sign...
        else
           ! We save initial electronic charge density
           initial_problem%rho_elec = problem%rho_elec
        end if

        ! Calculate dielectric on full grid
        call calc_dielectric_medium(initial_problem%rho_elec,problem%eps_full)
        initial_problem%eps_full = problem%eps_full

        ! Calculate dielectric on half grids
        call interpolate_dielectric_medium(initial_problem%rho_elec, &
             problem%eps_half)
        initial_problem%eps_half = problem%eps_half

        ! ...and also calculate the surface area and cavitation energy
        ! of the current cavity
        call calc_cavitation(initial_problem)

        problem%cavity_surface_area = initial_problem%cavity_surface_area
        problem%cavitation_energy   = initial_problem%cavitation_energy
        problem%cavity_volume       = initial_problem%cavity_volume

        ! Remember we've done that
        have_initial_eps = .true.

      end if

    case('SELF_CONSISTENT')

      ! We can go ahead and use current charge density to update dielectric
      call calc_dielectric_medium(problem%rho_elec,problem%eps_full)
      call interpolate_dielectric_medium(problem%rho_elec,problem%eps_half)

      ! Calculate the surface area and cavitation energy of the current cavity
      ! jd: In the SCF case, also calculate V_cav (eq. (14) in [1]).
      call calc_cavitation(problem)

    case default
      call utils_abort('Unrecognized pub_is_dielectric_model.')
    end select

    call utils_trace_out('initialize_solvation_problem')

  end subroutine initialize_solvation_problem
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine corrective_solvation_potential(solv_problem)
    !=========================================================================!
    ! A stub for future implementation of the corrective solvation potential. !
    !=========================================================================!

    use is_poisson, only: POISSON_PROBLEM
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(POISSON_PROBLEM), intent(inout) :: solv_problem

    ! -----------------------------------------------------------------------

    call utils_abort('Corrective solvation potential not implemented in onetep')
    solv_problem = solv_problem ! jd: Kill 'solv_problem usused' compiler warn.

  end subroutine corrective_solvation_potential
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine direct_solvation_potential(solv_problem)
    !=========================================================================!
    ! This subroutine uses the solvation type data to calculate the           !
    ! Hartree energy and potential in the presence of implicit solvent by     !
    ! directly solving the Poisson equation in dielectric.                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   solv_problem, intent(inout): the POISSON_PROBLEM in question.         !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal starting 07/2009                               !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use multigrid_methods, only: multigrid_calculate_hartree
    use rundat, only: pub_is_bulk_permittivity

    implicit none

    ! jd: Arguments
    type(POISSON_PROBLEM), intent(inout) :: solv_problem

    !------------------------------------------------------------------------

    call utils_trace_in('direct_solvation_potential')

    solv_problem%E_Hartree = multigrid_calculate_hartree( &
         solv_problem%phi_mol, &                            ! output
         solv_problem%rho_elec + solv_problem%rho_ion, &    ! input
         uniform_eps = pub_is_bulk_permittivity, &          ! input
         eps_full = solv_problem%eps_full, &                ! input
         eps_half = solv_problem%eps_half)                  ! input

    call utils_trace_out('direct_solvation_potential')

  end subroutine direct_solvation_potential
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine calc_dielectric_medium(den,dielectric)
    !=========================================================================!
    ! For a given electronic density (den) this subroutine will calculate the !
    ! dielectric functional (dielectric) for later use in solving the Poisson !
    ! equation in the presence of a dielectric.                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   den,              intent=in,  the input electron density              !
    !   dielectric        intent=out, the dielectric we seek                  !
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use constants, only: DP
    use rundat, only: pub_is_bulk_permittivity, pub_is_dielectric_function
    use utils, only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP), dimension(pub_fine_grid%ld1,pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12), intent(in)  :: den
    real(kind=DP), dimension(pub_fine_grid%ld1,pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12), intent(out) :: dielectric

    ! Internal variables
    logical :: use_fattebert

    !------------------------------------------------------------------------

    call utils_trace_in('calc_dielectric_medium')

    ! jd: Deal with the select-case first, so that we don't unnecessarily
    ! compare strings inside a tight loop.
    use_fattebert = .true. ! jd: Silence the compiler
    select case(pub_is_dielectric_function)
    case('FGF')
       use_fattebert = .true.
    case('GAUSSIAN')
       use_fattebert = .false.
    case default
       call utils_abort('Unrecognized dielectric fn in calc_dielectric_medium')
    end select

    ! jd: Fill the dielectric array, using elemental functions
    if(use_fattebert) then
       dielectric = eps_function_fattebert_gygi(den)
    else
       dielectric = eps_function_gaussian(den)
    end if

    call utils_trace_out('calc_dielectric_medium')

  end subroutine calc_dielectric_medium
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine interpolate_dielectric_medium(den,eps_half)
    !=========================================================================!
    ! For a given electronic density (den) this subroutine will use FFTs to   !
    ! interpolate the density onto 'half grids' and then calculates the       !
    ! dielectric functional (dielectric) for later use in solving the Poisson !
    ! equation in the presence of a dielectric.                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   den,              intent=in,  the input electron density              !
    !   eps_half          intent=out, the dielectric functional we seek       !
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: cell_grid_recip_pt, pub_fine_grid
    use comms, only: pub_my_node_id
    use constants, only: DP
    use fourier, only: fourier_apply_cell_forward, fourier_apply_cell_backward
    use geometry, only: magnitude
    use rundat, only: pub_is_dielectric_function, pub_is_bulk_permittivity
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    real(kind=DP), dimension(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12),  intent(in)  :: den
    real(kind=DP), dimension(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12,3),intent(out) :: eps_half

    ! Internal variables
    integer                                           :: ig1, ig2, ig3, dir
    complex(kind=DP), allocatable, dimension(:,:,:)   :: rho_recip
    complex(kind=DP), allocatable, dimension(:,:,:,:) :: phases
    real(kind=DP), allocatable, dimension(:,:,:,:)    :: rho_interp
    real(kind=DP), dimension(3)                       :: half_grids
    logical                                           :: use_fattebert
    integer                                           :: ierr ! jd: Error flag
    real(kind=DP)                                     :: d1f, d2f, d3f
    real(kind=DP)                                     :: gvec(3)

    !------------------------------------------------------------------------

    call utils_trace_in('interpolate_dielectric_medium')

    ! Allocate memory
    allocate(rho_recip(pub_fine_grid%ld3, pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs23),stat=ierr)
    call utils_alloc_check('interpolate_dielectric_medium','rho_recip',ierr)
    allocate(phases(pub_fine_grid%ld3, pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs23,3),stat=ierr)
    call utils_alloc_check('interpolate_dielectric_medium','phases',ierr)
    allocate(rho_interp(pub_fine_grid%ld1, pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12,3),stat=ierr)
    call utils_alloc_check('interpolate_dielectric_medium','rho_interp',ierr)

    ! jd: Deal with the select-case first, so that we don't unnecessarily
    ! compare strings inside a tight loop.
    use_fattebert = .true. ! jd: Silence the compiler
    select case(pub_is_dielectric_function)
    case('FGF')
       use_fattebert = .true.
    case('GAUSSIAN')
       use_fattebert = .false.
    case default
       call utils_abort('Unrecognized dielectric fn in calc_dielectric_medium')
    end select

    ! Store half-grid distances
    d1f = magnitude(pub_fine_grid%da1)
    d2f = magnitude(pub_fine_grid%da2)
    d3f = magnitude(pub_fine_grid%da3)
    half_grids(1) = d1f * 0.5_DP
    half_grids(2) = d2f * 0.5_DP
    half_grids(3) = d3f * 0.5_DP

    ! Calculate the needed phase shifts to get the density on the half grids
    ! NB: As long as pub_recip_grid orders indices differently from eps_half,
    !     there is no way this can be done cache-efficiently.
    do ig1 = 1, pub_fine_grid%num_slabs23
       do ig2 = 1, pub_fine_grid%ld2
          do ig3 = 1, pub_fine_grid%ld3
             call cell_grid_recip_pt(gvec,ig1 + &
                  pub_fine_grid%first_slab23(pub_my_node_id) - 1,ig2,ig3, &
                  pub_fine_grid)
             do dir = 1, 3
                phases(ig3,ig2,ig1,dir) = exp(cmplx(0.0_DP, &
                     gvec(dir) * half_grids(dir)))
             end do
          end do
       end do
    end do

    ! Now transform density to reciprocal space
    call fourier_apply_cell_forward(den,rho_recip,pub_fine_grid)

    ! Compute product of rho(G)*phase shift, transform to real space to get
    ! density in rho_interp
    do dir=1, 3
       phases(:,:,:,dir) = phases(:,:,:,dir) * rho_recip(:,:,:)
       call fourier_apply_cell_backward(rho_interp(:,:,:,dir), &
            phases(:,:,:,dir),pub_fine_grid)
    end do

    ! Initialize dielectric to bulk @removeme
    eps_half(:,:,:,:) = pub_is_bulk_permittivity

    ! Fill in eps_half using the interpolated density
    if(use_fattebert) then
       eps_half = eps_function_fattebert_gygi(rho_interp)
    else
       eps_half = eps_function_gaussian(rho_interp)
    end if

    ! jd: Clean up
    deallocate(rho_interp,stat=ierr)
    call utils_dealloc_check('interpolate_dielectric_medium','rho_interp',ierr)
    deallocate(phases,stat=ierr)
    call utils_dealloc_check('interpolate_dielectric_medium','phases',ierr)
    deallocate(rho_recip,stat=ierr)
    call utils_dealloc_check('interpolate_dielectric_medium','rho_recip',ierr)

    call utils_trace_out('interpolate_dielectric_medium')

  end subroutine interpolate_dielectric_medium
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine calc_cavitation(problem)
    !=========================================================================!
    ! This subroutine will take a given electron density stored in 'den' to   !
    ! calculate the surface area of the associated dielectric cavity and the  !
    ! cavitation energy.  We base this calculation upon the relationship      !
    ! between the density and the dielectric functional.  The calculation is  !
    ! done by computing the difference in volume between two nearby density   !
    ! isovalues (giving the volume of a 'shell' of charge) and dividing       !
    ! by the local thickness of this charge film or shell - giving the surface!
    ! area of the associated dielectric cavity.                               !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   solv_problem, intent(inout): the POISSON_PROBLEM in question.         !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    ! Added a volume calculation, Jacek Dziedzic 06/2010                      !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: pub_on_root, comms_reduce
    use constants, only: DP
    use finite_differences, only: finite_difference_gradient
    use multigrid_methods, only: mg
    use rundat, only: pub_is_density_threshold, pub_is_surface_thickness, &
         pub_is_solvent_surface_tension, pub_is_discretization_order, &
         pub_is_dielectric_model, pub_is_include_cavitation
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_sanity_check

    implicit none

    ! Arguments
    type(POISSON_PROBLEM), intent(inout) :: problem

    ! Internal variables
    real(kind=DP), dimension(:,:,:,:), allocatable :: grad_den
    real(kind=DP), dimension(:,:,:,:,:), allocatable :: second_der
    real(kind=DP) :: lcl_mod_grad_den
    real(kind=DP) :: lcl_den
    ! Density thresholds which are used to define isosurfaces
    real(kind=DP) :: dthresh_a,dthresh_b
    real(kind=DP) :: smooth_a,smooth_b,smooth_0    !smoothed density isosurfaces
    integer       :: i1, i2, i3
    integer       :: ierr
    real(kind=DP) :: cavity_surface_area, cavity_volume
    logical       :: want_scf
    real(kind=DP) :: V_cav_here
    integer       :: i, j ! jd: x, y or z

    !------------------------------------------------------------------------

    call utils_trace_in('calc_cavitation')

    want_scf = (pub_is_dielectric_model == 'SELF_CONSISTENT')

    ! jd: If cavitation is off, don't bother with all this
    if(.not. pub_is_include_cavitation) then
       problem%cavitation_energy = 0D0
       problem%cavity_surface_area = 0D0
       problem%cavity_volume = 0D0
       if(want_scf) problem%V_cav = 0D0
       call utils_trace_out('calc_cavitation')
       return
    end if

    ! Allocate memory for density gradient calculation
    allocate(grad_den(pub_fine_grid%ld1, pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12,3),stat=ierr)
    call utils_alloc_check('calc_cavitation','grad_den',ierr)

    ! Use finite differences to calculate density gradient
    call finite_difference_gradient(grad_den, &              ! out
         problem%rho_elec,pub_is_discretization_order,mg)    ! in

    ! jd: Allocate memory for second derivatives, needed only in the
    !     self-consistent version
    if(want_scf) then
       allocate(second_der(pub_fine_grid%ld1, pub_fine_grid%ld2, &
            pub_fine_grid%max_slabs12,3,3),stat=ierr)
       call utils_alloc_check('calc_cavitation','second_der',ierr)
       do i=1, 3
          call finite_difference_gradient(second_der(:,:,:,:,i), &   ! out
               grad_den(:,:,:,i),pub_is_discretization_order,mg)     ! in
       end do
    end if

    ! Use density_threshold and surface_thickness parameters to define two
    ! isosurfaces for finite difference calculation
    dthresh_a = pub_is_density_threshold - 0.5_DP * pub_is_surface_thickness
    dthresh_b = pub_is_density_threshold + 0.5_DP * pub_is_surface_thickness

    ! jd: Calculate the cavity surface area
    ! @improveme: not cache friendly wrt grad_den
    ! @improveme: not cache friendly wrt second_der
    cavity_surface_area = 0.0_DP
    cavity_volume = 0.0_DP
    if(want_scf) problem%V_cav = 0.0_DP ! jd: Take care of pad between pt and ld

    do i3=1, pub_fine_grid%num_my_slabs12
       do i2=1, pub_fine_grid%n2
          do i1=1, pub_fine_grid%n1

             ! Calculate local value of |grad den|
             lcl_mod_grad_den = sqrt(grad_den(i1,i2,i3,1)**2 + &
                  grad_den(i1,i2,i3,2)**2 + grad_den(i1,i2,i3,3)**2)

             ! And local value of density
             lcl_den = problem%rho_elec(i1,i2,i3)

             smooth_a = smoothing_function(lcl_den,dthresh_a)
             smooth_b = smoothing_function(lcl_den,dthresh_b)
             smooth_0 = smoothing_function(lcl_den,pub_is_density_threshold)

             cavity_surface_area = cavity_surface_area + (smooth_a-smooth_b) * &
                  lcl_mod_grad_den / pub_is_surface_thickness

             if(want_scf) then
                V_cav_here = 0D0

                ! jd: Take into account only points where the magnitude of the
                !     gradient is not infinitesimally small. This avoids trouble
                !     in the corner case where a centre of a symmetrical mole-
                !     cule lies exactly on a grid point, yielding something like
                !     1E-12 in the gradient, which then makes V_cav_here explode
                !     The term should then be zero.
                if(abs(lcl_mod_grad_den) > 1D-8) then
                   do i=1, 3
                      do j=1, 3
                         V_cav_here = V_cav_here + &
                              (grad_den(i1,i2,i3,j) * grad_den(i1,i2,i3,i) * &
                              second_der(i1,i2,i3,j,i)) / &
                              (lcl_mod_grad_den * lcl_mod_grad_den)
                      end do
                      V_cav_here = V_cav_here - second_der(i1,i2,i3,i,i)
                   end do

                   V_cav_here = V_cav_here / lcl_mod_grad_den * &
                        (smooth_a-smooth_b) * &
                        pub_is_solvent_surface_tension/pub_is_surface_thickness
                else
                   ! jd: No-op -- if lcl_mod_grad_den is zero, we avoid /0
                   !     by letting V_cav_here remain zero.
                end if

                problem%V_cav(i1,i2,i3) = V_cav_here
             end if

             cavity_volume = cavity_volume + smooth_0

          end do
       end do
    end do

    call comms_reduce('SUM',cavity_surface_area)
    call comms_reduce('SUM',cavity_volume)

    ! Now normalize the surface area and the volume
    cavity_surface_area = cavity_surface_area * pub_fine_grid%weight
    cavity_volume = cavity_volume * pub_fine_grid%weight

    ! Update the problem
    problem%cavitation_energy = &
         pub_is_solvent_surface_tension * cavity_surface_area
    problem%cavity_surface_area = cavity_surface_area
    problem%cavity_volume = cavity_volume

    ! jd: Clean up
    if(want_scf) then
       call utils_sanity_check(problem%V_cav(:,:,:),'V_cav calc_cavitation')
       deallocate(second_der,stat=ierr)
       call utils_dealloc_check('calc_cavitation','second_der',ierr)
    end if
    deallocate(grad_den,stat=ierr)
    call utils_dealloc_check('calc_cavitation','grad_den',ierr)

    call utils_trace_out('calc_cavitation')

  end subroutine calc_cavitation
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine calc_diel_correction(solv_problem)
    !=========================================================================!
    ! Given the electrostatic potential computed in the presence of a         !
    ! dielectric we compute the contribution to the Kohn-Sham potential       !
    ! arising from the dielectric depending on the charge density. This is    !
    ! necessary to ensure that the charge density and dielectric can respond  !
    ! to each other within the SCF procedure such that we achieve a tunable   !
    ! dielectric cavity for our solvated electronic structure calculation     !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   solv_problem, intent(inout): the POISSON_PROBLEM in question.         !
    ! ----------------------------------------------------------------------- !
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use constants, only: DP, PI
    use is_poisson, only: POISSON_PROBLEM
    use rundat, only: pub_is_dielectric_model, pub_fine_grid_scale, &
         pub_is_multigrid_nlevels, pub_is_discretization_order
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert

    implicit none

    ! Arguments
    type(POISSON_PROBLEM), intent(inout) :: solv_problem

    ! Internal variables
    real(kind=DP), dimension(:,:,:), allocatable :: deps_drho
    real(kind=DP), parameter :: inv_eight_pi = 1.0_DP/(8.0_DP * PI)
    integer :: ierr

    !------------------------------------------------------------------------

    call utils_trace_in('calc_diel_pot_and_corr')

    call utils_assert(pub_is_dielectric_model == 'SELF_CONSISTENT', &
         'Internal error in is_solvation_mod, calc_diel_correction')

    ! Allocate and get the dielectric medium derivative
    allocate(deps_drho(pub_fine_grid%ld1, pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('calc_diel_pot_and_corr','deps_drho',ierr)

    deps_drho = dielectric_functional_deriv(solv_problem%rho_elec)

    solv_problem%V_eps = -inv_eight_pi * deps_drho

    deallocate(deps_drho,stat=ierr)
    call utils_dealloc_check('calc_diel_pot_and_corr','deps_drho',ierr)

    ! Need to make sure the global grad_phi_sqd is set up
    call get_grad_phi_sqd(solv_problem)

    solv_problem%V_eps = solv_problem%V_eps * solv_problem%grad_phi_sqd

    call utils_trace_out('calc_diel_pot_and_corr')

  end subroutine calc_diel_correction
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine get_grad_phi_sqd(solv_problem)
    !=========================================================================!
    ! Compute the gradient of the solvation potential and store in global     !
    ! array which is needed in both the energy evaluation and the dielectric  !
    ! potential. Uses global flag have_grad_phi_sqd                           !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   solv_problem, intent(inout): the POISSON_PROBLEM in question.         !
    ! ----------------------------------------------------------------------- !
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use constants, only: DP
    use finite_differences, only: finite_difference_mod_grad_sqd
    use multigrid_methods, only: mg
    use is_poisson, only: POISSON_PROBLEM
    use rundat, only: pub_is_discretization_order
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(POISSON_PROBLEM), intent(inout) :: solv_problem

    !------------------------------------------------------------------------

    call utils_trace_in('get_grad_phi_sqd')

    if (.not. solv_problem%have_grad_phi_sqd) then
       ! jd: Get the gradient
       call finite_difference_mod_grad_sqd(solv_problem%grad_phi_sqd, &
            solv_problem%phi_mol, pub_is_discretization_order,mg)

       ! Lastly, set the flag so we don't have to repeat this process again
       solv_problem%have_grad_phi_sqd = .true.

    end if

    call utils_trace_out('get_grad_phi_sqd')

  end subroutine get_grad_phi_sqd
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine check_elec_energy_gradient(initial)
    !=========================================================================!
    ! Takes the given charge density and compares the numerical electrostatic !
    ! energy gradient to the analytic.                                        !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   initial, intent(inout): the initial POISSON_PROBLEM.                  !
    ! ----------------------------------------------------------------------- !
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: pub_on_root, comms_reduce
    use constants, only: DP, VERBOSE, stdout
    use is_poisson, only: POISSON_PROBLEM, allocate_poisson_problem, &
         deallocate_poisson_problem
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_output_detail, pub_is_dielectric_model, &
         pub_is_solvation_method, pub_is_include_cavitation
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort

    implicit none

    ! Arguments
    type(POISSON_PROBLEM), intent(inout) :: initial

    ! Internal variables
    type(POISSON_PROBLEM) :: varied
    real(kind=DP), allocatable, dimension(:,:,:) :: gradient
    real(kind=DP)                   :: actual,predicted,ratio,Q_varied
    integer :: ierr

    real(kind=DP), parameter :: alpha = 1.0D-5 ! jd: step length

    !------------------------------------------------------------------------

    call utils_trace_in('check_elec_energy_gradient')

    ! Allocate memory
    call allocate_poisson_problem(varied,pub_is_solvation_method)
    allocate(gradient(pub_fine_grid%ld1, pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('check_elec_energy_gradient','gradient',ierr)

    ! Now copy functional derivative into gradient (V_eps = 0 for fixed eps)
    gradient = initial%phi_mol

    if(pub_is_dielectric_model == 'SELF_CONSISTENT') then
       gradient = gradient + initial%V_eps
       if(pub_is_include_cavitation) then
          gradient = gradient + initial%V_cav
       end if
    end if

    ! Now use the gradient to define a new density
    varied%rho_elec = initial%rho_elec + alpha * gradient
    varied%rho_ion = initial%rho_ion ! direct_solvation_potential relies on this

    ! Correctly evaluate the dielectric
    select case (pub_is_dielectric_model)
    case('FIX_INITIAL','GAUSSIAN_IONS')
       ! hhh: Use the same eps as initial problem
       varied%eps_full = initial%eps_full
       varied%eps_half = initial%eps_half
    case('SELF_CONSISTENT') ! Must evaluate the new eps for the varied problem
       call calc_dielectric_medium(varied%rho_elec,varied%eps_full)
       call interpolate_dielectric_medium(varied%rho_elec,varied%eps_half)
    case default
       call utils_abort('Unrecognized pub_is_dielectric_model')
    end select

    ! And now we need to solve the varied solvation problem
    select case(pub_is_solvation_method)
    case ('DIRECT')
       ! Directly solve the Poisson equation in dielectric to find the solv. pot
       call direct_solvation_potential(varied)
    case ('CORRECTIVE')
       ! Use the corrective potential method to get the solvation potential
       call corrective_solvation_potential(varied)
    case default
      call utils_abort('Unrecognized pub_is_solvation_method')
    end select

    ! Compute the numerical energy gradient
    actual     = (varied%E_Hartree - initial%E_Hartree)/alpha

    ! Now we calculate the analytic gradient which is just
    ! int gradient * gradient dr

    ! We also calculate the total charge in rho_varied while we are at it
    predicted = integrals_product_on_grid(pub_fine_grid, &
         gradient, gradient)
    Q_varied = integrals_product_on_grid(pub_fine_grid, &
         varied%rho_elec)

    ! Calculate the ratio actual/predicted ('badness')
    ratio = actual/predicted

    ! Output
    if(pub_on_root .and. pub_output_detail == VERBOSE) then
      write(stdout,*) 'Electrostatic energy gradient check:'
      write(stdout,*) '  Q_varied      ', Q_varied
      write(stdout,*) '  Step length   ', alpha
      write(stdout,*) '  E_varied      ', varied%E_Hartree
      write(stdout,*) '  E_initial     ', initial%E_Hartree
      write(stdout,*) '  predicted     ', predicted
      write(stdout,*) '  actual        ', actual
      write(stdout,*) '  badness       ',ratio
    end if

    ! Store energy gradient in initial solvation problem
    initial%energy_gradient_ratio_rho = ratio
    initial%energy_gradient_ratio_eps = 0.0_DP ! jd: Removed functionality

    ! Remember
    initial%have_energy_gradient=.true.

    ! Clean up
    call deallocate_poisson_problem(varied,pub_is_solvation_method)
    deallocate(gradient,stat=ierr)
    call utils_dealloc_check('check_elec_energy_gradient','gradient',ierr)

    call utils_trace_out('check_elec_energy_gradient')

  end subroutine check_elec_energy_gradient
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine write_solvation_data(solv_problem)
    !=========================================================================!
    ! Writes out grid data for the solvation problem.                         !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   solv_problem, intent(inout): the POISSON_PROBLEM in question.         !
    ! ----------------------------------------------------------------------- !
    ! Written by Hatem H Helal                                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    !=========================================================================!

    use is_poisson, only: POISSON_PROBLEM, write_poisson_problem
    use rundat, only: pub_is_dielectric_model, pub_is_solvation_method, &
         pub_is_solvation_output_detail

    implicit none

    ! Arguments
    type(POISSON_PROBLEM), intent(in)  :: solv_problem

    !------------------------------------------------------------------------

    call utils_trace_in('write_solvation_data')

    ! Write out all potentials associated with solvation current problem
    if(pub_is_solvation_output_detail /= 'NONE') then
       call write_poisson_problem(solv_problem,'current', &
            pub_is_solvation_method, pub_is_dielectric_model)

       ! Also write out initial data in initial_problem if we have
       ! a fixed eps or gaussian ions
       if(have_initial_eps) then
          call write_poisson_problem(initial_problem,'initial','INITIAL', &
               'NONE')
       end if
    end if

    call utils_trace_out('write_solvation_data')

  end subroutine write_solvation_data
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



  !===========================================================================!
  !                           Private functions                               !
  !===========================================================================!

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  elemental function eps_function_fattebert_gygi(rho_r)
    !=========================================================================!
    ! Smooth dielectric function, proposed by Jean-Luc Fattebert and          !
    ! Francois Gygi.  See J Comp Chem 23 pgs 662-666                          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rho_r,            intent=in,  the electronic charge density at r      !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal      07/09/2007                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    !=========================================================================!

    use constants, only: DP
    use rundat, only: pub_is_bulk_permittivity, pub_is_density_threshold, &
         pub_is_solvation_beta

    implicit none

    ! Return type
    real(kind=DP) :: eps_function_fattebert_gygi

    ! Arguments
    real(kind=DP), intent(in) :: rho_r

    ! Internal variables
    real(kind=DP) :: powbeta, rho

    !------------------------------------------------------------------------

    rho = rho_r
    if(rho < 0.0_DP) rho = 0.0_DP ! jd: Filter out ringing and noise

    powbeta = (rho/pub_is_density_threshold)**(2.0_DP*pub_is_solvation_beta)

    eps_function_fattebert_gygi = &
         1.0_DP + ((pub_is_bulk_permittivity-1.0_DP)/2.0_DP) * &
         ( 1.0_DP + (1.0_DP-powbeta)/(1.0_DP+powbeta) )

  end function eps_function_fattebert_gygi
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  elemental function eps_function_gaussian(rho_r)
    !=========================================================================!
    ! Smooth, Gaussian-shaped dielectric function depending on only one param.!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rho_r,            intent=in,  the electronic charge density at r      !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal      21/09/2009                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    !=========================================================================!

    use constants, only: DP
    use rundat, only: pub_is_density_threshold, pub_is_bulk_permittivity

    implicit none

    ! Return type
    real(kind=DP) :: eps_function_gaussian

    ! Arguments
    real(kind=DP), intent(in) :: rho_r

    ! Internal variables
    real(kind=DP) :: rho, rho_frac

    !------------------------------------------------------------------------

    rho = rho_r
    if(rho < 0.0_DP) rho = 0.0_DP ! jd: Filter out ringing and noise

    rho_frac = rho/pub_is_density_threshold

    eps_function_gaussian = 1.0_DP + &
         (pub_is_bulk_permittivity-1.0_DP) * exp(-rho_frac**2.0_DP)

  end function eps_function_gaussian
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  elemental function dielectric_functional_deriv(rho_r)
    !=========================================================================!
    ! Computes the functional derivative of the dielectric                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! rho_r, intent =in,  the electronic charge density at r                  !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal 2/11/2007                                      !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    !=========================================================================!

    use constants, only: DP
    use rundat, only: pub_is_density_threshold, pub_is_bulk_permittivity, &
         pub_is_solvation_beta, pub_is_dielectric_function

    implicit none

    ! Return type
    real(kind=DP) :: dielectric_functional_deriv

    ! Arguments
    real(kind=DP), intent(in) :: rho_r

    ! Internal variables
    real(kind=DP)             :: num,denom,rho,rho_frac,twobeta

    !------------------------------------------------------------------------

    rho = rho_r
    if(rho < 0.0_DP) rho = 0.0_DP ! jd: Filter out ringing and noise

    rho_frac = rho/pub_is_density_threshold

    twobeta  = 2.0_DP*pub_is_solvation_beta

    num = (1.0_DP-pub_is_bulk_permittivity)*(twobeta*rho_frac**(twobeta-1.0_DP))

    denom = pub_is_density_threshold*(1.0_DP+rho_frac**twobeta)**2.0_DP

    dielectric_functional_deriv = num/denom

  end function dielectric_functional_deriv
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  function smoothing_function(rho_r,thresh)
    !=========================================================================!
    ! Smooth threshold function used for computing cavity volume and surface  !
    ! cf. [1]                                                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rho_r,            intent=in,  the electronic charge density at r      !
    !   thresh,           intent=in,  the threshold parameter                 !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal      24/08/2007                                !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    !=========================================================================!

    use constants, only: DP
    use rundat, only: pub_is_solvation_beta

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: rho_r
    real(kind=DP), intent(in) :: thresh

    ! Internal variables
    real(kind=DP) :: powbeta, rho

    ! Return type
    real(kind=DP) :: smoothing_function

    !------------------------------------------------------------------------

    rho = rho_r
    if(rho < 0.0_DP) rho = 0.0_DP ! jd: Filter out ringing and noise

    powbeta = (rho/thresh)**(2.0_DP*pub_is_solvation_beta)

    smoothing_function = 0.5_DP*( (powbeta-1.0_DP) / (powbeta+1.0_DP) + 1.0_DP )

  end function smoothing_function
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module is_solvation
