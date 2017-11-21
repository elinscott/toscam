!=============================================================================!
!                         I S_ S M E A R E D _ I O N S                        !
!=============================================================================!
! This module adds the smeared ion representation due to                      !
! [1] Scherlis et al, described in J. Chem. Phys. 124 (2006).                 !
!-----------------------------------------------------------------------------!
! Written by Jacek Dziedzic in May 2010, based largely on the abovementioned  !
! paper (Appendix A) and a similar module for CASTEP due to Hatem H Helal.    !
!-----------------------------------------------------------------------------!

module is_smeared_ions

  use constants, only: DP
  use utils, only: utils_trace_in, utils_trace_out

  implicit none
  private

  public :: smeared_ion_apply_Vloc_corr
  public :: smeared_ion_calc_E_self
  public :: smeared_ion_calc_E_smeared
  public :: smeared_ion_exit
  public :: smeared_ion_hartree
  public :: smeared_ion_initialise

  real(kind=DP), public, save :: smeared_ion_E_self
  real(kind=DP), public, save :: smeared_ion_E_smeared

  real(kind=DP), public, allocatable, save :: rho_ion(:,:,:)


contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_initialise(elements)
    !=========================================================================!
    ! Initializes the smeared ions:                                           !
    ! - generates the smeared ion density,                                    !
    ! - calculates the self-interaction term (constant),                      !
    ! - calculates the nonself-interactin term (constant until the ions move).!
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   elements(input): The array describing the ions.                       !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    !=========================================================================!

    use ion, only: ELEMENT
    use cell_grid, only: pub_fine_grid
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check

    implicit none

    ! jd: Arguments
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! jd: Internal variables
    integer :: ierr ! jd: Error flag

    !------------------------------------------------------------------------

    call utils_trace_in('smeared_ion_initialise')

    ! jd: Allocate the array that will hold smeared ion density
    allocate(rho_ion(pub_fine_grid%ld1, pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('smeared_ion_initialise','rho_ion',ierr)

    ! jd: Populate this array
    call smeared_ion_generate_density(elements)

    ! jd: Compute constant terms to energy
    smeared_ion_E_self = smeared_ion_calc_E_self(elements)
    ! jd: @note: must be invalidated when ionic positions change!
    smeared_ion_E_smeared = smeared_ion_calc_E_smeared(elements)

    call utils_trace_out('smeared_ion_initialise')

  end subroutine smeared_ion_initialise
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_exit
    !=========================================================================!
    ! Cleans up after the smeared ions:                                       !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    !=========================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! jd: Internal variables
    integer :: ierr ! jd: Error flag

    !------------------------------------------------------------------------

    ! jd: Allocate the array that will hold smeared ion density
    deallocate(rho_ion,stat=ierr)
    call utils_dealloc_check('smeared_ion_exit','rho_ion',ierr)

  end subroutine smeared_ion_exit
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function smeared_ion_calc_E_self(elements)
    !=========================================================================!
    ! Calculates the self-interaction energy of smeared ions.                 !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   elements(input): The array describing the ions.                       !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    !=========================================================================!

    use comms, only: comms_reduce, pub_on_root, pub_my_node_id
    use constants, only: DP, pi, stdout
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_first_atom_on_node, pub_num_atoms_on_node
    use rundat, only: pub_is_smeared_ion_width
    use simulation_cell, only: pub_cell

    implicit none

    ! jd: Arguments
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! jd: Internal variables
    integer :: my_first_at, my_last_at, I      ! jd: Indices and bounds
    real(kind=DP) :: inv_sigma_I               ! jd: Inverse of Gaussian spread
    real(kind=DP) :: Z_I                       ! jd: Core charge
    real(kind=DP) :: E_self                    ! jd: Accumulated quantity

    !------------------------------------------------------------------------

    call utils_trace_in('smeared_ion_calc_E_self')

    inv_sigma_I = 1.0_DP / pub_is_smeared_ion_width

    my_first_at = pub_first_atom_on_node(pub_my_node_id)
    my_last_at = my_first_at + pub_num_atoms_on_node(pub_my_node_id) - 1

    E_self = 0.0_DP

    ! jd: Loop over my atoms
    do I = my_first_at, my_last_at
       Z_I = real(elements(I)%ion_charge, kind=DP)
       E_self = E_self + Z_I*Z_I * inv_sigma_I
    end do

    ! jd: Apply prefactor and reduce
    E_self = E_self / (-sqrt(2.0_DP * pi))
    call comms_reduce('SUM', E_self)

    if(pub_on_root) then
       write(stdout,'(/a,f0.6/)') 'Smeared ion self correction:     ', E_self
    end if

    call utils_trace_out('smeared_ion_calc_E_self')

    smeared_ion_calc_E_self = E_self

  end function smeared_ion_calc_E_self
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function smeared_ion_calc_E_smeared(elements)
    !=========================================================================!
    ! Calculates the non-self-interaction energy of smeared ions.             !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   elements(input): The array describing the ions.                       !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    !=========================================================================!

    use comms, only: comms_reduce, pub_on_root, pub_my_node_id
    use constants, only: DP, stdout
    use geometry, only: POINT, magnitude, OPERATOR(-)
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_first_atom_on_node, pub_num_atoms_on_node
    use rundat, only: pub_is_smeared_ion_width
    use simulation_cell, only: pub_cell
    use utils, only: utils_erf

    implicit none

    ! jd: Arguments
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! jd: Internal variables
    integer :: my_first_at, my_last_at, I, J         ! jd: Indices and bounds
    real(kind=DP) :: sigma_I, sigma_J, sigma_IJ      ! jd: Gaussian spreads
    real(kind=DP) :: Z_I, Z_J                        ! jd: Charges
    real(kind=DP) :: R_IJ                            ! jd: Distance betw. cores
    real(kind=DP) :: E_smeared                       ! jd: Accumulated quantity

    !------------------------------------------------------------------------

    call utils_trace_in('smeared_ion_calc_E_smeared')

    my_first_at = pub_first_atom_on_node(pub_my_node_id)
    my_last_at = my_first_at + pub_num_atoms_on_node(pub_my_node_id) - 1

    E_smeared = 0.0_DP

    ! jd: Loop over my atoms
    do I = my_first_at, my_last_at
       Z_I = real(elements(I)%ion_charge, kind=DP)
       sigma_I = pub_is_smeared_ion_width

       ! jd: Loop over all other atoms
       do J = 1, pub_cell%nat

          if (J == I) cycle

          ! jd: Find charge and distance of other atom
          Z_J = real(elements(J)%ion_charge,kind=DP)
          R_IJ = magnitude(elements(I)%centre - elements(J)%centre)

          ! jd: Calculate sqrt(sigma_i^2 * sigma_j^2)
          sigma_J = pub_is_smeared_ion_width
          sigma_IJ = sqrt(sigma_I*sigma_I + sigma_J*sigma_J)

          ! jd: Add Gaussian smeared charge contribution to energy
          E_smeared = E_smeared + Z_i*Z_J / R_IJ * utils_erf(R_IJ / sigma_IJ)

       end do

    end do

    ! jd: E_smeared has a minus at the front and the sum runs over I<J.
    !     Since I/=J was excluded, we just have to divide by two...
    E_smeared = E_smeared * (-0.5_DP)
    call comms_reduce('SUM', E_smeared)

    if(pub_on_root) then
       write(stdout,'(/a,f0.6)') 'Smeared ion non-self correction: ', E_smeared
    end if

    smeared_ion_calc_E_smeared = E_smeared

    call utils_trace_out('smeared_ion_calc_E_smeared')

  end function smeared_ion_calc_E_smeared
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_apply_vloc_corr(elements, v_loc_on_fine)
    !=========================================================================!
    ! Applies the smeared-ion correction to vloc.                             !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   elements(input):      The array describing the ions.                  !
    !   v_loc_on_fine(inout): Vloc on the fine grid.                          !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: pub_my_node_id
    use constants, only: DP, pi
    use geometry, only: POINT, magnitude, OPERATOR(-), OPERATOR(+), OPERATOR(*)
    use ion, only: ELEMENT
    use rundat, only: pub_is_smeared_ion_width
    use simulation_cell, only: pub_cell
    use utils, only: utils_erf

    implicit none

    ! jd: Arguments
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    real(kind=DP), intent(inout) :: v_loc_on_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12)

    ! jd: Local variables and parameters
    real(kind=DP), parameter :: sqrt_pi = 1.7724538509055_DP
    real(kind=DP) :: Z_I, d, sigma_I, term, accum
    integer :: i1, i2, i3, islab12, I
    type(POINT) :: R_I, r

    !------------------------------------------------------------------------

    call utils_trace_in('smeared_ion_apply_vloc_corr')

    ! jd: Add the correction to v_loc_on_fine
    do islab12 = 1, pub_fine_grid%num_my_slabs12
       i3 = pub_fine_grid%first_slab12(pub_my_node_id) + islab12 - 1

       do i2 = 1, pub_fine_grid%n2
          do i1 = 1, pub_fine_grid%n1

             r = &
                  real((i1-1),kind=DP) * pub_fine_grid%da1 + &
                  real((i2-1),kind=DP) * pub_fine_grid%da2 + &
                  real((i3-1),kind=DP) * pub_fine_grid%da3

             accum = 0.0_DP
             do I = 1, pub_cell%nat

                R_I = elements(I)%centre
                Z_I = real(elements(I)%ion_charge,kind=DP)
                sigma_I = pub_is_smeared_ion_width

                d = magnitude(r-R_I)

                ! jd: Calculate erf(d/sigma_I)/d, considering d==0 separately
                if (d < 1E-10_DP) then
                   ! jd: lim_{d->0} erf(d/sigma_I)/d = 2/(sqrt(pi)*sigma_I)
                   term = 2.0_DP/sigma_I/sqrt_pi
                else
                   term = utils_erf(d/sigma_I) / d
                end if

                ! jd: It's faster to accumulate in a variable first
                accum = accum + Z_I * term

             end do ! over I

             ! jd: Subtract the whole thing from current point in v_loc_fine
             v_loc_on_fine(i1,i2,islab12) = v_loc_on_fine(i1,i2,islab12) + accum

          end do ! over i1
       end do ! over i2
    end do ! over 12-slabs

    call utils_trace_out('smeared_ion_apply_vloc_corr')

  end subroutine smeared_ion_apply_vloc_corr
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_generate_density(elements)
    !=========================================================================!
    ! Generates the smeared-ion density, stores it into rho_ion.              !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   elements(input):      The array describing the ions.                  !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: pub_my_node_id
    use constants, only: pi
    use geometry, only: POINT, magnitude, OPERATOR(-), OPERATOR(+), OPERATOR(*)
    use ion, only: ELEMENT
    use rundat, only: pub_is_smeared_ion_width
    use simulation_cell, only: pub_cell

    implicit none

    ! jd: Arguments
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! jd: Local variables
    real(kind=DP) :: prefactor
    integer :: islab12, i1, i2, i3, I
    type(POINT) :: r, R_I
    real(kind=DP) :: sigma_I, Z_I, d, accum

    !------------------------------------------------------------------------

    call utils_trace_in('smeared_ion_generate_density')

    prefactor = -1.0_DP/pi**1.5_DP

    rho_ion = 0.0_DP ! jd: Take care of padding between pt and ld

    ! jd: Generate smeared ion density on every point of fine grid
    do islab12 = 1, pub_fine_grid%num_my_slabs12
       i3 = pub_fine_grid%first_slab12(pub_my_node_id) + islab12 - 1

       do i2 = 1, pub_fine_grid%n2
          do i1 = 1, pub_fine_grid%n1

             r = &
                  real((i1-1),kind=DP) * pub_fine_grid%da1 + &
                  real((i2-1),kind=DP) * pub_fine_grid%da2 + &
                  real((i3-1),kind=DP) * pub_fine_grid%da3

             accum = 0.0_DP
             do I = 1, pub_cell%nat

                R_I = elements(I)%centre
                Z_I = real(elements(I)%ion_charge,kind=DP)
                sigma_I = pub_is_smeared_ion_width

                d = magnitude(r-R_I)

                ! jd: Calculate density value for d, accumulate
                accum = accum + Z_I / sigma_I**3 * exp(-d*d/(sigma_I*sigma_I))

             end do ! over I

             ! jd: Store
             rho_ion(i1,i2,islab12) = prefactor * accum

          end do ! over i1
       end do ! over i2
    end do ! over 12-slabs

    call utils_trace_out('smeared_ion_generate_density')

  end subroutine smeared_ion_generate_density
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine smeared_ion_hartree(phi, rho_elec, E_Hartree)
    !=========================================================================!
    ! Calculates the Hartree energy in the smeared-ion representation.        !
    ! This is the molecular Hartree energy as described in [1], Appendix.     !
    ! ----------------------------------------------------------------------- !
    ! Arguments:                                                              !
    !   phi(input):      Molecular electrostatic potential on the fine grid.  !
    !   rho_elec(input): *Electronic* charge density on the fine grid.        !
    !                    (smeared-ion density is stored by the module.)       !
    !   E_Hartree(output), optional: Resultant energy.                        !
    ! ----------------------------------------------------------------------- !
    ! Written by Jacek Dziedzic, 05/2010                                      !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: pub_on_root
    use constants, only: stdout
    use multigrid_methods, only: multigrid_calculate_hartree, &
         multigrid_electrostatic_energy, multigrid_prepare_bound_cond
    use integrals, only: integrals_product_on_grid
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    ! jd: Arguments
    ! jd: Electronic density in real space:
    real(kind=DP), intent(inout) :: rho_elec(pub_fine_grid%ld1, &
         pub_fine_grid%ld2,pub_fine_grid%max_slabs12)

    ! jd: Resultant molecular potential, from the MG
    real(kind=DP), intent(out) :: phi(pub_fine_grid%ld1,pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12)

    ! jd: Calculated Hartree energy, if required
    real(kind=DP), intent(out), optional :: E_Hartree

    ! jd: Local variables
    real(kind=DP), allocatable :: rho_tot(:,:,:)
    real(kind=DP) :: E_Hartree_loc

    integer :: ierr ! jd: Error flag

    !------------------------------------------------------------------------

    call utils_trace_in('smeared_ion_hartree')

    ! jd: Allocate the array that will hold smeared ion density
    allocate(rho_tot(pub_fine_grid%ld1,pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('smeared_ion_hartree','rho_tot',ierr)

    ! jd: Prepare a molecular density
    rho_tot = rho_ion + rho_elec

    ! jd: Calculate Hartree potential and energy
    E_Hartree_loc = multigrid_calculate_hartree(phi, rho_tot)
    if(present(E_Hartree)) E_Hartree = E_Hartree_loc

    ! jd: Clean up
    deallocate(rho_tot,stat=ierr)
    call utils_dealloc_check('smeared_ion_hartree','rho_tot',ierr)

    call utils_trace_out('smeared_ion_hartree')

  end subroutine smeared_ion_hartree
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module is_smeared_ions
