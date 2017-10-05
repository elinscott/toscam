! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                           Hartree module                       !
!----------------------------------------------------------------!
! Originally written by Chris-Kriton Skylaris in June 2000.      !
! Rewritten by Peter Haynes, 28/6/04                             !
! Updated for variable grids by Nicholas Hine, June 2010.        !
!================================================================!

module hartree

  use constants, only: DP
  implicit none

  private

  public :: hartree_on_grid
  public :: hartree_via_multigrid

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hartree_on_grid(hartree_potential,density,grid)

    use cell_grid, only: GRID_INFO
    use constants, only: DP, UP, DN, PI
    use fourier, only: fourier_apply_cell_forward, fourier_apply_cell_backward
    use pbc_corrections, only: pbc_corr_hartree
    use simulation_cell, only: pub_cell
    use rundat, only: pub_mt_cutoff
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
use geometry
use comms, only: comms_reduce, pub_on_root, pub_my_node_id
    implicit none

    ! Arguments
    ! Grid definition
    type(GRID_INFO), intent(inout) :: grid
    ! Density in real space:
    real(kind=DP), intent(inout) :: density(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
    ! Potential in real space:
    real(kind=DP), intent(out) :: hartree_potential(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)

    ! Local variables
    integer :: ierr                                ! Error flag
    real(kind=DP), parameter :: fourpi = 4.0_DP*PI ! 4 pi
    complex(kind=DP), allocatable :: zwork(:,:,:)  ! Workspace

    ! Start timer
    call timer_clock('hartree_on_grid',1)

    ! Allocate workspace
    allocate(zwork(grid%ld3,grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('hartree_on_grid','zwork',ierr)

    if (pub_cell%num_spins == 2) density(:,:,:,1) = &
         density(:,:,:,1) + density(:,:,:,2)

    ! Fourier transform density to reciprocal space
    call fourier_apply_cell_forward(density(:,:,:,1),zwork,grid)

    ! jd: Apply the Martyna-Tuckerman correction, if requested
    if(pub_mt_cutoff == 0.0_DP) then
      ! Multiply by 4pi/g^2 to obtain the usual Hartree potential
      zwork = zwork * grid%coulomb_recip(:,:,:) * fourpi
    else
      ! jd: Use the PBC-corrected formula
      call pbc_corr_hartree(zwork,grid)
    end if

    ! Fourier transform potential to real space
    call fourier_apply_cell_backward(hartree_potential(:,:,:,1),zwork,grid)

    if (pub_cell%num_spins == 2) then
       hartree_potential(:,:,:,DN) = hartree_potential(:,:,:,UP)
       density(:,:,:,1) = density(:,:,:,1) - density(:,:,:,2)
    end if

    ! Deallocate workspace
    deallocate(zwork,stat=ierr)
    call utils_dealloc_check('hartree_on_grid','zwork',ierr)

    ! Stop timer
    call timer_clock('hartree_on_grid',2)

  end subroutine hartree_on_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hartree_via_multigrid(phi, rho_elec, E_Hartree, E_cav)
    !=========================================================================!
    ! This is a front-end subroutine for the calculation of the Hartree energy!
    ! and Hartree potential using a multigrid solver. There are three,        !
    ! mutually-exclusive cases when this is in order:                         !
    !   a) An open BC calculation, without smeared ions, in vacuum.           !
    !   b) An open BC calculation, with smeared ions, in vacuum.              !
    !   c) An open BC calculation, with smeared ions, in solvent.             !
    ! In a) the Hartree energy and potential are only due to electrons.       !
    ! In b) and c) the Hartree energy and potential are due to the whole      !
    !       molecule, i.e. electrons and smeared ions.                        !
    ! With b) and c), because of its molecular nature, E_Hartree cannot be    !
    ! calculated in the usual fashion, by integrating the Hartree potential   !
    ! with the electronic density, hence this subroutine calculates it on its !
    ! own. Futhermore, in the subcase of c) with a self-consistently changing !
    ! dielectric cavity, extra terms appear in the gradient, which are most   !
    ! conveniently treated by adding them to the potential. This subroutine   !
    ! does this.                                                              !
    ! For convenience, in c) this subroutine also calculates the cavitation   !
    ! energy, because it's easiest to do this at the same time. In cases a)   !
    ! and b) this term will be zero.                                          !
    ! Because the multigrid only works with the fine grid, there is no need   !
    ! to pass a grid argument to this subroutine.                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   phi,              intent=out, the calculated Hartree potential        !
    !   rho_elec,         intent=inout,  the input _electronic_ density       !
    !                     NB: inout, as it's temporarily summed over spins    !
    !   E_Hartree, (opt)  intent=out, calculated Hartree energy, if requested !
    !   E_cav,     (opt)  intent=out, cavitation energy, if requested         !
    ! Q. Why is E_cav calculated here?                                        !
    ! A. Because current_problem lives only for the duration of this routine. !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 10/12/2010.                                  !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use constants, only: DP, DN, UP
    use is_smeared_ions, only: smeared_ion_hartree
    use is_solvation, only: implicit_solvent_hartree
    use multigrid_methods, only: multigrid_calculate_hartree, &
         multigrid_initialise
    use rundat, only: pub_is_smeared_ion_rep, pub_is_implicit_solvent, &
         pub_multigrid_hartree
    use simulation_cell, only: pub_cell
    use utils, only: utils_abort

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)   :: phi(pub_fine_grid%ld1, pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(inout) :: rho_elec(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(out), optional :: E_Hartree
    real(kind=DP), intent(out), optional :: E_cav

    ! jd: Local variables
    real(kind=DP) :: E_Hartree_loc

    ! ------------------------------------------------------------------------

    if(present(E_Hartree)) E_Hartree = 0D0
    if(present(E_cav)) E_cav = 0D0

    ! jd: Use total electronic density, if we have two spins
    if (pub_cell%num_spins == 2) rho_elec(:,:,:,UP) = &
         rho_elec(:,:,:,UP) + rho_elec(:,:,:,DN)

    ! jd: Need to initialise multigrid first, because the hartree routines
    !     dimension arrays with variables that are set up in the initialiser
    call multigrid_initialise

    ! jd: Case c)
    if(pub_is_implicit_solvent) then
       call implicit_solvent_hartree(phi, rho_elec(:,:,:,UP), E_Hartree, E_cav)

    ! jd: Case b)
    else if(pub_is_smeared_ion_rep) then
       call smeared_ion_hartree(phi, rho_elec(:,:,:,UP), E_Hartree)

    ! jd: Case a)
    else if(pub_multigrid_hartree) then
       E_Hartree_loc = multigrid_calculate_hartree(phi, rho_elec(:,:,:,UP))
       if(present(E_Hartree)) E_Hartree = E_Hartree_loc
    else
       call utils_abort('Internal error in hartree_via_multigrid.')
    end if

    ! jd: Undo the spin summation, update second spin component of potential
    if (pub_cell%num_spins == 2) then
       phi(:,:,:,DN) = phi(:,:,:,UP)
       rho_elec(:,:,:,UP) = rho_elec(:,:,:,UP) - rho_elec(:,:,:,DN)
    end if

  end subroutine hartree_via_multigrid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module hartree

