! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*- !
!================================================================!
!                                                                !
!             Hubbard DFT+U initialisation module                !
!                                                                !
!----------------------------------------------------------------!
! Written by David O'Regan in April 2009 based on an older       !
! development version of DFT+U begun in February 2008.           !
!----------------------------------------------------------------!
! Extended in November 2011 by Gilberto Teobaldi to accomodate   !
! information for CONSTRAINED_DFT simulations                    !
!================================================================!
! This module contains all information needed for DFT+U          !
! calculations which is specific to each species of Hubbard atom.!
! This includes all information in the HUBBARD block of the      !
! input file, and subsequently species-specific information on   !
! on hydrogenic projectors used for the DFT+U occupancy matrix   !
!================================================================!


module hubbard_init

  use constants

  implicit none

  private

  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ddor: Type for containing DFT+U information for each species
  type HUBBARD_SPECIES

     ! ddor: Name of this Hubbard species in input file
     character(len=4) :: hub_species

     ! ddor: Angular momentum of Hubbard projector. This is used to select 
     !       the most appropriate projectors if self-consistent NGWF projectors are used.   
     integer :: hub_ang_mom

     ! ddor: Hubbard U parameter
     real(kind=DP) :: hub_u

     ! ddor: Effective charge Z/(a_mu / a_0) for given species, atomic units
     real(kind=DP) :: hub_charge

     ! ddor: Hubbard alpha parameter used for determination of the U parameter
     !       consistent with the screening in the system using linear-response.
     real(kind=DP) :: hub_alpha

     ! ddor:  An energy splitting between up and down spin channels which
     !        is used to break local spin-symmetry, say for antiferromagnetic systems.
     !        The potential is lowered by hub_spin_splitting/2 for the 
     !        projection of the spin up density matrix onto the same orbitals used 
     !        as Hubbard projectors, and raised correspondingly for spin down.
     !        Note that this is not a penalty functional for a given expected
     !        local magnetic moment, but may be used to break magnetic symmetry in any
     !        system initially before eg by setting U and alpha to 
     !        zero in an uncorrelated material and selecting those orbitals 
     !        which are most responsible for magnetism as projectors. 
     real(kind=DP) :: hub_spin_splitting

! gibo: === CONSTRAINED_DFT EXTENSIONS: START ============================
     ! Name of this constrained species in input file
     character(len=4) :: cdft_species

     ! Angular momentum of Hubbard projector. This is used to select 
     ! the most appropriate projectors if self-consistent NGWF projectors are used.   
     integer :: cdft_ang_mom

     ! Charge constrain parameter
     real(kind=DP) :: cdft_u_charge, cdft_u_charge_up, cdft_u_charge_down

     ! Spin constrain parameter
     real(kind=DP) :: cdft_u_spin

     ! target value of electrons_up
     real(kind=DP) :: cdft_target_up

     ! target value of electrons_down
     real(kind=DP) :: cdft_target_down

     ! Spin target value
     real(kind=DP) :: cdft_target_spin

     ! logical for donor-atoms
     logical :: cdft_donor

     ! logical for acceptor-atoms
     logical :: cdft_acceptor
! gibo:  === CONSTRAINED_DFT EXTENSIONS: END   ============================

  end type HUBBARD_SPECIES
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  public :: HUBBARD_SPECIES
  public :: hubbard_init_species
  public :: hubbard_init_species_exit

  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  ! ddor: An array num_hub_species long, ordered as the Hubbard species 
  !      are arranged in "block hubbard"
  type(HUBBARD_SPECIES), public, allocatable, dimension(:) :: h_species

  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hubbard_init_species

    !==========================================================! 
    ! This subroutine allocates some memory for                !
    ! the h_species array.                                     !
    !----------------------------------------------------------!
    ! Modified for DFT+U module by David O'Regan in April 2009 !
    !==========================================================!

    use comms, only : pub_on_root,comms_abort
    use simulation_cell, only: pub_cell
    use utils, only : utils_alloc_check
    use rundat, only : pub_output_detail

    implicit none

    ! ddor: Local variables
    integer :: ierr

    ! ddor: Allocate h_species type(hubbard_species) array
    if (.not. allocated(h_species)) then
       allocate(h_species(1:pub_cell%num_hub_species),stat=ierr)
       call utils_alloc_check('hubbard_init_species','h_species',ierr)
    else
       if (pub_on_root) write(stdout,'(a)') &
            'Error in hubbard_init_species: &
            &h_species already allocated'
       call comms_abort
    end if

  end subroutine hubbard_init_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hubbard_init_species_exit

    !==========================================================! 
    ! This subroutine deallocates some memory allocated for    !
    ! the h_species array.                                     !
    !----------------------------------------------------------!
    ! Written for DFT+U module by David O'Regan in April 2009  !
    !==========================================================!

    use comms, only : pub_on_root,comms_abort
    use utils, only : utils_dealloc_check
    use rundat, only : pub_output_detail

    implicit none

    ! ddor: Local variables
    integer :: ierr

    if (allocated(h_species)) then
       deallocate(h_species,stat=ierr)
       call utils_dealloc_check('hubbard_init_species_exit', &
            'h_species',ierr)
    end if

  end subroutine hubbard_init_species_exit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module hubbard_init

