! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*- !
!================================================================!
!                                                                !
!             Hubbard DFT+U build module                         !
!                                                                !
!----------------------------------------------------------------!
! Written by David O'Regan in April 2009 and maintained          !
! subsequently by David O'Regan and Nicholas Hine.               !                                                  !
! Re-worked to use the new HUBBARD_MODEL type by Nicholas Hine   !
! in 2011.                                                       !
! Addition of Constrained DFT code by Gilberto Teobaldi, 2011.   !
!================================================================!
! This module contains all information needed for DFT+U and cDFT !
! calculations.                                                  !
!                                                                !
! It defines the HUBBARD_MODEL type which contains all data      !
! relating to an instance of the Hubbard model, defining the     !
! correlated sites, the projectors used to define the correlated !
! manifold, the overlaps between NGWFs and Hubbard projectors    !
! of both hydrogenic and self-consistent NGWF projector type,    !
! the tensorial correction matrix needed to compute the          !
! occupation matrix when non-orthogonal projectors are used, and !
! the Hubbard Hamiltonian for the system in the projector        !
! representation.                                                !
!                                                                !
! Computation of the DFT+U energy and Hamiltonian is             !
! carried out in this module, as is mixing of the projectors     !
! to obtain projector-self-consistency.                          !
!                                                                !
! For details, see the following papers:                         !
!                                                                !
!  Linear-scaling DFT+U with full local orbital optimization     !
!  D. D. O'Regan, N. D. M. Hine, M. C. Payne and A. A. Mostofi,  !
!  Phys. Rev. B 85, 085107 (2012).                               !
!                                                                !
!  Optimised Projections for the Ab Initio Simulation of Large   !
!  and Strongly Correlated Systems,                              !
!  D. D. O'Regan, (Springer, Berlin, Heidelberg, 2012) 1st Ed.,  !
!  Springer Theses XVI, 225p., ISBN 978-3-642-23237-4.           !
!                                                                !
!  Subspace representations in ab initio methods for strongly    !
!  correlated systems,                                           !
!  D. D. O'Regan, M. C. Payne and A. A. Mostofi,                 !
!  Phys. Rev. B 83, 245124 (2011)                                !
!                                                                !
!  Projector self-consistent DFT+U using non-orthogonal          !
!  generalized Wannier functions                                 !
!  D. D. O'Regan, N. D. M. Hine, M. C. Payne and A. A. Mostofi,  !
!  Phys. Rev. B 82, 081102(R) (2010).                            !
!                                                                !
!================================================================!

module hubbard_build

!CW
#ifdef GPU_SPEEDUP
  use fortran_cuda
#endif
  use mpi
!END CW

  use constants, only: DP, stdout, VERBOSE
  use sparse, only: SPAM3
  use projectors, only: PROJECTOR_SET

  implicit none

  private

!CW
    TYPE array_of_matrices
      real(8),allocatable :: occupancy(:,:,:)       
     end TYPE 
!END CW


  ! ndmh: "Container" type to store all information relating to an instance
  ! ndmh: of a Hubbard model projected from DFT+U orbitals
  type HUBBARD_MODEL

    ! Original atom number of each Hubbard atom
    integer, allocatable       :: orig(:)

    ! The hubbard_species index of each Hubbard atom
    integer, allocatable       :: species_number(:)

    ! A list of the NGWFs with greatest overlap with localised
    ! hydrogenic orbitals on this atom. Ordered according to
    ! the selection criterion.
    integer, allocatable       :: localisedngwfs(:,:) !(2l+1,nat_hub)

    ! ddor: The total occupancy for each site and spin
    real(kind=DP), allocatable :: population(:,:)

    ! ddor: The DFT+U energy for each site and spin
    real(kind=DP), allocatable :: energy(:,:)

    ! ddor: The total projection of the current set of NGWF projectors onto
    !       hydrogenic orbitals for each site
    real(kind=DP), allocatable :: current_projection(:)

    ! gibo: the cDFT gradient component for each site and spin
    real(kind=DP), allocatable :: cdft_gradient(:,:)

    !gibo: complementary arrays for atom-resolved (not species-) cDFT U-pot.
    real(kind=DP), allocatable :: cdft_u_charge_up(:)
    real(kind=DP), allocatable :: cdft_u_charge_down(:)
    real(kind=DP), allocatable :: cdft_u_spin(:)

    ! ndmh: projector set type
    type(PROJECTOR_SET)        :: projectors

    ! ddor: arrays to store NGWFs, projectors and energies for projector
    !       self-consistency
    real(kind=DP), allocatable :: consistency_ngwfs(:)
    real(kind=DP), allocatable :: consistency_projs(:)
    real(kind=DP), allocatable :: consistency_energies(:)
    real(kind=DP), allocatable :: consistency_projections(:)

    type(SPAM3) :: o_matrix    !ddor: Hubbard projector overlap matrix
    type(SPAM3) :: u_matrix    !ddor: Hubbard U, alpha
    type(SPAM3) :: up_matrix   !      and spin-splitting parameters
    type(SPAM3) :: down_matrix !      in block-diagonal form
    type(SPAM3), allocatable :: occupancy_matrix(:)
    type(SPAM3), allocatable :: projector_ham(:)

    ! ddor: The final NGWF overlap matrix of a Hubbard projector
    !       iteration is needed for the next
    !       in order to construct the covariant metric on the correlated sites
    type(SPAM3) :: consistency_overlap

    ! ddor: The current iteration number of  the Hubbard projector
    !       self-consistency cycle.
    integer :: consistency_iteration

    ! ddor: True if we wish to carry out DFT+U with hydrogenic projectors
    !       on the first HUBBARDSCF iteration, false for plain DFT.
    logical :: dftu_on_first_hubbardscf

    ! ddor: Make a copy of maxit_ngwf_cg for restart of HUBBARDSCF
    integer :: store_maxit_ngwf_cg

  end type HUBBARD_MODEL

  public :: HUBBARD_MODEL

  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  public :: hubbard_model_init
  public :: hubbard_model_exit
  public :: hubbard_build_matrices
  public :: hubbard_build_matrices_exit
  public :: cdft_build_matrices
  public :: hubbard_build_consist
  public :: hubbard_build_consist_exit
  public :: hubbard_spin_splitting_zero
  public :: hubbard_ham_matrix
  public :: hubbard_projector_ham
  public :: hubbard_energy_total
  public :: hubbard_energy_info
  public :: cdft_energy_total
  public :: cdft_energy_info
  public :: hubbard_projection_mtx
  public :: hubbard_projector_consistency
  public :: hubbard_test_convergence
  public :: hubbard_species_proj
  public :: hubbard_species_exit_proj
  public :: hubbard_calculate_forces
!CW
  public :: hubbard_dmft_interface
!END CW

  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hubbard_model_init(hub,elements)

    !==========================================================!
    ! This subroutine allocates memory for the HUBBARD_MODEL   !
    ! type.                                                    !
    !----------------------------------------------------------!
    ! Written by Nicholas Hine in November 2011, based on code !
    ! by David O'Regan written in April 2009                   !
    !==========================================================!

    use comms, only: pub_on_root,comms_abort
    use hubbard_init, only: h_species
    use ion, only: ELEMENT
    use rundat, only: pub_cdft, &
        pub_cdft_atom_charge, pub_cdft_atom_spin, &
        pub_cdft_group_charge, pub_cdft_group_spin, &
        pub_cdft_group_charge_diff, pub_cdft_group_spin_diff
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! Local Variables
    integer :: ierr
    integer :: dummy_h, iat_orig
    integer :: channels, sp

    allocate(hub%orig(pub_cell%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_model_init','hub%orig',ierr)
    allocate(hub%species_number(pub_cell%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_model_init','hub%species_number',ierr)

    ! Set up orig and species_number arrays
    dummy_h = 0
    do sp=1,pub_cell%num_hub_species
       do iat_orig=1,pub_cell%nat ! The total number of atoms in the cell
          if (elements(iat_orig)%species_id == h_species(sp)%hub_species) then
             dummy_h = dummy_h + 1
             ! ddor: Original atom number of Hubbard atom
             hub%orig(dummy_h) = iat_orig
             hub%species_number(dummy_h) = sp
          end if
       end do
    end do

    if (dummy_h .ne. pub_cell%nat_hub) then
       if (pub_on_root) write(stdout,*) &
            'Error in hubbard_model_init: final dummy_h', &
            dummy_h,'does not match pub_cell%nat_hub',pub_cell%nat_hub
       call comms_abort
    end if

    channels = maxval(2*h_species(:)%hub_ang_mom + 1)
    allocate(hub%population(pub_cell%num_spins,pub_cell%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_model_init','hub%population',ierr)
    allocate(hub%energy(pub_cell%num_spins,pub_cell%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_model_init','hub%energy',ierr)
    allocate(hub%localisedngwfs(channels,pub_cell%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_model_init','hub%localisedngwfs',ierr)
    allocate(hub%current_projection(pub_cell%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_model_init','hub%current_projection',ierr)
    ! gibo: for cDFT, allocate cDFT-gradient components
    if (pub_cdft) then
       allocate(hub%cdft_gradient(pub_cell%num_spins,pub_cell%nat_hub), &
            stat=ierr)
       call utils_alloc_check('hubbard_model_init','hub%cdft_gradient',ierr)

       ! and atom specific potentials/target arrays
       if (pub_cdft_atom_charge .OR. pub_cdft_group_charge .OR. &
            pub_cdft_group_charge_diff) then

          allocate(hub%cdft_u_charge_up(pub_cell%nat_hub),stat=ierr)
          call utils_alloc_check('hubbard_model_init',&
               'hub%cdft_u_charge_up',ierr)
          allocate(hub%cdft_u_charge_down(pub_cell%nat_hub),stat=ierr)
          call utils_alloc_check('hubbard_model_init',&
               'hub%cdft_u_charge_down',ierr)

          ! initialise atom-specific U-potentials from species-specific input
          ! Why? To be able to optimise different cDFT U-potentials for cDFT-atoms
          ! of the same species (thence targeted population) in different
          ! chemical environments...
          dummy_h = 0
           do sp=1,pub_cell%num_hub_species
             do iat_orig=1,pub_cell%nat ! The total number of atoms in the cell
                if (elements(iat_orig)%species_id &
                        == h_species(sp)%hub_species) then
                    dummy_h = dummy_h + 1
                   hub%cdft_u_charge_up(dummy_h)   =  &
                        h_species(sp)%cdft_u_charge_up
                   hub%cdft_u_charge_down(dummy_h) =  &
                        h_species(sp)%cdft_u_charge_down
                end if
             end do
          end do

       else if (pub_cdft_atom_spin .OR. pub_cdft_group_spin .OR. &
               pub_cdft_group_spin_diff) then

          allocate(hub%cdft_u_spin(pub_cell%nat_hub),stat=ierr)
          call utils_alloc_check('hubbard_model_init',&
               'hub%cdft_u_spin',ierr)

          ! initialise atom-specific U-potentials from species-specific input
          ! Why? To be able to optimise different cDFT U-potentials for cDFT-atoms
          ! of the same species (thence targeted population) in different
          ! chemical environments...
          dummy_h = 0
          do sp=1,pub_cell%num_hub_species
             do iat_orig=1,pub_cell%nat ! The total number of atoms in the cell
                if (elements(iat_orig)%species_id == h_species(sp)%hub_species) then
                   dummy_h = dummy_h + 1
                   hub%cdft_u_spin(dummy_h) =  h_species(sp)%cdft_u_spin
                end if
             end do
          end do

       end if
    end if

  end subroutine hubbard_model_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_model_exit(hub)

    !==========================================================!
    ! This subroutine deallocates memory for the HUBBARD_MODEL !
    ! type.                                                    !
    !----------------------------------------------------------!
    ! Written by Nicholas Hine in November 2011, based on code !
    ! by David O'Regan written in April 2009                   !
    !==========================================================!

    use rundat, only: pub_cdft,pub_cdft_atom_charge, pub_cdft_atom_spin, &
         pub_cdft_group_charge, pub_cdft_group_spin, &
         pub_cdft_group_charge_diff, pub_cdft_group_spin_diff        
    use simulation_cell, only: pub_cell
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub

    ! Local Variables
    integer :: ierr

    ! Deallocate allocatable arrays in HUBBARD_MODEL
    if (pub_cdft) then
       deallocate(hub%cdft_gradient,stat=ierr)
       call utils_dealloc_check('hubbard_model_exit','hub%cdft_gradient',ierr)

       ! and atom specific potentials/target arrays
       if (pub_cdft_atom_charge .OR. pub_cdft_group_charge .OR. &
           pub_cdft_group_charge_diff) then

         deallocate(hub%cdft_u_charge_up,stat=ierr)
         call utils_dealloc_check('hubbard_model_exit',&
              'hub%cdft_u_charge_up',ierr)
         deallocate(hub%cdft_u_charge_down,stat=ierr)
         call utils_dealloc_check('hubbard_model_exit',&
              'hub%cdft_u_charge_down',ierr)

       else if (pub_cdft_atom_spin .OR. pub_cdft_group_spin .OR. &
               pub_cdft_group_spin_diff) then

          deallocate(hub%cdft_u_spin,stat=ierr)
          call utils_dealloc_check('hubbard_model_exit',&
               'hub%cdft_u_spin',ierr)
       end if
    end if
    deallocate(hub%current_projection,stat=ierr)
    call utils_dealloc_check('hubbard_model_exit','hub%current_projection',ierr)
    deallocate(hub%localisedngwfs,stat=ierr)
    call utils_dealloc_check('hubbard_model_exit','hub%localisedngwfs',ierr)
    deallocate(hub%energy,stat=ierr)
    call utils_dealloc_check('hubbard_model_exit','hub%energy',ierr)
    deallocate(hub%population,stat=ierr)
    call utils_dealloc_check('hubbard_model_exit','hub%population',ierr)
    deallocate(hub%species_number,stat=ierr)
    call utils_dealloc_check('hubbard_model_init','hub%species_number',ierr)
    deallocate(hub%orig,stat=ierr)
    call utils_dealloc_check('hubbard_model_init','hub%orig',ierr)

  end subroutine hubbard_model_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_build_matrices(hub,hub_proj_basis)

    !==========================================================!
    ! This subroutine creates DFT+U SPAM3 matrices             !
    !----------------------------------------------------------!
    ! Written for DFT+U module by David O'Regan in August 2009 !
    ! originally as hubbard_build_spam.                        !
    !==========================================================!

    use comms, only: pub_my_node_id, pub_on_root
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only:  pub_num_hub_atoms_on_node, &
         pub_hub_atoms_on_node, pub_distr_atom, pub_first_atom_on_node
    use rundat, only: pub_aug, pub_cdft, pub_cdft_hubbard, &
         task, pub_hubbard_restart, pub_hubbard_atomsolve
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, &
         sparse_scale, sparse_put_element
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    integer :: sp, hat_on_node, hat, theatom
    integer :: hub_proj
    integer :: is, ierr
    real(kind=DP) :: half_u, potential_up, potential_down

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_build_matrices'
#endif

    ! ndmh: Create projector overlap matrix
    hub%o_matrix%structure = 'G'
    call sparse_create(hub%o_matrix)

    ! ndmh: Create U and alpha matrices
    hub%u_matrix%structure = 'G'
    call sparse_create(hub%u_matrix)
    call sparse_create(hub%up_matrix,hub%u_matrix)
    call sparse_create(hub%down_matrix,hub%u_matrix)

    allocate(hub%occupancy_matrix(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('hubbard_build_matrices','occupancy_matrix',ierr)
    allocate(hub%projector_ham(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('hubbard_build_matrices','projector_ham',ierr)
    do is=1,pub_cell%num_spins
       hub%occupancy_matrix(is)%structure = 'G'
       call sparse_create(hub%occupancy_matrix(is))
       hub%projector_ham(is)%structure = 'G'
       call sparse_create(hub%projector_ham(is))
    end do

    ! gibo: Hubbard potentials only for pure Hubbard-only .OR. cdft_hubbard simulations
    if ((.not.pub_cdft) .or. (pub_cdft_hubbard)) then

       ! ddor: Loop over Hubbard atoms on my node
       do hat_on_node=1,pub_num_hub_atoms_on_node(pub_my_node_id)

          hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
          sp = hub%species_number(hat)
          theatom = pub_distr_atom(hub%orig(hat))

          half_u = h_species(sp)%hub_u * 0.5_DP
          potential_up = h_species(sp)%hub_alpha - &
               ( h_species(sp)%hub_spin_splitting * 0.5_DP )
          potential_down = h_species(sp)%hub_alpha + &
               ( h_species(sp)%hub_spin_splitting * 0.5_DP )

          ! ddor: Loop over Hubbard projectors on my node
          do hub_proj = hub_proj_basis%first_on_atom(theatom), &
               hub_proj_basis%first_on_atom(theatom) + &
               hub_proj_basis%num_on_atom(theatom) - 1

             call sparse_put_element(half_u,&
                  hub%u_matrix,hub_proj,hub_proj)
             call sparse_put_element(potential_up,&
                  hub%up_matrix,hub_proj,hub_proj)
             call sparse_put_element(potential_down,&
                  hub%down_matrix,hub_proj,hub_proj)

          end do

       end do

    end if

    ! gibo: for constrained-DFT create and add cDFT terms
    if (pub_cdft) call cdft_build_matrices(hub,hub_proj_basis)

    if ((task == 'HUBBARDSCF') .or. pub_hubbard_restart &
         & .or. pub_hubbard_atomsolve) then
       hub%consistency_overlap%structure = 'S'
       call sparse_create(hub%consistency_overlap)
    endif

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_build_matrices'
#endif

  end subroutine hubbard_build_matrices



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_build_matrices_exit(hub)

    !==========================================================!
    ! This subroutine destroys DFT+U SPAM3 matrices            !
    !----------------------------------------------------------!
    ! Written for DFT+U module by David O'Regan in August 2009 !
    ! originally as hubbard_build_spam_exit.                   !
    !==========================================================!

    use rundat, only: task, pub_hubbard_restart, pub_hubbard_atomsolve
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_destroy
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub

    ! Local Variables
    integer :: ierr
    integer :: is

    if (task == 'HUBBARDSCF' .or. pub_hubbard_restart &
         & .or. pub_hubbard_atomsolve) then
       call sparse_destroy(hub%consistency_overlap)
       call hubbard_build_consist_exit(hub)
    endif

    call sparse_destroy(hub%down_matrix)
    call sparse_destroy(hub%up_matrix)
    call sparse_destroy(hub%u_matrix)
    call sparse_destroy(hub%o_matrix)
    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(hub%occupancy_matrix(is))
       call sparse_destroy(hub%projector_ham(is))
    end do
    deallocate(hub%projector_ham,stat=ierr)
    call utils_dealloc_check('hubbard_build_matrices_exit', &
         'projector_ham',ierr)

  end subroutine hubbard_build_matrices_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cdft_build_matrices(hub,hub_proj_basis)

    !==========================================================!
    ! This subroutine creates CDFT SPAM3 matrices              !
    !----------------------------------------------------------!
    ! Adapted to CDFT from hubbard_build_spam in November 2011 !
    ! by Gilberto Teobaldi.                                    !
    !==========================================================!

    use comms, only: pub_my_node_id, pub_on_root
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only:  pub_num_hub_atoms_on_node, &
         pub_hub_atoms_on_node, pub_distr_atom, pub_first_atom_on_node
    use rundat, only: pub_aug, pub_cdft, pub_cdft_hubbard, &
          pub_cdft_atom_charge, pub_cdft_atom_spin, &
          pub_cdft_group_charge, pub_cdft_group_spin, &
          pub_cdft_group_charge_diff, pub_cdft_group_spin_diff
    use sparse, only: SPAM3, sparse_create, sparse_destroy,&
         sparse_scale, sparse_put_element, sparse_axpy

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    type(SPAM3) :: cdft_potential_up
    type(SPAM3) :: cdft_potential_down
    integer :: sp, hat_on_node, hat, theatom
    integer :: hub_proj
    real(kind=DP) :: u_charge_up, u_charge_down
    real(kind=DP) :: u_spin_up, u_spin_down

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &cdft_build_matrices'
#endif

    ! gibo: create local spam matrixes for CDFT UP/DOWN potential
    call sparse_create(cdft_potential_up,   hub%up_matrix)
    call sparse_create(cdft_potential_down, hub%down_matrix)


!====select CDFT-mode and iterate...
    ! gibo: ATOM/GROUP-CHARGE-[DIFF] constrained run
    CDFT_MODE: if ((pub_cdft_atom_charge .OR. pub_cdft_group_charge).OR.&
                   pub_cdft_group_charge_diff) then
!====select CDFT-mode and iterate...

    ! ddor: Loop over Hubbard atoms on my node
    cdft_atoms_1 : do hat_on_node = 1, &
         pub_num_hub_atoms_on_node(pub_my_node_id)

       hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
       sp = hub%species_number(hat)
       theatom = pub_distr_atom(hub%orig(hat))

       u_charge_up   = hub%cdft_u_charge_up(hat)
       u_charge_down = hub%cdft_u_charge_down(hat)

       ! ddor: Loop over Hubbard projectors on my node
       hubprojs_cdft_1: do hub_proj = hub_proj_basis%first_on_atom(theatom), &
            hub_proj_basis%first_on_atom(theatom) + &
            hub_proj_basis%num_on_atom(theatom) - 1

          call sparse_put_element(u_charge_up,&
               cdft_potential_up,hub_proj,hub_proj)
          call sparse_put_element(u_charge_down,&
               cdft_potential_down,hub_proj,hub_proj)

       end do hubprojs_cdft_1

    end do cdft_atoms_1


!====select CDFT-mode and iterate...
    ! gibo: ATOM/GROUP-SPIN-[DIFF] constrained run
    else if ((pub_cdft_atom_spin .OR. pub_cdft_group_spin).OR.&
                   pub_cdft_group_spin_diff) then
!====select CDFT-mode and iterate...

    ! ddor: Loop over Hubbard atoms on my node
    cdft_atoms_2 : do hat_on_node = 1, &
         pub_num_hub_atoms_on_node(pub_my_node_id)

       hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
       sp = hub%species_number(hat)
       theatom = pub_distr_atom(hub%orig(hat))

       ! MIND THE SIGNS: Ecdft_spin = Us*[Nup -Ndown - DN_target]
       ! thus V_up = Us, and V_down = -Us
       u_spin_up   =  hub%cdft_u_spin(hat)
       u_spin_down = -hub%cdft_u_spin(hat)

       ! ddor: Loop over Hubbard projectors on my node
       hubprojs_cdft_2: do hub_proj = hub_proj_basis%first_on_atom(theatom), &
            hub_proj_basis%first_on_atom(theatom) + &
            hub_proj_basis%num_on_atom(theatom) - 1

          call sparse_put_element(u_spin_up,&
               cdft_potential_up,hub_proj,hub_proj)
          call sparse_put_element(u_spin_down,&
               cdft_potential_down,hub_proj,hub_proj)

       end do hubprojs_cdft_2

    end do cdft_atoms_2

!====select CDFT-mode and iterate...
    end if CDFT_MODE
!====select CDFT-mode and iterate...

    ! gibo: add cdft_potential_UP/DOWN to hub_UP/DOWN_matrix
    call sparse_axpy(hub%up_matrix,   cdft_potential_up,   1.0_DP)
    call sparse_axpy(hub%down_matrix, cdft_potential_down, 1.0_DP)

    ! gibo: destroy local spam matrixes
    call sparse_destroy(cdft_potential_up)
    call sparse_destroy(cdft_potential_down)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &cdft_build_matrices'
#endif

  end subroutine cdft_build_matrices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_build_consist(hub,hub_proj_basis,ngwf_basis)

    !==========================================================!
    ! This subroutine allocates arrays necessary for carrying  !
    ! out self-consistency over Hubbard projectors             !
    !----------------------------------------------------------!
    ! Written for DFT+U module by David O'Regan in Oct 2009    !
    !==========================================================!

    use function_basis, only: FUNC_BASIS
    use utils, only: utils_alloc_check
    use simulation_cell, only: pub_cell
    use rundat, only: pub_hub_conv_win

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    integer :: ierr

    ! ddor: initialise matrices used for DFT+U
    !       with self-consistent projectors
    allocate(hub%consistency_energies(pub_hub_conv_win),stat=ierr)
    call utils_alloc_check('hubbard_build_consist', &
         'hub%consistency_energies',ierr)
    hub%consistency_energies = 0.0_DP
    allocate(hub%consistency_projections(pub_hub_conv_win),stat=ierr)
    call utils_alloc_check('hubbard_build_consist', &
         'hub%consistency_projections',ierr)
    hub%consistency_projections = 0.0_DP
    allocate(hub%consistency_ngwfs(&
         &ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('hubbard_build_consist', &
         'hub%consistency_ngwfs',ierr)
    allocate(hub%consistency_projs(&
         &hub_proj_basis%n_ppds*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('hubbard_build_consist', &
         'hub%consistency_projs',ierr)

    hub%consistency_iteration = 1
    hub%dftu_on_first_hubbardscf = .true.

  end subroutine hubbard_build_consist


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_build_consist_exit(hub)

    !==========================================================!
    ! This subroutine deallocates arrays for carrying          !
    ! out self-consistency over Hubbard projectors             !
    !----------------------------------------------------------!
    ! Written for DFT+U module by David O'Regan in Oct 2009    !
    !==========================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub

    ! Local Variables
    integer :: ierr

    ! ddor: Deallocate arrays used for DFT+U with self-consistent projectors
    deallocate(hub%consistency_projs,stat=ierr)
    call utils_dealloc_check('hubbard_build_consist_exit', &
         'hub%consistency_projs',ierr)
    deallocate(hub%consistency_ngwfs,stat=ierr)
    call utils_dealloc_check('hubbard_build_consist_exit', &
         'hub%consistency_ngwfs',ierr)
    deallocate(hub%consistency_projections,stat=ierr)
    call utils_dealloc_check('hubbard_build_consist_exit', &
         'hub%consistency_projections',ierr)
    deallocate(hub%consistency_energies,stat=ierr)
    call utils_dealloc_check('hubbard_build_consist_exit', &
         'hub%consistency_energies',ierr)

  end subroutine hubbard_build_consist_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_spin_splitting_zero(hub,hub_proj_basis)

    !==========================================================!
    ! This subroutine sets any spin-splitting back to zero.    !
    !----------------------------------------------------------!
    ! Written for DFT+U module by David O'Regan in August 2009 !
    !==========================================================!

    use comms, only: pub_my_node_id, pub_on_root
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: pub_num_hub_atoms_on_node, &
         pub_hub_atoms_on_node, pub_distr_atom
    use sparse, only: SPAM3, sparse_create, &
         sparse_scale, sparse_put_element

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    integer :: sp, hat_on_node, hat, theatom
    integer :: hub_proj
    real(kind=DP) :: potential_up, potential_down

    if (SUM(ABS(h_species(:)%hub_spin_splitting)) .gt. 1.0e-10_DP) then

       if (pub_on_root) write(stdout,*) &
            'WARNING : Removing DFT+U spin-splitting potential.'

       call sparse_scale(hub%up_matrix,0.0_DP)
       call sparse_scale(hub%down_matrix,0.0_DP)

       h_species(:)%hub_spin_splitting = 0.0_DP

       ! ddor: Loop over Hubbard atoms on my node
       do hat_on_node=1,pub_num_hub_atoms_on_node(pub_my_node_id)

          hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
          sp = hub%species_number(hat)
          theatom = pub_distr_atom(hub%orig(hat))

          potential_up = h_species(sp)%hub_alpha
          potential_down = h_species(sp)%hub_alpha

          ! ddor: Loop over Hubbard projectors on my node
          do hub_proj = hub_proj_basis%first_on_atom(theatom), &
               hub_proj_basis%first_on_atom(theatom) + &
               hub_proj_basis%num_on_atom(theatom) - 1

             call sparse_put_element(potential_up,&
                  hub%up_matrix,hub_proj,hub_proj)
             call sparse_put_element(potential_down,&
                  hub%down_matrix,hub_proj,hub_proj)

          end do

       end do

    end if

  end subroutine hubbard_spin_splitting_zero


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_ham_matrix(hub,hubbard_ham,hub_overlap,hub_overlap_t)

    !==========================================================!
    ! This subroutine calculates the Hubbard DFT+U Hamiltonian !
    ! in SPAM3 form, including the components due to the alpha !
    ! and spin-splitting parameters, used for determining the  !
    ! value of U which is  consistent with the screening in    !
    ! the system and breaking magnetic symmetry, respectively. !
    !----------------------------------------------------------!
    ! Written by David O'Regan in August 2009                  !
    ! Modified by Nicholas Hine in April 2011 to separate      !
    ! creation of projector Ham from creation of Ham           !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: pub_on_root
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(in) :: hub
    type(SPAM3), intent(inout) :: hubbard_ham(pub_cell%num_spins)
    type(SPAM3), intent(in) :: hub_overlap
    type(SPAM3), intent(in) :: hub_overlap_t

    ! Local Variables
    type(SPAM3) :: v_hamiltonian_buffer
    integer :: is

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_ham_matrix'
#endif

    call timer_clock('hubbard_ham_matrix',1)

    call sparse_create(v_hamiltonian_buffer,hub_overlap)                !VG=V

    do is = 1, pub_cell%num_spins

       ! ddor: V [ U/2 (1-2N) + a]  blanked to 'VG'
       call sparse_product(v_hamiltonian_buffer,&
            hub_overlap,hub%projector_ham(is))
       ! ddor: Hamiltonian
       call sparse_product(hubbard_ham(is),&
            v_hamiltonian_buffer,hub_overlap_t)

    end do

    call sparse_destroy(v_hamiltonian_buffer)

    call timer_clock('hubbard_ham_matrix',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_ham_matrix'
#endif

  end subroutine hubbard_ham_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_projector_ham(hub,kernel,hub_overlap, &
       hub_overlap_t)

    !==========================================================!
    ! Calculates the DFT+U Hamiltonian in the Hubbard projector!
    ! basis for a given spin channel and for the quadratic     !
    ! functional.                                              !
    !----------------------------------------------------------!
    ! Written by David O'Regan in September 2009.              !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: comms_abort, pub_on_root
    use hubbard_init, only: h_species
    use rundat, only: pub_hub_functional, pub_cdft, pub_cdft_hubbard
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, &
         sparse_product, sparse_scale, sparse_axpy, sparse_copy
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(SPAM3), intent(in) :: kernel(pub_cell%num_spins)
    type(SPAM3), intent(in) :: hub_overlap
    type(SPAM3), intent(in) :: hub_overlap_t

    ! Local Variables
    type(SPAM3) :: site_buffer
    type(SPAM3) :: occupancy_buffer, occupancy_buffer2, occupancy_buffer3
    integer :: is

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_projector_ham'
#endif

    call timer_clock('hubbard_projector_ham',1)

    ! gibo: calculate Hubbard terms only for pure Hubbard-only
    !       .OR. cdft_hubbard simulations
    ! *** WARNING: this if-construct is closed on line 996
    HUBBARD_CDFT: if ((.not.pub_cdft) .OR. (pub_cdft_hubbard) ) then

    call sparse_create(site_buffer,hub_overlap_t,kernel(1))     ! WK
    occupancy_buffer%structure = 'G'
    call sparse_create(occupancy_buffer)                         !WKV=G

    SPIN_loop: do is=1,pub_cell%num_spins

       ! Calculate Hubbard term G
       ! WK
       call sparse_product(site_buffer,hub_overlap_t,kernel(is))
       ! WKV blanked to 'G'
       call sparse_product(occupancy_buffer,site_buffer,hub_overlap)

       ! ddor: The Hamiltonian for the quadratic energy penalty functional
       if (pub_hub_functional == 1 ) then

          ! (1 - 2 WKV) blanked to 'G'
          call sparse_scale(occupancy_buffer,-2.0_DP,1.0_DP)

          ! ddor: The Hamiltonian for the quartic energy penalty functional
       elseif (pub_hub_functional == 2) then

          ! Create temporary matrices with 'G' structure
          call sparse_create(occupancy_buffer2,occupancy_buffer)    !G
          call sparse_create(occupancy_buffer3,occupancy_buffer)    !G

          ! (WKV)^2
          call sparse_product(occupancy_buffer2,occupancy_buffer, &
               occupancy_buffer2)
          ! -(WKV) (1 - WKV)
          call sparse_axpy(occupancy_buffer2,occupancy_buffer,-1.0_DP)
          ! (1 - 2 WKV)
          call sparse_copy(occupancy_buffer3,occupancy_buffer)
          call sparse_scale(occupancy_buffer3,-2.0_DP,1.0_DP)
          ! -(WKV) (1 - WKV) (1 - 2 WKV)
          call sparse_product(occupancy_buffer,occupancy_buffer2, &
               occupancy_buffer3)
          ! 2 (WKV) (1 - WKV) (1 - 2 WKV)
          call sparse_scale(occupancy_buffer,-2.0_DP)

          ! Destroy temporary matrices
          call sparse_destroy(occupancy_buffer3)
          call sparse_destroy(occupancy_buffer2)

       else
          write(stdout,*) 'Invalid value for pub_hub_functional, aborting.'
          call comms_abort
       end if

       ! U/2 H_DFT
       call sparse_product(hub%projector_ham(is),hub%u_matrix,occupancy_buffer)

       ! Add contribution due to hub_alpha or spin-splitting if necessary
       ! gibo:*******  WARNING:
       ! gibo:*******  for (pub_cdft_hubbard), add also the cDFT terms
       !               [in hub_up/down_matrix]
       ! gibo:*******  (see also sbrtne hubbard_build_spam &
       !               cdft_build_matrices)
       ! gibo:*******  WARNING:
       if ( ( ( SUM(ABS(h_species(:)%hub_spin_splitting)) + &
            &SUM(ABS(h_species(:)%hub_alpha)) ) .gt. 0.0_DP ) &
            & .OR. (pub_cdft_hubbard)) then
          if (is == 2) then
             call sparse_axpy(hub%projector_ham(is),hub%down_matrix,1.0_DP)
          else
             call sparse_axpy(hub%projector_ham(is),hub%up_matrix,1.0_DP)
          end if
       end if

    end do SPIN_loop

    call sparse_destroy(occupancy_buffer)
    call sparse_destroy(site_buffer)

    ! gibo: for CDFT-only simulations, add cDFT contribution to 'Hubbard' hamiltonian
    ! [in hub_up/down_matrix following call to hubbard_build_spam & cdft_build_matrices]
    elseif ((pub_cdft) .AND. (.not. pub_cdft_hubbard)) then

         do is=1,pub_cell%num_spins
           if (is == 2) then
             call sparse_copy(hub%projector_ham(is),hub%down_matrix)
           else
             call sparse_copy(hub%projector_ham(is),hub%up_matrix)
           end if
         end do

    end if HUBBARD_CDFT

    call timer_clock('hubbard_projector_ham',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_projector_ham'
#endif

  end subroutine hubbard_projector_ham


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_energy_total(hub,hubbard_e, &
       pur_denskern, hub_proj_basis, hub_overlap, hub_overlap_t)

    !==========================================================!
    ! Calculates the Hubbard DFT+U energy contribution to the  !
    ! total energy, and writes out the occupancy matrix for    !
    ! each correlated site.                                    !
    !----------------------------------------------------------!
    ! Written by David O'Regan in August 2009                  !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: pub_on_root,comms_abort, &
         pub_my_node_id, comms_barrier, comms_reduce
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: pub_num_hub_atoms_on_node, &
         pub_hub_atoms_on_node, pub_distr_atom
    use rundat, only: pub_hub_functional, pub_output_detail, &
         pub_cdft, pub_cdft_hubbard
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_axpy, sparse_copy, &
         sparse_product, sparse_create, sparse_destroy, &
         sparse_get_element, sparse_trace
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(out) :: hubbard_e  ! Hubbard DFT+U energy
    type(SPAM3), intent(in) :: pur_denskern(pub_cell%num_spins)
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(SPAM3), intent(in) :: hub_overlap, hub_overlap_t

    ! Local Variables
    type(SPAM3) :: site_buffer
    type(SPAM3) :: energy_buffer
    type(SPAM3) :: energy_buffer2
    integer :: channels
    integer :: hat_on_node, hat, sp, theatom
    integer :: is
    integer :: i,ihub
    real(kind=DP) :: occupancyelement, energyelement
    real(kind=DP) :: u,a,s

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering hubbard_energy_total'
#endif

    call timer_clock('hubbard_energy_total',1)

    ! gibo: Initialise energy and occupancies here not to mess with
    !       cdft_energy_total for CDFT
    hubbard_e = 0.0_DP
    hub%energy = 0.0_DP
    hub%population = 0.0_DP

    ! gibo: calculate Hubbard terms only for pure Hubbard-only
    !       .OR. cdft_hubbard simulations
    HUBBARD_CDFT_out: if ((.not.pub_cdft) .OR. (pub_cdft_hubbard) ) then

       ! ddor: Sanity check for spin-symmetry breaking
       if ( (sum(abs(h_species(:)%hub_spin_splitting)) .ne. 0.0_DP) .and. &
            & (pub_cell%num_spins .eq. 1) ) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in hubbard_energy_total: &
               &An attempt has been made to break spin symmetry &
               &while performing an unpolarised calculation.'
          call comms_abort
       end if

       call sparse_create(site_buffer,hub_overlap_t,pur_denskern(1))
       ! ddor: Blank occupancy matrix and products thereof to
       !       block-diagonal matrices in order to prevent
       !       Hubbard sites with null overlap interfering

       energy_buffer%structure = 'G'
       call sparse_create(energy_buffer)
       energy_buffer2%structure = 'G'
       call sparse_create(energy_buffer2)

       spin: do is=1,pub_cell%num_spins

          ! WK
          call sparse_product(site_buffer,hub_overlap_t,pur_denskern(is))
          ! WKV blanked to 'G'
          call sparse_product(hub%occupancy_matrix(is),site_buffer,hub_overlap)
          ! WKV.WKV blanked to 'G'
          call sparse_product(energy_buffer2,hub%occupancy_matrix(is), &
               hub%occupancy_matrix(is))

          call sparse_copy(energy_buffer,hub%occupancy_matrix(is))

          ! ddor: N - N^2
          call sparse_axpy(energy_buffer,energy_buffer2,-1.0_DP)

          ! ddor: Use quartic penalty functional instead of quadratic
          ! ddor: ( N - N^2 )^2
          if (pub_hub_functional == 2) then

             call sparse_copy(energy_buffer2,energy_buffer)
             call sparse_product(energy_buffer,energy_buffer2,energy_buffer2)

          end if

          ! ddor: Loop over Hubbard atoms on my node
          hubatoms: do hat_on_node=1,pub_num_hub_atoms_on_node(pub_my_node_id)

             hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
             theatom = pub_distr_atom(hub%orig(hat))
             channels = hub_proj_basis%num_on_atom(theatom)
             sp = hub%species_number(hat)

             if ( channels .ne. (2 * h_species(sp)%hub_ang_mom + 1) ) then
                write(stdout,'(a)') &
                     'Error in hubbard_energy_total: &
                     &hub_proj_basis%num_on_atom and hub_ang_mom mismatch'
                call comms_abort
             end if

             ! ddor: Local copies of Hubbard parameters
             u = h_species(sp)%hub_u
             a = h_species(sp)%hub_alpha
             s = h_species(sp)%hub_spin_splitting

             do i=1,channels
                ihub = hub_proj_basis%first_on_atom(theatom) + i - 1
                call sparse_get_element(energyelement,energy_buffer, &
                     ihub,ihub)
                hub%energy(is,hat) = hub%energy(is,hat) &
                     + energyelement
                call sparse_get_element(occupancyelement, &
                     hub%occupancy_matrix(is),ihub,ihub)
                hub%population(is,hat) = hub%population(is,hat) &
                     + occupancyelement
             end do

             ! ddor: For spin-polarised systems
             !       Calculate energy
             !       U/2 Tr[N - N^2] + alpha Tr[N] +
             !       spin_spitting/2 Tr[N_down - N_up]
             !       That is (U/2 + alpha + spin_splitting/2) Tr[N_down] +
             !            (U/2 + alpha - spin_splitting/2) Tr[N_up] - U/2 Tr[N^2]
             if (pub_cell%num_spins == 2) then
                if (is == 1) then
                   hub%energy(is,hat) = &
                        ( ((-0.5_DP*s)+a)* hub%population(is,hat) ) + &
                        ( (u*0.5_DP)* hub%energy(is,hat) )
                elseif (is == 2) then
                   hub%energy(is,hat) = &
                        ( ((0.5_DP*s)+a)* hub%population(is,hat) ) + &
                        ( (u*0.5_DP)* hub%energy(is,hat) )
                end if
                ! ddor: For un-spin-polarised systems
                !       Calculate energy U Tr[N - N^2] + 2*alpha Tr[N]
                !       That is (U + 2*alpha) Tr[N] - U Tr[N^2]
                !       Since we are dealing with a one-electron density matrix
             elseif  (pub_cell%num_spins == 1) then
                hub%energy(is,hat) = &
                     ( (2.0_DP*a)* hub%population(is,hat) ) + &
                     (u* hub%energy(is,hat) )
             end if

             hubbard_e = hubbard_e + hub%energy(is,hat)

          end do hubatoms

       end do spin

       ! gibo: calculate DFT+U occupancy matrices and energies for pure
       ! gibo: Hubbard-only simulation
       HUBBARD_CDFT_in: if (.not.pub_cdft) then

          call comms_reduce('SUM',hub%energy)
          call comms_reduce('SUM',hub%population)

          ! ddor: Write out the DFT+U occupancy matrices and energies
          call comms_reduce('SUM', hubbard_e )
          if (pub_output_detail == VERBOSE) &
               call hubbard_energy_info(hub,hub_proj_basis)

       end if HUBBARD_CDFT_in

       call sparse_destroy(energy_buffer2)
       call sparse_destroy(energy_buffer)
       call sparse_destroy(site_buffer)

    end if HUBBARD_CDFT_out

    ! gibo: for cDFT runs, calculate cDFT energy energy and add it to hubbard_e
    CDFT_ENERGY: if (pub_cdft) then

       call cdft_energy_total(hub,hubbard_e, &
            pur_denskern, hub_proj_basis, hub_overlap, hub_overlap_t)

    end if CDFT_ENERGY

    call timer_clock('hubbard_energy_total',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_energy_total'
#endif

  end subroutine hubbard_energy_total


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_energy_info(hub,hub_proj_basis)

    !==========================================================!
    ! Writes out occupancies and energies of Hubbard sites.    !
    !----------------------------------------------------------!
    ! Written by David O'Regan in November 2009                !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: comms_abort, comms_bcast, pub_on_root, &
         pub_my_node_id
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: pub_node_of_atom, pub_distr_atom
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_get_block
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(in) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    integer :: species, this_species_atom ! ddor: Loop over species
    integer :: jj,hat,project,spin,theatom, node
    integer :: ierr
    real(kind=DP) :: local_population, local_moment, local_energy
    real(kind=DP),allocatable :: occupancy_block(:,:,:)
    character(10) :: fmt, tmp

    allocate(occupancy_block(hub_proj_basis%max_on_atom, &
         hub_proj_basis%max_on_atom,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('hubbard_energy_info','occupancy_block',ierr)

    if(pub_my_node_id==0) write(*,*) 'hubbard_energy_info'

    do species=1,pub_cell%num_hub_species

       this_species_atom = 0

       atomloop: do hat=1,pub_cell%nat_hub

          if (hub%species_number(hat) == species) then

             this_species_atom = this_species_atom + 1
             if (this_species_atom .gt. 2) exit atomloop

             theatom = pub_distr_atom(hub%orig(hat))
             project = hub_proj_basis%num_on_atom(theatom)
             node    = pub_node_of_atom(hub%orig(hat))

             if (node==pub_my_node_id) then
                do spin=1,pub_cell%num_spins
                   call sparse_get_block(occupancy_block(:,:,spin), &
                        hub%occupancy_matrix(spin),theatom,theatom)
                end do
             end if
             call comms_bcast(node,occupancy_block)

             if (.not.pub_on_root) cycle

             write(stdout,'(/a)')'########################################&
                  &########################################'

!CW
             if(node==pub_my_node_id) then
               call diago_occupancy(h_atoms_occupancy(hub_proj_basis,hub,hat,1),project,hat)
             endif
!END CW

             ! ddor: TELL ME ABOUT THIS ATOM
             write(stdout,'(a,i6,a,a)') 'DFT+U information on atom ', &
                  this_species_atom,' of Hubbard species ',&
                  h_species(species)%hub_species

             do spin=1,pub_cell%num_spins

                write(stdout,'(a)')'########################################&
                     &########################################'

                if (pub_cell%num_spins == 2) then
                   ! ddor: WRITE OUT INFORMATION ON OCCUPANCY MATRIX
                   write(stdout,'(a,2(i6,a))') 'Occupancy matrix of &
                        &Hubbard site ',hat,' and spin ',spin,' is '
                else

                end if

                write(tmp,'(i6)') project
                write(fmt,'(3a)') '(',trim(adjustl(tmp)),'f12.8)'
                do jj=1,project
                   write(stdout,fmt) occupancy_block(jj,1:project,spin)
                end do

                if (MAXVAL(occupancy_block(:,:,spin)) .gt. 1.0_DP) then
                   write(stdout,'(a,2(i6,a))') 'WARNING: OCCUPANCY MATRIX &
                        &of Hubbard site ',&
                        hat,' and spin ',spin,' exceeeds 1.'
                end if
                if (hub%energy(spin,hat) .lt. 0.0_DP) then
                   write(stdout,'(a,2(i6,a))') 'WARNING: DFT+U ENERGY &
                        &of Hubbard site ',&
                        hat,' and spin ',spin,' is negative.'
                end if

             end do

             write(stdout,'(a)')'########################################&
                  &########################################'

             ! ddor: WRITE OUT INFORMATION ON OCCUPANCY OF SITE
             if (pub_cell%num_spins==2) then
                local_population = sum(hub%population(:,hat))
                local_energy = hub%energy(1,hat) + &
                     hub%energy(2,hat)
             else
                local_population = 2.0_DP*hub%population(1,hat)
                local_energy = hub%energy(1,hat)
             end if

             write(stdout,'(a,i6,a,f12.8,a)') 'Total occupancy &
                  &of Hubbard site ', hat,' is       ', &
                  &local_population, ' e'

             if (pub_cell%num_spins==2) then
                local_moment =  &
                     &hub%population(1,hat) - hub%population(2,hat)
                write(stdout,'(a,i6,a,f12.8,a)') 'Local magnetic moment &
                     &of Hubbard site ', hat,' is ', local_moment,' mu_B'
             end if

             ! ddor: WRITE OUT INFORMATION ON DFT+U ENERGY
             write(stdout,'(a,i6,a,f12.8,a)') 'DFT+U energy of Hubbard &
                  &site ',hat,' is          ', local_energy,' Ha'

             write(stdout,'(a/)')'########################################&
                  &########################################'

          end if

       end do atomloop

    end do

    deallocate(occupancy_block,stat=ierr)
    call utils_dealloc_check('hubbard_energy_info','occupancy_block',ierr)

  end subroutine hubbard_energy_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cdft_energy_total(hub,hubbard_e, &
       pur_denskern, hub_proj_basis, hub_overlap, hub_overlap_t)

    !==========================================================!
    ! Calculates the cDFT energy contribution to the Hubbard   !
    ! (thence total) energy,                                   !
    ! and writes out the occupancy matrix for                  !
    ! each correlated site.                                    !
    !----------------------------------------------------------!
    ! Adapted from hubbard_energy_total by gibo in Nov. 2011   !
    !==========================================================!

    use comms, only: pub_on_root,comms_abort, &
         pub_my_node_id, comms_barrier, comms_reduce
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: pub_num_hub_atoms_on_node, &
         pub_hub_atoms_on_node, pub_distr_atom
    use rundat, only: pub_hub_functional, pub_output_detail, &
         pub_cdft, pub_cdft_hubbard, &
         pub_cdft_atom_charge, pub_cdft_atom_spin, &
         pub_cdft_group_charge, pub_cdft_group_spin, &
         pub_cdft_group_charge_diff, pub_cdft_group_spin_diff, &
         pub_cdft_group_charge_target, pub_cdft_group_spin_target, &
         pub_cdft_group_charge_diff_target, pub_cdft_group_spin_diff_target, &
         pub_cdft_group_u, pub_cdft_group_diff_u, &
         pub_cdft_type
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_axpy, sparse_copy, &
         sparse_product, sparse_create, sparse_destroy, &
         sparse_get_element, sparse_trace
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(inout) :: hubbard_e  ! Hubbard DFT+U energy
    type(SPAM3), intent(in) :: pur_denskern(pub_cell%num_spins)
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(SPAM3), intent(in) :: hub_overlap, hub_overlap_t

    ! Local Variables
    type(SPAM3) :: site_buffer

    integer :: channels
    integer :: hat_on_node, hat, sp, theatom
    integer :: is
    integer :: i,ihub
    real(kind=DP) :: occupancyelement
    real(kind=DP) :: cdft_e  ! cDFT energy
    real(kind=DP) :: u_charge_up, u_charge_down
    real(kind=DP) :: u_spin
    real(kind=DP) :: n_t_up, n_t_down, delta_n_t
    real(kind=DP) :: cdft_penalty
    real(kind=DP) :: pop_acceptor, pop_donor
    real(kind=DP) :: spinmom_acceptor, spinmom_donor
    real(kind=DP) :: pop_group, spinmom_group

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &cdft_energy_total'
#endif

    call timer_clock('cdft_energy_total',1)

    call sparse_create(site_buffer,hub_overlap_t,pur_denskern(1))

    ! initialise cDFT-energy to zero
    cdft_e = 0.0_DP

    ! initialise cdft_gradient to zero
    hub%cdft_gradient = 0.0_DP


    ! For cDFT+U runs [pub_cdft_hubbard=.T.],
    ! we have already calculated the populations in hubbard_energy_total
    ! so, no need to recalculate them here...

    CDFT_HUBBARD: if (.not.pub_cdft_hubbard) then

    ! Calculate population first,
    ! then work out energy depending on given cdft-mode.
    ! Slightly less efficient (we repeat the spin_ and hubatoms_ loops) for
    ! cdft_atom_charge = T, yet save a lot of repeated coding

    spin: do is = 1, pub_cell%num_spins

       ! WK
       call sparse_product(site_buffer,hub_overlap_t,pur_denskern(is))
       ! WKV blanked to 'G'
       call sparse_product(hub%occupancy_matrix(is),site_buffer,hub_overlap)

       ! ddor: Loop over Hubbard atoms on my node
       hubatoms : do hat_on_node=1,pub_num_hub_atoms_on_node(pub_my_node_id)

          hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
          theatom = pub_distr_atom(hub%orig(hat))
          channels = hub_proj_basis%num_on_atom(theatom)
          sp = hub%species_number(hat)

          if ( channels .ne. (2 * h_species(sp)%hub_ang_mom + 1) ) then
             write(stdout,'(a)') &
                  'Error in cdft_energy_total: &
                  &hub_atom_hub_projs and hub_ang_mom mismatch'
             call comms_abort
          end if

          do i=1,channels
             ihub = hub_proj_basis%first_on_atom(theatom) + i - 1
             call sparse_get_element(occupancyelement,hub%occupancy_matrix(is), &
                  ihub,ihub)
             hub%population(is,hat) = hub%population(is,hat) &
                  + occupancyelement
          end do

       end do hubatoms

    end do spin
 
    endif CDFT_HUBBARD

    ! Once the cDFT populations is calculated,
    ! deal with them according to the selected cDFT-mode

!====select CDFT-mode and iterate...
    CDFT_MODE: SELECT CASE (pub_cdft_type)

    ! gibo: ATOM-CHARGE constrained run: pub_cdft_type=1
    case(1)

    ! gibo: ATOM-CHARGE constrained run
    !CDFT_MODE: if (pub_cdft_atom_charge) then
!====select CDFT-mode and iterate...

       ! ddor: Loop over Hubbard atoms on my node
       hubatoms_1 : do hat_on_node=1,pub_num_hub_atoms_on_node(pub_my_node_id)

          hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
          sp = hub%species_number(hat)

          ! gibo: local copies of CONSTRAINED_DFT CHARGE parameters
          u_charge_up    = hub%cdft_u_charge_up(hat)
          u_charge_down  = hub%cdft_u_charge_down(hat)

          n_t_up      = h_species(sp)%cdft_target_up
          n_t_down    = h_species(sp)%cdft_target_down

       spin_1: do is = 1, pub_cell%num_spins

          if (is == 1) then

           ! U_q (Tr[n_up] - n_t_up)
           cdft_penalty = hub%population(is,hat) - n_t_up
           hub%energy(is,hat) =  hub%energy(is,hat) + u_charge_up * cdft_penalty

           cdft_e = cdft_e + hub%energy(is,hat)

          elseif (is == 2) then

           ! U_q (Tr[n_down] - n_t_down)
           cdft_penalty = hub%population(is,hat) - n_t_down
           hub%energy(is,hat) =  hub%energy(is,hat) + u_charge_down * cdft_penalty

           cdft_e = cdft_e + hub%energy(is,hat)

          end if

           ! store the gradient
           hub%cdft_gradient(is,hat) = cdft_penalty

        end do spin_1

       end do hubatoms_1


!====select CDFT-mode and iterate...
    ! gibo: ATOM-SPIN constrained run: pub_cdft_type=2
    case(2)

!    ! gibo: ATOM-SPIN constrained run
!    elseif (pub_cdft_atom_spin) then
!====select CDFT-mode and iterate...


    ! gibo: add spin penalization  U_s * (Tr[Nup] - Tr[Ndown] - Delta_n_t)
    hubatoms_cdft_2 : do hat_on_node=1,pub_num_hub_atoms_on_node(pub_my_node_id)

        hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)

        sp = hub%species_number(hat)

        ! gibo: local copies of CONSTRAINED_DFT SPIN parameters
        u_spin      = hub%cdft_u_spin(hat)
        delta_n_t   = h_species(sp)%cdft_target_spin

        ! hub%population(is,hat) has been calculated in SPIN_2: do-end do loop above
        cdft_penalty = hub%population(1,hat) - hub%population(2,hat) - delta_n_t

        ! distribute equally the spin-penalization between the UP- and DOWN-spins...
        do is = 1, pub_cell%num_spins

          hub%energy(is,hat) =  hub%energy(is,hat) + 0.5_DP * u_spin * cdft_penalty

          cdft_e = cdft_e + hub%energy(is,hat)

           ! store the gradient 
           ! [equally distributed between the UP- and DOWN- spin channels]
           hub%cdft_gradient(is,hat) = 0.5_DP * cdft_penalty

        end do

    end do hubatoms_cdft_2

!====select CDFT-mode and iterate...
    ! gibo: GROUP-CHARGE constrained run: pub_cdft_type=3
    case(3)

!    ! gibo: GROUP-CHARGE constrained run
!    elseif (pub_cdft_group_charge) then
!====select CDFT-mode and iterate...

    ! gibo: initialise total population of acceptor and donor (to zero)
    pop_group = 0._DP

       ! ddor: Loop over Hubbard atoms on my node
       hubatoms_3 : do hat_on_node=1,pub_num_hub_atoms_on_node(pub_my_node_id)

          hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)

        ! gibo: loop over spins
        spin_3: do is = 1, pub_cell%num_spins

          ! for group-charge-constrained cDFT runs, keep track of total populations
          pop_group = pop_group + hub%population(is,hat)

         end do spin_3

        end do hubatoms_3


    ! collect and sum group population
    call comms_reduce('SUM', pop_group)

    ! calculate group-constrained energy
    cdft_penalty = pop_group - pub_cdft_group_charge_target
    cdft_e = pub_cdft_group_u * cdft_penalty

    ! collect and sum hubbard_e
    call comms_reduce('SUM', hubbard_e )

    ! add cDFT energy to Hubbard energy
    hubbard_e = hubbard_e + cdft_e

    ! split equally e_cdft among CDFT_atoms...
    cdft_e = cdft_e/real(pub_cell%nat_hub,kind=DP)

    ! split equally gradient (cdft_penalty) among CDFT_atoms...
    cdft_penalty = cdft_penalty/real(pub_cell%nat_hub,kind=DP)

    do hat_on_node=1,pub_num_hub_atoms_on_node(pub_my_node_id)

       hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)

       !...and spin channels (0.5 factor)
       do is = 1, pub_cell%num_spins
          hub%energy(is,hat) =  hub%energy(is,hat) + 0.5_DP*cdft_e
          hub%cdft_gradient(is,hat) =  0.5_DP*cdft_penalty
       end do

    end do

!====select CDFT-mode and iterate...
    ! gibo: GROUP-SPIN constrained run: pub_cdft_type=4
    case(4)

!    ! gibo: GROUP-SPIN constrained run
!    elseif (pub_cdft_group_spin) then
!====select CDFT-mode and iterate...

    ! gibo: initialise total population of acceptor and donor (to zero)
    spinmom_group = 0._DP

       ! ddor: Loop over Hubbard atoms on my node
       hubatoms_4 : do hat_on_node =1,pub_num_hub_atoms_on_node(pub_my_node_id)

          hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)

          ! for group-spin-constrained cDFT runs, keep track of total spin_moment
          spinmom_group = spinmom_group + &
               hub%population(1,hat) - hub%population(2,hat)

        end do hubatoms_4


    ! collect and sum group spin-moment
    call comms_reduce('SUM', spinmom_group)

    ! calculate group-constrained energy
    cdft_penalty = spinmom_group - pub_cdft_group_spin_target
    cdft_e = pub_cdft_group_u * cdft_penalty

    ! collect and sum hubbard_e
    call comms_reduce('SUM', hubbard_e )

    ! add cDFT energy to Hubbard energy
    hubbard_e = hubbard_e + cdft_e

    ! split equally e_cdft among CDFT_atoms...
    cdft_e = cdft_e/real(pub_cell%nat_hub,kind=DP)

    ! split equally gradient (cdft_penalty) among CDFT_atoms...
    cdft_penalty = cdft_penalty/real(pub_cell%nat_hub,kind=DP)

    do hat_on_node=1,pub_num_hub_atoms_on_node(pub_my_node_id)

      hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)

      !...and spin channels (0.5 factor)
      do is = 1, pub_cell%num_spins
       hub%energy(is,hat) =  hub%energy(is,hat) + 0.5_DP*cdft_e
       hub%cdft_gradient(is,hat) =  0.5_DP*cdft_penalty
      end do

    end do

!====select CDFT-mode and iterate...
!    ! gibo: GROUP-CHARGE_DIFFERENCE constrained run:  pub_cdft_type=5
     case(5)

!    ! gibo: GROUP-CHARGE_DIFFERENCE constrained run
!    elseif (pub_cdft_group_charge_diff) then
!====select CDFT-mode and iterate...

    ! gibo: initialise total population of acceptor and donor (to zero)
    pop_acceptor = 0._DP
    pop_donor    = 0._DP

       ! ddor: Loop over Hubbard atoms on my node
       hubatoms_5 : do hat_on_node = 1, pub_num_hub_atoms_on_node(pub_my_node_id)

          hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)

          sp = hub%species_number(hat)

      ! gibo: loop over spin
      spin_5: do is = 1, pub_cell%num_spins

          ! for group-charge-constrained cDFT runs, keep track of total
          ! populations
          if (h_species(sp)%cdft_acceptor) &
              pop_acceptor = pop_acceptor + hub%population(is,hat)

          if (h_species(sp)%cdft_donor)    &
              pop_donor = pop_donor + hub%population(is,hat)

        end do spin_5

        end do hubatoms_5

    ! collect and sum acceptor/donor populations
    call comms_reduce('SUM', pop_acceptor)
    call comms_reduce('SUM', pop_donor)

    ! calculate group-constrained energy
    ! gibo: penalty inverted to be consistent with -ve (+ve) potentials at
    !       acceptor (donor) sites [pub_cdft_group_diff_u has been set +ve in 
    !       rundat_blocks_mod.F90]
    cdft_penalty = pub_cdft_group_charge_diff_target - pop_acceptor + pop_donor
    !cdft_penalty = pop_acceptor - pop_donor - pub_cdft_group_charge_diff_target
    cdft_e = pub_cdft_group_diff_u * cdft_penalty

    ! collect and sum hubbard_e
    call comms_reduce('SUM', hubbard_e )

    ! add cDFT energy to Hubbard energy
    hubbard_e = hubbard_e + cdft_e

    ! split equally e_cdft among CDFT_atoms...
    cdft_e = cdft_e/real(pub_cell%nat_hub,kind=DP)

    ! split equally gradient (cdft_penalty) among CDFT_atoms...
    cdft_penalty = cdft_penalty/real(pub_cell%nat_hub,kind=DP)

    do hat_on_node = 1, pub_num_hub_atoms_on_node(pub_my_node_id)

      hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
      sp = hub%species_number(hat)

      !...and spin channels (0.5 factor)
      do is = 1, pub_cell%num_spins
       hub%energy(is,hat) =  hub%energy(is,hat) + 0.5_DP*cdft_e
       hub%cdft_gradient(is,hat) =  0.5_DP*cdft_penalty
       ! gibo: for acceptors, invert sign of gradient to be consistent
       !      with (+ve) sign of pub_cdft_group_diff_u
       if (h_species(sp)%cdft_acceptor) &
          hub%cdft_gradient(is,hat) = -hub%cdft_gradient(is,hat)
      end do

    end do


!====select CDFT-mode and iterate...
    ! gibo: GROUP-SPIN_DIFFERENCE constrained run: pub_cdft_type=6
    case(6)

!    ! gibo: GROUP-SPIN_DIFFERENCE constrained run
!    elseif (pub_cdft_group_spin_diff) then
!====select CDFT-mode and iterate...

    ! gibo: initialise total population of acceptor and donor (to zero)
    spinmom_acceptor = 0._DP
    spinmom_donor    = 0._DP

    ! work out spin-pipulation difference (local spin_moment) for acceptor and donor
    hubatoms_cdft_6 : do hat_on_node = 1, pub_num_hub_atoms_on_node(pub_my_node_id)

        hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)

        sp = hub%species_number(hat)

         if (h_species(sp)%cdft_acceptor) &
             spinmom_acceptor = spinmom_acceptor + &
                                hub%population(1,hat) - hub%population(2,hat)

         if (h_species(sp)%cdft_donor)    &
              spinmom_donor = spinmom_donor + &
                              hub%population(1,hat) - hub%population(2,hat)

    end do hubatoms_cdft_6

    ! collect and sum acceptor/donor spin_moments
    call comms_reduce('SUM', spinmom_acceptor)
    call comms_reduce('SUM', spinmom_donor)

    ! calculate group-constrained energy
    ! gibo: penalty inverted to be consistent with -ve (+ve) potentials at
    !       acceptor (donor) sites [pub_cdft_group_diff_u has been set +ve in 
    !       rundat_blocks_mod.F90]
    cdft_penalty = pub_cdft_group_spin_diff_target - spinmom_acceptor + spinmom_donor
    !cdft_penalty = spinmom_acceptor - spinmom_donor - pub_cdft_group_spin_diff_target
    cdft_e = pub_cdft_group_diff_u * cdft_penalty

    ! collect and sum hubbard_e
    call comms_reduce('SUM', hubbard_e )

    ! add cDFT energy to Hubbard energy
    hubbard_e = hubbard_e + cdft_e

    ! split equally e_cdft among CDFT_atoms...
    cdft_e = cdft_e/real(pub_cell%nat_hub,kind=DP)

    ! split equally gradient (cdft_penalty) among CDFT_atoms...
    cdft_penalty = cdft_penalty/real(pub_cell%nat_hub,kind=DP)

    do hat_on_node = 1, pub_num_hub_atoms_on_node(pub_my_node_id)

       hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
       sp = hub%species_number(hat)

       !...and spin channels (0.5 factor)
       do is = 1, pub_cell%num_spins
       hub%energy(is,hat) =  hub%energy(is,hat) + 0.5_DP*cdft_e
       hub%cdft_gradient(is,hat) =  0.5_DP*cdft_penalty
       ! gibo: for acceptors, invert sign of gradient to be consistent
       !      with (+ve) sign of pub_cdft_group_diff_u
       !MIND: since magn-mom acceptors are DOWN-repulsive
       !      gradient for is=2 needs to be inverted again...
       !      the code below take care of that

       ! invert only UP-gradient for acceptor-atoms
       if ((h_species(sp)%cdft_acceptor).AND.(is==1)) &
          hub%cdft_gradient(is,hat) = -hub%cdft_gradient(is,hat)

       ! invert only DOWN-gradient for donor-atoms
       if ((h_species(sp)%cdft_donor).AND.(is==2)) &
          hub%cdft_gradient(is,hat) = -hub%cdft_gradient(is,hat)
      end do

    end do

!====select CDFT-mode and iterate...
    ! gibo: we somehow lost info on cDFT-mode, shut everything down...
    case default
      write(stdout,'(a)') &
           'Error in cdft_energy_total: unassigned pub_cdft_type run'
      write(stdout,'(a)') &
           'check cdft_atom/group_charge/spin[_diff] keywords in input file'
      call comms_abort
!====select CDFT-mode and iterate...

!====select CDFT-mode and iterate...
    end select CDFT_MODE
    !end if CDFT_MODE
!====select CDFT-mode and iterate...


    ! gibo: for ATOM-CHARGE/SPIN-constrained simulation,
    ! update hubbard_e and write out the cDFT(+U) occupancy matrices and energies
    ATOMS: if (pub_cdft_atom_charge .OR. pub_cdft_atom_spin) then

      ! add cDFT energy to Hubbard energy
      hubbard_e = hubbard_e + cdft_e


      ! gibo: MIND that at this point hubbard_e contains also cdft_e (see above)
      call comms_reduce('SUM', hubbard_e )

    end if ATOMS

    ! Write out the cDFT(+U) occupancy matrices and energies=== START
    call comms_reduce('SUM',hub%energy)!
    call comms_reduce('SUM',hub%population)

    if (pub_output_detail == VERBOSE) &
         call cdft_energy_info(hub,hub_proj_basis)
    ! Write out the cDFT(+U) occupancy matrices and energies=== END


    ! gibo: for cdft_group_charge/spin[_difference] print out
    ! info on total population/spin-moment [difference]
    TOT_POP: if ((.not.pub_cdft_atom_charge) .OR. (.not.pub_cdft_atom_spin)) then
    ROOT: if (pub_on_root) then

       write(stdout,'(/a)')'########################################&
                          &########################################'

       if (pub_cdft_group_charge) then
         write(stdout,'(a)') 'cDFT information on GROUP-CHARGE constrained run'
         write(stdout,'(a)')'########################################&
                            &########################################'
         write(stdout,'(2(a,f12.8))') 'Total group population:               ',&
              &pop_group,&
              &' e    / Target value: ', pub_cdft_group_charge_target

       elseif (pub_cdft_group_spin) then
         write(stdout,'(a)') 'cDFT information on GROUP-SPIN constrained run'
         write(stdout,'(a)')'########################################&
                            &########################################'
         write(stdout,'(2(a,f12.8))') 'Total group magnetic-moment:          ',&
              &spinmom_group,&
              &' mu_B / Target value: ', pub_cdft_group_spin_target

       elseif (pub_cdft_group_charge_diff) then
         write(stdout,'(a)') 'cDFT information on GROUP-CHARGE-DIFFERENCE constrained run'
         write(stdout,'(a)')'########################################&
                            &########################################'
         write(stdout,'(a,f12.8,a)')    'Total ACCEPTOR population:            ',&
              &pop_acceptor, ' e'
         write(stdout,'(a,f12.8,a)')    'Total DONOR population:               ',&
              &pop_donor,    ' e'
         write(stdout,'(2(a,f12.8))') 'ACCEPTOR-DONOR population diff.:      ',&
              &pop_acceptor - pop_donor,&
              &' e    / Target value: ', pub_cdft_group_charge_diff_target

       elseif (pub_cdft_group_spin_diff) then
         write(stdout,'(a)') 'cDFT information on GROUP-SPIN-DIFFERENCE constrained run'
         write(stdout,'(a)')'########################################&
                            &########################################'
         write(stdout,'(a,f12.8,a)')    'Total ACCEPTOR magnetic-moment:       ',&
              &spinmom_acceptor,' mu_B'
         write(stdout,'(a,f12.8,a)')    'Total DONOR magnetic-moment:          ',&
              &spinmom_donor,   ' mu_B'
         write(stdout,'(2(a,f12.8))') 'ACCEPTOR-DONOR magnetic-moment diff.: ',&
              &spinmom_acceptor - spinmom_donor,&
              &' mu_B / Target value: ', pub_cdft_group_spin_diff_target

       end if

       write(stdout,'(a/)')'########################################&
                          &########################################'

    end if ROOT
    end if TOT_POP

    ! gibo: destroy local SPAM3 matrixes
    call sparse_destroy(site_buffer)

    call timer_clock('cdft_energy_total',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &cdft_energy_total'
#endif

  end subroutine cdft_energy_total


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cdft_energy_info(hub,hub_proj_basis)

    !==========================================================!
    ! Writes out occupancies and energies of Hubbard sites.    !
    !----------------------------------------------------------!
    ! Written by David O'Regan in November 2009                !
    ! Adapted to cDFT by Gilberto Teobaldi in November 2011    !
    !==========================================================!

    use comms, only: comms_abort, comms_bcast, pub_on_root, &
         pub_my_node_id
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use parallel_strategy, only: pub_node_of_atom, pub_distr_atom
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_get_block
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(in) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    integer :: species, this_species_atom
    integer :: jj,hat,project,spin,theatom, node
    integer :: ierr
    real(kind=DP) :: local_population, local_moment, local_energy
    real(kind=DP) :: local_population_up, local_population_down
    real(kind=DP),allocatable :: occupancy_block(:,:,:)
    character(10) :: fmt, tmp

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &cdft_energy_info'
#endif

    allocate(occupancy_block(hub_proj_basis%max_on_atom, &
         hub_proj_basis%max_on_atom,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('hubbard_energy_info','occupancy_block',ierr)

    write(*,*) 'cdft_energy_info'

    do species=1,pub_cell%num_hub_species

       this_species_atom = 0

       atomloop: do hat=1,pub_cell%nat_hub

          if (hub%species_number(hat) == species) then

             this_species_atom = this_species_atom + 1

             !if (this_species_atom .gt. 2) exit atomloop
! 4 GIBO: add logical to print all the populations of the cDFT sites one the calculation is converged...
             !if ( (this_species_atom .gt. 2) .AND. (.not. cdft_done) ) exit atomloop
! 4 GIBO: add logical to print all the populations of the cDFT sites one the calculation is converged...

             theatom = pub_distr_atom(hub%orig(hat))
             project = hub_proj_basis%num_on_atom(theatom)
             node = pub_node_of_atom(hub%orig(hat))

             if (node==pub_my_node_id) then
                do spin=1,pub_cell%num_spins
                   call sparse_get_block(occupancy_block(:,:,spin), &
                        hub%occupancy_matrix(spin),theatom,theatom)
                end do
             end if
             call comms_bcast(node,occupancy_block)

             if (.not.pub_on_root) cycle

             write(stdout,'(/a)')'########################################&
                  &########################################'


             ! ddor: TELL ME ABOUT THIS ATOM
             write(stdout,'(a,i6,a,a)') 'cDFT information on &
                  &atom ',this_species_atom,' of constrained species ',&
                  h_species(species)%hub_species

             if (pub_cell%num_spins == 2) then

                do spin=1, pub_cell%num_spins

                   write(stdout,'(a)')'########################################&
                        &########################################'

                   ! ddor: WRITE OUT INFORMATION ON OCCUPANCY MATRIX
                   write(stdout,'(a,2(i6,a))') 'Occupancy matrix of &
                        &cDFT site ', hat,' and spin ',spin,' is '

                   write(tmp,'(i6)') project
                   write(fmt,'(3a)') '(',trim(adjustl(tmp)),'f12.8)'
                   do jj=1,project
                      write(stdout,fmt) occupancy_block(jj,1:project,spin)
                   end do

                   if (MAXVAL(occupancy_block(:,:,spin)) .gt. 1.0_DP) then
                      write(stdout,'(a,2(i6,a))') 'WARNING: OCCUPANCY MATRIX &
                           &of cDFT site ',&
                           hat,' and spin ',spin,' exceeeds 1.'
                   end if

                end do

                write(stdout,'(a)')'########################################&
                     &########################################'

                ! ddor: WRITE OUT INFORMATION ON OCCUPANCY OF SITE
                local_population_up   =  hub%population(1,hat)
                local_population_down =  hub%population(2,hat)

                local_population = hub%population(1,hat) &
                     + hub%population(2,hat)
                local_moment = hub%population(1,hat) &
                     - hub%population(2,hat)
                local_energy = hub%energy(1,hat) &
                     + hub%energy(2,hat)
                write(stdout,'(a,i6,a,f12.8,a)') 'Total occupancy &
                     &of cDFT site       ',  hat,' is ', &
                     &local_population, ' e'
                write(stdout,'(a,i6,a,f12.8,a)') 'UP occupancy &
                     &of cDFT site          ',  hat,' is ', &
                     &local_population_up, ' e'
                write(stdout,'(a,i6,a,f12.8,a)') 'DOWN occupancy &
                     &of cDFT site        ',  hat,' is ', &
                     &local_population_down, ' e'
!                ! gibo: print out local gradient component
!                !      (left here for debugging)
!                write(stdout,'(a,i6,a,f12.8)')   'UP gradient    &
!                     &of cDFT site        ',  hat,' is ', &
!                     &cdft_gradient(atom,1)
!                write(stdout,'(a,i6,a,f12.8)')   'DOWN gradient  &
!                     &of cDFT site        ',  hat,' is ', &
!                     &cdft_gradient(atom,2)
                write(stdout,'(a,i6,a,f12.8,a)') 'Local magnetic moment &
                     &of cDFT site ',  hat,' is ', local_moment,' mu_B'

                ! ddor: WRITE OUT INFORMATION ON DFT+U ENERGY
                write(stdout,'(a,i6,a,f12.8,a)') 'cDFT energy of cDFT &
                     &site           ', hat,' is ', local_energy,' Ha'

                write(stdout,'(a/)')'########################################&
                     &########################################'

             end if

          end if

       end do atomloop

    end do

    deallocate(occupancy_block,stat=ierr)
    call utils_dealloc_check('hubbard_energy_info','occupancy_block',ierr)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &cdft_energy_info'
#endif

  end subroutine cdft_energy_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_ngwf_select(hub,ngwf_basis,hub_proj_basis, &
       hub_consist_matrix,hub_consist_tmatrix)

    !==========================================================!
    ! This subroutine determines which NGWFs have the greatest !
    ! projection onto a set of localised hydrogenic orbitals   !
    ! which are chosen by the user.                            !
    !----------------------------------------------------------!
    ! Written by David O'Regan in April 2009.                  !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: pub_my_node_id, comms_abort, comms_barrier, &
         comms_reduce, comms_bcast, pub_on_root, pub_root_node_id
    use function_basis, only: FUNC_BASIS
    use hubbard_init, only: h_species
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_num_hub_atoms_on_node, &
         pub_hub_atoms_on_node, pub_distr_atom
!CW
    use rundat, only: pub_output_detail,pub_dmft_order_proj,pub_dmft_switch_off_proj_order
!END CW
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_get_element
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use wrappers, only: wrappers_invert_sym_matrix

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(SPAM3), intent(in) :: hub_consist_matrix
    type(SPAM3), intent(in) :: hub_consist_tmatrix

    ! Local Variables
    integer, allocatable, dimension(:) :: ll
    integer :: ierr
!CW
    integer :: ijk,mm(1),jj
    logical :: ngwf_order
    integer :: jjj
!END CW
    integer :: hat_on_node, hat, theatom
    integer :: olaps, projs
    integer :: ngwf_we_are_on, hub_proj_we_are_on
    integer :: ngwf_count, hub_proj_count
    integer :: kk
    real(kind=DP) :: overlap_element, transpose_overlap_element
    real(kind=DP), allocatable :: ngwfsquareovlp(:,:,:)
    real(kind=DP), allocatable :: ngwfsquareoverlapsum(:)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_ngwf_select'
#endif

    call timer_clock('hubbard_ngwf_select',1)

    ! Allocate workspace
    allocate(ll(1),stat=ierr)
    call utils_alloc_check('hubbard_ngwf_select','ll',ierr)
    allocate(ngwfsquareovlp(ngwf_basis%max_on_atom,hub_proj_basis%max_on_atom, &
         pub_cell%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_ngwf_select','ngwfsquareovlp',ierr)
    allocate(ngwfsquareoverlapsum(ngwf_basis%max_on_atom),stat=ierr)
    call utils_alloc_check('hubbard_ngwf_select','ngwfsquareoverlapsum',ierr)

    ngwfsquareovlp(:,:,:) = 0.0_DP

    ! ddor: Loop over Hubbard atoms on my node
    do hat_on_node = 1, pub_num_hub_atoms_on_node(pub_my_node_id)

       hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)

       ! Atom number of Hubbard atom
       theatom = pub_distr_atom(hub%orig(hat))
       ! ddor: Number of NGWFs on this Hubbard atom
       olaps = ngwf_basis%num_on_atom(theatom)
       ! ddor: Number of Hubbard projectors on this Hubbard atom
       projs = hub_proj_basis%num_on_atom(theatom)

       !==============================================================
       ! LOCATE THE MOST LOCALISED NGWFS
       !==============================================================

       ! ddor: Initially set output to zero
       hub%localisedngwfs(:,hat_on_node) = 0
       ngwfsquareoverlapsum = 0.0_DP

       ! ddor: First do a loop over local row ngwfs
       !       in the hydrogenic-ngwf overlap matrix
       !       to see which NGWFs \psi on the Hubbard atom
       !       maximise the
       !       sum_{localised atomic projectors \phi} | < \phi | \psi > |**2

       ! ddor: loop over the NGWFs on this Hubbard atom
       do ngwf_count = 1, olaps

          ngwf_we_are_on = ngwf_count + ngwf_basis%first_on_atom(theatom) - 1

          do hub_proj_count = 1, projs

             hub_proj_we_are_on = hub_proj_count + &
                  hub_proj_basis%first_on_atom(theatom) - 1

             ! ddor: Add up the square of the overlap between
             !       the current NGWF and
             !       the contravariant projectors.
             !       As well as correcting for nonorthogonality,
             !       this fixes the normalisation of the NGWF projectors.
             call sparse_get_element(overlap_element,hub_consist_matrix,&
                  ngwf_we_are_on, hub_proj_we_are_on)
             call sparse_get_element(transpose_overlap_element,&
                  hub_consist_tmatrix,hub_proj_we_are_on, ngwf_we_are_on)

             ! ddor: The elements on the trace of the
             !       \sum_i <NGWF_|\phi_i><\phi_i|NGWF^> matrix
             ngwfsquareovlp(ngwf_count,hub_proj_count,hat) = &
                  overlap_element * transpose_overlap_element
             ngwfsquareoverlapsum(ngwf_count) = &
                  ngwfsquareoverlapsum(ngwf_count) + &
                  ngwfsquareovlp(ngwf_count,hub_proj_count,hat)

          end do

       end do

       ! ddor: Find out which NGWFs produce the greatest sum of
       !       square overlaps with the atomic projectors
       projectortest: do kk = 1, projs

!CW
          ! ddor: find the location of the maximum overlap**2
          ll = MAXLOC(ngwfsquareoverlapsum(1:olaps) - (/( pub_dmft_order_proj*dble(ijk),ijk=1,olaps )/) )
!END CW
          ! ddor: take note of the NGWF number of the projector
          hub%localisedngwfs(kk,hat_on_node) = ll(1) + &
               ngwf_basis%first_on_atom(theatom) - 1

          ! ddor: destroy this overlap
          ngwfsquareoverlapsum(ll(1)) = 0.0_DP

       end do projectortest

!CW
!CORRECTIONS DAVID OREGAN
       ! ddor: Assuming that projectors are always grouped together,
       !       re-sort them to preserve magnetic quantum number 
       inquire(file='mask_ngwfs_force',exist=ngwf_order)
       if(ngwf_order) then
         open(unit=115,file='mask_ngwfs_force')
          mm(1) = MINVAL(hub%localisedngwfs(1:projs,hat_on_node))
          do jjj = 1, projs
            read(115,*) jj
            if(pub_my_node_id==0) write(*,*) 'FORCING PROJECTION ',jj
            hub%localisedngwfs(jjj,hat_on_node) = mm(1) + jj - 1
          end do
         close(115)
       endif
       if(.not.pub_dmft_switch_off_proj_order)then
         mm(1) = MINVAL(hub%localisedngwfs(1:projs,hat_on_node))
         do jj = 1, projs
            hub%localisedngwfs(jj,hat_on_node) = mm(1) + jj - 1
         end do
       endif
!END CW
    end do ! ddor: End loop over Hubbard atoms on my node

    ! ddor: Reduce, Allow other nodes to catch up
    call comms_reduce('SUM',ngwfsquareovlp)

    call comms_barrier
    ! ddor: Write out the | < \phi | \psi > |**2 matrix
    if (pub_on_root .and. pub_output_detail == VERBOSE) &
         &call internal_ngwf_character_info

    ! Deallocate temporary storage
    deallocate(ngwfsquareoverlapsum,stat=ierr)
    call utils_dealloc_check('hubbard_ngwf_select','ngwfsquareoverlapsum',ierr)
    deallocate(ngwfsquareovlp,stat=ierr)
    call utils_dealloc_check('hubbard_ngwf_select','ngwfsquareovlp',ierr)
    deallocate(ll,stat=ierr)
    call utils_dealloc_check('hubbard_ngwf_select','ll',ierr)

    call timer_clock('hubbard_ngwf_select',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_ngwf_select'
#endif

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_ngwf_character_info

      implicit none

      integer :: species, this_species_atom ! ddor: Loop over species
      integer :: project, jj
      character(10) :: fmt, tmp

      do species=1,pub_cell%num_hub_species

         this_species_atom = 0

         atomloop : do hat=1,pub_cell%nat_hub

            if (hub%species_number(hat) == species) then

               this_species_atom = this_species_atom + 1
               theatom = pub_distr_atom(hub%orig(hat))
               if (this_species_atom .gt. 2) exit atomloop

               write(stdout,'(/a)')'########################################&
                    &########################################'

               ! ddor: TELL ME ABOUT THIS ATOM
               write(stdout,'(a,i6,a,a)') 'DFT+U information on &
                    &atom ',this_species_atom,' of Hubbard species ',&
                    h_species(species)%hub_species

               ! ddor: Number of NGWFs on this Hubbard atom
               olaps = ngwf_basis%num_on_atom(theatom)

               ! ddor: Number of Hubbard projectors on this Hubbard atom
               project = hub_proj_basis%num_on_atom(theatom)

               write(stdout,'(a)')'########################################&
                    &########################################'

               ! ddor: WRITE OUT INFORMATION ON | < \phi | \psi > |**2 MATRIX
               write(stdout,'(a)') 'Using the on-atom contravariant metric, &
                    &the NGWF hydrogenic character'
               write(stdout,'(a,i6,a)') '| < \phi | \psi > |**2 &
                    &matrix of Hubbard site ',hat,' is '

               ! ndmh: write header line
               write(stdout,'(a)',advance='no') '  m_l ='
               do jj=-(project-1)/2,(project-1)/2
                  write(stdout,'(i5,7x)',advance='no') jj
               end do
               write(stdout,'(a)') '  SUM'

               ! ndmh: write format string
               write(tmp,'(i6)') project+1
               write(fmt,'(3a)') '(',trim(adjustl(tmp)),'f12.8)'

               ! ndmh: write matrix
               do jj=1,olaps
                  write(stdout,fmt) ngwfsquareovlp(jj,1:project,hat), &
                       SUM(ngwfsquareovlp(jj,1:project,hat))
               end do

               write(stdout,'(a,f12.8)') 'Total hydrogenic projection of &
                    &NGWFs on this atom           ', &
                    &SUM(ngwfsquareovlp(1:olaps,1:project,hat))

               write(stdout,'(a/)')'########################################&
                    &########################################'

            end if

         end do atomloop

      end do

    end subroutine internal_ngwf_character_info

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine hubbard_ngwf_select


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_tensorial_correction(hub,ngwf_basis,proj_basis, &
       consistency,inv_overlap,loc_o_matrix)

    !==========================================================!
    ! This subroutine calculates for a given set of NGWFs the  !
    ! appropriate tensorial correction matrix, which is stored !
    ! in hubbard_o_matrix or hubbard_consist_o_matrix.         !
    !----------------------------------------------------------!
    ! Written by David O'Regan in April 2009.                  !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: pub_my_node_id, comms_abort, &
         pub_total_num_nodes, pub_on_root
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: pub_num_hub_atoms_on_node, &
         pub_hub_atoms_on_node, pub_distr_atom
    use rundat, only: pub_output_detail, pub_hub_tensor_corr
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_get_element, sparse_put_element, &
         sparse_scale
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use wrappers, only: wrappers_invert_sym_matrix

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    logical, intent(in) :: consistency
    type(SPAM3), intent(in) :: inv_overlap
    type(SPAM3), intent(inout) :: loc_o_matrix

    ! Local Variables
    integer :: hat_on_node, hat, theatom
    integer :: iii, jjj ! Used for tensorial correction
    integer :: num_on_atom_row, num_on_atom_col
    integer :: ierr
    real(kind=DP), allocatable :: o_matrix(:,:)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_tensorial_correction'
#endif

    call timer_clock('hubbard_tensorial_correction',1)

    if (consistency) then
       num_on_atom_row = ngwf_basis%max_on_atom
    else
       num_on_atom_row = proj_basis%max_on_atom
    end if

    allocate(o_matrix(num_on_atom_row,num_on_atom_row),stat=ierr)
    call utils_alloc_check('hubbard_tensorial_correction','o_matrix',ierr)

    call sparse_scale(hub%o_matrix,0.0_DP)

    ! ddor: Loop over Hubbard atoms on my node
    do hat_on_node = 1, pub_num_hub_atoms_on_node(pub_my_node_id)

       hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)

       theatom = pub_distr_atom(hub%orig(hat)) !Atom number of Hubbard atom

       !==============================================================
       ! CALCULATE THE CONTRAVARIANT METRIC ON EACH HUBBARD SITE
       !==============================================================

       ! ddor: If using nonorthogonal Hubbard projectors,
       !       calculate the necessary tensorial correction matrix
       !       tensorial: if (tensor) then

       o_matrix(:,:) = 0.0_DP

       if (consistency) then
          num_on_atom_row = ngwf_basis%num_on_atom(theatom)
          num_on_atom_col = ngwf_basis%num_on_atom(theatom)
       else
          num_on_atom_row = proj_basis%num_on_atom(theatom)
          num_on_atom_col = proj_basis%num_on_atom(theatom)
       end if

       ! ddor: Loop over column Hubbard projectors
       tensor_col_ngwf_loop1: do jjj = 1, num_on_atom_col

          ! ddor: Loop over row Hubbard projectors
          tensor_row_ngwf_loop1: do iii = 1,  num_on_atom_row

             ! ddor: The tensorial correction for all the NGWFs on an atom
             if (consistency) then

                ! ddor: ONLY RETRIEVE ELEMENTS FROM hub%consistency_overlap WHEN
                !       FORMING PROJECTOR DUALS IN THE CORRELATED SUBSPACE
                if (pub_hub_tensor_corr == 1) then

                   call sparse_get_element(&
                        o_matrix(iii,jjj), hub%consistency_overlap,&
                        iii + ngwf_basis%first_on_atom(theatom) - 1,&
                        jjj + ngwf_basis%first_on_atom(theatom) - 1)

                   ! ddor: WHEN FORMING PROJECTOR DUALS IN THE SIMULATION CELL,
                   !       USE FULL OVERLAP INVERSE FOR THE TENSOR CORRECTION.
                else if (pub_hub_tensor_corr == 2) then

                   call sparse_get_element(&
                        o_matrix(iii,jjj), inv_overlap,&
                        iii + ngwf_basis%first_on_atom(theatom) - 1,&
                        jjj + ngwf_basis%first_on_atom(theatom) - 1)

                   ! ddor: WHEN NOT PERFORMING TENSORIAL CORRECTION,
                   !       JUST USE A KRONECKER DELTA FUNCTION,
                   !       WITH A SCALING FACTOR WHICH ASSUMES ORTHOGONALITY
                   !       (TOTAL NUMBER OF CORRELATED ELECTRONS NOT PRESERVED).
                else if (pub_hub_tensor_corr == 3) then

                   if (iii .eq. jjj) then

                      call sparse_get_element(&
                           o_matrix(iii,jjj), hub%consistency_overlap,&
                           iii + ngwf_basis%first_on_atom(theatom) - 1,&
                           jjj + ngwf_basis%first_on_atom(theatom) - 1)

                      o_matrix(iii,jjj) = o_matrix(iii,jjj)**(-1.0_DP)
                   else
                      o_matrix(iii,jjj) = 0.0_DP
                   end if

                else if (pub_hub_tensor_corr .gt. 3) then
                   if (pub_on_root) write(stdout,'(a)') &
                        &'Invalid value for pub_hub_tensor_corr,&
                        & aborting.'
                   call comms_abort

                end if

             else

                ! ddor: ONLY RETRIEVE ELEMENTS FROM hub%consistency_overlap WHEN
                !       FORMING PROJECTOR DUALS IN THE CORRELATED SUBSPACE
                if (pub_hub_tensor_corr == 1) then

                   call sparse_get_element(&
                        o_matrix(iii,jjj), hub%consistency_overlap,&
                        hub%localisedngwfs(iii,hat_on_node),&
                        hub%localisedngwfs(jjj,hat_on_node))

                   ! ddor: WHEN FORMING PROJECTOR DUALS IN THE SIMULATION CELL,
                   !       USE FULL OVERLAP INVERSE FOR THE TENSOR CORRECTION.
                else if (pub_hub_tensor_corr == 2) then

                   call sparse_get_element(&
                        o_matrix(iii,jjj), inv_overlap,&
                        hub%localisedngwfs(iii,hat_on_node),&
                        hub%localisedngwfs(jjj,hat_on_node))

                   ! ddor: WHEN NOT PERFORMING TENSORIAL CORRECTION,
                   !       JUST USE A KRONECKER DELTA FUNCTION.
                else if (pub_hub_tensor_corr == 3) then

                   if (iii .eq. jjj) then

                      call sparse_get_element(&
                           o_matrix(iii,jjj), hub%consistency_overlap,&
                           hub%localisedngwfs(iii,hat_on_node),&
                           hub%localisedngwfs(jjj,hat_on_node))

                      o_matrix(iii,jjj) = o_matrix(iii,jjj)**(-1.0_DP)
                   else
                      o_matrix(iii,jjj) = 0.0_DP
                   end if

                else if (pub_hub_tensor_corr .gt. 3) then
                   write(stdout,*) 'Invalid value for pub_hub_tensor_corr, &
                        &aborting.'
                   call comms_abort

                end if

             end if

          end do tensor_row_ngwf_loop1

       end do tensor_col_ngwf_loop1

       !write(stdout,*) 'Covariant metric on correlated space on site ',&
       !     hat,' is '
       !do mm=1,num_on_atom_col
       !   write(*,*) hub%o_matrix(mm,1:num_on_atom_row)
       !end do

       ! ddor: ONLY INVERT PROJECTOR OVERLAP MATRIX WHEN
       !       FORMING PROJECTOR DUALS IN THE CORRELATED SUBSPACE
       if (pub_hub_tensor_corr == 1) then
          ! ddor: If we want the tensorial correction for all
          !       the NGWFs on an atom
          if (consistency) then
             ! ddor: Perform matrix inverse to find the
             !       contravariant metric on the correlated subspace
             call wrappers_invert_sym_matrix(o_matrix(:,:),&
                  ngwf_basis%num_on_atom(theatom) )
          else

             ! ddor: Perform matrix inverse to find the
             !       contravariant metric on the correlated subspace
             call wrappers_invert_sym_matrix(o_matrix(:,:),&
                  proj_basis%num_on_atom(theatom) )
          end if
       end if

       ! ddor: Loop over column Hubbard projectors
       tensor_col_ngwf_loop2: do jjj = 1,  num_on_atom_col

          ! ddor: Loop over row Hubbard projectors
          tensor_row_ngwf_loop2: do iii = 1, num_on_atom_row

             ! ddor: If we want the tensorial correction for all
             !       the NGWFs on an atom
             if (consistency) then

                ! ddor: Get the corresponding element of old S
                !       and put it in hub_consist_o_matrix
                call sparse_put_element(&
                     o_matrix(iii,jjj), loc_o_matrix,&
                     iii + ngwf_basis%first_on_atom(theatom) - 1,&
                     jjj + ngwf_basis%first_on_atom(theatom) - 1)

             else

                ! ddor: Get the corresponding element of old S
                !       and put it in hub%o_matrix
                call sparse_put_element(&
                     o_matrix(iii,jjj), loc_o_matrix,&
                     iii + proj_basis%first_on_atom(theatom) - 1,&
                     jjj + proj_basis%first_on_atom(theatom) - 1)

             end if

          end do tensor_row_ngwf_loop2

       end do tensor_col_ngwf_loop2

       !write(stdout,*) 'Contravariant metric on correlated space on site ',&
       !     hat,' is '
       !do mm=1,num_on_atom_col
       !   write(*,*) o_matrix(mm,1:num_on_atom_row)
       !end do

    end do ! ddor: End loop over Hubbard atoms on my node

    deallocate(o_matrix,stat=ierr)
    call utils_dealloc_check('hubbard_tensorial_correction','o_matrix',ierr)

    call timer_clock('hubbard_tensorial_correction',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_tensorial_correction'
#endif

  end subroutine hubbard_tensorial_correction


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_ngwf_on_grid(hub,ngwf_basis,hub_proj_basis)

    !==========================================================!
    ! This subroutine copies the NGWFs which have been chosen  !
    ! as Hubbard projectors into an array which is indexed in  !
    ! the same order as how the <NGWF|projector> matrix        !
    ! elements are sorted.                                     !
    !----------------------------------------------------------!
    ! Written by David O'Regan in September 2009.              !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: pub_my_node_id, pub_on_root, &
         comms_abort
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_brappd_ketppd
    use parallel_strategy, only: pub_first_atom_on_node, &
         pub_num_atoms_on_node, pub_num_hub_atoms_on_node, &
         pub_hub_atoms_on_node, pub_distr_atom
    use simulation_cell, only: pub_cell
    use sparse, only : SPAM3, sparse_create, sparse_destroy
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis

    ! Local Variables
    integer :: ngwf_count
    integer :: ngwf_offset
    integer :: hubbard_offset
    integer :: loc_iat
    integer :: iat
    integer :: hat_on_node, hat
    integer :: ngwf, ngwf_skip
    integer :: ngwf_fix_up, ngwf_fix_down
    integer :: hub_proj
    integer :: index, projector_index, ngwf_index

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_ngwf_on_grid'
#endif

    call timer_clock('hubbard_ngwf_on_grid',1)

    hub%consistency_projs = 0.0_DP

    ngwf_count = 1
    ngwf_offset = 1
    hubbard_offset = 1

    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1

       do hat_on_node=1,pub_num_hub_atoms_on_node(pub_my_node_id)

          hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
          if (iat .eq. pub_distr_atom(hub%orig(hat)) ) then

             do hub_proj=1,hub_proj_basis%num_on_atom(iat)

                do ngwf=ngwf_basis%first_on_atom(iat),&
                     ngwf_basis%first_on_atom(iat)+ngwf_basis%num_on_atom(iat)-1

                   if (hub%localisedngwfs(hub_proj,hat_on_node) .eq. ngwf) then

                      do ngwf_fix_up=1,(ngwf - ngwf_basis%first_on_atom(iat))

                         ngwf_offset = ngwf_offset + &
                              pub_cell%n_pts*ngwf_basis%spheres(&
                              &ngwf_count)%n_ppds_sphere

                         ngwf_count = ngwf_count + 1

                      end do

                      do index=1,pub_cell%n_pts*&
                           &ngwf_basis%spheres(ngwf_count)%n_ppds_sphere

                         projector_index = index + hubbard_offset - 1
                         ngwf_index = index + ngwf_offset - 1

                         hub%consistency_projs(projector_index) = &
                              &hub%consistency_ngwfs(ngwf_index)

                      end do

                      hubbard_offset = hubbard_offset + &
                           &pub_cell%n_pts*&
                           &ngwf_basis%spheres(ngwf_count)%n_ppds_sphere

                      do ngwf_fix_down =1, &
                           (ngwf - ngwf_basis%first_on_atom(iat))

                         ngwf_count = ngwf_count - 1

                         ngwf_offset = ngwf_offset - &
                              pub_cell%n_pts*ngwf_basis%spheres(&
                              &ngwf_count)%n_ppds_sphere

                      end do

                   end if

                end do

             end do

          end if

       end do

       do ngwf_skip=1,ngwf_basis%num_on_atom(iat)

          if (ngwf_count .lt. ngwf_basis%num) then
             ngwf_offset = ngwf_offset + &
                  pub_cell%n_pts*ngwf_basis%spheres(ngwf_count)%n_ppds_sphere
          elseif (ngwf_count .gt. ngwf_basis%num) then
             if (pub_on_root) write(stdout,'(a)') &
                  &'Something wrong in hubbard_ngwf_on_grid, aborting'
             call comms_abort
          end if

          ngwf_count = ngwf_count + 1

       end do

    end do

    call timer_clock('hubbard_ngwf_on_grid',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_ngwf_on_grid'
#endif

  end subroutine hubbard_ngwf_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_projection_mtx(tensor,consistency, &
       hub_overlap, hub_overlap_t, &
       hub_o_matrix, hub_consist_matrix, hub_consist_tmatrix)

    !==========================================================!
    ! This matrix takes the  <NGWF|projector> matrix and the   !
    ! tensorial correction matrix if (tensor == .true.) and    !
    ! from then calculates the SPAM3 Hubbard projection.       !
    !----------------------------------------------------------!
    ! Written by David O'Regan in April 2009.                  !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: pub_on_root
    use hubbard_init, only: h_species
    use sparse, only: SPAM3, sparse_create, sparse_destroy, &
         sparse_transpose, sparse_product
    use timer, only: timer_clock

    implicit none

    ! Arguments
    ! ddor: Apply tensorial correction for non-orthogonal Hubbard projectors
    logical, intent(in) :: tensor
    ! ddor: Apply tensorial correction for non-orthogonal Hubbard projectors
    !       for the sake of determining which NGWFs to use as projectors
    logical, intent(in) :: consistency
    type(SPAM3), intent(in) :: hub_overlap
    type(SPAM3), intent(inout) :: hub_overlap_t
    type(SPAM3), intent(in) :: hub_o_matrix
    type(SPAM3), optional, intent(in) :: hub_consist_matrix
    type(SPAM3), optional, intent(inout) :: hub_consist_tmatrix

    ! Local Variables
    type(SPAM3) :: temp_overlap_t

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_projection_mtx'
#endif

    call timer_clock('hubbard_projection_mtx',1)

    ! ddor: generate hub_consist_tmatrix instead of hub_overlap_t
    if (consistency) then

       ! ddor: Carry out tensorial correction for nonorthogonal projectors
       if (tensor) then
          ! ddor: Raise projector index in <Projector_NGWF> matrix
          !       Here the NGWF index is raised, not the projector, so the
          !       correction matrix appears on the right.
          call sparse_create(temp_overlap_t,hub_consist_tmatrix)
          call sparse_transpose(temp_overlap_t,hub_consist_matrix)
          call sparse_product(hub_consist_tmatrix,temp_overlap_t,hub_o_matrix)
          call sparse_destroy(temp_overlap_t)
       else
          call sparse_transpose(hub_consist_tmatrix,hub_consist_matrix)
       end if

    else

       ! ddor: Carry out tensorial correction for nonorthogonal projectors
       if (tensor) then
          ! ddor: Raise projector index in <Projector_NGWF> matrix
          call sparse_create(temp_overlap_t,hub_overlap_t)
          call sparse_transpose(temp_overlap_t,hub_overlap)
          call sparse_product(hub_overlap_t,hub_o_matrix,temp_overlap_t)
          call sparse_destroy(temp_overlap_t)
       else
          call sparse_transpose(hub_overlap_t,hub_overlap)
       end if

    end if

    call timer_clock('hubbard_projection_mtx',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_projection_mtx'
#endif

  end subroutine hubbard_projection_mtx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_projector_character(hub,hub_proj_basis, &
       hub_consist_matrix)

    !==========================================================!
    ! This subroutine recomputes the hydrogenic character of   !
    ! NGWF DFT+U projectors using the updated                  !
    ! tensorial correction.                                    !
    !----------------------------------------------------------!
    ! Written by David O'Regan in October 2009.                !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: pub_my_node_id, comms_abort, comms_barrier, &
         comms_reduce, pub_on_root
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: pub_num_hub_atoms_on_node, &
         pub_hub_atoms_on_node, pub_distr_atom
    use rundat, only: pub_output_detail
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_get_element, sparse_product
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(SPAM3), intent(in) :: hub_consist_matrix

    ! Local Variables
    integer :: hat_on_node, hat, theatom
    integer :: projs
    integer :: ierr
    integer :: hub_proj_row_we_are_on, hub_proj_col_we_are_on
    integer :: hub_proj_count
    integer :: hub_hydro_we_are_on
    integer :: hub_tensorial_index_row, hub_tensorial_index_col
    integer :: hub_hydro_count, hub_tensorial_count
    real(kind=DP) :: overlap_element
    real(kind=DP) :: transpose_overlap_element, tensorial_element
    real(kind=DP), allocatable :: ngwfsquareovlp(:,:,:)


#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_projector_character'
#endif

    call timer_clock('hubbard_projector_character',1)

    allocate(ngwfsquareovlp(hub_proj_basis%max_on_atom, &
         hub_proj_basis%max_on_atom,pub_cell%nat_hub),stat=ierr)
    call utils_alloc_check('hubbard_ngwf_select','ngwfsquareovlp',ierr)
    ngwfsquareovlp = 0.0_DP

!CW BUG CORRECTED
    hub%current_projection=0.0_DP
!END CW

    ! ddor: Loop over Hubbard atoms on my node
    do hat_on_node = 1, pub_num_hub_atoms_on_node(pub_my_node_id)

       hat = pub_hub_atoms_on_node(hat_on_node,pub_my_node_id)
       theatom = pub_distr_atom(hub%orig(hat)) !Atom number of Hubbard atom
       ! ddor: Number of Hubbard projectors on this Hubbard atom
       projs = hub_proj_basis%num_on_atom(theatom)

       ! ddor: loop over the row projectors on this Hubbard atom
       do hub_proj_count = 1, projs

          ! ddor: The rows of <NGWF_|phi> and columns of <NGWF^|NGWF^>
          !       are indexed differently
          hub_proj_row_we_are_on = &
               hub%localisedngwfs(hub_proj_count,hat_on_node)
          hub_proj_col_we_are_on = hub_proj_count + &
               ( hub_proj_basis%first_on_atom(theatom) - 1 )

          ! ddor: loop over the hydrogenic projectors on this Hubbard atom
          do hub_hydro_count = 1, projs

             hub_hydro_we_are_on = hub_hydro_count + &
                  hub_proj_basis%first_on_atom(theatom) - 1

             ! ddor: Add up the square of the overlap between
             !       the current NGWF and
             !       the contravariant projectors.
             !       As well as correcting for nonorthogonality,
             !       this fixes the normalisation of the NGWF projectors.
             call sparse_get_element(overlap_element,hub_consist_matrix,&
                  hub_proj_row_we_are_on, hub_hydro_we_are_on)

             ! ddor: Loop over rows of the contravariant matrix on the site
             do hub_tensorial_count = 1, projs

                hub_tensorial_index_col = &
                     hub%localisedngwfs(hub_tensorial_count,hat_on_node)

                hub_tensorial_index_row = &
                     hub_tensorial_count + &
                     hub_proj_basis%first_on_atom(theatom) - 1

                ! ddor: Hydrogenic rows, NGWF projector columns,
                !       but use the transpose
                call sparse_get_element(transpose_overlap_element,&
                     hub_consist_matrix,&
                     hub_tensorial_index_col, hub_hydro_we_are_on)

                call sparse_get_element(tensorial_element,&
                     hub%o_matrix,&
                     hub_tensorial_index_row, hub_proj_col_we_are_on)

                ! ddor: The elements (alpha, i)
                !                 { <NGWF_alpha|\phi_i>    }
                !       sum_beta  { <\phi_i|NGWF_beta>     }
                !                 { <NGWF^beta|NGWF^alpha> }
                ngwfsquareovlp(hub_proj_count,&
                     hub_hydro_count,hat) = &
                     ngwfsquareovlp(hub_proj_count,&
                     hub_hydro_count,hat) + &
                     overlap_element * transpose_overlap_element * &
                     tensorial_element

             end do

          end do

       end do

       hub%current_projection(hat) = &
            SUM(ngwfsquareovlp(1:projs,1:projs,hat))

    end do ! ddor: End loop over Hubbard atoms on my node

    ! ddor: Reduce, Allow other nodes to catch up
    call comms_reduce('SUM',ngwfsquareovlp)
    call comms_reduce('SUM',hub%current_projection)

    call comms_barrier

    ! ddor: Write out the | < \phi | \psi > |**2 matrix
    if (pub_on_root .and. pub_output_detail == VERBOSE) &
         &call internal_ngwf_character_info2

    deallocate(ngwfsquareovlp,stat=ierr)
    call utils_dealloc_check('hubbard_ngwf_select','ngwfsquareovlp',ierr)

    call timer_clock('hubbard_projector_character',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_projector_character'
#endif

  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_ngwf_character_info2

      implicit none

      ! Local Variables
      integer :: project, jj
      character(10) :: fmt, tmp

      do hat=1,pub_cell%nat_hub

         ! ddor: Number of Hubbard projectors on this Hubbard atom
         theatom = pub_distr_atom(hub%orig(hat))
         project = hub_proj_basis%num_on_atom(theatom)

         write(stdout,'(/a)')'########################################&
              &########################################'

         ! ddor: WRITE OUT INFORMATION ON | < \phi | \psi > |**2 MATRIX
         write(stdout,'(a)') 'Using the subspace &
              &contravariant metric, &
              &the NGWF hydrogenic character'
         write(stdout,'(a,i6,a)') '| < \phi | \psi > |**2 &
              &matrix of Hubbard site ',hat,' is '

         ! ndmh: write header line
         write(stdout,'(a)',advance='no') '  m_l ='
         do jj=-(project-1)/2,(project-1)/2
            write(stdout,'(i5,7x)',advance='no') jj
         end do
         write(stdout,'(a)') '  SUM'

         ! ndmh: write format string
         write(tmp,'(i6)') project+1
         write(fmt,'(3a)') '(',trim(adjustl(tmp)),'f12.8)'

         ! ndmh: write matrix
         do jj=1,project
            write(stdout,fmt) ngwfsquareovlp(jj,1:project,hat), &
                 SUM(ngwfsquareovlp(jj,1:project,hat))
         end do

         write(stdout,'(a,f12.8)') 'Total hydrogenic projection of &
              &NGWFs on this correlated site', hub%current_projection(hat)

         write(stdout,'(a/)')'########################################&
              &########################################'

      end do

    end subroutine internal_ngwf_character_info2

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine hubbard_projector_character


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_projector_consistency( &
       elements, ngwfs_on_grid, ngwf_basis, projector_basis, hub_proj_basis, &
       hub, overlap, inv_overlap, hub_overlap, hub_overlap_t, &
       sp_overlap, hub_proj_paw_overlap, hub_ngwf_paw_overlap)

    !==========================================================!
    !  Depending on the projector iteration we are on,         !
    !  carry out some projector and metric mixing,             !
    !  and generate the new projector-NGWF overlap matrix.     !
    !----------------------------------------------------------!
    ! Written by David O'Regan in July 2009.                   !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use augmentation, only: augmentation_overlap
    use comms, only: pub_on_root
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT
    use ion, only: element
    use integrals, only: integrals_brappd_ketppd
    use parallel_strategy, only: pub_elements_on_node
    use paw, only: paw_projectors
    use projectors, only: projectors_func_ovlp_box
    use restart, only : restart_ngwfs_tightbox_output, &
         restart_ngwfs_tightbox_input
    use rundat, only: pub_hub_proj_mixing, pub_output_detail, &
         write_tightbox_ngwfs, maxit_ngwf_cg, &
         pub_hubbard_atomsolve, pub_aug
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_scale, sparse_axpy,&
         sparse_copy, sparse_create, sparse_destroy
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(ELEMENT), intent(in)       :: elements(pub_cell%nat)
    type(FUNC_BASIS), intent(in)    :: ngwf_basis
    type(FUNC_BASIS), intent(in)    :: projector_basis
    type(FUNC_BASIS), intent(inout) :: hub_proj_basis
    real(kind=DP), intent(in)       :: ngwfs_on_grid(&
         &ngwf_basis%n_ppds*pub_cell%n_pts)
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(SPAM3), intent(in)         :: overlap
    type(SPAM3), intent(in)         :: inv_overlap
    type(SPAM3), intent(inout)      :: hub_overlap
    type(SPAM3), intent(inout)      :: hub_overlap_t
    ! ddor: Needed for PAW+U
    type(SPAM3), intent(in)         :: sp_overlap
    type(SPAM3), intent(inout)      :: hub_proj_paw_overlap
    type(SPAM3), intent(inout)      :: hub_ngwf_paw_overlap

    ! Local Variables
    ! ddor: Contravariant <projector|projector> matrix for tensorial correction
    !       In this case it is for all the NGWFs on a Hubbard atom, so that
    !       the most localised ones can be chosen in tensorially consistent way.
    type(SPAM3) :: hub_consist_o_matrix
    ! ddor: <NGWF|projector> Matrix used to determine which NGWFs
    !        to use as projectors
    type(SPAM3) :: hub_consist_matrix
    ! ddor: Corrected <projector|NGWF> Matrix used to determine
    !       which NGWFs to use as projectors
    type(SPAM3) :: hub_consist_tmatrix

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_projector_consistency'
#endif

    call timer_clock('hubbard_projector_consistency',1)

    ! ddor: On first projector optimisation iteration do a plain
    !       ONETEP calculation, unless we have projectors stored in
    !       tightbox representation, in which case we go straight to
    !       DFT+U with NGWF projectors.

    if (( hub%consistency_iteration > 1) .or. &
        & (pub_hub_proj_mixing .lt. 0.0_DP) .or. &
         & (pub_hubbard_atomsolve) ) then

       ! ddor: If necessary, mix in some of previous projectors
       !       to aid convergence
       if ( ( ( pub_hub_proj_mixing .eq. 0.0_DP ) .and. &
            & (.not. (pub_hubbard_atomsolve .and. &
            & (hub%consistency_iteration >= 2) ) ) ) &
            & .or. ((pub_hub_proj_mixing .gt. 0.0_DP).and. &
            & ( hub%consistency_iteration == 2 )) &
            & .or. ((pub_hub_proj_mixing .lt. 0.0_DP).and. &
            & ( hub%consistency_iteration >= 2 )) ) then

          ! ddor: Nothing to mix
          hub%consistency_ngwfs = ngwfs_on_grid

          ! ddor: The old overlap matrix is needed to construct the
          !       covariant metric on the correlated sites
          call sparse_copy(hub%consistency_overlap,overlap)

          ! ddor: Restore normal value
!CW
! FOR DMFT
!         if (pub_hub_proj_mixing .lt. 0.0_DP) then
          if ((pub_hub_proj_mixing .lt. 0.0_DP).and.(pub_hub_proj_mixing>-1.99_DP)) then
!END CW
             maxit_ngwf_cg = hub%store_maxit_ngwf_cg
          end if

       else if (( pub_hub_proj_mixing .lt. 0.0_DP ) .and. &
            ( hub%consistency_iteration .eq. 1)) then

          ! ddor: If mixing is negative, take this to
          !       mean read Hubbard NGWF projectors from
          !       file in tightbox representation
          call restart_ngwfs_tightbox_input(&
               hub%consistency_ngwfs,ngwf_basis,elements,'tightbox_hub_projs')

          ! ddor: Overlap of Hubbard NGWF projectors from file
          call integrals_brappd_ketppd(hub%consistency_overlap, & !output
               hub%consistency_ngwfs, ngwf_basis, &               !input
               hub%consistency_ngwfs, ngwf_basis)                 !input

          ! ddor: Augment the <NGWF|NGWF> matrix for PAW+U or USP calculations
          if (pub_aug) then

             ! ddor: Calculate the overlap matrix of the Hubbard
             ! ddor: and PAW projectors
             call projectors_func_ovlp_box(hub_ngwf_paw_overlap, &
                  hub%consistency_ngwfs,ngwf_basis,projector_basis, &
                  paw_projectors)

             ! ndmh: Calculate the augmentation of the overlap matrix due to the
             ! ndmh: augmentation region part of the overlap operator
             call augmentation_overlap(hub%consistency_overlap,sp_overlap,&
                  hub_ngwf_paw_overlap)

          end if

          ! ddor: Do no NGWF iterations on first HUBBARDSCF iteration
          !       when restarting.
          hub%store_maxit_ngwf_cg = maxit_ngwf_cg
!CW
! added "if(...)" and "end if"
         if(pub_hub_proj_mixing>-1.99_DP)then
           maxit_ngwf_cg = 0
         endif
!END CW
       else

          ! ddor: Mix the projectors for iterations > 2
          hub%consistency_ngwfs = &
               ( ABS(pub_hub_proj_mixing) * hub%consistency_ngwfs ) + &
               ( ( 1.0_DP - ABS(pub_hub_proj_mixing) ) * ngwfs_on_grid )

          ! ddor: Mix the convariant metrics for iterations > 2
          !       Cheaper than calculating
          !       <hub%consistency_ngwfs|hub%consistency_ngwfs>
          call sparse_scale(hub%consistency_overlap,ABS(pub_hub_proj_mixing))
          call sparse_axpy(hub%consistency_overlap,&
               overlap,1.0_DP-ABS(pub_hub_proj_mixing))

       end if

       call sparse_create(hub_consist_matrix,hub_overlap)
       call sparse_create(hub_consist_tmatrix,hub_overlap_t)

       ! ddor: Renew hubbard_consistency_matrix
       !       (the overlap of NGWFs and atomic orbitals)
       !       needed to select which NGWFs have most localised character
       !       on this iteration.
       call projectors_func_ovlp_box(hub_consist_matrix, &
            ngwfs_on_grid,ngwf_basis,hub_proj_basis,hub%projectors)

       hub_consist_o_matrix%structure = 'D'
       call sparse_create(hub_consist_o_matrix)

       ! ddor: Find the tensorial correction matrix for all candidate
       !       NGWF projectors on Hubbard atoms
       call hubbard_tensorial_correction(hub,ngwf_basis,hub_proj_basis,.true.,&
            inv_overlap, hub_consist_o_matrix)

       ! ddor: With hub_consist_matrix and the tensorial correction matrix,
       !       calculate the Hubbard projection operator
       !       for all candidate projectors
       call hubbard_projection_mtx(.true.,.true.,hub_overlap,hub_overlap_t, &
            hub_consist_o_matrix, hub_consist_matrix, hub_consist_tmatrix)

       call sparse_destroy(hub_consist_o_matrix)

       ! ddor: Select the most localised NGWFs to use as Hubbard
       !       projectors on each correlated site.
       !       Using the hub_overlap from previous Hubbard iteration
       !       (e.g. hubbard_consistency_overlap) calculate
       !       the new tensorial correction matrix, the contravariant
       !       covariant metric on each correlated site.
       call hubbard_ngwf_select(hub, ngwf_basis, hub_proj_basis,  &
            hub_consist_matrix, hub_consist_tmatrix)

       ! ddor: Move the Hubbard NGWF projectors on the grid
       !       to an array containing the minimal set.
       call hubbard_ngwf_on_grid(hub,ngwf_basis, hub_proj_basis)

       ! ddor: Compute new NGWF-projector overlap matrix
       call integrals_brappd_ketppd(hub_overlap, &              !input-output
            ngwfs_on_grid, ngwf_basis, hub%consistency_projs, & !input
            hub_proj_basis)                                     !input

       ! ddor: Augment the <NGWF|proj> matrix for PAW+U or USP calculations
       if (pub_aug) then

          ! ddor: Calculate the overlap matrix of the Hubbard
          ! ddor: and PAW projectors
          call projectors_func_ovlp_box(hub_proj_paw_overlap, &
               hub%consistency_projs,hub_proj_basis,&
               projector_basis,paw_projectors)

          ! ndmh: Calculate the augmentation of the overlap matrix due to the
          ! ndmh: augmentation region part of the overlap operator
          call augmentation_overlap(hub_overlap,sp_overlap,&
               hub_proj_paw_overlap)

       end if

       ! ddor: write new Hubbard NGWF projectors to file in tightbox rep.
       !       write out all NGWFs, so we won't need projector tightboxes,
       !       or write out hub%consistency_overlap.
!CW
!       if (( hub%consistency_iteration .gt. 1) .and. write_tightbox_ngwfs) then
!FOR DMFT
!pub_hub_proj_mixing=-2 : read tightbox_hub but does not overwrite it
        if ( (pub_hub_proj_mixing>-1.99) .and. (hub%consistency_iteration .gt. 0) .and. write_tightbox_ngwfs) then
!END CW
          call restart_ngwfs_tightbox_output(&
               hub%consistency_ngwfs,ngwf_basis,elements,'tightbox_hub_projs')
       end if

       ! ddor: Find the tensorial correction matrix for
       !       the chosen NGWF Hubbard projectors
       call hubbard_tensorial_correction(hub,ngwf_basis,hub_proj_basis,&
            .false.,inv_overlap,hub%o_matrix)

       ! ddor: Recalculate the hydrogenic character of the NGWF projectors
       !       using the appropriate tensorial correction
       call hubbard_projector_character(hub,hub_proj_basis, &
            hub_consist_matrix)

       call sparse_destroy(hub_consist_tmatrix)
       call sparse_destroy(hub_consist_matrix)

       ! ddor: With hub_overlap and the tensorial correction matrix,
       !       calculate the Hubbard projection operator for each site
       call hubbard_projection_mtx(.true.,.false.,hub_overlap,hub_overlap_t, &
            hub%o_matrix)

    end if

    call timer_clock('hubbard_projector_consistency',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_projector_consistency'
#endif

  end subroutine hubbard_projector_consistency


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function hubbard_test_convergence(hub,total_energy)

    !==========================================================!
    ! When carrying out a DFT+U calculation with               !
    ! self-consistent NGWF projectors, this function tests     !
    ! the convergence by seeing if the energy difference       !
    ! between ground-state energies with is less than          !
    ! hub_energy_tol over hub_conv_win projector optimisation  !
    ! steps.                                                   !
    !                                                          !
    ! The self-consistent scheme also halts if the sum of the  !
    ! projections onto hydrogenic orbitals of the NGWF Hubbard !
    ! projectors remains constant over hub_conv_win projector  !
    ! optimisation steps.                                      !
    !----------------------------------------------------------!
    ! Written by David O'Regan in April 2009.                  !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use comms, only: pub_on_root
    use rundat, only: pub_hub_energy_tol,pub_hub_conv_win
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    real(kind=DP), intent(in) :: total_energy

    ! Local Variables
    integer ::  hub_energies_flag,  hub_projections_flag
    integer ::  hub_iterations_counter

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_test_convergence'
#endif

    ! ddor: Zero if the energy convergence criterion has been met
    hub_energies_flag = 0
    ! ddor: Zero if the hydrogenic projection convergence criterion has been met
    hub_projections_flag = 0

    ! ddor: Check the differences of historical energies and projections
    do hub_iterations_counter = 1, (pub_hub_conv_win - 1)

       ! ddor: Energies
       if (ABS( hub%consistency_energies(hub_iterations_counter) - &
            hub%consistency_energies(hub_iterations_counter + 1 ) )&
            & .gt. pub_hub_energy_tol ) then
          ! ddor: The difference of energies exceeds the threshold
          hub_energies_flag = 1
       end if
       ! ddor: Move energies down in the array
       hub%consistency_energies(hub_iterations_counter) = &
            hub%consistency_energies(hub_iterations_counter + 1 )

       ! ddor: Projections
       if (ABS( hub%consistency_projections(hub_iterations_counter) - &
            hub%consistency_projections(hub_iterations_counter + 1 ) )&
            & .gt. 1.0e-8_DP ) then
          ! ddor: The difference of projections exceeds the threshold
          hub_projections_flag = 1
       end if
       ! ddor: Move projections down in the array
       hub%consistency_projections(hub_iterations_counter) = &
            hub%consistency_projections(hub_iterations_counter + 1 )

    end do

    ! ddor: fill the top of the array with the current energy
    hub%consistency_energies(pub_hub_conv_win) = total_energy

    ! ddor: fill the top of the array with the current
    !       sum of hydrogenic projections of Hubbard NGWFs of all sites
    if (hub%consistency_iteration .gt. 1) then
       hub%consistency_projections(pub_hub_conv_win) = &
            SUM(hub%current_projection(1:pub_cell%nat_hub))
    end if

    ! ddor: Check the difference between current and previous energy
    if (ABS( hub%consistency_energies(pub_hub_conv_win - 1 ) - &
         hub%consistency_energies(pub_hub_conv_win) ) .gt. &
         &pub_hub_energy_tol ) then
       hub_energies_flag = 1
    end if

    ! ddor: Check the difference between current and previous projection
    if (ABS( hub%consistency_projections(pub_hub_conv_win - 1 ) - &
         hub%consistency_projections(pub_hub_conv_win ) ) .gt. &
         & 1.0e-8_DP ) then
       hub_projections_flag = 1
    end if

    ! ddor: If we still think that the energy has converged, and we've
    !       done sufficiently many iterations, then pass the convergence test
    if ( ( hub_projections_flag .eq. 0) &
         & .and. (hub%consistency_iteration .gt. pub_hub_conv_win+1)) then
       hubbard_test_convergence = .true.
       ! ddor: Write out a congratulatory message
       if (pub_on_root) then
          write(stdout,'(/a)')'########################################&
               &########################################'
          write(stdout,'(a,i6)')'DFT+U projectors optimised on iteration ', &
               hub%consistency_iteration
          if (hub_energies_flag .eq. 0) then
             write(stdout,'(a,f12.8,a,i6,a)') 'Energy tolerance of ', &
                  pub_hub_energy_tol,' maintained over ',pub_hub_conv_win, &
                  ' iterations'
          end if
          if (hub_projections_flag .eq. 0) then
             write(stdout,'(a,f12.8,a,i6,a)') 'Projector character &
                  &tolerance of ', 1.0e-8_DP,' maintained over ',&
                  pub_hub_conv_win, ' iterations'
          end if
          write(stdout,'(a/)')'########################################&
               &########################################'
       end if
    else
       ! ddor: Otherwise, keep going.
       hubbard_test_convergence = .false.
    end if

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &hubbard_test_convergence'
#endif

  end function hubbard_test_convergence


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_species_proj(hub,elements)

    !=================================================================!
    ! This subroutine allocates all the memory for                    !
    ! h_proj data on hydrogenic projectors used by DFT+U and          !
    ! calculates this data for each Hubbard species.                  !
    ! This subroutine allocates and initialises hub_fftbox_proj_recip !
    ! for each h_species element. Each such hub_fftbox_proj_recip is  !
    ! a full non-local projector in reciprocal space in the FFTbox.   !
    !-----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/1/2004.                  !
    ! Tidied up by Peter D. Haynes on 1/7/2004.                       !
    ! Modified for DFT+U module by David O'Regan in July 2008.        !
    ! Modified by Nicholas Hine in November 2011 to use the           !
    ! HUBBARD_MODEL type                                              !
    !=================================================================!

    use comms, only: pub_on_root,comms_abort
    use constants, only: DP, PI
    use fourier, only: fourier_apply_box
    use geometry, only: POINT, OPERATOR(.dot.)
    use hubbard_init, only: h_species
    use ion, only: element
    use parallel_strategy, only: pub_distr_atom
    use projectors, only: projectors_allocate_set, projector_in_box_recip, &
         projectors_init_fftbox_recip, projectors_exit_fftbox_recip
    use rundat, only: pub_output_detail
    use simulation_cell, only: pub_cell, pub_fftbox
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! Local Variables
    integer :: n1, n2, n3, ld1, ld2
    integer :: ps
    integer :: proj_count, hat
    integer :: iat, orig_iat
    real(kind=DP) :: spacing,deltar ! Grid variables
    type(POINT) :: kpt_loc

    ! ndmh: removed hubbard_proj type, so these variables need to be local
    ! ddor: number of points in the real radial grid
    integer :: n_r_rad_pts
    ! ddor: number of points in the reciprocal radial grid
    integer :: n_g_rad_pts, max_n_g_rad_pts
    ! ddor: maximum radial reciprocal k-vector up to which
    !      the hydrogenic projectors are defined
    real(kind=DP) :: gmax
    ! ddor: hydrogenic projector radius in real space
    real(kind=DP) :: rmax

    if (pub_on_root  .and. pub_output_detail == VERBOSE) &
         write(stdout,'(a)') &
         '... Hubbard projector initialisation'

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &hubbard_species_proj'
#endif

    call timer_clock('hubbard_species_proj',1)

    ! Local copies of FFT box dimensions
    n1 = pub_fftbox%total_pt1
    n2 = pub_fftbox%total_pt2
    n3 = pub_fftbox%total_pt3
    ld1 = pub_fftbox%total_ld1
    ld2 = pub_fftbox%total_ld2

    kpt_loc%x = 0.0_DP ; kpt_loc%y = 0.0_DP ; kpt_loc%z = 0.0_DP

    ! ndmh: find max value of recip grid points
    max_n_g_rad_pts = -1
    spacing = min(pub_fftbox%d1 ,pub_fftbox%d2, pub_fftbox%d3)
    do ps=1,pub_cell%num_hub_species
       do hat=1,pub_cell%nat_hub
          orig_iat = hub%orig(hat)
          if (hub%species_number(hat) == ps) exit
       end do
       rmax = elements(orig_iat)%radius
       n_r_rad_pts = int(100.0_DP*rmax/spacing)
       deltar = rmax / real(n_r_rad_pts - 1,kind=DP)
       n_g_rad_pts = nint(rmax/deltar) + 1
       max_n_g_rad_pts = max(max_n_g_rad_pts,n_g_rad_pts)
    end do

    ! ndmh: count unique projectors
    hub%projectors%n_proj_species = pub_cell%num_hub_species
    proj_count = 0
    do ps=1,pub_cell%num_hub_species
       proj_count = proj_count + h_species(ps)%hub_ang_mom*2 + 1
    end do

    ! ndmh: allocate hub%projectors type
    call projectors_allocate_set(hub%projectors,1,max_n_g_rad_pts)
    hub%projectors%normalise = .true.

    ! ndmh: set up entries in hub%projectors type
    proj_count = 1
    do ps=1,pub_cell%num_hub_species
       hub%projectors%species_num_proj(ps) = h_species(ps)%hub_ang_mom*2 + 1

       ! ndmh: find Gmax for this species
       do hat=1,pub_cell%nat_hub
          orig_iat = hub%orig(hat)
          if (hub%species_number(hat) == ps) exit
       end do
       rmax = elements(orig_iat)%radius
       n_r_rad_pts = int(100.0_DP*rmax/spacing)
       deltar = rmax / real(n_r_rad_pts - 1,kind=DP)
       n_g_rad_pts = nint(rmax/deltar) + 1

       ! ndmh: calculate projector in reciprocal space for this species
       call hubbard_internal_radial_proj(h_species(ps))
       hub%projectors%n_rad_pts(ps) = n_g_rad_pts
       hub%projectors%gmax(ps) = gmax
       hub%projectors%num_shells(ps) = 1
       hub%projectors%ang_mom(1,ps) = h_species(ps)%hub_ang_mom
       hub%projectors%species_first_proj(ps) = proj_count
       proj_count = proj_count + hub%projectors%species_num_proj(ps)
    end do

    ! ndmh: copy species number, atom centre and NGWF radius for Hubbard atoms
    hub%projectors%proj_species(:) = -1
    hub%projectors%proj_centre(:)%X = 0.0_DP
    hub%projectors%proj_centre(:)%Y = 0.0_DP
    hub%projectors%proj_centre(:)%Z = 0.0_DP
    hub%projectors%proj_max_radius(:) = -999.999_DP
    do hat=1,pub_cell%nat_hub
       orig_iat = hub%orig(hat)
       iat = pub_distr_atom(orig_iat)
       hub%projectors%proj_species(iat) = hub%species_number(hat)
       hub%projectors%proj_max_radius(iat) = elements(orig_iat)%radius
       hub%projectors%proj_centre(iat) = elements(orig_iat)%centre
    end do

    call timer_clock('hubbard_species_proj',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving hubbard_species_proj'
#endif

  contains


    subroutine hubbard_internal_radial_proj(hd_species)

      !==================================================================!
      ! This module builds the radial part of hydrogenic wavefunctions   !
      ! in reciprocal space.                                             !
      !------------------------------------------------------------------!
      ! Written by David D. O'Regan in April 2009.                       !
      !==================================================================!

      use comms, only: pub_on_root
      use hubbard_init, only: HUBBARD_SPECIES

      implicit none

      ! Arguments
      type(HUBBARD_SPECIES), intent(in) :: hd_species

      ! Local variables
      integer :: g_rad_pt                       ! Grid point counters
      real(kind=DP) :: gnorm                    ! Norm in recip space
      real(kind=DP) :: delg                     ! Grid variables
      real(kind=DP) :: gradius, denom, gath
      real(kind=DP) :: proj_max

      ! ddor: Set parameters for radial grid
      delg = PI / rmax
      proj_max = 0.0_DP
      gnorm = 0.0_DP
      gmax = 0.0_DP

      gath = ( hd_species%hub_charge )**( 0.5_DP*( &
           &real(2*hd_species%hub_ang_mom + 5,kind=DP) ) )

      ! ddor: Construct radial component in real space
      do g_rad_pt=1,n_g_rad_pts

         gradius =  real(g_rad_pt-1,DP) * delg

         denom = ( (hd_species%hub_charge)**2.0_DP + &
              (gradius * real(hd_species%hub_ang_mom + 1))**2.0_DP )**&
              (-1.0_DP*real(hd_species%hub_ang_mom) - 2.0_DP)

         select case (hd_species%hub_ang_mom)

         case(0) ! 1s

            hub%projectors%rad_proj_recip(g_rad_pt,1,ps) = &
                 SQRT(32.0_DP / PI) * gath * denom

         case(1) ! 2p

            hub%projectors%rad_proj_recip(g_rad_pt,1,ps) = &
                 (128.0_DP / SQRT(3.0_DP*PI) ) * gradius * gath * denom

         case(2) ! 3d

            hub%projectors%rad_proj_recip(g_rad_pt,1,ps) = &
                 1728.0_DP * SQRT( 3.0_DP / (5.0_DP*PI) ) * &
                 &(gradius**2.0_DP) * gath * denom

         case(3) ! 4f

            hub%projectors%rad_proj_recip(g_rad_pt,1,ps) = &
                 131072.0_DP * SQRT( 2.0_DP / (35.0_DP*PI) ) * &
                 &(gradius**3.0_DP) * gath * denom

         case default
            if (pub_on_root) write(stdout,*) 'Radial or angular quantum',&
                 ' number too high in', &
                 ' hubbard_internal_get_radial_proj. Only up to r=4,l=3',&
                 ' atomic radial functions are currently supported.'
            call comms_abort
         end select

         gnorm = gnorm + delg*(gradius * &
              hub%projectors%rad_proj_recip(g_rad_pt,1,ps))**2.0_DP

         ! ddor: Set core radius if magnitude of
         !       projector is less then 0.01% of
         !       maximum value and cumulative norm is at least 99.99%.
         proj_max = max(proj_max,abs(hub%projectors%rad_proj_recip(g_rad_pt,1,ps)))
         if (g_rad_pt > 1) then
            if ( (abs(hub%projectors%rad_proj_recip(g_rad_pt,1,ps)/proj_max) < &
                 &1.0e-4_DP) .and. (gnorm > 0.9999_DP) )then
               gmax = gradius ! The maximum necessary G
               n_g_rad_pts = g_rad_pt
               exit
            end if
         end if

      end do

!CW
 if(g_rad_pt>size(hub%projectors%rad_proj_recip,1))then
  write(*,*) 'Rad proj recip : out-bound argument'
  g_rad_pt=size(hub%projectors%rad_proj_recip,1)
 endif
!END

      ! ddor: Sanity check on Fourier-Bessel transform
      if  ( (abs(hub%projectors%rad_proj_recip(g_rad_pt,1,ps)/proj_max) > &
           &1.0e-4_DP) .or. (gnorm < 0.9999_DP) ) then
         gmax = real(n_g_rad_pts-1,DP) * delg
         if (pub_on_root) write(stdout,*) 'Fourier-Bessel normalisation',&
              ' not reached in', &
              ' hubbard_internal_get_radial_proj. Normalisation',&
              ' attained was.',gnorm*100,'%',' gmax set to',gmax
      end if

    end subroutine hubbard_internal_radial_proj

  end subroutine hubbard_species_proj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_species_exit_proj(hub)

    !==========================================================!
    ! This subroutine deallocates memory allocated for         !
    ! data on projecors.                                       !
    !----------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/1/2004.           !
    ! Tidied up by Peter D. Haynes on 1/7/2004.                !
    ! Modified for DFT+U module by David O'Regan in July 2008  !
    ! Modified by Nicholas Hine in November 2011 to use the    !
    ! HUBBARD_MODEL type                                       !
    !==========================================================!

    use projectors, only: projectors_deallocate_set

    implicit none

    ! Arguments
    type(HUBBARD_MODEL), intent(inout) :: hub

    call projectors_deallocate_set(hub%projectors)

  end subroutine hubbard_species_exit_proj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine hubbard_calculate_forces(hub_forces, &
       ngwfs_on_grid,ngwf_basis,hub_proj_basis,hub,pur_denskern, &
       hub_overlap,hub_overlap_t)

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the Hubbard DFT+U correction.                                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !    hub_forces      : output : DFT+U forces                              !
    !    ngwfs_on_grid   : input  : all of the tight boxes                    !
    !    ngwf_basis      : input  : all of the ngwf spheres                   !
    !    hub_proj_basis  : input  : all of the projector spheres              !
    !    pur_denskern    : input  : purified density kernel SPAM3             !
    !-------------------------------------------------------------------------!
    ! Written by David D. O'Regan in September 2009 based on                  !
    ! pseudo_nl_calculate_forces of the ODG.                                  !
    ! Modified by Nicholas Hine in November 2011 to use the                   !
    ! HUBBARD_MODEL type                                                      !
    !=========================================================================!

    use comms, only: comms_barrier, comms_reduce, pub_on_root, pub_my_node_id
    use function_basis, only: FUNC_BASIS
    use integrals, only : integrals_grad
    use ion, only: ELEMENT
    use simulation_cell, only: pub_cell
    use parallel_strategy, only: pub_first_atom_on_node, pub_orig_atom
    use projectors, only: projectors_func_grad_ovlp_box
    use rundat, only: task, pub_hubbard_restart, pub_hubbard_atomsolve
    use sparse, only: SPAM3, sparse_axpy, sparse_create, sparse_destroy, &
         sparse_get_element, sparse_product, sparse_transpose, &
         sparse_scale, sparse_rms_element
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: hub_forces(1:3,pub_cell%nat)
    type(FUNC_BASIS) :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(FUNC_BASIS) :: hub_proj_basis
    type(HUBBARD_MODEL) :: hub
    type(SPAM3), intent(in) :: pur_denskern(pub_cell%num_spins)
    type(SPAM3), intent(in) :: hub_overlap
    type(SPAM3), intent(in) :: hub_overlap_t

    ! Local Variables
    type(SPAM3) :: siGp_overlap(3),iGps_overlap(3),iGps_overlap_buffer
    type(SPAM3) :: sp_overlap_hub,hub_ps_overlap
    type(SPAM3) :: hub_force_mat(3)
    type(SPAM3) :: kh,rkh
    integer :: is
    integer :: cart
    integer :: iat, orig_iat
    integer :: atom_proj, global_proj
    real(kind=DP) :: hub_force

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering hubbard_calculate_forces'
#endif

    ! Start timer
    call timer_clock('hubbard_calculate_forces',1)

    ! Initialise
    hub_forces  = 0.0_DP

    ! Create result matrices to hold DFT+U forces
    do cart=1,3
       hub_force_mat%structure = 'G'
       call sparse_create(hub_force_mat(cart))
       call sparse_scale(hub_force_mat(cart),0.0_DP)
    end do

    ! ddor: G = U/2 [ 1 - 2 (Occupancy matrix in projector representation) ]
    !       Call it hubbard_term

    ! Create matrices to hold <NGWF_a|Proj_i> G and G <Proj_i|NGWF_a>
    call sparse_create(hub_ps_overlap,hub_overlap_t)
    call sparse_create(sp_overlap_hub,hub_overlap)

    ! Create matrices to hold <phi|iG*proj> overlap matrices in each direction
    do cart=1,3
       call sparse_create(siGp_overlap(cart),hub_overlap)
    end do

    ! Calculate <phi|iG*proj_> overlap matrix
    if (((task == 'HUBBARDSCF') .and. (hub%consistency_iteration > 1)) &
         & .or. pub_hubbard_restart .or. pub_hubbard_atomsolve) then
       call integrals_grad(siGp_overlap,ngwfs_on_grid,ngwf_basis,&
            hub%consistency_projs,hub_proj_basis)
    else
       do cart=1,3
          call projectors_func_grad_ovlp_box(siGp_overlap(cart), &
               ngwfs_on_grid,ngwf_basis,hub_proj_basis,hub%projectors,cart)
       end do
    end if

    ! Transpose them to get <iG*proj_|phi> overlap matrix
    ! This is the gradient of projectors, not duals of projectors.
    do cart=1,3
       call sparse_create(iGps_overlap(cart),hub_overlap_t)
       ! ddor: If using NGWF projectors, tensorially correct projector gradient
       if (((task == 'HUBBARDSCF') .and. (hub%consistency_iteration > 1)) &
            & .or. pub_hubbard_restart .or. pub_hubbard_atomsolve) then
          call sparse_create(iGps_overlap_buffer,hub_overlap_t)
          call sparse_transpose(iGps_overlap_buffer,siGp_overlap(cart))
          ! The gradient commutes with the tensorial correction so that
          ! <iG*proj^|phi> = O^^ <iG*proj_|phi>
          call sparse_product(iGps_overlap(cart),&
               hub%o_matrix,iGps_overlap_buffer)
          call sparse_destroy(iGps_overlap_buffer)
          ! ddor: Otherwise the transpose will do
       else
          call sparse_transpose(iGps_overlap(cart),siGp_overlap(cart))
       end if
    end do

    spinloop: do is= 1, pub_cell%num_spins

       ! Calculate the matrices <NGWF_a|Proj_i> H^_ and H^_ <Proj^i|NGWF_a>
       call sparse_product(sp_overlap_hub,hub_overlap,hub%projector_ham(is)) !V
       call sparse_product(hub_ps_overlap,hub%projector_ham(is),hub_overlap_t) !W

       ! Create temporary matrices kh and rkh
       call sparse_create(kh,pur_denskern(is),sp_overlap_hub)
       call sparse_create(rkh,hub_force_mat(1))

       ! Loop over Cartesian co-ordinates
       do cart=1,3
          ! Calculate H^_ <proj^j|phi_b> K^ab <phi_a|iG.proj_i>
          call sparse_product(kh,pur_denskern(is),siGp_overlap(cart))
          call sparse_product(rkh,hub_ps_overlap,kh)
          call sparse_axpy(hub_force_mat(cart),rkh,1.0_DP)

          ! Calculate <iG.proj^j|phi_b> K^ab <phi_a|proj_i> H^_
          call sparse_product(kh,pur_denskern(is),sp_overlap_hub)
          call sparse_product(rkh,iGps_overlap(cart),kh)
          call sparse_axpy(hub_force_mat(cart),rkh,1.0_DP)

          ! Force mat has projector tensor format F^m_m'
       end do

       ! Destroy temporary matrices
       call sparse_destroy(rkh)
       call sparse_destroy(kh)

    end do spinloop

    ! Destroy temporary matrices
    do cart=3,1,-1
       call sparse_destroy(iGps_overlap(cart))
       call sparse_destroy(siGp_overlap(cart))

    end do
    call sparse_destroy(sp_overlap_hub)
    call sparse_destroy(hub_ps_overlap)

    ! ddor: Force depends only on centre of Hubbard projector, that
    !       is the position of the Hubbard atom. There is no contribution
    !       to the force on atoms overlapping with Hubbard projectors.

    ! Loop over atoms
    do iat=pub_first_atom_on_node(pub_my_node_id), &
         pub_first_atom_on_node(pub_my_node_id + 1) - 1

       ! Find atom number in input file order
       orig_iat = pub_orig_atom(iat)

       ! Loop over Hubbard projectors on this atom
       do atom_proj=1,hub_proj_basis%num_on_atom(iat)
          global_proj = hub_proj_basis%first_on_atom(iat) + atom_proj - 1

          ! Loop over Cartesian co-ordinates
          do cart=1,3

             ! Find contribution of this projector to force on this atom
             ! from diagonal elements of hub_force_mat for this coordinate
             call sparse_get_element(hub_force,hub_force_mat(cart), &
                  global_proj, global_proj)
             hub_forces(cart,orig_iat) = hub_forces(cart,orig_iat) + hub_force

          end do  ! cart

       end do  ! atom_proj

    end do   ! loc_iat

    ! ddor: Put spin degeneracy into DFT+U force
    if (pub_cell%num_spins .eq. 2) then
       hub_forces = hub_forces * 2.0_DP
    end if

    ! Reduce result across nodes
    call comms_barrier
    call comms_reduce('SUM',hub_forces)

    ! Destroy temporary matrices
    do cart=1,3
       call sparse_destroy(hub_force_mat(cart))
    end do

    ! Stop timer
    call timer_clock('hubbard_calculate_forces',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving hubbard_calculate_forces'
#endif

  end subroutine hubbard_calculate_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CW
#include "hubbard_build_mod_dmft.h"
#include "hubbard_build_mod_gen_routines.h"
!END CW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module hubbard_build

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$  subroutine hubbard_ngwf_ovlp(overlap,hub_overlap,ngwf_basis,hub_proj_basis)
!!$
!!$    !==========================================================!
!!$    ! This calculates the updated  <NGWF|projector> matrix,    !
!!$    ! using as projectors the NGWFs on Hubbard atoms with the  !
!!$    ! greatest projection onto a set of hydrogenic orbitals    !
!!$    ! which are chosen by the user.                            !
!!$    ! It must be preceded by a call to integrals_brappd_ketppd !
!!$    ! to calculate the required sparse projection matrix.      !
!!$    !----------------------------------------------------------!
!!$    ! Written by David O'Regan in April 2009.                  !
!!$    !==========================================================!
!!$
!!$    use comms, only: pub_my_node_id, comms_abort, &
!!$         pub_on_root, comms_barrier
!!$    use function_basis, only: FUNC_BASIS
!!$    use hubbard_init, only: pub_cell%nat_hub
!!$    use parallel_strategy, only: pub_node_of_atom, &
!!$         pub_distr_atom, pub_num_hub_atoms_on_node, pub_hub_atoms_on_node
!!$    use sparse, only: SPAM3, sparse_get_element, sparse_put_element, &
!!$         sparse_scale, sparse_element_exists
!!$    use timer, only: timer_clock
!!$
!!$    implicit none
!!$
!!$    ! ddor: Arguments
!!$    type(SPAM3), intent(in) :: overlap
!!$    type(FUNC_BASIS), intent(in) :: ngwf_basis
!!$    type(FUNC_BASIS), intent(in) :: hub_proj_basis
!!$    type(SPAM3), intent(inout) :: hub_overlap
!!$
!!$    ! ddor: local variables
!!$    integer :: hub_atom, hat, theatom
!!$    integer :: hub_proj_we_are_on
!!$    integer :: row_ngwf
!!$    integer :: my_first_ngwf, my_last_ngwf
!!$    integer :: ll, nn
!!$    real(kind=DP) :: overlap_element
!!$    logical :: overlap_test
!!$
!!$#ifdef DEBUG
!!$    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
!!$         &hubbard_ngwf_ovlp'
!!$#endif
!!$
!!$    call timer_clock('hubbard_ngwf_ovlp',1)
!!$
!!$    ! ddor: Initially set output to zero (just the part on my node)
!!$    call sparse_scale(hub_overlap,0.0_DP)
!!$
!!$    my_first_ngwf = ngwf_basis%first_on_node(pub_my_node_id)
!!$    my_last_ngwf = ngwf_basis%first_on_node(pub_my_node_id) + &
!!$         ngwf_basis%num_on_node(pub_my_node_id) - 1
!!$
!!$    ! ddor: Loop over Hubbard atoms on my node
!!$    !       This works since we only consider Hubbard projectors
!!$    !       which are NGWFs on the Hubbard atom, thus on the same node.
!!$    hubatoms: do hub_atom = 1, pub_num_hub_atoms_on_node(pub_my_node_id)
!!$
!!$       hat = pub_hub_atoms_on_node(hub_atom,pub_my_node_id)
!!$       theatom = pub_distr_atom(hub%orig(hat)) !Atom number of Hubbard atom
!!$
!!$       ! ddor: Loop over Hubbard projectors on this atom
!!$       do nn = 1, hub_proj_basis%num_on_atom(theatom)
!!$
!!$          ! ddor: find the location of the corresponding projector
!!$          hub_proj_we_are_on = hub_proj_basis%first_on_atom(theatom) + nn - 1
!!$
!!$          ! ddor: find the location of NGWFs producing the greatest sum of
!!$          !       square overlaps with the atomic projectors
!!$          !       in terms of new atom numbers
!!$          ll = hub%localisedngwfs(nn)
!!$
!!$          !write(stdout, *)'Assigning NGWF number ',ll,&
!!$          !     ' which has original number ',orig_col_atom,&
!!$          !     ' as projector ',nn,' on this Hubbard atom'
!!$
!!$          ! ddor: Forcing only projectors on the Hubbard atom to be allowed
!!$          if (( theatom .eq. ngwf_basis%atom_of_func(ll) ) .and. &
!!$               & ( theatom .eq. &
!!$               & hub_proj_basis%atom_of_func(hub_proj_we_are_on) )) then
!!$
!!$             !ddor: Loop over NGWFs on my node
!!$             row_ngwf_loop: do row_ngwf = my_first_ngwf,my_last_ngwf
!!$
!!$                !ddor: Find out if they overlap with the Hubbard projector
!!$                overlap_test = sparse_element_exists(&
!!$                     &hub_overlap,row_ngwf,hub_proj_we_are_on)
!!$
!!$                ! Only proceed if NGWF overlaps with current projector
!!$                if (overlap_test) then
!!$
!!$                   ! ddor: Get the corresponding element of the projection
!!$                   ! ddor: matrix and put it in hub_overlap
!!$                   call sparse_get_element(overlap_element,overlap,row_ngwf,ll)
!!$
!!$                   ! ddor: Set the convention for the internal Hubbard matrices
!!$                   call sparse_put_element(overlap_element,&
!!$                        hub_overlap,row_ngwf,hub_proj_we_are_on)
!!$
!!$                end if
!!$
!!$             end do row_ngwf_loop
!!$
!!$          else ! ddor: The projector lies on a different atom
!!$             write(stdout, *) &
!!$                  &'WARNING: This Hubbard projector belongs to&
!!$                  & another atom though by construction&
!!$                  & it should not, aborting.'
!!$             call comms_abort
!!$
!!$          end if
!!$
!!$       end do
!!$
!!$    end do hubatoms ! ddor: End loop over Hubbard atoms
!!$
!!$    !ddor: Collect all elements across nodes
!!$    call comms_barrier
!!$
!!$    call timer_clock('hubbard_ngwf_ovlp',2)
!!$
!!$#ifdef DEBUG
!!$    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
!!$         &hubbard_ngwf_ovlp'
!!$#endif
!!$
!!$  end subroutine hubbard_ngwf_ovlp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$
!!$  subroutine hubbard_ngwf_grad_test(pur_denskern, hub_projector_ham, &
!!$       ngwf_basis, hub_proj_basis, hub_overlap, hub_overlap_t, &
!!$       hub_o_matrix, change_in_ngwfs, step)
!!$
!!$    !==========================================================!
!!$    ! This calculates the first order change in DFT+U energy   !
!!$    ! when the NGWFs are altered.                              !
!!$    !----------------------------------------------------------!
!!$    ! Written by David O'Regan in October 2009                 !
!!$    !==========================================================!
!!$
!!$    use comms, only: pub_on_root
!!$    use function_basis, only: FUNC_BASIS
!!$    use ion, only: ELEMENT
!!$    use integrals, only : integrals_brappd_ketppd
!!$    use projectors, only: projectors_func_ovlp_box
!!$    use simulation_cell, only: pub_cell
!!$    use sparse, only: sparse_create, sparse_destroy, &
!!$         SPAM3, sparse_product, sparse_trace, &
!!$         sparse_transpose, sparse_copy, sparse_scale
!!$    use rundat, only : task
!!$
!!$    implicit none
!!$
!!$    ! ddor: Arguments:
!!$
!!$    ! ddor: The density matrix
!!$    type(SPAM3), intent(in) :: pur_denskern(pub_cell%num_spins)
!!$    type(SPAM3), intent(in) :: hub_projector_ham(pub_cell%num_spins)
!!$    !type(ELEMENT), intent(in)       :: elements(pub_cell%nat)
!!$    type(FUNC_BASIS), intent(in)    :: ngwf_basis
!!$    type(FUNC_BASIS), intent(in) :: hub_proj_basis
!!$    real(kind=DP), intent(in)       :: change_in_ngwfs(&
!!$         &ngwf_basis%n_ppds*pub_cell%n_pts)
!!$    real(kind=DP), intent(in)       :: step
!!$    type(SPAM3), intent(in) :: hub_overlap
!!$    type(SPAM3), intent(in) :: hub_overlap_t
!!$    type(SPAM3), intent(in) :: hub_o_matrix
!!$
!!$    ! ddor: Local buffers
!!$    type(SPAM3) :: kernel_buf
!!$    type(SPAM3) :: hub_delta_matrix
!!$    type(SPAM3) :: hub_delta_tmatrix, hub_delta_tmatrix_buf
!!$    type(SPAM3) :: v_hamiltonian_buffer
!!$    type(SPAM3) :: v_hamiltonian_buffer_w
!!$    real(kind=DP) :: delta_e
!!$    integer :: is
!!$
!!$#ifdef DEBUG
!!$    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
!!$         &hubbard_ngwf_grad_test'
!!$#endif
!!$
!!$    call sparse_create(kernel_buf,pur_denskern(1))
!!$    call sparse_create(hub_delta_matrix,hub_overlap)
!!$    call sparse_create(hub_delta_tmatrix_buf,hub_overlap_t)
!!$    call sparse_create(hub_delta_tmatrix,hub_overlap_t)
!!$    call sparse_create(v_hamiltonian_buffer,hub_overlap)                  !VG=V
!!$    call sparse_create(v_hamiltonian_buffer_w,&
!!$         hub_overlap,hub_overlap_t)                                      !VGW=VW
!!$
!!$    if (task == 'HUBBARDSCF') then
!!$
!!$       ! < psi_n | f_b >
!!$       call integrals_brappd_ketppd(hub_delta_tmatrix_buf, &      !input-output
!!$            hub%consistency_projs, hub_proj_basis, change_in_ngwfs, & !input
!!$            ngwf_basis)                                               !input
!!$
!!$       ! O^mn < psi_n | f_b > =  < psi^m | f_b >
!!$       call sparse_product(hub_delta_tmatrix,hub_o_matrix,hub_delta_tmatrix_buf)
!!$
!!$       call sparse_transpose(hub_delta_matrix,hub_delta_tmatrix_buf)
!!$
!!$    else
!!$
!!$       call projectors_func_ovlp_box(hub_delta_matrix, &
!!$            change_in_ngwfs,ngwf_basis,hub_proj_basis,hub%projectors)
!!$
!!$       ! < psi_n | f_b >
!!$       call sparse_transpose(hub_delta_tmatrix,hub_delta_matrix)
!!$
!!$    end if
!!$
!!$    delta_e = 0.0_DP
!!$
!!$    spin: do is = 1, pub_cell%num_spins
!!$
!!$       call sparse_copy(kernel_buf,pur_denskern(is))
!!$
!!$       if (pub_cell%num_spins .eq. 1) then
!!$          call sparse_scale(kernel_buf,0.5_DP)
!!$       end if
!!$
!!$       ! ddor: Calculate Hubbard Hamiltonian U/2 (P - 2 PKP) + aP
!!$       !       This subroutine expects a kernel normalised without spin-
!!$       !       degeneracy factor, and it produces a Hamiltonian for same.
!!$       !call hubbard_projector_ham(hamiltonian_buffer,kernel_buf,hub_overlap, &
!!$       !     hub_overlap_t,is)
!!$
!!$       ! <f_a | \phi_m> H^m_n
!!$       call sparse_product(v_hamiltonian_buffer,&
!!$            hub_delta_matrix,hub_projector_ham(is))
!!$
!!$       ! <f_a | \phi_m> H^m_n < psi^n | \phi_b >
!!$       call sparse_product(v_hamiltonian_buffer_w,&
!!$            v_hamiltonian_buffer,hub_overlap_t)
!!$
!!$       ! Tr [ <f_a | \phi_m> H^m_n < psi^n | \phi_b > K^bc ]
!!$       delta_e = delta_e + sparse_trace(v_hamiltonian_buffer_w,kernel_buf)
!!$
!!$       ! <\phi_a | \phi_m> H^m_n
!!$       call sparse_product(v_hamiltonian_buffer,&
!!$            hub_overlap,hub_projector_ham(is))
!!$
!!$       ! <\phi_a | \phi_m> H^m_n < psi^n | f_b >
!!$       call sparse_product(v_hamiltonian_buffer_w,&
!!$            v_hamiltonian_buffer,hub_delta_tmatrix)
!!$       !v_hamiltonian_buffer,hub_delta_tmatrix_buf)
!!$
!!$       ! Tr [ <\phi_a | \phi_m> H^m_n < psi^n | f_b > K^bc ]
!!$       delta_e = delta_e + sparse_trace(v_hamiltonian_buffer_w,kernel_buf)
!!$
!!$    end do spin
!!$
!!$    delta_e = delta_e * step
!!$
!!$    if (pub_cell%num_spins .eq. 1) then
!!$       delta_e = delta_e * 2.0_DP
!!$    end if
!!$
!!$    if (pub_on_root) write(*,*) 'Expected change in DFT+U energy, step length', &
!!$         delta_e, step
!!$
!!$    call sparse_destroy(v_hamiltonian_buffer_w)
!!$    call sparse_destroy(v_hamiltonian_buffer)
!!$    call sparse_destroy(hub_delta_tmatrix)
!!$    call sparse_destroy(hub_delta_tmatrix_buf)
!!$    call sparse_destroy(hub_delta_matrix)
!!$    call sparse_destroy(kernel_buf)
!!$
!!$#ifdef DEBUG
!!$    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
!!$         &hubbard_ngwf_grad_test'
!!$#endif
!!$
!!$    !fftbox_batch(:,:,:,:,3,:) = 0.0_DP
!!$    !fftbox_batch(:,:,:,:,4,:) = 0.0_DP
!!$
!!$    !call hubbard_ngwf_grad_test(pur_denskern, elements, &
!!$    !     ngwf_basis, hub_proj_basis, direction_on_grid, trial_length)
!!$    !! This trick will always miss the minus sign
!!$    !write(*,*) 'alpha G dot G', trial_length * &
!!$    !     wrappers_ddot(ngwf_basis%n_ppds*pub_cell%n_pts, &
!!$    !     contra_grad_on_grid,1,cov_grad_on_grid,1)
!!$    !   call electronic_energy_components(&
!!$    !        ewald_energy, overlap, pur_denskern,kinet,nonlocpot, &
!!$    !        localpseudo_fine,core_density_fine, &
!!$    !        ngwfs_on_grid,ngwf_basis,hub_proj_basis, &
!!$    !        elements,iteration, hfexchange)
!!$
!!$  end subroutine hubbard_ngwf_grad_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$  subroutine hubbard_kernel_grad_test(hub_projector_ham, &
!!$       hub_overlap, hub_overlap_t, change_in_kernel, step)
!!$
!!$    !==========================================================!
!!$    ! This calculates the first order change in DFT+U energy   !
!!$    ! when the kernel is altered.                              !
!!$    !----------------------------------------------------------!
!!$    ! Written by David O'Regan in October 2009                 !
!!$    !==========================================================!
!!$
!!$    use comms, only: pub_on_root
!!$    use simulation_cell, only: pub_cell
!!$    use sparse, only: sparse_create, sparse_destroy, &
!!$         SPAM3, sparse_product, sparse_trace
!!$    use timer, only: timer_clock
!!$
!!$    implicit none
!!$
!!$    ! ddor: Arguments:
!!$
!!$    ! ddor: The density matrix
!!$    type(SPAM3), intent(in) :: change_in_kernel(pub_cell%num_spins)
!!$    real(kind=DP), intent(in)       :: step
!!$    type(SPAM3), intent(in) :: hub_projector_ham(pub_cell%num_spins)
!!$    type(SPAM3), intent(in) :: hub_overlap
!!$    type(SPAM3), intent(in) :: hub_overlap_t
!!$
!!$    ! ddor: Local buffers
!!$    type(SPAM3) :: v_hamiltonian_buffer
!!$    type(SPAM3) :: hamiltonian
!!$    real(kind=DP) :: delta_e
!!$    integer :: is
!!$
!!$#ifdef DEBUG
!!$    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
!!$         &hubbard_kernel_grad_test'
!!$#endif
!!$
!!$    call sparse_create(v_hamiltonian_buffer,hub_overlap)                !VG=V
!!$    hamiltonian%structure = 'H'
!!$    call sparse_create(hamiltonian)
!!$
!!$    delta_e = 0.0_DP
!!$
!!$    spin: do is = 1, pub_cell%num_spins
!!$
!!$       ! ddor: Calculate Hubbard Hamiltonian U/2 (P - 2 PKP) + aP
!!$       !       This subroutine expects a kernel normalised without spin-
!!$       !       degeneracy factor, and it produces a Hamiltonian for same.
!!$       !call hubbard_projector_ham(hamiltonian_buffer,pur_denskern(is), &
!!$       !     hub_overlap,hub_overlap_t,is)
!!$
!!$       ! ddor: V [ U/2 (1-2N) + a]  blanked to 'VG'
!!$       call sparse_product(v_hamiltonian_buffer,&
!!$            hub_overlap,hub_projector_ham(is))
!!$       ! ddor: Hamiltonian
!!$       call sparse_product(hamiltonian,&
!!$            v_hamiltonian_buffer,hub_overlap_t)
!!$
!!$       delta_e = delta_e + sparse_trace(hamiltonian,change_in_kernel(is))
!!$
!!$    end do spin
!!$
!!$    delta_e = delta_e * step
!!$
!!$    if (pub_cell%num_spins .eq. 1) then
!!$       delta_e = delta_e * 2.0_DP
!!$    end if
!!$
!!$    if (pub_on_root) write(*,*) 'Expected change in DFT+U energy, step length', &
!!$         delta_e, step
!!$
!!$    call sparse_destroy(hamiltonian)
!!$    call sparse_destroy(v_hamiltonian_buffer)
!!$
!!$#ifdef DEBUG
!!$    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
!!$         &hubbard_kernel_grad_test'
!!$#endif
!!$
!!$    line_search_step = 0.0001_DP
!!$    call hubbard_kernel_grad_test(aux, con_direction, line_search_step)
!!$    expected = 0.0_DP
!!$    do is = 1, pub_cell%num_spins
!!$       expected = expected + sparse_trace(con_gradient(is),co_gradient(is))
!!$    end do
!!$    ! This trick will always miss the minus sign
!!$    write(*,*) 'alpha K dot K', line_search_step * expected * real(pub_cell%num_spins,kind=DP)
!!$
!!$    ! Update the auxiliary matrix
!!$    !if (pub_lnv_check_trial_steps) then
!!$    !    call internal_kernel_protected_step(trial_aux,line_search_step, &
!!$    !         aux,con_direction)
!!$    !    do is=1,pub_cell%num_spins
!!$    !       call sparse_copy(aux(is),trial_aux(is))
!!$    !    end do
!!$    !else
!!$    do is=1,pub_cell%num_spins
!!$       !      call sparse_axpy(aux(is),con_direction(is),line_search_step)
!!$       call sparse_axpy(aux(is),con_gradient(is),line_search_step)
!!$    end do
!!$    !end if
!!$
!!$    internal_energy = hubbard_energy
!!$    !call sparse_axpy(ham(is),hubbard_ham(is),1.0_DP)
!!$    call sparse_copy(ham(is),hubbard_ham(is))
!!$
!!$  end subroutine hubbard_kernel_grad_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

