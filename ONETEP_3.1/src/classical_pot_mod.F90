! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                   Classical potential module                   !
!                                                                !
! This module implements classical potential approaches.         !
!----------------------------------------------------------------!
! Written by Chris-Kriton Skylaris and Alvaro Ruiz Serrano       !
! on 21/04/2009.                                                 !
!================================================================!
module classical_pot

  use ion, only: element

  implicit none

  private

  public :: classical_pot_ii_energy
  public :: classical_pot_init
  public :: classical_pot_struct_fac
  public :: classical_pot_dealloc
  public :: classical_pot_recip
  public :: classical_pot_ii_forces

  type(ELEMENT), save, public, allocatable, dimension(:) :: classical_elements

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine classical_pot_init(nat_classical)

    use utils, only: utils_alloc_check
    implicit none

    integer, intent(in) :: nat_classical
    integer :: ierr

    allocate(classical_elements(nat_classical),stat=ierr)
    call utils_alloc_check('classical_pot_init', &
         'classical_elements',ierr)

  end subroutine classical_pot_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine classical_pot_struct_fac(struct_fac_classical,grid)

    !===========================================================!
    ! This subroutine generates a "structure factor" for a      !
    ! collection of point charges ("classical atoms").          !
    !-----------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris and Alvaro Ruiz Serrano  !
    ! on 21/04/2009 starting from the                           !
    ! pseudo_make_structure_factor subroutine.                  !
    !===========================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_node_id, comms_abort, pub_on_root, comms_barrier
    use constants, only: DP, stdout
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    complex(kind=DP), intent(out) :: struct_fac_classical(grid%ld3,grid%ld2, &
         grid%max_slabs23)

    ! Local variables
    integer :: i3,i2,islab23
    integer :: atom
    real(kind=DP) :: gvec(3),gdotr

    ! Start timer
    call timer_clock('classical_pot_struct_fac',1)

    if (pub_on_root) write(stdout,'(/a,i6,a)',advance='no') &
         'Calculating structure factor for ', pub_cell%nat_classical, &
         ' classical atoms ...'

    struct_fac_classical = 0.0_DP

    ! Loop over reciprocal space on this node
    do islab23=1,grid%num_slabs23      ! along b1
       do i2=1,grid%n2                 ! along b2
          ! ndmh: initialise this slab of the structure factor array
          struct_fac_classical(:,i2,islab23) = (0.0_DP,0.0_DP)
          do i3=1,grid%n3              ! along b3

             call cell_grid_recip_pt(gvec,islab23 + &
                  grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)

             ! Loop over atoms in cell
             do atom=1,pub_cell%nat_classical
                gdotr = gvec(1) * classical_elements(atom)%centre%x + &
                     gvec(2) * classical_elements(atom)%centre%y + &
                     gvec(3) * classical_elements(atom)%centre%z

                struct_fac_classical(i3,i2,islab23) = &
                     struct_fac_classical(i3,i2,islab23) + &
                     classical_elements(atom)%ion_charge * &
                     exp(cmplx(0.0_DP,-gdotr,kind=DP))

             end do   ! loop over atoms

          end do      ! loop along b3
       end do         ! loop along b2
    end do            ! loop along b1


    ! Stop timer
    call timer_clock('classical_pot_struct_fac',2)

    call comms_barrier
    if (pub_on_root) write(stdout,'(a)') '... done'

    ! kaw: removed this as deallocation is now done in the forces routine,
    ! in hindsight the deallocation should be done both here and there and
    ! then executed depending on the calculation type...
    ! Perhaps call the routine after the call to this one rather than inside?
    ! Are the calculation type variable available here? Would we even want
    ! them to be?

    ! cks: classical elements are no longer needed
    ! call classical_pot_dealloc

  end subroutine classical_pot_struct_fac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine classical_pot_dealloc

    use utils, only: utils_dealloc_check
    implicit none

    integer :: ierr

    deallocate(classical_elements,stat=ierr)
    call utils_dealloc_check('classical_pot_dealloc', &
         'classical_elements',ierr)

  end subroutine classical_pot_dealloc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine classical_pot_recip(fine_complex, & ! output
         struct_fac_classical, grid) ! input

    !===============================================================!
    ! This subroutine generates in reciprocal space the total local !
    ! Coulomb potential due to a collection of point charges.       !
    !-------------------------------------------------------------- !
    ! Written by Chris-Kriton Skylaris and Alvaro Ruiz Serrano      !
    ! on 21/04/2009 starting from the                               !
    ! pseudopotentials_sum_local_rec subroutine.                    !
    ! Modified by Nicholas Hine in December 2009 to allow less than !
    ! one slab23 per node.                                          !
    ! Moved to classical_pot_mod by Nicholas Hine in July 2010.     !
    ! Gaussian smearing of point charges added by Chris-Kriton      !
    ! Skylaris on 24/9/2010.                                        !
    !===============================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_node_id, pub_total_num_nodes, comms_abort, &
         pub_on_root, comms_barrier
    use constants, only: DP, PI, stdout, VERBOSE
    use rundat, only: pub_coulomb_cutoff, pub_output_detail
    use services, only: services_flush
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in)   :: grid
    complex(kind=DP), intent(in) :: struct_fac_classical(grid%ld3,grid%ld2, &
         grid%max_slabs23)
    complex(kind=DP), intent(inout) :: fine_complex(grid%ld3,grid%ld2,&
         grid%max_slabs23)

    ! Local variables
    integer :: i3,i2,islab23        ! Reciprocal grid loop counters
    real(kind=DP) :: v_loc_value
    real(kind=DP) :: factor
    real(kind=DP) :: r_expo
    real(kind=DP) :: g_expo
    real(kind=DP) :: gvec(3), gsq



    factor = 4.0_DP * PI / grid%weight
    ! cks: spread the point charge to a Gaussian with halfwidth of 0.3 a0
    ! cks: the real-space exponent is ln(2)/(0.3^2)
    r_expo =log(2.0_DP)/0.09_DP
    ! kaw: Change to stupid number to allow comparison between QMMM and QM.
    !r_expo = 50000000.0

    g_expo =1.0_DP/(4.0_DP* r_expo)

    ! cks: report the embedding
    if (pub_on_root .and. pub_output_detail == VERBOSE) &
       write(stdout,'(a,f12.10,a)',advance='no') &
       'Calculating embedding potential with Gauss exp= ',r_expo,'  ...'


    ! Loop over reciprocal space grid on this node
    do islab23=1,grid%num_slabs23           ! along b1
       do i2=1,grid%n2                      ! along b2
          do i3=1,grid%n3                   ! along b3

             !cks: g^2
             call cell_grid_recip_pt(gvec,islab23 + &
                  grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)
             gsq = sum(gvec(1:3)**2)

             ! cks: v_loc_value
             ! cks: The g=0 term is zero for PBC while it is equal to
             ! cks: (4Pi/V)*(Rc^2/2) for cutoff Coulomb
             v_loc_value = -grid%coulomb_recip(i3,i2,islab23)*factor &
                  *exp(-g_expo*gsq)

             ! cks: multiply with "structure factor" to obtain Coulomb potential
             ! cks: for all point charges
             fine_complex(i3,i2,islab23) = fine_complex(i3,i2,islab23) + &
                  struct_fac_classical(i3,i2,islab23) * v_loc_value

          end do   ! b3
       end do      ! b2
    end do         ! b1


    ! G=0 element must be real
    if (pub_my_node_id==grid%node_slab23(1)) then

       ! cks: add also the G=0 term, which in the case of PBC
       ! cks: should be equal to -Pi*Q/(V*r_expo).
       ! cks: We actually add the negative of this because of the convention
       ! cks: in electronic structure codes that the potential of nuclei is negative
       if (.not. pub_coulomb_cutoff) then
          fine_complex(1,1,1) = fine_complex(1,1,1)  &
               +PI*struct_fac_classical(1,1,1)/(r_expo*grid%weight)
       endif


       if (aimag(fine_complex(1,1,1)) /= 0.0_DP) then
          write(stdout,'(a)') 'Error in classical_pot_recip: &
               &potential not real'
          call comms_abort
       end if

    end if

    ! Nyquist filter (fine grid is always going to be even)
    fine_complex(grid%n3/2+1,:,:) = (0.0_DP,0.0_DP)
    fine_complex(:,grid%n2/2+1,:) = (0.0_DP,0.0_DP)
    ! Nyquist filter for last slab
    if (pub_my_node_id==grid%node_slab23(grid%n1/2+1)) &
         fine_complex(:,:,grid%num_slabs23) = (0.0_DP,0.0_DP)

    call comms_barrier
    if (pub_on_root) write(stdout,'(a)') '... done'
    call services_flush


  end subroutine classical_pot_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine classical_pot_ii_energy(ii_energy, & ! output
       elements)

    !=====================================================================!
    ! Calculates the direct ion-ion contribution to the energy. Replaces  !
    ! ewald_calculate_energy for finite systems                           !
    !---------------------------------------------------------------------!
    ! Arguments:                                                          !
    ! elements (input): list of all the atoms in the system               !
    ! ii_energy (output): the ion-ion energy                              !
    !---------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 10/01/2011, starting from       !
    ! cutoff_coulomb_ii_energy which was written by Nick Hine.            !
    !=====================================================================!

    use comms, only: pub_total_num_nodes, pub_my_node_id, comms_reduce
    use constants, only: DP
    use ion, only: element
    use geometry, only: magnitude, OPERATOR(-), POINT
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    real(kind =DP), intent(out) :: ii_energy
    type(ELEMENT), intent(in) :: elements(:)


    type(ELEMENT), allocatable, dimension(:) :: all_elements
    integer :: all_atoms
    integer :: iatom
    integer :: jatom
    integer :: ierr
    integer :: iatom_start
    integer :: iatom_end
    integer :: nc_atoms
    real(kind =DP) :: qi, qj
    real(kind =DP) :: r_dist
    type(POINT) :: rij

    all_atoms =pub_cell%nat +pub_cell%nat_classical

    all_atoms =pub_cell%nat +pub_cell%nat_classical

    allocate(all_elements(pub_cell%nat +pub_cell%nat_classical),stat=ierr)
    call utils_alloc_check('classical_pot_ii_energy', &
         'all_elements',ierr)

    ! cks: stick ions and classical atoms in one big array
    all_elements(1: pub_cell%nat) = elements(1: pub_cell%nat)
    all_elements(pub_cell%nat+1 : pub_cell%nat +pub_cell%nat_classical)= &
         classical_elements(1: pub_cell%nat_classical)


    ! cks: parallelisation of outer loop
    nc_atoms    =all_atoms/pub_total_num_nodes
    iatom_start =1 +pub_my_node_id*nc_atoms
    iatom_end   =(pub_my_node_id+1)*nc_atoms
    if (pub_my_node_id == (pub_total_num_nodes -1) ) then
       iatom_end =iatom_end +mod(all_atoms, pub_total_num_nodes)
    endif

    ii_energy =0.0_DP
    do iatom = iatom_start, iatom_end

       qi = all_elements(iatom)%ion_charge

       ! Loop over all other atoms
       do jatom = 1, all_atoms

          if (jatom == iatom ) cycle

          ! Find charge and distance of other atom
          qj = all_elements(jatom)%ion_charge
          rij = all_elements(iatom)%centre - all_elements(jatom)%centre

          r_dist =magnitude(rij)

          ! cks: avoid division by zero in the case of overlapping real atoms
          ! cks: and embedding (classical) atoms
          if (r_dist > tiny(1.0_DP)) then
             ! Add point charge contribution to energy
             ii_energy = ii_energy + qi*qj/magnitude(rij)
          endif


       enddo
    enddo

    ! cks: sum effort of each core
    call comms_reduce('SUM',ii_energy)

    ! remove double counting
    ii_energy = ii_energy * 0.5_DP

!    if (pub_on_root) then
       write(*,*)
       write(*,*)    'Ion-Ion Energy        : ',ii_energy
!    end if
    deallocate(all_elements,stat=ierr)
    call utils_dealloc_check('classical_pot_ii_energy', &
         'all_elements',ierr)


  end subroutine classical_pot_ii_energy



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine classical_pot_ii_forces(ii_forces,elements)

    !=====================================================================!
    ! Calculates direct QM ion-classical ion contribution to the forces   !
    ! acting on the QM ions.                                              !
    !---------------------------------------------------------------------!
    ! Arguments:                                                          !
    ! elements (input): list of all the atoms in the system               !
    ! ii_forces (output): the contribution to the forces acting on the QM !
    ! ions that arises from the classical ions.                           !
    !---------------------------------------------------------------------!
    ! Written by Karl Wilkinson on 12/12/2011, starting from              !
    ! classical_pot_ii_force which was written by Chris-Kriton Skylaris   !
    ! (based upon cutoff_coulomb_ii_energy which was written by Nick Hine)!
    ! and cutoff_coulomb_ii_forces written by  Nicholas Hine.             !
    !=====================================================================!

    use comms, only: pub_total_num_nodes, pub_my_node_id, comms_reduce
    use constants, only: DP
    use ion, only: element
    use geometry, only: magnitude, OPERATOR(-), POINT
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    real (kind=DP), dimension(1:3,1:pub_cell%nat), intent(out) :: ii_forces
    type(ELEMENT), intent(in) :: elements(:)


    type(ELEMENT), allocatable, dimension(:) :: all_elements
    integer :: iatom
    integer :: jatom
    integer :: ierr
    integer :: iatom_start
    integer :: iatom_end
    integer :: nc_atoms
    integer :: all_atoms
    real(kind=DP) :: qi, qj
    real(kind=DP) :: r_dist
    type(POINT) :: rij
    real(kind=DP) :: rij3
    real(kind=DP) :: fij(3)        ! force between particles i and j

    ! kaw: Place ions and classical atoms in one big array
    all_atoms =pub_cell%nat +pub_cell%nat_classical
    allocate(all_elements(pub_cell%nat +pub_cell%nat_classical),stat=ierr)
    call utils_alloc_check('classical_pot_ii_forces', &
         'all_elements',ierr)
    all_elements(1: pub_cell%nat) = elements(1:pub_cell%nat)
    all_elements(pub_cell%nat+1:all_atoms)=classical_elements(1:pub_cell%nat_classical)


    ! kaw: parallelisation of outer loop
    nc_atoms    =pub_cell%nat/pub_total_num_nodes
    iatom_start =1 + pub_my_node_id*nc_atoms
    iatom_end   =(pub_my_node_id+1)*nc_atoms
    if (pub_my_node_id == (pub_total_num_nodes -1) ) then
       iatom_end =iatom_end +mod(pub_cell%nat, pub_total_num_nodes)
    endif

    ii_forces(:,:) = 0.0d0
    ! kaw: loop over all QM atoms for this node
    do iatom = iatom_start, iatom_end
       qi = all_elements(iatom)%ion_charge

       ! kaw: Loop over all classical atoms
       do jatom = pub_cell%nat + 1, pub_cell%nat_classical + pub_cell%nat

          ! Find charge and distance of other atom
          qj = all_elements(jatom)%ion_charge
          rij = all_elements(iatom)%centre - all_elements(jatom)%centre
          rij3 =magnitude(rij)**3

          ! Find force between atoms i and j
          fij(1) = (qi*qj/rij3)*rij%X
          fij(2) = (qi*qj/rij3)*rij%Y
          fij(3) = (qi*qj/rij3)*rij%Z

          ! Add contribution of jatom to force on iatom
          ii_forces(:,iatom) = ii_forces(:,iatom) + fij

       enddo
    enddo

    ! cks: sum effort of each core
    call comms_reduce('SUM',ii_forces)

    deallocate(all_elements,stat=ierr)
    call utils_dealloc_check('classical_pot_ii_forces', &
         'all_elements',ierr)

    call classical_pot_dealloc

  end subroutine classical_pot_ii_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module classical_pot
