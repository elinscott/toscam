! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Nicholas D.M. Hine
!
!   Thomas Young Centre
!   Imperial College London
!   Exhibition Road
!
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module augmentation

  use constants, only: DP, PI, stdout

  implicit none

  private

  ! Max number of augmentation functions on any atom
  integer :: max_proj_tot

  integer, parameter :: pad_pts = 4

  public :: aug_projector_denskern
  public :: aug_nonlocal_mat
  public :: aug_nonlocal_commutator_mat
  public :: augmentation_overlap
  public :: augmentation_pos
  public :: augmentation_grad
  public :: augmentation_density_on_grid
  public :: augmentation_box_init
  public :: augmentation_density_forces
  public :: augmentation_screen_dij
  public :: aug_nl_calculate_forces

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_box_init(grid)

    !==================================================================!
    ! This subroutine calculates the size of the augmentation box used !
    ! in the calculation of nhat and other quantities.                 !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  grid (in) : the whole-cell grid for which the augmentation      !
    !              box is required.                                    !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30/06/10.                            !
    ! Moved to augmentation_mod by Nicholas Hine on 14/02/11.          !
    !==================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce, pub_on_root
    use geometry, only: geometry_magnitude
    use paw, only: paw_projectors
    use pseudopotentials, only: nlps_projectors
    use rundat, only: pub_paw, pub_usp
    use simulation_cell, only: pub_aug_box_n1, pub_aug_box_n2, pub_aug_box_n3
    use utils, only: utils_abort

    ! Arguments
    type(GRID_INFO), intent(in) :: grid

    ! Local variables
    real(kind=DP) :: max_rad

    ! Find size of augmentation density box
    if (pub_paw) then
       max_rad = maxval(paw_projectors%proj_max_radius(:))
       max_proj_tot = maxval(paw_projectors%species_num_proj(:))
    else if (pub_usp) then
       max_rad = maxval(nlps_projectors%proj_max_radius(:))
       max_proj_tot = maxval(nlps_projectors%species_num_proj(:))
    else
       ! why are we here?
       call utils_abort('Error in augmentation_box_init: Augmentation box &
            &requested, but no augmentation present')
       max_rad = 0.0_DP ! suppresses warning
    end if

    ! ndmh: find minimum sizes, check they are smaller than cell, and make odd
    pub_aug_box_n1 = 2.0_DP*max_rad / geometry_magnitude(grid%da1) + pad_pts
    pub_aug_box_n1 = min(grid%n1,pub_aug_box_n1)
    pub_aug_box_n1 = pub_aug_box_n1 - modulo(pub_aug_box_n1,2) + 1
    call comms_reduce('MAX',pub_aug_box_n1)
    pub_aug_box_n2 = 2.0_DP*max_rad / geometry_magnitude(grid%da2) + pad_pts
    pub_aug_box_n2 = min(grid%n2,pub_aug_box_n2)
    pub_aug_box_n2 = pub_aug_box_n2 - modulo(pub_aug_box_n2,2) + 1
    call comms_reduce('MAX',pub_aug_box_n2)
    pub_aug_box_n3 = 2.0_DP*max_rad / geometry_magnitude(grid%da3) + pad_pts
    pub_aug_box_n3 = min(grid%n3,pub_aug_box_n3)
    pub_aug_box_n3 = pub_aug_box_n3 - modulo(pub_aug_box_n3,2) + 1
    call comms_reduce('MAX',pub_aug_box_n3)

    ! Report size of aug box
    if (pub_on_root) then
       write(stdout,'(a)') ''
       write(stdout,'(a)') '============================= Charge Augmentation =============================='
       write(stdout,'(3(a,i4))') '                           Aug box size: ',&
            pub_aug_box_n1,' x',pub_aug_box_n2,' x',pub_aug_box_n3
       write(stdout,'(a)') '================================================================================'
    end if

  end subroutine augmentation_box_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine aug_projector_denskern(proj_denskern,denskern,sp_overlap)

    !====================================================================!
    ! This subroutine creates the reduced density matrix for each        !
    ! atomic site, in a diagonal matrix of size nproj x nproj. This      !
    ! is rho_ij in the PAW and USP formalisms, where \rho_{ij} is the    !
    ! occupancy of partial wave i,j.                                     !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  proj_denskern (inout) : The projector density kernel rho^ij       !
    !  denskern (in) : The NGWF density kernel K^ab                      !
    !  sp_overlap (in) : The NGWF-Projector overlap matrix <phi_a|p_i>.  !
    !--------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.             !
    ! Moved to projectors_mod by Nicholas Hine, 11/02/11.                !
    !====================================================================!

    use comms, only: pub_on_root
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product, &
         sparse_transpose, sparse_scale, sparse_transpose_structure
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: proj_denskern(pub_cell%num_spins)
    type(SPAM3), intent(in) :: sp_overlap
    type(SPAM3), intent(in) :: denskern(pub_cell%num_spins)

    ! Local Variables
    type(SPAM3) :: ps_overlap
    type(SPAM3) :: ksp
    integer :: is

    ! Transpose <phi_a|p^i> to get <p^j|phi_b>
    call sparse_transpose_structure(ps_overlap%structure,sp_overlap)
    call sparse_create(ps_overlap)
    call sparse_transpose(ps_overlap,sp_overlap)
    call sparse_create(ksp,denskern(1),sp_overlap)

    ! Calculate <p^i|\hat{\rho}|p^j>, the projector density kernel:
    ! rho^ij = <p^i|phi_a>.(K^ab.<phi_b|p^j>)
    do is=1,pub_cell%num_spins
       call sparse_product(ksp,denskern(is),sp_overlap)
       call sparse_product(proj_denskern(is),ps_overlap,ksp)
    end do

    call sparse_destroy(ksp)
    call sparse_destroy(ps_overlap)

  end subroutine aug_projector_denskern


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_overlap(overlap,bra_proj_overlap,ket_proj_overlap)

    !====================================================================!
    ! This subroutine augments an overlap matrix by adding the           !
    ! contribution from the overlap operator S = 1 + |p_i>O_ij<p_j|.     !
    ! The matrix returned is                                             !
    !      <bra_a|S|ket_b> = <bra_a|ket_b> + <bra_a|p_i>O_ij<p_j|ket_b>  !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  overlap (inout) : On entry: a general overlap matrix between two  !
    !      set of functions <bra_a|ket_b>.                               !
    !      On exit: the overlap matrix augmented with the projector part !
    !  bra_proj_overlap (in) : The Bra-Projector overlap matrix.         !
    !  ket_proj_overlap (in) : The Ket-Projector overlap matrix.         !
    !--------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 28/05/10.             !
    ! Moved to augmentation_mod by Nicholas Hine 14/02/11.               !
    ! Adapted to cope with augmenting overlaps between different sets of !
    ! functions in bras and kets.                                        !
    !====================================================================!

    use comms, only: pub_on_root
    use pseudopotentials, only: pseudo_aug_Q_matrix
    use paw, only: paw_projector_overlap
    use rundat, only: pub_paw, pub_usp
    use sparse, only: SPAM3, sparse_axpy, sparse_create, sparse_destroy, &
         sparse_product, sparse_put_element, sparse_transpose, &
         sparse_transpose_structure
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: overlap
    type(SPAM3), intent(in) :: bra_proj_overlap
    type(SPAM3), intent(in), optional :: ket_proj_overlap

    ! Local Variables
    type(SPAM3) :: proj_overlap
    type(SPAM3) :: proj_ket_overlap
    type(SPAM3) :: bra_proj_overlap_O

    ! Check inputs
    if (present(ket_proj_overlap)) then
       if (ket_proj_overlap%iscmplx.neqv.bra_proj_overlap%iscmplx) then
          call utils_abort('Error in augmentation_overlap: bra and ket overlap &
               &matrices must be both real or both complex')
       end if
    end if

    ! Create temporary matrices to store <bra_a|p_i>, <bra_a|p_i> O_ij  and O_ij
    if (present(ket_proj_overlap)) then
       call sparse_transpose_structure(proj_ket_overlap%structure, &
            ket_proj_overlap)
    else
       call sparse_transpose_structure(proj_ket_overlap%structure, &
            bra_proj_overlap)
    end if
    call sparse_create(proj_ket_overlap,iscmplx=bra_proj_overlap%iscmplx)
    call sparse_create(bra_proj_overlap_O,bra_proj_overlap)
    proj_overlap%structure = 'E'
    call sparse_create(proj_overlap,iscmplx=bra_proj_overlap%iscmplx)

    if (.not.present(ket_proj_overlap)) then
       ! If functions of the bras are the same as the functions of the kets, then just
       ! transpose <bra_a|p_i> to get <p_i|ket_a>
       call sparse_transpose(proj_ket_overlap,bra_proj_overlap)
    else
       ! If they are different, transpose the supplied <ket_a|p_i> to get <p_i|ket_a>
       call sparse_transpose(proj_ket_overlap,ket_proj_overlap)
    end if

    ! Get projector overlap matrix O_ij
    if (pub_paw) then
       call paw_projector_overlap(proj_overlap)
    else if (pub_usp) then
       call pseudo_aug_Q_matrix(proj_overlap)
    end if

    ! Create product <bra_a|p_i> O_ij
    call sparse_product(bra_proj_overlap_O,bra_proj_overlap,proj_overlap)

    ! Create NGWF-sized proj_overlap matrix <bra_a|p_i>O_ij<p_j|ket_b>
    call sparse_destroy(proj_overlap)
    call sparse_create(proj_overlap,overlap)
    call sparse_product(proj_overlap,bra_proj_overlap_O,proj_ket_overlap)

    ! Add it to the unaugmented overlap matrix <bra_a|ket_b>
    call sparse_axpy(overlap,proj_overlap,1.0_DP)

    ! Clean up temporary matrices
    call sparse_destroy(proj_overlap)
    call sparse_destroy(bra_proj_overlap_O)
    call sparse_destroy(proj_ket_overlap)

  end subroutine augmentation_overlap


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_pos(r_elements,proj_basis,bra_proj_overlap, &
       ket_proj_overlap,axis)

    !====================================================================!
    ! This subroutine augments the matrix for the position operator in   !
    ! the NGWF representation, by adding to the r_elements matrix the    !
    ! contribution |p_i>R_ij<p_j| from the sphere terms.                 !
    ! where in PAW, R_ij = <phi_i|r|phi_j>-<tphi_i|r|tphi_j>.            !
    ! The matrices returned (one for each cartesian direction) are       !
    !   <bra_a|r|ket_b> = <bra_a|r|ket_b> + <bra_a|p_i>R_ij<p_j|ket_b>   !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  r_elements(3) (inout) : On entry: position operator matrices      !
    !      between two set of functions <bra_a|r|ket_b>.                 !
    !      On exit: the pos matrix augmented with the projector part     !
    !  proj_basis (in) : FUNC_BASIS type describing projectors.          !
    !  bra_proj_overlap (in) : The Bra-Projector overlap matrix.         !
    !  ket_proj_overlap (in) : The Ket-Projector overlap matrix.         !
    !  axis (in,optional)    : axis specifier (all axes if not present)  !
    !--------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 02/12/11.             !
    !====================================================================!

    use comms, only: pub_on_root, pub_my_node_id
    use function_basis, only: FUNC_BASIS
    use pseudopotentials, only: pseudo_aug_Q_matrix
    use parallel_strategy, only: pub_first_atom_on_node, pub_elements_on_node
    use paw, only: paw_projector_overlap, paw_position_operator
    use rundat, only: pub_paw, pub_usp
    use sparse, only: SPAM3, sparse_axpy, sparse_create, sparse_destroy, &
         sparse_product, sparse_put_element, &
         sparse_transpose, sparse_transpose_structure, sparse_get_element
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: r_elements(3)
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(SPAM3), intent(in) :: bra_proj_overlap
    type(SPAM3), intent(in), optional :: ket_proj_overlap
    integer, intent(in), optional :: axis

    ! Local Variables
    type(SPAM3) :: proj_overlap
    type(SPAM3) :: r_sphere(3)
    type(SPAM3) :: proj_r_elements
    type(SPAM3) :: proj_ket_overlap
    type(SPAM3) :: bra_proj_overlap_R
    integer :: xyz, axmin, axmax
    integer :: loc_iproj, iproj, jproj
    integer :: iat, loc_iat
    real(kind=DP) :: R_atom(3), r_el, o_el

    ! Check inputs
    if (present(ket_proj_overlap)) then
       if (ket_proj_overlap%iscmplx.neqv.bra_proj_overlap%iscmplx) then
          call utils_abort('Error in augmentation_overlap: bra and ket overlap &
               &matrices must be both real or both complex')
       end if
    end if
    if (present(axis)) then
       axmin = axis
       axmax = axis
    else
       axmin = 1
       axmax = 3
    end if

    ! Create temporary matrices to store <bra_a|p_i>, <bra_a|p_i> R_ij  and R_ij
    if (present(ket_proj_overlap)) then
       call sparse_transpose_structure(proj_ket_overlap%structure, &
            ket_proj_overlap)
    else
       call sparse_transpose_structure(proj_ket_overlap%structure, &
            bra_proj_overlap)
    end if
    call sparse_create(proj_ket_overlap,iscmplx=bra_proj_overlap%iscmplx)
    call sparse_create(bra_proj_overlap_R,bra_proj_overlap)
    proj_overlap%structure = 'E'
    call sparse_create(proj_overlap,iscmplx=bra_proj_overlap%iscmplx)
    do xyz=axmin,axmax
       r_sphere(xyz)%structure = 'E'
       call sparse_create(r_sphere(xyz),iscmplx=bra_proj_overlap%iscmplx)
    end do

    if (.not.present(ket_proj_overlap)) then
       ! If functions of the bras are the same as the functions of the kets, then just
       ! transpose <bra_a|p_i> to get <p_i|ket_a>
       call sparse_transpose(proj_ket_overlap,bra_proj_overlap)
    else
       ! If they are different, transpose the supplied <ket_a|p_i> to get <p_i|ket_a>
       call sparse_transpose(proj_ket_overlap,ket_proj_overlap)
    end if

    ! Get projector position operator R_ij(1:3) and overlap matrix O_ij
    if (pub_paw) then
       call paw_projector_overlap(proj_overlap)
       ! ndmh: calculate sphere part of position operator
       call paw_position_operator(r_sphere)
    else if (pub_usp) then
       call pseudo_aug_Q_matrix(proj_overlap)
       call utils_abort('Error in augmentation_pos: USP position operator &
            &does not exist')
    end if

    ! ndmh: cycle over projectors on this node, applying correction to
    ! ndmh: r_sphere to move it to the atom centre
    do loc_iproj=1,proj_basis%num_on_node(pub_my_node_id)
       iproj = loc_iproj + proj_basis%first_on_node(pub_my_node_id) - 1
       iat = proj_basis%atom_of_func(iproj)
       loc_iat = iat - pub_first_atom_on_node(pub_my_node_id) + 1
       do jproj=proj_basis%first_on_atom(iat), &
            proj_basis%first_on_atom(iat)+proj_basis%num_on_atom(iat)-1

          ! Extract overlap element
          call sparse_get_element(o_el,proj_overlap,jproj,iproj)

          R_atom(1) = pub_elements_on_node(loc_iat)%centre%x
          R_atom(2) = pub_elements_on_node(loc_iat)%centre%y
          R_atom(3) = pub_elements_on_node(loc_iat)%centre%z

          ! Extract element from r_sphere and shift by R_atom*o_el
          ! ddor: get elements for one direction only if axis is specified
          do xyz=axmin,axmax
             call sparse_get_element(r_el,r_sphere(xyz),jproj,iproj)
             r_el = R_atom(xyz) * o_el + r_el
             call sparse_put_element(r_el,r_sphere(xyz),jproj,iproj)
          end do

       end do
    end do

    call sparse_create(proj_r_elements,r_elements(axmin))

    do xyz=axmin,axmax

       ! Create product <bra_a|p_i> R_ij
       call sparse_product(bra_proj_overlap_R,bra_proj_overlap,r_sphere(xyz))

       ! Create NGWF-sized r_elements matrix <bra_a|p_i>R_ij<p_j|ket_b>
       call sparse_product(proj_r_elements,bra_proj_overlap_R,proj_ket_overlap)

       ! Add it to the unaugmented overlap matrix <bra_a|ket_b>
       call sparse_axpy(r_elements(xyz),proj_r_elements,1.0_DP)

    end do

    ! Clean up temporary matrices
    call sparse_destroy(proj_r_elements)
    do xyz=axmax,axmin,-1
       call sparse_destroy(r_sphere(xyz))
    end do
    call sparse_destroy(proj_overlap)
    call sparse_destroy(bra_proj_overlap_R)
    call sparse_destroy(proj_ket_overlap)

  end subroutine augmentation_pos


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_grad(grad_elements,bra_proj_overlap, &
       ket_proj_overlap,axis)

    !====================================================================!
    ! This subroutine augments the matrix for the position operator in   !
    ! the NGWF representation, by adding to the r_elements matrix the    !
    ! contribution |p_i>R_ij<p_j| from the sphere terms.                 !
    ! where in PAW, R_ij = <phi_i|r|phi_j>-<tphi_i|r|tphi_j>.            !
    ! The matrices returned (one for each cartesian direction) are       !
    !   <bra_a|r|ket_b> = <bra_a|r|ket_b> + <bra_a|p_i>R_ij<p_j|ket_b>   !
    !--------------------------------------------------------------------!
    ! Arguments:                                                         !
    !  grad_elements(3) (inout) : On entry: grad operator matrices       !
    !      between two set of functions <bra_a|nabla|ket_b>.             !
    !      On exit: the grad matrix augmented with the projector part    !
    !  bra_proj_overlap (in) : The Bra-Projector overlap matrix.         !
    !  ket_proj_overlap (in) : The Ket-Projector overlap matrix.         !
    !  axis (in,optional)    : axis specifier (all axes if not present)  !
    !--------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 02/12/11.             !
    !====================================================================!

    use comms, only: pub_on_root, pub_my_node_id
    use function_basis, only: FUNC_BASIS
    use pseudopotentials, only: pseudo_aug_Q_matrix
    use parallel_strategy, only: pub_first_atom_on_node, pub_elements_on_node
    use paw, only: paw_projector_overlap, paw_grad_operator
    use rundat, only: pub_paw, pub_usp
    use sparse, only: SPAM3, sparse_axpy, sparse_create, sparse_destroy, &
         sparse_product, sparse_transpose, sparse_transpose_structure
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: grad_elements(3)
    type(SPAM3), intent(in) :: bra_proj_overlap
    type(SPAM3), intent(in), optional :: ket_proj_overlap
    integer, intent(in), optional :: axis

    ! Local Variables
    type(SPAM3) :: proj_overlap
    type(SPAM3) :: grad_sphere(3)
    type(SPAM3) :: proj_grad_elements
    type(SPAM3) :: proj_ket_overlap
    type(SPAM3) :: bra_proj_overlap_grad
    integer :: xyz, axmin, axmax

    ! Check inputs
    if (present(ket_proj_overlap)) then
       if (ket_proj_overlap%iscmplx.neqv.bra_proj_overlap%iscmplx) then
          call utils_abort('Error in augmentation_overlap: bra and ket overlap &
               &matrices must be both real or both complex')
       end if
    end if
    if (present(axis)) then
       axmin = axis
       axmax = axis
    else
       axmin = 1
       axmax = 3
    end if

    ! Create temporary matrices to store <bra_a|p_i>, <bra_a|p_i> grad_ij
    if (present(ket_proj_overlap)) then
       call sparse_transpose_structure(proj_ket_overlap%structure, &
            ket_proj_overlap)
    else
       call sparse_transpose_structure(proj_ket_overlap%structure, &
            bra_proj_overlap)
    end if
    call sparse_create(proj_ket_overlap,iscmplx=bra_proj_overlap%iscmplx)
    call sparse_create(bra_proj_overlap_grad,bra_proj_overlap)
    proj_overlap%structure = 'E'
    call sparse_create(proj_overlap,iscmplx=bra_proj_overlap%iscmplx)
    do xyz=axmin,axmax
       grad_sphere(xyz)%structure = 'E'
       call sparse_create(grad_sphere(xyz),iscmplx=bra_proj_overlap%iscmplx)
    end do

    if (.not.present(ket_proj_overlap)) then
       ! If functions of the bras are the same as the functions of the kets, 
       ! then just transpose <bra_a|p_i> to get <p_i|ket_a>
       call sparse_transpose(proj_ket_overlap,bra_proj_overlap)
    else
       ! If they are different, transpose the supplied <ket_a|p_i> 
       ! to get <p_i|ket_a>
       call sparse_transpose(proj_ket_overlap,ket_proj_overlap)
    end if

    ! Get projector position operator grad_ij(1:3) and overlap matrix O_ij
    if (pub_paw) then
       call paw_projector_overlap(proj_overlap)
       ! ndmh: calculate sphere part of grad operator
       call paw_grad_operator(grad_sphere)
    else if (pub_usp) then
       call pseudo_aug_Q_matrix(proj_overlap)
       call utils_abort('Error in augmentation_pos: USP grad operator &
            &does not exist')
    end if

    call sparse_create(proj_grad_elements,grad_elements(axmin))

    do xyz=axmin,axmax

       ! Create product <bra_a|p_i> grad_ij
       call sparse_product(bra_proj_overlap_grad,bra_proj_overlap, &
            grad_sphere(xyz))

       ! Create NGWF-sized r_elements matrix <bra_a|p_i>grad_ij<p_j|ket_b>
       call sparse_product(proj_grad_elements,bra_proj_overlap_grad, &
            proj_ket_overlap)

       ! Add it to the unaugmented overlap matrix <bra_a|ket_b>
       call sparse_axpy(grad_elements(xyz),proj_grad_elements,1.0_DP)

    end do

    ! Clean up temporary matrices
    call sparse_destroy(proj_grad_elements)
    do xyz=axmax,axmin,-1
       call sparse_destroy(grad_sphere(xyz))
    end do
    call sparse_destroy(proj_overlap)
    call sparse_destroy(bra_proj_overlap_grad)
    call sparse_destroy(proj_ket_overlap)

  end subroutine augmentation_grad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine aug_nonlocal_mat(nonlocpot,dijhat,rho_ij,sp_overlap, &
       paw_sphere_energies,show_matrices)

    !==================================================================!
    ! This subroutine creates the PAW nonlocal matrix given by         !
    !  V^nl_ab = <phi_a|p_i> D_ij <p_j|phi_b>                          !
    ! where D_ij are the nonlocal energies given by                    !
    !  D_ij = \hat{D}_ij + D^1_ij - \tilde{D}^1_ij                     !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  nonlocpot (inout) : Nonlocal potential matrix V^nl_ab           !
    !  dijhat (in) : Augmentation functions screened with locpot       !
    !  rho_ij (in) : Projector density kernel rho_ij                   !
    !  sp_overlap (in) : NGWF-Projector overlap matrix <phi_a|p_i>     !
    !  paw_sphere_energies (inout) : Sphere energy contribution terms  !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.           !
    !==================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_my_node_id, pub_on_root
    use constants, only: VERBOSE, paw_en_size, paw_en_dijhat
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: pub_num_atoms_on_node
    use paw, only: paw_nonlocal_energies
    use pseudopotentials, only: pseudo_get_dij
    use rundat, only: pub_paw, pub_usp, pub_paw_output_detail
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product, &
         sparse_trace, sparse_axpy, sparse_transpose, sparse_convert, &
         sparse_transpose_structure, sparse_show_matrix
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: nonlocpot(pub_cell%num_spins)
    type(SPAM3), intent(in) :: dijhat(pub_cell%num_spins)
    type(SPAM3), intent(in) :: rho_ij(pub_cell%num_spins)
    type(SPAM3), intent(in) :: sp_overlap
    real(kind=DP), intent(inout), optional :: paw_sphere_energies(paw_en_size)
    logical, intent(in), optional :: show_matrices

    ! Local Variables
    type(SPAM3), allocatable :: dij(:)
    type(SPAM3) :: sp_overlap_dij
    type(SPAM3) :: ps_overlap
    integer :: ierr
    integer :: is
    logical :: iscmplx
    logical :: loc_show_matrices

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering aug_nonlocal_mat'
#endif

    ! Start Timer
    call timer_clock('aug_nonlocal_mat',1)

    ! Optional argument
    loc_show_matrices = .false.
    if (present(show_matrices)) loc_show_matrices = show_matrices

    ! Allocate arrays of matrices
    allocate(dij(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('aug_nonlocal_mat','dij',ierr)
    do is=1,pub_cell%num_spins
       call sparse_create(dij(is),rho_ij(is))
    end do

    ! Calculate the nonlocal energies D_ij, and the sphere energies if
    ! required
    if (pub_paw) then
       if (present(paw_sphere_energies)) then
          call paw_nonlocal_energies(dij,rho_ij, &
               paw_sphere_energies,loc_show_matrices)
       else
          call paw_nonlocal_energies(dij,rho_ij, &
               show_matrices=loc_show_matrices)
       end if

       do is=1,pub_cell%num_spins
          if ((loc_show_matrices).and.(pub_paw_output_detail>=VERBOSE)) then
             if (pub_on_root) write(stdout,'(a,i4)') 'dijhat', is
             call sparse_show_matrix(dijhat(is),show_elems=.true.)
          end if
          call sparse_axpy(dij(is),dijhat(is),1.0_DP)
          if (present(paw_sphere_energies)) then
             paw_sphere_energies(paw_en_dijhat) = &
                  paw_sphere_energies(paw_en_dijhat) &
                  + sparse_trace(rho_ij(is),dijhat(is))
          end if
       end do

       if ((loc_show_matrices).and.(pub_paw_output_detail>=VERBOSE)) then
          do is=1,pub_cell%num_spins
             if (pub_on_root) write(stdout,'(a,i4)') 'dij total', is
             call sparse_show_matrix(dij(is),show_elems=.true.)
          end do
       end if

    else if (pub_usp) then

       do is=1,pub_cell%num_spins
          call pseudo_get_dij(dij(is))
       end do

    end if

    ! Create the matrix structures (set as real or complex depending on
    ! whether sp_overlap is real or complex) for sp and ps matrices
    iscmplx = sp_overlap%iscmplx
    call sparse_create(sp_overlap_dij,sp_overlap,iscmplx=iscmplx)
    call sparse_transpose_structure(ps_overlap%structure,sp_overlap)
    call sparse_create(ps_overlap,iscmplx=iscmplx)

    ! Create the transpose of the NGWF-projector overlap matrix
    call sparse_transpose(ps_overlap,sp_overlap)

    ! Loop over spins
    do is=1,pub_cell%num_spins

       ! Calculate the matrix <NGWF_a|Proj_i> * D_ij
       call sparse_product(sp_overlap_dij,sp_overlap,dij(is))

       ! Calculate the matrix \sum_i (<NGWF_a|Proj_i>D_ij<Proj_j|NGWF_b>)
       call sparse_product(nonlocpot(is),sp_overlap_dij,ps_overlap)

    end do

    ! Clean up temporary matrices
    call sparse_destroy(ps_overlap)
    call sparse_destroy(sp_overlap_dij)
    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(dij(is))
    end do
    deallocate(dij,stat=ierr)
    call utils_dealloc_check('aug_nonlocal_mat','dij',ierr)

    ! Start Timer
    call timer_clock('aug_nonlocal_mat',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving aug_nonlocal_mat'
#endif

  end subroutine aug_nonlocal_mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine aug_nonlocal_commutator_mat(nonlocpot_com, proj_basis, &
       dijhat, rho_ij,  ngwf_basis, ngwfs_on_grid, sp_overlap, delta_in)

    !==================================================================!
    ! This subroutine calculates the commutator between the nonlocal   !
    ! potential and the position operator for the 3 Cartesian          !
    ! directions.                                                      !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  nonlocpot (inout) : Nonlocal potential, position op commutator  !
    !  dijhat (in) : Augmentation functions screened with locpot       !
    !  rho_ij (in) : Projector density kernel rho_ij                   !
    !  sp_overlap (in) : NGWF-Projector overlap matrix <phi_a|p_i>     !
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 02/12/11.           !
    !==================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_my_node_id, pub_on_root
    use constants, only: VERBOSE, paw_en_size, paw_en_dijhat
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_num_atoms_on_node
    use paw, only: paw_projectors, paw_nonlocal_energies, paw_species_init_proj
    use projectors, only: PROJECTOR_SET, projectors_commutator, &
         projectors_deallocate_set
    use pseudopotentials, only: pseudo_get_dij
    use rundat, only: pub_paw, pub_usp, pub_paw_output_detail
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_axpy
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: nonlocpot_com(pub_cell%num_spins)
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(SPAM3), intent(in) :: dijhat(pub_cell%num_spins)
    type(SPAM3), intent(in) :: rho_ij(pub_cell%num_spins)
    type(SPAM3), intent(in) :: sp_overlap
    real(kind=DP), intent(in) :: delta_in ! finite difference shift

    ! Local Variables
    type(SPAM3), allocatable :: dij(:)
    integer :: ierr
    integer :: is

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &aug_nonlocal_commutator_mat'
#endif

    ! Start Timer
    call timer_clock('aug_nonlocal_commutator_mat',1)

    ! Allocate arrays of matrices
    allocate(dij(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('aug_nonlocal_commutator_mat','dij',ierr)
    do is=1,pub_cell%num_spins
       call sparse_create(dij(is),rho_ij(is))
    end do

    ! Calculate the nonlocal energies D_ij
    if (pub_paw) then
       call paw_nonlocal_energies(dij,rho_ij)
       do is=1,pub_cell%num_spins
          call sparse_axpy(dij(is),dijhat(is),1.0_DP)
       end do
    else if (pub_usp) then
       do is=1,pub_cell%num_spins
          call pseudo_get_dij(dij(is))
       end do
    end if

    do is=1,pub_cell%num_spins,1
       call projectors_commutator(nonlocpot_com(is), proj_basis, &
            ngwf_basis, ngwfs_on_grid, sp_overlap, paw_projectors, &
            delta_in, dij(is))
    end do

    ! Deallocate nonlocal energies
    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(dij(is))
    end do
    deallocate(dij,stat=ierr)
    call utils_dealloc_check('aug_nonlocal_commutator_mat','dij',ierr)

    ! Start Timer
    call timer_clock('aug_nonlocal_commutator_mat',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &aug_nonlocal_commutator_mat'
#endif

  end subroutine aug_nonlocal_commutator_mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_density_on_grid(nhat_den_grad,grid,denskern, &
       sp_overlap)

    !==================================================================!
    ! This subroutine creates the compensation density \hat{n}(r) on   !
    ! the simulation cell fine grid, and also the gradient of the      !
    ! augmentation density in each cartesian direction.                !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  grid (in) : Grid definition for the fine grid.                  !
    !  nhat (inout) : The compensation density \hat{n}(r) on the grid. !
    !  denskern (in) : The NGWF density kernel K^ab.                   !
    !  sp_overlap (in) : The NGWF-Projector overlap matrix <phi_a|p_i>.!
    !------------------------------------------------------------------!
    ! This subroutine was written by Nicholas Hine 24/05/10.           !
    ! Moved to projectors_mod by Nicholas Hine on 11/02/11.            !
    !==================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_deposit_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: comms_reduce, pub_on_root, pub_my_node_id
    use constants, only: max_spins, VERBOSE
    use function_basis, only: FUNC_BASIS
    use geometry, only: geometry_magnitude
    use parallel_strategy, only: pub_elements_on_node, pub_max_atoms_on_node, &
         pub_first_atom_on_node, pub_num_atoms_on_node
    use paw, only: paw_atom_aug_den
    use pseudopotentials, only: pseudo_atom_aug_den
    use rundat, only: pub_paw, pub_usp, pub_output_detail, pub_aug_den_dim, &
         pub_aug_funcs_recip
    use simulation_cell, only: pub_cell, pub_aug_box_n1, pub_aug_box_n2, &
         pub_aug_box_n3
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_get_block
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check
    use xc, only: pub_xc_gradient_corrected

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: nhat_den_grad(grid%ld1, &
         grid%ld2, grid%max_slabs12, pub_cell%num_spins, 0:pub_aug_den_dim)
    type(SPAM3), intent(in) :: denskern(pub_cell%num_spins)
    type(SPAM3), intent(in) :: sp_overlap

    ! Local Variables
    type(SPAM3), allocatable :: rho_ij(:)
    character(20) :: fmt,tmp
    integer :: iat, loc_iat
    integer :: ierr
    integer :: is
    integer :: isp
    integer :: box_n1,box_n2,box_n3
    integer :: box_start1,box_start2,box_start3
    integer :: cart
    real(kind=DP) :: total_nhat(max_spins)
    real(kind=DP) :: total_nhat_targ(max_spins)
    real(kind=DP),allocatable :: rho_ij_block(:,:,:)
    real(kind=DP), allocatable :: atom_nhat(:,:,:,:)
    real(kind=DP), allocatable :: atom_grad_nhat(:,:,:,:,:)
    real(kind=DP), allocatable :: atom_aug_func(:,:,:)
    real(kind=DP), allocatable :: atom_aug_func_grad(:,:,:,:)
    real(kind=DP), allocatable :: buffer(:,:,:)
    complex(kind=DP), allocatable :: atom_nhat_recip(:,:,:,:)
    complex(kind=DP), allocatable :: atom_grad_nhat_recip(:,:,:,:,:)
    complex(kind=DP), allocatable :: atom_aug_func_recip(:,:,:)
    complex(kind=DP), allocatable :: atom_aug_func_grad_recip(:,:,:,:)
    logical :: i_have_box

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &augmentation_density_on_grid'
#endif

    ! Start Timer
    call timer_clock('augmentation_density_on_grid',1)

    ! Find size of box
    box_n1 = pub_aug_box_n1
    box_n2 = pub_aug_box_n2
    box_n3 = pub_aug_box_n3

    ! Create projector density kernel
    allocate(rho_ij(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('augmentation_density_on_grid','rho_ij',ierr)
    do is=1,pub_cell%num_spins
       rho_ij(is)%structure = 'E'
       call sparse_create(rho_ij(is))
    end do

    call aug_projector_denskern(rho_ij,denskern,sp_overlap)

    ! Allocate workspace
    allocate(rho_ij_block(max_proj_tot,max_proj_tot, &
         pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('augmentation_density_on_grid','rho_ij_block',ierr)
    allocate(atom_nhat(box_n1,box_n2,box_n3,pub_cell%num_spins), &
         stat=ierr)
    call utils_alloc_check('augmentation_density_on_grid','atom_nhat',ierr)
    allocate(atom_grad_nhat(box_n1,box_n2,box_n3,pub_cell%num_spins,3), &
         stat=ierr)
    call utils_alloc_check('augmentation_density_on_grid','atom_grad_nhat', &
         ierr)
    allocate(buffer(box_n1,box_n2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('augmentation_density_on_grid','buffer',ierr)

    if (pub_aug_funcs_recip) then
       allocate(atom_aug_func_recip(box_n1,box_n2,box_n3),stat=ierr)
       call utils_alloc_check('augmentation_density_on_grid', &
            'atom_aug_func_recip',ierr)
       allocate(atom_aug_func_grad_recip(box_n1,box_n2,box_n3,3),stat=ierr)
       call utils_alloc_check('augmentation_density_on_grid', &
            'atom_aug_func_grad_recip',ierr)
       allocate(atom_nhat_recip(box_n1,box_n2,box_n3,pub_cell%num_spins), &
            stat=ierr)
       call utils_alloc_check('augmentation_density_on_grid', &
            'atom_nhat_recip',ierr)
       allocate(atom_grad_nhat_recip(box_n1,box_n2,box_n3,pub_cell%num_spins, &
            3),stat=ierr)
       call utils_alloc_check('augmentation_density_on_grid', &
            'atom_grad_nhat_recip',ierr)
    else
       allocate(atom_aug_func(box_n1,box_n2,box_n3),stat=ierr)
       call utils_alloc_check('augmentation_density_on_grid', &
            'atom_aug_func',ierr)
       allocate(atom_aug_func_grad(box_n1,box_n2,box_n3,3),stat=ierr)
       call utils_alloc_check('augmentation_density_on_grid', &
            'atom_aug_func_grad',ierr)
    end if

    total_nhat = 0.0_DP
    total_nhat_targ = 0.0_DP

    !if (pub_on_root) write(stdout,'(a)')' iat lup mup ipt ipw  li  mi jpt &
    !     &jpw  lj  mj            nLij             Gij           rhoij         &
    !     &    qij             qLM'
    do loc_iat=1,pub_max_atoms_on_node

       if (loc_iat<=pub_num_atoms_on_node(pub_my_node_id)) then
          iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1
          isp = pub_elements_on_node(loc_iat)%pspecies_number

          ! Find where box for this atom is located in simulation cell
          call cell_grid_box_start_wrt_atom( &
               box_start1, box_start2, box_start3, &
               pub_elements_on_node(loc_iat)%centre, box_n1, box_n2, box_n3, &
               grid)

          ! Get block of \rho_ij for this atom
          do is=1,pub_cell%num_spins
             call sparse_get_block(rho_ij_block(:,:,is),rho_ij(is),iat,iat)
          end do

          ! Call appropriate routine to generate the augmentation density for
          ! this atom in the augmentation box using rhoij
          if (pub_paw) then
             if (pub_aug_funcs_recip) then
                call paw_atom_aug_den(atom_nhat,atom_grad_nhat,total_nhat, &
                     total_nhat_targ,rho_ij_block, &
                     isp,pub_elements_on_node(loc_iat)%centre,grid, &
                     box_n1,box_n2,box_n3,box_start1,box_start2,box_start3, &
                     atom_aug_func_recip=atom_aug_func_recip, &
                     atom_aug_func_grad_recip=atom_aug_func_grad_recip, &
                     atom_nhat_recip=atom_nhat_recip, &
                     atom_grad_nhat_recip=atom_grad_nhat_recip)
             else
                call paw_atom_aug_den(atom_nhat,atom_grad_nhat,total_nhat, &
                     total_nhat_targ,rho_ij_block, &
                     isp,pub_elements_on_node(loc_iat)%centre,grid, &
                     box_n1,box_n2,box_n3,box_start1,box_start2,box_start3, &
                     atom_aug_func_real=atom_aug_func, &
                     atom_aug_func_grad_real=atom_aug_func_grad)
             end if
          else if (pub_usp) then
             call pseudo_atom_aug_den(atom_nhat,atom_grad_nhat,atom_aug_func, &
                  atom_aug_func_grad,total_nhat,total_nhat_targ,rho_ij_block, &
                  isp,pub_elements_on_node(loc_iat)%centre,grid, &
                  box_n1,box_n2,box_n3,box_start1,box_start2,box_start3)
          end if

          i_have_box = .true.
       else
          ! Nothing to deposit on this node
          i_have_box = .false.
       end if

       ! Deposit this box to the simulation cell if present, or just wait for
       ! data from other nodes if no box
       do is=1,pub_cell%num_spins
          call cell_grid_deposit_box(nhat_den_grad(:,:,:,is,0), &
               atom_nhat(:,:,:,is), buffer, grid, &
               box_n1, box_n2, box_n3, box_n1, box_n2, &
               box_start1, box_start2, box_start3, i_have_box, .false.)
          ! Deposit gradient to whole-cell arrays
          if (pub_xc_gradient_corrected) then
             do cart=1,3
                call cell_grid_deposit_box(nhat_den_grad(:,:,:,is,cart), &
                     atom_grad_nhat(:,:,:,is,cart), buffer, grid, &
                     box_n1, box_n2, box_n3, box_n1, box_n2, &
                     box_start1, box_start2, box_start3, i_have_box, .false.)
             end do
          end if
       end do
    end do

    ! Get sum of compensation density over all nodes, compare to target
    call comms_reduce('SUM',total_nhat(:))
    call comms_reduce('SUM',total_nhat_targ(:))
    if (any(abs(total_nhat_targ(1:pub_cell%num_spins) &
         - total_nhat(1:pub_cell%num_spins)) > 1.0d-15)) then
       if (pub_on_root.and.(pub_output_detail>=VERBOSE).and. &
            (.not.pub_aug_funcs_recip)) then
          write(tmp,'(i5)') pub_cell%num_spins
          write(fmt,'(3a)') '(a,',trim(adjustl(tmp)),'f14.8)'
          write(stdout,fmt) 'Total Compensation Charge Target     &
               &    : ',total_nhat_targ(1:pub_cell%num_spins)
          write(stdout,fmt) 'Total Compensation Charge on Regular &
               &Grid: ',total_nhat(1:pub_cell%num_spins)
       end if
    end if

    ! Deallocate temporary arrays and matrices
    if (pub_aug_funcs_recip) then
       deallocate(atom_grad_nhat_recip,stat=ierr)
       call utils_dealloc_check('augmentation_density_on_grid', &
            'atom_grad_nhat_recip',ierr)
       deallocate(atom_nhat_recip,stat=ierr)
       call utils_dealloc_check('augmentation_density_on_grid', &
            'atom_nhat_recip',ierr)
       deallocate(atom_aug_func_grad_recip,stat=ierr)
       call utils_dealloc_check('augmentation_density_on_grid', &
            'atom_aug_func_grad_recip',ierr)
       deallocate(atom_aug_func_recip,stat=ierr)
       call utils_dealloc_check('augmentation_density_on_grid', &
            'atom_aug_func_recip',ierr)
    else
       deallocate(atom_aug_func_grad,stat=ierr)
       call utils_dealloc_check('augmentation_density_on_grid', &
            'atom_aug_func_grad',ierr)
       deallocate(atom_aug_func,stat=ierr)
       call utils_dealloc_check('augmentation_density_on_grid', &
            'atom_aug_func',ierr)
    end if

    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('augmentation_density_on_grid','buffer',ierr)
    deallocate(atom_grad_nhat,stat=ierr)
    call utils_dealloc_check('augmentation_density_on_grid','atom_grad_nhat', &
         ierr)
    deallocate(atom_nhat,stat=ierr)
    call utils_dealloc_check('augmentation_density_on_grid','atom_nhat',ierr)
    deallocate(rho_ij_block,stat=ierr)
    call utils_dealloc_check('augmentation_density_on_grid','rho_ij_block',ierr)

    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(rho_ij(is))
    end do

    deallocate(rho_ij,stat=ierr)
    call utils_dealloc_check('augmentation_density_on_grid','rho_ij',ierr)

    ! Stop Timer
    call timer_clock('augmentation_density_on_grid',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving augmentation_density_on_grid'
#endif

  end subroutine augmentation_density_on_grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_screen_dij(dij,locpot,grid)

    !=====================================================================!
    ! This subroutine calculates the contribution from the interaction of !
    ! the augmentation charge with the effective potential                !
    !     \hat{D}_ij = \sum_LM \int \tilde{v}_eff(r) \hat{Q}_ij^LM (r) dr !
    ! to the nonlocal term Dij of the Hamiltonian, where \tilde{v}_eff(r) !
    ! is the effective potential resulting from the smooth part of the    !
    ! density, v_H[\tilde{n}+\hat{n}+\tilde{n}_Zc]                        !
    !         +v_xc[\tilde{n}+\hat{n}+\tilde{n}_c]                        !
    !---------------------------------------------------------------------!
    ! Arguments:                                                          !
    !   dij(inout)    : Compensation density contribution to nonlocal     !
    !                   energy term \hat{D}_ij.                           !
    !   locpot(in)    : Local potential on fine grid.                     !
    !   grid(in)      : Grid definition for fine grid.                    !
    !---------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30 May 2010.                            !
    !=====================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_extract_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: pub_my_node_id, pub_on_root
    use constants, only: max_spins
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use geometry, only: geometry_magnitude
    use parallel_strategy, only: pub_num_atoms_on_node, &
         pub_first_atom_on_node, pub_elements_on_node, pub_max_atoms_on_node
    use paw, only: paw_atom_aug_integrals
    use pseudopotentials, only: pseudo_atom_aug_integrals
    use rundat, only: pub_paw, pub_usp, pub_aug_funcs_recip
    use simulation_cell, only: pub_cell, pub_aug_box_n1, pub_aug_box_n2, &
         pub_aug_box_n3
    use sparse, only: SPAM3, sparse_get_block, sparse_put_block
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, &
         utils_sanity_check

    implicit none

    ! Arguments
    type(GRID_INFO),intent(in) :: grid
    type(SPAM3),intent(inout) :: dij(:)
    real(kind=DP),intent(in) :: locpot(grid%ld1,grid%ld2, &
         grid%max_slabs12,size(dij))

    ! Local Variables
    integer :: loc_iat, iat
    integer :: isp
    integer :: is
    integer :: ierr
    integer :: box_n1
    integer :: box_n2
    integer :: box_n3
    integer :: box_start1
    integer :: box_start2
    integer :: box_start3
    logical :: i_need_box
    integer :: num_spins
    real(kind=DP), allocatable :: dij_at(:,:,:)
    real(kind=DP), allocatable :: locpot_box(:,:,:,:)
    real(kind=DP), allocatable :: atom_aug_func(:,:,:)
    real(kind=DP), allocatable :: buffer(:,:,:)
    complex(kind=DP), allocatable :: atom_aug_func_recip(:,:,:)
    complex(kind=DP), allocatable :: locpot_box_recip(:,:,:,:)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &augmentation_screen_dij'
    call utils_sanity_check(locpot(:,:,:,1),'locpot')
#endif

    ! Start Timer
    call timer_clock('augmentation_screen_dij',1)

    if (size(locpot,4)/=size(dij,1)) then
       call utils_abort('Error in augmentation_screen_dij: Inconsistent &
            &sizes of locpot and dij')
    end if
    num_spins = size(dij,1)
    isp = 0

    ! Find size of box
    box_n1 = pub_aug_box_n1
    box_n2 = pub_aug_box_n2
    box_n3 = pub_aug_box_n3

    ! Allocate temporary arrays
    allocate(dij_at(max_proj_tot,max_proj_tot, &
         num_spins),stat=ierr)
    call utils_alloc_check('augmentation_screen_dij','dij_at',ierr)
    allocate(locpot_box(box_n1,box_n2,box_n3,num_spins),stat=ierr)
    call utils_alloc_check('augmentation_screen_dij','locpot_box',ierr)
    allocate(buffer(box_n1,box_n2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('augmentation_screen_dij','buffer',ierr)
    if (pub_aug_funcs_recip) then
       allocate(atom_aug_func_recip(box_n1,box_n2,box_n3),stat=ierr)
       call utils_alloc_check('augmentation_screen_dij','atom_aug_func_recip',ierr)
       allocate(locpot_box_recip(box_n1,box_n2,box_n3,num_spins),stat=ierr)
       call utils_alloc_check('augmentation_screen_dij','locpot_box_recip',ierr)
    else
       allocate(atom_aug_func(box_n1,box_n2,box_n3),stat=ierr)
       call utils_alloc_check('augmentation_screen_dij','atom_aug_func',ierr)
    end if

    ! Loop over atoms
    do loc_iat=1,pub_max_atoms_on_node

       ! Only need to extract if there is an atom left on this node
       if (loc_iat<=pub_num_atoms_on_node(pub_my_node_id)) then
          iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1
          isp = pub_elements_on_node(loc_iat)%pspecies_number

          ! Find where box for this atom is located in simulation cell
          call cell_grid_box_start_wrt_atom( &
               box_start1, box_start2, box_start3, &
               pub_elements_on_node(loc_iat)%centre, box_n1, box_n2, box_n3, &
               grid)

          i_need_box = .true.
       else
          i_need_box = .false.
       end if

       ! Extract tightbox of data from effective potential over simulation
       ! cell for this atom
       do is=1,num_spins
          call cell_grid_extract_box(locpot_box(:,:,:,is), &
               buffer, locpot(:,:,:,is), grid, &
               box_n1, box_n2, box_n3, box_n1, box_n2, &
               box_start1, box_start2, box_start3, i_need_box, .false.)
       end do  ! is
#ifdef DEBUG
       call utils_sanity_check(locpot_box(:,:,:,1),'locpot_box')
#endif

       ! Only need to take overlap if there is an atom left on this node
       if (loc_iat<=pub_num_atoms_on_node(pub_my_node_id)) then

          ! Get block of dij from SPAM3 matrix
          !do is=1,num_spins
          !   call sparse_get_block(dij_at(:,:,is),dij(is),iat,iat)
          !end do  ! is
          dij_at(:,:,:) = 0.0_DP

          if (pub_aug_funcs_recip) then
             do is=1,num_spins
                locpot_box_recip(:,:,:,is) = locpot_box(:,:,:,is)
                call fourier_apply_box('F','B',locpot_box_recip(:,:,:,is), &
                     aug=.true.)
             end do
          end if

          ! Call appropriate routine to calculate the integral of the
          ! augmentation function for this atom with the local potential
          ! in the augmentation box
          if (pub_paw) then
             if (pub_aug_funcs_recip) then
                call paw_atom_aug_integrals(dij_at,num_spins,isp, &
                     pub_elements_on_node(loc_iat)%centre,grid, &
                     box_n1,box_n2,box_n3,box_start1,box_start2,box_start3, &
                     locpot_box_recip=locpot_box_recip, &
                     atom_aug_func_recip=atom_aug_func_recip)
             else
                call paw_atom_aug_integrals(dij_at,num_spins,isp, &
                     pub_elements_on_node(loc_iat)%centre,grid, &
                     box_n1,box_n2,box_n3,box_start1,box_start2,box_start3, &
                     locpot_box,atom_aug_func)
             end if
          else if (pub_usp) then
             ! NOT DONE YET
             call pseudo_atom_aug_integrals(locpot_box,atom_aug_func, &
                  dij_at,num_spins,isp,pub_elements_on_node(loc_iat)%centre, &
                  grid,box_n1,box_n2,box_n3,box_start1,box_start2,box_start3)
          end if

          ! Put block of dij into SPAM3 matrix
          do is=1,num_spins
             call sparse_put_block(dij_at(:,:,is),dij(is),iat,iat)
          end do  ! is

       end if  ! loc_iat < loc_nat
    end do  ! loc_iat

    ! Deallocate temporary arrays
    if (pub_aug_funcs_recip) then
       deallocate(locpot_box_recip,stat=ierr)
       call utils_dealloc_check('augmentation_screen_dij','locpot_box_recip',ierr)
       deallocate(atom_aug_func_recip,stat=ierr)
       call utils_dealloc_check('augmentation_screen_dij','atom_aug_func_recip',ierr)
    else
       deallocate(atom_aug_func,stat=ierr)
       call utils_dealloc_check('augmentation_screen_dij','atom_aug_func',ierr)
    end if
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('augmentation_screen_dij','buffer',ierr)
    deallocate(locpot_box,stat=ierr)
    call utils_dealloc_check('augmentation_screen_dij','locpot_box',ierr)
    deallocate(dij_at,stat=ierr)
    call utils_dealloc_check('augmentation_screen_dij','dij_at',ierr)

    ! Stop Timer
    call timer_clock('augmentation_screen_dij',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving augmentation_screen_dij'
#endif

  end subroutine augmentation_screen_dij


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine aug_nl_calculate_forces(nlps_forces, &
       ngwfs_on_grid,ngwf_basis,proj_basis, &
       sp_overlap,inv_overlap,pur_denskern,ham,dijhat)

    !=========================================================================!
    ! This subroutine calculates the contribution to the ionic forces coming  !
    ! from the nonlocal projectors.                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  nlps_forces     : output : nonlocal forces                             !
    !  sp_overlap      : input  : ngwf-projector overlap matrix               !
    !  ngwfs_on_grid   : input  : NGWF data in ppd format                     !
    !  ngwf_basis      : input  : Function basis description for NGWFs        !
    !  proj_basis      : input  : Function basis description for projectors   !
    !  pur_denskern    : input  : purified density kernel SPAM3               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 08/06/2010.                                 !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_barrier, comms_reduce, pub_on_root, pub_my_node_id
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use pseudopotentials, only: pseudo_get_dij, pseudo_aug_Q_matrix, &
         nlps_projectors
    use parallel_strategy, only: pub_first_atom_on_node, pub_orig_atom
    use paw, only: paw_projectors, paw_projector_overlap, paw_nonlocal_energies
    use projectors, only: projectors_func_grad_ovlp_box
    use rundat, only: pub_paw, pub_usp
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_axpy, sparse_create, sparse_destroy, &
         sparse_get_element, sparse_product, sparse_transpose, sparse_scale, &
         sparse_transpose_structure
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: nlps_forces(1:3,pub_cell%nat)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(SPAM3), intent(in) :: sp_overlap
    type(SPAM3), intent(in) :: inv_overlap
    type(SPAM3), intent(in) :: pur_denskern(pub_cell%num_spins)
    type(SPAM3), intent(in) :: ham(pub_cell%num_spins)
    type(SPAM3), intent(in) :: dijhat(pub_cell%num_spins)

    ! Local Variables
    type(SPAM3) :: nl_force_mat(3)
    type(SPAM3) :: siGp_overlap,iGps_overlap,ps_overlap
    type(SPAM3) :: rkq,dij_rkq,kh,khsq,oij
    type(SPAM3),allocatable :: dij(:)
    type(SPAM3),allocatable :: rho_ij(:)
    type(SPAM3),allocatable :: khs(:)
    integer :: cart, is
    integer :: iat, orig_iat
    integer :: atom_proj, global_proj
    integer :: ierr
    real(kind=DP) :: proj_force

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering aug_nl_calculate_forces'
#endif

    ! Start timer
    call timer_clock('aug_nl_calculate_forces',1)

    ! Initialise
    nlps_forces  = 0.0_DP

    ! Create result matrices to hold nonlocal forces, projector-by-projector
    do cart=1,3
       nl_force_mat(cart)%structure = 'E'
       call sparse_create(nl_force_mat(cart))
    end do

    ! Create matrix arrays and structures
    oij%structure = 'E'
    call sparse_create(oij)
    allocate(dij(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('aug_nl_calculate_forces','dij',ierr)
    allocate(rho_ij(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('aug_nl_calculate_forces','rho_ij',ierr)
    do is=1,pub_cell%num_spins
       dij(is)%structure = 'E'
       call sparse_create(dij(is))
    end do
    do is=1,pub_cell%num_spins
       call sparse_create(rho_ij(is),dij(is))
    end do

    if (pub_paw) then

       ! Get projector overlap matrix
       call paw_projector_overlap(oij)

       ! Get projector density kernel
       call aug_projector_denskern(rho_ij,pur_denskern,sp_overlap)

       ! Get nonlocal energies in diagonal matrix
       call paw_nonlocal_energies(dij,rho_ij)
       do is=1,pub_cell%num_spins
          call sparse_axpy(dij(is),dijhat(is),1.0_DP)
       end do

    else if (pub_usp) then

       ! Get projector overlap matrix
       call pseudo_aug_Q_matrix(oij)

       ! Get nonlocal energies in diagonal matrix
       do is=1,pub_cell%num_spins
          call pseudo_get_dij(dij(is))
          call sparse_axpy(dij(is),dijhat(is),1.0_DP)
       end do

    end if

    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(rho_ij(is))
    end do
    deallocate(rho_ij,stat=ierr)
    call utils_dealloc_check('aug_nl_calculate_forces','rho_ij',ierr)

    ! Create matrices to hold <phi|iG*proj> overlap matrix and transpose
    call sparse_transpose_structure(ps_overlap%structure,sp_overlap)
    call sparse_create(ps_overlap)
    call sparse_transpose(ps_overlap,sp_overlap)
    call sparse_create(siGp_overlap,sp_overlap)
    call sparse_create(iGps_overlap,ps_overlap)

    ! Create temporary matrices rkq and dij_rkq
    call sparse_create(rkq,nl_force_mat(1))
    call sparse_create(dij_rkq,rkq)

    ! Calculate K H S^-1
    allocate(khs(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('aug_nl_calculate_forces','khs',ierr)
    call sparse_create(kh,pur_denskern(1),ham(1))
    call sparse_destroy(kh)
    do is=1,pub_cell%num_spins
       call sparse_create(khs(is),kh,inv_overlap)
    end do
    call sparse_create(kh)
    do is=1,pub_cell%num_spins
       call sparse_product(kh,pur_denskern(is),ham(is))
       call sparse_product(khs(is),kh,inv_overlap)
    end do
    call sparse_destroy(kh)
    call sparse_create(khsq,khs(1),sp_overlap)

    ! Loop over Cartesian directions
    do cart=1,3

       ! Calculate <phi|iG*proj> overlap matrix
       if (pub_paw) then
          call projectors_func_grad_ovlp_box(siGp_overlap, &
               ngwfs_on_grid,ngwf_basis,proj_basis,paw_projectors,cart)
       else if (pub_usp) then
          call projectors_func_grad_ovlp_box(siGp_overlap, &
               ngwfs_on_grid,ngwf_basis,proj_basis,nlps_projectors,cart)
       end if

       ! Transpose it to get <iG*proj|phi> overlap matrix
       call sparse_transpose(iGps_overlap,siGp_overlap)

       do is=1,pub_cell%num_spins

          ! Calculate D_ij <proj_j|phi_a> K^ab <phi_b|iG.proj_i>
          call sparse_product(khsq,pur_denskern(is),siGp_overlap)
          call sparse_product(rkq,ps_overlap,khsq)
          call sparse_product(dij_rkq,dij(is),rkq)
          call sparse_axpy(nl_force_mat(cart),dij_rkq,1.0_DP)

          ! Calculate <iG.proj_i|phi_a> K^ab <phi_b|proj_j> D_ji
          call sparse_product(khsq,pur_denskern(is),sp_overlap)
          call sparse_product(rkq,iGps_overlap,khsq)
          call sparse_product(dij_rkq,rkq,dij(is))
          call sparse_axpy(nl_force_mat(cart),dij_rkq,1.0_DP)

          ! Calculate O_ij <proj_j|phi_c>K^cd H_da S^ab <phi_b|iG.proj_i>
          call sparse_product(khsq,khs(is),siGp_overlap)
          call sparse_product(rkq,ps_overlap,khsq)
          call sparse_product(dij_rkq,oij,rkq)
          call sparse_axpy(nl_force_mat(cart),dij_rkq,-1.0_DP)

          ! Calculate <iG.proj_i|phi_a> K^ag H_gd S^da <phi_b|proj_j> O_ji
          call sparse_product(khsq,khs(is),sp_overlap)
          call sparse_product(rkq,iGps_overlap,khsq)
          call sparse_product(dij_rkq,rkq,oij)
          call sparse_axpy(nl_force_mat(cart),dij_rkq,-1.0_DP)

       end do
    end do

    ! Destroy temporary matrices and deallocate arrays
    call sparse_destroy(khsq)
    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(khs(is))
    end do
    deallocate(khs,stat=ierr)
    call utils_dealloc_check('aug_nl_calculate_forces','khs',ierr)
    call sparse_destroy(dij_rkq)
    call sparse_destroy(rkq)
    call sparse_destroy(iGps_overlap)
    call sparse_destroy(siGp_overlap)
    call sparse_destroy(ps_overlap)
    call sparse_destroy(oij)
    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(dij(is))
    end do
    deallocate(dij,stat=ierr)
    call utils_dealloc_check('aug_nl_calculate_forces','dij',ierr)

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
    call timer_clock('aug_nl_calculate_forces',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving aug_nl_calculate_forces'
#endif

  end subroutine aug_nl_calculate_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine augmentation_density_forces(nhat_forces,denskern,sp_overlap, &
       locpot,grid)

    !=========================================================================!
    ! This subroutine calculates the contribution to the forces on each atom  !
    ! from the interaction of the compensation density nhat with the          !
    ! effective potential.                                                    !
    !  F_I = -d/dR_I(\int v_eff(r) n_hat(r) dr)                               !
    !      = \int v_eff(r) \sum_ijLM d/dR_I(g_L(r)S_LM(r)) rho_ij q_ij^LM dr  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 30 May 2010.                                !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_extract_box, &
         cell_grid_box_start_wrt_atom
    use comms, only: comms_barrier, comms_reduce, pub_my_node_id, pub_on_root
    use constants, only: max_spins
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use geometry, only: geometry_magnitude
    use parallel_strategy, only: pub_num_atoms_on_node, &
         pub_first_atom_on_node, pub_elements_on_node, pub_max_atoms_on_node, &
         pub_orig_atom
    use paw, only: paw_atom_aug_force
    use pseudopotentials, only: pseudo_atom_aug_force
    use rundat, only: pub_paw, pub_usp, pub_aug_funcs_recip
    use simulation_cell, only: pub_cell, pub_aug_box_n1, pub_aug_box_n2, &
         pub_aug_box_n3
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_get_block
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, utils_flush ! TEMP

    implicit none

    ! Arguments
    real(kind=DP),intent(out) :: nhat_forces(3,pub_cell%nat)
    type(GRID_INFO),intent(in) :: grid
    real(kind=DP),intent(in) :: locpot(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
    type(SPAM3),intent(in) :: denskern(pub_cell%num_spins)
    type(SPAM3),intent(in) :: sp_overlap

    ! Local Variables
    type(SPAM3), allocatable :: rho_ij(:)
    integer :: loc_iat, orig_iat, iat
    integer :: isp
    integer :: is
    integer :: ierr
    integer :: box_n1
    integer :: box_n2
    integer :: box_n3
    integer :: box_start1
    integer :: box_start2
    integer :: box_start3
    logical :: i_need_box
    integer :: num_spins
    real(kind=DP), allocatable :: buffer(:,:,:)
    real(kind=DP), allocatable :: rhoij_at(:,:,:)
    real(kind=DP), allocatable :: locpot_box(:,:,:,:)
    real(kind=DP), allocatable :: atom_grad_aug_func(:,:,:,:)
    complex(kind=DP), allocatable :: locpot_box_recip(:,:,:,:)
    complex(kind=DP), allocatable :: atom_grad_aug_func_recip(:,:,:,:)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering augmentation_density_forces'
#endif

    ! Start Timer
    call timer_clock('augmentation_density_forces',1)

    num_spins = pub_cell%num_spins
    isp = 0

    ! Find size of box
    box_n1 = pub_aug_box_n1
    box_n2 = pub_aug_box_n2
    box_n3 = pub_aug_box_n3

    ! Create projector density kernel
    allocate(rho_ij(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('augmentation_density_forces','rho_ij',ierr)

    do is=1,pub_cell%num_spins
       rho_ij(is)%structure = 'E'
       call sparse_create(rho_ij(is))
    end do  ! is

    call aug_projector_denskern(rho_ij,denskern,sp_overlap)

    ! Allocate temporary arrays
    allocate(rhoij_at(max_proj_tot,max_proj_tot,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('augmentation_density_forces','rhoij_at',ierr)
    allocate(locpot_box(box_n1,box_n2,box_n3,num_spins),stat=ierr)
    call utils_alloc_check('augmentation_density_forces','locpot_box',ierr)
    allocate(buffer(box_n1,box_n2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('augmentation_density_forces','buffer',ierr)
    if (pub_aug_funcs_recip) then
       allocate(locpot_box_recip(box_n1,box_n2,box_n3,num_spins),stat=ierr)
       call utils_alloc_check('augmentation_density_forces', &
            'locpot_box_recip',ierr)
       allocate(atom_grad_aug_func_recip(box_n1,box_n2,box_n3,3),stat=ierr)
       call utils_alloc_check('augmentation_density_forces', &
            'atom_grad_aug_func_recip',ierr)
    else
       allocate(atom_grad_aug_func(box_n1,box_n2,box_n3,3),stat=ierr)
       call utils_alloc_check('augmentation_density_forces', &
            'atom_grad_aug_func',ierr)
    end if

    ! Initialise
    nhat_forces = 0.0_DP

    ! Loop over atoms
    do loc_iat=1,pub_max_atoms_on_node

       ! Only need to extract if there is an atom left on this node
       if (loc_iat<=pub_num_atoms_on_node(pub_my_node_id)) then

          ! Find atom number in input file order and species number
          iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1
          isp = pub_elements_on_node(loc_iat)%pspecies_number
          orig_iat = pub_orig_atom(iat)

          ! Find where box for this atom is located in simulation cell
          call cell_grid_box_start_wrt_atom( &
               box_start1, box_start2, box_start3, &
               pub_elements_on_node(loc_iat)%centre, box_n1, box_n2, box_n3, &
               grid)

          i_need_box = .true.
       else
          i_need_box = .false.
          orig_iat = -1  ! suppress warning
       end if ! loc_iat < loc_nat

       ! Extract tightbox of data from effective potential over simulation
       ! cell for this atom
       do is=1,num_spins
          call cell_grid_extract_box(locpot_box(:,:,:,is), &
               buffer, locpot(:,:,:,is), grid, &
               box_n1, box_n2, box_n3, box_n1, box_n2, &
               box_start1, box_start2, box_start3, i_need_box, .false.)
       end do  ! is

       ! Fourier transform the locpot box if required
       if (pub_aug_funcs_recip) then
          do is=1,num_spins
             locpot_box_recip(:,:,:,is) = locpot_box(:,:,:,is)
             call fourier_apply_box('F','B',locpot_box_recip(:,:,:,is), &
                  aug=.true.)
          end do
       end if

       ! Only need to calculate force if there is an atom left on this node
       if (loc_iat<=pub_num_atoms_on_node(pub_my_node_id)) then

          ! Get blocks of projector density matrix for this atom
          do is=1,pub_cell%num_spins
             call sparse_get_block(rhoij_at(:,:,is),rho_ij(is),iat,iat)
          end do  ! is

          ! Call appropriate routine to calculate the contribution to the
          ! augmentation density force for this atom using the augmentation box
          if (pub_paw) then
             if (pub_aug_funcs_recip) then
                call paw_atom_aug_force(nhat_forces(:,orig_iat), &
                     rhoij_at,num_spins,isp, &
                     pub_elements_on_node(loc_iat)%centre, &
                     grid,box_n1,box_n2,box_n3,box_start1,box_start2, &
                     box_start3,locpot_box_recip=locpot_box_recip, &
                     atom_grad_aug_func_recip=atom_grad_aug_func_recip)
             else
                call paw_atom_aug_force(nhat_forces(:,orig_iat), &
                     rhoij_at,num_spins,isp, &
                     pub_elements_on_node(loc_iat)%centre, &
                     grid,box_n1,box_n2,box_n3,box_start1,box_start2, &
                     box_start3,locpot_box_real=locpot_box, &
                     atom_grad_aug_func_real=atom_grad_aug_func)
             end if
          else if (pub_usp) then
             call pseudo_atom_aug_force(nhat_forces(:,orig_iat),locpot_box, &
                  atom_grad_aug_func,rhoij_at,num_spins,isp, &
                  pub_elements_on_node(loc_iat)%centre,grid,box_n1,box_n2,box_n3, &
                  box_start1,box_start2,box_start3)
          end if

       end if  ! loc_iat < loc_nat
    end do  ! loc_iat

    ! Reduce result across nodes
    call comms_barrier
    call comms_reduce('SUM',nhat_forces,3*pub_cell%nat)

    ! Deallocate temporary arrays
    if (pub_aug_funcs_recip) then
       deallocate(atom_grad_aug_func_recip,stat=ierr)
       call utils_dealloc_check('augmentation_density_forces', &
            'atom_grad_aug_func_recip',ierr)
       deallocate(locpot_box_recip,stat=ierr)
       call utils_dealloc_check('augmentation_density_forces', &
            'locpot_box_recip',ierr)
    else
       deallocate(atom_grad_aug_func,stat=ierr)
       call utils_dealloc_check('augmentation_density_forces', &
            'atom_grad_aug_func',ierr)
    end if
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('augmentation_density_forces','buffer',ierr)
    deallocate(locpot_box,stat=ierr)
    call utils_dealloc_check('augmentation_density_forces','locpot_box',ierr)
    deallocate(rhoij_at,stat=ierr)
    call utils_dealloc_check('augmentation_density_forces','rhoij_at',ierr)

    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(rho_ij(is))
    end do  ! is

    deallocate(rho_ij,stat=ierr)
    call utils_dealloc_check('augmentation_density_forces','rho_ij',ierr)

    ! Stop Timer
    call timer_clock('augmentation_density_forces',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving augmentation_density_forces'
#endif

  end subroutine augmentation_density_forces

end module augmentation
