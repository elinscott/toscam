! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-

module nonsc_forces

  !=====================================================================!
  ! This module contains the subroutines to calculate the correction to !
  ! the ONETEP forces due to non full self-consistency of the NGWF      !
  ! optimisation.                                                       !
  ! If the NGWFs are optimised to a tight threshold (as in a standard   !
  ! calculation), then the contribution of the non self-consistent force!
  ! should be residual. However, if the NGWFs are not totally optimised !
  ! this force should be taken into account.                            !
  !---------------------------------------------------------------------!
  ! *** Please note that this capability is still under development.    !
  !---------------------------------------------------------------------!
  ! This module was created by Alvaro Ruiz Serrano in November 2010.    !
  !=====================================================================!

  implicit none

  public :: nonsc_forces_ngwfs_calc
  public :: nonsc_forces_overlap_deriv


contains

  subroutine nonsc_forces_ngwfs_calc(nonsc_forces,&
       ngwfs_on_grid, contra_grad_on_grid, ngwf_basis)

    !======================================================================!
    ! This subroutine calculates the forces due to non self-consistency of !
    ! the outer loop. These forces behave like Pulay forces as if the NGWF !
    ! functions were the ultimate basis set.                               !
    ! See: J. Chem. Phys. vol 121, no 13, October 2004.                    !
    !======================================================================!
    !  Arguments:                                                          !
    !    nonsc_forces (out)  : NGWF non self-consistent forces in          !
    !                              cartesian coordinates array             !
    !    ngwfs_on_grid  (in) : current set of NGWFS on grid in PPD format  !
    !    contra_gradient (in): contravariant NGFW gradient on grid (PPD)   !
    !    ngwf_basis (in)     : the basis type for the NGWFs                !
    !======================================================================!
    ! Written by Alvaro Ruiz Serrano in January 2010.                      !
    ! Clean-up by Alvaro Ruiz Serrano in December 2010.                    !
    !======================================================================!

    use constants, only: DP
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_grad
    use kinetic, only: kinetic_grad_to_func_batch
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_scale
    use timer, only: timer_clock

    implicit none

    ! Arguments
    real(kind=DP),   intent(out) :: nonsc_forces(:,:)
    type(FUNC_BASIS),intent(in) :: ngwf_basis
    real(kind=DP),intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP),intent(in) :: contra_grad_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)

    ! Local Variables
    ! ars: nabla matrices. The diagonal contains the forces
    type(SPAM3) :: nabla(3)
    integer :: ndim


    ! Start timer
    call timer_clock('nonsc_forces_ngwfs_calc',1)

    ! ars: create nablas
    do ndim = 1, 3
       nabla(ndim)%structure = 'D'
       call sparse_create(nabla(ndim))
    end do

    ! ars: calculate nabla = <contra_grad|d(ngwf)/dx>
    call integrals_grad(nabla, contra_grad_on_grid, ngwf_basis, ngwfs_on_grid, &
         ngwf_basis)

    ! ars: at this point, the deriv of the NGWFs should be scaled according to
    !      d(ngwf)/dx = - d(ngwf)/dX_alpha
    !      and then, when extracting the forces,
    !      F = - sum_alpha <contra_grad|d(ngwf)/dX_alpha>
    ! ---> avoid multiply twice times -1

    ! ars: extract nonsc_forces forces
    call nonsc_forces_extract_forces(nonsc_forces, nabla, ngwf_basis)

    ! ars: destroy nablas
    do ndim = 1, 3
       call sparse_destroy(nabla(ndim))
    end do

    ! Stop timer
    call timer_clock('nonsc_forces_ngwfs_calc',2)


  end subroutine nonsc_forces_ngwfs_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nonsc_forces_overlap_deriv(ov_deriv, ngwfs_on_grid, ngwf_basis)

    !======================================================================!
    ! This subroutine calculates the derivative of the overlap matrix wrt  !
    ! the atomic positions in NGWF representation and stores the data in   !
    ! SPAM3 omega matrices corresponding to X, Y and Z coordinates.        !
    !======================================================================!
    !  Arguments:                                                          !
    !    ov_deriv (inout) : deriv of the overlap matrix wrt atpos elements !
    !    ngwfs_on_grid(in): current set of NGWFS on grid in PPD format     !
    !    ngwf_basis (in)  : The basis type for the NGWFs                   !
    !======================================================================!
    ! Written by Alvaro Ruiz Serrano in November 2009                      !
    ! Clean-up by Alvaro Ruiz Serrano in December 2010.                    !
    !======================================================================!

    use constants, only: DP
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_grad
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_scale
    use timer, only: timer_clock


    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: ov_deriv(3)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)

    ! Local Variables
    integer :: ndim


    ! Start timer
    call timer_clock('nonsc_forces_overlap_deriv',1)

    ! ars: calculate ov_deriv = d/dx S_{alpha,beta}
    call integrals_grad(ov_deriv, ngwfs_on_grid, ngwf_basis, ngwfs_on_grid, &
         ngwf_basis)

    ! ars: scale according to
    !      d(ngwf)/dx = - d(ngwf)/dX_alpha
    do ndim = 1, 3
       call sparse_scale(ov_deriv(ndim), -1.0_DP)
    end do

    ! ars: add transpose to diagonal blocks
    call nonsc_forces_diagonal_blocks(ov_deriv, ngwf_basis)

    ! Stop timer
    call timer_clock('nonsc_forces_overlap_deriv',2)


  end subroutine nonsc_forces_overlap_deriv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nonsc_forces_extract_forces(nonsc_forces, nabla, ngwf_basis)

    !======================================================================!
    ! This subroutine extracts the nonsc_forces forces from the diagonal   !
    ! terms of nablaXYZ matrices:                                          !
    ! F_atom = -sum(NGWFs_on atom) * nabla_atomblock(NGWF, NGWF)           !
    !======================================================================!
    !  Arguments:                                                          !
    !    nonsc_forces (out)  : nonsc_forces in cartesian coordinates array !
    !    nabla        (in)   : matrices whose diagonal elements contain the!
    !                          forces. The rest of the elements are junk.  !
    !======================================================================!
    ! Written by Alvaro Ruiz Serrano in January 2010.                      !
    ! Clean-up by Alvaro Ruiz Serrano in December 2010.                    !
    !======================================================================!

    use comms, only: pub_my_node_id, pub_on_root, comms_reduce, comms_barrier
    use constants, only : DP
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: pub_first_atom_on_node,&
         pub_orig_atom, pub_num_atoms_on_node, pub_elements_on_node
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_get_block, SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none


    ! ars: arguments
    real(kind=DP), intent(out) :: nonsc_forces(:,:)
    type(SPAM3), intent(in) :: nabla(3)
    type(FUNC_BASIS), intent(in) :: ngwf_basis


    ! ars: local variables
    real(kind=DP), allocatable :: block(:,:)

    integer :: iat, loc_iat, ingwf, ndim, ierr
    real(kind=DP) :: block_tr


    ! ars: initialise
    nonsc_forces(:,:) = 0.0_DP

    allocate(block(ngwf_basis%max_on_atom,ngwf_basis%max_on_atom),stat=ierr)
    call utils_alloc_check('nonsc_forces_extract_forces','block',ierr)

    ! ars: extract forces on each atom
    iat = pub_first_atom_on_node(pub_my_node_id)
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       do ndim = 1, 3
          block(:,:) = 0.0_DP
          block_tr = 0.0_DP
          call sparse_get_block(block,nabla(ndim),iat,iat)

          do ingwf=1,ngwf_basis%num_on_atom(iat)
             block_tr = block_tr + block(ingwf,ingwf)
          end do

          nonsc_forces(ndim,pub_orig_atom(iat)) = block_tr
       end do

       iat = iat + 1
    end do

    ! ars: deallocate memory
    deallocate(block, stat=ierr)
    call utils_dealloc_check('nonsc_forces_extract_forces','block',ierr)

    ! ars: add factors and distribute over the nodes
    nonsc_forces(:,:) = nonsc_forces(:,:)/pub_cell%weight
    call comms_reduce('SUM', nonsc_forces)


  end subroutine nonsc_forces_extract_forces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nonsc_forces_diagonal_blocks(ov_deriv, ngwf_basis)

    !======================================================================!
    ! This subroutine extracts the diagonal atom-blocks of ov_deriv and    !
    ! adds the transpose of themselves to calculate the correct derivative !
    ! of the overlap matrix wrt the centre of the NGWF spheres.            !
    !                                                                      !
    ! ***Note: in the vast majority of calculations, this routine will set !
    ! the diagonal blocks of ov_deriv to zero or almost zero.              !
    !======================================================================!
    !  Arguments:                                                          !
    !   ov_deriv (inout): deriv of the overlap wrt the centre of the sphere!
    !======================================================================!
    ! Originally written by Alvaro Ruiz Serrano in November 2009.          !
    ! Clean-up by Alvaro Ruiz Serrano in December 2010.                    !
    !======================================================================!


    use comms, only: pub_my_node_id
    use constants, only : DP
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: pub_first_atom_on_node, pub_num_atoms_on_node
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_get_block, sparse_put_block, SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check


    implicit none

    ! ars: <<arguments>>
    type(SPAM3), intent(inout) :: ov_deriv(3)
    type(FUNC_BASIS), intent(in) :: ngwf_basis

    ! ars: <<local variables>>
    real(kind=DP), allocatable :: block(:,:)
    integer :: atom, no_ngwfs, local_atom, ierr, ndim


    ! ars: loop over the atoms on this node
    do local_atom = 1,pub_num_atoms_on_node(pub_my_node_id)
       atom = pub_first_atom_on_node(pub_my_node_id) + local_atom - 1

       ! ars: define block
       no_ngwfs = ngwf_basis%num_on_atom(atom)
       allocate(block(no_ngwfs, no_ngwfs), stat=ierr)
       call utils_alloc_check('nonsc_forces_diagonal_blocks','block',ierr)

       do ndim = 1, 3

          ! ars: init block
          block(:,:) = 0.0_DP

          ! ars: extract block corresponding to this atom
          call sparse_get_block(block, ov_deriv(ndim), atom, atom)

          ! ars: add the transpose
          block(:,:) = block(:,:) + transpose(block(:,:))

          ! ars: put block back in omega SPAM3
          call sparse_put_block(block, ov_deriv(ndim), atom, atom)

       end do

       ! ars: deallocate and start new atom
       deallocate(block, stat=ierr)
       call utils_dealloc_check('nonsc_forces_diagonal_blocks','block',ierr)

    end do


  end subroutine nonsc_forces_diagonal_blocks


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module nonsc_forces
