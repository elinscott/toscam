! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                  Density kernel utility module                 !
!                                                                !
! This module implements several utilities for the density       !
! kernel
!----------------------------------------------------------------!
! Written by Peter Haynes, 18/11/04                              !
! Modified by Nicholas Hine, December 2007 to reduce number of   !
! sparse algebra operations                                      !
!================================================================!

module kernel

  use constants, only: DP
  use sparse, only: SPAM3

  implicit none

  ! Re-usable quantities
  type(SPAM3), allocatable, public :: pub_ks(:)  ! Product of kernel and overlap matrix
  type(SPAM3), allocatable, public :: pub_ksk(:) ! Product of kernel and sk

  logical, public :: ks_valid,ksk_valid          ! Flags for validity of pub_ks, pub_ksk

  private

  logical, public :: pub_kernel_workspace_allocated

  public :: kernel_rms_err
  public :: kernel_occupancy_bounds
  public :: kernel_fix
  public :: kernel_normalise
  public :: kernel_rescale
  public :: kernel_purify
  public :: kernel_rms_commutator
  public :: kernel_init_core_ham

  public :: kernel_workspace_allocate
  public :: kernel_workspace_deallocate

  public :: kernel_workspace_invalidate
  public :: kernel_validate_ks
  public :: kernel_validate_ksk

  public :: kernel_from_vecs_asc2

  public :: kernel_basis_transform
  public :: kernel_basis_update
  public :: kernel_christoffel

contains

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_workspace_allocate(denskern,overlap)

    !=======================================================================!
    ! This subroutine allocates the shared pub_ks and pub_ksk workspaces    !
    !=======================================================================!
    ! Arguments:                                                            !
    ! denskern (inout)        : Density kernel in SPAM3 format              !
    ! overlap (input)         : Overlap matrix in SPAM3 format              !
    !=======================================================================!
    ! Written by Nick Hine, November 2007                                   !
    !=======================================================================!

    use simulation_cell, only : pub_cell
    use sparse, only: sparse_create
    use utils, only : utils_alloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: denskern(pub_cell%num_spins)
    type(SPAM3), intent(in) :: overlap

    ! Local Variables
    integer :: is,ierr

    ! Allocate structures for sparse matrices and initialise
    allocate(pub_ks(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('kernel_workspace_allocate', &
         'pub_ks',ierr)
    allocate(pub_ksk(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('kernel_workspace_allocate', &
         'pub_ksk',ierr)
    do is=1,pub_cell%num_spins
       call sparse_create(pub_ks(is),denskern(is),overlap)
       call sparse_create(pub_ksk(is),pub_ks(is),denskern(is))
    end do

    pub_kernel_workspace_allocated = .true.
    call kernel_workspace_invalidate()

  end subroutine kernel_workspace_allocate

  subroutine kernel_workspace_deallocate()

    !=======================================================================!
    ! This subroutine deallocates the shared pub_ks and pub_ksk workspaces  !
    !=======================================================================!
    ! Arguments:                                                            !
    !  None                                                                 !
    !=======================================================================!
    ! Written by Nick Hine, November 2007                                   !
    !=======================================================================!

    use simulation_cell, only : pub_cell
    use sparse, only: sparse_destroy
    use utils, only : utils_dealloc_check

    implicit none

    ! Local Variables
    integer :: is,ierr

    pub_kernel_workspace_allocated=.false.

    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(pub_ksk(is))
       call sparse_destroy(pub_ks(is))
    enddo

    deallocate(pub_ksk,stat=ierr)
    call utils_dealloc_check('kernel_workspace_deallocate', &
         'pub_ksk',ierr)
    deallocate(pub_ks,stat=ierr)
    call utils_dealloc_check('kernel_workspace_deallocate', &
         'pub_ks',ierr)

  end subroutine kernel_workspace_deallocate

  subroutine kernel_workspace_invalidate()

    !=======================================================================!
    ! This subroutine sets the flags of the workspace arrays to invalid     !
    !=======================================================================!
    ! Arguments:                                                            !
    !  None                                                                 !
    !=======================================================================!
    ! Written by Nick Hine, November 2007                                   !
    !=======================================================================!

    implicit none

    ks_valid = .false.
    ksk_valid = .false.

  end subroutine kernel_workspace_invalidate

  subroutine kernel_validate_ks(denskern,overlap)

    !=======================================================================!
    ! This subroutine recalculates ks if required and sets the flag of the  !
    ! workspace array to valid                                              !
    !=======================================================================!
    ! Arguments:                                                            !
    ! denskern (input)        : Density kernel in SPAM3 format              !
    ! overlap (input)         : Overlap matrix in SPAM3 format              !
    !=======================================================================!
    ! Written by Nick Hine, November 2007                                   !
    !=======================================================================!

    use simulation_cell, only : pub_cell
    use sparse, only: SPAM3, sparse_product

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: denskern(pub_cell%num_spins)
    type(SPAM3), intent(in) :: overlap

    ! Local Variables
    integer :: is

    if (ks_valid)then
       ! No need to recalculate it
       return
    endif

    do is=1,pub_cell%num_spins
       call sparse_product(pub_ks(is),denskern(is),overlap)
    enddo

    ks_valid=.true.

  end subroutine kernel_validate_ks

  subroutine kernel_validate_ksk(denskern,overlap)

    !=======================================================================!
    ! This subroutine recalculates ksk if required and sets the flag of the !
    ! workspace array to valid                                              !
    !=======================================================================!
    ! Arguments:                                                            !
    ! denskern (input)        : Density kernel in SPAM3 format              !
    ! overlap (input)         : Overlap matrix in SPAM3 format              !
    !=======================================================================!
    ! Written by Nick Hine, November 2007                                   !
    !=======================================================================!

    use simulation_cell, only : pub_cell
    use sparse, only: SPAM3, sparse_product

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: denskern(pub_cell%num_spins)
    type(SPAM3), intent(in) :: overlap

    ! Local Variables
    integer :: is

    if (ksk_valid)then
       ! No need to recalculate it
       return
    endif

    if (.not.(ks_valid))then
       call kernel_validate_ks(denskern,overlap)
    endif

    do is=1,pub_cell%num_spins
       call sparse_product(pub_ksk(is),pub_ks(is),denskern(is))
    enddo

    ksk_valid=.true.

  end subroutine kernel_validate_ksk

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  real(kind=DP) function kernel_rms_err(denskern, overlap)

    !=======================================================================!
    ! This subroutine returns the rms occupancy error of a given density    !
    ! kernel by calculating the value of the penalty functional.            !
    !                                                                       !
    ! See P.D. Haynes and M.C. Payne, Phys. Rev. B 59, 12173 (1999).        !
    !=======================================================================!
    ! Arguments:                                                            !
    ! denskern (inout)        : Density kernel in SPAM3 format              !
    ! overlap (input)         : Overlap matrix in SPAM3 format              !
    !=======================================================================!
    ! Written by Nicholas Hine, September 2009.                             !
    !=======================================================================!

    use constants, only: DP, max_spins
    use simulation_cell, only : pub_cell
    use sparse, only: sparse_axpy, sparse_copy, sparse_create, sparse_destroy, &
         sparse_num_rows, sparse_product, sparse_trace, sparse_scale

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: denskern(pub_cell%num_spins)  ! Density kernel
    type(SPAM3), intent(in) :: overlap                       ! Overlap matrix

    ! Local Variables
    type(SPAM3) :: ks, ksks
    real(kind=DP) :: pen
    integer :: is
    logical :: deallocate_workspace


    ! Create workspace if required
    deallocate_workspace = .false.
    if (.not.pub_kernel_workspace_allocated) then
       deallocate_workspace = .true.
       call kernel_workspace_allocate(denskern,overlap)
    end if

    ! Initialise
    pen = 0.0_DP

    ! Create temporary matrices
    call sparse_create(ks,pub_ks(1))
    call sparse_create(ksks,ks,ks)

    ! ndmh: recalculate K.S if it has become invalid
    call kernel_validate_ks(denskern,overlap)

    do is=1,pub_cell%num_spins

       ! Calculate ks := 1 - K.S
       call sparse_copy(ks,pub_ks(is))
       call sparse_scale(ks,-1.0_DP,1.0_DP)

       ! Calculate P = Tr[KSKS(1-KS)(1-KS)]
       call sparse_product(ksks,pub_ks(is),ks)
       pen = pen + sparse_trace(ksks,ksks)

    end do

    pen = pen / pub_cell%num_spins
    kernel_rms_err = sqrt(abs(pen)/sparse_num_rows(denskern(1)))

    call sparse_destroy(ksks)
    call sparse_destroy(ks)

    ! Deallocate workspace if we created it on entry
    if (deallocate_workspace) call kernel_workspace_deallocate


  end function kernel_rms_err

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_occupancy_bounds(max_occ, min_occ, denskern, overlap)

    !==========================================================================!
    ! This subroutine returns the bounds on the occupancies of a given density !
    ! kernel by calculating the maximum and minimum eigenvalues of KS.         !
    !==========================================================================!
    ! Arguments:                                                               !
    ! denskern (inout)        : Density kernel in SPAM3 format                 !
    ! overlap (input)         : Overlap matrix in SPAM3 format                 !
    !==========================================================================!
    ! Written by Nicholas Hine, September 2009.                                !
    !==========================================================================!

    use constants, only: DP, max_spins
    use simulation_cell, only : pub_cell
    use sparse, only: sparse_copy, sparse_create, sparse_destroy, &
         sparse_extremal_eigenvalue, sparse_scale

    implicit none

    ! Arguments
    real(kind=DP),intent(out) :: max_occ(max_spins)         ! Max KS eigenvalue
    real(kind=DP),intent(out) :: min_occ(max_spins)         ! Min KS eigenvalue
    type(SPAM3), intent(in) :: denskern(pub_cell%num_spins) ! Density kernel
    type(SPAM3), intent(in) :: overlap                      ! Overlap matrix

    ! Local Variables
    integer :: is
    logical :: deallocate_workspace

    ! Create workspace if required
    deallocate_workspace = .false.
    if (.not.pub_kernel_workspace_allocated) then
       deallocate_workspace = .true.
       call kernel_workspace_allocate(denskern,overlap)
    end if

    ! Recalculate K.S if required
    call kernel_validate_ks(denskern,overlap)

    ! Find maximum eigenvalues of KS and 1-KS
    do is=1,pub_cell%num_spins
       call sparse_extremal_eigenvalue(pub_ks(is),overlap,max_occ(is),0.001_DP)
       call sparse_scale(pub_ks(is),-1.0_DP,1.0_DP)
       call sparse_extremal_eigenvalue(pub_ks(is),overlap,min_occ(is),0.001_DP)
       min_occ(is) = 1.0_DP - min_occ(is)
       ! Restore ks to avoid recalculating it next time
       call sparse_scale(pub_ks(is),-1.0_DP,1.0_DP)
    end do
    if (pub_cell%num_spins == 1) then
       max_occ(2) = max_occ(1)
       min_occ(2) = min_occ(1)
    end if

    ! Deallocate workspace if we created it on entry
    if (deallocate_workspace) call kernel_workspace_deallocate

  end subroutine kernel_occupancy_bounds

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_fix(denskern, overlap, inv_overlap, miniter, const_ne)

    !=======================================================================!
    ! This subroutine iteratively improves the idempotency of the density   !
    ! kernel (for a fixed set of NGWFs) using the penalty functional.       !
    !                                                                       !
    ! See P.D. Haynes and M.C. Payne, Phys. Rev. B 59, 12173 (1999).        !
    !                                                                       !
    !=======================================================================!
    ! Arguments:                                                            !
    ! denskern (inout)        : Density kernel in SPAM3 format              !
    ! overlap (input)         : Overlap matrix in SPAM3 format              !
    ! inv_overlap (input)     : Inverse overlap matrix in SPAM3 format      !
    ! miniter (input)         : Minimum number of iterations to do          !
    ! const_ne (input)        : Flag to conserve number of electrons        !
    !=======================================================================!
    ! Written by Peter Haynes, November 2004                                !
    ! Improvements by Chris-Kriton Skylaris, June 2009                      !
    ! Speed improvements by Nicholas Hine, November 2009                    !
    !=======================================================================!

    use comms, only: pub_on_root, comms_bcast, pub_root_node_id
    use constants, only: stdout, DP, max_spins, VERBOSE
    use rundat, only: pub_output_detail, pub_kerfix, pub_maxit_kernel_fix
    use simulation_cell, only : pub_cell
    use sparse, only: SPAM3, sparse_extremal_eigenvalue, sparse_copy, &
         sparse_num_rows, sparse_scale, sparse_trace
    use timer, only: timer_clock

    implicit none

    ! Arguments: input/output
    type(SPAM3), intent(inout)     :: denskern(:)  ! Density kernel

    ! Arguments: input only
    type(SPAM3), intent(in)        :: overlap      ! Overlap matrix
    type(SPAM3), intent(in)        :: inv_overlap  ! Inverse overlap matrix
    integer, intent(in), optional  :: miniter      ! Min num iterations
    logical, intent(in), optional  :: const_ne     ! Flag for constant Ne

    ! Local variables
    integer :: maxiter
    integer :: iter,miniter_local
    integer :: nspins
    integer :: is
    integer :: n_occ(max_spins)
    real(kind=DP), parameter :: minbound = -0.36602540378443864676_DP
    real(kind=DP), parameter :: maxbound = 1.36602540378443864676_DP
    real(kind=DP) :: spin_fac
    real(kind=DP) :: pen
    real(kind=DP) :: rms_err,old_rms_err
    real(kind=DP) :: steplen
    real(kind=DP) :: ne(max_spins)
    real(kind=DP) :: min_occ(max_spins), max_occ(max_spins)
    logical :: local_const_ne
    logical :: converged,temp_workspace
    type(SPAM3) :: ksc,ksks,ksgs,gsgs
    type(SPAM3), allocatable :: dir(:)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering kernel_fix'
#endif

    ! Start timer
    call timer_clock('kernel_fix',1)

    ! cks: intialise variables
    maxiter = pub_maxit_kernel_fix
    ne     =-1.0_DP
    n_occ  = 0.0_DP
    rms_err     =-100.0_DP
    old_rms_err =-100.0_DP


    if (pub_on_root .and. pub_output_detail == VERBOSE) &
         write(stdout,'(/a)') '=================== Penalty &
         &functional idempotency correction =================='

    ! Local copy of miniter
    if (present(miniter)) then
       miniter_local = miniter
    else
       miniter_local = 1
    end if

    ! Local copy of const_ne
    if (present(const_ne)) then
       local_const_ne = const_ne
    else
       local_const_ne = .true.
    end if

    ! Establish number of spins
    nspins = size(denskern)
    spin_fac = 2.0_DP / nspins

    ! cks: determine electron numbers
    call kernel_validate_ks(denskern, overlap)
    if (local_const_ne .and. (pub_kerfix == 2) ) then
       do is=1, nspins
          n_occ(is) =nint( sparse_trace(pub_ks(is)) )
       enddo
    endif


    ! Allocate local workspace
    call internal_workspace_alloc

    ! Set convergence flag
    converged = .false.

    ! Set bounds
    min_occ(:) = 0.0_DP ; max_occ(:) = 1.0_DP

    ! Loop over iterations
    do iter=1,maxiter

       old_rms_err = rms_err

       ! Calculate penalty functional value and (contravariant) gradient
       call internal_penalty_value_grad

       ! cks: Check for convergence only after having done at least
       ! cks: one iteration. This prevents early exit in cases of little
       ! cks: or no denskern truncation which would otherwise hinder
       ! cks: LNV convergence.
       if (iter > 1) then
          ! Check for convergence (i.e. idempotency)
          converged = (pen/sparse_num_rows(denskern(1)) < &
               epsilon(1.0_DP) * 100.0_DP)
       end if

       call comms_bcast(pub_root_node_id,converged)
       if (converged) exit

       ! ndmh: check that rms_err is not stuck (usually indicates kernel
       ! ndmh: with occupation numbers that do not sum to Ne)
       !if ((iter>1).and.(rms_err>old_rms_err*0.9)) then
       !   if (pub_on_root) write(stdout,'(a,a\a)') 'WARNING in kernel_fix:', &
       !        'RMS occupancy error not reducing:', &
       !        'Possibly indicates wrong occupation number'
       !end if
       ! old_rms_err = rms_err

       ! Make sure gradient preserves electron number if required
       if (local_const_ne .and. (pub_kerfix ==1)) call internal_correct_gradient

       ! Take optimal length step in this direction
       call internal_optimal_step(steplen)

       ! Print out details
       if (pub_on_root .and. pub_output_detail == VERBOSE) then
          write(stdout,'(a,i3,a,f9.6,a,f7.4,a,f12.4,a)') &
               ' RMS occupancy error at iteration ',iter,': ', &
               rms_err,' (step ',steplen,' Ne1=',ne(1),')'
          if (nspins == 2) write(stdout,'(a,f12.4)') &
               '                                      &
               &                         Ne2=',ne(2)
       endif


       ! cks: rescale density kernel to correct ne
       if (local_const_ne .and. (pub_kerfix ==2)) &
            call kernel_rescale(denskern,overlap,n_occ,.false.,.true.)


       ! Recalculate K.S
       call kernel_validate_ks(denskern,overlap)

       ! Check bounds for convergence
       call kernel_occupancy_bounds(max_occ,min_occ,denskern,overlap)
       if (minval(min_occ(:)) > minbound .and. &
            maxval(max_occ(:)) < maxbound .and. iter >= miniter_local) exit

    end do

    ! ndmh: evaluate kernel occupancy bounds and kernel RMS occupancy error
    ! ndmh: if this was not already done in the loop
    ! ndmh: 04/06/10 - also evaluate Ne here as this was not done in loop
    if (maxiter==0) then
       call kernel_occupancy_bounds(max_occ,min_occ,denskern,overlap)
       rms_err = kernel_rms_err(denskern,overlap)
       old_rms_err = rms_err
       do is=1,pub_cell%num_spins
          ne(is) = sparse_trace(pub_ks(is))
       end do
    end if

    ! cks: fix occupancies by diagonalisation if it is a hopeless case
    ! cks: and pub_kerfix == 2
    if ( ((rms_err > 0.001_DP) .and. (rms_err > old_rms_err*0.9_DP) ) &
         .and. (pub_kerfix == 2)  )then

       if (pub_on_root .and. pub_output_detail == VERBOSE) &
            write(stdout,'(a)') &
            '               Reset of occupancies applied'

       call kernel_reset_occupancies(denskern, &  ! input-output
            overlap, n_occ )
!CW
       ! cks: denskern has changed: invalidate ks and ksk workspaces
! BUG CW : already done in reset_occupancies above
       !call kernel_workspace_invalidate()
!END CW
    endif



    ! Print results
    ! ndmh: only recalculate pen if displayed value > 0
    if ((.not. converged).and.(rms_err>0.0000005_DP)) &
         call internal_penalty_value_grad
    if (pub_on_root .and. pub_output_detail == VERBOSE) then
         write(stdout,'(a,f9.6,a,f12.4,a)') &
         '            Final RMS occupancy error   : ', &
         rms_err,' (Ne1=',ne(1),')'
         if (nspins == 2) write(stdout,'(a,f12.4)') &
              '                                     &
              &                Ne2=',ne(2)
      endif

    if (pub_on_root .and. pub_output_detail == VERBOSE) then
       if (nspins == 1) then
          write(stdout,'(a,2(f8.4,a))') &
               '               Final occupancy bounds   : [',min_occ(1), &
               ',',max_occ(1),']'
       else
          do is=1,nspins
             write(stdout,'(a,i1,a,2(f7.3,a))') &
                  '            Final occupancy bounds spin ',is,': [', &
                  min_occ(is),',',max_occ(is),']'
          end do
       end if
    end if

    ! Deallocate local workspace
    call internal_workspace_dealloc

    if (pub_on_root .and. pub_output_detail == VERBOSE) &
         write(stdout,'(a/)') '==============================&
         &=================================================='

    ! Stop timer
    call timer_clock('kernel_fix',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving kernel_fix'
#endif

    return

  contains

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine internal_workspace_alloc

      use sparse, only: sparse_create
      use utils, only: utils_alloc_check

      implicit none

      ! Local variable
      integer :: ierr    ! Error flag

      ! Allocate a temporary workspace if we need one
      if(.not.pub_kernel_workspace_allocated)then
          call kernel_workspace_allocate(denskern,overlap)
          temp_workspace=.true.
      else
          temp_workspace=.false.
      endif

      call sparse_create(ksc,pub_ks(1))
      call sparse_create(ksks,pub_ks(1),pub_ks(1))
      call sparse_create(ksgs,ksks)
      call sparse_create(gsgs,ksks)

      allocate(dir(nspins),stat=ierr)
      call utils_alloc_check('internal_workspace_alloc (kernel_fix)','dir',ierr)
      do is=1,nspins
         call sparse_create(dir(is),denskern(is))
      end do

    end subroutine internal_workspace_alloc

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine internal_workspace_dealloc

      use sparse, only: sparse_destroy
      use utils, only: utils_dealloc_check

      implicit none

      ! Local variable
      integer :: ierr   ! Error flag

      ! Deallocate workspace
      do is=nspins,1,-1
         call sparse_destroy(dir(is))
      end do
      deallocate(dir,stat=ierr)
      call utils_dealloc_check('internal_workspace_dealloc (kernel_fix)', &
           'dir',ierr)
      call sparse_destroy(gsgs)
      call sparse_destroy(ksgs)
      call sparse_destroy(ksks)
      call sparse_destroy(ksc)

      ! ndmh: destroy temporary workspace if we created one
      if(temp_workspace)then
         call kernel_workspace_deallocate()
      endif

    end subroutine internal_workspace_dealloc

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !--------------------------------------------------------------------!
    ! The following subroutine calculates the penalty functional value:  !
    !                                                                    !
    !   P = Tr[KSKS(1-KS)(1-KS)]                                         !
    !                                                                    !
    ! and gradient:                                                      !
    !                                                                    !
    !   G = 2(1-KS)(1-2KS)K                                              !
    !                                                                    !
    !--------------------------------------------------------------------!

    subroutine internal_penalty_value_grad

      use sparse, only: sparse_trace, sparse_product, sparse_axpy

      implicit none

      ! Initialise
      ne = 0.0_DP ; pen = 0.0_DP

      ! ndmh: calculate K.S and K.S.K again if they have become invalid
      call kernel_validate_ks(denskern,overlap)
      call kernel_validate_ksk(denskern,overlap)

      do is=1,nspins

         ! Calculate ne = Tr[KS]
         ne(is) = sparse_trace(pub_ks(is))

         ! Calculate ksc := 1 - K.S
         call sparse_copy(ksc,pub_ks(is))
         call sparse_scale(ksc,-1.0_DP,1.0_DP)

         ! Calculate P = Tr[KSKS(1-KS)(1-KS)]
         call sparse_product(ksks,pub_ks(is),ksc)
         pen = pen + sparse_trace(ksks,ksks)

         ! ndmh: a faster way to calculate dir, which also doesn't need to use
         ! ndmh: an extra kskc matrix (saving memory)

         ! ndmh: Calculate G = 2(2.KS.KSK - 3.KSK + K) in dir
         call sparse_product(dir(is),pub_ks(is),pub_ksk(is))
         call sparse_scale(dir(is),2.0_DP)
         call sparse_axpy(dir(is),pub_ksk(is),-3.0_DP)
         call sparse_axpy(dir(is),denskern(is),1.0_DP)
         call sparse_scale(dir(is),2.0_DP)

      end do

      ne = ne * spin_fac
      pen = pen / nspins
      rms_err = sqrt(abs(pen)/sparse_num_rows(denskern(1)))

    end subroutine internal_penalty_value_grad

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !--------------------------------------------------------------------!
    ! The following subroutine ensures that the gradient preserves the   !
    ! correct number of electrons i.e.                                   !
    !                                                                    !
    !   Tr[DS] = 0                                                       !
    !                                                                    !
    !--------------------------------------------------------------------!
    ! Spin polarisation improvements, Chris-Kriton Skylaris, June 2009   !
    !--------------------------------------------------------------------!

    subroutine internal_correct_gradient

      use comms, only: comms_abort
      use sparse, only: sparse_trace, sparse_axpy

      implicit none

      real(kind=DP) :: trds,trsinvs,fac

      ! Calculate Tr(Sinv.S)
      trsinvs = sparse_trace(inv_overlap,overlap)
      if (abs(trsinvs) < epsilon(1.0_DP)*sparse_num_rows(denskern(1))) then
         if (pub_on_root) write(stdout,'(a)') &
              'Error in internal_correct_gradient (kernel_fix): &
              &Tr(Sinv.S) = 0'
         call comms_abort
      end if


      do is=1,nspins

         ! Calculate Tr(D.S)
         trds = sparse_trace(dir(is),overlap)

         fac = -trds / trsinvs

         ! Correct gradient
         call sparse_axpy(dir(is),inv_overlap,fac)

      end do



    end subroutine internal_correct_gradient

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !--------------------------------------------------------------------!
    ! The following subroutine takes an optimal step in the direction G, !
    ! whose length is calculated from the analytic quartic form of the   !
    ! penalty functional                                                 !
    !--------------------------------------------------------------------!

    subroutine internal_optimal_step(steplen)

      use sparse, only: sparse_is_dense, sparse_axpy

      implicit none

      real(kind=DP), intent(out) :: steplen

      real(kind=DP) :: a(4)   ! quartic coefficients

      ! a(1) = 2 Tr[KS(1-KS)(1-2KS)GS]
      ! a(2) = Tr[(1-2KS)(1-2KS)GSGS] - 2 Tr[KSGS(1-KS)GS]
      ! a(3) = -2 Tr[(1-2KS)GSGSGS]
      ! a(4) = Tr[GSGSGSGS]

      ! ndmh: for gradient below threshold, set step to -0.5 if kernel is
      ! ndmh: dense (reliable and much faster than explicit evaluation of
      ! ndmh: actual coefficients in dense cases... may be ok for
      ! ndmh: much sparser kernels and higher thresholds).
      if ((rms_err < 0.000005_DP) .and. sparse_is_dense(denskern(1))) then

          steplen = -0.50000_DP

      ! ndmh: else calculate step explicitly
      else

         ! Calculate coefficients of quartic
         call internal_quartic_coeffs(a)

         ! Calculate optimal step length
         call internal_quartic_min(a,steplen)

      end if

      ! Take optimal step
      do is=1,nspins
         call sparse_axpy(denskern(is),dir(is),steplen)
      end do
      ! ndmh: denskern has changed: invalidate ks and ksk workspaces
      call kernel_workspace_invalidate()

    end subroutine internal_optimal_step

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine internal_quartic_coeffs(a)

      use sparse, only: sparse_trace, sparse_product

      implicit none

      real(kind=DP), intent(out) :: a(4)
      real(kind=DP) :: tr_ksgs, tr_ksksgs, tr_ksgsgs, tr_ksksgsgs, tr_ksksksgs, &
           tr_gsgs, tr_gsgsgs, tr_gsgsgsgs, tr_ksgsgsgs, tr_ksgsksgs


      ! ndmh: New way of finding coefficients involving far fewer
      ! ndmh: sparse_product calls (at cost of 1 extra matrix stored)
      a = 0.0_DP
      do is=1,nspins

         ! Calculate G.S in ksc
         call sparse_product(ksc,dir(is),overlap)

         ! Calculate KS.KS
         call sparse_product(ksks,pub_ks(is),pub_ks(is))

         ! Calculate KS.GS
         call sparse_product(ksgs,pub_ks(is),ksc)

          ! Calculate GS.GS
         call sparse_product(gsgs,ksc,ksc)

         ! Calculate traces involving the above
         tr_ksgs = sparse_trace(ksgs)
         tr_gsgs = sparse_trace(gsgs)
         tr_ksksgs = sparse_trace(pub_ks(is),ksgs)
         tr_ksgsgs = sparse_trace(ksgs,ksc)
         tr_gsgsgs = sparse_trace(gsgs,ksc)
         tr_ksksksgs = sparse_trace(ksks,ksgs)
         tr_ksksgsgs = sparse_trace(ksks,gsgs)
         tr_ksgsksgs = sparse_trace(ksgs,ksgs)
         tr_ksgsgsgs = sparse_trace(ksgs,gsgs)
         tr_gsgsgsgs = sparse_trace(gsgs,gsgs)

         ! a(1) = 2 Tr[KS(1-KS)(1-2KS)GS]
         a(1) = a(1) + 2.0_DP*(tr_ksgs - 3.0_DP*tr_ksksgs + 2.0_DP*tr_ksksksgs)

         ! a(2) = Tr[(1-2KS)(1-2KS)GSGS] - 2 Tr[KSGS(1-KS)GS]
         a(2) = a(2) + tr_gsgs - 4.0_DP*tr_ksgsgs + 4.0_DP*tr_ksksgsgs - &
              2.0_DP*(tr_ksgsgs - tr_ksgsksgs)

         ! a(3) = -2 Tr[(1-2KS)GSGSGS]
         a(3) = a(3) - 2.0_DP*(tr_gsgsgs - 2.0_DP*tr_ksgsgsgs)

         ! a(4) = Tr[GSGSGSGS]
         a(4) = a(4) + tr_gsgsgsgs

      end do

      a = a / nspins

    end subroutine internal_quartic_coeffs

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !--------------------------------------------------------------------!
    ! The following subroutine finds the global minimum of               !
    !   a1*x + a2*x**2 + a3*x**3 + a4*x**4 = 0                           !
    !--------------------------------------------------------------------!

    subroutine internal_quartic_min(a,x)

      use comms, only: comms_abort
      use constants, only: DP, PI, stdout

      implicit none

      ! Arguments
      real(kind=DP), intent(in) :: a(4)
      real(kind=DP), intent(out) :: x

      ! Local variables
      integer :: i
      logical :: foundmin
      real(kind=DP), parameter :: THIRD = 1.0_DP / 3.0_DP
      real(kind=DP) :: aa,bb,cc
      real(kind=DP) :: q,r,theta,rtq,y,z
      real(kind=DP) :: x1(3),q1(3),qmin

      ! If a4 vanishes then the gradient must be zero in which case choose
      ! step length to correspond to purification transformation
      if (abs(a(4)) < epsilon(1.0_DP)) then
         x = -0.5_DP
         return
      end if

      ! Differentiate and write result as  x**3 + a*x**2 + b*x + c = 0
      aa = 0.75_DP * a(3) / a(4)
      bb = 0.5_DP * a(2) / a(4)
      cc = 0.25_DP * a(1) / a(4)

      ! Solve cubic equation
      q = (aa*aa - 3.0_DP*bb)/9.0_DP
      r = (2*aa*aa*aa - 9.0_DP*aa*bb + 27.0_DP*cc)/54.0_DP
      if (r*r < q*q*q) then   ! three real roots
         rtq = sqrt(q)
         theta = acos(r/(rtq*q))
         x1(1) = -2.0_DP * rtq * cos(theta*THIRD) - aa*THIRD
         x1(2) = -2.0_DP * rtq * cos((theta+2.0_DP*PI)*THIRD) - aa*THIRD
         x1(3) = -2.0_DP * rtq * cos((theta-2.0_DP*PI)*THIRD) - aa*THIRD
         foundmin = .false.
         qmin = 0.0_DP ! qoh: Initialise to prevent compiler warning
         do i=1,3
            if (a(2)+3.0_DP*x1(i)*(a(3)+2.0_DP*a(4)*x1(i)) > 0.0_DP) then
               q1(i) = x1(i)*(a(1)+x1(i)*(a(2)+x1(i)*(a(3)+x1(i)*a(4))))
               if (foundmin) then
                  if (q1(i) < qmin) then
                     x = x1(i)
                     qmin = q1(i)
                  end if
               else
                  x = x1(i)
                  qmin = q1(i)
                  foundmin = .true.
               end if
            end if
         end do
         if (.not. foundmin) then
            if (pub_on_root) write(stdout,'(a)') &
                 'Error in internal_quartic_min (kernel_mod.F90): &
                 &no minimum found'
            call comms_abort
         end if
      else   ! only one root -> must be a minimum
         y = -sign(1.0_DP,r)*(abs(r)+sqrt(r*r-q*q*q))**THIRD
         if (abs(y) < epsilon(1.0_DP)) then
            z = 0.0_DP
         else
            z = q / y
         end if
         x = y + z - aa*THIRD
      end if

    end subroutine internal_quartic_min

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  end subroutine kernel_fix

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_normalise(denskern,overlap,inv_overlap,n_occ)

    !=======================================================================!
    ! This subroutine normalises the density kernel so that it satisfies a  !
    ! normalisation constraint, by adding an amount of the contravariant    !
    ! electron number gradient (inverse overlap matrix)                     !
    !-----------------------------------------------------------------------!
    ! Written by Peter Haynes, November 2004                                !
    ! Modified to support spin polarisation, Peter Haynes, July 2006        !
    ! Fixed bug that occured when the number of up and down electrons       !
    ! differed, Chris-Kriton Skylaris, 17 June 2009.                        !
    !=======================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, VERBOSE, stdout, max_spins
    use rundat, only : pub_output_detail
    use sparse, only: SPAM3, sparse_trace, sparse_axpy

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: denskern(:)  ! Density kernel
    type(SPAM3), intent(in)    :: overlap      ! Overlap matrix
    type(SPAM3), intent(in)    :: inv_overlap  ! Inverse overlap matrix
    integer, intent(in)        :: n_occ(max_spins)

    ! Local variables
    integer :: nspins                  ! Number of spins
    integer :: is                      ! Spin counter
    real(kind=DP) :: corrected_ne      ! Desired electron number
    real(kind=DP) :: uncorrected_ne    ! Initial uncorrected electron number
    real(kind=DP) :: rate_of_change_ne ! Rate of change of electron number as
                                       ! gradient is added
    real(kind=DP) :: final_ne          ! Final electron number
    real(kind=DP) :: step              ! Amount by which to correct kernel
    real(kind=DP) :: spin_fac          ! Spin polarisation factor

    ! Establish spin polarisation status
    nspins = size(denskern)
    spin_fac = 2.0_DP / nspins

    ! Calculate rate of change of electron number
    rate_of_change_ne = sparse_trace(inv_overlap,overlap)


    do is=1,nspins

       ! Calculate current electron number
       uncorrected_ne = sparse_trace(denskern(is),overlap)

       if (pub_output_detail == VERBOSE .and. pub_on_root) &
            write(stdout,'(a,i1,a,f16.8)') '   Initial electron number for spin ', is,': ', &
            uncorrected_ne*spin_fac

       ! Do the correction
       corrected_ne = real(n_occ(is), kind=DP)
       step = (corrected_ne - uncorrected_ne) / rate_of_change_ne

       call sparse_axpy(denskern(is),inv_overlap,step)

    enddo


    ! ndmh: mark sk and ks workspaces as invalid
    call kernel_workspace_invalidate()

    ! Check the correction if required
    if (pub_output_detail == VERBOSE) then

       do is=1,nspins
          final_ne = sparse_trace(denskern(is),overlap)

          final_ne = final_ne * spin_fac
          if (pub_on_root) write(stdout,'(a,i1,a,f16.8)') &
               '     Final electron number for spin ', is,': ',final_ne
       enddo

    end if

  end subroutine kernel_normalise


  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_rescale(denskern,overlap,n_occ,can_rescale_ks,silent)

    !=======================================================================!
    ! This subroutine rescales the density kernel so that it satisfies a    !
    ! normalisation constraint                                              !
    !-----------------------------------------------------------------------!
    ! Written by Peter Haynes, November 2004                                !
    ! Spin polarised by Peter Haynes, July 2006                             !
    ! Fixed bug that occured when the number of up and down electrons       !
    ! differed, Chris-Kriton Skylaris, 5 May 2009.                          !
    !=======================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout, max_spins, VERBOSE
    use rundat, only : pub_output_detail
    use sparse, only: SPAM3, sparse_scale, sparse_trace

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: denskern(:)  ! Density kernel
    type(SPAM3), intent(in)    :: overlap      ! Overlap matrix
    integer, intent(in)        :: n_occ(max_spins)
    logical, intent(in),optional :: can_rescale_ks
    logical, intent(in),optional :: silent     ! Whether to supress output

    ! Local variables
    integer :: nspins               ! Number of spins
    integer :: is                   ! Spin counter
    real(kind=DP) :: spin_fac       ! Spin degeneracy factor
    real(kind=DP) :: corrected_ne   ! Desired electron number
    real(kind=DP) :: uncorrected_ne ! Initial uncorrected electron number
    real(kind=DP) :: factor         ! Scaling factor
    real(kind=DP) :: final_ne       ! Final electron number
    logical       :: loc_silent     ! Whether to supress output

    ! ndmh: optional argument
    loc_silent = .false.
    if (present(silent)) then
       loc_silent = silent
    end if

    ! Determine number of spins
    nspins = size(denskern)
    spin_fac = 2.0_DP / nspins

    do is=1,nspins

       corrected_ne = n_occ(is)

       ! Calculate current electron number
       uncorrected_ne = sparse_trace(denskern(is),overlap)

       if (pub_output_detail == VERBOSE .and. pub_on_root .and. (.not.loc_silent)) &
            write(stdout,'(a,i1,a,f16.8)') '    Initial Ne',&
       is,': ', uncorrected_ne*spin_fac

       ! ndmh: protection against divide-by-zero in one spin channel
       if (uncorrected_ne == 0.0_DP) then
          factor = 0.0_DP
       else
          factor = corrected_ne / uncorrected_ne
       end if

       ! Rescale
       call sparse_scale(denskern(is),factor)

       ! ndmh: allow rescaling of ks rather than recalculation
       if(present(can_rescale_ks))then
          if(ks_valid.and.can_rescale_ks)then

             call sparse_scale(pub_ks(is),factor)

             call kernel_workspace_invalidate()
             ks_valid=.true.
          else
             call kernel_workspace_invalidate()
          endif
       else
          call kernel_workspace_invalidate()
       endif

    end do

    ! Check the correction if required
    if (pub_output_detail == VERBOSE .and. (.not.loc_silent)) then

       do is=1,nspins
          final_ne = sparse_trace(denskern(is),overlap)
          if (pub_on_root) write(stdout,'(a,i1,a,f16.8)') &
               '      Final Ne',is,': ',final_ne*spin_fac
       enddo

    end if

  end subroutine kernel_rescale

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_purify(pur_denskern, denskern, overlap, inv_overlap, &
       n_occ, fixed_denskern)

    !=======================================================================!
    ! This subroutine accepts a density kernel and overlap matrix in        !
    ! SPAM3 format and returns a purified density kernel with the given     !
    ! sparsity pattern.  The subroutine is robust in the sense that there   !
    ! are no truncations in matrix multiplications.                         !
    !-----------------------------------------------------------------------!
    ! Written by Peter Haynes on 19 November 2004                           !
    ! Based on earlier versions by Chris-Kriton Skylaris and Peter Haynes   !
    ! Modified by Chris-Kriton Skylaris on 28/03/2005 so that it conserves  !
    ! Ne during kernel_fix.                                                 !
    !-----------------------------------------------------------------------!

#ifdef DEBUG
    use comms, only : pub_on_root
#endif
    use constants, only: stdout, max_spins
    use sparse, only: SPAM3, sparse_copy, sparse_create, sparse_destroy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments:
    type(SPAM3), intent(inout) :: pur_denskern(:) ! Purified density kernel
    type(SPAM3), intent(inout) :: denskern(:)     ! Original density kernel
    type(SPAM3), intent(in) :: overlap            ! Overlap matrix
    type(SPAM3), intent(in) :: inv_overlap        ! Inverse overlap matrix
    integer, intent(in) :: n_occ(max_spins)       ! Number of occupied bands
    logical, intent(in), optional :: fixed_denskern ! Set this to true if
                                                    ! denskern must be fixed

    ! Local variables
    integer :: nspins
    integer :: is
    integer :: ierr
    logical :: temp_workspace
    logical :: temp_denskern_copy
    logical :: loc_fixed_denskern
    type(SPAM3), allocatable :: denskern_backup(:) ! Backup of input denskern
                                                   ! in case it must be altered

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering kernel_purify'
#endif

    ! Start timer
    call timer_clock('kernel_purify',1)

    if (present(fixed_denskern)) then
       loc_fixed_denskern = fixed_denskern
    else
       loc_fixed_denskern = .false.
    end if

    ! Establish number of spins
    nspins = size(denskern)

    ! Allocate temporary workspace if we need one
    if (.not.pub_kernel_workspace_allocated) then
       temp_workspace=.true.
       call kernel_workspace_allocate(denskern,overlap)
    else
       temp_workspace=.false.
    endif
    temp_denskern_copy = .false.

    ! Check that purification will be stable
    call internal_check_denskern

    ! Now do the purification
    call internal_purify

    ! Restore input kernel if we altered it
    if (temp_denskern_copy.and.loc_fixed_denskern) then
       do is=nspins,is,-1
         call sparse_copy(denskern(is),denskern_backup(is))
         call sparse_destroy(denskern_backup(is))
       end do
       deallocate(denskern_backup,stat=ierr)
       call utils_dealloc_check('internal_check_denskern (kernel_purify)', &
            'denskern_backup',ierr)
    end if

    ! Deallocate workspace if we created a temporary one
    if (temp_workspace) then
       call kernel_workspace_deallocate()
    endif

    ! Stop timer
    call timer_clock('kernel_purify',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving kernel_purify'
#endif

  contains

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine internal_check_denskern

      use constants, only: DP, max_spins
      use sparse, only: sparse_extremal_eigenvalue, sparse_copy, sparse_scale

      implicit none

      ! Local variables
      real(kind=DP), parameter :: maxocc = 1.36602540378443864676_DP
      real(kind=DP), parameter :: minocc = -0.36602540378443864676_DP
      real(kind=DP) :: max_occ(max_spins),min_occ(max_spins)

      ! Calculate K.S product
      call kernel_validate_ks(denskern,overlap)

      ! Check bounds
      call kernel_occupancy_bounds(max_occ,min_occ,denskern,overlap)

      ! If the density kernel is stable
      if (maxval(max_occ) < maxocc .and. minval(min_occ) > minocc) then

         return

      else

         if (loc_fixed_denskern) then
            allocate(denskern_backup(nspins),stat=ierr)
            call utils_alloc_check('internal_check_denskern (kernel_purify)', &
                 'denskern_backup',ierr)
            do is=1,nspins
               call sparse_create(denskern_backup(is),denskern(is))
               call sparse_copy(denskern_backup(is),denskern(is))
            end do
            temp_denskern_copy = .true.
         end if

         ! cks: It is massively important to normalise to Ne before the fix
         ! cks: to avoid obtaining wrong occupancies
         call kernel_rescale(denskern, overlap, n_occ, can_rescale_ks=.true.)

         ! Fix the kernel
         call kernel_fix(denskern,overlap,inv_overlap)

         ! Re-calculate product K.S
         call kernel_validate_ks(denskern,overlap)

      end if

    end subroutine internal_check_denskern

    subroutine internal_purify

      use constants, only: DP
      use sparse, only: sparse_copy, sparse_scale, sparse_product

      implicit none

      ! Calculate product K.S.K
      call kernel_validate_ksk(denskern,overlap)

      do is=1,nspins

         ! Put 3I - 2K.S in pub_ks
         call sparse_scale(pub_ks(is),-2.0_DP,3.0_DP)

         ! Calculate product 3K.S.K - 2K.S.K.S.K in pur_denskern
         call sparse_product(pur_denskern(is),pub_ks(is),pub_ksk(is))

      end do

      ! ndmh: pub_ks was broken in calculating this, so invalidate it
      ks_valid = .false.

    end subroutine internal_purify

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  end subroutine kernel_purify

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  real(kind=DP) function kernel_rms_commutator(denskern,ham,overlap)

    !====================================================================!
    ! This function returns the rms value of the commutator between the  !
    ! Hamiltonial matrix and the density kernel: HKS-SKH.                !
    !--------------------------------------------------------------------!
    ! Written by Peter Haynes, November 2004                             !
    ! Based on an original version by Arash Mostofi, July 2002           !
    ! Modifications for faster calculation of commutator by Nicholas     !
    ! Hine in July 2009 and November 2010.                               !
    !====================================================================!

    use constants, only: DP
    use simulation_cell, only : pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_copy, sparse_product, &
         sparse_transpose, sparse_axpy, sparse_destroy, &
         sparse_num_element, sparse_num_rows, sparse_rms_element

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: denskern(:)   ! Density kernel
    type(SPAM3), intent(in) :: ham(:)        ! Hamiltonian
    type(SPAM3), intent(in) :: overlap       ! Overlap matrix

    ! Local variables
    type(SPAM3) :: skh,hks
    type(SPAM3) :: sk,kh
    integer :: nspins
    integer :: is
    real(kind=DP) :: rms

    ! Establish number of spins
    nspins = size(denskern)

    ! ndmh: faster way of calculating commutator

    ! ndmh: construct matrices with HKH structure
    ! ndmh: (a destroyed matrix is fine for structure code use)
    call sparse_create(kh,denskern(1),ham(1))
    call sparse_destroy(kh)
    call sparse_create(skh,ham(1),kh)
    call sparse_create(hks,skh)

    ! ndmh: create temporary for SK
    call sparse_create(sk,overlap,denskern(1))

    ! ndmh: re-evaluate ks if required
    call kernel_validate_ks(denskern,overlap)

    kernel_rms_commutator = 0.0_DP
    do is=1,nspins

       ! ndmh: faster way of calculating commutator
       ! Calculate S.K
       call sparse_transpose(sk,pub_ks(is))

       ! Calculate (S.K).H and transpose to get HKS
       call sparse_product(skh,sk,ham(is))
       call sparse_transpose(hks,skh)

       ! Calculate the commutator SKH-HKS
       call sparse_axpy(skh,hks,-1.0_DP)

       rms = sparse_rms_element(skh)

       kernel_rms_commutator = kernel_rms_commutator + rms * rms

    end do

    kernel_rms_commutator = sqrt(kernel_rms_commutator * &
         sparse_num_element(skh)/real(sparse_num_rows(denskern(1)) &
         *nspins, kind=DP))

    ! Destroy structures
    call sparse_destroy(sk)
    call sparse_destroy(hks)
    call sparse_destroy(skh)

  end function kernel_rms_commutator

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_init_core_ham(denskern, &  ! input-output
       coreham, overlap, n_occ)                ! input

    !======================================================================!
    ! This subroutine initialises the density matrix to the density        !
    ! that is obtained from the orbitals of the core-Hamiltonian matrix    !
    ! in row-indexed sparse matrix format.                                 !
    !----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 3/7/2001 for the ONES code.      !
    ! Adapted for ONETEP by Chris-Kriton Skylaris on 12/1/2004.            !
    ! Modified for parallel SPAM3 by Peter Haynes, July 2006               !
    ! Moved to kernel_mod by Nicholas Hine, August 2008                    !
    ! Reduced memory usage by re-using buffer, Nicholas Hine, October 2009 !
    ! Re-written for general dense matrices by Nicholas Hine, February 2010!
    !======================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use dense, only: DEM, dense_create, dense_convert, dense_destroy, &
         dense_eigensolve, dense_product
    use sparse, only: SPAM3, sparse_num_rows, sparse_scale
    use simulation_cell, only: pub_cell
    use utils, only : utils_alloc_check, utils_dealloc_check

    implicit none


    type(SPAM3), intent(inout) :: denskern
    type(SPAM3), intent(in)    :: coreham
    type(SPAM3), intent(in)    :: overlap
    integer, intent(in)       :: n_occ


    ! cks: internal declarations
    integer :: ierr
    integer :: num
    type(DEM) :: coreham_dens, overlap_dens, denskern_dens
    type(DEM) :: eigenvecs_dens
    real(kind=DP), allocatable :: eigenvalues(:)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering kernel_init_core_ham'
#endif

    ! Do nothing if no states are required
    if (n_occ==0) then
       call sparse_scale(denskern,0.0_DP)
       return
    end if

    ! ndmh: allocate storage for list of eigenvalues
    num = sparse_num_rows(denskern)
    allocate(eigenvalues(num),stat=ierr)
    call utils_alloc_check('kernel_init_core_ham','eigenvalues',ierr)

    ! ndmh: create temporary dense matrices for eigenvecs, Ham and overlap
    call dense_create(eigenvecs_dens,num,num)!n_occ)
    call dense_create(coreham_dens,num,num)
    call dense_create(overlap_dens,num,num)

    ! ndmh: convert sparse Ham and overlap to dense
    call dense_convert(coreham_dens,coreham)
    call dense_convert(overlap_dens,overlap)

    ! ndmh: solve the generalised eigenvalue problem.
    ! ndmh: On return from this subroutine eigenvecs_dens will contain
    !       the eigenvectors of coreham_dens
    call dense_eigensolve(n_occ,eigenvalues,coreham_dens,overlap_dens,1, &
         eigenvecs_dens)

    ! ndmh: remove temporary matrices that are no longer required
    call dense_destroy(overlap_dens)
    call dense_destroy(coreham_dens)

    ! ndmh: create density kernel from eigenvectors by summing
    ! ndmh: first n_occ eigenvectors multiplied by transposed matrix
    call dense_create(denskern_dens,num,num)
    call dense_product(denskern_dens,eigenvecs_dens,eigenvecs_dens, &
         transpose_amat=.false.,transpose_bmat=.true.,first_k=1,last_k=n_occ)

    ! ndmh: convert dense version of denskern back to sparse
    call dense_convert(denskern,denskern_dens)

    ! ndmh: deallocate temporary memory
    call dense_destroy(denskern_dens)
    call dense_destroy(eigenvecs_dens)
    deallocate(eigenvalues,stat=ierr)
    call utils_dealloc_check('kernel_init_core_ham','eigenvalues',ierr)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving kernel_init_core_ham'
#endif

  end subroutine kernel_init_core_ham

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_reset_occupancies(denskern, &  ! input-output
       overlap, n_occ)                ! input

    !====================================================================!
    ! This subroutine manually resets the occupancies of the density     !
    ! kernel using diagonalisation.                                      !
    !--------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in June 2009.                     !
    !====================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, max_spins, stdout
    use dense, only: DEM, dense_create, dense_convert, dense_destroy, &
         dense_eigensolve, dense_product, dense_scale
    use sparse, only: SPAM3, sparse_num_rows
    use simulation_cell, only: pub_cell
    use utils, only : utils_alloc_check, utils_dealloc_check

    implicit none


    type(SPAM3), intent(inout) :: denskern(:)
    type(SPAM3), intent(in)    :: overlap
    integer, intent(in)        :: n_occ(max_spins)

    ! cks: internal declarations
    integer :: ierr
    integer :: num
    integer :: is    ! spin counter
    integer :: nspins
    real(kind=DP), dimension(:), allocatable :: eigenvalues
    type(DEM) :: denskern_dens, overlap_dens, eigs_dens

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &kernel_reset_occupancies'
#endif

    ! Establish number of spins
    nspins = size(denskern)

    num = sparse_num_rows(denskern(1))

    allocate(eigenvalues(num),stat=ierr)
    call utils_alloc_check('kernel_reset_occupancies','eigenvalues',ierr)

    ! ndmh: create temporary dense matrices
    call dense_create(denskern_dens,num,num)
    call dense_create(overlap_dens,num,num)
    call dense_create(eigs_dens,num,num)

    ! cks: loop over spins
    do is=1,nspins

       ! ndmh: expand the density kernel to a dense matrix
       call dense_convert(denskern_dens,denskern(is))

       ! ndmh: multiply by -1 so that we get the highest occupancy eigenvalues
       call dense_scale(denskern_dens,-1.0_DP)

       ! ndmh: expand the overlap matrix to a dense matrix
       call dense_convert(overlap_dens,overlap)

       ! ndmh: solve the generalised eigenvalue problem
       call dense_eigensolve(n_occ(is),eigenvalues,denskern_dens,overlap_dens, &
            2,eigs_dens)

       ! ndmh: construct a density kernel from the first n_occ(is) eigenvectors
       call dense_product(denskern_dens,eigs_dens,eigs_dens, &
            transpose_bmat=.true.,first_k=1,last_k=n_occ(is))

       ! ndmh: convert back to sparse matrix
       call dense_convert(denskern(is),denskern_dens)

    end do

    ! ndmh: clean up temporary matrices
    call dense_destroy(eigs_dens)
    call dense_destroy(overlap_dens)
    call dense_destroy(denskern_dens)

    deallocate(eigenvalues,stat=ierr)
    call utils_dealloc_check('kernel_reset_occupancies','eigenvalues',ierr)

    ! cks: denskern has changed: invalidate ks and ksk workspaces
    call kernel_workspace_invalidate()

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &kernel_reset_occupancies'
#endif

  end subroutine kernel_reset_occupancies

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

  subroutine kernel_from_vecs_asc(density_sq, &
       orbitals,num,n_occ)

    !================================================================!
    ! This subroutine constructs an idempotent density kernel from a !
    ! set of eigenvectors of the density kernel (assumes that the    !
    ! eigenvectors are ordered in ascending order of eigenvalue      !
    ! (occupancy).                                                   !
    !----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in June 2009.                 !
    !================================================================!

    use constants, only: DP, stdout
    use comms, only: comms_abort, pub_on_root

    implicit none

    integer, intent(in) :: num, n_occ
    real(kind=DP), intent(out) :: density_sq(num,num)
    real(kind=DP), intent(in) :: orbitals(num,num)

    ! cks: internal variable declarations
    integer :: row, col, occ_count
    real(kind=DP) :: el

    if (n_occ > num) then
       if (pub_on_root) write(stdout,'(a)') &
            'Error in kernel_construct_from_eigvecs: &
            &number of occupied orbitals exceeds number of NGWFs'
       call comms_abort
    end if

    ! cks: initialise
    density_sq = 0.0_DP

    do col=1,num
       do row=1,num

          ! cks: loop over the last n_occ (highest occupancy) orbitals
          el = 0.0_DP
          do occ_count = num, num-n_occ +1, -1
             el = el + orbitals(row,occ_count) * orbitals(col,occ_count)
          end do

          density_sq(row,col) = el

       end do
    end do


  end subroutine kernel_from_vecs_asc

 subroutine kernel_from_vecs_asc2(density_sq, &
       orbitals,num,n_occ)

    !================================================================!
    ! This subroutine constructs an idempotent density kernel from a !
    ! set of eigenvectors of the density kernel (assumes that the    !
    ! eigenvectors are ordered in ascending order of eigenvalue      !
    ! (occupancy).                                                   !
    !----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in June 2009.                 !
    !================================================================!

    use constants, only: DP, stdout
    use comms, only: comms_abort, pub_on_root

    implicit none

    integer, intent(in) :: num, n_occ
    real(kind=DP), intent(out) :: density_sq(num,num)
    real(kind=DP), intent(in) :: orbitals(num,num)

    ! cks: internal variable declarations
    integer :: row, col, occ_count
    real(kind=DP) :: el

    if (n_occ > num) then
       if (pub_on_root) write(stdout,'(a)') &
            'Error in kernel_construct_from_eigvecs: &
            &number of occupied orbitals exceeds number of NGWFs'
       call comms_abort
    end if

    ! cks: initialise
    density_sq = 0.0_DP

    do col=1,num
       do row=1,num

          ! cks: loop over the last n_occ (highest occupancy) orbitals
          el = 0.0_DP
          do occ_count = 1, n_occ
             el = el + orbitals(row,occ_count) * orbitals(col,occ_count)
          end do

          density_sq(row,col) = el

       end do
    end do


  end subroutine kernel_from_vecs_asc2


  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

    subroutine kernel_christoffel(ngwf_basis, &
         new_ngwfs_on_grid, old_ngwfs_on_grid, new_denskern, & 
         old_denskern, new_inv_overlap, old_inv_overlap, & 
         new_sp_overlap, old_sp_overlap)

      !==============================================================!
      ! This routine updates the density kernel with terms from the  !
      ! Christoffel symbols which preserve the completeness of the   !
      ! basis.                                                       !
      !--------------------------------------------------------------!
      ! Written by David O'Regan in June 2010                        !
      ! Modifications by Nicholas Hine, November 2010.               !
      ! PAW added by David O'Regan 29/3/11                           !
      ! Moved from ngwf_cg_mod and modified accordingly by           ! 
      ! Simon Dubois 20/06/11                                        !
      !==============================================================!

      use augmentation, only: augmentation_overlap
      use function_basis, only: FUNC_BASIS
      use integrals, only: integrals_brappd_ketppd
      use rundat, only: pub_aug
      use simulation_cell, only: pub_cell
      use sparse, only: SPAM3, sparse_create,  sparse_transpose, &
           sparse_product, sparse_destroy, sparse_axpy, sparse_copy, &
           sparse_scale
      use utils, only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Arguments
      type(FUNC_BASIS), intent(in) :: ngwf_basis
      real(kind=DP), intent(in) :: new_ngwfs_on_grid(ngwf_basis%n_ppds * &
           pub_cell%n_pts)
      real(kind=DP), intent(in) :: old_ngwfs_on_grid(ngwf_basis%n_ppds * &
           pub_cell%n_pts)
      type(SPAM3), intent(inout) :: new_denskern(pub_cell%num_spins)
      type(SPAM3), intent(in) :: old_denskern(pub_cell%num_spins)
      type(SPAM3), intent(in) :: new_inv_overlap
      type(SPAM3), intent(in) :: old_inv_overlap
      type(SPAM3), intent(in) :: new_sp_overlap
      type(SPAM3), intent(in) :: old_sp_overlap

      ! Local Variables
      type(SPAM3) :: step, tc_step, step_sp_overlap
      type(SPAM3) :: ktmp
      real(kind=DP), allocatable :: step_on_grid(:)
      logical, parameter :: using_old_inv_overlap = .true.
      integer :: ierr, is

      allocate(step_on_grid(ngwf_basis%n_ppds * pub_cell%n_pts),stat=ierr)
      call utils_alloc_check('kernel_christoffel',&
           &'step_on_grid',ierr)

      ! Create temporary matrix structures
      step%structure = 'S'
      call sparse_create(step)
      call sparse_create(tc_step,step,new_inv_overlap)
      call sparse_create(ktmp,old_denskern(1))

      ! Calculate change |dphi> in NGWFs
      step_on_grid = new_ngwfs_on_grid - old_ngwfs_on_grid

      ! Calculate overlap of the change with the old NGWFs
      call integrals_brappd_ketppd(step, &
           step_on_grid,ngwf_basis,old_ngwfs_on_grid,ngwf_basis)

      ! ndmh: Apply the augmentated overlap operator to this matrix
      if (pub_aug) then

         ! ndmh: <proj_i|dphi_a> = <proj_i|new_phi_a> - <proj_i|old_phi_a>
         call sparse_create(step_sp_overlap,old_sp_overlap)
         call sparse_copy(step_sp_overlap,new_sp_overlap)
         call sparse_axpy(step_sp_overlap,old_sp_overlap,-1.0_DP)

        ! ndmh: Calculate the augmentation of the "step" matrix
         call augmentation_overlap(step,old_sp_overlap,step_sp_overlap)

         ! ndmh: Clean up <proj_i|dphi_a> (no longer needed)
         call sparse_destroy(step_sp_overlap)

      end if

      ! Raise an index in the step-NGWF overlap
      if (using_old_inv_overlap) then
         ! ddor: Ideally, the inverse overlap for the point at
         !       which the gradient is computed is used.
         ! ddor: Here, we have the old NGWFs, but translated
         !       to their new positions.
         call sparse_product(tc_step,step,old_inv_overlap)
      else
         call sparse_product(tc_step,step,new_inv_overlap)
      endif

      do is=1,pub_cell%num_spins
         ! Compute -K <g|phi> S^^
         call sparse_product(new_denskern(is),old_denskern(is),tc_step)
         call sparse_scale(new_denskern(is),-1.0_DP)
         ! Transpose to get -S^^ <phi|g> K
         call sparse_transpose(ktmp,new_denskern(is))
         ! Add this to what we got before
         call sparse_axpy(new_denskern(is),ktmp,1.0_DP)
         ! Add the previous kernel to the calculated change to get new kernel
         call sparse_axpy(new_denskern(is),old_denskern(is),1.0_DP)
      enddo

      ! Destroy temporary matrix structures
      call sparse_destroy(ktmp)
      call sparse_destroy(tc_step)
      call sparse_destroy(step)

      deallocate(step_on_grid,stat=ierr)
      call utils_dealloc_check('kernel_christoffel',&
           &'step_on_grid',ierr)

    end subroutine kernel_christoffel

  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////

    subroutine kernel_basis_transform(ngwf_basis, &
         new_ngwfs_on_grid, old_ngwfs_on_grid, new_denskern, &
         old_denskern, new_inv_overlap, new_sp_overlap, &
         old_sp_overlap)

      !==============================================================!
      ! Transforms the density kernel to its representation in terms !
      ! of new NGWFs                                                 !
      !--------------------------------------------------------------!
      ! Written by Nicholas Hine on 22/10/09                         !
      ! PAW added by David O'Regan 29/3/11                           !
      ! Moved from ngwf_cg_mod and modified accordingly by           ! 
      ! Simon Dubois 20/06/11                                        !
      !==============================================================!

      use augmentation, only: augmentation_overlap
      use function_basis, only: FUNC_BASIS
      use integrals, only: integrals_brappd_ketppd
      use rundat, only: pub_aug
      use simulation_cell, only: pub_cell
      use sparse, only: SPAM3, sparse_create,  sparse_transpose, &
           sparse_product, sparse_destroy, sparse_copy, sparse_axpy

      implicit none

      ! Arguments
      type(FUNC_BASIS), intent(in) :: ngwf_basis
      real(kind=DP), intent(in) :: new_ngwfs_on_grid(ngwf_basis%n_ppds * &
           pub_cell%n_pts)
      real(kind=DP), intent(in) :: old_ngwfs_on_grid(ngwf_basis%n_ppds * &
           pub_cell%n_pts)
      type(SPAM3), intent(in) :: new_inv_overlap
      type(SPAM3), intent(in) :: new_sp_overlap
      type(SPAM3), intent(in) :: old_sp_overlap
      type(SPAM3), intent(inout) :: new_denskern(pub_cell%num_spins)
      type(SPAM3), intent(in) :: old_denskern(pub_cell%num_spins)

      ! Local Variables
      type(SPAM3) :: overlap_old_new
      type(SPAM3) :: sk,ks,ksk
      integer :: is

      ! Create temporary matrix structures
      overlap_old_new%structure = 'S'
      call sparse_create(overlap_old_new)
      sk%structure = 'SK'
      call sparse_create(sk)
      ks%structure = 'KS'
      call sparse_create(ks)
      ksk%structure = 'KSK'
      call sparse_create(ksk)

      ! Calculate overlap of old NGWFs with new NGWFs
      call integrals_brappd_ketppd(overlap_old_new, &
           old_ngwfs_on_grid,ngwf_basis,new_ngwfs_on_grid,ngwf_basis)

      if (pub_aug) then
         ! ddor: Calculate the augmentation of the overlap matrix
         call augmentation_overlap(overlap_old_new,old_sp_overlap, &
              new_sp_overlap)
      end if

      ! Calculate sk_a^b = <f_a|f'_c>.(S^-1)^cb
      call sparse_product(sk,overlap_old_new,new_inv_overlap)
      ! Transpose to find ks^a_b = (S^-1)^ac.<f'_c|f_b>
      call sparse_transpose(ks,sk)

      ! Calculate density kernel transformed to new NGWF basis
      ! K'^ah = (S^-1)^ab.<f'_b|f_g> K^ge <f_e|f'_z>.(S^-1)^zh
      do is=1,pub_cell%num_spins
         call sparse_product(ksk,ks,old_denskern(is))
         call sparse_product(new_denskern(is),ksk,sk)
      end do

      ! Destroy temporary matrix structures
      call sparse_destroy(ksk)
      call sparse_destroy(ks)
      call sparse_destroy(sk)
      call sparse_destroy(overlap_old_new)

    end subroutine kernel_basis_transform


  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////


    subroutine kernel_basis_update(new_denskern, old_denskern, &
         new_overlap, old_overlap, new_inv_overlap)

      !==============================================================!
      ! Update the density kernel such that                          !
      ! K_new*S_new = K_old*S_old                                    !
      !--------------------------------------------------------------!
      ! Written by Simon Dubois on 20/06/11                          !
      !==============================================================!

      use simulation_cell, only: pub_cell
      use sparse, only: SPAM3, sparse_create,  &
           sparse_product, sparse_destroy, sparse_copy, sparse_axpy

      implicit none

      ! Arguments
      type(SPAM3), intent(inout) :: new_denskern(pub_cell%num_spins)
      type(SPAM3), intent(in) :: old_denskern(pub_cell%num_spins)
      type(SPAM3), intent(in) :: new_overlap
      type(SPAM3), intent(in) :: old_overlap
      type(SPAM3), intent(in) :: new_inv_overlap
      
      ! Local Variables
      type(SPAM3) :: ds,ks,ksk
      integer :: is

      ! Create temporary matrix structures
      ds%structure = 'S'
      call sparse_create(ds)
      ks%structure = 'KS'
      call sparse_create(ks)
      ksk%structure = 'KSK'
      call sparse_create(ksk)

      ! Calculate ds
      call sparse_copy(ds,new_overlap)
      call sparse_axpy(ds,old_overlap,-1.0_DP)

      ! Calculate ks, ksk and new_denskern
      do is=1,pub_cell%num_spins
         call sparse_product(ks,old_denskern(is),ds)
         call sparse_product(ksk,ks,new_inv_overlap)
         call sparse_copy(new_denskern(is),old_denskern(is))
         call sparse_axpy(new_denskern(is),ksk,-1.0_DP)
      end do

      ! Destroy temporary matrix structures
      call sparse_destroy(ksk)
      call sparse_destroy(ks)
      call sparse_destroy(ds)

    end subroutine kernel_basis_update


  !/////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////////




end module kernel
