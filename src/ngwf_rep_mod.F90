! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-   !
!==================================================================!
!                                                                  !
!                  NGWF Representation module                      !
!                                                                  !
! This module contains routines associated with the use of the     !
! NGWF_REP and NGWF_HAM types. NGWF_REP is a container to simplify !
! the use of a given NGWF representation, its associated matrices  !
! (overlap, inverse overlap, projector overlap, and the parts of   !
! the Hamiltonian that do not depend on the density). NGWF_HAM is  !
! a container to store the matrices associated with the            !
! Hamiltonian in a given NGWF representation.                      !
!------------------------------------------------------------------!
! This module was created by Nicholas Hine in May 2010.            !
!==================================================================!

module ngwf_representation

   use constants, only: DP, PI, SQRT_PI, NORMAL, VERBOSE, stdout, max_spins
   use sparse, only: SPAM3

   implicit none

   private

   type NGWF_REP

      ! The NGWFs themselves, in PPDs on the standard grid
      real(kind=DP), allocatable     :: ngwfs_on_grid(:)

      ! Overlap matrix (augmented, in PAW or USP formalisms)
      type(SPAM3) :: overlap
      ! Original overlap <phi_a|phi_b> before augmentation
      type(SPAM3) :: ngwf_overlap
      ! Overlap between conduction NGWFs and valence NGWFs
      type(SPAM3) :: cross_overlap, cross_overlap_tr
      ! Inverse overlap matrix (and whether it is initialised)
      type(SPAM3) :: inv_overlap
      logical :: inv_overlap_init
      ! NGWF-projector overlap matrix
      type(SPAM3) :: sp_overlap

      ! NGWF-Hubbard projector overlap matrix
      type(SPAM3) :: hub_overlap, hub_overlap_t

      ! ddor: Hubbard projector - PAW projector overlaps for PAW+U
      type(SPAM3) :: hub_proj_paw_overlap, hub_ngwf_paw_overlap

      ! Occupation number of each spin channel
      integer :: n_occ(max_spins)

      ! Parts of Hamiltonian indep of density kernel
      type(SPAM3) :: kinet
      type(SPAM3) :: nonlocpot(1) ! only for nonlocal NCPPs

      ! Postfix string for the structure code of matrices of this rep
      character(len=9) :: postfix

   end type NGWF_REP

   type NGWF_HAM

      ! Components of the Hamiltonian dependent on the density kernel
      type(SPAM3), allocatable :: nonlocpot(:) ! only in PAW/USP
      type(SPAM3), allocatable :: hubbard_ham(:)
      type(SPAM3), allocatable :: hfexchange(:)
      type(SPAM3), allocatable :: lhxc(:)
      type(SPAM3), allocatable :: ham(:)
      type(SPAM3), allocatable :: cond_non_proj_ham(:) ! only in conduction
      type(SPAM3), allocatable :: dijhat(:)
      real(kind=DP) :: cond_shift ! only used in conduction

      ! jd: Metric matrix for HFx
      type(SPAM3) :: full_vmatrix

   end type NGWF_HAM

   public :: NGWF_REP
   public :: NGWF_HAM

   public :: ngwf_rep_create
   public :: ngwf_rep_destroy
   public :: ngwf_ham_create
   public :: ngwf_ham_destroy

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_rep_create(rep,postfix,elements)

    !=========================================================================!
    ! This subroutine allocates memory for the arrays in the NGWF_REP type,   !
    ! apart from the density-dependent matrices of the Hamiltonian (which are !
    ! done separately to save memory).                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rep (inout) : Container for NGWFs and matrices to be allocated.       !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 07/10/2010                                  !
    !=========================================================================!

    use comms, only: pub_on_root
    use ion, only: ELEMENT
    use rundat, only: pub_any_nl_proj, pub_paw, pub_aug, pub_output_detail, &
         pub_hubbard
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_create, sparse_destroy
    use sparse_initialise, only: sparse_init_rep
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(NGWF_REP), intent(inout) :: rep
    type(ELEMENT), intent(in) :: elements(:)     ! List of atoms in cell
    character(*), intent(in) :: postfix

    ! Set postfix, if present
    rep%postfix = postfix

    ! Initialise matrices for this representation
    call sparse_init_rep(elements,rep%postfix)

    ! Banner before allocating all these matrices so it is clear where
    ! it goes wrong if it runs out of memory
    if ((pub_on_root).and.(pub_output_detail>=NORMAL)) then
       if (postfix=='c') then
          if (pub_on_root) write(stdout,'(a)',advance='no') &
               'Conduction NGWF sparse matrix initialisation ...'
       else if (postfix=='j') then
          if (pub_on_root) write(stdout,'(a)',advance='no') &
               'Joint Valence + Conduction NGWF sparse matrix initialisation ...'
       else if (postfix=='a') then
          if (pub_on_root) write(stdout,'(a)',advance='no') &
               'Auxiliary NGWF sparse matrix initialisation ...'
       else
          if (pub_on_root) write(stdout,'(a)',advance='no') &
               'Sparse matrix initialisation ...'
       end if
    end if

    ! Create the matrix structure for overlap
    rep%overlap%structure = 'S'//rep%postfix
    call sparse_create(rep%overlap)

    ! Create the matrix structure for direct NGWF-only overlap
    ! We do not actually need the matrix to be allocated though, as
    ! it will only be used for its structure code and library entry
    if (postfix/='a') then
       if (pub_aug) then
          rep%ngwf_overlap%structure = 'O'//rep%postfix
       else
          rep%ngwf_overlap%structure = 'S'//rep%postfix
       end if
    else
       ! ndmh: for auxiliary NGWFs, we only care about diagonal blocks
       ! ndmh: in kinet, lhxc etc
       if (pub_aug) then
          rep%ngwf_overlap%structure = 'O'//rep%postfix
       else
          rep%ngwf_overlap%structure = 'S'//rep%postfix
       end if
    end if
    call sparse_create(rep%ngwf_overlap)
    call sparse_destroy(rep%ngwf_overlap)

    ! Create the matrix structure for inverse overlap
    rep%inv_overlap%structure = 'K'//rep%postfix
    call sparse_create(rep%inv_overlap)

    ! Create the NGWF-projector overlap matrix if required
    if (pub_any_nl_proj.or.pub_paw) then
       rep%sp_overlap%structure = 'Q'//rep%postfix
       call sparse_create(rep%sp_overlap)
    end if

    ! Create the NGWF-projector overlap matrix if required
    if (pub_hubbard) then
       rep%hub_overlap%structure = 'V'//rep%postfix
       call sparse_create(rep%hub_overlap)
       rep%hub_overlap_t%structure = 'W'//rep%postfix
       call sparse_create(rep%hub_overlap_t)
       if (pub_aug) then ! ddor: PAW+U
          rep%hub_proj_paw_overlap%structure = 'Y'//rep%postfix
          call sparse_create(rep%hub_proj_paw_overlap)
          rep%hub_ngwf_paw_overlap%structure = 'Q'//rep%postfix
          call sparse_create(rep%hub_ngwf_paw_overlap)
       endif
    end if

    ! If we have norm-conserving PPs and nonlocal projectors, create
    ! nonlocal potential matrix
    if (pub_any_nl_proj.and.(.not.pub_aug)) then
       rep%nonlocpot%structure = 'H'//rep%postfix
       call sparse_create(rep%nonlocpot(1))
    end if

    ! Create the kinetic part of the hamiltonian
    call sparse_create(rep%kinet,rep%ngwf_overlap)

    ! ndmh: ensure inverse overlap is (re-)initialised for this rep
    rep%inv_overlap_init = .false.

    ! lr408: If this is a conduction representation, create the cross overlap 
    ! lr408: matrix
    if ((rep%postfix=='c').or.(rep%postfix=='j').or.(rep%postfix=='a')) then
       rep%cross_overlap%structure = 'T'//rep%postfix
       call sparse_create(rep%cross_overlap)
    end if
    if ((rep%postfix=='a')) then
       rep%cross_overlap_tr%structure = 'U'//rep%postfix
       call sparse_create(rep%cross_overlap_tr)
    end if

  end subroutine ngwf_rep_create


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_rep_destroy(rep)

    !=========================================================================!
    ! This subroutine allocates memory for the arrays in the NGWF_REP type,   !
    ! apart from the density-dependent matrices of the Hamiltonian (which are !
    ! done separately to save memory).                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   rep (inout) : Container for NGWFs and matrices to be deallocated.     !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 07/10/2010                                  !
    !=========================================================================!

    use rundat, only: pub_any_nl_proj, pub_aug, pub_paw, pub_hubbard
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_destroy
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(NGWF_REP), intent(inout) :: rep

    if ((rep%postfix=='c').or.(rep%postfix=='j')) then
       call sparse_destroy(rep%cross_overlap)
    end if

    ! cks: memory deallocation for sparse matrices
    if (pub_hubbard) then
       if (pub_aug) then ! ddor: PAW+U
          call sparse_destroy(rep%hub_ngwf_paw_overlap)
          call sparse_destroy(rep%hub_proj_paw_overlap)
       endif
       call sparse_destroy(rep%hub_overlap_t)
       call sparse_destroy(rep%hub_overlap)
    end if
    if (pub_any_nl_proj.or.pub_paw) then
       call sparse_destroy(rep%sp_overlap)
    end if
    if (pub_any_nl_proj.and.(.not.pub_aug)) then
       call sparse_destroy(rep%nonlocpot(1))
    end if

    rep%inv_overlap_init = .false.
    call sparse_destroy(rep%inv_overlap)
    call sparse_destroy(rep%kinet)
    call sparse_destroy(rep%overlap)

  end subroutine ngwf_rep_destroy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_ham_create(ham,rep,lhxc_nl_only)

    !=========================================================================!
    ! This subroutine allocates memory for the arrays in the NGWF_HAM type    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   ham (inout) : Container for Hamiltonian matrices to be allocated.     !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 07/10/2010                                  !
    !=========================================================================!

    use rundat, only: pub_any_nl_proj, pub_aug, pub_hubbard, pub_usehfx, &
         cond_init_shift
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_create
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(inout) :: ham
    type(NGWF_REP), intent(in) :: rep
    logical, intent(in), optional :: lhxc_nl_only

    ! Local Variables
    integer :: is
    integer :: ierr
    logical :: loc_usehfx, loc_hubbard
    logical :: loc_cond

    ! Deal with optional argument (used to suppress hfx and hubbard parts)
    if (present(lhxc_nl_only)) then
       if (lhxc_nl_only) then
          loc_usehfx = .false.
          loc_hubbard = .false.
       else
          loc_usehfx = pub_usehfx
          loc_hubbard = pub_hubbard
       end if
    else
       loc_usehfx = pub_usehfx
       loc_hubbard = pub_hubbard
    end if

    if (rep%postfix=='c') then
       loc_cond = .true.
    else
       loc_cond = .false.
    end if

    ! Allocate matrices for local potential and final Hamiltonian
    allocate(ham%ham(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('ngwf_ham_create','ham',ierr)
    allocate(ham%lhxc(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('ngwf_ham_create','lhxc',ierr)

    ! ndmh: in PAW, nonlocpot is density-dependent
    if (pub_aug) then
       allocate(ham%nonlocpot(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','nonlocpot',ierr)
       allocate(ham%dijhat(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','dijhat',ierr)
    else
       allocate(ham%nonlocpot(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','nonlocpot',ierr)
       allocate(ham%dijhat(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','dijhat',ierr)
    end if

    ! qoh: Allocate Hartree-Fock exchange sparse matrix if necessary
    if (loc_usehfx) then
       allocate(ham%hfexchange(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','hfexchange',ierr)
    else
       allocate(ham%hfexchange(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','hfexchange',ierr)
    end if

    ! ddor: Allocate DFT+U Hamiltonian if necessary
    if (loc_hubbard) then
       allocate(ham%hubbard_ham(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','hubbard_ham',ierr)
    else
       allocate(ham%hubbard_ham(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','hubbard_ham',ierr)
    endif

    ! lr408: Allocate conduction projected Hamiltonian if necessary -
    ! lr408: only needed if this is for the cond rep
    if (loc_cond) then
       allocate(ham%cond_non_proj_ham(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','cond_non_proj_ham',ierr)
       ham%cond_shift = cond_init_shift
    else
       allocate(ham%cond_non_proj_ham(0),stat=ierr)
       call utils_alloc_check('ngwf_ham_create','cond_non_proj_ham',ierr)
    end if


    do is=1,pub_cell%num_spins
       if (loc_cond) then
          ham%ham(is)%structure = 'L'//rep%postfix
       else
          ham%ham(is)%structure = 'H'//rep%postfix
       end if
       call sparse_create(ham%ham(is))
       call sparse_create(ham%lhxc(is),rep%ngwf_overlap)
       if (pub_aug) then
          ham%nonlocpot%structure = 'H'//rep%postfix
          call sparse_create(ham%nonlocpot(is))
          ham%dijhat%structure = 'E'
          call sparse_create(ham%dijhat(is))
       end if
       if (loc_usehfx) then
          ham%hfexchange(is)%structure = 'K'//rep%postfix
          call sparse_create(ham%hfexchange(is))
       end if
       if (loc_hubbard) then
          ham%hubbard_ham(is)%structure = 'H'//rep%postfix
          call sparse_create(ham%hubbard_ham(is))
       end if
       if (loc_cond) then
          ham%cond_non_proj_ham(is)%structure = 'H'//rep%postfix
          call sparse_create(ham%cond_non_proj_ham(is))
       end if
    end do

    if (pub_hubbard) then


    end if

  end subroutine ngwf_ham_create


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_ham_destroy(ham)

    !=========================================================================!
    ! This subroutine allocates memory for the arrays in the NGWF_REP type,   !
    ! apart from the density-dependent matrices of the Hamiltonian (which are !
    ! done separately to save memory).                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   ham (inout) : Container for Hamiltonian matrices to be deallocated.   !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 07/10/2010                                  !
    !=========================================================================!

    use rundat, only: pub_aug, pub_hubbard, pub_usehfx
    use simulation_cell, only: pub_cell
    use sparse, only: sparse_destroy
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(NGWF_HAM), intent(inout) :: ham

    ! Local Variables
    integer :: is
    integer :: ierr

    ! Deallocate storage in matrices

    do is=pub_cell%num_spins,1,-1
       if (size(ham%cond_non_proj_ham)>0) &
          call sparse_destroy(ham%cond_non_proj_ham(is))
       if (pub_hubbard.and.size(ham%hubbard_ham)>0) then
          call sparse_destroy(ham%hubbard_ham(is))
       end if
       if (pub_usehfx.and.size(ham%hfexchange)>0) &
          call sparse_destroy(ham%hfexchange(is))
       if (pub_aug) then
          call sparse_destroy(ham%dijhat(is))
          call sparse_destroy(ham%nonlocpot(is))
       end if
       call sparse_destroy(ham%lhxc(is))
       call sparse_destroy(ham%ham(is))
    end do

    ! Deallocate all matrix arrays
    if (allocated(ham%cond_non_proj_ham)) then
       deallocate(ham%cond_non_proj_ham,stat=ierr)
       call utils_dealloc_check('ngwf_ham_create','cond_non_proj_ham',ierr)
    end if
    if (allocated(ham%hubbard_ham)) then
       deallocate(ham%hubbard_ham,stat=ierr)
       call utils_dealloc_check('ngwf_ham_create','hubbard_ham',ierr)
    end if
    if (allocated(ham%hfexchange)) then
       deallocate(ham%hfexchange,stat=ierr)
       call utils_dealloc_check('ngwf_ham_create','hfexchange',ierr)
    end if
    if (allocated(ham%nonlocpot)) then
       deallocate(ham%nonlocpot,stat=ierr)
       call utils_dealloc_check('ngwf_ham_create','nonlocpot',ierr)
    end if
    if (allocated(ham%dijhat)) then
       deallocate(ham%dijhat,stat=ierr)
       call utils_dealloc_check('ngwf_ham_create','dijhat',ierr)
    end if
    if (allocated(ham%lhxc)) then
       deallocate(ham%lhxc,stat=ierr)
       call utils_dealloc_check('ngwf_ham_create','lhxc',ierr)
    end if
    if (allocated(ham%ham)) then
       deallocate(ham%ham,stat=ierr)
       call utils_dealloc_check('ngwf_ham_create','ham',ierr)
    end if

  end subroutine ngwf_ham_destroy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module ngwf_representation
