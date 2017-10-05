! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=================================================================!
!                                                                 !
!                    Bandstructure Module                         !
!                                                                 !
! This module calculates the bandstructure using the NGWFs as a   !
! basis.                                                          !
!-----------------------------------------------------------------!
! Written by Peter Haynes on 10/10/2006 - "k.p" style             !
! Added "tight-binding" style on 16/1/2010                        !
! This file created by Nicholas Hine 09/03/2010 by moving routine !
! from properties_mod.                                            !
! Modified to use NGWF_REP and NGWF_HAM by Nicholas Hine in       !
! October 2010.                                                   !
!=================================================================!

module bandstructure

  implicit none

  private

  public :: bandstructure_calculate

contains

!CW: added ,hub argument
  subroutine bandstructure_calculate(ham, rep, ngwf_basis, proj_basis, elements, &
       ham_type,hub,hub_proj_basis)
!END CW

    !======================================================================!
    ! This subroutine calculates the bandstructure using the NGWFs as a    !
    ! basis.                                                               !
    !----------------------------------------------------------------------!
    ! Written by Peter Haynes on 10/10/2006 - "k.p" style                  !
    ! Added "tight-binding" style on 16/1/2010                             !
    !======================================================================!

    use comms, only: pub_on_root, pub_total_num_nodes, pub_my_node_id, &
         comms_abort, comms_bcast, comms_barrier
    use constants, only: DP, TWO_PI, max_spins, stdout, stderr
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(.DOT.), operator(-), operator(*), &
         geometry_magnitude, operator(+)
    use integrals, only: integrals_grad
    use ion, only: element
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use parallel_strategy, only: pub_num_atoms_on_node, &
         pub_first_atom_on_node, pub_orig_atom
    use projectors, only: projectors_func_ovlp_box, projectors_deallocate_set
    use pseudopotentials, only: pseudopotentials_nonlocal_mat, &
         pseudo_species_init_proj, nlps_projectors
    use rundat, only: pub_any_nl_proj, pub_bs_unfold, pub_bs_kpoint_path_end, &
         pub_bs_num_eigenvalues, pub_bs_kpoint_path_spacing, pub_bs_method, &
         pub_bs_kpoint_path_length, pub_bs_kpoint_path_start, pub_rootname, &
         pub_cond_calculate
!CW
    use rundat, only : pub_dmft_nkpoints
    use hubbard_build, only: HUBBARD_MODEL
    use sparse, only: sparse_convert,sparse_num_cols,sparse_num_rows,sparse_get_element,sparse_put_element
!END CW
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_copy, &
         sparse_axpy, sparse_convert, sparse_get_block, sparse_put_block, &
         sparse_index_length, sparse_generate_index
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit

    implicit none

    ! Arguments
!CW
     TYPE(HUBBARD_MODEL),optional      :: hub
     type(FUNC_BASIS), optional        :: hub_proj_basis
!END CW
    type(NGWF_HAM), intent(in) :: ham              ! Hamiltonian matrices
    type(NGWF_REP), intent(in) :: rep              ! NGWF Representation
    type(FUNC_BASIS), intent(in) :: ngwf_basis     ! NGWF basis type
    type(FUNC_BASIS), intent(in) :: proj_basis     ! Projector basis type
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    ! lr408: optional string denoting which basis for conduction calculations
    character(len=*), intent(in) :: ham_type 

    ! Local variables
    type(SPAM3), allocatable :: ham0(:)   ! k-independent part of Hamiltonian
    type(SPAM3), allocatable :: nl_k(:)   ! Nonlocal matrix at k-point
    type(SPAM3), allocatable :: ham_ks(:) ! Hamiltonian for k-point and spin
    type(SPAM3), allocatable :: over_ks(:) ! Overlap for k-point and spin
!CW
    type(SPAM3), allocatable :: hubT_ks(:) ! projectors Pk
!END CW
    type(SPAM3) :: grad(3)                ! Grad matrix
    type(POINT) :: kstart                 ! k-point at start of path segment
    type(POINT) :: kend                   ! k-point at end of path segment
    type(POINT) :: kseg                   ! path segment
    type(POINT) :: kptdiff                ! difference between k-points
    type(POINT) :: atdisp                 ! atom displacement
    type(POINT), allocatable :: kpt(:)    ! k-point
    integer :: ierr                       ! Error flag
    integer :: nspins                     ! Number of spins
    integer :: is                         ! Spin counter
    integer :: dim                        ! Dimension counter
    integer :: ipath                      ! Bandstructure path counter
    integer :: nsegpts                    ! Number of points along path segment
    integer :: isegpt                     ! Path segment point counter
    integer :: nkpts                      ! Total number of kpts
    integer :: ikpt,jkpt                  ! k-point counter
    integer :: ingwf,jngwf                ! NGWF counters
    integer :: max_ngwfs_on_atom          ! Max number of NGWFs on any atom
    integer :: nlower,nupper              ! Limits on NGWF bounds for output
    integer :: nkpts_red                  ! Reduced list of unique k-points
    integer :: node                       ! Node counter
    integer :: lwork                      ! Workspace length
    integer :: iat,jat                    ! Atom counters
    integer :: orig_iat, orig_jat         ! Atom counters in input file order
    integer :: loc_iat                    ! Atom counter on local node
    integer :: jdx                        ! Sparse matrix counters
    integer :: ispec                      ! Species of atom iat
    integer :: wrapped                    ! Flag for wrapping
    integer :: bands_unit                 ! I/O unit for .bands file
    integer :: agr_unit                   ! I/O unit for .agr file
    integer :: hamidxlen,overidxlen       ! Sparse index lengths
    integer :: max_degeneracy             ! Maximum degeneracy
    integer, allocatable :: uni_entry(:)  ! Pointer to unique entry
    integer, allocatable :: translate(:,:)    ! Translations for unfolding
    integer, allocatable :: hamidx(:)     ! Sparse Hamiltonian index
    integer, allocatable :: overidx(:)    ! Sparse overlap index
    logical :: unique                     ! Flag for reducing k-points
    logical :: unfold                     ! Flag for unfolding Brillouin zone
    logical :: found                      ! Flag for finding equivalent atom
    logical, allocatable :: flag(:)       ! Flag when finding equivalent atoms
    real(kind=DP), parameter :: tol = 1.0e-8_DP
    real(kind=DP), parameter :: recip_twopi = 1.0_DP / TWO_PI
    real(kind=DP) :: seglen               ! Path segment length
    real(kind=DP) :: pathlen              ! Total path length
    real(kind=DP) :: distsq               ! Squared distance
    real(kind=DP) :: phase                ! Translation phase
    real(kind=DP) :: kcart(3)             ! k-point in Cartesians
    real(kind=DP) :: kfrac(3)             ! k-point in fractional coordinates
    real(kind=DP) :: trvec(3)             ! Translation vector in real-space
    real(kind=DP) :: trpos(3)             ! Translated atom position
    real(kind=DP) :: dispfrac(3)          ! Displacement vector (fractional)
    real(kind=DP) :: dispabs(3)           ! Displacement vector (absolute)
    real(kind=DP) :: efermi               ! Fermi energies
    real(kind=DP), allocatable :: ov_sq(:,:)  ! Overlap in dense format
    real(kind=DP), allocatable :: dwork(:,:)  ! Diagonalization workspace
    real(kind=DP), allocatable :: en_ks(:)    ! Eigenvalues for k-point/spin
    real(kind=DP), allocatable :: kpt_ks(:,:) ! Unfolded k-point for k-point/sp
    real(kind=DP), allocatable :: eig(:,:,:)  ! All eigenvalues
    real(kind=DP), allocatable :: kunfold(:,:,:,:)  ! Unfolded k-point
    real(kind=DP), allocatable :: posfrac(:,:)   ! Fractional atomic positions
    complex(kind=DP) :: i2pi                     ! 2i*pi
    complex(kind=DP) :: zphase                   ! Complex translation phase
    complex(kind=DP), allocatable :: hm_sq(:,:)  ! Hamiltonian in dense format
    complex(kind=DP), allocatable :: zov_sq(:,:) ! Overlap in dense format
    complex(kind=DP), allocatable :: bmat(:,:)   ! Diagonalization workspace
!CW
    complex(8),allocatable        :: bmat_rec(:,:),bmatT_rec(:,:)
    integer(4)                    :: hub_proj,nproj,ii1,ii2,ngwfs_on_atom
    complex(8)                    :: cblock
!END CW
    complex(kind=DP), allocatable :: zwork(:)    ! Diagonalization workspace
    complex(kind=DP), allocatable :: trwfn(:)    ! Translated wavefunction
    complex(kind=DP), allocatable :: block(:,:)  ! Sparse matrix block
    complex(kind=DP), allocatable :: psitrpsi(:,:)  ! Inner products
    complex(kind=DP), allocatable :: zlambda(:)  ! Complex eigenvalues
    character(len=90) :: filename                ! Output filename
!CW
    integer  ::  i1,i2,i3,iproj,ingwfs,jp
    real(8)  ::  t1,t2,t3
!END CW

    ! ndmh: SPAM3 approach to nonlocal pseudopotentials
    type(SPAM3) :: sp_overlap
    type(SPAM3) :: sp_overlap_rc

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering bandstructure_calculate'
#endif

    !pdh: define constant to avoid problem with gfortran
    i2pi = cmplx(0.0_DP,-TWO_PI,kind=DP)

    ! Check details of unfolding
    unfold = (pub_bs_unfold(1) > 0 .and. pub_bs_unfold(2) > 0 .and. &
         pub_bs_unfold(3) > 0)

!CW
    if(abs(pub_dmft_nkpoints)>1)then
      unfold=.false.
    endif
!END CW

    if (unfold) pub_bs_unfold = max(pub_bs_unfold,1)

    ! Count the number of k-points needed
    nkpts = 0
    pathlen = 0.0_DP
    do ipath=1,pub_bs_kpoint_path_length
       kstart = pub_bs_kpoint_path_start(1,ipath) * pub_cell%b1 + &
            pub_bs_kpoint_path_start(2,ipath) * pub_cell%b2 + &
            pub_bs_kpoint_path_start(3,ipath) * pub_cell%b3
       kend = pub_bs_kpoint_path_end(1,ipath) * pub_cell%b1 + &
            pub_bs_kpoint_path_end(2,ipath) * pub_cell%b2 + &
            pub_bs_kpoint_path_end(3,ipath) * pub_cell%b3
       kseg = kend - kstart
       seglen = geometry_magnitude(kseg)
       pathlen = pathlen + seglen
       nsegpts = max(int(seglen / pub_bs_kpoint_path_spacing),1)
       nkpts = nkpts + nsegpts + 1
    end do
!CW
    if(abs(pub_dmft_nkpoints)>1)then
     jp=abs(pub_dmft_nkpoints)

     if(pub_dmft_nkpoints<0)then
      nkpts=0
      do i1=-(jp-1)/2,(jp-1)/2
       do i2=-(jp-1)/2,(jp-1)/2
        do i3=-(jp-1)/2,(jp-1)/2
         if(  i3>0  .or.  (i3==0  .and.  ( (i2>0).or.(i2==0.and.i1>=0)) ) )then
          nkpts=nkpts+1
         endif
        enddo
       enddo
      enddo
      write(*,*) 'THERE ARE [X] K POINTS WITH THE INVERSION SYMMETRY : ', nkpts
      allocate(kpt(nkpts))
      ikpt=0
      !inversion symmetry only
      do i1=-(jp-1)/2,(jp-1)/2
       do i2=-(jp-1)/2,(jp-1)/2
        do i3=-(jp-1)/2,(jp-1)/2
         if(  i3>0  .or.  (i3==0  .and.  ( (i2>0).or.(i2==0.and.i1>=0)  ) ) )then
          ikpt=ikpt+1 
          t1=dble(i1) /  dble(jp) 
          t2=dble(i2) /  dble(jp) 
          t3=dble(i3) /  dble(jp) 
          kpt(ikpt)=t1*pub_cell%b1+t2*pub_cell%b2+t3*pub_cell%b3
         endif
        enddo
       enddo
      enddo
     else
       nkpts=jp**3
       allocate(kpt(nkpts),stat=ierr)
       ikpt=0
       do i1=1,jp
        do i2=1,jp
         do i3=1,jp
          ikpt=ikpt+1
          t1=dble(i1-1)  /  dble(2*jp-1)
          t2=dble(i2-1)  /  dble(2*jp-1)
          t3=dble(i3-1)  /  dble(2*jp-1)
          kpt(ikpt)=t1*pub_cell%b1+t2*pub_cell%b2+t3*pub_cell%b3
         enddo
        enddo
       enddo
     endif
     goto 40
    endif
!END CW

    ! Compile a complete list of k-points
    allocate(kpt(nkpts),stat=ierr)
    call utils_alloc_check('bandstructure_calculate','kpt',ierr)
    ikpt = 0
    do ipath=1,pub_bs_kpoint_path_length
       kstart = pub_bs_kpoint_path_start(1,ipath) * pub_cell%b1 + &
            pub_bs_kpoint_path_start(2,ipath) * pub_cell%b2 + &
            pub_bs_kpoint_path_start(3,ipath) * pub_cell%b3
       kend = pub_bs_kpoint_path_end(1,ipath) * pub_cell%b1 + &
            pub_bs_kpoint_path_end(2,ipath) * pub_cell%b2 + &
            pub_bs_kpoint_path_end(3,ipath) * pub_cell%b3
       kseg = kend - kstart
       seglen = geometry_magnitude(kseg)
       nsegpts = max(int(seglen / pub_bs_kpoint_path_spacing),1)
       do isegpt=0,nsegpts
          ikpt = ikpt + 1
          kpt(ikpt) = kstart + (isegpt / real(nsegpts,kind=DP)) * kseg
       end do
    end do
!CW
    40 continue
!ENDCW
    ! Now compile a reduced list of unique k-points avoiding repetition
    ! (e.g. at start and end of consecutive segments)
    allocate(uni_entry(nkpts),stat=ierr)
    call utils_alloc_check('bandstructure_calculate','uni_entry',ierr)
    nkpts_red = 0
    do ikpt=1,nkpts
       unique = .true.
       do jkpt=1,ikpt-1
          kptdiff = kpt(ikpt) - kpt(jkpt)
          if (geometry_magnitude(kptdiff) < 1.0e-6_DP) then
             unique = .false.
             uni_entry(ikpt) = uni_entry(jkpt)
             exit
          end if
       end do
       if (unique) then
          nkpts_red = nkpts_red + 1
          uni_entry(ikpt) = nkpts_red
!CW
          write(*,*) 'K POINTS : ', uni_entry(ikpt) 
!END CW
       end if
    end do

    ! Allocate workspace for diagonalisation and unfolding
    nspins = pub_cell%num_spins
    allocate(eig(ngwf_basis%num,nkpts,nspins),stat=ierr)
    call utils_alloc_check('bandstructure_calculate','eig',ierr)
    allocate(en_ks(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('bandstructure_calculate','en_ks',ierr)
    allocate(hm_sq(ngwf_basis%num,ngwf_basis%num),stat=ierr)
    call utils_alloc_check('bandstructure_calculate','hm_sq',ierr)
    allocate(dwork(ngwf_basis%num,3),stat=ierr)
    call utils_alloc_check('bandstructure_calculate','dwork',ierr)
    allocate(zwork(ngwf_basis%num*2),stat=ierr)
    call utils_alloc_check('bandstructure_calculate','zwork',ierr)
    allocate(bmat(ngwf_basis%num,ngwf_basis%num),stat=ierr)
    call utils_alloc_check('bandstructure_calculate','bmat',ierr)
!CW
    if(abs(pub_dmft_nkpoints)>1)then
       iproj  = NINT(sparse_num_rows(rep%hub_overlap_t))
       ingwfs = NINT(sparse_num_cols(rep%hub_overlap_t))
       allocate(bmatT_rec(iproj,ingwfs))
       allocate(bmat_rec(ingwfs,iproj))
    endif
!END CW   
    allocate(zov_sq(ngwf_basis%num,ngwf_basis%num),stat=ierr)
    call utils_alloc_check('bandstructure_calculate','zov_sq',ierr)
    if (unfold) then
       allocate(kpt_ks(ngwf_basis%num,3),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','kpt_ks',ierr)
       allocate(translate(pub_cell%nat,3),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','translate',ierr)
       allocate(kunfold(ngwf_basis%num,3,nkpts_red,nspins),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','kunfold',ierr)
       allocate(flag(pub_cell%nat),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','flag',ierr)
       allocate(posfrac(3,pub_cell%nat),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','posfrac',ierr)
       allocate(trwfn(ngwf_basis%num),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','trwfn',ierr)
       max_degeneracy = 9
       allocate(psitrpsi(max_degeneracy,max_degeneracy),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','psitrpsi',ierr)
       allocate(zlambda(max_degeneracy),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','zlambda',ierr)
    end if

    ! Prepare translations etc for unfolding
    if (unfold) then
       do iat=1,pub_cell%nat
          ! ndmh: elements is in input-file order, so use pub_orig_atom
          orig_iat = pub_orig_atom(iat)
          posfrac(1,iat) = elements(orig_iat)%centre .DOT. pub_cell%b1
          posfrac(2,iat) = elements(orig_iat)%centre .DOT. pub_cell%b2
          posfrac(3,iat) = elements(orig_iat)%centre .DOT. pub_cell%b3
       end do
       posfrac = posfrac * recip_twopi
       posfrac = modulo(posfrac,1.0_DP)
       do dim=1,3
          if (pub_bs_unfold(dim) == 1) cycle
          trvec = 0.0_DP
          trvec(dim) = 1.0_DP / pub_bs_unfold(dim)
          flag = .false.
          do iat=1,pub_cell%nat
             orig_iat = pub_orig_atom(iat)
             ispec = elements(orig_iat)%species_number
             trpos = posfrac(:,iat) + trvec
             wrapped = 1
             if (trpos(dim) >= 1.0_DP) then
                trpos(dim) = trpos(dim) - 1.0_DP
                wrapped = -1
             end if
             found = .false.
             do jat=1,pub_cell%nat
                ! ndmh: elements is in input-file order
                orig_jat = pub_orig_atom(jat)
                if (flag(jat) .or. elements(orig_jat)%species_number /= ispec) cycle
                dispfrac = posfrac(:,jat) - trpos
                dispabs(1) = dispfrac(1)*pub_cell%a1%x + &
                     dispfrac(2)*pub_cell%a2%x + dispfrac(3)*pub_cell%a3%x
                dispabs(2) = dispfrac(1)*pub_cell%a1%y + &
                     dispfrac(2)*pub_cell%a2%y + dispfrac(3)*pub_cell%a3%y
                dispabs(3) = dispfrac(1)*pub_cell%a1%z + &
                     dispfrac(2)*pub_cell%a2%z + dispfrac(3)*pub_cell%a3%z
                distsq = dispabs(1)*dispabs(1) + dispabs(2)*dispabs(2) + &
                     dispabs(3)*dispabs(3)
                if (distsq < tol) then
                   flag(jat) = .true.
                   translate(iat,dim) = jat * wrapped
                   found = .true.
                   exit
                end if
             end do
             if (.not. found) then
                write(stderr,'(a,i5,a,3(f6.3,a))') &
                     'Error in bandstructure_calculate: &
                     &no atom equivalent to atom',iat, &
                     ' found when translating by (',trvec(1),',',trvec(2), &
                     ',',trvec(3),')'
                call comms_abort
             end if
          end do
       end do
    end if

!CW
    if(abs(pub_dmft_nkpoints)>1) pub_bs_method='TB'
    if(pub_my_node_id==0) write(*,*) 'BAND STRUCTURE WITH METHOD : ', pub_bs_method
!END CW

    ! Choose method
    select case (pub_bs_method)

    ! pdh: there is almost certainly more scope for code-sharing below
    ! pdh: and also for modularity that could assist future developments
    ! pdh: but the different approaches can be optimised rather
    ! pdh: differently and one may well be deprecated so for the time
    ! pdh: being separate the code for clarity

    !==================================!
    !                                  !
    !    TIGHT-BINDING STYLE METHOD    !
    !                                  !
    !==================================!

    case ('TB')

       ! Diagonalization workspace
       lwork = -1
       if (unfold) then
          call zhegv(1,'V','L',ngwf_basis%num,hm_sq,ngwf_basis%num,zov_sq, &
               ngwf_basis%num,en_ks,zwork,lwork,dwork(1,1),ierr)
       else
          call zhegv(1,'N','L',ngwf_basis%num,hm_sq,ngwf_basis%num,zov_sq, &
               ngwf_basis%num,en_ks,zwork,lwork,dwork(1,1),ierr)
       end if
       if (ierr /= 0) then
          lwork = 2*ngwf_basis%num
       else
          lwork = nint(real(zwork(1),kind=DP))
          deallocate(zwork,stat=ierr)
          call utils_dealloc_check('bandstructure_calculate','zwork',ierr)
          allocate(zwork(lwork),stat=ierr)
          call utils_alloc_check('bandstructure_calculate','zwork',ierr)
       end if

       ! Allocate arrays for this method
       hamidxlen = sparse_index_length(ham%ham(1))
       allocate(hamidx(hamidxlen),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','hamidx',ierr)
       call sparse_generate_index(hamidx,ham%ham(1))
       overidxlen = sparse_index_length(rep%overlap)
       allocate(overidx(overidxlen),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','overidx',ierr)
       call sparse_generate_index(overidx,rep%overlap)
       allocate(ham_ks(0:pub_total_num_nodes-1),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','ham_ks',ierr)
       allocate(over_ks(0:pub_total_num_nodes-1),stat=ierr)
!CW
       if(abs(pub_dmft_nkpoints)>1)then
        if(.not.present(hub)) then
           write(*,*) 'ERROR hub argument should be present for k-point DMFT'
           stop
        endif 
        allocate(hubT_ks(0:pub_total_num_nodes-1),stat=ierr)
         do node=0,pub_total_num_nodes-1
           call sparse_create(hubT_ks(node),rep%hub_overlap_t, iscmplx=.true. )
         end do
       endif 
!END CW
       call utils_alloc_check('bandstructure_calculate','over_ks',ierr)
       do node=0,pub_total_num_nodes-1
           call sparse_create(ham_ks(node),ham%ham(1),iscmplx=.true.)
           call sparse_create(over_ks(node),rep%overlap,iscmplx=.true.)
       end do
       max_ngwfs_on_atom = maxval(ngwf_basis%num_on_atom)
       allocate(block(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','block',ierr)
       block = (0.0_DP,0.0_DP)

       ! Loop over spins
       do is=1,nspins

          ! Main loop over unique k-points
          do ikpt=1,nkpts_red,pub_total_num_nodes
             do node=0,pub_total_num_nodes-1
                if (ikpt+node > nkpts_red) exit

                ! Find k-point in original list that matches
                do jkpt=1,nkpts
                   if (uni_entry(jkpt) == ikpt+node) exit
                end do
                if (jkpt > nkpts) then
                   write(stdout,*) 'Error in bandstructure_calculate: &
                        &malformed list of unique k-points'
                   call comms_abort
                end if

                kfrac(1) = kpt(jkpt) .DOT. pub_cell%a1
                kfrac(2) = kpt(jkpt) .DOT. pub_cell%a2
                kfrac(3) = kpt(jkpt) .DOT. pub_cell%a3
                kfrac = kfrac * recip_twopi

                ! Create overlap matrix for this k-point and spin
                call sparse_copy(over_ks(node),rep%overlap)
                do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
                   iat = loc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
                   ! ndmh: elements is in input-file order
                   orig_iat = pub_orig_atom(iat)
                   do jdx=overidx(loc_iat),overidx(loc_iat+1)-1
                      jat = overidx(jdx)
                      orig_jat = pub_orig_atom(jat)
                      atdisp = elements(orig_iat)%centre - &
                           elements(orig_jat)%centre
                      dispfrac(1) = atdisp .DOT. pub_cell%b1
                      dispfrac(2) = atdisp .DOT. pub_cell%b2
                      dispfrac(3) = atdisp .DOT. pub_cell%b3
                      dispfrac = dispfrac * recip_twopi
                      zphase = (1.0_DP,0.0_DP)
                      do dim=1,3
                         if (dispfrac(dim) > 0.5_DP) &
                              zphase = zphase * exp(i2pi * kfrac(dim))
                         if (dispfrac(dim) < -0.5_DP) &
                              zphase = zphase * exp(-i2pi * kfrac(dim))
                      end do

                      if (zphase /= (1.0_DP,0.0_DP)) then
                         call sparse_get_block(block,over_ks(node),jat,iat)
                         block = block * zphase
                         call sparse_put_block(block,over_ks(node),jat,iat)
                      end if
                   end do
                end do
                call sparse_convert(bmat,over_ks(node))

!CW
!               if (node == pub_my_node_id) zov_sq = bmat
                if (node == pub_my_node_id) then
                   zov_sq = bmat
                   if(abs(pub_dmft_nkpoints)>1)then
                     open(unit=5001,file='store_ovk_'//trim(adjustl(toString( ikpt+pub_my_node_id )))//'_'//trim(adjustl(toString( is ))),form='unformatted')
                     write(5001) zov_sq
                     close(5001)
                   endif
                 endif
!END CW

!CW
              if(abs(pub_dmft_nkpoints)>1.and.is==1)then

                call sparse_copy(hubT_ks(node),rep%hub_overlap_t)

                do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
                   iat = loc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
                   orig_iat = pub_orig_atom(iat)
                   do jdx=overidx(loc_iat),overidx(loc_iat+1)-1
                      jat      = overidx(jdx)
                      orig_jat = pub_orig_atom(jat)

                      atdisp = elements(orig_iat)%centre - elements(orig_jat)%centre
                      dispfrac(1) = atdisp .DOT. pub_cell%b1
                      dispfrac(2) = atdisp .DOT. pub_cell%b2
                      dispfrac(3) = atdisp .DOT. pub_cell%b3
                      dispfrac    = dispfrac * recip_twopi
                      zphase = (1.0_DP,0.0_DP)
                      do dim=1,3
                         if (dispfrac(dim) >  0.5_DP) &
                              zphase = zphase * exp( i2pi * kfrac(dim))
                         if (dispfrac(dim) < -0.5_DP) &
                              zphase = zphase * exp(-i2pi * kfrac(dim))
                      end do

                      if (zphase /= (1.0_DP,0.0_DP)) then
                        nproj=hub_proj_basis%num_on_atom(jat)
                        if(nproj>0)then
                          do hub_proj=1,nproj
                           do ngwfs_on_atom=1,ngwf_basis%num_on_atom(iat)

                              ii1 = hub_proj_basis%first_on_atom(jat) + hub_proj      -1
                              ii2 =     ngwf_basis%first_on_atom(iat) + ngwfs_on_atom -1  
                              !jat,iat in overlap
                              call sparse_get_element(cblock,hubT_ks(node),ii1,ii2)
                              cblock=cblock *       zphase
                              call sparse_put_element(cblock,hubT_ks(node),ii1,ii2)
                           enddo
                          enddo
                        else
                           !......!
                        endif

                      end if
                   end do
                end do
                call sparse_convert(bmatT_rec,hubT_ks(node))
                bmat_rec=transpose(conjg(bmatT_rec))

                if (node == pub_my_node_id) then
                     open(unit=5001,file='store_hubTk_'//trim(adjustl(toString( ikpt+pub_my_node_id ))),form='unformatted')
                     write(5001) bmatT_rec
                     close(5001)
                     open(unit=5001,file='store_hubk_'//trim(adjustl(toString( ikpt+pub_my_node_id ))),form='unformatted')
                     write(5001) bmat_rec
                     close(5001)
                endif

              endif
!END CW

                ! Create Hamiltonian matrix for this k-point and spin
                call sparse_copy(ham_ks(node),ham%ham(is))
                do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
                   iat = loc_iat + pub_first_atom_on_node(pub_my_node_id) - 1
                   ! ndmh: elements is in input-file order
                   orig_iat = pub_orig_atom(iat)
                   do jdx=hamidx(loc_iat),hamidx(loc_iat+1)-1
                      jat = hamidx(jdx)
                      orig_jat = pub_orig_atom(jat)
                      atdisp = elements(orig_iat)%centre - &
                           elements(orig_jat)%centre
                      dispfrac(1) = atdisp .DOT. pub_cell%b1
                      dispfrac(2) = atdisp .DOT. pub_cell%b2
                      dispfrac(3) = atdisp .DOT. pub_cell%b3
                      dispfrac = dispfrac * recip_twopi
                      zphase = (1.0_DP,0.0_DP)
                      do dim=1,3
                         if (dispfrac(dim) >  0.5_DP) &
                              zphase = zphase * exp( i2pi * kfrac(dim))
                         if (dispfrac(dim) < -0.5_DP) &
                              zphase = zphase * exp(-i2pi * kfrac(dim))
                      end do
                      if (zphase /= (1.0_DP,0.0_DP)) then
                         call sparse_get_block(block,ham_ks(node),jat,iat)
                         block = block * zphase
                         call sparse_put_block(block,ham_ks(node),jat,iat)
                      end if
                   end do
                end do
                call sparse_convert(bmat,ham_ks(node))

!CW
             !  if (node == pub_my_node_id) hm_sq = bmat
                if (node == pub_my_node_id) then
                   hm_sq = bmat
                   if(abs(pub_dmft_nkpoints)>1)then
                     open(unit=5001,file='store_hamk_'//trim(adjustl(toString( ikpt+pub_my_node_id )))//'_'//trim(adjustl(toString( is ))),form='unformatted')
                     write(5001) hm_sq
                     close(5001)
                   endif
                 endif
!END CW

             end do

!CW
             if(abs(pub_dmft_nkpoints)>1) cycle
!END CW

             ! Diagonalise different k-points on different nodes
             if (ikpt+pub_my_node_id <= nkpts_red) then

                ! Diagonalise Hamiltonian
                if (unfold) then
                   bmat = zov_sq
                   call zhegv(1,'V','L',ngwf_basis%num,hm_sq,ngwf_basis%num, &
                        bmat,ngwf_basis%num,en_ks,zwork,lwork,dwork(1,1),ierr)
                else
                   call zhegv(1,'N','L',ngwf_basis%num,hm_sq,ngwf_basis%num, &
                        zov_sq,ngwf_basis%num,en_ks,zwork,lwork,dwork(1,1),ierr)
                end if
                if (ierr /= 0) then
                   write(stderr,'(a,i6)') 'Error in bandstructure_calculate: &
                        &zhegv failed with code ',ierr
                   call comms_abort
                end if

                if (unfold) then

                   ! Recalculate fractional coordinates of relevant k-point
                   do jkpt=1,nkpts
                      if (uni_entry(jkpt) == ikpt+pub_my_node_id) exit
                   end do
                   kfrac(1) = kpt(jkpt) .DOT. pub_cell%a1
                   kfrac(2) = kpt(jkpt) .DOT. pub_cell%a2
                   kfrac(3) = kpt(jkpt) .DOT. pub_cell%a3
                   kfrac = kfrac * recip_twopi

                   call internal_unfold

                end if

             end if

             ! Communicate results between nodes
             call comms_barrier

             do node=0,pub_total_num_nodes-1
                if (ikpt+node > nkpts_red) exit

                if (node == pub_my_node_id) dwork(:,1) = en_ks
                call comms_bcast(node,dwork(1,1),ngwf_basis%num)
                do jkpt=1,nkpts
                   if (uni_entry(jkpt) == ikpt+node) &
                        eig(:,jkpt,is) = dwork(:,1)
                end do

                if (unfold) then
                   if (node == pub_my_node_id) dwork = kpt_ks
                   call comms_bcast(node,dwork(1,1),3*ngwf_basis%num)
                   do jkpt=1,nkpts
                      if (uni_entry(jkpt) == ikpt+node) &
                           kunfold(:,:,jkpt,is) = dwork
                   end do
                end if

             end do

          end do    ! unique k-points

       end do   ! spins

       ! Deallocate arrays for this method
       deallocate(block,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','block',ierr)
       do node=pub_total_num_nodes-1,0,-1
           call sparse_destroy(over_ks(node))
           call sparse_destroy(ham_ks(node))
       end do
!CW
       if(abs(pub_dmft_nkpoints)>1)then
         do node=0,pub_total_num_nodes-1
           call sparse_destroy( hubT_ks(node) )
         end do
        deallocate(hubT_ks,stat=ierr)
       endif
!END CW


       deallocate(over_ks,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','over_ks',ierr)
       deallocate(ham_ks,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','ham_ks',ierr)
       deallocate(overidx,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','overidx',ierr)
       deallocate(hamidx,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','hamidx',ierr)

    !==================================!
    !                                  !
    !         K.P STYLE METHOD         !
    !                                  !
    !==================================!

    case ('KP')

       ! Diagonalization workspace
       lwork = -1
       if (unfold) then
          call zheev('V','L',ngwf_basis%num,hm_sq,ngwf_basis%num,en_ks, &
               zwork,lwork,dwork(1,1),ierr)
       else
          call zheev('N','L',ngwf_basis%num,hm_sq,ngwf_basis%num,en_ks, &
               zwork,lwork,dwork(1,1),ierr)
       end if
       if (ierr /= 0) then
          lwork = 2*ngwf_basis%num
       else
          lwork = nint(real(zwork(1),kind=DP))
          deallocate(zwork,stat=ierr)
          call utils_dealloc_check('bandstructure_calculate','zwork',ierr)
          allocate(zwork(lwork),stat=ierr)
          call utils_alloc_check('bandstructure_calculate','zwork',ierr)
       end if

       ! Prepare k-independent part of Hamiltonian
       allocate(ham0(nspins),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','ham0',ierr)
       do is=1,nspins
          call sparse_create(ham0(is),ham%ham(is),iscmplx=.true.)
          call sparse_copy(ham0(is),ham%ham(is))
          if (pub_any_nl_proj) call sparse_axpy(ham0(is),rep%nonlocpot(1),-1.0_DP)
       end do

       ! Calculate grad matrix elements
       do dim=1,3
          call sparse_create(grad(dim),rep%overlap)
       end do
       call integrals_grad(grad, rep%ngwfs_on_grid, ngwf_basis, &
            rep%ngwfs_on_grid, ngwf_basis)

       ! Initialise nonlocal matrix and Hamiltonian
       if (pub_any_nl_proj) then
          allocate(nl_k(0:pub_total_num_nodes-1),stat=ierr)
          call utils_alloc_check('bandstructure_calculate','nl_k',ierr)
       end if
       allocate(ham_ks(0:pub_total_num_nodes-1),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','ham_ks',ierr)
       do node=0,pub_total_num_nodes-1
          if (pub_any_nl_proj) &
               call sparse_create(nl_k(node),rep%nonlocpot(1),iscmplx=.true.)
          call sparse_create(ham_ks(node),ham%ham(1),iscmplx=.true.)
       end do

       ! Convert overlap matrix to dense format and obtain
       ! Cholesky decomposition
       allocate(ov_sq(ngwf_basis%num,ngwf_basis%num),stat=ierr)
       call utils_alloc_check('bandstructure_calculate','ov_sq',ierr)
       call sparse_convert(ov_sq,rep%overlap)
       zov_sq = ov_sq

       call dpotrf('L',ngwf_basis%num,ov_sq,ngwf_basis%num,ierr)
       if (ierr /= 0) then
          write(stderr,'(a,i6)') 'Error in bandstructure_calculate: &
               &dpotrf failed with code ',ierr
          call comms_abort
       end if

       ! ndmh: create temporary matrices for real/complex sp overlaps
       if (pub_any_nl_proj) then
          call sparse_create(sp_overlap,rep%sp_overlap,iscmplx=.true.)
          call sparse_create(sp_overlap_rc,rep%sp_overlap)
       end if

       ! Loop over spins
       do is=1,nspins

          ! Main loop over unique k-points
          do ikpt=1,nkpts_red,pub_total_num_nodes
             do node=0,pub_total_num_nodes-1
                if (ikpt+node > nkpts_red) exit

!CW
                write(*,*) 'NODE / K POINT : ', node,ikpt+node
!ENDCW

                ! Find k-point in original list that matches
                do jkpt=1,nkpts
                   if (uni_entry(jkpt) == ikpt+node) exit
                end do
                if (jkpt > nkpts) then
                   write(stdout,*) 'Error in bandstructure_calculate: &
                        &malformed list of unique k-points'
                   call comms_abort
                end if

                kcart(1) = kpt(jkpt)%x ; kcart(2) = kpt(jkpt)%y
                kcart(3) = kpt(jkpt)%z

                if (node == pub_my_node_id) then
                   kfrac(1) = kpt(jkpt) .DOT. pub_cell%a1
                   kfrac(2) = kpt(jkpt) .DOT. pub_cell%a2
                   kfrac(3) = kpt(jkpt) .DOT. pub_cell%a3
                   kfrac = kfrac * recip_twopi
                end if

                ! Create Hamiltonian at this k-point and spin
                call sparse_copy(ham_ks(node),ham0(is))

                ! Calculate the nonlocal matrix for this k-point
                if (pub_any_nl_proj) then

                   ! ndmh: calculate real part of <ngwf|projector> overlap matrix
                   ! ndmh: at this kpt
                   call projectors_func_ovlp_box(sp_overlap_rc, &
                        rep%ngwfs_on_grid,ngwf_basis,proj_basis, nlps_projectors, &
                        kshift=kpt(jkpt),swap_rc=.false.)

                   ! ndmh: add real part to complex sp_overlap matrix
                   call sparse_copy(sp_overlap,sp_overlap_rc)

                   ! ndmh: calculate imaginary part of <ngwf|projector> overlap matrix
                   call projectors_func_ovlp_box(sp_overlap_rc, &
                        rep%ngwfs_on_grid,ngwf_basis,proj_basis, nlps_projectors, &
                        kshift=kpt(jkpt),swap_rc=.true.)

                   ! ndmh: add imaginary part to complex sp_overlap matrix
                   call sparse_axpy(sp_overlap,sp_overlap_rc,(0.0_DP,1.0_DP))

                   ! ndmh: calculate nonlocal pseudopotential matrix at this kpt
                   call pseudopotentials_nonlocal_mat(nl_k(node),sp_overlap)

                   call sparse_axpy(ham_ks(node),nl_k(node),1.0_DP)

                end if
                do dim=1,3
                   call sparse_axpy(ham_ks(node),grad(dim), &
                        cmplx(0.0_DP,-kcart(dim),kind=DP))
                end do
                call sparse_axpy(ham_ks(node),rep%overlap, &
                     0.5_DP*(kcart(1)*kcart(1) + kcart(2)*kcart(2) + &
                     kcart(3)*kcart(3)))
                call sparse_convert(bmat,ham_ks(node))

                if (node == pub_my_node_id) hm_sq = bmat
             end do


             ! Diagonalise different k-points on different nodes
             if (ikpt+pub_my_node_id <= nkpts_red) then

                ! Copy overlap Cholesky factorization
                do ingwf=1,ngwf_basis%num
                   do jngwf=ingwf,ngwf_basis%num
                      bmat(jngwf,ingwf) = &
                           cmplx(ov_sq(jngwf,ingwf),0.0_DP,kind=DP)
                   end do
                end do

                ! Convert Hamiltonian to orthogonal form
                call zhegst(1,'L',ngwf_basis%num,hm_sq,ngwf_basis%num,bmat, &
                     ngwf_basis%num,ierr)
                if (ierr /= 0) then
                   write(stderr,'(a,i6)') 'Error in bandstructure_calculate: &
                        &zhegst failed with code ',ierr
                   call comms_abort
                end if

                ! Diagonalise Hamiltonian
                if (unfold) then
                   call zheev('V','L',ngwf_basis%num,hm_sq,ngwf_basis%num,en_ks, &
                        zwork,lwork,dwork(1,1),ierr)
                else
                   call zheev('N','L',ngwf_basis%num,hm_sq,ngwf_basis%num,en_ks, &
                        zwork,lwork,dwork(1,1),ierr)
                end if
                if (ierr /= 0) then
                   write(stderr,'(a,i6)') 'Error in bandstructure_calculate: &
                        &zheev failed with code ',ierr
                   call comms_abort
                end if

                if (unfold) then

                   ! Recalculate fractional coordinates of relevant k-point
                   do jkpt=1,nkpts
                      if (uni_entry(jkpt) == ikpt+pub_my_node_id) exit
                   end do
                   kfrac(1) = kpt(jkpt) .DOT. pub_cell%a1
                   kfrac(2) = kpt(jkpt) .DOT. pub_cell%a2
                   kfrac(3) = kpt(jkpt) .DOT. pub_cell%a3
                   kfrac = kfrac * recip_twopi

                   ! Recover eigenvectors
                   call ztrsm('L','L','T','N',ngwf_basis%num,ngwf_basis%num, &
                        1.0_DP,bmat,ngwf_basis%num,hm_sq,ngwf_basis%num)

                   call internal_unfold

                end if
             end if

             ! Communicate results between nodes
             call comms_barrier

             do node=0,pub_total_num_nodes-1
                if (ikpt+node > nkpts_red) exit

                if (node == pub_my_node_id) dwork(:,1) = en_ks
                call comms_bcast(node,dwork(1,1),ngwf_basis%num)
                do jkpt=1,nkpts
                   if (uni_entry(jkpt) == ikpt+node) &
                        eig(:,jkpt,is) = dwork(:,1)
                end do

                if (unfold) then
                   if (node == pub_my_node_id) dwork = kpt_ks
                   call comms_bcast(node,dwork)
                   do jkpt=1,nkpts
                      if (uni_entry(jkpt) == ikpt+node) &
                           kunfold(:,:,jkpt,is) = dwork
                   end do
                end if

             end do

          end do    ! unique k-points

       end do   ! spins

       ! Deallocate arrays for this method
       if (pub_any_nl_proj) then
          call sparse_destroy(sp_overlap_rc)
          call sparse_destroy(sp_overlap)
       end if
       deallocate(ov_sq,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','ov_sq',ierr)
       do node=pub_total_num_nodes-1,0,-1
          call sparse_destroy(ham_ks(node))
          if (pub_any_nl_proj) call sparse_destroy(nl_k(node))
       end do
       deallocate(ham_ks,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','ham_ks',ierr)
       if (pub_any_nl_proj) then
          deallocate(nl_k,stat=ierr)
          call utils_dealloc_check('bandstructure_calculate','nl_k',ierr)
       end if
       do dim=3,1,-1
          call sparse_destroy(grad(dim))
       end do
       do is=1,nspins
          call sparse_destroy(ham0(is))
       end do
       deallocate(ham0,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','ham0',ierr)

    end select

    ! Write results to output file
    ! ndmh: write all bands if bs_num_eigenvalues is left at default value
    if (pub_bs_num_eigenvalues < 0) pub_bs_num_eigenvalues = maxval(rep%n_occ(:))
    if (pub_on_root .and. pub_bs_num_eigenvalues > 0) then
       write(stdout,'(/a)') &
            ' + ============================================================= +'
       write(stdout,'(a)') &
            ' +                     Electronic energies                       +'
       write(stdout,'(a)') &
            ' +                     -------------------                       +'
       write(stdout,'(a)') &
            ' +                                                               +'
       if (unfold) then
          write(stdout,'(a)') ' +   Band number   Energy in Ha          &
               &Unfolded k-point        +'
       else
          write(stdout,'(a)') ' +   Band number   Energy in Ha          &
               &                        +'
       end if
       write(stdout,'(a)') &
            ' + ============================================================= +'
       write(stdout,'(a)') &
            ' +                                                               +'
       do is=1,nspins
          nlower = max(rep%n_occ(is)-pub_bs_num_eigenvalues+1,1)
          nupper = min(rep%n_occ(is)+pub_bs_num_eigenvalues,ngwf_basis%num)
          ikpt = 0
          do ipath=1,pub_bs_kpoint_path_length
             kstart = pub_bs_kpoint_path_start(1,ipath) * pub_cell%b1 + &
                  pub_bs_kpoint_path_start(2,ipath) * pub_cell%b2 + &
                  pub_bs_kpoint_path_start(3,ipath) * pub_cell%b3
             kend = pub_bs_kpoint_path_end(1,ipath) * pub_cell%b1 + &
                  pub_bs_kpoint_path_end(2,ipath) * pub_cell%b2 + &
                  pub_bs_kpoint_path_end(3,ipath) * pub_cell%b3
             kseg = kend - kstart
             seglen = geometry_magnitude(kseg)
             nsegpts = max(int(seglen / pub_bs_kpoint_path_spacing),1)
             do isegpt=0,nsegpts
                ikpt = ikpt + 1
                kfrac = pub_bs_kpoint_path_start(:,ipath) + &
                     (isegpt / real(nsegpts,kind=DP)) * &
                     (pub_bs_kpoint_path_end(:,ipath) - &
                     pub_bs_kpoint_path_start(:,ipath))
                write(stdout,'(a)') ' + --------------------------------------&
                     &----------------------- +'
                write(stdout,'(a,i1,a,i5,a,f8.5,1x,f8.5,1x,f8.5,a)') &
                     ' +  Spin=',is,' kpt=',ikpt,' (',kfrac,')                +'
                write(stdout,'(a)') ' + --------------------------------------&
                     &----------------------- +'
                write(stdout,'(a)') ' +                                       &
                     &                        +'
                do ingwf=nlower,nupper
                   if (unfold) then
                      kfrac = kunfold(ingwf,:,ikpt,is)
                      write(stdout,'(a,i12,1x,f15.6,5x,a,3f8.4,a,4x,a)') &
                           ' +',ingwf,eig(ingwf,ikpt,is),'(', &
                           kfrac,')','+'
                   else
                      write(stdout,'(a,i12,1x,f15.6,35x,a)') ' +',ingwf, &
                           eig(ingwf,ikpt,is),'+'
                   end if
                   if (ingwf == rep%n_occ(is)) write(stdout,'(a)') &
                     ' +         ... Fermi level ...                          &
                     &         +'
                end do
                write(stdout,'(a)') ' +                                       &
                     &                        +'
             end do
          end do
       end do
       write(stdout,'(a)') &
            ' + ============================================================= +'
    end if

    ! Write out .bands file
    if (pub_on_root) then

       ! Open output file
       bands_unit = utils_unit()
       if (.not.pub_cond_calculate) then
          write(filename,'(2a)') trim(pub_rootname),'_BS.bands'
       else
          write(filename,'(4a)') trim(pub_rootname),'_',trim(ham_type),'_BS.bands'
       end if
       open(bands_unit,file=trim(filename),iostat=ierr)
       if (ierr /= 0) then
          write(stderr,'(3a,i6)') &
               'Error in bandstructure_calculate: opening "', &
               trim(filename),'" failed with code ',ierr
          call comms_abort
       end if

       write(bands_unit,'(a,i6)') 'Number of k-points',nkpts
       write(bands_unit,'(a,i2)') 'Number of spin components',nspins
       write(bands_unit,'(a,2i8)') 'Number of electrons ',rep%n_occ(1:nspins)
       write(bands_unit,'(a,2i8)') 'Number of eigenvalues ',ngwf_basis%num
       if (rep%n_occ(1)>0) then
          efermi = maxval(eig(rep%n_occ(1),:,1))
       else
          efermi = 0_DP
       endif
       if (nspins == 1) then
          write(bands_unit,'(a,i8)') 'Number of eigenvalues ',ngwf_basis%num
          write(bands_unit,'(a,f12.6)') 'Fermi energy (in atomic units) ', &
               efermi
       else
          write(bands_unit,'(a,2i8)') 'Number of eigenvalues ',ngwf_basis%num, &
               ngwf_basis%num
          if (rep%n_occ(2)>0) then
             efermi = max(efermi,maxval(eig(rep%n_occ(2),:,2)))
          else
             efermi = 0_DP
          endif
          write(bands_unit,'(a,2f12.6)') 'Fermi energies (in atomic units) ', &
               efermi,efermi
       end if
       write(bands_unit,'(a)') 'Unit cell vectors'
       write(bands_unit,'(3f12.6)') pub_cell%a1%x,pub_cell%a1%y,pub_cell%a1%z
       write(bands_unit,'(3f12.6)') pub_cell%a2%x,pub_cell%a2%y,pub_cell%a2%z
       write(bands_unit,'(3f12.6)') pub_cell%a3%x,pub_cell%a3%y,pub_cell%a3%z
       ikpt = 0
       do ipath=1,pub_bs_kpoint_path_length
          kstart = pub_bs_kpoint_path_start(1,ipath) * pub_cell%b1 + &
               pub_bs_kpoint_path_start(2,ipath) * pub_cell%b2 + &
               pub_bs_kpoint_path_start(3,ipath) * pub_cell%b3
          kend = pub_bs_kpoint_path_end(1,ipath) * pub_cell%b1 + &
               pub_bs_kpoint_path_end(2,ipath) * pub_cell%b2 + &
               pub_bs_kpoint_path_end(3,ipath) * pub_cell%b3
          kseg = kend - kstart
          seglen = geometry_magnitude(kseg)
          nsegpts = max(int(seglen / pub_bs_kpoint_path_spacing),1)
          do isegpt=0,nsegpts
             ikpt = ikpt + 1
             kfrac = pub_bs_kpoint_path_start(:,ipath) + &
                  (isegpt / real(nsegpts,kind=DP)) * &
                  (pub_bs_kpoint_path_end(:,ipath) - &
                  pub_bs_kpoint_path_start(:,ipath))
             write(bands_unit,'(a,i6,4f12.8)') 'K-point',ikpt,kfrac,1.0_DP/nkpts
             do is=1,nspins
                write(bands_unit,'(a,i2)') 'Spin component',is
                do ingwf=1,ngwf_basis%num
                   write(bands_unit,'(f14.8)') eig(ingwf,ikpt,is)
                end do
             end do
          end do
       end do

       ! Close output file
       close(bands_unit,iostat=ierr)
       if (ierr /= 0) then
          write(stderr,'(3a,i6)') &
               'Error in bandstructure_calculate: closing "', &
               trim(filename),'" failed with code ',ierr
          call comms_abort
       end if

    end if

    ! Write out .agr file
    if (pub_on_root) then

       ! Open output file
       agr_unit = utils_unit()
       write(filename,'(2a)') trim(pub_rootname),'_BS.agr'
       open(agr_unit,file=trim(filename),iostat=ierr)
       if (ierr /= 0) then
          write(stderr,'(3a,i6)') &
               'Error in bandstructure_calculate: opening "', &
               trim(filename),'" failed with code ',ierr
          call comms_abort
       end if

       write(agr_unit,'(a)') '@version 50112'
       write(agr_unit,'(a)') '@with g0'
       write(agr_unit,'(a,f12.6)') '@    world xmin ',0.0_DP
       write(agr_unit,'(a,f12.6)') '@    world xmax ',pathlen
       write(agr_unit,'(a,f12.6)') '@    world ymin ', &
            real(floor(minval(eig)),kind=DP)
       write(agr_unit,'(a,f12.6)') '@    world ymax ', &
            real(ceiling(maxval(eig)),kind=DP)
       write(agr_unit,'(a)') '@    yaxis  tick on'
       write(agr_unit,'(a,f12.6)') '@    yaxis  tick major ',1.0_DP
       write(agr_unit,'(a)') '@    yaxis  tick minor off'
       write(agr_unit,'(a)') '@    xaxis  tick major off'
       write(agr_unit,'(a)') '@    xaxis  tick minor off'
       write(agr_unit,'(a)') '@    xaxis  ticklabel on'
       write(agr_unit,'(a)') '@    xaxis  ticklabel type spec'
       write(agr_unit,'(a)') '@    xaxis  tick type spec'
       write(agr_unit,'(a,i2,a,f12.6)') '@    xaxis  tick ',0,', ',0.0_DP
       write(agr_unit,'(a,2(i2,a))') '@    xaxis  ticklabel ',0,', "\6k\4\s', &
            0,'"'
       seglen = 0.0_DP
       do ipath=1,pub_bs_kpoint_path_length
          kstart = pub_bs_kpoint_path_start(1,ipath) * pub_cell%b1 + &
               pub_bs_kpoint_path_start(2,ipath) * pub_cell%b2 + &
               pub_bs_kpoint_path_start(3,ipath) * pub_cell%b3
          kend = pub_bs_kpoint_path_end(1,ipath) * pub_cell%b1 + &
               pub_bs_kpoint_path_end(2,ipath) * pub_cell%b2 + &
               pub_bs_kpoint_path_end(3,ipath) * pub_cell%b3
          kseg = kend - kstart
          seglen = seglen + geometry_magnitude(kseg)
          write(agr_unit,'(a,i2,a,f12.6)') '@    xaxis  tick ',ipath,', ',seglen
          write(agr_unit,'(a,2(i2,a))') '@    xaxis  ticklabel ',ipath, &
               ', "\6k\4\s',ipath,'"'
       end do
       write(agr_unit,'(a,i3)') '@    xaxis  tick spec ', &
            pub_bs_kpoint_path_length+1

       do is=1,nspins
          do ingwf=1,ngwf_basis%num
             write(agr_unit,'(a,i3)') '@TYPE xy'
             ikpt = 0
             pathlen = 0.0_DP
             do ipath=1,pub_bs_kpoint_path_length
                kstart = pub_bs_kpoint_path_start(1,ipath) * pub_cell%b1 + &
                     pub_bs_kpoint_path_start(2,ipath) * pub_cell%b2 + &
                     pub_bs_kpoint_path_start(3,ipath) * pub_cell%b3
                kend = pub_bs_kpoint_path_end(1,ipath) * pub_cell%b1 + &
                     pub_bs_kpoint_path_end(2,ipath) * pub_cell%b2 + &
                     pub_bs_kpoint_path_end(3,ipath) * pub_cell%b3
                kseg = kend - kstart
                seglen = geometry_magnitude(kseg)
                nsegpts = max(int(seglen / pub_bs_kpoint_path_spacing),1)
                do isegpt=0,nsegpts
                   ikpt = ikpt + 1
                   if (isegpt > 0) pathlen = pathlen + seglen / nsegpts
                   write(agr_unit,'(2f12.6)') pathlen,eig(ingwf,ikpt,is)
                end do
             end do
             write(agr_unit,'(a,i3)') '&'
          end do
       end do

       ! Close output file
       close(agr_unit,iostat=ierr)
       if (ierr /= 0) then
          write(stderr,'(3a,i6)') &
               'Error in bandstructure_calculate: closing "', &
               trim(filename),'" failed with code ',ierr
          call comms_abort
       end if

    end if

    ! Destroy workspace
    deallocate(zwork,stat=ierr)
    call utils_dealloc_check('bandstructure_calculate','zwork',ierr)
    if (unfold) then
       deallocate(zlambda,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','zlambda',ierr)
       deallocate(psitrpsi,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','psitrpsi',ierr)
       deallocate(trwfn,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','trwfn',ierr)
       deallocate(posfrac,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','posfrac',ierr)
       deallocate(flag,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','flag',ierr)
       deallocate(kunfold,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','kunfold',ierr)
       deallocate(translate,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','translate',ierr)
       deallocate(kpt_ks,stat=ierr)
       call utils_dealloc_check('bandstructure_calculate','kpt_ks',ierr)
    end if
    deallocate(zov_sq,stat=ierr)
    call utils_dealloc_check('bandstructure_calculate','zov_sq',ierr)
!CW
    if(abs(pub_dmft_nkpoints)>1)then
      deallocate(bmat_rec,bmatT_rec)
    endif
!END CW  
    deallocate(bmat,stat=ierr)
    call utils_dealloc_check('bandstructure_calculate','bmat',ierr)
    deallocate(dwork,stat=ierr)
    call utils_dealloc_check('bandstructure_calculate','dwork',ierr)
    deallocate(hm_sq,stat=ierr)
    call utils_dealloc_check('bandstructure_calculate','hm_sq',ierr)
    deallocate(en_ks,stat=ierr)
    call utils_dealloc_check('bandstructure_calculate','en_ks',ierr)
    deallocate(eig,stat=ierr)
    call utils_dealloc_check('bandstructure_calculate','eig',ierr)
    deallocate(uni_entry,stat=ierr)
    call utils_dealloc_check('bandstructure_calculate','uni_entry',ierr)
    deallocate(kpt,stat=ierr)
    call utils_dealloc_check('bandstructure_calculate','kpt',ierr)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving bandstructure_calculate'
#endif

  contains

!CW
  character(3) function toString(i)
   integer :: i
    if(i>999)then
     write(*,*) 'ERROR function toString 1' ;stop
    endif
    write(toString,'(i3)') i
   end function
!END CW

    subroutine internal_unfold

      ! Local variables
      integer :: iwfn,jwfn,kwfn             ! Wavefunction counters
      integer :: degen                      ! Degeneracy
      real(kind=DP), parameter :: DEGEN_TOL = 1.0e-6_DP
      complex(kind=DP) :: dotp              ! Dot product

      ! Loop over eigenstates
      iwfn = 0
      degen = 1
      do
         iwfn = iwfn + degen
         if (iwfn > ngwf_basis%num) exit

         ! Count degeneracy of this state
         do jwfn=iwfn+1,ngwf_basis%num
            if (en_ks(jwfn) - en_ks(iwfn) > DEGEN_TOL) exit
         end do
         degen = jwfn - iwfn

         ! Reallocate memory if necessary
         if (degen > max_degeneracy) then
            deallocate(zlambda,stat=ierr)
            call utils_dealloc_check( &
                 'internal_unfold (bandstructure_calculate)','zlambda',ierr)
            deallocate(psitrpsi,stat=ierr)
            call utils_dealloc_check( &
                 'internal_unfold (bandstructure_calculate)','psitrpsi',ierr)
            allocate(psitrpsi(degen,degen),stat=ierr)
            call utils_alloc_check( &
                 'internal_unfold (bandstructure_calculate)','psitrpsi',ierr)
            allocate(zlambda(degen),stat=ierr)
            call utils_alloc_check( &
                 'internal_unfold (bandstructure_calculate)','zlambda',ierr)
            max_degeneracy = degen
         end if

         ! Attempt displacements along each direction
         do dim=1,3

            if (pub_bs_unfold(dim) == 1) then
               kpt_ks(iwfn:iwfn+degen-1,dim) = 0.0_DP
            else

               ! Loop over degenerate states
               do jwfn=1,degen

                  ! Apply translation operator
                  do ingwf=1,ngwf_basis%num
                     iat = ngwf_basis%atom_of_func(ingwf)
                     jat = translate(iat,dim)
                     wrapped = sign(1,jat)
                     jat = abs(jat)
                     jngwf = ngwf_basis%first_on_atom(jat) + ingwf - &
                          ngwf_basis%first_on_atom(iat)
                     trwfn(jngwf) = hm_sq(ingwf,iwfn+jwfn-1)
                     if (wrapped == -1) trwfn(jngwf) = trwfn(jngwf) * &
                          exp(i2pi*kfrac(dim))
                  end do

                  ! Calculate expectation value (inner product)
                  do kwfn=1,degen
                     psitrpsi(kwfn,jwfn) = (0.0_DP,0.0_DP)
                     do jngwf=1,ngwf_basis%num
                        dotp = (0.0_DP,0.0_DP)
                        do ingwf=1,ngwf_basis%num
                           dotp = dotp + conjg(hm_sq(ingwf,iwfn+kwfn-1)) * &
                                zov_sq(ingwf,jngwf)
                        end do
                        psitrpsi(kwfn,jwfn) = psitrpsi(kwfn,jwfn) + &
                             dotp * trwfn(jngwf)
                     end do
                  end do

               end do

               ! Calculate complex phases
               if (degen == 1) then
                  zlambda(1) = psitrpsi(1,1)
               else
                  call zgeev('N','N',degen,psitrpsi(1,1),max_degeneracy, &
                       zlambda,zphase,1,zphase,1,zwork,lwork, &
                       dwork(1,1),ierr)
                  if (ierr /= 0) then
                     write(stderr,'(a,i6)') 'Error in internal_unfold &
                          &(bandstructure_calculate): zgeev failed with code ',ierr
                     call comms_abort
                  end if
               end if

               ! Calculate unfolded k-points
               do jwfn=1,degen
                  phase = atan2(aimag(zlambda(jwfn)), &
                       real(zlambda(jwfn),kind=DP))
                  kpt_ks(iwfn+jwfn-1,dim) = -phase * recip_twopi
               end do

            end if

         end do

      end do

    end subroutine internal_unfold

  end subroutine bandstructure_calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module bandstructure
