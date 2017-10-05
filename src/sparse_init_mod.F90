! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-   !

module sparse_initialise

  use constants, only: DP, LONG, stdout

  implicit none

  public :: sparse_init_rep

  private

contains

  !===========================================================================!
  ! This subroutine performs all the initialisation routines involving the    !
  ! sparse matrices for a given set of NGWFs (eg valence, conduction, joint). !
  !---------------------------------------------------------------------------!
  ! Arguments:                                                                !
  !   elements (input) : A list of atoms in the cell                          !
  !   suffix (input)   : The suffix for this NGWF representation.             !
  !---------------------------------------------------------------------------!
  ! SPAM3 version written by Nicholas Hine, May 2009.                         !
  !---------------------------------------------------------------------------!
  ! SPAM2 version originally written by Peter Haynes, March 2004.             !
  ! Revised for distributed data, June 2006.                                  !
  ! Indexing for dense matrices added by Nicholas Hine, Dec 2007              !
  ! Modified Hamiltonian sparsity for HF exchange, Quintin Hill Feb 2008      !
  ! Modified to include conduction matrix types by Laura Ratcliff, Oct 2010   !
  ! Massive re-write by Nicholas Hine, April 2011.                            !
  !===========================================================================!

  subroutine sparse_init_rep(elements,suffix)

    use comms, only: pub_on_root, comms_abort, pub_my_node_id, &
         pub_total_num_nodes, comms_barrier
    use constants, only: NORMAL
    use ion, only : ELEMENT
    use parallel_strategy, only: parallel_strategy_list_overlaps, &
         pub_num_overlaps, pub_overlap_list, pub_first_atom_on_node, &
         pub_atoms_on_node, pub_num_atoms_on_node
    use rundat, only: pub_usehfx, kernel_cutoff, pub_ovlp_for_nonlocal, &
         pub_any_nl_proj, pub_paw, pub_aug, pub_hubbard, pub_output_detail, &
         pub_cond_calculate, cond_kernel_cutoff, pub_hfx_nlpp_for_exchange, &
         pub_hfxsw, pub_eels_calculate
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, BLKS_NGWF, BLKS_COND, BLKS_JOINT, BLKS_PROJ, &
         BLKS_HUB_PROJ, BLKS_SW, BLKS_CORE, BLKS_AUX, sparse_fill_fac_denom, &
         sparse_num_elems_on_atom, sparse_count_ss, sparse_index_ss, &
         sparse_count_union, sparse_index_union, sparse_create, &
         sparse_destroy, sparse_num_element, pub_sparse_allow_new_matrices
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
!CW
    use parallel_strategy, only: pub_elements_on_node
    use sparse, only: my_first_blk,my_last_blk
!END CW

#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Argument
    type(ELEMENT), intent(in) :: elements(:)     ! List of atoms in cell
    character(*), intent(in) :: suffix

    ! Local variables
    integer :: ierr                    ! Error flag
    integer :: iat                     ! Atom loop counter
    integer :: node                    ! Node loop counter
    integer :: iatnode                 ! Atom on node loop counter
    integer :: iblk                    ! Block counter
    integer :: nzb,my_nze,my_nzb       ! Various sizes of current matrix
    integer(kind=LONG) :: nze          ! Number of elems of current matrix
    integer(kind=LONG) :: nze_q,nze_r  ! Number of elems of q and r matrices
    integer(kind=LONG) :: nze_b,nze_c  ! Number of elems of b and c matrices
    !integer(kind=LONG) :: nze_j,nze_p  ! Number of elems of b and c matrices
    integer :: olib                    ! Library entry for basic NGWF overlap
    integer :: slib                    ! Library entry for overlap
    integer :: klib                    ! Library entry for K
    integer :: qrlib,vwlib             ! Library entries for H intermediates
    integer :: hlib                    ! Library entry for H
    integer :: kslib,kskslib           ! Library entries for kernel products
    integer :: ksksklib
    integer :: tlib, tclib
    integer :: cond_proj_ham_lib
    integer :: cond_grad_lib1, cond_grad_lib2
    integer :: blks
    character(len=10) :: pscode
    real(kind=DP) :: fill_fac          ! Denominator for filling fractions
    integer,allocatable :: seg_nzb(:)  ! Number of nonzero blocks on each seg
    integer,allocatable :: seg_nze(:)  ! Number of nonzero elements on each seg
    type(ELEMENT), allocatable :: elems_sfc(:)   ! List of atoms in SFC order
    logical :: any_proj
!CW
    integer :: i,j
!END CW

    ! Start Timer (and initially synchronise all procs)
    call comms_barrier
    call timer_clock('sparse_init',1)
#ifdef ITC_TRACE
    call VTBEGIN(vt_sparse_init,vt_err)
#endif

    pub_sparse_allow_new_matrices = .true.

    any_proj = (pub_any_nl_proj.or.pub_paw)

    ! Check arguments
    if (size(elements) /= pub_cell%nat) then
       if (pub_on_root) write(stdout,'(a)') &
            'Error in sparse_init: inconsistency in arguments.'
       call comms_abort
    end if

    ! Allocate workspace for this routine
    allocate(elems_sfc(pub_cell%nat),stat=ierr)
    call utils_alloc_check('sparse_init','elems_sfc',ierr)
    allocate(seg_nzb(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('sparse_init','seg_nzb',ierr)
    allocate(seg_nze(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('sparse_init','seg_nze',ierr)

    ! Transfer atoms to space-filling curve order
    iblk = 1
    do node=0,pub_total_num_nodes-1
       do iatnode=1,pub_num_atoms_on_node(node)
          iat = pub_atoms_on_node(iatnode,node)
          elems_sfc(iblk) = elements(iat)
          iblk = iblk + 1
       end do
    end do

    ! Sparse matrix information header
    if ((pub_on_root).and.(pub_output_detail>=NORMAL)) then
       if (suffix=='c') then
          write(stdout,'(/a)') '===================== Conduction Sparse &
               &matrix information ====================='
       else if (suffix=='j') then
          write(stdout,'(/a)') '============== Joint Valence + Conduction &
               &Sparse matrix information ============'
       else if (suffix=='a') then
          write(stdout,'(/a)') '===================== Auxiliary NGWF Sparse &
               &matrix information ====================='
       else
          write(stdout,'(/a)') '=========================== Sparse matrix &
               &information =========================='
       end if
    end if

    !---------------------------------------------------------------------------
    ! NGWF overlaps and kernel
    !---------------------------------------------------------------------------

    if (suffix=='') then
       blks = BLKS_NGWF
       pscode = 'Region'
    else if (suffix=='c') then
       blks = BLKS_COND
       pscode = 'A'
    else if (suffix=='j') then
       blks = BLKS_JOINT
       pscode = 'A'
    else if (suffix=='a') then
       blks = BLKS_AUX
       pscode = 'Region'
    end if

    call internal_create_representation

    ! Create other matrices: ks,ksk,ksks
    if (suffix/='a') call internal_create_ksksk

    !---------------------------------------------------------------------------
    ! Hamiltonian matrices
    !---------------------------------------------------------------------------

    call internal_create_hamiltonian

    ! Calculate filling factor and display
    fill_fac = sparse_fill_fac_denom(blks,blks)
    if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
         write(stdout,'(20x,a,f7.2,a)') 'Hamiltonian matrix filling:', &
         fill_fac * nze,'%'

    call internal_create_ham_strucs
    if (any_proj) call internal_create_nl_strucs

    ! Create Hubbard structures if required
    if (pub_hubbard.and.(suffix=='')) call internal_create_hubbard_strucs

    if ((suffix=='c').or.(suffix=='j').or.(suffix=='a')) then
       call internal_create_cond_structures
    end if

    !---------------------------------------------------------------------------
    ! Block diagonal matrices
    !---------------------------------------------------------------------------

    ! One-time initialisation of projector-diagonal matrices (done for valence
    ! rep only)
    if (suffix=='') then

       ! Projector diagonal matrix
       if (pub_any_nl_proj.or.pub_paw) then
          do iblk=pub_first_atom_on_node(pub_my_node_id), &
               pub_first_atom_on_node(pub_my_node_id + 1) - 1
             if (sparse_num_elems_on_atom(iblk,BLKS_PROJ)>0) then
                pub_num_overlaps(iblk) = 1
                pub_overlap_list(1,iblk) = iblk
             else
                pub_num_overlaps(iblk) = 0
             end if
          end do
          call sparse_count_ss(BLKS_PROJ,BLKS_PROJ,nze,nzb,my_nze,my_nzb, &
               seg_nze,seg_nzb)
          call sparse_index_ss(BLKS_PROJ,BLKS_PROJ,'E',nze,nzb,my_nze,my_nzb, &
               seg_nze,seg_nzb)
       end if

       ! Hubbard Projector diagonal matrix
       if (pub_hubbard) then
          do iblk=pub_first_atom_on_node(pub_my_node_id), &
               pub_first_atom_on_node(pub_my_node_id + 1) - 1
             if (sparse_num_elems_on_atom(iblk,BLKS_HUB_PROJ)>0) then
                pub_num_overlaps(iblk) = 1
                pub_overlap_list(1,iblk) = iblk
             else
                pub_num_overlaps(iblk) = 0
             end if
          end do
          call sparse_count_ss(BLKS_HUB_PROJ,BLKS_HUB_PROJ,nze,nzb,my_nze, &
               my_nzb,seg_nze,seg_nzb)
          call sparse_index_ss(BLKS_HUB_PROJ,BLKS_HUB_PROJ,'G',nze,nzb,my_nze, &
               my_nzb,seg_nze,seg_nzb)
!CW
!FULL PROJ
        do iblk=my_first_blk,my_last_blk 
          if (sparse_num_elems_on_atom(iblk,BLKS_HUB_PROJ)>0) then
             j=1
             pub_overlap_list(1,iblk) = iblk

             do i=my_first_blk,my_last_blk 
               if (i/=iblk.and.sparse_num_elems_on_atom(i,BLKS_HUB_PROJ)>0.and. &
                             & sparse_num_elems_on_atom(i,BLKS_HUB_PROJ)== &
                             & sparse_num_elems_on_atom(iblk,BLKS_HUB_PROJ)) then
                 j=j+1
                 pub_overlap_list(j,iblk) = i
               endif
             enddo
             pub_num_overlaps(iblk) = j
             write(*,*) 'hubbard atom [x] has [y] buddies : ', iblk,j
             write(*,*) 'his pals are                     : ', pub_overlap_list(1:j,iblk)
             write(*,*) 'NUMS elements on atoms           : ', sparse_num_elems_on_atom(iblk,BLKS_HUB_PROJ)
          else
             pub_num_overlaps(iblk) = 0
          end if
       end do
       call sparse_count_ss(BLKS_HUB_PROJ,BLKS_HUB_PROJ,nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
       call sparse_index_ss(BLKS_HUB_PROJ,BLKS_HUB_PROJ,'FULL',nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb)
       if(pub_my_node_id==0)then
       open(unit=7454,file='positions_cor_atoms_in_onetep')
        do iblk=my_first_blk,my_last_blk
         if (sparse_num_elems_on_atom(iblk,BLKS_HUB_PROJ)>0) then
           write(7454,'(a,3f10.4)')  trim(adjustl(pub_elements_on_node(iblk-my_first_blk+1)%symbol)),             &
                                    &             pub_elements_on_node(iblk-my_first_blk+1)%centre%x*0.529177249, &
                                    &             pub_elements_on_node(iblk-my_first_blk+1)%centre%y*0.529177249, &
                                    &             pub_elements_on_node(iblk-my_first_blk+1)%centre%z*0.529177249
         endif
        enddo
       close(7454)
       endif
!END CW

       end if


       if (pub_eels_calculate) then ! lr408
          do iblk=pub_first_atom_on_node(pub_my_node_id), &
               pub_first_atom_on_node(pub_my_node_id + 1) - 1
             if (sparse_num_elems_on_atom(iblk,BLKS_PROJ)>0 .and. &
                  sparse_num_elems_on_atom(iblk,BLKS_CORE)>0) then
                pub_num_overlaps(iblk) = 1
                pub_overlap_list(1,iblk) = iblk
             else
                pub_num_overlaps(iblk) = 0
             end if
          end do
          call sparse_count_ss(BLKS_CORE,BLKS_PROJ,nze,nzb,my_nze,my_nzb, &
               seg_nze,seg_nzb)
          call sparse_index_ss(BLKS_CORE,BLKS_PROJ,'J',nze,nzb,my_nze,my_nzb, &
               seg_nze,seg_nzb,transpose_name='P')
          call sparse_count_ss(BLKS_PROJ,BLKS_CORE,nze,nzb,my_nze,my_nzb, &
               seg_nze,seg_nzb)
          call sparse_index_ss(BLKS_PROJ,BLKS_CORE,'P',nze,nzb,my_nze,my_nzb, &
               seg_nze,seg_nzb,transpose_name='J')

       end if

    end if

    ! NGWF diagonal matrix for this rep
    fill_fac = sparse_fill_fac_denom(blks,blks)
    do iblk=pub_first_atom_on_node(pub_my_node_id), &
         pub_first_atom_on_node(pub_my_node_id + 1) - 1
       pub_num_overlaps(iblk) = 1
       pub_overlap_list(1,iblk) = iblk
    end do
    call sparse_count_ss(blks,blks,nze,nzb,my_nze,my_nzb, &
         seg_nze,seg_nzb)
    call sparse_index_ss(blks,blks,'D'//suffix,nze,nzb,my_nze, &
         my_nzb,seg_nze,seg_nzb)
    if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
         write(stdout,'(18x,a,f7.2,a)') 'NGWF diagonal matrix filling:', &
         fill_fac * nze,'%'

    !---------------------------------------------------------------------------
    ! Display filling factors of other structures
    !---------------------------------------------------------------------------

    ! Calculate filling factor denominator for fraction
    fill_fac = sparse_fill_fac_denom(blks,blks)

    ! Find the number of nonzero elements in the patterns
    ! of commonly-used matrices and print results as part of filling
    ! factor summary
    if ((pub_on_root).and.(pub_output_detail>=NORMAL).and.(suffix/='a')) then
       write(stdout,'(29x,a,f7.2,2a)') 'KS matrix filling:', &
            fill_fac * sparse_num_element(kslib),'% '
       write(stdout,'(27x,a,f7.2,2a)') 'KSKS matrix filling:', &
            fill_fac * sparse_num_element(kskslib),'% '
    end if

    ! For valence rep only, create v matrix for SWHF
    if (pub_hfxsw) call internal_create_vmatrix

    ! Sparse matrix information footer
    if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
         write(stdout,'(a/)') '============================&
         &===================================================='

    pub_sparse_allow_new_matrices = .false.

    call comms_barrier

    ! Deallocate workspace for this routine
    deallocate(seg_nze,stat=ierr)
    call utils_dealloc_check('sparse_init','seg_nze',ierr)
    deallocate(seg_nzb,stat=ierr)
    call utils_dealloc_check('sparse_init','seg_nzb',ierr)
    deallocate(elems_sfc,stat=ierr)
    call utils_dealloc_check('sparse_init','elems_sfc',ierr)

    ! Stop Timer
    call timer_clock('sparse_init',2)
#ifdef ITC_TRACE
    call VTEND(vt_sparse_init,vt_err)
#endif

  contains


    !==========================================================================!
    ! This subroutine creates the structures of an NGWF representation: the    !
    ! overlap, density kernel, and the sp- and ps- overlaps with the nonlocal  !
    ! and Hubbard projectors.                                                  !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in April 2011 based on existing code.           !
    !==========================================================================!

    subroutine internal_create_representation

      implicit none

      ! Calculate filling factor denominator for fraction
      fill_fac = sparse_fill_fac_denom(blks,blks)

      !---------------------------------------------------------------------------
      ! Overlap matrix
      !---------------------------------------------------------------------------

      ! Obtain a list of overlaps between spherical regions
      call parallel_strategy_list_overlaps(elems_sfc,pscode,pscode)

      ! Count number of corresponding nonzero NGWF matrix elements
      call sparse_count_ss(blks,blks,nze,nzb,my_nze,my_nzb, &
           seg_nze,seg_nzb)

      ! Create sparse overlap matrix
      if (.not.pub_aug) then
         call sparse_index_ss(blks,blks,'S'//suffix,nze,nzb,my_nze, &
              my_nzb,seg_nze,seg_nzb,rlib=slib)
      else
         call sparse_index_ss(blks,blks,'O'//suffix,nze,nzb,my_nze, &
              my_nzb,seg_nze,seg_nzb,rlib=olib)
      end if

      if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
           write(stdout,'(24x,a,f7.2,a)') 'Overlap matrix filling:', &
           fill_fac * nze,'%'

      ! Synchronise
      call comms_barrier

      !---------------------------------------------------------------------------
      ! Density kernel
      !---------------------------------------------------------------------------

      ! Obtain a list of overlaps of fixed cutoff
      call parallel_strategy_list_overlaps(elems_sfc,'Fixed','Fixed', &
           0.5_DP * kernel_cutoff)

      ! Count number of corresponding nonzero NGWF matrix elements
      call sparse_count_ss(blks,blks,nze,nzb,my_nze,my_nzb, &
           seg_nze,seg_nzb)

      ! Create density kernel matrix
      call sparse_index_ss(blks,blks,'K'//suffix,nze,nzb,my_nze,my_nzb, &
           seg_nze,seg_nzb,rlib=klib)

      if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
           write(stdout,'(24x,a,f7.2,a)') 'Density kernel filling:', &
           fill_fac * nze,'%'

      ! Synchronise
      call comms_barrier

      !---------------------------------------------------------------------------
      ! Projector-NGWF matrices
      !---------------------------------------------------------------------------

      ! Set up these matrices only if there are any projectors
      if (any_proj) then

         ! Obtain a list of overlaps between spherical regions and ion cores
         call parallel_strategy_list_overlaps(elems_sfc,pscode,'Core')

         ! Count number of corresponding nonzero matrix elements
         call sparse_count_ss(blks,BLKS_PROJ,nze_r,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

         ! Create Projector-NGWF overlap matrix
         call sparse_index_ss(blks,BLKS_PROJ,'R'//suffix,nze_r,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb,transpose_name='Q'//suffix)

         ! Synchronise
         call comms_barrier

         ! Obtain a list of overlaps between ion cores and spherical regions
         call parallel_strategy_list_overlaps(elems_sfc,'Core',pscode)

         ! Count number of corresponding nonzero matrix elements
         call sparse_count_ss(BLKS_PROJ,blks,nze_q,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

         ! Create NGWF-Projector overlap matrix
         call sparse_index_ss(BLKS_PROJ,blks,'Q'//suffix,nze_q,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb,transpose_name='R'//suffix)

         ! Calculate filling factor denominator for fraction
         fill_fac = sparse_fill_fac_denom(blks,BLKS_PROJ)

         if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
              write(stdout,'(14x,a,f7.2,a)') 'NGWF-Proj overlap matrix &
              &filling:', fill_fac * nze_r,'%'

         ! Synchronise
         call comms_barrier

         if (nze_q /= nze_r) then
            if ((pub_on_root).and.(pub_output_detail>=NORMAL)) then
               write(stdout,'(14x,a,f7.2,a)') 'Proj-NGWF &
                    &overlap matrix filling:', fill_fac * nze_q,'%'
               write(stdout,'(a/a/)') 'Error in sparse_init: NGWF-Projector and &
                    &Projector-NGWF overlap matrices must','have same number &
                    &of nonzero elements.'
               call comms_abort
            end if
         end if

         ! Create matrix product of Q and R matrices as nonlocal psp matrix
         call internal_create_product(qrlib,'Q'//suffix,'R'//suffix)

         ! Augmented overlap 'S' is union of direct sphere-sphere S and QR
         if (pub_aug) then

            ! Count elements in union of O and QR
            call sparse_count_union(olib,qrlib,nze,nzb,my_nze,my_nzb, &
                 seg_nze,seg_nzb)

            ! S structure becomes this union
            call sparse_index_union(olib,qrlib,'S'//suffix,nze,nzb,my_nze, &
                 my_nzb,seg_nze,seg_nzb,rlib=slib)

            ! Calculate filling factor denominator for fraction
            fill_fac = sparse_fill_fac_denom(blks,blks)

            if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
                 write(stdout,'(20x,a,f7.2,a)') 'Aug Overlap matrix filling:', &
                 fill_fac * nze,'%'
         end if

      else
         if (pub_aug) then
            if (pub_on_root) write(stdout,'(a)') &
                 'Error in sparse_init: Cannot have no projectors in augmented &
                  &formalisms.'
            call comms_abort
         end if
         qrlib = 0
      end if



      ! Set up these matrices only if there are any projectors ! lr408
      if (pub_eels_calculate) then

         !---------------------------------------------------------------------------
         ! Core-NGWF matrices
         !---------------------------------------------------------------------------

         ! Obtain a list of overlaps between spherical regions and ion cores
         call parallel_strategy_list_overlaps(elems_sfc,pscode,'Laura')

         ! Count number of corresponding nonzero matrix elements
         call sparse_count_ss(blks,BLKS_CORE,nze_c,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

         ! Create Projector-NGWF overlap matrix
         call sparse_index_ss(blks,BLKS_CORE,'C'//suffix,nze_c,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb,transpose_name='B'//suffix)

         ! Synchronise
         call comms_barrier

         ! Obtain a list of overlaps between ion cores and spherical regions
         call parallel_strategy_list_overlaps(elems_sfc,'Laura',pscode)

         ! Count number of corresponding nonzero matrix elements
         call sparse_count_ss(BLKS_CORE,blks,nze_b,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

         ! Create NGWF-Projector overlap matrix
         call sparse_index_ss(BLKS_CORE,blks,'B'//suffix,nze_b,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb,transpose_name='C'//suffix)

         ! Calculate filling factor denominator for fraction
         fill_fac = sparse_fill_fac_denom(blks,BLKS_CORE)

         if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
              write(stdout,'(14x,a,f7.2,a)') 'NGWF-Core overlap matrix &
              &filling:', fill_fac * nze_c,'%'

         ! Synchronise
         call comms_barrier

         if (nze_b /= nze_c) then
            if ((pub_on_root).and.(pub_output_detail>=NORMAL)) then
               write(stdout,'(14x,a,f7.2,a)') 'Core-NGWF &
                    &overlap matrix filling:', fill_fac * nze_b,'%'
               write(stdout,'(a/a/)') 'Error in sparse_init: NGWF-Core and &
                    &Core-NGWF overlap matrices must','have same number &
                    &of nonzero elements.'
               call comms_abort
            end if
         end if



      end if

      !---------------------------------------------------------------------------
      ! Hubbard Projector-NGWF overlaps
      !---------------------------------------------------------------------------

      ! Create W and V matrices if required
      if (pub_hubbard) then

         ! Obtain a list of overlaps between spherical regions and ion cores
         call parallel_strategy_list_overlaps(elems_sfc,pscode, 'Region')

         ! Count number of corresponding nonzero matrix elements
         call sparse_count_ss(blks,BLKS_HUB_PROJ,nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

         ! Create NGWF-Hubbard Projector overlap matrix
         call sparse_index_ss(blks,BLKS_HUB_PROJ,'W'//suffix,nze,nzb,my_nze, &
              my_nzb,seg_nze,seg_nzb,transpose_name='V'//suffix)

         ! Obtain a list of overlaps between spherical regions and ion cores
         call parallel_strategy_list_overlaps(elems_sfc,'Region', pscode)

         ! Count number of corresponding nonzero matrix elements
         call sparse_count_ss(BLKS_HUB_PROJ,blks,nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

         ! Create Hubbard Projector-NGWF overlap matrix
         call sparse_index_ss(BLKS_HUB_PROJ,blks,'V'//suffix,nze,nzb,my_nze, &
              my_nzb,seg_nze,seg_nzb,transpose_name='W'//suffix)

         ! Calculate filling factor denominator for fraction
         fill_fac = sparse_fill_fac_denom(blks,BLKS_HUB_PROJ)

         if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
              write(stdout,'(6x,a,f7.2,a)') 'NGWF-Hubbard Proj overlap matrix &
              &filling:', fill_fac * nze,'%'

         ! Create matrix product of V and W matrices as Hubbard Hamiltonian
         call internal_create_product(vwlib,'V'//suffix,'W'//suffix)

         ! Calculate filling factor denominator for fraction
         fill_fac = sparse_fill_fac_denom(blks,blks)
         if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
              write(stdout,'(29x,a,f7.2,2a)') 'VW matrix filling:', &
              fill_fac * sparse_num_element(vwlib),'% '
      else
         vwlib = 0
      end if

      !---------------------------------------------------------------------------
      ! Hubbard Projector - Projector overlaps
      !---------------------------------------------------------------------------

      ! ddor: Create Y and X matrices if required, where these are akin
      ! ddor: to W and V, respectively, but with non-local projectors
      ! ddor: replacing NGWFs.
      if (pub_hubbard .and. any_proj) then

         ! Obtain a list of overlaps between spherical regions and ion cores
         call parallel_strategy_list_overlaps(elems_sfc,'Core', 'Region')

         ! Count number of corresponding nonzero matrix elements
         call sparse_count_ss(BLKS_PROJ,BLKS_HUB_PROJ,nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

         ! Create NGWF-Hubbard Projector overlap matrix
         call sparse_index_ss(BLKS_PROJ,BLKS_HUB_PROJ,'Y'//suffix,nze,nzb,&
              my_nze, my_nzb,seg_nze,seg_nzb,transpose_name='X'//suffix)

         ! Obtain a list of overlaps between spherical regions and ion cores
         call parallel_strategy_list_overlaps(elems_sfc,'Region', 'Core')

         ! Count number of corresponding nonzero matrix elements
         call sparse_count_ss(BLKS_HUB_PROJ,BLKS_PROJ,nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

         ! Create Hubbard Projector-NGWF overlap matrix
         call sparse_index_ss(BLKS_HUB_PROJ,BLKS_PROJ,'X'//suffix,nze,nzb, &
              my_nze, my_nzb,seg_nze,seg_nzb,transpose_name='Y'//suffix)

      end if

    end subroutine internal_create_representation


    !==========================================================================!
    ! This subroutine creates the structures of the Hamiltonian, as the union  !
    ! of various structures that might contribute nonzero elements.            !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in April 2011 based on existing code.           !
    !==========================================================================!

    subroutine internal_create_hamiltonian

      implicit none

      ! Start by copying the overlap matrix into Hamiltonian
      call sparse_count_union(slib,slib,nze,nzb,my_nze,my_nzb, &
           seg_nze,seg_nzb)
      call sparse_index_union(slib,slib,'H'//suffix,nze,nzb,my_nze,my_nzb, &
           seg_nze,seg_nzb,rlib=hlib)

      if (any_proj) then

         ! Count elements in union of H and QR
         call sparse_count_union(hlib,qrlib,nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)
         ! Hamiltonian structure becomes this union
         call sparse_index_union(hlib,qrlib,'H'//suffix,nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

      end if

      if (pub_hubbard) then

         ! Count elements in union of H and VW
         call sparse_count_union(hlib,vwlib,nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)
         ! Hamiltonian structure becomes this union
         call sparse_index_union(hlib,vwlib,'H'//suffix,nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

      end if

      if (pub_usehfx .and. .not. pub_hfx_nlpp_for_exchange) then

        ! Count elements in union of H and KSKSK
        call sparse_count_union(hlib,ksksklib,nze,nzb,my_nze,my_nzb, &
             seg_nze,seg_nzb)
        ! Hamiltonian structure becomes this union
        call sparse_index_union(hlib,ksksklib,'H'//suffix,nze,nzb,my_nze, &
             my_nzb,seg_nze,seg_nzb)

      end if


    end subroutine internal_create_hamiltonian


    !==========================================================================!
    ! This subroutine creates a structure by counting overlaps for a given     !
    ! pair of region types and a pair of block lists.                          !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in April 2011 based on existing code.           !
    !==========================================================================!

    subroutine internal_create_structure(col_blks,row_blks,row_ovlp_type, &
         col_ovlp_type,name,nze,nzb,my_nze,my_nzb,seg_nze,seg_nzb,elems_sfc)

      implicit none

      ! Arguments
      integer, intent(in) :: col_blks
      integer, intent(in) :: row_blks
      character, intent(in) :: col_ovlp_type
      character, intent(in) :: row_ovlp_type
      character(*), intent(in) :: name
      integer(kind=LONG), intent(out) :: nze
      integer, intent(out) :: nzb
      integer, intent(out) :: my_nze, my_nzb
      integer, intent(out) :: seg_nze(0:pub_total_num_nodes)
      integer, intent(out) :: seg_nzb(0:pub_total_num_nodes)
      type(ELEMENT), intent(in) :: elems_sfc(pub_cell%nat)

      ! Obtain a list of overlaps between spherical regions
      call parallel_strategy_list_overlaps(elems_sfc,row_ovlp_type, &
           col_ovlp_type)

      ! Count number of corresponding nonzero NGWF matrix elements
      call sparse_count_ss(col_blks,row_blks,nze,nzb,my_nze,my_nzb, &
           seg_nze,seg_nzb)

      ! Create sparse overlap matrix
      call sparse_index_ss(col_blks,row_blks,name,nze,nzb,my_nze,my_nzb, &
           seg_nze,seg_nzb)

    end subroutine internal_create_structure


    !==========================================================================!
    ! This subroutine creates the structures of various matrices associated    !
    ! with the calculation of conduction NGWFs                                 !
    !--------------------------------------------------------------------------!
    ! Created by Nicholas Hine in April 2011, based on code by Laura Ratcliff, !
    ! from June 2010.                                                          !
    !==========================================================================!

    subroutine internal_create_cond_structures

      !------------------------------------------
      ! Cross overlap matrices with Valence NGWFs
      !------------------------------------------

      ! For auxiliary NGWFs, we only want the diagonal blocks. For Conduction or
      ! joint NGWFs, we need the full cross overlap matrix
      if (suffix=='a') then
         do iblk=pub_first_atom_on_node(pub_my_node_id), &
              pub_first_atom_on_node(pub_my_node_id + 1) - 1
            pub_num_overlaps(iblk) = 1
            pub_overlap_list(1,iblk) = iblk
         end do
      else
         ! Obtain a list of overlaps between valence NGWF region and other region
         call parallel_strategy_list_overlaps(elems_sfc,'Region',pscode)
      end if

      ! Count number of corresponding nonzero matrix elements
      call sparse_count_ss(BLKS_NGWF,blks,nze,nzb,my_nze,my_nzb, &
           seg_nze,seg_nzb)

      ! Create NGWF-COND overlap matrix
      call sparse_index_ss(BLKS_NGWF,blks,'U'//suffix,nze,nzb,my_nze, &
           my_nzb,seg_nze,seg_nzb,transpose_name='T'//suffix,rlib=tlib)

      if (suffix=='a') then
         do iblk=pub_first_atom_on_node(pub_my_node_id), &
              pub_first_atom_on_node(pub_my_node_id + 1) - 1
            pub_num_overlaps(iblk) = 1
            pub_overlap_list(1,iblk) = iblk
         end do
      else
         ! Obtain a list of overlaps between conduction NGWF region and NGWF region
         call parallel_strategy_list_overlaps(elems_sfc,pscode,'Region')
      end if

      ! Count number of corresponding nonzero matrix elements
      call sparse_count_ss(blks,BLKS_NGWF,nze,nzb,my_nze,my_nzb, &
           seg_nze,seg_nzb)

      ! Create COND-NGWF blks matrix
      call sparse_index_ss(blks,BLKS_NGWF,'T'//suffix,nze,nzb,my_nze, &
           my_nzb,seg_nze,seg_nzb,transpose_name='U'//suffix,rlib=tclib)

      ! Calculate filling factor denominator for fraction
      fill_fac = sparse_fill_fac_denom(BLKS_NGWF,blks)

      if ((pub_on_root).and.(pub_output_detail>=NORMAL)) &
           write(stdout,'(14x,a,f7.2,a)') 'NGWF-COND &
           &overlap matrix filling:', fill_fac * nze,'%'

      ! Rest is not relevant to joint matrices - only during COND optimisation
      if (suffix/='c') return

      !-------------------------------
      ! Conduction-projector overlaps
      !-------------------------------

      ! Create COND-PROJ and PROJ-COND overlaps
      if (any_proj) then

         ! Create matrix product of Q and R matrices as nonlocal psp matrix
         !call internal_create_product(tlib,'Qc','R')

         ! Create matrix product of Q and R matrices as nonlocal psp matrix
         !call internal_create_product(tlib,'Q','Rc')
      end if

      call internal_cond_ham_products

      !---------------------------------------------------------------------------
      ! Projected conduction Hamiltonian
      !---------------------------------------------------------------------------

      call sparse_count_union(cond_proj_ham_lib,hlib,nze,nzb,my_nze,my_nzb, &
           seg_nze,seg_nzb)

      call sparse_index_union(cond_proj_ham_lib,hlib,'L'//suffix,nze,nzb, &
           my_nze,my_nzb,seg_nze,seg_nzb)

      fill_fac = sparse_fill_fac_denom(BLKS_COND,BLKS_COND)

      ! Find the number of nonzero elements in the patterns
      ! of commonly-used matrices and print results as part of filling
      ! factor summary
      if ((pub_on_root).and.(pub_output_detail>=NORMAL)) then
         write(stdout,'(8x,a,f7.2,2a)') 'Proj. Cond. Hamiltonian matrix &
              &filling:',fill_fac * nze,'% '
      end if

      ! Gradient term
      call sparse_count_union(cond_grad_lib1,cond_grad_lib2,nze,nzb,my_nze, &
           my_nzb,seg_nze,seg_nzb)

      call sparse_index_union(cond_grad_lib1,cond_grad_lib2,'M'//suffix,nze, &
           nzb,my_nze,my_nzb,seg_nze,seg_nzb)

      call internal_cond_extra_products

    end subroutine internal_create_cond_structures


    !==========================================================================!
    ! This subroutine creates the filling structures of the projected          !
    ! conduction Hamiltonian and associated gradient structure                 !
    !--------------------------------------------------------------------------!
    ! Written by Laura Ratcliff, June 2010.                                    !
    !==========================================================================!

    subroutine internal_cond_ham_products

      implicit none

      ! Local variables
      type(SPAM3) :: k, h, u, t, s, m
      type(SPAM3) :: kt, kh, ks, uk
      type(SPAM3) :: ukh, uks
      type(SPAM3) :: kskt, khkt
      type(SPAM3) :: ksktm, ukhkt, ukskt, khktm

      k%structure = 'K'         ! valence kernel
      h%structure = 'H'         ! valence hamiltonian
      u%structure = 'U'//suffix ! cross overlap <\chi_\alpha|\phi_\beta>
      t%structure = 'T'//suffix ! cross overlap <\phi_\alpha|\chi_\beta>
      s%structure = 'S'         ! valence overlap
      m%structure = 'K'//suffix ! cond kernel

      ! Create the matrices and immediately destroy them to free up the data
      ! part... the structures will be kept in each case (this allows contiguous
      ! blocks of memory to remain after the creation of the structures)

      ! lr408: Products for projected conduction Hamiltonian
      call sparse_create(m)
      call sparse_destroy(m)
      call sparse_create(s)
      call sparse_destroy(s)
      call sparse_create(k)
      call sparse_destroy(k)
      call sparse_create(h)
      call sparse_destroy(h)
      call sparse_create(u)
      call sparse_destroy(u)
      call sparse_create(t)
      call sparse_destroy(t)
      call sparse_create(kt,k,t)
      call sparse_destroy(kt)
      call sparse_create(uk,u,k)
      call sparse_destroy(uk)
      call sparse_create(ukh,uk,h)
      call sparse_destroy(ukh)
      call sparse_create(uks,uk,s)
      call sparse_destroy(uks)
      call sparse_create(kh,k,h)
      call sparse_destroy(kh)
      call sparse_create(khkt,kh,kt)
      call sparse_destroy(khkt)
      call sparse_create(ukhkt,u,khkt,rlib=cond_proj_ham_lib)
      call sparse_destroy(ukhkt)
      call sparse_create(ukskt,uks,kt)
      call sparse_destroy(ukskt)

      ! lr408:  products for projected conduction Hamiltonian gradient
      call sparse_create(khktm,khkt,m,rlib=cond_grad_lib1)
      call sparse_destroy(khktm)
      call sparse_create(ks,k,s)
      call sparse_destroy(ks)
      call sparse_create(kskt,ks,kt)
      call sparse_destroy(kskt)
      call sparse_create(ksktm,kskt,m,rlib=cond_grad_lib2)
      call sparse_destroy(ksktm)

    end subroutine internal_cond_ham_products


    !==========================================================================!
    ! This subroutine creates the structures associated with products of the   !
    ! projected conduction Hamiltonian and other matrices.                     !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in November 2011.                               !
    !==========================================================================!

    subroutine internal_cond_extra_products

      implicit none

      ! Local variables
      type(SPAM3) :: k, l, t, s, kv, m
      type(SPAM3) :: kl, ks, lk, tk, ms
      type(SPAM3) :: skl, lks, lkl, klk, kvtk

      k%structure = 'K'//suffix ! cond kernel
      l%structure = 'L'//suffix ! projected cond hamiltonian
      t%structure = 'T'//suffix ! cross overlap <\phi_\alpha|\chi_\beta>
      s%structure = 'S'//suffix ! cond overlap
      m%structure = 'M'//suffix ! cond gradient term
      kv%structure = 'K'        ! valence kernel

      ! Create the matrices and immediately destroy them to free up the data
      ! part... the structures will be kept in each case (this allows contiguous
      ! blocks of memory to remain after the creation of the structures)

      ! ndmh: Products involving projected conduction Hamiltonian
      call sparse_create(k)
      call sparse_destroy(k)
      call sparse_create(l)
      call sparse_destroy(l)
      call sparse_create(k)
      call sparse_destroy(k)
      call sparse_create(s)
      call sparse_destroy(s)
      call sparse_create(kl,k,l)
      call sparse_destroy(kl)
      call sparse_create(lk,l,k)
      call sparse_destroy(lk)
      call sparse_create(ks,k,s)
      call sparse_destroy(ks)
      call sparse_create(skl,s,kl)
      call sparse_destroy(skl)
      call sparse_create(lks,l,ks)
      call sparse_destroy(lks)
      call sparse_create(lkl,l,kl)
      call sparse_destroy(lkl)
      call sparse_create(klk,k,lk)
      call sparse_destroy(klk)

      ! ndmh: Products involving gradient term m
      call sparse_create(m)
      call sparse_destroy(m)
      call sparse_create(ms,m,s)
      call sparse_destroy(ms)

      ! ndmh: other cond matrices
      call sparse_create(t)
      call sparse_destroy(t)
      call sparse_create(kv)
      call sparse_destroy(kv)
      call sparse_create(tk,t,k)
      call sparse_destroy(tk)
      call sparse_create(kvtk,kv,tk)
      call sparse_destroy(kvtk)

    end subroutine internal_cond_extra_products

    !==========================================================================!
    ! This subroutine creates the structures of the ks, ... ,ksksk matrices    !
    ! so that their filling factors can be displayed at init time              !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine, July 2008.                                     !
    ! Modified by Quintin Hill, 02/02/2009 to include LSLSL for HF exchange.   !
    ! Updated for SPAM3 by Nicholas Hine, May 2009.                            !
    ! Modified for extra memory contiguousness by Nicholas Hine in Jan 2010.   !
    ! Modified to allow reuse for conduction matrices by Laura Ratcliff        !
    ! in Oct 2010.                                                             !
    !==========================================================================!

    subroutine internal_create_ksksk

      implicit none

      ! Local variables
      type(SPAM3) :: k,s,ks,ksk,sk,sks,ksks,ksksk,kk

      k%structure = 'K'//suffix
      s%structure = 'S'//suffix

      ! Create the matrices and immediately destroy them to free up the data
      ! part... the structures will be kept in each case (this allows contiguous
      ! blocks of memory to remain after the creation of the structures)
      call sparse_create(k)
      call sparse_destroy(k)
      call sparse_create(s)
      call sparse_destroy(s)
      call sparse_create(ks,k,s,rlib=kslib)
      call sparse_destroy(ks)
      call sparse_create(ksk,ks,k)
      call sparse_destroy(ksk)
      call sparse_create(ksks,ksk,s,rlib=kskslib)
      call sparse_destroy(ksks)
      call sparse_create(sk,s,k)
      call sparse_destroy(sk)
      call sparse_create(sks,sk,s)
      call sparse_destroy(sks)
      call sparse_create(kk,k,k)
      call sparse_destroy(kk)

      ! With Hartree-Fock exchange, sparsity pattern of Hamiltonian is KSKSK
      if (pub_usehfx) then
         call sparse_create(ksksk,ksks,k,rlib=ksksklib)
         call sparse_destroy(ksksk)
       end if

    end subroutine internal_create_ksksk


    !==========================================================================!
    ! This subroutine creates the structures of all required matrices relating !
    ! to the Hamiltonian, to prevent creation of new structures at run time.   !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in Jan 2010.                                    !
    ! Modified to allow reuse for conduction matrices by Laura Ratcliff        !
    ! in Oct 2010.                                                             !
    !==========================================================================!

    subroutine internal_create_ham_strucs

      implicit none

      ! Local variables
      type(SPAM3) :: k,s,h,kh,hk,hkh,khk,skh,hks

      k%structure = 'K'//suffix
      h%structure = 'H'//suffix
      s%structure = 'S'//suffix

      ! Create the matrices and immediately destroy them to free up the data
      ! part... the structures will be kept in each case (this allows contiguous
      ! blocks of memory to remain after the creation of the structures)
      call sparse_create(k)
      call sparse_destroy(k)
      call sparse_create(s)
      call sparse_destroy(s)
      call sparse_create(h)
      call sparse_destroy(h)
      call sparse_create(kh,k,h)
      call sparse_destroy(kh)
      call sparse_create(hk,h,k)
      call sparse_destroy(hk)
      call sparse_create(hkh,h,kh)
      call sparse_destroy(hkh)
      call sparse_create(khk,kh,k)
      call sparse_destroy(khk)
      call sparse_create(skh,s,kh)
      call sparse_destroy(skh)
      call sparse_create(hks,hk,s)
      call sparse_destroy(hks)

    end subroutine internal_create_ham_strucs


    !==========================================================================!
    ! This subroutine creates the structures of all required matrices relating !
    ! to nonlocal psps, to prevent creation of new structures at run time.     !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in Jan 2010.                                    !
    ! Modified to allow reuse for conduction matrices by Laura Ratcliff        !
    ! in Oct 2010.                                                             !
    !==========================================================================!

    subroutine internal_create_nl_strucs

      implicit none

      ! Local variables
      type(SPAM3) :: r,q,k,s,rk,rks,kq

      r%structure = 'R'//suffix
      q%structure = 'Q'//suffix
      k%structure = 'K'//suffix
      s%structure = 'S'//suffix

      ! Create the matrices and immediately destroy them to free up the data
      ! part... the structures will be kept in each case (this allows contiguous
      ! blocks of memory to remain after the creation of the structures)
      call sparse_create(r)
      call sparse_destroy(r)
      call sparse_create(q)
      call sparse_destroy(q)
      call sparse_create(k)
      call sparse_destroy(k)
      call sparse_create(s)
      call sparse_destroy(s)
      call sparse_create(rk,r,k)
      call sparse_destroy(rk)
      call sparse_create(rks,rk,s)
      call sparse_destroy(rks)
      call sparse_create(kq,k,q)
      call sparse_destroy(kq)

    end subroutine internal_create_nl_strucs

    !==========================================================================!
    ! This subroutine creates the structures of all required matrices relating !
    ! to Hubbard projectors, to prevent creation of new structures at run time.!
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in Jan 2010.                                    !
    !==========================================================================!

    subroutine internal_create_hubbard_strucs

      implicit none

      ! Local variables
      type(SPAM3) :: v,w,k,s,wk,wks,kv

      v%structure = 'V'//suffix
      w%structure = 'W'//suffix
      k%structure = 'K'//suffix
      s%structure = 'S'//suffix

      ! Create the matrices and immediately destroy them to free up the data
      ! part... the structures will be kept in each case (this allows contiguous
      ! blocks of memory to remain after the creation of the structures)
      call sparse_create(v)
      call sparse_destroy(v)
      call sparse_create(w)
      call sparse_destroy(w)
      call sparse_create(k)
      call sparse_destroy(k)
      call sparse_create(s)
      call sparse_destroy(s)
      call sparse_create(wk,w,k)
      call sparse_destroy(wk)
      call sparse_create(wks,wk,s)
      call sparse_destroy(wks)
      call sparse_create(kv,k,v)
      call sparse_destroy(kv)

    end subroutine internal_create_hubbard_strucs

    !==========================================================================!
    ! This subroutine creates the structure of a matrix product from two given !
    ! structure codes. Matrix product indexing is normally done automatically  !
    ! upon creation of a new matrix with an unknown code, but for the          !
    ! creation of the Hamiltonian matrix we need the patterns of some extra    !
    ! matrices such as QR and VW, so this routine can be used to force the     !
    ! necessary call to sparse_create.                                         !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine, May 2009 as internal_create_qr.                !
    ! Generalised to any matrix product by Nicholas Hine in September 2009.    !
    ! Modified for extra memory contiguousness by Nicholas Hine in Jan 2010.   !
    !==========================================================================!

    subroutine internal_create_product(newlib,astruc,bstruc)

      implicit none

      ! Arguments
      integer, intent(out) :: newlib
      character(len=*), intent(in) :: astruc,bstruc

      ! Local variables
      type(SPAM3) :: amat,bmat,cmat

      amat%structure = astruc
      bmat%structure = bstruc

      ! Create temporary matrices and immediately destroy them to free up the
      ! memory allocate to the data before creating new structures
      call sparse_create(amat)
      call sparse_destroy(amat)
      call sparse_create(bmat)
      call sparse_destroy(bmat)

      ! Create the product and store its library entry
      call sparse_create(cmat,amat,bmat,rlib=newlib)

      ! Destroy temporary matrix
      call sparse_destroy(cmat)

      ! Synchronise
      call comms_barrier

    end subroutine internal_create_product

    !==========================================================================!
    ! This subroutine creates the structure of the vmatrix used in the         !
    ! hf_exchange module.  The vmatrix is the SW-SW two electron integral      !
    ! metric matrix.  This matrix has the same atom block sparsity pattern as  !
    ! the purified density kernel (but with different sized blocks).           !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill, October 2009.                                   !
    !==========================================================================!

    subroutine internal_create_vmatrix

      implicit none

      type(SPAM3) :: k,s,ks,ksk,ksksk
      integer :: mslib, mksksklib              ! Library indicies
      integer :: mqrlib, mrlib, mqlib          ! Library indicies

      ! qoh: SW-SW overlap matrix

      ! Obtain a list of overlaps between spherical regions
      call parallel_strategy_list_overlaps(elems_sfc,'Region','Region')

      ! qoh: Count number of corresponding nonzero SW matrix elements
      call sparse_count_ss(BLKS_SW,BLKS_SW,nze,nzb,my_nze,my_nzb, &
           seg_nze,seg_nzb)
      ! qoh: Create SW-SW overlap matrix
      call sparse_index_ss(BLKS_SW,BLKS_SW,'M',nze,nzb,my_nze,my_nzb, &
           seg_nze,seg_nzb,rlib=mslib)
      
      if (pub_hfx_nlpp_for_exchange.and.(pub_any_nl_proj.or.pub_paw).and. &
           (.not.pub_ovlp_for_nonlocal)) then

         ! Obtain a list of overlaps between spherical regions and ion cores
         call parallel_strategy_list_overlaps(elems_sfc,'Region','Core')

         ! Count number of corresponding nonzero matrix elements
         call sparse_count_ss(BLKS_SW,BLKS_PROJ,nze_r,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

         ! Create Projector-NGWF overlap matrix
         call sparse_index_ss(BLKS_SW,BLKS_PROJ,'MR',nze_r,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb,rlib=mrlib)

         ! Obtain a list of overlaps between ion cores and spherical regions
         call parallel_strategy_list_overlaps(elems_sfc,'Core','Region')

         ! Count number of corresponding nonzero matrix elements
         call sparse_count_ss(BLKS_PROJ,BLKS_SW,nze_q,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

         ! Create SW-Projector overlap matrix
         call sparse_index_ss(BLKS_PROJ,BLKS_SW,'MQ',nze_q,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb,rlib=mqlib)

         ! Create matrix product of Q and R matrices as nonlocal psp matrix
         call internal_create_product(mqrlib,'MQ','MR')

         ! Count elements in union of QR and S
         call sparse_count_union(mslib,mqrlib,nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

         ! Create Hamiltonian matrix
         call sparse_index_union(mslib,mqrlib,'M',nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

      else

         ! qoh: SW-SW two electron integral metric matrix with denskern sparsity

         ! Obtain a list of overlaps of fixed cutoff
         call parallel_strategy_list_overlaps(elems_sfc,'Fixed','Fixed', &
              0.5_DP * kernel_cutoff)

         ! qoh: Count number of corresponding nonzero SW matrix elements
         call sparse_count_ss(BLKS_SW,BLKS_SW,nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)
         ! qoh: Create SW-SW two electron integral metric
         call sparse_index_ss(BLKS_SW,BLKS_SW,'MK',nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

         k%structure = 'MK'
         s%structure = 'M'

         ! qoh: Create the matrices, and immediately destroy them so we don't
         ! qoh: waste memory on the unused data arrays.
         call sparse_create(k)
         call sparse_destroy(k)
         call sparse_create(s)
         call sparse_destroy(s)
         call sparse_create(ks,k,s)
         call sparse_destroy(ks)
         call sparse_create(ksk,ks,k)
         call sparse_destroy(ksk)
         call sparse_create(ksksk,ks,ksk,rlib=mksksklib)
         call sparse_destroy(ksksk)

         ! Count elements in union of MS and MKMSMKMSMK
         call sparse_count_union(mslib,mksksklib,nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)
         ! Create union of MS and MKMSMKMSMK as SW two electron integral metric
         call sparse_index_union(mslib,mksksklib,'M',nze,nzb,my_nze,my_nzb, &
              seg_nze,seg_nzb)

      end if
      
    end subroutine internal_create_vmatrix

  end subroutine sparse_init_rep

end module sparse_initialise
