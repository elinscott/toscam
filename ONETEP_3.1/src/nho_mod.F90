! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!               Natural Hybrid Orbitals module                   !
!                                                                !
!----------------------------------------------------------------!
! This module is an implementation of the method by              !
! J. P. Foster and F. Weinhold,                                  !
! J. Am. Chem. Soc. vol.102 (1980) p.7211                        !
! to generate *nonorthogonal* natural hybrid orbitals from       !
! NGWFs in calculations with the ONETEP program.                 !
!----------------------------------------------------------------!
! Written by Chris-Kriton Skylaris, 17/02/2006                   !
! Modifications by Nicholas Hine for SPAM3, function_basis and   !
! use of various public arrays in October 2009.                  !
!================================================================!



module nho

  use constants, only: DP
  implicit none

  private


  ! cks: defined type to hold the hybrid and reverse transform for an atom
  type ATOMIC_HYBRID

     ! cks: counter to the number of hybrids generated so far
     integer :: count

     ! cks: linear transformation from AOs to hybrids on atom
     real(kind=DP), pointer, dimension(:,:) :: mix

     ! cks: inverse transform from hybrids to AOs on atom
     real(kind=DP), pointer, dimension(:,:) :: inv_mix

  end type ATOMIC_HYBRID


  public :: nho_generate


contains 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine nho_generate(ao_to_hy, hy_to_ao, &  ! output
       ngwfs_on_grid, prev_direction_on_grid,  &  ! in-out
       kernel, overlap, elements, ngwf_basis) ! input

    !================================================================!
    ! This subroutine converts NGWFs to nonorthogonal natural hybrid !
    ! orbitals (NNHOs) and also returns the transform matrices for   !
    ! NGWF->NNHO and back.                                           !
    !----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris, 17/02/2006                   !
    !================================================================!

    use comms, only: pub_my_node_id, pub_on_root, comms_barrier
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use ion, only: element
    use parallel_strategy, only: pub_first_atom_on_node, &
           pub_num_atoms_on_node, pub_orig_atom
    use sparse, only: SPAM3, sparse_convert
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    type(SPAM3), intent(inout) :: ao_to_hy ! NGWF to hybrid transform
    type(SPAM3), intent(inout) :: hy_to_ao ! Hybrid to NGWF transform
    type(FUNC_BASIS), intent(in) :: ngwf_basis ! NGWF basis
    real(kind=DP), intent(inout) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts) ! NGWFs/Hybrids
    real(kind=DP), intent(inout) :: prev_direction_on_grid(:) ! cov search direction
    type(SPAM3), intent(in) :: kernel  ! Density kernel for hybrid generation
    type(SPAM3), intent(in) :: overlap ! Overlap matrix for hybrid generation
    type(ELEMENT), intent(in) :: elements(pub_cell%nat) ! elements of all nodes

    ! cks: <<< local variables>>>
    type(ATOMIC_HYBRID), allocatable, dimension(:) :: hybrids ! hybrids for this node
    real(kind=DP), allocatable :: full_k(:,:) ! full square density kernel
    real(kind=DP), allocatable :: full_s(:,:) ! full square overlap

    ! cks: Start timer
    call timer_clock('nho_generate',1)

    ! cks: print entry message
    if (pub_on_root) write(stdout,'(a)',advance ='no') &
         'Constructing Nonorthogonal Natural Hybrid Orbitals from current NGWFs....'

    ! cks: allocate space for the hybrids array
    call internal_hybrids_allocate

    ! cks: convert SPAM3 kernel to full square
    call sparse_convert(full_k, kernel)

    ! cks: multiply by 2.0_DP to account for spin degeneracy
    full_k =2.0_DP*full_k

    ! cks: convert SPAM3 overlap to full square
    call sparse_convert(full_s, overlap)

    ! cks: find AO -> hybrid transform for each atom on pub_my_node_id
    call nho_find_transforms(hybrids, &        
         ngwf_basis, full_k, full_s, elements)

    ! cks: apply hybrid transform to the covariant ngwf search direction
    call nho_build_hybrids(prev_direction_on_grid,    &
         hybrids, ngwf_basis, .false.)

    ! cks: apply transform on each atom to build hybrids
    call nho_build_hybrids(ngwfs_on_grid,       &
         hybrids, ngwf_basis, .true.)


    ! cks: create AO <-> hy transform matrices
    call nho_transform_matrices(ao_to_hy, hy_to_ao, &  ! output
         ngwf_basis,hybrids)


    ! cks: deallocate memory used for hybrids
    call internal_hybrids_deallocate

    ! cks: print exit message
    call comms_barrier
    if (pub_on_root) write(stdout,*)' done'

    ! cks: Stop timer
    call timer_clock('nho_generate',2)


  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_hybrids_allocate

      !==============================================!
      ! Memory allocation for the hybrids.           !
      !----------------------------------------------!
      ! Written by Chris-Kriton Skylaris, 17/02/2006.!
      ! Tidied up by Nicholas Hine, 09/10/2009.      !
      !==============================================!

      implicit none

      integer :: fc
      integer :: latom, atom
      integer :: ierr ! memory allocation error flag
      integer :: row


      allocate(hybrids(pub_num_atoms_on_node(pub_my_node_id)), stat=ierr)
      call utils_alloc_check('internal_hybrids_allocate (nho)','hybrids',ierr)

      ! cks: initialise
      hybrids(:)%count =0


      do latom =1,pub_num_atoms_on_node(pub_my_node_id)

         atom = latom + pub_first_atom_on_node(pub_my_node_id) - 1         
         fc = ngwf_basis%num_on_atom(atom)
         
         allocate(hybrids(latom)%mix(fc, fc), stat=ierr)
         call utils_alloc_check('internal_hybrids_allocate (nho)',&
              'hybrids(latom)%mix',ierr)

         allocate(hybrids(latom)%inv_mix(fc, fc), stat=ierr)
         call utils_alloc_check('internal_hybrids_allocate (nho)',&
              'hybrids(latom)%inv_mix',ierr)

         ! cks: initialise to diagonal matrix
         hybrids(latom)%mix(:,:)=0.0_DP
         ! cks: initialise to diagonal matrix
         hybrids(latom)%inv_mix(:,:)=0.0_DP
         do row=1,fc
            hybrids(latom)%mix(row, row)=1.0_DP
            hybrids(latom)%inv_mix(row, row)=1.0_DP
         enddo
      enddo

      ! cks: allocate full square density kernel
      allocate(full_k(ngwf_basis%num, ngwf_basis%num), stat=ierr)
      call utils_alloc_check('internal_hybrids_allocate (nho)',&
           'full_k',ierr)

      ! cks: allocate full square overlap matrix
      allocate(full_s(ngwf_basis%num, ngwf_basis%num), stat=ierr)
      call utils_alloc_check('internal_hybrids_allocate (nho)',&
           'full_s',ierr)

    end subroutine internal_hybrids_allocate


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_hybrids_deallocate

      !==============================================!
      ! Memory de-allocation for the hybrids.        !
      !----------------------------------------------!
      ! Written by Chris-Kriton Skylaris, 17/02/2006.!
      !==============================================!

      implicit none

      integer :: latom
      integer :: ierr



      ! cks: deallocate full square overlap matrix
      deallocate(full_s, stat=ierr)
      call utils_dealloc_check('internal_hybrids_deallocate (nho)',&
           'full_s',ierr)

      ! cks: deallocate full density kernel matrix
      deallocate(full_k, stat=ierr)
      call utils_dealloc_check('internal_hybrids_deallocate (nho)',&
           'full_k',ierr)

      ! cks: free hybrids memory
      do latom =1,pub_num_atoms_on_node(pub_my_node_id)
         deallocate(hybrids(latom)%mix, stat=ierr)
         call utils_dealloc_check('internal_hybrids_deallocate (nho)',&
              'hybrids(latom)%mix',ierr)

         deallocate(hybrids(latom)%inv_mix, stat=ierr)
         call utils_dealloc_check('internal_hybrids_deallocate (nho)',&
              'hybrids(latom)%inv_mix',ierr)
      enddo


      deallocate(hybrids, stat=ierr)
      call utils_dealloc_check('internal_hybrids_deallocate (nho)',&
           'hybrids',ierr)

    end subroutine internal_hybrids_deallocate

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  end subroutine nho_generate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nho_transform_matrices(ao_to_hy, hy_to_ao, & ! output
       ngwf_basis, hybrids)                               ! input

    !=============================================================!
    ! This subroutine generates the NGWF->Hybrid transform matrix !
    ! and the matrix of the reverse transform. These matrices     !
    ! consist simply of a diagonal of non-symmetric atomic blocks.!
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 17/02/2006.             !
    !=============================================================!

    use comms, only: comms_reduce, pub_my_node_id
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: pub_num_atoms_on_node, pub_first_atom_on_node
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_convert
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none


    type(SPAM3), intent(inout) :: ao_to_hy 
    type(SPAM3), intent(inout) :: hy_to_ao
    type(ATOMIC_HYBRID), intent(in) :: hybrids(:)
    type(FUNC_BASIS), intent(in) :: ngwf_basis

    ! cks: <<< local variables >>
    real(kind=DP), allocatable :: full_ao_to_hy(:,:)
    real(kind=DP), allocatable :: full_hy_to_ao(:,:)
    integer :: my_atoms_start
    integer :: my_atoms_end
    integer :: start_ngwf
    integer :: end_ngwf
    integer :: this_atom
    integer :: row
    integer :: col
    integer :: latom
    integer :: ierr


    ! cks: allocate full AO->hy matrix
    allocate(full_ao_to_hy(ngwf_basis%num, ngwf_basis%num), stat=ierr)
    call utils_alloc_check('nho_transform_matrices','full_ao_to_hy',ierr)

    ! cks: allocate full hy->AO matrix
    allocate(full_hy_to_ao(ngwf_basis%num, ngwf_basis%num), stat=ierr)
    call utils_alloc_check('nho_transform_matrices','full_hy_to_ao',ierr)

    ! cks: initialise
    full_hy_to_ao =0.0_DP
    full_ao_to_hy =0.0_DP

    my_atoms_start = pub_first_atom_on_node(pub_my_node_id)
    my_atoms_end = pub_first_atom_on_node(pub_my_node_id+1) - 1

    ! cks: build AO->hybrid transform matrix
    ! cks: and its inverse

    ! cks: loop over diagonal atom blocks
    latom = 0
    do this_atom=my_atoms_start,my_atoms_end
       latom = latom +1

       start_ngwf = ngwf_basis%first_on_atom(this_atom)
       end_ngwf = start_ngwf + ngwf_basis%num_on_atom(this_atom) - 1

       do row=start_ngwf,end_ngwf

          do col=start_ngwf,end_ngwf

             full_hy_to_ao(row,col) =&
                  hybrids(latom)%inv_mix(row-start_ngwf+1, col-start_ngwf+1 )

             full_ao_to_hy(row, col) =&
                  hybrids(latom)%mix(row-start_ngwf+1, col-start_ngwf+1 )        
          enddo
       enddo

    enddo


    ! cks: add the effort of each node
    call comms_reduce('SUM', full_ao_to_hy)
    call comms_reduce('SUM', full_hy_to_ao)


    ! cks: convert to sparse format
    call sparse_convert(ao_to_hy, full_ao_to_hy)
    call sparse_convert(hy_to_ao, full_hy_to_ao)



    ! cks: deallocate full hy->AO matrix
    deallocate(full_hy_to_ao, stat=ierr)
    call utils_dealloc_check('nho_transform_matrices','full_hy_to_ao',ierr)

    ! cks: deallocate full AO->hy matrix
    deallocate(full_ao_to_hy, stat=ierr)
    call utils_dealloc_check('nho_transform_matrices','full_ao_to_hy',ierr)

  end subroutine nho_transform_matrices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine nho_find_transforms(hybrids, &                             ! in-out
       ngwf_basis, full_k, full_s, elements)                            ! input

    !=============================================================!
    ! This subroutine generates the NGWF->Hybrid transform and    !
    ! its inverse for each atom and stores it in the hybrids      !
    ! defined type. For each atom the hybrids due to lone pairs   !
    ! are first extracted and the density kernel is depleted from !
    ! them after. Then hybrids due to "bonded" atom pairs are     !
    ! extracted.                                                  !
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 17/02/2006.             !
    !=============================================================!

    use comms, only: comms_barrier, pub_my_node_id
    use constants, only: DP, stdout, VERBOSE
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, geometry_distance, operator(*), operator(+)
    use ion, only: element
    use parallel_strategy, only: pub_orig_atom, pub_first_atom_on_node, &
         pub_num_atoms_on_node
    use rundat, only: pub_output_detail
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_convert
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none


    type(ATOMIC_HYBRID), intent(inout) :: hybrids(:) ! hybrid transforms
    type(FUNC_BASIS), intent(in) :: ngwf_basis ! NGWF basis
    real(kind=DP), intent(inout) :: full_k(:,:) ! full square denskern
    real(kind=DP), intent(in) :: full_s(:,:)    ! full square overlap
    type(ELEMENT),  intent(in) :: elements(pub_cell%nat)    ! in input file order

    ! cks: <<<local variables>>>
    real(kind=DP), allocatable:: mini_k(:,:) ! mini denskern for hybrid extraction
    real(kind=DP), allocatable:: mini_s(:,:) ! mini overlap for hybrid extracti
    real(kind=DP), allocatable:: hybrid_occs(:) ! mini denskern occupancies
    real(kind=DP) :: dist ! distance between atom pairs
    real(kind=DP) :: per_dist ! distance between periodic image atom pair
    real(kind=DP), parameter :: bond_threshold =5.5_DP ! threshold for chemical bonds
    real(kind=DP) :: sec_threshold ! actual bond threshold
    type(POINT) :: that_pc ! periodic centre on that atom
    integer :: latom ! local atom counter
    integer :: this_atom ! hybrid-extraction atom
    integer :: that_atom ! atom "bonded" to this_atom  
    integer :: this_orig ! in original input file order
    integer :: that_orig ! in original input file order
    integer :: this_num  ! number of NGWFs on this_atom
    integer :: that_num  ! number of NGWFs on that_atom
    integer :: block_num ! size of this-that matrix block
    integer :: ierr      ! memory allocation error flag
    integer :: row_start    ! start of row extraction
    integer :: row_midend   ! mid-end of row extraction
    integer :: row_midstart ! mid-start of row extraction
    integer :: row_end      ! mid-end of row extraction
    integer :: my_atoms_start ! global atom counting for this node
    integer :: my_atoms_end   ! global atom counting for this node
    integer :: loc1 ! periodic a1-image index
    integer :: loc2 ! periodic a2-image index
    integer :: loc3 ! periodic a3-image index
    integer :: row  ! row counter
    logical :: deplete_k    ! remove current hybrid vector from denskern
    type(SPAM3) :: k_spam3  ! spam3 of the entire density kernel


    ! cks: global limits according to space-filling atom distribution
    my_atoms_start = pub_first_atom_on_node(pub_my_node_id)
    my_atoms_end = pub_first_atom_on_node(pub_my_node_id+1) - 1


    ! cks: ========= LONE-PAIR HYBRIDS ======================================

    deplete_k =.true.
    that_atom =huge(1)
    latom =0
    ! cks: loop over all atoms on this node 
    do this_atom=my_atoms_start, my_atoms_end
       this_orig =pub_orig_atom(this_atom)
       this_num  = ngwf_basis%num_on_atom(this_atom)

       latom =latom +1

       ! cks: proceed if there is more than one NGWF on current atom
       if (this_num > 1) then

          block_num =this_num

          allocate(mini_k(block_num, block_num), stat=ierr)
          call utils_alloc_check('nho_find_transforms','mini_k',ierr)
          allocate(mini_s(block_num, block_num), stat=ierr)
          call utils_alloc_check('nho_find_transforms','mini_s',ierr)
          allocate(hybrid_occs(block_num), stat=ierr)
          call utils_alloc_check('nho_find_transforms','hybrid_occs',ierr)

          ! cks: set block extraction limits from square matrix
          row_start = ngwf_basis%first_on_atom(this_atom)
          row_midend = row_start + ngwf_basis%num_on_atom(this_atom) - 1
          row_midstart = 1
          row_end = 0


          ! cks: initialise mini_k and mini_s
          call internal_extract_blocks

          ! cks: create and store hybrids on this_atom and that_atom
          call internal_hybrid_from_block


          deallocate(hybrid_occs, stat=ierr)
          call utils_dealloc_check('nho_find_transforms','hybrid_occs',ierr)
          deallocate(mini_k, stat=ierr)
          call utils_dealloc_check('nho_find_transforms','mini_k',ierr)
          deallocate(mini_s, stat=ierr)
          call utils_dealloc_check('nho_find_transforms','mini_s',ierr)

       endif


    end do

    ! cks: ===== END LONE-PAIR HYBRIDS ======================================


    

    ! cks: +++++ copy depleted diagonal kernel blocks to all nodes ++++++++++
    call comms_barrier

   
    ! cks: create the entire kernel in spam3
    ! cks: For the diagonal blocks, only the ones
    ! cks: on atoms of pub_my_node_id are correct on the full square kernel,
    ! cks: this problem is corrected by going to SPAM3
    k_spam3%structure ='K'
    call sparse_create(k_spam3) 
    call sparse_convert(k_spam3, full_k)

    ! cks: put kernel with correct diagonal atom-atom blocks,
    ! cks: depleted from extraction of lone-pair hybrids,
    ! cks: back to full square format
    call sparse_convert(full_k, k_spam3)

    ! cks: deallocate the spam3 structure
    call sparse_destroy(k_spam3)

    call comms_barrier
    ! cks: ++ END copy depleted diagonal kernel blocks to all nodes +++++++++




    ! cks: =========== BONDING HYBRIDS ======================================
    deplete_k =.false.
    latom =0
    ! cks: loop over pairs of atoms that could be connected by a chemical bond 
    do this_atom =my_atoms_start, my_atoms_end
       this_orig =pub_orig_atom(this_atom)
       this_num = ngwf_basis%num_on_atom(this_atom)

       sec_threshold =bond_threshold
       if (elements(this_orig)%atomic_number <= 10) then
          sec_threshold =sec_threshold -1.5_DP 
       end if


       latom = latom + 1

       ! cks: proceed if there is more than one NGWF on current atom
       if (this_num > 1) then

          ! cks: find all atoms within bonding distance from this atom
          do that_atom= 1, pub_cell%nat 
             that_orig =pub_orig_atom(that_atom)
             that_num  = ngwf_basis%num_on_atom(that_atom)


             ! cks: take account of periodic boundary conditions
             dist =huge(1.0_DP)
             do loc1 = -1, 1
                do loc2 = -1, 1
                   do loc3 = -1, 1

                      that_pc =elements(that_orig)%centre &
                           +real(loc1, kind=DP)*pub_cell%a1 &
                           +real(loc2, kind=DP)*pub_cell%a2 &
                           +real(loc3, kind=DP)*pub_cell%a3

                      per_dist =geometry_distance(elements(this_orig)%centre, that_pc)

                      dist= min(dist, per_dist)

                   enddo
                enddo
             enddo


             ! cks: look for a hybrid if this and that are within bonding distance
             ! cks: but are not the same
             if ( (dist <= sec_threshold) .and. (this_atom /= that_atom) ) then

                block_num = this_num + that_num

                allocate(mini_k(block_num, block_num), stat=ierr)
                call utils_alloc_check('nho_find_transforms','mini_k',ierr)
                allocate(mini_s(block_num, block_num), stat=ierr)
                call utils_alloc_check('nho_find_transforms','mini_s',ierr)
                allocate(hybrid_occs(block_num), stat=ierr)
                call utils_alloc_check('nho_find_transforms','hybrid_occs',ierr)

                ! cks: set block extraction limits from square matrices
                if (this_atom < that_atom) then
                   row_start    = ngwf_basis%first_on_atom(this_atom)
                   row_midend   = row_start + ngwf_basis%num_on_atom(this_atom) - 1
                   row_midstart = ngwf_basis%first_on_atom(that_atom)
                   row_end      = row_midstart + ngwf_basis%num_on_atom(that_atom) - 1
                else
                   row_start    = ngwf_basis%first_on_atom(that_atom)
                   row_midend   = row_start + ngwf_basis%num_on_atom(that_atom) - 1
                   row_midstart = ngwf_basis%first_on_atom(this_atom)
                   row_end      = row_midstart + ngwf_basis%num_on_atom(this_atom) - 1
                endif

                ! cks: initialise mini_k and mini_s
                call internal_extract_blocks

                ! cks: create and store hybrids on this_atom and that_atom
                call internal_hybrid_from_block

                deallocate(hybrid_occs, stat=ierr)
                call utils_dealloc_check('nho_find_transforms','hybrid_occs',&
                     ierr)
                deallocate(mini_k, stat=ierr)
                call utils_dealloc_check('nho_find_transforms','mini_k',ierr)
                deallocate(mini_s, stat=ierr)
                call utils_dealloc_check('nho_find_transforms','mini_s',ierr)

             endif

          enddo
       end if
    end do

    ! cks: ======= END BONDING HYBRIDS ======================================





    ! cks: ======= ROBUSTNESS MEASURES ======================================
    ! cks: undo the hubrid transforms on atoms that did not manage to produce
    ! cks: as many hybrids as original NGWFs.
    
    latom = 0
    ! cks: loop over atoms of this proc 
    do this_atom =my_atoms_start, my_atoms_end
       this_orig =pub_orig_atom(this_atom)
       this_num = ngwf_basis%num_on_atom(this_atom)
       
       latom = latom + 1
       
       if ( (hybrids(latom)%count .ne. this_num) .and. (this_num > 1) ) then
          if (pub_output_detail == VERBOSE) then             
             write(stdout,'(a,i5)')'WARNING: Wrong number of hybrids on atom:', this_orig
             write(stdout,'(a,i5,a,i3,a,i3,a)') &
                  'hcount=',hybrids(latom)%count,', this_num=',this_num,&
                  'proc=',pub_my_node_id,'Resetting to NGWFs'
          endif
          ! cks: undo hybrid transform
          hybrids(latom)%mix=0.0_DP
          do row=1, this_num
             hybrids(latom)%mix(row, row)=1.0_DP
          enddo 
          
       endif
       
    enddo
    ! cks: === END ROBUSTNESS MEASURES ======================================
       

  contains


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    subroutine internal_extract_blocks

      !=============================================================!
      ! This subroutine extracts mini denskern and overlap blocks   !
      ! from the entire full square matrices.                       !
      !-------------------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 17/02/2006.             !
      !=============================================================!

      implicit none

      ! cks: local variables
      integer :: mrow
      integer :: mcol
      integer :: row
      integer :: col


      ! cks: extract left half
      mrow=0
      do row= row_start, row_midend
         mrow= mrow +1
         mcol =0
         do col= row_start, row_midend
            mcol= mcol +1
            mini_k(mrow, mcol) = full_k(row, col)
            mini_s(mrow, mcol) = full_s(row, col)
         enddo
         do col= row_midstart, row_end
            mcol= mcol +1
            mini_k(mrow, mcol) = full_k(row, col)
            mini_s(mrow, mcol) = full_s(row, col)
         enddo
      enddo


      ! cks: extract right half
      do row= row_midstart, row_end
         mrow= mrow +1
         mcol =0
         do col= row_start, row_midend
            mcol= mcol +1
            mini_k(mrow, mcol) = full_k(row, col)
            mini_s(mrow, mcol) = full_s(row, col)
         enddo
         do col= row_midstart, row_end
            mcol= mcol +1
            mini_k(mrow, mcol) = full_k(row, col)
            mini_s(mrow, mcol) = full_s(row, col)
         enddo
      enddo

    end subroutine internal_extract_blocks


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    subroutine internal_hybrid_from_block

      !=============================================================!
      ! This subroutine solves the eigenvalue problem               !
      !                            KSC=Cf                           !
      ! for the blocks extracted from the large square matrices and !
      ! stores hybrid trasforms C on an atom depending on the       ! 
      ! occupancy eigenvalue f.                                     !
      ! It also depletes the spectral expansion of the full square  ! 
      ! K from current "lone-pair" eigenvectors, upon request.      !
      !-------------------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 17/02/2006.             !
      !=============================================================!

      use wrappers, only: wrappers_dsygv_lt_2
      implicit none

      ! cks: <<< local variables >>>
      integer :: val     ! occupancy eigenvalue counter
      integer :: hcount  ! hybrid counter on current atom
      integer :: rr      ! row counter
      integer :: cc      ! column counter
      integer :: r_local ! local (block) row counter
      integer :: c_local ! local (block) column counter
      real(kind=DP), parameter ::two_els =1.75_DP ! definition of electron "pair"


      ! cks: diagonalise mini_k
      call wrappers_dsygv_lt_2(mini_k, hybrid_occs, mini_s, block_num)

      ! cks: store bonding hybrid for current atom 
      ! cks: loop over occupancies in descending order
      val_loop: do val=block_num, 1, -1


         ! cks: store hybrids if this coresponds to electron pair
         if (hybrid_occs(val) >= two_els) then

            ! cks: increase by one the index of stored hybrids on this atom
            ! cks: and extract current hybrid
            hybrids(latom)%count =hybrids(latom)%count +1
            hcount = hybrids(latom)%count

            if (hcount .le. this_num) then

               if (this_atom < that_atom) then              
                  hybrids(latom)%mix( 1: this_num, hcount) =&
                       mini_k(1: this_num, val)
               else
                  hybrids(latom)%mix( 1: this_num, hcount) =&
                       mini_k(that_num+1:  block_num, val)
               endif
            else 
               if (pub_output_detail == VERBOSE) then
                  write(stdout,'(a,i5,a,i3,a,i4,a,f6.3)')  &
                       ' hcount=',hcount,', this_num=',this_num,&
                       ' proc=',pub_my_node_id,' occ=', hybrid_occs(val)
                  write(stdout,*)'WARNING: More hybrids than NGWFs on atom', this_orig
               endif
               ! cks: undo hybrid transform if more hybrids than NGWFs
               hybrids(latom)%mix=0.0_DP
               do rr=1, this_num
                  hybrids(latom)%mix(rr,rr)=1.0_DP
               enddo
               exit val_loop
            endif


            ! cks: remove lone pair orbitals from full_k
            if (deplete_k) then

               r_local =0
               do rr =row_start, row_midend
                  r_local =r_local +1

                  c_local =0
                  do cc =row_start, row_midend
                     c_local =c_local +1

                     full_k(rr,cc) =full_k(rr,cc) &
                          -hybrid_occs(val)*mini_k(r_local, val)*mini_k(c_local, val)

                  enddo
               enddo

            endif

         endif


      enddo val_loop



    end subroutine internal_hybrid_from_block


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine nho_find_transforms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nho_build_hybrids(ngwfs_on_grid,       &
       hybrids, ngwf_basis, generate_hao_trans)

    !=============================================================!
    ! This subroutine constructs the hybrids and stores them in   !
    ! the ngwfs_on_grid array that contained the original NGWFs.  !
    ! It also constructs and stores in the hybrids type the       !
    ! reverse transform from hybrids back to the original NGWFs.  !
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 17/02/2006.             !
    !=============================================================!

    use comms, only: pub_my_node_id
    use function_basis, only: FUNC_BASIS
    use ion, only: element
    use parallel_strategy, only: pub_orig_atom, pub_num_atoms_on_node, &
         pub_first_atom_on_node
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check
    use wrappers, only: wrappers_invert_sym_matrix
    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(inout) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts) ! NGWFs -> hybrids storage
    type(ATOMIC_HYBRID), intent(inout) :: hybrids(pub_num_atoms_on_node(pub_my_node_id)) ! transforms
    logical, intent(in) :: generate_hao_trans ! generate hybrid to NGWF transformations


    ! cks: <<<local variables>>>
    real(kind=DP), allocatable :: holap(:,:) ! hybrid overlap on current atom
    real(kind=DP), allocatable :: buffer_on_grid(:) ! temporary NGWF storage
    real(kind=DP) :: bracket ! bracket integral
    integer :: this_atom ! counter for current global space-filling atom
    integer :: this_num  ! number of NGWFs on current atom
    integer :: row  ! row counter
    integer :: col  ! column counter
    integer :: nho_pt_start  ! local point start of current nho
    integer :: nho_pt_end    ! local point end of current nho
    integer :: ngwf_pt_start ! local point start of current ngwf
    integer :: ngwf_pt_end   ! local point end of current ngwf
    integer :: ierr          ! memory allocation error flag 
    integer :: start_nho     ! first nho on atom wrt local space-fill order
    integer :: end_nho       ! last nho on atom wrt local space-fill order
    integer :: my_atoms_start! first atom on this node in global space-fill count
    integer :: my_atoms_end  ! last atom on this node in global space-fill count
    integer :: latom         ! local atom counter


    my_atoms_start = pub_first_atom_on_node(pub_my_node_id)
    my_atoms_end   = pub_first_atom_on_node(pub_my_node_id+1) - 1

    allocate(buffer_on_grid(pub_cell%n_pts*ngwf_basis%n_ppds), stat=ierr)
    call utils_alloc_check('nho_build_hybrids','buffer_on_grid',ierr)

    ! cks: initialise
    buffer_on_grid(:) = ngwfs_on_grid(:)
    ngwfs_on_grid(:) = 0.0_DP


    ! cks: loop over all atoms on this node
    latom =0
    this_atom_loop: do this_atom= my_atoms_start, my_atoms_end

       this_num = ngwf_basis%num_on_atom(this_atom)
       latom =latom +1

       ! cks: first and last NWGF of this_atom wrt counting 
       ! cks: on pub_my_node_id
       start_nho = ngwf_basis%first_on_atom(this_atom) &
            - ngwf_basis%first_on_node(pub_my_node_id) + 1
       end_nho = start_nho + ngwf_basis%num_on_atom(this_atom) - 1


       ! cks: construct hybrids for this atom
       do col=0, this_num -1

          nho_pt_start = ngwf_basis%spheres(start_nho +col)%offset
          nho_pt_end = ngwf_basis%spheres(start_nho +col)%offset  &
               +pub_cell%n_pts*ngwf_basis%spheres(start_nho +col)%n_ppds_sphere -1

          ! cks: loop over NGWF "basis" to construct hybrid
          do row=start_nho,end_nho 

             ngwf_pt_start = ngwf_basis%spheres(row)%offset
             ngwf_pt_end = ngwf_basis%spheres(row)%offset  &
                  +pub_cell%n_pts*ngwf_basis%spheres(row)%n_ppds_sphere -1

             ! cks: accumulate the hybrid
             ngwfs_on_grid(nho_pt_start: nho_pt_end)= &
                  ngwfs_on_grid(nho_pt_start: nho_pt_end) &    
                  +hybrids(latom)%mix(row -start_nho +1, col+1) &
                  *buffer_on_grid(ngwf_pt_start: ngwf_pt_end)
          enddo
       enddo



       ! cks: ======== HYBRID to NGWF transforms ======================

       if (generate_hao_trans) then

          allocate(holap(this_num, this_num), stat=ierr)
          call utils_alloc_check('nho_build_hybrids','holap',ierr)

          ! cks: store inverse transform
          ! cks: loop over hybrids
          do row=0, this_num -1

             nho_pt_start = ngwf_basis%spheres(start_nho +row)%offset
             nho_pt_end = ngwf_basis%spheres(start_nho +row)%offset  &
                  +pub_cell%n_pts*ngwf_basis%spheres(start_nho +row)%n_ppds_sphere -1

             ! cks: loop over NGWFs 
             do col= start_nho, end_nho 
                ngwf_pt_start = ngwf_basis%spheres(col)%offset
                ngwf_pt_end = ngwf_basis%spheres(col)%offset  &
                     +pub_cell%n_pts*ngwf_basis%spheres(col)%n_ppds_sphere -1

                bracket = pub_cell%weight &
                     * sum(ngwfs_on_grid(nho_pt_start: nho_pt_end) &
                     * buffer_on_grid(ngwf_pt_start: ngwf_pt_end))

                ! cks: hybrid-ngwf overlap matrix on this_atom
                hybrids(latom)%inv_mix(row+1, col -start_nho +1) = bracket

                bracket=pub_cell%weight*sum(ngwfs_on_grid(nho_pt_start: nho_pt_end)&
                     *ngwfs_on_grid(ngwf_pt_start: ngwf_pt_end) )

                ! cks: hybrid-hybrid overlap matrix on this atom
                holap(row+1, col -start_nho +1)=bracket
             enddo
          enddo


          ! cks: invert hybrid overlap atomic matrix
          if (this_num > 1) then

             call wrappers_invert_sym_matrix(holap(1:this_num, 1:this_num), this_num)

             ! cks: create and store complete hybrid->NGWF transform
             hybrids(latom)%inv_mix(1:this_num, 1:this_num) =&
                  matmul(holap(1:this_num, 1:this_num), &
                  hybrids(latom)%inv_mix(1:this_num, 1:this_num)  )
          else
             hybrids(latom)%inv_mix(1,1) =hybrids(latom)%inv_mix(1,1)/holap(1,1)
          endif


          deallocate(holap, stat=ierr)
          call utils_dealloc_check('nho_build_hybrids','holap',ierr)

       end if

       ! cks: ===== END HYBRID to NGWF transforms ======================

    enddo this_atom_loop



    deallocate(buffer_on_grid, stat=ierr)
    call utils_dealloc_check('nho_build_hybrids','buffer_on_grid',ierr)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine nho_build_hybrids


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end module nho




