! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!============================================================================!
!                                                                            !
!               Electronic optimisation restart module                       !
!                                                                            !
!   The subroutines in this file were written by Chris-Kriton Skylaris       !
!   They are mainly subroutines for storing/retrieving from disk information !
!   in binary format that is used for continuining/restarting a ONETEP       !
!   calculation.                                                             !
!----------------------------------------------------------------------------!
!   Written by Chris-Kriton Skylaris, February 2004                          !
!                                                                            !
!   TCM Group, Cavendish laboratory                                          !
!   Madingley Road                                                           !
!   Cambridge CB3 0HE                                                        !
!   UK                                                                       !
!============================================================================!


module restart

  use constants, only: dp
  
  implicit none
  
  private
  
  !=== Subroutines ==========================================================!
  
  public :: restart_ngwfs_tightbox_input
  public :: restart_ngwfs_tightbox_output
  
  public :: restart_sph_waves_input
  public :: restart_sph_waves_output
  
  ! pdh: new density kernel routines to use SPAM3 format
  ! ndmh: adapted for SPAM3
  public :: restart_kernel_write
  public :: restart_kernel_read
  
  ! smmd: Subroutines and variables used to store/retrieve the denskern
  ! and ngwfs on/from memory
  public :: restart_kernel_store
  public :: restart_kernel_retrieve
  public :: restart_kernel_compose
  public :: restart_ngwfs_tightbox_store
  public :: restart_ngwfs_tightbox_compose
  public :: restart_ngwfs_tightbox_retrieve
  
  public :: restart_store_create
  public :: restart_store_reset
  public :: restart_store_destroy
  public :: restart_store_new_item
  
  !=== Type definitions for chained calculations  ===========================!
  
  type, public ::  store_item
  
  ! Label
  integer                         :: label
  
  ! Atomic coordinates
  real(kind=DP), pointer          :: coord(:,:)
  
  ! NGWFs
  integer, pointer                :: ngwfs_nodes(:)
  real(kind=DP), pointer          :: ngwfs_orig(:,:)
  real(kind=DP), pointer          :: ngwfs(:,:,:,:)
  integer                         :: tightbox_size(3)
  
  ! Density kernel
  integer, pointer                :: dk_idx(:)
  real(kind=DP), pointer          :: dk(:)
  
  end type
  
  !=== Variables used for chained calculations ==============================!
  
  type(store_item), allocatable, save  :: store(:)
  integer, save                        :: store_size
  integer, save                        :: store_nitem
  integer, allocatable, save           :: store_pointer(:)
  
  logical, public, save   :: store_tightbox_ngwfs = .false.
  logical, public, save   :: store_denskern = .false.
  logical, public, save   :: store_coordinates = .false.
  logical, public, save   :: retrieve_tightbox_ngwfs = .false.
  logical, public, save   :: retrieve_denskern = .false.
  logical, public, save   :: retrieve_coordinates = .false.
  
contains
  
!============================================================================!
!============================================================================!
!============================================================================!
  
  subroutine restart_store_create(max_items)
    
    !========================================================================!
    ! This subroutine allocates the global arrays for storage of             !
    ! coordinates, NGWFs and density kernel                                  !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 06/12/2010                            !
    !        modified by SMM Dubois on 06/04/2011                            !
    !========================================================================!
    
    use comms, only : pub_on_root
    use constants, only: DP, stdout
    use utils, only: utils_alloc_check
    
    implicit none
    
    ! Arguments
    integer, intent(in) :: max_items
    
    ! Local Variables
    integer  :: iel, ierr
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering restart_store_create'
#endif
    
    ! Define store size 
    store_size = max_items
    store_nitem = 0
    
    ! Allocate and initialise global arrays 
    allocate(store_pointer(store_size),stat=ierr)
    call utils_alloc_check('restart_store_create','store_pointer',ierr)
    store_pointer(:) = 0
    
    allocate(store(store_size),stat=ierr)
    call utils_alloc_check('restart_store_create','store',ierr)
    do iel = 1, store_size
       store(iel)%label = 0
       nullify(store(iel)%coord)
       nullify(store(iel)%ngwfs_nodes,store(iel)%ngwfs_orig,store(iel)%ngwfs)
       nullify(store(iel)%dk_idx,store(iel)%dk)
    enddo
    
    ! Store is ready to receive items
    store_coordinates = .true.
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving restart_store_create'
#endif
    
    return
  
  end subroutine restart_store_create
  
!============================================================================!
!============================================================================!
!============================================================================!
  
  subroutine restart_store_reset()
  
    !========================================================================!
    ! This subroutine deallocates the storage arrays                         !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 06/12/2010                            !
    !        modified by SMM Dubois on 06/04/2011                            !
    !========================================================================!
    
    use comms, only : pub_on_root
    use constants, only: DP, stdout
    use utils, only: utils_alloc_check, utils_dealloc_check
    
    implicit none
    
    ! Local Variables
    integer  :: iel, ierr
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering restart_store_reset'
#endif
    
    ! For each items in store, deallocate arrays and nullify pointers
    do iel = 1, min(store_size,store_nitem)
    if (associated(store(iel)%coord)) then
      deallocate(store(iel)%coord,stat=ierr)
      call utils_dealloc_check('restart_store_reset','store%coord',ierr)
      nullify(store(iel)%coord)
    endif
    
    if (associated(store(iel)%ngwfs_nodes)) then
      deallocate(store(iel)%ngwfs_nodes,stat=ierr)
      call utils_dealloc_check('restart_store_reset','store%ngwfs_nodes',ierr)
      nullify(store(iel)%ngwfs_nodes)
    endif
      
    if (associated(store(iel)%ngwfs_orig)) then
      deallocate(store(iel)%ngwfs_orig,stat=ierr)
      call utils_dealloc_check('restart_store_reset','store%ngwfs_orig',ierr)
      nullify(store(iel)%ngwfs_orig)
    endif
      
    if (associated(store(iel)%ngwfs)) then
      deallocate(store(iel)%ngwfs,stat=ierr)
      call utils_dealloc_check('restart_store_reset','store%ngwfs',ierr)
      nullify(store(iel)%ngwfs)
    endif
      
    if (associated(store(iel)%dk_idx)) then
      deallocate(store(iel)%dk_idx,stat=ierr)
      call utils_dealloc_check('restart_store_reset','store%dk_idx',ierr)
      nullify(store(iel)%dk_idx)
    endif
      
    if (associated(store(iel)%dk)) then
      deallocate(store(iel)%dk,stat=ierr)
      call utils_dealloc_check('restart_store_reset','store%dk',ierr)
      nullify(store(iel)%dk)
    endif
    enddo
    
    ! Update number of items to be coherent
    store_nitem = 0
    
    ! Store is resetted, ready to receive new items only
    store_tightbox_ngwfs = .false.
    store_denskern = .false.
    retrieve_tightbox_ngwfs = .false.
    retrieve_denskern = .false.
    retrieve_coordinates = .false.
    
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving restart_store_reset'
#endif
    
    return
  
  end subroutine restart_store_reset
  
!============================================================================!
!============================================================================!
!============================================================================!
  subroutine restart_store_destroy()
  
    !========================================================================!
    ! This subroutine deallocates the storage arrays                         !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 06/12/2010                            !
    !        modified by SMM Dubois on 06/04/2011                            !
    !========================================================================!
    
    use comms, only : pub_on_root
    use constants, only: DP, stdout
    use utils, only: utils_alloc_check, utils_dealloc_check
    
    implicit none
    
    ! Local Variables
    integer  :: ierr
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering restart_store_destroy'
#endif
    ! Deallocate all internal arrays
    call restart_store_reset() 
    
    ! Deallocate global arrays
    deallocate(store,stat=ierr)
    call utils_dealloc_check('restart_store_destroy','store',ierr)
    deallocate(store_pointer,stat=ierr)
    call utils_dealloc_check('restart_store_destroy','store_pointer',ierr)
    
    ! Update store size and number of items to be coherent
    store_size = 0
    
    ! Store is closed, no item storage/delivery 
    store_coordinates = .false.
    
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving restart_store_destroy'
#endif
    
    return
  
  end subroutine restart_store_destroy
  
!============================================================================!
!============================================================================!
!============================================================================!
  
  subroutine restart_store_new_item(coordinates,label)
  
    !========================================================================!
    ! This subroutine creates a new item in store                            !
    ! if (store_nitem .lt. store_size)                                       !
    !    item is placed in the first slot available                          !
    ! elseif (store_nitem .eq. store_size)                                   !
    !    item replace the oldest item in store                               !
    ! endif                                                                  !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 06/12/2010                            !
    !        modified by SMM Dubois on 06/04/2011                            !
    !========================================================================!
    
    use comms, only : pub_on_root
    use constants, only: DP, stdout
    use rundat, only: mix_dkn_num, mix_ngwfs_num
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check
    
    implicit none
    
    ! Arguments
    real(kind=DP), intent(in)      :: coordinates(3,pub_cell%nat)
    integer, optional, intent(out) :: label
    
    ! Local Variables
    integer  :: iel, ierr
    integer  :: storage_idx
    integer  :: storage_label
    
    !======================================================================!
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering restart_store_new_item'
#endif
    
    ! Determine storage index :
    ! If store is not full then place new item in the first slot available
    ! else replace the oldest item
    if (store_nitem .lt. store_size) then
       store_nitem = store_nitem + 1
       storage_idx = store_nitem
    else
       storage_idx = store_pointer(store_nitem)
    endif
    
    ! Update store pointer so that it lists the store items by order 
    ! of entry time (newest first)
    if (store_nitem .gt. 1 .and. store_size .gt. 1) then
       do iel = store_size, 2, -1
          store_pointer(iel) = store_pointer(iel-1)
       enddo
    endif
    store_pointer(1) = storage_idx
    
    ! If required free space for storage
    if (associated(store(storage_idx)%dk)) then
       deallocate(store(storage_idx)%dk_idx,stat=ierr)
       call utils_dealloc_check('restart_kernel_store','store%dk_idx',ierr)
       nullify(store(storage_idx)%dk_idx)
       deallocate(store(storage_idx)%dk,stat=ierr)
       nullify(store(storage_idx)%dk)
       call utils_dealloc_check('restart_kernel_store','store%dk',ierr)
    endif

    if (associated(store(storage_idx)%ngwfs_nodes)) then
       deallocate(store(storage_idx)%ngwfs,stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_store',&
         'store%ngwfs',ierr)
       nullify(store(storage_idx)%ngwfs)
       deallocate(store(storage_idx)%ngwfs_orig,stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_store',&
         'store%ngwfs_orig',ierr)
       nullify(store(storage_idx)%ngwfs_orig)
       deallocate(store(storage_idx)%ngwfs_nodes,stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_store',&
       'store%ngwfs_nodes',ierr)
       nullify(store(storage_idx)%ngwfs_nodes)
    endif

    if (associated(store(storage_idx)%coord)) then
       deallocate(store(storage_idx)%coord,stat=ierr)
       call utils_dealloc_check('restart_store_new_item','store%coord',ierr)
       nullify(store(storage_idx)%coord)
    endif
    
    ! Store coordinates
    allocate(store(storage_idx)%coord(3,pub_cell%nat),stat=ierr)
    call utils_alloc_check('restart_store_new_item','store%coord',ierr)
    store(storage_idx)%coord(:,:) = coordinates(:,:)
    
    ! Determine label given to new item
    store(storage_idx)%label = 0
    storage_label = 0
    do iel = 1, store_size
       if (store(iel)%label .gt. storage_label) storage_label = store(iel)%label
    enddo 
    store(storage_idx)%label = storage_label + 1
    if (present(label)) label = store(storage_idx)%label
   
    ! Store is ready to receive NGWFs and density kernel
    if (mix_ngwfs_num .gt. 0) store_tightbox_ngwfs = .true.
    if (mix_dkn_num .gt. 0) store_denskern = .true.
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving restart_store_new_item'
#endif
    
    return
  
  end subroutine restart_store_new_item
  
!============================================================================!
!============================================================================!

  subroutine restart_sph_waves_output(funcs_on_grid, fbasis, elements, &
     file_extension)

    !=======================================================================!
    ! Write functions on grid in spherical waves representation on a file so!
    ! calculations can be restarted.                                        !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Originally written by Alvaro Ruiz Serrano in January 2009.            !
    ! Revisited by Alvaro Ruiz Serrano in November 2011 for re-structuring  !
    ! and bug-fixing.                                                       !
    !=======================================================================!

    use constants, only: DP, stdout,NORMAL
    use comms, only: comms_abort, pub_on_root,pub_root_node_id,pub_my_node_id,&
         comms_reduce, comms_barrier, comms_free
    use function_basis, only: FUNC_BASIS, function_basis_ppds_to_sph_waves
    use geometry, only: POINT
    use ion, only : ELEMENT
    use parallel_strategy, only: pub_distr_atom, pub_node_of_atom
    use rundat, only : pub_rootname, pub_output_detail, &
         pub_write_max_l,pub_extra_n_sw
    use simulation_cell, only : pub_cell, pub_tb_recip_grid
    use spherical_wave, only: sw_init,sw_bessel_zeros, sw_bessel_zeros_init
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check,utils_unit


    implicit none

    ! ars: << arguments >>
    type(FUNC_BASIS), intent(in) :: fbasis
    real(kind =DP),   intent(in) :: funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)
    type(ELEMENT),    intent(in) :: elements(pub_cell%nat)
    character(len=*), intent(in) :: file_extension ! ddor: Extension of filename to write

    ! ars: << local variables>>
    character(len=256) :: output_file ! output file name buffer
    integer :: output_unit            ! fortran output unit number
    integer :: distr_func, atom_func  ! ngwf counters
    integer :: orig_atom              ! global atom counter in input file order
    integer :: ll, mm, nn             ! angular momentum counters
    integer :: ierr                   ! memory allocation error flag
    integer :: maxn                   ! ars: maximum number of Bessel zeros into the biggest sphere on each node
    integer :: num_swcoeff_on_file    ! ars; Number of SW coeffcients (size of the SW basis set)
    real(kind=DP), allocatable :: sw_coeffs (:,:,:) ! ars: spherical waves coefficients
    integer, allocatable, dimension(:,:) :: maxln(:,:) ! ars: maximum n for a given atom and l number
    real(kind=DP) :: alpha_init,alpha,radius ! ars: calculation of maxn


    call timer_clock("restart_sph_waves_output", 1)

    !---------------------------------------------------------------------------
    ! ars: header
    !---------------------------------------------------------------------------
    if (pub_on_root) then
       write(output_file,'(a80)') trim(pub_rootname)//"."//trim(file_extension)
       output_file = adjustl(output_file)
       
       ! ddor: Writeout depends on the filename extension.
       if (trim(file_extension) .eq. 'sw_ngwfs') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Writing NGWFs to file "', trim(output_file),'"....'
       elseif  (trim(file_extension) .eq. 'sw_hub_projs') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Writing HUBBARD PROJs to file "', trim(output_file),'"....'
       else if (trim(file_extension) .eq. 'sw_ngwfs_cond') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Writing Conduction NGWFs to file "', trim(output_file),'"....'
       else
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Writing functions to file "', trim(output_file),'"....'
       endif
    end if

    !---------------------------------------------------------------------------
    ! ars: initialisation of SW structures
    !---------------------------------------------------------------------------


    ! ars : Initialise sw module using approximate method.
    maxn = ceiling(maxval(fbasis%spheres(:)%radius) / &
         max(pub_cell%d1,pub_cell%d2,pub_cell%d3))
    call comms_reduce('MAX',maxn)
    call sw_bessel_zeros_init(2*maxn,pub_write_max_l)

    allocate(maxln(pub_cell%nat,0:pub_write_max_l), stat=ierr)
    call utils_alloc_check('restart_sph_waves_output','maxn',ierr)

    ! ars: calculate maxn for each l
    maxln = 1
    do orig_atom=1,pub_cell%nat
       alpha_init = max(pub_cell%d1,pub_cell%d2,pub_cell%d3)/&
                        elements(orig_atom)%radius

       do ll=0, pub_write_max_l
          nn = 1
          alpha = alpha_init

          do while(alpha.ge.alpha_init)
             nn = nn + 1
             ! ars : Security statement that re-initialises sw_bessel_zeros
             ! ars : if maxn_node increases in more than 2 times its
             ! ars :  previous value
             if (nn.gt.maxn+5) then
                maxn=nn
                call sw_bessel_zeros_init(2*maxn,pub_write_max_l)
             endif
             alpha = 1 - sw_bessel_zeros(nn-1,ll)/sw_bessel_zeros(nn,ll)

          enddo
          maxln(orig_atom,ll)=(nn-1)+pub_extra_n_sw
       enddo
    enddo

    ! ars : Maximum value of n within all nodes.
    maxn=maxval(maxln)

    ! ars: init SW module
    radius = maxval(fbasis%spheres(:)%radius)**2 * 2 &
         * maxval(pub_tb_recip_grid(5,:,:,:))
    call comms_reduce('MAX',radius)
    call sw_init(pub_write_max_l,maxn,radius)

    ! ars : calculate the total number of SW coefficients
    num_swcoeff_on_file = 0
    do orig_atom =1, pub_cell%nat
       do atom_func=1,fbasis%num_on_atom(pub_distr_atom(orig_atom))
          do ll =0, pub_write_max_l
             num_swcoeff_on_file = num_swcoeff_on_file + (2*ll+1)*maxln(orig_atom,ll)
          end do
       enddo
    end do

    ! ars : Check if any value of max_nl_matrix is less than or equal to zero as
    ! ars : pub_extra_n_sw can be negative
    if (any(maxln.le.0)) then
       if (pub_on_root) then
          write(stdout,'(/a17,i3)') 'pub_extra_n_sw = ',pub_extra_n_sw
          write(stdout,'(a32,i3)') 'Minimum n number = ', minval(maxln)
          write(stdout,'(a44)') 'Error: negative n number for spherical waves. ONETEP Stops'
       endif
       call comms_abort
    endif


    !---------------------------------------------------------------------------
    ! ars: end initialisation
    !---------------------------------------------------------------------------


    !---------------------------------------------------------------------------
    ! ars: write on disk
    !---------------------------------------------------------------------------

    output_unit =utils_unit()
    if (pub_on_root) then

       open(unit=output_unit, form="unformatted" ,file=trim(output_file), &
            action="write" )

       write(output_unit,err=140) fbasis%num          ! Number of NGWFs
       write(output_unit,err=140) num_swcoeff_on_file ! Number of SW coeffs
       write(output_unit,err=140) pub_write_max_l     ! Max_l
       write(output_unit,err=140) maxn                ! Max_n (global)
       do orig_atom=1,pub_cell%nat                           ! Element, Atom,l,max_n
          do ll=0, pub_write_max_l
             write(output_unit, err=140) &
                elements(orig_atom)%species_id,elements(orig_atom)%symbol,&
                orig_atom,ll,maxln(orig_atom,ll)
          enddo
       enddo
    endif

    ! ars : allocate sw coeffficient
    allocate(&
      sw_coeffs(0:pub_write_max_l,-pub_write_max_l:pub_write_max_l,1:maxn),&
      stat=ierr)
    call utils_alloc_check('restart_sph_waves_output','sw_coeffs',ierr)

    ! cks: Loop over all atoms in the order they appear in the input file
    do orig_atom =1, pub_cell%nat

       ! ars: loop over functions on atom
       do atom_func= 1,fbasis%num_on_atom(pub_distr_atom(orig_atom))

          distr_func = fbasis%first_on_atom(pub_distr_atom(orig_atom)) + &
               atom_func - 1

          ! cks: this essentially converts the non-blocking sends of
          ! cks: comms_send to blocking sends. Since parallel performance
          ! cks: is not the issue here, this makes the code simpler.
          call comms_barrier

          call function_basis_ppds_to_sph_waves(sw_coeffs, &
               maxln(orig_atom,:), maxn, distr_func, pub_root_node_id, &
               funcs_on_grid, fbasis)


          ! &&&&&&&&&&& WRITE THE BUFFER &&&&&&&&&&&&&&&&&&&&&&&
          if (pub_on_root) then
             do ll = 0, pub_write_max_l
                do mm = -ll,ll
                   do nn = 1, maxln(orig_atom,ll)

                      ! Atom,NGWF,l,m,n,swcoeff
                      write(output_unit,err=140) &
                           elements(orig_atom)%species_id,&
                           elements(orig_atom)%symbol,&
                           orig_atom, atom_func,ll,mm,nn,sw_coeffs(ll,mm,nn)
                      
                   end do
                end do
             end do
          end if
          ! &&&&&&& END WRITE THE BUFFER &&&&&&&&&&&&&&&&&&&&&&&

       end do
    end do

    ! cks: root node closes output file
    if (pub_on_root) then
       close(unit =output_unit)

       if (pub_output_detail >= NORMAL ) write(stdout,*)' done'

    endif


    ! cks: don't go away just yet! Wait for all communications to finish!
    ! cks: also wait before deallocating uni_tbox which is used in non-blocking
    ! cks: send!
    call comms_barrier
    ! cks: clean up send stack too
    call comms_free

    ! ars: deallocate and finish
    deallocate(sw_coeffs, stat=ierr)
    call utils_dealloc_check('restart_sph_waves_output','sw_coeffs',ierr)
    deallocate(maxln, stat=ierr)
    call utils_dealloc_check('restart_sph_waves_output','max_n',ierr)


    call timer_clock("restart_sph_waves_output", 2)

    return

140 if (pub_on_root) write(stdout,*)'Problem writing to file in&
                      & restart_sph_waves_output. ONETEP stops'
    call comms_abort


  end subroutine restart_sph_waves_output


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine restart_sph_waves_input(funcs_on_grid,fbasis,file_extension)


    !=======================================================================!
    ! Read functions on grid in spherical waves representation on a file so !
    ! calculations can be restarted.                                        !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Originally written by Alvaro Ruiz Serrano in January 2009.            !
    ! Revisited by Alvaro Ruiz Serrano in November 2011 for re-structuring  !
    ! and bug-fixing.                                                       !
    !=======================================================================!


    use constants, only: DP, stdout, NORMAL
    use comms, only: comms_bcast, comms_barrier,comms_abort, comms_free, &
         pub_root_node_id, pub_my_node_id, pub_on_root
    use basis, only : SPHERE, FUNCTION_TIGHT_BOX
    use function_basis, only: FUNC_BASIS,function_basis_sph_waves_to_ppds
    use parallel_strategy, only: pub_distr_atom, pub_node_of_atom
    use rundat, only : pub_rootname, pub_read_max_l, pub_output_detail
    use simulation_cell, only : pub_cell, pub_fftbox, pub_tb_recip_grid
    use spherical_wave, only: sw_init,sw_bessel_zeros, sw_recp_generate_in_tb
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check,utils_unit


    implicit none

    ! ars: Arguments
    type(FUNC_BASIS), intent(in   ) :: fbasis
    real(kind=DP),    intent(inout) :: funcs_on_grid(fbasis%n_ppds*pub_cell%n_pts)
    character(len=*), intent(in   ) :: file_extension ! ddor: Extension of filename to write

    ! ars: Variables
    character(len=256) :: read_file  ! read file name buffer
    integer :: input_unit   ! fortran input unit number
    integer :: ierr         ! memory allocation error flag
    integer :: local_ngwf   ! ngwf counter on this proc
    integer :: orig_atom    ! global atom counter in input file order
    integer :: atom_ngwf,distr_ngwf ! ngwf of current atom counter
    integer :: ll,mm,nn,ii
    integer :: num_funcs_on_file    ! total number of ngwfs read from file
    integer :: maxl                 ! maximum value of l within all funcs
    integer :: maxn                 ! maximum value of n within all funcs
    character(len=2) :: read_atom_symbol
    character(len=4) :: read_species_id
    integer :: read_atom,read_ngwf,read_ll,read_mm,read_nn

    integer :: num_swcoeffs, num_swcoeffs_perfunc,num_swcoeffs_on_file
    integer :: skip_lines ! skip a number of lines if required
    real(kind=DP) :: max_radius
    integer, allocatable :: maxln(:,:)    ! ars: n zeros array in serial
    real(kind=DP), allocatable :: sw_coeffs(:,:,:) ! ars: spherical waves coefficients


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Initialisations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    call timer_clock("restart_sph_waves_input", 1)

    !---------------------------------------------------------------------------
    ! ars: header
    !---------------------------------------------------------------------------
    if (pub_on_root) then
       write(read_file,'(a80)') trim(pub_rootname)//"."//trim(file_extension)
       read_file = adjustl(read_file)
       
       ! ddor: Writeout depends on the filename extension.
       if (trim(file_extension) .eq. 'sw_ngwfs') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Reading NGWFs from file "', trim(read_file),'"....'
       elseif  (trim(file_extension) .eq. 'sw_hub_projs') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Reading HUBBARD PROJs from file "', trim(read_file),'"....'
       else if (trim(file_extension) .eq. 'sw_ngwfs_cond') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Reading Conduction NGWFs from file "', trim(read_file),'"....'
       else
          if (pub_output_detail >= NORMAL ) write(stdout,'(/,3a)',advance ='no') &
               'Reading functions from file "', trim(read_file),'"....'
       endif
    end if

    ! ars : radius to initialise sw module.
    ! ars : this must be the same than the used in restart_ngwfs_swcoeff_output

    max_radius = maxval(fbasis%spheres(:)%radius)**2 * 2 &
         * maxval(pub_tb_recip_grid(5,:,:,:))

    ! ars : initialise checking variables
    num_swcoeffs = 0
    num_swcoeffs_perfunc = 0

    ! ars : this is the array that comes out from this subroutine
    funcs_on_grid(:) = 0.0_DP

    !---------------------------------------------------------------------------
    ! ars: read headers from disk and safe-check
    !---------------------------------------------------------------------------

    ! cks: get a unit number that is free
    input_unit =utils_unit()

    ! cks: root node opens input file, reads basic info
    if (pub_on_root) then

       open(unit =input_unit, form="unformatted" ,file=trim(read_file), &
            action="read", status='old')


       ! ars : read information on file and broadcast
       read(input_unit,err =150) num_funcs_on_file    ! Number of funcs saved on file
       read(input_unit,err =150) num_swcoeffs_on_file ! Number of SW coeffs
       read(input_unit,err =150) maxl                 ! Max_l
       read(input_unit,err =150) maxn                 ! Max_n

    endif
    call comms_bcast(pub_root_node_id, num_funcs_on_file)
    call comms_bcast(pub_root_node_id, num_swcoeffs_on_file)
    call comms_bcast(pub_root_node_id, maxl)
    call comms_bcast(pub_root_node_id, maxn)


    ! allocate max_n, read and broadcast
    allocate(maxln(pub_cell%nat,0:maxl),stat=ierr)
    call utils_alloc_check('restart_sph_waves_input','maxln',ierr)
    maxln=0.0_DP
    if (pub_on_root) then
       ! ars: read from file and fill max_nl_matrix_serial
       do orig_atom=1,pub_cell%nat
          do ll=0, maxl
             read(input_unit,err=150) read_species_id, &
                read_atom_symbol,read_atom,read_ll,maxln(read_atom,read_ll)
          enddo
       enddo
    endif
    call comms_bcast(pub_root_node_id, maxln)

    ! ars: check if the information has been read correctly
    if(num_funcs_on_file.ne.fbasis%num) then
       if (pub_on_root) then
          write(stdout,*)'Number of functions on file = ', num_funcs_on_file, &
               ', Number of functions on basis set = ',fbasis%num
          write(stdout,*)'The number of basis set functions is different &
            &to the number of functions on file'
          write(stdout,*)'Problem with reading spherical waves file. ONETEP stops.'
       end if
       call comms_abort
    endif

    if (pub_read_max_l.gt.maxl) then
       if (pub_on_root) then
          write(stdout,*)'max l on file = ', maxl
          write(stdout,*)'read_max_l = ', pub_read_max_l
          write(stdout,*)'The angular momentum in the input file is higher &
               &than in the file'
          write(stdout,*)'Maximum angular momentum has been reset to: ',&
               maxl
       endif
       pub_read_max_l = maxl
    endif

    ! ars : initialise spherical_wave module
    call sw_init(maxl,maxn+1,max_radius)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End initialisations ~~~~~~~~~~~~~~~~~~~~~~~~~~!



    !---------------------------------------------------------------------------
    ! ars: init SW module and read coeffs from disk
    !---------------------------------------------------------------------------

    ! ars : allocate arrays
    allocate(sw_coeffs(0:pub_read_max_l,-pub_read_max_l:pub_read_max_l,maxn),stat=ierr)
    call utils_alloc_check ('restart_sph_waves_input','sw_coeffs',ierr)
    sw_coeffs = 0.0_DP

    ! cks: Loop over all atoms in the order they appear in the input file
    do orig_atom =1, pub_cell%nat

       ! ars: loop over functions on atom
       do atom_ngwf =1,fbasis%num_on_atom(pub_distr_atom(orig_atom))

          distr_ngwf = fbasis%first_on_atom(pub_distr_atom(orig_atom)) + &
               atom_ngwf - 1


          ! cks: this essentially converts the non-blocking sends of
          ! cks: comms_send to blocking sends. Since parallel performance
          ! cks: is not the issue here, this makes the code simpler.
          call comms_barrier

          ! ars : n_swcoeffs_ngwf becomes zero after the ngwf loop
          num_swcoeffs_perfunc=0
          ! ars : if pub_read_max_l < read_max_l then the program "skip_liness" extra lines
          skip_lines=0
          do ll = 0,maxl
             skip_lines=skip_lines+maxln(orig_atom,ll)*(2*ll+1)
          enddo

          ! ars : loops over l, m, n
          do ll = 0,pub_read_max_l
             do mm = -ll,ll
                do nn = 1, maxln(orig_atom,ll)


                   ! &&&&&&&&&&& READ INTO BUFFER &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

                   if (pub_on_root) then
                      read(input_unit, err=150)&
                           read_species_id,read_atom_symbol,read_atom,read_ngwf,&
                           read_ll,read_mm,read_nn,sw_coeffs(read_ll,read_mm,read_nn)

                      num_swcoeffs = num_swcoeffs + 1
                      num_swcoeffs_perfunc = num_swcoeffs_perfunc + 1
                   endif
                   ! &&&&&&& END READ INTO BUFFER &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


                enddo  !!! End Loop 5 !!!
             enddo  !!! End Loop 4 !!!
          enddo  !!! End Loop 3 !!!

          if (pub_on_root) then
             do ii=1,skip_lines-num_swcoeffs_perfunc
                read(input_unit, err=150) ! ars: read nothing
             enddo
          end if

          ! ars: create functions in PPD representation
          call function_basis_sph_waves_to_ppds(fbasis,funcs_on_grid,sw_coeffs,&
               maxln(orig_atom,:),pub_read_max_l,maxn,distr_ngwf,pub_root_node_id)

       enddo  !!! End Loop 2 !!!
    enddo  !!! End Loop 1 !!!


    ! cks: root node closes output file
    if (pub_on_root) then
       close(unit=input_unit)

       ! ars : print info about the reading process
       if (pub_output_detail.gt.NORMAL) then
          write(stdout,'(a/)')' done'
          write(stdout,'(a,13x,i0)') "Maximum l of the SW read = ", pub_read_max_l
          write(stdout,'(a,17x,i0)') "Number of SW coefficients read = ",&
               num_swcoeffs
          write(stdout,'(a,i0)') "Number of SW coefficients in the file = ",&
               num_swcoeffs_on_file
       end if
    endif

    if (num_swcoeffs > num_swcoeffs_on_file) then
       write(stdout,'(a)')'The number of coefficients read is bigger than&
            & the number of coefficients in the file. ONETEP Stops'
       call comms_abort
    endif


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Deallocate arrays ~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    deallocate(maxln, stat=ierr)
    call utils_dealloc_check('restart_sph_waves_input','maxln',ierr)
    deallocate(sw_coeffs,stat=ierr)
    call utils_dealloc_check('restart_sph_waves_input','swcoeff',ierr)
    !~~~~~~~~~~~~~~~~~~~~~~~~~ End deallocate arrays ~~~~~~~~~~~~~~~~~~~~~~~~~~~!


    ! cks: close comms
    call comms_barrier
    call comms_free


    call timer_clock("restart_sph_waves_input", 2)

    return

150 if (pub_on_root) write(stdout,*)&
         'Problem reading from file in restart_sph_waves_input. ONETEP stops'
    call comms_abort

  end subroutine restart_sph_waves_input


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine restart_ngwfs_tightbox_output(ngwfs_on_grid, ngwf_basis, &
       elements, file_extension)

    !=======================================================================!
    ! This subroutine outputs the current NGWFs in "universal tightbox"     !
    ! representation to a .tightbox_ngwfs file. The order of the NGWFs is   !
    ! the one in which they appear in the input file and is not affected    !
    ! by the use of the space-filling curve. Also the .tightbox_ngwfs       !
    ! remains usable if the positions of the atoms change slightly          !
    ! after the restart.                                                    !
    ! ndmh 18/12/2009: In cases where the FFTbox coincides with the cell    !
    ! but the individual tightboxes do not need to, a new temporary NGWF    !
    ! basis is constructed, for which the FFTbox only coincides with the    !
    ! cell along direction i if L_i < 2*R_NGWF + \delta_i.                  !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 12/3/2004.                        !
    ! Modified by Chris-Kriton Skylaris on 13/4/2004 so that it also outputs!
    ! the coordinates of the centre of each NGWF in (real) numbers of grid  !
    ! points from the origin of its tightbox.                               !
    ! Modified by D. O'Regan in November 2009 for different file extensions.!
    ! Modified by Nicholas Hine in December 2009, to allow reduced file     !
    ! sizes by only reading in tightboxes even if FFTbox coincides with     !
    ! cell in some directions.                                              !
    !=======================================================================!

    use comms, only: comms_abort, comms_reduce, comms_barrier, comms_free, &
         pub_my_node_id, pub_on_root, pub_root_node_id
    use constants, only: DP, stdout, NORMAL
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
         function_basis_deallocate, function_basis_distribute, &
         function_basis_copy_spheres, function_basis_init_tight_boxes, &
         function_basis_ppds_to_tightbox
    use geometry, only: magnitude
    use ion, only : ELEMENT
    use parallel_strategy, only: pub_distr_atom, pub_node_of_atom, &
         pub_elements_on_node
    use rundat, only : pub_rootname, pub_output_detail
    use simulation_cell, only : pub_cell, pub_fftbox
    use timer, only: timer_clock
    use utils, only : utils_alloc_check, utils_dealloc_check, utils_unit

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind =DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds *pub_cell%n_pts) ! NGWFs on this proc
    type(ELEMENT), intent(in) :: elements(pub_cell%nat) ! elements of all proc (in input file order)
    character(len=*), intent(in) :: file_extension ! ddor: Extension of filename to write

    ! Local variables
    character(len=256) :: output_file  ! outputr file name buffer
    real(kind=DP), allocatable, dimension(:,:,:) :: uni_tbox  ! universal tightbox
    real(kind=DP) :: tb_orig1   ! a1-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: tb_orig2   ! a2-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: tb_orig3   ! a3-position of NGWF wrt tightbox in (real number of) grid points
    integer :: maxt_n1     ! a1-maximum tightbox points
    integer :: maxt_n2     ! a2-maximum tightbox points
    integer :: maxt_n3     ! a3-maximum tightbox points
    integer :: orig_atom   ! global atom counter in input file order
    integer :: atom_ngwf   ! ngwf of current atom counter
    integer :: output_unit ! fortran output unit number
    integer :: distr_ngwf  ! global ngwf counter in spacefil order
    real(kind=DP) :: sendbuf(3) ! comms send pack buffer
    real(kind=DP) :: recvbuf(3) ! comms recv pack buffer
    real(kind=DP) :: max_radius ! Maximum NGWF radius
    integer :: ierr        ! memory allocation error flag
    type(FUNC_BASIS) :: tmp_basis ! ndmh: for dealing with FFTbox coinciding with cell
    logical :: orig_coin1, orig_coin2, orig_coin3 ! ndmh: whether the fftbox
                           ! originally coincided with cell in 1,2,3 directions

    call timer_clock("restart_ngwfs_tightbox_output", 1)

    ! cks: find maximum number of points for universal tightbox that has odd
    ! cks: number of points in each dimension
    maxt_n1 = ngwf_basis%maxtight_pts1
    maxt_n2 = ngwf_basis%maxtight_pts2
    maxt_n3 = ngwf_basis%maxtight_pts3

    ! ndmh: if the FFTbox coincides with the cell in any direction, but the
    ! ndmh: tightbox does not necessarily need to, we can define a
    ! ndmh: temporary NGWF basis where the tightboxes are as small as they
    ! ndmh: can be, so that we do not get simulation-cell-sized tightboxes
    ! ndmh: written out to files unnecessarily
    orig_coin1 = pub_fftbox%coin1
    orig_coin2 = pub_fftbox%coin2
    orig_coin3 = pub_fftbox%coin3
    if (orig_coin1.or.orig_coin2.or.orig_coin3) then

       max_radius = maxval(ngwf_basis%spheres(:)%radius)
       call comms_reduce('MAX',max_radius)

       ! ndmh: temporarily override coin1,2,3 variables if possible
       if (magnitude(pub_fftbox%a1) > 2.0_DP*max_radius + pub_fftbox%d1) &
            pub_fftbox%coin1 = .false.
       if (magnitude(pub_fftbox%a2) > 2.0_DP*max_radius + pub_fftbox%d2) &
            pub_fftbox%coin2 = .false.
       if (magnitude(pub_fftbox%a3) > 2.0_DP*max_radius + pub_fftbox%d3) &
            pub_fftbox%coin3 = .false.

    end if

    ! ndmh: if any of the coin1,2,3 variables have been changed...
    if ((orig_coin1.neqv.pub_fftbox%coin1).or. &
         (orig_coin2.neqv.pub_fftbox%coin2).or. &
         (orig_coin3.neqv.pub_fftbox%coin3)) then

       ! ndmh: create a temporary function basis where tightboxes
       ! ndmh: do not coincide with cell if possible
       call function_basis_allocate(tmp_basis,ngwf_basis%num, &
            'tmp_'//ngwf_basis%name)
       call function_basis_distribute(tmp_basis,elements)
       call function_basis_copy_spheres(tmp_basis,ngwf_basis)
       call function_basis_init_tight_boxes(tmp_basis)

       ! ndmh: set a universal tightbox size which depends on the actual size
       ! ndmh: of the temporary NGWF tightboxes
       maxt_n1 = maxval(tmp_basis%tight_boxes(:)%tight_pts1)
       maxt_n2 = maxval(tmp_basis%tight_boxes(:)%tight_pts2)
       maxt_n3 = maxval(tmp_basis%tight_boxes(:)%tight_pts3)
       call comms_reduce('MAX',maxt_n1)
       call comms_reduce('MAX',maxt_n2)
       call comms_reduce('MAX',maxt_n3)

    end if

    ! cks: allocate universal tightbox buffer
    allocate(uni_tbox(maxt_n1,maxt_n2,maxt_n3),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_output','uni_tbox',ierr)

    ! cks: get a unit number that is free
    output_unit = utils_unit()

    ! cks: root node allocates buffer, opens output file, write basic info
    if (pub_on_root) then

       !  write(stdout,*)'output_unit=', output_unit

       write(output_file,'(a80)') trim(pub_rootname)//"."//trim(file_extension)
       output_file = adjustl(output_file)

       ! ddor: Writeout depends on the filename extension.
       if (trim(file_extension) .eq. 'tightbox_ngwfs') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(3a)',advance ='no') &
               'Writing NGWFs to file "', trim(output_file),'"....'
       elseif  (trim(file_extension) .eq. 'tightbox_hub_projs') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(3a)',advance ='no') &
               'Writing HUBBARD PROJs to file "', trim(output_file),'"....'
       else if (trim(file_extension) .eq. 'tightbox_ngwfs_cond') then
          if (pub_output_detail >= NORMAL ) write(stdout,'(3a)',advance ='no') &
               'Writing Conduction NGWFs to file "', trim(output_file),'"....'
       endif

       open(unit=output_unit, form="unformatted" ,file=trim(output_file), &
            action="write" )

       write(output_unit, err =140) ngwf_basis%num
       write(output_unit, err =140) maxt_n1
       write(output_unit, err =140) maxt_n2
       write(output_unit, err =140) maxt_n3

    end if

    ! cks: Loop over all atoms in the order they appear in the input file
    atom_loop: do orig_atom =1, pub_cell%nat


       ngwfs_on_atom_loop: do atom_ngwf= &
            1,ngwf_basis%num_on_atom(pub_distr_atom(orig_atom))


          distr_ngwf = ngwf_basis%first_on_atom(pub_distr_atom(orig_atom)) + &
               atom_ngwf - 1

          ! cks: this essentially converts the non-blocking sends of
          ! cks: comms_send to blocking sends. Since parallel performance
          ! cks: is not the issue here, this makes the code simpler.
          call comms_barrier

          if ((orig_coin1.neqv.pub_fftbox%coin1).or. &
            (orig_coin2.neqv.pub_fftbox%coin2).or. &
            (orig_coin3.neqv.pub_fftbox%coin3)) then

             call function_basis_ppds_to_tightbox(uni_tbox, &
                  tb_orig1,tb_orig2,tb_orig3,distr_ngwf,pub_root_node_id, &
                  maxt_n1,maxt_n2,maxt_n3,ngwfs_on_grid,tmp_basis, &
                  sendbuf,recvbuf)
          else
             call function_basis_ppds_to_tightbox(uni_tbox, &
                  tb_orig1,tb_orig2,tb_orig3,distr_ngwf,pub_root_node_id, &
                  maxt_n1,maxt_n2,maxt_n3,ngwfs_on_grid,ngwf_basis, &
                  sendbuf,recvbuf)
          end if


          ! &&&&&&&&&&& WRITE THE BUFFER &&&&&&&&&&&&&&&&&&&&&&&
          if (pub_on_root) then

             write(output_unit, err=140) tb_orig1, tb_orig2, tb_orig3

             write(output_unit, err=140) uni_tbox(1:maxt_n1,1:maxt_n2,1:maxt_n3)

          end if
          ! &&&&&&& END WRITE THE BUFFER &&&&&&&&&&&&&&&&&&&&&&&


       enddo ngwfs_on_atom_loop
    enddo atom_loop




    ! cks: root node closes output file
    if (pub_on_root) then
       close(unit =output_unit)

       if (pub_output_detail >= NORMAL ) write(stdout,*)' done'

    endif


    ! cks: don't go away just yet! Wait for all communications to finish!
    ! cks: also wait before deallocating uni_tbox which is used in non-blocking
    ! cks: send!
    call comms_barrier
    ! cks: clean up send stack too
    call comms_free


    ! cks: deallocate universal tightbox buffer
    deallocate(uni_tbox,stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_output','uni_tbox',ierr)


    ! ndmh: Clean up temporary basis if we created one, and reset coin1,2,3
    if ((orig_coin1.neqv.pub_fftbox%coin1).or. &
         (orig_coin2.neqv.pub_fftbox%coin2).or. &
         (orig_coin3.neqv.pub_fftbox%coin3)) then

       call function_basis_deallocate(tmp_basis)
       pub_fftbox%coin1 = orig_coin1
       pub_fftbox%coin2 = orig_coin2
       pub_fftbox%coin3 = orig_coin3

    end if

    call timer_clock("restart_ngwfs_tightbox_output", 2)

    return

140 if (pub_on_root) write(stdout,*)'Problem writing to file in &
         &restart_ngwfs_tightbox_output. ONETEP stops'
    call comms_abort


  end subroutine restart_ngwfs_tightbox_output




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine restart_ngwfs_tightbox_input(ngwfs_on_grid, ngwf_basis, &
       elements, file_extension)

    !=======================================================================!
    ! This subroutine reads the current NGWFs in "universal tightbox"       !
    ! representation from a .tightbox_ngwfs file. The order of the NGWFs    !
    ! in the .tightbox_ngwfs file is the one in which they appear in the    !
    ! input file and is not affected by the use of the space-filling curve. !
    ! Most importantly the "universal tightbox" storage allows              !
    ! the .tightbox_ngwfs file of a previous atomic                         !
    ! configuration to be read, i.e. one can initialise the NGWFs with      !
    ! the converged NGWFs of a previous calculation where the atoms         !
    ! were in different positions.                                          !
    ! ndmh 18/12/2009: In cases where the FFTbox coincides with the cell    !
    ! but the individual tightboxes do not need to, a new temporary NGWF    !
    ! basis is constructed, for which the FFTbox only coincides with the    !
    ! cell along direction i if L_i < 2*R_NGWF + \delta_i.                  !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 12/3/2004                         !
    ! Modified by Chris-Kriton Skylaris on 13/4/2004 so that the read NGWFs !
    ! are shifted to the centres of their corresponding atoms.              !
    ! Fix by Chris-Kriton Skylaris on 16/02/2005 of bug related to          !
    ! initialisation of tb_start1/2/3 reported by Victor Milman.            !
    ! Modified by D. O'Regan in November 2009 for different file extensions.!
    ! Modified by Nicholas Hine in December 2009, to allow reduced file     !
    ! sizes by only reading in tightboxes even if FFTbox coincides with     !
    ! cell in some directions.                                              !
    !=======================================================================!

    use basis, only: basis_ket_start_wrt_fftbox
    use comms, only: comms_abort, comms_barrier, comms_bcast, comms_free, &
         comms_recv, comms_reduce, comms_send, pub_my_node_id, pub_on_root, &
         pub_root_node_id
    use constants, only: DP, stdout
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
         function_basis_deallocate, function_basis_distribute, &
         function_basis_copy_spheres, function_basis_init_tight_boxes, &
         function_basis_tightbox_to_ppds
    use geometry, only: magnitude
    use ion, only : ELEMENT
    use parallel_strategy, only: pub_distr_atom, pub_node_of_atom, &
         pub_elements_on_node
    use rundat, only: pub_rootname
    use simulation_cell, only : pub_cell, pub_fftbox
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(out) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts) ! NGWFs on this proc
    type(ELEMENT), intent(in) :: elements(pub_cell%nat) ! elements of all proc (in input file order)
    character(len=*), intent(in) :: file_extension ! ddor: Extension of filename to read

    ! Local Variables
    character(len=256) :: read_file  ! read file name buffer
    real(kind=DP), allocatable, dimension(:,:,:) :: uni_tbox  ! universal tightbox
    complex(kind=DP), allocatable, dimension(:,:,:) :: fftbox_complex ! complex fftbox space
    complex(kind=DP), allocatable, dimension(:,:,:) :: fftbox_complex_shifted ! complex fftbox space
    real(kind=DP), allocatable, dimension(:,:,:) :: fftbox_buffer     ! real fftbox space
    real(kind=DP) :: max_radius          ! globaly maximum ngwf radius
    real(kind=DP) :: read_tb_orig1 ! read a1-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: read_tb_orig2 ! read a2-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: read_tb_orig3 ! read a3-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: read_weight   ! read grid point weight
    integer :: n1, n2, n3, ld1, ld2 ! fftbox dimensions
    integer :: orig_atom    ! global atom counter in input file order
    integer :: atom_ngwf    ! ngwf of current atom counter
    integer :: input_unit   ! fortran input unit number
    integer :: read_num     ! total number of ngwfs read from file
    integer :: maxt_n1      ! a1-maximum tightbox points
    integer :: maxt_n2      ! a2-maximum tightbox points
    integer :: maxt_n3      ! a3-maximum tightbox points
    integer :: read_maxt_n1 ! a1-maximum tightbox points read from file
    integer :: read_maxt_n2 ! a2-maximum tightbox points read from file
    integer :: read_maxt_n3 ! a3-maximum tightbox points read from file
    integer :: tb_start1    ! a1-grid point start of universal tightbox wrt fftbox
    integer :: tb_start2    ! a2-grid point start of universal tightbox wrt fftbox
    integer :: tb_start3    ! a3-grid point start of universal tightbox wrt fftbox
    integer :: distr_ngwf   ! global ngwf counter in spacefill order
    real(kind=DP) :: sendbuf(3) ! comms send pack buffer
    real(kind=DP) :: recvbuf(3) ! comms recv pack buffer
    integer :: ierr         ! memory allocation error flag

    ! ndmh:
    type(FUNC_BASIS) :: tmp_basis
    logical :: orig_coin1,orig_coin2,orig_coin3 ! ndmh: whether the fftbox
                           ! originally coincided with cell in 1,2,3 directions

    call timer_clock("restart_ngwfs_tightbox_input", 1)


    ! cks: ********** INITIALISATIONS ***************************************
    ngwfs_on_grid = 0.0_DP

    ! cks: find maximum number of (odd numbers) points for universal tightbox
    max_radius = maxval(ngwf_basis%spheres(:)%radius)
    call comms_reduce('MAX',max_radius)
    maxt_n1 = int(2.0_DP*max_radius / pub_cell%d1) + 1
    maxt_n1 = maxt_n1 + 1 - mod(maxt_n1,2)
    maxt_n2 = int(2.0_DP*max_radius / pub_cell%d2) + 1
    maxt_n2 = maxt_n2 + 1 - mod(maxt_n2,2)
    maxt_n3 = int(2.0_DP*max_radius / pub_cell%d3) + 1
    maxt_n3 = maxt_n3 + 1 - mod(maxt_n3,2)

    n1 = pub_fftbox%total_pt1
    n2 = pub_fftbox%total_pt2
    n3 = pub_fftbox%total_pt3
    ld1 = pub_fftbox%total_ld1
    ld2 = pub_fftbox%total_ld2

    read_maxt_n1 =0
    read_maxt_n2 =0
    read_maxt_n3 =0
    ! cks: ****** END INITIALISATIONS ***************************************

    ! cks: get a unit number that is free
    input_unit = utils_unit()

    ! cks: root node opens input file, reads basic info
    if (pub_on_root) then

       ! write(stdout,*)'input_unit=', input_unit

       write(read_file,'(a80)') trim(pub_rootname)//"."//trim(file_extension)
       read_file = adjustl(read_file)

       write(stdout,*)' '
       ! ddor: Writeout depends on the filename extension.
       if (trim(file_extension) .eq. 'tightbox_ngwfs') then
          write(stdout,'(3a)',advance ='no') &
               'Reading NGWFs from file "', trim(read_file),'" ...'
       elseif (trim(file_extension) .eq. 'tightbox_hub_projs') then
          write(stdout,'(3a)',advance ='no') &
               'Reading HUBBARD PROJs from file "', trim(read_file),'" ...'
       elseif (trim(file_extension) .eq. 'tightbox_ngwfs_cond') then
          write(stdout,'(3a)',advance ='no') &
               'Reading conduction NGWFs from file "', trim(read_file),'" ...'
       endif

       open(unit =input_unit, form="unformatted" ,file=trim(read_file), &
            action="read", status='old' )
       read(input_unit, err =150) read_num
       read(input_unit, err =150) read_maxt_n1
       read(input_unit, err =150) read_maxt_n2
       read(input_unit, err =150) read_maxt_n3

       ! cks: check that the correct number of NGWFs will be read from file
       if (read_num.ne.ngwf_basis%num) then

          write(stdout,*)'Problem with reading ',trim(file_extension),&
               &' file. ONETEP stops.'
          write(stdout,*)'read_num=', read_num, ' ngwf_basis%num=', ngwf_basis%num
          write(stdout,*)'read_maxt_n1=', read_maxt_n1, ' maxt_n1=', maxt_n1
          write(stdout,*)'read_maxt_n2=', read_maxt_n2, ' maxt_n2=', maxt_n2
          write(stdout,*)'read_maxt_n3=', read_maxt_n3, ' maxt_n3=', maxt_n3

          call comms_abort
       endif

    endif

    ! ndmh: broadcast the read_maxt's from the root to all other nodes, since
    ! ndmh: they are needed to determine whether to initialise temporary basis
    call comms_bcast(pub_root_node_id, read_maxt_n1)
    call comms_bcast(pub_root_node_id, read_maxt_n2)
    call comms_bcast(pub_root_node_id, read_maxt_n3)

    ! ndmh: if the FFTbox coincides with the cell in any direction, but the
    ! ndmh: tightbox does not necessarily need to, we may need to define a
    ! ndmh: temporary NGWF basis where the tightboxes are as small as they
    ! ndmh: can be, unless the file we are reading in was written with
    ! ndmh: tightboxes coinciding with cell
    orig_coin1 = pub_fftbox%coin1
    orig_coin2 = pub_fftbox%coin2
    orig_coin3 = pub_fftbox%coin3
    if (orig_coin1.or.orig_coin2.or.orig_coin3) then

       ! ndmh: temporarily override coin1,2,3 variables if possible (as long
       ! ndmh: as the written values are not for tightbox==fftbox)
       if ((magnitude(pub_fftbox%a1) > 2.0_DP*max_radius + &
            pub_fftbox%d1*pub_cell%n_pt1).and. &
            (read_maxt_n1 /= pub_fftbox%total_pt1)) pub_fftbox%coin1 = .false.
       if ((magnitude(pub_fftbox%a2) > 2.0_DP*max_radius + &
            pub_fftbox%d2*pub_cell%n_pt1).and. &
            (read_maxt_n2 /= pub_fftbox%total_pt2)) pub_fftbox%coin2 = .false.
       if ((magnitude(pub_fftbox%a3) > 2.0_DP*max_radius + &
            pub_fftbox%d3*pub_cell%n_pt1).and. &
            (read_maxt_n3 /= pub_fftbox%total_pt3)) pub_fftbox%coin3 = .false.

#ifdef DEBUG
       if ((orig_coin1.neqv.pub_fftbox%coin1).and.pub_on_root) write(stdout,*) &
            'DEBUG: Temporarily overriding pub_fftbox%coin1'
       if ((orig_coin2.neqv.pub_fftbox%coin2).and.pub_on_root) write(stdout,*) &
            'DEBUG: Temporarily overriding pub_fftbox%coin2'
       if ((orig_coin3.neqv.pub_fftbox%coin3).and.pub_on_root) write(stdout,*) &
            'DEBUG: Temporarily overriding pub_fftbox%coin3'
#endif
    end if

    ! ndmh: if any of the coin1,2,3 variables have been changed...
    if ((orig_coin1.neqv.pub_fftbox%coin1).or. &
         (orig_coin2.neqv.pub_fftbox%coin2).or. &
         (orig_coin3.neqv.pub_fftbox%coin3)) then

       ! ndmh: create a temporary function basis where tightboxes
       ! ndmh: do not coincide with cell if possible
       call function_basis_allocate(tmp_basis,ngwf_basis%num, &
            'tmp_'//ngwf_basis%name)
       call function_basis_distribute(tmp_basis,elements)
       call function_basis_copy_spheres(tmp_basis, ngwf_basis)
       call function_basis_init_tight_boxes(tmp_basis)

       ! ndmh: set a universal tightbox size which depends on the actual size
       ! ndmh: of the temporary NGWF tightboxes
       maxt_n1 = maxval(tmp_basis%tight_boxes(:)%tight_pts1)
       maxt_n2 = maxval(tmp_basis%tight_boxes(:)%tight_pts2)
       maxt_n3 = maxval(tmp_basis%tight_boxes(:)%tight_pts3)
       call comms_reduce('MAX',maxt_n1)
       call comms_reduce('MAX',maxt_n2)
       call comms_reduce('MAX',maxt_n3)

    end if

    ! cks: NGWFs will be placed in the centre of the fftbox
    ! cks: or left where they are in the simulation cell
    ! cks: (when the simulation cell and FFT-box coincide)
    call basis_ket_start_wrt_fftbox(tb_start1, tb_start2, tb_start3, n1, n2, n3)

    ! ndmh: check tightbox will not extend out of FFTbox if the coin flag has
    ! ndmh: been overridden. If so, try to adjust tb_start to compensate. If
    ! ndmh: that fails, print an error and quit.
    if ((tb_start1 + maxt_n1 > n1).and.(.not.pub_fftbox%coin1)) then
#ifdef DEBUG
       if (pub_on_root) write(stdout,*) 'DEBUG: Attempting to adjust tb_start1'
#endif
       if (n1 - maxt_n1 > 0) then
          tb_start1 = (n1 - maxt_n1)/2
       else
          write(stdout,'(a,i4,a)') 'Error in restart_ngwfs_tightbox_input: &
               &Tightbox size ',maxt_n1, 'along 1-axis'
          write(stdout,'(a,i4)') 'is not compatible with FFT-box size ',n1
          call comms_abort
       end if
    end if
    if ((tb_start2 + maxt_n2 > n2).and.(.not.pub_fftbox%coin2)) then
#ifdef DEBUG
       if (pub_on_root) write(stdout,*) 'DEBUG: Attempting to adjust tb_start2'
#endif
       if (n2 - maxt_n2 > 0) then
          tb_start2 = (n2 - maxt_n2)/2
       else
          write(stdout,'(a,i4,a)') 'Error in restart_ngwfs_tightbox_input: &
               &Tightbox size ',maxt_n2, 'along 2-axis'
          write(stdout,'(a,i4)') 'is not compatible with FFT-box size ',n2
          call comms_abort
       end if
    end if
    if ((tb_start3 + maxt_n3 > n3).and.(.not.pub_fftbox%coin3)) then
#ifdef DEBUG
       if (pub_on_root) write(stdout,*) 'DEBUG: Attempting to adjust tb_start3'
#endif
       if (n3 - maxt_n3 > 0) then
          tb_start3 = (n3 - maxt_n3)/2
       else
          write(stdout,'(a,i4,a)') 'Error in restart_ngwfs_tightbox_input: &
               &Tightbox size ',maxt_n3, 'along 3-axis'
          write(stdout,'(a,i4)') 'is not compatible with FFT-box size ',n3
          call comms_abort
       end if
    end if

    ! cks: allocate universal tightbox buffer
    maxt_n1 =max(min(maxt_n1,n1), read_maxt_n1)
    maxt_n2 =max(min(maxt_n2,n2), read_maxt_n2)
    maxt_n3 =max(min(maxt_n3,n3), read_maxt_n3)

    ! cks: send the new maxt's from the root to all other nodes
    call comms_bcast(pub_root_node_id, maxt_n1)
    call comms_bcast(pub_root_node_id, maxt_n2)
    call comms_bcast(pub_root_node_id, maxt_n3)

    ! cks: send read_weight to all nodes
    call comms_bcast(pub_root_node_id, read_weight)

    ! ndmh: allocate temporary storage
    allocate(uni_tbox(maxt_n1,maxt_n2,maxt_n3),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_input','uni_tbox',ierr)
    allocate(fftbox_complex(ld1,ld2,n3),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_input','fftbox_complex',ierr)
    allocate(fftbox_complex_shifted(ld1,ld2,n3),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_input', &
         'fftbox_complex_shifted',ierr)
    allocate(fftbox_buffer(ld1,ld2,n3),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_input','fftbox_buffer',ierr)

    ! cks: Loop over all atoms in the order they appear in the input file
    atom_loop: do orig_atom =1, pub_cell%nat

       ngwfs_on_atom_loop: do atom_ngwf=1, &
            ngwf_basis%num_on_atom(pub_distr_atom(orig_atom))

          distr_ngwf = ngwf_basis%first_on_atom(pub_distr_atom(orig_atom)) + &
               atom_ngwf - 1

          ! cks: this essentially converts the non-blocking sends of
          ! cks: comms_send to blocking sends. Since parallel performance
          ! cks: is not the issue here, this makes the code simpler.
          call comms_barrier

          ! cks: initialise
          uni_tbox =0.0_DP

          ! &&&&&&&&&&& READ INTO BUFFER &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          if (pub_on_root) then

             read(input_unit, err=150) read_tb_orig1, read_tb_orig2, &
                  read_tb_orig3

             read(input_unit, err=150) uni_tbox(1:read_maxt_n1, &
                  1:read_maxt_n2, 1:read_maxt_n3)

          endif
          ! &&&&&&& END READ INTO BUFFER &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

          if ((orig_coin1.neqv.pub_fftbox%coin1).or. &
            (orig_coin2.neqv.pub_fftbox%coin2).or. &
            (orig_coin3.neqv.pub_fftbox%coin3)) then

             call function_basis_tightbox_to_ppds(uni_tbox, &
                  read_tb_orig1,read_tb_orig2,read_tb_orig3, &
                  distr_ngwf,pub_root_node_id,maxt_n1,maxt_n2,maxt_n3, &
                  ngwfs_on_grid,tmp_basis,tb_start1,tb_start2,tb_start3, &
                  fftbox_complex,fftbox_complex_shifted,fftbox_buffer, &
                  sendbuf,recvbuf)

          else

             call function_basis_tightbox_to_ppds(uni_tbox, &
                  read_tb_orig1,read_tb_orig2,read_tb_orig3, &
                  distr_ngwf,pub_root_node_id,maxt_n1,maxt_n2,maxt_n3, &
                  ngwfs_on_grid,ngwf_basis,tb_start1,tb_start2,tb_start3, &
                  fftbox_complex,fftbox_complex_shifted,fftbox_buffer, &
                  sendbuf,recvbuf)

          end if

       enddo ngwfs_on_atom_loop
    enddo atom_loop

    ! cks: root node closes output file
    if (pub_on_root) then
       close(unit=input_unit)
       write(stdout,*)' done'
    endif

    ! cks: don't go away just yet! Wait for all communications to finish!
    call comms_barrier
    ! cks: clean up send stack too
    call comms_free

    deallocate(fftbox_buffer,stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_input', &
         'fftbox_buffer',ierr)
    deallocate(fftbox_complex_shifted,stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_input', &
         'fftbox_complex_shifted',ierr)
    deallocate(fftbox_complex,stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_input', &
         'fftbox_complex',ierr)
    deallocate(uni_tbox,stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_input', &
         'uni_tbox',ierr)

    ! ndmh: Clean up temporary basis if we created one, and reset coin1,2,3
    if ((orig_coin1.neqv.pub_fftbox%coin1).or. &
         (orig_coin2.neqv.pub_fftbox%coin2).or. &
         (orig_coin3.neqv.pub_fftbox%coin3)) then

       call function_basis_deallocate(tmp_basis)
       pub_fftbox%coin1 = orig_coin1
       pub_fftbox%coin2 = orig_coin2
       pub_fftbox%coin3 = orig_coin3

    end if


    call timer_clock("restart_ngwfs_tightbox_input", 2)

    return

150 if (pub_on_root) write(stdout,*)'Problem reading from file in restart_ngwfs_tightbox_input. ONETEP stops'
    call comms_abort


  end subroutine restart_ngwfs_tightbox_input


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine restart_kernel_read(denskern,read_cond,read_dmft, read_xc, &
      read_hartree)

    !======================================================================!
    ! This subroutine reads in the density kernel from an external binary  !
    ! file.                                                                !
    ! Its operation is not affected by changes in:                         !
    !  1) internal ordering of atoms                                       !
    !  2) number of processors                                             !
    !  3) density kernel cutoff                                            !
    !  4) atomic positions                                                 !
    !----------------------------------------------------------------------!
    ! Written by Peter Haynes 15/2/2007                                    !
    ! Modified by Nicholas Hine 11/08/2008                                 !
    !======================================================================!

    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout
    use rundat, only : pub_rootname
    use sparse, only: SPAM3, sparse_read
    use timer, only: timer_clock

    implicit none

    ! Argument
    type(SPAM3), intent(inout) :: denskern(:)
    logical, optional, intent(in) :: read_cond

    ! Local variables
    character(len=80) :: filename
    logical :: fileexists
!CW
    logical,optional :: read_dmft, read_xc, read_hartree
!END CW
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering restart_kernel_read'
#endif

    call timer_clock("restart_kernel_read", 1)

    ! Establish filename from input file
    ! lr408: Check if we are reading in a conduction kernel
    if (present(read_cond)) then
       if (read_cond) then
          write(filename,'(2a)') trim(pub_rootname),'.dkn_cond'
       else
          write(filename,'(2a)') trim(pub_rootname),'.dkn'
       end if
!CW
    elseif(present(read_dmft))then
       if (read_dmft) then
          write(filename,'(2a)') trim(pub_rootname),'.dkn_dmft'
       else
          write(filename,'(2a)') trim(pub_rootname),'.dkn'
       end if
    elseif(present(read_xc))then
       if (read_xc) then
          write(filename,'(2a)') trim(pub_rootname),'.xc'
       else
          write(filename,'(2a)') trim(pub_rootname),'.dkn'
       end if
    elseif(present(read_hartree))then
       if (read_hartree) then
          write(filename,'(2a)') trim(pub_rootname),'.hartree'
       else
          write(filename,'(2a)') trim(pub_rootname),'.dkn'
       end if

!END CW
    else
       write(filename,'(2a)') trim(pub_rootname),'.dkn'
    end if

    if (pub_on_root) write(stdout,'(/3a)',advance='no') &
         'Reading density kernel from file "', trim(filename),'" ...'

    ! Check that the file exists
    ! ndmh: only root node needs to be able to see the file
    if (pub_on_root) then
       inquire(file=filename,exist=fileexists)
    else
       fileexists = .true.
    end if

    if (fileexists) then

       ! Read density kernel from this file
       call sparse_read(denskern,trim(filename))

    else

       if (pub_on_root) write(stdout,'(/a/)') ' File not found, quitting'
       call comms_abort

    end if

    if (pub_on_root) write(stdout,'(a)') ' done'

    call timer_clock("restart_kernel_read", 2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving restart_kernel_read'
#endif

  end subroutine restart_kernel_read



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!CW
  subroutine restart_kernel_write(denskern, write_cond, &
       & write_hartree, write_xc, write_dmft)
!END CW
    !======================================================================!
    ! This subroutine writes the density kernel to an external binary      !
    ! file.                                                                !
    ! Its operation is not affected by changes in:                         !
    !  1) internal ordering of atoms                                       !
    !  2) number of processors                                             !
    !  3) density kernel cutoff                                            !
    !  4) atomic positions                                                 !
    !----------------------------------------------------------------------!
    ! Written by Peter Haynes 15/2/2007                                    !
    !======================================================================!

    use comms, only : pub_on_root
    use constants, only: DP,stdout, normal
    use rundat, only : pub_rootname, pub_output_detail
    use sparse, only: SPAM3, sparse_write
    use timer, only: timer_clock

    implicit none

    ! Argument
    type(SPAM3), intent(inout) :: denskern(:)
    logical, optional, intent(in) :: write_cond
    logical, optional, intent(in) :: write_hartree ! ddor-28feb17
    logical, optional, intent(in) :: write_xc ! ddor-28feb17
!CW
    logical, optional, intent(in) :: write_dmft 
!END CW
    ! Local variables
    character(len=80) :: filename

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering restart_kernel_write'
#endif

    call timer_clock("restart_kernel_write", 1)

    ! Establish filename from input file
    ! lr408: Check if we are writing out a conduction kernel
    if (present(write_cond)) then
       if (write_cond) then
          write(filename,'(2a)') trim(pub_rootname),'.dkn_cond'
       else
          write(filename,'(2a)') trim(pub_rootname),'.dkn'
       end if
!    else
!       write(filename,'(2a)') trim(pub_rootname),'.dkn'
!    end if
!CW
     elseif (present(write_xc)) then
         if (write_xc) then
          write(filename,'(2a)') trim(pub_rootname),'.xc'
         else
          write(filename,'(2a)') trim(pub_rootname),'.dkn'
         end if
     elseif (present(write_hartree)) then
         if (write_hartree) then
          write(filename,'(2a)') trim(pub_rootname),'.hartree'
         else
          write(filename,'(2a)') trim(pub_rootname),'.dkn'
         end if
     elseif (present(write_dmft)) then
       if (write_dmft) then
          write(filename,'(2a)') trim(pub_rootname),'.dkn_dmft'
       else
          write(filename,'(2a)') trim(pub_rootname),'.dkn'
       end if
    else
          write(filename,'(2a)') trim(pub_rootname),'.dkn'
    end if
!END CW

    if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
       write(stdout,'(/3a)',advance='no') &
            'Writing density kernel to file "', trim(filename),'" ...'
    end if

    ! Write density kernel to this file
    call sparse_write(denskern,trim(filename))

    if (pub_on_root .and. (pub_output_detail >= NORMAL )) then
       write(stdout,'(a)') ' done'
    end if

    call timer_clock("restart_kernel_write", 2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving restart_kernel_write'
#endif

  end subroutine restart_kernel_write

  
!============================================================================!
!============================================================================!
!============================================================================!
  
  
  subroutine restart_kernel_store(denskern,label)
  
    !======================================================================!
    ! This subroutine store the density kernel for chained calculation     !
    !    (input)  density kernel to be stored [SPAM3]                      !
    !    (output) density kernel stored in unsegmented format              !
    ! If (label == 0 or not present)                                       !
    !   store kernel in most recent item of the store                      !
    ! elseif (label=-i < 0)                                                !
    !   store kernel in the (i+1)th item of the store                      !
    ! elseif (label=i > 0)                                                 !
    !   store kernel in item labelled i                                    !
    ! endif                                                                !
    !----------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois 25/2/2010                              !
    !        modified by SMM Dubois (02/12/2010)                           !
    !        modified by SMM Dubois (06/04/2011)                           !
    !======================================================================!
    
    
    use comms, only : pub_on_root, pub_my_node_id
    use constants, only: DP, stdout
    use timer, only: timer_clock
    use rundat, only: mix_dkn_num
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
    use sparse, only: SPAM3, sparse_convert_unsegment_real
    
    implicit none
    
    ! Argument
    type(SPAM3), intent(in) :: denskern(:)     ! Density kernel to be stored
    integer, optional, intent(inout) :: label
    
    ! Internal variables
    integer  :: storage_label 
    integer  :: dklen, idxlen, storage_idx
    integer  :: ierr
    
    !======================================================================!
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering restart_kernel_store'
#endif
    
    call timer_clock("restart_kernel_store", 1)
    
    if (.not. present(label)) then
       storage_label = 0
    else
       storage_label = label
    endif 
    
    ! Find storage index
    call restart_storage_index(storage_idx,storage_label)
    
    ! Return label
    if (present(label)) label = store(storage_idx)%label
    if (pub_on_root) write(stdout,'(a,i5,a)') &
       'restart_kernel_store : storage_idx  (', storage_idx, ')'
    
    ! smmd: If required free space for storage
    if (associated(store(storage_idx)%dk)) then
       deallocate(store(storage_idx)%dk_idx,stat=ierr)
       call utils_dealloc_check('restart_kernel_store','store%dk_idx',ierr)
       nullify(store(storage_idx)%dk_idx)
       deallocate(store(storage_idx)%dk,stat=ierr)
       nullify(store(storage_idx)%dk)
       call utils_dealloc_check('restart_kernel_store','store%dk',ierr)
    endif
    
    ! smmd: allocate the storage arrays
    call sparse_convert_unsegment_real(denskern,idxlen,store(storage_idx)%dk_idx, &
       dklen,store(storage_idx)%dk,.true.)
    
    allocate(store(storage_idx)%dk_idx(idxlen),stat=ierr)
    call utils_alloc_check('restart_kernel_store','store%dk_idx',ierr)
    allocate(store(storage_idx)%dk(dklen),stat=ierr)
    call utils_alloc_check('restart_kernel_store','store%dk',ierr)
    
    ! smmd: store denskern and update storage index
    call sparse_convert_unsegment_real(denskern,idxlen,store(storage_idx)%dk_idx, &
       dklen,store(storage_idx)%dk,.false.)
    
    ! Store contains at least one density kernel
    retrieve_denskern = .true.
    
    call timer_clock("restart_kernel_store", 2)
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving restart_kernel_store'
#endif
  
  end subroutine restart_kernel_store
  
  
!============================================================================!
!============================================================================!
!============================================================================!
  
  
  subroutine restart_kernel_retrieve(denskern,label)
  
    !======================================================================!
    ! This subroutine store the density kernel for chained calculation     !
    ! This operation is not affected by changes in atomic positions        !
    !    (input)  density kernel stored in unsegmented format              !
    !    (output) density kernel within denskern in SPAM3 format           !
    ! If (label == 0 or not present)                                       !
    !   retrieve most recent kernel in store                               !
    ! elseif (label=-i < 0)                                                !
    !   retrieve (i+1)th kernel in store                                   !
    ! elseif (label=i > 0)                                                 !
    !   retrieve kernel labelled i                                         !
    ! endif                                                                !
    !----------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois 25/2/2010, modified (02/12/2010).      !
    !======================================================================!
    
    use comms, only : pub_on_root
    use constants, only: DP, stdout
    use timer, only: timer_clock
    use simulation_cell, only : pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_axpy, &
      sparse_destroy, sparse_convert_segment_real
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_abort
    
    implicit none
    
    ! Argument
    type(SPAM3), intent(inout)     :: denskern(:)
    integer, intent(inout), optional  :: label
    
    ! Internal variables
    integer   :: storage_label, storage_idx
    
    !======================================================================!
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering restart_kernel_retrieve'
#endif
    
    call timer_clock("restart_kernel_retrieve", 1)
    
    ! Find storage index
    if (.not. present(label)) then
       storage_label = 0
    else
       storage_label = label
    endif 
    call restart_retrieve_index(storage_idx, storage_label, 'denskern')
    
    ! Return label
    if (present(label)) label = store(storage_idx)%label
    if (pub_on_root) write(stdout,'(a,i5,a)') &
       'restart_kernel_retrieve : retrieve kernel from store  (',storage_idx, ')'
    
    ! Retrieve density kernel from store
    call sparse_convert_segment_real(denskern, &
       store(storage_idx)%dk_idx, store(storage_idx)%dk)
    
    call timer_clock("restart_kernel_retrieve", 2)
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving restart_kernel_retrieve'
#endif
  
  end subroutine restart_kernel_retrieve
  
  
!============================================================================!
!============================================================================!
!============================================================================!
  
  
  subroutine restart_kernel_compose(denskern, rep, ngwf_basis, proj_basis, &
       elements)
  
    !======================================================================!
    ! This subroutine store the density kernel for chained calculation     !
    ! This operation is not affected by changes in atomic positions        !
    !                                                                      !
    !----------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois 25/2/2010, modified (02/12/2010).      !
    !======================================================================!
    
    use comms, only : pub_on_root
    use constants, only: DP, stdout
    use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
      dense_invert
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_brappd_ketppd
    use ion, only: ELEMENT
    use kernel, only: kernel_basis_transform, kernel_basis_update, &
      kernel_christoffel
    use ngwf_representation, only: NGWF_REP
    use projectors, only: projectors_func_ovlp_box
    use pseudopotentials, only: nlps_projectors
    use rundat, only: pub_aug, pub_realspace_projectors, &
      mix_dkn_type, mix_ngwfs_type, mix_dkn_num
    use services, only: services_rms_fit
    use simulation_cell, only : pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_axpy, &
      sparse_destroy, sparse_convert_segment_real, &
      sparse_scale, sparse_copy
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    
    implicit none
    
    ! Argument
    type(SPAM3), intent(inout)     :: denskern(:)
    type(FUNC_BASIS), intent(in)   :: ngwf_basis
    type(FUNC_BASIS), intent(in)   :: proj_basis
    type(NGWF_REP), intent(in)     :: rep
    type(ELEMENT), intent(in)      :: elements(pub_cell%nat)
    
    ! Local Variables
    real(kind=DP), allocatable  :: linear_coeffs(:)
    real(kind=DP), allocatable  :: old_coords(:,:,:)
    real(kind=DP), allocatable  :: old_ngwfs_on_grid(:)
    type(SPAM3), allocatable    :: old_denskern(:)
    type(SPAM3) :: old_overlap
    type(SPAM3) :: old_sp_overlap
    type(SPAM3) :: old_inv_overlap
    type(DEM)   :: old_inv_overlap_dens
    
    ! Internal variables
    integer       :: mix_num, storage_idx, label
    integer       :: istep, is, ierr
    
    !======================================================================!
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering restart_kernel_compose'
#endif
    
    call timer_clock("restart_kernel_compose", 1)
    
    !=============================================================!
    ! No mixing : new denskern = last denskern stored
    !=============================================================!
    if (mix_dkn_type == 0) then
    
       call restart_kernel_retrieve(denskern)
    
    !=============================================================!
    ! Linear mixing :  
    ! algorithm based on Arias et al, PRB 45, 1538 (1992)
    !=============================================================!
    elseif (mix_dkn_type == 1) then
    
       ! Number of density kernel to be mixed
       mix_num = max(1,min(mix_dkn_num, store_nitem-1))
       
       ! Mixing coefficients
       allocate(linear_coeffs(mix_num),stat=ierr)
       call utils_alloc_check('restart_kernel_compose','linear_coeffs',ierr)
       allocate(old_coords(3,pub_cell%nat,mix_num),stat=ierr)
       call utils_alloc_check('restart_kernel_compose','old_coords',ierr)
       
       if (mix_num == 1) then

          ! Interpolation coefficients set to 1.0_dp
          linear_coeffs(:) = 1.0_dp

       else

          ! Set of coordinates to be interpolated
          do istep = 1, mix_num
             label = istep - 1
             call restart_retrieve_index(storage_idx,label,'denskern')
             old_coords(:,:,istep) = store(storage_idx)%coord(:,:)
          enddo

          ! Interpolation coefficients
          call services_rms_fit(linear_coeffs,mix_num,3,pub_cell%nat,&
             old_coords,store(store_pointer(1))%coord)

       endif
       
       ! Linear combination of density kernels
       allocate(old_denskern(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('restart_kernel_compose','old_denskern',ierr)
       do is=1,pub_cell%num_spins
          call sparse_create(old_denskern(is),denskern(is))
       end do
       
       do istep = 1, mix_num
          label = istep - 1
          call restart_kernel_retrieve(old_denskern,label=label)
          do is=1,pub_cell%num_spins
             call sparse_axpy(denskern(is),old_denskern(is),linear_coeffs(istep))
          enddo 
       enddo 
       
       do is=1,pub_cell%num_spins
          call sparse_destroy(old_denskern(is))
       end do
       
       deallocate(old_denskern,stat=ierr)
       call utils_dealloc_check('restart_kernel_compose','old_denskern',ierr)
       deallocate(linear_coeffs,stat=ierr)
       call utils_dealloc_check('restart_kernel_compose','linear_coeffs',ierr)
    
    !=============================================================!
    ! Apply christoffel corrections to last denskern stored
    !=============================================================!
    elseif (retrieve_tightbox_ngwfs .and. mix_dkn_type == 2) then
    
       ! Retrieve density kernel from store
       allocate(old_denskern(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('restart_kernel_compose','old_denskern',ierr)
       do is=1,pub_cell%num_spins
         call sparse_create(old_denskern(is),denskern(is))
       end do
       call restart_kernel_retrieve(old_denskern)
       
       ! Retrieve old ngwfs from store
       allocate(old_ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
       call utils_alloc_check('restart_kernel_compose', &
              'old_ngwfs_on_grid',ierr)
       call restart_ngwfs_tightbox_retrieve(old_ngwfs_on_grid, ngwf_basis, &
              elements)
       
       ! Create sparse arrays for christoffel corrections to density kernel
       call sparse_create(old_sp_overlap,rep%sp_overlap)
       call sparse_create(old_overlap,rep%overlap)
       call sparse_create(old_inv_overlap,rep%inv_overlap)
       call dense_create(old_inv_overlap_dens,ngwf_basis%num,ngwf_basis%num)
       
       ! If required, calculate overlap of old NGWFs with projectors
       if (pub_aug) then
          if (.not.pub_realspace_projectors) then
             call projectors_func_ovlp_box(old_sp_overlap, &
                old_ngwfs_on_grid,ngwf_basis,&
                proj_basis,nlps_projectors)
          else
             call integrals_brappd_ketppd(old_sp_overlap,&
                 old_ngwfs_on_grid, &
                 ngwf_basis,nlps_projectors%projs_on_grid,proj_basis)
          end if
       endif
       
       ! Calculate overlap and inv_overlap of old NGWFs
       call integrals_brappd_ketppd(old_overlap, old_ngwfs_on_grid, &
           ngwf_basis, old_ngwfs_on_grid, ngwf_basis)
       call dense_convert(old_inv_overlap_dens, old_overlap)
       call dense_invert(old_inv_overlap_dens)
       call dense_convert(old_inv_overlap, old_inv_overlap_dens)
       
       ! Apply christoffel corrections to density kernel
       call kernel_christoffel(ngwf_basis, rep%ngwfs_on_grid, &
           old_ngwfs_on_grid, denskern, old_denskern, rep%inv_overlap, &
           old_inv_overlap, rep%sp_overlap, old_sp_overlap)
       
       ! Deallocate temp arrays
       call dense_destroy(old_inv_overlap_dens)
       call sparse_destroy(old_inv_overlap)
       call sparse_destroy(old_overlap)
       call sparse_destroy(old_sp_overlap)
       do is=pub_cell%num_spins,1,-1
          call sparse_destroy(old_denskern(is))
       end do
       
       deallocate(old_ngwfs_on_grid,stat=ierr)
       call utils_dealloc_check('restart_kernel_compose', &
              'start_ngwfs_on_grid',ierr)
       deallocate(old_denskern,stat=ierr)
       call utils_dealloc_check('restart_kernel_compose','start_denskern',ierr)
    
    !=============================================================!
    ! Express last denskern stored in terms of new NGWFs 
    !=============================================================!
    elseif (retrieve_tightbox_ngwfs .and. mix_dkn_type == 3) then
    
       ! Retrieve density kernel from store
       allocate(old_denskern(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('restart_kernel_compose','old_denskern',ierr)
       do is=1,pub_cell%num_spins
         call sparse_create(old_denskern(is),denskern(is))
       end do
       call restart_kernel_retrieve(old_denskern)
       
       ! Retrieve old ngwfs from store
       allocate(old_ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
       call utils_alloc_check('restart_kernel_compose', &
              'old_ngwfs_on_grid',ierr)
       call restart_ngwfs_tightbox_retrieve(old_ngwfs_on_grid, ngwf_basis, &
              elements)
       
       ! Create sparse arrays for christoffel corrections to density kernel
       call sparse_create(old_sp_overlap,rep%sp_overlap)
       
       ! If required, calculate overlap of old NGWFs with projectors
       if (pub_aug) then
         if (.not.pub_realspace_projectors) then
            call projectors_func_ovlp_box(old_sp_overlap, &
               old_ngwfs_on_grid,ngwf_basis,&
               proj_basis,nlps_projectors)
         else
            call integrals_brappd_ketppd(old_sp_overlap,&
               old_ngwfs_on_grid, &
               ngwf_basis,nlps_projectors%projs_on_grid,proj_basis)
         end if
       endif
       
       ! Apply basis transformation to density kernel
       call kernel_basis_transform(ngwf_basis, rep%ngwfs_on_grid, &
          old_ngwfs_on_grid, denskern, old_denskern, rep%inv_overlap, &
          rep%sp_overlap, old_sp_overlap)
       
       ! Deallocate temp arrays
       call sparse_destroy(old_sp_overlap)
       do is=pub_cell%num_spins,1,-1
          call sparse_destroy(old_denskern(is))
       end do
       
       deallocate(old_ngwfs_on_grid,stat=ierr)
       call utils_dealloc_check('restart_kernel_compose', &
              'start_ngwfs_on_grid',ierr)
       deallocate(old_denskern,stat=ierr)
       call utils_dealloc_check('restart_kernel_compose','start_denskern',ierr)
    
    !=============================================================!
    ! Apply corrections to last denskern stored in order to
    ! recover (K.S)_new = (K.S)_old
    !=============================================================!
    elseif (retrieve_tightbox_ngwfs .and. mix_dkn_type == 4) then
    
       ! Retrieve density kernel from store
       allocate(old_denskern(pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('restart_kernel_compose','old_denskern',ierr)
       do is=1,pub_cell%num_spins
          call sparse_create(old_denskern(is),denskern(is))
       end do
       call restart_kernel_retrieve(old_denskern)
       
       ! Retrieve old ngwfs from store
       allocate(old_ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
       call utils_alloc_check('restart_kernel_compose', &
              'old_ngwfs_on_grid',ierr)
       call restart_ngwfs_tightbox_retrieve(old_ngwfs_on_grid, ngwf_basis, &
              elements)
       
       ! Calculate overlap and inv_overlap of old NGWFs
       call sparse_create(old_overlap,rep%overlap)
       call integrals_brappd_ketppd(old_overlap, old_ngwfs_on_grid, &
           ngwf_basis, old_ngwfs_on_grid, ngwf_basis)
       
       ! Apply basis KS corrections to density kernel
       call kernel_basis_update(denskern, old_denskern, &
         rep%overlap, old_overlap, rep%inv_overlap)
       
       ! Deallocate temp arrays
       call sparse_destroy(old_overlap)
       do is=pub_cell%num_spins,1,-1
          call sparse_destroy(old_denskern(is))
       end do
       
       deallocate(old_denskern,stat=ierr)
       call utils_dealloc_check('restart_kernel_compose','old_denskern',ierr)
    
    endif
    
    call timer_clock("restart_kernel_compose", 2)
    
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving restart_kernel_compose'
#endif
  
  end subroutine restart_kernel_compose
  
  
!============================================================================!
!============================================================================!
!============================================================================!
  
  
  subroutine restart_ngwfs_tightbox_store(ngwfs_on_grid, ngwf_basis, &
            elements, label)
  
  
    !========================================================================!
    ! This subroutine store the current NGWFs in a "universal tightbox"      !
    ! representation in list_ngwfs. The order of the NGWFs in list_ngwfs     !
    ! is the one in which they appear in the input file and is not affected  !
    ! by the use of the space-filling curve.                                 !
    !                                                                        !
    ! ndmh 18/12/2009: In cases where the FFTbox coincides with the cell     !
    ! but the individual tightboxes do not need to, a new temporary NGWF     !
    ! basis is constructed, for which the FFTbox only coincides with the     !
    ! cell along direction i if L_i < 2*R_NGWF + \delta_i.                   !
    !                                                                        !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 18/2/2010 based on the                !
    ! restart_ngwfs_tightbox_write subroutine originally written             !
    ! by Chris-Kriton Skylaris on 12/3/2004.                                 !
    !========================================================================!
    
    use comms, only: comms_abort, comms_reduce, comms_barrier, comms_free, &
     pub_my_node_id, pub_on_root, pub_root_node_id
    use constants, only: DP, stdout, NORMAL
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
     function_basis_deallocate, function_basis_distribute, &
     function_basis_copy_spheres, function_basis_init_tight_boxes, &
     function_basis_ppds_to_tightbox
    use geometry, only: magnitude
    use ion, only : ELEMENT
    use parallel_strategy, only: pub_distr_atom, pub_node_of_atom, &
     pub_elements_on_node
    use rundat, only : pub_rootname, pub_output_detail
    use simulation_cell, only : pub_cell, pub_fftbox
    use timer, only: timer_clock
    use utils, only : utils_alloc_check, utils_dealloc_check, &
     utils_unit, utils_abort
    
    implicit none
    
    ! Arguments
    type(FUNC_BASIS), intent(in)  :: ngwf_basis
    real(kind =DP), intent(in)    :: ngwfs_on_grid(ngwf_basis%n_ppds *pub_cell%n_pts) ! NGWFs on this proc
    type(ELEMENT), intent(in)     :: elements(pub_cell%nat) ! elements of all proc (in input file order)
    integer, optional, intent(inout)      :: label
    
    ! Local variables
    real(kind=DP), allocatable, dimension(:,:,:) :: uni_tbox  ! universal tightbox
    real(kind=DP) :: tb_orig1   ! a1-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: tb_orig2   ! a2-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: tb_orig3   ! a3-position of NGWF wrt tightbox in (real number of) grid points
    integer :: maxt_n1     ! a1-maximum tightbox points
    integer :: maxt_n2     ! a2-maximum tightbox points
    integer :: maxt_n3     ! a3-maximum tightbox points
    integer :: orig_atom   ! global atom counter in input file order
    integer :: atom_ngwf   ! ngwf of current atom counter
    integer :: distr_ngwf  ! global ngwf counter in spacefil order
    real(kind=DP) :: sendbuf(3) ! comms send pack buffer
    real(kind=DP) :: recvbuf(3) ! comms recv pack buffer
    real(kind=DP) :: max_radius ! Maximum NGWF radius
    integer :: glob_ngwf, loc_ngwf  ! global/local NGWFs counter
    integer :: storage_idx, storage_label, storage_node   
    type(FUNC_BASIS) :: tmp_basis ! ndmh: for dealing with FFTbox coinciding with cell
    logical :: orig_coin1, orig_coin2, orig_coin3 ! ndmh: whether the fftbox
    		   ! originally coincided with cell in 1,2,3 directions
    integer :: ierr
    
    call timer_clock("restart_ngwfs_tightbox_store", 1)
    
    ! cks: find maximum number of points for universal tightbox that has odd
    ! cks: number of points in each dimension
    maxt_n1 = ngwf_basis%maxtight_pts1
    maxt_n2 = ngwf_basis%maxtight_pts2
    maxt_n3 = ngwf_basis%maxtight_pts3
    
    ! ndmh: if the FFTbox coincides with the cell in any direction, but the
    ! ndmh: tightbox does not necessarily need to, we can define a
    ! ndmh: temporary NGWF basis where the tightboxes are as small as they
    ! ndmh: can be, so that we do not get simulation-cell-sized tightboxes
    ! ndmh: written out to files unnecessarily
    orig_coin1 = pub_fftbox%coin1
    orig_coin2 = pub_fftbox%coin2
    orig_coin3 = pub_fftbox%coin3
    if (orig_coin1.or.orig_coin2.or.orig_coin3) then
    
    max_radius = maxval(ngwf_basis%spheres(:)%radius)
    call comms_reduce('MAX',max_radius)
    
    ! ndmh: temporarily override coin1,2,3 variables if possible
    if (magnitude(pub_fftbox%a1) > 2.0_DP*max_radius + pub_fftbox%d1) &
        pub_fftbox%coin1 = .false.
    if (magnitude(pub_fftbox%a2) > 2.0_DP*max_radius + pub_fftbox%d2) &
        pub_fftbox%coin2 = .false.
    if (magnitude(pub_fftbox%a3) > 2.0_DP*max_radius + pub_fftbox%d3) &
        pub_fftbox%coin3 = .false.
    
    end if
    
    ! ndmh: if any of the coin1,2,3 variables have been changed...
    if ((orig_coin1.neqv.pub_fftbox%coin1).or. &
        (orig_coin2.neqv.pub_fftbox%coin2).or. &
        (orig_coin3.neqv.pub_fftbox%coin3)) then
    
       ! ndmh: create a temporary function basis where tightboxes
       ! ndmh: do not coincide with cell if possible
       call function_basis_allocate(tmp_basis,ngwf_basis%num, &
           'tmp_'//ngwf_basis%name)
       call function_basis_distribute(tmp_basis,elements)
       call function_basis_copy_spheres(tmp_basis,ngwf_basis)
       call function_basis_init_tight_boxes(tmp_basis)
       
       ! ndmh: set a universal tightbox size which depends on the actual size
       ! ndmh: of the temporary NGWF tightboxes
       maxt_n1 = maxval(tmp_basis%tight_boxes(:)%tight_pts1)
       maxt_n2 = maxval(tmp_basis%tight_boxes(:)%tight_pts2)
       maxt_n3 = maxval(tmp_basis%tight_boxes(:)%tight_pts3)
       call comms_reduce('MAX',maxt_n1)
       call comms_reduce('MAX',maxt_n2)
       call comms_reduce('MAX',maxt_n3)
    
    end if
    
    ! cks: allocate universal tightbox buffer
    allocate(uni_tbox(maxt_n1,maxt_n2,maxt_n3),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_store','uni_tbox',ierr)
    
    ! Find storage index
    if (.not. present(label)) then
       storage_label = 0
    else
       storage_label = label
    endif 
    call restart_storage_index(storage_idx,storage_label)
    
    ! Return label
    if (present(label)) label = store(storage_idx)%label
    if (pub_on_root) write(stdout,'(a,i5,a)') &
          'restart_ngwfs_store : storage_idx  (', storage_idx,')'
    
    ! Store contains at least one density kernel
    retrieve_tightbox_ngwfs = .true.
    
    ! smmd: If required free space for storage
    if (associated(store(storage_idx)%ngwfs_nodes)) then
       deallocate(store(storage_idx)%ngwfs,stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_store',&
         'store%ngwfs',ierr)
       nullify(store(storage_idx)%ngwfs)
       deallocate(store(storage_idx)%ngwfs_orig,stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_store',&
         'store%ngwfs_orig',ierr)
       nullify(store(storage_idx)%ngwfs_orig)
       deallocate(store(storage_idx)%ngwfs_nodes,stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_store',&
         'store%ngwfs_nodes',ierr)
       nullify(store(storage_idx)%ngwfs_nodes)
    endif
    
    ! smmd: allocate the storage arrays
    allocate(store(storage_idx)%ngwfs(maxt_n1, &
    maxt_n2,maxt_n3,ngwf_basis%node_num),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_store', &
    'store%ngwfs',ierr)
    allocate(store(storage_idx)%ngwfs_orig(3, &
    ngwf_basis%node_num),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_store', &
    'store%ngwfs_orig',ierr)
    allocate(store(storage_idx)%ngwfs_nodes( &
    ngwf_basis%num),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_nodes', &
    'store%ngwfs_nodes',ierr)
    
    ! smmd: Store tightbox_size
    store(storage_idx)%tightbox_size(:) = &
          (/ maxt_n1, maxt_n2, maxt_n3 /)
    
    ! smmd: Initialise the NGWFs counting indexes
    glob_ngwf = 0
    loc_ngwf = 0
    
    ! cks: Loop over all atoms in the order they appear in the input file
    atom_loop: do orig_atom =1, pub_cell%nat
    
       ngwfs_on_atom_loop: do atom_ngwf= &
           1,ngwf_basis%num_on_atom(pub_distr_atom(orig_atom))
       
         ! smmd: global ngwf index (before distribution)
         glob_ngwf = glob_ngwf + 1
       
         ! smmd: global ngwf index (after distribution)
         distr_ngwf = ngwf_basis%first_on_atom(pub_distr_atom(orig_atom)) + &
              atom_ngwf - 1
       
         ! smmd: storage node
         storage_node = ngwf_basis%node_of_func(distr_ngwf)
       
         ! cks: this essentially converts the non-blocking sends of
         ! cks: comms_send to blocking sends. Since parallel performance
         ! cks: is not the issue here, this makes the code simpler.
         call comms_barrier
       
         if ((orig_coin1.neqv.pub_fftbox%coin1).or. &
           (orig_coin2.neqv.pub_fftbox%coin2).or. &
           (orig_coin3.neqv.pub_fftbox%coin3)) then
       
            call function_basis_ppds_to_tightbox(uni_tbox, &
               tb_orig1,tb_orig2,tb_orig3,distr_ngwf,storage_node, &
               maxt_n1,maxt_n2,maxt_n3,ngwfs_on_grid,tmp_basis, &
               sendbuf,recvbuf)
         else
            call function_basis_ppds_to_tightbox(uni_tbox, &
               tb_orig1,tb_orig2,tb_orig3,distr_ngwf,storage_node, &
               maxt_n1,maxt_n2,maxt_n3,ngwfs_on_grid,ngwf_basis, &
               sendbuf,recvbuf)
         end if
       
         ! Store current NGWF on storage_node
         if (pub_my_node_id==storage_node) then
       
            loc_ngwf = loc_ngwf + 1
            store(storage_idx)%ngwfs(:,:,:,loc_ngwf) = &
               uni_tbox(1:maxt_n1,1:maxt_n2,1:maxt_n3)
            store(storage_idx)%ngwfs_orig(:,loc_ngwf) = &
              (/ tb_orig1, tb_orig2, tb_orig3 /)
         end if
       
         ! Keep the identity of storage_node in memory
         store(storage_idx)%ngwfs_nodes(glob_ngwf) = storage_node
       
       enddo ngwfs_on_atom_loop
    enddo atom_loop
    
    
    ! cks: don't go away just yet! Wait for all communications to finish!
    ! cks: also wait before deallocating uni_tbox which is used in non-blocking
    ! cks: send!
    call comms_barrier
    ! cks: clean up send stack too
    call comms_free
    
    
    ! cks: deallocate universal tightbox buffer
    deallocate(uni_tbox,stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_store','uni_tbox',ierr)
    
    
    ! ndmh: Clean up temporary basis if we created one, and reset coin1,2,3
    if ((orig_coin1.neqv.pub_fftbox%coin1).or. &
        (orig_coin2.neqv.pub_fftbox%coin2).or. &
        (orig_coin3.neqv.pub_fftbox%coin3)) then
    
       call function_basis_deallocate(tmp_basis)
       pub_fftbox%coin1 = orig_coin1
       pub_fftbox%coin2 = orig_coin2
       pub_fftbox%coin3 = orig_coin3
    
    end if
    
    call timer_clock("restart_ngwfs_tightbox_store", 2)
    
    return
  
  end subroutine restart_ngwfs_tightbox_store
  
  
!============================================================================!
!============================================================================!
!============================================================================!
  
  
  subroutine restart_ngwfs_tightbox_compose(ngwfs_on_grid, ngwf_basis, &
  elements)
  
    !========================================================================!
    ! This subroutine creates new NGWFs from data kept in store. New ngwfs   !
    ! are in "universal tightbox" representation from previous stored.       !
    !                                                                        !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 20/06/2011                            !
    !========================================================================!
    
    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use ion, only : ELEMENT
    use rundat, only: pub_rootname, mix_ngwfs_num, &
        mix_ngwfs_type, mix_ngwfs_coeff
    use services, only : services_polynomial_step, services_rms_fit
    use simulation_cell, only : pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit
    
    implicit none
    
    ! Arguments
    type(ELEMENT), intent(in)     :: elements(pub_cell%nat)
    type(FUNC_BASIS), intent(in)  :: ngwf_basis
    real(kind=DP), intent(out)    :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    
    ! Local variables
    real(kind=DP), allocatable  :: buffer_on_grid(:)
    real(kind=DP), allocatable  :: linear_coeffs(:), generalized_coeffs(:,:)
    real(kind=DP), allocatable  :: old_ngwfs_on_grid(:,:)
    real(kind=DP), allocatable  :: old_coords(:,:,:), coord(:,:)
    real(kind=DP), allocatable  :: old_steps(:)
    integer  :: orig_atom
    integer  :: mix_num, label, storage_idx
    integer  :: istep, ierr
    real(kind=DP) :: step
    
    ! Preliminaries 
    if (mix_ngwfs_type .gt. 0) then
       mix_num = max(1,min(mix_ngwfs_num, store_nitem-1))
    endif
    
    !=============================================================!
    ! No mixing : new NGWFs = last NGWFs stored
    !=============================================================!
    if (mix_ngwfs_type == 0 .or. mix_num == 1) then
    
       call restart_ngwfs_tightbox_retrieve(ngwfs_on_grid, ngwf_basis, elements)
    
    !=============================================================!
    ! Linear mixing :  
    ! algorithm based on Arias et al, PRB 45, 1538 (1992)
    !=============================================================!
    elseif (mix_ngwfs_type == 1 .and. mix_num .gt. 1) then
    
       ! Mixing coefficients
       allocate(linear_coeffs(mix_num),stat=ierr)
       call utils_alloc_check('restart_ngwfs_tightbox_compose','linear_coeffs',ierr)
       allocate(old_coords(3,pub_cell%nat,mix_num),stat=ierr)
       call utils_alloc_check('restart_ngwfs_tightbox_compose','old_coords',ierr)
       
       ! Set of coordinates to be interpolated
       do istep = 1, mix_num
         label = istep - 1
         call restart_retrieve_index(storage_idx,label,'ngwfs')
         old_coords(:,:,istep) = store(storage_idx)%coord(:,:)
       enddo
       ! Interpolation coefficients
       call services_rms_fit(linear_coeffs,mix_num,3,pub_cell%nat,&
              old_coords,store(store_pointer(1))%coord)
       
       deallocate(old_coords, stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_compose','old_coords',ierr)
       
       ! Linear combination of NGWFs 
       allocate(buffer_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts), stat=ierr)
       call utils_alloc_check('restart_ngwfs_tightbox_compose','data_on_grid',ierr)
       
       ngwfs_on_grid(:) = 0.0_dp
       ngwfs_loop1: do istep = 1, mix_num
    
          buffer_on_grid(:) = 0.0_dp
          label = istep - 1
          call restart_ngwfs_tightbox_retrieve(buffer_on_grid(:), ngwf_basis, &
                           elements, label=label)
          ngwfs_on_grid(:) = ngwfs_on_grid(:) + linear_coeffs(istep)*buffer_on_grid(:)
       
       enddo ngwfs_loop1
    
       deallocate(linear_coeffs, stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_compose','linear_coeffs',ierr)
       deallocate(buffer_on_grid, stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_compose','buffer_on_grid',ierr)
    
    !=============================================================!
    ! Generalized linear mixing :  
    ! for each atom, a different set of coefficients is computed
    ! as a function of its local environemet 
    !=============================================================!
    elseif (mix_ngwfs_type == 2 .and. mix_num.gt.1) then
    
       ! Mixing coefficients
       allocate(generalized_coeffs(mix_num,pub_cell%nat),stat=ierr)
       call utils_alloc_check('restart_ngwfs_tightbox_compose','generalized_coeffs',ierr)
       allocate(linear_coeffs(mix_num),stat=ierr)
       call utils_alloc_check('restart_ngwfs_tightbox_compose','linear_coeffs',ierr)
       
       ! Set of coordinates
       allocate(old_coords(3,pub_cell%nat,mix_num),stat=ierr)
       call utils_alloc_check('restart_ngwfs_tightbox_compose','old_coords',ierr)
       allocate(coord(3,pub_cell%nat),stat=ierr)
       call utils_alloc_check('restart_ngwfs_tightbox_compose','coord',ierr)
       do istep = 1, mix_num
          label = istep - 1
          call restart_retrieve_index(storage_idx,label,'ngwfs')
          old_coords(:,:,istep) = store(storage_idx)%coord(:,:)
       enddo
       coord(:,:) = store(store_pointer(1))%coord(:,:)
       
       ! Loop over all atoms 
       atom_loop: do orig_atom = 1, pub_cell%nat

          ! Compute weigthed atomic coordinates
          call generalized_coordinates(mix_num,old_coords,coord,coord(:,orig_atom))

          ! Interpolation coefficients
          call services_rms_fit(generalized_coeffs(:,orig_atom),mix_num, &
                3,pub_cell%nat,old_coords,coord) 

       enddo atom_loop

       deallocate(old_coords, stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_compose','old_coords',ierr)
       deallocate(coord, stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_compose','coord',ierr)
       
       ! Linear combination of NGWFs 
       allocate(buffer_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts), stat=ierr)
       call utils_alloc_check('restart_ngwfs_tightbox_compose','data_on_grid',ierr)
       ngwfs_on_grid(:) = 0.0_dp
       
       ngwfs_loop2: do istep = 1, mix_num
       
         label = istep-1
         buffer_on_grid(:) = 0.0_dp
         linear_coeffs(:) = generalized_coeffs(istep,:)          
         call restart_ngwfs_tightbox_retrieve(buffer_on_grid(:), ngwf_basis, &
               elements, label=label, coeffs=linear_coeffs)
         ngwfs_on_grid(:) = ngwfs_on_grid(:) + buffer_on_grid(:)
       
       enddo ngwfs_loop2
       
       deallocate(buffer_on_grid, stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_compose','buffer_on_grid',ierr)
       deallocate(linear_coeffs, stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_compose','linear_coeffs',ierr)
       deallocate(generalized_coeffs, stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_compose','generalized_coeffs',ierr)
  
  
    !=============================================================!
    ! Polynomial extrapolation in real space
    !=============================================================!
    elseif (mix_ngwfs_type == 3 .and. mix_num.gt.1) then
  
       ngwfs_on_grid(:) = 0.0_dp
  
       ! Loop over the ngwfs to be mixed
       allocate(old_steps(mix_num), stat=ierr)
       call utils_alloc_check('restart_ngwfs_tightbox_compose','time_ref',ierr)
       allocate(buffer_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts), stat=ierr)
       call utils_alloc_check('restart_ngwfs_tightbox_compose','buffer_on_grid',ierr)
       allocate(old_ngwfs_on_grid(mix_num,ngwf_basis%n_ppds*pub_cell%n_pts), stat=ierr)
       call utils_alloc_check('restart_ngwfs_tightbox_compose','buffer_poly',ierr)
       buffer_on_grid(:) = 0.0_dp
       old_ngwfs_on_grid(:,:) = 0.0_dp
  
       ngwfs_loop3: do istep = 1, mix_num
       
          label = istep - 1 
          buffer_on_grid(:) = 0.0_dp
          call restart_ngwfs_tightbox_retrieve(buffer_on_grid, ngwf_basis, &
                   elements, label=label)
  
          old_ngwfs_on_grid(istep,:) = buffer_on_grid(:)
          old_steps(istep) = real(istep)
       
       enddo ngwfs_loop3
  
       step = real(mix_num+1)
       call services_polynomial_step(ngwfs_on_grid,mix_num,ngwf_basis%n_ppds*pub_cell%n_pts,&
              old_ngwfs_on_grid,old_steps,step)
    
       ngwfs_on_grid(:) = (1.0_dp-mix_ngwfs_coeff)*old_ngwfs_on_grid(1,:) + &
                     mix_ngwfs_coeff*ngwfs_on_grid(:)
       
       deallocate(buffer_on_grid, stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_compose','buffer_poly',ierr)
       deallocate(old_ngwfs_on_grid, stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_compose','buffer_poly',ierr)
       deallocate(old_steps, stat=ierr)
       call utils_dealloc_check('restart_ngwfs_tightbox_compose','time_ref',ierr)
  
    endif 
  
  end subroutine restart_ngwfs_tightbox_compose
  
!============================================================================!
!============================================================================!
!============================================================================!

  subroutine generalized_coordinates(ncoords,old_coords,coord,ref_point)

    use constants, only: DP, stdout
    use fourier, only: fourier_apply_box
    use simulation_cell, only : pub_cell
    use rundat, only: mix_local_smear, mix_local_length
    use services, only: services_rationalise_coords
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none
    
    ! Arguments
    integer, intent(in)  :: ncoords 
    real(kind=DP), intent(inout) :: coord(3,pub_cell%nat) 
    real(kind=DP), intent(inout) :: old_coords(3,pub_cell%nat,ncoords) 
    real(kind=DP), intent(in) :: ref_point(3) 

    ! Local variables
    real(kind=DP) :: center_shift(3), dist
    real(kind=DP) :: atom_weight(pub_cell%nat)
    integer       :: ic, iat

    ! Center point of simulation cell 
    center_shift(1) = (pub_cell%a1%x + pub_cell%a2%x + pub_cell%a3%x)/2.0d0
    center_shift(2) = (pub_cell%a1%y + pub_cell%a2%y + pub_cell%a3%y)/2.0d0
    center_shift(3) = (pub_cell%a1%z + pub_cell%a2%z + pub_cell%a3%z)/2.0d0
    
    ! Center coordinates around ref_point and rationalise 
    do iat = 1, pub_cell%nat
       coord(:,iat) = coord(:,iat) - ref_point(:) &
                                + center_shift(:)
    enddo
    call services_rationalise_coords(pub_cell%nat,coord)
    do ic = 1, ncoords  
       do iat = 1, pub_cell%nat
          old_coords(:,iat,ic) = old_coords(:,iat,ic) - ref_point(:) &
                                   + center_shift(:)
       enddo
       call services_rationalise_coords(pub_cell%nat,old_coords(:,:,ic))
    enddo

    ! Compute atomic weight 
    atom_weight(:) = 0.0d0    
    do iat = 1, pub_cell%nat
       dist = sqrt(sum((coord(:,iat)-ref_point(:))**2))
       if (mix_local_smear.le.0.0d0) then
          if (dist.lt.mix_local_length) atom_weight(iat) = 1.0d0  
       else
          atom_weight(iat) = (exp((dist-mix_local_length)/mix_local_smear)+1)**(-1)
       endif
    enddo

    ! Apply atomic weight to coordinates
    do iat = 1, pub_cell%nat
       if (atom_weight(iat).gt.0.0d0) then
          coord(:,iat) = coord(:,iat)*atom_weight(iat)
       endif 
    enddo
    do ic = 1, ncoords  
       do iat = 1, pub_cell%nat
          if (atom_weight(iat).gt.0.0d0) then
             old_coords(:,iat,ic) = old_coords(:,iat,ic)*atom_weight(iat)
          endif 
       enddo
    enddo

  end subroutine generalized_coordinates

!============================================================================!
!============================================================================!
!============================================================================!


  subroutine restart_ngwfs_tightbox_retrieve(ngwfs_on_grid, ngwf_basis, &
       elements, label, coeffs)

    !========================================================================!
    ! This subroutine retrieve the  NGWFs in "universal tightbox"            !
    ! representation from list_ngwfs. The order of the NGWFs in list_ngwfs   !
    ! is the one in which they appear in the input file and is not affected  !
    ! by the use of the space-filling curve.                                 !
    !                                                                        !
    ! ndmh 18/12/2009: In cases where the FFTbox coincides with the cell     !
    ! but the individual tightboxes do not need to, a new temporary NGWF     !
    ! basis is constructed, for which the FFTbox only coincides with the     !
    ! cell along direction i if L_i < 2*R_NGWF + \delta_i.                   !
    !                                                                        !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois on 18/2/2010 based on the                !
    ! restart_ngwfs_tightbox_read subroutine originally written              !
    ! by Chris-Kriton Skylaris on 12/3/2004.                                 !
    !========================================================================!

    use basis, only: basis_ket_start_wrt_fftbox
    use comms, only: comms_abort, comms_barrier, comms_bcast, comms_free, &
         comms_recv, comms_reduce, comms_send, pub_my_node_id, pub_on_root, &
         pub_root_node_id
    use constants, only: DP, stdout
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS, function_basis_allocate, &
         function_basis_deallocate, function_basis_distribute, &
         function_basis_copy_spheres, function_basis_init_tight_boxes, &
         function_basis_tightbox_to_ppds
    use geometry, only: magnitude
    use ion, only : ELEMENT
    use parallel_strategy, only: pub_distr_atom, pub_node_of_atom, &
         pub_elements_on_node
    use rundat, only: pub_rootname
    use simulation_cell, only : pub_cell, pub_fftbox
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_unit, utils_abort

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(out) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts) ! NGWFs on this proc
    type(ELEMENT), intent(in) :: elements(pub_cell%nat) ! elements of all proc (in input file order)
    real(kind=DP), intent(in), optional :: coeffs(pub_cell%nat) 
    integer, intent(inout), optional  :: label

    ! Local Variables
    real(kind=DP), allocatable, dimension(:,:,:) :: uni_tbox  ! universal tightbox
    complex(kind=DP), allocatable, dimension(:,:,:) :: fftbox_complex ! complex fftbox space
    complex(kind=DP), allocatable, dimension(:,:,:) :: fftbox_complex_shifted ! complex fftbox space
    real(kind=DP), allocatable, dimension(:,:,:) :: fftbox_buffer     ! real fftbox space
    real(kind=DP) :: max_radius          ! globaly maximum ngwf radius
    real(kind=DP) :: read_tb_orig1 ! read a1-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: read_tb_orig2 ! read a2-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: read_tb_orig3 ! read a3-position of NGWF wrt tightbox in (real number of) grid points
    real(kind=DP) :: read_weight   ! read grid point weight
    integer :: n1, n2, n3, ld1, ld2 ! fftbox dimensions
    integer :: orig_atom    ! global atom counter in input file order
    integer :: atom_ngwf    ! ngwf of current atom counter
    integer :: maxt_n1      ! a1-maximum tightbox points
    integer :: maxt_n2      ! a2-maximum tightbox points
    integer :: maxt_n3      ! a3-maximum tightbox points
    integer :: read_maxt_n1 ! a1-maximum tightbox points read from file
    integer :: read_maxt_n2 ! a2-maximum tightbox points read from file
    integer :: read_maxt_n3 ! a3-maximum tightbox points read from file
    integer :: tb_start1    ! a1-grid point start of universal tightbox wrt fftbox
    integer :: tb_start2    ! a2-grid point start of universal tightbox wrt fftbox
    integer :: tb_start3    ! a3-grid point start of universal tightbox wrt fftbox
    integer :: distr_ngwf   ! global ngwf counter in spacefill order
    real(kind=DP) :: sendbuf(3) ! comms send pack buffer
    real(kind=DP) :: recvbuf(3) ! comms recv pack buffer
    integer :: ierr         ! memory allocation error flag

    ! smmd:
    real(kind=DP), allocatable :: buffer_on_grid(:)
    integer :: glob_ngwf, loc_ngwf  ! global/local NGWFs counter
    integer :: storage_label, storage_idx
    integer :: storage_node

    ! ndmh:
    type(FUNC_BASIS) :: tmp_basis
    logical :: orig_coin1,orig_coin2,orig_coin3 ! ndmh: whether the fftbox
                           ! originally coincided with cell in 1,2,3 directions

    call timer_clock("restart_ngwfs_tightbox_retrieve", 1)

    ! Find storage index
    if (.not. present(label)) then
       storage_label = 0
    else
       storage_label = label
    endif 
    call restart_retrieve_index(storage_idx,storage_label,'ngwfs')

    ! Return label
    if (present(label)) label = store(storage_idx)%label
    if (pub_on_root) write(stdout,'(a,i5,a)') &
       'restart_ngwfs_tightbox_retrieve : retrieve ngwfs from store  (',storage_idx, ')'


    ! cks: ********** INITIALISATIONS ***************************************
    ngwfs_on_grid = 0.0_DP

    ! cks: find maximum number of (odd numbers) points for universal tightbox
    max_radius = maxval(ngwf_basis%spheres(:)%radius)
    call comms_reduce('MAX',max_radius)
    maxt_n1 = int(2.0_DP*max_radius / pub_cell%d1) + 1
    maxt_n1 = maxt_n1 + 1 - mod(maxt_n1,2)
    maxt_n2 = int(2.0_DP*max_radius / pub_cell%d2) + 1
    maxt_n2 = maxt_n2 + 1 - mod(maxt_n2,2)
    maxt_n3 = int(2.0_DP*max_radius / pub_cell%d3) + 1
    maxt_n3 = maxt_n3 + 1 - mod(maxt_n3,2)

    n1 = pub_fftbox%total_pt1
    n2 = pub_fftbox%total_pt2
    n3 = pub_fftbox%total_pt3
    ld1 = pub_fftbox%total_ld1
    ld2 = pub_fftbox%total_ld2

    read_maxt_n1 = store(storage_idx)%tightbox_size(1)
    read_maxt_n2 = store(storage_idx)%tightbox_size(2)
    read_maxt_n3 = store(storage_idx)%tightbox_size(3)
    ! cks: ****** END INITIALISATIONS ***************************************


    ! ndmh: if the FFTbox coincides with the cell in any direction, but the
    ! ndmh: tightbox does not necessarily need to, we may need to define a
    ! ndmh: temporary NGWF basis where the tightboxes are as small as they
    ! ndmh: can be, unless the file we are reading in was written with
    ! ndmh: tightboxes coinciding with cell
    orig_coin1 = pub_fftbox%coin1
    orig_coin2 = pub_fftbox%coin2
    orig_coin3 = pub_fftbox%coin3
    if (orig_coin1.or.orig_coin2.or.orig_coin3) then

       ! ndmh: temporarily override coin1,2,3 variables if possible (as long
       ! ndmh: as the written values are not for tightbox==fftbox)
       if ((magnitude(pub_fftbox%a1) > 2.0_DP*max_radius + &
            pub_fftbox%d1*pub_cell%n_pt1).and. &
            (read_maxt_n1 /= pub_fftbox%total_pt1)) pub_fftbox%coin1 = .false.
       if ((magnitude(pub_fftbox%a2) > 2.0_DP*max_radius + &
            pub_fftbox%d2*pub_cell%n_pt1).and. &
            (read_maxt_n2 /= pub_fftbox%total_pt2)) pub_fftbox%coin2 = .false.
       if ((magnitude(pub_fftbox%a3) > 2.0_DP*max_radius + &
            pub_fftbox%d3*pub_cell%n_pt1).and. &
            (read_maxt_n3 /= pub_fftbox%total_pt3)) pub_fftbox%coin3 = .false.

#ifdef DEBUG
       if ((orig_coin1.neqv.pub_fftbox%coin1).and.pub_on_root) write(stdout,*) &
            'DEBUG: Temporarily overriding pub_fftbox%coin1'
       if ((orig_coin2.neqv.pub_fftbox%coin2).and.pub_on_root) write(stdout,*) &
            'DEBUG: Temporarily overriding pub_fftbox%coin2'
       if ((orig_coin3.neqv.pub_fftbox%coin3).and.pub_on_root) write(stdout,*) &
            'DEBUG: Temporarily overriding pub_fftbox%coin3'
#endif
    end if

    ! ndmh: if any of the coin1,2,3 variables have been changed...
    if ((orig_coin1.neqv.pub_fftbox%coin1).or. &
         (orig_coin2.neqv.pub_fftbox%coin2).or. &
         (orig_coin3.neqv.pub_fftbox%coin3)) then

       ! ndmh: create a temporary function basis where tightboxes
       ! ndmh: do not coincide with cell if possible
       call function_basis_allocate(tmp_basis,ngwf_basis%num, &
            'tmp_'//ngwf_basis%name)
       call function_basis_distribute(tmp_basis,elements)
       call function_basis_copy_spheres(tmp_basis, ngwf_basis)
       call function_basis_init_tight_boxes(tmp_basis)

       ! ndmh: set a universal tightbox size which depends on the actual size
       ! ndmh: of the temporary NGWF tightboxes
       maxt_n1 = maxval(tmp_basis%tight_boxes(:)%tight_pts1)
       maxt_n2 = maxval(tmp_basis%tight_boxes(:)%tight_pts2)
       maxt_n3 = maxval(tmp_basis%tight_boxes(:)%tight_pts3)
       call comms_reduce('MAX',maxt_n1)
       call comms_reduce('MAX',maxt_n2)
       call comms_reduce('MAX',maxt_n3)

    end if

    ! cks: NGWFs will be placed in the centre of the fftbox
    ! cks: or left where they are in the simulation cell
    ! cks: (when the simulation cell and FFT-box coincide)
    call basis_ket_start_wrt_fftbox(tb_start1, tb_start2, tb_start3, n1, n2, n3)

    ! ndmh: check tightbox will not extend out of FFTbox if the coin flag has
    ! ndmh: been overridden. If so, try to adjust tb_start to compensate. If
    ! ndmh: that fails, print an error and quit.
    if ((tb_start1 + maxt_n1 > n1).and.(.not.pub_fftbox%coin1)) then
#ifdef DEBUG
       if (pub_on_root) write(stdout,*) 'DEBUG: Attempting to adjust tb_start1'
#endif
       if (n1 - maxt_n1 > 0) then
          tb_start1 = (n1 - maxt_n1)/2
       else
          write(stdout,'(a,i4,a)') 'Error in restart_ngwfs_tightbox_retrieve: &
               &Tightbox size ',maxt_n1, 'along 1-axis'
          write(stdout,'(a,i4)') 'is not compatible with FFT-box size ',n1
          call comms_abort
       end if
    end if
    if ((tb_start2 + maxt_n2 > n2).and.(.not.pub_fftbox%coin2)) then
#ifdef DEBUG
       if (pub_on_root) write(stdout,*) 'DEBUG: Attempting to adjust tb_start2'
#endif
       if (n2 - maxt_n2 > 0) then
          tb_start2 = (n2 - maxt_n2)/2
       else
          write(stdout,'(a,i4,a)') 'Error in restart_ngwfs_tightbox_retrieve: &
               &Tightbox size ',maxt_n2, 'along 2-axis'
          write(stdout,'(a,i4)') 'is not compatible with FFT-box size ',n2
          call comms_abort
       end if
    end if
    if ((tb_start3 + maxt_n3 > n3).and.(.not.pub_fftbox%coin3)) then
#ifdef DEBUG
       if (pub_on_root) write(stdout,*) 'DEBUG: Attempting to adjust tb_start3'
#endif
       if (n3 - maxt_n3 > 0) then
          tb_start3 = (n3 - maxt_n3)/2
       else
          write(stdout,'(a,i4,a)') 'Error in restart_ngwfs_tightbox_retrieve: &
               &Tightbox size ',maxt_n3, 'along 3-axis'
          write(stdout,'(a,i4)') 'is not compatible with FFT-box size ',n3
          call comms_abort
       end if
    end if

    ! cks: allocate universal tightbox buffer
    maxt_n1 =max(min(maxt_n1,n1), read_maxt_n1)
    maxt_n2 =max(min(maxt_n2,n2), read_maxt_n2)
    maxt_n3 =max(min(maxt_n3,n3), read_maxt_n3)

    ! cks: send the new maxt's from the root to all other nodes
    call comms_bcast(pub_root_node_id, maxt_n1)
    call comms_bcast(pub_root_node_id, maxt_n2)
    call comms_bcast(pub_root_node_id, maxt_n3)

    ! cks: send read_weight to all nodes
    call comms_bcast(pub_root_node_id, read_weight)

    ! ndmh: allocate temporary storage
    allocate(uni_tbox(maxt_n1,maxt_n2,maxt_n3),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_retrieve','uni_tbox',ierr)
    allocate(fftbox_complex(ld1,ld2,n3),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_retrieve','fftbox_complex',ierr)
    allocate(fftbox_complex_shifted(ld1,ld2,n3),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_retrieve', &
         'fftbox_complex_shifted',ierr)
    allocate(fftbox_buffer(ld1,ld2,n3),stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_retrieve','fftbox_buffer',ierr)
    allocate(buffer_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts), stat=ierr)
    call utils_alloc_check('restart_ngwfs_tightbox_retrieve','buffer_on_grid',ierr)

    ! Initialisation
    ngwfs_on_grid(:) = 0.d0
    loc_ngwf = 0
    glob_ngwf = 0

    ! cks: Loop over all atoms in the order they appear in the input file
    atom_loop: do orig_atom =1, pub_cell%nat

       ! Initialise buffer
       buffer_on_grid(:) = 0.0d0

       ngwfs_on_atom_loop: do atom_ngwf=1, &
            ngwf_basis%num_on_atom(pub_distr_atom(orig_atom))

          ! smmd: global ngwfs index before distribution
          glob_ngwf = glob_ngwf + 1

          ! smmd: global ngwfs index after distribution
          distr_ngwf = ngwf_basis%first_on_atom(pub_distr_atom(orig_atom)) + &
               atom_ngwf - 1

          ! smmd: storage node
          storage_node = store(storage_idx)%ngwfs_nodes(glob_ngwf)

          ! cks: this essentially converts the non-blocking sends of
          ! cks: comms_send to blocking sends. Since parallel performance
          ! cks: is not the issue here, this makes the code simpler.
          call comms_barrier

          ! cks: initialise
          uni_tbox = 0.0_DP

          ! smmd: if storage node, read data into buffer
          if (pub_my_node_id == storage_node) then

             loc_ngwf = loc_ngwf + 1

             read_tb_orig1 = store(storage_idx)%ngwfs_orig(1,loc_ngwf)
             read_tb_orig2 = store(storage_idx)%ngwfs_orig(2,loc_ngwf)
             read_tb_orig3 = store(storage_idx)%ngwfs_orig(3,loc_ngwf)

             uni_tbox(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3) = &
                 store(storage_idx)%ngwfs(1:read_maxt_n1,1:read_maxt_n2,1:read_maxt_n3,loc_ngwf)

          endif

          if ((orig_coin1.neqv.pub_fftbox%coin1).or. &
            (orig_coin2.neqv.pub_fftbox%coin2).or. &
            (orig_coin3.neqv.pub_fftbox%coin3)) then

             call function_basis_tightbox_to_ppds(uni_tbox, &
                  read_tb_orig1,read_tb_orig2,read_tb_orig3, &
                  distr_ngwf,storage_node,maxt_n1,maxt_n2,maxt_n3, &
                  buffer_on_grid,tmp_basis,tb_start1,tb_start2,tb_start3, &
                  fftbox_complex,fftbox_complex_shifted,fftbox_buffer, &
                  sendbuf,recvbuf)

          else

             call function_basis_tightbox_to_ppds(uni_tbox, &
                  read_tb_orig1,read_tb_orig2,read_tb_orig3, &
                  distr_ngwf,storage_node,maxt_n1,maxt_n2,maxt_n3, &
                  buffer_on_grid,ngwf_basis,tb_start1,tb_start2,tb_start3, &
                  fftbox_complex,fftbox_complex_shifted,fftbox_buffer, &
                  sendbuf,recvbuf)

          end if

       enddo ngwfs_on_atom_loop

       if (present(coeffs)) then
          ngwfs_on_grid(:) = ngwfs_on_grid(:) + &
                    buffer_on_grid(:) * coeffs(orig_atom)
       else
          ngwfs_on_grid(:) = ngwfs_on_grid(:) + buffer_on_grid(:)
       endif

    enddo atom_loop

    ! cks: don't go away just yet! Wait for all communications to finish!
    call comms_barrier
    ! cks: clean up send stack too
    call comms_free


    deallocate(buffer_on_grid, stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_retrieve','buffer_on_grid',ierr)
    deallocate(fftbox_buffer,stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_retrieve', &
         'fftbox_buffer',ierr)
    deallocate(fftbox_complex_shifted,stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_retrieve', &
         'fftbox_complex_shifted',ierr)
    deallocate(fftbox_complex,stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_retrieve', &
         'fftbox_complex',ierr)
    deallocate(uni_tbox,stat=ierr)
    call utils_dealloc_check('restart_ngwfs_tightbox_retrieve', &
         'uni_tbox',ierr)

    ! ndmh: Clean up temporary basis if we created one, and reset coin1,2,3
    if ((orig_coin1.neqv.pub_fftbox%coin1).or. &
         (orig_coin2.neqv.pub_fftbox%coin2).or. &
         (orig_coin3.neqv.pub_fftbox%coin3)) then

       call function_basis_deallocate(tmp_basis)
       pub_fftbox%coin1 = orig_coin1
       pub_fftbox%coin2 = orig_coin2
       pub_fftbox%coin3 = orig_coin3

    end if

    call timer_clock("restart_ngwfs_tightbox_retrieve", 2)

    return

  end subroutine restart_ngwfs_tightbox_retrieve


!============================================================================!
!============================================================================!
!============================================================================!

  subroutine restart_storage_index(idx,label)

    !========================================================================!
    ! Written by Simon M.-M. Dubois 12/09/2011                               !
    !========================================================================!

    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout
    use utils, only: utils_abort

    implicit none

    ! Argument
    integer, intent(in)     :: label    
    integer, intent(out)    :: idx   

    ! Internal variables
    integer                 :: idk

    !======================================================================!

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
       'DEBUG: Entering restart_storage_index'
#endif

    ! Find storage index
    idx = 0
    if (label == 0) then
       idx = store_pointer(1)
    elseif (label .lt. 0) then
       if (abs(label) .ge. store_size) then
          if (pub_on_root) call utils_abort('Error in restart_storage_index: &
                   &storage label incompatible with store size')
       endif 
       idx = store_pointer(1+abs(label))
    elseif (label .gt. 0) then
       do idk = 1, store_size
          if (store(idk)%label == label) then
             idx = idk
             exit
          endif
       enddo
    endif
    if (idx == 0) then
       if (pub_on_root) call utils_abort('Error in restart_storage_index: &
                &storage label incompatible with storage elements ')
    endif 

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
       'DEBUG: Leave restart_storage_index'
#endif

    return

  end subroutine restart_storage_index

!============================================================================!
!============================================================================!
!============================================================================!

  subroutine restart_retrieve_index(idx,label,task)

    !========================================================================!
    ! Written by Simon M.-M. Dubois 12/09/2011                               !
    !========================================================================!

    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout
    use utils, only: utils_abort

    implicit none

    ! Argument
    integer, intent(in)     :: label    
    integer, intent(out)         :: idx   
    character(len=*), intent(in)   :: task   

    ! Internal variables
    integer                 :: idk, iel
    logical                 :: dkn_idx, ngwfs_idx

    !======================================================================!

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
       'DEBUG: Entering restart_retrieve_index'
#endif

    dkn_idx = .false.
    ngwfs_idx = .false.
    if (task == 'denskern') then
       dkn_idx = .true.
    elseif (task == 'ngwfs') then
       ngwfs_idx = .true.
    endif

    ! Find storage index
    idx = 0
    if (label == 0) then
       do idk = 1, store_size
          if ((associated(store(store_pointer(idk))%dk) .and. dkn_idx) .or. &
               (associated(store(store_pointer(idk))%ngwfs) .and. ngwfs_idx)) then
             idx = store_pointer(idk)
             exit
          endif
       enddo
    elseif (label .lt. 0) then
       iel = 0
       do idk = 1, store_size
          if ((associated(store(store_pointer(idk))%dk) .and. dkn_idx) &
             .or.  (associated(store(store_pointer(idk))%ngwfs) .and. ngwfs_idx) &
             .and. iel == abs(label)) then
             idx = store_pointer(idk)
             exit
          elseif ((associated(store(store_pointer(idk))%dk) .and. dkn_idx) &
             .or.  (associated(store(store_pointer(idk))%ngwfs) .and. ngwfs_idx)) then
             iel = iel + 1
          endif
       enddo
    elseif (label .gt. 0) then
       do idk = 1, store_size
          if ((associated(store(store_pointer(idk))%dk) .and. dkn_idx) &
             .or.  (associated(store(store_pointer(idk))%ngwfs) .and. ngwfs_idx) &
             .and. store(idk)%label == label) then
             idx = idk
             exit
          endif
       enddo
    endif

    if (idx == 0) then
       if (pub_on_root) call utils_abort('Error in restart_retrieve_index: &
                &storage label incompatible with storage elements ')
    endif 

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
       'DEBUG: Leave restart_retrieve_index'
#endif

    return

  end subroutine restart_retrieve_index

!============================================================================!
!============================================================================!
!============================================================================!

end module restart
