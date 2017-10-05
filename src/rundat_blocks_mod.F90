! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:

module rundat_blocks

  !=====================================================================!
  ! This module contains the subroutines that read the simulation cell  !
  ! and species parameters from the input file in form of data blocks.  !
  !                                                                     !
  ! The routines in this module are adapted from the former routines    !
  ! services_read_system_input and services_read_system_input_old.      !
  !---------------------------------------------------------------------!
  ! This module was created by Alvaro Ruiz Serrano in November 2010.    !
  !=====================================================================!

  implicit none

  public :: rundat_blocks_exec

contains



  subroutine rundat_blocks_exec(elements,nat,tag)

    !==============================================================!
    ! Driver routine that reads esdf_block from the input file and !
    ! initialises various parameters of the simulation cell and    !
    ! the ppd grid.                                                !
    !--------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano on 17/11/2010.                !
    !==============================================================!

    use ion, only: ELEMENT
    use rundat, only: pub_old_input

    implicit none

    ! ars: arguments
    integer,                    intent(in   ) :: nat
    type(ELEMENT),              intent(inout) :: elements(nat)
    character(len=*), optional, intent(in   ) :: tag

    ! aam: Initialise elements and most of pub_cell
    if (pub_old_input) then
       call rundat_blocks_old(elements,nat)
    else
       call rundat_blocks_new(elements, nat, tag)
    endif

  end subroutine rundat_blocks_exec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rundat_blocks_new(elements,nat,tag)

    !==============================================================!
    ! Driver routine that reads esdf_block from the input file     !
    ! using new input file syntax and initialises various          !
    ! parameters of the simulation cell and the ppd grid.          !
    !--------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano on 17/11/2010.                !
    !==============================================================!

    use geometry, only: POINT
    use ion, only: ELEMENT
    use rundat, only: pub_is_implicit_solvent
    use simulation_cell, only: pub_cell
    use utils, only: utils_assert

    implicit none

    ! ars: arguments
    integer,                    intent(in   ) :: nat
    type(ELEMENT),              intent(inout) :: elements(nat)
    character(len=*), optional, intent(in   ) :: tag

    ! ars: local variables
    type(POINT) :: a1,a2,a3
    type(POINT) :: a1_pad,a2_pad,a3_pad
    integer     :: num, num_cond, num_aux, class_nat

    ! ars: read blocks from input file
    call rundat_blocks_read(a1, a2, a3, a1_pad, a2_pad, a3_pad, &
         num, num_cond, num_aux, class_nat, elements, nat, tag)

    ! ars: initialise grid
    call rundat_blocks_init_grid(a1, a2, a3, a1_pad, a2_pad, a3_pad, &
         elements, num, num_cond, num_aux, nat, class_nat)

    ! ars: check compatibility
    call utils_assert(.not. ((pub_cell%nat_classical .gt. 0) &
         .and. pub_is_implicit_solvent),'is_implicit_solvent cannot be used &
         &simultaneously with classical atoms, sorry.')

  end subroutine rundat_blocks_new


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine rundat_blocks_read(a1, a2, a3, a1_pad, a2_pad, a3_pad, &
       num, num_cond, num_aux, class_nat, elements, nat, tag)

    !=============================================================!
    ! This subroutine reads the input file using the esdf module  !
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2001.                   !
    ! Modified by Peter D. Haynes in 2004 so that the root node   !
    ! reads and communicates across all nodes.                    !
    ! Modified by Chris-Kriton Skylaris on 09/11/2004 so that it  !
    ! determines grids by calling ppd_strategy rather than        !
    ! reading them from input.                                    !
    !=============================================================!
    ! Rewritten by Arash A Mostofi, July 2005                     !
    ! NGWF plot block added by Chris-Kriton Skylaris on 13/03/2006!
    ! Hubbard DFT+U block added by D.D. O'Regan on 15/04/2009     !
    !=============================================================!
    ! Adapted and cleaned up from services_read_system_input by   !
    ! Alvaro Ruiz Serrano, 16/11/2010.                            !
    !=============================================================!


    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use geometry, only: POINT
    use ion, only: ELEMENT
    use utils, only: utils_alloc_check, utils_dealloc_check


    implicit none

    ! ars: arguments
    type(POINT),                intent(  out) :: a1,a2,a3
    type(POINT),                intent(  out) :: a1_pad,a2_pad,a3_pad
    integer,                    intent(  out) :: num,num_cond,num_aux,class_nat
    integer,                    intent(in   ) :: nat
    type(ELEMENT),              intent(inout) :: elements(nat)
    character(len=*), optional, intent(in   ) :: tag


    ! ars: ---- global variables ----
    integer :: nsp, hub_nsp, cdft_nsp, ierr
    integer :: num_hub_proj
    character(len=4), pointer, dimension(:) :: id_list
    character(len=4), allocatable :: all_species_id(:)


#ifdef DEBUG
    if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering rundat_blocks_read'
#endif



    ! ars: ============= READ LATTICE BLOCKS =============
    call internal_lattice_cart()
    call internal_padded_lattice_cart()
    ! ars: =========== END READ LATTICE BLOCKS ===========




    ! ars: ============ READ ATOMIC POSITIONS ============
    call internal_positions_abs()

    ! ars: allocate global arrays
    if (pub_on_root) then
       allocate(id_list(1:nsp),stat=ierr)
       call utils_alloc_check('rundat_blocks_read','id_list',ierr)
    end if
    allocate(all_species_id(1:nsp), stat=ierr)
    call utils_alloc_check('rundat_blocks_read','all_species_id',ierr)
    ! ars: ========== END READ ATOMIC POSITIONS ==========




    ! smmd: ============ READ ATOMIC VELOCITIES ==========
    call internal_velocities()
    ! smmd: ========== END READ ATOMIC VELOCITIES ========

    ! smmd: ======== READ THERMOSTAT DEFINITION ==========
    call internal_thermostat()
    ! smmd: ====== END READ THERMOSTAT DEFINITION ========

    ! smmd: =========== READ TRANSPORT SETUP =============
    call internal_etrans_setup()
    ! smmd: ========= END READ TRANSPORT SETUP ===========

    ! ars: ======= READ SPECIES - MANDATORY BLOCKS =======
    call internal_species()
    call internal_species_pot()
    ! ars: ===== END READ SPECIES - MANDATORY BLOCKS =====




    ! ars: ======= READ SPECIES - OPTIONAL BLOCKS ========
    call internal_species_cond()
    call internal_species_aux()
    call internal_species_core_wf()
    call internal_species_atomic_set('')
    call internal_species_atomic_set('_cond')
    call internal_species_atomic_set('_aux')
    call internal_species_constraints()
    call internal_species_ngwf_plot()
    call internal_species_ldos_groups()

    ! lpl: NBO block
    call internal_species_write_nbo()

    ! ars: ===== END READ SPECIES - OPTIONAL BLOCKS ======




    ! ars: ============= READ HUBBARD BLOCK ==============
    call internal_hubbard()
    ! ars: =========== END READ HUBBARD BLOCK ============

    ! gibo: ======== READ CONSTRAINED_DFT BLOCK ======
    call internal_cdft()
    ! gibo:  ======== READ CONSTRAINED_DFT BLOCK ======


    ! ars: ======== READ CLASSICAL CHARGES BLOCK =========
    call internal_classical_atoms()
    ! ars: ====== END READ CLASSICAL CHARGES BLOCK =======

    !qoh: VDW parameter override block
    call internal_vdwparam_override()

    ! ars: deallocate global arrays
    if (pub_on_root) then
       deallocate(id_list,stat=ierr)
       call utils_dealloc_check('rundat_blocks_read','id_list',ierr)
    endif
    deallocate(all_species_id, stat=ierr)
    call utils_dealloc_check('rundat_blocks_read','all_species_id',ierr)


#ifdef DEBUG
    ! ars: check data if #DEBUG
    call internal_debug()

    if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving rundat_blocks_read'
#endif

    return

  contains


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!----------------------------    DATA BLOCKS    ----------------------------!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!




    !-----------------------------
    !--------- %block lattice_cart
    !-----------------------------

    subroutine internal_lattice_cart()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_reduce, &
                      esdf_convfac

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_len_unit
      !-------------
      ! ars: buffers
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: real
      real(kind=dp) :: lenconfac
      ! ars: integer
      integer :: num_lines

      ! pdh: allocate buffers
      allocate(dbuf(max(nat*4,9)),stat=ierr) !ddor: Altered for DFT+U
      call utils_alloc_check('internal_lattice_cart','dbuf',ierr)
      dbuf = 0.0_DP

      ! aam: Root node only: read information from
      !      <lattice_cart> block of input file
      if (pub_on_root) then

         if (esdf_block('lattice_cart',num_lines)) then
            ! aam: 3/6/09 added ability to read optional  unit string
            !      ("bohr" or "ang") in first line of lattice_cart block
            !      (default is bohr)
            ! smmd: use esdf_convfac to deal with the optional unit string

            lenconfac=0.0_dp ! jd: kills compiler warning
            if (num_lines == 4) then
               read(block_data(1),'(a)') dummy_len_unit
               read(block_data(2),*) dbuf(1:3)
               read(block_data(3),*) dbuf(4:6)
               read(block_data(4),*) dbuf(7:9)
               lenconfac = esdf_convfac(dummy_len_unit,'bohr')
            else if (num_lines == 3) then
               read(block_data(1),*) dbuf(1:3)
               read(block_data(2),*) dbuf(4:6)
               read(block_data(3),*) dbuf(7:9)
               lenconfac=1.0_dp
            else
               write(stdout,'(/a)') 'Error in internal_lattice_cart: &
                   &malformed lattice_cart vector specification'
               call comms_abort
            end if
         end if
         dbuf(1:9)=lenconfac*dbuf(1:9)
      end if

      ! pdh: broadcast and unpack
      call comms_bcast(pub_root_node_id,dbuf,9)
      a1%x = dbuf(1) ; a1%y = dbuf(2) ; a1%z = dbuf(3)
      a2%x = dbuf(4) ; a2%y = dbuf(5) ; a2%z = dbuf(6)
      a3%x = dbuf(7) ; a3%y = dbuf(8) ; a3%z = dbuf(9)

      ! pdh: deallocate buffers
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_lattice_cart','dbuf',ierr)

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') &
           'DEBUG: Read <lattice_cart> block'
#endif

    end subroutine internal_lattice_cart



    !------------------------------------
    !--------- %block padded_lattice_cart
    !------------------------------------

    subroutine internal_padded_lattice_cart

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_reduce, &
                      esdf_convfac
      use rundat, only: pub_coulomb_cutoff

      implicit none


      !-------------
      ! ars: dummies
      character(len=80) :: dummy_len_unit
      !-------------
      ! ars: buffers
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: real
      real(kind=dp) :: lenconfac
      ! ars: integer
      integer :: num_lines


      ! ndmh: Root node only: read information from
      !      <padded_lattice_cart> block of input file
      if (pub_coulomb_cutoff) then

         ! pdh: allocate buffers
         allocate(dbuf(max(nat*4,9)),stat=ierr) !ddor: Altered for DFT+U
         call utils_alloc_check('internal_padded_lattice_cart','dbuf',ierr)

         if (pub_on_root) then
            lenconfac=1.0_dp ! ars: optional units
            if (esdf_block('padded_lattice_cart',num_lines)) then
               ! ndmh: 2/9/09 added ability to read optional unit string
               !      ("bohr" or "ang") in first line of padded_lattice_cart block
               !      (default is bohr)
               ! smmd: use esdf_convfac to deal with the optional unit string
               if (num_lines == 4) then
                  read(block_data(1),'(a)') dummy_len_unit
                  read(block_data(2),*) dbuf(1:3)
                  read(block_data(3),*) dbuf(4:6)
                  read(block_data(4),*) dbuf(7:9)
                  lenconfac = esdf_convfac(dummy_len_unit,'bohr')
               else if (num_lines == 3) then
                  read(block_data(1),*) dbuf(1:3)
                  read(block_data(2),*) dbuf(4:6)
                  read(block_data(3),*) dbuf(7:9)
                  lenconfac=1.0_dp
               else
                  write(stdout,'(/a)') 'Error in cutoff_coulomb_pad_cell: &
                       &malformed padded cell lattice vector specification'
                  call comms_abort
               end if
            else
               dbuf(1) = -1.0_DP; dbuf(2) =  0.0_DP; dbuf(3) =  0.0_DP
               dbuf(4) =  0.0_DP; dbuf(5) = -1.0_DP; dbuf(6) =  0.0_DP
               dbuf(7) =  0.0_DP; dbuf(8) =  0.0_DP; dbuf(9) = -1.0_DP
            end if
            dbuf(1:9)=lenconfac*dbuf(1:9) ! ars: optional units
         end if

         ! ndmh: broadcast and unpack
         call comms_bcast(pub_root_node_id,dbuf,9)
         a1_pad%x = dbuf(1) ; a1_pad%y = dbuf(2) ; a1_pad%z = dbuf(3)
         a2_pad%x = dbuf(4) ; a2_pad%y = dbuf(5) ; a2_pad%z = dbuf(6)
         a3_pad%x = dbuf(7) ; a3_pad%y = dbuf(8) ; a3_pad%z = dbuf(9)

         ! pdh: allocate buffers
         deallocate(dbuf,stat=ierr) !ddor: Altered for DFT+U
         call utils_dealloc_check('internal_padded_lattice_cart','dbuf',ierr)

#ifdef DEBUG
         if (pub_on_root) write(stdout,'(a)') &
              'DEBUG: Read <padded_lattice_cart> block'
#endif

      end if

    end subroutine internal_padded_lattice_cart


    !------------------------------
    !--------- %block positions_abs
    !------------------------------

    subroutine internal_positions_abs()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_reduce, &
                      esdf_convfac, esdf_line_divide

      implicit none

      !--------------
      ! ars: dummies
      integer :: dummy_nat
      character(len=80) :: dummy_len_unit
      !--------------
      ! ars: buffers
      integer, allocatable :: ibuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      character(len=4),allocatable :: species_id(:)
      !--------------
      ! ars: integers
      integer           :: ishift, row, ii
      ! ars: real
      real(kind=dp)     :: lenconfac
      ! ars: character
      character(len=30) :: struc_tag
      !--------------
      ! smmd: allow for specification of groups of atoms
      integer, allocatable :: group_buf(:)
      integer :: maxp
      character(len=80) :: ctemp
      character(len=80) :: cbuf(4)

      ! check whether tag is provided
      if (present(tag)) then
         struc_tag = 'positions_abs_'//trim(tag)
      else
         struc_tag = 'positions_abs'
      end if


      ! pdh: allocate buffers
      allocate(dbuf(max(nat*4,9)),stat=ierr) !ddor: Altered for DFT+U
      call utils_alloc_check('internal_positions_abs','dbuf',ierr)
      allocate(ibuf(nat*4),stat=ierr)
      call utils_alloc_check('internal_positions_abs','ibuf',ierr)
      allocate(group_buf(nat),stat=ierr)
      call utils_alloc_check('internal_positions_abs','group_buf',ierr)
      allocate(species_id(nat),stat=ierr)
      call utils_alloc_check('internal_positions_abs','species_id',ierr)


      ishift=0 ! qoh: initialise on all processors

      ! aam: Root node only: read absolute cartesian positions from
      !      <positions_abs> block of input file
      if (pub_on_root) then
         lenconfac=1.0_dp
         nsp = 0
         if (esdf_block(struc_tag,dummy_nat)) then
            ! aam: 3/6/09 added ability to read optional  unit string
            !      ("bohr" or "ang") in first line of positions_abs block
            !      (default is bohr)
            ! smmd: use esdf_convfac to deal with the optional unit string
            if (dummy_nat == nat+1) then
               read(block_data(1),'(a)') dummy_len_unit
               lenconfac = esdf_convfac(dummy_len_unit,'bohr')
               ishift=1
            elseif (dummy_nat.ne.nat) then
               write(stdout,'(/a)') 'Error in internal_positions_abs: &
                    &mismatching numbers of atoms in'
               write(stdout,'(a)') 'positions_abs specification'
               call comms_abort
            endif
            if (nat > 0) nsp = 1
            do row=1,nat

               ! smmd: modified in order to allow for groups definition
               maxp = 5
               call esdf_line_divide(maxp,cbuf,block_data(row+ishift)) 
               if (maxp == 4) then
                   read(block_data(row+ishift),*) species_id(row),dbuf(row*3-2:row*3)
                  group_buf(row) = 1
               elseif (maxp == 5) then
                  read(block_data(row+ishift),*) species_id(row),dbuf(row*3-2:row*3),group_buf(row)
               endif         

               do ii=1,4
                  ibuf((row-1)*4+ii) = iachar(species_id(row)(ii:ii))
               end do

               ! aam: count number of distinct species
               sp_count: do ii=1,row-1
                  if (species_id(row) == species_id(ii)) then
                     exit sp_count
                  else
                     if (ii == row-1) nsp = nsp + 1
                  end if
               end do sp_count
            end do
            dbuf(1:3*nat) = lenconfac*dbuf(1:3*nat)
         else
            write(stdout,'(/a)') 'Error in internal_positions_abs: &
                 &<positions_abs> block not found in input file'
            call comms_abort
         end if
      end if

      ! pdh: broadcast and unpack positions
      call comms_bcast(pub_root_node_id,dbuf,nat*3)
      call comms_bcast(pub_root_node_id,group_buf,nat)

      do row=1,nat
         elements(row)%centre%x = dbuf(row*3-2)
         elements(row)%centre%y = dbuf(row*3-1)
         elements(row)%centre%z = dbuf(row*3)
         elements(row)%group_id = group_buf(row)
      end do

      ! pdh: broadcast and unpack species_id's
      call comms_bcast(pub_root_node_id,ibuf)
      do row=1,nat
         do ii=1,4
            elements(row)%species_id(ii:ii) = achar(ibuf((row-1)*4+ii))
         end do
      end do

      ! pdh: deallocate buffers
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_positions_abs','dbuf',ierr)
      deallocate(ibuf,stat=ierr)
      call utils_dealloc_check('internal_positions_abs','ibuf',ierr)
      deallocate(group_buf,stat=ierr)
      call utils_dealloc_check('internal_positions_abs','group_buf',ierr)
      deallocate(species_id,stat=ierr)
      call utils_dealloc_check('internal_positions_abs','species_id',ierr)

      ! pdh: broadcast number of species
      call comms_bcast(pub_root_node_id,nsp)


#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') 'DEBUG: Read <positions_abs> block'
#endif

    end subroutine internal_positions_abs


    !------------------------------
    !--------- %block velocities
    !------------------------------

    subroutine internal_velocities()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_reduce, &
                      esdf_convfac
      use rundat, only: md_init_velocities

      implicit none

      !--------------
      ! smmd: dummies
      integer :: dummy_nat
      character(len=80) :: dummy_vel_unit
      !--------------
      ! smmd: buffers
      real(kind=DP), allocatable :: dbuf(:)
      !--------------
      ! smmd: integers
      integer           :: ishift, row
      ! smmd: real
      real(kind=dp)     :: velconfac

      ! pdh: allocate buffers
      allocate(dbuf(nat*4),stat=ierr)
      call utils_alloc_check('internal_velocities','dbuf',ierr)
      dbuf(:) = 0.d0

      ishift=0 ! qoh: initialise on all processors
      ! aam: Root node only: read absolute cartesian positions from
      !      <positions_abs> block of input file
      if (pub_on_root) then
         velconfac=1.0_dp
         if (esdf_block('velocities',dummy_nat)) then
            ! smmd: use esdf_convfac to deal with the optional unit string
            if (dummy_nat == nat+1) then
               read(block_data(1),'(a)') dummy_vel_unit
               velconfac = esdf_convfac(dummy_vel_unit,'auv')
               ishift=1
            elseif (dummy_nat.ne.nat) then
               write(stdout,'(/a)') 'Error in internal_velocities: &
                    &mismatching numbers of atoms in'
               write(stdout,'(a)') 'velocities specification'
               call comms_abort
            endif
            do row=1,nat
               read(block_data(row+ishift),*) dbuf(row*3-2:row*3)
            end do
            dbuf(1:3*nat) = velconfac*dbuf(1:3*nat)
            md_init_velocities = .false.
         else
            md_init_velocities = .true.
         end if
      end if

      call comms_bcast(pub_root_node_id,md_init_velocities)
      call comms_bcast(pub_root_node_id,dbuf,nat*3)
      do row=1,nat
         elements(row)%ion_velocity(1:3) = dbuf(row*3-2:row*3)
      end do

      ! smmd: deallocate buffers
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_velocities','dbuf',ierr)


#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') 'DEBUG: Read <velocities> block'
#endif

    end subroutine internal_velocities



    !------------------------------
    !--------- %block etrans_setup
    !------------------------------

    subroutine internal_etrans_setup()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_line_divide
      use transport, only: pub_nleads, info_leads, info_device

      implicit none

      !--------------
      integer, allocatable :: ibuf(:)
      character(len=80)    :: line, cjunk(5)
      character(len=1), allocatable :: cbuf(:)
      integer           :: row, nrow, il, maxp
      logical           :: etrans

      ! smmd: allocate buffers
      if (pub_on_root) then
         etrans = esdf_block('etrans_setup',nrow)
         if (etrans .and. nrow .lt. 3) then
            write(stdout,'(/a)') 'Error in internal_etrans_setup: &
                 &incorrect number of lines in'
            write(stdout,'(a)') 'transport setup specification'
         endif
      endif
      call comms_bcast(pub_root_node_id,etrans)
      call comms_bcast(pub_root_node_id,nrow)

      if (etrans) then
         if (nrow .ge. 3) then
            allocate(ibuf(nrow*4),stat=ierr)
            call utils_alloc_check('internal_etrans_setup','ibuf',ierr)
            allocate(cbuf(nrow),stat=ierr)
            call utils_alloc_check('internal_etrans_setup','cbuf',ierr)
         else
            call comms_abort
         endif
         
         ! smmd: Root node only: read transport setup description from 
         !      <etrans_setup> block of input file
         if (pub_on_root) then
            do row = 1, nrow
               line = block_data(row)
         
               maxp = 5 
               call esdf_line_divide(maxp,cjunk,line)
               cbuf(row) = trim(adjustl(cjunk(1)))
               if (cbuf(row)=='D') then
                  do il=2,3
                     read(cjunk(il),*) ibuf((row-1)*4+il-1)
                  enddo
               else
                  do il=2,5
                     read(cjunk(il),*) ibuf((row-1)*4+il-1)
                  enddo
               endif
         
            end do
         end if
         
         call comms_bcast(pub_root_node_id,ibuf,nrow*4)
         call comms_bcast(pub_root_node_id,cbuf,nrow)
         
         pub_nleads = nrow-1
         if (pub_on_root) write(stdout,*) 'TRC setup : nleads = ', pub_nleads
         allocate(info_leads(pub_nleads),stat=ierr)
         call utils_alloc_check('internal_etrans_setup','pub_leads',ierr)
        
         il = 0 
         do row = 1, nrow
            if (cbuf(row).eq.'D') then
               info_device%atms(1:2) = ibuf(row*4-3:row*4-2)
               info_device%type = cbuf(row)
               if (pub_on_root) write(stdout,*) 'TRC setup : Device ==', info_device%atms(1:2), info_device%type
            else
               il = il + 1
               info_leads(il)%atms(1:4) = ibuf(row*4-3:row*4)
               info_leads(il)%type = cbuf(row)
               if (pub_on_root) write(stdout,*) 'TRC setup : Leads ==', il, info_leads(il)%atms(1:2), info_leads(il)%type
            endif
         end do
         
         ! smmd: deallocate buffers
         deallocate(ibuf,stat=ierr)
         call utils_dealloc_check('internal_etrans_setup','ibuf',ierr)
         deallocate(cbuf,stat=ierr)
         call utils_dealloc_check('internal_etrans_setup','cbuf',ierr)
      endif

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') 'DEBUG: Read <etrans_setup> block'
#endif

    end subroutine internal_etrans_setup



    !------------------------------------
    !--------- %block thermostat
    !------------------------------------

    subroutine internal_thermostat()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_convfac, esdf_reduce
      use md_thermostat, only: md_thermo, md_thermo_num
      use rundat, only: md_delta_t

      implicit none

      !-------------
      real(kind=DP)    :: io_temp, io_tgrad
      real(kind=DP)    :: io_freq, io_mix, io_damp
      integer          :: io_kind, io_start, io_stop
      integer          :: io_group, io_nhc, io_nhi
      logical          :: io_upd
      integer          :: nthermo
      real(kind=dp)    :: confac
      !-------------
      integer, allocatable :: ibuf(:)
      logical, allocatable :: lbuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      integer :: row, nrow
      integer :: ith, ierr
      integer :: nhcl
      character (len=80) :: cjunk
      character (len=80) :: cbuf(80)
      character (len=4) :: thermotype(4)
      logical :: thermo
      !-------------
      !data thermotype /'none','andersen','langevin','nosehoover'/

      ! smmd: allocate buffers
      if (pub_on_root) thermo = esdf_block('thermostat',nrow)
      call comms_bcast(pub_root_node_id,nrow)

      if (nrow .ge. 1) then
         allocate(lbuf(nrow),stat=ierr) 
         call utils_alloc_check('internal_thermostat','lbuf',ierr)
         allocate(dbuf(nrow*5),stat=ierr) 
         call utils_alloc_check('internal_thermostat','dbuf',ierr)
         allocate(ibuf(nrow*6),stat=ierr)
         call utils_alloc_check('internal_thermostat','ibuf',ierr)
         dbuf(:) = 0.0_dp
         ibuf(:) = 0
      endif

      ! smmd: Root node only: read information from <thermostat> 
      ! block of input file
      ierr = 0
      if (pub_on_root .and. nrow .ge. 1) then
         if (esdf_block('thermostat',nrow)) then
            nthermo = 0
            row = 1
            do while (row .le. nrow)

               ! Parse each user defined thermostat
               cjunk = block_data(row)

               if (index(cjunk,'none').gt.0 .or. index(cjunk,'andersen').gt.0 .or.&
                   index(cjunk,'langevin').gt.0 .or. index(cjunk,'nosehoover').gt.0) then

                  ! Read mandatory parameters
                  call thermostat_read_params(cjunk, &
                           io_temp,io_kind,io_start,io_stop)
                  nthermo = nthermo + 1
                  row = row + 1

                  ! Default parameters
                  io_tgrad  = 0.0_dp 
                  io_freq   = 0.1_dp/md_delta_t 
                  io_mix    = 1.0_dp
                  io_damp   = 0.2_dp
                  io_group  = 0 
                  io_nhc    = 0
                  io_nhi    = 20
                  io_upd    = .false.

                  ! User defined parameters
                  options_loop : do 
                     cjunk = block_data(row)
                     if (index(cjunk,'none').gt.0 .or. index(cjunk,'andersen').gt.0 .or.&
                         index(cjunk,'langevin').gt.0 .or. index(cjunk,'nosehoover').gt.0 .or.&
                         row .gt. nrow) then
                        exit options_loop
                     else
                        ! Look for optional parameters
                        call thermostat_read_options(cjunk,io_tgrad, &
                                 io_freq,io_mix,io_damp,io_group,io_nhc,io_nhi,io_upd)
                        row = row + 1
                     endif
                  enddo options_loop

               else
                  write(stdout,'(/a)') 'Error in internal_thermostat: &
                      &wrong thermostat definition in input file'
                  call comms_abort
               endif
               
               ! Load buffer arrays
               dbuf((nthermo-1)*5+1:nthermo*5) = &
                   (/io_temp,io_tgrad,io_freq,io_mix,io_damp/)
               ibuf((nthermo-1)*6+1:nthermo*6) = &
                   (/io_kind,io_start,io_stop,io_group,io_nhc,io_nhi/)
               lbuf(nthermo) = io_upd

            enddo
         endif
      endif

      if (nrow .ge. 1) then
         ! smmd: Broadcast buffers
         call comms_bcast(pub_root_node_id,nthermo)
         call comms_bcast(pub_root_node_id,dbuf,nrow*5)
         call comms_bcast(pub_root_node_id,ibuf,nrow*6)
         call comms_bcast(pub_root_node_id,lbuf,nrow)

         ! smmd: allocate public arrays in thermostat_mod 
         md_thermo_num = nthermo
         allocate(md_thermo(md_thermo_num),stat=ierr) 
         call utils_alloc_check('internal_thermostat','md_thermo',ierr)
         
         ! smmd: process thermostat parameters 
         do ith=1,md_thermo_num
         
            md_thermo(ith)%type   = ibuf((ith-1)*6+1)
            md_thermo(ith)%start  = ibuf((ith-1)*6+2)
            md_thermo(ith)%stop   = ibuf((ith-1)*6+3)
            md_thermo(ith)%group  = ibuf((ith-1)*6+4)
            md_thermo(ith)%nhc_length = ibuf((ith-1)*6+5)
            md_thermo(ith)%nhc_integ_nstep = ibuf((ith-1)*6+6)
         
            md_thermo(ith)%Tinit = dbuf((ith-1)*5+1)
            md_thermo(ith)%Tgrad = dbuf((ith-1)*5+2)
            md_thermo(ith)%freq  = dbuf((ith-1)*5+3)
            md_thermo(ith)%mix   = dbuf((ith-1)*5+4)
            md_thermo(ith)%damp  = dbuf((ith-1)*5+5)

            md_thermo(ith)%nhc_upd = lbuf(ith)

         enddo  
         
         ! smmd: deallocate buffers
         deallocate(ibuf,stat=ierr)
         call utils_dealloc_check('internal_thermostat','ibuf',ierr)
         deallocate(dbuf,stat=ierr)
         call utils_dealloc_check('internal_thermostat','dbuf',ierr)
         deallocate(lbuf,stat=ierr)
         call utils_dealloc_check('internal_thermostat','lbuf',ierr)
      endif

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') &
           'DEBUG: Read <thermostat> block'
#endif

    end subroutine internal_thermostat


    !------------------------
    !--------- %block_species
    !------------------------

    subroutine internal_species

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_reduce, &
                      esdf_convfac
      use rundat, only: pub_ngwf_halo

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_len_unit
      character(len=4)  :: dummy_id
      character(len=2)  :: dummy_symb
      integer :: dummy_nsp, dummy_nfunc, dummy_Z
      real(kind=dp) :: dummy_rad
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      integer, allocatable :: ibuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: integers
      integer :: ishift, row, ii
      ! ars: logical
      logical :: element_found
      ! ars: real
      real(kind=dp) :: lenconfac



      ! pdh: allocate buffers
      allocate(dbuf(max(nat*4,9)),stat=ierr) !ddor: Altered for DFT+U
      call utils_alloc_check('internal_species','dbuf',ierr)
      allocate(ibuf(nat*4),stat=ierr)
      call utils_alloc_check('internal_species','ibuf',ierr)
      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species','cbuf',ierr)

      ! aam: Root node only: read information from <species> block of input file
      if (pub_on_root) then
         lenconfac=1.0_dp
         ishift=0
         if (esdf_block('species',dummy_nsp)) then
            ! aam: 3/6/09 added ability to read optional  unit string
            !      ("bohr" or "ang") in first line of species block
            !      (default is bohr)
            ! smmd: use esdf_convfac to deal with the optional unit string
            if (dummy_nsp == nsp+1) then
               read(block_data(1),'(a)') dummy_len_unit
               lenconfac = esdf_convfac(dummy_len_unit,'bohr')
               ishift=1
            endif
            if ( (dummy_nsp.ne.nsp) .and. (dummy_nsp.ne.nsp+1) ) then
               write(stdout,'(/a)') 'Error in internal_species: &
                    &mismatching numbers of species in <species> specification'
               call comms_abort
            endif
            ! aam: check for multiply defined species
            call internal_check_species(nsp,'<species>') !ddor
            do row=1,nsp
               read(block_data(row+ishift),*) dummy_id,dummy_symb, &
                    ibuf(row*2-1:row*2),dbuf(row)
               cbuf(row)(1:4) = dummy_id
               cbuf(row)(5:6) = dummy_symb
               ! ars: save species_id's in buffer
               all_species_id(row)(1:4) = cbuf(row)(1:4)
            end do
            dbuf(1:nsp)=lenconfac*dbuf(1:nsp)
         else
            write(stdout,'(/a)') 'Error in internal_species: &
                 &<species> block not found in input file'
            call comms_abort
         end if
      end if

      ! pdh: broadcast species information
      call comms_bcast(pub_root_node_id,ibuf,nsp*2)
      call comms_bcast(pub_root_node_id,dbuf,nsp)
      do row=1,nsp
         call comms_bcast(pub_root_node_id,cbuf(row),6)
      end do

      ! extract element information from species
      num = 0 ! initialise counter of NGWFs
      element_found = .false. ! qoh: Initialise to avoid compiler warnings
      do row=1,nsp
         element_found = .false. ! aam: error check flag
         dummy_id = cbuf(row)(1:4)
         dummy_symb = cbuf(row)(5:6)
         dummy_Z = ibuf(row*2-1)
         dummy_nfunc = ibuf(row*2)
         dummy_rad = dbuf(row)
         do ii=1,nat
            if (elements(ii)%species_id == dummy_id) then
               elements(ii)%symbol         = dummy_symb
               elements(ii)%atomic_number  = dummy_Z
               elements(ii)%nfunctions     = dummy_nfunc
               elements(ii)%species_number = row
               elements(ii)%radius         = dummy_rad
               num = num + dummy_nfunc  ! aam: count total number of NGWFS
               element_found = .true.
            end if
         end do
         if (.not. element_found) then
            if (pub_on_root) write(stdout,'(/a)') &
                 'Error in internal_species: &
                 &mismatching species in <species> and <positions_abs> &
                 &blocks of input file'
            call comms_abort
         end if
      end do

      ! cks: 21/04/2006 Warning! If there is a halo, the element radius
      ! cks: is the NGWF radius plus the halo length and so it is not
      ! cks: equal to the NGWF-sphere radius (which does not include the halo)
      if (pub_ngwf_halo > 0.0_DP) then
         ! cks: increase element radius by halo length
         do ii=1,nat
            elements(ii)%radius = elements(ii)%radius + pub_ngwf_halo
         end do
      end if

      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_species','cbuf',ierr)
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_species','dbuf',ierr)
      deallocate(ibuf,stat=ierr)
      call utils_dealloc_check('internal_species','ibuf',ierr)


#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') 'DEBUG: Read <species> block'
#endif

    end subroutine internal_species


    !----------------------------
    !--------- %block species_pot
    !----------------------------

    subroutine internal_species_pot()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast, comms_abort
      use esdf, only: block_data, esdf_block

      implicit none

      !-------------
      ! ars: dummies
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      character(len=64) :: dummy_pseudoname
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      !-------------
      ! ars: integers
      integer :: row, ii
      ! ars: logical
      logical :: element_found


      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species_pot','cbuf',ierr)

      ! aam: Root node only: read information from
      !      <species_pot> block of input file
      if (pub_on_root) then
         if (esdf_block('species_pot',dummy_nsp)) then
            if (nsp == dummy_nsp)  then
               ! aam: check for multiply defined species
               call internal_check_species(nsp,'<species_pot>')
               do row=1,nsp
                  read(block_data(row),*) dummy_id, dummy_pseudoname
                  cbuf(row)(1:4) = dummy_id
                  cbuf(row)(5:68) = dummy_pseudoname
               end do
            else
               write(stdout,'(/a)') 'Error in internal_species_pot: &
                    &mismatching numbers of species in'
               write(stdout,'(a)') '<species_pot> block of input file'
               call comms_abort
            end if
         else
            write(stdout,'(/a)') 'Error in internal_species_pot: &
                 &<species_pot> block not found in input file'
            call comms_abort
         end if
      end if

      ! pdh: broadcast species information
      do row=1,nsp
         call comms_bcast(pub_root_node_id,cbuf(row))
      end do

      ! extract element information from species
      do row=1,nsp
         element_found = .false. ! aam: error check flag
         dummy_id = cbuf(row)(1:4)
         dummy_pseudoname = cbuf(row)(5:68)
         do ii=1,nat
            if (elements(ii)%species_id == dummy_id) then
               elements(ii)%pseudo_name = dummy_pseudoname
               element_found = .true.
            end if
         end do
         if (.not. element_found) then
            if (pub_on_root) write(stdout,'(/a)') &
                 'Error in internal_species_pot: &
                 &mismatching species in <species_pot> and <positions_abs> &
                 &blocks of input file'
            call comms_abort
         end if
      end do


      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_species_pot','cbuf',ierr)


#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') 'DEBUG: Read <species_pot> block'
#endif

    end subroutine internal_species_pot


    !----------------------------
    !--------- %block species_core_wf
    !----------------------------

    subroutine internal_species_core_wf()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast, comms_abort
      use esdf, only: block_data, esdf_block

      implicit none

      ! Local Variables
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      character(len=64) :: dummy_corewfname
      character(len=68), allocatable :: cbuf(:)
      integer :: row, ii
      logical :: element_found

      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species_core_wf','cbuf',ierr)

      ! ndmh: Root node only: read information from
      !       <species_core_wf> block of input file
      if (pub_on_root) then
         if (esdf_block('species_core_wf',dummy_nsp)) then
            if (nsp == dummy_nsp)  then
               ! aam: check for multiply defined species
               call internal_check_species(nsp,'<species_core_wf>')
               do row=1,nsp
                  read(block_data(row),*) dummy_id, dummy_corewfname
                  cbuf(row)(1:4) = dummy_id
                  cbuf(row)(5:68) = dummy_corewfname
               end do
            else
               write(stdout,'(/a)') 'Error in internal_species_core_wf: &
                    &mismatching numbers of species in'
               write(stdout,'(a)') '<species_core_wf> block of input file'
               call comms_abort
            end if
         else
            !write(stdout,'(/a/)',advance='no') &
            !     '  <species_core_wf> block not found'
            do row=1,nsp
               cbuf(row)(1:4) = all_species_id(row)(1:4)
               cbuf(row)(5:68) = 'NONE'//repeat(' ',60)
            end do
         end if
      end if

      ! ndmh: broadcast species information
      do row=1,nsp
         call comms_bcast(pub_root_node_id,cbuf(row))
      end do

      ! extract element information from species
      do row=1,nsp
         element_found = .false. ! aam: error check flag
         dummy_id = cbuf(row)(1:4)
         dummy_corewfname = cbuf(row)(5:68)
         do ii=1,nat
            if (elements(ii)%species_id == dummy_id) then
               elements(ii)%core_wf_name = dummy_corewfname
               element_found = .true.
            end if
         end do
         if (.not. element_found) then
            if (pub_on_root) write(stdout,'(/a)') &
                 'Error in internal_species_core_wf: &
                 &mismatching species in <species_core_wf> and <positions_abs> &
                 &blocks of input file'
            call comms_abort
         end if
      end do


      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_species_core_wf','cbuf',ierr)


#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') 'DEBUG: Read <species_core_wf> block'
#endif

    end subroutine internal_species_core_wf


    !-----------------------------
    !--------- %block species_cond
    !-----------------------------

    subroutine internal_species_cond()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_reduce, &
                      esdf_convfac
      use rundat, only: pub_cond_calculate, pub_ngwf_halo

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_len_unit
      character(len=4)  :: dummy_id
      character(len=2)  :: dummy_symb
      integer :: dummy_nsp, dummy_nfunc, dummy_Z
      real(kind=dp) :: dummy_rad
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      integer, allocatable :: ibuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: integers
      integer :: ishift, row, ii
      ! ars: logical
      logical :: element_found
      ! ars: real
      real(kind=dp) :: lenconfac


      if (pub_cond_calculate) then

         ! pdh: allocate buffers
         allocate(dbuf(max(nat*4,9)),stat=ierr) !ddor: Altered for DFT+U
         call utils_alloc_check('internal_species_cond','dbuf',ierr)
         allocate(ibuf(nat*4),stat=ierr)
         call utils_alloc_check('internal_species_cond','ibuf',ierr)
         allocate(cbuf(nsp),stat=ierr)
         call utils_alloc_check('internal_species_cond','cbuf',ierr)


         ! ndmh: Root node only: read information from <species_cond> block of input
         ! ndmh: file for conduction band NGWFs
         ! lr408: If this isn't a conduction calculation, don't bother checking
         ! lr408: for conduction species block
         if (pub_on_root) then
            lenconfac=1.0_dp
            ishift=0
            if (esdf_block('species_cond',dummy_nsp)) then
               ! aam: 3/6/09 added ability to read optional  unit string
               !      ("bohr" or "ang") in first line of species block
               !      (default is bohr)
               ! smmd: use esdf_convfac to deal with the optional unit string
               if (dummy_nsp == nsp+1) then
                  read(block_data(1),'(a)') dummy_len_unit
                  lenconfac = esdf_convfac(dummy_len_unit,'bohr')
                  ishift=1
               endif
               if ( (dummy_nsp.ne.nsp) .and. (dummy_nsp.ne.nsp+1) ) then
                  write(stdout,'(/a)') 'Error in internal_species_cond: &
                       &mismatching numbers of species in'
                  write(stdout,'(a)') '<species_cond> specification'
                  call comms_abort
               endif
               ! aam: check for multiply defined species
               call internal_check_species(nsp,'<species_cond>') !ddor
               do row=1,nsp
                  read(block_data(row+ishift),*) dummy_id,dummy_symb, &
                       ibuf(row*2-1:row*2),dbuf(row)
                  cbuf(row)(1:4) = dummy_id
                  cbuf(row)(5:6) = dummy_symb
               end do
               dbuf(1:nsp)=lenconfac*dbuf(1:nsp)
            else
               write(stdout,'(/a)') 'Error in internal_species_cond: &
                    &conduction calculation specified in task, but '
               write(stdout,'(a)') 'no species_cond block has been specified'
               call comms_abort
            end if
         end if



         ! ndmh: broadcast species information
         call comms_bcast(pub_root_node_id,ibuf,nsp*2)
         call comms_bcast(pub_root_node_id,dbuf,nsp)
         do row=1,nsp
            call comms_bcast(pub_root_node_id,cbuf(row),6)
         end do

         ! ndmh: extract element information from species
         num_cond = 0 ! ndmh: initialise counter of conduction NGWFs
         element_found = .false.
         do row=1,nsp
            element_found = .false.
            dummy_id = cbuf(row)(1:4)
            dummy_symb = cbuf(row)(5:6)
            dummy_Z = ibuf(row*2-1)
            dummy_nfunc = ibuf(row*2)
            dummy_rad = dbuf(row)
            do ii=1,nat
               if (elements(ii)%species_id == dummy_id) then
                  ! ndmh: check Z and symbol match
                  if (elements(ii)%symbol/=dummy_symb) then
                     if (pub_on_root) then
                        write(stdout,'(a)') 'Error in services_read_system&
                             &_input: mismatching element symbols in'
                        write(stdout,'(a)') '<species_cond> and <species> blocks'
                     end if
                     call comms_abort
                  end if
                  if (elements(ii)%atomic_number /= dummy_Z) then
                     if (pub_on_root) then
                        write(stdout,'(a)') 'Error in services_read_system_&
                             &input: mismatching atomic numbers in'
                        write(stdout,'(a)') '<species_cond> and <species> blocks'
                     end if
                     call comms_abort
                  end if
                  ! ndmh: copy info to elements array
                  elements(ii)%nfunctions_cond = dummy_nfunc
                  elements(ii)%radius_cond     = dummy_rad
                  num_cond = num_cond + dummy_nfunc
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               if (pub_on_root) then
                  write(stdout,'(/a)') 'Error in internal_species_cond: &
                       &mismatching species in'
                  write(stdout,'(a)') '<species_cond> and <positions_abs> &
                       &blocks of input file'
               end if
               call comms_abort
            end if
         end do

         ! cks: 21/04/2006 Warning! If there is a halo, the element radius
         ! cks: is the NGWF radius plus the halo length and so it is not
         ! cks: equal to the NGWF-sphere radius (which does not include the halo)
         if (pub_ngwf_halo > 0.0_DP) then
            ! cks: increase element radius by halo length
            do ii=1,nat
               elements(ii)%radius_cond = elements(ii)%radius_cond + pub_ngwf_halo
            end do
         end if

         ! pdh: deallocate buffers
         deallocate(cbuf,stat=ierr)
         call utils_dealloc_check('internal_species_cond','cbuf',ierr)
         deallocate(dbuf,stat=ierr)
         call utils_dealloc_check('internal_species_cond','dbuf',ierr)
         deallocate(ibuf,stat=ierr)
         call utils_dealloc_check('internal_species_cond','ibuf',ierr)

#ifdef DEBUG
         if (pub_on_root) write(stdout,'(a)') 'DEBUG: Read <species_cond> block'
#endif

      else
         ! ndmh: no conduction NGWF species block is present, so set all related
         ! ndmh: entries to zero
         num_cond = 0
         do ii=1,nat
            elements(ii)%radius_cond = 0.0_DP
            elements(ii)%nfunctions_cond = 0
         end do
      end if


    end subroutine internal_species_cond


    !-----------------------------
    !--------- %block species_aux
    !-----------------------------

    subroutine internal_species_aux()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block, esdf_reduce, &
                      esdf_convfac
      use rundat, only: pub_use_aux_ngwfs

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_len_unit
      character(len=4)  :: dummy_id
      character(len=2)  :: dummy_symb
      integer :: dummy_nsp, dummy_nfunc, dummy_Z
      real(kind=dp) :: dummy_rad
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      integer, allocatable :: ibuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: integers
      integer :: ishift, row, ii
      ! ars: logical
      logical :: element_found
      ! ars: real
      real(kind=dp) :: lenconfac

      if (pub_on_root) then
         ! ndmh: activate auxiliary NGWFs if block is present
         pub_use_aux_ngwfs = esdf_block('species_aux',dummy_nsp)
      end if
      call comms_bcast(pub_root_node_id,pub_use_aux_ngwfs)

      if (pub_use_aux_ngwfs) then

         ! pdh: allocate buffers
         allocate(dbuf(max(nat*4,9)),stat=ierr)
         call utils_alloc_check('internal_species_aux','dbuf',ierr)
         allocate(ibuf(nat*4),stat=ierr)
         call utils_alloc_check('internal_species_aux','ibuf',ierr)
         allocate(cbuf(nsp),stat=ierr)
         call utils_alloc_check('internal_species_aux','cbuf',ierr)


         ! ndmh: Root node only: read information from <species_aux> block of input
         ! ndmh: file for auxiliary NGWFs
         if (pub_on_root) then
            lenconfac=1.0_dp
            ishift=0
            if (esdf_block('species_aux',dummy_nsp)) then
               ! aam: 3/6/09 added ability to read optional  unit string
               !      ("bohr" or "ang") in first line of species block
               !      (default is bohr)
               ! smmd: use esdf_convfac to deal with the optional unit string
               if (dummy_nsp == nsp+1) then
                  read(block_data(1),'(a)') dummy_len_unit
                  lenconfac = esdf_convfac(dummy_len_unit,'bohr')
                  ishift=1
               endif
               if ( (dummy_nsp.ne.nsp) .and. (dummy_nsp.ne.nsp+1) ) then
                  write(stdout,'(/a)') 'Error in internal_species_aux: &
                       &mismatching numbers of species in'
                  write(stdout,'(a)') '<species_aux> specification'
                  call comms_abort
               endif
               ! aam: check for multiply defined species
               call internal_check_species(nsp,'<species_aux>') !ddor
               do row=1,nsp
                  read(block_data(row+ishift),*) dummy_id,dummy_symb, &
                       ibuf(row*2-1:row*2),dbuf(row)
                  cbuf(row)(1:4) = dummy_id
                  cbuf(row)(5:6) = dummy_symb
               end do
               dbuf(1:nsp)=lenconfac*dbuf(1:nsp)
            else
               write(stdout,'(/a)') 'Error in internal_species_aux: &
                    &auxiliary NGWFs required due to task, but '
               write(stdout,'(a)') 'no species_aux block has been specified'
               call comms_abort
            end if
         end if



         ! ndmh: broadcast species information
         call comms_bcast(pub_root_node_id,ibuf,nsp*2)
         call comms_bcast(pub_root_node_id,dbuf,nsp)
         do row=1,nsp
            call comms_bcast(pub_root_node_id,cbuf(row),6)
         end do

         ! ndmh: extract element information from species
         num_aux = 0 ! ndmh: initialise counter of auxiliary NGWFs
         element_found = .false.
         do row=1,nsp
            element_found = .false.
            dummy_id = cbuf(row)(1:4)
            dummy_symb = cbuf(row)(5:6)
            dummy_Z = ibuf(row*2-1)
            dummy_nfunc = ibuf(row*2)
            dummy_rad = dbuf(row)
            do ii=1,nat
               if (elements(ii)%species_id == dummy_id) then
                  ! ndmh: check Z and symbol match
                  if (elements(ii)%symbol/=dummy_symb) then
                     if (pub_on_root) then
                        write(stdout,'(a)') 'Error in services_read_system&
                             &_input: mismatching element symbols in'
                        write(stdout,'(a)') '<species_aux> and <species> blocks'
                     end if
                     call comms_abort
                  end if
                  if (elements(ii)%atomic_number /= dummy_Z) then
                     if (pub_on_root) then
                        write(stdout,'(a)') 'Error in services_read_system_&
                             &input: mismatching atomic numbers in'
                        write(stdout,'(a)') '<species_aux> and <species> blocks'
                     end if
                     call comms_abort
                  end if
                  ! ndmh: copy info to elements array
                  elements(ii)%nfunctions_aux = dummy_nfunc
                  num_aux = num_aux + dummy_nfunc
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               if (pub_on_root) then
                  write(stdout,'(/a)') 'Error in internal_species_aux: &
                       &mismatching species in'
                  write(stdout,'(a)') '<species_aux> and <positions_abs> &
                       &blocks of input file'
               end if
               call comms_abort
            end if
         end do

         ! pdh: deallocate buffers
         deallocate(cbuf,stat=ierr)
         call utils_dealloc_check('internal_species_aux','cbuf',ierr)
         deallocate(dbuf,stat=ierr)
         call utils_dealloc_check('internal_species_aux','dbuf',ierr)
         deallocate(ibuf,stat=ierr)
         call utils_dealloc_check('internal_species_aux','ibuf',ierr)

#ifdef DEBUG
         if (pub_on_root) write(stdout,'(a)') 'DEBUG: Read <species_aux> block'
#endif

      else
         ! ndmh: no auxiliary NGWF species block is present, so set all related
         ! ndmh: entries to zero/false
         pub_use_aux_ngwfs = .false.
         num_aux = 0
         do ii=1,nat
            elements(ii)%nfunctions_aux = 0
         end do
      end if


    end subroutine internal_species_aux





    !-----------------------------------
    !--------- %block species_atomic_set
    !-----------------------------------

    subroutine internal_species_atomic_set(which_type)

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block

      implicit none

      ! Arguments
      character(*), intent(in) :: which_type

      !-------------
      ! ars: dummies
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      character(len=64) :: dummy_ngwf_set
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      !-------------
      ! ars: integers
      integer :: row, ii
      ! ars: logical
      logical :: element_found
      logical :: copy_val_ngwfs


      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species_atomic_set','cbuf',ierr)

      ! aam: Root node only: read information from
      !      <species_atomic_set> block of input file
      if (pub_on_root) then
         if (esdf_block('species_atomic_set'//trim(which_type),dummy_nsp)) then
            if (nsp == dummy_nsp)  then
               ! aam: check for multiply defined species
               call internal_check_species(nsp,'<species_atomic_set'// &
                    trim(which_type)//'>')
               do row=1,nsp
                  read(block_data(row),*) dummy_id,dummy_ngwf_set
                  cbuf(row)(1:4) = dummy_id
                  cbuf(row)(5:68) = dummy_ngwf_set
               end do
            else
               write(stdout,'(/a)') 'Error in internal_species_atomic_set: &
                    &mismatching numbers of species in'
               write(stdout,'(a)') '<species_atomic_set'// &
                    trim(which_type)//'> block of input file'
               call comms_abort
            end if
            copy_val_ngwfs = .false.
         else
            copy_val_ngwfs = .true.
            if (trim(which_type)=='') then
               write(stdout,'(/a/)',advance='no') &
                    '  <species_atomic_set> block not found:&
                    & NGWF initialisation set to AUTO'
               ! cks: Use automatic NGWF initialisation if species_atomic_set
               ! cks: is missing
               do row=1,nsp
                  cbuf(row)(1:4) = all_species_id(row)(1:4)
                  cbuf(row)(5:68) = 'AUTO'//repeat(' ',60)
               end do
               copy_val_ngwfs = .false.
            end if
         end if
      end if

      ! pdh: broadcast species information
      call comms_bcast(pub_root_node_id,copy_val_ngwfs)
      do row=1,nsp
         call comms_bcast(pub_root_node_id,cbuf(row))
      end do

      do row=1,nsp
         element_found = .false. ! aam: error check flag
         ! ndmh: for cond/aux NGWFs, if we are copying the NGWF sets from
         ! ndmh: the valence NGWF sets, do so now and exit
         if (copy_val_ngwfs) then
            do ii=1,nat
               if (trim(which_type)=='_cond') then
                  elements(ii)%cond_ngwf_set = elements(ii)%ngwf_set
               else if (trim(which_type)=='_aux') then
                  elements(ii)%aux_ngwf_set = elements(ii)%ngwf_set
               end if
            end do
         else
            dummy_id = cbuf(row)(1:4)
            dummy_ngwf_set = cbuf(row)(5:68)
            do ii=1,nat
               if (elements(ii)%species_id == dummy_id) then
                  if (trim(which_type)=='') then
                     elements(ii)%ngwf_set = dummy_ngwf_set
                  else if (trim(which_type)=='_cond') then
                     elements(ii)%cond_ngwf_set = dummy_ngwf_set
                  else if (trim(which_type)=='_aux') then
                     elements(ii)%aux_ngwf_set = dummy_ngwf_set
                  else
                     write(stdout,'(/a)') 'Error in internal_species_&
                          &atomic_set: unexpected type of NGWF set'
                  end if
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               if (pub_on_root) then
                  write(stdout,'(/a)') 'Error in internal_species_atomic_set: &
                       &mismatching species in'
                  write(stdout,'(a)') '<species_atomic_set> and &
                       &<positions_abs> blocks of input file'
               end if
               call comms_abort
            end if
         end if
      end do

      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_species_atomic_set','cbuf',ierr)

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_atomic_set> block'
#endif

    end subroutine internal_species_atomic_set






    !------------------------------------
    !--------- %block species_constraints
    !------------------------------------

    subroutine internal_species_constraints()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block

      implicit none


      !-------------
      ! ars: dummies
      character(len=5)  :: dummy_contype
      real(kind=dp)     :: dummy_convec(3)
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: integers
      integer :: row, ii, ierr
      ! ars: logical
      logical :: element_found



      ! pdh: allocate buffers
      allocate(dbuf(max(nat*4,9)),stat=ierr) !ddor: Altered for DFT+U
      call utils_alloc_check('internal_species_constraints','dbuf',ierr)
      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species_constraints','cbuf',ierr)

      ! aam: Root node only: read information from
      !      <species_constraints> block of input file
      ! pdh: set error flag
      ierr = 0
      if (pub_on_root) then
         if (esdf_block('species_constraints',dummy_nsp)) then
            if (dummy_nsp > nsp) then ! aam: error check
               write(stdout,'(/a)') 'Error in internal_species_constraints: &
                    &too many species in'
               write(stdout,'(a)') '<species_constraints> block &
                    &of input file'
               ierr = 1
            else
               ! aam: check for multiply defined species
               call internal_check_species(dummy_nsp,'<species_constraints>')
               do row=1,dummy_nsp ! aam: dummy_nsp is not necessarily equal to nsp
                  if (index(block_data(row),'FIXED') /= 0 .or. &
                       index(block_data(row),'NONE') /= 0) then
                     ! aam: only species id and constraint type required here
                     read(block_data(row),*) dummy_id, dummy_contype
                     dbuf(row*3-2:row*3) = 0.0_DP
                  else if (index(block_data(row),'LINE') /= 0 .or. &
                       index(block_data(row),'PLANE') /= 0) then
                     ! aam: also require the constraint vector here
                     read(block_data(row),*) dummy_id, dummy_contype, &
                          dbuf(row*3-2:row*3)
                  else
                     write(stdout,'(/a)') 'Error in internal_species_constraints: &
                          &unknown constraint identifier in'
                     write(stdout,'(a)') '<species_constraints> block of input &
                          &file'
                     ierr = 2
                     exit
                  end if
                  cbuf(row)(1:4) = dummy_id
                  cbuf(row)(5:9) = dummy_contype
               end do
            end if
         end if
      end if

      ! pdh: broadcast error flag for a clean shutdown if necessary
      call comms_bcast(pub_root_node_id,ierr)
      if (ierr /= 0) call comms_abort

      ! pdh: broadcast
      call comms_bcast(pub_root_node_id,dummy_nsp)
      do row=1,dummy_nsp
         call comms_bcast(pub_root_node_id,cbuf(row),9)
      end do
      call comms_bcast(pub_root_node_id,dbuf,dummy_nsp*3)

      ! aam: initialise to default values
      do ii=1,nat
         elements(ii)%ion_constraint_type = 'NONE'
         elements(ii)%ion_constraint      = 0.0_DP
      end do

      ! pdh: extract constraint information
      do row=1,dummy_nsp
         element_found = .false. ! aam: error check flag
         dummy_id = cbuf(row)(1:4)
         dummy_contype = cbuf(row)(5:9)
         dummy_convec = dbuf(row*3-2:row*3)
         do ii=1,nat
            if (elements(ii)%species_id == dummy_id) then
               elements(ii)%ion_constraint_type = dummy_contype
               elements(ii)%ion_constraint      = dummy_convec
               element_found = .true.
            end if
         end do
         if (.not. element_found) then
            if (pub_on_root) then
               write(stdout,'(/a)') 'Error in internal_species_constraints: &
                    &mismatching species in'
               write(stdout,'(a)') '<species_constraints> and &
                    &<positions_abs> blocks of input file'
            end if
            call comms_abort
         end if
      end do

      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_species_constraints','cbuf',ierr)
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_species_constraints','dbuf',ierr)

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_constraints> block'
#endif


    end subroutine internal_species_constraints



    !----------------------------------
    !--------- %block_species_ngwf_plot
    !----------------------------------

    subroutine internal_species_ngwf_plot()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block

      implicit none

      !-------------
      ! ars: dummies
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      !-------------
      ! ars: integers
      integer ::  row, ii, ierr
      ! ars: logical
      logical :: element_found



      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_species_ngwf_plot','cbuf',ierr)

      ! cks: Root node only: read information from
      !      <species_ngwf_plot> block of input file
      ! pdh: set error flag
      ierr = 0
      if (pub_on_root) then
         if (esdf_block('species_ngwf_plot',dummy_nsp)) then
            if (dummy_nsp > nsp) then ! aam: error check
               write(stdout,'(/a)') 'Error in internal_species_ngwf_plot: &
                    &too many species in <species_ngwf_plot>'
               write(stdout,'(a)') 'block of input file'
               ierr = 1
            else
               ! aam: check for multiply defined species
               call internal_check_species(dummy_nsp,'<species_ngwf_plot>')
               do row=1,dummy_nsp ! aam: dummy_nsp is not necessarily equal to nsp
                  read(block_data(row),*) dummy_id
                  cbuf(row)(1:4) = dummy_id
               end do
            end if
         end if
      end if

      ! pdh: broadcast error flag for a clean shutdown if necessary
      call comms_bcast(pub_root_node_id,ierr)
      if (ierr /= 0) call comms_abort

      ! pdh: broadcast
      call comms_bcast(pub_root_node_id,dummy_nsp)
      do row=1,dummy_nsp
         call comms_bcast(pub_root_node_id,cbuf(row),4)
      end do

      ! cks: initialise to default values
      do ii=1,nat
         elements(ii)%ngwf_plot = .false.
      end do

      ! pdh: extract information
      do row=1,dummy_nsp
         element_found = .false. ! aam: error check flag
         dummy_id = cbuf(row)(1:4)
         do ii=1,nat
            if (elements(ii)%species_id == dummy_id) then
               elements(ii)%ngwf_plot = .true.
               element_found = .true.
            end if
         end do
         if (.not. element_found) then
            if (pub_on_root) then
               write(stdout,'(/a)') 'Error in internal_species_ngwf_plot: &
                    &mismatching species in'
               write(stdout,'(a)') '<species_ngwf_plot> and <positions_abs> &
                    &blocks of input file'
            end if
            call comms_abort
         endif
      end do

      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_species_ngwf_plot','cbuf',ierr)

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_ngwf_plot> block'
#endif

    end subroutine internal_species_ngwf_plot


    !------------------------------------
    !--------- %block species_ldos_groups
    !------------------------------------

    subroutine internal_species_ldos_groups()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block
      use rundat, only: pub_ldos_ngroups, pub_ldos_group_nsp, pub_ldos_groups

      implicit none

      !-------------
      ! ars: dummies
      character(len=80) :: dummy_ldos_group
      character(len=4)  :: dummy_id
      integer :: dummy_nsp
      !-------------
      ! ars: integers
      integer :: row, ii, jj
      ! ars: logical
      logical :: element_found



      ! ndmh: Root node only: read information from
      !      <species_ldos_groups> block of input file
      ! ndmh: check if the block exists
      if (pub_on_root) then
         if (esdf_block('species_ldos_groups',pub_ldos_ngroups)) then
            if (pub_ldos_ngroups > nat) then
               write(stdout,'(/a)') 'Error in internal_species_ldos_groups: &
                    &more groups than atoms defined in <species_ldos_groups>'
               call comms_abort
            end if
         else
            pub_ldos_ngroups = 0
            dummy_nsp = 0
         end if
      end if

      ! ndmh: allocate array for numbers of species in each group
      call comms_bcast(pub_root_node_id,pub_ldos_ngroups)
      if (allocated(pub_ldos_group_nsp)) then
         deallocate(pub_ldos_group_nsp,stat=ierr)
         call utils_dealloc_check('internal_species_ldos_groups', &
              'pub_ldos_group_nsp',ierr)
      end if
      allocate(pub_ldos_group_nsp(pub_ldos_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_ldos_groups', &
           'pub_ldos_group_nsp',ierr)

      ! ndmh: count number of species in each group
      if (pub_on_root.and.(pub_ldos_ngroups > 0)) then
         dummy_nsp = 0
         ! loop over lines in block (LDOS groups)
         do row=1,pub_ldos_ngroups
            dummy_ldos_group = trim(adjustl(block_data(row)))
            pub_ldos_group_nsp(row) = 1
            ! ndmh: count number of spaces (distinct species) in this group
            jj = 0
            do ii=1,len_trim(dummy_ldos_group)
               if ((dummy_ldos_group(ii:ii)==' ').and. &
                    (dummy_ldos_group(ii+1:ii+1)/=' ')) then
                  if (ii > jj + 1) then
                     pub_ldos_group_nsp(row) = pub_ldos_group_nsp(row) + 1
                     jj = ii
                  else
                     jj = ii
                  end if
               end if
            end do
         end do
      end if

      ! ndmh: allocate array for storing LDOS group species ID's
      call comms_bcast(pub_root_node_id,pub_ldos_group_nsp)
      dummy_nsp = maxval(pub_ldos_group_nsp)
      call comms_bcast(pub_root_node_id,dummy_nsp)
      if (allocated(pub_ldos_groups)) then
         deallocate(pub_ldos_groups,stat=ierr)
         call utils_dealloc_check('internal_species_ldos_groups','pub_ldos_groups',&
              ierr)
      end if
      allocate(pub_ldos_groups(dummy_nsp,pub_ldos_ngroups),stat=ierr)
      call utils_alloc_check('internal_species_ldos_groups','pub_ldos_groups',ierr)
      pub_ldos_groups(:,:) = ''

      ! ndmh: now fill the pub_ldos_groups array with these strings
      if (pub_on_root.and.(pub_ldos_ngroups > 0)) then
         do row=1,pub_ldos_ngroups
            dummy_ldos_group = trim(adjustl(block_data(row)))
            ii = 1
            jj = len_trim(dummy_ldos_group)
            do dummy_nsp=1,pub_ldos_group_nsp(row)
               if (index(dummy_ldos_group(ii:jj),' ') > 0) then
                  jj = ii + index(dummy_ldos_group(ii:jj),' ') - 1
                  do
                     if (dummy_ldos_group(jj+1:jj+1)/=' ') exit
                     jj = jj + 1
                  end do
               end if
               pub_ldos_groups(dummy_nsp,row) = &
                    adjustl(trim(dummy_ldos_group(ii:jj)))
               ii = jj + 1
               jj = len_trim(dummy_ldos_group)
            end do
         end do
      end if

      ! ndmh: broadcast LDOS groups from root node
      do row=1,pub_ldos_ngroups
         do dummy_nsp=1,pub_ldos_group_nsp(row)
            call comms_bcast(pub_root_node_id,pub_ldos_groups(dummy_nsp,row))
         end do
      end do

      ! ndmh: check LDOS group species match elements array species id's
      do row=1,pub_ldos_ngroups
         do ii=1,pub_ldos_group_nsp(row)
            element_found = .false.
            dummy_id = pub_ldos_groups(ii,row)
            do jj=1,nat
               if (elements(jj)%species_id == dummy_id) then
                  element_found = .true.
               end if
            end do
            if (.not. element_found) then
               write(stdout,'(/a)') 'Error in internal_species_ldos_groups: &
                    &mismatching species in'
               write(stdout,'(a)') '<species_ldos_groups> and &
                    &<positions_abs> blocks of input file'
               call comms_abort
            endif
         end do
      end do


#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') &
           'DEBUG: Read <species_ldos_groups> block'
#endif

    end subroutine internal_species_ldos_groups


    !------------------------
    !--------- %block hubbard
    !------------------------

    subroutine internal_hubbard()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use constants, only: HARTREE_IN_EVS
      use esdf, only: block_data, esdf_block
      use hubbard_init, only: hubbard_init_species, h_species
      use rundat, only: pub_hubbard, pub_hubbard_atomsolve, pub_hubbard_restart
      use simulation_cell, only: pub_cell

      implicit none

      !-------------
      ! ars: dummies
      character(len=4)  :: dummy_id
      integer           :: dummy_l, dummy_h
      real(kind=dp)     :: dummy_u, dummy_c, dummy_a, dummy_s !ddor
      integer           :: dummy_hub_nsp !ddor
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      integer, allocatable :: ibuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: integers
      integer :: row, ii
      ! ars: logical
      logical :: element_found



      ! pdh: allocate buffers
      allocate(dbuf(max(nat*4,9)),stat=ierr) !ddor: Altered for DFT+U
      call utils_alloc_check('internal_hubbard','dbuf',ierr)
      allocate(ibuf(nat*4),stat=ierr)
      call utils_alloc_check('internal_hubbard','ibuf',ierr)
      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_hubbard','cbuf',ierr)

      ! ddor: Root node only: read DFT+U information from <hubbard> block of
      ! ddor: input file
      hub_nsp = 0 ! Number of Hubbard species in cell
      pub_hubbard = .false.
      element_found = .false. ! error check flag
      if (pub_on_root) then
         if (esdf_block('hubbard',dummy_hub_nsp)) then
            if ( (dummy_hub_nsp .le. nsp) .and. ( dummy_hub_nsp > 0 ) ) then
               ! aam: check for multiply defined species
               call internal_check_species(dummy_hub_nsp,'<hubbard>')
               hub_nsp = dummy_hub_nsp ! The number of Hubbard species
               pub_hubbard = .true. ! Turn on DFT+U if Hubbard atom found
               do row=1,dummy_hub_nsp
                  read(block_data(row),*) dummy_id, &
                       ibuf(row),dbuf(row*4-3:row*4)
                  cbuf(row)(1:4) = dummy_id
               end do
            else
               write(stdout,'(/a)') 'Error in internal_hubbard: &
                    &mismatching number of species in'
               write(stdout,'(a)') '<hubbard> block of input file'
               call comms_abort
            end if
         else
            if (pub_hubbard) then
               write(stdout,'(/a)') 'Error in internal_hubbard: &
                    &<hubbard> block not found in input'
               write(stdout,'(a)') 'file yet DFT+U feature is activated'
               call comms_abort
            endif
         end if
      end if

      ! ddor: broadcast flag for DFT+U calculation
      call comms_bcast(pub_root_node_id,pub_hubbard)
      ! ddor: broadcast number of Hubbard species
      call comms_bcast(pub_root_node_id,hub_nsp)
      pub_cell%num_hub_species = hub_nsp

      if (pub_hubbard) then

         ! ddor: broadcast species information
         call comms_bcast(pub_root_node_id,ibuf,hub_nsp)
         call comms_bcast(pub_root_node_id,dbuf,hub_nsp*4)
         do row=1,hub_nsp
            call comms_bcast(pub_root_node_id,cbuf(row),4)
         end do

         ! ddor: allocate public DFT+U h_species type
         call hubbard_init_species

         ! ddor: extract Hubbard element information to h_species
         dummy_h = 0 ! initialise counter of Hubbard atoms
         num_hub_proj = 0 ! initialise counter of Hubbard projectors
         ! Set some defaults for non-Hubbard atoms
         do row=1,hub_nsp
            dummy_id = cbuf(row)(1:4)
            dummy_l = ibuf(row)
            dummy_u = dbuf(row*4-3)
            dummy_c = dbuf(row*4-2)
            dummy_a = dbuf(row*4-1)
            dummy_s = dbuf(row*4)
            do ii=1,nat ! The total number of atoms in the cell
               if (elements(ii)%species_id == dummy_id) then
                  ! If it's a Hubbard atom then...
                  ! Store namelist of Hubbard atoms
                  h_species(row)%hub_species = dummy_id(1:4)
                  ! Angular momentum channel of projector
                  h_species(row)%hub_ang_mom = dummy_l
                  ! Hubbard U parameter (eV)
                  h_species(row)%hub_u       = dummy_u / HARTREE_IN_EVS
                  ! Effective charge for radial function
                  h_species(row)%hub_charge  = ABS(dummy_c)
                  ! Hubbard alpha parameter (eV)
                  h_species(row)%hub_alpha   = dummy_a / HARTREE_IN_EVS
                  ! Hubbard spin Zeeman splitting (eV)
                  h_species(row)%hub_spin_splitting = dummy_s / HARTREE_IN_EVS
                  ! ddor: Trigger fireball NGWFs by using negative Z parameter
                  if ( dummy_c < 0.0_DP ) then
                     pub_hubbard_atomsolve = .true.
                  endif

                  if ( ( h_species(row)%hub_ang_mom < 0 ) .or. &
                       &( h_species(row)%hub_ang_mom > 3 ) ) then
                     if (pub_on_root) then
                        write(stdout,'(/a)') 'Error in &
                             &internal_hubbard: unphysical value of '
                        write(stdout,'(a)') 'angular momentum for projector in &
                             &<hubbard> block'
                     end if
                     call comms_abort
                  endif

                  ! ddor: count total number of Hubbard atoms
                  dummy_h = dummy_h + 1
                  num_hub_proj = num_hub_proj + 2 * h_species(row)%hub_ang_mom + 1
                  element_found = .true.

               end if
            end do

           if (pub_hubbard_restart .and. pub_hubbard_atomsolve) then
               if (pub_on_root) then
                  write(stdout,'(/a)') 'Error in &
                       &internal_hubbard: Cannot have negative values '
                  write(stdout,'(a)') 'for DFT+U charge Z with PAW or &
                       &negative projector mixing parameter in &
                       &<hubbard> block'
                end if
                call comms_abort
            endif

            call comms_bcast(pub_root_node_id, pub_hubbard_atomsolve)
            call comms_bcast(pub_root_node_id, h_species(row)%hub_species)
            call comms_bcast(pub_root_node_id, h_species(row)%hub_ang_mom)
            call comms_bcast(pub_root_node_id, h_species(row)%hub_u)
            call comms_bcast(pub_root_node_id, h_species(row)%hub_charge)
            call comms_bcast(pub_root_node_id, h_species(row)%hub_alpha)
            call comms_bcast(pub_root_node_id, &
                 & h_species(row)%hub_spin_splitting)
         end do

         if (.not. element_found) then
            if (pub_on_root) write(stdout,'(/a)') &
                 'Error in internal_hubbard: &
                 &mismatching species in <hubbard> and <positions_abs> &
                 &blocks of input file'
            call comms_abort
         else
            pub_cell%nat_hub = dummy_h ! The number of Hubbard atoms in the cell
            pub_cell%num_hub_proj = num_hub_proj
         end if

         ! ddor: broadcast number of Hubbard atoms in the cell
         call comms_bcast(pub_root_node_id,pub_cell%nat_hub)
         call comms_bcast(pub_root_node_id,pub_cell%num_hub_proj)

#ifdef DEBUG
         if (pub_on_root) write(stdout,'(a)') 'DEBUG: Read <hubbard> block'
#endif
      end if

      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_hubbard','cbuf',ierr)
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_hubbard','dbuf',ierr)
      deallocate(ibuf,stat=ierr)
      call utils_dealloc_check('internal_hubbard','ibuf',ierr)

    end subroutine internal_hubbard


    !------------------------
    !--------- %block constrained_dft
    !------------------------

! IDEA: to (internally) parse the constrained_dft block as a modified version of 
! hubbard_init (written by ddor)
!
    subroutine internal_cdft()

      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use constants, only: HARTREE_IN_EVS
      use esdf, only: block_data, esdf_block
      use hubbard_init, only: hubbard_init_species, h_species
      use rundat, only: pub_hubbard, pub_hubbard_atomsolve, pub_hubbard_restart, &
           pub_cdft,pub_cdft_atom_charge, pub_cdft_atom_spin, &
           pub_cdft_group_charge, pub_cdft_group_spin, &
           pub_cdft_group_charge_diff, pub_cdft_group_spin_diff, &
           pub_cdft_group_u, pub_cdft_group_diff_u, &
           pub_spin_polarised, pub_cdft_guru
      use simulation_cell, only: pub_cell

      implicit none

      !-------------
      ! ars: dummies
      character(len=4)  :: dummy_id
      integer           :: dummy_l, dummy_h
      real(kind=dp)     :: dummy_u, dummy_c, dummy_a, dummy_s
      integer           :: dummy_hub_nsp
      real(kind=dp)     :: dummy_u_q, dummy_u_s, dummy_s_target
      real(kind=dp)     :: dummy_u_q_up, dummy_u_q_down, dummy_n_up, dummy_n_down
      integer           :: dummy_cdft_nsp
      !-------------
      ! ars: buffers
      character(len=68), allocatable :: cbuf(:)
      integer, allocatable :: ibuf(:)
      real(kind=DP), allocatable :: dbuf(:)
      !-------------
      ! ars: integers
      integer :: row, ii
      ! ars: logical
      logical :: element_found
      ! gibo: smallest representable real
      real(KIND=DP) :: eps


      ! pdh: allocate buffers
      allocate(dbuf(max(nat*8,9)),stat=ierr) ! gibo: modified for cDFT 04.10.11
      call utils_alloc_check('internal_cdft','dbuf',ierr)
      allocate(ibuf(nat*8),stat=ierr)  ! gibo: modified for cDFT 04.10.11
      call utils_alloc_check('internal_cdft','ibuf',ierr)
      allocate(cbuf(nsp),stat=ierr)
      call utils_alloc_check('internal_cdft','cbuf',ierr)

      ! Root node only: read DFT+U information from <constrained_dft> block of input file
      element_found = .false. ! error check flag

      cdft_nsp = 0 ! Number of cDFT(=Hubbard) species in cell
      pub_cdft = .false.

      if (pub_on_root) then
         if (esdf_block('constrained_dft',dummy_cdft_nsp)) then
            if ( (dummy_cdft_nsp .le. nsp) .and. ( dummy_cdft_nsp > 0 ) ) then
               ! aam: check for multiply defined species
               call internal_check_species(dummy_cdft_nsp,'<constrained_dft>')
               cdft_nsp = dummy_cdft_nsp ! The number of constrained species
               pub_cdft = .true. ! Turn on cDFT

              !gibo: make sure the <hubbard> and <constrained_cdft> blocks 
              !      are *NOT* simultanously defined in input file
              !MIND: internal_hubbard is called before internal_cdft therefore
              !      pub_hubbard=T if <hubbard> block in .dat file...
              if ( (pub_cdft) .AND. (pub_hubbard) ) then
               write(stdout,'(/a)') 'Error in internal_cdft:'
               write(stdout,'(a)') '<hubbard> and <constrained_dft> blocks &
                                 &cannot be simultanously present in input file'
               write(stdout,'(a)') 're-run with either <hubbard> or &
                                 &<constrained_dft> block in input file'
               call comms_abort
              endif

            ! *** CRUCIAL **** if (pub_cdft) then turn on also DFT+U...
               hub_nsp = cdft_nsp ! number of Hub species=number of cDFT species
               pub_hubbard = .true. ! Turn on also DFT+U 
            ! *** CRUCIAL **** if (pub_cdft) then turn on also DFT+U...

               do row=1,dummy_cdft_nsp
                  read(block_data(row),*) dummy_id, &
                       ibuf(row),dbuf(row*8-7:row*8)  ! gibo: modified 04.10.11
                  cbuf(row)(1:4) = dummy_id
               end do
            else
               write(stdout,'(/a)') 'Error in internal_cdft: &
                    &mismatching number of species in'
               write(stdout,'(a)') '<constrained_dft> block of input file'
               call comms_abort
            end if
         else
            if (pub_cdft) then
               write(stdout,'(/a)') 'Error in internal_cdft: &
                    &<constrained_dft> block not found in input'
               write(stdout,'(a)') 'file yet cDFT feature is activated'
               write(stdout,'(a)') '[pub_cdft=.TRUE.]'
               call comms_abort
            endif
         end if
      end if

      ! gibo: broadcast flag for cDFT calculation
       call comms_bcast(pub_root_node_id,pub_cdft)

      ! gibo: to avoid messing with pure Hubbard calculations...
      if (pub_cdft) then 
         ! ddor: broadcast flag for DFT+U calculation
         call comms_bcast(pub_root_node_id,pub_hubbard)
         ! ddor: broadcast number of Hubbard species
         call comms_bcast(pub_root_node_id,hub_nsp)
         pub_cell%num_hub_species  = hub_nsp
         call comms_bcast(pub_root_node_id,cdft_nsp)
         pub_cell%num_cdft_species = cdft_nsp
      endif 

      ! gibo: follow ddor's DFT+U framework as much as possible... 
      CHECK_CDFT: if (pub_cdft) then

         ! If required, prepare to (hopefully) help unexperienced users
         if (.not.pub_cdft_guru) then 
            eps = epsilon(1.0_DP) 
            if (pub_on_root) then
               write(stdout,'(/a)') 'WARNING in &
                    &internal_cdft: '
               write(stdout,'(a)') 'Fail-safe initialisation of &
                    & cDFT (non-zero) |U-potentials| to 1 eV'
               write(stdout,'(a)') 'If unhappy about this, &
                    &and you really know what you are doing, &
                    &set "cdft_guru: T" in input file'
            end if
         endif
         ! gibo: turn spin-polarisation on for constrained-DFT runs
         pub_spin_polarised = .true.
         call comms_bcast(pub_root_node_id, pub_spin_polarised)

         if (pub_on_root) write(stdout,'(a)') &
         'WARNING in internal_cdft: turn spin-polarisation on for &
          &constrained-DFT run'

         ! ddor: broadcast species information
         call comms_bcast(pub_root_node_id,ibuf,hub_nsp)
         call comms_bcast(pub_root_node_id,dbuf,hub_nsp*8)
         do row=1,hub_nsp
            call comms_bcast(pub_root_node_id,cbuf(row),8)
         end do

         ! ddor: allocate public DFT+U h_species type
         call hubbard_init_species

         ! ddor: extract Hubbard element information to h_species
         dummy_h = 0 ! initialise counter of Hubbard atoms
         num_hub_proj = 0 ! initialise counter of Hubbard projectors
         ! Set some defaults for non-Hubbard atoms
         do row=1,hub_nsp
            dummy_id = cbuf(row)(1:4)
            dummy_l = ibuf(row)

           ! different U_q for alpha and beta spin-channels 04.10.11
            dummy_c         = dbuf(row*8-7)
            dummy_u         = dbuf(row*8-6)
            dummy_u_q_up    = dbuf(row*8-5)
            dummy_u_q_down  = dbuf(row*8-4)
            dummy_u_s       = dbuf(row*8-3)
            dummy_n_up      = dbuf(row*8-2)
            dummy_n_down    = dbuf(row*8-1)
            dummy_s_target  = dbuf(row*8)
           ! different U_q for alpha and beta spin-channels 04.10.11

           ! Unless the user claims to be a cDFT guru [pub_cdft_guru = .TRUE.]
           ! [who should know about the best cDFT U-potentials initialisation]
           ! initialise the (non-zero) cDFT U-potentials to 1 eV
            if (.not.pub_cdft_guru) then
                if (ABS(dummy_u_q_up) > eps) &
                        dummy_u_q_up = dummy_u_q_up/ABS(dummy_u_q_up)
                if (ABS(dummy_u_q_down) > eps) &
                        dummy_u_q_down = dummy_u_q_down/ABS(dummy_u_q_down)
                if (ABS(dummy_u_s) > eps) &
                        dummy_u_s = dummy_u_s/ABS(dummy_u_s)
            endif

            ATOMS_CELL: do ii=1,nat ! The total number of atoms in the cell
               if (elements(ii)%species_id == dummy_id) then
                  ! gibo: If it is a cDFT (=Hubbard) atom then...
                  ! Store namelist of Hubbard atoms
                  h_species(row)%hub_species = dummy_id(1:4)

                ! modified ==== start
                  ! Angular momentum channel of projector
                  h_species(row)%cdft_ang_mom        = dummy_l
                  ! Copy cdft angular momentum into DFT-U counterpart
                  h_species(row)%hub_ang_mom         = dummy_l 

                  ! Effective charge for radial function
                  h_species(row)%hub_charge          = ABS(dummy_c)

                  ! Hubbard U parameter (eV)
                  h_species(row)%hub_u               = dummy_u / HARTREE_IN_EVS

                ! different U_q for alpha and beta spin-channels 04.10.11
                  h_species(row)%cdft_u_charge_up    = dummy_u_q_up / HARTREE_IN_EVS

                  h_species(row)%cdft_u_charge_down  = dummy_u_q_down / HARTREE_IN_EVS
                ! different U_q for alpha and beta spin-channels 04.10.11

                  ! SPIN constrain (U) parameter (eV)
                  h_species(row)%cdft_u_spin         = dummy_u_s / HARTREE_IN_EVS

                  ! TARGET: number of electrons_UP
                  h_species(row)%cdft_target_up      = dummy_n_up

                  ! TARGET: number of electrons_DOWN
                  h_species(row)%cdft_target_down    = dummy_n_down

                  ! TARGET: spin population (N_alpha-N_beta)
                  h_species(row)%cdft_target_spin    = dummy_s_target 

                  ! pub_hubbard_atomsolve = .true.
                  ! ddor: Trigger fireball NGWFs by using negative Z parameter
                  if ( dummy_c < 0.0_DP ) then
                     pub_hubbard_atomsolve = .true.
                  endif
                ! modified ==== end

                  if ( ( h_species(row)%hub_ang_mom < 0 ) .or. &
                       &( h_species(row)%hub_ang_mom > 3 ) ) then
                     if (pub_on_root) then
                        write(stdout,'(/a)') 'Error in &
                             &internal_cdft: unphysical value of '
                        write(stdout,'(a)') 'angular momentum for projector in &
                             &<constrained_dft> block'
                     end if
                     call comms_abort
                  endif

                  ! If constraining group CHARGE difference, 
                  ! decide whether a given atom is an acceptor or a donor...
                  if (pub_cdft_group_charge_diff) then

                   if (dummy_u_q_up > 0._DP) then
                     h_species(row)%cdft_donor    = .TRUE.
                     h_species(row)%cdft_acceptor = .FALSE.

                   elseif (dummy_u_q_up < 0._DP) then
                     h_species(row)%cdft_donor    = .FALSE.
                     h_species(row)%cdft_acceptor = .TRUE.

                   endif

                  endif

                  ! If constraining group CHARGE (difference),
                  ! avoid silly inputs like (U_up > 0, U_down <0 and 
                  ! cdft_group_charge/spin_[diff] active)
                  if (pub_cdft_group_charge_diff .OR. pub_cdft_group_charge) then

                   if ( dummy_u_q_up*dummy_u_q_down < 0._DP) then
                     if (pub_on_root) then
                        write(stdout,'(/a)') 'Error in &
                             &internal_cdft: '
                        write(stdout,'(a)') 'cdft_group_charge[_diff] active, &
                             &yet different signs for U_up and U_down of the &
                             &same atom in <constrained_dft> block'
                     end if
                     call comms_abort
                   endif

                  endif
                  ! gibo. decide whether given atom is an acceptor or a donor===END

                  ! If constraining group SPIN difference, 
                  ! decide whether a given atom is an acceptor or a donor...
                  if (pub_cdft_group_spin_diff) then

                   if (dummy_u_s > 0._DP) then
                     h_species(row)%cdft_donor    = .TRUE.
                     h_species(row)%cdft_acceptor = .FALSE.

                   elseif (dummy_u_s < 0._DP) then
                     h_species(row)%cdft_donor    = .FALSE.
                     h_species(row)%cdft_acceptor = .TRUE.

                   endif

                  endif

                  ! ddor: count total number of Hubbard atoms
                  ! gibo: number of Hubbard atoms = number of cDFT atoms
                  dummy_h = dummy_h + 1
                  num_hub_proj = num_hub_proj + 2 * h_species(row)%hub_ang_mom + 1
                  element_found = .true.

               end if

            enddo ATOMS_CELL

           ! 4 GIBO: TO BE CHECKED AFTER PAW GO AHEAD FROM NICK ==== START
           !: cDFT within PAW is something we fancy...
           if (pub_hubbard_restart .and. pub_hubbard_atomsolve) then
               if (pub_on_root) then
                  write(stdout,'(/a)') 'Error in &
                       &internal_cdft: Cannot have negative values '
                  write(stdout,'(a)') 'for DFT+U charge Z with PAW or &
                       &negative projector mixing parameter in &
                       &<hubbard> block'
                end if
                call comms_abort
            endif
           ! 4 GIBO: TO BE CHECKED AFTER PAW GO AHEAD FROM NICK ==== END

            call comms_bcast(pub_root_node_id, pub_hubbard_atomsolve)
            call comms_bcast(pub_root_node_id, h_species(row)%hub_species)
            call comms_bcast(pub_root_node_id, h_species(row)%hub_ang_mom)
            call comms_bcast(pub_root_node_id, h_species(row)%hub_u)
            call comms_bcast(pub_root_node_id, h_species(row)%hub_charge)

            ! different Uq for alpha and beta spin 04.10.11
            call comms_bcast(pub_root_node_id, h_species(row)%cdft_u_charge_up)
            call comms_bcast(pub_root_node_id, h_species(row)%cdft_u_charge_down)
            ! different Uq for alpha and beta spin 04.10.11

            call comms_bcast(pub_root_node_id, h_species(row)%cdft_u_spin)


            call comms_bcast(pub_root_node_id, h_species(row)%cdft_target_up)
            call comms_bcast(pub_root_node_id, h_species(row)%cdft_target_down)

            call comms_bcast(pub_root_node_id, h_species(row)%cdft_target_spin)

            call comms_bcast(pub_root_node_id, h_species(row)%cdft_donor)
            call comms_bcast(pub_root_node_id, h_species(row)%cdft_acceptor)


            ! gibo: to avoid problems, initialise (unused) Hubbard variables 
            !      and send them around
            h_species(row)%hub_alpha = 0._DP
            h_species(row)%hub_spin_splitting = 0._DP
            call comms_bcast(pub_root_node_id, h_species(row)%hub_alpha)
            call comms_bcast(pub_root_node_id, &
                 & h_species(row)%hub_spin_splitting)

         end do

         ! gibo: if <constrained_dft> block is present
         ![we are inside the CHECK_CDFT IF-THEN construct]
         ! make sure at least one cDFT-mode is activated
         if ( .not.(pub_cdft_atom_charge .OR. &
                pub_cdft_atom_spin  .OR. pub_cdft_group_charge      .OR. &
                pub_cdft_group_spin .OR. pub_cdft_group_charge_diff .OR. &
                pub_cdft_group_spin_diff ) ) then
           if (pub_on_root) then
              write(*,'(a)') 'Error in internal_cdft: no cDFT mode selected, &
                &despite the presence of  <constrained_dft> block, aborting'
              write(*,'(a)') 'Select *ONE* from available cDFT-modes:'
              write(*,'(a)') '[1] pub_cdft_atom_charge       : T'
              write(*,'(a)') '[2] pub_cdft_atom_spin         : T'
              write(*,'(a)') '[3] pub_cdft_group_charge      : T'
              write(*,'(a)') '[4] pub_cdft_group_spin        : T'
              write(*,'(a)') '[5] pub_cdft_group_charge_diff : T'
              write(*,'(a)') '[6] pub_cdft_group_spin_diff   : T'
           endif
           call comms_abort
         endif

         !gibo: for group-charge_diff runs, make sure all the acceptor and donor 
         !      atoms have the same constraining potentials
         if (pub_cdft_group_charge_diff) then

          do ii=1,hub_nsp

           do row = ii+1, hub_nsp

             if (h_species(ii)%cdft_donor .AND. h_species(row)%cdft_donor ) then

               if ( (h_species(ii)%cdft_u_charge_up   .NE. &
                     h_species(row)%cdft_u_charge_up) .OR. &
                    (h_species(ii)%cdft_u_charge_down .NE. &
                     h_species(row)%cdft_u_charge_down) ) then

                if (pub_on_root) write(stdout,'(/a)') &
                     'Error in internal_cdft: &
                     &mismatching U_charge values in <constrained_dft> block &
                     &for pub_cdft_group_charge_diff simulation [DONOR], &
                     &aborting' 
                call comms_abort

               endif

             endif

             if (h_species(ii)%cdft_acceptor .AND. h_species(row)%cdft_acceptor ) then

               if ( (h_species(ii)%cdft_u_charge_up   .NE. &
                     h_species(row)%cdft_u_charge_up) .OR. &
                    (h_species(ii)%cdft_u_charge_down .NE. &
                     h_species(row)%cdft_u_charge_down) ) then

                if (pub_on_root) write(stdout,'(/a)') &
                     'Error in internal_cdft: &
                     &mismatching U_charge values in <constrained_dft> block &
                     &for pub_cdft_group_charge_diff simulation [ACCEPTOR], &
                     &aborting'
                call comms_abort

               endif

             endif



           enddo
          enddo

          ! gibo: since all acceptor and donor atoms have the same constraining 
          !       potentials, set pub_cdft_charge_group_diff_u 
          !       [from 1st hub_species], and broadcast as needed in 
          !       sbrtne cdft_energy_total [hubbard_build_mod.F90]
          ! gibo: use ABS and act accordingly in 'sbrtne cdft_energy_total'
          ! [alternative: re-loop over hub_nsp and single out acc. and don. species]
          pub_cdft_group_diff_u = ABS( h_species(1)%cdft_u_charge_up )
          call comms_bcast(pub_root_node_id,pub_cdft_group_diff_u)

         endif

         ! gibo: for group-spin_diff runs,
         !      make sure that all the acceptor and donor atoms have the same 
         !      constraining potentials
         if (pub_cdft_group_spin_diff) then

          do ii=1,hub_nsp

           do row = ii+1, hub_nsp

             if (h_species(ii)%cdft_donor .AND. h_species(row)%cdft_donor ) then

               if (h_species(ii)%cdft_u_spin .NE. h_species(row)%cdft_u_spin) then

                if (pub_on_root) write(stdout,'(/a)') &
                     'Error in internal_cdft: &
                     &mismatching U_spin values in <constrained_dft> block for &
                     &pub_cdft_group_spin_diff simulation [DONOR]'
                call comms_abort

               endif

             endif

             if (h_species(ii)%cdft_acceptor .AND. h_species(row)%cdft_acceptor ) then

               if (h_species(ii)%cdft_u_spin .NE. h_species(row)%cdft_u_spin) then

                if (pub_on_root) write(stdout,'(/a)') &
                     'Error in internal_cdft: &
                     &mismatching U_spin values in <constrained_dft> block for &
                     &pub_cdft_group_spin_diff simulation [ACCEPTOR]'
                call comms_abort

               endif

             endif

           enddo
          enddo

          ! gibo: since all acceptor and donor atoms have the same constraining 
          !       potentials, set pub_cdft_charge_group_diff_u 
          !       [from 1st hub_species], and broadcast as needed in 
          !       sbrtne cdft_energy_total [hubbard_build_mod.F90]
          ! gibo: use ABS and act accordingly in 'sbrtne cdft_energy_total'
          ! [alternative: re-loop over hub_nsp and single out acc. and don. species]
          pub_cdft_group_diff_u = ABS( h_species(1)%cdft_u_spin )
          call comms_bcast(pub_root_node_id,pub_cdft_group_diff_u)

         endif



         ! gibo: for group-charge runs,
         ! make sure all the group-atoms have the same constraining potentials
         if (pub_cdft_group_charge) then
          do ii=1,hub_nsp

           do row = ii+1, hub_nsp

               if ( (h_species(ii)%cdft_u_charge_up   .NE. &
                     h_species(row)%cdft_u_charge_up) .OR. &
                    (h_species(ii)%cdft_u_charge_down .NE. &
                     h_species(row)%cdft_u_charge_down) ) then

                if (pub_on_root) write(stdout,'(/a)') &
                     'Error in internal_cdft: &
                     &mismatching U_charge values in <constrained_dft> block &
                     &for pub_cdft_group_charge simulation, aborting'
                call comms_abort

               endif

           enddo
          enddo

          ! gibo: since all the group-atoms have the same constraining potentials,
          !       set pub_cdft_charge_group_u (from 1st hubbard_atom, we need it 
          !       in cdft_energy_total) and broadcast
          pub_cdft_group_u = h_species(1)%cdft_u_charge_up 
          call comms_bcast(pub_root_node_id,pub_cdft_group_u)

         endif

         ! gibo: for group-spin runs,
         !       make sure all the group-atoms have the same constraining
         !       potentials
         if (pub_cdft_group_spin) then
          do ii=1,hub_nsp

           do row = ii+1, hub_nsp

               if (h_species(ii)%cdft_u_spin .NE. h_species(row)%cdft_u_spin) then

                if (pub_on_root) write(stdout,'(/a)') &
                     'Error in internal_cdft: &
                     &mismatching U_spin values in <constrained_dft> block for &
                     &pub_cdft_group_spin simulation, aborting'
                call comms_abort

               endif

           enddo
          enddo

          ! gibo: since all the group-atoms have the same constraining potentials,
          !       set pub_cdft_charge_group_u (from 1st hubbard_atom, we need it 
          !       in cdft_energy_total) and broadcast
          pub_cdft_group_u = h_species(1)%cdft_u_spin 

          call comms_bcast(pub_root_node_id,pub_cdft_group_u)

         endif

         if (.not. element_found) then
            if (pub_on_root) write(stdout,'(/a)') &
                 'Error in internal_cdft: &
                 &mismatching species in <constrained_dft> and <positions_abs> &
                 &blocks of input file, aborting'
            call comms_abort
         else
            pub_cell%nat_hub = dummy_h ! The number of Hubbard atoms in the cell
            pub_cell%num_hub_proj = num_hub_proj
         end if

         ! ddor: broadcast number of Hubbard atoms in the cell
         call comms_bcast(pub_root_node_id,pub_cell%nat_hub)
         call comms_bcast(pub_root_node_id,pub_cell%num_hub_proj)

#ifdef DEBUG
         if (pub_on_root) write(stdout,'(a)') 'DEBUG: Read <constrained_dft> block'
#endif
      end if CHECK_CDFT

      ! pdh: deallocate buffers
      deallocate(cbuf,stat=ierr)
      call utils_dealloc_check('internal_cdft','cbuf',ierr)
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('internal_cdft','dbuf',ierr)
      deallocate(ibuf,stat=ierr)
      call utils_dealloc_check('internal_cdft','ibuf',ierr)

    end subroutine internal_cdft



    !-------------------------------
    !--------- %block classical_info
    !-------------------------------

    subroutine internal_classical_atoms()

      use classical_pot, only: classical_elements,  classical_pot_init
      use comms, only: pub_root_node_id, comms_abort, comms_bcast
      use esdf, only: block_data, esdf_block

      implicit none

      !-------------
      ! ars: buffers
      integer, allocatable :: ibuf_class(:)
      real(kind=DP), allocatable :: dbuf_class(:)
      ! ars: integer
      integer :: row, ii
      ! ars: character
      character(len=4),allocatable  :: species_id_classical(:)
      ! ars: real
      real(kind=DP), allocatable    :: classical_charge(:)
      ! ars: logical
      logical :: class_nats_exist


      ! cks: initialise number of classical atoms
      class_nat =0
      if (pub_on_root) then
         class_nats_exist =esdf_block('classical_info',class_nat)
      endif

      ! cks: make sure all cores have info about classical nats
      call comms_bcast(pub_root_node_id,class_nat)
      call comms_bcast(pub_root_node_id,class_nats_exist)

      ! cks: allocate classical arrays only if there are classical atoms
      if (class_nats_exist) then

         call classical_pot_init(class_nat)

         allocate(species_id_classical(class_nat),stat=ierr)
         call utils_alloc_check('internal_classical_atoms',&
              'species_id_classical',ierr)
         allocate(classical_charge(class_nat),stat=ierr)
         call utils_alloc_check('internal_classical_atoms','classical_charge',&
              ierr)

         allocate(dbuf_class(3*class_nat),stat=ierr)
         call utils_alloc_check('internal_classical_atoms','dbuff_class',&
              ierr)

         allocate(ibuf_class(4*class_nat),stat=ierr)
         call utils_alloc_check('internal_classical_atoms','ibuf_class',&
              ierr)

         ! ars: Root node only: read absolute cartesian positions from
         !      <classical_info> block of input file
         if (pub_on_root) then

            do row=1,class_nat
               read(block_data(row),*) species_id_classical(row),&
                    dbuf_class(row*3-2:row*3),classical_charge(row)
               do ii=1,4
                  ibuf_class((row-1)*4+ii) = iachar(species_id_classical(row)(ii:ii))
               end do

            end do
         end if


         ! cks: broadcast positions and charges
         call comms_bcast(pub_root_node_id,dbuf_class,class_nat*3)
         call comms_bcast(pub_root_node_id,classical_charge,class_nat)

         ! cks: broadcast and unpack species_id_classical
         call comms_bcast(pub_root_node_id,ibuf_class)
         do row=1,class_nat
            do ii=1,4
               species_id_classical(row)(ii:ii) = achar(ibuf_class((row-1)*4+ii))
            end do
         end do

         ! cks: unpack and store ids, coordinates and charges
         do row=1,class_nat
            classical_elements(row)%symbol = species_id_classical(row)(1:2)
            classical_elements(row)%centre%x = dbuf_class(row*3-2)
            classical_elements(row)%centre%y = dbuf_class(row*3-1)
            classical_elements(row)%centre%z = dbuf_class(row*3)
            classical_elements(row)%ion_charge = classical_charge(row)
         end do

         ! cks: deallocate classical arrays
         deallocate(ibuf_class,stat=ierr)
         call utils_dealloc_check('internal_classical_atoms',&
              'ibuf_class',ierr)

         deallocate(dbuf_class, stat=ierr)
         call utils_dealloc_check('internal_classical_atoms',&
              'dbuff_class',ierr)

         deallocate(classical_charge,stat=ierr)
         call utils_dealloc_check('internal_classical_atoms',&
              'classical_charge',ierr)

         deallocate(species_id_classical,stat=ierr)
         call utils_dealloc_check('internal_classical_atoms',&
              'species_id_classical',ierr)
      endif

    end subroutine internal_classical_atoms


 !--------------------------------------------------------------------------!
 ! lpl: NBO Partial matrix atom list & NGWF labelling (16/06/2011)          !
 !      Modified from %block_species_ngwf_plot in this module               !
 !      Reads: %block nbo_write_species                                     !
 !             %block nbo_species_ngwflabel                                 !
 !             %block nbo_list_plotnbo                                      !
 !--------------------------------------------------------------------------!

    subroutine internal_species_write_nbo()

      use comms,  only: pub_root_node_id, comms_abort, comms_bcast
      use esdf,   only: block_data, esdf_block
      use rundat, only: pub_nbo_write_species, pub_nbo_ngwf_label, &
                        pub_nbo_list_plotnbo, pub_plot_nbo

      implicit none

      ! lpl: Local variables
      character(len=4) :: dummy_id
      character(len=256) :: dummy_ngwf_label
      integer :: dummy_nsp_cbuf1   ! lpl: For <nbo_write_species>
      integer :: dummy_nsp_cbuf2   ! lpl: For <nbo_species_ngwflabel>
      integer :: nnbo, nbo_it

      ! lpl: Local list buffers
      ! lpl: For <nbo_write_species>
      character(len=68), allocatable :: cbuf1(:)
      ! lpl: For <nbo_species_ngwflabel>
      character(len=256),  allocatable :: cbuf2(:)

      ! lpl: Misc. scalars
      integer ::  row, ii, ierr

      ! lpl: Checking flags
      logical :: element_found
      logical :: skip_species_list   ! lpl: For <nbo_write_species>
      logical :: nbolist_exists, skip

      ierr = 0

      allocate(cbuf1(nsp),stat=ierr)
      call utils_alloc_check('internal_species_write_nbo','cbuf1',ierr)
      allocate(cbuf2(nsp),stat=ierr)
      call utils_alloc_check('internal_species_write_nbo','cbuf2',ierr)

      ! lpl: Allocate necessary arrays
      if (allocated(pub_nbo_write_species)) then
         deallocate(pub_nbo_write_species,stat=ierr)
         call utils_dealloc_check('internal_species_write_nbo', &
              'pub_nbo_write_species',ierr)
      end if
      allocate(pub_nbo_write_species(nat),stat=ierr)
      call utils_alloc_check('internal_species_write_nbo', &
           'pub_nbo_write_species',ierr)

      if (allocated(pub_nbo_ngwf_label)) then
         deallocate(pub_nbo_ngwf_label,stat=ierr)
         call utils_dealloc_check('internal_species_write_nbo', &
              'pub_nbo_ngwf_label',ierr)
      end if
      allocate(pub_nbo_ngwf_label(nsp),stat=ierr)
      call utils_alloc_check('internal_species_write_nbo', &
           'pub_nbo_ngwf_label',ierr)

      if (allocated(pub_nbo_list_plotnbo)) then
         deallocate(pub_nbo_list_plotnbo,stat=ierr)
         call utils_dealloc_check('internal_species_write_nbo', &
              'pub_nbo_list_plotnbo',ierr)
      end if

      ! lpl: Get NBO list block if plot_nbo = T
      nnbo = 0
      if(pub_plot_nbo) then

         if(pub_on_root) nbolist_exists &
              = esdf_block('nbo_list_plotnbo',nnbo)

         call comms_bcast(pub_root_node_id,nbolist_exists)
         call comms_bcast(pub_root_node_id,nnbo)

         if(nbolist_exists .and. nnbo > 0) then
            allocate(pub_nbo_list_plotnbo(0:nnbo),stat=ierr)
            call utils_alloc_check('internal_species_write_nbo', &
                 'pub_nbo_list_plotnbo',ierr)

            pub_nbo_list_plotnbo(0) = 0 ! idx 0 stores #NBOs
            pub_nbo_list_plotnbo(1:nnbo) = -1 ! lpl: -1 for error checking
         else
            write(stdout,'(a)') 'ERROR: Error reading block &
                 &nbo_list_plotnbo'
            ierr=1
         end if

         if(pub_on_root) then ! lpl: Read block

            do row=1,nnbo
               read(block_data(row),*) nbo_it
               if(nbo_it < 0) then
                  write(stdout,'(a)') 'ERROR: nbo_it < 0'
                  ierr = 1
               end if

               ! lpl: Check for multiply defined NBOs
               if(row == 1) then ! Read 1st element
                  pub_nbo_list_plotnbo(0) = 1
                  pub_nbo_list_plotnbo(1) = nbo_it
               else ! Check against read elements

                  skip = .false.
                  do ii=1,pub_nbo_list_plotnbo(0)-1
                     if(nbo_it == pub_nbo_list_plotnbo(ii)) then
                        write(stdout,'(a,I5,a)') 'WARNING: NBO # ', &
                             nbo_it,' multiply defined.'
                        skip = .true.
                        exit
                     end if
                  end do ! END  do col=1,row-1
                  if(.not. skip) then
                     pub_nbo_list_plotnbo(0) = & ! idx 0 stores # NBOs
                          pub_nbo_list_plotnbo(0) + 1
                     pub_nbo_list_plotnbo(pub_nbo_list_plotnbo(0)) &
                          = nbo_it
                  end if

               end if
            end do ! END do row=1,nnbo

         end if ! END if(pub_on_root)

      end if ! END if(pub_plot_nbo)

      call comms_bcast(pub_root_node_id,ierr)
      if(ierr /= 0) call comms_abort

      ! lpl: Broadcast pub_nbo_list_plotnbo if present & allocated
      if (allocated(pub_nbo_list_plotnbo)) &
           call comms_bcast(pub_root_node_id,pub_nbo_list_plotnbo)

      ! lpl: initialise info arrays to default
      pub_nbo_write_species = .false.
      do row=1,nsp
         pub_nbo_ngwf_label(row)(1:4) = all_species_id(row)(1:4)
         pub_nbo_ngwf_label(row)(5:)  = 'AUTO'
      end do

      ! lpl: Read info from <nbo_write_species> and <nbo_species_ngwflabel>
      !      blocks in input file on root node only. If blocks are
      !      not present default values are assumed
      ierr = 0
      if (pub_on_root) then

         ! lpl: <nbo_write_species>
         skip_species_list = .false.
         if (esdf_block('nbo_write_species',dummy_nsp_cbuf1)) then
            if (dummy_nsp_cbuf1 > nsp) then ! aam: error check
               write(stdout,'(/a)') 'Error in internal_species_write_nbo: &
                    &too many species in <nbo_write_species>'
               write(stdout,'(a)')  'block of input file'
               ierr = 1
            else
               ! aam: check for multiply defined species
               call internal_check_species(dummy_nsp_cbuf1, &
                    '<nbo_write_species>')
               ! aam: dummy_nsp is not necessarily equal to nsp
               do row=1,dummy_nsp_cbuf1
                  read(block_data(row),*) dummy_id
                  cbuf1(row)(1:4) = dummy_id
               end do
            end if
         else
            ! lpl: If block isn't specified, skip reading species list
            skip_species_list = .true.
         end if

         ! lpl: <nbo_species_ngwflabel>
         if (esdf_block('nbo_species_ngwflabel',dummy_nsp_cbuf2)) then
            if (dummy_nsp_cbuf2 > nsp) then ! aam: error check
               write(stdout,'(/a)') 'Error in internal_species_write_nbo: &
                    &too many species in <nbo_species_ngwflabel>'
               write(stdout,'(a)')  'block of input file'
               ierr = 1
            else
               ! aam: check for multiply defined species
               call internal_check_species(dummy_nsp_cbuf2, &
                    '<nbo_species_ngwflabel>')
               ! aam: dummy_nsp is not necessarily equal to nsp
               do row=1,dummy_nsp_cbuf2 
                  read(block_data(row),*) dummy_id,dummy_ngwf_label
                  cbuf2(row)(1:4) = dummy_id
                  cbuf2(row)(5:)  = dummy_ngwf_label
               end do
            end if
         end if

      end if

      ! pdh: broadcast error flag for a clean shutdown if necessary
      call comms_bcast(pub_root_node_id,ierr)
      if (ierr /= 0) call comms_abort

      ! lpl: broadcast string buffers & indicators
      call comms_bcast(pub_root_node_id,dummy_nsp_cbuf1)
      call comms_bcast(pub_root_node_id,dummy_nsp_cbuf2)
      do row=1,dummy_nsp_cbuf1
         call comms_bcast(pub_root_node_id,cbuf1(row),4)
      end do
      do row=1,dummy_nsp_cbuf2
         call comms_bcast(pub_root_node_id,cbuf2(row))
      end do
      call comms_bcast(pub_root_node_id,skip_species_list)

      ! lpl: extract information from <nbo_write_species>
      if(.not. skip_species_list) then
         do row=1,dummy_nsp_cbuf1

            element_found = .false. ! aam: error check flag
            dummy_id = cbuf1(row)(1:4)

            ! lpl: <nbo_write_species> loop
            do ii=1,nat
               if (elements(ii)%species_id == dummy_id) then
                  pub_nbo_write_species(ii) = .true.
                  element_found = .true.
               end if
            end do

            if (.not. element_found) then
               if (pub_on_root) then
                  write(stdout,'(/a)') 'Error in internal_nbo_write_species: &
                       &mismatching species in'
                  write(stdout,'(a)') '<internal_species_write_nbo> and &
                       &<positions_abs> blocks of input file'
               end if
               call comms_abort
            end if
         end do
      else
      ! lpl: If no nbo_write_species block is specified, include all atoms
      !      in output
         pub_nbo_write_species = .true.
      end if

      ! lpl: extract information from <nbo_species_ngwflabel>
      do row=1,dummy_nsp_cbuf2

         element_found = .false. ! aam: error check flag
         dummy_id = cbuf2(row)(1:4)
         dummy_ngwf_label = cbuf2(row)(5:)

         ! lpl: Check that this species exists
         do ii=1,nat
            if (elements(ii)%species_id == dummy_id) element_found = .true.
         end do

         if (.not. element_found) then
            if (pub_on_root) then
               write(stdout,'(/a)') 'Error in internal_species_write_nbo: &
                    &mismatching species in'
               write(stdout,'(a)') '<nbo_species_ngwflabel> and &
                    &<positions_abs> blocks of input file'
            end if
            call comms_abort
         end if

         ! lpl: Overwrite default ('AUTO') with use-specified configuration
         do ii=1,nsp
            if(pub_nbo_ngwf_label(ii)(1:4) == dummy_id) then
               pub_nbo_ngwf_label(ii)(5:)  = dummy_ngwf_label
            end if
         end do

      end do

      ! lpl: Deallocate buffers
      deallocate(cbuf1,stat=ierr)
      call utils_dealloc_check('internal_species_write_nbo','cbuf1',ierr)
      deallocate(cbuf2,stat=ierr)
      call utils_dealloc_check('internal_species_write_nbo','cbuf2',ierr)

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') &
           'DEBUG: Read <nbo_write_species> block'
#endif

    end subroutine internal_species_write_nbo

 !--------------------------------------------------------------------------!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !--------------------------------------------------------------------------!


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!---------------------------   END DATA BLOCKS    --------------------------!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!





    !-----------------------
    !--------- check_species
    !-----------------------
    subroutine internal_check_species(nrecords,blockname)

      use comms, only: comms_abort
      use esdf, only: block_data

      implicit none

      integer, intent(in) :: nrecords
      character(len=*), intent(in) :: blockname

      ! <<<local variables>>>
      integer :: pp,qq

      do pp=1,nrecords
         read(block_data(pp),*) id_list(pp)
      enddo
      do pp=1,nrecords
         do qq=pp+1,nrecords
            if (id_list(pp).eq.id_list(qq)) then
               write(stdout,'(/a,a,a)') 'Error in rundat_blocks_read: &
                    &multiply defined species in ',trim(blockname),' block &
                    &of input file'
               call comms_abort
            endif
         enddo
      enddo

    end subroutine internal_check_species

    subroutine internal_vdwparam_override()

      !==================================================================!
      ! This subroutine reads the VDW parameter override block from the  !
      ! input file and dumps it into vdwparam_override in the            !
      ! vdwcorrection module.                                            !
      !------------------------------------------------------------------!
      ! Written by Quintin Hill in December 2010.                        !
      !==================================================================!

      use comms, only: pub_on_root
      use esdf, only: block_data, esdf_block
      use vdwcorrection, only: vdwcorrection_override_alloc

      implicit none

      integer :: num_rows ! number of rows

      if (pub_on_root) then
         if (esdf_block('vdw_params',num_rows)) then
            call vdwcorrection_override_alloc(block_data(1:num_rows),num_rows)
         end if
      end if

    end subroutine internal_vdwparam_override

#ifdef DEBUG

    subroutine internal_debug()

      use comms, only: comms_barrier
      use hubbard_init, only: h_species
      use rundat, only: pub_hubbard, pub_cond_calculate
      use simulation_cell, only: pub_cell

      implicit none

      integer :: ii, row

      call comms_barrier

      if (pub_on_root) then
         write(stdout,'(/a,i6)') 'DEBUG: No of Atoms:',nat
         write(stdout,'(a,i5)') 'DEBUG: No of Species: ',nsp
         if (pub_hubbard) then
            write(stdout,'(/a,i6)') 'DEBUG: No of Hubbard Atoms:',pub_cell%nat_hub
            write(stdout,'(a,i5)') 'DEBUG: No of Hubbard Species: ',hub_nsp
         endif
         write(stdout,*) 'DEBUG: =======================================&
              &============================'
         do ii=1,nat
            write(stdout,*) 'DEBUG: Atom       :',ii
            write(stdout,*) 'DEBUG: Symbol     : ',elements(ii)%symbol
            write(stdout,*) 'DEBUG: Z          :',elements(ii)%atomic_number
            write(stdout,*) 'DEBUG: NGWFs      :',elements(ii)%nfunctions
            write(stdout,*) 'DEBUG: Radius     :',elements(ii)%radius
            write(stdout,*) 'DEBUG: Pseudopot  : ',trim(elements(ii)%pseudo_name)
            write(stdout,*) 'DEBUG: Fireball   : ',trim(elements(ii)%ngwf_set)
            write(stdout,'(a,3f10.5)') ' DEBUG: Centre     :', &
                 elements(ii)%centre%x,&
                 elements(ii)%centre%y,elements(ii)%centre%z
            write(stdout,*) 'DEBUG: Constraint : ', &
                 elements(ii)%ion_constraint_type
            write(stdout,'(a,3f10.5)') ' DEBUG: Constraint :', &
                 elements(ii)%ion_constraint(1),&
                 elements(ii)%ion_constraint(2),elements(ii)%ion_constraint(3)
            if (pub_hubbard) then
               do row=1,hub_nsp
                  if (elements(ii)%species_id == h_species(row)%hub_species) then
                     write(stdout,*) 'DEBUG: Hubbard projector angular momentum :',&
                          h_species(row)%hub_ang_mom
                     write(stdout,*) 'DEBUG: Hubbard U parameter (Ha)           :',&
                          h_species(row)%hub_u
                     write(stdout,*) 'DEBUG: Hubbard projector effective charge :',&
                          h_species(row)%hub_charge
                     write(stdout,*) 'DEBUG: Hubbard alpha parameter (Ha)       :',&
                          h_species(row)%hub_alpha
                     write(stdout,*) 'DEBUG: Hubbard spin splitting (Ha)        :',&
                          h_species(row)%hub_spin_splitting
                  endif
               enddo
            endif
            if (pub_cond_calculate) then
               write(stdout,*) 'DEBUG: Conduction NGWFs   :',elements(ii)%nfunctions_cond
               write(stdout,*) 'DEBUG: Conduction radius  :',elements(ii)%radius_cond
            end if
            write(stdout,*) 'DEBUG: ====================================&
                 &==============================='
         end do
      end if

    end subroutine internal_debug

#endif

  end subroutine rundat_blocks_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine thermostat_read_params(line,temp,tkind,tstart,tstop)
 
    use constants, only: dp,stdout
    use comms,     only: comms_abort, pub_on_root 
    use esdf,      only: esdf_line_divide, esdf_convfac, esdf_reduce
  
    implicit none
  
    ! Arguments
    character(len=80), intent(in) :: line
    real(kind=dp), intent(out)  :: temp
    integer, intent(out)        :: tkind, tstart, tstop
    
    ! Local variables
    character(len=80) :: cbuf(5), cjunk
    character(len=20) :: runit
    real(kind=dp)     :: confac
    integer           :: maxp

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') &
           'DEBUG: Enter thermostat_read_params'
#endif

    temp = 0.0_dp
    tkind = 0 
    tstart = 0
    tstop = 0
 
    maxp = 5 
    call esdf_line_divide(maxp,cbuf,line)

    if (maxp .ne. 5) then 
       write(stdout,'(/a)') 'Error in thermostat_read_params: &
             &mandatory parameters not found in input file'
       call comms_abort
    endif

    ! Read thermostat range
    read(cbuf(1),*) tstart
    read(cbuf(2),*) tstop
  
    ! Process thermostat type
    cjunk = esdf_reduce(cbuf(3))
    select case (cjunk)
    case ('none') 
       tkind = 0
    case ('andersen') 
       tkind = 1
    case ('langevin') 
       tkind = 2
    case ('nosehoover') 
       tkind = 3
    end select
  
    ! Read thermostat temperature
    read(cbuf(4),*) temp
    read(cbuf(5),*) runit
    confac = esdf_convfac(runit,'hartree')
    temp = temp*confac
  

  end subroutine thermostat_read_params

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine thermostat_read_options(line,tgrad,freq,mix,damp,group,nhc,nhi,upd)
  
    use constants, only: dp,stdout
    use comms, only: pub_on_root 
    use esdf, only: esdf_line_divide, esdf_convfac, esdf_reduce
    use rundat, only: md_delta_t
  
    implicit none
  
    ! Arguments 
    character(len=80), intent(in) :: line
    real(kind=dp), intent(inout) :: tgrad, freq, mix, damp
    integer, intent(inout) :: group, nhc, nhi
    logical, intent(inout) :: upd
    
    ! Local variables
    character(len=80) :: cbuf(16), cjunk
    character(len=20) :: runit
    real(kind=dp)     :: confac
    integer           :: ip, maxp

    maxp = 16
    call esdf_line_divide(maxp,cbuf,line)

    if (pub_on_root) then 
       write(stdout,'(a,i4)') 'maxp ', maxp
       do ip = 1, maxp
          write(stdout,'(i4,3x,a)') ip,cbuf(ip)
       enddo
    endif
    ip = 1
    do while (ip.lt.maxp)
       cjunk = esdf_reduce(cbuf(ip))
       select case (cjunk)
       case ('tgrad')
          read(cbuf(ip+1),*) tgrad
          read(cbuf(ip+2),*) runit
          confac = esdf_convfac(runit,'hartree')
          tgrad = tgrad*confac/md_delta_t 
          ip = ip + 3
       case ('freq') 
          read(cbuf(ip+1),*) freq
          freq = freq/md_delta_t
          ip = ip + 2
       case ('mix')
          read(cbuf(ip+1),*) mix
          ip = ip + 2
       case ('damp')
          read(cbuf(ip+1),*) damp
          ip = ip + 2
       case ('group')
          read(cbuf(ip+1),*) group
          ip = ip + 2
       case ('nchain')
          read(cbuf(ip+1),*) nhc
          ip = ip + 2
       case ('nstep')
          read(cbuf(ip+1),*) nhi
          ip = ip + 2
       case ('update')
          read(cbuf(ip+1),*) upd
          ip = ip + 2
       end select
    enddo
  
  end subroutine thermostat_read_options

  !------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine rundat_blocks_init_grid(a1, a2, a3, a1_pad, a2_pad, a3_pad, &
       elements, num, num_cond, num_aux, nat, class_nat)

    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use geometry, only: point
    use ion, only: element
    use ppd_strategy, only: ppd_strategy_determine
    use rundat, only: pub_coulomb_cutoff, psinc_spacing
    use simulation_cell, only: pub_cell, pub_padded_cell, &
         simulation_cell_initialise, simulation_cell_add_padding
    use services, only: services_rationalise_coords

    implicit none

    ! Arguments
    type(POINT),   intent(in   ) :: a1, a2, a3
    type(POINT),   intent(inout) :: a1_pad, a2_pad, a3_pad
    integer,       intent(in   ) :: nat, num, num_cond, num_aux, class_nat
    type(ELEMENT), intent(inout) :: elements(nat)

    ! Local Variables
    real(kind=DP) :: d1, d2, d3
    real(kind=DP) :: r_abs(3,nat)
    integer       :: iat
    integer       :: n_pt1, n_pt2, n_pt3


    ! cks: determine grids and ppds based on KE-cutoff input parameter
    call ppd_strategy_determine(d1,d2,d3,n_pt1,n_pt2,n_pt3, &
         a1,a2,a3,elements)

    ! Initialise pub_cell
    call simulation_cell_initialise(pub_cell,elements,num,num_cond,num_aux, &
         nat,class_nat,a1,a2,a3,n_pt1,n_pt2,n_pt3,d1,d2,d3)

    ! ndmh: initialise pub_padded_cell
    if (pub_coulomb_cutoff) then

       ! ndmh: apply override to psinc spacing so that spacing in padded
       ! ndmh: cell matches that of original cell
       psinc_spacing(1) = d1
       psinc_spacing(2) = d2
       psinc_spacing(3) = d3

       ! ndmh: check padding
       call simulation_cell_add_padding(d1,d2,d3,n_pt1,n_pt2,n_pt3, &
            a1,a2,a3,a1_pad,a2_pad,a3_pad)

       ! ndmh: determine grids and ppds of pub_padded_cell
       call ppd_strategy_determine(d1, d2, d3, n_pt1, n_pt2, n_pt3, &
            a1_pad, a2_pad, a3_pad, elements)

       ! ndmh: initialise pub_padded_cell
       call simulation_cell_initialise(pub_padded_cell, elements, num, num_cond, &
            num_aux, nat, class_nat, a1_pad, a2_pad, a3_pad, n_pt1, n_pt2, n_pt3, &
            d1, d2, d3)

    end if

    ! smmd: ensure all ions are contained within the cell
    do iat=1,nat
       r_abs(1,iat) = elements(iat)%centre%x
       r_abs(2,iat) = elements(iat)%centre%y
       r_abs(3,iat) = elements(iat)%centre%z
    enddo

    call services_rationalise_coords(nat,r_abs)

    do iat=1,nat
       elements(iat)%centre%x =  r_abs(1,iat)
       elements(iat)%centre%y =  r_abs(2,iat)
       elements(iat)%centre%z =  r_abs(3,iat)
    enddo

#ifdef DEBUG
         if (pub_on_root) write(stdout,'(a)') 'DEBUG: exit init_blocks_init_grid'
#endif

  end subroutine rundat_blocks_init_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rundat_blocks_old(elements,nat)

    !=============================================================!
    ! This subroutine reads information about the cell from the   !
    ! input file using the esdf module.                           !
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2001.                   !
    ! Modified by Peter D. Haynes in 2004 so that the root node   !
    ! reads and communicates across all nodes.                    !
    ! Modified by Chris-Kriton Skylaris on 09/11/2004 so that it  !
    ! determines grids by calling ppd_strategy rather than        !
    ! reading them from input.                                    !
    !=============================================================!
    ! Adapted and cleaned up from services_read_system_input_old  !
    ! by Alvaro Ruiz Serrano, 17/11/2010.                         !
    !=============================================================!

    use comms, only: comms_abort, comms_bcast, pub_on_root, pub_root_node_id
    use constants, only: DP, stdout
    use esdf, only: block_data, esdf_block
    use geometry, only: point
    use ion, only: element
    use ppd_strategy, only: ppd_strategy_determine
    use rundat, only: pub_ngwf_halo
    use simulation_cell, only: pub_cell, simulation_cell_initialise

    implicit none

    integer, intent(in) :: nat
    type(element), intent(out) :: elements(nat)

    ! cks: internal
    integer :: row,dummy_nat,num_lines,num  !, ppd_1d_npoints
    integer :: n_pt1, n_pt2, n_pt3
    type(POINT) :: a1,a2,a3
    real(kind=DP) :: d1, d2, d3

    if (pub_on_root) then
       dummy_nat=-1
       if (esdf_block('geometry',dummy_nat)) then
          if (nat.eq.dummy_nat)  then

             num=0
             do row=1,nat

                read(block_data(row),*) elements(row)%symbol,&
                     elements(row)%atomic_number,elements(row)%pseudo_name, &
                     elements(row)%centre%x,elements(row)%centre%y,&
                     elements(row)%centre%z,elements(row)%nfunctions,&
                     elements(row)%radius

                ! cks: 21/04/2006 Warning! If there is a halo, the element radius is
                ! cks: the NGWF radius plus the halo length and so it is not equal to
                ! cks: the NGWF-sphere radius (which does not include the halo)
                if (pub_ngwf_halo > 0.0_DP) then
                   ! cks: increase element radius by halo length
                   elements(row)%radius =elements(row)%radius +pub_ngwf_halo
                endif

                ! cks: this counts the total number of NGWFs by
                !      adding the number of functions on every atom.
                num=num+elements(row)%nfunctions

             enddo
          else
             write(stdout,'(a)') 'Error in services_read_system_input: &
                  &mismatching numbers of atoms in geometry specification'
             call comms_abort
          end if
       else
          write(stdout,'(a)') 'Error in services_read_system_input: &
               &no geometry specification'
          call comms_abort
       end if
    end if

    do row=1,nat
       call comms_bcast(pub_root_node_id,elements(row)%symbol)
       call comms_bcast(pub_root_node_id,elements(row)%atomic_number)
       call comms_bcast(pub_root_node_id,elements(row)%pseudo_name)
       call comms_bcast(pub_root_node_id,elements(row)%centre%x)
       call comms_bcast(pub_root_node_id,elements(row)%centre%y)
       call comms_bcast(pub_root_node_id,elements(row)%centre%z)
       call comms_bcast(pub_root_node_id,elements(row)%nfunctions)
       call comms_bcast(pub_root_node_id,elements(row)%radius)
    end do
    call comms_bcast(pub_root_node_id,num)

    if (pub_on_root) then
       dummy_nat=-1
       if (esdf_block('atomic_sets',dummy_nat)) then
          if (nat.eq.dummy_nat)  then
             do row=1,nat
                read(block_data(row),*) elements(row)%ngwf_set
             enddo
          else
             write(stdout,'(a)') 'Error in services_read_system_input: &
                  &mismatching numbers of atoms in atomic sets specification'
             call comms_abort
          end if
       else
          write(stdout,'(a)') 'Error in services_read_system_input: &
               &no atomic sets specification'
          call comms_abort
       end if
    end if

    do row=1,nat
       call comms_bcast(pub_root_node_id,elements(row)%ngwf_set,64)
    end do

    if (pub_on_root) then
       num_lines=-1
       if (esdf_block('lattice_cart',num_lines)) then
          if (num_lines /= 3) then
             write(stdout,'(a)') 'Error in services_read_system_input: &
                  &malformed lattice vector specification'
             call comms_abort
          end if

          read(block_data(1),*) a1%x,a1%y,a1%z
          read(block_data(2),*) a2%x,a2%y,a2%z
          read(block_data(3),*) a3%x,a3%y,a3%z

       end if
    end if

    call comms_bcast(pub_root_node_id,a1%x)
    call comms_bcast(pub_root_node_id,a1%y)
    call comms_bcast(pub_root_node_id,a1%z)
    call comms_bcast(pub_root_node_id,a2%x)
    call comms_bcast(pub_root_node_id,a2%y)
    call comms_bcast(pub_root_node_id,a2%z)
    call comms_bcast(pub_root_node_id,a3%x)
    call comms_bcast(pub_root_node_id,a3%y)
    call comms_bcast(pub_root_node_id,a3%z)

#ifdef DEBUG
    ! DIAGNOSTICS
    if (pub_on_root) then
       write(stdout,*) 'No of Atoms:   ',nat
       write(stdout,*) 'No of Species:  N/A'
       write(stdout,*) '========================================================================='
       do row=1,nat
          write(stdout,*) 'Atom      :',row
          write(stdout,*) 'ID        : N/A'
          write(stdout,*) 'Symbol    :',elements(row)%symbol
          write(stdout,*) 'Z         :',elements(row)%atomic_number
          write(stdout,*) 'NGWFs     :',elements(row)%nfunctions
          write(stdout,*) 'Radius    :',elements(row)%radius
          write(stdout,*) 'Pseudopot :',elements(row)%pseudo_name
          write(stdout,*) 'Fireball  :',elements(row)%ngwf_set
          write(stdout,*) 'Centre    :',elements(row)%centre%x,&
               elements(row)%centre%y,elements(row)%centre%z
          write(stdout,*) '========================================================================='
       enddo
    endif
#endif

    ! cks: determine grids and ppds based on KE-cutoff input parameter
    call ppd_strategy_determine(d1, d2, d3, n_pt1, n_pt2, n_pt3, &
         a1, a2, a3, elements)

    ! Initialise pub_cell
    call simulation_cell_initialise(pub_cell,elements,num,0,nat,0,0,&
         a1,a2,a3,n_pt1,n_pt2,n_pt3,d1,d2,d3)


  end subroutine rundat_blocks_old

end module rundat_blocks
