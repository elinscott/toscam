! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                    T S S E A R C H                                          !
!=============================================================================!
!                                                                             !
! $Id: tssearch_mod.F90,v 1.11 2009/09/23 16:56:58 cks22 Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! This module performs a transition state search starting from a reactant,    !
! product and intermediate (if required) structures                           !
!                                                                             !
!-----------------------------------------------------------------------------!
! Written from  "Module Specification for New Plane Wave Code",  M. Segall,   !
!    P. Lindan, M. Probert, C. Pickard, P. Hasnip, S. Clark and M. Payne      !
!                           Copyright 1999/2000                               !
!-----------------------------------------------------------------------------!
! Written by Niri Govind, v0.1, 19/11/2001                                    !
!-----------------------------------------------------------------------------!
! modification information                                                    !
!=============================================================================!

module tssearch
  use simulation_cell, only: castep_cell

  implicit none                                 !Impose strong typing

  private                                       !Everything is private ...

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: tssearch_run
  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!


  integer, save :: num_ionic_constraints

  ! Cell data in CASTEP format
  type(castep_cell)  :: current_cell

contains

  subroutine tssearch_run(elements,output_file)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mdl, intent=inout, the model to be optimised                          !
    !   output_file, intent=in, the filename to which results will be written.!
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   mdl has already been properly read and initialised                    !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, 20/11/2001                                !
    !=========================================================================!

    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout
    use ion, only   : element
    use rundat, only: tssearch_method
    use simulation_cell, only: pub_cell, castep_model, copy_to_castep_cell, &
         castep_model_dealloc, castep_cell_dealloc, castep_model_nullify, &
         castep_model_alloc, castep_cell_copy, castep_cell_2_elements
    use tssearch_utils, only: tssearch_utils_io_abort, constrain_ions, &
         tssearch_utils_initialize

    implicit none

    ! Arguments
    type(ELEMENT), intent(inout) :: elements(pub_cell%nat)
    character(len=*),  intent(in)    :: output_file

    ! <<< local variables >>>
    type(castep_model) :: mdl

    integer :: ndim  !the number of dimensions in the search space
                   !so for the augmented-Hessian BFGS method used here,
                   !we have ndim=9+3*num_ions as we seek to optimize the
                   !cell (9 vars) and the ions (3*num_ions d.of.f) simultaneously.

    ! qoh: Actually initialise ndim!
    ndim = 9 + 3*pub_cell%nat

    ! fill CASTEP-style cell data structures from pub_cell and "elements"
    call copy_to_castep_cell(current_cell,elements)

    ! Consistency check
    if (current_cell%num_ions<=0) then
       if (pub_on_root) write(stdout,'(a)') &
            'Error in ONETEP: current_cell uninitialised in tssearch_run'
       call comms_abort
    endif

        ! Nullify pointers in mdl
    call castep_model_nullify(mdl)

    ! Allocate model data
    call castep_model_alloc(mdl,ndim,constrain_ions,.false.)

    call castep_cell_copy(current_cell,mdl%cell)

    !set the mdl_ptr and output_file_ptr
    call tssearch_utils_initialize(mdl,elements) 

    !Analyse constraints to see what is being asked to relax
    constrain_ions=.false.

    num_ionic_constraints=min(3*mdl%cell%num_ions,num_ionic_constraints)  !catch any sillies
    num_ionic_constraints=max(0,num_ionic_constraints)                    !catch any sillies
    if (num_ionic_constraints>0) constrain_ions=.true.
    if (num_ionic_constraints==3*mdl%cell%num_ions) call tssearch_utils_io_abort &
       & ('ERROR: Cannot perform transition state search run when all atoms are fixed')

    !Write main header
    if (pub_on_root) then
      write(stdout,1)
      write(stdout,2) 'Starting ONETEP Transition State Search'
      write(stdout,1)
      write(stdout,'(/a,a,a)') &
            ' Transition State Search: output file is "',trim(output_file),'"'
    end if


    !Now select transition state search method
    select case (tssearch_method)
    case ('LSTQST')
       call tssearch_LSTQST(mdl,output_file)
    case default
       call tssearch_utils_io_abort('Error in tssearch_run - unrecognised transition state search request='//tssearch_method//'.')
    end select

    ! Make sure that the elements ONETEP array is updated
    call castep_cell_2_elements(mdl%cell,elements)

    ! Deallocate memory
    call castep_model_dealloc(mdl,.false.)
    call castep_cell_dealloc(current_cell)

1   format(80('='))
2   format(19('<'),1x,a,1x,20('>'))

    return
  end subroutine tssearch_run

  subroutine tssearch_LSTQST(mdl,output_file)
    !=========================================================================!
    ! Initializes the reactant, product and intermediate cells                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   I/O                                                                   !
    !   Wavefunction                                                          !
    !   Density                                                               !
    !   Cell                                                                  !
    !   Parameters                                                            !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, version 0.1, 21/11/01                           !
    !=========================================================================!

    use comms,    only : pub_on_root,comms_bcast,pub_root_node_id
    use constants,only : DP
    use esdf,     only : esdf_init, esdf_close
    use ion,      only : element
    use rundat,   only : pub_rootname, tssearch_lstqst_protocol, &
         tssearch_cg_max_iter, tssearch_force_tol, &
         tssearch_disp_tol, tssearch_qst_max_iter
    use rundat_blocks, only: rundat_blocks_exec
    use simulation_cell, only: castep_model, castep_cell_dealloc, &
         castep_cell_alloc, castep_cell_copy, castep_cell_frac_to_cart, &
         copy_to_castep_cell
    use tssearch_algor,  only: tssearch_algor_QST, tssearch_algor_lst
    use tssearch_utils,  only: tssearch_utils_io_abort, &
         tssearch_utils_get_coordinates, tssearch_utils_fix_coordinates, &
         tssearch_utils_check_coord, tssearch_utils_deallocate

    implicit none

    ! Arguments

    type(castep_model), intent(inout) :: mdl ! The model to be initialised
    character(len=*),  intent(in)     :: output_file

    ! Local Variables
    type(ELEMENT),allocatable,dimension(:) :: elements_tmp
    integer :: ndim_atoms
    integer :: nat
    integer :: nat_reac
    integer :: nat_prod
    integer :: nat_diff
    integer :: ierr
    real (kind=dp), allocatable, dimension(:,:,:) :: xbuff
    real (kind=dp), allocatable, dimension(:)     :: x_reac,x_prod,x_intm,xclean1,xclean2
    real (kind=dp), allocatable, dimension(:)     :: g_reac,g_prod,g_intm
    real (kind=dp)                                :: e_reac,e_prod,e_intm
    real (kind=dp)                                :: p12,p13,p23
    integer,  allocatable, dimension(:)           :: idfix
    logical                                       :: tsfound,period,impose
    integer                                       :: neg,nfix
    integer                                       :: mtsmethod,nstrucs
    integer                                       :: icontinue

    ! First initialize the cells
    !NB we properly nullify any pointers before first usage
    call castep_cell_dealloc(mdl%orig_cell)
    call castep_cell_dealloc(mdl%cell)
    call castep_cell_dealloc(mdl%reac_cell)
    call castep_cell_dealloc(mdl%prod_cell)
    call castep_cell_dealloc(mdl%intm_cell)

    ! Now, allocate the original, model cell, reactant, product and intermediate cells
    ! and copy some dummy contents

    call castep_cell_alloc(mdl%orig_cell)
    call castep_cell_alloc(mdl%cell)
    call castep_cell_alloc(mdl%reac_cell)
    call castep_cell_alloc(mdl%prod_cell)
    call castep_cell_alloc(mdl%intm_cell)

    call castep_cell_copy(current_cell,mdl%orig_cell)
    call castep_cell_copy(current_cell,mdl%cell)
    call castep_cell_copy(current_cell,mdl%reac_cell)
    call castep_cell_copy(current_cell,mdl%prod_cell)
    call castep_cell_copy(current_cell,mdl%intm_cell)

    ! Copy relevant coordinates into appropriate locations
    nstrucs = 0

    !check model has all the bits needed
    if (.not.associated(mdl%forces)) then
       allocate(mdl%forces(1:3,1:mdl%cell%max_ions_in_species,1:mdl%cell%num_species),stat=ierr)
       if (ierr/=0) call tssearch_utils_io_abort('Error in allocating mdl%forces in tssearch_LSTQST')
       mdl%forces=0.0_dp
    end if

    ! Reactant structure - already copied (take current one)
    nstrucs = nstrucs + 1
    nat_reac  = mdl%reac_cell%num_ions

    !Calculate number of degrees of freedom
    nat = nat_reac
    ndim_atoms = 3*nat

    ! prepare a spare array for reading configurations
    allocate(elements_tmp(nat),stat=ierr)
    if (ierr /= 0) call tssearch_utils_io_abort('Error in allocating elements_tmp in tssearch_LSTQST')

    if (pub_on_root) call esdf_init(trim(pub_rootname)//'.dat',ierr)

    ! Product structure
    nstrucs = nstrucs + 1
    call rundat_blocks_exec(elements_tmp,nat,'PRODUCT')
    call copy_to_castep_cell(mdl%prod_cell,elements_tmp)
    nat_prod  = mdl%prod_cell%num_ions

    ! Intermediate structure
    select case (tssearch_lstqst_protocol)
    case ('QST/OPTIMIZATION')
        call rundat_blocks_exec(elements_tmp,nat,'INTERMEDIATE')
        call copy_to_castep_cell(mdl%intm_cell,elements_tmp)
        nstrucs = nstrucs + 1
    end select
    
    if (pub_on_root) call esdf_close()

    !Check number of atoms in reactant and product models
    nat_diff = nat_reac - nat_prod
    if (nat_diff /=0) ierr = 1
    if (ierr/=0) call tssearch_utils_io_abort('Error in tssearch_LSTQST - Unequal number of atoms &
            & in reactant and product')

    !Allocate reactant, product, intermediate arrays   
    !Reactant
    allocate (x_reac(1:ndim_atoms),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating x_reac in tssearch_LSTQST')
    x_reac=0.0_dp
    allocate (g_reac(1:ndim_atoms),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating g_reac in tssearch_LSTQST')
    g_reac=0.0_dp

    !Product
    allocate (x_prod(1:ndim_atoms),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating x_prod in tssearch_LSTQST')
    x_prod=0.0_dp
    allocate (g_prod(1:ndim_atoms),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating g_prod in tssearch_LSTQST')
    g_prod=0.0_dp

    !Intermediate
    allocate (x_intm(1:ndim_atoms),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating x_intm in tssearch_LSTQST')
    x_intm=0.0_dp
    allocate (g_intm(1:ndim_atoms),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating g_intm in tssearch_LSTQST')
    g_intm=0.0_dp

    !Arrays for clean structures
    allocate (xclean1(1:ndim_atoms),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating xclean1 in tssearch_LSTQST')
    xclean1=0.0_dp
    allocate (xclean2(1:ndim_atoms),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating xclean2 in tssearch_LSTQST')
    xclean2=0.0_dp

    !Buffer variable coordinates
    allocate(xbuff(1:3,1:mdl%cell%max_ions_in_species,1:mdl%cell%num_species),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating xbuff in tssearch_LSTQST')
    xbuff=0.0_dp

    !Fixed coordinate information
    allocate (idfix(1:ndim_atoms),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating idfix in tssearch_LSTQST')
    idfix=0
    nfix=0

    !Reactant
    call castep_cell_frac_to_cart(mdl%reac_cell,xbuff)
    call tssearch_utils_get_coordinates(mdl,xbuff,x_reac)
    xclean1 = x_reac      !preserve a clean set of reactant coordinates

    !Product
    call castep_cell_frac_to_cart(mdl%prod_cell,xbuff)
    call tssearch_utils_get_coordinates(mdl,xbuff,x_prod)
    xclean2 = x_prod      !preserve a clean set of product coordinates

    !Intermediate
    call castep_cell_frac_to_cart(mdl%intm_cell,xbuff)
    call tssearch_utils_get_coordinates(mdl,xbuff,x_intm)

    !get fixed coordinate information
    call tssearch_utils_fix_coordinates(nfix,idfix)

    !check the consistency of the structures with fixed atom information and in relation to each other
    call tssearch_utils_check_coord(ndim_atoms,x_reac,x_prod,p12,idfix)

    !set necessary flags before calling main drivers
    impose = .FALSE.        !no superposition for periodic models
    period = .TRUE.         !periodic boundary condition flag
    tsfound = .FALSE.       !transition state flag
    neg = 0                 !energy-gradient counter

    !set the method flags
    select case (tssearch_lstqst_protocol)
    case ('LSTMAXIMUM')
            mtsmethod = 1
            if (nstrucs < 2) call tssearch_utils_io_abort('Error in tssearch_LSTQST - & 
            & Insufficient number of structures provided. Please check the seed.dat file')
    case ('HALGREN-LIPSCOMB')
            mtsmethod = 12
            if (nstrucs < 2) call tssearch_utils_io_abort('Error in tssearch_LSTQST - & 
            & Insufficient number of structures provided. Please check the seed.dat file')
    case ('LST/OPTIMIZATION')
            mtsmethod = 12
            if (nstrucs < 2) call tssearch_utils_io_abort('Error in tssearch_LSTQST - & 
            & Insufficient number of structures provided. Please check the seed.dat file')
    case ('COMPLETELSTQST')
            mtsmethod = 123
            if (nstrucs < 2) call tssearch_utils_io_abort('Error in tssearch_LSTQST - & 
            & Insufficient number of structures provided. Please check the seed.dat file')
    case ('QST/OPTIMIZATION')
            mtsmethod = 3
            if (nstrucs /= 3) call tssearch_utils_io_abort('Error in tssearch_LSTQST - & 
            & Insufficient number of structures provided. Please check the seed.dat file')
    !check the consistency of the structures with fixed atom information and in relation to each other
            call tssearch_utils_check_coord(ndim_atoms,x_reac,x_intm,p13,idfix)
            call tssearch_utils_check_coord(ndim_atoms,x_prod,x_intm,p23,idfix)
    end select

    !select driver
    icontinue=0
    if ( mtsmethod .eq. 1  .or. mtsmethod .eq. 12 .or.  mtsmethod .eq. 123 ) then
        if(pub_on_root)then
           ! Do algorithimic things only on the root node 
           ! all other nodes are "on alert" in the  tssearch_dummy
           ! where the are waiting for the signal to join electronic calculations
           call tssearch_algor_LST(mtsmethod, nat, ndim_atoms, period, impose, & 
                                x_reac, g_reac, e_reac, x_prod, g_prod, e_prod, &
                                x_intm, g_intm, e_intm, xclean1,xclean2, nfix, idfix,&
                                tsfound,neg,output_file,&
                                tssearch_cg_max_iter,tssearch_force_tol,&
                                tssearch_disp_tol)
           ! signal to let other nodes to run:
           call comms_bcast(pub_root_node_id,icontinue)

        else
           call tssearch_dummy(ndim_atoms,x_reac, g_reac)
        endif
        ! let's inform all the nodes about current status:
        call comms_bcast(pub_root_node_id,tsfound)

    end if

    if (.not.tsfound ) then
       if ( mtsmethod .eq. 123 .or.  mtsmethod .eq. 3   ) then
        if(pub_on_root)then
           ! Do algorithimic things only on the root node 
           ! all other nodes are "on alert" in the  tssearch_dummy
           ! where the are waiting for the signal to join electronic calculations
           call tssearch_algor_QST(mtsmethod, nat, ndim_atoms, period, impose, & 
                                x_reac, g_reac, e_reac, x_prod, g_prod, e_prod, &
                                x_intm, g_intm, e_intm, xclean1,xclean2, nfix, idfix,&
                                tsfound,neg,output_file,tssearch_qst_max_iter,&
                                tssearch_cg_max_iter,tssearch_force_tol,&
                                tssearch_disp_tol)
           ! signal to let other nodes to run:
           call comms_bcast(pub_root_node_id,icontinue)
        else
           call tssearch_dummy(ndim_atoms,x_reac, g_reac)
        endif

       end if
    end if

    !Deallocate local arrays
    deallocate (x_reac,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating x_reac in tssearch_LSTQST')
    deallocate (g_reac,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating f_reac in tssearch_LSTQST')
    deallocate (x_prod,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating x_prod in tssearch_LSTQST')
    deallocate (g_prod,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating f_prod in tssearch_LSTQST')
    deallocate (x_intm,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating x_intm in tssearch_LSTQST')
    deallocate (g_intm,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating f_intm in tssearch_LSTQST')
    deallocate (xbuff,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating xbuff in tssearch_LSTQST')
    deallocate (idfix,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating idfix in tssearch_LSTQST')
    deallocate (xclean1,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating xclean1 in tssearch_LSTQST')
    deallocate (xclean2,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating xclean2 in tssearch_LSTQST')
    deallocate(elements_tmp,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating elements_tmp in tssearch_LSTQST')

    !all done
    !Deallocate public arrays in tssearch_utils
    call tssearch_utils_deallocate()

    return
  end subroutine tssearch_LSTQST

  subroutine tssearch_dummy(ndim,x, g)
    !=========================================================================!
    ! This routine is called on each non root node. An execution stops here   !
    ! and program wait for the data from the root node to start force         !
    ! calculations for the geometry defined by the TSSEARCH algorithm (root   !
    ! node). When data are ready variable icontinue=1 is sent to each node.   !
    ! After everything is done icontinue=0 is sent and routine return control ! 
    ! the main flow                                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:(see below)                                                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used: constants,comms,tssearch_utils                            !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   mdl has already been properly read and initialised                    !
    !-------------------------------------------------------------------------!
    ! Written by Alexander Perlov 20/Dec/2006                                 !
    !=========================================================================!
    use constants, only: DP
    use comms, only : comms_bcast,pub_root_node_id
    use tssearch_utils, only : tssearch_utils_energy_gradient
    implicit none
    ! arguments 
    integer, intent(in)           ::ndim         ! number of atom dimensions
    real(kind=dp), intent(inout)  ::x(:)         ! Cartesian coords array
    real(kind=dp), intent(inout)  ::g(:)         ! Gradients array


    !local variables
    real(kind=dp)                 ::energy   ! Total energy
    integer                       :: icontinue
    call comms_bcast(pub_root_node_id,icontinue)
    do while (icontinue == 1) 
       call tssearch_utils_energy_gradient(ndim,x,g,energy)  
       call comms_bcast(pub_root_node_id,icontinue)
    end do
  end  subroutine tssearch_dummy
  
end module tssearch
