! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                T S S E A R C H _ U T I L S                                  !
!=============================================================================!
!                                                                             !
! $Id: tssearch_utils_mod.F90,v 1.10 2009/06/23 09:28:18 cks22 Exp $              !
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

module tssearch_utils
  use constants, only : dp
  use ion, only: element
  use simulation_cell, only: castep_model
  implicit none                                 !Impose strong typing

  private                                       !Everything is private ...

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: tssearch_utils_get_forces
  public :: tssearch_utils_energy_gradient
  public :: tssearch_utils_cart_to_frac
  public :: tssearch_utils_write_summary
  public :: tssearch_utils_write_forces
  public :: tssearch_utils_write
  public :: tssearch_utils_write_trajectory
  public :: tssearch_utils_get_coordinates
  public :: tssearch_utils_fix_coordinates
  public :: tssearch_utils_check_coord
  public :: tssearch_utils_initialize
  public :: tssearch_utils_get_images
  public :: tssearch_utils_io_abort
  public :: tssearch_utils_deallocate

  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  type(castep_model),pointer,public,save  :: mdl_ptr
  real(kind=dp), public,save              :: alength  ! lattice parameters: lengths
  real(kind=dp), public,save              :: blength
  real(kind=dp), public,save              :: clength
  real(kind=dp), public,save              :: alpha    ! lattice parameters: angles
  real(kind=dp), public,save              :: beta
  real(kind=dp), public,save              :: gamma
  integer, public,save                    :: ncell
  integer, public,save                    :: nxcell
  integer, public,save                    :: nycell
  integer, public,save                    :: nzcell
  logical, public,save                    :: other_images
  logical, public,save                    :: root_write      ! are we on root node
  integer, dimension(:,:), allocatable, public,save  :: icell
  real(kind=dp), public,save              :: alpha_cos
  real(kind=dp), public,save              :: beta_sin
  real(kind=dp), public,save              :: beta_cos
  real(kind=dp), public,save              :: gamma_sin
  real(kind=dp), public,save              :: gamma_cos
  real(kind=dp), public,save              :: beta_term
  real(kind=dp), public,save              :: gamma_term
  real(kind=dp), public,save              :: energy_conv
  real(kind=dp), public,save              :: force_conv
  real(kind=dp), public,save              :: length_conv
  character (len=30), public,save         :: length_label_ts
  character (len=30), public,save         :: energy_label_ts
  character (len=30), public,save         :: force_label_ts
  logical, public,save                    :: scf_fail
  integer, public,save                    :: ts_stdout
  logical, dimension(:,:,:), allocatable, public,save :: atom_move_ts
  logical, public, save                   :: constrain_ions ! set in tssearch

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  type(element),pointer,dimension(:),save   :: elements_ptr

contains

  subroutine tssearch_utils_initialize(mdl,elements)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, 28/11/2001                                !
    !=========================================================================!
    use comms, only: pub_on_root
    use constants, only: stdout
    use ion, only: element
    use simulation_cell, only: castep_model, castep_cell_cart_lattice_to_abc, &
         pub_cell

    implicit none

    ! Arguments
    type(castep_model), intent(inout), target :: mdl
    type(ELEMENT), intent(inout), target      :: elements(pub_cell%nat)

    !point mdl_ptr to structure
    mdl_ptr => mdl

    !point elements_ptr to structure
    elements_ptr => elements

    !define output logical unit
    ts_stdout = stdout

    !get the lattice parameters (lengths and angles)
    call castep_cell_cart_lattice_to_abc(mdl%cell,alength,blength,clength,alpha,&
                       beta,gamma)

    !get energy, force and length units
    length_conv = 1.0_dp ! io_atomic_to_unit(1.0_dp,length_unit)
    force_conv  = 1.0_dp ! io_atomic_to_unit(1.0_dp,force_unit)
    energy_conv = 1.0_dp ! io_atomic_to_unit(1.0_dp,energy_unit)

    !get energy, force and length labels
    energy_label_ts     = 'Ha'         ! call io_unit_label(energy_unit,energy_label)
    force_label_ts      = 'Ha/Bohr'    ! call io_unit_label(force_unit,force_label)
    length_label_ts     = 'Bohr'       ! call io_unit_label(length_unit,di_length_label)

    !get the root node id. same as on_root
    !output is written ONLY on root node
    root_write = pub_on_root

    return
  end subroutine tssearch_utils_initialize

  subroutine tssearch_utils_get_images(rcut)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, 30/11/2001                                !
    !=========================================================================!
    use constants, only : dp,pi

    implicit none

    ! Arguments
    real(kind=dp), intent(in) :: rcut
 
    ! Local variables
    integer :: i,j,k,ii,ierr
    real(kind=dp) :: xbox,ybox,zbox
    real(kind=dp) :: xcut,ycut,zcut
    real(kind=dp) :: imagemax
    real(kind=dp)::alpha_rad,beta_rad,gamma_rad  ! angles in radians

    !setup some lattice related quantities
    xbox  = alength
    ybox  = blength
    zbox  = clength

    !change angles to radians
    alpha_rad=alpha*pi/180.0_dp
    beta_rad=beta*pi/180.0_dp
    gamma_rad=gamma*pi/180.0_dp

    alpha_cos  = cos(alpha_rad)
    beta_sin   = sin(beta_rad)
    beta_cos   = cos(beta_rad)
    gamma_sin  = sin(gamma_rad)
    gamma_cos  = cos(gamma_rad)
    beta_term  = (alpha_cos - beta_cos*gamma_cos) / gamma_sin
    gamma_term = sqrt(1.0_dp - beta_cos**2 - beta_term**2)

    !project non-cubic box onto cubic
    xcut = xbox * beta_sin * gamma_sin
    ycut = ybox * gamma_sin
    zcut = zbox * beta_sin

    !maximum number of images determined by smallest cell dimension
    imagemax = min(xcut,ycut,zcut)

    !cells in each dimension
    nxcell = int(rcut/xcut)
    nycell = int(rcut/ycut)
    nzcell = int(rcut/zcut)

    !number of cells
    ncell = (2*nxcell + 1) * (2*nycell + 1) * (2*nzcell + 1) - 1

    !set images flag
    if (rcut.le.imagemax.or.ncell.eq.0) then
        other_images = .false.
    else
        other_images = .true.
    end if

    !allocate icell
    allocate(icell(3,ncell),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating icell in &
                          & tssearch_utils_image_setup')
    !initialize icell
    do ii=1,ncell
       icell(1,ii) = 0
       icell(2,ii) = 0
       icell(3,ii) = 0
    end do

    !fill icell array
    ii = 0
    do k = -nzcell,nzcell,1
    do j = -nycell,nycell,1
    do i = -nxcell,nxcell,1
      if (k.ne.0 .or. j.ne.0 .or. i.ne.0) then
          ii = ii + 1
          icell(1,ii) = i
          icell(2,ii) = j
          icell(3,ii) = k
      end if
    end do
    end do
    end do
    
    return
  end subroutine tssearch_utils_get_images

  subroutine tssearch_utils_get_forces(mdl,elements)
    !=========================================================================!
    ! Find ground state energy, forces and stresses for given model.          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mdl (inout) : The model to be checked.                                !
    !   elements (inout) : structure on ONETEP convention.                    !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman 27/03/06                                       !
    !=========================================================================!
#ifdef debug 
    use constants,        only: stdout
#endif
    use energy_and_force, only: energy_and_force_calculate
    use forces, only: forces_apply_constraints
    use ion, only: element
    use simulation_cell, only: castep_model, pub_cell, castep_cell_2_elements
    
    implicit none

    ! Arguments
    type(castep_model), intent(inout) :: mdl ! The model to check for the grounds state
    type(element), intent(inout)      :: elements(1:pub_cell%nat)

    !local variables
    integer       :: ii,is

#ifdef debug
    write (stdout,*) 'starting tssearch_utils_get_forces'
#endif


    !-------------------------------------------------------------------------!
    !Now we must find the ground state wvfn for this model                    !
    !-------------------------------------------------------------------------!
    if (.not.mdl%found_ground_state) then

       ! aam: Update ionic coordinates in elements with new ones from current_cell
       call castep_cell_2_elements(mdl%cell,elements)

       ! aam: Calculate new ground state energy and forces
       call energy_and_force_calculate(mdl%total_energy,mdl%forces,elements)

       ! aam: In the presence of ionic constraints keep copy of unconstrained forces
       if (constrain_ions) then 
          do ii=1,pub_cell%nat                          ! Each ion has its own species
             do is=1,mdl%cell%max_ions_in_species   ! ie max_ions_species = 1 
                mdl%orig_forces(1,is,ii) = mdl%forces(1,is,ii)
                mdl%orig_forces(2,is,ii) = mdl%forces(2,is,ii)
                mdl%orig_forces(3,is,ii) = mdl%forces(3,is,ii)
             enddo
          enddo
       endif

       ! aam: Apply ionic constraints
       if (constrain_ions) call forces_apply_constraints(mdl%forces,elements)

       mdl%found_ground_state = .true.
    end if

    return
  end subroutine tssearch_utils_get_forces

  subroutine tssearch_utils_get_coordinates(model,xin,xout)
    !==========================================================================!
    ! Purpose:                                                                 !
    ! This routine just packs ionic Cartesian coords                           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! model(in): needed just for indexing of atoms                             !
    ! xin(in): Cartesian coordinates in expanded form                          !
    ! xout(out): array which returns packed Cartesian coordinates              !
    !--------------------------------------------------------------------------!
    ! Key internal variables: none                                             !
    !--------------------------------------------------------------------------!
    ! Written by Niri Govind, version 0.1, 26/11/01                            !
    !==========================================================================!

    use simulation_cell, only: castep_model
    
    implicit none

    type(castep_model), intent(in)               :: model
    real(kind=dp), dimension(:,:,:), intent(in)  :: xin   ! Input Cartesian coords array
    real(kind=dp), dimension(:),intent(out)      :: xout  ! Packed Cartesian coords array
    integer::i,i0
    integer::ispec
    integer::iatom

    !Extract coordinates
    i0 = 0
    do ispec=1,model%reac_cell%num_species
       do iatom=1,model%reac_cell%num_ions_in_species(ispec)
          do i=1,3
              i0 = i0 + 1
              xout(i0) = xin(i,iatom,ispec)
          end do
       end do
    end do

    return
  end subroutine tssearch_utils_get_coordinates

  subroutine tssearch_utils_fix_coordinates(nfix,idfix)
    !==========================================================================!
    ! Purpose:                                                                 !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !--------------------------------------------------------------------------!
    ! Parent module variables used: castep_cell                                  !
    !--------------------------------------------------------------------------!
    ! Modules used: io                                                         !
    !--------------------------------------------------------------------------!
    ! Key internal variables: none                                             !
    !--------------------------------------------------------------------------!
    ! Necessary conditions:                                                    !
    !--------------------------------------------------------------------------!
    ! Written by Niri Govind, version 0.1, 12/12/01                            !
    !==========================================================================!

    implicit none

    ! Arguments
    integer, intent(out) :: idfix(:)
    integer, intent(out) :: nfix

    ! Local variables
    integer :: nsp,ni,i,i0,iatom
    logical :: fixed
    integer :: ierr

    nfix = 0
    idfix = 0

    !allocate logical move array and initialize
    allocate(atom_move_ts(1:3,1:mdl_ptr%cell%max_ions_in_species,1:mdl_ptr%cell%num_species),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating atom_move_ts in &
                  &  tssearch_utils_fix_coordinates')
    atom_move_ts = .TRUE.

    ! get fixed coordinate information
    ! this only assumes diagonal components as we look for individual coordinates only
    i0 = 0
    iatom = 0

    do nsp = 1,mdl_ptr%cell%num_species
      do ni = 1,mdl_ptr%cell%num_ions_in_species(nsp)
        iatom = iatom + 1
        fixed = (elements_ptr(iatom)%ion_constraint_type .eq. 'FIXED')
        if (fixed) nfix = nfix + 1

        do i = 1,3
          i0 = i0 + 1
          if (fixed) then
             idfix(i0) = 1
             atom_move_ts(i,ni,nsp) = .FALSE.
          end if
        end do
      end do
    end do

    return
  end subroutine tssearch_utils_fix_coordinates

  subroutine tssearch_utils_check_coord(ndim,x1,x2,p12,idfix)
    !==========================================================================!
    ! Purpose:                                                                 !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !--------------------------------------------------------------------------!
    ! Parent module variables used: castep_cell                                  !
    !--------------------------------------------------------------------------!
    ! Modules used: io                                                         !
    !--------------------------------------------------------------------------!
    ! Key internal variables: none                                             !
    !--------------------------------------------------------------------------!
    ! Necessary conditions:                                                    !
    !--------------------------------------------------------------------------!
    ! Written by Niri Govind, version 0.1, 12/12/01                            !
    !==========================================================================!
    use constants, only : dp, stdout

    implicit none

    ! Arguments
    real(kind=dp), intent(inout) :: x1(:),x2(:)
    integer, intent(inout) :: idfix(:),ndim

    ! Local variables
    real(kind=dp) :: deltax,p12
    real(kind=dp) :: xdiff,ydiff,zdiff
    real (kind=dp), allocatable, dimension(:) :: xr,yr,zr
    real (kind=dp), allocatable, dimension(:) :: xp,yp,zp
    
    integer :: i,ierr,nat

    !get number of atoms
    nat = ndim/3

    !allocate temporary variable coordinates
    allocate(xr(nat),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating xr in tssearch_utils_check_coord')
    allocate(yr(nat),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating yr in tssearch_utils_check_coord')
    allocate(zr(nat),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating zr in tssearch_utils_check_coord')
    allocate(xp(nat),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating xp in tssearch_utils_check_coord')
    allocate(yp(nat),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating yp in tssearch_utils_check_coord')
    allocate(zp(nat),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating zp in tssearch_utils_check_coord')

    deltax = 1.0e-03_dp
    p12 = 0.0_dp
    do i=1,ndim
        p12 = p12 + (x1(i)-x2(i))*(x1(i)-x2(i))
    enddo
    p12 = sqrt(p12)

    !check fixed coordinate consistency
    do i=1,nat
      !reactant
      xr(i) = x1(3*i-2)
      yr(i) = x1(3*i-1)
      zr(i) = x1(3*i)

      !product
      xp(i) = x2(3*i-2)
      yp(i) = x2(3*i-1)
      zp(i) = x2(3*i)

      !coordinate differences
      xdiff = dabs(xp(i)-xr(i))
      ydiff = dabs(yp(i)-yr(i))
      zdiff = dabs(zp(i)-zr(i))

      if ( (xdiff >= deltax .and. idfix(3*i-2) == 1) .or. &
           (ydiff >= deltax .and. idfix(3*i-1) == 1) .or. &
           (zdiff >= deltax .and. idfix(3*i) == 1)  ) then
            ierr = 1
            write(stdout,'(a)') 'Error: Fixed coordinate mismatch'
            write(stdout,'(a,3(3x,f10.4))') &
             'Fixed Reactant Atom: ',xr(i)*length_conv,yr(i)*length_conv,zr(i)*length_conv
            write(stdout,'(a,3(3x,f10.4))') &
             'Fixed Product Atom:  ',xp(i)*length_conv,yp(i)*length_conv,zp(i)*length_conv
            write(stdout,'(a,f10.4,1x,a)') &
             'Fixed coordinates should match in both structures to within ',deltax,length_label_ts
            write(stdout,'(a)') 'Please check coordinates in the cell file' 
            if (ierr/=0) call tssearch_utils_io_abort('Error detected in tssearch_utils_check_coord')
      end if
    end do

    !check p vectors for structures
    if (p12 <= deltax) then
       ierr = 1
       write(stdout,'(a)') 'Error: Reactant and Product structures may be the same or very close to each other' 
       write(stdout,'(a)') 'Please check coordinates in the cell file' 
       if (ierr/=0) call tssearch_utils_io_abort('Error detected in tssearch_utils_check_coord')
    end if

    !all done
    deallocate(xr,yr,zr,xp,yp,zp)

    return
  end subroutine tssearch_utils_check_coord

  subroutine tssearch_utils_cart_to_frac(r_cart,r_frac)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   Constants                                                             !
    !   Cell                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, 27/11/2001                                !
    !=========================================================================!
    use constants, only: DP
    use simulation_cell, only: castep_cell, castep_cell_dealloc, &
         castep_cell_cart_to_frac, castep_cell_nullify, castep_cell_alloc, &
         castep_cell_copy

    implicit none
  
    ! Arguments

    real(kind=dp),dimension(:),intent(in)  :: r_cart    ! Cartesian coords array
    real(kind=dp),dimension(:),intent(out) :: r_frac    ! Fractional coords array

    ! Local variables

    real (kind=dp) :: xtmp(3)  ! temporary array
    type(castep_cell) :: tmp_cell ! temporary cell
    integer :: i,i0
    integer :: ispec
    integer :: iatom

    ! allocate and initialize tmp_cell
    call castep_cell_nullify(tmp_cell)
    call castep_cell_alloc(tmp_cell)
    call castep_cell_copy(mdl_ptr%cell,tmp_cell)
    tmp_cell%ionic_positions(:,:,:)=0.0_dp
    xtmp=0.0_dp

    !fill temporary array 
    i0 = 0
    do ispec=1,tmp_cell%num_species
       do iatom=1,tmp_cell%num_ions_in_species(ispec)
          do i=1,3
              i0 = i0 + 1
              xtmp(i) = r_cart(i0)
          end do
          ! positions are in cartesian coords - now convert to fractional
          call castep_cell_cart_to_frac(tmp_cell,xtmp,tmp_cell%ionic_positions(:,iatom,ispec))
       end do
    end do
    
    !pass back to main routine as fractional
    i0 = 0
    do ispec=1,tmp_cell%num_species
      do iatom=1,tmp_cell%num_ions_in_species(ispec)
          do i=1,3
            i0 = i0 + 1
            r_frac(i0) = tmp_cell%ionic_positions(i,iatom,ispec)
           end do
      end do
    end do

    !all done
    !get rid of tmp_cell, xtmp
    call castep_cell_dealloc(tmp_cell)

    return
  end subroutine tssearch_utils_cart_to_frac

  subroutine tssearch_utils_energy_gradient(ndim,r,gradient,energy)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   Constants                                                             !
    !   Parameters                                                            !
    !   Cell                                                                  !
    !   Model                                                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, 27/11/2001                                !
    !=========================================================================!
    use comms, only : comms_bcast, pub_on_root, pub_root_node_id
    use constants, only: DP
    use simulation_cell, only: castep_cell, castep_model_cell_changed, &
         castep_cell_nullify, castep_cell_alloc, castep_cell_copy, &
         castep_cell_cart_to_frac, castep_cell_dealloc

    implicit none
  
    ! Arguments

    integer, intent(in)::ndim         ! number of atom dimensions
    real(kind=dp), intent(inout)::r(:)       ! Cartesian coords array
    real(kind=dp), intent(out)::energy       ! Total energy
    real(kind=dp), intent(out)::gradient(:)  ! Gradients array

    ! Local variables

    real (kind=dp) xtmp(3)  ! temporary array
    type(castep_cell) :: tmp_cell ! temporary cell
    integer :: i,i0
    integer :: ispec
    integer :: iatom
    integer :: icontinue
    
    if(pub_on_root)then
       icontinue=1
       call comms_bcast(pub_root_node_id,icontinue)
    endif
    call comms_bcast(pub_root_node_id,r,ndim)

    ! allocate and initialize tmp_cell
    call castep_cell_nullify(tmp_cell)
    call castep_cell_alloc(tmp_cell)
    call castep_cell_copy(mdl_ptr%cell,tmp_cell)
    tmp_cell%ionic_positions(:,:,:)=0.0_dp
    xtmp=0.0_dp

    !update coordinates
    i0 = 0
    do ispec=1,tmp_cell%num_species
       do iatom=1,tmp_cell%num_ions_in_species(ispec)
          do i=1,3
              i0 = i0 + 1
              xtmp(i) = r(i0)
          end do
          !positions are in cartesian coords - now convert to fractional
          call castep_cell_cart_to_frac(tmp_cell,xtmp,tmp_cell%ionic_positions(:,iatom,ispec))
       end do
    end do

    !copy the information
    call castep_cell_copy(tmp_cell,mdl_ptr%cell)

    !... and set the appropriate flags in model
    call castep_model_cell_changed(mdl_ptr)

    !call energy-gradient calculator
    scf_fail = .FALSE.
    call tssearch_utils_get_forces(mdl_ptr,elements_ptr) 

    !extract total energy and gradients
    energy = mdl_ptr%total_energy

    ! get the gradients
    i0 = 0
    do ispec=1,tmp_cell%num_species
       do iatom=1,tmp_cell%num_ions_in_species(ispec)
          do i=1,3
              i0 = i0 + 1
              gradient(i0) = mdl_ptr%forces(i,iatom,ispec)
              gradient(i0) = -gradient(i0)             ! note sign flip
          end do
       end do
    end do

    !all done
    !get rid of tmp_cell, xtmp
    call castep_cell_dealloc(tmp_cell)

    return
  end subroutine tssearch_utils_energy_gradient

  subroutine tssearch_utils_write_summary(nmode,file_unit,energy_reac,energy_prod,&
                  energy_intm,pathvalue)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   Constants                                                             !
    !   Parameters                                                            !
    !   Cell                                                                  !
    !   Model                                                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, 06/12/2001                                !
    !=========================================================================!

    implicit none
  
    ! Arguments

    integer, intent(in)            :: nmode           ! summary mode !3= Transition State
    integer, intent(in)            :: file_unit       ! output file unit
    real(kind=dp), intent(in)      :: energy_reac     ! total energy of reactant !also used as local variable
    real(kind=dp), intent(in)      :: energy_prod     ! total energy of product
    real(kind=dp), intent(in)      :: energy_intm     ! total energy of intermediate
    real(kind=dp), intent(in)      :: pathvalue       ! pathvalue

    !LST maximum summary
    if (nmode == 1) then
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,'(a)')             ' LST Maximum Found!'
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of reactant:           ',energy_reac*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of product:            ',energy_prod*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of LST maximum:        ',energy_intm*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5)')       ' Location of LST maximum:      ',pathvalue
      write(file_unit,'(a,f18.5,1x,a)')  ' Barrier from reactant:        ',(energy_intm-energy_reac)*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Barrier from product:         ',(energy_intm-energy_prod)*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of reaction:           ',(energy_prod-energy_reac)*energy_conv,energy_label_ts
      !write(file_unit,'(a,i6)')       ' Cumulative number of energy/gradient calls: ',neg
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,*)
    end if

    !QST maximum summary
    if (nmode == 2) then
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,'(a)')             ' QST Maximum Found!'
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of reactant:           ',energy_reac*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of product:            ',energy_prod*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of QST maximum:        ',energy_intm*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5)')       ' Location of QST maximum:      ',pathvalue
      write(file_unit,'(a,f18.5,1x,a)')  ' Barrier from reactant:        ',(energy_intm-energy_reac)*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Barrier from product:         ',(energy_intm-energy_prod)*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of reaction:           ',(energy_prod-energy_reac)*energy_conv,energy_label_ts
      !write(file_unit,'(a,i6)')       ' Cumulative number of energy/gradient calls: ',neg
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,*)
    end if

    !Transition state summary
    if (nmode == 3) then
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,'(a)')             ' Transition State Found!'
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of reactant:           ',energy_reac*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of product:            ',energy_prod*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of transition state:   ',energy_intm*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5)')       ' Location of transition state: ',pathvalue
      write(file_unit,'(a,f18.5,1x,a)')  ' Barrier from reactant:        ',(energy_intm-energy_reac)*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Barrier from product:         ',(energy_intm-energy_prod)*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of reaction:           ',(energy_prod-energy_reac)*energy_conv,energy_label_ts
      !write(file_unit,'(a,i6)')       ' Cumulative number of energy/gradient calls: ',neg
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,*)
    end if

    !Halgren-Lipscomb or LST/Optimization summary
    if (nmode == 4) then
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of reactant:           ',energy_reac*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of product:            ',energy_prod*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of intermediate:       ',energy_intm*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5)')       ' Location of intermediate:     ',pathvalue
      write(file_unit,'(a,f18.5,1x,a)')  ' Barrier from reactant:        ',(energy_intm-energy_reac)*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Barrier from product:         ',(energy_intm-energy_prod)*energy_conv,energy_label_ts
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy of reaction:           ',(energy_prod-energy_reac)*energy_conv,energy_label_ts
      !write(file_unit,'(a,i6)')       ' Cumulative number of energy/gradient calls: ',neg
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,*)
    end if

    !General summary
    if (nmode == 5) then
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,'(a,f18.5)')       ' Reactant                      '
      write(file_unit,'(a,f18.5,1x,a)')  ' Path coordinate:              ',pathvalue
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy:                       ',energy_reac*energy_conv,energy_label_ts
      !write(file_unit,'(a,i6)')       ' Cumulative number of energy/gradient calls: ',neg
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,*)
    end if

    if (nmode == 6) then
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,'(a,f18.5)')       ' Product                       '
      write(file_unit,'(a,f18.5,1x,a)')  ' Path coordinate:              ',pathvalue
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy:                       ',energy_reac*energy_conv,energy_label_ts
      !write(file_unit,'(a,i6)')       ' Cumulative number of energy/gradient calls: ',neg
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,*)
    end if

    if (nmode == 7) then
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,'(a,f18.5)')       ' Intermediate                  '
      write(file_unit,'(a,f18.5,1x,a)')  ' Path coordinate:              ',pathvalue
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy:                       ',energy_reac*energy_conv,energy_label_ts
      !write(file_unit,'(a,i6)')       ' Cumulative number of energy/gradient calls: ',neg
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,*)
    end if

    if (nmode == 8) then
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,'(a,f18.5,1x,a)')  ' Path coordinate:              ',pathvalue
      write(file_unit,'(a,f18.5,1x,a)')  ' Energy:                       ',energy_reac*energy_conv,energy_label_ts
      !write(file_unit,'(a,i6)')       ' Cumulative number of energy/gradient calls: ',neg
      write(file_unit,'(a)') ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(file_unit,*)
    end if

    return
  end subroutine tssearch_utils_write_summary

  subroutine tssearch_utils_write_forces(ndim,file_unit,r,gradient)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   Constants                                                             !
    !   Parameters                                                            !
    !   Cell                                                                  !
    !   Model                                                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, 05/04/2002                                !
    !=========================================================================!
    use constants, only: DP
    use simulation_cell, only: castep_cell,castep_cell_nullify,&
         castep_cell_copy, castep_cell_alloc, castep_cell_dealloc

    implicit none
 
    ! Arguments

    integer, intent(in)            :: ndim         ! number of atom dimensions
    integer, intent(in)            :: file_unit    ! output file unit
    real(kind=dp), intent(in)      :: r(:)         ! Cartesian coords array
    real(kind=dp), intent(in)      :: gradient(:)  ! Gradients array

    ! Local variables

    real (kind=dp), allocatable, dimension(:,:,:) :: xtmp  ! temporary array
    real (kind=dp), allocatable, dimension(:,:,:) :: gtmp  ! temporary array
    type(castep_cell) :: tmp_cell ! temporary cell
    integer :: i,i0,j
    integer :: ispec
    integer :: iatom
    integer :: ierr
    real(kind=dp) :: rfrac(ndim) ! Fractional coords array

    ! allocate and initialize tmp_cell
    call castep_cell_nullify(tmp_cell)
    call castep_cell_alloc(tmp_cell)
    call castep_cell_copy(mdl_ptr%cell,tmp_cell)
    tmp_cell%ionic_positions(:,:,:)=0.0_dp

    !allocate and update temporary variable coordinates
    allocate(xtmp(1:3,1:tmp_cell%max_ions_in_species,1:tmp_cell%num_species),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating xtmp in tssearch_utils_write_forces')
    xtmp=0.0_dp
    allocate(gtmp(1:3,1:tmp_cell%max_ions_in_species,1:tmp_cell%num_species),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating gtmp in tssearch_utils_write_forces')
    gtmp=0.0_dp

    !convert from Cartesian to Fractionals
    rfrac=0.0_dp
    call tssearch_utils_cart_to_frac(r,rfrac)

    !transfer variables
    i0 = 0
    do ispec=1,tmp_cell%num_species
       do iatom=1,tmp_cell%num_ions_in_species(ispec)
          do i=1,3
              i0 = i0 + 1
              xtmp(i,iatom,ispec) = rfrac(i0)
              gtmp(i,iatom,ispec) = -gradient(i0)*force_conv    !forces = -gradient
          end do
       end do
    end do

    !write coordinates and forces to the output file
    write(file_unit,*)'           ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(file_unit,'(a)')'                             Coordinates: Fractional components'
    write(file_unit,*)'           ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(file_unit,*)'           +  Element    Atom                                                 +'
    write(file_unit,*)'           +            Number           x          y          z              +'
    write(file_unit,*)'           ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    do i=1,tmp_cell%num_species
     do j=1,tmp_cell%num_ions_in_species(i)
            write(file_unit,3)tmp_cell%species_symbol(i), &
                 & i,xtmp(1,j,i),xtmp(2,j,i),xtmp(3,j,i),"        "
     end do
    end do
    write(file_unit,*)'           ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(file_unit,'(a,a)')'                           Forces: Cartesian components: ',force_label_ts
    write(file_unit,*)'           ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(file_unit,*)'           +  Element    Atom                                                 +'
    write(file_unit,*)'           +            Number           x          y          z              +'
    write(file_unit,*)'           ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    do i=1,tmp_cell%num_species
     do j=1,tmp_cell%num_ions_in_species(i)
          if (atom_move_ts(1,j,i)) then
            write(file_unit,3)tmp_cell%species_symbol(i), &
                 &  i,gtmp(1,j,i),gtmp(2,j,i),gtmp(3,j,i),"        "
          else
            write(file_unit,3)tmp_cell%species_symbol(i), &
                 &  i,gtmp(1,j,i),gtmp(2,j,i),gtmp(3,j,i),"(cons'd)"
          end if    
     end do
    end do
    write(file_unit,*)'           ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(file_unit,*)

    !all done
    !get rid of tmp_cell, xtmp, gtmp
    call castep_cell_dealloc(tmp_cell)
    deallocate (xtmp,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating xtmp in tssearch_utils_write_forces')
    deallocate (gtmp,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating gtmp in tssearch_utils_write_forces')

3   format(11x,' +',a6,1x,i8,7x,3f11.6,1x,a8,'  +')

    return
  end subroutine tssearch_utils_write_forces

  subroutine tssearch_utils_write(r,gradient,energy,neg,output_file,&
                       status_name,status_number,pathvalue)
    !=========================================================================!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   Constants                                                             !
    !   Parameters                                                            !
    !   Cell                                                                  !
    !   Model                                                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, 06/12/2001                                !
    !=========================================================================!
    use constants, only: DP
    use simulation_cell, only: castep_cell, castep_cell_nullify, &
         castep_cell_alloc, castep_cell_copy, castep_cell_dealloc

    implicit none
  
    ! Arguments

    integer, intent(in)            :: neg          ! number of energy-gradient calls
    real(kind=dp), intent(in)      :: r(:)         ! Cartesian coords array
    real(kind=dp), intent(in)      :: energy       ! Total energy
    real(kind=dp), intent(in)      :: gradient(:)  ! Gradients array
    character(len=*), intent(in)   :: output_file  !designated output file
    character(len=*), intent(in)   :: status_name
    integer, intent(in)            :: status_number
    real(kind=dp), intent(in)      :: pathvalue

    ! Local variables

    real (kind=dp), allocatable, dimension(:,:,:) :: xtmp  ! temporary array
    real (kind=dp), allocatable, dimension(:,:,:) :: gtmp  ! temporary array
    type(castep_cell) :: tmp_cell ! temporary cell
    integer :: i,i0
    integer :: ispec
    integer :: iatom
    integer :: ierr

    ! allocate and initialize tmp_cell
    call castep_cell_nullify(tmp_cell)
    call castep_cell_alloc(tmp_cell)
    call castep_cell_copy(mdl_ptr%cell,tmp_cell)
    tmp_cell%ionic_positions(:,:,:)=0.0_dp

    !allocate and update temporary variable coordinates
    allocate(xtmp(1:3,1:tmp_cell%max_ions_in_species,1:tmp_cell%num_species),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating xtmp in tssearch_utils_write')
    xtmp=0.0_dp
    allocate(gtmp(1:3,1:tmp_cell%max_ions_in_species,1:tmp_cell%num_species),stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in allocating gtmp in tssearch_utils_write')
    gtmp=0.0_dp

    !transfer variables
    i0 = 0
    do ispec=1,tmp_cell%num_species
       do iatom=1,tmp_cell%num_ions_in_species(ispec)
          do i=1,3
              i0 = i0 + 1
              xtmp(i,iatom,ispec) = r(i0)
              gtmp(i,iatom,ispec) = -gradient(i0)    !forces = -gradient
          end do
       end do
    end do

    !call trajectory writer to write out the frame
    !note that the forces are written out to the trajectory file
    call tssearch_utils_write_trajectory(mdl_ptr,xtmp,gtmp,energy,neg,output_file,&
                   status_name,status_number,pathvalue)

    !all done
    !get rid of tmp_cell, xtmp, gtmp
    call castep_cell_dealloc(tmp_cell)
    deallocate (xtmp,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating xtmp in tssearch_utils_write')
    deallocate (gtmp,stat=ierr)
    if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating gtmp in tssearch_utils_write')

    return
  end subroutine tssearch_utils_write

  subroutine tssearch_utils_write_trajectory(mdl,pos,grad,energy,neg,output_file,&
                status_name,status_number,pathvalue)
    !=========================================================================!
    ! Output relevant data to specified trajectory file                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   on_root                                                               !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, 06/12/2001                                !
    !=========================================================================!
    use comms, only: comms_abort
    use simulation_cell, only: castep_model
    use utils, only : utils_unit

    implicit none

    ! Arguments
    type(castep_model), intent(inout):: mdl
    real (kind=dp),intent(in)        :: pos(:,:,:) !  cartesian coordinate array
    real (kind=dp),intent(in)        :: grad(:,:,:) ! gradient array
    real(kind=dp), intent(in)        :: energy       ! Total energy
    integer                          :: neg
    character(len=*),  intent(in)    :: output_file  !designated output file
    character(len=*),  intent(in)    :: status_name  
    integer, intent(in)              :: status_number
    real(kind=dp), intent(in)        :: pathvalue  

    ! Local variables
    integer                          :: i,iatom,ispec

    !for filehandling we need
    integer                            :: out_unit, io_status
    character(len=12)                  :: action,form,stat,position,access


    action   = 'WRITE'
    form     = 'FORMATTED'
    stat     = 'UNKNOWN'
    position = 'REWIND'
    access   = 'SEQUENTIAL'
    
    if (neg==1) then
       !create a new output file
       stat     ='REPLACE'
    else
       !append to existing output file
       position ='APPEND'
    end if
    
    if (root_write) then
       
       ! Find next available unit specifier
       out_unit = utils_unit()
       
       ! open output file
       open(unit=out_unit,iostat=io_status,file=output_file,status=stat,&
            access=access,form=form,position=position,action=action)

       if (io_status/=0) call tssearch_utils_io_abort &
         ('Error in tssearch_utils_write_trajectory: failed to open file '//trim(output_file))

      !now write out all the relevant data (root node only)
       write (out_unit,1) status_name, status_number, pathvalue
       !in eV
       !write (out_unit,2) io_atomic_to_unit(energy,energy_unit),io_atomic_to_unit(energy,energy_unit)
       !in Hartrees
       write (out_unit,2) energy,energy

       !cell stuff
       !cell vectors
       do i=1,3
             !in angstroms
             !write (out_unit,5) io_atomic_to_unit(mdl%cell%real_lattice(i,1),length_unit), &
             !     & io_atomic_to_unit(mdl%cell%real_lattice(i,2),length_unit), &
             !     & io_atomic_to_unit(mdl%cell%real_lattice(i,3),length_unit)
             !in atomic units
             write (out_unit,5) mdl%cell%real_lattice(i,:)
       end do

       !write all the ion stuff
       !ionic positions
       do ispec=1,mdl%cell%num_species
             do iatom=1,mdl%cell%num_ions_in_species(ispec)
                !in angstroms
                !write (out_unit,3) mdl%cell%species_symbol(ispec), iatom, &
                !     & io_atomic_to_unit(pos(1,iatom,ispec),length_unit), &
                !     & io_atomic_to_unit(pos(2,iatom,ispec),length_unit), &
                !     & io_atomic_to_unit(pos(3,iatom,ispec),length_unit)
                !in atomic units
                write (out_unit,3) mdl%cell%species_symbol(ispec),iatom, &
                      & pos(1,iatom,ispec),pos(2,iatom,ispec),pos(3,iatom,ispec)
             end do
       end do

       !ionic forces
       do ispec=1,mdl%cell%num_species
             do iatom=1,mdl%cell%num_ions_in_species(ispec)
                !in units of eV/A
                !write (out_unit,4) mdl%cell%species_symbol(ispec),iatom,  &
                !     & io_atomic_to_unit(grad(1,iatom,ispec),force_unit), &
                !     & io_atomic_to_unit(grad(2,iatom,ispec),force_unit), &
                !     & io_atomic_to_unit(grad(3,iatom,ispec),force_unit)
                !in Hartree/Bohr
                write (out_unit,4) mdl%cell%species_symbol(ispec),iatom, &
                      & grad(1,iatom,ispec),grad(2,iatom,ispec),grad(3,iatom,ispec)
             end do
       end do

       !blank line to signal end of this iter (datablock size variable)
       write (out_unit,*) ' '

       close (unit=out_unit)

    end if

1    format(1x,a3,1x,i4,3x,es18.8e3)
2    format(9x,2(3x,es18.8e3),21x,'  <-- E')
3    format(1x,a3,1x,i4,3(3x,es18.8e3),'  <-- R')
4    format(1x,a3,1x,i4,3(3x,es18.8e3),'  <-- F')
5    format(9x,3(3x,es18.8e3),'  <-- h')

    return
  end subroutine tssearch_utils_write_trajectory

  subroutine tssearch_utils_io_abort(message)
    !=========================================================================!
    ! Description:                                                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, Dec 17 2001                               !
    !=========================================================================!
    use comms, only : comms_abort
    use constants, only: stdout

    implicit none

    character(len=*), intent(in) :: message
    !-------------------------------------------------------------------------!

    if (root_write) then
       write (stdout,'(a)') message
    end if
    call comms_abort


    return
  end subroutine tssearch_utils_io_abort

  subroutine tssearch_utils_deallocate
    !=========================================================================!
    ! Description: Deallocate all public variables allocated in tssearch_utils!
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Niri Govind, v0.1, May 15 2002                               !
    !=========================================================================!

    implicit none

    ! Local variables

    integer :: ierr

    if(allocated(icell))then
      deallocate(icell,stat=ierr)
      if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating icell in tssearch_utils_deallocate')  
    endif
    
    if(allocated(icell))then
       deallocate(atom_move_ts,stat=ierr)
       if (ierr/=0) call tssearch_utils_io_abort('Error in deallocating atom_move_ts in tssearch_utils_deallocate')  
    endif
    
    return
  end subroutine tssearch_utils_deallocate


end module tssearch_utils
