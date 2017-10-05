! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-


!============================================================================!
!                                                                            !
!                        Molecular Dynamics Module                           !
!                                                                            !
!----------------------------------------------------------------------------!
! Written by by Simon M.-M. Dubois (May 2010) based on the original          !
! module created by Arash A. Mostofi (May 2006)                              !
!============================================================================!

module md_ions

  use constants, only: dp, stdout

  implicit none

  private

  ! Set of classical particles.
  type, public :: subset 
     integer              :: label
     integer              :: num_atom
     integer, allocatable :: list_atom(:)
  end type subset
  
  ! MD global variables
  real(kind=dp), allocatable, public, save :: r_t(:,:)
  real(kind=dp), allocatable, public, save :: v_t(:,:)        
  real(kind=dp), allocatable, public, save :: a_t(:,:)
  real(kind=dp), allocatable, public, save :: ion_mass(:)

  ! MD properties
  real(kind=dp), public, save :: md_volume
  real(kind=dp), public, save :: md_total_polarisation(3)        
  real(kind=dp), public, save :: md_ionic_polarisation(3)        
  real(kind=dp), public, save :: md_elec_polarisation(3)        

  ! Subsets of atoms
  integer, public, save :: group_num
  type(subset), allocatable, public, save :: groups(:)

  ! Module parameters
  real(kind=dp), parameter :: electron_mass_si = 9.1093826e-31_dp
  real(kind=dp), parameter :: avogadro_si = 6.0221415e23_dp

  ! List of subroutines and functions :
  public :: kinetic_energy
  public :: group_kinetic_energy
  public :: md_allocate_struct
  public :: md_deallocate_struct
  public :: md_initialise_groups
  public :: md_initialise_ion_masses
  public :: md_initialise_positions
  public :: md_initialise_velocities

  contains
       
!==============================================================================!
!==============================================================================!

    function kinetic_energy()
    
      use simulation_cell, only: pub_cell

      implicit none
    
      real(kind=dp) :: kinetic_energy
      integer :: iat, idir
      
      kinetic_energy = 0.0_dp
      do iat = 1, pub_cell%nat
         do idir = 1, 3
            kinetic_energy = kinetic_energy + ion_mass(iat)*v_t(idir,iat)**2
         enddo
      enddo
      kinetic_energy = 0.5_dp*kinetic_energy
    
      return
    
    end function kinetic_energy

!==============================================================================!
!==============================================================================!

    function group_kinetic_energy(group)
    
      implicit none
    
      integer, intent(in) :: group

      real(kind=dp) :: group_kinetic_energy
      integer       :: nat, iat, jat, idir
      
      nat = groups(group)%num_atom
      group_kinetic_energy = 0.0_dp
      do jat = 1, nat
         iat = groups(group)%list_atom(jat)
         do idir = 1, 3
            group_kinetic_energy = group_kinetic_energy + ion_mass(iat)*v_t(idir,iat)**2
         enddo
      enddo
      group_kinetic_energy = 0.5_dp*group_kinetic_energy
    
      return
    
    end function group_kinetic_energy

!==============================================================================!
!==============================================================================!

    subroutine md_allocate_struct(elements)

      use ion,             only: element
      use simulation_cell, only: pub_cell
      use utils,           only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Arguments
      type(ELEMENT), intent(in) :: elements(pub_cell%nat)

      ! Variables
      integer, allocatable :: group_nat(:)
      integer, allocatable :: group_label(:)
      integer, allocatable :: group_idx(:)
      integer :: ngroup
      integer :: iat, ig
      integer :: ierr

      ! Allocate global arrays 
      allocate(r_t(3,pub_cell%nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct (md_main)','r_t',ierr)
      allocate(v_t(3,pub_cell%nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct (md_main)','v_t',ierr)
      allocate(a_t(3,pub_cell%nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct (md_main)','a_t',ierr)
      allocate(ion_mass(pub_cell%nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct (md_main)','ion_mass',ierr)

      r_t=0.0_dp; v_t=0.0_dp; a_t=0.0_dp; ion_mass=0.0_dp

      ! Determine number of subset
      allocate(group_nat(pub_cell%nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct','group_nat',ierr)
      allocate(group_idx(pub_cell%nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct','group_idx',ierr)
      allocate(group_label(pub_cell%nat),stat=ierr)
      call utils_alloc_check('md_allocate_struct','group_label',ierr)

      ! Determine number of groups
      ngroup = 1
      group_label(1) = elements(1)%group_id
      group_idx(1) = 1 
      group_nat(1) = 1 

      group_label(2:pub_cell%nat) = -1
      group_idx(2:pub_cell%nat) = -1
      group_nat(2:pub_cell%nat) = -1

      do iat = 2, pub_cell%nat 
         do ig = 1, ngroup
            if (elements(iat)%group_id == group_label(ig)) then
               group_nat(ig) = group_nat(ig) + 1
               group_idx(iat) = ig
               exit
            endif
         enddo
         if (group_idx(iat) .lt. 0) then
            ngroup = ngroup + 1
            group_label(ngroup) = elements(iat)%group_id
            group_nat(ngroup) = 1
            group_idx(iat) = ngroup 
         endif
  
      enddo
      group_num = ngroup

      ! Allocate atoms
      allocate(groups(group_num),stat=ierr)
      call utils_alloc_check('md_allocate_struct','group_num',ierr)
     
      ! Allocate arrays for each group of atoms
      do ig = 1, group_num
         groups(ig)%label = group_label(ig)
         groups(ig)%num_atom = group_nat(ig)

         allocate(groups(ig)%list_atom(group_nat(ig)),stat=ierr)
         call utils_alloc_check('md_allocate_struct','groups%list_atom',ierr)
      enddo

      deallocate(group_nat,stat=ierr)
      call utils_dealloc_check('md_allocate_struct','group_nat',ierr)
      deallocate(group_idx,stat=ierr)
      call utils_dealloc_check('md_allocate_struct','group_idx',ierr)
      deallocate(group_label,stat=ierr)
      call utils_dealloc_check('md_allocate_struct','group_label',ierr)

      return

    end subroutine md_allocate_struct

!==============================================================================!
!==============================================================================!

    subroutine md_deallocate_struct()

      use utils,           only: utils_dealloc_check

      implicit none

      integer :: ierr, ig

      ! Deallocate global arrays 
      deallocate(ion_mass,stat=ierr)
      call utils_dealloc_check('md_deallocate_struct','ion_mass',ierr)
      deallocate(a_t,stat=ierr)
      call utils_dealloc_check('md_deallocate_struct','a_t',ierr)
      deallocate(v_t,stat=ierr)
      call utils_dealloc_check('md_deallocate_struct','v_t',ierr)
      deallocate(r_t,stat=ierr)
      call utils_dealloc_check('md_deallocate_struct','r_t',ierr)


      ! dellocate arrays for each group of atoms
      do ig = 1, group_num
         deallocate(groups(ig)%list_atom,stat=ierr)
         call utils_dealloc_check('md_deallocate_struct','groups%list_atom',ierr)
      enddo
      deallocate(groups,stat=ierr)
      call utils_dealloc_check('md_deallocate_struct','groups',ierr)

      return

    end subroutine md_deallocate_struct

!==============================================================================!
!==============================================================================!

    subroutine md_initialise_groups(elements)

      use ion,              only: element
      use simulation_cell,  only: pub_cell

      implicit none

      ! Arguments
      type(ELEMENT), intent(in)  :: elements(pub_cell%nat) 
     
      ! Local variables
      integer :: iat, jat, ig

      ! Create subset of atoms
      do ig = 1, group_num
         iat = 0
         do jat=1,pub_cell%nat
            if (elements(jat)%group_id == groups(ig)%label) then
               iat = iat + 1
               groups(ig)%list_atom(iat) = jat
            endif
         enddo
      enddo

    end subroutine md_initialise_groups

!==============================================================================!
!==============================================================================!

    subroutine md_initialise_ion_masses(elements)


      use constants,        only: dp,stdout,periodic_table_mass
      use ion,              only: element
      use simulation_cell,  only: pub_cell
      
      implicit none

      ! Arguments
      type(ELEMENT), intent(in)  :: elements(pub_cell%nat) 
     
      ! Local variables
      integer :: iat
      integer :: atom_Z

      ! Initialise ionic masses (convert masses to atomic units)
      do iat=1,pub_cell%nat
         atom_Z = elements(iat)%atomic_number
         ion_mass(iat) = periodic_table_mass(atom_Z)*1e-3_dp/avogadro_si/electron_mass_si
      enddo
      
    end subroutine md_initialise_ion_masses

!==============================================================================!
!==============================================================================!

    subroutine md_initialise_positions(elements)


      use constants,        only: dp,stdout
      use ion,              only: element
      use simulation_cell,  only: pub_cell
      
      implicit none

      ! Arguments
      type(ELEMENT), intent(in)  :: elements(pub_cell%nat) 
     
      ! Local variables
      integer :: iat

      ! Initialise coordinates
      do iat=1,pub_cell%nat
         r_t(1,iat) = elements(iat)%centre%x
         r_t(2,iat) = elements(iat)%centre%y
         r_t(3,iat) = elements(iat)%centre%z
      enddo

    end subroutine md_initialise_positions

!==============================================================================!
!==============================================================================!

    subroutine md_initialise_velocities(group, temp, elements)


      use comms,            only: pub_on_root, comms_bcast, pub_root_node_id
      use constants,        only: dp, stdout
      use ion,              only: element
      use rundat,           only: md_init_velocities, md_restart
      use services,         only: services_maxboltzdist
      use simulation_cell,  only: pub_cell
      
      implicit none

      ! Arguments
      integer, intent(in)         :: group
      real(kind=dp), intent(in)   :: temp
      type(ELEMENT), intent(in)   :: elements(pub_cell%nat) 
     
      ! Local variables
      integer :: ndof(3)
      real(kind=dp) :: sigma, tscale(3)
      real(kind=dp) :: mtot, vtot(3), ptot(3), mv2tot(3)
      integer :: iat, jat, idir

      ! User defined velocities 
      if ((.not. md_init_velocities) .or. md_restart) then
         do iat = 1, pub_cell%nat
            jat = groups(group)%list_atom(iat)
            v_t(:,jat) = elements(iat)%ion_velocity(:)
         enddo
     
      ! Velocities from a Maxwell-Boltzmann distribution at 
      ! desired temperature
      else 
         ndof = 0
         ptot = 0.0_dp
         mtot = 0.0_dp
         do iat = 1, groups(group)%num_atom
            jat = groups(group)%list_atom(iat)
         
            sigma = sqrt(temp/ion_mass(jat))
            do idir = 1, 3
               ndof(idir) = ndof(idir) + 1
               v_t(idir,jat) = sigma*services_maxboltzdist()
               ptot(idir) = ptot(idir) + ion_mass(jat)*v_t(idir,jat)
            enddo
         
            mtot = mtot + ion_mass(jat) 
         enddo
         
         ! Effective velocity of total mass
         if (mtot .gt. tiny(1.0_dp)) then
            vtot(:) = ptot(:)/mtot
         else
            vtot = 0.0_dp
         endif
         
         ! Correct for centre of mass drift
         ! and compute total KE
         mv2tot = 0.0_dp
         do iat = 1, groups(group)%num_atom
            jat = groups(group)%list_atom(iat)
         
            v_t(:,jat) = v_t(:,jat) - vtot(:)
            do idir = 1, 3
               mv2tot(idir) = mv2tot(idir) + ion_mass(jat)*v_t(idir,jat)**2
            enddo
         
         enddo
         
         ! Scale factor to fix temperature
         tscale = 0.0_dp
         do idir = 1, 3
            if (ndof(idir).gt.0) then
               mv2tot(idir)=mv2tot(idir)/real(ndof(idir),kind=dp)
            else
               mv2tot(idir)=0.0_dp
            endif
            if (mv2tot(idir).gt.tiny(1.0_dp)) tscale(idir)=sqrt(temp/mv2tot(idir))
         enddo
         
         ! Rescale velocities
         do iat = 1, groups(group)%num_atom
            jat = groups(group)%list_atom(iat)
         
            do idir = 1, 3
               v_t(idir,jat) = v_t(idir,jat)*tscale(idir)
            enddo
         enddo

      endif

      ! Make sure all nodes have the same velocity
      call comms_bcast(pub_root_node_id,v_t(:,:))

      return

    end subroutine md_initialise_velocities


!==============================================================================!
!==============================================================================!

 
end module md_ions
    
    
    
    
    
