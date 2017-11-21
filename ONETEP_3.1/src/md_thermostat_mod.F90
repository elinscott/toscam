! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*- !

!============================================================================!
!                                                                            !
!                 Molcular Dynamics Thermostat Module                        !
!                                                                            !
! Thermostat subroutines for a system of classical particles coupled to an   !
! thermal bath. Includes ANDERSON/LANGEVIN/NOSE-HOOVER/NOSE-HOOVER-CHAIN     !
! schemes for the description of the coupling.                               !
!                                                                            !
! Calls to thermostat subroutines should be interleaved with the velocity    !
! Verlet algorithm as follow :                                               !
!                                                                            ! 
!  thermo1                                                                   !
!     v(t+dt/2) = v(t) + a(t)*dt/2                                           !
!  thermo2                                                                   !
!     r(t+dt) = r(t) + v(t+dt/2)*dt                                          !
!     call energy_and_force_calculate()                                      !
!  thermo3                                                                   !
!     v(t+dt) = v(t+dt/2) + a(t+dt)*dt/2                                     !
!  thermo4                                                                   !
!                                                                            !
!----------------------------------------------------------------------------!
! Written by Simon M.-M. Dubois (May 2011)                                   !
!============================================================================!



module md_thermostat

  use constants, only: dp, stdout

  implicit none

  private

  ! Type definition for the thermostat.
  type, public :: thermostat

     ! General parameters
     integer         :: type        ! 0 = NONE, 1 = ANDERSEN,
                                    ! 2 = LANGEVIN, 3 = NOSE_HOOVER
     real(kind=DP)   :: T           ! Thermostat temperature
     real(kind=DP)   :: Tinit       ! Thermostat temperature
     real(kind=DP)   :: Tgrad       ! Temperature gradient
     integer         :: group       ! Subset of atoms to which the thermostat
                                    ! is linked
     integer         :: start  ! Time window (expressed in md steps) in
     integer         :: stop   ! which the thermostat is defined

     ! Andersen :
     real(kind=DP)   :: freq        ! Collision frequency
     real(kind=DP)   :: mix         ! Mixing coeff for a softer rescaling

     ! Langevin :
     real(kind=DP)   :: damp        ! Langevin damping parameter

     ! Nose-Hoover chain :
     integer         :: nhc_length         ! #-elements in the NH chain
     integer         :: nhc_integ_nstep    ! #-steps to integrate NHC dynamics
     logical         :: nhc_upd            
     real(kind=DP), pointer :: nhc_x(:)    ! Thermostat positions
     real(kind=DP), pointer :: nhc_v(:)    ! Thermostat velocities
     real(kind=DP), pointer :: nhc_g(:)    ! Thermostat accelerations
     real(kind=DP), pointer :: nhc_mass(:) ! Thermostat masses

  end type thermostat

  ! Module variables
  integer, public :: md_thermo_num
  type(thermostat), public, allocatable :: md_thermo(:)

  ! Module parameters
  real(kind=DP), parameter :: aut2fs = 0.0241888468_dp ! atomic unit of time --> fs

  ! List of subroutines
  public  :: md_velocity_verlet_thermostat
  public  :: md_velocity_verlet_thermo1
  public  :: md_velocity_verlet_thermo2
  public  :: md_velocity_verlet_thermo3
  public  :: md_velocity_verlet_thermo4
  public  :: md_select_thermostat
  public  :: md_update_thermostat
  public  :: md_initialise_thermostats
  public  :: md_restore_thermostat
  public  :: md_backup_thermostat
  public  :: md_write_thermostat
  private  :: md_integrate_nhc

  contains


!=============================================================================!
!=============================================================================!

    subroutine md_velocity_verlet_thermostat(pos,time)

      !==============================================================! 
      !                                                              !
      !--------------------------------------------------------------!
      ! Written by Simon M.-M. Dubois (May 2011)                     !
      !==============================================================!

      use comms,           only: pub_on_root
      use constants,       only: dp,stdout
      use md_ions,         only: groups,group_num
      use rundat,          only: md_delta_t,mts_delta_t,mts_xi
      use simulation_cell, only: pub_cell
      use utils,           only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Argument      
      integer, intent(in) :: pos
      real(kind=DP), intent(in) :: time

      ! Variables
      real(kind=DP) :: dt
      integer       :: thermo, group
           
 
      !==============================================================!

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_velocity_verlet_thermostat'
#endif

      if (mts_xi) then
         dt = mts_delta_t
      else
         dt = md_delta_t
      endif

      do group = 1, group_num

         ! Select and initialise the adequate thermostat
         call md_select_thermostat(thermo,group,time)
         if (pub_on_root) write(stdout,'(a,2(i4,2x),f16.8)') &
                 'velocity_verlet_thermostat (group,thermo,time) ', group, thermo, time 
         call md_update_thermostat(thermo,group,time) 


         ! Apply thermostat on group
         if (pos == 1) then 
            call md_velocity_verlet_thermo1(thermo, group, dt)
         elseif (pos == 2) then 
            call md_velocity_verlet_thermo2(thermo, group, dt)
         elseif (pos == 3) then 
            call md_velocity_verlet_thermo3(thermo, group, dt)
         elseif (pos == 4) then 
            call md_velocity_verlet_thermo4(thermo, group, dt)
         endif

      enddo

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_velocity_verlet_thermostat'
#endif

      return


    end subroutine md_velocity_verlet_thermostat

!=============================================================================!
!=============================================================================!

    subroutine md_velocity_verlet_thermo1(thermo,group,dt)

      !==============================================================! 
      !                                                              !
      !--------------------------------------------------------------!
      ! Written by Simon M.-M. Dubois (May 2011)                     !
      !==============================================================!

      use comms,           only: pub_root_node_id, comms_bcast, pub_on_root
      use constants,       only: dp,stdout
      use md_ions,         only: r_t, v_t, a_t, ion_mass, groups, group_num
      use services,        only: services_maxboltzdist
      use simulation_cell, only: pub_cell
      use utils,           only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Argument      
      integer, intent(in) :: thermo
      integer, intent(in) :: group
      real(kind=DP), intent(in) :: dt

      ! Variables
      real(kind=DP) :: dt2
      real(kind=DP) :: scale_v, scale_a
      real(kind=DP) :: scale_random, random_a 
      integer  :: iat, jat, idir 
           
 
      !==============================================================!
     
      ! No thermostat fount 
      if (thermo .le. 0) then

         return 

      ! Thermostat = NONE 
      elseif (md_thermo(thermo)%type == 0) then 
        
         return 

      ! Thermostat = ANDERSEN 
      elseif (md_thermo(thermo)%type == 1) then 

         return   

      ! Thermostat = LANGEVIN 
      elseif (md_thermo(thermo)%type == 2) then 
    
         ! Scale velocities and accelerations according to the damping
         ! factor introduced in Langevin dynamics equation.
         dt2 = 0.5_dp*dt
         scale_v = exp(-md_thermo(thermo)%damp*dt2)
         scale_a = (1-scale_v)/(md_thermo(thermo)%damp*dt2)
         
         !! smmdebug
         !if (pub_on_root) write(stdout,'(a,2(f12.8,2x))') &
         !     'Langevin thermostat (step1) : rescale factors are ', scale_v, scale_a 

         do iat = 1, groups(group)%num_atom
            jat = groups(group)%list_atom(iat)

            !! smmdebug
            !if (pub_on_root) write(stdout,'(a,i,3x,3(f12.8,2x),a)',advance='no') ' v',iat, v_t(:,jat), '-->' 

            v_t(:,jat) = scale_v*v_t(:,jat)

            !! smmdebug
            !if (pub_on_root) write(stdout,'(3(f12.8,2x))') v_t(:,jat)
            !if (pub_on_root) write(stdout,'(a,i,3x,3(f12.8,2x),a)',advance='no') ' a',iat, a_t(:,jat), '-->' 

            a_t(:,jat) = scale_a*a_t(:,jat)

            !! smmdebug
            !if (pub_on_root) write(stdout,'(3(f12.8,2x))') a_t(:,jat)

         enddo


         ! Add random fluctuating accelerations.
         scale_random = md_thermo(thermo)%T*(1-scale_v**2)/dt2**2

         !! smmdebug
         !if (pub_on_root) write(stdout,'(a)') ' Langevin thermostat (step1) : random accelerations...'

         do iat = 1, groups(group)%num_atom
            jat = groups(group)%list_atom(iat)

            !! smmdebug
            !if (pub_on_root) write(stdout,'(a,i,3x,3(f12.8,2x),a)',advance='no') ' a',iat, a_t(:,jat), '-->' 

            do idir = 1,3
               random_a = sqrt(scale_random/ion_mass(jat))*services_maxboltzdist()
               a_t(idir,jat) = a_t(idir,jat) + random_a               
            enddo

            !! smmdebug
            !if (pub_on_root) write(stdout,'(3(f12.8,2x))') a_t(:,jat) 

         enddo 

         ! Make sure all nodes have the same velocity and
         ! acceleration
         call comms_bcast(pub_root_node_id,v_t(:,:))
         call comms_bcast(pub_root_node_id,a_t(:,:))

      ! Thermostat = NOSE-HOOVER 
      elseif (md_thermo(thermo)%type == 3) then 

         !! smmdebug
         !if (pub_on_root) write(stdout,'(a)') &
         !        'Nose-Hoover thermostat (step1) : integrate nhc equation of motion... '

         ! Integrate the dynamics of the Nose-Hoover chain. Rescale the
         ! atomic velocities according to the coupling with the thermostat 
         ! variables
         call md_integrate_nhc(thermo,group,dt)

      endif

      return

    end subroutine md_velocity_verlet_thermo1

!=============================================================================!
!=============================================================================!

    subroutine md_velocity_verlet_thermo2(thermo,group,dt)

      !==============================================================! 
      !                                                              !
      !--------------------------------------------------------------!
      ! Written by Simon M.-M. Dubois (May 2011)                     !
      !==============================================================!

      use comms,           only: pub_root_node_id, pub_on_root, comms_bcast
      use constants,       only: dp, stdout, two_pi
      use md_ions,         only: r_t, v_t, a_t, ion_mass, groups, group_num
      use services,        only: services_maxboltzdist
      use simulation_cell, only: pub_cell
      use utils,           only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Argument      
      integer, intent(in) :: thermo
      integer, intent(in) :: group
      real(kind=DP), intent(in) :: dt

      ! Variables
      real(kind=DP) :: dt2, ran, sigma, pcol, mix1, mix2, vtmp
      integer  :: iat, jat, idir 
      
      !==============================================================!


      ! No thermostat found 
      if (thermo .le. 0) then

         return 

      ! Thermostat = NONE 
      elseif (md_thermo(thermo)%type == 0) then 
        
         return 

      ! Thermostat = ANDERSEN 
      elseif (md_thermo(thermo)%type == 1) then 

         ! Replace the velocities of some atoms with velocities derived 
         ! from the Boltzmann distribution at the thermostat temperature
         dt2 = 0.5_dp*dt
         pcol = 1.0_dp - exp(-md_thermo(thermo)%freq*dt2)
         mix1 = sqrt(1-md_thermo(thermo)%mix**2) 
         mix2 = md_thermo(thermo)%mix
         
         !! smmdebug
         !if (pub_on_root) write(stdout,'(a,2(f12.8,2x))') &
         !        'Andersen thermostat (step2) : pcol, mix ', pcol, mix2

         do iat = 1, groups(group)%num_atom
            jat = groups(group)%list_atom(iat)

            !! smmdebug
            !if (pub_on_root) write(stdout,'(a,i,3x,3(f12.8,2x),a)',advance='no') ' v',iat, v_t(:,jat), '-->' 

            call random_number(ran)

            if (ran.lt.pcol) then
               sigma = sqrt(md_thermo(thermo)%T/ion_mass(jat))
               do idir = 1,3
                  vtmp = sigma*services_maxboltzdist() 
                  v_t(idir,jat) = mix1*v_t(idir,jat)+mix2*vtmp 
               enddo
            endif

            !! smmdebug
            !if (pub_on_root) write(stdout,'(a,f12.8,a,3(f12.8,2x))') '(',ran,')  ',v_t(:,jat)

         enddo

         ! Make sure all nodes have the same velocity
         call comms_bcast(pub_root_node_id,v_t(:,:))

      ! Thermostat = LANGEVIN
      elseif (md_thermo(thermo)%type == 2) then 

         return

      ! Thermostat = NOSE-HOOVER 
      elseif (md_thermo(thermo)%type == 3) then 
 
         return

      endif

    end subroutine md_velocity_verlet_thermo2

!=============================================================================!
!=============================================================================!

    subroutine md_velocity_verlet_thermo3(thermo,group,dt)

      !==============================================================! 
      !                                                              !
      !--------------------------------------------------------------!
      ! Written by Simon M.-M. Dubois (May 2011)                     !
      !==============================================================!

      use comms,           only: pub_root_node_id, comms_bcast, pub_on_root
      use constants,       only: dp,stdout
      use md_ions,         only: r_t, v_t, a_t, ion_mass, groups, group_num
      use services,        only: services_maxboltzdist
      use simulation_cell, only: pub_cell
      use utils,           only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Argument      
      integer, intent(in) :: thermo
      integer, intent(in) :: group
      real(kind=DP), intent(in) :: dt

      ! Variables
      real(kind=DP) :: dt2
      real(kind=DP) :: scale_v, scale_a
      real(kind=DP) :: scale_random, random_a 
      integer  :: iat, jat, idir 
      
      !==============================================================!

      ! No thermostat found 
      if (thermo .le. 0) then

         return 

      ! Thermostat = NONE 
      elseif (md_thermo(thermo)%type == 0) then 
        
         return 

      ! Thermostat = ANDERSEN 
      elseif (md_thermo(thermo)%type == 1) then 

         return   

      ! Thermostat = LANGEVIN 
      elseif (md_thermo(thermo)%type == 2) then 

         ! Scale velocities and accelerations according to the damping
         ! factor introduced in Langevin dynamics equation.
         dt2 = 0.5_dp*dt
         scale_v = exp(-md_thermo(thermo)%damp*dt2)
         scale_a = (1-scale_v)/(md_thermo(thermo)%damp*dt2)
         
         !! smmdebug
         !if (pub_on_root) write(stdout,'(a,2(f12.8,2x))') &
         !     'Langevin thermostat (step3) : rescale factors are ', scale_v, scale_a 

         do iat = 1,groups(group)%num_atom
            jat = groups(group)%list_atom(iat)

            !! smmdebug
            !if (pub_on_root) write(stdout,'(a,i,3x,3(f12.8,2x),a)',advance='no') ' v',iat, v_t(:,jat), '-->' 

            v_t(:,jat) = scale_v*v_t(:,jat)

            !! smmdebug
            !if (pub_on_root) write(stdout,'(3(f12.8,2x))') v_t(:,jat)
            !if (pub_on_root) write(stdout,'(a,i,3x,3(f12.8,2x),a)',advance='no') ' a',iat, a_t(:,jat), '-->' 

            a_t(:,jat) = scale_a*a_t(:,jat)

            !! smmdebug
            !if (pub_on_root) write(stdout,'(3(f12.8,2x))') a_t(:,jat)

         enddo

         ! Add random fluctuating accelerations.
         scale_random = md_thermo(thermo)%T*(1-scale_v**2)/dt2**2

         !! smmdebug
         !if (pub_on_root) write(stdout,'(a)') ' Langevin thermostat (step3) : random accelerations...'

         do iat = 1, groups(group)%num_atom
            jat = groups(group)%list_atom(iat)

            !! smmdebug
            !if (pub_on_root) write(stdout,'(a,i,3x,3(f12.8,2x),a)',advance='no') ' a',iat, a_t(:,jat), '-->' 

            do idir = 1,3
               random_a = sqrt(scale_random/ion_mass(jat))*services_maxboltzdist()
               a_t(idir,jat) = a_t(idir,jat) + random_a
            enddo

            !! smmdebug
            !if (pub_on_root) write(stdout,'(3(f12.8,2x))') a_t(:,jat)

         enddo 
         
         ! Make sure all nodes have the same velocity and
         ! acceleration
         call comms_bcast(pub_root_node_id,v_t(:,:))
         call comms_bcast(pub_root_node_id,a_t(:,:))
         
      ! Thermostat = NOSE-HOOVER
      elseif (md_thermo(thermo)%type == 3) then 

         return

      endif

    end subroutine md_velocity_verlet_thermo3

!==============================================================================!
!==============================================================================!

    subroutine md_velocity_verlet_thermo4(thermo,group,dt)

      !==============================================================! 
      !                                                              !
      !--------------------------------------------------------------!
      ! Written by Simon M.-M. Dubois (May 2011)                     !
      !==============================================================!

      use constants,       only: dp,stdout
      use comms,           only: pub_root_node_id, comms_bcast, pub_on_root
      use md_ions,         only: r_t, v_t, a_t, ion_mass, groups, group_num
      use services,        only: services_maxboltzdist
      use simulation_cell, only: pub_cell
      use utils,           only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Argument      
      integer, intent(in) :: thermo
      integer, intent(in) :: group
      real(kind=DP), intent(in) :: dt

      ! Variables
      real(kind=DP) :: dt2, ran, sigma, pcol, mix1, mix2, vtmp
      integer  :: iat, jat, idir 
      
      !==============================================================!

      ! No thermostat found 
      if (thermo .le. 0) then

         return 

      ! Thermostat = NONE 
      elseif (md_thermo(thermo)%type == 0) then 
        
         return 

      ! Thermostat = ANDERSEN 
      elseif (md_thermo(thermo)%type == 1) then 

         ! Replace the velocities of some atoms with velocities derived 
         ! from the Boltzmann distribution at the thermostat temperature
         dt2 = 0.5_dp*dt
         pcol = 1.0_dp - exp(-md_thermo(thermo)%freq*dt2)
         mix1 = sqrt(1-md_thermo(thermo)%mix**2) 
         mix2 = md_thermo(thermo)%mix
         
         !! smmdebug
         !if (pub_on_root) write(stdout,'(a,2(f12.8,2x))') &
         !        'Andersen thermostat (step2) : pcol, mix ', pcol, mix2

         do iat = 1, groups(group)%num_atom
            jat = groups(group)%list_atom(iat)

            !! smmdebug
            !if (pub_on_root) write(stdout,'(a,i,3x,3(f12.8,2x),a)',advance='no') ' v',iat, v_t(:,jat), '-->' 

            call random_number(ran)

            if (ran.lt.pcol) then
               sigma = sqrt(md_thermo(thermo)%T/ion_mass(jat))
               do idir = 1,3
                  vtmp = sigma*services_maxboltzdist() 
                  v_t(idir,jat) = mix1*v_t(idir,jat)+mix2*vtmp 
               enddo
            endif

            !! smmdebug
            !if (pub_on_root) write(stdout,'(3(f12.8,2x))') v_t(:,jat)

         enddo

         ! Make sure all nodes have the same velocity
         call comms_bcast(pub_root_node_id,v_t(:,:))

      ! Thermostat = LANGEVIN 
      elseif (md_thermo(thermo)%type == 2) then 
   
         return 

      ! Thermostat = NOSE-HOOVER 
      elseif (md_thermo(thermo)%type == 3) then 


         !! smmdebug
         !if (pub_on_root) write(stdout,'(a)') &
         !        'Nose-Hoover thermostat (step4) : integrate nhc equation of motion... '

         ! Integrate the dynamics of the Nose-Hoover chain. Rescale the
         ! atomic velocities according to the coupling wthermo the thermostat 
         ! variables
         call md_integrate_nhc(thermo,group,dt)

      endif

    end subroutine md_velocity_verlet_thermo4

!==============================================================================!
!==============================================================================!

    subroutine md_integrate_nhc(thermo,group,dt)

      !================================================================! 
      ! This subroutine propagates a system of particles coupled       !
      ! to a single Nose-Hoover chain of thermostat according to the   !
      ! NHC part of the extended Liouville evolution operator.         !
      !                                                                !
      ! This is part of the Nose-Hoover approach to NVT extended by    !
      ! Martyna, Tuckermann and Klein to a chain of thermostat.        !
      !                                                                !
      !----------------------------------------------------------------!
      ! Implemented by Simon M.-M. Dubois (May 2011)                   !
      !================================================================!

      use comms,           only: pub_on_root
      use constants,       only: dp,stdout
      use md_ions,         only: v_t, group_num, groups, group_kinetic_energy
      use simulation_cell, only: pub_cell
      use utils,           only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Argument      
      integer, intent(in) :: thermo
      integer, intent(in) :: group
      real(kind=DP), intent(in) :: dt

      ! Variables
      integer    :: nth
      integer    :: nstep
      integer    :: ndof
      real(kind=DP) :: wdt2, wdt4, wdt8
      real(kind=DP) :: akin
      real(kind=DP) :: scale, aa
      integer    :: is, inhc, iat, jat

      !====================================================================!


      ! Initialise NHC local parameters
      nth   = md_thermo(thermo)%nhc_length
      nstep = md_thermo(thermo)%nhc_integ_nstep
      ndof  = 3*groups(group)%num_atom

      scale = 1.0_dp
      wdt2  = dt/(2.0_dp*nstep)
      wdt4  = dt/(4.0_dp*nstep)
      wdt8  = dt/(8.0_dp*nstep)

      ! Get the group kinetic energy
      akin = 2.0_dp*group_kinetic_energy(group)
      
      !! smmdebug
      !if (pub_on_root) write(stdout,'(a,<nth>(f12.8,2x))') ' nhc : mass  ',md_thermo(thermo)%nhc_mass(:) 
      !if (pub_on_root) write(stdout,'(a,<nth>(f12.8,2x))') ' nhc : x_init  ',md_thermo(thermo)%nhc_x(:) 
      !if (pub_on_root) write(stdout,'(a,<nth>(f12.8,2x))') ' nhc : v_init  ',md_thermo(thermo)%nhc_v(:) 
      !if (pub_on_root) write(stdout,'(a,<nth>(f12.8,2x))') ' nhc : g_init  ',md_thermo(thermo)%nhc_g(:) 

      ! Compute force on first thermostat in NHC
      md_thermo(thermo)%nhc_g(1) = (akin - ndof*md_thermo(thermo)%T)/md_thermo(thermo)%nhc_mass(1) 
      
      ! Multiple time step procedure for the propagation 
      do is = 1, md_thermo(thermo)%nhc_integ_nstep
  
         ! Update the thermostat velocities
         md_thermo(thermo)%nhc_v(nth) =  md_thermo(thermo)%nhc_v(nth) + wdt4*md_thermo(thermo)%nhc_g(nth)
         do inhc = 1, nth-1
            aa = exp(-wdt8*md_thermo(thermo)%nhc_v(nth+1-inhc))
            md_thermo(thermo)%nhc_v(nth+1-inhc) = md_thermo(thermo)%nhc_v(nth+1-inhc)*aa*aa &
                  + wdt4*md_thermo(thermo)%nhc_g(nth+1-inhc)*aa
         enddo 

         ! Update the scaling factor for particle velocities
         aa = exp(-wdt2*md_thermo(thermo)%nhc_v(1))
         scale = scale*aa

         ! Update the forces on first thermostat in chain
         md_thermo(thermo)%nhc_g(1) = (scale*scale*akin - ndof * md_thermo(thermo)%T) &
                  / md_thermo(thermo)%nhc_mass(1)

         ! Update the thermostat positions
         do inhc = 1, nth
            md_thermo(thermo)%nhc_x(inhc) = md_thermo(thermo)%nhc_x(inhc) + wdt2 * md_thermo(thermo)%nhc_v(inhc)
         enddo

         ! Update the thermostat velocities and forces
         do inhc = 1, nth-1
            aa = exp(-wdt8*md_thermo(thermo)%nhc_v(inhc+1))
            md_thermo(thermo)%nhc_v(inhc) = md_thermo(thermo)%nhc_v(inhc) * aa * aa &
                  + wdt4 * md_thermo(thermo)%nhc_g(inhc) * aa
            md_thermo(thermo)%nhc_g(inhc) = (md_thermo(thermo)%nhc_mass(inhc) &
                  * md_thermo(thermo)%nhc_v(inhc) * md_thermo(thermo)%nhc_v(inhc) &
                  - md_thermo(thermo)%T) / md_thermo(thermo)%nhc_mass(inhc+1)
         enddo 
         md_thermo(thermo)%nhc_v(nth) =  md_thermo(thermo)%nhc_v(nth) & 
                  + wdt4*md_thermo(thermo)%nhc_g(nth)
      enddo

      ! Update the particle velocities according to the scaling factor
      do iat = 1, groups(group)%num_atom
         jat = groups(group)%list_atom(iat)

         !! smmdebug
         !if (pub_on_root) write(stdout,'(a,i,3x,3(f12.8,2x),a)',advance='no') ' v',iat, v_t(:,jat), '-->' 

         v_t(:,jat) = scale  * v_t(:,jat)

         !! smmdebug
         !if (pub_on_root) write(stdout,'(3(f12.8,2x))') v_t(:,jat)

      enddo
      
      !! smmdebug
      !if (pub_on_root) write(stdout,'(a,<nth>(f12.8,2x))') ' nhc : x_end  ',md_thermo(thermo)%nhc_x(:) 
      !if (pub_on_root) write(stdout,'(a,<nth>(f12.8,2x))') ' nhc : v_end  ',md_thermo(thermo)%nhc_v(:) 
      !if (pub_on_root) write(stdout,'(a,<nth>(f12.8,2x))') ' nhc : g_end  ',md_thermo(thermo)%nhc_g(:) 

      ! Stability check
      if (maxval(abs(md_thermo(thermo)%nhc_v)).gt. 1.0E8_dp) then
         if (pub_on_root) then
            write(stdout,*) 'Warning : dynamic of Nose-Hoover chain going unstable '
            write(stdout,*) 'Action  : increase the thermostat masses  ' 
            write(stdout,*) '   or/and decrease the step length for NHC integration ' 
            write(stdout,*) '   or/and decrease the setp lenght for MD integration. '
         endif
      endif
 
      return
      
    end subroutine md_integrate_nhc
    
       
!==============================================================================!
!==============================================================================!

    subroutine md_initialise_thermostats()

      !================================================================! 
      !                                                                !
      !----------------------------------------------------------------!
      ! Implemented by Simon M.-M. Dubois (August 2011)                   !
      !================================================================!

      use constants,       only: stdout, dp
      use utils,           only: utils_alloc_check

      implicit none

      ! Local variables
      integer :: thermo, nhcl
      integer :: ierr
      
      !====================================================================!

      ! Loop over all thermostats
      thermo_loop : do thermo = 1, md_thermo_num

         md_thermo(thermo)%T = md_thermo(thermo)%Tinit

         if (md_thermo(thermo)%type == 3) then

            nhcl = md_thermo(thermo)%nhc_length

            ! Allocate dynamical variables
            allocate(md_thermo(thermo)%nhc_x(nhcl),stat=ierr)
            call utils_alloc_check('md_initialise_thermostats','nhc_x',ierr)
            allocate(md_thermo(thermo)%nhc_v(nhcl),stat=ierr)
            call utils_alloc_check('md_initialise_thermostats','nhc_v',ierr)
            allocate(md_thermo(thermo)%nhc_g(nhcl),stat=ierr)
            call utils_alloc_check('md_initialise_thermostats','nhc_g',ierr)
            allocate(md_thermo(thermo)%nhc_mass(nhcl),stat=ierr)
            call utils_alloc_check('md_initialise_thermostats','nhc_mass',ierr)

            ! Initialise dynamical variables
            md_thermo(thermo)%nhc_x = 0.0_dp
            md_thermo(thermo)%nhc_v = 0.0_dp
            md_thermo(thermo)%nhc_g = 0.0_dp
            md_thermo(thermo)%nhc_mass = 0.0_dp

         endif 

      enddo thermo_loop

      return

    end subroutine md_initialise_thermostats
 
!==============================================================================!
!==============================================================================!

    subroutine md_update_thermostat(thermo,group,time)

      !================================================================! 
      !                                                                !
      !----------------------------------------------------------------!
      ! Implemented by Simon M.-M. Dubois (May 2011)                   !
      !================================================================!

      use constants,       only: stdout, two_pi
      use comms,           only: pub_on_root
      use md_ions,         only: groups
      use rundat,          only: md_delta_t, mts_delta_t
      use utils,           only: utils_alloc_check

      implicit none

      ! Argument      
      integer, intent(in)        :: thermo
      integer, intent(in)        :: group
      real(kind=DP), intent(in)  :: time

      ! Local variables
      real(kind=DP)  :: time_start
      integer :: ndof
      
      !====================================================================!

      if (thermo .gt. 0) then

         time_start = (md_thermo(thermo)%start-1) * md_delta_t

         ! Update the thermostat temperature  
         md_thermo(thermo)%T = md_thermo(thermo)%Tinit + &
            (time - time_start) * md_thermo(thermo)%Tgrad / md_delta_t
         
         if (pub_on_root) write(stdout,'(a,2(f8.4,2x))') &
                 'md_update_thermostat (T,Tinit) ', md_thermo(thermo)%T, md_thermo(thermo)%Tinit 

         ! If required, update the NHC masses parameters  
         if (md_thermo(thermo)%type ==  3) then 

            ndof = groups(group)%num_atom*3
            if (md_thermo(thermo)%nhc_upd) then
               md_thermo(thermo)%nhc_mass(:) = md_thermo(thermo)%T/(two_pi*md_thermo(thermo)%freq)**2
            else
               md_thermo(thermo)%nhc_mass(:) = md_thermo(thermo)%Tinit/(two_pi*md_thermo(thermo)%freq)**2
            endif
            md_thermo(thermo)%nhc_mass(1) = ndof*md_thermo(thermo)%nhc_mass(1)
         endif

      endif

      return

    end subroutine md_update_thermostat
 
!==============================================================================!
!==============================================================================!

 
    subroutine md_select_thermostat(thermo,group,time)
   
      !================================================================! 
      !                                                                !
      !----------------------------------------------------------------!
      ! Implemented by Simon M.-M. Dubois (August 2011)                !
      !================================================================!

      use comms,           only: pub_on_root
      use md_ions,         only: groups 
      use rundat,          only: md_delta_t
    
      implicit none
    
      ! Argument      
      integer, intent(in)             :: group
      integer, intent(out)            :: thermo
      real(kind=DP), intent(in)       :: time
    
      ! Local variables
      real(kind=DP) :: time_start, time_stop
      integer       :: ith
      logical       :: found
      
      found = .false.
      
      ! Look for the thermostat corresponding to group
      do ith = 1, md_thermo_num
         time_start = (md_thermo(ith)%start-1) * md_delta_t + tiny(1.0_dp)
         time_stop = (md_thermo(ith)%stop + 0.5_dp) * md_delta_t
         if ((md_thermo(ith)%group .eq. groups(group)%label .or. &
             md_thermo(ith)%group .eq. 0 ) .and. time_start .lt. time .and. & 
             time_stop .ge. time) then
            thermo = ith
            found = .true.
            exit
         endif
      enddo
    
      ! If not found, report an error message
      ! and use Newton's equation of motion
      if (.not. found) then
         thermo = -1
         if (pub_on_root) then
            write(stdout,'(/a)') 'Error in md_select_thermostat: '
            write(stdout,'(a,i4.4,a,f16.8,a)') 'No suitable thermostat found for group ', &
                  groups(group)%label, ' at time ', time, ' aut !' 
         endif
      endif

    end subroutine md_select_thermostat
 
!==============================================================================!
!==============================================================================!

    subroutine md_restore_thermostat()

      use constants,       only: stdout,dp
      use comms,           only: pub_on_root,pub_root_node_id, &
             comms_abort, comms_bcast
      use rundat,          only: pub_rootname
      use simulation_cell, only: pub_cell
      use utils,           only: utils_unit, utils_close_unit_check, &
             utils_alloc_check, utils_dealloc_check

      implicit none

      ! Local variables
      real(kind=dp), allocatable :: dbuf(:)   
      integer           :: dbuf_length, nhcl
      integer           :: io_unit, io_stat
      integer           :: ith, il, idx
      integer           :: ierr
      character(len=80) :: filename


#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_restore_thermostat'
#endif

      ! Compute buffer length
      dbuf_length = 0
      do ith = 1, md_thermo_num
         dbuf_length = dbuf_length + 1
         if (md_thermo(ith)%nhc_length .ge. 1) then
            dbuf_length = dbuf_length + 4*md_thermo(ith)%nhc_length
         endif
      enddo

      ! Allocate buffer
      allocate(dbuf(dbuf_length),stat=ierr)
      call utils_alloc_check('md_restore_thermostat','dbuf',ierr)

      if (pub_on_root) then       

         ! Find available unit specifier
         io_unit = utils_unit()
         write(filename,'(a,a)') trim(pub_rootname),'.thermo.restart'

         ! Open .thermo file
         open(unit=io_unit,iostat=io_stat,file=filename,&
              access='SEQUENTIAL',form='UNFORMATTED',position='REWIND', &
              action='READ')

         ! Read data in buffer 
         if (io_stat == 0) then 

            idx = 0
            thermo_loop1 : do ith = 1, md_thermo_num

               idx = idx + 1
               read(io_unit) dbuf(idx) 

               nhcl = md_thermo(ith)%nhc_length
               if (md_thermo(ith)%type ==3 .and. nhcl .ge. 1) then
                  do il = 1, 4
                     read(io_unit) dbuf(idx+1:idx+nhcl) 
                     idx = idx + nhcl 
                  enddo
               endif

            enddo thermo_loop1

            close(unit=io_unit,iostat=io_stat)
            call utils_close_unit_check('md_restore_thermostat',filename,io_stat)

         else
            write(stdout,'(3a,/,a)') 'WARNING : md_restore_thermostat failed to open "', &
              trim(filename),'"','Initialise thermostat dynamic variables to default...'

         endif
      endif
     
      ! Broadcast buffer to all nodes 
      call comms_bcast(pub_root_node_id,io_stat)
      call comms_bcast(pub_root_node_id,dbuf,dbuf_length)

      ! Read buffer into thermostat arrays
      if (io_stat == 0) then

         idx = 0
         thermo_loop2 : do ith = 1, md_thermo_num

            idx = idx + 1 
            md_thermo(ith)%T = dbuf(idx) 

            nhcl = md_thermo(ith)%nhc_length
            if (nhcl .ge. 1) then
               md_thermo(ith)%nhc_x(:) = dbuf(idx+1:idx+nhcl)  
               md_thermo(ith)%nhc_v(:) = dbuf(idx+nhcl+1:idx+2*nhcl)  
               md_thermo(ith)%nhc_g(:) = dbuf(idx+2*nhcl+1:idx+3*nhcl)  
               md_thermo(ith)%nhc_mass(:) = dbuf(idx+3*nhcl+1:idx+4*nhcl)  
               idx = idx + 4*nhcl
            endif

         enddo thermo_loop2
      else

         thermo_loop3 : do ith = 1, md_thermo_num

            md_thermo(ith)%T = md_thermo(ith)%Tinit

            nhcl = md_thermo(ith)%nhc_length
            if (nhcl .ge. 1) then
               md_thermo(ith)%nhc_x(:) = 0.0_dp
               md_thermo(ith)%nhc_v(:) = 0.0_dp
               md_thermo(ith)%nhc_g(:) = 0.0_dp
               md_thermo(ith)%nhc_mass(:) = 0.0_dp
            endif

         enddo thermo_loop3
      endif

      ! Allocate buffer
      deallocate(dbuf,stat=ierr)
      call utils_dealloc_check('md_restore_thermostat','dbuf',ierr)

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_restore_thermostat'
#endif

      return

    end subroutine md_restore_thermostat

!==============================================================================!
!==============================================================================!

    subroutine md_backup_thermostat()

      use constants,       only: dp
      use comms,           only: pub_on_root
      use rundat,          only: pub_rootname
      use utils,           only: utils_unit,utils_open_unit_check, &
            utils_close_unit_check
      
      implicit none
    
      ! Local variables
      character(len=80) :: filename
      integer           :: fileunit, iostatus
      integer           :: ith

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_backup_thermostat'
#endif

      if (pub_on_root) then       
         fileunit = utils_unit()
         write(filename,'(a,a)') trim(pub_rootname),'.thermo.restart'
         open(unit=fileunit,file=filename,iostat=iostatus,status='UNKNOWN', &
                form='UNFORMATTED',position='APPEND',action='WRITE')
         call utils_open_unit_check('md_backup_thermostat',filename,iostatus)
         
         thermo_loop : do ith = 1, md_thermo_num
            write(fileunit) md_thermo(ith)%T
            if (md_thermo(ith)%nhc_length .ge. 1) then
               write(fileunit) md_thermo(ith)%nhc_x(:)
               write(fileunit) md_thermo(ith)%nhc_v(:)
               write(fileunit) md_thermo(ith)%nhc_g(:)
               write(fileunit) md_thermo(ith)%nhc_mass(:)
            endif
         enddo thermo_loop

         close(unit=fileunit,iostat=iostatus)
         call utils_close_unit_check('md_backup_thermostat',filename,iostatus)

      endif

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_backup_thermostat'
#endif

      return

    end subroutine md_backup_thermostat

!==============================================================================!
!==============================================================================!

    subroutine md_write_thermostat(time)

      use constants,       only: dp
      use comms,           only: pub_on_root
      use rundat,          only: pub_rootname
      use utils,           only: utils_unit,utils_open_unit_check, &
            utils_close_unit_check
      
      implicit none
    
      ! Argument
      real(kind=dp)     :: time

      ! Local variables
      character(len=80) :: filename
      character(len=12) :: tkind(4)
      character(len=80) :: position, stat
      character(len=5)  :: tmp_length
      character(len=80) :: nhc_format
      integer           :: fileunit, iostatus
      integer           :: ith

      data tkind /'None','Andersen','Langevin','Nose-Hoover'/

      if (time == 0.0_dp) then
         position = 'REWIND'
         stat = 'REPLACE'
      else
         position = 'APPEND'
         stat = 'UNKNOWN'
      endif 

      if (pub_on_root) then       
         fileunit = utils_unit()
         write(filename,'(a,a)') trim(pub_rootname),'.thermo'
         open(unit=fileunit,file=filename,iostat=iostatus,status=stat, &
                form='FORMATTED',position=position,action='WRITE')
         call utils_open_unit_check('md_write_thermostat',filename,iostatus)
         
         write(fileunit,1) time*aut2fs, 'fs'

         thermo_loop : do ith = 1, md_thermo_num
            write(fileunit,2) '       ', ith
            write(fileunit,3) '       ', tkind(md_thermo(ith)%type + 1)
            write(fileunit,4) 'T     :', md_thermo(ith)%T
            write(fileunit,6) 'Start :', md_thermo(ith)%start, 'Stop  :',md_thermo(ith)%stop
            write(fileunit,5) 'Tinit :', md_thermo(ith)%Tinit, 'Tgrad :',md_thermo(ith)%Tgrad
            write(fileunit,5) 'Freq  :', md_thermo(ith)%freq,  'Damp  :',md_thermo(ith)%damp 
            write(fileunit,6) 'NHC   :', md_thermo(ith)%nhc_length, 'Int   :',md_thermo(ith)%nhc_integ_nstep
            if (md_thermo(ith)%type == 3 .and. md_thermo(ith)%nhc_length .ge. 1) then
               write(tmp_length,'(i5)') md_thermo(ith)%nhc_length
               write(nhc_format,'(a,i4,a)') '(9x,',trim(adjustl(tmp_length)),'(3x,es18.8e3))'
               write(fileunit,nhc_format) md_thermo(ith)%nhc_x(:)
               write(fileunit,nhc_format) md_thermo(ith)%nhc_v(:)
               write(fileunit,nhc_format) md_thermo(ith)%nhc_g(:)
               write(fileunit,nhc_format) md_thermo(ith)%nhc_mass(:)
            endif
         enddo thermo_loop
         write(fileunit,7) 

         close(unit=fileunit,iostat=iostatus)
         call utils_close_unit_check('md_write_thermostat',filename,iostatus)

      endif

1     format(/,12x,es18.8e3,2x,a)
2     format(/,4x,a,1x,i4) 
3     format(4x,a,4x,a12) 
4     format(4x,a,1x,es18.8e3)
5     format(4x,2(a,1x,es18.8e3,3x))
6     format(4x,2(a,1x,i4,17x))
7     format(/,58('-'))

      return

    end subroutine md_write_thermostat

!==============================================================================!
!==============================================================================!




 
end module md_thermostat
    
    
    
    
    
