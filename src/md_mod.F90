! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-


!============================================================================!
!                                                                            !
!                        Molecular Dynamics Module                           !
!                                                                            !
!----------------------------------------------------------------------------!
! Written by by Simon M.-M. Dubois (May 2010) based on the original          !
! module created by Arash A. Mostofi (May 2006)                              !
!============================================================================!

module md

  use constants, only: dp, stdout

  implicit none

  private

  ! MD status 
  real(kind=dp), save :: md_time
  real(kind=dp), save :: md_energy
  real(kind=dp), save :: md_kinetic_energy
  real(kind=dp), save :: md_potential_energy  
  real(kind=dp), save :: md_actual_temp

  ! MTS integration parameters
  real(kind=dp), save :: mts_ngwfs_threshold
  real(kind=dp), save :: mts_lnv_threshold
  real(kind=dp), save :: ref_ngwfs_threshold
  real(kind=dp), save :: ref_lnv_threshold

  ! Module parameters
  real(kind=DP), parameter :: aut2fs = 0.0241888468_dp ! atomic unit of time --> fs
  real(kind=DP), parameter :: ha2k   = 3.1577462e5_dp  ! hartree --> kelvin
  character (len=30), parameter :: time_label = 'fs'
  character (len=30), parameter :: energy_label = 'Ha'
  character (len=30), parameter :: temperature_label = 'K'

  ! List of subroutines and functions :
  public  :: md_main

  contains

!==============================================================================!
!==============================================================================!

    subroutine md_main(total_energy,forces,elements)

      use constants,       only: dp,stdout
      use comms,           only: pub_on_root
      use ion,             only: element
      use md_ions,         only: md_allocate_struct, md_deallocate_struct
      use simulation_cell, only: pub_cell

      implicit none

      ! Arguments
      type(ELEMENT),    intent(inout) :: elements(pub_cell%nat)
      real(kind=DP),    intent(inout) :: total_energy
      real(kind=DP),    intent(inout) :: forces(1:3,pub_cell%nat)

      ! Internal variables

      ! Banner
      if (pub_on_root) write(stdout,1) ' Starting ONETEP Molecular Dynamics '

      ! Allocate module arrays
      call md_allocate_struct(elements)

      ! Velocity verlet integration algorithm
      call md_velocity_verlet(total_energy,forces,elements)

      ! Deallocate module arrays
      call md_deallocate_struct()

      if (pub_on_root) write(stdout,1) '  End of ONETEP Molecular Dynamics  '

1     format(/,80('x'),/,11('x'),a47,11x,11('x'),/,80('x'),/)

      return

    end subroutine md_main

!==============================================================================!
!==============================================================================!

    subroutine md_velocity_verlet(total_energy,forces,elements)

    !====================================================================! 
    ! Velocity Verlet algorithm :                                        !
    ! Integrate the equation of motion for a system of classical         !
    ! particles, if required coupled with one of the following           ! 
    ! thermostat (ANDERSON/LANGEVIN/NOSE-HOOVER/NOSE-HOOVER-CHAIN)       !
    !                                                                    !
    ! Calls to thermostat subroutine are interleaved with the usual      !
    ! velocity verlet algorythm:                                         !
    !                                                                    !
    ! thermo (step1)                                                     !
    !    v(t+dt/2) = v(t) + a(t)*dt/2                                    !
    ! thermo (step2)                                                     !
    !    r(t+dt) = r(t) + v(t+dt/2)*dt                                   !
    !    call energy_and_force_calculate()                               !
    ! thermo (step3)                                                     !
    !    v(t+dt) = v(t+dt/2) + a(t+dt)*dt/2                              !
    ! thermo (step4)                                                     !
    !                                                                    !
    !--------------------------------------------------------------------!
    ! Implemented by Simon M.-M. Dubois (May 2011)                       !
    !====================================================================!

      use comms,            only: comms_abort,pub_on_root,pub_my_node_id
      use constants,        only: dp,stdout
      use ion,              only: element
      use md_ions,          only: a_t, v_t, r_t, ion_mass, &
           groups, group_num, md_initialise_groups, md_initialise_velocities, &
           md_initialise_positions, md_initialise_ion_masses, &
           kinetic_energy
      use md_thermostat,    only: md_select_thermostat, md_thermo, &
           md_restore_thermostat, md_backup_thermostat, md_write_thermostat, &
           md_velocity_verlet_thermostat, md_initialise_thermostats
      use restart,          only: restart_store_new_item, &
           restart_store_create, restart_store_destroy, restart_store_reset
      use rundat,           only: md_num_iter,md_delta_t, mts_delta_t, &
           pub_write_positions, md_init_velocities, mix_dkn_num, &
           mix_ngwfs_num, mts_nstep, mts_xi, ngwf_threshold_orig, &
           lnv_threshold_orig, mts_ngwf_threshold, mts_lnv_threshold, &
           pub_md_properties, md_restart, md_reset_dkn_ngwfs
      use simulation_cell,  only: pub_cell
      use utils,            only: utils_alloc_check, utils_dealloc_check, &
           utils_unit

      implicit none

      ! Argument      
      type(ELEMENT),    intent(inout) :: elements(pub_cell%nat)
      real(kind=DP),    intent(inout) :: total_energy
      real(kind=DP),    intent(inout) :: forces(1:3,pub_cell%nat)

      ! Internal variables
      real(kind=DP) :: tgroup
      integer :: md_iter, mts_step
      integer :: ndof
      integer :: iat, ig, ith

    !====================================================================!

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_velocity_verlet'
#endif

      ! Initialise MD
      md_time = 0.0_dp; ndof = 3*pub_cell%nat

      ! Initialise multi-time-step integration scheme
      if (pub_on_root) write(stdout,'(a)',advance='no') &
          ' Initialising multi time-step (MTS) integration scheme...'
      call init_mts_scheme()
      if (pub_on_root) write(stdout,*) 'done'

      ! Initialise ion masses
      if (pub_on_root) write(stdout,'(a)',advance='no') &
          ' Initialising ion masses...'
      call md_initialise_ion_masses(elements)
      if (pub_on_root) write(stdout,*) 'done'

      ! Initialise groups
      if (pub_on_root) write(stdout,'(a)',advance='no') &
          ' Initialising ion groups...'
      call md_initialise_groups(elements)
      if (pub_on_root) write(stdout,*) 'done'

      if (md_restart) then 
         ! Read positions and velocities
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             ' MD structure from previous calculation...'
         call md_restore_trajectory(elements)
         if (pub_on_root) write(stdout,*) 'done'
      endif
      
      if (md_restart) then 
         ! Read thermostat parameters 
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             ' MD thermostats from previous calculation...'
         call md_restore_thermostat()
         if (pub_on_root) write(stdout,*) 'done'
      endif
      
      if (.not. md_restart) then
         ! Initialise positions 
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             ' Initialising MD structure...'
         call md_initialise_positions(elements)
         if (pub_on_root) write(stdout,*) 'done'
         
         ! Initialise velocities 
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             ' Initialising MD velocities ...'
         do ig = 1, group_num
            call md_select_thermostat(ith,ig,mts_delta_t)
            tgroup = md_thermo(ith)%T 
            call md_initialise_velocities(ig, tgroup, elements) 
         enddo
         if (pub_on_root) write(stdout,*) 'done'

         ! Initialise positions 
         if (pub_on_root) write(stdout,'(a)',advance='no') &
             ' Initialising MD thermostats...'
         call md_initialise_thermostats()
         if (pub_on_root) write(stdout,*) 'done'
         
      endif

      ! Write initial coordinates
      if (pub_write_positions) call md_write_coordinates(elements)

      ! Store the particle positions (required for extrapolation 
      ! of NGWFS/denskern at subsequent MD steps)
      if (mix_dkn_num .gt. 0  .or. mix_ngwfs_num .gt. 0) then
         call restart_store_create(max(mix_ngwfs_num, mix_dkn_num)+1)
         call restart_store_new_item(r_t)
      endif 

      ! Calculate energy and forces 
      call md_energy_and_force(total_energy,forces,elements,'normal')
     
      ! Calculate acceleration
      if (pub_on_root) write(stdout,*) 'Calculate accelerations... '
      do iat=1,pub_cell%nat
         a_t(:,iat) = forces(:,iat)/ion_mass(iat)
      enddo

      ! Calculate energies and temperature
      if (pub_on_root) write(stdout,*) 'Calculate energies and temperature... '
      md_potential_energy = total_energy
      md_kinetic_energy = kinetic_energy()
      md_energy = md_kinetic_energy+md_potential_energy
      md_actual_temp = 2.0_dp*md_kinetic_energy/ndof

      ! Output MD info
      call md_write_md_info()

      ! Write initial state
      if (.not. md_restart) then
         if (pub_on_root) write(stdout,'(a)',advance='no') ' Writing initial state...'
         call md_write_trajectory(elements)
         call md_write_thermostat(md_time)
         if (pub_on_root) write(stdout,*) 'done'
      endif

      ! Write backup files
      if (pub_on_root) write(stdout,'(a)',advance='no') ' Writing MD restart files...'
      call md_backup_trajectory()        
      call md_backup_thermostat()        
      if (pub_on_root) write(stdout,*) 'done'
           
      md_loop : do md_iter=1,md_num_iter

         if (pub_on_root) then
            write(stdout,1) 'Starting MD iteration :', md_iter
         endif 

         mts_loop : do mts_step = 1, mts_nstep

            if (pub_on_root .and. mts_nstep .gt. 1) then
               write(stdout,2) ' Multi-time step integration :', mts_step
            endif

            md_time = md_time + mts_delta_t
         
            ! Scale particle velocities, 
            ! update thermostat velocities and positions
            call md_velocity_verlet_thermostat(1,md_time)
          
            ! v(t+dt/2) = v(t) + 0.5 * dt * a(t)
            do iat=1,pub_cell%nat
               v_t(:,iat) = v_t(:,iat) + 0.5_dp * mts_delta_t * a_t(:,iat)
            enddo
            
            ! Scale particle velocities, 
            ! update thermostat velocities and positions
            call md_velocity_verlet_thermostat(2,md_time)

            ! r(t+dt) = r(t) + dt * v(t+dt/2)
            do iat=1,pub_cell%nat
               r_t(:,iat) = r_t(:,iat) + md_delta_t * v_t(:,iat)
            enddo
         
            ! Apply constraints
            ! call shake() 
         
            ! Store particle coordinates (required for extrapolation 
            ! of NGWFS/denskern at subsequent MD steps)
            if (mix_dkn_num .gt. 0 .or. mix_ngwfs_num .gt. 0) then
               
               ! If required, reset the extrapolation scheme
               if (mts_step == 1 .and. &
                   modulo(md_iter,md_reset_dkn_ngwfs) == 0) then
                  call restart_store_reset()
               endif

               call restart_store_new_item(r_t)

            endif 
         
            ! Update elements 
            call md_update_coordinates(elements)
         
            ! Write coordinates
            if (pub_write_positions) call md_write_coordinates(elements)
        
            ! Calculate energy and forces 
            if (mts_nstep.gt.1 .and. mts_step == mts_nstep) then
               call md_energy_and_force(total_energy,forces,elements,'high')
            elseif (mts_nstep.gt.1) then
               call md_energy_and_force(total_energy,forces,elements,'low')
            else 
               call md_energy_and_force(total_energy,forces,elements,'normal')
            endif

            ! Calculate acceleration
            do iat=1,pub_cell%nat
               a_t(:,iat) = forces(:,iat)/ion_mass(iat)
            enddo
            
            ! Scale particle velocities, 
            ! update thermostat velocities and positions
            call md_velocity_verlet_thermostat(3,md_time)

            ! v(t+dt) = v(t+dt/2) + 0.5 * dt * a(t+dt)
            do iat=1,pub_cell%nat
               v_t(:,iat) = v_t(:,iat) + 0.5_dp * mts_delta_t * a_t(:,iat)
            enddo
         
            ! Apply constraints
            ! call shake() 
         
            ! Scale particle velocities, update thermostat velocities and positions
            call md_velocity_verlet_thermostat(4,md_time)
         
            ! Energies and temperature
            md_potential_energy = total_energy
            md_kinetic_energy = kinetic_energy()
            md_energy = md_kinetic_energy + md_potential_energy
            md_actual_temp = 2.0_dp * md_kinetic_energy/ndof
         
            ! Output MD info
            call md_write_md_info()
         
            ! Write trajectory
            if (pub_on_root) write(stdout,'(a)',advance='no') ' Writing MD trajectory...'
            call md_write_trajectory(elements)        
            call md_write_thermostat(md_time)        
            if (pub_on_root) write(stdout,*) 'done'
           
            ! Write backup files
            if (pub_on_root) write(stdout,'(a)',advance='no') ' Writing MD restart files...'
            call md_backup_trajectory()        
            call md_backup_thermostat()        
            if (pub_on_root) write(stdout,*) 'done'
           
         enddo mts_loop 

      enddo md_loop


      ! Deallocate extrapolation arrays
      if (mix_dkn_num .gt. 0  .or. mix_ngwfs_num .gt. 0) then
         call restart_store_destroy()
      endif 

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_velocity_verlet'
#endif
      return

1     format(/,80('x'),/,11('x'),a37,1x,i9,11x,11('x'),/,80('x'),/) 
2     format(/,10x,60('x'),/,10x,'x',a37,1x,i9,11x,'x',/,10x,60('x'),/) 

    end subroutine md_velocity_verlet

!==============================================================================!
!==============================================================================!

    subroutine md_energy_and_force(total_energy,forces,elements,mode)
 
    !====================================================================! 
    ! Buffer layer between MD integration subroutines and                !
    ! energy_and_force_calculate subroutine.                             !
    !                                                                    !
    !--------------------------------------------------------------------!
    ! Implemented by Simon M.-M. Dubois (May 2011)                       !
    !====================================================================!

      use comms,            only: comms_abort,comms_barrier,pub_on_root
      use constants,        only: dp,stdout
      use energy_and_force, only: energy_and_force_calculate
      use ion,              only: ELEMENT
      use simulation_cell,  only: pub_cell
      use rundat,           only: mts_nstep
      use utils,            only: utils_alloc_check, utils_dealloc_check

      implicit none

      ! Argumentas
      type(ELEMENT), intent(inout)   :: elements(pub_cell%nat)
      real(kind=DP), intent(out)     :: total_energy
      real(kind=DP), intent(out)     :: forces(1:3,1:pub_cell%nat)
      character(len=*), intent(in)   :: mode
     
      ! Variables
      real(kind=DP), allocatable :: forces_mts(:,:)
      real(kind=DP) :: fc(3)
      integer       :: iat, ierr

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_energy_and_force'
#endif


      ! Compute energy and forces
      
      if (mode == 'normal') then

         call comms_barrier()
         call energy_and_force_calculate(total_energy,forces,elements)

      ! Turn to mts convergence parameters for the low-accuracy
      ! mts steps
      elseif (mode == 'low') then

         call comms_barrier()
         call activate_mts_param()
         call energy_and_force_calculate(total_energy,forces,elements)
         call deactivate_mts_param()

      ! mts apply corrections to the forces
      elseif (mode == 'high') then

         ! Start with a low-accuracy estimation of the forces
         allocate(forces_mts(3,pub_cell%nat),stat=ierr)
         call utils_alloc_check('md_energy_and_force','forces_mts',ierr)

         call comms_barrier()
         call activate_mts_param()
         call energy_and_force_calculate(total_energy,forces_mts,elements)
         call deactivate_mts_param()

         ! Turn on the high accuracy parameters 
         call energy_and_force_calculate(total_energy,forces,elements)
         
         ! Apply corrections to forces
         do iat = 1,pub_cell%nat
            fc(:) = real(mts_nstep)*(forces(:,iat) - forces_mts(:,iat))
            forces(:,iat) = forces_mts(:,iat) + fc(:)
         enddo
         
         deallocate(forces_mts,stat=ierr)
         call utils_dealloc_check('md_energy_and_force','forces_mts',ierr)

      endif
 
#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_energy_and_force'
#endif

      return 

    end subroutine md_energy_and_force

!==============================================================================!
!==============================================================================!


    subroutine md_update_coordinates(elements)

      use comms,           only: pub_on_root
      use constants,       only: dp
      use ion,             only: element
      use md_ions,         only: r_t, v_t
      use simulation_cell, only: pub_cell
      use services,        only: services_rationalise_coords

      implicit none

      ! Arguments
      type(element), intent(inout) :: elements(pub_cell%nat)

      ! Local variables
      integer :: iat
     
#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_update_coordinates'
#endif

      call services_rationalise_coords(pub_cell%nat,r_t)
      do iat=1,pub_cell%nat
         elements(iat)%centre%x =  r_t(1,iat)
         elements(iat)%centre%y =  r_t(2,iat)
         elements(iat)%centre%z =  r_t(3,iat)
         elements(iat)%ion_velocity(:) =  v_t(:,iat)
      enddo

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_update_coordinates'
#endif

      return

    end subroutine md_update_coordinates


!==============================================================================!
!==============================================================================!


    subroutine md_write_coordinates(elements)


      use comms,           only: pub_on_root
      use constants,       only: stdout
      use ion,             only: element
      use simulation_cell, only: pub_cell

      implicit none

      type(element), intent(in) :: elements(pub_cell%nat)

      ! Local variables
      integer :: iat

      if (pub_on_root) then
         write(stdout,*)
         write(stdout,*)'                           -------------------------------'
         write(stdout,*)'                                     Cell Contents'
         write(stdout,*)'                           -------------------------------'
         write(stdout,*)         
         write(stdout,21)
         write(stdout,*)'           x  Element    Atom     Absolute co-ordinates of atoms (bohr) x'
         write(stdout,*)'           x            Number         x           y           z        x'
         write(stdout,22)
         do iat=1,pub_cell%nat
            write(stdout,20) elements(iat)%symbol, iat, elements(iat)%centre%x, &
                             elements(iat)%centre%y, elements(iat)%centre%z
         end do
         write(stdout,21)
         write(stdout,*)
      endif

20    format(11x,' x',a6,1x,i8,5x,3f12.6,'    x')
21    format(1x,'           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
22    format(1x,'           x------------------------------------------------------------x')

      return

    end subroutine md_write_coordinates


!==============================================================================!
!==============================================================================!


    subroutine md_write_md_info()


      use comms,     only: pub_on_root
      use constants, only: stdout

      implicit none

      if (pub_on_root) then
         write(stdout,*)
         write(stdout,21)

         write(stdout,1)'MD INFO                        '
         write(stdout,1)' '
         write(stdout,2)'            time:',md_time*aut2fs,trim(time_label)
         write(stdout,1)' '
         write(stdout,2)'Potential Energy:',md_potential_energy,trim(energy_label)
         write(stdout,2)'Kinetic   Energy:',md_kinetic_energy,trim(energy_label)
         write(stdout,2)'Total     Energy:',md_energy,trim(energy_label)
         write(stdout,1)' '
         write(stdout,2)'     Temperature:',md_actual_temp*ha2k,trim(temperature_label)
         write(stdout,1)' '
         write(stdout,21)
         write(stdout,*)
      endif

1     format(11('x'),a58                ,11('x'))
2     format(11('x'),a20,f14.6,1x,a20,3x,11('x'))
21    format(80('x')) 

      return

    end subroutine md_write_md_info



!==============================================================================!
!==============================================================================!

    subroutine md_restore_trajectory(elements)


      use constants,       only: stdout,dp
      use comms,           only: pub_on_root,pub_root_node_id, &
             comms_abort, comms_bcast
      use ion,             only: element
      use md_ions,         only: r_t, v_t
      use rundat,          only: pub_rootname, md_restart
      use simulation_cell, only: pub_cell
      use utils,           only: utils_unit,utils_open_unit_check, &
             utils_close_unit_check
      
      implicit none

      ! Argument 
      type(element), intent(inout) :: elements(pub_cell%nat)

      ! Local variables
      integer           :: io_unit,io_stat
      integer           :: idir,iat
      real(kind=DP)  :: time_tmp
      real(kind=DP)  :: r0_tmp(3,pub_cell%nat)
      real(kind=DP)  :: v0_tmp(3,pub_cell%nat)
      character(len=80) :: filename
   
#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_restore_trajectory'
#endif

      if (pub_on_root) then       
      
         ! Find available unit specifier
         io_unit = utils_unit()
         write(filename,'(a,a)') trim(pub_rootname),'.md.restart'
         
         open(unit=io_unit,iostat=io_stat,file=filename,&
              access='SEQUENTIAL',form='UNFORMATTED',position='REWIND', &
              action='READ')

         if (io_stat == 0) then

            ! Time, ionic positions, velocities and forces
            read(io_unit) time_tmp
            read(io_unit) ((r0_tmp(idir,iat),idir=1,3),iat=1,pub_cell%nat)
            read(io_unit) ((v0_tmp(idir,iat),idir=1,3),iat=1,pub_cell%nat)
            
            close(unit=io_unit,iostat=io_stat)
            call utils_close_unit_check('md_restore_trajectory',filename,io_stat) 
            
            if (time_tmp .gt. 0.0_dp) then
               md_time = time_tmp
               r_t = r0_tmp 
               v_t = v0_tmp 
            endif

         else
            write(stdout,'(/,3a,/,a)') ' WARNING : md_restore_trajectory failed to open "', &
              trim(filename),'"',' Initialise positions and velocities from input file...'
            md_restart = .false.
         endif

      endif

      call comms_bcast(pub_root_node_id,md_restart)
      call comms_bcast(pub_root_node_id,md_time)
      call comms_bcast(pub_root_node_id,r_t(:,:))
      call comms_bcast(pub_root_node_id,v_t(:,:))

      ! Update time, ionic positions and velocities
      if (md_restart) then
          call md_update_coordinates(elements)
      endif

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_restore_trajectory'
#endif

      return

    end subroutine md_restore_trajectory

!==============================================================================!
!==============================================================================!

    subroutine md_backup_trajectory()

      use constants,       only: stdout
      use comms,           only: pub_on_root, comms_abort
      use ion,             only: element
      use md_ions,         only: r_t, v_t, a_t, ion_mass, &
             md_ionic_polarisation, md_elec_polarisation, &
             md_total_polarisation
      use rundat,          only: pub_rootname
      use simulation_cell, only: pub_cell
      use utils,           only: utils_unit, utils_open_unit_check, &
             utils_close_unit_check
      
      implicit none

      ! Local variables
      integer           :: io_unit,io_stat
      integer           :: iat,idir
      character(len=80) :: filename

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Entering md_backup_trajectory'
#endif

      if (pub_on_root) then       
      
         ! Find available unit specifier
         io_unit = utils_unit()
         write(filename,'(a,a)') trim(pub_rootname),'.md.restart'
         
         open(unit=io_unit,iostat=io_stat,file=filename,status='REPLACE',&
              form='UNFORMATTED',position='REWIND',action='WRITE')
         call utils_open_unit_check('md_backup',filename,io_stat) 
       
         ! output is in atomic units except md_time (in fs) and temperature (in kelvin)
         write(io_unit) md_time

         ! ionic positions and velocities
         write(io_unit) ((r_t(idir,iat),idir=1,3),iat=1,pub_cell%nat)
         write(io_unit) ((v_t(idir,iat),idir=1,3),iat=1,pub_cell%nat)
         
         close(io_unit,iostat=io_stat)
         call utils_open_unit_check('md_backup',filename,io_stat) 

      endif

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(/a)') &
         'DEBUG: Leaving md_backup_trajectory'
#endif

      return

    end subroutine md_backup_trajectory

!==============================================================================!
!==============================================================================!

    subroutine md_write_trajectory(elements)


      use constants,       only: stdout
      use comms,           only: pub_on_root, comms_abort
      use ion,             only: element
      use md_ions,         only: r_t, v_t, a_t, ion_mass, &
             md_ionic_polarisation, md_elec_polarisation, &
             md_total_polarisation
      use rundat,          only: pub_rootname
      use simulation_cell, only: pub_cell
      use utils,           only: utils_unit, utils_open_unit_check, &
             utils_close_unit_check
      
      implicit none

      type(element),    intent(in) :: elements(pub_cell%nat)

      ! Local variables
      integer           :: io_unit,io_stat
      integer           :: ispec,iat
      character(len=12) :: action,form,stat,position,access
      character(len=80) :: filename

      action   = 'WRITE'
      form     = 'FORMATTED'
      stat     = 'UNKNOWN'
      position = 'REWIND'
      access   = 'SEQUENTIAL'
      
      if (md_time==0.0_dp) then
         stat     ='REPLACE'        !create a new output file
      else
         position ='APPEND'         !append to existing output file
      end if
    
 
      if (pub_on_root) then       
      
         ! Find available unit specifier
         io_unit = utils_unit()
         write(filename,'(a,a)') trim(pub_rootname),'.md'
         
         open(unit=io_unit,iostat=io_stat,file=filename,status=stat,&
              access=access,form=form,position=position,action=action)
         call utils_open_unit_check('md_write_trajectory',filename,io_stat) 
       
         ! output is in atomic units except md_time (in fs) and temperature (in kelvin)
         write(io_unit,1) md_time*aut2fs
         write(io_unit,2) md_potential_energy,md_energy,md_kinetic_energy
         write(io_unit,3) md_actual_temp*ha2k

         ! cell vectors
         write(io_unit,4) pub_cell%a1%x, pub_cell%a1%y, pub_cell%a1%z
         write(io_unit,4) pub_cell%a2%x, pub_cell%a2%y, pub_cell%a2%z
         write(io_unit,4) pub_cell%a3%x, pub_cell%a3%y, pub_cell%a3%z

         ispec=1

         ! ionic positions
         do iat=1,pub_cell%nat
            write(io_unit,7) elements(iat)%symbol,ispec,r_t(:,iat)
         enddo
         
         ! ionic velocities
         do iat=1,pub_cell%nat
            write(io_unit,8) elements(iat)%symbol,ispec,v_t(:,iat)
         end do
         
         ! ionic forces
         do iat=1,pub_cell%nat
            write(io_unit,9) elements(iat)%symbol,ispec,ion_mass(iat)*a_t(:,iat)
         end do
        
         ! Polarisation
         write(io_unit,10) md_elec_polarisation(:) 
         write(io_unit,10) md_ionic_polarisation(:) 
         write(io_unit,10) md_total_polarisation(:) 
 
         ! blank line to signal end of iteration
         write(io_unit,*) ' '
         
         close(io_unit,iostat=io_stat)
         call utils_open_unit_check('md_write_trajectory',filename,io_stat) 

      endif

1     format(12x,es18.8e3)
2     format(9x,3(3x,es18.8e3),'  <-- E')
3     format(12x,es18.8,T73,   '  <-- T')
4     format(9x,3(3x,es18.8e3),'  <-- h')
7     format(1x,a3,1x,i4,3(3x,es18.8e3),'  <-- R')
8     format(1x,a3,1x,i4,3(3x,es18.8e3),'  <-- V')
9     format(1x,a3,1x,i4,3(3x,es18.8e3),'  <-- F')
10    format(9x,3(3x,es18.8e3),'  <-- P')

      return

    end subroutine md_write_trajectory

!==============================================================================!
!==============================================================================!

    subroutine init_mts_scheme()

      use comms,     only: pub_on_root
      use constants, only: stdout
      use rundat,    only: md_delta_t, mts_nstep, mts_delta_t, &
             ngwf_threshold_orig, lnv_threshold_orig, &
             mts_ngwf_threshold, mts_lnv_threshold

      implicit none

      mts_delta_t = md_delta_t/real(mts_nstep)

      if (mts_nstep .gt. 1) then
         ref_ngwfs_threshold  = ngwf_threshold_orig
         ref_lnv_threshold    = lnv_threshold_orig
         mts_ngwfs_threshold = mts_ngwf_threshold
         mts_lnv_threshold   = mts_lnv_threshold
      else
         ref_ngwfs_threshold  = ngwf_threshold_orig
         ref_lnv_threshold    = lnv_threshold_orig
         mts_ngwfs_threshold = ngwf_threshold_orig
         mts_lnv_threshold   = lnv_threshold_orig
      endif 

      return

    end subroutine init_mts_scheme

!==============================================================================!
!==============================================================================!

    subroutine activate_mts_param()

      use rundat,  only: ngwf_threshold_orig, lnv_threshold_orig

      implicit none

      ngwf_threshold_orig = mts_ngwfs_threshold 
      lnv_threshold_orig  = mts_lnv_threshold

    end subroutine activate_mts_param

!==============================================================================!
!==============================================================================!

    subroutine deactivate_mts_param()

      use rundat,  only: ngwf_threshold_orig, lnv_threshold_orig

      implicit none

      ngwf_threshold_orig = ref_ngwfs_threshold 
      lnv_threshold_orig  = ref_lnv_threshold

    end subroutine deactivate_mts_param

!==============================================================================!
!==============================================================================!
!
!    subroutine md_step_properties(time,temp,task)
!
!
!      use comms,            only: pub_on_root
!      use constants,        only: dp, stdout
!      use rundat,           only: pub_rootname, md_formatted_output, &
!            md_pcf_calculate, md_vibra_calculate
!      use utils,            only: utils_unit, utils_open_unit_check, &
!            utils_close_unit_check
!      
!      implicit none
!
!      ! Arguments
!      real(kind=DP), intent(in) :: time 
!      real(kind=DP), intent(in) :: temp
!      character(len=*), intent(in) :: task 
!
!      ! Local variables
!      character(len=80) :: file_name
!      character(len=12) :: action,form,stat,position,access
!      integer           :: file_unit
!      integer           :: ierr, idir, iat
!
!      
!      action   = 'WRITE'
!      stat     = 'UNKNOWN'
!      position = 'APPEND'
!      access   = 'SEQUENTIAL'
!
!      if (task == 'init') then
!         stat = 'REPLACE'
!         position = 'REWIND'        
!      endif
!
!      if (md_formatted_output) then
!         form = 'FORMATTED'  
!      else
!         form = 'UNFORMATTED'  
!      end if
!
!      if (pub_on_root) then       
!
!         file_unit = utils_unit()
!         write(file_name,'(a,a)') trim(pub_rootname),'_data.md'
!         
!         open(unit=file_unit,iostat=ierr,file=file_name,status=stat,&
!              access=access,form=form,position=position,action=action)
!         call utils_open_unit_check('md_init_properties',file_name,ierr)
!         
!         if (task == 'init') then
!            ! MD parameters
!            write(file_unit,*) md_num_iter, mts_nstep
!            write(file_unit,*) md_delta_t, mts_delta_t 
!         
!            ! MD properties parameters
!            write(file_unit,*) md_pcf_calculate, md_vibra_calculate, md_irspec_calculate
!         
!            ! Computation cell 
!            write(file_unit,*) pub_cell%nat
!            write(file_unit,*) pub_cell%a1%x, pub_cell%a1%y, pub_cell%a1%z
!            write(file_unit,*) pub_cell%a2%x, pub_cell%a2%y, pub_cell%a2%z
!            write(file_unit,*) pub_cell%a3%x, pub_cell%a3%y, pub_cell%a3%z
!         
!         
!         elseif (task == 'append')
!
!            ! Time and temperature         
!            write(file_unit,*) time, temp
!         
!            ! Positions
!            if (md_pcf_calculate) write(file_unit,*) & 
!                   ((r_t(idir,iat),idir=1,3),iat=1,pub_cell%nat) 
!        
!            ! Velocities 
!            if (md_vibra_calculate) write(file_unit,*) & 
!                   ((v_t(idir,iat),idir=1,3),iat=1,pub_cell%nat) 
!         
!         endif
!         
!         close(file_unit,iostat=ierr)
!         call utils_close_unit_check('md_init_properties',file_name,ierr)
!
!      endif
!         
!    end subroutine md_step_properties
!
!==============================================================================!
!==============================================================================!

end module md
