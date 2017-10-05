! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                                 O N E T E P                                 !
!=============================================================================!
!                                                                             !
! The main ONETEP program.                                                    !
!                                                                             !
!-----------------------------------------------------------------------------!
!                                                                             !
!                 The ONETEP code is written and maintained by                !
!         Chris-Kriton Skylaris, Arash A. Mostofi, Nicholas D.M. Hine         !
!                           and Peter D. Haynes                               !
!                                                                             !
!-----------------------------------------------------------------------------!
! Arash A. Mostofi, Version 0.01, 01/11/2004                                  !
! Originally written by Chris-Kriton Skylaris in 2000.                        !
! Improved and parallelised by Chris-Kriton Skylaris in November 2003.        !
! Improvements by Peter D. Haynes in 2004.                                    !
!-----------------------------------------------------------------------------!

program onetep

  use comms, only: comms_abort, comms_barrier, comms_bcast, comms_exit, &
       comms_init, comms_groups_init, pub_on_root, pub_root_node_id, &
       pub_total_num_nodes
  use constants, only: stdout, DP
  use energy_and_force, only: energy_and_force_calculate
  use esdf, only: block_data, esdf_block, esdf_close, esdf_init, esdf_reduce, &
       esdf_stdout_dump, esdf_warnout
  use ewald, only: ewald_exit
  use geometry_optimiser
  use hubbard_init, only: hubbard_init_species_exit
  use ion, only: ELEMENT
#ifdef ACCELRYS
  use license
#endif
  use md, only: md_main
  use phonon, only: phonon_main
  use rundat, only: get_rundat, rundat_exit, pub_do_properties, pub_old_input, &
       pub_rootname, pub_write_forces, pub_write_params, pub_write_xyz, &
       read_denskern, read_tightbox_ngwfs, write_denskern, &
       write_tightbox_ngwfs, task, geom_continuation, pub_hubbard, &
       pub_cond_calculate
  use rundat_blocks, only: rundat_blocks_exec
  use services, only: services_flush, services_write_xyz
  use spherical_wave, only: sw_exit
  use timer, only: timer_clock
  use tssearch, only: tssearch_run
#ifdef DEBUG_ARRAYS
  use utils, only: utils_alloc_check, utils_dealloc_check, &
                    utils_init_array_checker, utils_array_checker_report
#else
  use utils, only: utils_alloc_check, utils_dealloc_check
#endif
#ifdef ITC_TRACE
  use vt
#endif
  use vdwcorrection, only: vdwcorrection_override_dealloc

  implicit none

  type(ELEMENT), dimension(:), allocatable   :: elements
  real(kind=DP), dimension(:,:), allocatable :: total_forces
  real(kind=DP)                              :: total_energy
  character(len=80)                          :: input_file,output_file
  character(len=80)                          :: dummy_string
  integer                                    :: nat,ierr
#ifdef ACCELRYS
  integer                                    :: num_licenses
  character (len=3)                          :: msversion="5.0"
#endif
  character (len=4)                          :: copyrightyear="2010"
  character (len=12)                         :: compilation_date
  logical                                    :: final_properties

  ! Initialise communications module
  call comms_init

  ! pdh: extract input file name from command line
  call internal_arguments

#ifdef __DATE__
    compilation_date=__DATE__
    copyrightyear=compilation_date(8:11)
#endif

#ifdef ACCELRYS
#ifdef MS_VERSION
    write(msversion,"(f3.1)") MS_VERSION
#endif
#endif

  ! pdh: print out something friendly so you know it's working...
  if (pub_on_root) call internal_welcome_banner

#ifdef ACCELRYS
  ! checkout licenses
  ! from 4.2 release there is no increase in number of checked out
  ! licenses with the number of processes; one license per job
  num_licenses = 1
  call license_checkout(num_licenses,ierr)
  if (ierr /= LIC_SUCCESS) then
     if (pub_on_root) write (stdout,'(a)') 'Error in ONETEP: &
          &failed to check out ONETEP licences'
     call comms_abort
  end if

  ! Set traps for signals
  call license_setup_traps()

#endif
#ifdef ITC_TRACE
  call vt_init_symbols()
#endif
#ifdef DEBUG_ARRAYS
  call utils_init_array_checker
#endif

  ! initialisation of the timer subroutine
  call timer_clock('total_time',0)

  ! pdh: read input file
  if (pub_on_root) call esdf_init(input_file,ierr)
  call comms_bcast(pub_root_node_id,ierr)
  if (ierr /= 0) then
     if (pub_on_root) write(stdout,'(3a)') 'Error in ONETEP: reading file "', &
          trim(input_file),'" failed'
     call comms_abort
  end if

  ! cks: read run-time parameters from file
  call comms_barrier
  if (pub_on_root) write(stdout,'(3a)',advance='no') &
       'Reading parameters from file "',trim(input_file),'" ...'
  call get_rundat
  call comms_barrier
  if (pub_on_root) write(stdout,'(a)') '... done'

  ! check for geometry block, and extract number of atoms
  if (pub_on_root) then
     if (pub_old_input) then  ! aam: flag for old-style input file
        if (.not. esdf_block('geometry',nat)) then
           write(stdout,'(a)') &
                'Error in ONETEP: no geometry block found in input file'
        end if
     else
        if (.not. esdf_block('positions_abs',nat)) then
           write(stdout,'(a)') &
                'Error in ONETEP: no positions_abs block found in input file'
        end if
        ! aam: 3/6/09 added ability to read optional  unit string
        !      ("bohr" or "ang") in first line of positions_abs block
        !      (default is bohr)
        read(block_data(1),*) dummy_string
        dummy_string=esdf_reduce(dummy_string)
        if ((index(dummy_string,"ang").ne.0).or.(index(dummy_string,"bohr").ne.0) ) nat=nat-1
     end if
  end if

  ! cks: broadcast number of atoms to all nodes
  call comms_bcast(pub_root_node_id,nat)
  if (nat == 0) call comms_abort
#ifdef DEBUG
  if (pub_on_root) write(stdout,'(a,i6,a)') 'Input file contains',nat,' atoms'
#endif

  allocate(elements(nat),stat=ierr)
  call utils_alloc_check('ONETEP (main program)','elements', ierr)

  call comms_barrier
  if (pub_on_root) write(stdout,'(3a)',advance='no') &
       'Reading geometry and species blocks from file "',trim(input_file), &
       '" ...'

  ! aam: Initialise elements and most of pub_cell
  call rundat_blocks_exec(elements,nat)

  call comms_barrier
  if (pub_on_root) write(stdout,'(a)') '... done'


  ! cks: write coordinates in xyz file
  ! qoh: but not if we are doing a geometry optimistation continuation
  if (pub_write_xyz .and. .not. geom_continuation) call services_write_xyz( &
       elements, pub_rootname, "Initial atomic coordinates")


  if (pub_on_root) then

     if (pub_write_params) then
        write(stdout,'(a)') '---------------------------------------&
             &----------------------------------------'
        write(stdout,'(a)') '---------------------------- &
             &RUN-TIME PARAMETERS ------------------------------'
        write(stdout,'(a)') '---------------------------------------&
             &----------------------------------------'

        ! ndmh: write parameters to stdout
        call esdf_stdout_dump

        write(stdout,'(a)') '---------------------------------------&
             &----------------------------------------'
        write(stdout,'(a/)') '---------------------------------------&
             &----------------------------------------'
     else
        write(stdout,*)
     end if

     ! cks: Close ESDF subroutines
     call esdf_warnout
     call esdf_close

  end if

  call services_flush

  ! Initialise groups of cores in communications module
  call comms_groups_init()

  ! Allocate forces
  allocate(total_forces(1:3,1:nat),stat=ierr)
  call utils_alloc_check('ONETEP (main program)','total_forces', ierr)

  call services_flush

  call comms_barrier

  ! aam: Determine what to do next based on value of "task"
  select case (task)

  case ('SINGLEPOINT')

     call energy_and_force_calculate(total_energy,total_forces,elements)

  case ('PROPERTIES')

     call energy_and_force_calculate(total_energy,total_forces,elements)

  case ('PHONON')

     call phonon_main(total_energy,total_forces,elements)

  case ('GEOMETRYOPTIMIZATION')

     write(output_file,'(a,a)') trim(pub_rootname),'.geom'
     if (pub_do_properties) then
        final_properties = .true.
        pub_do_properties = .false.
     else
        final_properties = .false.
     end if
     call geometry_optimise(total_energy,total_forces,elements,output_file)
     if (final_properties) then
        pub_do_properties = .true.
        call energy_and_force_calculate(total_energy,total_forces,elements)
     end if

  case ('MOLECULARDYNAMICS')

     write(output_file,'(a,a)') trim(pub_rootname),'.md'
     call md_main(total_energy,total_forces,elements)
#ifdef ACCELRYS
     call energy_and_force_calculate(total_energy,total_forces,elements, &
          properties_only=.true.)
#endif

  case ('TRANSITIONSTATESEARCH')

     write(output_file,'(a,a)') trim(pub_rootname),'.ts'
     call tssearch_run(elements,output_file)
#ifdef ACCELRYS
     call energy_and_force_calculate(total_energy,total_forces,elements, &
          properties_only=.true.)
#endif

  case ('FORCETEST')

     write(output_file,'(a,a)') trim(pub_rootname),'.forcetest'
     call internal_test_forces

  case ('HUBBARDSCF')

     write(output_file,'(a,a)') trim(pub_rootname),'.hubbardscf'
     call energy_and_force_calculate(total_energy,total_forces,elements)

  case ('TDDFT')

     write(output_file,'(a,a)') trim(pub_rootname),'.tddft'
     call energy_and_force_calculate(total_energy,total_forces,elements)

  case ('COND','PROPERTIES_COND')

     write(output_file,'(a,a)') trim(pub_rootname),'.cond'
     pub_cond_calculate = .true.
     call energy_and_force_calculate(total_energy,total_forces,elements)

  case default

     if (pub_on_root) write(stdout,'(a,a)') &
          'Error in ONETEP: illegal value for task, ', task
     call comms_abort

  end select

  !call ewald_exit(elements, total_forces)

  call rundat_exit
  call vdwcorrection_override_dealloc
  call sw_exit

  ! Deallocate forces and elements
  deallocate(total_forces,stat=ierr)
  call utils_dealloc_check('ONETEP (main program)','total_forces', ierr)

  if (pub_hubbard) then
     ! ddor: deallocate hubbard_species type array in DFT+U
     call hubbard_init_species_exit
  end if

  deallocate(elements,stat=ierr)
  call utils_dealloc_check('ONETEP (main program)','elements', ierr)

  ! shutdown of the timer subroutine
  call timer_clock('total_time',3)

#ifdef DEBUG_ARRAYS
  call utils_array_checker_report
#endif

  ! say goodbye
  if (pub_on_root) call internal_farewell

#ifdef ACCELRYS
  ! check the license back in
  call license_checkin(ierr)
#endif

  ! clean up comms module
  call comms_exit

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine internal_arguments

    !===================================================!
    ! This subroutine deals with command line arguments !
    !===================================================!

    implicit none

    ! Local variables
    integer :: num_args, len_root
#ifndef INTGETARG
    integer, external :: iargc
    external getarg
#endif

    if (pub_on_root) num_args = iargc()
    call comms_bcast(pub_root_node_id,num_args)

    if (num_args == 0) then
       input_file = 'onetep.dat'
    else
       if (pub_on_root) call getarg(1,input_file)
       call comms_bcast(pub_root_node_id,input_file)
    end if

    len_root = index(input_file,'.dat') - 1
    if (len_root > 0) then
       pub_rootname = input_file(1:len_root)
    else
       pub_rootname = input_file
    end if
    write(input_file,'(a80)') trim(pub_rootname)//'.dat'
    input_file = adjustl(input_file)
#ifdef ACCELRYS
    write(output_file,'(a80)') trim(pub_rootname)//'.onetep'
    output_file = adjustl(output_file)
    if (pub_on_root)then
       open (unit=stdout,file=output_file)
    end if
#endif
  end subroutine internal_arguments

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine internal_welcome_banner

    !========================================================!
    ! Banner printed at the beginning of ONETEP calculation. !
    !--------------------------------------------------------!
    ! Updated by Chris-Kriton Skylaris on 14/10/2004.        !
    !========================================================!

    implicit none

    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone

    call date_and_time(date,time,zone)

    write(stdout,*) '+---------------------------------------------------------------+'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|        ####### #     # ####### ####### ####### ######         |'
    write(stdout,*) '|        #     # ##    # #          #    #       #     #        |'
    write(stdout,*) '|        #     # # #   # #          #    #       #     #        |'
    write(stdout,*) '|        #     # #  #  # #####      #    #####   ######         |'
    write(stdout,*) '|        #     # #   # # #          #    #       #              |'
    write(stdout,*) '|        #     # #    ## #          #    #       #              |'
    write(stdout,*) '|        ####### #     # #######    #    ####### #              |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|        Linear-Scaling Ab Initio Total Energy Program          |'
    write(stdout,*) '|                                                               |'
#ifdef ACCELRYS
    write(stdout,*) '|                Materials Studio version '//msversion//'                   |'
#else
    write(stdout,*) '|          Release for academic collaborators of ODG            |'
    write(stdout,*) '|                                             Version 3.1.15.3  |'
#endif
    write(stdout,*) '|                                                               |'
    write(stdout,*) '+---------------------------------------------------------------+'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  Authors:                                                     |'
    write(stdout,*) '|  Peter D. Haynes, Nicholas D. M. Hine, Arash. A. Mostofi,     |'
    write(stdout,*) '|  Mike C. Payne and Chris-Kriton Skylaris                      |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  Contributors:                                                |'
    write(stdout,*) '|  J. Aarons, P. W. Avraam, S. J. Clark, F. Corsetti,           |'
    write(stdout,*) '|  N. Corsini, O. Dieguez, S. M. M. Dubois, J. Dziedzic,        |'
    write(stdout,*) '|  H. H. Helal, Q. O. Hill, D. D. O`Regan, C. J. Pickard,       |'
    write(stdout,*) '|  M. I. J. Probert, L. Ratcliff, M. Robinson, A. Ruiz Serrano, |'
    write(stdout,*) '|  and G. Teobaldi                                              |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|                                   Copyright (c) 2004-'//&
         & copyrightyear//'     |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  Please cite:                                                 |'
    write(stdout,*) '|  "Introducing ONETEP: Linear-scaling density functional       |'
    write(stdout,*) '|   simulations on parallel computers".                         |'
    write(stdout,*) '|   C.-K. Skylaris, P. D. Haynes, A. A. Mostofi, M. C. Payne.   |'
    write(stdout,*) '|   J. Chem. Phys. 122 084119 (2005).                           |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|          in all publications arising from your use of ONETEP. |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '+---------------------------------------------------------------+'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|   ONETEP is based on developments described in the following  |'
    write(stdout,*) '|   publications:                                               |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  "Nonorthogonal generalized Wannier function pseudopotential  |'
    write(stdout,*) '|   plane-wave method".                                         |'
    write(stdout,*) '|   C.-K. Skylaris, A. A. Mostofi, P. D. Haynes, O. Dieguez,    |'
    write(stdout,*) '|   M. C. Payne.                                                |'
    write(stdout,*) '|   Phys. Rev. B 66 035119 (2002).                              |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  "Preconditioned iterative minimization for linear-scaling    |'
    write(stdout,*) '|   electronic structure calculations".                         |'
    write(stdout,*) '|   A. A. Mostofi, P. D. Haynes, C.-K. Skylaris, M. C. Payne.   |'
    write(stdout,*) '|   J. Chem. Phys. 119(17), pp.8842-8848 (2003).                |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  "Total-energy calculations on a real space grid with         |'
    write(stdout,*) '|   localized functions and a plane-wave basis".                |'
    write(stdout,*) '|   A. A. Mostofi, C.-K. Skylaris, P. D. Haynes, M. C. Payne.   |'
    write(stdout,*) '|   Comput. Phys. Commun. 147, pp.788-802 (2002).               |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  "Accurate kinetic energy evaluation in electronic structure  |'
    write(stdout,*) '|   calculations with localized functions on real space grids"  |'
    write(stdout,*) '|   C.-K. Skylaris, A. A. Mostofi, P. D. Haynes, C. J. Pickard, |'
    write(stdout,*) '|   M. C. Payne.                                                |'
    write(stdout,*) '|   Comput. Phys. Commun. 140, pp.315-322 (2001).               |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '|  "Accurate ionic forces and geometry optimization in linear-  |'
    write(stdout,*) '|   scaling density-functional theory with local orbitals"      |'
    write(stdout,*) '|   N. D. M. Hine, M. Robinson, P. D. Haynes, C.-K. Skylaris,   |'
    write(stdout,*) '|   M. C. Payne, and A. A. Mostofi.                             |'
    write(stdout,*) '|   Phys. Rev. B 83 195102 (2011).                              |'
    write(stdout,*) '|                                                               |'
    write(stdout,*) '+---------------------------------------------------------------+'

#ifdef ACCELRYS
#ifdef __DATE__
#ifndef PLATFORM
#define PLATFORM "MS Windows"
#endif
#ifdef debug
#define DEBUG " DEBUG "
#else
#define DEBUG " "
#endif
    write(stdout,*)
    write(stdout,*) "This",DEBUG, "version was compiled for ",PLATFORM, &
             & " on ", __DATE__
    write(stdout,*)
#endif
#endif

    write(stdout,'(/a,2(a2,a1),a4,1x,a2,a1,a2,a2,a5,a1/)') 'Job started: ', &
         date(7:8),'-',date(5:6),'-',date(1:4),time(1:2),':',time(3:4),&
         ' (',zone,')'

    ! Write out how many processors are available
    if (pub_total_num_nodes == 1) then
       write(stdout,'(a/)') 'Running on 1 processor.'
    else
       write(stdout,'(a,i6,a/)') 'Running on',pub_total_num_nodes,' processors.'
    end if

    call services_flush

  end subroutine internal_welcome_banner

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine internal_farewell

    implicit none

    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone

    call date_and_time(date,time,zone)

    write(stdout,'(/a,2(a2,a1),a4,1x,a2,a1,a2,a2,a5,a1/)') 'Job completed: ', &
         date(7:8),'-',date(5:6),'-',date(1:4),time(1:2),':',time(3:4),&
         ' (',zone,')'

    call services_flush

  end subroutine internal_farewell

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine internal_test_forces

    !=========================================================================!
    ! Test the calculated forces against numerical finite differences         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! total_energy -                                                          !
    ! total_forces -                                                          !
    ! elements     -                                                          !
    ! nat          - number of atoms                                          !
    ! ierr         - error flag                                               !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! numerical_forces                                                        !
    ! analytical_forces                                                       !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   elements has been properly read and initialised                       !
    !   nat properly initialised                                              !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, v1.0, 08/03/2005                            !
    !=========================================================================!

    use rundat, only: maxit_pen, pub_devel_code
    use utils, only: utils_abort

    implicit none

    integer :: atom,i
    real(kind=DP) :: average_force(1:3),fd_delta_vec(1:3),fd_delta,Eplus,Eminus,Ecentral
    real(kind=DP), allocatable, dimension(:,:) :: numerical_forces
    real(kind=DP), allocatable, dimension(:,:) :: analytical_forces
    character(len=200) :: forcetest_devel_code ! ars
    integer :: fd_type, start_pos, stop_pos, test_pos ! ars
    logical :: fd_read_dkn, fd_read_tb_ngwfs ! ars
    logical :: convgd


    if (pub_on_root) then
       write(stdout,1)
       write(stdout,2) 'Starting ONETEP Force Test'
       write(stdout,1)
    endif


    ! ars: set flags
    fd_type = 1
    fd_delta = 1e-4_dp
    fd_read_dkn = .true.
    fd_read_tb_ngwfs = .true.
    if (pub_on_root) then
       forcetest_devel_code=pub_devel_code
       if (len_trim(forcetest_devel_code)>0) then
          start_pos=index(forcetest_devel_code,'FORCETEST:')
          stop_pos=index(forcetest_devel_code,':FORCETEST')
          if (stop_pos<=0) stop_pos=len_trim(forcetest_devel_code) !missing end so go to end of string
          if (start_pos>0) then

             ! ars: set finite differences scheme
             test_pos=index(forcetest_devel_code,'TYPE=')
             if (test_pos>start_pos.and.test_pos<stop_pos) then
                test_pos=test_pos+len('TYPE=')
                read(forcetest_devel_code(test_pos:test_pos+ &
                     & index(forcetest_devel_code(test_pos:stop_pos),':')-2),*) fd_type
             end if
             ! ars: set FD step
             test_pos=index(forcetest_devel_code,'DELTA=')
             if (test_pos>start_pos.and.test_pos<stop_pos) then
                test_pos=test_pos+len('DELTA=')
                read(forcetest_devel_code(test_pos:test_pos+ &
                     & index(forcetest_devel_code(test_pos:stop_pos),':')-2),*) fd_delta
             end if
             ! ars: set read_desnkern flag
             if (index(forcetest_devel_code(start_pos:stop_pos),'READ_DKN=T')>0) then
                fd_read_dkn=.true.
             else if (index(forcetest_devel_code(start_pos:stop_pos),'READ_DKN=F')>0) then
                fd_read_dkn=.false.
             end if
             ! ars: set read_tightbox_ngwfs flag
             if (index(forcetest_devel_code(start_pos:stop_pos),'READ_TB_NGWFS=T')>0) then
                fd_read_tb_ngwfs=.true.
             else if (index(forcetest_devel_code(start_pos:stop_pos),'READ_TB_NGWFS=F')>0) then
                fd_read_tb_ngwfs=.false.
             end if

          end if
       end if
    end if
    call comms_bcast(pub_root_node_id,fd_type)
    call comms_bcast(pub_root_node_id,fd_delta)
    call comms_bcast(pub_root_node_id,fd_read_dkn)
    call comms_bcast(pub_root_node_id,fd_read_tb_ngwfs)

    ! Allocate forces arrays
    allocate(numerical_forces(1:3,1:nat),stat=ierr)
    call utils_alloc_check('internal_test_forces (onetep.F90)',&
         'numerical_forces', ierr)
    allocate(analytical_forces(1:3,1:nat),stat=ierr)
    call utils_alloc_check('internal_test_forces (onetep.F90)',&
         'analytical_forces', ierr)


    ! Set flags
    if (fd_read_dkn) then
       write_denskern = .TRUE.
       if (pub_on_root) write(stdout,*) 'Force Test: Setting write_denskern = TRUE'
    endif
    if (fd_read_tb_ngwfs) then
       write_tightbox_ngwfs = .TRUE.
       if (pub_on_root) write(stdout,*) 'Force Test: Setting write_tightbox_ngwfs = TRUE'
    endif

    if (.not.pub_write_forces) then
       pub_write_forces = .TRUE.
       if (pub_on_root) write(stdout,*) 'Force Test: Setting write_forces = TRUE'
    endif

    ! Calculate total energy and forces for initial configuration
    call energy_and_force_calculate(Ecentral,total_forces,elements, &
         return_converged=convgd)
    if (.not.convgd) then
       call utils_abort('Error in internal_test_forces: Initial NGWF &
            &optimisation did not converge.')
    end if

    ! Copy
    do atom=1,nat
       analytical_forces(:,atom) = total_forces(:,atom)
    enddo

    ! ars: set read dkn and ngwfs after the central calculation
    write_denskern = .FALSE.
    if (pub_on_root) write(stdout,*) 'Force Test: Setting write_denskern = FALSE'
    write_tightbox_ngwfs = .FALSE.
    if (pub_on_root) write(stdout,*) 'Force Test: Setting write_tightbox_ngwfs = FALSE'
    if(fd_read_dkn) then
       maxit_pen = 0
       if (pub_on_root) write(stdout,*) 'Force Test: Setting maxit_pen = 0'
       read_denskern = .TRUE.
       if (pub_on_root) write(stdout,*) 'Force Test: Setting read_denskern = TRUE'
    endif
    if (fd_read_tb_ngwfs) then
       read_tightbox_ngwfs = .TRUE.
       if (pub_on_root) write(stdout,*) 'Force Test: Setting read_tightbox_ngwfs = TRUE'
    endif

    ! Calculate numerical forces
    do atom=1,nat                 ! loop over ions
       do i=1,3                   ! loop over Cartesian direction

          fd_delta_vec(:)=0.0_dp
          fd_delta_vec(i)=fd_delta

          ! ars: print banner
          if(pub_on_root) then
             write(stdout,*) "~~~~~~~~~~~~~~~~~~~~ FORCETEST ~~~~~~~~~~~~~~~~~~~~"
             write(stdout,*) "Moving atom: ", atom, "along coordinate: ", i
             write(stdout,*) "fd_delta_vec = ", fd_delta_vec(:)
             if(fd_type.eq.1) then
                write(stdout,*) "Using central finite differences"
             elseif(fd_type.eq.2) then
                write(stdout,*) "Using forward finite differences"
             elseif(fd_type.eq.3) then
                write(stdout,*) "Using backward finite differences"
             else ! ars: ==> abort
                call utils_abort("Unkonwn finite differences method for task FORCETEST")
             end if
             write(stdout,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
          end if


          if(fd_type.eq.1) then ! ars: ==> central differences

             ! Move atom in positive direction
             elements(atom)%centre%x = elements(atom)%centre%x + fd_delta_vec(1)
             elements(atom)%centre%y = elements(atom)%centre%y + fd_delta_vec(2)
             elements(atom)%centre%z = elements(atom)%centre%z + fd_delta_vec(3)

             ! Calculate new energy
             call energy_and_force_calculate(Eplus,total_forces,elements)

             ! Move atom in negative direction
             elements(atom)%centre%x = elements(atom)%centre%x - 2.0_dp*fd_delta_vec(1)
             elements(atom)%centre%y = elements(atom)%centre%y - 2.0_dp*fd_delta_vec(2)
             elements(atom)%centre%z = elements(atom)%centre%z - 2.0_dp*fd_delta_vec(3)

             ! Calculate new energy
             call energy_and_force_calculate(Eminus,total_forces,elements)

             ! Calculate force by finite-difference
             numerical_forces(i,atom) = (Eminus-Eplus)/(2.0_dp*fd_delta)

             ! Move atom back to original position
             elements(atom)%centre%x = elements(atom)%centre%x + fd_delta_vec(1)
             elements(atom)%centre%y = elements(atom)%centre%y + fd_delta_vec(2)
             elements(atom)%centre%z = elements(atom)%centre%z + fd_delta_vec(3)

          elseif(fd_type.eq.2) then ! ars: ==> forward differences

             ! Move atom in positive direction
             elements(atom)%centre%x = elements(atom)%centre%x + fd_delta_vec(1)
             elements(atom)%centre%y = elements(atom)%centre%y + fd_delta_vec(2)
             elements(atom)%centre%z = elements(atom)%centre%z + fd_delta_vec(3)

             ! Calculate new energy
             call energy_and_force_calculate(Eplus,total_forces,elements)
             if (pub_on_root) write(stdout,*) "Eplus = ", Eplus
             if (pub_on_root) write(stdout,*) "fd_delta = ", fd_delta

             ! Calculate force by finite-difference
             numerical_forces(i,atom) = (Ecentral-Eplus)/fd_delta

             ! Move atom back to original position
             elements(atom)%centre%x = elements(atom)%centre%x - fd_delta_vec(1)
             elements(atom)%centre%y = elements(atom)%centre%y - fd_delta_vec(2)
             elements(atom)%centre%z = elements(atom)%centre%z - fd_delta_vec(3)

          elseif(fd_type.eq.3) then ! ars: ==> backward differences

             ! Move atom in negative direction
             elements(atom)%centre%x = elements(atom)%centre%x - fd_delta_vec(1)
             elements(atom)%centre%y = elements(atom)%centre%y - fd_delta_vec(2)
             elements(atom)%centre%z = elements(atom)%centre%z - fd_delta_vec(3)

             ! Calculate new energy
             call energy_and_force_calculate(Eminus,total_forces,elements)

             ! Calculate force by finite-difference
             numerical_forces(i,atom) = (Eminus-Ecentral)/fd_delta

             ! Move atom back to original position
             elements(atom)%centre%x = elements(atom)%centre%x + fd_delta_vec(1)
             elements(atom)%centre%y = elements(atom)%centre%y + fd_delta_vec(2)
             elements(atom)%centre%z = elements(atom)%centre%z + fd_delta_vec(3)

          end if
       enddo
    enddo

    ! Calculate average force
    average_force=0.0_dp
    do atom=1,nat
       average_force(:) = average_force(:) + numerical_forces(:,atom)/nat
    enddo
    if (pub_on_root) write(stdout,*) "average_force(1) = ", average_force(1)

    ! Subtract average force
    do atom=1,nat
       numerical_forces(:,atom) = numerical_forces(:,atom) - average_force(:)
    enddo

    ! Write numerical forces to file
    call internal_output_forces(analytical_forces,numerical_forces)

    ! Deallocate force arrays
    deallocate(analytical_forces,stat=ierr)
    call utils_dealloc_check('internal_test_forces (onetep.F90)',&
         'numerical_forces', ierr)
    deallocate(numerical_forces,stat=ierr)
    call utils_dealloc_check('internal_test_forces (onetep.F90)',&
         'analytical_forces', ierr)

1   format(80('='))
2   format(26('<'),1x,a,1x,26('>'))

    return

  end subroutine internal_test_forces


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  subroutine internal_output_forces(analytical_forces,numerical_forces)
    !=========================================================================!
    ! Output calculated forces and numerical finite difference forces to file !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! elements                                                                !
    ! nat          - number of atoms                                          !
    ! output_file  - name of file to which forces will be written             !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, v1.0, 08/03/2005                            !
    !=========================================================================!

    use utils, only: utils_unit
    implicit none

    real(kind=DP), intent(in) :: analytical_forces(1:3,1:nat)
    real(kind=DP), intent(in) :: numerical_forces(1:3,1:nat)

    integer :: atom,write_unit,ios

    if (pub_on_root) then

       ! Find available unit specifier
       write_unit = utils_unit()

       ! Open file
       open(unit=write_unit,file=output_file,iostat=ios,&
            form='FORMATTED',action='WRITE')
       if (ios/=0) then
          write(stdout,'(3a)') &
               'Error in ONETEP: failed to open file "',trim(output_file),'"'
          call comms_abort
       endif

       ! Write analytical forces to file
       write(write_unit,'(a)') ' '
       write(write_unit,'(a)') '******************* Analytical Forces ********************'
       write(write_unit,'(a)') '*                                                        *'
       write(write_unit,'(a)') '* Element  Atom         Cartesian components (Eh/a)      *'
       write(write_unit,'(a)') '* ------------------------------------------------------ *'
       write(write_unit,'(a)') '*                       x            y            z      *'
       write(write_unit,'(a)') '*                                                        *'
       do atom=1,nat
          write(write_unit,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
               elements(atom)%symbol,' ',atom,'   ', &
               analytical_forces(1,atom),analytical_forces(2,atom),analytical_forces(3,atom),' *'
       enddo
       write(write_unit,'(a)') '*                                                        *'
       write(write_unit,'(a)') '**********************************************************'

       ! Write numerical forces to file
       write(write_unit,'(a)') ' '
       write(write_unit,'(a)') '******************** Numerical Forces ********************'
       write(write_unit,'(a)') '*                                                        *'
       write(write_unit,'(a)') '* Element  Atom         Cartesian components (Eh/a)      *'
       write(write_unit,'(a)') '* ------------------------------------------------------ *'
       write(write_unit,'(a)') '*                       x            y            z      *'
       write(write_unit,'(a)') '*                                                        *'
       do atom=1,nat
          write(write_unit,'(a1,a6,a1,i6,a3,3f13.8,a2)') '*', &
               elements(atom)%symbol,' ',atom,'   ', &
               numerical_forces(1,atom),numerical_forces(2,atom),numerical_forces(3,atom),' *'
       enddo
       write(write_unit,'(a)') '*                                                        *'
       write(write_unit,'(a)') '**********************************************************'

       ! Close file
       close(write_unit,iostat=ios)
       if (ios/=0) then
          write(stdout,'(3a)') 'Error in ONETEP: failed to close file "',trim(output_file),'"'
          call comms_abort
       endif

    endif


    return

  end subroutine internal_output_forces


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end program onetep

