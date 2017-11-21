! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                    G E O M E T R Y                                          !
!=============================================================================!
!                                                                             !
! $Id: geometry_optimiser_mod.F90,v 1.37 2009/09/23 16:56:54 cks22 Exp $
!                                                                             !
!-----------------------------------------------------------------------------!
! This module performs a geometry optimisation on the specified model.        !
! It may be used to modify ionic positions and/or cell parameters.            !
!                                                                             !
!-----------------------------------------------------------------------------!
! Written from  "Module Specification for New Plane Wave Code",  M. Segall,   !
!    P. Lindan, M. Probert, C. Pickard, P. Hasnip, S. Clark and M. Payne      !
!                           Copyright 1999/2000                               !
!-----------------------------------------------------------------------------!
! Written by Matt Probert, v1.0, 11/05/2002                                   !
!-----------------------------------------------------------------------------!
! Modified for ONETEP by Arash A Mostofi                                      !
! Cleaned up by Alvaro Ruiz Serrano in September 2009.                        !
!-----------------------------------------------------------------------------!
!
! $Log: geometry_optimiser_mod.F90,v $
! Revision 1.37  2009/09/23 16:56:54  cks22
!
! cks_2009_09_23:
! -> Alvaro's clean up of the geometry_optimiser
! -> Fix geometry optimisation continuations
! -> Add XLYP dispersion parameters
! -> Faster calculation of spherical harmonics
! -> Various minor fixes and cleanups
!
! Revision 1.36  2009/08/14 10:39:03  ndmh3
!
! ndmh_2009_08_14: Addition of SPAM3, function_basis and new function sum/integrals
! routines. Cleanup of whitespace, capitalisation and indentation. Bugfixes to
! TDDFT and DFT+U. Minor format changes, rearrangement of atom counting information
! to show projector count. New approach to nonlocal potential matrices in
! norm_conserv_pseudo, to go with SPAM3. Improvements to Peano Space Filling Curve
! distribution in parallel_strategy_mod. Rearrangement of ngwf_gradient_mod to
! group routines relating to construction of coefficients for function sums.
!
! Revision 1.35  2009/06/23 09:28:13  cks22
!
! cks_2009_06_23: Bug fixes for spin polarisation, an additional kernel fix approach,
! BLYP parameters and periodic boundary conditions for the VdW module, clean-up of
! various subroutines.
!
! Revision 1.34  2009/05/11 14:32:40  ndmh3
!
! ndmh_2009_05_11: Removed old-style usage of output_detail, replaced with
! pub_output_detail throughout. Fixed memory allocation problem when using
! NNHO's with dense_matrices: T. Tidied up some formatting in properties_mod.
! Fixed a badly formatted error message in energy_and_force_mod.
! qoh: fixed bug in spherical_wave_mod introduced in 2.3.13
!
! Revision 1.33  2009/04/22 15:11:54  cks22
!
! cks_2009_04_22: several minor bug fixes and the inclusion of "classical" point
! charges to the external potential.
!
! Revision 1.32  2009/02/05 03:20:59  cks22
!
! cks_2009_02_05: Hartree-Fock exchange (slow algorithm, using the FFT box equal
! to the simulation cell). Also a bug fix by Nick Hine in
! restart_ngwfs_tightbox_input to follow the F95 standard and fix wrong reading
! of NGWFs with pg compilers.
!
! Revision 1.19  2009/02/02 17:53:02  cks
!
! cks_2009_02_02: Version 2.2.16
!
! Revision 1.31  2009/01/26 16:49:33  ndmh3
! ndmh_2009_01_26: fixed step sizes in geometry_optimiser_mod by fixing bug
! in inverse Hessian initialisation involving units of mass (1 m_e /= 1 amu!)
!
! Revision 1.5  2008/12/01 14:08:43  vmilman
! Merge with the updated ODG version 2.2.14
!
! Revision 1.29  2008/10/28 11:09:56  cks22
!
! cks_2008_10_28: Trivial changes to
! 1) Eliminate "cell" as an argument and replace it with "pub_cell".
! 2) Fix a recent bug in the initialisation of tight boxes which made the code
! crash when the number of ppds was not the same along a1 and a2.
! 3) Make the code run without crashing when the error bounds chacking is used
! with pgf90.
!
! Revision 1.28  2008/09/22 18:43:45  ndmh3
!
! ndmh_22_09_2008: addendums to previous checkin, cache-optimisation of
! norm_conserv_pseudo_mod structure factor code, various cleanups and minor
! fixes, adaptable output during penalty and lnv for very large systems (no
! change for small systems).
!
! Revision 1.27  2008/08/27 15:08:48  cks22
!
! cks_2008_08_27: Restored working of empirical dispersion in the code.
!Also added the input flag geom_reuse_dk_ngwfs (true by default) which when
!set to false it prevents reading the density kernel and ngwfs from the
!previous geometry step as this can sometimes cause (mysterious) ngwf
!convergence problems.
!
! Revision 1.26  2008/08/13 15:34:14  ndmh3
!
! ndmh_2008_08_13: file omitted from previous check-in
!
! Revision 1.25  2008/07/31 17:48:24  cks22
!
! cks_2008_07_31: (1) Trivial changes in spherical_mod and (2) incorporation of
! empirical dispersion forces and energy in order to be able to use dispersion
! contributions in geometry optimisation.
!
! Revision 1.12  2008/07/31 17:41:16  cks
!
! cks_2008_07_31: VdW forces and geopt.
!
! Revision 1.24  2008/07/02 16:39:22  pdh1001
! Minor changes to satisfy gfortran.
!
! Revision 1.23  2008/03/30 14:56:57  cks22
!
! cks_2008_03_30: Memory deallocation bug fixes and write_xyz now writes xyz
! history during optimisation.
!
! Revision 1.22  2008/02/15 10:38:17  pdh1001
! Several minor bug fixes.
! Added NGWF analysis functionality.
! Added change from Victor to population analysis output.
!
! Revision 1.21  2007/08/19 08:54:54  cks22
!
! cks_2007_08_19: Extension to use simulation cells smaller than the 6 NGWF
! radii rule, along
! one or more lattice vector directions. The minimum simulation cell is now
!defined by the twice
! the maximum NGWF radius.
!
! Revision 1.3  2007/08/17 17:17:14  cks
!
! cks_2007_08_17: Merged flex and current versions.
!
! Revision 1.20  2007/02/17 00:08:40  pdh1001
! Added capability to write sparse matrices in SPAM3 format to a binary file
! in a versioned, extensible format.
!
! Revision 1.19  2007/01/31 14:10:45  pdh1001
! Victor's changes to TS search and delocalised internals, plus Windows fix.
!
! Revision 1.18  2006/12/11 11:08:38  pdh1001
! Minor changes to avoid continuation bug: number of spins now included in
! .continuation file and SPAM%matrix initialised to zero after allocation.
!
! Revision 1.17  2006/11/16 13:23:55  pdh1001
! Bug fix for spin polarisation in backup, restore and continuation file routine
!
! Revision 1.16  2006/09/08 16:32:46  pdh1001
! FFTw version 3 interface added (compile with -DFFTW3 instead of -DFFTW)
! Energy components printed in BRIEF mode
! Other superficial changes
!
! Revision 1.15  2006/05/15 08:20:35  cks22
!
! cks_2006_05_15: added properties module, brief output level, Victor's
! latest TS changes, ngwf initialisation up to element 53.
!
! Revision 1.14  2006/04/21 15:42:09  cks22
!
! cks_2004_04_21: Added Victor's transition state search files.
! Also added ngwf_halo, a new experimentlal feature which can take into account
! small but non-zero matrix elements.
!
! Revision 1.13  2006/02/10 16:45:07  pdh1001
! Merged with source from Victor Milman (mainly to incorporate delocalized
! internals) and some minor bug fixes.
!
! Revision 1.12  2005/09/30 13:22:06  aam24
! Small change to supress inverse Hessian printout.
!
! Revision 1.11  2005/09/27 10:23:45  aam24
! New logical keywords "geom_print_inv_hessian" and "verbose_ewald_forces" for
! controlling whether to write to standard output the inverse hessian
! during a geometry optimisation run and debugging information on the
! Ewald forces, respectively.
!
! Revision 1.10  2005/09/13 15:17:54  aam24
! Parallel bug fix in subroutine geom_opt_continuation_read for geometry
! optimiser continuation facility.
!
! Revision 1.9  2005/08/31 14:39:29  aam24
! Small change to make task=GeometryOptimization work with old_input=TRUE.
!
! Revision 1.8  2005/08/04 14:10:32  aam24
! Output absolute cartesian (rather than fractional) coordinates in geometry
! optimisation.
!
! Revision 1.7  2005/08/04 11:02:58  aam24
! First implementation of ionic constraints.
!
! Revision 1.6  2005/07/26 10:24:18  aam24
! Geometry optimiser now able to cope with fractional positions outside
! the interval [0,1]. Minor changes to norm_conserve_pseudo, electronic,
! energy_and_force and onetep.F90.
!
! Revision 1.5  2005/06/08 10:16:10  aam24
! Minor structural changes to geometry_optimiser and onetep
! modules. Average residual force on atoms is now explicitly
! removed. Addition of new "task" and subroutines for testing
! forces. Addition of a function in energy_and_force module to count the
! number of psinc functions enclosed in NGWF spheres (you need to
! uncomment the relevant call).
!
! Revision 1.4  2005/02/28 17:10:26  aam24
! Fixing changes in onetep.F90 that were accidentally overwritten in the
! last revision. Very minor change in geometry_optimiser module.
!
! Revision 1.3  2005/02/23 18:32:52  aam24
! Additional geometry optimisation parameters in rundat and esdf_key,
! including logical keyword "geom_continuation" which allows one to
! resume a previous geometry optimisation. Corresponding changes to
! geometry_optimiser module and onetep.F90. Addition of two "physical"
! types in esdf_mod.F90 to deal with units of "ha/bohr**3" for
! rundat parameter "geom_modulus_est" (estimate of the bulk modulus) and
! "ha/bohr" for "geom_force_tol" (convergence tolerance for the forces).
!
! Revision 1.2  2005/01/27 17:49:22  aam24
! Small changes to geometry_optimiser module to get around an apparent
! Intel compiler bug.
!
! Revision 1.1  2004/12/02 11:59:09  aam24
! Addition of geometry optimisation module with corresponding minor
! changes to some other modules. Geometry optimisation is invoked with
! the 'task' rundat parameter. Addition of two periodic table arrays
! to constants_mod.F90.
!
! Revision 1.60  2003/11/20 17:20:22  mijp2
! catch case of non-default fixed_npw with non-varying cell [bug 03324vymq01]
!
! Revision 1.59  2003/11/03  17:29:24  mijp2
! made frac_symm_ops allocatable to fix bug 03304vymq05
!
! Revision 1.58  2003/10/31  14:58:10  mijp2
! added constraints bug-fix in deloc-internals by Victor
!
! Revision 1.57  2003/09/30  09:30:52  mijp2
! minor forchk update
!
! Revision 1.56  2003/08/06 13:53:27  mijp2
! fixed latest restart bug
!
! Revision 1.55  2003/08/05  14:55:41  mijp2
! fixed bug 03148vymq02
!
! Revision 1.54  2003/07/28  17:00:33  mijp2
! added fixed_npw
! forchk neatening up
!
! Revision 1.53  2003/07/16  16:55:37  mijp2
! fixed variable cell, variable #pw restart bug
! fixed nb=0 bug
!
! Revision 1.52  2003/05/24  18:13:03  mijp2
! fixed bug in backup/restore with const Ecut cell changes
!
! Revision 1.51  2003/05/22  09:57:14  mijp2
! better purification of noise in frac_symm_ops
! changed convergence path in variable cell with const_pw=.false.
!
! Revision 1.50  2003/05/12  21:01:29  mijp2
! corrected intent() errors
!
! Revision 1.49  2003/05/07  15:02:10  mijp2
! added missing initialisation for last_write_time
!
! Revision 1.48  2003/04/30  22:11:40  mijp2
! Added saves to module variables.
! Added delocalised internals minimiser (Victor)
! Changed to fixed Ecutoff not fixed #plane waves
!
! Revision 1.47  2003/04/11  16:52:41  mijp2
! fixed run_time restart bug
!
! Revision 1.46  2003/03/20  11:35:04  mijp1
! FORCHK caught a type error in unsymm_to_symm
!
! Revision 1.44  2003/02/07  09:37:02  mijp2
! fixed parallel bug in timed backups
!
! Revision 1.43  2003/01/29  21:29:00  mijp2
! added timed backups using backup_interval
!
! Revision 1.42  2002/11/15  11:06:21  mijp2
! renamed %nbands -> %nbands_max and %nbands_occ -> %nbands
!
! Revision 1.41  2002/08/23 13:22:52  mijp2
! added off-diagonal external pressure (ie shear) capability
!
! Revision 1.40  2002/07/29  15:12:42  mijp2
! tweaked previous RCS log message as it causes pgf90 problems!
!
! Revision 1.39  2002/07/29  11:04:48  mijp2
! added frequency_label,
! deleted off-diagonal external_pressures (bug 02128s1wq01),
! clarified cell_constraints warning,
! added some extra checks on lambda to ensure minimum cell vector lengths,
! fixed bug in hs_vec calculation in geom_update_inv_Hessian_params,
! minor neatening of output,
! explicit dE/ion in convergence table,
! only do bisection if trying to guarantee downhill only steps.
! Phew!
!
! Revision 1.38  2002/07/16  12:42:02  mijp2
! extended geom_symm_delta to operate on ionic displacements as well as cell
! strains
!
! Revision 1.37  2002/07/12  16:24:16  mijp2
! fixed bug 02159vymq01 - BFGS sometimes violates cell constraints
! by the addition of new private routine geom_symm_delta which
! ensures that the strain obeys symmetry before the cell is updated.
!
! Revision 1.36  2002/07/09  14:20:48  mijp2
! fixed bug 02168vymq01 (final output forces showed zeroed constraints)
! changed output to only print cell info if cell changed and v.v. ions
!
! Revision 1.35  2002/05/16  16:08:27  mijp2
! fixed subtle restart bug [02136vymq01]
!
! Revision 1.34  2002/05/14  16:08:09  mijp2
! corrected dE bug in geom_converged [bug 02134vymq01]
!
! Revision 1.33  2002/05/13  12:30:33  mijp2
! tweaked uphill / quad heuristics to fix bug
!                    02119g1eq03
!
! Revision 1.32  2002/05/12  00:09:42  mijp2
! numerous minor bug fixes and tweaks:
! improved output messages to reduce user confusion [bugs 02120s1wq01 and
! 02119g1eq03]
! revised geom_check_lambda [bug 02114s1wq01]
! updated force output with constraints [bug 02114vymq02]
! deleted fake_forces as no longer needed
! improved backtracking if need to go uphill
! tweaked quad minimisation heuristic
! no longer rationalise ionic coords in .castep output
! sanitised final messages
! removed windows on |F|max, |dR|max and Smax convergence
!
! Revision 1.31  2002/04/24  17:09:02  mijp2
! bug fix in random noise
!
! Revision 1.30  2002/04/19  12:26:40  mijp2
! updated FBSC routine to be in sync with castep.f90 v1.30
!
! Revision 1.29  2002/04/18  16:30:56  mijp2
! deleted cell_update_symmetry (now in cell module) and fixed up call-logic
! beefed-up constraints & symmetry checks for cell and ions
! minor FBSC bug fixed
! catch instability in inv_Hessian for long runs and correct it
!
! Revision 1.28  2002/04/16  16:51:01  mijp2
! random noise bug fix
! changed order of basis_initialise and model_cell_changed to sort out beta-phi
! bug
! correctly updating symmetry if cell vectors change
! tweaked termination conditions if try to bisect and in the force/stress noise
! floor
! catch system that is fully constrained by symmetry
! correctly update the |dR|max field if backtrack a step
!
! Revision 1.27  2002/04/03  13:05:02  mijp2
! undid fudge factor
! tweaked error messages
!
! Revision 1.26  2002/03/04  17:26:37  mijp2
! changed BFGS backup file to scratch type so not visible to user
! made FBSC more efficient
!
! Revision 1.25  2002/02/18  16:20:37  mijp2
! added geom_update_params to catch user changes to modulus/frequency
! estimates and update inverse Hessian appropriately
! tweaked o/p some more
!
! Revision 1.24  2002/02/15  01:25:11  mijp2
! added nullify methods
! added DM stuff
! changed output formats
! neatened output and verbosity
! neatened up convergence output
!
! Revision 1.23  2002/02/07  17:41:41  mijp2
! fixed format error
!
! Revision 1.22  2002/02/07  17:21:49  mijp2
! allowed for on-the-fly unit changes in output
!
! Revision 1.21  2002/02/07  12.38.53  mijp2
! added check for no relaxation required
!
! Revision 1.20  2002/02/04  15:06:51  mijp2
! changed LAPACK usage for portability between v2 and v3
!
! Revision 1.19  2002/01/30  19:00:54  mijp2
! changed FBSC usage as problems with doing it mid-minimization
!
! Revision 1.18  2002/01/30  17:20:03  mijp2
! changed force & stress output
! added fermi energy to fbsc
! changed wave functions to S versions
!
! Revision 1.17  2002/01/17  17:40:03  mijp2
! fixed transposed cell vectors, and added efermi
!
! Revision 1.16  2001/12/19 10:54:43  mijp2
! fixed <-- E formatting in trajectory file
!
! Revision 1.15  2001/12/11  15:50:32  mijp2
! rewrote perfectly good F90 in geom_get_forces to workaround pgf90 compiler
! bug with -fast optimisation
!
! Revision 1.14  2001/11/28  09:42:57  mijp2
! fixed weird SGI bug with anint function
!
! Revision 1.13  2001/11/24  14:33:39  mijp2
! added file_maxpath
! changed .geom file format - now includes a header, cell vectors even for
! fixed cell calculation, non-rationalized positions, and everything in AU
! tweaks to output
! restart bug fixed
!
! Revision 1.12  2001/11/08  23:56:56  mijp2
! re-introduced finite basis correction calculation as in main
! deleted geom_ucase
! only o/p modulus and frequency if different to initial values
!
! Revision 1.11  2001/10/31  22:55:14  mijp2
! commented out one offending line missed in last checkin
!
! Revision 1.10  2001/10/31  22:36:42  mijp2
! removed finite basis set correction calculation - now in main
!
! Revision 1.9  2001/10/12  22:27:19  mijp2
! updated format statement
!
! Revision 1.8  2001/09/29  00:00:56  mijp2
! parallel IO bugfix, trajectory output bugfix
! improved inv_Hess reinitialise
! improved check_lambda
! improved IO
!
! Revision 1.7  2001/07/23  17:02:33  mijp2
! improved output
! fixed bug in |dR| calculation
!
! Revision 1.6  2001/07/19  10:50:51  mijp1
! added bisection search for line minimization failure situation
!  replaced quadratic lambda calculation algorithm
! tidied up logic for backtracking/further minimization
! implemented opt_strategy for fast backtracking
!
! Revision 1.5  2001/07/15  00:41:56  mijp2
! minor tweaks to logic and backup/backtracking
! some improvements to output formatting
!
! Revision 1.4  2001/07/03  02:35:42  mijp2
! finite basis corrections
!
! Revision 1.2  2001/06/27  22:06:05  mijp2
! added trajectory file output
!
! Revision 1.1  2001/06/15  14:31:18  mijp2
! Initial revision
!
!                                                                             !
!=============================================================================!

module geometry_optimiser

  use constants, only: DP
  use simulation_cell, only: castep_cell

  implicit none                                 !Impose strong typing

  private                                       !Everything is private ...

     public :: geometry_optimise

     ! Cell data in CASTEP format
     type(castep_cell)  :: current_cell

     ! Fake constants
     integer, parameter :: file_maxpath = 256

     ! Constants and parameters
     real(kind=DP), parameter :: geom_stress_tol      = 0.1_dp
     real(kind=DP), dimension(6), parameter :: external_pressure=0.0_dp
     logical, save :: constrain_ions
     logical, save :: fix_all_ions
     integer, save :: num_ionic_constraints
     integer, save :: ndim  !the number of dimensions in the search space
     !so for the augmented-Hessian BFGS method used here,
     !we have ndim=9+3*num_ions as we seek to optimize the
     !cell (9 vars) and the ions (3*num_ions d.of.f) simultaneously.

     integer, save :: ndof  !the number of degrees of freedom in the problem
     ! - may be less than ndim due to symmetry and/or constraints

     integer, parameter :: shuffle_forwards  =  1
     integer, parameter :: shuffle_none      =  0
     integer, parameter :: shuffle_backwards = -1

     integer, parameter, dimension(1:3,1:3) :: &
          unit_matrix=reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))

     integer, save            :: reset_iter

     real(kind=DP), save      :: old_frequency_est
     real(kind=DP), save      :: new_frequency_est, new_modulus_est

     !Set up the pretty unit labels
     character(len=*), parameter :: length_label    = 'Bohr'
     character(len=*), parameter :: energy_label    = 'Ha'
     character(len=*), parameter :: force_label     = 'Ha/Bohr'
     character(len=*), parameter :: pressure_label  = 'Ha/Bohr**3'
     character(len=*), parameter :: frequency_label = '1/aut'

     character(len=file_maxpath) :: cont_file

   contains

     subroutine geometry_optimise(total_energy,total_forces,elements,output_file)

       !=========================================================================!
       ! Find ground state ionic positions and/or cell vectors for given model   !
       ! Also sets the private flags constrain_ions, symm_ions, etc.             !
       !-------------------------------------------------------------------------!
       ! Arguments:                                                              !
       !   mdl, intent=inout, the model to be optimised                          !
       !   output_file, intent=in, the filename to which results will be written.!
       !-------------------------------------------------------------------------!
       ! Parent module variables used:                                           !
       !   model_file is set to checkpoint parameter                             !
       !   constrain_ions, symm_ions, constrain_cell, and symm_cell are set here !
       !   for all minimizers.                                                   !
       !   ndof and ndim are set here for all minimizers                         !
       !   pretty printing labels are set for all minimizers                     !
       !   on_root flag is set as appropriate for output on parallel machines    !
       !-------------------------------------------------------------------------!
       ! Modules used:                                                           !
       !   md for damped-md optimiser                                            !
       !   model for type definitions                                            !
       !   parameters to define the job etc                                      !
       !   cell to setup constraints flags                                       !
       !   io for error messages                                                 !
       !   comms to set on_root flag                                             !
       !-------------------------------------------------------------------------!
       ! Key Internal Variables:                                                 !
       !-------------------------------------------------------------------------!
       ! Necessary conditions:                                                   !
       !   mdl has already been properly read and initialised                    !
       !-------------------------------------------------------------------------!
       ! Written by Matt Probert, v1.0, 11/05/2002                               !
       !-------------------------------------------------------------------------!
       ! Modified for ONETEP by Arash A Mostofi                                  !
       !=========================================================================!

       use classical_pot, only: classical_elements, classical_pot_ii_energy, classical_pot_ii_forces
       use comms, only: pub_on_root, comms_abort, comms_barrier
       use constants, only: dp, stdout, NORMAL
       use ion, only: element
       use geometry, only: operator(.DOT.)
       use energy_and_force, only: energy_and_force_calculate
       use forces, only: forces_apply_constraints
       use rundat, only: write_denskern, write_tightbox_ngwfs, read_denskern, &
            read_tightbox_ngwfs, pub_geom_reuse_dk_ngwfs, maxit_pen, &
            geom_method, pub_rootname, pub_old_input, geom_continuation, &
            pub_write_positions
       use simulation_cell, only: castep_model, pub_cell, &
            copy_to_castep_cell, castep_model_nullify, castep_cell_copy, &
            castep_model_alloc, castep_cell_2_elements, castep_model_dealloc, &
            castep_cell_dealloc

       implicit none

       character(len=*), intent(in) :: output_file
       type(ELEMENT), intent(inout) :: elements(pub_cell%nat)
       real(kind=DP), intent(inout) :: total_energy
       real(kind=DP), intent(inout) :: total_forces(1:3,pub_cell%nat)

       ! <<< local variables >>>
       type(castep_model) :: mdl
       integer :: ion_i, ierr, ion_in_species

       ! aam: total number of ionic constraints
       num_ionic_constraints = 0
       do ion_i=1,pub_cell%nat
          ! aam: old input format does not support ionic constraints
          if (pub_old_input) elements(ion_i)%ion_constraint_type='NONE'
          select case (elements(ion_i)%ion_constraint_type)
          case ('NONE')  ; continue
          case ('LINE')  ; num_ionic_constraints = num_ionic_constraints + 1
          case ('PLANE') ; num_ionic_constraints = num_ionic_constraints + 2
          case ('FIXED') ; num_ionic_constraints = num_ionic_constraints + 3
          case default
             if (pub_on_root) then
                write(stdout,'(a,i6,a)') 'Error in geometry_optimise: &
             &illegal value for elements(',ion_i,')%ion_constraint_type'
                call comms_abort
             endif
          end select
       enddo

       ! ars 0 <= num_ionic_constraints <= 3*pub_cell%nat
       ! aam: Analyse constraints
       num_ionic_constraints=min(3*pub_cell%nat,num_ionic_constraints)  !catch any sillies
       num_ionic_constraints=max(0,num_ionic_constraints)               !catch any sillies
       if (num_ionic_constraints>0) constrain_ions=.true.
       ! ars : correct initialisation of fix_all_ions
       if (num_ionic_constraints==3*pub_cell%nat) then
          fix_all_ions=.true.
       else
          fix_all_ions=.false.
       end if

       ndof = 3*pub_cell%nat - num_ionic_constraints ! aam: cell is fixed ; no symmetry

          ndim = 9 + 3*pub_cell%nat             !we always allow for both cell & ions here regardless
                                                !of constraints to allow for changes upon restart, etc.

       !Make sure there is something here to do
       if (fix_all_ions) then
          if (pub_on_root) write (stdout,'(a)') &
          'WARNING - there is nothing to optimise - skipping relaxation'
          return
       end if

       ! aam: Write a little banner
       if (pub_on_root) then
          write(stdout,'(a)') ''
          if (geom_continuation) then
             write(stdout,1)
        write(stdout,2) 'Resuming previous ONETEP Geometry Optimisation'
             write(stdout,1)
          else
             write(stdout,1)
             write(stdout,3) 'Starting ONETEP Geometry Optimisation'
             write(stdout,1)
          endif
          write(stdout,'(/a,a,a)') &
               ' Geometry Optimisation: output file is "',trim(output_file),'"'
       endif

       if (pub_on_root) write(stdout,*) ''

        ! fill CASTEP-style cell data structures from pub_cell and "elements"
        ! kaw: Directly analogous to cks changes for energy calculation
       if (pub_cell%nat_classical > 0) then
        call copy_to_castep_cell(current_cell,elements,classical_elements)
       else
        call copy_to_castep_cell(current_cell,elements)
       end if

       ! Consistency check
       if (current_cell%num_ions<=0) then
          if (pub_on_root) write(stdout,'(a)') &
      'Error in ONETEP: current_cell uninitialised in geometry_optimise'
          call comms_abort
       endif

       ! Nullify pointers in mdl
       call castep_model_nullify(mdl)

       ! Allocate model data
       call castep_model_alloc(mdl,ndim,constrain_ions,.false.)

       ! Over-ride the density kernel and tightbox_ngwf write flags if necessary
       if (.not.write_denskern .and. pub_geom_reuse_dk_ngwfs) then
          write_denskern=.TRUE.
          if (pub_on_root) write(stdout,*) &
      'Geometry Optimisation: write_denskern parameter over-ridden to T'
       endif
       if (.not.write_tightbox_ngwfs .and. pub_geom_reuse_dk_ngwfs) then
          write_tightbox_ngwfs=.TRUE.
          if (pub_on_root) write(stdout,*) &
              'Geometry Optimisation: write_tightbox_ngwfs parameter over-ridden &
              &to TRUE'
       endif
       ! End over-rides


       !**************************************************************!
       !                     Copy model data                          !
       !**************************************************************!

       ! Name of continuation file
       write(cont_file,'(a,a)') trim(pub_rootname),'.continuation'

       if (geom_continuation) then

          ! Read continuation data into model and current_cell.
          ! NB This subroutine overwrites current_cell%ionic_positions
          call geom_opt_continuation_read(mdl,elements,cont_file)

          call castep_cell_2_elements(current_cell,elements)

       else

  !mdl%cell is either not being allocated correctly for classical_atoms or is being used incorrectly
          call castep_cell_copy(current_cell,mdl%cell)

          if (pub_write_positions) call geom_output(mdl)

          ! If this is not a continuation of a previous geometry optimisation:
          ! Find ground state energy and forces for initial configuration
          call energy_and_force_calculate(total_energy,total_forces,elements)

          ! Total energy
          mdl%total_energy = total_energy

          ! Forces
          do ion_i=1,pub_cell%nat                                   ! Each ion has its own species
             do ion_in_species=1,current_cell%max_ions_in_species   ! ie max_ions_species = 1
                mdl%forces(1,ion_in_species,ion_i) = total_forces(1,ion_i)
                mdl%forces(2,ion_in_species,ion_i) = total_forces(2,ion_i)
                mdl%forces(3,ion_in_species,ion_i) = total_forces(3,ion_i)
                if (constrain_ions) then ! aam: keep copy of unconstrained forces
                   mdl%orig_forces(1,ion_in_species,ion_i) = total_forces(1,ion_i)
                   mdl%orig_forces(2,ion_in_species,ion_i) = total_forces(2,ion_i)
                   mdl%orig_forces(3,ion_in_species,ion_i) = total_forces(3,ion_i)
                endif
             enddo
          enddo

          ! BFGS data
          mdl%bfgs_iteration     = 0
          if (pub_on_root) then
             mdl%bfgs_inv_Hessian   = 0.0_dp
          end if

       endif

       mdl%found_ground_state = .true.
       mdl%stress             = 0.0_dp
       mdl%strain             = 0.0_dp

       ! aam: apply ionic constraints
       if (constrain_ions) call forces_apply_constraints(mdl%forces,elements)

       call castep_cell_copy(current_cell,mdl%cell)
       call castep_cell_copy(current_cell,mdl%orig_cell)
       call castep_cell_copy(current_cell,mdl%ref_cell)

       !**************************************************************!
       !                    Model data copy completed                 !
       !**************************************************************!


       ! Over-ride the density kernel and tightbox_ngwf read flags if necessary
       ! This is required for back-tracking in the BFGS algorithm
       ! cks, 2008/08/26: Added the ability to prevent the usage of the
       ! cks: density kernel and NGWFs from the previous geometry step
       ! cks: as this can occasionally cause convergence problems in the energy.
       if (.not.read_denskern .and.  pub_geom_reuse_dk_ngwfs) then
          read_denskern=.TRUE.
          if (pub_on_root) write(stdout,*) &
       'Geometry Optimisation: read_denskern parameter over-ridden to T'
       endif
       if (.not.read_tightbox_ngwfs .and.  pub_geom_reuse_dk_ngwfs) then
          read_tightbox_ngwfs=.TRUE.
          if (pub_on_root) write(stdout,*) &
               'Geometry Optimisation: read_tightbox_ngwfs parameter over-ridden &
               &to TRUE'
       endif
       ! cks: do not do penalty if NGWFs and kernel are re-used
       if (pub_geom_reuse_dk_ngwfs) maxit_pen = 0
       ! End over-rides

       ! aam: output initial atomic positions in fractional coordinates
       if (pub_write_positions) call geom_output(mdl)
       call comms_barrier

       !Now select geometry minimizer (might want to add others in the future)
       select case (geom_method)
       case ('CARTESIAN')

          call geom_BFGS(mdl,elements,output_file)

       case ('DELOCALIZED')

          call geom_DELOC(mdl,elements,output_file,ierr)

          if (ierr/=0) then
             geom_method = 'CARTESIAN'
             call geom_BFGS(mdl,elements,output_file)
          endif

       case default
          write (stdout,'(a)') 'Error in geom_optimise - unrecognised geometry&
               &optimisation request='//geom_method//'.'
          call comms_abort
       end select

       ! Deallocate memory

       call castep_model_dealloc(mdl,.false.)

       call castep_cell_dealloc(current_cell)

   1   format(80('='))
   2   format(16('<'),1x,a,1x,16('>'))
   3   format(20('<'),1x,a,1x,21('>'))

       return
     end subroutine geometry_optimise


     !---------------------------------------------------------------------------!
     !                P R I V A T E   R O U T I N E S                            !
     !---------------------------------------------------------------------------!

     subroutine geom_BFGS(mdl,elements,output_file)
       !=========================================================================!
       ! Find ground state ionic positions and/or cell vectors for given model.  !
       ! Based upon "Relaxation of Crystals with the Quasi-Newton Method" by     !
       !   B.G. Pfrommer et al, J.Comp.Phys. v131, p233-240 (1997).              !
       !                                                                         !
       ! We are in the constant pressure, adiabatic ensemble, and so work with   !
       !   the enthalpy: H=Etot+p^ext.V and seek to minimize this.               !
       ! Minimizing wrt ionic coords -> equilibrium condition: force=0           !
       !   and min. wrt cell strains -> equilibrium condition: stress+pressure=0 !
       ! That is, we use the following sign convention for stress & pressure:    !
       !   A positive external pressure tends to compress the system, and        !
       !   a positive internal stress   tends to expand   the system.            !
       !                                                                         !
       ! NB We use the ref_cell as the strain reference if doing variable cell.  !
       !                                                                         !
       ! NB As of version 1.48 the minimiser can be either constant Ecut or      !
       !    constant #PW, and as of version 1.54 this can be set at run-time via !
       !    the fixed_npw parameter.                                             !
       !    Constant Ecut is most physical and works best if long way from eqm   !
       !   - but can get trapped into finding wrong minimum (eg Fe).             !
       !    Constant PW can give smoother convergence near the minimum but is    !
       !    physical. Might be tempted to try to hybridise these schemes, but I  !
       !    found you can get nasty discontinuities if switch mid-minimisation   !
       !     => better to set it once before minimisation starts and not change. !
       !                                                                         !
       ! NB If variable cell with fixed_npw=F then we must abandon our strictly  !
       !    downhill in energy strategy. As we have a changing number of plane   !
       !    waves then there will be discontinuities in the energy as we step    !
       !    around the cell vectors and so an E-based search will go wild. Hence !
       !    the search direction must be based upon reducing F & S rather than E.!
       !-------------------------------------------------------------------------!
       ! Arguments:                                                              !
       !   mdl, intent=inout, the model to be optimised                          !
       !   output_file, intent=in, the filename to which results will be written.!
       !-------------------------------------------------------------------------!
       ! Parent module variables used:                                           !
       !   on_root used for diagnostics                                          !
       !-------------------------------------------------------------------------!
       ! Modules used:                                                           !
       !   model for type definitions                                            !
       !   parameters                                                            !
       !   cell for cell vectors and constraints                                 !
       !   io for output and unit conversion routines                            !
       !-------------------------------------------------------------------------!
       ! Key Internal Variables:                                                 !
       !   x_vec       =  ndim   vector of strains & fractional coords           !
       !   f_vec       =  ndim   vector of stresses & forces (without the metric)!
       !   delta_vec   =  ndim   vector of x_vec update: new_x=old_x+lambda*delta!
       !   lambda      =  scalar = line minimization parameter along the proposed!
       !                  BFGS update direction delta                            !
       !   trial_OK    =  logical flag to say whether or not trial step has gone !
       !                  up or down in enthalpy. Used if need to backtrack.     !
       !                  ditto for line_OK and quad_OK.                         !
       !                                                                         !
       !-------------------------------------------------------------------------!
       ! Necessary conditions:                                                   !
       !   mdl has already been properly read and initialised                    !
       !-------------------------------------------------------------------------!
       ! Written by Matt Probert, v1.0, 11/05/2002                               !
       !-------------------------------------------------------------------------!
       ! Modified for ONETEP by Arash A Mostofi, 2004                            !
       !=========================================================================!

       use constants, only: dp, stdout, NORMAL
       use comms, only: pub_on_root, pub_root_node_id, comms_abort, comms_barrier, &
            comms_bcast
       use ion, only: element
       use rundat, only: geom_energy_tol,geom_frequency_est, pub_output_detail, &
            geom_max_iter,geom_modulus_est,geom_backup_iter,pub_rootname, &
            pub_write_positions,geom_reset_dk_ngwfs_iter,pub_geom_reuse_dk_ngwfs, &
            read_denskern,read_tightbox_ngwfs
       use simulation_cell, only: castep_cell, castep_model, pub_cell, &
            castep_cell_nullify, castep_cell_alloc, castep_cell_copy, &
            castep_cell_dealloc
       use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
                        utils_open_unit_check, utils_close_unit_check, &
                        utils_read_check, utils_write_check, utils_flush

       implicit none

       type(castep_model), intent(inout) :: mdl
       type(element),      intent(inout) :: elements(1:pub_cell%nat)
       character(len=*),   intent(in)    :: output_file
       type(castep_cell) :: init_cell
       real(kind=DP), allocatable, dimension(:)     :: x_vec, old_x_vec !strains and fractional coords
       real(kind=DP), allocatable, dimension(:)     :: f_vec, old_f_vec !derivatives of enthalpy wrt x_vec
       real(kind=DP), allocatable, dimension(:)     :: delta_vec        !update x_vec

       type(castep_cell)                            :: backup_cell
       real(kind=DP)                                :: backup_total_energy,backup_enthalpy
       real(kind=DP), allocatable, dimension(:,:,:) :: backup_forces
       real(kind=DP), dimension(1:6)                :: backup_stress
       real(kind=DP), dimension(1:3,1:3)            :: backup_strain
       real(kind=DP), allocatable, dimension(:)     :: backup_x_vec, backup_f_vec
       integer                                      :: backup_dk_unit,backup_ngwf_unit !backup_unit
       character (len=file_maxpath)                 :: backup_dk_file,backup_ngwf_file !backup_file

       integer                                      :: io_status

       real(kind=DP), allocatable, dimension(:)     :: trial_f_vec, line_f_vec, quad_f_vec

       real(kind=DP) :: enthalpy, old_enthalpy, trial_enthalpy, line_enthalpy, quad_enthalpy

       logical       :: abort_BFGS  !have we got stuck?
       real(kind=DP) :: fudge       !amount of uphill ignored if monotonic_E=false

       logical       :: trial_OK, line_OK, quad_OK
       logical       :: done_trial, done_line, done_quad
       logical       :: revert_old

       real(kind=DP) :: dE
       real(kind=DP) :: Fmax   !max Force
       real(kind=DP) :: dRmax  !max displacement
       real(kind=DP) :: Smax   !max Stress
       logical       :: converged, converged_dE, converged_Fmax, converged_dRmax, converged_Smax

       real(kind=DP) :: lambda_trial,lambda_line,lambda_quad,lambda_backtrack
       logical       :: lambda_diff

       real(kind=DP) :: old_f_dot_delta, trial_f_dot_delta, line_f_dot_delta, quad_f_dot_delta
       real(kind=DP) :: df_dot_delta, df_dot_delta_dlambda

       integer       :: i,ierr
       logical       :: re_initialise

#ifdef DEBUG
       if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_BFGS'
#endif

       !-------------------------------------------------------------------------!
       !Initialize data structures required                                      !
       !-------------------------------------------------------------------------!

       abort_BFGS=.false.
       fudge=geom_energy_tol*real(mdl%cell%num_ions,kind=dp) !absolute fudge in enthalpy for uphill steps

       old_frequency_est=geom_frequency_est             !we store the starting values as 'old' to test o/p at end
       new_frequency_est=geom_frequency_est             !we store the working value as 'new' s.t. can cope with
       new_modulus_est  =geom_modulus_est               !automatic and user updates to the estimates

       !NB we properly nullify any pointers before first usage

       call castep_cell_nullify(init_cell)

       call castep_cell_nullify(backup_cell)


       !We need a local cell as that at which the inv_Hessian was last initialised ...

       call castep_cell_alloc(init_cell)

       call castep_cell_copy(mdl%cell,init_cell)


       !Then we set up the ndim arrays
       allocate (x_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','x_vec',ierr)
       x_vec=0.0_dp

       allocate (old_x_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','old_x_vec',ierr)
       old_x_vec=0.0_dp

       allocate (f_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','f_vec',ierr)
       f_vec=0.0_dp

       allocate (old_f_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','old_f_vec',ierr)
       old_f_vec=0.0_dp

       allocate (delta_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','delta_vec',ierr)
       delta_vec=0.0_dp

       !then we setup the backup arrays as needed

       call castep_cell_alloc(backup_cell)


       ! AAM: Set file names for NGWF and density kernel backup files
       write(backup_ngwf_file,'(2a)')&
            trim(pub_rootname),'.tightbox_ngwfs.bfgs_back'

       write(backup_dk_file,'(2a)') trim(pub_rootname),'.dkn.bfgs_back'


       if (pub_on_root) write(stdout,'(/3a)') &
            ' BFGS: continuation file name is "',trim(cont_file),'"'


       !kaw: We only need forces for the QM ions so the final dimension of this was too large
       !allocate(backup_forces(1:3,1:mdl%cell%max_ions_in_species,1:mdl%cell%num_species),stat=ierr)
       allocate(backup_forces(1:3,1:mdl%cell%max_ions_in_species,1:pub_cell%nat),stat=ierr)
       call utils_alloc_check('geom_BFGS','backup_forces',ierr)
       backup_forces=0.0_dp

       allocate (backup_x_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','backup_x_vec',ierr)
       backup_x_vec=0.0_dp

       allocate (backup_f_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','backup_f_vec',ierr)
       backup_f_vec=0.0_dp

       allocate (trial_f_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','trial_f_vec',ierr)
       trial_f_vec=0.0_dp

       allocate (line_f_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','line_f_vec',ierr)
       line_f_vec=0.0_dp

       allocate (quad_f_vec(1:ndim),stat=ierr)
       call utils_alloc_check('geom_BFGS','quad_f_vec',ierr)
       quad_f_vec=0.0_dp

       !-------------------------------------------------------------------------!
       !Sort out start from restart                                              !
       !-------------------------------------------------------------------------!

       !If this is not a restart, then need to initialise the inv_Hessian
       if (mdl%bfgs_iteration==0) then

          call geom_inv_Hessian_initialize(init_cell,mdl)
          mdl%strain=0.0_dp

          !Put the strains wrt original cell and fractional atomic coords wrt current_cell into x_vec
          call geom_mdl_to_xvec(mdl,x_vec)

          !Get stresses and forces from mdl and put into f_vec
          call geom_mdl_to_fvec(mdl,f_vec,enthalpy)

          !start up the output file
          call geom_write_trajectory(mdl,enthalpy,output_file)

       else
          if (pub_on_root) write (stdout,*) 'BFGS: Restarting previous &
     &                                       geometry optimization.'

          !Catch the possibility that model_continuation/reuse has decided to nuke the inv_Hessian
          re_initialise = .false.
          if (pub_on_root) then
             if (all(abs(mdl%bfgs_inv_Hessian)<tiny(0.0_dp))) re_initialise = .true.
          end if
          call comms_bcast(pub_root_node_id,re_initialise)
          if (re_initialise) then
             call castep_cell_copy(mdl%cell,init_cell)
             call geom_inv_Hessian_initialize(init_cell,mdl)
          end if

          !Put the strains wrt original cell and fractional atomic coords wrt current_cell into x_vec
          call geom_mdl_to_xvec(mdl,x_vec)

          !Get stresses and forces from mdl and put into f_vec
          call geom_mdl_to_fvec(mdl,f_vec,enthalpy)

       end if

       !now prepare to relax the structure ...
       !calculate delta_vec (the update step)
       !NB Here and elsewhere we use matmul rather than BLAS for clarity,
       !   as efficiency is not an issue here (matrix sizes are either 3*3 or ndim*ndim)
       !   and electronic_minimisation always dominates.
       if (pub_on_root) then
          delta_vec=matmul(mdl%bfgs_inv_Hessian,f_vec)
       end if
       call comms_bcast(pub_root_node_id, delta_vec)

       !store old x_vec, f_vec and enthalpy
       old_x_vec=x_vec
       old_f_vec=f_vec
       old_enthalpy=enthalpy

       !NB We reset bfgs_iteration here as the geom_max_iter refers to number of iterations
       !   THIS pass and not in total.
       mdl%bfgs_iteration=0
       reset_iter=mdl%bfgs_iteration

       !initialise convergence tests (impossible to converge here with convergence windows)
       call geom_converged(mdl,old_x_vec,enthalpy,shuffle_forwards,dE,Fmax,dRmax,Smax, &
            & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax)
       !29/09/05: Actually it is possible to converge here now - we no longer use windows on
       !          Fmax, dRmax and Smax, and if the initial configuration passes all those but
       !          not dE(as that still has a window) then I think it is OK to call this
       !          structure converged and save some CPU time!

       !tell user where we are starting from
       call geom_output_converged(mdl%bfgs_iteration,enthalpy,dE,Fmax,dRmax, &
            converged_dE,converged_Fmax,converged_dRmax)

       !-------------------------------------------------------------------------!
       !BEGIN main loop over successive configurations as we minimize cell & ions!
       !-------------------------------------------------------------------------!
       main_loop: do

       !check to see if we have done enough iterations yet?
       !or perhaps already converged?
       !NB Test is done this way around s.t. we don't need to undo anything at end
       if (mdl%bfgs_iteration>=geom_max_iter .or. converged) exit main_loop

       mdl%bfgs_iteration=mdl%bfgs_iteration+1

       !banner to output
       if (pub_on_root) write (stdout,4) 'Starting BFGS iteration ',mdl%bfgs_iteration,' ...'

       ! aam: write geometry optimiser continuation data?
       if (mod(mdl%bfgs_iteration,geom_backup_iter)==0) &

            call geom_opt_continuation_write(mdl,elements,cont_file)

       ! ndmh: reset dk and NGWFs this iteration?
       if ((geom_reset_dk_ngwfs_iter>0).and.pub_geom_reuse_dk_ngwfs) then
          if (mod(mdl%bfgs_iteration,geom_reset_dk_ngwfs_iter)==0) then
             if (read_denskern) then
                read_denskern=.FALSE.
                if (pub_on_root) write(stdout,*) &
                     'Geometry Optimisation: Resetting: read_denskern &
                     &parameter set to FALSE'
             endif
             if (read_tightbox_ngwfs) then
                read_tightbox_ngwfs=.FALSE.
                if (pub_on_root) write(stdout,*) &
                     'Geometry Optimisation: Resetting: read_tightbox_ngwfs &
                     &parameter set to FALSE'
             endif
           end if
       end if


       !store key data in case of backtrack (including old x_vec, f_vec etc)
       call geom_BFGS_backup

       !set accept/reject flags
       done_trial=.false.
       trial_OK  =.false.
       done_line =.false.
       line_OK   =.false.
       done_quad =.false.
       quad_OK   =.false.

       !----------------------------------------------------------------------!
       !BEGIN trial move                                                      !
       !----------------------------------------------------------------------!

       !store F^old.delta for lambda line minimisation
       old_f_dot_delta=dot_product(f_vec,delta_vec)
       call geom_BFGS_write_summary

       !propose a trial move with lambda=1
       lambda_trial=1.0_dp

       !check trial lambda does not lead to too large a displacement wrt old structure ...
       call geom_check_lambda(lambda_trial,delta_vec)


       !generate trial structure
       x_vec=old_x_vec+lambda_trial*delta_vec             !x=x+dx - includes strains & ionic coords

       !diagnostic
       if (pub_output_detail>NORMAL.and.pub_on_root) then
          write (stdout,7)    'BFGS: trial lambda =',lambda_trial
          write (stdout,10)   'BFGS: trial:','i','xvec','old_xvec','delta_vec'
!'
          do i=1,ndim
             write (stdout,1) 'BFGS: trial:',i,x_vec(i),old_x_vec(i),delta_vec(i)
          end do
       end if

       !update model and current_cell for new cell & positions
       call geom_xvec_to_mdl(x_vec,mdl)

       !tell user what's about to happen
       if (pub_on_root)&
            write (stdout,5) 'BFGS: starting iteration',mdl%bfgs_iteration, &
            & ' with trial guess (lambda=',lambda_trial
       if (pub_write_positions) call geom_output(mdl)  ! info on current_cell to stdout


       !evaluate trial structure
       call geom_get_forces(mdl,elements)


       !get stresses and forces from mdl and put into f_vec
       call geom_mdl_to_fvec(mdl,f_vec,enthalpy)
       trial_f_vec=f_vec

       !store F^trial.delta for lambda line minimisation
       trial_f_dot_delta=dot_product(f_vec,delta_vec)

       !setup backtrack logic
       trial_enthalpy=enthalpy
       done_trial=.true.
       trial_OK=.false.

       if (abs(trial_f_dot_delta)<=abs(old_f_dot_delta)) trial_OK=.true.


       !inform user about what is going on ...
       call geom_BFGS_write_summary

       !test convergence
       call geom_converged(mdl,old_x_vec,enthalpy,shuffle_forwards,dE,Fmax,dRmax,Smax, &
            & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax)

       !tell user what just happened
       if (converged) then
          call geom_output_converged(mdl%bfgs_iteration,enthalpy,dE,Fmax,dRmax,&
               converged_dE,converged_Fmax,converged_dRmax)
          exit main_loop
       end if

       ! ndmh: stop resetting dk and NGWFs if we just did so above
       if ((geom_reset_dk_ngwfs_iter>0).and.pub_geom_reuse_dk_ngwfs) then
          if (mod(mdl%bfgs_iteration,geom_reset_dk_ngwfs_iter)==0) then
             if (.not.read_denskern) then
                read_denskern=.TRUE.
                if (pub_on_root) write(stdout,*) &
                     'Geometry Optimisation: read_denskern parameter &
                     &set to TRUE'
             endif
             if (.not.read_tightbox_ngwfs) then
                read_tightbox_ngwfs=.TRUE.
                if (pub_on_root) write(stdout,*) &
                     'Geometry Optimisation: read_tightbox_ngwfs parameter &
                     &set to TRUE'
             endif
          end if
       end if

       !line minimisation for best lambda is safe -> lambda_line
       !NB If the energy surface is quadratic, then lambda=1 is optimal, but
       !   if we start some way from the minimum this is unlikely, so we do
       !   a line minimization instead.
       !BUT for efficiency we keep the trial structure if lambda_line~1
       !   (unless it has gone uphill) to prevent needless calls to geom_get_forces!!!
       if (abs(trial_f_dot_delta-old_f_dot_delta)<epsilon(old_f_dot_delta)) then
          if (pub_on_root) then
             write (stdout,*) 'BFGS: Looks like this system is as converged as possible.'
             write (stdout,*) '      Maybe your geometry convergence tolerances are too tight?'
          end if

          !test convergence
          call geom_converged(mdl,old_x_vec,enthalpy,shuffle_none,dE,Fmax,dRmax,Smax, &
               & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax)
          !spoof the main flag
          converged=.true.

          !report the result
          call geom_output_converged(mdl%bfgs_iteration,enthalpy,dE,Fmax,dRmax,&
               converged_dE,converged_Fmax,converged_dRmax)
          exit main_loop
       end if

       !so we know that abs(df_dot_delta)>epsilon ...
       df_dot_delta=trial_f_dot_delta-old_f_dot_delta
       lambda_line=-lambda_trial*old_f_dot_delta/df_dot_delta

       !check line minimisation does not lead to too large a displacement wrt old structure ...
       if ((1.0_dp-lambda_trial)>epsilon(lambda_trial).and.lambda_line>lambda_trial) then
          !lambda_trial has been truncated and line>trial => line will be truncated too => nothing will change
          !so this is a time NOT to check the proposed new lambda and let a big displacment go through
          if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
             'BFGS: trial lambda was truncated for safety so NOT truncating new&
             & lambda'
       else
          call geom_check_lambda(lambda_line,delta_vec)
       end if


       !best line lambda not close enough to trial lambda (Pfrommer tolerance)
       !or trial structure has gone uphill from end of previous iteration and lambda_line<>lambda_trial
       ! => evaluate new trial structure
       lambda_diff=abs(lambda_line-lambda_trial)>epsilon(lambda_line)
       if ((abs(lambda_line-lambda_trial)>0.6_dp) .or. &
            & (((trial_enthalpy-fudge)>old_enthalpy).and.lambda_diff)) then

          !-------------------------------------------------------------------!
          !BEGIN line-minimisation move                                       !
          !-------------------------------------------------------------------!

          !store key data in case of backtrack (including trial x_vec, f_vec etc)
          if (trial_OK) call geom_BFGS_backup

          !setup new state of system
          x_vec=old_x_vec+lambda_line*delta_vec

          !diagnostic
          if (pub_output_detail>NORMAL.and.pub_on_root) then
             write (stdout,7)    'BFGS: line minimisation lambda=',lambda_line
             write (stdout,10)   'BFGS: line :','i','xvec','old_xvec','delta_vec'
             do i=1,ndim
                write (stdout,1) 'BFGS: line :',i,x_vec(i),old_x_vec(i),delta_vec(i)
             end do
          end if

          !update model and current_cell for new cell & positions
          call geom_xvec_to_mdl(x_vec,mdl)

          !tell user what's about to happen
          if (pub_on_root) write (stdout,5) 'BFGS: improving iteration',mdl%bfgs_iteration, &
               & ' with line minimization (lambda=',lambda_line
          if (pub_write_positions) call geom_output(mdl) !info on current_cell to stdout

          !evaluate line minimisation structure
          call geom_get_forces(mdl,elements)

          !get stresses and forces from mdl and put into f_vec
          call geom_mdl_to_fvec(mdl,f_vec,enthalpy)
          line_f_vec=f_vec

          !store F^line.delta for quad minimisation
          line_f_dot_delta=dot_product(f_vec,delta_vec)

          !setup backtrack logic
          line_enthalpy=enthalpy
          done_line=.true.
          line_OK=.false.

             if (abs(line_f_dot_delta)<=abs(trial_f_dot_delta)) then
                trial_OK=.false.
                line_OK=.true.
             end if

          !inform user about what is going on ...
          call geom_BFGS_write_summary

          !test convergence
          call geom_converged(mdl,old_x_vec,enthalpy,shuffle_none,dE,Fmax,dRmax,Smax, &
               & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax)

          !tell user what just happened
          if (converged) then
             call geom_output_converged(mdl%bfgs_iteration,enthalpy,dE,Fmax,dRmax,&
                  converged_dE,converged_Fmax,converged_dRmax)
             exit main_loop
          end if

          !Is quadratic interpolation required?
          !We actually do this using simple Newton-Raphson and assuming quadratic form for (F.delta) vs lambda
          !This rarely makes a big difference, except when the line minimisation step has failed to go downhill
          ! whereupon it is worth doing everything possible to avoid backtracking to old structure
          ! and resetting the inv_Hessian.
          !Thats also why I've added the monotonic_E flag - sometimes going uphill gets around a barrier
          !if all else has failed!

          !we know that lambda_trial-lambda_line>epsilon and lambda_line<>0 from if blocks above
          df_dot_delta_dlambda=(trial_f_dot_delta-line_f_dot_delta)/(lambda_trial-lambda_line)
          if (lambda_line<lambda_trial) then
             !can do two-sided estimate (averaging gradients) to get a best value
             df_dot_delta_dlambda=0.5_dp*(df_dot_delta_dlambda+(line_f_dot_delta-old_f_dot_delta)/lambda_line)
          else
             !cannot do any more as one sided estimate is all that is possible I'm afraid!
          end if
          if (abs(df_dot_delta_dlambda)>tiny(df_dot_delta_dlambda)) then
             lambda_quad=lambda_line-line_f_dot_delta/df_dot_delta_dlambda
          else
             !this implies that the function is linear and so lambda_line must be the best answer!
             lambda_quad=lambda_line
          end if

          !check quadratic lambda does not lead to too large a displacement wrt old structure ...
          if ((1.0_dp-lambda_trial)>epsilon(lambda_trial).and.lambda_quad>lambda_trial) then
             !lambda_trial has been truncated and quad>trial => quad will be truncated too => nothing will change
             !so this is a time NOT to check the proposed new lambda
             if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
                'BFGS: trial lambda was truncated for safety so NOT truncating &
                &new lambda'
          else
             !compare actual to expected value of lamba_line (not saved)
             if ((lambda_line+lambda_trial*old_f_dot_delta/(trial_f_dot_delta-old_f_dot_delta))<epsilon(lambda_line) &
                  & .and.lambda_quad>lambda_line) then
                !lambda_line has been truncated and quad>line => quad will be truncated too => nothing will change
                !so this is a time NOT to check the proposed new lambda
                if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
                'BFGS: line lambda was truncated for safety so NOT truncating &
                &new lambda'
             else
                call geom_check_lambda(lambda_quad,delta_vec)
             end if
          end if

          !best quadratic lambda not close enough to line minimisation lambda (my tolerance: difference > 30%)
          !or line minimisation structure has gone uphill from trial structure and lamba_quad<>lambda_line
          lambda_diff=abs(lambda_quad-lambda_line)>epsilon(lambda_quad)
          if (  ((trial_enthalpy-fudge)>old_enthalpy).and. &
               &   ((line_enthalpy-fudge)>old_enthalpy).and.lambda_diff)  then

             !----------------------------------------------------------------!
             !BEGIN quadratic-minimisation move                               !
             !----------------------------------------------------------------!

             !store key data in case of backtrack (including line minimisation x_vec, f_vec etc)
             if (line_OK) call geom_BFGS_backup

             !set up new state of system
             x_vec=old_x_vec+lambda_quad*delta_vec

             !diagnostic
             if (pub_output_detail>NORMAL.and.pub_on_root) then
                write (stdout,7)    'BFGS: quad lambda =',lambda_quad
                write (stdout,10)   'BFGS: quad :','i','xvec','old_xvec','delta_vec'
                do i=1,ndim
                   write (stdout,1) 'BFGS: quad :',i,x_vec(i),old_x_vec(i),delta_vec(i)
                end do
             end if

             !update model and current_cell for new cell & positions
             call geom_xvec_to_mdl(x_vec,mdl)

             !tell user what's about to happen
             if (pub_on_root)write(stdout,5) 'BFGS: improving iteration',mdl%bfgs_iteration, &
                  & ' with quad minimization (lambda=',lambda_quad
             if (pub_write_positions) call geom_output(mdl)              !info on current_cell to stdout

             !evaluate quad structure
             call geom_get_forces(mdl,elements)

             !get stresses and forces from mdl and put into f_vec
             call geom_mdl_to_fvec(mdl,f_vec,enthalpy)
             quad_f_vec=f_vec

             !store F^quad.delta for neat output
             quad_f_dot_delta=dot_product(f_vec,delta_vec)

             !setup backtrack logic
             quad_enthalpy=enthalpy
             done_quad=.true.
             quad_OK=.false.

                if (abs(quad_f_dot_delta)<=abs(line_f_dot_delta)) then
                   trial_OK=.false.
                   line_OK=.false.
                   quad_OK=.true.
                end if

             !inform user about what is going on ...
             call geom_BFGS_write_summary

             !test convergence
             call geom_converged(mdl,old_x_vec,enthalpy,shuffle_none,dE,Fmax,dRmax,Smax, &
                  & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax)

             !tell user what just happened
             if (converged) then
                call geom_output_converged(mdl%bfgs_iteration,enthalpy,dE,Fmax,dRmax,&
                     converged_dE,converged_Fmax,converged_dRmax)
                exit main_loop
             end if

             !----------------------------------------------------------------!
             !END quadratic-minimisation move                                 !
             !----------------------------------------------------------------!

          else !lambda_quad close enough to lambda_line so no need to do any more
             if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
                'BFGS: line minimization close enough - no need to evaluate &
                &quadratic step'

             !----------------------------------------------------------------!
             !END line-minimisation move                                      !
             !----------------------------------------------------------------!

          end if

       else !lambda_line close enough to trial (lambda=1) so no need to change anything
          if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
     'BFGS: trial structure close enough - no need to evaluate line step'

          !-------------------------------------------------------------------!
          !END trial move                                                     !
          !-------------------------------------------------------------------!

       end if


       !----------------------------------------------------------------------!
       !BEGIN backtrack code                                                  !
       !----------------------------------------------------------------------!

       !check if we need to backtrack the proposed structure
       lambda_backtrack=-999.9_dp    !set to impossible value
       revert_old =.false.
       if (done_quad) then           !we may or may not have done a quad step, so check first
          if (quad_OK) then
             if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
                'BFGS: accepting quad step - no need to backtrack'
          else                       !quad   step has gone uphill so reject and revert to line step
             if (line_OK) then       !line step went downhill so accept line step
                if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
                   'BFGS: backtracking from quad to line'
                lambda_backtrack=lambda_line
             else                    !line step has gone uphill so reject and revert to trial step
                if (trial_OK) then   !trial  step went downhill so accept trial step
                   if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
                      'BFGS: backtracking from quad to trial'
                   lambda_backtrack=lambda_trial
                else                 !trial  step has gone uphill so reject and revert to old structure
                   if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
                      'BFGS: backtracking from quad to old'
                   lambda_backtrack=0.0_dp
                   revert_old=.true.
                end if
             end if
          end if
       else if (done_line) then    !we may or may not have done a line step, so check next
          if (line_OK) then
             if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
                'BFGS: accepting line step - no need to backtrack'
          else                       !line step has gone uphill so reject and revert to trial step
             if (trial_OK) then      !trial  step went downhill so accept trial step
                if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
                   'BFGS: backtracking from line to trial'
                lambda_backtrack=lambda_trial
             else                    !trial  step has gone uphill so reject and revert to old structure
                if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
                   'BFGS: backtracking from line to old'
                lambda_backtrack=0.0_dp
                revert_old=.true.
             end if
          end if
       else                          !we always do a trial step, so check last.
          if (trial_OK) then
             if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
                'BFGS: accepting trial step - no need to backtrack'
          else                       !trial  step has gone uphill so reject and revert to old structure
             if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,*) &
                'BFGS: backtracking from trial to old'
             lambda_backtrack=0.0_dp
             revert_old=.true.
          end if
       end if

       !warn user we are about to backtrack
       if (revert_old) then

          if (pub_on_root) write (stdout,*) 'BFGS: reverting to earlier &
     &                                       configuration'
          abort_BFGS=.false.
          call geom_BFGS_revert_old !calculate new value for abort_BFGS
          if (abort_BFGS) then
             exit main_loop     !no point continuing
          else
             cycle main_loop    !successful backtrack, do fresh iteration
          end if

       else                     !not revert_old

          if (lambda_backtrack>-999.0_dp) then !ie we have assigned a definite value in check-logic above

             x_vec=old_x_vec+lambda_backtrack*delta_vec
             call geom_xvec_to_mdl(x_vec,mdl)

             !tell user what's about to happen
             if (pub_on_root) write (stdout,*) 'BFGS: reverting to &
                                                earlier configuration'
             if (pub_write_positions) call geom_output(mdl)              !info on current_cell to stdout

             !NB Presume above logic is consistent with geom_BFGS_backup logic s.t. we get the right restore ...
             call geom_BFGS_restore

          end if                !lambda_backtrack>0

       end if                   !revert_old

       !----------------------------------------------------------------------!
       !END backtrack code                                                    !
       !----------------------------------------------------------------------!

       !----------------------------------------------------------------------!
       !BEGIN finish off this iteration                                       !
       !----------------------------------------------------------------------!

       !update inv_Hessian using BFGS
       call geom_update_inv_Hessian(old_x_vec,x_vec,old_f_vec,f_vec,init_cell,mdl)

       !apply inv_Hessian to get new delta
       if (pub_on_root) then
          delta_vec=matmul(mdl%bfgs_inv_Hessian,f_vec)
       end if
       call comms_bcast(pub_root_node_id,delta_vec)

       !test convergence
       call geom_converged(mdl,old_x_vec,enthalpy,shuffle_none,dE,Fmax,dRmax,Smax, &
            & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax)

       !tell user what just happened
       call geom_output_converged(mdl%bfgs_iteration,enthalpy,dE,Fmax,dRmax,&
            & converged_dE,converged_Fmax,converged_dRmax)

       !are we done yet?
       if (converged) exit main_loop

       !update output file
       call geom_write_trajectory(mdl,enthalpy,output_file)

       !store old x_vec, f_vec and enthalpy
       old_x_vec=x_vec
       old_f_vec=f_vec
       old_enthalpy=enthalpy


       !----------------------------------------------------------------------!
       !END finish off this iteration                                         !
       !----------------------------------------------------------------------!

      end do main_loop
      !-------------------------------------------------------------------------!
      !END main loop over successive configurations as we minimize cell & ions  !
      !-------------------------------------------------------------------------!

      !-------------------------------------------------------------------------!
      !BEGIN final output and tidy up                                           !
      !-------------------------------------------------------------------------!

      !write final entry to the output file
      if (converged) call geom_write_trajectory(mdl,enthalpy,output_file)


      !Actually, we need to UNSET the found_forces flag if we have any ionic
      !constraints s.t. check_forces_stresses in main will recalculate the
      !unconstrained forces and so produce non-zero outputs with firstd_output_forces
      !and then reset the found_forces flag before the call to model_write in main.

      if (pub_on_root) then
       if (converged) then
          write(stdout,*) 'BFGS: Geometry optimization completed successfully.'
       else
          write(stdout,*) 'BFGS: Geometry optimization failed to converge &
            &after',mdl%bfgs_iteration,'steps.'
       end if
       write (stdout,11) 'BFGS: Final Configuration:'
       call geom_output(mdl)
       write (stdout,2) 'BFGS: Final Enthalpy     =',enthalpy,trim(energy_label)
      end if

      if (mdl%bfgs_iteration-reset_iter>1) call geom_update_inv_Hessian_params(init_cell,mdl)
      if (pub_on_root) then
         if (.not.fix_all_ions) then
            if (abs((new_frequency_est-old_frequency_est)/old_frequency_est)>0.01_dp) then !warn if 1% change
               write (stdout,2) 'BFGS: Final <frequency>  =',new_frequency_est,trim(frequency_label)
            end if
         end if
      end if

    ! aam: store as cont_file whether converged or not
    call geom_opt_continuation_write(mdl,elements,cont_file)

    !clean up at the end
    call castep_cell_dealloc(backup_cell)

    ! ===============================================================================
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<< Remove backup files >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ===============================================================================
    if (pub_on_root) then
       ! NGWF backup file first....
       ! Find next available unit specifier
       backup_ngwf_unit = utils_unit()
       ! tell user what's going on
       write(stdout,'(/3a)',advance='no') &
            ' Removing NGWF backup file "',trim(backup_ngwf_file),'" ... '
       ! open file
       open(unit=backup_ngwf_unit,file=backup_ngwf_file,iostat=io_status,&
            form='UNFORMATTED',action='READWRITE')
       call utils_open_unit_check('geom_BFGS','backup_ngwf_unit',io_status)

       ! delete file
       close(unit=backup_ngwf_unit,iostat=io_status,status='DELETE')
       call utils_close_unit_check('geom_BFGS','backup_ngwf_unit',io_status)

       ! success
       write(stdout,'(a)') 'done'
       ! ... then density kernel backup file...
       ! Find next available unit specifier
       backup_dk_unit = utils_unit()
       ! tell user what's going on
       write(stdout,'(3a)',advance='no') &
            ' Removing density kernel backup file "',trim(backup_dk_file),'" ... '
       ! open file
       open(unit=backup_dk_unit,file=backup_dk_file,iostat=io_status,&
            form='UNFORMATTED',action='READWRITE')
       call utils_open_unit_check('geom_BFGS','backup_dk_unit',io_status)

       ! delete file
       close(unit=backup_dk_unit,iostat=io_status,status='DELETE')
       call utils_close_unit_check('geom_BFGS','backup_dk_unit', io_status)

       ! success
       write(stdout,'(a)') 'done'
    endif
    ! =================================================================================
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<< End remove backup files >>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! =================================================================================


    deallocate (trial_f_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','trial_f_vec',ierr)

    deallocate (line_f_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','line_f_vec',ierr)

    deallocate (quad_f_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','quad_f_vec',ierr)

    deallocate(backup_forces,stat=ierr)
    call utils_dealloc_check('geom_BFGS','backup_forces',ierr)

    deallocate (backup_x_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','backup_x_vec',ierr)

    deallocate (backup_f_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','backup_f_vec',ierr)

    call castep_cell_dealloc(init_cell)

    deallocate (x_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','x_vec',ierr)

    deallocate (old_x_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','old_x_vec',ierr)

    deallocate (f_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','f_vec',ierr)

    deallocate (old_f_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','old_f_vec',ierr)

    deallocate (delta_vec,stat=ierr)
    call utils_dealloc_check('geom_BFGS','delta_vec',ierr)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_BFGS'
#endif

    !-------------------------------------------------------------------------!
    !Its the end of the BFGS code as we know it and I feel fine ...           !
    !-------------------------------------------------------------------------!

1   format(1x,a,i4,4f15.5)
2   format(1x,a,es17.8e3,1x,a)
4   format(/,80('='),/,1x,a,i4,a,/,80('='))
5   format(/,80('-'),/,1x,a,i4,a,f12.6,')',/,80('-'))
7   format(1x,a,f11.6)
10  format(1x,a,a4,4a15)
11  format(/,80('='),/,1x,a,/,80('='))

    return

  contains

    subroutine geom_BFGS_revert_old
      !=========================================================================!
      ! 'contain'ed subroutine to implement backtracking from current to        !
      ! previous configuration as done in several places and in same way        !
      !-------------------------------------------------------------------------!
      ! Arguments:                                                              !
      !   none - access all geom_BFGS variables directly as in same scope       !
      !-------------------------------------------------------------------------!
      ! Parent module variables used:                                           !
      !   same as geom_BFGS as in same scope                                    !
      !-------------------------------------------------------------------------!
      ! Modules used:                                                           !
      !   same as geom_BFGS as in same scope                                    !
      !-------------------------------------------------------------------------!
      ! Key Internal Variables:                                                 !
      !   reset_iter = value of mdl%bfgs_iteration last time we did this        !
      !                so can trap consecutive backtracks and abort if required !
      !   monotonic_E= can turn off if we get desperate                         !
      !   abort_BFGS = give up altogether                                       !
      !-------------------------------------------------------------------------!
      ! Necessary conditions:                                                   !
      !-------------------------------------------------------------------------!
      ! Written by Matt Probert, v1.0, 11/05/2002                               !
      !-------------------------------------------------------------------------!
      ! Modified for ONETEP by Arash A Mostofi, 2004                            !
      !=========================================================================!


      use constants, only: dp, stdout

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') &
           'DEBUG: Entering geom_BFGS_revert_old'
#endif

      !backtrack model
      x_vec=old_x_vec
      f_vec=old_f_vec
      enthalpy=old_enthalpy
      call geom_xvec_to_mdl(old_x_vec,mdl)
      call geom_fvec_to_mdl(old_f_vec,mdl)

      !NB geom_xvec_to_mdl calls model_cell_changed
      !=> invalidates model flags so we need to reset them here
      !=> but we have not restored the wavefunction etc


      !before we throw away the old inv_Hessian, lets try to extract some useful stuff out of it
      if (mdl%bfgs_iteration-reset_iter>1) call geom_update_inv_Hessian_params(init_cell,mdl)

      call castep_cell_copy(mdl%cell,init_cell)
      call geom_inv_Hessian_initialize(init_cell,mdl)

      !apply inv_Hessian to get new delta
      if (pub_on_root) then
         delta_vec=matmul(mdl%bfgs_inv_Hessian,f_vec)
      end if
      call comms_bcast(pub_root_node_id,delta_vec)

      !test convergence
      call geom_converged(mdl,old_x_vec,enthalpy,shuffle_backwards,dE,Fmax,dRmax,Smax, &
           & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax)

      !tell user what just happened
      call geom_output_converged(mdl%bfgs_iteration,enthalpy,dE,Fmax,dRmax,&
           converged_dE,converged_Fmax,converged_dRmax)

      !abort or continue?
      if (mdl%bfgs_iteration-reset_iter==1) then !whoops - we've reseting after a previous reset!
         !Test if we are able to toggle off-diag element in inv_Hessian on/off ..

            !We're in trouble ...
            if (pub_on_root) then
               write (stdout,*) 'BFGS: Warning - Repeated consecutive reset of inverse Hessian'
               write (stdout,*) 'BFGS:           without satisfying convergence criteria which'
               write (stdout,*) 'BFGS:           looks like BFGS has run out of search directions.'
            end if

            if (pub_on_root) write (stdout,*) 'BFGS: Warning - Terminating BFGS loop.'
            abort_BFGS=.true.

      else                            !not consecutive reset
         reset_iter=mdl%bfgs_iteration
      end if

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving geom_BFGS_revert_old'
#endif

      return
    end subroutine geom_BFGS_revert_old

    subroutine geom_BFGS_backup
      !=======================================================================!
      ! 'contain'ed subroutine to store key data to speed up backtracking     !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      !   none - access all geom_BFGS variables directly as in same scope     !
      !-----------------------------------------------------------------------!
      ! Parent module variables used:                                         !
      !   same as geom_BFGS as in same scope                                  !
      !-----------------------------------------------------------------------!
      ! Modules used:                                                         !
      !   same as geom_BFGS as in same scope                                  !
      !-----------------------------------------------------------------------!
      ! Key Internal Variables:                                               !
      !-----------------------------------------------------------------------!
      ! Necessary conditions:                                                 !
      !-----------------------------------------------------------------------!
      ! Written by Matt Probert, v1.0, 11/05/2002                             !
      !-----------------------------------------------------------------------!
      ! Modified for ONETEP by Arash A Mostofi, 2004                          !
      !=======================================================================!

      use comms, only: comms_abort, comms_barrier, pub_on_root
      use constants, only: DP, stderr
      use utils, only: utils_binary_copy

      implicit none

      integer :: read_unit
      integer :: atom,ngwf,dim
      integer :: num,tn1,tn2,tn3,tn(3)
#ifdef ACCELRYS
      integer :: row           ! vm: required for Intel bug workaround
#endif
      real(kind=DP) :: orig1,orig2,orig3
      real(kind=DP), allocatable, dimension(:,:,:) :: uni_tbox
      character(len=file_maxpath) :: read_file

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_BFGS_backup'
#endif

      !always backup mdl%cell in memory - small arrays so not an issue here
      call castep_cell_copy(mdl%cell,backup_cell)

      !backup energies
      backup_total_energy = mdl%total_energy
      backup_enthalpy = enthalpy

      if (pub_on_root) then

         !**********************************************************!
         ! AAM: read current tightbox_ngwfs and copy to backup file !
         !**********************************************************!

         ! Find next available unit specifier
         backup_ngwf_unit = utils_unit()

         ! Open backup file
         open(unit=backup_ngwf_unit,file=backup_ngwf_file,iostat=io_status, &
              form='UNFORMATTED',action='WRITE')
         call utils_open_unit_check('geom_BFGS_backup','backup_ngwf_file',io_status)
         rewind(backup_ngwf_unit)

         ! Find available unit for reading current NGWFs
         read_unit = utils_unit()

         ! Set up name of current NGWF file to be read
         write(read_file,'(2a)') trim(pub_rootname),'.tightbox_ngwfs'

         ! Tell user what is going on
         write(stdout,'(3a)',advance='no') ' Backing up NGWFs to file "', &
              trim(backup_ngwf_file),'" ...'

         ! Open current NGWF file to be read
         open(unit=read_unit,file=read_file,form='UNFORMATTED',status='OLD', &
              action='READ',iostat=io_status)
         call utils_open_unit_check('geom_BFGS_backup','read_file',io_status)
         rewind(read_unit)

         ! Read number of NGWFs
         read(read_unit,iostat=io_status) num
         call utils_read_check('geom_BFGS_backup','num',io_status)

         ! Sanity check
         if (num /= pub_cell%num_ngwfs) then
            write(stderr,'(2(a,i6),a)') 'Error in geom_BFGS_backup: num (', &
                 num,') not equal to pub_cell%num_ngwfs (',pub_cell%num_ngwfs,')'
            call comms_abort
         end if

         ! Write number of NGWFs to backup file
         write(backup_ngwf_unit,iostat=io_status) num
         call utils_write_check('geom_BFGS_backup','num',io_status)

         ! Loop over dimensions
         do dim=1,3

            ! Read number of points in dimension dim of universal tightbox
            read(read_unit,iostat=io_status) tn(dim)
            call utils_read_check('geom_BFGS_backup','tn',io_status)

            ! Write number of points in first dimension of universal tightbox
            write(backup_ngwf_unit,iostat=io_status) tn(dim)
            call utils_write_check('geom_BFGS_backup','backup_ngwf_unit',io_status)

         end do

         ! Allocate universal tightbox
         tn1 = tn(1) ; tn2 = tn(2) ; tn3 = tn(3)
         allocate(uni_tbox(tn1,tn2,tn3),stat=ierr)
         call utils_alloc_check('geom_BFGS_backup','uni_tbox',ierr)
         uni_tbox = 0.0_dp

         atom_loop: do atom=1,pub_cell%nat

            ngwfs_on_atom_loop: do ngwf=1,elements(atom)%nfunctions

               ! Read origin
               read(read_unit,iostat=io_status) orig1, orig2, orig3
               call utils_read_check('geom_BFGS_backup','orig1, orig2, orig3',io_status)

               ! Write origin
               write(backup_ngwf_unit,iostat=io_status) orig1, orig2, orig3
               call utils_write_check('geom_BFGS_backup','orig1, orig2, orig3',io_status)

               ! Read tightbox
               read(read_unit,iostat=io_status) uni_tbox(1:tn1,1:tn2,1:tn3)
               call utils_read_check('geom_BFGS_backup','uni_tbox',io_status)

               ! Write tightbox
               write(backup_ngwf_unit,iostat=io_status) uni_tbox(1:tn1,1:tn2,1:tn3)
               call utils_write_check('geom_BFGS_backup','uni_tbox',io_status)

            end do ngwfs_on_atom_loop

         end do atom_loop

         ! Close NGWF backup file
         close(backup_ngwf_unit,iostat=io_status)
         call utils_close_unit_check('geom_BFGS_backup','backup_ngwf_unit',io_status)

         ! Close current NGWF file
         close(read_unit,iostat=io_status)
         call utils_close_unit_check('geom_BFGS_backup','read_unit',io_status)

         write(stdout,'(a)') ' done'

         ! Deallocate memory
         deallocate(uni_tbox,stat=ierr)
         call utils_dealloc_check('geom_BFGS_backup','uni_tbox',ierr)

         !***************************************************************!
         ! AAM: read current density kernel from file and backup to file !
         !***************************************************************!

         write(stdout,'(3a)',advance='no') &
              ' Backing up density kernel to file "', &
              trim(backup_dk_file), '" ...'

         ! Find available unit specifier for density kernel backup
         backup_dk_unit = utils_unit()

         ! Open density kernel backup file
         open(unit=backup_dk_unit,file=backup_dk_file,iostat=io_status, &
              form='UNFORMATTED',action='WRITE')
         call utils_open_unit_check('geom_BFGS_backup','backup_dk_file',io_status)
         rewind(backup_dk_unit)

         ! Find available unit for current density kernel file to be read
         read_unit = utils_unit()

         ! Name of current density kernel file to be read
         write(read_file,'(2a)') trim(pub_rootname),'.dkn'

         ! Open current density kernel file
         open(unit=read_unit,file=read_file,form='UNFORMATTED', &
              status='OLD',action='READ',iostat=io_status)
         call utils_open_unit_check('geom_BFGS_backup','read_file',io_status)
         rewind(read_unit)

         ! Copy data from one file to the other
         call utils_binary_copy(read_unit,backup_dk_unit)

         ! Close density kernel backup file
         close(backup_dk_unit,iostat=io_status)
         call utils_close_unit_check('geom_BFGS_backup','backup_dk_unit',io_status)

         ! Close current density kernel file
         close(read_unit,iostat=io_status)
         call utils_close_unit_check('geom_BFGS_backup','read_unit',io_status)

         write(stdout,'(a)') ' done'

      end if

      !backup ionic stuff
      backup_forces = mdl%forces
      backup_stress = mdl%stress
      backup_strain = mdl%strain
      backup_x_vec = x_vec
      backup_f_vec = f_vec

      call comms_barrier

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_BFGS_backup'
#endif

      return
    end subroutine geom_BFGS_backup

    subroutine geom_BFGS_restore
      !=======================================================================!
      ! 'contain'ed subroutine to restore key data to speed up backtracking   !
      !-----------------------------------------------------------------------!
      ! Arguments:                                                            !
      !   none - access all geom_BFGS variables directly as in same scope     !
      !-----------------------------------------------------------------------!
      ! Parent module variables used:                                         !
      !   same as geom_BFGS as in same scope                                  !
      !-----------------------------------------------------------------------!
      ! Modules used:                                                         !
      !   same as geom_BFGS as in same scope                                  !
      !-----------------------------------------------------------------------!
      ! Key Internal Variables:                                               !
      !-----------------------------------------------------------------------!
      ! Necessary conditions:                                                 !
      !-----------------------------------------------------------------------!
      ! Written by Matt Probert, v1.0, 11/05/2002                             !
      !-----------------------------------------------------------------------!
      ! Modified for ONETEP by Arash A Mostofi, 2004                          !
      !=======================================================================!

      use constants, only: DP, stdout, stderr
      use utils, only: utils_binary_copy

      implicit none

      integer :: write_unit
      integer :: atom,ngwf,dim
      integer :: num,tn1,tn2,tn3,tn(3)
#ifdef ACCELRYS
      integer :: row           ! vm: required for Intel bug workaround
#endif
      real(kind=DP) :: orig1,orig2,orig3
      real(kind=DP), allocatable, dimension(:,:,:) :: uni_tbox
      character(len=file_maxpath) :: write_file

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_BFGS_restore'
#endif

      !restore mdl%cell from memory
      call castep_cell_copy(backup_cell,mdl%cell)
      !and make sure that current_cell is in sync
      call castep_cell_copy(mdl%cell,current_cell)

      ! AAM: bring orig_cell in sync too
      call castep_cell_copy(mdl%cell,mdl%orig_cell)

      !restore energies
      mdl%total_energy = backup_total_energy
      enthalpy = backup_enthalpy

      if (pub_on_root) then

         !*******************************************************************!
         !                AAM: start with tightbox_ngwfs                     !
         !*******************************************************************!

         ! Find available unit specifier for NGWF backup
         backup_ngwf_unit = utils_unit()

         ! Open backup file
         open(unit=backup_ngwf_unit,file=backup_ngwf_file,iostat=io_status, &
              form='UNFORMATTED',action='READ')
         call utils_open_unit_check('geom_BFGS_backup','backup_ngwf_file',io_status)
         rewind(backup_ngwf_unit)

         ! Tell user what's going on
         write(stdout,'(3a)',advance='no') &
              'Restoring NGWFs from file "', trim(backup_ngwf_file),'" ...'

         ! Find a free unit specifier
         write_unit = utils_unit()

         ! Filename to be written
         write(write_file,'(2a)') trim(pub_rootname),'.tightbox_ngwfs'

         ! Open file to be written
         open(unit=write_unit,file=write_file,form='UNFORMATTED', &
              action='WRITE',iostat=io_status)
         call utils_open_unit_check('geom_BFGS_backup','write_file',io_status)
         rewind(write_unit)

         ! Read number of NGWFs (backup_ngwf_unit already open)
         read(backup_ngwf_unit,iostat=io_status) num
         call utils_read_check('geom_BFGS_restore','num',io_status)

         if (num /= pub_cell%num_ngwfs) then
            write(stderr,'(2(a,i6),a)') 'Error in geom_BFGS_restore: num (', &
                 num,') not equal to pub_cell%num_ngwfs (',pub_cell%num_ngwfs,')'
            call comms_abort
         end if

         ! Write number of NGWFs to write_file
         write(write_unit,iostat=io_status) num
         call utils_write_check('geom_BFGS_restore','num',io_status)

         ! Loop over dimensions
         do dim=1,3

            ! Read number of points in dimension dim of universal tightbox
            read(backup_ngwf_unit,iostat=io_status) tn(dim)
            call utils_read_check('geom_BFGS_restore','tn',io_status)

            ! Write number of points in dimension dim of universal tightbox
            write(write_unit,iostat=io_status) tn(dim)
            call utils_write_check('geom_BFGS_restore','tn',io_status)

         end do

         ! Allocate universal tightbox
         tn1 = tn(1) ; tn2 = tn(2) ; tn3 = tn(3)
         allocate(uni_tbox(tn1,tn2,tn3),stat=ierr)
         call utils_alloc_check('geom_BFGS_restore','uni_tbox',ierr)
         uni_tbox = 0.0_dp

         atom_loop: do atom=1,pub_cell%nat

            ngwfs_on_atom_loop: do ngwf=1,elements(atom)%nfunctions

               ! Read origin
               read(backup_ngwf_unit,iostat=io_status) orig1, orig2, orig3
               call utils_read_check('geom_BFGS_restore','orig1, orig2, orig3',io_status)

               ! Write origin
               write(write_unit,iostat=io_status) orig1, orig2, orig3
               call utils_write_check('geom_BFGS_restore','orig1, orig2, orig3',io_status)

               ! Read tightbox
               read(backup_ngwf_unit,iostat=io_status) uni_tbox(1:tn1,1:tn2,1:tn3)
               call utils_read_check('geom_BFGS_restore','uni_tbox',io_status)

               ! Write tightbox
               write(write_unit,iostat=io_status) uni_tbox(1:tn1,1:tn2,1:tn3)
               call utils_write_check('geom_BFGS_restore','uni_tbox',io_status)

            end do ngwfs_on_atom_loop

         end do atom_loop

         close(backup_ngwf_unit,iostat=io_status)
         call utils_close_unit_check('geom_BFGS_restore','backup_ngwf_unit',io_status)

         close(write_unit,iostat=io_status)
         call utils_close_unit_check('geom_BFGS_restore','write_unit',io_status)

         write(stdout,'(a)') ' done'

         deallocate(uni_tbox,stat=ierr)
         call utils_dealloc_check('geom_BFGS_restore','uni_tbox',ierr)

         !********************************************************************!
         !        AAM: read density kernel from backup file and restore       !
         !********************************************************************!

         ! Name of current density kernel file to be written
         write(write_file,'(2a)') trim(pub_rootname),'.dkn'
         write(stdout,'(3a)',advance='no') &
              'Restoring density kernel to file "', trim(write_file),'" ...'

         ! Find available unit specifier for density kernel backup
         backup_dk_unit = utils_unit()

         ! Open density kernel backup file
         open(unit=backup_dk_unit,file=backup_dk_file,iostat=io_status, &
              form='UNFORMATTED',action='READ')
         call utils_open_unit_check('geom_BFGS_restore','backup_dk_file',io_status)
         rewind(backup_dk_unit)

         ! Find available unit
         write_unit = utils_unit()

         ! Open file
         open(unit=write_unit,file=write_file,form='UNFORMATTED', &
              action='WRITE',iostat=io_status)
         call utils_open_unit_check('geom_BFGS_restore','write_file',io_status)
         rewind(write_unit)

         ! Copy data from one file to the other
         call utils_binary_copy(backup_dk_unit,write_unit)

         close(backup_dk_unit,iostat=io_status)
         call utils_close_unit_check('geom_BFGS_restore','backup_dk_unit',io_status)

         close(write_unit,iostat=io_status)
         call utils_close_unit_check('geom_BFGS_restore','write_unit',io_status)

         write(stdout,'(a)') ' done'

      end if

      !restore ionic stuff
      mdl%forces = backup_forces
      mdl%stress = backup_stress
      mdl%strain = backup_strain
      x_vec = backup_x_vec
      f_vec = backup_f_vec

      !restore model flags
      mdl%found_ground_state = .true.

      !NB No need to call wave_Sorthonormalise as wvfn was Sorthonormal when stored

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_BFGS_restore'
#endif

      return
    end subroutine geom_BFGS_restore

    subroutine geom_BFGS_write_summary()
      !=========================================================================!
      ! contain'ed subroutine to write out the BFGS minimisation summary.       !
      !-------------------------------------------------------------------------!
      ! Arguments:                                                              !
      !    none                                                                 !
      !-------------------------------------------------------------------------!
      ! Parent module variables used:                                           !
      !    none                                                                 !
      !-------------------------------------------------------------------------!
      ! Modules used:                                                           !
      !    none                                                                 !
      !-------------------------------------------------------------------------!
      ! Key Internal Variables:                                                 !
      !    none                                                                 !
      !-------------------------------------------------------------------------!
      ! Necessary conditions:                                                   !
      !-------------------------------------------------------------------------!
      ! Written by Matt Probert, v1.0, 11/05/2002                               !
      !-------------------------------------------------------------------------!
      ! Modified for ONETEP by Arash A Mostofi, 2004                            !
      !=========================================================================!


      use constants, only: dp, stdout
      use services, only: services_flush

      !local vars
      character(len=80) :: data_string
      character(len=80) :: divider_string
      character(len=80) :: label_string
      integer           :: string_index
      integer           :: len_label, len_lambda, len_fdelta, len_enthalpy

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') &
           'DEBUG: Entering geom_BFGS_write_summary'
#endif

      !initialise strings
      data_string     = repeat(' ',len(data_string))
      divider_string  = repeat(' ',len(divider_string))
      label_string    = repeat(' ',len(label_string))
      !                  12345678901234 2345678901234 2345678901234 23456789012345678
      divider_string  = '+------------+-------------+-------------+-----------------+'
      label_string    = '|    Step    |   lambda    |   F.delta   |    enthalpy     |'
      !e.g.              | trial step    -100.12345    -100.12345    -123456.123456

      if(pub_on_root) then

         !header
         write(stdout,*) ' '
         write(stdout,5) divider_string
         write(stdout,5) label_string
         write(stdout,5) divider_string

         !previous                                                    '123456789012'
         string_index=1
         len_label   =14
         len_lambda  =14
         len_fdelta  =14
         len_enthalpy=18
         write(data_string(string_index:string_index+len_label),   1) '  previous  '
         string_index=string_index+len_label
         write(data_string(string_index:string_index+len_lambda),  2) 0.0_dp
         string_index=string_index+len_lambda
         write(data_string(string_index:string_index+len_fdelta),  3) old_f_dot_delta
         string_index=string_index+len_fdelta
         write(data_string(string_index:string_index+len_enthalpy),4) old_enthalpy
         string_index=string_index+len_enthalpy
         !cross check all is well (NB strings start from 1 not 0 ...)
         if (string_index>71) write (stdout,*) 'geom_BFGS_write_summary: format problem?'
         write(stdout,5) data_string

         !trial                                                          '123456789012'
         if (done_trial) then
            string_index=1
            write(data_string(string_index:string_index+len_label),   1) ' trial step '
            string_index=string_index+len_label
            write(data_string(string_index:string_index+len_lambda),  2) lambda_trial
            string_index=string_index+len_lambda
            write(data_string(string_index:string_index+len_fdelta),  3) trial_f_dot_delta
            string_index=string_index+len_fdelta
            write(data_string(string_index:string_index+len_enthalpy),4) trial_enthalpy
            string_index=string_index+len_enthalpy
            write(stdout,5) data_string
         end if

         !line                                                           '123456789012'
         if (done_line) then
            string_index=1
            write(data_string(string_index:string_index+len_label),   1) '  line step '
            string_index=string_index+len_label
            write(data_string(string_index:string_index+len_lambda),  2) lambda_line
            string_index=string_index+len_lambda
            write(data_string(string_index:string_index+len_fdelta),  3) line_f_dot_delta
            string_index=string_index+len_fdelta
            write(data_string(string_index:string_index+len_enthalpy),4) line_enthalpy
            string_index=string_index+len_enthalpy
            write(stdout,5) data_string
         end if

         !quad                                                           '123456789012'
         if (done_quad) then
            string_index=1
            write(data_string(string_index:string_index+len_label),   1) '  quad step '
            string_index=string_index+len_label
            write(data_string(string_index:string_index+len_lambda),  2) lambda_quad
            string_index=string_index+len_lambda
            write(data_string(string_index:string_index+len_fdelta),  3) quad_f_dot_delta
            string_index=string_index+len_fdelta
            write(data_string(string_index:string_index+len_enthalpy),4) quad_enthalpy
            string_index=string_index+len_enthalpy
            write(stdout,5) data_string
         end if

         !footer
         write(stdout,5) divider_string
         write(stdout,*) ' '

      end if          !pub_on_root
      call services_flush

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') &
           'DEBUG: Leaving geom_BFGS_write_summary'
#endif

1     format('|',a11,1x,'|')        !1+11+1+1=14
2     format(1x,f11.5,1x,'|')       !1+11+1+1=14
3     format(1x,f11.6,1x,'|')       !1+11+1+1=14
4     format(1x,f15.6,1x,'|')       !1+15+1+1=18
      !14+14+14+18=60
5     format(1x,a60,' <-- min BFGS')

    end subroutine geom_BFGS_write_summary
  end subroutine geom_BFGS

  subroutine geom_DELOC(mdl,elements,output_file,status)
    !=========================================================================!
    ! Find ground state ionic positions for given model using delocalized     !
    ! internal coordinates. This technology is based on the DMol3             !
    ! implementation which is described in:                                   !
    ! J. Andzelm, R. D. King-Smith, G. Fitzgerald,                            !
    ! Chem. Phys. Lett. 335 (2001) 321-326                                    !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mdl, intent=inout, the model to be optimised                          !
    !   output_file, intent=in, the filename to which results will be written.!
    !   status, intent=out                                                    !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   model                                                                 !
    !   parameters                                                            !
    !   comms                                                                 !
    !   geometry_deloc_utils                                                  !
    !   geometry_deloc_algor                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   mdl has already been properly read and initialised                    !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !=========================================================================!

    use simulation_cell, only: castep_model, pub_cell, castep_cell_nullify, &
             castep_cell_alloc, castep_cell_copy, castep_cell_dealloc
    use comms, only: comms_bcast, pub_root_node_id, pub_on_root
    use constants, only: dp, stdout
    use geometry_deloc_utils, only: deloc_utils_io_abort,deloc_utils_initialize,&
         deloc_utils_mdl_to_internals, deloc_utils_update_params, &
         deloc_utils_energy_gradient, deloc_utils_output_converged, &
         deloc_utils_internals_to_mdl, deloc_utils_deallocate,deloc_utils_nullify
    use geometry_deloc_algor, only: geom_algor_deloc
    use ion, only: element
    use rundat, only: geom_backup_iter, pub_output_detail, &
         geom_frequency_est, geom_modulus_est, geom_max_iter, &
         pub_write_positions
    use services, only: services_flush
    use utils, only: utils_alloc_check, utils_dealloc_check


    implicit none

    type(castep_model), intent(inout)               :: mdl
    type(element),      intent(inout)               :: elements(1:pub_cell%nat)
    character(len=*), intent(in)                    :: output_file ! trajectory
    integer,intent(out)                             :: status
    !-------------------------------------------------------------------------!

    type (castep_cell) :: init_cell !the current_cell when last initialised inv_Hessian
    integer :: ierr,fake_iter
    logical :: converged

    real(kind=DP) :: enthalpy
    real(kind=DP) :: dE
    real(kind=DP) :: Fmax   !max Force
    real(kind=DP) :: dRmax  !max displacement
    real(kind=DP) :: Smax   !max Stress
    logical       :: first_step=.true.
    logical       :: converged_dE, converged_Fmax, converged_dRmax, converged_Smax
    real(kind=DP), allocatable, dimension(:)     :: old_x_vec !old strains and fractional coords
    !only to be able to use geometry_converged
    !-------------------------------------------------------------------------!
    status = 0

    old_frequency_est=geom_frequency_est             !we store the starting values as 'old' to test o/p at end
    new_frequency_est=geom_frequency_est             !we store the working value as 'new' s.t. can cope with
    new_modulus_est  =geom_modulus_est               !automatic and user updates to the estimates

    allocate(old_x_vec(1:ndim),stat=ierr)
    call utils_alloc_check('geom_DELOC','old_x_vec',ierr)

    if (pub_on_root) write(stdout,'(/a,a,a)') &
         ' DI:   continuation file name is "',trim(cont_file),'"'

    !NB we properly nullify any pointers before first usage

    call castep_cell_nullify(init_cell)


    call castep_cell_alloc(init_cell)

    call castep_cell_copy(mdl%cell,init_cell)

    !-------------------------------------------------------------------------!
    !Sort out start from restart                                              !
    !-------------------------------------------------------------------------!

    !If this is not a restart, then need to initialise the inv_Hessian
    if (mdl%bfgs_iteration==0) then

       call geom_inv_Hessian_initialize(init_cell,mdl)
       mdl%strain=0.0_dp

       call geom_get_forces(mdl,elements)

       enthalpy = mdl%total_energy
       !start up the output file
       call geom_write_trajectory(mdl,enthalpy,output_file)

    else

       if (pub_on_root)write (stdout,*) 'Delocalized Internals: Restarting previous geometry optimization.'

       !Catch the possibility that model_continuation/reuse has decided to nuke the inv_Hessian
       if (all(abs(mdl%bfgs_inv_Hessian)<tiny(0.0_dp))) then
          call castep_cell_copy(mdl%cell,init_cell)
          call geom_inv_Hessian_initialize(init_cell,mdl)
       end if

    end if

    ! initialize geometry_deloc_utils module
    call deloc_utils_initialize(mdl,pub_output_detail)

    !NB We reset bfgs_iteration here as the geom_max_iter refers to number of iterations
    !   THIS pass and not in total.
    mdl%bfgs_iteration=0
    reset_iter=mdl%bfgs_iteration

    ! fill symmetry information
    call deloc_utils_mdl_to_internals(mdl,elements)

    ! generate initial delocalized coordinates
    fake_iter = 0
    if (pub_on_root) then
       call geom_algor_DELOC(.true.,fake_iter,status)
    end if
    call comms_bcast(pub_root_node_id,status)
    if (status/=0) go to 999

    !Put the strains wrt original cell and fractional atomic coords wrt current_cell into old_x_vec
    call geom_mdl_to_xvec(mdl,old_x_vec)

    main_loop : do

       ! check to see if we have done enough iterations yet?
       !NB Test is done this way around s.t. we don't need to undo anything at end
       if (mdl%bfgs_iteration>=geom_max_iter) exit main_loop

       !checkpoint and update as appropriate
       call deloc_utils_update_params

       ! get energy and gradients
       call deloc_utils_energy_gradient                ! SCF energy/gradient

       if (first_step) then
          !initialise convergence tests (impossible to converge here with convergence windows)
          call geom_converged(mdl,old_x_vec,enthalpy,shuffle_forwards,dE,Fmax,dRmax,Smax, &
               & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax)

          !tell user where we are starting from
          call deloc_utils_output_converged(mdl%bfgs_iteration,enthalpy,dE,Fmax,dRmax, &
                converged_dE,converged_Fmax,converged_dRmax)

          first_step = .false.

       endif

       mdl%bfgs_iteration = mdl%bfgs_iteration + 1             ! next geometry cycle

       !banner to output
       if (pub_on_root) write (stdout,4) 'Starting DI iteration ',mdl%bfgs_iteration,' ...'

       !store the current model?

       !Put the strains wrt original cell and fractional atomic coords wrt current_cell into old_x_vec
       call geom_mdl_to_xvec(mdl,old_x_vec)

       call deloc_utils_mdl_to_internals(mdl,elements)

       !info on current_cell to stdout
       ! make a step of geometry optimization

       if (pub_on_root) then
          ! Do all the things on the root node (i.e. look for  new coordinates):
          call geom_algor_DELOC(.false.,mdl%bfgs_iteration,status)
       end if
       call comms_bcast(pub_root_node_id,status)
       if (status/=0) go to 999

       ! update model data
       call deloc_utils_internals_to_mdl(mdl)
       call castep_cell_copy(mdl%cell,current_cell)

       !tell user what's happening
       if (pub_on_root) write (stdout,5) 'DI: structure has been updated on iteration',mdl%bfgs_iteration

       !tell user what's the current structure
       if (pub_write_positions) call geom_output(mdl)

       ! find the forces and energy
       call geom_get_forces(mdl,elements)

       call services_flush

       ! save trajectory data
       enthalpy = mdl%total_energy
       call geom_write_trajectory(mdl,enthalpy,output_file)

       !test convergence
       call geom_converged(mdl,old_x_vec,enthalpy,shuffle_forwards,dE,Fmax,dRmax,Smax, &
            & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax)

       !tell user what just happened (stop iterations if converged)
       call deloc_utils_output_converged(mdl%bfgs_iteration,enthalpy,dE,Fmax,dRmax, &
             converged_dE,converged_Fmax,converged_dRmax)

       if (converged) exit main_loop


       !checkpoint and update as appropriate
       call deloc_utils_update_params
       if (mod(mdl%bfgs_iteration,geom_backup_iter)==0) &
            call geom_opt_continuation_write(mdl,elements,cont_file)

    enddo main_loop

    if (pub_on_root) then
       if (converged) then
          write(stdout,*) 'DI: Geometry optimization completed successfully.'
       else
          write(stdout,*) 'DI: Geometry optimization failed to converge after',mdl%bfgs_iteration,'steps.'
       end if
       write (stdout,11) 'DI: Final Configuration:'
       call geom_output(mdl)
       write (stdout,2) 'DI: Final Enthalpy     =',enthalpy,trim(energy_label)
    end if


999 continue

    call deloc_utils_deallocate
    call deloc_utils_nullify

    deallocate(old_x_vec,stat=ierr)
    call utils_dealloc_check('geom_DELOC','old_x_vec',ierr)

    call castep_cell_dealloc(init_cell)

2   format(1x,a,es17.8e3,1x,a)
4   format(/,80('='),/,1x,a,i4,a,/,80('='))
5   format(/,80('-'),/,1x,a,i4,/,80('-'))
11  format(/,80('='),/,1x,a,/,80('='))

    return

  end subroutine geom_DELOC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine geom_xvec_to_mdl(x_vec,mdl)
    !=========================================================================!
    ! Updates mdl and current_cell for current strains and fractional ionic   !
    ! positions in x_vec                                                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   x_vec, intent=in, the ndim vector of search parameters                !
    !   mdl,intent=inout, the model to be updated for new cell & coords       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !   fixed_npw to see which kind of basis re-initialisation is required    !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   model for type definitions                                            !
    !   cell to update components of cell                                     !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout, two_pi
    use simulation_cell, only: castep_model, castep_cell_copy, &
         castep_model_cell_changed

    implicit none

    real(kind=DP), dimension(1:ndim), intent(in)  :: x_vec
    type(castep_model), intent(inout)             :: mdl

    !local variables
    integer                                       :: i,j,iatom,ispec

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_xvec_to_mdl'
#endif

    !Only need to update the fractional atomic coords if allowed to move
    if (.not.fix_all_ions) then
       !NB We no longer apply PBCs here and so mdl%cell is allowed to become
       !   unrationalized, i.e. fractional coords>1 etc.
       i=10
       do ispec=1,mdl%cell%num_species
          do iatom=1,mdl%cell%num_ions_in_species(ispec)
             do j=1,3
                mdl%cell%ionic_positions(j,iatom,ispec)=x_vec(i+j-1)
             end do
             i=i+3
          end do
       end do
    end if

    !Not much to do for the basis here:
    !Just synchronize mdl%cell with current_cell ...
    call castep_cell_copy(mdl%cell,current_cell)

    !... and update model flags and S-orthogonalise wvfn etc
    call castep_model_cell_changed(mdl)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_xvec_to_mdl'
#endif

    return
  end subroutine geom_xvec_to_mdl

  subroutine geom_fvec_to_mdl(f_vec,mdl)
    !=========================================================================!
    ! Updates mdl and current_cell for current stresses and fractional forces !
    ! positions in f_vec                                                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   f_vec, intent=in, the ndim vector of fractional forces                !
    !   mdl,intent=inout, the model to be updated                             !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   model for type definitions                                            !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: dp, stdout
    use simulation_cell, only: castep_model, castep_cell_frac_to_cart
    use wrappers, only: wrappers_invert_sym_matrix

    implicit none

    real(kind=DP), dimension(1:ndim), intent(in)  :: f_vec
    type(castep_model), intent(inout)             :: mdl

    !local variables
    real(kind=DP), dimension(1:3,1:3)             :: metric, inv_metric
    real(kind=DP), dimension(1:3)                 :: frac_forces
    integer                                       :: i,iatom,ispec

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_fvec_to_mdl'
#endif

    !NB WE REQUIRE that geom_xvec_to_mdl has been called immediately before this routine
    !   s.t. mdl%cell and mdl%strain are in sync with x_vec

    !calculate metric in same form as inv_Hessian, metric=h^T.h
    metric=matmul(mdl%cell%real_lattice,transpose(mdl%cell%real_lattice))
    inv_metric=metric ; call wrappers_invert_sym_matrix(inv_metric,3)

    if (.not.fix_all_ions) then
       i=10
       do ispec=1,mdl%cell%num_species
          do iatom=1,mdl%cell%num_ions_in_species(ispec)

             !convert forces from fractional form
             frac_forces=matmul(inv_metric,f_vec(i:i+2))

             call castep_cell_frac_to_cart(current_cell,frac_forces,mdl%forces(:,iatom,ispec))

             !next ion
             i=i+3
          end do
       end do
    end if

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_fvec_to_mdl'
#endif

    return
  end subroutine geom_fvec_to_mdl

  subroutine geom_mdl_to_xvec(mdl,x_vec)
    !=========================================================================!
    ! Extract strains wrt original cell and fractional ionic coords from mdl  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mdl,   intent=in, the model to extract strain & coords from           !
    !   x_vec, intent=out,the ndim vector of search parameters                !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   io for error messages                                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use constants, only: dp, stdout, NORMAL
    use comms, only: pub_on_root
    use rundat, only: pub_output_detail
    use simulation_cell, only: castep_model

    implicit none

    type(castep_model), intent(in)                :: mdl
    real(kind=DP), dimension(1:ndim), intent(out) :: x_vec

    !local variables
    integer                                       :: i,j,iatom,ispec

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_mdl_to_xvec'
#endif

    !Pack strain into x_vec using same order as stresses into f_vec
    x_vec(1:9)=reshape(mdl%strain,(/9/))

    !... and pack the fractional atomic coords into x_vec
    i=10



    do ispec=1,mdl%cell%num_species
       do iatom=1,mdl%cell%num_ions_in_species(ispec)
          do j=1,3
             x_vec(i+j-1)=mdl%cell%ionic_positions(j,iatom,ispec)
          end do
          i=i+3
       end do
    end do

    if (pub_output_detail>NORMAL.and.pub_on_root) then
       do i=1,ndim
          write (stdout,99) i,x_vec(i)
       end do
    end if
99  format('x_vec(',i2,')=',f10.5)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_mdl_to_xvec'
#endif

    return
  end subroutine geom_mdl_to_xvec

  subroutine geom_mdl_to_fvec(mdl,f_vec,enthalpy)
    !=========================================================================!
    ! Extract stresses and forces from mdl                                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mdl,   intent=in, the model to extract strain & coords from           !
    !   f_vec,         intent=out,the ndim vector of gradients                !
    !   enthalpy,      intent=out,the enthalpy=Etot+pV of the system          !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   cell for external pressure and cart_to_frac conversion                !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   stress = total unsymmetrised stress, including external pressure      !
    !   strain = difference between current cell & original cell              !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: dp, stdout, NORMAL, two_pi
    use rundat, only: pub_output_detail
    use simulation_cell, only: castep_model

    implicit none

    type(castep_model), intent(in)                :: mdl
    real(kind=DP), dimension(1:ndim), intent(out) :: f_vec
    real(kind=DP),    intent(out)                 :: enthalpy


    !local variables
    real(kind=DP), dimension(1:3,1:3)             :: metric
    real(kind=DP), dimension(1:3)                 :: frac_forces
    real(kind=DP)                                 :: pressure
    real(kind=DP)                                 :: shear_energy

    integer                                       :: i,iatom,ispec

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_mdl_to_fvec'
#endif

    pressure=(external_pressure(1)+external_pressure(2)+external_pressure(3))/3.0_dp

    !and now include the ionic forces (Pfrommer eq 5) ...
    !calculate metric in same form as inv_Hessian, metric=h^T.h
    metric=matmul(mdl%cell%real_lattice,transpose(mdl%cell%real_lattice))

    ! jd: Takes care of the part of f_vec that would otherwise remain
    !     uninitialized. Since inv_Hessian is not necessarily diagonal,
    !     leaving garbage in f_vec(1:9) is risky.
    f_vec(1:9) = 0.0_DP

    if (.not.fix_all_ions) then
       i=10
       do ispec=1,mdl%cell%num_species
          do iatom=1,mdl%cell%num_ions_in_species(ispec)

             !convert forces from model to fractional form
             !NB we DO NOT use cell_cart_to_frac as this imposes PBCs on frac
             !   which is disastrous for forces!!!
             frac_forces=matmul(mdl%cell%recip_lattice,mdl%forces(:,iatom,ispec))/two_pi

             !store metric*fractional_forces in f_vec
             f_vec(i:i+2)=matmul(metric,frac_forces)

             !next ion
             i=i+3
          end do
       end do
    end if

    enthalpy=mdl%total_energy+pressure*mdl%cell%volume

    !Calculate elastic shear energy contribution (off diagonal pressure) according to Kittel (3rd ed) eq 4.14
    !NB this requires strain wrt EQUILIBRIUM => only correct if ref_cell is equilibrated ...
    !NB if external pressure is shear, then we ought to have no symmetry ops and so cannot presume
    !   that strain is symmetric => cannot do factor of 0.5 & 2 cancellation ...
    shear_energy=0.0_dp
    !pdh: remove any() structure which upsets gfortran
    if (abs(external_pressure(4))>epsilon(1.0_dp).or.abs(external_pressure(5))>epsilon(1.0_dp).or. &
         &abs(external_pressure(6))>epsilon(1.0_dp)) then
       shear_energy=shear_energy+(mdl%stress(4)+external_pressure(4))*(mdl%strain(2,3)+mdl%strain(3,2))
       shear_energy=shear_energy+(mdl%stress(5)+external_pressure(5))*(mdl%strain(3,1)+mdl%strain(1,3))
       shear_energy=shear_energy+(mdl%stress(6)+external_pressure(6))*(mdl%strain(1,2)+mdl%strain(2,1))
       shear_energy=0.5_dp*mdl%cell%volume*shear_energy
       enthalpy=enthalpy+shear_energy
    end if

    if (pub_output_detail>NORMAL.and.pub_on_root) then
       write (stdout,*)  'Enthalpy contributions:'
       write (stdout,98) 'total energy',mdl%total_energy,trim(energy_label)
       write (stdout,98) 'diagonal pressure term',pressure*mdl%cell%volume,trim(energy_label)
       write (stdout,98) 'off-diagonal pressure term',shear_energy,trim(energy_label)
       do i=1,ndim
          write (stdout,99) i,f_vec(i)
       end do
    end if

98  format(a30,'=',1x,es17.8e3,1x,a)
99  format('f_vec(',i2,')=',f10.5)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_mdl_to_fvec'
#endif

    return
  end subroutine geom_mdl_to_fvec

  subroutine geom_check_lambda(lambda,delta_vec)
    !=========================================================================!
    ! Check proposed value of lambda does not lead to too large a displacement!
    ! and if it does, reduce it until tolerance satisfied.                    !
    ! Large values of lambda in themselves are not a problem - only if cause  !
    ! an unreasonably large change in cell or ionic coords.                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   lambda, intent=inout, the lambda to be tested/updated                 !
    !   delta_vec, intent=in, the displacement vector to check against        !
    !   old_x_vec, intent=in, the reference structure to check delta against  !
    !   mdl, intent=in, the model for the proposed cell vector change         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   on_root for diagnostics                                               !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: dp, stdout, NORMAL
    use rundat, only: pub_output_detail
    use simulation_cell, only: castep_cell_frac_to_cart

    implicit none
    real(kind=DP), intent(inout)                 :: lambda
    real(kind=DP), intent(in), dimension(1:ndim) :: delta_vec

    !local vars
    real(kind=DP), parameter                     :: delta_ions_tol=0.5_dp !max allowed mag. change in Bohr (au)
    real(kind=DP)                                :: lambda_new
    real(kind=DP)                                :: delta_mag
    real(kind=DP), dimension(1:3)                :: delta_frac, delta_cart
    integer                                      :: i,j
    logical                                      :: warn_ions

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_check_lambda'
#endif

    !convert lambda*delta into 3D ionic change
    if (.not.fix_all_ions) then

       delta_mag=-1.0_dp
       do i=10,ndim,3     !check appropriate displacements
          do j=1,3
             delta_frac(j)=lambda*delta_vec(i+j-1)  !this is the proposed displacement wrt that at end of last iter
          end do                                    !NOT wrt last minimization !?

          !convert from fractional to cartesian
          call castep_cell_frac_to_cart(current_cell,delta_frac,delta_cart)

          !compare cartesian length to max allowed displacement (delta_ions_tol)
          delta_mag=max(delta_mag,sqrt(dot_product(delta_cart,delta_cart)))
       end do

       if (delta_mag>delta_ions_tol) then
          !NB delta_mag >=0 but lambda may (if .not.monotonic_E) be < 0 so must be careful here
          lambda_new=lambda*delta_ions_tol/delta_mag
          if (pub_output_detail>NORMAL.and.pub_on_root) then
             write (stdout,1) 'geom_check_lambda: Warning - BFGS step would result in ionic displacement >',delta_ions_tol,' au'
             write (stdout,2) 'geom_check_lambda: Warning - reducing lambda from ',lambda,' to ',lambda_new,' for safety'
          end if
          lambda=lambda_new
       end if

       !warn user if overall ionic coords have changed by a lot
       warn_ions=.false.
       do i=10,ndim,3     !check appropriate displacements
          do j=1,3
             delta_frac(j)=lambda*delta_vec(i+j-1)  !this is the checked & scaled displacement
          end do
          if (maxval(delta_frac)>0.1_dp) warn_ions=.true.
       end do
       if (warn_ions .and. pub_on_root) then
          write (stdout,*) 'geom_check_lambda: Warning - fractional ionic displacement > 10%'
          write (stdout,*) 'geom_check_lambda:         - suggest you reconsider the initial ionic positions'
       end if

    end if

1   format(1x,a,f10.3,a)
2   format(1x,a,f10.6,a,f10.6,a)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_check_lambda'
#endif

    return
  end subroutine geom_check_lambda

  subroutine geom_inv_Hessian_initialize(init_cell,mdl)
    !=========================================================================!
    ! Initialize the inverse Hessian matrix for good convergence to gs ions   !
    ! as in B.G. Pfrommer et al, J.Comp.Phys. v131, p233-240 (1997).          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   toggle,    intent=in,  switch to determine symm break or not          !
    !   init_cell, intent=in,  the cell to derive inv_Hessian from            !
    !   mdl,     intent=inout, model containing suitably initialized inv_Hess !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   parameters to define the bulk modulus and average phonon frequency    !
    !   io for error messages and unit conversion routines                    !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   init_cell must be valid                                               !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: dp, stdout, NORMAL
    use rundat, only: geom_print_inv_hessian, pub_output_detail
    use simulation_cell, only: castep_cell, castep_model
    use wrappers, only: wrappers_invert_sym_matrix

    implicit none

    type(castep_cell), intent(in)        :: init_cell
    type(castep_model), intent(inout)    :: mdl

    !local variables
    real(kind=DP), dimension(1:3,1:3) :: metric
    real(kind=DP) :: avg_mass
    real(kind=DP) :: inv_3VB
    integer :: i,j,k,l,ispec

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering geom_inv_Hessian_initialize'
#endif

    !Initialize ionic part of inv_Hessian using average phonon frequency
    !3x3 blocks=(metric^-1)/(mass*omega^2), where metric=h^T.h
    metric=matmul(init_cell%real_lattice,transpose(init_cell%real_lattice))
    call wrappers_invert_sym_matrix(metric,3)

    if (.not. pub_on_root) return

    !Initialize most of the matrix with a lot of zeros
    mdl%bfgs_inv_Hessian=0.0_dp

    !Initialize stress/strain part of inv_Hessian using bulk modulus
    !first 9 diagonal elements=1/(3*Vol*B)
    inv_3VB=1.0_dp/(3.0_dp*init_cell%volume*new_modulus_est)
    do i=1,9
       mdl%bfgs_inv_Hessian(i,i)=inv_3VB
    end do

    if (pub_output_detail>NORMAL .and. pub_on_root) then
       write(stdout,*) 'inv_3VB=           ',inv_3VB
       write(stdout,*) 'init_cell%volume=  ',init_cell%volume
       write(stdout,*) 'new_modulus_est=   ',new_modulus_est
       write(stdout,*) 'new_frequency_est= ',new_frequency_est
    endif


    avg_mass=0.0_dp
    do ispec=1,init_cell%num_species
       avg_mass=avg_mass+init_cell%species_mass(ispec)*init_cell%num_ions_in_species(ispec)
    end do

    ! ndmh: 26/01/09 very important fix: atomic units of mass are m_e so
    ! ndmh: multiply by (1 a.m.u.) / (1 m_e) = 1822.888485
    avg_mass=1822.888485_dp*avg_mass/init_cell%num_ions

    if (pub_output_detail>NORMAL .and. pub_on_root) write(stdout,*) 'avg_mass= ',avg_mass

    metric=metric/(avg_mass*new_frequency_est**2) !(metric^-1)/(mass*omega^2)

    do i=1,init_cell%num_ions
       k=10+3*(i-1)
       mdl%bfgs_inv_Hessian(k:k+2,k:k+2)=metric(1:3,1:3)
    end do

    !Finally, if symmetry is not to be strictly enforced, then we may add terms
    !to mix forces & stresses so ensure symmetry breaking is possible if required
    !NB We add these terms asymmetrically - we want to break symm for ions but not cell!


    if (geom_print_inv_hessian.and.pub_on_root) then
       l=ndim/5
       do i=1,ndim
          do j=1,5*l,5
             write (stdout,99) (i,j+k,mdl%bfgs_inv_Hessian(i,j+k:j+k),k=0,4)
          end do
          if (5*l<ndim) write (stdout,99) (i,j,mdl%bfgs_inv_Hessian(i,j:j),j=5*l,ndim)
       end do
    end if

99  format('inv_Hessian:',5('(',i2,',',i2,')=',f15.10,:))

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving geom_inv_Hessian_initialize'
#endif

    return
  end subroutine geom_inv_Hessian_initialize

  subroutine geom_update_inv_Hessian(old_x_vec,x_vec,old_f_vec,f_vec,init_cell,mdl)
    !=========================================================================!
    ! Update the inverse Hessian matrix using the standard BFGS algorithm     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   inv_Hessian, intent=inout, BFGS updated (Hessian)^-1                  !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !   reset_iter and previous_reset_iter may be updated                     !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use comms, only: pub_on_root, comms_abort
    use constants, only: dp, stdout, NORMAL
    use rundat, only: geom_print_inv_hessian, pub_output_detail
    use simulation_cell, only: castep_cell, castep_model, castep_cell_copy

    implicit none

    real(kind=DP), dimension(1:ndim),        intent(in)    :: old_x_vec
    real(kind=DP), dimension(1:ndim),        intent(in)    :: x_vec
    real(kind=DP), dimension(1:ndim),        intent(in)    :: old_f_vec
    real(kind=DP), dimension(1:ndim),        intent(in)    :: f_vec
    type(castep_cell),                       intent(inout) :: init_cell
    type(castep_model),                      intent(inout) :: mdl

    !local variables
    real(kind=DP), dimension(1:ndim) :: dx_vec, df_vec, u_vec, h_df_vec
    real(kind=DP)                    :: df_dot_h_df, dx_dot_df
    integer                          :: i,j,k,l

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering geom_update_inv_Hessian'
#endif

    if (.not.pub_on_root) return

    dx_vec=x_vec-old_x_vec
    df_vec=f_vec-old_f_vec
    dx_dot_df=dot_product(dx_vec,df_vec)

    !diagnostic
    if (pub_on_root.and.pub_output_detail>NORMAL) then
       if (dx_dot_df>tiny(dx_dot_df)) then
          write (stdout,*) 'BFGS: dx.df=',dx_dot_df,' NOT OK!'
       else
          write (stdout,*) 'BFGS: dx.df=',dx_dot_df,' OK'
       end if
    end if

    !To ensure the inv_Hessian stays positive-definite we may need to reset it rather than update it
    !if we get too much error accumulation due to non-exact line minimisation,
    !which will happen after too many steps or a bad update, so:
    if (((mdl%bfgs_iteration-reset_iter)>ndof).or.(dx_dot_df>tiny(dx_dot_df))) then   !dx.df ought to be <0
       if (pub_on_root.and.pub_output_detail>NORMAL) then
          write (stdout,*) 'BFGS: resetting inv_Hessian to prevent error accumulation'
          if ((mdl%bfgs_iteration-reset_iter)>ndof) write (stdout,*) '      due to number of steps taken'
          if (dx_dot_df>tiny(dx_dot_df))            write (stdout,*) '      due to bad update'
       end if
       call geom_update_inv_Hessian_params(init_cell,mdl)
       call castep_cell_copy(mdl%cell,init_cell)
       call geom_inv_Hessian_initialize(init_cell,mdl)
       reset_iter=mdl%bfgs_iteration
       return
    end if

    !All is OK so we can continue to do the update:

    h_df_vec=matmul(mdl%bfgs_inv_Hessian,df_vec)
    df_dot_h_df=dot_product(df_vec,h_df_vec)

    !this ought to be impossible with the above positive-definite test, but for safety:
    if (abs(dx_dot_df)<tiny(dx_dot_df)) then
       if (pub_on_root) write(stdout,*) 'Error in ONETEP: in geom_update_inv_Hessian - dx_dot_df->0'
       call comms_abort
    endif
    if (abs(df_dot_h_df)<tiny(df_dot_h_df)) then
       if (pub_on_root) write(stdout,*) 'Error in ONETEP: in geom_update_inv_Hessian - df_dot_h_df->0'
       call comms_abort
    endif

    u_vec=(dx_vec/dx_dot_df) - (h_df_vec/df_dot_h_df)

    do i=1,ndim
       do j=1,ndim
          mdl%bfgs_inv_Hessian(j,i)=mdl%bfgs_inv_Hessian(j,i) - dx_vec(i)*dx_vec(j)/dx_dot_df &
               &                   - h_df_vec(i)*h_df_vec(j)/df_dot_h_df &
               &                   + df_dot_h_df*u_vec(i)*u_vec(j)
       end do
    end do

    if (geom_print_inv_hessian.and.pub_on_root) then
       l=ndim/5
       do i=1,ndim
          do j=1,5*l,5
             write (stdout,99) (i,j+k,mdl%bfgs_inv_Hessian(i,j+k:j+k),k=0,4)
          end do
          if (5*l<ndim) write (stdout,99) (i,j,mdl%bfgs_inv_Hessian(i,j:j),j=5*l,ndim)
       end do
    end if

99  format('inv_Hessian:',5('(',i2,',',i2,')=',f15.10,:))

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving geom_update_inv_Hessian'
#endif

    return
  end subroutine geom_update_inv_Hessian

  subroutine geom_update_inv_Hessian_params(init_cell,mdl,final)
    !=========================================================================!
    ! Update the parameters used to initialize the inverse Hessian matrix,    !
    ! i.e. geom_modulus_est and geom_frequency_est, by analysing the current  !
    ! inv_Hessian.                                                            !
    ! This follows Pfrommer section 3 closely and uses same notation.         !
    !                                                                         !
    ! NB Cannot call this routine WITHOUT then reseting inv_Hessian as it will!
    !    update geom_modulus_est and/or geom_frequency_est which on the next  !
    !    call will then be out of sync with what has actually been used!      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   init_cell,   intent=in, cell when inv_Hessian last initialized        !
    !   mdl,         intent=in, current state of the system                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ndim and on_root                                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   geom_modulus_est & geom_frequency_est in parameters are recalculated  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    ! v2.0 now includes Householder transformation s.t. can handle mixed case !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    ! Version with allocatable arrays rather than automatics adapted from     !
    ! CASTEP 5.0.1 by Nicholas D.M. Hine, Jan 2010.                           !
    !=========================================================================!

    use comms, only: pub_on_root, pub_root_node_id, comms_abort, comms_bcast
    use constants, only: dp, stdout, NORMAL
    use rundat, only: pub_output_detail
    use simulation_cell, only: castep_cell, castep_model
    use utils, only: utils_alloc_check, utils_dealloc_check
    use wrappers, only: wrappers_invert_sym_matrix

    implicit none

    type(castep_cell),            intent(in)    :: init_cell
    type(castep_model),           intent(inout) :: mdl
    logical, optional,            intent(in)    :: final

    !local variables
    real(kind=DP), dimension(:,:), allocatable :: D_matrix    !Pfrommer eqn 11
    real(kind=DP), dimension(:,:), allocatable :: Hp_matrix   !Pfrommer eqn 12
    real(kind=DP), dimension(:,:), allocatable :: Hpp_matrix  !symmetrised Hp
    real(kind=DP), dimension(:,:), allocatable :: Hpp_update  !Hpp-H0pp
    real(kind=DP), dimension(:,:), allocatable :: H0pp_matrix !symmetrised H0
    real(kind=DP), dimension(:,:), allocatable :: tmp_Hess_matrix
    real(kind=DP), dimension(:),   allocatable :: sigma       !SVD singular values
    real(kind=DP), dimension(:,:), allocatable :: Y_matrix    !useful left SVD vectors
    real(kind=DP), dimension(:,:), allocatable :: right_vecs  !unwanted SVD vectors
    real(kind=DP), dimension(:),   allocatable :: u1_vec, e1_vec !Householder vectors
    real(kind=DP), dimension(:,:), allocatable :: Q_matrix    !Householder transformation matrix
    real(kind=DP), dimension(:,:), allocatable :: s_vec       !pure ionic  SVD vectors
    real(kind=DP), dimension(:,:), allocatable :: mu_matrix   !Pfrommer eqn 15
    real(kind=DP), dimension(1:3,1:3)          :: tmp_strain
    real(kind=DP), dimension(1:6,1:6)          :: e_vec       !pure strain SVD vectors

    real(kind=DP), dimension(:,:), allocatable :: Abar        !Hessian in reduced coords, Pfrommer eqn 16
    real(kind=DP), dimension(:,:), allocatable :: Abar_e      !strain part of Hessian in reduced coords
    real(kind=DP), dimension(:,:), allocatable :: Abar_s      !ionic  part of Hessian in reduced coords
    real(kind=DP), dimension(:,:), allocatable :: Abar_es     !mixed  part of Hessian in reduced coords

    real (kind=DP), dimension(:,:), allocatable :: S_matrix    !overlap matrix of ionic displacement vectors, Pfrommer eqn 20
    real (kind=DP), dimension(:,:,:),allocatable:: hs_vec      !basis vectors for ionic displacements
    real (kind=DP), dimension(:),   allocatable :: omega2_vec  !phonon frequencies^2 (eigenvalues of Pfrommer eqn 20)
    real (kind=DP), dimension(:,:), allocatable :: Bbar        !projected elastic stiffness Pfrommer eqn 21

    logical                                     :: found_pure_strain
    logical                                     :: found_pure_ionic
    integer                                     :: m_space     !number of sampled d.of.f
    integer                                     :: me_space    !number of sampled d.of.f in strain (epsilon) space
    integer                                     :: ndof_cell_max !max allowed d.of.f for cell variables
    integer                                     :: ms_space    !number of sampled d.of.f in ionic  (frac s)  space
    integer                                     :: ndof_ions_max !max allowed d.of.f for ion variables

    integer :: i,j,k,iatom,ispec
    integer :: ierr
    integer :: ndim3, ndof_max

    real(kind=DP) :: trace_ei, trace_ej !Traces of e_vec in Pfrommer eqn 23
    real(kind=DP) :: dummy,max_B
    real (kind=DP) :: alpha              !norm of vectors in Householder
    logical        :: update_OK, local_final

    !LAPACK requires
    integer :: work_size, info
    real(kind=DP), dimension(:), allocatable :: work

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering geom_update_inv_Hessian_params'
#endif

    if (.not.pub_on_root) return

    if (present(final)) then           !don't print analysis results here if final call
       local_final=final
    else
       local_final=.false.
    end if

    !initialize symmetric array size
    ndim3=ndim-3
    ! qoh: Initialize to avoid compiler warning
    update_ok = .false.

    !First get I+strain from init_cell to current_cell
    !call algor_invert(3,init_cell%real_lattice,inverse=tmp_strain) !using tmp_strain as temp ref_cell^-1
    !h'=(1+e)h0 => (h^T)'=(h0^T)(1+e)^T => (1+e)^T=(h0^T)^-1.(h^T)'
    !tmp_strain=transpose(matmul(tmp_strain,mdl%cell%real_lattice))
    tmp_strain = unit_matrix

    !now allocate arrays needed to generate Hpp_update
    allocate(D_matrix(1:ndim,1:ndim),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','D_matrix',ierr)
    allocate(Hp_matrix(1:ndim,1:ndim),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','Hp_matrix',ierr)
    allocate(tmp_Hess_matrix(1:ndim,1:ndim),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','tmp_Hess_matrix',ierr)
    allocate(Hpp_matrix(1:ndim3,1:ndim3),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','Hpp_matrix',ierr)
    allocate(H0pp_matrix(1:ndim3,1:ndim3),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','H0pp_matrix',ierr)
    allocate(Hpp_update(1:ndim3,1:ndim3),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','Hpp_update',ierr)

    !Now correct inv_Hessian for finite strain (Pfrommer eqn 12)
    D_matrix=0.0_DP
    do i=1,9,3
       D_matrix(i:i+2,i:i+2)=tmp_strain     !Pfrommer eqn 11
    end do
    do i=10,ndim
       D_matrix(i,i)=1.0_DP
    end do

    if (pub_on_root.and.pub_output_detail>NORMAL) write(stdout,*) 'Inverting D_matrix...'
    call wrappers_invert_sym_matrix(D_matrix,ndim)

    !Hp now describes change in enthalpy around current configuration, not original configuration
    Hp_matrix=matmul(mdl%bfgs_inv_Hessian,D_matrix)         !Hp=inv_H*(D^-1)
    Hp_matrix=matmul(transpose(D_matrix),Hp_matrix)         !Hp=(D^-1)^T*inv_H*(D^-1)

    !now restrict Hp to be symmetric in strain -> Hpp
    call unsymm_to_symm(Hp_matrix,Hpp_matrix)

    !now calculate original inv_Hessian in Hp_matrix as we have finished with that now
    !(backup current inv_Hessian in tmp_Hess_matrix)
    !NB must NOT have changed new_frequency_est or new_modulus_est between the 'real' call in main program
    !to get inv_Hessian and this call!
    tmp_Hess_matrix=mdl%bfgs_inv_Hessian
    call geom_inv_Hessian_initialize(init_cell,mdl)
    Hp_matrix=mdl%bfgs_inv_Hessian
    mdl%bfgs_inv_Hessian=tmp_Hess_matrix

    !Hp now describes change in enthalpy around current configuration, not original configuration
    Hp_matrix=matmul(Hp_matrix,D_matrix)                    !Hp=inv_H*(D^-1)
    Hp_matrix=matmul(transpose(D_matrix),Hp_matrix)         !Hp=(D^-1)^T*inv_H*(D^-1)

    !and symmetrize H0pp
    call unsymm_to_symm(Hp_matrix,H0pp_matrix)

    !now calculate update step
    Hpp_update=Hpp_matrix-H0pp_matrix

    !deallocate the ones we've finished with
    deallocate(D_matrix,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','D_matrix',ierr)
    deallocate(Hp_matrix,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','Hp_matrix',ierr)
    deallocate(tmp_Hess_matrix,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','tmp_Hess_matrix',ierr)
    !now allocate SVD and Householder stuff
    allocate(Y_matrix(1:ndim3,1:ndim3),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','Y_matrix',ierr)
    allocate(Q_matrix(1:ndim3,1:ndim3),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','Q_matrix',ierr)
    allocate(sigma(1:ndim3),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','sigma',ierr)
    allocate(right_vecs(1:ndim3,1:1),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','right_vecs',ierr)
    allocate(u1_vec(1:ndim3),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','u1_vec',ierr)
    allocate(e1_vec(1:ndim3),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','e1_vec',ierr)

    !get LAPACK optimal work size for dgesvd work array
    work_size=5*ndim3
    allocate(work(1:work_size),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','work',ierr)

    !now do SVD decomposition of update step (Hpp_update) using LAPACK
    !ie. Hpp_update=Y*sigma*right_vecs^T, where
    !    Y and right_vecs are orthogonal matrices containing singular vectors of Hpp_update
    !    and sigma is diagonal being the singular values of Hpp_update,
    !    and we overwrite Hpp_update with the left_vecs found
    !NB Left_vecs are the sensible ones, corresponding to non-zero sigma, whilst right_vecs span the null-space
    info=0
    call dgesvd('O','N',ndim3,ndim3,Hpp_update,ndim3,sigma,Y_matrix,ndim3,right_vecs,ndim3,work,work_size,info)
    !JOBU='O' => left_vecs returned in Hpp_update : 'A' => all left_vecs returned in Y_matrix
    !JOBVT='N' => right_vecs not calculated : 'A' => all right_vecs returned in right_vecs
    deallocate(work,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','work',ierr)

    !check LAPACK is OK
    !call comms_bcast(pub_root_node_id,info,1)
    if (info/=0) then
       if (pub_on_root) write (stdout,*) 'BFGS: Warning - SVD failed - skipping rest of inverse Hessian analysis'
       return
    end if

    !parallel synch
    !call comms_bcast(pub_root_node_id,sigma(1),ndim3)
    !call comms_bcast(pub_root_node_id,Hpp_update(1,1),ndim3*ndim3)

    !diagnostic
    !if ((verbose.or.iprint>2).and.pub_on_root) call geom_print_vector('SVD sigma:',ndim3,sigma)

    !now test the singular values in sigma to see which col-vecs in Hpp_update are to be put into Y for usage
    !NB sigma is guaranteed to be positive definite and stored in max->min size order, but
    !   may well contain zero values that we want to discard.

       ndof_cell_max=0

    if (fix_all_ions) then
       ndof_ions_max=0
    else
       ndof_ions_max=ndim3-6                   !(-6 for cell) but may be less due to symmetry
       if (constrain_ions) then
          ndof_ions_max=max(ndof_ions_max-num_ionic_constraints,1)
       else
          ndof_ions_max=max(ndof_ions_max-3,1) !fixed com even if not specified!
       end if
    end if

    ndof_max=ndof_cell_max+ndof_ions_max

    !This is for JOBU='O' and JOBVT='N' in dgesvd call above so left_vecs are returned in Hpp_update and not Y_matrix
    i=1
    Y_matrix=0.0_DP               !zero whole array to prevent NaN in unused cols
    Y_matrix(:,1)=Hpp_update(:,1) !there will always be at least one singular value if anything has changed
    do j=2,ndim3                  !and the sigma array has been sorted into descending order
       update_OK=(sigma(j)>epsilon(sigma(1)))                                         !parallel paranoia
       !call comms_bcast(pub_root_node_id,update_OK,1)                                 !parallel paranoia
       if (update_OK) then                                                            !non-zero singular value
          if (i<ndof_max) then
             i=i+1                                                                    ! => non-degenerate  basis vectors
             Y_matrix(:,i)=Hpp_update(:,j)                                            !    with non-trivial singular values
          end if
       end if
    end do
    m_space=i
    !so Y_matrix has m_space columns of length ndim3, each being an e-vec of Hpp_update

    !diagnostic
    !if ((verbose.or.iprint>2).and.pub_on_root) call geom_print_matrix('SVD Y:',ndim3,m_space,Y_matrix)

    !Now make new linear combinations of the e-vecs s.t. are either pure strain or pure ionic
    !using Householder transformations
    !if (.not.fix_all_cell) then
    if (.false.) then
       pure_strain: do i=1,min(6,m_space)

          !choose a column to reflect
          u1_vec=Y_matrix(:,i)
          alpha=dot_product(u1_vec,u1_vec)
          if (alpha>0.0_DP) then
             alpha=sign(sqrt(alpha),-u1_vec(1))
          else
             cycle pure_strain
          end if

          !choose reflect axis
          e1_vec=0.0_DP
          e1_vec(i)=1.0_DP
          u1_vec=u1_vec-alpha*e1_vec
          alpha=dot_product(u1_vec,u1_vec)
          if (alpha>0.0_DP) then
             alpha=sqrt(alpha)
          else
             cycle pure_strain
          end if
          u1_vec=u1_vec/alpha

          !construct Householder reflection matrix
          Q_matrix=0.0_DP
          do j=1,ndim3
             Q_matrix(j,j)=1.0_DP
          end do
          Q_matrix=Q_matrix-2*spread(u1_vec,dim=2,ncopies=ndim3)*spread(u1_vec,dim=1,ncopies=ndim3)

          !and apply it cumulatively to the original matrix
          Y_matrix=matmul(Q_matrix,Y_matrix)

       end do pure_strain

       !diagnostic
       !if ((verbose.or.iprint>2).and.pub_on_root) call geom_print_matrix('SVD QY:',ndim3,m_space,Y_matrix)
    end if

    !finished Householder
    deallocate(Hpp_update,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','Hpp_update',ierr)
    deallocate(Q_matrix,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','Q_matrix',ierr)
    deallocate(sigma,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','sigma',ierr)
    deallocate(right_vecs,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','right_vecs',ierr)
    deallocate(u1_vec,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','u1_vec',ierr)
    deallocate(e1_vec,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','e1_vec',ierr)
    !now allocate remainder
    allocate(s_vec(1:ndim3,1:ndim3),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','s_vec',ierr)
    allocate(mu_matrix(1:ndim3,1:ndim3),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','mu_matrix',ierr)

    !now sort the col-vecs s.t. strains are first, followed by pure ionic
    me_space=0                   !sampled dof in strain space
    ms_space=0                   !sampled dof in ionic  space
    e_vec=0.0_DP
    s_vec=0.0_DP
    m_space_loop: do i=1,m_space
       found_pure_ionic =.false.
       found_pure_strain =.false.

       if (.not.fix_all_ions) then
          found_pure_ionic =.true.
          do j=1,6
             if (abs(Y_matrix(j,i))>0.000001_DP) found_pure_ionic=.false.
          end do
          if (found_pure_ionic) then
             ms_space=ms_space+1
             if (ms_space>ndim3) then
                if (pub_on_root) write (stdout,*) 'BFGS: inverse Hessian analysis error - too many ionic dof found - skipping rest'
                exit m_space_loop
             end if
             s_vec(1:ndim3-6,ms_space)=Y_matrix(7:ndim3,i)
          end if
       end if

       if (.not.found_pure_ionic) then !catch double counting by giving s_vecs precedence
          found_pure_strain =.true.
          do j=7,ndim3
             if (abs(Y_matrix(j,i))>0.000001_DP) found_pure_strain=.false.
          end do
          if (found_pure_strain) then
             me_space=me_space+1
             if (me_space>6) then
                if (pub_on_root) write (stdout,*) 'BFGS: inverse Hessian analysis error - too many cell dof found - skipping rest'
                exit m_space_loop
             end if
             e_vec(1:6,me_space)=Y_matrix(1:6,i) !symmetric strain tensor in CASTEP order
          end if
       end if

    end do m_space_loop

    !catch exit cases
    !call comms_bcast(pub_root_node_id,ms_space,1)
    if (ms_space>ndim3) then
      return
    end if
    !call comms_bcast(pub_root_node_id,me_space,1)
    if (me_space>6) then
      return
    end if

    !count total number of pure e-vecs
    m_space=me_space+ms_space

    !diagnostic
    if ((pub_output_detail>NORMAL).and.pub_on_root) then
       write (stdout,*) 'geom_update_inv_Hessian_param: SVD analysis'
       write (stdout,*) '                               number of pure strain e-vecs found=',me_space
       write (stdout,*) '                               number of pure ionic  e-vecs found=',ms_space
       if (m_space<=0) write (stdout,*) '                               no pure e-vecs found - skipping rest of analysis'
    end if

    !catch bad case - nothing to do!
    if (m_space<=0) then
      return
    end if

    !now sort pure vecs back into Y and throw away any non-pure ones
    Y_matrix=0.0_DP
    do i=1,me_space
       Y_matrix(1:6,i)=e_vec(1:6,i)
    end do
    do i=1,ms_space
       Y_matrix(7:ndim3,me_space+i)=s_vec(1:ndim3-6,i)
    end do

    !now create mu^1/2 matrix (Pfrommer eq 15)
    !NB mu is diagonal, so mu^1/2 is just sqrt of diagonal elements
    mu_matrix=0.0_DP
    do i=1,6
       mu_matrix(i,i)=1.0_DP
    end do
    ispec=1
    iatom=1
    do i=7,ndim3,3
       mu_matrix(i,i)    =sqrt(mdl%cell%species_mass(ispec))
       mu_matrix(i+1,i+1)=sqrt(mdl%cell%species_mass(ispec))
       mu_matrix(i+2,i+2)=sqrt(mdl%cell%species_mass(ispec))
       if (iatom<mdl%cell%num_ions_in_species(ispec)) then
          iatom=iatom+1
       else
          iatom=1
          ispec=ispec+1
       end if
    end do

    !now create the Hessian in reduced coords (Pfrommer eqn 16)
    allocate (Abar(1:m_space,1:m_space),stat=ierr)
    call utils_alloc_check('geom_update_inv_Hessian_params','Abar',ierr)


    !we use H0pp as a working array for the multiplications:
    !Y is (ndim3)^2 but we only want (ndim3*m_space) - other RH columns are zero after sorting
    !mu is (ndim3)^2 and Hpp is (ndim3)^2
    !so the result - Abar - is (m_space)^2
    !but need (ndim3)^2 array at intermediate stages => reuse H0pp
    H0pp_matrix=matmul(mu_matrix,Y_matrix)               !ndim3^2 but only want ndim3*m_space
    H0pp_matrix=matmul(Hpp_matrix,H0pp_matrix)           !ndim3^2 but only want ndim3*m_space
    H0pp_matrix=matmul(mu_matrix,H0pp_matrix)            !ndim3^2 but only want ndim3*m_space
    H0pp_matrix=matmul(transpose(Y_matrix),H0pp_matrix)  !ndim3^2 but only want m_space*m_space
    !A=Y^T.mu^1/2.Hpp.mu^1/2.Y

    Abar(1:m_space,1:m_space)=H0pp_matrix(1:m_space,1:m_space) !extract only the chunk that we want
    if (pub_on_root.and.pub_output_detail>NORMAL) write(stdout,*) 'Inverting Abar...'
    call wrappers_invert_sym_matrix(Abar,m_space)

    !diagnostic
    !if ((verbose.or.iprint>2).and.pub_on_root) call geom_print_matrix('SVD Abar:',m_space,m_space,Abar)

    !clean up
    deallocate(Y_matrix,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','Y_matrix',ierr)
    deallocate(Hpp_matrix,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','Hpp_matrix',ierr)

    if (ms_space>0) then
       !now calculate phonon modes ...
       allocate (S_matrix(1:ms_space,1:ms_space),stat=ierr)
       call utils_alloc_check('geom_update_inv_Hessian_params','S_matrix',ierr)

       allocate (hs_vec(1:3,1:(ndim-9)/3,1:ms_space),stat=ierr)
       call utils_alloc_check('geom_update_inv_Hessian_params','hs_vec',ierr)

       !now create the Hessian sub-matrices in reduced coords (Pfrommer eqn 16)
       allocate (Abar_s(1:ms_space,1:ms_space),stat=ierr)
       call utils_alloc_check('geom_update_inv_Hessian_params','Abar_s',ierr)

       Abar_s(1:ms_space,1:ms_space)=Abar(me_space+1:m_space,me_space+1:m_space)

       !get non-fractional displacement for each atom under each s-vec eigenvector
       do i=1,ms_space
          do j=1,(ndim-9)/3
             k=j*3-2
             hs_vec(:,j,i)=matmul(transpose(mdl%cell%real_lattice),s_vec(k:k+2,i))
          end do
       end do

       !now construct ionic displacement overlap matrix
       do i=1,ms_space
          do j=1,ms_space
             S_matrix(i,j)=0.0_DP
             do k=1,(ndim-9)/3
                S_matrix(i,j)=S_matrix(i,j)+dot_product(hs_vec(:,k,i),hs_vec(:,k,j))
             end do
          end do
       end do

       !now use LAPACK to solve the generalized eigenvalue problem (Pfrommer eqn 20): (Abar_s - omega*S)*x=0
       !get LAPACK optimal block size for dsygv work array

       allocate(omega2_vec(1:ms_space),stat=ierr)
       call utils_alloc_check('geom_update_inv_Hessian_params','omega2_vec',ierr)

       !Would like to use LAPACK intrinsics to get optimal block size for dgetri work array
       !but this is not portable between different systems & LAPACK versions
       !so we do this instead:
       work_size=max(1,3*ms_space-1)

       !now allocate optimal size work array
       allocate(work(1:work_size),stat=ierr)
       call utils_alloc_check('geom_update_inv_Hessian_params','work',ierr)

       !assume Abar_s and S_matrix are symmetric, then
       info=0
       call dsygv(1,'N','U',ms_space,Abar_s,ms_space,S_matrix,ms_space,omega2_vec,work,work_size,info)
       !NB contents of Abar_s are destroyed by LAPACK call

       !finally, calculate <frequency> from generalized eigenvalues - is it really worth all this hassle???
       !NB Eigenvalues come from atomic units, where omega is an energy but when IO converts to/from energy it uses
       !   E=hbar.omega and so no two_pi factors here
       !NB We only update new_frequency_est if all is OK: freq(H2 vib)=132 THz so accept (0<freq<396 THz)
       dummy = (396.0_dp*1.5177e-4_dp)**2

       update_OK=((info==0) .and. (minval(omega2_vec)>0.0_DP) .and. (maxval(omega2_vec)<dummy))
       !call comms_bcast(pub_root_node_id,update_OK,1)

       if (update_OK) then
          new_frequency_est=0.0_DP
          do i=1,ms_space
             new_frequency_est=new_frequency_est+sqrt(omega2_vec(i))
          end do
          new_frequency_est=new_frequency_est/ms_space
          if (pub_on_root.and..not.local_final) write (stdout,12) 'BFGS: updated estimated <frequency>  =', &
               new_frequency_est,trim(frequency_label)
               !& io_atomic_to_unit(new_frequency_est,frequency_unit),trim(frequency_label)
       else
          if ((pub_output_detail>NORMAL).and.pub_on_root) then
             write (stdout,*) 'geom_update_inv_Hessian_param: skipping <frequency> update because'
             if (info/=0) write (stdout,*)                   '                               dsygv failed to converge'
             if (minval(omega2_vec)>0.0_DP) write (stdout,*) '                               min omega**2<0'
             if (maxval(omega2_vec)<dummy)  write (stdout,*) '                               max omega**2 too high'
          end if
       end if

       deallocate(omega2_vec,stat=ierr)
       call utils_dealloc_check('geom_update_inv_Hessian_params','omega2_vec',ierr)

       deallocate (work,stat=ierr)
       call utils_dealloc_check('geom_update_inv_Hessian_params','work',ierr)

       deallocate (hs_vec,stat=ierr)
       call utils_dealloc_check('geom_update_inv_Hessian_params','hs_vec',ierr)

       deallocate(S_matrix,stat=ierr)
       call utils_dealloc_check('geom_update_inv_Hessian_params','S_matrix',ierr)

       deallocate (Abar_s,stat=ierr)
       call utils_dealloc_check('geom_update_inv_Hessian_params','Abar_s',ierr)

    end if


    if (me_space>0) then
       !now calculate bulk modulus

       allocate (Bbar(1:me_space,1:me_space),stat=ierr)
       call utils_alloc_check('geom_update_inv_Hessian_params','Bbar',ierr)

       !now create the Hessian sub-matrices in reduced coords (Pfrommer eqn 16)
       allocate (Abar_e(1:me_space,1:me_space),stat=ierr)
       call utils_alloc_check('geom_update_inv_Hessian_params','Abar_e',ierr)
       Abar_e(1:me_space,1:me_space)=Abar(1:me_space,1:me_space)

       allocate (Abar_s(1:ms_space,1:ms_space),stat=ierr)
       call utils_alloc_check('geom_update_inv_Hessian_params','Abar_s',ierr)
       Abar_s(1:ms_space,1:ms_space)=Abar(me_space+1:m_space,me_space+1:m_space)

       allocate (Abar_es(1:me_space,1:ms_space),stat=ierr)
       call utils_alloc_check('geom_update_inv_Hessian_params','Abar_es',ierr)
       Abar_es(1:me_space,1:ms_space)=Abar(1:me_space,me_space+1:m_space)

       if (ms_space>0) then
          !Bbar=Abar_e-(Abar_es.Abar_s^-1.Abar_es^T) = (me_space)^2

          if (pub_on_root.and.pub_output_detail>NORMAL) write(stdout,*) 'Inverting Abar_s...'
          call wrappers_invert_sym_matrix(Abar_s,ms_space)

          !use H0pp=(ndim3)^2 to work in
          H0pp_matrix=0.0_DP
          H0pp_matrix(1:ms_space,1:me_space)=matmul(Abar_s,transpose(Abar_es))
          H0pp_matrix(1:me_space,1:me_space)=matmul(Abar_es,H0pp_matrix(1:ms_space,1:me_space))
          Bbar=Abar_e-H0pp_matrix(1:me_space,1:me_space)
       else
          Bbar=Abar_e
       end if

       if (pub_on_root.and.pub_output_detail>NORMAL) write(stdout,*) 'Inverting Bbar...'
       call wrappers_invert_sym_matrix(Bbar,me_space)

       !finally, do Pfrommer eqn 23 to get a new bulk modulus - is it really worth all this hassle???
       dummy=0.0_DP
       do i=1,me_space
          trace_ei=0.0_DP
          do k=1,3
             trace_ei=trace_ei+e_vec(k,i)   !e_vec is set of symmetric strain tensors in CASTEP order
          end do
          do j=1,me_space
             trace_ej=0.0_DP
             do k=1,3
                trace_ej=trace_ej+e_vec(k,j)   !e_vec is set of symmetric strain tensors in CASTEP order
             end do
             dummy=dummy+trace_ei*Bbar(i,j)*trace_ej
          end do
       end do

       !NB we only update new_modulus_est if all is OK (0<B<3000 GPa)
       dummy=1.0_DP/(mdl%cell%volume*dummy)
       if (pub_on_root) then
          !max_B=io_unit_to_atomic(3000.0_DP,'GPa')           !only works on root node
          max_B = 3000000.0_DP ! FIXME
          update_OK=(dummy>0.0_DP).and.(dummy<max_B)
       end if
       !call comms_bcast(pub_root_node_id,update_OK,1)

       if (update_OK) then
          new_modulus_est=dummy
          if (pub_on_root.and..not.local_final) write (stdout,12) 'BFGS: updated estimated bulk modulus =', &
               new_modulus_est,trim(pressure_label)
               !& io_atomic_to_unit(new_modulus_est,pressure_unit),trim(pressure_label)
       else
          if ((pub_output_detail>NORMAL).and.pub_on_root) then
             write (stdout,*) 'geom_update_inv_Hessian_param: skipping bulk modulus update because'
             if (dummy<0) then
                write (stdout,*) '                               update B <0'
             else
                write (stdout,*) '                               update B >3000 GPa'
             end if
          end if
       end if

       deallocate (Abar_es,stat=ierr)
       call utils_dealloc_check('geom_update_inv_Hessian_params','Abar_es',ierr)

       deallocate (Abar_s,stat=ierr)
       call utils_dealloc_check('geom_update_inv_Hessian_params','Abar_s',ierr)

       deallocate (Abar_e,stat=ierr)
       call utils_dealloc_check('geom_update_inv_Hessian_params','Abar_e',ierr)

       deallocate (Bbar,stat=ierr)
       call utils_dealloc_check('geom_update_inv_Hessian_params','Bbar',ierr)

    end if


    !clean up at the end
    deallocate (Abar,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','Abar',ierr)
    deallocate(H0pp_matrix,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','H0pp_matrix',ierr)
    deallocate(s_vec,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','s_vec',ierr)
    deallocate(mu_matrix,stat=ierr)
    call utils_dealloc_check('geom_update_inv_Hessian_params','mu_matrix',ierr)

    !parallel sync
    !call comms_bcast(pub_root_node_id,new_frequency_est,1)
    !call comms_bcast(pub_root_node_id,new_modulus_est,1)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving geom_update_inv_Hessian_params'
#endif

12  format(1x,a,f10.5,1x,a)

    return

  contains

    subroutine unsymm_to_symm(in_matrix,out_matrix)

      !Convert in_matrix (9+3N)^2 to out_matrix (6+3N)^2 by transforming
      !unsymmetrised strain space to symmetrised strain space

      !unsymm strain = 9x9 -> (e_11, e_12, e_13, e_21, e_22, e_23, e_31, e_32, e_33) square
      !  symm strain = 6x6 -> (e_11, e_22, e_33, e_23, e_31, e_12)                   square
      !where symm e_ij=0.5*(un_symm e_ij + un_symm e_ji) and using CASTEP order

      implicit none

      real(kind=DP), intent(in),  dimension(1:ndim,1:ndim)     :: in_matrix
      real(kind=DP), intent(out), dimension(1:ndim-3,1:ndim-3) :: out_matrix

      !local vars
      integer :: i,j,k


      !First, do the strain-strain bits ...
      !e11
      out_matrix(1,1)=in_matrix(1,1)
      out_matrix(2,1)=in_matrix(5,1)
      out_matrix(3,1)=in_matrix(9,1)
      out_matrix(4,1)=(in_matrix(6,1)+in_matrix(8,1))/2.0_DP
      out_matrix(5,1)=(in_matrix(3,1)+in_matrix(7,1))/2.0_DP
      out_matrix(6,1)=(in_matrix(2,1)+in_matrix(4,1))/2.0_DP
      !e22
      out_matrix(1,2)=in_matrix(1,5)
      out_matrix(2,2)=in_matrix(5,5)
      out_matrix(3,2)=in_matrix(9,5)
      out_matrix(4,2)=(in_matrix(6,5)+in_matrix(8,5))/2.0_DP
      out_matrix(5,2)=(in_matrix(3,5)+in_matrix(7,5))/2.0_DP
      out_matrix(6,2)=(in_matrix(2,5)+in_matrix(4,5))/2.0_DP
      !e33
      out_matrix(1,3)=in_matrix(1,9)
      out_matrix(2,3)=in_matrix(5,9)
      out_matrix(3,3)=in_matrix(9,9)
      out_matrix(4,3)=(in_matrix(6,9)+in_matrix(8,9))/2.0_DP
      out_matrix(5,3)=(in_matrix(3,9)+in_matrix(7,9))/2.0_DP
      out_matrix(6,3)=(in_matrix(2,9)+in_matrix(4,9))/2.0_DP
      !e23
      out_matrix(1,4)=(in_matrix(1,6)+in_matrix(1,8))/2.0_DP
      out_matrix(2,4)=(in_matrix(5,6)+in_matrix(5,8))/2.0_DP
      out_matrix(3,4)=(in_matrix(9,6)+in_matrix(9,8))/2.0_DP
      out_matrix(4,4)=(in_matrix(6,6)+in_matrix(8,6)+in_matrix(6,8)+in_matrix(8,8))/4.0_DP
      out_matrix(5,4)=(in_matrix(3,6)+in_matrix(7,6)+in_matrix(3,8)+in_matrix(7,8))/4.0_DP
      out_matrix(6,4)=(in_matrix(2,6)+in_matrix(4,6)+in_matrix(2,8)+in_matrix(4,8))/4.0_DP
      !e31
      out_matrix(1,5)=(in_matrix(1,3)+in_matrix(1,7))/2.0_DP
      out_matrix(2,5)=(in_matrix(5,3)+in_matrix(5,7))/2.0_DP
      out_matrix(3,5)=(in_matrix(9,3)+in_matrix(9,7))/2.0_DP
      out_matrix(4,5)=(in_matrix(6,3)+in_matrix(8,3)+in_matrix(6,7)+in_matrix(8,7))/4.0_DP
      out_matrix(5,5)=(in_matrix(3,3)+in_matrix(7,3)+in_matrix(3,7)+in_matrix(7,7))/4.0_DP
      out_matrix(6,5)=(in_matrix(2,3)+in_matrix(4,3)+in_matrix(2,7)+in_matrix(4,7))/4.0_DP
      !e12
      out_matrix(1,6)=(in_matrix(1,2)+in_matrix(1,4))/2.0_DP
      out_matrix(2,6)=(in_matrix(5,2)+in_matrix(5,4))/2.0_DP
      out_matrix(3,6)=(in_matrix(9,2)+in_matrix(9,4))/2.0_DP
      out_matrix(4,6)=(in_matrix(6,2)+in_matrix(8,2)+in_matrix(6,4)+in_matrix(8,4))/4.0_DP
      out_matrix(5,6)=(in_matrix(3,2)+in_matrix(7,2)+in_matrix(3,4)+in_matrix(7,4))/4.0_DP
      out_matrix(6,6)=(in_matrix(2,2)+in_matrix(4,2)+in_matrix(2,4)+in_matrix(4,4))/4.0_DP


      !... then do the strain-ion bits ...
      do j=7,ndim-3
         k=j+3
         out_matrix(1,j)=in_matrix(1,k)
         out_matrix(2,j)=in_matrix(5,k)
         out_matrix(3,j)=in_matrix(9,k)
         out_matrix(4,j)=0.5_DP*(in_matrix(6,k)+in_matrix(8,k))
         out_matrix(5,j)=0.5_DP*(in_matrix(3,k)+in_matrix(7,k))
         out_matrix(6,j)=0.5_DP*(in_matrix(2,k)+in_matrix(4,k))
      end do

      !... then do the ion-strain bits ...
      do i=7,ndim-3
         k=i+3
         out_matrix(i,1)=in_matrix(k,1)
         out_matrix(i,2)=in_matrix(k,5)
         out_matrix(i,3)=in_matrix(k,9)
         out_matrix(i,4)=0.5_DP*(in_matrix(k,6)+in_matrix(k,8))
         out_matrix(i,5)=0.5_DP*(in_matrix(k,3)+in_matrix(k,7))
         out_matrix(i,6)=0.5_DP*(in_matrix(k,2)+in_matrix(k,4))
      end do

      !... and finally do the ion-ion bits
      out_matrix(7:ndim3,7:ndim3)=in_matrix(10:ndim,10:ndim)

      return
    end subroutine unsymm_to_symm

  end subroutine geom_update_inv_Hessian_params


  subroutine geom_get_forces(mdl,elements)
    !=========================================================================!
    ! Find ground state energy, forces and stresses for given model.          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mdl, intent=inout, the model to be evaluated                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   fix_all_ions, constrain_ions, symm_ions, fix_all_cell, constrain_cell,!
    !   and symm_cell are used to determine which firstd calls we need.       !
    !   pub_output_detail, stdout and on_root for diagnostics                            !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   model for type definitions and checkpointing                          !
    !   firstd for forces and stresses                                        !
    !   electronic for electron minimiser                                     !
    !   wave for adding random noise to highest occupied band                 !
    !   parameters to define the job etc                                      !
    !   cell for cell vectors and constraints                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   mdl has already been properly read and initialised                    !
    !   current_cell must be valid                                            !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: dp, stdout
    use energy_and_force, only : energy_and_force_calculate
    use forces, only: forces_apply_constraints
    use ion, only: element
    use simulation_cell, only: castep_model, pub_cell, castep_cell_2_elements

    implicit none

    type(castep_model), intent(inout) :: mdl
    type(element),      intent(inout) :: elements(1:pub_cell%nat)

    !local variables
    integer       :: ii,is

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_get_forces'
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
             do is=1,current_cell%max_ions_in_species   ! ie max_ions_species = 1
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

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_get_forces'
#endif

    return
  end subroutine geom_get_forces

  subroutine geom_converged(mdl,old_x_vec,enthalpy,shuffle,dE,Fmax,dRmax,Smax, &
       & converged,converged_dE,converged_Fmax,converged_dRmax,converged_Smax)
    !=========================================================================!
    ! Check convergence of current configuration of the system.               !
    ! Need to satisfy convergence in enthalpy for a given number of iterations!
    ! and SIMULTANEOUSLY get force, displacement and stress < tolerances.     !
    ! Complicated by possibility of user changing length of history during    !
    ! run and possibility of backtracking to old configurations.              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mdl,    intent=in, the current state of the system                    !
    !   old_x_vec,    =in, the previous fractional coords for calculating dR  !
    !   enthalpy,     =in, the current enthalpy                               !
    !                                                                         !
    !   shuffle, intent=in,integer flag to update, or not change histories    !
    !   dE,        intent=out, width of enthalpy window                       !
    !   Fmax,          "     , max force on any ion                           !
    !   dRmax,         "     , max displacement of any ion                    !
    !   Smax,          "     , max component of cell stress                   !
    !   converged,     "     , logical flag to signal model has converged     !
    !   converged_dE,  "     , logical flag to signal dE     "     "          !
    !   converged_Fmax,"     , logical flag to signal max force    "          !
    !   converged_dRmax,"    , logical flag to signal max disp     "          !
    !   converged_Smax, "    , logical flag to signal stress       "          !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_output_detail, stdout and on_root for diagnostics                 !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !    parameters                                                           !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !  current_cell must be in sync with current model                        !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    ! kaw: Gives access to the number of classical atoms.
    use comms, only: pub_on_root,comms_abort
    use constants, only: dp, stdout, NORMAL
    use rundat, only: geom_convergence_win, geom_disp_tol, &
         geom_energy_tol, geom_force_tol, pub_output_detail, geom_continuation
    use simulation_cell, only: castep_model, castep_cell_frac_to_cart, pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    type(castep_model),               intent(in)  :: mdl
    real(kind=DP), dimension(1:ndim), intent(in)  :: old_x_vec
    real(kind=DP),                    intent(in)  :: enthalpy
    integer,                          intent(in)  :: shuffle

    real(kind=DP),                    intent(out) :: dE
    real(kind=DP),                    intent(out) :: Fmax   !max Force
    real(kind=DP),                    intent(out) :: dRmax  !max displacement
    real(kind=DP),                    intent(out) :: Smax   !max Stress
    logical,                          intent(out) :: converged
    logical,                          intent(out) :: converged_dE
    logical,                          intent(out) :: converged_Fmax
    logical,                          intent(out) :: converged_dRmax
    logical,                          intent(out) :: converged_Smax

    !local variables
    integer, save                                        :: history_iter=0
    integer, parameter                                   :: max_history_length=10

    real(kind=DP), save, dimension(1:max_history_length) :: history_E    !NB E not dE
    real(kind=DP), save, dimension(1:max_history_length) :: history_Fmax
    real(kind=DP), save, dimension(1:max_history_length) :: history_dRmax
    real(kind=DP), save, dimension(1:max_history_length) :: history_Smax

    real(kind=DP), allocatable, dimension(:,:,:)         :: old_frac
    real(kind=DP), dimension(1:3)                        :: dR, dR_frac

    integer                                              :: i,j,iatom,ispec
    integer                                              :: num_shuffle,ierr

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_converged'
#endif

    ! ars: catch errors in tag 'shuffle'
    select case (shuffle)
    case (shuffle_forwards, shuffle_none, shuffle_backwards)
    case default
       if (pub_on_root) write(stdout,*) 'Error in ONETEP: unrecognised shuffle in geom_converged:',shuffle
       call comms_abort
    end select

    !cross-check history lengths
    !why ever would you want such a long convergence window? You must be joking ....
    if (geom_convergence_win>max_history_length) then
       if (pub_on_root) write (stdout,*) 'geom_converged: truncating geom_convergence_win to ',max_history_length-1
       !NB we need the -1 to allow for backtracking
       geom_convergence_win=max_history_length-1
    end if
    !... and a window of 1 will give all sorts of problems
    if (geom_convergence_win<2) then
       if (pub_on_root) write (stdout,*) 'geom_converged: increasing geom_convergence_win to 2'
       geom_convergence_win=2
    end if

    !initialize various bits
    allocate (old_frac(1:3,1:mdl%cell%max_ions_in_species,1:mdl%cell%num_species),stat=ierr)
    call utils_alloc_check('geom_converged','old_frac',ierr)

    i=10
    do ispec=1,mdl%cell%num_species
       do iatom=1,mdl%cell%num_ions_in_species(ispec)
          do j=1,3
             old_frac(j,iatom,ispec)=old_x_vec(i+j-1)
          end do
          i=i+3
       end do
    end do

    if (history_iter==0) then
       history_E    =0.0_dp
       history_Fmax =0.0_dp
       history_dRmax=0.0_dp
       history_Smax =0.0_dp
    end if

    !calculate convergence criteria (dE has to wait until done any shuffling)

    !calculate ionic max |F| and |dR| values
    Fmax =0.0_dp
    dRmax=0.0_dp
    if (.not.fix_all_ions) then
    ! kaw: Allows for classical atoms.
       do ispec=1,mdl%cell%num_species - pub_cell%nat_classical
          do iatom=1,mdl%cell%num_ions_in_species(ispec)
             Fmax=max(Fmax,dot_product(mdl%forces(:,iatom,ispec),mdl%forces(:,iatom,ispec)))
             dR_frac=geom_opt_min_image_frac(mdl%cell%ionic_positions(:,iatom,ispec),old_frac(:,iatom,ispec))
             call castep_cell_frac_to_cart(current_cell,dR_frac,dR)
             dRmax=max(dRmax,dot_product(dR,dR))
          end do
       end do
       Fmax =sqrt(Fmax)
       dRmax=sqrt(dRmax)
    end if

    !calculate cell max |S| values
    Smax=0.0_dp

    !update history list? forwards/backwards/none?
    select case (shuffle)
    case (shuffle_forwards)       !first attempt at new configuation this iter

       !update history_iter
       history_iter=history_iter+1

       !now shuffle histories forwards, i.e. normal order in time
       num_shuffle=max(min(history_iter,max_history_length)-1,0)
       do i=num_shuffle,1,-1
          history_E(i+1)    =history_E(i)
          history_Fmax(i+1) =history_Fmax(i)
          history_dRmax(i+1)=history_dRmax(i)
          history_Smax(i+1) =history_Smax(i)
       end do

       !assign new values
       history_E(1)    =enthalpy/mdl%cell%num_ions  !enthalpy    per atom
       history_Fmax(1) =Fmax                        !max |F|  on any atom
       history_dRmax(1)=dRmax                       !max |dR| of any atom
       history_Smax(1) =Smax                        !max S component

    case (shuffle_none)           !subsequent attempts at new configuration this iter
       !assign new values, overwriting previous attempt at this iter, without updating history_iter
       history_E(1)    =enthalpy/mdl%cell%num_ions  !enthalpy    per atom
       history_Fmax(1) =Fmax                        !max |F|  on any atom
       history_dRmax(1)=dRmax                       !max |dR| of any atom
       history_Smax(1) =Smax                        !max S component

    case (shuffle_backwards)      !backtracking from current configuration to old one
       !now shuffle histories backwards as we are backtracking
       !and assign the old values to the working versions!
       !backtrack history_iter
       history_iter=history_iter-1

       num_shuffle=min(history_iter,max_history_length)-1
       do i=1,num_shuffle
          history_E(i)    =history_E(i+1)
          history_Fmax(i) =history_Fmax(i+1)
          history_dRmax(i)=history_dRmax(i+1)
          history_Smax(i) =history_Smax(i+1)
       end do

       !assign old values
       Fmax =history_Fmax(1)                        !max |F|  on any atom
       dRmax=history_dRmax(1)                       !max |dR| of any atom
       Smax =history_Smax(1)                        !max S component

    case default
       if (pub_on_root) &
            write(stdout,*) 'Error in ONETEP: unrecognised shuffle in geom_converged: ',shuffle
       call comms_abort
    end select

    !Calculate energy convergence
    !NB Energy is |max-min| enthalpy/atom over window, unlike others, as we don't know the min value!
    if (history_iter>=geom_convergence_win) then
       dE=maxval(history_E(1:geom_convergence_win))-minval(history_E(1:geom_convergence_win))
    else
       dE=maxval(history_E(1:history_iter))-minval(history_E(1:history_iter))
    end if

    !test convergence
    if (history_iter>=geom_convergence_win) then
       converged_dE=(dE<=geom_energy_tol)
    else
       converged_dE=.false.
    end if
    converged_Fmax =( Fmax<=geom_force_tol)
    if (history_iter>1) then
       converged_dRmax=(dRmax<=geom_disp_tol)
    else
       converged_dRmax=.false.
    end if

    converged_Smax =( Smax<=geom_stress_tol)

    !spoof flags where irrevelant
    if (fix_all_ions) then
       converged_Fmax =.true.
       converged_dRmax=.true.
    end if
    converged_Smax =.true.

    !special case: initial configuration satisfies F and S so don't want to do
    !any more just to satisfy dE or dR
    ! qoh: However if we are doing a continuation don't stop just because the
    ! qoh: forces are converged.
    if (mdl%bfgs_iteration==0 .and. converged_Fmax .and. converged_Smax &
         .and. .not. geom_continuation) then
       converged_dE=.true.
       converged_dRmax=.true.
    end if

    !combine all the different parts into an overall convergence flag
    converged=(converged_dE.and.converged_Fmax.and.converged_dRmax.and.converged_Smax)

    !inform user
    if (pub_output_detail>NORMAL.and.pub_on_root) then
       write (stdout,99)    'geom_converged: history_iter=',history_iter
       do i=1,min(history_iter,max_history_length)
          write (stdout,98) 'geom_converged: history_E    (',i,')=',history_E(i),trim(energy_label)
          write (stdout,98) 'geom_converged: history_Fmax (',i,')=',history_Fmax(i),trim(force_label)
          write (stdout,98) 'geom_converged: history_dRmax(',i,')=',history_dRmax(i),trim(length_label)
          write (stdout,98) 'geom_converged: history_Smax (',i,')=',history_Smax(i),trim(pressure_label)
       end do
    end if

    !clean up at end
    deallocate (old_frac,stat=ierr)
    call utils_dealloc_check('geom_converged','old_frac',ierr)

98  format (1x,a,i4,a,1x,f20.6,1x,a)
99  format (1x,a,i4)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_converged'
#endif

    return
  end subroutine geom_converged

  subroutine geom_output_converged(bfgs_iteration,enthalpy,dE,Fmax,dRmax, &
        converged_dE,converged_Fmax,converged_dRmax)
    !=========================================================================!
    ! Output relevant convergence data to stdout                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
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
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use comms, only: pub_on_root
    use constants, only: dp, stdout
    use rundat, only:  geom_disp_tol, geom_energy_tol, geom_force_tol
    use services, only: services_flush

    implicit none

    integer      , intent(in) :: bfgs_iteration
    real(kind=DP), intent(in) :: enthalpy

    real(kind=DP), intent(in) :: dE
    real(kind=DP), intent(in) :: Fmax   !max Force
    real(kind=DP), intent(in) :: dRmax  !max displacement

    logical,        intent(in) :: converged_dE
    logical,        intent(in) :: converged_Fmax
    logical,        intent(in) :: converged_dRmax

    !local vars
    character(len=80) :: data_string
    character(len=80) :: divider_string
    character(len=80) :: label_string
    integer           :: string_index
    integer           :: len_label, len_value, len_tol, len_unit, len_flag

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)')'DEBUG: Entering geom_output_converged'
#endif

    !initialise strings
    data_string     = repeat(' ',len(data_string))
    divider_string  = repeat(' ',len(divider_string))
    label_string    = repeat(' ',len(label_string))
    !                  1234567890123 23456789012345678 23456789012345678 234567890123 23456
    divider_string  = '+-----------+-----------------+-----------------+------------+-----+'
    label_string    = '| Parameter |      value      |    tolerance    |    units   | OK? |'

    !write out all the relevant data (root node only)
    if (pub_on_root) then

       !explain what we are testing ...
       write (stdout,1) 'BFGS: finished iteration',bfgs_iteration,' with enthalpy=', &
            & enthalpy,trim(energy_label)

       !header
       write(stdout,*) ' '
       write(stdout,7) divider_string
       write(stdout,7) label_string
       write(stdout,7) divider_string

       !dE                                                       '1234567890'
       string_index=1
       len_label   =13
       len_value   =18
       len_tol     =18
       len_unit    =13
       len_flag    = 6
       write(data_string(string_index:string_index+len_label),2) '  dE/ion  '
       string_index=string_index+len_label
       write(data_string(string_index:string_index+len_value),3) dE
       string_index=string_index+len_value
       write(data_string(string_index:string_index+len_tol),  4) geom_energy_tol
       string_index=string_index+len_tol
       write(data_string(string_index:string_index+len_unit), 5) trim(energy_label)
       string_index=string_index+len_unit
       if (converged_dE) then
          write(data_string(string_index:string_index+len_flag), 6) 'Yes'
       else
          write(data_string(string_index:string_index+len_flag), 6) 'No '
       end if
       string_index=string_index+len_flag
       !cross check all is well (NB strings start from 1 not 0 ...)
       if (string_index>71) write (stdout,*) 'geom_output_converged: format problem?'
       write(stdout,7) data_string

       !|F|max                                                      '1234567890'
       if (.not.fix_all_ions) then
          string_index=1
          write(data_string(string_index:string_index+len_label),2) '  |F|max  '
          string_index=string_index+len_label
          write(data_string(string_index:string_index+len_value),3) Fmax
          string_index=string_index+len_value
          write(data_string(string_index:string_index+len_tol),  4) geom_force_tol
          string_index=string_index+len_tol
          write(data_string(string_index:string_index+len_unit), 5) trim(force_label)
          string_index=string_index+len_unit
          if (converged_Fmax) then
             write(data_string(string_index:string_index+len_flag), 6) 'Yes'
          else
             write(data_string(string_index:string_index+len_flag), 6) 'No '
          end if
          write(stdout,7) data_string
       end if

       !|dR|max                                                     '1234567890'
       if (.not.fix_all_ions) then
          string_index=1
          write(data_string(string_index:string_index+len_label),2) '  |dR|max '
          string_index=string_index+len_label
          write(data_string(string_index:string_index+len_value),3) dRmax
          string_index=string_index+len_value
          write(data_string(string_index:string_index+len_tol),  4) geom_disp_tol
          string_index=string_index+len_tol
          write(data_string(string_index:string_index+len_unit), 5) trim(length_label)
          string_index=string_index+len_unit
          if (converged_dRmax) then
             write(data_string(string_index:string_index+len_flag), 6) 'Yes'
          else
             write(data_string(string_index:string_index+len_flag), 6) 'No '
          end if
          write(stdout,7) data_string
       end if

       !footer
       write(stdout,7) divider_string
       write(stdout,*) ' '

       !back to all nodes
    end if

    call services_flush

1   format(1x,a,i4,1x,a,es21.12e3,1x,a)

2   format('|',a10,    1x,'|')       !1+10+1+1=13
3   format(1x,es15.6e3,1x,'|')       !1+15+1+1=18
4   format(1x,es15.6e3,1x,'|')       !1+15+1+1=18
    !NB In IO the phys_unit is defined as 30 characters long but max. is 10 in practice
5   format(1x,a10,     1x,'|')       !1+10+1+1=13
6   format(1x,a3,      1x,'|')       !1+3+1+1=6
    !13+18+18+13+6=68
7   format(1x,a68,' <-- BFGS')

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)')'DEBUG: Leaving geom_output_converged'
#endif

    return
  end subroutine geom_output_converged

     subroutine geom_output(mdl)
     !=========================================================================!
     ! Output relevant parts of model to stdout. Based on cell_output.         !
     ! NB unlike md_output this is does NOT output energies etc as called at   !
     !    intermediate stages of calculation as well as at end.                !
     !-------------------------------------------------------------------------!
     ! Arguments:                                                              !
     !   mdl, intent=in, the model to be written                               !
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
     ! Written by Matt Probert, v1.0, 11/05/2002                               !
     !-------------------------------------------------------------------------!
     ! Modified for ONETEP by Arash A Mostofi, 2004                            !
     !=========================================================================!

     use comms, only: pub_on_root
     use constants, only: dp, stdout
     use simulation_cell, only: castep_model, castep_cell_frac_to_cart

     implicit none

     type(castep_model), intent(in) :: mdl

     !local variables
     integer                        :: i,j,k
     real(kind=DP)                  :: cart(3)

#ifdef DEBUG
      if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_output'
#endif

     ! Do it all on root node
     if(pub_on_root)then


       if (.not.fix_all_ions) then

          write(stdout,*)
          write(stdout,*)'                          -------------------&
     &------------'
          write(stdout,*)'                                    Cell Cont&
     &ents'
          write(stdout,*)'                          -------------------&
     &------------'
          write(stdout,*)

          ! write the ionic positions and species in cell
          write(stdout,7)
          write(stdout,*)'           x  Element    Atom      Absolute  &
     &    co-ordinates of atoms     x'
          write(stdout,*)'           x            Number        x      &
     &     y           z       x'
          write(stdout,8)


          do i=1,mdl%cell%num_species
             do j=1,mdl%cell%num_ions_in_species(i)
!             do j=1,1
                !NB Model is no longer rationalised and we don't bother to make it so here!
                !                write(stdout,3) mdl%cell%species_symbol(i), j, (mdl%cell%ionic_positions(k,j,i),k=1,3)
                !               write(stdout,3) mdl%cell%species_symbol(i), i, (mdl%cell%ionic_positions(k,j,i),k=1,3)
                ! aam: convert to absolute cartesian co-ordinates before outputting
                call castep_cell_frac_to_cart(current_cell, &
     &               mdl%cell%ionic_positions(:,j,i),cart)
                write(stdout,3) mdl%cell%species_symbol(i), i, &
     &(cart(k),k=1,3)
             end do
          end do
          write(stdout,7)
3         format(11x,' x',a6,1x,i8,4x,3f12.6,'   x')
7         format(1x,'           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx&
     &xxxxxxxxxxxxxxxxxxxxx')
8         format(1x,'           x--------------------------------------&
     &--------------------x')

       end if

       write(stdout,*)

       ! Back to all nodes
     endif

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_output'
#endif

    return
  end subroutine geom_output

  subroutine geom_write_trajectory(mdl,enthalpy,output_file)
    !=========================================================================!
    ! Output relevant geom_data in AU to specified trajectory file.           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mdl, intent=in, the model to be written                               !
    !   enthalpy, intent=in, the current enthalpy                             !
    !   output_file, intent=in, the filename to which results will be written.!
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
    ! Written by Matt Probert, v1.0, 11/05/2002                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    ! kaw: Gives access to the number of classical atoms.
    use comms, only: pub_on_root
    use constants, only: dp, stdout
    use simulation_cell, only: castep_model, castep_cell_frac_to_cart, &
                               pub_cell
    use utils, only : utils_unit,utils_open_unit_check

    implicit none

    type(castep_model), intent(in)  :: mdl
    real(kind=DP),    intent(in)    :: enthalpy
    character(len=*),  intent(in)   :: output_file

    !local variables
    integer                         :: i,iatom,ispec
    real(kind=DP), dimension(1:3)   :: frac_pos, abs_pos
    character(len=12)               :: action,form,stat,position,access

    integer                         :: kk

    !for filehandling we need
    integer                         :: out_unit, io_status

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering geom_write_trajectory'
#endif

    action   = 'WRITE'
    form     = 'FORMATTED'
    stat     = 'UNKNOWN'
    position = 'REWIND'
    access   = 'SEQUENTIAL'

    if (mdl%bfgs_iteration==0) then
       !create a new output file
       stat     ='REPLACE'
    else
       !append to existing output file
       position ='APPEND'
    end if

    if (pub_on_root) then

       ! Find next available unit specifier
       out_unit = utils_unit()

       open(unit=out_unit,iostat=io_status,file=output_file,status=stat,&
            access=access,form=form,position=position,action=action)
       call utils_open_unit_check('geom_write_trajectory','output_file',io_status)



       !now write out all the relevant data in atomic units (root node only)
       write (out_unit,1) mdl%bfgs_iteration
       write (out_unit,2) mdl%total_energy,enthalpy

       !cell vectors
       do i=1,3
          write (out_unit,4) mdl%cell%real_lattice(i,:)
       end do

       !ion stuff
       if (.not.fix_all_ions) then

          !ionic positions
          do ispec=1,mdl%cell%num_species
             do iatom=1,mdl%cell%num_ions_in_species(ispec)
                frac_pos(1:3)=mdl%cell%ionic_positions(1:3,iatom,ispec)
                ! rationalise fractional co-ordinates before
                ! converting to cartesians and writing to file
                do kk=1,3
                   if (frac_pos(kk).lt.0.0_dp) then
                      frac_pos(kk)=frac_pos(kk)+1.0_dp
                   else if (frac_pos(kk).gt.1.0_dp) then
                      frac_pos(kk)=frac_pos(kk)-1.0_dp
                   else
                      continue
                   endif
                enddo
                call castep_cell_frac_to_cart(current_cell,frac_pos,abs_pos)
                write (out_unit,7) mdl%cell%species_symbol(ispec),iatom, abs_pos(:)
             end do
          end do

          !ionic forces
    ! kaw: modified this to allow for classical atoms.
          do ispec=1,mdl%cell%num_species - pub_cell%nat_classical
             do iatom=1,mdl%cell%num_ions_in_species(ispec)
                write (out_unit,9)  mdl%cell%species_symbol(ispec),iatom, mdl%forces(:,iatom,ispec)
             end do
          end do

       end if

       !blank line to signal end of this iter (datablock size variable)
       write (out_unit,*) ' '

       close (unit=out_unit)

    end if

    !NB Try to keep format numbers/layout same here as in md and fit into 80-char page-width
1   format(12x,i18)
2   format(9x,2(3x,es18.8e3),21x,'  <-- E')
4   format(9x,3(3x,es18.8e3),    '  <-- h')
7   format(1x,a3,1x,i4,3(3x,es18.8e3),'  <-- R')
9   format(1x,a3,1x,i4,3(3x,es18.8e3),'  <-- F')

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving geom_write_trajectory'
#endif

    return
  end subroutine geom_write_trajectory


  function geom_opt_min_image_frac(f1,f2)
    !=========================================================================!
    ! Handles minimum image convention for periodic boundary conditions by    !
    ! applying translation vectors t of Bravais lattice to minimize r12, where!
    ! f12 = |f1 - f2 - t| for input f1 & f2 vectors, for any Bravais lattice. !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   f1,   intent=in, fractional input vector1                             !
    !   f2,   intent=in, fractional input vector2                             !
    ! Returns:                                                                !
    !   fractional vector f12=f1-f2                                           !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use constants, only: dp
    use simulation_cell, only:castep_cell_frac_to_cart, castep_cell_cart_to_frac

    implicit none

    real(kind=DP),dimension(3), intent(in) :: f1,f2    ! Fractional coords
    real(kind=DP), dimension(3) :: geom_opt_min_image_frac

    ! <<< local variables >>>
    real(kind=DP), dimension(3) :: r1,r2,r12,f12

    ! Get absolute Cartesian co-ordinates
    call castep_cell_frac_to_cart(current_cell,f1,r1)
    call castep_cell_frac_to_cart(current_cell,f2,r2)

    ! Get minimum image
    r12 = geom_opt_min_image_cart(r1,r2)

    ! Set to fractional
    call castep_cell_cart_to_frac(current_cell,r12,f12)

    ! Set function
    geom_opt_min_image_frac=f12

    return
  end function geom_opt_min_image_frac


  function geom_opt_min_image_cart(r1,r2)
    !=========================================================================!
    ! Handles minimum image convention for periodic boundary conditions by    !
    ! applying translation vectors t of Bravais lattice to minimize r12, where!
    ! r12 = |r1 - r2 - t| for input r1 & r2 vectors, for any Bravais lattice. !
    ! Routine taken from SCAMPI (mijp's path integral MD code)                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   r1,   intent=in, Cartesian input vector1                              !
    !   r2,   intent=in, Cartesian input vector2                              !
    ! Returns:                                                                !
    !   Cartesian vector r12=r1-r2                                            !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! current_cell                                                            !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    ! parameters, io                                                          !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   transvec = the 8 different translation vectors = 't' in above         !
    !            = 0, a, b, c, a+b, a+c, b+c, a+b+c                           !
    !   r12      = the basic r1-r2 vector in the unit cell                    !
    !   mod_r12  = the 8 different values of |r12-t| in the unit cell         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.2, 05/01/2001                               !
    !-------------------------------------------------------------------------!
    ! Modified for ONETEP by Arash A Mostofi, 2004                            !
    !=========================================================================!

    use constants, only: dp
    use simulation_cell, only: castep_cell_cart_to_frac

    implicit none

    real(kind=DP), dimension(1:3), intent(in) :: r1, r2 !Cartesian
    real(kind=DP), dimension(1:3)             :: geom_opt_min_image_cart

    !local variables
    real(kind=DP), dimension(1:3)             :: r12
    real(kind=DP), dimension(1:3)             :: fract_r12
    real(kind=DP), dimension(1:8,1:3)         :: transvec
    real(kind=DP), dimension(1:8)             :: mod_r12
    real(kind=DP), dimension(1:3)             :: min_r12
    real(kind=DP)                             :: min_mod
    integer                                   :: i,j

    !Calculate basic r12=r1-r2 displacement vector ...
    r12=r1-r2

    call castep_cell_cart_to_frac(current_cell,r12,fract_r12)

    !Find origin of cell containing r12 ...
    do i=1,3
       if (fract_r12(i)>=0.0_dp) then
          fract_r12(i)=aint(fract_r12(i),kind=dp)
       else
          fract_r12(i)=aint(fract_r12(i)-1.0_dp,kind=dp)
       end if
    end do

    !... and map r12 s.t. new r12(i)= new fract_r12(i)*real_lattice(i),
    !where 0 <= new fract_r12 < 1 ...
    do i=1,3
       do j=1,3
          r12(i)=r12(i)-fract_r12(j)*current_cell%real_lattice(j,i)
       end do
    end do

    !... and now calculate translation vectors from each of eight corners
    !(of parallelepiped that contains r12) to r12 ...

    transvec(1,:)=r12(:)
    transvec(2,:)=r12(:)-current_cell%real_lattice(1,:)
    transvec(3,:)=r12(:)-current_cell%real_lattice(2,:)
    transvec(4,:)=r12(:)-current_cell%real_lattice(3,:)
    transvec(5,:)=r12(:)-current_cell%real_lattice(1,:)-current_cell%real_lattice(2,:)
    transvec(6,:)=r12(:)-current_cell%real_lattice(1,:)-current_cell%real_lattice(3,:)
    transvec(7,:)=r12(:)-current_cell%real_lattice(2,:)-current_cell%real_lattice(3,:)
    transvec(8,:)=r12(:)-current_cell%real_lattice(1,:)-current_cell%real_lattice(2,:)-current_cell%real_lattice(3,:)

    mod_r12=0.0_dp
    do i=1,3
       mod_r12(:)=mod_r12(:)+transvec(:,i)**2
    end do

    !... and finally chose the one with minimum length ...
    min_mod=mod_r12(1)
    min_r12(:)=transvec(1,:)
    do i=2,8
       if (mod_r12(i)<min_mod) then
          min_mod=mod_r12(i)
          min_r12(:)=transvec(i,:)
       end if
    end do

    ! Set function
    geom_opt_min_image_cart=min_r12

    return
  end function geom_opt_min_image_cart




  subroutine geom_opt_continuation_read(mdl,elements,filename)
    !=========================================================================!
    ! Read all relevant data from file for geometry optimisation continuation !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mdl, intent=inout, the model                                          !
    !   elements, intent=in, the element data array                           !
    !   filename, intent=in, the filename from which continuation             !
    !                           data will be read.                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !    model type                                                           !
    !    ndim                                                                 !
    !    pub_cell                                                             !
    !    current_cell                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !    rundat for geom_modulus_est, geom_frequency_est                      !
    !    comms for pub_on_root                                                !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   mdl has already been properly read and initialised                    !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, 08/02/2005                                  !
    !=========================================================================!

    use comms, only: pub_on_root, comms_abort, comms_barrier, comms_bcast, pub_root_node_id
    use constants, only: DP, stderr, stdout
    use ion, only: element
    use rundat, only: geom_modulus_est,geom_frequency_est, pub_rootname
    use simulation_cell, only: castep_model, pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
                     utils_open_unit_check, utils_close_unit_check, &
                     utils_read_check, utils_write_check, utils_binary_copy

    implicit none

    type(castep_model), intent(inout)       :: mdl
    type(element), intent(in)               :: elements(pub_cell%nat)
    character(len=file_maxpath), intent(in) :: filename

    ! <<< local variables >>>
    integer :: read_unit,write_unit
    integer :: ios,ierr
    integer :: atom,ngwf,dim
    integer :: ionidx
    integer :: num,tn1,tn2,tn3,tn(3)
#ifdef ACCELRYS
    integer :: row           ! vm: required for Intel bug workaround
#endif
    real(kind=DP) :: orig1,orig2,orig3
    real(kind=DP), allocatable, dimension(:,:,:) :: uni_tbox
    character(len=file_maxpath) :: ngwf_file,dk_file

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering geom_opt_continuation_read'
#endif

    if (pub_on_root) then

       ! Find available unit specifier for reading continuation file
       read_unit = utils_unit()

       ! Open continuation file
       open(unit=read_unit,file=filename,iostat=ios,&
            form='UNFORMATTED',action='READ',status='OLD')
       call utils_open_unit_check('geom_opt_continuation_read','filename',ios)
       rewind(read_unit)

       ! Tell user what is going on
       write(stdout,'(3a)',advance ='no') &
            ' Reading continuation file "',trim(filename),'" ...'

       ! Read bfgs_iteration
       read(read_unit,iostat=ios) mdl%bfgs_iteration
       call utils_read_check('geom_opt_continuation_read','mdl%bfgs_iteration',ios)

       ! Read total_energy
       read(read_unit,iostat=ios) mdl%total_energy
       call utils_read_check('geom_opt_continuation_read','mdl%total_energy',ios)

       ! Read number of spins
       read(read_unit,iostat=ios) pub_cell%num_spins
       call utils_read_check('geom_opt_continuation_read','pub_cell%num_spins',ios)

       ! Read ionic positions in fractional co-ordinates into current_cell
       ! NB it is assumed that each ion has its own species
       do ionidx=1,pub_cell%nat
          read(read_unit,iostat=ios) current_cell%ionic_positions(1:3,1,ionidx)
          call utils_read_check('geom_opt_continuation_read','current_cell%ionic_positions',ios)
       end do

       ! Read ionic forces
       ! NB it is assumed that each ion has its own species
       do ionidx=1,pub_cell%nat
          read(read_unit,iostat=ios) mdl%forces(1:3,1,ionidx)
          call utils_read_check('geom_opt_continuation_read','mdl%forces',ios)
       end do

       ! If constrain_ions is TRUE, then mdl%orig_forces will have been
       ! allocated (in castep_model_alloc) and we copy mdl%forces to it
       ! NB it is assumed that each ion has its own species
       if (constrain_ions) then
          do ionidx=1,pub_cell%nat
             mdl%orig_forces(1:3,1,ionidx) = mdl%forces(1:3,1,ionidx)
          end do
       end if

       ! Read inverse Hessian
       !  read(read_unit,iostat=ios) mdl%bfgs_inv_Hessian(1:ndim,1:ndim)
       do dim = 1,ndim
          read(read_unit,iostat=ios) mdl%bfgs_inv_Hessian(dim,1:ndim)
       enddo
       call utils_read_check('geom_BFGS_restore','mdl%bfgs_inv_Hessian',ios)

       ! Read geom_modulus_est
       read(read_unit,iostat=ios) geom_modulus_est
       call utils_read_check('geom_BFGS_restore','geom_modulus_est',ios)

       ! Read geom_frequency_est
       read(read_unit,iostat=ios) geom_frequency_est
       call utils_read_check('geom_BFGS_restore','geom_frequency_est',ios)

       !******************************************************************!
       ! AAM: read ngwfs from continuation file and write to              !
       !      .tightbox_ngwfs file                                        !
       !******************************************************************!

       ! Find available unit for writing ngwfs
       write_unit = utils_unit()

       ! Set up name of current NGWF file to be written
       write(ngwf_file,'(2a)') trim(pub_rootname),'.tightbox_ngwfs'

       ! Open current NGWF file to be written
       open(unit=write_unit,file=ngwf_file,form='UNFORMATTED', &
            action='WRITE',iostat=ios)
       call utils_open_unit_check('geom_opt_continuation_read',ngwf_file,ios)
       rewind(write_unit)

       ! Read number of NGWFs
       read(read_unit,iostat=ios) num
       call utils_read_check('geom_opt_continuation_read','num',ios)

       ! Sanity check
       if (num /= pub_cell%num_ngwfs) then
          write(stderr,'(2(a,i6),a)') 'Error in geom_opt_continuation_read: &
               &num (',num,') not equal to pub_cell%num_ngwfs (',pub_cell%num_ngwfs,')'
          call comms_abort
       end if

       ! Write number of NGWFs to file
       write(write_unit,iostat=ios) num
       call utils_write_check('geom_opt_continuation_read','num',ios)

       ! Loop over dimensions
       do dim=1,3

          ! Read number of points in dimension dim of universal tightbox
          read(read_unit,iostat=ios) tn(dim)
          call utils_read_check('geom_opt_continuation_read','tn',ios)

          ! Write number of points in dimension dim of universal tightbox
          write(write_unit,iostat=ios) tn(dim)
          call utils_write_check('geom_opt_continuation_read','tn',ios)

       end do

       ! Allocate universal tightbox
       tn1 = tn(1) ; tn2 = tn(2) ; tn3 = tn(3)
       allocate(uni_tbox(tn1,tn2,tn3),stat=ierr)
       call utils_alloc_check('geom_opt_continuation_read','uni_tbox',ierr)
       uni_tbox = 0.0_dp

       atom_loop: do atom=1,pub_cell%nat

          ngwfs_on_atom_loop: do ngwf=1,elements(atom)%nfunctions

             ! Read origin
             read(read_unit,iostat=ios) orig1, orig2, orig3
             call utils_read_check('geom_opt_continuation_read','orig1, orig2, orig3',ios)

             ! Write origin
             write(write_unit,iostat=ios) orig1, orig2, orig3
             call utils_write_check('geom_opt_continuation_read','orig1, orig2, orig3',ios)

             ! Read tightbox
             read(read_unit,iostat=ios) uni_tbox(1:tn1,1:tn2,1:tn3)
             call utils_read_check('geom_opt_continuation_read','uni_tbox',ios)

             ! Write tightbox
             write(write_unit,iostat=ios) uni_tbox(1:tn1,1:tn2,1:tn3)
             call utils_write_check('geom_opt_continuation_read','uni_tbox',ios)

          end do ngwfs_on_atom_loop

       end do atom_loop

       ! Close current NGWF file
       close(write_unit,iostat=ios)
       call utils_close_unit_check('geom_opt_continuation_read','write_unit',ios)

       ! Deallocate memory
       deallocate(uni_tbox,stat=ierr)
       call utils_dealloc_check('geom_opt_continuation_read','uni_tbox',ierr)

       !******************************************************************!
       !================   End read/write tightbox_ngwfs  ================!
       !******************************************************************!

       !******************************************************************!
       ! AAM: read current density kernel from continuation file and      !
       !      write .dkn file                                             !
       !******************************************************************!

       ! Find available unit for current density kernel file to be written
       write_unit = utils_unit()

       ! Name of current density kernel file to be written
       write(dk_file,'(2a)') trim(pub_rootname),'.dkn'

       ! Open density kernel file
       open(unit=write_unit,file=dk_file,form='UNFORMATTED', &
            action='WRITE',iostat=ios)
       call utils_open_unit_check('geom_opt_continuation_read','dk_file',ios)
       rewind(write_unit)

       ! Read in density kernel file
       call utils_binary_copy(read_unit,write_unit)

       ! Close current density kernel file
       close(write_unit,iostat=ios)
       call utils_close_unit_check('geom_opt_continuation_read','write_unit',ios)

       !******************************************************************!
       !============= END read/write current density kernel ==============!
       !******************************************************************!

       ! Close read_unit
       close(read_unit,iostat=ios)
       call utils_close_unit_check('geom_opt_continuation_read','read_unit',ios)

       write(stdout,'(a)') 'done'

    end if

    call comms_barrier

    ! Broadcast data to all nodes
    call comms_bcast(pub_root_node_id,mdl%bfgs_iteration)
    call comms_bcast(pub_root_node_id,mdl%total_energy)
    call comms_bcast(pub_root_node_id,mdl%forces)
    !call comms_bcast(pub_root_node_id,mdl%bfgs_inv_Hessian)
    call comms_bcast(pub_root_node_id,current_cell%ionic_positions)
    call comms_bcast(pub_root_node_id,geom_modulus_est)
    call comms_bcast(pub_root_node_id,geom_frequency_est)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving geom_opt_continuation_read'
#endif

    return
  end subroutine geom_opt_continuation_read

  subroutine geom_opt_continuation_write(mdl,elements,output_file)
    !=========================================================================!
    ! Write all relevant data to file for geometry optimisation continuation  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   mdl, intent=in, the model                                             !
    !   elements, intent=in, the element data array                           !
    !   output_file, intent=in, the filename to which continuation            !
    !                           data will be written.                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !    model type                                                           !
    !    ndim                                                                 !
    !    pub_cell                                                             !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !    comms for pub_on_root                                                !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   mdl has already been properly read and initialised                    !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, 08/02/2005                                  !
    ! Coordinate output in xyz format by Chris-Kriton Skylaris, 28/03/2008    !
    !=========================================================================!

    use comms, only: pub_on_root, comms_abort
    use constants, only: DP, stdout, stderr
    use ion, only: element
    use rundat, only: pub_write_xyz, pub_rootname
    use services, only: services_write_xyz
    use simulation_cell, only: castep_model, pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
                     utils_open_unit_check, utils_close_unit_check, &
                     utils_read_check, utils_write_check, utils_binary_copy

    implicit none

    type(castep_model), intent(in)          :: mdl
    type(element), intent(in)               :: elements(pub_cell%nat)
    character(len=file_maxpath), intent(in) :: output_file

    ! <<< local variables >>>
    integer :: output_unit,read_unit,ios,ierr,ionidx
    integer :: atom,ngwf,dim
    integer :: num,tn1,tn2,tn3,tn(3)
#ifdef ACCELRYS
    integer :: row           ! vm: required for Intel bug workaround
#endif
    real(kind=DP) :: orig1,orig2,orig3
    real(kind=DP), allocatable, dimension(:,:,:) :: uni_tbox
    character(len=file_maxpath) :: read_file
    character(len=file_maxpath) :: xyz_title_line ! title line of coords in xyz file

#ifdef DEBUG
    if (pub_on_root) &
         write(stdout,'(a)') 'DEBUG: Entering geom_opt_continuation_write'
#endif

    if (pub_on_root) then

       ! cks: -------- XYZ OUTPUT ----------------------------------------
       if (pub_write_xyz) then
          write(xyz_title_line,'(a, i5)')"Geometry optimisation iteration: ", mdl%bfgs_iteration
          call services_write_xyz(elements, pub_rootname, trim(xyz_title_line) )
       endif
       ! cks: ---- END XYZ OUTPUT ----------------------------------------



       ! Find available unit specifier for writing continuation file
       output_unit = utils_unit()

       ! Open output_file
       open(unit=output_unit,file=output_file,iostat=ios,&
            form='UNFORMATTED',action='WRITE')
       call utils_open_unit_check('geom_opt_continuation_write',output_file,ios)
       rewind(output_unit)

       ! Tell user what is going on
       write(stdout,'(/3a)',advance='no') ' Writing continuation file "', &
            trim(output_file),'" ...'

       ! Write bfgs_iteration
       write(output_unit,iostat=ios) mdl%bfgs_iteration
       call utils_write_check('geom_opt_continuation_write','mdl%bfgs_iteration',ios)

       ! Write total_energy
       write(output_unit,iostat=ios) mdl%total_energy
       call utils_write_check('geom_opt_continuation_write','mdl%total_energy',ios)

       ! Write number of spins
       write(output_unit,iostat=ios) pub_cell%num_spins
       call utils_write_check('geom_opt_continuation_write','pub_cell%num_spins',ios)

       ! Write ionic positions in fractional co-ordinates
       ! NB it is assumed that each ion has its own species
       do ionidx=1,mdl%cell%num_ions
          write(output_unit,iostat=ios) mdl%cell%ionic_positions(1:3,1,ionidx)
          call utils_write_check('geom_opt_continuation_write','mdl%cell%ionic_positions',ios)
       end do

       ! Write the UNCONSTRAINED ionic forces so that on continuation different
       ! ionic constraints may be used
       ! NB it is assumed that each ion has its own species
       !do ionidx=1,mdl%cell%num_ions
       ! kaw: We only need forces for the QM ions
       do ionidx=1,pub_cell%nat
          if (constrain_ions) then
             write(output_unit,iostat=ios) mdl%orig_forces(1:3,1,ionidx)
             call utils_write_check('geom_opt_continuation_write','mdl%orig_forces',ios)
          else
             write(output_unit,iostat=ios) mdl%forces(1:3,1,ionidx)
             call utils_write_check('geom_opt_continuation_write','mdl%forces',ios)
          end if
       end do

       ! Write inverse Hessian
       ! write(output_unit,iostat=ios) mdl%bfgs_inv_Hessian(1:ndim,1:ndim)
       do dim=1,ndim
          write(output_unit,iostat=ios) mdl%bfgs_inv_Hessian(dim,1:ndim)
       end do
       call utils_write_check('geom_opt_continuation_write','mdl%bfgs_inv_Hessian',ios)

       ! Write the latest estimate of the bulk modulus
       write(output_unit,iostat=ios) new_modulus_est
       call utils_write_check('geom_opt_continuation_write','new_modulus_est',ios)

       ! Write the latest estimate of the phonon frequency
       write(output_unit,iostat=ios) new_frequency_est
       call utils_write_check('geom_opt_continuation_write','new_frequency_est',ios)

       !******************************************************************!
       ! AAM: read current tightbox_ngwfs from file and write to          !
       !      continuation file                                           !
       !******************************************************************!

       ! Find available unit for reading current NGWFs
       read_unit = utils_unit()

       ! Set up name of current NGWF file to be read
       write(read_file,'(a,a)') trim(pub_rootname),'.tightbox_ngwfs'

       ! Open current NGWF file to be read
       open(unit=read_unit,file=read_file,form='UNFORMATTED',status='OLD', &
            action='READ',iostat=ios)
       call utils_open_unit_check('geom_opt_continuation_write',read_file,ios)
       rewind(read_unit)

       ! Read number of NGWFs
       read(read_unit,iostat=ios) num
       call utils_read_check('geom_opt_continuation_write','num',ios)


       ! Sanity check
       if (num /= pub_cell%num_ngwfs) then
          write(stderr,'(2(a,i6),a)') 'Error in geom_opt_continuation_write: &
               &num (',num,') not equal to pub_cell%num_ngwfs (',pub_cell%num_ngwfs,')'
          call comms_abort
       end if

       ! Write number of NGWFs to backup file
       write(output_unit,iostat=ios) num
       call utils_write_check('geom_opt_continuation_write','num',ios)

       ! Loop over dimensions
       do dim=1,3

          ! Read number of points in dimension dim of universal tightbox
          read(read_unit,iostat=ios) tn(dim)
          call utils_read_check('geom_opt_continuation_write','tn',ios)

          ! Write number of points in first dimension of universal tightbox
          write(output_unit,iostat=ios) tn(dim)
          call utils_write_check('geom_opt_continuation_write','tn',ios)

       end do

       ! Allocate universal tightbox
       tn1 = tn(1) ; tn2 = tn(2) ; tn3 = tn(3)
       allocate(uni_tbox(tn1,tn2,tn3),stat=ierr)
       call utils_alloc_check('geom_opt_continuation_write','uni_tbox',ierr)
       uni_tbox = 0.0_dp

       atom_loop: do atom=1,pub_cell%nat

          ngwfs_on_atom_loop: do ngwf=1,elements(atom)%nfunctions

             ! Read origin
             read(read_unit,iostat=ios) orig1, orig2, orig3
             call utils_read_check('geom_opt_continuation_write','orig1, orig2, orig3',ios)

             ! Write origin
             write(output_unit,iostat=ios) orig1, orig2, orig3
             call utils_write_check('geom_opt_continuation_write','orig1, orig2, orig3',ios)

             ! Read tightbox
             read(read_unit,iostat=ios) uni_tbox(1:tn1,1:tn2,1:tn3)
             call utils_read_check('geom_opt_continuation_write','uni_tbox',ios)

             ! Write tightbox
             write(output_unit,iostat=ios) uni_tbox(1:tn1,1:tn2,1:tn3)
             call utils_write_check('geom_opt_continuation_write','uni_tbox',ios)

          end do ngwfs_on_atom_loop

       end do atom_loop

       ! Close current NGWF file
       close(read_unit,iostat=ios)
       call utils_close_unit_check('geom_opt_continuation_write','read_unit',ios)

       ! Deallocate memory
       deallocate(uni_tbox,stat=ierr)
       call utils_dealloc_check('geom_opt_continuation_write','uni_tbox',ierr)

       !******************************************************************!
       !================   End read/write tightbox_ngwfs  ================!
       !******************************************************************!

       !******************************************************************!
       ! AAM: read current density kernel from file and write to          !
       !      continuation file                                           !
       !******************************************************************!

       ! Name of current density kernel file to be read

       write(read_file,'(2a)') trim(pub_rootname),'.dkn'

       ! Find available unit for current density kernel file to be read
       read_unit = utils_unit()

       ! Open current density kernel file
       open(unit=read_unit,file=read_file,form='UNFORMATTED', &
            status='OLD',action='READ',iostat=ios)
       call utils_open_unit_check('geom_opt_continuation_write',read_file,ios)
       rewind(read_unit)

       ! Read in density kernel file
       call utils_binary_copy(read_unit,output_unit)

       ! Close current density kernel file
       close(read_unit,iostat=ios)
       call utils_close_unit_check('geom_opt_continuation_write','read_unit',ios)

       !******************************************************************!
       !============= END read/write current density kernel ==============!
       !******************************************************************!

       ! Close output_unit
       close(output_unit,iostat=ios)
       call utils_close_unit_check('geom_opt_continuation_write','output_unit',ios)

       write(stdout,'(a)') 'done'

    end if

#ifdef DEBUG
    if (pub_on_root) &
         write(stdout,'(a)') 'DEBUG: Leaving geom_opt_continuation_write'
#endif

    return
  end subroutine geom_opt_continuation_write

end module geometry_optimiser
