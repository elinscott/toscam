! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!                          E W A L D                                          !
!=============================================================================!
!                                                                             !
! $Id: ewald.F90,v 1.69 2009/10/20 14:02:02 kr264 Exp $
!-----------------------------------------------------------------------------!
! This module calculates the ion-ion contribution to the total energy,        !
! (and its derivatives) using the Ewald sum technique.                        !
!                                                                             !
!-----------------------------------------------------------------------------!
! Written from  "Module Specification for New Plane Wave Code",  M. Segall,   !
!    P. Lindan, M. Probert, C. Pickard, P. Hasnip, S. Clark and M. Payne      !
!                           Copyright 1999/2000                               !
!-----------------------------------------------------------------------------!
! Written by Matt Probert, v0.1, 06/07/2000                                   !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log: ewald.F90,v $
! Revision 1.69  2009/10/20 14:02:02  kr264
! Changes to comms_reassign_strategy to allow idle node pickup
!
! Revision 1.68.100.1  2009/10/11 10:21:48  kr264
! Updated saved-state checks to allow for old data differing on previously idle nodes.
!
! Revision 1.68  2009/09/02 14:37:17  kr264
! Replaced call to IO_ABORT on ALLOCATE failure with
! call to specific handler routine IO_ALLOCATE_ABORT
!
! Revision 1.67  2009/08/17 10:02:08  kr264
! Fixed bug in Ewald_calc_e2_dipole_dipole where the Born charge tensors were transposed.
!
! Revision 1.66  2009/08/05 11:02:39  pjh1003
! Added support for band-parallelism.
!
! Revision 1.65  2009/06/29 22:59:59  mijp2
! slackened tolerances due to memory costs
!
! Revision 1.64  2009/05/26 13:29:20  mijp2
! incorporated Keiths skewness suggestion and re-optimised parameters
!
! Revision 1.63  2009/05/22 13:51:39  mijp2
! tightened default precision in eta calculation
!
! Revision 1.62  2008/10/01 16:37:31  mijp2
! added EW3D correction via devel_code
!
! Revision 1.61  2008/09/24 14:09:42  mijp2
! logic switched in defaults for ewald_calc_storage
!
! Revision 1.60  2008/08/29 13:50:53  kr264
! Fixed character argument length bug (from A.Perlov)
!
! Revision 1.59  2008/08/24 23:00:50  mijp2
! added ewald_calc_storage
!
! Revision 1.58  2008/04/28 14:06:14  mijp2
! re-optimised precision and r2r after 1.57 bugfix
!
! Revision 1.57  2008/02/25 14:15:18  mijp2
! fixed cut-offs problem which caused molecule-in-a-box bug [08049vymq02]
!
! Revision 1.56  2008/01/18 15:14:51  mijp2
! eliminated redundant array in ewald_calculate_forces
!
! Revision 1.55  2008/01/08 15:53:49  mijp2
! tidied trace
!
! Revision 1.54  2007/05/17 16:35:26  kr264
! Fixed incorrect trace call.
!
! Revision 1.53  2007/05/17 14:05:21  pjh1003
! Merged trace branch.
!
! Revision 1.52  2007/04/19 17:20:49  mijp2
! tighted precision to fix castep_2c ARTS test
!
! Revision 1.51  2007/03/27 21:51:25  mijp2
! deleted redundant E2 routines
!
! Revision 1.50  2007/03/26 12:30:18  mijp2
! reworked eta algorithm so now scales O(N^1.5)
!
! Revision 1.49  2006/10/12 12:17:57  mijp2
! fixed minor typo
!
! Revision 1.48  2006/10/11 08:25:53  mijp2
! changes to dipole-dipole routine (partial fix of 06283vymq01) by Keith
!
! Revision 1.47  2006/09/06 14:52:33  mijp2
! added missing variable to debug block
!
! Revision 1.46  2006/09/04 19:10:41  mijp2
! removed g95 warnings
!
! Revision 1.45  2006/09/02 12:10:36  mijp2
! updated debug and allocate error messages
!
! Revision 1.44  2006/08/10 14:48:46  mijp2
! added parallel desync check in ewald_calculate_energy
!
!Revision 1.43  2005/11/01  15:54:04  mijp2
!minor buglet fixed [Keith]
!
! Revision 1.42  2005/08/30 16:21:19  mijp2
! forchk tidy up
!
! Revision 1.41  2005/08/30 15:11:10  mijp2
! E2 updates from Keith
!
! Revision 1.40.100.4  2005/08/26 09:29:07  kr264
! Changed arg E2_matrix to ewald_calculate_E2_dipole_dipole to fullsquare dimensions as individual 3x3 blocks are NOT symmetric.
! Added additional overloaded version of ewald_calculate_E2 with square E2_matrix arg.
! Deleted unused ewald_calculate_E2_nq version.
!
! Revision 1.40.100.3  2005/08/24 08:26:22  kr264
! Modified ion->nsp,ni mapping for E2 to use new
! cell%ion_pack_* variables.
!
! Revision 1.40.100.2  2005/08/08 17:15:03  kr264
! Fixed a number of errors in e2_dipole_dipole.
! Result is not eta-invariand and Z^2/e scale invariant.
!
! Revision 1.40.100.1  2005/07/05 16:47:36  kr264
! Merged phonon_supercell_3 and main branches.
! Normalized q in ew_dipole_dipole into range [-0.5,0.5) in response to convergence issues.
!
! Revision 1.40  2005/04/06 16:07:27  mijp2
! bug fix for 05091vymq01 - tighter tolerances and better handling of very skewed cells
!
!Revision 1.39  2005/02/23  15:18:32  mijp2
!fixed debug bug in calculate_E2
!
!Revision 1.38  2005/01/06  17:00:19  mijp2
!fixed bug in E2 caused by checkin 1.36
!
! Revision 1.37  2004/11/26 08:59:22  mijp2
! tightened convergence tolerance
!
!Revision 1.36  2004/08/10  23:05:17  mijp2
!fixed way-out-of-cell real-space bug [04223e1kq01]
!
!Revision 1.35  2004/07/27  22:53:00  mijp2
!VCA now works
!
!Revision 1.34  2004/07/16  10:21:07  mijp2
!made VCA changes
!
!Revision 1.33  2004/06/11  16:13:53  mijp2
!reordered to prevent memory leak in redundant call situation
!
!Revision 1.32  2004/02/06  00:41:33  mijp2
!added checks for changins number of ions for supercell calculations (thanks Barbara)
!
! Revision 1.31  2003/04/30 22:10:32  mijp2
! Added saves to module variables
!
!Revision 1.30  2003/03/26  14:16:28  mijp2
!added missing deallocates in ewald_calculate_E2
!
!Revision 1.29  2003/03/20  10:31:08  mijp2
!FORCHK caught a missing explicit cmplx conversion and a set of redundant arrays
!
!Revision 1.28  2002/12/19  09:12:42  mijp2
!neater, shorter, faster.
!
!Revision 1.27  2002/12/19  08:07:57  mijp2
!changed spec of ewald_calculate_E2 (KR)
!
! Revision 1.26  2002/11/21 09:47:37  mijp2
! *** empty log message ***
!
!Revision 1.25  2002/10/03  13:05:21  mijp2
!fixed bug 02275vymq02 (corrected zeroing of reuse sum arrays in parallel)
!
!Revision 1.24  2002/08/05  08:44:49  mijp2
!fixed parallel bug in LR - thanks to Keith Refson
!
!Revision 1.23  2002/07/29  11:03:36  mijp2
!used g_epsilon for consistency within LR
!
!Revision 1.22  2002/06/27  14:30:21  mijp2
!Fixed ewald_calculate_E2 errors (Keith Refson):
!E2 errors which cropped up for non-cubic cells and bug which computed incorrect results for off-diaginal imaginary terms.  (Real-space part was CC).
!
! Revision 1.21  2002/04/23  12:01:06  mijp2
! bug fixes in ewald_calculate_E2 (Keith Refson)
! re-use bug fixes
! redundant variables in stress & E2 eliminated
!
!Revision 1.20  2002/04/03  13:06:01  mijp2
!tweaked error messages
!
!Revision 1.19  2002/03/04  17:08:24  mijp2
!~20% speedup in force routine
!
!Revision 1.18  2002/02/12  13:02:27  mijp2
!use cmplx_0 throughout
!
!Revision 1.17  2002/01/23  18:19:30  mijp2
!changed iprint values
!
!Revision 1.16  2001/11/17  03:18:27  mijp2
!initial LR calculation
!
!Revision 1.15  2001/10/12  22:26:47  mijp2
!deleted timers
!
!Revision 1.14  2001/09/25  11:21:16  mijp2
!replaced private erfc function with algor call
!
!Revision 1.13  2001/06/27  22:05:07  mijp2
!neatened up iprint usage and redundant variables
!
!Revision 1.12  2001/03/21  10:47:37  mijp2
!revised ordering of recip_lattice in accordance with spec change
!
!Revision 1.11  2001/03/01  15:55:03  mijp2
!typo introduced in last revision now fixed
!
!Revision 1.10  2001/03/01  11:45:20  mijp2
!type conversions simplified
!
!Revision 1.9  2001/02/27  14:59:01  mijp2
!RCS command correction
!
!                                                                             !
!=============================================================================!

module ewald
  use constants, only : DP

  implicit none                                 !Impose strong typing
  private                                       !Everything is private ...

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: ewald_calculate_energy
  public :: ewald_calculate_forces
  public :: ewald_exit
  !public :: ewald_calculate_stress
  !public :: ewald_calculate_E2
  !public :: ewald_calc_storage

  !interface ewald_calculate_E2
  !   module procedure ewald_calc_E2_qpoint_full
  !   module procedure ewald_calc_E2_dipole_dipole
  !end interface
  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  real(kind=dp), save :: eta                                     !Ewald convergence parameter (dimensions L^-1)
  real(kind=dp), save :: elect=1.0_dp                            !e^2/(4*pi*epsilon0)=1.0 in atomic units
  !real(kind=dp) :: elect=1.0000086637_dp                        !correct for incorrect castep factor
  !(0.001% accuracy in energy). NB value in MSI CASTEP v4.2 is out-of-date with current SI constants.

  real(kind=dp), save :: real_cutoff                             !radius in real space (dimensionless due to eta)

  !NB old cut-off = 4.0 gives 0.000023%    error in Ewald energy/atom for 64-atom FCC lattice w.r.t. Madelung
  !v1.37:new cut-off gives  0.00000000068% error by setting real_cutoff to 10.0 & leaving recip_cutoff at 4.0
  !The 4.0 value gives a sub-meV error for most systems but might become noticable for larger ones
  !hence new default is 10.0 as of 26/11/04 and the extra cost is negligible (1 sec for 64 atoms on my Alpha)

  !v1.40: Seems that for highly skewed systems even (10,4) is not good enough to always get meV accuracy
  !e.g. Li-primitive cell+symmetry with redefined lattice using (/1,n,0 /0,1,0/ 0,0,1/) gives
  !different answers for all n in v1.39 (skewing effects choice of eta with old algorithm)
  !and -2.353 eV change in Ewald energy (n=0 -> n=9)
  !Using min() rather than geometric mean for eta gives different answers for n>2 and max error now -0.006 eV
  !Increasing recip_cutoff to 10 increases all Ewald energies by 6*10^-7 eV and no change in max error
  !Setting both cutoffs to 20 makes no difference in Ewald energies for small n and reduces max error to
  !6*10^-5 eV but DOES have an impact on runtime (takes 64-atom FCC from 5 -> 20 secs) which I guess
  !is negligible in the overall scheme of things.

  !v1.50: Now got a smarter system: real_cutoff is no longer hard-wired but can adapt to the situation
  !with a new algorithm for eta that changes according to the level of precision (convergence of the Ewald sum)
  !that is required. After some experimentation, the defaults have been set. The devel_code is left in
  !in case a pathological case is found in the future that requires higher precision.
  !As it stands, should get accuracy ~10^-16 with precision=36.
  !Algorithm follows Fincham (CCP5 newsletter 1993)
  !and now that eta varies with num_atoms the overall
  !algorithm has a scaling~O(N^3/2) not O(N^2) as before.
  !Values set in ewald_calculate_num_cells ...

  !v1.57 Fixed typo in Fincham algorithm
  !v1.58 Consequently re-optimised value of precision and real-to-reciprocal space calculation time,
  !for some large systems (too fast for small ones!) using trace_sections and ewald_test.f90
  !So, for 579-atom alpha-quartz system, time to calculate Ewald energy+force+stress has gone from
  !1974 secs (v1.49) to 79 secs (v1.58) on my Alpha, with 5*10^(-8)% change in Ewald energy.
  !v1.63 Revisited again! Another user has found a very skewed system that gives a different answer
  !depending on skew (even with old algorithm) so decided to tighten defaults up to precision=1000
  !and stop this nonsense, and then let "power users" turn it down if they need to.
  !v1.64 Bah! Forgot memory implications. Reset precision=100

  real(kind=dp), parameter :: g_epsilon =  1.0e-10_dp

  real(kind=dp), parameter :: real_min_len    = g_epsilon        !avoid self-self interactions
  real(kind=dp), parameter :: recip_min_len   = g_epsilon        !avoid self-self interactions
  !NB we now use g_epsilon from constants for consistency with rest of the code (particular for LR bits)

  real(kind=dp), save :: real_cutoff_sq
  real(kind=dp), save :: real_min_len_sq
  real(kind=dp), save :: recip_cutoff                 !radius in recip space (dimensionless due to eta)
  real(kind=dp), save :: recip_cutoff_sq
  real(kind=dp), save :: recip_min_len_sq

  integer, save :: max_real_cells_x                   !max number of real-space cells in x-direction
  integer, save :: max_real_cells_y                   !max number of real-space cells in y-direction
  integer, save :: max_real_cells_z                   !max number of real-space cells in z-direction

  integer, save :: num_real_cells_x                   !total number of real-space cells in x-direction
  integer, save :: num_real_cells_y                   !total number of real-space cells in y-direction
  integer, save :: num_real_cells_z                   !total number of real-space cells in z-direction

  integer, save :: num_real_cells                     !total number of real-space cells

  integer, save :: max_recip_cells_x                  !max number of reciprocal-space cells in x-direction
  integer, save :: max_recip_cells_y                  !max number of reciprocal-space cells in y-direction
  integer, save :: max_recip_cells_z                  !max number of reciprocal-space cells in z-direction

  integer, save :: num_recip_cells_x                  !total number of reciprocal-space cells in x-direction
  integer, save :: num_recip_cells_y                  !total number of reciprocal-space cells in y-direction
  integer, save :: num_recip_cells_z                  !total number of reciprocal-space cells in z-direction

  integer, save :: num_recip_cells=0                  !total number of reciprocal-space cells

  !lookup tables
  integer,        allocatable, dimension(:,:), save :: recip_nxyz !packed -> unpacked cell indices
  real (kind=dp), allocatable, dimension(:,:), save :: real_dr    !real space cell displacements
  real (kind=dp), allocatable, dimension(:,:), save :: recip_dg   !recip space cell displacements

  !3Dcorrection
  !logical, save :: EW3DC=.false.

contains
  subroutine ewald_calculate_energy(elements, ewald_energy, finished, &
       classical_elements)
    !=========================================================================!
    ! Calculate the Ewald energy of current_cell                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   ewald_energy, intent=out, the Ewald energy in atomic units            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   all of them!                                                          !
    !   If lattice_changed, then will call ewald_calculate_num_cells to get   !
    !      new optimal number of cells in each direction, real & reciprocal.  !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   cell for current_cell information                                     !
    !   comms for parallelization to speed up where possible.                 !
    !   io for output and unit conversion routines                            !
    !   parameters for verbosity and units                                    !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   The real space sum is split into two contributions - self & non-self. !
    !   The self-term is independant of the atomic coords and only depends on !
    !       the lattice vectors - hence we can store it and reuse as long as  !
    !       the cell has not changed.                                         !
    !   The non-self term has to be calculated afresh every time, unless the  !
    !       atoms have not moved between calls!                               !
    !   Ditto for reciprocal space sum, with the addition that there is an    !
    !       array of G-vector expressions that can be stored and reused as    !
    !       long as the cell has not changed.                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   current_cell must have valid data for atomic coords, cell vectors and !
    !       ionic charges.                                                    !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.1, 05/07/2000                               !
    !=========================================================================!
    use comms, only: comms_abort, comms_free, comms_reduce, pub_on_root, &
         pub_my_node_id, pub_total_num_nodes
    use constants, only: dp,stdout, NORMAL, two_pi, stderr, PI, cmplx_i, cmplx_0
    use ion, only: element
    use rundat, only: pub_output_detail
    use simulation_cell, only: pub_cell, castep_cell,castep_cell_dealloc,&
         copy_to_castep_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_erfc

    implicit none

    type(element), dimension(1:pub_cell%nat), intent(in) :: elements
    real (kind=dp), intent(out) :: ewald_energy
    logical, intent(in) :: finished
    type(element), dimension(1:pub_cell%nat_classical), &
         intent(in), optional :: classical_elements

    type(castep_cell) :: current_cell

    !local bits of the energy
    real   (kind=dp)       :: real_ewald_energy            !total real-space sum contribution
    real   (kind=dp)       :: recip_ewald_energy           !total recip-space sum contribution
    real   (kind=dp)       :: dE                           !working contribution

    !reusable bits
    logical, save                                    :: lattice_changed  !save time by not repeating bits needlessly
    real(kind=dp), dimension(1:3,1:3), save          :: old_real_lattice !store old lattice to check if changed
    real(kind=dp), dimension(:), allocatable, save   :: recip_gauss      !actually exp(-G*G)/(G*G), stored for reuse.
    logical, save                                    :: ions_moved       !save time by not recomputing if nothing moved!
    real(kind=dp), dimension(:,:), allocatable, save :: old_ion_position !store old ionic coords in case of really dumb usage.
    real(kind=dp), save                              :: old_ewald_energy !store old answer        "  "    "    "     "    "
    integer, save                                    :: old_num_ions     !store old number of ions
    logical, save                                    :: first_pass=.true.
    logical, save                                    :: num_ions_changed !save time by not recomputing if number of ions has not changed!

    !lengths
    real   (kind=dp) :: real_len                           !dimensionless length=|r1-r2-R!*eta
    real   (kind=dp) :: real_len_sq                        !real_len**2
    real   (kind=dp) :: recip_len_sq                       ![dimensionless length=|G|/(2*eta)] **2

    real   (kind=dp) :: xdiff,ydiff,zdiff                  !(r1-r2-R)x,y,z length
    real   (kind=dp) :: gx,gy,gz                           ! Gx/y/z recip-vec length
    real   (kind=dp), dimension(1:3) :: a_diff             !displacement vec in arbitary space
    real   (kind=dp) :: rx0,ry0,rz0                        !offset origin for real-space sums
    real   (kind=dp) :: gx0,gy0,gz0                        !offset origin for recip-space sums

    !cell phases
    complex(kind=dp), dimension(:), allocatable :: exp_recip_phase_x !exp(i(G.dR)x)
    complex(kind=dp), dimension(:), allocatable :: exp_recip_phase_y
    complex(kind=dp), dimension(:), allocatable :: exp_recip_phase_z
    complex(kind=dp) :: recip_phase,exp_recip_phase

    !ionic stuff
    real(kind=dp), dimension(:), allocatable :: ion_charge !packed array of charge of each ion in current cell
    !integer :: iion,imix
    !real(kind=dp) :: mix_charge
    real(kind=dp) :: tot_ion_charge                        !sum of (ion_charge)
    real(kind=dp) :: tot_ion_charge_sq                     !sum of (ion_charge)**2
    real(kind=dp) :: charge_ij                             !ion_charge(i)*ion_charge(j) - real for convenience
    real(kind=dp), dimension(:,:), allocatable :: ion_position !packed array of positions of each ion in current cell

    !parallelization stuff
    integer       :: num_active_nodes                      !<> pub_total_num_nodes if we have any idle nodes!

    !miscellaneous scale factors
    real(kind=dp) :: scale_real_energy                     !converts from dimensionless lengths and other units to Hartrees.
    real(kind=dp) :: scale_recip_energy                    ! ditto
    real(kind=dp) :: inv_two_eta_sq                        !1/(2*eta)**2
    complex(kind=dp) :: two_i_pi                           !2*sqrt(-1)*pi

    !miscellaneous loop indices, etc.
    integer :: ispec,iatom                                 !unpacked array indices
    integer :: ion_i,ion_j                                 !packed array indices
    integer :: nx,ny,nz                                    !loop indices in x/y/z directions
    integer :: nxyz                                        !combined loop index for x & y & z directions
    integer :: ng                                          !index for packed arrays of reciprocal lattice vectors
    integer :: idim,jdim                                   !3D loop indices

    integer :: ierr
    !integer                         :: status

    !pretty printing
    character(len=*), parameter :: energy_label = 'Hartree'

    ! qoh: Clean up by deallocating saved arrays
    if (finished)  then
       if (allocated(recip_gauss)) then
          deallocate(recip_gauss,stat=ierr)
          call utils_dealloc_check('ewald_calculate_energy','recip_gauss',ierr)
       end if
       if (allocated(old_ion_position)) then
          deallocate(old_ion_position,stat=ierr)
          call utils_dealloc_check('ewald_calculate_energy','old_ion_position',&
               ierr)
       end if
       ! ndmh: deallocate arrays allocated in ewald_calculate_num_cells
       if (allocated(recip_nxyz)) then
          deallocate(recip_nxyz,stat=ierr)
          call utils_dealloc_check('ewald_calculate_num_cells','recip_nxyz',&
               ierr)
       end if
       if (allocated(recip_dg)) then
          deallocate(recip_dg,stat=ierr)
          call utils_dealloc_check('ewald_calculate_num_cells','recip_dg',&
               ierr)
       end if
       if (allocated(real_dr)) then
          deallocate(real_dr,stat=ierr)
          call utils_dealloc_check('ewald_calculate_num_cells','real_dr',&
               ierr)
       end if
       ! ndmh: flag to reallocate arrays next time
       first_pass = .true.
       return
    end if

    ! ndmh: start timer
    call timer_clock("ewald_calculate_energy",1)

    !First copy data from ONES structures into CASTEP cell structure
    !This is an inefficient hack, but it will be the most reliable way...
    if (pub_cell%nat_classical.ne.0) then
      call copy_to_castep_cell(current_cell,elements,classical_elements)
    else
      call copy_to_castep_cell(current_cell,elements)
    end if

    if (current_cell%num_ions<=0) then
       write (stderr,'(a)') 'Error in ewald_calculate_energy - &
            &current_cell uninitialised'
       call comms_abort
    end if

    !Set up parallelization over all active nodes
    num_active_nodes=pub_total_num_nodes
    !setup packed arrays for charge and position for efficiency and parallelisation
    allocate(ion_charge(1:current_cell%num_ions),stat=ierr)
    call utils_alloc_check('ewald_calculate_energy','ion_charge',ierr)
    allocate(ion_position(1:3,1:current_cell%num_ions),stat=ierr)
    call utils_alloc_check('ewald_calculate_energy','ion_position',ierr)

    !check if cell has changed (or first pass)
    lattice_changed=.false.
    num_ions_changed=.false.
    if (first_pass) then
       old_num_ions=current_cell%num_ions
       old_real_lattice=current_cell%real_lattice
       allocate(old_ion_position(1:3,1:current_cell%num_ions),stat=ierr)
       call utils_alloc_check('ewald_calculate_energy','old_ion_position',ierr)
       old_ion_position=0.0_dp
       first_pass=.false.
       lattice_changed=.true.
       ions_moved=.true.
       num_ions_changed=.true.
       old_ewald_energy=0.0_dp
    else
       do idim=1,3
          do jdim=1,3
             if (old_real_lattice(jdim,idim)/=current_cell%real_lattice(jdim,idim)) lattice_changed=.true.
          end do
       end do
    end if

    ! Need to synchronise as parallel distribution may have changed and old values differ on previously idle nodes

    !call comms_reduce_kp(lattice_changed,1,'OR')
    !call comms_reduce_bnd(lattice_changed,1,'OR')
    !call comms_reduce_gv(lattice_changed,1,'OR')
    call comms_reduce('OR',lattice_changed)

    !calculate packed array of ion charges and positions and total sums
    !tot_ion_charge=0.0_dp
    !tot_ion_charge_sq=0.0_dp
    !do iion=1,current_cell%num_ions
    !   ispec=current_cell%ion_pack_species(iion)
    !   ion_charge(iion)=current_cell%ionic_charge(ispec)
    !end do
    !and allow for VCA possibility
    !=> assign appropriate weighted charge to each compound atom
    !do imix = 1,current_cell%num_mixture_atoms
    !   mix_charge = 0.0_dp
    !   do ion_i = 1,current_cell%indx_mixture_atoms(imix)
    !      iion  = current_cell%list_mixture_atoms(ion_i,imix)
    !      iatom = current_cell%ion_pack_index(iion)
    !      ispec = current_cell%ion_pack_species(iion)
    !      mix_charge = mix_charge + current_cell%ionic_charge(ispec)*current_cell%mixture_weight(iatom,ispec)
    !   end do
    !   do ion_i = 1,current_cell%indx_mixture_atoms(imix)
    !      iion  = current_cell%list_mixture_atoms(ion_i,imix)
    !      if (ion_i==1) then
    !         ion_charge(iion) = mix_charge
    !      else
    !         ion_charge(iion) = 0.0_dp !label this atom as a VCA fake
    !      end if
    !   end do
    !end do
    !do ion_i=1,current_cell%num_ions
    !   iatom = current_cell%ion_pack_index(ion_i)
    !   ispec = current_cell%ion_pack_species(ion_i)
    !   tot_ion_charge    =tot_ion_charge+ion_charge(ion_i)
    !   tot_ion_charge_sq =tot_ion_charge_sq+ion_charge(ion_i)**2
    !   do idim=1,3
    !      ion_position(idim,ion_i)=current_cell%ionic_positions(idim,iatom,ispec) &
    !           & -nint(current_cell%ionic_positions(idim,iatom,ispec))
    !   end do
    !end do

    !calculate packed array of ion charges and positions and total sums
    tot_ion_charge=0.0_dp
    tot_ion_charge_sq=0.0_dp
    ion_i=1
    do ispec=1,current_cell%num_species
       do iatom=1,current_cell%num_ions_in_species(ispec)
          ion_charge(ion_i) =current_cell%ionic_charge(ispec)
          tot_ion_charge    =tot_ion_charge+ion_charge(ion_i)
          tot_ion_charge_sq =tot_ion_charge_sq+ion_charge(ion_i)**2
          do idim=1,3
             ion_position(idim,ion_i)=current_cell%ionic_positions(idim,iatom,ispec)
          end do
          ion_i=ion_i+1
       end do
    end do

    !check to make sure that the ions have actually moved
    ions_moved=.false.
    num_ions_changed=.false.
    if(current_cell%num_ions/=old_num_ions) then
       deallocate(old_ion_position,stat=ierr)
       call utils_dealloc_check('ewald_calculate_energy','old_ion_position',ierr)
       allocate(old_ion_position(1:3,1:current_cell%num_ions),stat=ierr)
       call utils_alloc_check('ewald_calculate_energy','old_ion_position',ierr)
       old_ion_position=0.0_dp
       num_ions_changed=.true.
    else
       do ion_i=1,current_cell%num_ions
          do idim=1,3
             if (old_ion_position(idim,ion_i)/=ion_position(idim,ion_i)) ions_moved=.true.
          end do
       end do
    end if

    ! Need to synchronise as parallel distribution may have changed and old values differ on previously idle nodes

    !call comms_reduce_kp(ions_moved,1,'OR')
    !call comms_reduce_bnd(ions_moved,1,'OR')
    !call comms_reduce_gv(ions_moved,1,'OR')
    call comms_reduce('OR',ions_moved)

    !check if we actually need to do a calculation or not
    if ((.not.lattice_changed).and.(.not.ions_moved).and.(.not.num_ions_changed)) then
       ewald_energy=old_ewald_energy
       if (pub_output_detail>NORMAL.and.pub_on_root) then
          write (stdout,*) 'NB: reusing old value of ewald_energy as nothing has changed!'
       end if
       deallocate(ion_charge,stat=ierr)
       call utils_dealloc_check('ewald_calculate_energy','ion_charge',ierr)
       deallocate(ion_position,stat=ierr)
       call utils_dealloc_check('ewald_calculate_energy','ion_position',ierr)
       !call trace_exit('ewald_calculate_energy',status)
       call castep_cell_dealloc(current_cell)
       call timer_clock("ewald_calculate_energy",2)
       return
    end if

    !Calculate eta and number of cells needed for real and reciprocal space sums
    if (lattice_changed.or.num_ions_changed) call ewald_calculate_num_cells(current_cell)

    !... so can now dimension allocatable arrays
    if (lattice_changed.or.num_ions_changed) then
       if (allocated(recip_gauss)) then
         deallocate(recip_gauss,stat=ierr)
         call utils_dealloc_check('ewald_calculate_energy','recip_gauss',ierr)
       end if
       allocate(recip_gauss(1:num_recip_cells),stat=ierr)
       call utils_alloc_check('ewald_calculate_energy','recip_gauss',ierr)
    end if     !lattice changed

    if (.not.allocated(recip_gauss)) then
       allocate(recip_gauss(1:num_recip_cells),stat=ierr)
       call utils_alloc_check('ewald_calculate_energy','recip_gauss',ierr)
    end if
    allocate(exp_recip_phase_x(1:num_recip_cells_x),stat=ierr)
    call utils_alloc_check('ewald_calculate_energy','exp_recip_phase_x',ierr)
    allocate(exp_recip_phase_y(1:num_recip_cells_y),stat=ierr)
    call utils_alloc_check('ewald_calculate_energy','exp_recip_phase_y',ierr)
    allocate(exp_recip_phase_z(1:num_recip_cells_z),stat=ierr)
    call utils_alloc_check('ewald_calculate_energy','exp_recip_phase_z',ierr)

    !and initialize to zero
    ewald_energy=0.0_dp
    real_ewald_energy=0.0_dp
    recip_ewald_energy=0.0_dp
    exp_recip_phase_x=cmplx_0
    exp_recip_phase_y=cmplx_0
    exp_recip_phase_z=cmplx_0
    if (lattice_changed.or.num_ions_changed) then
       recip_gauss=0.0_dp
    end if

    !miscellaneous scaling factors
    two_i_pi=two_pi*cmplx_i
    inv_two_eta_sq=1.0_dp/(2.0_dp*eta)**2
    scale_real_energy =eta*elect
    scale_recip_energy=(pi*elect)/(current_cell%volume*eta**2)

    !=========================================================================!
    !calculate real-space contribution to Ewald energy
    !=========================================================================!

    !VCA - had to abandon the splitting into self and non-self terms as no longer have equivalent ions
    !      for a given species so simplest to just ignore this minor speedup
    !      and just remember to allow for double counting of the self-term in the general sum

    !Now do the ion(i)->ion(j) interaction and implicitly fill in the (j,i) terms by symmetry using charge_ij factor ,..
    !parallelizing outer loop over all active nodes
    do ion_i=pub_my_node_id+1,current_cell%num_ions,num_active_nodes
       do ion_j=ion_i,current_cell%num_ions

          charge_ij=ion_charge(ion_i)*ion_charge(ion_j)
          if (ion_i==ion_j) charge_ij=charge_ij/2.0_dp

          !initialise the distances in the x,y,z directions between ion(i) and the periodic images of ion_j
          !in the corner cell of the block of cells included in the real-space part of the Ewald sum
          a_diff(1)=ion_position(1,ion_i)-ion_position(1,ion_j)+real(max_real_cells_x,kind=dp)
          a_diff(2)=ion_position(2,ion_i)-ion_position(2,ion_j)+real(max_real_cells_y,kind=dp)
          a_diff(3)=ion_position(3,ion_i)-ion_position(3,ion_j)+real(max_real_cells_z,kind=dp)

          rx0=dot_product(current_cell%real_lattice(:,1),a_diff(:))
          ry0=dot_product(current_cell%real_lattice(:,2),a_diff(:))
          rz0=dot_product(current_cell%real_lattice(:,3),a_diff(:))

          do nxyz=1,num_real_cells
             xdiff=rx0-real_dr(1,nxyz)
             ydiff=ry0-real_dr(2,nxyz)
             zdiff=rz0-real_dr(3,nxyz)

             if (abs(xdiff)<real_min_len) xdiff=0.0_dp    !catch underflow
             if (abs(ydiff)<real_min_len) ydiff=0.0_dp
             if (abs(zdiff)<real_min_len) zdiff=0.0_dp

             real_len_sq=(xdiff**2+ydiff**2+zdiff**2)*(eta**2)
             !NB Might be tempted to remove this if-block on a vector machine, thinking that the cost of the extra
             !calculations would be offset by vectorizing this inner loop.
             !But, whilst this may be computationally more efficient, it would be PHYSICALLY WRONG!
             !We need a spherical cutoff to get reliable convergence of the sum,
             !and to get a proper treatment of the surface term. mijp2
             if (real_len_sq>real_min_len_sq.and.real_len_sq<real_cutoff_sq) then
                real_len=sqrt(real_len_sq)
!                dE=charge_ij*utils_custom_erfc(real_len)/real_len
                dE=charge_ij*utils_erfc(real_len)/real_len
                real_ewald_energy=real_ewald_energy+dE
             end if

          end do

          !move onto next ion_j
       end do

       !move onto next ion_i
    end do

    !and gather up over all the nodes and convert to atomic units
    !call comms_reduce_gv(real_ewald_energy,1,'SUM')
    !call comms_reduce_bnd(real_ewald_energy,1,'SUM')
    !call comms_reduce_kp(real_ewald_energy,1,'SUM')
    call comms_reduce('SUM',real_ewald_energy)
    real_ewald_energy=real_ewald_energy*scale_real_energy
    if (pub_output_detail>NORMAL.and.pub_on_root) then
       write (stdout,*)
       write (stdout,99) 'Ewald: real-space energy =', real_ewald_energy
!         & io_atomic_to_unit(real_ewald_energy,energy_unit),trim(energy_label)
    end if

    !... finished real space term

    !=========================================================================!
    !calculate reciprocal-space contribution to Ewald energy
    !=========================================================================!

    !First we precompute various factors ...
    if (lattice_changed.or.num_ions_changed) then
       a_diff(1)=real(max_recip_cells_x,kind=dp)
       a_diff(2)=real(max_recip_cells_y,kind=dp)
       a_diff(3)=real(max_recip_cells_z,kind=dp)

       gx0=-dot_product(current_cell%recip_lattice(:,1),a_diff(:))
       gy0=-dot_product(current_cell%recip_lattice(:,2),a_diff(:))
       gz0=-dot_product(current_cell%recip_lattice(:,3),a_diff(:))

       !parallelizing outer loop over all active nodes
       !NB This way of doing things with the single loop and the lookup table is slightly faster than the previous triple loop
       !   AND, more significantly, it makes a longer loop => better load balancing with more processors. mijp.
       do ng=pub_my_node_id+1,num_recip_cells,num_active_nodes
          gx=gx0+recip_dg(1,ng)
          gy=gy0+recip_dg(2,ng)
          gz=gz0+recip_dg(3,ng)

          recip_len_sq=(gx**2+gy**2+gz**2)*inv_two_eta_sq
          !NB Might be tempted to remove this if-block on a vector machine, thinking that the cost of the extra
          !calculations would be offset by vectorizing this inner loop.
          !But, whilst this may be computationally more efficient, it would be PHYSICALLY WRONG!
          !We need a spherical cutoff to get reliable convergence of the sum,
          !and to get a proper treatment of the surface term. mijp2
          if(recip_len_sq>recip_min_len_sq.and.recip_len_sq<recip_cutoff_sq) then
             recip_gauss(ng)=exp(-recip_len_sq)/recip_len_sq
          else
             recip_gauss(ng)=0.0_dp
          endif

       end do

       !and we need recip_gauss on all nodes, so gather it up ...
       !call comms_reduce_gv(recip_gauss(1),num_recip_cells,'SUM')
       !call comms_reduce_bnd(recip_gauss(1),num_recip_cells,'SUM')
       !call comms_reduce_kp(recip_gauss(1),num_recip_cells,'SUM')
       call comms_reduce('SUM',recip_gauss)

    end if !lattice_changed.or.num_ions_changed

    !VCA - had to abandon the splitting into self and non-self terms as no longer have equivalent ions
    !      for a given species so simplest to just ignore this minor speedup
    !      and just remember to allow for double counting of the self-term in the general sum

    !... and now do the reciprocal space sum over ion(i)->ion(j) contributions
    do ion_i=pub_my_node_id+1,current_cell%num_ions,num_active_nodes
       do ion_j=ion_i,current_cell%num_ions

          charge_ij=ion_charge(ion_i)*ion_charge(ion_j)
          if (ion_i==ion_j) charge_ij=charge_ij/2.0_dp

          !precompute x-direction phase factors=exp(i.Gx.(r1x-r2x))
          recip_phase=(ion_position(1,ion_i)-ion_position(1,ion_j))*two_i_pi
          exp_recip_phase=exp(recip_phase)
          exp_recip_phase_x(1)=exp(-recip_phase*max_recip_cells_x)
          do nx=2,num_recip_cells_x
             exp_recip_phase_x(nx)=exp_recip_phase_x(nx-1)*exp_recip_phase
          end do

          !precompute y-direction phase factors
          recip_phase=(ion_position(2,ion_i)-ion_position(2,ion_j))*two_i_pi
          exp_recip_phase=exp(recip_phase)
          exp_recip_phase_y(1)=exp(-recip_phase*max_recip_cells_y)
          do ny=2,num_recip_cells_y
             exp_recip_phase_y(ny)=exp_recip_phase_y(ny-1)*exp_recip_phase
          end do

          !precompute z-direction phase factors
          recip_phase=(ion_position(3,ion_i)-ion_position(3,ion_j))*two_i_pi
          exp_recip_phase=exp(recip_phase)
          exp_recip_phase_z(1)=exp(-recip_phase*max_recip_cells_z)
          do nz=2,num_recip_cells_z
             exp_recip_phase_z(nz)=exp_recip_phase_z(nz-1)*exp_recip_phase
          end do

          do ng=1,num_recip_cells
             nx=recip_nxyz(1,ng)
             ny=recip_nxyz(2,ng)
             nz=recip_nxyz(3,ng)
             exp_recip_phase=exp_recip_phase_x(nx)*exp_recip_phase_y(ny)*exp_recip_phase_z(nz)
             recip_ewald_energy=recip_ewald_energy+recip_gauss(ng)*real(exp_recip_phase,kind=dp)*charge_ij
          end do

          !move onto the next ion_j
       end do

       !move onto the next ion_i
    end do

    !and gather up over all the nodes and convert to atomic units
    !call comms_reduce_gv(recip_ewald_energy,1,'SUM')
    !call comms_reduce_bnd(recip_ewald_energy,1,'SUM')
    !call comms_reduce(recip_ewald_energy,1,'SUM')
    call comms_reduce('SUM',recip_ewald_energy)
    recip_ewald_energy=recip_ewald_energy*scale_recip_energy

    if (pub_output_detail>NORMAL.and.pub_on_root) write (stdout,99) 'Ewald: reciprocal-space energy =', &
         recip_ewald_energy
!         & io_atomic_to_unit(recip_ewald_energy,energy_unit),trim(energy_label)

    !... finished reciprocal space term

    !... finally add real and reciprocal space terms and the other bits
    !(due to uniform charge-neutralising background and boundary conditions)
    ewald_energy=real_ewald_energy+recip_ewald_energy &
         & -(0.5_dp*elect*pi*tot_ion_charge**2/(current_cell%volume*eta**2)) &
         & -(eta*elect*tot_ion_charge_sq/sqrt(pi))

    if (pub_output_detail>NORMAL.and.pub_on_root) then
       write (stdout,99) 'Ewald: (sum_charge)**2 Energy =', &
            & 0.5_dp*elect*pi*tot_ion_charge**2/(current_cell%volume*eta**2)
!            & io_atomic_to_unit(0.5_dp*elect*pi*tot_ion_charge**2/(current_cell%volume*eta**2),energy_unit),trim(energy_label)
       write (stdout,99) 'Ewald: (sum_charge**2) Energy =', &
            & eta*elect*tot_ion_charge_sq/sqrt(pi)
!            & io_atomic_to_unit(eta*elect*tot_ion_charge_sq/sqrt(pi),energy_unit),trim(energy_label)
       write (stdout,99) 'Ewald Energy =', &
            & ewald_energy
!            & io_atomic_to_unit(ewald_energy,energy_unit),trim(energy_label)
    else if (pub_on_root) then
       write(stdout,'(f20.6,1x,a)') ewald_energy,trim(energy_label)
    end if

    !store result for next time ...
    old_ewald_energy=ewald_energy
    old_num_ions=current_cell%num_ions
    old_ion_position=ion_position
    old_real_lattice=current_cell%real_lattice

    !... clean up temporary arrays ...
    deallocate(exp_recip_phase_x,stat=ierr)
    call utils_dealloc_check('ewald_calculate_energy','exp_recip_phase_x',ierr)
    deallocate(exp_recip_phase_y,stat=ierr)
    call utils_dealloc_check('ewald_calculate_energy','exp_recip_phase_y',ierr)
    deallocate(exp_recip_phase_z,stat=ierr)
    call utils_dealloc_check('ewald_calculate_energy','exp_recip_phase_z',ierr)

    deallocate(ion_charge,stat=ierr)
    call utils_dealloc_check('ewald_calculate_energy','ion_charge',ierr)
    deallocate(ion_position,stat=ierr)
    call utils_dealloc_check('ewald_calculate_energy','ion_position',ierr)

    call castep_cell_dealloc(current_cell)

99  format(A40,F15.6)

    ! Free up comms memory
    call comms_free

    !add in 3D correction?
    !if (EW3DC) call ewald_calculate_EW3DC_energy(ewald_energy)

    !call trace_exit('ewald_calculate_energy',status)
    call timer_clock("ewald_calculate_energy",2)

    return
  end subroutine ewald_calculate_energy


  subroutine ewald_calculate_forces(elements, ewald_forces, finished)
    !=========================================================================!
    ! Calculate the Ewald forces on the ions in the current_cell              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   ewald_forces, intent=out, the Ewald forces in atomic units            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   all of them!                                                          !
    !   If lattice_changed or num_ions_changed, then will call                !
    !   ewald_calculate_num_cells to get the new optimal number of cells in   !
    !   each direction, real & reciprocal.                                    !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   cell for current_cell information                                     !
    !   comms for parallelization to speed up where possible.                 !
    !   io for output and unit conversion routines                            !
    !   parameters for verbosity and units                                    !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   We use packed arrays of ionic positions/forces here to get long loops !
    !       which we then parallelize.  Using unpacked arrays (ions,species)  !
    !       is not so efficient in this respect.                              !
    !   Of course, the result is returned in unpacked form!                   !
    !                                                                         !
    !   With the reciprocal space sum, there is an array of G-vector          !
    !       expressions that can be stored and reused aslong as the cell has  !
    !       not changed.                                                      !
    !NB There is no self-force, and so we do not split force up, unlike energy!
    !   However we do make use of inversion symmetry, calculating Fij as a sum!
    !       over cells and then filling in Fji=-Fij as a transpose at end.    !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   current_cell must have valid data for atomic coords, cell vectors and !
    !       ionic charges.                                                    !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.1, 05/07/2000                               !
    ! Modified for ONES code by Arash Mostofi.                                !
    !=========================================================================!
    use comms, only: comms_abort, comms_free, comms_reduce, pub_on_root, &
         pub_my_node_id, pub_total_num_nodes
    use constants, only: DP, stdout, two_pi, NORMAL, PI, stderr, cmplx_i, cmplx_0
    use ion, only: element
    use rundat, only: pub_output_detail, pub_is_smeared_ion_rep
    use simulation_cell, only: pub_cell, castep_cell,castep_cell_dealloc,&
         copy_to_castep_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_erfc
    !kaw: Allows access to coordinates of classical atoms
    use classical_pot, only: classical_elements

    implicit none

    type(ELEMENT), dimension(1:pub_cell%nat), intent(in) :: elements
    real (kind=dp), dimension(1:3,1:pub_cell%nat), intent(out) :: ewald_forces
    logical, intent(in) :: finished

    type(castep_cell) :: current_cell

    !local bits of the force
    !real (kind=dp), dimension(:,:), allocatable       :: tmp_ewald_force    !real or recip-space, Fxij etc.
    !real (kind=dp), dimension(1:3)                    :: tmp_force          !speedup array
    real (kind=dp)                                    :: dF                 !working contribution

    !reusable bits
    !logical, save                                     :: lattice_changed    !save time by not repeating bits needlessly
    !real(kind=dp), dimension(1:3,1:3), save           :: old_real_lattice   !store old lattice to check if changed
    !logical, save                                     :: ions_moved         !save time by not recomputing if nothing moved!
    !real(kind=dp), dimension(:,:), allocatable, save  :: old_ion_position   !store old ionic coords in case of really dumb usage.
    !real(kind=dp), dimension(:,:), allocatable, save  :: old_ewald_forces   !store old answer        "  "    "    "     "    "
    !integer, save                                     :: old_num_ions       !store old number of ions
    !logical, save                                     :: first_pass=.true.
    !logical, save                                     :: num_ions_changed   !save time by not recomputing if number of ions has not changed!

    !lengths
    real (kind=dp) :: real_len                           !dimensionless length=|r1-r2-R!*eta
    real (kind=dp) :: real_len_sq                        !real_len**2
    real (kind=dp) :: recip_len_sq                       ![dimensionless length=|G|/(2*eta)] **2

    real (kind=dp) :: xdiff,ydiff,zdiff                  !(r1-r2-R)x,y,z length
    real (kind=dp) :: gx,gy,gz                           ! Gx/y/z recip-vec length
    real (kind=dp), dimension(1:3) :: a_diff             !displacement vec in arbitary space
    real (kind=dp) :: rx0,ry0,rz0                        !offset origin for real-space sums
    real (kind=dp) :: gx0,gy0,gz0                        !offset origin for recip-space sums

    ! pre-calculated values of Gx*exp(-G*G)/(G*G),etc.
    real (kind=dp), dimension(:,:), allocatable :: recip_gauss_force

    !cell phases
    complex(kind=dp), dimension(:), allocatable :: exp_recip_phase_x !exp(i(G.dR)x)
    complex(kind=dp), dimension(:), allocatable :: exp_recip_phase_y
    complex(kind=dp), dimension(:), allocatable :: exp_recip_phase_z
    complex(kind=dp) :: recip_phase,exp_recip_phase

    !ionic stuff
    real(kind=dp), dimension(:), allocatable :: ion_charge !packed array of charge of each ion in current cell
    real(kind=dp) :: charge_ij                             !ion_charge(i)*ion_charge(j) - real for convenience
    real(kind=dp), dimension(:,:), allocatable :: ion_position !packed array of positions of each ion in current cell

    !parallelization stuff
    integer       :: num_active_nodes                      !<> pub_total_num_nodes if we have any idle nodes!

    !miscellaneous scale factors
    real(kind=dp) :: scale_real_force                      !converts from dimensionless lengths and other units to Hartrees.
    real(kind=dp) :: scale_recip_force                     ! ditto
    real(kind=dp) :: two_by_root_pi                        !2/sqrt(pi)
    real(kind=dp) :: inv_two_eta_sq                        !1/(2*eta)**2
    complex(kind=dp) :: two_i_pi                           !2*sqrt(-1)*pi

    !miscellaneous loop indices, etc.
    integer :: ispec,iatom                                 !unpacked array indices
    integer :: ion_i,ion_j,ion_k                           !packed array indices
    integer :: nx,ny,nz                                    !loop indices in x/y/z directions
    integer :: nxyz                                        !combined loop index for x & y & z directions
    integer :: ng                                          !index for packed arrays of reciprocal lattice vectors
    integer :: idim                                       !3D loop index

    integer :: ierr

    !pretty printing
    !character(len=*), parameter :: force_label = 'Hartree/bohr'

    !Set up parallelization over all active nodes
    !integer                         :: status

    ! qoh: Clean up by deallocating saved arrays
    if (finished)  then
    !   if (allocated(old_ion_position)) then
    !      deallocate(old_ion_position,stat=ierr)
    !      call utils_dealloc_check('ewald_calculate_forces','old_ion_position',&
    !           ierr)
    !   end if
    !   if (allocated(old_ewald_forces)) then
    !      deallocate(old_ewald_forces,stat=ierr)
    !      call utils_dealloc_check('ewald_calculate_forces','old_ewald_forces',&
    !           ierr)
    !   end if
       ! ndmh: deallocate arrays allocated in ewald_calculate_num_cells
       if (allocated(recip_nxyz)) then
          deallocate(recip_nxyz,stat=ierr)
          call utils_dealloc_check('ewald_calculate_num_cells','recip_nxyz',&
               ierr)
       end if
       if (allocated(recip_dg)) then
          deallocate(recip_dg,stat=ierr)
          call utils_dealloc_check('ewald_calculate_num_cells','recip_dg',&
               ierr)
       end if
       if (allocated(real_dr)) then
          deallocate(real_dr,stat=ierr)
          call utils_dealloc_check('ewald_calculate_num_cells','real_dr',&
               ierr)
       end if
       ! ndmh: flag to reallocate arrays next time
       !first_pass = .true.
       return
    end if


    call timer_clock("ewald_calculate_forces",1)

    ! jd: Sanity check
    call utils_assert(.not. pub_is_smeared_ion_rep,'Ewald forces not supported &
         &with is_smeared_ion_rep, because the derivatives of the smeared-ion &
         &correction to Ewald are not implemented.')

!    call utils_assert(pub_cell%nat_classical == 0,'Ewald forces are not &
!         &implemented for classical atoms (sparkles).')

    !First copy data from ONES structures into CASTEP cell structure
    !This is an inefficient hack, but it will be the most reliable way...
    if (pub_cell%nat_classical.ne.0) then
      call copy_to_castep_cell(current_cell,elements,classical_elements)
    else
      call copy_to_castep_cell(current_cell,elements)
    end if


    if (current_cell%num_ions<=0) then
       write (stderr,'(a)') 'Error in ewald_calculate_forces - &
            &current_cell uninitialised'
       call comms_abort
    end if

    !Set up parallelization over all active nodes
    num_active_nodes=pub_total_num_nodes


    !call trace_entry('ewald_calculate_forces',status)

    !Set up pretty unit labels
    !call io_unit_label(force_unit,force_label)
    !pub_on_root=(pub_my_node_id.eq.root_node_id)

    !setup packed arrays for charge and position for efficiency and parallelisation
    allocate(ion_charge(1:current_cell%num_ions),stat=ierr)
    call utils_alloc_check('ewald_calculate_forces','ion_charge',ierr)
    allocate(ion_position(1:3,1:current_cell%num_ions),stat=ierr)
    call utils_alloc_check('ewald_calculate_forces','ion_position',ierr)

    !check if cell has changed (or first pass)
    !lattice_changed=.false.
    !num_ions_changed=.false.
    !if (first_pass) then
       !old_num_ions=current_cell%num_ions
       !old_real_lattice=current_cell%real_lattice
       !allocate(old_ion_position(1:3,1:current_cell%num_ions),stat=ierr)
       !call utils_alloc_check('ewald_calculate_forces','old_ion_position',ierr)
!       allocate(old_ewald_forces(1:3,1:current_cell%max_ions_in_species,1:current_cell%num_species),stat=ierr)
       !allocate(old_ewald_forces(1:3,1:current_cell%num_ions),stat=ierr)
       !call utils_alloc_check('ewald_calculate_forces','old_ewald_forces',ierr)
       !old_ion_position=0.0_dp
       !old_ewald_forces=0.0_dp
       !first_pass=.false.
       !lattice_changed=.true.
       !ions_moved=.true.
       !num_ions_changed=.true.
    !else
       !do idim=1,3
       !   do jdim=1,3
       !      if (old_real_lattice(jdim,idim)/=current_cell%real_lattice(jdim,idim)) lattice_changed=.true.
       !   end do
       !end do
       !if(current_cell%num_ions/=old_num_ions) num_ions_changed=.true.
    !end if

    ! Need to synchronise as parallel distribution may have changed and old values differ on previously idle nodes

    !call comms_reduce_kp(lattice_changed,1,'OR')
    !call comms_reduce_bnd(lattice_changed,1,'OR')
    !call comms_reduce_gv(lattice_changed,1,'OR')
    !call comms_reduce('OR',lattice_changed)



    !calculate packed array of ion charges and positions, this code uses the castep cell as the
    !source of the data


    ion_i=1
    do ispec=1,current_cell%num_species
       do iatom=1,current_cell%num_ions_in_species(ispec)
          ion_charge(ion_i) =current_cell%ionic_charge(ispec) !*current_cell%mixture_weight(iatom,ispec)
          !ion_charge(ion_i) =current_cell%ionic_charge(ispec)

          !NB We do NOT do the compound atom trick here (unlike in energy & stress) as
          ! a) it is not necessary and
          ! b) it complicates things, as resulting force (whilst correct) is only non-zero on the
          !    compound atoms => need to re-split it at end for rest of CASTEP to re-combine it in MD etc!
          ! => simplest to leave this alone.
          do idim=1,3
             ion_position(idim,ion_i)=current_cell%ionic_positions(idim,iatom,ispec) &
               & -nint(current_cell%ionic_positions(idim,iatom,ispec))
          end do
          ion_i=ion_i+1
       end do
    end do

    !check to make sure that the ions have actually moved
    !ions_moved=.false.
    !num_ions_changed=.false.
    !if(current_cell%num_ions/=old_num_ions) then
    !   deallocate(old_ion_position,stat=ierr)
    !   call utils_dealloc_check('ewald_calculate_forces','old_ion_position',ierr)
    !   allocate(old_ion_position(1:3,1:current_cell%num_ions),stat=ierr)
    !   call utils_alloc_check('ewald_calculate_forces','old_ion_position',ierr)
    !   deallocate(old_ewald_forces,stat=ierr)
    !   call utils_dealloc_check('ewald_calculate_forces','old_ewald_forces',ierr)
    !   !allocate(old_ewald_forces(1:3,1:current_cell%max_ions_in_species,1:current_cell%num_species),stat=ierr)
    !   allocate(old_ewald_forces(1:3,1:current_cell%num_ions),stat=ierr)
    !   call utils_alloc_check('ewald_calculate_forces','old_ewald_forces',ierr)
    !   old_ion_position=0.0_dp
    !   num_ions_changed=.true.
    !else
    !   do ion_i=1,current_cell%num_ions
    !      do idim=1,3
    !         if (old_ion_position(idim,ion_i)/=ion_position(idim,ion_i)) ions_moved=.true.
    !      end do
    !   end do
    !end if

    ! Need to synchronise as parallel distribution may have changed and old values differ on previously idle nodes

    !call comms_reduce_kp(ions_moved,1,'OR')
    !call comms_reduce_bnd(ions_moved,1,'OR')
    !call comms_reduce_gv(ions_moved,1,'OR')
    !call comms_reduce('OR',ions_moved)

    !check if we actually need to do a calculation or not
    !if ((.not.lattice_changed).and.(.not.ions_moved).and.(.not.num_ions_changed)) then
    !   ewald_forces=old_ewald_forces
    !   if (pub_output_detail>NORMAL.and.pub_on_root) then
    !      write (stdout,*) 'NB: reusing old value of ewald_forces as nothing has changed!'
    !   end if
    !   deallocate(ion_charge,stat=ierr)
    !   call utils_dealloc_check('ewald_calculate_forces','ion_charge',ierr)
    !   deallocate(ion_position,stat=ierr)
    !   call utils_dealloc_check('ewald_calculate_forces','ion_position',ierr)
    !   !call trace_exit('ewald_calculate_forces',status)
    !   call castep_cell_dealloc(current_cell)
    !   call timer_clock("ewald_calculate_energy",2)
    !   return
    !end if

    !Calculate eta and number of cells needed for real and reciprocal space sums
    !if (lattice_changed.or.num_ions_changed)
    call ewald_calculate_num_cells(current_cell)

    !... so can now dimension allocatable arrays
    !if (lattice_changed.or.num_ions_changed) then
    !   if (allocated(recip_gauss_force)) then
    !      deallocate(recip_gauss_force,stat=ierr)
    !      call utils_dealloc_check('ewald_calculate_forces','recip_gauss_force',ierr)
    !   end if
    !end if

    ! ndmh: fixed insane memory allocation for temp forces!
    !allocate(tmp_ewald_force(1:3,1:current_cell%num_ions,1:current_cell%num_ions),stat=ierr)
    !call utils_alloc_check('ewald_calculate_forces','tmp_ewald_force',ierr)

    !and initialize to zero
    ewald_forces=0.0_dp
    !tmp_ewald_force=0.0_dp
    !tmp_force = 0.0_dp

    !miscellaneous scaling factors
    two_by_root_pi=2.0_dp/sqrt(pi)
    inv_two_eta_sq=1.0_dp/(2.0_dp*eta)**2
    two_i_pi=two_pi*cmplx_i
    scale_real_force =eta**3*elect
    scale_recip_force=(pi*elect)/(current_cell%volume*eta**2)


    !=========================================================================!
    !calculate real-space contribution to Ewald force
    !=========================================================================!
    call timer_clock("ewald_calculate_forces_realspace",1)
    !NB There is no contribution to the total real-space force, from ion(i) interacting with its own images
    !   and so there is no need to split this term off and calculate separately.

    !Now do the ion(i)->ion(j) interaction  (will fill in the ion(j)->ion(i) bit by symmetry at end)
    !parallelizing outer loop over all active nodes


    !kaw: Outer loop should only be over the QM atoms
    do ion_i=pub_my_node_id+1,current_cell%num_ions - pub_cell%nat_classical,num_active_nodes
       !do ion_j=ion_i+1,current_cell%num_ions
       do ion_j=1,current_cell%num_ions

          charge_ij=ion_charge(ion_i)*ion_charge(ion_j)

          !initialise the distances in the x,y,z directions between ion(i) and the periodic images of ion_j
          !in the corner cell of the block of cells included in the real-space part of the Ewald sum
          a_diff(1)=ion_position(1,ion_i)-ion_position(1,ion_j)+real(max_real_cells_x,kind=dp)
          a_diff(2)=ion_position(2,ion_i)-ion_position(2,ion_j)+real(max_real_cells_y,kind=dp)
          a_diff(3)=ion_position(3,ion_i)-ion_position(3,ion_j)+real(max_real_cells_z,kind=dp)

          rx0=dot_product(current_cell%real_lattice(:,1),a_diff(:))
          ry0=dot_product(current_cell%real_lattice(:,2),a_diff(:))
          rz0=dot_product(current_cell%real_lattice(:,3),a_diff(:))

          do nxyz=1,num_real_cells
             xdiff=rx0-real_dr(1,nxyz)
             ydiff=ry0-real_dr(2,nxyz)
             zdiff=rz0-real_dr(3,nxyz)

             if (abs(xdiff)<real_min_len) xdiff=0.0_dp !catch numerical underflow inaccurarices
             if (abs(ydiff)<real_min_len) ydiff=0.0_dp
             if (abs(zdiff)<real_min_len) zdiff=0.0_dp

             !make length dimensionless using eta
             real_len_sq=(xdiff**2+ydiff**2+zdiff**2)*(eta**2)
             !NB Might be tempted to remove this if-block on a vector machine, thinking that the cost of the extra
             !calculations would be offset by vectorizing this inner loop.
             !But, whilst this may be computationally more efficient, it would be PHYSICALLY WRONG!
             !We need a spherical cutoff to get reliable convergence of the sum. mijp2
             if (real_len_sq>real_min_len_sq.and.real_len_sq<real_cutoff_sq) then
                real_len=sqrt(real_len_sq)
                dF=charge_ij*(exp(-real_len_sq)*two_by_root_pi+ &
                     utils_erfc(real_len)/real_len)/real_len_sq*scale_real_force
                ewald_forces(1,ion_i) = ewald_forces(1,ion_i)+dF*xdiff
                ewald_forces(2,ion_i) = ewald_forces(2,ion_i)+dF*ydiff
                ewald_forces(3,ion_i) = ewald_forces(3,ion_i)+dF*zdiff
             end if
          end do

          !move onto next ion_j
       end do

       !move onto next ion_i
    end do

    !now fill in transpose by inversion symmetry
    !do ion_i=2,current_cell%num_ions
    !   do ion_j=1,ion_i-1
    !      do idim=1,3
    !         tmp_ewald_force(idim,ion_j,ion_i)=-tmp_ewald_force(idim,ion_i,ion_j)
    !      end do
    !   end do
    !end do

    !and gather up over all the nodes and convert to atomic units
    !call comms_reduce_gv(tmp_ewald_force(1,1,1),3*current_cell%num_ions**2,'SUM')
    !call comms_reduce_bnd(tmp_ewald_force(1,1,1),3*current_cell%num_ions**2,'SUM')
    !call comms_reduce_kp(tmp_ewald_force(1,1,1),3*current_cell%num_ions**2,'SUM')
    !call comms_reduce('SUM',tmp_ewald_force)
    !tmp_ewald_force=tmp_ewald_force*scale_real_force

    call timer_clock("ewald_calculate_forces_realspace",2)
    !print real-space force (user units)
!    if (pub_output_detail>NORMAL.and.pub_on_root) then
!       write (stdout,99) 'Ewald: real-space force:'
!       do ion_i=1,current_cell%num_ions
!          do ion_j=1,current_cell%num_ions
!             write (stdout,97) ion_j,ion_i,(tmp_ewald_force(idim,ion_j,ion_i),idim=1,3)
!             write (stdout,97) ion_j,ion_i,(io_atomic_to_unit(tmp_ewald_force(idim,ion_j,ion_i),force_unit), &
!                  &  idim=1,3),trim(force_label)
!          end do
!       end do
!    end if
    !... finished real space term

    !... store the real space contribution in packed array ...
    ion_k=0
    !do ispec=1,current_cell%num_species
    !   do ion_i=1,current_cell%num_ions_in_species(ispec)
    !      ion_k=ion_k+1
    !      do ion_j=1,current_cell%num_ions
    !         do idim=1,3
    !            ewald_forces(idim,ion_i,ispec)=ewald_forces(idim,ion_i,ispec)+ &
    !                 & tmp_ewald_force(idim,ion_j,ion_k)
    !         end do
    !      end do
    !   end do
    !end do
    !do ion_i=1,current_cell%num_ions
    !   do ion_j=1,current_cell%num_ions
    !      do idim=1,3
    !         ewald_forces(idim,ion_i)=ewald_forces(idim,ion_i)+ &
    !              & tmp_ewald_force(idim,ion_i)
    !      end do
    !   end do
    !end do


    !=========================================================================!
    !calculate reciprocal-space contribution to Ewald force
    !=========================================================================!
    call timer_clock("ewald_calculate_forces_recipspace",1)

    allocate(exp_recip_phase_x(1:num_recip_cells_x),stat=ierr)
    call utils_alloc_check('ewald_calculate_forces','exp_recip_phase_x',ierr)
    allocate(exp_recip_phase_y(1:num_recip_cells_y),stat=ierr)
    call utils_alloc_check('ewald_calculate_forces','exp_recip_phase_y',ierr)
    allocate(exp_recip_phase_z(1:num_recip_cells_z),stat=ierr)
    call utils_alloc_check('ewald_calculate_forces','exp_recip_phase_z',ierr)
    allocate(recip_gauss_force(1:3,1:num_recip_cells),stat=ierr)
    call utils_alloc_check('ewald_calculate_forces','recip_gauss_force',ierr)
    exp_recip_phase_x=cmplx_0
    exp_recip_phase_y=cmplx_0
    exp_recip_phase_z=cmplx_0
    recip_gauss_force=0.0_dp
    !if (lattice_changed.or.num_ions_changed) then
    !   recip_gauss_force=0.0_dp
    !end if

    !reset tmp array
    !tmp_ewald_force=0.0_dp

    !First we precompute various factors ...
    !if (lattice_changed.or.num_ions_changed) then

       a_diff(1)=real(max_recip_cells_x,kind=dp)
       a_diff(2)=real(max_recip_cells_y,kind=dp)
       a_diff(3)=real(max_recip_cells_z,kind=dp)

       gx0=-dot_product(current_cell%recip_lattice(:,1),a_diff(:))
       gy0=-dot_product(current_cell%recip_lattice(:,2),a_diff(:))
       gz0=-dot_product(current_cell%recip_lattice(:,3),a_diff(:))

       !parallelizing outer loop over all active nodes
       !NB This way of doing things with the single loop and the lookup table is slightly faster than the previous triple loop
       !   AND, more significantly, it makes a longer loop => better load balancing with more processors. mijp.
       do ng=pub_my_node_id+1,num_recip_cells,num_active_nodes
          gx=gx0+recip_dg(1,ng)
          gy=gy0+recip_dg(2,ng)
          gz=gz0+recip_dg(3,ng)

          !calculate the dimensionless length of the reciprocal lattice vector
          recip_len_sq=(gx**2+gy**2+gz**2)*inv_two_eta_sq
          !NB Might be tempted to remove this if-block on a vector machine, thinking that the cost of the extra
          !calculations would be offset by vectorizing this inner loop.
          !But, whilst this may be computationally more efficient, it would be PHYSICALLY WRONG!
          !We need a spherical cutoff to get reliable convergence of the sum. mijp2
          if(recip_len_sq>recip_min_len_sq.and.recip_len_sq<recip_cutoff_sq) then
             dF=exp(-recip_len_sq)/recip_len_sq
             recip_gauss_force(1,ng)=gx*dF
             recip_gauss_force(2,ng)=gy*dF
             recip_gauss_force(3,ng)=gz*dF
          else
             recip_gauss_force(:,ng)=0.0_dp
          endif

       end do

       !and we need recip_gauss on all nodes, so gather it up ...
       !call comms_reduce_gv(recip_gauss_force(1,1),3*num_recip_cells,'SUM')
       !call comms_reduce_bnd(recip_gauss_force(1,1),3*num_recip_cells,'SUM')
       !call comms_reduce_kp(recip_gauss_force(1,1),3*num_recip_cells,'SUM')
       call comms_reduce('SUM',recip_gauss_force)

    !end if !lattice_changed.or.num_ions_changed


    !Now do the ion(i)->ion(j) interaction  (will fill in the ion(j)->ion(i) bit by symmetry at end)
    !kaw: See above
    do ion_i=pub_my_node_id+1,current_cell%num_ions - pub_cell%nat_classical,num_active_nodes
       !do ion_j=ion_i+1,current_cell%num_ions
       do ion_j=1,current_cell%num_ions

          charge_ij=ion_charge(ion_i)*ion_charge(ion_j)

          !precompute x-direction phase factors=exp(i.Gx.(r1x-r2x))
          recip_phase=(ion_position(1,ion_i)-ion_position(1,ion_j))*two_i_pi
          exp_recip_phase=exp(recip_phase)
          exp_recip_phase_x(1)=exp(-recip_phase*max_recip_cells_x)
          do nx=2,num_recip_cells_x
             exp_recip_phase_x(nx)=exp_recip_phase_x(nx-1)*exp_recip_phase
          end do

          !precompute y-direction phase factors
          recip_phase=(ion_position(2,ion_i)-ion_position(2,ion_j))*two_i_pi
          exp_recip_phase=exp(recip_phase)
          exp_recip_phase_y(1)=exp(-recip_phase*max_recip_cells_y)
          do ny=2,num_recip_cells_y
             exp_recip_phase_y(ny)=exp_recip_phase_y(ny-1)*exp_recip_phase
          end do

          !precompute z-direction phase factors
          recip_phase=(ion_position(3,ion_i)-ion_position(3,ion_j))*two_i_pi
          exp_recip_phase=exp(recip_phase)
          exp_recip_phase_z(1)=exp(-recip_phase*max_recip_cells_z)
          do nz=2,num_recip_cells_z
             exp_recip_phase_z(nz)=exp_recip_phase_z(nz-1)*exp_recip_phase
          end do

          !speedup
          !tmp_force(:)=tmp_ewald_force(:,ion_i)

          do ng=1,num_recip_cells
             nx=recip_nxyz(1,ng)
             ny=recip_nxyz(2,ng)
             nz=recip_nxyz(3,ng)

             exp_recip_phase=exp_recip_phase_x(nx)*exp_recip_phase_y(ny)*exp_recip_phase_z(nz)
             dF=aimag(exp_recip_phase)*charge_ij*scale_recip_force
             !do idim=1,3
             !   tmp_ewald_force(idim,ion_i)=tmp_ewald_force(idim,ion_i)+recip_gauss_force(idim,ng)*dF
             !end do
             ewald_forces(:,ion_i)=ewald_forces(:,ion_i)+recip_gauss_force(:,ng)*dF

          end do

          !speedup
          !tmp_ewald_force(:,ion_i)=tmp_force(:)

          !move onto the next ion_j
       end do

       !move onto the next ion_i
    end do

    !now fill in transpose by inversion symmetry
    !do ion_i=2,current_cell%num_ions
    !   do ion_j=1,ion_i-1
    !      do idim=1,3
    !         tmp_ewald_force(idim,ion_j,ion_i)=-tmp_ewald_force(idim,ion_i,ion_j)
    !      end do
    !   end do
    !end do

    !and gather up over all the nodes and convert to atomic units
    !call comms_reduce_gv(tmp_ewald_force(1,1,1),3*current_cell%num_ions**2,'SUM')
    !call comms_reduce_bnd(tmp_ewald_force(1,1,1),3*current_cell%num_ions**2,'SUM')
    !call comms_reduce_kp(tmp_ewald_force(1,1,1),3*current_cell%num_ions**2,'SUM')
    !call comms_reduce('SUM',tmp_ewald_force)
    !tmp_ewald_force=tmp_ewald_force*scale_recip_force
    call comms_reduce('SUM',ewald_forces)
    call timer_clock("ewald_calculate_forces_recipspace",2)
    !print recip-space force (user units)
!    if (pub_output_detail>NORMAL.and.pub_on_root) then
!       write (stdout,99) 'Ewald: recip-space force:'
!       do ion_i=1,current_cell%num_ions
!          do ion_j=1,current_cell%num_ions
!             write (stdout,97) ion_j,ion_i,(tmp_ewald_force(idim,ion_j,ion_i),idim=1,3)
!             write (stdout,97) ion_j,ion_i,(io_atomic_to_unit(tmp_ewald_force(idim,ion_j,ion_i),force_unit), &
!                  & idim=1,3), trim(force_label)
!          end do
!       end do
!    end if
    !... finished reciprocal space term

    !... finally add reciprocal space contribution (already have real-space in packed array) ...
    !ion_k=0
    !do ispec=1,current_cell%num_species
    !   do ion_i=1,current_cell%num_ions_in_species(ispec)
    !      ion_k=ion_k+1
    !      do ion_j=1,current_cell%num_ions
    !         do idim=1,3
    !            ewald_forces(idim,ion_i,ispec)=ewald_forces(idim,ion_i,ispec)+ &
    !                 & tmp_ewald_force(idim,ion_j,ion_k)
    !         end do
    !      end do
    !   end do
    !end do
    !do ion_i=1,current_cell%num_ions
    !   do ion_j=1,current_cell%num_ions
    !      do idim=1,3
    !         ewald_forces(idim,ion_i)=ewald_forces(idim,ion_i)+ &
    !              & tmp_ewald_force(idim,ion_i)
    !      end do
    !   end do
    !end do

!    if (pub_output_detail>NORMAL.and.pub_on_root) then
!       write (stdout,99) 'Ewald force:'
!       do ispec=1,current_cell%num_species
!          do ion_i=1,current_cell%num_ions
!             write (stdout,96) ispec,ion_i,(ewald_forces(idim,ion_i),idim=1,3)
!             write (stdout,96) ispec,ion_i,(io_atomic_to_unit(ewald_forces(idim,ion_i,ispec),force_unit),idim=1,3), &
!                  & trim(force_label)
!          end do
!       end do
!    end if

    !store result for next time ...
    !old_ewald_forces=ewald_forces
    !old_num_ions=current_cell%num_ions
    !old_ion_position=ion_position
    !old_real_lattice=current_cell%real_lattice

    !... clean up temporary arrays ...
    !deallocate(tmp_ewald_force,stat=ierr)
    !call utils_dealloc_check('ewald_calculate_forces','tmp_ewald_force',ierr)
    deallocate(recip_gauss_force,stat=ierr)
    call utils_dealloc_check('ewald_calculate_forces','recip_gauss_force',ierr)
    deallocate(exp_recip_phase_x,stat=ierr)
    call utils_dealloc_check('ewald_calculate_forces','exp_recip_phase_x',ierr)
    deallocate(exp_recip_phase_y,stat=ierr)
    call utils_dealloc_check('ewald_calculate_forces','exp_recip_phase_y',ierr)
    deallocate(exp_recip_phase_z,stat=ierr)
    call utils_dealloc_check('ewald_calculate_forces','exp_recip_phase_z',ierr)

    deallocate(ion_charge,stat=ierr)
    call utils_dealloc_check('ewald_calculate_forces','ion_charge',ierr)
    deallocate(ion_position,stat=ierr)
    call utils_dealloc_check('ewald_calculate_forces','ion_position',ierr)

    call castep_cell_dealloc(current_cell)

    ! Free up comms memory
    call comms_free

!99  format(A40)
!97  format('   Atom #',I3,' Atom #',I3,' Fx:',G12.6,' Fy:',G12.6,' Fz:',G12.6,1X,A)
!96  format('Species #',I3,' Atom #',I3,' Fx:',G12.6,' Fy:',G12.6,' Fz:',G12.6,1X,A)

    !add in 3D correction?
    !if (EW3DC) call ewald_calculate_EW3DC_forces(ewald_forces)

    !call trace_exit('ewald_calculate_forces',status)

    call timer_clock("ewald_calculate_forces",2)

    return
  end subroutine ewald_calculate_forces

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ewald_exit(elements, ewald_forces)

    !=========================================================================!
    ! Clean up by deallocating saved arrays in the ewald module.              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   ewald_forces, intent=out, just for compatibility with the subroutines !
    !       called.                                                           !
    !-------------------------------------------------------------------------!
    ! Written by Quintin Hill 18/03/2008                                      !
    !=========================================================================!

    use simulation_cell, only: pub_cell
    use ion, only: element

    implicit none

    type(ELEMENT), dimension(1:pub_cell%nat), intent(in) :: elements
    real (kind=dp), dimension(1:3,1:pub_cell%nat), intent(out) :: ewald_forces
    real (kind=dp) :: ewald_energy
    logical, parameter :: finished = .true.

    call ewald_calculate_energy(elements, ewald_energy, finished)
    call ewald_calculate_forces(elements, ewald_forces, finished)

  end subroutine ewald_exit

#if 0
  subroutine ewald_calc_storage(ram,disc,type)
    !=========================================================================!
    ! Estimates memory required for the various ewald_calculate_* routines.   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   ram,   intent=out, the additional RAM required for specified routine  !
    !   disc,  intent=out, the additional disc space required                 !
    !   type, optional intent=in, the type of Ewald calculation to report     !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   cell for number of atoms etc                                          !
    !   io for output routines                                                !
    !   algor for sizeof routines                                             !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   cell_read must have already been called                               !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v1.0, 25/08/2008                               !
    !=========================================================================!
    use algor,      only: algor_sizeof
    use parameters, only: calculate_stress

    use trace, only: trace_entry, trace_exit

    implicit none

    real(kind=dp),     intent(out)         :: ram
    real(kind=dp),     intent(out)         :: disc
    character(len=*), optional, intent(in) :: type

    !local variables
    character(len=5) :: local_type
    !integer          :: status
    integer          :: num_ions, num_spec, max_ions

    !call trace_entry('ewald_calc_storage',status)

    if (present(type)) then
       local_type=type                !value passed in
    else
       if (calculate_stress) then
          local_type="EFS"            !Can be E, F, S, 2 or combinations
       else
          local_type="EF"             !Can be E, F, S, 2 or combinations
       end if
    endif

    ram  = 0.0_dp
    disc = 0.0_dp
    num_ions=current_cell%num_ions
    num_spec=current_cell%num_species
    max_ions=current_cell%max_ions_in_species

    !make sure Ewald is properly initialised
    if (num_recip_cells<1) call ewald_calculate_num_cells

    !check allowed values for type
    if (index(local_type,'E')>0) then
       !add contribution due to ewald_calculate_energy
       ram = ram + real(num_ions,dp)* real(4*algor_sizeof(1.0_dp),dp)         !ion_charge & ion_position
       ram = ram + real(num_ions,dp)* real(3*algor_sizeof(1.0_dp),dp)         !old_ion_position  *saved*
       ram = ram + real(num_recip_cells,dp)* real(algor_sizeof(1.0_dp),dp)    !recip_gauss       *saved*
       ram = ram + real(num_recip_cells_x,dp)* real(algor_sizeof(cmplx_i),dp) !exp_recip_phase_x
       ram = ram + real(num_recip_cells_y,dp)* real(algor_sizeof(cmplx_i),dp) !exp_recip_phase_y
       ram = ram + real(num_recip_cells_z,dp)* real(algor_sizeof(cmplx_i),dp) !exp_recip_phase_z
    end if
    if (index(local_type,'F')>0) then
       !add contribution due to ewald_calculate_force
       ram = ram + real(num_ions,dp)* real(4*algor_sizeof(1.0_dp),dp)         !ion_charge & ion_position
       ram = ram + real(num_ions,dp)* real(3*algor_sizeof(1.0_dp),dp)         !old_ion_position *saved*
       ram = ram + real(max_ions,dp)* real(3*num_spec*algor_sizeof(1.0_dp),dp)!old_ewald_forces *saved*
       ram = ram + real(num_recip_cells,dp)* real(3*algor_sizeof(1.0_dp),dp)  !recip_gauss_force *saved*
       ram = ram + real(num_ions,dp)*real(num_ions,dp)* real(3*algor_sizeof(1.0_dp),dp) !tmp_ewald_force
       ram = ram + real(num_recip_cells_x,dp)* real(algor_sizeof(cmplx_i),dp) !exp_recip_phase_x
       ram = ram + real(num_recip_cells_y,dp)* real(algor_sizeof(cmplx_i),dp) !exp_recip_phase_y
       ram = ram + real(num_recip_cells_z,dp)* real(algor_sizeof(cmplx_i),dp) !exp_recip_phase_z
    end if
    if (index(local_type,'S')>0) then
       !add contribution due to ewald_calculate_stress
       ram = ram + real(num_ions,dp)* real(4*algor_sizeof(1.0_dp),dp)         !ion_charge & ion_position
       ram = ram + real(num_ions,dp)* real(3*algor_sizeof(1.0_dp),dp)         !old_ion_position *saved*
       ram = ram + real(6*algor_sizeof(1.0_dp),dp)                            !old_ewald_stress *saved*
       ram = ram + real(num_recip_cells,dp)* real(6*algor_sizeof(1.0_dp),dp)  !recip_gauss_stress *saved*
       ram = ram + real(num_recip_cells_x,dp)* real(algor_sizeof(cmplx_i),dp) !exp_recip_phase_x
       ram = ram + real(num_recip_cells_y,dp)* real(algor_sizeof(cmplx_i),dp) !exp_recip_phase_y
       ram = ram + real(num_recip_cells_z,dp)* real(algor_sizeof(cmplx_i),dp) !exp_recip_phase_z
    end if
    if (index(local_type,'2')>0) then
       !add contribution due to ewald_calculate_E2
       ram = ram + real(num_ions,dp)*real(num_ions,dp)* real(6*algor_sizeof(1.0_dp),dp)  !real_diag_E2 *saved*
       ram = ram + real(num_ions,dp)*real(num_ions,dp)* real(6*algor_sizeof(1.0_dp),dp)  !recip_diag_E2 *saved*
       ram = ram + real(num_ions,dp)*real(num_ions,dp)* real(6*algor_sizeof(cmplx_i),dp) !real_E2
       ram = ram + real(num_ions,dp)*real(num_ions,dp)* real(6*algor_sizeof(cmplx_i),dp) !recip_E2
       ram = ram + real(num_real_cells,dp)* real(algor_sizeof(cmplx_i),dp)    !exp_real_phase *saved*
       ram = ram + real(num_ions,dp)* real(3*algor_sizeof(1.0_dp),dp)         !old_ion_position *saved*
       ram = ram + real(num_ions,dp)* real(3*algor_sizeof(1.0_dp),dp)         !ion_position
    end if

    ! Convert to megabytes
    ram  = ram/(1024.0_dp**2)
    disc = disc/(1024.0_dp**2)

    !call trace_exit('ewald_calc_storage',status)

    return
  end subroutine ewald_calc_storage
#endif
  !---------------------------------------------------------------------------!
  !                P R I V A T E   R O U T I N E S                            !
  !---------------------------------------------------------------------------!
  subroutine ewald_calculate_num_cells(current_cell)
    !=========================================================================!
    ! Calculate the optimum value of eta for the current cell, and the        !
    ! optimum number of cells for the real and reciprocal space sums.         !
    ! Also setup the lookup tables for moving between the different cells.    !
    ! Depends only on size of current_cell.                                   !
    ! Also reads devel_code to set various private options, such as reverting !
    ! to old eta algorithm or setting EW3DC on/off.                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! If the current_cell lattice vectors change, then:                       !
    !   eta (convergence parmameter) is recalculated                          !
    !   max_real_cells_x,y,z          "      "                                !
    !   num_real_cells_x,y,z          "      "                                !
    !   num_real_cells                "      "                                !
    !   max_recip_cells_x,y,z         "      "                                !
    !   num_recip_cells_x,y,z         "      "                                !
    !   num_recip_cells               "      "                                !
    !   real_dr                       "      "                                !
    !   recip_dg                      "      "                                !
    !   recip_nxyz                    "      "                                !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   cell for current_cell lattice vectors (real and reciprocal space)     !
    !   io   for abort                                                        !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   a1=length of first real-space lattice vector                          !
    !   b1=   '   "   "    recip-space   "      "                             !
    !   etc.                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   current_cell must be valid                                            !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.1, 05/07/2000                               !
    !=========================================================================!
    use comms, only: comms_abort, comms_bcast, pub_on_root, pub_root_node_id
    use constants, only: DP, stderr, stdout, PI, NORMAL
    use rundat, only: pub_devel_code, pub_output_detail
    use simulation_cell, only: castep_cell
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    type(castep_cell), intent(in) :: current_cell

    !local variables
    real(kind=dp) :: min_length
    real(kind=dp) :: a1,a2,a3,a
    real(kind=dp) :: b1,b2,b3,b

    integer       :: nx,ny,nz,nxyz,ng
    integer       :: ierr

    !devel code stuff
    character (len=80) :: ewald_devel_code
    integer            :: ewald_start_pos, ewald_stop_pos, ewald_test_pos

    !eta calculation
    real(kind=dp), save  :: precision = 36.0_dp
    !convergence => neglecting real space term ~exp(-precision): 16->7dp, 25->11dp, 36->16dp, 49->21dp

    !NB Old algorithm (pre v1.50) had hard-wired cutoffs and hence implicit precision.
    !Having cut-off=4.0 is equivalent to precision = 16, v1.37 cutoff -> 10, v1.40 cutoff -> 20 => precision = 400
    !With v1.50 introduced this more flexible scheme which enables control of precision and
    !does away with hard-wired cutoffs.

    logical, save        :: old_eta=.false.      !T => use algorithm pre v1.50
    real (kind=dp), save :: real_to_recip=1.0_dp !approx real/recip calc time

    real (kind=dp)       :: mean_a, real_skewness

    !integer              :: status

    !call trace_entry('ewald_calculate_num_cells',status)

    !Check local copy of devel_code string and act accordingly
    if (pub_on_root) then
       ewald_devel_code=pub_devel_code
       if (len_trim(ewald_devel_code)>0) then
          ewald_start_pos=index(ewald_devel_code,'EWALD:')
          ewald_stop_pos=index(ewald_devel_code,':ENDEWALD')
          if (ewald_stop_pos<=0) ewald_stop_pos=len_trim(ewald_devel_code) !missing end so go to end of string
          if (ewald_start_pos>0) then

             !found an EWALD section so parse what follows

             !look for a setting of "precision"
             ewald_test_pos=index(ewald_devel_code,'PREC=')
             if (ewald_test_pos>ewald_start_pos.and.ewald_test_pos<ewald_stop_pos) then
                ewald_test_pos=ewald_test_pos+len('PREC=')
                !read next character string to : as the convergence window width
                read(ewald_devel_code(ewald_test_pos:ewald_test_pos+ &
                     & index(ewald_devel_code(ewald_test_pos:ewald_stop_pos),':')-2),*) precision
                if (pub_on_root.and.(pub_output_detail>=NORMAL)) then
                   write(stdout,'(a40,f15.6)') 'Ewald: Precision = ', precision
                end if
             end if

             !look for setting of old or new algorithm for eta
             if (index(ewald_devel_code(ewald_start_pos:ewald_stop_pos),'OLD=T')>0) then
                old_eta=.true.
             else if (index(ewald_devel_code(ewald_start_pos:ewald_stop_pos),'OLD=F')>0) then
                old_eta=.false.
             end if

             !look for setting of real_to_recip ratio
             ewald_test_pos=index(ewald_devel_code,'R2R=')
             if (ewald_test_pos>ewald_start_pos.and.ewald_test_pos<ewald_stop_pos) then
                ewald_test_pos=ewald_test_pos+len('R2R=')
                !read next character string to : as the convergence window width
                read(ewald_devel_code(ewald_test_pos:ewald_test_pos+ &
                     & index(ewald_devel_code(ewald_test_pos:ewald_stop_pos),':')-2),*) real_to_recip
             end if

             !look for turning on Ewald3D correction
             !if (index(ewald_devel_code(ewald_start_pos:ewald_stop_pos),'EW3DC=T')>0) then
             !   EW3DC=.true.
             !else if (index(ewald_devel_code(ewald_start_pos:ewald_stop_pos),'EW3DC=F')>0 ) then
             !   EW3DC=.false.
             !end if

         end if
       end if
    end if
    !and make sure all nodes in sync
    call comms_bcast(pub_root_node_id,precision)
    call comms_bcast(pub_root_node_id,old_eta)
    call comms_bcast(pub_root_node_id,real_to_recip)
    !call comms_bcast(pub_root_node_id,EW3DC)

    !get some basics about the cell
    a1=sqrt(current_cell%real_lattice(1,1)**2+current_cell%real_lattice(1,2)**2+current_cell%real_lattice(1,3)**2)
    a2=sqrt(current_cell%real_lattice(2,1)**2+current_cell%real_lattice(2,2)**2+current_cell%real_lattice(2,3)**2)
    a3=sqrt(current_cell%real_lattice(3,1)**2+current_cell%real_lattice(3,2)**2+current_cell%real_lattice(3,3)**2)
    mean_a=(a1*a2*a3)**(1.0_dp/3.0_dp)
    real_skewness=max(a1,a2,a3)/mean_a

    if (old_eta) then
       !deprecated original method with hard-wired lengths - kept purely for reversion testing
       !min_length=(a1*a2*a3)**(1.0_dp/3.0_dp)
       min_length=min(a1,a2,a3)
       eta=sqrt(pi)/min_length
       real_cutoff = 20.0_dp
       recip_cutoff = 20.0_dp
    else
       !new alternative approach - eta now depends on number of atoms and cutoff goes with precision
       !NB At a decent level of precision, results & timings pretty robust wrt real_to_recip (timing ratio)
       !NB A geometric mean based algorithm for eta should be unbiased by skew and most reliable & robust
       eta=(sqrt(pi)/mean_a)*(real_to_recip*current_cell%num_ions)**(1.0_dp/6.0_dp)
       !NB A min(length) [or max()] based algorithm for eta will biase real:recip cutoffs and will therefore
       !give underconvergence for very skew cells.
       !eta=(sqrt(pi)/min(a1,a2,a3))*(real_to_recip*current_cell%num_ions)**(1.0_dp/6.0_dp)
       !NB A volume based algorithm for eta will not discriminate between two cells of same volume
       !where one has been skewed wrt other. This sounds good, as final answer should be invariant to skewing,
       !but need to put more sums into shorter real-space length when skewed so actually this is bad.
       !eta=sqrt(pi)*(real_to_recip*current_cell%num_ions/(current_cell%volume**2))**(1.0_dp/6.0_dp)
       real_cutoff = sqrt(precision)*real_skewness
       recip_cutoff = 2.0_dp*sqrt(precision)*real_skewness
    end if
    call comms_bcast(pub_root_node_id,eta)
    call comms_bcast(pub_root_node_id,real_cutoff)
    call comms_bcast(pub_root_node_id,recip_cutoff)

    !set the other bits to be in sync
    real_cutoff_sq  = real_cutoff**2
    real_min_len_sq = real_min_len**2
    recip_cutoff_sq = recip_cutoff**2
    recip_min_len_sq= recip_min_len**2

    !calculate new number of cells in real space ...
    a=sqrt((real_cutoff/eta)**2 &
         & + (abs(current_cell%real_lattice(1,1))+abs(current_cell%real_lattice(2,1))+abs(current_cell%real_lattice(3,1)))**2 &
         & + (abs(current_cell%real_lattice(1,2))+abs(current_cell%real_lattice(2,2))+abs(current_cell%real_lattice(3,2)))**2 &
         & + (abs(current_cell%real_lattice(1,3))+abs(current_cell%real_lattice(2,3))+abs(current_cell%real_lattice(3,3)))**2)

    max_real_cells_x=int(a/a1)+1
    max_real_cells_y=int(a/a2)+1
    max_real_cells_z=int(a/a3)+1

    !paranoid sync
    call comms_bcast(pub_root_node_id,max_real_cells_x)
    call comms_bcast(pub_root_node_id,max_real_cells_y)
    call comms_bcast(pub_root_node_id,max_real_cells_z)

    num_real_cells_x=2*max_real_cells_x+1
    num_real_cells_y=2*max_real_cells_y+1
    num_real_cells_z=2*max_real_cells_z+1

    num_real_cells=num_real_cells_x*num_real_cells_y*num_real_cells_z

    !calculate new number of cells in reciprocal space ...
    b1=sqrt(current_cell%recip_lattice(1,1)**2+current_cell%recip_lattice(1,2)**2+current_cell%recip_lattice(1,3)**2)
    b2=sqrt(current_cell%recip_lattice(2,1)**2+current_cell%recip_lattice(2,2)**2+current_cell%recip_lattice(2,3)**2)
    b3=sqrt(current_cell%recip_lattice(3,1)**2+current_cell%recip_lattice(3,2)**2+current_cell%recip_lattice(3,3)**2)

    b=sqrt((recip_cutoff*2.0_dp*eta)**2 &               !factors of pi cancel
         & +(abs(current_cell%recip_lattice(1,1))+abs(current_cell%recip_lattice(2,1))+abs(current_cell%recip_lattice(3,1)))**2 &
         & +(abs(current_cell%recip_lattice(1,2))+abs(current_cell%recip_lattice(2,2))+abs(current_cell%recip_lattice(3,2)))**2 &
         & +(abs(current_cell%recip_lattice(1,3))+abs(current_cell%recip_lattice(2,3))+abs(current_cell%recip_lattice(3,3)))**2)

    max_recip_cells_x=int(b/b1)+1
    max_recip_cells_y=int(b/b2)+1
    max_recip_cells_z=int(b/b3)+1

    !paranoid sync
    call comms_bcast(pub_root_node_id,max_recip_cells_x)
    call comms_bcast(pub_root_node_id,max_recip_cells_y)
    call comms_bcast(pub_root_node_id,max_recip_cells_z)

    num_recip_cells_x=2*max_recip_cells_x+1
    num_recip_cells_y=2*max_recip_cells_y+1
    num_recip_cells_z=2*max_recip_cells_z+1

    num_recip_cells=num_recip_cells_x*num_recip_cells_y*num_recip_cells_z

    !setup lookup tables
    if (allocated(real_dr)) then
       deallocate(real_dr,stat=ierr)
       call utils_dealloc_check('ewald_calculate_num_cells','real_dr',ierr)
    end if
    allocate(real_dr(1:3,1:num_real_cells),stat=ierr)
    call utils_alloc_check('ewald_calculate_num_cells','real_dr',ierr)

    nxyz=0
    do nz=1,num_real_cells_z
       do ny=1,num_real_cells_y
          do nx=1,num_real_cells_x
             nxyz=nxyz+1
             real_dr(1,nxyz)=(nx-1)*current_cell%real_lattice(1,1)+ &
                  & (ny-1)*current_cell%real_lattice(2,1)+(nz-1)*current_cell%real_lattice(3,1)
             real_dr(2,nxyz)=(nx-1)*current_cell%real_lattice(1,2)+ &
                  & (ny-1)*current_cell%real_lattice(2,2)+(nz-1)*current_cell%real_lattice(3,2)
             real_dr(3,nxyz)=(nx-1)*current_cell%real_lattice(1,3)+ &
                  & (ny-1)*current_cell%real_lattice(2,3)+(nz-1)*current_cell%real_lattice(3,3)
          end do
       end do
    end do

    if (allocated(recip_dg)) then
       deallocate(recip_dg,stat=ierr)
       call utils_dealloc_check('ewald_calculate_num_cells','recip_dg',ierr)
    end if
    allocate(recip_dg(1:3,1:num_recip_cells),stat=ierr)
    call utils_alloc_check('ewald_calculate_num_cells','recip_dg',ierr)

    if (allocated(recip_nxyz)) then
       deallocate(recip_nxyz,stat=ierr)
       call utils_dealloc_check('ewald_calculate_num_cells','recip_nxyz',ierr)
    end if
    allocate(recip_nxyz(1:3,1:num_recip_cells),stat=ierr)
    call utils_alloc_check('ewald_calculate_num_cells','recip_nxyz',ierr)

    ng=0
    do nz=1,num_recip_cells_z
       do ny=1,num_recip_cells_y
          do nx=1,num_recip_cells_x
             ng=ng+1
             recip_nxyz(1,ng)=nx
             recip_nxyz(2,ng)=ny
             recip_nxyz(3,ng)=nz
             recip_dg(1,ng)=(nx-1)*current_cell%recip_lattice(1,1)+ &
                  & (ny-1)*current_cell%recip_lattice(2,1)+(nz-1)*current_cell%recip_lattice(3,1)
             recip_dg(2,ng)=(nx-1)*current_cell%recip_lattice(1,2)+ &
                  & (ny-1)*current_cell%recip_lattice(2,2)+(nz-1)*current_cell%recip_lattice(3,2)
             recip_dg(3,ng)=(nx-1)*current_cell%recip_lattice(1,3)+ &
                  & (ny-1)*current_cell%recip_lattice(2,3)+(nz-1)*current_cell%recip_lattice(3,3)
          end do
       end do
    end do

    !call trace_exit('ewald_calculate_num_cells',status)

    return
  end subroutine ewald_calculate_num_cells

end module ewald
