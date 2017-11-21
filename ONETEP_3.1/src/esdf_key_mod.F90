! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!
! Module to hold keyword list. This must be updated as
! new keywords are brought into existence.
!
! The 'LABEL' is the LABEL as used in calling the esdf routines
! 'TYP' defines the TYPe, with the following syntax. It is 3 characters
! long.
! The first indicates:
!  I - integer
!  S - single
!  D - double
!  P - physical
!  T - string (text)
!  E - defined (exists)
!  L - boolean (logical)
!  B - block
! The second is always a colon (:)
! The third indicates the "level" of the keyword
!  B - Basic
!  I - Intermediate
!  E - Expert
!  D - Dummy
!
! 'DSCRPT' is a description of the variable. It should contain a (short) title
! enclosed between *! ... !*, and then a more detailed description of the
! variable.
!
module esdf_key

  implicit none

  type KW_TYPE
     Character(30)   :: LABEL
     Character(3)    :: TYP
!    Character(2500) :: DSCRPT
!CW : 3000->300, for 32-bit portability
     Character(300)  :: DSCRPT
  end type KW_TYPE
!CW

!CW
#define i 352
  integer, parameter :: numkw=423
! integer, parameter :: numkw=352
!END CW

  type(KW_TYPE) :: KW(numkw)

  ! Now define the keywords
  DATA KW(1)%LABEL /"LATTICE_CART"/
  DATA KW(1)%TYP /"B:B"/
  DATA KW(1)%DSCRPT/"*! The simulation cell lattice vectors !*"/

  DATA KW(2)%LABEL /"PSINC_SPACING"/
  DATA KW(2)%TYP /"T:E"/
  DATA KW(2)%DSCRPT /"*! PSINC grid spacing in atomic units !*"/

  DATA KW(3)%LABEL /"GEOMETRY"/
  DATA KW(3)%TYP /"B:B"/
  DATA KW(3)%DSCRPT/"*! input geometry !*"/

  DATA KW(4)%LABEL /"ATOMIC_SETS"/
  DATA KW(4)%TYP /"B:I"/
  DATA KW(4)%DSCRPT /"*! atomic set names for each atom !*"/

  DATA KW(5)%LABEL /"PPD_NPOINTS"/
  DATA KW(5)%TYP /"T:E"/
  DATA KW(5)%DSCRPT /"*! PPD edge length in grid points for &
       &each lattice direction !*"/

  DATA KW(6)%LABEL /"XC_FUNCTIONAL"/
  DATA KW(6)%TYP /"T:B"/
  DATA KW(6)%DSCRPT /"*! Exchange-correlation functional !*"/

  DATA KW(7)%LABEL /"KERNEL_CUTOFF"/
  DATA KW(7)%TYP /"D:B"/
  DATA KW(7)%DSCRPT /"*! Density kernel radius in Bohr !*"/

  DATA KW(8)%LABEL /"MAXIT_NGWF_CG"/
  DATA KW(8)%TYP /"I:I"/
  DATA KW(8)%DSCRPT /"*! Max number of NGWF conjugate gradients (CG) &
       &iterations !*"/

  DATA KW(9)%LABEL /"MAXIT_LNV"/
  DATA KW(9)%TYP /"I:I"/
  DATA KW(9)%DSCRPT /"*! Max number of LNV iterations !*"/

  DATA KW(10)%LABEL /"LNV_THRESHOLD_ORIG"/
  DATA KW(10)%TYP /"D:I"/
  DATA KW(10)%DSCRPT /"*! LNV convergence threshold !*"/

  DATA KW(11)%LABEL /"NGWF_THRESHOLD_ORIG"/
  DATA KW(11)%TYP /"D:I"/
  DATA KW(11)%DSCRPT /"*! NGWF convergence threshold !*"/

  DATA KW(12)%LABEL /"OVLP_FOR_NONLOCAL"/
  DATA KW(12)%TYP /"L:E"/
  DATA KW(12)%DSCRPT /"*! Overlap sparsity for nonlocal !*"/

  DATA KW(13)%LABEL /"EXACT_LNV"/
  DATA KW(13)%TYP /"L:I"/
  DATA KW(13)%DSCRPT /"*! Use original LNV algorithm !*"/

  DATA KW(14)%LABEL /"PRECOND_RECIP"/
  DATA KW(14)%TYP /"L:E"/
  DATA KW(14)%DSCRPT /"*! Recip-space KE preconditioning !*"/

  DATA KW(15)%LABEL /"PRECOND_SCHEME"/
  DATA KW(15)%TYP /"T:E"/
  DATA KW(15)%DSCRPT /"*! Recip-space preconditioning scheme !*&
       & BG = Bowler-Gillan method;&
       & MAURI = Mauri method;&
       & TETER = Teter-Allen-Payne method"/

  DATA KW(16)%LABEL /"K_ZERO"/
  DATA KW(16)%TYP /"D:E"/
  DATA KW(16)%DSCRPT /"*! KE preconditioning parameter !*"/

  DATA KW(17)%LABEL /"PRECOND_REAL"/
  DATA KW(17)%TYP /"L:E"/
  DATA KW(17)%DSCRPT /"*! Real-space KE preconditioning !*"/

  DATA KW(18)%LABEL /"OUTPUT_DETAIL"/
  DATA KW(18)%TYP /"T:B"/
  DATA KW(18)%DSCRPT /"*! Level of output detail !*&
       &BRIEF, NORMAL or VERBOSE"/

  DATA KW(19)%LABEL /"KERNEL_UPDATE"/
  DATA KW(19)%TYP /"L:E"/
  DATA KW(19)%DSCRPT /"*! Update density kernel during NGWF line search !*"/

  DATA KW(20)%LABEL /"MAXIT_PEN"/
  DATA KW(20)%TYP /"I:I"/
  DATA KW(20)%DSCRPT /"*! Max number of penalty functional iterations !*"/

  DATA KW(21)%LABEL /"PEN_PARAM"/
  DATA KW(21)%TYP /"D:I"/
  DATA KW(21)%DSCRPT /"*! Penalty functional parameter !*"/

  DATA KW(22)%LABEL /"R_PRECOND"/
  DATA KW(22)%TYP /"D:E"/
  DATA KW(22)%DSCRPT /"*! Radial cut-off for real-space preconditioner !*"/

  DATA KW(23)%LABEL /"WRITE_NGWF_PLOT"/
  DATA KW(23)%TYP /"L:B"/
  DATA KW(23)%DSCRPT /"*! Write NGWFs in plotting format !*"/

  DATA KW(24)%LABEL /"WRITE_NGWF_GRAD_PLOT"/
  DATA KW(24)%TYP /"L:B"/
  DATA KW(24)%DSCRPT /"*! Write NGWF Gradients in plotting format !*"/

  DATA KW(25)%LABEL /"MAXIT_HOTELLING"/
  DATA KW(25)%TYP /"I:I"/
  DATA KW(25)%DSCRPT /"*! Number of Hotelling iteration per NGWF change !*"/

  DATA KW(26)%LABEL /"LOCPOT_INT_BATCH_SIZE"/
  DATA KW(26)%TYP /"I:E"/
  DATA KW(26)%DSCRPT /"*! Size of NGWF in fftbox batch for locpot integral evaluation !*"/

  DATA KW(27)%LABEL /"KINETIC_INT_BATCH_SIZE"/
  DATA KW(27)%TYP /"I:E"/
  DATA KW(27)%DSCRPT /"*! Size of NGWF in fftbox batch for kinetic integral evaluation !*"/

  DATA KW(28)%LABEL /"DENSITY_BATCH_SIZE"/
  DATA KW(28)%TYP /"I:E"/
  DATA KW(28)%DSCRPT /"*! Size of NGWF in fftbox batch for charge density evaluation !*"/

  DATA KW(29)%LABEL /"NGWF_GRAD_BATCH_SIZE"/
  DATA KW(29)%TYP /"I:E"/
  DATA KW(29)%DSCRPT /"*! Size of NGWF in fftbox batch for NGWF gradient calculation !*"/

  DATA KW(30)%LABEL /"USE_SPACE_FILLING_CURVE"/
  DATA KW(30)%TYP /"L:E"/
  DATA KW(30)%DSCRPT /"*! Re-arrange atoms according to space-filling curve !*"/

  DATA KW(31)%LABEL /"COREHAM_DENSKERN_GUESS"/
  DATA KW(31)%TYP /"L:E"/
  DATA KW(31)%DSCRPT /"*! Initial guess for density kernel from core Hamiltonian !*"/

  DATA KW(32)%LABEL /"WRITE_DENSKERN"/
  DATA KW(32)%TYP /"L:B"/
  DATA KW(32)%DSCRPT /"*! Write density kernel restart information !*"/

  DATA KW(33)%LABEL /"READ_DENSKERN"/
  DATA KW(33)%TYP /"L:B"/
  DATA KW(33)%DSCRPT /"*! Read density kernel restart information !*"/

  DATA KW(34)%LABEL /"NGWF_CG_TYPE"/
  DATA KW(34)%TYP /"T:E"/
  DATA KW(34)%DSCRPT /"*! Type of CG coefficient for NGWF optimisation !*&
      & NGWF_POLAK = Polak-Ribbiere formula; &
      & NGWF_FLETCHER = Fletcher-Reeves formula."/

  DATA KW(35)%LABEL /"LNV_CG_TYPE"/
  DATA KW(35)%TYP /"T:E"/
  DATA KW(35)%DSCRPT /"*! Type of CG coefficient for LNV denskern optimisation !*&
       & LNV_POLAK = Polak-Ribbiere formula; &
       & LNV_FLETCHER = Fletcher-Reeves formula."/

  DATA KW(36)%LABEL /"WRITE_TIGHTBOX_NGWFS"/
  DATA KW(36)%TYP /"L:B"/
  DATA KW(36)%DSCRPT /"*! Write in universal tightbox NGWFs restart information !*"/

  DATA KW(37)%LABEL /"READ_TIGHTBOX_NGWFS"/
  DATA KW(37)%TYP /"L:B"/
  DATA KW(37)%DSCRPT /"*! Read in universal tightbox NGWFs restart information !*"/

  DATA KW(38)%LABEL /"WRITE_DENSITY_PLOT"/
  DATA KW(38)%TYP /"L:B"/
  DATA KW(38)%DSCRPT /"*! Write the charge density in plotting format !*"/

  DATA KW(39)%LABEL /"FFTBOX_PREF"/
  DATA KW(39)%TYP /"T:I"/
  DATA KW(39)%DSCRPT /"*! Preferred FFT box dimensions !*"/

  DATA KW(40)%LABEL /"TIMINGS_LEVEL"/
  DATA KW(40)%TYP /"I:I"/
  DATA KW(40)%DSCRPT /"*! Level of timings output: 0(none) 1(nodes summary), 2(node details) !*"/

  DATA KW(41)%LABEL /"CONSTANT_EFIELD"/
  DATA KW(41)%TYP /"T:B"/
  DATA KW(41)%DSCRPT /"*! Cartesian coordinates of constant electric field vector !*"/

  DATA KW(42)%LABEL /"CHARGE"/
  DATA KW(42)%TYP /"D:B"/ ! pa: changed from I to D
  DATA KW(42)%DSCRPT /"*! The total charge of the system !*"/

  DATA KW(43)%LABEL /"WRITE_FORCES"/
  DATA KW(43)%TYP /"L:B"/
  DATA KW(43)%DSCRPT /"*! Write ionic forces !*"/

  DATA KW(44)%LABEL /"WRITE_POSITIONS"/
  DATA KW(44)%TYP /"L:B"/
  DATA KW(44)%DSCRPT /"*! Write ionic positions each geometry or MD step !*"/

  DATA KW(45)%LABEL /"TASK"/
  DATA KW(45)%TYP /"T:B"/
  DATA KW(45)%DSCRPT /"*! Type of calculation !*"/

  DATA KW(46)%LABEL /"CUTOFF_ENERGY"/
  DATA KW(46)%TYP /"P:B"/
  DATA KW(46)%DSCRPT /"*! Plane wave kinetic energy cutoff !*"/

  DATA KW(47)%LABEL /"NUM_EIGENVALUES"/
  DATA KW(47)%TYP /"I:I"/
  DATA KW(47)%DSCRPT /"*! Number of energy and occupancy eigenvalues to print &
       &below and above the Fermi level !*"/

  DATA KW(48)%LABEL /"OLD_LNV"/
  DATA KW(48)%TYP /"L:E"/
  DATA KW(48)%DSCRPT /"*! Use LNV algorithm backwards compatible pre Dec 2004 !*"/

  DATA KW(49)%LABEL /"SMOOTH_PROJECTORS"/
  DATA KW(49)%TYP /"D:E"/
  DATA KW(49)%DSCRPT /"*! Halfwidth of Gaussian filter for nonlocal projectors !*"/

  DATA KW(50)%LABEL /"SMOOTH_LOC_PSPOT"/
  DATA KW(50)%TYP /"D:E"/
  DATA KW(50)%DSCRPT /"*! Halfwidth of Gaussian filter for local pseudopotential !*"/

  DATA KW(51)%LABEL /"OCC_MIX"/
  DATA KW(51)%TYP /"D:E"/
  DATA KW(51)%DSCRPT /"*! Mix fraction of occupancy preconditioned NGWF cov grad !*"/

  DATA KW(52)%LABEL /"NGWF_CG_MAX_STEP"/
  DATA KW(52)%TYP /"D:E"/
  DATA KW(52)%DSCRPT /"*! Maximum length of trial step for NGWF optimisation line search !*"/

  DATA KW(53)%LABEL /"LNV_CG_MAX_STEP"/
  DATA KW(53)%TYP /"D:E"/
  DATA KW(53)%DSCRPT /"*! Maximum length of trial step for kernel optimisation line search !*"/

  DATA KW(54)%LABEL /"GEOM_MODULUS_EST"/
  DATA KW(54)%TYP /"P:I"/
  DATA KW(54)%DSCRPT /"*! The estimated bulk modulus !*"/

  DATA KW(55)%LABEL /"GEOM_FREQUENCY_EST"/
  DATA KW(55)%TYP /"P:I"/
  DATA KW(55)%DSCRPT /"*! The estimated average phonon frequency at the gamma point !*"/

  DATA KW(56)%LABEL /"GEOM_ENERGY_TOL"/
  DATA KW(56)%TYP /"P:I"/
  DATA KW(56)%DSCRPT /"*! Geometry optimization energy convergence tolerance.&
       & The difference between max and min energies over&
       & geom_convergence_win iterations must be less than this !*"/

  DATA KW(57)%LABEL /"GEOM_FORCE_TOL"/
  DATA KW(57)%TYP /"P:I"/
  DATA KW(57)%DSCRPT /"*! Geometry optimization force convergence tolerance !*"/

  DATA KW(58)%LABEL /"GEOM_DISP_TOL"/
  DATA KW(58)%TYP /"P:I"/
  DATA KW(58)%DSCRPT /"*! Geometry optimization displacement convergence tolerance !*"/

  DATA KW(59)%LABEL /"GEOM_MAX_ITER"/
  DATA KW(59)%TYP /"I:B"/
  DATA KW(59)%DSCRPT /"*! Maximum number of geometry optimization iterations !*"/

  DATA KW(60)%LABEL /"GEOM_CONVERGENCE_WIN"/
  DATA KW(60)%TYP /"I:I"/
  DATA KW(60)%DSCRPT /"*! Geometry optimization convergence tolerance window.&
       & The geometry optimization convergence criteria must all be met for&
       & geom_convergence_win iterations before acceptance !*"/

  DATA KW(61)%LABEL /"GEOM_CONTINUATION"/
  DATA KW(61)%TYP /"L:I"/
  DATA KW(61)%DSCRPT /"*! Read information for continuation of a previous geometry optimisation !*"/

  DATA KW(62)%LABEL /"GEOM_BACKUP_ITER"/
  DATA KW(62)%TYP /"I:I"/
  DATA KW(62)%DSCRPT /"*! Number of geometry optimisation iterations between backups of all data for continuation !*"/

  DATA KW(63)%LABEL /"DELTA_E_CONV"/
  DATA KW(63)%TYP /"L:E"/
  DATA KW(63)%DSCRPT /"*! Use consecutive energy gains as a criterion for NGWF convergence !*"/

  DATA KW(64)%LABEL /"PRINT_QC"/
  DATA KW(64)%TYP /"L:E"/
  DATA KW(64)%DSCRPT /"*! Print Quality Control information !*"/

  DATA KW(65)%LABEL /"MINIT_LNV"/
  DATA KW(65)%TYP /"I:I"/
  DATA KW(65)%DSCRPT /"*! Min number of LNV iterations !*"/

  DATA KW(66)%LABEL /"MAX_RESID_HOTELLING"/
  DATA KW(66)%TYP /"D:E"/
  DATA KW(66)%DSCRPT /"*! Max allowed value in Hotelling residual !*"/

  DATA KW(67)%LABEL /"CUBE_FORMAT"/
  DATA KW(67)%TYP /"L:B"/
  DATA KW(67)%DSCRPT /"*! Allow .cube format for plot outputs !*"/

  DATA KW(68)%LABEL /"GRD_FORMAT"/
  DATA KW(68)%TYP /"L:B"/
  DATA KW(68)%DSCRPT /"*! Allow .grd format for plot outputs !*"/

  DATA KW(69)%LABEL /"OLD_INPUT"/
  DATA KW(69)%TYP /"L:I"/
  DATA KW(69)%DSCRPT /"*! Read old-format input file !*"/

  DATA KW(70)%LABEL /"POSITIONS_ABS"/
  DATA KW(70)%TYP /"B:B"/
  DATA KW(70)%DSCRPT /"*! Cartesian positions for each atom !*"/

  DATA KW(71)%LABEL /"SPECIES"/
  DATA KW(71)%TYP /"B:B"/
  DATA KW(71)%DSCRPT /"*! Species information (symbol, atomic number, number of NGWFs, NGWF radius) !*"/

  DATA KW(72)%LABEL /"SPECIES_POT"/
  DATA KW(72)%TYP /"B:B"/
  DATA KW(72)%DSCRPT /"*! Pseudopotential name for each species !*"/

  DATA KW(73)%LABEL /"SPECIES_ATOMIC_SET"/
  DATA KW(73)%TYP /"B:I"/
  DATA KW(73)%DSCRPT /"*! Atomic set name for each species !*"/

  DATA KW(74)%LABEL /"SPECIES_CONSTRAINTS"/
  DATA KW(74)%TYP /"B:B"/
  DATA KW(74)%DSCRPT /"*! Ionic constraints for each species !*"/

  DATA KW(75)%LABEL /"ODD_PSINC_GRID"/
  DATA KW(75)%TYP /"L:E"/
  DATA KW(75)%DSCRPT /"*! Force odd number of points in simcell psinc grid !*"/

  DATA KW(76)%LABEL /"GEOM_PRINT_INV_HESSIAN"/
  DATA KW(76)%TYP /"L:E"/
  DATA KW(76)%DSCRPT /"*! Write inverse Hessian to standard output !*"/

  DATA KW(77)%LABEL /"VERBOSE_EWALD_FORCES"/
  DATA KW(77)%TYP /"L:E"/
  DATA KW(77)%DSCRPT /"*! Write debugging information for Ewald forces to standard output !*"/

  DATA KW(78)%LABEL /"GEOM_METHOD"/
  DATA KW(78)%TYP /"T:B"/
  DATA KW(78)%DSCRPT /"*! Method for geometry optimization !*&
       & CARTESIAN = BFGS from CASTEP; &
       & DELOCALIZED = DELOCALIZED INTERNALS from CASTEP."/

  DATA KW(79)%LABEL /"NNHO"/
  DATA KW(79)%TYP /"L:B"/
  DATA KW(79)%DSCRPT /"*! Initialise NGWFs to nonorthogonal natural hybrid orbitals !*"/

  DATA KW(80)%LABEL /"ELEC_CG_MAX"/
  DATA KW(80)%TYP /"I:E"/
  DATA KW(80)%DSCRPT /"*! Number of NGWF iterations to reset CG !*"/

  DATA KW(81)%LABEL /"CHECK_ATOMS"/
  DATA KW(81)%TYP /"L:E"/
  DATA KW(81)%DSCRPT /"*! Check atoms on top of each other !*"/

  DATA KW(82)%LABEL /"SPECIES_NGWF_PLOT"/
  DATA KW(82)%TYP /"B:B"/
  DATA KW(82)%DSCRPT /"*! Species whose NGWFs to plot !*"/

  DATA KW(83)%LABEL /"TSSEARCH_METHOD"/
  DATA KW(83)%TYP /"T:I"/
  DATA KW(83)%DSCRPT /"*! Specifies method to be used for TS search (e.g., LSTQST !*"/

  DATA KW(84)%LABEL /"TSSEARCH_LSTQST_PROTOCOL"/
  DATA KW(84)%TYP /"T:I"/
  DATA KW(84)%DSCRPT /"*! Specifies LSTQST protocol !*"/

  DATA KW(85)%LABEL /"TSSEARCH_QST_MAX_ITER"/
  DATA KW(85)%TYP /"I:E"/
  DATA KW(85)%DSCRPT /"*! Specifies maximum number of QST steps !*"/

  DATA KW(86)%LABEL /"TSSEARCH_CG_MAX_ITER"/
  DATA KW(86)%TYP /"I:E"/
  DATA KW(86)%DSCRPT /"*! Specifies maximum number of CG steps !*"/

  DATA KW(87)%LABEL /"TSSEARCH_FORCE_TOL"/
  DATA KW(87)%TYP /"P:I"/
  DATA KW(87)%DSCRPT /"*! Force tolerance for TS search !*"/

  DATA KW(88)%LABEL /"TSSEARCH_DISP_TOL"/
  DATA KW(88)%TYP /"P:I"/
  DATA KW(88)%DSCRPT /"*! Displacement tolerance for TS search !*"/

  DATA KW(89)%LABEL /"POSITIONS_ABS_PRODUCT"/
  DATA KW(89)%TYP /"B:I"/
  DATA KW(89)%DSCRPT /"*! Cartesian positions for each atom in the product (TS search)!*"/

  DATA KW(90)%LABEL /"POSITIONS_ABS_INTERMEDIATE"/
  DATA KW(90)%TYP /"B:I"/
  DATA KW(90)%DSCRPT /"*! Cartesian positions for each atom in the intermediate structure (TS search)!*"/

  DATA KW(91)%LABEL /"NGWF_HALO"/
  DATA KW(91)%TYP /"D:E"/
  DATA KW(91)%DSCRPT /"*! Halo extension to NGWF radii !*"/

  DATA KW(92)%LABEL /"HOMO_DENS_PLOT"/
  DATA KW(92)%TYP /"I:B"/
  DATA KW(92)%DSCRPT /"*! Number of squared MOs to plot from HOMO and lower !*"/

  DATA KW(93)%LABEL /"LUMO_DENS_PLOT"/
  DATA KW(93)%TYP /"I:B"/
  DATA KW(93)%DSCRPT /"*! Number of squared MOs to plot from LUMO and higher !*"/

  DATA KW(94)%LABEL /"HOMO_PLOT"/
  DATA KW(94)%TYP /"I:B"/
  DATA KW(94)%DSCRPT /"*! Number of MOs to plot from HOMO and lower !*"/

  DATA KW(95)%LABEL /"LUMO_PLOT"/
  DATA KW(95)%TYP /"I:B"/
  DATA KW(95)%DSCRPT /"*! Number of MOs to plot from LUMO and higher !*"/

  DATA KW(96)%LABEL /"DOS_SMEAR"/
  DATA KW(96)%TYP /"P:I"/
  DATA KW(96)%DSCRPT /"*! Half width of smearing Gaussians for DOS !*"/

  DATA KW(97)%LABEL /"LDOS_SMEAR"/
  DATA KW(97)%TYP /"P:I"/
  DATA KW(97)%DSCRPT /"*! Lorentzian broadening of LDOS !*"/

  DATA KW(98)%LABEL /"MD_DELTA_T"/
  DATA KW(98)%TYP /"P:B"/
  DATA KW(98)%DSCRPT /"*! Molecular dynamics time step !*"/

  DATA KW(99)%LABEL /"MD_NUM_ITER"/
  DATA KW(99)%TYP /"I:B"/
  DATA KW(99)%DSCRPT /"*! Maximum number of molecular dynamics iterations !*"/

  DATA KW(100)%LABEL /"MTS_XI"/
  DATA KW(100)%TYP /"L:E"/
  DATA KW(100)%DSCRPT/"*! Internal thermostat in the multiple time-step scheme !*"/ 
             
  DATA KW(101)%LABEL /"MTS_NSTEP"/
  DATA KW(101)%TYP /"I:E"/
  DATA KW(101)%DSCRPT/"*! Number of time steps in the multiple time-step scheme !*"/ 
             
  DATA KW(102)%LABEL /"MTS_NGWF_THRESHOLD"/
  DATA KW(102)%TYP /"D:E"/
  DATA KW(102)%DSCRPT/"*! NGWF convergence threshold for the mts correction !*"/ 
             
  DATA KW(103)%LABEL /"MTS_LNV_THRESHOLD"/
  DATA KW(103)%TYP /"D:E"/
  DATA KW(103)%DSCRPT/"*! LNV convergence threshold for the mts correction !*"/ 

  DATA KW(104)%LABEL /"SPIN_POLARISED"/
  DATA KW(104)%TYP /"L:B"/
  DATA KW(104)%DSCRPT /"*! Switch for spin polarisation !*"/

  DATA KW(105)%LABEL /"SPIN_POLARIZED"/
  DATA KW(105)%TYP /"L:B"/
  DATA KW(105)%DSCRPT /"*! Switch for spin polarization !*"/

  DATA KW(106)%LABEL /"SPIN"/
  DATA KW(106)%TYP /"I:B"/
  DATA KW(106)%DSCRPT /"*! Total spin of system !*"/

  DATA KW(107)%LABEL /"MAXIT_PALSER_MANO"/
  DATA KW(107)%TYP /"I:I"/
  DATA KW(107)%DSCRPT /"*! Maximum number of iterations for Palser-Manolopoulos scheme !*"/

  DATA KW(108)%LABEL /"MAXIT_KERNEL_FIX"/
  DATA KW(108)%TYP /"I:I"/
  DATA KW(108)%DSCRPT /"*! Maximum # iterations of Penalty Functional idempotency correction per LNV step !*"/

  DATA KW(109)%LABEL /"DO_PROPERTIES"/
  DATA KW(109)%TYP /"L:B"/
  DATA KW(109)%DSCRPT /"*! Allow calculation of properties !*"/

  DATA KW(110)%LABEL /"POPN_CALCULATE"/
  DATA KW(110)%TYP /"L:B"/
  DATA KW(110)%DSCRPT /"*! Allow population analysis !*"/

  DATA KW(111)%LABEL /"POPN_BOND_CUTOFF"/
  DATA KW(111)%TYP /"P:B"/
  DATA KW(111)%DSCRPT /"*! Bond length cutoff for population analysis !*"/

  DATA KW(112)%LABEL /"DEVEL_CODE"/
  DATA KW(112)%TYP /"T:E"/
  DATA KW(112)%DSCRPT /"*! For development code only !*"/

  DATA KW(113)%LABEL /"WRITE_XYZ"/
  DATA KW(113)%TYP /"L:B"/
  DATA KW(113)%DSCRPT /"*! Output coordinates in .xyz file !*"/

  DATA KW(114)%LABEL /"WRITE_PARAMS"/
  DATA KW(114)%TYP /"L:B"/
  DATA KW(114)%DSCRPT /"*! Output runtime parameters at startup !*"/

  DATA KW(115)%LABEL /"DENSE_THRESHOLD"/
  DATA KW(115)%TYP /"D:E"/
  DATA KW(115)%DSCRPT /"*! Threshold for matrix segment filling for segment  to be dense!*"/

  DATA KW(116)%LABEL /"NGWF_ANALYSIS"/
  DATA KW(116)%TYP /"L:B"/
  DATA KW(116)%DSCRPT /"*! Perform NGWF analysis !*"/

  DATA KW(117)%LABEL /"SPREAD_CALCULATE"/
  DATA KW(117)%TYP /"L:B"/
  DATA KW(117)%DSCRPT /"*! Calculate spread of NGWFs !*"/

  DATA KW(118)%LABEL /"DISPERSION"/
  DATA KW(118)%TYP /"I:B"/
  DATA KW(118)%DSCRPT /"*! Select dispersion correction !*"/

  DATA KW(119)%LABEL /"PADDED_LATTICE_CART"/
  DATA KW(119)%TYP /"B:I"/
  DATA KW(119)%DSCRPT/"*! The simulation cell lattice vectors for the padded cell !*"/

  DATA KW(120)%LABEL /"COULOMB_CUTOFF_TYPE"/
  DATA KW(120)%TYP /"T:I"/
  DATA KW(120)%DSCRPT/"*! Type of cutoff coulomb interaction: NONE, SPHERE, CYLINDER, SLAB, WIRE !*"/

  DATA KW(121)%LABEL /"COULOMB_CUTOFF_RADIUS"/
  DATA KW(121)%TYP /"P:I"/
  DATA KW(121)%DSCRPT/"*! Radius of sphere or cylinder for cutoff coulomb interaction !*"/

  DATA KW(122)%LABEL /"COULOMB_CUTOFF_LENGTH"/
  DATA KW(122)%TYP /"P:I"/
  DATA KW(122)%DSCRPT/"*! Length of cylinder or width of slab for cutoff coulomb interaction !*"/

  DATA KW(123)%LABEL /"COULOMB_CUTOFF_WRITE_INT"/
  DATA KW(123)%TYP /"L:I"/
  DATA KW(123)%DSCRPT/"*! Write real-space cutoff Coulomb interaction scalarfield !*"/

  DATA KW(124)%LABEL /"LOCPOT_SCHEME"/
  DATA KW(124)%TYP /"T:E"/
  DATA KW(124)%DSCRPT/"*! Scheme for evaluating local potential matrix elements !*&
      & FULL = Calculate matrix and symmetrize;&
      & LOWER = Calculate lower triangle only and expand;&
      & ALTERNATE = Calculate alternating elements from both triangles and expand !*"/

  DATA KW(125)%LABEL /"LNV_CHECK_TRIAL_STEPS"/
  DATA KW(125)%TYP /"L:I"/
  DATA KW(125)%DSCRPT /"*! Check stability of kernel at each trial step during LNV !*"/

  DATA KW(126)%LABEL /"BS_KPOINT_PATH"/
  DATA KW(126)%TYP /"B:B"/
  DATA KW(126)%DSCRPT /"*! K-point path for bandstructure calculation !*"/

  DATA KW(127)%LABEL /"BS_KPOINT_PATH_SPACING"/
  DATA KW(127)%TYP /"P:I"/
  DATA KW(127)%DSCRPT /"*! K-point spacing along bandstructure path !*"/

  DATA KW(128)%LABEL /"BS_NUM_EIGENVALUES"/
  DATA KW(128)%TYP /"I:I"/
  DATA KW(128)%DSCRPT /"*! Num of energy and occ. eigenvalues to print below and above the Fermi level from a bs calc !*"/

  DATA KW(129)%LABEL /"BS_UNFOLD"/
  DATA KW(129)%TYP /"T:I"/
  DATA KW(129)%DSCRPT /"*! Number of times to unfold Brillouin zone in each lattice direction !*"/

  DATA KW(130)%LABEL /"BS_METHOD"/
  DATA KW(130)%TYP /"T:I"/
  DATA KW(130)%DSCRPT /"*! Method to use: 'KP' or 'TB' !*"/

  DATA KW(131)%LABEL /"GEOM_REUSE_DK_NGWFS"/
  DATA KW(131)%TYP /"L:E"/
  DATA KW(131)%DSCRPT /"*! Re-use density kernel and NGWFs during geometry optimisation steps !*"/

  DATA KW(132)%LABEL /"WRITE_CONVERGED_DK_NGWFS"/
  DATA KW(132)%TYP /"L:I"/
  DATA KW(132)%DSCRPT /"*! Only write Density Kernel and NGWFs upon convergence of NGWF optimisation !*"/

  DATA KW(133)%LABEL /"WRITE_SW_NGWFS"/
  DATA KW(133)%TYP /"L:B"/
  DATA KW(133)%DSCRPT /"*! Write NGWFs restart information in spherical waves representation !*"/

  DATA KW(134)%LABEL /"READ_SW_NGWFS"/
  DATA KW(134)%TYP /"L:B"/
  DATA KW(134)%DSCRPT /"*! Read NGWFs restart information in spherical waves representation !*"/

  DATA KW(135)%LABEL /"WRITE_MAX_L"/
  DATA KW(135)%TYP /"I:I"/
  DATA KW(135)%DSCRPT /"*! Maximum angular momentum number when writing in SW representation !*"/

  DATA KW(136)%LABEL /"READ_MAX_L"/
  DATA KW(136)%TYP /"I:I"/
  DATA KW(136)%DSCRPT /"*! Maximum angular momentum number when reading in SW representation !*"/

  DATA KW(137)%LABEL /"EXTRA_N_SW"/
  DATA KW(137)%TYP /"I:I"/
  DATA KW(137)%DSCRPT /"*! Maximum number of zeros of the sph Bessel function for the SW !*"/

  DATA KW(138)%LABEL /"CLASSICAL_INFO"/
  DATA KW(138)%TYP /"B:B"/
  DATA KW(138)%DSCRPT /"*! Include classical atoms Coulomb interaction in the calculation !*"/

  DATA KW(139)%LABEL /"HUBBARD_MAX_ITER"/
  DATA KW(139)%TYP /"I:I"/
  DATA KW(139)%DSCRPT /"*! Maximum number of DFT+U projector optimisation steps, 0 for none !*"/

  DATA KW(140)%LABEL /"HUBBARD_ENERGY_TOL"/
  DATA KW(140)%TYP /"P:I"/
  DATA KW(140)%DSCRPT /"*! Energy tolerance when using DFT+U projector optimisation !*"/

  DATA KW(141)%LABEL /"HUBBARD_CONV_WIN"/
  DATA KW(141)%TYP /"I:I"/
  DATA KW(141)%DSCRPT /"*! Energy convergence window when using DFT+U projector optimisation !*"/

  DATA KW(142)%LABEL /"HUBBARD_PROJ_MIXING"/
  DATA KW(142)%TYP /"D:I"/
  DATA KW(142)%DSCRPT /"*! Proportion of old Hubbard projector to mix with new !*"/

  DATA KW(143)%LABEL /"HUBBARD"/
  DATA KW(143)%TYP /"B:I"/
  DATA KW(143)%DSCRPT /"*! Hubbard species info (symb., ang. mom., U parameter (eV), effective charge, alpha parameter (eV)) !*"/

  DATA KW(144)%LABEL /"HUBBARD_FUNCTIONAL"/
  DATA KW(144)%TYP /"I:I"/
  DATA KW(144)%DSCRPT /"*! DFT+U energy functional to use !*"/

  DATA KW(145)%LABEL /"HUBBARD_TENSOR_CORR"/
  DATA KW(145)%TYP /"I:I"/
  DATA KW(145)%DSCRPT /"*! DFT+U projector tensorial correction to use !*"/

  DATA KW(146)%LABEL /"HUBBARD_NGWF_SPIN_THRESHOLD"/
  DATA KW(146)%TYP /"D:I"/
  DATA KW(146)%DSCRPT /"*! NGWF RMS  gradient threshold at which to switch off DFT+U spin-splitting !*"/

  DATA KW(147)%LABEL /"POLARISATION_CALCULATE"/
  DATA KW(147)%TYP /"L:B"/
  DATA KW(147)%DSCRPT /"*! Allow calculation of polarisation !*"/

  DATA KW(148)%LABEL /"KERFIX"/
  DATA KW(148)%TYP /"I:E"/
  DATA KW(148)%DSCRPT /"*! Density kernel fixing approach !*"/

  DATA KW(149)%LABEL /"PAW"/
  DATA KW(149)%TYP /"L:E"/
  DATA KW(149)%DSCRPT /"*! Uses a PAW construction to find correct core densities/wavefunctions !*"/

  DATA KW(150)%LABEL /"aug_grid_factor"/
  DATA KW(150)%TYP /"D:E"/
  DATA KW(150)%DSCRPT /"*! Sets the scale of the augmentation grid for PAW calculations !*"/

  DATA KW(151)%LABEL /"DO_TDDFT"/
  DATA KW(151)%TYP /"L:E"/
  DATA KW(151)%DSCRPT /"*! Allow Time-Dependent DFT calculation !*"/

  DATA KW(152)%LABEL /"TDDFT_MAXIMUM_ENERGY"/
  DATA KW(152)%TYP /"P:E"/
  DATA KW(152)%DSCRPT /"*! Desired maximum of spectrum from TDDFT (in Hartree) !*"/

  DATA KW(153)%LABEL /"TDDFT_RESOLUTION"/
  DATA KW(153)%TYP /"P:E"/
  DATA KW(153)%DSCRPT /"*! Desired resolution of spectrum from TDDFT (in Hartree) !*"/

  DATA KW(154)%LABEL /"TDDFT_PROPAGATION_METHOD"/
  DATA KW(154)%TYP /"T:E"/
  DATA KW(154)%DSCRPT /"*! Method used to integrate von Neumann equation eg. RUNGEKUTTA or CRANKNICHOLSON!*"/

  DATA KW(155)%LABEL /"TDDFT_SPARSITY_LEVEL"/
  DATA KW(155)%TYP /"I:E"/
  DATA KW(155)%DSCRPT /"*! Matrix sparsity when computing propagators e.g. 0(recommended),1,2,3 !*"/

  DATA KW(156)%LABEL /"TDDFT_TAMMDANCOFF"/
  DATA KW(156)%TYP /"L:E"/
  DATA KW(156)%DSCRPT /"*! Invoke Tamm-Dancoff decoupling approximation !*"/

  DATA KW(157)%LABEL /"TDDFT_DIPOLE_KICK_STRENGTH"/
  DATA KW(157)%TYP /"T:E"/
  DATA KW(157)%DSCRPT /"*! Maximum allowed phase shift in TDDFT delta-kick, units of PI !*"/

  DATA KW(158)%LABEL /"TDDFT_XC_FUNCTIONAL"/
  DATA KW(158)%TYP /"T:E"/
  DATA KW(158)%DSCRPT /"*! Exchange-correlation functional for TDDFT !*&
       & LDA = Adiabatic Perdew-Zunger LDA;&
       & NONE = Random Phase Approximation."/

  DATA KW(159)%LABEL /"TDDFT_HAMILTONIAN_MIXING"/
  DATA KW(159)%TYP /"I:E"/
  DATA KW(159)%DSCRPT /"*! Order of polynomial extrapolation to H(t + half Delta t) 0,1,2 !*"/

  DATA KW(160)%LABEL /"TDDFT_DAMPING"/
  DATA KW(160)%TYP /"P:E"/
  DATA KW(160)%DSCRPT /"*! Energy smearing when Fourier transforming for frequency-dependent dipole moment !*"/

  DATA KW(161)%LABEL /"TDDFT_ENFORCED_IDEMPOTENCY"/
  DATA KW(161)%TYP /"L:E"/
  DATA KW(161)%DSCRPT /"*! Project out at each timestep that part of change to denskern not respecting idempotency to 1st order !*"/

  DATA KW(162)%LABEL /"TDDFT_MAXIT_HOTELLING"/
  DATA KW(162)%TYP /"I:E"/
  DATA KW(162)%DSCRPT /"*! Number of Hotelling iteration per propagation step in Crank-Nicholson propagator !*"/

  DATA KW(163)%LABEL /"TDDFT_MAX_RESID_HOTELLING"/
  DATA KW(163)%TYP /"D:E"/
  DATA KW(163)%DSCRPT /"*! Max allowed value in Hotelling residual for Crank-Nicholson propagator !*"/

  DATA KW(164)%LABEL /"TDDFT_INV_OVERLAP_EXACT"/
  DATA KW(164)%TYP /"L:E"/
  DATA KW(164)%DSCRPT /"*! Renew inverse overlap with O N^3 algorithm before beginning TDDFT !*"/

  DATA KW(165)%LABEL /"PBC_CORRECTION_CUTOFF"/
  DATA KW(165)%TYP /"D:E"/
  DATA KW(165)%DSCRPT /"*! alpha*L cutoff parameter for Martyna-Tuckerman PBC correction !*"/

  DATA KW(166)%LABEL /"SPECIES_LDOS_GROUPS"/
  DATA KW(166)%TYP /"B:I"/
  DATA KW(166)%DSCRPT /"*! Species groups for Local density of states calculation !*"/

  DATA KW(167)%LABEL /"DX_FORMAT"/
  DATA KW(167)%TYP /"L:B"/
  DATA KW(167)%DSCRPT /"*! Allow .dx format for plot outputs !*"/

  DATA KW(168)%LABEL /"DX_FORMAT_DIGITS"/
  DATA KW(168)%TYP /"I:I"/
  DATA KW(168)%DSCRPT /"*! Number of significant digits in .dx output !*"/

  DATA KW(169)%LABEL /"DX_FORMAT_COARSE"/
  DATA KW(169)%TYP /"L:I"/
  DATA KW(169)%DSCRPT /"*! Output only points on the coarse grid !*"/

  DATA KW(170)%LABEL /"NGWF_DIIS_SIZE"/
  DATA KW(170)%TYP /"I:E"/
  DATA KW(170)%DSCRPT /"*! Max size of the DIIS iterative subspace !*"/

  DATA KW(171)%LABEL /"MAXIT_NGWF_DIIS"/
  DATA KW(171)%TYP /"I:E"/
  DATA KW(171)%DSCRPT /"*! Max number of NGWF DIIS steps !*"/

  DATA KW(172)%LABEL /"NGWF_DIIS_MAX_STEP"/
  DATA KW(172)%TYP /"D:E"/
  DATA KW(172)%DSCRPT /"*! Max length of trial step for NGWF DIIS optimisation !*"/

  DATA KW(173)%LABEL /"NGWF_DIIS_RESET"/
  DATA KW(173)%TYP /"I:E"/
  DATA KW(173)%DSCRPT /"*! Max number of iteration for the resetting of DIIS optimisation of the NGWFs !*"/

  DATA KW(174)%LABEL /"PAW_OUTPUT_DETAIL"/
  DATA KW(174)%TYP /"T:B"/
  DATA KW(174)%DSCRPT /"*! Level of output detail for PAW: BRIEF, NORMAL or VERBOSE !*"/

  DATA KW(175)%LABEL /"FINE_GRID_SCALE"/
  DATA KW(175)%TYP /"D:B"/
  DATA KW(175)%DSCRPT /"*! Ratio of size of fine grid to standard grid !*"/

  DATA KW(176)%LABEL /"SPECIES_COND"/
  DATA KW(176)%TYP /"B:B"/
  DATA KW(176)%DSCRPT /"*! Species information for Conduction NGWFs (symbol, atomic number, number of NGWFs, NGWF radius) !*"/

  DATA KW(177)%LABEL /"LIBXC_X_FUNC_ID"/
  DATA KW(177)%TYP /"I:I"/
  DATA KW(177)%DSCRPT /"*! Functional identifier for which exchange functional to use with LIBXC - see LIBXC documentation !*"/

  DATA KW(178)%LABEL /"LIBXC_C_FUNC_ID"/
  DATA KW(178)%TYP /"I:I"/
  DATA KW(178)%DSCRPT /"*! Functional identifier for which correlation functional to use with LIBXC - see LIBXC documentation !*"/

  DATA KW(179)%LABEL /"IS_DENSITY_THRESHOLD"/
  DATA KW(179)%TYP /"D:I"/
  DATA KW(179)%DSCRPT /"*! Implicit solvent: rho_0 parameter in Fattebert-Gygi functional !*"/

  DATA KW(180)%LABEL /"IS_SOLVATION_BETA"/
  DATA KW(180)%TYP /"D:I"/
  DATA KW(180)%DSCRPT /"*! Implicit solvent: beta parameter in Fattebert-Gygi functional !*"/

  DATA KW(181)%LABEL /"IS_BULK_PERMITTIVITY"/
  DATA KW(181)%TYP /"D:B"/
  DATA KW(181)%DSCRPT /"*! Implicit solvent: eps_inf parameter in Fattebert-Gygi functional !*"/

  DATA KW(182)%LABEL /"IS_MULTIGRID_ERROR_TOL"/
  DATA KW(182)%TYP /"D:I"/
  DATA KW(182)%DSCRPT /"*! Implicit solvent: stop criterion for the multigrid solver !*"/

  DATA KW(183)%LABEL /"IS_SMEARED_ION_WIDTH"/
  DATA KW(183)%TYP /"D:I"/
  DATA KW(183)%DSCRPT /"*! Implicit solvent: Width of Gaussian smearing in bohr !*"/

  DATA KW(184)%LABEL /"IS_IMPLICIT_SOLVENT"/
  DATA KW(184)%TYP /"L:B"/
  DATA KW(184)%DSCRPT /"*! Use implicit solvent? !*"/

  DATA KW(185)%LABEL /"IS_CHECK_SOLV_ENERGY_GRAD"/
  DATA KW(185)%TYP /"L:I"/
  DATA KW(185)%DSCRPT /"*! Implicit solvent: sanity check the energy gradient !*"/

  DATA KW(186)%LABEL /"IS_SMEARED_ION_REP"/
  DATA KW(186)%TYP /"L:I"/
  DATA KW(186)%DSCRPT /"*! Implicit solvent: use smeared ions for electrostatics !*"/

  DATA KW(187)%LABEL /"IS_SOLVATION_METHOD"/
  DATA KW(187)%TYP /"T:I"/
  DATA KW(187)%DSCRPT /"*! Implicit solvent: direct or corrective method !*"/

  DATA KW(188)%LABEL /"IS_DIELECTRIC_MODEL"/
  DATA KW(188)%TYP /"T:I"/
  DATA KW(188)%DSCRPT /"*! Implicit solvent: how the cavity is determined !*"/

  DATA KW(189)%LABEL /"IS_DIELECTRIC_FUNCTION"/
  DATA KW(189)%TYP /"T:I"/
  DATA KW(189)%DSCRPT /"*! Implicit solvent: how the cavity is determined !*"/

  DATA KW(190)%LABEL /"IS_SOLVATION_OUTPUT_DETAIL"/
  DATA KW(190)%TYP /"T:I"/
  DATA KW(190)%DSCRPT /"*! Implicit solvent: controls extra output !*"/

  DATA KW(191)%LABEL /"IS_SURFACE_THICKNESS"/
  DATA KW(191)%TYP /"D:E"/
  DATA KW(191)%DSCRPT /"*! Implicit solvent: thickness used for SA calculation !*"/

  DATA KW(192)%LABEL /"IS_SOLVENT_SURFACE_TENSION"/
  DATA KW(192)%TYP /"P:B"/
  DATA KW(192)%DSCRPT /"*! Implicit solvent: surface tension of solvent !*"/

  DATA KW(193)%LABEL /"IS_MULTIGRID_MAX_ITERS"/
  DATA KW(193)%TYP /"I:E"/
  DATA KW(193)%DSCRPT /"*! Implicit solvent: max number of iterations in multigrid solver !*"/

  DATA KW(194)%LABEL /"IS_MULTIGRID_NLEVELS"/
  DATA KW(194)%TYP /"I:E"/
  DATA KW(194)%DSCRPT /"*! Implicit solvent: number of levels in multigrid solver !*"/

  DATA KW(195)%LABEL /"IS_MULTIGRID_DEFECT_ERROR_TOL"/
  DATA KW(195)%TYP /"D:I"/
  DATA KW(195)%DSCRPT /"*! Implicit solvent: stop criterion for defect correction !*"/

  DATA KW(196)%LABEL /"IS_DISCRETIZATION_ORDER"/
  DATA KW(196)%TYP /"I:I"/
  DATA KW(196)%DSCRPT /"*! Implicit solvent: discretization order (2nd, 4th, ...) for the PB solver !*"/

  DATA KW(197)%LABEL /"IS_INCLUDE_CAVITATION"/
  DATA KW(197)%TYP /"L:B"/
  DATA KW(197)%DSCRPT/"*! Implicit solvent: include cavitation energy !*"/

  DATA KW(198)%LABEL /"IS_BC_COARSENESS"/
  DATA KW(198)%TYP /"I:I"/
  DATA KW(198)%DSCRPT/"*! Open BCs: controls boundary condition coarse-graining !*"/

  DATA KW(199)%LABEL /"IS_BC_SURFACE_COARSENESS"/
  DATA KW(199)%TYP /"I:I"/
  DATA KW(199)%DSCRPT/"*! Open BCs: controls boundary condition coarse-graining !*"/

  DATA KW(200)%LABEL /"IS_BC_THRESHOLD"/
  DATA KW(200)%TYP /"D:E"/
  DATA KW(200)%DSCRPT/"*! Open BCs: controls boundary condition coarse-graining !*"/

  DATA KW(201)%LABEL /"OPENBC_PSPOT"/
  DATA KW(201)%TYP /"L:I"/
  DATA KW(201)%DSCRPT/"*! Force open BCs in local pseudopotential !*"/

  DATA KW(202)%LABEL /"OPENBC_PSPOT_FINETUNE_F"/
  DATA KW(202)%TYP /"I:E"/
  DATA KW(202)%DSCRPT/"*! Open BCs in local pseudo, fineness parameter!*"/

  DATA KW(203)%LABEL /"OPENBC_PSPOT_FINETUNE_NPTSX"/
  DATA KW(203)%TYP /"I:E"/
  DATA KW(203)%DSCRPT/"*! Open BCs in local pseudo, npts_x parameter!*"/

  DATA KW(204)%LABEL /"OPENBC_PSPOT_FINETUNE_ALPHA"/
  DATA KW(204)%TYP /"D:E"/
  DATA KW(204)%DSCRPT/"*! Open BCs in local pseudo, alpha parameter!*"/

  DATA KW(205)%LABEL /"OPENBC_ION_ION"/
  DATA KW(205)%TYP /"L:I"/
  DATA KW(205)%DSCRPT/"*! Force open BCs in ion-ion energy !*"/

  DATA KW(206)%LABEL /"COND_NUM_STATES"/
  DATA KW(206)%TYP /"I:B"/
  DATA KW(206)%DSCRPT /"*! Number of conduction states to be optimised for !*"/ ! lr408

  DATA KW(207)%LABEL /"COND_READ_DENSKERN"/
  DATA KW(207)%TYP /"L:B"/
  DATA KW(207)%DSCRPT /"*! Read in conduction density kernel !*"/ ! lr408

  DATA KW(208)%LABEL /"COND_READ_TIGHTBOX_NGWFS"/
  DATA KW(208)%TYP /"L:B"/
  DATA KW(208)%DSCRPT /"*! Read in universal tightbox conduction NGWFs !*"/ ! lr408

  DATA KW(209)%LABEL /"COND_KERNEL_CUTOFF"/
  DATA KW(209)%TYP /"D:B"/
  DATA KW(209)%DSCRPT /"*! Conduction density kernel radius in Bohr !*"/ ! lr408

  DATA KW(210)%LABEL /"COND_INIT_SHIFT"/
  DATA KW(210)%TYP /"D:B"/
  DATA KW(210)%DSCRPT /"*! Initial shifting factor for projected conduction Hamiltonian !*"/ ! lr408

  DATA KW(211)%LABEL /"COND_SHIFT_BUFFER"/
  DATA KW(211)%TYP /"D:B"/
  DATA KW(211)%DSCRPT /"*! Additional buffer for updating projected Hamiltonian shift !*"/ ! lr408

  DATA KW(212)%LABEL /"COND_FIXED_SHIFT"/
  DATA KW(212)%TYP /"L:B"/
  DATA KW(212)%DSCRPT /"*! Fixed projected conduction Hamiltonian shift !*"/ ! lr408

  DATA KW(213)%LABEL /"COND_CALC_MAX_EIGEN"/
  DATA KW(213)%TYP /"L:B"/
  DATA KW(213)%DSCRPT /"*! Calculate maximum conduction Hamiltonian eigenvalue !*"/ ! lr408

  DATA KW(214)%LABEL /"COND_CALC_OPTICAL_SPECTRA"/
  DATA KW(214)%TYP /"L:B"/
  DATA KW(214)%DSCRPT /"*! Calculate matrix elements for optical absorption spectra !*"/ ! lr408

  DATA KW(215)%LABEL /"COND_SPEC_CALC_MOM_MAT_ELS"/
  DATA KW(215)%TYP /"L:B"/
  DATA KW(215)%DSCRPT /"*! Calculate momentum matrix elements (default true otherwise use position) !*"/ ! lr408

  DATA KW(216)%LABEL /"COND_SPEC_CALC_NONLOC_COMM"/
  DATA KW(216)%TYP /"L:B"/
  DATA KW(216)%DSCRPT /"*! Calculate nonlocal commutator for momentum matrix elements (default true) !*"/ ! lr408

  DATA KW(217)%LABEL /"COND_SPEC_CONT_DERIV"/
  DATA KW(217)%TYP /"L:B"/
  DATA KW(217)%DSCRPT /"*! Calculate non-local commuator using continuous deriv in k-space (default true) !*"/ ! lr408

  DATA KW(218)%LABEL /"COND_SPEC_NONLOC_COMM_SHIFT"/
  DATA KW(218)%TYP /"D:B"/
  DATA KW(218)%DSCRPT /"*! Finite difference shift for non-local commutator if calculating using finite difference !*"/ ! lr408

  DATA KW(219)%LABEL /"KERNEL_DIIS"/
  DATA KW(219)%TYP /"L:I"/
  DATA KW(219)%DSCRPT /"*! Use kernel DIIS in the inner loop !*"/

  DATA KW(220)%LABEL /"KERNEL_DIIS_MAX"/
  DATA KW(220)%TYP /"I:I"/
  DATA KW(220)%DSCRPT /"*! Max number of kernels saved during kernel DIIS !*"/

  DATA KW(221)%LABEL /"KERNEL_DIIS_MAXIT"/
  DATA KW(221)%TYP /"I:I"/
  DATA KW(221)%DSCRPT /"*! Max number of kernel DIIS iterations !*"/

  DATA KW(222)%LABEL /"KERNEL_DIIS_TYPE"/
  DATA KW(222)%TYP /"T:A"/
  DATA KW(222)%DSCRPT /"*! Density mixing method !*"/

  DATA KW(223)%LABEL /"KERNEL_DIIS_THRES"/
  DATA KW(223)%TYP /"D:A"/
  DATA KW(223)%DSCRPT /"*! Density mixing convergence threshold !*"/

  DATA KW(224)%LABEL /"KERNEL_DIIS_LITER"/
  DATA KW(224)%TYP /"I:A"/
  DATA KW(224)%DSCRPT /"*! Number of linear mixing iterations !*"/

  DATA KW(225)%LABEL /"KERNEL_DIIS_C_IN"/
  DATA KW(225)%TYP /"D:A"/
  DATA KW(225)%DSCRPT /"*! Coefficient for the input kernel in linear mixing DIIS !*"/

  DATA KW(226)%LABEL /"KERNEL_DIIS_CRITERIA"/
  DATA KW(226)%TYP /"T:A"/
  DATA KW(226)%DSCRPT /"*! Density mixing convergence criteria !*"/

  DATA KW(227)%LABEL /"KERNEL_CHRISTOFFEL_UPDATE"/
  DATA KW(227)%TYP /"L:I"/
  DATA KW(227)%DSCRPT /"*! Update density kernel during NGWF line search !*"/

  DATA KW(228)%LABEL /"ELEC_ENERGY_TOL"/
  DATA KW(228)%TYP /"P:I"/
  DATA KW(228)%DSCRPT /"*! Tolerance on total energy change during NGWF optimisation !*"/

  DATA KW(229)%LABEL /"ELEC_FORCE_TOL"/
  DATA KW(229)%TYP /"P:I"/
  DATA KW(229)%DSCRPT /"*! Tolerance on max force change during NGWF optimisation !*"/

  DATA KW(230)%LABEL /"NGWF_MAX_GRAD"/
  DATA KW(230)%TYP /"D:I"/
  DATA KW(230)%DSCRPT /"*! Maximum permissible value of NGWF Gradient for convergence !*"/

  DATA KW(231)%LABEL /"NONSC_FORCES"/
  DATA KW(231)%TYP /"L:A"/
  DATA KW(231)%DSCRPT /"*! Include non self-consistent forces due to NGWF optimisation !*"/

  DATA KW(232)%LABEL /"GEOM_RESET_DK_NGWFS_ITER"/
  DATA KW(232)%TYP /"I:I"/
  DATA KW(232)%DSCRPT /"*! Number of geom iterations between resets of kernel and NGWFs !*"/

  DATA KW(233)%LABEL /"VELOCITIES"/
  DATA KW(233)%TYP /"B:B"/
  DATA KW(233)%DSCRPT /"*! Initial velocities for each atom !*"/

  DATA KW(234)%LABEL /"MIX_DENSKERN_NUM"/
  DATA KW(234)%TYP /"I:E"/
  DATA KW(234)%DSCRPT /"*! Number of coefficients used to build new guess for dkn !*"/

  DATA KW(235)%LABEL /"MIX_NGWFS_NUM"/
  DATA KW(235)%TYP /"I:E"/
  DATA KW(235)%DSCRPT /"*! Number of coefficients used to build new guess for NGWFS !*"/

  DATA KW(236)%LABEL /"MIX_DENSKERN_TYPE"/
  DATA KW(236)%TYP /"I:E"/
  DATA KW(236)%DSCRPT /"*! Type of mixing used to build new guess for dkn !*"/
             
  DATA KW(237)%LABEL /"MIX_NGWFS_TYPE"/
  DATA KW(237)%TYP /"I:E"/
  DATA KW(237)%DSCRPT /"*! Type of mixing used to build new guess for NGWFS !*"/
             
  DATA KW(238)%LABEL /"MIX_NGWFS_COEFF"/
  DATA KW(238)%TYP /"D:E"/
  DATA KW(238)%DSCRPT /"*! Coefficient for mixing extrapolated NGWFS with NGWFS at last step !*"/
             
  DATA KW(239)%LABEL /"MIX_LOCAL_LENGTH"/
  DATA KW(239)%TYP /"P:E"/
  DATA KW(239)%DSCRPT /"*! Max radius for local mixing of NGWFs !*"/
             
  DATA KW(240)%LABEL /"MIX_LOCAL_SMEAR"/
  DATA KW(240)%TYP /"P:E"/
  DATA KW(240)%DSCRPT /"*! Radial smearing for local mixing of NGWFs !*"/
             
  DATA KW(241)%LABEL /"OPENBC_HARTREE"/
  DATA KW(241)%TYP /"L:I"/
  DATA KW(241)%DSCRPT/"*! Force open BCs in Hartree potential !*"/
             
  DATA KW(242)%LABEL /"VDW_PARAMS"/
  DATA KW(242)%TYP /"B:E"/
  DATA KW(242)%DSCRPT /"*! Replacement VDW parameters (atomic number, c6coeff, radzero, neff) !*"/
             
  DATA KW(243)%LABEL /"VDW_DCOEFF"/
  DATA KW(243)%TYP /"D:E"/
  DATA KW(243)%DSCRPT /"*! Replacement VDW damping coefficient !*"/
             
  DATA KW(244)%LABEL /"COMMS_GROUP_SIZE"/
  DATA KW(244)%TYP /"I:I"/
  DATA KW(244)%DSCRPT /"*! Number of nodes in a group (determines comms efficiency) !*"/
             
  DATA KW(245)%LABEL /"DBL_GRID_SCALE"/
  DATA KW(245)%TYP /"D:B"/
  DATA KW(245)%DSCRPT /"*! Ratio of charge density / potential working grid to standard grid (1 or 2 only) !*"/
             
  DATA KW(246)%LABEL /"REALSPACE_PROJECTORS"/
  DATA KW(246)%TYP /"L:I"/
  DATA KW(246)%DSCRPT /"*! Whether to evaluate and store projectors in real space !*"/
             
  DATA KW(247)%LABEL /"EVEN_PSINC_GRID"/
  DATA KW(247)%TYP /"L:E"/
  DATA KW(247)%DSCRPT /"*! Force even number of points in simcell psinc grid !*"/
             
  DATA KW(248)%LABEL /"MAXIT_LNV_COARSE"/
  DATA KW(248)%TYP /"I:E"/
  DATA KW(248)%DSCRPT /"*! Max number of low-accuracy LNV iterations !*"/
             
  DATA KW(249)%LABEL /"COND_PLOT_JOINT_ORBITALS"/
  DATA KW(249)%TYP /"L:B"/
  DATA KW(249)%DSCRPT /"*! Plot orbitals in the joint basis following a conduction calculation !*"/ ! lr408
             
  DATA KW(250)%LABEL /"COND_PLOT_VC_ORBITALS"/
  DATA KW(250)%TYP /"L:B"/
  DATA KW(250)%DSCRPT /"*! Plot orbitals in the val and cond NGWF basis sets after a cond calc !*"/ ! lr408
             
  DATA KW(251)%LABEL /"COND_NUM_EXTRA_STATES"/
  DATA KW(251)%TYP /"I:I"/
  DATA KW(251)%DSCRPT /"*! Number of extra conduction states for initial 'preconditioning' !*"/ ! lr408
             
  DATA KW(252)%LABEL /"COND_NUM_EXTRA_ITS"/
  DATA KW(252)%TYP /"I:I"/
  DATA KW(252)%DSCRPT /"*! Number of NGWF iterations with extra conduction states for 'preconditioning' !*"/ ! lr408
             
  DATA KW(253)%LABEL /"WRITE_NGWF_RADIAL"/
  DATA KW(253)%TYP /"I:B"/
  DATA KW(253)%DSCRPT/"*! Write NGWFs radial distributions !*"/
             
  DATA KW(254)%LABEL /"WRITE_NGWF_GRAD_RADIAL"/
  DATA KW(254)%TYP /"I:B"/
  DATA KW(254)%DSCRPT/"*! Write NGWFs gradients radial distributions !*"/
             
  DATA KW(255)%LABEL /"WRITE_RADIAL_STEP"/
  DATA KW(255)%TYP /"P:B"/
  DATA KW(255)%DSCRPT/"*! Define the grid step used for writing radial distributions !*"/
             
  DATA KW(256)%LABEL /"WRITE_RADIAL_SMEAR"/
  DATA KW(256)%TYP /"P:B"/
  DATA KW(256)%DSCRPT/"*! Define the gaussian smearing used for writing radial distributions !*"/
             
  DATA KW(257)%LABEL /"SMOOTH_SCHEME"/
  DATA KW(257)%TYP /"T:E"/
  DATA KW(257)%DSCRPT/"*! Smoothing scheme for the NGWF gradients at the edges !*"/
             
  DATA KW(258)%LABEL /"R_SMOOTH"/
  DATA KW(258)%TYP /"P:E"/
  DATA KW(258)%DSCRPT/"*! Radius of the unshaved NGWF gradients !*"/
             
  DATA KW(259)%LABEL /"K_SMOOTH"/
  DATA KW(259)%TYP /"D:E"/
  DATA KW(259)%DSCRPT/"*! Characteristic wavelength of the smoothing function !*"/
             
  DATA KW(260)%LABEL /"AUG_FUNCS_RECIP"/
  DATA KW(260)%TYP /"L:E"/
  DATA KW(260)%DSCRPT/"*! Construct Augmentation functions in recip space (T) or real (F) !*"/
             
  DATA KW(261)%LABEL /"INITIAL_DENS_REALSPACE"/
  DATA KW(261)%TYP /"L:E"/
  DATA KW(261)%DSCRPT/"*! Construct initial density in real space from atomsolver density !*"/
             
  ! lpl: NB62O output
  DATA KW(262)%LABEL /"WRITE_NBO"/
  DATA KW(262)%TYP /"L:B"/
  DATA KW(262)%DSCRPT /"*! Performs Natural Population Analysis and writes &
       &a FILE.47 input for GENNBO !*"/
           
  DATA KW(263)%LABEL /"NBO_WRITE_SPECIES"/
  DATA KW(263)%TYP /"B:B"/
  DATA KW(263)%DSCRPT /"*! List of atoms to be included in the output to &
       &GENNBO FILE.47 !*"/
           
  DATA KW(264)%LABEL /"NBO_SPECIES_NGWFLABEL"/
  DATA KW(264)%TYP /"B:I"/
  DATA KW(264)%DSCRPT /"*! User-specified NGWF (false) lm-label and &
       &NMB/NRBs for GENNBO !*"/
           
  DATA KW(265)%LABEL /"NBO_INIT_LCLOWDIN"/
  DATA KW(265)%TYP /"L:E"/
  DATA KW(265)%DSCRPT /"*! Performs atom-centered Lowdin symmetric &
       &orthogonalization in generating the NAOs. !*"/
           
  DATA KW(266)%LABEL /"NBO_WRITE_LCLOWDIN"/
  DATA KW(266)%TYP /"L:E"/
  DATA KW(266)%DSCRPT /"*! Write a GENNBO FILE.47 containing all the atoms &
       &in the atom-centered Lowdin basis to satisfy the strict &
       &lm-orthogonality requirement in GENNBO !*"/
           
  DATA KW(267)%LABEL /"NBO_WRITE_NPACOMP"/
  DATA KW(267)%TYP /"L:B"/
  DATA KW(267)%DSCRPT /"*! Writes individual NAO population into the &
       &standard output. !*"/
           
  DATA KW(268)%LABEL /"NBO_SCALE_DM"/
  DATA KW(268)%TYP /"L:E"/
  DATA KW(268)%DSCRPT /"*! Scales density matrix in the FILE.47 output to &
       &achieve charge integrality (Required for proper GENNBO &
       &functionality). !*"/

  DATA KW(269)%LABEL /"NLPP_FOR_EXCHANGE"/
  DATA KW(269)%TYP /"L:E"/
  DATA KW(269)%DSCRPT /"*! Give exchange matrix same sparsity as non-local pseudopotential matrix !*"/

  DATA KW(270)%LABEL /"RADIAL_SEGMENTS"/
  DATA KW(270)%TYP /"I:E"/
  DATA KW(270)%DSCRPT /"*! Number of radial segments in vmatrix evaluation !*"/

  DATA KW(271)%LABEL /"ANGULAR_SEGMENTS"/
  DATA KW(271)%TYP /"I:E"/
  DATA KW(271)%DSCRPT /"*! Number of angular segments in vmatrix evaluation !*"/

  DATA KW(272)%LABEL /"HFX_INTEGRATION"/
  DATA KW(272)%TYP /"I:E"/
  DATA KW(272)%DSCRPT /"*! Variant of the numerical pointwise approach in HF exchange !*"/

  DATA KW(273)%LABEL /"MD_PROPERTIES"/
  DATA KW(273)%TYP /"L:I"/
  DATA KW(273)%DSCRPT/"*! Compute vibrational and IR spectra from MD !*"/ 

  DATA KW(274)%LABEL /"MD_RESTART"/
  DATA KW(274)%TYP /"L:I"/
  DATA KW(274)%DSCRPT/"*! Restart MD from backup files !*"/ 

  DATA KW(275)%LABEL /"THERMOSTAT"/
  DATA KW(275)%TYP /"B:B"/
  DATA KW(275)%DSCRPT/"*! Thermostat for MD in NVT ensemble !*"/ 

  DATA KW(276)%LABEL /"MD_RESET_DENSKERN_NGWFS"/
  DATA KW(276)%TYP /"I:I"/
  DATA KW(276)%DSCRPT/"*! Reset mixing scheme for initial guess of NGWFs and density kernel !*"/ 

  DATA KW(277)%LABEL /"NBO_AOPNAO_SCHEME"/
  DATA KW(277)%TYP /"T:E"/
  DATA KW(277)%DSCRPT /"*! AO to PNAO scheme to use in generating NAOs (for testing purposes). !*"/

  DATA KW(278)%LABEL /"NBO_LIST_PLOTNBO"/
  DATA KW(278)%TYP /"B:I"/
  DATA KW(278)%DSCRPT /"*! List of NBOs to be plotted according to GENNBO &
       &output indices. !*"/

  DATA KW(279)%LABEL /"PLOT_NBO"/
  DATA KW(279)%TYP /"L:I"/
  DATA KW(279)%DSCRPT /"*! Plot NBO's orbitals from FILE.xx as defined by nbo_plot_orbtype. !*"/

  DATA KW(280)%LABEL /"NBO_PLOT_ORBTYPE"/
  DATA KW(280)%TYP /"T:I"/
  DATA KW(280)%DSCRPT /"*! Type of GENNBO-generated orbital to plot. !*"/

  DATA KW(281)%LABEL /"NBO_PNAO_ANALYSIS"/
  DATA KW(281)%TYP /"L:I"/
  DATA KW(281)%DSCRPT /"*! S/P/D/F CHARACTER ANALYSIS ON PNAO. !*"/

  DATA KW(282)%LABEL /"PHONON_FINITE_DISP"/
  DATA KW(282)%TYP /"P:E"/
  DATA KW(282)%DSCRPT /"*! Amplitude of the ionic perturbation to be used in a finite displacement phonon calculation !*"/

  DATA KW(283)%LABEL /"PHONON_FMAX"/
  DATA KW(283)%TYP /"P:E"/
  DATA KW(283)%DSCRPT /"*! Maximum force allowed on the unperturbed system for a phonon calculation !*"/

  DATA KW(284)%LABEL /"PHONON_FARMING_TASK"/
  DATA KW(284)%TYP /"I:I"/
  DATA KW(284)%DSCRPT /"*! Operation to perform for phonon calc (for task farming or post-proc. of dynamical matrix) !*"/

  DATA KW(285)%LABEL /"PHONON_DISP_LIST"/
  DATA KW(285)%TYP /"B:I"/
  DATA KW(285)%DSCRPT /"*! List of displacements to perform for phonon calculation !*"/

  DATA KW(286)%LABEL /"PHONON_TMIN"/
  DATA KW(286)%TYP /"P:B"/
  DATA KW(286)%DSCRPT /"*! Lower bound of temperature range for computation of vibrational thermodynamic quantities !*"/

  DATA KW(287)%LABEL /"PHONON_TMAX"/
  DATA KW(287)%TYP /"P:B"/
  DATA KW(287)%DSCRPT /"*! Upper bound of temperature range for computation of vibrational thermodynamic quantities !*"/

  DATA KW(288)%LABEL /"PHONON_DELTAT"/
  DATA KW(288)%TYP /"P:B"/
  DATA KW(288)%DSCRPT /"*! Temperature step for computation of vibrational thermodynamic quantities !*"/

  DATA KW(289)%LABEL /"PHONON_MIN_FREQ"/
  DATA KW(289)%TYP /"P:E"/
  DATA KW(289)%DSCRPT /"*! Discard phonon frequencies smaller than this value for computation of vibrational &
       &thermodynamic quantities !*"/

  DATA KW(290)%LABEL /"GEOM_LBFGS"/
  DATA KW(290)%TYP /"L:B"/
  DATA KW(290)%DSCRPT /"*! Whether to perform LBFGS rather than BFGS in a Geometry Optimization !*"/

  DATA KW(291)%LABEL /"GEOM_LBFGS_MAX_UPDATES"/
  DATA KW(291)%TYP /"I:I"/
  DATA KW(291)%DSCRPT /"*! Number of LBFGS update vectors to store !*"/

  DATA KW(292)%LABEL /"GEOM_LBFGS_BLOCK_LENGTH"/
  DATA KW(292)%TYP /"I:E"/
  DATA KW(292)%DSCRPT /"*! How many updates to store before reallocation in an unbounded LBFGS calculation !*"/

  DATA KW(293)%LABEL /"SPECIES_CORE_WF"/
  DATA KW(293)%TYP /"B:I"/
  DATA KW(293)%DSCRPT /"*! Core Wavefunction filename for each species !*"/

  DATA KW(294)%LABEL /"COND_SPEC_PRINT_MAT_ELS"/
  DATA KW(294)%TYP /"L:B"/
  DATA KW(294)%DSCRPT /"*! Write optical matrix elements to file !*"/ 

  DATA KW(295)%LABEL /"COND_SPEC_SCISS_OP"/
  DATA KW(295)%TYP /"P:I"/
  DATA KW(295)%DSCRPT /"*! Scissor operator for JDOS and imag. diel. fn. !*"/ 

  DATA KW(296)%LABEL /"COND_SPEC_OPT_SMEAR"/
  DATA KW(296)%TYP /"P:I"/
  DATA KW(296)%DSCRPT /"*! Half width of smearing Gaussians for JDOS and imag. diel. fn. !*"/ 

  DATA KW(297)%LABEL /"COND_CALC_EELS"/
  DATA KW(297)%TYP /"L:B"/
  DATA KW(297)%DSCRPT /"*! Calculate matrix elements for electron energy loss spectra (EELS) !*"/ ! lr408

  DATA KW(298)%LABEL /"EXTERNAL_PRESSURE"/
  DATA KW(298)%TYP /"P:B"/
  DATA KW(298)%DSCRPT/"*! External applied pressure !*"/

  DATA KW(299)%LABEL /"SMOOTHING_FACTOR"/
  DATA KW(299)%TYP /"D:B"/
  DATA KW(299)%DSCRPT/"*! Smoothing factor in volume term !*"/

  DATA KW(300)%LABEL /"ISOSURFACE_CUTOFF"/
  DATA KW(300)%TYP /"D:B"/
  DATA KW(300)%DSCRPT/"*! Isosurface cutoff in volume term !*"/

  DATA KW(301)%LABEL /"CONSTRAINED_DFT"/
  DATA KW(301)%TYP /"B:I"/
  DATA KW(301)%DSCRPT /"*! Constrained_DFT species info [1:symb., 2:l=angul. mom., 3. Ionic_charge, 4: U_occ. (eV), 5:U_q_up (eV), &
           &6:U_q_down (eV), 7:U_spin (eV), 8: Targeted N_up, 9: Targeted N_down, 10:targeted N_up - N_down !*"/

  DATA KW(302)%LABEL /"CI_CDFT"/
  DATA KW(302)%TYP /"L:I"/
  DATA KW(302)%DSCRPT /"*! Perform a CONFIGURATION-INTERACTION CDFT simulation !*"/

  DATA KW(303)%LABEL /"CDFT_ATOM_CHARGE"/
  DATA KW(303)%TYP /"L:I"/
  DATA KW(303)%DSCRPT /"*! Perform an ATOM-CHARGE-constrained CDFT simulation !*"/

  DATA KW(304)%LABEL /"CDFT_ATOM_SPIN"/
  DATA KW(304)%TYP /"L:I"/
  DATA KW(304)%DSCRPT /"*! Perform an ATOM-SPIN-constrained CDFT simulation !*"/

  DATA KW(305)%LABEL /"CDFT_GROUP_CHARGE"/
  DATA KW(305)%TYP /"L:I"/
  DATA KW(305)%DSCRPT /"*! Perform a GROUP-CHARGE-constrained CDFT simulation !*"/

  DATA KW(306)%LABEL /"CDFT_GROUP_SPIN"/
  DATA KW(306)%TYP /"L:I"/
  DATA KW(306)%DSCRPT /"*! Perform a GROUP-SPIN-constrained CDFT simulation !*"/

  DATA KW(307)%LABEL /"CDFT_GROUP_CHARGE_DIFF"/
  DATA KW(307)%TYP /"L:I"/
  DATA KW(307)%DSCRPT /"*! Perform a GROUP-CHARGE-DIFFERENCE-constrained CDFT simulation !*"/

  DATA KW(308)%LABEL /"CDFT_GROUP_SPIN_DIFF"/
  DATA KW(308)%TYP /"L:I"/
  DATA KW(308)%DSCRPT /"*! Perform a GROUP-SPIN-DIFFERENCE-constrained CDFT simulation !*"/

  DATA KW(309)%LABEL /"CDFT_HUBBARD"/
  DATA KW(309)%TYP /"L:I"/
  DATA KW(309)%DSCRPT /"*! Perform a constrained-DFT+U simulation !*"/

  DATA KW(310)%LABEL /"CDFT_GROUP_CHARGE_TARGET"/
  DATA KW(310)%TYP /"D:I"/
  DATA KW(310)%DSCRPT /"*! Targeted group-CHARGE for GROUP-CHARGE-constrained cDFT !*"/

  DATA KW(311)%LABEL /"CDFT_GROUP_SPIN_TARGET"/
  DATA KW(311)%TYP /"D:I"/
  DATA KW(311)%DSCRPT /"*! Targeted group-SPIN for GROUP-SPIN-constrained cDFT !*"/

  DATA KW(312)%LABEL /"CDFT_GROUP_CHARGE_DIFF_TARGET"/
  DATA KW(312)%TYP /"D:I"/
  DATA KW(312)%DSCRPT /"*! Targeted CHARGE difference (acceptor-donor) for  GROUP-CHARGE-DIFFERENCE-constrained cDFT !*"/

  DATA KW(313)%LABEL /"CDFT_GROUP_SPIN_DIFF_TARGET"/
  DATA KW(313)%TYP /"D:I"/
  DATA KW(313)%DSCRPT /"*! Targeted SPIN difference (acceptor-donor) for  GROUP-SPIN-DIFFERENCE-constrained cDFT !*"/

  DATA KW(314)%LABEL /"CDFT_GROUP_U"/
  DATA KW(314)%TYP /"P:I"/
  DATA KW(314)%DSCRPT /"*! Constraining potential (eV) for GROUP-CHARGE/SPIN-constrained CDFT simulation !*"/

  DATA KW(315)%LABEL /"CDFT_GROUP_DIFF_U"/
  DATA KW(315)%TYP /"P:I"/
  DATA KW(315)%DSCRPT /"*! Constraining potential (eV) for GROUP-CHARGE/SPIN-DIFFERENCE-constrained CDFT simulation !*"/

  DATA KW(316)%LABEL /"MAXIT_CDFT_U_CG"/
  DATA KW(316)%TYP /"I:I"/
  DATA KW(316)%DSCRPT /"*! Max number of cdFT-U conjugate gradients (CG) &
       &iterations !*"/

  DATA KW(317)%LABEL /"CDFT_CG_TYPE"/
  DATA KW(317)%TYP /"T:E"/
  DATA KW(317)%DSCRPT /"*! Type of CG coefficient for CDFT U-optimisation !*&
      & NGWF_POLAK = Polak-Ribbiere formula; &
      & NGWF_FLETCHER = Fletcher-Reeves formula."/

  DATA KW(318)%LABEL /"CDFT_CG_THRESHOLD"/
  DATA KW(318)%TYP /"D:I"/
  DATA KW(318)%DSCRPT /"*! RMS gradient convergence threshold for U-pot. in CDFT !*"/

  DATA KW(319)%LABEL /"CDFT_CG_MAX"/
  DATA KW(319)%TYP /"I:E"/
  DATA KW(319)%DSCRPT /"*! Number of U-opt iterations to reset CG !*"/

  DATA KW(320)%LABEL /"CDFT_MAX_GRAD"/
  DATA KW(320)%TYP /"D:I"/
  DATA KW(320)%DSCRPT /"*! Maximum permissible value of CDFT U-Gradient for convergence !*"/

  DATA KW(321)%LABEL /"CDFT_ELEC_ENERGY_TOL"/
  DATA KW(321)%TYP /"P:I"/
  DATA KW(321)%DSCRPT /"*! Tolerance on total energy change during CDFT optimisation !*"/

  DATA KW(322)%LABEL /"CDFT_CG_MAX_STEP"/
  DATA KW(322)%TYP /"D:E"/
  DATA KW(322)%DSCRPT /"*! Maximum length of trial step for cDFT optimisation line search !*"/

  DATA KW(323)%LABEL /"CDFT_GURU"/
  DATA KW(323)%TYP /"L:I"/
  DATA KW(323)%DSCRPT /"*! Let the user signal she/he does not need helpt with the cDFT U-initialisation !*"/

  DATA KW(324)%LABEL /"CDFT_CONTINUATION"/
  DATA KW(324)%TYP /"L:I"/
  DATA KW(324)%DSCRPT /"*! Logical to restart (from the *.cdft file) a cDFT U-optimisation !*"/

  DATA KW(325)%LABEL /"ETRANS_CALCULATE"/
  DATA KW(325)%TYP /"L:I"/
  DATA KW(325)%DSCRPT /"*! Compute electronic transport coefficients !*"/
             
  DATA KW(326)%LABEL /"ETRANS_SAME_LEADS"/
  DATA KW(326)%TYP /"L:I"/
  DATA KW(326)%DSCRPT /"*! Use same description for all leads !*"/
             
  DATA KW(327)%LABEL /"ETRANS_BULK"/
  DATA KW(327)%TYP /"L:I"/
  DATA KW(327)%DSCRPT /"*! Compute leads self-energies on-the-fly !*"/
             
  DATA KW(328)%LABEL /"ETRANS_ECMPLX"/
  DATA KW(328)%TYP /"P:I"/
  DATA KW(328)%DSCRPT /"*! Imaginary part of the energy !*"/
             
  DATA KW(329)%LABEL /"ETRANS_SETUP"/
  DATA KW(329)%TYP /"B:I"/
  DATA KW(329)%DSCRPT /"*! Transport setup description !*"/
             
  DATA KW(330)%LABEL /"ETRANS_EMAX"/
  DATA KW(330)%TYP /"P:I"/
  DATA KW(330)%DSCRPT /"*! Highest energy for the calculation of transmission coefficients !*"/
              
  DATA KW(331)%LABEL /"ETRANS_EMIN"/
  DATA KW(331)%TYP /"P:I"/
  DATA KW(331)%DSCRPT /"*! Lowest energy for the calculation of transmission coefficients !*"/
             
  DATA KW(332)%LABEL /"ETRANS_ENUM"/
  DATA KW(332)%TYP /"I:I"/
  DATA KW(332)%DSCRPT /"*! Number of energy points for the calculation of transmission coefficients !*"/
              
  DATA KW(333)%LABEL /"ETRANS_SOURCE"/
  DATA KW(333)%TYP /"I:I"/
  DATA KW(333)%DSCRPT /"*! Index of the lead used as source !*"/
            
  DATA KW(334)%LABEL /"ETRANS_DRAIN"/
  DATA KW(334)%TYP /"I:I"/
  DATA KW(334)%DSCRPT /"*! Index of the lead used as drain !*"/
            
  DATA KW(335)%LABEL /"ETRANS_WRITE_SETUP"/
  DATA KW(335)%TYP /"L:I"/
  DATA KW(335)%DSCRPT /"*! Write hamiltonian corresponding to transport setup !*"/

  DATA KW(336)%LABEL /"SPECIES_ATOMIC_SET_COND"/
  DATA KW(336)%TYP /"B:I"/
  DATA KW(336)%DSCRPT /"*! Atomic set description for each species, for initialising Conduction NGWFs !*"/

  DATA KW(337)%LABEL /"SPECIES_ATOMIC_SET_AUX"/
  DATA KW(337)%TYP /"B:I"/
  DATA KW(337)%DSCRPT /"*! Atomic set description for each species, for initialising Auxiliary NGWFs !*"/

  DATA KW(338)%LABEL /"SPECIES_AUX"/
  DATA KW(338)%TYP /"B:B"/
  DATA KW(338)%DSCRPT /"*! Species information for Auxiliary NGWFs (symbol, atomic number, number of NGWFs, NGWF radius) !*"/

  DATA KW(339)%LABEL /"HFX_NLPP_FOR_EXCHANGE"/
  DATA KW(339)%TYP /"L:E"/
  DATA KW(339)%DSCRPT /"*! Give exchange matrix same sparsity as non-local pseudopotential matrix !*"/

  DATA KW(340)%LABEL /"HFX_RADIAL_SEGMENTS"/
  DATA KW(340)%TYP /"I:E"/
  DATA KW(340)%DSCRPT /"*! Number of radial segments in vmatrix evaluation !*"/

  DATA KW(341)%LABEL /"HFX_ANGULAR_SEGMENTS"/
  DATA KW(341)%TYP /"I:E"/
  DATA KW(341)%DSCRPT /"*! Number of angular segments in vmatrix evaluation !*"/

  DATA KW(342)%LABEL /"HFX_CHEB_ORDER"/
  DATA KW(342)%TYP /"I:E"/
  DATA KW(342)%DSCRPT /"*! Order for Chebyshev interpolation !*"/

  DATA KW(343)%LABEL /"HFX_CHEB_INTERVALS"/
  DATA KW(343)%TYP /"I:E"/
  DATA KW(343)%DSCRPT /"*! Number of intervals Chebyshev interpolation !*"/

  DATA KW(344)%LABEL /"HFX_CHEB_A_BATCHSIZE"/
  DATA KW(344)%TYP /"I:E"/
  DATA KW(344)%DSCRPT /"*! Number of SW expansions buffered !*"/

  DATA KW(345)%LABEL /"HFX_CHEB_B_BATCHSIZE"/
  DATA KW(345)%TYP /"I:E"/
  DATA KW(345)%DSCRPT /"*! Number of SW pot expansions buffered !*"/

  DATA KW(346)%LABEL /"HFX_READ_VMATRIX"/
  DATA KW(346)%TYP /"L:B"/
  DATA KW(346)%DSCRPT /"*! Read restart information for the V matrix !*"/

  DATA KW(347)%LABEL /"HFX_WRITE_VMATRIX"/
  DATA KW(347)%TYP /"L:B"/
  DATA KW(347)%DSCRPT /"*! Write restart information for the V matrix !*"/

  DATA KW(348)%LABEL /"HFX_MAX_L"/
  DATA KW(348)%TYP /"I:I"/
  DATA KW(348)%DSCRPT /"*! Maximum quantum number l used for SW expansion !*"/

  DATA KW(349)%LABEL /"HFX_MAX_ZEROS"/
  DATA KW(349)%TYP /"I:I"/
  DATA KW(349)%DSCRPT /"*! Maximum number of Bessel zeros used for SW expansion !*"/

  DATA KW(350)%LABEL /"HFX_DEBUG"/
  DATA KW(350)%TYP /"L:E"/
  DATA KW(350)%DSCRPT /"*! Debug HFx? !*"/

  DATA KW(351)%LABEL /"HFX_READ_XMATRIX"/
  DATA KW(351)%TYP /"L:E"/
  DATA KW(351)%DSCRPT /"*! Read restart information for the X matrix !*"/

  DATA KW(352)%LABEL /"HFX_WRITE_XMATRIX"/
  DATA KW(352)%TYP /"L:E"/
  DATA KW(352)%DSCRPT /"*! Write restart information for the X matrix !*"/

!CW
  DATA KW(i+1)%LABEL /"DMFT_POINTS"/
  DATA KW(i+1)%TYP /"I:I"/
  DATA KW(i+1)%DSCRPT /"*! Number of DMFT energy points on real or matsubara axis !*"/
  
  DATA KW(i+2)%LABEL /"DMFT_TEMP"/
  DATA KW(i+2)%TYP /"P:I"/
  DATA KW(i+2)%DSCRPT /"*! Temperature (in Hartree) for DFT+DMFT!*"/

  DATA KW(i+3)%LABEL /"DMFT_EMIN"/
  DATA KW(i+3)%TYP /"P:I"/
  DATA KW(i+3)%DSCRPT /"*! Min energy on real axis (in Hartree) for DFT+DMFT!*"/

  DATA KW(i+4)%LABEL /"DMFT_EMAX"/
  DATA KW(i+4)%TYP /"P:I"/
  DATA KW(i+4)%DSCRPT /"*! Max energy (in Hartree) for DFT+DMFT!*"/

  DATA KW(i+5)%LABEL /"DMFT_CHEM_SHIFT"/
  DATA KW(i+5)%TYP /"P:I"/
  DATA KW(i+5)%DSCRPT /"*! Shift of chemical potential (chemical doping, in Hartree) for DFT+DMFT!*"/

  DATA KW(i+6)%LABEL /"DMFT_SMEAR_SHIFT"/
  DATA KW(i+6)%TYP /"P:I"/
  DATA KW(i+6)%DSCRPT /"*! Frequency dependent Smearing (in Hartree) for DFT+DMFT!*"/

  DATA KW(i+7)%LABEL /"DMFT_SMEAR"/
  DATA KW(i+7)%TYP /"P:I"/
  DATA KW(i+7)%DSCRPT /"*! Smearing (in Hartree) for DFT+DMFT!*"/

  DATA KW(i+8)%LABEL /"DMFT_PARAMAGNETIC"/
  DATA KW(i+8)%TYP /"L:I"/
  DATA KW(i+8)%DSCRPT /"*! Imposes the paramagnetic state in DFT+DMFT!*"/

  DATA KW(i+9)%LABEL /"DMFT_CUTOFF_TAIL"/
  DATA KW(i+9)%TYP /"P:I"/
  DATA KW(i+9)%DSCRPT /"*!ATOMIC ENERGY THRESHOLD FOR DOUBLE/SINGLE PREC FOR THE GREEN-FUNCTION !*"/

  DATA KW(i+10)%LABEL /"DMFT_ROTATE_GREEN"/
  DATA KW(i+10)%TYP /"L:I"/
  DATA KW(i+10)%DSCRPT /"*!IF TRUE IT WILL WRITE OUT THE GREEN FUNCTION IN THE ROTATED LOCAL ATOMIC BASIS !*"/

  DATA KW(i+11)%LABEL /"DMFT_CUTOFF_SMALL"/
  DATA KW(i+11)%TYP /"P:I"/
  DATA KW(i+11)%DSCRPT /"*!ATOMIC ENERGY THRESHOLD FOR DOUBLE/SINGLE PREC FOR THE GREEN-FUNCTION !*"/

  DATA KW(i+12)%LABEL /"DMFT_FREE_GREEN_FREQU"/
  DATA KW(i+12)%TYP /"P:I"/
  DATA KW(i+12)%DSCRPT /"*!ATOMIC ENERGY THRESHOLD FOR THE FREE GREEN-FUNCTION !*"/

  DATA KW(i+13)%LABEL /"DMFT_PLOT_REAL_SPACE"/
  DATA KW(i+13)%TYP /"L:I"/
  DATA KW(i+13)%DSCRPT /"*!IF TRUE IT WILL WRITE OUT THE REAL SPACE DMFT QUANTITIES !*"/

  DATA KW(i+14)%LABEL /"DMFT_OPTICS"/
  DATA KW(i+14)%TYP /"L:I"/
  DATA KW(i+14)%DSCRPT /"*!IF TRUE IT WILL COMPUTE THE OPTICAL CONDUCTIVITY FROM THE DMFT GREEN FUNCTION !*"/

  DATA KW(i+15)%LABEL /"DMFT_OPTICS_WINDOW"/
  DATA KW(i+15)%TYP /"P:I"/
  DATA KW(i+15)%DSCRPT /"*!WINDOW OF ENERGY AROUND FERMI ENERGY FOR OPTICAL CONDUCTIVITY !*"/

  DATA KW(i+16)%LABEL /"DMFT_OPTICS_I1"/
  DATA KW(i+16)%TYP /"I:I"/
  DATA KW(i+16)%DSCRPT /"*! first direction for optical conductivity current-current correlator !*"/

  DATA KW(i+17)%LABEL /"DMFT_OPTICS_I2"/
  DATA KW(i+17)%TYP /"I:I"/
  DATA KW(i+17)%DSCRPT /"*! second direction for optical conductivity current-current correlator !*"/

  DATA KW(i+18)%LABEL /"DMFT_DOS_MIN"/
  DATA KW(i+18)%TYP /"P:I"/
  DATA KW(i+18)%DSCRPT /"*! Min of energy window on real axis (in Hartree) for DOS calculations !*"/

  DATA KW(i+19)%LABEL /"DMFT_DOS_MAX"/
  DATA KW(i+19)%TYP /"P:I"/
  DATA KW(i+19)%DSCRPT /"*! Max of energy window on real axis (in Hartree) for DOS calculations !*"/

  DATA KW(i+20)%LABEL /"DMFT_SMEAR_T"/
  DATA KW(i+20)%TYP /"P:I"/
  DATA KW(i+20)%DSCRPT /"*! Frequency dependent Smearing (in Hartree) for DFT+DMFT!*"/

  DATA KW(i+21)%LABEL /"DMFT_SMEAR_W"/
  DATA KW(i+21)%TYP /"P:I"/
  DATA KW(i+21)%DSCRPT /"*! Frequency dependent Smearing (in Hartree) for DFT+DMFT!*"/

  DATA KW(i+22)%LABEL /"DMFT_SMEAR_ETA"/
  DATA KW(i+22)%TYP /"P:I"/
  DATA KW(i+22)%DSCRPT /"*! Frequency dependent Smearing (in Hartree) for DFT+DMFT!*"/

  DATA KW(i+23)%LABEL /"DMFT_PLOT_REAL_SPACE_SIGMA"/
  DATA KW(i+23)%TYP /"L:I"/
  DATA KW(i+23)%DSCRPT /"*!IF TRUE IT WILL ONLY WRITE OUT THE REAL SPACE DMFT QUANTITIES FOR SIGMA !*"/

  DATA KW(i+24)%LABEL /"DMFT_NORM_PROJ"/
  DATA KW(i+24)%TYP /"I:I"/
  DATA KW(i+24)%DSCRPT /"*!IF =1,2 IT WILL ORTHONORMALIZE THE HUBBARD PROJECTORS U^+ U , IF =3 U U^+FOR THE DMFT !*"/

  DATA KW(i+25)%LABEL /"DMFT_PLOT_ALL_PROJ"/
  DATA KW(i+25)%TYP /"L:I"/
  DATA KW(i+25)%DSCRPT /"*! IF TRUE WILL PLOT ALL DIAOGNAL CORRELATED ELEMENTS EVEN THE ONE NOT INCLUDED IN THE CLUSTER CALCULATION (MASK_DIMER) !*"/

  DATA KW(i+26)%LABEL /"DMFT_INTEGRATE_GREEN"/
  DATA KW(i+26)%TYP /"L:I"/
  DATA KW(i+26)%DSCRPT /"*! IF TRUE WILL PLOT THE INTEGRATED DENSITY (both in real space and per orbital) !*"/

  DATA KW(i+27)%LABEL /"DMFT_OPTICS_X1"/
  DATA KW(i+27)%TYP /"D:B"/
  DATA KW(i+27)%DSCRPT /"*! UNITARY VECTOR (ex componant) WHICH DEFINES THE OPTICAL DIRECTION !*"/

  DATA KW(i+28)%LABEL /"DMFT_OPTICS_Y1"/
  DATA KW(i+28)%TYP /"D:B"/
  DATA KW(i+28)%DSCRPT /"*! UNITARY VECTOR (ey componant) WHICH DEFINES THE OPTICAL DIRECTION !*"/

  DATA KW(i+29)%LABEL /"DMFT_OPTICS_Z1"/
  DATA KW(i+29)%TYP /"D:B"/
  DATA KW(i+29)%DSCRPT /"*! UNITARY VECTOR (ez componant) WHICH DEFINES THE OPTICAL DIRECTION !*"/

  DATA KW(i+30)%LABEL /"DMFT_E_PLOT"/
  DATA KW(i+30)%TYP /"P:I"/
  DATA KW(i+30)%DSCRPT /"*! Energy level to plot the real space upfolded functions from the DMFT !*"/

  DATA KW(i+31)%LABEL /"DMFT_IN_BOHR"/
  DATA KW(i+31)%TYP /"L:I"/
  DATA KW(i+31)%DSCRPT /"*! Should be true if the rotations-dimer mask files in dmft are in bohr notations !*"/

  DATA KW(i+32)%LABEL /"DMFT_ORDER_PROJ"/
  DATA KW(i+32)%TYP /"D:I"/
  DATA KW(i+32)%DSCRPT /"*! small number to enforce the orbital order of the NGWFS when the overlap with projectors are identical !*"/

  DATA KW(i+33)%LABEL /"DMFT_NMU_LOOP"/
  DATA KW(i+33)%TYP /"I:I"/
  DATA KW(i+33)%DSCRPT /"*! number of iterations for the Newtons method used to find the chemical potential!*"/

  DATA KW(i+34)%LABEL /"DMFT_SC"/
  DATA KW(i+34)%TYP /"L:I"/
  DATA KW(i+34)%DSCRPT /"*! if true will run the self consistent ONETEP+DMFT, if false, one shot ONETEP+DMFT is used*"/

  DATA KW(i+35)%LABEL /"DMFT_KERNEL"/
  DATA KW(i+35)%TYP /"I:I"/
  DATA KW(i+35)%DSCRPT /"*! 1:writes dmft density kernel, 0:does not calculate it, -1:writes the purified dmft density kernel !*"/

  DATA KW(i+36)%LABEL /"DMFT_SPOIL_KERNEL"/
  DATA KW(i+36)%TYP /"L:I"/
  DATA KW(i+36)%DSCRPT /"*! if false will NOT update LHXC potential*"/

  DATA KW(i+37)%LABEL /"DMFT_FULLY_SC"/
  DATA KW(i+37)%TYP /"L:I"/
  DATA KW(i+37)%DSCRPT /"*! if true uses the self energy in the ONETEP kernel/NGWF optimization, in the energy functional*"/

  DATA KW(i+38)%LABEL /"DMFT_KERNEL_MIX"/
  DATA KW(i+38)%TYP /"D:I"/
  DATA KW(i+38)%DSCRPT /"*! mixing of DMFT and DFT kernel for DFT+DMFT*"/

  DATA KW(i+39)%LABEL /"DMFT_PURIFY_SC"/
  DATA KW(i+39)%TYP /"L:I"/
  DATA KW(i+39)%DSCRPT /"*! if true uses the purified kernel for DFT_DMFT in the property module*"/

  DATA KW(i+40)%LABEL /"DMFT_KS_SHIFT"/
  DATA KW(i+40)%TYP /"L:I"/
  DATA KW(i+40)%DSCRPT /"*! if true adds the correction to the energy in the SC dmft minimization during kernel optimization where the re-occupation of the energy level by the DMFT density kernel is taken into account*"/

  DATA KW(i+41)%LABEL /"DMFT_FULLY_SC_H"/
  DATA KW(i+41)%TYP /"L:I"/
  DATA KW(i+41)%DSCRPT /"*! if true the H used in the DFT is the KS H built from 1 shot DMFT*"/

  DATA KW(i+42)%LABEL /"DMFT_SWITCH_OFF_PROJ_ORDER"/
  DATA KW(i+42)%TYP /"L:I"/
  DATA KW(i+42)%DSCRPT /"*! for compatibility with older version of ONETEP-DMFT, if true switches off the natural order of the projections*"/

  DATA KW(i+43)%LABEL /"DMFT_MU_DIFF_MAX"/
  DATA KW(i+43)%TYP /"D:I"/
  DATA KW(i+43)%DSCRPT /"*! max diff from target when imposing DMFT chemical potential*"/

  DATA KW(i+44)%LABEL /"DMFT_EMBED_ITER"/
  DATA KW(i+44)%TYP /"I:I"/
  DATA KW(i+44)%DSCRPT /"*! number of iterations for the embedding procedure within DMFT !*"/

  DATA KW(i+45)%LABEL /"DMFT_EMBED_MIX"/
  DATA KW(i+45)%TYP /"D:I"/
  DATA KW(i+45)%DSCRPT /"*! mixing of internal convergence for the embedding procedure within DMFT !*"/

  DATA KW(i+46)%LABEL /"DMFT_NKPOINTS"/
  DATA KW(i+46)%TYP /"I:I"/
  DATA KW(i+46)%DSCRPT /"*! number of K points for averaging the lattice Green Function !*"/

  DATA KW(i+47)%LABEL /"DMFT_MU_ORDER"/
  DATA KW(i+47)%TYP /"I:I"/
  DATA KW(i+47)%DSCRPT /"*! order of the Newton Method used to fix the chemical potential (HouseHolder general form) !*"/

  DATA KW(i+48)%LABEL /"DMFT_SKIP_ENERGY"/
  DATA KW(i+48)%TYP /"L:I"/
  DATA KW(i+48)%DSCRPT /"*! if true the energy is not computed along the Newton steepest descent used to find the chemical potential*"/

  DATA KW(i+49)%LABEL /"DMFT_KPOINTS_SYM"/
  DATA KW(i+49)%TYP /"L:I"/
  DATA KW(i+49)%DSCRPT /"*! if true and using k points, additional cubic symmetry are used to reduce the number of k points, and not only the k inversion symmetry. Note: this is NOT valid when computing the DMFT density kernel*"/

  DATA KW(i+50)%LABEL /"DMFT_KPOINTS_KERNEL_GAMMA"/
  DATA KW(i+50)%TYP /"L:I"/
  DATA KW(i+50)%DSCRPT /"*! if true and using k points, the kpoints are used to fix the chemical potential, but only the Gamma point is used to output the dmft density kernel*"/

  DATA KW(i+51)%LABEL /"DMFT_LIN_SCALING"/
  DATA KW(i+51)%TYP /"L:I"/
  DATA KW(i+51)%DSCRPT /"*! uses a linear scaling inversion scheme for the Green function - as implemented, this provides a good estimate of the projected GF required for the DMFT, but the kernel would be badly approximated, so it cannot be used as is when dmft_kernel/=0*"/

  DATA KW(i+52)%LABEL /"DMFT_WIN"/
  DATA KW(i+52)%TYP /"D:I"/
  DATA KW(i+52)%DSCRPT /"*!window of energies around the chemical potential considered for the GF inversion *"/

  DATA KW(i+53)%LABEL /"DMFT_NVAL"/
  DATA KW(i+53)%TYP /"I:I"/
  DATA KW(i+53)%DSCRPT /"*!max number of eigenstates in this energy window used for the GF inversion *"/

  DATA KW(i+54)%LABEL /"DMFT_SCALING_ITER"/
  DATA KW(i+54)%TYP /"I:I"/
  DATA KW(i+54)%DSCRPT /"*! *"/

  DATA KW(i+55)%LABEL /"DMFT_SCALING_MAXSPACE"/
  DATA KW(i+55)%TYP /"I:I"/
  DATA KW(i+55)%DSCRPT /"*! *"/

  DATA KW(i+56)%LABEL /"DMFT_SCALING_TOL"/
  DATA KW(i+56)%TYP /"D:I"/
  DATA KW(i+56)%DSCRPT /"*! *"/

  DATA KW(i+57)%LABEL /"DMFT_SCALING_CUTOFF"/
  DATA KW(i+57)%TYP /"D:I"/
  DATA KW(i+57)%DSCRPT /"*! *"/

  DATA KW(i+58)%LABEL /"DMFT_SCALING_NMPI"/
  DATA KW(i+58)%TYP /"I:I"/
  DATA KW(i+58)%DSCRPT /"*! *"/

  DATA KW(i+59)%LABEL /"DMFT_SCALING_METH"/
  DATA KW(i+59)%TYP /"I:I"/
  DATA KW(i+59)%DSCRPT /"*! *"/
 
  DATA KW(i+60)%LABEL /"DMFT_SCALING_TAIL"/
  DATA KW(i+60)%TYP /"D:I"/
  DATA KW(i+60)%DSCRPT /"*! *"/ 

  DATA KW(i+61)%LABEL /"DMFT_SPLIT"/
  DATA KW(i+61)%TYP /"L:I"/
  DATA KW(i+61)%DSCRPT /"*!  *"/

  DATA KW(i+62)%LABEL /"DMFT_INVERT_OVERLAP"/
  DATA KW(i+62)%TYP /"L:I"/
  DATA KW(i+62)%DSCRPT /"*!  *"/

  DATA KW(i+63)%LABEL /"DMFT_LOCAL_SCRATCH"/
  DATA KW(i+63)%TYP /"L:I"/
  DATA KW(i+63)%DSCRPT /"*!  *"/

  DATA KW(i+64)%LABEL /"DMFT_SPLITK"/
  DATA KW(i+64)%TYP /"L:I"/
  DATA KW(i+64)%DSCRPT /"*!  *"/

  DATA KW(i+65)%LABEL /"DMFT_SCALING_CUTOFF_H"/
  DATA KW(i+65)%TYP /"D:I"/
  DATA KW(i+65)%DSCRPT /"*! *"/

  DATA KW(i+66)%LABEL /"WRITE_POLARISATION_PLOT"/
  DATA KW(i+66)%TYP /"L:I"/
  DATA KW(i+66)%DSCRPT /"*! *"/

  DATA KW(i+67)%LABEL /"DMFT_IMPOSE_SAME_COEFFS"/
  DATA KW(i+67)%TYP /"L:I"/
  DATA KW(i+67)%DSCRPT /"*! *"/

  DATA KW(i+68)%LABEL /"DMFT_IMPOSE_CHEM_SPIN"/
  DATA KW(i+68)%TYP /"L:I"/
  DATA KW(i+68)%DSCRPT /"*! *"/

  DATA KW(i+69)%LABEL /"DMFT_NBO"/
  DATA KW(i+69)%TYP /"L:I"/
  DATA KW(i+69)%DSCRPT /"*! *"/

  DATA KW(i+70)%LABEL /"DMFT_DOPING"/
  DATA KW(i+70)%TYP /"D:I"/
  DATA KW(i+70)%DSCRPT /"*! additional charge or doping when imposing the chemical potential*"/

  DATA KW(i+71)%LABEL /"DMFT_IGNORE_TYPE"/
  DATA KW(i+71)%TYP /"L:I"/
  DATA KW(i+71)%DSCRPT /"*!  *"/

!END CW

end module esdf_key
