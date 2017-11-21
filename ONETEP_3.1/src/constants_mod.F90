! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by Chris-Kriton Skylaris
!
!   February 2000
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module constants

  integer, parameter :: DP = kind(1.0d0)               ! Double precision real type
  integer, parameter :: LONG_R = 18                    ! jd: Needed elsewhere
  integer, parameter :: LONG=selected_int_kind(LONG_R) ! Long integer type

  real(kind=DP), parameter :: PI=3.141592653589793238462643383279502884197_DP
  real(kind=DP), parameter :: SQRT_PI=1.772453850905516027298167483341145182798_DP
  real(kind=dp), parameter :: TWO_PI=2.0_DP * pi

#ifndef OLD_CONSTANTS

  ! ndmh 11/02/2011:
  ! For consistency with CASTEP pseudopotentials, the following conversion
  ! factors have been obtained by calling CASTEP's internal conversion routines
  ! as follows:
  !   ANGSTROM = io_atomic_to_unit(1.0_dp,'ang')
  !   HARTREE_IN_EVS = io_atomic_to_unit(1.0_dp,'eV')
  ! based on the CODATA 2002 values (CASTEP's default choice). This produced
  !   electron_mass_si*speed_light_si*fine_structure_si/hbar_si*1e-10_dp
  !   elementary_charge_si/(fine_structure_si**2*electron_mass_si*speed_light_si**2)
  ! where these were derived from the following definitions (CODATA2002):
  !real(kind=dp), parameter, public :: speed_light_si = 299792458.0_dp
  !real(kind=dp), parameter, public :: planck_si = 6.6260693e-34_dp
  !real(kind=dp), parameter, public :: elementary_charge_si = 1.60217653e-19_dp
  !real(kind=dp), parameter, public :: electron_mass_si = 9.1093826e-31_dp
  !real(kind=dp), parameter, public :: avogadro_si = 6.0221415e23_dp
  !real(kind=dp), parameter, public :: pi=3.141592653589793238462643383279502884197_dp
  !real(kind=dp), parameter, public :: two_pi=2.0_dp*pi
  !real(kind=dp), parameter, public :: hbar_si = planck_si/two_pi
  !real(kind=dp), parameter, public :: mu_0_si = 4.0_dp*pi*1e-7_dp
  !real(kind=dp), parameter, public :: epsilon_0_si = 1.0_dp/(mu_0_si*speed_light_si**2)
  !real(kind=dp), parameter, public :: fine_structure_si = elementary_charge_si**2
  !                                 / (4.0_dp*pi*epsilon_0_si*hbar_si*speed_light_si)

  ! The value of one Angstrom in terms of Bohr radii
  real(kind=DP), parameter :: ANGSTROM=1.889726134583548707935_DP

  ! cks: The value of one Hartree in terms of electron-volts
  real(kind=DP), parameter :: HARTREE_IN_EVS=27.2113846081672_DP

#else

  ! ndmh: equivalent values in terms of the old versions of the CASTEP constants
  real(kind=DP), parameter :: ANGSTROM=1.889726313_DP
  real(kind=DP), parameter :: HARTREE_IN_EVS=27.2116529_DP

#endif

  ! pdh: square root of minus one
  complex(kind=DP), parameter :: cmplx_i = (0.0_DP,1.0_DP)

  ! ndmh: complex 0
  complex(kind=DP), parameter :: cmplx_0 = (0.0_DP,0.0_DP)
  complex(kind=DP), parameter :: cmplx_1 = (1.0_DP,0.0_DP)

  ! pdh: units for standard output and standard error
  integer, parameter :: stdout = 6
  integer, parameter :: stderr = stdout

  ! ndmh: sizes of datatypes, for size estimation
  integer, parameter :: logical_size = 1
  integer, parameter :: char_size = 1
  integer, parameter :: int_size = 4
  integer, parameter :: real_size = 8
  integer, parameter :: cmplx_size = 16

  ! cks: output level integer constants
  integer, parameter :: BRIEF   =0
  integer, parameter :: NORMAL  =1
  integer, parameter :: VERBOSE =2

  ! pdh: spin polarisation
  integer, parameter :: max_spins = 2
  integer, parameter :: UP = 1
  integer, parameter :: DN = 2

  ! ndmh: PAW energy component labels
  integer, parameter, public :: paw_en_size     = 9
  integer, parameter, public :: paw_en_dijhat   = 1
  integer, parameter, public :: paw_en_dij0     = 2
  integer, parameter, public :: paw_en_ehart    = 3
  integer, parameter, public :: paw_en_exc      = 4
  integer, parameter, public :: paw_en_exc_dc   = 5
  integer, parameter, public :: paw_en_etxc     = 6
  integer, parameter, public :: paw_en_etxc_dc  = 7
  integer, parameter, public :: paw_en_dijxc    = 8
  integer, parameter, public :: paw_en_exc_core = 9

  ! aam: The symbols of the elements in the periodic table
  character(len=2), parameter, dimension(109) :: periodic_table_name= (/ &
       & 'H ',                                                                                'He', &
       & 'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &
       & 'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &
       & 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
       & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
       & 'Cs','Ba', &
       & 'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu', &
       & 'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
       & 'Fr','Ra', &
       & 'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr', &
       & 'Rf','Db','Sg','Bh','Hs','Mt' /)

  ! aam: The atomic masses of the elements in the periodic table
  real(kind=dp), parameter, dimension(109) :: periodic_table_mass = (/ &
     & 1.00794_dp,                                                                                                    4.00260_dp, &
     & 6.941_dp, 9.012187_dp,                            10.811_dp, 12.0107_dp, 14.00674_dp, 15.9994_dp, 18.99840_dp, 20.1797_dp, &
     & 22.98977_dp, 24.3050_dp,                           26.98154_dp, 28.0855_dp, 30.97376_dp, 32.066_dp, 35.4527_dp, 39.948_dp, &
     & 39.0983_dp, 40.078_dp, &
     &   44.95591_dp, 47.867_dp, 50.9415_dp, 51.9961_dp, 54.93805_dp, 55.845_dp, 58.93320_dp, 58.6934_dp, 63.546_dp, 65.39_dp, &
     &                                                           69.723_dp, 72.61_dp, 74.92160_dp, 78.96_dp, 79.904_dp, 83.80_dp, &
     & 85.4678_dp, 87.62_dp, &
     &   88.90585_dp, 91.224_dp, 92.90638_dp, 95.94_dp, 98.0_dp, 101.07_dp, 102.90550_dp, 106.42_dp, 107.8682_dp, 112.411_dp, &
     &                                                    114.818_dp, 118.710_dp, 121.760_dp, 127.60_dp, 126.90447_dp, 131.29_dp, &
     & 132.90545_dp, 137.327_dp, &
     &   138.9055_dp, 140.116_dp, 140.90765_dp, 144.24_dp, 145.0_dp, 150.36_dp, 151.964_dp, 157.25_dp, 158.92534_dp, &
     &                                    162.50_dp, 164.93032_dp, 167.26_dp, 168.93421_dp, 173.04_dp, 174.967_dp, &
     &   178.49_dp, 180.9479_dp, 183.84_dp, 186.207_dp, 190.23_dp, 192.217_dp, 195.078_dp, 196.96655_dp, 200.59_dp, &
     &                                                      204.3833_dp, 207.2_dp, 208.98038_dp, 209.0_dp, 210.0_dp, 222.0_dp, &
     & 223.0_dp, 226.0_dp, &
     &   227.0_dp, 232.0381_dp, 231.03588_dp, 238.0289_dp, 237.0_dp, 244.0_dp, 243.0_dp, 247.0_dp, 247.0_dp, 251.0_dp, 252.0_dp, &
     &                                                                           257.0_dp, 258.0_dp, 259.0_dp, 262.0_dp, &
     &   261.0_dp, 262.0_dp, 263.0_dp, 264.0_dp, 265.0_dp, 268.0_dp /)

end module constants


