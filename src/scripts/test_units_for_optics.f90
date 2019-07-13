program ttest
   implicit none

   real(kind=DP), parameter      ::  me_in_ev = 510998.91d0
   real(kind=DP), parameter      ::  small_scale_derivative = 0.005d0, hartree_in_eV = 27.211396132d0, eV_in_hartree = 0.036749308824d0
   real(kind=DP), parameter      ::  pi = 3.14159265358979323846264338327950288419716939937510d0
   real(kind=DP), parameter      ::  one_angstrom_in_bohr_radius = 1.889725989d0
   real(kind=DP), parameter      ::  hbar_divided_electron_mass_in_meter2_per_second = 0.000115767622d0
   real(kind=DP), parameter      ::  hbar_divided_by_second_in_ev = 6.58211899131*1.d-16
   real(kind=DP), parameter      ::  hbar_divided_by_me_in_Angstrom_square_divided_by_second = 1.15767635048*1.d16
   real(kind=DP), parameter      ::  one_over_Angstrom_in_one_over_cm = 100000000.d0
   complex(kind=DP), parameter   ::  imi = CMPLX(0.d0, 1.0d0, 8)
   real(kind=DP), parameter      ::  MAX_EXP = 700.d0
   real(kind=DP), parameter      ::  MIN_EXP = -700.d0
   integer                ::  rank, size2, ierr
   real(kind=DP), parameter      ::  Angst_in_ev = 1968.d0
   real(kind=DP), parameter      ::  Angst_in_cm = 1.d-8
   real(kind=DP), parameter      ::  hbar_divided_e2_in_ohm = 25812.d0/2.d0/pi
   real(kind=DP), parameter      ::  cm_moins_1_in_ev = 0.000123980262342235
   real(kind=DP), parameter      ::  Kelvin_in_ev = 8.617343d-5
   real(kind=DP), parameter      ::  ev_in_Kelvin = 11604.505008098
   real(kind=DP)                ::  units_for_optics, volume

   volume = 3171.23691267275

   units_for_optics = pi/Volume
   units_for_optics = units_for_optics/(hartree_in_eV**2)
   units_for_optics = units_for_optics/hbar_divided_e2_in_ohm
   units_for_optics = units_for_optics*(hbar_divided_by_second_in_ev**2)
   units_for_optics = units_for_optics*(hbar_divided_by_me_in_Angstrom_square_divided_by_second**2)
   units_for_optics = units_for_optics*(one_over_Angstrom_in_one_over_cm)
   units_for_optics = units_for_optics*(one_angstrom_in_bohr_radius**2)

   write (*, *) units_for_optics
   write (*, *) 'old scaling : ', 1.05904848572949*one_angstrom_in_bohr_radius**2

end program
