program ttest
implicit none

 real(8),parameter      ::  me_in_ev = 510998.91d0
 real(8),parameter      ::  small_scale_derivative=0.005d0, hartree_in_eV=27.211396132d0, eV_in_hartree=0.036749308824d0
 real(8),parameter      ::  pi       =  3.14159265358979323846264338327950288419716939937510d0
 real(8),parameter      ::  one_angstrom_in_bohr_radius=1.889725989d0
 real(8),parameter      ::  hbar_divided_electron_mass_in_meter2_per_second = 0.000115767622d0
 real(8),parameter      ::  hbar_divided_by_second_in_ev=6.58211899131*1.d-16
 real(8),parameter      ::  hbar_divided_by_me_in_Angstrom_square_divided_by_second=1.15767635048*1.d16
 real(8),parameter      ::  one_over_Angstrom_in_one_over_cm=100000000.d0
 complex(8),parameter   ::  imi      =  CMPLX( 0.d0 , 1.0d0, 8 )
 real(8),parameter      ::  MAX_EXP  =  700.d0
 real(8),parameter      ::  MIN_EXP  = -700.d0
 integer                ::  rank,size2,ierr
 real(8),parameter      ::  Angst_in_ev             = 1968.d0
 real(8),parameter      ::  Angst_in_cm             = 1.d-8
 real(8),parameter      ::  hbar_divided_e2_in_ohm  = 25812.d0/2.d0/pi
 real(8),parameter      ::  cm_moins_1_in_ev        = 0.000123980262342235
 real(8),parameter      ::  Kelvin_in_ev            = 8.617343d-5
 real(8),parameter      ::  ev_in_Kelvin            = 11604.505008098
 real(8)                ::  units_for_optics,volume

    volume=    3171.23691267275     

    units_for_optics =                pi /  Volume
    units_for_optics = units_for_optics  / (hartree_in_eV**2)
    units_for_optics = units_for_optics  / hbar_divided_e2_in_ohm
    units_for_optics = units_for_optics  * (hbar_divided_by_second_in_ev**2)
    units_for_optics = units_for_optics  * (hbar_divided_by_me_in_Angstrom_square_divided_by_second**2)
    units_for_optics = units_for_optics  * (one_over_Angstrom_in_one_over_cm)
    units_for_optics = units_for_optics  * (one_angstrom_in_bohr_radius**2)

    write(*,*) units_for_optics
    write(*,*) 'old scaling : ', 1.05904848572949*one_angstrom_in_bohr_radius**2

end program
