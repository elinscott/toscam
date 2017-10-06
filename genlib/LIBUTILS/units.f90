module units

 use linalg

 !------------------------------------------------------!
 ! relations:                                           ! 
 ! phi_0=h/2/e                                          !
 ! e^2/(4*pi*Epsilon0*c*hbar)                           !
 ! Neff = int_cutoff(opt_cond) *  2 me V / e^2 hbar     !
 ! hbar/e^2 = 4000 Ohm                                  !
 !------------------------------------------------------!

  real(8),parameter,private :: Angst_in_ev             = 1968.d0
  real(8),parameter,private :: Angst_in_cm             = 1.d-8
  real(8),parameter,private :: me_in_ev                = 511000.d0
  real(8),parameter,private :: hbar_divided_e2_in_ohm  = 25812.d0/2.d0/pi
  real(8),parameter,private :: mu_b_in_ev_by_tesla     = 5.7*1.d-5
  real(8),parameter,private :: cm_moins_1_in_ev        = 0.000123980262342235 
  real(8),parameter,private :: Kelvin_in_ev            = 8.617343d-5  
  real(8),parameter,private :: ev_in_Kelvin            = 11604.505008098

contains

 !-------------------------------------------------------------------------------------------------------------------------!
 ! sigma  =  e^2/hbar * 1/c /(Nk*b*a) * pi * Sum_k Int_E (fermi(E)-fermi(E+omega))/omega * Trace(vk rokE vk rok(E+omega))  !
 ! sigma  =  spin *sigma                                                                                                   !  
 ! j_ki   =  derivate H_k                                                                                                  !
 ! N/V * densofstates                                                                                                      !
 ! densofstate = 1/N * Sum_k , normalized                                                                                  !
 ! (1,0) dans la base : (1,+-1)                                                                                            !
 ! (1,0) = (v1+v2)/2                                                                                                       !
 ! (0,1) = (v1-v2)/2                                                                                                       !
 ! BZ : k_x,k_y : norme 2pi/sqrt(2)                                                                                        !
 !-------------------------------------------------------------------------------------------------------------------------!

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 real(8) function optical_conductivity_conv_in_ohm_by_cm_(with_spin,dist_plane_angstrom,VBZ,aa,bb,as,bs)
 implicit none
 logical ::  with_spin
 real(8) ::  dist_plane_angstrom,VBZ,aa(3),bb(3),as(3),bs(3)
  optical_conductivity_conv_in_ohm_by_cm_ = 1.d0/(dist_plane_angstrom*Angst_in_cm)/hbar_divided_e2_in_ohm * &
                                           & pi/norme(vecprod(aa,bb)) * VBZ/(4.d0*pi**2) 
  if(with_spin) optical_conductivity_conv_in_ohm_by_cm_=optical_conductivity_conv_in_ohm_by_cm_*2.d0
 end function

                !--------------------------------------------------!

 real(8) function optical_conductivity_conv_in_ohm_by_cm(with_spin,dist_plane_angstrom,VBZ,aa,bb,as,bs)
 implicit none
 logical ::  with_spin
 real(8) ::  dist_plane_angstrom,VBZ,aa(3),bb(3),as(3),bs(3)
  optical_conductivity_conv_in_ohm_by_cm = 1.d0/(dist_plane_angstrom*Angst_in_cm)/hbar_divided_e2_in_ohm * pi/norme(vecprod(aa,bb)) 
  if(with_spin) optical_conductivity_conv_in_ohm_by_cm=optical_conductivity_conv_in_ohm_by_cm*2.d0
 end function


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 real(8) function optical_conductivity_Neff_(with_spin,a,b,VBZ,aa,bb,as,bs)
 implicit none
 logical ::  with_spin
 real(8) ::  a,b,VBZ,aa(3),bb(3),as(3),bs(3)
   optical_conductivity_Neff_ = a*b/Angst_in_ev/Angst_in_ev / pi * 2.d0 * me_in_ev * norme(vecprod(as,bs))
   optical_conductivity_Neff_ = optical_conductivity_Neff_     * pi/norme(vecprod(aa,bb)) * VBZ/(4.d0*pi**2)
   if(with_spin) optical_conductivity_Neff_=optical_conductivity_Neff_*2.d0
 end function

                !--------------------------------------------------!

 real(8) function optical_conductivity_Neff(with_spin,a,b,VBZ,aa,bb,as,bs)
 implicit none
 logical ::  with_spin
 real(8) ::  a,b,VBZ,aa(3),bb(3),as(3),bs(3)
   optical_conductivity_Neff = a*b/Angst_in_ev/Angst_in_ev / pi * 2.d0 * me_in_ev * norme(vecprod(as,bs))
   optical_conductivity_Neff = optical_conductivity_Neff * pi/norme(vecprod(aa,bb)) 
   if(with_spin) optical_conductivity_Neff=optical_conductivity_Neff*2.d0
 end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

end module
