! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!==============================================================================!
!                                                                              !
!  ONETEP exchange-correlation functionals module: xc_funcs_mod.F90            !
!                                                                              !
!  The subroutines in this file were written by Peter Haynes, Quintin Hill,    !
!  Nicholas Hine, and Jacek Dziedzic, between February 2004 and October 2011.  !
!                                                                              !
!  Addition of hybrid functionals, BLYP and XLYP, by Quintin Hill in           !
!  February/March 2009 under the supervision of Chris-Kriton Skylaris.         !
!                                                                              !
!==============================================================================!

module xc_funcs

  use constants, only: DP

  implicit none

  real(kind=DP), parameter :: dentol = 1.0e-15_DP
  real(kind=DP), parameter :: THIRD = 1.0_DP / 3.0_DP
  real(kind=DP), parameter :: FTHRD = 4.0_DP / 3.0_DP

  private

  ! Public subroutines
  public :: xc_lda_x_point
  public :: xc_lsda_x_point
  public :: xc_capz_c_point
  public :: xc_capz_c_point_sp
  public :: xc_vwn_c_point
  public :: xc_vwn_c_point_sp
  public :: xc_pw92_c_point
  public :: xc_pw92_c_point_sp
  public :: xc_pbe_x_point
  public :: xc_rpbe_x_point
  public :: xc_revpbe_x_point
  public :: xc_wc_x_point
  public :: xc_pbe0_x_point
  public :: xc_pbe_c_point
  public :: xc_pbe_x_point_sp
  public :: xc_rpbe_x_point_sp
  public :: xc_revpbe_x_point_sp
  public :: xc_wc_x_point_sp
  public :: xc_pbe0_x_point_sp
  public :: xc_pbe_c_point_sp
  public :: xc_pw91_x_point
  public :: xc_pw91_c_point
  public :: xc_pw91_x_point_sp
  public :: xc_pw91_c_point_sp
  public :: xc_b88_x_point
  public :: xc_b88_x_point_sp
  public :: xc_b1_x_point
  public :: xc_b3_x_point
  public :: xc_x_x_point
  public :: xc_x3_x_point
  public :: xc_b1_x_point_sp
  public :: xc_b3_x_point_sp
  public :: xc_x_x_point_sp
  public :: xc_x3_x_point_sp
  public :: xc_lyp_c_point
  public :: xc_lyp_c_point_sp
  public :: xc_b3lyp_c_point
  public :: xc_b3lyp_c_point_sp
  public :: xc_b3pw91_c_point
  public :: xc_b3pw91_c_point_sp
  public :: xc_x3lyp_c_point
  public :: xc_x3lyp_c_point_sp
  public :: xc_none_c_point
  public :: xc_none_c_point_sp

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines for evaluating the energy and potential contibution of each
! point.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lda_x_point(den,x_energy,x_pot)

    !==========================================================================!
    ! This subroutine calculates the local density appoximation (LDA) exchange !
    ! energy and potential at a point given the density at that point.         !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in March 2009 based on existing xc_capz          !
    ! subroutine written by Peter Haynes in February 2004.                     !
    !==========================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den      ! Density
    real(kind=DP), intent(out) :: x_energy ! Exchange energy
    real(kind=DP), intent(out) :: x_pot    ! Exchange potential

    real(kind=DP) :: rs
    real(kind=DP) :: epsx
    real(kind=DP), parameter :: xf = -0.458165293283142893475554_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    if (den > 0.0_DP) then   ! Positive charge density

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       epsx = xf / rs
       x_energy = den * epsx
       x_pot = FTHRD * epsx
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
    end if

  end subroutine xc_lda_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lsda_x_point(den1,den2,x_energy,x_pot1,x_pot2)

    !==========================================================================!
    ! This subroutine calculates the local spin density appoximation (LSDA)    !
    ! exchange energy and potential at a point given the density at that point !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in March 2009 based on existing xc_capz          !
    ! subroutine written by Peter Haynes in February 2004.                     !
    !==========================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den1
    real(kind=DP), intent(in)  :: den2
    real(kind=DP), intent(out) :: x_energy
    real(kind=DP), intent(out) :: x_pot1
    real(kind=DP), intent(out) :: x_pot2

    real(kind=DP) :: epsx1
    real(kind=DP) :: epsx2
    real(kind=DP), parameter :: xfp = -0.930525736349100025002010_DP

    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    ! Check for Positive charge densities
    if (den1 >= 0.0_DP) then
       epsx1 = xfp * den1**THIRD
       x_energy = den1 * epsx1
       x_pot1 = FTHRD * epsx1
    end if
    if (den2 >= 0.0_DP) then
       epsx2 = xfp * den2**THIRD
       x_energy = x_energy + den2 * epsx2
       x_pot2 = FTHRD * epsx2
    end if

  end subroutine xc_lsda_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_capz_c_point(den,c_energy,c_pot)

    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the exchange-correlation energy and potential, according to the local     !
    ! density approximation (LDA), from the Monte Carlo data by Ceperley and    !
    ! Alder, the Gell-Mann & Brueckner expansion, and the parametrisation by    !
    ! Perdew & Zunger.                                                          !
    !                                                                           !
    ! References:                                                               !
    !  D M Ceperley & B J Alder, Phys. Rev. Lett. 45, 566 (1980)                !
    !  M Gell-Mann & K Brueckner, Phys. Rev. 106, 364 (1957)                    !
    !     see also N W Ashcroft & N D Mermin 'Solid State Physics',             !
    !              International Edition, p. 336                                !
    !  J P Perdew & A Zunger, Phys. Rev. B 23, 5048 (1981)                      !
    !---------------------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.                        !
    ! Split into point contribution by Quintin Hill on 09/03/2009.              !
    !===========================================================================!


    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den
    real(kind=DP), intent(out) :: c_energy
    real(kind=DP), intent(out) :: c_pot

    real(kind=DP) :: rs           ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs       ! Square root of rs
    real(kind=DP) :: logrs        ! Natural logarithm of rs
    real(kind=DP) :: epsc         ! Correlation energy per particle

    ! Constants
    real(kind=DP), parameter :: SSXTH = 7.0_DP / 6.0_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    ! Constants for spin-unpolarised case:
    real(kind=DP), parameter :: g = -0.1423_DP
    real(kind=DP), parameter :: b1 = 1.0529_DP
    real(kind=DP), parameter :: b2 = 0.3334_DP
    real(kind=DP), parameter :: aa = 0.0311_DP
    real(kind=DP), parameter :: bb = -0.048_DP
    real(kind=DP), parameter :: cc = 0.0020_DP
    real(kind=DP), parameter :: dd = -0.0116_DP

    if (den > 0.0_DP) then   ! Positive charge density

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius

       if (rs >= 1.0_DP) then  ! Low charge density

          sqrtrs = sqrt(rs)
          epsc = g / (1.0_DP + b1 * sqrtrs + b2 * rs)
          c_energy = den * epsc
          c_pot = epsc * epsc * (1.0_DP + SSXTH * b1 * sqrtrs + &
               FTHRD * b2 * rs) / g

       else  ! High charge density

          logrs = log(rs)
          epsc = aa * logrs + bb + &
               (cc * logrs + dd) * rs
          c_energy = den * epsc
          c_pot = epsc - THIRD * (aa + &
               (cc * logrs + cc + dd) * rs)

       end if

    else   ! Negative charge density
       c_energy = 0.0_DP
       c_pot = 0.0_DP

    end if

  end subroutine xc_capz_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_capz_c_point_sp(den,den1,den2,c_energy,c_pot1,c_pot2)

    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the exchange-correlation energy and potential, according to the local     !
    ! density approximation (LDA), from the Monte Carlo data by Ceperley and    !
    ! Alder, the Gell-Mann & Brueckner expansion, and the parametrisation by    !
    ! Perdew & Zunger.                                                          !
    !                                                                           !
    ! References:                                                               !
    !  D M Ceperley & B J Alder, Phys. Rev. Lett. 45, 566 (1980)                !
    !  M Gell-Mann & K Brueckner, Phys. Rev. 106, 364 (1957)                    !
    !     see also N W Ashcroft & N D Mermin 'Solid State Physics',             !
    !              International Edition, p. 336                                !
    !  J P Perdew & A Zunger, Phys. Rev. B 23, 5048 (1981)                      !
    !---------------------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.                        !
    ! Split into point contribution by Quintin Hill on 09/03/2009.              !
    !===========================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den
    real(kind=DP), intent(in)  :: den1
    real(kind=DP), intent(in)  :: den2
    real(kind=DP), intent(out) :: c_energy
    real(kind=DP), intent(out) :: c_pot1
    real(kind=DP), intent(out) :: c_pot2

    real(kind=DP) :: zeta         ! Spin polarisation
    real(kind=DP) :: fzeta, dfdz  ! Function of zeta and its derivative
    real(kind=DP) :: rs           ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs       ! Square root of rs
    real(kind=DP) :: logrs        ! Natural logarithm of rs
    real(kind=DP) :: epsc         ! Correlation energy per particle
    real(kind=DP) :: epscp        ! Polarised correlation energy per particle
    real(kind=DP) :: vcu,vcp   ! Unpolarised and polarised correlation potentials

    ! Constants
    real(kind=DP), parameter :: SSXTH = 7.0_DP / 6.0_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    ! Constants for spin-unpolarised case:
    real(kind=DP), parameter :: g = -0.1423_DP
    real(kind=DP), parameter :: b1 = 1.0529_DP
    real(kind=DP), parameter :: b2 = 0.3334_DP
    real(kind=DP), parameter :: aca = 0.0311_DP
    real(kind=DP), parameter :: bca = -0.048_DP
    real(kind=DP), parameter :: cca = 0.0020_DP
    real(kind=DP), parameter :: dca = -0.0116_DP

    ! Constants for spin-polarised case:
    real(kind=DP), parameter :: gp = -0.0843_DP
    real(kind=DP), parameter :: b1p = 1.3981_DP
    real(kind=DP), parameter :: b2p = 0.2611_DP
    real(kind=DP), parameter :: acap = 0.01555_DP
    real(kind=DP), parameter :: bcap = -0.0269_DP
    real(kind=DP), parameter :: ccap = 0.0007_DP
    real(kind=DP), parameter :: dcap = -0.0048_DP

    ! Positive charge densities
    if (den > 0.0_DP) then

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius

       if (rs >= 1.0_DP) then    ! Low charge density

          sqrtrs = sqrt(rs)
          epsc = g / (1.0_DP + b1 * sqrtrs + b2 * rs)
          vcu = epsc * epsc * (1.0_DP + SSXTH * b1 * sqrtrs + &
               FTHRD * b2 * rs) / g
          epscp = gp / (1.0_DP + b1p * sqrtrs + b2p * rs)
          vcp = epscp * epscp * (1.0_DP + SSXTH * b1p * sqrtrs + &
               FTHRD * b2p * rs) / gp

       else   ! High charge density

          logrs = log(rs)
          epsc = aca * logrs + bca + (cca * logrs + dca) * rs
          vcu = epsc - &
               (aca + (cca * logrs + cca + dca) * rs) * THIRD
          epscp = acap * logrs + bcap + (ccap * logrs + dcap) * rs
          vcp = epscp - &
               (acap + (ccap * logrs + ccap + dcap) * rs) * THIRD

       end if

       zeta = (den1 - den2) / den   ! Polarisation
       !pdh: check that -1 <= zeta <= 1
       zeta = min(zeta,1.0_DP)
       zeta = max(zeta,-1.0_DP)

       fzeta = ((1.0_DP + zeta)**FTHRD + &
            (1.0_DP - zeta)**FTHRD - 2.0_DP) / &
            (2.0_DP**FTHRD - 2.0_DP)
       dfdz = FTHRD * ((1.0_DP + zeta)**THIRD - &
            (1.0_DP - zeta)**THIRD) / &
            (2.0_DP**FTHRD - 2.0_DP)

       c_pot1 = vcu + fzeta * (vcp - vcu) + &
            (epscp - epsc) * (1.0_DP - zeta) * dfdz
       c_pot2 = vcu + fzeta * (vcp - vcu) - &
            (epscp - epsc) * (1.0_DP + zeta) * dfdz

       c_energy = den * (epsc + fzeta * (epscp - epsc))

    else  ! Negative charge density
       c_energy = 0.0_DP
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP
    end if

  end subroutine xc_capz_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vwn_c_point(den,c_energy,c_pot)
    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the correlation energy and potential, according to the local density      !
    ! approximation (LDA), from the parameterisation by Vosko, Wilk and Nusair. !
    !                                                                           !
    ! Spin-unpolarized version.                                                 !
    !                                                                           !
    ! Reference:                                                                !
    !  [1] S H Vosko, L Wilk & M Nusair, Can. J. Phys. 58, 1200 (1980)          !
    !  [2] S H Vosko and L Wilk, Phys. Rev. B 22, 3812 (1980)                   !
    !                                                                           !
    !---------------------------------------------------------------------------!
    ! Written in 2011/05/31 by Jacek Dziedzic, basing on a template by Peter    !
    ! Haynes, modified by Quintin Hill.                                         !
    !===========================================================================!
    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den
    real(kind=DP), intent(out) :: c_energy
    real(kind=DP), intent(out) :: c_pot

    real(kind=DP) :: rs ! Wigner-Seitz radius corresponding to den

    ! Constants
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: minus_fourpi_over_nine = -4.0_DP/9.0_DP * PI
    real(kind=DP), parameter :: MIN_DENS = 5.0D-13

    ! Local variables
    real(kind=DP) :: c_eps     ! Correlation energy density
    real(kind=DP) :: drs_drho  ! Derivative of rs wrt electronic density

    ! -------------------------------------------------------------------------

    ! We're only concerned about positive densities.
    ! Furthermore, to get perfect agreement with LIBXC, we ignore densities
    ! below MIN_DENS just like they do
    if (den > MIN_DENS) then

       rs = (TONFPI / den)**THIRD ! the Wigner-Seitz radius

       ! Calculate energy and potential (eq. (4.4) in [1] and its derivative)
       call xc_vwn_eps_c_helper(c_eps,c_pot,rs,1)

       ! Energy density -> energy
       c_energy = c_eps * den

       ! drs/drho
       drs_drho = minus_fourpi_over_nine * rs**4.0_DP

       ! Resulting potential (functional derivative of E_c (not eps_c) wrt rho)
       ! (v_c = eps_c + rho * de_c/drho), where de_c/drho = de_c/drs * drs/drho
       c_pot = c_eps + den * c_pot * drs_drho

    else
       c_energy = 0.0_DP
       c_pot = 0.0_DP
    end if

  end subroutine xc_vwn_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vwn_c_point_sp(den,den1,den2,c_energy,c_pot1,c_pot2)
    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the correlation energy and potential, according to the local density      !
    ! approximation (LDA), from the parameterisation by Vosko, Wilk and Nusair. !
    !                                                                           !
    ! Spin-polarized version.                                                   !
    !                                                                           !
    ! Reference:                                                                !
    !  [1] S H Vosko, L Wilk & M Nusair, Can. J. Phys. 58, 1200 (1980)          !
    !  [2] S H Vosko and L Wilk, Phys. Rev. B 22, 3812 (1980)                   !
    !                                                                           !
    !---------------------------------------------------------------------------!
    ! Written in 2011/05/31 by Jacek Dziedzic, basing on a template by Peter    !
    ! Haynes, modified by Quintin Hill.                                         !
    !===========================================================================!
    use constants, only: DP, PI
    use utils, only: utils_abort
    implicit none

    ! Arguments
    real(kind=DP), intent(in)  :: den, den1, den2
    real(kind=DP), intent(out) :: c_energy
    real(kind=DP), intent(out) :: c_pot1, c_pot2

    ! Constants
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: d2f_at_zero = &
         1.7099209341613656176_DP ! f''(0) = 4/(9 * (-1 + 2**(1/3))), cf. [1]
    real(kind=DP), parameter :: minus_fourpi_over_nine = -4.0_DP/9.0_DP * PI
    real(kind=DP), parameter :: MIN_DENS = 5.0D-13

    ! Local variables
    real(kind=DP) :: c_eps     ! Correlation energy density
    real(kind=DP) :: rs        ! Wigner-Seitz radius
    real(kind=DP) :: alpha_c   ! See discussion below (9) in [2]
    real(kind=DP) :: eps_f     ! Ferromagnetic contribution to eps
    real(kind=DP) :: eps_p     ! Paramagnetic contribution to eps
    real(kind=DP) :: pot_c, pot_f, pot_p ! Corresponding 1st derivatives
    real(kind=DP) :: eps_delta ! (7) in [2]
    real(kind=DP) :: pot_delta ! (7) in [2] with epsilons replaced by derivs
    real(kind=DP) :: beta      ! (8) in [2]
    real(kind=DP) :: beta_pot  ! (8) in [2] with epsilons replaced by derivs
    real(kind=DP) :: f_of_zeta ! (4) in [2]
    real(kind=DP) :: fprime_of_zeta ! d/dzeta of (4) in [2]
    real(kind=DP) :: zeta      ! Polarization
    real(kind=DP) :: zeta3     ! zeta**3
    real(kind=DP) :: zeta4     ! zeta**4
    real(kind=DP) :: de_drs    ! deps/drs
    real(kind=DP) :: de_dzeta, de_dzeta1, de_dzeta2 ! deps/dzeta and its terms
    real(kind=DP) :: c_pot     ! Resulting potential, before zeta correction
    real(kind=DP) :: drs_drho  ! Derivative of rs wrt electronic density

    ! -------------------------------------------------------------------------

    ! We're only concerned about positive densities.
    ! Furthermore, to get perfect agreement with LIBXC, we ignore densities
    ! below MIN_DENS just like they do
    if (den > MIN_DENS) then

       rs = (TONFPI / den)**THIRD ! the Wigner-Seitz radius

       ! Calculate spin-stiffness contribution, alpha_c(r_s)
       ! see discussion below (9) in [2].
       call xc_vwn_eps_c_helper(alpha_c,pot_c,rs,3)

       ! Calculate paramagnetic contribution eps_p
       ! see discussion below (9) in [2].
       call xc_vwn_eps_c_helper(eps_p,pot_p,rs,1)

       ! Calculate ferromagnetic contribution eps_f
       ! see discussion below (9) in [2].
       call xc_vwn_eps_c_helper(eps_f,pot_f,rs,2)

       ! This should never happen. Since this is executed many times, an
       ! explicit check is better than an assert.
       if(alpha_c == 0.0_DP) then
          call utils_abort("Division by zero in xc_vwn_c_point_sp")
       end if

       ! (8) in [2]
       beta = d2f_at_zero * (eps_f - eps_p) / alpha_c - 1.0_DP
       beta_pot = d2f_at_zero * (pot_f - pot_p) / pot_c - 1.0_DP

       ! Polarization
       zeta = (den1 - den2) / den
       zeta = min(zeta,1.0_DP)
       zeta = max(zeta,-1.0_DP)

       ! zeta**3 and zeta**4
       zeta3 = zeta*zeta*zeta
       zeta4 = zeta*zeta3

       ! (4) in [2]
       f_of_zeta = 0.5_DP * ((1+zeta)**(4.0_DP/3.0_DP) + &
            (1-zeta)**(4.0_DP/3.0_DP) - 2.0_DP) / (2.0_DP**THIRD - 1.0_DP)

       ! 1st derivative of the above
       fprime_of_zeta = ( -4.0_DP/3.0_DP*(1.0_DP-zeta)**THIRD + &
            4.0_DP/3.0_DP*(1.0_DP+zeta)**THIRD ) / &
            (2.0_DP * (-1.0_DP + 2.0_DP**THIRD))

       ! Corrections to energy density and potential due to polarization:
       ! - delta_eps_c
       eps_delta = alpha_c * f_of_zeta/d2f_at_zero * (1.0_DP + beta * zeta4)

       ! - delta_v_c
       pot_delta = pot_c * f_of_zeta/d2f_at_zero * (1.0_DP + beta_pot * zeta4)

       ! Corrected energy density: eps_paramagnetic + delta_eps
       c_eps = eps_p + eps_delta

       ! Corrected derivative of c_eps wrt rs
       de_drs = pot_p + pot_delta

       ! derivative of c_eps wrt zeta (= derivative of eps_delta wrt zeta)
       ! term1 = d/dzeta [ alpha_c * f(zeta)/f''(0) * (1-zeta**4) ]
       ! term2 = d/dzeta [ f(zeta) * zeta**4 * (eps_f - eps_p) ]
       de_dzeta1 = alpha_c / d2f_at_zero * (fprime_of_zeta * (1.0_DP - zeta4) +&
            f_of_zeta * (-4.0_DP*zeta3))
       de_dzeta2 = (fprime_of_zeta * zeta4 + 4.0_DP * zeta3 * f_of_zeta) * &
            (eps_f - eps_p)

       de_dzeta = de_dzeta1 + de_dzeta2

       ! Energy density -> energy
       c_energy = c_eps * den

       ! drs/drho
       drs_drho = minus_fourpi_over_nine * rs**4.0_DP

       ! Resulting potential (functional derivative of E_c (not eps_c) wrt rho)
       ! (v_c = eps_c + rho * de_c/drho), where de_c/drho = de_c/drs * drs/drho
       c_pot = c_eps + den * de_drs * drs_drho

       ! This is the spin-polarized case, take zeta into account
       c_pot1 = c_pot - (zeta - 1.0_DP) * de_dzeta
       c_pot2 = c_pot - (zeta + 1.0_DP) * de_dzeta

    else
       c_energy = 0.0_DP
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP
    end if

  end subroutine xc_vwn_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vwn_eps_c_helper(c_eps, c_pot, rs, paramset)
    !===========================================================================!
    ! Given the Wigner-Seitz radius, rs, calculates the correlation energy      !
    ! density, according to expression (4.4) in [1]. The first derivative       !
    ! wrt rs is also calculated.                                                !
    !---------------------------------------------------------------------------!
    ! Arguments:                                                                !
    !   c_eps    (out): Returned energy density.                                !
    !   c_pot    (out): Returned potential.                                     !
    !   rs       (in):  Wigner-Seitz radius (a simple function of density)      !
    !   paramset (in):  Selects the VWN parameter set:                          !
    !                   1 - paramagnetic case                                   !
    !                   2 - ferromagnetic case                                  !
    !                   3 - calculation of alpha_c (cf. text below (9) in [2].  !
    !---------------------------------------------------------------------------!
    ! Reference:                                                                !
    !  [1] S H Vosko, L Wilk & M Nusair, Can. J. Phys. 58, 1200 (1980)          !
    !  [2] S H Vosko and L Wilk, Phys. Rev. B 22, 3812 (1980)                   !
    !---------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 01/06/2011.                                    !
    !===========================================================================!
    use constants, only: DP, PI
    use utils, only: utils_abort
    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: c_eps
    real(kind=DP), intent(out) :: c_pot
    real(kind=DP), intent(in)  :: rs
    integer                    :: paramset

    ! jd: Parameters of the method, following [1], [2]
    !     Note A(:) are divided by 2 because of Ry->Ha conversion.
    !     This also applies to A(3), which becomes 1/(6 pi^2) rather than
    !     1/(3 pi^2)
    real(kind=DP), parameter :: As(3) = &
         (/ 0.0310907_DP, 0.01554535_DP, -0.01688686394038963_DP /)
    real(kind=DP), parameter :: bs(3) = (/ 3.72744_DP, 7.06042_DP, 1.13107_DP /)
    real(kind=DP), parameter :: cs(3) = (/ 12.9352_DP, 18.0578_DP, 13.0045_DP /)
    real(kind=DP), parameter :: x0s(3) = &
         (/ -0.10498_DP, -0.325_DP, -0.0047584_DP /)
    real(kind=DP), parameter :: Qs(3) = &
         (/ 6.15199082_DP, 4.73092691_DP, 7.123108918_DP /) ! sqrt(4c-b2)
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    ! jd: Internal variables
    real(kind=DP) :: x                 ! jd: =rs^(1/2), cf. [1], below [4.3]
    real(kind=DP) :: X_of_x, X_of_x0   ! jd:            cf. [1], below [4.4]
    real(kind=DP) :: den               ! jd: the density corresponding to rs
    real(kind=DP) :: ln1, ln2, arctan1 ! jd: Temporaries
    real(kind=DP) :: A,b,c,x0,Q        ! jd: Parameters for current paramset

    ! -------------------------------------------------------------------------

    if(paramset < 1 .or. paramset > 3) then
       call utils_abort("Bad paramset in xc_vwn_eps_c_helper")
    endif

    A = As(paramset)
    b = bs(paramset)
    c = cs(paramset)
    x0 = x0s(paramset)
    Q = Qs(paramset)

    x = sqrt(rs) ! rs>0 is guaranteed by the caller
    den = TONFPI / (rs*rs*rs)
    X_of_x = x*x + b*x + c
    X_of_x0 = x0*x0 + b*x0 + c

    ln1 = log(x*x/X_of_x)
    ln2 = log((x-x0)*(x-x0)/X_of_x)
    arctan1 = atan(Q/(2.0_DP*x+b))

    ! Calculate energy density according to (4.4)
    c_eps = A* &
         (ln1 + 2.0_DP*b/Q * arctan1 - b*x0/X_of_x0 * &
         (ln2 + 2.0_DP*(b+2.0_DP*x0)/Q * arctan1))

    ! Derivative of c_eps wrt rs
    c_pot = (A*c*x-A * (c+b*x)*x0) / (rs*(c+b*x+rs)*(x-x0))

  end subroutine xc_vwn_eps_c_helper


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw92_c_point(den,c_energy,c_pot)

    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the exchange-correlation energy and potential, according to the local     !
    ! density approximation (LDA), from the Monte Carlo data by Ceperley and    !
    ! Alder, the Gell-Mann & Brueckner expansion, and the parametrisation by    !
    ! Perdew & Wang.                                                            !
    !                                                                           !
    ! References:                                                               !
    !  D M Ceperley & B J Alder, Phys. Rev. Lett. 45, 566 (1980)                !
    !  M Gell-Mann & K Brueckner, Phys. Rev. 106, 364 (1957)                    !
    !     see also N W Ashcroft & N D Mermin 'Solid State Physics',             !
    !              International Edition, p. 336                                !
    !  J P Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992)                     !
    !---------------------------------------------------------------------------!
    ! Originally written by Nicholas Hine, March 2011.                          !
    !===========================================================================!


    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den
    real(kind=DP), intent(out) :: c_energy
    real(kind=DP), intent(out) :: c_pot

    real(kind=DP) :: rs           ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs       ! Square root of rs
    real(kind=DP) :: ecunif,eurs

    ! Constants
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    ! Constants for spin-unpolarised case:

    tol: if (den > dentol) then   ! Positive charge density
       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,ecunif,eurs)
       c_energy = den * ecunif
       c_pot = ecunif-THIRD*rs*eurs

    else
       c_energy = 0.0_DP
       c_pot = 0.0_DP
    end if tol

  end subroutine xc_pw92_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw92_c_point_sp(den,den1,den2,c_energy,c_pot1,c_pot2)

    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the exchange-correlation energy and potential, according to the local     !
    ! density approximation (LDA), from the Monte Carlo data by Ceperley and    !
    ! Alder, the Gell-Mann & Brueckner expansion, and the parametrisation by    !
    ! Perdew & Zunger.                                                          !
    !                                                                           !
    ! References:                                                               !
    !  D M Ceperley & B J Alder, Phys. Rev. Lett. 45, 566 (1980)                !
    !  M Gell-Mann & K Brueckner, Phys. Rev. 106, 364 (1957)                    !
    !     see also N W Ashcroft & N D Mermin 'Solid State Physics',             !
    !              International Edition, p. 336                                !
    !  J P Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992)                     !
    !---------------------------------------------------------------------------!
    ! Originally written by Nicholas Hine, March 2011.                          !
    !===========================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den
    real(kind=DP), intent(in)  :: den1
    real(kind=DP), intent(in)  :: den2
    real(kind=DP), intent(out) :: c_energy
    real(kind=DP), intent(out) :: c_pot1
    real(kind=DP), intent(out) :: c_pot2

    real(kind=DP) :: zeta         ! Spin polarisation
    real(kind=DP) :: rs           ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs       ! Square root of rs
    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: ecunif               ! LDA correlation energy per particle
    real(kind=DP) :: ecd,ecz,eu,eurs,ecrs
    real(kind=DP) :: eprs,ep,alfm,alfrsm
    real(kind=DP) :: eczet
    real(kind=DP) :: zeta3,zeta4          ! zeta**3 and zeta**4
    real(kind=DP) :: fzeta,dfdz           ! Function of zeta and derivative

    ! Constants
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: zetamin = -1.0_DP + 1.0e-10_DP
    real(kind=DP), parameter :: zetamax = 1.0_DP - 1.0e-10_DP
    real(kind=DP), parameter :: gaminv = 1.92366105093153631981_DP
    real(kind=DP), parameter :: fzzinv = 1.0_DP / (8.0_DP * gaminv / 9.0_DP)

    tol: if (den > dentol .and. den1 > 0.0_DP .and. &
         den2 > 0.0_DP) then

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vectors

       zeta = (den1 - den2) / den   ! Polarisation
       zeta = max(zeta,zetamin)
       zeta = min(zeta,zetamax)

       ! PW91 correlation:
       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,eu,eurs)
       call xc_pw91_eq10(0.01554535_DP,0.20548_DP,14.1189_DP, &
            6.1977_DP,3.3662_DP,0.62517_DP,rs,sqrtrs,ep,eprs)
       call xc_pw91_eq10(0.0168869_DP,0.11125_DP,10.357_DP, &
            3.6231_DP,0.88026_DP,0.49671_DP,rs,sqrtrs,alfm,alfrsm)

       fzeta = gaminv * &
            ((1.0_DP+zeta)**FTHRD+(1.0_DP-zeta)**FTHRD-2.0_DP)
       zeta3 = zeta*zeta*zeta
       zeta4 = zeta3*zeta

       ecunif = eu*(1.0_DP-fzeta*zeta4) + ep*fzeta*zeta4 - &
            alfm*fzeta*(1.0_DP-zeta4)*fzzinv
       ecrs = eurs*(1.0_DP-fzeta*zeta4) + eprs*fzeta*zeta4 - &
            alfrsm*fzeta*(1.0_DP-zeta4)*fzzinv
       dfdz = gaminv * FTHRD * &
            ((1.0_DP+zeta)**THIRD-(1.0_DP-zeta)**THIRD)
       eczet = (4.0_DP*zeta3*fzeta+zeta4*dfdz) * &
            (ep-eu+alfm*fzzinv)-dfdz*alfm*fzzinv

       ecd = ecunif-THIRD*rs*ecrs
       ecz = eczet * den

       c_energy = den * ecunif
       c_pot1 = ecd+(1.0_DP-zeta)*ecz/den
       c_pot2 = ecd-(1.0_DP+zeta)*ecz/den

    else

       ! qoh: Set all values to zero for very small densities
       c_energy = 0.0_DP
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP

    end if tol

  end subroutine xc_pw92_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbe_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the PBE exchange at a point       !
    ! in the spin unpolarised case.  (As in original paper.)       !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fxpbe,fs             ! PBE exchange variables
    real(kind=DP) :: exunif               ! LDA exchange per particle
    real(kind=DP) :: p0,p1

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: ul = um / uk


    tol: if (den > dentol) then

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       exunif = ax*den**THIRD

       ! qoh: PBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s*s
       p1 = 1.0_DP / p0
       fxpbe = 1.0_DP + uk - uk * p1
       fs = 2.0_DP * um * s * p1 * p1

       x_energy =  den * exunif * fxpbe
       x_pot = FTHRD * exunif*(fxpbe - fs*s)
       x_dfdmgd = 0.5_DP * ax * fs * TTPI23

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_pbe_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_rpbe_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the RPBE exchange at a point      !
    ! in the spin unpolarised case.                                !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    ! Uses  the RPBE exchange enhancement  factor and derivative   !
    !   Hammer, Hansen & Norskov, Phys. Rev. B 59, 7413 (1999)     !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fxpbe,fs             ! PBE exchange variables
    real(kind=DP) :: exunif               ! LDA exchange per particle
    real(kind=DP) :: p0

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: ul = um / uk

    tol: if (den > dentol) then

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       exunif = ax*den**THIRD
       ! qoh: RPBE exchange enhancement factor and derivative
       p0 = exp(-ul*s*s)
       fxpbe = 1.0_DP + uk - uk * p0
       fs = 2.0_DP * um * s * p0

       x_energy =  den * exunif * fxpbe
       x_pot = FTHRD * exunif*(fxpbe - fs*s)
       x_dfdmgd = 0.5_DP * ax * fs * TTPI23

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_rpbe_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_revpbe_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the revPBE exchange at a point    !
    ! in the spin unpolarised case.                                !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    ! Uses the revised PBE exchange enhancement factor and         !
    ! derivative                                                   !
    ! Zhang & Yang, Phys. Rev. Lett. 80, 890 (1998)                !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fxpbe,fs             ! PBE exchange variables
    real(kind=DP) :: exunif               ! LDA exchange per particle
    real(kind=DP) :: p0,p1

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 1.245_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: ul = um / uk

    tol: if (den > dentol) then

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       exunif = ax*den**THIRD
       ! qoh: revPBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s*s
       p1 = 1.0_DP / p0
       fxpbe = 1.0_DP + uk - uk * p1
       fs = 2.0_DP * um * s * p1 * p1

       x_energy =  den * exunif * fxpbe
       x_pot = FTHRD * exunif*(fxpbe - fs*s)
       x_dfdmgd = 0.5_DP * ax * fs * TTPI23

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_revpbe_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_wc_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the WC exchange at a point        !
    ! in the spin unpolarised case.                                !
    ! Z Wu and R E Cohen, Phys. Rev. B 73, 235116 (2006)           !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    ! Adapted for Wu-Cohen functional by Nicholas Hine, 03/05/2010.!
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: s                    ! Dimensionless density gradient
    real(kind=DP) :: fxwc,fs              ! WC exchange variables
    real(kind=DP) :: exunif               ! LDA exchange per particle
    real(kind=DP) :: p1
    real(kind=DP) :: x

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: c  = 0.007937469335162_DP

    tol: if (den > dentol) then

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless gradient

       exunif = ax*den**THIRD

       ! ndmh: WC x function
       x = 10.0_DP/81.0_DP*s*s+(um-10.0_DP/81.0_DP)*s*s*exp(-s*s)+log(1.0_DP+c*s**4)
       ! ndmh: WC exchange enhancement factor and derivative
       ! ndmh: fxwc = 1 + uk - uk/(1+x/uk) (Eqn 4 of Ref)
       p1 = 1.0_DP / (1.0_DP + x/uk)
       fxwc = 1.0_DP + uk - uk * p1
       fs = (2.0_DP*s*(10.0_dp/81.0_dp+(um-10.0_dp/81.0_dp)*(1.0_DP-s*s)* &
            exp(-s*s)+2.0_DP*c*s*s/(1.0_DP+c*s**4)))*p1*p1

       x_energy =  den * exunif * fxwc
       x_pot = FTHRD * exunif*(fxwc - fs*s)
       x_dfdmgd = 0.5_DP * ax * fs * TTPI23

    else if (den > 0.0_DP) then   ! Low density: exchange only
       exunif = ax*den**THIRD
       x_energy =  den * exunif
       x_pot = FTHRD * exunif
       x_dfdmgd = 0.0_DP
    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_wc_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbe0_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the PBE0 exchange at a point      !
    ! in the spin unpolarised case.                                !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    ! Hybrid functional form from:                                 !
    !  Adamo, C., Cossi, Maurizio, and Barone, Vincenzo:           !
    !    J.Mol. Struc. (Theochem) 493, 145 (1999)                  !
    ! PBE0 = PBE + 1/4 (HF_X - PBE_X)                              !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                       !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP), parameter   :: pbe_fac = 0.75_DP ! PBE0 factor (1 - 1/4)

    call xc_pbe_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)
    x_energy = pbe_fac * x_energy
    x_pot    = pbe_fac * x_pot
    x_dfdmgd = pbe_fac * x_dfdmgd

  end subroutine xc_pbe0_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbe_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the PBE correlation at a point    !
    ! in the spin unpolarised case.                                !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy      ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    real(kind=DP) :: rs                   ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs               ! Square root of rs
    real(kind=DP) :: kf                   ! Fermi wave-vectors
    real(kind=DP) :: ks                   ! Thomas-Fermi screening wave-vector
    real(kind=DP) :: t                    ! Dimensionless density gradient
    real(kind=DP) :: ecunif               ! LDA correlation per particle
    real(kind=DP) :: q1,q4                ! LDA correlation variables
    real(kind=DP) :: t2,bt2,q5,q6,q7      ! PBE correlation variables
    real(kind=DP) :: h,h1,ha,ht,pon,b
    real(kind=DP) :: eurs,aec,ecn,ect

    ! Constants
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: FRONPI = 4.0_DP / PI

    ! PBE correlation constants
    real(kind=dp), parameter :: gamma = 0.03109069086965489503494086371273_DP
    real(kind=dp), parameter :: beta = 0.06672455060314922_DP
    real(kind=dp), parameter :: delta = beta / gamma

    tol: if (den > dentol) then

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       ks = sqrt(FRONPI * kf)       ! Thomas-Fermi wave-vector
       t = mgd / (2.0_DP*den*ks)    ! Dimensionless gradient

       ! PBE correlation:
       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,ecunif,eurs)
       pon = -ecunif/gamma
       q1 = exp(pon)
       b = delta/(q1-1.0_DP)
       t2 = t*t
       bt2 = b*t2
       q4 = 1.0_DP+bt2
       q5 = 1.0_DP+bt2+bt2*bt2
       h = gamma*log(1.0_DP+delta*q4*t2/q5)
       ect = ecunif + h
       c_energy = den * ect
       q6 = 1.0_DP/((gamma*q5+beta*t2*q4)*q5)
       q7 = gamma*q6
       ha = -beta*t2*t2*t2*b*(1.0_DP+q4)*q7
       aec = b*b*q1/beta
       ht = 2.0_DP*beta*t*(1.0_DP+2.0_DP*bt2)*q7
       h1 = (rs*ha*aec*eurs+3.5_DP*t*ht)*THIRD
       ecn = THIRD*rs*eurs
       c_pot = ect - (ecn + h1)
       c_dfdmgd = 0.5_DP * ht / ks

    else
       c_energy = 0.0_DP
       c_pot = 0.0_DP
       c_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_pbe_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbe_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the PBE exchange at a point       !
    ! in the spin polarised case.  (As in original paper.)         !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: kf1,kf2              ! Fermi wave-vectors
    real(kind=DP) :: s1,s2                ! Dimensionless density gradients
    real(kind=DP) :: fs1,fs2              ! PBE exchange variables
    real(kind=DP) :: fxpbe1,fxpbe2
    real(kind=DP) :: exunif1,exunif2      ! LDA exchange per particle
    real(kind=DP) :: p0,p1

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: SXPISQ = 6.0_DP * PI * PI

    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: ul = um / uk

    ! qoh: Initialise values to zero to deal with small densities easily
    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    !qoh: Spin 1:
    tol1: if (den1 > dentol*0.5_DP) then
       kf1 = (SXPISQ * den1)**THIRD

       s1 = mgd1 / (den1*kf1*2.0_DP) ! Dimensionless gradient

       exunif1 = ax*(2.0_DP*den1)**THIRD
       ! qoh: PBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s1*s1
       p1 = 1.0_DP / p0
       fxpbe1 = 1.0_DP + uk - uk * p1
       fs1 = 2.0_DP * um * s1 * p1 * p1

       x_energy = den1 * exunif1 * fxpbe1
       x_pot1 = FTHRD*exunif1*(fxpbe1 - fs1*s1)
       x_dfdmgd1 = 0.5_DP * ax * fs1 * TTPI23

    end if tol1

    !qoh:  Spin 2:
    tol2: if (den2 > 0.5_DP*dentol) then

       kf2 = (SXPISQ * den2)**THIRD
       s2 = mgd2 / (den2*kf2*2.0_DP) !Dimensionless gradient
       exunif2 = ax*(2.0_DP*den2)**THIRD
       ! qoh: PBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s2*s2
       p1 = 1.0_DP / p0
       fxpbe2 = 1.0_DP + uk - uk * p1
       fs2 = 2.0_DP * um * s2 * p1 * p1

       x_energy = x_energy + den2 * exunif2 * fxpbe2
       x_pot2 = FTHRD*exunif2*(fxpbe2 - fs2*s2)
       x_dfdmgd2 = 0.5_DP * ax * fs2 * TTPI23

    end if tol2

  end subroutine xc_pbe_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_rpbe_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the RPBE exchange at a point      !
    ! in the spin polarised case.                                  !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    ! Uses  the RPBE exchange enhancement  factor and derivative   !
    !   Hammer, Hansen & Norskov, Phys. Rev. B 59, 7413 (1999)     !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: kf1,kf2              ! Fermi wave-vectors
    real(kind=DP) :: s1,s2                ! Dimensionless density gradients
    real(kind=DP) :: fs1,fs2              ! PBE exchange variables
    real(kind=DP) :: fxpbe1,fxpbe2
    real(kind=DP) :: exunif1,exunif2      ! LDA exchange per particle
    real(kind=DP) :: p0

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: SXPISQ = 6.0_DP * PI * PI

    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: ul = um / uk

    ! qoh: Initialise values to zero to deal with small densities easily
    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    !qoh: Spin 1:
    tol1: if (den1 > dentol*0.5_DP) then
       kf1 = (SXPISQ * den1)**THIRD

       s1 = mgd1 / (den1*kf1*2.0_DP) ! Dimensionless gradient

       exunif1 = ax*(2.0_DP*den1)**THIRD
       ! qoh: RPBE exchange enhancement factor and derivative
       p0 = exp(-ul*s1*s1)
       fxpbe1 = 1.0_DP + uk - uk * p0
       fs1 = 2.0_DP * um * s1 * p0

       x_energy = den1 * exunif1 * fxpbe1
       x_pot1 = FTHRD*exunif1*(fxpbe1 - fs1*s1)
       x_dfdmgd1 = 0.5_DP * ax * fs1 * TTPI23

    end if tol1

    !qoh:  Spin 2:
    tol2: if (den2 > 0.5_DP*dentol) then

       kf2 = (SXPISQ * den2)**THIRD
       s2 = mgd2 / (den2*kf2*2.0_DP) !Dimensionless gradient
       exunif2 = ax*(2.0_DP*den2)**THIRD
       ! qoh: RPBE exchange enhancement factor and derivative
       p0 = exp(-ul*s2*s2)
       fxpbe2 = 1.0_DP + uk - uk * p0
       fs2 = 2.0_DP * um * s2 * p0

       x_energy = x_energy + den2 * exunif2 * fxpbe2
       x_pot2 = FTHRD*exunif2*(fxpbe2 - fs2*s2)
       x_dfdmgd2 = 0.5_DP * ax * fs2 * TTPI23

    end if tol2

  end subroutine xc_rpbe_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_revpbe_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the revPBE exchange at a point    !
    ! in the spin polarised case.                                  !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    ! Uses the revised PBE exchange enhancement factor and         !
    ! derivative                                                   !
    ! Zhang & Yang, Phys. Rev. Lett. 80, 890 (1998)                !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: kf1,kf2              ! Fermi wave-vectors
    real(kind=DP) :: s1,s2                ! Dimensionless density gradients
    real(kind=DP) :: fs1,fs2              ! PBE exchange variables
    real(kind=DP) :: fxpbe1,fxpbe2
    real(kind=DP) :: exunif1,exunif2      ! LDA exchange per particle
    real(kind=DP) :: p0,p1

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: SXPISQ = 6.0_DP * PI * PI

    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 1.245_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: ul = um / uk

    ! qoh: Initialise values to zero to deal with small densities easily
    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    !qoh: Spin 1:
    tol1: if (den1 > dentol*0.5_DP) then
       kf1 = (SXPISQ * den1)**THIRD

       s1 = mgd1 / (den1*kf1*2.0_DP) ! Dimensionless gradient

       exunif1 = ax*(2.0_DP*den1)**THIRD
       ! qoh: revPBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s1*s1
       p1 = 1.0_DP / p0
       fxpbe1 = 1.0_DP + uk - uk * p1
       fs1 = 2.0_DP * um * s1 * p1 * p1

       x_energy = den1 * exunif1 * fxpbe1
       x_pot1 = FTHRD*exunif1*(fxpbe1 - fs1*s1)
       x_dfdmgd1 = 0.5_DP * ax * fs1 * TTPI23

    end if tol1

    !qoh:  Spin 2:
    tol2: if (den2 > 0.5_DP*dentol) then

       kf2 = (SXPISQ * den2)**THIRD
       s2 = mgd2 / (den2*kf2*2.0_DP) !Dimensionless gradient
       exunif2 = ax*(2.0_DP*den2)**THIRD
       ! qoh: revPBE exchange enhancement factor and derivative
       p0 = 1.0_DP + ul*s2*s2
       p1 = 1.0_DP / p0
       fxpbe2 = 1.0_DP + uk - uk * p1
       fs2 = 2.0_DP * um * s2 * p1 * p1

       x_energy = x_energy + den2 * exunif2 * fxpbe2
       x_pot2 = FTHRD*exunif2*(fxpbe2 - fs2*s2)
       x_dfdmgd2 = 0.5_DP * ax * fs2 * TTPI23

    end if tol2

  end subroutine xc_revpbe_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_wc_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the WC exchange at a point        !
    ! in the spin polarised case.                                  !
    ! Z Wu and R E Cohen, Phys. Rev. B 73, 235116 (2006)           !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    ! Adapted for Wu-Cohen functional by Nicholas Hine, 03/05/2010.!
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: kf1,kf2              ! Fermi wave-vectors
    real(kind=DP) :: s1,s2                ! Dimensionless density gradients
    real(kind=DP) :: fs1,fs2              ! PBE exchange variables
    real(kind=DP) :: fxwc1,fxwc2
    real(kind=DP) :: exunif1,exunif2      ! LDA exchange per particle
    real(kind=DP) :: p1
    real(kind=DP) :: x1,x2

    ! Constants
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP
    real(kind=DP), parameter :: SXPISQ = 6.0_DP * PI * PI

    ! Exchange constants
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: uk = 0.8040_DP
    real(kind=DP), parameter :: um = 0.2195149727645171_DP
    real(kind=DP), parameter :: c  = 0.007937469335162_DP

    ! qoh: Initialise values to zero to deal with small densities easily
    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    ! ndmh: Spin 1:
    tol1: if (den1 > dentol*0.5_DP) then
       kf1 = (SXPISQ * den1)**THIRD

       s1 = mgd1 / (den1*kf1*2.0_DP) ! Dimensionless gradient
       exunif1 = ax*(2.0_DP*den1)**THIRD

       ! ndmh: WC x function
       x1 = 10.0_DP/81.0_DP*s1*s1+(um-10.0_DP/81.0_DP)*s1*s1*exp(-s1*s1)+log(1+c*s1**4)
       ! ndmh: WC exchange enhancement factor and derivative
       ! ndmh: fxwc = 1 + uk - uk/(1+x/uk) (Eqn 4 of Ref)
       p1 = 1.0_DP / (1.0_DP + x1/uk)
       fxwc1 = 1.0_DP + uk - uk * p1
       fs1 = (2.0_DP*s1*(10.0_dp/81.0_dp+(um-10.0_dp/81.0_dp)*(1.0_DP-s1*s1)* &
            exp(-s1*s1)+2.0_DP*c*s1*s1/(1.0_DP+c*s1**4)))*p1*p1

       x_energy = den1 * exunif1 * fxwc1
       x_pot1 = FTHRD*exunif1*(fxwc1 - fs1*s1)
       x_dfdmgd1 = 0.5_DP * ax * fs1 * TTPI23

    end if tol1

    ! ndmh:  Spin 2:
    tol2: if (den2 > 0.5_DP*dentol) then

       kf2 = (SXPISQ * den2)**THIRD
       s2 = mgd2 / (den2*kf2*2.0_DP) !Dimensionless gradient
       exunif2 = ax*(2.0_DP*den2)**THIRD

       ! ndmh: WC x function
       x2 = 10.0_DP/81.0_DP*s2*s2+(um-10.0_DP/81.0_DP)*s2*s2*exp(-s2*s2)+log(1+c*s2**4)
       ! ndmh: WC exchange enhancement factor and derivative
       ! ndmh: fxwc = 1 + uk - uk/(1+x/uk) (Eqn 4 of Ref)
       p1 = 1.0_DP / (1.0_DP + x2/uk)
       fxwc2 = 1.0_DP + uk - uk * p1
       fs2 = (2.0_DP*s2*(10.0_dp/81.0_dp+(um-10.0_dp/81.0_dp)*(1.0_DP-s2*s2)* &
            exp(-s2*s2)+2.0_DP*c*s2*s2/(1.0_DP+c*s2**4)))*p1*p1

       x_energy = x_energy + den2 * exunif2 * fxwc2
       x_pot2 = FTHRD*exunif2*(fxwc2 - fs2*s2)
       x_dfdmgd2 = 0.5_DP * ax * fs2 * TTPI23

    end if tol2

  end subroutine xc_wc_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbe0_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the PBE0 exchange at a point      !
    ! in the spin polarised case.                                  !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    ! Hybrid functional form from:                                 !
    !  Adamo, C., Cossi, Maurizio, and Barone, Vincenzo:           !
    !    J.Mol. Struc. (Theochem) 493, 145 (1999)                  !
    ! PBE0 = PBE + 1/4 (HF_X - PBE_X)                              !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                       !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP), parameter   :: xfac = 0.75_DP ! PBE0 factor (1 - 1/4)

    call xc_pbe_x_point_sp(den1,den2,mgd1,mgd2,&
         x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)
    x_energy  = xfac * x_energy
    x_pot1    = xfac * x_pot1
    x_pot2    = xfac * x_pot2
    x_dfdmgd1 = xfac * x_dfdmgd1
    x_dfdmgd2 = xfac * x_dfdmgd2

  end subroutine xc_pbe0_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pbe_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==============================================================!
    ! This subroutine calculates the PBE correlation at a point    !
    ! in the spin polarised case.                                  !
    !  J P Perdew, K Burke & M Ernzerhof,                          !
    !    Phys. Rev. Lett. 77, 3865 (1996)                          !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 11/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|

    real(kind=DP) :: rs                   ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs               ! Square root of rs
    real(kind=DP) :: kf                   ! Fermi wave-vectors
    real(kind=DP) :: ks                   ! Thomas-Fermi screening wave-vector
    real(kind=DP) :: t                    ! Dimensionless density gradient
    real(kind=DP) :: ecunif               ! LDA correlation per particle
    real(kind=DP) :: q1,q2,q3,q4          ! LDA correlation variables
    real(kind=DP) :: t2,q5,q6,q7,q8       ! PBE correlation variables
    real(kind=DP) :: h,ht,pon,b,hb
    real(kind=DP) :: eurs
    real(kind=DP) :: eprs,ep,alfm,alfrsm
    real(kind=DP) :: ecd,ecz,eu,ecrs
    real(kind=DP) :: ec1d,ec1z,eczet
    real(kind=DP) :: g2,g3,g4,be,bg,gz
    real(kind=DP) :: b2,t3,t4,g
    real(kind=DP) :: zeta                 ! Spin polarisation
    real(kind=DP) :: zeta3,zeta4          ! zeta**3 and zeta**4
    real(kind=DP) :: fzeta,dfdz           ! Function of zeta and derivative
    real(kind=DP) :: mgd12                ! mgd1 + mgd2 to avoid compiler warning

    ! Constants
    real(kind=DP), parameter :: SSXTH = 7.0_DP / 6.0_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: FRONPI = 4.0_DP / PI
    real(kind=DP), parameter :: TWTHRD = 2.0_DP / 3.0_DP

    ! PBE correlation constants
    real(kind=dp), parameter :: gamma = 0.03109069086965489503494086371273_DP
    real(kind=dp), parameter :: beta = 0.06672455060314922_DP
    real(kind=dp), parameter :: delta = beta / gamma
    real(kind=DP), parameter :: zetamin = -1.0_DP + 1.0e-10_DP
    real(kind=DP), parameter :: zetamax = 1.0_DP - 1.0e-10_DP

    ! PW91 correlation parameters
    real(kind=DP), parameter :: gaminv = 1.92366105093153631981_DP
    real(kind=DP), parameter :: fzzinv = 1.0_DP / (8.0_DP * gaminv / 9.0_DP)

    ! qoh: Avoid compiler warning by making use of mgd1 and mgd2
    mgd12 = mgd1 + mgd2

    tol: if (den > dentol .and. den1 > 0.0_DP .and. den2 > 0.0_DP) then

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vectors
       ks = sqrt(FRONPI * kf)       ! Thomas-Fermi wave-vector

       zeta = (den1 - den2) / den   ! Polarisation
       zeta = max(zeta,zetamin)
       zeta = min(zeta,zetamax)
       g = 0.5_DP*((1.0_DP-zeta)**TWTHRD+(1.0_DP+zeta)**TWTHRD)

       t = mgd / (den*ks*g*2.0_DP)! Dimensionless gradient

       ! PBE correlation:
       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,eu,eurs)
       call xc_pw91_eq10(0.01554535_DP,0.20548_DP, &
            14.1189_DP,6.1977_DP,3.3662_DP,0.62517_DP,rs,sqrtrs, &
            ep,eprs)
       call xc_pw91_eq10(0.0168869_DP,0.11125_DP,10.357_DP, &
            3.6231_DP,0.88026_DP,0.49671_DP,rs,sqrtrs, &
            alfm,alfrsm)

       fzeta = gaminv * &
            ((1.0_DP+zeta)**FTHRD+(1.0_DP-zeta)**FTHRD-2.0_DP)
       zeta3 = zeta*zeta*zeta
       zeta4 = zeta3*zeta

       ecunif = eu*(1.0_DP-fzeta*zeta4) + ep*fzeta*zeta4 - &
            alfm*fzeta*(1.0_DP-zeta4)*fzzinv
       ecrs = eurs*(1.0_DP-fzeta*zeta4) + eprs*fzeta*zeta4 - &
            alfrsm*fzeta*(1.0_DP-zeta4)*fzzinv
       dfdz = gaminv * FTHRD * &
            ((1.0_DP+zeta)**THIRD-(1.0_DP-zeta)**THIRD)
       eczet = (4.0_DP*zeta3*fzeta+zeta4*dfdz) * &
            (ep-eu+alfm*fzzinv)-dfdz*alfm*fzzinv
       ecd = ecunif - THIRD*rs*ecrs
       ecz = eczet * den

       g2 = g*g
       g3 = g2*g
       g4 = g2*g2
       pon = -ecunif/(g3*gamma)
       ep = exp(pon)
       b = delta/(ep-1.0_DP)
       b2 = b*b
       t2 = t*t
       t3 = t2*t
       t4 = t2*t2
       q1 = t2 + b*t4
       q2 = 1.0_DP+b*t2+b2*t4
       q3 = q1/q2
       q4 = q2*q2 + delta*q1*q2
       q5 = beta*g3/q4
       q6 = 2.0_DP*t+4.0_DP*b*t3
       q7 = -b*t4*t2*(2.0_DP+b*t2)
       q8 = ep*b2/beta
       h = gamma*g3*log(1.0_DP + delta*q3)
       ht = q5*q6
       hb = q5*q7
       be = q8/g3
       bg = -3.0_DP*ecunif*q8/g4
       gz = ((1.0_DP+zeta)**(-THIRD)-(1.0_DP-zeta)**(-THIRD))*THIRD
       ec1d = h - THIRD*rs*ecrs*be*hb - SSXTH*t*ht
       ec1z = den * (3.0_DP*h*gz/g + hb*(bg*gz + be*eczet) - ht*t/g*gz)

       c_energy = den * (ecunif + h)
       c_pot1 = ecd+ec1d+(1.0_DP-zeta)*(ecz+ec1z)/den
       c_pot2 = ecd+ec1d-(1.0_DP+zeta)*(ecz+ec1z)/den

       ! qoh: Correlation gradient dependent part
       ! ndmh: separated c_dfdmgd from c_dfdmgd1
       c_dfdmgd = 0.5_DP * ht / (g*ks)
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP
    else
       c_energy = 0.0_DP ! ddor
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP
       c_dfdmgd = 0.0_DP
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP
    end if tol

  end subroutine xc_pbe_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw91_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the PW91 exchange at a point      !
    ! in the spin unpolarised case.                                !
    !    J P Perdew & Y Wang, Phys. Rev. B 45, 13244 (1992)        !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution  by Quintin Hill on 10/03/2009.!
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: kf                   ! Fermi wave-vectors
    real(kind=DP) :: s                    ! Dimensionless density gradients
    real(kind=DP) :: fac,ss,s3,s4,f,fs    ! PW91 exchange variables
    real(kind=DP) :: p0,p1,p2,p3,p4,p5,p6

    ! Constants
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP

    ! PW91 exchange constants
    real(kind=DP), parameter :: a1 = 0.19645_DP
    real(kind=DP), parameter :: a2 = 0.27430_DP
    real(kind=DP), parameter :: a3 = 0.15084_DP
    real(kind=DP), parameter :: a4 = 100.0_DP
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: a = 7.7956_DP
    real(kind=DP), parameter :: b1 = 0.004_DP

    tol: if (den > dentol) then   ! Positive charge density

       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       s = mgd / (2.0_DP*den*kf)    ! Dimensionless

       fac = ax*den**THIRD
       ss = s*s
       s3 = ss*s
       s4 = ss*ss
       p0 = 1.0_DP/sqrt(1.0_DP+a*a*ss)
       p1 = log(a*s+1.0_DP/p0)
       p2 = exp(-a4*ss)
       p3 = 1.0_DP/(1.0_DP+a1*s*p1+b1*s4)
       p4 = 1.0_DP+a1*s*p1+(a2-a3*p2)*ss
       f = p3*p4
       x_energy = den * fac * f

       p5 = 2.0_DP*(s*(a2-a3*p2)+a3*a4*s3*p2-2.0_DP*b1*s3)
       p6 = (a1*(p1+a*s*p0)+4.0_DP*b1*s3)*((a2-a3*p2)*ss-b1*s4)
       fs = (p5*p3-p6*p3*p3)
       x_pot = FTHRD*fac*(f-s*fs)

       x_dfdmgd = 0.5_DP * ax * fs * TTPI23

    else if (den > 0.0_DP) then   ! Low density: exchange only

       fac = ax*den**THIRD
       x_energy = den * fac
       x_pot = FTHRD * fac
       x_dfdmgd = 0.0_DP

    else
       x_energy = 0.0_DP
       x_pot = 0.0_DP
       x_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_pw91_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw91_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the PW91 correlation at a point   !
    ! in the spin unpolarised case.                                !
    !    J P Perdew & Y Wang, Phys. Rev. B 45, 13244 (1992)        !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 10/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy      ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    real(kind=DP) :: rs                   ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs               ! Square root of rs
    real(kind=DP) :: kf                   ! Fermi wave-vectors
    real(kind=DP) :: ks                   ! Thomas-Fermi screening wave-vector
    real(kind=DP) :: t                    ! Dimensionless density gradient
    real(kind=DP) :: ecunif               ! LDA correlation energy per particle
    real(kind=DP) :: ecd,eurs,ec1d        ! PW91 correlation variables
    real(kind=DP) :: bet,delt,pon,b,b2,t2,t4,t6,rs2
    real(kind=DP) :: rs3,q4,q5,q6,q7,cc,r0,r1,coeff,r2,r3
    real(kind=DP) :: h0,h1,h,q8,h0t,h0b,h0rs,h1t,ccrs,r1rs
    real(kind=DP) :: h1rs,ht,hrs!,g,h0z,h1z,hz,gz,g3,g4

    ! Constants
    real(kind=DP), parameter :: SSXTH = 7.0_DP / 6.0_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: FRONPI = 4.0_DP / PI
    real(kind=DP), parameter :: THSVTH = 3.0_DP / 7.0_DP

    ! PW91 correlation parameters
    real(kind=DP), parameter :: a4 = 100.0_DP
    real(kind=DP), parameter :: xnu = 15.75592_dp
    real(kind=DP), parameter :: cc0 = 0.004235_dp
    real(kind=DP), parameter :: cx = -0.001667212_dp
    real(kind=DP), parameter :: alf = 0.09_dp
    real(kind=DP), parameter :: c1 = 0.002568_dp
    real(kind=DP), parameter :: c2 = 0.023266_dp
    real(kind=DP), parameter :: c3 = 0.000007389_dp
    real(kind=DP), parameter :: c4 = 8.723_dp
    real(kind=DP), parameter :: c5 = 0.472_dp
    real(kind=DP), parameter :: c6 = 0.07389_dp

    tol: if (den > dentol) then   ! Positive charge density
       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vector
       ks = sqrt(FRONPI * kf)       ! Thomas-Fermi wave-vector

       t = mgd / (2.0_DP*den*ks)    ! gradients

       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,ecunif,eurs)
       ecd = ecunif-THIRD*rs*eurs
       bet = xnu*cc0
       delt = 2.0_DP*alf/bet
       pon = -delt*ecunif/bet
       b = delt/(exp(pon)-1.0_DP)
       b2 = b*b
       t2 = t*t
       t4 = t2*t2
       t6 = t4*t2
       rs2 = rs*rs
       rs3 = rs2*rs
       q4 = 1.0_DP+b*t2
       q5 = 1.0_DP+b*t2+b2*t4
       q6 = c1+c2*rs+c3*rs2
       q7 = 1.0_DP+c4*rs+c5*rs2+c6*rs3
       cc = -cx + q6/q7
       r0 = ks*ks/(kf*kf)
       r1 = a4*r0
       coeff = cc-cc0-THSVTH*cx
       r2 = xnu*coeff
       r3 = exp(-r1*t2)
       h0 = (bet/delt)*log(1.0_DP+delt*q4*t2/q5)
       h1 = r3*r2*t2
       h = h0+h1
       q8 = q5*q5+delt*q4*q5*t2
       h0t = 2.0_DP*bet*t*(1.0_DP+2.0_DP*b*t2)/q8
       h0b = -bet*t6*(2.0_DP*b+b2*t2)/q8
       h0rs = h0b*b*eurs*(b+delt)/bet
       h1t = 2.0_DP*r3*r2*t*(1.0_DP-r1*t2)
       ccrs = (c2+2.0_DP*c3*rs)/q7 - &
            q6*(c4+2.0_DP*c5*rs+3.0_DP*c6*rs2)/(q7*q7)
       r1rs = 100.0_DP*r0/rs
       h1rs = xnu*t2*r3*(ccrs - coeff*t2*r1rs)
       ht = h0t + h1t
       hrs = h0rs + h1rs
       ec1d = h-THIRD*rs*hrs-SSXTH*t*ht

       c_energy = den * (ecunif + h)
       c_pot = ecd + ec1d
       c_dfdmgd = 0.5_DP * ht / ks

    else
       c_energy = 0.0_DP
       c_pot = 0.0_DP
       c_dfdmgd = 0.0_DP
    end if tol

  end subroutine xc_pw91_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw91_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the PW91 exchange at a point      !
    ! in the spin polarised case.                                  !
    !    J P Perdew & Y Wang, Phys. Rev. B 45, 13244 (1992)        !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 10/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: kf1,kf2              ! Fermi wave-vectors
    real(kind=DP) :: s1,s2                ! Dimensionless density gradients
    real(kind=DP) :: fac,ss,s3,s4,f,fs    ! PW91 exchange variables
    real(kind=DP) :: p0,p1,p2,p3,p4,p5,p6

    ! Constants
    real(kind=DP), parameter :: SXPISQ = 6.0_DP * PI * PI
    real(kind=DP), parameter :: TTPI23 = 0.32324091934799096266_DP

    ! PW91 exchange constants
    real(kind=DP), parameter :: a1 = 0.19645_DP
    real(kind=DP), parameter :: a2 = 0.27430_DP
    real(kind=DP), parameter :: a3 = 0.15084_DP
    real(kind=DP), parameter :: a4 = 100.0_DP
    real(kind=DP), parameter :: ax = -0.738558766382022405884230032680836_DP
    real(kind=DP), parameter :: a = 7.7956_DP
    real(kind=DP), parameter :: b1 = 0.004_DP

    ! qoh: Initialise values to zero to deal with small densities easily
    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    tol1: if (den1 > dentol*0.5_DP) then
       ! Spin 1:
       kf1 = (SXPISQ * den1)**THIRD
       s1 = mgd1 / (den1*kf1*2.0_DP) ! Dimensionless gradients
       fac = ax*(2.0_DP*den1)**THIRD
       ss = s1*s1
       s3 = ss*s1
       s4 = ss*ss
       p0 = 1.0_DP/sqrt(1.0_DP+a*a*ss)
       p1 = log(a*s1+1.0_DP/p0)
       p2 = exp(-a4*ss)
       p3 = 1.0_DP/(1.0_DP+a1*s1*p1+b1*s4)
       p4 = 1.0_DP+a1*s1*p1+(a2-a3*p2)*ss
       f = p3*p4
       x_energy = den1 * fac * f
       p5 = 2.0_DP*(s1*(a2-a3*p2)+a3*a4*s3*p2-2.0_DP*b1*s3)
       p6 = (a1*(p1+a*s1*p0)+4.0_DP*b1*s3)*((a2-a3*p2)*ss-b1*s4)
       fs = (p5*p3-p6*p3*p3)
       x_pot1 = FTHRD * fac * (f - s1*fs)
       x_dfdmgd1 = 0.5_DP * ax * fs * TTPI23

    end if tol1

    tol2: if (den2 > 0.5_DP*dentol) then

       ! Spin 2:
       kf2 = (SXPISQ * den2)**THIRD
       s2 = mgd2 / (den2*kf2*2.0_DP) ! Dimensionless gradients
       fac = ax*(2.0_DP*den2)**THIRD
       ss = s2*s2
       s3 = ss*s2
       s4 = ss*ss
       p0 = 1.0_DP/sqrt(1.0_DP+a*a*ss)
       p1 = log(a*s2+1.0_DP/p0)
       p2 = exp(-a4*ss)
       p3 = 1.0_DP/(1.0_DP+a1*s2*p1+b1*s4)
       p4 = 1.0_DP+a1*s2*p1+(a2-a3*p2)*ss
       f = p3*p4
       x_energy = x_energy + den2 * fac * f
       p5 = 2.0_DP*(s2*(a2-a3*p2)+a3*a4*s3*p2-2.0_DP*b1*s3)
       p6 = (a1*(p1+a*s2*p0)+4.0_DP*b1*s3)*((a2-a3*p2)*ss-b1*s4)
       fs = (p5*p3-p6*p3*p3)
       x_pot2 = FTHRD * fac * (f - s2*fs)
       x_dfdmgd2 = 0.5_DP * ax * fs * TTPI23

    end if tol2

  end subroutine xc_pw91_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw91_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==============================================================!
    ! This subroutine calculates the PW91 correlation at a point   !
    ! in the spin polarised case.                                  !
    !    J P Perdew & Y Wang, Phys. Rev. B 45, 13244 (1992)        !
    !--------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.           !
    ! Split into point contribution by Quintin Hill on 10/03/2009. !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad_n_2|

    real(kind=DP) :: rs                   ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs               ! Square root of rs
    real(kind=DP) :: kf                   ! Fermi wave-vector
    real(kind=DP) :: ks                   ! Thomas-Fermi screening wave-vector
    real(kind=DP) :: t                    ! Dimensionless density gradient
    real(kind=DP) :: ecunif               ! LDA correlation energy per particle
    real(kind=DP) :: ecd,eu,eurs,ec1d     ! PW91 correlation variables
    real(kind=DP) :: bet,delt,pon,b,b2,t2,t4,t6,rs2
    real(kind=DP) :: rs3,q4,q5,q6,q7,cc,r0,r1,coeff,r2,r3
    real(kind=DP) :: h0,h1,h,q8,h0t,h0b,h0rs,h1t,ccrs,r1rs
    real(kind=DP) :: h1rs,ht,hrs,g,h0z,h1z,hz,gz,g3,g4
    real(kind=DP) :: ecz,ecrs,eczet,ec1z
    real(kind=DP) :: ep,eprs,alfm,alfrsm
    real(kind=DP) :: zeta                 ! Spin polarisation
    real(kind=DP) :: zeta3,zeta4          ! zeta**3 and zeta**4
    real(kind=DP) :: fzeta,dfdz           ! Function of zeta and derivative
    real(kind=DP) :: mgd12                ! mgd1 + mgd 2 avoids compiler warning

    ! Constants
    real(kind=DP), parameter :: SSXTH = 7.0_DP / 6.0_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: THPISQ = 3.0_DP * PI * PI
    real(kind=DP), parameter :: FRONPI = 4.0_DP / PI
    real(kind=DP), parameter :: THSVTH = 3.0_DP / 7.0_DP
    real(kind=DP), parameter :: TWTHRD = 2.0_DP / 3.0_DP

    ! PW91 correlation parameters
    real(kind=DP), parameter :: a4 = 100.0_DP
    real(kind=DP), parameter :: xnu = 15.75592_dp
    real(kind=DP), parameter :: cc0 = 0.004235_dp
    real(kind=DP), parameter :: cx = -0.001667212_dp
    real(kind=DP), parameter :: alf = 0.09_dp
    real(kind=DP), parameter :: c1 = 0.002568_dp
    real(kind=DP), parameter :: c2 = 0.023266_dp
    real(kind=DP), parameter :: c3 = 0.000007389_dp
    real(kind=DP), parameter :: c4 = 8.723_dp
    real(kind=DP), parameter :: c5 = 0.472_dp
    real(kind=DP), parameter :: c6 = 0.07389_dp
    real(kind=DP), parameter :: gaminv = 1.92366105093153631981_DP
    real(kind=DP), parameter :: fzzinv = 1.0_DP / (8.0_DP * gaminv / 9.0_DP)
    real(kind=DP), parameter :: zetamin = -1.0_DP + 1.0e-10_DP
    real(kind=DP), parameter :: zetamax = 1.0_DP - 1.0e-10_DP

    ! qoh: Avoid compiler warning by making use of mgd1 and mgd2
    mgd12 = mgd1 + mgd2

    tol: if (den > dentol .and. den1 > 0.0_DP .and. &
         den2 > 0.0_DP) then

       rs = (TONFPI / den)**THIRD   ! Wigner-Seitz radius
       kf = (THPISQ * den)**THIRD   ! Fermi wave-vectors
       ks = sqrt(FRONPI * kf)       ! Thomas-Fermi wave-vector

       zeta = (den1 - den2) / den   ! Polarisation
       zeta = max(zeta,zetamin)
       zeta = min(zeta,zetamax)
       g = 0.5_DP*((1.0_dp-zeta)**TWTHRD+(1.0_DP+zeta)**TWTHRD)

       t = mgd / (den*ks*g*2.0_DP) ! Dimensionless gradient

       ! PW91 correlation:
       sqrtrs = sqrt(rs)
       call xc_pw91_eq10(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,eu,eurs)
       call xc_pw91_eq10(0.01554535_DP,0.20548_DP,14.1189_DP, &
            6.1977_DP,3.3662_DP,0.62517_DP,rs,sqrtrs,ep,eprs)
       call xc_pw91_eq10(0.0168869_DP,0.11125_DP,10.357_DP, &
            3.6231_DP,0.88026_DP,0.49671_DP,rs,sqrtrs,alfm,alfrsm)

       fzeta = gaminv * &
            ((1.0_DP+zeta)**FTHRD+(1.0_DP-zeta)**FTHRD-2.0_DP)
       zeta3 = zeta*zeta*zeta
       zeta4 = zeta3*zeta

       ecunif = eu*(1.0_DP-fzeta*zeta4) + ep*fzeta*zeta4 - &
            alfm*fzeta*(1.0_DP-zeta4)*fzzinv
       ecrs = eurs*(1.0_DP-fzeta*zeta4) + eprs*fzeta*zeta4 - &
            alfrsm*fzeta*(1.0_DP-zeta4)*fzzinv
       dfdz = gaminv * FTHRD * &
            ((1.0_DP+zeta)**THIRD-(1.0_DP-zeta)**THIRD)
       eczet = (4.0_DP*zeta3*fzeta+zeta4*dfdz) * &
            (ep-eu+alfm*fzzinv)-dfdz*alfm*fzzinv

       ecd = ecunif-THIRD*rs*ecrs
       ecz = eczet * den

       bet = xnu*cc0
       delt = 2.0_DP*alf/bet
       gz = ((1.0_DP+zeta)**(-THIRD)-(1.0_DP-zeta)**(-THIRD))*THIRD
       g3 = g*g*g
       g4 = g3*g
       pon = -delt*ecunif/(g3*bet)
       b = delt/(exp(pon)-1.0_DP)
       b2 = b*b
       t2 = t*t
       t4 = t2*t2
       t6 = t4*t2
       rs2 = rs*rs
       rs3 = rs2*rs
       q4 = 1.0_DP+b*t2
       q5 = 1.0_DP+b*t2+b2*t4
       q6 = c1+c2*rs+c3*rs2
       q7 = 1.0_DP+c4*rs+c5*rs2+c6*rs3
       cc = -cx + q6/q7
       r0 = ks*ks/(kf*kf)
       r1 = a4*r0*g4
       coeff = cc-cc0-THSVTH*cx
       r2 = xnu*coeff*g3
       r3 = exp(-r1*t2)
       h0 = g3*(bet/delt)*log(1.0_DP+delt*q4*t2/q5)
       h1 = r3*r2*t2
       h = h0+h1
       q8 = q5*q5+delt*q4*q5*t2
       h0t = 2.0_DP*bet*t*g3*(1.0_DP+2.0_DP*b*t2)/q8
       h0b = -bet*t6*g3*(2.0_DP*b+b2*t2)/q8
       h0rs = h0b*b*ecrs*(b+delt)/(bet*g3)
       h0z = 3.0_DP*h0*gz/g + &
            h0b*b*(b+delt)*(eczet-3.0_DP*ecunif*gz/g)/(g3*bet)
       h1t = 2.0_DP*r3*r2*t*(1.0_DP-r1*t2)
       ccrs = (c2+2.0_DP*c3*rs)/q7 - &
            q6*(c4+2.0_DP*c5*rs+3.0_DP*c6*rs2)/(q7*q7)
       r1rs = r1/rs
       h1rs = xnu*t2*r3*g3*(ccrs - coeff*t2*r1rs)
       h1z = h1*(3.0_DP-4.0_DP*r1*t2)*gz/g
       ht = h0t + h1t
       hrs = h0rs + h1rs
       hz = h0z + h1z
       ec1d = h-THIRD*rs*hrs-SSXTH*t*ht
       ec1z = den*(hz-ht*t*gz/g)

       c_energy = den * (ecunif + h)
       c_pot1 = ecd+ec1d+(1.0_DP-zeta)*(ecz+ec1z)/den
       c_pot2=  ecd+ec1d-(1.0_DP+zeta)*(ecz+ec1z)/den

       ! qoh: Correlation gradient dependent par
       ! ndmh: separated c_dfdmgd from c_dfdmgd1t
       c_dfdmgd = 0.5_DP * ht / (g*ks)
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP
    else

       ! qoh: Set all values to zero for very small densities
       c_energy = 0.0_DP
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP
       c_dfdmgd = 0.0_DP
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP

    end if tol

  end subroutine xc_pw91_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw91_eq10(A,alpha1,beta1,beta2,beta3,beta4,rs,sqrtrs,G,dGdrs)

    !---------------------------------------------------------!
    ! This subroutine evaluates Eq. 10 of:                    !
    !  J P Perdew & Y Wang, Phys. Rev. B 45, 13244 (1992)     !
    ! for G(rs,A,alpha1,beta1,beta2,beta3,beta4,p=1)          !
    ! and its derivative with respect to rs.                  !
    !---------------------------------------------------------!

    implicit none

    !------------------------
    ! INPUT/OUTPUT variables:
    !------------------------

    real(kind=DP), intent(in)  :: A,alpha1,beta1,beta2,beta3,beta4,rs,sqrtrs
    real(kind=DP), intent(out) :: G,dGdrs

    !------------------------
    ! Local variables:
    !------------------------

    real(kind=DP) :: q0,rs32,rs2,q1,q2,q3

    q0 = -2.0_DP*A*(1.0_DP+alpha1*rs)
    rs32 = sqrtrs*rs
    rs2 = rs*rs
    q1 = 2.0_DP*A*(beta1*sqrtrs+beta2*rs+beta3*rs32+beta4*rs2)
    q2 = log(1.0_DP+1.0_DP/q1)
    G = q0*q2
    q3 = A*(beta1/sqrtrs+2.0_DP*beta2+3.0_DP*beta3*sqrtrs+4.0_DP*beta4*rs)
    dGdrs = -2.0_DP*A*alpha1*q2-q0*q3/(q1*q1+q1)

  end subroutine xc_pw91_eq10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b88_x_point(denin,mgdin,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the Becke 88 exchange at a point  !
    ! in the spin unpolarised case.                                !
    !  A.D. Becke                                                  !
    !  Density-functional exchange-energy approximation with       !
    !  correct asymptotic behaviour                                !
    !  Phys. Rev. A38 (1988) 3098-3100                             !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 19/02/2009.                       !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: denin    ! The charge density at a point
    real(kind=dp), intent(in)  :: mgdin    ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! Exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: den ! Density
    real(kind=DP) :: xsigma ! x_{sigma} in paper
    real(kind=DP) :: mgd, mgdsq ! Modulus of the the density gradient and square
    real(kind=DP) :: crden, densq ! Cube root of den and den squared
    real(kind=DP) :: crdenbq ! den ^ (4/3)
    real(kind=DP) :: asinhx ! Inverse hyperbolic sine of x
    real(kind=DP) :: dasinhx  ! Derivative wrt x of asinh(x)
    real(kind=DP) :: edenom, edenomsq ! Energy denominator and its square
    real(kind=DP) :: dnedenom ! Energy denominator density derivative
    real(kind=DP) :: dgedenom ! Energy denominator density gradient derivative

    ! qoh: Parameters:
    real(kind=DP), parameter :: beta = 0.0042_DP ! As in paper
    real(kind=DP), parameter :: sixbeta = 6.0_DP * beta
    real(kind=DP), parameter :: ldafac =  -0.73855876638202234_DP
    real(kind=DP), parameter :: ftldafac = fthrd*ldafac
    real(kind=DP), parameter :: crtwo = 1.2599210498948732_DP ! spin factor
    real(kind=DP), parameter :: rcrtwo = 1.0_DP / crtwo ! spin factor

    den = max(0.0_dp,denin)
    tol: if(den.gt.dentol) then
       mgd = max(0.0_dp,mgdin)
       crden = den**THIRD
       crdenbq = crden*den
       mgdsq = mgd**2
       xsigma = crtwo * mgd / crdenbq
       asinhx = log(xsigma+sqrt(1.0_DP+xsigma**2))
       edenom = rcrtwo*crdenbq+sixbeta*mgd*asinhx
       x_energy = ldafac * crdenbq - beta * mgdsq /edenom

       edenomsq = edenom**2
       densq = den**2
       dasinhx = 1.0_DP / sqrt(1.0_dp+xsigma**2)
       dnedenom = fthrd *(rcrtwo * crden - sixbeta* mgdsq * dasinhx * crtwo / &
            (crden * densq))
       dgedenom = sixbeta* mgd * dasinhx * crtwo / crdenbq + sixbeta*asinhx

       x_pot = ftldafac * crden + (beta * mgdsq * dnedenom)/ edenomsq
       x_dfdmgd = - 2.0_DP * beta * mgd / edenom &
            + (beta * mgdsq * dgedenom) / edenomsq
    else
       ! qoh: Set all values to zero for very small densities
       x_energy = 0.0_dp
       x_pot = 0.0_dp
       x_dfdmgd = 0.0_dp
    endif tol

  end subroutine xc_b88_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b88_x_point_sp(den1in,den2in,mgd1in,mgd2in,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the Becke 88 exchange at a point  !
    ! in the spin polarised case.                                  !
    !  A.D. Becke                                                  !
    !  Density-functional exchange-energy approximation with       !
    !  correct asymptotic behaviour                                !
    !  Phys. Rev. A38 (1988) 3098-3100                             !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 23/02/2009.                       !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1in   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2in   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1in   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2in   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: den ! Density
    real(kind=DP) :: xsigma ! x_{sigma} in paper
    real(kind=DP) :: mgd, mgdsq ! Modulus of the the density gradient and square
    real(kind=DP) :: crden, densq ! Cube root of den and den squared
    real(kind=DP) :: crdenbq ! den ^ (4/3)
    real(kind=DP) :: asinhx ! Inverse hyperbolic sine of x
    real(kind=DP) :: dasinhx  ! Derivative wrt x of asinh(x)
    real(kind=DP) :: edenom, edenomsq ! Energy denominator and its square
    real(kind=DP) :: dnedenom ! Energy denominator density derivative
    real(kind=DP) :: dgedenom ! Energy denominator density gradient derivative
    ! qoh: Parameters:
    real(kind=DP), parameter :: beta = 0.0042_DP
    real(kind=DP), parameter :: sixbeta = 6.0_DP * beta
    real(kind=DP), parameter :: ldafac = -0.930525736349100025002010_DP
    real(kind=DP), parameter :: ftldafac = fthrd*ldafac

    ! qoh: Initialise values to zero to deal with small densities easily
    x_energy = 0.0_DP
    x_pot1 = 0.0_DP
    x_pot2 = 0.0_DP
    x_dfdmgd1 = 0.0_DP
    x_dfdmgd2 = 0.0_DP

    den = max(0.0_dp,den1in)
    spin1: if (den.gt.0.5_DP*dentol) then
       crden = den**THIRD
       crdenbq = crden*den
       densq = den**2
       mgd = max(0.0_dp,mgd1in)
       mgdsq = mgd**2
       xsigma = mgd / crdenbq
       asinhx = log(xsigma+sqrt(1.0_DP+xsigma**2))
       edenom = crdenbq+sixbeta*mgd*asinhx
       edenomsq = edenom**2

       x_energy = ldafac * crdenbq - beta * mgdsq /edenom

       dasinhx = 1.0_DP / sqrt(1.0_dp+xsigma**2)
       dnedenom = fthrd *(crden - sixbeta* mgdsq * dasinhx / (crden*densq))
       dgedenom = sixbeta* mgd * dasinhx / crdenbq + sixbeta*asinhx

       x_pot1 = ftldafac * crden + (beta * mgdsq * dnedenom)/ edenomsq
       x_dfdmgd1 = - 2.0_DP * beta * mgd / edenom &
            + (beta * mgdsq * dgedenom) / edenomsq
    end if spin1

    den = max(0.0_dp,den2in)
    spin2: if (den.gt.0.5_DP*dentol) then
       mgd = max(0.0_dp,mgd2in)
       crden = den**THIRD
       crdenbq = crden*den
       mgdsq = mgd**2
       xsigma = mgd / crdenbq
       asinhx = log(xsigma+sqrt(1.0_DP+xsigma**2))
       edenom = crdenbq+sixbeta*mgd*asinhx

       !qoh: Add to x_energy from spin 1
       x_energy = x_energy + ldafac * crdenbq - beta * mgdsq /edenom

       edenomsq = edenom**2
       densq = den**2
       dasinhx = 1.0_DP / sqrt(1.0_dp+xsigma**2)
       dnedenom = fthrd *(crden - sixbeta* mgdsq * dasinhx / (crden*densq))
       dgedenom = sixbeta* mgd * dasinhx / crdenbq + sixbeta*asinhx

       x_pot2 = ftldafac * crden + (beta * mgdsq * dnedenom)/ edenomsq
       x_dfdmgd2 = - 2.0_DP * beta * mgd / edenom &
            + (beta * mgdsq * dgedenom) / edenomsq
    end if spin2

  end subroutine xc_b88_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b1_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !===================================================================!
    ! This subroutine calculates the B1{LYP,PW91} exchange at a         !
    ! point in the spin unpolarised case.                               !
    ! B1{LYP,PW91}_X = 1/4 HF_X + 3/4 B88_X                             !
    !  C. Adamo and V. Barone, Chemical Physics Letters 274, 242 (1997) !
    !-------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                            !
    !===================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP), parameter :: b88_fac = 0.75_DP ! B1{LYP,PW91} a_x

    call xc_b88_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)
    x_energy = x_energy*b88_fac
    x_pot = x_pot * b88_fac
    x_dfdmgd = x_dfdmgd * b88_fac

  end subroutine xc_b1_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b3_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==================================================================!
    ! This subroutine calculates the B3{LYP,PW91} exchange at a point  !
    ! in the spin unpolarised case.                                    !
    ! B3{LYP,PW91} = (1 - a_0 - a_x ) LDA_X + a_0 HF_X + a_x B88_X     !
    !                + (1 - a_c) VWN_C + a_c {LYP,PW91}_C              !
    ! a_0 = 0.2, a_x = 0.72, a_c = 0.81
    ! A. D. Becke, J. Chem. Phys. 98, 5648 (1993)                      !
    ! P. J. Stephens, F. J. Devlin, C. F. Chabalowski, and             !
    !   M. J. Frisch, J. Phys. Chem. 98, 11623 (1994)                  !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                           !
    !==================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot

    real(kind=DP), parameter :: b88_fac = 0.72_DP ! B3{LYP,PW91} a_x
    real(kind=DP), parameter :: lda_fac = 0.08_DP ! B3{LYP,PW91} 1 - a_0 - a_x

    call xc_lda_x_point(den,lda_energy,lda_pot)
    call xc_b88_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)
    x_energy = x_energy*b88_fac + lda_energy*lda_fac
    x_pot = x_pot * b88_fac + lda_pot * lda_fac
    x_dfdmgd = x_dfdmgd * b88_fac

  end subroutine xc_b3_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_x_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the X exchange at a point         !
    ! in the spin unpolarised case.                                !
    ! X = 0.722 B88_X + 0.347 PW91_X - 0.069 LDA_X                 !
    ! X. Xu and W. Goddard, PNAS 101, 2673 (2004)                  !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 12/03/2009.                       !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot
    real(kind=DP) :: pw91_energy
    real(kind=DP) :: pw91_pot
    real(kind=DP) :: pw91_dfdmgd

    real(kind=DP), parameter :: b88_fac  =  0.722_DP
    real(kind=DP), parameter :: pw91_fac =  0.347_DP
    real(kind=DP), parameter :: lda_fac  = -0.069_DP

    call xc_lda_x_point(den,lda_energy,lda_pot)
    call xc_pw91_x_point(den,mgd,pw91_energy,pw91_pot,pw91_dfdmgd)
    call xc_b88_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    x_energy = x_energy*b88_fac + lda_energy*lda_fac + pw91_energy * pw91_fac
    x_pot = x_pot * b88_fac + lda_pot * lda_fac + pw91_pot * pw91_fac
    x_dfdmgd = x_dfdmgd * b88_fac + pw91_dfdmgd * pw91_fac

   end subroutine xc_x_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_x3_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the X3LYP exchange at a point     !
    ! in the spin unpolarised case.                                !
    ! X3LYP = a_0 HF_X + (1 - a_0 - a_x) LDA_X                     !
    !         + a_x ( 0.765 B88_X + 0.235 PW91_X) + a_c VWN_C      !
    !         + (1 - a_C) LYP_C                                    !
    ! a_0 = 0.218, a_x = 0.709, a_c = 0.129.                       !
    ! X. Xu and W. Goddard, PNAS 101, 2673 (2004)                  !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 12/03/2009.                       !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot    ! Local X pot contribution df_{x}/dn
    real(kind=dp), intent(out) :: x_dfdmgd ! Derivative df_{x}/d|grad n|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot
    real(kind=DP) :: pw91_energy
    real(kind=DP) :: pw91_pot
    real(kind=DP) :: pw91_dfdmgd

    real(kind=DP), parameter :: b88_fac = 0.765_DP * 0.709_DP
    real(kind=DP), parameter :: pw91_fac = 0.235_DP * 0.709_DP
    real(kind=DP), parameter :: lda_fac = 0.073_DP

    call xc_lda_x_point(den,lda_energy,lda_pot)
    call xc_pw91_x_point(den,mgd,pw91_energy,pw91_pot,pw91_dfdmgd)
    call xc_b88_x_point(den,mgd,x_energy,x_pot,x_dfdmgd)

    x_energy = x_energy*b88_fac + lda_energy*lda_fac + pw91_energy * pw91_fac
    x_pot = x_pot * b88_fac + lda_pot * lda_fac + pw91_pot * pw91_fac
    x_dfdmgd = x_dfdmgd * b88_fac + pw91_dfdmgd * pw91_fac

  end subroutine xc_x3_x_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b1_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !===================================================================!
    ! This subroutine calculates the B1{LYP,PW91} exchange at a         !
    ! point in the spin polarised case.                                 !
    ! B1{LYP,PW91}_X = 1/4 HF_X + 3/4 B88_X                             !
    !  C. Adamo and V. Barone, Chemical Physics Letters 274, 242 (1997) !
    !-------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                            !
    !===================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP), parameter :: b88_fac = 0.75_DP ! B1{LYP,PW91} a_x

    call xc_b88_x_point_sp(den1,den2,mgd1,mgd2,&
         x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)
    x_energy = x_energy*b88_fac
    x_pot1 = x_pot1 * b88_fac
    x_pot2 = x_pot2 * b88_fac
    x_dfdmgd1 = x_dfdmgd1 * b88_fac
    x_dfdmgd2 = x_dfdmgd2 * b88_fac

  end subroutine xc_b1_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b3_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==================================================================!
    ! This subroutine calculates the B3{LYP,PW91} exchange at a point  !
    ! in the spin polarised case.                                      !
    ! B3{LYP,PW91} = (1 - a_0 - a_x ) LDA_X + a_0 HF_X + a_x B88_X     !
    !                + (1 - a_c) VWN_C + a_c {LYP,PW91}_C              !
    ! a_0 = 0.2, a_x = 0.72, a_c = 0.81                                !
    ! A. D. Becke, J. Chem. Phys. 98, 5648 (1993)                      !
    ! P. J. Stephens, F. J. Devlin, C. F. Chabalowski, and             !
    !   M. J. Frisch, J. Phys. Chem. 98, 11623 (1994)                  !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                           !
    !==================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot1
    real(kind=DP) :: lda_pot2

    real(kind=DP), parameter :: b88_fac = 0.72_DP ! B3{LYP,PW91} a_x
    real(kind=DP), parameter :: lda_fac = 0.08_DP ! B3{LYP,PW91} 1 - a_0 - a_x

    call xc_lsda_x_point(den1,den2,lda_energy,lda_pot1, lda_pot2)
    call xc_b88_x_point_sp(den1,den2,mgd1,mgd2,&
         x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)
    x_energy = x_energy*b88_fac + lda_energy*lda_fac
    x_pot1 = x_pot1 * b88_fac + lda_pot1 * lda_fac
    x_pot2 = x_pot2 * b88_fac + lda_pot2 * lda_fac
    x_dfdmgd1 = x_dfdmgd1 * b88_fac
    x_dfdmgd2 = x_dfdmgd2 * b88_fac

  end subroutine xc_b3_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_x_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the X exchange at a point         !
    ! in the spin polarised case.                                  !
    ! X = 0.722 B88_X + 0.347 PW91_X - 0.069 * LDA_X               !
    ! X. Xu and W. Goddard, PNAS 101, 2673 (2004)                  !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 12/03/2009.                       !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot1
    real(kind=DP) :: lda_pot2
    real(kind=DP) :: pw91_energy
    real(kind=DP) :: pw91_pot1
    real(kind=DP) :: pw91_pot2
    real(kind=DP) :: pw91_dfdmgd1
    real(kind=DP) :: pw91_dfdmgd2

    real(kind=DP), parameter :: b88_fac  =  0.722_DP
    real(kind=DP), parameter :: pw91_fac =  0.347_DP
    real(kind=DP), parameter :: lda_fac  = -0.069_DP

    call xc_lsda_x_point(den1,den2,lda_energy,lda_pot1, lda_pot2)
    call xc_pw91_x_point_sp(den1,den2,mgd1,mgd2,&
         pw91_energy,pw91_pot1,pw91_pot2,pw91_dfdmgd1,pw91_dfdmgd2)
    call xc_b88_x_point_sp(den1,den2,mgd1,mgd2,&
         x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    x_energy = x_energy*b88_fac + lda_energy*lda_fac + pw91_energy * pw91_fac
    x_pot1 = x_pot1 * b88_fac + lda_pot1 * lda_fac + pw91_pot1 * pw91_fac
    x_pot2 = x_pot2 * b88_fac + lda_pot2 * lda_fac + pw91_pot2 * pw91_fac
    x_dfdmgd1 = x_dfdmgd1 * b88_fac + pw91_dfdmgd1 * pw91_fac
    x_dfdmgd2 = x_dfdmgd2 * b88_fac + pw91_dfdmgd2 * pw91_fac

  end subroutine xc_x_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_x3_x_point_sp(den1,den2,mgd1,mgd2,&
       x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    !==============================================================!
    ! This subroutine calculates the X3LYP exchange at a point     !
    ! in the spin polarised case.                                  !
    ! X3LYP = a_0 HF_X + (1 - a_0 - a_x) LDA_X                     !
    !         + a_x ( 0.765 B88_X + 0.235 PW91_X) + a_c VWN_C      !
    !         + (1 - a_C) LYP_C                                    !
    ! a_0 = 0.218, a_x = 0.709, a_c = 0.129.                       !
    ! X. Xu and W. Goddard, PNAS 101, 2673 (2004)                  !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 12/03/2009.                       !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: x_energy ! The exchange energy
    real(kind=dp), intent(out) :: x_pot1  ! Local X pot contribution df_{x}/dn s1
    real(kind=dp), intent(out) :: x_pot2  ! Local X pot contribution df_{x}/dn s2
    real(kind=dp), intent(out) :: x_dfdmgd1 ! Derivative df_{x}/d|grad n| spin 1
    real(kind=dp), intent(out) :: x_dfdmgd2 ! Derivative df_{x}/d|grad n| spin 2

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot1
    real(kind=DP) :: lda_pot2
    real(kind=DP) :: pw91_energy
    real(kind=DP) :: pw91_pot1
    real(kind=DP) :: pw91_pot2
    real(kind=DP) :: pw91_dfdmgd1
    real(kind=DP) :: pw91_dfdmgd2

    real(kind=DP), parameter :: b88_fac = 0.765_DP * 0.709_DP
    real(kind=DP), parameter :: pw91_fac = 0.235_DP * 0.709_DP
    real(kind=DP), parameter :: lda_fac = 0.073_DP

    call xc_lsda_x_point(den1,den2,lda_energy,lda_pot1, lda_pot2)
    call xc_pw91_x_point_sp(den1,den2,mgd1,mgd2,&
         pw91_energy,pw91_pot1,pw91_pot2,pw91_dfdmgd1,pw91_dfdmgd2)
    call xc_b88_x_point_sp(den1,den2,mgd1,mgd2,&
         x_energy,x_pot1,x_pot2,x_dfdmgd1,x_dfdmgd2)

    x_energy = x_energy*b88_fac + lda_energy*lda_fac + pw91_energy * pw91_fac
    x_pot1 = x_pot1 * b88_fac + lda_pot1 * lda_fac + pw91_pot1 * pw91_fac
    x_pot2 = x_pot2 * b88_fac + lda_pot2 * lda_fac + pw91_pot2 * pw91_fac
    x_dfdmgd1 = x_dfdmgd1 * b88_fac + pw91_dfdmgd1 * pw91_fac
    x_dfdmgd2 = x_dfdmgd2 * b88_fac + pw91_dfdmgd2 * pw91_fac

  end subroutine xc_x3_x_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lyp_c_point(denin,mgdin,c_energy,c_pot,c_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the LYP correlation at a point    !
    ! in the spin unpolarised case.                                !
    !     C. Lee, W. Yang, and R.G. Parr                           !
    !     Development of the Colle-Salvetti correlation-energy     !
    !     formula into a functional of the electron density        !
    !     Phys. Rev. B37 (1988) 785-789                            !
    !  Using reformulation (to avoid Laplacian) given by:          !
    !     B. Miehlich, A. Savin, H. Stoll and H. Preuss            !
    !     Results obtained with the correlation energy density     !
    !     functionals of becke and Lee, Yang and Parr              !
    !     Chem. Phys. Lett. 157 (1989) 200-206                     !
    !  Note: There is a mistake in Handy's formula in ESQC book.   !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 24/02/2009                        !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: denin    ! The charge density at a point
    real(kind=dp), intent(in)  :: mgdin    ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy      ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    real(kind=DP) :: den, rden ! Density and its reciprocal
    real(kind=DP) :: mgd, mgdsq ! Modulus of density gradient and it's square
    real(kind=DP) :: crden ! Cube root of den
    real(kind=DP) :: rcrden ! Reciprocal cube root of den
    real(kind=DP) :: densq, rdensq ! den squared and its reciprocal
    !real(kind=DP) :: rcrdensq ! Reciprocal cube root of den squared
    real(kind=DP) :: redenom, redenomsq ! Reciprocal of E denominator and square
    real(kind=DP) :: expo ! Exponential term
    real(kind=DP) :: omega ! As in paper
    real(kind=DP) :: delta ! As in paper
    real(kind=DP) :: bracket
    real(kind=DP) :: dndelta ! Density derivative of density
    real(kind=DP) :: dnomega ! Density derivative of omega

    ! qoh: Parameters:
    real(kind=DP), parameter :: aa = 0.04918_DP ! As in paper
    real(kind=DP), parameter :: bb = 0.132_DP   ! As in paper
    real(kind=DP), parameter :: ab = aa * bb
    real(kind=DP), parameter :: cc = 0.25330_DP ! As in paper
    real(kind=DP), parameter :: dd = 0.3490_DP  ! As in paper
    real(kind=DP), parameter :: ETHRD = 8.0_DP / 3.0_DP
    real(kind=DP), parameter :: cf = 2.8712340001881911_DP ! 0.3*(3*PI^2)^(2/3)
    real(kind=DP), parameter :: se72 = 7.0_DP/72.0_DP
    real(kind=DP), parameter :: o24 = 1.0_DP/24.0_DP

    den = max(0.0_dp,denin)
    tol: if(den.gt.dentol) then
       ! qoh: Calculate useful powers of den
       densq = den**2
       rdensq = 1.0_DP / densq
       rden = 1.0_DP / den
       crden = den**THIRD
       rcrden = 1/crden
       !rcrdensq = rcrden**2
       ! qoh: Initialise density gradient and square
       mgd = max(0.0_dp,mgdin)
       mgdsq = mgd**2

       redenom = 1/( 1.0_dp+dd*rcrden )
       redenomsq = redenom**2
       expo = exp(-cc*rcrden)

       omega = expo * crden * redenom * rdensq * rdensq
       delta = cc*rcrden + dd*rcrden * redenom

       bracket = CF*densq*den*rcrden - mgdsq*(o24+se72*delta)
       dnomega = ( THIRD * rdensq * rdensq * redenom * expo * rden ) * &
            (cc - 11.0_DP*crden + dd * redenom )
       dndelta = - THIRD * rden * rcrden * ( cc + dd * redenomsq )

       c_energy = - aa*den*redenom - ab * omega*densq * bracket

       c_pot = - aa * ( redenom + THIRD*rcrden*dd*redenomsq ) &
            - ab * ( dnomega*densq + 2.0_DP*omega*den) * bracket &
            - ab * omega*densq  * (ETHRD*CF*densq*rcrden &
            - mgdsq * (se72*dndelta))

       !qoh: Note: No minus sign since we have merged brackets
       c_dfdmgd = ab * omega * densq * mgd * 2.0_DP * ( o24 + se72 * delta)
    else
       ! qoh: Set all values to zero for very small densities
       c_energy = 0.00_DP
       c_pot = 0.0_DP
       c_dfdmgd = 0.0_DP
    endif tol

  end subroutine xc_lyp_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lyp_c_point_sp(denin,den1in,den2in,mgdin,mgd1in,mgd2in,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==============================================================!
    ! This subroutine calculates the LYP correlation at a point    !
    ! in the spin polarised case.                                  !
    !     C. Lee, W. Yang, and R.G. Parr                           !
    !     Development of the Colle-Salvetti correlation-energy     !
    !     formula into a functional of the electron density        !
    !     Phys. Rev. B37 (1988) 785-789                            !
    !  Using reformulation (to avoid Laplacian) given by:          !
    !     B. Miehlich, A. Savin, H. Stoll and H. Preuss            !
    !     Results obtained with the correlation energy density     !
    !     functionals of becke and Lee, Yang and Parr              !
    !     Chem. Phys. Lett. 157 (1989) 200-206                     !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 04/03/2009                        !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: denin    ! Density at point
    real(kind=dp), intent(in)  :: den1in   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2in   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgdin    ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd1in   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2in   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|

    ! qoh: Local variables
    real(kind=DP) :: den, rden   ! Total density and its reciprocal
    real(kind=DP) :: den1!, rden1 ! Spin 1 density and its reciprocal
    real(kind=DP) :: den2!, rden2 ! Spin 2 density and its reciprocal
    real(kind=DP) :: mgd, mgdsq  ! Mod density gradient and square
    real(kind=DP) :: mgd1, mgd1sq ! Mod spin 1 density gradient and square
    real(kind=DP) :: mgd2, mgd2sq ! Mod spin 2 density gradient and square
    real(kind=DP) :: crden ! Cube root of den
    real(kind=DP) :: rcrden ! Reciprocal cube root of den
    real(kind=DP) :: densq, rdensq ! den squared and its reciprocal
    real(kind=DP) :: crden1sq ! Cube root of den spin 1 squared
    real(kind=DP) :: den1sq ! den spin 1 squared
    real(kind=DP) :: crden2sq ! Cube root of den spin 2 squared
    real(kind=DP) :: den2sq ! den spin 2 squared
    real(kind=DP) :: redenom, redenomsq ! Reciprocal of E denominator and square
    real(kind=DP) :: expo ! Exponential term
    !real(kind=DP) :: rcrdensq ! Reciprocal of the cube root of den squared
    real(kind=DP) :: omega ! As in paper
    real(kind=DP) :: delta ! As in paper
    real(kind=DP) :: dndelta ! delta density derivative
    real(kind=DP) :: dnomega ! omega density derivative
    real(kind=DP) :: inbracket ! Inner bracket [] from paper
    real(kind=DP) :: outbracket ! Outer bracket {} from paper
    real(kind=DP) :: dn1inbracket ! Spin 1 density derivative of inner bracket
    real(kind=DP) :: dn2inbracket ! Spin 2 density derivative of inner bracket

    ! qoh: Parameters:
    real(kind=DP), parameter :: aa = 0.04918_DP ! As in paper
    real(kind=DP), parameter :: bb = 0.132_DP   ! As in paper
    real(kind=DP), parameter :: ab = aa * bb
    real(kind=DP), parameter :: cc = 0.25330_DP ! As in paper
    real(kind=DP), parameter :: dd = 0.3490_DP  ! As in paper
    real(kind=DP), parameter :: TTHRD = 2.0_DP / 3.0_DP
    real(kind=DP), parameter :: ETHRD = 8.0_DP / 3.0_DP
    real(kind=DP), parameter :: cfp = 36.462398978764767_DP ! C_F * 2^{11/3}
    real(kind=DP), parameter :: on9 = 1.0_DP/9.0_DP
    real(kind=DP), parameter :: tw9 = 2.0_DP/9.0_DP
    real(kind=DP), parameter :: el9 = 11.0_DP/9.0_DP
    real(kind=DP), parameter :: tt9 = 22.0_DP/9.0_DP
    real(kind=DP), parameter :: on18 = 1.0_DP/18.0_DP
    real(kind=DP), parameter :: se18 = 7.0_DP/18.0_DP
    real(kind=DP), parameter :: fs18 = 47.0_DP/18.0_DP

    den = max(0.0_dp,denin)
    den1 = max(0.0_dp,den1in)
    den2 = max(0.0_dp,den2in)
    tol: if(den1.gt.dentol .and. den2.gt.dentol) then

       ! qoh: Calculate useful powers of den[ 12]
       densq = den**2
       rdensq = 1.0_DP / densq
       rden = 1.0_DP / den
       crden = den**THIRD
       rcrden = 1/crden
       !rcrdensq = rcrden**2
       den1sq = den1**2
       crden1sq = den1**TTHRD
       den2sq = den2**2
       crden2sq = den2**TTHRD
       ! qoh: Modulus of density gradient and square initialised
       mgd = max(0.0_dp,mgdin)
       mgdsq = mgd**2
       mgd1 = max(0.0_dp,mgd1in)
       mgd1sq = mgd1**2
       mgd2 = max(0.0_dp,mgd2in)
       mgd2sq = mgd2**2

       redenom = 1/( 1.0_dp+dd*rcrden )
       redenomsq = redenom**2
       expo = exp(-cc*rcrden)

       omega = expo * crden * redenom * rdensq * rdensq
       delta = cc*rcrden + dd*rcrden * redenom

       inbracket = CFP*(den1sq*crden1sq + den2sq*crden2sq)  &
            + mgdsq*(fs18-se18*delta) - (2.5_DP - on18*delta)*(mgd1sq+mgd2sq) &
            - (on9*delta -el9)*(den1*rden*mgd1sq+den2*rden*mgd2sq)

       outbracket = den1*den2 * inbracket - TTHRD*densq*mgdsq&
            + (TTHRD*densq-den1sq)*mgd2sq +  (TTHRD*densq-den2sq)*mgd1sq

       dnomega = ( THIRD * rdensq * rdensq * redenom * expo * rden ) * &
            (cc - 11.0_DP*crden + dd * redenom )
       dndelta = - THIRD * rden * rcrden * ( cc + dd * redenomsq )

       dn1inbracket = CFP*ETHRD*den1*crden1sq &
            - dndelta*(se18*mgdsq - on18*(mgd1sq + mgd2sq)) &
            - dndelta*on9*rden*(den1*mgd1sq + den2*mgd2sq) &
            - (on9*delta-el9)*den2*rdensq*(mgd1sq - mgd2sq)

       dn2inbracket =  CFP*ETHRD*den2*crden2sq &
            - dndelta*(se18*mgdsq - on18*(mgd1sq + mgd2sq)) &
            - dndelta*on9*rden*(den1*mgd1sq + den2*mgd2sq) &
            - (on9*delta-el9)*den1*rdensq*(mgd2sq - mgd1sq)

       c_energy = - aa*4.0_DP*den1*den2*rden*redenom - ab * omega * outbracket

       ! qoh: Calculate potential
       c_pot1 = -aa * 4.0_DP*(redenom *(den2sq*rdensq)&
            + THIRD*redenomsq*dd*rcrden*den1*den2 * rdensq)&
            - ab*dnomega*outbracket -ab*omega*den2*inbracket &
            - ab * omega*den1*den2*dn1inbracket &
            + ab*omega*FTHRD*den*(mgdsq-mgd1sq-mgd2sq)&
            + ab*omega*2.0_DP*den1*mgd2sq

       c_pot2 = -aa * 4.0_DP*(redenom *(den1sq*rdensq)&
            + THIRD*redenomsq*dd*rcrden*den1*den2 * rdensq)&
            - ab*dnomega*outbracket -ab*omega*den1*inbracket &
            - ab*omega*den1*den2* dn2inbracket &
            + ab*omega*FTHRD*den*(mgdsq-mgd1sq-mgd2sq)&
            + ab*omega*2.0_DP*den2*mgd1sq

       c_dfdmgd = 0.0_DP

       c_dfdmgd12 = - ab * omega * ( den1 * den2 *(fs18 - se18 * delta) &
            - TTHRD*densq) * (mgdsq - mgd1sq - mgd2sq)

       c_dfdmgd1 = -ab * omega * mgd1 * (den1 * den2 *((tw9 - TTHRD*delta) &
            - (tw9*delta - tt9) * den1 * rden) - 2.0_DP * den2sq )

       c_dfdmgd2 = -ab * omega * mgd2 * ( den1 * den2 *((tw9 - TTHRD*delta) &
            - (tw9*delta - tt9) * den2 * rden) - 2.0_DP * den1sq )

     else
       ! qoh: Set all values to zero for very small densities
       c_energy = 0.0_DP
       c_pot1 = 0.0_DP
       c_pot2 = 0.0_DP
       c_dfdmgd = 0.0_DP
       c_dfdmgd1 = 0.0_DP
       c_dfdmgd2 = 0.0_DP
       c_dfdmgd12 = 0.0_DP
    endif tol

  end subroutine xc_lyp_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b3lyp_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==================================================================!
    ! This subroutine calculates the B3LYP correlation at a point      !
    ! in the spin unpolarised case.                                    !
    ! B3LYP = (1 - a_0 - a_x ) LDA_X + a_0 HF_X + a_x B88_X            !
    !                + (1 - a_c) VWN_C + a_c LYP_C                     !
    ! A. D. Becke, J. Chem. Phys. 98, 5648 (1993)                      !
    ! P. J. Stephens, F. J. Devlin, C. F. Chabalowski, and             !
    !   M. J. Frisch, J. Phys. Chem. 98, 11623 (1994)                  !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                           !
    !==================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy      ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot

    real(kind=DP), parameter :: lda_fac = 0.19_DP ! 1- a_c
    real(kind=DP), parameter :: lyp_fac = 0.81_DP ! a_c

    call xc_vwn_c_point(den,lda_energy,lda_pot)
    call xc_lyp_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)
    c_energy = c_energy*lyp_fac + lda_fac * lda_energy
    c_pot = c_pot * lyp_fac + lda_fac * lda_pot
    c_dfdmgd = c_dfdmgd * lyp_fac

   end subroutine xc_b3lyp_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b3lyp_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==================================================================!
    ! This subroutine calculates the B3LYP correlation at a point      !
    ! in the spin polarised case.                                      !
    ! B3LYP = (1 - a_0 - a_x ) LDA_X + a_0 HF_X + a_x B88_X            !
    !                + (1 - a_c) VWN_C + a_c LYP_C                     !
    ! A. D. Becke, J. Chem. Phys. 98, 5648 (1993)                      !
    ! P. J. Stephens, F. J. Devlin, C. F. Chabalowski, and             !
    !   M. J. Frisch, J. Phys. Chem. 98, 11623 (1994)                  !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                           !
    !==================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den   ! Density at point
    real(kind=dp), intent(in)  :: den1   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd    ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd1   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot1
    real(kind=DP) :: lda_pot2

    real(kind=DP), parameter :: lda_fac = 0.19_DP ! 1- a_c
    real(kind=DP), parameter :: lyp_fac = 0.81_DP ! a_c

    call xc_vwn_c_point_sp(den,den1,den2,lda_energy,lda_pot1, lda_pot2)
    call xc_lyp_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
         c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)
    c_energy = c_energy*lyp_fac + lda_fac * lda_energy
    c_pot1 = c_pot1 * lyp_fac + lda_fac * lda_pot1
    c_pot2 = c_pot2 * lyp_fac + lda_fac * lda_pot2
    c_dfdmgd = c_dfdmgd * lyp_fac
    c_dfdmgd1 = c_dfdmgd1 * lyp_fac
    c_dfdmgd2 = c_dfdmgd2 * lyp_fac
    c_dfdmgd12 = c_dfdmgd12 * lyp_fac

  end subroutine xc_b3lyp_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b3pw91_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==================================================================!
    ! This subroutine calculates the B3PW91 correlation at a point     !
    ! in the spin unpolarised case.                                    !
    ! B3PW91 = (1 - a_0 - a_x ) LDA_X + a_0 HF_X + a_x B88_X           !
    !                + (1 - a_c) VWN_C + a_c PW91_C                    !
    ! A. D. Becke, J. Chem. Phys. 98, 5648 (1993)                      !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                           !
    !==================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy      ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot

    real(kind=DP), parameter :: lda_fac  = 0.19_DP ! 1- a_c
    real(kind=DP), parameter :: pw91_fac = 0.81_DP ! a_c

    call xc_vwn_c_point(den,lda_energy,lda_pot)
    call xc_pw91_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)
    c_energy = c_energy*pw91_fac + lda_fac * lda_energy
    c_pot = c_pot * pw91_fac + lda_fac * lda_pot
    c_dfdmgd = c_dfdmgd * pw91_fac

   end subroutine xc_b3pw91_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_b3pw91_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==================================================================!
    ! This subroutine calculates the B3PW91 correlation at a point     !
    ! in the spin polarised case.                                      !
    ! B3PW91 = (1 - a_0 - a_x ) LDA_X + a_0 HF_X + a_x B88_X           !
    !                + (1 - a_c) VWN_C + a_c PW91_C                    !
    ! A. D. Becke, J. Chem. Phys. 98, 5648 (1993)                      !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 11/03/2009.                           !
    !==================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den   ! Density at point
    real(kind=dp), intent(in)  :: den1   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd    ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd1   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad_n_2|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot1
    real(kind=DP) :: lda_pot2

    real(kind=DP), parameter :: lda_fac  = 0.19_DP ! 1- a_c
    real(kind=DP), parameter :: pw91_fac = 0.81_DP ! a_c

    call xc_vwn_c_point_sp(den,den1,den2,lda_energy,lda_pot1, lda_pot2)
    call xc_pw91_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
         c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)
    c_energy = c_energy*pw91_fac + lda_fac * lda_energy
    c_pot1 = c_pot1 * pw91_fac + lda_fac * lda_pot1
    c_pot2 = c_pot2 * pw91_fac + lda_fac * lda_pot2
    c_dfdmgd = c_dfdmgd * pw91_fac
    c_dfdmgd1 = c_dfdmgd1 * pw91_fac
    c_dfdmgd2 = c_dfdmgd2 * pw91_fac
    c_dfdmgd12 = c_dfdmgd12 * pw91_fac

  end subroutine xc_b3pw91_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_x3lyp_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==============================================================!
    ! This subroutine calculates the X3LYP correlation at a point  !
    ! in the spin unpolarised case.                                !
    ! X3LYP = a_0 HF_X + (1 - a_0 - a_x) LDA_X                     !
    !         + a_x ( 0.765 B88_X + 0.235 PW91_X) + a_c VWN_C      !
    !         + (1 - a_C) LYP_C                                    !
    ! a_0 = 0.218, a_x = 0.709, a_c = 0.129.                       !
    ! X. Xu and W. Goddard, PNAS 101, 2673 (2004)                  !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 12/03/2009.                       !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy      ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot

    real(kind=DP), parameter :: lda_fac = 0.129_DP ! a_c
    real(kind=DP), parameter :: lyp_fac = 0.871_DP ! 1 - a_c

    call xc_vwn_c_point(den,lda_energy,lda_pot)
    call xc_lyp_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)
    c_energy = c_energy*lyp_fac + lda_fac * lda_energy
    c_pot = c_pot * lyp_fac + lda_fac * lda_pot
    c_dfdmgd = c_dfdmgd * lyp_fac

  end subroutine xc_x3lyp_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_x3lyp_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==============================================================!
    ! This subroutine calculates the X3LYP correlation at a point  !
    ! in the spin polarised case.                                  !
    ! X3LYP = a_0 HF_X + (1 - a_0 - a_x) LDA_X                     !
    !         + a_x ( 0.765 B88_X + 0.235 PW91_X) + a_c VWN_C      !
    !         + (1 - a_C) LYP_C                                    !
    ! a_0 = 0.218, a_x = 0.709, a_c = 0.129.                       !
    ! X. Xu and W. Goddard, PNAS 101, 2673 (2004)                  !
    !--------------------------------------------------------------!
    ! Written by Quintin Hill on 12/03/2009.                       !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den   ! Density at point
    real(kind=dp), intent(in)  :: den1   ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2   ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd    ! mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd1   ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2   ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1  ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2  ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|

    real(kind=DP) :: lda_energy
    real(kind=DP) :: lda_pot1
    real(kind=DP) :: lda_pot2

    real(kind=DP), parameter :: lda_fac = 0.129_DP ! a_c
    real(kind=DP), parameter :: lyp_fac = 0.871_DP ! 1 - a_c

    call xc_vwn_c_point_sp(den,den1,den2,lda_energy,lda_pot1, lda_pot2)
    call xc_lyp_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
         c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)
    c_energy = c_energy*lyp_fac + lda_fac * lda_energy
    c_pot1 = c_pot1 * lyp_fac + lda_fac * lda_pot1
    c_pot2 = c_pot2 * lyp_fac + lda_fac * lda_pot2
    c_dfdmgd = c_dfdmgd * lyp_fac
    c_dfdmgd1 = c_dfdmgd1 * lyp_fac
    c_dfdmgd2 = c_dfdmgd2 * lyp_fac
    c_dfdmgd12 = c_dfdmgd12 * lyp_fac

  end subroutine xc_x3lyp_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_none_c_point(den,mgd,c_energy,c_pot,c_dfdmgd)

    !==============================================================!
    ! This subroutine is a dummy for adding no correlation.        !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! The charge density at a point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(d) at the same point
    real(kind=dp), intent(out) :: c_energy ! The C energy
    real(kind=dp), intent(out) :: c_pot    ! Local C pot contribution dfxc/dn
    real(kind=dp), intent(out) :: c_dfdmgd ! Derivative  dfxc/d|grad n|

    if (den==mgd) c_energy = 0.0_DP
    c_energy = 0.0_DP
    c_pot = 0.0_DP
    c_dfdmgd = 0.0_DP

  end subroutine xc_none_c_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_none_c_point_sp(den,den1,den2,mgd,mgd1,mgd2,&
       c_energy,c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)

    !==============================================================!
    ! This subroutine is a dummy for adding no correlation.        !
    !==============================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=dp), intent(in)  :: den      ! Density at point
    real(kind=dp), intent(in)  :: den1     ! Density at point for spin 1
    real(kind=dp), intent(in)  :: den2     ! Density at point for spin 2
    real(kind=dp), intent(in)  :: mgd1     ! Spin1 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd2     ! Spin2 mod_grad(n) at the same point
    real(kind=dp), intent(in)  :: mgd      ! mod_grad(n) at the same point
    real(kind=dp), intent(out) :: c_energy ! The correlation energy
    real(kind=dp), intent(out) :: c_pot1   ! Local C pot contribution df_{c}/dn s1
    real(kind=dp), intent(out) :: c_pot2   ! Local C pot contribution df_{c}/dn s2
    real(kind=dp), intent(out) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=dp), intent(out) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n_1|
    real(kind=dp), intent(out) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n_2|
    real(kind=dp), intent(out) :: c_dfdmgd12  ! Derivative df_{c}/d|grad n_1.grad n_2|

    if ((den==den1).and.(den==den2).and.(mgd==mgd1).and.(mgd==mgd2)) &
         c_energy = 0.0_DP
    c_energy = 0.0_DP
    c_pot1 = 0.0_DP
    c_pot2 = 0.0_DP
    c_dfdmgd = 0.0_DP
    c_dfdmgd1 = 0.0_DP
    c_dfdmgd2 = 0.0_DP
    c_dfdmgd12 = 0.0_DP

  end subroutine xc_none_c_point_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module
