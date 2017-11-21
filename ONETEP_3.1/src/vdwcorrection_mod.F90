! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!          Van der Waals empirical correction module             !
!                                                                !
! This module contains subroutines for calculating a Van der     !
! Waal energy correction as a sum of damped London potentials.   !
!----------------------------------------------------------------!
! Written by Quintin Hill in 2007/8 with assistance from         !
! Chris-Kriton Skylaris.                                         !
! Forces added July 2008 by Quintin Hill                         !
!================================================================!

module vdwcorrection

  use constants, only: DP

  implicit none

  private

  ! qoh: Structure to store parameters in

  TYPE VDWPARAMETERS
     !qoh: C_6 coefficients
     REAL(kind=DP), DIMENSION(109) :: C6COEFF
     !qoh: R_0 values
     REAL(kind=DP), DIMENSION(109) :: RADZERO
     !qoh: damping coefficient
     REAL(kind=DP), DIMENSION(3) :: DCOEFF
     !qoh: Effective number of electrons
     REAL(kind=DP), DIMENSION(109) :: NEFF
  END TYPE VDWPARAMETERS

  public :: vdwcorrection_calculate_energy
  public :: vdwcorrection_calculate_forces
  public :: vdwcorrection_override_alloc
  public :: vdwcorrection_override_dealloc

  ! qoh: Share array of Van der Waals parameters accross module
  type(VDWPARAMETERS) :: vdwparams ! Array of parameters
  real(kind=DP), save, public :: pub_dispersion_energy
  ! The dispersion correction energy
  character(len=80), save, allocatable :: vdwparam_override(:)

  real(kind=DP), parameter :: cutoff   = 20.0_DP ! Cutoff distance in Bohr
  real(kind=DP), parameter :: cutoffsq = 400.0_DP ! Square of cutoff distance

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdwcorrection_initializeparams

    !==================================================================!
    ! This subroutine populates the array of Van der Waals parameters  !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in 2007, with some modifications in 2008 !
    ! Some unoptimised parameters (e.g. for P) were taken from Halgren.!
    ! Journal of the American Chemical Society 114(20), 7827-7843, 1992!
    !==================================================================!

    use rundat, only: xc_functional, pub_dispersion, pub_libxc_c_func_id, &
         pub_libxc_x_func_id
#ifdef LIBXC
    use libxc_funcs_m
#endif

    implicit none

    ! ndmh: Local Variables
    character(len=80) :: xc_func

    ! ndmh: set up xc_func, overriding parameter for libxc functionals
    xc_func = xc_functional
#ifdef LIBXC
    if (xc_func=='LIBXC') then
       if ((pub_libxc_x_func_id==XC_GGA_X_B88) .and. &
            (pub_libxc_c_func_id==XC_GGA_C_LYP)) then
          xc_func = 'BLYP'
       else if ((pub_libxc_x_func_id==XC_GGA_X_PBE) .and. &
            (pub_libxc_c_func_id==XC_GGA_C_PBE)) then
          xc_func = 'PBE'
       else if ((pub_libxc_x_func_id==XC_GGA_X_PW91) .and. &
            (pub_libxc_c_func_id==XC_GGA_C_PW91)) then
          xc_func = 'PW91'
       else if ((pub_libxc_x_func_id==XC_GGA_X_PBE_R) .and. &
            (pub_libxc_c_func_id==XC_GGA_C_PBE)) then
          xc_func = 'REVPBE'
       else if ((pub_libxc_x_func_id==XC_GGA_X_RPBE) .and. &
            (pub_libxc_c_func_id==XC_GGA_C_PBE)) then
          xc_func = 'RPBE'
       else if ((pub_libxc_x_func_id==XC_GGA_XC_XLYP)) then
          xc_func = 'XLYP'
       end if
    end if
#endif

    ! qoh: Initialise the vdwparams array
    vdwparams%c6coeff = 0.0_DP
    vdwparams%neff    = 0.0_DP
    ! qoh: Setting to 1 prevents division by zero
    vdwparams%radzero = 1.0_DP
    ! qoh: Populate the vdwparams array

    select case (xc_func)

    case('PBE')

       pbedamp: select case (pub_dispersion)

       case(1) pbedamp

          vdwparams%dcoeff(1)=3.2607_DP
          vdwparams%dcoeff(2)=3.5400_DP
          vdwparams%dcoeff(3)=23.0000_DP
          vdwparams%c6coeff(1)=2.9239_DP
          vdwparams%c6coeff(6)=27.3561_DP
          vdwparams%c6coeff(7)=19.5089_DP
          vdwparams%c6coeff(8)=11.7697_DP
          vdwparams%c6coeff(15)=190.5033_DP
          vdwparams%c6coeff(16)=161.0368_DP
          vdwparams%radzero(1)=2.9635_DP
          vdwparams%radzero(2)=3.0000_DP
          vdwparams%radzero(3)=3.8000_DP
          vdwparams%radzero(4)=3.8000_DP
          vdwparams%radzero(5)=3.8000_DP
          vdwparams%radzero(6)=3.3610_DP
          vdwparams%radzero(7)=3.5136_DP
          vdwparams%radzero(8)=3.7294_DP
          vdwparams%radzero(9)=3.8000_DP
          vdwparams%radzero(10)=3.8000_DP
          vdwparams%radzero(11)=4.8000_DP
          vdwparams%radzero(12)=4.8000_DP
          vdwparams%radzero(13)=4.8000_DP
          vdwparams%radzero(14)=4.8000_DP
          vdwparams%radzero(15)=4.8000_DP
          vdwparams%radzero(16)=5.1033_DP

       case(2) pbedamp

          vdwparams%dcoeff(1)=3.0000_DP
          vdwparams%dcoeff(2)=3.5400_DP
          vdwparams%dcoeff(3)=23.0000_DP
          vdwparams%c6coeff(1)=2.8450_DP
          vdwparams%c6coeff(6)=27.3200_DP
          vdwparams%c6coeff(7)=19.4800_DP
          vdwparams%c6coeff(8)=11.7600_DP
          vdwparams%c6coeff(15)=190.5033_DP
          vdwparams%c6coeff(16)=161.0232_DP
          vdwparams%radzero(1)=3.0000_DP
          vdwparams%radzero(2)=3.0000_DP
          vdwparams%radzero(3)=3.8000_DP
          vdwparams%radzero(4)=3.8000_DP
          vdwparams%radzero(5)=3.8000_DP
          vdwparams%radzero(6)=3.8000_DP
          vdwparams%radzero(7)=3.8000_DP
          vdwparams%radzero(8)=3.8000_DP
          vdwparams%radzero(9)=3.8000_DP
          vdwparams%radzero(10)=3.8000_DP
          vdwparams%radzero(11)=4.8000_DP
          vdwparams%radzero(12)=4.8000_DP
          vdwparams%radzero(13)=4.8000_DP
          vdwparams%radzero(14)=4.8000_DP
          vdwparams%radzero(15)=4.8000_DP
          vdwparams%radzero(16)=5.5998_DP

       case(3) pbedamp

          vdwparams%dcoeff(1)=3.0000_DP
          vdwparams%dcoeff(2)=3.5400_DP
          vdwparams%dcoeff(3)=23.0025_DP
          vdwparams%c6coeff(1)=2.9865_DP
          vdwparams%c6coeff(6)=27.3784_DP
          vdwparams%c6coeff(7)=19.5223_DP
          vdwparams%c6coeff(8)=11.7733_DP
          vdwparams%c6coeff(15)=190.5033_DP
          vdwparams%c6coeff(16)=161.0388_DP
          vdwparams%radzero(1)=2.8996_DP
          vdwparams%radzero(2)=3.0000_DP
          vdwparams%radzero(3)=3.8000_DP
          vdwparams%radzero(4)=3.8000_DP
          vdwparams%radzero(5)=3.8000_DP
          vdwparams%radzero(6)=3.0001_DP
          vdwparams%radzero(7)=3.2659_DP
          vdwparams%radzero(8)=3.6630_DP
          vdwparams%radzero(9)=3.8000_DP
          vdwparams%radzero(10)=3.8000_DP
          vdwparams%radzero(11)=4.8000_DP
          vdwparams%radzero(12)=4.8000_DP
          vdwparams%radzero(13)=4.8000_DP
          vdwparams%radzero(14)=4.8000_DP
          vdwparams%radzero(15)=4.8000_DP
          vdwparams%radzero(16)=4.7948_DP

       end select pbedamp

    case('RPBE')

       rpbedamp: select case (pub_dispersion)

       case(1) rpbedamp

          vdwparams%dcoeff(1)=3.4967_DP
          vdwparams%dcoeff(2)=3.5400_DP
          vdwparams%dcoeff(3)=23.0000_DP
          vdwparams%c6coeff(1)=3.2054_DP
          vdwparams%c6coeff(6)=27.4608_DP
          vdwparams%c6coeff(7)=19.5765_DP
          vdwparams%c6coeff(8)=11.7926_DP
          vdwparams%c6coeff(15)=190.5033_DP
          vdwparams%c6coeff(16)=161.0465_DP
          vdwparams%radzero(1)=2.8062_DP
          vdwparams%radzero(2)=3.0000_DP
          vdwparams%radzero(3)=3.8000_DP
          vdwparams%radzero(4)=3.8000_DP
          vdwparams%radzero(5)=3.8000_DP
          vdwparams%radzero(6)=3.0375_DP
          vdwparams%radzero(7)=3.2484_DP
          vdwparams%radzero(8)=3.6109_DP
          vdwparams%radzero(9)=3.8000_DP
          vdwparams%radzero(10)=3.8000_DP
          vdwparams%radzero(11)=4.8000_DP
          vdwparams%radzero(12)=4.8000_DP
          vdwparams%radzero(13)=4.8000_DP
          vdwparams%radzero(14)=4.8000_DP
          vdwparams%radzero(15)=4.8000_DP
          vdwparams%radzero(16)=4.0826_DP

       case(2) rpbedamp

          vdwparams%dcoeff(1)=3.0000_DP
          vdwparams%dcoeff(2)=3.9308_DP
          vdwparams%dcoeff(3)=23.0000_DP
          vdwparams%c6coeff(1)=3.2290_DP
          vdwparams%c6coeff(6)=27.4707_DP
          vdwparams%c6coeff(7)=19.5922_DP
          vdwparams%c6coeff(8)=11.8007_DP
          vdwparams%c6coeff(15)=190.5033_DP
          vdwparams%c6coeff(16)=161.0388_DP
          vdwparams%radzero(1)=2.9406_DP
          vdwparams%radzero(2)=3.0000_DP
          vdwparams%radzero(3)=3.8000_DP
          vdwparams%radzero(4)=3.8000_DP
          vdwparams%radzero(5)=3.8000_DP
          vdwparams%radzero(6)=3.4751_DP
          vdwparams%radzero(7)=3.6097_DP
          vdwparams%radzero(8)=3.7393_DP
          vdwparams%radzero(9)=3.8000_DP
          vdwparams%radzero(10)=3.8000_DP
          vdwparams%radzero(11)=4.8000_DP
          vdwparams%radzero(12)=4.8000_DP
          vdwparams%radzero(13)=4.8000_DP
          vdwparams%radzero(14)=4.8000_DP
          vdwparams%radzero(15)=4.8000_DP
          vdwparams%radzero(16)=4.8007_DP

       case(3) rpbedamp

          vdwparams%dcoeff(1)=3.0000_DP
          vdwparams%dcoeff(2)=3.5400_DP
          vdwparams%dcoeff(3)=23.0033_DP
          vdwparams%c6coeff(1)=3.0210_DP
          vdwparams%c6coeff(6)=27.3845_DP
          vdwparams%c6coeff(7)=19.5208_DP
          vdwparams%c6coeff(8)=11.7723_DP
          vdwparams%c6coeff(15)=190.5033_DP
          vdwparams%c6coeff(16)=161.0437_DP
          vdwparams%radzero(1)=2.8667_DP
          vdwparams%radzero(2)=3.0000_DP
          vdwparams%radzero(3)=3.8000_DP
          vdwparams%radzero(4)=3.8000_DP
          vdwparams%radzero(5)=3.8000_DP
          vdwparams%radzero(6)=2.9914_DP
          vdwparams%radzero(7)=3.2931_DP
          vdwparams%radzero(8)=3.6758_DP
          vdwparams%radzero(9)=3.8000_DP
          vdwparams%radzero(10)=3.8000_DP
          vdwparams%radzero(11)=4.8000_DP
          vdwparams%radzero(12)=4.8000_DP
          vdwparams%radzero(13)=4.8000_DP
          vdwparams%radzero(14)=4.8000_DP
          vdwparams%radzero(15)=4.8000_DP
          vdwparams%radzero(16)=3.9912_DP

       end select rpbedamp

    case('REVPBE')

       revpbedamp: select case (pub_dispersion)

       case(1) revpbedamp
          vdwparams%dcoeff(1)=3.4962_DP
          vdwparams%dcoeff(2)=3.5400_DP
          vdwparams%dcoeff(3)=23.0000_DP
          vdwparams%c6coeff(1)=3.2167_DP
          vdwparams%c6coeff(6)=27.4616_DP
          vdwparams%c6coeff(7)=19.5760_DP
          vdwparams%c6coeff(8)=11.7927_DP
          vdwparams%c6coeff(15)=190.5033_DP
          vdwparams%c6coeff(16)=161.0475_DP
          vdwparams%radzero(1)=2.7985_DP
          vdwparams%radzero(2)=3.0000_DP
          vdwparams%radzero(3)=3.8000_DP
          vdwparams%radzero(4)=3.8000_DP
          vdwparams%radzero(5)=3.8000_DP
          vdwparams%radzero(6)=3.0398_DP
          vdwparams%radzero(7)=3.2560_DP
          vdwparams%radzero(8)=3.6122_DP
          vdwparams%radzero(9)=3.8000_DP
          vdwparams%radzero(10)=3.8000_DP
          vdwparams%radzero(11)=4.8000_DP
          vdwparams%radzero(12)=4.8000_DP
          vdwparams%radzero(13)=4.8000_DP
          vdwparams%radzero(14)=4.8000_DP
          vdwparams%radzero(15)=4.8000_DP
          vdwparams%radzero(16)=3.9811_DP

       case(2) revpbedamp

          vdwparams%dcoeff(1)=3.0000_DP
          vdwparams%dcoeff(2)=3.8282_DP
          vdwparams%dcoeff(3)=23.0000_DP
          vdwparams%c6coeff(1)=3.1160_DP
          vdwparams%c6coeff(6)=27.4184_DP
          vdwparams%c6coeff(7)=19.5513_DP
          vdwparams%c6coeff(8)=11.7867_DP
          vdwparams%c6coeff(15)=190.5033_DP
          vdwparams%c6coeff(16)=161.0389_DP
          vdwparams%radzero(1)=2.9540_DP
          vdwparams%radzero(2)=3.0000_DP
          vdwparams%radzero(3)=3.8000_DP
          vdwparams%radzero(4)=3.8000_DP
          vdwparams%radzero(5)=3.8000_DP
          vdwparams%radzero(6)=3.5630_DP
          vdwparams%radzero(7)=3.6696_DP
          vdwparams%radzero(8)=3.7581_DP
          vdwparams%radzero(9)=3.8000_DP
          vdwparams%radzero(10)=3.8000_DP
          vdwparams%radzero(11)=4.8000_DP
          vdwparams%radzero(12)=4.8000_DP
          vdwparams%radzero(13)=4.8000_DP
          vdwparams%radzero(14)=4.8000_DP
          vdwparams%radzero(15)=4.8000_DP
          vdwparams%radzero(16)=4.7980_DP

       case(3) revpbedamp

          vdwparams%dcoeff(1)=3.0000_DP
          vdwparams%dcoeff(2)=3.5400_DP
          vdwparams%dcoeff(3)=23.0034_DP
          vdwparams%c6coeff(1)=3.0258_DP
          vdwparams%c6coeff(6)=27.3854_DP
          vdwparams%c6coeff(7)=19.5210_DP
          vdwparams%c6coeff(8)=11.7724_DP
          vdwparams%c6coeff(15)=190.5033_DP
          vdwparams%c6coeff(16)=161.0435_DP
          vdwparams%radzero(1)=2.8623_DP
          vdwparams%radzero(2)=3.0000_DP
          vdwparams%radzero(3)=3.8000_DP
          vdwparams%radzero(4)=3.8000_DP
          vdwparams%radzero(5)=3.8000_DP
          vdwparams%radzero(6)=2.9923_DP
          vdwparams%radzero(7)=3.2952_DP
          vdwparams%radzero(8)=3.6750_DP
          vdwparams%radzero(9)=3.8000_DP
          vdwparams%radzero(10)=3.8000_DP
          vdwparams%radzero(11)=4.8000_DP
          vdwparams%radzero(12)=4.8000_DP
          vdwparams%radzero(13)=4.8000_DP
          vdwparams%radzero(14)=4.8000_DP
          vdwparams%radzero(15)=4.8000_DP
          vdwparams%radzero(16)=3.9910_DP

       end select revpbedamp

    case('PW91')

       pw91damp: select case (pub_dispersion)

       case(1) pw91damp

          vdwparams%dcoeff(1)=3.2106_DP
          vdwparams%dcoeff(2)=3.5400_DP
          vdwparams%dcoeff(3)=23.0000_DP
          vdwparams%c6coeff(1)=2.8701_DP
          vdwparams%c6coeff(6)=27.3422_DP
          vdwparams%c6coeff(7)=19.5030_DP
          vdwparams%c6coeff(8)=11.7670_DP
          vdwparams%c6coeff(15)=190.5033_DP
          vdwparams%c6coeff(16)=161.0346_DP
          vdwparams%radzero(1)=3.0013_DP
          vdwparams%radzero(2)=3.0000_DP
          vdwparams%radzero(3)=3.8000_DP
          vdwparams%radzero(4)=3.8000_DP
          vdwparams%radzero(5)=3.8000_DP
          vdwparams%radzero(6)=3.4423_DP
          vdwparams%radzero(7)=3.5445_DP
          vdwparams%radzero(8)=3.7444_DP
          vdwparams%radzero(9)=3.8000_DP
          vdwparams%radzero(10)=3.8000_DP
          vdwparams%radzero(11)=4.8000_DP
          vdwparams%radzero(12)=4.8000_DP
          vdwparams%radzero(13)=4.8000_DP
          vdwparams%radzero(14)=4.8000_DP
          vdwparams%radzero(15)=4.8000_DP
          vdwparams%radzero(16)=5.4111_DP

       case(2) pw91damp

          vdwparams%dcoeff(1)=3.0000_DP
          vdwparams%dcoeff(2)=3.3102_DP
          vdwparams%dcoeff(3)=23.0000_DP
          vdwparams%c6coeff(1)=2.5834_DP
          vdwparams%c6coeff(6)=27.2743_DP
          vdwparams%c6coeff(7)=19.4633_DP
          vdwparams%c6coeff(8)=11.7472_DP
          vdwparams%c6coeff(15)=190.5033_DP
          vdwparams%c6coeff(16)=161.0233_DP
          vdwparams%radzero(1)=3.0759_DP
          vdwparams%radzero(2)=3.0000_DP
          vdwparams%radzero(3)=3.8000_DP
          vdwparams%radzero(4)=3.8000_DP
          vdwparams%radzero(5)=3.8000_DP
          vdwparams%radzero(6)=3.9644_DP
          vdwparams%radzero(7)=3.8390_DP
          vdwparams%radzero(8)=3.8330_DP
          vdwparams%radzero(9)=3.8000_DP
          vdwparams%radzero(10)=3.8000_DP
          vdwparams%radzero(11)=4.8000_DP
          vdwparams%radzero(12)=4.8000_DP
          vdwparams%radzero(13)=4.8000_DP
          vdwparams%radzero(14)=4.8000_DP
          vdwparams%radzero(15)=4.8000_DP
          vdwparams%radzero(16)=5.6250_DP

       case(3) pw91damp

          vdwparams%dcoeff(1)=3.0000_DP
          vdwparams%dcoeff(2)=3.5400_DP
          vdwparams%dcoeff(3)=23.0004_DP
          vdwparams%c6coeff(1)=2.9252_DP
          vdwparams%c6coeff(6)=27.3608_DP
          vdwparams%c6coeff(7)=19.5150_DP
          vdwparams%c6coeff(8)=11.7706_DP
          vdwparams%c6coeff(15)=190.5033_DP
          vdwparams%c6coeff(16)=161.0381_DP
          vdwparams%radzero(1)=2.9543_DP
          vdwparams%radzero(2)=3.0000_DP
          vdwparams%radzero(3)=3.8000_DP
          vdwparams%radzero(4)=3.8000_DP
          vdwparams%radzero(5)=3.8000_DP
          vdwparams%radzero(6)=3.0785_DP
          vdwparams%radzero(7)=3.3119_DP
          vdwparams%radzero(8)=3.6858_DP
          vdwparams%radzero(9)=3.8000_DP
          vdwparams%radzero(10)=3.8000_DP
          vdwparams%radzero(11)=4.8000_DP
          vdwparams%radzero(12)=4.8000_DP
          vdwparams%radzero(13)=4.8000_DP
          vdwparams%radzero(14)=4.8000_DP
          vdwparams%radzero(15)=4.8000_DP
          vdwparams%radzero(16)=4.9014_DP

       end select pw91damp

    case('BLYP')
       blypdamp: select case (pub_dispersion)

       case(1) blypdamp
          vdwparams%dcoeff(1)=3.4942_DP
          vdwparams%c6coeff(1)=3.2260_DP
          vdwparams%c6coeff(6)=27.4628_DP
          vdwparams%c6coeff(7)=19.5741_DP
          vdwparams%c6coeff(8)=11.7910_DP
          vdwparams%c6coeff(16)=161.0455_DP
          vdwparams%radzero(1)=2.7898_DP
          vdwparams%radzero(6)=3.0362_DP
          vdwparams%radzero(7)=3.2664_DP
          vdwparams%radzero(8)=3.6239_DP
          vdwparams%radzero(16)=4.1775_DP

       case(2) blypdamp
          vdwparams%dcoeff(2)=3.9215_DP
          vdwparams%c6coeff(1)=3.2514_DP
          vdwparams%c6coeff(6)=27.4690_DP
          vdwparams%c6coeff(7)=19.5853_DP
          vdwparams%c6coeff(8)=11.7973_DP
          vdwparams%c6coeff(16)=161.0389_DP
          vdwparams%radzero(1)=2.9339_DP
          vdwparams%radzero(6)=3.4803_DP
          vdwparams%radzero(7)=3.6233_DP
          vdwparams%radzero(8)=3.7464_DP
          vdwparams%radzero(16)=4.7968_DP

       case(3) blypdamp

          vdwparams%dcoeff(3)=23.0035_DP
          vdwparams%c6coeff(1)=3.0293_DP
          vdwparams%c6coeff(6)=27.3861_DP
          vdwparams%c6coeff(7)=19.5207_DP
          vdwparams%c6coeff(8)=11.7721_DP
          vdwparams%c6coeff(16)=161.0430_DP
          vdwparams%radzero(1)=2.8570_DP
          vdwparams%radzero(6)=2.9890_DP
          vdwparams%radzero(7)=3.3017_DP
          vdwparams%radzero(8)=3.6800_DP
          vdwparams%radzero(16)=4.0893_DP

       end select blypdamp

    case('XLYP')
       xlypdamp: select case (pub_dispersion)

       case(1) xlypdamp
          vdwparams%dcoeff(1)=3.4939_DP
          vdwparams%c6coeff(1)=3.2150_DP
          vdwparams%c6coeff(6)=27.4630_DP
          vdwparams%c6coeff(7)=19.5753_DP
          vdwparams%c6coeff(8)=11.7908_DP
          vdwparams%c6coeff(16)=161.0455_DP
          vdwparams%radzero(1)=2.7993_DP
          vdwparams%radzero(6)=3.0332_DP
          vdwparams%radzero(7)=3.2581_DP
          vdwparams%radzero(8)=3.6244_DP
          vdwparams%radzero(16)=4.1801_DP

       case(2) xlypdamp
          vdwparams%dcoeff(2)=4.0506_DP
          vdwparams%c6coeff(1)=3.4379_DP
          vdwparams%c6coeff(6)=27.5653_DP
          vdwparams%c6coeff(7)=19.6598_DP
          vdwparams%c6coeff(8)=11.8190_DP
          vdwparams%c6coeff(16)=161.0374_DP
          vdwparams%radzero(1)=2.9208_DP
          vdwparams%radzero(6)=3.3614_DP
          vdwparams%radzero(7)=3.5352_DP
          vdwparams%radzero(8)=3.7252_DP
          vdwparams%radzero(16)=4.9034_DP

       case(3) xlypdamp

          vdwparams%dcoeff(3)=23.0034_DP
          vdwparams%c6coeff(1)=3.0238_DP
          vdwparams%c6coeff(6)=27.3853_DP
          vdwparams%c6coeff(7)=19.5206_DP
          vdwparams%c6coeff(8)=11.7719_DP
          vdwparams%c6coeff(16)=161.0439_DP
          vdwparams%radzero(1)=2.8619_DP
          vdwparams%radzero(6)=2.9877_DP
          vdwparams%radzero(7)=3.2994_DP
          vdwparams%radzero(8)=3.6813_DP
          vdwparams%radzero(16)=3.9901_DP

       end select xlypdamp

    case default

       ! qoh: The damping coefficient for the three damping functions
       vdwparams%dcoeff(1)=3.0000_DP
       vdwparams%dcoeff(2)=3.5400_DP
       vdwparams%dcoeff(3)=23.0000_DP

       ! qoh:  Values from Wu and Yang in hartree/bohr ^ 6
       vdwparams%c6coeff(1)=2.8450_DP
       vdwparams%c6coeff(6)=27.3200_DP
       vdwparams%c6coeff(7)=19.4800_DP
       vdwparams%c6coeff(8)=11.7600_DP
       ! qoh: C6 from Halgren
       vdwparams%c6coeff(16)=161.0388_DP

       ! qoh: vdw radii from Elstner in Angstrom
       vdwparams%radzero(1)=3.0000_DP
       vdwparams%radzero(2)=3.0000_DP
       vdwparams%radzero(3)=3.8000_DP
       vdwparams%radzero(4)=3.8000_DP
       vdwparams%radzero(5)=3.8000_DP
       vdwparams%radzero(6)=3.8000_DP
       vdwparams%radzero(7)=3.8000_DP
       vdwparams%radzero(8)=3.8000_DP
       vdwparams%radzero(16)=4.8000_DP

    end select

    ! qoh: Unoptimised parameters from Halgren
    vdwparams%c6coeff(9)=6.2413_DP
    vdwparams%c6coeff(15)=190.5033_DP
    vdwparams%c6coeff(17)=103.5612_DP
    vdwparams%c6coeff(35)=201.8972_DP

    vdwparams%radzero(9)=3.09_DP
    vdwparams%radzero(15)=4.8000_DP
    vdwparams%radzero(17)=4.09_DP
    vdwparams%radzero(35)=4.33_DP

    ! qoh: Array containing the Neff from Wu and Yang
    vdwparams%neff(1)=0.5300_DP
    vdwparams%neff(2)=0.0000_DP
    vdwparams%neff(3)=0.0000_DP
    vdwparams%neff(4)=0.0000_DP
    vdwparams%neff(5)=0.0000_DP
    vdwparams%neff(6)=2.0200_DP
    vdwparams%neff(7)=2.5200_DP
    vdwparams%neff(8)=2.6500_DP
    ! qoh: Neff from Halgren
    vdwparams%neff(9)=3.48_DP
    vdwparams%neff(15)=4.5_DP
    vdwparams%neff(16)=4.8_DP
    vdwparams%neff(17)=5.10_DP
    vdwparams%neff(35)=6.00_DP

    call internal_param_override()

  contains

    subroutine internal_param_override()

    !==================================================================!
    ! This subroutine overrides the parameters set above with those    !
    ! specified in the input file.                                     !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in December 2010.                        !
    !==================================================================!

      use constants, only: DP, VERBOSE, stdout
      use rundat, only: pub_output_detail, pub_vdw_dcoeff

      implicit none

      real(kind=DP) :: c6coeff
      real(kind=DP) :: radzero
      real(kind=DP) :: neff
      integer :: nzatom ! atomic number
      integer :: row ! current row
      integer :: ios ! I/O status flag

      if (allocated(vdwparam_override)) then
         do row=1,size(vdwparam_override)
            ! qoh: Initialise variables
            nzatom=0
            c6coeff=-1.0_DP
            radzero=-1.0_DP
            neff=-1.0_DP

            ! qoh: read a row of the override block
            read (vdwparam_override(row),*,iostat=ios) nzatom, c6coeff, &
                 radzero, neff
            ! qoh: Verify content of the row
            if ((nzatom .gt. 109 .or. nzatom .lt. 1) .or. ios .ne. 0) then
               if (pub_output_detail == VERBOSE) then
                  write(stdout,*) "WARNING: Invalid line in VDW_PARAMS block &
                       &ignoring. Line is:"
                  write(stdout,*) vdwparam_override(row)
               end if
               cycle
            end if
            if (c6coeff .ge. 0.0_DP) vdwparams%c6coeff(nzatom) = c6coeff
            if (radzero .gt. epsilon(1.0_DP)) &
                 vdwparams%radzero(nzatom) = radzero
            if (neff .ge. 0.0_DP) vdwparams%neff(nzatom) = neff
         end do
      end if

      if (pub_vdw_dcoeff .ge. 0.0_DP) &
           vdwparams%dcoeff(pub_dispersion) = pub_vdw_dcoeff

    end subroutine internal_param_override

  end subroutine vdwcorrection_initializeparams

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdwcorrection_calculate_energy(elements) !input

    !==================================================================!
    ! This subroutine calculates the dispersion correction to the      !
    ! total energy.                                                    !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in 2007, with some modifications in 2008 !
    ! Modified by Quintin Hill on 29/05/2009 to account for periodic   !
    ! boundary conditions.                                             !
    !==================================================================!

    use comms,           only: pub_on_root, comms_bcast, pub_root_node_id
    use constants,       only: DP, stdout
    use geometry,        only: point, operator(.DOT.), operator(+), &
         operator(*), operator(-)
    use ion,             only: element
    use rundat,          only: pub_dispersion
    use simulation_cell, only: pub_cell

    implicit none

    type(ELEMENT), intent(in)  :: elements(pub_cell%nat)

    ! Internal variables
    integer       :: atom1, atom2 ! atom counters for loops
    integer       :: a1_neighbour,a2_neighbour, a3_neighbour ! Periodic images
    real(kind=DP) :: distance     ! Distance between pairs of atoms
    real(kind=DP) :: sqdist       ! Square distance between pairs of atoms
    real(kind=DP) :: c6coeff      ! The c6coefficient of the pair
    real(kind=DP) :: damping      ! The damping for the pair
    real(kind=DP) :: distfromco   ! distance - cutoff
    real(kind=DP) :: smoothco     ! cutoff smoothing
    type(point)   :: displacement ! Vector between pairs of atoms

    pub_dispersion_energy = 0.0_DP

    if (pub_dispersion /= 0 .and. pub_on_root) then

       call vdwcorrection_warnings(elements)

       call vdwcorrection_initializeparams

       ! qoh: Loop over all distinct pairs of atoms
       at1: do atom1=1,pub_cell%nat
          at2: do atom2=1,atom1-1

             ! qoh: Calculate c6 coefficient
             c6coeff = vdwcorrection_c6(elements(atom1)%atomic_number,&
                  elements(atom2)%atomic_number)

             ! qoh: Loop over periodic images of atom2 in adjacent cells
             a1: do a1_neighbour = -1,1
                a2: do a2_neighbour = -1,1
                   a3: do a3_neighbour = -1,1

                      displacement =  elements(atom1)%centre &
                           - elements(atom2)%centre &
                           - real(a1_neighbour,kind=DP)*pub_cell%a1 &
                           - real(a2_neighbour,kind=DP)*pub_cell%a2 &
                           - real(a3_neighbour,kind=DP)*pub_cell%a3
                      sqdist = displacement .DOT. displacement

                      within_cutoff: if (sqdist < cutoffsq) then

                         distance = sqrt(sqdist)
                         distfromco = distance - cutoff
                         smoothco = 1.0_DP - exp(-1.0_DP * distfromco**2)

                         ! qoh : Get damping function
                         damping = vdwcorrection_damping(&
                              elements(atom1)%atomic_number,&
                              elements(atom2)%atomic_number,distance)

                         ! qoh: distance**6 = sqdist**3

                         pub_dispersion_energy = pub_dispersion_energy &
                              - (c6coeff * damping / sqdist**3)*smoothco

                      end if within_cutoff
                   end do a3
                end do a2
             end do a1

          enddo at2
       enddo at1

       write(stdout,'(a, e15.8,1x,a)') &
            'Dispersion Correction Energy: ',pub_dispersion_energy, 'Hartree'

    end if

    call comms_bcast(pub_root_node_id, pub_dispersion_energy)

  end subroutine vdwcorrection_calculate_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdwcorrection_calculate_forces(vdw_forces,elements)

    !==================================================================!
    ! This subroutine calculates the dispersion correction to the      !
    ! total energy.                                                    !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in 2007, with some modifications in 2008 !
    ! Modified by Quintin Hill on 29/05/2009 to account for periodic   !
    ! boundary conditions.                                             !
    !==================================================================!

    use comms,           only: pub_on_root, comms_bcast, pub_root_node_id
    use geometry,        only: point, operator(.DOT.), operator(+), &
         operator(*), operator(-)
    use ion,             only: element
    use rundat,          only: pub_dispersion
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments

    real(kind=DP), intent(out) :: vdw_forces(3,pub_cell%nat)
    type(ELEMENT), intent(in)  :: elements(pub_cell%nat)

    ! Internal variables

    integer       :: atom1, atom2 ! atom counters for loops
    integer       :: a1_neighbour,a2_neighbour, a3_neighbour ! Periodic images
    real(kind=DP) :: distance     ! Distance between pairs of atoms
    real(kind=DP) :: sqdist       ! Square distance between pairs of atoms
    real(kind=DP) :: c6coeff      ! The c6coefficient of the pair
    real(kind=DP) :: damping      ! The damping for the pair
    real(kind=DP) :: dampingdrv   ! The damping derivative for the pair
    real(kind=DP) :: drvcommon    ! The common part of the derivative
    real(kind=DP) :: distfromco   ! distance - cutoff
    real(kind=DP) :: smoothco     ! cutoff smoothing
    real(kind=DP) :: smoothcodrv  ! cutoff smoothing derivative
    type(point)   :: displacement ! Vector between pairs of atoms

    vdw_forces = 0.0_DP

    if (pub_dispersion /= 0 .and. pub_on_root) then

       call vdwcorrection_initializeparams

       ! qoh: Loop over all atom pairs
       at1: do atom1=1,pub_cell%nat
          at2: do atom2=1,pub_cell%nat

             if ( atom1 .ne. atom2) then

                ! qoh: Calculate c6 coefficient
                c6coeff = vdwcorrection_c6(elements(atom2)%atomic_number,&
                     elements(atom1)%atomic_number)

                ! qoh: Loop over periodic images of atom2 in adjacent cells
                a1: do a1_neighbour = -1,1
                   a2: do a2_neighbour = -1,1
                      a3: do a3_neighbour = -1,1

                         displacement =  elements(atom1)%centre &
                              - elements(atom2)%centre &
                              - real(a1_neighbour,kind=DP)*pub_cell%a1 &
                              - real(a2_neighbour,kind=DP)*pub_cell%a2 &
                              - real(a3_neighbour,kind=DP)*pub_cell%a3
                         sqdist = displacement .DOT. displacement

                         within_cutoff: if (sqdist < cutoffsq) then

                            distance   = sqrt(sqdist)
                            distfromco = distance - cutoff
                            smoothco = 1.0_DP - exp(-1.0_DP * distfromco**2)
                            smoothcodrv = 2.0_DP*distfromco*(1.0_DP - smoothco)

                            ! qoh : Get damping function
                            damping = vdwcorrection_damping(&
                                 elements(atom1)%atomic_number,&
                                 elements(atom2)%atomic_number,distance)

                            dampingdrv = &
                                 vdwcorrection_drvdamping(&
                                 elements(atom2)%atomic_number,&
                                 elements(atom1)%atomic_number,distance)

                            ! qoh: distance**6 = sqdist**3

                            drvcommon = (c6coeff/sqdist**3)*(smoothco*&
                                 (dampingdrv-6.0_DP*damping/sqdist) &
                                 + smoothcodrv*damping)

                            ! cks: This is the *negative* of the derivative
                            ! cks: of the dispersion energy w.r.t. atomic
                            ! cks: coordinates.
                            vdw_forces(1,atom1) = vdw_forces(1,atom1) &
                                 +drvcommon*displacement%x
                            vdw_forces(2,atom1) = vdw_forces(2,atom1) &
                                 +drvcommon*displacement%y
                            vdw_forces(3,atom1) = vdw_forces(3,atom1) &
                                 +drvcommon* displacement%z

                         end if within_cutoff
                      end do a3
                   end do a2
                end do a1

             end if

          enddo at2
       enddo at1
    end if

    ! cks: broadcast VDW forces from root to all nodes
    call comms_bcast(pub_root_node_id, vdw_forces, 3*pub_cell%nat)

  end subroutine vdwcorrection_calculate_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdwcorrection_warnings(elements)

    !==================================================================!
    ! This subroutine warns about the use of unoptimised or unavailable!
    ! dispersion parameters.                                           !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill on 13/02/2009.                           !
    ! Quintin Hill added functional check on 23/02/2009.               !
    !==================================================================!

    use constants,       only: periodic_table_name, stdout, VERBOSE
    use ion,             only: element
    use rundat,          only: xc_functional, pub_output_detail
    use simulation_cell, only: pub_cell

    implicit none

    type(element), intent(in)    :: elements(pub_cell%nat)

    logical, save                :: already_warned = .false.
    integer                      :: atcnum ! atomic number
    ! qoh: Atomic numbers of elements with unoptimised dispersion parameters:
    integer, parameter           :: unoptimised(4) = (/9,15,17,35/)
    ! qoh: Atomic numbers of elements with optimised dispersion parameters:
    integer, parameter           :: optimised(5) = (/1,6,7,8,16/)
    ! qoh: XC functionals with optimised parameters:
    character(len=6), parameter :: xcfoptimised(6) &
         = (/'BLYP  ','PBE   ','PW91  ','REVPBE','RPBE  ','XLYP  '/)

    if ( pub_output_detail /= VERBOSE .or. already_warned ) return
    ! qoh: Loop over all elements and check if each present in this calculation
    do atcnum=1,109
       elpresent: if (any(elements%atomic_number == atcnum) &
            .and. any(xcfoptimised == xc_functional)) then
          if (any(unoptimised == atcnum)) then
             write(stdout,'(a,a2)') 'WARNING: Unoptimised dispersion &
                  &parameters used for ', periodic_table_name(atcnum)
          elseif (.not. any(optimised == atcnum)) then
             write(stdout,'(a,a2)') 'WARNING: No dispersion parameters &
                  &available for ', periodic_table_name(atcnum)
          end if
       end if elpresent
    end do

    if (.not. any(xcfoptimised == xc_functional)) &
         write(stdout,'(2a)') 'WARNING: No optimised dispersion parameters &
         &available for ', trim(xc_functional)

    already_warned = .true.

  end subroutine vdwcorrection_warnings

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdwcorrection_override_alloc(override_block,nrows)

    !==================================================================!
    ! This subroutine allocates an array to hold parameter overides.   !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in December 2010.                        !
    !==================================================================!

    use utils, only: utils_alloc_check

    implicit none

    integer, intent(in) :: nrows ! Number of rows of override data
    character(len=80), intent(in) :: override_block(nrows) ! Override data
    integer             :: ierr ! error flag

    allocate(vdwparam_override(nrows), stat=ierr)
    call utils_alloc_check('vdwcorrection_overide_alloc','vdwparam_override',&
         ierr)

    vdwparam_override = override_block

  end subroutine vdwcorrection_override_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vdwcorrection_override_dealloc()

    !==================================================================!
    ! This subroutine deallocates the array to hold parameter overides.!
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in December 2010.                        !
    !==================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    integer             :: ierr ! error flag

    if (allocated(vdwparam_override)) then
       deallocate(vdwparam_override, stat=ierr)
       call utils_dealloc_check('vdwcorrection_overide_alloc',&
            'vdwparam_override',ierr)
    end if

  end subroutine vdwcorrection_override_dealloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function vdwcorrection_c6(nzatom1,nzatom2)

    !==================================================================!
    ! This function calculates a heteroatomic C_6 coefficient from     !
    ! homoatomic C_6 coefficients using the formula given in Elstner's !
    ! paper (J. Chem. Phys. 114(12), 5149-5155). (Corrected error in   !
    ! the formula appearing in this paper.)                            !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in 2007, with some modifications in 2008 !
    !==================================================================!

    implicit none

    real(kind=DP)       :: vdwcorrection_c6
    integer, intent(in) :: nzatom1 ! Atomic number of atom 1
    integer, intent(in) :: nzatom2 ! Atomic number of atom 2
    real(kind=DP)       :: c61     ! c6 coefficient of atom 1
    real(kind=DP)       :: c62     ! c6 coefficient of atom 2
    real(kind=DP)       :: ne1     ! Effective number of electrons for atom 1
    real(kind=DP)       :: ne2     ! Effective number of electrons for atom 2
    real(kind=DP), parameter :: third = 1.0_DP/3.0_DP

    ! qoh: Set up shorthands

    c61 = vdwparams%c6coeff(nzatom1)
    c62 = vdwparams%c6coeff(nzatom2)
    ne1 = vdwparams%neff(nzatom1)
    ne2 = vdwparams%neff(nzatom2)

    if (c61 .lt. epsilon(1.0_DP) .or. c62 .lt.epsilon(1.0_DP) .or. &
         ne1 .lt. epsilon(1.0_DP) .or. ne2 .lt. epsilon(1.0_DP)) then
       vdwcorrection_c6=0.0_DP
    else
       vdwcorrection_c6 = 2.0_DP * (c61**2 * c62**2 * ne1 * ne2)**(third)/&
            (((c61 * ne2**2)**(third)) + ((c62* ne1**2)**(third)))
    end if

  end function vdwcorrection_c6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function vdwcorrection_damping(nzatom1,nzatom2,separation)

    !==================================================================!
    ! This function calculates the damping function specified by       !
    ! pub_dispersion:                                                  !
    ! (1) Damping function from Elstner                                !
    !     (J. Chem. Phys. 114(12), 5149-5155).                         !
    ! (2) First damping function from Wu and Yang (I)                  !
    !     (J. Chem. Phys. 116(2), 515-524).                            !
    ! (3) Second damping function from Wu and Yang (II)                !
    !     (J. Chem. Phys. 116(2), 515-524).                            !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in 2007, with some modifications in 2008 !
    ! Merged into a single function in July 2008                       !
    !==================================================================!

    use rundat, only: pub_dispersion
    implicit none

    real(kind=DP)             :: vdwcorrection_damping
    integer, intent(in)       :: nzatom1 ! Atomic number of atom 1
    integer, intent(in)       :: nzatom2 ! Atomic number of atom 2
    real(kind=DP), intent(in) :: separation
    real(kind=DP)             :: radzero
    real(kind=DP)             :: expo ! Exponent
    integer, parameter        :: mexpo(2) = (/4,2/)
    integer, parameter        :: nexpo(2) = (/7,3/)
    integer                   :: mexp
    integer                   :: nexp

    radzero = vdwcorrection_radzero(nzatom1,nzatom2)

    select case (pub_dispersion)
    case(1,2)

       mexp = mexpo(pub_dispersion)
       nexp = nexpo(pub_dispersion)

       expo = -vdwparams%dcoeff(pub_dispersion)*(separation/radzero)**nexp
       vdwcorrection_damping = ( 1.0_DP - exp(expo))**mexp

    case(3)

       expo = -vdwparams%dcoeff(3)*((separation/radzero)-1.0_DP)
       vdwcorrection_damping = 1.0_DP/( 1.0_DP + exp(expo))

    case default
       vdwcorrection_damping = 1.0_DP
    end select

  end function vdwcorrection_damping


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function vdwcorrection_drvdamping(nzatom1,nzatom2,separation)

    !==================================================================!
    ! This function calculates the derivative with respect to atomic   !
    !cooridinates of the damping function specified by pub_dispersion: !
    ! (1) Damping funtion from Elstner                                 !
    !     (J. Chem. Phys. 114(12), 5149-5155).                         !
    ! (2) First damping function from Wu and Yang (I)                  !
    !     (J. Chem. Phys. 116(2), 515-524).                            !
    ! (3) Second damping function from Wu and Yang (II)                !
    !     (J. Chem. Phys. 116(2), 515-524).                            !
    ! Note: For simplicity the (r_{A,i}-r_{B_i}) (where i is x,y or z) !
    !       term is omitted here and included in the calling routine.  !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in July 2008.                            !
    !==================================================================!

    use rundat, only: pub_dispersion
    implicit none

    real(kind=DP)             :: vdwcorrection_drvdamping
    integer,       intent(in) :: nzatom1 ! Atomic number of atom 1
    integer,       intent(in) :: nzatom2 ! Atomic number of atom 2
    real(kind=DP), intent(in) :: separation
    real(kind=DP)             :: radzero
    real(kind=DP)             :: expo ! Exponent
    integer, parameter        :: mexpo(2) = (/4,2/)
    integer, parameter        :: nexpo(2) = (/7,3/)
    integer                   :: mexp
    integer                   :: nexp

    radzero = vdwcorrection_radzero(nzatom1,nzatom2)

    select case (pub_dispersion)
    case(1,2)

       mexp = mexpo(pub_dispersion)
       nexp = nexpo(pub_dispersion)

       expo = -vdwparams%dcoeff(pub_dispersion)*(separation/radzero)**nexp

       vdwcorrection_drvdamping = real((mexp*nexp),kind=DP)*exp(expo)* &
            vdwparams%dcoeff(pub_dispersion)*( 1.0_DP - exp(expo))**(mexp-1)*&
            separation**(nexp-2)/radzero**nexp

    case(3)

       expo = -vdwparams%dcoeff(3)*((separation/radzero)-1.0_DP)

       vdwcorrection_drvdamping = vdwparams%dcoeff(3)*exp(expo)/&
            (separation*radzero*( 1.0_DP + exp(expo))**2)

    case default
       vdwcorrection_drvdamping = 0.0_DP
    end select

  end function vdwcorrection_drvdamping

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function vdwcorrection_radzero(nzatom1,nzatom2)

    !==================================================================!
    ! Function to calculate the R_0 for an atom pair. Uses expression  !
    ! found in Elstner's paper (J. Chem. Phys. 114(12), 5149-5155).    !
    !------------------------------------------------------------------!
    ! Written by Quintin Hill in 2007, with some modifications in 2008 !
    !==================================================================!

    use constants, only: DP, ANGSTROM
    implicit none

    real(kind=DP) :: vdwcorrection_radzero
    integer, intent(in) :: nzatom1 ! Atomic number of atom 1
    integer, intent(in) :: nzatom2 ! Atomic number of atom 2

    vdwcorrection_radzero = ANGSTROM*(vdwparams%radzero(nzatom1)**3 + &
         vdwparams%radzero(nzatom2)**3)/&
         (vdwparams%radzero(nzatom1)**2 + vdwparams%radzero(nzatom2)**2)

  end function vdwcorrection_radzero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module vdwcorrection
