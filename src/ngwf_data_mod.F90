! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                     NGWFs data module                          !
!                                                                !
! This module contains data for NGWF initialisation.             !
!----------------------------------------------------------------!
! Written by Chris-Kriton Skylaris on 06/09/2006.                !
!================================================================!


module ngwf_data

  use constants, only: DP

  implicit none
  private


  ! CKS: PUBLIC FUNCTIONS AND SUBROUTINES
  public :: ngwf_data_cocos_and_expos01_20
  public :: ngwf_data_cocos_and_expos21_40
  public :: ngwf_data_cocos_and_expos41_60
  public :: ngwf_data_cocos_and_expos61_80
  public :: ngwf_data_cocos_and_expos81_103


  ! cks: Define type for Gaussian-Type-Orbitals (GTO) set of each atom
  type GTO_SET

     ! cks: number of shells of Gaussian atomic set
     integer :: nshells

     ! cks: angular momentum of each shell
     integer, pointer, dimension(:) :: angmom

     ! cks: exponents of Gaussian primitives
     real(kind =DP), pointer, dimension(:,:) :: expo

     ! cks: contraction coefficients of Gaussian primitives
     real(kind =DP), pointer, dimension(:,:) :: coco

  end type GTO_SET

  ! cks: make gto_set type available outside the module
  public GTO_SET

  ! cks: number of elements supported by built-in STO-3G basis
  integer, parameter :: pub_num_gtatoms =103

  public :: pub_num_gtatoms

#ifdef OLD_CONSTANTS
  ! ndmh: compile with -DOLD_CONSTANTS to use previous definitions of
  ! ndmh: STO-3G in terms of single-precision reals.
  integer, parameter :: SP=kind(1.0)
#define DP SP
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_data_cocos_and_expos01_20(gbasis)

    !===============================================================!
    ! Initialisation of Gaussian basis set exponents and contraction!
    ! coefficients. All the data of the STO-3G basis set are        !
    ! included for each element of the periodic table as well as    !
    ! the one "polarisation" function of the 6-31G* set.            !
    !---------------------------------------------------------------!
    ! This subroutine was created by Chris-Kriton Skylaris on       !
    ! 31/03/2006 but the actual STO-3G and 6-31G* parameters        !
    ! which appear below were typed by Shyong Chen, in March 2006.  !
    !===============================================================!

    implicit none

    type(GTO_SET), intent(inout) :: gbasis(1: pub_num_gtatoms)


    ! cks: all data in this subroutine was typed by Shyong Chen

    ! ====begin H =========================================
    gbasis(1)%angmom(1)=0
    gbasis(1)%expo(1,1)=3.425250910_DP
    gbasis(1)%expo(1,2)=0.623913730_DP
    gbasis(1)%expo(1,3)=0.168855400_DP
    gbasis(1)%coco(1,1)=0.15432897_DP
    gbasis(1)%coco(1,2)=0.53532814_DP
    gbasis(1)%coco(1,3)=0.44463454_DP

    gbasis(1)%angmom(2)=1
    gbasis(1)%expo(2,1)=0.161277800_DP
    gbasis(1)%coco(2,1)=1.00000000_DP
    ! ====end H ============================================


    ! ====begin He ============================================
    gbasis(2)%angmom(1)=0
    gbasis(2)%expo(1,1)=6.362421390_DP
    gbasis(2)%expo(1,2)=1.158923000_DP
    gbasis(2)%expo(1,3)=0.313649790_DP
    gbasis(2)%coco(1,1)=0.15432897_DP
    gbasis(2)%coco(1,2)=0.53532814_DP
    gbasis(2)%coco(1,3)=0.44463454_DP

    gbasis(2)%angmom(2)=1
    gbasis(2)%expo(2,1)=0.297964000_DP
    gbasis(2)%coco(2,1)=1.00000000_DP
    ! ====end He ============================================


    ! ====begin Li ==========================================
    gbasis(3)%angmom(1)=0
    gbasis(3)%expo(1,1)=16.119575000_DP
    gbasis(3)%expo(1,2)=2.936200700_DP
    gbasis(3)%expo(1,3)=0.794650500_DP
    gbasis(3)%coco(1,1)=0.15432897_DP
    gbasis(3)%coco(1,2)=0.53532814_DP
    gbasis(3)%coco(1,3)=0.44463454_DP

    gbasis(3)%angmom(2)=0
    gbasis(3)%expo(2,1)=0.636289700_DP
    gbasis(3)%expo(2,2)=0.147860100_DP
    gbasis(3)%expo(2,3)=0.048088700_DP
    gbasis(3)%coco(2,1)=-0.09996723_DP
    gbasis(3)%coco(2,2)=0.39951283_DP
    gbasis(3)%coco(2,3)=0.70011547_DP

    gbasis(3)%angmom(3)=1
    gbasis(3)%expo(3,1)=0.636289700_DP
    gbasis(3)%expo(3,2)=0.147860100_DP
    gbasis(3)%expo(3,3)=0.048088700_DP
    gbasis(3)%coco(3,1)=0.15591627_DP
    gbasis(3)%coco(3,2)=0.60768372_DP
    gbasis(3)%coco(3,3)=0.39195739_DP

    gbasis(3)%angmom(4)=2
    gbasis(3)%expo(4,1)=0.20000000_DP
    gbasis(3)%coco(4,1)=1.00000000_DP
    ! ====end Li ==========================================


    ! ====begin Be ==========================================
    gbasis(4)%angmom(1)=0
    gbasis(4)%expo(1,1)=30.167871000_DP
    gbasis(4)%expo(1,2)=5.495115300_DP
    gbasis(4)%expo(1,3)=1.487192700_DP
    gbasis(4)%coco(1,1)=0.15432897_DP
    gbasis(4)%coco(1,2)=0.53532814_DP
    gbasis(4)%coco(1,3)=0.44463454_DP

    gbasis(4)%angmom(2)=0
    gbasis(4)%expo(2,1)=1.314833100_DP
    gbasis(4)%expo(2,2)=0.305538900_DP
    gbasis(4)%expo(2,3)=0.099370700_DP
    gbasis(4)%coco(2,1)=-0.09996723_DP
    gbasis(4)%coco(2,2)=0.39951283_DP
    gbasis(4)%coco(2,3)=0.70011547_DP

    gbasis(4)%angmom(3)=1
    gbasis(4)%expo(3,1)=1.314833100_DP
    gbasis(4)%expo(3,2)=0.305538900_DP
    gbasis(4)%expo(3,3)=0.099370700_DP
    gbasis(4)%coco(3,1)=0.15591627_DP
    gbasis(4)%coco(3,2)=0.60768372_DP
    gbasis(4)%coco(3,3)=0.39195739_DP

    gbasis(4)%angmom(4)=2
    gbasis(4)%expo(4,1)=0.40000000_DP
    gbasis(4)%coco(4,1)=1.00000000_DP
    ! ====end Be ==========================================

    ! ====begin B ==========================================
    gbasis(5)%angmom(1)=0
    gbasis(5)%expo(1,1)=48.791113000_DP
    gbasis(5)%expo(1,2)=8.887362200_DP
    gbasis(5)%expo(1,3)=2.405267000_DP
    gbasis(5)%coco(1,1)=0.15432897_DP
    gbasis(5)%coco(1,2)=0.53532814_DP
    gbasis(5)%coco(1,3)=0.44463454_DP

    gbasis(5)%angmom(2)=0
    gbasis(5)%expo(2,1)=2.236956100_DP
    gbasis(5)%expo(2,2)=0.519820500_DP
    gbasis(5)%expo(2,3)=0.169061800_DP
    gbasis(5)%coco(2,1)=-0.09996723_DP
    gbasis(5)%coco(2,2)=0.39951283_DP
    gbasis(5)%coco(2,3)=0.70011547_DP

    gbasis(5)%angmom(3)=1
    gbasis(5)%expo(3,1)=2.236956100_DP
    gbasis(5)%expo(3,2)=0.519820500_DP
    gbasis(5)%expo(3,3)=0.169061800_DP
    gbasis(5)%coco(3,1)=0.15591627_DP
    gbasis(5)%coco(3,2)=0.60768372_DP
    gbasis(5)%coco(3,3)=0.39195739_DP

    gbasis(5)%angmom(4)=2
    gbasis(5)%expo(4,1)=0.60000000_DP
    gbasis(5)%coco(4,1)=1.00000000_DP
    ! ====end B ==========================================


    ! ====begin C ==========================================
    gbasis(6)%angmom(1)=0
    gbasis(6)%expo(1,1)=71.616837000_DP
    gbasis(6)%expo(1,2)=13.045096000_DP
    gbasis(6)%expo(1,3)=3.530512200_DP
    gbasis(6)%coco(1,1)=0.15432897_DP
    gbasis(6)%coco(1,2)=0.53532814_DP
    gbasis(6)%coco(1,3)=0.44463454_DP

    gbasis(6)%angmom(2)=0
    gbasis(6)%expo(2,1)=2.941249400_DP
    gbasis(6)%expo(2,2)=0.683483100_DP
    gbasis(6)%expo(2,3)=0.222289900_DP
    gbasis(6)%coco(2,1)=-0.09996723_DP
    gbasis(6)%coco(2,2)=0.39951283_DP
    gbasis(6)%coco(2,3)=0.70011547_DP

    gbasis(6)%angmom(3)=1
    gbasis(6)%expo(3,1)=2.941249400_DP
    gbasis(6)%expo(3,2)=0.683483100_DP
    gbasis(6)%expo(3,3)=0.222289900_DP
    gbasis(6)%coco(3,1)=0.15591627_DP
    gbasis(6)%coco(3,2)=0.60768372_DP
    gbasis(6)%coco(3,3)=0.39195739_DP

    gbasis(6)%angmom(4)=2
    gbasis(6)%expo(4,1)=0.80000000_DP
    gbasis(6)%coco(4,1)=1.00000000_DP
    ! ====end C ==========================================


    ! ====begin N ==========================================
    gbasis(7)%angmom(1)=0
    gbasis(7)%expo(1,1)=99.106169000_DP
    gbasis(7)%expo(1,1)=18.052312000_DP
    gbasis(7)%expo(1,1)=4.885660200_DP
    gbasis(7)%coco(1,1)=0.15432897_DP
    gbasis(7)%coco(1,2)=0.53532814_DP
    gbasis(7)%coco(1,3)=0.44463454_DP

    gbasis(7)%angmom(2)=0
    gbasis(7)%expo(2,1)=3.780455900_DP
    gbasis(7)%expo(2,2)=0.878496600_DP
    gbasis(7)%expo(2,3)=0.285714400_DP
    gbasis(7)%coco(2,1)=-0.09996723_DP
    gbasis(7)%coco(2,2)=0.39951283_DP
    gbasis(7)%coco(2,3)=0.70011547_DP

    gbasis(7)%angmom(3)=1
    gbasis(7)%expo(3,1)=3.780455900_DP
    gbasis(7)%expo(3,2)=0.878496600_DP
    gbasis(7)%expo(3,3)=0.285714400_DP
    gbasis(7)%coco(3,1)=0.15591627_DP
    gbasis(7)%coco(3,2)=0.60768372_DP
    gbasis(7)%coco(3,3)=0.39195739_DP

    gbasis(7)%angmom(4)=2
    gbasis(7)%expo(4,1)=0.80000000_DP
    gbasis(7)%coco(4,1)=1.00000000_DP
    ! ====end N ==========================================


    ! ====begin O ==========================================
    gbasis(8)%angmom(1)=0
    gbasis(8)%expo(1,1)=130.709320000_DP
    gbasis(8)%expo(1,2)=23.808861000_DP
    gbasis(8)%expo(1,3)=6.443608300_DP
    gbasis(8)%coco(1,1)=0.15432897_DP
    gbasis(8)%coco(1,2)=0.53532814_DP
    gbasis(8)%coco(1,3)=0.44463454_DP

    gbasis(8)%angmom(2)=0
    gbasis(8)%expo(2,1)=5.033151300_DP
    gbasis(8)%expo(2,2)=1.169596100_DP
    gbasis(8)%expo(2,3)=0.380389000_DP
    gbasis(8)%coco(2,1)=-0.09996723_DP
    gbasis(8)%coco(2,2)=0.39951283_DP
    gbasis(8)%coco(2,3)=0.70011547_DP

    gbasis(8)%angmom(3)=1
    gbasis(8)%expo(3,1)=5.033151300_DP
    gbasis(8)%expo(3,2)=1.169596100_DP
    gbasis(8)%expo(3,3)=0.380389000_DP
    gbasis(8)%coco(3,1)=0.15591627_DP
    gbasis(8)%coco(3,2)=0.60768372_DP
    gbasis(8)%coco(3,3)=0.39195739_DP

    gbasis(8)%angmom(4)=2
    gbasis(8)%expo(4,1)=0.80000000_DP
    gbasis(8)%coco(4,1)=1.00000000_DP
    ! ====end O ==========================================


    ! ====begin F ==========================================
    gbasis(9)%angmom(1)=0
    gbasis(9)%expo(1,1)=166.679130000_DP
    gbasis(9)%expo(1,2)=30.360812000_DP
    gbasis(9)%expo(1,3)=8.216820700_DP
    gbasis(9)%coco(1,1)=0.15432897_DP
    gbasis(9)%coco(1,2)=0.53532814_DP
    gbasis(9)%coco(1,3)=0.44463454_DP

    gbasis(9)%angmom(2)=0
    gbasis(9)%expo(2,1)=6.464803200_DP
    gbasis(9)%expo(2,2)=1.502281200_DP
    gbasis(9)%expo(2,3)=0.488588500_DP
    gbasis(9)%coco(2,1)=-0.09996723_DP
    gbasis(9)%coco(2,2)=0.39951283_DP
    gbasis(9)%coco(2,3)=0.70011547_DP

    gbasis(9)%angmom(3)=1
    gbasis(9)%expo(3,1)=6.464803200_DP
    gbasis(9)%expo(3,2)=1.502281200_DP
    gbasis(9)%expo(3,3)=0.488588500_DP
    gbasis(9)%coco(3,1)=0.15591627_DP
    gbasis(9)%coco(3,2)=0.60768372_DP
    gbasis(9)%coco(3,3)=0.39195739_DP

    gbasis(9)%angmom(4)=2
    gbasis(9)%expo(4,1)=0.80000000_DP
    gbasis(9)%coco(4,1)=1.00000000_DP
    ! ====end F ==========================================


    ! ====begin Ne ==========================================
    gbasis(10)%angmom(1)=0
    gbasis(10)%expo(1,1)=207.015610000_DP
    gbasis(10)%expo(1,2)=37.708151000_DP
    gbasis(10)%expo(1,3)=10.205297000_DP
    gbasis(10)%coco(1,1)=0.15432897_DP
    gbasis(10)%coco(1,2)=0.53532814_DP
    gbasis(10)%coco(1,3)=0.44463454_DP

    gbasis(10)%angmom(2)=0
    gbasis(10)%expo(2,1)=8.246315100_DP
    gbasis(10)%expo(2,2)=1.916266200_DP
    gbasis(10)%expo(2,3)=0.623229300_DP
    gbasis(10)%coco(2,1)=-0.09996723_DP
    gbasis(10)%coco(2,2)=0.39951283_DP
    gbasis(10)%coco(2,3)=0.70011547_DP

    gbasis(10)%angmom(3)=1
    gbasis(10)%expo(3,1)=8.246315100_DP
    gbasis(10)%expo(3,2)=1.916266200_DP
    gbasis(10)%expo(3,3)=0.623229300_DP
    gbasis(10)%coco(3,1)=0.15591627_DP
    gbasis(10)%coco(3,2)=0.60768372_DP
    gbasis(10)%coco(3,3)=0.39195739_DP

    gbasis(10)%angmom(4)=2
    gbasis(10)%expo(4,1)=0.80000000_DP
    gbasis(10)%coco(4,1)=1.00000000_DP
    ! ====end Ne ==========================================


    ! ====begin Na ==========================================
    gbasis(11)%angmom(1)=0
    gbasis(11)%expo(1,1)=250.772430000_DP
    gbasis(11)%expo(1,2)=45.678511000_DP
    gbasis(11)%expo(1,3)=12.362388000_DP
    gbasis(11)%coco(1,1)=0.15432897_DP
    gbasis(11)%coco(1,2)=0.53532814_DP
    gbasis(11)%coco(1,3)=0.44463454_DP

    gbasis(11)%angmom(2)=0
    gbasis(11)%expo(2,1)=12.040193000_DP
    gbasis(11)%expo(2,2)=2.797881900_DP
    gbasis(11)%expo(2,3)=0.909958000_DP
    gbasis(11)%coco(2,1)=-0.09996723_DP
    gbasis(11)%coco(2,2)=0.39951283_DP
    gbasis(11)%coco(2,3)=0.70011547_DP

    gbasis(11)%angmom(3)=1
    gbasis(11)%expo(3,1)=12.040193000_DP
    gbasis(11)%expo(3,2)=2.797881900_DP
    gbasis(11)%expo(3,3)=0.909958000_DP
    gbasis(11)%coco(3,1)=0.15591627_DP
    gbasis(11)%coco(3,2)=0.60768372_DP
    gbasis(11)%coco(3,3)=0.39195739_DP

    gbasis(11)%angmom(4)=0
    gbasis(11)%expo(4,1)=1.478740600_DP
    gbasis(11)%expo(4,2)=0.412564900_DP
    gbasis(11)%expo(4,3)=0.161475100_DP
    gbasis(11)%coco(4,1)=-0.21962037_DP
    gbasis(11)%coco(4,2)=0.22559543_DP
    gbasis(11)%coco(4,3)=0.90039843_DP

    gbasis(11)%angmom(5)=1
    gbasis(11)%expo(5,1)=1.478740600_DP
    gbasis(11)%expo(5,2)=0.412564900_DP
    gbasis(11)%expo(5,3)=0.161475100_DP
    gbasis(11)%coco(5,1)=0.01058760_DP
    gbasis(11)%coco(5,2)=0.59516701_DP
    gbasis(11)%coco(5,3)=0.46200101_DP

    gbasis(11)%angmom(6)=2
    gbasis(11)%expo(6,1)=0.17500000_DP
    gbasis(11)%coco(6,1)=1.00000000_DP
    ! ====end Na ==========================================


    ! ====begin Mg ==========================================
    gbasis(12)%angmom(1)=0
    gbasis(12)%expo(1,1)=299.237400000_DP
    gbasis(12)%expo(1,2)=54.506470000_DP
    gbasis(12)%expo(1,3)=14.751580000_DP
    gbasis(12)%coco(1,1)=0.15432897_DP
    gbasis(12)%coco(1,2)=0.53532814_DP
    gbasis(12)%coco(1,3)=0.44463454_DP

    gbasis(12)%angmom(2)=0
    gbasis(12)%expo(2,1)=15.121820000_DP
    gbasis(12)%expo(2,2)=3.513987000_DP
    gbasis(12)%expo(2,3)=1.142857000_DP
    gbasis(12)%coco(2,1)=-0.09996723_DP
    gbasis(12)%coco(2,2)=0.39951283_DP
    gbasis(12)%coco(2,3)=0.70011547_DP

    gbasis(12)%angmom(3)=1
    gbasis(12)%expo(3,1)=15.121820000_DP
    gbasis(12)%expo(3,2)=3.513987000_DP
    gbasis(12)%expo(3,3)=1.142857000_DP
    gbasis(12)%coco(3,1)=0.15591628_DP
    gbasis(12)%coco(3,2)=0.60768372_DP
    gbasis(12)%coco(3,3)=0.39195739_DP

    gbasis(12)%angmom(4)=0
    gbasis(12)%expo(4,1)=1.395448000_DP
    gbasis(12)%expo(4,2)=0.389326000_DP
    gbasis(12)%expo(4,3)=0.152380000_DP
    gbasis(12)%coco(4,1)=-0.21962037_DP
    gbasis(12)%coco(4,2)=0.22559543_DP
    gbasis(12)%coco(4,3)=0.90039843_DP

    gbasis(12)%angmom(5)=1
    gbasis(12)%expo(5,1)=1.395448000_DP
    gbasis(12)%expo(5,2)=0.389326000_DP
    gbasis(12)%expo(5,3)=0.152380000_DP
    gbasis(12)%coco(5,1)=0.01058760_DP
    gbasis(12)%coco(5,2)=0.59516701_DP
    gbasis(12)%coco(5,3)=0.46200101_DP

    gbasis(12)%angmom(6)=2
    gbasis(12)%expo(6,1)=0.17500000_DP
    gbasis(12)%coco(6,1)=1.00000000_DP
    ! ====end Mg ==========================================


    ! ====begin Al ==========================================
    gbasis(13)%angmom(1)=0
    gbasis(13)%expo(1,1)=351.421476700_DP
    gbasis(13)%expo(1,2)=64.011860670_DP
    gbasis(13)%expo(1,3)=17.324107610_DP
    gbasis(13)%coco(1,1)=0.15432897_DP
    gbasis(13)%coco(1,2)=0.53532814_DP
    gbasis(13)%coco(1,3)=0.44463454_DP

    gbasis(13)%angmom(2)=0
    gbasis(13)%expo(2,1)=18.899396210_DP
    gbasis(13)%expo(2,2)=4.391813233_DP
    gbasis(13)%expo(2,3)=1.428353970_DP
    gbasis(13)%coco(2,1)=-0.09996723_DP
    gbasis(13)%coco(2,2)=0.39951283_DP
    gbasis(13)%coco(2,3)=0.70011547_DP

    gbasis(13)%angmom(3)=1
    gbasis(13)%expo(3,1)=18.899396210_DP
    gbasis(13)%expo(3,2)=4.391813233_DP
    gbasis(13)%expo(3,3)=1.428353970_DP
    gbasis(13)%coco(3,1)=0.15591628_DP
    gbasis(13)%coco(3,2)=0.60768372_DP
    gbasis(13)%coco(3,3)=0.39195739_DP

    gbasis(13)%angmom(4)=0
    gbasis(13)%expo(4,1)=1.395448293_DP
    gbasis(13)%expo(4,2)=0.389326532_DP
    gbasis(13)%expo(4,3)=0.152379766_DP
    gbasis(13)%coco(4,1)=-0.21962037_DP
    gbasis(13)%coco(4,2)=0.22559543_DP
    gbasis(13)%coco(4,3)=0.90039843_DP

    gbasis(13)%angmom(5)=1
    gbasis(13)%expo(5,1)=1.395448293_DP
    gbasis(13)%expo(5,2)=0.389326532_DP
    gbasis(13)%expo(5,3)=0.152379766_DP
    gbasis(13)%coco(5,1)=0.01058760_DP
    gbasis(13)%coco(5,2)=0.59516701_DP
    gbasis(13)%coco(5,3)=0.46200101_DP

    gbasis(13)%angmom(6)=2
    gbasis(13)%expo(6,1)=0.32500000_DP
    gbasis(13)%coco(6,1)=1.00000000_DP
    ! ====end Al ==========================================


    ! ====begin Si ==========================================
    gbasis(14)%angmom(1)=0
    gbasis(14)%expo(1,1)=407.797551400_DP
    gbasis(14)%expo(1,2)=74.280833050_DP
    gbasis(14)%expo(1,3)=20.103292290_DP
    gbasis(14)%coco(1,1)=0.15432897_DP
    gbasis(14)%coco(1,2)=0.53532814_DP
    gbasis(14)%coco(1,3)=0.44463454_DP

    gbasis(14)%angmom(2)=0
    gbasis(14)%expo(2,1)=23.193656060_DP
    gbasis(14)%expo(2,2)=5.389706871_DP
    gbasis(14)%expo(2,3)=1.752899952_DP
    gbasis(14)%coco(2,1)=-0.09996723_DP
    gbasis(14)%coco(2,2)=0.39951283_DP
    gbasis(14)%coco(2,3)=0.70011547_DP

    gbasis(14)%angmom(3)=1
    gbasis(14)%expo(3,1)=23.193656060_DP
    gbasis(14)%expo(3,2)=5.389706871_DP
    gbasis(14)%expo(3,3)=1.752899952_DP
    gbasis(14)%coco(3,1)=0.15591628_DP
    gbasis(14)%coco(3,2)=0.60768372_DP
    gbasis(14)%coco(3,3)=0.39195739_DP

    gbasis(14)%angmom(4)=0
    gbasis(14)%expo(4,1)=1.478740622_DP
    gbasis(14)%expo(4,2)=0.412564880_DP
    gbasis(14)%expo(4,3)=0.161475098_DP
    gbasis(14)%coco(4,1)=-0.21962037_DP
    gbasis(14)%coco(4,2)=0.22559543_DP
    gbasis(14)%coco(4,3)=0.90039843_DP

    gbasis(14)%angmom(5)=1
    gbasis(14)%expo(5,1)=1.478740622_DP
    gbasis(14)%expo(5,2)=0.412564880_DP
    gbasis(14)%expo(5,3)=0.161475098_DP
    gbasis(14)%coco(5,1)=0.01058760_DP
    gbasis(14)%coco(5,2)=0.59516701_DP
    gbasis(14)%coco(5,3)=0.46200101_DP

    gbasis(14)%angmom(6)=2
    gbasis(14)%expo(6,1)=0.45000000_DP
    gbasis(14)%coco(6,1)=1.00000000_DP
    ! ====end Si ==========================================


    ! ====begin P ==========================================
    gbasis(15)%angmom(1)=0
    gbasis(15)%expo(1,1)=468.365637800_DP
    gbasis(15)%expo(1,2)=85.313385590_DP
    gbasis(15)%expo(1,3)=23.089131560_DP
    gbasis(15)%coco(1,1)=0.15432897_DP
    gbasis(15)%coco(1,2)=0.53532814_DP
    gbasis(15)%coco(1,3)=0.44463454_DP

    gbasis(15)%angmom(2)=0
    gbasis(15)%expo(2,1)=28.032639580_DP
    gbasis(15)%expo(2,2)=6.514182577_DP
    gbasis(15)%expo(2,3)=2.118614352_DP
    gbasis(15)%coco(2,1)=-0.09996723_DP
    gbasis(15)%coco(2,2)=0.39951283_DP
    gbasis(15)%coco(2,3)=0.70011547_DP

    gbasis(15)%angmom(3)=1
    gbasis(15)%expo(3,1)=28.032639580_DP
    gbasis(15)%expo(3,2)=6.514182577_DP
    gbasis(15)%expo(3,3)=2.118614352_DP
    gbasis(15)%coco(3,1)=0.15591628_DP
    gbasis(15)%coco(3,2)=0.60768372_DP
    gbasis(15)%coco(3,3)=0.39195739_DP

    gbasis(15)%angmom(4)=0
    gbasis(15)%expo(4,1)=1.743103231_DP
    gbasis(15)%expo(4,2)=0.486321377_DP
    gbasis(15)%expo(4,3)=0.190342891_DP
    gbasis(15)%coco(4,1)=-0.21962037_DP
    gbasis(15)%coco(4,2)=0.22559543_DP
    gbasis(15)%coco(4,3)=0.90039843_DP

    gbasis(15)%angmom(5)=1
    gbasis(15)%expo(5,1)=1.743103231_DP
    gbasis(15)%expo(5,2)=0.486321377_DP
    gbasis(15)%expo(5,3)=0.190342891_DP
    gbasis(15)%coco(5,1)=0.01058760_DP
    gbasis(15)%coco(5,2)=0.59516701_DP
    gbasis(15)%coco(5,3)=0.46200101_DP

    gbasis(15)%angmom(6)=2
    gbasis(15)%expo(6,1)=0.55000000_DP
    gbasis(15)%coco(6,1)=1.00000000_DP
    ! ====end P ==========================================


    ! ====begin S ==========================================
    gbasis(16)%angmom(1)=0
    gbasis(16)%expo(1,1)=533.125735900_DP
    gbasis(16)%expo(1,2)=97.109518300_DP
    gbasis(16)%expo(1,3)=26.281625420_DP
    gbasis(16)%coco(1,1)=0.15432897_DP
    gbasis(16)%coco(1,2)=0.53532814_DP
    gbasis(16)%coco(1,3)=0.44463454_DP

    gbasis(16)%angmom(2)=0
    gbasis(16)%expo(2,1)=33.329751730_DP
    gbasis(16)%expo(2,2)=7.745117521_DP
    gbasis(16)%expo(2,3)=2.518952599_DP
    gbasis(16)%coco(2,1)=-0.09996723_DP
    gbasis(16)%coco(2,2)=0.39951283_DP
    gbasis(16)%coco(2,3)=0.70011547_DP

    gbasis(16)%angmom(3)=1
    gbasis(16)%expo(3,1)=33.329751730_DP
    gbasis(16)%expo(3,2)=7.745117521_DP
    gbasis(16)%expo(3,3)=2.518952599_DP
    gbasis(16)%coco(3,1)=0.15591628_DP
    gbasis(16)%coco(3,2)=0.60768372_DP
    gbasis(16)%coco(3,3)=0.39195739_DP

    gbasis(16)%angmom(4)=0
    gbasis(16)%expo(4,1)=2.029194274_DP
    gbasis(16)%expo(4,2)=0.566140052_DP
    gbasis(16)%expo(4,3)=0.221583379_DP
    gbasis(16)%coco(4,1)=-0.21962037_DP
    gbasis(16)%coco(4,2)=0.22559543_DP
    gbasis(16)%coco(4,3)=0.90039843_DP

    gbasis(16)%angmom(5)=1
    gbasis(16)%expo(5,1)=2.029194274_DP
    gbasis(16)%expo(5,2)=0.566140052_DP
    gbasis(16)%expo(5,3)=0.221583379_DP
    gbasis(16)%coco(5,1)=0.01058760_DP
    gbasis(16)%coco(5,2)=0.59516701_DP
    gbasis(16)%coco(5,3)=0.46200101_DP

    gbasis(16)%angmom(6)=2
    gbasis(16)%expo(6,1)=0.65000000_DP
    gbasis(16)%coco(6,1)=1.00000000_DP
    ! ====end S ============================================


    ! ====begin Cl ========================================
    gbasis(17)%angmom(1)=0
    gbasis(17)%expo(1,1)=601.345613600_DP
    gbasis(17)%expo(1,2)=109.535854200_DP
    gbasis(17)%expo(1,3)=29.644676860_DP
    gbasis(17)%coco(1,1)=0.15432897_DP
    gbasis(17)%coco(1,2)=0.53532814_DP
    gbasis(17)%coco(1,3)=0.44463454_DP

    gbasis(17)%angmom(2)=0
    gbasis(17)%expo(2,1)=38.960418890_DP
    gbasis(17)%expo(2,2)=9.053563477_DP
    gbasis(17)%expo(2,3)=2.944499834_DP
    gbasis(17)%coco(2,1)=-0.09996723_DP
    gbasis(17)%coco(2,2)=0.39951283_DP
    gbasis(17)%coco(2,3)=0.70011547_DP

    gbasis(17)%angmom(3)=1
    gbasis(17)%expo(3,1)=38.960418890_DP
    gbasis(17)%expo(3,2)=9.053563477_DP
    gbasis(17)%expo(3,3)=2.944499834_DP
    gbasis(17)%coco(3,1)=0.15591628_DP
    gbasis(17)%coco(3,2)=0.60768372_DP
    gbasis(17)%coco(3,3)=0.39195739_DP

    gbasis(17)%angmom(4)=0
    gbasis(17)%expo(4,1)=2.129386495_DP
    gbasis(17)%expo(4,2)=0.594093427_DP
    gbasis(17)%expo(4,3)=0.232524141_DP
    gbasis(17)%coco(4,1)=-0.21962037_DP
    gbasis(17)%coco(4,2)=0.22559543_DP
    gbasis(17)%coco(4,3)=0.90039843_DP

    gbasis(17)%angmom(5)=1
    gbasis(17)%expo(5,1)=2.129386495_DP
    gbasis(17)%expo(5,2)=0.594093427_DP
    gbasis(17)%expo(5,3)=0.232524141_DP
    gbasis(17)%coco(5,1)=0.01058760_DP
    gbasis(17)%coco(5,2)=0.59516701_DP
    gbasis(17)%coco(5,3)=0.46200101_DP

    gbasis(17)%angmom(6)=2
    gbasis(17)%expo(6,1)=0.75000000_DP
    gbasis(17)%coco(6,1)=1.00000000_DP
    ! ====end Cl ==========================================


    ! ====begin Ar ==========================================
    gbasis(18)%angmom(1)=0
    gbasis(18)%expo(1,1)=674.446518400_DP
    gbasis(18)%expo(1,2)=122.851275300_DP
    gbasis(18)%expo(1,3)=33.248349450_DP
    gbasis(18)%coco(1,1)=0.15432897_DP
    gbasis(18)%coco(1,2)=0.53532814_DP
    gbasis(18)%coco(1,3)=0.44463454_DP

    gbasis(18)%angmom(2)=0
    gbasis(18)%expo(2,1)=45.164243920_DP
    gbasis(18)%expo(2,2)=10.495199000_DP
    gbasis(18)%expo(2,3)=3.413364448_DP
    gbasis(18)%coco(2,1)=-0.09996723_DP
    gbasis(18)%coco(2,2)=0.39951283_DP
    gbasis(18)%coco(2,3)=0.70011547_DP

    gbasis(18)%angmom(3)=1
    gbasis(18)%expo(3,1)=45.164243920_DP
    gbasis(18)%expo(3,2)=10.495199000_DP
    gbasis(18)%expo(3,3)=3.413364448_DP
    gbasis(18)%coco(3,1)=0.15591628_DP
    gbasis(18)%coco(3,2)=0.60768372_DP
    gbasis(18)%coco(3,3)=0.39195739_DP

    gbasis(18)%angmom(4)=0
    gbasis(18)%expo(4,1)=2.621366518_DP
    gbasis(18)%expo(4,2)=0.731354605_DP
    gbasis(18)%expo(4,3)=0.286247236_DP
    gbasis(18)%coco(4,1)=-0.21962037_DP
    gbasis(18)%coco(4,2)=0.22559543_DP
    gbasis(18)%coco(4,3)=0.90039843_DP

    gbasis(18)%angmom(5)=1
    gbasis(18)%expo(5,1)=2.621366518_DP
    gbasis(18)%expo(5,2)=0.731354605_DP
    gbasis(18)%expo(5,3)=0.286247236_DP
    gbasis(18)%coco(5,1)=0.01058760_DP
    gbasis(18)%coco(5,2)=0.59516701_DP
    gbasis(18)%coco(5,3)=0.46200101_DP

    gbasis(18)%angmom(6)=2
    gbasis(18)%expo(6,1)=0.85000000_DP
    gbasis(18)%coco(6,1)=1.00000000_DP
    ! ====end Ar ==========================================


    ! ====begin K ==========================================
    gbasis(19)%angmom(1)=0
    gbasis(19)%expo(1,1)=771.510368100_DP
    gbasis(19)%expo(1,2)=140.531576600_DP
    gbasis(19)%expo(1,3)=38.033328990_DP
    gbasis(19)%coco(1,1)=0.15432897_DP
    gbasis(19)%coco(1,2)=0.53532814_DP
    gbasis(19)%coco(1,3)=0.44463454_DP

    gbasis(19)%angmom(2)=0
    gbasis(19)%expo(2,1)=52.402039790_DP
    gbasis(19)%expo(2,2)=12.177107100_DP
    gbasis(19)%expo(2,3)=3.960373165_DP
    gbasis(19)%coco(2,1)=-0.09996723_DP
    gbasis(19)%coco(2,2)=0.39951283_DP
    gbasis(19)%coco(2,3)=0.70011547_DP

    gbasis(19)%angmom(3)=1
    gbasis(19)%expo(3,1)=52.402039790_DP
    gbasis(19)%expo(3,2)=12.177107100_DP
    gbasis(19)%expo(3,3)=3.960373165_DP
    gbasis(19)%coco(3,1)=0.15591628_DP
    gbasis(19)%coco(3,2)=0.60768372_DP
    gbasis(19)%coco(3,3)=0.39195739_DP

    gbasis(19)%angmom(4)=0
    gbasis(19)%expo(4,1)=3.651583985_DP
    gbasis(19)%expo(4,2)=1.018782663_DP
    gbasis(19)%expo(4,3)=0.398744630_DP
    gbasis(19)%coco(4,1)=-0.21962037_DP
    gbasis(19)%coco(4,2)=0.22559543_DP
    gbasis(19)%coco(4,3)=0.90039843_DP

    gbasis(19)%angmom(5)=1
    gbasis(19)%expo(5,1)=3.651583985_DP
    gbasis(19)%expo(5,2)=1.018782663_DP
    gbasis(19)%expo(5,3)=0.398744630_DP
    gbasis(19)%coco(5,1)=0.01058760_DP
    gbasis(19)%coco(5,2)=0.59516701_DP
    gbasis(19)%coco(5,3)=0.46200101_DP

    gbasis(19)%angmom(6)=0
    gbasis(19)%expo(6,1)=0.503982251_DP
    gbasis(19)%expo(6,2)=0.186001147_DP
    gbasis(19)%expo(6,3)=0.082140067_DP
    gbasis(19)%coco(6,1)=-0.30884412_DP
    gbasis(19)%coco(6,2)=0.01960641_DP
    gbasis(19)%coco(6,3)=1.13103444_DP

    gbasis(19)%angmom(7)=1
    gbasis(19)%expo(7,1)=0.503982251_DP
    gbasis(19)%expo(7,2)=0.186001147_DP
    gbasis(19)%expo(7,3)=0.082140067_DP
    gbasis(19)%coco(7,1)=-0.12154686_DP
    gbasis(19)%coco(7,2)=0.57152276_DP
    gbasis(19)%coco(7,3)=0.54989495_DP

    gbasis(19)%angmom(8)=2
    gbasis(19)%expo(8,1)=0.20000000_DP
    gbasis(19)%coco(8,1)=1.00000000_DP
    ! ====end K ==========================================


    ! ====begin Ca ==========================================
    gbasis(20)%angmom(1)=0
    gbasis(20)%expo(1,1)=854.032495100_DP
    gbasis(20)%expo(1,2)=155.563085100_DP
    gbasis(20)%expo(1,3)=42.101441790_DP
    gbasis(20)%coco(1,1)=0.15432897_DP
    gbasis(20)%coco(1,2)=0.53532814_DP
    gbasis(20)%coco(1,3)=0.44463454_DP

    gbasis(20)%angmom(2)=0
    gbasis(20)%expo(2,1)=59.560299440_DP
    gbasis(20)%expo(2,2)=13.840532700_DP
    gbasis(20)%expo(2,3)=4.501370797_DP
    gbasis(20)%coco(2,1)=-0.09996723_DP
    gbasis(20)%coco(2,2)=0.39951283_DP
    gbasis(20)%coco(2,3)=0.70011547_DP

    gbasis(20)%angmom(3)=1
    gbasis(20)%expo(3,1)=59.560299440_DP
    gbasis(20)%expo(3,2)=13.840532700_DP
    gbasis(20)%expo(3,3)=4.501370797_DP
    gbasis(20)%coco(3,1)=0.15591628_DP
    gbasis(20)%coco(3,2)=0.60768372_DP
    gbasis(20)%coco(3,3)=0.39195739_DP

    gbasis(20)%angmom(4)=0
    gbasis(20)%expo(4,1)=4.374706256_DP
    gbasis(20)%expo(4,2)=1.220531941_DP
    gbasis(20)%expo(4,3)=0.477707930_DP
    gbasis(20)%coco(4,1)=-0.21962037_DP
    gbasis(20)%coco(4,2)=0.22559543_DP
    gbasis(20)%coco(4,3)=0.90039843_DP

    gbasis(20)%angmom(5)=1
    gbasis(20)%expo(5,1)=4.374706256_DP
    gbasis(20)%expo(5,2)=1.220531941_DP
    gbasis(20)%expo(5,3)=0.477707930_DP
    gbasis(20)%coco(5,1)=0.01058760_DP
    gbasis(20)%coco(5,2)=0.59516701_DP
    gbasis(20)%coco(5,3)=0.46200101_DP

    gbasis(20)%angmom(6)=0
    gbasis(20)%expo(6,1)=0.455848976_DP
    gbasis(20)%expo(6,2)=0.168236941_DP
    gbasis(20)%expo(6,3)=0.074295207_DP
    gbasis(20)%coco(6,1)=-0.30884412_DP
    gbasis(20)%coco(6,2)=0.01960641_DP
    gbasis(20)%coco(6,3)=1.13103444_DP

    gbasis(20)%angmom(7)=1
    gbasis(20)%expo(7,1)=0.455848976_DP
    gbasis(20)%expo(7,2)=0.168236941_DP
    gbasis(20)%expo(7,3)=0.074295207_DP
    gbasis(20)%coco(7,1)=-0.12154686_DP
    gbasis(20)%coco(7,2)=0.57152276_DP
    gbasis(20)%coco(7,3)=0.54989495_DP

    gbasis(20)%angmom(8)=2
    gbasis(20)%expo(8,1)=0.20000000_DP
    gbasis(20)%coco(8,1)=1.00000000_DP
    ! ====end Ca ==========================================

  end subroutine ngwf_data_cocos_and_expos01_20

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_data_cocos_and_expos21_40(gbasis)

    !===============================================================!
    ! Initialisation of Gaussian basis set exponents and contraction!
    ! coefficients. All the data of the STO-3G basis set are        !
    ! included for each element of the periodic table as well as    !
    ! the one "polarisation" function of the 6-31G* set.            !
    !---------------------------------------------------------------!
    ! This subroutine was created by Chris-Kriton Skylaris on       !
    ! 31/03/2006 but the actual STO-3G and 6-31G* parameters        !
    ! which appear below were typed by Shyong Chen, in March 2006.  !
    !===============================================================!

    implicit none

    type(GTO_SET), intent(inout) :: gbasis(1: pub_num_gtatoms)


    ! cks: all data in this subroutine was typed by Shyong Chen

    ! ====begin Sc ==========================================
    gbasis(21)%angmom(1)=0
    gbasis(21)%expo(1,1)=941.662425000_DP
    gbasis(21)%expo(1,2)=171.524986200_DP
    gbasis(21)%expo(1,3)=46.421355160_DP
    gbasis(21)%coco(1,1)=0.15432897_DP
    gbasis(21)%coco(1,2)=0.53532814_DP
    gbasis(21)%coco(1,3)=0.44463454_DP

    gbasis(21)%angmom(2)=0
    gbasis(21)%expo(2,1)=67.176687710_DP
    gbasis(21)%expo(2,2)=15.610417540_DP
    gbasis(21)%expo(2,3)=5.076992278_DP
    gbasis(21)%coco(2,1)=-0.09996723_DP
    gbasis(21)%coco(2,2)=0.39951283_DP
    gbasis(21)%coco(2,3)=0.70011547_DP

    gbasis(21)%angmom(3)=1
    gbasis(21)%expo(3,1)=67.176687710_DP
    gbasis(21)%expo(3,2)=15.610417540_DP
    gbasis(21)%expo(3,3)=5.076992278_DP
    gbasis(21)%coco(3,1)=0.15591628_DP
    gbasis(21)%coco(3,2)=0.60768372_DP
    gbasis(21)%coco(3,3)=0.39195739_DP

    gbasis(21)%angmom(4)=0
    gbasis(21)%expo(4,1)=4.698159231_DP
    gbasis(21)%expo(4,2)=1.433088313_DP
    gbasis(21)%expo(4,3)=0.552930024_DP
    gbasis(21)%coco(4,1)=-0.22776350_DP
    gbasis(21)%coco(4,2)=0.21754360_DP
    gbasis(21)%coco(4,3)=0.91667696_DP

    gbasis(21)%angmom(5)=1
    gbasis(21)%expo(5,1)=4.698159231_DP
    gbasis(21)%expo(5,2)=1.433088313_DP
    gbasis(21)%expo(5,3)=0.552930024_DP
    gbasis(21)%coco(5,1)=0.00495151_DP
    gbasis(21)%coco(5,2)=0.57776647_DP
    gbasis(21)%coco(5,3)=0.48464604_DP

    gbasis(21)%angmom(6)=2
    gbasis(21)%expo(6,1)=0.551700068_DP
    gbasis(21)%expo(6,2)=0.168286106_DP
    gbasis(21)%expo(6,3)=0.064930011_DP
    gbasis(21)%coco(6,1)=0.21976795_DP
    gbasis(21)%coco(6,2)=0.65554736_DP
    gbasis(21)%coco(6,3)=0.28657326_DP

    gbasis(21)%angmom(7)=0
    gbasis(21)%expo(7,1)=0.630932838_DP
    gbasis(21)%expo(7,2)=0.232853898_DP
    gbasis(21)%expo(7,3)=0.102830736_DP
    gbasis(21)%coco(7,1)=-0.30884412_DP
    gbasis(21)%coco(7,2)=0.01960641_DP
    gbasis(21)%coco(7,3)=1.13103444_DP

    gbasis(21)%angmom(8)=1
    gbasis(21)%expo(8,1)=0.630932838_DP
    gbasis(21)%expo(8,2)=0.232853898_DP
    gbasis(21)%expo(8,3)=0.102830736_DP
    gbasis(21)%coco(8,1)=-0.12154686_DP
    gbasis(21)%coco(8,2)=0.57152276_DP
    gbasis(21)%coco(8,3)=0.54989495_DP

    gbasis(21)%angmom(9)=3
    gbasis(21)%expo(9,1)=0.80000000_DP
    gbasis(21)%coco(9,1)=1.00000000_DP
    ! ====end Sc ==========================================


    ! ====begin Ti ==========================================
    gbasis(22)%angmom(1)=0
    gbasis(22)%expo(1,1)=1033.571245000_DP
    gbasis(22)%expo(1,2)=188.266292600_DP
    gbasis(22)%expo(1,3)=50.952206010_DP
    gbasis(22)%coco(1,1)=0.15432897_DP
    gbasis(22)%coco(1,2)=0.53532814_DP
    gbasis(22)%coco(1,3)=0.44463454_DP

    gbasis(22)%angmom(2)=0
    gbasis(22)%expo(2,1)=75.251204600_DP
    gbasis(22)%expo(2,2)=17.486761620_DP
    gbasis(22)%expo(2,3)=5.687237606_DP
    gbasis(22)%coco(2,1)=-0.09996723_DP
    gbasis(22)%coco(2,2)=0.39951283_DP
    gbasis(22)%coco(2,3)=0.70011547_DP

    gbasis(22)%angmom(3)=1
    gbasis(22)%expo(3,1)=75.251204600_DP
    gbasis(22)%expo(3,2)=17.486761620_DP
    gbasis(22)%expo(3,3)=5.687237606_DP
    gbasis(22)%coco(3,1)=0.15591628_DP
    gbasis(22)%coco(3,2)=0.60768372_DP
    gbasis(22)%coco(3,3)=0.39195739_DP

    gbasis(22)%angmom(4)=0
    gbasis(22)%expo(4,1)=5.395535474_DP
    gbasis(22)%expo(4,2)=1.645810296_DP
    gbasis(22)%expo(4,3)=0.635004777_DP
    gbasis(22)%coco(4,1)=-0.22776350_DP
    gbasis(22)%coco(4,2)=0.21754360_DP
    gbasis(22)%coco(4,3)=0.91667696_DP

    gbasis(22)%angmom(5)=1
    gbasis(22)%expo(5,1)=5.395535474_DP
    gbasis(22)%expo(5,2)=1.645810296_DP
    gbasis(22)%expo(5,3)=0.635004777_DP
    gbasis(22)%coco(5,1)=0.00495151_DP
    gbasis(22)%coco(5,2)=0.57776647_DP
    gbasis(22)%coco(5,3)=0.48464604_DP

    gbasis(22)%angmom(6)=2
    gbasis(22)%expo(6,1)=1.645981194_DP
    gbasis(22)%expo(6,2)=0.502076728_DP
    gbasis(22)%expo(6,3)=0.193716810_DP
    gbasis(22)%coco(6,1)=0.21976795_DP
    gbasis(22)%coco(6,2)=0.65554736_DP
    gbasis(22)%coco(6,3)=0.28657326_DP

    gbasis(22)%angmom(7)=0
    gbasis(22)%expo(7,1)=0.712264025_DP
    gbasis(22)%expo(7,2)=0.262870220_DP
    gbasis(22)%expo(7,3)=0.116086261_DP
    gbasis(22)%coco(7,1)=-0.30884412_DP
    gbasis(22)%coco(7,2)=0.01960641_DP
    gbasis(22)%coco(7,3)=1.13103444_DP

    gbasis(22)%angmom(8)=1
    gbasis(22)%expo(8,1)=0.712264025_DP
    gbasis(22)%expo(8,2)=0.262870220_DP
    gbasis(22)%expo(8,3)=0.116086261_DP
    gbasis(22)%coco(8,1)=-0.12154686_DP
    gbasis(22)%coco(8,2)=0.57152276_DP
    gbasis(22)%coco(8,3)=0.54989495_DP

    gbasis(22)%angmom(9)=3
    gbasis(22)%expo(9,1)=0.80000000_DP
    gbasis(22)%coco(9,1)=1.00000000_DP
    ! ====end Ti ==========================================


    ! ====begin V ==========================================
    gbasis(23)%angmom(1)=0
    gbasis(23)%expo(1,1)=1130.762517000_DP
    gbasis(23)%expo(1,2)=205.969804100_DP
    gbasis(23)%expo(1,3)=55.743467110_DP
    gbasis(23)%coco(1,1)=0.15432897_DP
    gbasis(23)%coco(1,2)=0.53532814_DP
    gbasis(23)%coco(1,3)=0.44463454_DP

    gbasis(23)%angmom(2)=0
    gbasis(23)%expo(2,1)=83.783850110_DP
    gbasis(23)%expo(2,2)=19.469564930_DP
    gbasis(23)%expo(2,3)=6.332106784_DP
    gbasis(23)%coco(2,1)=-0.09996723_DP
    gbasis(23)%coco(2,2)=0.39951283_DP
    gbasis(23)%coco(2,3)=0.70011547_DP

    gbasis(23)%angmom(3)=1
    gbasis(23)%expo(3,1)=83.783850110_DP
    gbasis(23)%expo(3,2)=19.469564930_DP
    gbasis(23)%expo(3,3)=6.332106784_DP
    gbasis(23)%coco(3,1)=0.15591628_DP
    gbasis(23)%coco(3,2)=0.60768372_DP
    gbasis(23)%coco(3,3)=0.39195739_DP

    gbasis(23)%angmom(4)=0
    gbasis(23)%expo(4,1)=6.141151276_DP
    gbasis(23)%expo(4,2)=1.873246881_DP
    gbasis(23)%expo(4,3)=0.722756883_DP
    gbasis(23)%coco(4,1)=-0.22776350_DP
    gbasis(23)%coco(4,2)=0.21754360_DP
    gbasis(23)%coco(4,3)=0.91667696_DP

    gbasis(23)%angmom(5)=1
    gbasis(23)%expo(5,1)=6.141151276_DP
    gbasis(23)%expo(5,2)=1.873246881_DP
    gbasis(23)%expo(5,3)=0.722756883_DP
    gbasis(23)%coco(5,1)=0.00495151_DP
    gbasis(23)%coco(5,2)=0.57776647_DP
    gbasis(23)%coco(5,3)=0.48464604_DP

    gbasis(23)%angmom(6)=2
    gbasis(23)%expo(6,1)=2.964817927_DP
    gbasis(23)%expo(6,2)=0.904363968_DP
    gbasis(23)%expo(6,3)=0.348931734_DP
    gbasis(23)%coco(6,1)=0.21976795_DP
    gbasis(23)%coco(6,2)=0.65554736_DP
    gbasis(23)%coco(6,3)=0.28657326_DP

    gbasis(23)%angmom(7)=0
    gbasis(23)%expo(7,1)=0.712264025_DP
    gbasis(23)%expo(7,2)=0.262870220_DP
    gbasis(23)%expo(7,3)=0.116086261_DP
    gbasis(23)%coco(7,1)=-0.30884412_DP
    gbasis(23)%coco(7,2)=0.01960641_DP
    gbasis(23)%coco(7,3)=1.13103444_DP

    gbasis(23)%angmom(8)=1
    gbasis(23)%expo(8,1)=0.712264025_DP
    gbasis(23)%expo(8,2)=0.262870220_DP
    gbasis(23)%expo(8,3)=0.116086261_DP
    gbasis(23)%coco(8,1)=-0.12154686_DP
    gbasis(23)%coco(8,2)=0.57152276_DP
    gbasis(23)%coco(8,3)=0.54989495_DP

    gbasis(23)%angmom(9)=3
    gbasis(23)%expo(9,1)=0.80000000_DP
    gbasis(23)%coco(9,1)=1.00000000_DP
    ! ====end V ==========================================


    ! ====begin Cr ==========================================
    gbasis(24)%angmom(1)=0
    gbasis(24)%expo(1,1)=1232.320450000_DP
    gbasis(24)%expo(1,2)=224.468708200_DP
    gbasis(24)%expo(1,3)=60.749992510_DP
    gbasis(24)%coco(1,1)=0.15432897_DP
    gbasis(24)%coco(1,2)=0.53532814_DP
    gbasis(24)%coco(1,3)=0.44463454_DP

    gbasis(24)%angmom(2)=0
    gbasis(24)%expo(2,1)=92.774624230_DP
    gbasis(24)%expo(2,2)=21.558827490_DP
    gbasis(24)%expo(2,3)=7.011599810_DP
    gbasis(24)%coco(2,1)=-0.09996723_DP
    gbasis(24)%coco(2,2)=0.39951283_DP
    gbasis(24)%coco(2,3)=0.70011547_DP

    gbasis(24)%angmom(3)=1
    gbasis(24)%expo(3,1)=92.774624230_DP
    gbasis(24)%expo(3,2)=21.558827490_DP
    gbasis(24)%expo(3,3)=7.011599810_DP
    gbasis(24)%coco(3,1)=0.15591628_DP
    gbasis(24)%coco(3,2)=0.60768372_DP
    gbasis(24)%coco(3,3)=0.39195739_DP

    gbasis(24)%angmom(4)=0
    gbasis(24)%expo(4,1)=6.899488096_DP
    gbasis(24)%expo(4,2)=2.104563782_DP
    gbasis(24)%expo(4,3)=0.812006134_DP
    gbasis(24)%coco(4,1)=-0.22776350_DP
    gbasis(24)%coco(4,2)=0.21754360_DP
    gbasis(24)%coco(4,3)=0.91667696_DP

    gbasis(24)%angmom(5)=1
    gbasis(24)%expo(5,1)=6.899488096_DP
    gbasis(24)%expo(5,2)=2.104563782_DP
    gbasis(24)%expo(5,3)=0.812006134_DP
    gbasis(24)%coco(5,1)=0.00495151_DP
    gbasis(24)%coco(5,2)=0.57776647_DP
    gbasis(24)%coco(5,3)=0.48464604_DP

    gbasis(24)%angmom(6)=2
    gbasis(24)%expo(6,1)=4.241479241_DP
    gbasis(24)%expo(6,2)=1.293786360_DP
    gbasis(24)%expo(6,3)=0.499182999_DP
    gbasis(24)%coco(6,1)=0.21976795_DP
    gbasis(24)%coco(6,2)=0.65554736_DP
    gbasis(24)%coco(6,3)=0.28657326_DP

    gbasis(24)%angmom(7)=0
    gbasis(24)%expo(7,1)=0.754778054_DP
    gbasis(24)%expo(7,2)=0.278560571_DP
    gbasis(24)%expo(7,3)=0.123015285_DP
    gbasis(24)%coco(7,1)=-0.30884412_DP
    gbasis(24)%coco(7,2)=0.01960641_DP
    gbasis(24)%coco(7,3)=1.13103444_DP

    gbasis(24)%angmom(8)=1
    gbasis(24)%expo(8,1)=0.754778054_DP
    gbasis(24)%expo(8,2)=0.278560571_DP
    gbasis(24)%expo(8,3)=0.123015285_DP
    gbasis(24)%coco(8,1)=-0.12154686_DP
    gbasis(24)%coco(8,2)=0.57152276_DP
    gbasis(24)%coco(8,3)=0.54989495_DP

    gbasis(24)%angmom(9)=3
    gbasis(24)%expo(9,1)=0.80000000_DP
    gbasis(24)%coco(9,1)=1.00000000_DP
    ! ====end Cr ==========================================


    ! ====begin Mn ==========================================
    gbasis(25)%angmom(1)=0
    gbasis(25)%expo(1,1)=1337.153266000_DP
    gbasis(25)%expo(1,2)=243.564136500_DP
    gbasis(25)%expo(1,3)=65.917960620_DP
    gbasis(25)%coco(1,1)=0.15432897_DP
    gbasis(25)%coco(1,2)=0.53532814_DP
    gbasis(25)%coco(1,3)=0.44463454_DP

    gbasis(25)%angmom(2)=0
    gbasis(25)%expo(2,1)=102.022002100_DP
    gbasis(25)%expo(2,2)=23.707719230_DP
    gbasis(25)%expo(2,3)=7.710486098_DP
    gbasis(25)%coco(2,1)=-0.09996723_DP
    gbasis(25)%coco(2,2)=0.39951283_DP
    gbasis(25)%coco(2,3)=0.70011547_DP

    gbasis(25)%angmom(3)=1
    gbasis(25)%expo(3,1)=102.022002100_DP
    gbasis(25)%expo(3,2)=23.707719230_DP
    gbasis(25)%expo(3,3)=7.710486098_DP
    gbasis(25)%coco(3,1)=0.15591628_DP
    gbasis(25)%coco(3,2)=0.60768372_DP
    gbasis(25)%coco(3,3)=0.39195739_DP

    gbasis(25)%angmom(4)=0
    gbasis(25)%expo(4,1)=7.701960922_DP
    gbasis(25)%expo(4,2)=2.349343572_DP
    gbasis(25)%expo(4,3)=0.906449787_DP
    gbasis(25)%coco(4,1)=-0.22776350_DP
    gbasis(25)%coco(4,2)=0.21754360_DP
    gbasis(25)%coco(4,3)=0.91667696_DP

    gbasis(25)%angmom(5)=1
    gbasis(25)%expo(5,1)=7.701960922_DP
    gbasis(25)%expo(5,2)=2.349343572_DP
    gbasis(25)%expo(5,3)=0.906449787_DP
    gbasis(25)%coco(5,1)=0.00495151_DP
    gbasis(25)%coco(5,2)=0.57776647_DP
    gbasis(25)%coco(5,3)=0.48464604_DP

    gbasis(25)%angmom(6)=2
    gbasis(25)%expo(6,1)=5.426950461_DP
    gbasis(25)%expo(6,2)=1.655392868_DP
    gbasis(25)%expo(6,3)=0.638702032_DP
    gbasis(25)%coco(6,1)=0.21976795_DP
    gbasis(25)%coco(6,2)=0.65554736_DP
    gbasis(25)%coco(6,3)=0.28657326_DP

    gbasis(25)%angmom(7)=0
    gbasis(25)%expo(7,1)=0.670982286_DP
    gbasis(25)%expo(7,2)=0.247634663_DP
    gbasis(25)%expo(7,3)=0.109358078_DP
    gbasis(25)%coco(7,1)=-0.30884412_DP
    gbasis(25)%coco(7,2)=0.01960641_DP
    gbasis(25)%coco(7,3)=1.13103444_DP

    gbasis(25)%angmom(8)=1
    gbasis(25)%expo(8,1)=0.670982286_DP
    gbasis(25)%expo(8,2)=0.247634663_DP
    gbasis(25)%expo(8,3)=0.109358078_DP
    gbasis(25)%coco(8,1)=-0.12154686_DP
    gbasis(25)%coco(8,2)=0.57152276_DP
    gbasis(25)%coco(8,3)=0.54989495_DP

    gbasis(25)%angmom(9)=3
    gbasis(25)%expo(9,1)=0.80000000_DP
    gbasis(25)%coco(9,1)=1.00000000_DP
    ! ====end Mn ==========================================


    ! ====begin Fe ==========================================
    gbasis(26)%angmom(1)=0
    gbasis(26)%expo(1,1)=1447.400411000_DP
    gbasis(26)%expo(1,2)=263.645791600_DP
    gbasis(26)%expo(1,3)=71.352840190_DP
    gbasis(26)%coco(1,1)=0.15432897_DP
    gbasis(26)%coco(1,2)=0.53532814_DP
    gbasis(26)%coco(1,3)=0.44463454_DP

    gbasis(26)%angmom(2)=0
    gbasis(26)%expo(2,1)=111.919489100_DP
    gbasis(26)%expo(2,2)=26.007682360_DP
    gbasis(26)%expo(2,3)=8.458505490_DP
    gbasis(26)%coco(2,1)=-0.09996723_DP
    gbasis(26)%coco(2,2)=0.39951283_DP
    gbasis(26)%coco(2,3)=0.70011547_DP

    gbasis(26)%angmom(3)=1
    gbasis(26)%expo(3,1)=111.919489100_DP
    gbasis(26)%expo(3,2)=26.007682360_DP
    gbasis(26)%expo(3,3)=8.458505490_DP
    gbasis(26)%coco(3,1)=0.15591628_DP
    gbasis(26)%coco(3,2)=0.60768372_DP
    gbasis(26)%coco(3,3)=0.39195739_DP

    gbasis(26)%angmom(4)=0
    gbasis(26)%expo(4,1)=8.548569754_DP
    gbasis(26)%expo(4,2)=2.607586250_DP
    gbasis(26)%expo(4,3)=1.006087840_DP
    gbasis(26)%coco(4,1)=-0.22776350_DP
    gbasis(26)%coco(4,2)=0.21754360_DP
    gbasis(26)%coco(4,3)=0.91667696_DP

    gbasis(26)%angmom(5)=1
    gbasis(26)%expo(5,1)=8.548569754_DP
    gbasis(26)%expo(5,2)=2.607586250_DP
    gbasis(26)%expo(5,3)=1.006087840_DP
    gbasis(26)%coco(5,1)=0.00495151_DP
    gbasis(26)%coco(5,2)=0.57776647_DP
    gbasis(26)%coco(5,3)=0.48464604_DP

    gbasis(26)%angmom(6)=2
    gbasis(26)%expo(6,1)=6.411803475_DP
    gbasis(26)%expo(6,2)=1.955804428_DP
    gbasis(26)%expo(6,3)=0.754610151_DP
    gbasis(26)%coco(6,1)=0.21976795_DP
    gbasis(26)%coco(6,2)=0.65554736_DP
    gbasis(26)%coco(6,3)=0.28657326_DP

    gbasis(26)%angmom(7)=0
    gbasis(26)%expo(7,1)=0.592115681_DP
    gbasis(26)%expo(7,2)=0.218527925_DP
    gbasis(26)%expo(7,3)=0.096504236_DP
    gbasis(26)%coco(7,1)=-0.30884412_DP
    gbasis(26)%coco(7,2)=0.01960641_DP
    gbasis(26)%coco(7,3)=1.13103444_DP

    gbasis(26)%angmom(8)=1
    gbasis(26)%expo(8,1)=0.592115681_DP
    gbasis(26)%expo(8,2)=0.218527925_DP
    gbasis(26)%expo(8,3)=0.096504236_DP
    gbasis(26)%coco(8,1)=-0.12154686_DP
    gbasis(26)%coco(8,2)=0.57152276_DP
    gbasis(26)%coco(8,3)=0.54989495_DP

    gbasis(26)%angmom(9)=3
    gbasis(26)%expo(9,1)=0.80000000_DP
    gbasis(26)%coco(9,1)=1.00000000_DP
    ! ====end Fe ==========================================


    ! ====begin Co ==========================================
    gbasis(27)%angmom(1)=0
    gbasis(27)%expo(1,1)=1557.298704000_DP
    gbasis(27)%expo(1,2)=283.663902900_DP
    gbasis(27)%expo(1,3)=76.770522340_DP
    gbasis(27)%coco(1,1)=0.15432897_DP
    gbasis(27)%coco(1,2)=0.53532814_DP
    gbasis(27)%coco(1,3)=0.44463454_DP

    gbasis(27)%angmom(2)=0
    gbasis(27)%expo(2,1)=121.834474100_DP
    gbasis(27)%expo(2,2)=28.311711640_DP
    gbasis(27)%expo(2,3)=9.207847321_DP
    gbasis(27)%coco(2,1)=-0.09996723_DP
    gbasis(27)%coco(2,2)=0.39951283_DP
    gbasis(27)%coco(2,3)=0.70011547_DP

    gbasis(27)%angmom(3)=1
    gbasis(27)%expo(3,1)=121.834474100_DP
    gbasis(27)%expo(3,2)=28.311711640_DP
    gbasis(27)%expo(3,3)=9.207847321_DP
    gbasis(27)%coco(3,1)=0.15591628_DP
    gbasis(27)%coco(3,2)=0.60768372_DP
    gbasis(27)%coco(3,3)=0.39195739_DP

    gbasis(27)%angmom(4)=0
    gbasis(27)%expo(4,1)=9.480851678_DP
    gbasis(27)%expo(4,2)=2.891961952_DP
    gbasis(27)%expo(4,3)=1.115808827_DP
    gbasis(27)%coco(4,1)=-0.22776350_DP
    gbasis(27)%coco(4,2)=0.21754360_DP
    gbasis(27)%coco(4,3)=0.91667696_DP

    gbasis(27)%angmom(5)=1
    gbasis(27)%expo(5,1)=9.480851678_DP
    gbasis(27)%expo(5,2)=2.891961952_DP
    gbasis(27)%expo(5,3)=1.115808827_DP
    gbasis(27)%coco(5,1)=0.00495151_DP
    gbasis(27)%coco(5,2)=0.57776647_DP
    gbasis(27)%coco(5,3)=0.48464604_DP

    gbasis(27)%angmom(6)=2
    gbasis(27)%expo(6,1)=7.664527389_DP
    gbasis(27)%expo(6,2)=2.337925151_DP
    gbasis(27)%expo(6,3)=0.902044205_DP
    gbasis(27)%coco(6,1)=0.21976795_DP
    gbasis(27)%coco(6,2)=0.65554736_DP
    gbasis(27)%coco(6,3)=0.28657326_DP

    gbasis(27)%angmom(7)=0
    gbasis(27)%expo(7,1)=0.592115681_DP
    gbasis(27)%expo(7,2)=0.218527925_DP
    gbasis(27)%expo(7,3)=0.096504236_DP
    gbasis(27)%coco(7,1)=-0.30884412_DP
    gbasis(27)%coco(7,2)=0.01960641_DP
    gbasis(27)%coco(7,3)=1.13103444_DP

    gbasis(27)%angmom(8)=1
    gbasis(27)%expo(8,1)=0.592115681_DP
    gbasis(27)%expo(8,2)=0.218527925_DP
    gbasis(27)%expo(8,3)=0.096504236_DP
    gbasis(27)%coco(8,1)=-0.12154686_DP
    gbasis(27)%coco(8,2)=0.57152276_DP
    gbasis(27)%coco(8,3)=0.54989495_DP

    gbasis(27)%angmom(9)=3
    gbasis(27)%expo(9,1)=0.80000000_DP
    gbasis(27)%coco(9,1)=1.00000000_DP
    ! ====end Co ==========================================


    ! ====begin Ni ==========================================
    gbasis(28)%angmom(1)=0
    gbasis(28)%expo(1,1)=1679.771028000_DP
    gbasis(28)%expo(1,2)=305.972389600_DP
    gbasis(28)%expo(1,3)=82.808069430_DP
    gbasis(28)%coco(1,1)=0.15432897_DP
    gbasis(28)%coco(1,2)=0.53532814_DP
    gbasis(28)%coco(1,3)=0.44463454_DP

    gbasis(28)%angmom(2)=0
    gbasis(28)%expo(2,1)=132.858889900_DP
    gbasis(28)%expo(2,2)=30.873548780_DP
    gbasis(28)%expo(2,3)=10.041036270_DP
    gbasis(28)%coco(2,1)=-0.09996723_DP
    gbasis(28)%coco(2,2)=0.39951283_DP
    gbasis(28)%coco(2,3)=0.70011547_DP

    gbasis(28)%angmom(3)=1
    gbasis(28)%expo(3,1)=132.858889900_DP
    gbasis(28)%expo(3,2)=30.873548780_DP
    gbasis(28)%expo(3,3)=10.041036270_DP
    gbasis(28)%coco(3,1)=0.15591628_DP
    gbasis(28)%coco(3,2)=0.60768372_DP
    gbasis(28)%coco(3,3)=0.39195739_DP

    gbasis(28)%angmom(4)=0
    gbasis(28)%expo(4,1)=10.330743350_DP
    gbasis(28)%expo(4,2)=3.151206003_DP
    gbasis(28)%expo(4,3)=1.215833241_DP
    gbasis(28)%coco(4,1)=-0.22776350_DP
    gbasis(28)%coco(4,2)=0.21754360_DP
    gbasis(28)%coco(4,3)=0.91667696_DP

    gbasis(28)%angmom(5)=1
    gbasis(28)%expo(5,1)=10.330743350_DP
    gbasis(28)%expo(5,2)=3.151206003_DP
    gbasis(28)%expo(5,3)=1.215833241_DP
    gbasis(28)%coco(5,1)=0.00495151_DP
    gbasis(28)%coco(5,2)=0.57776647_DP
    gbasis(28)%coco(5,3)=0.48464604_DP

    gbasis(28)%angmom(6)=2
    gbasis(28)%expo(6,1)=8.627722755_DP
    gbasis(28)%expo(6,2)=2.631730438_DP
    gbasis(28)%expo(6,3)=1.015403419_DP
    gbasis(28)%coco(6,1)=0.21976795_DP
    gbasis(28)%coco(6,2)=0.65554736_DP
    gbasis(28)%coco(6,3)=0.28657326_DP

    gbasis(28)%angmom(7)=0
    gbasis(28)%expo(7,1)=0.630932838_DP
    gbasis(28)%expo(7,2)=0.232853898_DP
    gbasis(28)%expo(7,3)=0.102830736_DP
    gbasis(28)%coco(7,1)=-0.30884412_DP
    gbasis(28)%coco(7,2)=0.01960641_DP
    gbasis(28)%coco(7,3)=1.13103444_DP

    gbasis(28)%angmom(8)=1
    gbasis(28)%expo(8,1)=0.630932838_DP
    gbasis(28)%expo(8,2)=0.232853898_DP
    gbasis(28)%expo(8,3)=0.102830736_DP
    gbasis(28)%coco(8,1)=-0.12154686_DP
    gbasis(28)%coco(8,2)=0.57152276_DP
    gbasis(28)%coco(8,3)=0.54989495_DP

    gbasis(28)%angmom(9)=3
    gbasis(28)%expo(9,1)=0.80000000_DP
    gbasis(28)%coco(9,1)=1.00000000_DP
    ! ====end Ni ==========================================


    ! ====begin Cu ==========================================
    gbasis(29)%angmom(1)=0
    gbasis(29)%expo(1,1)=1801.806730000_DP
    gbasis(29)%expo(1,2)=328.201345000_DP
    gbasis(29)%expo(1,3)=88.824092280_DP
    gbasis(29)%coco(1,1)=0.15432897_DP
    gbasis(29)%coco(1,2)=0.53532814_DP
    gbasis(29)%coco(1,3)=0.44463454_DP

    gbasis(29)%angmom(2)=0
    gbasis(29)%expo(2,1)=144.121218400_DP
    gbasis(29)%expo(2,2)=33.490671730_DP
    gbasis(29)%expo(2,3)=10.892205880_DP
    gbasis(29)%coco(2,1)=-0.09996723_DP
    gbasis(29)%coco(2,2)=0.39951283_DP
    gbasis(29)%coco(2,3)=0.70011547_DP

    gbasis(29)%angmom(3)=1
    gbasis(29)%expo(3,1)=144.121218400_DP
    gbasis(29)%expo(3,2)=33.490671730_DP
    gbasis(29)%expo(3,3)=10.892205880_DP
    gbasis(29)%coco(3,1)=0.15591628_DP
    gbasis(29)%coco(3,2)=0.60768372_DP
    gbasis(29)%coco(3,3)=0.39195739_DP

    gbasis(29)%angmom(4)=0
    gbasis(29)%expo(4,1)=11.307754020_DP
    gbasis(29)%expo(4,2)=3.449225397_DP
    gbasis(29)%expo(4,3)=1.330818388_DP
    gbasis(29)%coco(4,1)=-0.22776350_DP
    gbasis(29)%coco(4,2)=0.21754360_DP
    gbasis(29)%coco(4,3)=0.91667696_DP

    gbasis(29)%angmom(5)=1
    gbasis(29)%expo(5,1)=11.307754020_DP
    gbasis(29)%expo(5,2)=3.449225397_DP
    gbasis(29)%expo(5,3)=1.330818388_DP
    gbasis(29)%coco(5,1)=0.00495151_DP
    gbasis(29)%coco(5,2)=0.57776647_DP
    gbasis(29)%coco(5,3)=0.48464604_DP

    gbasis(29)%angmom(6)=2
    gbasis(29)%expo(6,1)=9.647911930_DP
    gbasis(29)%expo(6,2)=2.942920654_DP
    gbasis(29)%expo(6,3)=1.135470278_DP
    gbasis(29)%coco(6,1)=0.21976795_DP
    gbasis(29)%coco(6,2)=0.65554736_DP
    gbasis(29)%coco(6,3)=0.28657326_DP

    gbasis(29)%angmom(7)=0
    gbasis(29)%expo(7,1)=0.630932838_DP
    gbasis(29)%expo(7,2)=0.232853898_DP
    gbasis(29)%expo(7,3)=0.102830736_DP
    gbasis(29)%coco(7,1)=-0.30884412_DP
    gbasis(29)%coco(7,2)=0.01960641_DP
    gbasis(29)%coco(7,3)=1.13103444_DP

    gbasis(29)%angmom(8)=1
    gbasis(29)%expo(8,1)=0.630932838_DP
    gbasis(29)%expo(8,2)=0.232853898_DP
    gbasis(29)%expo(8,3)=0.102830736_DP
    gbasis(29)%coco(8,1)=-0.12154686_DP
    gbasis(29)%coco(8,2)=0.57152276_DP
    gbasis(29)%coco(8,3)=0.54989495_DP

    gbasis(29)%angmom(9)=3
    gbasis(29)%expo(9,1)=0.80000000_DP
    gbasis(29)%coco(9,1)=1.00000000_DP
    ! ====end Cu ==========================================


    ! ====begin Zn ==========================================
    gbasis(30)%angmom(1)=0
    gbasis(30)%expo(1,1)=1929.432301000_DP
    gbasis(30)%expo(1,2)=351.448502100_DP
    gbasis(30)%expo(1,3)=95.115680210_DP
    gbasis(30)%coco(1,1)=0.15432897_DP
    gbasis(30)%coco(1,2)=0.53532814_DP
    gbasis(30)%coco(1,3)=0.44463454_DP

    gbasis(30)%angmom(2)=0
    gbasis(30)%expo(2,1)=155.841675500_DP
    gbasis(30)%expo(2,2)=36.214253910_DP
    gbasis(30)%expo(2,3)=11.777999340_DP
    gbasis(30)%coco(2,1)=-0.09996723_DP
    gbasis(30)%coco(2,2)=0.39951283_DP
    gbasis(30)%coco(2,3)=0.70011547_DP

    gbasis(30)%angmom(3)=1
    gbasis(30)%expo(3,1)=155.841675500_DP
    gbasis(30)%expo(3,2)=36.214253910_DP
    gbasis(30)%expo(3,3)=11.777999340_DP
    gbasis(30)%coco(3,1)=0.15591628_DP
    gbasis(30)%coco(3,2)=0.60768372_DP
    gbasis(30)%coco(3,3)=0.39195739_DP

    gbasis(30)%angmom(4)=0
    gbasis(30)%expo(4,1)=12.281527440_DP
    gbasis(30)%expo(4,2)=3.746257327_DP
    gbasis(30)%expo(4,3)=1.445422541_DP
    gbasis(30)%coco(4,1)=-0.22776350_DP
    gbasis(30)%coco(4,2)=0.21754360_DP
    gbasis(30)%coco(4,3)=0.91667696_DP

    gbasis(30)%angmom(5)=1
    gbasis(30)%expo(5,1)=12.281527440_DP
    gbasis(30)%expo(5,2)=3.746257327_DP
    gbasis(30)%expo(5,3)=1.445422541_DP
    gbasis(30)%coco(5,1)=0.00495151_DP
    gbasis(30)%coco(5,2)=0.57776647_DP
    gbasis(30)%coco(5,3)=0.48464604_DP

    gbasis(30)%angmom(6)=2
    gbasis(30)%expo(6,1)=10.947370770_DP
    gbasis(30)%expo(6,2)=3.339297018_DP
    gbasis(30)%expo(6,3)=1.288404602_DP
    gbasis(30)%coco(6,1)=0.21976795_DP
    gbasis(30)%coco(6,2)=0.65554736_DP
    gbasis(30)%coco(6,3)=0.28657326_DP

    gbasis(30)%angmom(7)=0
    gbasis(30)%expo(7,1)=0.889713885_DP
    gbasis(30)%expo(7,2)=0.328360379_DP
    gbasis(30)%expo(7,3)=0.145007406_DP
    gbasis(30)%coco(7,1)=-0.30884412_DP
    gbasis(30)%coco(7,2)=0.01960641_DP
    gbasis(30)%coco(7,3)=1.13103444_DP

    gbasis(30)%angmom(8)=1
    gbasis(30)%expo(8,1)=0.889713885_DP
    gbasis(30)%expo(8,2)=0.328360379_DP
    gbasis(30)%expo(8,3)=0.145007406_DP
    gbasis(30)%coco(8,1)=-0.12154686_DP
    gbasis(30)%coco(8,2)=0.57152276_DP
    gbasis(30)%coco(8,3)=0.54989495_DP

    gbasis(30)%angmom(9)=3
    gbasis(30)%expo(9,1)=0.80000000_DP
    gbasis(30)%coco(9,1)=1.00000000_DP
    ! ====end Zn ==========================================


    ! ====begin Ga ==========================================
    gbasis(31)%angmom(1)=0
    gbasis(31)%expo(1,1)=2061.424532000_DP
    gbasis(31)%expo(1,2)=375.491051700_DP
    gbasis(31)%expo(1,3)=101.622532400_DP
    gbasis(31)%coco(1,1)=0.15432897_DP
    gbasis(31)%coco(1,2)=0.53532814_DP
    gbasis(31)%coco(1,3)=0.44463454_DP

    gbasis(31)%angmom(2)=0
    gbasis(31)%expo(2,1)=167.761868000_DP
    gbasis(31)%expo(2,2)=38.984250280_DP
    gbasis(31)%expo(2,3)=12.678888130_DP
    gbasis(31)%coco(2,1)=-0.09996723_DP
    gbasis(31)%coco(2,2)=0.39951283_DP
    gbasis(31)%coco(2,3)=0.70011547_DP

    gbasis(31)%angmom(3)=1
    gbasis(31)%expo(3,1)=167.761868000_DP
    gbasis(31)%expo(3,2)=38.984250280_DP
    gbasis(31)%expo(3,3)=12.678888130_DP
    gbasis(31)%coco(3,1)=0.15591628_DP
    gbasis(31)%coco(3,2)=0.60768372_DP
    gbasis(31)%coco(3,3)=0.39195739_DP

    gbasis(31)%angmom(4)=0
    gbasis(31)%expo(4,1)=12.615055200_DP
    gbasis(31)%expo(4,2)=3.847993927_DP
    gbasis(31)%expo(4,3)=1.484675684_DP
    gbasis(31)%coco(4,1)=-0.22776350_DP
    gbasis(31)%coco(4,2)=0.21754360_DP
    gbasis(31)%coco(4,3)=0.91667696_DP

    gbasis(31)%angmom(5)=1
    gbasis(31)%expo(5,1)=12.615055200_DP
    gbasis(31)%expo(5,2)=3.847993927_DP
    gbasis(31)%expo(5,3)=1.484675684_DP
    gbasis(31)%coco(5,1)=0.00495151_DP
    gbasis(31)%coco(5,2)=0.57776647_DP
    gbasis(31)%coco(5,3)=0.48464604_DP

    gbasis(31)%angmom(6)=2
    gbasis(31)%expo(6,1)=12.615055200_DP
    gbasis(31)%expo(6,2)=3.847993927_DP
    gbasis(31)%expo(6,3)=1.484675684_DP
    gbasis(31)%coco(6,1)=0.21976795_DP
    gbasis(31)%coco(6,2)=0.65554736_DP
    gbasis(31)%coco(6,3)=0.28657326_DP

    gbasis(31)%angmom(7)=0
    gbasis(31)%expo(7,1)=0.798524374_DP
    gbasis(31)%expo(7,2)=0.294705714_DP
    gbasis(31)%expo(7,3)=0.130145151_DP
    gbasis(31)%coco(7,1)=-0.30884412_DP
    gbasis(31)%coco(7,2)=0.01960641_DP
    gbasis(31)%coco(7,3)=1.13103444_DP

    gbasis(31)%angmom(8)=1
    gbasis(31)%expo(8,1)=0.798524374_DP
    gbasis(31)%expo(8,2)=0.294705714_DP
    gbasis(31)%expo(8,3)=0.130145151_DP
    gbasis(31)%coco(8,1)=-0.12154686_DP
    gbasis(31)%coco(8,2)=0.57152276_DP
    gbasis(31)%coco(8,3)=0.54989495_DP

    gbasis(31)%angmom(9)=2
    gbasis(31)%expo(9,1)=0.20700000000_DP
    gbasis(31)%coco(9,1)=1.00000000_DP
    ! ====end Ga ==========================================


    ! ====begin Ge ==========================================
    gbasis(32)%angmom(1)=0
    gbasis(32)%expo(1,1)=2196.384229000_DP
    gbasis(32)%expo(1,2)=400.074129200_DP
    gbasis(32)%expo(1,3)=108.275672600_DP
    gbasis(32)%coco(1,1)=0.15432897_DP
    gbasis(32)%coco(1,2)=0.53532814_DP
    gbasis(32)%coco(1,3)=0.44463454_DP

    gbasis(32)%angmom(2)=0
    gbasis(32)%expo(2,1)=180.389038000_DP
    gbasis(32)%expo(2,2)=41.918533040_DP
    gbasis(32)%expo(2,3)=13.633207950_DP
    gbasis(32)%coco(2,1)=-0.09996723_DP
    gbasis(32)%coco(2,2)=0.39951283_DP
    gbasis(32)%coco(2,3)=0.70011547_DP

    gbasis(32)%angmom(3)=1
    gbasis(32)%expo(3,1)=180.389038000_DP
    gbasis(32)%expo(3,2)=41.918533040_DP
    gbasis(32)%expo(3,3)=13.633207950_DP
    gbasis(32)%coco(3,1)=0.15591628_DP
    gbasis(32)%coco(3,2)=0.60768372_DP
    gbasis(32)%coco(3,3)=0.39195739_DP

    gbasis(32)%angmom(4)=0
    gbasis(32)%expo(4,1)=14.196656190_DP
    gbasis(32)%expo(4,2)=4.330432640_DP
    gbasis(32)%expo(4,3)=1.670815538_DP
    gbasis(32)%coco(4,1)=-0.22776350_DP
    gbasis(32)%coco(4,2)=0.21754360_DP
    gbasis(32)%coco(4,3)=0.91667696_DP

    gbasis(32)%angmom(5)=1
    gbasis(32)%expo(5,1)=14.196656190_DP
    gbasis(32)%expo(5,2)=4.330432640_DP
    gbasis(32)%expo(5,3)=1.670815538_DP
    gbasis(32)%coco(5,1)=0.00495151_DP
    gbasis(32)%coco(5,2)=0.57776647_DP
    gbasis(32)%coco(5,3)=0.48464604_DP

    gbasis(32)%angmom(6)=2
    gbasis(32)%expo(6,1)=14.196656190_DP
    gbasis(32)%expo(6,2)=4.330432640_DP
    gbasis(32)%expo(6,3)=1.670815538_DP
    gbasis(32)%coco(6,1)=0.21976795_DP
    gbasis(32)%coco(6,2)=0.65554736_DP
    gbasis(32)%coco(6,3)=0.28657326_DP

    gbasis(32)%angmom(7)=0
    gbasis(32)%expo(7,1)=0.985832560_DP
    gbasis(32)%expo(7,2)=0.363834215_DP
    gbasis(32)%expo(7,3)=0.160673025_DP
    gbasis(32)%coco(7,1)=-0.30884412_DP
    gbasis(32)%coco(7,2)=0.01960641_DP
    gbasis(32)%coco(7,3)=1.13103444_DP

    gbasis(32)%angmom(8)=1
    gbasis(32)%expo(8,1)=0.985832560_DP
    gbasis(32)%expo(8,2)=0.363834215_DP
    gbasis(32)%expo(8,3)=0.160673025_DP
    gbasis(32)%coco(8,1)=-0.12154686_DP
    gbasis(32)%coco(8,2)=0.57152276_DP
    gbasis(32)%coco(8,3)=0.54989495_DP

    gbasis(32)%angmom(9)=2
    gbasis(32)%expo(9,1)=0.24600000000_DP
    gbasis(32)%coco(9,1)=1.00000000_DP
    ! ====end Ge ==========================================


    ! ====begin As ==========================================
    gbasis(33)%angmom(1)=0
    gbasis(33)%expo(1,1)=2337.065673000_DP
    gbasis(33)%expo(1,2)=425.699429800_DP
    gbasis(33)%expo(1,3)=115.210879000_DP
    gbasis(33)%coco(1,1)=0.15432897_DP
    gbasis(33)%coco(1,2)=0.53532814_DP
    gbasis(33)%coco(1,3)=0.44463454_DP

    gbasis(33)%angmom(2)=0
    gbasis(33)%expo(2,1)=193.197053500_DP
    gbasis(33)%expo(2,2)=44.894840400_DP
    gbasis(33)%expo(2,3)=14.601195480_DP
    gbasis(33)%coco(2,1)=-0.09996723_DP
    gbasis(33)%coco(2,2)=0.39951283_DP
    gbasis(33)%coco(2,3)=0.70011547_DP

    gbasis(33)%angmom(3)=1
    gbasis(33)%expo(3,1)=193.197053500_DP
    gbasis(33)%expo(3,2)=44.894840400_DP
    gbasis(33)%expo(3,3)=14.601195480_DP
    gbasis(33)%coco(3,1)=0.15591628_DP
    gbasis(33)%coco(3,2)=0.60768372_DP
    gbasis(33)%coco(3,3)=0.39195739_DP

    gbasis(33)%angmom(4)=0
    gbasis(33)%expo(4,1)=15.871635840_DP
    gbasis(33)%expo(4,2)=4.841354819_DP
    gbasis(33)%expo(4,3)=1.867945198_DP
    gbasis(33)%coco(4,1)=-0.22776350_DP
    gbasis(33)%coco(4,2)=0.21754360_DP
    gbasis(33)%coco(4,3)=0.91667696_DP

    gbasis(33)%angmom(5)=1
    gbasis(33)%expo(5,1)=15.871635840_DP
    gbasis(33)%expo(5,2)=4.841354819_DP
    gbasis(33)%expo(5,3)=1.867945198_DP
    gbasis(33)%coco(5,1)=0.00495151_DP
    gbasis(33)%coco(5,2)=0.57776647_DP
    gbasis(33)%coco(5,3)=0.48464604_DP

    gbasis(33)%angmom(6)=2
    gbasis(33)%expo(6,1)=15.871635840_DP
    gbasis(33)%expo(6,2)=4.841354819_DP
    gbasis(33)%expo(6,3)=1.867945198_DP
    gbasis(33)%coco(6,1)=0.21976795_DP
    gbasis(33)%coco(6,2)=0.65554736_DP
    gbasis(33)%coco(6,3)=0.28657326_DP

    gbasis(33)%angmom(7)=0
    gbasis(33)%expo(7,1)=1.107681464_DP
    gbasis(33)%expo(7,2)=0.408804124_DP
    gbasis(33)%expo(7,3)=0.180532211_DP
    gbasis(33)%coco(7,1)=-0.30884412_DP
    gbasis(33)%coco(7,2)=0.01960641_DP
    gbasis(33)%coco(7,3)=1.13103444_DP

    gbasis(33)%angmom(8)=1
    gbasis(33)%expo(8,1)=1.107681464_DP
    gbasis(33)%expo(8,2)=0.408804124_DP
    gbasis(33)%expo(8,3)=0.180532211_DP
    gbasis(33)%coco(8,1)=-0.12154686_DP
    gbasis(33)%coco(8,2)=0.57152276_DP
    gbasis(33)%coco(8,3)=0.54989495_DP

    gbasis(33)%angmom(9)=2
    gbasis(33)%expo(9,1)=0.29300000000_DP
    gbasis(33)%coco(9,1)=1.00000000_DP
    ! ====end As ==========================================


    ! ====begin Se ==========================================
    gbasis(34)%angmom(1)=0
    gbasis(34)%expo(1,1)=2480.626814000_DP
    gbasis(34)%expo(1,2)=451.849270800_DP
    gbasis(34)%expo(1,3)=122.288046400_DP
    gbasis(34)%coco(1,1)=0.15432897_DP
    gbasis(34)%coco(1,2)=0.53532814_DP
    gbasis(34)%coco(1,3)=0.44463454_DP

    gbasis(34)%angmom(2)=0
    gbasis(34)%expo(2,1)=206.157878000_DP
    gbasis(34)%expo(2,2)=47.906657270_DP
    gbasis(34)%expo(2,3)=15.580731800_DP
    gbasis(34)%coco(2,1)=-0.09996723_DP
    gbasis(34)%coco(2,2)=0.39951283_DP
    gbasis(34)%coco(2,3)=0.70011547_DP

    gbasis(34)%angmom(3)=1
    gbasis(34)%expo(3,1)=206.157878000_DP
    gbasis(34)%expo(3,2)=47.906657270_DP
    gbasis(34)%expo(3,3)=15.580731800_DP
    gbasis(34)%coco(3,1)=0.15591628_DP
    gbasis(34)%coco(3,2)=0.60768372_DP
    gbasis(34)%coco(3,3)=0.39195739_DP

    gbasis(34)%angmom(4)=0
    gbasis(34)%expo(4,1)=17.639994140_DP
    gbasis(34)%expo(4,2)=5.380760465_DP
    gbasis(34)%expo(4,3)=2.076064666_DP
    gbasis(34)%coco(4,1)=-0.22776350_DP
    gbasis(34)%coco(4,2)=0.21754360_DP
    gbasis(34)%coco(4,3)=0.91667696_DP

    gbasis(34)%angmom(5)=1
    gbasis(34)%expo(5,1)=17.639994140_DP
    gbasis(34)%expo(5,2)=5.380760465_DP
    gbasis(34)%expo(5,3)=2.076064666_DP
    gbasis(34)%coco(5,1)=0.00495151_DP
    gbasis(34)%coco(5,2)=0.57776647_DP
    gbasis(34)%coco(5,3)=0.48464604_DP

    gbasis(34)%angmom(6)=2
    gbasis(34)%expo(6,1)=17.639994140_DP
    gbasis(34)%expo(6,2)=5.380760465_DP
    gbasis(34)%expo(6,3)=2.076064666_DP
    gbasis(34)%coco(6,1)=0.21976795_DP
    gbasis(34)%coco(6,2)=0.65554736_DP
    gbasis(34)%coco(6,3)=0.28657326_DP

    gbasis(34)%angmom(7)=0
    gbasis(34)%expo(7,1)=1.214644297_DP
    gbasis(34)%expo(7,2)=0.448280136_DP
    gbasis(34)%expo(7,3)=0.197965235_DP
    gbasis(34)%coco(7,1)=-0.30884412_DP
    gbasis(34)%coco(7,2)=0.01960641_DP
    gbasis(34)%coco(7,3)=1.13103444_DP

    gbasis(34)%angmom(8)=1
    gbasis(34)%expo(8,1)=1.214644297_DP
    gbasis(34)%expo(8,2)=0.448280136_DP
    gbasis(34)%expo(8,3)=0.197965235_DP
    gbasis(34)%coco(8,1)=-0.12154686_DP
    gbasis(34)%coco(8,2)=0.57152276_DP
    gbasis(34)%coco(8,3)=0.54989495_DP

    gbasis(34)%angmom(9)=2
    gbasis(34)%expo(9,1)=0.33800000000_DP
    gbasis(34)%coco(9,1)=1.00000000_DP
    ! ====end Se ==========================================


    ! ====begin Br ==========================================
    gbasis(35)%angmom(1)=0
    gbasis(35)%expo(1,1)=2629.997471000_DP
    gbasis(35)%expo(1,2)=479.057322400_DP
    gbasis(35)%expo(1,3)=129.651607000_DP
    gbasis(35)%coco(1,1)=0.15432897_DP
    gbasis(35)%coco(1,2)=0.53532814_DP
    gbasis(35)%coco(1,3)=0.44463454_DP

    gbasis(35)%angmom(2)=0
    gbasis(35)%expo(2,1)=219.835025500_DP
    gbasis(35)%expo(2,2)=51.084932220_DP
    gbasis(35)%expo(2,3)=16.614405460_DP
    gbasis(35)%coco(2,1)=-0.09996723_DP
    gbasis(35)%coco(2,2)=0.39951283_DP
    gbasis(35)%coco(2,3)=0.70011547_DP

    gbasis(35)%angmom(3)=1
    gbasis(35)%expo(3,1)=219.835025500_DP
    gbasis(35)%expo(3,2)=51.084932220_DP
    gbasis(35)%expo(3,3)=16.614405460_DP
    gbasis(35)%coco(3,1)=0.15591628_DP
    gbasis(35)%coco(3,2)=0.60768372_DP
    gbasis(35)%coco(3,3)=0.39195739_DP

    gbasis(35)%angmom(4)=0
    gbasis(35)%expo(4,1)=19.501731090_DP
    gbasis(35)%expo(4,2)=5.948649577_DP
    gbasis(35)%expo(4,3)=2.295173940_DP
    gbasis(35)%coco(4,1)=-0.22776350_DP
    gbasis(35)%coco(4,2)=0.21754360_DP
    gbasis(35)%coco(4,3)=0.91667696_DP

    gbasis(35)%angmom(5)=1
    gbasis(35)%expo(5,1)=19.501731090_DP
    gbasis(35)%expo(5,2)=5.948649577_DP
    gbasis(35)%expo(5,3)=2.295173940_DP
    gbasis(35)%coco(5,1)=0.00495151_DP
    gbasis(35)%coco(5,2)=0.57776647_DP
    gbasis(35)%coco(5,3)=0.48464604_DP

    gbasis(35)%angmom(6)=2
    gbasis(35)%expo(6,1)=19.501731090_DP
    gbasis(35)%expo(6,2)=5.948649577_DP
    gbasis(35)%expo(6,3)=2.295173940_DP
    gbasis(35)%coco(6,1)=0.21976795_DP
    gbasis(35)%coco(6,2)=0.65554736_DP
    gbasis(35)%coco(6,3)=0.28657326_DP

    gbasis(35)%angmom(7)=0
    gbasis(35)%expo(7,1)=1.396037488_DP
    gbasis(35)%expo(7,2)=0.515225632_DP
    gbasis(35)%expo(7,3)=0.227529071_DP
    gbasis(35)%coco(7,1)=-0.30884412_DP
    gbasis(35)%coco(7,2)=0.01960641_DP
    gbasis(35)%coco(7,3)=1.13103444_DP

    gbasis(35)%angmom(8)=1
    gbasis(35)%expo(8,1)=1.396037488_DP
    gbasis(35)%expo(8,2)=0.515225632_DP
    gbasis(35)%expo(8,3)=0.227529071_DP
    gbasis(35)%coco(8,1)=-0.12154686_DP
    gbasis(35)%coco(8,2)=0.57152276_DP
    gbasis(35)%coco(8,3)=0.54989495_DP

    gbasis(35)%angmom(9)=2
    gbasis(35)%expo(9,1)=0.38900000000_DP
    gbasis(35)%coco(9,1)=1.00000000_DP
    ! ====end Br ==========================================


    ! ====begin Kr ==========================================
    gbasis(36)%angmom(1)=0
    gbasis(36)%expo(1,1)=2782.160055000_DP
    gbasis(36)%expo(1,2)=506.773927000_DP
    gbasis(36)%expo(1,3)=137.152801900_DP
    gbasis(36)%coco(1,1)=0.15432897_DP
    gbasis(36)%coco(1,2)=0.53532814_DP
    gbasis(36)%coco(1,3)=0.44463454_DP

    gbasis(36)%angmom(2)=0
    gbasis(36)%expo(2,1)=233.951411800_DP
    gbasis(36)%expo(2,2)=54.365276810_DP
    gbasis(36)%expo(2,3)=17.681275330_DP
    gbasis(36)%coco(2,1)=-0.09996723_DP
    gbasis(36)%coco(2,2)=0.39951283_DP
    gbasis(36)%coco(2,3)=0.70011547_DP

    gbasis(36)%angmom(3)=1
    gbasis(36)%expo(3,1)=233.951411800_DP
    gbasis(36)%expo(3,2)=54.365276810_DP
    gbasis(36)%expo(3,3)=17.681275330_DP
    gbasis(36)%coco(3,1)=0.15591628_DP
    gbasis(36)%coco(3,2)=0.60768372_DP
    gbasis(36)%coco(3,3)=0.39195739_DP

    gbasis(36)%angmom(4)=0
    gbasis(36)%expo(4,1)=21.456846710_DP
    gbasis(36)%expo(4,2)=6.545022156_DP
    gbasis(36)%expo(4,3)=2.525273021_DP
    gbasis(36)%coco(4,1)=-0.22776350_DP
    gbasis(36)%coco(4,2)=0.21754360_DP
    gbasis(36)%coco(4,3)=0.91667696_DP

    gbasis(36)%angmom(5)=1
    gbasis(36)%expo(5,1)=21.456846710_DP
    gbasis(36)%expo(5,2)=6.545022156_DP
    gbasis(36)%expo(5,3)=2.525273021_DP
    gbasis(36)%coco(5,1)=0.00495151_DP
    gbasis(36)%coco(5,2)=0.57776647_DP
    gbasis(36)%coco(5,3)=0.48464604_DP

    gbasis(36)%angmom(6)=2
    gbasis(36)%expo(6,1)=21.456846710_DP
    gbasis(36)%expo(6,2)=6.545022156_DP
    gbasis(36)%expo(6,3)=2.525273021_DP
    gbasis(36)%coco(6,1)=0.21976795_DP
    gbasis(36)%coco(6,2)=0.65554736_DP
    gbasis(36)%coco(6,3)=0.28657326_DP

    gbasis(36)%angmom(7)=0
    gbasis(36)%expo(7,1)=1.590049336_DP
    gbasis(36)%expo(7,2)=0.586828205_DP
    gbasis(36)%expo(7,3)=0.259149523_DP
    gbasis(36)%coco(7,1)=-0.30884412_DP
    gbasis(36)%coco(7,2)=0.01960641_DP
    gbasis(36)%coco(7,3)=1.13103444_DP

    gbasis(36)%angmom(8)=1
    gbasis(36)%expo(8,1)=1.590049336_DP
    gbasis(36)%expo(8,2)=0.586828205_DP
    gbasis(36)%expo(8,3)=0.259149523_DP
    gbasis(36)%coco(8,1)=-0.12154686_DP
    gbasis(36)%coco(8,2)=0.57152276_DP
    gbasis(36)%coco(8,3)=0.54989495_DP

    gbasis(36)%angmom(9)=2
    gbasis(36)%expo(9,1)=0.44300000000_DP
    gbasis(36)%coco(9,1)=1.00000000_DP
    ! ====end Kr ==========================================


    ! ====begin Rb ==========================================
    gbasis(37)%angmom(1)=0
    gbasis(37)%expo(1,1)=2938.601529000_DP
    gbasis(37)%expo(1,2)=535.269937000_DP
    gbasis(37)%expo(1,3)=144.864934000_DP
    gbasis(37)%coco(1,1)=0.15432897_DP
    gbasis(37)%coco(1,2)=0.53532814_DP
    gbasis(37)%coco(1,3)=0.44463454_DP

    gbasis(37)%angmom(2)=0
    gbasis(37)%expo(2,1)=248.507037000_DP
    gbasis(37)%expo(2,2)=57.747691000_DP
    gbasis(37)%expo(2,3)=18.781341000_DP
    gbasis(37)%coco(2,1)=-0.09996723_DP
    gbasis(37)%coco(2,2)=0.39951283_DP
    gbasis(37)%coco(2,3)=0.70011547_DP

    gbasis(37)%angmom(3)=1
    gbasis(37)%expo(3,1)=248.507037000_DP
    gbasis(37)%expo(3,2)=57.747691000_DP
    gbasis(37)%expo(3,3)=18.781341000_DP
    gbasis(37)%coco(3,1)=0.15591628_DP
    gbasis(37)%coco(3,2)=0.60768372_DP
    gbasis(37)%coco(3,3)=0.39195739_DP

    gbasis(37)%angmom(4)=0
    gbasis(37)%expo(4,1)=23.505340970_DP
    gbasis(37)%expo(4,2)=7.169878201_DP
    gbasis(37)%expo(4,3)=2.766361909_DP
    gbasis(37)%coco(4,1)=-0.22776350_DP
    gbasis(37)%coco(4,2)=0.21754360_DP
    gbasis(37)%coco(4,3)=0.91667696_DP

    gbasis(37)%angmom(5)=1
    gbasis(37)%expo(5,1)=23.505340970_DP
    gbasis(37)%expo(5,2)=7.169878201_DP
    gbasis(37)%expo(5,3)=2.766361909_DP
    gbasis(37)%coco(5,1)=0.00495151_DP
    gbasis(37)%coco(5,2)=0.57776647_DP
    gbasis(37)%coco(5,3)=0.48464604_DP

    gbasis(37)%angmom(6)=2
    gbasis(37)%expo(6,1)=23.505340970_DP
    gbasis(37)%expo(6,2)=7.169878201_DP
    gbasis(37)%expo(6,3)=2.766361909_DP
    gbasis(37)%coco(6,1)=0.21976795_DP
    gbasis(37)%coco(6,2)=0.65554736_DP
    gbasis(37)%coco(6,3)=0.28657326_DP

    gbasis(37)%angmom(7)=0
    gbasis(37)%expo(7,1)=2.247796820_DP
    gbasis(37)%expo(7,2)=0.829578393_DP
    gbasis(37)%expo(7,3)=0.366350565_DP
    gbasis(37)%coco(7,1)=-0.30884412_DP
    gbasis(37)%coco(7,2)=0.01960641_DP
    gbasis(37)%coco(7,3)=1.13103444_DP

    gbasis(37)%angmom(8)=1
    gbasis(37)%expo(8,1)=2.247796820_DP
    gbasis(37)%expo(8,2)=0.829578393_DP
    gbasis(37)%expo(8,3)=0.366350565_DP
    gbasis(37)%coco(8,1)=-0.12154686_DP
    gbasis(37)%coco(8,2)=0.57152276_DP
    gbasis(37)%coco(8,3)=0.54989495_DP

    gbasis(37)%angmom(9)=0
    gbasis(37)%expo(9,1)=0.486993992_DP
    gbasis(37)%expo(9,2)=0.262216156_DP
    gbasis(37)%expo(9,3)=0.115825488_DP
    gbasis(37)%coco(9,1)=-0.38426426_DP
    gbasis(37)%coco(9,2)=-0.19725674_DP
    gbasis(37)%coco(9,3)=1.37549551_DP

    gbasis(37)%angmom(10)=1
    gbasis(37)%expo(10,1)=0.486993992_DP
    gbasis(37)%expo(10,2)=0.262216156_DP
    gbasis(37)%expo(10,3)=0.115825488_DP
    gbasis(37)%coco(10,1)=-0.34816915_DP
    gbasis(37)%coco(10,2)=0.62903237_DP
    gbasis(37)%coco(10,3)=0.66628327_DP

    gbasis(37)%angmom(11)=2
    gbasis(37)%expo(11,1)=0.27432280711E-01_DP
    gbasis(37)%coco(11,1)=1.00000000_DP
    ! ====end Rb ==========================================


    ! ====begin Sr ==========================================
    gbasis(38)%angmom(1)=0
    gbasis(38)%expo(1,1)=3100.983951000_DP
    gbasis(38)%expo(1,2)=564.848097800_DP
    gbasis(38)%expo(1,3)=152.869938900_DP
    gbasis(38)%coco(1,1)=0.15432897_DP
    gbasis(38)%coco(1,2)=0.53532814_DP
    gbasis(38)%coco(1,3)=0.44463454_DP

    gbasis(38)%angmom(2)=0
    gbasis(38)%expo(2,1)=263.501900700_DP
    gbasis(38)%expo(2,2)=61.232174930_DP
    gbasis(38)%expo(2,3)=19.914603720_DP
    gbasis(38)%coco(2,1)=-0.09996723_DP
    gbasis(38)%coco(2,2)=0.39951283_DP
    gbasis(38)%coco(2,3)=0.70011547_DP

    gbasis(38)%angmom(3)=1
    gbasis(38)%expo(3,1)=263.501900700_DP
    gbasis(38)%expo(3,2)=61.232174930_DP
    gbasis(38)%expo(3,3)=19.914603720_DP
    gbasis(38)%coco(3,1)=0.15591628_DP
    gbasis(38)%coco(3,2)=0.60768372_DP
    gbasis(38)%coco(3,3)=0.39195739_DP

    gbasis(38)%angmom(4)=0
    gbasis(38)%expo(4,1)=25.578866920_DP
    gbasis(38)%expo(4,2)=7.802369707_DP
    gbasis(38)%expo(4,3)=3.010396794_DP
    gbasis(38)%coco(4,1)=-0.22776350_DP
    gbasis(38)%coco(4,2)=0.21754360_DP
    gbasis(38)%coco(4,3)=0.91667696_DP

    gbasis(38)%angmom(5)=1
    gbasis(38)%expo(5,1)=25.578866920_DP
    gbasis(38)%expo(5,2)=7.802369707_DP
    gbasis(38)%expo(5,3)=3.010396794_DP
    gbasis(38)%coco(5,1)=0.00495151_DP
    gbasis(38)%coco(5,2)=0.57776647_DP
    gbasis(38)%coco(5,3)=0.48464604_DP

    gbasis(38)%angmom(6)=2
    gbasis(38)%expo(6,1)=25.578866920_DP
    gbasis(38)%expo(6,2)=7.802369707_DP
    gbasis(38)%expo(6,3)=3.010396794_DP
    gbasis(38)%coco(6,1)=0.21976795_DP
    gbasis(38)%coco(6,2)=0.65554736_DP
    gbasis(38)%coco(6,3)=0.28657326_DP

    gbasis(38)%angmom(7)=0
    gbasis(38)%expo(7,1)=2.461032403_DP
    gbasis(38)%expo(7,2)=0.908275734_DP
    gbasis(38)%expo(7,3)=0.401104140_DP
    gbasis(38)%coco(7,1)=-0.30884412_DP
    gbasis(38)%coco(7,2)=0.01960641_DP
    gbasis(38)%coco(7,3)=1.13103444_DP

    gbasis(38)%angmom(8)=1
    gbasis(38)%expo(8,1)=2.461032403_DP
    gbasis(38)%expo(8,2)=0.908275734_DP
    gbasis(38)%expo(8,3)=0.401104140_DP
    gbasis(38)%coco(8,1)=-0.12154686_DP
    gbasis(38)%coco(8,2)=0.57152276_DP
    gbasis(38)%coco(8,3)=0.54989495_DP

    gbasis(38)%angmom(9)=0
    gbasis(38)%expo(9,1)=0.437080480_DP
    gbasis(38)%expo(9,2)=0.235340816_DP
    gbasis(38)%expo(9,3)=0.103954177_DP
    gbasis(38)%coco(9,1)=-0.38426426_DP
    gbasis(38)%coco(9,2)=-0.19725674_DP
    gbasis(38)%coco(9,3)=1.37549551_DP

    gbasis(38)%angmom(10)=1
    gbasis(38)%expo(10,1)=0.437080480_DP
    gbasis(38)%expo(10,2)=0.235340816_DP
    gbasis(38)%expo(10,3)=0.103954177_DP
    gbasis(38)%coco(10,1)=-0.34816915_DP
    gbasis(38)%coco(10,2)=0.62903237_DP
    gbasis(38)%coco(10,3)=0.66628327_DP

    gbasis(38)%angmom(11)=2
    gbasis(38)%expo(11,1)=0.36655000000E-01_DP
    gbasis(38)%coco(11,1)=0.44818300000_DP
    ! ====end Sr ==========================================


    ! ====begin Y ==========================================
    gbasis(39)%angmom(1)=0
    gbasis(39)%expo(1,1)=3266.026869000_DP
    gbasis(39)%expo(1,2)=594.910871000_DP
    gbasis(39)%expo(1,3)=161.006099000_DP
    gbasis(39)%coco(1,1)=0.15432897_DP
    gbasis(39)%coco(1,2)=0.53532814_DP
    gbasis(39)%coco(1,3)=0.44463454_DP

    gbasis(39)%angmom(2)=0
    gbasis(39)%expo(2,1)=277.937724000_DP
    gbasis(39)%expo(2,2)=64.586750000_DP
    gbasis(39)%expo(2,3)=21.005616000_DP
    gbasis(39)%coco(2,1)=-0.09996723_DP
    gbasis(39)%coco(2,2)=0.39951283_DP
    gbasis(39)%coco(2,3)=0.70011547_DP

    gbasis(39)%angmom(3)=1
    gbasis(39)%expo(3,1)=277.937724000_DP
    gbasis(39)%expo(3,2)=64.586750000_DP
    gbasis(39)%expo(3,3)=21.005616000_DP
    gbasis(39)%coco(3,1)=0.15591628_DP
    gbasis(39)%coco(3,2)=0.60768372_DP
    gbasis(39)%coco(3,3)=0.39195739_DP

    gbasis(39)%angmom(4)=0
    gbasis(39)%expo(4,1)=30.671326000_DP
    gbasis(39)%expo(4,2)=8.557222000_DP
    gbasis(39)%expo(4,3)=3.349239000_DP
    gbasis(39)%coco(4,1)=-0.22776350_DP
    gbasis(39)%coco(4,2)=0.21754360_DP
    gbasis(39)%coco(4,3)=0.91667696_DP

    gbasis(39)%angmom(5)=1
    gbasis(39)%expo(5,1)=30.671326000_DP
    gbasis(39)%expo(5,2)=8.557222000_DP
    gbasis(39)%expo(5,3)=3.349239000_DP
    gbasis(39)%coco(5,1)=0.00495151_DP
    gbasis(39)%coco(5,2)=0.57776647_DP
    gbasis(39)%coco(5,3)=0.48464604_DP

    gbasis(39)%angmom(6)=2
    gbasis(39)%expo(6,1)=5.660043000_DP
    gbasis(39)%expo(6,2)=1.774715000_DP
    gbasis(39)%expo(6,3)=0.691295000_DP
    gbasis(39)%coco(6,1)=0.21976795_DP
    gbasis(39)%coco(6,2)=0.65554736_DP
    gbasis(39)%coco(6,3)=0.28657326_DP

    gbasis(39)%angmom(7)=0
    gbasis(39)%expo(7,1)=2.667688000_DP
    gbasis(39)%expo(7,2)=0.984544000_DP
    gbasis(39)%expo(7,3)=0.434785000_DP
    gbasis(39)%coco(7,1)=-0.33061006_DP
    gbasis(39)%coco(7,2)=0.05761095_DP
    gbasis(39)%coco(7,3)=1.15578745_DP

    gbasis(39)%angmom(8)=1
    gbasis(39)%expo(8,1)=2.667688000_DP
    gbasis(39)%expo(8,2)=0.984544000_DP
    gbasis(39)%expo(8,3)=0.434785000_DP
    gbasis(39)%coco(8,1)=-0.12839276_DP
    gbasis(39)%coco(8,2)=0.58520476_DP
    gbasis(39)%coco(8,3)=0.54394420_DP

    gbasis(39)%angmom(9)=2
    gbasis(39)%expo(9,1)=2.128212000_DP
    gbasis(39)%expo(9,2)=0.962594000_DP
    gbasis(39)%expo(9,3)=0.472861000_DP
    gbasis(39)%coco(9,1)=0.12506621_DP
    gbasis(39)%coco(9,2)=0.66867856_DP
    gbasis(39)%coco(9,3)=0.30524682_DP

    gbasis(39)%angmom(10)=0
    gbasis(39)%expo(10,1)=0.207424000_DP
    gbasis(39)%expo(10,2)=0.111685000_DP
    gbasis(39)%expo(10,3)=0.049333000_DP
    gbasis(39)%coco(10,1)=-0.38426426_DP
    gbasis(39)%coco(10,2)=-0.19725674_DP
    gbasis(39)%coco(10,3)=1.37549551_DP

    gbasis(39)%angmom(11)=1
    gbasis(39)%expo(11,1)=0.207424000_DP
    gbasis(39)%expo(11,2)=0.111685000_DP
    gbasis(39)%expo(11,3)=0.049333000_DP
    gbasis(39)%coco(11,1)=-0.34816915_DP
    gbasis(39)%coco(11,2)=0.62903237_DP
    gbasis(39)%coco(11,3)=0.66628327_DP

    gbasis(39)%angmom(12)=3
    gbasis(39)%expo(12,1)=0.26515000000_DP
    gbasis(39)%coco(12,1)=1.0000000000_DP
    ! ====end Y ==========================================


    ! ====begin Zr ==========================================
    gbasis(40)%angmom(1)=0
    gbasis(40)%expo(1,1)=3435.348677000_DP
    gbasis(40)%expo(1,2)=625.753049800_DP
    gbasis(40)%expo(1,3)=169.353195800_DP
    gbasis(40)%coco(1,1)=0.15432897_DP
    gbasis(40)%coco(1,2)=0.53532814_DP
    gbasis(40)%coco(1,3)=0.44463454_DP

    gbasis(40)%angmom(2)=0
    gbasis(40)%expo(2,1)=293.783029200_DP
    gbasis(40)%expo(2,2)=68.268857970_DP
    gbasis(40)%expo(2,3)=22.203151440_DP
    gbasis(40)%coco(2,1)=-0.09996723_DP
    gbasis(40)%coco(2,2)=0.39951283_DP
    gbasis(40)%coco(2,3)=0.70011547_DP

    gbasis(40)%angmom(3)=1
    gbasis(40)%expo(3,1)=293.783029200_DP
    gbasis(40)%expo(3,2)=68.268857970_DP
    gbasis(40)%expo(3,3)=22.203151440_DP
    gbasis(40)%coco(3,1)=0.15591628_DP
    gbasis(40)%coco(3,2)=0.60768372_DP
    gbasis(40)%coco(3,3)=0.39195739_DP

    gbasis(40)%angmom(4)=0
    gbasis(40)%expo(4,1)=30.732931030_DP
    gbasis(40)%expo(4,2)=9.374523538_DP
    gbasis(40)%expo(4,3)=3.616982618_DP
    gbasis(40)%coco(4,1)=-0.22776350_DP
    gbasis(40)%coco(4,2)=0.21754360_DP
    gbasis(40)%coco(4,3)=0.91667696_DP

    gbasis(40)%angmom(5)=1
    gbasis(40)%expo(5,1)=30.732931030_DP
    gbasis(40)%expo(5,2)=9.374523538_DP
    gbasis(40)%expo(5,3)=3.616982618_DP
    gbasis(40)%coco(5,1)=0.00495151_DP
    gbasis(40)%coco(5,2)=0.57776647_DP
    gbasis(40)%coco(5,3)=0.48464604_DP

    gbasis(40)%angmom(6)=2
    gbasis(40)%expo(6,1)=30.732931030_DP
    gbasis(40)%expo(6,2)=9.374523538_DP
    gbasis(40)%expo(6,3)=3.616982618_DP
    gbasis(40)%coco(6,1)=0.21976795_DP
    gbasis(40)%coco(6,2)=0.65554736_DP
    gbasis(40)%coco(6,3)=0.28657326_DP

    gbasis(40)%angmom(7)=0
    gbasis(40)%expo(7,1)=2.827607815_DP
    gbasis(40)%expo(7,2)=1.101055827_DP
    gbasis(40)%expo(7,3)=0.484687486_DP
    gbasis(40)%coco(7,1)=-0.33061006_DP
    gbasis(40)%coco(7,2)=0.05761095_DP
    gbasis(40)%coco(7,3)=1.15578745_DP

    gbasis(40)%angmom(8)=1
    gbasis(40)%expo(8,1)=2.827607815_DP
    gbasis(40)%expo(8,2)=1.101055827_DP
    gbasis(40)%expo(8,3)=0.484687486_DP
    gbasis(40)%coco(8,1)=-0.12839276_DP
    gbasis(40)%coco(8,2)=0.58520476_DP
    gbasis(40)%coco(8,3)=0.54394420_DP

    gbasis(40)%angmom(9)=2
    gbasis(40)%expo(9,1)=0.486993992_DP
    gbasis(40)%expo(9,2)=0.262216156_DP
    gbasis(40)%expo(9,3)=0.115825488_DP
    gbasis(40)%coco(9,1)=0.12506621_DP
    gbasis(40)%coco(9,2)=0.66867856_DP
    gbasis(40)%coco(9,3)=0.30524682_DP

    gbasis(40)%angmom(10)=0
    gbasis(40)%expo(10,1)=0.887830189_DP
    gbasis(40)%expo(10,2)=0.345716474_DP
    gbasis(40)%expo(10,3)=0.152185243_DP
    gbasis(40)%coco(10,1)=-0.38426426_DP
    gbasis(40)%coco(10,2)=-0.19725674_DP
    gbasis(40)%coco(10,3)=1.37549551_DP

    gbasis(40)%angmom(11)=1
    gbasis(40)%expo(11,1)=0.887830189_DP
    gbasis(40)%expo(11,2)=0.345716474_DP
    gbasis(40)%expo(11,3)=0.152185243_DP
    gbasis(40)%coco(11,1)=-0.34816915_DP
    gbasis(40)%coco(11,2)=0.62903237_DP
    gbasis(40)%coco(11,3)=0.66628327_DP

    gbasis(40)%angmom(12)=3
    gbasis(40)%expo(12,1)=0.39261000000_DP
    gbasis(40)%coco(12,1)=1.0000000000_DP
    ! ====end Zr ==========================================


  end subroutine ngwf_data_cocos_and_expos21_40

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_data_cocos_and_expos41_60(gbasis)

    !===============================================================!
    ! Initialisation of Gaussian basis set exponents and contraction!
    ! coefficients. All the data of the STO-3G basis set are        !
    ! included for each element of the periodic table as well as    !
    ! the one "polarisation" function of the 6-31G* set.            !
    !---------------------------------------------------------------!
    ! This subroutine was created by Chris-Kriton Skylaris on       !
    ! 31/03/2006 but the actual STO-3G and 6-31G* parameters        !
    ! which appear below were typed by Shyong Chen, in March 2006.  !
    !===============================================================!

    implicit none

    type(GTO_SET), intent(inout) :: gbasis(1: pub_num_gtatoms)


    ! cks: all data in this subroutine was typed by Shyong Chen

    ! ====begin Nb ==========================================
    gbasis(41)%angmom(1)=0
    gbasis(41)%expo(1,1)=3610.742864000_DP
    gbasis(41)%expo(1,2)=657.701320100_DP
    gbasis(41)%expo(1,3)=177.999644500_DP
    gbasis(41)%coco(1,1)=0.15432897_DP
    gbasis(41)%coco(1,2)=0.53532814_DP
    gbasis(41)%coco(1,3)=0.44463454_DP

    gbasis(41)%angmom(2)=0
    gbasis(41)%expo(2,1)=310.067572800_DP
    gbasis(41)%expo(2,2)=72.053035690_DP
    gbasis(41)%expo(2,3)=23.433883480_DP
    gbasis(41)%coco(2,1)=-0.09996723_DP
    gbasis(41)%coco(2,2)=0.39951283_DP
    gbasis(41)%coco(2,3)=0.70011547_DP

    gbasis(41)%angmom(3)=1
    gbasis(41)%expo(3,1)=310.067572800_DP
    gbasis(41)%expo(3,2)=72.053035690_DP
    gbasis(41)%expo(3,3)=23.433883480_DP
    gbasis(41)%coco(3,1)=0.15591628_DP
    gbasis(41)%coco(3,2)=0.60768372_DP
    gbasis(41)%coco(3,3)=0.39195739_DP

    gbasis(41)%angmom(4)=0
    gbasis(41)%expo(4,1)=33.019978580_DP
    gbasis(41)%expo(4,2)=10.072145940_DP
    gbasis(41)%expo(4,3)=3.886147028_DP
    gbasis(41)%coco(4,1)=-0.22776350_DP
    gbasis(41)%coco(4,2)=0.21754360_DP
    gbasis(41)%coco(4,3)=0.91667696_DP

    gbasis(41)%angmom(5)=1
    gbasis(41)%expo(5,1)=33.019978580_DP
    gbasis(41)%expo(5,2)=10.072145940_DP
    gbasis(41)%expo(5,3)=3.886147028_DP
    gbasis(41)%coco(5,1)=0.00495151_DP
    gbasis(41)%coco(5,2)=0.57776647_DP
    gbasis(41)%coco(5,3)=0.48464604_DP

    gbasis(41)%angmom(6)=2
    gbasis(41)%expo(6,1)=33.019978580_DP
    gbasis(41)%expo(6,2)=10.072145940_DP
    gbasis(41)%expo(6,3)=3.886147028_DP
    gbasis(41)%coco(6,1)=0.21976795_DP
    gbasis(41)%coco(6,2)=0.65554736_DP
    gbasis(41)%coco(6,3)=0.28657326_DP

    gbasis(41)%angmom(7)=0
    gbasis(41)%expo(7,1)=3.144798430_DP
    gbasis(41)%expo(7,2)=1.224568208_DP
    gbasis(41)%expo(7,3)=0.539057940_DP
    gbasis(41)%coco(7,1)=-0.33061006_DP
    gbasis(41)%coco(7,2)=0.05761095_DP
    gbasis(41)%coco(7,3)=1.15578745_DP

    gbasis(41)%angmom(8)=1
    gbasis(41)%expo(8,1)=3.144798430_DP
    gbasis(41)%expo(8,2)=1.224568208_DP
    gbasis(41)%expo(8,3)=0.539057940_DP
    gbasis(41)%coco(8,1)=-0.12839276_DP
    gbasis(41)%coco(8,2)=0.58520476_DP
    gbasis(41)%coco(8,3)=0.54394420_DP

    gbasis(41)%angmom(9)=2
    gbasis(41)%expo(9,1)=1.344878866_DP
    gbasis(41)%expo(9,2)=0.523688859_DP
    gbasis(41)%expo(9,3)=0.230529125_DP
    gbasis(41)%coco(9,1)=0.12506621_DP
    gbasis(41)%coco(9,2)=0.66867856_DP
    gbasis(41)%coco(9,3)=0.30524682_DP

    gbasis(41)%angmom(10)=0
    gbasis(41)%expo(10,1)=0.486993992_DP
    gbasis(41)%expo(10,2)=0.262216156_DP
    gbasis(41)%expo(10,3)=0.115825488_DP
    gbasis(41)%coco(10,1)=-0.38426426_DP
    gbasis(41)%coco(10,2)=-0.19725674_DP
    gbasis(41)%coco(10,3)=1.37549551_DP

    gbasis(41)%angmom(11)=1
    gbasis(41)%expo(11,1)=0.486993992_DP
    gbasis(41)%expo(11,2)=0.262216156_DP
    gbasis(41)%expo(11,3)=0.115825488_DP
    gbasis(41)%coco(11,1)=-0.34816915_DP
    gbasis(41)%coco(11,2)=0.62903237_DP
    gbasis(41)%coco(11,3)=0.66628327_DP

    gbasis(41)%angmom(12)=3
    gbasis(41)%expo(12,1)=0.52270000000_DP
    gbasis(41)%coco(12,1)=1.0000000000_DP
    ! ====end Nb ==========================================


    ! ====begin Mo ==========================================
    gbasis(42)%angmom(1)=0
    gbasis(42)%expo(1,1)=3788.666115000_DP
    gbasis(42)%expo(1,2)=690.110262300_DP
    gbasis(42)%expo(1,3)=186.770769100_DP
    gbasis(42)%coco(1,1)=0.15432897_DP
    gbasis(42)%coco(1,2)=0.53532814_DP
    gbasis(42)%coco(1,3)=0.44463454_DP

    gbasis(42)%angmom(2)=0
    gbasis(42)%expo(2,1)=326.430956700_DP
    gbasis(42)%expo(2,2)=75.855534200_DP
    gbasis(42)%expo(2,3)=24.670574010_DP
    gbasis(42)%coco(2,1)=-0.09996723_DP
    gbasis(42)%coco(2,2)=0.39951283_DP
    gbasis(42)%coco(2,3)=0.70011547_DP

    gbasis(42)%angmom(3)=1
    gbasis(42)%expo(3,1)=326.430956700_DP
    gbasis(42)%expo(3,2)=75.855534200_DP
    gbasis(42)%expo(3,3)=24.670574010_DP
    gbasis(42)%coco(3,1)=0.15591628_DP
    gbasis(42)%coco(3,2)=0.60768372_DP
    gbasis(42)%coco(3,3)=0.39195739_DP

    gbasis(42)%angmom(4)=0
    gbasis(42)%expo(4,1)=35.469481290_DP
    gbasis(42)%expo(4,2)=10.819322340_DP
    gbasis(42)%expo(4,3)=4.174430912_DP
    gbasis(42)%coco(4,1)=-0.22776350_DP
    gbasis(42)%coco(4,2)=0.21754360_DP
    gbasis(42)%coco(4,3)=0.91667696_DP

    gbasis(42)%angmom(5)=1
    gbasis(42)%expo(5,1)=35.469481290_DP
    gbasis(42)%expo(5,2)=10.819322340_DP
    gbasis(42)%expo(5,3)=4.174430912_DP
    gbasis(42)%coco(5,1)=0.00495151_DP
    gbasis(42)%coco(5,2)=0.57776647_DP
    gbasis(42)%coco(5,3)=0.48464604_DP

    gbasis(42)%angmom(6)=2
    gbasis(42)%expo(6,1)=35.469481290_DP
    gbasis(42)%expo(6,2)=10.819322340_DP
    gbasis(42)%expo(6,3)=4.174430912_DP
    gbasis(42)%coco(6,1)=0.21976795_DP
    gbasis(42)%coco(6,2)=0.65554736_DP
    gbasis(42)%coco(6,3)=0.28657326_DP

    gbasis(42)%angmom(7)=0
    gbasis(42)%expo(7,1)=3.496895188_DP
    gbasis(42)%expo(7,2)=1.361672861_DP
    gbasis(42)%expo(7,3)=0.599411746_DP
    gbasis(42)%coco(7,1)=-0.33061006_DP
    gbasis(42)%coco(7,2)=0.05761095_DP
    gbasis(42)%coco(7,3)=1.15578745_DP

    gbasis(42)%angmom(8)=1
    gbasis(42)%expo(8,1)=3.496895188_DP
    gbasis(42)%expo(8,2)=1.361672861_DP
    gbasis(42)%expo(8,3)=0.599411746_DP
    gbasis(42)%coco(8,1)=-0.12839276_DP
    gbasis(42)%coco(8,2)=0.58520476_DP
    gbasis(42)%coco(8,3)=0.54394420_DP

    gbasis(42)%angmom(9)=2
    gbasis(42)%expo(9,1)=1.702112315_DP
    gbasis(42)%expo(9,2)=0.662793713_DP
    gbasis(42)%expo(9,3)=0.291763424_DP
    gbasis(42)%coco(9,1)=0.12506621_DP
    gbasis(42)%coco(9,2)=0.66867856_DP
    gbasis(42)%coco(9,3)=0.30524682_DP

    gbasis(42)%angmom(10)=0
    gbasis(42)%expo(10,1)=0.512962508_DP
    gbasis(42)%expo(10,2)=0.276198597_DP
    gbasis(42)%expo(10,3)=0.122001777_DP
    gbasis(42)%coco(10,1)=-0.38426426_DP
    gbasis(42)%coco(10,2)=-0.19725674_DP
    gbasis(42)%coco(10,3)=1.37549551_DP

    gbasis(42)%angmom(11)=1
    gbasis(42)%expo(11,1)=0.512962508_DP
    gbasis(42)%expo(11,2)=0.276198597_DP
    gbasis(42)%expo(11,3)=0.122001777_DP
    gbasis(42)%coco(11,1)=-0.34816915_DP
    gbasis(42)%coco(11,2)=0.62903237_DP
    gbasis(42)%coco(11,3)=0.66628327_DP

    gbasis(42)%angmom(12)=3
    gbasis(42)%expo(12,1)=0.65545000000_DP
    gbasis(42)%coco(12,1)=1.0000000000_DP
    ! ====end Mo ==========================================


    ! ====begin Tc ==========================================
    gbasis(43)%angmom(1)=0
    gbasis(43)%expo(1,1)=3970.868257000_DP
    gbasis(43)%expo(1,2)=723.298609800_DP
    gbasis(43)%expo(1,3)=195.752831100_DP
    gbasis(43)%coco(1,1)=0.15432897_DP
    gbasis(43)%coco(1,2)=0.53532814_DP
    gbasis(43)%coco(1,3)=0.44463454_DP

    gbasis(43)%angmom(2)=0
    gbasis(43)%expo(2,1)=343.584632300_DP
    gbasis(43)%expo(2,2)=79.841679520_DP
    gbasis(43)%expo(2,3)=25.966992190_DP
    gbasis(43)%coco(2,1)=-0.09996723_DP
    gbasis(43)%coco(2,2)=0.39951283_DP
    gbasis(43)%coco(2,3)=0.70011547_DP

    gbasis(43)%angmom(3)=1
    gbasis(43)%expo(3,1)=343.584632300_DP
    gbasis(43)%expo(3,2)=79.841679520_DP
    gbasis(43)%expo(3,3)=25.966992190_DP
    gbasis(43)%coco(3,1)=0.15591628_DP
    gbasis(43)%coco(3,2)=0.60768372_DP
    gbasis(43)%coco(3,3)=0.39195739_DP

    gbasis(43)%angmom(4)=0
    gbasis(43)%expo(4,1)=38.089919830_DP
    gbasis(43)%expo(4,2)=11.618639620_DP
    gbasis(43)%expo(4,3)=4.482832367_DP
    gbasis(43)%coco(4,1)=-0.22776350_DP
    gbasis(43)%coco(4,2)=0.21754360_DP
    gbasis(43)%coco(4,3)=0.91667696_DP

    gbasis(43)%angmom(5)=1
    gbasis(43)%expo(5,1)=38.089919830_DP
    gbasis(43)%expo(5,2)=11.618639620_DP
    gbasis(43)%expo(5,3)=4.482832367_DP
    gbasis(43)%coco(5,1)=0.00495151_DP
    gbasis(43)%coco(5,2)=0.57776647_DP
    gbasis(43)%coco(5,3)=0.48464604_DP

    gbasis(43)%angmom(6)=2
    gbasis(43)%expo(6,1)=38.089919830_DP
    gbasis(43)%expo(6,2)=11.618639620_DP
    gbasis(43)%expo(6,3)=4.482832367_DP
    gbasis(43)%coco(6,1)=0.21976795_DP
    gbasis(43)%coco(6,2)=0.65554736_DP
    gbasis(43)%coco(6,3)=0.28657326_DP

    gbasis(43)%angmom(7)=0
    gbasis(43)%expo(7,1)=3.829752708_DP
    gbasis(43)%expo(7,2)=1.491285854_DP
    gbasis(43)%expo(7,3)=0.656467704_DP
    gbasis(43)%coco(7,1)=-0.33061006_DP
    gbasis(43)%coco(7,2)=0.05761095_DP
    gbasis(43)%coco(7,3)=1.15578745_DP

    gbasis(43)%angmom(8)=1
    gbasis(43)%expo(8,1)=3.829752708_DP
    gbasis(43)%expo(8,2)=1.491285854_DP
    gbasis(43)%expo(8,3)=0.656467704_DP
    gbasis(43)%coco(8,1)=-0.12839276_DP
    gbasis(43)%coco(8,2)=0.58520476_DP
    gbasis(43)%coco(8,3)=0.54394420_DP

    gbasis(43)%angmom(9)=2
    gbasis(43)%expo(9,1)=2.101373228_DP
    gbasis(43)%expo(9,2)=0.818263843_DP
    gbasis(43)%expo(9,3)=0.360201758_DP
    gbasis(43)%coco(9,1)=0.12506621_DP
    gbasis(43)%coco(9,2)=0.66867856_DP
    gbasis(43)%coco(9,3)=0.30524682_DP

    gbasis(43)%angmom(10)=0
    gbasis(43)%expo(10,1)=0.461699983_DP
    gbasis(43)%expo(10,2)=0.248596896_DP
    gbasis(43)%expo(10,3)=0.109809621_DP
    gbasis(43)%coco(10,1)=-0.38426426_DP
    gbasis(43)%coco(10,2)=-0.19725674_DP
    gbasis(43)%coco(10,3)=1.37549551_DP

    gbasis(43)%angmom(11)=1
    gbasis(43)%expo(11,1)=0.461699983_DP
    gbasis(43)%expo(11,2)=0.248596896_DP
    gbasis(43)%expo(11,3)=0.109809621_DP
    gbasis(43)%coco(11,1)=-0.34816915_DP
    gbasis(43)%coco(11,2)=0.62903237_DP
    gbasis(43)%coco(11,3)=0.66628327_DP

    gbasis(43)%angmom(12)=3
    gbasis(43)%expo(12,1)=0.79085000000_DP
    gbasis(43)%coco(12,1)=1.0000000000_DP
    ! ====end Tc ==========================================


    ! ====begin Ru ==========================================
    gbasis(44)%angmom(1)=0
    gbasis(44)%expo(1,1)=4159.274210000_DP
    gbasis(44)%expo(1,2)=757.616989400_DP
    gbasis(44)%expo(1,3)=205.040723900_DP
    gbasis(44)%coco(1,1)=0.15432897_DP
    gbasis(44)%coco(1,2)=0.53532814_DP
    gbasis(44)%coco(1,3)=0.44463454_DP

    gbasis(44)%angmom(2)=0
    gbasis(44)%expo(2,1)=360.798656100_DP
    gbasis(44)%expo(2,2)=83.841848430_DP
    gbasis(44)%expo(2,3)=27.267971270_DP
    gbasis(44)%coco(2,1)=-0.09996723_DP
    gbasis(44)%coco(2,2)=0.39951283_DP
    gbasis(44)%coco(2,3)=0.70011547_DP

    gbasis(44)%angmom(3)=1
    gbasis(44)%expo(3,1)=360.798656100_DP
    gbasis(44)%expo(3,2)=83.841848430_DP
    gbasis(44)%expo(3,3)=27.267971270_DP
    gbasis(44)%coco(3,1)=0.15591628_DP
    gbasis(44)%coco(3,2)=0.60768372_DP
    gbasis(44)%coco(3,3)=0.39195739_DP

    gbasis(44)%angmom(4)=0
    gbasis(44)%expo(4,1)=40.717516780_DP
    gbasis(44)%expo(4,2)=12.420140440_DP
    gbasis(44)%expo(4,3)=4.792076302_DP
    gbasis(44)%coco(4,1)=-0.22776350_DP
    gbasis(44)%coco(4,2)=0.21754360_DP
    gbasis(44)%coco(4,3)=0.91667696_DP

    gbasis(44)%angmom(5)=1
    gbasis(44)%expo(5,1)=40.717516780_DP
    gbasis(44)%expo(5,2)=12.420140440_DP
    gbasis(44)%expo(5,3)=4.792076302_DP
    gbasis(44)%coco(5,1)=0.00495151_DP
    gbasis(44)%coco(5,2)=0.57776647_DP
    gbasis(44)%coco(5,3)=0.48464604_DP

    gbasis(44)%angmom(6)=2
    gbasis(44)%expo(6,1)=40.717516780_DP
    gbasis(44)%expo(6,2)=12.420140440_DP
    gbasis(44)%expo(6,3)=4.792076302_DP
    gbasis(44)%coco(6,1)=0.21976795_DP
    gbasis(44)%coco(6,2)=0.65554736_DP
    gbasis(44)%coco(6,3)=0.28657326_DP

    gbasis(44)%angmom(7)=0
    gbasis(44)%expo(7,1)=4.197516371_DP
    gbasis(44)%expo(7,2)=1.634491118_DP
    gbasis(44)%expo(7,3)=0.719507014_DP
    gbasis(44)%coco(7,1)=-0.33061006_DP
    gbasis(44)%coco(7,2)=0.05761095_DP
    gbasis(44)%coco(7,3)=1.15578745_DP

    gbasis(44)%angmom(8)=1
    gbasis(44)%expo(8,1)=4.197516371_DP
    gbasis(44)%expo(8,2)=1.634491118_DP
    gbasis(44)%expo(8,3)=0.719507014_DP
    gbasis(44)%coco(8,1)=-0.12839276_DP
    gbasis(44)%coco(8,2)=0.58520476_DP
    gbasis(44)%coco(8,3)=0.54394420_DP

    gbasis(44)%angmom(9)=2
    gbasis(44)%expo(9,1)=2.390895761_DP
    gbasis(44)%expo(9,2)=0.931002417_DP
    gbasis(44)%expo(9,3)=0.409829556_DP
    gbasis(44)%coco(9,1)=0.12506621_DP
    gbasis(44)%coco(9,2)=0.66867856_DP
    gbasis(44)%coco(9,3)=0.30524682_DP

    gbasis(44)%angmom(10)=0
    gbasis(44)%expo(10,1)=0.413135485_DP
    gbasis(44)%expo(10,2)=0.222447917_DP
    gbasis(44)%expo(10,3)=0.098259157_DP
    gbasis(44)%coco(10,1)=-0.38426426_DP
    gbasis(44)%coco(10,2)=-0.19725674_DP
    gbasis(44)%coco(10,3)=1.37549551_DP

    gbasis(44)%angmom(11)=1
    gbasis(44)%expo(11,1)=0.413135485_DP
    gbasis(44)%expo(11,2)=0.222447917_DP
    gbasis(44)%expo(11,3)=0.098259157_DP
    gbasis(44)%coco(11,1)=-0.34816915_DP
    gbasis(44)%coco(11,2)=0.62903237_DP
    gbasis(44)%coco(11,3)=0.66628327_DP

    gbasis(44)%angmom(12)=3
    gbasis(44)%expo(12,1)=0.94314000000_DP
    gbasis(44)%coco(12,1)=1.0000000000_DP
    ! ====end Ru ==========================================


    ! ====begin Rh ==========================================
    gbasis(45)%angmom(1)=0
    gbasis(45)%expo(1,1)=4350.077794000_DP
    gbasis(45)%expo(1,2)=792.372100500_DP
    gbasis(45)%expo(1,3)=214.446813300_DP
    gbasis(45)%coco(1,1)=0.15432897_DP
    gbasis(45)%coco(1,2)=0.53532814_DP
    gbasis(45)%coco(1,3)=0.44463454_DP

    gbasis(45)%angmom(2)=0
    gbasis(45)%expo(2,1)=378.433426400_DP
    gbasis(45)%expo(2,2)=87.939789810_DP
    gbasis(45)%expo(2,3)=28.600748990_DP
    gbasis(45)%coco(2,1)=-0.09996723_DP
    gbasis(45)%coco(2,2)=0.39951283_DP
    gbasis(45)%coco(2,3)=0.70011547_DP

    gbasis(45)%angmom(3)=1
    gbasis(45)%expo(3,1)=378.433426400_DP
    gbasis(45)%expo(3,2)=87.939789810_DP
    gbasis(45)%expo(3,3)=28.600748990_DP
    gbasis(45)%coco(3,1)=0.15591628_DP
    gbasis(45)%coco(3,2)=0.60768372_DP
    gbasis(45)%coco(3,3)=0.39195739_DP

    gbasis(45)%angmom(4)=0
    gbasis(45)%expo(4,1)=43.521794550_DP
    gbasis(45)%expo(4,2)=13.275534540_DP
    gbasis(45)%expo(4,3)=5.122113939_DP
    gbasis(45)%coco(4,1)=-0.22776350_DP
    gbasis(45)%coco(4,2)=0.21754360_DP
    gbasis(45)%coco(4,3)=0.91667696_DP

    gbasis(45)%angmom(5)=1
    gbasis(45)%expo(5,1)=43.521794550_DP
    gbasis(45)%expo(5,2)=13.275534540_DP
    gbasis(45)%expo(5,3)=5.122113939_DP
    gbasis(45)%coco(5,1)=0.00495151_DP
    gbasis(45)%coco(5,2)=0.57776647_DP
    gbasis(45)%coco(5,3)=0.48464604_DP

    gbasis(45)%angmom(6)=2
    gbasis(45)%expo(6,1)=43.521794550_DP
    gbasis(45)%expo(6,2)=13.275534540_DP
    gbasis(45)%expo(6,3)=5.122113939_DP
    gbasis(45)%coco(6,1)=0.21976795_DP
    gbasis(45)%coco(6,2)=0.65554736_DP
    gbasis(45)%coco(6,3)=0.28657326_DP

    gbasis(45)%angmom(7)=0
    gbasis(45)%expo(7,1)=4.540857408_DP
    gbasis(45)%expo(7,2)=1.768186338_DP
    gbasis(45)%expo(7,3)=0.778359979_DP
    gbasis(45)%coco(7,1)=-0.33061006_DP
    gbasis(45)%coco(7,2)=0.05761095_DP
    gbasis(45)%coco(7,3)=1.15578745_DP

    gbasis(45)%angmom(8)=1
    gbasis(45)%expo(8,1)=4.540857408_DP
    gbasis(45)%expo(8,2)=1.768186338_DP
    gbasis(45)%expo(8,3)=0.778359979_DP
    gbasis(45)%coco(8,1)=-0.12839276_DP
    gbasis(45)%coco(8,2)=0.58520476_DP
    gbasis(45)%coco(8,3)=0.54394420_DP

    gbasis(45)%angmom(9)=2
    gbasis(45)%expo(9,1)=2.779066094_DP
    gbasis(45)%expo(9,2)=1.082153932_DP
    gbasis(45)%expo(9,3)=0.476366825_DP
    gbasis(45)%coco(9,1)=0.12506621_DP
    gbasis(45)%coco(9,2)=0.66867856_DP
    gbasis(45)%coco(9,3)=0.30524682_DP

    gbasis(45)%angmom(10)=0
    gbasis(45)%expo(10,1)=0.413135485_DP
    gbasis(45)%expo(10,2)=0.222447917_DP
    gbasis(45)%expo(10,3)=0.098259157_DP
    gbasis(45)%coco(10,1)=-0.38426426_DP
    gbasis(45)%coco(10,2)=-0.19725674_DP
    gbasis(45)%coco(10,3)=1.37549551_DP

    gbasis(45)%angmom(11)=1
    gbasis(45)%expo(11,1)=0.413135485_DP
    gbasis(45)%expo(11,2)=0.222447917_DP
    gbasis(45)%expo(11,3)=0.098259157_DP
    gbasis(45)%coco(11,1)=-0.34816915_DP
    gbasis(45)%coco(11,2)=0.62903237_DP
    gbasis(45)%coco(11,3)=0.66628327_DP

    gbasis(45)%angmom(12)=3
    gbasis(45)%expo(12,1)=1.0949900000_DP
    gbasis(45)%coco(12,1)=1.0000000000_DP
    ! ====end Rh ==========================================


    ! ====begin Pd ==========================================
    gbasis(46)%angmom(1)=0
    gbasis(46)%expo(1,1)=4545.160269000_DP
    gbasis(46)%expo(1,2)=827.906616800_DP
    gbasis(46)%expo(1,3)=224.063840200_DP
    gbasis(46)%coco(1,1)=0.15432897_DP
    gbasis(46)%coco(1,2)=0.53532814_DP
    gbasis(46)%coco(1,3)=0.44463454_DP

    gbasis(46)%angmom(2)=0
    gbasis(46)%expo(2,1)=396.488943300_DP
    gbasis(46)%expo(2,2)=92.135503650_DP
    gbasis(46)%expo(2,3)=29.965325350_DP
    gbasis(46)%coco(2,1)=-0.09996723_DP
    gbasis(46)%coco(2,2)=0.39951283_DP
    gbasis(46)%coco(2,3)=0.70011547_DP

    gbasis(46)%angmom(3)=1
    gbasis(46)%expo(3,1)=396.488943300_DP
    gbasis(46)%expo(3,2)=92.135503650_DP
    gbasis(46)%expo(3,3)=29.965325350_DP
    gbasis(46)%coco(3,1)=0.15591628_DP
    gbasis(46)%coco(3,2)=0.60768372_DP
    gbasis(46)%coco(3,3)=0.39195739_DP

    gbasis(46)%angmom(4)=0
    gbasis(46)%expo(4,1)=46.419450970_DP
    gbasis(46)%expo(4,2)=14.159412110_DP
    gbasis(46)%expo(4,3)=5.463141383_DP
    gbasis(46)%coco(4,1)=-0.22776350_DP
    gbasis(46)%coco(4,2)=0.21754360_DP
    gbasis(46)%coco(4,3)=0.91667696_DP

    gbasis(46)%angmom(5)=1
    gbasis(46)%expo(5,1)=46.419450970_DP
    gbasis(46)%expo(5,2)=14.159412110_DP
    gbasis(46)%expo(5,3)=5.463141383_DP
    gbasis(46)%coco(5,1)=0.00495151_DP
    gbasis(46)%coco(5,2)=0.57776647_DP
    gbasis(46)%coco(5,3)=0.48464604_DP

    gbasis(46)%angmom(6)=2
    gbasis(46)%expo(6,1)=46.419450970_DP
    gbasis(46)%expo(6,2)=14.159412110_DP
    gbasis(46)%expo(6,3)=5.463141383_DP
    gbasis(46)%coco(6,1)=0.21976795_DP
    gbasis(46)%coco(6,2)=0.65554736_DP
    gbasis(46)%coco(6,3)=0.28657326_DP

    gbasis(46)%angmom(7)=0
    gbasis(46)%expo(7,1)=4.919104589_DP
    gbasis(46)%expo(7,2)=1.915473830_DP
    gbasis(46)%expo(7,3)=0.843196295_DP
    gbasis(46)%coco(7,1)=-0.33061006_DP
    gbasis(46)%coco(7,2)=0.05761095_DP
    gbasis(46)%coco(7,3)=1.15578745_DP

    gbasis(46)%angmom(8)=1
    gbasis(46)%expo(8,1)=4.919104589_DP
    gbasis(46)%expo(8,2)=1.915473830_DP
    gbasis(46)%expo(8,3)=0.843196295_DP
    gbasis(46)%coco(8,1)=-0.12839276_DP
    gbasis(46)%coco(8,2)=0.58520476_DP
    gbasis(46)%coco(8,3)=0.54394420_DP

    gbasis(46)%angmom(9)=2
    gbasis(46)%expo(9,1)=3.025977448_DP
    gbasis(46)%expo(9,2)=1.178299934_DP
    gbasis(46)%expo(9,3)=0.518690532_DP
    gbasis(46)%coco(9,1)=0.12506621_DP
    gbasis(46)%coco(9,2)=0.66867856_DP
    gbasis(46)%coco(9,3)=0.30524682_DP

    gbasis(46)%angmom(10)=0
    gbasis(46)%expo(10,1)=0.437080480_DP
    gbasis(46)%expo(10,2)=0.235340816_DP
    gbasis(46)%expo(10,3)=0.103954177_DP
    gbasis(46)%coco(10,1)=-0.38426426_DP
    gbasis(46)%coco(10,2)=-0.19725674_DP
    gbasis(46)%coco(10,3)=1.37549551_DP

    gbasis(46)%angmom(11)=1
    gbasis(46)%expo(11,1)=0.437080480_DP
    gbasis(46)%expo(11,2)=0.235340816_DP
    gbasis(46)%expo(11,3)=0.103954177_DP
    gbasis(46)%coco(11,1)=-0.34816915_DP
    gbasis(46)%coco(11,2)=0.62903237_DP
    gbasis(46)%coco(11,3)=0.66628327_DP

    gbasis(46)%angmom(12)=3
    gbasis(46)%expo(12,1)=1.2462900000_DP
    gbasis(46)%coco(12,1)=1.0000000000_DP
    ! ====end Pd ==========================================


    ! ====begin Ag ==========================================
    gbasis(47)%angmom(1)=0
    gbasis(47)%expo(1,1)=4744.521634000_DP
    gbasis(47)%expo(1,2)=864.220538300_DP
    gbasis(47)%expo(1,3)=233.891804500_DP
    gbasis(47)%coco(1,1)=0.15432897_DP
    gbasis(47)%coco(1,2)=0.53532814_DP
    gbasis(47)%coco(1,3)=0.44463454_DP

    gbasis(47)%angmom(2)=0
    gbasis(47)%expo(2,1)=414.965206900_DP
    gbasis(47)%expo(2,2)=96.428989950_DP
    gbasis(47)%expo(2,3)=31.361700350_DP
    gbasis(47)%coco(2,1)=-0.09996723_DP
    gbasis(47)%coco(2,2)=0.39951283_DP
    gbasis(47)%coco(2,3)=0.70011547_DP

    gbasis(47)%angmom(3)=1
    gbasis(47)%expo(3,1)=414.965206900_DP
    gbasis(47)%expo(3,2)=96.428989950_DP
    gbasis(47)%expo(3,3)=31.361700350_DP
    gbasis(47)%coco(3,1)=0.15591628_DP
    gbasis(47)%coco(3,2)=0.60768372_DP
    gbasis(47)%coco(3,3)=0.39195739_DP

    gbasis(47)%angmom(4)=0
    gbasis(47)%expo(4,1)=49.410486050_DP
    gbasis(47)%expo(4,2)=15.071773140_DP
    gbasis(47)%expo(4,3)=5.815158634_DP
    gbasis(47)%coco(4,1)=-0.22776350_DP
    gbasis(47)%coco(4,2)=0.21754360_DP
    gbasis(47)%coco(4,3)=0.91667696_DP

    gbasis(47)%angmom(5)=1
    gbasis(47)%expo(5,1)=49.410486050_DP
    gbasis(47)%expo(5,2)=15.071773140_DP
    gbasis(47)%expo(5,3)=5.815158634_DP
    gbasis(47)%coco(5,1)=0.00495151_DP
    gbasis(47)%coco(5,2)=0.57776647_DP
    gbasis(47)%coco(5,3)=0.48464604_DP

    gbasis(47)%angmom(6)=2
    gbasis(47)%expo(6,1)=49.410486050_DP
    gbasis(47)%expo(6,2)=15.071773140_DP
    gbasis(47)%expo(6,3)=5.815158634_DP
    gbasis(47)%coco(6,1)=0.21976795_DP
    gbasis(47)%coco(6,2)=0.65554736_DP
    gbasis(47)%coco(6,3)=0.28657326_DP

    gbasis(47)%angmom(7)=0
    gbasis(47)%expo(7,1)=5.290230450_DP
    gbasis(47)%expo(7,2)=2.059988316_DP
    gbasis(47)%expo(7,3)=0.906811928_DP
    gbasis(47)%coco(7,1)=-0.33061006_DP
    gbasis(47)%coco(7,2)=0.05761095_DP
    gbasis(47)%coco(7,3)=1.15578745_DP

    gbasis(47)%angmom(8)=1
    gbasis(47)%expo(8,1)=5.290230450_DP
    gbasis(47)%expo(8,2)=2.059988316_DP
    gbasis(47)%expo(8,3)=0.906811928_DP
    gbasis(47)%coco(8,1)=-0.12839276_DP
    gbasis(47)%coco(8,2)=0.58520476_DP
    gbasis(47)%coco(8,3)=0.54394420_DP

    gbasis(47)%angmom(9)=2
    gbasis(47)%expo(9,1)=3.283395668_DP
    gbasis(47)%expo(9,2)=1.278537254_DP
    gbasis(47)%expo(9,3)=0.562815247_DP
    gbasis(47)%coco(9,1)=0.12506621_DP
    gbasis(47)%coco(9,2)=0.66867856_DP
    gbasis(47)%coco(9,3)=0.30524682_DP

    gbasis(47)%angmom(10)=0
    gbasis(47)%expo(10,1)=0.437080480_DP
    gbasis(47)%expo(10,2)=0.235340816_DP
    gbasis(47)%expo(10,3)=0.103954177_DP
    gbasis(47)%coco(10,1)=-0.38426426_DP
    gbasis(47)%coco(10,2)=-0.19725674_DP
    gbasis(47)%coco(10,3)=1.37549551_DP

    gbasis(47)%angmom(11)=1
    gbasis(47)%expo(11,1)=0.437080480_DP
    gbasis(47)%expo(11,2)=0.235340816_DP
    gbasis(47)%expo(11,3)=0.103954177_DP
    gbasis(47)%coco(11,1)=-0.34816915_DP
    gbasis(47)%coco(11,2)=0.62903237_DP
    gbasis(47)%coco(11,3)=0.66628327_DP

    gbasis(47)%angmom(12)=3
    gbasis(47)%expo(12,1)=1.3971100000_DP
    gbasis(47)%coco(12,1)=1.0000000000_DP
    ! ====end Ag ==========================================


    ! ====begin Cd ==========================================
    gbasis(48)%angmom(1)=0
    gbasis(48)%expo(1,1)=4950.261905000_DP
    gbasis(48)%expo(1,2)=901.696385600_DP
    gbasis(48)%expo(1,3)=244.034231300_DP
    gbasis(48)%coco(1,1)=0.15432897_DP
    gbasis(48)%coco(1,2)=0.53532814_DP
    gbasis(48)%coco(1,3)=0.44463454_DP

    gbasis(48)%angmom(2)=0
    gbasis(48)%expo(2,1)=433.446938500_DP
    gbasis(48)%expo(2,2)=100.723746900_DP
    gbasis(48)%expo(2,3)=32.758488610_DP
    gbasis(48)%coco(2,1)=-0.09996723_DP
    gbasis(48)%coco(2,2)=0.39951283_DP
    gbasis(48)%coco(2,3)=0.70011547_DP

    gbasis(48)%angmom(3)=1
    gbasis(48)%expo(3,1)=433.446938500_DP
    gbasis(48)%expo(3,2)=100.723746900_DP
    gbasis(48)%expo(3,3)=32.758488610_DP
    gbasis(48)%coco(3,1)=0.15591628_DP
    gbasis(48)%coco(3,2)=0.60768372_DP
    gbasis(48)%coco(3,3)=0.39195739_DP

    gbasis(48)%angmom(4)=0
    gbasis(48)%expo(4,1)=52.592792350_DP
    gbasis(48)%expo(4,2)=16.042478000_DP
    gbasis(48)%expo(4,3)=6.189686744_DP
    gbasis(48)%coco(4,1)=-0.22776350_DP
    gbasis(48)%coco(4,2)=0.21754360_DP
    gbasis(48)%coco(4,3)=0.91667696_DP

    gbasis(48)%angmom(5)=1
    gbasis(48)%expo(5,1)=52.592792350_DP
    gbasis(48)%expo(5,2)=16.042478000_DP
    gbasis(48)%expo(5,3)=6.189686744_DP
    gbasis(48)%coco(5,1)=0.00495151_DP
    gbasis(48)%coco(5,2)=0.57776647_DP
    gbasis(48)%coco(5,3)=0.48464604_DP

    gbasis(48)%angmom(6)=2
    gbasis(48)%expo(6,1)=52.592792350_DP
    gbasis(48)%expo(6,2)=16.042478000_DP
    gbasis(48)%expo(6,3)=6.189686744_DP
    gbasis(48)%coco(6,1)=0.21976795_DP
    gbasis(48)%coco(6,2)=0.65554736_DP
    gbasis(48)%coco(6,3)=0.28657326_DP

    gbasis(48)%angmom(7)=0
    gbasis(48)%expo(7,1)=5.674851796_DP
    gbasis(48)%expo(7,2)=2.209757875_DP
    gbasis(48)%expo(7,3)=0.972740857_DP
    gbasis(48)%coco(7,1)=-0.33061006_DP
    gbasis(48)%coco(7,2)=0.05761095_DP
    gbasis(48)%coco(7,3)=1.15578745_DP

    gbasis(48)%angmom(8)=1
    gbasis(48)%expo(8,1)=5.674851796_DP
    gbasis(48)%expo(8,2)=2.209757875_DP
    gbasis(48)%expo(8,3)=0.972740857_DP
    gbasis(48)%coco(8,1)=-0.12839276_DP
    gbasis(48)%coco(8,2)=0.58520476_DP
    gbasis(48)%coco(8,3)=0.54394420_DP

    gbasis(48)%angmom(9)=2
    gbasis(48)%expo(9,1)=3.642963976_DP
    gbasis(48)%expo(9,2)=1.418551290_DP
    gbasis(48)%expo(9,3)=0.624449770_DP
    gbasis(48)%coco(9,1)=0.12506621_DP
    gbasis(48)%coco(9,2)=0.66867856_DP
    gbasis(48)%coco(9,3)=0.30524682_DP

    gbasis(48)%angmom(10)=0
    gbasis(48)%expo(10,1)=0.594915098_DP
    gbasis(48)%expo(10,2)=0.320325000_DP
    gbasis(48)%expo(10,3)=0.141493186_DP
    gbasis(48)%coco(10,1)=-0.38426426_DP
    gbasis(48)%coco(10,2)=-0.19725674_DP
    gbasis(48)%coco(10,3)=1.37549551_DP

    gbasis(48)%angmom(11)=1
    gbasis(48)%expo(11,1)=0.594915098_DP
    gbasis(48)%expo(11,2)=0.320325000_DP
    gbasis(48)%expo(11,3)=0.141493186_DP
    gbasis(48)%coco(11,1)=-0.34816915_DP
    gbasis(48)%coco(11,2)=0.62903237_DP
    gbasis(48)%coco(11,3)=0.66628327_DP

    gbasis(48)%angmom(12)=3
    gbasis(48)%expo(12,1)=1.5981300000_DP
    gbasis(48)%coco(12,1)=1.0000000000_DP
    ! ====end Cd ==========================================


    ! ====begin In ==========================================
    gbasis(49)%angmom(1)=0
    gbasis(49)%expo(1,1)=5158.224714000_DP
    gbasis(49)%expo(1,2)=939.577070700_DP
    gbasis(49)%expo(1,3)=254.286223100_DP
    gbasis(49)%coco(1,1)=0.15432897_DP
    gbasis(49)%coco(1,2)=0.53532814_DP
    gbasis(49)%coco(1,3)=0.44463454_DP

    gbasis(49)%angmom(2)=0
    gbasis(49)%expo(2,1)=452.331322300_DP
    gbasis(49)%expo(2,2)=105.112071600_DP
    gbasis(49)%expo(2,3)=34.185707990_DP
    gbasis(49)%coco(2,1)=-0.09996723_DP
    gbasis(49)%coco(2,2)=0.39951283_DP
    gbasis(49)%coco(2,3)=0.70011547_DP

    gbasis(49)%angmom(3)=1
    gbasis(49)%expo(3,1)=452.331322300_DP
    gbasis(49)%expo(3,2)=105.112071600_DP
    gbasis(49)%expo(3,3)=34.185707990_DP
    gbasis(49)%coco(3,1)=0.15591628_DP
    gbasis(49)%coco(3,2)=0.60768372_DP
    gbasis(49)%coco(3,3)=0.39195739_DP

    gbasis(49)%angmom(4)=0
    gbasis(49)%expo(4,1)=55.975397690_DP
    gbasis(49)%expo(4,2)=17.074280440_DP
    gbasis(49)%expo(4,3)=6.587788204_DP
    gbasis(49)%coco(4,1)=-0.22776350_DP
    gbasis(49)%coco(4,2)=0.21754360_DP
    gbasis(49)%coco(4,3)=0.91667696_DP

    gbasis(49)%angmom(5)=1
    gbasis(49)%expo(5,1)=55.975397690_DP
    gbasis(49)%expo(5,2)=17.074280440_DP
    gbasis(49)%expo(5,3)=6.587788204_DP
    gbasis(49)%coco(5,1)=0.00495151_DP
    gbasis(49)%coco(5,2)=0.57776647_DP
    gbasis(49)%coco(5,3)=0.48464604_DP

    gbasis(49)%angmom(6)=2
    gbasis(49)%expo(6,1)=55.975397690_DP
    gbasis(49)%expo(6,2)=17.074280440_DP
    gbasis(49)%expo(6,3)=6.587788204_DP
    gbasis(49)%coco(6,1)=0.21976795_DP
    gbasis(49)%coco(6,2)=0.65554736_DP
    gbasis(49)%coco(6,3)=0.28657326_DP

    gbasis(49)%angmom(7)=0
    gbasis(49)%expo(7,1)=5.048549180_DP
    gbasis(49)%expo(7,2)=1.965878882_DP
    gbasis(49)%expo(7,3)=0.865384724_DP
    gbasis(49)%coco(7,1)=-0.33061006_DP
    gbasis(49)%coco(7,2)=0.05761095_DP
    gbasis(49)%coco(7,3)=1.11557874_DP

    gbasis(49)%angmom(8)=1
    gbasis(49)%expo(8,1)=5.048549180_DP
    gbasis(49)%expo(8,2)=1.965878882_DP
    gbasis(49)%expo(8,3)=0.865384724_DP
    gbasis(49)%coco(8,1)=-0.12839276_DP
    gbasis(49)%coco(8,2)=0.58520476_DP
    gbasis(49)%coco(8,3)=0.54394420_DP

    gbasis(49)%angmom(9)=2
    gbasis(49)%expo(9,1)=5.048549180_DP
    gbasis(49)%expo(9,2)=1.965878882_DP
    gbasis(49)%expo(9,3)=0.865384724_DP
    gbasis(49)%coco(9,1)=0.12506621_DP
    gbasis(49)%coco(9,2)=0.66867856_DP
    gbasis(49)%coco(9,3)=0.30524682_DP

    gbasis(49)%angmom(10)=0
    gbasis(49)%expo(10,1)=0.566923061_DP
    gbasis(49)%expo(10,2)=0.305253019_DP
    gbasis(49)%expo(10,3)=0.134835626_DP
    gbasis(49)%coco(10,1)=-0.38426426_DP
    gbasis(49)%coco(10,2)=-0.19725674_DP
    gbasis(49)%coco(10,3)=1.37549551_DP

    gbasis(49)%angmom(11)=1
    gbasis(49)%expo(11,1)=0.566923061_DP
    gbasis(49)%expo(11,2)=0.305253019_DP
    gbasis(49)%expo(11,3)=0.134835626_DP
    gbasis(49)%coco(11,1)=-0.34816915_DP
    gbasis(49)%coco(11,2)=0.62903237_DP
    gbasis(49)%coco(11,3)=0.66628327_DP

    gbasis(49)%angmom(12)=2
    gbasis(49)%expo(12,1)=0.18000000000_DP
    gbasis(49)%coco(12,1)=1.0000000000_DP
    ! ====end In ==========================================


    ! ====begin Sn ==========================================
    gbasis(50)%angmom(1)=0
    gbasis(50)%expo(1,1)=5370.466413000_DP
    gbasis(50)%expo(1,2)=978.237161100_DP
    gbasis(50)%expo(1,3)=264.749152200_DP
    gbasis(50)%coco(1,1)=0.15432897_DP
    gbasis(50)%coco(1,2)=0.53532814_DP
    gbasis(50)%coco(1,3)=0.44463454_DP

    gbasis(50)%angmom(2)=0
    gbasis(50)%expo(2,1)=472.051532200_DP
    gbasis(50)%expo(2,2)=109.694624300_DP
    gbasis(50)%expo(2,3)=35.676096360_DP
    gbasis(50)%coco(2,1)=-0.09996723_DP
    gbasis(50)%coco(2,2)=0.39951283_DP
    gbasis(50)%coco(2,3)=0.70011547_DP

    gbasis(50)%angmom(3)=1
    gbasis(50)%expo(3,1)=472.051532200_DP
    gbasis(50)%expo(3,2)=109.694624300_DP
    gbasis(50)%expo(3,3)=35.676096360_DP
    gbasis(50)%coco(3,1)=0.15591628_DP
    gbasis(50)%coco(3,2)=0.60768372_DP
    gbasis(50)%coco(3,3)=0.39195739_DP

    gbasis(50)%angmom(4)=0
    gbasis(50)%expo(4,1)=59.151411880_DP
    gbasis(50)%expo(4,2)=18.043066000_DP
    gbasis(50)%expo(4,3)=6.961575790_DP
    gbasis(50)%coco(4,1)=-0.22776350_DP
    gbasis(50)%coco(4,2)=0.21754360_DP
    gbasis(50)%coco(4,3)=0.91667696_DP

    gbasis(50)%angmom(5)=1
    gbasis(50)%expo(5,1)=59.151411880_DP
    gbasis(50)%expo(5,2)=18.043066000_DP
    gbasis(50)%expo(5,3)=6.961575790_DP
    gbasis(50)%coco(5,1)=0.00495151_DP
    gbasis(50)%coco(5,2)=0.57776647_DP
    gbasis(50)%coco(5,3)=0.48464604_DP

    gbasis(50)%angmom(6)=2
    gbasis(50)%expo(6,1)=59.151411880_DP
    gbasis(50)%expo(6,2)=18.043066000_DP
    gbasis(50)%expo(6,3)=6.961575790_DP
    gbasis(50)%coco(6,1)=0.21976795_DP
    gbasis(50)%coco(6,2)=0.65554736_DP
    gbasis(50)%coco(6,3)=0.28657326_DP

    gbasis(50)%angmom(7)=0
    gbasis(50)%expo(7,1)=5.583138529_DP
    gbasis(50)%expo(7,2)=2.174045204_DP
    gbasis(50)%expo(7,3)=0.957020051_DP
    gbasis(50)%coco(7,1)=-0.33061006_DP
    gbasis(50)%coco(7,2)=0.05761095_DP
    gbasis(50)%coco(7,3)=1.11557874_DP

    gbasis(50)%angmom(8)=1
    gbasis(50)%expo(8,1)=5.583138529_DP
    gbasis(50)%expo(8,2)=2.174045204_DP
    gbasis(50)%expo(8,3)=0.957020051_DP
    gbasis(50)%coco(8,1)=-0.12839276_DP
    gbasis(50)%coco(8,2)=0.58520476_DP
    gbasis(50)%coco(8,3)=0.54394420_DP

    gbasis(50)%angmom(9)=2
    gbasis(50)%expo(9,1)=5.583138529_DP
    gbasis(50)%expo(9,2)=2.174045204_DP
    gbasis(50)%expo(9,3)=0.957020051_DP
    gbasis(50)%coco(9,1)=0.12506621_DP
    gbasis(50)%coco(9,2)=0.66867856_DP
    gbasis(50)%coco(9,3)=0.30524682_DP

    gbasis(50)%angmom(10)=0
    gbasis(50)%expo(10,1)=0.623581642_DP
    gbasis(50)%expo(10,2)=0.335760162_DP
    gbasis(50)%expo(10,3)=0.148311168_DP
    gbasis(50)%coco(10,1)=-0.38426426_DP
    gbasis(50)%coco(10,2)=-0.19725674_DP
    gbasis(50)%coco(10,3)=1.37549551_DP

    gbasis(50)%angmom(11)=1
    gbasis(50)%expo(11,1)=0.623581642_DP
    gbasis(50)%expo(11,2)=0.335760162_DP
    gbasis(50)%expo(11,3)=0.148311168_DP
    gbasis(50)%coco(11,1)=-0.34816915_DP
    gbasis(50)%coco(11,2)=0.62903237_DP
    gbasis(50)%coco(11,3)=0.66628327_DP

    gbasis(50)%angmom(12)=2
    gbasis(50)%expo(12,1)=0.20500000000_DP
    gbasis(50)%coco(12,1)=1.0000000000_DP
    ! ====end Sn ==========================================


    ! ====begin Sb ==========================================
    gbasis(51)%angmom(1)=0
    gbasis(51)%expo(1,1)=5586.987002000_DP
    gbasis(51)%expo(1,2)=1017.676657000_DP
    gbasis(51)%expo(1,3)=275.423018900_DP
    gbasis(51)%coco(1,1)=0.15432897_DP
    gbasis(51)%coco(1,2)=0.53532814_DP
    gbasis(51)%coco(1,3)=0.44463454_DP

    gbasis(51)%angmom(2)=0
    gbasis(51)%expo(2,1)=492.192488800_DP
    gbasis(51)%expo(2,2)=114.374949400_DP
    gbasis(51)%expo(2,3)=37.198283360_DP
    gbasis(51)%coco(2,1)=-0.09996723_DP
    gbasis(51)%coco(2,2)=0.39951283_DP
    gbasis(51)%coco(2,3)=0.70011547_DP

    gbasis(51)%angmom(3)=1
    gbasis(51)%expo(3,1)=492.192488800_DP
    gbasis(51)%expo(3,2)=114.374949400_DP
    gbasis(51)%expo(3,3)=37.198283360_DP
    gbasis(51)%coco(3,1)=0.15591628_DP
    gbasis(51)%coco(3,2)=0.60768372_DP
    gbasis(51)%coco(3,3)=0.39195739_DP

    gbasis(51)%angmom(4)=0
    gbasis(51)%expo(4,1)=62.521797750_DP
    gbasis(51)%expo(4,2)=19.071141120_DP
    gbasis(51)%expo(4,3)=7.358239131_DP
    gbasis(51)%coco(4,1)=-0.22776350_DP
    gbasis(51)%coco(4,2)=0.21754360_DP
    gbasis(51)%coco(4,3)=0.91667696_DP

    gbasis(51)%angmom(5)=1
    gbasis(51)%expo(5,1)=62.521797750_DP
    gbasis(51)%expo(5,2)=19.071141120_DP
    gbasis(51)%expo(5,3)=7.358239131_DP
    gbasis(51)%coco(5,1)=0.00495151_DP
    gbasis(51)%coco(5,2)=0.57776647_DP
    gbasis(51)%coco(5,3)=0.48464604_DP

    gbasis(51)%angmom(6)=2
    gbasis(51)%expo(6,1)=62.521797750_DP
    gbasis(51)%expo(6,2)=19.071141120_DP
    gbasis(51)%expo(6,3)=7.358239131_DP
    gbasis(51)%coco(6,1)=0.21976795_DP
    gbasis(51)%coco(6,2)=0.65554736_DP
    gbasis(51)%coco(6,3)=0.28657326_DP

    gbasis(51)%angmom(7)=0
    gbasis(51)%expo(7,1)=6.120693149_DP
    gbasis(51)%expo(7,2)=2.383366187_DP
    gbasis(51)%expo(7,3)=1.049163663_DP
    gbasis(51)%coco(7,1)=-0.33061006_DP
    gbasis(51)%coco(7,2)=0.05761095_DP
    gbasis(51)%coco(7,3)=1.11557874_DP

    gbasis(51)%angmom(8)=1
    gbasis(51)%expo(8,1)=6.120693149_DP
    gbasis(51)%expo(8,2)=2.383366187_DP
    gbasis(51)%expo(8,3)=1.049163663_DP
    gbasis(51)%coco(8,1)=-0.12839276_DP
    gbasis(51)%coco(8,2)=0.58520476_DP
    gbasis(51)%coco(8,3)=0.54394420_DP

    gbasis(51)%angmom(9)=2
    gbasis(51)%expo(9,1)=6.120693149_DP
    gbasis(51)%expo(9,2)=2.383366187_DP
    gbasis(51)%expo(9,3)=1.049163663_DP
    gbasis(51)%coco(9,1)=0.12506621_DP
    gbasis(51)%coco(9,2)=0.66867856_DP
    gbasis(51)%coco(9,3)=0.30524682_DP

    gbasis(51)%angmom(10)=0
    gbasis(51)%expo(10,1)=0.652922693_DP
    gbasis(51)%expo(10,2)=0.351558503_DP
    gbasis(51)%expo(10,3)=0.155289573_DP
    gbasis(51)%coco(10,1)=-0.38426426_DP
    gbasis(51)%coco(10,2)=-0.19725674_DP
    gbasis(51)%coco(10,3)=1.37549551_DP

    gbasis(51)%angmom(11)=1
    gbasis(51)%expo(11,1)=0.652922693_DP
    gbasis(51)%expo(11,2)=0.351558503_DP
    gbasis(51)%expo(11,3)=0.155289573_DP
    gbasis(51)%coco(11,1)=-0.34816915_DP
    gbasis(51)%coco(11,2)=0.62903237_DP
    gbasis(51)%coco(11,3)=0.66628327_DP

    gbasis(51)%angmom(12)=2
    gbasis(51)%expo(12,1)=0.23060000000_DP
    gbasis(51)%coco(12,1)=1.0000000000_DP
    ! ====end Sb ==========================================


    ! ====begin Te ==========================================
    gbasis(52)%angmom(1)=0
    gbasis(52)%expo(1,1)=5810.061591000_DP
    gbasis(52)%expo(1,2)=1058.309972000_DP
    gbasis(52)%expo(1,3)=286.419979700_DP
    gbasis(52)%coco(1,1)=0.15432897_DP
    gbasis(52)%coco(1,2)=0.53532814_DP
    gbasis(52)%coco(1,3)=0.44463454_DP

    gbasis(52)%angmom(2)=0
    gbasis(52)%expo(2,1)=512.754192000_DP
    gbasis(52)%expo(2,2)=119.153047100_DP
    gbasis(52)%expo(2,3)=38.752269000_DP
    gbasis(52)%coco(2,1)=-0.09996723_DP
    gbasis(52)%coco(2,2)=0.39951283_DP
    gbasis(52)%coco(2,3)=0.70011547_DP

    gbasis(52)%angmom(3)=1
    gbasis(52)%expo(3,1)=512.754192000_DP
    gbasis(52)%expo(3,2)=119.153047100_DP
    gbasis(52)%expo(3,3)=38.752269000_DP
    gbasis(52)%coco(3,1)=0.15591628_DP
    gbasis(52)%coco(3,2)=0.60768372_DP
    gbasis(52)%coco(3,3)=0.39195739_DP

    gbasis(52)%angmom(4)=0
    gbasis(52)%expo(4,1)=65.985562270_DP
    gbasis(52)%expo(4,2)=20.127699700_DP
    gbasis(52)%expo(4,3)=7.765892279_DP
    gbasis(52)%coco(4,1)=-0.22776350_DP
    gbasis(52)%coco(4,2)=0.21754360_DP
    gbasis(52)%coco(4,3)=0.91667696_DP

    gbasis(52)%angmom(5)=1
    gbasis(52)%expo(5,1)=65.985562270_DP
    gbasis(52)%expo(5,2)=20.127699700_DP
    gbasis(52)%expo(5,3)=7.765892279_DP
    gbasis(52)%coco(5,1)=0.00495151_DP
    gbasis(52)%coco(5,2)=0.57776647_DP
    gbasis(52)%coco(5,3)=0.48464604_DP

    gbasis(52)%angmom(6)=2
    gbasis(52)%expo(6,1)=65.985562270_DP
    gbasis(52)%expo(6,2)=20.127699700_DP
    gbasis(52)%expo(6,3)=7.765892279_DP
    gbasis(52)%coco(6,1)=0.21976795_DP
    gbasis(52)%coco(6,2)=0.65554736_DP
    gbasis(52)%coco(6,3)=0.28657326_DP

    gbasis(52)%angmom(7)=0
    gbasis(52)%expo(7,1)=6.707956921_DP
    gbasis(52)%expo(7,2)=2.612043655_DP
    gbasis(52)%expo(7,3)=1.149828048_DP
    gbasis(52)%coco(7,1)=-0.33061006_DP
    gbasis(52)%coco(7,2)=0.05761095_DP
    gbasis(52)%coco(7,3)=1.11557874_DP

    gbasis(52)%angmom(8)=1
    gbasis(52)%expo(8,1)=6.707956921_DP
    gbasis(52)%expo(8,2)=2.612043655_DP
    gbasis(52)%expo(8,3)=1.149828048_DP
    gbasis(52)%coco(8,1)=-0.12839276_DP
    gbasis(52)%coco(8,2)=0.58520476_DP
    gbasis(52)%coco(8,3)=0.54394420_DP

    gbasis(52)%angmom(9)=2
    gbasis(52)%expo(9,1)=6.707956921_DP
    gbasis(52)%expo(9,2)=2.612043655_DP
    gbasis(52)%expo(9,3)=1.149828048_DP
    gbasis(52)%coco(9,1)=0.12506621_DP
    gbasis(52)%coco(9,2)=0.66867856_DP
    gbasis(52)%coco(9,3)=0.30524682_DP

    gbasis(52)%angmom(10)=0
    gbasis(52)%expo(10,1)=0.701271348_DP
    gbasis(52)%expo(10,2)=0.377591265_DP
    gbasis(52)%expo(10,3)=0.166788702_DP
    gbasis(52)%coco(10,1)=-0.38426426_DP
    gbasis(52)%coco(10,2)=-0.19725674_DP
    gbasis(52)%coco(10,3)=1.37549551_DP

    gbasis(52)%angmom(11)=1
    gbasis(52)%expo(11,1)=0.701271348_DP
    gbasis(52)%expo(11,2)=0.377591265_DP
    gbasis(52)%expo(11,3)=0.166788702_DP
    gbasis(52)%coco(11,1)=-0.34816915_DP
    gbasis(52)%coco(11,2)=0.62903237_DP
    gbasis(52)%coco(11,3)=0.66628327_DP

    gbasis(52)%angmom(12)=2
    gbasis(52)%expo(12,1)=0.25000000000_DP
    gbasis(52)%coco(12,1)=1.0000000000_DP
    ! ====end Te ==========================================


    ! ====begin I ==========================================
    gbasis(53)%angmom(1)=0
    gbasis(53)%expo(1,1)=6035.183623000_DP
    gbasis(53)%expo(1,2)=1099.316231000_DP
    gbasis(53)%expo(1,3)=297.517873700_DP
    gbasis(53)%coco(1,1)=0.15432897_DP
    gbasis(53)%coco(1,2)=0.53532814_DP
    gbasis(53)%coco(1,3)=0.44463454_DP

    gbasis(53)%angmom(2)=0
    gbasis(53)%expo(2,1)=533.736641800_DP
    gbasis(53)%expo(2,2)=124.028917100_DP
    gbasis(53)%expo(2,3)=40.338053280_DP
    gbasis(53)%coco(2,1)=-0.09996723_DP
    gbasis(53)%coco(2,2)=0.39951283_DP
    gbasis(53)%coco(2,3)=0.70011547_DP

    gbasis(53)%angmom(3)=1
    gbasis(53)%expo(3,1)=533.736641800_DP
    gbasis(53)%expo(3,2)=124.028917100_DP
    gbasis(53)%expo(3,3)=40.338053280_DP
    gbasis(53)%coco(3,1)=0.15591628_DP
    gbasis(53)%coco(3,2)=0.60768372_DP
    gbasis(53)%coco(3,3)=0.39195739_DP

    gbasis(53)%angmom(4)=0
    gbasis(53)%expo(4,1)=69.542705450_DP
    gbasis(53)%expo(4,2)=21.212741750_DP
    gbasis(53)%expo(4,3)=8.184535234_DP
    gbasis(53)%coco(4,1)=-0.22776350_DP
    gbasis(53)%coco(4,2)=0.21754360_DP
    gbasis(53)%coco(4,3)=0.91667696_DP

    gbasis(53)%angmom(5)=1
    gbasis(53)%expo(5,1)=69.542705450_DP
    gbasis(53)%expo(5,2)=21.212741750_DP
    gbasis(53)%expo(5,3)=8.184535234_DP
    gbasis(53)%coco(5,1)=0.00495151_DP
    gbasis(53)%coco(5,2)=0.57776647_DP
    gbasis(53)%coco(5,3)=0.48464604_DP

    gbasis(53)%angmom(6)=2
    gbasis(53)%expo(6,1)=69.542705450_DP
    gbasis(53)%expo(6,2)=21.212741750_DP
    gbasis(53)%expo(6,3)=8.184535234_DP
    gbasis(53)%coco(6,1)=0.21976795_DP
    gbasis(53)%coco(6,2)=0.65554736_DP
    gbasis(53)%coco(6,3)=0.28657326_DP

    gbasis(53)%angmom(7)=0
    gbasis(53)%expo(7,1)=7.295991196_DP
    gbasis(53)%expo(7,2)=2.841021154_DP
    gbasis(53)%expo(7,3)=1.250624506_DP
    gbasis(53)%coco(7,1)=-0.33061006_DP
    gbasis(53)%coco(7,2)=0.05761095_DP
    gbasis(53)%coco(7,3)=1.11557874_DP

    gbasis(53)%angmom(8)=1
    gbasis(53)%expo(8,1)=7.295991196_DP
    gbasis(53)%expo(8,2)=2.841021154_DP
    gbasis(53)%expo(8,3)=1.250624506_DP
    gbasis(53)%coco(8,1)=-0.12839276_DP
    gbasis(53)%coco(8,2)=0.58520476_DP
    gbasis(53)%coco(8,3)=0.54394420_DP

    gbasis(53)%angmom(9)=2
    gbasis(53)%expo(9,1)=7.295991196_DP
    gbasis(53)%expo(9,2)=2.841021154_DP
    gbasis(53)%expo(9,3)=1.250624506_DP
    gbasis(53)%coco(9,1)=0.12506621_DP
    gbasis(53)%coco(9,2)=0.66867856_DP
    gbasis(53)%coco(9,3)=0.30524682_DP

    gbasis(53)%angmom(10)=0
    gbasis(53)%expo(10,1)=0.790036458_DP
    gbasis(53)%expo(10,2)=0.425385789_DP
    gbasis(53)%expo(10,3)=0.187900384_DP
    gbasis(53)%coco(10,1)=-0.38426426_DP
    gbasis(53)%coco(10,2)=-0.19725674_DP
    gbasis(53)%coco(10,3)=1.37549551_DP

    gbasis(53)%angmom(11)=1
    gbasis(53)%expo(11,1)=0.790036458_DP
    gbasis(53)%expo(11,2)=0.425385789_DP
    gbasis(53)%expo(11,3)=0.187900384_DP
    gbasis(53)%coco(11,1)=-0.34816915_DP
    gbasis(53)%coco(11,2)=0.62903237_DP
    gbasis(53)%coco(11,3)=0.66628327_DP

    gbasis(53)%angmom(12)=2
    gbasis(53)%expo(12,1)=0.30900000000_DP
    gbasis(53)%coco(12,1)=1.0000000000_DP
    ! ====end I ==========================================


    ! cks: all data in this subroutine from here onwards was typed by Mark Robinson

    ! ====begin Xe ==========================================
    gbasis(54)%angmom(1) =0
    gbasis(54)%angmom(2) =0
    gbasis(54)%angmom(3) =1
    gbasis(54)%angmom(4) =0
    gbasis(54)%angmom(5) =1

    ! 3d:
    gbasis(54)%angmom(6)=2
    gbasis(54)%expo(6,1)=135.60300038_DP
    gbasis(54)%coco(6,1)=0.81873543290E-03_DP
    gbasis(54)%expo(6,2)=38.727062692_DP
    gbasis(54)%coco(6,2)=0.60897654151E-02_DP
    gbasis(54)%expo(6,3)=15.377328089_DP
    gbasis(54)%coco(6,3)=-0.92782985596E-02_DP

    ! 4s:
    gbasis(54)%angmom(7)=0
    gbasis(54)%expo(7,1)=7.295991196_DP
    gbasis(54)%expo(7,2)=2.841021154_DP
    gbasis(54)%expo(7,3)=1.250624506_DP
    gbasis(54)%coco(7,1)=-0.33061006_DP
    gbasis(54)%coco(7,2)=0.05761095_DP
    gbasis(54)%coco(7,3)=1.11557874_DP

    ! 4p:
    gbasis(54)%angmom(8)=1
    gbasis(54)%expo(8,1)=193.81432545_DP
    gbasis(54)%coco(8,1)=0.95394802497E-03_DP
    gbasis(54)%expo(8,2)=21.725228086_DP
    gbasis(54)%coco(8,2)=0.57393353332E-01_DP
    gbasis(54)%expo(8,3)=9.8891605641_DP
    gbasis(54)%coco(8,3)=-0.27974266640_DP

    ! 5s:
    gbasis(54)%angmom(9)=0
    gbasis(54)%expo(9,1)=6420.2481656_DP
    gbasis(54)%coco(9,1)=0.25092173886E-03_DP
    gbasis(54)%expo(9,2)=983.54530664_DP
    gbasis(54)%coco(9,2)=0.16251948178E-02_DP
    gbasis(54)%expo(9,3)=219.43881364_DP
    gbasis(54)%coco(9,3)=0.46037106451E-02_DP

    ! 4d:
    gbasis(54)%angmom(10)=2
    gbasis(54)%expo(10,1)=0.58050830139_DP
    gbasis(54)%coco(10,1)=1.0000000000_DP

    ! 5p:
    gbasis(54)%angmom(11)=1
    gbasis(54)%expo(11,1)=13.960683826_DP
    gbasis(54)%coco(11,1)=-0.50950157807E-01_DP
    gbasis(54)%expo(11,2)=4.0928947097_DP
    gbasis(54)%coco(11,2)=0.36669211800_DP
    gbasis(54)%expo(11,3)=2.2546815460_DP
    gbasis(54)%coco(11,3)=0.72619456861_DP

    ! 6s:
    gbasis(54)%angmom(12)=0
    gbasis(54)%expo(12,1)=2.6257211320_DP
    gbasis(54)%coco(12,1)=1.0000000000_DP

    ! 6p:
    gbasis(54)%angmom(13)=1
    gbasis(54)%expo(13,1)=0.52321923128_DP
    gbasis(54)%coco(13,1)=1.0000000000_DP


    ! ====begin Cs ==========================================
    gbasis(55)%angmom(1) =0
    gbasis(55)%angmom(2) =0
    gbasis(55)%angmom(3) =1
    gbasis(55)%angmom(4) =0
    gbasis(55)%angmom(5) =1
    gbasis(55)%angmom(6) =2
    gbasis(55)%angmom(7) =0
    gbasis(55)%angmom(8) =1
    gbasis(55)%angmom(9) =2

    ! 5s:
    gbasis(55)%angmom(10)=0
    gbasis(55)%expo(10,1)=5.8778113443_DP
    gbasis(55)%coco(10,1)=0.12859994983_DP
    gbasis(55)%expo(10,2)=4.3631538286_DP
    gbasis(55)%coco(10,2)=-0.34632569725_DP
    gbasis(55)%expo(10,3)=1.8048475155_DP
    gbasis(55)%coco(10,3)=0.69930637051_DP

    ! 5p:
    gbasis(55)%angmom(11)=1
    gbasis(55)%expo(11,1)=4.2751856154_DP
    gbasis(55)%coco(11,1)=0.45723074174E-01_DP
    gbasis(55)%expo(11,2)=1.9656663360_DP
    gbasis(55)%coco(11,2)=-0.25019961976_DP
    gbasis(55)%expo(11,3)=0.47689195212_DP
    gbasis(55)%coco(11,3)=0.55660850066_DP

    ! 6s:
    gbasis(55)%angmom(12)=0
    gbasis(55)%expo(12,1)=0.16384858778_DP
    gbasis(55)%coco(12,1)=1.0000000000_DP

    ! 6p:
    gbasis(55)%angmom(13)=1
    gbasis(55)%expo(13,1)=0.91450850296E-01_DP
    gbasis(55)%coco(13,1)=1.0000000000_DP

    ! 6d:
    gbasis(55)%angmom(14)=2
    gbasis(55)%expo(14,1)=0.13300000000_DP
    gbasis(55)%coco(14,1)=1.0000000000_DP


    ! ====begin Ba ==========================================
    gbasis(56)%angmom(1) =0
    gbasis(56)%angmom(2) =0
    gbasis(56)%angmom(3) =1
    gbasis(56)%angmom(4) =0
    gbasis(56)%angmom(5) =1
    gbasis(56)%angmom(6) =2
    gbasis(56)%angmom(7) =0
    gbasis(56)%angmom(8) =1
    gbasis(56)%angmom(9) =2

    ! 5s:
    gbasis(56)%angmom(10)=0
    gbasis(56)%expo(10,1)=2.3961900000_DP
    gbasis(56)%coco(10,1)=-5.9554766710_DP
    gbasis(56)%expo(10,2)=2.2433050000_DP
    gbasis(56)%coco(10,2)=6.6805648241_DP
    gbasis(56)%expo(10,3)=0.71740200000_DP
    gbasis(56)%coco(10,3)=-0.57055347037_DP

    ! 5p:
    gbasis(56)%angmom(11)=1
    gbasis(56)%expo(11,1)=2.9267420000_DP
    gbasis(56)%coco(11,1)=0.79815276764_DP
    gbasis(56)%expo(11,2)=2.5207180000_DP
    gbasis(56)%coco(11,2)=-1.0707548832_DP
    gbasis(56)%expo(11,3)=0.67578180854_DP
    gbasis(56)%coco(11,3)=0.34940781448_DP

    ! 6s:
    gbasis(56)%angmom(12)=0
    gbasis(56)%expo(12,1)=0.40933264188E-01_DP
    gbasis(56)%coco(12,1)=1.0000000000_DP

    ! 6p:
    gbasis(56)%angmom(13)=1
    gbasis(56)%expo(13,1)=0.33039565349E-01_DP
    gbasis(56)%coco(13,1)=1.0000000000_DP

    ! 6d:
    gbasis(56)%angmom(14)=2
    gbasis(56)%expo(14,1)=0.96631500000_DP
    gbasis(56)%coco(14,1)=-0.90893800000_DP
    gbasis(56)%expo(14,2)=0.89382800000_DP
    gbasis(56)%coco(14,2)=0.94724000000_DP
    gbasis(56)%expo(14,3)=0.27319500000_DP
    gbasis(56)%coco(14,3)=0.32205700000_DP


    ! ====begin La ==========================================
    gbasis(57)%angmom(1) =0
    gbasis(57)%angmom(2) =0
    gbasis(57)%angmom(3) =1
    gbasis(57)%angmom(4) =0
    gbasis(57)%angmom(5) =1
    gbasis(57)%angmom(6) =2
    gbasis(57)%angmom(7) =0
    gbasis(57)%angmom(8) =1
    gbasis(57)%angmom(9) =2

    ! 5s:
    gbasis(57)%angmom(10)=0
    gbasis(57)%expo(10,1)=5.0873990000_DP
    gbasis(57)%coco(10,1)=-0.44173707857_DP
    gbasis(57)%expo(10,2)=4.2709780000_DP
    gbasis(57)%coco(10,2)=0.85812489720_DP

    ! 5p:
    gbasis(57)%angmom(11)=1
    gbasis(57)%expo(11,1)=3.0251610000_DP
    gbasis(57)%coco(11,1)=-0.28580266655_DP
    gbasis(57)%expo(11,2)=2.3820950000_DP
    gbasis(57)%coco(11,2)=0.51718364644_DP
    gbasis(57)%expo(11,3)=0.67246658081_DP
    gbasis(57)%coco(11,3)=-0.44823115995_DP

    ! 6s:
    gbasis(57)%angmom(12)=0
    gbasis(57)%expo(12,1)=0.48954926461_DP
    gbasis(57)%coco(12,1)=1.0000000000_DP

    ! 4f:
    gbasis(57)%angmom(13)=3
    gbasis(57)%expo(13,1)=0.25683000000_DP
    gbasis(57)%coco(13,1)=1.0000000000_DP

    ! 5d:
    gbasis(57)%angmom(14)=2
    gbasis(57)%expo(14,1)=1.2667182018_DP
    gbasis(57)%coco(14,1)=-0.17445579381_DP
    gbasis(57)%expo(14,2)=0.89049882672_DP
    gbasis(57)%coco(14,2)=0.25137236211_DP
    gbasis(57)%expo(14,3)=0.32947291040_DP
    gbasis(57)%coco(14,3)=0.44604068240_DP

    ! 6p:
    gbasis(57)%angmom(15)=1
    gbasis(57)%expo(15,1)=0.25100000000E-01_DP
    gbasis(57)%coco(15,1)=1.0000000000_DP


    ! ====begin Ce ==========================================
    gbasis(58)%angmom(1) =0
    gbasis(58)%angmom(2) =0
    gbasis(58)%angmom(3) =1
    gbasis(58)%angmom(4) =0
    gbasis(58)%angmom(5) =1
    gbasis(58)%angmom(6) =2
    gbasis(58)%angmom(7) =0
    gbasis(58)%angmom(8) =1
    gbasis(58)%angmom(9) =2

    ! 5s:
    gbasis(58)%angmom(10)=0
    gbasis(58)%expo(10,1)=5.0873990000_DP
    gbasis(58)%coco(10,1)=-0.44173707857_DP
    gbasis(58)%expo(10,2)=4.2709780000_DP
    gbasis(58)%coco(10,2)=0.85812489720_DP

    ! 5p:
    gbasis(58)%angmom(11)=1
    gbasis(58)%expo(11,1)=3.0251610000_DP
    gbasis(58)%coco(11,1)=-0.28580266655_DP
    gbasis(58)%expo(11,2)=2.3820950000_DP
    gbasis(58)%coco(11,2)=0.51718364644_DP
    gbasis(58)%expo(11,3)=0.67246658081_DP
    gbasis(58)%coco(11,3)=-0.44823115995_DP

    ! 6s:
    gbasis(58)%angmom(12)=0
    gbasis(58)%expo(12,1)=0.9018000_DP
    gbasis(58)%coco(12,1)=1.0000000_DP
    gbasis(58)%expo(12,2)=0.3558000_DP
    gbasis(58)%coco(12,2)=1.0000000_DP
    gbasis(58)%expo(12,3)=0.2580000_DP
    gbasis(58)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(58)%angmom(13)=3
    gbasis(58)%expo(13,1)=54.7400000_DP
    gbasis(58)%coco(13,1)=1.0000000_DP
    gbasis(58)%expo(13,2)=18.2400000_DP
    gbasis(58)%coco(13,2)=1.0000000_DP
    gbasis(58)%expo(13,3)=6.9560000_DP
    gbasis(58)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(58)%angmom(14)=2
    gbasis(58)%expo(14,1)=3.6790000_DP
    gbasis(58)%coco(14,1)=1.0000000_DP
    gbasis(58)%expo(14,2)=0.5340000_DP
    gbasis(58)%coco(14,2)=1.0000000_DP
    gbasis(58)%expo(14,3)=0.2463000_DP
    gbasis(58)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(58)%angmom(15)=1
    gbasis(58)%expo(15,1)=1.8540000_DP
    gbasis(58)%coco(15,1)=1.0000000_DP
    gbasis(58)%expo(15,2)=0.2073000_DP
    gbasis(58)%coco(15,2)=1.0000000_DP
    gbasis(58)%expo(15,3)=0.1229000_DP
    gbasis(58)%coco(15,3)=1.0000000_DP


    ! ====begin Pr ==========================================
    gbasis(59)%angmom(1) =0
    gbasis(59)%angmom(2) =0
    gbasis(59)%angmom(3) =1
    gbasis(59)%angmom(4) =0
    gbasis(59)%angmom(5) =1
    gbasis(59)%angmom(6) =2
    gbasis(59)%angmom(7) =0
    gbasis(59)%angmom(8) =1
    gbasis(59)%angmom(9) =2
    gbasis(59)%angmom(10)=0
    gbasis(59)%angmom(11)=1

    ! 6s:
    gbasis(59)%angmom(12)=0
    gbasis(59)%expo(12,1)=2.5680000_DP
    gbasis(59)%coco(12,1)=1.0000000_DP
    gbasis(59)%expo(12,2)=0.3517000_DP
    gbasis(59)%coco(12,2)=1.0000000_DP
    gbasis(59)%expo(12,3)=0.3036000_DP
    gbasis(59)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(59)%angmom(13)=3
    gbasis(59)%expo(13,1)=59.5500000_DP
    gbasis(59)%coco(13,1)=1.0000000_DP
    gbasis(59)%expo(13,2)=18.2700000_DP
    gbasis(59)%coco(13,2)=1.0000000_DP
    gbasis(59)%expo(13,3)=7.2230000_DP
    gbasis(59)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(59)%angmom(14)=2
    gbasis(59)%expo(14,1)=4.4590000_DP
    gbasis(59)%coco(14,1)=1.0000000_DP
    gbasis(59)%expo(14,2)=0.5570000_DP
    gbasis(59)%coco(14,2)=1.0000000_DP
    gbasis(59)%expo(14,3)=0.2596000_DP
    gbasis(59)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(59)%angmom(15)=1
    gbasis(59)%expo(15,1)=0.9555000_DP
    gbasis(59)%coco(15,1)=1.0000000_DP
    gbasis(59)%expo(15,2)=0.6262000_DP
    gbasis(59)%coco(15,2)=1.0000000_DP
    gbasis(59)%expo(15,3)=0.4320000_DP
    gbasis(59)%coco(15,3)=1.0000000_DP


    ! ====begin Nd ==========================================
    gbasis(60)%angmom(1) =0
    gbasis(60)%angmom(2) =0
    gbasis(60)%angmom(3) =1
    gbasis(60)%angmom(4) =0
    gbasis(60)%angmom(5) =1
    gbasis(60)%angmom(6) =2
    gbasis(60)%angmom(7) =0
    gbasis(60)%angmom(8) =1
    gbasis(60)%angmom(9) =2
    gbasis(60)%angmom(10)=0
    gbasis(60)%angmom(11)=1

    ! 6s:
    gbasis(60)%angmom(12)=0
    gbasis(60)%expo(12,1)=2.5290000_DP
    gbasis(60)%coco(12,1)=1.0000000_DP
    gbasis(60)%expo(12,2)=0.3766000_DP
    gbasis(60)%coco(12,2)=1.0000000_DP
    gbasis(60)%expo(12,3)=0.3073000_DP
    gbasis(60)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(60)%angmom(13)=3
    gbasis(60)%expo(13,1)=62.6400000_DP
    gbasis(60)%coco(13,1)=1.0000000_DP
    gbasis(60)%expo(13,2)=19.7100000_DP
    gbasis(60)%coco(13,2)=1.0000000_DP
    gbasis(60)%expo(13,3)=7.8250000_DP
    gbasis(60)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(60)%angmom(14)=2
    gbasis(60)%expo(14,1)=4.3030000_DP
    gbasis(60)%coco(14,1)=1.0000000_DP
    gbasis(60)%expo(14,2)=0.5902000_DP
    gbasis(60)%coco(14,2)=1.0000000_DP
    gbasis(60)%expo(14,3)=0.2741000_DP
    gbasis(60)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(60)%angmom(15)=1
    gbasis(60)%expo(15,1)=2.1700000_DP
    gbasis(60)%coco(15,1)=1.0000000_DP
    gbasis(60)%expo(15,2)=0.3192000_DP
    gbasis(60)%coco(15,2)=1.0000000_DP
    gbasis(60)%expo(15,3)=0.1009000_DP
    gbasis(60)%coco(15,3)=1.0000000_DP


  end subroutine ngwf_data_cocos_and_expos41_60

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwf_data_cocos_and_expos61_80(gbasis)

    !===============================================================!
    ! Initialisation of Gaussian basis set exponents and contraction!
    ! coefficients. All the data of the STO-3G basis set are        !
    ! included for each element of the periodic table as well as    !
    ! the one "polarisation" function of the 6-31G* set.            !
    !---------------------------------------------------------------!
    ! This subroutine was created by Chris-Kriton Skylaris on       !
    ! 31/03/2006 but the actual STO-3G and 6-31G* parameters        !
    ! which appear below were typed by Mark Robinson, in April 2007 !
    !===============================================================!

    implicit none

    type(GTO_SET), intent(inout) :: gbasis(1: pub_num_gtatoms)


    ! cks: all data in this subroutine was typed by Mark Robinson


    ! ====begin Pm ==========================================
    gbasis(61)%angmom(1) =0
    gbasis(61)%angmom(2) =0
    gbasis(61)%angmom(3) =1
    gbasis(61)%angmom(4) =0
    gbasis(61)%angmom(5) =1
    gbasis(61)%angmom(6) =2
    gbasis(61)%angmom(7) =0
    gbasis(61)%angmom(8) =1
    gbasis(61)%angmom(9) =2
    gbasis(61)%angmom(10)=0
    gbasis(61)%angmom(11)=1

    ! 6s:
    gbasis(61)%angmom(12)=0
    gbasis(61)%expo(12,1)=2.3480000_DP
    gbasis(61)%coco(12,1)=1.0000000_DP
    gbasis(61)%expo(12,2)=0.3778000_DP
    gbasis(61)%coco(12,2)=1.0000000_DP
    gbasis(61)%expo(12,3)=0.3058000_DP
    gbasis(61)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(61)%angmom(13)=3
    gbasis(61)%expo(13,1)=70.8100000_DP
    gbasis(61)%coco(13,1)=1.0000000_DP
    gbasis(61)%expo(13,2)=23.1300000_DP
    gbasis(61)%coco(13,2)=1.0000000_DP
    gbasis(61)%expo(13,3)=9.1960000_DP
    gbasis(61)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(61)%angmom(14)=2
    gbasis(61)%expo(14,1)=4.2670000_DP
    gbasis(61)%coco(14,1)=1.0000000_DP
    gbasis(61)%expo(14,2)=0.6316000_DP
    gbasis(61)%coco(14,2)=1.0000000_DP
    gbasis(61)%expo(14,3)=0.2996000_DP
    gbasis(61)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(61)%angmom(15)=1
    gbasis(61)%expo(15,1)=2.4190000_DP
    gbasis(61)%coco(15,1)=1.0000000_DP
    gbasis(61)%expo(15,2)=0.3014000_DP
    gbasis(61)%coco(15,2)=1.0000000_DP
    gbasis(61)%expo(15,3)=0.1093000_DP
    gbasis(61)%coco(15,3)=1.0000000_DP


    ! ====begin Sm ==========================================
    gbasis(62)%angmom(1) =0
    gbasis(62)%angmom(2) =0
    gbasis(62)%angmom(3) =1
    gbasis(62)%angmom(4) =0
    gbasis(62)%angmom(5) =1
    gbasis(62)%angmom(6) =2
    gbasis(62)%angmom(7) =0
    gbasis(62)%angmom(8) =1
    gbasis(62)%angmom(9) =2
    gbasis(62)%angmom(10)=0
    gbasis(62)%angmom(11)=1

    ! 6s:
    gbasis(62)%angmom(12)=0
    gbasis(62)%expo(12,1)=21.5900000_DP
    gbasis(62)%coco(12,1)=1.0000000_DP
    gbasis(62)%expo(12,2)=0.3692000_DP
    gbasis(62)%coco(12,2)=1.0000000_DP
    gbasis(62)%expo(12,3)=0.2862000_DP
    gbasis(62)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(62)%angmom(13)=3
    gbasis(62)%expo(13,1)=76.0000000_DP
    gbasis(62)%coco(13,1)=1.0000000_DP
    gbasis(62)%expo(13,2)=25.7200000_DP
    gbasis(62)%coco(13,2)=1.0000000_DP
    gbasis(62)%expo(13,3)=10.2600000_DP
    gbasis(62)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(62)%angmom(14)=2
    gbasis(62)%expo(14,1)=4.2730000_DP
    gbasis(62)%coco(14,1)=1.0000000_DP
    gbasis(62)%expo(14,2)=0.6648000_DP
    gbasis(62)%coco(14,2)=1.0000000_DP
    gbasis(62)%expo(14,3)=0.3131000_DP
    gbasis(62)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(62)%angmom(15)=1
    gbasis(62)%expo(15,1)=2.7330000_DP
    gbasis(62)%coco(15,1)=1.0000000_DP
    gbasis(62)%expo(15,2)=0.2053000_DP
    gbasis(62)%coco(15,2)=1.0000000_DP
    gbasis(62)%expo(15,3)=0.1191000_DP
    gbasis(62)%coco(15,3)=1.0000000_DP


    ! ====begin Eu ==========================================
    gbasis(63)%angmom(1) =0
    gbasis(63)%angmom(2) =0
    gbasis(63)%angmom(3) =1
    gbasis(63)%angmom(4) =0
    gbasis(63)%angmom(5) =1
    gbasis(63)%angmom(6) =2
    gbasis(63)%angmom(7) =0
    gbasis(63)%angmom(8) =1
    gbasis(63)%angmom(9) =2
    gbasis(63)%angmom(10)=0
    gbasis(63)%angmom(11)=1

    ! 6s:
    gbasis(63)%angmom(12)=0
    gbasis(63)%expo(12,1)=16.4800000_DP
    gbasis(63)%coco(12,1)=1.0000000_DP
    gbasis(63)%expo(12,2)=0.4074000_DP
    gbasis(63)%coco(12,2)=1.0000000_DP
    gbasis(63)%expo(12,3)=0.2855000_DP
    gbasis(63)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(63)%angmom(13)=3
    gbasis(63)%expo(13,1)=74.3400000_DP
    gbasis(63)%coco(13,1)=1.0000000_DP
    gbasis(63)%expo(13,2)=24.3900000_DP
    gbasis(63)%coco(13,2)=1.0000000_DP
    gbasis(63)%expo(13,3)=9.6550000_DP
    gbasis(63)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(63)%angmom(14)=2
    gbasis(63)%expo(14,1)=4.8290000_DP
    gbasis(63)%coco(14,1)=1.0000000_DP
    gbasis(63)%expo(14,2)=0.6861000_DP
    gbasis(63)%coco(14,2)=1.0000000_DP
    gbasis(63)%expo(14,3)=0.3226000_DP
    gbasis(63)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(63)%angmom(15)=1
    gbasis(63)%expo(15,1)=2.8170000_DP
    gbasis(63)%coco(15,1)=1.0000000_DP
    gbasis(63)%expo(15,2)=0.2188000_DP
    gbasis(63)%coco(15,2)=1.0000000_DP
    gbasis(63)%expo(15,3)=0.1140000_DP
    gbasis(63)%coco(15,3)=1.0000000_DP


    ! ====begin Gd ==========================================
    gbasis(64)%angmom(1) =0
    gbasis(64)%angmom(2) =0
    gbasis(64)%angmom(3) =1
    gbasis(64)%angmom(4) =0
    gbasis(64)%angmom(5) =1
    gbasis(64)%angmom(6) =2
    gbasis(64)%angmom(7) =0
    gbasis(64)%angmom(8) =1
    gbasis(64)%angmom(9) =2
    gbasis(64)%angmom(10)=0
    gbasis(64)%angmom(11)=1

    ! 6s:
    gbasis(64)%angmom(12)=0
    gbasis(64)%expo(12,1)=17.4000000_DP
    gbasis(64)%coco(12,1)=1.0000000_DP
    gbasis(64)%expo(12,2)=0.4215000_DP
    gbasis(64)%coco(12,2)=1.0000000_DP
    gbasis(64)%expo(12,3)=0.2946000_DP
    gbasis(64)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(64)%angmom(13)=3
    gbasis(64)%expo(13,1)=79.0300000_DP
    gbasis(64)%coco(13,1)=1.0000000_DP
    gbasis(64)%expo(13,2)=25.7300000_DP
    gbasis(64)%coco(13,2)=1.0000000_DP
    gbasis(64)%expo(13,3)=10.1900000_DP
    gbasis(64)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(64)%angmom(14)=2
    gbasis(64)%expo(14,1)=4.8620000_DP
    gbasis(64)%coco(14,1)=1.0000000_DP
    gbasis(64)%expo(14,2)=0.7236000_DP
    gbasis(64)%coco(14,2)=1.0000000_DP
    gbasis(64)%expo(14,3)=0.3426000_DP
    gbasis(64)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(64)%angmom(15)=1
    gbasis(64)%expo(15,1)=2.8130000_DP
    gbasis(64)%coco(15,1)=1.0000000_DP
    gbasis(64)%expo(15,2)=0.1974000_DP
    gbasis(64)%coco(15,2)=1.0000000_DP
    gbasis(64)%expo(15,3)=0.1504000_DP
    gbasis(64)%coco(15,3)=1.0000000_DP


    ! ====begin Tb ==========================================
    gbasis(65)%angmom(1) =0
    gbasis(65)%angmom(2) =0
    gbasis(65)%angmom(3) =1
    gbasis(65)%angmom(4) =0
    gbasis(65)%angmom(5) =1
    gbasis(65)%angmom(6) =2
    gbasis(65)%angmom(7) =0
    gbasis(65)%angmom(8) =1
    gbasis(65)%angmom(9) =2
    gbasis(65)%angmom(10)=0
    gbasis(65)%angmom(11)=1

    ! 6s:
    gbasis(65)%angmom(12)=0
    gbasis(65)%expo(12,1)=17.5500000_DP
    gbasis(65)%coco(12,1)=1.0000000_DP
    gbasis(65)%expo(12,2)=0.6583000_DP
    gbasis(65)%coco(12,2)=1.0000000_DP
    gbasis(65)%expo(12,3)=0.2507000_DP
    gbasis(65)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(65)%angmom(13)=3
    gbasis(65)%expo(13,1)=85.7900000_DP
    gbasis(65)%coco(13,1)=1.0000000_DP
    gbasis(65)%expo(13,2)=27.7600000_DP
    gbasis(65)%coco(13,2)=1.0000000_DP
    gbasis(65)%expo(13,3)=10.9800000_DP
    gbasis(65)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(65)%angmom(14)=2
    gbasis(65)%expo(14,1)=4.9970000_DP
    gbasis(65)%coco(14,1)=1.0000000_DP
    gbasis(65)%expo(14,2)=0.7526000_DP
    gbasis(65)%coco(14,2)=1.0000000_DP
    gbasis(65)%expo(14,3)=0.3470000_DP
    gbasis(65)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(65)%angmom(15)=1
    gbasis(65)%expo(15,1)=2.7960000_DP
    gbasis(65)%coco(15,1)=1.0000000_DP
    gbasis(65)%expo(15,2)=0.2043000_DP
    gbasis(65)%coco(15,2)=1.0000000_DP
    gbasis(65)%expo(15,3)=0.1503000_DP
    gbasis(65)%coco(15,3)=1.0000000_DP


    ! ====begin Dy ==========================================
    gbasis(66)%angmom(1) =0
    gbasis(66)%angmom(2) =0
    gbasis(66)%angmom(3) =1
    gbasis(66)%angmom(4) =0
    gbasis(66)%angmom(5) =1
    gbasis(66)%angmom(6) =2
    gbasis(66)%angmom(7) =0
    gbasis(66)%angmom(8) =1
    gbasis(66)%angmom(9) =2
    gbasis(66)%angmom(10)=0
    gbasis(66)%angmom(11)=1

    ! 6s:
    gbasis(66)%angmom(12)=0
    gbasis(66)%expo(12,1)=14.8900000_DP
    gbasis(66)%coco(12,1)=1.0000000_DP
    gbasis(66)%expo(12,2)=0.8004000_DP
    gbasis(66)%coco(12,2)=1.0000000_DP
    gbasis(66)%expo(12,3)=0.2464000_DP
    gbasis(66)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(66)%angmom(13)=3
    gbasis(66)%expo(13,1)=91.0600000_DP
    gbasis(66)%coco(13,1)=1.0000000_DP
    gbasis(66)%expo(13,2)=30.8800000_DP
    gbasis(66)%coco(13,2)=1.0000000_DP
    gbasis(66)%expo(13,3)=11.7900000_DP
    gbasis(66)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(66)%angmom(14)=2
    gbasis(66)%expo(14,1)=5.9300000_DP
    gbasis(66)%coco(14,1)=1.0000000_DP
    gbasis(66)%expo(14,2)=0.8184000_DP
    gbasis(66)%coco(14,2)=1.0000000_DP
    gbasis(66)%expo(14,3)=0.3938000_DP
    gbasis(66)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(66)%angmom(15)=1
    gbasis(66)%expo(15,1)=4.8570000_DP
    gbasis(66)%coco(15,1)=1.0000000_DP
    gbasis(66)%expo(15,2)=0.2386000_DP
    gbasis(66)%coco(15,2)=1.0000000_DP
    gbasis(66)%expo(15,3)=0.1759000_DP
    gbasis(66)%coco(15,3)=1.0000000_DP


    ! ====begin Ho ==========================================
    gbasis(67)%angmom(1) =0
    gbasis(67)%angmom(2) =0
    gbasis(67)%angmom(3) =1
    gbasis(67)%angmom(4) =0
    gbasis(67)%angmom(5) =1
    gbasis(67)%angmom(6) =2
    gbasis(67)%angmom(7) =0
    gbasis(67)%angmom(8) =1
    gbasis(67)%angmom(9) =2
    gbasis(67)%angmom(10)=0
    gbasis(67)%angmom(11)=1

    ! 6s:
    gbasis(67)%angmom(12)=0
    gbasis(67)%expo(12,1)=11.0700000_DP
    gbasis(67)%coco(12,1)=1.0000000_DP
    gbasis(67)%expo(12,2)=0.4329000_DP
    gbasis(67)%coco(12,2)=1.0000000_DP
    gbasis(67)%expo(12,3)=0.3406000_DP
    gbasis(67)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(67)%angmom(13)=3
    gbasis(67)%expo(13,1)=101.6000000_DP
    gbasis(67)%coco(13,1)=1.0000000_DP
    gbasis(67)%expo(13,2)=34.3800000_DP
    gbasis(67)%coco(13,2)=1.0000000_DP
    gbasis(67)%expo(13,3)=12.9400000_DP
    gbasis(67)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(67)%angmom(14)=2
    gbasis(67)%expo(14,1)=5.3680000_DP
    gbasis(67)%coco(14,1)=1.0000000_DP
    gbasis(67)%expo(14,2)=0.8922000_DP
    gbasis(67)%coco(14,2)=1.0000000_DP
    gbasis(67)%expo(14,3)=0.3919000_DP
    gbasis(67)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(67)%angmom(15)=1
    gbasis(67)%expo(15,1)=2.5840000_DP
    gbasis(67)%coco(15,1)=1.0000000_DP
    gbasis(67)%expo(15,2)=0.5200000_DP
    gbasis(67)%coco(15,2)=1.0000000_DP
    gbasis(67)%expo(15,3)=0.2019000_DP
    gbasis(67)%coco(15,3)=1.0000000_DP


    ! ====begin Er ==========================================
    gbasis(68)%angmom(1) =0
    gbasis(68)%angmom(2) =0
    gbasis(68)%angmom(3) =1
    gbasis(68)%angmom(4) =0
    gbasis(68)%angmom(5) =1
    gbasis(68)%angmom(6) =2
    gbasis(68)%angmom(7) =0
    gbasis(68)%angmom(8) =1
    gbasis(68)%angmom(9) =2
    gbasis(68)%angmom(10)=0
    gbasis(68)%angmom(11)=1

    ! 6s:
    gbasis(68)%angmom(12)=0
    gbasis(68)%expo(12,1)=18.0900000_DP
    gbasis(68)%coco(12,1)=1.0000000_DP
    gbasis(68)%expo(12,2)=0.5156000_DP
    gbasis(68)%coco(12,2)=1.0000000_DP
    gbasis(68)%expo(12,3)=0.3522000_DP
    gbasis(68)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(68)%angmom(13)=3
    gbasis(68)%expo(13,1)=104.5000000_DP
    gbasis(68)%coco(13,1)=1.0000000_DP
    gbasis(68)%expo(13,2)=35.8800000_DP
    gbasis(68)%coco(13,2)=1.0000000_DP
    gbasis(68)%expo(13,3)=13.5500000_DP
    gbasis(68)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(68)%angmom(14)=2
    gbasis(68)%expo(14,1)=6.3310000_DP
    gbasis(68)%coco(14,1)=1.0000000_DP
    gbasis(68)%expo(14,2)=0.8450000_DP
    gbasis(68)%coco(14,2)=1.0000000_DP
    gbasis(68)%expo(14,3)=0.3822000_DP
    gbasis(68)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(68)%angmom(15)=1
    gbasis(68)%expo(15,1)=2.4370000_DP
    gbasis(68)%coco(15,1)=1.0000000_DP
    gbasis(68)%expo(15,2)=0.2513000_DP
    gbasis(68)%coco(15,2)=1.0000000_DP
    gbasis(68)%expo(15,3)=0.1799000_DP
    gbasis(68)%coco(15,3)=1.0000000_DP


    ! ====begin Tm ==========================================
    gbasis(69)%angmom(1) =0
    gbasis(69)%angmom(2) =0
    gbasis(69)%angmom(3) =1
    gbasis(69)%angmom(4) =0
    gbasis(69)%angmom(5) =1
    gbasis(69)%angmom(6) =2
    gbasis(69)%angmom(7) =0
    gbasis(69)%angmom(8) =1
    gbasis(69)%angmom(9) =2
    gbasis(69)%angmom(10)=0
    gbasis(69)%angmom(11)=1

    ! 6s:
    gbasis(69)%angmom(12)=0
    gbasis(69)%expo(12,1)=13.2500000_DP
    gbasis(69)%coco(12,1)=1.0000000_DP
    gbasis(69)%expo(12,2)=0.5350000_DP
    gbasis(69)%coco(12,2)=1.0000000_DP
    gbasis(69)%expo(12,3)=0.3677000_DP
    gbasis(69)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(69)%angmom(13)=3
    gbasis(69)%expo(13,1)=106.5000000_DP
    gbasis(69)%coco(13,1)=1.0000000_DP
    gbasis(69)%expo(13,2)=35.6900000_DP
    gbasis(69)%coco(13,2)=1.0000000_DP
    gbasis(69)%expo(13,3)=13.5400000_DP
    gbasis(69)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(69)%angmom(14)=2
    gbasis(69)%expo(14,1)=6.4840000_DP
    gbasis(69)%coco(14,1)=1.0000000_DP
    gbasis(69)%expo(14,2)=0.8741000_DP
    gbasis(69)%coco(14,2)=1.0000000_DP
    gbasis(69)%expo(14,3)=0.3826000_DP
    gbasis(69)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(69)%angmom(15)=1
    gbasis(69)%expo(15,1)=2.4040000_DP
    gbasis(69)%coco(15,1)=1.0000000_DP
    gbasis(69)%expo(15,2)=0.2581000_DP
    gbasis(69)%coco(15,2)=1.0000000_DP
    gbasis(69)%expo(15,3)=0.1840000_DP
    gbasis(69)%coco(15,3)=1.0000000_DP


    ! ====begin Yb ==========================================
    gbasis(70)%angmom(1) =0
    gbasis(70)%angmom(2) =0
    gbasis(70)%angmom(3) =1
    gbasis(70)%angmom(4) =0
    gbasis(70)%angmom(5) =1
    gbasis(70)%angmom(6) =2
    gbasis(70)%angmom(7) =0
    gbasis(70)%angmom(8) =1
    gbasis(70)%angmom(9) =2
    gbasis(70)%angmom(10)=0
    gbasis(70)%angmom(11)=1

    ! 6s:
    gbasis(70)%angmom(12)=0
    gbasis(70)%expo(12,1)=15.0100000_DP
    gbasis(70)%coco(12,1)=1.0000000_DP
    gbasis(70)%expo(12,2)=0.4727000_DP
    gbasis(70)%coco(12,2)=1.0000000_DP
    gbasis(70)%expo(12,3)=0.3908000_DP
    gbasis(70)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(70)%angmom(13)=3
    gbasis(70)%expo(13,1)=125.5000000_DP
    gbasis(70)%coco(13,1)=1.0000000_DP
    gbasis(70)%expo(13,2)=40.5500000_DP
    gbasis(70)%coco(13,2)=1.0000000_DP
    gbasis(70)%expo(13,3)=14.9800000_DP
    gbasis(70)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(70)%angmom(14)=2
    gbasis(70)%expo(14,1)=7.0050000_DP
    gbasis(70)%coco(14,1)=1.0000000_DP
    gbasis(70)%expo(14,2)=0.9656000_DP
    gbasis(70)%coco(14,2)=1.0000000_DP
    gbasis(70)%expo(14,3)=0.4348000_DP
    gbasis(70)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(70)%angmom(15)=1
    gbasis(70)%expo(15,1)=1.4250000_DP
    gbasis(70)%coco(15,1)=1.0000000_DP
    gbasis(70)%expo(15,2)=0.2992000_DP
    gbasis(70)%coco(15,2)=1.0000000_DP
    gbasis(70)%expo(15,3)=0.2278000_DP
    gbasis(70)%coco(15,3)=1.0000000_DP


    ! ====begin Lu ==========================================
    gbasis(71)%angmom(1) =0
    gbasis(71)%angmom(2) =0
    gbasis(71)%angmom(3) =1
    gbasis(71)%angmom(4) =0
    gbasis(71)%angmom(5) =1
    gbasis(71)%angmom(6) =2
    gbasis(71)%angmom(7) =0
    gbasis(71)%angmom(8) =1
    gbasis(71)%angmom(9) =2
    gbasis(71)%angmom(10)=0
    gbasis(71)%angmom(11)=1

    ! 6s:
    gbasis(71)%angmom(12)=0
    gbasis(71)%expo(12,1)=5.4060000_DP
    gbasis(71)%coco(12,1)=1.0000000_DP
    gbasis(71)%expo(12,2)=0.5869000_DP
    gbasis(71)%coco(12,2)=1.0000000_DP
    gbasis(71)%expo(12,3)=0.4803000_DP
    gbasis(71)%coco(12,3)=1.0000000_DP

    ! 4f:
    gbasis(71)%angmom(13)=3
    gbasis(71)%expo(13,1)=120.7000000_DP
    gbasis(71)%coco(13,1)=1.0000000_DP
    gbasis(71)%expo(13,2)=40.2300000_DP
    gbasis(71)%coco(13,2)=1.0000000_DP
    gbasis(71)%expo(13,3)=15.3400000_DP
    gbasis(71)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(71)%angmom(14)=2
    gbasis(71)%expo(14,1)=9.1770000_DP
    gbasis(71)%coco(14,1)=1.0000000_DP
    gbasis(71)%expo(14,2)=0.8922000_DP
    gbasis(71)%coco(14,2)=1.0000000_DP
    gbasis(71)%expo(14,3)=0.3919000_DP
    gbasis(71)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(71)%angmom(15)=1
    gbasis(71)%expo(15,1)=1.1320000_DP
    gbasis(71)%coco(15,1)=1.0000000_DP
    gbasis(71)%expo(15,2)=0.6640000_DP
    gbasis(71)%coco(15,2)=1.0000000_DP
    gbasis(71)%expo(15,3)=0.3620000_DP
    gbasis(71)%coco(15,3)=1.0000000_DP


    ! ====begin Hf ==========================================
    gbasis(72)%angmom(1) =0
    gbasis(72)%angmom(2) =0
    gbasis(72)%angmom(3) =1
    gbasis(72)%angmom(4) =0
    gbasis(72)%angmom(5) =1
    gbasis(72)%angmom(6) =2
    gbasis(72)%angmom(7) =0
    gbasis(72)%angmom(8) =1
    gbasis(72)%angmom(9) =2
    gbasis(72)%angmom(10)=0
    gbasis(72)%angmom(11)=1

    ! 4f:
    gbasis(72)%angmom(12)=3
    gbasis(72)%expo(12,1)=0.31547000000_DP
    gbasis(72)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(72)%angmom(13)=2
    gbasis(72)%expo(13,1)=3.9820623229_DP
    gbasis(72)%coco(13,1)=-0.42446023944E-01_DP
    gbasis(72)%expo(13,2)=1.3077987772_DP
    gbasis(72)%coco(13,2)=0.18409232564_DP
    gbasis(72)%expo(13,3)=0.53272310298_DP
    gbasis(72)%coco(13,3)=0.41484423516_DP

    ! 6s:
    gbasis(72)%angmom(14)=0
    gbasis(72)%expo(14,1)=2.55500_DP
    gbasis(72)%coco(14,1)=1.0_DP

    ! 6p:
    gbasis(72)%angmom(15)=1
    gbasis(72)%expo(15,1)=6.7265310000_DP
    gbasis(72)%coco(15,1)=0.65076205466_DP
    gbasis(72)%expo(15,2)=5.9599790000_DP
    gbasis(72)%coco(15,2)=-0.83713508620_DP
    gbasis(72)%expo(15,3)=1.0921117054_DP
    gbasis(72)%coco(15,3)=0.50570506286_DP


    ! ====begin Ta ==========================================
    gbasis(73)%angmom(1) =0
    gbasis(73)%angmom(2) =0
    gbasis(73)%angmom(3) =1
    gbasis(73)%angmom(4) =0
    gbasis(73)%angmom(5) =1
    gbasis(73)%angmom(6) =2
    gbasis(73)%angmom(7) =0
    gbasis(73)%angmom(8) =1
    gbasis(73)%angmom(9) =2
    gbasis(73)%angmom(10)=0
    gbasis(73)%angmom(11)=1

    ! 4f:
    gbasis(73)%angmom(12)=3
    gbasis(73)%expo(12,1)=0.37386000000_DP
    gbasis(73)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(73)%angmom(13)=2
    gbasis(73)%expo(13,1)=3.9738796278_DP
    gbasis(73)%coco(13,1)=-0.52799310714E-01_DP
    gbasis(73)%expo(13,2)=1.4528884813_DP
    gbasis(73)%coco(13,2)=0.18558319471_DP
    gbasis(73)%expo(13,3)=0.61042908544_DP
    gbasis(73)%coco(13,3)=0.42959071631_DP

    ! 6s:
    gbasis(73)%angmom(14)=0
    gbasis(73)%expo(14,1)=5.2740000_DP
    gbasis(73)%coco(14,1)=1.0_DP

    ! 6p:
    gbasis(73)%angmom(15)=1
    gbasis(73)%expo(15,1)=7.4188720000_DP
    gbasis(73)%coco(15,1)=0.26979695152_DP
    gbasis(73)%expo(15,2)=5.6984100000_DP
    gbasis(73)%coco(15,2)=-0.46968874449_DP
    gbasis(73)%expo(15,3)=1.1777211960_DP
    gbasis(73)%coco(15,3)=0.50905100155_DP


    ! ====begin W ==========================================
    gbasis(74)%angmom(1) =0
    gbasis(74)%angmom(2) =0
    gbasis(74)%angmom(3) =1
    gbasis(74)%angmom(4) =0
    gbasis(74)%angmom(5) =1
    gbasis(74)%angmom(6) =2
    gbasis(74)%angmom(7) =0
    gbasis(74)%angmom(8) =1
    gbasis(74)%angmom(9) =2
    gbasis(74)%angmom(10)=0
    gbasis(74)%angmom(11)=1

    ! 4f:
    gbasis(74)%angmom(12)=3
    gbasis(74)%expo(12,1)=0.43199000000_DP
    gbasis(74)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(74)%angmom(13)=2
    gbasis(74)%expo(13,1)=4.0131231332_DP
    gbasis(74)%coco(13,1)=-0.63416491808E-01_DP
    gbasis(74)%expo(13,2)=1.6237452450_DP
    gbasis(74)%coco(13,2)=0.18452733732_DP
    gbasis(74)%expo(13,3)=0.69187452392_DP
    gbasis(74)%coco(13,3)=0.44121676843_DP

    ! 6s:
    gbasis(74)%angmom(14)=0
    gbasis(74)%expo(14,1)=2.8533000_DP
    gbasis(74)%coco(14,1)=1.0_DP

    ! 6p:
    gbasis(74)%angmom(15)=1
    gbasis(74)%expo(15,1)=7.2496570000_DP
    gbasis(74)%coco(15,1)=0.45857265530_DP
    gbasis(74)%expo(15,2)=6.0848760000_DP
    gbasis(74)%coco(15,2)=-0.66720686321_DP
    gbasis(74)%expo(15,3)=1.2523777812_DP
    gbasis(74)%coco(15,3)=0.52090351989_DP


    ! ====begin Re ==========================================
    gbasis(75)%angmom(1) =0
    gbasis(75)%angmom(2) =0
    gbasis(75)%angmom(3) =1
    gbasis(75)%angmom(4) =0
    gbasis(75)%angmom(5) =1
    gbasis(75)%angmom(6) =2
    gbasis(75)%angmom(7) =0
    gbasis(75)%angmom(8) =1
    gbasis(75)%angmom(9) =2
    gbasis(75)%angmom(10)=0
    gbasis(75)%angmom(11)=1

    ! 4f:
    gbasis(75)%angmom(12)=3
    gbasis(75)%expo(12,1)=0.48989000000_DP
    gbasis(75)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(75)%angmom(13)=2
    gbasis(75)%expo(13,1)=4.0616274916_DP
    gbasis(75)%coco(13,1)=-0.76667295205E-01_DP
    gbasis(75)%expo(13,2)=1.8487943286_DP
    gbasis(75)%coco(13,2)=0.18067495215_DP
    gbasis(75)%expo(13,3)=0.78702048429_DP
    gbasis(75)%coco(13,3)=0.45160766937_DP

    ! 6s:
    gbasis(75)%angmom(14)=0
    gbasis(75)%expo(14,1)=2.676_DP
    gbasis(75)%coco(14,1)=1.0_DP

    ! 6p:
    gbasis(75)%angmom(15)=1
    gbasis(75)%expo(15,1)=7.4040050000_DP
    gbasis(75)%coco(15,1)=0.54510937865_DP
    gbasis(75)%expo(15,2)=6.3502060000_DP
    gbasis(75)%coco(15,2)=-0.76096056536_DP
    gbasis(75)%expo(15,3)=1.3311713455_DP
    gbasis(75)%coco(15,3)=0.53252276765_DP


    ! ====begin Os ==========================================
    gbasis(76)%angmom(1) =0
    gbasis(76)%angmom(2) =0
    gbasis(76)%angmom(3) =1
    gbasis(76)%angmom(4) =0
    gbasis(76)%angmom(5) =1
    gbasis(76)%angmom(6) =2
    gbasis(76)%angmom(7) =0
    gbasis(76)%angmom(8) =1
    gbasis(76)%angmom(9) =2
    gbasis(76)%angmom(10)=0
    gbasis(76)%angmom(11)=1

    ! 4f:
    gbasis(76)%angmom(12)=3
    gbasis(76)%expo(12,1)=0.55065000000_DP
    gbasis(76)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(76)%angmom(13)=2
    gbasis(76)%expo(13,1)=3.5808150000_DP
    gbasis(76)%coco(13,1)=-0.42773932262_DP
    gbasis(76)%expo(13,2)=3.1961650000_DP
    gbasis(76)%coco(13,2)=0.44501403992_DP
    gbasis(76)%expo(13,3)=1.0550226354_DP
    gbasis(76)%coco(13,3)=0.43689459722_DP

    ! 6s:
    gbasis(76)%angmom(14)=0
    gbasis(76)%expo(14,1)=2.2173_DP
    gbasis(76)%coco(14,1)=1.0_DP

    ! 6p:
    gbasis(76)%angmom(15)=1
    gbasis(76)%expo(15,1)=7.9362790000_DP
    gbasis(76)%coco(15,1)=0.33899369659_DP
    gbasis(76)%expo(15,2)=6.3036410000_DP
    gbasis(76)%coco(15,2)=-0.56379438345_DP
    gbasis(76)%expo(15,3)=1.4351579676_DP
    gbasis(76)%coco(15,3)=0.52360040577_DP


    ! ====begin Ir ==========================================
    gbasis(77)%angmom(1) =0
    gbasis(77)%angmom(2) =0
    gbasis(77)%angmom(3) =1
    gbasis(77)%angmom(4) =0
    gbasis(77)%angmom(5) =1
    gbasis(77)%angmom(6) =2
    gbasis(77)%angmom(7) =0
    gbasis(77)%angmom(8) =1
    gbasis(77)%angmom(9) =2
    gbasis(77)%angmom(10)=0
    gbasis(77)%angmom(11)=1

    ! 4f:
    gbasis(77)%angmom(12)=3
    gbasis(77)%expo(12,1)=0.61008000000_DP
    gbasis(77)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(77)%angmom(13)=2
    gbasis(77)%expo(13,1)=3.6994850000_DP
    gbasis(77)%coco(13,1)=-0.54840155965_DP
    gbasis(77)%expo(13,2)=3.3617430000_DP
    gbasis(77)%coco(13,2)=0.56835425804_DP
    gbasis(77)%expo(13,3)=1.1177053710_DP
    gbasis(77)%coco(13,3)=0.44304381513_DP

    ! 6s:
    gbasis(77)%angmom(14)=0
    gbasis(77)%expo(14,1)=2.1794_DP
    gbasis(77)%coco(14,1)=1.0_DP

    ! 6p:
    gbasis(77)%angmom(15)=1
    gbasis(77)%expo(15,1)=8.6697960000_DP
    gbasis(77)%coco(15,1)=0.21469187346_DP
    gbasis(77)%expo(15,2)=6.2456140000_DP
    gbasis(77)%coco(15,2)=-0.44863071815_DP
    gbasis(77)%expo(15,3)=1.5657036153_DP
    gbasis(77)%coco(15,3)=0.50653393257_DP


    ! ====begin Pt ==========================================
    gbasis(78)%angmom(1) =0
    gbasis(78)%angmom(2) =0
    gbasis(78)%angmom(3) =1
    gbasis(78)%angmom(4) =0
    gbasis(78)%angmom(5) =1
    gbasis(78)%angmom(6) =2
    gbasis(78)%angmom(7) =0
    gbasis(78)%angmom(8) =1
    gbasis(78)%angmom(9) =2
    gbasis(78)%angmom(10)=0
    gbasis(78)%angmom(11)=1

    ! 4f:
    gbasis(78)%angmom(12)=3
    gbasis(78)%expo(12,1)=0.66813000000_DP
    gbasis(78)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(78)%angmom(13)=2
    gbasis(78)%expo(13,1)=4.6299536825_DP
    gbasis(78)%coco(13,1)=-0.87774450596E-01_DP
    gbasis(78)%expo(13,2)=2.1980241252_DP
    gbasis(78)%coco(13,2)=0.21158360681_DP
    gbasis(78)%expo(13,3)=0.93629991261_DP
    gbasis(78)%coco(13,3)=0.46533857641_DP

    ! 6s:
    gbasis(78)%angmom(14)=0
    gbasis(78)%expo(14,1)=2.3654_DP
    gbasis(78)%coco(14,1)=1.0_DP

    ! 6p:
    gbasis(78)%angmom(15)=1
    gbasis(78)%expo(15,1)=8.1000000000_DP
    gbasis(78)%coco(15,1)=0.72955608128_DP
    gbasis(78)%expo(15,2)=7.2000000000_DP
    gbasis(78)%coco(15,2)=-0.95441807252_DP
    gbasis(78)%expo(15,3)=1.5588402917_DP
    gbasis(78)%coco(15,3)=0.57140490320_DP


    ! ====begin Au ==========================================
    gbasis(79)%angmom(1) =0
    gbasis(79)%angmom(2) =0
    gbasis(79)%angmom(3) =1
    gbasis(79)%angmom(4) =0
    gbasis(79)%angmom(5) =1
    gbasis(79)%angmom(6) =2
    gbasis(79)%angmom(7) =0
    gbasis(79)%angmom(8) =1
    gbasis(79)%angmom(9) =2
    gbasis(79)%angmom(10)=0
    gbasis(79)%angmom(11)=1

    ! 4f:
    gbasis(79)%angmom(12)=3
    gbasis(79)%expo(12,1)=0.72482000000_DP
    gbasis(79)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(79)%angmom(13)=2
    gbasis(79)%expo(13,1)=4.1439490000_DP
    gbasis(79)%coco(13,1)=-0.37099566643_DP
    gbasis(79)%expo(13,2)=3.5682570000_DP
    gbasis(79)%coco(13,2)=0.40197233762_DP
    gbasis(79)%expo(13,3)=1.2345757130_DP
    gbasis(79)%coco(13,3)=0.46001988624_DP

    ! 6s:
    gbasis(79)%angmom(14)=0
    gbasis(79)%expo(14,1)=2.4651_DP
    gbasis(79)%coco(14,1)=1.0_DP

    ! 6p:
    gbasis(79)%angmom(15)=1
    gbasis(79)%expo(15,1)=8.6096650000_DP
    gbasis(79)%coco(15,1)=0.50053018599_DP
    gbasis(79)%expo(15,2)=7.3353260000_DP
    gbasis(79)%coco(15,2)=-0.72681584494_DP
    gbasis(79)%expo(15,3)=1.6575296365_DP
    gbasis(79)%coco(15,3)=0.57315511417_DP


    ! ====begin Hg ==========================================
    gbasis(80)%angmom(1) =0
    gbasis(80)%angmom(2) =0
    gbasis(80)%angmom(3) =1
    gbasis(80)%angmom(4) =0
    gbasis(80)%angmom(5) =1
    gbasis(80)%angmom(6) =2
    gbasis(80)%angmom(7) =0
    gbasis(80)%angmom(8) =1
    gbasis(80)%angmom(9) =2
    gbasis(80)%angmom(10)=0
    gbasis(80)%angmom(11)=1

    ! 4f:
    gbasis(80)%angmom(12)=3
    gbasis(80)%expo(12,1)=0.79569000000_DP
    gbasis(80)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(80)%angmom(13)=2
    gbasis(80)%expo(13,1)=4.2800000000_DP
    gbasis(80)%coco(13,1)=-0.32033414305_DP
    gbasis(80)%expo(13,2)=3.4800000000_DP
    gbasis(80)%coco(13,2)=0.37748470203_DP
    gbasis(80)%expo(13,3)=1.2535869926_DP
    gbasis(80)%coco(13,3)=0.48342275878_DP

    ! 6s:
    gbasis(80)%angmom(14)=0
    gbasis(80)%expo(14,1)=2.4957_DP
    gbasis(80)%coco(14,1)=1.0_DP

    ! 6p:
    gbasis(80)%angmom(15)=1
    gbasis(80)%expo(15,1)=9.7729900000_DP
    gbasis(80)%coco(15,1)=0.21210940502_DP
    gbasis(80)%expo(15,2)=7.1690950000_DP
    gbasis(80)%coco(15,2)=-0.44026342153_DP
    gbasis(80)%expo(15,3)=1.7559567931_DP
    gbasis(80)%coco(15,3)=0.58406314464_DP



  end subroutine ngwf_data_cocos_and_expos61_80



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine ngwf_data_cocos_and_expos81_103(gbasis)

    !===============================================================!
    ! Initialisation of Gaussian basis set exponents and contraction!
    ! coefficients. All the data of the STO-3G basis set are        !
    ! included for each element of the periodic table as well as    !
    ! the one "polarisation" function of the 6-31G* set.            !
    !---------------------------------------------------------------!
    ! This subroutine was created by Chris-Kriton Skylaris on       !
    ! 31/03/2006 but the actual STO-3G and 6-31G* parameters        !
    ! which appear below were typed by Mark Robinson, in April 2007 !
    !===============================================================!

    implicit none

    type(GTO_SET), intent(inout) :: gbasis(1: pub_num_gtatoms)


    ! cks: all data in this subroutine was typed by Mark Robinson


    ! ====begin Tl ==========================================
    gbasis(81)%angmom(1) =0
    gbasis(81)%angmom(2) =0
    gbasis(81)%angmom(3) =1
    gbasis(81)%angmom(4) =0
    gbasis(81)%angmom(5) =1
    gbasis(81)%angmom(6) =2
    gbasis(81)%angmom(7) =0
    gbasis(81)%angmom(8) =1
    gbasis(81)%angmom(9) =2

    ! 5s:
    gbasis(81)%angmom(10)=0
    gbasis(81)%expo(10,1)=45.356924655_DP
    gbasis(81)%coco(10,1)=0.50938157397E-02_DP
    gbasis(81)%expo(10,2)=20.278843336_DP
    gbasis(81)%coco(10,2)=-0.14805704314_DP
    gbasis(81)%expo(10,3)=14.420612394_DP
    gbasis(81)%coco(10,3)=0.36768453564_DP

    ! 5p:
    gbasis(81)%angmom(11)=1
    gbasis(81)%expo(11,1)=9.0323058024_DP
    gbasis(81)%coco(11,1)=0.22738839085_DP
    gbasis(81)%expo(11,2)=7.5356419671_DP
    gbasis(81)%coco(11,2)=-0.35117975889_DP
    gbasis(81)%expo(11,3)=1.0521421998_DP
    gbasis(81)%coco(11,3)=0.25690085604_DP

    ! 6s:
    gbasis(81)%angmom(12)=0
    gbasis(81)%expo(12,1)=0.17643545481_DP
    gbasis(81)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(81)%angmom(13)=2
    gbasis(81)%expo(13,1)=9.4421940760_DP
    gbasis(81)%coco(13,1)=0.58250776986E-01_DP
    gbasis(81)%expo(13,2)=7.4410003381_DP
    gbasis(81)%coco(13,2)=-0.11789199301_DP
    gbasis(81)%expo(13,3)=1.9831617772_DP
    gbasis(81)%coco(13,3)=0.34232057636_DP

    ! 6p:
    gbasis(81)%angmom(14)=1
    gbasis(81)%expo(14,1)=0.16675605056_DP
    gbasis(81)%coco(14,1)=1.0000000000_DP

    ! 6d:
    gbasis(81)%angmom(15)=2
    gbasis(81)%expo(15,1)=0.15500000000_DP
    gbasis(81)%coco(15,1)=1.0000000000_DP


    ! ====begin Pb ==========================================
    gbasis(82)%angmom(1) =0
    gbasis(82)%angmom(2) =0
    gbasis(82)%angmom(3) =1
    gbasis(82)%angmom(4) =0
    gbasis(82)%angmom(5) =1
    gbasis(82)%angmom(6) =2
    gbasis(82)%angmom(7) =0
    gbasis(82)%angmom(8) =1
    gbasis(82)%angmom(9) =2

    ! 5s:
    gbasis(82)%angmom(10)=0
    gbasis(82)%expo(10,1)=55.727520333_DP
    gbasis(82)%coco(10,1)=0.29681153058E-02_DP
    gbasis(82)%expo(10,2)=20.302912827_DP
    gbasis(82)%coco(10,2)=-0.14352634362_DP
    gbasis(82)%expo(10,3)=14.558628750_DP
    gbasis(82)%coco(10,3)=0.37554188507_DP

    ! 5p:
    gbasis(82)%angmom(11)=1
    gbasis(82)%expo(11,1)=9.5531412199_DP
    gbasis(82)%coco(11,1)=0.20542943990_DP
    gbasis(82)%expo(11,2)=7.5327685129_DP
    gbasis(82)%coco(11,2)=-0.36407590058_DP
    gbasis(82)%expo(11,3)=1.1185539567_DP
    gbasis(82)%coco(11,3)=0.32217481686_DP

    ! 6s:
    gbasis(82)%angmom(12)=0
    gbasis(82)%expo(12,1)=0.20418933370_DP
    gbasis(82)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(82)%angmom(13)=2
    gbasis(82)%expo(13,1)=13.494505735_DP
    gbasis(82)%coco(13,1)=0.11571042371E-01_DP
    gbasis(82)%expo(13,2)=7.0609770411_DP
    gbasis(82)%coco(13,2)=-0.66746645076E-01_DP
    gbasis(82)%expo(13,3)=2.1581835375_DP
    gbasis(82)%coco(13,3)=0.33978760942_DP

    ! 6p:
    gbasis(82)%angmom(14)=1
    gbasis(82)%expo(14,1)=0.20027494801_DP
    gbasis(82)%coco(14,1)=1.0000000000_DP

    ! 6d:
    gbasis(82)%angmom(15)=2
    gbasis(82)%expo(15,1)=0.17500000000_DP
    gbasis(82)%coco(15,1)=1.0000000000_DP


    ! ====begin Bi ==========================================
    gbasis(83)%angmom(1) =0
    gbasis(83)%angmom(2) =0
    gbasis(83)%angmom(3) =1
    gbasis(83)%angmom(4) =0
    gbasis(83)%angmom(5) =1
    gbasis(83)%angmom(6) =2
    gbasis(83)%angmom(7) =0
    gbasis(83)%angmom(8) =1
    gbasis(83)%angmom(9) =2

    ! 5s:
    gbasis(83)%angmom(10)=0
    gbasis(83)%expo(10,1)=175.19632102_DP
    gbasis(83)%coco(10,1)=0.80278879028E-03_DP
    gbasis(83)%expo(10,2)=19.820147783_DP
    gbasis(83)%coco(10,2)=-0.14624461250_DP
    gbasis(83)%expo(10,3)=14.960110102_DP
    gbasis(83)%coco(10,3)=0.37435719885_DP

    ! 5p:
    gbasis(83)%angmom(11)=1
    gbasis(83)%expo(11,1)=10.672046293_DP
    gbasis(83)%coco(11,1)=0.14391258793_DP
    gbasis(83)%expo(11,2)=7.1992185915_DP
    gbasis(83)%coco(11,2)=-0.37994608603_DP
    gbasis(83)%expo(11,3)=1.5359697528_DP
    gbasis(83)%coco(11,3)=0.53268032065_DP

    ! 6s:
    gbasis(83)%angmom(12)=0
    gbasis(83)%expo(12,1)=0.23405391503_DP
    gbasis(83)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(83)%angmom(13)=2
    gbasis(83)%expo(13,1)=15.834273510_DP
    gbasis(83)%coco(13,1)=0.85062993231E-02_DP
    gbasis(83)%expo(13,2)=7.0956896087_DP
    gbasis(83)%coco(13,2)=-0.65954380016E-01_DP
    gbasis(83)%expo(13,3)=2.3363137007_DP
    gbasis(83)%coco(13,3)=0.33707163406_DP

    ! 6p:
    gbasis(83)%angmom(14)=1
    gbasis(83)%expo(14,1)=0.23769418731_DP
    gbasis(83)%coco(14,1)=1.0000000000_DP

    ! 6d:
    gbasis(83)%angmom(15)=2
    gbasis(83)%expo(15,1)=0.20000000000_DP
    gbasis(83)%coco(15,1)=1.0000000000_DP


    ! ====begin Po ==========================================
    gbasis(84)%angmom(1) =0
    gbasis(84)%angmom(2) =0
    gbasis(84)%angmom(3) =1
    gbasis(84)%angmom(4) =0
    gbasis(84)%angmom(5) =1
    gbasis(84)%angmom(6) =2
    gbasis(84)%angmom(7) =0
    gbasis(84)%angmom(8) =1
    gbasis(84)%angmom(9) =2

    ! 5s:
    gbasis(84)%angmom(10)=0
    gbasis(84)%expo(10,1)=278.13783916_DP
    gbasis(84)%coco(10,1)=0.91690497839E-03_DP
    gbasis(84)%expo(10,2)=19.659291105_DP
    gbasis(84)%coco(10,2)=-0.14713009283_DP
    gbasis(84)%expo(10,3)=15.016229647_DP
    gbasis(84)%coco(10,3)=0.38790413631_DP

    ! 5p:
    gbasis(84)%angmom(11)=1
    gbasis(84)%expo(11,1)=10.445945973_DP
    gbasis(84)%coco(11,1)=0.23631841405_DP
    gbasis(84)%expo(11,2)=8.0293101888_DP
    gbasis(84)%coco(11,2)=-0.45017818996_DP
    gbasis(84)%expo(11,3)=1.3588798023_DP
    gbasis(84)%coco(11,3)=0.44160334658_DP

    ! 6s:
    gbasis(84)%angmom(12)=0
    gbasis(84)%expo(12,1)=0.26329582133_DP
    gbasis(84)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(84)%angmom(13)=2
    gbasis(84)%expo(13,1)=18.292670150_DP
    gbasis(84)%coco(13,1)=0.72030963652E-02_DP
    gbasis(84)%expo(13,2)=7.1695578539_DP
    gbasis(84)%coco(13,2)=-0.66160544754E-01_DP
    gbasis(84)%expo(13,3)=2.5029225180_DP
    gbasis(84)%coco(13,3)=0.34363449378_DP

    ! 6p:
    gbasis(84)%angmom(14)=1
    gbasis(84)%expo(14,1)=0.27242942777_DP
    gbasis(84)%coco(14,1)=1.0000000000_DP

    ! 6d:
    gbasis(84)%angmom(15)=2
    gbasis(84)%expo(15,1)=0.21500000000_DP
    gbasis(84)%coco(15,1)=1.0000000000_DP


    ! ====begin At ==========================================
    gbasis(85)%angmom(1) =0
    gbasis(85)%angmom(2) =0
    gbasis(85)%angmom(3) =1
    gbasis(85)%angmom(4) =0
    gbasis(85)%angmom(5) =1
    gbasis(85)%angmom(6) =2
    gbasis(85)%angmom(7) =0
    gbasis(85)%angmom(8) =1
    gbasis(85)%angmom(9) =2

    ! 5s:
    gbasis(85)%angmom(10)=0
    gbasis(85)%expo(10,1)=303.23376176_DP
    gbasis(85)%coco(10,1)=0.63504969788E-03_DP
    gbasis(85)%expo(10,2)=20.289449240_DP
    gbasis(85)%coco(10,2)=-0.13823165350_DP
    gbasis(85)%expo(10,3)=15.383783126_DP
    gbasis(85)%coco(10,3)=0.38333193283_DP

    ! 5p:
    gbasis(85)%angmom(11)=1
    gbasis(85)%expo(11,1)=10.722023134_DP
    gbasis(85)%coco(11,1)=0.26194870008_DP
    gbasis(85)%expo(11,2)=8.4823338656_DP
    gbasis(85)%coco(11,2)=-0.46607228264_DP
    gbasis(85)%expo(11,3)=1.3474291294_DP
    gbasis(85)%coco(11,3)=0.42119853309_DP

    ! 6s:
    gbasis(85)%angmom(12)=0
    gbasis(85)%expo(12,1)=0.29269461143_DP
    gbasis(85)%coco(12,2)=1.0000000000_DP

    ! 5d:
    gbasis(85)%angmom(13)=2
    gbasis(85)%expo(13,1)=20.799819604_DP
    gbasis(85)%coco(13,1)=0.65227965381E-02_DP
    gbasis(85)%expo(13,2)=7.1664616770_DP
    gbasis(85)%coco(13,2)=-0.67780885217E-01_DP
    gbasis(85)%expo(13,3)=2.6934397048_DP
    gbasis(85)%coco(13,3)=0.33954991335_DP

    ! 6p:
    gbasis(85)%angmom(14)=1
    gbasis(85)%expo(14,1)=0.31036518470_DP
    gbasis(85)%coco(14,1)=1.0000000000_DP

    ! 6d:
    gbasis(85)%angmom(15)=2
    gbasis(85)%expo(15,1)=0.23500000000_DP
    gbasis(85)%coco(15,1)=1.0000000000_DP


    ! ====begin Rn ==========================================
    gbasis(86)%angmom(1) =0
    gbasis(86)%angmom(2) =0
    gbasis(86)%angmom(3) =1
    gbasis(86)%angmom(4) =0
    gbasis(86)%angmom(5) =1
    gbasis(86)%angmom(6) =2
    gbasis(86)%angmom(7) =0
    gbasis(86)%angmom(8) =1
    gbasis(86)%angmom(9) =2

    ! 5s:
    gbasis(86)%angmom(10)=0
    gbasis(86)%expo(10,1)=5187.7442468_DP
    gbasis(86)%coco(10,1)=0.16097034848E-03_DP
    gbasis(86)%expo(10,2)=768.14491207_DP
    gbasis(86)%coco(10,2)=0.10508214789E-02_DP
    gbasis(86)%expo(10,3)=159.35182065_DP
    gbasis(86)%coco(10,3)=0.27375438093E-02_DP

    ! 5p:
    gbasis(86)%angmom(11)=1
    gbasis(86)%expo(11,1)=194.68918723_DP
    gbasis(86)%coco(11,1)=0.26265814046E-03_DP
    gbasis(86)%expo(11,2)=11.210013525_DP
    gbasis(86)%coco(11,2)=0.19952302621_DP
    gbasis(86)%expo(11,3)=7.5298316489_DP
    gbasis(86)%coco(11,3)=-0.45994479561_DP

    ! 6s:
    gbasis(86)%angmom(12)=0
    gbasis(86)%expo(12,1)=2.0850437166_DP
    gbasis(86)%coco(12,1)=1.0000000000_DP

    ! 5d:
    gbasis(86)%angmom(13)=2
    gbasis(86)%expo(13,1)=75.302723665_DP
    gbasis(86)%coco(13,1)=0.73590238429E-03_DP
    gbasis(86)%expo(13,2)=17.965913830_DP
    gbasis(86)%coco(13,2)=0.89490827389E-02_DP
    gbasis(86)%expo(13,3)=7.3741298033_DP
    gbasis(86)%coco(13,3)=-0.76641663693E-01_DP

    ! 6p:
    gbasis(86)%angmom(14)=1
    gbasis(86)%expo(14,1)=0.47444492594_DP
    gbasis(86)%coco(14,1)=1.0000000000_DP

    ! 6d:
    gbasis(86)%angmom(15)=2
    gbasis(86)%expo(15,1)=0.46860482223_DP
    gbasis(86)%coco(15,1)=1.0000000000_DP

    ! 6f:
    gbasis(86)%angmom(16)=3
    gbasis(86)%expo(16,1)=2.1000000000_DP
    gbasis(86)%coco(16,1)=1.0000000000_DP


    ! ====begin Fr ==========================================
    gbasis(87)%angmom(1) =0
    gbasis(87)%angmom(2) =0
    gbasis(87)%angmom(3) =1
    gbasis(87)%angmom(4) =0
    gbasis(87)%angmom(5) =1
    gbasis(87)%angmom(6) =2
    gbasis(87)%angmom(7) =0
    gbasis(87)%angmom(8) =1
    gbasis(87)%angmom(9) =2
    gbasis(87)%angmom(10)=0
    gbasis(87)%angmom(11)=1
    gbasis(87)%angmom(12)=3
    gbasis(87)%angmom(13)=0

    ! 5d:
    gbasis(87)%angmom(14)=2
    gbasis(87)%expo(14,1)=0.8220000_DP
    gbasis(87)%coco(14,1)=1.0000000_DP
    gbasis(87)%expo(14,2)=0.6784000_DP
    gbasis(87)%coco(14,2)=1.0000000_DP
    gbasis(87)%expo(14,3)=0.1289000_DP
    gbasis(87)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(87)%angmom(15)=1
    gbasis(87)%expo(15,1)=0.9503000_DP
    gbasis(87)%coco(15,1)=1.0000000_DP
    gbasis(87)%expo(15,2)=0.4237000_DP
    gbasis(87)%coco(15,2)=1.0000000_DP
    gbasis(87)%expo(15,3)=0.1567000_DP
    gbasis(87)%coco(15,3)=1.0000000_DP

    ! 7s:
    gbasis(87)%angmom(16)=0
    gbasis(87)%expo(16,1)=0.7475000_DP
    gbasis(87)%coco(16,1)=1.0000000_DP
    gbasis(87)%expo(16,2)=0.5133000_DP
    gbasis(87)%coco(16,2)=1.0000000_DP
    gbasis(87)%expo(16,3)=0.2220000_DP
    gbasis(87)%coco(16,3)=1.0000000_DP


    ! ====begin Ra ==========================================
    gbasis(88)%angmom(1) =0
    gbasis(88)%angmom(2) =0
    gbasis(88)%angmom(3) =1
    gbasis(88)%angmom(4) =0
    gbasis(88)%angmom(5) =1
    gbasis(88)%angmom(6) =2
    gbasis(88)%angmom(7) =0
    gbasis(88)%angmom(8) =1
    gbasis(88)%angmom(9) =2
    gbasis(88)%angmom(10)=0
    gbasis(88)%angmom(11)=1
    gbasis(88)%angmom(12)=3
    gbasis(88)%angmom(13)=0

    ! 5d:
    gbasis(88)%angmom(14)=2
    gbasis(88)%expo(14,1)=1.6150000_DP
    gbasis(88)%coco(14,1)=1.0000000_DP
    gbasis(88)%expo(14,2)=0.4550000_DP
    gbasis(88)%coco(14,2)=1.0000000_DP
    gbasis(88)%expo(14,3)=0.1690000_DP
    gbasis(88)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(88)%angmom(15)=1
    gbasis(88)%expo(15,1)=1.1518000_DP
    gbasis(88)%coco(15,1)=1.0000000_DP
    gbasis(88)%expo(15,2)=0.4229000_DP
    gbasis(88)%coco(15,2)=1.0000000_DP
    gbasis(88)%expo(15,3)=0.1640000_DP
    gbasis(88)%coco(15,3)=1.0000000_DP

    ! 7s:
    gbasis(88)%angmom(16)=0
    gbasis(88)%expo(16,1)=0.8573000_DP
    gbasis(88)%coco(16,1)=1.0000000_DP
    gbasis(88)%expo(16,2)=0.5490000_DP
    gbasis(88)%coco(16,2)=1.0000000_DP
    gbasis(88)%expo(16,3)=0.1998000_DP
    gbasis(88)%coco(16,3)=1.0000000_DP


    ! ====begin Ac ==========================================
    gbasis(89)%angmom(1) =0
    gbasis(89)%angmom(2) =0
    gbasis(89)%angmom(3) =1
    gbasis(89)%angmom(4) =0
    gbasis(89)%angmom(5) =1
    gbasis(89)%angmom(6) =2
    gbasis(89)%angmom(7) =0
    gbasis(89)%angmom(8) =1
    gbasis(89)%angmom(9) =2
    gbasis(89)%angmom(10)=0
    gbasis(89)%angmom(11)=1
    gbasis(89)%angmom(12)=0

    ! 4f:
    gbasis(89)%angmom(13)=3
    gbasis(89)%expo(13,1)=0.1110000_DP
    gbasis(89)%coco(13,1)=1.0000000_DP
    gbasis(89)%expo(13,2)=0.0239000_DP
    gbasis(89)%coco(13,2)=1.0000000_DP
    gbasis(89)%expo(13,3)=0.0080000_DP
    gbasis(89)%coco(13,3)=1.0000000_DP

    ! 5d:
    gbasis(89)%angmom(14)=2
    gbasis(89)%expo(14,1)=1.3667000_DP
    gbasis(89)%coco(14,1)=1.0000000_DP
    gbasis(89)%expo(14,2)=0.3422000_DP
    gbasis(89)%coco(14,2)=1.0000000_DP
    gbasis(89)%expo(14,3)=0.1329000_DP
    gbasis(89)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(89)%angmom(15)=1
    gbasis(89)%expo(15,1)=0.9402000_DP
    gbasis(89)%coco(15,1)=1.0000000_DP
    gbasis(89)%expo(15,2)=0.6756000_DP
    gbasis(89)%coco(15,2)=1.0000000_DP
    gbasis(89)%expo(15,3)=0.2585000_DP
    gbasis(89)%coco(15,3)=1.0000000_DP

    ! 7s:
    gbasis(89)%angmom(16)=0
    gbasis(89)%expo(16,1)=0.9636000_DP
    gbasis(89)%coco(16,1)=1.0000000_DP
    gbasis(89)%expo(16,2)=0.5266000_DP
    gbasis(89)%coco(16,2)=1.0000000_DP
    gbasis(89)%expo(16,3)=0.1513000_DP
    gbasis(89)%coco(16,3)=1.0000000_DP


    ! ====begin Th ==========================================
    gbasis(90)%angmom(1) =0
    gbasis(90)%angmom(2) =0
    gbasis(90)%angmom(3) =1
    gbasis(90)%angmom(4) =0
    gbasis(90)%angmom(5) =1
    gbasis(90)%angmom(6) =2
    gbasis(90)%angmom(7) =0
    gbasis(90)%angmom(8) =1
    gbasis(90)%angmom(9) =2
    gbasis(90)%angmom(10)=0
    gbasis(90)%angmom(11)=1
    gbasis(90)%angmom(12)=2
    gbasis(90)%angmom(13)=0

    ! 4f:
    gbasis(90)%angmom(14)=3
    gbasis(90)%expo(14,1)=3.6871000_DP
    gbasis(90)%coco(14,1)=1.0000000_DP
    gbasis(90)%expo(14,2)=1.4405000_DP
    gbasis(90)%coco(14,2)=1.0000000_DP
    gbasis(90)%expo(14,3)=0.5334000_DP
    gbasis(90)%coco(14,3)=1.0000000_DP

    ! 6p:
    gbasis(90)%angmom(15)=1
    gbasis(90)%expo(15,1)=1.0955000_DP
    gbasis(90)%coco(15,1)=1.0000000_DP
    gbasis(90)%expo(15,2)=0.7035000_DP
    gbasis(90)%coco(15,2)=1.0000000_DP
    gbasis(90)%expo(15,3)=0.2839000_DP
    gbasis(90)%coco(15,3)=1.0000000_DP

    ! 7s:
    gbasis(90)%angmom(16)=0
    gbasis(90)%expo(16,1)=0.8559000_DP
    gbasis(90)%coco(16,1)=1.0000000_DP
    gbasis(90)%expo(16,2)=0.6788000_DP
    gbasis(90)%coco(16,2)=1.0000000_DP
    gbasis(90)%expo(16,3)=0.1687000_DP
    gbasis(90)%coco(16,3)=1.0000000_DP

    ! 6d:
    gbasis(90)%angmom(17)=2
    gbasis(90)%expo(17,1)=2.1232000_DP
    gbasis(90)%coco(17,1)=1.0000000_DP
    gbasis(90)%expo(17,2)=0.3507000_DP
    gbasis(90)%coco(17,2)=1.0000000_DP
    gbasis(90)%expo(17,3)=0.1351000_DP
    gbasis(90)%coco(17,3)=1.0000000_DP


    ! ====begin Pa ==========================================
    gbasis(91)%angmom(1) =0
    gbasis(91)%angmom(2) =0
    gbasis(91)%angmom(3) =1
    gbasis(91)%angmom(4) =0
    gbasis(91)%angmom(5) =1
    gbasis(91)%angmom(6) =2
    gbasis(91)%angmom(7) =0
    gbasis(91)%angmom(8) =1
    gbasis(91)%angmom(9) =2
    gbasis(91)%angmom(10)=0
    gbasis(91)%angmom(11)=1
    gbasis(91)%angmom(12)=3
    gbasis(91)%angmom(13)=2
    gbasis(91)%angmom(14)=0

    ! 6p:
    gbasis(91)%angmom(15)=1
    gbasis(91)%expo(15,1)=1.1821000_DP
    gbasis(91)%coco(15,1)=1.0000000_DP
    gbasis(91)%expo(15,2)=0.7366000_DP
    gbasis(91)%coco(15,2)=1.0000000_DP
    gbasis(91)%expo(15,3)=0.2963000_DP
    gbasis(91)%coco(15,3)=1.0000000_DP

    ! 7s:
    gbasis(91)%angmom(16)=0
    gbasis(91)%expo(16,1)=0.9652000_DP
    gbasis(91)%coco(16,1)=1.0000000_DP
    gbasis(91)%expo(16,2)=0.6098000_DP
    gbasis(91)%coco(16,2)=1.0000000_DP
    gbasis(91)%expo(16,3)=0.1636000_DP
    gbasis(91)%coco(16,3)=1.0000000_DP

    ! 5f:
    gbasis(91)%angmom(17)=3
    gbasis(91)%expo(17,1)=3.9144000_DP
    gbasis(91)%coco(17,1)=1.0000000_DP
    gbasis(91)%expo(17,2)=1.5449000_DP
    gbasis(91)%coco(17,2)=1.0000000_DP
    gbasis(91)%expo(17,3)=0.5737000_DP
    gbasis(91)%coco(17,3)=1.0000000_DP

    ! 6d:
    gbasis(91)%angmom(18)=2
    gbasis(91)%expo(18,1)=1.7027000_DP
    gbasis(91)%coco(18,1)=1.0000000_DP
    gbasis(91)%expo(18,2)=0.3748000_DP
    gbasis(91)%coco(18,2)=1.0000000_DP
    gbasis(91)%expo(18,3)=0.1417000_DP
    gbasis(91)%coco(18,3)=1.0000000_DP


    ! ====begin U  ==========================================
    gbasis(92)%angmom(1) =0
    gbasis(92)%angmom(2) =0
    gbasis(92)%angmom(3) =1
    gbasis(92)%angmom(4) =0
    gbasis(92)%angmom(5) =1
    gbasis(92)%angmom(6) =2
    gbasis(92)%angmom(7) =0
    gbasis(92)%angmom(8) =1
    gbasis(92)%angmom(9) =2
    gbasis(92)%angmom(10)=0
    gbasis(92)%angmom(11)=1
    gbasis(92)%angmom(12)=3
    gbasis(92)%angmom(13)=2
    gbasis(92)%angmom(14)=0

    ! 6p:
    gbasis(92)%angmom(15)=1
    gbasis(92)%expo(15,1)=1.4248000_DP
    gbasis(92)%coco(15,1)=1.0000000_DP
    gbasis(92)%expo(15,2)=0.6453000_DP
    gbasis(92)%coco(15,2)=1.0000000_DP
    gbasis(92)%expo(15,3)=0.2711000_DP
    gbasis(92)%coco(15,3)=1.0000000_DP

    ! 7s:
    gbasis(92)%angmom(16)=0
    gbasis(92)%expo(16,1)=0.9978000_DP
    gbasis(92)%coco(16,1)=1.0000000_DP
    gbasis(92)%expo(16,2)=0.7281000_DP
    gbasis(92)%coco(16,2)=1.0000000_DP
    gbasis(92)%expo(16,3)=0.2132000_DP
    gbasis(92)%coco(16,3)=1.0000000_DP

    ! 5f:
    gbasis(92)%angmom(17)=3
    gbasis(92)%expo(17,1)=4.3777000_DP
    gbasis(92)%coco(17,1)=1.0000000_DP
    gbasis(92)%expo(17,2)=1.7970000_DP
    gbasis(92)%coco(17,2)=1.0000000_DP
    gbasis(92)%expo(17,3)=0.7050000_DP
    gbasis(92)%coco(17,3)=1.0000000_DP

    ! 6d:
    gbasis(92)%angmom(18)=2
    gbasis(92)%expo(18,1)=2.1505000_DP
    gbasis(92)%coco(18,1)=1.0000000_DP
    gbasis(92)%expo(18,2)=0.3844000_DP
    gbasis(92)%coco(18,2)=1.0000000_DP
    gbasis(92)%expo(18,3)=0.1419000_DP
    gbasis(92)%coco(18,3)=1.0000000_DP


    ! ====begin Np ==========================================
    gbasis(93)%angmom(1) =0
    gbasis(93)%angmom(2) =0
    gbasis(93)%angmom(3) =1
    gbasis(93)%angmom(4) =0
    gbasis(93)%angmom(5) =1
    gbasis(93)%angmom(6) =2
    gbasis(93)%angmom(7) =0
    gbasis(93)%angmom(8) =1
    gbasis(93)%angmom(9) =2
    gbasis(93)%angmom(10)=0
    gbasis(93)%angmom(11)=1
    gbasis(93)%angmom(12)=3
    gbasis(93)%angmom(13)=2
    gbasis(93)%angmom(14)=0

    ! 6p:
    gbasis(93)%angmom(15)=1
    gbasis(93)%expo(15,1)=1.4945000_DP
    gbasis(93)%coco(15,1)=1.0000000_DP
    gbasis(93)%expo(15,2)=0.6844000_DP
    gbasis(93)%coco(15,2)=1.0000000_DP
    gbasis(93)%expo(15,3)=0.2952000_DP
    gbasis(93)%coco(15,3)=1.0000000_DP

    ! 7s:
    gbasis(93)%angmom(16)=0
    gbasis(93)%expo(16,1)=1.0012000_DP
    gbasis(93)%coco(16,1)=1.0000000_DP
    gbasis(93)%expo(16,2)=0.7687000_DP
    gbasis(93)%coco(16,2)=1.0000000_DP
    gbasis(93)%expo(16,3)=0.2064000_DP
    gbasis(93)%coco(16,3)=1.0000000_DP

    ! 5f:
    gbasis(93)%angmom(17)=3
    gbasis(93)%expo(17,1)=4.3677000_DP
    gbasis(93)%coco(17,1)=1.0000000_DP
    gbasis(93)%expo(17,2)=1.7779000_DP
    gbasis(93)%coco(17,2)=1.0000000_DP
    gbasis(93)%expo(17,3)=0.6802000_DP
    gbasis(93)%coco(17,3)=1.0000000_DP

    ! 6d:
    gbasis(93)%angmom(18)=2
    gbasis(93)%expo(18,1)=2.7871000_DP
    gbasis(93)%coco(18,1)=1.0000000_DP
    gbasis(93)%expo(18,2)=0.3930000_DP
    gbasis(93)%coco(18,2)=1.0000000_DP
    gbasis(93)%expo(18,3)=0.1438000_DP
    gbasis(93)%coco(18,3)=1.0000000_DP


    ! ====begin Pu ==========================================
    gbasis(94)%angmom(1) =0
    gbasis(94)%angmom(2) =0
    gbasis(94)%angmom(3) =1
    gbasis(94)%angmom(4) =0
    gbasis(94)%angmom(5) =1
    gbasis(94)%angmom(6) =2
    gbasis(94)%angmom(7) =0
    gbasis(94)%angmom(8) =1
    gbasis(94)%angmom(9) =2
    gbasis(94)%angmom(10)=0
    gbasis(94)%angmom(11)=1
    gbasis(94)%angmom(12)=3
    gbasis(94)%angmom(13)=2
    gbasis(94)%angmom(14)=0

    ! 6p:
    gbasis(94)%angmom(15)=1
    gbasis(94)%expo(15,1)=1.5611000_DP
    gbasis(94)%coco(15,1)=1.0000000_DP
    gbasis(94)%expo(15,2)=0.7242000_DP
    gbasis(94)%coco(15,2)=1.0000000_DP
    gbasis(94)%expo(15,3)=0.3081000_DP
    gbasis(94)%coco(15,3)=1.0000000_DP

    ! 7s:
    gbasis(94)%angmom(16)=0
    gbasis(94)%expo(16,1)=1.0817000_DP
    gbasis(94)%coco(16,1)=1.0000000_DP
    gbasis(94)%expo(16,2)=0.7457000_DP
    gbasis(94)%coco(16,2)=1.0000000_DP
    gbasis(94)%expo(16,3)=0.2115000_DP
    gbasis(94)%coco(16,3)=1.0000000_DP

    ! 5f:
    gbasis(94)%angmom(17)=3
    gbasis(94)%expo(17,1)=4.4495000_DP
    gbasis(94)%coco(17,1)=1.0000000_DP
    gbasis(94)%expo(17,2)=1.8334000_DP
    gbasis(94)%coco(17,2)=1.0000000_DP
    gbasis(94)%expo(17,3)=0.7127000_DP
    gbasis(94)%coco(17,3)=1.0000000_DP

    ! 6d:
    gbasis(94)%angmom(18)=2
    gbasis(94)%expo(18,1)=2.6903000_DP
    gbasis(94)%coco(18,1)=1.0000000_DP
    gbasis(94)%expo(18,2)=0.4088000_DP
    gbasis(94)%coco(18,2)=1.0000000_DP
    gbasis(94)%expo(18,3)=0.1434000_DP
    gbasis(94)%coco(18,3)=1.0000000_DP


    ! ====begin Am ==========================================
    gbasis(95)%angmom(1) =0
    gbasis(95)%angmom(2) =0
    gbasis(95)%angmom(3) =1
    gbasis(95)%angmom(4) =0
    gbasis(95)%angmom(5) =1
    gbasis(95)%angmom(6) =2
    gbasis(95)%angmom(7) =0
    gbasis(95)%angmom(8) =1
    gbasis(95)%angmom(9) =2
    gbasis(95)%angmom(10)=0
    gbasis(95)%angmom(11)=1
    gbasis(95)%angmom(12)=3
    gbasis(95)%angmom(13)=2
    gbasis(95)%angmom(14)=0

    ! 5d:
    gbasis(95)%angmom(15)=2
    gbasis(95)%expo(15,1)=11.5158890_DP
    gbasis(95)%coco(15,1)=1.0000000_DP
    gbasis(95)%expo(15,2)=1.9585446_DP
    gbasis(95)%coco(15,2)=1.0000000_DP
    gbasis(95)%expo(15,3)=1.3284298_DP
    gbasis(95)%coco(15,3)=1.0000000_DP

    ! 6p:
    gbasis(95)%angmom(16)=1
    gbasis(95)%expo(16,1)=0.0433588_DP
    gbasis(95)%coco(16,1)=1.0000000_DP
    gbasis(95)%expo(16,2)=0.0152267_DP
    gbasis(95)%coco(16,2)=1.0000000_DP

    ! 5f:
    gbasis(95)%angmom(17)=3
    gbasis(95)%expo(17,1)=5.3774034_DP
    gbasis(95)%coco(17,1)=1.0000000_DP
    gbasis(95)%expo(17,2)=2.4281957_DP
    gbasis(95)%coco(17,2)=1.0000000_DP
    gbasis(95)%expo(17,3)=1.1237942_DP
    gbasis(95)%coco(17,3)=1.0000000_DP


    ! ====begin Cm ==========================================
    gbasis(96)%angmom(1) =0
    gbasis(96)%angmom(2) =0
    gbasis(96)%angmom(3) =1
    gbasis(96)%angmom(4) =0
    gbasis(96)%angmom(5) =1
    gbasis(96)%angmom(6) =2
    gbasis(96)%angmom(7) =0
    gbasis(96)%angmom(8) =1
    gbasis(96)%angmom(9) =2
    gbasis(96)%angmom(10)=0
    gbasis(96)%angmom(11)=1
    gbasis(96)%angmom(12)=3
    gbasis(96)%angmom(13)=2
    gbasis(96)%angmom(14)=0

    ! 5d:
    gbasis(96)%angmom(15)=2
    gbasis(96)%expo(15,1)=11.0265260_DP
    gbasis(96)%coco(15,1)=1.0000000_DP
    gbasis(96)%expo(15,2)=2.3919811_DP
    gbasis(96)%coco(15,2)=1.0000000_DP
    gbasis(96)%expo(15,3)=1.3199383_DP
    gbasis(96)%coco(15,3)=1.0000000_DP

    ! 6p:
    gbasis(96)%angmom(16)=1
    gbasis(96)%expo(16,1)=0.0426398_DP
    gbasis(96)%coco(16,1)=1.0000000_DP
    gbasis(96)%expo(16,2)=0.0154094_DP
    gbasis(96)%coco(16,2)=1.0000000_DP

    ! 5f:
    gbasis(96)%angmom(17)=3
    gbasis(96)%expo(17,1)=5.3398412_DP
    gbasis(96)%coco(17,1)=1.0000000_DP
    gbasis(96)%expo(17,2)=2.4037829_DP
    gbasis(96)%coco(17,2)=1.0000000_DP
    gbasis(96)%expo(17,3)=1.1019327_DP
    gbasis(96)%coco(17,3)=1.0000000_DP


    ! ====begin Bk ==========================================
    gbasis(97)%angmom(1) =0
    gbasis(97)%angmom(2) =0
    gbasis(97)%angmom(3) =1
    gbasis(97)%angmom(4) =0
    gbasis(97)%angmom(5) =1
    gbasis(97)%angmom(6) =2
    gbasis(97)%angmom(7) =0
    gbasis(97)%angmom(8) =1
    gbasis(97)%angmom(9) =2
    gbasis(97)%angmom(10)=0
    gbasis(97)%angmom(11)=1
    gbasis(97)%angmom(12)=3
    gbasis(97)%angmom(13)=2
    gbasis(97)%angmom(14)=0

    ! 5d:
    gbasis(97)%angmom(15)=2
    gbasis(97)%expo(15,1)=16.7548840_DP
    gbasis(97)%coco(15,1)=1.0000000_DP
    gbasis(97)%expo(15,2)=2.2932401_DP
    gbasis(97)%coco(15,2)=1.0000000_DP
    gbasis(97)%expo(15,3)=1.2577831_DP
    gbasis(97)%coco(15,3)=1.0000000_DP

    ! 6p:
    gbasis(97)%angmom(16)=1
    gbasis(97)%expo(16,1)=0.0439084_DP
    gbasis(97)%coco(16,1)=1.0000000_DP
    gbasis(97)%expo(16,2)=0.0156209_DP
    gbasis(97)%coco(16,2)=1.0000000_DP

    ! 5f:
    gbasis(97)%angmom(17)=3
    gbasis(97)%expo(17,1)=5.7033444_DP
    gbasis(97)%coco(17,1)=1.0000000_DP
    gbasis(97)%expo(17,2)=2.5741569_DP
    gbasis(97)%coco(17,2)=1.0000000_DP
    gbasis(97)%expo(17,3)=1.1748984_DP
    gbasis(97)%coco(17,3)=1.0000000_DP


    ! ====begin Cf ==========================================
    gbasis(98)%angmom(1) =0
    gbasis(98)%angmom(2) =0
    gbasis(98)%angmom(3) =1
    gbasis(98)%angmom(4) =0
    gbasis(98)%angmom(5) =1
    gbasis(98)%angmom(6) =2
    gbasis(98)%angmom(7) =0
    gbasis(98)%angmom(8) =1
    gbasis(98)%angmom(9) =2
    gbasis(98)%angmom(10)=0
    gbasis(98)%angmom(11)=1
    gbasis(98)%angmom(12)=3
    gbasis(98)%angmom(13)=2
    gbasis(98)%angmom(14)=0

    ! 5d:
    gbasis(98)%angmom(15)=2
    gbasis(98)%expo(15,1)=20.7339090_DP
    gbasis(98)%coco(15,1)=1.0000000_DP
    gbasis(98)%expo(15,2)=2.0937082_DP
    gbasis(98)%coco(15,2)=1.0000000_DP
    gbasis(98)%expo(15,3)=1.5363060_DP
    gbasis(98)%coco(15,3)=1.0000000_DP

    ! 6p:
    gbasis(98)%angmom(16)=1
    gbasis(98)%expo(16,1)=0.0482689_DP
    gbasis(98)%coco(16,1)=1.0000000_DP
    gbasis(98)%expo(16,2)=0.0163072_DP
    gbasis(98)%coco(16,2)=1.0000000_DP

    ! 5f:
    gbasis(98)%angmom(17)=3
    gbasis(98)%expo(17,1)=6.0047068_DP
    gbasis(98)%coco(17,1)=1.0000000_DP
    gbasis(98)%expo(17,2)=2.7325189_DP
    gbasis(98)%coco(17,2)=1.0000000_DP
    gbasis(98)%expo(17,3)=1.3102299_DP
    gbasis(98)%coco(17,3)=1.0000000_DP


    ! ====begin Es ==========================================
    gbasis(99)%angmom(1) =0
    gbasis(99)%angmom(2) =0
    gbasis(99)%angmom(3) =1
    gbasis(99)%angmom(4) =0
    gbasis(99)%angmom(5) =1
    gbasis(99)%angmom(6) =2
    gbasis(99)%angmom(7) =0
    gbasis(99)%angmom(8) =1
    gbasis(99)%angmom(9) =2
    gbasis(99)%angmom(10)=0
    gbasis(99)%angmom(11)=1
    gbasis(99)%angmom(12)=3
    gbasis(99)%angmom(13)=2
    gbasis(99)%angmom(14)=0

    ! 5d:
    gbasis(99)%angmom(15)=2
    gbasis(99)%expo(15,1)=14.9607040_DP
    gbasis(99)%coco(15,1)=1.0000000_DP
    gbasis(99)%expo(15,2)=2.1642731_DP
    gbasis(99)%coco(15,2)=1.0000000_DP
    gbasis(99)%expo(15,3)=1.5640496_DP
    gbasis(99)%coco(15,3)=1.0000000_DP

    ! 6p:
    gbasis(99)%angmom(16)=1
    gbasis(99)%expo(16,1)=0.0489932_DP
    gbasis(99)%coco(16,1)=1.0000000_DP
    gbasis(99)%expo(16,2)=0.0167286_DP
    gbasis(99)%coco(16,2)=1.0000000_DP

    ! 5f:
    gbasis(99)%angmom(17)=3
    gbasis(99)%expo(17,1)=6.1699240_DP
    gbasis(99)%coco(17,1)=1.0000000_DP
    gbasis(99)%expo(17,2)=2.8361203_DP
    gbasis(99)%coco(17,2)=1.0000000_DP
    gbasis(99)%expo(17,3)=1.3227024_DP
    gbasis(99)%coco(17,3)=1.0000000_DP


    ! ====begin Fm ==========================================
    gbasis(100)%angmom(1) =0
    gbasis(100)%angmom(2) =0
    gbasis(100)%angmom(3) =1
    gbasis(100)%angmom(4) =0
    gbasis(100)%angmom(5) =1
    gbasis(100)%angmom(6) =2
    gbasis(100)%angmom(7) =0
    gbasis(100)%angmom(8) =1
    gbasis(100)%angmom(9) =2
    gbasis(100)%angmom(10)=0
    gbasis(100)%angmom(11)=1
    gbasis(100)%angmom(12)=3
    gbasis(100)%angmom(13)=2
    gbasis(100)%angmom(14)=0

    ! 5d:
    gbasis(100)%angmom(15)=2
    gbasis(100)%expo(15,1)=14.3727930_DP
    gbasis(100)%coco(15,1)=1.0000000_DP
    gbasis(100)%expo(15,2)=2.6279910_DP
    gbasis(100)%coco(15,2)=1.0000000_DP
    gbasis(100)%expo(15,3)=1.5581072_DP
    gbasis(100)%coco(15,3)=1.0000000_DP

    ! 6p:
    gbasis(100)%angmom(16)=1
    gbasis(100)%expo(16,1)=0.0486447_DP
    gbasis(100)%coco(16,1)=1.0000000_DP
    gbasis(100)%expo(16,2)=0.0166478_DP
    gbasis(100)%coco(16,2)=1.0000000_DP

    ! 5f:
    gbasis(100)%angmom(17)=3
    gbasis(100)%expo(17,1)=6.3362744_DP
    gbasis(100)%coco(17,1)=1.0000000_DP
    gbasis(100)%expo(17,2)=2.8970511_DP
    gbasis(100)%coco(17,2)=1.0000000_DP
    gbasis(100)%expo(17,3)=1.3538754_DP
    gbasis(100)%coco(17,3)=1.0000000_DP


    ! ====begin Md ==========================================
    gbasis(101)%angmom(1) =0
    gbasis(101)%angmom(2) =0
    gbasis(101)%angmom(3) =1
    gbasis(101)%angmom(4) =0
    gbasis(101)%angmom(5) =1
    gbasis(101)%angmom(6) =2
    gbasis(101)%angmom(7) =0
    gbasis(101)%angmom(8) =1
    gbasis(101)%angmom(9) =2
    gbasis(101)%angmom(10)=0
    gbasis(101)%angmom(11)=1
    gbasis(101)%angmom(12)=3
    gbasis(101)%angmom(13)=2
    gbasis(101)%angmom(14)=0

    ! 5d:
    gbasis(101)%angmom(15)=2
    gbasis(101)%expo(15,1)=22.0297190_DP
    gbasis(101)%coco(15,1)=1.0000000_DP
    gbasis(101)%expo(15,2)=2.1642731_DP
    gbasis(101)%coco(15,2)=1.0000000_DP
    gbasis(101)%expo(15,3)=1.5673249_DP
    gbasis(101)%coco(15,3)=1.0000000_DP

    ! 6p:
    gbasis(101)%angmom(16)=1
    gbasis(101)%expo(16,1)=0.0516183_DP
    gbasis(101)%coco(16,1)=1.0000000_DP
    gbasis(101)%expo(16,2)=0.0169365_DP
    gbasis(101)%coco(16,2)=1.0000000_DP

    ! 5f:
    gbasis(101)%angmom(17)=3
    gbasis(101)%expo(17,1)=6.8582197_DP
    gbasis(101)%coco(17,1)=1.0000000_DP
    gbasis(101)%expo(17,2)=3.1413322_DP
    gbasis(101)%coco(17,2)=1.0000000_DP
    gbasis(101)%expo(17,3)=1.4417434_DP
    gbasis(101)%coco(17,3)=1.0000000_DP


    ! ====begin No ==========================================
    gbasis(102)%angmom(1) =0
    gbasis(102)%angmom(2) =0
    gbasis(102)%angmom(3) =1
    gbasis(102)%angmom(4) =0
    gbasis(102)%angmom(5) =1
    gbasis(102)%angmom(6) =2
    gbasis(102)%angmom(7) =0
    gbasis(102)%angmom(8) =1
    gbasis(102)%angmom(9) =2
    gbasis(102)%angmom(10)=0
    gbasis(102)%angmom(11)=1
    gbasis(102)%angmom(12)=3
    gbasis(102)%angmom(13)=2
    gbasis(102)%angmom(14)=0

    ! 5d:
    gbasis(102)%angmom(15)=2
    gbasis(102)%expo(15,1)=20.2777110_DP
    gbasis(102)%coco(15,1)=1.0000000_DP
    gbasis(102)%expo(15,2)=2.1438433_DP
    gbasis(102)%coco(15,2)=1.0000000_DP
    gbasis(102)%expo(15,3)=1.8043526_DP
    gbasis(102)%coco(15,3)=1.0000000_DP

    ! 6p:
    gbasis(102)%angmom(16)=1
    gbasis(102)%expo(16,1)=0.0501206_DP
    gbasis(102)%coco(16,1)=1.0000000_DP
    gbasis(102)%expo(16,2)=0.0163136_DP
    gbasis(102)%coco(16,2)=1.0000000_DP

    ! 5f:
    gbasis(102)%angmom(17)=3
    gbasis(102)%expo(17,1)=7.3730206_DP
    gbasis(102)%coco(17,1)=1.0000000_DP
    gbasis(102)%expo(17,2)=3.2225813_DP
    gbasis(102)%coco(17,2)=1.0000000_DP
    gbasis(102)%expo(17,3)=1.3980402_DP
    gbasis(102)%coco(17,3)=1.0000000_DP


    ! ====begin Lr ==========================================
    gbasis(103)%angmom(1) =0
    gbasis(103)%angmom(2) =0
    gbasis(103)%angmom(3) =1
    gbasis(103)%angmom(4) =0
    gbasis(103)%angmom(5) =1
    gbasis(103)%angmom(6) =2
    gbasis(103)%angmom(7) =0
    gbasis(103)%angmom(8) =1
    gbasis(103)%angmom(9) =2
    gbasis(103)%angmom(10)=0
    gbasis(103)%angmom(11)=1
    gbasis(103)%angmom(12)=3
    gbasis(103)%angmom(13)=2
    gbasis(103)%angmom(14)=0

    ! 5d:
    gbasis(103)%angmom(15)=2
    gbasis(103)%expo(15,1)=18.9171780_DP
    gbasis(103)%coco(15,1)=1.0000000_DP
    gbasis(103)%expo(15,2)=2.5830972_DP
    gbasis(103)%coco(15,2)=1.0000000_DP
    gbasis(103)%expo(15,3)=1.7557520_DP
    gbasis(103)%coco(15,3)=1.0000000_DP

    ! 6p:
    gbasis(103)%angmom(16)=1
    gbasis(103)%expo(16,1)=0.0563985_DP
    gbasis(103)%coco(16,1)=1.0000000_DP
    gbasis(103)%expo(16,2)=0.0187504_DP
    gbasis(103)%coco(16,2)=1.0000000_DP

    ! 5f:
    gbasis(103)%angmom(17)=3
    gbasis(103)%expo(17,1)=7.5803955_DP
    gbasis(103)%coco(17,1)=1.0000000_DP
    gbasis(103)%expo(17,2)=3.4452695_DP
    gbasis(103)%coco(17,2)=1.0000000_DP
    gbasis(103)%expo(17,3)=1.5715476_DP
    gbasis(103)%coco(17,3)=1.0000000_DP


  end subroutine ngwf_data_cocos_and_expos81_103

#ifdef OLD_CONSTANTS
#undef DP
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end module ngwf_data



