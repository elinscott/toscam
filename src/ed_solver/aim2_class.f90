MODULE AIM2_class

   use masked_matrix_class, only: masked_matrix_type

   implicit none

   private

   !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   !$$ANDERSON IMPURITY MODEL CLASS$$
   !$$GENERIC QUADRATIC PART ONLY$$
   !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   !-------------------------------------!
   ! AIM2 = IMPURITY[QUADRATIC] + BATH   !
   !-------------------------------------!

   TYPE AIM2_type
      ! SPIN (IF NECESSARY)
      INTEGER :: spin = 0
      ! TOTAL NUMBER OF SITES
      INTEGER :: Ns = 0
      ! NUMBER OF SITES IN IMPURITY
      INTEGER :: Nc = 0
      ! NUMBER OF SITES IN BATH
      INTEGER :: Nb = 0
      ! SUPERCONDUCTING FLAG
      LOGICAL :: SUPER = .false.
      ! IMPURITY ENERGY  MATRIX
      TYPE(masked_matrix_type) :: Ec
      ! BATH     ENERGY  MATRICES
      TYPE(masked_matrix_type) :: Eb
      TYPE(masked_matrix_type) :: Vbc
      ! IMPURITY ORBITALS ORDERED IN FULL AIM ORBITAL BASIS
      INTEGER          :: IMPnorbs = 0
      INTEGER, POINTER :: IMPiorb(:) => NULL()
      ! BATH     ORBITALS ORDERED IN FULL AIM ORBITAL BASIS
      INTEGER          :: BATHnorbs = 0
      INTEGER, POINTER :: BATHiorb(:) => NULL()
   END TYPE

   public :: AIM2_type
   public :: delete_AIM2
   public :: new_AIM2

contains

   subroutine new_AIM2(AIM2, AIM, SPIN)

      use aim_class, only: aim_type
      use bath_class, only: nambu_eb, nambu_vbc
      use genvar, only: DBL
      use impurity_class, only: nambu_ec
      use masked_matrix_class, only: new_masked_matrix

      implicit none

      TYPE(AIM2_type), INTENT(INOUT) :: AIM2
      TYPE(AIM_type), INTENT(IN)     :: AIM
      INTEGER, OPTIONAL, INTENT(IN)  :: SPIN
      INTEGER :: spin_, shift_spin, site
      INTEGER :: BATHbds(2), IMPbds(2)

      CALL delete_AIM2(AIM2)

      AIM2%Nc = AIM%Nc
      AIM2%Nb = AIM%bath%Nb
      AIM2%Ns = AIM%Ns
      AIM2%SUPER = AIM%bath%SUPER

      ! NUMBER OF ORBITALS IN IMPURITY/BATH

      IF (AIM2%SUPER) THEN ! SUPERCONDUCTING : Sz sector
         AIM2%IMPnorbs = AIM2%Nc*2
         AIM2%BATHnorbs = AIM2%Nb*2
      ELSE               ! NORMAL : (npart, spin) sector
         AIM2%IMPnorbs = AIM2%Nc
         AIM2%BATHnorbs = AIM2%Nb
      END IF

      ! ORDER IMPURITY/BATH ORBITALS IN THE FULL BASIS

      IF (AIM2%SUPER) THEN
         ! IMPURITY ORBITALS IN NAMBU NOTATION: spin up first
         ! WE KEEP spin up FIRST IN FULL AIM BASIS, BUT WE START WE BATH
         ! ORBITALS
         ALLOCATE (AIM2%IMPiorb(AIM2%IMPnorbs))
         AIM2%IMPiorb = (/AIM%IMPiorb(:, 1), AIM%IMPiorb(:, 2)/)
         ! BATH ORBITALS IN NAMBU ORDER
         ALLOCATE (AIM2%BATHiorb(AIM2%BATHnorbs))
         AIM2%BATHiorb = (/AIM%BATHiorb(:, 1), AIM%BATHiorb(:, 2)/)
      ELSE
         ! FOR A GIVEN SPIN WE JUST START WITH BATH ORBITALS
         ALLOCATE (AIM2%IMPiorb(SIZE(AIM%IMPiorb, 1)))
         AIM2%IMPiorb = AIM%IMPiorb(:, 1)
         ! BATH ORBITALS FOR A GIVEN SPIN
         ALLOCATE (AIM2%BATHiorb(SIZE(AIM%BATHiorb, 1)))
         AIM2%BATHiorb = AIM%BATHiorb(:, 1)
      END IF

      ! GENERIC QUADRATIC ENERGY MATRICES

      IF (AIM2%SUPER) THEN ! NAMBU
         IF (PRESENT(SPIN)) STOP "ERROR IN new_AIM2: NO NEED TO SPECIFY SPIN &
              &DIRECTION!!"
         CALL Nambu_Ec(AIM2%Ec, AIM%impurity%Ec)
         CALL Nambu_Eb(AIM2%Eb, AIM%bath%Eb, AIM%bath%Pb)
         CALL Nambu_Vbc(AIM2%Vbc, AIM%bath%Vbc, AIM%bath%PVbc)
      ELSE ! FIXED SPIN
         IF (.NOT. PRESENT(SPIN)) STOP "ERROR IN new_AIM2: SHOULD SPECIFY A SPIN &
              &DIRECTION!!"
         AIM2%spin = SPIN
         ! these tmp matrices enforce tabulation of sector up & down in the
         ! NAMBU basis
         BATHbds = (/(SPIN - 1)*AIM2%Nb + 1, (SPIN - 1)*AIM2%Nb + AIM2%Nb/)
         IMPbds = (/(SPIN - 1)*AIM2%Nc + 1, (SPIN - 1)*AIM2%Nc + AIM2%Nc/)
         CALL new_masked_matrix(AIM2%Ec, AIM%impurity%Ec(SPIN))
         CALL new_masked_matrix(AIM2%Eb, AIM%bath%Eb(SPIN))
         CALL new_masked_matrix(AIM2%Vbc, AIM%bath%Vbc(SPIN))
      END IF

      ! WARNING: HERE THE MASKS STAND FOR NON-ZERO MATRIX ELEMENT ... AVOIDS
      ! UNNECESSARY TABULATIONS

      AIM2%Ec%rc%MASK%mat = .true.
      WHERE (AIM2%Ec%rc%mat == 0.0_DBL)
      AIM2%Ec%rc%MASK%mat = .false.
      END WHERE
      AIM2%Eb%rc%MASK%mat = .true.
      WHERE (AIM2%Eb%rc%mat == 0.0_DBL)
      AIM2%Eb%rc%MASK%mat = .false.
      END WHERE
      AIM2%Vbc%rc%MASK%mat = .true.
      WHERE (AIM2%Vbc%rc%mat == 0.0_DBL)
      AIM2%Vbc%rc%MASK%mat = .false.
      END WHERE

   end subroutine

   subroutine copy_AIM2(AIM2OUT, AIM2IN)

      ! CREATE AIM2OUT BY COPYING EXISTING AIM2IN

      use masked_matrix_class, only: copy_masked_matrix

      implicit none

      TYPE(AIM2_type), INTENT(INOUT) :: AIM2OUT
      TYPE(AIM2_type), INTENT(IN)    :: AIM2IN

      AIM2OUT%Nc = AIM2IN%Nc
      AIM2OUT%Nb = AIM2IN%Nb
      AIM2OUT%Ns = AIM2IN%Ns
      AIM2OUT%IMPnorbs = AIM2IN%IMPnorbs
      AIM2OUT%BATHnorbs = AIM2IN%BATHnorbs
      AIM2OUT%IMPiorb = AIM2IN%IMPiorb
      AIM2OUT%BATHiorb = AIM2IN%BATHiorb
      AIM2OUT%SUPER = AIM2IN%SUPER
      CALL copy_masked_matrix(AIM2OUT%Ec, AIM2IN%Ec)
      CALL copy_masked_matrix(AIM2OUT%Eb, AIM2IN%Eb)
      CALL copy_masked_matrix(AIM2OUT%Vbc, AIM2IN%Vbc)
   end subroutine

   subroutine delete_AIM2(AIM2)

      ! CREATE AIM2OUT BY COPYING EXISTING AIM2IN

      use masked_matrix_class, only: delete_masked_matrix

      implicit none

      TYPE(AIM2_type), INTENT(INOUT) :: AIM2

      IF (ASSOCIATED(AIM2%IMPiorb)) DEALLOCATE (AIM2%IMPiorb)
      IF (ASSOCIATED(AIM2%BATHiorb)) DEALLOCATE (AIM2%BATHiorb)
      CALL delete_masked_matrix(AIM2%Ec)
      CALL delete_masked_matrix(AIM2%Eb)
      CALL delete_masked_matrix(AIM2%Vbc)
   end subroutine

end module
