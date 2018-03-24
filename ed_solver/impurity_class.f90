MODULE impurity_class

   !$$$$$$$$$$$$$$$$$$$$
   !$$ IMPURITY CLASS $$
   !$$$$$$$$$$$$$$$$$$$$

   ! use common_def,              only: dump_message
   use genvar,                  only: DBL! , cspin, log_unit, rank
   use globalvar_ed_solver,     only: istati
   use linalg,                  only: arrayi
   ! use matrix,                  only: diag, write_array
   use masked_matrix_class, only: masked_matrix_type, masked_real_matrix_type!, &
   !      new_masked_matrix, write_masked_matrix
   ! use masked_matrix_class_mod, only: new_masked_real_matrix
   !use latticerout   , only : web, T1_T2_connect_unitcells
   !use quantum_hamilt, only : Hamiltonian

   implicit none

   private

   TYPE symmetry
      ! ======= groupe ponctuel ======== !
      integer, dimension(:, :), allocatable    :: rot, trans, alltrans, flip, &
                                                  invtrans
      real(8), dimension(:, :), allocatable    :: ttrans
      integer                                  :: nmax, nnmax !max numb of sym.
      real(8)                                  :: BZVEC(3), BZVEC2(3)
      complex(8), dimension(:), allocatable    :: rotchar, flipchar, transchar
      complex(8)                               :: chirot, chiflip
      integer                                  :: ntrans, nrot, nflip
      real(8), dimension(:, :, :), allocatable :: rotmat
      real(8)                                  ::  centerrot(3), diagonale(3)
      integer                                  ::  flag(20), flageff(20)
      integer, dimension(:, :, :), allocatable ::  symtab
      integer, dimension(:, :), allocatable    ::  sympoint, whichk
      type(arrayi), dimension(:), allocatable  ::  tabid
      integer                                  ::  tabidli
      integer, dimension(:), allocatable       ::  tabidcol
      ! ======== groupe ponctuel ======== !
   END TYPE

   TYPE unitcell
      integer                                   :: q1, q1small, q1mag, ncor, npk
      real(8)                                   :: band_plot_shift
      real(8), dimension(:,:), allocatable      :: insidepos, pk
      real(8), dimension(:,:,:), allocatable    :: clock
      integer, dimension(:,:), allocatable      :: half
      integer, dimension(:), allocatable        :: nneigh
      integer, dimension(:,:,:), allocatable    :: invlink
      real(8)                                   :: a(3), b(3), c(3), as(3), &
                                                   bs(3), cs(3), VBZ
      integer, dimension(:), allocatable        :: smallsite, smallsite2, &
                                                   orbitals
      integer, dimension(:,:), allocatable      :: phase2, phase3
      integer                                   :: struct_nm, struct_nlink
      integer, dimension(:,:), allocatable      :: struct
      logical                                   :: FLAG_MATRIX_ELEMENTS
      complex(8), dimension(:,:,:), allocatable :: matkxky, matkxkyshift
      real(8), dimension(:,:), allocatable      :: dispersion
      complex(8), dimension(:,:,:), allocatable :: eigenstates
      integer                                   :: k_zero
      real(8)                                   :: reciproc_base(3, 3), P(3, 3), &
                                                   T1(3), T2(3), T3(3)
      real(8)                                   :: reciproc_boundary(3, 3)
      real(8), dimension(:), allocatable        :: xr, yr, zr, xrpp, yrpp, zrpp
      integer                                   :: nL, nLz
      integer                                   :: nfermiNnumb = 2000
      real(8), dimension(:,:), allocatable      :: nfermiN, Akwfree
      real(8), dimension(:), allocatable        :: polygonx, polygony
      integer                                   :: Npoly
      integer, dimension(:,:), allocatable      :: site_connected_to_bath
      integer, dimension(:,:,:), allocatable    :: site_connected_to_cor
      logical, dimension(:), allocatable        :: sitecor
      integer, dimension(:), allocatable        :: reduced_space
      integer                                   :: mbath = 0, mdata = 0, &
                                                   nmanybody = 0
      logical                                   :: FLAG_ONLY_NN
      integer                                   :: FLAG_MAG_FIELD
      integer                                   :: FLAG_HORSDIAG_SIGMA
   END TYPE


   TYPE web

      type(unitcell)                            :: cell
      type(symmetry)                            :: sym

      !web:
      logical                                   :: symmetries
      real(8), dimension(:,:), allocatable      :: nlink_
      integer                                   :: mini, maxi, Ncu, maxconn
      integer, dimension(:,:,:,:), allocatable  :: struct_link
      real(8), dimension(:,:,:), allocatable    :: struct_pos
      integer                                   :: plaquettecentre(100), plaquetteN
      integer, dimension(:), allocatable        :: latticemapi
      real(8)                                   :: P(3, 3), teta(2)
      real(8)                                   :: T1(3), T2(3), T3(3)
      logical                                   :: FLAG_NO_SYM_AT_ALL
      logical                                   :: FLAG_OPEN_BC, &
                                                   FLAG_OPEN_BUT_PERIODIC_LINKED, &
                                                   FLAG_UNIT_CELL_INSIDE, &
                                                   FLAG_PLOT_3D
      real(8), dimension(:,:), allocatable      :: x, xma
      real(8)                                   :: randvec(3)
      integer, dimension(:), allocatable        :: site, site2, site3
      integer                                   :: open1, open2, open3, nL, nLz
      complex(8), dimension(:,:), allocatable   :: phase, phase4
      integer                                   :: centresite
      integer, dimension(:,:), allocatable      :: vecin
      real(8), dimension(:), allocatable        :: latticemap
      real(8), dimension(:,:), allocatable      :: distbord, periodist, &
                                                   periodist_open_bc
      real(8), dimension(:), allocatable        :: distances, angles
      real(8), dimension(:,:), allocatable      :: distancesp
      integer, dimension(:), allocatable        :: bordure, whichbord
      integer, dimension(:,:), allocatable      :: nombrevoisin
      integer                                   :: maxd, maxd_real, maxvoisins, &
                                                   bordcount
      integer, dimension(:,:), allocatable      :: maxdistsite
      integer, dimension(:,:,:), allocatable    :: longvoisin
      real(8)                                   :: centre(3)
      integer, dimension(:,:), allocatable      :: ineigh, cadran, direc, &
                                                   centerCu, centreangle
      integer, dimension(:,:), allocatable      :: links
      integer                                   :: nlinks, N
      integer, dimension(:), allocatable        :: conv_unitcell_site_, &
                                                   centredist, k_inv
      logical, dimension(:), allocatable        :: major_site
      integer                                   :: normalisation, nlink_cell
      complex(8), dimension(:,:), allocatable   :: unitary_transform

      !--------------------------------!
      ! copy of the unitcell variables !
      !--------------------------------!

      real(8), dimension(:,:), allocatable      :: dispersion
      real(8)                                   :: reciproc_base(3, 3)
      real(8)                                   :: reciproc_boundary(3, 3)
      real(8), dimension(:), allocatable        :: xr, yr, zr, xrpp, yrpp, zrpp
      complex(8), dimension(:,:,:), allocatable :: matkxky, matkxkyshift
      integer                                   :: q1, q1small, q1mag
      real(8), dimension(:,:), allocatable      :: insidepos
      real(8), dimension(:,:,:), allocatable    :: clock
      integer, dimension(:,:), allocatable      :: half
      integer, dimension(:), allocatable        :: nneigh
      real(8)                                   :: a(3), b(3), c(3)
      integer, dimension(:), allocatable        :: smallsite, smallsite2
      integer, dimension(:,:), allocatable      :: phase2, phase3
      integer                                   :: struct_nm, struct_nlink
      integer, dimension(:,:), allocatable      :: struct

   END TYPE


   TYPE Hamiltonian
      character(22)                                 :: hamiltonian_title
      integer                                       :: RANGE_TETA = 1
      logical                                       :: spin_sector_not_coupled, &
                                                       NO_DC = .false., &
                                                       DC_NO_SAME_ORBITALS = &
                                                       .false.
      integer                                       :: ordertype, ncor, mbath, &
                                                       mdata ! 1 = para, 2 =
      ! spinupdn, 3 = bcs, 4 = currents
      character(22)                                 :: orderlabel(100)
      complex(8), dimension(:,:), allocatable       :: teta, delta, tetadn, &
                                                       deltap, sigma_hf
      complex(8), dimension(:,:,:,:,:), allocatable :: teta_, tetadn_
      complex(8), dimension(:,:), allocatable       :: delta_mat, teta_mat, &
                                                       teta_mat_dn, teta_mat0
      complex(8), dimension(:,:), allocatable       :: delta_mat_p
      real(8), dimension(:,:), allocatable          :: rrrrsign, rrrrrsign
      real(8), dimension(:), allocatable            :: rrsign, rsign
      real(8), dimension(:), allocatable            :: eps_mat
      logical, dimension(:), allocatable            :: full_proj
      real(8), dimension(:), allocatable            :: eps, dU
      complex(8), dimension(:,:), allocatable       :: epsk
      integer, dimension(:,:), allocatable          :: type
      real(8), dimension(:,:), allocatable          :: Vrep, Jterm, field, &
                                                       field_cell
      integer                                       :: subE, ntype, subDiag, q1
      logical, dimension(:), allocatable            :: diagornot, forbidden
      character(20), dimension(:), allocatable      :: label
      integer                                       :: nH
      character(22)                                 :: labH(100)
      real(8)                                       :: min(100), max(100), hund, &
                                                       disorder, a_Angstrom, &
                                                       b_Angstrom, c_Angstrom, &
                                                       dist_plane
      real(8)                                       :: cor1, cor2, cor3, cor4, &
                                                       quanta, temperature
      real(8)                                       :: rrrsign
   end type

   TYPE(masked_matrix_type), PUBLIC, SAVE :: Eccc

   TYPE impurity_type
      ! IMPURITY TYPE
      ! NUMBER OF SITES
      INTEGER :: Nc      = 0
      ! NUMBER OF 1-PARTICLE ORBITALS
      INTEGER :: norbs   = 0
      ! SIZE OF REDUCED HILBERT SPACE OF THE IMPURITY
      INTEGER :: nstates = 0
      ! RANK OF ORBITALS
      INTEGER, POINTER :: iorb(:, :) => NULL()
      ! IMPURITY QUADRATIC ENERGY
      TYPE(masked_matrix_type), POINTER :: Ec(:) => NULL() ! in (site, site)
      ! basis for a given spin
      ! IMPURITY QUARTIC ENERGY
      TYPE(masked_real_matrix_type) :: U  ! in (site, site) basis
   END TYPE


   public :: define_impurity
   public :: delete_impurity
   public :: hamiltonian
   public :: impurity_type
   public :: nambu_ec
   public :: symmetry
   public :: update_impurity
   public :: web

contains

   subroutine T1_T2_connect_unitcells(xx, i, j, n1, n2, n3, cadran)

      use matrix, only: decomposevec
      use linalg, only: cells

      implicit none

      type(web)         :: xx
      integer           :: i, j, k, l, m, n1, n2, n3
      real(8)           :: base(3, 3), Tdecomp(3), vv(3), aa, bb, cc
      logical, optional :: cadran

      if(.not.present(cadran))then
         base(1, :) = xx%T1
         base(2, :) = xx%T2
         base(3, :) = xx%T3
         vv        = xx%xma(j, :) - xx%xma(i, :)
         call decomposevec(Tdecomp, vv, base)
         n1 = NINT(Tdecomp(1))
         n2 = NINT(Tdecomp(2))
         n3 = NINT(Tdecomp(3))
      else
         call cells(xx%cadran(i, j), aa, bb, cc)
         n1 = NINT(aa)
         n2 = NINT(bb)
         n3 = NINT(cc)
      endif

      return
   end subroutine

   subroutine new_impurity(impurity, Nc, IMASKE, IMASKU)

      use genvar,                  only: cspin
      use linalg,                  only: ramp
      use masked_matrix_class_mod, only: new_masked_real_matrix
      use masked_matrix_class,     only: new_masked_matrix

      implicit none

      TYPE(impurity_type), INTENT(INOUT) :: impurity
      INTEGER, INTENT(IN)                :: Nc
      INTEGER, OPTIONAL, INTENT(IN)      :: IMASKU(Nc, Nc), IMASKE(Nc, Nc, 2)
      INTEGER :: spin

      CALL delete_impurity(impurity)

      ! NUMBER OF SITES
      impurity%Nc    = Nc

      ! NUMBER OF 1-PARTICLE ORBITALS
      impurity%norbs = Nc * 2

      ! NUMBER OF IMPURITY STATES
      impurity%nstates = 2**impurity%norbs

      ! ORDERING OF ORBITALS WITH INCREASING RANK = |(site, up) > |(site,
      ! down) >

      IF(ASSOCIATED(impurity%iorb)) DEALLOCATE(impurity%iorb, STAT = istati)
      ALLOCATE(impurity%iorb(Nc, 2))

      CALL ramp(impurity%iorb(:, 1))

      impurity%iorb(:, 2) = impurity%iorb(:, 1) + Nc

      IF(PRESENT(IMASKE))THEN
         ! QUADRATIC ENERGY
         if(ASSOCIATED(impurity%Ec)) DEALLOCATE(impurity%Ec, STAT = istati)
         ALLOCATE(impurity%Ec(SIZE(IMASKE, 3)))
         DO spin = 1, SIZE(IMASKE, 3)
            CALL new_masked_matrix(impurity%Ec(spin), "Ec(sz = " // &
                 TRIM(cspin(spin)) // ")", Nc, Nc, IMASK = IMASKE(:, :, spin), &
                 IS_HERM = .true.)
         ENDDO
      ENDIF

      IF(PRESENT(IMASKU))THEN
         ! QUARTIC ENERGY
         CALL new_masked_real_matrix(impurity%U, "U", Nc, Nc, IMASK = IMASKU, &
              IS_HERM = .true.)
      ENDIF

   end subroutine

   subroutine delete_impurity(IMP)

      use masked_matrix_class, only: delete_masked_matrix

      implicit none

      TYPE(impurity_type), INTENT(INOUT) :: IMP
      INTEGER :: spin

      IF(ASSOCIATED(IMP%iorb)) DEALLOCATE(IMP%iorb, STAT = istati)
      IF(ASSOCIATED(IMP%Ec))THEN
         DO spin = 1, SIZE(IMP%Ec)
            CALL delete_masked_matrix(IMP%Ec(spin))
         ENDDO
         DEALLOCATE(IMP%Ec, STAT = istati)
      ENDIF

   end subroutine

   subroutine copy_impurity(IMPOUT, IMPIN)

      use masked_matrix_class_mod, only: copy_masked_real_matrix
      use masked_matrix_class,     only: copy_masked_matrix

      implicit none

      TYPE(impurity_type), INTENT(INOUT) :: IMPOUT
      TYPE(impurity_type), INTENT(IN)    :: IMPIN
      INTEGER :: spin


      IF(.NOT.ASSOCIATED(IMPIN%Ec)) STOP "ERROR IN copy_impurity: INPUT ISNT &
           &ALLOCATED!"
      IF(.NOT.ASSOCIATED(IMPOUT%Ec)) STOP "ERROR IN copy_impurity: OUTPUT ISNT &
           &ALLOCATED!"
      IMPOUT%Nc      = IMPIN%Nc
      IMPOUT%norbs   = IMPIN%norbs
      IMPOUT%nstates = IMPIN%nstates
      DO spin = 1, SIZE(IMPIN%Ec)
         CALL copy_masked_matrix(IMPOUT%Ec(spin), IMPIN%Ec(spin))
      ENDDO
      CALL copy_masked_real_matrix(IMPOUT%U, IMPIN%U)

   end subroutine

   subroutine define_impurity(impurity, mmu, impurity_, Himp, Eimp)

      ! $$ READ IMPURITY PARAMETERS

      use genvar,                  only: log_unit, rank
      use globalvar_ed_solver,     only: uumatrix
      use masked_matrix_class_mod, only: fill_masked_real_matrix, &
           new_masked_real_matrix, test_masked_real_matrix_symmetric
      use matrix,                  only: write_array

      implicit none

      TYPE(impurity_type), INTENT(INOUT) :: impurity
      type(web)              :: impurity_
      type(hamiltonian)      :: Himp
      integer                :: jj, ff, k, i, j, kkk, ijk, ki, kj, kki, kkj, ii
      REAL(DBL)              :: rval, mmu, val(impurity_%N**2)
      INTEGER                :: mu, Nc, spin, iind_, iind
      INTEGER, ALLOCATABLE   :: IMASKU(:,:)
#ifdef _complex
      COMPLEX(DBL), OPTIONAL :: Eimp(:,:)
#else
      REAL(DBL), OPTIONAL    :: Eimp(:,:)
#endif

      CALL delete_impurity(impurity)
      Nc = impurity_%N

      write(log_unit, *) '...... starting impurity problem with number of site &
           &: ', Nc
      CALL new_impurity(impurity, Nc)

      ! NON-QUADRATIC ENERGY

      if(allocated(IMASKU)) deallocate(IMASKU, STAT = istati)
      ALLOCATE(IMASKU(Nc, Nc))

      IMASKU = 0
      kkk = 0
      val = 0.d0

      do jj = 1, Nc
         kkk           = kkk + 1
         IMASKU(jj, jj) = kkk
         if(.not.allocated(UUmatrix))then
            val(kkk)      = Himp%dU(impurity_%site(jj))
         else
            val(kkk)      = UUmatrix(jj, jj)
         endif
      enddo

      kkk = jj-1

      do jj = 1, Nc
         if(.not.allocated(UUmatrix))then
            do i = 1, impurity_%nneigh(impurity_%site(jj))
               if(impurity_%cadran(jj, i) == 5)then
                  kkk                               = kkk + 1
                  if(kkk > size(val)) then
                     do j = 1, Nc
                        do ii = 1, impurity_%nneigh(impurity_%site(j))
                           if(impurity_%cadran(j, ii) == 5)then
                              write(log_unit, *) ' site, neighbor : ', j, &
                                   impurity_%ineigh(j, ii)
                           endif
                        enddo
                     enddo
                     write(*, *) 'size val   : ', size(val)
                     write(*, *) 'impurity%N : ', impurity_%N
                     write(*, *) 'Nc         : ', Nc
                     stop 'error build impurity Vrep matrix, too much &
                          &elements, a geometry problem?'
                  endif
                  IMASKU(jj, impurity_%ineigh(jj, i)) = kkk
                  val(kkk) = Himp%Vrep(impurity_%site(jj), i)
               endif
            enddo
         else
            do i = jj + 1, Nc
               kkk = kkk + 1
               if(kkk > size(val)) then
                  write(*, *) 'stop error in building impurity'
                  stop
               endif
               IMASKU(jj, i) = kkk
               val(kkk)     = UUmatrix(jj, i)
            enddo
         endif
      enddo

      CALL new_masked_real_matrix(impurity%U, "U", Nc, Nc, IMASK = IMASKU, &
           IS_HERM = .true.)

      write(145 + rank, *) ' ===== DIAGONALIZING IMP WITH ONSITE REPULSION &
           &===== '
      DO iind = 1, impurity%U%MASK%nind
         rval = val(iind)
         CALL fill_masked_real_matrix(impurity%U, iind, rval)
      ENDDO
      call write_array(impurity%U%mat, ' Coulomb repulsion ', unit = 145 + &
           rank, short = .true.)
      write(145 + rank, *) ' ==================================================='


      ! TEST HERMITICITY
      write(log_unit, *) 'test hermiticity impurity U'
      CALL test_masked_real_matrix_symmetric(impurity%U)

      ! QUADRATIC ENERGY

      if(associated(impurity%Ec)) deallocate(impurity%Ec, STAT = istati)
      ALLOCATE(impurity%Ec(2))

      call update_impurity(impurity_%N, impurity, mmu, impurity_, Himp, Eimp = &
           Eimp)

      if(allocated(IMASKU)) deallocate(IMASKU, STAT = istati)

   end subroutine

   subroutine update_impurity(Nc, impurity, mmu, impurity_, Himp, Eimp)

      use genvar,              only: cspin, log_unit
      use masked_matrix_class, only: clean_redundant_imask, &
           delete_masked_matrix, fill_masked_matrix, new_masked_matrix, &
           test_masked_matrix_hermitic
      use matrix,              only: write_array

      implicit none

      TYPE(impurity_type), INTENT(INOUT) :: impurity
      integer                            :: Nc
      type(web)                          :: impurity_
      type(hamiltonian)                  :: Himp
      integer                            :: jj, ff, k, i, j
#ifdef _complex
      COMPLEX(DBL)                       :: val(Nc**2*4)
      COMPLEX(DBL), OPTIONAL             :: Eimp(:,:)
#else
      REAL(DBL)                          :: val(Nc**2*4)
      REAL(DBL), OPTIONAL                :: Eimp(:,:)
#endif
      REAL(DBL)                          :: rval, mmu
      INTEGER                            :: mu, spin, iind_, iind
      INTEGER                            :: IMASKE(Nc, Nc, 2)


      IMASKE = 0
      k = 0
      val = 0.d0

      if(.not.present(Eimp))then

         do i = 1, impurity_%N
            k = k + 1
            if(allocated(Himp%eps)) val(k) = -mmu + Himp%eps(impurity_%site(i))
            IMASKE(i, i, 1:2) = k
            if(allocated(impurity_%nneigh))then
               do j = 1, impurity_%nneigh(impurity_%site(i))
                  if(allocated(impurity_%ineigh) .and. impurity_%cadran(i, j) &
                       == 5)then
                     k = k + 1
                     val(k) = Himp%teta(impurity_%site(i), j)
                     IMASKE(i, impurity_%ineigh(i, j), 1:2) = k
                  endif
               enddo
            endif
         enddo

      else

         if(Nc /= impurity_%N)then
            write(*, *) 'something wrong in impurity class, size Nc'
            stop
         endif
         if(size(Eimp, 1) /= 2*impurity_%N)then
            write(*, *) 'something wrong in impurity class, size Eimp and Nc &
                 &do not match'
            write(*, *) 'shape Eimp : ', shape(Eimp)
            write(*, *) 'impurity%N : ', impurity_%N
            stop
         endif

         do i = 1, impurity_%N
            do j = 1, impurity_%N
               if(abs(Eimp(i, j)) > 1.d-7 .or. i == j)then
                  k = k + 1
                  val(k) =   Eimp(i, j)
                  if(i == j) val(k) =   val(k) - mmu
                  IMASKE(i, j, 1) = k
               endif
               if(abs(Eimp(j + Nc, i + Nc)) > 1.d-7 .or. i == j)then
                  k = k + 1
                  val(k) = - Eimp(j + Nc, i + Nc) ! Eimp in the supra form,
                  ! here we want the TB form
                  if(i == j) val(k) =   val(k) - mmu
                  IMASKE(i, j, 2) = k
               endif
            enddo
         enddo

      endif

      DO spin = 1, 2
         CALL new_masked_matrix(impurity%Ec(spin), "Ec(sz = " // &
              TRIM(cspin(spin)) // ")", Nc, Nc, IMASK = IMASKE(:, :, spin), &
              IS_HERM = .true.)
      ENDDO
      CALL clean_redundant_imask(impurity%Ec)

      DO spin = 1, SIZE(impurity%Ec)
         do i = 1, k
            CALL fill_masked_matrix(impurity%Ec(spin), i, val(i))
         enddo
         call write_array( impurity%Ec(spin)%rc%mat, ' Ec(spin) ', unit = &
              log_unit, short = .true. )
         call write_array( IMASKE(:, :, spin), ' mask spin ', unit = log_unit )
      enddo

      if(present(Eimp))then
         call write_array( Eimp, ' real Eimp ', unit = log_unit, short = .true.)
      endif

      ! TEST HERMITICITY
      DO spin = 1, SIZE(impurity%Ec)
         write(log_unit, *) 'test hermiticity Ec spin ', spin
         CALL test_masked_matrix_hermitic(impurity%Ec(spin))
      ENDDO

      CALL delete_masked_matrix(Eccc)
      CALL Nambu_Ec(Eccc, impurity%Ec)

   end subroutine

   subroutine write_impurity(impurity, UNIT)

      use common_def,              only: dump_message
      use genvar,                  only: log_unit
      use masked_matrix_class,     only: write_masked_matrix
      use masked_matrix_class_mod, only: write_masked_real_matrix

      implicit none

      TYPE(impurity_type), INTENT(IN) :: impurity
      INTEGER, OPTIONAL, INTENT(IN)   :: UNIT
      INTEGER :: unit_, spin


      IF(.NOT.ASSOCIATED(impurity%U%mat)) STOP "ERROR IN write_impurity: INPUT &
           &ISNT ALLOCATED!"

      CALL dump_message(UNIT = UNIT, TEXT = "################")
      CALL dump_message(UNIT = UNIT, TEXT = "### IMPURITY ###")
      CALL dump_message(UNIT = UNIT, TEXT = "################")

      unit_ = log_unit
      IF(PRESENT(UNIT)) unit_ = UNIT

      WRITE(unit_, '(a, I0)') "# Nb of sites in the impurity : Nc = ", &
           impurity%Nc

      DO spin = 1, SIZE(impurity%Ec)
         write(unit_, *) ' =================================== '
         CALL write_masked_matrix(impurity%Ec(spin), UNIT = UNIT, SHORT = .true.)
      ENDDO

      write(unit_, *) ' =================================== '
      CALL write_masked_real_matrix(impurity%U, UNIT = UNIT, SHORT = .true.)

   end subroutine

   subroutine Nambu_Ec(EcNambu, Ec)

      ! EMBED IMPURITY QUADRATIC ENERGY IN NAMBU MATRICES

      use globalvar_ed_solver, only: energy_global_shift
      use masked_matrix_class, only: masked_matrix_type, new_masked_matrix
      use matrix,              only: diag

      implicit none

      TYPE(masked_matrix_type), INTENT(INOUT) :: EcNambu
      TYPE(masked_matrix_type), INTENT(IN)    :: Ec(:)
      INTEGER :: Nc ! for clarity only

      IF(SIZE(Ec) == 0) STOP "ERROR IN Nambu_Ec: INPUT Ec ISNT ALLOCATED!"

      Nc = Ec(1)%rc%n1
      CALL new_masked_matrix(EcNambu, "EcNambu", Nc*2, Nc*2, IS_HERM = .true.)

      ! UPPER LEFT BLOCK (SPIN UP)
      EcNambu%rc%mat(   1:Nc,     1:Nc)   =             Ec(1)%rc%mat

      ! LOWER RIGHT BLOCK (SPIN DOWN)
      EcNambu%rc%mat(Nc + 1:Nc*2, Nc + 1:Nc*2) = - TRANSPOSE(Ec(2)%rc%mat)

      energy_global_shift =  sum(diag( Ec(2)%rc%mat ))

   end subroutine

   function average_chem_pot(impurity)

      use matrix, only: new_diag

      implicit none

      TYPE(impurity_type), INTENT(IN) :: impurity
      REAL(DBL) :: average_chem_pot
      INTEGER   :: spin, nspin
      LOGICAL   :: is_diag(impurity%Nc, impurity%Nc)

      CALL new_diag(is_diag, impurity%Nc)

      ! COMPUTE AVERAGE CHEMICAL POTENTIAL
      nspin = SIZE(impurity%Ec)

      DO spin = 1, nspin
         average_chem_pot = average_chem_pot + SUM(impurity%Ec(spin)%rc%mat, &
              is_diag)
      ENDDO
      average_chem_pot = average_chem_pot / ( nspin * impurity%Nc )

   end function

   subroutine shift_average_chem_pot(new_chem_pot, impurity)

      use matrix, only: new_diag

      implicit none

      REAL(DBL), INTENT(IN) :: new_chem_pot
      TYPE(impurity_type) :: impurity
      REAL(DBL)           :: mean_chem_pot
      INTEGER             :: spin, nspin
      LOGICAL             :: is_diag(impurity%Nc, impurity%Nc)

      CALL new_diag(is_diag, impurity%Nc)

      ! SHIFT AVERAGE CHEMICAL POTENTIAL TO new_chem_pot
      nspin = SIZE(impurity%Ec)

      ! COMPUTE OLD AVERAGE CHEMICAL POTENTIAL
      mean_chem_pot = average_chem_pot(impurity)

      DO spin = 1, nspin
         WHERE(is_diag)
         impurity%Ec(spin)%rc%mat = impurity%Ec(spin)%rc%mat + new_chem_pot - &
              mean_chem_pot
         END WHERE
      ENDDO

   end subroutine

end module
