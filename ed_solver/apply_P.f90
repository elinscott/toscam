MODULE apply_P

   use genvar, only: DBL

   implicit none

   private

   !-----------------------------------------------------!
   ! CLASS TO APPLY PERMUTATION OPERATORS ON EIGENSTATE  !
   !-----------------------------------------------------!

   integer, allocatable :: CYC3(:, :, :), CYC4(:, :, :)
   integer, allocatable :: IMPiorbupdo(:, :), IMPiorbsz(:, :)
   integer              :: nCYC3, nCYC4, Ns, Ns2
   logical, parameter   :: force_reset_list = .false.

   public :: apply_p3
   public :: apply_p4
   public :: init_apply_p
   public :: p_sector

contains

   subroutine init_apply_P(triplets, quadruplets, AIM)

      ! READ MULTIPLETS

      use aim_class, only: aim_type

      implicit none

      INTEGER, INTENT(IN)        :: triplets(:,:), quadruplets(:,:)
      TYPE(AIM_type), INTENT(IN) :: AIM
      INTEGER :: icyc, site, Nc

      if(ANY(shape(triplets) > 0))then
         nCYC3 = SIZE(triplets, 1)
         !CEDRIC november 20 2012
         if(size(triplets, 2) /= 3)then
            write(*, *) 'shape triplets', shape(triplets)
            return
         endif
         ALLOCATE(CYC3(nCYC3, 3, 2))
         CYC3(:, :, 1) = triplets
         DO icyc = 1, nCYC3
            CYC3(icyc, :, 2) = (/triplets(icyc, 1), (triplets(icyc, &
                 3-(site-1)), site = 1, 3-1)/)
         ENDDO
      endif

      if(ANY(shape(quadruplets) > 0))then
         nCYC4 = SIZE(quadruplets, 1)
         ALLOCATE(CYC4(nCYC4, 4, 2))
         CYC4(:, :, 1) = quadruplets
         DO icyc = 1, nCYC4
            CYC4(icyc, :, 2) = (/quadruplets(icyc, 1), (quadruplets(icyc, &
                 4-(site-1)), site = 1, 4-1)/)
         ENDDO
      endif

      ! ORDERING OF ORBITALS
      Nc  = AIM%Nc
      Ns  = AIM%Ns
      Ns2 = Ns*2
      ALLOCATE(IMPiorbupdo(Nc, 2))
      IMPiorbupdo(:, 1) = AIM%IMPiorb(:, 1)
      IMPiorbupdo(:, 2) = AIM%IMPiorb(:, 1) ! RESCALED: WE NEVER WORK IN THE
      ! FULL (nup, ndo) BASIS
      ALLOCATE(IMPiorbsz(Nc, 2))
      IMPiorbsz  (:, 1) = AIM%IMPiorb(:, 1)
      IMPiorbsz  (:, 2) = AIM%IMPiorb(:, 2)

   end subroutine

   subroutine P_sector(Csec, pm, sector_in)

      use sector_class, only: delete_sector, new_sector, sector_type

      implicit none

      TYPE(sector_type), INTENT(INOUT) :: Csec
      TYPE(sector_type), INTENT(IN)    :: sector_in
      CHARACTER(LEN = 1), INTENT(IN)   :: pm

      CALL delete_sector(Csec)
      CALL new_sector(Csec, sector_in)
   end subroutine

   subroutine apply_P3(Ces, pm, MASK, es, esrank)

      use eigen_sector_class, only: eigensector_type

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: Ces
      CHARACTER(LEN = 1), INTENT(IN)        :: pm
      LOGICAL, INTENT(IN)                   :: MASK(:)
      TYPE(eigensector_type), INTENT(IN)    :: es
      INTEGER, INTENT(IN)                   :: esrank

      CALL apply_Pn(Ces, 3, pm, MASK, es, esrank)
   end subroutine

   subroutine apply_P4(Ces, pm, MASK, es, esrank)

      use eigen_sector_class, only: eigensector_type

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: Ces
      CHARACTER(LEN = 1), INTENT(IN)        :: pm
      LOGICAL, INTENT(IN)                   :: MASK(:)
      TYPE(eigensector_type), INTENT(IN)    :: es
      INTEGER, INTENT(IN)                   :: esrank

      CALL apply_Pn(Ces, 4, pm, MASK, es, esrank)
   end subroutine

   subroutine apply_CHI(Ces, pm, MASK, es, esrank)

      use eigen_class,        only: add_eigen, delete_eigenlist
      use sector_class,       only: delete_sector, sector_type
      use eigen_sector_class, only: delete_eigensector, eigensector_type, &
           new_eigensector

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: Ces
      CHARACTER(LEN = 1), INTENT(IN)        :: pm
      LOGICAL, INTENT(IN)                   :: MASK(:)
      TYPE(eigensector_type), INTENT(IN)    :: es
      INTEGER, INTENT(IN)                   :: esrank
      TYPE(eigensector_type) :: P3p, P3m
      TYPE(sector_type)      :: Psecp, Psecm
      INTEGER                :: ieigen

      if(force_reset_list) CALL delete_eigenlist(Ces%lowest)
      ! APPLY P3 +
      CALL P_sector(Psecp, ' + ', es%sector)
      CALL new_eigensector(P3p, Psecp)
      CALL apply_Pn(P3p, 3, ' + ', MASK, es, esrank)
      ! APPLY P3-
      CALL P_sector(Psecm, '-', es%sector)
      CALL new_eigensector(P3m, Psecm)
      CALL apply_Pn(P3m, 3, '-', MASK, es, esrank)
      ! SUBSTRACT
      DO ieigen = 1, P3p%lowest%neigen
         P3p%lowest%eigen(ieigen)%vec%rc = P3p%lowest%eigen(ieigen)%vec%rc - &
              P3m%lowest%eigen(ieigen)%vec%rc
         CALL add_eigen(P3p%lowest%eigen(ieigen), Ces%lowest)
      ENDDO
      CALL delete_eigensector(P3p)
      CALL delete_eigensector(P3m)
      CALL delete_sector(Psecp)
      CALL delete_sector(Psecm)
   end subroutine

   subroutine apply_Pn(Ces, ncyc, pm, MASK, es, esrank)

      ! $$ COMPUTE P^-[site]|0 >, P^ + [site]|0 >

      use eigen_class,           only: add_eigen, delete_eigen, delete_eigenlist, &
           eigen_type, new_eigen
      use fermion_sector2_class, only: rankupdo
      use fermion_ket_class,     only: fermion_ket_type, new_ket
      use eigen_sector_class,    only: eigensector_type
      use sector_class,          only: dimen_func

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: Ces
      INTEGER, INTENT(IN)                   :: ncyc
      CHARACTER(LEN = 1), INTENT(IN)        :: pm
      LOGICAL, DIMENSION(:), INTENT(IN)     :: MASK
      TYPE(eigensector_type), INTENT(IN)    :: es
      INTEGER, INTENT(IN)                   :: esrank
      TYPE(eigen_type)                    :: eigen_out
      TYPE(eigen_type), POINTER           :: eigen_in => NULL()
      TYPE(fermion_ket_type)              :: ket_in, ket_in_up, ket_in_do
      INTEGER                             :: istate, jstate, iup, ido, jup, jdo
      INTEGER                             :: icyc, jcyc, iP, nP, iprod, nprod
      INTEGER, ALLOCATABLE                :: multiplet(:,:)
      TYPE(fermion_ket_type), ALLOCATABLE :: Pket(:), Pket_up(:), Pket_do(:)

      eigen_in => es%lowest%eigen(esrank)
      if(force_reset_list) CALL delete_eigenlist(Ces%lowest)
      CALL new_eigen(eigen_out, dimen_func(Ces%sector))

      eigen_out%val       = eigen_in%val
      eigen_out%converged = eigen_in%converged

      nprod = 2**ncyc

      SELECT CASE(ncyc)
      CASE(3)
         nP = nCYC3
         ALLOCATE(multiplet(nP, 3))
         IF(pm == ' + ') multiplet = CYC3(:, :, 1)
         IF(pm == '-') multiplet = CYC3(:, :, 2)
      CASE(4)
         nP = nCYC4
         ALLOCATE(multiplet(nP, 4))
         IF(pm == ' + ') multiplet = CYC4(:, :, 1)
         IF(pm == '-') multiplet = CYC4(:, :, 2)
      END SELECT

      ! WE ONLY NEED TO CONSIDER A SUBSET OF PERMUTATIONS

      DO iP = 1, nP
         IF(MASK(iP))THEN

            ! FIRST WE CREATE THE OUTPUT VECTOR IN THE RELEVANT SECTOR

            eigen_out%rank   = iP
            eigen_out%vec%rc = 0.0_DBL

            ! THEN WE SELECT RELEVANT MULTIPLET
            ! THEN PARSE THE INPUT SECTOR TO APPLY RELEVANT PERMUTATION OPERATOR

            IF(ASSOCIATED(es%sector%sz))THEN ! Sz-SECTOR
               ALLOCATE(Pket(nprod))
               DO istate = 1, es%sector%sz%dimen
                  IF(eigen_in%vec%rc(istate) /= 0.0_DBL)THEN
                     CALL new_ket(ket_in, es%sector%sz%state(istate), Ns2)
                     CALL Psz(Pket, multiplet(iP, :), ket_in)
                     DO iprod = 1, nprod
                        IF(.NOT.Pket(iprod)%is_nil)THEN
                           jstate = Ces%sector%sz%rank(Pket(iprod)%state)
                           eigen_out%vec%rc(jstate) = eigen_out%vec%rc(jstate) &
                                + eigen_in%vec%rc(istate) * &
                                Pket(iprod)%fermion_sign
                        ENDIF

                     ENDDO
                  ENDIF
               ENDDO
            ELSE IF(ASSOCIATED(es%sector%updo))THEN ! (nup, ndo) SECTOR
               ALLOCATE(Pket_up(0:nprod-1), Pket_do(0:nprod-1))
               DO ido = 1, es%sector%updo%down%dimen
                  CALL new_ket(ket_in_do, es%sector%updo%down%state(ido), Ns)
                  CALL Pupdo(Pket_do, multiplet(iP, :), 2, ket_in_do)
                  DO iup = 1, es%sector%updo%up%dimen
                     istate = rankupdo(iup, ido, es%sector%updo)
                     IF(eigen_in%vec%rc(istate) /= 0.0_DBL)THEN
                        CALL new_ket(ket_in_up, es%sector%updo%up%state(iup), &
                             Ns)
                        CALL Pupdo(Pket_up, multiplet(iP, :), 1, ket_in_up)
                        ! NOW IT REMAINS TO FORM THE DIRECT PRODUCT
                        DO icyc = 0, nprod-1
                           jcyc = nprod-1-icyc
                           IF((.NOT.Pket_up(jcyc)%is_nil) .AND. &
                                (.NOT.Pket_do(icyc)%is_nil))THEN
                              jup = &
                                   Ces%sector%updo%up%rank(Pket_up(jcyc)%state)
                              jdo = &
                                   Ces%sector%updo%down%rank(Pket_do(icyc)%state)
                              jstate = rankupdo(jup, jdo, es%sector%updo)
                              eigen_out%vec%rc(jstate) = &
                                   eigen_out%vec%rc(jstate) + &
                                   Pket_up(jcyc)%fermion_sign * &
                                   Pket_do(icyc)%fermion_sign * &
                                   eigen_in%vec%rc(istate)
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF

            CALL add_eigen(eigen_out, Ces%lowest)
         ENDIF
      ENDDO

      IF(ALLOCATED(Pket))    DEALLOCATE(Pket)
      IF(ALLOCATED(Pket_up)) DEALLOCATE(Pket_up)
      IF(ALLOCATED(Pket_do)) DEALLOCATE(Pket_do)

      DEALLOCATE(multiplet)
      CALL delete_eigen(eigen_out)

   end subroutine

   subroutine Pupdo(Pket, multi, spin, ket)

      use fermion_ket_class, only: fermion_ket_type, hop, new_ket

      implicit none

      INTEGER, INTENT(IN)                   :: multi(:)
      TYPE(fermion_ket_type), INTENT(INOUT) :: Pket(0:2**SIZE(multi)-1)
      TYPE(fermion_ket_type), INTENT(IN)    :: ket
      INTEGER, INTENT(IN)                   :: spin
      INTEGER :: prod(0:2**SIZE(multi)-1, 0:SIZE(multi))
      INTEGER :: nprod, ioldgen, ioldgenmin, ioldgenmax
      INTEGER :: igen, icyc, ncyc

      CALL new_ket(Pket(0), ket)

      ncyc       = SIZE(multi)
      prod       = 0
      nprod      = 0
      ioldgenmin = 0
      ioldgenmax = 0

      DO igen = 1, ncyc
         DO ioldgen = ioldgenmin, ioldgenmax
            DO icyc = prod(ioldgen, igen-1) + 1, ncyc
               nprod = nprod + 1
               prod(nprod, 0:igen) = (/prod(ioldgen, 0:igen-1), icyc/)
               CALL hop(Pket(nprod), IMPiorbupdo(multi(MOD(icyc, ncyc) + 1), &
                    spin), IMPiorbupdo(multi(icyc), spin), Pket(ioldgen))
            ENDDO
         ENDDO
         ioldgenmin = ioldgenmax + 1
         ioldgenmax = nprod
      ENDDO

      ! now ket(:)%state = (1, X1, X2, X3, X2X1, X3X1, X3X2, X3X2X1)|0 > if
      ! ncyc = 3
      ! now ket(:)%state = (1, X1, X2, X3, X4, X2X1, X3X1, X4X1, X3X2, X4X2,
      ! X4X3, X3X2X1, X4X2X1, X4X3X1, X4X3X2, X4X3X2X1)|0 > if ncyc = 4

   end subroutine

   subroutine Psz(Pket, multi, ket)

      use fermion_ket_class, only: copy_ket, create, destroy, &
           fermion_ket_type, new_ket

      implicit none

      INTEGER, INTENT(IN)                   :: multi(:)
      TYPE(fermion_ket_type), INTENT(INOUT) :: Pket(2**SIZE(multi))
      TYPE(fermion_ket_type), INTENT(IN)    :: ket
      TYPE(fermion_ket_type) :: prodket(0:2**(SIZE(multi) + 1)-2)
      INTEGER                :: iprod, nprod, ioldgen, ioldgenmin, ioldgenmax
      INTEGER                :: igen, icyc, ncyc

      CALL new_ket(prodket(0), ket)
      ncyc  = SIZE(multi)
      nprod = 0
      DO icyc = 1, ncyc
         ioldgenmin = 2**(icyc-1) - 1
         ioldgenmax = 2** icyc    - 2
         ! SPIN UP
         DO ioldgen = ioldgenmin, ioldgenmax
            nprod = nprod + 1
            CALL destroy(prodket(nprod), IMPiorbsz(multi(icyc), 1), &
                 prodket(ioldgen))
            CALL create(prodket(nprod), IMPiorbsz(multi(MOD(icyc, ncyc) + 1), &
                 1), prodket(nprod))
         ENDDO
         ! SPIN DOWN - NAMBU
         DO ioldgen = ioldgenmin, ioldgenmax
            nprod = nprod + 1
            CALL create(prodket(nprod), IMPiorbsz(multi(icyc), 2), &
                 prodket(ioldgen))
            CALL destroy(prodket(nprod), IMPiorbsz(multi(MOD(icyc, ncyc) + 1), &
                 2), prodket(nprod))
         ENDDO
      ENDDO
      DO iprod = 1, 2**ncyc
         CALL copy_ket(Pket(iprod), prodket(2**ncyc-2 + iprod))
      ENDDO
   end subroutine

end module
