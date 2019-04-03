MODULE apply_C

   use genvar, only: log_unit

   implicit none

   private

   !--------------------------------------------------------------------!
   ! CLASS TO APPLY BASIS CREATION/DESTRUCTION OPERATORS ON EIGENSTATE  !
   !--------------------------------------------------------------------!

   INTEGER, POINTER   :: IMPiorb(:, :) => NULL(), &
                         AIMIMPiorbupdo(:, :) => NULL()
   INTEGER, POINTER   :: AIMIMPiorbsz_updo(:, :) => NULL(), &
                         AIMIMPiorbsz(:) => NULL()
   INTEGER            ::  Nc

   logical, parameter :: force_reset_list = .false.

   public :: apply_cdo
   public :: apply_cup
   public :: apply_n_cdo
   public :: apply_n_cup
   public :: cdo_sector
   public :: cup_sector
   public :: init_apply_c

contains

   subroutine init_apply_C(AIM)

      use aim_class, only: aim_type

      implicit none

      TYPE(AIM_type), INTENT(IN) :: AIM

      Nc = AIM%Nc
      IMPiorb => AIM%impurity%iorb ! not used
      IF(ASSOCIATED(AIMIMPiorbupdo)) DEALLOCATE(AIMIMPiorbupdo)
      ALLOCATE(AIMIMPiorbupdo(AIM%Nc, 2))
      AIMIMPiorbupdo(:, 1) = AIM%IMPiorb(:, 1)
      AIMIMPiorbupdo(:, 2) = AIM%IMPiorb(:, 1) ! rescale ( we never work in
      ! full (up, do) basis!)
      IF(ASSOCIATED(AIMIMPiorbsz_updo)) DEALLOCATE(AIMIMPiorbsz_updo)
      ALLOCATE(AIMIMPiorbsz_updo(AIM%Nc, 2))
      AIMIMPiorbsz_updo(:, 1) = AIM%IMPiorb(:, 1)
      AIMIMPiorbsz_updo(:, 2) = AIM%IMPiorb(:, 2)
      IF(ASSOCIATED(AIMIMPiorbsz)) DEALLOCATE(AIMIMPiorbsz)
      ALLOCATE(AIMIMPiorbsz(AIM%Nc*2))
      AIMIMPiorbsz = (/AIM%IMPiorb(:, 1), AIM%IMPiorb(:, 2)/)
   end subroutine

   subroutine Csz_sector(Csec, pm, sector_in)

      use sector_class, only: delete_sector, norbs__, npart_func, sector_type
      use fermion_hilbert_class, only: new_fermion_sector

      implicit none

      TYPE(sector_type), INTENT(INOUT) :: Csec
      TYPE(sector_type), INTENT(IN)    :: sector_in
      CHARACTER(LEN = 1), INTENT(IN)   :: pm

      CALL delete_sector(Csec)
      ALLOCATE(Csec%sz)
      IF(pm == '+') CALL new_fermion_sector(Csec%sz, npart_func(sector_in) + &
           1, norbs__(sector_in), SZ = .true.)
      IF(pm == '-') CALL new_fermion_sector(Csec%sz, npart_func(sector_in)-1, &
           norbs__(sector_in), SZ = .true.)
   end subroutine

   subroutine Cup_sector(Csec, pm, sector_in)

      use sector_class, only: sector_type

      implicit none

      TYPE(sector_type), INTENT(INOUT) :: Csec
      TYPE(sector_type), INTENT(IN)    :: sector_in
      CHARACTER(LEN = 1), INTENT(IN)   :: pm

#ifdef DEBUG
      write(log_unit, '(a)') "DEBUG: entering apply_C_Cup_sector"
#endif

      IF     (ASSOCIATED(sector_in%sz))  THEN
         CALL Csz_sector_updo(Csec, pm, 1, sector_in)
      ELSE IF(ASSOCIATED(sector_in%updo))THEN
         CALL    Cupdo_sector(Csec, pm, 1, sector_in)
      ENDIF


#ifdef DEBUG
      write(log_unit, '(a)') "DEBUG: leaving apply_C_Cup_sector"
#endif
   end subroutine

   subroutine Cdo_sector(Csec, pm, sector_in)

      use sector_class, only: sector_type

      implicit none

      TYPE(sector_type), INTENT(INOUT) :: Csec
      TYPE(sector_type), INTENT(IN)    :: sector_in
      CHARACTER(LEN = 1), INTENT(IN)   :: pm

#ifdef DEBUG
      write(log_unit, '(a)') "DEBUG: entering apply_C_Cdo_sector"
#endif

      IF     (ASSOCIATED(sector_in%sz))  THEN
         CALL Csz_sector_updo(Csec, pm, 2, sector_in)
      ELSE IF(ASSOCIATED(sector_in%updo))THEN
         CALL    Cupdo_sector(Csec, pm, 2, sector_in)
      ENDIF

#ifdef DEBUG
      write(log_unit, '(a)') "DEBUG: leaving apply_C_Cdo_sector"
#endif

   end subroutine

   subroutine apply_Cup(Ces, pm, MASK, es, esrank)

      use eigen_sector_class, only: eigensector_type

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: Ces
      CHARACTER(LEN = 1), INTENT(IN)        :: pm
      TYPE(eigensector_type), INTENT(IN)    :: es
      LOGICAL, INTENT(IN)                   :: MASK(:)
      INTEGER, INTENT(IN)                   :: esrank

#ifdef DEBUG
      write(log_unit, '(a)') "DEBUG: entering apply_C_apply_Cup"
#endif

      IF     (ASSOCIATED(es%sector%sz))  THEN
         CALL apply_Csz_updo(Ces, pm, 1, MASK, es, esrank)
      ELSE IF(ASSOCIATED(es%sector%updo))THEN
         CALL    apply_Cupdo(Ces, pm, 1, MASK, es, esrank)
      ENDIF

#ifdef DEBUG
      write(log_unit, '(a)') "DEBUG: leaving apply_C_apply_Cup"
#endif

   end subroutine

   subroutine apply_N_Cup(Ces, pm, MASK, es, esrank)

      use eigen_sector_class, only: eigensector_type

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: Ces
      CHARACTER(LEN = 1), INTENT(IN)        :: pm
      TYPE(eigensector_type), INTENT(IN)    :: es
      LOGICAL, INTENT(IN)                   :: MASK(:)
      INTEGER, INTENT(IN)                   :: esrank

      IF     (ASSOCIATED(es%sector%sz))  THEN
         CALL apply_Csz_N_updo(Ces, pm, 1, MASK, es, esrank)
      ELSE IF(ASSOCIATED(es%sector%updo))THEN
         write(*, *) 'N_CUP not working with up dn tensor product basis'
         stop
      ENDIF
   end subroutine

   subroutine apply_Cdo(Ces, pm, MASK, es, esrank)

      use eigen_sector_class, only: eigensector_type

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: Ces
      CHARACTER(LEN = 1), INTENT(IN)        :: pm
      LOGICAL, INTENT(IN)                   :: MASK(:)
      TYPE(eigensector_type), INTENT(IN)    :: es
      INTEGER, INTENT(IN)                   :: esrank

#ifdef DEBUG
      write(log_unit, '(a)') "DEBUG: entering apply_C_apply_Cdo"
#endif

      IF     (ASSOCIATED(es%sector%sz))  THEN
         CALL apply_Csz_updo(Ces, pm, 2, MASK, es, esrank)
      ELSE IF(ASSOCIATED(es%sector%updo))THEN
         CALL    apply_Cupdo(Ces, pm, 2, MASK, es, esrank)
      ENDIF

#ifdef DEBUG
      write(log_unit, '(a)') "DEBUG: leaving apply_C_apply_Cdo"
#endif

   end subroutine

   subroutine apply_N_Cdo(Ces, pm, MASK, es, esrank)

      use eigen_sector_class, only: eigensector_type

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: Ces
      CHARACTER(LEN = 1), INTENT(IN)        :: pm
      LOGICAL, INTENT(IN)                   :: MASK(:)
      TYPE(eigensector_type), INTENT(IN)    :: es
      INTEGER, INTENT(IN)                   :: esrank

      IF     (ASSOCIATED(es%sector%sz))  THEN
         CALL apply_Csz_N_updo(Ces, pm, 2, MASK, es, esrank)
      ELSE IF(ASSOCIATED(es%sector%updo))THEN
         write(*, *) 'cor hopping N Cup not working in up dn tensor product &
              &basis'
         stop
      ENDIF
   end subroutine

   subroutine Csz_sector_updo(Csec, pm, spin, sector_in)

      use sector_class, only: delete_sector, norbs__, npart_func, sector_type
      use fermion_hilbert_class, only: new_fermion_sector

      implicit none

      TYPE(sector_type), INTENT(INOUT) :: Csec
      TYPE(sector_type), INTENT(IN)    :: sector_in
      CHARACTER(LEN = 1), INTENT(IN)   :: pm
      INTEGER, INTENT(IN)              :: spin

      CALL delete_sector(Csec)

      ALLOCATE(Csec%sz)

      SELECT CASE(spin)
      CASE(1)
         IF(pm == '+') CALL new_fermion_sector(Csec%sz, &
              npart_func(sector_in) + 1, norbs__(sector_in), SZ = .true.)
         IF(pm == '-') CALL new_fermion_sector(Csec%sz, &
              npart_func(sector_in)-1, norbs__(sector_in), SZ = .true.)
      CASE(2) ! NAMBU
         IF(pm == '+') CALL new_fermion_sector(Csec%sz, &
              npart_func(sector_in)-1, norbs__(sector_in), SZ = .true.)
         IF(pm == '-') CALL new_fermion_sector(Csec%sz, npart_func(sector_in) &
              + 1, norbs__(sector_in), SZ = .true.)
      END SELECT

   end subroutine

   subroutine Cupdo_sector(Csec, pm, spin, sector_in)

      use sector_class, only: delete_sector, norbs__, npart_func, sector_type
      use fermion_sector2_class, only: new_fermion_sector2

      implicit none

      TYPE(sector_type), INTENT(INOUT) :: Csec
      TYPE(sector_type), INTENT(IN)    :: sector_in
      CHARACTER(LEN = 1), INTENT(IN)   :: pm
      INTEGER, INTENT(IN)              :: spin
      INTEGER :: nup, ndo, Ns, ii

      CALL delete_sector(Csec)

      Ns  = norbs__(sector_in) / 2
      nup = npart_func(sector_in, 1)
      ndo = npart_func(sector_in, 2)

      ALLOCATE(Csec%updo)

      SELECT CASE(spin)
      CASE(1) ! SPIN UP
         IF(pm == '+') then
            CALL new_fermion_sector2(Csec%updo, nup + 1, ndo, Ns)
         endif
         IF(pm == '-') then
            CALL new_fermion_sector2(Csec%updo, nup-1, ndo, Ns)
         endif
      CASE(2) ! SPIN DOWN
         IF(pm == '+') then
            CALL new_fermion_sector2(Csec%updo, nup, ndo + 1, Ns)
         endif
         IF(pm == '-') then
            CALL new_fermion_sector2(Csec%updo, nup, ndo-1, Ns)
         endif

      END SELECT

   end subroutine

   subroutine apply_Csz(Ces, pm, MASK, es, esrank)

      ! $$ COMPUTE c[site]|0>, c^+[site]|0>

      use eigen_class, only: add_eigen, delete_eigen, delete_eigenlist, &
           eigen_type, new_eigen
      use eigen_sector_class, only: eigensector_type
      use fermion_ket_class,  only: create, destroy, fermion_ket_type, new_ket
      use genvar,             only: DBL
      use sector_class,       only: dimen_func

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: Ces
      CHARACTER(LEN = 1), INTENT(IN)        :: pm
      LOGICAL, INTENT(IN)                   :: MASK(:)
      TYPE(eigensector_type), INTENT(IN)    :: es
      INTEGER, INTENT(IN)                   :: esrank
      TYPE(eigen_type)          :: eigen_out
      TYPE(eigen_type), POINTER :: eigen_in => NULL()
      TYPE(fermion_ket_type)    :: ket_in, ket_out
      INTEGER                   :: istate, jstate, iorb

      eigen_in => es%lowest%eigen(esrank)
      if(force_reset_list) CALL delete_eigenlist(Ces%lowest)
      CALL new_eigen(eigen_out, dimen_func(Ces%sector))

      eigen_out%val       = eigen_in%val
      eigen_out%converged = eigen_in%converged

      ! WE ONLY NEED TO CONSIDER A SUBSET OF ORBITALS
      DO iorb = 1, SIZE(MASK)
         IF(MASK(iorb))THEN

            ! FIRST WE CREATE THE OUTPUT VECTOR IN THE RELEVANT SECTOR

            eigen_out%rank   = iorb
            eigen_out%vec%rc = 0.0_DBL

            ! THEN PARSE THE INPUT SECTOR TO APPLY RELEVANT
            ! CREATION/ANNIHILATION OPERATOR

            DO istate = 1, es%sector%sz%dimen
               IF(eigen_in%vec%rc(istate) /= 0.0_DBL)THEN
                  CALL new_ket(ket_in, es%sector%sz%state(istate), &
                       es%sector%sz%norbs)
                  IF(pm == '+') CALL create(ket_out, AIMIMPiorbsz(iorb), &
                       ket_in) ! |Ceigen >= c^ + |eigen >
                  IF(pm == '-') CALL destroy(ket_out, AIMIMPiorbsz(iorb), &
                       ket_in) ! |Ceigen >= c |eigen >

                  IF(.NOT.ket_out%is_nil)THEN
                     jstate = Ces%sector%sz%rank(ket_out%state)
                     eigen_out%vec%rc(jstate) = eigen_out%vec%rc(jstate) + &
                          ket_out%fermion_sign * eigen_in%vec%rc(istate)
                  ENDIF

               ENDIF
            ENDDO

            ! THEN WE UPDATE THE LIST OF OUTPUT VECTORS

            CALL add_eigen(eigen_out, Ces%lowest)
         ENDIF
      ENDDO
      CALL delete_eigen(eigen_out)

   end subroutine

   subroutine apply_Csz_updo(Ces, pm, spin, MASK, es, esrank)

      ! $$ COMPUTE c[site]|0 >, c^ + [site]|0 >

      use eigen_class, only: add_eigen, delete_eigen, delete_eigenlist, &
           eigen_type, new_eigen
      use fermion_ket_class,  only: create, destroy, fermion_ket_type, new_ket
      use eigen_sector_class, only: eigensector_type
      use genvar,             only: DBL
      use sector_class,       only: dimen_func

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: Ces
      CHARACTER(LEN = 1), INTENT(IN)        :: pm
      INTEGER, INTENT(IN)                   :: spin
      LOGICAL, INTENT(IN)                   :: MASK(:)
      TYPE(eigensector_type), INTENT(IN)    :: es
      INTEGER, INTENT(IN)                   :: esrank
      TYPE(eigen_type)          :: eigen_out
      TYPE(eigen_type), POINTER :: eigen_in => NULL()
      TYPE(fermion_ket_type)    :: ket_in, ket_out
      INTEGER                   :: istate, jstate, site

      eigen_in => es%lowest%eigen(esrank)
      if(force_reset_list) CALL delete_eigenlist(Ces%lowest)
      CALL new_eigen(eigen_out, dimen_func(Ces%sector))

      eigen_out%val        = eigen_in%val
      eigen_out%converged  = eigen_in%converged

      ! WE ONLY NEED TO CONSIDER A SUBSET OF ORBITALS
      DO site = 1, SIZE(MASK)
         IF(MASK(site))THEN

            ! FIRST WE CREATE THE OUTPUT VECTOR IN THE RELEVANT SECTOR

            eigen_out%rank   = site
            eigen_out%vec%rc = 0.0_DBL

            ! THEN PARSE THE INPUT SECTOR TO APPLY RELEVANT
            ! CREATION/ANNIHILATION OPERATOR

            DO istate = 1, es%sector%sz%dimen
               IF(eigen_in%vec%rc(istate) /= 0.0_DBL)THEN
                  CALL new_ket(ket_in, es%sector%sz%state(istate), &
                       es%sector%sz%norbs)
                  SELECT CASE(spin)
                  CASE(1)
                     IF(pm == '+') CALL create(ket_out, &
                          AIMIMPiorbsz_updo(site, spin), ket_in) ! |Ceigen >=
                     ! c^ + |eigen >
                     IF(pm == '-') CALL destroy(ket_out, &
                          AIMIMPiorbsz_updo(site, spin), ket_in) ! |Ceigen >=
                     ! c |eigen >
                  CASE(2) ! NAMBU
                     IF(pm == '+') CALL destroy(ket_out, &
                          AIMIMPiorbsz_updo(site, spin), ket_in) ! |Ceigen >=
                     ! c |eigen >
                     IF(pm == '-') CALL create(ket_out, &
                          AIMIMPiorbsz_updo(site, spin), ket_in) ! |Ceigen >=
                     ! c^ + |eigen >
                  END SELECT
                  IF(.NOT.ket_out%is_nil)THEN
                     jstate = Ces%sector%sz%rank(ket_out%state)
                     eigen_out%vec%rc(jstate) = eigen_out%vec%rc(jstate) + &
                          ket_out%fermion_sign * eigen_in%vec%rc(istate)
                  ENDIF
               ENDIF
            ENDDO

            ! THEN WE UPDATE THE LIST OF OUTPUT VECTORS

            CALL add_eigen(eigen_out, Ces%lowest)
         ENDIF
      ENDDO
      CALL delete_eigen(eigen_out)

   end subroutine

   subroutine apply_Csz_N_updo(Ces, pm, spin, MASK, es, esrank)

      ! $$ COMPUTE c[site]|0 >, c^ + [site]|0 >

      use eigen_class,        only: add_eigen, delete_eigen, delete_eigenlist, &
           eigen_type, new_eigen
      use fermion_ket_class,  only: create, destroy, fermion_ket_type, new_ket, &
           ni__
      use eigen_sector_class, only: eigensector_type
      use genvar,             only: DBL
      use sector_class,       only: dimen_func

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: Ces
      CHARACTER(LEN = 1), INTENT(IN)        :: pm
      INTEGER, INTENT(IN)                   :: spin
      LOGICAL, INTENT(IN)                   :: MASK(:)
      TYPE(eigensector_type), INTENT(IN)    :: es
      INTEGER, INTENT(IN)                   :: esrank
      TYPE(eigen_type)          :: eigen_out
      TYPE(eigen_type), POINTER :: eigen_in => NULL()
      TYPE(fermion_ket_type)    :: ket_in, ket_out
      INTEGER                   :: istate, jstate, site, goforit, NUP, NDN, i, &
                                   nn
      INTEGER                   :: tt1, tt2



      eigen_in => es%lowest%eigen(esrank)
      if(force_reset_list) CALL delete_eigenlist(Ces%lowest)
      CALL new_eigen(eigen_out, dimen_func(Ces%sector))

      eigen_out%val        = eigen_in%val
      eigen_out%converged  = eigen_in%converged

      ! WE ONLY NEED TO CONSIDER A SUBSET OF ORBITALS

      DO site = 1, SIZE(MASK)
         IF(MASK(site))THEN

            ! FIRST WE CREATE THE OUTPUT VECTOR IN THE RELEVANT SECTOR

            eigen_out%rank   = site
            eigen_out%vec%rc = 0.0_DBL

            ! THEN PARSE THE INPUT SECTOR TO APPLY RELEVANT
            ! CREATION/ANNIHILATION OPERATOR

            DO istate = 1, es%sector%sz%dimen
               IF(eigen_in%vec%rc(istate) /= 0.0_DBL)THEN
                  CALL new_ket(ket_in, es%sector%sz%state(istate), &
                       es%sector%sz%norbs)

                  nn = es%sector%sz%norbs
                  tt1 = AIMIMPiorbsz_updo(site, 1)
                  tt2 = AIMIMPiorbsz_updo(site, 2)
                  NUP = ni__(tt1, ket_in)
                  NDN = ni__(tt2, ket_in)
                  goforit = 0
                  if(spin == 1)then !UP
                     if(NDN == 0) goforit = 1 !Nambu
                  else            !DN
                     if(NUP == 1) goforit = 1
                  endif
                  if(goforit == 1)then
                     SELECT CASE(spin)
                     CASE(1)
                        IF(pm == '+') CALL create(ket_out, &
                             AIMIMPiorbsz_updo(site, spin), ket_in) ! |Ceigen
                        ! >= c^ + |eigen >
                        IF(pm == '-') CALL destroy(ket_out, &
                             AIMIMPiorbsz_updo(site, spin), ket_in) ! |Ceigen
                        ! >= c |eigen >
                     CASE(2) ! NAMBU
                        IF(pm == '+') CALL destroy(ket_out, &
                             AIMIMPiorbsz_updo(site, spin), ket_in) ! |Ceigen
                        ! >= c |eigen >
                        IF(pm == '-') CALL create(ket_out, &
                             AIMIMPiorbsz_updo(site, spin), ket_in) ! |Ceigen
                        ! >= c^ + |eigen >
                     END SELECT
                  else
                     ket_out%is_nil = .true.
                  endif

                  IF(.NOT.ket_out%is_nil)THEN
                     jstate = Ces%sector%sz%rank(ket_out%state)
                     eigen_out%vec%rc(jstate) = eigen_out%vec%rc(jstate) + &
                          ket_out%fermion_sign * eigen_in%vec%rc(istate)
                  ENDIF

               ENDIF
            ENDDO

            ! THEN WE UPDATE THE LIST OF OUTPUT VECTORS

            CALL add_eigen(eigen_out, Ces%lowest)
         ENDIF
      ENDDO
      CALL delete_eigen(eigen_out)

   end subroutine

   subroutine apply_Cupdo(Ces, pm, spin, MASK, es, esrank)

      ! $$ COMPUTE c[site, spin]|0 >, c^ + [site, spin]|0 >

      use eigen_class,           only: add_eigen, delete_eigen, delete_eigenlist, &
           eigen_type, new_eigen
      use eigen_sector_class,    only: eigensector_type
      use fermion_hilbert_class, only: fermion_sector_type
      use fermion_ket_class,     only: create, destroy, fermion_ket_type, new_ket
      use fermion_sector2_class, only: tabrankupdo
      use genvar,                only: DBL, iproc, size2
      use globalvar_ed_solver,   only: use_transpose_trick_mpi
      use linalg,                only: ramp
      use mpirout,               only: mpibarrier, scatter_it
      use sector_class,          only: dimen_func

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: Ces
      CHARACTER(LEN = 1), INTENT(IN)        :: pm
      INTEGER, INTENT(IN)                   :: spin
      LOGICAL, INTENT(IN)                   :: MASK(:)
      TYPE(eigensector_type), INTENT(IN)    :: es
      INTEGER, INTENT(IN)                   :: esrank
      TYPE(fermion_sector_type), POINTER :: Csec => NULL(), sec => NULL()
      TYPE(eigen_type)                   :: eigen_out
      TYPE(eigen_type), POINTER          :: eigen_in => NULL()
      TYPE(fermion_ket_type)             :: ket_in, ket_out
      INTEGER                            :: istate, jstate, site, dim_stride, &
                                            iup, ido
      INTEGER, ALLOCATABLE               :: myramp(:), tabistate(:), &
                                            tabjstate(:)
      INTEGER                            :: ii
#ifdef _complex
      COMPLEX(8), ALLOCATABLE            :: tempp(:)
#else
      REAL(8), ALLOCATABLE               :: tempp(:)
#endif

#ifdef DEBUG
      write(log_unit, '(a)') "DEBUG: entering apply_C_apply_Cupdo"
#endif

      eigen_in => es%lowest%eigen(esrank)
      if(force_reset_list) CALL delete_eigenlist(Ces%lowest)

      SELECT CASE(spin)
      CASE(1) ! SPIN UP
         sec  =>  es%sector%updo%up
         Csec => Ces%sector%updo%up
         dim_stride = es%sector%updo%down%dimen
      CASE(2) ! SPIN DOWN
         sec  =>  es%sector%updo%down
         Csec => Ces%sector%updo%down
         dim_stride = es%sector%updo%up%dimen
      END SELECT

      if(USE_TRANSPOSE_TRICK_MPI)then
         CALL new_eigen(eigen_out, &
              Ces%sector%updo%up%dimen*Ces%sector%updo%down%dimen/size2)
      else
         CALL new_eigen(eigen_out, dimen_func(Ces%sector))
      endif

      eigen_out%val         =  eigen_in%val
      eigen_out%converged   =  eigen_in%converged
      eigen_out%dim_space = &
           Ces%sector%updo%up%dimen*Ces%sector%updo%down%dimen
      eigen_out%dim_sector  =  eigen_out%dim_space/size2

      call collect_on_rank0(es%lowest%eigen(esrank))

      ALLOCATE(myramp(dim_stride))
      CALL ramp(myramp, dim_stride)
      ALLOCATE(tabistate(dim_stride), tabjstate(dim_stride))

      ! WE ONLY NEED TO CONSIDER A SUBSET OF ORBITALS
      DO site = 1, SIZE(MASK)
         IF(MASK(site))THEN

            ! FIRST WE CREATE THE OUTPUT VECTOR IN THE RELEVANT SECTOR

            eigen_out%rank   = site
            eigen_out%vec%rc = 0.0_DBL

            ! THEN PARSE THE RELEVANT SPIN SECTOR TO APPLY RELEVANT
            ! CREATION/ANNIHILATION OPERATOR

            call collect_on_rank0(eigen_out)
            if(iproc /= 1 .and. USE_TRANSPOSE_TRICK_MPI) goto 35

            DO istate = 1, sec%dimen
               CALL new_ket(ket_in, sec%state(istate), sec%norbs)
               IF(pm == '+') CALL create(ket_out, AIMIMPiorbupdo(site, &
                    spin), ket_in) ! |Ceigen >= c^ + |eigen >
               IF(pm == '-') CALL destroy(ket_out, AIMIMPiorbupdo(site, spin), &
                    ket_in) ! |Ceigen >= c |eigen >
               IF(.NOT.ket_out%is_nil)THEN
                  jstate = Csec%rank(ket_out%state)
                  SELECT CASE(spin)
                  CASE(1) ! SPIN UP
                     CALL tabrankupdo(tabistate, istate, myramp, es%sector%updo)
                     CALL tabrankupdo(tabjstate, jstate, myramp, &
                          Ces%sector%updo)
                  CASE(2) ! SPIN DOWN
                     CALL tabrankupdo(tabistate, myramp, istate, es%sector%updo)
                     CALL tabrankupdo(tabjstate, myramp, jstate, &
                          Ces%sector%updo)
                  END SELECT
                  eigen_out%vec%rc(tabjstate) = eigen_out%vec%rc(tabjstate) + &
                       ket_out%fermion_sign * eigen_in%vec%rc(tabistate)
               ENDIF
            ENDDO

            35   continue
            call mpibarrier
            call scatter_rank0(eigen_out)

            CALL add_eigen(eigen_out, Ces%lowest)

         ENDIF
      ENDDO

      CALL delete_eigen(eigen_out)

      IF(ALLOCATED(myramp))    DEALLOCATE(myramp)
      IF(ALLOCATED(tabistate)) DEALLOCATE(tabistate)
      IF(ALLOCATED(tabjstate)) DEALLOCATE(tabjstate)

      call scatter_rank0(es%lowest%eigen(esrank))

#ifdef DEBUG
      write(log_unit, '(a)') "DEBUG: leaving apply_C_apply_Cupdo"
#endif

   contains

      subroutine collect_on_rank0(eigen)

         use mpirout,             only: mpibarrier, mpigather_on_masternode
         use globalvar_ed_solver, only: use_transpose_trick_mpi
         use genvar,              only: iproc
         use eigen_class,         only: eigen_type

         implicit none

         TYPE(eigen_type) :: eigen

         if(USE_TRANSPOSE_TRICK_MPI)then
            call mpibarrier
            if(iproc == 1) then
               allocate(tempp(eigen%dim_space))
            else
               allocate(tempp(1))
            endif
            call mpigather_on_masternode(tempp, eigen%vec%rc)
            if(iproc == 1) then
               deallocate(eigen%vec%rc)
               allocate(eigen%vec%rc(size(tempp)))
               eigen%vec%rc = tempp
            endif
            deallocate(tempp)
         endif
      end subroutine

      subroutine scatter_rank0(eigen)

         use globalvar_ed_solver, only: use_transpose_trick_mpi
         use genvar,              only: iproc
         use mpirout,             only: mpibarrier, scatter_it
         use eigen_class,         only: eigen_type

         implicit none

         TYPE(eigen_type) :: eigen

         if(USE_TRANSPOSE_TRICK_MPI)then
            call mpibarrier
            allocate(tempp(eigen%dim_sector))
            call scatter_it(eigen%vec%rc, tempp)
            if(iproc == 1) then
               deallocate(eigen%vec%rc)
               allocate(eigen%vec%rc(size(tempp)))
            endif
            eigen%vec%rc = tempp
            deallocate(tempp)
         endif
      end subroutine

   end subroutine

end module
