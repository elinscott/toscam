MODULE fermion_Hilbert_class

   ! use specialFunction,      only : factorial_rec, combination
   ! use mpirout,              only : split
   ! use readable_vec_class

   IMPLICIT NONE

   private

   ! BREAK DOWN HILBERT SPACE 2^Ns WITH norbs FERMION ORBITALS INTO ORTHOGONAL
   ! SECTORS CONTAINING 0 <= N <= norbs FERMIONS !

   ! TOTAL HILBERT SPACE OF SPINLESS FERMIONS

   TYPE fermion_Hilbert_type
      INTEGER                            :: norbs     = 0        ! total number of 1-particle orbitals (# sites)
      INTEGER(8)                         :: dimen     = 0        ! dimension of total Hilbert space
      INTEGER                            :: nsectors  = 0        ! total number of sectors
      INTEGER                            :: max_dimen = 2**30+1  ! max dimension for keeping into memory the Hilbert%rank array
      TYPE(fermion_sector_type), POINTER :: sector(:) => NULL()  ! array of sectors (dim=nsectors)
      INTEGER,                   POINTER :: rank(:)   => NULL()  ! 1 <= state <= dimen=2^norbs
   END TYPE

   ! FERMION SECTOR TYPE

   TYPE fermion_sector_type
      LOGICAL            :: SZ        = .false.      !
      LOGICAL            :: not_physical = .false.
      INTEGER            :: norbs     = 0            ! total number of 1-particle fermion orbitals (# sites)
      INTEGER            :: npart     = 0            ! # of occupied orbitals (# spinless fermions)
      INTEGER            :: nstates   = 0            ! # of 1-fermion states (=2**norbs)
      CHARACTER(LEN=100) :: title     = '\0'         ! writable spin label: npart or Sz
      INTEGER            :: dimen     = 0            ! dimension of the sector
      INTEGER, POINTER   :: state(:)     => NULL()   ! dimen states of the sector in integer representation 
      INTEGER, POINTER   :: rank(:)      => NULL()   ! points towards Hilbert%rank !
      INTEGER, POINTER   :: chunk(:)     => NULL()   ! size of the chunk of the sector actually stored on processors
      INTEGER, POINTER   :: istatemin(:) => NULL()   ! bounds of the chunk
      INTEGER, POINTER   :: istatemax(:) => NULL()   !
   END TYPE

   INTERFACE new_fermion_sector
      MODULE PROCEDURE new_fermion_sector_from_scratch
      MODULE PROCEDURE new_fermion_sector_from_old
   END INTERFACE


   LOGICAL                     :: HILBERT_SPACE_SPLITED_AMONG_NODES, &
                                  FLAG_FORBID_SPLITTING = .false.
   INTEGER,            private :: istati
   !BUG September 6th 2012, Hilbert with save attribute, to compile with
   ! gfortran
   TYPE(fermion_Hilbert_type), private, save        :: Hilbert


contains

   subroutine new_fermion_sector_from_scratch(SEC, npart, norbs, SZ)

      use common_def,            only: c2s, i2c
      use quantum_algebra,       only: fermion

      implicit none

      TYPE(fermion_sector_type), INTENT(INOUT) :: SEC
      INTEGER, INTENT(IN)                      :: norbs ! # of orbitals
      INTEGER, INTENT(IN)                      :: npart ! # of spinless fermions in this sector
      LOGICAL, OPTIONAL, INTENT(IN)            :: SZ
      LOGICAL, SAVE :: first_call = .true.

      CALL delete_fermion_sector(SEC)

      ! IF THIS IS THE FIRST CALL EVER WE NEED TO CREATE THE FULL HILBERT
      ! SPACE I.E. *ALL* THE SECTORS (npart, norbs)

      IF(first_call)THEN
         CALL new_fermion_Hilbert_space(Hilbert, norbs, SZ = SZ)
         first_call = .false.
      ENDIF

      IF(norbs /= Hilbert%norbs) then
         write(*, *) "ERROR IN new_fermion_sector: norbs = " // &
              c2s(i2c(norbs)) // " BUT Hilbert%norbs = " // &
              c2s(i2c(Hilbert%norbs))
         stop
      endif

      IF(PRESENT(SZ)) sec%SZ = SZ

      sec%norbs     =  norbs
      sec%npart     =  npart

      IF(npart < 0)      sec%npart = 0
      if(npart > norbs)  sec%npart = norbs
      sec%not_physical = npart < 0 .or. npart > norbs

      sec%nstates   =  Hilbert%sector(sec%npart)%nstates
      sec%title     =  TRIM(Hilbert%sector(sec%npart)%title)
      sec%dimen     =  Hilbert%sector(sec%npart)%dimen

      ! HILBERT SPACE WAS ALREADY CREATED: JUST NEED TO POINT TO HILBERT
      ! ARRAYS

      sec%state => Hilbert%sector(sec%npart)%state
      sec%rank  => Hilbert%rank

      CALL split_fermion_sector(sec)

      if(.not.associated(SEC%istatemin)) stop 'error istatemin not associated &
           &new fermion sector'
      if(.not.associated(SEC%istatemax)) stop 'error istatemax not associated &
           &new fermion sector'

   end subroutine

   subroutine new_fermion_sector_from_old(SECOUT, SECIN)

      implicit none

      TYPE(fermion_sector_type), INTENT(INOUT) :: SECOUT
      TYPE(fermion_sector_type), INTENT(IN)    :: SECIN

      CALL delete_fermion_sector(SECOUT)
      CALL new_fermion_sector_from_scratch(SECOUT, SECIN%npart, SECIN%norbs, &
           SZ = SECIN%SZ)
      CALL copy_fermion_sector(SECOUT, SECIN)
   end subroutine

   subroutine copy_fermion_sector(SECOUT, SECIN)

      implicit none

      TYPE(fermion_sector_type), INTENT(INOUT) :: SECOUT
      TYPE(fermion_sector_type), INTENT(IN)    :: SECIN

      SECOUT%SZ      =  SECIN%SZ
      SECOUT%norbs   =  SECIN%norbs
      SECOUT%npart   =  SECIN%npart
      SECOUT%nstates =  SECIN%nstates
      SECOUT%title   =  TRIM(SECIN%title)
      SECOUT%dimen   =  SECIN%dimen
      SECOUT%state   => Hilbert%sector(SECOUT%npart)%state
      SECOUT%rank    => Hilbert%rank
      CALL split_fermion_sector(SECOUT)
   end subroutine

   subroutine delete_fermion_sector(SEC)

      implicit none

      TYPE(fermion_sector_type), INTENT(INOUT) :: SEC

      IF(ASSOCIATED(SEC%state))     NULLIFY(SEC%state)
      IF(ASSOCIATED(SEC%rank))      NULLIFY(SEC%rank)
      IF(ASSOCIATED(SEC%chunk))     DEALLOCATE(SEC%chunk, STAT = istati)
      IF(ASSOCIATED(SEC%istatemin)) DEALLOCATE(SEC%istatemin, STAT = istati)
      IF(ASSOCIATED(SEC%istatemax)) DEALLOCATE(SEC%istatemax, STAT = istati)
   end subroutine

   subroutine write_raw_fermion_sector(sec, UNIT)

      implicit none

      TYPE(fermion_sector_type), INTENT(IN) :: sec
      INTEGER, INTENT(IN)                   :: UNIT

      WRITE(UNIT, *) sec%norbs
      WRITE(UNIT, *) sec%npart
   end subroutine

   subroutine read_raw_fermion_sector(SEC, UNIT, SZ)

      implicit none

      TYPE(fermion_sector_type), INTENT(INOUT) :: SEC
      INTEGER, INTENT(IN)                      :: UNIT
      LOGICAL, OPTIONAL, INTENT(IN)            :: SZ
      INTEGER :: norbs, npart

      CALL delete_fermion_sector(SEC)
      READ(UNIT, *) norbs
      READ(UNIT, *) npart
      CALL new_fermion_sector(SEC, npart, norbs, SZ = SZ)
   end subroutine

   subroutine split_fermion_sector(SEC)

      ! DISTRIBUTE THE SECTOR OVER nproc >= 1 PROCESSORS

      use genvar,  only: iproc, log_unit, messages3, no_mpi, nproc, size2
      use mpirout, only: split

      implicit none

      TYPE(fermion_sector_type), INTENT(INOUT) :: SEC
      INTEGER :: jproc

      IF(ASSOCIATED(SEC%chunk))     DEALLOCATE(SEC%chunk)
      IF(ASSOCIATED(SEC%istatemin)) DEALLOCATE(SEC%istatemin)
      IF(ASSOCIATED(SEC%istatemax)) DEALLOCATE(SEC%istatemax)

      ALLOCATE(SEC%chunk(nproc), SEC%istatemin(nproc), SEC%istatemax(nproc))

      ! SPLIT_SECTOR_ON_NODES

      HILBERT_SPACE_SPLITED_AMONG_NODES = size2 > 1 .and. &
           .not.FLAG_FORBID_SPLITTING .and. .not.no_mpi

      if(.not.HILBERT_SPACE_SPLITED_AMONG_NODES)then
         SEC%chunk     = SEC%dimen
         SEC%istatemin = 1
         SEC%istatemax = SEC%dimen
      else
         CALL split(SEC%dimen, SEC%istatemin, SEC%istatemax, SEC%chunk)
         if(messages3)then
            write(log_unit, *) 'the hilbert space was splitted among nodes'
            write(log_unit, *) ' chunk, statemin, statemax on my node is : ', &
                 SEC%chunk(iproc), SEC%istatemin(iproc), SEC%istatemax(iproc)
         endif
      endif

      RETURN
   end subroutine

   function is_in_fermion_sector(state, SEC) RESULT(is_in_fs)

      implicit none

      TYPE(fermion_sector_type), INTENT(IN) :: SEC
      INTEGER, INTENT(IN)                   :: state
      LOGICAL :: is_in_fs
      INTEGER :: iorb, npart

      is_in_fs = .false.
      npart    = COUNT((/(BTEST(state, iorb-1), iorb = 1, SEC%norbs)/))
      IF(npart == SEC%npart) is_in_fs = .true.
   end function

   subroutine allocate_fermion_sector(SEC, npart, norbs, SZ)

      use specialfunction,       only: combination
      use common_def,            only: c2s, i2c

      implicit none

      TYPE(fermion_sector_type), INTENT(INOUT) :: SEC
      INTEGER, INTENT(IN)                      :: norbs ! # of orbitals
      INTEGER, INTENT(IN)                      :: npart ! # of spinless
      ! fermions in this sector
      LOGICAL, OPTIONAL, INTENT(IN)            :: SZ

      CALL delete_fermion_sector(SEC)

      SEC%norbs   = norbs
      SEC%npart   = npart
      SEC%nstates = 2**norbs
      SEC%title   = "Nocc = " // c2s(i2c(SEC%npart))

      IF(PRESENT(SZ))THEN
         IF(SZ) SEC%title = "Sz*2 = " // c2s(i2c(npart-norbs/2))
      ENDIF

      ! FIND DIMENSION OF SECTOR

      SEC%dimen = combination(norbs, npart)

      ! ALLOCATE ARRAY OF ALL STATES (INTEGER REP.) IN THE SECTOR

      ALLOCATE(SEC%state(SEC%dimen))

      ! IT JUST REMAINS TO ALLOCATE RANK (LATER ON) AND DISTRIBUTE THE SECTOR
      ! OVER nproc PROCESSORS

      CALL split_fermion_sector(SEC)

   end subroutine

   subroutine delete_fermion_Hilbert_space(Hilbert)

      implicit none

      TYPE(fermion_Hilbert_type), INTENT(INOUT) :: Hilbert
      INTEGER :: npart

      IF(ASSOCIATED(Hilbert%sector))THEN
         DO npart = 0, Hilbert%norbs
            CALL delete_fermion_sector(Hilbert%sector(npart))
            DEALLOCATE(Hilbert%sector)
         ENDDO
      ENDIF
      IF(ASSOCIATED(Hilbert%rank) .and. Hilbert%dimen < Hilbert%max_dimen) &
           DEALLOCATE(Hilbert%rank)
   end subroutine

   subroutine new_fermion_Hilbert_space(Hilbert, norbs, SZ)

      use fermion_ket_class,     only: fermion_ket_type, new_ket, noccupied
      use common_def,            only: c2s, dump_message, i2c

      implicit none

      TYPE(fermion_Hilbert_type), INTENT(INOUT) :: Hilbert
      INTEGER, INTENT(IN)                       :: norbs
      LOGICAL, OPTIONAL, INTENT(IN)             :: SZ
      INTEGER                :: state, npart, iorb, dimen_check(0:norbs)
      TYPE(fermion_ket_type) :: ket

      CALL delete_fermion_Hilbert_space(Hilbert)
      Hilbert%norbs = norbs

      ! FIND DIMENSION OF THE HILBERT SPACE
      Hilbert%dimen = 2**norbs

      ! FIND NUMBER OF SECTORS
      Hilbert%nsectors = norbs + 1

      ! FIND LIST OF SECTORS AND THEIR DIMENSION
      ALLOCATE(Hilbert%sector(0:norbs))
      DO npart = 0, norbs
         CALL allocate_fermion_sector(Hilbert%sector(npart), npart, norbs, SZ &
              = SZ)
      ENDDO

      if(Hilbert%dimen < Hilbert%max_dimen) then
         ALLOCATE(Hilbert%rank(0:Hilbert%dimen-1))
      else
         write(*, *) 'SIZE of requested hilbert space : ', Hilbert%dimen
         write(*, *) 'MAX size allowed                : ', Hilbert%max_dimen
         STOP 'fermion_hilbert, array too large, you will have to implement a &
              &bissection instead of bookkeeping the rank array'
      endif

      ! FIND THE BASIS STATES IN EVERY SECTOR
      IF(norbs > KIND(state)*8)THEN
         CALL dump_message(TEXT = "ERROR IN make_fermion: norbs = " // &
              c2s(i2c(norbs)) // " > nbits = " // c2s(i2c(KIND(state)*8)))
         CALL dump_message(TEXT = "CHANGE INTEGER KIND TO REPRESENT STATES! &
              &(KIND = " // c2s(i2c(KIND(state))) // ")")
         STOP
      ENDIF

      dimen_check = 0
      DO state = 0, Hilbert%dimen-1
         ! Span the entire subspace with norbs spinless fermion orbitals

         ! find # particles (occupied orbitals)
         CALL new_ket(ket, state, norbs)
         npart = noccupied(ket)

         ! update/check dimension
         dimen_check(npart) = dimen_check(npart) + 1
         IF(dimen_check(npart) > Hilbert%sector(npart)%dimen)THEN
            CALL dump_message(TEXT = "ERROR in new_fermion_Hilbert_space: &
                 &sector npart = " // c2s(i2c(npart)) // " dim_check = " // &
                 c2s(i2c(dimen_check(npart))) // " > dim = " // &
                 c2s(i2c(Hilbert%sector(npart)%dimen)))
            STOP
         ENDIF

         ! update relevant sector
         Hilbert%sector(npart)%state(dimen_check(npart)) = state
         Hilbert%rank(state)                             = dimen_check(npart)

      ENDDO

      ! final check of dimension
      DO npart = 0, norbs
         IF(dimen_check(npart) /= Hilbert%sector(npart)%dimen)THEN
            CALL dump_message(TEXT = "ERROR in new_fermion_Hilbert_space: &
                 &npart = " // c2s(i2c(npart)) // " dimen_check = " // &
                 c2s(i2c(dimen_check(npart))) // " < dim = " // &
                 c2s(i2c(Hilbert%sector(npart)%dimen)))
            STOP
         ENDIF
      ENDDO

      ! NOW ARRAY rank IN EVERY SECTOR POINTS TO Hilbert%rank
      DO npart = 0, norbs
         Hilbert%sector(npart)%rank => Hilbert%rank
      ENDDO

   end subroutine

end module
