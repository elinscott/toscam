MODULE fermion_sector2_class

   use fermion_hilbert_class, only: fermion_sector_type
   use genvar, only: log_unit

   IMPLICIT NONE

   private

   !-----------------------------------------------------!
   ! DIRECT PRODUCT OF 2 FERMION SECTORS, LIKE (nup, ndo) !
   !-----------------------------------------------------!

   TYPE fermion_sector2_type
      LOGICAL                   :: NAMBU = .false.
      LOGICAL                   :: not_physical = .false.
      TYPE(fermion_sector_type) :: up    ! nup sector
      TYPE(fermion_sector_type) :: down  ! ndo sector
      CHARACTER(LEN=100)      :: title = '\0' ! (nup, ndo) in writable form
      INTEGER                   :: norbs = 0    ! total # of orbitals
      INTEGER                   :: npart = 0    ! # of occupied orbitals
      INTEGER                   :: nstates = 0    ! total # of 1-particle states
      INTEGER                   :: dimen = 0    ! dimension of the sector
      INTEGER                   :: strideup = 0, stridedown = 0 ! stride between consecutive up/do states
      INTEGER, POINTER          :: chunk(:) => NULL() ! size of the chunk of the sector actually stored on the local processors
      INTEGER, POINTER          :: istatemin(:) => NULL() ! bounds of the local chunks in the full basis
      INTEGER, POINTER          :: istatemax(:) => NULL()
   END TYPE

   INTERFACE tabrankupdo
      MODULE PROCEDURE rankupdo_strideup
      MODULE PROCEDURE rankupdo_stridedo
   END INTERFACE

   INTERFACE new_fermion_sector2
      MODULE PROCEDURE new_fermion_sector2_from_scratch
      MODULE PROCEDURE new_fermion_sector2_from_old
   END INTERFACE

   public :: copy_fermion_sector2
   public :: delete_fermion_sector2
   public :: fermion_sector2_type
   public :: is_in_fermion_sector2
   public :: new_fermion_sector2
   public :: rank2
   public :: rankupdo
   public :: read_raw_fermion_sector2
   public :: state2
   public :: tabrankupdo
   public :: write_raw_fermion_sector2

contains

   subroutine new_fermion_sector2_from_scratch(SEC, nup, ndo, Ns, NAMBU)

      use genvar, only: messages4
      use common_def, only: c2s, i2c
      use quantum_algebra, only: fermion
      use fermion_hilbert_class, only: new_fermion_sector

      implicit none

      TYPE(fermion_sector2_type), INTENT(INOUT) :: SEC ! this returns a pointer to be passed to 'to_sector'
      INTEGER, INTENT(IN)                       :: Ns ! # of sites
      INTEGER, INTENT(IN)                       :: nup, ndo ! # of particles with spin up, down
      LOGICAL, OPTIONAL, INTENT(IN)             :: NAMBU

      SEC%not_physical = .false.

      IF (nup + ndo < 0 .OR. nup + ndo > Ns*2) then
         write (*, *) 'DANGER nup + ndn < 0 or < Ns'
         SEC%not_physical = .true.
      ENDIF

      CALL delete_fermion_sector2(SEC)

      SEC%title = "(nup, ndo) = ("//c2s(i2c(nup))//", "//c2s(i2c(ndo)) &
                  //")" ! note: this is not Ns-do!

      SEC%NAMBU = .false.
      IF (PRESENT(NAMBU)) SEC%NAMBU = NAMBU

      ! FIRST BUILD THE nup AND ndo SECTORS SEPARATELY

      if (messages4) write (*, *) 'build new up and dn sectors'

      ! spin up
      CALL new_fermion_sector(SEC%up, nup, Ns)
      SEC%up%title = "nup = "//c2s(i2c(nup))

      ! spin down
      CALL new_fermion_sector(SEC%down, ndo, Ns)
      SEC%down%title = "ndo = "//c2s(i2c(ndo))

      if (messages4) write (*, *) 'split among nodes'

      if (.not. associated(SEC%up%istatemin)) stop 'error new fermion sector2 : &
           &sec%up%istatemin not associated'

      ! THEN BUILD THE DIRECT PRODUCT

      SEC%norbs = SEC%up%norbs + SEC%down%norbs
      SEC%npart = SEC%up%npart + SEC%down%npart
      SEC%nstates = SEC%up%nstates*SEC%down%nstates
      SEC%dimen = SEC%up%dimen*SEC%down%dimen

      SEC%not_physical = SEC%not_physical .or. SEC%up%not_physical .or. &
                         SEC%down%not_physical

      ! DISTRIBUTE OVER nproc >= 1 PROCESSORS

      CALL split_fermion_sector2(SEC)

      if (messages4) write (*, *) 'done now return from new fermion 2 from &
           &scratch'

   end subroutine

   subroutine new_fermion_sector2_from_old(SECOUT, SECIN)

      implicit none

      TYPE(fermion_sector2_type), INTENT(INOUT) :: SECOUT
      TYPE(fermion_sector2_type), INTENT(IN)    :: SECIN

      CALL delete_fermion_sector2(SECOUT)
      CALL new_fermion_sector2(SECOUT, SECIN%up%npart, SECIN%down%npart, &
                               SECIN%up%norbs, NAMBU=SECIN%NAMBU)
      CALL copy_fermion_sector2(SECOUT, SECIN)
   end subroutine

   subroutine copy_fermion_sector2(SECOUT, SECIN)

      use fermion_hilbert_class, only: copy_fermion_sector

      implicit none

      TYPE(fermion_sector2_type), INTENT(INOUT) :: SECOUT
      TYPE(fermion_sector2_type), INTENT(IN)    :: SECIN

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: entering fermion_sector2_class_copy_&
           &fermion_sector2"
#endif DEBUG

      SECOUT%NAMBU = SECIN%NAMBU
      SECOUT%title = TRIM(SECIN%title)
      SECOUT%dimen = SECIN%dimen
      SECOUT%norbs = SECIN%norbs
      SECOUT%npart = SECIN%npart
      SECOUT%nstates = SECIN%nstates
      CALL copy_fermion_sector(SECOUT%up, SECIN%up)
      CALL copy_fermion_sector(SECOUT%down, SECIN%down)
      CALL split_fermion_sector2(SECOUT)

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: leaving fermion_sector2_class_copy_&
           &fermion_sector2"
#endif DEBUG

   end subroutine

   subroutine delete_fermion_sector2(sector2)

      use fermion_hilbert_class, only: delete_fermion_sector

      implicit none

      TYPE(fermion_sector2_type), INTENT(INOUT) :: sector2

      CALL delete_fermion_sector(sector2%up)
      CALL delete_fermion_sector(sector2%down)
      IF (ASSOCIATED(sector2%chunk)) DEALLOCATE (sector2%chunk)
      IF (ASSOCIATED(sector2%istatemin)) DEALLOCATE (sector2%istatemin)
      IF (ASSOCIATED(sector2%istatemax)) DEALLOCATE (sector2%istatemax)
   end subroutine

   subroutine split_fermion_sector2(SEC)

      ! DISTRIBUTE THE SECTOR OVER nproc >= 1 PROCESSORS
      ! ORDERING : (iup = [1, dimup], ido = 1), (iup = [1, dimup], ido = 2),
      ! ..., (iup = [1, dimup], ido = dimdo)

      use common_def, only: create_seg_fault
      use genvar, only: messages4, nproc
      use quantum_algebra, only: fermion
      use mpirout, only: split

      implicit none

      TYPE(fermion_sector2_type), INTENT(INOUT) :: SEC
      INTEGER :: jproc

      if (messages4) write (*, *) 'test if sectors are allocated'

      if (.not. associated(SEC%up%istatemin)) then
         write (*, *) 'error split fermion sector 2 SEC%up%istatemin not &
              &associated'
         call create_seg_fault
      endif
      if (.not. associated(SEC%up%istatemax)) then
         write (*, *) 'error split fermion sector 2 SEC%up%istatemax not &
              &associated'
         call create_seg_fault
      endif

      SEC%up%istatemin = 1
      SEC%up%istatemax = SEC%up%dimen
      SEC%up%chunk = SEC%up%dimen

      SEC%strideup = 1             ! STORAGE up   IS CONTIGUOUS
      SEC%stridedown = SEC%up%dimen  ! STORAGE down IS NOT

      IF (ASSOCIATED(SEC%chunk)) DEALLOCATE (SEC%chunk)
      IF (ASSOCIATED(SEC%istatemin)) DEALLOCATE (SEC%istatemin)
      IF (ASSOCIATED(SEC%istatemax)) DEALLOCATE (SEC%istatemax)

      if (messages4) write (*, *) 'allocate sectors'

      ALLOCATE (SEC%chunk(nproc), SEC%istatemin(nproc), SEC%istatemax(nproc))

      if (messages4) write (*, *) 'define istatemin istatemax'

      DO jproc = 1, nproc
         SEC%istatemin(jproc) = rankupdo(SEC%up%istatemin(jproc), &
                                         SEC%down%istatemin(jproc), SEC)
         SEC%istatemax(jproc) = rankupdo(SEC%up%istatemax(jproc), &
                                         SEC%down%istatemax(jproc), SEC)
      ENDDO

      SEC%chunk = SEC%up%chunk*SEC%down%chunk

      if (messages4) write (*, *) 'done, now return from split'

   end subroutine

   function is_in_fermion_sector2(state, sector2) RESULT(is_in_fs2)

      use fermion_hilbert_class, only: is_in_fermion_sector

      implicit none

      TYPE(fermion_sector2_type), INTENT(IN) :: sector2
      INTEGER, INTENT(IN)                    :: state
      LOGICAL :: is_in_fs2
      INTEGER :: stateup, statedo, nstates_sz
      INTEGER :: iorb

      is_in_fs2 = .false.
      nstates_sz = 2**sector2%down%norbs
      stateup = MOD(state, nstates_sz)
      IF (is_in_fermion_sector(stateup, sector2%up)) THEN
         statedo = (state - stateup)/nstates_sz
         IF (is_in_fermion_sector(statedo, sector2%down)) is_in_fs2 = .true.
      ENDIF
   end function

   function rankupdo(iup, ido, sector2) RESULT(ranki)

      ! FINDS THE RANK OF (iup, ido) IN FULL DIRECT PRODUCT BASIS

      implicit none

      INTEGER, INTENT(IN)                    :: iup, ido
      TYPE(fermion_sector2_type), INTENT(IN) :: sector2
      INTEGER :: ranki

      ranki = iup + (ido - 1)*sector2%up%dimen
   end function

   subroutine rankupdo_stridedo(rank, iup, tabido, sector2)

      ! FINDS THE RANK OF (iup, :) IN FULL DIRECT PRODUCT BASIS

      implicit none

      INTEGER, INTENT(INOUT)                 :: rank(:)
      INTEGER, INTENT(IN)                    :: iup, tabido(:)
      TYPE(fermion_sector2_type), INTENT(IN) :: sector2

      rank = iup + (tabido - 1)*sector2%up%dimen
   end subroutine

   subroutine rankupdo_strideup(rank, tabiup, ido, sector2)

      ! FINDS THE RANK OF (:, ido) IN FULL DIRECT PRODUCT BASIS

      implicit none

      INTEGER, INTENT(INOUT)                 :: rank(:)
      INTEGER, INTENT(IN)                    :: tabiup(:), ido
      TYPE(fermion_sector2_type), INTENT(IN) :: sector2

      rank = tabiup + (ido - 1)*sector2%up%dimen
   end subroutine

   function rankupdochunk(iup, ido, sector2) RESULT(rank)

      ! FINDS THE RANK OF (iup, ido) IN LOCAL CHUNK OF THE DIRECT PRODUCT BASIS

      use genvar, only: iproc

      implicit none

      INTEGER, INTENT(IN)                    :: iup, ido
      TYPE(fermion_sector2_type), INTENT(IN) :: sector2
      INTEGER :: rank
      INTEGER :: iuploc, idoloc

      iuploc = iup - sector2%up%istatemin(iproc) + 1
      idoloc = ido - sector2%down%istatemin(iproc) + 1
      rank = rankupdo(iuploc, idoloc, sector2)
   end function

   function stateupdo(iup, ido, sector2) RESULT(state)

      ! FINDS THE STATE (iup, ido)

      implicit none

      INTEGER, INTENT(IN)                    :: iup, ido
      TYPE(fermion_sector2_type), INTENT(IN) :: sector2
      INTEGER :: state

      state = sector2%up%state(iup) + sector2%down%state(ido)* &
              sector2%up%nstates
   end function

   function rank2(state, sector2)

      ! FIND THE RANK OF STATE IN FULL DIRECT PRODUCT BASIS

      implicit none

      INTEGER, INTENT(IN)                    :: state
      TYPE(fermion_sector2_type), INTENT(IN) :: sector2
      INTEGER :: rank2
      INTEGER :: stateup, statedo

      stateup = MOD(state, sector2%down%nstates)
      statedo = (state - stateup)/sector2%up%nstates
      rank2 = rankupdo(sector2%up%rank(stateup), sector2%down%rank(statedo), &
                       sector2)
   end function

   function state2(rank, sector2)

      ! FIND THE STATE rank IN FULL DIRECT PRODUCT BASIS

      implicit none

      INTEGER, INTENT(IN)                    :: rank
      TYPE(fermion_sector2_type), INTENT(IN) :: sector2
      INTEGER :: state2
      INTEGER :: iup, ido

      iup = MOD(rank - 1, sector2%up%dimen) + 1
      ido = (rank - iup)/sector2%up%dimen + 1
      state2 = stateupdo(iup, ido, sector2)
   end function

   subroutine write_raw_fermion_sector2(sector2, UNIT)

      use fermion_hilbert_class, only: write_raw_fermion_sector

      implicit none

      TYPE(fermion_sector2_type), INTENT(IN) :: sector2
      INTEGER, INTENT(IN)                    :: UNIT

      CALL write_raw_fermion_sector(sector2%up, UNIT)
      CALL write_raw_fermion_sector(sector2%down, UNIT)
   end subroutine

   subroutine read_raw_fermion_sector2(sector2, UNIT, NAMBU)

      use fermion_hilbert_class, only: delete_fermion_sector, &
         fermion_sector_type, read_raw_fermion_sector

      implicit none

      TYPE(fermion_sector2_type), INTENT(INOUT) :: sector2
      INTEGER, INTENT(IN)                       :: UNIT
      LOGICAL, OPTIONAL, INTENT(IN)             :: NAMBU
      TYPE(fermion_sector_type) :: up, down

      CALL delete_fermion_sector2(sector2)
      CALL read_raw_fermion_sector(up, UNIT)
      CALL read_raw_fermion_sector(down, UNIT)
      CALL new_fermion_sector2(sector2, up%npart, down%npart, up%norbs, NAMBU &
                               =NAMBU)
      CALL delete_fermion_sector(up)
      CALL delete_fermion_sector(down)
   end subroutine

end module
