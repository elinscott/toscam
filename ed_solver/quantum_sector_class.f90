MODULE sector_class

   use fermion_hilbert_class, only: fermion_sector_type
   use fermion_sector2_class, only: fermion_sector2_type

   IMPLICIT NONE

   private

   !----------------------!
   ! GENERIC SECTOR TYPE  !
   !----------------------!

   TYPE sector_type
      TYPE(fermion_sector_type),  POINTER :: sz => NULL() ! plain Sz fermionic sector
      TYPE(fermion_sector2_type), POINTER :: updo => NULL() ! (nup, ndo) as a direct product of 2 plain fermionic sectors
   END TYPE

   INTERFACE to_sector
      MODULE PROCEDURE   sz_to_sector
      MODULE PROCEDURE updo_to_sector
   END INTERFACE

   public :: chunk_func
   public :: copy_sector
   public :: delete_sector
   public :: dimen_func
   public :: equal_sector
   public :: is_in_sector
   public :: istatemax__
   public :: istatemin__
   public :: new_sector
   public :: norbs__
   public :: npart_func
   public :: rank_func
   public :: read_raw_sector
   public :: sector_type
   public :: state_func
   public :: title_sector
   public :: write_raw_sector

contains

   function sz_to_sector(sz_sector) RESULT(SEC)

      use fermion_hilbert_class, only: fermion_sector_type

      implicit none

      TYPE(sector_type)                 :: SEC
      TYPE(fermion_sector_type), TARGET :: sz_sector ! new_fermion_sector
      ! sends a pointer here

      SEC%sz => sz_sector
   end function

   function updo_to_sector(updo_sector) RESULT(SEC)

      use fermion_sector2_class, only: fermion_sector2_type

      implicit none

      TYPE(sector_type)                  :: SEC
      TYPE(fermion_sector2_type), TARGET :: updo_sector ! new_fermion_sector2
      ! sends a pointer here

      SEC%updo => updo_sector
   end function

   subroutine new_sector(SECOUT, SECIN)

      use fermion_hilbert_class, only: new_fermion_sector
      use fermion_sector2_class, only: new_fermion_sector2

      implicit none

      TYPE(sector_type), INTENT(INOUT) :: SECOUT
      TYPE(sector_type), INTENT(IN)    :: SECIN

      IF     (ASSOCIATED(SECIN%updo))THEN
         ALLOCATE(SECOUT%updo)
         CALL new_fermion_sector2(SECOUT%updo, SECIN%updo)
      ELSE IF(ASSOCIATED(SECIN%sz))  THEN
         ALLOCATE(SECOUT%sz)
         CALL new_fermion_sector(SECOUT%sz, SECIN%sz)
      ENDIF
   end subroutine

   subroutine copy_sector(SECOUT, SECIN)

      use fermion_hilbert_class, only: copy_fermion_sector
      use fermion_sector2_class, only: copy_fermion_sector2

      implicit none

      TYPE(sector_type), INTENT(INOUT) :: SECOUT
      TYPE(sector_type), INTENT(IN)    :: SECIN

      IF     (ASSOCIATED(SECIN%updo))THEN
         ALLOCATE(SECOUT%updo)
         CALL copy_fermion_sector2(SECOUT%updo, SECIN%updo)
      ELSE IF(ASSOCIATED(SECIN%sz))  THEN
         ALLOCATE(SECOUT%sz)
         CALL copy_fermion_sector(SECOUT%sz, SECIN%sz)
      ENDIF
   end subroutine

   subroutine delete_sector(SEC)

      use fermion_sector2_class, only: delete_fermion_sector2
      use fermion_hilbert_class, only: delete_fermion_sector

      implicit none

      TYPE(sector_type), INTENT(INOUT) :: SEC

      IF     (ASSOCIATED(SEC%updo))THEN
         CALL delete_fermion_sector2(SEC%updo)
         NULLIFY(SEC%updo)
      ELSE IF(ASSOCIATED(SEC%sz))  THEN
         CALL delete_fermion_sector( SEC%sz)
         NULLIFY(SEC%sz)
      ENDIF
   end subroutine

   function is_in_sector(state, SEC)

      use fermion_sector2_class, only: is_in_fermion_sector2
      use fermion_hilbert_class, only: is_in_fermion_sector

      implicit none

      INTEGER, INTENT(IN)           :: state
      TYPE(sector_type), INTENT(IN) :: SEC
      LOGICAL :: is_in_sector

      IF     (ASSOCIATED(SEC%updo))THEN
         is_in_sector = is_in_fermion_sector2(state, SEC%updo)
      ELSE IF(ASSOCIATED(SEC%sz))  THEN
         is_in_sector = is_in_fermion_sector (state, SEC%sz)
      ENDIF
   end function

   function equal_sector(SEC1, SEC2)

      implicit none

      TYPE(sector_type), INTENT(IN) :: SEC1, SEC2
      LOGICAL :: equal_sector

      equal_sector = .false.
      IF (ASSOCIATED(SEC1%updo)) THEN ! fermion_sector2 defined by (nup, ndo,
      ! Ns)
         IF(.NOT.ASSOCIATED(SEC2%updo)) STOP "ERROR IN equal_sector: &
              &INCONSISTENT INPUT SECTORS!"
         IF( SEC1%updo%up%npart == SEC2%updo%up%npart .AND. &
              SEC1%updo%down%npart == SEC2%updo%down%npart.AND. &
              SEC1%updo%up%norbs == SEC2%updo%up%norbs .AND. &
              SEC1%updo%down%norbs == SEC2%updo%down%norbs) equal_sector = &
              .true.
      ELSE IF (ASSOCIATED(SEC1%sz)) THEN ! fermion_sector defined by (npart,
      ! norbs)
         IF(.NOT.ASSOCIATED(SEC2%sz)) STOP "ERROR IN equal_sector: &
              &INCONSISTENT INPUT SECTORS!"
         IF( SEC1%sz%npart == SEC2%sz%npart .AND. SEC1%sz%norbs == &
              SEC2%sz%norbs) equal_sector = .false.
      ENDIF
   end function

   subroutine write_raw_sector(SEC, UNIT)

      use fermion_hilbert_class, only: write_raw_fermion_sector
      use fermion_sector2_class, only: write_raw_fermion_sector2

      implicit none

      TYPE(sector_type), INTENT(IN) :: SEC
      INTEGER, INTENT(IN)           :: UNIT

      IF     (ASSOCIATED(SEC%updo))THEN
         WRITE(UNIT, *) 'updo'
         CALL write_raw_fermion_sector2(SEC%updo, UNIT)
      ELSE IF(ASSOCIATED(SEC%sz))  THEN
         WRITE(UNIT, *) 'sz'
         CALL write_raw_fermion_sector(SEC%sz, UNIT)
      ENDIF
   end subroutine

   subroutine read_raw_sector(SEC, UNIT, SZNAMBU)

      use fermion_hilbert_class, only: fermion_sector_type, &
           read_raw_fermion_sector
      use fermion_sector2_class, only: fermion_sector2_type, &
           read_raw_fermion_sector2

      implicit none

      TYPE(sector_type), INTENT(INOUT) :: SEC
      INTEGER, INTENT(IN)              :: UNIT
      LOGICAL, OPTIONAL, INTENT(IN)    :: SZNAMBU
      CHARACTER(LEN = 100)       :: sectype
      TYPE(fermion_sector2_type) :: updo
      TYPE(fermion_sector_type)  :: sz

      CALL delete_sector(SEC)
      ! READS AN EMPTY SECTOR
      READ(UNIT, *) sectype
      IF(INDEX(sectype, 'updo') /= 0)THEN
         CALL read_raw_fermion_sector2(updo, UNIT, NAMBU = SZNAMBU)
         SEC = updo_to_sector(updo)
      ELSE IF(INDEX(sectype, 'sz') /= 0)THEN
         CALL read_raw_fermion_sector(sz, UNIT, SZ = SZNAMBU)
         SEC =   sz_to_sector(sz)
      ENDIF
   end subroutine

   pure function dimen_func(SEC)

      implicit none

      TYPE(sector_type), INTENT(IN) :: SEC
      INTEGER :: dimen_func

      IF     (ASSOCIATED(SEC%updo))THEN
         dimen_func = SEC%updo%dimen
      ELSE IF(ASSOCIATED(SEC%sz))  THEN
         dimen_func =   SEC%sz%dimen
      ENDIF
   end function

   function title_sector(SEC)

      implicit none

      TYPE(sector_type), INTENT(IN) :: SEC
      CHARACTER(LEN = 100) :: title_sector

      IF     (ASSOCIATED(SEC%updo))THEN
         title_sector = TRIM(SEC%updo%title)
      ELSE IF(ASSOCIATED(SEC%sz))  THEN
         title_sector = TRIM(SEC%sz%title)
      ENDIF
   end function

   function state_func(rank, SEC)

      use fermion_sector2_class, only: state2

      implicit none

      TYPE(sector_type), INTENT(IN) :: SEC
      INTEGER, INTENT(IN)           :: rank
      INTEGER :: state_func

      IF     (ASSOCIATED(SEC%updo))THEN
         state_func = state2(rank, SEC%updo)
      ELSE IF(ASSOCIATED(SEC%sz))  THEN
         state_func = SEC%sz%state(rank)
      ENDIF
   end function

   function rank_func(state, SEC)

      use fermion_sector2_class, only: rank2

      implicit none

      TYPE(sector_type), INTENT(IN) :: SEC
      INTEGER, INTENT(IN)           :: state
      INTEGER :: rank_func

      IF     (ASSOCIATED(SEC%updo))THEN
         rank_func = rank2(state, SEC%updo)
      ELSE IF(ASSOCIATED(SEC%sz))  THEN
         rank_func = SEC%sz%rank(state)
      ENDIF
   end function

   function norbs__(SEC)

      implicit none

      TYPE(sector_type), INTENT(IN) :: SEC
      INTEGER :: norbs__

      IF     (ASSOCIATED(SEC%updo))THEN
         norbs__ = SEC%updo%norbs
      ELSE IF(ASSOCIATED(SEC%sz))  THEN
         norbs__ =   SEC%sz%norbs
      ENDIF
   end function

   function npart_func(SEC, SPIN)

      implicit none

      TYPE(sector_type), INTENT(IN) :: SEC
      INTEGER, OPTIONAL, INTENT(IN) :: SPIN
      INTEGER :: npart_func

      IF     (ASSOCIATED(SEC%updo))THEN
         IF(PRESENT(SPIN))THEN
            IF(SPIN == 1) npart_func = SEC%updo%up%npart
            IF(SPIN == 2) npart_func = SEC%updo%down%npart
         ELSE
            npart_func = SEC%updo%npart
         ENDIF
      ELSE IF(ASSOCIATED(SEC%sz))  THEN
         IF(PRESENT(SPIN)) STOP "ERROR IN npart: SPIN CAN'T BE SPECIFIED IN Sz &
              &SECTOR!"
         npart_func =   SEC%sz%npart
      ENDIF
   end function

   function istatemin__(SEC)

      use genvar, only: iproc

      implicit none

      TYPE(sector_type), INTENT(IN) :: SEC
      INTEGER :: istatemin__

      IF     (ASSOCIATED(SEC%updo))THEN
         istatemin__ = SEC%updo%istatemin(iproc)
      ELSE IF(ASSOCIATED(SEC%sz))  THEN
         istatemin__ =   SEC%sz%istatemin(iproc)
      ENDIF
   end function

   function istatemax__(SEC)

      use genvar, only: iproc

      implicit none

      TYPE(sector_type), INTENT(IN) :: SEC
      INTEGER :: istatemax__

      IF     (ASSOCIATED(SEC%updo))THEN
         istatemax__ = SEC%updo%istatemax(iproc)
      ELSE IF(ASSOCIATED(SEC%sz))  THEN
         istatemax__ =   SEC%sz%istatemax(iproc)
      ENDIF
   end function

   function chunk_func(SEC)

      use genvar, only: iproc

      implicit none

      TYPE(sector_type), INTENT(IN) :: SEC
      INTEGER :: chunk_func

      IF     (ASSOCIATED(SEC%updo))THEN
         chunk_func = SEC%updo%chunk(iproc)
      ELSE IF(ASSOCIATED(SEC%sz))  THEN
         chunk_func =   SEC%sz%chunk(iproc)
      ENDIF
   end function

end module
