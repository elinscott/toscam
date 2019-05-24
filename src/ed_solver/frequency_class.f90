MODULE frequency_class

   use genvar, only: DBL

   !$$$$$$$$$$$$$$$$$$$$$
   !$$FREQUENCY CLASS$$
   !$$$$$$$$$$$$$$$$$$$$$

   implicit none

   private

   INTEGER, PRIVATE                                  :: istati

   TYPE freq_type
      INTEGER               :: Nw = 0
      INTEGER               :: iwprint = 0
      CHARACTER(LEN=9)    :: title = '\0' ! FERMIONIC OR BOSONIC OR RETARDED OR ADVANCED
      REAL(DBL)             :: beta = 0  ! IF MATSUBARA
      REAL(DBL)             :: wmax = 0  ! IF REAL FREQ.
      REAL(DBL)             :: wmin = 0  ! IF REAL FREQ.
      REAL(DBL)             :: width = 0  ! IF REAL FREQ.
      COMPLEX(DBL), POINTER :: vec(:) => NULL() ! frequency array
      INTEGER, POINTER :: iwmin(:) => NULL(), iwmax(:) => NULL(), chunk(:) => NULL()
   END TYPE

   INTERFACE new_freq
      MODULE PROCEDURE new_rfreq_from_scratch
      MODULE PROCEDURE new_ifreq_from_scratch
      MODULE PROCEDURE new_freq_from_old
   END INTERFACE

   REAL(8), PRIVATE :: rdelta_frequ_eta1, rdelta_frequ_w0, rdelta_frequ_T, &
                       rdelta_frequ_eta2

   public :: copy_frequency
   public :: delete_frequency
   public :: freq_type
   public :: init_frequency_module
   public :: new_freq
   public :: read_raw_freq
   public :: write_raw_freq

contains

   subroutine init_frequency_module(rdelta_frequ_eta1_, rdelta_frequ_w0_, &
                                    rdelta_frequ_T_, rdelta_frequ_eta2_)

      implicit none

      real(8)           :: rdelta_frequ_eta1_, rdelta_frequ_w0_, rdelta_frequ_T_
      real(8), optional :: rdelta_frequ_eta2_

      rdelta_frequ_eta1 = rdelta_frequ_eta1_
      rdelta_frequ_w0 = rdelta_frequ_w0_
      rdelta_frequ_T = rdelta_frequ_T_
      if (present(rdelta_frequ_eta2_)) then
         rdelta_frequ_eta2 = rdelta_frequ_eta2_
      else
         rdelta_frequ_eta2 = 0.d0
      endif

   end subroutine

   subroutine new_rfreq_from_scratch(FREQ, Nw, wmin, wmax, width, sense)

      use genvar, only: ADVANCED, imi, iwprint, RETARDED
      use linalg, only: dexpc, ramp

      implicit none

      TYPE(freq_type), INTENT(INOUT) :: FREQ
      INTEGER, INTENT(IN)            :: Nw
      REAL(DBL), INTENT(IN)          :: width, wmax, wmin
      CHARACTER(LEN=9), INTENT(IN) :: sense
      REAL(DBL) :: dw
      INTEGER   :: iw, ramp_freq(Nw), i

      if (Nw == 0) stop 'error new_rfreq_from_scratch: 0 frequency'

      CALL delete_frequency(FREQ)

      !real freq. step

      dw = (wmax - wmin)/dble(Nw - 1)
      FREQ%title = sense
      ALLOCATE (FREQ%vec(Nw))
      CALL ramp(ramp_freq)
      FREQ%vec = wmin + dble(ramp_freq - 1)*dw

      SELECT CASE (sense)
      CASE (RETARDED)
         do i = 1, Nw
            FREQ%vec(i) = FREQ%vec(i) + imi*(width + rdelta_frequ_eta1/( &
                                             DEXPc(-(real(FREQ%vec(i)) - rdelta_frequ_w0)/rdelta_frequ_T) + &
                                             1.d0) + rdelta_frequ_eta1/(DEXPc((real(FREQ%vec(i)) + &
                                                                               rdelta_frequ_w0)/rdelta_frequ_T) + 1.d0))
            if (abs(rdelta_frequ_eta2) > 1.d-3) then
               FREQ%vec(i) = FREQ%vec(i) + imi*(rdelta_frequ_eta2* &
                                                rdelta_frequ_eta1/(DEXPc(-(real(FREQ%vec(i)) - &
                                                                           rdelta_frequ_eta2*rdelta_frequ_w0)/(rdelta_frequ_eta2* &
                                                                                     rdelta_frequ_T)) + 1.d0) + rdelta_frequ_eta2* &
                                                rdelta_frequ_eta1/(DEXPc((real(FREQ%vec(i)) + &
                                                                          rdelta_frequ_eta2*rdelta_frequ_w0)/(rdelta_frequ_eta2 &
                                                                                                          *rdelta_frequ_T)) + 1.d0))
            endif
         enddo
      CASE (ADVANCED)
         do i = 1, Nw
            FREQ%vec(i) = FREQ%vec(i) - imi*(width + rdelta_frequ_eta1/( &
                                             DEXPc(-(real(FREQ%vec(i)) - rdelta_frequ_w0)/rdelta_frequ_T) + &
                                             1.d0) + rdelta_frequ_eta1/(DEXPc((real(FREQ%vec(i)) + &
                                                                               rdelta_frequ_w0)/rdelta_frequ_T) + 1.d0))
            if (abs(rdelta_frequ_eta2) > 1.d-3) then
               FREQ%vec(i) = FREQ%vec(i) - imi*(rdelta_frequ_eta2* &
                                                rdelta_frequ_eta1/(DEXPc(-(real(FREQ%vec(i)) - &
                                                                           rdelta_frequ_eta2*rdelta_frequ_w0)/(rdelta_frequ_eta2 &
                                                                                    *rdelta_frequ_T)) + 1.d0) + rdelta_frequ_eta2* &
                                                rdelta_frequ_eta1/(DEXPc((real(FREQ%vec(i)) + &
                                                                          rdelta_frequ_eta2*rdelta_frequ_w0)/(rdelta_frequ_eta2 &
                                                                                                          *rdelta_frequ_T)) + 1.d0))
            endif
         enddo
      END SELECT

      FREQ%iwprint = MIN(iwprint, Nw)
      FREQ%Nw = Nw
      FREQ%wmin = wmin
      FREQ%wmax = wmax
      FREQ%width = width

      CALL split_frequency(FREQ)

      return
   end subroutine

   subroutine new_ifreq_from_scratch(FREQ, Nw, beta, stat)

      use genvar, only: BOSONIC, FERMIONIC, imi, iwprint, pi, pi2
      use linalg, only: ramp

      implicit none

      TYPE(freq_type), INTENT(INOUT) :: FREQ
      INTEGER, INTENT(IN)            :: Nw
      REAL(DBL), INTENT(IN)          :: beta
      CHARACTER(LEN=9), INTENT(IN) :: stat
      INTEGER :: iw, ramp_freq(Nw)

      if (Nw == 0) stop 'error new_ifreq_from_scratch: 0 frequency'

      CALL delete_frequency(FREQ)

      FREQ%title = stat
      ALLOCATE (FREQ%vec(Nw))
      CALL ramp(ramp_freq)
      SELECT CASE (stat)
      CASE (BOSONIC)
         FREQ%vec = imi*pi2/beta*ramp_freq
      CASE (FERMIONIC)
         FREQ%vec = imi*pi2/beta*ramp_freq - imi*pi/beta
      END SELECT
      FREQ%iwprint = MIN(iwprint, Nw)
      FREQ%Nw = Nw
      FREQ%beta = beta
      CALL split_frequency(FREQ)

   end subroutine

   subroutine new_freq_from_old(FREQOUT, FREQIN)

      use genvar, only: ADVANCED, BOSONIC, FERMIONIC, RETARDED

      implicit none

      TYPE(freq_type), INTENT(INOUT) :: FREQOUT
      TYPE(freq_type), INTENT(IN)    :: FREQIN

      CALL delete_frequency(FREQOUT)

      SELECT CASE (FREQIN%title)
      CASE (RETARDED, ADVANCED)
         CALL new_rfreq_from_scratch(FREQOUT, FREQIN%Nw, FREQIN%wmin, &
                                     FREQIN%wmax, FREQIN%width, FREQIN%TITLE)
      CASE (FERMIONIC, BOSONIC)
         CALL new_ifreq_from_scratch(FREQOUT, FREQIN%Nw, FREQIN%beta, &
                                     FREQIN%TITLE)
      END SELECT

      CALL copy_frequency(FREQOUT, FREQIN)

   end subroutine

   subroutine split_frequency(FREQ)

      use genvar, only: nproc
      use mpi_mod, only: split

      implicit none

      TYPE(freq_type), INTENT(INOUT) :: FREQ

      ! DISTRIBUTE THE FREQUENCIES OVER nproc >= 1 PROCESSORS

      IF (ASSOCIATED(FREQ%chunk)) DEALLOCATE (FREQ%chunk, STAT=istati)
      IF (ASSOCIATED(FREQ%iwmin)) DEALLOCATE (FREQ%iwmin, STAT=istati)
      IF (ASSOCIATED(FREQ%iwmax)) DEALLOCATE (FREQ%iwmax, STAT=istati)
      ALLOCATE (FREQ%chunk(nproc), FREQ%iwmin(nproc), FREQ%iwmax(nproc))
      if (FREQ%Nw <= 0) stop 'error split_frequency, Nw = 0'
      CALL split(FREQ%Nw, FREQ%iwmin, FREQ%iwmax, FREQ%chunk)

   end subroutine

   subroutine copy_frequency(FREQOUT, FREQIN)

      implicit none

      TYPE(freq_type), INTENT(INOUT) :: FREQOUT
      TYPE(freq_type), INTENT(IN)    :: FREQIN

      if (.not. ASSOCIATED(FREQOUT%vec)) ALLOCATE (FREQOUT%vec(FREQIN%Nw))
      if (.not. ASSOCIATED(FREQOUT%chunk)) &
         ALLOCATE (FREQOUT%chunk(size(FREQIN%chunk)))
      if (.not. ASSOCIATED(FREQOUT%iwmin)) &
         ALLOCATE (FREQOUT%iwmin(size(FREQIN%iwmin)))
      if (.not. ASSOCIATED(FREQOUT%iwmax)) &
         ALLOCATE (FREQOUT%iwmax(size(FREQIN%iwmax)))

      FREQOUT%title = FREQIN%title
      FREQOUT%iwprint = FREQIN%iwprint
      FREQOUT%Nw = FREQIN%Nw
      FREQOUT%width = FREQIN%width
      FREQOUT%wmin = FREQIN%wmin
      FREQOUT%wmax = FREQIN%wmax
      FREQOUT%beta = FREQIN%beta
      FREQOUT%vec = FREQIN%vec
      FREQOUT%iwmin = FREQIN%iwmin
      FREQOUT%iwmax = FREQIN%iwmax
      FREQOUT%chunk = FREQIN%chunk

   end subroutine

   subroutine delete_frequency(FREQ)

      use genvar, only: log_unit, messages3

      implicit none

      TYPE(freq_type), INTENT(INOUT) :: FREQ

      if (messages3) write (log_unit, *) '--- > delete frequency vec < ---'
      IF (ASSOCIATED(FREQ%vec)) DEALLOCATE (FREQ%vec, STAT=istati)
      if (messages3) write (log_unit, *) 'delete frequency iwmin'
      IF (ASSOCIATED(FREQ%iwmin)) DEALLOCATE (FREQ%iwmin, STAT=istati)
      if (messages3) write (log_unit, *) 'delete frequency iwmax'
      IF (ASSOCIATED(FREQ%iwmax)) DEALLOCATE (FREQ%iwmax, STAT=istati)
      if (messages3) write (log_unit, *) 'delete frequency chunk'
      IF (ASSOCIATED(FREQ%chunk)) DEALLOCATE (FREQ%chunk, STAT=istati)
      if (messages3) write (log_unit, *) 'done return......'
   end subroutine

   subroutine write_raw_freq(FREQ, UNIT)

      ! WRITE FREQUENCIES TO FILE IN RAW (UNFORMATTED) FORM

      use genvar, only: ADVANCED, BOSONIC, FERMIONIC, pi2, RETARDED

      implicit none

      TYPE(freq_type), INTENT(IN) :: FREQ
      INTEGER, INTENT(IN)         :: UNIT

      WRITE (UNIT, *) FREQ%title
      WRITE (UNIT, *) FREQ%Nw
      SELECT CASE (FREQ%title)
      CASE (FERMIONIC, BOSONIC) ! MATSUBARA FREQ.
         WRITE (UNIT, *) pi2/(AIMAG(FREQ%vec(2) - FREQ%vec(1)))
      CASE (RETARDED, ADVANCED) ! REAL FREQ.
         WRITE (UNIT, *) FREQ%wmin, FREQ%wmax
         WRITE (UNIT, *) ABS(AIMAG(FREQ%vec(1)))
      END SELECT
      WRITE (UNIT, *) FREQ%vec

   end subroutine

   subroutine read_raw_freq(FREQ, UNIT)

      ! CREATE NEW CORRELATION FROM RAW FILE

      use genvar, only: ADVANCED, BOSONIC, FERMIONIC, RETARDED

      implicit none

      TYPE(freq_type), INTENT(INOUT) :: FREQ
      INTEGER, INTENT(IN)            :: UNIT
      CHARACTER(LEN=9) :: title
      INTEGER            :: Nw
      REAL(DBL)          :: wmin, wmax, width, beta

      CALL delete_frequency(FREQ)
      ! READ RAW CORRELATION FROM INPUT FILE
      READ (UNIT, *) title
      READ (UNIT, *) Nw
      if (Nw == 0) stop 'error read raw frequency: 0 frequency'

      SELECT CASE (title)
      CASE (FERMIONIC, BOSONIC) ! MATSUBARA FREQ.
         READ (UNIT, *) beta
         CALL new_ifreq_from_scratch(FREQ, Nw, beta, title)
      CASE (RETARDED, ADVANCED) ! REAL FREQ.
         READ (UNIT, *) wmin, wmax
         READ (UNIT, *) width
         CALL new_rfreq_from_scratch(FREQ, Nw, wmin, wmax, width, title)
      END SELECT
      READ (UNIT, *) FREQ%vec
   end subroutine

end module
