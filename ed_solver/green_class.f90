MODULE green_class

   use genvar,              only: DBL, log_unit
   use correl_class,        only: correl_type
   use masked_matrix_class, only: masked_matrix_type

   IMPLICIT NONE

   private

   LOGICAL, PARAMETER :: ALLOCATE_ALL = .true.

   TYPE green_type
      ! DYNAMICAL CORRELATIONS OF OPERATORS A, B
      INTEGER                  :: N  = 0
      INTEGER                  :: Nw = 0
      CHARACTER(LEN = 100)     :: title = '\0'
      LOGICAL                  :: compute(2, 2) = RESHAPE((/.false., .false., &
                                  .false., .false./), (/2, 2/)) ! independent green's fctns
      TYPE(correl_type)        :: correl (2, 2) ! green's fctn in all channels
      TYPE(masked_matrix_type) :: correlstat(2, 2) ! static green's fctn in
      ! all channels particle, hole
#ifdef _complex
      ! mean of part & hole operators
      COMPLEX(DBL), POINTER    :: Amean(:, :) => NULL(), Bmean(:, :) => NULL()
#else
      ! mean of part & hole operators
      REAL(DBL),    POINTER    :: Amean(:, :) => NULL(), Bmean(:, :) => NULL()
#endif
   END TYPE

   INTERFACE new_green
      MODULE PROCEDURE new_green_rfreq
      MODULE PROCEDURE new_green_ifreq
   END INTERFACE

   public :: copy_green
   public :: green_type
   public :: new_green
   public :: pad_green
   public :: write_green

contains

   subroutine new_green_rfreq(GREEN, compute, title, N, Nw, wmin, wmax, width, &
        STAT, IMASK, AB, force_compute)

      use genvar,              only: dbl, log_unit, messages3, pm
      use correl_class,        only: correl2vec, new_correl
      use masked_matrix_class, only: masked_matrix2vec, new_masked_matrix

      implicit none

      TYPE(green_type), INTENT(INOUT) :: GREEN
      LOGICAL, INTENT(IN)             :: compute(2, 2)
      CHARACTER(LEN = *), INTENT(IN)  :: title
      INTEGER, INTENT(IN)             :: N
      CHARACTER(LEN = 9), INTENT(IN)  :: STAT
      INTEGER, INTENT(IN)             :: Nw
      REAL(DBL), INTENT(IN)           :: width, wmin, wmax
      INTEGER, OPTIONAL, INTENT(IN)   :: IMASK(N, N)
      LOGICAL, OPTIONAL, INTENT(IN)   :: AB
      INTEGER           :: ipm, jpm
      LOGICAL           :: IS_HERM, AB_
      LOGICAL, OPTIONAL :: force_compute(2, 2)


      AB_ = .false.
      IF(PRESENT(AB)) AB_ = AB ! TRUE IF < A(z)B >, FALSE IF < A(z)A >

      if(messages3) write(log_unit, *) 'delete green function'
      CALL delete_green(GREEN)

      if(messages3) write(log_unit, *) 'done, now define main variables...'

      GREEN%title = TRIM(title(1:MIN(LEN_TRIM(title), 100)))
      GREEN%N     = N
      GREEN%Nw    = Nw

      ! ELIMINATE REDUNDANCIES => DEDUCE THEM BY SYMMETRY: < a(z)b >= <
      ! A(-z*)B > *

      IF(compute(1, 1) .OR. compute(2, 2))THEN ! COMPUTE ONLY < a(z)*b > BY
      ! DEFAULT
         GREEN%compute(1, 1) = .false.
         GREEN%compute(2, 2) = .true.
      ENDIF
      IF(compute(1, 2) .OR. compute(2, 1))THEN ! COMPUTE ONLY < a(z)*B > BY
      ! DEFAULT
         GREEN%compute(1, 2) = .false.
         GREEN%compute(2, 1) = .true.
      ENDIF

      if(present(force_compute)) GREEN%compute = force_compute

      if(messages3) write(log_unit, *) 'build new arrays'

      DO ipm = 1, 2
         DO jpm = 1, 2
#ifdef _complex
            IS_HERM = .false.
#else
            IS_HERM = .false.
            IF(ipm /= jpm .AND. (.NOT.AB_)) IS_HERM = .true.
#endif
            ! EQUAL-TIME
            CALL new_masked_matrix(GREEN%correlstat(ipm, jpm), &
                 TRIM(GREEN%title) // pm(ipm) // pm(jpm) // "_stat", N, N, &
                 IMASK = IMASK, IS_HERM = IS_HERM)
            CALL masked_matrix2vec(GREEN%correlstat(ipm, jpm)) ! dummy to
            ! create vec of ind elemts

            ! DYNAMIC
            IF(GREEN%compute(ipm, jpm) .OR. GREEN%compute(3-ipm, 3-jpm) .OR. &
                 ALLOCATE_ALL)THEN
               CALL new_correl(GREEN%correl(ipm, jpm), TRIM(GREEN%title) // &
                    pm(ipm) // pm(jpm), N, Nw, wmin, wmax, width, STAT, IMASK &
                    = IMASK, AB = AB)
               CALL correl2vec(GREEN%correl(ipm, jpm)) ! dummy to create vec
               ! of ind elemts
            ENDIF

         ENDDO
      ENDDO

      if(messages3) write(log_unit, *) '...new frequ rfrequ....'
      IF(N /= 0)THEN
         if(associated(GREEN%Amean)) deallocate(GREEN%Amean)
         if(associated(GREEN%Bmean)) deallocate(GREEN%Bmean)
         ALLOCATE(GREEN%Amean(N, 2), GREEN%Bmean(N, 2))
         GREEN%Amean = 0.0_DBL
         GREEN%Bmean = 0.0_DBL
      ENDIF
      if(messages3) write(log_unit, *) ' .. allocation  done .. '

   end subroutine

   subroutine new_green_ifreq(GREEN, compute, title, N, Nw, beta, STAT, IMASK, &
        AB, force_compute)

      use genvar,              only: dbl, log_unit, messages3, pm
      use correl_class,        only: correl2vec, new_correl
      use masked_matrix_class, only: masked_matrix2vec, new_masked_matrix

      implicit none

      TYPE(green_type), INTENT(INOUT) :: GREEN
      LOGICAL, INTENT(IN)             :: compute(2, 2)
      CHARACTER(LEN = *), INTENT(IN)  :: title
      INTEGER, INTENT(IN)             :: N
      CHARACTER(LEN = 9), INTENT(IN)  :: STAT
      INTEGER, INTENT(IN)             :: Nw
      REAL(DBL), INTENT(IN)           :: beta
      INTEGER, OPTIONAL, INTENT(IN)   :: IMASK(N, N)
      LOGICAL, OPTIONAL, INTENT(IN)   :: AB
      INTEGER           :: ipm, jpm
      LOGICAL           :: IS_HERM, AB_
      LOGICAL, OPTIONAL :: force_compute(2, 2)

      AB_ = .false.
      IF(PRESENT(AB)) AB_ = AB ! TRUE IF < A(z)B >, FALSE IF < A(z)A >

      if(messages3) write(log_unit, *) 'delete green function'

      CALL delete_green(GREEN)

      GREEN%title   = TRIM(title(1:MIN(LEN_TRIM(title), 100)))
      GREEN%N       = N
      GREEN%Nw      = Nw

      ! ELIMINATE REDUNDANCIES => DEDUCE THEM BY SYMMETRY: < a(z)b >= <
      ! A(-z*)B > *

      IF(compute(1, 1) .OR. compute(2, 2))THEN ! COMPUTE ONLY < a(z)*b > BY
      ! DEFAULT
         GREEN%compute(1, 1) = .false.
         GREEN%compute(2, 2) = .true.
      ENDIF
      IF(compute(1, 2) .OR. compute(2, 1))THEN ! COMPUTE ONLY < a(z)*B > BY
      ! DEFAULT
         GREEN%compute(1, 2) = .false.
         GREEN%compute(2, 1) = .true.
      ENDIF

      if(present(force_compute)) GREEN%compute = force_compute

      if(messages3) write(log_unit, *) 'build new arrays'

      DO ipm = 1, 2
         DO jpm = 1, 2
#ifdef _complex
            IS_HERM = .false.
#else
            IS_HERM = .false.
            IF(ipm /= jpm .AND. (.NOT.AB_)) IS_HERM = .true.
#endif
            ! EQUAL-TIME
            CALL new_masked_matrix(GREEN%correlstat(ipm, jpm), &
                 TRIM(GREEN%title) // pm(ipm) // pm(jpm) // "_stat", N, N, &
                 IMASK = IMASK, IS_HERM = IS_HERM)
            CALL masked_matrix2vec(GREEN%correlstat(ipm, jpm)) ! dummy to
            ! create vec of ind elemts

            ! DYNAMIC
            IF(GREEN%compute(ipm, jpm) .OR. GREEN%compute(3-ipm, 3-jpm) .OR. &
                 ALLOCATE_ALL)THEN
               CALL new_correl(GREEN%correl(ipm, jpm), TRIM(GREEN%title) // &
                    pm(ipm) // pm(jpm), N, Nw, beta, STAT, IMASK = IMASK, AB = &
                    AB)
               CALL correl2vec(GREEN%correl(ipm, jpm)) ! dummy to create vec
               ! of ind elemts
            ENDIF

         ENDDO
      ENDDO

      if(messages3) write(log_unit, *) '...new frequ rfrequ....'
      IF(N /= 0)THEN
         if(associated(GREEN%Amean)) deallocate(GREEN%Amean)
         if(associated(GREEN%Bmean)) deallocate(GREEN%Bmean)
         ALLOCATE(GREEN%Amean(N, 2), GREEN%Bmean(N, 2))
         GREEN%Amean = 0.0_DBL
         GREEN%Bmean = 0.0_DBL
      ENDIF
      if(messages3) write(log_unit, *) '...allocation done...'

   end subroutine

   subroutine delete_green(GREEN)

      use correl_class,        only: delete_correl
      use genvar,              only: log_unit, messages3
      use masked_matrix_class, only: delete_masked_matrix

      implicit none

      TYPE(green_type), INTENT(INOUT) :: GREEN
      INTEGER :: ipm, jpm

      if(messages3) write(log_unit, *) '...delete correlations contained in &
           &green...'

      DO ipm = 1, 2
         DO jpm = 1, 2
            if(messages3) write(log_unit, *) 'delete dynamical part'
            IF(GREEN%compute(ipm, jpm) .OR. GREEN%compute(3-ipm, 3-jpm)) CALL &
                 delete_correl(GREEN%correl(ipm, jpm))
            if(messages3) write(log_unit, *) 'delete equal time correlations'
            CALL delete_masked_matrix(GREEN%correlstat(ipm, jpm))
            if(messages3) write(log_unit, *) 'done......'
         ENDDO
      ENDDO

      if(messages3) write(log_unit, *) '...delete green function s averages...'
      IF(ASSOCIATED(GREEN%Amean)) DEALLOCATE(GREEN%Amean)
      IF(ASSOCIATED(GREEN%Bmean)) DEALLOCATE(GREEN%Bmean)

   end subroutine

   subroutine copy_green(GREENOUT, GREENIN)

      use genvar,              only: log_unit, messages3
      use masked_matrix_class, only: copy_masked_matrix
      use correl_class,        only: copy_correl

      implicit none

      TYPE(green_type), INTENT(IN)    :: GREENIN
      TYPE(green_type), INTENT(INOUT) :: GREENOUT
      INTEGER :: ipm, jpm

      GREENOUT%title   = GREENIN%title
      GREENOUT%N       = GREENIN%N
      GREENOUT%Nw      = GREENIN%Nw
      GREENOUT%compute = GREENIN%compute

      DO ipm = 1, 2
         DO jpm = 1, 2
            ! DYNAMIC
            if(messages3) write(log_unit, *) ' ----- > copy_green : correl '
            IF(GREENOUT%compute(ipm, jpm) .OR. GREENOUT%compute(3-ipm, 3-jpm)) &
                 then
               CALL copy_correl(GREENOUT%correl(ipm, jpm), GREENIN%correl(ipm, &
                    jpm))
            ENDIF
            if(messages3) write(log_unit, *) ' ----- > copy_green : mask '
            ! EQUAL-TIME
            CALL copy_masked_matrix(GREENOUT%correlstat(ipm, jpm), &
                 GREENIN%correlstat(ipm, jpm))
         ENDDO
      ENDDO

      GREENOUT%Amean = GREENIN%Amean
      GREENOUT%Bmean = GREENIN%Bmean

   end subroutine

   subroutine pad_green(GREENOUT, GREENIN)

      ! PAD REDUNDANT MATRIX ELEMENTS IF NECESSARY
      ! USING < a(z)*b >= (-1)^F < A(-z*)*B > *
      ! FIRST PAD IN (ipm, jpm) SECTOR USING EXTERNAL GREEN'S FUNCTION IF
      ! NECESSARY

      use masked_matrix_class, only: pad_masked_matrix
      use correl_class,        only: pad_correl

      implicit none

      TYPE(green_type), INTENT(INOUT)        :: GREENOUT
      TYPE(green_type), INTENT(IN), OPTIONAL :: GREENIN
      INTEGER :: ipm, jpm, mipm, mjpm, iw, miw

#ifdef DEBUG
      write(log_unit, *) "DEBUG: entering green_class_pad_green"
#endif

      DO ipm = 1, 2
         DO jpm = 1, 2
            IF(PRESENT(GREENIN))THEN
               GREENOUT%Amean = GREENIN%Amean ! WARNING: OVERRIDE!... SO PAD
               ! GREENOUT WITH GREENIN BEFORE COMPUTING GREENIN!
               GREENOUT%Bmean = GREENIN%Bmean ! WARNING: OVERRIDE!... SO PAD
               ! GREENOUT WITH GREENIN BEFORE COMPUTING GREENIN!
               ! EQUAL-TIME
               CALL pad_masked_matrix(GREENOUT%correlstat(ipm, jpm), &
                    GREENIN%correlstat(ipm, jpm))
               ! DYNAMIC
               IF( (GREENOUT%compute(ipm, jpm) .OR. GREENOUT%compute(3-ipm, &
                    3-jpm)) .AND. ( GREENIN%compute(ipm, jpm) .OR. &
                    GREENIN%compute(3-ipm, 3-jpm))) CALL &
                    pad_correl(GREENOUT%correl(ipm, jpm), GREENIN%correl(ipm, &
                    jpm))
            ELSE
               ! EQUAL-TIME
               CALL pad_masked_matrix(GREENOUT%correlstat(ipm, jpm))
               ! DYNAMIC
               IF( (GREENOUT%compute(ipm, jpm) .OR. GREENOUT%compute(3-ipm, &
                    3-jpm)) .AND. (GREENIN%compute(ipm, jpm) .OR. &
                    GREENIN%compute(3-ipm, 3-jpm))) CALL &
                    pad_correl(GREENOUT%correl(ipm, jpm))
            ENDIF
         ENDDO
      ENDDO

#ifdef DEBUG
      write(log_unit, *) "DEBUG: leaving green_class_pad_green"
#endif

   end subroutine

   subroutine write_green(GREEN)

      use correl_class, only: write_correl

      implicit none

      TYPE(green_type), INTENT(IN) :: GREEN
      INTEGER :: ipm, jpm

      DO ipm = 1, 2
         DO jpm = 1, 2
            IF(GREEN%compute(ipm, jpm) .OR. GREEN%compute(3-ipm, 3-jpm))THEN
               CALL write_correl(GREEN%correl(ipm, jpm))
            ENDIF
         ENDDO
      ENDDO
   end subroutine

   subroutine glimpse_green(GREEN, UNIT)

      use correl_class, only: glimpse_correl

      implicit none

      TYPE(green_type), INTENT(IN) :: GREEN
      INTEGER, OPTIONAL :: UNIT
      INTEGER           :: ipm, jpm

      DO ipm = 1, 2
         DO jpm = 1, 2
            IF(GREEN%compute(ipm, jpm) .OR. GREEN%compute(3-ipm, 3-jpm))THEN
               CALL glimpse_correl(GREEN%correl(ipm, jpm), UNIT)
            ENDIF
         ENDDO
      ENDDO
   end subroutine

   subroutine slice_green(GREENOUT, GREENIN, rbounds, cbounds)

      use common_def,          only: c2s, i2c
      use correl_class,        only: slice_correl
      use masked_matrix_class, only: slice_masked_matrix

      implicit none

      TYPE(green_type), INTENT(INOUT) :: GREENOUT
      TYPE(green_type), INTENT(IN)    :: GREENIN
      INTEGER, INTENT(IN)             :: rbounds(2), cbounds(2)
      INTEGER :: ipm, jpm, n1slice, n2slice

      CALL delete_green(GREENOUT)
      n1slice = rbounds(2)-rbounds(1) + 1
      n2slice = cbounds(2)-cbounds(1) + 1

      IF(n1slice /= n2slice) STOP "ERROR IN slice_green: SQUARE SLICE &
           &REQUIRED!"
      IF(n1slice <= 0) STOP "ERROR IN slice_green: NON-NULL DIMENSIONS &
           &REQUIRED!"

      GREENOUT%title = TRIM(GREENIN%title) // "_" // c2s(i2c(rbounds(1))) // &
           "_" // c2s(i2c(rbounds(2))) // "_" // c2s(i2c(cbounds(1))) // "_" &
           // c2s(i2c(cbounds(2)))
      GREENOUT%compute  = GREENIN%compute
      GREENOUT%Nw       = GREENIN%Nw
      GREENOUT%N        = n1slice

      DO ipm = 1, 2
         DO jpm = 1, 2
            IF(GREENIN%compute(ipm, jpm) .OR. GREENIN%compute(3-ipm, 3-jpm))THEN
               ! DYNAMIC
               CALL slice_correl(GREENOUT%correl(ipm, jpm), &
                    GREENIN%correl(ipm, jpm), rbounds, cbounds)
               ! EQUAL-TIME
               CALL slice_masked_matrix(GREENOUT%correlstat(ipm, jpm), &
                    GREENIN%correlstat(ipm, jpm), rbounds, cbounds)
            ENDIF
         ENDDO
      ENDDO

   end subroutine

end module
