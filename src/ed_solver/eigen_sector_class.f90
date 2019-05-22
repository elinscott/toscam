MODULE eigen_sector_class

   use eigen_class, only: eigenlist_type
   use genvar, only: DBL, DP, log_unit
   use sector_class, only: sector_type

   implicit none

   private

   type eigensector_type
      type(sector_type)    :: sector
      type(eigenlist_type) :: lowest
   end type

   interface new_eigensector
      module procedure copy_eigensector
      module procedure new_eigensector_from_scratch
   end interface

   type eigensectorlist_type
      integer                         :: nsector = 0
      type(eigensector_type), pointer :: es(:) => null()
   end type

   interface new_rcvec_from_file
      module procedure new_rcvec_from_file_c, new_rcvec_from_file_r
   end interface

   public :: add_eigensector
   public :: copy_eigensector
   public :: copy_eigensectorlist
   public :: delete_eigensector
   public :: delete_eigensectorlist
   public :: eigensector_type
   public :: eigensectorlist_type
   public :: filter_eigensector
   public :: gsenergy
   public :: nambu_energy_shift
   public :: new_eigensector
   public :: not_commensurate_sector
   public :: not_commensurate_sector_
   public :: nlowest
   public :: partition
   public :: print_eigensectorlist
   public :: rank_sector_in_list
   public :: read_raw_eigensectorlist
   public :: write_raw_eigensectorlist

contains

   subroutine new_rcvec_from_file_r(unit, vec, sector)

      use sector_class, only: rank_func, sector_type

      implicit none

      integer           :: unit, istate, jstate, stateup, statedo
      TYPE(sector_type) :: sector
      real(8)           :: coeff, vec(:)

      DO istate = 1, size(vec)
         READ (unit, *) stateup, statedo, coeff
         jstate = stateup + IBITS(NOT(statedo), 0, 12)*(2**12)
         vec(rank_func(jstate, sector)) = coeff
      ENDDO
   end subroutine

   subroutine new_rcvec_from_file_c(unit, vec, sector)

      use sector_class, only: rank_func, sector_type

      implicit none

      integer           :: unit, istate, jstate, stateup, statedo
      TYPE(sector_type) :: sector
      complex(8)        :: coeff, vec(:)

      DO istate = 1, size(vec)
         READ (unit, *) stateup, statedo, coeff
         jstate = stateup + IBITS(NOT(statedo), 0, 12)*(2**12)
         vec(rank_func(jstate, sector)) = coeff
      ENDDO
   end subroutine

   subroutine new_eigensectorlist(list, nsector)

      use genvar, only: log_unit, messages3

      implicit none

      TYPE(eigensectorlist_type), INTENT(INOUT) :: list
      INTEGER, INTENT(IN)                       :: nsector

      if (messages3) write (log_unit, *) 'delete eigensector'
      CALL delete_eigensectorlist(list)
      list%nsector = nsector
      if (messages3) write (log_unit, *) 'new eigensector list nsector --- > ', &
         nsector

      IF (list%nsector > 0) then
         if (associated(list%es)) deallocate (list%es)
         if (messages3) write (log_unit, *) 'allocate new list....'
         ALLOCATE (list%es(nsector))
      ELSE
         if (messages3) write (log_unit, *) 'create list with 0 eigensector'
         if (associated(list%es)) deallocate (list%es)
         list%es => NULL()
      ENDIF

   end subroutine

   subroutine remove_eigensector(sector, list)

      use sector_class, only: sector_type, title_sector

      implicit none

      TYPE(eigensectorlist_type) :: list
      TYPE(sector_type)          :: sector
      INTEGER                    :: isector, jsector, truerank
      TYPE(eigensectorlist_type) :: tmp

      truerank = rank_sector_in_list(sector, list)
      IF (truerank == 0) then
         write (*, *) "ERROR IN remove_eigensector: SECTOR "// &
            TRIM(title_sector(sector))//" ISNT IN LIST!"
         stop
      endif
      CALL copy_eigensectorlist(tmp, list)
      CALL new_eigensectorlist(list, list%nsector - 1)
      jsector = 0
      DO isector = 1, tmp%nsector
         IF (isector /= truerank) THEN
            jsector = jsector + 1
            CALL copy_eigensector(list%es(jsector), tmp%es(isector))
         ENDIF
      ENDDO
      CALL delete_eigensectorlist(tmp)

   end subroutine

   subroutine add_eigensector(es, list)

      implicit none

      TYPE(eigensectorlist_type), INTENT(INOUT) :: list
      TYPE(eigensector_type), INTENT(IN)        :: es
      INTEGER                    :: isector
      TYPE(eigensectorlist_type) :: tmp

      CALL copy_eigensectorlist(tmp, list)
      CALL new_eigensectorlist(list, list%nsector + 1)
      DO isector = 1, tmp%nsector
         CALL copy_eigensector(list%es(isector), tmp%es(isector))
      ENDDO
      CALL copy_eigensector(list%es(list%nsector), es)
      CALL delete_eigensectorlist(tmp)

   end subroutine

   subroutine new_eigensector_from_scratch(es, sector_in)

      ! NEW EMPTY EIGENSECTOR

      use sector_class, only: copy_sector, sector_type

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: es
      TYPE(sector_type), INTENT(IN)         :: sector_in

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: entering eigen_sector_class_new_eigensector_from_&
           &scratch"
#endif DEBUG

      CALL delete_eigensector(es)
      CALL copy_sector(es%sector, sector_in)

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: leaving eigen_sector_class_new_eigensector_from_&
           &scratch"
#endif DEBUG

   end subroutine

   subroutine copy_eigensector(esout, esin)

      use eigen_class, only: copy_eigenlist
      use sector_class, only: copy_sector

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: esout
      TYPE(eigensector_type), INTENT(IN)    :: esin

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: entering eigen_sector_class_copy_eigensector"
#endif DEBUG

      CALL delete_eigensector(esout)
      CALL copy_sector(esout%sector, esin%sector)
      CALL copy_eigenlist(esout%lowest, esin%lowest)

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: leaving eigen_sector_class_copy_eigensector"
#endif DEBUG

   end subroutine

   subroutine copy_eigensectorlist(listout, listin)

      implicit none

      TYPE(eigensectorlist_type), INTENT(INOUT) :: listout
      TYPE(eigensectorlist_type), INTENT(IN)    :: listin
      INTEGER :: isector

      CALL new_eigensectorlist(listout, listin%nsector)
      if (listin%nsector == 0) return
      DO isector = 1, listin%nsector
         CALL copy_eigensector(listout%es(isector), listin%es(isector))
      ENDDO

   end subroutine

   subroutine delete_eigensector(es)

      use eigen_class, only: delete_eigenlist
      use sector_class, only: delete_sector

      implicit none

      TYPE(eigensector_type), INTENT(INOUT) :: es

      CALL delete_sector(es%sector)
      CALL delete_eigenlist(es%lowest)
   end subroutine

   subroutine delete_eigensectorlist(list)

      use genvar, only: log_unit, messages3

      implicit none

      TYPE(eigensectorlist_type), INTENT(INOUT) :: list
      INTEGER :: isector

      DO isector = 1, list%nsector
         CALL delete_eigensector(list%es(isector))
      ENDDO
      if (messages3) write (log_unit, *) 'delete list'
      IF (ASSOCIATED(list%es)) DEALLOCATE (list%es)
      list%nsector = 0
      list%es => NULL()
      if (messages3) write (log_unit, *) '---done---'

   end subroutine

   subroutine filter_eigensector(list, window)

      ! FILTER OUT ALL EIGENSTATES OUTSIDE window IN ALL EIGENSECTORS OF list

      use common_def, only: utils_assert
      use eigen_class, only: filter_eigen
      use genvar, only: log_unit, messages4
      use stringmanip, only: toString

      implicit none

      TYPE(eigensectorlist_type), INTENT(IN) :: list
      REAL(DBL), INTENT(IN)                  :: window(2)
      INTEGER :: isector, nn

      isector = 1
      nn = size(list%es(:))

      DO
         if (messages4) then
            write (log_unit, *) '.....start to filter sector : ', isector, &
               window
            write (log_unit, *) '.....ISECTOR, NN, SIZE OF LIST : ', isector, &
               nn, size(list%es(:)), list%nsector
         endif
         if (isector > nn) exit
         if (list%es(isector)%lowest%neigen > 0) CALL &
            filter_eigen(list%es(isector)%lowest, window)
         IF (list%es(isector)%lowest%neigen == 0) then
            if (messages4) write (log_unit, *) '----remove empty sector----', &
               isector
            CALL remove_eigensector(list%es(isector)%sector, list)
            if (messages4) write (log_unit, *) 'done, now remains .. sectors : &
                 &', nn - 1
            nn = nn - 1
            call utils_assert(nn == list%nsector, 'Error in &
                 &filter_eigensector: nn = '//trim(adjustl(toString(nn))) &
                 //' is not equal to list%nsector = ' &
                 //trim(adjustl(toString(list%nsector))))
         else
            isector = isector + 1
         endif
      ENDDO

   end subroutine

   subroutine nambu_energy_shift(list, back)

      use genvar, only: log_unit, strongstop
      use globalvar_ed_solver, only: energy_global_shift, energy_global_shift2

      implicit none

      TYPE(eigensectorlist_type) :: list
      INTEGER                    :: isector, nn, k, kk
      logical, optional          :: back

      nn = size(list%es(:))

      DO isector = 1, nn

         if (ASSOCIATED(list%es(isector)%lowest%eigen)) then

            kk = size(list%es(isector)%lowest%eigen(:))

            do k = 1, kk
               if (.not. present(back)) then
                  list%es(isector)%lowest%eigen(k)%val = &
                     list%es(isector)%lowest%eigen(k)%val + &
                     energy_global_shift + energy_global_shift2
               else
                  list%es(isector)%lowest%eigen(k)%val = &
                     list%es(isector)%lowest%eigen(k)%val - &
                     energy_global_shift - energy_global_shift2
               endif
            enddo

         else

            write (*, *) 'SECTOR IN NAMBU NOT ASSOCIATED'
            write (*, *) 'SECTOR : ', isector, '  /  ', nn

            if (strongstop) stop 'critical'

         endif

      ENDDO

      write (log_unit, *) ' --- the ground state energy was shifted by --- '
      write (log_unit, *) ' ---           shift due to mu              --- '
      write (log_unit, *) energy_global_shift
      write (log_unit, *) ' ---       shift due to epsilons bath       --- '
      write (log_unit, *) energy_global_shift2

   end subroutine

   function GSenergy(list)

      ! FIND THE LOWEST EIGENVALUE IN list

      use eigen_class, only: min_eigen
      use genvar, only: huge_
      implicit none

      TYPE(eigensectorlist_type), INTENT(IN) :: list
      REAL(DBL) :: GSenergy
      INTEGER   :: isector
      REAL(DBL) :: E0

      GSenergy = huge_
      DO isector = 1, list%nsector
         E0 = min_eigen(list%es(isector)%lowest) ! GS ENERGY OF THE LOCAL SECTOR
         IF (E0 < GSenergy) GSenergy = E0
      ENDDO

   end function

   logical function not_commensurate_sector(list, isector)

      use genvar, only: size2

      implicit none

      TYPE(eigensectorlist_type), INTENT(IN) :: list
      INTEGER :: isector

      if (.not. associated(list%es(isector)%sector%updo)) then
         not_commensurate_sector = .true.
         return
      endif

      not_commensurate_sector = ANY(list%es(isector)%sector%updo%up%chunk /= &
                                    list%es(isector)%sector%updo%up%chunk(1)) .or. &
         ANY(list%es(isector)%sector%updo%down%chunk /= &
             list%es(isector)%sector%updo%down%chunk(1)) .or. &
         mod(list%es(isector)%sector%updo%up%chunk(1), size2) /= 0
   end function

   logical function not_commensurate_sector_(sector)

      use genvar, only: size2
      use sector_class, only: sector_type

      implicit none

      TYPE(sector_type), INTENT(IN) :: sector

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: entering eigen_sector_class_not_commensurate_sector_"
#endif

      if (.not. associated(sector%updo)) then
         not_commensurate_sector_ = .true.
         return
      endif

      not_commensurate_sector_ = ANY(sector%updo%up%chunk /= &
                                     sector%updo%up%chunk(1)) .or. ANY(sector%updo%down%chunk /= &
                                                               sector%updo%down%chunk(1)) .or. mod(sector%updo%up%chunk(1), size2) &
         /= 0

#ifdef DEBUG
      write (log_unit, '(a)') "DEBUG: leaving eigen_sector_class_not_commensurate_sector_"
#endif

   end function

   function partition(beta, list) RESULT(Zpart)

      use eigen_class, only: sum_boltz

      implicit none

      REAL(DBL), INTENT(IN)                  :: beta
      TYPE(eigensectorlist_type), INTENT(IN) :: list
      REAL(DBL) :: Zpart
      INTEGER   :: isector
      REAL(DBL) :: E0

      E0 = GSenergy(list)
      Zpart = 0.0_DP
      DO isector = 1, list%nsector
         Zpart = Zpart + sum_boltz(list%es(isector)%lowest, beta, E0)
      ENDDO

   end function

   function average_energy(beta, list) RESULT(EE)

      use eigen_class, only: sum_boltz, sum_boltz_times_e

      implicit none

      REAL(DBL), INTENT(IN)                  :: beta
      TYPE(eigensectorlist_type), INTENT(IN) :: list
      REAL(DBL) :: Zpart, EE
      INTEGER   :: isector
      REAL(DBL) :: E0

      E0 = GSenergy(list)
      Zpart = 0.0_DP
      EE = 0.0_DP
      DO isector = 1, list%nsector
         Zpart = Zpart + sum_boltz(list%es(isector)%lowest, beta, E0)
         EE = EE + sum_boltz_times_E(list%es(isector)%lowest, beta, E0)
      ENDDO
      EE = EE/Zpart

   end function

   function nlowest(list)

      implicit none

      TYPE(eigensectorlist_type), INTENT(IN) :: list
      INTEGER :: nlowest
      INTEGER :: isector

      nlowest = 0
      DO isector = 1, list%nsector
         nlowest = nlowest + list%es(isector)%lowest%neigen
      ENDDO
   end function

   subroutine print_eigensectorlist(list, UNIT)

      use common_def, only: c2s, dump_message, i2c
      use sector_class, only: dimen_func, title_sector
      use eigen_class, only: print_eigenlist

      implicit none

      TYPE(eigensectorlist_type), INTENT(IN) :: list
      INTEGER, OPTIONAL, INTENT(IN)          :: UNIT
      INTEGER :: isector

      DO isector = 1, list%nsector
         CALL dump_message(UNIT=UNIT, TEXT="# SECTOR "// &
                           TRIM(title_sector(list%es(isector)%sector))//" [dim = "// &
                           c2s(i2c(dimen_func(list%es(isector)%sector)))//"]")
         CALL print_eigenlist(list%es(isector)%lowest, UNIT=UNIT)
      ENDDO
   end subroutine

   subroutine write_raw_eigensectorlist(list, FILEOUT)

      use common_def, only: close_safe, open_safe
      use eigen_class, only: write_raw_eigenlist
      use string5, only: get_unit
      use sector_class, only: write_raw_sector

      implicit none

      TYPE(eigensectorlist_type), INTENT(IN) :: list
      CHARACTER(LEN=*), INTENT(IN)         :: FILEOUT
      INTEGER :: isector
      INTEGER :: UNIT

      CALL open_safe(UNIT, FILEOUT, 'UNKNOWN', 'WRITE', get_unit=.true.)
      DO isector = 1, list%nsector
         WRITE (UNIT, *) .true. ! there is a new item
         CALL write_raw_sector(list%es(isector)%sector, UNIT)
         CALL write_raw_eigenlist(list%es(isector)%lowest, UNIT)
      ENDDO
      WRITE (UNIT, *) .false. ! last item
      CALL close_safe(UNIT)
   end subroutine

   subroutine read_raw_eigensectorlist(list, FILEIN)

      use common_def, only: close_safe, open_safe
      use sector_class, only: read_raw_sector
      use eigen_class, only: read_raw_eigenlist

      implicit none

      TYPE(eigensectorlist_type), INTENT(INOUT) :: list
      CHARACTER(LEN=*), INTENT(IN)            :: FILEIN
      TYPE(eigensector_type) :: es
      INTEGER                :: UNIT
      LOGICAL                :: new_item

      CALL open_safe(UNIT, FILEIN, 'UNKNOWN', 'READ', get_unit=.true.)

      DO
         READ (UNIT, *) new_item
         IF (.NOT. new_item) EXIT
         CALL read_raw_sector(es%sector, UNIT)
         CALL read_raw_eigenlist(es%lowest, UNIT)
         CALL add_eigensector(es, list)
      ENDDO
      CALL close_safe(UNIT)

   end subroutine

   function rank_sector_in_list(sector, list)

      use sector_class, only: equal_sector, sector_type

      implicit none

      TYPE(sector_type), INTENT(IN)          :: sector
      TYPE(eigensectorlist_type), INTENT(IN) :: list
      INTEGER :: rank_sector_in_list
      INTEGER :: isector

      rank_sector_in_list = 0
      DO isector = 1, list%nsector
         IF (equal_sector(sector, list%es(isector)%sector)) THEN
            rank_sector_in_list = isector
            EXIT
         ENDIF
      ENDDO
   end function

end module
