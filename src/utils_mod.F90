! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=============================================================================!
!                                Utils module                                 !
!=============================================================================!
! The module for miscellaneous utilities in ONETEP. Routines in this module   !
! can only depend on the constants, comms and io modules.                     !
!-----------------------------------------------------------------------------!
!                                                                             !
!                 The ONETEP code is written and maintained by                !
!          Chris-Kriton Skylaris, Arash A. Mostofi and Peter D. Haynes        !
!                                                                             !
!-----------------------------------------------------------------------------!
! This file written by Peter Haynes, 22 March 2005.                           !
!=============================================================================!

module utils

  implicit none

  private

#ifdef DEBUG_ARRAYS
  type array_info
     character(len=36) :: name
     character(len=72) :: routine
     integer :: allocated
     integer :: max_allocated
  end type array_info

  integer, parameter :: max_arrays = 2000
  integer :: num_arrays
  type(array_info) :: arrays(max_arrays)
#endif

  public :: utils_init_array_checker
  public :: utils_array_checker_report

  interface utils_heapsort
     module procedure utils_heapsort_integer
     module procedure utils_heapsort_real
  end interface

  interface utils_assert
     module procedure utils_assert_with_string
     module procedure utils_assert_with_integers
     module procedure utils_assert_with_real
  end interface

  interface utils_use_var
     module procedure utils_use_var_character
     module procedure utils_use_var_integer
     module procedure utils_use_var_real
  end interface

  public :: utils_alloc_check
  public :: utils_dealloc_check
  public :: utils_flush
  public :: utils_unit
  public :: utils_binary_copy
  public :: utils_heapsort
  public :: utils_open_unit_check
  public :: utils_close_unit_check
  public :: utils_read_check
  public :: utils_write_check
  public :: utils_assert
  public :: utils_abort
  public :: utils_trace_in
  public :: utils_trace_out
  public :: utils_sanity_check
  public :: utils_erf
  public :: utils_erfc
  public :: utils_custom_erfc
  public :: utils_dump_array3D_to_file
  public :: utils_isnan
  public :: utils_use_var

  integer, private, save :: nindent = 0

contains

  subroutine utils_init_array_checker()

#ifdef DEBUG_ARRAYS
    num_arrays = 0
    arrays(:)%name = ''
    arrays(:)%routine = ''
    arrays(:)%allocated = 0
    arrays(:)%max_allocated = 0
#endif

  end subroutine utils_init_array_checker

  subroutine utils_array_checker_report()

    use comms, only: pub_on_root
    use constants, only: stdout
    implicit none

#ifdef DEBUG_ARRAYS
    integer :: iarray

    if (pub_on_root) then
       do iarray=1,num_arrays
          if (arrays(iarray)%allocated>0) &
               write(stdout,'(i5,i5,i5,1x,a36,a72)') iarray, &
               arrays(iarray)%allocated,arrays(iarray)%max_allocated, &
               arrays(iarray)%name,arrays(iarray)%routine
       end do
    end if
#endif

  end subroutine utils_array_checker_report

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !============================================================================!
  ! This subroutine checks whether an array has been allocated successfully.   !
  ! On failure it writes an error message and aborts.                          !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   routine (in) : Name of the subroutine which called this routine          !
  !   array   (in) : Name of the array allocated                               !
  !   status  (in) : Status flag returned by allocate                          !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2005.                                       !
  !============================================================================!

  subroutine utils_alloc_check(routine,array,status)

    use comms, only: pub_on_root, comms_abort
    use constants, only: stderr, stdout
    implicit none

    ! Arguments
    character(len=*), intent(in) :: routine    ! Routine name
    character(len=*), intent(in) :: array      ! Array name
    integer, intent(in)          :: status     ! Status flag

    ! Locals
#ifdef DEBUG_ARRAYS
    integer :: iarray
#endif

    if (status /= 0) then
       write(stderr,'(3a/3a,i6)') 'Error in ',routine,':', &
            '  allocating ',array,' failed with code ',status
       call comms_abort
    end if

#ifdef DEBUG_ARRAYS
    do iarray=1,num_arrays
       ! ndmh: if current array name and parent has been used before
       if  ((array == arrays(iarray)%name) .and. &
            (routine == arrays(iarray)%routine)) then
          ! ndmh: increment allocated count for this array
          arrays(iarray)%allocated = arrays(iarray)%allocated + 1
          arrays(iarray)%max_allocated = max(arrays(iarray)%allocated, &
               arrays(iarray)%max_allocated)
          return
       end if
    end do

    ! ndmh: make sure you don't exceed allocated memory
    if (num_arrays == max_arrays) then
       if (pub_on_root) write(stdout,'(a,i4,a)') 'Error in utils_alloc_check: &
            &maximum number of arrays (',max_arrays, ') exceeded'
       call comms_abort
    end if

    num_arrays = num_arrays + 1
    arrays(num_arrays)%name = array
    arrays(num_arrays)%routine = routine
    arrays(num_arrays)%allocated = 1
    arrays(iarray)%max_allocated = max(arrays(iarray)%allocated, &
         arrays(iarray)%max_allocated)
#endif

  end subroutine utils_alloc_check

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !============================================================================!
  ! This subroutine checks whether an array has been deallocated successfully. !
  ! On failure it writes an error message and aborts.                          !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   routine (in) : Name of the subroutine which called this routine          !
  !   array   (in) : Name of the array deallocated                             !
  !   status  (in) : Status flag returned by deallocate                        !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2005.                                       !
  !============================================================================!

  subroutine utils_dealloc_check(routine,array,status)

    use comms, only: pub_on_root, comms_abort
    use constants, only: stderr, stdout
    implicit none

    ! Arguments
    character(len=*), intent(in) :: routine    ! Routine name
    character(len=*), intent(in) :: array      ! Array name
    integer, intent(in)          :: status     ! Status flag

    ! Locals
#ifdef DEBUG_ARRAYS
    integer :: iarray
    logical :: any_routine
    character(72) :: routine_no_de
#endif

    if (status /= 0) then
       write(stderr,'(3a/3a,i6)') 'Error in ',routine,':', &
            '  deallocating ',array,' failed with code ',status
       call comms_abort
    end if

#ifdef DEBUG_ARRAYS


    routine_no_de = ''
    iarray = index(routine,'dealloc')
    if (iarray>0) then
       routine_no_de(1:iarray) = routine(1:iarray)
       routine_no_de(iarray:(len_trim(routine)-2)) = &
            routine((iarray+2):len_trim(routine))
    else
       routine_no_de = routine
    end if
    any_routine = .false.
    iarray = index(routine,'_exit')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'_free')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'_destroy')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'_close')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'classical_pot_dealloc')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'basis_sphere_deallocate')
    if (iarray>0) any_routine = .true.
    iarray = index(routine,'sparse_segments_dealloc')
    if (iarray>0) any_routine = .true.

    do iarray=1,num_arrays
       ! ndmh: find array in list of allocated arrays
       if ((array == arrays(iarray)%name) .and. &
            ((routine_no_de == arrays(iarray)%routine).or. &
            any_routine)) then
          ! ndmh: decrement number of allocated versions
          if (arrays(iarray)%allocated > 0) then
             arrays(iarray)%allocated = arrays(iarray)%allocated - 1
             return
          end if
       end if
    end do

    if (pub_on_root) write(stdout,'(3a/3a/)') 'WARNING in &
         &utils_dealloc_check: matching tag "',trim(array),'" not found', &
         'when deallocating in routine "',trim(routine),'"'
#endif

  end subroutine utils_dealloc_check

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !============================================================================!
  ! Flushes output (stdout by default)                                         !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes 12 November 2004                                   !
  ! Modified by Jacek Dziedzic on 08/06/2010 to accept an optional parameter,  !
  ! which, if .true., abandons the MPI reduction at the end. This allows one to!
  ! call this routine on certain nodes only, from places where services_flush  !
  ! is not accessible due to circular dependencies.                            !
  ! Default behaviour remains unchanged.                                       !
  !============================================================================!

  subroutine utils_flush(unit, dont_reduce_error_flag)

    use comms, only : pub_on_root, comms_reduce
    use constants, only: stdout

    implicit none

    ! Arguments
    integer, optional, intent(in) :: unit
    logical, optional, intent(in) :: dont_reduce_error_flag

    ! Local variables
    integer :: unit_local     ! Local copy of argument
    integer :: ierr           ! Error flag
#ifdef SUN
    integer :: flush
#endif
#ifndef INTFLUSH
    external flush
#endif

    if (present(unit)) then
       unit_local = unit
    else
       unit_local = stdout
    end if

    ! Set error flag to success
    ierr = 0

#ifdef SUN
    ierr = flush(unit_local)
#else
#ifdef SGI
    ! SGI introduced a second argument starting from the v7.4
    ! of their runtime environment. This form is still acceptable
    ! by the older environments - the return value is ignored in those
    ! older cases, and the code exits on error.
    ! This second argument is also allowed under UNICOS as an
    ! optional argument.
    call flush(unit_local,ierr)
#else
!CW
#ifndef FLUSH_EXT_NAME
          call flush(unit_local)
#else
          call flush_(unit_local)
#endif
!END CW
#endif
#endif

    ! jd: Prevent reduction, if asked to, so that flushing on selected
    !     nodes becomes possible
    if (present(dont_reduce_error_flag)) then
       if (dont_reduce_error_flag) return
    end if

    ! Compare error flags
    call comms_reduce('SUM',ierr)
    if (ierr /= 0 .and. pub_on_root) write(stdout,'(a)') &
         'WARNING in utils_flush: flush failed on at least one processor'

  end subroutine utils_flush

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !============================================================================!
  ! Finds a free unit for I/O                                                  !
  !----------------------------------------------------------------------------!
  ! Amalgamation of restart_find_unit by Chris-Kriton Skylaris and esdf_unit   !
  ! by Chris Pickard.                                                          !
  !============================================================================!

  integer function utils_unit()

    use comms, only : comms_abort
    use constants, only: stderr

    ! Local variables
    integer :: trial_unit
    integer :: ierr
    logical :: ex
    logical :: op

    utils_unit = -1
    do trial_unit=1,99
       inquire(unit=trial_unit,exist=ex,opened=op,iostat=ierr)
       if (ierr /= 0) then
          write(stderr,'(a,i2,a,i6)') &
               'Error in utils_unit: inquiring about unit ',trial_unit, &
               ' failed with code ',ierr
          call comms_abort
       end if
       if (ex .and. (.not. op)) then
          utils_unit = trial_unit
          exit
       end if
    end do

    if (utils_unit == -1) then
       write(stderr,'(a)') 'Error in utils_unit: no I/O units available'
       call comms_abort
    end if

  end function utils_unit

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !============================================================================!
  ! Copies data from/to binary files in the extensible format used for the     !
  ! density kernel.                                                            !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, February 2007.                                    !
  !============================================================================!

  subroutine utils_binary_copy(iunit,ounit)

    use comms, only : comms_abort
    use constants, only: stderr, DP

    implicit none

    ! Arguments
    integer, intent(in) :: iunit              ! Input I/O unit
    integer, intent(in) :: ounit              ! Output I/O unit

    ! Local variables
    integer :: ierr                           ! Error flag
    integer :: ilen                           ! Record length
    character(len=12) :: tag                  ! Record tag
    integer, allocatable :: ibuf(:)           ! I/O buffer for integer data
    real(kind=DP), allocatable :: dbuf(:)     ! I/O buffer for real data
    complex(kind=DP), allocatable :: zbuf(:)  ! I/O buffer for complex data


    ! Indefinite loop until endfile record found
    do
       ! Read tag for record and check for endfile
       read(iunit) tag
       if (tag(1:7) == 'ENDFILE') then
          write(ounit) tag
          write(ounit) 1
          write(ounit) 0
          exit
       end if


       ! Check record type is valid
       if (tag(12:12) /= 'I' .and. tag(12:12) /= 'D' .and. &
            tag(12:12) /= 'Z') then
          write(stderr,'(5a)') 'Error in utils_binary_copy: &
               &malformed record tag'
          call comms_abort
       end if
       write(ounit) tag


       ! Read record length
       read(iunit) ilen
       if (ilen <= 0) then
          write(stderr,'(5a)') 'Error in utils_binary_copy: &
               &malformed record length'
          call comms_abort
       end if
       write(ounit) ilen


       ! Allocate workspace for this type as required
       if (tag(12:12) == 'I') then
          if (allocated(ibuf)) then
             if (size(ibuf) < ilen) then
                deallocate(ibuf,stat=ierr)
                call utils_dealloc_check('utils_binary_copy','ibuf',ierr)
             end if
          end if
          if (.not. allocated(ibuf)) then
             allocate(ibuf(ilen),stat=ierr)
             call utils_alloc_check('utils_binary_copy','ibuf',ierr)
          end if
          read(iunit) ibuf(1:ilen)
          write(ounit) ibuf(1:ilen)
          deallocate(ibuf,stat=ierr)
          call utils_dealloc_check('utils_binary_copy','ibuf',ierr)
       else if (tag(12:12) == 'D') then
          if (allocated(dbuf)) then
             if (size(dbuf) < ilen) then
                deallocate(dbuf,stat=ierr)
                call utils_dealloc_check('utils_binary_copy','dbuf',ierr)
             end if
          end if
          if (.not. allocated(dbuf)) then
             allocate(dbuf(ilen),stat=ierr)
             call utils_alloc_check('utils_binary_copy','dbuf',ierr)
          end if
          read(iunit) dbuf(1:ilen)
          write(ounit) dbuf(1:ilen)
          deallocate(dbuf,stat=ierr)
          call utils_dealloc_check('utils_binary_copy','dbuf',ierr)
       else if (tag(12:12) == 'Z') then
          if (allocated(zbuf)) then
             if (size(zbuf) < ilen) then
                deallocate(zbuf,stat=ierr)
                call utils_dealloc_check('utils_binary_copy','zbuf',ierr)
             end if
          end if
          if (.not. allocated(zbuf)) then
             allocate(zbuf(ilen),stat=ierr)
             call utils_alloc_check('utils_binary_copy','zbuf',ierr)
          end if
          read(iunit) zbuf(1:ilen)
          write(ounit) zbuf(1:ilen)
          deallocate(zbuf,stat=ierr)
          call utils_dealloc_check('utils_binary_copy','zbuf',ierr)

       end if

       ! End loop over file
    end do

    return

  end subroutine utils_binary_copy

  !--------------------------------------------------------------------------

  subroutine utils_heapsort_integer(n,key,idx)

    !=======================================================================!
    ! This subroutine creates an index based on a list of keys, which are   !
    ! sorted by the heap sort algorithm, which is guaranteed to complete    !
    ! after O(n log n) operations.                                          !
    !-----------------------------------------------------------------------!
    ! Arguments:             !
    ! n (input)          : Number of items in list to be indexed            !
    ! key (input/output) : Key according to which the items are sorted      !
    ! idx (output)       : The resulting index                              !
    !-----------------------------------------------------------------------!
    ! Written by Peter Haynes, 10/7/03                               !
    !=======================================================================!

    use constants, only: LONG

    implicit none

    ! Arguments
    integer, intent(in) :: n
    integer(kind=LONG), intent(inout) :: key(n)
    integer, intent(out) :: idx(n)

    ! Local variables
    integer :: i,j                 ! Loop variables over list
    integer :: ic,is               ! Counters for heap creation/selection
    integer :: temp_idx            ! Temporary copy of index
    integer(kind=LONG) :: temp_key ! Temporary copy of key

    ! If the list is empty, there's nothing to do!
    if (n == 0) return

    ! If there is only one item in the list, there's no sorting to be done!
    if (n == 1) then
       idx(1) = 1
       return
    end if

    ! Initialise the sorting index
    do i=1,n
       idx(i) = i
    end do

    ! Do the heap sort
    ic = n/2+1
    is = n
    do
       if (ic > 1) then ! in heap creation phase
          ic = ic - 1
          temp_key = key(ic) ; temp_idx = idx(ic)
       else ! in heap selection phase
          temp_key = key(is) ; temp_idx = idx(is)
          key(is) = key(1) ; idx(is) = idx(1)
          is = is - 1
          if (is == 1) then
             key(1) = temp_key ; idx(1) = temp_idx
             exit
          end if
       end if
       ! Sift down temporary copy to correct level in heap
       i = ic
       j = 2 * ic
       do while (j <= is)
          if (j < is) then
             if (key(j) < key(j+1)) j = j + 1
          end if
          if (temp_key < key(j)) then ! demote temporary copy
             key(i) = key(j) ; idx(i) = idx(j)
             i = j
             j = 2 * j
          else ! found level for temporary copy
             j = is + 1
          end if
       end do
       key(i) = temp_key ; idx(i) = temp_idx
    end do

  end subroutine utils_heapsort_integer

  !--------------------------------------------------------------------------

  subroutine utils_heapsort_real(n,key,idx)

    !=======================================================================!
    ! This subroutine creates an index based on a list of keys, which are   !
    ! sorted by the heap sort algorithm, which is guaranteed to complete    !
    ! after O(n log n) operations.                                          !
    !-----------------------------------------------------------------------!
    ! Arguments:             !
    ! n (input)          : Number of items in list to be indexed            !
    ! key (input/output) : Key according to which the items are sorted      !
    ! idx (output)       : The resulting index                              !
    !-----------------------------------------------------------------------!
    ! Written by Peter Haynes, 10/7/03                               !
    !=======================================================================!

    use constants, only: DP

    implicit none

    ! Arguments
    integer, intent(in) :: n
    real(kind=DP), intent(inout) :: key(n)
    integer, intent(out) :: idx(n)

    ! Local variables
    integer :: i,j                 ! Loop variables over list
    integer :: ic,is               ! Counters for heap creation/selection
    integer :: temp_idx            ! Temporary copy of index
    real(kind=DP) :: temp_key      ! Temporary copy of key

    ! If the list is empty, there's nothing to do!
    if (n == 0) return

    ! If there is only one item in the list, there's no sorting to be done!
    if (n == 1) then
       idx(1) = 1
       return
    end if

    ! Initialise the sorting index
    do i=1,n
       idx(i) = i
    end do

    ! Do the heap sort
    ic = n/2+1
    is = n
    do
       if (ic > 1) then ! in heap creation phase
          ic = ic - 1
          temp_key = key(ic) ; temp_idx = idx(ic)
       else ! in heap selection phase
          temp_key = key(is) ; temp_idx = idx(is)
          key(is) = key(1) ; idx(is) = idx(1)
          is = is - 1
          if (is == 1) then
             key(1) = temp_key ; idx(1) = temp_idx
             exit
          end if
       end if
       ! Sift down temporary copy to correct level in heap
       i = ic
       j = 2 * ic
       do while (j <= is)
          if (j < is) then
             if (key(j) < key(j+1)) j = j + 1
          end if
          if (temp_key < key(j)) then ! demote temporary copy
             key(i) = key(j) ; idx(i) = idx(j)
             i = j
             j = 2 * j
          else ! found level for temporary copy
             j = is + 1
          end if
       end do
       key(i) = temp_key ; idx(i) = temp_idx
    end do

  end subroutine utils_heapsort_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_open_unit_check(routine,file_name,ios_status)

    !==========================================================================!
    ! This subroutine checks whether an unit has been opened successfully.     !
    ! On failure it writes an error message and aborts.                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   routine   (in) : Name of the subroutine which called this routine      !
    !   file_name (in) : Name of the file to be opened                         !
    !   ios_status(in) : Status flag returned by open                          !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano, July 2009.                               !
    !==========================================================================!

    use comms, only: comms_abort
    use constants, only: stderr
    implicit none

    ! Arguments
    character(len=*), intent(in) :: routine    ! Routine name
    character(len=*), intent(in) :: file_name  ! File name
    integer, intent(in)          :: ios_status ! Status flag

    if (ios_status /= 0) then
       write(stderr,'(/3a/3a,i6)') 'Error in ',routine,':', &
            '  failed to open file "',file_name,'" with code ',ios_status
       call comms_abort
    end if

  end subroutine utils_open_unit_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_close_unit_check(routine,unit_name,ios_status)

    !==========================================================================!
    ! This subroutine checks whether an unit has been closed successfully.     !
    ! On failure it writes an error message and aborts.                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   routine   (in) : Name of the subroutine which called this routine      !
    !   unit_name (in) : Name of the unit to be closed                         !
    !   ios_status(in) : Status flag returned by open                          !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano, July 2009.                               !
    !==========================================================================!

    use comms, only: comms_abort
    use constants, only: stderr
    implicit none

    ! Arguments
    character(len=*), intent(in) :: routine    ! Routine name
    character(len=*), intent(in) :: unit_name  ! File name
    integer, intent(in)          :: ios_status ! Status flag

    if (ios_status /= 0) then
       write(stderr,'(/3a/3a,i6)') 'Error in ',routine,':', &
            '  failed to close unit "',unit_name,'" with code ',ios_status
       call comms_abort
    end if

  end subroutine utils_close_unit_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_read_check(routine,var_name,ios_status)

    !==========================================================================!
    ! This subroutine checks whether a value has been readed successfully.     !
    ! On failure it writes an error message and aborts.                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   routine   (in) : Name of the subroutine which called this routine      !
    !   var_name  (in) : Name of the variable to be read                       !
    !   ios_status(in) : Status flag returned by open                          !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano, July 2009.                               !
    !==========================================================================!

    use comms, only: comms_abort
    use constants, only: stderr
    implicit none

    ! Arguments
    character(len=*), intent(in) :: routine    ! Routine name
    character(len=*), intent(in) :: var_name   ! Variable name
    integer, intent(in)          :: ios_status ! Status flag

    if (ios_status /= 0) then
       write(stderr,'(/3a/3a,i6)') 'Error in ',routine,':', &
            '  failed to read variable "',var_name,'" with code ',ios_status
       call comms_abort
    end if

  end subroutine utils_read_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_write_check(routine,var_name,ios_status)

    !==========================================================================!
    ! This subroutine checks whether a value has been written successfully.    !
    ! On failure it writes an error message and aborts.                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   routine   (in) : Name of the subroutine which called this routine      !
    !   var_name  (in) : Name of the variable to be read                       !
    !   ios_status(in) : Status flag returned by open                          !
    !--------------------------------------------------------------------------!
    ! Written by Alvaro Ruiz Serrano, July 2009.                               !
    !==========================================================================!

    use comms, only: comms_abort
    use constants, only: stderr
    implicit none

    ! Arguments
    character(len=*), intent(in) :: routine    ! Routine name
    character(len=*), intent(in) :: var_name   ! Variable name
    integer, intent(in)          :: ios_status ! Status flag

    if (ios_status /= 0) then
       write(stderr,'(/3a/3a,i6)') 'Error in ',routine,':', &
            '  failed to write variable "',var_name,'" with code ',ios_status
       call comms_abort
    end if

  end subroutine utils_write_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function utils_custom_erfc(x)
    !=========================================================================!
    ! Calculate an accurage approximation to the complemenary error function, !
    !    erfc(x)=(2/sqrt(pi))*integral(x->infinity)[exp(-t^2)]dt.             !
    !    Based upon parameterization given in NSWC Mathematics Library.       !
    !    NB This is machine portable and not dependent on a system function.  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   x, intent=in, argument of erfc to evaluate.                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   none.                                                                 !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.1, 25/09/2001                               !
    !=========================================================================!

    use constants, only: DP

    implicit none

    real(kind=dp), intent(in) :: x
    real(kind=dp)             :: utils_custom_erfc

    !expansion parameters
    real(kind=dp), dimension(1:5) :: a=(/  0.771058495001320E-04_dp, -0.133733772997339E-02_dp,  &
         &  0.323076579225834E-01_dp,  0.479137145607681E-01_dp,  0.128379167095513E+00_dp /)
    real(kind=dp), dimension(1:3) :: b=(/  0.301048631703895E-02_dp,  0.538971687740286E-01_dp,  &
         &  0.375795757275549E+00_dp /)
    real(kind=dp), dimension(1:8) :: p=(/ -1.36864857382717E-07_dp,  5.64195517478974E-01_dp, &
         &  7.21175825088309E+00_dp,  4.31622272220567E+01_dp,  1.52989285046940E+02_dp,  &
         &  3.39320816734344E+02_dp,  4.51918953711873E+02_dp,  3.00459261020162E+02_dp /)
    real(kind=dp), dimension(1:8) :: q=(/  1.00000000000000E+00_dp,  1.27827273196294E+01_dp, &
         &  7.70001529352295E+01_dp, 2.77585444743988E+02_dp,  6.38980264465631E+02_dp, &
         &  9.31354094850610E+02_dp, 7.90950925327898E+02_dp,  3.00459260956983E+02_dp /)
    real(kind=dp), dimension(1:5) :: r=(/  2.10144126479064E+00_dp,  2.62370141675169E+01_dp, &
         &  2.13688200555087E+01_dp, 4.65807828718470E+00_dp,  2.82094791773523E-01_dp /)
    real(kind=dp), dimension(1:4) :: s=(/ 9.41537750555460E+01_dp, 1.87114811799590E+02_dp, &
         &  9.90191814623914E+01_dp, 1.80124575948747E+01_dp /)
    real(kind=dp)                 :: c = .564189583547756_dp

    ! Local variables
    real(kind=dp) :: ax, x2, t, top, bot

    ax=abs(x)
    x2=x*x
    if (ax<0.5_dp) then
       t=x2
       top=((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0_dp
       bot=((  b(1)*t + b(2))*t + b(3)) * t + 1.0_dp
       utils_custom_erfc=0.5_dp + (0.5_dp-x*(top/bot))

    else if (ax<4.0_dp) then
       top=((((((p(1) *ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax + p(6))*ax  &
            &  + p(7))*ax + p(8)
       bot=((((((q(1) *ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax + q(6))*ax  &
            &  + q(7))*ax + q(8)
       utils_custom_erfc=exp(-x2)*top/bot
       if (x<0.0_dp) utils_custom_erfc=2.0_dp - utils_custom_erfc

    else
       if (x<=-5.6_dp) then
          utils_custom_erfc=2.0_dp    !large negative x limit

       else if (x>100.0_dp) then
          utils_custom_erfc=0.0_dp    !large positive x limit

       else
          t=1.0_dp/x2
          top=(((r(1)*t + r(2))*t + r(3))*t + r(4)) * t + r(5)
          bot=(((s(1)*t + s(2))*t + s(3))*t + s(4)) * t + 1.0_dp
          utils_custom_erfc=exp(-x2)*(c-t*top/bot)/ax
          if (x<0.0_dp) utils_custom_erfc=2.0_dp - utils_custom_erfc

       end if
    end if

    return
  end function utils_custom_erfc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function utils_erf(x)
    !==========================================================================!
    ! Returns erf(x), calculated from derf(x) which is either vendor-provided, !
    ! or linked from libc. Alternatively, iff CUSTOM_ERF is defined, it uses   !
    ! utils_custom_erfc, defined above, to compute erf(x).                     !
    !--------------------------------------------------------------------------!
    ! Both versions are accurate to at least 3D-15, utils_custom_erfc() is     !
    ! several times slower.                                                    !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, November 2010.                                !
    !==========================================================================!

    use constants, only: DP

    implicit none

    real(kind=DP) :: utils_erf

#ifdef __PGI
    real(kind=DP) :: derf ! PGI compiler requires this to be declared
#endif

    ! Arguments
    real(kind=DP), intent(in) :: x

    ! -----------------------------------------------------------------------

#ifndef CUSTOM_ERF
    utils_erf = derf(x)
#else
    utils_erf = 1.0_DP - utils_custom_erfc(x)
#endif

  end function utils_erf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function utils_erfc(x)
    !==========================================================================!
    ! Returns erfc(x), calculated from derfc(x), which is either               !
    ! vendor-provided, or linked from libc. Alternatively, iff CUSTOM_ERF is   !
    ! defined, it uses utils_custom_erfc, defined above.                       !
    !--------------------------------------------------------------------------!
    ! Both versions are accurate to at least 3D-15, utils_custom_erfc() is     !
    ! several times slower.                                                    !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, November 2010.                                !
    !==========================================================================!

    use constants, only: DP

    implicit none

    real(kind=DP) :: utils_erfc

#ifdef __PGI
    real(kind=DP) :: derfc ! PGI compiler requires this to be declared
#endif

    ! Arguments
    real(kind=DP), intent(in) :: x

    ! -----------------------------------------------------------------------

#ifndef CUSTOM_ERF
    utils_erfc = derfc(x)
#else
    utils_erfc = utils_custom_erfc(x)
#endif

  end function utils_erfc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_abort(error_message)

    !==========================================================================!
    ! Prints out 'error_message' to stderr and stdout, flushes (both on the    !
    ! root node),  then calls comms_abort.                                     !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   error_message (input): The error message to print before aborting.     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2010.                                   !
    !==========================================================================!

    use constants, only: stderr, stdout
    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    character(len=*), intent(in) :: error_message

    ! -----------------------------------------------------------------------

    if (pub_on_root) then
       write(stdout,'(a)') error_message
       write(stderr,'(a)') error_message
       call utils_flush(stderr,.true.)
    end if

    call comms_abort

  end subroutine utils_abort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_assert_with_integers(condition,error_message, &
       value_1_to_print,value_2_to_print,value_3_to_print,value_4_to_print)

    !==========================================================================!
    ! Asserts that 'condition' is .true., if not -- prints out 'error_message' !
    ! to stderr, followed by 'value_to_print', flushes and calls comms_abort.  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   condition (input): The condition to check.                             !
    !   error_message (input): The error message to print if assertion fails.  !
    !   value_to_print (input): The value to print if assertion fails.         !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Avoid using conditions with side effects (such as function calls).     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2010.                                   !
    !==========================================================================!

    ! Arguments
    logical, intent(in) :: condition
    character(len=*), intent(in) :: error_message
    integer, intent(in) :: value_1_to_print
    integer, intent(in), optional :: value_2_to_print
    integer, intent(in), optional :: value_3_to_print
    integer, intent(in), optional :: value_4_to_print

    ! Local variables
    character(len=256) :: value_1_as_string, value_2_as_string, &
         value_3_as_string, value_4_as_string

    ! -----------------------------------------------------------------------

    if(condition) then
       return ! jd: Incur no extra overhead
    else
       write(value_1_as_string,'(i0)') value_1_to_print
       value_2_as_string=''
       value_3_as_string=''
       value_4_as_string=''
       if(present(value_2_to_print)) &
            write(value_2_as_string,'(i0)') value_2_to_print
       if(present(value_3_to_print)) &
            write(value_3_as_string,'(i0)') value_3_to_print
       if(present(value_4_to_print)) &
            write(value_4_as_string,'(i0)') value_4_to_print
       call utils_assert_with_string(condition,error_message, &
            trim(value_1_as_string)//' '//&
            trim(value_2_as_string)//' '//&
            trim(value_3_as_string)//' '//&
            trim(value_4_as_string))
    end if

  end subroutine utils_assert_with_integers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_assert_with_real(condition,error_message,value_to_print)

    !==========================================================================!
    ! Asserts that 'condition' is .true., if not -- prints out 'error_message' !
    ! to stderr, followed by 'value_to_print', flushes and calls comms_abort.  !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   condition (input): The condition to check.                             !
    !   error_message (input): The error message to print if assertion fails.  !
    !   value_to_print (input): The value to print if assertion fails.         !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Avoid using conditions with side effects (such as function calls).     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2010.                                   !
    !==========================================================================!

    use constants, only: DP

    ! Arguments
    logical, intent(in) :: condition
    character(len=*), intent(in) :: error_message
    real(kind=DP), intent(in) :: value_to_print

    ! Local variables
    character(len=256) :: value_as_string

    ! -----------------------------------------------------------------------

    if(condition) then
       return ! jd: Incur no extra overhead
    else
       write(value_as_string,'(e20.12)') value_to_print
       call utils_assert_with_string(condition,error_message, &
            trim(value_as_string))
    end if

  end subroutine utils_assert_with_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_assert_with_string(condition,error_message,optional_string)

    !==========================================================================!
    ! Asserts that 'condition' is .true., if not -- prints out 'error_message' !
    ! to stderr, optionally followed by 'optional_string', flushes and         !
    ! calls comms_abort.                                                       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   condition (input): The condition to check.                             !
    !   error_message (input): The error message to print if assertion fails.  !
    !   optional_string (input): The optional string to print if assertion     !
    !                            fails.                                        !
    !--------------------------------------------------------------------------!
    ! Caveats:                                                                 !
    !   Avoid using conditions with side effects (such as function calls).     !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2010.                                   !
    !==========================================================================!

    use constants, only: stderr
    use comms, only: comms_abort, pub_my_node_id

    implicit none

    ! Arguments
    logical, intent(in) :: condition
    character(len=*), intent(in) :: error_message
    character(len=*), intent(in), optional :: optional_string

    ! -----------------------------------------------------------------------

    if(condition) return

    write(stderr,'(/a)') &
         '#####################################################################'
    write(stderr,'(a,i0)',advance='no') 'Assertion failed (on node #',pub_my_node_id
    write(stderr,'(a)') '): '
    write(stderr,'(a)') error_message
    if(present(optional_string)) write(stderr,'(a)') optional_string
    write(stderr,'(a)') &
         '#####################################################################'

    call utils_flush(stderr,.true.)

    call comms_abort

  end subroutine utils_assert_with_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_trace_in(msg)
    !=======================================================================!
    ! Enters a node in the call graph, outputs the caller name to stdout.   !
    !-----------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 06/2010.                                 !
    !=======================================================================!

    use comms, only: pub_my_node_id
    use constants, only: stdout

    implicit none

    character(len=*), intent(in) :: msg
    integer :: i

    ! -----------------------------------------------------------------------

#ifdef TRACE
    write(stdout,*)
    do i=1,nindent
       write(stdout,'(a)',advance='no') ' '
    end do

#ifdef MPI
    write(stdout,'(3a,i0,a)',advance='no') '-->> ', msg, ' (node #', &
         pub_my_node_id,')'
#else
    write(stdout,*) '-->> ', msg
#endif

    nindent = nindent+4
    call utils_flush(stdout,.true.)

#else
    i=iachar(msg(1:1)) ! jd: Pointless, but kills warning about i and msg unused
#endif

  end subroutine utils_trace_in

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_trace_out(msg)
    !=======================================================================!
    ! Exits a node in the call graph, outputs the caller name to stdout.    !
    !-----------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 06/2010.                                 !
    !=======================================================================!

    use constants, only: stdout
    use comms, only: pub_my_node_id

    implicit none

    character(len=*), intent(in) :: msg
    integer :: i

    ! -----------------------------------------------------------------------

#ifdef TRACE
    nindent = nindent-4
    write(stdout,*)
    do i=1,nindent
       write(stdout,'(a)',advance='no') ' '
    end do

#ifdef MPI
    write(stdout,'(3a,i0,a)',advance='no') '<<-- ', msg, ' (node #', &
         pub_my_node_id,')'
#else
    write(stdout,*) '<<-- ', msg
#endif

    call utils_flush(stdout,.true.)

#else
    i=iachar(msg(1:1)) ! jd: Pointless, but kills warning about i and msg unused
#endif

  end subroutine utils_trace_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function utils_isnan(x)
    !=========================================================================!
    ! Returns .true. iff the double precision argument x is NaN.              !
    ! Allows abstracting away from the platform-specific mess.                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   x (input) : the double precision value to check against NaN.          !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 01/2011.                                   !
    !=========================================================================!

    ! jd: This implementation was found to work with the following compilers:
    !     - Intel's ifort 10.1, 11.1
    !     - GNU fortran 4.3, 4.4
    !     - PGI fortran 7.1, 10.9, 11.4
    !     - Sun Pro fortran 95 8.5
    !
    ! jd: It does not work (fails to compile) with
    !     - GNU fortran 4.1.2
    !
    ! If you find this to be the case, add "-DNO_ISNAN" to your compile options.
    ! This will cause isnan to always return false.

    use constants, only: DP

#if defined(NO_ISNAN)
    real(kind=DP), intent(in) :: x

    utils_isnan = .false.
    return
#else
#if defined(__SUNPRO_F90) || defined(__SUNPRO_F95)
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: ieee_features
#endif

    real(kind=DP), intent(in) :: x

    ! jd: PGI uses an externally linked isnand function
#if defined(__PGI)
    logical isnand
    utils_isnan = isnand(x)
#else
    ! jd: SUN does it the modern way, with ieee_is_nan,
    !     requiring ieee_* intrinsics
#if defined(__SUNPRO_F90) || defined(__SUNPRO_F95)
    utils_isnan = ieee_is_nan(x)
#else
    ! jd: other compilers use C's isnan
    utils_isnan = isnan(x)
#endif
#endif
#endif

  end function utils_isnan

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_sanity_check(arr,arrname,l1,l2,l3)
    !=========================================================================!
    ! Scans a 3D array of double precision numbers for NaNs, Inf, -Inf and    !
    ! ridiculously large values.                                              !
    ! If any are found, the location and value of the first one is printed and!
    ! utils_abort is called.                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   arr (input)     : The array to check                                  !
    !   arrname (input) : The name of the array, written out if check fails.  !
    !   l1,l2,l3 (input, optional): Lower bounds of arr, only needed if the   !
    !                               array is not dimensioned from 1.          !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic in 06/2010.                                   !
    !=========================================================================!

    use comms, only: pub_my_node_id
    use constants, only: stderr, DP

    implicit none

    ! jd: Input arguments
    real(kind=DP), intent(in)     :: arr(:,:,:)
    character(len=*), intent(in)  :: arrname
    integer, intent(in), optional :: l1,l2,l3

#ifdef DEBUG
    integer :: i,j,k       ! jd: Indices
    integer :: ii,jj,kk    ! jd: Indices
    real(kind=DP) :: value ! jd: arr(i,j,k)
    logical :: failed      ! jd: .true. iff arr(i,j,k) failed sanity check
    integer :: myrank      ! jd: Alias for pub_my_node_id

    !------------------------------------------------------------------------


    myrank = pub_my_node_id ! jd: Local copy is accessible from debugger
    do k=lbound(arr,3), ubound(arr,3)
       do j=lbound(arr,2), ubound(arr,2)
          do i=lbound(arr,1), ubound(arr,1)

             value = arr(i,j,k)
             failed = .false.

             ! jd: Detect infinities and ridiculously large values
!             if(value+1D10 == value) then
!                write(stderr,*) 'Infinity or very large value in ', arrname
!                failed = .true.
!             end if
             ! jd: Detect -infinities and ridiculously large negative values
!             if(value-1D10 == value) then
!                write(stderr,*) '-Inf or very large negative value in ', arrname
!                failed = .true.
!             end if
             if(utils_isnan(value)) then
                write(stderr,*) 'NaN detected in ', arrname
                failed = .true.
             end if

             if(failed) then
                write(stderr,*) 'Element at ', i,' ',j,' ',k,' is :',arr(i,j,k)
                if(present(l1) .or. present(l2) .or. present(l3)) then
                   ii=i
                   jj=j
                   kk=k
                   if(present(l1)) ii=ii+l1-1
                   if(present(l2)) jj=jj+l2-1
                   if(present(l3)) kk=kk+l3-1
                   write(stderr,*) 'NB: The array might not be dimensioned from&
                        & 1 -- actually the element is at ',ii,' ',jj,' ',kk
                end if
                write(stderr,*) 'node: ',myrank
                call utils_flush(stderr,.true.)
                call utils_abort('Sanity check failed.')
             end if

          end do
       end do
    end do
#endif

    ! jd: Suppress 'variable unused' warnings when DEBUG not defined
    if(present(l1)) call utils_use_var(l1)
    if(present(l2)) call utils_use_var(l2)
    if(present(l3)) call utils_use_var(l3)
    call utils_use_var(arrname)
    call utils_use_var(arr(lbound(arr,1),lbound(arr,2),lbound(arr,3)))

  end subroutine utils_sanity_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_dump_array3D_to_file(arr, filename, output_format)
    !==========================================================================!
    ! Writes the contents of a 3D local (non-distributed) array to a file      !
    ! for debugging purposes. The output filename is constructed by suffixing  !
    ! the passed filename with the node number and a global running counter.   !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   arr (input): the 3D array to write out.                                !
    !   filename (input): The basename of the file where arr will be written.  !
    !   output_format (input, optional): The format to use for output.         !
    !                                    If omitted, it defaults to (e20.12).  !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, November 2010.                                !
    !==========================================================================!

    use constants, only: DP
    use comms, only: pub_my_node_id

    implicit none

    ! jd: Input arguments
    real(kind=DP), intent(in)              :: arr(:,:,:)
    character(len=*), intent(in)           :: filename
    character(len=*), intent(in), optional :: output_format

    !------------------------------------------------------------------------

    ! jd: Local variables
    integer :: outunit
    integer, save :: file_counter = 0
    character(len=200) :: outfilename
    character(len=200) :: suffix1, suffix2
    integer :: i,j,k
    logical :: at_least_one

    ! jd: Construct the full filename
    write(suffix1,'(i0)') pub_my_node_id
    write(suffix2,'(i0)') file_counter
    outfilename=trim(trim(filename//'_node'//suffix1)//'_'//suffix2)

    ! jd: Open the file
    outunit = utils_unit()
    open(outunit, file=outfilename, err=10)

    ! jd: Output the contents of the array
    at_least_one = .false.
    do k=lbound(arr,3), ubound(arr,3)
       do j=lbound(arr,2), ubound(arr,2)
          do i=lbound(arr,1), ubound(arr,1)
             at_least_one = .true.
             if(present(output_format)) then
                write(outunit,output_format) arr(i,j,k)
             else
                write(outunit,'(e20.12)') arr(i,j,k)
             end if
          end do
       end do
    end do

    ! jd: Close the file, increment running counter
    close(outunit, err=20)
    file_counter = file_counter + 1

    return

    ! jd: I/O error handling
10  call utils_abort('Error during creation of file: '//outfilename)
20  call utils_abort('Error during closing of file: '//outfilename)

  end subroutine utils_dump_array3D_to_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_use_var_integer(x)
    !==========================================================================!
    ! Pretends to use its argument, so that compiler warnings about unused     !
    ! variables can be suppressed.                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   x (input): An integer value to be 'used'.                              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2011.                                   !
    !==========================================================================!

    integer, intent(in) :: x
    integer :: dummy
    dummy = x

  end subroutine utils_use_var_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_use_var_real(x)
    !==========================================================================!
    ! Pretends to use its argument, so that compiler warnings about unused     !
    ! variables can be suppressed.                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   x (input): An integer value to be 'used'.                              !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2011.                                   !
    !==========================================================================!
    use constants, only: DP

    real(kind=DP), intent(in) :: x
    real(kind=DP) :: dummy
    dummy = x

  end subroutine utils_use_var_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine utils_use_var_character(x)
    !==========================================================================!
    ! Pretends to use its argument, so that compiler warnings about unused     !
    ! variables can be suppressed.                                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   x (input): A character value to be 'used'.                             !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2011.                                   !
    !==========================================================================!

    character(len=1), intent(in) :: x
    character(len=1) :: dummy
    dummy = x

  end subroutine utils_use_var_character

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module utils
