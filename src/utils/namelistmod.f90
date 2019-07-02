module namelistmod

   use StringManip
   use string4, only: s_to_i4, chrctp
   use string5, only: s_eqi, get_unit, s_to_r8
   use common_def, only: get_address, put_address, open_safe
   use genvar

   IMPLICIT NONE

   private
   public :: look_for_command_line_argument
   public :: look_for_namelist_in_file
   public :: namelist_init
   public :: namelist_set
   public :: putel_in_namelist
   public :: putel_in_namelist_char
   public :: putel_in_namelist_i
   public :: putel_in_namelist_ivec
   public :: putel_in_namelist_l
   public :: putel_in_namelist_r

   TYPE namelist_member
      integer(8)             ::  address
      integer                ::  type, nsize
      character(200)         ::  a, a0
      logical                ::  l, l0
      integer                ::  i, i0
      real(kind=DP)                ::  r, r0
      integer(4), allocatable ::  ivec(:), ivec0(:)
      real(kind=DP), allocatable    ::  v(:), v0(:)
      complex(kind=DP), allocatable ::  vc(:), vc0(:)
      complex(kind=DP)             ::  c, c0
      character(400)         ::  comment, label
   END TYPE

!BUG September 6th 2012, corrected, indeed namelist_set is public, so have to be all types contained in namelist_set
!                        so following line was added
   public :: namelist_member

   TYPE namelist_set
      integer                           :: nel, cc
      TYPE(namelist_member), allocatable :: members(:)
      character(200)                    :: tree(100)
      character(200)                    :: namelist_name
      integer                           :: ntree = 0
   END TYPE

   ! INTERFACE get_info
   !   MODULE PROCEDURE get_infoi,get_infor,get_infoc,get_infol
   ! END INTERFACE

   INTERFACE putel_in_namelist
      MODULE PROCEDURE putel_in_namelist_i, putel_in_namelist_r, putel_in_namelist_c, &
                       & putel_in_namelist_l, putel_in_namelist_rr, putel_in_namelist_rvec, putel_in_namelist_char, &
                         putel_in_namelist_cvec, putel_in_namelist_ivec
   END INTERFACE

   logical, private :: wait
   logical, private :: menu_parsed = .false.

contains

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine build_namelist_tree(nm)
      TYPE(namelist_set) :: nm
      integer            :: i
      write (*, *) '====== start to build simulation tree in local directory ======'
      if (nm%ntree == 0) return
      do i = 1, nm%ntree
#ifndef NO_SYS_CALL
         if (rank == 0) call system("mkdir  "//TRIM(ADJUSTL(nm%tree(i))))
#endif
      enddo
      write (*, *) 'tree built....'
   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine look_for_namelist_in_file(nm, filename, nosync)
      TYPE(namelist_set) :: nm
      integer            :: i, jjj
      character(200)     :: value
      character(*)      :: filename
      logical, optional   :: nosync
      call get_unit(jjj)
      open (jjj, file=trim(adjustl(filename)))
      i = 0
      do
         i = i + 1
         read (jjj, '(a200)', end=10, err=11) value
11       continue
         if (rank == 0 .and. testing) write (*, *) 'in file : ', i, TRIM(ADJUSTL(value))
         call update_data_in_namelist(nm, TRIM(ADJUSTL(value)))
      enddo
10    close (jjj)

      if (.not. present(nosync)) then
      if (size2 > 1 .and. .not. no_mpi) then
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      endif
      write (*, *) 'file sparsed ... return - RANK / NOMPI / SIZE / FILENAME ', rank, no_mpi, size2, trim(adjustl(filename))
      endif

   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine look_for_command_line_argument(nm, nosync)
      TYPE(namelist_set) :: nm
      integer            :: i, num, len, status
      character(200)     :: value
      integer, parameter  :: nchar = 200
      logical, optional   :: nosync

      if (.not. present(nosync) .and. .not. no_mpi) then
         if (rank == 0) num = COMMAND_ARGUMENT_COUNT()
      else
         num = COMMAND_ARGUMENT_COUNT()
      endif

      if (rank == 0) then
         write (*, *) "==========================================="
         write (*, *) '  number of command line arguments : ', num
         write (*, *) '  Usage : '
         write (*, *) '         help       -  to print the full list of arguments'
         write (*, *) '         summary    -  to print a shortlist of argument definitions'
         write (*, *) '         build_tree -  build local tree of directories'
         write (*, *) '         wait       -  wait a few seconds to read this menu'
      endif

      if (.not. present(nosync)) then
      if (size2 > 1 .and. .not. no_mpi) then
         call MPI_BCAST(num, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      endif
      endif

      wait = .false.

      do i = 1, num
         if (present(nosync) .or. no_mpi) then
            write (*, *) 'WARNING : command line argument and no-sync or no_mpi'
            write (*, *) '          present nosync ? ', present(nosync)
            write (*, *) '          no_mpi         ? ', no_mpi
            write (*, *) '          command line arg are not always well read by childrens'
         endif
         if (.not. present(nosync) .and. .not. no_mpi) then
            if (rank == 0) call GET_COMMAND_ARGUMENT(i, value, len, status)
            if (rank == 0) write (*, *) 'argument : ', i, TRIM(ADJUSTL(value))
            if (size2 > 1 .and. .not. no_mpi) then
               call MPI_BCAST(value(1:nchar), nchar, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
               call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            endif
         else
            call GET_COMMAND_ARGUMENT(i, value, len, status)
            write (*, *) 'argument : ', i, TRIM(ADJUSTL(value))
         endif

         SELECT CASE (TRIM(ADJUSTL(value)))

         CASE ("help")
            wait = .true.
            call display_menu(nm)
            goto 11
         CASE ("summary")
            wait = .true.
            call display_menu_short(nm)
            goto 11
         CASE ("build_tree")
            if (.not. menu_parsed) call build_namelist_tree(nm)
         CASE ("wait")
            wait = .true.
         END SELECT
#ifndef NO_SYS_CALL
         if (wait) call system("sleep 3")
#endif
         write (*, *) 'update data in namelist '
         call update_data_in_namelist(nm, TRIM(ADJUSTL(value)))
11       continue

      enddo
      if (rank == 0) write (*, *) "==========================================="
      menu_parsed = .true.

      if (.not. present(nosync)) then
      if (size2 > 1 .and. .not. no_mpi) then
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      endif
      endif

      write (*, *) 'command line parsed ... return - RANK / NOMPI / SIZE ', rank, no_mpi, size2

   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine update_data_in_namelist(nm, label)
      TYPE(namelist_set) :: nm
      integer            :: i, j, length, ierror
      character(*)      :: label
      logical            :: lll
      integer            :: iii, jjj, jjj0
      real(kind=DP)            :: rrr
      complex(kind=DP)         :: ccc
      real(kind=DP)            :: vec(1000), ttt
      complex(kind=DP)         :: vecc(1000)
      integer(4)         :: ivec(1000), ss, iss, ok, ok2

      ok = 0; ok2 = 0

      do i = 1, len(label)
         if (label(i:i) == '=') then
            ok2 = 1
            exit
         endif
      enddo

      do j = 1, nm%cc

         if (s_eqi(label(1:i - 1), nm%members(j)%label)) then

            if (messages3) write (*, *) 'VARIABLE FOUND IN FILE : ', label(1:i - 1)

            if (nm%members(j)%type == 1) then
               lll = label(i + 1:len(label)) == ".true."
               call put_address(nm%members(j)%address, lll)
               nm%members(j)%l = lll
            elseif (nm%members(j)%type == 2) then
               call s_to_r8(label(i + 1:len(label)), rrr, ierror, length)
               iii = NINT(rrr)
               call put_address(nm%members(j)%address, iii)
               nm%members(j)%i = iii
            elseif (nm%members(j)%type == 3) then
               call s_to_r8(label(i + 1:len(label)), rrr, ierror, length)
               call put_address(nm%members(j)%address, rrr)
               nm%members(j)%r = rrr
            elseif (nm%members(j)%type == 4) then
               call chrctp(label(i + 1:len(label)), ccc, ierror, length)
               call put_address(nm%members(j)%address, ccc)
               nm%members(j)%c = ccc
            elseif (nm%members(j)%type == 5) then
               jjj = i; jjj0 = 0; iii = 1; vec = 0.
               do
                  jjj = jjj + 1
                  if (label(jjj:jjj) == '[' .or. label(jjj:jjj) == ',' .or. label(jjj:jjj) == ']') then
                     if (iii > nm%members(j)%nsize) &
                   & stop 'error put_name_in_namelist, size in namelist is not matching the size read to file/keyboard'
                     if (iii > 998) stop 'error put_name_in_namelist, size of temporary array too small'
                     if (label(jjj:jjj) /= '[') call s_to_r8(label(jjj0:jjj - 1), vec(iii), ierror, length)
                     if (label(jjj:jjj) == ',') iii = iii + 1
                     if (label(jjj:jjj) == ']') exit
                     jjj0 = jjj + 1
                  endif
               enddo
               if (iii /= nm%members(j)%nsize) then
                  write (*, *) 'found .. elements in vec : ', iii, vec(1:iii)
                  write (*, *) 'and should be.. elements : ', nm%members(j)%nsize
                  stop
               endif
               call put_address(nm%members(j)%address, vec(1:iii))
               nm%members(j)%v(1:iii) = vec(1:iii)
            elseif (nm%members(j)%type == 6) then
               ss = nm%members(j)%nsize
               do iss = 1, ss; nm%members(j)%a(iss:iss) = " "; enddo
               if (len(label) - i > ss) stop 'error update_data_in_namelist character : value too long for variable'
               nm%members(j)%a(1:len(label) - i) = label(i + 1:len(label))
               call put_address(nm%members(j)%address, nm%members(j)%a(1:ss))
            elseif (nm%members(j)%type == 7) then
               jjj = i; jjj0 = 0; iii = 1; vecc = 0.d0
               do
                  jjj = jjj + 1
                  if (label(jjj:jjj) == '[' .or. label(jjj:jjj) == ',' .or. label(jjj:jjj) == ']' .or. label(jjj:jjj) == 'I') then
                     if (iii > 998) stop 'error put_name_in_namelist, size of temporary array too small'
                     if (label(jjj:jjj) /= '[') then
                        call s_to_r8(label(jjj0:jjj - 1), ttt, ierror, length)
                        if (iii > nm%members(j)%nsize) then
                           write (*, *) 'iii : ', iii
                           write (*, *) 'member%nsize in namelist : ', nm%members(j)%nsize
                           stop 'error put_name_in_namelist, size in namelist is not matching the size read to file/keyboard'
                        endif
                        if (label(jjj:jjj) == 'I') then
                           vecc(iii) = vecc(iii) + ttt
                        else
                           vecc(iii) = vecc(iii) + imi*ttt
                        endif
                        if (label(jjj:jjj) == ',') iii = iii + 1
                     endif
                     if (label(jjj:jjj) == ']') exit
                     jjj0 = jjj + 1
                  endif
               enddo
               if (iii /= nm%members(j)%nsize) then
                  write (*, *) 'found ... elements in vec  : ', iii, vecc(1:iii)
                  write (*, *) 'and should be ... elements : ', nm%members(j)%nsize
                  stop
               endif
               call put_address(nm%members(j)%address, vecc(1:iii))
               nm%members(j)%vc(1:iii) = vecc(1:iii)
            elseif (nm%members(j)%type == 8) then
               jjj = i; jjj0 = 0; iii = 1; ivec = 0
               do
                  jjj = jjj + 1
                  if (label(jjj:jjj) == '[' .or. label(jjj:jjj) == ',' .or. label(jjj:jjj) == ']') then
                     if (iii > nm%members(j)%nsize) then
                        write (*, *) 'size indicated in namelist : ', nm%members(j)%nsize
                        write (*, *) ' iii  is now               : ', iii
                        stop 'error put_name_in_namelist, size in namelist is not matching the size read to file/keyboard'
                     endif
                     if (iii > 998) stop 'error put_name_in_namelist, size of temporary array too small'
                     if (label(jjj:jjj) /= '[') call s_to_i4(label(jjj0:jjj - 1), ivec(iii), ierror, length)
                     if (label(jjj:jjj) == ',') iii = iii + 1
                     if (label(jjj:jjj) == ']') exit
                     jjj0 = jjj + 1
                  endif
               enddo
               if (iii /= nm%members(j)%nsize) then
                  write (*, *) 'found .. elements in vec : ', iii, ivec(1:iii)
                  write (*, *) 'and should be.. elements : ', nm%members(j)%nsize
                  stop
               endif
               call put_address(nm%members(j)%address, ivec(1:iii))
               nm%members(j)%ivec(1:iii) = ivec(1:iii)
            endif

            ok = 1

         endif

      enddo

      if (ok == 0 .and. ok2 > 0) then
         if (testing) then
            write (*, *) 'WARNING : VARIABLE IN FILE DOES NOT CORRESPOND TO A NAMELIST VARIABLE -> ', label(1:i - 1)
            write (*, *) 'INSIDE NAMELIST CALLED : ', nm%namelist_name
         endif
         if (strongstop .or. testing) then
            write (*, *) ' list of all possible variables : '
            do j = 1, nm%cc; write (*, *) nm%members(j)%label; enddo
            !stop 'critical error'
         endif
      endif

      return
   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine display_menu(nm)
      TYPE(namelist_set) :: nm
      integer            :: ii, kk
      if (rank /= 0) return
      do ii = 1, nm%cc
         write (*, *) '==============================='
         write (*, *) 'VARIABLE      = ', TRIM(ADJUSTL(nm%members(ii)%label))
         write (*, *) 'MEM ADDRESS   = ', nm%members(ii)%address
         write (*, *) 'TYPE          = ', nm%members(ii)%type
         if (nm%members(ii)%type == 1) then
            write (*, *) 'VALUE         = ', nm%members(ii)%l
            write (*, *) 'DEFAULT VALUE = ', nm%members(ii)%l0
         elseif (nm%members(ii)%type == 2) then
            write (*, *) 'VALUE         = ', nm%members(ii)%i
            write (*, *) 'DEFAULT VALUE = ', nm%members(ii)%i0
         elseif (nm%members(ii)%type == 3) then
            write (*, *) 'VALUE         = ', nm%members(ii)%r
            write (*, *) 'DEFAULT VALUE = ', nm%members(ii)%r0
         elseif (nm%members(ii)%type == 4) then
            write (*, *) 'VALUE         = ', nm%members(ii)%c
            write (*, *) 'DEFAULT VALUE = ', nm%members(ii)%c0
         elseif (nm%members(ii)%type == 5) then
            write (*, '(a17,200f10.4)') 'VALUE         = ', (nm%members(ii)%v(kk), kk=1, nm%members(ii)%nsize)
            write (*, '(a17,200f10.4)') 'DEFAULT VALUE = ', (nm%members(ii)%v0(kk), kk=1, nm%members(ii)%nsize)
         elseif (nm%members(ii)%type == 6) then
            write (*, *) 'VALUE         = ', nm%members(ii)%a
            write (*, *) 'DEFAULT VALUE = ', nm%members(ii)%a0
         elseif (nm%members(ii)%type == 7) then
            write (*, '(a17,200f10.4)') 'VALUE         = ', (nm%members(ii)%vc(kk), kk=1, nm%members(ii)%nsize)
            write (*, '(a17,200f10.4)') 'DEFAULT VALUE = ', (nm%members(ii)%vc0(kk), kk=1, nm%members(ii)%nsize)
         elseif (nm%members(ii)%type == 8) then
            write (*, '(a17,200i5)') 'VALUE         = ', (nm%members(ii)%ivec(kk), kk=1, nm%members(ii)%nsize)
            write (*, '(a17,200i5)') 'DEFAULT VALUE = ', (nm%members(ii)%ivec0(kk), kk=1, nm%members(ii)%nsize)
         endif
         write (*, *) 'COMMENT       = ', TRIM(ADJUSTL(nm%members(ii)%comment))
         write (*, *) '==============================='
      enddo
#ifndef NO_SYS_CALL
      if (wait) call system("sleep 5")
#endif
   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine display_menu_short(nm)
      TYPE(namelist_set) :: nm
      integer            :: ii, kk, unit_
      integer, save       :: count = 0
      character(22)      :: filename

      if (rank /= 0) return

      count = count + 1
      filename = 'namelist_number_'//TRIM(ADJUSTL(toString(count)))
      CALL open_safe(unit_, TRIM(ADJUSTL(filename)), "UNKNOWN", "WRITE", get_unit=.true.)

      call dump_(6)
      call dump_(unit_)
#ifndef NO_SYS_CALL
      if (wait) call system("sleep 5")
#endif
      close (unit_)

   contains

      subroutine dump_(unit)
         integer :: unit
         do ii = 1, nm%cc
            write (unit, *)
            write (unit, *) '=================================================================='
            write (unit, '(a55)', ADVANCE='NO') nm%members(ii)%label
            if (nm%members(ii)%type == 1) then
               write (unit, '(2l4)', ADVANCE='NO') nm%members(ii)%l, nm%members(ii)%l0
            elseif (nm%members(ii)%type == 2) then
               write (unit, '(2i12)', ADVANCE="NO") nm%members(ii)%i, nm%members(ii)%i0
            elseif (nm%members(ii)%type == 3) then
               write (unit, '(2f10.4)', ADVANCE="NO") nm%members(ii)%r, nm%members(ii)%r0
            elseif (nm%members(ii)%type == 4) then
               write (unit, '(4f10.4)', ADVANCE="NO") nm%members(ii)%c, nm%members(ii)%c0
            elseif (nm%members(ii)%type == 5) then
               write (unit, '(200f10.4)', ADVANCE="NO") (nm%members(ii)%v(kk), kk=1, nm%members(ii)%nsize)
               write (unit, '(200f10.4)', ADVANCE="NO") (nm%members(ii)%v0(kk), kk=1, nm%members(ii)%nsize)
            elseif (nm%members(ii)%type == 6) then
               write (unit, '(2a20)', ADVANCE="NO") nm%members(ii)%a, nm%members(ii)%a0
            elseif (nm%members(ii)%type == 7) then
               write (unit, '(200f10.4)', ADVANCE="NO") (nm%members(ii)%vc(kk), kk=1, nm%members(ii)%nsize)
               write (unit, '(200f10.4)', ADVANCE="NO") (nm%members(ii)%vc0(kk), kk=1, nm%members(ii)%nsize)
            elseif (nm%members(ii)%type == 8) then
               write (unit, '(200i5)', ADVANCE="NO") (nm%members(ii)%ivec(kk), kk=1, nm%members(ii)%nsize)
               write (unit, '(200i5)', ADVANCE="NO") (nm%members(ii)%ivec0(kk), kk=1, nm%members(ii)%nsize)
            endif
            write (unit, '(a)', ADVANCE="NO") "    "//TRIM(ADJUSTL(nm%members(ii)%comment))
         enddo
         write (unit, *)
         write (unit, *) '=================================================================='
      end subroutine

   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!
!
! subroutine get_infoi(nm,i)
!  TYPE(namelist_set) :: nm
!  integer            :: i
!  integer            :: ii
!  integer(8)         :: j
!  j=get_address(i)
!  do ii=1,nm%cc
!   if(nm%members(ii)%address==j) exit
!  enddo
!  if(rank/=0) return
!  write(*,*) 'VARIABLE      = ' , nm%members(ii)%label
!  write(*,*) 'DEFAULT VALUE = ' , nm%members(ii)%i0
!  write(*,*) 'COMMENT       = ' , nm%members(ii)%comment
! end subroutine
!
!
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
!
!
! subroutine get_infor(nm,i)
!  TYPE(namelist_set) :: nm
!  real(kind=DP)            :: i
!  integer            :: ii
!  integer(8)         :: j
!  if(rank/=0) return
!  j=get_address(i)
!  do ii=1,nm%cc
!   if(nm%members(ii)%address==j) exit
!  enddo
!  write(*,*) 'VARIABLE      = ' , nm%members(ii)%label
!  write(*,*) 'DEFAULT VALUE = ' , nm%members(ii)%i0
!  write(*,*) 'COMMENT       = ' , nm%members(ii)%comment
! end subroutine
!
!
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
!
!
! subroutine get_infoc(nm,i)
!  TYPE(namelist_set) :: nm
!  complex(kind=DP)         :: i
!  integer            :: ii
!  integer(8)         :: j
!  if(rank/=0) return
!  j=get_address(i)
!  do ii=1,nm%cc
!   if(nm%members(ii)%address==j) exit
!  enddo
!  write(*,*) 'VARIABLE      = ' , nm%members(ii)%label
!  write(*,*) 'DEFAULT VALUE = ' , nm%members(ii)%i0
!  write(*,*) 'COMMENT       = ' , nm%members(ii)%comment
! end subroutine
!
!
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
!
!
! subroutine get_infol(nm,i)
!  TYPE(namelist_set) :: nm
!  logical            :: i
!  integer            :: ii
!  integer(8)         :: j
!  if(rank/=0) return
!  j=get_address(i)
!  do ii=1,nm%cc
!   if(nm%members(ii)%address==j) exit
!  enddo
!  write(*,*) 'VARIABLE      = ' , nm%members(ii)%label
!  write(*,*) 'DEFAULT VALUE = ' , nm%members(ii)%i0
!  write(*,*) 'COMMENT       = ' , nm%members(ii)%comment
! end subroutine
!
!
!
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine namelist_init(nm, nel, tree, name_of_namelist)
      TYPE(namelist_set)      :: nm
      integer                 :: nel, i
      character(*), optional   :: tree(:)
      character(*), optional   :: name_of_namelist

      nm%nel = nel
      if (allocated(nm%members)) deallocate (nm%members)
      allocate (nm%members(nel))
      nm%cc = 0
      if (present(tree)) then
         nm%ntree = size(tree)
         do i = 1, nm%ntree
            nm%tree(i) = tree(i)
         enddo
      else
         nm%ntree = 0
      endif
      if (present(name_of_namelist)) then
         nm%namelist_name = name_of_namelist
      else
         nm%namelist_name = 'default namelist'
      endif

   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine putel_in_namelist_i(nm, i, label, i0, comment)
      TYPE(namelist_set) :: nm
      integer            :: i, i0
      character(*)      :: comment, label
      if (nm%cc > nm%nel) stop 'error, namelist is full'
      nm%cc = nm%cc + 1

      nm%members(nm%cc)%address = get_address(i)
      nm%members(nm%cc)%type = 2
      nm%members(nm%cc)%i = i0
      nm%members(nm%cc)%i0 = i0
      i = i0
      nm%members(nm%cc)%comment = comment
      nm%members(nm%cc)%label = label
   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine putel_in_namelist_r(nm, i, label, i0, comment)
      TYPE(namelist_set) :: nm
      real(kind=DP)            :: i, i0
      character(*)      :: comment, label
      if (nm%cc > nm%nel) stop 'error, namelist is full'
      nm%cc = nm%cc + 1

      nm%members(nm%cc)%address = get_address(i)
      nm%members(nm%cc)%type = 3
      nm%members(nm%cc)%r = i0
      nm%members(nm%cc)%r0 = i0
      i = i0
      nm%members(nm%cc)%comment = comment
      nm%members(nm%cc)%label = label
   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine putel_in_namelist_rvec(nm, i, label, i0, comment)
      TYPE(namelist_set) :: nm
      real(kind=DP)            :: i(:), i0(:)
      character(*)      :: comment, label
      if (nm%cc > nm%nel) stop 'error, namelist is full'
      if (size(i) /= size(i0)) stop 'error in namelist, vectorial elements: size and size of default values are different'
      nm%cc = nm%cc + 1
      nm%members(nm%cc)%address = get_address(i)
      nm%members(nm%cc)%type = 5
      if (allocated(nm%members(nm%cc)%v)) deallocate (nm%members(nm%cc)%v)
      allocate (nm%members(nm%cc)%v(size(i)))
      if (allocated(nm%members(nm%cc)%v0)) deallocate (nm%members(nm%cc)%v0)
      allocate (nm%members(nm%cc)%v0(size(i0)))
      nm%members(nm%cc)%v0 = i0
      nm%members(nm%cc)%v = i0
      nm%members(nm%cc)%nsize = size(i)
      i = i0
      nm%members(nm%cc)%comment = comment
      nm%members(nm%cc)%label = label
   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine putel_in_namelist_ivec(nm, i, label, i0, comment)
      TYPE(namelist_set) :: nm
      integer(4)         :: i(:), i0(:)
      character(*)      :: comment, label
      if (nm%cc > nm%nel) stop 'error, namelist is full'
      if (size(i) /= size(i0)) stop 'error in namelist, vectorial elements: size and size of default values are different'
      nm%cc = nm%cc + 1
      nm%members(nm%cc)%address = get_address(i)
      nm%members(nm%cc)%type = 8
      if (allocated(nm%members(nm%cc)%ivec)) deallocate (nm%members(nm%cc)%ivec)
      allocate (nm%members(nm%cc)%ivec(size(i)))
      if (allocated(nm%members(nm%cc)%ivec0)) deallocate (nm%members(nm%cc)%ivec0)
      allocate (nm%members(nm%cc)%ivec0(size(i0)))
      nm%members(nm%cc)%ivec0 = i0
      nm%members(nm%cc)%ivec = i0
      nm%members(nm%cc)%nsize = size(i)
      i = i0
      nm%members(nm%cc)%comment = comment
      nm%members(nm%cc)%label = label
   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine putel_in_namelist_cvec(nm, i, label, i0, comment)
      TYPE(namelist_set) :: nm
      complex(kind=DP)         :: i(:), i0(:)
      character(*)      :: comment, label
      if (nm%cc > nm%nel) stop 'error, namelist is full'
      if (size(i) /= size(i0)) stop 'error in namelist, vectorial elements: size and size of default values are different'
      nm%cc = nm%cc + 1
      nm%members(nm%cc)%address = get_address(i)
      nm%members(nm%cc)%type = 7
      if (allocated(nm%members(nm%cc)%vc)) deallocate (nm%members(nm%cc)%vc)
      allocate (nm%members(nm%cc)%vc(size(i)))
      if (allocated(nm%members(nm%cc)%vc0)) deallocate (nm%members(nm%cc)%vc0)
      allocate (nm%members(nm%cc)%vc0(size(i0)))
      nm%members(nm%cc)%vc0 = i0
      nm%members(nm%cc)%vc = i0
      nm%members(nm%cc)%nsize = size(i)
      i = i0
      nm%members(nm%cc)%comment = comment
      nm%members(nm%cc)%label = label
   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine putel_in_namelist_char(nm, i, label, i0, comment)
      TYPE(namelist_set) :: nm
      character(LEN=*)   :: i, i0
      character(LEN=*)   :: comment, label
      if (nm%cc > nm%nel) stop 'error, namelist is full'
      nm%cc = nm%cc + 1
      nm%members(nm%cc)%address = get_address(i)
      nm%members(nm%cc)%type = 6
      nm%members(nm%cc)%a = i0
      nm%members(nm%cc)%a0 = i0
      nm%members(nm%cc)%nsize = len(i)
      if (len(i0) > len(i)) stop 'error putel_in_namelist_char, default character string too long'
      i = i0
      nm%members(nm%cc)%comment = comment
      nm%members(nm%cc)%label = label
   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine putel_in_namelist_rr(nm, i, ib, label, labelb, i0, i0b, comment)
      TYPE(namelist_set) :: nm
      real(kind=DP)            :: i, i0, ib, i0b
      character(*)      :: comment, label, labelb

      if (nm%cc > nm%nel) stop 'error, namelist is full'
      nm%cc = nm%cc + 1

      nm%members(nm%cc)%address = get_address(i)
      nm%members(nm%cc)%type = 3
      nm%members(nm%cc)%r = i0
      nm%members(nm%cc)%r0 = i0
      i = i0
      nm%members(nm%cc)%comment = comment
      nm%members(nm%cc)%label = label

      nm%cc = nm%cc + 1
      if (nm%cc > nm%nel) stop 'error, namelist is full'
      nm%members(nm%cc)%address = get_address(ib)
      nm%members(nm%cc)%type = 3
      nm%members(nm%cc)%r = i0b
      nm%members(nm%cc)%r0 = i0b
      ib = i0b
      nm%members(nm%cc)%comment = comment
      nm%members(nm%cc)%label = labelb

   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine putel_in_namelist_c(nm, i, label, i0, comment)
      TYPE(namelist_set) :: nm
      complex(kind=DP)         :: i, i0
      character(*)      :: comment, label
      if (nm%cc > nm%nel) stop 'error, namelist is full'
      nm%cc = nm%cc + 1

      nm%members(nm%cc)%address = get_address(i)
      nm%members(nm%cc)%type = 4
      nm%members(nm%cc)%c = i0
      nm%members(nm%cc)%c0 = i0
      i = i0
      nm%members(nm%cc)%comment = comment
      nm%members(nm%cc)%label = label
   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine putel_in_namelist_l(nm, i, label, i0, comment)
      TYPE(namelist_set) :: nm
      logical            :: i, i0
      character(*)      :: comment, label
      if (nm%cc > nm%nel) stop 'error, namelist is full'
      nm%cc = nm%cc + 1

      nm%members(nm%cc)%address = get_address(i)
      nm%members(nm%cc)%type = 1
      nm%members(nm%cc)%l = i0
      nm%members(nm%cc)%l0 = i0
      i = i0
      nm%members(nm%cc)%comment = comment
      nm%members(nm%cc)%label = label
   end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

end module

