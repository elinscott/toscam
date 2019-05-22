module string5

   private
   public :: ch_cap
   public :: ch_eqi
   public :: get_unit
   public :: s_eqi
   public :: s_to_r8

contains

   subroutine ch_cap(c)

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
      implicit none

      character c
      integer itemp

      itemp = ichar(c)

      if (97 <= itemp .and. itemp <= 122) then
         c = char(itemp - 32)
      end if

      return
   end subroutine

   function ch_eqi(c1, c2)

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for
!equality.
!
!  Examples:
!
!    CH_EQI ( 'A', 'a' ) is TRUE.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
      implicit none

      character c1
      character c1_cap
      character c2
      character c2_cap
      logical ch_eqi

      c1_cap = c1
      c2_cap = c2

      call ch_cap(c1_cap)
      call ch_cap(c2_cap)

      if (c1_cap == c2_cap) then
         ch_eqi = .true.
      else
         ch_eqi = .false.
      end if

      return
   end function

   subroutine ch_to_digit(c, digit)

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
      implicit none

      character c
      integer digit

      if (lle('0', c) .and. lle(c, '9')) then

         digit = ichar(c) - 48

      else if (c == ' ') then

         digit = 0

      else

         digit = -1

      end if

      return
   end subroutine

   subroutine get_unit(iunit)

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN
!    unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
      implicit none

      integer i
      integer ios
      integer iunit
      logical lopen

      iunit = 0

      do i = 1, 99

         if (i /= 5 .and. i /= 6 .and. i /= 9) then

            inquire (unit=i, opened=lopen, iostat=ios)

            if (ios == 0) then
               if (.not. lopen) then
                  iunit = i
                  return
               end if
            end if

         end if

      end do

      return
   end subroutine

! subroutine s_cat ( s1, s2, s3 )
!
! !*****************************************************************************80
! !
! !! S_CAT concatenates two strings to make a third string.
! !
! !  Modified:
! !
! !    18 September 2000
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, character ( len = * ) S1, the "prefix" string.
! !
! !    Input, character ( len = * ) S2, the "postfix" string.
! !
! !    Output, character ( len = * ) S3, the string made by
! !    concatenating S1 and S2, ignoring any trailing blanks.
! !
!   implicit none
!
!   character ( len = * ) s1
!   character ( len = * ) s2
!   character ( len = * ) s3
!
!   if ( s1 == ' ' .and. s2 == ' ' ) then
!     s3 = ' '
!   else if ( s1 == ' ' ) then
!     s3 = s2
!   else if ( s2 == ' ' ) then
!     s3 = s1
!   else
!     s3 = trim ( s1 ) // trim ( s2 )
!   end if
!
!   return
! end subroutine
!
!
!
!
!
!
!
!
!
   function s_eqi(s1, s2)

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Examples:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
      implicit none

      character c1
      character c2
      integer i
      integer len1
      integer len2
      integer lenc
      logical s_eqi
      character(len=*) s1
      character(len=*) s2

      len1 = len(s1)
      len2 = len(s2)
      lenc = min(len1, len2)

      s_eqi = .false.

      do i = 1, lenc

         c1 = s1(i:i)
         c2 = s2(i:i)
         call ch_cap(c1)
         call ch_cap(c2)

         if (c1 /= c2) then
            return
         end if

      end do

      do i = lenc + 1, len1
         if (s1(i:i) /= ' ') then
            return
         end if
      end do

      do i = lenc + 1, len2
         if (s2(i:i) /= ' ') then
            return
         end if
      end do

      s_eqi = .true.

      return
   end function

   subroutine s_to_r8(s, dval, ierror, length)

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    The routine will read as many characters as possible until it
!    reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Examples:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
      implicit none

      ! logical ch_eqi
      character c
      real(kind=8) dval
      integer ierror
      integer ihave
      integer isgn
      integer iterm
      integer jbot
      integer jsgn
      integer jtop
      integer length
      integer nchar
      integer ndig
      real(kind=8) rbot
      real(kind=8) rexp
      real(kind=8) rtop
      character(len=*) s

      nchar = len_trim(s)

      ierror = 0
      dval = 0.0D+00
      length = -1
      isgn = 1
      rtop = 0
      rbot = 1
      jsgn = 1
      jtop = 0
      jbot = 1
      ihave = 1
      iterm = 0

      do

         length = length + 1

         if (nchar < length + 1) then
            exit
         end if

         c = s(length + 1:length + 1)
!
!  Blank character.
!
         if (c == ' ') then

            if (ihave == 2) then

            else if (ihave == 6 .or. ihave == 7) then
               iterm = 1
            else if (1 < ihave) then
               ihave = 11
            end if
!
!  Comma.
!
         else if (c == ',' .or. c == ';') then

            if (ihave /= 1) then
               iterm = 1
               ihave = 12
               length = length + 1
            end if
!
!  Minus sign.
!
         else if (c == '-') then

            if (ihave == 1) then
               ihave = 2
               isgn = -1
            else if (ihave == 6) then
               ihave = 7
               jsgn = -1
            else
               iterm = 1
            end if
!
!  Plus sign.
!
         else if (c == '+') then

            if (ihave == 1) then
               ihave = 2
            else if (ihave == 6) then
               ihave = 7
            else
               iterm = 1
            end if
!
!  Decimal point.
!
         else if (c == '.') then

            if (ihave < 4) then
               ihave = 4
            else if (6 <= ihave .and. ihave <= 8) then
               ihave = 9
            else
               iterm = 1
            end if
!
!  Scientific notation exponent marker.
!
         else if (ch_eqi(c, 'E') .or. ch_eqi(c, 'D')) then

            if (ihave < 6) then
               ihave = 6
            else
               iterm = 1
            end if
!
!  Digit.
!
         else if (ihave < 11 .and. lle('0', c) .and. lle(c, '9')) then

            if (ihave <= 2) then
               ihave = 3
            else if (ihave == 4) then
               ihave = 5
            else if (ihave == 6 .or. ihave == 7) then
               ihave = 8
            else if (ihave == 9) then
               ihave = 10
            end if

            call ch_to_digit(c, ndig)

            if (ihave == 3) then
               rtop = 10.0D+00*rtop + real(ndig, kind=8)
            else if (ihave == 5) then
               rtop = 10.0D+00*rtop + real(ndig, kind=8)
               rbot = 10.0D+00*rbot
            else if (ihave == 8) then
               jtop = 10*jtop + ndig
            else if (ihave == 10) then
               jtop = 10*jtop + ndig
               jbot = 10*jbot
            end if
!
!  Anything else is regarded as a terminator.
!
         else
            iterm = 1
         end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
         if (iterm == 1) then
            exit
         end if

      end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
      if (iterm /= 1 .and. length + 1 == nchar) then
         length = nchar
      end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
      if (ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7) then
         ierror = ihave
         write (*, '(a)') ' '
         write (*, '(a)') 'S_TO_R8 - Serious error!'
         write (*, '(a)') '  Illegal or nonnumeric input:'
         write (*, '(a)') '    '//trim(s)
         return
      end if
!
!  Number seems OK.  Form it.
!
      if (jtop == 0) then
         rexp = 1.0D+00
      else
         if (jbot == 1) then
            rexp = 10.0D+00**(jsgn*jtop)
         else
            rexp = 10.0D+00**(real(jsgn*jtop, kind=8) &
                              /real(jbot, kind=8))
         end if
      end if

      dval = real(isgn, kind=8)*rexp*rtop/rbot

      return
   end subroutine
!
!
!
!
!
!
!
! subroutine timestamp ( )
!
! !*****************************************************************************80
! !
! !! TIMESTAMP prints the current YMDHMS date as a time stamp.
! !
! !  Example:
! !
! !    31 May 2001   9:45:54.872 AM
! !
! !  Modified:
! !
! !    06 August 2005
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    None
! !
!   implicit none
!
!   character ( len = 8 ) ampm
!   integer d
!   integer h
!   integer m
!   integer mm
!   character ( len = 9 ), parameter, dimension(12) :: month = (/ &
!     'January  ', 'February ', 'March    ', 'April    ', &
!     'May      ', 'June     ', 'July     ', 'August   ', &
!     'September', 'October  ', 'November ', 'December ' /)
!   integer n
!   integer s
!   integer values(8)
!   integer y
!
!   call date_and_time ( values = values )
!
!   y = values(1)
!   m = values(2)
!   d = values(3)
!   h = values(5)
!   n = values(6)
!   s = values(7)
!   mm = values(8)
!
!   if ( h < 12 ) then
!     ampm = 'AM'
!   else if ( h == 12 ) then
!     if ( n == 0 .and. s == 0 ) then
!       ampm = 'Noon'
!     else
!       ampm = 'PM'
!     end if
!   else
!     h = h - 12
!     if ( h < 12 ) then
!       ampm = 'PM'
!     else if ( h == 12 ) then
!       if ( n == 0 .and. s == 0 ) then
!         ampm = 'Midnight'
!       else
!         ampm = 'AM'
!       end if
!     end if
!   end if
!
!   write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
!     d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
!
!   return
! end subroutine
!
!
!
!
!
! subroutine word_next_read ( s, word, done )
!
! !*****************************************************************************80
! !
! !! WORD_NEXT_READ "reads" words from a string, one at a time.
! !
! !  Special cases:
! !
! !    The following characters are considered to be a single word,
! !    whether surrounded by spaces or not:
! !
! !      " ( ) { } [ ]
! !
! !    Also, if there is a trailing comma on the word, it is stripped off.
! !    This is to facilitate the reading of lists.
! !
! !  Modified:
! !
! !    23 May 2001
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    Input, character ( len = * ) S, a string, presumably containing
! !    words
! !    separated by spaces.
! !
! !    Output, character ( len = * ) WORD.
! !    If DONE is FALSE, then WORD contains the "next" word read.
! !    If DONE is TRUE, then WORD is blank, because there was no more to
! !    read.
! !
! !    Input/output, logical DONE.
! !    On input with a fresh string, set DONE to TRUE.
! !    On output, the routine sets DONE:
! !      FALSE if another word was read,
! !      TRUE if no more words could be read.
! !
!   implicit none
!
!   logical done
!   integer ilo
!   integer, save :: lenc = 0
!   integer, save :: next = 1
!   character ( len = * ) s
!   character, parameter :: TAB = char ( 9 )
!   character ( len = * ) word
! !
! !  We "remember" LENC and NEXT from the previous call.
! !
! !  An input value of DONE = TRUE signals a new line of text to examine.
! !
!   if ( done ) then
!
!     next = 1
!     done = .false.
!     lenc = len_trim ( s )
!
!     if ( lenc <= 0 ) then
!       done = .true.
!       word = ' '
!       return
!     end if
!
!   end if
! !
! !  Beginning at index NEXT, search the string for the next nonblank,
! !  which signals the beginning of a word.
! !
!   ilo = next
! !
! !  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
! !
!   do
!
!     if ( lenc < ilo ) then
!       word = ' '
!       done = .true.
!       next = lenc + 1
!       return
!     end if
! !
! !  If the current character is blank, skip to the next one.
! !
!     if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
!       exit
!     end if
!
!     ilo = ilo + 1
!
!   end do
! !
! !  ILO is the index of the next nonblank character in the string.
! !
! !  If this initial nonblank is a special character,
! !  then that's the whole word as far as we're concerned,
! !  so return immediately.
! !
!   if ( s(ilo:ilo) == '"' .or. &
!        s(ilo:ilo) == '(' .or. &
!        s(ilo:ilo) == ')' .or. &
!        s(ilo:ilo) == '{' .or. &
!        s(ilo:ilo) == '}' .or. &
!        s(ilo:ilo) == '[' .or. &
!        s(ilo:ilo) == ']' ) then
!
!     word = s(ilo:ilo)
!     next = ilo + 1
!     return
!
!   end if
! !
! !  Now search for the last contiguous character that is not a
! !  blank, TAB, or special character.
! !
!   next = ilo + 1
!
!   do while ( next <= lenc )
!
!     if ( s(next:next) == ' ' ) then
!       exit
!     else if ( s(next:next) == TAB ) then
!       exit
!     else if ( s(next:next) == '"' ) then
!       exit
!     else if ( s(next:next) == '(' ) then
!       exit
!     else if ( s(next:next) == ')' ) then
!       exit
!     else if ( s(next:next) == '{' ) then
!       exit
!     else if ( s(next:next) == '}' ) then
!       exit
!     else if ( s(next:next) == '[' ) then
!       exit
!     else if ( s(next:next) == ']' ) then
!       exit
!     end if
!
!     next = next + 1
!
!   end do
!
!   if ( s(next-1:next-1) == ',' ) then
!     word = s(ilo:next-2)
!   else
!     word = s(ilo:next-1)
!   end if
!
!   return
! end subroutine

end module
