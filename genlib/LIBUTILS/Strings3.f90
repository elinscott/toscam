! bof
! **********************************************************************
! Fortran 95 module character_functions

! **********************************************************************
! Source Control Strings

! $Id: charfunc.f90 1.8 2003/06/01 13:10:41Z Dan Release $

! **********************************************************************
!  Copyright 2003 Purple Sage Computing Solutions, Inc.

! **********************************************************************
! a set of parameters and functions useful when using character entities

! **********************************************************************
! Summary of License

!   This library is free software; you can redistribute it and/or
!   modify it under the terms of the GNU Library General Public
!   License as published by the Free Software Foundation; either
!   version 2 of the License, or (at your option) any later version.

!   This library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!   Library General Public License for more details.

!   You should have received a copy of the GNU Library General Public
!   License along with this library; if not, write to the Free
!   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

! To report bugs, suggest enhancements, etc. to the Authors,
! Contact:
!    Purple Sage Computing Solutions, Inc.
!                               send email to dnagle@erols.com
!                                   or fax to 703 471 0684 (USA)
!                                  or mail to 12142 Purple Sage Ct.
!                                             Reston, VA 20194-5621 USA

! **********************************************************************

!  character_functions constants

!     ascii_???= defines an integer parameter for each ascii character

!     esc_newline= ascii newline
!     esc_cr= ascii carriage return
!     esc_tab= ascii horizontal tab
!     esc_vtab= ascii vertical tab
!     esc_newpage= ascii form feed
!     esc_escape= ascii escape
!     esc_eof= ascii end of file

!  character_functions operators

!     .eqic.
!     .neic.
!     .ltic.
!     .leic.
!     .geic.
!     .gtic. compare ignoring case
!     .lexeq.
!     .lexne.
!     .lexlt.
!     .lexle.
!     .lexge.
!     .lexgt. compare punct, digit, alphabetic (ignore case)

!  character_functions library

!     isascii() true if ascii character
!     isupper() true if ascii upper case
!     islower() true if ascii lower case
!     isprint() true if ascii printing
!     isgraph() true if ascii graphic
!     isalpha() true if ascii alphabetic
!     isalnum() true if ascii alfanumeric
!     isdigit() true if ascii numeric
!     isxdigit() true if ascii hex digit 0-9, A-F, a-f
!     ispunct() true if ascii punctuation
!     isspace() true if ascii space, tab, cr, newline, vtab, ff
!     iscntl() true if ascii control

!     toupper() character upper case if lower case else noop
!     tolower() character lower case if upper case else noop

!     isevenp() true if character has even parity
!     isoddp() true if character has odd parity
!     ismarkp() true if character has mark parity
!     isclrp() true if character has clear parity

!     toevenp() character with even parity added
!     tooddp() character with odd parity added
!     tomarkp() character with mark parity added
!     toclrp() character with clear parity added

!     rbn(), rnb() convert trailing blanks to/from nulls
!     bcomp() blank compress ignoring quoted strings
!     bclc() to blank compressed lower case
!     bcuc() to blank compressed upper case

!     tr() translates character substrings, arrays by translation table

!     compare() return index of first character to differ

!     find_matching() returns index of matching (), [], or {}
!     unquote_string() returns a quoted string without the quotes

! **********************************************************************

!  character_functions

module character_functions

use genvar

! **********************************************************************
! declare all variables

implicit none                                                        ! turn off implicit typing

! **********************************************************************
! export all names

private                                                              ! turn off implicit exporting

integer,parameter :: substring_not_found=0

! **********************************************************************
! save all names

save                                                                 ! force static module data

! **********************************************************************

!  RCS strings

! **********************************************************************

character( len= *), public, parameter :: character_functions_rcs_id = &
      '$Id: charfunc.f90 1.8 2003/06/01 13:10:41Z Dan Release $'

! **********************************************************************

!  ascii character codes (default integers)

integer, public, parameter :: ascii_nul = 0                          ! null

integer, public, parameter :: ascii_soh = 1                          ! start of heading
integer, public, parameter :: ascii_stx = 2                          ! start of text
integer, public, parameter :: ascii_etx = 3                          ! end of text
integer, public, parameter :: ascii_eot = 4                          ! end of transmission
integer, public, parameter :: ascii_enq = 5                          ! enquiry
integer, public, parameter :: ascii_ack = 6                          ! acknowledge
integer, public, parameter :: ascii_bel = 7                          ! bell
integer, public, parameter :: ascii_bs = 8                           ! backspace
integer, public, parameter :: ascii_ht = 9                           ! horizontal tab
integer, public, parameter :: ascii_lf = 10                          ! line feed
integer, public, parameter :: ascii_vt = 11                          ! vertical tab
integer, public, parameter :: ascii_ff = 12                          ! form feed
integer, public, parameter :: ascii_cr = 13                          ! carriage return
integer, public, parameter :: ascii_s0 = 14                          ! shift out
integer, public, parameter :: ascii_si = 15                          ! shift in
integer, public, parameter :: ascii_dle = 16                         ! data link escape
integer, public, parameter :: ascii_dc1 = 17                         ! device control 1
integer, public, parameter :: ascii_dc2 = 18                         ! device control 2
integer, public, parameter :: ascii_dc3 = 19                         ! device control 3
integer, public, parameter :: ascii_dc4 = 20                         ! device control 4
integer, public, parameter :: ascii_nak = 21                         ! negative acknowledge
integer, public, parameter :: ascii_syn = 22                         ! synchronous idle
integer, public, parameter :: ascii_etb = 23                         ! end of transmission block
integer, public, parameter :: ascii_can = 24                         ! cancel
integer, public, parameter :: ascii_em = 25                          ! end of medium
integer, public, parameter :: ascii_sub = 26                         ! substitute
integer, public, parameter :: ascii_esc = 27                         ! escape
integer, public, parameter :: ascii_fs = 28                          ! file separator
integer, public, parameter :: ascii_gs = 29                          ! group separator
integer, public, parameter :: ascii_rs = 30                          ! record separator
integer, public, parameter :: ascii_us = 31                          ! unit separator

integer, public, parameter :: ascii_sp = 32                          ! space
integer, public, parameter :: ascii_ep = 33                          ! ! exclamation point
integer, public, parameter :: ascii_qtm = 34                         ! " quotation marks
integer, public, parameter :: ascii_ns = 35                          ! # number sign
integer, public, parameter :: ascii_cs = 36                          ! $ currency symbol
integer, public, parameter :: ascii_pct = 37                         ! % percent
integer, public, parameter :: ascii_amp = 38                         ! & ampersand
integer, public, parameter :: ascii_apo = 39                         ! ' apostrophe
integer, public, parameter :: ascii_lp = 40                          ! ( left parenthesis
integer, public, parameter :: ascii_rp = 41                          ! ) right parenthesis
integer, public, parameter :: ascii_ast = 42                         ! * asterisk
integer, public, parameter :: ascii_pls = 43                         ! + plus
integer, public, parameter :: ascii_com = 44                         ! , comma
integer, public, parameter :: ascii_mns = 45                         ! - minus
integer, public, parameter :: ascii_prd = 46                         ! . period
integer, public, parameter :: ascii_sl = 47                          ! / slash

integer, public, parameter :: ascii_0 = 48                           ! 0 zero
integer, public, parameter :: ascii_1 = 49                           ! 1 one
integer, public, parameter :: ascii_2 = 50                           ! 2 two
integer, public, parameter :: ascii_3 = 51                           ! 3 three
integer, public, parameter :: ascii_4 = 52                           ! 4 four
integer, public, parameter :: ascii_5 = 53                           ! 5 five
integer, public, parameter :: ascii_6 = 54                           ! 6 six
integer, public, parameter :: ascii_7 = 55                           ! 7 seven
integer, public, parameter :: ascii_8 = 56                           ! 8 eight
integer, public, parameter :: ascii_9 = 57                           ! 9 nine

integer, public, parameter :: ascii_cl = 58                          ! : colon
integer, public, parameter :: ascii_scl = 59                         ! ; semicolon
integer, public, parameter :: ascii_lt = 60                          ! < less than
integer, public, parameter :: ascii_eq = 61                          ! = equals
integer, public, parameter :: ascii_gt = 62                          ! > greater than
integer, public, parameter :: ascii_qm = 63                          ! ? question mark
integer, public, parameter :: ascii_ats = 64                         ! @ at sign

integer, public, parameter :: ascii_uca = 65                         ! upper case A
integer, public, parameter :: ascii_ucb = 66                         ! upper case B
integer, public, parameter :: ascii_ucc = 67                         ! upper case C
integer, public, parameter :: ascii_ucd = 68                         ! upper case D
integer, public, parameter :: ascii_uce = 69                         ! upper case E
integer, public, parameter :: ascii_ucf = 70                         ! upper case F
integer, public, parameter :: ascii_ucg = 71                         ! upper case G
integer, public, parameter :: ascii_uch = 72                         ! upper case H
integer, public, parameter :: ascii_uci = 73                         ! upper case I
integer, public, parameter :: ascii_ucj = 74                         ! upper case J
integer, public, parameter :: ascii_uck = 75                         ! upper case K
integer, public, parameter :: ascii_ucl = 76                         ! upper case L
integer, public, parameter :: ascii_ucm = 77                         ! upper case M
integer, public, parameter :: ascii_ucn = 78                         ! upper case N
integer, public, parameter :: ascii_uco = 79                         ! upper case O
integer, public, parameter :: ascii_ucp = 80                         ! upper case P
integer, public, parameter :: ascii_ucq = 81                         ! upper case Q
integer, public, parameter :: ascii_ucr = 82                         ! upper case R
integer, public, parameter :: ascii_ucs = 83                         ! upper case S
integer, public, parameter :: ascii_uct = 84                         ! upper case T
integer, public, parameter :: ascii_ucu = 85                         ! upper case U
integer, public, parameter :: ascii_ucv = 86                         ! upper case V
integer, public, parameter :: ascii_ucw = 87                         ! upper case W
integer, public, parameter :: ascii_ucx = 88                         ! upper case X
integer, public, parameter :: ascii_ucy = 89                         ! upper case Y
integer, public, parameter :: ascii_ucz = 90                         ! upper case Z

integer, public, parameter :: ascii_lb = 91                          ! [ left bracket
integer, public, parameter :: ascii_bsl = 92                         ! \ backslash
integer, public, parameter :: ascii_rb = 93                          ! ] right bracket
integer, public, parameter :: ascii_crt = 94                         ! ^ caret
integer, public, parameter :: ascii_und = 95                         ! _ underscore
integer, public, parameter :: ascii_gra = 96                         ! ` grave accent

integer, public, parameter :: ascii_lca = 97                         ! lower case a
integer, public, parameter :: ascii_lcb = 98                         ! lower case b
integer, public, parameter :: ascii_lcc = 99                         ! lower case c
integer, public, parameter :: ascii_lcd = 100                        ! lower case d
integer, public, parameter :: ascii_lce = 101                        ! lower case e
integer, public, parameter :: ascii_lcf = 102                        ! lower case f
integer, public, parameter :: ascii_lcg = 103                        ! lower case g
integer, public, parameter :: ascii_lch = 104                        ! lower case h
integer, public, parameter :: ascii_lci = 105                        ! lower case i
integer, public, parameter :: ascii_lcj = 106                        ! lower case j
integer, public, parameter :: ascii_lck = 107                        ! lower case k
integer, public, parameter :: ascii_lcl = 108                        ! lower case l
integer, public, parameter :: ascii_lcm = 109                        ! lower case m
integer, public, parameter :: ascii_lcn = 110                        ! lower case n
integer, public, parameter :: ascii_lco = 111                        ! lower case o
integer, public, parameter :: ascii_lcp = 112                        ! lower case p
integer, public, parameter :: ascii_lcq = 113                        ! lower case q
integer, public, parameter :: ascii_lcr = 114                        ! lower case r
integer, public, parameter :: ascii_lcs = 115                        ! lower case s
integer, public, parameter :: ascii_lct = 116                        ! lower case t
integer, public, parameter :: ascii_lcu = 117                        ! lower case u
integer, public, parameter :: ascii_lcv = 118                        ! lower case v
integer, public, parameter :: ascii_lcw = 119                        ! lower case w
integer, public, parameter :: ascii_lcx = 120                        ! lower case x
integer, public, parameter :: ascii_lcy = 121                        ! lower case y
integer, public, parameter :: ascii_lcz = 122                        ! lower case z

integer, public, parameter :: ascii_lbr = 123                        ! { left brace
integer, public, parameter :: ascii_vl = 124                         ! | vertical line
integer, public, parameter :: ascii_rbr = 125                        ! } right brace
integer, public, parameter :: ascii_tld = 126                        ! ~ tilde

integer, public, parameter :: ascii_del = 127                        ! delete

! **********************************************************************public,

!  public, escape character

! **********************************************************************

!  ascii newline

character( len= *), public, parameter :: esc_newline = achar( ascii_lf)

!  ascii tab

character( len= *), public, parameter :: esc_tab = achar( ascii_ht)

!  ascii newpage

character( len= *), public, parameter :: esc_newpage = achar( ascii_ff)

!  ascii escape

character( len= *), public, parameter :: esc_escape = achar( ascii_esc)

!  ascii eof

character( len= *), public, parameter :: esc_eof = achar( ascii_sub)

! **********************************************************************

!  local data

! **********************************************************************

!  difference between uppercase and lowercase (in ascii)

integer, parameter :: change_case = 32

!  biases for the compare (ignoring case) and the compare (lexical) functions

integer, parameter :: bias_digit = 256

integer, parameter :: bias_letter = 512

!  parity bit (2^7) Z'80'

integer, parameter :: parity_bit = 128

!  lower seven bits set Z'7f'

integer, parameter :: ascii_bits = 127

!  true if parity is even

logical, dimension( 0: 127), parameter :: even_parity = (/ &
                                   .true.,  &                        ! 0000000
                                   .false., &                        ! 0000001
                                   .false., &                        ! 0000010
                                   .true.,  &                        ! 0000011
                                   .false., &                        ! 0000100
                                   .true.,  &                        ! 0000101
                                   .true.,  &                        ! 0000110
                                   .false., &                        ! 0000111
                                   .false., &                        ! 0001000
                                   .true.,  &                        ! 0001001
                                   .true.,  &                        ! 0001010
                                   .false., &                        ! 0001011
                                   .true.,  &                        ! 0001100
                                   .false., &                        ! 0001101
                                   .false., &                        ! 0001110
                                   .true.,  &                        ! 0001111
                                   .false., &                        ! 0010000
                                   .true.,  &                        ! 0010001
                                   .true.,  &                        ! 0010010
                                   .false., &                        ! 0010011
                                   .true.,  &                        ! 0010100
                                   .false., &                        ! 0010101
                                   .false., &                        ! 0010110
                                   .true.,  &                        ! 0010111
                                   .true.,  &                        ! 0011000
                                   .false., &                        ! 0011001
                                   .false., &                        ! 0011010
                                   .true.,  &                        ! 0011011
                                   .false., &                        ! 0011100
                                   .true.,  &                        ! 0011101
                                   .true.,  &                        ! 0011110
                                   .false., &                        ! 0011111
                                   .false., &                        ! 0100000
                                   .true.,  &                        ! 0100001
                                   .true.,  &                        ! 0100010
                                   .false., &                        ! 0100011
                                   .true.,  &                        ! 0100100
                                   .false., &                        ! 0100101
                                   .false., &                        ! 0100110
                                   .true.,  &                        ! 0100111
                                   .true.,  &                        ! 0101000
                                   .false., &                        ! 0101001
                                   .false., &                        ! 0101010
                                   .true.,  &                        ! 0101011
                                   .false., &                        ! 0101100
                                   .true.,  &                        ! 0101101
                                   .true.,  &                        ! 0101110
                                   .false., &                        ! 0101111
                                   .true.,  &                        ! 0110000
                                   .false., &                        ! 0110001
                                   .false., &                        ! 0110010
                                   .true.,  &                        ! 0110011
                                   .false., &                        ! 0110100
                                   .true.,  &                        ! 0110101
                                   .true.,  &                        ! 0110110
                                   .false., &                        ! 0110111
                                   .false., &                        ! 0111000
                                   .true.,  &                        ! 0111001
                                   .true.,  &                        ! 0111010
                                   .false., &                        ! 0111011
                                   .true.,  &                        ! 0111100
                                   .false., &                        ! 0111101
                                   .false., &                        ! 0111110
                                   .true.,  &                        ! 0111111
                                   .false., &                        ! 1000000
                                   .true.,  &                        ! 1000001
                                   .true.,  &                        ! 1000010
                                   .false., &                        ! 1000011
                                   .true.,  &                        ! 1000100
                                   .false., &                        ! 1000101
                                   .false., &                        ! 1000110
                                   .true.,  &                        ! 1000111
                                   .true.,  &                        ! 1001000
                                   .false., &                        ! 1001001
                                   .false., &                        ! 1001010
                                   .true.,  &                        ! 1001011
                                   .false., &                        ! 1001100
                                   .true.,  &                        ! 1001101
                                   .true.,  &                        ! 1001110
                                   .false., &                        ! 1001111
                                   .true.,  &                        ! 1010000
                                   .false., &                        ! 1010001
                                   .false., &                        ! 1010010
                                   .true.,  &                        ! 1010011
                                   .false., &                        ! 1010100
                                   .true.,  &                        ! 1010101
                                   .true.,  &                        ! 1010110
                                   .false., &                        ! 1010111
                                   .false., &                        ! 1011000
                                   .true.,  &                        ! 1011001
                                   .true.,  &                        ! 1011010
                                   .false., &                        ! 1011011
                                   .true.,  &                        ! 1011100
                                   .false., &                        ! 1011101
                                   .false., &                        ! 1011110
                                   .true.,  &                        ! 1011111
                                   .true.,  &                        ! 1100000
                                   .false., &                        ! 1100001
                                   .false., &                        ! 1100010
                                   .true.,  &                        ! 1100011
                                   .false., &                        ! 1100100
                                   .true.,  &                        ! 1100101
                                   .true.,  &                        ! 1100110
                                   .false., &                        ! 1100111
                                   .false., &                        ! 1101000
                                   .true.,  &                        ! 1101001
                                   .true.,  &                        ! 1101010
                                   .false., &                        ! 1101011
                                   .true.,  &                        ! 1101100
                                   .false., &                        ! 1101101
                                   .false., &                        ! 1101110
                                   .true.,  &                        ! 1101111
                                   .false., &                        ! 1110000
                                   .true.,  &                        ! 1110001
                                   .true.,  &                        ! 1110010
                                   .false., &                        ! 1110011
                                   .true.,  &                        ! 1110100
                                   .false., &                        ! 1110101
                                   .false., &                        ! 1110110
                                   .true.,  &                        ! 1110111
                                   .true.,  &                        ! 1111000
                                   .false., &                        ! 1111001
                                   .false., &                        ! 1111010
                                   .true.,  &                        ! 1111011
                                   .false., &                        ! 1111100
                                   .true.,  &                        ! 1111101
                                   .true.,  &                        ! 1111110
                                   .false. /)                        ! 1111111

!  characters used repeatedly

character( len= *), parameter :: blank = achar( ascii_sp)

character( len= *), parameter :: single_quote = achar( ascii_apo)

character( len= *), parameter :: double_quote = achar( ascii_qtm)

!  null string to initialize output strings

character( len= *), parameter :: null_string = ''

! **********************************************************************

!  character_functions library

! **********************************************************************

!  declare specific functions implementing the .eqic., etc operators

public :: operator( .eqic.)                                          ! operator name

interface operator( .eqic.)
   module procedure ascii_eqic
end interface

public :: operator( .neic.)                                          ! operator name

interface operator( .neic.)
   module procedure ascii_neic
end interface

public :: operator( .ltic.)                                          ! operator name

interface operator( .ltic.)
   module procedure ascii_ltic
end interface

public :: operator( .leic.)                                          ! operator name

interface operator( .leic.)
   module procedure ascii_leic
end interface

public :: operator( .geic.)                                          ! operator name

interface operator( .geic.)
   module procedure ascii_geic
end interface

public :: operator( .gtic.)                                          ! operator name

interface operator( .gtic.)
   module procedure ascii_gtic
end interface

! **********************************************************************

!  declare specific functions implementing the .lexeq., etc operators

public :: operator( .lexeq.)                                         ! operator name

interface operator( .lexeq.)
   module procedure ascii_lexeq
end interface

public :: operator( .lexne.)                                         ! operator name

interface operator( .lexne.)
   module procedure ascii_lexne
end interface

public :: operator( .lexlt.)                                         ! operator name

interface operator( .lexlt.)
   module procedure ascii_lexlt
end interface

public :: operator( .lexle.)                                         ! operator name

interface operator( .lexle.)
   module procedure ascii_lexle
end interface

public :: operator( .lexge.)                                         ! operator name

interface operator( .lexge.)
   module procedure ascii_lexge
end interface

public :: operator( .lexgt.)                                         ! operator name

interface operator( .lexgt.)
   module procedure ascii_lexgt
end interface

! **********************************************************************

!  character utility functions

! **********************************************************************

!  declare specific functions implementing the isascii() function

public :: isascii                                                    ! generic name

interface isascii
   !module procedure byte_isascii
   !module procedure short_isascii
    module procedure int_isascii
   !module procedure long_isascii
    module procedure character_isascii
end interface

! **********************************************************************

!  character translation functions

! **********************************************************************

!  declare specific functions implementing the tr() subroutine

public :: tr                                                         ! generic name

interface tr
   module procedure array_tr
   module procedure substring_tr
end interface

! **********************************************************************

!  make the other functions public

! **********************************************************************

public :: isupper                                                    ! export names
public :: islower
public :: isprint
public :: isgraph
public :: isalpha
public :: isalnum
public :: isdigit
public :: isxdigit
public :: ispunct
public :: isspace
public :: iscntl

public :: toupper
public :: tolower

public :: isevenp
public :: isoddp
public :: ismarkp
public :: isclrp

public :: toevenp
public :: tooddp
public :: tomarkp
public :: toclrp

public :: rbn
public :: rnb
public :: bcomp
public :: bclc
public :: bcuc

public :: compare_string
public :: replace_substring
public :: string_cat

public :: matching_pair
public :: unquote_string

! **********************************************************************

!  character_functions module procedures

! **********************************************************************

contains

! **********************************************************************

!  .eqic., .neic., .ltic., .leic., .geic., .gtic. compare ignoring case

! **********************************************************************

!  ascii_eqic()

elemental logical function ascii_eqic( string_a, string_b)

!  ascii_eqic() interface

character( len= *), intent( in) :: string_a, string_b

! **********************************************************************

!  ascii_eqic() local

   integer( kind= int_k) :: int_a                                    ! from string_a

   integer( kind= int_k) :: int_b                                    ! from string_b

   integer( kind= int_k) :: sort_a                                   ! bias applied to a

   integer( kind= int_k) :: sort_b                                   ! bias applied to b

   integer( kind= int_k) :: len_i                                    ! loop variable

! **********************************************************************

!  ascii_eqic() text

continue                                                             ! .eqic.

! ----------------------------------------------------------------------

!  loop thru each character

   comp_ab: do len_i = 1, min( len( string_a), len( string_b))

! ----------------------------------------------------------------------

!  boost character from a?

      int_a = iachar( string_a( len_i: len_i) )                      ! get ascii character code

      bias_a: select case( int_a)                                    ! which character

      case( ascii_uca: ascii_ucz) bias_a                             ! upper case

         sort_a = change_case                                        ! boost

      case default bias_a                                            ! otherwise

         sort_a = 0                                                  ! ignore

      end select bias_a                                              ! which character

! ----------------------------------------------------------------------

!  boost character from b?

      int_b = iachar( string_b( len_i: len_i) )                      ! get ascii character code

      bias_b: select case( int_b)                                    ! which character

      case( ascii_uca: ascii_ucz) bias_b                             ! upper case

         sort_b = change_case                                        ! boost

      case default bias_b                                            ! otherwise

         sort_b = 0                                                  ! ignore

      end select bias_b                                              ! which character

! ----------------------------------------------------------------------

!  compare a ? b

      comp_eq: if( int_a + sort_a /= int_b + sort_b )then            ! biased comparison

         ascii_eqic = .false.                                        ! not equal

         return                                                      ! .eqic.

      endif comp_eq                                                  ! biased comparison

   enddo comp_ab                                                     ! do each character

! ----------------------------------------------------------------------

!  last chance comparison

   ascii_eqic = len( string_a) == len( string_b)                     ! lengths differ?

return                                                               ! .eqic.

! **********************************************************************

!  ascii_eqic()

end function ascii_eqic

! **********************************************************************

!  ascii_neic()

elemental logical function ascii_neic( string_a, string_b)

!  ascii_neic() interface

character( len= *), intent( in) :: string_a, string_b

! **********************************************************************

!  ascii_neic() local

   integer( kind= int_k) :: int_a                                    ! from string_a

   integer( kind= int_k) :: int_b                                    ! from string_b

   integer( kind= int_k) :: sort_a                                   ! bias applied to a

   integer( kind= int_k) :: sort_b                                   ! bias applied to b

   integer( kind= int_k) :: len_i                                    ! loop variable

! **********************************************************************

!  ascii_neic()

continue                                                             ! .neic.

! ----------------------------------------------------------------------

!  loop thru each character

   comp_ab: do len_i = 1, min( len( string_a), len( string_b))

! ----------------------------------------------------------------------

!  boost character from a?

      int_a = iachar( string_a( len_i: len_i) )                      ! get ascii character code

      bias_a: select case( int_a)                                    ! which character

      case( ascii_uca: ascii_ucz) bias_a                             ! upper case

         sort_a = change_case                                        ! boost

      case default bias_a                                            ! otherwise

         sort_a = 0                                                  ! ignore

      end select bias_a                                              ! which character

! ----------------------------------------------------------------------

!  boost character from b?

      int_b = iachar( string_b( len_i: len_i) )                      ! get ascii character code

      bias_b: select case( int_b)                                    ! which character

      case( ascii_uca: ascii_ucz) bias_b                             ! upper case

         sort_b = change_case                                        ! boost

      case default bias_b                                            ! otherwise

         sort_b = 0                                                  ! ignore

      end select bias_b                                              ! which character

! ----------------------------------------------------------------------

!  compare a ? b

      comp_ne: if( int_a + sort_a /= int_b + sort_b )then            ! biased comparison

         ascii_neic = .true.                                         ! a is not b

         return                                                      ! .neic.

      endif comp_ne                                                  ! biased comparison

! ----------------------------------------------------------------------

   enddo comp_ab

   ascii_neic = len( string_a) /= len( string_b)                     ! lengths differ?

! ----------------------------------------------------------------------

return                                                               ! .neic.

! **********************************************************************

!  ascii_neic()

end function ascii_neic

! **********************************************************************

!  ascii_ltic()

elemental logical function ascii_ltic( string_a, string_b)

!  ascii_ltic() interface

character( len= *), intent( in) :: string_a, string_b

! **********************************************************************

!  ascii_ltic() local

   integer( kind= int_k) :: int_a                                    ! from string_a

   integer( kind= int_k) :: int_b                                    ! from string_b

   integer( kind= int_k) :: sort_a                                   ! bias applied to a

   integer( kind= int_k) :: sort_b                                   ! bias applied to b

   integer( kind= int_k) :: len_i                                    ! loop variable

! **********************************************************************

!  ascii_ltic()

continue                                                             ! .ltic.

! ----------------------------------------------------------------------

!  loop thru each character

   comp_ab: do len_i = 1, min( len( string_a), len( string_b))

! ----------------------------------------------------------------------

      int_a = iachar( string_a( len_i: len_i) )                      ! get ascii character code

      bias_a: select case( int_a)                                    ! which character

      case( ascii_uca: ascii_ucz) bias_a                             ! upper case

         sort_a = change_case                                        ! boost

      case default bias_a                                            ! otherwise

         sort_a = 0                                                  ! ignore

      end select bias_a                                              ! which character

! ----------------------------------------------------------------------

      int_b = iachar( string_b( len_i: len_i) )                      ! get ascii character code

      bias_b: select case( int_b)                                    ! which character

      case( ascii_uca: ascii_ucz) bias_b                             ! upper case

         sort_b = change_case                                        ! boost

      case default bias_b                                            ! otherwise

         sort_b = 0                                                  ! ignore

      end select bias_b                                              ! which character

! ----------------------------------------------------------------------

!  compare a ? b

      comp_lt: if( int_a + sort_a < int_b + sort_b )then             ! biased comparison

         ascii_ltic = .true.                                         ! a less than b

         return                                                      ! .ltic.

      elseif( int_a + sort_a > int_b + sort_b )then comp_lt

         ascii_ltic = .false.                                        ! a greater than b

         return                                                      ! .ltic.

      endif comp_lt                                                  ! biased comparison

   enddo comp_ab

! ----------------------------------------------------------------------

   ascii_ltic = len( string_a) < len( string_b)                      ! lengths differ?

! ----------------------------------------------------------------------

return                                                               ! .ltic.

! **********************************************************************

!  ascii_ltic()

end function ascii_ltic

! **********************************************************************

!  ascii_leic()

elemental logical function ascii_leic( string_a, string_b)

!  ascii_leic() interface

character( len= *), intent( in) :: string_a, string_b

! **********************************************************************

!  ascii_leic() local

   integer( kind= int_k) :: int_a                                    ! from string_a

   integer( kind= int_k) :: int_b                                    ! from string_b

   integer( kind= int_k) :: sort_a                                   ! bias applied to a

   integer( kind= int_k) :: sort_b                                   ! bias applied to b

   integer( kind= int_k) :: len_i                                    ! loop variable

! **********************************************************************

!  ascii_leic()

continue                                                             ! .leic.

! ----------------------------------------------------------------------

!  loop thru each character

   comp_ab: do len_i = 1, min( len( string_a), len( string_b))

! ----------------------------------------------------------------------

      int_a = iachar( string_a( len_i: len_i) )                      ! get ascii character code

      bias_a: select case( int_a)                                    ! which character

      case( ascii_uca: ascii_ucz) bias_a                             ! upper case

         sort_a = change_case                                        ! boost

      case default bias_a                                            ! otherwise

         sort_a = 0                                                  ! ignore

      end select bias_a                                              ! which character

! ----------------------------------------------------------------------

      int_b = iachar( string_b( len_i: len_i) )                      ! get ascii character code

      bias_b: select case( int_b)                                    ! which character

      case( ascii_uca: ascii_ucz) bias_b                             ! upper case

         sort_b = change_case                                        ! boost

      case default bias_b                                            ! otherwise

         sort_b = 0                                                  ! ignore

      end select bias_b                                              ! which character

! ----------------------------------------------------------------------

      comp_le: if( int_a + sort_a < int_b + sort_b )then

         ascii_leic = .true.

         return                                                      ! .leic.

      elseif( int_a + sort_a > int_b + sort_b )then comp_le

         ascii_leic = .false.

         return                                                      ! .leic.

      endif comp_le

   enddo comp_ab

! ----------------------------------------------------------------------

   ascii_leic = len( string_a) <= len( string_b)                     ! lengths differ?

! ----------------------------------------------------------------------

return                                                               ! .leic.

! **********************************************************************

!  ascii_leic()

end function ascii_leic

! **********************************************************************

!  ascii_geic()

elemental logical function ascii_geic( string_a, string_b)

!  ascii_geic() interface

character( len= *), intent( in) :: string_a, string_b

! **********************************************************************

!  ascii_geic() local

   integer( kind= int_k) :: int_a                                    ! from string_a

   integer( kind= int_k) :: int_b                                    ! from string_b

   integer( kind= int_k) :: sort_a                                   ! bias applied to a

   integer( kind= int_k) :: sort_b                                   ! bias applied to b

   integer( kind= int_k) :: len_i                                    ! loop variable

! **********************************************************************

!  ascii_geic()

continue                                                             ! .geic.

! ----------------------------------------------------------------------

!  loop thru each character

   comp_ab: do len_i = 1, min( len( string_a), len( string_b))

! ----------------------------------------------------------------------

      int_a = iachar( string_a( len_i: len_i) )                      ! get ascii character code

      bias_a: select case( int_a)                                    ! which character

      case( ascii_uca: ascii_ucz) bias_a                             ! upper case

         sort_a = change_case                                        ! boost

      case default bias_a                                            ! otherwise

         sort_a = 0                                                  ! ignore

      end select bias_a                                              ! which character

! ----------------------------------------------------------------------

      int_b = iachar( string_b( len_i: len_i) )                      ! get ascii character code

      bias_b: select case( int_b)                                    ! which character

      case( ascii_uca: ascii_ucz) bias_b                             ! upper case

         sort_b = change_case                                        ! boost

      case default bias_b                                            ! otherwise

         sort_b = 0                                                  ! ignore

      end select bias_b                                              ! which character

! ----------------------------------------------------------------------

      comp_ge: if( int_a + sort_a > int_b + sort_b )then

         ascii_geic = .true.

         return                                                      ! .geic.

      elseif( int_a + sort_a < int_b + sort_b )then comp_ge

         ascii_geic = .false.

         return                                                      ! .geic.

      endif comp_ge

   enddo comp_ab

! ----------------------------------------------------------------------

   ascii_geic = len( string_a) >= len( string_b)                     ! lengths differ?

! ----------------------------------------------------------------------

return                                                               ! .geic.

! **********************************************************************

!  ascii_geic()

end function ascii_geic

! **********************************************************************

!  ascii_gtic()

elemental logical function ascii_gtic( string_a, string_b)

!  ascii_gtic() interface

character( len= *), intent( in) :: string_a, string_b

! **********************************************************************

!  ascii_gtic() local

   integer( kind= int_k) :: int_a                                    ! from string_a

   integer( kind= int_k) :: int_b                                    ! from string_b

   integer( kind= int_k) :: sort_a                                   ! bias applied to a

   integer( kind= int_k) :: sort_b                                   ! bias applied to b

   integer( kind= int_k) :: len_i                                    ! loop variable

! **********************************************************************

!  ascii_gtic()

continue                                                             ! .gtic.

! ----------------------------------------------------------------------

!  loop thru each character

   comp_ab: do len_i = 1, min( len( string_a), len( string_b))

! ----------------------------------------------------------------------

      int_a = iachar( string_a( len_i: len_i) )                      ! get ascii character code

      bias_a: select case( int_a)                                    ! which character

      case( ascii_uca: ascii_ucz) bias_a                             ! upper case

         sort_a = change_case                                        ! boost

      case default bias_a                                            ! otherwise

         sort_a = 0                                                  ! ignore

      end select bias_a                                              ! which character

! ----------------------------------------------------------------------

      int_b = iachar( string_b( len_i: len_i) )                      ! get ascii character code

      bias_b: select case( int_b)                                    ! which character

      case( ascii_uca: ascii_ucz) bias_b                             ! upper case

         sort_b = change_case                                        ! boost

      case default bias_b                                            ! otherwise

         sort_b = 0                                                  ! ignore

      end select bias_b                                              ! which character

! ----------------------------------------------------------------------

      comp_gt: if( int_a + sort_a > int_b + sort_b )then

         ascii_gtic = .true.

         return                                                      ! .gtic.

      endif comp_gt

   enddo comp_ab

! ----------------------------------------------------------------------

   ascii_gtic = len( string_a) > len( string_b)                      ! lengths differ?

! ----------------------------------------------------------------------

return                                                               ! .gtic.

! **********************************************************************

!  ascii_gtic()

end function ascii_gtic

! **********************************************************************

!  .lexeq., .lexne., .lexlt., .lexle., .lexge., .lexgt. lexical compare

! **********************************************************************

!  ascii_lexeq()

elemental logical function ascii_lexeq( string_a, string_b)

!  ascii_lexeq() interface

character( len= *), intent( in) :: string_a, string_b

! **********************************************************************

!  ascii_lexeq() local

   integer( kind= int_k) :: int_a                                    ! from string_a

   integer( kind= int_k) :: int_b                                    ! from string_b

   integer( kind= int_k) :: sort_a                                   ! bias applied to a

   integer( kind= int_k) :: sort_b                                   ! bias applied to b

   integer( kind= int_k) :: len_i                                    ! loop variable

! **********************************************************************

!  ascii_lexeq()

continue                                                             ! .lexeq.

! ----------------------------------------------------------------------

!  loop thru each character

   comp_ab: do len_i = 1, min( len( string_a), len( string_b))

! ----------------------------------------------------------------------

      int_a = iachar( string_a( len_i: len_i) )                      ! get ascii character code

      bias_a: select case( int_a)                                    ! which character

      case( ascii_0: ascii_9) bias_a                                 ! digit

         sort_a = bias_digit                                         ! boost

      case( ascii_uca: ascii_ucz) bias_a                             ! upper case

         sort_a = bias_letter + change_case                          ! boost

      case( ascii_lca: ascii_lcz) bias_a                             ! lower case

         sort_a = bias_letter                                        ! boost

      case default bias_a                                            ! otherwise

         sort_a = 0                                                  ! ignore

      end select bias_a                                              ! which character

! ----------------------------------------------------------------------

      int_b = iachar( string_b( len_i: len_i) )                      ! get ascii character code

      bias_b: select case( int_b)                                    ! which character

      case( ascii_0: ascii_9) bias_b                                 ! digit

         sort_b = bias_digit                                         ! boost

      case( ascii_uca: ascii_ucz) bias_b                             ! upper case

         sort_b = bias_letter + change_case                          ! boost

      case( ascii_lca: ascii_lcz) bias_b                             ! lower case

         sort_b = bias_letter                                        ! boost

      case default bias_b                                            ! otherwise

         sort_b = 0                                                  ! ignore

      end select bias_b                                              ! which character

! ----------------------------------------------------------------------

      comp_eq: if( int_a + sort_a /= int_b + sort_b )then

         ascii_lexeq = .false.

         return                                                      ! .lexeq.

      endif comp_eq

   enddo comp_ab

! ----------------------------------------------------------------------

   ascii_lexeq = len( string_a) == len( string_b)                    ! lengths differ?

! ----------------------------------------------------------------------

return                                                               ! .lexeq.

! **********************************************************************

!  ascii_lexeq()

end function ascii_lexeq

! **********************************************************************

!  ascii_lexne()

elemental logical function ascii_lexne( string_a, string_b)

!  ascii_lexne() interface

character( len= *), intent( in) :: string_a, string_b

! **********************************************************************

!  ascii_lexne() local

   integer( kind= int_k) :: int_a                                    ! from string_a

   integer( kind= int_k) :: int_b                                    ! from string_b

   integer( kind= int_k) :: sort_a                                   ! bias applied to a

   integer( kind= int_k) :: sort_b                                   ! bias applied to b

   integer( kind= int_k) :: len_i                                    ! loop variable

! **********************************************************************

!  ascii_lexne()

continue                                                             ! .lexne.

! ----------------------------------------------------------------------

!  loop thru each character

   comp_ab: do len_i = 1, min( len( string_a), len( string_b))

! ----------------------------------------------------------------------

      int_a = iachar( string_a( len_i: len_i) )                      ! get ascii character code

      bias_a: select case( int_a)                                    ! which character

      case( ascii_0: ascii_9) bias_a                                 ! digit

         sort_a = bias_digit                                         ! boost

      case( ascii_uca: ascii_ucz) bias_a                             ! upper case

         sort_a = bias_letter + change_case                          ! boost

      case( ascii_lca: ascii_lcz) bias_a                             ! lower case

         sort_a = bias_letter                                        ! boost

      case default bias_a                                            ! otherwise

         sort_a = 0                                                  ! ignore

      end select bias_a                                              ! which character

! ----------------------------------------------------------------------

      int_b = iachar( string_b( len_i: len_i) )                      ! get ascii character code

      bias_b: select case( int_b)                                    ! which character

      case( ascii_0: ascii_9) bias_b                                 ! digit

         sort_b = bias_digit                                         ! boost

      case( ascii_uca: ascii_ucz) bias_b                             ! upper case

         sort_b = bias_letter + change_case                          ! boost

      case( ascii_lca: ascii_lcz) bias_b                             ! lower case

         sort_b = bias_letter                                        ! boost

      case default bias_b                                            ! otherwise

         sort_b = 0                                                  ! ignore

      end select bias_b                                              ! which character

! ----------------------------------------------------------------------

      comp_ne: if( int_a + sort_a /= int_b + sort_b )then

         ascii_lexne = .true.

         return                                                      ! .lexne.

      endif comp_ne

   enddo comp_ab

! ----------------------------------------------------------------------

   ascii_lexne = len( string_a) /= len( string_b)                    ! lengths differ?

! ----------------------------------------------------------------------

return                                                               ! .lexne.

! **********************************************************************

!  ascii_lexne()

end function ascii_lexne

! **********************************************************************

!  ascii_lexlt()

elemental logical function ascii_lexlt( string_a, string_b)

!  ascii_lexlt() interface

character( len= *), intent( in) :: string_a, string_b

! **********************************************************************

!  ascii_lexlt() local

   integer( kind= int_k) :: int_a                                    ! from string_a

   integer( kind= int_k) :: int_b                                    ! from string_b

   integer( kind= int_k) :: sort_a                                   ! bias applied to a

   integer( kind= int_k) :: sort_b                                   ! bias applied to b

   integer( kind= int_k) :: len_i                                    ! loop variable

! **********************************************************************

!  ascii_lexlt()

continue                                                             ! .lexlt.

! ----------------------------------------------------------------------

!  loop thru each character

   comp_ab: do len_i = 1, min( len( string_a), len( string_b))

! ----------------------------------------------------------------------

      int_a = iachar( string_a( len_i: len_i) )                      ! get ascii character code

      bias_a: select case( int_a)                                    ! which character

      case( ascii_0: ascii_9) bias_a                                 ! digit

         sort_a = bias_digit                                         ! boost

      case( ascii_uca: ascii_ucz) bias_a                             ! upper case

         sort_a = bias_letter + change_case                          ! boost

      case( ascii_lca: ascii_lcz) bias_a                             ! lower case

         sort_a = bias_letter                                        ! boost

      case default bias_a                                            ! otherwise

         sort_a = 0                                                  ! ignore

      end select bias_a                                              ! which character

! ----------------------------------------------------------------------

      int_b = iachar( string_b( len_i: len_i) )                      ! get ascii character code

      bias_b: select case( int_b)                                    ! which character

      case( ascii_0: ascii_9) bias_b                                 ! digit

         sort_b = bias_digit                                         ! boost

      case( ascii_uca: ascii_ucz) bias_b                             ! upper case

         sort_b = bias_letter + change_case                          ! boost

      case( ascii_lca: ascii_lcz) bias_b                             ! lower case

         sort_b = bias_letter                                        ! boost

      case default bias_b                                            ! otherwise

         sort_b = 0                                                  ! ignore

      end select bias_b                                              ! which character

! ----------------------------------------------------------------------

      comp_lt: if( int_a + sort_a < int_b + sort_b )then

         ascii_lexlt = .true.

         return                                                      ! .lexlt.

      elseif( int_a + sort_a > int_b + sort_b )then comp_lt

         ascii_lexlt = .false.

         return                                                      ! .lexlt.

      endif comp_lt

   enddo comp_ab

! ----------------------------------------------------------------------

   ascii_lexlt = len( string_a) < len( string_b)                     ! lengths differ?

! ----------------------------------------------------------------------

return                                                               ! .lexlt.

! **********************************************************************

!  ascii_lexlt()

end function ascii_lexlt

! **********************************************************************

!  ascii_lexle()

elemental logical function ascii_lexle( string_a, string_b)

!  ascii_lexle() interface

character( len= *), intent( in) :: string_a, string_b

! **********************************************************************

!  ascii_lexle() local

   integer( kind= int_k) :: int_a                                    ! from string_a

   integer( kind= int_k) :: int_b                                    ! from string_b

   integer( kind= int_k) :: sort_a                                   ! bias applied to a

   integer( kind= int_k) :: sort_b                                   ! bias applied to b

   integer( kind= int_k) :: len_i                                    ! loop variable

! **********************************************************************

!  ascii_lexle()

continue                                                             ! .lexle.

! ----------------------------------------------------------------------

!  loop thru each character

   comp_ab: do len_i = 1, min( len( string_a), len( string_b))

! ----------------------------------------------------------------------

      int_a = iachar( string_a( len_i: len_i) )                      ! get ascii character code

      bias_a: select case( int_a)                                    ! which character

      case( ascii_0: ascii_9) bias_a                                 ! digit

         sort_a = bias_digit                                         ! boost

      case( ascii_uca: ascii_ucz) bias_a                             ! upper case

         sort_a = bias_letter + change_case                          ! boost

      case( ascii_lca: ascii_lcz) bias_a                             ! lower case

         sort_a = bias_letter                                        ! boost

      case default bias_a                                            ! otherwise

         sort_a = 0                                                  ! ignore

      end select bias_a                                              ! which character

! ----------------------------------------------------------------------

      int_b = iachar( string_b( len_i: len_i) )                      ! get ascii character code

      bias_b: select case( int_b)                                    ! which character

      case( ascii_0: ascii_9) bias_b                                 ! digit

         sort_b = bias_digit                                         ! boost

      case( ascii_uca: ascii_ucz) bias_b                             ! upper case

         sort_b = bias_letter + change_case                          ! boost

      case( ascii_lca: ascii_lcz) bias_b                             ! lower case

         sort_b = bias_letter                                        ! boost

      case default bias_b                                            ! otherwise

         sort_b = 0                                                  ! ignore

      end select bias_b                                              ! which character

! ----------------------------------------------------------------------

      comp_le: if( int_a + sort_a < int_b + sort_b )then

         ascii_lexle = .true.

         return                                                      ! .lexle.

      elseif( int_a + sort_b > int_b + sort_b )then comp_le

         ascii_lexle = .false.

         return                                                      ! .lexle.

      endif comp_le

   enddo comp_ab

! ----------------------------------------------------------------------

   ascii_lexle = len( string_a) <= len( string_b)                    ! lengths differ?

! ----------------------------------------------------------------------

return                                                               ! .lexle.

! **********************************************************************

!  ascii_lexle()

end function ascii_lexle

! **********************************************************************

!  ascii_lexge()

elemental logical function ascii_lexge( string_a, string_b)

!  ascii_lexge() interface

character( len= *), intent( in) :: string_a, string_b

! **********************************************************************

!  ascii_lexge() local

   integer( kind= int_k) :: int_a                                    ! from string_a

   integer( kind= int_k) :: int_b                                    ! from string_b

   integer( kind= int_k) :: sort_a                                   ! bias applied to a

   integer( kind= int_k) :: sort_b                                   ! bias applied to b

   integer( kind= int_k) :: len_i                                    ! loop variable

! **********************************************************************

!  ascii_lexge()

continue                                                             ! .lexge.

! ----------------------------------------------------------------------

!  loop thru each character

   comp_ab: do len_i = 1, min( len( string_a), len( string_b))

! ----------------------------------------------------------------------

      int_a = iachar( string_a( len_i: len_i) )                      ! get ascii character code

      bias_a: select case( int_a)                                    ! which character

      case( ascii_0: ascii_9) bias_a                                 ! digit

         sort_a = bias_digit                                         ! boost

      case( ascii_uca: ascii_ucz) bias_a                             ! upper case

         sort_a = bias_letter + change_case                          ! boost

      case( ascii_lca: ascii_lcz) bias_a                             ! lower case

         sort_a = bias_letter                                        ! boost

      case default bias_a                                            ! otherwise

         sort_a = 0                                                  ! ignore

      end select bias_a                                              ! which character

! ----------------------------------------------------------------------

      int_b = iachar( string_b( len_i: len_i) )                      ! get ascii character code

      bias_b: select case( int_b)                                    ! which character

      case( ascii_0: ascii_9) bias_b                                 ! digit

         sort_b = bias_digit                                         ! boost

      case( ascii_uca: ascii_ucz) bias_b                             ! upper case

         sort_b = bias_letter + change_case                          ! boost

      case( ascii_lca: ascii_lcz) bias_b                             ! lower case

         sort_b = bias_letter                                        ! boost

      case default bias_b                                            ! otherwise

         sort_b = 0                                                  ! ignore

      end select bias_b                                              ! which character

! ----------------------------------------------------------------------

      comp_ge: if( int_a + sort_a > int_b + sort_b )then

         ascii_lexge = .true.

         return                                                      ! .lexge.

      elseif( int_a + sort_a < int_b + sort_b )then comp_ge

         ascii_lexge = .false.

         return                                                      ! .lexge.

      endif comp_ge

   enddo comp_ab

! ----------------------------------------------------------------------

   ascii_lexge = len( string_a) >= len( string_b)                    ! lengths differ?

! ----------------------------------------------------------------------

return                                                               ! .lexge.

! **********************************************************************

!  ascii_lexge()

end function ascii_lexge

! **********************************************************************

!  ascii_lexgt()

elemental logical function ascii_lexgt( string_a, string_b)

!  ascii_lexgt() interface

character( len= *), intent( in) :: string_a, string_b

! **********************************************************************

!  ascii_lexgt() local

   integer( kind= int_k) :: int_a                                    ! from string_a

   integer( kind= int_k) :: int_b                                    ! from string_b

   integer( kind= int_k) :: sort_a                                   ! bias applied to a

   integer( kind= int_k) :: sort_b                                   ! bias applied to b

   integer( kind= int_k) :: len_i                                    ! loop variable

! **********************************************************************

!  ascii_lexgt()

continue                                                             ! .lexgt.

! ----------------------------------------------------------------------

!  loop thru each character

   comp_ab: do len_i = 1, min( len( string_a), len( string_b))

! ----------------------------------------------------------------------

      int_a = iachar( string_a( len_i: len_i) )                      ! get ascii character code

      bias_a: select case( int_a)                                    ! which character

      case( ascii_0: ascii_9) bias_a                                 ! digit

         sort_a = bias_digit                                         ! boost

      case( ascii_uca: ascii_ucz) bias_a                             ! upper case

         sort_a = bias_letter + change_case                          ! boost

      case( ascii_lca: ascii_lcz) bias_a                             ! lower case

         sort_a = bias_letter                                        ! boost

      case default bias_a                                            ! otherwise

         sort_a = 0                                                  ! ignore

      end select bias_a                                              ! which character

! ----------------------------------------------------------------------

      int_b = iachar( string_b( len_i: len_i) )                      ! get ascii character code

      bias_b: select case( int_b)                                    ! which character

      case( ascii_0: ascii_9) bias_b                                 ! digit

         sort_b = bias_digit                                         ! boost

      case( ascii_uca: ascii_ucz) bias_b                             ! upper case

         sort_b = bias_letter + change_case                          ! boost

      case( ascii_lca: ascii_lcz) bias_b                             ! lower case

         sort_b = bias_letter                                        ! boost

      case default bias_b                                            ! otherwise

         sort_b = 0                                                  ! ignore

      end select bias_b                                              ! which character

! ----------------------------------------------------------------------

      comp_gt: if( int_a + sort_a > int_b + sort_b )then

         ascii_lexgt = .true.

         return                                                      ! .lexgt.

      endif comp_gt

   enddo comp_ab

! ----------------------------------------------------------------------

   ascii_lexgt = len( string_a) > len( string_b)                     ! lengths differ?

! ----------------------------------------------------------------------

return                                                               ! .lexgt.

! **********************************************************************

!  ascii_lexgt()

end function ascii_lexgt

! **********************************************************************

!  isascii(): true if input is ascii character

! **********************************************************************

!  byte_isascii()

elemental logical function byte_isascii( bch)

!  byte_isascii() interface

integer(kind=byte_k), intent( in) :: bch

! **********************************************************************

!  byte_isascii() text

continue                                                             ! isascii()

! ----------------------------------------------------------------------

   ch_range: select case( bch)                                       ! which character

   case( ascii_nul: ascii_del) ch_range                              ! seven bit ascii

      byte_isascii = .true.                                          ! found it

   case default ch_range                                             ! otherwise

      byte_isascii = .false.                                         ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! isascii()

! **********************************************************************

!  byte_isascii()

end function byte_isascii

! **********************************************************************

!  short_isascii()

elemental logical function short_isascii( sch)

!  short_isascii() interface

integer( kind= short_k), intent( in) :: sch

! **********************************************************************

!  short_isascii() text

continue                                                             ! isascii()

! ----------------------------------------------------------------------

   ch_range: select case( sch)                                       ! which character

   case( ascii_nul: ascii_del) ch_range                              ! seven bit ascii

      short_isascii = .true.                                         ! found it

   case default ch_range                                             ! otherwise

      short_isascii = .false.                                        ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! isascii()

!  short_isascii()

end function short_isascii

! **********************************************************************

!  int_isascii()

elemental logical function int_isascii( ich)

! int_isascii() interface

integer( kind= int_k), intent( in) :: ich

! **********************************************************************

!  int_isascii() text

continue                                                             ! isascii()

! ----------------------------------------------------------------------

   ch_range: select case( ich)                                       ! which character

   case( ascii_nul: ascii_del) ch_range                              ! seven bit ascii

      int_isascii = .true.                                           ! found it

   case default ch_range                                             ! otherwise

      int_isascii = .false.                                          ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! isascii()

! **********************************************************************

!  int_isascii()

end function int_isascii

! **********************************************************************

!  long_isascii()

elemental logical function long_isascii( lch)

! long_isascii() interface

integer( kind= long_k), intent( in) :: lch

! **********************************************************************

!  long_isascii() text

continue                                                             ! isascii()

! ----------------------------------------------------------------------

   ch_range: select case( lch)                                       ! which character

   case( ascii_nul: ascii_del) ch_range                              ! seven bit ascii

      long_isascii = .true.                                          ! found it

   case default ch_range                                             ! otherwise

      long_isascii = .false.                                         ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! isascii()

! **********************************************************************

!  long_isascii()

end function long_isascii

! **********************************************************************

!  character_isascii()

elemental logical function character_isascii( ch)

!  character_isascii() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  character_isascii() local

   integer :: ich

! **********************************************************************

!  character_isascii() text

continue                                                             ! isascii()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case( ich)                                       ! which character

   case( ascii_nul: ascii_del) ch_range                              ! seven bit ascii

      character_isascii = .true.                                     ! found it

   case default ch_range                                             ! otherwise

      character_isascii = .false.                                    ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! isascii()

! **********************************************************************

!  character_isascii()

end function character_isascii

! **********************************************************************

!  isupper(): true if input ascii character is uppercase

elemental logical function isupper( ch)

!  isupper() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  isupper() local

   integer :: ich

! **********************************************************************

!  isupper() text

continue                                                             ! isupper()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case( ich)                                       ! which character

   case( ascii_uca: ascii_ucz) ch_range                              ! upper case

      isupper = .true.                                               ! found it

   case default ch_range                                             ! otherwise

      isupper = .false.                                              ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! isupper()

! **********************************************************************

!  isupper()

end function isupper

! **********************************************************************

!  islower(): true if input ascii character is lowercase

elemental logical function islower( ch)

!  islower() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  islower() local

   integer :: ich

! **********************************************************************

!  islower() text

continue                                                             ! islower()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case( ich)                                       ! which character

   case( ascii_lca: ascii_lcz) ch_range                              ! lower case

      islower = .true.                                               ! found it

   case default ch_range                                             ! otherwise

      islower = .false.                                              ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! islower()

! **********************************************************************

!  islower()

end function islower

! **********************************************************************

!  isprint(): true if input ascii character is printable

elemental logical function isprint( ch)

!  isprint() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  isprint() local

   integer :: ich

! **********************************************************************

!  isprint() text

continue                                                             ! isprint()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case( ich)                                       ! which character

   case( ascii_sp: ascii_tld) ch_range                               ! graphic character

      isprint = .true.                                               ! found it

   case default ch_range                                             ! otherwise

      isprint = .false.                                              ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! isprint()

! **********************************************************************

!  isprint()

end function isprint

! **********************************************************************

!  isgraph(): true if input ascii character is printable

elemental logical function isgraph( ch)

!  isgraph() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  isgraph() local

   integer :: ich

! **********************************************************************

!  isgraph() text

continue                                                             ! isgraph()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case( ich)                                       ! which character

   case( ascii_ep: ascii_tld) ch_range                               ! graphic character

      isgraph = .true.                                               ! found it

   case default ch_range                                             ! otherwise

      isgraph = .false.                                              ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! isgraph()

! **********************************************************************

!  isgraph()

end function isgraph

! **********************************************************************

!  isalpha(): true if input ascii character is alphabetic

elemental logical function isalpha( ch)

!  isalpha() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  isalpha() local

   integer :: ich

! **********************************************************************

!  isalpha() text

continue                                                             ! isalpha()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case( ich)                                       ! which character

!  upper case letter or lower case letter

   case( ascii_uca: ascii_ucz, ascii_lca: ascii_lcz) ch_range

      isalpha = .true.                                               ! found it

   case default ch_range                                             ! otherwise

      isalpha = .false.                                              ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! isalpha()

! **********************************************************************

!  isalpha()

end function isalpha

! **********************************************************************

!  isalnum(): true if input ascii character is alphanumeric

elemental logical function isalnum( ch)

!  isalnum() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  isalnum() local

   integer :: ich

! **********************************************************************

!  isalnum() text

continue                                                             ! isalnum()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case( ich)                                       ! which character

!  digit or upper case letter or lower case letter

   case( ascii_0: ascii_9, ascii_uca: ascii_ucz, ascii_lca: ascii_lcz) ch_range

      isalnum = .true.                                               ! found it

   case default ch_range                                             ! otherwise

      isalnum = .false.                                              ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! isalnum()

! **********************************************************************

!  isalnum()

end function isalnum

! **********************************************************************

!  isdigit(): true if input ascii character is numeric

elemental logical function isdigit( ch)

!  isdigit() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  isdigit() local

   integer :: ich

! **********************************************************************

!  isdigit() text

continue                                                             ! isdigit()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case( ich)                                       ! which character

   case( ascii_0: ascii_9) ch_range                                  ! digit

      isdigit = .true.                                               ! found it

   case default ch_range                                             ! otherwise

      isdigit = .false.                                              ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! isdigit()

! **********************************************************************

!  isdigit()

end function isdigit

! **********************************************************************

!  isxdigit(): true if input ascii character is numeric

elemental logical function isxdigit( ch)

!  isxdigit() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  isxdigit() local

   integer :: ich

! **********************************************************************

!  isxdigit() text

continue                                                             ! isxdigit()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case( ich)                                       ! which character

   case( ascii_0: ascii_9, ascii_uca: ascii_ucz, ascii_lca: ascii_lcz) ch_range

      isxdigit = .true.                                               ! found it

   case default ch_range                                             ! otherwise

      isxdigit = .false.                                              ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! isxdigit()

! **********************************************************************

!  isxdigit()

end function isxdigit

! **********************************************************************

!  ispunct(): true if input ascii character is punctuation

elemental logical function ispunct( ch)

!  ispunct() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  ispunct() local

   integer :: ich

! **********************************************************************

!  ispunct() text

continue                                                             ! ispunct()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case ( ich)                                      ! which character

!  punctuation character

   case ( ascii_ep: ascii_sl, ascii_cl: ascii_ats, ascii_lb: ascii_gra, ascii_lbr: ascii_tld) ch_range

         ispunct = .true.                                            ! found it

   case default ch_range                                             ! otherwise

         ispunct = .false.                                           ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! ispunct()

! **********************************************************************

!  ispunct()

end function ispunct

! **********************************************************************

!  isspace(): true if input ascii character is numeric

elemental logical function isspace( ch)

!  isspace() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  isspace() local

   integer :: ich

! **********************************************************************

!  isspace() text

continue                                                             ! isspace()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case( ich)                                       ! which character

   case( ascii_ht, ascii_lf, ascii_sp) ch_range                      ! white space character

      isspace = .true.                                               ! found it

   case default ch_range                                             ! otherwise

      isspace = .false.                                              ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! isspace()

! **********************************************************************

!  isspace()

end function isspace

! **********************************************************************

!  iscntl(): true if input ascii character is control

elemental logical function iscntl( ch)

!  iscntl() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  iscntl() local

   integer :: ich

! **********************************************************************

!  iscntl() text

continue                                                             ! iscntl()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case( ich)                                       ! which character

   case( ascii_soh: ascii_us) ch_range                               ! control character

      iscntl = .true.                                                ! found it

   case default ch_range                                             ! otherwise

      iscntl = .false.                                               ! found it not

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! iscntl()

! **********************************************************************

!  iscntl()

end function iscntl

! **********************************************************************

!  toupper(): if lowercase, change to uppercase

elemental character( len= 1) function toupper( ch)

!  toupper() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  toupper() local

   integer :: ich

! **********************************************************************

!  toupper() text

continue                                                             ! toupper()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case( ich)                                       ! which character

   case( ascii_lca: ascii_lcz) ch_range                              ! lower case

      toupper = achar( ich - change_case)                            ! change to upper case

   case default ch_range                                             ! otherwise

      toupper = achar( ich)                                          ! ignore

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! toupper()

! **********************************************************************

!  toupper()

end function toupper

! **********************************************************************

!  tolower(): if uppercase, change to lowercase

elemental character( len= 1) function tolower( ch)

!  tolower() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  tolower() local

   integer :: ich

! **********************************************************************

!  tolower() text

continue                                                             ! tolower()

! ----------------------------------------------------------------------

   ich = iachar( ch)                                                 ! get ascii character code

   ch_range: select case( ich)                                       ! which character

   case( ascii_uca: ascii_ucz) ch_range                              ! upper case

      tolower = achar( ich + change_case)                            ! change to lower case

   case default ch_range                                             ! otherwise

      tolower = achar( ich)                                          ! ignore

   end select ch_range                                               ! which character

! ----------------------------------------------------------------------

return                                                               ! tolower()

! **********************************************************************

!  tolower()

end function tolower

! **********************************************************************

!     isevenp()/isoddp()/ismarkp()/isclrp() characterize parity

! **********************************************************************

!  isevenp(): true if character has even parity

elemental logical function isevenp( ch)

! isevenp() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  isevenp() local

   integer :: ich

! **********************************************************************

!  isevenp() text

continue                                                             ! isevenp()

! ----------------------------------------------------------------------

   ich = iand( ichar( ch), ascii_bits)                               ! all but high bit

   parity: if( even_parity( ich) )then                               ! even parity

      isevenp = iand( ich, parity_bit) == 0                          ! high bit clear

   else parity                                                       ! otherwise

      isevenp = iand( ich, parity_bit) == parity_bit                 !  parity bit set

   endif parity                                                      ! even parity

! ----------------------------------------------------------------------

return                                                               ! isevenp()

! **********************************************************************

!  isevenp()

end function isevenp

! **********************************************************************

!  isoddp(): true if character has even parity

elemental logical function isoddp( ch)

!  isoddp() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  isoddp() local

   integer( kind= int_k) :: ich

! **********************************************************************

!  isoddp() text

continue                                                             ! isoddp()

! ----------------------------------------------------------------------

   ich = iand( ichar( ch), ascii_bits)                               ! all but high bit

   parity: if( even_parity( ich) )then                               ! even parity

      isoddp = iand( ich, parity_bit) == parity_bit                  !  parity bit set

   else parity                                                       ! otherwise

      isoddp = iand( ich, parity_bit) == 0                           ! parity bit clear

   endif parity                                                      ! even parity

! ----------------------------------------------------------------------

return                                                               ! isoddp()

! **********************************************************************

!  isoddp()

end function isoddp

! **********************************************************************

!  ismarkp(): true if character has even parity

elemental logical function ismarkp( ch)

character( len= 1), intent( in) :: ch

! **********************************************************************

!  ismarkp()

continue                                                             ! ismarkp()

! ----------------------------------------------------------------------

   ismarkp = iand( ichar( ch), parity_bit) == parity_bit             ! parity bit set

! ----------------------------------------------------------------------

return                                                               ! ismarkp()

! **********************************************************************

!  ismarkp()

end function ismarkp

! **********************************************************************

!  isclrp(): true if character has even parity

elemental logical function isclrp( ch)

character( len= 1), intent( in) :: ch

! **********************************************************************

!  isclrp()

continue                                                             ! isclrp()

! ----------------------------------------------------------------------

   isclrp = iand( ichar( ch), parity_bit) == 0                       ! parity bit clear

! ----------------------------------------------------------------------

return                                                               ! isclrp()

! **********************************************************************

!  isclrp()

end function isclrp

! **********************************************************************

!     toevenp()/tooddp()/tomarkp()/toclrp() set parity

! **********************************************************************

!  toevenp(): true if character has even parity

elemental character( len= 1) function toevenp( ch)

!  toevenp() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  toevenp() local

   integer( kind= int_k) :: ich

! **********************************************************************

!  toevenp()

continue                                                             ! toevenp()

! ----------------------------------------------------------------------

   ich = iand( ichar( ch), ascii_bits)                               ! all but high bit

   parity: if( even_parity( ich) )then                               ! even parity

      toevenp = achar( ich)                                          ! already ok

   else parity                                                       ! otherwise

      toevenp = char( ior( ich, parity_bit))                         ! set parity bit

   endif parity                                                      ! even parity

! ----------------------------------------------------------------------

return                                                               ! toevenp()

! **********************************************************************

!  toevenp()

end function toevenp

! **********************************************************************

!  tooddp(): true if character has even parity

elemental character( len= 1) function tooddp( ch)

!  tooddp() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  tooddp() local

   integer( kind= int_k) :: ich

! **********************************************************************

!  tooddp() text

continue                                                             ! tooddp()

! ----------------------------------------------------------------------

   ich = iand( iachar( ch), ascii_bits)                              ! all but high bit

   parity: if( even_parity( ich) )then                               ! even parity

      tooddp = char( ior( ich, parity_bit))                          ! set parity bit

   else parity                                                       ! otherwise

      tooddp = char( ich)                                            ! ignore

   endif parity                                                      ! even parity

! ----------------------------------------------------------------------

return                                                               ! tooddp()

! **********************************************************************

!  tooddp()

end function tooddp

! **********************************************************************

!  tomarkp(): true if character has even parity

elemental character( len= 1) function tomarkp( ch)

!  tomarkp() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  tomarkp() text

continue                                                             ! tomarkp()

! ----------------------------------------------------------------------

   tomarkp = char( ior( ichar( ch), parity_bit))                     ! set parity bit

! ----------------------------------------------------------------------

return                                                               ! tomarkp()

! **********************************************************************

!  tomarkp()

end function tomarkp

! **********************************************************************

!  toclrp(): true if character has even parity

elemental character( len= 1) function toclrp( ch)

!  toclrp() interface

character( len= 1), intent( in) :: ch

! **********************************************************************

!  toclrp() text

continue                                                             ! toclrp()

! ----------------------------------------------------------------------

   toclrp = char( iand( ichar( ch), ascii_bits) )                    ! clear high bit

! ----------------------------------------------------------------------

return                                                               ! toclrp()

! **********************************************************************

!  toclrp()

end function toclrp

! **********************************************************************

!  rbn()/rnb(): convert trailing blanks to/from nulls

! **********************************************************************

!  rbn(): trailing blanks to nulls

elemental subroutine rbn( string)

!  rbn() interface

character( len= *), intent( inout) :: string

! **********************************************************************

!  rnb() local

   integer :: len_i

! **********************************************************************

!  rbn() text

continue                                                             ! rbn()

! ----------------------------------------------------------------------

   trail: do len_i = len_trim( string)+1, len( string)               ! trailing blanks if any

      string( len_i: len_i) = achar( ascii_nul)                      ! to nulls

   enddo trail                                                       ! trailing blanks if any

! ----------------------------------------------------------------------

return                                                               ! rbn()

! **********************************************************************

!  rbn()

end subroutine rbn

! **********************************************************************

!  rnb(): trailing nulls to blanks

elemental subroutine rnb( string)

!  rnb() interface

character( len= *), intent( inout) :: string

! **********************************************************************

!  rnb() local

   integer :: len_i

! **********************************************************************

!  rnb() text

continue                                                             ! rnb()

! ----------------------------------------------------------------------

   trail: do len_i = len( string), 1, -1                             ! last to first

      null_quit: if( string( len_i: len_i) == achar( ascii_nul) )then

         string( len_i: len_i) = achar( ascii_sp)                    ! null to blank

      else null_quit                                                 ! null or quit

         exit                                                        ! quit at first non-null

      endif null_quit                                                ! null or quit

   enddo trail                                                       ! last to first

! ----------------------------------------------------------------------

return                                                               ! rnb()

! **********************************************************************

!  rnb()

end subroutine rnb

! **********************************************************************

!  bcomp(): blank compress ignoring quoted strings

! **********************************************************************

!  bcomp()

elemental subroutine bcomp( outstr, instr)

!  bcomp() interface

character( len= *), intent( in) :: instr

character( len= *), intent( out) :: outstr

! **********************************************************************

!  bcomp() local

   integer :: len_to_do

   integer :: in_ptr

   integer :: out_ptr

   integer :: index_ptr

! **********************************************************************

!  bcomp() text

continue

! ----------------------------------------------------------------------

!  do each character in instr

   len_to_do = len_trim( instr)

!  initialize outstr and character pointers

   outstr = null_string

   in_ptr = 1
   out_ptr = 1

! ----------------------------------------------------------------------

!  loop thru each character

   each_char: do while( in_ptr <= len_to_do)                         ! do instr

      non_blanks: if( instr( in_ptr: in_ptr) /= blank )then
 
! ----------------------------------------------------------------------

!  check non blank characters for single quote

         sq_string: if( instr( in_ptr: in_ptr) == single_quote )then

            index_ptr = index( instr( in_ptr+1: ), single_quote)

            no_sq: if( index_ptr > substring_not_found )then

!  if single quoted string, copy as is

               outstr( out_ptr: out_ptr+index_ptr) = instr( in_ptr: in_ptr+index_ptr)

               in_ptr = in_ptr + index_ptr                           ! advance in pointer
               out_ptr = out_ptr + index_ptr                         ! advance out pointer

               cycle each_char                                       ! get next character

            endif no_sq

         endif sq_string

! ----------------------------------------------------------------------

!  check non blank characters for double quote

         dq_string: if( instr( in_ptr: in_ptr) == double_quote )then

            index_ptr = index( instr( in_ptr+1: ), double_quote)

            no_dq: if( index_ptr > substring_not_found )then

!  if double quoted string, copy as is

               outstr( out_ptr: out_ptr+index_ptr) = instr( in_ptr: in_ptr+index_ptr)

               in_ptr = in_ptr + index_ptr                           ! advance in pointer
               out_ptr = out_ptr + index_ptr                         ! advance out pointer

               cycle each_char                                       ! get next character

            endif no_dq

         endif dq_string

! ----------------------------------------------------------------------

!  copy non blank character

         outstr( out_ptr: out_ptr) = instr( in_ptr: in_ptr)

         out_ptr = out_ptr + 1                                       ! advance out pointer

         if_exit: if( out_ptr > len( outstr) )then

            exit each_char                                           ! truncate or blank fill

         endif if_exit

      endif non_blanks

! ----------------------------------------------------------------------

!  increment

      in_ptr = in_ptr + 1                                            ! advance in pointer

   enddo each_char                                                   ! do instr

!  truncate or blank fill

   to_end: if( out_ptr <= len( outstr) )then

      outstr( out_ptr: ) = blank

   endif to_end

! ----------------------------------------------------------------------

return                                                               ! bcomp()

! **********************************************************************

!  bcomp()

end subroutine bcomp

! **********************************************************************

!  bclc(): blank compress to lower case ignoring quoted strings

! **********************************************************************

!  bclc()

elemental subroutine bclc( outstr, instr)

!  bclc() interface

character( len= *), intent( in) :: instr

character( len= *), intent( out) :: outstr

! **********************************************************************

!  bclc() local

   integer :: len_to_do

   integer :: in_ptr

   integer :: out_ptr

   integer :: index_ptr

   integer :: ich

! **********************************************************************

!  bclc() text

continue

! ----------------------------------------------------------------------

!  do each character in instr

   len_to_do = len_trim( instr)

!  initialize outstr and character pointers

   outstr = null_string

   in_ptr = 1
   out_ptr = 1

! ----------------------------------------------------------------------

!  loop thru each character

   each_char: do while( in_ptr <= len_to_do)

      non_blanks: if( instr( in_ptr: in_ptr) /= blank )then
 
! ----------------------------------------------------------------------

!  check non blank characters for single quote

         sq_string: if( instr( in_ptr: in_ptr) == single_quote )then

            index_ptr = index( instr( in_ptr+1: ), single_quote)

            no_sq: if( index_ptr > substring_not_found )then

!  if single quoted string, copy as is

               outstr( out_ptr: out_ptr+index_ptr) = instr( in_ptr: in_ptr+index_ptr)

               in_ptr = in_ptr + index_ptr
               out_ptr = out_ptr + index_ptr

               cycle each_char

            endif no_sq

         endif sq_string

! ----------------------------------------------------------------------

!  check non blank characters for double quote

         dq_string: if( instr( in_ptr: in_ptr) == double_quote )then

            index_ptr = index( instr( in_ptr+1: ), double_quote)

            no_dq: if( index_ptr > substring_not_found )then

!  if double quoted string, copy as is

               outstr( out_ptr: out_ptr+index_ptr) = instr( in_ptr: in_ptr+index_ptr)

               in_ptr = in_ptr + index_ptr
               out_ptr = out_ptr + index_ptr

               cycle each_char

            endif no_dq

         endif dq_string

! ----------------------------------------------------------------------

!  convert to lower case

         ich = iachar( instr( in_ptr: in_ptr) )                      ! get ascii character code

         select case( ich)                                           ! which character

         case( ascii_uca: ascii_ucz)                                 ! upper case

            outstr( out_ptr: out_ptr) = achar( ich + change_case)

         case default                                                ! otherwise

            outstr( out_ptr: out_ptr) = achar( ich)

         end select                                                  ! which character

! ----------------------------------------------------------------------

!  increment

         out_ptr = out_ptr + 1

         if_exit: if( out_ptr > len( outstr) )then

            exit each_char

         endif if_exit

      endif non_blanks

! ----------------------------------------------------------------------

!  increment

      in_ptr = in_ptr + 1

   enddo each_char                                                   ! do instr

   to_end: if( out_ptr <= len( outstr) )then

      outstr( out_ptr: ) = blank

   endif to_end

! ----------------------------------------------------------------------

return                                                               ! bclc()

! **********************************************************************

!  bclc()

end subroutine bclc

! **********************************************************************

!  bcuc(): blank compress to lower case ignoring quoted strings

! **********************************************************************

!  bcuc()

elemental subroutine bcuc( outstr, instr)

character( len= *), intent( in) :: instr

character( len= *), intent( out) :: outstr

! **********************************************************************

!  bcuc() local

   integer :: len_to_do

   integer :: in_ptr

   integer :: out_ptr

   integer :: index_ptr

   integer :: ich

! **********************************************************************

continue

! ----------------------------------------------------------------------

!  do each character in instr

   len_to_do = len_trim( instr)

!  initialize outstr and character pointers

   outstr = null_string

   in_ptr = 1
   out_ptr = 1

! ----------------------------------------------------------------------

!  loop thru each character

   each_char: do while( in_ptr <= len_to_do)

      non_blanks: if( instr( in_ptr: in_ptr) /= blank )then
 
! ----------------------------------------------------------------------

!  check non blank characters for single quote

         sq_string: if( instr( in_ptr: in_ptr) == single_quote )then

            index_ptr = index( instr( in_ptr+1: ), single_quote)

            no_sq: if( index_ptr > substring_not_found )then

!  if single quoted string, copy as is

               outstr( out_ptr: out_ptr+index_ptr) = instr( in_ptr: in_ptr+index_ptr)

               in_ptr = in_ptr + index_ptr
               out_ptr = out_ptr + index_ptr

               cycle each_char

            endif no_sq

         endif sq_string

! ----------------------------------------------------------------------

!  check non blank characters for double quote

         dq_string: if( instr( in_ptr: in_ptr) == double_quote )then

            index_ptr = index( instr( in_ptr+1: ), double_quote)

            no_dq: if( index_ptr > substring_not_found )then

!  if double quoted string, copy as is

               outstr( out_ptr: out_ptr+index_ptr) = instr( in_ptr: in_ptr+index_ptr)

               in_ptr = in_ptr + index_ptr
               out_ptr = out_ptr + index_ptr

               cycle each_char

            endif no_dq

         endif dq_string

! ----------------------------------------------------------------------

!  convert to upper case

         ich = iachar( instr( in_ptr: in_ptr) )                      ! get ascii character code

         select case( ich)                                           ! which character

         case( ascii_lca: ascii_lcz)                                 ! lower case

            outstr( out_ptr: out_ptr) = achar( ich - change_case)

         case default                                                ! otherwise

            outstr( out_ptr: out_ptr) = achar( ich)

         end select                                                  ! which character

! ----------------------------------------------------------------------

!  increment

         out_ptr = out_ptr + 1

         if_exit: if( out_ptr > len( outstr) )then

            exit each_char                                           ! truncate or blank fill

         endif if_exit

      endif non_blanks

! ----------------------------------------------------------------------

!  increment

      in_ptr = in_ptr + 1

   enddo each_char                                                   ! do instr

! truncate or blank fill

   to_end: if( out_ptr <= len( outstr) )then

      outstr( out_ptr: ) = blank

   endif to_end

! ----------------------------------------------------------------------

return                                                               ! bcuc()

! **********************************************************************

!  bcuc()

end subroutine bcuc

! **********************************************************************

!  tr(): translate string in place via translation table

! **********************************************************************

!  array_tr(): tr() for array of char*1(i)

subroutine array_tr( string, offset, nchars, ttable)

character( len= 1), dimension( *), intent( inout) :: string

integer( kind= int_k), intent( in) :: offset, nchars

character( len= 1), dimension( 0: 255), intent( in) :: ttable

! **********************************************************************

!  array_tr() local

   integer :: len_i

! **********************************************************************

!  array_tr()

continue                                                             ! tr()

! ----------------------------------------------------------------------

   each_char: do len_i = offset, offset + nchars                     ! loop thru characters

      string( len_i) = ttable( ichar( string( len_i)) )              ! translate

   enddo each_char                                                   ! loop thru characters

! ----------------------------------------------------------------------

return                                                               ! tr()

! **********************************************************************

!  array_tr()

end subroutine array_tr

! **********************************************************************

!  substring_tr(): tr() for substring char( i: j)

subroutine substring_tr( string, offset, nchars, ttable)

character( len= *), intent( inout) :: string

integer( kind= int_k), intent( in) :: offset, nchars

character( len= 1), dimension( 0: 255), intent( in) :: ttable

! **********************************************************************

!  substring_tr() local

   integer :: len_i

! **********************************************************************

!  substring_tr()

continue                                                             ! tr()

! ----------------------------------------------------------------------

   each_char: do len_i = offset, offset + nchars                     ! loop thru characters

      string( len_i: len_i) = ttable( ichar( string( len_i: len_i)) )

   enddo each_char                                                   ! loop thru characters

! ----------------------------------------------------------------------

return                                                               ! tr()

! **********************************************************************

!  substring_tr()

end subroutine substring_tr

! **********************************************************************

!  compare_string() location of first difference between two strings

elemental integer function compare_string( string_a, string_b)

! **********************************************************************

!  compare_string() interface

character( len= *), intent( in) :: string_a                          ! the first string

character( len= *), intent( in) :: string_b                          ! the second string

! **********************************************************************

!  compare_string() local

! ----------------------------------------------------------------------

!  loop thru the strings

   integer :: i                                                      ! loop index

!  loop length

   integer :: loop_len                                               ! length of shorter string

! **********************************************************************

!  compare_string() text

continue                                                             ! compare_string()

!  initialize

   compare_string = 0                                                ! true if match found

   loop_len = min( len( string_a), len( string_b))

! ----------------------------------------------------------------------

!  scan thru the quoted string

   scan_string: do i = 1, loop_len

      mismatch: if( string_a( i: i) /= string_b( i: i) )then

         compare_string = i

         exit scan_string

      endif mismatch

   enddo scan_string                                                 ! loop thru quoted string

! ----------------------------------------------------------------------

return                                                               ! compare_string()

!  compare_string()

! **********************************************************************

end function compare_string

! **********************************************************************

!  replace_substring() replace all occurances of a target substring by a replacement

elemental subroutine replace_substring( string, substring, repstring)

! **********************************************************************

!  replace_substring() interface

character( len= *), intent( inout) :: string                         ! the string containing the substrings

character( len= *), intent( in) :: substring                         ! the target substring

character( len= *), intent( in) :: repstring                         ! the replacement string

! **********************************************************************

!  replace_substring() local

! ----------------------------------------------------------------------

!  search the string

   integer :: i                                                      ! loop index

! **********************************************************************

!  replace_substring() text

continue                                                             ! replace_substring()

!  initialize

   i = index( string, substring)                                     ! try to find substring

! ----------------------------------------------------------------------

!  scan thru the quoted string

   scan_string: do while( i > substring_not_found)                   ! find all substrings

      string = string( 1: i - 1) // repstring // string( i + len( substring): )

      i = index( string, substring)                                  ! find next substring, if any

   enddo scan_string                                                 ! loop thru quoted string

! ----------------------------------------------------------------------

return                                                               ! replace_substring()

!  replace_substring()

! **********************************************************************

end subroutine replace_substring

! **********************************************************************

!  string_cat() copy string_b into string_a following all nonblanks

elemental subroutine string_cat( string_a, string_b)

! **********************************************************************

!  string_cat() interface

character( len= *), intent( inout) :: string_a                       ! the string containing the substrings

character( len= *), intent( in) :: string_b                          ! the target substring

! **********************************************************************

!  Note that call string_cat( foo // ' ', bar) is the same as
!  call string_cat( foo, bar).  To place a blank between the two strings,
!  use call string_cat( foo, ' ' // bar).

! **********************************************************************

!  string_cat() text

continue                                                             ! string_cat()

! ----------------------------------------------------------------------

!  concatenate string_b after len_trim( string_a)

   string_a = trim( string_a) // string_b

! ----------------------------------------------------------------------

return                                                               ! string_cat()

!  string_cat()

! **********************************************************************

end subroutine string_cat

! **********************************************************************

!  matching_pair() true if successfully found matching (), [], or {}

logical function matching_pair( string, start, found)

! **********************************************************************

!  matching_pair() interface

character( len= *), intent( in) :: string                            ! the string to be searched

integer, intent( in) :: start                                        ! location to start the search

integer, intent( out) :: found                                       ! location of matching paren

! **********************************************************************

!  matching_pair() constants

   character( len= *), parameter :: open_paren = '('
   character( len= *), parameter :: close_paren = ')'

   character( len= *), parameter :: open_brace = '['
   character( len= *), parameter :: close_brace = ']'

   character( len= *), parameter :: open_curly = '{'
   character( len= *), parameter :: close_curly = '}'

! **********************************************************************

!  matching_pair() local

   integer :: step                                                   ! local back

! ----------------------------------------------------------------------

!  counters and pointers

   integer :: level                                                  ! level of parenthesis encountered

   integer :: iptr                                                   ! point to character in string

   integer :: end_point                                              ! either end of string

   character( len= 1) :: another                                     ! starting character

   character( len= 1) :: matching                                    ! target of search

! **********************************************************************

!  matching_pair() text

continue                                                             ! matching_pair()

!  not correct yet

   matching_pair = .false.                                           ! true if match found

! ----------------------------------------------------------------------

!  initialize before scanning string

   match_direction: select case( string( start: start) )             ! which search?

!  seek forward matching close parenthesis

   case( open_paren) match_direction                                 ! got (, find )

      step = 1                                                       ! go forward

      end_point = len_trim( string)                                  ! to end

      another = open_paren                                           ! got (

      matching = close_paren                                         ! find )

!  seek backward matching open parenthesis

   case( close_paren) match_direction                                ! got ), find (

      step = -1                                                      ! go backward

      end_point = 1                                                  ! to beginning

      another = close_paren                                          ! got )

      matching = open_paren                                          ! find (

!  seek forward matching close brace

   case( open_brace) match_direction                                 ! got [, find ]

      step = 1                                                       ! go forward

      end_point = len_trim( string)                                  ! to end

      another = open_brace                                           ! got [

      matching = close_brace                                         ! find ]

!  seek backward matching open brace

   case( close_brace) match_direction                                ! got ], find [

      step = -1                                                      ! go backward

      end_point = 1                                                  ! to beginning

      another = close_brace                                          ! got ]

      matching = open_brace                                          ! find [

!  seek forward matching close curly brace

   case( open_curly) match_direction                                 ! got {, find }

      step = 1                                                       ! go forward

      end_point = len_trim( string)                                  ! to end

      another = open_curly                                           ! got {

      matching = close_curly                                         ! find }

!  seek backward matching open curly brace

   case( close_curly) match_direction                                ! got }, find {

      step = -1                                                      ! go backward

      end_point = 1                                                  ! to beginning

      another = close_curly                                          ! got }

      matching = open_curly                                          ! find {

!  nothing to seek

   case default match_direction                                      ! not right

      return                                                         ! quit

   end select match_direction                                        ! which search?

!  matching character has the same level

   level = 0                                                         ! count nesting level

! ----------------------------------------------------------------------

!  scan through string

   search: do iptr = start + step, end_point, step                   ! look thru string

!  count nesting levels

      levels: if( string( iptr: iptr) == another )then               ! nested

         level = level + 1                                           ! increment counter

      elseif( string( iptr: iptr) == matching )then levels           ! end

!  if matching character has the same level

         found_match: if( level == 0 )then                           ! match?

!  found target

            found = iptr                                             ! mark

            matching_pair = .true.                                   ! eureka

            exit search                                              ! quit

         endif found_match                                           ! match?

!  count level

         level = level - 1                                           ! decrement counter

      endif levels                                                   ! nesting

!  continue scan

   enddo search                                                      ! look thru string

! ----------------------------------------------------------------------

return                                                               ! matching_pair()

!  matching_pair()

! **********************************************************************

end function matching_pair

! **********************************************************************

!  unquote_string() true if extracts string from between quotes

logical function unquote_string( quoted_str, unquoted_str, esc, quoted_len, unquoted_len)

! **********************************************************************

!  unquote_string() interface

character( len= *), intent( in) :: quoted_str                        ! the string to be unquoted

character( len= *), intent( out) :: unquoted_str                     ! the unquoted string

character( len= 1), optional, intent( in) :: esc                     ! escape character

integer, intent( out) :: quoted_len                                  ! length scanned in quoted string

integer, intent( out) :: unquoted_len                                ! length of unquoted string

! **********************************************************************

!  unquote_string() constants

   character( len= *), parameter :: single_quote = "'"

   character( len= *), parameter :: double_quote = '"'

! **********************************************************************

!  unquote_string() local

! ----------------------------------------------------------------------

!  either single quote or double quote

   character( len= 1) :: quote                                       ! whichever quote is to be used

! **********************************************************************

!  unquote_string() text

continue                                                             ! unquote_string()

!  not correct yet

   unquote_string = .false.                                          ! true if match found

   quote = ''                                                        ! no quote yet

! ----------------------------------------------------------------------

!  which quote is the first quote (if either)

   which_quote: select case( quoted_str( 1: 1) )                     ! which quote

! ----------------------------------------------------------------------

!  string delimited by single quote

   case( single_quote) which_quote                                   ! use single quote

      quote = single_quote                                           ! quote is single quote

      quoted_len = 2                                                 ! one character in quoted string

      unquoted_len = 0                                               ! no characters in unquoted string

! ----------------------------------------------------------------------

!  string delimited by double quote

   case( double_quote) which_quote                                   ! use double quote

      quote = double_quote                                           ! quote is double quote

      quoted_len = 2                                                 ! one character in quoted string

      unquoted_len = 0                                               ! no characters in unquoted string

! ----------------------------------------------------------------------

!  string delimited by neither quote- nothing to do

   case default which_quote                                          ! neither

      quoted_len = 0                                                 ! no quote

      unquoted_len = 0                                               ! no unquoted string

      unquoted_str = null_string                                     ! null output

      return                                                         ! quit

   end select which_quote

! ----------------------------------------------------------------------

!  is there an escape character?

   got_esc: if( present( esc) )then                                  ! process woth or without esc

! ----------------------------------------------------------------------

!  scan thru the quoted string

      scan_string_esc: do while( quoted_len <= len_trim( quoted_str) )

!  if find one matching quote

         next_char_esc: if( quoted_str( quoted_len: quoted_len) == quote )then

!  check next character

            quoted_len = quoted_len + 1                              ! next character in quoted string

!  check for a pair of quotes

            next_quote_esc: if( quoted_str( quoted_len: quoted_len) == quote )then

!  two consequetive quotes represents one quote

               unquoted_len = unquoted_len + 1                       ! next character in unquoted

               unquoted_str( unquoted_len: unquoted_len) = quoted_str( quoted_len: quoted_len)

!  quote followed by any other character

            else next_quote_esc                                      ! find two adjacent quotes

!  one quote is the matching end quote

               exit scan_string_esc                                  ! quit

            endif next_quote_esc                                     ! find two adjacent quotes

! ----------------------------------------------------------------------

!  character is the escape character

         elseif( quoted_str( quoted_len: quoted_len) == esc )then next_char_esc

!  get next character

            quoted_len = quoted_len + 1                              ! next character in quoted string

            unquoted_len = unquoted_len + 1                          ! next character in unquoted

            unquoted_str( unquoted_len: unquoted_len) = quoted_str( quoted_len: quoted_len)

! ----------------------------------------------------------------------

!  character is not a quote nor an escape

         else next_char_esc                                          ! find a quote

            unquoted_len = unquoted_len + 1                          ! next character in unquoted

            unquoted_str( unquoted_len: unquoted_len) = quoted_str( quoted_len: quoted_len)

         endif next_char_esc                                         ! find a quote

!  check next character

         quoted_len = quoted_len + 1                                 ! next character in quoted string

      enddo scan_string_esc                                          ! loop thru quoted string

! ----------------------------------------------------------------------

   else got_esc                                                      ! process woth or without esc

! ----------------------------------------------------------------------

!  scan thru the quoted string

      scan_string: do while( quoted_len <= len_trim( quoted_str) )

!  if find one matching quote

         next_char: if( quoted_str( quoted_len: quoted_len) == quote )then

!  check next character

            quoted_len = quoted_len + 1                              ! next character in quoted string

!  check for a pair of quotes

            next_quote: if( quoted_str( quoted_len: quoted_len) == quote )then

!  two consequetive quotes represents one quote

               unquoted_len = unquoted_len + 1                       ! next character in unquoted

               unquoted_str( unquoted_len: unquoted_len) = quoted_str( quoted_len: quoted_len)

!  quote followed by any other character

            else next_quote                                          ! find two adjacent quotes

!  one quote is the matching end quote

               exit scan_string                                      ! quit

            endif next_quote                                         ! find two adjacent quotes

! ----------------------------------------------------------------------

!  character is not a matching quote

         else next_char                                              ! find a quote

            unquoted_len = unquoted_len + 1                          ! next character in unquoted

            unquoted_str( unquoted_len: unquoted_len) = quoted_str( quoted_len: quoted_len)

         endif next_char                                             ! find a quote

!  check next character

         quoted_len = quoted_len + 1                                 ! next character in quoted string

      enddo scan_string                                              ! loop thru quoted string

! ----------------------------------------------------------------------

   endif got_esc                                                     ! process woth or without esc

   unquote_string = .true.

! ----------------------------------------------------------------------

return                                                               ! unquote_string()

!  unquote_string()

! **********************************************************************

end function unquote_string

! **********************************************************************

!  character_functions

! **********************************************************************

! $Id: charfunc.f90 1.8 2003/06/01 13:10:41Z Dan Release $

end module character_functions                                       ! eof


