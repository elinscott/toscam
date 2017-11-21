module StringManip

private
public :: toString
public :: strint2
! INTERFACE putinString
!  MODULE PROCEDURE putinStringstring,putinStringreal, &
!                   & putinStringint,putinStringarray, &
!                   & putinStringreals,putinStringComp
! END INTERFACE
INTERFACE toString
 MODULE PROCEDURE StringToInt,StringToReal,StringToComp,StringToRealr,StringToLog
END INTERFACE

! INTERFACE ASSIGNMENT (=)
!  MODULE PROCEDURE putinStringreal,                   &
!                 & putinStringint,putinStringarray,   &
!                 & putinStringreals,putinStringComp
! END INTERFACE
! 
contains


!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
! 
! function string_vec(k1,k2,decb,dec)
! real(8)                   :: k1,k2,dec_
! integer                   :: dec,i,decb
! character(dec+decb+2)     :: aa1,aa2
! character(7+dec*2+decb*2) :: string_vec
!  dec_=1.d0
!  do i=1,dec
!   dec_=dec_*10.d0
!  enddo
!  aa1=TRIM(ADJUSTL(toString(dble(NINT(real(k1*dec_)))/dec_)))
!  aa2=TRIM(ADJUSTL(toString(dble(NINT(real(k2*dec_)))/dec_)))
!  string_vec="("//aa1//","//aa2//")"
! end function
! 
! 
! !***********************************************
! !***********************************************
! !***********************************************
! 
! subroutine define_safe_array(filename,filename2)
! implicit none
! character*(*)   :: filename
! character*(*)   :: filename2
! integer,save    :: ER=0
! integer         :: i
! 
! interface
!  logical function ch_is_space(ch)
!   character :: ch
!  end function
! end interface
! 
!  do i=1,len(filename)
!   if(i<=len(filename2))then 
!    if( CH_IS_SPACE(filename(i:i) ) ) then
!      filename2(i:i)=" "
!    else
!      filename2(i:i)=filename(i:i)
!    endif
!   endif
!  enddo
! 
!  call S_CONTROL_BLANK(filename2)
!  call S_BLANKS_DELETE(filename2)
!  call S_TAB_BLANK(filename2)
!  call S_B2U(filename2)
! 
!  do i=1,len(filename2)
!    if( CH_IS_SPACE(filename2(i:i) ) ) then
!      filename2(i:i)=" "
!    else
!      filename2(i:i)=filename2(i:i)
!    endif
!  enddo
! 
!  if(len(filename)<1) then
!   ER=ER+1
!   filename2=toString(ER)//'no valid name'
!  endif
! 
! return
! end subroutine
! 
! !***********************************************
! !***********************************************
! !***********************************************
! 
! ! Nom       :NUMTOASCII
! ! Fonction  :Transforme en CHARACTER
! 
!       CHARACTER(LEN=10) FUNCTION NUMTOASCII(I)
!       IMPLICIT NONE
!       CHARACTER(LEN=10)            :: TMP
!       INTEGER,INTENT(IN)            :: I
!       INTEGER                        :: J,K
!       J=I
!       K=1
!       NUMTOASCII(1:10)='          '
!       TMP=NUMTOASCII
!       DO WHILE (J.GT.0.OR.K.EQ.1)
!           TMP(K:K)=CHAR(ICHAR('0')+MOD(J,10))
!           K=K+1
!           J=J/10
!       END DO
!       J=INDEX(TMP,' ')
!       DO K=1,J-1
!           NUMTOASCII(K:K)=TMP(J-K:J-K)
!       END DO
!       END FUNCTION NUMTOASCII
! 
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   !----------------------!

character*30 function StringToLog(i)
logical :: i
!BUG September 6th, corrected
!call initString(StringToInt)
call initString(StringToLog)
if(i) then
 write(StringToLog,'(a1)') 'T'
else
 write(StringToLog,'(a1)') 'F'
endif
StringToLog=TRIM(ADJUSTL(StringToLog))
end function

character*30 function StringToInt(i)
integer :: i
call initString(StringToInt)
write(StringToInt,'(i8)') i
StringToInt=TRIM(ADJUSTL(StringToInt))
end function

character*30 function StringToReal(i)
real(8) :: i
call initString(StringToReal)
if(abs(i)>=1.d-8)then
 write(StringToReal,'(f15.8)') i
 StringToReal=TRIM(ADJUSTL(StringToReal))
else
 write(StringToReal,'(f17.15)') i
 StringToReal=TRIM(ADJUSTL(StringToReal))
endif
end function

character*30 function StringToRealr(i)
real(4) :: i
call initString(StringToRealr)
write(StringToRealr,'(f12.4)') i
StringToRealr=TRIM(ADJUSTL(StringToRealr))
end function

character*30 function StringToComp(i)
complex(8) :: i
call initString(StringToComp)
write(StringToComp,'(2f10.3)') real(i),aimag(i)
StringToComp=TRIM(ADJUSTL(StringToComp))
end function
! 
!    !----------------------!
! 
! !***********************************************
! !***********************************************
! !***********************************************
! 
!    subroutine HowManyLineString(string,k)
!     implicit none
!    character*(*) string 
!     integer :: i,j,k,n
!     k=0
!    do i=1,len(string)
!       if(string(i:i).eq.';') then
!          k=k+1
!         endif
!     enddo
!    return
!    end subroutine
! 
! !***********************************************
! !***********************************************
! !***********************************************
! 
!    subroutine NiethLineString(string,n,string2)
!     implicit none
!    character*(*) string 
!     character*(*) string2
!     integer :: i,j,k,n
!     i=0
!     j=0
!     k=0
!    do i=1,len(string)
!       if(string(i:i).ne.';') then
!          k=k+1
!          string2(k:k)=string(i:i)
!         else 
!          j=j+1
!          if(j==n)goto 102
!          k=0
!         endif
!     enddo
!     if(j/=n) k=0
!     102 continue
!     do i=k+1,len(string2)
!      string2(i:i)=' '
!     enddo
!     string2=ADJUSTL(string2)
!      return
!    end subroutine
! 
! !***********************************************
! !***********************************************
! !***********************************************
! 
!    subroutine reformString(string,arrayd)
!     implicit none
!    character*(*) string 
!    character*(*) arrayd
!    integer i,last
!    do last = len(string),1,-1
!       if (string(last:last).ne. ' ') goto 102
!     enddo
! 102   continue
!    do i = 1, min(len(arrayd),last)
!       if (string(i:i)==CHAR(0)) goto 101
!       arrayd(i:i)=string(i:i)
!     enddo
! 101   continue
! 
!    return
!    end subroutine
! 
!***********************************************
!***********************************************
!***********************************************

subroutine initString(STRING)
character*(*) :: STRING
integer :: i
 do i=1,len(STRING)
  STRING(i:i) = " "
 enddo
end subroutine
! 
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! 
! subroutine finalizeString(STRING)
! character*(*) :: STRING
! integer :: aa
! aa=LEN(STRING)
! if(LEN_TRIM(STRING)==aa)then
!  aa=aa-1
!  write(*,*) 'trouble finalize string'
! endif
! STRING(aa+1:aa+1)=CHAR(0)
! end subroutine
! 
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! 
! subroutine putinStringstring(STRING,sput)
! character*(*),intent(inout) :: STRING
! character*(*),intent(in)  :: sput
! integer :: aa
! aa=LEN(STRING)
! if(LEN_TRIM(STRING)+LEN(sput)>aa) then
! write(*,*) 'trouble to put texte in String'
! return
! endif
! write(STRING,*) STRING(1:LEN_TRIM(STRING))//sput
! end subroutine
! 
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! 
! subroutine putinStringreal(STRING,numb)
! character*(*),intent(inout) :: STRING
! real(8),intent(in) :: numb
! character*40  :: sput
! write(sput,*) numb
! write(STRING,*) STRING(1:LEN_TRIM(STRING))//' '//&
!  & sput(1:LEN_TRIM(sput))
! end subroutine
! 
! !***********************************************
! !***********************************************
! !***********************************************
! 
! subroutine putinStringreals(STRING,numb)
! character*(*),intent(inout) :: STRING
! real,intent(in) :: numb
! character*40  :: sput
! write(sput,'(f10.3)') numb
! write(STRING,*) STRING(1:LEN_TRIM(STRING))//' '//&
!  & sput(1:LEN_TRIM(sput))
! end subroutine
! 
! !***********************************************
! !***********************************************
! !***********************************************
! 
! subroutine putinStringComp(STRING,numb)
! character*(*),intent(inout) :: STRING
! complex(8),intent(in) :: numb
! character*40  :: sput
! write(sput,'(2f10.3)') numb
! write(STRING,*) STRING(1:LEN_TRIM(STRING))//' '//&
!  & sput(1:LEN_TRIM(sput))
! end subroutine
! 
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! 
! subroutine putinStringint(STRING,numb)
! character*(*),intent(inout) :: STRING
! integer,intent(in) :: numb
! character*40  :: sput
! write(sput,*) numb
! write(STRING,*) STRING(1:LEN_TRIM(STRING))//' '// & 
!  & sput(1:LEN_TRIM(sput))
! end subroutine
! 
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! 
! subroutine putinStringarray(STRING,array)
! character*(*),intent(inout) :: STRING
! character*(*), dimension(:),intent(in) :: array
! integer :: i
! 
! do i=1,size(array(:))
! if(i==1) write(STRING,*) STRING(1:LEN_TRIM(STRING))//' '//array(i)
! if(i>1 ) write(STRING,*) STRING(1:LEN_TRIM(STRING))//array(i)
! enddo
! end subroutine
! 
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! 
! subroutine PrintString(STRING)
! character*(*) :: STRING
! write(*,*) STRING(1:LEN_TRIM(STRING))
! end subroutine
! 
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! 
! subroutine LineBreakString(STRING)
! character*(*) :: STRING
! STRING=STRING(1:LEN_TRIM(STRING))//" ;  "
! end subroutine
! 
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! 
!       SUBROUTINE FORMTQ(STRING,NCHAR,NNN)
!       CHARACTER*(*) STRING
!       READ(*,200) STRING
!  200  FORMAT(A)
!       NNN=0
!       DO 10 I=1,NCHAR
!   10    IF (STRING(I:I).NE.' ') NNN=NNN+1
!       END subroutine
! 
! !***********************************************************
! !***********************************************************
! !***********************************************************
! !***********************************************************
! !***********************************************************
! 
!       SUBROUTINE LOGQYN(S,D,L)
!       LOGICAL      L,L2
!       CHARACTER*1  D,D2,A
!       CHARACTER*45 STRING
!       CHARACTER    S*(*)
!       IF (D.EQ.'Y') THEN
!         L=.TRUE.
!         D2='N'
!         L2=.FALSE.
!       ELSEIF (D.EQ.'N') THEN
!         L=.FALSE.
!         D2='Y'
!         L2=.TRUE.
!       ELSE
!         WRITE(*,*)' Default should be Y or N !'
!         RETURN
!       ENDIF
!       CALL STPACK(STRING,S,45)
!    1  WRITE(*,100) STRING,D
!  100  FORMAT(A45,' Y/N (def=',A1,')  : ')
!       CALL FORMTQ(A,1,NNN)
!       IF (NNN.EQ.0) RETURN
!       IF (A.EQ.'y' .OR. A.EQ.'T' .OR. A.EQ.'t') A='Y'
!       IF (A.EQ.'n' .OR. A.EQ.'F' .OR. A.EQ.'f') A='N'
!       IF (A.EQ.D) THEN
!         RETURN
!       ELSEIF (A.EQ.D2) THEN
!         L=L2
!         RETURN
!       ENDIF
!       GOTO 1
!       END subroutine
! 
! !***********************************************************
! !***********************************************************
! !***********************************************************
! !***********************************************************
! !***********************************************************
! 
! 
!       SUBROUTINE STPACK(STRING,S,N)
!       CHARACTER STRING*(*),S*(*)
!       DO 10 I=1,N
!         STRING(I:I)=S(I:I)
!         IF (S(I:I).EQ.'?') GOTO 20
!   10  CONTINUE
!   20  DO 30 J=I+1,N
!   30    STRING(J:J)=' '
!       END subroutine
! 
! !***********************************************************
! !***********************************************************
! !***********************************************************
! !***********************************************************
! !***********************************************************
! 
!       SUBROUTINE FINDNC(STRING,NCHARS,NC)
!       CHARACTER STRING(*)
!       integer i,NC,NCHARS
! 
!       DO 10 I=1,NCHARS
!   10    IF (STRING(I).NE.' ') GOTO 1
!    1  IMIN=I
!       DO 20 I=NCHARS,1,-1
!   20    IF (STRING(I).NE.' ') GOTO 2
!    2  IMAX=I
!       IF (IMIN.LE.IMAX) THEN
!         NC=1
!         J=IMIN
!         DO 30 I=IMIN,IMAX
!           IF (J.GT.IMAX) RETURN
!           IF (STRING(J).EQ.' ' .OR. STRING(J).EQ.',') THEN
!             NC=NC+1
!    3        IF (STRING(J+1).EQ.' ' .OR. STRING(J+1).EQ.',') THEN
!               J=J+1
!               GOTO 3
!             ENDIF
!           ENDIF
!           J=J+1
!   30    CONTINUE
!       ELSE
!         NC=0
!       ENDIF
!       END subroutine
! 
! !***********************************************************
! !***********************************************************
! !***********************************************************
! !***********************************************************
! !***********************************************************
! 
!       INTEGER FUNCTION CTOI (S, I)
!       CHARACTER*(*) S
!       INTEGER I
! 
! ! Attempt to read an integer from a character string, and return
! ! the result. No attempt is made to avoid integer overflow. A valid
! ! integer is any sequence of decimal digits.
! ! Returns:
! !  CTOI            : the value of the integer; if the first character
! !                    read is not a decimal digit, the value returned
! !                    is zero.
! ! Arguments:
! !  S      (input)  : character string to be parsed.
! !  I      (in/out) : on input, I is the index of the first character
! !                    in S to be examined; on output, either it points
! !                    to the next character after a valid integer, or
! !                    it is equal to LEN(S)+1.
!       INTEGER K
!       CHARACTER*1 DIGITS(0:9)
!       DATA  DIGITS/'0','1','2','3','4','5','6','7','8','9'/
!       CTOI = 0
!    10 IF (I.GT.LEN(S)) RETURN
!       IF (S(I:I).EQ.' ') THEN
!          I = I+1
!          GOTO 10
!       END IF
!       DO K=0,9
!           IF (S(I:I).EQ.DIGITS(K)) GOTO 30
!       enddo
!       RETURN
!    30 CTOI = CTOI*10 + K
!       I = I+1
!       GOTO 10
!       END function
! 
! !***********************************************************
! !***********************************************************
! !***********************************************************
! !***********************************************************
! !***********************************************************
! 
! !     ---------------------------------
!       subroutine GetName ( ic , fname )
! !     ---------------------------------
! !     A sample routine to construct a name for frame #ic
! !     can be used for generating movie sequences in a cycle.
! !     input  : ic - number of the frame
! !     output : fname - filename
!       implicit none
!       integer ic
!       character*15 fname
!       integer i1, i2, i3, i4
!       character*3 prefix
!       character*8 suffix
!       i1 = ic / 1000
!       i2 = (ic - (ic/1000)*1000)/100
!       i3 = (ic - (ic/100)*100) / 10
!       i4 = (ic - (ic/10)*10)
!       prefix = 'f1_'
!       suffix = '.gif/gif'
!       fname = prefix//d2ch(i1)//d2ch(i2)//d2ch(i3)//d2ch(i4)//suffix
!       return
!       end subroutine
! 
! !***********************************************************
! !***********************************************************
! !***********************************************************
! !***********************************************************
! !***********************************************************
! 
! !     --------------------------------
!       character function d2ch ( idig )
! !     --------------------------------
! !     Converts a digit 0 <= idig <= 9 into a corresponding character
!       implicit none
!       integer idig
!       if ( (idig .ge. 0).and.(idig.le.9) ) then
!         d2ch = char(48 + idig)
!       else
!         write(*,*) 'd2ch error : idig =', idig,' is not a digit'
!       endif
!       return
!       end function
! 
! !***************************************
! !***************************************
! !***************************************
! !***************************************
! 
! !convert string ch to a real(8) (positive only)
! Function StrDbl(ch)
! real(8) :: StrDbl
! character*(*) ch
! integer i, ipos
! real(8)  x
! 
! ipos=0
! do i=1,LEN_TRIM(ch)
!   if (ch(i:i).eq.'.') ipos=i
! end do
! if(ipos==0) then
! x=0.
! ipos=LEN_TRIM(ch)
! do i=1,ipos
!   x=x+(IACHAR(ch(i:i))-48)*10**(ipos-i)
! end do
! StrDbl=x
! return
! endif
! 
! x=0.d0
! do i=1, ipos-1
!   x=x+(IACHAR(ch(i:i))-48)*10**(ipos-1-i)
! end do
! do i=ipos+1, LEN_TRIM(ch)
!   x=x+REAL(IACHAR(ch(i:i))-48)/(10**(i-ipos))
! end do
! StrDbl=x
! END Function StrDbl
! 
! !***************************************
! !***************************************
! !***************************************
! !***************************************
! 
! !convert string ch to an integer*4 (positive only)
! INTEGER(4) Function StrInt(ch)
! character*(*) ch
! integer i, ilen
! integer*4 j
! j=0
! ilen=LEN_TRIM(ch)
! do i=1,ilen+1
!   if (i.eq.ilen) then
!     j=j+IACHAR(ch(i:i))-48
!   else
!     j=j+INT((IACHAR(ch(i:i))-48)*10**(ilen-i))
!   end if
! end do
! StrInt=j
! END Function
! 
! !***************************************
! !***************************************
! !***************************************
! !***************************************
! 
! !convert string ch to an integer (positive only)
INTEGER Function StrInt2(ch)
character*(*) ch
integer i,j,ilen
j=0;
ilen=LEN_TRIM(ch)
do i=1, ilen !+1 BUG corrected
  if (i.eq.ilen) then
    j=j+IACHAR(ch(i:i))-48
  else
    j=j+INT((IACHAR(ch(i:i))-48)*10**(ilen-i))
  end if
end do
StrInt2=j
END Function
! 
! !***************************************
! !***************************************
! !***************************************
! !***************************************

end module
