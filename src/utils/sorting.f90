module sorting
   use genvar, only: SP, DP
   use linalg
   use random

   implicit none
   private
   public :: fastsearchreal
   public :: group_data_rrr
   public :: qsort_adj_array
   public :: qsort_array

! INTERFACE askinset
!    MODULE PROCEDURE askinseti,askinsetr,askinsets,askinsetc,askinseta
! END INTERFACE
!
   INTERFACE qsort_adj_array
      MODULE PROCEDURE qsort_adj_array_c, qsort_adj_array_r, qsort_adj_array_c_s, qsort_adj_array_rs &
                    & , qsort_adj_array_veccs, qsort_adj_array_vecrs, qsort_adj_array_vecc, qsort_adj_array_vecr,&
                    &  qsort_adj_array_i
   END INTERFACE

   INTERFACE qsort_array
      MODULE PROCEDURE qsort_array_r, qsort_array_rs, qsort_array_i
   END INTERFACE
!
! INTERFACE group_data
!  MODULE PROCEDURE group_data_rr,group_data_r,group_data_c
! END INTERFACE
!
   INTERFACE fastsearchreal
      MODULE PROCEDURE fastsearchreal_r, fastsearchreal_d
   END INTERFACE
!
! INTERFACE searchwitherror
!  MODULE PROCEDURE searchwitherror_r,searchwitherror_f
! END INTERFACE
!
! INTERFACE maxval_set
!  MODULE PROCEDURE maxval_setr,maxval_setrr
! END INTERFACE
!
! INTERFACE minval_set
!  MODULE PROCEDURE minval_setr,minval_setrr
! END INTERFACE
!
! INTERFACE SSORT
!  MODULE PROCEDURE SSORT_,SSORT__
! END INTERFACE

contains
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
!  real(kind=DP) function maxval_setr(tab,set)
!  implicit none
!  real(kind=DP)  :: tab(:,:),r1,r2
!  integer :: set(:),j
!   r1=-1.d30
!   do j=1,size(tab(:,1))
!    if(askinset(j,set))then
!      r2=maxval(tab(j,:))
!      if(r2>r1) r1=r2
!    endif
!   enddo
!   maxval_setr=r1
!  end function
!
!  real(kind=DP) function maxval_setrr(tab,set)
!  implicit none
!  real  :: tab(:,:),r1,r2
!  integer :: set(:),j
!   r1=-1.d30
!   do j=1,size(tab(:,1))
!    if(askinset(j,set))then
!      r2=maxval(tab(j,:))
!      if(r2>r1) r1=r2
!    endif
!   enddo
!   maxval_setrr=r1
!  end function
!
!  real(kind=DP) function minval_setr(tab,set)
!  implicit none
!  real(kind=DP)  :: tab(:,:),r1,r2
!  integer :: set(:),j
!   r1=1.d30
!   do j=1,size(tab(:,1))
!    if(askinset(j,set))then
!      r2=minval(tab(j,:))
!      if(r2<r1) r1=r2
!    endif
!   enddo
!   minval_setr=r1
!  end function
!
!  real(kind=DP) function minval_setrr(tab,set)
!  implicit none
!  real  :: tab(:,:),r1,r2
!  integer :: set(:),j
!   r1=1.d30
!   do j=1,size(tab(:,1))
!    if(askinset(j,set))then
!      r2=minval(tab(j,:))
!      if(r2<r1) r1=r2
!    endif
!   enddo
!   minval_setrr=r1
!  end function
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
!       ! Nom : SORT3
!       SUBROUTINE SORT3(N,np,RA,RB,WKSP,IWKSP)
!       implicit real(kind=DP)(a-h,o-z)
!       DIMENSION RA(Np),RB(Np,Np),WKSP(Np),IWKSP(Np)
!       CALL INDEXX(N,np,RA,IWKSP)
!       DO 11 J=1,N
!       WKSP(J)=RA(J)
! 11    CONTINUE
!       DO 12 J=1,N
!       RA(J)=WKSP(IWKSP(J))
! 12    CONTINUE
!       do i=1,n
!       DO 13 J=1,N
!       WKSP(J)=RB(i,J)
! 13    CONTINUE
!       DO 14 J=1,N
!       RB(i,J)=WKSP(IWKSP(J))
! 14    CONTINUE
!       enddo
!       RETURN
!
!      contains
!
!     !----------------!
!     !----------------!
!     !----------------!
!     !----------------!
!     !----------------!
!
!       ! Nom : INDEXX
!       SUBROUTINE INDEXX(N,np,ARRIN,INDX)
!       implicit real(kind=DP)(a-h,o-z)
!       DIMENSION ARRIN(Np),INDX(Np)
!       DO 11 J=1,N
!       INDX(J)=J
! 11    CONTINUE
!       L=N/2+1
!       IR=N
! 10    CONTINUE
!       IF(L.GT.1)THEN
!       L=L-1
!       INDXT=INDX(L)
!       Q=ARRIN(INDXT)
!       ELSE
!       INDXT=INDX(IR)
!       Q=ARRIN(INDXT)
!       INDX(IR)=INDX(1)
!       IR=IR-1
!       IF(IR.EQ.1)THEN
!       INDX(1)=INDXT
!       RETURN
!       ENDIF
!       ENDIF
!       I=L
!       J=L+L
! 20    IF(J.LE.IR)THEN
!       IF(J.LT.IR)THEN
!       IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
!       ENDIF
!       IF(Q.LT.ARRIN(INDX(J)))THEN
!       INDX(I)=INDX(J)
!       I=J
!       J=J+J
!       ELSE
!       J=IR+1
!       ENDIF
!       GO TO 20
!       ENDIF
!       INDX(I)=INDXT
!       GO TO 10
!       return
!       END subroutine
!
!     !----------------!
!     !----------------!
!     !----------------!
!     !----------------!
!     !----------------!
!
!    end subroutine
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
!       ! Nom : SORT3
!       SUBROUTINE SORT3q(N,np,RA,RB,WKSP,IWKSP)
!       implicit real(16)(a-h,o-z)
!       DIMENSION RA(Np),RB(Np,Np),WKSP(Np),IWKSP(Np)
!       CALL INDEXX(N,np,RA,IWKSP)
!       DO 11 J=1,N
!       WKSP(J)=RA(J)
! 11    CONTINUE
!       DO 12 J=1,N
!       RA(J)=WKSP(IWKSP(J))
! 12    CONTINUE
!       do i=1,n
!       DO 13 J=1,N
!       WKSP(J)=RB(i,J)
! 13    CONTINUE
!       DO 14 J=1,N
!       RB(i,J)=WKSP(IWKSP(J))
! 14    CONTINUE
!       enddo
!       RETURN
!
!      contains
!
!     !----------------!
!     !----------------!
!     !----------------!
!     !----------------!
!     !----------------!
!
!       ! Nom : INDEXX
!       SUBROUTINE INDEXX(N,np,ARRIN,INDX)
!       implicit real(16)(a-h,o-z)
!       DIMENSION ARRIN(Np),INDX(Np)
!       DO 11 J=1,N
!       INDX(J)=J
! 11    CONTINUE
!       L=N/2+1
!       IR=N
! 10    CONTINUE
!       IF(L.GT.1)THEN
!       L=L-1
!       INDXT=INDX(L)
!       Q=ARRIN(INDXT)
!       ELSE
!       INDXT=INDX(IR)
!       Q=ARRIN(INDXT)
!       INDX(IR)=INDX(1)
!       IR=IR-1
!       IF(IR.EQ.1)THEN
!       INDX(1)=INDXT
!       RETURN
!       ENDIF
!       ENDIF
!       I=L
!       J=L+L
! 20    IF(J.LE.IR)THEN
!       IF(J.LT.IR)THEN
!       IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
!       ENDIF
!       IF(Q.LT.ARRIN(INDX(J)))THEN
!       INDX(I)=INDX(J)
!       I=J
!       J=J+J
!       ELSE
!       J=IR+1
!       ENDIF
!       GO TO 20
!       ENDIF
!       INDX(I)=INDXT
!       GO TO 10
!       return
!       END subroutine
!
!     !----------------!
!     !----------------!
!     !----------------!
!     !----------------!
!     !----------------!
!
!    end subroutine
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! function ind_cycle(i,siz,n,ind)
! implicit none
! integer          :: k,kk,siz,i,n,ind_cycle(siz),s1
! integer,optional :: ind(:)
!
! !ind:forbidden indices
!
! if(present(ind)) s1=size(ind)
! k=0; kk=0
!
! do
!  if(present(ind))then
!   if(askinset(i+kk,ind)) goto 38
!  endif
!  if(i+kk<=n) then
!   ind_cycle(k+1)=i+kk
!  else
!   ind_cycle(k+1)=i+kk-n
!  endif
!  k=k+1
!  if(k>siz-1)exit
!  38 continue
!  kk=kk+1
! enddo
!
! return
! end function
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! logical function askinseti(a,set)
! implicit none
! integer,intent(in)::a,set(:)
! integer::sizei,i
! sizei=size(set)
! askinseti=.false.
! do i=1,sizei
! if(a==set(i))then
!  askinseti=.true.
!  return
! endif
! enddo
! end function
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! logical function askinsets(a,set)
! implicit none
! real(kind=SP),intent(in)::a,set(:)
! integer::sizei,i
! sizei=size(set)
! askinsets=.false.
! do i=1,sizei
! if(a==set(i))then
!  askinsets=.true.
!  return
! endif
! enddo
! end function
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! logical function askinsetr(a,set)
! implicit none
! real(kind=DP),intent(in)::a,set(:)
! integer::sizei,i
! sizei=size(set)
! askinsetr=.false.
! do i=1,sizei
! if(a==set(i))then
!  askinsetr=.true.
!  return
! endif
! enddo
! end function
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! logical function askinsetc(a,set)
! implicit none
! complex(kind=DP),intent(in)::a,set(:)
! integer::sizei,i
! sizei=size(set)
! askinsetc=.false.
! do i=1,sizei
! if(a==set(i))then
!  askinsetc=.true.
!  return
! endif
! enddo
! end function
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! logical function askinseta(a,set)
! implicit none
! character*(*),intent(in)::a,set(:)
! integer::sizei,i
! sizei=size(set)
! askinseta=.false.
! do i=1,sizei
! if(a==set(i))then
!  askinseta=.true.
!  return
! endif
! enddo
! end function
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! integer function find_vec_in_array(vec,table)
! implicit none
! real(kind=DP),dimension(:,:),intent(in)::table
! real(kind=DP),dimension(:) :: vec
! integer :: size1,i,jj
! size1=size(table(:,1))
! find_vec_in_array=0
!  do jj=1,size1
!   if(norme(vec-table(jj,:))<1.d-3) then
!    find_vec_in_array=jj
!    return
!   endif
!  enddo
! return
! end function
!
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! subroutine SSORT__(n,ra)
!       implicit none
!       integer          :: n,l,ir,i,j
!       real(kind=DP) :: ra(:)
!       real(kind=DP) :: rra
!       l=n/2+1
!       ir=n
! 10      continue
!       if (l.gt.1) then
!             l=l-1
!             rra=ra(l)
!       else
!             rra=ra(ir)
!             ra(ir)=ra(1)
!             ir=ir-1
!             if (ir.eq.1) then
!                   ra(1)=rra
!                   return
!             end if
!       end if
!       i=l
!       j=l+l
! 20      if (j.le.ir) then
!             if (j.lt.ir) then
!                   if (ra(j).lt.ra(j+1)) j=j+1
!             end if
!             if (rra.lt.ra(j)) then
!                   ra(i)=ra(j)
!                   i=j
!                   j=j+j
!             else
!                   j=ir+1
!             end if
!             go to 20
!       end if
!       ra(i)=rra
!       go to 10
!
! return
! end subroutine
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! subroutine ordertab(dim,A,Idisorder,Iorder,Jdisorder, Jorder, numbcross)
! implicit none
! integer,intent(in)::  dim
! integer, intent(in)  ::  numbcross
! integer, dimension(numbcross), intent(in)  ::  Jdisorder
! integer, dimension(numbcross), intent(out)::  Jorder
! real(kind=DP), dimension(numbcross,dim), intent(in)  ::  Idisorder
! real(kind=DP), dimension(numbcross,dim), intent(out) ::  Iorder
! real(kind=DP), dimension(numbcross) ::  norm
! real(kind=DP)   ::  normeval,max
! integer  ::  temp
! integer::  i,j
! real(kind=DP), dimension(dim), intent(in)::  A
!
! do i=1,numbcross
! normeval=norme(Idisorder(i,:) - A(:))
! norm(i)=normeval
! enddo
! do j=1,numbcross
! max=0
! do i=1, numbcross
! if(norm(i)>=max) then
! max=norm(i)
! temp=i
!  endif
! enddo
! norm(temp)=0
! Iorder(numbcross-j+1,:) = Idisorder(temp,:)
! Jorder(numbcross-j+1) = Jdisorder(temp)
! enddo
! end subroutine
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
!  SUBROUTINE SORT (SIZ, ARR, numb)
!  implicit none
!   INTEGER*4   SIZ, i,j,numb(SIZ),TMP2
!   real(kind=DP)      ARR(SIZ), TMP
!   do i=1,SIZ
!    numb(i)=i
!   enddo
!   do i=1,SIZ-1
!    do j=i+1,SIZ
!     if(ARR(i) > ARR(j))then
!      TMP    = ARR(i)
!      ARR(i) = ARR(j)
!      ARR(j) = TMP
!      TMP2   = numb(i)
!      numb(i)= numb(j)
!      numb(j)= TMP2
!     endif
!    enddo
!   enddo
! end subroutine
!
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   function compare(f, g)
      implicit none
      real(kind=DP) :: f, g
      integer :: compare
      if (f < g) then
         compare = -1
      else
         compare = 1
      endif
   end function

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   function compare_r(f, g)
      real :: f, g
      integer :: compare_r
      if (f < g) then
         compare_r = -1
      else
         compare_r = 1
      endif
   end function

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!
! subroutine group_data_c(n1,array,array2,lbin)
! implicit none
! integer     :: n1
! real(kind=SP)     :: array(n1),dist,temp,good
! complex     :: array2(n1)
! integer     :: lbin,i,j,k,l,m,n,siz,count
! complex     :: mean(n1)
! integer     :: howmany(n1)
!
! siz=n1; temp=0.;count=0; good=0.0; howmany=0; mean=0.
!
! do i=1,siz
!  temp=array(i)
!  dist=abs(temp-good)
!  if(dist>3.d-3.or.count==0)then
!   count=count+1
!   good=temp
!   howmany(count)=1
!   mean(count)=array2(i)
!  else
!   howmany(count)=howmany(count)+1
!   mean(count)=mean(count)+array2(i)
!  endif
! enddo
!
! lbin=count
! do i=1,lbin
! temp=dble(howmany(i))
!  if(abs(temp)>1.d-4) then
!   array2(i)=mean(i)/temp
!  else
!   array2(i)=0.
!  endif
! enddo
!
! return
! end subroutine
!
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   subroutine group_data_rrr(n1, array, array2, array2b, lbin)
      implicit none
      integer     :: n1
      real(kind=DP)     :: array(n1), dist, temp, good
      real(kind=DP)     :: array2(n1), array2b(n1)
      integer     :: lbin, i, siz, count
      real(kind=DP)     :: mean(n1)
      integer     :: howmany(n1)

      siz = n1; temp = 0.; count = 0; good = 0.d0; howmany = 0; mean = 0.

      do i = 1, siz
         temp = array(i); dist = abs(temp - good)
         if (dist > 3.d-3 .or. count == 0) then
            count = count + 1
            howmany(count) = 1
            mean(count) = array2(i)
            array2b(count) = array(i)
            good = temp
         else
            howmany(count) = howmany(count) + 1
            mean(count) = mean(count) + array2(i)
         endif
      enddo

      lbin = count
      do i = 1, lbin
         temp = dble(howmany(i))
         if (abs(temp) > 1.d-4) then
            array2(i) = mean(i)/temp
         else
            array2(i) = 0.
         endif
      enddo

      return
   end subroutine

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!
! subroutine group_data_rr(n1,array,array2,lbin)
! implicit none
! integer     :: n1
! real(kind=DP)     :: array(n1),dist,temp,good
! real(kind=DP)     :: array2(n1)
! integer     :: lbin,i,j,k,l,m,n,siz,count
! real(kind=DP)     :: mean(n1)
! integer     :: howmany(n1)
! siz=n1; temp=0.;count=0; good=0.d0; howmany=0; mean=0.
! do i=1,siz
!  temp=array(i)
!  dist=abs(temp-good)
!  if(dist>3.d-3.or.count==0)then
!   count=count+1
!   howmany(count)=1
!   mean(count)=array2(i)
!   good=temp
!  else
!   howmany(count)=howmany(count)+1
!   mean(count)=mean(count)+array2(i)
!  endif
! enddo
! lbin=count
! do i=1,lbin
! temp=dble(howmany(i))
!  if(abs(temp)>1.d-4) then
!   array2(i)=mean(i)/temp
!  else
!   array2(i)=0.
!  endif
! enddo
! return
! end subroutine
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! subroutine group_data_r(n1,array,array2,lbin)
! implicit none
! integer     :: n1
! real(kind=SP)     :: array(n1),dist,temp,good
! real(kind=SP)     :: array2(n1)
! integer     :: lbin,i,j,k,l,m,n,siz,count
! real(kind=SP)     :: mean(n1)
! integer     :: howmany(n1)
! siz=n1; temp=0.;count=0; good=0.d0; howmany=0; mean=0.
! do i=1,siz
!  temp=array(i)
!  dist=abs(temp-good)
!  if(dist>3.d-3.or.count==0)then
!   count=count+1
!   howmany(count)=1
!   mean(count)=array2(i)
!   good=temp
!  else
!   howmany(count)=howmany(count)+1
!   mean(count)=mean(count)+array2(i)
!  endif
! enddo
!
! lbin=count
! do i=1,lbin
! temp=dble(howmany(i))
!  if(abs(temp)>1.d-4) then
!   array2(i)=mean(i)/temp
!  else
!   array2(i)=0.
!  endif
! enddo
! return
! end subroutine
!
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   !-------------!

   subroutine qsort_adj_array_i(array, order)
      implicit none
      integer, dimension(:)            :: array
      integer, dimension(size(array))  :: backup
      integer, dimension(size(array))  :: order
      integer                          :: j
      backup = array
      do j = 1, size(array)
         array(j) = backup(order(j))
      enddo
   end subroutine

   !-------------!

   subroutine qsort_adj_array_veccs(array, order)
      implicit none
      complex(kind=SP), dimension(:, :)                          :: array
      complex(kind=SP), dimension(size(array, 1), size(array, 2))  :: backup
      integer, dimension(size(array, 2))                   :: order
      integer                                             :: j
      backup = array
      do j = 1, size(array, 2)
         array(:, j) = backup(:, order(j))
      enddo
   end subroutine

   !-------------!

   subroutine qsort_adj_array_vecrs(array, order)
      implicit none
      real(kind=SP), dimension(:, :)                          :: array
      real(kind=SP), dimension(size(array, 1), size(array, 2))  :: backup
      integer, dimension(size(array, 2))                :: order
      integer                                          :: j
      backup = array
      do j = 1, size(array, 2)
         array(:, j) = backup(:, order(j))
      enddo
   end subroutine

   !-------------!

   subroutine qsort_adj_array_vecc(array, order)
      implicit none
      complex(kind=DP), dimension(:, :)                          :: array
      complex(kind=DP), dimension(size(array, 1), size(array, 2))  :: backup
      integer, dimension(size(array, 2))                :: order
      integer                                          :: j
      backup = array
      do j = 1, size(array, 2)
         array(:, j) = backup(:, order(j))
      enddo
   end subroutine

   !-------------!

   subroutine qsort_adj_array_vecr(array, order)
      implicit none
      real(kind=DP), dimension(:, :)                          :: array
      real(kind=DP), dimension(size(array, 1), size(array, 2))  :: backup
      integer, dimension(size(array, 2))                :: order
      integer                                          :: j
      backup = array
      do j = 1, size(array, 2)
         array(:, j) = backup(:, order(j))
      enddo
   end subroutine

   !-------------!

   subroutine qsort_adj_array_r(array, order)
      implicit none
      real(kind=DP), dimension(:)            :: array
      real(kind=DP), dimension(size(array))  :: backup
      integer, dimension(size(array)) :: order
      integer                         :: j
      backup = array
      do j = 1, size(array)
         array(j) = backup(order(j))
      enddo
   end subroutine

   !-------------!

   subroutine qsort_adj_array_rs(array, order)
      implicit none
      real, dimension(:)            :: array
      real, dimension(size(array))  :: backup
      integer, dimension(size(array)) :: order
      integer                         :: j
      backup = array
      do j = 1, size(array)
         array(j) = backup(order(j))
      enddo
   end subroutine

   !-------------!

   subroutine qsort_adj_array_c(array, order)
      implicit none
      complex(kind=DP), dimension(:)            :: array
      complex(kind=DP), dimension(size(array))  :: backup
      integer, dimension(size(array))     :: order
      integer                             :: j
      backup = array
      do j = 1, size(array)
         array(j) = backup(order(j))
      enddo
   end subroutine

   !-------------!

   subroutine qsort_adj_array_c_s(array, order)
      implicit none
      complex, dimension(:)            :: array
      complex, dimension(size(array))  :: backup
      integer, dimension(size(array))  :: order
      integer                          :: j
      backup = array
      do j = 1, size(array)
         array(j) = backup(order(j))
      enddo
   end subroutine

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   subroutine qsort_array_r(array, order2)
      implicit none
      real(kind=DP), dimension(:)                    :: array
      real(kind=DP), dimension(size(array))          :: backup
      integer, dimension(size(array))          :: order
      integer, dimension(size(array)), optional :: order2
      integer                                 :: i
      do i = 1, size(order)
         order(i) = i
      enddo
      call qsort_sort(array, order, 1, size(array))
      do i = 1, size(order)
         backup(i) = array(order(i))
      enddo
      array = backup
      if (present(order2)) order2 = order
   end subroutine

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   subroutine qsort_array_rs(array, order2)
      implicit none
      real, dimension(:)                       :: array
      real, dimension(size(array))             :: backup
      integer, dimension(size(array))          :: order
      integer, dimension(size(array)), optional :: order2
      integer                                 :: i
      do i = 1, size(order)
         order(i) = i
      enddo
      call qsort_sort_r(array, order, 1, size(array))
      do i = 1, size(order)
         backup(i) = array(order(i))
      enddo
      array = backup
      if (present(order2)) order2 = order
   end subroutine

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   subroutine qsort_array_i(array_, order2)
      implicit none
      integer, dimension(:)                    :: array_
      real, dimension(size(array_))            :: array
      real, dimension(size(array))             :: backup
      integer, dimension(size(array))          :: order
      integer, dimension(size(array)), optional :: order2
      integer                                 :: i
      do i = 1, size(order)
         order(i) = i
      enddo
      array = float(array_)
      call qsort_sort_r(array, order, 1, size(array))
      do i = 1, size(order)
         backup(i) = array(order(i))
      enddo
      array = backup
      if (present(order2)) order2 = order
      array_ = NINT(array)
   end subroutine

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   recursive subroutine qsort_sort(array, order, left, right)
      implicit none
      real(kind=DP), dimension(:)         :: array
      integer, dimension(:)         :: order
      integer                       :: left
      integer                       :: right
      integer                       :: i
      integer                       :: last
      if (left .ge. right) return
      call qsort_swap(order, left, qsort_rand(left, right))
      last = left
      do i = left + 1, right
         if (compare(array(order(i)), array(order(left))) .lt. 0) then
            last = last + 1
            call qsort_swap(order, last, i)
         endif
      enddo
      call qsort_swap(order, left, last)
      call qsort_sort(array, order, left, last - 1)
      call qsort_sort(array, order, last + 1, right)

   end subroutine qsort_sort

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   recursive subroutine qsort_sort_r(array, order, left, right)
      implicit none
      real, dimension(:)             :: array
      integer, dimension(:)         :: order
      integer                       :: left
      integer                       :: right
      integer                       :: i
      integer                       :: last
      if (left .ge. right) return
      call qsort_swap(order, left, qsort_rand(left, right))
      last = left
      do i = left + 1, right
         if (compare_r(array(order(i)), array(order(left))) .lt. 0) then
            last = last + 1
            call qsort_swap(order, last, i)
         endif
      enddo
      call qsort_swap(order, left, last)
      call qsort_sort_r(array, order, left, last - 1)
      call qsort_sort_r(array, order, last + 1, right)
   end subroutine

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   subroutine qsort_swap(order, first, second)
      implicit none
      integer, dimension(:)         :: order
      integer                       :: first, second
      integer                       :: tmp
      tmp = order(first)
      order(first) = order(second)
      order(second) = tmp
   end subroutine

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   integer function qsort_rand(lower, upper)
      implicit none
      integer                       :: lower, upper
      real                          :: r
      call random_number_wrapper(r)
      qsort_rand = lower + nint(r*(upper - lower))
   end function qsort_rand

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!
!     recursive function qfindmm(a, b, l, r) result(res)
!       !bissection recursive
!       implicit none
!       integer(4)                           :: res
!       integer(4), intent(in)               :: l,r
!       real(kind=DP), intent(in)                   :: b
!       real(kind=DP), dimension(:)                 :: a
!       integer(4)                           :: i,j,m
!       real(kind=DP)                               :: x
!       i=l
!       j=r
!       m=(l+r)/2
!       x=a(m)
!       if(x.eq.b) then
!        res=m
!       else if(x.gt.b) then
!        res=qfindmm(a,b,l,m-1)
!       else if(x.lt.b) then
!        res=qfindmm(a,b,m+1,r)
!       endif
!     end function
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
!   subroutine fastsearchrealwithdeg(xx,tab,jmin,jmax,err)
!   implicit none
!   integer          :: i,j,k,l,siz,jmin,jmax
!   real(kind=DP)          :: xx,tab(:),err2
!   real(kind=DP),optional :: err
!
!   err2=1.d-3
!   if(present(err)) err2=err
!
!   siz=size(tab(:))
!   i=searchwitherror(xx,tab,err2)
!
!   jmin=i; jmax=i
!
!   do
!    jmin=jmin-1
!    if(jmin<1) exit
!    if(abs(tab(jmin)-xx)>err2) exit
!   enddo
!   jmin=max(1,jmin)
!
!   do
!    jmax=jmax+1
!    if(jmax>siz) exit
!    if(abs(tab(jmax)-xx)>err2) exit
!   enddo
!   jmax=min(size(tab),jmax)
!
!   return
!   end subroutine
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
   function fastsearchreal_r(xx, tab)
      implicit none
      !bissection non recursive
      !LES BIN COMMENCE A GAUCHE : gauche_bin_i=tab(xi)
      integer :: i1, i2, is, siz
      real(kind=SP) :: xx, tab(:)
      integer :: fastsearchreal_r
      siz = size(tab)
      is = siz/2
      i1 = 1
      i2 = siz
      fastsearchreal_r = 0
      if (tab(1) > xx) then
         fastsearchreal_r = siz
         write (*, *) 'element pas dans le tableau'
         return
      endif
      if (tab(siz) <= xx) then
         fastsearchreal_r = siz
         return
      endif
      do
         !*********************************
         if (tab(is) <= xx) then
            i1 = is
            is = i1 + max(1, (i2 - i1)/2)
            goto 28
         endif
         !*********************************
         if (tab(is) > xx) then
            i2 = is
            is = i1 + (i2 - i1)/2
            goto 28
         endif
         !**********************************
28       continue
         if (is == siz .and. tab(is) <= xx) then
            fastsearchreal_r = is
            return
         endif
         if (is + 1 <= siz) then
         if (tab(is) <= xx .and. tab(is + 1) > xx) then
            fastsearchreal_r = is
            return
         endif
         endif
      enddo
   end function

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************

   function fastsearchreal_d(xx, tab)
      implicit none
      !bissection non recursive
      !LES BIN COMMENCE A GAUCHE : gauche_bin_i=tab(xi)
      integer :: i1, i2, is, siz
      real(kind=DP)  :: xx, tab(:)
      integer :: fastsearchreal_d
      siz = size(tab)
      is = siz/2
      i1 = 1
      i2 = siz
      fastsearchreal_d = 0
      if (tab(1) > xx) then
         fastsearchreal_d = siz
         write (*, *) 'element pas dans le tableau'
         return
      endif
      if (tab(siz) <= xx) then
         fastsearchreal_d = siz
         return
      endif
      do
         !*********************************
         if (tab(is) <= xx) then
            i1 = is
            is = i1 + max(1, (i2 - i1)/2)
            goto 28
         endif
         !*********************************
         if (tab(is) > xx) then
            i2 = is
            is = i1 + (i2 - i1)/2
            goto 28
         endif
         !**********************************
28       continue
         if (is == siz .and. tab(is) <= xx) then
            fastsearchreal_d = is
            return
         endif
         if (is + 1 <= siz) then
         if (tab(is) <= xx .and. tab(is + 1) > xx) then
            fastsearchreal_d = is
            return
         endif
         endif
      enddo
   end function

!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!***********************************************
!
!    function searchwitherror_r(xx,tab,error)
!    implicit none
!
!    !LES BIN COMMENCE A GAUCHE : gauche_bin_i=tab(xi)
!
!    integer                         :: i1,i2,is,siz,jj
!    real(kind=DP),intent(in),dimension(:) :: tab
!    real(kind=DP),intent(in)              :: xx,error
!    integer                         :: searchwitherror_r
!    siz=size(tab(:))
!    is=siz/2
!    i1=1
!    i2=siz
!    searchwitherror_r=0
!    if(tab(1)-error>xx)then
!     searchwitherror_r=0
!     return
!    endif
!    if(tab(siz)-error<=xx)then
!     searchwitherror_r=siz
!     return
!    endif
!      do jj=1,100000
!      !*********************************
!      if(tab(is)-error<xx) then
!       i1=is
!       is=i1+max(1,(i2-i1)/2)
!       goto 28
!      endif
!      !*********************************
!      if(tab(is)-error>xx)then
!       i2=is
!       is=i1+(i2-i1)/2
!       goto 28
!      endif
!      !**********************************
!      28 continue
!      if(is==siz.and.tab(is)-error<=xx)then
!       searchwitherror_r=is
!       return
!      endif
!      if(is+1<=siz)then
!      if(tab(is)-error<=xx.and.tab(is+1)-error>xx)then
!       searchwitherror_r=is
!       return
!      endif
!      endif
!      enddo
!      write(*,*) 'elements not found in search'
!      write(*,*) 'elements : ', xx
!      write(*,*) 'tab : ', tab(:)
!      is=1
!      return
!    end function
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
!    function searchwitherror_f(xx,tab,error,calls,func)
!    implicit none
!
!    !LES BIN COMMENCE A GAUCHE : gauche_bin_i=tab(xi)
!
!    integer                         :: i1,i2,is,siz,jj,calls
!    real(kind=DP),intent(in),dimension(:) :: tab
!    real(kind=DP),intent(in)              :: xx,error
!    integer                         :: searchwitherror_f
!
!    interface
!     real(kind=DP) function func(x,i)
!     implicit none
!     real(kind=DP)      :: x
!     integer      :: i
!     end function
!    end interface
!
!    calls=NINT(func(0.d0,-2))
!
!    siz=size(tab(:))
!    is=siz/2
!    i1=1
!    i2=siz
!    searchwitherror_f=0
!
!      if(func(tab(1),1)-error>xx)then
!        searchwitherror_f=0
!        goto 31
!      endif
!      if(func(tab(siz),1)-error<=xx)then
!        searchwitherror_f=siz
!        goto 31
!      endif
!
!      do jj=1,100000
!
!      !*********************************
!      if(func(tab(is),1)-error<xx) then
!       i1=is
!       is=i1+max(1,(i2-i1)/2)
!       goto 28
!      endif
!      !*********************************
!      if(func(tab(is),1)-error>xx)then
!       i2=is
!       is=i1+(i2-i1)/2
!       goto 28
!      endif
!      !**********************************
!      28 continue
!
!      if(is==siz.and.func(tab(is),1)-error<=xx)then
!       searchwitherror_f=is
!       goto 31
!      endif
!      if(is+1<=siz)then
!       if(func(tab(is),1)-error<=xx.and.func(tab(is+1),1)-error>xx)then
!        searchwitherror_f=is
!        goto 31
!       endif
!      endif
!
!
!      enddo
!
!
!      write(*,*) 'elements not found in search'
!      write(*,*) 'elements : ', xx
!      is=1
!
!
! 31   continue
!      calls=NINT(func(0.d0,-1))
!      return
!
!    end function
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
!     recursive subroutine qsortmm(a,order,l,r)
!       implicit none
!       integer*4,dimension(:)                 :: order
!       real(kind=DP), dimension(:)                   :: a
!       integer*4                              :: i,j,l,r,temp
!       real(kind=DP)                                 :: x, y
!       i = l
!       j = r
!       x = a((l+r)/2)
!       do while(i.le.j)
!         do while (a(i).lt.x)
!           i = i+1
!         enddo
!         do while (x.lt.a(j))
!           j = j-1
!         enddo
!         if (i.le.j) then
!           y = a(i)
!           a(i) = a(j)
!           a(j) = y
!           temp=order(i)
!           order(i)=order(j)
!           order(j)=temp
!           i = i+1
!           j = j-1
!         endif
!       enddo
!       if(l.lt.j) call qsortmm(a,order,l,j)
!       if(i.lt.r) call qsortmm(a,order,i,r)
!
!      return
!     end subroutine
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !     -----------------------------
!       SUBROUTINE indexx(n,arr,indx)
! !     -----------------------------
!       INTEGER n,indx(n),M,NSTACK
!       REAL arr(n)
!       PARAMETER (M=7,NSTACK=50)
!       INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
!       REAL a
!       do 11 j=1,n
!         indx(j)=j
! 11    continue
!       jstack=0
!       l=1
!       ir=n
! 1     if(ir-l.lt.M)then
!         do 13 j=l+1,ir
!           indxt=indx(j)
!           a=arr(indxt)
!           do 12 i=j-1,1,-1
!             if(arr(indx(i)).le.a)goto 2
!             indx(i+1)=indx(i)
! 12        continue
!           i=0
! 2         indx(i+1)=indxt
! 13      continue
!         if(jstack.eq.0)return
!         ir=istack(jstack)
!         l=istack(jstack-1)
!         jstack=jstack-2
!       else
!         k=(l+ir)/2
!         itemp=indx(k)
!         indx(k)=indx(l+1)
!         indx(l+1)=itemp
!         if(arr(indx(l+1)).gt.arr(indx(ir)))then
!           itemp=indx(l+1)
!           indx(l+1)=indx(ir)
!           indx(ir)=itemp
!         endif
!         if(arr(indx(l)).gt.arr(indx(ir)))then
!           itemp=indx(l)
!           indx(l)=indx(ir)
!           indx(ir)=itemp
!         endif
!         if(arr(indx(l+1)).gt.arr(indx(l)))then
!           itemp=indx(l+1)
!           indx(l+1)=indx(l)
!           indx(l)=itemp
!         endif
!         i=l+1
!         j=ir
!         indxt=indx(l)
!         a=arr(indxt)
! 3       continue
!           i=i+1
!         if(arr(indx(i)).lt.a)goto 3
! 4       continue
!           j=j-1
!         if(arr(indx(j)).gt.a)goto 4
!         if(j.lt.i)goto 5
!         itemp=indx(i)
!         indx(i)=indx(j)
!         indx(j)=itemp
!         goto 3
! 5       indx(l)=indx(j)
!         indx(j)=indxt
!         jstack=jstack+2
!         if(jstack.gt.NSTACK) write(*,*)  'NSTACK too small in indexx'
!         if(ir-i+1.ge.j-l)then
!           istack(jstack)=ir
!           istack(jstack-1)=i
!           ir=j-1
!         else
!           istack(jstack)=j-1
!           istack(jstack-1)=l
!           l=i
!         endif
!       endif
!       goto 1
!       END subroutine
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
!       SUBROUTINE SORTLI(N,RA1,IA1,IA2,IA3)
!       REAL RA1, RRA1
!       INTEGER L, N, IR, I, J, IRA1, IRA2, IRA3, IA1, IA2, IA3
!       DIMENSION RA1(*), IA1(*), IA2(*) , IA3(*)
!       L=N/2+1
!       IR=N
!  10   CONTINUE
!       IF(L.GT.1)THEN
!          L=L-1
!          RRA1=RA1(L)
!          IRA1=IA1(L)
!          IRA2=IA2(L)
!          IRA3=IA3(L)
!       ELSE
!          RRA1=RA1(IR)
!          IRA1=IA1(IR)
!          IRA2=IA2(IR)
!          IRA3=IA3(IR)
!          RA1(IR)=RA1(1)
!          IA1(IR)=IA1(1)
!          IA2(IR)=IA2(1)
!          IA3(IR)=IA3(1)
!          IR=IR-1
!          IF(IR.EQ.1)THEN
!             RA1(1)=RRA1
!             IA1(1)=IRA1
!             IA2(1)=IRA2
!             IA3(1)=IRA3
!             RETURN
!          ENDIF
!       ENDIF
!       I=L
!       J=L+L
!  20   IF(J.LE.IR)THEN
!          IF(J.LT.IR)THEN
!             IF(RA1(J).LT.RA1(J+1))J=J+1
!          ENDIF
!          IF(RRA1.LT.RA1(J))THEN
!             RA1(I)=RA1(J)
!             IA1(I)=IA1(J)
!             IA2(I)=IA2(J)
!             IA3(I)=IA3(J)
!             I=J
!             J=J+J
!          ELSE
!             J=IR+1
!          ENDIF
!          GO TO 20
!       ENDIF
!       RA1(I)=RRA1
!       IA1(I)=IRA1
!       IA2(I)=IRA2
!       IA3(I)=IRA3
!       GO TO 10
!
!       END subroutine
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
!       SUBROUTINE SORTPP(N,ITYPE,RA1,RA2,RA3)
!       REAL RA1, RA2, RA3, RRA1, RRA2, RRA3
!       INTEGER ITYPE(*), L, N, IR, I, J, IRRA1
!       DIMENSION RA1(*), RA2(*), RA3(*)
!       L=N/2+1
!       IR=N
!  10   CONTINUE
!       IF(L.GT.1)THEN
!          L=L-1
!          RRA1=RA1(L)
!          IRRA1=ITYPE(L)
!          RRA2=RA2(L)
!          RRA3=RA3(L)
!       ELSE
!          RRA1=RA1(IR)
!          IRRA1=ITYPE(IR)
!          RRA2=RA2(IR)
!          RRA3=RA3(IR)
!          RA1(IR)=RA1(1)
!          ITYPE(IR)=ITYPE(1)
!          RA2(IR)=RA2(1)
!          RA3(IR)=RA3(1)
!          IR=IR-1
!          IF(IR.EQ.1)THEN
!             RA1(1)=RRA1
!             ITYPE(1)=IRRA1
!             RA2(1)=RRA2
!             RA3(1)=RRA3
!             RETURN
!          ENDIF
!       ENDIF
!       I=L
!       J=L+L
!  20   IF(J.LE.IR)THEN
!          IF(J.LT.IR)THEN
!             IF(RA1(J).LT.RA1(J+1))J=J+1
!          ENDIF
!          IF(RRA1.LT.RA1(J))THEN
!             RA1(I)=RA1(J)
!             ITYPE(I)=ITYPE(J)
!             RA2(I)=RA2(J)
!             RA3(I)=RA3(J)
!             I=J
!             J=J+J
!          ELSE
!             J=IR+1
!          ENDIF
!          GO TO 20
!       ENDIF
!       RA1(I)=RRA1
!       ITYPE(I)=IRRA1
!       RA2(I)=RRA2
!       RA3(I)=RRA3
!       GO TO 10
!       END subroutine
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
!       SUBROUTINE SSORT_ (X, IY, N)
!       IMPLICIT NONE
! !
! !***PURPOSE  Sort an array and make the same interchanges in
! !            an auxiliary array.  The array is sorted in
! !            decreasing order.
! !
! !   Description of Parameters
! !      X - array of values to be sorted. Array will not be modified
! !      IY - array to be carried with X (all swaps are done on this array
! !            X(IY(J)) is sorted in increasing order.
! !      N - number of values in array X to be sorted
! !
! !     .. Scalar Arguments ..
!       INTEGER N
! !     .. Array Arguments ..
!       REAL(8) X ( * )
!       INTEGER IY ( * )
! !     .. Local Scalars ..
!       REAL TEMP
!       INTEGER I, J, JMAX, ITEMP
!       JMAX = N - 1
!       DO 200 I = 1, N - 1
!          ITEMP = 99999
!          DO 100 J = 1, JMAX
!             IF (X (IY (J) ) .LT.X (IY (J + 1) ) ) GOTO 100
!             ITEMP = IY (J)
!             IY (J) = IY (J + 1)
!             IY (J + 1) = ITEMP
!   100    END DO
!          IF (ITEMP.EQ.99999) GOTO 300
!          JMAX = JMAX - 1
!   200 END DO
!   300 RETURN
!       END SUBROUTINE
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
!
! SUBROUTINE creorder( e,isrt,flg )
! !!!----------------------------------------------------------------------!!!
! !!!  subroutine creorder( e,isrt,flg )                                   !!!
! !!!  This subroutine reorders four complex energies, e(i), in such       !!!
! !!!  a way that any identical energies are in the low index array        !!!
! !!!  slots.  There are a total of 15 different possibilities that        !!!
! !!!  can occur and to distinguish between these we set up a unique       !!!
! !!!  array of prime numbers and then use the product of the array        !!!
! !!!  elements to identify which case occurs.  Note that the way we       !!!
! !!!  check for identical elements could in principle result in more      !!!
! !!!  than 15 possibilities if say element 1 and 2 are close within       !!!
! !!!  tolerance and element 2 and 3 are close within tolerance it         !!!
! !!!  could occur that element 1 and 3 would not be counted as equal      !!!
! !!!  if they were not within the tolerance limit.  This case is          !!!
! !!!  identified by a negative value for the flag, flg and should be      !!!
! !!!  treated (as error) in the calling routine.                          !!!
! !!!                                                                      !!!
! !!!  Inputs:  e     - energy points (overwritten)                        !!!
! !!!  Outputs: isrt  - pointer array                                      !!!
! !!!           flg   - 1 if all energies different                        !!!
! !!!                   2 if two energies are identical                    !!!
! !!!                   3 if two pairs of energies are identical           !!!
! !!!                   4 if three energies are identical                  !!!
! !!!                   5 if all energies are identical                    !!!
! !!!                  -1 if no order is found                             !!!
! !!!----------------------------------------------------------------------!!!
!   !---------- Passed variables ----------
!   INTEGER    :: flg,isrt(4)
!   complex(kind=DP) :: e(4)
!   !---------- Local parameters ----------
!   real(kind=DP), parameter :: eps = 1.d-6
!   !---------- Local variables ----------
!   INTEGER    :: p,q,idx
!   real(kind=DP)    :: rez,imz
!   complex(kind=DP) :: z,w,e_new(4)
!
!   !---------- Set up prime number matrix ----------
!   idx = 1
!   DO p = 1,3
!      DO q = p+1,4
!         z = e(p)-e(q)
!         rez = dabs( dble(z) )
!         imz = dabs( aimag(z) )
!         IF( rez.LE.eps.AND.imz.LE.eps ) THEN
!            idx = idx*prime(p,q)
!         ENDIF
!      ENDDO
!   ENDDO
!   !---------- Find flg ----------
!   flg = 1
!   IF(idx.GT.1) flg = flg+1
!   IF(idx.GT.13) flg = flg+1
!   IF(idx.GT.35) flg = flg+1
!   IF(idx.GT.1001) flg = flg+1
!   !---------- Find new ordering ----------
!   IF(idx.EQ.1.OR.idx.EQ.2.OR.idx.EQ.42.OR.idx.EQ.30030) THEN
!      CALL setmap(isrt,1,2,3,4)
!   ELSE IF(idx.EQ.3) THEN
!      CALL setmap(isrt,1,3,2,4)
!   ELSE IF(idx.EQ.5) THEN
!      CALL setmap(isrt,1,4,2,3)
!   ELSE IF(idx.EQ.7.OR.idx.EQ.35) THEN
!      CALL setmap(isrt,2,3,1,4)
!   ELSE IF(idx.EQ.11.OR.idx.EQ.33) THEN
!      CALL setmap(isrt,2,4,1,3)
!   ELSE IF(idx.EQ.13.OR.idx.EQ.26) THEN
!      CALL setmap(isrt,3,4,1,2)
!   ELSE IF(idx.EQ.195) THEN
!      CALL setmap(isrt,1,3,4,2)
!   ELSE IF(idx.EQ.110) THEN
!      CALL setmap(isrt,1,2,4,3)
!   ELSE IF(idx.EQ.1001) THEN
!      CALL setmap(isrt,2,3,4,1)
!   ELSE
!      flg = -1
!      CALL setmap(isrt,1,2,3,4)
!   ENDIF
!   !---------- Reorder energies ----------
!   DO p = 1,4
!      e_new(p) = e(isrt(p))
!   ENDDO
!   DO p = 1,4
!      e(p) = e_new(p)
!   ENDDO
!   RETURN
! END SUBROUTINE
!
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
!
! FUNCTION prime(p,q)
! !!!------------------------------------------------------------------------!!!
! !!!  This function associates a prime number to pairs of small integers.   !!!
! !!!  The association is the following:                                     !!!
! !!!  (1,2) = 2, (1,3) = 3, (1,4) = 5                                       !!!
! !!!             (2,3) = 7, (2,4) = 11                                      !!!
! !!!                        (3,4) = 13                                      !!!
! !!!------------------------------------------------------------------------!!!
!   INTEGER :: prime
!   INTEGER :: p,q
!   prime = p*q
!   IF(prime.EQ.4.OR.prime.EQ.6.OR.prime.EQ.12) THEN
!      prime = prime+1
!      RETURN
!   ELSE IF(prime.EQ.8) THEN
!      prime = prime+3
!      RETURN
!   ENDIF
!   RETURN
! END FUNCTION prime
!
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
!
! SUBROUTINE setmap(v,a,b,c,d)
! !!!-----------------------------------------------------------!!!
! !!!                                                           !!!
! !!!  This subroutine takes in one vector, v, of integers and  !!!
! !!!  four numbers a,b,c,d and assigns them the elements of    !!!
! !!!  the vector v, which is of length four.                   !!!
! !!!                                                           !!!
! !!!-----------------------------------------------------------!!!
!   INTEGER v(4),a,b,c,d
!   v(1) = a
!   v(2) = b
!   v(3) = c
!   v(4) = d
!   RETURN
!
! END SUBROUTINE
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! SUBROUTINE Permute_TEV_g(ndim, zp, indx)
!   IMPLICIT NONE
!   complex(kind=DP), intent(inout) :: zp(ndim,4)
!   INTEGER, intent(in)       :: ndim
!   INTEGER, intent(in)       :: indx(ndim,4)
!   complex(kind=DP) :: eval(ndim)
!   INTEGER     :: i, p
!   DO i=1,4
!      DO p=1,ndim
!         eval(p) = zp(indx(p,i),i)
!      ENDDO
!      zp(:,i) = eval
!   ENDDO
! END SUBROUTINE
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! SUBROUTINE InversePermute_g(ndim, indx, wgh)
!   IMPLICIT NONE
!   complex(kind=DP), intent(inout) :: wgh(ndim,4)
!   INTEGER, intent(in)       :: ndim
!   INTEGER, intent(in)       :: indx(ndim,4)
!   complex(kind=DP)                :: eval(ndim)
!   INTEGER                   :: i, p
!   DO i=1,4
!      DO p=1,ndim
!         eval(indx(p,i)) = wgh(p,i)
!      ENDDO
!      wgh(:,i) = eval
!   ENDDO
! END SUBROUTINE
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! SUBROUTINE eig_order_real_part(ev, idxarr, ndim)
!   IMPLICIT NONE
! !!!-----------------------------------------------------------------!!!
! !!! This routine sorts complex eigenvalues of a matrix according to !!!
! !!! its real parts with the smallest in the first slot and reorders !!!
! !!! the matrices of left (row) and right (column) eigenvectors in a !!!
! !!! corresponding manner.                                           !!!
! !!!-----------------------------------------------------------------!!!
!   !---------- Passed variables ----------
!   complex(kind=DP), intent(in) :: ev(ndim)         ! Array of eigenvalues
!   INTEGER, intent(out)   :: idxarr(ndim)     ! Index array which gives proper order
!   INTEGER :: ndim                            ! Dimension of matrices
!   !---------- Parameters ----------
!   real(kind=DP), PARAMETER :: maxval = 1000.d0
!   !---------- Local variables ----------
!   LOGICAL, ALLOCATABLE :: sorted(:)
!   real(kind=DP),  ALLOCATABLE :: sortonr(:)
!   INTEGER :: p
!   INTEGER :: q
!   INTEGER :: idx
!   real(kind=DP)  :: min
!   !---------- Allocate dynamic memory storage ----------
!   ALLOCATE(sortonr(ndim), sorted(ndim))
!   !---------- Initialize arrays ----------
!   idxarr = 0
!   sortonr = DBLE(ev)
!   sorted = .FALSE.
!   !---------- Create index array for real value ----------
!   sorted = .FALSE.
!   DO p = 1,ndim
!      min = maxval
!      DO q = 1,ndim
!         IF(.NOT.sorted(q).AND.min.GT.sortonr(q)) THEN
!            min = sortonr(q)
!            idx = q
!         ENDIF
!      ENDDO
!      idxarr(p) = idx
!      sorted(idx) = .TRUE.
!   ENDDO
!   DEALLOCATE(sortonr, sorted)
!   RETURN
! END SUBROUTINE eig_order_real_part
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! SUBROUTINE permute_eigensystem(idxarr, ev, evl, evr, ndim)
!   IMPLICIT NONE
!   !---------- Passed variables ----------
!   INTEGER, intent(in)       :: idxarr(ndim)     ! Index array which gives proper order
!   complex(kind=DP), intent(inout) :: ev(ndim)         ! Array of eigenvalues
!   complex(kind=DP), intent(inout) :: evl(ndim,ndim)   ! Matrix of left eigenvectors  (row)
!   complex(kind=DP), intent(inout) :: evr(ndim,ndim)   ! Matrix of right eigenvectors (column)
!   INTEGER :: ndim                               ! Dimension of matrices
!   !---------- Local variables ------------------
!   INTEGER :: p
!   complex(kind=DP), ALLOCATABLE :: eval(:)
!   complex(kind=DP), ALLOCATABLE :: evec(:,:)
!   ALLOCATE(eval(ndim), evec(ndim,ndim))
!   !---------- Permute the eigenvalues ----------
!   DO p = 1,ndim
!      eval(p) = ev(idxarr(p))
!   ENDDO
!   ev = eval
!   !---------- Permute the right eigenvectors ----------
!   DO p = 1,ndim
!      evec(:,p) = evr(:,idxarr(p))
!   ENDDO
!   evr = evec
!   !---------- Permute the left eigenvectors ----------
!   DO p = 1,ndim
!      evec(p,:) = evl(idxarr(p),:)
!   ENDDO
!   evl = evec
!   !---------- Deallocate dynamic memory storage ----------
!   DEALLOCATE(eval, evec)
!   RETURN
! END SUBROUTINE permute_eigensystem
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! SUBROUTINE permute_eigenvals(idxarr, ev, ndim)
!   IMPLICIT NONE
!   !---------- Passed variables ----------
!   INTEGER, intent(in)       :: idxarr(ndim)     ! Index array which gives proper order
!   complex(kind=DP), intent(inout) :: ev(ndim)         ! Array of eigenvalues
!   INTEGER :: ndim                               ! Dimension of matrices
!   !---------- Local variables ------------------
!   INTEGER :: p
!   complex(kind=DP), ALLOCATABLE :: eval(:)
!   ALLOCATE(eval(ndim))
!   !---------- Permute the eigenvalues ----------
!   DO p = 1,ndim
!      eval(p) = ev(idxarr(p))
!   ENDDO
!   ev = eval
!   !---------- Deallocate dynamic memory storage ----------
!   DEALLOCATE(eval)
!   RETURN
! END SUBROUTINE permute_eigenvals
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
!
! ! sorting subroutine
!
!       SUBROUTINE EIGSRT (D, V, N, NP)
!       IMPLICIT none
!       INTEGER N, NP, i, k, j
!       REAL(8) D (NP), V (NP, NP)
!       REAL(8) P
!
!       DO 13 I = 1, N - 1
!          K = I
!          P = D (I)
!          DO 11 J = I + 1, N
!             IF (D (J) .GE.P) THEN
!                K = J
!                P = D (J)
!             ENDIF
!    11    END DO
!          IF (K.NE.I) THEN
!             D (K) = D (I)
!             D (I) = P
!             DO 12 J = 1, N
!                P = V (J, I)
!                V (J, I) = V (J, K)
!                V (J, K) = P
!    12       END DO
!          ENDIF
!    13 END DO
!       RETURN
!       END SUBROUTINE
!
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************
! !***********************************************

end module

