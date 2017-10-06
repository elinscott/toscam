module listemod
use sorting 

!********************************************************!
  type llist                     
      integer             :: index        
      complex(8)          :: comp !contenu liste
      real(8)             :: rr   !comparaison liste     
      type(llist),pointer :: aft   
      type(llist),pointer :: bef
  end type
!********************************************************!

type(llist),dimension(:),allocatable,target  :: adress
integer                                      :: ntot,nmax


contains


!=================================================!
!=================================================!
!=================================================!

 subroutine initlist(nel)
 implicit none
 integer :: nel
  if(allocated(adress))deallocate(adress)
  allocate(adress(nel))
  nmax=nel;ntot=0
 end subroutine

!=================================================!
!=================================================!
!=================================================!

 subroutine addel(rr,comp)
 implicit none
 real(8)     :: rr
 complex(8)  :: comp
 integer     :: i,j
 
 ntot=ntot+1

 if(ntot==1)then
  adress(nmax)%bef=>adress(1)
  adress(nmax)%aft=>null()
  adress(nmax)%index=nmax
  adress(nmax)%rr=1.d20
  adress(nmax)%comp=0.
  adress(1)%bef=>null()
  adress(1)%aft=>adress(nmax)
  adress(1)%index=1
  adress(1)%rr=rr
  adress(1)%comp=comp
  return
 else
  i=findn(rr) !element gauche
  adress(ntot)%aft=>adress(i)%aft
  adress(ntot)%bef=>adress(i)
  adress(ntot)%index=ntot
  adress(i)%aft=>adress(ntot)
  adress(adress(i)%aft%index)%bef=>adress(i)
  adress(ntot)%rr=rr
  adress(ntot)%comp=comp
 endif

 end subroutine

!=================================================!
!=================================================!
!=================================================!

 function findn(rr)
 implicit none
 real(8) :: rr
 integer :: findn

 !bissection sur la liste
  findn=1

 end function

!=================================================!
!=================================================!
!=================================================!

subroutine test_liste
integer    :: i,j,k,l
complex(8) :: rr

rr=(1.d0,1.d0)
call initlist(1000)
do i=1,10
 call addel(dble(-i),rr)
enddo

do i=1,10
 write(*,*) i,adress(i)%index,adress(i)%rr
enddo

end subroutine


!=================================================!
!=================================================!
!=================================================!

end module



