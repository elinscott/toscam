module combinatory

real(8) :: LeviCivita(3,3,3)

contains

!***************************************
!***************************************
!***************************************
!***************************************

subroutine defineLevi
implicit none
integer :: i,j,k

LeviCivita=0.d0
do i=1,3
 do j=1,3
  do k=1,3
   if(i==j.or.i==k.or.j==k)   LeviCivita(i,j,k)= 0.d0
   if(i==1.and.j==2.and.k==3) LeviCivita(i,j,k)= 1.d0
   if(i==1.and.j==3.and.k==2) LeviCivita(i,j,k)=-1.d0
   if(i==2.and.j==3.and.k==1) LeviCivita(i,j,k)= 1.d0
   if(i==2.and.j==1.and.k==3) LeviCivita(i,j,k)=-1.d0
   if(i==3.and.j==1.and.k==2) LeviCivita(i,j,k)= 1.d0
   if(i==3.and.j==2.and.k==1) LeviCivita(i,j,k)=-1.d0
  enddo
 enddo
enddo

end subroutine

!***************************************
!***************************************
!***************************************
!***************************************

 ! return greater common divisor of two long integers
   integer Function GCD(a, b)
    integer a, b, ta, tb
    integer r
    ta=a; tb=b
    a=abs(a)
    b=abs(b)
    if (a.eq.0.or.b.eq.0) then
    GCD=1
    return
    end if
    if (a<b) then
      r=a; a=b; b=r
    end if
    r=1
    do while (r>0)
      r=MOD(a, b)
      a=b; b=r
    end do
  GCD=a; a=ta; b=tb
  end Function GCD

!***************************************
!***************************************
!***************************************
!***************************************



end module
