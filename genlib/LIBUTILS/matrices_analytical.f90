module analytic_matrix

 use matrix
 use random

contains

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine invert_analitically_coef_matrix___(out,term,bottom,nL,np,matin)
 implicit none
  integer             :: i,j,k,nL,np,l1,l2
  complex(8)          :: bottom(np,nL),out(np),term(np),denom
  complex(8),optional :: matin(np,nL,3,3)

  if(present(matin))then
   do i=1,np
    do k=1,nL
     bottom(i,k) =     - abs(matin(i,k,1,2))**2   * matin(i,k,2,2) - &
                   &     abs(matin(i,k,1,3))**2   * matin(i,k,2,2) &
                   &       + matin(i,k,1,2)       * matin(i,k,2,3) *  conjg(matin(i,k,1,3))  + &
                   &         matin(i,k,1,3) * conjg(matin(i,k,1,2))*  conjg(matin(i,k,2,3))
     denom       =  matin(i,k,2,2)**2 -abs(matin(i,k,2,3))**2

     if(abs(denom)<1.d-17)  denom=1.d-17

     bottom(i,k) = bottom(i,k) / denom + matin(i,k,1,1)
    enddo
   enddo
  endif

   do i=1,np
    out(i) =  sum( CMPLX(1.d0,0.d0,kind=8)  /  ( bottom(i,:) +  term(i) ) ) / dble(nL)
   enddo

 return
 end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine invert_analitically_coef_matrix2___(out,term,bottom,a1,b1,b2,nL,np,matin)
 implicit none

  integer             :: i,j,k,nL,np,l1,l2
  complex(8)          :: bottom(np,nL),out(np),term(np),denom,a1(np,nL),b1(np,nL),b2(np,nL)
  complex(8),optional :: matin(np,nL,3,3)

  if(present(matin))then
   do i=1,np
    do k=1,nL

     bottom(i,k) =     - abs(matin(i,k,1,2))**2   * matin(i,k,2,2) - &
                   &     abs(matin(i,k,1,3))**2   * matin(i,k,2,2) &
                   &       + matin(i,k,1,2)       * matin(i,k,2,3) *  conjg(matin(i,k,1,3))  + &
                   &         matin(i,k,1,3) * conjg(matin(i,k,1,2))*  conjg(matin(i,k,2,3)) 
     denom       =  matin(i,k,2,2)**2 -abs(matin(i,k,2,3))**2

     if(abs(denom)<1.d-14)then
       write(*,*) 'ERROR invmat analatycally bad behavior'
       denom=1.d-14
     endif

     bottom(i,k) =  bottom(i,k) / denom + matin(i,k,1,1) 
     b1(i,k)     =      matin(i,k,2,2)/denom
     a1(i,k)     =      matin(i,k,1,1)
     b2(i,k)     = (abs(matin(i,k,1,3))**2)/denom
    enddo
   enddo
  endif

  do i=1,np
    out(i) =  sum( (  (term(i)+a1(i,:))*b1(i,:) - b2(i,:)    )  /  ( bottom(i,:) +  term(i) ) ) / dble(nL)
  enddo

 return
 end subroutine


!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine invert_analitically_coef_matrix2_nosum(out,term,bottom,a1,b1,b2,nL,np)
 implicit none
  integer             :: i,j,k,nL,np,l1,l2
  complex(8)          :: bottom(np,nL),out(nL,np),term(np),denom,a1(np,nL),b1(np,nL),b2(np,nL)
  do j=1,nL
  do i=1,np
    out(j,i) =   (  (term(i)+a1(i,j))*b1(i,j) - b2(i,j)    )  /  ( bottom(i,j) +  term(i) ) 
  enddo
  enddo
 return
 end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!

 subroutine invert_analitically_coef_matrix_no_sum(out,term,bottom,nL,np)
 implicit none
  integer             :: i,j,k,nL,np,l1,l2
  complex(8)          :: bottom(np,nL),out(nL,np),term(np),denom
  do j=1,nL
   do i=1,np
     out(j,i) =   1.d0 /  ( bottom(i,j) +  term(i) )  
   enddo
  enddo
 return
 end subroutine

!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!
!*****************************************************!


end module
