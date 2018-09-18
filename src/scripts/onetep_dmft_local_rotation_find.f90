program band
!---------------!
  !use dmftmod
  use genvar
  use init_and_close_my_sim
  ! use fortran_cuda
  use linalg
!  use DLPLOT
!  use plotlib
!  use plot2d
  USE minimization_wrapping
!---------------!
implicit none

real(8)       :: rot(3,3),Rd(5,5),search_step,dist_min,dist_max,mat(5,5),tmp(10),test(9)
integer       :: ii,i,j,k,bb(5)
character(13) :: FIT_METH
logical       :: flip

 !###############################
  call initialize_my_simulation
  testing=.false.
  fast_invmat=.true.
  use_cuda_routines=.false.
  use_cula_routines=.false.
  force_invmat_single_prec=.false.
  use_openmp_invmat=.false.
  diag_use_LU_instead_of_pivot=.false.
  flag_use_invmat_jordan=.true.
  flag_use_invmat_jordan_real=.true.
  enable_mpi_dot=.false.
 !###############################

  flip=.true.

  open(unit=30,file='matrix.txt')
  do i=1,5
   read(30,*) (tmp(j),j=1,10) 
   do j=1,5
    mat(i,j)=tmp(2*j-1)
   enddo
   write(*,'(a,100f10.3)') 'matrix',(mat(i,j),j=1,5)
  enddo
  close(30)

  open(unit=19812,file='mask_ngwfs')
  read(19812,*) (bb(ii),ii=1,5)
  write(*,*) 'NGWFs order : ', bb
  close(19812)


  search_step=0.001
  dist_min=0.000000001
  dist_max=0.000000001
  FIT_METH='MINIMIZE'

  test=0.d0 
  k=0
  do i=1,3
   do j=1,3
    k=k+1
    if(i==j) test(k)=1.
   enddo
  enddo

 call minimize_func_wrapper(distance_func,test,9,FIT_METH,1000000000,dist_min,dist_max,search_step,.false.)
 
 call get_rd
 do i=1,5
   write(*,'(a,100f10.3)') 'Rd',(Rd(i,j),j=1,5)
 enddo

 write(*,*) 'new vectors'
 k=0
 do i=1,3
  do j=1,3
   k=k+1
   rot(i,j)=test(k)
  enddo
 enddo
 do i=1,3
  rot(:,i)=rot(:,i)/norme(rot(:,i))
  write(*,'(i3,200f10.3)') i,(rot(j,i),j=1,3)
 enddo

 do i=1,3
 do j=1,3
  write(*,'(2i3,2f10.4)') i,j,scalprod(rot(:,i),rot(:,j))/norme(rot(:,i))/norme(rot(:,j))
 enddo
 enddo

 call finalize_my_simulation

contains

 SUBROUTINE rotate_our_notation(Rd,R)
 IMPLICIT NONE
 real(8), intent(in)  :: R(3,3)
 real(8), intent(out) :: Rd(5,5)
 real(8)              :: s3,r1,r2,r3,r4,r5,r6,r7,r8,r9

   s3=sqrt(3.d0)

   r1=R(1,1)
   r2=R(1,2)
   r3=R(1,3)
   r4=R(2,1)
   r5=R(2,2)
   r6=R(2,3)
   r7=R(3,1)
   r8=R(3,2)
   r9=R(3,3)

   Rd(1,1)=r2*r4+r1*r5
   Rd(1,2)=r3*r5+r2*r6
   Rd(1,3)=s3*r3*r6
   Rd(1,4)=r3*r4+r1*r6
   Rd(1,5)=2.d0*r1*r4+r3*r6

   Rd(2,1)=r5*r7+r4*r8
   Rd(2,2)=r6*r8+r5*r9
   Rd(2,3)=s3*r6*r9
   Rd(2,4)=r6*r7+r4*r9
   Rd(2,5)= 2.*r4*r7+r6*r9

   Rd(3,1)= s3*r7*r8
   Rd(3,2)= s3*r8*r9
   Rd(3,3)= (3.* r9**2 -1.d0)/2.d0
   Rd(3,4)= s3*r7*r9
   Rd(3,5)= s3*(2.*r7**2 + r9**2 - 1.d0)/2.d0

   Rd(4,1)= r2*r7+r1*r8
   Rd(4,2)=  r3*r8+r2*r9
   Rd(4,3)=  s3*r3*r9
   Rd(4,4)= r3*r7+r1*r9
   Rd(4,5)= 2.*r1*r7 +r3*r9

   Rd(5,1)=r1*r2-r4*r5
   Rd(5,2)=r2*r3-r5*r6
   Rd(5,3)=s3*(r3**2-r6**2)/2.d0
   Rd(5,4)=r1*r3-r4*r6
   Rd(5,5)=(2.*r1**2 + r3**2 - 2.*r4**2 - r6**2)/2.d0

 end subroutine

 SUBROUTINE distance_func(dist,n,vec)
    REAL(8), INTENT(OUT)   :: dist
    INTEGER, INTENT(IN)    :: n
    REAL(8), INTENT(IN)    :: vec(n)
    REAL(8)                :: vec_(3,3),Rd(5,5),tmp(5,5)
    integer                :: i,k,j
      k=0
      do i=1,3
       do j=1,3
        k=k+1
        vec_(i,j)=vec(k)
       enddo
      enddo
      call rotate_our_notation(Rd,vec_)
      tmp=Rd
      if(flip)then
      do i=1,5
       do j=1,5
        Rd(i,j)=tmp(bb(i),bb(j))
       enddo
      enddo
      endif
      dist=sum(abs(Rd-mat)**2)/25.d0
 END SUBROUTINE

 subroutine get_rd
 implicit none
 integer :: i,j,k
 real(8) :: vec_(3,3),tmp(5,5)

      k=0
      do i=1,3
       do j=1,3
        k=k+1
        vec_(i,j)=test(k)
       enddo
      enddo
      call rotate_our_notation(Rd,vec_)
      tmp=Rd
      if(flip)then
      do i=1,5
       do j=1,5
        Rd(i,j)=tmp(bb(i),bb(j))
       enddo
      enddo
      endif

 end subroutine

end program







