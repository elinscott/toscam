  program example
 
  use minimization_wrapping_module
  implicit none 
  integer,parameter  :: N=200
  real(8)            :: phi(N),Orot(5,5),Otarget(5,5),rot_vec(3,3),param(N,N)
  integer            :: i,j,k,l,m,u(3)
  real(8)            :: distance,vv(3),dist(N,N,N),rot1(3,3),rot2(3,3),rot3(3,3)
  real(8)            :: ang1,ang2,ang3,diagO(5),Orot_onetep(3,3),Otilde(3,3)


   write(*,*) 'reading file Orot_onetep.txt, containing the rotation'
   write(*,*) 'matrix used in ONETEP'

   open(file='Orot_onetep.txt',unit=20)
   do i=1,3
   read(20,*) (Orot_onetep(i,j),j=1,3)
   write(*,'(5f15.5)') (Orot_onetep(i,j),j=1,3)
   enddo
   close(20)


   open(file='Otarget.txt',unit=20)
   write(*,*) 'target matrix'
   do i=1,5
   read(20,*) (Otarget(i,j),j=1,5)
   write(*,'(5f15.5)') (Otarget(i,j),j=1,5)
   enddo
   close(20)

   do i=1,N
      phi(i)=-acos(-1.d0)+2.d0*dble(i-1)/dble(N-1)*acos(-1.d0)  
   enddo

   do k=1,N
   do l=1,N
   do m=1,N
   ang1=phi(k)
   ang2=phi(l)
   ang3=phi(m)

   rot_vec=0
   rot1=0.
   rot2=0.
   rot3=0.
   rot1(3,3)=1.
   rot1(1,1)= cos(ang1)
   rot1(1,2)= sin(ang1)
   rot1(2,1)=-sin(ang1)
   rot1(2,2)= cos(ang1)

   rot2(2,2)=1.
   rot2(1,1)= cos(ang2)
   rot2(1,3)= sin(ang2)
   rot2(3,1)=-sin(ang2)
   rot2(3,3)= cos(ang2)

   rot3(1,1)=1.
   rot3(2,2)= cos(ang3)
   rot3(2,3)= sin(ang3)
   rot3(3,2)=-sin(ang3)
   rot3(3,3)= cos(ang3)

   rot_vec=MATMUL(MATMUL(rot1,rot2),rot3)

  ! Orot=dspace_rot(rot_vec)
    call rotate_our_notation(Orot,rot_vec)
    diagO=[1,2,3,4,5]
    call rearrange_columns_to_identity(5,Orot,diagO)
 

    distance=0 
    do i=1,5
    do j=1,5
      distance = distance + abs( Orot(i,j) - Otarget(i,j) ) 
    enddo
    enddo

    dist(k,l,m)=distance
    !write(102,*) k,l,distance
    enddo
    enddo 
    enddo
   
 
    write(*,*) 'minimal distance:', minval(dist)
    
    u=minloc(dist)
  
    ang1=phi(u(1))
    ang2=phi(u(2))
    ang3=phi(u(3))
    rot_vec=0
    rot1=0.
    rot2=0.
    rot3=0.
    rot1(3,3)=1.
    rot1(1,1)= cos(ang1)
    rot1(1,2)= sin(ang1)
    rot1(2,1)=-sin(ang1)
    rot1(2,2)= cos(ang1)
    rot2(2,2)=1.
    rot2(1,1)= cos(ang2)
    rot2(1,3)= sin(ang2)
    rot2(3,1)=-sin(ang2)
    rot2(3,3)= cos(ang2)
    rot3(1,1)=1.
    rot3(2,2)= cos(ang3)
    rot3(2,3)= sin(ang3)
    rot3(3,2)=-sin(ang3)
    rot3(3,3)= cos(ang3)
    rot_vec=MATMUL(MATMUL(rot1,rot2),rot3)

    write(*,*) 'final rotation matrix: '
    do i=1,3
     write(*,'(3f15.5)') (rot_vec(i,j),j=1,3)
    enddo
    write(*,*)
    write(*,*) 'final rotation angles: ', ang1,ang2,ang3
   
    call rotate_our_notation(Orot,rot_vec)
    diagO=[1,2,3,4,5]
    call rearrange_columns_to_identity(5,Orot,diagO)

    distance=0
    do i=1,5
      write(*,'(5f10.5)') (Orot(i,j),j=1,5)
    enddo
    write(*,*) 'target was:'
    do i=1,5
      write(*,'(5f10.5)') (Otarget(i,j),j=1,5)
    enddo
    write(*,*) 're-ordering is:', NINT(diagO)
    write(*,*) 'Determinant  : ',  det(rot_vec)
    write(*,*) 'O O^T        : ',  matmul(rot_vec,transpose(rot_vec))


    Otilde= matmul(Orot_onetep,rot_vec)

   write(*,*) 'NEW MATRIX FOR ONETEP:'
   do i=1,3
   write(*,'(5f15.5)') (Otilde(i,j),j=1,3)
   enddo
   close(20)



end program example 
