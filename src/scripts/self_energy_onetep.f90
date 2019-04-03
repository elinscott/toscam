program selfit
implicit none

integer                  :: i,j,k,l,kk
character*20             :: filename
real*8,allocatable       :: mat(:,:)
 
write(*,*) 'name of file:'
read(*,*) filename

open(unit=1000,file=filename)

l=0
do
 read(1000,*,end=66)
 l=l+1
enddo
66 continue

allocate(mat(l,2*l))

rewind(1000)
do i=1,l
 read(1000,*) (mat(i,2*kk-1),mat(i,2*kk),kk=1,l)
enddo

write(*,*) '===========REAL PART=========='
do i=1,l
 write(*,'(400f4.1)') (mat(i,kk),kk=1,2*l,2)
enddo
write(*,*) '===========IM PART=========='
do i=1,l
 write(*,'(400f4.1)') (mat(i,kk),kk=2,2*l,2)
enddo


close(1000)



end program
