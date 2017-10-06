program histogram
use sorting
use StringManip

implicit none

real(8),allocatable :: mat(:)
integer,allocatable :: order(:)
integer             :: i,j,u,v,LLL,LL,col
real(8)             :: cutoff

 write(*,*) 'Please enter cutoff for keeping relevant states'
 read(*,*) cutoff
 write(*,*) 'Please enter the angular momentum (L=2 or 3,1 not implemented yet)'
 read(*,*) LLL

 select case(LLL)
 case(2)
  LL=13
  col=12
 case(3)
  LL=17
  col=16
 case default
  write(*,*) 'ERROR, L not implemented'
  stop 
 end select

 open(unit=10,file='Probability.dat')
 read(10,*)
 j=0
 do
   read(10,*,end=25)
   j=j+1
 enddo
 25 continue
 write(*,*) 'THERE ARE [x] STATES : ', j

 if(allocated(mat)) deallocate(mat,order)
 allocate(mat(j),order(j))

 rewind(10)
 read(10,*)
 do i=1,j
  read(10,*) u,v,mat(i)
 enddo

 mat=-mat
 call qsort_array(mat,order)
 mat=-mat

 write(*,*) 'Probabilities are : ', mat
 
 j=0
 do
  j=j+1
  if(mat(j)<cutoff) exit
 enddo
 write(*,*) 'We are keeping [x] states : ', j

 write(*,*) 'Probabilities : ', mat(1:j)
 write(*,*) 'states labels : ', order(1:j)

 do i=1,j
  write(*,*) 'State/Prob : ', order(i),mat(i)
  call system(" echo "//adjustl(trim(toString(mat(i))))//" `onetep.dmft.histogram.script.out "// &
                      & adjustl(trim(toString(order(i))))//"  "//adjustl(trim(toString(LL)))//"  "//adjustl(trim(toString(col)))//" |tail -1` ")
  write(*,*) 
 enddo

 write(*,*) '========== PIE CHART ============'
 do i=1,j
  call system(" echo "//adjustl(trim(toString(mat(i))))//" `onetep.dmft.histogram.script.out "// &
                      & adjustl(trim(toString(order(i))))//"  "//adjustl(trim(toString(LL)))//"  "//adjustl(trim(toString(col)))//" |tail -1` ")
  write(*,*)
 enddo

 close(10)

end program
