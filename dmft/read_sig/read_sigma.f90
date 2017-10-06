program readsigma
integer :: chan,i,j
complex(8),allocatable :: mat1(:,:),mat2(:,:)
write(*,*) 'ENTER CHANNELS'
read(*,*) chan
allocate(mat1(chan,chan),mat2(chan,chan))
open(unit=10,file='v1',form='unformatted')
i=0
do 
i=i+1
read(10,end=60) mat1
write(11,*) (real(mat1(j,j)),j=1,5)
enddo
60 continue
close(10);
write(*,*) 'total frequ : ',i-1

open(unit=10,file='v2',form='unformatted')
i=0
do
i=i+1
read(10,end=70) mat1
write(12,*) (real(mat1(j,j)),j=1,5)
enddo
70 continue
close(10);
write(*,*) 'total frequ : ',i-1

open(unit=10,file='v12',form='unformatted')
i=0
do
i=i+1
read(10,end=80) mat1
write(13,*) abs(mat1(1,4)), abs(mat1(4,1))
enddo
80 continue
close(10);
write(*,*) 'total frequ : ',i-1



end program
