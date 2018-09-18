program readoverlap
use matrix
use linalg
implicit none
real(8),allocatable :: overlap(:,:),overlap_inv(:,:)
integer             :: i1,i2,i,j,k,l,m,n,n1,n2,block,nblock

 write(*,*) 'ENTER SIZE OF OVERLAP MATRIX (n1,n2) :'
 read(*,*) n1,n2
 write(*,*) 'ENTER BLOCK SIZES'
 read(*,*) block
 if(mod(n2,block)/=0)then
  write(*,*) 'ERROR n2 should be a multiple of block'
  write(*,*) 'N2    : ', n2
  write(*,*) 'Block : ', block
  write(*,*) 'MOD   : ', mod(n2,block)
  stop
 endif

 allocate(overlap(n1,n2),overlap_inv(n2,n1))
 open(unit=30,file='store_hub_overlap_matrix',form='unformatted')
 read(30) overlap
 close(30)

 open(unit=30,file='store_hub_overlap_tmatrix',form='unformatted')
 read(30) overlap_inv
 close(30)
 
 nblock=n2/block

 write(*,*) 'SCANNING OVERLAP MATRIX, LINES'
 do i=1,n1
  write(*,*) ' line [x] : ', i, sum(abs( overlap(i,:))**2)
 enddo
 write(*,*) 'SCANNING OVERLAP MATRIX, COLUMNS'
 do i=1,n2
  write(*,*) ' columns [x] : ', i, sum(abs( overlap(:,i))**2)
 enddo


 do i=1,nblock
 i1=(i-1)*block+1
 i2= i*block
 write(*,*) 'BLOCK : ', i
 call write_array( matmul(overlap_inv(i1:i2,:),overlap(:,i1:i2)),' O^T * O   ',short=.true.,unit=6)
 enddo

 deallocate(overlap,overlap_inv)

end program
