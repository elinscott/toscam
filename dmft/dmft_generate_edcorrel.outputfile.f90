program generateedcorrel
use StringManip ,only : StrInt2,toString
implicit none
integer                :: ii,i,j,k,l,ien,spins
integer                :: num,len__,status__,paramagnetic
character(200)         :: value(2000),filename
integer,allocatable    :: dd(:,:)
logical                :: check_mask

 num=COMMAND_ARGUMENT_COUNT()
 if(num>size(value))then
  write(*,*) 'ERROR size value too small'
  stop
 endif
 do i=1,num
   call GET_COMMAND_ARGUMENT(i, value(i), len__, status__)
   write(*,*) 'MY ARGUMENTS : ', value(i)
 enddo

 if(num/=3)then
  write(*,*) 'ERROR generating ed.correl is expecting 3 argument'
  stop
 endif

 k     = StrInt2(TRIM(ADJUSTL(value(1))))
 spins = StrInt2(TRIM(ADJUSTL(value(2))))
 filename=trim(adjustl(value(3)))

 write(*,*) '================ SIZE OF CORREL : ',k

open(unit=140,file=trim(adjustl(filename)))
allocate(dd(2*k,2*k))
dd=0.0
do i=1,2*k
 dd(i,i)=i
enddo
inquire(file='mask_sym_green_ed',exist=check_mask)
if(check_mask)then
open(unit=2222,file='mask_sym_green_ed')
 do i=1,2*k
  read(2222,*) (dd(i,j),j=1,2*k)
 enddo
close(2222)
endif

write(140,*) '###############################################'
write(140,*) '### ELECTRONIC CORRELATIONS ON THE IMPURITY ###'
write(140,*) '###############################################'
write(140,*) '########################################'
write(140,*) '###     MASK FOR GREENS FUNCTION    ###'
write(140,*) '### [ <c[up]*C[up]>  <c[up]*c[do]> ] ###'   
write(140,*) '### [ <C[do]*C[up]>  <c[do]*C[do]> ] ###'
write(140,*) '########################################'
do i=1,2*k
 write(140,'(200i4)') (dd(i,j),j=1,2*k)
enddo
write(140,*) '###########################################'
write(140,*) '### MASKS FOR SPIN/DENSITY CORRELATIONS ###'
write(140,*) '###########################################'
write(140,*) 'F' 
write(140,*) 'F' 
write(140,*) 'F' 

dd=0
if(k==1.or.spins==0)then
 do i=1,k
  dd(i,i)=i
 enddo
else 
 ii=0
 do i=1,k
  ii=ii+1
  dd(ii,ii)=ii
 enddo
 do i=1,k-1
  do j=i+1,k
    ii=ii+1
    dd(i,j)=ii
    dd(j,i)=ii
  enddo
 enddo
endif

do i=1,k
 write(140,'(200i4)') (dd(i,j),j=1,k) 
enddo

write(140,*) 2.  ,'# wmax = real freq. range [-wmax,wmax] for Sz'
write(140,*) 2.  ,'# wmax = real freq. range [-wmax,wmax] for S+-'
write(140,*) 30. ,'# wmax = real freq. range [-wmax,wmax] for N'
write(140,*) '##################################'
write(140,*) '### MASK FOR Pijk CORRELATIONS ###'
write(140,*) '##################################'
write(140,*) 'F' 
write(140,*) 0 
write(140,*) '# LIST OF TRIPLETS'
write(140,*) '# MASK OF CORRELATIONS'
write(140,*) 4.
write(140,*) '###################################'
write(140,*) '### MASK FOR Pijkl CORRELATIONS ###'
write(140,*) '###################################'
write(140,*) 'F' 
write(140,*) 0
write(140,*) 
write(140,*) 
write(140,*) 4.
write(140,*) '##################################################'
write(140,*) '### LIST OF PROJECTION VECTORS IN BASIS |.ud2> ###'
write(140,*) '##################################################'
write(140,*) 0 

close(140)
end program
