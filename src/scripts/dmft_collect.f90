program collect
use strings    , only : replace_in_string,string,assignment (=)
implicit none
integer                :: i,j,ien,pub_dmft_points,channels,atom,channelsb,atomb
real(8)                :: mmu,ttemp
integer                :: num,len__,status__,LL,paramagnetic
character(200)         :: value(4000)
character(200)         :: filename
complex(8)             :: frequ
complex(8),allocatable :: green_(:,:,:),green_temp(:,:),frequ_(:)
real(8),allocatable    :: occupation_matrix(:,:,:)
type(string)           :: cc_

 num=COMMAND_ARGUMENT_COUNT()
 if(num>size(value))then
  write(*,*) 'ERROR size value too small'
  stop
 endif
 do i=1,num
   call GET_COMMAND_ARGUMENT(i, value(i), len__, status__)
   write(*,*) 'MY ARGUMENTS',value(i)
 enddo

 do i=1,num
  open(unit=10001,file=TRIM(ADJUSTL(value(i))),form='unformatted')
  read(10001,end=66) ien,pub_dmft_points,ttemp,channels,channelsb,atom,atomb
  rewind(10001)
  write(*,*) 'ien,npoints,ttemp,channels,atom : ', ien,pub_dmft_points,ttemp,channels,channelsb,atom,atomb
  if(.not.allocated(green_)) then
       allocate(frequ_(pub_dmft_points),green_(pub_dmft_points,channels,channelsb), &
                                    & green_temp(channels,channelsb),occupation_matrix(channels,channels,2))
       green_=0.d0
  endif
  write(*,*) 'start loop'
  j=0
  do
   read(10001,end=66) ien,pub_dmft_points,ttemp,channels,channelsb,atom,atomb
   read(10001,end=66) occupation_matrix(:,:,1),occupation_matrix(:,:,2)
   read(10001,end=66) mmu,frequ,green_temp
     j=j+1
     green_(ien,:,:)   = green_(ien,:,:)+green_temp
     if(abs(frequ_(ien))>1.d-12)then
        if(abs(frequ_(ien)-frequ)>1.d-12)then
         write(*,*) 'ERROR IN COLLECT GREEN , FREQUENCIES DO NOT MATCH'
         write(*,*) 'IN ARRAY WE HAD  : ', frequ_(ien)
         write(*,*) 'NEW FREQUENCY IS : ', frequ
         write(*,*) 'THEY SHOULD BE THE SAME - might be due to having several k points in the same green_output files' 
         stop
        endif
     endif
     frequ_(ien)       = frequ
  enddo
  66 continue
  write(*,*) 'number of frequencies : ', j
  write(*,*) 'file/tot : ', i,num
  close(10001)
 enddo

 do i=1,num
  call system( " rm "//TRIM(ADJUSTL(value(i))) )
 enddo

  write(*,*) 'TEMP = ', ttemp
  write(*,*) 'TARGET FILE'
  filename=TRIM(ADJUSTL(value(1)))
  do i=LEN(filename),1,-1
   if(filename(i:i)=="k")then
     exit
   else
     filename(i:i)=" "
   endif
  enddo
  cc_=TRIM(ADJUSTL(filename))
  call replace_in_string(cc_,'_rank','','last')
  filename=cc_

  write(*,*) 'TARGET GREEN FUNCTION FILE = ',TRIM(ADJUSTL(filename))
  write(*,*) 'writing file      : ',mmu
  write(*,*) 'pub_dmft_points   : ',pub_dmft_points
  write(*,*) 'ttemp             : ',ttemp
  write(*,*) 'channels          : ',channels
  write(*,*) 'channelsb         : ',channelsb
  write(*,*) 'atom              : ',atom
  write(*,*) 'atomb             : ',atomb
  write(*,*) 'occupation matrix : ',occupation_matrix(:,:,1)
  !value(1) was already removed above at line 51
  call system("rm "//TRIM(ADJUSTL(filename))//" 2> /dev/null  ") 
  write(*,*) 'opening file'
  open(unit=2000,file=TRIM(ADJUSTL(filename)),form='unformatted')
  write(*,*) 'start loop'
  do j=1,pub_dmft_points
   write(2000) j,pub_dmft_points,ttemp,channels,channelsb,atom,atomb
   write(2000) occupation_matrix(:,:,1),occupation_matrix(:,:,2)
   write(2000) mmu,frequ_(j),green_(j,:,:)
  enddo
  close(2000)

stop
end program
