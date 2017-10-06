program projectorback

use strings    , only : replace_in_string,string,assignment (=)
use StringManip ,only : StrInt2,toString
use matrix,      only : write_array

implicit none

integer                :: kk_,i,j,k,l,pub_dmft_points
real(8)                :: mmu,ttemp
integer                :: num,len__,status__,paramagnetic
character(200)         :: value(2000)
character(200)         :: filename
complex(8),allocatable :: green_2(:,:),green_(:,:,:),green_temp(:,:)
type(string)           :: cc_
logical,allocatable    :: mask_proj(:)
integer                :: nproj
integer                :: k_,l_,iii,kk,ll,newatom
character(2000)        :: tt,uu
integer                :: nnn,kkk,lll
integer,allocatable    :: dimer_table(:,:)
logical                :: check,check_test
integer                :: ndimer,ii,atom_number


  write(*,*) '====================================='
  write(*,*) '       projecting back sigma         '
  write(*,*) '====================================='


  num=COMMAND_ARGUMENT_COUNT()
  if(num>size(value))then
   write(*,*) 'ERROR size value too small'
   stop
  endif
  do i=1,num
    call GET_COMMAND_ARGUMENT(i, value(i), len__, status__)
    write(*,*) 'MY ARGUMENTS : ', value(i)
  enddo
  if(num/=1)then
    write(*,*) 'ERROR combining dimer script is expecting 1 argument'
    stop
  endif

  INQUIRE(file='compound_tot_orb',EXIST=check)
  if(.not.check)then
   write(*,*) 'ERROR compounds_tot file does not exist'
   stop
  endif
  open(unit=444,file='compound_tot_orb')
  nnn=0

  call get_atom_number(atom_number,value(1))

  write(*,*) 'THIS IS ATOM : ', atom_number

  do ii=1,atom_number
   read(444,*,end=32) nnn
  enddo

  32 continue
  write(*,*) 'COMPOUNDS WITH [x] total elements : ', nnn
  close(444)

  INQUIRE(file='mask_dimer',EXIST=check)
  if(.not.check) then
   write(*,*) 'ERROR doing single site DMFT'
   stop
  endif

  i=0
  open(unit=313,file='mask_dimer')
  do
   read(313,*,end=54)
   i=i+1
  enddo
  54 continue
  allocate(dimer_table(i,2))
  dimer_table=0
  rewind(313)
  do j=1,i
   read(313,*) dimer_table(j,1),dimer_table(j,2)
   write(*,*) 'DIMER TABLE = ', dimer_table(j,1:2)
  enddo
  close(313)

  call process_inputs

  if(newatom==0)then
   ndimer=1
  else
   ndimer=2
  endif
  write(*,*) 'NDIMER = ', ndimer


  open(unit=1010,file='mask_projections')
  write(*,*) 'reading  mask'
  allocate(mask_proj(nnn*ndimer))
  do j=1,atom_number
   read(1010,*)  (mask_proj(i),i=1,nnn*ndimer)
  enddo
  write(*,*)    (mask_proj(i),i=1,nnn*ndimer)
  nproj=0
  do i=1,nnn*ndimer
    if(mask_proj(i)) nproj=nproj+1
  enddo
  nproj=nproj/ndimer
  write(*,*) 'THERE ARE [x] KEPT FUNCTIONS : ', nproj
  close(1010)


  write(*,*) 'scanning file .... : ', TRIM(ADJUSTL(value(1)))

  INQUIRE(file=TRIM(ADJUSTL(value(1))),exist=check_test)

  if(.not.check_test)then
   write(*,*) 'argument input file does not exist, file name is : ', TRIM(ADJUSTL(value(1))) 
   stop
  endif

  if(.not.allocated(green_temp)) allocate(green_temp(ndimer*nproj,ndimer*nproj))
  open(unit=10001,file=TRIM(ADJUSTL(value(1))),form='unformatted')
  j=0
  do
    read(10001,end=63) green_temp
    call write_array(green_temp,' GREEN TEMP ', short=.true.,unit=6)
    j=j+1
  enddo
  63 continue
  write(*,*) ' number of frequencies -first round : ', j
  pub_dmft_points=j
  rewind(10001)
  if(j==0)then
   write(*,*) 'ERROR FILE IS EMPTY : ', value(1)
   stop
  endif

  if(.not.allocated(green_)) allocate(green_(pub_dmft_points,ndimer*nproj,ndimer*nproj))

  green_=0.
  j=0
  do
     read(10001,end=66)  green_temp
     j=j+1
     green_(j,:,:) = green_temp
  enddo
  66 continue
  write(*,*) 'number of frequencies -second round: ', j
  if(j/=pub_dmft_points)then
   write(*,*) 'ERROR error in project back, number of frequencies do not match'
   stop
  endif
  close(10001)

  allocate(green_2(nnn,nnn))

  do i=1,ndimer**2
    write(*,*) 'output files are : ', TRIM(ADJUSTL(value(i)))
    write(*,*) 'erasing these files first...'
    call system("rm "//TRIM(ADJUSTL(value(i))))
  enddo

  write(*,*) 'now generates them again ...'
  do i=1,ndimer**2
  open(unit=2000,file=TRIM(ADJUSTL(value(i))),form='unformatted')
  write(*,*) 'WRITING FILE : ', TRIM(ADJUSTL(value(i))) 
  do j=1,pub_dmft_points
        green_2=0.
     !-------------------------------!
      select case(i) 
      case(1)
        kk=0;ll=0;kkk=0;lll=0;
      case(3)
        kk=nproj;ll=nproj;kkk=nnn;lll=nnn;
      case(2)
        kk=0;ll=nproj;kkk=0;lll=nnn;
      case(4)
        kk=nproj;ll=0;kkk=nnn;lll=0;
      end select
      k_=0
      do k=1,nnn
       if(mask_proj(k+kkk))then
        k_=k_+1
        l_=0
        do l=1,nnn
            if(mask_proj(l+lll))then
              l_=l_+1
              green_2(k,l) = green_(j,k_+kk,l_+ll)
            endif
        enddo
       endif
      enddo
     !-------------------------------!
     write(2000) green_2
  enddo
  close(2000)
  enddo

contains

!===================================================================================================!

subroutine process_inputs
implicit none

  filename=TRIM(ADJUSTL(value(1)))
  write(*,*) 'TARGET COMBINED PROJECTED GREEN FUNCTION FILE = ',TRIM(ADJUSTL(filename))
  write(*,*) 'pub_dmft_points   : ',pub_dmft_points

  value(1)=TRIM(ADJUSTL(filename))
  value(2)=TRIM(ADJUSTL(filename))//"_dimer"

  do i=LEN(filename),1,-1
   if(filename(i:i)=="t")then
     exit
   endif
  enddo
  i=i+1

  do j=i,2000
   if(filename(j:j)=="_")then
     exit
   endif
  enddo
  j=j-1

  uu=TRIM(ADJUSTL(filename(j+1:LEN(filename))))
  write(*,*) 'ATOM NUMBER : ', filename(i:j) 

  newatom=StrInt2(TRIM(ADJUSTL(filename(i:j))))

  do kk_=1,size(dimer_table,1)
   if(dimer_table(kk_,1)==newatom)then
    newatom=dimer_table(kk_,2) 
    exit
   endif
  enddo

  write(*,*) 'NEW ATOM : ', newatom
  cc_=TRIM(ADJUSTL(value(1)))
  write(*,*) 'replace .... = ',TRIM(ADJUSTL(filename(i:LEN(filename))))
  write(*,*) '...with      = ',TRIM(ADJUSTL(toString(newatom)))//TRIM(ADJUSTL(uu))
  call replace_in_string(cc_,TRIM(ADJUSTL(filename(i:LEN(filename)))),TRIM(ADJUSTL(toString(newatom)))//TRIM(ADJUSTL(uu)),'last')
  value(3)=cc_
  value(4)=TRIM(ADJUSTL(value(3)))//"_dimer"
  write(*,*) 'reference file is : ', TRIM(ADJUSTL(filename)) 

  write(*,*) 'value(1) : ', trim(adjustl(value(1)))
  write(*,*) 'value(2) : ', trim(adjustl(value(2)))
  write(*,*) 'value(3) : ', trim(adjustl(value(3)))
  write(*,*) 'value(4) : ', trim(adjustl(value(4)))


end subroutine

!===================================================================================================!

subroutine get_atom_number(newatom,filename)
implicit none
integer       :: newatom,j,i
character*(*) :: filename
  do i=LEN(filename),1,-1
   if(filename(i:i)=="t")then
     exit
   endif
  enddo
  i=i+1
  do j=i,2000
   if(filename(j:j)=="_")then
     exit
   endif
  enddo
  j=j-1
  write(*,*) 'ATOM NUMBER : ', filename(i:j)
  newatom=StrInt2(TRIM(ADJUSTL(filename(i:j))))
end subroutine
!===================================================================================================!

end program
