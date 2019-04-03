!**************************************!
!**************************************!
!**************************************!

module variables

use namelistmod
use StringManip, only : tostring
use linalg, only : norm
character(20),allocatable :: leftsymbol(:),rightsymbol(:),devicesymbol(:)
real(8),allocatable       ::  leftcoord(:,:),rightcoord(:,:),devicecoord(:,:)
real(8)                   :: leftcarth(3,3),rightcarth(3,3),devicecarth(3,3)
character(20),allocatable :: labels(:)
integer,allocatable       :: orbitals(:)
integer                   :: nlab
character(20),allocatable :: leftshiftedsymbol(:),rightshiftedsymbol(:)
real(8),allocatable       ::  leftshiftedcoord(:,:),rightshiftedcoord(:,:)

integer                   :: d3,t1,t2,s,t
real(8)                   :: DDD(3),cutoff_same_site
integer,allocatable       :: device_in_LR_shifted(:),LR_in_LR_shifted(:)
real(8)                   :: v1(3),v2(3)
integer,allocatable       :: sites_in_device(:),sites_in_lr(:)


contains

 subroutine set_var
  type(namelist_set) :: nm
   call namelist_init(nm,400,name_of_namelist='stitch_oxides')
   call putel_in_namelist(nm,d3,'d3',1,'Definition:heterostructure is stacked along d3')
   call look_for_namelist_in_file(nm,'var.txt')
   call look_for_command_line_argument(nm)
 end subroutine

end module

program stitch
use variables
implicit none
integer :: i,j,k,l,r,d

 
 call set_var

 DDD=0.0
 DDD(3)=0.
 cutoff_same_site=0.3

 l=0
 open(unit=100,file='stitch_left')
 read(100,*) nlab
 allocate(labels(nlab),orbitals(nlab))
 read(100,*) (labels(i),i=1,nlab)
 read(100,*) (orbitals(i),i=1,nlab)
 do i=1,3
   read(100,*) (leftcarth(i,j),j=1,3)
 enddo
 do 
  read(100,*,end=50) 
  l=l+1
 enddo
50 continue
 rewind(100)
  read(100,*);read(100,*);read(100,*)
  do i=1,3; read(100,*); enddo
 allocate(leftcoord(l,3),leftsymbol(l))
 allocate(leftshiftedcoord(l,3),leftshiftedsymbol(l))
 do i=1,l
   read(100,*) leftsymbol(i),(leftcoord(i,j),j=1,3)
 enddo
 close(100)
 open(unit=100,file='stitch_left_shifted')
  read(100,*);read(100,*);read(100,*)
  do i=1,3; read(100,*); enddo
  do i=1,l
   read(100,*) leftshiftedsymbol(i),(leftshiftedcoord(i,j),j=1,3)
  enddo
 close(100)

 do i=1,nlab
  write(*,*) 'ORBITAL WITH LABEL : ', trim(adjustl(labels(i)))
  write(*,*) 'HAS [x] ORBITALS   : ', orbitals(i) 
 enddo

 r=0
 open(unit=100,file='stitch_right')
 read(100,*);read(100,*);read(100,*)
 do i=1,3
   read(100,*) (rightcarth(i,j),j=1,3)
 enddo
 do
  read(100,*,end=51)
  r=r+1
 enddo
51 continue
 rewind(100)
  read(100,*);read(100,*);read(100,*)
  do i=1,3; read(100,*); enddo
 allocate(rightcoord(r,3),rightsymbol(r))
 allocate(rightshiftedcoord(r,3),rightshiftedsymbol(r))
 do i=1,r
   read(100,*) rightsymbol(i),(rightcoord(i,j),j=1,3)
 enddo
 close(100)
 open(unit=100,file='stitch_right_shifted')
  read(100,*);read(100,*);read(100,*)
  do i=1,3; read(100,*); enddo
  do i=1,r
   read(100,*) rightshiftedsymbol(i),(rightshiftedcoord(i,j),j=1,3)
  enddo
 close(100)

 d=0
 open(unit=100,file='stitch_device')
 read(100,*);read(100,*);read(100,*)
 do i=1,3
   read(100,*) (devicecarth(i,j),j=1,3)
 enddo
 do
  read(100,*,end=52)
  d=d+1
 enddo
52 continue
 rewind(100)
  read(100,*);read(100,*);read(100,*)
  do i=1,3; read(100,*); enddo
 allocate(devicecoord(d,3),devicesymbol(d))
 do i=1,d
   read(100,*) devicesymbol(i),(devicecoord(i,j),j=1,3)
 enddo
 close(100)

 do i=1,l
  leftcoord(i,:)=leftcoord(i,:)-leftcarth(d3,:)
  leftshiftedcoord(i,:)=leftshiftedcoord(i,:)-leftcarth(d3,:) - DDD
 enddo
 do i=1,r
  rightcoord(i,:)=rightcoord(i,:)+devicecarth(d3,:)
  rightshiftedcoord(i,:)=rightshiftedcoord(i,:)+devicecarth(d3,:) + DDD
 enddo

 open(unit=100,file='stitch_full.xyz')
 write(100,*) l+d+r
 write(100,*) (trim(adjustl(labels(i)))//" ",i=1,nlab)
  do i=1,l
    write(100,'(a,3f17.8)') leftsymbol(i),(leftcoord(i,j),j=1,3)
  enddo   
  do i=1,d
    write(100,'(a,3f17.8)') devicesymbol(i),(devicecoord(i,j),j=1,3)
  enddo
  do i=1,r
    write(100,'(a,3f17.8)') rightsymbol(i),(rightcoord(i,j),j=1,3)
  enddo
 close(100)

 open(unit=100,file='stitch_full_shifted.xyz')
 write(100,*) l+d+r
 write(100,*) (trim(adjustl(labels(i)))//" ",i=1,nlab)
  do i=1,l
    write(100,'(a,3f17.8)') leftshiftedsymbol(i),(leftshiftedcoord(i,j),j=1,3)
  enddo
  do i=1,d
    write(100,'(a,3f17.8)') devicesymbol(i),(devicecoord(i,j),j=1,3)
  enddo
  do i=1,r
    write(100,'(a,3f17.8)') rightshiftedsymbol(i),(rightshiftedcoord(i,j),j=1,3)
  enddo
 close(100)

!#####################################################################!

 !LEFT 
 k=0
 do i=1,d
  do j=1,l
   if( norm(devicecoord(i,:) - leftshiftedcoord(j,:)) < cutoff_same_site ) then
       write(*,*) ' POINT [x] IN DEVICE WILL BE CONNECTED TO L ', i
       write(*,*) ' POINT [x] IN L WILL BE A TARGET            ', j
       k=k+1
   endif
  enddo
 enddo

 write(*,*) 'THERE ARE [x] SOURCES IN D AND [x] TARGETS IN L : ', k

 if(allocated(sites_in_device))      deallocate(sites_in_device);      allocate(sites_in_device(k)) 
 if(allocated(sites_in_lr))          deallocate(sites_in_lr)    ;      allocate(sites_in_lr(k))
 if(allocated(device_in_LR_shifted)) deallocate(device_in_LR_shifted); allocate(device_in_LR_shifted(d))
 if(allocated(LR_in_LR_shifted))     deallocate(LR_in_LR_shifted) ;    allocate(LR_in_LR_shifted(l))
 write(*,*) 'building dictionary device/LR_shifted'
 device_in_LR_shifted=0;LR_in_LR_shifted=0
 k=0
 do i=1,d
  do j=1,l
   if( norm(devicecoord(i,:) - leftshiftedcoord(j,:)) < cutoff_same_site ) then
       if(trim(adjustl(devicesymbol(i)))/=trim(adjustl(leftshiftedsymbol(j))))then
         write(*,*) 'ERROR SITE MISMATCH RIGHT 1'
         write(*,*) trim(adjustl(devicesymbol(i))),trim(adjustl(leftshiftedsymbol(j)))
         stop
       endif
       k=k+1
       device_in_LR_shifted(i)=j
       sites_in_device(k)=i
       sites_in_lr(k)=j
   endif
  enddo
 enddo
 write(*,*) '======THERE ARE [x] SITES AT THE INTERFACE =======', k

 write(*,*) 'building dictionnary LR/LR_shifted'
 do i=1,l
  do j=1,l
   if( norm(leftcoord(i,:) - leftshiftedcoord(j,:)) < cutoff_same_site ) then
       LR_in_LR_shifted(i)=j
   endif
  enddo
 enddo

 do t1=1,k
 do t2=1,k
   i=sites_in_device(t1)
   j=sites_in_lr(t2)
   v1 =  leftcoord(j,:) - devicecoord(i,:)
   v2 =  leftshiftedcoord(LR_in_LR_shifted(j),:) - leftshiftedcoord( device_in_LR_shifted(i),: )
   if(norm(v1-v2)>cutoff_same_site)then
    write(*,*) 'ERROR'
    stop
   endif
   v2 =  leftcoord(LR_in_LR_shifted(j),:) - leftcoord( device_in_LR_shifted(i),: )
   if(norm(v1-v2)>cutoff_same_site)then
    write(*,*) 'ERROR2'
    stop
   endif
 enddo
 enddo

 open(unit=103,file='connections_D_L_device_sites')
 open(unit=102,file='connections_D_L_lead_sites')
 open(unit=100,file='connections_D_L')
 write(103,*) k
 do t2=1,k
  i=sites_in_device(t2)
  j=sites_in_lr(t2)
  if(trim(adjustl(devicesymbol(i)))/=trim(adjustl(leftsymbol(j))))then
    write(*,*) 'ERROR SITE MISMATCH RIGHT 1'
    write(*,*) trim(adjustl(devicesymbol(i))),trim(adjustl(leftsymbol(j)))
    stop 
  endif
  write(102,*) j
  write(103,*) i, orb_num_orbitals( devicesymbol(i) )
 enddo

 do t1=1,k
 do t2=1,k
  i=sites_in_device(t1)
  j=sites_in_lr(t2)
  s=device_in_LR_shifted(i)   
  t=LR_in_LR_shifted(j)       
  !(i,j) = device site - connected to - L site
  !(s,t) = equivalent sites in L, to get H_st = H_ij
  if(trim(adjustl(devicesymbol(i)))/=trim(adjustl(leftshiftedsymbol(s))))then
    write(*,*) 'ERROR SITE MISMATCH LEFT 1'
    write(*,*) trim(adjustl(devicesymbol(i))),trim(adjustl(leftshiftedsymbol(s)))
    stop
  endif
  if( trim(adjustl(leftshiftedsymbol(t))) /= trim(adjustl(leftsymbol(j)))  )then
    write(*,*) 'ERROR SITE MISMATCH LEFT 2'
    write(*,*) trim(adjustl(leftsymbol(j))),trim(adjustl(leftshiftedsymbol(t)))
    stop
  endif
  if( trim(adjustl(leftsymbol(t))) /= trim(adjustl(leftsymbol(j)))  )then
    write(*,*) 'ERROR SITE MISMATCH LEFT 3'
    write(*,*) trim(adjustl(leftsymbol(t))),trim(adjustl(leftsymbol(j)))
    stop
  endif
  write(100,*) t1,t2,s,t
 enddo  
 enddo
 close(100)
 close(102)
 close(103)

!#####################################################################!

 !RIGHT
 k=0
 do i=1,d
  do j=1,r
   if( norm(devicecoord(i,:) - rightshiftedcoord(j,:)) < cutoff_same_site ) then
       write(*,*) ' POINT [x] IN DEVICE WILL BE CONNECTED TO R ', i
       write(*,*) ' POINT [x] IN R WILL BE A TARGET            ', j
       k=k+1
   endif
  enddo
 enddo

 write(*,*) 'THERE ARE [x] SOURCES IN D AND [x] TARGETS IN R : ', k

 if(allocated(sites_in_device))      deallocate(sites_in_device); allocate(sites_in_device(k))
 if(allocated(sites_in_lr))          deallocate(sites_in_lr)    ; allocate(sites_in_lr(k))
 if(allocated(device_in_LR_shifted)) deallocate(device_in_LR_shifted); allocate(device_in_LR_shifted(d))
 if(allocated(LR_in_LR_shifted))     deallocate(LR_in_LR_shifted) ; allocate(LR_in_LR_shifted(r))
 write(*,*) 'building dictionary device/LR_shifted'
 device_in_LR_shifted=0;LR_in_LR_shifted=0
 k=0
 do i=1,d
  do j=1,r
   if( norm(devicecoord(i,:) - rightshiftedcoord(j,:)) < cutoff_same_site ) then
       if(trim(adjustl(devicesymbol(i)))/=trim(adjustl(rightshiftedsymbol(j))))then
         write(*,*) 'ERROR SITE MISMATCH RIGHT 1'
         write(*,*) trim(adjustl(devicesymbol(i))),trim(adjustl(rightshiftedsymbol(j)))
         stop
       endif
       k=k+1
       device_in_LR_shifted(i)=j
       sites_in_device(k)=i
       sites_in_lr(k)=j
   endif
  enddo
 enddo
 write(*,*) 'building dictionnary LR/LR_shifted'
 do i=1,r
  do j=1,r
   if( norm(rightcoord(i,:) - rightshiftedcoord(j,:)) < cutoff_same_site ) then
       LR_in_LR_shifted(i)=j
   endif
  enddo
 enddo

 do t1=1,k
 do t2=1,k
   i=sites_in_device(t1)
   j=sites_in_lr(t2)
   v1 =  rightcoord(j,:) - devicecoord(i,:)
   v2 =  rightshiftedcoord(LR_in_LR_shifted(j),:) - rightshiftedcoord( device_in_LR_shifted(i),: )
   if(norm(v1-v2)>cutoff_same_site)then
    write(*,*) 'ERROR in building R'
    write(*,*) 'difference : ', norm(v1-v2)
    write(*,*) 'position in device : ', devicecoord(i,:)
    write(*,*) 'position in lead   : ', rightcoord(j,:)
    stop
   endif
   v2 = rightcoord(LR_in_LR_shifted(j),:) - rightcoord( device_in_LR_shifted(i),: )
   if(norm(v1-v2)>cutoff_same_site)then
    write(*,*) 'ERROR2 in building R'
    stop
   endif
 enddo
 enddo

 open(unit=103,file='connections_D_R_device_sites')
 open(unit=102,file='connections_D_R_lead_sites')
 open(unit=100,file='connections_D_R')
 write(103,*) k
 do t2=1,k
  i=sites_in_device(t2)
  j=sites_in_lr(t2)
  if(trim(adjustl(devicesymbol(i)))/=trim(adjustl(rightsymbol(j))))then
    write(*,*) 'ERROR SITE MISMATCH RIGHT 1'
    write(*,*) trim(adjustl(devicesymbol(i))),trim(adjustl(rightsymbol(j)))
    stop
  endif
  write(102,*) j
  write(103,*) i, orb_num_orbitals( devicesymbol(i) )
 enddo

 do t1=1,k
 do t2=1,k
  i=sites_in_device(t1)
  j=sites_in_lr(t2)
  s=device_in_LR_shifted(i)
  t=LR_in_LR_shifted(j)

  if(trim(adjustl(devicesymbol(i)))/=trim(adjustl(rightshiftedsymbol(s))))then
    write(*,*) 'ERROR SITE MISMATCH RIGHT 1'
    write(*,*) trim(adjustl(devicesymbol(i))),trim(adjustl(rightshiftedsymbol(s)))
    stop
  endif
  if(trim(adjustl(devicesymbol(i)))/=trim(adjustl(rightsymbol(s))))then
    write(*,*) 'ERROR SITE MISMATCH RIGHT 1b'
    write(*,*) trim(adjustl(devicesymbol(i))),trim(adjustl(rightsymbol(s)))
    stop
  endif
  if(trim(adjustl(rightsymbol(j)))/=trim(adjustl(rightshiftedsymbol(t))))then
    write(*,*) 'ERROR SITE MISMATCH RIGHT 2'
    write(*,*) trim(adjustl(rightsymbol(j))),trim(adjustl(rightshiftedsymbol(t)))
    stop
  endif
  if( trim(adjustl(rightsymbol(t))) /= trim(adjustl(rightsymbol(j)))  )then
    write(*,*) 'ERROR SITE MISMATCH RIGHT 3'
    write(*,*) trim(adjustl(rightsymbol(t))),trim(adjustl(rightsymbol(j)))
    stop
  endif

  !(i,j) = device site - connected to - L site
  !(s,t) = equivalent sites in L, to get H_st = H_ij
  write(100,*) t1,t2,s,t
 enddo
 enddo
 close(100)
 close(102)
 close(103)

!#####################################################################!

contains

  integer(4) function orb_num_orbitals(orb_label)
  implicit none 
  character*(*) :: orb_label
  integer       :: i,j,k
    do i=1,size(orbitals)
     if( trim(adjustl(orb_label)) == labels(i) )then
       orb_num_orbitals= orbitals(i) ; return ;
     endif
    enddo
    write(*,*) 'could not find orbital'
    stop
  end function

end program
