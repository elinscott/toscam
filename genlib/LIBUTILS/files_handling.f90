module file_handling

 use genvar
 use splines
 use common_def, only : close_safe,open_safe 

INTERFACE ASSIGNMENT (=)
 MODULE PROCEDURE  read_file_in_array_comp_,read_file_in_array_, read_file_in_array__
END INTERFACE

contains

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 integer function check_size_array_in_file(unit,imax)
 implicit none
 integer :: i,j,k,l,m,imax,unit
 real(8) :: x(imax)

  i=0
  do 
   i=i+1
   if(i>imax) goto 51
   rewind(unit)
   read(unit,end=51,err=51) x(1:i)
  enddo
51 continue

  i=i-1
  check_size_array_in_file=i

 end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine read_file_put_in_array_and_resize(filename,array)
 implicit none
 character*(*)                      :: filename
 real(8)                            :: array(:,:)
 real(8),dimension(:,:),allocatable :: temp
 integer                            :: i,j,k,l,m,s1,s2,unit_
 logical                            :: check

     s1=size(array(:,1))
     s2=size(array(1,:))
     INQUIRE(file=TRIM(ADJUSTL(filename)),EXIST=check)
     if(.not.check)then
      array=0.;return;
     endif
     CALL open_safe(unit_,TRIM(ADJUSTL(filename)),"UNKNOWN","READ",get_unit=.true.)

     k=0
     do
      read(unit_,*,end=30)
      k=k+1
     enddo
     30 continue
     if(allocated(temp))deallocate(temp)
     allocate(temp(k,s2))
     rewind(unit_)
     do i=1,k
      read(unit_,*,err=10,end=10)  (temp(i,l),l=1,s2)
     enddo
     10 continue

     do l=2,s2
      call resampleit(temp(:,1),temp(:,l),array(:,1),array(:,l),0.d0)
     enddo

     if(allocated(temp))deallocate(temp)

     call close_safe(unit_)

 return
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine read_file_in_array_(temp,filename)
 implicit none
  real(8),intent(inout)    :: temp(:,:)
  character*(*),intent(in) :: filename
  integer                  :: i,j,k,l,s2,npoints,unit_
  logical                  :: check

  INQUIRE(file=TRIM(ADJUSTL(filename)),EXIST=check)
   if(.not.check)then
   temp=0.;return;
  endif

  CALL open_safe(unit_,TRIM(ADJUSTL(filename)),"UNKNOWN","READ",get_unit=.true.)

  npoints=size(temp(:,1))
  s2=size(temp(1,:))

  do i=1,npoints
   read(unit_,*,err=30,end=30) (temp(i,j),j=1,s2)
  enddo
  30 continue

  call close_safe(unit_)

 return
 end subroutine

subroutine read_file_in_array_comp_(temp,filename)
 implicit none
  complex(8),intent(inout) :: temp(:,:)
  character*(*),intent(in) :: filename
  integer                  :: i,j,k,l,s2,npoints,unit_
  logical                  :: check

  INQUIRE(file=TRIM(ADJUSTL(filename)),EXIST=check)
   if(.not.check)then
   temp=0.;return;
  endif

  CALL open_safe(unit_,TRIM(ADJUSTL(filename)),"UNKNOWN","READ",get_unit=.true.)

  npoints=size(temp(:,1))
  s2=size(temp(1,:))

  do i=1,npoints
   read(unit_,*,err=30,end=30) (temp(i,j),j=1,s2)
  enddo
  30 continue

  call close_safe(unit_)

 return
 end subroutine


!**************************************************************************
!**************************************************************************

 subroutine read_file_in_array__(temp,filename)
 implicit none
  real,intent(inout)       :: temp(:,:)
!BUG September 6th,2012 : for compilation with gfortran
! character*(*)            :: filename
  character*(*),intent(in) :: filename
  integer                  :: i,j,k,l,s2,npoints,unit_
  logical                  :: check

  INQUIRE(file=TRIM(ADJUSTL(filename)),EXIST=check)
  if(.not.check)then
    temp=0.;return;
  endif

  CALL open_safe(unit_,TRIM(ADJUSTL(filename)),"UNKNOWN","READ",get_unit=.true.)

  npoints=size(temp(:,1))
  s2=size(temp(1,:))

  do i=1,npoints
   read(unit_,*,err=30,end=30) (temp(i,j),j=1,s2)
  enddo
  30 continue

  call close_safe(unit_)
 return
 end subroutine

!**************************************************************************
!**************************************************************************
 
 integer function length_of_file(filename)
 implicit none
 character*(*)                     :: filename
 integer                           :: i,j,k,npoints,s,unit_

  CALL open_safe(unit_,TRIM(ADJUSTL(filename)),"UNKNOWN","READ",get_unit=.true.)

  npoints=0
  do 
   read(unit_,*,err=30,end=30)
   npoints=npoints+1
  enddo
  30 continue
  length_of_file=npoints
  call close_safe(unit_)

 return
 end function

!**************************************************************************
!**************************************************************************
 
 subroutine write_end_file(ifile)
 implicit none
 integer :: ifile,i,j
  do i=1,10000
   write(ifile,err=10) i,pi
  enddo
  10 continue
  call close_safe(ifile)

 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

end module
