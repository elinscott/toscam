module smart_pointers

 USE namelistmod 
 USE common_def

IMPLICIT NONE
PRIVATE



TYPE , public :: pointer__
 integer                   :: istart=-1
 integer                   :: itot
 integer                   :: address_var,address
 real(8), pointer          :: a1,a2(:),a3(:,:),a4(:,:,:)
 complex(8),pointer        :: b1,b2(:),b3(:,:),b4(:,:,:) 
 integer, pointer          :: c1,c2(:),c3(:,:),c4(:,:,:)
 integer                   :: pointing_type(2)    ! first index =a,b,c, second index=1,2,3,4
 integer                   :: memory_type         ! 0=pointing to VOID 
                                                  ! 1=pointing to variable
                                                  ! 2=allocated pointer
                                                  ! 3=pointing to pointer
 character(22)             :: memory_type_label(0:10)
END TYPE



TYPE :: array_pointer
   type(pointer__), pointer :: array_pointer
END TYPE


TYPE, public, extends(pointer__) :: pointer_
   type(array_pointer) ::  POINTING_TO_ME(50) 
END TYPE 


INTERFACE point_to
 MODULE PROCEDURE  &
    &    POINTER_TO_VAR_1c,POINTER_TO_VAR_2c,POINTER_TO_VAR_3c,POINTER_TO_VAR_4c,POINTER_TO_VAR_1b,POINTER_TO_VAR_2b, &
    &    POINTER_TO_VAR_3b,POINTER_TO_VAR_4b,POINTER_TO_VAR_1,POINTER_TO_VAR_2,POINTER_TO_VAR_3,POINTER_TO_VAR_4
END INTERFACE 


INTERFACE ASSIGNMENT(=)
 MODULE PROCEDURE  &
    &    COPY_VAR_1c,COPY_VAR_2c,COPY_VAR_3c,COPY_VAR_4c,COPY_VAR_1b,COPY_VAR_2b, &
    & COPY_VAR_3b,COPY_VAR_4b,COPY_VAR_1,COPY_VAR_2,COPY_VAR_3,COPY_VAR_4 !,POINTER_TO_POINTER_, need intent(in)
END INTERFACE 

INTERFACE ASSIGNMENT(=)
 MODULE PROCEDURE  & 
    &  COPY_P_TO_VAR1c,COPY_P_TO_VAR2c,COPY_P_TO_VAR3c,COPY_P_TO_VAR4c, &
    &  COPY_P_TO_VAR1a,COPY_P_TO_VAR2a,COPY_P_TO_VAR3a,COPY_P_TO_VAR4a, &
    &  COPY_P_TO_VAR1b,COPY_P_TO_VAR2b,COPY_P_TO_VAR3b,COPY_P_TO_VAR4b  
END INTERFACE 

INTERFACE NULLIFY
 MODULE PROCEDURE KILL_POINTER_
END INTERFACE

PUBLIC  :: display_info_pointer,point_to,nullify,test_smart_pointers 
INTEGER :: istati 



CONTAINS

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine test_smart_pointers
  real(8)             :: aa
  real(8),allocatable :: a(:,:),b(:,:)
  TYPE(pointer_)      :: pp,pp2

   write(*,*) '#########################################'
   write(*,*) '#########################################'
   write(*,*) '#########################################'
   write(*,*) '#########################################'

   write(*,*) 'TEST 1 : pointing to VAR '

   allocate(a(2,3)); a=1.; allocate(b(2,3)); 
   writE(*,*) ' P1 point to double a=1'

   call point_to(pp,a)

   write(*,*) ' info on p1 .....'

   call display_info_pointer(pp) 

   write(*,*) 'p2 --> p1'
   pp2=pp

   write(*,*) 'info on p1.....'
   call display_info_pointer(pp)

   write(*,*) 'info on p2....'
   call display_info_pointer(pp2)

   write(*,*) 'kill p1 ....'
   call nullify(pp)

   write(*,*) 'info on p2....'
   call display_info_pointer(pp2)

   write(*,*) 'kill p2...'
   call nullify(pp2)

   write(*,*) 'info on p2 .....'
   call display_info_pointer(pp2)

   write(*,*) '#########################################'
   write(*,*) '#########################################'
   write(*,*) '#########################################'
   write(*,*) '#########################################'

   write(*,*) 'TEST 2 : COPYING VAR, allocate pointer'

   write(*,*) ' P1 point to double a=1'

   pp=a

   write(*,*) ' info on p1 .....'
   call display_info_pointer(pp)
   write(*,*) 'p2 --> p1'
   pp2=pp
   write(*,*) 'info on p1.....'
   call display_info_pointer(pp)
   write(*,*) 'info on p2....'

   b=pp
   write(*,*) ' b=pp2, b is : ', b

   call display_info_pointer(pp2)
   write(*,*) 'kill p1 ....'
   call nullify(pp)
   write(*,*) 'info on p2....'
   call display_info_pointer(pp2)
   write(*,*) 'kill p2...'
   call nullify(pp2)
   write(*,*) 'info on p2 .....'
   call display_info_pointer(pp2)

  stop 

 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine display_info_pointer(PP2)
 TYPE(pointer_) :: PP2
 integer        :: i,j,k

 if(PP2%istart==-1)then
  write(*,*) 'THIS POINTER WAS NEVER USED OR INITIALIAZED'
  return
 endif

 write(*,*) '========================================='
 write(*,*) ' POINTER INFO  : '
 write(*,*) ' pointing to memory address             : ', PP2%address
 write(*,*) ' how many pointers pointing to this guy : ', PP2%itot

 if(PP2%itot>0)then
 write(*,*) ' ----- '
 do i=1,PP2%itot
  write(*,*) 'pointing to me from address : ' , PP2%POINTING_TO_ME(i)%array_pointer%address
 enddo
 write(*,*) ' ----- '
 endif

 write(*,*) ' type of pointer                        : ', PP2%memory_type_label(PP2%memory_type)
 write(*,*) ' data pointed by the pointer            : '
 SELECT CASE( PP2%pointing_type(1) )
  CASE(1)
   SELECT CASE( PP2%pointing_type(2) )
    CASE(1); write(*,*) PP2%a1
    CASE(2); write(*,*) PP2%a2
    CASE(3); write(*,*) PP2%a3
    CASE(4); write(*,*) PP2%a4
   END SELECT
  CASE(2)
   SELECT CASE( PP2%pointing_type(2) )
    CASE(1); write(*,*) PP2%b1
    CASE(2); write(*,*) PP2%b2
    CASE(3); write(*,*) PP2%b3
    CASE(4); write(*,*) PP2%b4
   END SELECT
  CASE(3)
   SELECT CASE( PP2%pointing_type(2) )
    CASE(1); write(*,*) PP2%c1
    CASE(2); write(*,*) PP2%c2
    CASE(3); write(*,*) PP2%c3
    CASE(4); write(*,*) PP2%c4
   END SELECT
 END SELECT
 write(*,*) '========================================='

 end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine COPY_P_TO_VAR1c(a,PP)
 implicit none
  TYPE(pointer_),intent(in) :: PP
  integer(4),intent(inout)  :: a
   if(PP%memory_type==0 .or. PP%pointing_type(1)/=3 .or. PP%pointing_type(2)/=1 ) then
     write(*,*) 'error not good format'; return
   endif
   a=PP%c1
 end subroutine

 subroutine COPY_P_TO_VAR2c(a,PP)
 implicit none
  TYPE(pointer_),intent(in) :: PP
  integer(4),intent(inout)  :: a(:)
   if(PP%memory_type==0 .or. PP%pointing_type(1)/=3 .or. PP%pointing_type(2)/=2 ) then
     write(*,*) 'error not good format'; return
   endif
   a=PP%c2
 end subroutine

 subroutine COPY_P_TO_VAR3c(a,PP)
 implicit none
  TYPE(pointer_),intent(in) :: PP
  integer(4),intent(inout)  :: a(:,:)
   if(PP%memory_type==0 .or. PP%pointing_type(1)/=3 .or. PP%pointing_type(2)/=3 ) then
     write(*,*) 'error not good format'; return
   endif
   a=PP%c3
 end subroutine

 subroutine COPY_P_TO_VAR4c(a,PP)
 implicit none
  TYPE(pointer_),intent(in) :: PP
  integer(4),intent(inout)  :: a(:,:,:)
   if(PP%memory_type==0 .or. PP%pointing_type(1)/=3 .or. PP%pointing_type(2)/=4 ) then
     write(*,*) 'error not good format'; return
   endif
   a=PP%c4
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine COPY_P_TO_VAR1b(a,PP)
 implicit none
  TYPE(pointer_),intent(in) :: PP
  complex(8),intent(inout)  :: a
   if(PP%memory_type==0 .or. PP%pointing_type(1)/=2 .or. PP%pointing_type(2)/=1 ) then
     write(*,*) 'error not good format'; return
   endif
   a=PP%b1
 end subroutine

 subroutine COPY_P_TO_VAR2b(a,PP)
 implicit none
  TYPE(pointer_),intent(in) :: PP
  complex(8),intent(inout)  :: a(:)
   if(PP%memory_type==0 .or. PP%pointing_type(1)/=2 .or. PP%pointing_type(2)/=2 ) then
     write(*,*) 'error not good format'; return
   endif
   a=PP%b2
 end subroutine

 subroutine COPY_P_TO_VAR3b(a,PP)
 implicit none
  TYPE(pointer_),intent(in) :: PP
  complex(8),intent(inout)  :: a(:,:)
   if(PP%memory_type==0 .or. PP%pointing_type(1)/=2 .or. PP%pointing_type(2)/=3 ) then
     write(*,*) 'error not good format'; return
   endif
   a=PP%b3
 end subroutine

 subroutine COPY_P_TO_VAR4b(a,PP)
 implicit none
  TYPE(pointer_),intent(in) :: PP
  complex(8),intent(inout)  :: a(:,:,:)
   if(PP%memory_type==0 .or. PP%pointing_type(1)/=2 .or. PP%pointing_type(2)/=4 ) then
     write(*,*) 'error not good format'; return
   endif
   a=PP%b4
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine COPY_P_TO_VAR1a(a,PP)
 implicit none
  TYPE(pointer_),intent(in) :: PP
   real(8),intent(inout)    :: a
   if(PP%memory_type==0 .or. PP%pointing_type(1)/=1 .or. PP%pointing_type(2)/=1 ) then
    write(*,*) 'error not good format'; return
   endif
   a=PP%a1
 end subroutine

 subroutine COPY_P_TO_VAR2a(a,PP)
 implicit none
  TYPE(pointer_),intent(in) :: PP
  real(8),intent(inout)     :: a(:)
   if(PP%memory_type==0 .or. PP%pointing_type(1)/=1 .or. PP%pointing_type(2)/=2 ) then
    write(*,*) 'error not good format'; return
   endif
   a=PP%a2
 end subroutine

 subroutine COPY_P_TO_VAR3a(a,PP)
 implicit none
  TYPE(pointer_),intent(in) :: PP
  real(8),intent(inout)     :: a(:,:)
   if(PP%memory_type==0 .or. PP%pointing_type(1)/=1 .or. PP%pointing_type(2)/=3 ) then
    write(*,*) 'error not good format'; return
   endif
   a=PP%a3
 end subroutine

 subroutine COPY_P_TO_VAR4a(a,PP)
 implicit none
  TYPE(pointer_),intent(in) :: PP
  real(8),intent(inout)     :: a(:,:,:)
   if(PP%memory_type==0 .or. PP%pointing_type(1)/=1 .or. PP%pointing_type(2)/=4 ) then
    write(*,*) 'error not good format'; return
   endif
   a=PP%a4
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
 
 subroutine COPY_VAR_1c(PP,a)
 implicit none
  TYPE(pointer_),intent(inout) :: PP
  integer(4),intent(in) :: a
   call init_pointer__(PP%pointer__); PP%memory_type = 2 ; PP%pointing_type(1) = 3 ; PP%pointing_type(2) = 1; PP%c1=a
  PP%address_var=get_address(a)
 end subroutine
 subroutine COPY_VAR_2c(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  integer(4),intent(in)     :: a(:)
   IF(PP%memory_type/=2.or.PP%pointing_type(1)/=3.or.PP%pointing_type(2)/=2) call nullify_(PP%pointer__)
   if(.not.associated(PP%c2)) allocate(PP%c2(size(a,1)))
   if(ANY(shape(a)-shape(PP%c2)>0)) then ; deallocate(PP%c2,STAT=istati) ; allocate(PP%c2(size(a,1))); endif
   call init_pointer__(PP%pointer__); PP%memory_type = 2 ; PP%pointing_type(1) = 3 ; PP%pointing_type(2) = 2; PP%c2=a
  PP%address_var=get_address(a)
 end subroutine
 subroutine COPY_VAR_3c(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  integer(4),intent(in) :: a(:,:)
   IF(PP%memory_type/=2.or.PP%pointing_type(1)/=3.or.PP%pointing_type(2)/=3) call nullify_(PP%pointer__)
   if(.not.associated(PP%c3)) allocate(PP%c3(size(a,1),size(a,2)))
   if(ANY(shape(a)-shape(PP%c3)>0)) then ; deallocate(PP%c3,STAT=istati) ; allocate(PP%c3(size(a,1),size(a,2))); endif
   call init_pointer__(PP%pointer__); PP%memory_type = 2 ; PP%pointing_type(1) = 3 ; PP%pointing_type(2) = 3; PP%c3=a
  PP%address_var=get_address(a(1,1))
 end subroutine
 subroutine COPY_VAR_4c(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  integer(4),intent(in) :: a(:,:,:)
   IF(PP%memory_type/=2.or.PP%pointing_type(1)/=3.or.PP%pointing_type(2)/=4) call nullify_(PP%pointer__)
   if(.not.associated(PP%c4)) allocate(PP%c4(size(a,1),size(a,2),size(a,3)))
   if(ANY(shape(a)-shape(PP%c4)>0)) then ; deallocate(PP%c4,STAT=istati) ; allocate(PP%c4(size(a,1),size(a,2),size(a,3))); endif
   call init_pointer__(PP%pointer__); PP%memory_type = 2 ; PP%pointing_type(1) = 3 ; PP%pointing_type(2) = 4; PP%c4=a
  PP%address_var=get_address(a(1,1,1))
 end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
 
 subroutine POINTER_TO_VAR_1c(PP,a)
 implicit none
  TYPE(pointer_),intent(inout):: PP
  integer(4),intent(in),target :: a
   call init_pointer__(PP%pointer__); PP%memory_type = 1 ; PP%pointing_type(1) = 3 ; PP%pointing_type(2) = 1; PP%c1=>a
  PP%address_var=get_address(a)
 end subroutine
 subroutine POINTER_TO_VAR_2c(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  integer(4),intent(in),target :: a(:)
   call init_pointer__(PP%pointer__); PP%memory_type = 1 ; PP%pointing_type(1) = 3 ; PP%pointing_type(2) = 2; PP%c2=>a
  PP%address_var=get_address(a)
 end subroutine
 subroutine POINTER_TO_VAR_3c(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  integer(4),intent(in),target :: a(:,:)
   call init_pointer__(PP%pointer__); PP%memory_type = 1 ; PP%pointing_type(1) = 3 ; PP%pointing_type(2) = 3; PP%c3=>a
  PP%address_var=get_address(a(1,1))
 end subroutine
 subroutine POINTER_TO_VAR_4c(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  integer(4),intent(in),target :: a(:,:,:)
   call init_pointer__(PP%pointer__); PP%memory_type = 1 ; PP%pointing_type(1) = 3 ; PP%pointing_type(2) = 4; PP%c4=>a
  PP%address_var=get_address(a(1,1,1))
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine COPY_VAR_1b(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  COMPLEX(8),intent(in) :: a
   call init_pointer__(PP%pointer__); PP%memory_type = 2 ; PP%pointing_type(1) = 2 ; PP%pointing_type(2) = 1; PP%b1=a
  PP%address_var=get_address(a)
 end subroutine
 subroutine COPY_VAR_2b(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  COMPLEX(8),intent(in) :: a(:)
   IF(PP%memory_type/=2.or.PP%pointing_type(1)/=2.or.PP%pointing_type(2)/=2) call nullify_(PP%pointer__)
   if(.not.associated(PP%b2)) allocate(PP%b2(size(a,1)))
   if(ANY(shape(a)-shape(PP%b2)>0)) then ; deallocate(PP%b2,STAT=istati) ;  allocate(PP%b2(size(a,1))) ; endif
   call init_pointer__(PP%pointer__); PP%memory_type = 2 ; PP%pointing_type(1) = 2 ; PP%pointing_type(2) = 2; PP%b2=a
  PP%address_var=get_address(a)
 end subroutine
 subroutine COPY_VAR_3b(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  COMPLEX(8),intent(in) :: a(:,:)
   IF(PP%memory_type/=2.or.PP%pointing_type(1)/=2.or.PP%pointing_type(2)/=3) call nullify_(PP%pointer__)
   if(.not.associated(PP%b3))  allocate(PP%b3(size(a,1),size(a,2)))
   if(ANY(shape(a)-shape(PP%b3)>0)) then ; deallocate(PP%b3,STAT=istati) ; allocate(PP%b3(size(a,1),size(a,2))); endif
   call init_pointer__(PP%pointer__); PP%memory_type = 2 ; PP%pointing_type(1) = 2 ; PP%pointing_type(2) = 3; PP%b3=a
  PP%address_var=get_address(a(1,1))
 end subroutine
 subroutine COPY_VAR_4b(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  COMPLEX(8),intent(in) :: a(:,:,:)
   IF(PP%memory_type/=2.or.PP%pointing_type(1)/=2.or.PP%pointing_type(2)/=4) call nullify_(PP%pointer__)
   if(.not.associated(PP%b4)) allocate(PP%b4(size(a,1),size(a,2),size(a,3))) 
   if(ANY(shape(a)-shape(PP%b4)>0)) then ; deallocate(PP%b4,STAT=istati) ; allocate(PP%b4(size(a,1),size(a,2),size(a,3))) ; endif
   call init_pointer__(PP%pointer__); PP%memory_type = 2 ; PP%pointing_type(1) = 2 ; PP%pointing_type(2) = 4; PP%b4=a
  PP%address_var=get_address(a(1,1,1))
 end subroutine
 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine POINTER_TO_VAR_1b(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  COMPLEX(8),intent(in),target :: a
   call init_pointer__(PP%pointer__); PP%memory_type = 1 ; PP%pointing_type(1) = 2 ; PP%pointing_type(2) = 1; PP%b1=>a
  PP%address_var=get_address(a)
 end subroutine
 subroutine POINTER_TO_VAR_2b(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  COMPLEX(8),intent(in),target :: a(:)
   call init_pointer__(PP%pointer__); PP%memory_type = 1 ; PP%pointing_type(1) = 2 ; PP%pointing_type(2) = 2; PP%b2=>a
  PP%address_var=get_address(a)
 end subroutine
 subroutine POINTER_TO_VAR_3b(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  COMPLEX(8),intent(in),target :: a(:,:)
   call init_pointer__(PP%pointer__); PP%memory_type = 1 ; PP%pointing_type(1) = 2 ; PP%pointing_type(2) = 3; PP%b3=>a
  PP%address_var=get_address(a(1,1))
 end subroutine
 subroutine POINTER_TO_VAR_4b(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  COMPLEX(8),intent(in),target :: a(:,:,:)
   call init_pointer__(PP%pointer__); PP%memory_type = 1 ; PP%pointing_type(1) = 2 ; PP%pointing_type(2) = 4; PP%b4=>a
  PP%address_var=get_address(a(1,1,1))
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine COPY_VAR_1(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  real(8),intent(in) :: a
   call init_pointer__(PP%pointer__); PP%memory_type = 2 ; PP%pointing_type(1) = 1 ; PP%pointing_type(2) = 1; PP%a1=a
  PP%address_var=get_address(a)
 end subroutine
 subroutine COPY_VAR_2(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  real(8),intent(in) :: a(:)
   IF(PP%memory_type/=2.or.PP%pointing_type(1)/=1.or.PP%pointing_type(2)/=2) call nullify_(PP%pointer__) 
   if(.not.associated(PP%a2)) allocate(PP%a2(size(a,1))) 
   if(ANY(shape(a)-shape(PP%a2)>0)) then ; deallocate(PP%a2,STAT=istati) ; allocate(PP%a2(size(a,1))) ; endif
   call init_pointer__(PP%pointer__); PP%memory_type = 2 ; PP%pointing_type(1) = 1 ; PP%pointing_type(2) = 2; PP%a2=a
  PP%address_var=get_address(a)
 end subroutine
 subroutine COPY_VAR_3(PP,a)
 implicit none
  TYPE(pointer_) , intent(inout):: PP
  real(8),intent(in) :: a(:,:)
   IF(PP%memory_type/=2.or.PP%pointing_type(1)/=1.or.PP%pointing_type(2)/=3) call nullify_(PP%pointer__)
   if(.not.associated(PP%a3)) allocate(PP%a3(size(a,1),size(a,2)))   
   if(ANY(shape(a)-shape(PP%a3)>0)) then ; deallocate(PP%a3,STAT=istati) ; allocate(PP%a3(size(a,1),size(a,2))) ; endif
   call init_pointer__(PP%pointer__); PP%memory_type = 2 ; PP%pointing_type(1) = 1 ; PP%pointing_type(2) = 3; PP%a3=a
  PP%address_var=get_address(a(1,1))
 end subroutine
 subroutine COPY_VAR_4(PP,a)
 implicit none
  TYPE(pointer_),intent(inout) :: PP
  real(8),intent(in) :: a(:,:,:)
   IF(PP%memory_type/=2.or.PP%pointing_type(1)/=1.or.PP%pointing_type(2)/=4) call nullify_(PP%pointer__)
   if(.not.associated(PP%a4)) allocate(PP%a4(size(a,1),size(a,2),size(a,3))) 
   if(ANY(shape(a)-shape(PP%a4)>0)) then ; deallocate(PP%a4,STAT=istati) ; allocate(PP%a4(size(a,1),size(a,2),size(a,3))) ; endif
   call init_pointer__(PP%pointer__); PP%memory_type = 2 ; PP%pointing_type(1) = 1 ; PP%pointing_type(2) = 4; PP%a4=a
  PP%address_var=get_address(a(1,1,1))
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine POINTER_TO_VAR_1(PP,a)
 implicit none
  TYPE(pointer_),intent(inout) :: PP
  real(8),intent(in),target :: a
   call init_pointer__(PP%pointer__); PP%memory_type = 1 ; PP%pointing_type(1) = 1 ; PP%pointing_type(2) = 1; PP%a1=>a
  PP%address_var=get_address(a)
 end subroutine
 subroutine POINTER_TO_VAR_2(PP,a)
 implicit none
  TYPE(pointer_),intent(inout) :: PP
  real(8),intent(in),target :: a(:)
   call init_pointer__(PP%pointer__); PP%memory_type = 1 ; PP%pointing_type(1) = 1 ; PP%pointing_type(2) = 2; PP%a2=>a
  PP%address_var=get_address(a)
 end subroutine
 subroutine POINTER_TO_VAR_3(PP,a)
 implicit none
  TYPE(pointer_),intent(inout) :: PP
  real(8),intent(in),target :: a(:,:)
   call init_pointer__(PP%pointer__); PP%memory_type = 1 ; PP%pointing_type(1) = 1 ; PP%pointing_type(2) = 3; PP%a3=>a
  PP%address_var=get_address(a(1,1))
 end subroutine
 subroutine POINTER_TO_VAR_4(PP,a)
 implicit none
  TYPE(pointer_),intent(inout) :: PP
  real(8),intent(in),target :: a(:,:,:)
   call init_pointer__(PP%pointer__); PP%memory_type = 1 ; PP%pointing_type(1) = 1 ; PP%pointing_type(2) = 4; PP%a4=>a
  PP%address_var=get_address(a(1,1,1))
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine KILL_POINTER_(PP)
 TYPE(pointer_) :: PP
 integer        :: i,j,k,l,m
  if(PP%itot>0)then
   do i=1,PP%itot
    if(PP%POINTING_TO_ME(i)%array_pointer%memory_type/=3) stop 'error kill pointer memory type not consistent'
    CALL NULLIFY_(PP%POINTING_TO_ME(i)%array_pointer)
   enddo
  endif
  call NULLIFY_(PP%pointer__); PP%itot=0
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine POINTER_TO_POINTER_(PPfrom,PPto)
 implicit none
  TYPE(pointer_),target,intent(inout)  :: PPfrom
  TYPE(pointer_),intent(inout)         :: PPto

  call init_pointer__(PPfrom%pointer__); call init_pointer__(PPto%pointer__)
  PPto%itot=PPto%itot+1

  if(PPto%itot>size(PPto%POINTING_TO_ME(:))) stop 'error POINTER_TO_POINTER size is too big'

  PPto%POINTING_TO_ME(PPto%itot)%array_pointer => PPfrom%pointer__
  call nullify_(PPfrom%pointer__) 

  PPfrom%memory_type   = 3 
  PPfrom%pointing_type = PPto%pointing_type
 
  call point_to_(PPfrom%pointer__,PPto%pointer__) 
 end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine init_pointer__(PP)
  TYPE(pointer__) :: PP
   if(PP%istart>0) return
   PP%memory_type           =  0
   PP%pointing_type         =  0
   PP%itot                  =  0
   PP%istart                =  1
   PP%address               =  get_address(PP%istart)
   PP%memory_type_label(0)  = 'pointer undefined'
   PP%memory_type_label(1)  = 'pointing to variable'
   PP%memory_type_label(2)  = 'allocated pointer'
   PP%memory_type_label(3)  = 'pointing to pointer'
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


 subroutine deallocate_memory_space_pointer_(PP)
 TYPE(pointer__) :: PP
 SELECT  CASE( PP%pointing_type(1) )
  CASE(1)
   SELECT  CASE( PP%pointing_type(2) )
    CASE(1); IF(ASSOCIATED(PP%a1)) DEALLOCATE(PP%a1)
    CASE(2); IF(ASSOCIATED(PP%a2)) DEALLOCATE(PP%a2)
    CASE(3); IF(ASSOCIATED(PP%a3)) DEALLOCATE(PP%a3)
    CASE(4); IF(ASSOCIATED(PP%a4)) DEALLOCATE(PP%a4)
   END SELECT
  CASE(2)
   SELECT  CASE( PP%pointing_type(2) )
    CASE(1); IF(ASSOCIATED(PP%b1)) DEALLOCATE(PP%b1)
    CASE(2); IF(ASSOCIATED(PP%b2)) DEALLOCATE(PP%b2)
    CASE(3); IF(ASSOCIATED(PP%b3)) DEALLOCATE(PP%b3)
    CASE(4); IF(ASSOCIATED(PP%b4)) DEALLOCATE(PP%b4)
   END SELECT
  CASE(3)
   SELECT  CASE( PP%pointing_type(2) )
    CASE(1); IF(ASSOCIATED(PP%c1)) DEALLOCATE(PP%c1)
    CASE(2); IF(ASSOCIATED(PP%c2)) DEALLOCATE(PP%c2)
    CASE(3); IF(ASSOCIATED(PP%c3)) DEALLOCATE(PP%c3)
    CASE(4); IF(ASSOCIATED(PP%c4)) DEALLOCATE(PP%c4)
   END SELECT
 END SELECT
  PP%memory_type         = 0
  PP%pointing_type       = 0
 end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


 subroutine point_to_(PP1,PP2)
 TYPE(pointer__) :: PP1,PP2
 SELECT CASE( PP2%pointing_type(1) )
  CASE(1)
   SELECT CASE( PP2%pointing_type(2) )
    CASE(1); PP1%a1=>PP2%a1
    CASE(2); PP1%a2=>PP2%a2
    CASE(3); PP1%a3=>PP2%a3
    CASE(4); PP1%a4=>PP2%a4
   END SELECT
  CASE(2)
   SELECT CASE( PP2%pointing_type(2) )
    CASE(1); PP1%b1=>PP2%b1
    CASE(2); PP1%b2=>PP2%b2
    CASE(3); PP1%b3=>PP2%b3
    CASE(4); PP1%b4=>PP2%b4
   END SELECT
  CASE(3)
   SELECT CASE( PP2%pointing_type(2) )
    CASE(1); PP1%c1=>PP2%c1
    CASE(2); PP1%c2=>PP2%c2
    CASE(3); PP1%c3=>PP2%c3
    CASE(4); PP1%c4=>PP2%c4
   END SELECT
 END SELECT
 end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


 subroutine NULLIFY_POINTER_(PP)
 TYPE(pointer__) :: PP


 if(PP%pointing_type(2)<1.or.PP%pointing_type(2)>4)then
  write(*,*) ' NULLIFY POINTER CASE NOT DEFINED , pointing type : ' , PP%pointing_type(1:2)
  stop 'critical'
 endif

 if(PP%memory_type==1.or.PP%memory_type==3)then
  SELECT CASE( PP%pointing_type(1) )
  CASE(1)
   SELECT CASE( PP%pointing_type(2) )
    CASE(1); NULLIFY(PP%a1)  
    CASE(2); NULLIFY(PP%a2)
    CASE(3); NULLIFY(PP%a3)
    CASE(4); NULLIFY(PP%a4)
   END SELECT   
  CASE(2)
   SELECT CASE( PP%pointing_type(2) )
    CASE(1); NULLIFY(PP%b1)
    CASE(2); NULLIFY(PP%b2)
    CASE(3); NULLIFY(PP%b3)
    CASE(4); NULLIFY(PP%b4)
   END SELECT
  CASE(3)
   SELECT CASE( PP%pointing_type(2) )
    CASE(1); NULLIFY(PP%c1)
    CASE(2); NULLIFY(PP%c2)
    CASE(3); NULLIFY(PP%c3)
    CASE(4); NULLIFY(PP%c4)
   END SELECT
  CASE DEFAULT 
    write(*,*) 'pointing type is : ', PP%pointing_type(1)
    STOP 'nullify pointer case not defined'
  END SELECT
  PP%memory_type         = 0
  PP%pointing_type       = 0
 else
  stop 'error NULLIFY pointer and undefined memory zone'
 endif

 end subroutine


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine nullify_(PP)
  TYPE(pointer__) :: PP
  !(1)='pointing to variable'
  !(2)='allocated pointer'
  !(3)='pointing to pointer'

  SELECT CASE(PP%memory_type)
  CASE(0)
  CASE(1)
    CALL NULLIFY_POINTER_(PP)
  CASE(2)
    call deallocate_memory_space_pointer_(PP)
  CASE(3)
    CALL NULLIFY_POINTER_(PP)
  CASE DEFAULT
   stop 'nullify_pointer case not defined'
  END SELECT

 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************



end module

