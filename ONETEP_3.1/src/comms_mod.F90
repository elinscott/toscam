! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                       Communications module                    !
!----------------------------------------------------------------!
! This version uses the MPI library version 1.1                  !
!----------------------------------------------------------------!
! Written by Peter Haynes, 16/7/02                               !
!================================================================!

module comms

  use constants, only : DP, stdout

!CW
!#ifdef GPU_SPEEDUP
#if defined (GPU_SPEEDUP_WRAPPER) || defined (FFTW3GPU)
 use fortran_cuda
#endif
!END CW

  implicit none

  private

#ifdef MPI
#include "mpif.h"
#else
  integer, parameter :: MPI_COMM_NULL = 0
  integer, parameter :: MPI_COMM_WORLD = -1
  integer, parameter :: MPI_ANY_TAG = -1
#endif

  ! Information about parallelisation
  logical, public :: pub_comms_initialised = .false. ! Initialisation flag
  integer, public :: pub_total_num_nodes  ! Total number of nodes
  integer, public :: pub_my_node_id       ! My node ID
  integer, public :: pub_root_node_id     ! Root node ID
  logical, public :: pub_on_root          ! True only if my node is the root

  ! Information about groups of nodes
  logical, public :: pub_comms_groups_initialised
  integer, public :: pub_comms_group_size
  integer, public :: pub_num_comms_groups
  integer, public :: pub_my_comms_group
  integer, public :: pub_my_rank_in_group
  integer, public :: pub_first_node_in_group
  integer, public :: pub_last_node_in_group
  integer, public :: pub_my_rank_in_rank

  ! Publically accessible communicators
  integer, public :: pub_world_comm = MPI_COMM_NULL
  integer, public :: pub_group_comm = MPI_COMM_NULL
  integer, public :: pub_rank_comm = MPI_COMM_NULL

  ! Define null handle
#ifdef MPI
  integer, parameter, public :: pub_null_handle = MPI_REQUEST_NULL
#else
  integer, parameter, public :: pub_null_handle = 0
#endif

  ! jd: Define null handle
#ifdef MPI
  integer, parameter, public :: pub_null_proc = MPI_PROC_NULL
#else
  integer, parameter, public :: pub_null_proc = 0
#endif

  ! Public subroutines

  public :: comms_init
  public :: comms_groups_init
  public :: comms_free
  public :: comms_exit
  public :: comms_abort
  public :: comms_send
  public :: comms_recv
  public :: comms_irecv
  public :: comms_bcast
  public :: comms_reduce
  public :: comms_alltoall
  public :: comms_allgather
  public :: comms_barrier
  public :: comms_probe
  public :: comms_test
  public :: comms_testsome
  public :: comms_waitany
  public :: comms_wait

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  ! Interfaces for overloaded subroutines

  interface comms_send

     module procedure comms_send_integer_0
     module procedure comms_send_integer_1
     module procedure comms_send_integer_2
     module procedure comms_send_integer_3

     module procedure comms_send_real_0
     module procedure comms_send_real_1
     module procedure comms_send_real_2
     module procedure comms_send_real_3
     module procedure comms_send_real_4

     module procedure comms_send_complex_0
     module procedure comms_send_complex_1
     module procedure comms_send_complex_2
     module procedure comms_send_complex_3

     module procedure comms_send_character_0

     module procedure comms_send_logical_0
     module procedure comms_send_logical_1

  end interface


  interface comms_recv

     module procedure comms_recv_integer_0
     module procedure comms_recv_integer_1
     module procedure comms_recv_integer_2
     module procedure comms_recv_integer_3

     module procedure comms_recv_real_0
     module procedure comms_recv_real_1
     module procedure comms_recv_real_2
     module procedure comms_recv_real_3
     module procedure comms_recv_real_4

     module procedure comms_recv_complex_0
     module procedure comms_recv_complex_1
     module procedure comms_recv_complex_2
     module procedure comms_recv_complex_3

     module procedure comms_recv_character_0

     module procedure comms_recv_logical_0
     module procedure comms_recv_logical_1

  end interface

  interface comms_irecv

     module procedure comms_irecv_integer_0
     module procedure comms_irecv_integer_1
     module procedure comms_irecv_integer_2
     module procedure comms_irecv_integer_3

     module procedure comms_irecv_real_0
     module procedure comms_irecv_real_1
     module procedure comms_irecv_real_2
     module procedure comms_irecv_real_3
     module procedure comms_irecv_real_4

     module procedure comms_irecv_complex_0
     module procedure comms_irecv_complex_1
     module procedure comms_irecv_complex_2
     module procedure comms_irecv_complex_3

     module procedure comms_recv_character_0

     module procedure comms_irecv_logical_0
     module procedure comms_irecv_logical_1

  end interface


  interface comms_bcast

     module procedure comms_bcast_integer_0
     module procedure comms_bcast_integer_1
     module procedure comms_bcast_integer_2
     module procedure comms_bcast_integer_3

     module procedure comms_bcast_real_0
     module procedure comms_bcast_real_1
     module procedure comms_bcast_real_2
     module procedure comms_bcast_real_3

     module procedure comms_bcast_complex_0
     module procedure comms_bcast_complex_1
     module procedure comms_bcast_complex_2
     module procedure comms_bcast_complex_3

     module procedure comms_bcast_character_0
     module procedure comms_bcast_character_1

     module procedure comms_bcast_logical_0
     module procedure comms_bcast_logical_1

  end interface


  interface comms_reduce

     module procedure comms_reduce_integer_0
     module procedure comms_reduce_integer_1
     module procedure comms_reduce_integer_2
     module procedure comms_reduce_integer_3

     module procedure comms_reduce_real_0
     module procedure comms_reduce_real_1
     module procedure comms_reduce_real_2
     module procedure comms_reduce_real_3

     module procedure comms_reduce_complex_0
     module procedure comms_reduce_complex_1
     module procedure comms_reduce_complex_2
     module procedure comms_reduce_complex_3

     module procedure comms_reduce_logical_0
     module procedure comms_reduce_logical_1

     module procedure comms_reduce_integer_long_0

  end interface


  interface comms_alltoall

     module procedure comms_alltoall_integer_1
     module procedure comms_alltoall_integer_2
     module procedure comms_alltoall_integer_3

     module procedure comms_alltoall_real_1
     module procedure comms_alltoall_real_2
     module procedure comms_alltoall_real_3

     module procedure comms_alltoall_complex_1
     module procedure comms_alltoall_complex_2
     module procedure comms_alltoall_complex_3

     module procedure comms_alltoall_logical_1

  end interface

  interface comms_allgather

     module procedure comms_allgather_integer_0
     module procedure comms_allgather_integer_1
     module procedure comms_allgather_integer_2
     module procedure comms_allgather_integer_3

     module procedure comms_allgather_real_0
     module procedure comms_allgather_real_1
     module procedure comms_allgather_real_2
     module procedure comms_allgather_real_3

     module procedure comms_allgather_complex_0
     module procedure comms_allgather_complex_1
     module procedure comms_allgather_complex_2
     module procedure comms_allgather_complex_3

!     module procedure comms_allgather_logical_1

  end interface

  interface comms_copy

     module procedure comms_copy_integer_2
     module procedure comms_copy_integer_3

     module procedure comms_copy_real_2
     module procedure comms_copy_real_3

     module procedure comms_copy_complex_2
     module procedure comms_copy_complex_3

  end interface

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  ! Private data for module

  ! Communicator defined here
  integer :: communicator = MPI_COMM_NULL

  ! Default tag for MPI send / receive pairs
  integer, parameter :: default_tag = 0

  ! Upper bound for tags
  integer :: tag_ub

  ! Types defined here for comms module
#ifdef MPI
  integer, parameter :: integer_type = MPI_INTEGER
  integer, parameter :: real_type = MPI_DOUBLE_PRECISION
  integer, parameter :: complex_type = MPI_DOUBLE_COMPLEX
  integer, parameter :: character_type = MPI_CHARACTER
  integer, parameter :: logical_type = MPI_LOGICAL
  integer            :: integer_kind_long_type ! jd: Initialised in comms_init
#endif

  ! Private workspace
  integer, allocatable, dimension(:) :: iwork1
  real(kind=DP), allocatable, dimension(:) :: dwork1
  complex(kind=DP), allocatable, dimension(:) :: zwork1
  logical, allocatable, dimension(:) :: lwork1

  ! Stack variables for handles from non-blocking routines
  integer, parameter :: max_handles = 1024
  integer :: num_send_handles
  integer, dimension(max_handles) :: send_stack
  logical, dimension(max_handles) :: stack_flags
  integer, dimension(max_handles) :: completed_list

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_init

    !=========================================================================!
    ! This subroutine initialises the comms library and sets the public       !
    ! variables pub_total_num_nodes, pub_my_node_id, pub_root_node_id and     !
    ! pub_on_root.                                                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! None                                                                    !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! pub_total_num_nodes   : On exit, set to number of nodes in              !
    !                           MPI_COMM_WORLD                                !
    !   pub_my_node_id        : On exit, the id of the node in global group   !
    !   pub_root_node_id      : On exit, set to the id of the root node       !
    !   pub_on_root           : On exit, indicates whether this node is root  !
    !   pub_comms_initialised : On exit, set to .true.                        !
    !   communicator          : The group of currently active nodes.          !
    !                           On exit, set to MPI_COMM_WORLD                !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error    : Used as an error flag for MPI routines                       !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    use constants, only: LONG_R

    implicit none

    ! Additional #ifdef statements have been introduced for Accelrys builds
    ! on Linux where we use a confidential license key from HP to enable
    ! HP-MPI mechanism.

!CW
    integer :: mygpu
    logical :: check
!END CW

#ifdef MPI
    ! Local variables
    integer :: error                             ! Error flag
    logical :: flag                              ! Attribute flag
#endif

    ! If already initialised, return immediately
    if (pub_comms_initialised) return

#ifdef MPI

    call HPMPIINIT
    ! Initialise MPI
    call MPI_INIT(error)
    if (error /= MPI_SUCCESS) then
       write(stdout,'(a,i3)') 'Error in comms_init: &
            &MPI_INIT failed with code ',error
       stop
    end if
#endif

    ! Get the total number of nodes
#ifdef MPI
    call MPI_COMM_SIZE(MPI_COMM_WORLD,pub_total_num_nodes,error)
    if (error /= MPI_SUCCESS) then
       write(stdout,'(a,i3)') &
            'Error in comms_init: MPI_COMM_SIZE failed with code ',error
       call MPI_ABORT(MPI_COMM_WORLD,1,error)
    end if
#else
    pub_total_num_nodes = 1
#endif

    ! Get the id of the local node
#ifdef MPI
    call MPI_COMM_RANK(MPI_COMM_WORLD,pub_my_node_id,error)
    if (error /= MPI_SUCCESS) then
       write(stdout,'(a,i3)') &
            'Error in comms_init: MPI_COMM_RANK failed with code ',error
       call MPI_ABORT(MPI_COMM_WORLD,1,error)
    end if
#else
    pub_my_node_id = 0
#endif

!CW
#if defined (GPU_SPEEDUP_WRAPPER) || defined (FFTW3GPU)
  write(*,*) 'bringing online GPUs'
  call init_gpu_device
  INQUIRE(file='my_gpu',EXIST=check)
  if(.not.check)then
   call choose_gpu(pub_my_node_id)
   write(*,*) 'CPU [x] is choosing GPU [y] : ', pub_my_node_id,pub_my_node_id
  else
   open(unit=1326,file='my_gpu')
   read(1326,*) mygpu
   close(1326)
   if(pub_my_node_id/=0)then
     write(*,*) 'file mygpu is present, but it is not compatible with mpi'
     stop
   endif
   write(*,*) '-----------> USING GPU NUMBER (first start at 0): ',mygpu
   call choose_gpu(mygpu)
  endif
#endif
!END CW

    ! Choose node 0 as the root node
    pub_root_node_id = 0

    ! Set flag to indicate whether local node is root node
    pub_on_root = (pub_my_node_id == pub_root_node_id)

    ! Get the upper bound on tags
#ifdef MPI
    call MPI_ATTR_GET(MPI_COMM_WORLD,MPI_TAG_UB,tag_ub,flag,error)
    if (error /= MPI_SUCCESS) then
       write(stdout,'(a,i3)') &
            'Error in comms_init: MPI_COMM_GET_ATTR failed with code ',error
       call MPI_ABORT(MPI_COMM_WORLD,1,error)
    end if
    if (.not. flag) tag_ub = 32767
#endif

    ! Set the current communicator of active nodes to MPI_COMM_WORLD
    communicator = MPI_COMM_WORLD

    ! Initialise the handle stack for non-blocking routines
    num_send_handles = 0

    ! Set flag to indicate comms module has been initialised
    pub_comms_initialised = .true.

#ifdef MPI
    ! Perform a barrier synchronisation before returning
    call MPI_BARRIER(MPI_COMM_WORLD,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_init: MPI_BARRIER failed with code ',error
       call MPI_ABORT(MPI_COMM_WORLD,1,error)
    end if
#endif

  end subroutine comms_init

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_groups_init()

    !=========================================================================!
    ! This subroutine initialises the communicator for a group of nodes
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! None                                                                    !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_group_size    : On exit, set to number of nodes in          !
    !                             pub_group_comm                              !
    !   pub_my_rank_in_group    : On exit, the id of the node in local group  !
    !   pub_my_comms_group      : On exit, the group number of the local node !
    !   pub_num_comms_groups    : On exit, the number of comms groups         !
    !   pub_first_node_in_group : On exit, the first node in this node's group!
    !   pub_last_node_in_group  : On exit, the last node in this node's group !
    !   pub_comms_groups_initialised : On exit, set to .true.                 !
    !   pub_group_comm          : The communicator for the local nodes.       !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error    : Used as an error flag for MPI routines                       !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 03/11/2010.                                   !
    !=========================================================================!

    implicit none

    ! Local variables
    integer :: ierr
    integer :: node
#ifdef MPI
    integer :: group_world
    integer :: group
#endif
    integer, allocatable, dimension(:) :: nodes_in_group

    ! If already initialised, return immediately
    if (pub_comms_groups_initialised) return

    pub_comms_group_size = min(pub_total_num_nodes,pub_comms_group_size)

    ! Check the group size not greater than the total number of nodes
    if (pub_total_num_nodes<pub_comms_group_size) then
       if (pub_on_root) write(stdout,'(a)') 'Error in comms_groups_init: &
            &Comms group size is greater than total number of nodes'
       call comms_abort(.true.)
    end if
    ! Check the group size is a divisor of the total number of nodes
    if (modulo(pub_total_num_nodes,pub_comms_group_size)/=0) then
       if (pub_on_root) write(stdout,'(a)') 'Error in comms_groups_init: &
            &Comms group size not a divisor of number of nodes'
       call comms_abort(.true.)
    end if
    pub_num_comms_groups = pub_total_num_nodes / pub_comms_group_size

    ! Find out the group number of the local node
    pub_my_comms_group = int(pub_my_node_id / pub_comms_group_size)

    ! Create a list of the nodes in the local group
    allocate(nodes_in_group(pub_comms_group_size),stat=ierr)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_groups_init: &
            &allocating nodes_in_group failed with code ',ierr
       call comms_abort
    end if

    do node=1,pub_comms_group_size
       nodes_in_group(node) = pub_my_comms_group*pub_comms_group_size &
            + node - 1
    end do
    pub_first_node_in_group = nodes_in_group(1)
    pub_last_node_in_group = nodes_in_group(pub_comms_group_size)

#ifdef MPI
    ! Obtain the group reference for MPI_COMM_WORLD
    pub_world_comm = communicator
    call MPI_COMM_GROUP(pub_world_comm,group_world,ierr)

    ! Create the local group
    call MPI_GROUP_INCL(group_world,pub_comms_group_size, &
         nodes_in_group,group,ierr)

    ! Find information about the this node and the local group
    call MPI_GROUP_RANK(group,pub_my_rank_in_group,ierr)
    call MPI_GROUP_SIZE(group,pub_comms_group_size,ierr)

    ! Create a communicator for the group
    call MPI_COMM_CREATE(pub_world_comm,group,pub_group_comm,ierr)
#else
    pub_world_comm = communicator
    pub_group_comm = communicator
    pub_my_rank_in_group = 0
#endif

    ! Deallocate groups list
    deallocate(nodes_in_group,stat=ierr)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_groups_init: &
            &deallocating nodes_in_group failed with code ',ierr
       call comms_abort
    end if

    ! Create a list of the nodes in the rank group
    allocate(nodes_in_group(pub_num_comms_groups),stat=ierr)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_groups_init: &
            &allocating nodes_in_group failed with code ',ierr
       call comms_abort
    end if

    do node=1,pub_num_comms_groups
       nodes_in_group(node) = pub_my_rank_in_group + &
            pub_comms_group_size*(node-1)
    end do

#ifdef MPI
    ! Create the local group
    call MPI_GROUP_INCL(group_world,pub_num_comms_groups, &
         nodes_in_group,group,ierr)

    ! Find information about the this node and the local group
    call MPI_GROUP_RANK(group,pub_my_rank_in_rank,ierr)
    call MPI_GROUP_SIZE(group,pub_num_comms_groups,ierr)

    ! Create a communicator for the group
    call MPI_COMM_CREATE(pub_world_comm,group,pub_rank_comm,ierr)
#else
    pub_rank_comm = communicator
    pub_my_rank_in_rank = 0
#endif

    ! Deallocate groups list
    deallocate(nodes_in_group,stat=ierr)
    if (ierr /= 0) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_groups_init: &
            &deallocating nodes_in_group failed with code ',ierr
       call comms_abort
    end if

    ! Record completion of initialisation
    pub_comms_groups_initialised = .true.

  end subroutine comms_groups_init

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_probe(result,node,tag,comm)

    !=========================================================================!
    ! This subroutine probes to see if a message is waiting from a given node !
    ! with a given tag.                                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! result (output)      : True if there is a message waiting matching      !
    !                          the given criteria, false if there is not      !
    ! node (input)         : The node to check for a message from             !
    !   tag (output)         : Tag of the message required from node          !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error   : Used as an error flag for deallocate                          !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine 28/10/08                                       !
    !=========================================================================!

    implicit none

    ! Arguments
    logical,intent(out) :: result
    integer,intent(in) :: node
    integer,intent(in),optional :: tag
    integer,intent(in),optional :: comm

#ifdef MPI
    ! Local variables
    integer :: error ! Error Flag
    integer :: status(MPI_STATUS_SIZE)
    integer :: local_tag
    integer :: local_comm

    ! Make local copy of optional arguments
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Call MPI_IPROBE routine to execute test for incoming message
    call MPI_IPROBE(node,local_tag,local_comm,result,status,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i10)') &
            'Error in comms_probe: MPI_IPROBE failed with code ',error
       call comms_abort
    end if
#else
    result = .false.
    if (present(tag) .and. present(comm) .and. node == 1 ) return
#endif

  end subroutine comms_probe

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_test(result,handle)

    !=========================================================================!
    ! This subroutine test to see if a message with a given handle has        !
    ! completed yet.                                                          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! result (output)      : True if the message with the given handle has    !
    !                          completed, false if it still exists            !
    !   handle (input)       : Handle to test                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error   : Used as an error flag for deallocate                          !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine 05/08/09                                       !
    !=========================================================================!

    implicit none

    ! Arguments
    logical, intent(out) :: result
    integer, intent(in) :: handle

#ifdef MPI
    ! Local variables
    integer :: error ! Error Flag
    integer :: status(MPI_STATUS_SIZE)

    ! Call MPI_TEST routine to execute test for incoming message
    call MPI_TEST(handle,result,status,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i10)') &
            'Error in comms_test: MPI_TEST failed with code ',error
       call comms_abort
    end if
#else
    if (handle<0) result = .false.
    result = .true.
#endif

  end subroutine comms_test

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_testsome(nhandle,handles,nfreed,freed)

    !=========================================================================!
    ! This subroutine test to see if any message with one of a list of        !
    ! handle has completed yet.                                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! nhandles (input) : Number of handles to test                            !
    ! handles (input)  : A list of length nhandles to test                    !
    ! nfreed (output)  : Number of the handles asked about which have been    !
    !                    freed.                                               !
    ! freed (output)   : A list of length nfreed (though allocated to size    !
    !                    nhandles as this is the biggest it could be) giving  !
    !                    the indices of the handles which have been freed.    !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error   : Used as an error flag for deallocate                          !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine 05/08/09                                       !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: nhandle
    integer, intent(in) :: handles(nhandle)
    integer, intent(out) :: nfreed
    integer, intent(out) :: freed(nhandle)

#ifdef MPI
    ! Local variables
    integer :: error ! Error Flag
    integer :: status(MPI_STATUS_SIZE)

    ! Call MPI_TESTSOME routine to list completed operations
    call MPI_TESTSOME(nhandle,handles,nfreed,freed,status,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i10)') &
            'Error in comms_testsome: MPI_TESTSOME failed with code ',error
       call comms_abort
    end if
#else
    if ((nhandle < 0).or.(all(handles<0))) nfreed = 0
    nfreed = 0
    freed = 0
#endif

  end subroutine comms_testsome

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_waitany(count,handle_array,index)

    !=========================================================================!
    ! This subroutine waits for any of an array of handles passed in to       !
    ! complete, and then returns the index of the handle that finished first  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! count (input)        : Number of handles in handle_array                !
    ! handle_array (input) : An array of handles previously returned by       !
    !                          send or irecv calls                            !
    !   index (output)       : The index in handle_array of the handle that   !
    !                          returned (NB: several may have returned        !
    !                          simulateneously - do not rely on trapping them !
    !                          all as return values from this routine!)       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error   : Used as an error flag for deallocate                          !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine 06/11/08                                       !
    !=========================================================================!

    implicit none

    ! Arguments
    integer,intent(in) :: count
    integer,intent(in) :: handle_array(1:count)
    integer,intent(out),optional :: index

#ifdef MPI
    ! Local variables
    integer :: error ! Error Flag
    integer :: status(MPI_STATUS_SIZE)
    integer :: local_index

    ! Call MPI_WAITANY routine to wait for completed operation
    call MPI_WAITANY(count,handle_array,local_index,status,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i10)') &
            'Error in comms_waitany: MPI_WAITANY failed with code ',error
       call comms_abort
    end if
    if (present(index)) then
       index = local_index
    end if
#else
    if (present(index)) index = 0
    if (size(handle_array) == count) return
#endif

  end subroutine comms_waitany

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_wait(handle,free_stack)

    !=========================================================================!
    ! This subroutine waits for mpi operation associated with the specified   !
    ! handle to complete. If the operation is a send request, free_stack can  !
    ! be set to .true. to remove this entry from the send stack.              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! handle (input)     : The handle of the send or receive to wait for      !
    !                        completion of                                    !
    !   free_stack (input) : If waiting for a send, set this to true to mark  !
    !                        the send as completed in the send stack          !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error   : Used as a return value for mpi_wait                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine 06/11/08                                       !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: handle
    logical, optional, intent(in) :: free_stack

    ! Local variables

#ifdef MPI
    integer :: status(MPI_STATUS_SIZE)
    integer :: error ! Error Flag
    integer :: istack,jstack
    logical :: found_handle
#endif
    logical :: local_free_stack

    if (present(free_stack)) then
       local_free_stack = free_stack
    else
       local_free_stack = .false.
    end if

#ifdef MPI

    if (local_free_stack .and. (handle /= pub_null_handle)) then

       ! Find this handle in the stack
       found_handle = .false.
       do istack=1,num_send_handles
         if (handle==send_stack(istack)) then
            found_handle = .true.
            exit
         end if
       end do

       ! If the handle was not found in the send stack, then
       ! the send must already have completed, so return
       if (.not.found_handle) return

       ! Call MPI_WAIT routine to wait for send to complete
       call MPI_WAIT(send_stack(istack),status,error)

       ! Shuffle the rest of the handles down
       do jstack=istack,num_send_handles-1
          send_stack(jstack) = send_stack(jstack+1)
       end do
       num_send_handles = num_send_handles - 1
    else
       ! Call MPI_WAIT routine to wait for message to complete
       call MPI_WAIT(handle,status,error)
    end if

    if (error /= MPI_SUCCESS) then
       !if (pub_on_root)
       write(stdout,'(a,i10)') &
            'Error in comms_wait: MPI_WAIT failed with code ',error
       call comms_abort
    end if
#else
    if (handle == -1 ) return
#endif

  end subroutine comms_wait

!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_free

    !=========================================================================!
    ! This subroutine frees up workspace used by the comms library.           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! None                                                                    !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised                                                 !
    !   iwork1, iwork2, iwork3, dwork1 etc.                                   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error   : Used as an error flag for deallocate                          !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Local variables
    integer :: error                   ! Error flag

    ! If the comms library is not initialised, return immediately
    if (.not. pub_comms_initialised) return

    if (allocated(iwork1)) then
       deallocate(iwork1,stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_free: &
               &deallocating iwork1 failed with code ',error
          call comms_abort
       end if
    end if

    if (allocated(dwork1)) then
       deallocate(dwork1,stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_free: &
               &deallocating dwork1 failed with code ',error
          call comms_abort
       end if
    end if

    if (allocated(zwork1)) then
       deallocate(zwork1,stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_free: &
               &deallocating zwork1 failed with code ',error
          call comms_abort
       end if
    end if

    if (allocated(lwork1)) then
       deallocate(lwork1,stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_free: &
               &deallocating lwork1 failed with code ',error
          call comms_abort
       end if
    end if

    ! Empty the stack of outstanding non-blocking sends
    call comms_empty_send_stack

  end subroutine comms_free

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_exit

    !=========================================================================!
    ! This subroutine finalises the comms library and halts execution. Used   !
    ! for normal termination of execution.                                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! None                                                                    !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : On exit, set to .false.                       !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error   : Used as an error flag for MPI routines                        !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

#ifdef MPI
    ! Local variables
    integer :: error                   ! Error flag
#endif

    ! If the comms library is not initialised, return immediately
    if (.not. pub_comms_initialised) return

    ! Free up any workspace
    call comms_free

    ! Free memory associated with group comms
    call comms_groups_exit

#ifdef MPI
    ! Call a barrier to synchronise all nodes
    call MPI_BARRIER(MPI_COMM_WORLD,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_exit: MPI_BARRIER failed with code ',error
       call MPI_ABORT(MPI_COMM_WORLD,1,error)
    end if

    ! Finalise MPI
    call MPI_FINALIZE(error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_exit: MPI_FINALIZE failed with code ',error
       call MPI_ABORT(MPI_COMM_WORLD,1,error)
    end if
#endif

    ! Set current communicator to MPI_COMM_NULL
    communicator = MPI_COMM_NULL

    ! Set flag to indicate comms module is no longer initialised
    pub_comms_initialised = .false.

  end subroutine comms_exit

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_groups_exit

    !=========================================================================!
    ! This subroutine deallocates the storage associated with group comms.    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! None                                                                    !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_groups_initialised : On exit, set to .false.                !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error   : Used as an error flag for MPI routines                        !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 06/11/10                                      !
    !=========================================================================!

    implicit none

#ifdef MPI
    ! Local variables
    integer :: error                   ! Error flag
#endif

    ! If the groups are not already initialised, return immediately
    if (.not. pub_comms_groups_initialised) return

    ! Free the communicators created in comms_groups_init
#ifdef MPI
    call MPI_COMM_FREE(pub_rank_comm,error)
    if (error /= 0) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_groups_exit: &
            &freeing rank communicator failed with code ',error
       call comms_abort
    end if
    call MPI_COMM_FREE(pub_group_comm,error)
    if (error /= 0) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_groups_exit: &
            &freeing group communicator failed with code ',error
       call comms_abort
    end if
#endif

    pub_comms_groups_initialised = .false.

  end subroutine comms_groups_exit

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_abort(require_all)

    !=========================================================================!
    ! This subroutine aborts execution on all nodes. Used for abnormal        !
    ! termination of execution.                                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   require_all (in, optional): If .true., will wait on a barrier before  !
    !                               aborting.                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   error : Used as an error flag for MPI routines                        !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    logical, optional, intent(in) :: require_all

#ifdef MPI
    ! Local variables
    integer :: error ! Error flag
#ifdef SUN
    integer :: flush ! Sun flush function
#endif
#ifndef INTFLUSH
    external flush
#endif
#endif

    ! If the comms library has been initialised, call MPI_ABORT. Otherwise
    ! simply halt execution.
    if (pub_comms_initialised) then
#ifdef MPI
       if (pub_on_root) then
          write(stdout,'(a)') 'ONETEP execution aborted'
#ifdef SUN
          error = flush(stdout)
#else
#ifndef SGI
!CW
#ifndef FLUSH_EXT_NAME
          call flush(stdout)
#else
          call flush_(stdout)
#endif
!END CW
#else
         ! SGI introduced a second argument starting from the v7.4
         ! of their runtime environment. This form is still acceptable
         ! by the older environments - the return value is ignored in those
         ! older cases, and the code exits on error.
         ! This second argument is also allowed under UNICOS as an
         ! optional argument.
          call flush(stdout,error)
#endif
#endif
       end if

    ! SGI implementation of MPI_ABORT is faulty: the bug is in the
    ! runtime Fortran environment, and it will not be fixed until 7.4 version
    ! of MIPS Pro compilers. The problem is that MPI_ABORT kills all
    ! parents of the process, which means MaterialsStudio Gateway is
    ! killed, too. A temporary solution is to store PIDs of all MPI
    ! processes (reportpid call in license.F90), create a killfile at this
    ! point, and trust the Perl script to terminate all these PIDs.

       if(present(require_all)) then
          if(require_all) call comms_barrier()
       end if

#ifdef ACCELRYS
       write (stdout,*) ' '
       call comms_killfile()
#ifdef SGI
       call MPI_FINALIZE(error)
       stop
#endif
#endif

       call MPI_ABORT(MPI_COMM_WORLD,1,error)
       return
#else
       stop 'ONETEP execution aborted'
#endif
    else
       stop 'ONETEP execution aborted'
    end if

  end subroutine comms_abort

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_barrier(comm)

    !=========================================================================!
    ! This subroutine performs a barrier synchronisation across all active    !
    ! nodes.                                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised   !
    !   communicator          : The communicator corresponding to active      !
    !                           nodes in the current distribution strategy    !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, optional, intent(in) :: comm

    ! Local variables
    integer :: local_comm
#ifdef MPI
    integer :: error ! Error flag
#endif

    ! Check comms have been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_barrier: comms not initialised'
       call comms_abort
    end if

    ! Check for optional argument
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

#ifdef MPI
    ! Call MPI_BARRIER routine to execute barrier synchronisation
    call MPI_BARRIER(local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_barrier: MPI_BARRIER failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_barrier

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !=========================================================================!
  !                                                                         !
  !      Private subroutines                                                !
  !                                                                         !
  !=========================================================================!


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_integer_0(node,i_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the integer scalar form of comms_send. It sends the  !
    ! integer i_array, or optionally the first length integers pointed to by  !
    ! i_array to node.                                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   i_array (input) : The first element of the array of data to be sent.  !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    integer, target, intent(in) :: i_array   ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    integer, optional, intent(in) :: comm    ! The communicator to use
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_integer_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_integer_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(i_array,local_length,integer_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_integer_0: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
    if (node == i_array) return
#endif

  end subroutine comms_send_integer_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_integer_1(node,i_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the integer vector form of comms_send. It sends the  !
    ! first length integers of i_array to node, or the whole array if length  !
    ! is absent.                                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   i_array (input) : The array of data to be sent.                       !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    integer, target, intent(in) :: i_array(:)  ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_integer_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_integer_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_integer_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(i_array(1),local_length,integer_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_integer_1: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
#endif

  end subroutine comms_send_integer_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_integer_2(node,i_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the integer matrix form of comms_send. It sends the  !
    ! first length integers of i_array to node, or the whole array if length  !
    ! is absent.                                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   i_array (input) : The array of data to be sent.                       !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    integer, target, intent(in) :: i_array(:,:)  ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_integer_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_integer_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_integer_2: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(i_array(1,1),local_length,integer_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_integer_2: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
    if (node == 1) return
#endif

  end subroutine comms_send_integer_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_integer_3(node,i_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the integer rank 3 tensor form of comms_send. It     !
    ! sends the first length integers of i_array to node, or the whole array  !
    ! if length is absent.                                                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   i_array (input) : The array of data to be sent.                       !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    integer, target, intent(in) :: i_array(:,:,:)  ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_integer_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_integer_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_integer_3: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(i_array(1,1,1),local_length,integer_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_integer_3: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
    if (node == 1) return
#endif

  end subroutine comms_send_integer_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_real_0(node,d_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the (double precision) real scalar form of           !
    ! comms_send. It sends the real d_array, or optionally the first length   !
    ! real pointed to by d_array to node.                                     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   d_array (input) : The first element of the array of data to be sent.  !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    real(kind=DP), target, intent(in) :: d_array  ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_real_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) 'Error in comms_send_real_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(d_array,local_length,real_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_real_0: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
#endif

  end subroutine comms_send_real_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_real_1(node,d_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the (double precision) real vector form of           !
    ! comms_send. It sends the first length reals of d_array, or the whole    !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   d_array (input) : The array of data to be sent.                       !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    real(kind=DP), target, intent(in) :: d_array(:)  ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_real_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_real_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_real_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(d_array(1),local_length,real_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_real_1: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
#endif

  end subroutine comms_send_real_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_real_2(node,d_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the (double precision) real matrix form of           !
    ! comms_send. It sends the first length reals of d_array, or the whole    !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   d_array (input) : The array of data to be sent.                       !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    real(kind=DP), target, intent(in) :: d_array(:,:)  ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_real_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_real_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_real_2: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(d_array(1,1),local_length,real_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_real_2: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle

#else
    if (present(return_handle)) return_handle = -1
    if (node == size(d_array)) return
#endif

  end subroutine comms_send_real_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_real_3(node,d_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the (double precision) real rank 3 tensor form of    !
    ! comms_send. It sends the first length reals of d_array, or the whole    !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   d_array (input) : The array of data to be sent.                       !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    real(kind=DP), target, intent(in) :: d_array(:,:,:)  ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_real_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_real_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_real_3: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(d_array(1,1,1),local_length,real_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_real_3: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
#endif

  end subroutine comms_send_real_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_real_4(node,d_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the (double precision) real rank 3 tensor form of    !
    ! comms_send. It sends the first length reals of d_array, or the whole    !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   d_array (input) : The array of data to be sent.                       !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine based on code by Peter Haynes, 25/11/08        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    real(kind=DP), target, intent(in) :: d_array(:,:,:,:)  ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_real_4: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_real_4: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_real_4: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(d_array(1,1,1,1),local_length,real_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_real_4: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
#endif

  end subroutine comms_send_real_4

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_complex_0(node,z_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the (double precision) complex scalar form of        !
    ! comms_send. It sends the complex z_array, or optionally the first       !
    ! length complex pointed to by z_array to node.                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   z_array (input) : The first element of the array of data to be sent.  !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    complex(kind=DP), target, intent(in) :: z_array  ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_complex_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_complex_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(z_array,local_length,complex_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_complex_0: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
#endif

  end subroutine comms_send_complex_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_complex_1(node,z_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the (double precision) complex vector form of        !
    ! comms_send. It sends the first length complex of d_array, or the whole  !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   z_array (input) : The array of data to be sent.                       !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    complex(kind=DP), target, intent(in) :: z_array(:)  ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_complex_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_complex_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_complex_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(z_array(1),local_length,complex_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_complex_1: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle to caller if requested
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
    if (node == 1) return
#endif

  end subroutine comms_send_complex_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_complex_2(node,z_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the (double precision) complex matrix form of        !
    ! comms_send. It sends the first length complex of z_array, or the whole  !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   z_array (input) : The array of data to be sent.                       !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    complex(kind=DP), target, intent(in) :: z_array(:,:)  ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_complex_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_complex_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_complex_2: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(z_array(1,1),local_length,complex_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_complex_2: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
    if (node == 1) return
#endif

  end subroutine comms_send_complex_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_complex_3(node,z_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the (double precision) complex rank 3 tensor form of !
    ! comms_send. It sends the first length complex of z_array, or the whole  !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   z_array (input) : The array of data to be sent.                       !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    complex(kind=DP), target, intent(in) :: z_array(:,:,:) ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_complex_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_complex_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_complex_3: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(z_array(1,1,1),local_length,complex_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_complex_3: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
    if (node == 1) return
#endif

  end subroutine comms_send_complex_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_character_0(node,c_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the character scalar form of comms_send. It sends    !
    ! the character c_array, or optionally the first length characters        !
    ! pointed to by c_array to node.                                          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   c_array (input) : The first element of the array of data to be sent.  !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   character_type    : The MPI type for characters.                      !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(c_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    character, target, intent(in) :: c_array ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_character_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_character_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(c_array,local_length,character_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_character_0: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
    if (node == 1 .and. len(c_array) == 0) return
#endif

  end subroutine comms_send_character_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_logical_0(node,l_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the logical scalar form of comms_send. It sends the  !
    ! logical l_array, or optionally the first length logicals pointed to by  !
    ! l_array to node.                                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   l_array (input) : The first element of the array of data to be sent.  !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   logical_type      : The MPI type for logicals.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(l_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    logical, target, intent(in) :: l_array   ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_logical_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_logical_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(l_array,local_length,logical_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_logical_0: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
    if (node == 1 .and. l_array) return
#endif

  end subroutine comms_send_logical_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_send_logical_1(node,l_array,length,tag,return_handle,comm, &
       add_to_stack)

    !=========================================================================!
    ! This subroutine is the logical vector form of comms_send. It sends the  !
    ! first length logicals of l_array to node, or the whole array if length  !
    ! is absent.                                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)    : The node to which the data will be sent.              !
    !   l_array (input) : The array of data to be sent.                       !
    !   length (input)  : The (optional) number of elements to be sent.       !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   logical_type      : The MPI type for logicals.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(l_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The destination node for the data
    logical, target, intent(in) :: l_array(:)  ! The data to be sent
    integer, optional, intent(in) :: length  ! The number of elements to send
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: return_handle ! Handle for the send
    logical, optional, intent(in) :: add_to_stack

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: handle                        ! Non-blocking request handle
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_send_logical_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(l_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_logical_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(l_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_send_logical_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Check send stack
    call comms_free_send_stack

    ! Call MPI_ISEND to transfer data
    call MPI_ISEND(l_array(1),local_length,logical_type,node,local_tag, &
         local_comm,handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_send_logical_1: MPI_ISEND failed with code ',error
       call comms_abort
    end if

    ! Add handle to send stack if required
    if (present(add_to_stack)) then
       if (add_to_stack) call comms_add_send_stack(handle)
    else
       call comms_add_send_stack(handle)
    end if

    ! Return handle if needed
    if (present(return_handle)) return_handle = handle
#else
    if (present(return_handle)) return_handle = -1
    if (node == 1) return
#endif

  end subroutine comms_send_logical_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_integer_0(node,i_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the integer scalar form of comms_recv. It receives   !
    ! the integer i_array, or optionally the first length integers pointed to !
    ! by i_array from node.                                                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   i_array (output) : The first element of the array of data.            !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)      : The (optional) tag to attach to the message.       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    integer, target, intent(out) :: i_array  ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_integer_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_integer_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(i_array,local_length,integer_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_integer_0: MPI_RECV failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_recv_integer_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_integer_0(node,i_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the integer scalar form of comms_irecv. It receives  !
    ! the integer i_array, or optionally the first length integers pointed to !
    ! by i_array from node.                                                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   i_array (output) : The first element of the array of data.            !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)      : The (optional) tag to attach to the message.       !
    !   handle (output)  : The (optional) handle for the recv to return.      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 12/11/08                                      !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    integer, target, intent(out) :: i_array  ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! Return the handle to this recv

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle argument
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_integer_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_integer_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_IRECV to transfer data
       call MPI_IRECV(i_array,local_length,integer_type,node,local_tag, &
            local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_integer_0: MPI_IRECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_integer_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_integer_1(node,i_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the integer vector form of comms_recv. It receives   !
    ! the first length integers of i_array to node, or the whole array if     !
    ! length is absent.                                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   i_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)      : The (optional) tag to attach to the message.       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    integer, intent(inout) :: i_array(:)     ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_integer_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_integer_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_integer_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(i_array(1),local_length,integer_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_integer_1: MPI_RECV failed with code ',error
       call comms_abort
    end if

#else
    if (node == 1) return
#endif

  end subroutine comms_recv_integer_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_integer_1(node,i_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the integer vector form of comms_irecv. It receives  !
    ! the first length integers of i_array to node, or the whole array if     !
    ! length is absent.                                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   i_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)      : The (optional) tag to attach to the message.       !
    !   handle (output)  : The (optional) handle for the sent to return.      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 12/11/08                                      !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    integer, intent(inout) :: i_array(:)     ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! Return the handle to this recv

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle argument
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument


    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_integer_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_integer_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_integer_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_IRECV to open receive buffer
    call MPI_IRECV(i_array(1),local_length,integer_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_integer_1: MPI_IRECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_integer_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_integer_2(node,i_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the integer matrix form of comms_recv. It receives   !
    ! the first length integers of i_array to node, or the whole array if     !
    ! length is absent.                                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   i_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)      : The (optional) tag to attach to the message.       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    integer, intent(inout) :: i_array(:,:)   ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_integer_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_integer_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_integer_2: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(i_array(1,1),local_length,integer_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_integer_2: MPI_RECV failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_recv_integer_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_integer_2(node,i_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the integer matrix form of comms_irecv. It receives  !
    ! the first length integers of i_array to node, or the whole array if     !
    ! length is absent.                                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   i_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)      : The (optional) tag to attach to the message.       !
    !   handle (output)  : The (optional) handle for the sent to return.      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 12/11/08                                      !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    integer, intent(inout) :: i_array(:,:)   ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! Return the handle to this recv

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle argument
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument


    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_integer_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_integer_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_integer_2: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_IRECV to open receive buffer
    call MPI_IRECV(i_array(1,1),local_length,integer_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_integer_2: MPI_IRECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_integer_2

!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_integer_3(node,i_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the integer rank 3 tensor form of comms_recv. It     !
    ! receives the first length integers of i_array to node, or the whole     !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   i_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    integer, intent(inout) :: i_array(:,:,:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_integer_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_integer_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_integer_3: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(i_array(1,1,1),local_length,integer_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_integer_3: MPI_RECV failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_recv_integer_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_integer_3(node,i_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the integer rank 3 tensor form of comms_irecv. It    !
    ! receives the first length integers of i_array to node, or the whole     !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   i_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 13/12/10                                      !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    integer, intent(inout) :: i_array(:,:,:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle  ! The handle to return

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle argument
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_integer_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_integer_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_integer_3: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_IRECV to transfer data
    call MPI_IRECV(i_array(1,1,1),local_length,integer_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_integer_3: MPI_IRECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_integer_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_real_0(node,d_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real scalar form of           !
    ! comms_recv. It receives the real d_array, or optionally the first       !
    ! length reals pointed to by d_array to node.                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   d_array (output) : The first element of the array of data.            !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), target, intent(out) :: d_array  ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_real_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_real_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(d_array,local_length,real_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_real_0: MPI_RECV failed with code ',error
       call comms_abort
    end if

#else
    if (node == 1) return
#endif

  end subroutine comms_recv_real_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_real_0(node,d_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the (double precision) real scalar form of           !
    ! comms_irecv. It receives the real d_array, or optionally the first      !
    ! length reals pointed to by d_array to node.                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   d_array (output) : The first element of the array of data.            !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 13/12/10                                    !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), target, intent(out) :: d_array  ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! The handle to return

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle argument
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_real_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_real_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_IRECV to transfer data
    call MPI_IRECV(d_array,local_length,real_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_real_0: MPI_IRECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_real_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_real_1(node,d_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real vector form of           !
    ! comms_recv. It receives the first length reals of d_array, or the whole !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   d_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), intent(inout) :: d_array(:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_real_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_real_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_real_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(d_array(1),local_length,real_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_real_1: MPI_RECV failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_recv_real_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_real_1(node,d_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the (double precision) real vector form of           !
    ! comms_recv. It receives the first length reals of d_array, or the whole !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   d_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)      : The (optional) tag to attach to the message.       !
    !   handle (output)  : The (optional) handle for the sent to return.      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 06/11/08                                      !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), intent(inout) :: d_array(:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! The handle to return

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle argument
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_real_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_real_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_real_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_IRECV(d_array(1),local_length,real_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_real_1: MPI_IRECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_real_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_real_2(node,d_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real matrix form of           !
    ! comms_recv. It receives the first length reals of d_array, or the whole !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   d_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), intent(inout) :: d_array(:,:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_real_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_real_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_real_2: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(d_array(1,1),local_length,real_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_real_2: MPI_RECV failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_recv_real_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_real_2(node,d_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the (double precision) real matrix form of           !
    ! comms_recv. It receives the first length reals of d_array, or the whole !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   d_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)      : The (optional) tag to attach to the message.        !
    !   handle (output)  : The (optional) handle for the sent to return.      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 12/11/08                                      !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), intent(inout) :: d_array(:,:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! The handle to return

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle output
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_real_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_real_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_real_2: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_IRECV to transfer data
    call MPI_IRECV(d_array(1,1),local_length,real_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_real_2: MPI_RECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_real_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_real_3(node,d_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real rank 3 tensor form of    !
    ! comms_recv. It receives the first length reals of d_array, or the whole !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   d_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), intent(inout) :: d_array(:,:,:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_real_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_real_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_real_3: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(d_array(1,1,1),local_length,real_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_real_3: MPI_RECV failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_recv_real_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_real_3(node,d_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the (double precision) real rank 3 tensor form of    !
    ! comms_recv. It receives the first length reals of d_array, or the whole !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   d_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)      : The (optional) tag to attach to the message.       !
    !   handle (output)  : The (optional) handle for the sent to return.      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 25/11/08                                      !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), intent(inout) :: d_array(:,:,:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! The handle to return

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle output
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_real_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_real_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_real_3: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_IRECV to transfer data
    call MPI_IRECV(d_array(1,1,1),local_length,real_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_real_3: MPI_IRECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_real_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_real_4(node,d_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real rank 3 tensor form of    !
    ! comms_recv. It receives the first length reals of d_array, or the whole !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   d_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)      : The (optional) tag to attach to the message.       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 25/11/08                                      !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), intent(inout) :: d_array(:,:,:,:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_real_4: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_real_4: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_real_4: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(d_array(1,1,1,1),local_length,real_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_real_4: MPI_RECV failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_recv_real_4

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_real_4(node,d_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the (double precision) real rank 3 tensor form of    !
    ! comms_recv. It receives the first length reals of d_array, or the whole !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   d_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)      : The (optional) tag to attach to the message.       !
    !   handle (output)  : The (optional) handle for the sent to return.      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 25/11/08                                      !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), intent(inout) :: d_array(:,:,:,:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! The handle to return

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle output
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_real_4: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_real_4: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_real_4: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_IRECV(d_array(1,1,1,1),local_length,real_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_real_4: MPI_IRECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_real_4

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_complex_0(node,z_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex scalar form of        !
    ! comms_recv. It receives the complex z_array, or optionally the first    !
    ! length complex pointed to by z_array from node.                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   z_array (output) : The first element of the array of data.            !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    complex(kind=DP), target, intent(out) :: z_array  ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_complex_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_complex_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(z_array,local_length,complex_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_complex_0: MPI_RECV failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_recv_complex_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_complex_0(node,z_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the (double precision) complex scalar form of        !
    ! comms_irecv. It receives the complex z_array, or optionally the first    !
    ! length complex pointed to by z_array from node.                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   z_array (output) : The first element of the array of data.            !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_irecv pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 13/12/10.                                   !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    complex(kind=DP), target, intent(out) :: z_array  ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! The handle to return

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle argument
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_complex_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_complex_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_IRECV to transfer data
    call MPI_IRECV(z_array,local_length,complex_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_complex_0: MPI_IRECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_complex_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_complex_1(node,z_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex vector form of        !
    ! comms_recv. It receives the first length complex of d_array, or the     !
    ! whole array if length is absent.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   z_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    complex(kind=DP), intent(inout) :: z_array(:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_complex_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_complex_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_complex_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(z_array(1),local_length,complex_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_complex_1: MPI_RECV failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_recv_complex_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_complex_1(node,z_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the (double precision) complex vector form of        !
    ! comms_recv. It receives the first length complex of d_array, or the     !
    ! whole array if length is absent.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   z_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)      : The (optional) tag to attach to the message.       !
    !   handle (output)  : The (optional) handle for the sent to return.      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 12/11/08                                      !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    complex(kind=DP), intent(inout) :: z_array(:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! The handle to return

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle output
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_complex_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_complex_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_complex_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_IRECV to transfer data
    call MPI_IRECV(z_array(1),local_length,complex_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_complex_1: MPI_IRECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_complex_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_complex_2(node,z_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex matrix form of        !
    ! comms_recv. It receives the first length complex of z_array, or the     !
    ! whole array if length is absent.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   z_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    complex(kind=DP), intent(inout) :: z_array(:,:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_complex_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_complex_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_complex_2: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(z_array(1,1),local_length,complex_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_complex_2: MPI_RECV failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_recv_complex_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_complex_2(node,z_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the (double precision) complex matrix form of        !
    ! comms_recv. It receives the first length complex of z_array, or the     !
    ! whole array if length is absent.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   z_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 13/12/10.                                   !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    complex(kind=DP), intent(inout) :: z_array(:,:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! The handle to return

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle argument
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_complex_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_complex_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_complex_2: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_IRECV to transfer data
    call MPI_IRECV(z_array(1,1),local_length,complex_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_complex_2: MPI_IRECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_complex_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_complex_3(node,z_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex rank 3 tensor form of !
    ! comms_recv. It receives the first length complex of z_array, or the     !
    ! whole array if length is absent.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   z_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    complex(kind=DP), intent(inout) :: z_array(:,:,:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_complex_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_complex_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_complex_3: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(z_array(1,1,1),local_length,complex_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_complex_3: MPI_RECV failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_recv_complex_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_complex_3(node,z_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the (double precision) complex rank 3 tensor form of !
    ! comms_irecv. It receives the first length complex of z_array, or the    !
    ! whole array if length is absent.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   z_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    complex(kind=DP), intent(inout) :: z_array(:,:,:) ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! The handle to return

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle argument
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_complex_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_complex_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_complex_3: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_IRECV to transfer data
    call MPI_IRECV(z_array(1,1,1),local_length,complex_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_complex_3: MPI_IRECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_complex_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_character_0(node,c_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the character scalar form of comms_recv. It receives !
    ! the character c_array, or optionally the first length characters        !
    ! pointed to by c_array from node.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   c_array (output) : The first element of the array of data.            !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   character_type    : The MPI type for characters.                      !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(c_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none


    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    character, target, intent(out) :: c_array  ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_character_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_character_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(c_array,local_length,character_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_character_0: MPI_RECV failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_recv_character_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_logical_0(node,l_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the logical scalar form of comms_recv. It receives   !
    ! the logical l_array, or optionally the first length logicals pointed to !
    ! by l_array from node.                                                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   l_array (output) : The first element of the array of data.            !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   logical_type      : The MPI type for logicals.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(l_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    logical, target, intent(out) :: l_array  ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_logical_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_logical_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(l_array,local_length,logical_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_logical_0: MPI_RECV failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_recv_logical_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_irecv_logical_0(node,l_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the logical scalar form of comms_irecv. It receives  !
    ! the logical l_array, or optionally the first length logicals pointed to !
    ! by l_array from node.                                                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   l_array (output) : The first element of the array of data.            !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)      : The (optional) tag to attach to the message.       !
    !   handle (output)  : The (optional) handle for the sent to return.      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   logical_type      : The MPI type for logicals.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(l_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine, 06/11/08                                      !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    logical, target, intent(out) :: l_array  ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! The handle to return

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle output
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_logical_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_logical_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_IRECV(l_array,local_length,logical_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_logical_0: MPI_RECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_logical_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_recv_logical_1(node,l_array,length,tag,comm)

    !=========================================================================!
    ! This subroutine is the logical vector form of comms_recv. It receives   !
    ! the first length logicals of l_array from node, or the whole array if   !
    ! length is absent.                                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   l_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   logical_type      : The MPI type for logicals.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(l_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    logical, intent(inout) :: l_array(:)     ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_recv_logical_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(l_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_logical_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(l_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_recv_logical_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_RECV to transfer data
    call MPI_RECV(l_array(1),local_length,logical_type,node,local_tag, &
         local_comm,MPI_STATUS_IGNORE,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_recv_logical_1: MPI_RECV failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_recv_logical_1

  subroutine comms_irecv_logical_1(node,l_array,length,tag,comm,handle)

    !=========================================================================!
    ! This subroutine is the logical vector form of comms_irecv. It receives  !
    ! the first length logicals of l_array from node, or the whole array if   !
    ! length is absent.                                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be received.       !
    !   l_array (output) : The array of data to be received.                  !
    !   length (input)   : The (optional) number of elements to be received.  !
    !   tag (input)     : The (optional) tag to attach to the message.        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   logical_type      : The MPI type for logicals.                        !
    !   default_tag       : Default tag for the MPI_ISEND & MPI_RECV pairs.   !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(l_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    logical, intent(inout) :: l_array(:)     ! The data to be received
    integer, optional, intent(in) :: length  ! The number of elements to receive
    integer, optional, intent(in) :: tag     ! The tag to use
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(out) :: handle ! The handle to return

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: local_handle                  ! Local copy of handle argument
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_tag                     ! Local copy of tag argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_irecv_logical_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(l_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(tag)) then
       if (tag == MPI_ANY_TAG) then
          local_tag = MPI_ANY_TAG
       else
          local_tag = mod(tag,tag_ub+1)
       end if
    else
       local_tag = default_tag
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_logical_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(l_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_irecv_logical_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_IRECV to transfer data
    call MPI_IRECV(l_array(1),local_length,logical_type,node,local_tag, &
         local_comm,local_handle,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_irecv_logical_1: MPI_IRECV failed with code ',error
       call comms_abort
    end if

    if (present(handle)) handle = local_handle
#else
    if (present(handle)) handle = -1
    if (node == 1) return
#endif

  end subroutine comms_irecv_logical_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_integer_0(node,i_array,length,comm)

    !=========================================================================!
    ! This subroutine is the integer scalar form of comms_bcast. It           !
    ! broadcasts the integer i_array, or optionally the first length integers !
    ! pointed to by i_array from node.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be broadcast.      !
    !   i_array (in/out) : The first element of the array of data to be bcast.!
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    integer, target, intent(inout) :: i_array  ! The data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_integer_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_integer_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(i_array,local_length,integer_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_integer_0: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_bcast_integer_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_integer_1(node,i_array,length,comm)

    !=========================================================================!
    ! This subroutine is the integer vector form of comms_bcast. It           !
    ! broadcasts the first length integers of i_array from node, or the whole !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be broadcast.      !
    !   i_array (in/out) : The array of data to be broadcast.                 !
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    integer, intent(inout) :: i_array(:)     ! The data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_integer_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_integer_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_integer_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) then
       return
    end if

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(i_array(1),local_length,integer_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_integer_1: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_bcast_integer_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_integer_2(node,i_array,length,comm)

    !=========================================================================!
    ! This subroutine is the integer matrix form of comms_bcast. It           !
    ! broadcasts the first length integers of i_array from node, or the whole !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be broadcast.      !
    !   i_array (in/out) : The array of data to be broadcast.                 !
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    integer, intent(inout) :: i_array(:,:)   ! The data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_integer_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_integer_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_integer_2: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(i_array(1,1),local_length,integer_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_integer_2: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
    i_array = i_array
#endif

  end subroutine comms_bcast_integer_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_integer_3(node,i_array,length,comm)

    !=========================================================================!
    ! This subroutine is the integer rank 3 tensor form of comms_bcast. It    !
    ! broadcasts the first length integers of i_array from node, or the whole !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node to which the data will be broadcast.        !
    !   i_array (in/out) : The array of data to be broadcast.                 !
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    integer, intent(inout) :: i_array(:,:,:) ! The data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_integer_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_integer_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_integer_3: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(i_array(1,1,1),local_length,integer_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_integer_3: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_bcast_integer_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_real_0(node,d_array,length,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real scalar form of           !
    ! comms_bcast. It broadcasts the real d_array, or optionally the first    !
    ! length real pointed to by d_array from node.                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be broadcast.      !
    !   d_array (in/out) : The first element of the array of data to be bcast.!
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), target, intent(inout) :: d_array  ! The data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_real_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_real_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) then
       return
    end if

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(d_array,local_length,real_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_real_0: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_bcast_real_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_real_1(node,d_array,length,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real vector form of           !
    ! comms_bcast. It broadcasts the first length reals of d_array, or the    !
    ! whole array if length is absent.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be broadcast.      !
    !   d_array (in/out) : The array of data to be broadcast.                 !
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), intent(inout) :: d_array(:) ! The data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_real_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_real_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_real_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) then
       return
    end if

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(d_array(1),local_length,real_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_real_1: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_bcast_real_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_real_2(node,d_array,length,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real matrix form of           !
    ! comms_bcast. It broadcasts the first length reals of d_array, or the    !
    ! whole array if length is absent.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be broadcast.      !
    !   d_array (in/out) : The array of data to be broadcast.                 !
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), intent(inout) :: d_array(:,:) ! The data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_real_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_real_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_real_2: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(d_array(1,1),local_length,real_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_real_2: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
    d_array = d_array
#endif

  end subroutine comms_bcast_real_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_real_3(node,d_array,length,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real rank 3 tensor form of    !
    ! comms_bcast. It broadcasts the first length reals of d_array, or the    !
    ! whole array if length is absent.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be broadcast.      !
    !   d_array (in/out) : The array of data to be broadcast.                 !
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    real(kind=DP), intent(inout) :: d_array(:,:,:) ! The data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_real_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_real_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_real_3: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(d_array(1,1,1),local_length,real_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_real_3: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
    d_array = d_array
#endif

  end subroutine comms_bcast_real_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_complex_0(node,z_array,length,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex scalar form of        !
    ! comms_bcast. It broadcasts the complex z_array, or optionally the first !
    ! length complex pointed to by z_array from node.                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be broadcast.      !
    !   z_array (in/out) : The first element of the array of data to be bcast.!
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    complex(kind=DP), target, intent(inout) :: z_array  ! Data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_complex_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_complex_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) then
       return
    end if

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(z_array,local_length,complex_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_complex_0: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_bcast_complex_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_complex_1(node,z_array,length,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex vector form of        !
    ! comms_bcast. It broadcasts the first length complex of d_array, or the  !
    ! whole array if length is absent.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node to which the data will be broadcast.        !
    !   z_array (in/out) : The array of data to be broadcast.                 !
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    complex(kind=DP), intent(inout) :: z_array(:) ! The data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_complex_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_complex_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_complex_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) then
       return
    end if

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(z_array(1),local_length,complex_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_complex_1: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_bcast_complex_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_complex_2(node,z_array,length,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex matrix form of        !
    ! comms_bcast. It broadcasts the first length complex of z_array, or the  !
    ! whole array if length is absent.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node to which the data will be broadcast.        !
    !   z_array (in/out) : The array of data to be broadcast.                 !
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    complex(kind=DP), intent(inout) :: z_array(:,:) ! The data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_complex_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_complex_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_complex_2: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(z_array(1,1),local_length,complex_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_complex_2: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_bcast_complex_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_complex_3(node,z_array,length,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex rank 3 tensor form of !
    ! comms_bcast. It broadcasts the first length complex of z_array, or the  !
    ! whole array if length is absent.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node to which the data will be broadcast.        !
    !   z_array (in/out) : The array of data to be broadcast.                 !
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for double precision complex.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    complex(kind=DP), intent(inout) :: z_array(:,:,:) ! The data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_complex_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_complex_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_complex_3: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(z_array(1,1,1),local_length,complex_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_complex_3: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_bcast_complex_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_character_0(node,c_array,length,comm)

    !=========================================================================!
    ! This subroutine is the character scalar form of comms_bcast. It         !
    ! broadcasts the character c_array, or optionally the first length        !
    ! characters pointed to by c_array to node.                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node to which the data will be broadcast.        !
    !   c_array (in/out) : The first element of the array of data to be bcast.!
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   character_type    : The MPI type for characters.                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(c_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node                ! The source node for the data
    character(len=*), target, intent(inout) :: c_array  ! Data to be broadcast
    integer, optional, intent(in) :: length    ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_character_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = len(c_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_character_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    ! ifdef USE_CVF is due to the deficiency of MPICH-NT implementation - it
    ! does not handle correctly the length of the string as a hidden argument
#ifdef USE_CVF
    call MPI_BCAST_STR(c_array(1:1),local_length,character_type,node, &
         local_comm,error)
#else
    call MPI_BCAST(c_array(1:1),local_length,character_type,node, &
         local_comm,error)
#endif
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_character_0: MPI_BCAST failed with code ', &
            error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_bcast_character_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_character_1(node,c_array,length,comm)

    !=========================================================================!
    ! This subroutine is the character scalar form of comms_bcast. It         !
    ! broadcasts the character c_array, or optionally the first length        !
    ! characters pointed to by c_array to node.                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node to which the data will be broadcast.        !
    !   c_array (in/out) : The first element of the array of data to be bcast.!
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   character_type    : The MPI type for characters.                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(c_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node                ! The source node for the data
    character(len=*), target, intent(inout) :: c_array(:)  ! Data to be broadcast
    integer, optional, intent(in) :: length    ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_character_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = len(c_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_character_1: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    ! ifdef USE_CVF is due to the deficiency of MPICH-NT implementation - it
    ! does not handle correctly the length of the string as a hidden argument
#ifdef USE_CVF
    call MPI_BCAST_STR(c_array(1),local_length,character_type,node, &
         local_comm,error)
#else
    call MPI_BCAST(c_array(1),local_length,character_type,node, &
         local_comm,error)
#endif
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_character_1: MPI_BCAST failed with code ', &
            error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_bcast_character_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_logical_0(node,l_array,length,comm)

    !=========================================================================!
    ! This subroutine is the logical scalar form of comms_broadcast. It       !
    ! broadcasts the logical l_array, or optionally the first length logicals !
    ! pointed to by l_array from node.                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be broadcast.      !
    !   l_array (in/out) : The first element of the array of data to be bcast.!
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   logical_type      : The MPI type for logicals.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(l_array) >= length (can't check this explicitly)                 !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    logical, target, intent(inout) :: l_array  ! The data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_logical_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_logical_0: length < 0'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(l_array,local_length,logical_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_logical_0: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
    l_array = l_array
#endif

  end subroutine comms_bcast_logical_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_bcast_logical_1(node,l_array,length,comm)

    !=========================================================================!
    ! This subroutine is the logical vector form of comms_bcast. It           !
    ! broadcasts the first length logicals of l_array from node, or the whole !
    ! array if length is absent.                                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! node (input)     : The node from which the data will be broadcast.      !
    !   l_array (in/out) : The array of data to be broadcast.                 !
    !   length (input)   : The (optional) number of elements to be broadcast. !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   logical_type      : The MPI type for logicals.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(l_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 16/7/02                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: node              ! The source node for the data
    logical, intent(inout) :: l_array(:)     ! The data to be broadcast
    integer, optional, intent(in) :: length  ! The number of elements to bcast
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
#ifdef MPI
    integer :: error                         ! Error flag
#endif

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_bcast_logical_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(l_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_logical_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(l_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_bcast_logical_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Call MPI_BCAST to transfer data
    call MPI_BCAST(l_array(1),local_length,logical_type,node, &
         local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_bcast_logical_1: MPI_BCAST failed with code ',error
       call comms_abort
    end if
#else
    if (node == 1) return
#endif

  end subroutine comms_bcast_logical_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_integer_0(op,i_array,comm,root)

    !=========================================================================!
    ! This subroutine is the integer scalar form of comms_reduce. It performs !
    ! a global reduction of i_array across all nodes and the result is        !
    ! returned in i_array.                                                    !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'SUM'  - global sum                                                   !
    !   'PROD' - global product                                               !
    !   'MAX'  - global maximum                                               !
    !   'MIN'  - global minimum                                               !
    !   'AND'  - bitwise AND                                                  !
    !   'OR'   - bitwise OR                                                   !
    !   'XOR'  - bitwise XOR                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   i_array (in/out) : On entry, the data to be reduced. On exit, the     !
    !                      reduced data.                                      !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op       ! The reduction operation
    integer, intent(inout) :: i_array        ! On entry, data to be reduced.
                                             ! On exit the reduced data
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce

#ifdef MPI
    ! Local variables
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    integer :: iwork                         ! Local temporary variable
#endif
    integer :: local_comm                    ! Local copy of comm argument
    logical :: local_recv                    ! Is this node receiving data?


    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_integer_0: comms not initialised'
       call comms_abort
    end if

    ! Check arguments
    if (op /= 'SUM' .and. op /= 'PROD' .and. op /= 'MAX' .and. &
         op /= 'MIN' .and. op /= 'AND' .and. op /= 'OR' .and. &
         op /= 'XOR') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_integer_0: invalid operation "',op,'"'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(root)) then
       local_recv = (root==pub_my_node_id)
    else
       local_recv = .true.
    end if

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('SUM')  ; mpi_op = MPI_SUM
       case ('PROD') ; mpi_op = MPI_PROD
       case ('MAX')  ; mpi_op = MPI_MAX
       case ('MIN')  ; mpi_op = MPI_MIN
       case ('AND')  ; mpi_op = MPI_BAND
       case ('OR')   ; mpi_op = MPI_BOR
       case ('XOR')  ; mpi_op = MPI_BXOR
    end select

    ! Make a local copy of the data
    iwork = i_array

    ! Perform the reduction
    if (.not.present(root)) then
       call MPI_ALLREDUCE(iwork,i_array,1,integer_type,mpi_op, &
            local_comm,error)
    else
       call MPI_REDUCE(iwork,i_array,1,integer_type,mpi_op,root, &
            local_comm,error)
    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in &
            &comms_reduce_integer_0: MPI_(ALL)REDUCE failed with code ',error
       call comms_abort
    end if
#else
    i_array = i_array
#endif


  end subroutine comms_reduce_integer_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_integer_long_0(op,i_array,comm,root)

    !=========================================================================!
    ! This subroutine is the integer scalar form of comms_reduce. It performs !
    ! a global reduction of i_array across all nodes and the result is        !
    ! returned in i_array.                                                    !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'SUM'  - global sum                                                   !
    !   'PROD' - global product                                               !
    !   'MAX'  - global maximum                                               !
    !   'MIN'  - global minimum                                               !
    !   'AND'  - bitwise AND                                                  !
    !   'OR'   - bitwise OR                                                   !
    !   'XOR'  - bitwise XOR                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   i_array (in/out) : On entry, the data to be reduced. On exit, the     !
    !                      reduced data.                                      !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !-------------------------------------------------------------------------!
    ! Cloned by Jacek Dziedzic 06/06/2011 from comms_reduce_integer_0,        !
    !                                     written by Peter Haynes, 15/7/03    !
    !=========================================================================!

    use constants, only: LONG

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op           ! The reduction operation
    integer(kind=LONG), intent(inout) :: i_array ! On entry, data to be reduced.
                                             ! On exit the reduced data
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce

#ifdef MPI
    ! Local variables
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    integer(kind=LONG) :: iwork              ! Local temporary variable
#endif
    integer :: local_comm                    ! Local copy of comm argument
    logical :: local_recv                    ! Is this node receiving data?


    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_integer_0: comms not initialised'
       call comms_abort
    end if

    ! Check arguments
    if (op /= 'SUM' .and. op /= 'PROD' .and. op /= 'MAX' .and. &
         op /= 'MIN' .and. op /= 'AND' .and. op /= 'OR' .and. &
         op /= 'XOR') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_integer_0: invalid operation "',op,'"'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(root)) then
       local_recv = (root==pub_my_node_id)
    else
       local_recv = .true.
    end if

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('SUM')  ; mpi_op = MPI_SUM
       case ('PROD') ; mpi_op = MPI_PROD
       case ('MAX')  ; mpi_op = MPI_MAX
       case ('MIN')  ; mpi_op = MPI_MIN
       case ('AND')  ; mpi_op = MPI_BAND
       case ('OR')   ; mpi_op = MPI_BOR
       case ('XOR')  ; mpi_op = MPI_BXOR
    end select

    ! Make a local copy of the data
    iwork = i_array

    ! Perform the reduction
    if (.not.present(root)) then
       call MPI_ALLREDUCE(iwork,i_array,1,integer_kind_long_type,mpi_op, &
            local_comm,error)
    else
       call MPI_REDUCE(iwork,i_array,1,integer_kind_long_type,mpi_op,root, &
            local_comm,error)
    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in &
            &comms_reduce_integer_long_0: MPI_(ALL)REDUCE failed with code ',error
       call comms_abort
    end if
#else
    i_array = i_array
#endif


  end subroutine comms_reduce_integer_long_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_integer_1(op,i_array,length,comm,root)

    !=========================================================================!
    ! This subroutine is the integer vector form of comms_reduce. It performs !
    ! a global reduction of i_array (or optionally the first length integers  !
    ! in i_array across all nodes and the result is returned in i_array.      !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'SUM'  - global sum                                                   !
    !   'PROD' - global product                                               !
    !   'MAX'  - global maximum                                               !
    !   'MIN'  - global minimum                                               !
    !   'AND'  - bitwise AND                                                  !
    !   'OR'   - bitwise OR                                                   !
    !   'XOR'  - bitwise XOR                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   i_array (in/out) : On entry, the array of data to be reduced.         !
    !                      On exit, the reduced data.                         !
    !   length (input)   : The (optional) number of elements to reduce.       !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op       ! The reduction operation
    integer, intent(inout) :: i_array(:)     ! On entry, data to be reduced.
                                             ! On exit the reduced data
    integer, optional, intent(in) :: length  ! The number of elements to reduce
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce

    ! Local variables
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    logical :: local_recv                    ! Is this node receiving data?
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    integer :: iwork                         ! Local temporary variable
#endif


    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_integer_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(root)) then
       local_recv = (root==pub_my_node_id)
    else
       local_recv = .true.
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_integer_1: length < 0'
       call comms_abort
    end if
    if ((local_length > size(i_array)).and.local_recv) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_integer_1: length exceeds array size'
       call comms_abort
    end if
    if (present(root)) then
       if ((root<0).or.(root>pub_total_num_nodes)) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_integer_1: invalid root specified'
       end if
    end if
    if (op /= 'SUM' .and. op /= 'PROD' .and. op /= 'MAX' .and. &
         op /= 'MIN' .and. op /= 'AND' .and. op /= 'OR' .and. &
         op /= 'XOR') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_integer_1: invalid operation "',op,'"'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) then
       return
    end if

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('SUM')  ; mpi_op = MPI_SUM
       case ('PROD') ; mpi_op = MPI_PROD
       case ('MAX')  ; mpi_op = MPI_MAX
       case ('MIN')  ; mpi_op = MPI_MIN
       case ('AND')  ; mpi_op = MPI_BAND
       case ('OR')   ; mpi_op = MPI_BOR
       case ('XOR')  ; mpi_op = MPI_BXOR
    end select

    ! MPI needs separate send and receive buffers, so single element is a
    ! special case:
    if (local_length == 1) then

       ! Make a local copy of the data
       iwork = i_array(1)

       ! Perform the reduction
       if (.not.present(root)) then
          call MPI_ALLREDUCE(iwork,i_array(1),1,integer_type,mpi_op, &
               local_comm,error)
       else
          call MPI_REDUCE(iwork,i_array(1),1,integer_type,mpi_op,root, &
               local_comm,error)
       end if

    else

       ! Allocate workspace if required
       if (allocated(iwork1)) then
          if (size(iwork1) < local_length) then
             deallocate(iwork1,stat=error)
             if (error /= 0) then
                if (pub_on_root) write(stdout,'(a,i3)') &
                     'Error in comms_reduce_integer_1: &
                     &deallocating iwork1 failed with code ',error
                call comms_abort
             end if
          end if
       end if
       if (.not. allocated(iwork1)) then
          allocate(iwork1(local_length),stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce_&
                  &integer_1: allocating iwork1 failed with code ',error
             call comms_abort
          end if
       end if

       ! Copy data into local array
       iwork1(1:local_length) = i_array(1:local_length)

       ! Perform reduction
       if (.not.present(root)) then
          call MPI_ALLREDUCE(iwork1(1),i_array(1),local_length,integer_type, &
               mpi_op,local_comm,error)
       else
          call MPI_REDUCE(iwork1(1),i_array(1),local_length,integer_type, &
               mpi_op,root,local_comm,error)
       end if

    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce_&
            &integer_1: MPI_(ALL)REDUCE failed with code ',error
       call comms_abort
    end if

#endif


  end subroutine comms_reduce_integer_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_integer_2(op,i_array,length,comm,root)

    !=========================================================================!
    ! This subroutine is the integer matrix form of comms_reduce. It performs !
    ! a global reduction of i_array (or optionally the first length integers  !
    ! in i_array across all nodes and the result is returned in i_array.      !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'SUM'  - global sum                                                   !
    !   'PROD' - global product                                               !
    !   'MAX'  - global maximum                                               !
    !   'MIN'  - global minimum                                               !
    !   'AND'  - bitwise AND                                                  !
    !   'OR'   - bitwise OR                                                   !
    !   'XOR'  - bitwise XOR                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   i_array (in/out) : On entry, the array of data to be reduced.         !
    !                      On exit, the reduced data.                         !
    !   length (input)   : The (optional) number of elements to reduce.       !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op       ! The reduction operation
    integer, intent(inout) :: i_array(:,:)   ! On entry, data to be reduced.
                                             ! On exit the reduced data
    integer, optional, intent(in) :: length  ! The number of elements to reduce
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce

    ! Local variables
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    integer :: iwork                         ! Local temporary variable
#endif

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_integer_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_integer_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_integer_2: length exceeds array size'
       call comms_abort
    end if
    if (present(root)) then
       if ((root<0).or.(root>pub_total_num_nodes)) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_integer_2: invalid root specified'
       end if
    end if
    if (op /= 'SUM' .and. op /= 'PROD' .and. op /= 'MAX' .and. &
         op /= 'MIN' .and. op /= 'AND' .and. op /= 'OR' .and. &
         op /= 'XOR') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_integer_2: invalid operation "',op,'"'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('SUM')  ; mpi_op = MPI_SUM
       case ('PROD') ; mpi_op = MPI_PROD
       case ('MAX')  ; mpi_op = MPI_MAX
       case ('MIN')  ; mpi_op = MPI_MIN
       case ('AND')  ; mpi_op = MPI_BAND
       case ('OR')   ; mpi_op = MPI_BOR
       case ('XOR')  ; mpi_op = MPI_BXOR
    end select

    ! MPI needs separate send and receive buffers, so single element is a
    ! special case:
    if (local_length == 1) then

       ! Make a local copy of the data
       iwork = i_array(1,1)

       ! Perform the reduction
       if (.not.present(root)) then
          call MPI_ALLREDUCE(iwork,i_array(1,1),1,integer_type,mpi_op, &
               local_comm,error)
       else
          call MPI_REDUCE(iwork,i_array(1,1),1,integer_type,mpi_op,root, &
               local_comm,error)
       end if

    else

       ! Allocate workspace if required
       if (allocated(iwork1)) then
          if (size(iwork1) < local_length) then
             deallocate(iwork1,stat=error)
             if (error /= 0) then
                if (pub_on_root) write(stdout,'(a,i3)') &
                     'Error in comms_reduce_integer_2: &
                     &deallocating iwork1 failed with code ',error
                call comms_abort
             end if
          end if
       end if
       if (.not. allocated(iwork1)) then
          allocate(iwork1(local_length),stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') &
                  'Error in comms_reduce_integer_2: &
                  &allocating iwork1 failed with code ',error
             call comms_abort
          end if
       end if

       ! Make a local copy of the data
       call comms_copy(local_length,i_array,iwork1)

       ! Perform reduction
       if (.not.present(root)) then
          call MPI_ALLREDUCE(iwork1(1),i_array(1,1),local_length,integer_type, &
               mpi_op,local_comm,error)
       else
          call MPI_REDUCE(iwork1(1),i_array(1,1),local_length,integer_type, &
               mpi_op,root,local_comm,error)
       end if

    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce&
            &_integer_2: MPI_(ALL)REDUCE failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_reduce_integer_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_integer_3(op,i_array,length,comm,root)

    !=========================================================================!
    ! This subroutine is the integer rank 3 tensor form of comms_reduce. It   !
    ! performs a global reduction of i_array (or optionally the first length  !
    ! integer in i_array across all nodes and the result is returned in       !
    ! i_array.                                                                !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'SUM'  - global sum                                                   !
    !   'PROD' - global product                                               !
    !   'MAX'  - global maximum                                               !
    !   'MIN'  - global minimum                                               !
    !   'AND'  - bitwise AND                                                  !
    !   'OR'   - bitwise OR                                                   !
    !   'XOR'  - bitwise XOR                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   i_array (in/out) : On entry, the array of data to be reduced.         !
    !                      On exit, the reduced data.                         !
    !   length (input)   : The (optional) number of elements to reduce.       !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(i_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op       ! The reduction operation
    integer, intent(inout) :: i_array(:,:,:) ! On entry, data to be reduced.
                                             ! On exit the reduced data
    integer, optional, intent(in) :: length  ! The number of elements to reduce
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce

    ! Local variables
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    integer :: iwork                         ! Local temporary variable
#endif

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_integer_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_integer_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_integer_3: length exceeds array size'
       call comms_abort
    end if
    if (present(root)) then
       if ((root<0).or.(root>pub_total_num_nodes)) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_integer_3: invalid root specified'
       end if
    end if
    if (op /= 'SUM' .and. op /= 'PROD' .and. op /= 'MAX' .and. &
         op /= 'MIN' .and. op /= 'AND' .and. op /= 'OR' .and. &
         op /= 'XOR') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_integer_3: invalid operation "',op,'"'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('SUM')  ; mpi_op = MPI_SUM
       case ('PROD') ; mpi_op = MPI_PROD
       case ('MAX')  ; mpi_op = MPI_MAX
       case ('MIN')  ; mpi_op = MPI_MIN
       case ('AND')  ; mpi_op = MPI_BAND
       case ('OR')   ; mpi_op = MPI_BOR
       case ('XOR')  ; mpi_op = MPI_BXOR
    end select

    ! MPI needs separate send and receive buffers, so single element is a
    ! special case:
    if (local_length == 1) then

       ! Make a local copy of the data
       iwork = i_array(1,1,1)

       ! Perform the reduction
       if (.not.present(root)) then
          call MPI_ALLREDUCE(iwork,i_array(1,1,1),1,integer_type,mpi_op, &
               local_comm,error)
       else
          call MPI_REDUCE(iwork,i_array(1,1,1),1,integer_type,mpi_op,root, &
               local_comm,error)
       end if

    else

       ! Allocate workspace if required
       if (allocated(iwork1)) then
          if (size(iwork1) < local_length) then
             deallocate(iwork1,stat=error)
             if (error /= 0) then
                if (pub_on_root) write(stdout,'(a,i3)') &
                     'Error in comms_reduce_integer_3: &
                     &deallocating iwork1 failed with code ',error
                call comms_abort
             end if
          end if
       end if
       if (.not. allocated(iwork1)) then
          allocate(iwork1(local_length),stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') &
                  'Error in comms_reduce_integer_3: &
                  &allocating iwork1 failed with code ',error
             call comms_abort
          end if
       end if

       ! Make a local copy of the data
       call comms_copy(local_length,i_array,iwork1)

       ! Perform reduction
       if (.not.present(root)) then
          call MPI_ALLREDUCE(iwork1(1),i_array(1,1,1),local_length,integer_type, &
               mpi_op,local_comm,error)
       else
          call MPI_REDUCE(iwork1(1),i_array(1,1,1),local_length,integer_type, &
               mpi_op,root,local_comm,error)
       end if

    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce&
            &_integer_3: MPI_(ALL)REDUCE failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_reduce_integer_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_real_0(op,d_array,comm,root)

    !=========================================================================!
    ! This subroutine is the (double precision) real scalar form of           !
    ! comms_reduce. It performs a global reduction of d_array across all      !
    ! nodes and the result is returned in d_array.                            !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'SUM'  - global sum                                                   !
    !   'PROD' - global product                                               !
    !   'MAX'  - global maximum                                               !
    !   'MIN'  - global minimum                                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   d_array (in/out) : On entry, the data to be reduced. On exit, the     !
    !                      reduced data.                                      !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for (double precision) reals.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op       ! The reduction operation
    real(kind=DP), intent(inout) :: d_array  ! On entry, data to be reduced.
                                             ! On exit the reduced data
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    real(kind=DP) :: dwork                   ! Local temporary variable
#endif
    integer :: local_comm                    ! Local copy of comm argument


    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_real_0: comms not initialised'
       call comms_abort
    end if

    ! Check arguments
    if (op /= 'SUM' .and. op /= 'PROD' .and. op /= 'MAX' .and. &
         op /= 'MIN') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_real_0: invalid operation "',op,'"'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('SUM')  ; mpi_op = MPI_SUM
       case ('PROD') ; mpi_op = MPI_PROD
       case ('MAX')  ; mpi_op = MPI_MAX
       case ('MIN')  ; mpi_op = MPI_MIN
    end select

    ! Make a local copy of the data
    dwork = d_array

    ! Perform the reduction
    if (.not.present(root)) then
       call MPI_ALLREDUCE(dwork,d_array,1,real_type,mpi_op, &
            local_comm,error)
    else
       call MPI_REDUCE(dwork,d_array,1,real_type,mpi_op,root, &
            local_comm,error)
    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce_real_0: &
            &MPI_(ALL)REDUCE failed with code ',error
       call comms_abort
    end if
#else
    d_array = d_array
#endif


  end subroutine comms_reduce_real_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_real_1(op,d_array,length,comm,root,d_array_src)

    !=========================================================================!
    ! This subroutine is the (double precision) real vector form of           !
    ! comms_reduce. It performs a global reduction of d_array (or optionally  !
    ! the first length reals in d_array across all nodes and the result is    !
    ! returned in d_array.                                                    !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'SUM'  - global sum                                                   !
    !   'PROD' - global product                                               !
    !   'MAX'  - global maximum                                               !
    !   'MIN'  - global minimum                                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   d_array (in/out) : On entry, the array of data to be reduced.         !
    !                      On exit, the reduced data.                         !
    !   length (input)   : The (optional) number of elements to reduce.       !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for (double precision) reals.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op       ! The reduction operation
    real(kind=DP), intent(inout) :: d_array(:) ! On entry, data to be reduced.
                                               ! On exit the reduced data
    integer, optional, intent(in) :: length  ! The number of elements to reduce
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce
    real(kind=DP), optional, intent(in) :: d_array_src(:) ! Optional: data to be reduced.

    ! Local variables
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    logical :: local_recv                    ! Is this node receiving data?
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    real(kind=DP) :: dwork                   ! Local temporary variable
#endif


    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_real_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(root)) then
       local_recv = (root==pub_my_node_id)
    else
       local_recv = .true.
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_real_1: length < 0'
       call comms_abort
    end if
    if ((local_length > size(d_array)).and.local_recv) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_real_1: length exceeds array size'
       call comms_abort
    end if
    if (op /= 'SUM' .and. op /= 'PROD' .and. op /= 'MAX' .and. &
         op /= 'MIN') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_real_1: invalid operation "',op,'"'
       call comms_abort
    end if
    if (present(root)) then
       if ((root<0).or.(root>pub_total_num_nodes)) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_real_1: invalid root specified'
       end if
    end if
    if (present(d_array_src)) then
       if (size(d_array_src)<local_length) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_real_1: length exceeds size of &
               &source array'
       end if
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) then
       return
    end if

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('SUM')  ; mpi_op = MPI_SUM
       case ('PROD') ; mpi_op = MPI_PROD
       case ('MAX')  ; mpi_op = MPI_MAX
       case ('MIN')  ; mpi_op = MPI_MIN
    end select

    ! MPI needs separate send and receive buffers, so single element is a
    ! special case:
    if (local_length == 1) then

       ! Make a local copy of the data
       dwork = d_array(1)

       if (.not.present(d_array_src)) then

          ! Make a local copy of the data
          dwork = d_array(1)

          ! Perform the reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(dwork,d_array(1),1,real_type,mpi_op, &
                  local_comm,error)
          else
             call MPI_REDUCE(dwork,d_array(1),1,real_type,mpi_op,root, &
                  local_comm,error)
          end if

       else

          ! Perform the reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(d_array_src(1),d_array(1),1,real_type,mpi_op, &
                  local_comm,error)
          else
             call MPI_REDUCE(d_array_src(1),d_array(1),1,real_type,mpi_op,root, &
                  local_comm,error)
          end if

       end if

    else

       if (.not.present(d_array_src)) then

          ! Allocate workspace if required
          if (allocated(dwork1)) then
             if (size(dwork1) < local_length) then
                deallocate(dwork1,stat=error)
                if (error /= 0) then
                   if (pub_on_root) write(stdout,'(a,i3)') &
                        'Error in comms_reduce_real_1: &
                        &deallocating dwork1 failed with code ',error
                   call comms_abort
                end if
             end if
          end if
          if (.not. allocated(dwork1)) then
             allocate(dwork1(local_length),stat=error)
             if (error /= 0) then
                if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce&
                     &_real_1: allocating dwork1 failed with code ',error
                call comms_abort
             end if
          end if

          ! Copy data into local array
          dwork1(1:local_length) = d_array(1:local_length)

          ! Perform reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(dwork1(1),d_array(1),local_length,real_type, &
                  mpi_op,local_comm,error)
          else
             call MPI_REDUCE(dwork1(1),d_array(1),local_length,real_type, &
                  mpi_op,root,local_comm,error)
          end if

       else

          ! Perform reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(d_array_src(1),d_array(1),local_length,real_type, &
                  mpi_op,local_comm,error)
          else
             call MPI_REDUCE(d_array_src(1),d_array(1),local_length,real_type, &
                  mpi_op,root,local_comm,error)
          end if

       end if

    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce_real_1:&
            &MPI_(ALL)REDUCE failed with code ',error
       call comms_abort
    end if
#else
    if (present(d_array_src)) then
       d_array(:) = d_array_src(:)
    end if
#endif


  end subroutine comms_reduce_real_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_real_2(op,d_array,length,comm,root,d_array_src)

    !=========================================================================!
    ! This subroutine is the (double precision) real matrix form of           !
    ! comms_reduce. It performs a global reduction of d_array (or optionally  !
    ! the first length reals in d_array across all nodes and the result is    !
    ! returned in d_array.                                                    !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'SUM'  - global sum                                                   !
    !   'PROD' - global product                                               !
    !   'MAX'  - global maximum                                               !
    !   'MIN'  - global minimum                                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   d_array (in/out) : On entry, the array of data to be reduced.         !
    !                      On exit, the reduced data.                         !
    !   length (input)   : The (optional) number of elements to reduce.       !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for (double precision) reals.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op       ! The reduction operation
    real(kind=DP), intent(inout) :: d_array(:,:) ! On entry, data to be reduced.
                                                 ! On exit the reduced data
    integer, optional, intent(in) :: length  ! The number of elements to reduce
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce
    real(kind=DP), optional, intent(in) :: d_array_src(:,:) ! On entry data to be reduced

    ! Local variables
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    logical :: local_recv                    ! Is this node receiving data?
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    real(kind=DP) :: dwork                   ! Local temporary variable
#endif

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_real_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(root)) then
       local_recv = (root==pub_my_node_id)
    else
       local_recv = .true.
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_real_2: length < 0'
       call comms_abort
    end if
    if ((local_length > size(d_array)).and.local_recv) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_real_2: length exceeds array size'
       call comms_abort
    end if
    if (op /= 'SUM' .and. op /= 'PROD' .and. op /= 'MAX' .and. &
         op /= 'MIN') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_real_2: invalid operation "',op,'"'
       call comms_abort
    end if
    if (present(root)) then
       if ((root<0).or.(root>pub_total_num_nodes)) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_real_2: invalid root specified'
       end if
    end if
    if (present(d_array_src)) then
       if (size(d_array_src)<local_length) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_real_2: length exceeds size of &
               &source array'
       end if
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('SUM')  ; mpi_op = MPI_SUM
       case ('PROD') ; mpi_op = MPI_PROD
       case ('MAX')  ; mpi_op = MPI_MAX
       case ('MIN')  ; mpi_op = MPI_MIN
    end select

    ! MPI needs separate send and receive buffers, so single element is a
    ! special case:
    if (local_length == 1) then

       if (.not.present(d_array_src)) then

          ! Make a local copy of the data
          dwork = d_array(1,1)

          ! Perform the reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(dwork,d_array(1,1),1,real_type,mpi_op, &
                  local_comm,error)
          else
             call MPI_REDUCE(dwork,d_array(1,1),1,real_type,mpi_op,root, &
                  local_comm,error)
          end if

       else

          ! Perform the reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(d_array_src(1,1),d_array(1,1),1,real_type,&
                  mpi_op,local_comm,error)
          else
             call MPI_REDUCE(d_array_src(1,1),d_array(1,1),1,real_type,mpi_op, &
                  root,local_comm,error)
          end if

       end if

    else

       if (.not.present(d_array_src)) then

          ! Allocate workspace if required
          if (allocated(dwork1)) then
             if (size(dwork1) < local_length) then
                deallocate(dwork1,stat=error)
                if (error /= 0) then
                   if (pub_on_root) write(stdout,'(a,i3)') &
                        'Error in comms_reduce_real_2: &
                        &deallocating dwork1 failed with code ',error
                   call comms_abort
                end if
             end if
          end if
          if (.not. allocated(dwork1)) then
             allocate(dwork1(local_length),stat=error)
             if (error /= 0) then
                if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce&
                     &_real_2: allocating dwork1 failed with code ',error
                call comms_abort
             end if
          end if

          ! Make a local copy of the data
          call comms_copy(local_length,d_array,dwork1)

          ! Perform reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(dwork1(1),d_array(1,1),local_length,real_type, &
                  mpi_op,local_comm,error)
          else
             call MPI_REDUCE(dwork1(1),d_array(1,1),local_length,real_type, &
                  mpi_op,root,local_comm,error)
          end if

       else

          ! Perform reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(d_array_src(1,1),d_array(1,1),local_length, &
                  real_type,mpi_op,local_comm,error)
          else
             call MPI_REDUCE(d_array_src(1,1),d_array(1,1),local_length, &
                  real_type,mpi_op,root,local_comm,error)
          end if

       end if

    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce_real_2: &
            &MPI_(ALL)REDUCE failed with code ',error
       call comms_abort
    end if
#else
    if (present(d_array_src)) then
       d_array(:,:) = d_array_src(:,:)
    end if
#endif

  end subroutine comms_reduce_real_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_real_3(op,d_array,length,comm,root,d_array_src)

    !=========================================================================!
    ! This subroutine is the (double precision) real rank 3 tensor form of    !
    ! comms_reduce. It performs a global reduction of d_array (or optionally  !
    ! the first length reals in d_array across all nodes and the result is    !
    ! returned in d_array.                                                    !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'SUM'  - global sum                                                   !
    !   'PROD' - global product                                               !
    !   'MAX'  - global maximum                                               !
    !   'MIN'  - global minimum                                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   d_array (in/out) : On entry, the array of data to be reduced.         !
    !                      On exit, the reduced data.                         !
    !   length (input)   : The (optional) number of elements to reduce.       !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for (double precision) reals.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(d_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op       ! The reduction operation
    real(kind=DP), intent(inout) :: d_array(:,:,:) ! On entry data to be reduced
                                             ! On exit the reduced data
    integer, optional, intent(in) :: length  ! The number of elements to reduce
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce
    real(kind=DP), optional, intent(in) :: d_array_src(:,:,:) ! On entry data to be reduced

    ! Local variables
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    logical :: local_recv                    ! Is this node receiving data?
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    real(kind=DP) :: dwork                   ! Local temporary variable
#endif

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_real_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(d_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(root)) then
       local_recv = (root==pub_my_node_id)
    else
       local_recv = .true.
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_real_3: length < 0'
       call comms_abort
    end if
    if ((local_length > size(d_array)).and.local_recv) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_real_3: length exceeds array size'
       call comms_abort
    end if
    if (op /= 'SUM' .and. op /= 'PROD' .and. op /= 'MAX' .and. &
         op /= 'MIN') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_real_3: invalid operation "',op,'"'
       call comms_abort
    end if
    if (present(root)) then
       if ((root<0).or.(root>pub_total_num_nodes)) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_real_3: invalid root specified'
       end if
    end if
    if (present(d_array_src)) then
       if (size(d_array_src)<local_length) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_real_3: length exceeds size of &
               &source array'
       end if
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('SUM')  ; mpi_op = MPI_SUM
       case ('PROD') ; mpi_op = MPI_PROD
       case ('MAX')  ; mpi_op = MPI_MAX
       case ('MIN')  ; mpi_op = MPI_MIN
    end select

    ! MPI needs separate send and receive buffers, so single element is a
    ! special case:
    if (local_length == 1) then

       if (.not.present(d_array_src)) then

          ! Make a local copy of the data
          dwork = d_array(1,1,1)

          ! Perform the reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(dwork,d_array(1,1,1),1,real_type,mpi_op, &
                  local_comm,error)
          else
             call MPI_REDUCE(dwork,d_array(1,1,1),1,real_type,mpi_op,root, &
                  local_comm,error)
          end if

       else

          ! Perform the reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(d_array_src(1,1,1),d_array(1,1,1),1,real_type, &
                  mpi_op,local_comm,error)
          else
             call MPI_REDUCE(d_array_src(1,1,1),d_array(1,1,1),1,real_type, &
                  mpi_op,root,local_comm,error)
          end if

       end if

    else

       if (.not.present(d_array_src)) then

          ! Allocate workspace if required
          if (allocated(dwork1)) then
             if (size(dwork1) < local_length) then
                deallocate(dwork1,stat=error)
                if (error /= 0) then
                   if (pub_on_root) write(stdout,'(a,i3)') &
                        'Error in comms_reduce_real_3: &
                        &deallocating dwork1 failed with code ',error
                   call comms_abort
                end if
             end if
          end if
          if (.not. allocated(dwork1)) then
             allocate(dwork1(local_length),stat=error)
             if (error /= 0) then
                if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce&
                     &_real_3: allocating dwork1 failed with code ',error
                call comms_abort
             end if
          end if

          ! Make a local copy of the data
          call comms_copy(local_length,d_array,dwork1)

          ! Perform reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(dwork1(1),d_array(1,1,1),local_length, &
                  real_type,mpi_op,local_comm,error)
          else
             call MPI_REDUCE(dwork1(1),d_array(1,1,1),local_length, &
                  real_type,mpi_op,root,local_comm,error)
          end if

       else

          ! Perform reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(d_array_src(1,1,1),d_array(1,1,1),local_length, &
                  real_type,mpi_op,local_comm,error)
          else
             call MPI_REDUCE(d_array_src(1,1,1),d_array(1,1,1),local_length, &
                  real_type,mpi_op,root,local_comm,error)
          end if

       end if

    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce&
            &_real_3: MPI_(ALL)REDUCE failed with code ',error
       call comms_abort
    end if
#else
    if (present(d_array_src)) then
       d_array(:,:,:) = d_array_src(:,:,:)
    end if
#endif

  end subroutine comms_reduce_real_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_complex_0(op,z_array,comm,root)

    !=========================================================================!
    ! This subroutine is the (double precision) complex scalar form of        !
    ! comms_reduce. It performs a global reduction of z_array across all      !
    ! nodes and the result is returned in z_array.                            !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'SUM'  - global sum                                                   !
    !   'PROD' - global product                                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   z_array (in/out) : On entry, The first element of the data to be      !
    !                      reduced. On exit, the first element of the reduced !
    !                      data.                                              !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for (double precision) complex.      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op         ! The reduction operation
    complex(kind=DP), intent(inout) :: z_array ! On entry, data to be reduced
                                               ! On exit the reduced data
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce

    ! Local variables
    integer :: local_comm                    ! Local copy of comm argument
    logical :: local_recv                    ! Is this node receiving data?
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    complex(kind=DP) :: zwork                ! Local temporary variable
#endif


    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_complex_0: comms not initialised'
       call comms_abort
    end if

    ! Check arguments
    if (op /= 'SUM' .and. op /= 'PROD') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_complex_0: invalid operation "',op,'"'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(root)) then
       local_recv = (root==pub_my_node_id)
    else
       local_recv = .true.
    end if

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('SUM')  ; mpi_op = MPI_SUM
       case ('PROD') ; mpi_op = MPI_PROD
    end select

    ! Make a local copy of the data
    zwork = z_array

    ! Perform the reduction
    if (.not.present(root)) then
       call MPI_ALLREDUCE(zwork,z_array,1,complex_type,mpi_op, &
            local_comm,error)
    else
       call MPI_REDUCE(zwork,z_array,1,complex_type,mpi_op,root, &
            local_comm,error)
    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce_complex&
            &_0: MPI_(ALL)REDUCE failed with code ',error
       call comms_abort
    end if
#else
    z_array = z_array
#endif


  end subroutine comms_reduce_complex_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_complex_1(op,z_array,length,comm,root,z_array_src)

    !=========================================================================!
    ! This subroutine is the (double precision) complex vector form of        !
    ! comms_reduce. It performs a global reduction of z_array (or optionally  !
    ! the first length complex in z_array across all nodes and the result is  !
    ! returned in z_array.                                                    !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'SUM'  - global sum                                                   !
    !   'PROD' - global product                                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   z_array (in/out) : On entry, the array of data to be reduced.         !
    !                      On exit, the reduced data.                         !
    !   length (input)   : The (optional) number of elements to reduce.       !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for (double precision) complex.      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op       ! The reduction operation
    complex(kind=DP), intent(inout) :: z_array(:) ! On entry data to be reduced.
                                                  ! On exit the reduced data
    integer, optional, intent(in) :: length  ! The number of elements to reduce
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce
    complex(kind=DP), optional, intent(in) :: z_array_src(:) ! Optional: data to be reduced.


    ! Local variables
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    logical :: local_recv                    ! Is this node receiving data?
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    complex(kind=DP) :: zwork                ! Local temporary variable
#endif


    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_complex_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(root)) then
       local_recv = (root==pub_my_node_id)
    else
       local_recv = .true.
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_complex_1: length < 0'
       call comms_abort
    end if
    if ((local_length > size(z_array)).and.local_recv) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_complex_1: length exceeds array size'
       call comms_abort
    end if
    if (op /= 'SUM' .and. op /= 'PROD') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_complex_1: invalid operation "',op,'"'
       call comms_abort
    end if
    if (present(root)) then
       if ((root<0).or.(root>pub_total_num_nodes)) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_complex_1: invalid root specified'
       end if
    end if
    if (present(z_array_src)) then
       if (size(z_array_src)<local_length) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_complex_1: length exceeds size of &
               &source array'
       end if
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) then
       return
    end if

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('SUM')  ; mpi_op = MPI_SUM
       case ('PROD') ; mpi_op = MPI_PROD
    end select

    ! MPI needs separate send and receive buffers, so single element is a
    ! special case:
    if (local_length == 1) then

       if (.not.present(z_array_src)) then

          ! Make a local copy of the data
          zwork = z_array(1)

          ! Perform the reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(zwork,z_array(1),1,complex_type,mpi_op, &
                  local_comm,error)
          else
              call MPI_REDUCE(zwork,z_array(1),1,complex_type,mpi_op,root, &
                  local_comm,error)
          end if

       else

          ! Perform the reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(z_array_src,z_array(1),1,complex_type,mpi_op, &
                  local_comm,error)
          else
             call MPI_REDUCE(z_array_src,z_array(1),1,complex_type,mpi_op,root,&
                  local_comm,error)
          end if

       end if

    else

       if (.not.present(z_array_src)) then

          ! Allocate workspace if required
          if (allocated(zwork1)) then
             if (size(zwork1) < local_length) then
                deallocate(zwork1,stat=error)
                if (error /= 0) then
                   if (pub_on_root) write(stdout,'(a,i3)') &
                        'Error in comms_reduce_complex_1: &
                         &deallocating zwork1 failed with code ',error
                   call comms_abort
                end if
             end if
          end if
          if (.not. allocated(zwork1)) then
             allocate(zwork1(local_length),stat=error)
             if (error /= 0) then
                if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce&
                     &_complex_1: allocating zwork1 failed with code ',error
                call comms_abort
             end if
          end if

          ! Copy data into local array
          zwork1(1:local_length) = z_array(1:local_length)

          ! Perform reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(zwork1(1),z_array(1),local_length,complex_type, &
                  mpi_op,local_comm,error)
          else
             call MPI_REDUCE(zwork1(1),z_array(1),local_length,complex_type, &
                  mpi_op,root,local_comm,error)
          end if

       else

          ! Perform reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(z_array_src(1),z_array(1),local_length,complex_type, &
                  mpi_op,local_comm,error)
          else
             call MPI_REDUCE(z_array_src(1),z_array(1),local_length,complex_type, &
                  mpi_op,root,local_comm,error)
          end if

       end if

    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce_complex_1:&
            &MPI_REDUCE failed with code ',error
       call comms_abort
    end if

#else
    if (present(z_array_src)) then
       z_array(:) = z_array_src(:)
    end if
#endif


  end subroutine comms_reduce_complex_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_complex_2(op,z_array,length,comm,root,z_array_src)

    !=========================================================================!
    ! This subroutine is the (double precision) complex matrix form of        !
    ! comms_reduce. It performs a global reduction of z_array (or optionally  !
    ! the first length complexs in z_array across all nodes and the result is !
    ! returned in z_array.                                                    !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'SUM'  - global sum                                                   !
    !   'PROD' - global product                                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   z_array (in/out) : On entry, the array of data to be reduced.         !
    !                      On exit, the reduced data.                         !
    !   length (input)   : The (optional) number of elements to reduce.       !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for (double precision) complex.      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op       ! The reduction operation
    complex(kind=DP), intent(inout) :: z_array(:,:) ! On entry data for reducing
                                                    ! On exit the reduced data
    integer, optional, intent(in) :: length  ! The number of elements to reduce
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce
    complex(kind=DP), optional, intent(in) :: z_array_src(:,:) ! Optional: data to be reduced.

    ! Local variables
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    logical :: local_recv                    ! Is this node receiving data?
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    complex(kind=DP) :: zwork                ! Local temporary variable
#endif

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_complex_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(root)) then
       local_recv = (root==pub_my_node_id)
    else
       local_recv = .true.
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_complex_2: length < 0'
       call comms_abort
    end if
    if ((local_length > size(z_array)).and.local_recv) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_complex_2: length exceeds array size'
       call comms_abort
    end if
    if (op /= 'SUM' .and. op /= 'PROD') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_complex_2: invalid operation "',op,'"'
       call comms_abort
    end if
    if (present(root)) then
       if ((root<0).or.(root>pub_total_num_nodes)) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_complex_2: invalid root specified'
       end if
    end if
    if (present(z_array_src)) then
       if (size(z_array_src)<local_length) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_complex_2: length exceeds size of &
               &source array'
       end if
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('SUM')  ; mpi_op = MPI_SUM
       case ('PROD') ; mpi_op = MPI_PROD
    end select

    ! MPI needs separate send and receive buffers, so single element is a
    ! special case:
    if (local_length == 1) then

       if (.not.present(z_array_src)) then

          ! Make a local copy of the data
          zwork = z_array(1,1)

          ! Perform the reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(zwork,z_array(1,1),1,complex_type,mpi_op, &
                  local_comm,error)
          else
             call MPI_REDUCE(zwork,z_array(1,1),1,complex_type,mpi_op,root, &
                  local_comm,error)
          end if

       else

          ! Perform the reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(z_array_src(1,1),z_array(1,1),1,complex_type, &
                  mpi_op,local_comm,error)
          else
             call MPI_REDUCE(z_array_src(1,1),z_array(1,1),1,complex_type, &
                  mpi_op,root,local_comm,error)
          end if

       end if

    else

       if (.not.present(z_array_src)) then

          ! Allocate workspace if required
          if (allocated(zwork1)) then
             if (size(zwork1) < local_length) then
                deallocate(zwork1,stat=error)
                if (error /= 0) then
                   if (pub_on_root) write(stdout,'(a,i3)') &
                        'Error in comms_reduce_complex_2: &
                        &deallocating zwork1 failed with code ',error
                   call comms_abort
                end if
             end if
          end if
          if (.not. allocated(zwork1)) then
             allocate(zwork1(local_length),stat=error)
             if (error /= 0) then
                if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce&
                     &_complex_2: allocating zwork1 failed with code ',error
                call comms_abort
             end if
          end if

          ! Make a local copy of the data
          call comms_copy(local_length,z_array,zwork1)

          ! Perform reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(zwork1(1),z_array(1,1),local_length,complex_type, &
                  mpi_op,local_comm,error)
          else
             call MPI_REDUCE(zwork1(1),z_array(1,1),local_length,complex_type, &
                  mpi_op,root,local_comm,error)
          end if

       else

          ! Perform reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(z_array_src(1,1),z_array(1,1),local_length, &
                  complex_type,mpi_op,local_comm,error)
          else
             call MPI_REDUCE(z_array_src(1,1),z_array(1,1),local_length, &
                  complex_type,mpi_op,root,local_comm,error)
          end if

       end if

    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce&
            &_complex_2: MPI_REDUCE failed with code ',error
       call comms_abort
    end if

#else
    if (present(z_array_src)) then
       z_array(:,:) = z_array_src(:,:)
    end if

#endif

  end subroutine comms_reduce_complex_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_complex_3(op,z_array,length,comm,root,z_array_src)

    !=========================================================================!
    ! This subroutine is the (double precision) complex rank 3 tensor form of !
    ! comms_reduce. It performs a global reduction of z_array (or optionally  !
    ! the first length complexs in z_array across all nodes and the result is !
    ! returned in z_array.                                                    !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'SUM'  - global sum                                                   !
    !   'PROD' - global product                                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   z_array (in/out) : On entry, the array of data to be reduced.         !
    !                      On exit, the reduced data.                         !
    !   length (input)   : The (optional) number of elements to reduce.       !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for (double precision) complex.      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(z_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op       ! The reduction operation
    complex(kind=DP), intent(inout) :: z_array(:,:,:) ! On entry data for reduce
                                                      ! On exit the reduced data
    integer, optional, intent(in) :: length  ! The number of elements to reduce
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce
    complex(kind=DP), optional, intent(in) :: z_array_src(:,:,:) ! Optional: data to be reduced.

    ! Local variables
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    logical :: local_recv                    ! Is this node receiving data?
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    complex(kind=DP) :: zwork                ! Local temporary variable
#endif

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_complex_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(z_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(root)) then
       local_recv = (root==pub_my_node_id)
    else
       local_recv = .true.
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_complex_3: length < 0'
       call comms_abort
    end if
    if ((local_length > size(z_array)).and.local_recv) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_complex_3: length exceeds array size'
       call comms_abort
    end if
    if (op /= 'SUM' .and. op /= 'PROD') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_complex_3: invalid operation "',op,'"'
       call comms_abort
    end if
    if (present(root)) then
       if ((root<0).or.(root>pub_total_num_nodes)) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_complex_3: invalid root specified'
       end if
    end if
    if (present(z_array_src)) then
       if (size(z_array_src)<local_length) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_complex_3: length exceeds size of &
               &source array'
       end if
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('SUM')  ; mpi_op = MPI_SUM
       case ('PROD') ; mpi_op = MPI_PROD
    end select

    ! MPI needs separate send and receive buffers, so single element is a
    ! special case:
    if (local_length == 1) then

       if (.not.present(z_array_src)) then

          ! Make a local copy of the data
          zwork = z_array(1,1,1)

          ! Perform the reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(zwork,z_array(1,1,1),1,complex_type,mpi_op, &
                  local_comm,error)
          else
             call MPI_REDUCE(zwork,z_array(1,1,1),1,complex_type,mpi_op,root, &
                  local_comm,error)
          end if

       else

          ! Perform the reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(z_array_src(1,1,1),z_array(1,1,1),1,complex_type,mpi_op, &
                  local_comm,error)
          else
             call MPI_REDUCE(z_array_src(1,1,1),z_array(1,1,1),1,complex_type,mpi_op,root, &
                  local_comm,error)
          end if

       end if

    else

       if (.not.present(z_array_src)) then

          ! Allocate workspace if required
          if (allocated(zwork1)) then
             if (size(zwork1) < local_length) then
                deallocate(zwork1,stat=error)
                if (error /= 0) then
                   if (pub_on_root) write(stdout,'(a,i3)') &
                        'Error in comms_reduce_complex_3: &
                        &deallocating zwork1 failed with code ',error
                   call comms_abort
                end if
             end if
          end if
          if (.not. allocated(zwork1)) then
             allocate(zwork1(local_length),stat=error)
             if (error /= 0) then
                if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce&
                     &_complex_3: allocating zwork1 failed with code ',error
                call comms_abort
             end if
          end if

          ! Make a local copy of the data
          call comms_copy(local_length,z_array,zwork1)

          ! Perform reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(zwork1(1),z_array(1,1,1),local_length, &
                  complex_type,mpi_op,local_comm,error)
          else
             call MPI_REDUCE(zwork1(1),z_array(1,1,1),local_length, &
                  complex_type,mpi_op,root,local_comm,error)
          end if

       else

          ! Perform reduction
          if (.not.present(root)) then
             call MPI_ALLREDUCE(z_array_src(1,1,1),z_array(1,1,1), &
                  local_length,complex_type,mpi_op,local_comm,error)
          else
             call MPI_REDUCE(z_array_src(1,1,1),z_array(1,1,1), &
                  local_length,complex_type,mpi_op,root,local_comm,error)
          end if

       end if
    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce&
            &_complex_3: MPI_(ALL)REDUCE failed with code ',error
       call comms_abort
    end if

#else
    if (present(z_array_src)) then
       z_array(:,:,:) = z_array_src(:,:,:)
    end if

#endif

  end subroutine comms_reduce_complex_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_logical_0(op,l_array,comm,root)

    !=========================================================================!
    ! This subroutine is the logical scalar form of comms_reduce. It performs !
    ! a global reduction of l_array across all nodes and the result is        !
    ! returned in l_array.                                                    !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'AND'  - logical AND                                                  !
    !   'OR'   - logical OR                                                   !
    !   'XOR'  - logical XOR                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   l_array (in/out) : On entry, The first element of the data to be      !
    !                      reduced. On exit, the first element of the reduced !
    !                      data.                                              !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   logical_type      : The MPI type for logicals.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op       ! The reduction operation
    logical, intent(inout) :: l_array        ! On entry, data to be reduced.
                                             ! On exit the reduced data
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce

#ifdef MPI
    ! Local variables
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    logical :: lwork                         ! Local temporary variable
#endif
    integer :: local_comm                    ! Local copy of comm argument
    logical :: local_recv                    ! Is this node receiving data?

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_logical_0: comms not initialised'
       call comms_abort
    end if

    ! Check arguments
    if (op /= 'AND' .and. op /= 'OR' .and. op /= 'XOR') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_logical_0: invalid operation "',op,'"'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(root)) then
       local_recv = (root==pub_my_node_id)
    else
       local_recv = .true.
    end if

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('AND')  ; mpi_op = MPI_LAND
       case ('OR')   ; mpi_op = MPI_LOR
       case ('XOR')  ; mpi_op = MPI_LXOR
    end select

    ! Make a local copy of the data
    lwork = l_array

    ! Perform the reduction
    if (.not.present(root)) then
       call MPI_ALLREDUCE(lwork,l_array,1,logical_type,mpi_op, &
            local_comm,error)
    else
       call MPI_REDUCE(lwork,l_array,1,logical_type,mpi_op,root, &
            local_comm,error)
    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce_logical_0: &
            &MPI_(ALL)REDUCE failed with code ',error
       call comms_abort
    end if
#else
    l_array = l_array
#endif

  end subroutine comms_reduce_logical_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_reduce_logical_1(op,l_array,length,comm,root)

    !=========================================================================!
    ! This subroutine is the logical vector form of comms_reduce. It performs !
    ! a global reduction of l_array (or optionally the first length logicals  !
    ! in l_array across all nodes and the result is returned in l_array.      !
    ! The reduction operation which will be performed is determined by the    !
    ! value of the argument op:                                               !
    !   'AND'  - logical AND                                                  !
    !   'OR'   - logical OR                                                   !
    !   'XOR'  - logical XOR                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   op (input)       : String as above denoting the reduction operation.  !
    !   l_array (in/out) : On entry, the array of data to be reduced.         !
    !                      On exit, the reduced data.                         !
    !   length (input)   : The (optional) number of elements to reduce.       !
    !   root (input)     : The (optional) root process on which the result is !
    !                      required. If not present, all processes will get   !
    !                      the result (ie allreduce will be used).            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   logical_type      : The MPI type for logicals.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0 (if present)                                              !
    !   size(l_array) >= length                                               !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 15/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    character(len=*), intent(in) :: op       ! The reduction operation
    logical, intent(inout) :: l_array(:)     ! On entry, data to be reduced.
                                             ! On exit the reduced data
    integer, optional, intent(in) :: length  ! The number of elements to reduce
    integer, optional, intent(in) :: comm    ! The communicator to use
    integer, optional, intent(in) :: root    ! The root process for the reduce

    ! Local variables
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    logical :: local_recv                    ! Is this node receiving data?
#ifdef MPI
    integer :: error                         ! Error flag
    integer :: mpi_op                        ! The MPI operation parameter
    logical :: lwork                         ! Local temporary variable
#endif

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_reduce_logical_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(l_array)
    end if
    if (present(comm)) then
       local_comm = comm
    else
       local_comm = communicator
    end if
    if (present(root)) then
       local_recv = (root==pub_my_node_id)
    else
       local_recv = .true.
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_logical_1: length < 0'
       call comms_abort
    end if
    if ((local_length > size(l_array)).and.local_recv) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_reduce_logical_1: length exceeds array size'
       call comms_abort
    end if
    if (op /= 'AND' .and. op /= 'OR' .and. op /= 'XOR') then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in comms_reduce_logical_1: invalid operation "',op,'"'
       call comms_abort
    end if
    if (present(root)) then
       if ((root<0).or.(root>pub_total_num_nodes)) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in comms_reduce_logical_1: invalid root specified'
       end if
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! Determine the operation to be performed
    select case (op)
       case ('AND')  ; mpi_op = MPI_LAND
       case ('OR')   ; mpi_op = MPI_LOR
       case ('XOR')  ; mpi_op = MPI_LXOR
    end select

    ! MPI needs separate send and receive buffers, so single element is a
    ! special case:
    if (local_length == 1) then

       ! Make a local copy of the data
       lwork = l_array(1)

       ! Perform the reduction
       if (.not.present(root)) then
          call MPI_ALLREDUCE(lwork,l_array(1),1,logical_type,mpi_op, &
               local_comm,error)
       else
          call MPI_REDUCE(lwork,l_array(1),1,logical_type,mpi_op,root, &
               local_comm,error)
       end if

    else

       ! Allocate workspace if required
       if (allocated(lwork1)) then
          if (size(lwork1) < local_length) then
             deallocate(lwork1,stat=error)
             if (error /= 0) then
                if (pub_on_root) write(stdout,'(a,i3)') &
                     'Error in comms_reduce_logical_1: &
                     &deallocating lwork1 failed with code ',error
                call comms_abort
             end if
          end if
       end if
       if (.not. allocated(lwork1)) then
          allocate(lwork1(local_length),stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') &
                  'Error in comms_reduce_logical_1: &
                  &allocating lwork1 failed with code ',error
             call comms_abort
          end if
       end if

       ! Copy data into local array
       lwork1(1:local_length) = l_array(1:local_length)

       ! Perform reduction
       if (.not.present(root)) then
          call MPI_ALLREDUCE(lwork1(1),l_array(1),local_length,logical_type, &
               mpi_op,local_comm,error)
       else
          call MPI_REDUCE(lwork1(1),l_array(1),local_length,logical_type, &
               mpi_op,root,local_comm,error)
       end if

    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_reduce&
            &_logical_1: MPI_(ALL)REDUCE failed with code ',error
       call comms_abort
    end if

#endif

  end subroutine comms_reduce_logical_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_alltoall_integer_1(i_array,length,comm)

    !=========================================================================!
    ! This subroutine is the integer vector form of comms_alltoall. Each node !
    ! sends length (default 1) integers to every other node according to the  !
    ! order in which they appear in i_array.                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   i_array (in/out) : On entry, the data to be sent by this node.        !
    !                      On exit, the data received.                        !
    !   length (input)   : The (optional) number of elements to send to each  !
    !                      node.                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0                                                           !
    !   size(i_array) >= length*pub_total_num_nodes                           !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 22/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(inout) :: i_array(:)    ! On entry, data to be sent.
                                            ! On exit the data received
    integer, optional, intent(in) :: length ! The number of elements to send to/
                                            ! receive from each node
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                        ! Error flag
#endif
    integer :: local_length                 ! Local copy of length
    integer :: local_comm                   ! Local copy of comm argument
    integer :: comm_size                    ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) &
            'Error in comms_alltoall_integer_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_integer_1: length < 0'
       call comms_abort
    end if
    if (size(i_array) < local_length*comm_size) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_integer_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! MPI needs separate send and receive buffers, so allocate workspace
    ! if required
    if (allocated(iwork1)) then
       if (size(iwork1) < local_length*comm_size) then
          deallocate(iwork1,stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') &
                  'Error in comms_alltoall_integer_1: &
                  &deallocating iwork1 failed with code ',error
             call comms_abort
          end if
       end if
    end if
    if (.not. allocated(iwork1)) then
       allocate(iwork1(local_length*comm_size),stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
               &_integer_1: allocating iwork1 failed with code ',error
          call comms_abort
       end if
    end if

    ! Copy data into local array
    iwork1(1:local_length*comm_size) = &
         i_array(1:local_length*comm_size)

    ! Perform alltoall
    call MPI_ALLTOALL(iwork1(1),local_length,integer_type,i_array(1), &
         local_length,integer_type,local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
            &_integer_1: MPI_ALLTOALL failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_alltoall_integer_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_alltoall_integer_2(i_array,length,comm)

    !=========================================================================!
    ! This subroutine is the integer matrix form of comms_alltoall. Each node !
    ! sends size(i_array,1) integers to every other node according to the     !
    ! order in which they appear in i_array.                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   i_array (in/out) : On entry, the data to be sent by this node.        !
    !                      On exit, the data received.                        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   size(i_array,2) >= pub_total_num_nodes                                !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 22/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(inout) :: i_array(:,:)  ! On entry, data to be sent.
                                            ! On exit the data received
    integer, optional, intent(in) :: comm   ! The communicator to use
    integer, optional, intent(in) :: length ! Amount of data to send/receive

    ! Local variables
#ifdef MPI
    integer :: error                        ! Error flag
#endif
    integer :: local_length                 ! Local length of comm argument
    integer :: local_comm                   ! Local copy of comm argument
    integer :: comm_size                    ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) &
            'Error in comms_alltoall_integer_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = size(i_array,1)
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (size(i_array) < comm_size*local_length) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_integer_2: array size too small.'
       call comms_abort
    end if

#ifdef MPI
    ! MPI needs separate send and receive buffers, so allocate workspace
    ! if required
    if (allocated(iwork1)) then
       if (size(iwork1) < local_length*comm_size) then
          deallocate(iwork1,stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') &
                  'Error in comms_alltoall_integer_2: &
                  &deallocating iwork1 failed with code ',error
             call comms_abort
          end if
       end if
    end if
    if (.not. allocated(iwork1)) then
       allocate(iwork1(local_length*comm_size),stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') &
               'Error in comms_alltoall_integer_2: &
               &allocating iwork1 failed with code ',error
          call comms_abort
       end if
    end if

    ! Make a local copy of the data
    call comms_copy(local_length*comm_size,i_array,iwork1)

    ! Perform alltoall
    call MPI_ALLTOALL(iwork1(1),local_length,integer_type,i_array(1,1),local_length, &
         integer_type,local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
            &_integer_2: MPI_ALLTOALL failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_alltoall_integer_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_alltoall_integer_3(i_array,comm)

    !=========================================================================!
    ! This subroutine is the integer rank 3 tensor form of comms_alltoall.    !
    ! Each node sends size(i_array,1)*size(i_array,2) integers to every other !
    ! node according to the order in which they appear in i_array.            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   i_array (in/out) : On entry, the data to be sent by this node.        !
    !                      On exit, the data received.                        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type      : The MPI type for integers.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   size(i_array,3) >= pub_total_num_nodes (or size of communicator)      !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 22/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(inout) :: i_array(:,:,:) ! On entry, data to be sent
                                             ! On exit the data received
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                        ! Error flag
#endif
    integer :: length                       ! Amount of data to send/receive
    integer :: local_comm                   ! Local copy of comm argument
    integer :: comm_size                    ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) &
            'Error in comms_alltoall_integer_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (size(i_array,3) < comm_size) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_integer_3: array size too small.'
       call comms_abort
    end if

    ! Get length - number of elements to send/receive to/from each node
    length = size(i_array,1)*size(i_array,2)

#ifdef MPI
    ! MPI needs separate send and receive buffers, so allocate workspace
    ! if required
    if (allocated(iwork1)) then
       if (size(iwork1) < length*comm_size) then
          deallocate(iwork1,stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') &
                  'Error in comms_alltoall_integer_3: &
                  &deallocating iwork1 failed with code ',error
             call comms_abort
          end if
       end if
    end if
    if (.not. allocated(iwork1)) then
       allocate(iwork1(length*comm_size),stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
               &_integer_3: allocating iwork1 failed with code ',error
          call comms_abort
       end if
    end if

    ! Make a local copy of the data using an external subroutine to avoid
    ! the problem that the ranks of iwork1 and i_array differ
    call comms_copy(length*comm_size,i_array,iwork1)

    ! Perform alltoall
    call MPI_ALLTOALL(iwork1(1),length,integer_type,i_array(1,1,1),length, &
         integer_type,local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
            &_integer_3: MPI_ALLTOALL failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_alltoall_integer_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_alltoall_real_1(d_array,length,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real vector form of           !
    ! comms_alltoall. Each node sends length (default 1) reals to every other !
    ! node according to the order in which they appear in d_array.            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   d_array (in/out) : On entry, the data to be sent by this node.        !
    !                      On exit, the data received.                        !
    !   length (input)   : The (optional) number of elements to send to each  !
    !                      node.                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for (double precision) reals.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0                                                           !
    !   size(d_array) >= length*pub_total_num_nodes                           !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 22/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: d_array(:) ! On entry, data to be sent.
                                            ! On exit the data received
    integer, optional, intent(in) :: length ! The number of elements to send to/
                                            ! receive from each node
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                        ! Error flag
#endif
    integer :: local_length                 ! Local copy of length
    integer :: local_comm                   ! Local copy of comm argument
    integer :: comm_size                    ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_alltoall_real_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_real_1: length < 0'
       call comms_abort
    end if
    if (size(d_array) < local_length*comm_size) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_real_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! MPI needs separate send and receive buffers, so allocate workspace
    ! if required
    if (allocated(dwork1)) then
       if (size(dwork1) < local_length*comm_size) then
          deallocate(dwork1,stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') &
                  'Error in comms_alltoall_real_1: &
                  &deallocating dwork1 failed with code ',error
             call comms_abort
          end if
       end if
    end if
    if (.not. allocated(dwork1)) then
       allocate(dwork1(local_length*comm_size),stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
               &_real_1: allocating dwork1 failed with code ',error
          call comms_abort
       end if
    end if

    ! Copy data into local array
    dwork1(1:local_length*comm_size) = &
         d_array(1:local_length*comm_size)

    ! Perform alltoall
    call MPI_ALLTOALL(dwork1(1),local_length,real_type,d_array(1), &
         local_length,real_type,local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
            &_real_1: MPI_ALLTOALL failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_alltoall_real_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_alltoall_real_2(d_array,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real matrix form of           !
    ! comms_alltoall. Each node sends size(d_array,1) reals to every other    !
    ! node according to the order in which they appear in d_array.            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   d_array (in/out) : On entry, the data to be sent by this node.        !
    !                      On exit, the data received.                        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for (double precision) reals.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   size(d_array,2) >= pub_total_num_nodes                                !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 22/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: d_array(:,:) ! On entry, data to be sent.
                                            ! On exit the data received
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                        ! Error flag
#endif
    integer :: length                       ! Amount of data to send/receive
    integer :: local_comm                   ! Local copy of comm argument
    integer :: comm_size                    ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_alltoall_real_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = comm_size
    end if

    ! Check arguments
    if (size(d_array,2) < comm_size) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_real_2: array size too small.'
       call comms_abort
    end if

    ! Get length - number of elements to send/receive to/from each node
    length = size(d_array,1)

#ifdef MPI
    ! MPI needs separate send and receive buffers, so allocate workspace
    ! if required
    if (allocated(dwork1)) then
       if (size(dwork1) < length*comm_size) then
          deallocate(dwork1,stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') &
                  'Error in comms_alltoall_real_2: &
                  &deallocating dwork1 failed with code ',error
             call comms_abort
          end if
       end if
    end if
    if (.not. allocated(dwork1)) then
       allocate(dwork1(length*comm_size),stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
               &_real_2: allocating dwork1 failed with code ',error
          call comms_abort
       end if
    end if

    ! Make a local copy of the data
    call comms_copy(length*comm_size,d_array,dwork1)

    ! Perform alltoall
    call MPI_ALLTOALL(dwork1(1),length,real_type,d_array(1,1),length, &
         real_type,local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
            &_real_2: MPI_ALLTOALL failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_alltoall_real_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_alltoall_real_3(d_array,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real rank 3 tensor form of    !
    ! comms_alltoall. Each node sends size(d_array,1)*size(d_array,2) reals   !
    ! to every other node according to the order in which they appear in      !
    ! d_array.                                                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   d_array (in/out) : On entry, the data to be sent by this node.        !
    !                      On exit, the data received.                        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for (double precision) reals.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   size(d_array,2) >= pub_total_num_nodes                                !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 22/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: d_array(:,:,:) ! On entry, data to be sent.
                                            ! On exit the data received
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                        ! Error flag
#endif
    integer :: length                       ! Amount of data to send/receive
    integer :: local_comm                   ! Local copy of comm argument
    integer :: comm_size                    ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_alltoall_real_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = comm_size
    end if

    ! Check arguments
    if (size(d_array,3) < comm_size) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_real_3: array size too small.'
       call comms_abort
    end if

    ! Get length - number of elements to send/receive to/from each node
    length = size(d_array,1)*size(d_array,2)

#ifdef MPI
    ! MPI needs separate send and receive buffers, so allocate workspace
    ! if required
    if (allocated(dwork1)) then
       if (size(dwork1) < length*comm_size) then
          deallocate(dwork1,stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') &
                  'Error in comms_alltoall_real_3: &
                  &deallocating dwork1 failed with code ',error
             call comms_abort
          end if
       end if
    end if
    if (.not. allocated(dwork1)) then
       allocate(dwork1(length*comm_size),stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
               &_real_3: allocating dwork1 failed with code ',error
          call comms_abort
       end if
    end if

    ! Make a local copy of the data using an external subroutine to avoid
    ! the problem that the ranks of iwork1 and i_array differ
    call comms_copy(length*comm_size,d_array,dwork1)

    ! Perform alltoall
    call MPI_ALLTOALL(dwork1(1),length,real_type,d_array(1,1,1),length, &
         real_type,local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
            &_real_3: MPI_ALLTOALL failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_alltoall_real_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_alltoall_complex_1(z_array,length,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex vector form of        !
    ! comms_alltoall. Each node sends length (default 1) complex to every     !
    ! other node according to the order in which they appear in z_array.      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   z_array (in/out) : On entry, the data to be sent by this node.        !
    !                      On exit, the data received.                        !
    !   length (input)   : The (optional) number of elements to send to each  !
    !                      node.                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for (double precision) complex.      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0                                                           !
    !   size(z_array) >= length*pub_total_num_nodes                           !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 22/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout) :: z_array(:) ! On entry, data to be sent.
                                            ! On exit the data received
    integer, optional, intent(in) :: length ! The number of elements to send to/
                                            ! receive from each node
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                        ! Error flag
#endif
    integer :: local_length                 ! Local copy of length
    integer :: local_comm                   ! Local copy of comm argument
    integer :: comm_size                    ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) &
            'Error in comms_alltoall_complex_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_complex_1: length < 0'
       call comms_abort
    end if
    if (size(z_array) < local_length*comm_size) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_complex_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! MPI needs separate send and receive buffers, so allocate workspace
    ! if required
    if (allocated(zwork1)) then
       if (size(zwork1) < local_length*comm_size) then
          deallocate(zwork1,stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') &
                  'Error in comms_alltoall_complex_1: &
                  &deallocating zwork1 failed with code ',error
             call comms_abort
          end if
       end if
    end if
    if (.not. allocated(zwork1)) then
       allocate(zwork1(local_length*comm_size),stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
               &_complex_1: allocating zwork1 failed with code ',error
          call comms_abort
       end if
    end if

    ! Copy data into local array
    zwork1(1:local_length*comm_size) = &
         z_array(1:local_length*comm_size)

    ! Perform alltoall
    call MPI_ALLTOALL(zwork1(1),local_length,complex_type,z_array(1), &
         local_length,complex_type,local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
            &_complex_1: MPI_ALLTOALL failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_alltoall_complex_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_alltoall_complex_2(z_array,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex matrix form of        !
    ! comms_alltoall. Each node sends size(z_array,1) complex to every other  !
    ! node according to the order in which they appear in z_array.            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   z_array (in/out) : On entry, the data to be sent by this node.        !
    !                      On exit, the data received.                        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for (double precision) complex.      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   size(z_array,2) >= pub_total_num_nodes                                !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 22/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout) :: z_array(:,:) ! On entry, data to be sent.
                                            ! On exit the data received
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                        ! Error flag
#endif
    integer :: length                       ! Amount of data to send/receive
    integer :: local_comm                   ! Local copy of comm argument
    integer :: comm_size                    ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) &
            'Error in comms_alltoall_complex_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (size(z_array,2) < comm_size) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_complex_2: array size too small.'
       call comms_abort
    end if

    ! Get length - number of elements to send/receive to/from each node
    length = size(z_array,1)

#ifdef MPI
    ! MPI needs separate send and receive buffers, so allocate workspace
    ! if required
    if (allocated(zwork1)) then
       if (size(zwork1) < length*comm_size) then
          deallocate(zwork1,stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') &
                  'Error in comms_alltoall_complex_2: &
                  &deallocating zwork1 failed with code ',error
             call comms_abort
          end if
       end if
    end if
    if (.not. allocated(zwork1)) then
       allocate(zwork1(length*comm_size),stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
               &_complex_2: allocating zwork1 failed with code ',error
          call comms_abort
       end if
    end if

    ! Make a local copy of the data
    call comms_copy(length*comm_size,z_array,zwork1)

    ! Perform alltoall
    call MPI_ALLTOALL(zwork1(1),length,complex_type,z_array(1,1),length, &
         complex_type,local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
            &_complex_2: MPI_ALLTOALL failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_alltoall_complex_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_alltoall_complex_3(z_array,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex rank 3 tensor form of !
    ! comms_alltoall. Each node sends size(z_array,1)*size(z_array,2) complex !
    ! to every other node according to the order in which they appear in      !
    ! z_array.                                                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   z_array (in/out) : On entry, the data to be sent by this node.        !
    !                      On exit, the data received.                        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type      : The MPI type for (double precision) complex.      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   size(z_array,2) >= pub_total_num_nodes                                !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 22/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout) :: z_array(:,:,:) ! On entry, data to send.
                                            ! On exit the data received
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                        ! Error flag
#endif
    integer :: length                       ! Amount of data to send/receive
    integer :: local_comm                   ! Local copy of comm argument
    integer :: comm_size                    ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) &
            'Error in comms_alltoall_complex_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (size(z_array,3) < comm_size) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_complex_3: array size too small.'
       call comms_abort
    end if

    ! Get length - number of elements to send/receive to/from each node
    length = size(z_array,1)*size(z_array,2)

#ifdef MPI
    ! MPI needs separate send and receive buffers, so allocate workspace
    ! if required
    if (allocated(zwork1)) then
       if (size(zwork1) < length*comm_size) then
          deallocate(zwork1,stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') &
                  'Error in comms_alltoall_complex_3: &
                  &deallocating zwork1 failed with code ',error
             call comms_abort
          end if
       end if
    end if
    if (.not. allocated(zwork1)) then
       allocate(zwork1(length*comm_size),stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
               &_complex_3: allocating zwork1 failed with code ',error
          call comms_abort
       end if
    end if

    ! Make a local copy of the data using an external subroutine to avoid
    ! the problem that the ranks of iwork1 and i_array differ
    call comms_copy(length*comm_size,z_array,zwork1)

    ! Perform alltoall
    call MPI_ALLTOALL(zwork1(1),length,complex_type,z_array(1,1,1),length, &
         complex_type,local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
            &_complex_3: MPI_ALLTOALL failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_alltoall_complex_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_alltoall_logical_1(l_array,length,comm)

    !=========================================================================!
    ! This subroutine is the logical vector form of comms_alltoall. Each node !
    ! sends length (default 1) logicals to every other node according to the  !
    ! order in which they appear in l_array.                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   l_array (in/out) : On entry, the data to be sent by this node.        !
    !                      On exit, the data received.                        !
    !   length (input)   : The (optional) number of elements to send to each  !
    !                      node.                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   logical_type      : The MPI type for logicals.                        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                                        !
    !   length >= 0                                                           !
    !   size(l_array) >= length*pub_total_num_nodes                           !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 22/7/03                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    logical, intent(inout) :: l_array(:)    ! On entry, data to be sent.
                                            ! On exit the data received
    integer, optional, intent(in) :: length ! The number of elements to send to/
                                            ! receive from each node
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                        ! Error flag
#endif
    integer :: local_length                 ! Local copy of length
    integer :: local_comm                   ! Local copy of comm argument
    integer :: comm_size                    ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) &
            'Error in comms_alltoall_logical_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_logical_1: length < 0'
       call comms_abort
    end if
    if (size(l_array) < local_length*comm_size) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_alltoall_logical_1: length exceeds array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) return

#ifdef MPI
    ! MPI needs separate send and receive buffers, so allocate workspace
    ! if required
    if (allocated(lwork1)) then
       if (size(lwork1) < local_length*comm_size) then
          deallocate(lwork1,stat=error)
          if (error /= 0) then
             if (pub_on_root) write(stdout,'(a,i3)') &
                  'Error in comms_alltoall_logical_1: &
                  &deallocating lwork1 failed with code ',error
             call comms_abort
          end if
       end if
    end if
    if (.not. allocated(lwork1)) then
       allocate(lwork1(local_length*comm_size),stat=error)
       if (error /= 0) then
          if (pub_on_root) write(stdout,'(a,i3)') &
               'Error in comms_alltoall_logical_1: &
               &allocating lwork1 failed with code ',error
          call comms_abort
       end if
    end if

    ! Copy data into local array
    lwork1(1:local_length*comm_size) = &
         l_array(1:local_length*comm_size)

    ! Perform alltoall
    call MPI_ALLTOALL(lwork1(1),local_length,logical_type,l_array(1), &
         local_length,logical_type,local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_alltoall&
            &_logical_1: MPI_ALLTOALL failed with code ',error
       call comms_abort
    end if
#endif

  end subroutine comms_alltoall_logical_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_allgather_integer_0(i_array_dest,i_array_src,length,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) integer scalar form of           !
    ! comms_allgather. It allgathers the scalar i_array_src to i_array_dest   !
    ! on all nodes, or optionally the first length double precision values    !
    ! pointed to by i_array_src are allgathered.                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   i_array_dest (out) : The array of data which will contain the result. !
    !   i_array_src (in)   : The source scalar.                               !
    !   length (input)   : The (optional) number of elements to be sent on    !
    !                      each node.                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type         : The MPI type for double precision integers.          !
    !   pub_total_num_nodes : The total number of nodes.                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                       (checked)        !
    !   length >= 0 (if present)                             (checked)        !
    !   size(i_array_src) >= length                          (checked)        !
    !   side(i_array_dest) >= length * pub_total_num_nodes   (checked)        !
    !   i_array_src and i_array_dest do not overlap          (NOT checked)    !
    !   length is the same on all nodes                      (NOT checked)    !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   1) If i_array_dest is larger than pub_total_num_nodes * length, the   !
    !      remaining elements remain unchanged, will not be zeroed.           !
    !   2) Note that this subroutine _requires_ length to be the same on all  !
    !      nodes, thus it's not suitable to work with the usual slab-distri-  !
    !      buted structures of onetep -- the number of slabs per core is      !
    !      not necessarily the same across cores in ONETEP.                   !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 08/11/2010, based on code by Jacek          !
    ! Dziedzic and Peter Haynes.                                              !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(inout) :: i_array_dest(:)! Receives the result
    integer, intent(in) :: i_array_src       ! The data to be gathered
    integer, optional, intent(in) :: length  ! The # of elements on each node
    integer, optional, intent(in) :: comm    ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    integer :: comm_size                     ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_allgather_integer_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_integer_0: length < 0'
       call comms_abort
    end if
    if (local_length * comm_size > size(i_array_dest)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_integer_0: combined length exceeds&
            & dest array size'
       call comms_abort
    end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) then
       return
    end if

#ifdef MPI
    ! Call MPI_ALLGATHER to transfer data
    call MPI_ALLGATHER(i_array_src,local_length,integer_type, &
         i_array_dest(1),local_length,integer_type,local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_allgather_integer_0: MPI_ALLGATHER failed with code ',&
            error
       call comms_abort
    end if
#else
    i_array_dest = i_array_src
#endif

  end subroutine comms_allgather_integer_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_allgather_integer_1(i_array_dest,i_array_src,length_src, &
       lengths_dest,displs_dest,comm)

    !=========================================================================!
    ! This subroutine is the integer scalar form of comms_allgather. It       !
    ! allgathers the scalar i_array_src to i_array_dest on all nodes, or      !
    ! optionally the first length integer values pointed to by i_array_src    !
    ! are allgathered.                                                        !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   i_array_dest (out) : The array of data which will contain the result. !
    !   i_array_src (in)   : The source array of data.                        !
    !   length (input)   : The (optional) number of elements to be sent on    !
    !                      each node.                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type     : The MPI type for integers.                         !
    !   pub_total_num_nodes : The total number of nodes.                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                       (checked)        !
    !   length >= 0 (if present)                             (checked)        !
    !   size(i_array_src) >= length                          (checked)        !
    !   side(i_array_dest) >= length * pub_total_num_nodes   (checked)        !
    !   i_array_src and i_array_dest do not overlap          (NOT checked)    !
    !   length is the same on all nodes                      (NOT checked)    !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   1) If i_array_dest is larger than pub_total_num_nodes * length, the   !
    !      remaining elements remain unchanged, will not be zeroed.           !
    !   2) Note that this subroutine _requires_ length to be the same on all  !
    !      nodes, thus it's not suitable to work with the usual slab-distri-  !
    !      buted structures of onetep -- the number of slabs per core is      !
    !      not necessarily the same across cores in ONETEP.                   !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 08/11/2010, based on code by Jacek          !
    ! Dziedzic and Peter Haynes.                                              !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(inout) :: i_array_dest(:) ! Receives the result
    integer, intent(in) :: i_array_src(:)     ! The data to be gathered
    integer, optional, intent(in) :: length_src     ! The # of elements on each node
    integer, optional, intent(in) :: displs_dest(:) ! The start indices for each node
    integer, optional, intent(in) :: lengths_dest(:)! The # of elements on each node
    integer, optional, intent(in) :: comm     ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                   ! Local copy of comm argument
    integer :: comm_size                    ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_allgather_integer_1: comms not &
            &initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length_src)) then
       local_length = length_src
    else
       local_length = size(i_array_src)
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_integer_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array_src)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_integer_1: length exceeds source array size'
       call comms_abort
    end if
    if (((present(lengths_dest)).and.(.not.present(displs_dest))).or. &
         ((.not.present(lengths_dest)).and.(present(displs_dest)))) then
       if (pub_on_root) then
          write(stdout,*) 'Error in comms_allgather_integer_3: inconsistent &
               &arguments passed to routine.'
          write(stdout,*) 'Both of displs_dest and lengths_dest (or &
               &neither) must be provided'
       end if
       call comms_abort
    end if
    !if (local_length * comm_size > size(i_array_dest)) then
    !   if (pub_on_root) write(stdout,*) &
    !        'Error in comms_allgather_integer_1: combined length exceeds&
    !        & dest array size'
    !   call comms_abort
    !end if

    ! If length is less than 1, no communication necessary
    ! ndmh: unless doing an allgatherv
    if ((local_length < 1).and.(.not.present(lengths_dest))) then
       return
    end if

#ifdef MPI
    ! Call MPI_ALLGATHER to transfer data
    if (present(lengths_dest)) then
       call MPI_ALLGATHERV(i_array_src(1),local_length,integer_type, &
            i_array_dest(1),lengths_dest(1),displs_dest(1),integer_type, &
            local_comm,error)
    else
       call MPI_ALLGATHER(i_array_src(1),local_length,integer_type, &
            i_array_dest(1),local_length,integer_type,local_comm,error)
    end if
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_allgather_integer_1: MPI_ALLGATHER(V) failed with code ',&
            error
       call comms_abort
    end if
#else
    i_array_dest = i_array_src
#endif

  end subroutine comms_allgather_integer_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_allgather_integer_2(i_array_dest,i_array_src,length_src, &
       lengths_dest,displs_dest,comm)

    !=========================================================================!
    ! This subroutine is the integer matrix form of comms_allgather. It       !
    ! allgathers the first length integers of i_array_src, or the whole       !
    ! array if length is absent, to i_array_dest on all nodes.                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   i_array_dest (out) : The array of data which will contain the result. !
    !   i_array_src (in)   : The source array of data.                        !
    !   length (input)   : The (optional) number of elements to be sent on    !
    !                      each node.                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   integer_type     : The MPI type for integers.                         !
    !   pub_total_num_nodes : The total number of nodes.                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                       (checked)        !
    !   length >= 0 (if present)                             (checked)        !
    !   size(i_array_src) >= length                          (checked)        !
    !   side(i_array_dest) >= length * pub_total_num_nodes   (checked)        !
    !   i_array_src and i_array_dest do not overlap          (NOT checked)    !
    !   length is the same on all nodes                      (NOT checked)    !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   1) If i_array_dest is larger than pub_total_num_nodes * length, the   !
    !      remaining elements remain unchanged, will not be zeroed.           !
    !   2) Note that this subroutine _requires_ length to be the same on all  !
    !      nodes, thus it's not suitable to work with the usual slab-distri-  !
    !      buted structures of onetep -- the number of slabs per core is      !
    !      not necessarily the same across cores in ONETEP.                   !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 08/11/2010, based on code by Jacek          !
    ! Dziedzic and Peter Haynes.                                              !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(inout) :: i_array_dest(:,:) ! Receives the result
    integer, intent(in) :: i_array_src(:,:)     ! The data to be gathered
    integer, optional, intent(in) :: length_src     ! The # of elements on each node
    integer, optional, intent(in) :: displs_dest(:) ! The start indices for each node
    integer, optional, intent(in) :: lengths_dest(:)! The # of elements on each node
    integer, optional, intent(in) :: comm     ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    integer :: comm_size                     ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_allgather_integer_2: comms not &
            &initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length_src)) then
       local_length = length_src
    else
       local_length = size(i_array_src)
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_integer_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array_src)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_integer_2: length exceeds source array size'
       call comms_abort
    end if
    if (((present(lengths_dest)).and.(.not.present(displs_dest))).or. &
         ((.not.present(lengths_dest)).and.(present(displs_dest)))) then
       if (pub_on_root) then
          write(stdout,*) 'Error in comms_allgather_integer_2: inconsistent &
               &arguments passed to routine.'
          write(stdout,*) 'Both of displs_dest and lengths_dest (or &
               &neither) must be provided'
       end if
       call comms_abort
    end if
    !if (local_length * comm_size > size(i_array_dest)) then
    !   if (pub_on_root) write(stdout,*) &
    !        'Error in comms_allgather_integer_2: combined length exceeds&
    !        & dest array size'
    !   call comms_abort
    !end if

    ! If length is less than 1, no communication necessary
    ! ndmh: unless doing an allgatherv
    if ((local_length < 1).and.(.not.present(lengths_dest))) then
       return
    end if

#ifdef MPI
    ! Call MPI_ALLGATHER to transfer data
    if (present(lengths_dest)) then
       call MPI_ALLGATHERV(i_array_src(1,1),local_length,integer_type, &
            i_array_dest(1,1),lengths_dest(1),displs_dest(1),integer_type, &
            local_comm,error)
    else
       call MPI_ALLGATHER(i_array_src(1,1),local_length,integer_type, &
            i_array_dest(1,1),local_length,integer_type,local_comm,error)
    end if
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_allgather_integer_2: MPI_ALLGATHER failed with &
            &code',error
       call comms_abort
    end if
#else
    i_array_dest = i_array_src
#endif

  end subroutine comms_allgather_integer_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_allgather_integer_3(i_array_dest,i_array_src,length_src, &
       lengths_dest,displs_dest,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real tensor form of           !
    ! comms_allgather. It allgathers the first length reals of d_array_src,   !
    ! or the whole array if length is absent, to d_array_dest on all nodes.   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   d_array_dest (out) : The array of data which will contain the result. !
    !   d_array_src (in)   : The source array of data.                        !
    !   length (input)   : The (optional) number of elements to be sent on    !
    !                      each node.                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   pub_total_num_nodes : The total number of nodes.                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                       (checked)        !
    !   length >= 0 (if present)                             (checked)        !
    !   size(d_array_src) >= length                          (checked)        !
    !   side(d_array_dest) >= length * pub_total_num_nodes   (checked)        !
    !   d_array_src and d_array_dest do not overlap          (NOT checked)    !
    !   length is the same on all nodes                      (NOT checked)    !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   1) If d_array_dest is larger than pub_total_num_nodes * length, the   !
    !      remaining elements remain unchanged, will not be zeroed.           !
    !   2) Note that this subroutine _requires_ length to be the same on all  !
    !      nodes, thus it's not suitable to work with the usual slab-distri-  !
    !      buted structures of onetep -- the number of slabs per core is      !
    !      not necessarily the same across cores in ONETEP.                   !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 08/11/2010, based on code by Jacek          !
    ! Dziedzic and Peter Haynes.                                              !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(inout) :: i_array_dest(:,:,:)   ! Receives the result
    integer, intent(in) :: i_array_src(:,:,:)       ! The data to be gathered
    integer, optional, intent(in) :: length_src     ! The # of elements on each node
    integer, optional, intent(in) :: displs_dest(:) ! The start indices for each node
    integer, optional, intent(in) :: lengths_dest(:)! The # of elements on each node
    integer, optional, intent(in) :: comm     ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    integer :: comm_size                     ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_allgather_integer_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length_src)) then
       local_length = length_src
    else
       local_length = size(i_array_src)
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_integer_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(i_array_src)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_integer_3: length exceeds source array size'
       call comms_abort
    end if
    if (((present(lengths_dest)).and.(.not.present(displs_dest))).or. &
         ((.not.present(lengths_dest)).and.(present(displs_dest)))) then
       if (pub_on_root) then
          write(stdout,*) 'Error in comms_allgather_integer_3: inconsistent &
               &arguments passed to routine.'
          write(stdout,*) 'Both of displs_dest and lengths_dest (or &
               &neither) must be provided'
       end if
       call comms_abort
    end if
    !if (local_length * comm_size > size(i_array_dest)) then
    !   if (pub_on_root) write(stdout,*) &
    !        'Error in comms_allgather_integer_3: combined length exceeds&
    !        & dest array size'
    !   call comms_abort
    !end if

    ! If length is less than 1, no communication necessary
    ! ndmh: unless doing an allgatherv
    if ((local_length < 1).and.(.not.present(lengths_dest))) then
       return
    end if

#ifdef MPI
    ! Call MPI_ALLGATHER to transfer data
    if (present(lengths_dest)) then
       call MPI_ALLGATHERV(i_array_src(1,1,1),local_length,integer_type, &
            i_array_dest(1,1,1),lengths_dest(1),displs_dest(1),integer_type, &
            local_comm,error)
    else
       call MPI_ALLGATHER(i_array_src(1,1,1),local_length,integer_type, &
            i_array_dest(1,1,1),local_length,integer_type,local_comm,error)
    end if
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_allgather_integer_3: MPI_ALLGATHER failed with code ',&
            error
       call comms_abort
    end if
#else
    i_array_dest = i_array_src
#endif


  end subroutine comms_allgather_integer_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_allgather_real_0(d_array_dest,d_array_src,length,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real scalar form of           !
    ! comms_allgather. It allgathers the scalar d_array_src to d_array_dest   !
    ! on all nodes, or optionally the first length double precision values    !
    ! pointed to by d_array_src are allgathered.                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   d_array_dest (out) : The array of data which will contain the result. !
    !   d_array_src (in)   : The source scalar.                               !
    !   length (input)   : The (optional) number of elements to be sent on    !
    !                      each node.                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   pub_total_num_nodes : The total number of nodes.                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                       (checked)        !
    !   length >= 0 (if present)                             (checked)        !
    !   size(d_array_src) >= length                          (checked)        !
    !   side(d_array_dest) >= length * pub_total_num_nodes   (checked)        !
    !   d_array_src and d_array_dest do not overlap          (NOT checked)    !
    !   length is the same on all nodes                      (NOT checked)    !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   1) If d_array_dest is larger than pub_total_num_nodes * length, the   !
    !      remaining elements remain unchanged, will not be zeroed.           !
    !   2) Note that this subroutine _requires_ length to be the same on all  !
    !      nodes, thus it's not suitable to work with the usual slab-distri-  !
    !      buted structures of onetep -- the number of slabs per core is      !
    !      not necessarily the same across cores in ONETEP.                   !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 28/07/2010 based on code by Peter Haynes.    !
    !=========================================================================!

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: d_array_dest(:) ! Receives the result
    real(kind=DP), intent(in) :: d_array_src        ! The data to be gathered
    integer, optional, intent(in) :: length  ! The # of elements on each node
    integer, optional, intent(in) :: comm     ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    integer :: comm_size                     ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_allgather_real_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional argument
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_real_0: length < 0'
       call comms_abort
    end if
    !if (local_length * comm_size > size(d_array_dest)) then
    !   if (pub_on_root) write(stdout,*) &
    !        'Error in comms_allgather_real_0: combined length exceeds&
    !        & dest array size'
    !   call comms_abort
    !end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) then
       return
    end if

#ifdef MPI
    ! Call MPI_ALLGATHER to transfer data
    call MPI_ALLGATHER(d_array_src,local_length,real_type, &
         d_array_dest(1),local_length,real_type,local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_allgather_real_0: MPI_ALLGATHER failed with code ',&
            error
       call comms_abort
    end if
#else
    d_array_dest = d_array_src
#endif

  end subroutine comms_allgather_real_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_allgather_real_1(d_array_dest,d_array_src,length_src, &
       lengths_dest,displs_dest,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real vector form of           !
    ! comms_allgather. It allgathers the first length reals of d_array_src,   !
    ! or the whole array if length is absent, to d_array_dest on all nodes.   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   d_array_dest (out) : The array of data which will contain the result. !
    !   d_array_src (in)   : The source array of data.                        !
    !   length (input)   : The (optional) number of elements to be sent on    !
    !                      each node.                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   pub_total_num_nodes : The total number of nodes.                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                       (checked)        !
    !   length >= 0 (if present)                             (checked)        !
    !   size(d_array_src) >= length                          (checked)        !
    !   side(d_array_dest) >= length * pub_total_num_nodes   (checked)        !
    !   d_array_src and d_array_dest do not overlap          (NOT checked)    !
    !   length is the same on all nodes                      (NOT checked)    !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   1) If d_array_dest is larger than pub_total_num_nodes * length, the   !
    !      remaining elements remain unchanged, will not be zeroed.           !
    !   2) Note that this subroutine _requires_ length to be the same on all  !
    !      nodes, thus it's not suitable to work with the usual slab-distri-  !
    !      buted structures of onetep -- the number of slabs per core is      !
    !      not necessarily the same across cores in ONETEP.                   !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 28/07/2010 based on code by Peter Haynes.    !
    !=========================================================================!

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: d_array_dest(:) ! Receives the result
    real(kind=DP), intent(in) :: d_array_src(:)     ! The data to be gathered
    integer, optional, intent(in) :: length_src     ! The # of elements on each node
    integer, optional, intent(in) :: displs_dest(:) ! The start indices for each node
    integer, optional, intent(in) :: lengths_dest(:)! The # of elements on each node
    integer, optional, intent(in) :: comm     ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    integer :: comm_size                     ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_allgather_real_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length_src)) then
       local_length = length_src
    else
       local_length = size(d_array_src)
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_real_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array_src)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_real_1: length exceeds source array size'
       call comms_abort
    end if
    if (((present(lengths_dest)).and.(.not.present(displs_dest))).or. &
         ((.not.present(lengths_dest)).and.(present(displs_dest)))) then
       if (pub_on_root) then
          write(stdout,*) 'Error in comms_allgather_real_1: inconsistent &
               &arguments passed to routine.'
          write(stdout,*) 'Both of displs_dest and lengths_dest (or &
               &neither) must be provided'
       end if
       call comms_abort
    end if
    !if (local_length * comm_size > size(d_array_dest)) then
    !   if (pub_on_root) write(stdout,*) &
    !        'Error in comms_allgather_real_1: combined length exceeds&
    !        & dest array size'
    !   call comms_abort
    !end if

    ! If length is less than 1, no communication necessary
    ! ndmh: unless doing an allgatherv
    if ((local_length < 1).and.(.not.present(lengths_dest))) then
       return
    end if

#ifdef MPI
    ! Call MPI_ALLGATHER to transfer data
    if (present(lengths_dest)) then
       call MPI_ALLGATHERV(d_array_src(1),local_length,real_type, &
            d_array_dest(1),lengths_dest(1),displs_dest(1),real_type, &
            local_comm,error)
    else
       call MPI_ALLGATHER(d_array_src(1),local_length,real_type, &
            d_array_dest(1),local_length,real_type,local_comm,error)
    end if
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_allgather_real_1: MPI_ALLGATHER failed with code ',&
            error
       call comms_abort
    end if
#else
    d_array_dest = d_array_src
#endif

  end subroutine comms_allgather_real_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_allgather_real_2(d_array_dest,d_array_src,length_src, &
       lengths_dest,displs_dest,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real matrix form of           !
    ! comms_allgather. It allgathers the first length reals of d_array_src,   !
    ! or the whole array if length is absent, to d_array_dest on all nodes.   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   d_array_dest (out) : The array of data which will contain the result. !
    !   d_array_src (in)   : The source array of data.                        !
    !   length (input)   : The (optional) number of elements to be sent on    !
    !                      each node.                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   pub_total_num_nodes : The total number of nodes.                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                       (checked)        !
    !   length >= 0 (if present)                             (checked)        !
    !   size(d_array_src) >= length                          (checked)        !
    !   side(d_array_dest) >= length * pub_total_num_nodes   (checked)        !
    !   d_array_src and d_array_dest do not overlap          (NOT checked)    !
    !   length is the same on all nodes                      (NOT checked)    !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   1) If d_array_dest is larger than pub_total_num_nodes * length, the   !
    !      remaining elements remain unchanged, will not be zeroed.           !
    !   2) Note that this subroutine _requires_ length to be the same on all  !
    !      nodes, thus it's not suitable to work with the usual slab-distri-  !
    !      buted structures of onetep -- the number of slabs per core is      !
    !      not necessarily the same across cores in ONETEP.                   !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 28/07/2010 based on code by Peter Haynes.    !
    !=========================================================================!

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: d_array_dest(:,:) ! Receives the result
    real(kind=DP), intent(in) :: d_array_src(:,:)     ! The data to be gathered
    integer, optional, intent(in) :: length_src     ! The # of elements on each node
    integer, optional, intent(in) :: displs_dest(:) ! The start indices for each node
    integer, optional, intent(in) :: lengths_dest(:)! The # of elements on each node
    integer, optional, intent(in) :: comm     ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    integer :: comm_size                     ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_allgather_real_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length_src)) then
       local_length = length_src
    else
       local_length = size(d_array_src)
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_real_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array_src)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_real_2: length exceeds source array size'
       call comms_abort
    end if
    if (((present(lengths_dest)).and.(.not.present(displs_dest))).or. &
         ((.not.present(lengths_dest)).and.(present(displs_dest)))) then
       if (pub_on_root) then
          write(stdout,*) 'Error in comms_allgather_real_2: inconsistent &
               &arguments passed to routine.'
          write(stdout,*) 'Both of displs_dest and lengths_dest (or &
               &neither) must be provided'
       end if
       call comms_abort
    end if
    !if (local_length * comm_size > size(d_array_dest)) then
    !   if (pub_on_root) write(stdout,*) &
    !        'Error in comms_allgather_real_2: combined length exceeds&
    !        & dest array size'
    !   call comms_abort
    !end if

    ! If length is less than 1, no communication necessary
    ! ndmh: unless doing an allgatherv
    if ((local_length < 1).and.(.not.present(lengths_dest))) then
       return
    end if

#ifdef MPI
    ! Call MPI_ALLGATHER to transfer data
    if (present(lengths_dest)) then
       call MPI_ALLGATHERV(d_array_src(1,1),local_length,real_type, &
            d_array_dest(1,1),lengths_dest(1),displs_dest(1),real_type, &
            local_comm,error)
    else
       call MPI_ALLGATHER(d_array_src(1,1),local_length,real_type, &
            d_array_dest(1,1),local_length,real_type,local_comm,error)
    end if

    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_allgather_real_2: MPI_ALLGATHER failed with code ',&
            error
       call comms_abort
    end if
#else
    d_array_dest = d_array_src
#endif

  end subroutine comms_allgather_real_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_allgather_real_3(d_array_dest,d_array_src,length_src, &
       lengths_dest,displs_dest,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) real tensor form of           !
    ! comms_allgather. It allgathers the first length reals of d_array_src,   !
    ! or the whole array if length is absent, to d_array_dest on all nodes.   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   d_array_dest (out) : The array of data which will contain the result. !
    !   d_array_src (in)   : The source array of data.                        !
    !   length (input)   : The (optional) number of elements to be sent on    !
    !                      each node.                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   real_type         : The MPI type for double precision reals.          !
    !   pub_total_num_nodes : The total number of nodes.                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                       (checked)        !
    !   length >= 0 (if present)                             (checked)        !
    !   size(d_array_src) >= length                          (checked)        !
    !   side(d_array_dest) >= length * pub_total_num_nodes   (checked)        !
    !   d_array_src and d_array_dest do not overlap          (NOT checked)    !
    !   length is the same on all nodes                      (NOT checked)    !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   1) If d_array_dest is larger than pub_total_num_nodes * length, the   !
    !      remaining elements remain unchanged, will not be zeroed.           !
    !   2) Note that this subroutine _requires_ length to be the same on all  !
    !      nodes, thus it's not suitable to work with the usual slab-distri-  !
    !      buted structures of onetep -- the number of slabs per core is      !
    !      not necessarily the same across cores in ONETEP.                   !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 28/07/2010 based on code by Peter Haynes.    !
    !=========================================================================!

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: d_array_dest(:,:,:) ! Receives the result
    real(kind=DP), intent(in) :: d_array_src(:,:,:)     ! The data to be gathered
    integer, optional, intent(in) :: length_src  ! The # of elements on each node
    integer, optional, intent(in) :: displs_dest(:) ! The start indices for each node
    integer, optional, intent(in) :: lengths_dest(:)! The # of elements on each node
    integer, optional, intent(in) :: comm        ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    integer :: comm_size                     ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_allgather_real_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length_src)) then
       local_length = length_src
    else
       local_length = size(d_array_src)
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_real_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(d_array_src)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_real_3: length exceeds source array size'
       call comms_abort
    end if
    if (((present(lengths_dest)).and.(.not.present(displs_dest))).or. &
         ((.not.present(lengths_dest)).and.(present(displs_dest)))) then
       if (pub_on_root) then
          write(stdout,*) 'Error in comms_allgather_real_3: inconsistent &
               &arguments passed to routine.'
          write(stdout,*) 'Both of displs_dest and lengths_dest (or &
               &neither) must be provided'
       end if
       call comms_abort
    end if
    !if (local_length * comm_size > size(d_array_dest)) then
    !   if (pub_on_root) write(stdout,*) &
    !        'Error in comms_allgather_real_3: combined length exceeds&
    !        & dest array size'
    !   call comms_abort
    !end if

    ! If length is less than 1, no communication necessary
    ! ndmh: unless doing an allgatherv
    if ((local_length < 1).and.(.not.present(lengths_dest))) then
       return
    end if

#ifdef MPI
    ! Call MPI_ALLGATHER to transfer data
    if (present(lengths_dest)) then
       call MPI_ALLGATHERV(d_array_src(1,1,1),local_length,real_type, &
            d_array_dest(1,1,1),lengths_dest(1),displs_dest(1),real_type, &
            local_comm,error)
    else
       call MPI_ALLGATHER(d_array_src(1,1,1),local_length,real_type, &
            d_array_dest(1,1,1),local_length,real_type,local_comm,error)
    end if
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_allgather_real_3: MPI_ALLGATHER(V) failed with code ',&
            error
       call comms_abort
    end if
#else
    d_array_dest = d_array_src
#endif

  end subroutine comms_allgather_real_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_allgather_complex_0(z_array_dest,z_array_src,length,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex scalar form of        !
    ! comms_allgather. It allgathers the scalar z_array_src to z_array_dest   !
    ! on all nodes, or optionally the first length double precision values    !
    ! pointed to by z_array_src are allgathered.                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   z_array_dest (out) : The array of data which will contain the result. !
    !   z_array_src (in)   : The source scalar.                               !
    !   length (input)   : The (optional) number of elements to be sent on    !
    !                      each node.                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator          : The communicator corresponding to active nodes!
    !                           in the current distribution strategy.         !
    !   complex_type          : The MPI type for double precision complexs.   !
    !   pub_total_num_nodes  : The total number of nodes.                     !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                       (checked)        !
    !   length >= 0 (if present)                             (checked)        !
    !   size(z_array_src) >= length                          (checked)        !
    !   side(z_array_dest) >= length * pub_total_num_nodes   (checked)        !
    !   z_array_src and z_array_dest do not overlap          (NOT checked)    !
    !   length is the same on all nodes                      (NOT checked)    !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   1) If z_array_dest is larger than pub_total_num_nodes * length, the   !
    !      remaining elements remain unchanged, will not be zeroed.           !
    !   2) Note that this subroutine _requires_ length to be the same on all  !
    !      nodes, thus it's not suitable to work with the usual slab-distri-  !
    !      buted structures of onetep -- the number of slabs per core is      !
    !      not necessarily the same across cores in ONETEP.                   !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 08/11/2010, based on code by Jacek          !
    ! Dziedzic and Peter Haynes.                                              !
    !=========================================================================!

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout) :: z_array_dest(:) ! Receives the result
    complex(kind=DP), intent(in) :: z_array_src        ! The data to be gathered
    integer, optional, intent(in) :: length  ! The # of elements on each node
    integer, optional, intent(in) :: comm     ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    integer :: comm_size                     ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_allgather_complex_0: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional argument
    if (present(length)) then
       local_length = length
    else
       local_length = 1
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_complex_0: length < 0'
       call comms_abort
    end if
    !if (local_length * comm_size > size(z_array_dest)) then
    !   if (pub_on_root) write(stdout,*) &
    !        'Error in comms_allgather_complex_0: combined length exceeds&
    !        & dest array size'
    !   call comms_abort
    !end if

    ! If length is less than 1, no communication necessary
    if (local_length < 1) then
       return
    end if

#ifdef MPI
    ! Call MPI_ALLGATHER to transfer data
    call MPI_ALLGATHER(z_array_src,local_length,complex_type, &
         z_array_dest(1),local_length,complex_type,local_comm,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_allgather_complex_0: MPI_ALLGATHER failed with code ',&
            error
       call comms_abort
    end if
#else
    z_array_dest = z_array_src
#endif

  end subroutine comms_allgather_complex_0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_allgather_complex_1(z_array_dest,z_array_src,length_src, &
       lengths_dest,displs_dest,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex vector form of           !
    ! comms_allgather. It allgathers the first length complexs of z_array_src,   !
    ! or the whole array if length is absent, to z_array_dest on all nodes.   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   z_array_dest (out) : The array of data which will contain the result. !
    !   z_array_src (in)   : The source array of data.                        !
    !   length (input)   : The (optional) number of elements to be sent on    !
    !                      each node.                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator      : The communicator corresponding to active nodes    !
    !                       in the current distribution strategy.             !
    !   complex_type         : The MPI type for double precision complexs.          !
    !   pub_total_num_nodes : The total number of nodes.                      !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                       (checked)        !
    !   length >= 0 (if present)                             (checked)        !
    !   size(z_array_src) >= length                          (checked)        !
    !   side(z_array_dest) >= length * pub_total_num_nodes   (checked)        !
    !   z_array_src and z_array_dest do not overlap          (NOT checked)    !
    !   length is the same on all nodes                      (NOT checked)    !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   1) If z_array_dest is larger than pub_total_num_nodes * length, the   !
    !      remaining elements remain unchanged, will not be zeroed.           !
    !   2) Note that this subroutine _requires_ length to be the same on all  !
    !      nodes, thus it's not suitable to work with the usual slab-distri-  !
    !      buted structures of onetep -- the number of slabs per core is      !
    !      not necessarily the same across cores in ONETEP.                   !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 08/11/2010, based on code by Jacek          !
    ! Dziedzic and Peter Haynes.                                              !
    !=========================================================================!

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout) :: z_array_dest(:) ! Receives the result
    complex(kind=DP), intent(in) :: z_array_src(:)     ! The data to be gathered
    integer, optional, intent(in) :: length_src     ! The # of elements on each node
    integer, optional, intent(in) :: displs_dest(:) ! The start indices for each node
    integer, optional, intent(in) :: lengths_dest(:)! The # of elements on each node
    integer, optional, intent(in) :: comm     ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    integer :: comm_size                     ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_allgather_complex_1: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length_src)) then
       local_length = length_src
    else
       local_length = size(z_array_src)
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_complex_1: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array_src)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_complex_1: length exceeds source array size'
       call comms_abort
    end if
    if (((present(lengths_dest)).and.(.not.present(displs_dest))).or. &
         ((.not.present(lengths_dest)).and.(present(displs_dest)))) then
       if (pub_on_root) then
          write(stdout,*) 'Error in comms_allgather_complex_1: inconsistent &
               &arguments passed to routine.'
          write(stdout,*) 'Both of displs_dest and lengths_dest (or &
               &neither) must be provided'
       end if
       call comms_abort
    end if
    !if (local_length * comm_size > size(z_array_dest)) then
    !   if (pub_on_root) write(stdout,*) &
    !        'Error in comms_allgather_complex_1: combined length exceeds&
    !        & dest array size'
    !   call comms_abort
    !end if

    ! If length is less than 1, no communication necessary
    ! ndmh: unless doing an allgatherv
    if ((local_length < 1).and.(.not.present(lengths_dest))) then
       return
    end if

#ifdef MPI
    ! Call MPI_ALLGATHER to transfer data
    if (present(lengths_dest)) then
       call MPI_ALLGATHERV(z_array_src(1),local_length,complex_type, &
            z_array_dest(1),lengths_dest(1),displs_dest(1),complex_type, &
            local_comm,error)
    else
       call MPI_ALLGATHER(z_array_src(1),local_length,complex_type, &
            z_array_dest(1),local_length,complex_type,local_comm,error)
    end if
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_allgather_complex_1: MPI_ALLGATHER(V) failed with code ',&
            error
       call comms_abort
    end if
#else
    z_array_dest = z_array_src
#endif

  end subroutine comms_allgather_complex_1

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_allgather_complex_2(z_array_dest,z_array_src,length_src, &
       lengths_dest,displs_dest,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex matrix form of        !
    ! comms_allgather. It allgathers the first length complexs of z_array_src,!
    ! or the whole array if length is absent, to z_array_dest on all nodes.   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   z_array_dest (out) : The array of data which will contain the result. !
    !   z_array_src (in)   : The source array of data.                        !
    !   length (input)   : The (optional) number of elements to be sent on    !
    !                      each node.                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator          : The communicator corresponding to active nodes!
    !                           in the current distribution strategy.         !
    !   complex_type          : The MPI type for double precision complexs.   !
    !   pub_total_num_nodes   : The total number of nodes.                    !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                       (checked)        !
    !   length >= 0 (if present)                             (checked)        !
    !   size(z_array_src) >= length                          (checked)        !
    !   side(z_array_dest) >= length * pub_total_num_nodes   (checked)        !
    !   z_array_src and z_array_dest do not overlap          (NOT checked)    !
    !   length is the same on all nodes                      (NOT checked)    !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   1) If z_array_dest is larger than pub_total_num_nodes * length, the   !
    !      remaining elements remain unchanged, will not be zeroed.           !
    !   2) Note that this subroutine _requires_ length to be the same on all  !
    !      nodes, thus it's not suitable to work with the usual slab-distri-  !
    !      buted structures of onetep -- the number of slabs per core is      !
    !      not necessarily the same across cores in ONETEP.                   !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 08/11/2010, based on code by Jacek          !
    ! Dziedzic and Peter Haynes.                                              !
    !=========================================================================!

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout) :: z_array_dest(:,:) ! Receives the result
    complex(kind=DP), intent(in) :: z_array_src(:,:)     ! The data to be gathered
    integer, optional, intent(in) :: length_src  ! The # of elements on each node
    integer, optional, intent(in) :: displs_dest(:) ! The start indices for each node
    integer, optional, intent(in) :: lengths_dest(:)! The # of elements on each node
    integer, optional, intent(in) :: comm        ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    integer :: comm_size                     ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_allgather_complex_2: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length_src)) then
       local_length = length_src
    else
       local_length = size(z_array_src)
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_complex_2: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array_src)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_complex_2: length exceeds source array size'
       call comms_abort
    end if
    if (((present(lengths_dest)).and.(.not.present(displs_dest))).or. &
         ((.not.present(lengths_dest)).and.(present(displs_dest)))) then
       if (pub_on_root) then
          write(stdout,*) 'Error in comms_allgather_complex_2: inconsistent &
               &arguments passed to routine.'
          write(stdout,*) 'Both of displs_dest and lengths_dest (or &
               &neither) must be provided'
       end if
       call comms_abort
    end if
    !if (local_length * comm_size > size(z_array_dest)) then
    !   if (pub_on_root) write(stdout,*) &
    !        'Error in comms_allgather_complex_2: combined length exceeds&
    !        & dest array size'
    !   call comms_abort
    !end if

    ! If length is less than 1, no communication necessary
    ! ndmh: unless doing an allgatherv
    if ((local_length < 1).and.(.not.present(lengths_dest))) then
       return
    end if

#ifdef MPI
    ! Call MPI_ALLGATHER to transfer data
    if (present(lengths_dest)) then
       call MPI_ALLGATHERV(z_array_src(1,1),local_length,complex_type, &
            z_array_dest(1,1),lengths_dest(1),displs_dest(1),complex_type, &
            local_comm,error)
    else
       call MPI_ALLGATHER(z_array_src(1,1),local_length,complex_type, &
            z_array_dest(1,1),local_length,complex_type,local_comm,error)
    end if
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_allgather_complex_2: MPI_ALLGATHER(V) failed with code ',&
            error
       call comms_abort
    end if
#else
    z_array_dest = z_array_src
#endif

  end subroutine comms_allgather_complex_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_allgather_complex_3(z_array_dest,z_array_src,length_src, &
       lengths_dest,displs_dest,comm)

    !=========================================================================!
    ! This subroutine is the (double precision) complex tensor form of        !
    ! comms_allgather. It allgathers the first length complexs of z_array_src,!
    ! or the whole array if length is absent, to z_array_dest on all nodes.   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   z_array_dest (out) : The array of data which will contain the result. !
    !   z_array_src (in)   : The source array of data.                        !
    !   length (input)   : The (optional) number of elements to be sent on    !
    !                      each node.                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   pub_comms_initialised : Flag to indicate comms has been initialised.  !
    !   communicator          : The communicator corresponding to active nodes!
    !                           in the current distribution strategy.         !
    !   complex_type          : The MPI type for double precision complexs.   !
    !   pub_total_num_nodes   : The total number of nodes.                    !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   pub_comms_initialised = .true.                       (checked)        !
    !   length >= 0 (if present)                             (checked)        !
    !   size(z_array_src) >= length                          (checked)        !
    !   side(z_array_dest) >= length * pub_total_num_nodes   (checked)        !
    !   z_array_src and z_array_dest do not overlap          (NOT checked)    !
    !   length is the same on all nodes                      (NOT checked)    !
    !-------------------------------------------------------------------------!
    ! Caveats:                                                                !
    !   1) If z_array_dest is larger than pub_total_num_nodes * length, the   !
    !      remaining elements remain unchanged, will not be zeroed.           !
    !   2) Note that this subroutine _requires_ length to be the same on all  !
    !      nodes, thus it's not suitable to work with the usual slab-distri-  !
    !      buted structures of onetep -- the number of slabs per core is      !
    !      not necessarily the same across cores in ONETEP.                   !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 08/11/2010, based on code by Jacek          !
    ! Dziedzic and Peter Haynes.                                              !
    !=========================================================================!

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout) :: z_array_dest(:,:,:) ! Receives the result
    complex(kind=DP), intent(in) :: z_array_src(:,:,:)     ! The data to be gathered
    integer, optional, intent(in) :: length_src  ! The # of elements on each node
    integer, optional, intent(in) :: displs_dest(:) ! The start indices for each node
    integer, optional, intent(in) :: lengths_dest(:)! The # of elements on each node
    integer, optional, intent(in) :: comm        ! The communicator to use

    ! Local variables
#ifdef MPI
    integer :: error                         ! Error flag
#endif
    integer :: local_length                  ! Local copy of length argument
    integer :: local_comm                    ! Local copy of comm argument
    integer :: comm_size                     ! Size of communicator

    ! Check comms has been initialised
    if (.not. pub_comms_initialised) then
       write(stdout,*) 'Error in comms_allgather_complex_3: comms not initialised'
       call comms_abort
    end if

    ! Make local copies of optional arguments
    if (present(length_src)) then
       local_length = length_src
    else
       local_length = size(z_array_src)
    end if
    if (present(comm)) then
       local_comm = comm
       comm_size = 1
#ifdef MPI
       call MPI_COMM_SIZE(local_comm,comm_size,error)
#endif
    else
       local_comm = communicator
       comm_size = pub_total_num_nodes
    end if

    ! Check arguments
    if (local_length < 0) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_complex_3: length < 0'
       call comms_abort
    end if
    if (local_length > size(z_array_src)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_allgather_complex_3: length exceeds source array size'
       call comms_abort
    end if
    if (((present(lengths_dest)).and.(.not.present(displs_dest))).or. &
         ((.not.present(lengths_dest)).and.(present(displs_dest)))) then
       if (pub_on_root) then
          write(stdout,*) 'Error in comms_allgather_real_2: inconsistent &
               &arguments passed to routine.'
          write(stdout,*) 'Both of displs_dest and lengths_dest (or &
               &neither) must be provided'
       end if
       call comms_abort
    end if
    !if (local_length * comm_size > size(z_array_dest)) then
    !   if (pub_on_root) write(stdout,*) &
    !        'Error in comms_allgather_complex_3: combined length exceeds&
    !        & dest array size'
    !   call comms_abort
    !end if

    ! If length is less than 1, no communication necessary
    ! ndmh: unless doing an allgatherv
    if ((local_length < 1).and.(.not.present(lengths_dest))) then
       return
    end if

#ifdef MPI
    ! Call MPI_ALLGATHER to transfer data
    if (present(lengths_dest)) then
       call MPI_ALLGATHERV(z_array_src(1,1,1),local_length,complex_type, &
            z_array_dest(1,1,1),lengths_dest(1),displs_dest(1),complex_type, &
            local_comm,error)
    else
       call MPI_ALLGATHER(z_array_src(1,1,1),local_length,complex_type, &
            z_array_dest(1,1,1),local_length,complex_type,local_comm,error)
    end if
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') &
            'Error in comms_allgather_complex_3: MPI_ALLGATHER(V) failed with code ',&
            error
       call comms_abort
    end if
#else
    z_array_dest = z_array_src
#endif

  end subroutine comms_allgather_complex_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_free_send_stack

    !=========================================================================!
    ! This subroutine checks outstanding non-blocking sends to see which have !
    ! completed. If the stack is full, the routine waits until a handle       !
    ! becomes available.                                                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !    None                                                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   max_handles       : The handle stack size.                            !
    !   num_send_handles  : The number of outstanding sends.                  !
    !   send_stack        : The stack of handles.                             !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! error : Used as an error flag for MPI routines.                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 13/11/03                                       !
    !=========================================================================!

    implicit none

#ifdef MPI
    ! Local variables
    integer, parameter :: WARN_THRESH = 1000 ! Stack full warning threshold
    integer :: istack,jstack                 ! Stack loop variables
    integer :: num_freed                     ! Number of completed entries
    integer :: error                         ! Error flag
    integer, save :: num_warns = -1          ! Stack full warning counter
    ! Space for MPI status structure - bug in some MPI implementations can
    ! cause crashes when MPI_STATUSES_IGNORE is used
    integer, dimension(MPI_STATUS_SIZE,max_handles) :: status
#endif

    ! If the stack is empty, return
    if (num_send_handles == 0) return

#ifdef MPI
    ! If the stack is full, we need to wait for one entry in stack to complete
    ! and remove this entry
    if (num_send_handles == max_handles) then

       ! Print warning if necessary
       if (num_warns == -1) then
          if (pub_on_root) write(stdout,'(a)') &
               &'WARNING in comms_free_send_stack: stack full'
       else if (num_warns == WARN_THRESH) then
          if (pub_on_root) write(stdout,'(a,i5,a)') &
               &'WARNING in comms_free_send_stack: stack full',num_warns, &
               &' times since last warning'
          num_warns = -1
       end if
       num_warns = num_warns + 1

       ! Wait for something to complete
       call MPI_WAITANY(max_handles,send_stack,istack,MPI_STATUS_IGNORE,error)
       if (error /= MPI_SUCCESS) then
          if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_free_send&
               &_stack: MPI_WAITANY failed with code ',error
          call comms_abort
       end if
       send_stack(istack:max_handles-1) = send_stack(istack+1:max_handles)
       num_send_handles = num_send_handles - 1
    end if

    ! If the stack is empty, return (only if max_handles == 1)
    if (num_send_handles == 0) then
       return
    end if

    ! Test all outstanding sends to check whether they have completed
    ! Could use MPI_STATUSES_IGNORE instead of status below, but this can
    ! cause some MPI implementations to crash
    call MPI_TESTSOME(num_send_handles,send_stack,num_freed,completed_list, &
         status,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_free_send&
            &_stack: MPI_TESTSOME failed with code ',error
       call comms_abort
    end if
    if (num_freed == MPI_UNDEFINED) then
       if (pub_on_root) write(stdout,'(a)') 'Error in comms_free_send_stack: &
            &corrupted send stack'
       call comms_abort
    end if

    ! Set up list of flags to indicate completed requests
    stack_flags = .false.
    do istack = 1,num_freed
       jstack = completed_list(istack)
       stack_flags(jstack) = .true.
    end do

    ! Shuffle stack to remove completed requests
    jstack = 0
    do istack = 1,num_send_handles
       if (.not. stack_flags(istack)) then
          jstack = jstack + 1
          send_stack(jstack) = send_stack(istack)
       end if
    end do
    num_send_handles = num_send_handles - num_freed
#endif

  end subroutine comms_free_send_stack

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_add_send_stack(handle)

    !=========================================================================!
    ! This subroutine adds a handle for a non-blocking send to the stack of   !
    ! outstanding sends.                                                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   handle (input)    : The handle to add to the stack                    !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   max_handles       : The handle stack size.                            !
    !   num_send_handles  : The number of outstanding sends.                  !
    !   send_stack        : The stack of handles.                             !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 13/11/03                                       !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: handle     ! The handle to add to the stack

    ! If the stack is full, we have a problem
    if (num_send_handles == max_handles) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_add_send_stack: &
            &stack full'
       call comms_abort
    end if

    ! Increment number of outstanding sends
    num_send_handles = num_send_handles + 1

    ! Add this handle to the stack
    send_stack(num_send_handles) = handle

  end subroutine comms_add_send_stack

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_empty_send_stack

    !=========================================================================!
    ! This subroutine empties the stack for outstanding non-blocking sends by !
    ! waiting until all have completed.                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   max_handles       : The handle stack size.                            !
    !   num_send_handles  : The number of outstanding sends.                  !
    !   send_stack        : The stack of handles.                             !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 13/11/03                                       !
    !=========================================================================!

    implicit none

#ifdef MPI
    ! Local variable
    integer :: error

    ! Space for MPI status structure - bug in some MPI implementations can
    ! cause crashes when MPI_STATUSES_IGNORE is used
    integer, dimension(MPI_STATUS_SIZE,max_handles) :: status
#endif

    ! If the stack is already empty, return
    if (num_send_handles == 0) return

    ! Wait for all outstanding sends to complete
#ifdef MPI
    ! Could use MPI_STATUSES_IGNORE instead of status below, but this can
    ! cause some MPI implementations to crash
    call MPI_WAITALL(num_send_handles,send_stack,status,error)
    if (error /= MPI_SUCCESS) then
       if (pub_on_root) write(stdout,'(a,i3)') 'Error in comms_empty_send&
            &_stack: MPI_WAITALL failed with code ',error
       call comms_abort
    end if
#endif

    ! Reset stack size to zero
    num_send_handles = 0

  end subroutine comms_empty_send_stack

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_copy_integer_2(length,source,dest)

    !=========================================================================!
    ! This performs a copy of integer data from one two-dimensional array to  !
    ! a one-dimensional array.                                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   length (input) : The number of elements to copy.                      !
    !   source (input) : The array from which the data will be copied.        !
    !   dest (output)  : The array into which the data will be copied.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   length >=0                                                            !
    !   size(source) >= length                                                !
    !   size(dest) >= length                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 21/6/04                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: length        ! The number of elements to copy
    integer, intent(in) :: source(:,:)   ! The array containing the data
    integer, intent(out) :: dest(:)      ! The array for the copy of the data

    ! Local variables
    integer :: i,j,k                     ! Loop variables
    integer :: size1                     ! Size of first dimension of source

    ! Check arguments
    if (length < 0) then
       if (pub_on_root) write(stdout,*) 'Error in comms_copy_integer_2:&
            & length < 0'
       call comms_abort
    end if
    if (length > size(source)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_copy_integer_2: length exceeds source array size'
       call comms_abort
    end if
    if (length > size(dest)) then
       if (pub_on_root) write(stdout,*) 'Error in comms_copy_integer_2: &
            &length exceeds destination array size'
       call comms_abort
    end if

    ! Copy the data
    size1 = size(source,1)
    k = 1
    do j=1,length / size1
       do i=1,size1
          dest(k) = source(i,j)
          k = k + 1
       end do
    end do
    do i=1,length - k + 1
       dest(k) = source(i,j)
       k = k + 1
    end do

  end subroutine comms_copy_integer_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_copy_integer_3(length,source,dest)

    !=========================================================================!
    ! This performs a copy of integer data from one three-dimensional array   !
    ! to a one-dimensional array.                                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   length (input) : The number of elements to copy.                      !
    !   source (input) : The array from which the data will be copied.        !
    !   dest (output)  : The array into which the data will be copied.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   length >=0                                                            !
    !   size(source) >= length                                                !
    !   size(dest) >= length                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 21/6/04                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: length        ! The number of elements to copy
    integer, intent(in) :: source(:,:,:) ! The array containing the data
    integer, intent(out) :: dest(:)      ! The array for the copy of the data

    ! Local variables
    integer :: i,j,k,l                   ! Loop variables
    integer :: size1,size2               ! Sizes of 1st two dimensions of source

    ! Check arguments
    if (length < 0) then
       if (pub_on_root) write(stdout,*) 'Error in comms_copy_integer_3:&
            & length < 0'
       call comms_abort
    end if
    if (length > size(source)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_copy_integer_3: length exceeds source array size'
       call comms_abort
    end if
    if (length > size(dest)) then
       if (pub_on_root) write(stdout,*) 'Error in comms_copy_integer_3: &
            &length exceeds destination array size'
       call comms_abort
    end if

    ! Copy the data
    size1 = size(source,1)
    size2 = size(source,2)
    l = 1
    do k=1,length / (size1 * size2)
       do j=1,size2
          do i=1,size1
             dest(l) = source(i,j,k)
             l = l + 1
          end do
       end do
    end do
    do j=1,length / size1 - (k-1) * size2
       do i=1,size1
          dest(l) = source(i,j,k)
          l = l + 1
       end do
    end do
    do i=1,length - l + 1
       dest(l) = source(i,j,k)
       l = l + 1
    end do

  end subroutine comms_copy_integer_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_copy_real_2(length,source,dest)

    !=========================================================================!
    ! This performs a copy of real data from one two-dimensional array to a   !
    ! one-dimensional array.                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   length (input) : The number of elements to copy.                      !
    !   source (input) : The array from which the data will be copied.        !
    !   dest (output)  : The array into which the data will be copied.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   length >=0                                                            !
    !   size(source) >= length                                                !
    !   size(dest) >= length                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 21/6/04                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: length            ! The number of elements to copy
    real(kind=DP), intent(in) :: source(:,:) ! The source of the data
    real(kind=DP), intent(out) :: dest(:)    ! The destination for the data

    ! Local variables
    integer :: i,j,k                     ! Loop variables
    integer :: size1                     ! Size of first dimension of source

    ! Check arguments
    if (length < 0) then
       if (pub_on_root) write(stdout,*) 'Error in comms_copy_real_2: length < 0'
       call comms_abort
    end if
    if (length > size(source)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_copy_real_2: length exceeds source array size'
       call comms_abort
    end if
    if (length > size(dest)) then
       if (pub_on_root) write(stdout,*) 'Error in comms_copy_real_2: &
            &length exceeds destination array size'
       call comms_abort
    end if

    ! Copy the data
    size1 = size(source,1)
    k = 1
    do j=1,length / size1
       do i=1,size1
          dest(k) = source(i,j)
          k = k + 1
       end do
    end do
    do i=1,length - k + 1
       dest(k) = source(i,j)
       k = k + 1
    end do

  end subroutine comms_copy_real_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_copy_real_3(length,source,dest)

    !=========================================================================!
    ! This performs a copy of real data from one three-dimensional array to a !
    ! one-dimensional array.                                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   length (input) : The number of elements to copy.                      !
    !   source (input) : The array from which the data will be copied.        !
    !   dest (output)  : The array into which the data will be copied.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   length >=0                                                            !
    !   size(source) >= length                                                !
    !   size(dest) >= length                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 21/6/04                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: length              ! The number of elements to copy
    real(kind=DP), intent(in) :: source(:,:,:) ! The source of the data
    real(kind=DP), intent(out) :: dest(:)      ! The destination for the data

    ! Local variables
    integer :: i,j,k,l                   ! Loop variables
    integer :: size1,size2               ! Sizes of 1st two dimensions of source

    ! Check arguments
    if (length < 0) then
       if (pub_on_root) write(stdout,*) 'Error in comms_copy_real_3: length < 0'
       call comms_abort
    end if
    if (length > size(source)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_copy_real_3: length exceeds source array size'
       call comms_abort
    end if
    if (length > size(dest)) then
       if (pub_on_root) write(stdout,*) 'Error in comms_copy_real_3: &
            &length exceeds destination array size'
       call comms_abort
    end if

    ! Copy the data
    size1 = size(source,1)
    size2 = size(source,2)
    l = 1
    do k=1,length / (size1 * size2)
       do j=1,size2
          do i=1,size1
             dest(l) = source(i,j,k)
             l = l + 1
          end do
       end do
    end do
    do j=1,length / size1 - (k-1) * size2
       do i=1,size1
          dest(l) = source(i,j,k)
          l = l + 1
       end do
    end do
    do i=1,length - l + 1
       dest(l) = source(i,j,k)
       l = l + 1
    end do

  end subroutine comms_copy_real_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_copy_complex_2(length,source,dest)

    !=========================================================================!
    ! This performs a copy of complex data from one two-dimensional array to  !
    ! a one-dimensional array.                                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   length (input) : The number of elements to copy.                      !
    !   source (input) : The array from which the data will be copied.        !
    !   dest (output)  : The array into which the data will be copied.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   length >=0                                                            !
    !   size(source) >= length                                                !
    !   size(dest) >= length                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 21/6/04                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: length               ! The number of elements to copy
    complex(kind=DP), intent(in) :: source(:,:) ! The source of the data
    complex(kind=DP), intent(out) :: dest(:)    ! The destination for the data

    ! Local variables
    integer :: i,j,k                     ! Loop variables
    integer :: size1                     ! Size of first dimension of source

    ! Check arguments
    if (length < 0) then
       if (pub_on_root) write(stdout,*) 'Error in comms_copy_complex_2:&
            & length < 0'
       call comms_abort
    end if
    if (length > size(source)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_copy_complex_2: length exceeds source array size'
       call comms_abort
    end if
    if (length > size(dest)) then
       if (pub_on_root) write(stdout,*) 'Error in comms_copy_complex_2: &
            &length exceeds destination array size'
       call comms_abort
    end if

    ! Copy the data
    size1 = size(source,1)
    k = 1
    do j=1,length / size1
       do i=1,size1
          dest(k) = source(i,j)
          k = k + 1
       end do
    end do
    do i=1,length - k + 1
       dest(k) = source(i,j)
       k = k + 1
    end do

  end subroutine comms_copy_complex_2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_copy_complex_3(length,source,dest)

    !=========================================================================!
    ! This performs a copy of complex data from one three-dimensional array   !
    ! to a one-dimensional array.                                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   length (input) : The number of elements to copy.                      !
    !   source (input) : The array from which the data will be copied.        !
    !   dest (output)  : The array into which the data will be copied.        !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   length >=0                                                            !
    !   size(source) >= length                                                !
    !   size(dest) >= length                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Peter Haynes, 21/6/04                                        !
    !=========================================================================!

    implicit none

    ! Arguments
    integer, intent(in) :: length                 ! Number of elements to copy
    complex(kind=DP), intent(in) :: source(:,:,:) ! The source of the data
    complex(kind=DP), intent(out) :: dest(:)      ! The destination for the data

    ! Local variables
    integer :: i,j,k,l                   ! Loop variables
    integer :: size1,size2               ! Sizes of 1st two dimensions of source

    ! Check arguments
    if (length < 0) then
       if (pub_on_root) write(stdout,*) 'Error in comms_copy_complex_3:&
            & length < 0'
       call comms_abort
    end if
    if (length > size(source)) then
       if (pub_on_root) write(stdout,*) &
            'Error in comms_copy_complex_3: length exceeds source array size'
       call comms_abort
    end if
    if (length > size(dest)) then
       if (pub_on_root) write(stdout,*) 'Error in comms_copy_complex_3: &
            &length exceeds destination array size'
       call comms_abort
    end if

    ! Copy the data
    size1 = size(source,1)
    size2 = size(source,2)
    l = 1
    do k=1,length / (size1 * size2)
       do j=1,size2
          do i=1,size1
             dest(l) = source(i,j,k)
             l = l + 1
          end do
       end do
    end do
    do j=1,length / size1 - (k-1) * size2
       do i=1,size1
          dest(l) = source(i,j,k)
          l = l + 1
       end do
    end do
    do i=1,length - l + 1
       dest(l) = source(i,j,k)
       l = l + 1
    end do

  end subroutine comms_copy_complex_3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine comms_killfile()

    !=========================================================================!
    ! This creates an empty file called killfile in the current directory     !
    !-------------------------------------------------------------------------!
    ! Arguments:                  none                                        !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, version 0.1, 26/11/04                         !
    !=========================================================================!

    implicit none

    integer :: ios, lunit

    if (pub_on_root) then
         lunit = 60
         open(unit=lunit,form='FORMATTED',status='UNKNOWN',&
         access='SEQUENTIAL',file='killfile',iostat=ios)

         write (lunit,*) 'Job terminated'
         close (lunit)

    end if

    return
  end subroutine comms_killfile

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module comms
#ifdef MPI
subroutine HPMPIINIT
#ifdef ACCELRYS
#include "initHPMPI.h"
!    if (hpmpi_ierr /= MPI_SUCCESS) then
!       write(stdout,'(a,i3)') 'Error in comms_init: &
!            &HP-MPI MPI_Initialized failed with code ',hpmpi_ierr
!       stop
!    end if
#endif
end subroutine HPMPIINIT
#endif
