 module openmpmod

  use lockmod

  !----------------------------------------------------!
  !                                                    !
  !   This is a suggested module to declare the OpenMP !
  !   functions and subroutine interfaces. This        !
  !   interface also declares a type for the OpenMP    !
  !   lock that hopefully is portable between systems. !
  !                                                    !
  !   A lock is in a user program is declared as:      !
  !      type(omp_lock_t) :: lock                      !
  !                                                    !
  !  The declaration is such that it is possible to    !
  !  construct allocatable arrays of locks             !
  !                                                    !
  !      type(omp_lock_t), allocatable :: lock_array(:)!
  !                                                    !
  ! Everything is provided as is and with no guarantees!
  !  Nils Smeds <smeds@pdc.kth.se>, 2000-08-28         !
  !                                                    !
  !----------------------------------------------------!

  private
  public :: init_openmp
  public :: omp_get_num_threads
  public :: omp_get_thread_num
  public :: openmp_split_array

 !############################################!
  INTERFACE

    subroutine omp_set_num_threads(NP)
      integer:: NP
    end subroutine omp_set_num_threads

    function omp_get_num_threads()
      integer:: omp_get_num_threads
    end function omp_get_num_threads

    function omp_get_max_threads()
      integer:: omp_get_max_threads
    end function omp_get_max_threads

    function omp_get_thread_num()
      integer:: omp_get_thread_num
    end function omp_get_thread_num

    function omp_get_num_procs()
      integer:: omp_get_num_procs
    end function omp_get_num_procs

    subroutine omp_set_dynamic(flag)
      logical:: flag
    end subroutine omp_set_dynamic

    function omp_get_dynamic()
      logical:: omp_get_dynamic
    end function omp_get_dynamic

    function omp_in_parallel()
      logical:: omp_in_parallel
    end function omp_in_parallel

    subroutine omp_set_nested(flag)
      logical:: flag
    end subroutine omp_set_nested

    function omp_get_nested()
      logical:: omp_get_nested
    end function omp_get_nested

    subroutine omp_init_lock(lock)
     use lockmod
      type(omp_lock_t) :: lock
    end subroutine omp_init_lock

    subroutine omp_destroy_lock(lock)
      use lockmod
      type(omp_lock_t) :: lock
    end subroutine omp_destroy_lock

    subroutine omp_unset_lock(lock)
      use lockmod
      type(omp_lock_t) :: lock
    end subroutine omp_unset_lock

    function omp_test_lock(lock)
      use lockmod
      type(omp_lock_t) :: lock
      logical :: omp_test_lock
    end function omp_test_lock

  end interface

 !############################################!


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
!**************************************************************************

  subroutine openmp_split_array(dimen,tabimin,tabimax)
  implicit none
    INTEGER(4), INTENT(IN) :: dimen
    INTEGER(4)             :: tabimin(:),tabimax(:),tabchunk(size(tabimin))
    INTEGER(4)             :: ipart,npart,pp,i

    if(dimen==0) then; 
      tabimin=0
      tabimax=0
      return
    endif

    npart = size(tabimin)
    if(npart>1000.or.npart<1)then
      write(*,*) 'more than 1000 omp processess or less than 1, something fishy'
      write(*,*) 'shape tabimin  : ', shape(tabimin)
      write(*,*) 'shape tabimax  : ', shape(tabimax)
      write(*,*) 'dimen          : ', dimen
      write(*,*) 'npart          : ', npart
      stop
    endif

#ifdef debug 
    write(*,*) 'splitting array of dimension : ' , dimen 
    write(*,*) 'in [x] openmp chuncks        : ' , npart
#endif

    if(npart>1)then
     tabchunk            =     0
     pp                  =     MOD(dimen,npart)
     tabchunk(1:npart)   =     (dimen-pp) / (npart)
      do i=1,pp
       tabchunk(i)=tabchunk(i)+1
      enddo
    else
     tabchunk(npart)     =   dimen
    endif

    if(sum(tabchunk)/=dimen)then
     write(*,*) 'error in split, tabchunk : ', tabchunk
     write(*,*) ' dimen                   : ', dimen
     stop
    endif

    DO ipart=1,npart
     if(ipart>1)then
      if(SUM(tabchunk(1:ipart-1))<dimen)then
        tabimin(ipart) = SUM(tabchunk(1:ipart-1)) + 1
      else
        tabimin(ipart) = 0
        tabimax(ipart) = 0
        cycle
      endif
     else
      tabimin(ipart) = 1
     endif
      tabimax(ipart) = SUM(tabchunk(1:ipart))
    ENDDO

    if(tabimax(npart)>dimen)then
     write(*,*) 'ERROR OPENMP SPLIT'
     write(*,*) 'last chunck min max : ', tabimin(npart),tabimax(npart)
     write(*,*) 'dimen : ', dimen
     write(*,*) 'max larger than dimen'
     stop 
    endif

  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

   subroutine init_openmp
   implicit none

     !$OMP PARALLEL PRIVATE(NTHREADS, TID)
        TID      = OMP_GET_THREAD_NUM()
        PROCS    = OMP_GET_NUM_PROCS()
        NTHREADS = OMP_GET_NUM_THREADS()
        MAXT     = OMP_GET_MAX_THREADS()
        INPAR    = OMP_IN_PARALLEL()
        DYNAMIC  = OMP_GET_DYNAMIC()
        NESTED   = OMP_GET_NESTED()
     !$OMP END PARALLEL

      IF (TID==0) THEN
        write(*,*) '========================================================'
        write(*,*)
        write(*,*) '  .............. INITIALIZE OPEN MP ..................  '
        write(*,*)
        write(*,*) 'Number of processors          = ', PROCS
        write(*,*) 'Number of threads             = ', NTHREADS
        write(*,*) 'Max threads                   = ', MAXT
        write(*,*) 'In parallel?                  = ', INPAR
        write(*,*) 'Dynamic threads enabled?      = ', DYNAMIC
        write(*,*) 'Nested parallelism supported? = ', NESTED
        write(*,*) '========================================================'
      END IF

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
!**************************************************************************
!**************************************************************************
!**************************************************************************
! 
!       subroutine dolooptest
!       implicit none
!       integer,parameter :: N=10
!       integer           :: i,j,k,l
!       real(16)          :: A(N),B(N),C(N),t
! 
!       call init_openmp
! 
!       C=0; A=(/( i+2, i=1,N )/); B=(/( i-1, i=1,N )/)
! 
!       call OMP_SET_NUM_THREADS(PROCS)
!      !$OMP PARALLEL SHARED(C) PRIVATE(I,J,t,l)
!      !$OMP DO SCHEDULE(DYNAMIC)
!       DO I = 1, N
!          t=1.
!          do J=1,4000
!           do l=1,1
!            t=t*sin(sqrt(2.d0/dble(j))+cos(-dble(j)))**8.d0
!           enddo
!          enddo
!          if(MOD(I,10)==0) write(*,*) I
!          C(I) = A(I) + B(I)
!       ENDDO
!      !$OMP END DO NOWAIT
!      !$OMP END PARALLEL
!       write(*,*) 'C = ' , C 
!  
!      return
!      end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************

end module
