module mpirout

 !--------------!
    use genvar
    use linalg
    use random
 !--------------!
   implicit none
   private
   public :: adim0r_
   public :: adim1s
   public :: distributed_memory_transpose_mat
   public :: distributed_memory_transpose_mat_c
   public :: distributed_memory_transpose_mat_r
   public :: mpi_stop
   public :: mpibarrier
   public :: mpibcast
   public :: mpibcastb_
   public :: mpibcastb___
   public :: mpibcastb____
   public :: mpigather_on_masternode
   public :: mpigather_on_masternode_r
   public :: mpisum
   public :: MPI_DOT_PRODUCT
   public :: MPI_DOT_PRODUCTc
   public :: MPI_DOT_PRODUCTr
   public :: scatter_it
   public :: scatter_it_r
   public :: split

!--------------------------------------------------------------------------!
INTERFACE mpigather_on_masternode
    MODULE PROCEDURE mpigather_on_masternode_c,mpigather_on_masternode_r
END INTERFACE
!--------------------------------------------------------------------------!
INTERFACE distributed_memory_transpose_mat
    MODULE PROCEDURE distributed_memory_transpose_mat_r,distributed_memory_transpose_mat_c
END INTERFACE
!--------------------------------------------------------------------------!
INTERFACE scatter_it
    MODULE PROCEDURE scatter_it_r,scatter_it_c
END INTERFACE
!--------------------------------------------------------------------------!
INTERFACE mpibcast
!BUG September 6th, 2012 : mpibcast__ and mpibcast_ leads to a confusion with mpibcastbcomp_ and mpibcastb_ respectively
!                            , and they are not used throughout the codes, so we removed it from the list
    MODULE PROCEDURE mpibcast___,mpibcast____,mpibcastb_,mpibcastb__,mpibcastb___,mpibcastb____, &
                   & mpibcast_2arrays_r,mpibcast_2arrays_c,mpibcast_no_disp_r,mpibcast_no_disp_c, &
                   & mpibcastbcomp_,mpibcastbcomp_2d,mpibcastb_2d
END INTERFACE
! !--------------------------------------------------------------------------!
! INTERFACE mpierror
!     MODULE PROCEDURE dim1,dim2,dim3,rdim1,rdim2,rdim3
! END INTERFACE
! !--------------------------------------------------------------------------!
! INTERFACE mpiaverage
!     MODULE PROCEDURE adim1,adim2,adim3,ardim0,ardim1,ardim2,ardim3
! END INTERFACE
!--------------------------------------------------------------------------!
INTERFACE mpisum
    MODULE PROCEDURE adim5r,adim5s,adim1s,adim2s,adim3s,ardim1s,ardim2s,ardim3s,adim4s,ardim4s,adim0,adim0r,adim0_,adim0r_,adim0i_
END INTERFACE
!--------------------------------------------------------------------------!
! INTERFACE mpi_max_func
!     MODULE PROCEDURE mpi_max_i,mpi_max_i_vec
! END INTERFACE
! !--------------------------------------------------------------------------!
! !SI CHAQUE CPU INDEPENDANT, CALCULE L ERREUR COMME DISPERSION STATISTIQUE
! INTERFACE mpiaverage_
!     MODULE PROCEDURE ardim1_mask,ardim2_mask,ardim3_mask,ardim4_mask,el_mask,ardim5_mask
! END INTERFACE
! !--------------------------------------------------------------------------!
! INTERFACE mpiswap
!     MODULE PROCEDURE mpiswap_i,mpiswap_r,mpiswap_c, &
!                    & mpiswap_ib,mpiswap_rb,mpiswap_cb, &
!                    & mpiswap_ibb,mpiswap_rbb,mpiswap_cbb
! END INTERFACE
! !--------------------------------------------------------------------------!
! INTERFACE mpi_read
!  MODULE PROCEDURE mpi_readr,mpi_readl,mpi_readi,mpi_readc,mpi_readrd,mpi_readid
! END INTERFACE
!--------------------------------------------------------------------------!
INTERFACE MPI_DOT_PRODUCT
 MODULE PROCEDURE MPI_DOT_PRODUCTi,MPI_DOT_PRODUCTr,MPI_DOT_PRODUCTc
END INTERFACE
!--------------------------------------------------------------------------!
! 
! 
! integer :: mpisend_incr

contains
! 
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! 
!   subroutine test_mpi_gather_scatter
!   implicit none
!   real(8) :: vec_(size2,10),vec(size2*10)
!   integer :: i
! 
!   if(rank==0)then
!   do i=1,10*size2
!    vec(i)=i
!   enddo
!   endif
! 
!   call scatter_it(vec,vec_(rank+1,:))
!   write(*,*) 'my rank, my vec', rank,vec_(rank+1,:)
!   stop
! 
!   end subroutine
! 
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  subroutine scatter_it_c(vec_to_scatter,target_vec)
   implicit none
   complex(8),dimension(:)  :: vec_to_scatter,target_vec

    if(size2==1.or.no_mpi)then
       target_vec=vec_to_scatter
       return
    endif
    call MPI_SCATTER(vec_to_scatter,size(target_vec),MPI_DOUBLE_COMPLEX,target_vec,size(target_vec), &
                       &  MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

  end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  subroutine scatter_it_r(vec_to_scatter,target_vec)
   implicit none
   real(8),dimension(:)    :: vec_to_scatter,target_vec

    if(size2==1.or.no_mpi)then
       target_vec=vec_to_scatter
       return
    endif
    call MPI_SCATTER(vec_to_scatter,size(target_vec),MPI_DOUBLE_PRECISION,target_vec,size(target_vec), &
                       &  MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
! 
!   subroutine test_mpi_distributed_transpose
!   implicit none
!   integer,parameter :: n1=4,n2_loc=2
!   integer           :: n2
!   integer           :: i,j,k,l,ii,jj
!   real(8)           :: mat_glob(n1,size2*n2_loc),transp(size2*n2_loc,n1)
!   real(8)           :: local_vec(n1*n2_loc)
! 
!   n2=size2*n2_loc
! 
!  !contains the full matrix
!   do i=1,n1
!    do j=1,n2
!     mat_glob(i,j)=i-j
!    enddo
!   enddo
! 
!  !contains the columns stored locally on a node 
!   do i=1,n1
!    do jj=1,n2_loc
!     j=rank*n2_loc+jj
!     local_vec(i+(jj-1)*n1)=i-j
!    enddo
!   enddo
! 
!   call mpibarrier
!   call distributed_memory_transpose_mat( local_vec , n1 , n2_loc , n2)  
!   call mpibarrier
! 
!   transp=transpose(mat_glob)
! 
!   call print_it(local_vec,n2,n1/size2)
! 
!   call mpibarrier
!   stop
! 
!   contains
! 
! 
!   subroutine print_it(mat_mpi,n1,n2)
!   implicit none
!   integer :: n1,n2,i
!   real(8) :: mat_mpi(n1,n2)
! 
!   do ii=1,size2
!    call mpibarrier
!     if(rank==ii-1)then
!      write(*,*) 'dumping info, my rank=',rank
!      do j=1,n2
!       jj = j + rank*n2
!       do i=1,n1
!         write(*,'(i3,i4,i4,f10.4,f10.4)') rank,i,jj,mat_mpi(i,j),transp(i,jj)
!       enddo
!      enddo
!     endif
!    enddo
! 
!    end subroutine
! 
!   end subroutine
! 
! 
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

    subroutine distributed_memory_transpose_mat_r( vec , nup , ndn , ndn_global)
    implicit none

   !-------------------------------------------------------!
   ! nup x ndn_global matrix distributed amongst the nodes !
   ! local chunk on one node is nup x ndn                  !
   !-------------------------------------------------------!

    real(8)    :: vec(:),vec2(size(vec))
    integer(4) :: nup,ndn
    integer(4) :: send_cnt(size2),send_disp(size2),rcv_cnt(size2),rcv_disp(size2),i,ndn_global
    integer(4) :: npack
    logical    :: verbose

    if(no_mpi)then
     write(*,*) 'ERROR distributed_memory_transpose_mat_r not yet coded for no_mpi'
     stop
    endif

    verbose=.false.

    if(size(vec)/=nup*ndn) then
      write(*,*) 'distributed memory matrix transpose, size do not match'
      write(*,*) 'input vector : ', size(vec)
      write(*,*) 'nup,ndn      : ', nup,ndn
      write(*,*) 'ndn global   : ', ndn_global 
      stop 'critical'
    endif
 
    if(mod(nup,size2)>0) then
       write(*,*) 'DANGER : using vec split among nodes require nup to be '
       write(*,*)         ' commensurate with the number of nodes '
       write(*,*) 'nup   = ', nup
       write(*,*) 'nodes = ', size2
       stop 'critical'
    endif

     if(verbose) write(*,*) 'call global transpose'
     call global_transpose(vec2,vec,nup,ndn)

     npack=nup/size2*ndn
     do i=1,size2
      send_cnt(i)  =       npack
      send_disp(i) = (i-1)*npack
      rcv_cnt(i)   =       npack
      rcv_disp(i)  = (i-1)*npack
     enddo

     if(verbose)write(*,*) 'broadcast all-to-all'
     call MPI_Alltoallv(vec2,send_cnt,send_disp,MPI_DOUBLE_PRECISION,vec,  &
               &   rcv_cnt,rcv_disp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

     if(verbose)write(*,*) 'local block transpose'
     call local_transpose_of_the_blocks(vec,vec2)

     if(verbose)write(*,*) 'global transpose'
     if(nup/size2*ndn_global/=size(vec)) then
       write(*,*) 'mpi transpose error, size of vector : ', size(vec)
       write(*,*) 'nup ndn ndn_global   : ', nup,ndn,ndn_global 
       write(*,*) 'number of nodes      : ', size2 
       write(*,*) 'nup/size2*ndn_global : ', nup/size2*ndn_global
       stop 'critical'
     endif
     call global_transpose(vec,vec2,nup/size2,ndn_global)

     return
     contains

       !--------------------!
       !--------------------!
       !--------------------!

      subroutine global_transpose(vec2,vec,nup,ndn)
      implicit none
      real(8) :: vec(nup,ndn),vec2(ndn,nup)
      integer :: nup,ndn
        vec2=transpose(vec) 
      end subroutine

       !--------------------!
       !--------------------!
       !--------------------!

      subroutine local_transpose_of_the_blocks(vec,vec2)
      implicit none
       real(8) :: vec(:),vec2(:)
       do i=1,size2
        call transpose_block_i(vec((i-1)*npack+1:i*npack),vec2((i-1)*npack+1:i*npack))
       enddo
      end subroutine

       !--------------------!
       !--------------------!
       !--------------------!

      subroutine transpose_block_i(vec,vec2)
      implicit none
       integer :: i1,i2,i
       real(8) :: vec(ndn,nup/size2),vec2(nup/size2,ndn)
       vec2=transpose(vec)
      end subroutine

       !--------------------!
       !--------------------!
       !--------------------!

    end subroutine


!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

    subroutine distributed_memory_transpose_mat_c( vec , nup , ndn , ndn_global)
    implicit none

   !-------------------------------------------------------!
   ! nup x ndn_global matrix distributed amongst the nodes !
   ! local chunk on one node is nup x ndn                  !
   !-------------------------------------------------------!

    complex(8) :: vec(:),vec2(size(vec))
    integer(4) :: nup,ndn
    integer(4) :: send_cnt(size2),send_disp(size2),rcv_cnt(size2),rcv_disp(size2),i,ndn_global
    integer(4) :: npack
    logical    :: verbose

    if(no_mpi)then
     write(*,*) 'ERROR distributed_memory_transpose_mat_c not yet coded for no_mpi'
     stop
    endif

    verbose=.false.

    if(size(vec)/=nup*ndn) then
      write(*,*) 'distributed memory matrix transpose, size do not match'
      write(*,*) 'input vector : ', size(vec)
      write(*,*) 'nup,ndn      : ', nup,ndn
      write(*,*) 'ndn global   : ', ndn_global 
      stop 'critical'
    endif
 
    if(mod(nup,size2)>0) then
       write(*,*) 'DANGER : using vec split among nodes require nup to be '
       write(*,*)         ' commensurate with the number of nodes '
       write(*,*) 'nup   = ', nup
       write(*,*) 'nodes = ', size2
       stop 'critical'
    endif

     if(verbose) write(*,*) 'call global transpose'
     call global_transpose(vec2,vec,nup,ndn)

     npack=nup/size2*ndn
     do i=1,size2
      send_cnt(i)  =       npack
      send_disp(i) = (i-1)*npack
      rcv_cnt(i)   =       npack
      rcv_disp(i)  = (i-1)*npack
     enddo

     if(verbose)write(*,*) 'broadcast all-to-all'
     call MPI_Alltoallv(vec2,send_cnt,send_disp,MPI_DOUBLE_COMPLEX,vec,  &
               &   rcv_cnt,rcv_disp,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)

     if(verbose)write(*,*) 'local block transpose'
     call local_transpose_of_the_blocks(vec,vec2)

     if(verbose)write(*,*) 'global transpose'
     if(nup/size2*ndn_global/=size(vec)) then
       write(*,*) 'mpi transpose error, size of vector : ', size(vec)
       write(*,*) 'nup ndn ndn_global   : ', nup,ndn,ndn_global 
       write(*,*) 'number of nodes      : ', size2 
       write(*,*) 'nup/size2*ndn_global : ', nup/size2*ndn_global
       stop 'critical'
     endif
     call global_transpose(vec,vec2,nup/size2,ndn_global)

     return
     contains

       !--------------------!
       !--------------------!
       !--------------------!

      subroutine global_transpose(vec2,vec,nup,ndn)
      implicit none
      complex(8) :: vec(nup,ndn),vec2(ndn,nup)
      integer    :: nup,ndn
        vec2=transpose(vec) 
      end subroutine

       !--------------------!
       !--------------------!
       !--------------------!

      subroutine local_transpose_of_the_blocks(vec,vec2)
      implicit none
       complex(8) :: vec(:),vec2(:)
       do i=1,size2
        call transpose_block_i(vec((i-1)*npack+1:i*npack),vec2((i-1)*npack+1:i*npack))
       enddo
      end subroutine

       !--------------------!
       !--------------------!
       !--------------------!

      subroutine transpose_block_i(vec,vec2)
      implicit none
       integer :: i1,i2,i
       complex(8) :: vec(ndn,nup/size2),vec2(nup/size2,ndn)
       vec2=transpose(vec)
      end subroutine

       !--------------------!
       !--------------------!
       !--------------------!

    end subroutine


!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
! 
! subroutine mpi_max_i(vec)
! implicit none
! integer    :: vec(:,:),vecm(size(vec,1),size(vec,2)),vecmm(size(vec,1),size(vec,2))
! integer    :: i,j,k,sizev
! 
! if(size2==1.or.no_mpi) then
!  return
! endif
! 
! vecmm=-99999999; vecm=0 ;sizev=size(vec,1)*size(vec,2)
! 
! do j=0,size2-1
!  if(rank==j) vecm=vec
!  call MPI_BCAST(vecm,sizev,MPI_INTEGER,j,MPI_COMM_WORLD,ierr)
!  where(vecm>vecmm) vecmm=vecm
! enddo
! 
! vec=vecmm
! 
! return
! end subroutine
! 
!    !-----------------------------!
! 
! subroutine mpi_max_i_vec(vec)
! implicit none
! integer    :: vec(:),vecm(size(vec)),vecmm(size(vec))
! integer    :: i,j,k,sizev
! 
! if(size2==1.or.no_mpi) return
! 
! vecmm=-99999999; vecm=0 ;sizev=size(vec)
! 
! do j=0,size2-1
!  if(rank==j) vecm=vec
!  call MPI_BCAST(vecm,sizev,MPI_INTEGER,j,MPI_COMM_WORLD,ierr)
!  where(vecm>vecmm) vecmm=vecm
! enddo
! 
! vec=vecmm
! 
! return
! end subroutine
! 
! 
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE split(dimen,tabimin,tabimax,tabchunk,nparts)

  ! GIVES THE RANGES TO SPLIT ARRAY OF DIMENSION dimen INTO nparts PIECES on each CPU
    implicit none
    INTEGER, INTENT(IN)           :: dimen
    INTEGER, INTENT(INOUT)        :: tabimin(:),tabimax(:),tabchunk(:)
    INTEGER, INTENT(IN), OPTIONAL :: nparts
    INTEGER                       :: ipart,npart,pp
    integer                       :: i

    if(dimen==0) then
     stop 'ERROR : dimension supplied to split is zero'
    endif

    npart = nproc

    !write(*,*) 'start splitting among nodes :', nproc

    IF(PRESENT(nparts)) npart = nparts

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

    !write(*,*)'rank,tabchunk : ' , rank, tabchunk

    if(sum(tabchunk)/=dimen)then
     write(*,*) 'error in split, tabchunk: ', tabchunk
     write(*,*) ' dimen                  : ', dimen
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
     write(*,*) 'ERROR MPI SPLIT'
     write(*,*) 'last chunck min max : ', tabimin(npart),tabimax(npart)
     write(*,*) 'dimen : ', dimen
     write(*,*) 'max larger than dimen'
     stop 
    endif

  END SUBROUTINE

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

 subroutine mpi_stop(string)
 character*(*) :: string
     write(*,*) ' my rank is : ', rank,size2
     write(*,*) ' I am done...'
     if(.not.no_mpi) call mpi_finalize(ierr)
     write(*,*) ' mpi termination successfull'
     write(*,*) ' message : ', string
     stop 
 end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
! 
!  subroutine test_mpi_dot
!  real(8)    :: v1(3),v2(3)
!  complex(8) :: vv1(3),vv2(3)
! 
!  call randomize_rank1(v1)
!  call randomize_rank1(v2)
!  call randomize_rank1(vv1)
!  call randomize_rank1(vv2)
! 
!  if(rank==0)then
!   write(*,*) 'v1,v2'
!   write(*,*) v1,v2
!   write(*,*) 'vv1,vv2'
!   write(*,*) vv1,vv2
!   write(*,*) 'dot prod'
!   write(*,*) DOT_PRODUCT(v1,v2)
!   write(*,*) DOT_PRODUCT(vv1,vv2)
!  endif
!  
!   write(*,*) 'mpi dot prod'
!   write(*,*) rank,MPI_DOT_PRODUCT(v1,v2)
!   write(*,*) rank,MPI_DOT_PRODUCT(vv1,vv2)
! 
!  call MPI_FINALIZE(ierr)
!  stop
!  end subroutine
! 
     !------------------------!

 integer function MPI_DOT_PRODUCTi(v1,v2,split)
 implicit none
 integer          :: v1(:),v2(:)
 integer          :: m
 integer          :: i,j,k,l
 logical,optional :: split

 k=min(size(v1),size(v2))
 MPI_DOT_PRODUCTi=0

 if((.not.MPIseparate.and..not.no_mpi).and.enable_mpi_dot)then
   if(size2==1)then
   !$OMP PARALLEL PRIVATE(i), SHARED(v1,v2,k), &
   !$OMP REDUCTION(+:MPI_DOT_PRODUCTi)
   !$OMP DO 
         do i=1,k
            MPI_DOT_PRODUCTi = MPI_DOT_PRODUCTi + v1(i)*v2(i)
         end do
   !$OMP END DO
   !$OMP END PARALLEL 
   else
    if(present(split))then
     if(split)then
      do i=1,k
        MPI_DOT_PRODUCTi=MPI_DOT_PRODUCTi+v1(i)*v2(i)
      enddo
      call MPI_ALLREDUCE(MPI_DOT_PRODUCTi,m,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      MPI_DOT_PRODUCTi=m
     return
     endif
    endif
     do i=rank+1,k,size2
      MPI_DOT_PRODUCTi=MPI_DOT_PRODUCTi+v1(i)*v2(i)
     enddo
     call MPI_ALLREDUCE(MPI_DOT_PRODUCTi,m,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
     MPI_DOT_PRODUCTi=m
    endif
 else
    MPI_DOT_PRODUCTi=DOT_PRODUCT(v1,v2) 
 endif
 end function

     !------------------------!

 real(8) function MPI_DOT_PRODUCTr(v1,v2,split)
 implicit none
 real(8)          :: v1(:),v2(:)
 real(8)          :: m
 integer          :: i,j,k,l
 logical,optional :: split

 k=min(size(v1),size(v2))
 MPI_DOT_PRODUCTr=0

 if((.not.MPIseparate.and..not.no_mpi).and.enable_mpi_dot)then
  if(k==0) stop 'error mpi dot product'
  if(size2==1)then
  !$OMP PARALLEL PRIVATE(i), SHARED(v1,v2,k), &
  !$OMP REDUCTION(+:MPI_DOT_PRODUCTr)
  !$OMP DO
        do i=1,k
           MPI_DOT_PRODUCTr = MPI_DOT_PRODUCTr + v1(i)*v2(i)
        end do
  !$OMP END DO
  !$OMP END PARALLEL
  else
    if(present(split))then
     if(split)then
      do i=1,k
        MPI_DOT_PRODUCTr=MPI_DOT_PRODUCTr+v1(i)*v2(i)
      enddo
      call MPI_ALLREDUCE(MPI_DOT_PRODUCTr,m,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      MPI_DOT_PRODUCTr=m
     return
     endif
    endif
    do i=rank+1,k,size2
     MPI_DOT_PRODUCTr=MPI_DOT_PRODUCTr+v1(i)*v2(i)
    enddo
    call MPI_ALLREDUCE(MPI_DOT_PRODUCTr,m,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    MPI_DOT_PRODUCTr=m
  endif
 else
  MPI_DOT_PRODUCTr=DOT_PRODUCT(v1,v2)
 endif
 end function

     !------------------------!

 complex(8) function MPI_DOT_PRODUCTc(v1,v2,split)
 implicit none
 complex(8)       :: v1(:),v2(:)
 complex(8)       :: m
 integer          :: i,j,k,l
 logical,optional :: split

 k=min(size(v1),size(v2))
 MPI_DOT_PRODUCTc=0

 if((.not.MPIseparate.and.no_mpi).and.enable_mpi_dot)then
   if(size2==1)then
   !$OMP PARALLEL PRIVATE(i), SHARED(v1,v2,k), &
   !$OMP REDUCTION(+:MPI_DOT_PRODUCTc)
   !$OMP DO
         do i=1,k
            MPI_DOT_PRODUCTc = MPI_DOT_PRODUCTc + conjg(v1(i))*v2(i)
         end do
   !$OMP END DO
   !$OMP END PARALLEL 
   else
     if(present(split))then
      if(split)then
       do i=1,k
         MPI_DOT_PRODUCTc=MPI_DOT_PRODUCTc+conjg(v1(i))*v2(i)
       enddo
       call MPI_ALLREDUCE(MPI_DOT_PRODUCTc,m,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
       MPI_DOT_PRODUCTc=m
      return
      endif
     endif
     do i=rank+1,k,size2
      MPI_DOT_PRODUCTc=MPI_DOT_PRODUCTc+conjg(v1(i))*v2(i)
     enddo
     call MPI_ALLREDUCE(MPI_DOT_PRODUCTc,m,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
     MPI_DOT_PRODUCTc=m
   endif
 else
  MPI_DOT_PRODUCTc=DOT_PRODUCT(v1,v2)
 endif
 end function

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine mpibarrier
    if(size2==1.or.no_mpi)return
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine efficient_gather
! !written by A.Camyani
! implicit none
! integer                                :: rank, num_proc, ierr, i, share
! real,dimension(16),target              :: a
! real,dimension(:,:),allocatable,target :: b
! real,dimension(:,:),pointer            :: c => null()
! 
!  if(no_mpi.or.size2==1) return
! 
!  a = (/ ( real(i) , i = 1, 16 ) /)
!  allocate(b(4,4))
!  b = reshape( a, (/ 4,4 /) )
! 
!  print'(a,i2,a)','I''am:',rank,' and b is:'
!  print'(4(2x,f5.2))', ( b(i,:) , i = 1, 4 )
!  print*,''
! 
!  share = size(b,2)/num_proc    ! Just caca, take care of uncommensurability
!  c => b( : , 1+share*rank:share*(rank+1) )
!  print'(a,i2,a)','I''am: ',rank,' and c is:'
!  print'(2(2x,f5.2))', ( c(i,:) , i = 1, 4 )
!  print*,''
!  c = c*real(rank+2)
!  call MPI_Gather(c,size(c),MPI_REAL,b,size(c),MPI_REAL,0,MPI_COMM_WORLD, ierr)
!  c => null()
! 
!  print'(a,i2,a)','I''am:',rank,' and b is:'
!  print'(4(2x,f5.2))', ( b(i,:) , i = 1, 4 )
!  print*,''
! 
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine test_gather
! implicit none
! integer                                :: i1,i2,i,v(2),w(2)
! real(8),dimension(22)                  :: a
! 
!  write(*,*) '-------------------'
!  write(*,*) '-------------------'
!  a=0.d0
!  if(rank==0)then
!  do i=1,10
!   a(rank*10+i)=dble(i)
!  enddo
!  else
!  do i=1,12
!   a(rank*10+i)=dble(i)
!  enddo
!  endif
!  v=[10,12]
!  w=[0,10]
!  if(rank==0) call efficient_gather_(a,v,w)
!  if(rank==1) call efficient_gather_(a,v,w)
! 
!  write(*,*) 'FINAL after gather is : ', a
!  write(*,*) '-------------------'
!  write(*,*) '-------------------'
! 
!  call MPI_finalize(ierr)
!  stop
! 
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine efficient_gather_(a,v,w)
! implicit none
! integer                                :: v(:),w(:),siz
! real(8),dimension(:),target            :: a
! real(8),dimension(:),pointer           :: c => null()
! 
!  if(no_mpi.or.size2==1) return
! 
!  if(v(rank+1)==0)then
!    write(*,*) '================================'
!    write(*,*) 'WARNING mpibcast___ empty buffer'
!    c => null()
!    siz=0
!  else
!    c => a(1+w(rank+1):w(rank+1)+v(rank+1))
!    siz=size(c)
!  endif
! 
!  where(v==0) w=0
! 
!  call MPI_AllGatherv(c,siz,MPI_DOUBLE_PRECISION,a,v,w,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
!  c => null()
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! 
 !---------------------!

subroutine mpibcastbcomp_2d(a,iii)
implicit none
complex(8),dimension(:,:)              :: a
integer,optional                       :: iii
integer                                :: jjj
  if(no_mpi.or.size2==1) return
  if(present(iii))then
   jjj=iii
  else
   jjj=0
  endif
  call MPI_BCAST(a,size(a,1)*size(a,2),MPI_DOUBLE_COMPLEX,jjj,MPI_COMM_WORLD,ierr)
end subroutine

 !---------------------!

subroutine mpibcastb_2d(a,iii)
implicit none
real(8),dimension(:,:)                 :: a
integer,optional                       :: iii
integer                                :: jjj
  if(no_mpi.or.size2==1) return
  if(present(iii))then
   jjj=iii
  else
   jjj=0
  endif
  call MPI_BCAST(a,size(a,1)*size(a,2),MPI_DOUBLE_PRECISION,jjj,MPI_COMM_WORLD,ierr)
end subroutine


 !---------------------!

subroutine mpibcastbcomp_(a,iii)
implicit none
complex(8),dimension(:)                :: a
integer,optional                       :: iii
integer                                :: jjj
  if(no_mpi.or.size2==1) return
  if(present(iii))then
   jjj=iii
  else
   jjj=0
  endif
  call MPI_BCAST(a(:),size(a),MPI_DOUBLE_COMPLEX,jjj,MPI_COMM_WORLD,ierr)
end subroutine

 !---------------------!

subroutine mpibcastb_(a,iii)
implicit none
real(8),dimension(:)                   :: a
integer,optional                       :: iii
integer                                :: jjj
  if(no_mpi.or.size2==1) return
  if(present(iii))then
   jjj=iii
  else
   jjj=0
  endif
  call MPI_BCAST(a(:),size(a),MPI_DOUBLE_PRECISION,jjj,MPI_COMM_WORLD,ierr)
end subroutine

 !---------------------!

subroutine mpibcastb__(a,iii)
implicit none
real(8)                                :: a
integer,optional                       :: iii
integer                                :: jjj
  if(no_mpi.or.size2==1) return
  if(present(iii))then
   jjj=iii
  else
   jjj=0
  endif
  call MPI_BCAST(a,1,MPI_DOUBLE_PRECISION,jjj,MPI_COMM_WORLD,ierr)
end subroutine

 !---------------------!

subroutine mpibcastb___(ii,iii)
implicit none
integer                                :: ii
integer,optional                       :: iii
integer                                :: jjj
  if(present(iii))then
   jjj=iii
  else
   jjj=0
  endif
  call MPI_BCAST(ii,1,MPI_INTEGER,jjj,MPI_COMM_WORLD,ierr)
end subroutine

 !---------------------!

subroutine mpibcastb____(a,iii)
implicit none
integer                                :: step
integer,dimension(:)                   :: a
integer,optional                       :: iii
integer                                :: jjj
   if(no_mpi.or.size2==1) return
   if(present(iii))then
    jjj=iii
   else
    jjj=0
   endif
   call MPI_BCAST(a(:),size(a),MPI_INTEGER,jjj,MPI_COMM_WORLD,ierr)
end subroutine

 !---------------------!

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine mpibcast_(a,step)
implicit none
integer                                :: i1,i2,step,j
real(8),dimension(:)                   :: a
 if(no_mpi.or.size2==1) return
 i1=1
 i2=step
 do i=0,size2-1
  call MPI_BCAST(a(i1:i2),i2-i1+1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
  if(i/=size2-2)then
   i1=i1+step
   i2=i2+step
  else
   i1=i1+step
   i2=size(a)
  endif
 enddo
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine mpibcast___(a,v,w)
implicit none
integer                                :: v(:),w(:),siz
real(8),dimension(:),target            :: a
real(8),dimension(:),pointer           :: c => null()

 if(size2==1.or.no_mpi) return

 if(v(rank+1)==0)then
   write(*,*) '================================'
   write(*,*) 'WARNING mpibcast___ empty buffer'
   c => null()
   siz=0
 else
   c => a(1+w(rank+1):w(rank+1)+v(rank+1))
   siz=size(c)
 endif

 where(v==0) w=0

 call MPI_AllGatherv(c,siz,MPI_DOUBLE_PRECISION,a,v,w,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
 c => null()
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine mpibcast_2arrays_r(a,source,v,w)
implicit none
integer                                :: v(:),w(:)
real(8),dimension(:)                   :: a,source
 if(size2==1.or.no_mpi)then
   a=source
   return
 endif
 call MPI_AllGatherv(source,size(source),MPI_DOUBLE_PRECISION,a,v,w,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine mpigather_on_masternode_r(a,source)
implicit none
real(8),dimension(:)   :: a,source
 if(size2==1.or.no_mpi)then
   a=source
   return
 endif
 call MPI_Gather(source,size(source),MPI_DOUBLE_PRECISION,a, &
              & size(source),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine mpigather_on_masternode_c(a,source)
implicit none
complex(8),dimension(:)  :: a,source
 if(size2==1.or.no_mpi)then
   a=source
   return
 endif
 call MPI_Gather(source,size(source),MPI_DOUBLE_COMPLEX,a,size(source), &
               & MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine mpibcast_no_disp_r(a,source,v)
implicit none
integer                                :: v(:),disp(size(v)),i
real(8),dimension(:)                   :: a,source
 if(size2==1.or.no_mpi)then
   a=source
   return
 endif
 disp=0
 do i=2,size2
  disp(i)=disp(i-1)+v(i-1)
 enddo
 call MPI_AllGatherv(source,size(source),MPI_DOUBLE_PRECISION,a,v,disp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine mpibcast_no_disp_c(a,source,v)
implicit none
integer                                   :: v(:),disp(size(v)),i
complex(8),dimension(:)                   :: a,source
 if(size2==1.or.no_mpi)then
   a=source
   return
 endif
 disp=0
 do i=2,size2
  disp(i)=disp(i-1)+v(i-1)
 enddo
 call MPI_AllGatherv(source,size(source),MPI_DOUBLE_COMPLEX,a,v,disp,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine mpibcast_2arrays_c(a,source,v,w)
implicit none
integer                                :: v(:),w(:)
complex(8),dimension(:)                :: a,source
 if(size2==1.or.no_mpi)then 
   a=source
   return
 endif
 call MPI_AllGatherv(source,size(source),MPI_DOUBLE_COMPLEX,a,v,w,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine mpibcast____(a,v,w)
implicit none
integer                                :: v(:),w(:),siz
complex(8),dimension(:),target         :: a
complex(8),dimension(:),pointer        :: c => null()

 if(size2==1.or.no_mpi)return

 if(v(rank+1)==0)then
   write(*,*) '================================'
   write(*,*) 'WARNING mpibcast___ empty buffer'
   c => null()
   siz=0
 else
   c => a(1+w(rank+1):w(rank+1)+v(rank+1))
   siz=size(c)
 endif

 where(v==0) w=0

 call MPI_AllGatherv(c,siz,MPI_DOUBLE_COMPLEX,a,v,w,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,ierr)
 c => null()
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
! 
! subroutine mpibcast__(a,step)
! implicit none
! integer                              :: i1,i2,step,j,i
! complex(8),dimension(:)              :: a
!  if(no_mpi.or.size2==1) return
!  i1=1
!  i2=step
!  do i=0,size2-1
!   call MPI_BCAST(a(i1:i2),i2-i1+1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
!   if(i/=size2-2)then
!    i1=i1+step
!    i2=i2+step
!   else
!    i1=i1+step
!    i2=size(a)
!   endif
!  enddo
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine mpi_readr(unit,number)
! implicit none
! integer :: unit,ier
! real(8)  :: number
!  if(rank==0) read(unit,*,end=66,err=66) number
!  66 continue
!  if(size2>1) call MPI_BCAST(number,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpi_readrd(unit,number,number2)
! implicit none
! integer :: unit
! real(8)  :: number,number2
!  if(rank==0) read(unit,*,end=66,err=66) number,number2
!  66 continue
!  if(no_mpi.or.size2==1) return
!  call MPI_BCAST(number,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(number2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpi_readi(unit,number)
! implicit none
! integer  :: unit
! integer  :: number
!  if(rank==0) read(unit,*,end=66,err=66) number
!  66 continue
!  if(size2>1.and..not.no_mpi) call MPI_BCAST(number,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpi_readid(unit,number,number2)
! implicit none
! integer  :: unit
! integer  :: number,number2
!  if(rank==0) read(unit,*,end=66,err=66) number,number2
!  66 continue
!  if(no_mpi.or.size2==1) return
!  call MPI_BCAST(number,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!  call MPI_BCAST(number2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpi_readl(unit,number)
! implicit none
! integer  :: unit
! logical  :: number
!  if(rank==0) read(unit,*,end=66,err=66) number
!  66 continue
!  if(size2>1.and..not.no_mpi) call MPI_BCAST(number,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpi_readc(unit,number)
! implicit none
! integer     :: unit
! complex(8)  :: number
!  if(rank==0) read(unit,*,end=66,err=66) number
!  66 continue
!  if(size2>1.and..not.no_mpi) call MPI_BCAST(number,1,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! 
! !swap array on a ring for the different nodes
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpiswap_ibb(array1,ii)
! implicit none
! integer :: array1,temp(0:size2-1)
! integer :: i,j,k,l,m,s1,s2,ii
!  if(no_mpi.or.size2==1) return
!  do i=0,size2-1
!   temp(i)=array1
!   CALL MPI_BCAST(temp(i),1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
!  enddo
!  if(ii==1)then !turn right
!    if(rank/=size2-1)then
!     array1=temp(rank+1)
!    else
!     array1=temp(0)
!    endif
!  else
!    if(rank/=0)then
!     array1=temp(rank-1)
!    else          !turn left
!     array1=temp(size2-1)
!    endif
!  endif
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpiswap_rbb(array1,ii)
! implicit none
! real(8) :: array1,temp(0:size2-1)
! integer :: i,j,k,l,m,s1,s2,ii
!  if(no_mpi.or.size2==1) return
!  do i=0,size2-1
!   temp(i)=array1
!   CALL MPI_BCAST(temp(i),1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
!  enddo
!  if(ii==1)then !turn right
!    if(rank/=size2-1)then
!     array1=temp(rank+1)
!    else
!     array1=temp(0)
!    endif
!  else
!    if(rank/=0)then
!     array1=temp(rank-1)
!    else          !turn left
!     array1=temp(size2-1)
!    endif
!  endif
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpiswap_cbb(array1,ii)
! implicit none
! complex(8) :: array1,temp(0:size2-1)
! integer :: i,j,k,l,m,s1,s2,ii
!  if(no_mpi.or.size2==1) return
!  do i=0,size2-1
!   temp(i)=array1
!   CALL MPI_BCAST(temp(i),1,MPI_DOUBLE_COMPLEX,i,MPI_COMM_WORLD,ierr)
!  enddo
!  if(ii==1)then !turn right
!    if(rank/=size2-1)then
!     array1=temp(rank+1)
!    else
!     array1=temp(0)
!    endif
!  else
!    if(rank/=0)then
!     array1=temp(rank-1)
!    else          !turn left
!     array1=temp(size2-1)
!    endif
!  endif
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpiswap_i(array1,ii)
! implicit none
! integer :: array1(:,:),temp(0:size2-1,size(array1(:,1)),size(array1(1,:)))
! integer :: i,j,k,l,m,s1,s2,ii
! if(no_mpi.or.size2==1) return
! s1=size(array1(:,1))
! s2=size(array1(1,:))
!  do i=0,size2-1
!   temp(i,:,:)=array1
!   CALL MPI_BCAST(temp(i,:,:),s1*s2,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
!  enddo
!  if(ii==1)then !turn right
!    if(rank/=size2-1)then
!     array1=temp(rank+1,:,:)
!    else          
!     array1=temp(0,:,:)
!    endif
!  else
!    if(rank/=0)then
!     array1=temp(rank-1,:,:)
!    else          !turn left
!     array1=temp(size2-1,:,:)
!    endif
!  endif
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpiswap_ib(array1,ii)
! implicit none
! integer :: array1(:),temp(0:size2-1,size(array1(:)))
! integer :: i,j,k,l,m,s1,s2,ii
!  if(no_mpi.or.size2==1) return
!  s1=size(array1(:))
!  do i=0,size2-1
!   temp(i,:)=array1
!   CALL MPI_BCAST(temp(i,:),s1,MPI_INTEGER,i,MPI_COMM_WORLD,ierr)
!  enddo
!  if(ii==1)then !turn right
!    if(rank/=size2-1)then
!     array1=temp(rank+1,:)
!    else
!     array1=temp(0,:)
!    endif
!  else
!    if(rank/=0)then
!     array1=temp(rank-1,:)
!    else          !turn left
!     array1=temp(size2-1,:)
!    endif
!  endif
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpiswap_r(array1,ii)
! implicit none
! real(8)  :: array1(:,:),temp(0:size2-1,size(array1(:,1)),size(array1(1,:)))
! integer :: i,j,k,l,m,s1,s2,ii
!  if(no_mpi.or.size2==1) return
!  s1=size(array1(:,1))
!  s2=size(array1(1,:))
!  do i=0,size2-1
!   temp(i,:,:)=array1
!   CALL MPI_BCAST(temp(i,:,:),s1*s2,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
!  enddo
!  if(ii==1)then !turn right
!    if(rank/=size2-1)then
!     array1=temp(rank+1,:,:)
!    else
!     array1=temp(0,:,:)
!    endif
!  else
!    if(rank/=0)then
!     array1=temp(rank-1,:,:)
!    else          !turn left
!     array1=temp(size2-1,:,:)
!    endif
!  endif
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpiswap_rb(array1,ii)
! implicit none
! real(8)  :: array1(:),temp(0:size2-1,size(array1(:)))
! integer :: i,j,k,l,m,s1,s2,ii
!  if(no_mpi.or.size2==1) return
!  s1=size(array1(:))
!  do i=0,size2-1
!   temp(i,:)=array1
!   CALL MPI_BCAST(temp(i,:),s1,MPI_DOUBLE_PRECISION,i,MPI_COMM_WORLD,ierr)
!  enddo
!  if(ii==1)then !turn right
!    if(rank/=size2-1)then
!     array1=temp(rank+1,:)
!    else
!     array1=temp(0,:)
!    endif
!  else
!    if(rank/=0)then
!     array1=temp(rank-1,:)
!    else          !turn left
!     array1=temp(size2-1,:)
!    endif
!  endif
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpiswap_c(array1,ii)
! implicit none
! complex(8) :: array1(:,:),temp(0:size2-1,size(array1(:,1)),size(array1(1,:)))
! integer    :: i,j,k,l,m,s1,s2,ii
!  if(no_mpi.or.size2==1) return
!  s1=size(array1(:,1))
!  s2=size(array1(1,:))
!  do i=0,size2-1
!   temp(i,:,:)=array1
!   CALL MPI_BCAST(temp(i,:,:),s1*s2,MPI_DOUBLE_COMPLEX,i,MPI_COMM_WORLD,ierr)
!  enddo
!  if(ii==1)then !turn right
!    if(rank/=size2-1)then
!     array1=temp(rank+1,:,:)
!    else
!     array1=temp(0,:,:)
!    endif
!  else
!    if(rank/=0)then
!     array1=temp(rank-1,:,:)
!    else          !turn left
!     array1=temp(size2-1,:,:)
!    endif
!  endif
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! subroutine mpiswap_cb(array1,ii)
! implicit none
! complex(8) :: array1(:),temp(0:size2-1,size(array1(:)))
! integer    :: i,j,k,l,m,s1,s2,ii
!  if(no_mpi.or.size2==1) return
!  s1=size(array1(:))
!  do i=0,size2-1
!   temp(i,:)=array1
!   CALL MPI_BCAST(temp(i,:),s1,MPI_DOUBLE_COMPLEX,i,MPI_COMM_WORLD,ierr)
!  enddo
!  if(ii==1)then !turn right
!    if(rank/=size2-1)then
!     array1=temp(rank+1,:)
!    else
!     array1=temp(0,:)
!    endif
!  else
!    if(rank/=0)then
!     array1=temp(rank-1,:)
!    else          !turn left
!     array1=temp(size2-1,:)
!    endif
!  endif
! end subroutine
! 
!  !-----------------------------!
!  !-----------------------------!
! 
! 
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine dim1(vec,Ivec)
! implicit none
! complex(8) :: vec(:),vecm(size(vec)),vec2m(size(vec)),Ivec(size(vec))
! integer    :: i,j,k,sizev
!  Ivec=0
!  vec2m=0
!  vecm=0
!  if(no_mpi.or.size2==1) return
!  sizev=size(vec)
!  call erase_divergence(vec)
!  call MPI_ALLREDUCE(vec,vecm,sizev,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr) 
!  call MPI_ALLREDUCE(vec**2,vec2m,sizev,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
!  vecm  = vecm/dble(size2)
!  vec2m = vec2m/dble(size2)
!  Ivec  = sqrt( abs (vecm**2-vec2m)/dble(size2) )
! return 
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine dim2(mat,Imat)
! implicit none
! complex(8) :: mat(:,:),matm(size(mat(:,1)),size(mat(1,:))), &
!                     & mat2m(size(mat(:,1)),size(mat(1,:))), & 
!                     &  Imat(size(mat(:,1)),size(mat(1,:)))
! integer    :: i,j,k,s1,s2
! 
! s1=size(mat(:,1))
! s2=size(mat(1,:))
! Imat=0.
! mat2m=0.
! matm=0.
! if(no_mpi.or.size2==1) return
! call erase_divergence(mat)
! call MPI_ALLREDUCE(mat,matm,s1*s2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
! call MPI_ALLREDUCE(mat**2,mat2m,s1*s2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
! matm=matm/dble(size2)
! mat2m=mat2m/dble(size2)
! Imat= sqrt (  abs (matm**2-mat2m) / dble(size2) )
! return
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine dim3(mat,Imat)
! implicit none
! complex(8) :: mat(:,:,:),matm(size(mat(:,1,1)),size(mat(1,:,1)),size(mat(1,1,:))), &
!             & mat2m(size(mat(:,1,1)),size(mat(1,:,1)),size(mat(1,1,:))), & 
!            &  Imat(size(mat(:,1,1)),size(mat(1,:,1)),size(mat(1,1,:)))
! integer :: i,j,k,s1,s2,s3
! 
! s1=size(mat(:,1,1))
! s2=size(mat(1,:,1))
! s3=size(mat(1,1,:))
! Imat=0.
! mat2m=0.
! matm=0.
! if(no_mpi.or.size2==1) return
! call erase_divergence(mat)
! call MPI_ALLREDUCE(mat,matm,s1*s2*s3,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
! call MPI_ALLREDUCE(mat**2,mat2m,s1*s2*s3,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
! matm=matm/dble(size2)
! mat2m=mat2m/dble(size2)
! Imat= sqrt (  abs (matm**2-mat2m) / dble(size2) )
! return
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine rdim1(vec,Ivec)
! implicit none
! real(8) :: vec(:),vecm(size(vec)),vec2m(size(vec)),Ivec(size(vec))
! integer :: i,j,k,sizev
!  Ivec=0
!  vec2m=0
!  vecm=0
!  if(no_mpi.or.size2==1) return
!  sizev=size(vec)
!  call erase_divergence(vec)
!  call MPI_ALLREDUCE(vec,vecm,sizev,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) 
!  call MPI_ALLREDUCE(vec**2,vec2m,sizev,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!  vecm=vecm/dble(size2)
!  vec2m=vec2m/dble(size2)
!  Ivec= sqrt( abs (vecm**2-vec2m)/dble(size2) )
! return
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine rdim2(mat,Imat)
! implicit none
! real(8) :: mat(:,:),matm(size(mat(:,1)),size(mat(1,:))), &
!                  & mat2m(size(mat(:,1)),size(mat(1,:))), & 
!                  &  Imat(size(mat(:,1)),size(mat(1,:)))
! integer :: i,j,k,s1,s2
! 
! s1=size(mat(:,1))
! s2=size(mat(1,:))
! Imat=0.
! mat2m=0.
! matm=0.
! if(no_mpi.or.size2==1) return
! call erase_divergence(mat)
! call MPI_ALLREDUCE(mat,matm,s1*s2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
! call MPI_ALLREDUCE(mat**2,mat2m,s1*s2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
! matm=matm/dble(size2)
! mat2m=mat2m/dble(size2)
! Imat= sqrt (  abs (matm**2-mat2m) / dble(size2) )
! return
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine rdim3(mat,Imat)
! implicit none
! real(8) :: mat(:,:,:),matm(size(mat(:,1,1)),size(mat(1,:,1)),size(mat(1,1,:))), &
!             & mat2m(size(mat(:,1,1)),size(mat(1,:,1)),size(mat(1,1,:))), & 
!            &  Imat(size(mat(:,1,1)),size(mat(1,:,1)),size(mat(1,1,:)))
! integer :: i,j,k,s1,s2,s3
! 
! s1=size(mat(:,1,1));s2=size(mat(1,:,1));s3=size(mat(1,1,:))
! Imat=0.;mat2m=0.;matm=0.
! if(no_mpi.or.size2==1) return
! call erase_divergence(mat)
! call MPI_ALLREDUCE(mat,matm,s1*s2*s3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
! call MPI_ALLREDUCE(mat**2,mat2m,s1*s2*s3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
! matm=matm/dble(size2)
! mat2m=mat2m/dble(size2)
! Imat= sqrt (  abs (matm**2-mat2m) / dble(size2) )
! return
! end subroutine
! 
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine adim0(vec,vecm)
implicit none
complex(8) :: vec,vecm
 if(size2==1.or.no_mpi) then
  vecm=vec
  return
 endif
 vecm=0
 call erase_divergence(vec)
 call MPI_ALLREDUCE(vec,vecm,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
 vecm=vecm/dble(size2)
return
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine adim0r(vec,vecm)
implicit none
real(8)    :: vec,vecm
if(size2==1.or.no_mpi) then
 vecm=vec
 return
endif
vecm=0
call erase_divergence(vec)
call MPI_ALLREDUCE(vec,vecm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
vecm=vecm/dble(size2)
return
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine adim0_(vec)
implicit none
complex(8) :: vec,vecm
if(size2==1.or.no_mpi) then
 return
endif
vecm=0
call erase_divergence(vec)
call MPI_ALLREDUCE(vec,vecm,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
vec=vecm
return
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine adim0r_(vec)
implicit none
real(8)    :: vec,vecm
if(size2==1.or.no_mpi) then
 return
endif
vecm=0
call erase_divergence(vec)
call MPI_ALLREDUCE(vec,vecm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
vec=vecm
return
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine adim0i_(vec)
implicit none
integer    :: vec,vecm
if(size2==1.or.no_mpi) then
 return
endif
vecm=0
call erase_divergence(vec)
call MPI_ALLREDUCE(vec,vecm,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
vec=vecm
return
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine adim1(vec,vecm2)
implicit none
complex(8) :: vec(:),vecm(size(vec)),vecm2(size(vec))
integer    :: sizev
if(size2==1.or.no_mpi) then
 vecm2=vec
 return
endif
vecm=0;sizev=size(vec)
call erase_divergence(vec)
call MPI_ALLREDUCE(vec,vecm,sizev,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
vecm2=vecm/dble(size2)
return
end subroutine

!*********************************************
!*********************************************
!*********************************************

subroutine adim2(mat,matm)
implicit none
complex(8) :: mat(:,:),matm(size(mat,1),size(mat,2))
integer    :: s1,s2
if(size2==1.or.no_mpi) then
 matm=mat
 return
endif
s1=size(mat(:,1));s2=size(mat(1,:));matm=0.
call erase_divergence(mat)
call MPI_ALLREDUCE(mat,matm,s1*s2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
matm=matm/dble(size2)
return
end subroutine

!*********************************************
!*********************************************
!*********************************************

subroutine adim3(mat,matm)
implicit none
complex(8) :: mat(:,:,:),matm(size(mat(:,1,1)),size(mat(1,:,1)),size(mat(1,1,:)))
integer    :: s1,s2,s3
if(size2==1.or.no_mpi)then
 matm=mat
 return
endif
s1=size(mat(:,1,1));s2=size(mat(1,:,1));s3=size(mat(1,1,:));matm=0.
call erase_divergence(mat)
call MPI_ALLREDUCE(mat,matm,s1*s2*s3,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
matm=matm/dble(size2)
return
end subroutine

!*********************************************
!*********************************************
!*********************************************

subroutine adim1s(vec)
implicit none
complex(8) :: vec(:),vecm(size(vec))
integer    :: sizev
if(size2==1.or.no_mpi)then
 return
endif
vecm=0;sizev=size(vec)
call erase_divergence(vec)
call MPI_ALLREDUCE(vec,vecm,sizev,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
vec=vecm
return
end subroutine

!*********************************************
!*********************************************
!*********************************************

subroutine adim2s(mat)
implicit none
complex(8) :: mat(:,:),matm(size(mat(:,1)),size(mat(1,:)))
integer    :: s1,s2
if(size2==1.or.no_mpi)then
 return
endif
s1=size(mat(:,1));s2=size(mat(1,:));matm=0.
call erase_divergence(mat)
call MPI_ALLREDUCE(mat,matm,s1*s2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
mat=matm
return
end subroutine

!*********************************************
!*********************************************
!*********************************************

subroutine adim3s(mat)
implicit none
complex(8) :: mat(:,:,:),matm(size(mat(:,1,1)),size(mat(1,:,1)),size(mat(1,1,:)))
integer    :: s1,s2,s3
if(size2==1.or.no_mpi)then
 return
endif
s1=size(mat(:,1,1));s2=size(mat(1,:,1));s3=size(mat(1,1,:));matm=0.
call erase_divergence(mat)
call MPI_ALLREDUCE(mat,matm,s1*s2*s3,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
mat=matm
return
end subroutine

!*********************************************
!*********************************************
!*********************************************

subroutine adim4s(mat)
implicit none
complex(8) :: mat(:,:,:,:),matm(size(mat(:,1,1,1)),size(mat(1,:,1,1)),size(mat(1,1,:,1)),size(mat(1,1,1,:)))
integer    :: s(4),p
if(size2==1.or.no_mpi)then
 return
endif
s=shape(mat)
p=s(1)*s(2)*s(3)*s(4)
matm=0.
call erase_divergence(mat)
call MPI_ALLREDUCE(mat,matm,p,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
mat=matm
return
end subroutine

!*********************************************
!*********************************************
!*********************************************

subroutine adim5r(mat)
implicit none
real(8) :: mat(:,:,:,:,:), &
     & matm(size(mat(:,1,1,1,1)),size(mat(1,:,1,1,1)),size(mat(1,1,:,1,1)),size(mat(1,1,1,:,1)),size(mat(1,1,1,1,:)))
integer    :: s(5),p
if(size2==1.or.no_mpi)then
 return
endif
s=shape(mat)
p=s(1)*s(2)*s(3)*s(4)*s(5)
matm=0.
call erase_divergence(mat)
call MPI_ALLREDUCE(mat,matm,p,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
mat=matm
return
end subroutine

!*********************************************
!*********************************************
!*********************************************

subroutine adim5s(mat)
implicit none
complex(8) :: mat(:,:,:,:,:), &
     & matm(size(mat(:,1,1,1,1)),size(mat(1,:,1,1,1)),size(mat(1,1,:,1,1)),size(mat(1,1,1,:,1)),size(mat(1,1,1,1,:)))
integer    :: s(5),p
if(size2==1.or.no_mpi)then
 return
endif
s=shape(mat)
p=s(1)*s(2)*s(3)*s(4)*s(5)
matm=0.
call erase_divergence(mat)
call MPI_ALLREDUCE(mat,matm,p,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
mat=matm
return
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine ardim0(vec,vecm2)
implicit none
real(8)  :: vec,vecm,vecm2
integer  :: sizev
 if(size2==1.or.no_mpi)then
   vecm2=vec
   return
 endif
 vecm=0;sizev=1
 call erase_divergence(vec)
 call MPI_ALLREDUCE(vec,vecm,sizev,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 vecm2=vecm/dble(size2)
return
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

subroutine ardim1(vec,vecm2)
implicit none
real(8)  :: vec(:),vecm(size(vec)),vecm2(size(vec))
integer  :: sizev
 if(size2==1.or.no_mpi)then
   vecm2=vec
   return
 endif
 vecm=0;sizev=size(vec)
 call erase_divergence(vec)
 call MPI_ALLREDUCE(vec,vecm,sizev,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 vecm2=vecm/dble(size2)
return
end subroutine

!*********************************************
!*********************************************
!*********************************************

subroutine ardim2(mat,matm)
implicit none
real(8) :: mat(:,:),matm(size(mat(:,1)),size(mat(1,:)))
integer :: s1,s2
if(size2==1.or.no_mpi)then
 matm=mat
 return
endif
 s1=size(mat(:,1));s2=size(mat(1,:));matm=0.
 call erase_divergence(mat)
 call MPI_ALLREDUCE(mat,matm,s1*s2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 matm=matm/dble(size2)
return
end subroutine

!*********************************************
!*********************************************
!*********************************************

subroutine ardim3(mat,matm)
implicit none
real(8)  :: mat(:,:,:),matm(size(mat(:,1,1)),size(mat(1,:,1)),size(mat(1,1,:)))
integer  :: s1,s2,s3
 if(size2==1.or.no_mpi)then
   matm=mat
   return
 endif
 s1=size(mat(:,1,1));s2=size(mat(1,:,1));s3=size(mat(1,1,:));matm=0.
 call erase_divergence(mat)
 call MPI_ALLREDUCE(mat,matm,s1*s2*s3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 matm=matm/dble(size2)
return
end subroutine

!*********************************************
!*********************************************

subroutine ardim1s(vec)
implicit none
real(8)  :: vec(:),vecm(size(vec))
integer :: sizev
if(size2==1.or.no_mpi)then
 return
endif
 vecm=0;sizev=size(vec)
 call erase_divergence(vec)
 call MPI_ALLREDUCE(vec,vecm,sizev,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 vec=vecm
return
end subroutine

!*********************************************
!*********************************************
!*********************************************

subroutine ardim2s(mat)
implicit none
real(8)  :: mat(:,:),matm(size(mat(:,1)),size(mat(1,:)))
integer :: s1,s2
 if(size2==1.or.no_mpi)then
  return
 endif
 s1=size(mat(:,1));s2=size(mat(1,:));matm=0.
 call erase_divergence(mat)
 call MPI_ALLREDUCE(mat,matm,s1*s2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 mat=matm
return
end subroutine

!*********************************************
!*********************************************
!*********************************************

subroutine ardim3s(mat)
implicit none
real(8)  :: mat(:,:,:),matm(size(mat(:,1,1)),size(mat(1,:,1)),size(mat(1,1,:)))
integer :: s1,s2,s3
 if(size2==1.or.no_mpi)then
  return
 endif
 s1=size(mat(:,1,1));s2=size(mat(1,:,1));s3=size(mat(1,1,:));matm=0.
 call erase_divergence(mat)
 call MPI_ALLREDUCE(mat,matm,s1*s2*s3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 mat=matm
return
end subroutine


!*********************************************
!*********************************************
!*********************************************

subroutine ardim4s(mat)
implicit none
real(8)  :: mat(:,:,:,:),matm(size(mat(:,1,1,1)),size(mat(1,:,1,1)),size(mat(1,1,:,1)),size(mat(1,1,1,:)))
integer  :: s(4),p
  if(size2==1.or.no_mpi)then
    return
  endif
  s=shape(mat)
  matm=0. 
  p=s(1)*s(2)*s(3)*s(4)
  call erase_divergence(mat)
  call MPI_ALLREDUCE(mat,matm,p,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  mat=matm
return
end subroutine

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
! 
! subroutine ardim1_mask(vec,vec2,mask,onl,onlyav)
! implicit none
! real(8),intent(inout) :: vec(:),vec2(size(vec))
! real(8)               :: vecm(size(vec))
! integer               :: i,j,k,sizev,count
! integer,intent(in)    :: mask
! integer               :: tot,maskm
! logical,optional      :: onl,onlyav
! 
!  if(size2==1.or.no_mpi) then
!   vec2=0.
!   return
!  endif
!  
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  vecm=0;sizev=size(vec);maskm=0.
!  call erase_divergence(vec)
!  call MPI_ALLREDUCE(mask,maskm,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!  call MPI_ALLREDUCE(vec*dble(mask),vecm,sizev,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!  if(.not.present(onlyav))then
!   vec2=0.
!   call MPI_ALLREDUCE((vec**2)*dble(mask),vec2,sizev,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!  endif
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! 
!  if(maskm>0)then
!    vecm=vecm/dble(maskm)
!    if(.not.present(onl)) vec = vecm
!    if(.not.present(onlyav))then
!      vec2 = vec2/dble(maskm)
!      vec2 = sqrt( abs (vecm**2-vec2)/dble(maskm) )
!    endif
!  endif
! 
! return
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine el_mask(vec,vec2,mask,onl,onlyav)
! implicit none
! real(8),intent(inout) :: vec,vec2
! real(8)               :: vecm
! integer               :: i,j,k,sizev,count
! integer,intent(in)    :: mask
! integer               :: tot,maskm
! logical,optional      :: onl,onlyav
! 
!  if(size2==1.or.no_mpi) then
!    vec2=0.
!    return
!  endif 
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! 
!  vecm=0;maskm=0.
!  call erase_divergence(vec)
!  call MPI_ALLREDUCE(mask,maskm,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!  call MPI_ALLREDUCE(vec*dble(mask),vecm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!  if(.not.present(onlyav))then
!   vec2=0.
!   call MPI_ALLREDUCE((vec**2)*dble(mask),vec2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!  endif
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! 
!  if(maskm>0)then
!    vecm=vecm/dble(maskm)
!    if(.not.present(onl)) vec  = vecm
!    if(.not.present(onlyav))then
!      vec2 = vec2/dble(maskm)
!      vec2 = sqrt( abs (vecm**2-vec2)/dble(maskm) )
!    endif
!  endif
! 
! return
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine ardim2_mask(mat,matm2,mask,onl,onlyav)
! implicit none
!  real(8),intent(inout) :: mat(:,:)
!  real(8),intent(inout) :: matm2(:,:)
!  real(8)               :: matm(size(mat(:,1)),size(mat(1,:)))
!  integer               :: i,j,k,s1,s2
!  integer,intent(in)    :: mask
!  integer               :: maskm
!  logical,optional      :: onl,onlyav
! 
!  if(size2==1.or.no_mpi) then
!   matm2=0.
!   return
!  endif
! 
!  if(ANY(shape(matm2)-shape(mat)>0) ) then
!    write(*,*) 'shape matm2 : ',shape(matm2)
!    write(*,*) 'shape mat   : ',shape(mat)
!    stop 'error ardim2_mask MPI_routines'
!  endif
! 
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  s1=size(mat,1);s2=size(mat,2); matm=0.; maskm=0; 
!  call erase_divergence(mat)
!  call MPI_ALLREDUCE(mask,maskm,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!  call MPI_ALLREDUCE(mat*dble(mask),matm,s1*s2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!  if(.not.present(onlyav))then
!   matm2=0.
!   call MPI_ALLREDUCE(mat**2*dble(mask),matm2,s1*s2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!  endif
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! 
!   if(maskm>0)then
!     matm=matm/dble(maskm)
!     if(.not.present(onl)) mat = matm
!     if(.not.present(onlyav))then
!       matm2 = matm2/dble(maskm)
!       matm2 = sqrt( abs ((matm**2-matm2)/dble(maskm)) )
!     endif
!   endif
! 
! return
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine ardim3_mask(mat,matm2,mask,onl,onlyav)
! implicit none
!  real(8),intent(inout) :: mat(:,:,:)
!  real(8),intent(inout) :: matm2(:,:,:)
!  real(8)               :: matm(size(mat,1),size(mat,2),size(mat,3))
!  integer               :: i,j,k,s(3),p
!  integer,intent(in)    :: mask
!  integer               :: maskm
!  logical,optional      :: onl,onlyav
! 
!  if(size2==1.or.no_mpi) then
!   matm2=0.
!   return
!  endif
!  if(norme(shape(mat)-shape(matm2))>1.d-4) stop 'error ardim3_mask, bad shapes'
!  s=shape(mat)
!  p=s(1)*s(2)*s(3)
! 
!   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!   matm=0.; maskm=0;
!   call erase_divergence(mat)
!   call MPI_ALLREDUCE(mask,maskm,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!   call MPI_ALLREDUCE(mat*dble(mask),matm,p,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!   if(.not.present(onlyav))then
!     matm2=0.
!     call MPI_ALLREDUCE(mat**2*dble(mask),matm2,p,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!   endif
!   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! 
!  if(maskm>0)then
!    matm=matm/dble(maskm)
!    if(.not.present(onl)) mat   =  matm
!    if(.not.present(onlyav))then
!      matm2 = matm2/dble(maskm)
!      matm2 = sqrt( abs (matm**2-matm2)/dble(maskm) )
!    endif
!  endif
! 
! return
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine ardim4_mask(rmat,matm2,mask,onl,onlyav)
! implicit none
!  real(8)               :: rmat(:,:,:,:)
!  real(8)               :: matm2(:,:,:,:)
!  real(8)               :: matm(size(rmat,1),size(rmat,2),size(rmat,3),size(rmat,4))
!  integer               :: i,j,k,s(4),p
!  integer               :: mask
!  integer               :: maskm
!  logical,optional      :: onl,onlyav
! 
!  if(size2==1.or.no_mpi) then
!   matm2=0.
!   return
!  endif
! 
!  if(ANY(shape(rmat)-shape(matm2)>0)) stop 'error ardim4_mask, bad shapes ...'
!  if(ANY(shape(rmat)-shape(matm )>0)) stop 'error ardim4_mask, bad shapes ...'
!  
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! 
!  s=shape(rmat); p=s(1)*s(2)*s(3)*s(4); matm=0.; maskm=0; 
!  
!  call erase_divergence(rmat)
!  call MPI_ALLREDUCE(mask,maskm,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!  call MPI_ALLREDUCE(rmat*dble(mask),matm,p,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
! 
!  if(.not.present(onlyav))then
!    matm2=0.d0
!    call MPI_ALLREDUCE((rmat**2.d0)*dble(mask),matm2,p,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!  endif
! 
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! 
!  if(maskm>0)then
!    matm=matm/dble(maskm)
!    if(.not.present(onl)) rmat = matm
!    if(.not.present(onlyav))then
!      matm2 = matm2/dble(maskm)
!      matm2 = sqrt( abs(matm**2-matm2)/dble(maskm) )
!    endif
!  endif
! 
! return
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! 
! subroutine ardim5_mask(rmat,matm2,mask,onl,onlyav)
! implicit none
!  real(8),intent(inout) :: rmat(:,:,:,:,:)
!  real(8),intent(inout) :: matm2(:,:,:,:,:)
!  real(8)               :: matm(size(rmat(:,1,1,1,1)),size(rmat(1,:,1,1,1)),size(rmat(1,1,:,1,1)), &
!                              & size(rmat(1,1,1,:,1)),size(rmat(1,1,1,1,:)))
!  integer               :: i,j,k,s(5),p
!  integer,intent(in)    :: mask
!  integer               :: maskm
!  logical,optional      :: onl,onlyav
! 
!  if(size2==1.or.no_mpi) then
!   matm2=0.
!   return
!  endif
! 
!  if(norme(shape(rmat)-shape(matm2))>1.d-4) stop 'error ardim5_mask, bad shapes'
!  
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  s=shape(rmat);matm=0.; maskm=0; 
!  p=s(1)*s(2)*s(3)*s(4)*s(5)
!  call erase_divergence(rmat)
!  call MPI_ALLREDUCE(mask,maskm,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
!  call MPI_ALLREDUCE(rmat*dble(mask),matm,p,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!  if(.not.present(onlyav))then
!   matm2=0.
!   call MPI_ALLREDUCE(rmat**2*dble(mask),matm2,p,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!  endif
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! 
!  if(maskm>0)then
!    matm=matm/dble(maskm)
!    if(.not.present(onl))rmat   =  matm
!    if(.not.present(onlyav))then
!      matm2 = matm2/dble(maskm)
!      matm2 = sqrt( abs(matm**2-matm2)/dble(maskm) )
!    endif
!  endif
! 
! return
! end subroutine
! 
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************
! !*********************************************

end module
