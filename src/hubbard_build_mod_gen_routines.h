!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine shoot_next_mu(np,muvec,Nvec,maxdiff,mixing,target_N,reshoot,chem_shift)
 use comms, only: pub_my_node_id
 implicit none
 integer           :: np,nmu_step
 real(8)           :: muvec(np),Nvec(np),chem_shift,maxdiff,mixing,target_N
 real(8)           :: mat(np,np)
 real(16)          :: vectorNmu(np),vectormu(np)
 real(16)          :: f(100),x(100),fp(100),x0,fp0
 real(16)          :: tt1,tt2,tt,vecy(np)
 integer           :: i,j,id
 logical           :: reshoot
 real(16)          :: yn,fpyn,zn,fzn

!receives     : vectorNmu(nmu_step-np+1:nmu_step)
!receives     :  vectormu(nmu_step-np+1:nmu_step)
!root finding : f=(N(mu) - target_N)
!Newton       : x_{n+1} - x_{n} = - f(x_n) / f'(x_n+1)

 nmu_step=np
 vectorNmu=Nvec
 vectorMu=muvec

 if(np==1)then
  if(vectorNmu(nmu_step)>target_N)then
   chem_shift=-0.003
  else
   chem_shift= 0.003
  endif
  return
 endif

 do i=1,np
  f(i)=vectorNmu(nmu_step-(np-i))-target_N
  x(i)= vectorMu(nmu_step-(np-i))
 enddo
 x0=x(np); fp0=f(np)

 if(abs(x(2)-x(1))<1.d-8)then
   chem_shift=0.0 ! you got there...
   return
 endif

 if(np==2)then
  fp(1) = (f(2)-f(1)) / (x(2)-x(1))
 endif

 if(np>2)then
   mat=0.0;vecy=0.0
   do i=1,np
    do j=1,np
      mat(j,i)=x(j)**(np-i)
    enddo
    vecy(i)=f(i)
   enddo
   call invert(np,mat); vecy=matmul(mat,vecy)
 endif

 if(np>2)then
  do i=1,np-1
    fp(i)=derivate(np,vecy,i,x0)
  enddo
  if(fp(1)<-1.d-4)then
    write(*,*) 'WARNING : negative derivate in N(mu) - only use two points'
    write(*,*) '          originally [x] points : ', np
    fp=0.0; fp(1) = ( f(np)-f(np-1) ) / ( x(np)-x(np-1) )
  endif 
endif

   if(pub_my_node_id==0)          write(*,*)  'Position0: ', x0
   if(pub_my_node_id==0)          write(*,*)  'function0: ', fp0
 do i=1,np
   if(pub_my_node_id==0)          write(*,*)  'Position : ', x(i)
   if(pub_my_node_id==0)          write(*,*)  'Function : ', f(i)
   if(i<np.and.pub_my_node_id==0) write(*,*)  'Der      : ', fp(i)
 enddo

 if(np==2)then
    if(abs(fp(1))<50.0)then
      chem_shift = sign(real(1.0,kind=16),fp(1)) * ( - fp0 ) * 0.01
    else
      chem_shift =           mixing * ( - fp0 )  / fp(1)
    endif
 endif

 if(np==3)then
   tt = fp0*fp(2) / (2.d0*fp(1))
   chem_shift = - fp0 / ( fp(1) - tt )
 endif

 if(np>=4)then
   tt2 = (fp0**2) * ( 3.d0*(fp(2)**2)-fp(1)*fp(3) ) / ( 6.d0* (fp(1)**4) )
   tt1 = 1.d0 + fp0*fp(2) / ( 2.d0*(fp(1)**2) ) + tt2
   chem_shift = - fp0 / fp(1) * tt1
 endif

 if(np>2 .and. ( (fp0>0. .and. chem_shift>0.) .or. (fp0<0. .and. chem_shift<0.) )  )then
   write(*,*) 'WARNING going in wrong direction, switch back to first derivative only'
   write(*,*) ' fp0                 : ', fp0
   write(*,*) ' obtained chem_shift : ', chem_shift 
   fp=0.0; fp(1) = (f(np)-f(np-1)) / (x(np)-x(np-1))
   if(abs(fp(1))<50.0)then
      chem_shift = sign(real(1.0,kind=16),fp(1)) * ( - fp0 ) * 0.01
    else
      chem_shift =             mixing * ( - fp0 )  / fp(1)
    endif
    write(*,*) 'corrected chem_shift : ', chem_shift
 endif

if(reshoot.and.np>2)then
  yn= x0 - fp0/fp(1)
  fpyn=derivate(np,vecy,1,yn)
  zn= x0 - 2.d0*fp0/( fp(1) + fpyn )
  fzn =derivate(np,vecy,0,zn)
  chem_shift = zn - (fpyn + fp(1))/(3.d0*fpyn - fp(1)) * fzn / fp(1)
  chem_shift = (chem_shift - x0)
 endif

 if(chem_shift> maxdiff) chem_shift= maxdiff
 if(chem_shift<-maxdiff) chem_shift=-maxdiff

 contains

   !------------------------!

  real(16) FUNCTION derivate(n,vec,order,x0)
  integer :: n,order,iorder,derorder
  real(16) :: x0,vec(n),coef
  integer :: i,j,k
  derivate=0.d0
  do i=1,np
   iorder=np-i
   coef=1.0
   do derorder=1,order
    coef=coef*dble(iorder-derorder+1)
   enddo
   if(iorder-order>=0) derivate=derivate+coef*(x0**(iorder-order))*vec(i)
  enddo
  END FUNCTION

   !------------------------!

  SUBROUTINE invert(n,A)
    integer     ::  n,piv(n),INFO
    real(8)     ::  A(n,n)
    real(8)     ::  WORK(n)
    CALL DGETRF(n,n,A,n,piv,INFO)
    CALL DGETRI(n,A,n,piv,WORK,n,INFO)
    if(INFO/=0)then
      write(*,*) 'ERROR : matrix has no inverse'
    endif
  END SUBROUTINE

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine sync_nodes_by_shared_nfsc(totsum,operation,mpi_rank,mpi_size,longdelay)
 !use IFPORT, only: system
 implicit none
 complex(8)       :: totsum,sumtemp
 integer          :: i,jj,mpi_rank,mpi_size,unit__,ressys
 character*(*)    :: operation
 logical          :: checkfile
 logical,optional :: longdelay

  if(mpi_size/=1)then

      write(*,*) 'ENTERING NFS SHARED RANK : ', mpi_rank
      write(*,*) 'SIZE                     : ', mpi_size

#ifdef DMFT_SAFE_NFS
      ressys=system("sleep 1.5")
#endif

      unit__=free_unit(trim(adjustl(operation))//trim(adjustl(tostring(mpi_rank))),unformatted=.true.)  
      write(unit__) totsum
      close(unit__)

      write(*,*) 'ENTERING LOOP'

      do
#ifdef DMFT_SAFE_NFS
         ressys=system("sleep 0.5")
#else
         ressys=system("sleep 0.1")
#endif
         ressys=system(" echo 'I am done'  > onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//trim(adjustl(tostring(mpi_rank))))

         if(mpi_rank==1) ressys=system("echo `ls -l onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//"*   2>&1 | grep -v 'No such' | wc -l` > onetep_nfs_operation_count_"//trim(adjustl(operation)))

         open(unit=unit__,file='onetep_nfs_operation_count_'//trim(adjustl(operation)))
         jj=0
         read(unit__,*,end=67) jj
         67 continue
         close(unit__)
         write(*,*) 'N PROCESS ARE DONE = ',jj
         if(jj==mpi_size) then
           write(*,*) 'ALL ARE DONE SCALAR'
           exit
         endif
      enddo

      write(*,*) 'EXITING LOOP'
#ifdef DMFT_SAFE_NFS
      ressys=system("sleep 2.5")
#endif
      if(present(longdelay)) ressys=system("sleep 2.5")

      write(*,*) 'SUM , mpi_size = ', mpi_size

      totsum=0.0
      do i=1,mpi_size
         inquire(file=trim(adjustl(operation))//trim(adjustl(tostring(i))),exist=checkfile)
         if(.not.checkfile)then
          write(*,*) 'error NFS share, file erased in the process'
          write(*,*) 'MPI_RANK : ', mpi_rank
          stop
         endif
         open(file=trim(adjustl(operation))//trim(adjustl(tostring(i))),unit=unit__,action='read',position='rewind',form='unformatted')
         read(unit__) sumtemp
         close(unit__)
         totsum=totsum+sumtemp
      enddo

      write(*,*) 'CLEANING'

      if(mpi_rank==1)then
#ifdef DMFT_SAFE_NFS
        ressys=system("sleep 2.5")
#else
        ressys=system("sleep 1.5")
#endif
        if(present(longdelay)) ressys=system("sleep 2.5")
        write(*,*) 'COMMAND LINE : ', " rm "//trim(adjustl(operation))//"*" 
        ressys=system(" rm "//trim(adjustl(operation))//"*"//"  > /dev/null 2>&1")
        ressys=system(" rm onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_*"//"  > /dev/null 2>&1")
        ressys=system(" rm onetep_nfs_operation_count_"//trim(adjustl(operation))//"  > /dev/null 2>&1")
      endif

    endif

    write(*,*) 'LEAVING NFS SHARE'

 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine sync_nodes_by_shared_nfs(totsum,operation,mpi_rank,mpi_size,longdelay)
 !use IFPORT, only: system
 implicit none
 real(8)          :: totsum,sumtemp
 integer          :: i,jj,mpi_rank,mpi_size,unit__,ressys
 character*(*)    :: operation
 logical          :: checkfile
 logical,optional :: longdelay

  if(mpi_size/=1)then

      write(*,*) 'ENTERING NFS SHARED RANK : ', mpi_rank
      write(*,*) 'SIZE                     : ', mpi_size

#ifdef DMFT_SAFE_NFS
      ressys=system("sleep 1.5")
#endif

      unit__=free_unit(trim(adjustl(operation))//trim(adjustl(tostring(mpi_rank))),unformatted=.true.)  
      write(unit__) totsum
      close(unit__)

      write(*,*) 'ENTERING LOOP'

      do
#ifdef DMFT_SAFE_NFS
         ressys=system("sleep 0.5")
#else
         ressys=system("sleep 0.1")
#endif
         ressys=system(" echo 'I am done' > onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//trim(adjustl(tostring(mpi_rank))))

         if(mpi_rank==1) ressys=system("echo `ls -l onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//"*   2>&1  | grep -v 'No such' | wc -l` > onetep_nfs_operation_count_"//trim(adjustl(operation)))

         open(unit=unit__,file='onetep_nfs_operation_count_'//trim(adjustl(operation)))
         jj=0
         read(unit__,*,end=67) jj
         67 continue
         close(unit__)
         write(*,*) 'N PROCESS ARE DONE = ',jj
         if(jj==mpi_size) then
           write(*,*) 'ALL ARE DONE SCALAR'
           exit
         endif
      enddo

      write(*,*) 'EXITING LOOP'
#ifdef DMFT_SAFE_NFS
      ressys=system("sleep 2.5")
#endif
      if(present(longdelay)) ressys=system("sleep 2.5") 

      write(*,*) 'SUM , mpi_size = ', mpi_size

      totsum=0.0
      do i=1,mpi_size
         inquire(file=trim(adjustl(operation))//trim(adjustl(tostring(i))),exist=checkfile)
         if(.not.checkfile)then
          write(*,*) 'error NFS share, file erased in the process'
          write(*,*) 'MPI_RANK : ', mpi_rank
          stop
         endif
         open(file=trim(adjustl(operation))//trim(adjustl(tostring(i))),unit=unit__,action='read',position='rewind',form='unformatted')
         read(unit__) sumtemp
         close(unit__)
         totsum=totsum+sumtemp
      enddo
      
      write(*,*) 'CLEANING'

      if(mpi_rank==1)then
#ifdef DMFT_SAFE_NFS
        ressys=system("sleep 2.5")
#else
        ressys=system("sleep 1.5")
#endif
        if(present(longdelay)) ressys=system("sleep 2.5")
        write(*,*) 'COMMAND LINE : ', " rm "//trim(adjustl(operation))//"*" 
        ressys=system(" rm "//trim(adjustl(operation))//"*"//"  > /dev/null 2>&1")
        ressys=system(" rm onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_*"//"  > /dev/null 2>&1")
        ressys=system(" rm onetep_nfs_operation_count_"//trim(adjustl(operation))//"  > /dev/null 2>&1")
      endif

    endif

    write(*,*) 'LEAVING NFS SHARE'

 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine sync_nodes_by_shared_nfs_mat(totsum,operation,mpi_rank,mpi_size,siz,longdelay,safe)
 !use IFPORT, only: system
 implicit none
 integer           ::  siz
 real(8)           ::  totsum(:,:),sumtemp(size(totsum,1),size(totsum,2))
 integer           ::  i,jj,mpi_rank,mpi_size,unit__,retval,iout,ressys
 character*(*)     ::  operation
 logical,optional  ::  longdelay,safe
 character(2000)   ::  executeop

  if(size(totsum,1)/=siz)then
     write(*,*) 'ERROR IN MATRIX SIZE NFS SUM - SYNC NODES NFS_MAT'
     write(*,*) 'size totsum : ', shape(totsum)
     write(*,*) 'entry size  : ', siz
     write(*,*) 'mpi_rank    : ', mpi_rank
     stop
  endif

  if(mpi_size/=1)then

      write(*,*) 'ENTERING NFS SHARED RANK : ', mpi_rank
      write(*,*) 'SIZE                     : ', mpi_size

#ifdef DMFT_SAFE_NFS
      call sleep(1)
#else
      if(present(safe)) &
&     call sleep(1)
#endif

      unit__=free_unit(trim(adjustl(operation))//trim(adjustl(tostring(mpi_rank))),unformatted=.true.)
      write(unit__) totsum
      close(unit__)

      write(*,*) 'ENTERING LOOP'

      do

         call sleep(1)
         executeop=" echo 'I am done'  > onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//trim(adjustl(tostring(mpi_rank)))
         write(*,*) 'ECHO TO NFS KERNEL : ',trim(adjustl(executeop)) 
         retval=system( trim(adjustl(executeop)) )
         write(*,*) 'RETURN VALUE COMMAND LINE : ', retval
         if(retval/=0)then
            unit__=free_unit("onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//trim(adjustl(tostring(mpi_rank))))
            write(unit__,*) "I am done" 
            close(unit__)
         endif

         if(mpi_rank==1) then
           executeop="echo `ls -l onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//&
           & "*   2>&1  | grep -v 'No such' | wc -l` > onetep_nfs_operation_count_"//trim(adjustl(operation))
           write(*,*) 'COUNT ECHO : ',trim(adjustl(executeop))
           retval=system( trim(adjustl( executeop )) )
           write(*,*) 'RETURN VALUE COMMAND LINE : ', retval
           if(retval/=0)then
             call countfiles("onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_",mpi_size,iout)
             unit__=free_unit("onetep_nfs_operation_count_"//trim(adjustl(operation)))
             write(unit__,*) iout 
             close(unit__) 
           endif
         endif

         open(unit=unit__,file='onetep_nfs_operation_count_'//trim(adjustl(operation)))
         jj=0
         read(unit__,*,end=67) jj
         67 continue
         close(unit__)
         write(*,*) 'N PROCESS ARE DONE = ',jj
         if(jj==mpi_size) then
           write(*,*) 'ALL ARE DONE MATRIX'
           exit
         endif
      enddo

      write(*,*) 'EXITING LOOP'
#ifdef DMFT_SAFE_NFS
         call sleep(2)
#else
      if(present(safe)) then
         call sleep(2)
      endif
#endif
      if(present(longdelay)) call sleep(1)

      totsum=0.0
      do i=1,mpi_size
         open(file=trim(adjustl(operation))//trim(adjustl(tostring(i))),unit=unit__,action='read',position='rewind',form='unformatted')
         read(unit__) sumtemp
         close(unit__)
         totsum=totsum+sumtemp
      enddo

      if(mpi_rank==1)then
       if(present(longdelay)) call sleep(2)

#ifdef DMFT_SAFE_NFS
        call sleep(2)
#else
       if(present(safe)) then
         call sleep(2)
       else
         call sleep(1)
       endif
#endif
       if(retval==0)then
        retval=system(" rm "//trim(adjustl(operation))//"*"//"  > /dev/null 2>&1")
        retval=system(" rm onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//"*"//"  > /dev/null 2>&1")
        retval=system(" rm onetep_nfs_operation_count_"//trim(adjustl(operation))//"  > /dev/null 2>&1")
       else
        do i=1,mpi_size
         call rm( trim(adjustl(operation))//trim(adjustl(tostring(i))))
         call rm( "onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//trim(adjustl(tostring(i))) )
         call rm( "onetep_nfs_operation_count_"//trim(adjustl(operation)) )
        enddo
       endif
      
      endif

    endif

 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------!

subroutine rm(filename1)
implicit none
character(*)                     :: filename1
character(len_trim(filename1)+1) :: fileit
integer                          :: unit__
  open(file=trim(adjustl(filename1)),unit=unit__)
  close(unit__,status='delete')
end subroutine

!-------------------------!

subroutine countfiles(filename1,mpi_size,out)
implicit none
character(*)                     :: filename1
integer                          :: i,mpi_size,out
logical                          :: check

 out=0
 do i=1,mpi_size
  inquire(file=trim(adjustl(filename1))//trim(adjustl(tostring(i))),exist=check)
  if(check) out=out+1 
 enddo

end subroutine

!-------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine sync_nodes_by_shared_nfs_recc(totsum,operation,mpi_rank,mpi_size,siz,siz2,longdelay)
 !use IFPORT, only: system
 implicit none
 integer          :: siz,siz2
 complex(8)       :: totsum(siz,siz2),sumtemp(siz,siz2)
 integer          :: i,jj,mpi_rank,mpi_size,unit__,ressys
 character*(*)    :: operation
 logical,optional :: longdelay

  if(mpi_size/=1)then

      write(*,*) 'ENTERING NFS SHARED RANK : ', mpi_rank
      write(*,*) 'SIZE                     : ', mpi_size

#ifdef DMFT_SAFE_NFS
      ressys=system("sleep 1.5")
#endif

      unit__=free_unit(trim(adjustl(operation))//trim(adjustl(tostring(mpi_rank))),unformatted=.true.)
      write(unit__) totsum
      close(unit__)

      write(*,*) 'ENTERING LOOP'

      do
#ifdef DMFT_SAFE_NFS
         ressys=system("sleep 0.5")
#else
         ressys=system("sleep 0.1")
#endif
         ressys=system(" echo 'I am done'  > onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//trim(adjustl(tostring(mpi_rank))))

         if(mpi_rank==1) ressys=system("echo `ls -l onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//"*   2>&1  | grep -v 'No such'| wc -l` > onetep_nfs_operation_count_"//trim(adjustl(operation)))

         open(unit=unit__,file='onetep_nfs_operation_count_'//trim(adjustl(operation)))
         jj=0
         read(unit__,*,end=67) jj
         67 continue
         close(unit__)
         write(*,*) 'N PROCESS ARE DONE = ',jj
         if(jj==mpi_size) then
           write(*,*) 'ALL ARE DONE MATRIX'
           exit
         endif
      enddo

      write(*,*) 'EXITING LOOP'
#ifdef DMFT_SAFE_NFS
      ressys=system("sleep 2.5")
#endif
      if(present(longdelay)) ressys=system("sleep 1.5")

      totsum=0.0
      do i=1,mpi_size
         open(file=trim(adjustl(operation))//trim(adjustl(tostring(i))),unit=unit__,action='read',position='rewind',form='unformatted')
         read(unit__) sumtemp
         close(unit__)
         totsum=totsum+sumtemp
      enddo

      if(mpi_rank==1)then
#ifdef DMFT_SAFE_NFS
       ressys=system("sleep 2.5")
#else
       if(present(longdelay)) ressys=system("sleep 1.5")
       ressys=system("sleep 1.5")
#endif
       ressys=system(" rm "//trim(adjustl(operation))//"*"//"  > /dev/null 2>&1")
       ressys=system(" rm onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//"*"//"  > /dev/null 2>&1")
       ressys=system(" rm onetep_nfs_operation_count_"//trim(adjustl(operation))//"  > /dev/null 2>&1")
      endif

    endif

 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine sync_nodes_by_shared_nfs_matc(totsum,operation,mpi_rank,mpi_size,siz,longdelay)
 !use IFPORT, only: system
 implicit none
 integer          :: siz
 complex(8)       :: totsum(siz,siz),sumtemp(siz,siz)
 integer          :: i,jj,mpi_rank,mpi_size,unit__,ressys
 character*(*)    :: operation
 logical,optional :: longdelay

  if(mpi_size/=1)then

      write(*,*) 'ENTERING NFS SHARED RANK : ', mpi_rank
      write(*,*) 'SIZE                     : ', mpi_size

#ifdef DMFT_SAFE_NFS
      ressys=system("sleep 1.5")
#endif

      unit__=free_unit(trim(adjustl(operation))//trim(adjustl(tostring(mpi_rank))),unformatted=.true.)
      write(unit__) totsum
      close(unit__)

      write(*,*) 'ENTERING LOOP'

      do
#ifdef DMFT_SAFE_NFS
         ressys=system("sleep 0.5")
#else
         ressys=system("sleep 0.1")
#endif
         ressys=system(" echo 'I am done'  > onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//trim(adjustl(tostring(mpi_rank))))

         if(mpi_rank==1) ressys=system("echo `ls -l onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//"*  2>&1  | grep -v 'No such'| wc -l` > onetep_nfs_operation_count_"//trim(adjustl(operation)))

         open(unit=unit__,file='onetep_nfs_operation_count_'//trim(adjustl(operation)))
         jj=0
         read(unit__,*,end=67) jj
         67 continue
         close(unit__)
         write(*,*) 'N PROCESS ARE DONE = ',jj
         if(jj==mpi_size) then
           write(*,*) 'ALL ARE DONE MATRIX'
           exit
         endif
      enddo

      write(*,*) 'EXITING LOOP'
#ifdef DMFT_SAFE_NFS
      ressys=system("sleep 2.5")
#endif
      if(present(longdelay)) ressys=system("sleep 1.5")

      totsum=0.0
      do i=1,mpi_size
         open(file=trim(adjustl(operation))//trim(adjustl(tostring(i))),unit=unit__,action='read',position='rewind',form='unformatted')
         read(unit__) sumtemp
         close(unit__)
         totsum=totsum+sumtemp
      enddo

      if(mpi_rank==1)then
#ifdef DMFT_SAFE_NFS
       ressys=system("sleep 2.5")
#else
       if(present(longdelay)) ressys=system("sleep 1.5")
       ressys=system("sleep 1.5")
#endif
       ressys=system(" rm "//trim(adjustl(operation))//"*"//"  > /dev/null 2>&1" )
       ressys=system(" rm onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//"*"//"  > /dev/null 2>&1")
       ressys=system(" rm onetep_nfs_operation_count_"//trim(adjustl(operation))//"  > /dev/null 2>&1")
      endif

    endif

 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine sync_nodes_by_shared_nfs_vec(totsum,operation,mpi_rank,mpi_size,siz)
 !use IFPORT, only: system
 implicit none
 integer       :: siz
 real(8)       :: totsum(siz),sumtemp(siz)
 integer       :: i,jj,mpi_rank,mpi_size,unit__,ressys
 character*(*) :: operation

  if(mpi_size/=1)then

      write(*,*) 'ENTERING NFS SHARED RANK : ', mpi_rank
      write(*,*) 'SIZE                     : ', mpi_size

#ifdef DMFT_SAFE_NFS
      ressys=system("sleep 1.5")
#endif

      unit__=free_unit(trim(adjustl(operation))//trim(adjustl(tostring(mpi_rank))),unformatted=.true.)
      write(unit__) totsum
      close(unit__)

      write(*,*) 'ENTERING LOOP'

      do
#ifdef DMFT_SAFE_NFS
         ressys=system("sleep 0.5")
#else
         ressys=system("sleep 0.1")
#endif
         ressys=system(" echo 'I am done' > onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//trim(adjustl(tostring(mpi_rank))))

         if(mpi_rank==1) ressys=system("echo `ls -l onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//"*  2>&1  | grep -v 'No such' | wc -l` > onetep_nfs_operation_count_"//trim(adjustl(operation)))

         open(unit=unit__,file='onetep_nfs_operation_count_'//trim(adjustl(operation)))
         jj=0
         read(unit__,*,end=67) jj
         67 continue
         close(unit__)
         write(*,*) 'N PROCESS ARE DONE = ',jj
         if(jj==mpi_size) then
           write(*,*) 'ALL ARE DONE MATRIX'
           exit
         endif
      enddo

      write(*,*) 'EXITING LOOP'

#ifdef DMFT_SAFE_NFS
      ressys=system("sleep 2.5")
#endif

      totsum=0.0
      do i=1,mpi_size
         open(file=trim(adjustl(operation))//trim(adjustl(tostring(i))),unit=unit__,action='read',position='rewind',form='unformatted')
         read(unit__) sumtemp
         close(unit__)
         totsum=totsum+sumtemp
      enddo

      if(mpi_rank==1)then
#ifdef DMFT_SAFE_NFS
       ressys=system("sleep 2.5")
#else
       ressys=system("sleep 1.5")
#endif
       ressys=system(" rm "//trim(adjustl(operation))//"*"//"  > /dev/null 2>&1")
       ressys=system(" rm onetep_nfs_operation_kernel_"//trim(adjustl(operation))//"_"//"*"//"  > /dev/null 2>&1")
       ressys=system(" rm onetep_nfs_operation_count_"//trim(adjustl(operation))//"  > /dev/null 2>&1")
      endif

    endif

 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 integer(4) function lines_in_file(unit_)
 implicit none
 integer       :: unit_
 rewind(unit_)
 lines_in_file=0
 do 
  read(unit_,*,end=66)
  lines_in_file=lines_in_file+1
 enddo
 66 continue
 rewind(unit_)
 end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine check_not_opened(unit_)
 implicit none
 logical       :: is_it_opened
 integer       :: unit_
   INQUIRE(unit=unit_,OPENED=is_it_opened)
   if(is_it_opened)then
      write(*,*) 'ERROR , unit = ', unit_
      write(*,*) 'was already opened, and not supposed to be'
      stop
   endif
 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine check_files_exists(filename)
 implicit none
 character*(*) :: filename
 logical       :: is_it_there
 INQUIRE(file=trim(adjustl(filename)),EXIST=is_it_there)
 if(.not.is_it_there)then
  write(*,*) 'ERROR, file is not present : ', trim(adjustl(filename))
  stop
 endif
 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 integer(4) function free_unit(filename,unformatted)
 implicit none
 integer       :: unit_,ios
 logical       :: is_it_opened
 character*(*) :: filename
 logical,optional :: unformatted
  unit_=20
  do
   unit_=unit_+1
   INQUIRE(unit=unit_,OPENED=is_it_opened,iostat=ios)
   if(.not.is_it_opened.and.ios==0)exit
  enddo
  free_unit=unit_
  if(.not.present(unformatted))then
  open(unit=free_unit,file=trim(adjustl(filename)))
  else
  open(unit=free_unit,file=trim(adjustl(filename)),form='unformatted')
  endif
 end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine mpisum(mat)
   use comms, only: pub_total_num_nodes
   implicit none
   real(8)  :: mat(:,:,:),matm(size(mat(:,1,1)),size(mat(1,:,1)),size(mat(1,1,:)))
   integer  :: i,j,k,s1,s2,s3,ierr
    if(pub_total_num_nodes==1)then
     return
    endif
    s1=size(mat(:,1,1));s2=size(mat(1,:,1));s3=size(mat(1,1,:));matm=0.
    call MPI_ALLREDUCE(mat,matm,s1*s2*s3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    mat=matm
   return
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function h_atoms_occupancy(hub_proj_basis,hub,hat,spin)
   use function_basis,    only: FUNC_BASIS
   use sparse,            only: sparse_get_block
   use parallel_strategy, only: pub_num_hub_atoms_on_node, pub_hub_atoms_on_node, pub_distr_atom
   implicit none
#define n_on_at hub_proj_basis%num_on_atom(pub_distr_atom(hub%orig(hat)))
   type(FUNC_BASIS), intent(in)    :: hub_proj_basis
   type(HUBBARD_MODEL), intent(in) :: hub
   integer                         :: hat
   REAL(8)                         :: h_atoms_occupancy_tmp(hub_proj_basis%max_on_atom,hub_proj_basis%max_on_atom)
   REAL(8)                         :: h_atoms_occupancy(n_on_at,n_on_at) 
   integer                         :: spin
     call sparse_get_block(h_atoms_occupancy_tmp,hub%occupancy_matrix(spin),pub_distr_atom(hub%orig(hat)),pub_distr_atom(hub%orig(hat)))
     if(n_on_at>hub_proj_basis%max_on_atom)then
       write(*,*) 'ERROR inconsistent with number of orbitals'
       write(*,*) 'hub_proj_basis%max_on_atom : ', hub_proj_basis%max_on_atom
       write(*,*) 'hub_proj_basis%num_on_atom : ', hub_proj_basis%num_on_atom
       stop
     endif
     h_atoms_occupancy(1:n_on_at,1:n_on_at) = h_atoms_occupancy_tmp(1:n_on_at,1:n_on_at)
   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function diag(cc)
   implicit none
   real(8) :: cc(:,:),diag(size(cc,1))
   integer :: i
      do i=1,size(cc,1)
        diag(i)=cc(i,i)
      enddo
   end function

   function compare_r(f,g)
   real(8) :: f,g
   integer :: compare_r
    if(f<g) then
     compare_r=-1
    else
     compare_r=1
    endif
   end function

   function drand1()
   implicit none
   real(8) :: drand1
   real(4) :: r
    call random_number(r)
    drand1=dble(r)
   end function

   integer function qsort_rand( lower, upper )
   implicit none
       integer                       :: lower, upper
       real(4)                       :: r
       r=drand1()
       qsort_rand =  lower + nint(r * (upper-lower))
   end function qsort_rand

   subroutine qsort_adj_array_r(array,order)
   implicit none
      real(8), dimension(:)            :: array
      real(8), dimension(size(array))  :: backup
      integer, dimension(size(array))  :: order
      integer                          :: i,j
      backup=array
      do j=1,size(array)
        array(j)=backup(order(j))
      enddo
   end subroutine

   subroutine qsort_adj_array_vecr(array,order)
   implicit none
      real(8), dimension(:,:)                          :: array
      real(8), dimension(size(array,1),size(array,2))  :: backup
      integer, dimension(size(array,2))                :: order
      integer                                          :: i,j
      backup=array
      do j=1,size(array,2)
        array(:,j)=backup(:,order(j))
      enddo
   end subroutine

   subroutine qsort_swap( order, first, second )
   implicit none
       integer, dimension(:)         :: order
       integer                       :: first, second
       integer                       :: tmp
       tmp           = order(first)
       order(first)  = order(second)
       order(second) = tmp
   end subroutine

   recursive subroutine qsort_sort_r( array, order, left, right )
   implicit none
       real(8),dimension(:)          :: array
       integer, dimension(:)         :: order
       integer                       :: left
       integer                       :: right
       integer                       :: i
       integer                       :: last
       if ( left .ge. right ) return
       call qsort_swap( order, left, qsort_rand(left,right) )
       last = left
       do i = left+1, right
           if ( compare_r(array(order(i)), array(order(left)) ) .lt. 0 ) then
               last = last + 1
               call qsort_swap( order, last, i )
           endif
       enddo
       call qsort_swap( order, left, last )
       call qsort_sort_r( array, order, left, last-1 )
       call qsort_sort_r( array, order, last+1, right )
   end subroutine

   subroutine qsort_array_r(array,order2)
   implicit none
       real(8),dimension(:)                    :: array
       real(8),dimension(size(array))          :: backup
       integer,dimension(size(array))          :: order
       integer,dimension(size(array)),optional :: order2
       integer                                 :: i
       do i=1,size(order)
         order(i)=i
       enddo
       call qsort_sort_r( array, order, 1, size(array) )
       do i=1,size(order)
          backup(i)=array(order(i))
       enddo
       array=backup
       if(present(order2)) order2=order
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure real(8) function DEXPc(rr)
  implicit none
   real(8),parameter  :: MAX_EXP  =  700.d0
   real(8),parameter  :: MIN_EXP  = -700.d0
   real(8),intent(in) :: rr 
   if(rr<MAX_EXP) then
    if(rr<MIN_EXP)then
     DEXPc=0.d0
    else
     DEXPc=EXP(rr)
    endif
   else
     DEXPc=EXP(MAX_EXP)
   endif
  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(8) function norm_vector(xx)
  implicit none
  real(8) :: xx(3)
     norm_vector=sqrt(sum(xx(:)**2))
  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function vecprod(x,y)
 implicit none
 real(8), dimension(3),intent(in) ::  x,y
 real(8), dimension(3)            ::  vecprod
  vecprod(1) = x(2)*y(3) -x(3)*y(2)
  vecprod(2) = x(3)*y(1) -x(1)*y(3)
  vecprod(3) = x(1)*y(2) -x(2)*y(1)
 end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function scalprod(xx,yy)
 implicit none
 integer                         ::  j
 real(8),dimension(:),intent(in) ::  xx,yy
 real(8)                         ::  scalprod
  scalprod=0.
  do j=1,size(xx)
   scalprod=scalprod+xx(j)*yy(j)
  enddo
 end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  character(3) function toString(i)
   integer :: i
    if(i>999)then
     write(*,*) 'ERROR function toString 2' ;stop
    endif
    write(toString,'(i3)') i
   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   INTEGER Function StrInt2(ch)
   character*(*) ch
   integer i,j,ilen
   j=0;
   ilen=LEN_TRIM(ch)
   do i=1, ilen !+1 BUG corrected
     if (i.eq.ilen) then
       j=j+IACHAR(ch(i:i))-48
     else
       j=j+INT((IACHAR(ch(i:i))-48)*10**(ilen-i))
     end if
   end do
   StrInt2=j
   END Function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_transpose(chan,rot_vec_angle,hat)
  implicit none
  integer    :: chan,i,j,jj,hat
  real(8)    :: rotation(chan,chan),rotation_t(chan,chan)
  real(8)    :: rot_vec_angle(3,3),ttt(3,3),tttt(chan,chan)

  rotation   = cmp_rot(           rot_vec_angle,  (chan-1)/2, chan ,hat)
  rotation_t = cmp_rot( transpose(rot_vec_angle), (chan-1)/2, chan ,hat) 

  ttt = matmul( rot_vec_angle , transpose(rot_vec_angle) )

  do i=1,3
   if(abs(ttt(i,i)-1.d0)>1.d-5)then
      write(*,*) 'ERROR, rotation transpose issues : ', i,ttt(i,i)
      stop
   endif
  enddo

 if(maxval(abs(rotation_t-transpose(rotation)))>1.d-5)then
    write(*,*) 'error check_transpose routine'
    write(*,*) 'max vlue : ' , maxval(abs(rotation_t-transpose(rotation)))
    stop
  endif

  tttt = matmul( rotation_t , rotation )

  do i=1,chan
   if(abs(tttt(i,i)-1.d0)>1.d-5)then
      write(*,*) ' ERROR, rotation in cubic harmonic issues : ', i, tttt(i,i)
      write(*,*) ' channels  : ', chan
      write(*,*) ' L         : ', (chan-1)/2
      write(*,*) ' product '
      do j=1,chan
        write(*,'(20f10.4)') (tttt(j,jj),jj=1,chan)
      enddo
      write(*,*) 'sum on lines'
      do j=1,chan
        write(*,*) 'line, sum square           : ', j , sum(rotation(j,:)**2)
        write(*,*) 'line, sum square transpose : ', j , sum(rotation(:,j)**2)
      enddo
      stop
   endif
  enddo

  return
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rotate_to_local_basis_green(occupancy1,occupancy2,green,chan,rot_vec_angle1,rot_vec_angle2,hat,hatb,rotateit)
  implicit none
  integer    :: chan,i,j,hat,hatb
  complex(8) :: green(:,:)
  logical    :: read_in_file,check,rotateit
  real(8)    :: occupancy1(:,:),occupancy2(:,:),diagdens(chan),rotation1(chan,chan),rotation2(chan,chan),ddiag(chan,chan)
  real(8)    :: rotation2_spin(2*chan,2*chan)
  real(8)    :: rot_vec_angle1(:,:),rot_vec_angle2(:,:)
 
  read_in_file=maxval(abs(rot_vec_angle1))>1.d-4.and.maxval(abs(rot_vec_angle2))>1.d-4

  if(.not.read_in_file)then
    call eigenvector_matrix(chan,occupancy1,diagdens,rotation1)
    call eigenvector_matrix(chan,occupancy2,diagdens,rotation2)
  else
    call check_transpose(chan,rot_vec_angle1,hat )
    call check_transpose(chan,rot_vec_angle2,hatb)
    rotation1=cmp_rot(rot_vec_angle1,(chan-1)/2,chan,hat)
    rotation2=cmp_rot(rot_vec_angle2,(chan-1)/2,chan,hatb)
  endif

  if(hat==hatb.and.rotateit)then
    inquire(FILE='mask_loc_rot_atom_pm_'//TRIM(ADJUSTL(toString( hat  ))),EXIST=check)
    if(.not.check)then
     open(unit=98071,FILE='mask_loc_rot_atom_pm_'//TRIM(ADJUSTL(toString( hat  ))),form='unformatted')
     write(98071) shape(rotation2)
     write(98071)       rotation2
     close(98071)
    endif
    inquire(FILE='mask_loc_rot_atom_spin_'//TRIM(ADJUSTL(toString( hat  ))),EXIST=check)
    rotation2_spin                              = 0.d0
    rotation2_spin(     1:  chan,     1:  chan) = rotation2
    rotation2_spin(chan+1:2*chan,chan+1:2*chan) = rotation2
    if(.not.check)then
     open(unit=98071,FILE='mask_loc_rot_atom_spin_'//TRIM(ADJUSTL(toString( hat  ))),form='unformatted')
     write(98071) shape(rotation2_spin)
     write(98071)       rotation2_spin
     close(98071)
    endif
  endif

  green=MATMUL(MATMUL(transpose(rotation1),green),rotation2)

  return
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rotate_back_to_carth_basis_sigma(occupancy1,occupancy2,sigma,chan,rot_vec_angle1,rot_vec_angle2,hat,hatb)
  implicit none
  integer    :: chan,i,hat,hatb
  complex(8) :: sigma(:,:)
  logical    :: read_in_file
  real(8)    :: occupancy1(:,:),occupancy2(:,:),diagdens(chan),rotation1(chan,chan),rotation2(chan,chan),ddiag(chan,chan)
  real(8)    :: rot_vec_angle1(:,:),rot_vec_angle2(:,:)

  read_in_file=maxval(abs(rot_vec_angle1))>1.d-4.and.maxval(abs(rot_vec_angle2))>1.d-4

  if(.not.read_in_file)then
    call eigenvector_matrix(chan,occupancy1,diagdens,rotation1)
    call eigenvector_matrix(chan,occupancy2,diagdens,rotation2)
  else
    rotation1=cmp_rot(rot_vec_angle1,(chan-1)/2,chan,hat)
    rotation2=cmp_rot(rot_vec_angle2,(chan-1)/2,chan,hatb)
  endif

  sigma=MATMUL(MATMUL(rotation1,sigma),transpose(rotation2))

  return
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine diago_occupancy(tt,chan,hat,rot_vec_angle,tmp_dens)
  implicit none
  integer             :: chan,i,hat
  real(8)             :: ttt,entropy,tmp(chan,chan),tt(:,:),tt_(size(tt,1),size(tt,2)),diagdens(chan),rotation(chan,chan)
  logical             :: read_in_file
  real(8),optional    :: rot_vec_angle(3,3),tmp_dens(chan,chan)
  real(8)             :: n,d,l(4)

  read_in_file=.false.

  if(maxval(abs(tt))<1.d-8)then
   return
  endif

  if(present(rot_vec_angle))then
    read_in_file=maxval(abs(rot_vec_angle))>1.d-4
  endif

  if(.not.read_in_file)then
    tt_=(tt+transpose(tt))/2.d0
    call eigenvector_matrix(chan,tt_,diagdens,rotation)
  else
    rotation=cmp_rot(rot_vec_angle,(chan-1)/2,chan,hat)
  endif

  tmp=MATMUL(MATMUL(transpose(rotation),tt),rotation)

  if(present(tmp_dens)) tmp_dens=tmp

  ttt=0.d0
  do i=1,chan
    ttt=ttt+tt(i,i)
  enddo

  do i=1,chan
    diagdens(i)=tmp(i,i)
  enddo

  write(*,*) '================================================='
  call write_array(tmp,' DIAGONAL DENSITIES ')
  write(*,*) 'TOTAL (trace rotated matrix)  = ', sum(diagdens)
  write(*,*) 'CHECK (trace original matrix) = ', ttt
  write(*,*) '================================================='


  tmp=tmp/sum(diagdens)
  entropy=0.d0
  do i=1,chan
    n=2.0*tmp(i,i)
    d=n**2/4.0
    l(1)= 1.0 - n  + d
    l(2)= n/2.0 -d
    l(3)= n/2.0 -d
    l(4)= d
    write(*,'(a,i3,2f15.3)') 'ORBITAL / SUM P / ENT : ', i, sum(l), -sum(l*log(l))
    entropy = entropy + sum(l*log(l))
  enddo
  write(*,*) '================================================='
  write(*,*) 'ENTROPY 4 MANY BODY STATES (U=0)  = ', -entropy
  write(*,*) '================================================='
 
   
  entropy = 1.d0 - sum(diag(matmul(tmp,tmp)))
  write(*,*) '================================================='
  write(*,*) 'NORMALIZATION                   = ', sum(diag(tmp))
  write(*,*) 'LINEAR ENTROPY                  = ', entropy
  write(*,*) 'WITH SPIN FACTOR 2              = ', 2.d0*entropy 
  write(*,*) '================================================='

  
  tmp=tmp*2.d0
  tmp=tmp/sum(diag(tmp))
  entropy = 1.d0 - sum(diag(matmul(tmp,tmp)))
  write(*,*) '================================================='
  write(*,*) 'NORMALIZATION                   = ', sum(diag(tmp))
  write(*,*) 'LINEAR ENTROPY FULL RHO (WITH S)= ', entropy
  write(*,*) '================================================='


  return
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function diago_occupancy_func(tt,chan,rotateit,hat,rot_vec_angle)
  implicit none
  integer             :: chan,i,hat
  real(8)             :: diagdens(chan),ttt,tmp(chan,chan),tt(:,:),tt_(size(tt,1),size(tt,2)),rotation(chan,chan)
  real(8)             :: diago_occupancy_func(chan,chan)
  logical             :: read_in_file,rotateit
  real(8),optional    :: rot_vec_angle(3,3)

  read_in_file=.false.

  if(present(rot_vec_angle))then
    read_in_file=maxval(abs(rot_vec_angle))>1.d-4
  endif

  if(.not.read_in_file)then
    tt_=(tt+transpose(tt))/2.d0
    call eigenvector_matrix(chan,tt_,diagdens,rotation)
  else
    rotation=cmp_rot(rot_vec_angle,(chan-1)/2,chan,hat)
  endif

  if(rotateit)then
    diago_occupancy_func=MATMUL(MATMUL(transpose(rotation),tt),rotation)
  else
    diago_occupancy_func=tt 
  endif

  return
  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_array(tt,title)
  implicit none
  character*(*) :: title
  integer :: chan,i,j,k
  real(8) :: tt(:,:)
  chan=size(tt,1)
  write(*,*) '==========================================='
  write(*,*) title
  do i=1,chan
   write(*,'(200f10.4)') (tt(i,j),j=1,chan)
  enddo
  write(*,*) '==========================================='
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eigenvector_matrix(lsize,mat,vaps,eigenvec)
  use comms, only: pub_my_node_id
  implicit none
  integer                   :: lsize,i
  real(8)                   :: WORK(3*lsize),RWORK(3*lsize)
  real(8)                   :: mat(lsize,lsize)
  real(8)                   :: eigenvec(lsize,lsize)
  real(8)                   :: rrr
  integer                   :: INFO
  real(8)                   :: vaps(lsize),assym
     if(lsize<1) stop 'error eigenvector_matrix, 0 dim'
#ifdef debug
     assym=maxval(abs( mat - transpose((mat)) ))
     if(assym>1.d-3) then
         if(pub_my_node_id==0) then
            write(*,*) 'ERROR eigenvector_matrix_r, matrix non symmetric'
            write(*,*) 'assymetry : ', assym
            write(*,*) 'lsize     : ', lsize
            if(lsize<100)then
             call write_array( mat, ' - MATRIX - ')
            endif
         endif
     endif
#endif
     eigenvec=mat
     call DSYEV('V','U',lsize,eigenvec,lsize,vaps,WORK,3*lsize,INFO)
     if(INFO/=0)then
      write(*,*) 'BAD eigenvalue calculation , info = :', INFO
      write(*,*) 'stop calculations...'
      call write_array( mat, 'MATRIX WAS : ')
      stop
     endif
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function cmp_rot(rot_vec,l,ndim,hat)
    IMPLICIT NONE
    integer                   :: ndim,i,j,aa(7),hat
    COMPLEX(8)                :: cmp_rot(ndim,ndim),cmp_tmp(ndim,ndim)
    real(8),     intent(in)   :: rot_vec(3,3)
    real(8)                   :: Rd(5,5),Rf(7,7),Rp(3,3)
    real(8)                   :: tmp(5,5)
    INTEGER                   :: ii,l,bb(7)
    logical                   :: schemerot,check_ngwfs
    
   !-----------------------------------------------------------------!
   ! transformed G/sigma is obtained by G' = ( cmp_rot^T G cmp_rot ) !
   !-----------------------------------------------------------------!

    schemerot=.false.

    inquire(file='mask_ngwfs',exist=check_ngwfs)
    if(check_ngwfs) then
      open(unit=19812,file='mask_ngwfs')
      do ii=1,hat-1
        read(19812,*,end=111)
      enddo
      read(19812,*) (bb(ii),ii=1,ndim)
      close(19812)
      if(.false.)then
       111 continue
       write(*,*) 'ERROR NOT ENOUGH LINES IN MASK_NGWFS'
       stop
      endif
    else
      bb=(/( ii, ii=1,ndim )/)
    endif

    IF(l==1)THEN
     !===========================!
     ! ONETEP IS WITH : 
     ! -1: py
     !  0: pz
     ! +1: px
     !===========================!
     !==============================!
     ! Rp-language
     ! b1 : y
     ! b2 : z
     ! b3 : x
     !==============================!
     call Rotate_p(Rp, rot_vec)
     aa(1)=1
     aa(2)=2
     aa(3)=3
     do i=1,3
      do j=1,3
        cmp_rot(i,j)=Rp(aa(i),aa(j))
      enddo
     enddo
     if(ndim/=3)then
        write(*,*) 'ERROR, ndim not matching, L=1'
        stop
     endif
    ELSEIF(l==2)THEN
     !===========================!
     ! ONETEP IS WITH :
     ! 1 IS m=-2 : dxy
     ! 2 IS m=-1 : dzy
     ! 3 IS m= 0 : d3z2-r2
     ! 4 IS m= 1 : dxz
     ! 5 IS m= 2 : dx2-y2
     !===========================!
     ! THIS ROUTINE WITH NOTATIONS :
     !     yz
     !     zx
     !     xy
     !   x2-y2
     !   3z2-1
     !===========================!
      if(schemerot)then
        CALL Rotate_d(Rd,rot_vec)
        aa(1)=3
        aa(2)=1
        aa(3)=5
        aa(4)=2
        aa(5)=4
        do i=1,5
         do j=1,5
           cmp_rot(i,j)=Rd(aa(i),aa(j))
         enddo
        enddo
      else
        call rotate_our_notation(Rd,rot_vec)
        cmp_rot=Rd
      endif
      if(ndim/=5)then
        write(*,*) 'ERROR, ndim not matching, L=2'
        stop
      endif
    ELSEIF(l==3)THEN
     call Rotate_f(Rf, rot_vec)
    !onetep i is j in Rf-language
     aa(1)=1
     aa(2)=2
     aa(3)=3
     aa(4)=4
     aa(5)=5
     aa(6)=6
     aa(7)=7
    !==============================!
    ! Onetep
    ! b1 fy(3x^2-y^2)
    ! b2 fxyz         
    ! b3 f y z^2        
    ! b4 f   z^3         
    ! b5 f x z^2        
    ! b6 fz(x2-y2)    
    ! b7 fx(x^2-3y^2)
    !==============================!
    ! Rf-language
    ! b1 fy(3x2-y2) 
    ! b2 xyz
    ! b3 y(4z2 -x2   -y2)
    ! b4 z(2z2 -3x2 -3y2)
    ! b5 x(4z2 -x2 - y2)
    ! b6 z(x2 - y2)
    ! b7 x(x2 - 3y2)
    !==============================!
     do i=1,7
      do j=1,7
        cmp_rot(i,j)=Rf(aa(i),aa(j))
      enddo
     enddo
     if(ndim/=7)then
        write(*,*) 'ERROR, ndim not matching, L=3'
        stop
      endif
    ELSEIF(l>3)THEN
       write(*,*) 'rotation for angular momentum l/=2 not implemented'
       stop
    ENDIF

    cmp_tmp=cmp_rot
    do i=1,ndim
     do j=1,ndim
      cmp_rot(i,j)=cmp_tmp(bb(i),bb(j))
     enddo
    enddo

  return
  END function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE rotate_our_notation(Rd,R)
 IMPLICIT NONE
 real(8), intent(in)  :: R(3,3)
 real(8), intent(out) :: Rd(5,5)
 real(8)              :: s3,r1,r2,r3,r4,r5,r6,r7,r8,r9
 
   s3=sqrt(3.d0)
   
   r1=R(1,1)
   r2=R(1,2)
   r3=R(1,3)
   r4=R(2,1)
   r5=R(2,2)
   r6=R(2,3)
   r7=R(3,1)
   r8=R(3,2)
   r9=R(3,3)
   
   Rd(1,1)=r2*r4+r1*r5
   Rd(1,2)=r3*r5+r2*r6
   Rd(1,3)=s3*r3*r6
   Rd(1,4)=r3*r4+r1*r6
   Rd(1,5)=2.d0*r1*r4+r3*r6
   
   Rd(2,1)=r5*r7+r4*r8
   Rd(2,2)=r6*r8+r5*r9
   Rd(2,3)=s3*r6*r9
   Rd(2,4)=r6*r7+r4*r9
   Rd(2,5)= 2.*r4*r7+r6*r9
   
   Rd(3,1)= s3*r7*r8
   Rd(3,2)= s3*r8*r9
   Rd(3,3)= (3.* r9**2 -1.d0)/2.d0
   Rd(3,4)= s3*r7*r9
   Rd(3,5)= s3*(2.*r7**2 + r9**2 - 1.d0)/2.d0
   
   Rd(4,1)= r2*r7+r1*r8
   Rd(4,2)=  r3*r8+r2*r9
   Rd(4,3)=  s3*r3*r9
   Rd(4,4)= r3*r7+r1*r9
   Rd(4,5)= 2.*r1*r7 +r3*r9
   
   Rd(5,1)=r1*r2-r4*r5
   Rd(5,2)=r2*r3-r5*r6
   Rd(5,3)=s3*(r3**2-r6**2)/2.d0
   Rd(5,4)=r1*r3-r4*r6
   Rd(5,5)=(2.*r1**2 + r3**2 - 2.*r4**2 - r6**2)/2.d0
 
 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Rotate_d(Rd, R)

! Rotation of the d electron cubic harmonic orbitals. Real space 
! rotation matrix is R(3,3). The cubic wave functions are defined as follows:

! b1 = 2yz/r^2 * a
! b2 = 2xz/r^2 * a
! b3 = 2xy/r^2 * a
! b4 = (x^2-y^2)/r^2 * a
! b5 = (3z^2-r^2)/r^2 * a

! Rotation changes them as follows 

! R.b1 =   (R(2,2)*R(3,3) + R(2,3)*R(3,2))*b1 + 
!          (R(2,1)*R(3,3) + R(2,3)*R(3,1))*b2 
!        + (R(2,1)*R(3,2) + R(2,2)*R(3,1))*b3
!        + (R(2,1)*R(3,1) - R(2,2)*R(3,2))*b4 + 
!          sqrt(3.)*R(2,3)*R(3,3)*b5

  IMPLICIT NONE
  real(8), intent(in)  :: R(3,3)
  real(8), intent(out) :: Rd(5,5)

  Rd(1,1) = R(2,2)*R(3,3) + R(2,3)*R(3,2);
  Rd(1,2) = R(2,1)*R(3,3) + R(2,3)*R(3,1);
  Rd(1,3) = R(2,1)*R(3,2) + R(2,2)*R(3,1);
  Rd(1,4) = R(2,1)*R(3,1) - R(2,2)*R(3,2);
  Rd(1,5) = sqrt(3.)*R(2,3)*R(3,3);
  Rd(2,1) = R(1,2)*R(3,3) + R(1,3)*R(3,2);
  Rd(2,2) = R(1,1)*R(3,3) + R(1,3)*R(3,1);
  Rd(2,3) = R(1,1)*R(3,2) + R(1,2)*R(3,1);
  Rd(2,4) = R(1,1)*R(3,1) - R(1,2)*R(3,2);
  Rd(2,5) = sqrt(3.)*R(1,3)*R(3,3);
  Rd(3,1) = R(1,2)*R(2,3) + R(1,3)*R(2,2);
  Rd(3,2) = R(1,1)*R(2,3) + R(1,3)*R(2,1);
  Rd(3,3) = R(1,1)*R(2,2) + R(1,2)*R(2,1);
  Rd(3,4) = R(1,1)*R(2,1) - R(1,2)*R(2,2);
  Rd(3,5) = sqrt(3.)*R(1,3)*R(2,3);
  Rd(4,1) = R(1,2)*R(1,3) - R(2,2)*R(2,3);
  Rd(4,2) = R(1,1)*R(1,3) - R(2,1)*R(2,3);
  Rd(4,3) = R(1,1)*R(1,2) - R(2,1)*R(2,2);
  Rd(4,4) = 0.5*(R(1,1)**2+R(2,2)**2) -0.5*(R(1,2)**2 + R(2,1)**2);
  Rd(4,5) = 0.5*sqrt(3.)*(R(1,3)**2-R(2,3)**2);
  Rd(5,1) = sqrt(3.)*R(3,2)*R(3,3);
  Rd(5,2) = sqrt(3.)*R(3,1)*R(3,3);
  Rd(5,3) = sqrt(3.)*R(3,1)*R(3,2);
  Rd(5,4) = 0.5*sqrt(3.)*(R(3,1)**2-R(3,2)**2);
  Rd(5,5) = -0.5 + 1.5*R(3,3)**2;

return
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Rotate_p(Rp, R)
implicit none
real(8) :: R(3,3),Rp(3,3)
real(8) :: r1,r2,r3,r4,r5,r6,r7,r8,r9

! Rotation of the p electron cubic harmonic orbitals. Real space
! rotation matrix is R(3,3). The cubic wave functions are defined as follows:

!  m=-1, y 
!  m= 0, z 
!  m= 1, x 

r1=R(1,1)
r2=R(1,2)
r3=R(1,3)
r4=R(2,1)
r5=R(2,2)
r6=R(2,3)
r7=R(3,1)
r8=R(3,2)
r9=R(3,3)

Rp(1,1)= r5 
Rp(1,2)= r6
Rp(1,3)= r4
Rp(2,1)= r8
Rp(2,2)= r9
Rp(2,3)= r7
Rp(3,1)= r2
Rp(3,2)= r3
Rp(3,3)= r1

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Rotate_f(Rf, R)
implicit none
real(8) :: R(3,3),Rf(7,7)
real(8) :: a(7),b(7),c(7),d(7),e(7),f(7),g(7)
real(8) :: r1,r2,r3,r4,r5,r6,r7,r8,r9

! Rotation of the f electron cubic harmonic orbitals. Real space
! rotation matrix is R(3,3). The cubic wave functions are defined as follows:

! b1 m=-3, fy(3x^2-y^2) 
! b2 m=-2, fxyz         
! b3 m=-1, fyz^2        
! b4 m= 0, fz^3         
! b5 m= 1, fxz^2        
! b6 m= 2, fz(x2-y2)    
! b7 m= 3, fx(x^2-3y^2) 

r1=R(1,1)
r2=R(1,2)
r3=R(1,3)
r4=R(2,1)
r5=R(2,2)
r6=R(2,3)
r7=R(3,1)
r8=R(3,2)
r9=R(3,3)

a(1) = ( 2.d0*r2*(4.d0*r1*r4+r3*r6)+4.d0*r5*( r1**2 - r4**2 )+ r5 * ( r3**2 - r6**2 ) ) / 4.d0
a(2) = ( 4.d0*(r2*r4*r7 + r1*r5*r7 + r1*r4*r8) + r3*r6*r8 + r3*r5*r9 +r2*r6*r9 )/sqrt(6.d0)
a(3) = sqrt(5.d0/3.d0) * (2.d0*r8*(4.d0*r4*r7 + r6*r9) + r5*(4.d0*r7**2 + r9**2 - 1.d0) ) /4.d0
a(4) = sqrt(5.d0/2.d0) * r8 * (4.d0*r7**2 + r9**2 - 1.d0 ) / 2.d0
a(5) = sqrt(5.d0/3.d0) * ( 2.d0*r8 * (4.d0*r1*r7+r3*r9 ) + r2*(4.d0*r7**2 +r9**2 - 1.d0) ) /4.d0
a(6) = ( 8.d0*r7*( r1*r2 - r4*r5) + 4.d0*r8*(r1**2-r4**2) + r8*(r3**2-r6**2) + 2.d0*r9*(r2*r3 - r5*r6) ) / ( 2.d0*sqrt(6.d0) )   
a(7) = ( r2*(r3**2-4.d0*r4**2 -r6**2 +4.d0*r1**2) - 2.d0*r5*( 4.d0*r1*r4 + r3*r6 ) ) /4.d0

b(1) = sqrt(3.d0/2.d0)* ( r3*(r2*r4 + r1*r5) + r6*(r1*r2 - r4*r5) )
b(2) = (r3*r5 + r2*r6)*r7 + (r3*r4 + r1*r6)*r8 +(r2*r4 + r1*r5)*r9
b(3) = sqrt(5.d0/2.d0) * ( r6*r7*r8 + r5*r7*r9 + r4*r8*r9 )
b(4) = sqrt(15.d0)*r7*r8*r9
b(5) = sqrt(5.d0/2.d0)*(r3*r7*r8 + r2*r7*r9 + r1*r8*r9 )
b(6) = r7*(r2*r3 - r5*r6) + r8*(r1*r3 - r4*r6) + r9*(r1*r2 - r4*r5)
b(7) = sqrt(3.d0/2.d0)*( r3*(r1*r2 - r4*r5) - r6*(r2*r4 + r1*r5) )

c(1) = sqrt(15.d0)*( r3**2*r5 + 2.d0*r2*r3*r6 - r5*r6**2) / 4.d0
c(2) = sqrt(5.d0/2.d0) * ( r3*r6*r8 + r3*r5*r9 + r2*r6*r9 )
c(3) = (10.d0*r6*r8*r9 + r5*(5.d0*r9**2 - 1.d0) ) / 4.d0 
c(4) = sqrt(3.d0/2.d0)*r8*( 5.d0*r9**2 - 1.d0 ) / 2.d0
c(5) = ( 10.d0*r3*r8*r9 + r2*(5.d0*r9**2 - 1.d0) ) / 4.d0
c(6) = sqrt(5.d0/2.d0)* ( r3*(r3*r8+2.d0*r2*r9) - r6*(r6*r8 + 2.d0*r5*r9) ) / 2.d0
c(7) = sqrt(15.d0)*( r2*(r3**2 - r6**2) - 2.d0*r3*r5*r6 ) / 4.d0

d(1) = -sqrt(5.d0/2.d0)*r6*(r6**2 - 3.d0*r3**2) / 2.d0 
d(2) =  sqrt(15.d0) * r3*r6*r9
d(3) =  sqrt(3.d0/2.d0)* r6 * (5.d0*r9**2 - 1.d0) / 2.d0
d(4) =  r9 * (5.d0*r9**2 - 3.d0) / 2.d0
d(5) =  sqrt(3.d0/2.d0) * r3 * (5.d0*r9**2 - 1.d0) / 2.d0 
d(6) =  sqrt(15.d0)*r9*(r3**2 - r6**2) / 2.d0 
d(7) =  sqrt(5.d0/2.d0) * r3 * ( r3**2 - 3.d0*r6**2 ) / 2.d0

e(1) =  sqrt(15.d0)*(r3**2*r4 + 2.d0*r1*r3*r6 - r4*r6**2) / 4.d0
e(2) =  sqrt(5.d0/2.d0) * (r3*r6*r7 + r3*r4*r9 + r1*r6*r9)
e(3) =  (10.d0 * r6*r7*r9 + r4*(5.d0*r9**2 - 1.d0))/4.d0
e(4) =  sqrt(3.d0/2.d0)*r7*( 5.d0*r9**2 - 1.d0 ) / 2.d0  
e(5) =  ( 10.d0*r3*r7*r9 + r1*(5.d0*r9**2 - 1.d0) ) / 4.d0 
e(6) =  sqrt(5.d0/2.d0) * ( r3**2*r7 + 2.d0*r1*r3*r9 - r6*(r6*r7 + 2.d0*r4*r9 ) ) / 2.d0
e(7) =  sqrt(15.d0) * ( r1*(r3**2 - r6**2 ) - 2.d0*r3*r4*r6 ) / 4.d0

f(1) =  sqrt(3.d0/2.d0) * ( 4.d0*r1*r3*r4 + r6*(2.d0*r1**2 - 2.d0*r4**2 - r6**2 + 3.d0*r3**2) ) / 2.d0
f(2) =  2.d0*r7*(r3*r4 + r1*r6) + r9*(2.d0*r1*r4 + 3.d0*r3*r6)
f(3) =  sqrt(5.d0/2.d0) * ( 4.d0*r4*r7*r9 + r6*(2.d0*r7**2 + 3.d0*r9**2 - 1.d0) ) / 2.d0
f(4) =  sqrt(15.d0) * r9 * ( 2.d0*r7**2 + r9**2 - 1.d0 ) / 2.d0
f(5) =  sqrt(5.d0/2.d0) * ( 4.d0*r1*r7*r9 + r3*(2.d0*r7**2 + 3.d0*r9**2 - 1.d0) ) / 2.d0
f(6) =  2.d0*r7*( r1*r3 - r4*r6 ) + 3.d0*r9*( r3**2 - r6**2 )/2.d0 + r9*(r1**2 - r4**2 )
f(7) =  sqrt(3.d0/2.d0) * ( r3*(2.d0*r1**2 + r3**2 - 2.d0*r4**2 - 3.d0*r6**2 ) - 4.d0*r1*r4*r6 ) / 2.d0

g(1) =  ( r4*(4.d0*r5**2 + r6**2 - 4.d0*r2**2 - r3**2 ) - r1*( 8.d0*r2*r5 + 2.d0*r3*r6 ) ) / 4.d0
g(2) = -(r3*(r6*r7 + r4*r9) + 4.d0*r2*( r5*r7 + r4*r8 ) + r1*( r6*r9 + 4.d0*r5*r8 ) ) / sqrt(6.d0)
g(3) = -sqrt(5.d0/3.d0) * ( 2.d0*r7*( 4.d0*r5*r8 + r6*r9 ) + r4*(4.d0*r8**2 + r9**2 - 1.d0) ) / 4.d0
g(4) = -sqrt(5.d0/2.d0) * r7 * ( 4.d0*r8**2 + r9**2 - 1.d0 ) / 2.d0
g(5) = -sqrt(5.d0/3.d0) * (2.d0*r7*(4.d0*r2*r8 + r3*r9 ) + r1*(4.d0*r8**2 + r9**2 - 1.d0) ) /4.d0
g(6) =  ( r7 * ( r6**2 - 4.d0*r2**2 - r3**2 + 4.d0*r5**2 ) - 8.d0*r8*( r1*r2 - r4*r5 ) - 2.d0*r9*( r1*r3 - r4*r6 ) ) / sqrt(24.d0)
g(7) =  ( 2.d0*r4 * (4.d0*r2*r5 + r3*r6 ) + r1*( r6**2 - 4.d0*r2**2 -r3**2 + 4.d0*r5**2 ) ) / 4.d0 

Rf(:,1)=a
Rf(:,2)=b
Rf(:,3)=c
Rf(:,4)=d
Rf(:,5)=e
Rf(:,6)=f
Rf(:,7)=g

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ch_cap ( c )
    implicit none
    character c
    integer itemp
    itemp = ichar ( c )
    if ( 97 <= itemp .and. itemp <= 122 ) then
      c = char ( itemp - 32 )
    end if
    return
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ch_to_digit ( c, digit )
  implicit none
  character c
  integer digit
  if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
    digit = ichar ( c ) - 48
  else if ( c == ' ' ) then
    digit = 0
  else
    digit = -1
  end if
  return
  end subroutine

  function ch_eqi ( c1, c2 )
  implicit none
  character c1
  character c1_cap
  character c2
  character c2_cap
  logical   ch_eqi
  c1_cap = c1
  c2_cap = c2
  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )
  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if
  return
  end function

  real(8) function StringToReal(s)
  character*(*) :: s
  integer       :: ierror,length
    call s_to_r8(s,StringToReal,ierror,length)
  end function

  subroutine s_to_r8 ( s, dval, ierror, length )
  implicit none

  character c
  real ( kind = 8 ) dval
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer length
  integer nchar
  integer ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s

  nchar = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length+1 == nchar ) then
    length = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

  return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
