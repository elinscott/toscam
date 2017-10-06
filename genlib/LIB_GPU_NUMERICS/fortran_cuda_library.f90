module fortran_cuda

 use fortran_cuda_routines
#ifdef MAGMA
 use magma
#endif

 implicit none
 !private

!-----------------------------------------------------------------------------------!
!                                  --------------                                   !
!                                      WARNING                                      !
!                                  --------------                                   !
!  in this module, the conversion from Fortran (array storage in columns) to C      !
!  (array storage in lines) is taken into account, but can be tricky. For           !
!  matrices diagonalization, I use the fact :                                       !
!  (A-1)A=Id, A_T (A-1)_T =1, (A-1)_T=(A_T)-1                                       !
!  so C receives A_T, diagonalize it and get (A_T)-1 inside C, and fortran          !
!  gets ((A_T)-1)_T, which is A_-1                                                  ! 
!  in MATMUL_cuda the call is done to a C routine which takes into account          !
!  the fortran column memory format.                                                !
!  MAGMA uses the Lapack (and hence Fortran) notations                              !
!-----------------------------------------------------------------------------------!


     INTERFACE cuda_array_of_inverse
       MODULE PROCEDURE cuda_array_of_inverse_r,cuda_array_of_inverse_c
     END INTERFACE
    
     INTERFACE cuda_array_of_inverse_collect
       MODULE PROCEDURE cuda_array_of_inverse_collect_c,cuda_array_of_inverse_collect_r
     END INTERFACE
    
     INTERFACE cuda_array_of_inverse_array
       MODULE PROCEDURE cuda_array_of_inverse_array_r,cuda_array_of_inverse_array_c
     END INTERFACE
    
     INTERFACE matmul_cuda
       MODULE PROCEDURE matmul_cuda_r,matmul_cuda_c
     END INTERFACE

     INTERFACE matmulcuda
       MODULE PROCEDURE matmulcuda_r,matmulcuda_c
     END INTERFACE

     INTERFACE matmul_my_cuda
       MODULE PROCEDURE matmul_my_cuda_r,matmul_my_cuda_c 
     END INTERFACE 

public :: cuda_array_of_inverse_r, cuda_array_of_inverse_c, cuda_array_of_inverse_collect_c, cuda_array_of_inverse_array_c, cuda_array_of_inverse_collect_r, cuda_array_of_inverse_array_r,&
& FFT_CUDA_real_4_,diago_cuda_it_c, diago_cuda_it_r, test_cuda, test_cuda_c, test_cuda_MATRIXMUL, matmul_cuda_r, matmul_cuda_c, matmulcuda_r_cublas, matmulcuda_c_cublas, matinv_magma_complex, &
& matinv_magma_double, matinv_sym_complex, matinv_ge_complex, matmulcuda_r,&
& matmulcuda_c, matmulcuda,matmul_cuda,cuda_array_of_inverse_array,cuda_array_of_inverse_collect,cuda_array_of_inverse,matmul_my_cuda,eigenvector_cuda_type_r,eigenvector_cuda_type_c

public :: eigenvector_cuda_r,eigenvector_gen_cuda_r

integer, public :: CUFFT_FORWARD = -1
integer, public :: CUFFT_INVERSE = 1
integer, public :: CUFFT_R2C = Z'2a' ! Real to Complex (interleaved)
integer, public :: CUFFT_C2R = Z'2c' ! Complex (interleaved) to Real
integer, public :: CUFFT_C2C = Z'29' ! Complex to Complex, interleaved
integer, public :: CUFFT_D2Z = Z'6a' ! Double to Double-Complex
integer, public :: CUFFT_Z2D = Z'6c' ! Double-Complex to Double
integer, public :: CUFFT_Z2Z = Z'69' ! Double-Complex to Double-Complex

interface 
subroutine cufftPlan1d_C2C(plan, nx, batch) bind(C,name='cufftPlan1d_C2C_') 
implicit none
integer    :: plan
integer    :: nx, batch
end subroutine
subroutine cufftPlan1d_Z2Z(plan, nx, batch) bind(C,name='cufftPlan1d_Z2Z_')      
implicit none
integer    :: plan
integer    :: nx, batch
end subroutine
subroutine cufftPlan1d_Z2D(plan, nx, batch) bind(C,name='cufftPlan1d_Z2D_')      
implicit none
integer    :: plan
integer    :: nx, batch
end subroutine
subroutine cufftPlan1d_D2Z(plan, nx, batch) bind(C,name='cufftPlan1d_D2Z_')      
implicit none
integer    :: plan
integer    :: nx, batch
end subroutine

subroutine cufftPlan2d_C2C(plan, nx, ny,batch) bind(C,name='cufftPlan2d_C2C_')
implicit none
integer    :: plan
integer    :: nx, ny,batch
end subroutine
subroutine cufftPlan2d_Z2Z(plan, nx, ny,batch) bind(C,name='cufftPlan2d_Z2Z_')
implicit none
integer    :: plan
integer    :: nx, ny,batch
end subroutine
subroutine cufftPlan2d_Z2D(plan, nx, ny,batch) bind(C,name='cufftPlan2d_Z2D_')
implicit none
integer    :: plan
integer    :: nx, ny,batch
end subroutine
subroutine cufftPlan2d_D2Z(plan, nx, ny,batch) bind(C,name='cufftPlan2d_D2Z_')
implicit none
integer    :: plan
integer    :: nx, ny,batch
end subroutine

subroutine cufftPlan3d_C2C(plan, nx, ny,nz, batch) bind(C,name='cufftPlan3d_C2C_')
implicit none
integer    :: plan
integer    :: nx, ny,nz ,batch
end subroutine
subroutine cufftPlan3d_Z2Z(plan, nx, ny,nz, batch) bind(C,name='cufftPlan3d_Z2Z_')
implicit none
integer    :: plan
integer    :: nx, ny,nz ,batch
end subroutine
subroutine cufftPlan3d_Z2D(plan, nx, ny,nz, batch) bind(C,name='cufftPlan3d_Z2D_')
implicit none
integer    :: plan
integer    :: nx, ny,nz ,batch
end subroutine
subroutine cufftPlan3d_D2Z(plan, nx, ny,nz, batch) bind(C,name='cufftPlan3d_D2Z_')
implicit none
integer    :: plan
integer    :: nx, ny,nz ,batch
end subroutine

subroutine cufftDestroy(plan) bind(C,name='cufftDestroy_') 
implicit none
integer    :: plan
end subroutine

subroutine cufftExecC2C(plan, idata, odata, direction,nx,batch) bind(C,name='cufftExecC2C_') 
implicit none
integer    :: direction
integer    :: plan,nx,batch
complex(4) :: idata(nx*batch),odata(nx*batch)
end subroutine

subroutine cufftExecZ2Z(plan, idata, odata, direction,nx,batch) bind(C,name='cufftExecZ2Z_') 
implicit none
integer    :: direction
integer    :: plan,nx,batch
complex(8) :: idata(nx*batch),odata(nx*batch) 
end subroutine 

subroutine cufftExecZ2D(plan, idata, odata, direction,nx,batch) bind(C,name='cufftExecZ2D_')
implicit none
integer    :: direction
integer    :: plan,nx,batch
real(8)    :: idata(2*nx*batch)
real(8)    :: odata(nx*batch)
end subroutine

subroutine cufftExecD2Z(plan, idata, odata, direction,nx,batch) bind(C,name='cufftExecD2Z_')
implicit none
integer    :: direction
integer    :: plan,nx,batch
real(8)    :: idata(nx*batch)
complex(8) :: odata(nx*batch)
end subroutine

end interface 

contains

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

    subroutine eigenvector_cuda_type_r(itype,lsize,mat,overlap,vaps,eigenvec,force_singleprec)
    implicit none
    integer                   :: lsize,itype,i
    real(8)                   :: mat(:,:)
    real(8)                   :: eigenvec(:,:),overlap(:,:)
    real(8)                   :: vaps(:)
    logical                   :: force_singleprec
    real(4),allocatable       :: vapsr(:)
    real(4),allocatable       :: eigenvecr(:,:),matr(:,:)
    real(4),allocatable       :: overlapr(:,:)

       interface
        subroutine magma_dsygvd_type(itype,n,mat,overlap,vaps,vec)
        implicit none
        integer      :: n,itype
        real(8)      :: vec(n,n),vaps(n),mat(n,n),overlap(n,n)
        end subroutine
        subroutine magma_ssygvd_type(itype,n,mat,overlap,vaps,vec)
        implicit none
        integer      :: n,itype
        real(4)      :: vec(n,n),vaps(n),mat(n,n),overlap(n,n)
        end subroutine
       end interface
      if(.not.force_singleprec)then 
       call magma_dsygvd_type(itype,lsize,mat,overlap,vaps,eigenvec)
      else
       i=lsize
       allocate(vapsr(i),matr(i,i),overlapr(i,i),eigenvecr(i,i))
       matr=mat;overlapr=overlap
       call magma_ssygvd_type(itype,lsize,matr,overlapr,vapsr,eigenvecr)
       eigenvec=eigenvecr; vaps=vapsr
       deallocate(vapsr,matr,overlapr,eigenvecr)
      endif
   end subroutine
 
                  !------------------------------------!

   subroutine eigenvector_cuda_type_c(itype,lsize,mat,overlap,vaps,eigenvec,force_singleprec)
    implicit none
    integer                   :: lsize,itype,i
    complex(8)                :: mat(:,:),eigenvec(:,:),overlap(:,:)
    real(8)                   :: vaps(:)
    logical                   :: force_singleprec
    real(4),allocatable       :: vapsr(:)
    complex(4),allocatable    :: eigenvecr(:,:),matr(:,:)
    complex(4),allocatable    :: overlapr(:,:)
       interface
        subroutine magma_zhegvd_type(itype,n,mat,overlap,vaps,vec)
        implicit none
        integer      :: n,itype
        real(8)      :: vaps(n)
        complex(8)   :: vec(n,n),mat(n,n),overlap(n,n)
        end subroutine
        subroutine magma_chegvd_type(itype,n,mat,overlap,vaps,vec)
        implicit none
        integer      :: n,itype
        real(4)      :: vaps(n)
        complex(4)   :: vec(n,n),mat(n,n),overlap(n,n)
        end subroutine
       end interface
      if(.not.force_singleprec)then
       call magma_zhegvd_type(itype,lsize,mat,overlap,vaps,eigenvec)
      else
       i=lsize
       allocate(vapsr(i),matr(i,i),overlapr(i,i),eigenvecr(i,i))
       matr=mat;overlapr=overlap
       call magma_chegvd_type(itype,lsize,matr,overlapr,vapsr,eigenvecr)
       eigenvec=eigenvecr; vaps=vapsr
       deallocate(vapsr,matr,overlapr,eigenvecr)
      endif
    end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

    subroutine eigenvector_cuda_r(lsize,mat,vaps,eigenvec,force_singleprec)
    implicit none
    integer                   :: lsize,i
    real(8)                   :: mat(:,:)
    real(8)                   :: eigenvec(:,:)
    real(8)                   :: vaps(:)
    real(4),allocatable       :: vapsr(:)
    real(4),allocatable       :: eigenvecr(:,:)
    logical                   :: force_singleprec

    interface
     subroutine magma_dsyevd(n,vec,vaps)
     implicit none
     integer      :: n
     real(8)      :: vec(n,n),vaps(n)
     end subroutine
     subroutine magma_ssyevd(n,vec,vaps)
     implicit none
     integer      :: n
     real(4)      :: vec(n,n),vaps(n)
     end subroutine
    end interface
       if(lsize<1) stop 'error eigenvector_matrix, 0 dim'
      if(.not.force_singleprec)then
       eigenvec=mat
       call magma_dsyevd(lsize,eigenvec,vaps)
      else
       allocate(eigenvecr(lsize,lsize),vapsr(lsize))
       eigenvecr=mat
       call magma_ssyevd(lsize,eigenvecr,vapsr)
       vaps=vapsr; eigenvec=eigenvecr
       deallocate(eigenvecr,vapsr)
      endif
    end subroutine

             !---------------------------------------------!

    subroutine eigenvector_gen_cuda_r(lsize,mat,overlap,vaps,eigenvec,force_singleprec)
    implicit none
    integer               :: lsize,i
    real(8)               :: mat(:,:)
    real(8)               :: eigenvec(:,:),overlap(:,:)
    real(8)               :: vaps(:)
    real(4),allocatable   :: vapsr(:)
    real(4),allocatable   :: eigenvecr(:,:)
    real(4),allocatable   :: overlapr(:,:)
    logical               :: force_singleprec

    interface
     subroutine magma_dsygvd(n,vec,over,vaps)
     implicit none
     integer      :: n
     real(8)      :: vec(n,n),vaps(n),over(n,n)
     end subroutine
     subroutine magma_ssygvd(n,vec,over,vaps)
     implicit none
     integer      :: n
     real(4)      :: vec(n,n),vaps(n),over(n,n)
     end subroutine
    end interface
       if(lsize<1) stop 'error eigenvector_gen_matrix, 0 dim'
      if(.not.force_singleprec)then
       eigenvec=mat
       call magma_dsygvd(lsize,eigenvec,overlap,vaps)
      else
       allocate( vapsr(lsize) ) 
       allocate( eigenvecr(lsize,lsize) )
       allocate( overlapr(lsize,lsize) )
       eigenvecr=mat
       overlapr=overlap
       call magma_ssygvd(lsize,eigenvecr,overlapr,vapsr)
       vaps=vapsr;eigenvec=eigenvecr
       deallocate(vapsr,eigenvecr,overlapr)
      endif
    end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine test_cufft
implicit none
integer,parameter :: nn=400
complex(8) :: a(nn)
complex(4) :: b(nn)
integer(4) :: plan

  call select_gpu

  a=1.d0
  plan=1

  call cufftPlan1d_Z2Z(plan, nn, 1)
  write(*,*) 'PLAN : ',plan
  call cufftExecZ2Z(plan, a, a, -1,nn,1)
  call cufftExecZ2Z(plan, a, a,  1,nn,1)
  call cufftDestroy(plan) 
  write(*,'(a,500f10.2)') 'OUTPUT IS (real) : ', real(a)
  write(*,'(a,500f10.2)') 'OUTPUT IS (im)   : ', aimag(a)

  b=1.d0
  plan=1

  call cufftPlan1d_C2C(plan, nn,1)
  write(*,*) 'PLAN : ',plan
  call cufftExecC2C(plan, b, b, -1,nn,1)
  call cufftExecC2C(plan, b, b,  1,nn,1)
  call cufftDestroy(plan)
  write(*,'(a,500f10.2)') 'OUTPUT IS (real) : ', real(b)
  write(*,'(a,500f10.2)') 'OUTPUT IS (im)   : ', aimag(b)

stop
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

 subroutine cuda_array_of_inverse_r(nnn,nfrequ,Eb,totsum,frequ )
 implicit none
 integer :: nnn,nfrequ
 real(8) :: Eb(nnn,nnn),totsum(nnn,nnn),frequ(nfrequ)

 interface 
  subroutine sum_of_inverse_frequ(nnn,nfrequ,Eb,totsum,frequ)
   implicit none
   integer :: nnn,nfrequ
   real(8) :: Eb(nnn*nnn),totsum(nnn*nnn),frequ(nfrequ)  
  end subroutine
 end interface
  
 call sum_of_inverse_frequ(nnn,nfrequ,Eb,totsum,frequ)

 end subroutine

           !------------------------------!
           !------------------------------!
           !------------------------------!
           !------------------------------!

 subroutine cuda_array_of_inverse_c(nnn,nfrequ,Eb,totsum,frequ,firstlast )
 implicit none
 integer    :: nnn,nfrequ,firstlast
 complex(8) :: Eb(nnn,nnn),totsum(nnn,nnn),frequ(nfrequ)

 interface
  subroutine sum_of_inverse_frequ_complex(nnn,nfrequ,Eb,totsum,frequ,firstlast)
   implicit none
   integer    :: nnn,nfrequ,firstlast
   complex(8) :: Eb(nnn*nnn),totsum(nnn*nnn),frequ(nfrequ)
  end subroutine
 end interface

 call sum_of_inverse_frequ_complex(nnn,nfrequ,Eb,totsum,frequ,firstlast)

 end subroutine

           !------------------------------!
           !------------------------------!
           !------------------------------!
           !------------------------------!

 subroutine cuda_array_of_inverse_collect_c(nnn,nfrequ,Eb,cc1,cc2,frequ,firstlast)
 implicit none
 integer    :: nnn,nfrequ,firstlast
 complex(8) :: Eb(nnn,nnn),frequ(nfrequ),cc(nnn,nnn,nfrequ)
 real(8)    :: cc1(nnn,nnn,nfrequ),cc2(nnn,nnn,nfrequ)
 interface
  subroutine sum_of_inverse_frequ_complex_collect(nnn,nfrequ,Eb,cc,frequ,firstlast)
   implicit none
   integer    :: nnn,nfrequ,firstlast
   complex(8) :: Eb(nnn*nnn),cc(nnn*nnn*nfrequ),frequ(nfrequ)
  end subroutine
 end interface
 call sum_of_inverse_frequ_complex_collect(nnn,nfrequ,Eb,cc,frequ,firstlast)
 cc1=real(cc)
 cc2=aimag(cc)
 end subroutine

           !------------------------------!
           !------------------------------!
           !------------------------------!
           !------------------------------!

 subroutine cuda_array_of_inverse_array_c(nnn,nfrequ,cc,firstlast)
 implicit none
 integer    :: nnn,nfrequ,firstlast
 complex(8) :: cc(nnn,nnn,nfrequ)
 interface
  subroutine sum_of_inverse_frequ_complex_array(nnn,nfrequ,cc,firstlast)
   implicit none
   integer    :: nnn,nfrequ,firstlast
   complex(8) :: cc(nnn*nnn*nfrequ)
  end subroutine
 end interface
 call sum_of_inverse_frequ_complex_array(nnn,nfrequ,cc,firstlast)
 end subroutine

           !------------------------------!
           !------------------------------!
           !------------------------------!
           !------------------------------!

 subroutine cuda_array_of_inverse_collect_r(nnn,nfrequ,Eb,cc,frequ,firstlast)
 implicit none
 integer    :: nnn,nfrequ,firstlast
 real(8)    :: Eb(nnn,nnn),frequ(nfrequ),cc(nnn,nnn,nfrequ)
 interface
  subroutine sum_of_inverse_frequ_collect(nnn,nfrequ,Eb,cc,frequ,firstlast)
   implicit none
   integer    :: nnn,nfrequ,firstlast
   real(8)    :: Eb(nnn*nnn),cc(nnn*nnn*nfrequ),frequ(nfrequ)
  end subroutine
 end interface
 call sum_of_inverse_frequ_collect(nnn,nfrequ,Eb,cc,frequ,firstlast)
 end subroutine

           !------------------------------!
           !------------------------------!
           !------------------------------!
           !------------------------------!

 subroutine cuda_array_of_inverse_array_r(nnn,nfrequ,cc,firstlast)
 implicit none
 integer    :: nnn,nfrequ,firstlast
 real(8)    :: cc(nnn,nnn,nfrequ)
 interface
  subroutine sum_of_inverse_frequ_array(nnn,nfrequ,cc,firstlast)
   implicit none
   integer    :: nnn,nfrequ,firstlast
   real(8)    :: cc(nnn*nnn*nfrequ)
  end subroutine
 end interface
 call sum_of_inverse_frequ_array(nnn,nfrequ,cc,firstlast)
 end subroutine

           !------------------------------!
           !------------------------------!
           !------------------------------!
           !------------------------------!


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

      subroutine FFT_CUDA_real_4_(NX,BATCH,data)
      implicit none
    !----------------------------------------! 
    ! TYPICAL VALUES ARE : NX=256,BATCH=10   !
    !----------------------------------------!
      integer    :: CUFFT_FORWARD, CUFFT_INVERSE
      parameter(CUFFT_FORWARD=-1, CUFFT_INVERSE=1)
      integer    :: CUFFT_R2C, CUFFT_C2R, CUFFT_C2C
      parameter(CUFFT_R2C=X"2a", CUFFT_C2R=X"2c", CUFFT_C2C=X"29")
      integer    :: cudaMemcpyHostToDevice, cudaMemcpyDeviceToHost
      parameter(cudaMemcpyHostToDevice=1, cudaMemcpyDeviceToHost=2)
      integer    :: NX,BATCH
      real       :: PI
      parameter (PI=3.14159)
      integer    :: plan
      integer    :: i,err
      complex(4) :: data(:)

      interface
       subroutine cufftplan1d(plan,nx,cucu,batch)
        implicit none
        integer :: batch,nx,plan,cucu
       end subroutine
       subroutine cufftexecc2c(plan,data,data2,cucu,nx,batch)
        implicit none
        integer    :: batch,nx,plan,cucu
        complex(4) :: data(nx*batch),data2(nx*batch)
       end subroutine
       subroutine cufftdestroy(plan)
        implicit none
        integer :: plan
       end subroutine
      end interface

       if(size(data,1)/=NX*BATCH) stop 'FFT_CUDA_real_4_ size of data is different from NX*BATCH, critical stop'

 !     Create a 1D FFT plan. 
       call cufftplan1d(plan, NX, CUFFT_C2C, BATCH)
 !     Use the CUFFT plan to transform the signal in place.
       call cufftexecc2c(plan, data, data, CUFFT_FORWARD,NX,BATCH)
 !     Inverse transform the signal in place.
       call cufftexecc2c(plan, data, data, CUFFT_INVERSE,NX,BATCH)
 !     Destroy the CUFFT plan.
       call cufftdestroy(plan)

      return
      end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine diago_cuda_it_c(lin,mat)
implicit none
 integer                :: lin,block
 complex(8)             :: mat(:,:),inv(lin,lin)
 complex(8),allocatable :: mat_(:,:),inv_(:,:)
 integer                :: i,j,k,l,lin_

 interface
  subroutine cuda_complex_invert(k,aa,inva,nlines)
    complex(8) :: aa(nlines,nlines),inva(nlines,nlines)
    integer    :: k,nlines
  end subroutine
 end interface

  block=16
  if(allocated(mat_)) deallocate(mat_); if(allocated(inv_)) deallocate(inv_)
                lin_= lin-mod(lin,block)
  if(lin_/=lin) lin_= (lin_/block+1)*block

  allocate(mat_(lin_,lin_),inv_(lin_,lin_))

  mat_(1:lin,1:lin)=mat(1:lin,1:lin)
  do i=1,lin_
   do j=lin+1,lin_
    mat_(i,j) = 0.
   enddo
  enddo
  do j=1,lin_
   do i=lin+1,lin_
    mat_(i,j) = 0.
   enddo
  enddo
  do i=lin+1,lin_
   do j=lin+1,lin_
    mat_(i,j) = 0.
   enddo
  enddo
  do i=lin+1,lin_
   mat_(i,i) = 1.d0
  enddo

  inv_=0.d0
  write(*,*) 'inverting on GPU with linear size : ', lin_
  call cuda_complex_invert(block,mat_,inv_,lin_)
  mat(1:lin,1:lin)=inv_(1:lin,1:lin)
  if(allocated(mat_)) deallocate(mat_); if(allocated(inv_)) deallocate(inv_)

return
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
!**************************************************************************
!**************************************************************************

subroutine diago_cuda_it_r(lin,mat)
implicit none
 integer             :: lin,block
 real(8)             :: mat(:,:),inv(lin,lin)
 real(8),allocatable :: mat_(:,:),inv_(:,:)
 integer             :: i,j,k,l,lin_

 interface
  subroutine cuda_invert_(k,aa,inva,nlines)
   real(8) :: aa(nlines,nlines),inva(nlines,nlines)
   integer :: k,nlines
  end subroutine
 end interface

  block=16
  if(allocated(mat_)) deallocate(mat_); if(allocated(inv_)) deallocate(inv_)
                lin_= lin-mod(lin,block)
  if(lin_/=lin) lin_= (lin_/block+1)*block

  allocate(mat_(lin_,lin_),inv_(lin_,lin_))

  mat_(1:lin,1:lin)=mat(1:lin,1:lin)
  do i=1,lin_
   do j=lin+1,lin_
    mat_(i,j) = 0.
   enddo
  enddo
  do j=1,lin_
   do i=lin+1,lin_
    mat_(i,j) = 0.
   enddo
  enddo
  do i=lin+1,lin_
   do j=lin+1,lin_
    mat_(i,j) = 0.
   enddo
  enddo
  do i=lin+1,lin_
   mat_(i,i) = 1.d0
  enddo

  inv_=0.d0
  call cuda_invert_(block,mat_,inv_,lin_)

  mat(1:lin,1:lin)=inv_(1:lin,1:lin)
  if(allocated(mat_)) deallocate(mat_); if(allocated(inv_)) deallocate(inv_)


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
!**************************************************************************
!**************************************************************************

 subroutine test_cuda
 implicit none

 integer,parameter :: lin=31
 real(8)           :: mat(lin,lin),inv(lin,lin)
 real(4)           :: dd
 integer           :: i,j,k,l

 do i=1,lin; do j=1,lin;  call random_number(dd); mat(i,j)=dd ; enddo; enddo ; 

 if(lin<20)then
 write(*,*) '============================'
 write(*,*) 'matrix'
 do i=1,lin
  do j=1,lin
    write(*,*) i,j,mat(i,j)
  enddo
 enddo
 write(*,*) '============================'
 else
  write(*,*) 'maxval matrix input : ', maxval(abs(mat))
 endif

 inv=mat
 call diago_cuda_it_r(lin,inv)

 if(lin<20)then
 write(*,*) '============================'
 write(*,*) 'inverse matrix'
 do i=1,lin
  do j=1,lin
    write(*,*) i,j,inv(i,j)
  enddo
 enddo
 write(*,*) '============================'
 else
  write(*,*) 'maxval inv output : ', maxval(abs(inv))
 endif

 write(*,*) 'check...........'
 write(*,*) maxval(abs(MATMUL(mat,inv))),minval(abs(MATMUL(mat,inv)))
 write(*,*) 'done............'

 stop
 end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

 subroutine test_cuda_c
 implicit none

 integer,parameter :: lin=4
 complex(8)        :: mat(lin,lin),inv(lin,lin)
 real(4)           :: dd,dd1,dd2
 integer           :: i,j,k,l

 do i=1,lin; do j=1,lin; call random_number(dd1);call random_number(dd2);  mat(i,j)=cmplx(dd1,dd2,8) ; enddo; enddo ;

 write(*,*) '============================'
 write(*,*) 'matrix'
 do i=1,lin
  do j=1,lin
    write(*,*) i,j,mat(i,j)
  enddo
 enddo
 write(*,*) '============================'

 inv=mat
 call diago_cuda_it_c(lin,mat)

 write(*,*) '============================'
 write(*,*) 'inverse matrix'
 do i=1,lin
  do j=1,lin
    write(*,*) i,j,mat(i,j)
  enddo
 enddo
 write(*,*) '============================'
                                                   
 write(*,*) 'check...........'
 write(*,*) maxval(abs(MATMUL(mat,inv))),minval(abs(MATMUL(mat,inv)))
 write(*,*) 'done............'

 stop
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

subroutine test_cuda_MATRIXMUL
IMPLICIT NONE
INTEGER,PARAMETER   :: dima=44,dimb=40,dimc=60
real(8),allocatable :: M(:,:),N(:,:),P(:,:),Q(:,:),Q2(:,:),Q3(:,:),Q4(:,:),MT(:,:),NT(:,:),Q2T(:,:)
REAL(4)             :: dd
INTEGER             :: I,J,K,clock,clock_ref,blocksiz
REAL(8)             :: one_over_clock_rate
INTEGER             :: clock_rate

allocate( M(dima,dimb),N(dimb,dimc),P(dima,dimc),Q(dima,dimc),Q2(dima,dimc),Q3(dima,dimc),Q4(dima,dimc) )

M=0.
DO i=1,dima
do j=1,dimb
  call random_number(dd)
  M(i,j)=dd
enddo
ENDDO
N=0.
DO i=1,dimb
do j=1,dimc
  call random_number(dd)
  N(i,j)=dd
enddo
ENDDO

Q =0.
Q2=0.
Q3=0.
Q4=0.

 call freeMem()

 CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate)
 one_over_clock_rate  = 1.d0 / DBLE(clock_rate)

 CALL SYSTEM_CLOCK(COUNT=clock_ref)
 P=MATMUL(M,N)
 CALL SYSTEM_CLOCK(COUNT=clock)
 write(*,*) 'CPU MATMUL CEDRIC LIBRARY TOOK : ',  DBLE(clock-clock_ref) * one_over_clock_rate

 CALL SYSTEM_CLOCK(COUNT=clock_ref)
 write(*,*) 'start GPU 1'
 call matmulcuda_r(M,N,Q,dima,dimb,dimc)
 write(*,*) 'done'
 CALL SYSTEM_CLOCK(COUNT=clock)
 write(*,*) 'GPU MATMUL CEDRIC LIBRARY TOOK : ',  DBLE(clock-clock_ref) * one_over_clock_rate

 CALL SYSTEM_CLOCK(COUNT=clock_ref)
 write(*,*) 'start GPU 2'
 blocksiz=4
 if(mod(dima,blocksiz)/=0)then; write(*,*) 'blocksiz not commensurable dima';stop;endif
 if(mod(dimb,blocksiz)/=0)then; write(*,*) 'blocksiz not commensurable dimb';stop;endif
 if(mod(dimc,blocksiz)/=0)then; write(*,*) 'blocksiz not commensurable dimc';stop;endif
 CALL matmul_cuda_r(blocksiz,M,N,Q2,dima,dimb,dimc)

 write(*,*) 'done'
 CALL SYSTEM_CLOCK(COUNT=clock)
 write(*,*) 'GPU MATMUL CEDRIC LIBRARY TOOK : ',  DBLE(clock-clock_ref) * one_over_clock_rate

 CALL SYSTEM_CLOCK(COUNT=clock_ref)
 write(*,*) 'start GPU 3, single prec'
 call matmulcuda_r_singleprec_bypass(M,N,Q3,dima,dimb,dimc)
 write(*,*) 'done'
 CALL SYSTEM_CLOCK(COUNT=clock)
 write(*,*) 'GPU MATMUL CEDRIC LIBRARY TOOK : ',  DBLE(clock-clock_ref) * one_over_clock_rate

 CALL SYSTEM_CLOCK(COUNT=clock_ref)
 write(*,*) 'start GPU 4, cublas'
 call matmulcuda_r_cublas(M,N,Q4,dima,dimb,dimc)
 write(*,*) 'done'
 CALL SYSTEM_CLOCK(COUNT=clock)
 write(*,*) 'GPU MATMUL CEDRIC LIBRARY TOOK : ',  DBLE(clock-clock_ref) * one_over_clock_rate

 write(*,*) '============================='
 write(*,*) 'matmulcuda_r      Q     : ', maxval(abs(Q))
 write(*,*) 'matmul_cuda_r     Q2    : ', maxval(abs(Q2))
 write(*,*) 'magma single prec Q3    : ', maxval(abs(Q3))
 write(*,*) 'error matmulcuda_r      : ', maxval(abs(P-Q))
 write(*,*) 'error matmul_cuda_r     : ', maxval(abs(P-Q2))
 write(*,*) 'error magma single prec : ', maxval(abs(P-Q3))
 write(*,*) 'error cublas            : ', maxval(abs(P-Q4))
 write(*,*) '============================='

 deallocate(M,N,P,Q,Q2,Q3)
 
STOP
END subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine matmul_cuda_r(block,M,N,Q,dima,dimb,dimc)
implicit none
INTEGER           :: DIMA,DIMB,DIMC
REAL(8)           :: M(DIMA*DIMB),N(DIMB*DIMC),Q(DIMA*DIMC)
INTEGER           :: I,J,K,block
interface
  subroutine matmul_gpu_fortran(block,M,N,Q,a,b,c)
  implicit none
    integer :: a,b,c,block
    real(8) :: M(a*b),N(b*c),Q(a*c)
   end subroutine
  end interface
  CALL matmul_gpu_fortran(block,M,N,Q,DIMA,DIMB,DIMC)
end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine matmul_cuda_c(block,M,N,Q,dima,dimb,dimc)
implicit none
INTEGER           :: DIMA,DIMB,DIMC
complex(8)        :: M(DIMA*DIMB),N(DIMB*DIMC),Q(DIMA*DIMC)
INTEGER           :: I,J,K,block
interface
  subroutine matmul_gpu_fortran_c(block,M,N,Q,a,b,c)
  implicit none
    integer    :: a,b,c,block
    complex(8) :: M(a*b),N(b*c),Q(a*c)
   end subroutine
  end interface
  CALL matmul_gpu_fortran_c(block,M,N,Q,DIMA,DIMB,DIMC)
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

subroutine matmulcuda_r(M,N,Q,dima,dimb,dimc)
implicit none
INTEGER           :: DIMA,DIMB,DIMC
REAL(8)           :: M(DIMA*DIMB),N(DIMB*DIMC),Q(DIMA*DIMC)
interface
  subroutine magma_matrix_multiply(M,N,Q,a,b,c)
  implicit none
    integer :: a,b,c
    real(8) :: M(a*b),N(b*c),Q(a*c)
   end subroutine
  end interface
  CALL magma_matrix_multiply(M,N,Q,DIMA,DIMB,DIMC)
end subroutine

             !-----------------------!

subroutine matmulcuda_c(M,N,Q,dima,dimb,dimc)
implicit none
INTEGER           :: DIMA,DIMB,DIMC
COMPLEX(8)           :: M(DIMA*DIMB),N(DIMB*DIMC),Q(DIMA*DIMC)
interface
  subroutine magma_complex_matrix_multiply(M,N,Q,a,b,c)
  implicit none
    integer :: a,b,c
    complex(8) :: M(a*b),N(b*c),Q(a*c)
   end subroutine
  end interface
  CALL magma_complex_matrix_multiply(M,N,Q,DIMA,DIMB,DIMC)
end subroutine

             !-----------------------!

subroutine matmulcuda_r_cublas(M,N,Q,dima,dimb,dimc)
implicit none
INTEGER           :: DIMA,DIMB,DIMC
REAL(8)           :: M(DIMA*DIMB),N(DIMB*DIMC),Q(DIMA*DIMC)
interface
  subroutine cublas_matrix_multiply(M,N,Q,a,b,c)
  implicit none
    integer :: a,b,c
    real(8) :: M(a*b),N(b*c),Q(a*c)
   end subroutine
  end interface
  CALL cublas_matrix_multiply(M,N,Q,DIMA,DIMB,DIMC)
end subroutine

             !-----------------------!

subroutine matmulcuda_c_cublas(M,N,Q,dima,dimb,dimc)
implicit none
INTEGER           :: DIMA,DIMB,DIMC
COMPLEX(8)           :: M(DIMA*DIMB),N(DIMB*DIMC),Q(DIMA*DIMC)
interface
  subroutine cublas_complex_matrix_multiply(M,N,Q,a,b,c)
  implicit none
    integer :: a,b,c
    complex(8) :: M(a*b),N(b*c),Q(a*c)
   end subroutine
  end interface
  CALL cublas_complex_matrix_multiply(M,N,Q,DIMA,DIMB,DIMC)
end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine matinv_magma_complex_cxx(n,mat)
implicit none
integer              :: n
complex(8)           :: mat(:,:)
integer              :: i,j,k
interface
 subroutine invert_magma_ge(nn,a)
 implicit none
   integer    :: nn
   complex(8) :: a(nn,nn)
 end subroutine
end interface
call invert_magma_ge(n,mat)
end subroutine

subroutine matinv_magma_double_cxx(n,mat)
implicit none
integer              :: n
real(8)              :: mat(:,:)
integer              :: i,j,k
interface
 subroutine invert_magma_double_ge(nn,a)
 implicit none
   integer    :: nn
   real(8) :: a(nn,nn)
 end subroutine
end interface
call invert_magma_double_ge(n,mat)
end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine matinv_sym_complex(n,mat)
implicit none
integer              :: n
complex(8)           :: mat(:,:)
integer              :: i,j,k
real(8)              :: A(n,n),Am(n,n),r0(n,n),C(n,n),y11(n,n),y11m(n,n),y01(n,n)
logical,parameter    :: use_cpu_matmul=.false. 
interface
 subroutine invert_spd(nn,a,inva)
 implicit none
   integer :: nn
   real(8) :: a(nn,nn),inva(nn,nn)
 end subroutine
end interface

! (A+iC)^-1 = y11^-1 - i * y01
!        r0  =  A^-1 * C 
!        y11 = (C*ro+A)
!        y01 = r0*(y11^-1)

  A=real(mat)
  call invert_spd(n,A,Am)
  C=aimag(mat)
  if(use_cpu_matmul)then
    r0  = MATMUL(Am,C)
    y11 = MATMUL(C,r0)+A
    call invert_spd(n,y11,y11m)
    y01 = MATMUL(r0,y11m)
  else
    call matmulcuda_r_cublas(Am,C,r0,n,n,n)
    call matmulcuda_r_cublas(C,r0,y11,n,n,n)  
    y11=y11+A
    call invert_spd(n,y11,y11m)
    call matmulcuda_r_cublas(r0,y11m,y01,n,n,n)
  endif
 do i=1,n
  do k=1,n
    mat(i,k)= CMPLX(y11m(i,k),-y01(i,k),KIND=8)
  enddo
 enddo

return
end subroutine

     !----------------------!

subroutine matinv_ge_complex(n,mat)
implicit none
integer              :: n
complex(8)           :: mat(:,:)
integer              :: i,j,k
real(8)              :: A(n,n),Am(n,n),r0(n,n),C(n,n),y11(n,n),y11m(n,n),y01(n,n)
interface
 subroutine invert_ge(nn,a,inva)
 implicit none
   integer :: nn
   real(8) :: a(nn,nn),inva(nn,nn)
 end subroutine
end interface

! (A+iC)^-1 = y11^-1 - i * y01
!        r0  =  A^-1 * C
!        y11 = (C*ro+A)
!        y01 = r0*(y11^-1)

  A=real(mat)
  Am=A; call magma_fortran_double_(n,Am)
  C=aimag(mat)
  call matmulcuda_r(Am,C,r0,n,n,n)
  call matmulcuda_r(C,r0,y11,n,n,n)
  y11=y11+A
  y11m=y11; call magma_fortran_double_(n,y11m)
  call matmulcuda_r(r0,y11m,y01,n,n,n)
  do i=1,n
   do k=1,n
    mat(i,k)=CMPLX(y11m(i,k),-y01(i,k),KIND=8)
   enddo
  enddo

return
end subroutine

     !----------------------!

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
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine matmul_my_cuda_r(M,N,Q,dima,dimb,dimc)
implicit none
INTEGER             :: DIMA,DIMB,DIMC,block
real(8)             :: M(DIMA,DIMB),N(DIMB,DIMC),Q(DIMA,DIMC)
real(8),allocatable :: M_(:,:),N_(:,:),Q_(:,:)
INTEGER             :: I,J,K,i_,j_
logical,parameter   :: verbose=.true.

  if(allocated(M_)) deallocate(M_); if(allocated(N_)) deallocate(N_);
  if(allocated(Q_)) deallocate(Q_)
  i = size(M,1); i_= i-mod(i,16); if(i_/=i) i_= (i_/16+1)*16
  j = size(M,2); j_= j-mod(j,16); if(j_/=j) j_= (j_/16+1)*16
  if(verbose) write(*,*) 'M allocate array ,i ,j : ', i ,j
  if(verbose) write(*,*) 'M allocate array ,i_,j_: ', i_,j_
  allocate(M_(i_,j_));
  if(verbose) write(*,*) 'allocated,padd with zeros'

  M_=0.d0
  M_(1:i,1:j)=M(1:i,1:j)

  if(verbose) write(*,*) 'reshape matrix N'
  i = size(N,1); i_= i-mod(i,16); if(i_/=i)i_= (i_/16+1)*16
  j = size(N,2); j_= j-mod(j,16); if(j_/=j)j_= (j_/16+1)*16
  if(verbose) write(*,*) 'N allocate array ,i ,j : ', i ,j
  if(verbose) write(*,*) 'N allocate array ,i_,j_: ', i_,j_
  allocate(N_(i_,j_));

  N_=0.d0
  N_(1:i,1:j)=N(1:i,1:j)

  i=size(M_,1);j=size(M_,2);k=size(N_,2)
  if(allocated(Q_))deallocate(Q_); allocate(Q_(i,k))

  if(verbose)then
    write(*,*) 'CALLING MATMUL CUDA, dims : ',i,j,k
    write(*,*) 'shape M ',shape(M_)
    write(*,*) 'shape N ',shape(N_)
    write(*,*) 'shape Q ',shape(Q_)
    write(*,*) 'ijk     ',i,j,k
  endif

  Q_=0.d0
  CALL matmul_cuda_r(16,M_,N_,Q_,i,j,k)
  if(verbose) write(*,*) '...done...'
  Q(1:dima,1:dimc)=Q_(1:dima,1:dimc)

  if(allocated(M_)) deallocate(M_); if(allocated(N_)) deallocate(N_);
  if(allocated(Q_)) deallocate(Q_)

return
end subroutine

  !-------------------------------------------------------------------------------!

subroutine matmul_my_cuda_c(M,N,Q,dima,dimb,dimc)
implicit none
INTEGER                :: DIMA,DIMB,DIMC,block
complex(8)             :: M(DIMA,DIMB),N(DIMB,DIMC),Q(DIMA,DIMC)
complex(8),allocatable :: M_(:,:),N_(:,:),Q_(:,:)
INTEGER                :: I,J,K,i_,j_
logical,parameter      :: verbose=.true.

  if(allocated(M_)) deallocate(M_); if(allocated(N_)) deallocate(N_);
  if(allocated(Q_)) deallocate(Q_)
  i = size(M,1); i_= i-mod(i,16); if(i_/=i)i_= (i_/16+1)*16
  j = size(M,2); j_= j-mod(j,16); if(j_/=j)j_= (j_/16+1)*16
  if(verbose) write(*,*) 'M allocate array ,i ,j : ', i ,j
  if(verbose) write(*,*) 'M allocate array ,i_,j_: ', i_,j_
  allocate(M_(i_,j_));
  if(verbose) write(*,*) 'allocated,padd with zeros'

  M_=0.d0
  M_(1:i,1:j)=M(1:i,1:j)

  if(verbose) write(*,*) 'reshape matrix N'
  i = size(N,1); i_= i-mod(i,16); if(i_/=i)i_= (i_/16+1)*16
  j = size(N,2); j_= j-mod(j,16); if(j_/=j)j_= (j_/16+1)*16
  if(verbose) write(*,*) 'N allocate array ,i ,j : ', i ,j
  if(verbose) write(*,*) 'N allocate array ,i_,j_: ', i_,j_
  allocate(N_(i_,j_));

  N_=0.d0
  N_(1:i,1:j)=N(1:i,1:j)

  i=size(M_,1);j=size(M_,2);k=size(N_,2)
  allocate(Q_(i,k))
  if(verbose)then
  write(*,*) 'CALLING MATMUL CUDA, dims : ',i,j,k
  write(*,*) 'shape M ',shape(M_)
  write(*,*) 'shape N ',shape(N_)
  write(*,*) 'shape Q ',shape(Q_)
  write(*,*) 'ijk     ',i,j,k
  endif
  CALL matmul_cuda_c(16,M_,N_,Q_,i,j,k)
  if(verbose) write(*,*) '...done...'
  Q(1:dima,1:dimc)=Q_(1:dima,1:dimc)
  if(allocated(M_)) deallocate(M_); if(allocated(N_)) deallocate(N_);
  if(allocated(Q_)) deallocate(Q_)

return
end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

subroutine invert_my_matrix_nvcc(nb,n,a)
implicit none
interface
 subroutine invert_my_matrix(nb,n,a)
  implicit none
 integer ::nb,n
 complex(8) :: a(n,n)
 end subroutine
end interface
integer ::n ,nb
complex(8) :: a(n,n)
call invert_my_matrix(nb,n,a)
end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

      subroutine matinv_magma_complex(nn,AA)
      implicit none
      integer                       :: nn
      complex(8)                    :: AA(:,:)
      complex(8), allocatable       :: A(:,:), X(:,:)
      integer,    allocatable       :: ipiv(:)
      integer                       :: i, n, info, lda
      integer                       :: nrhs
      real(kind=8)                  :: flops, t
      integer                       :: k,j,tstart(2), tend(2)

      !call cublas_init()
      n = nn; lda = n; nrhs= n
      allocate(A(lda,n)); allocate(X(lda,nrhs)); allocate(ipiv(n));

      A = AA
      if(size(A,1)/=n) stop 'wrong dimension'
      do i=1,n
       do j=1,n
        if(i==j)then
         X(i,j)=1.d0
        else
         X(i,j)=0.d0
        endif
       enddo
      enddo

#ifdef MAGMA
!---- Call magma LU ----------!
      call magmaf_zgetrf(n, n, A, lda, ipiv, info)
      if ( info .ne. 0 )  then
         write(*,*) "Info zgetrf: ", info
      end if

!---- Call solve -------------!
      call zgetrs('n', n, nrhs, A, lda, ipiv, X, lda, info)
      if ( info .ne. 0 )  then
         write(*,*) "Info zgetrs: ", info
      end if
      AA=X
      deallocate(A,X,ipiv)
     !call cublas_shutdown()
#endif

 105  format((a35,es10.3))
      end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


      subroutine matinv_magma_double(nn,AA)
      implicit none
      integer :: nn
      real(8) :: AA(:,:)
      double precision, allocatable :: work(:)
      double precision, allocatable :: A(:), X(:)
      integer,    allocatable       :: ipiv(:)
      integer                       :: i, n, info, lda
      integer                       :: nrhs
      real(kind=8)                  :: flops, t
      integer                       :: k,j,tstart(2), tend(2)
     
#ifdef MAGMA 
      !call cublas_init()
      n = nn; lda = n; nrhs= n
      allocate(A(lda*n)); allocate(X(lda*nrhs)); allocate(ipiv(n)); allocate(work(n))

      k=0
      do j=1,n
       do i=1,n
         k=k+1
         A(k) = AA(i,j)
       end do
      enddo
      k=0
      do j=1,n
       do i=1,n
         k=k+1
         X(k) =0
         if(i==j) X(k)=1.d0
       end do
      enddo

!---- Call magma LU ----------!
      call magmaf_dgetrf(n, n, A, lda, ipiv, info)
      if ( info .ne. 0 )  then
         write(*,*) "Info : ", info
      end if

!---- Call solve -------------!
      call dgetrs('n', n, nrhs, A, lda, ipiv, X, lda, info)
      if ( info .ne. 0 )  then
         write(*,*) "Info : ", info
      end if
     
      k=0
      do j=1,n
       do i=1,n
         k=k+1
         AA(i,j)=X(k)
       end do
      enddo
 
      deallocate(A,   X, ipiv, work)
      !call cublas_shutdown()
#endif

 105  format((a35,es10.3))
      end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


subroutine Gpu_infc_LU(h_A,lda,ipiv_gpu)
   implicit none
   integer :: Np
   integer(kind=4), intent(in) :: lda
   integer(kind=4), dimension(:), intent(out) :: ipiv_gpu
   complex(kind=8), dimension(:,:), intent(inout) :: h_A
   integer(kind=4) :: m, n, lda_t, i, j
   character*100 :: var
   integer(kind=4) :: stat, info
   real(kind=8) :: end=0, start=0
   complex(kind=8),allocatable,dimension(:,:) :: h_A2
   integer, dimension(4) :: seed

#ifdef MAGMA

   Np=size(h_A,2)
   allocate(h_A2(Np,Np))
   seed(1)=124; seed(2)=352; seed(3)=753; seed(4)=977
   do i = 1, Np
    call zlarnv(2,seed,Np,h_A(:,i))
   end do
   h_A2=h_A
   m = Np; n = Np; lda_t = lda
   write(*,*) 'Matrix Size: N, N*M = ', Np, Np*Np
   call magmaf_zgetrf(m, n, h_A, lda_t, ipiv_gpu, info) ! Compute on gpu
   call zgetrf(Np,Np,h_A2,lda,ipiv_gpu,info)
   write(*,*) 'h_A,h_A2 =',maxval(abs(h_A-h_A2))
   if (info <  0) then
     write (6,*) 'magma: the i-th argument had an illegal value',INFO
     if(info .eq. -7) write(6,*) 'internal GPU memory allocation failed.'
   else if(info > 0) then
     write (6,*) 'magma: U(i,i) is exactly zero.',INFO
   end if

#endif

end subroutine 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  subroutine magma_fortran_comp(nn,aa,bb)
  integer    :: nn
  complex(8) :: aa(nn,nn),bb(nn,nn)
  interface
   subroutine magma_inv_routine_comp_all_gpu(nn,aa,bb)
   integer    :: nn
   complex(8) :: aa(nn,nn),bb(nn,nn)
   end subroutine
  end interface
#ifdef MAGMA
  call magma_inv_routine_comp_all_gpu(nn,aa,bb)
#endif
  end subroutine

  subroutine magma_fortran_comp_(nn,aa)
  integer    :: nn,i
  complex(8) :: aa(nn,nn),bb(nn,nn)
  interface
   subroutine magma_inv_routine_comp_all_gpu(nn,aa,bb)
   integer    :: nn
   complex(8) :: aa(nn,nn),bb(nn,nn)
   end subroutine
  end interface
   bb=0.d0; do i=1,size(aa,1); bb(i,i)=1.d0 ; enddo
#ifdef MAGMA
   call magma_inv_routine_comp_all_gpu(nn,aa,bb)
#endif
   aa=bb
  end subroutine

 subroutine matmulcuda_c_singleprec_bypass(M,N,Q,dima,dimb,dimc)
 implicit none
  INTEGER     :: DIMA,DIMB,DIMC
  COMPLEX(8)  :: M(DIMA*DIMB),N(DIMB*DIMC),Q(DIMA*DIMC)
  COMPLEX(4)  :: M_(DIMA*DIMB),N_(DIMB*DIMC),Q_(DIMA*DIMC)
  interface
   subroutine magma_complex_matrix_multiply_single(M,N,Q,a,b,c)
    implicit none
    integer :: a,b,c
    complex(4) :: M(a*b),N(b*c),Q(a*c)
   end subroutine
  end interface
  M_=M; N_=N
#ifdef MAGMA
  CALL magma_complex_matrix_multiply_single(M_,N_,Q_,DIMA,DIMB,DIMC)
#endif
  Q=Q_
 end subroutine

 subroutine matmulcuda_r_singleprec_bypass(M,N,Q,dima,dimb,dimc)
 implicit none
  INTEGER     :: DIMA,DIMB,DIMC
  real(8)  :: M(DIMA*DIMB),N(DIMB*DIMC),Q(DIMA*DIMC)
  real(4)  :: M_(DIMA*DIMB),N_(DIMB*DIMC),Q_(DIMA*DIMC)
  interface
   subroutine magma_real_matrix_multiply_single(M,N,Q,a,b,c)
    implicit none
    integer :: a,b,c
    real(4) :: M(a*b),N(b*c),Q(a*c)
   end subroutine
  end interface
  M_=M; N_=N
#ifdef MAGMA
  CALL magma_real_matrix_multiply_single(M_,N_,Q_,DIMA,DIMB,DIMC)
#endif
  Q=Q_
 end subroutine

 subroutine magma_fortran_comp_singleprec_bypass(nn,aa)
  integer    :: nn,i
  complex(8) :: aa(nn,nn)
  complex(4) :: bb_single(nn,nn),aa_single(nn,nn)
  interface
   subroutine magma_inv_routine_comp_all_gpu_single(nn,aa,bb)
   integer    :: nn
   complex(4) :: aa(nn,nn),bb(nn,nn)
   end subroutine
  end interface
   aa_single=aa
   bb_single=0.d0; do i=1,size(aa,1); bb_single(i,i)=1.d0 ; enddo
#ifdef MAGMA
   call magma_inv_routine_comp_all_gpu_single(nn,aa_single,bb_single)
#endif
   aa=bb_single
  end subroutine

 subroutine magma_fortran_double_singleprec_bypass(nn,aa)
  integer    :: nn,i
  real(8) :: aa(nn,nn)
  real(4) :: bb_single(nn,nn),aa_single(nn,nn)
  interface
   subroutine magma_inv_routine_double_all_gpu_single(nn,aa,bb)
   integer    :: nn
   real(4) :: aa(nn,nn),bb(nn,nn)
   end subroutine
  end interface
   aa_single=aa
   bb_single=0.d0; do i=1,size(aa,1); bb_single(i,i)=1.d0 ; enddo
#ifdef MAGMA
   call magma_inv_routine_double_all_gpu_single(nn,aa_single,bb_single)
#endif
   aa=bb_single
  end subroutine

  subroutine magma_fortran_double_(nn,aa)
  integer    :: nn,i
  real(8) :: aa(nn,nn),bb(nn,nn)
  interface
   subroutine magma_inv_routine_double_all_gpu(nn,aa,bb)
   integer    :: nn
   real(8) :: aa(nn,nn),bb(nn,nn)
   end subroutine
  end interface
   bb=0.d0; do i=1,size(aa,1); bb(i,i)=1.d0 ; enddo
#ifdef MAGMA
   call magma_inv_routine_double_all_gpu(nn,aa,bb)
#endif
   aa=bb
  end subroutine

  subroutine magma_fortran_double(nn,aa,bb)
   integer    :: nn
   real(8) :: aa(nn,nn),bb(nn,nn)
  interface
   subroutine magma_inv_routine_double_all_gpu(nn,aa,bb)
   integer    :: nn
   real(8) :: aa(nn,nn),bb(nn,nn)
   end subroutine
  end interface
#ifdef MAGMA
  call magma_inv_routine_double_all_gpu(nn,aa,bb)
#endif
  end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


subroutine testing_magma_inv
implicit none
integer                       :: i,kkk,k,kk,ss
complex(8),allocatable        :: aa(:,:),bb(:,:),cc(:,:)
integer,allocatable           :: ipiv(:)

do kkk=1,4
   write(*,*) '========================================'
   call freeMem()
   ss=kkk*20*32
   write(*,*) 'size : ', ss
   if(allocated(aa)) deallocate(aa,bb,ipiv,cc)
   allocate(aa(ss,ss),bb(ss,ss),ipiv(ss),cc(ss,ss))
  
   bb=0.
   do i=1,ss
   do k=1,ss
    aa(i,k)=cmplx(1.d0-2.d0*drand1(),1.d0-2.d0*drand1())
   enddo
   !aa(i,i)=aa(i,i)+100.
   enddo
   !aa=aa+transpose(aa)
   
   write(*,*) 'matrix A : ',maxval(abs(aa))
   bb=aa
   write(*,*) 'start'
   call magma_fortran_comp_(ss,aa)
   write(*,*) 'done, now : C=A*B'
   call matmulcuda_c(aa,bb,cc,ss,ss,ss)
   write(*,*) 'C=A*B  : ',maxval(abs(cc)),minval(abs(cc))
 
enddo
 
 contains
 function drand1()
 implicit none
 real(8) :: drand1
 real(4) :: r
  call random_number(r)
  drand1=dble(r)
 end function

end subroutine   

subroutine testing_magma_inv_double
implicit none
integer                       :: i,kkk,k,kk,ss
real(8),allocatable           :: aa(:,:),bb(:,:),cc(:,:)
integer,allocatable           :: ipiv(:)

do kkk=1,4
   write(*,*) '========================================'
   call freeMem()
   ss=kkk*20*32
   write(*,*) 'size : ', ss
   if(allocated(aa)) deallocate(aa,bb,ipiv,cc)
   allocate(aa(ss,ss),bb(ss,ss),ipiv(ss),cc(ss,ss))

   bb=0.
   do i=1,ss
   do k=1,ss
    aa(i,k)=1.d0-2.d0*drand1()
   enddo
   enddo
   !aa=aa+transpose(aa)

   write(*,*) 'matrix A : ',maxval(abs(aa))
   bb=aa
   write(*,*) 'start'
   call magma_fortran_double_(ss,aa)
   write(*,*) 'done, now : C=A*B'
   call matmulcuda_r(aa,bb,cc,ss,ss,ss)
   write(*,*) 'C=A*B  : ',maxval(abs(cc)),minval(abs(cc))
 enddo

 contains
 function drand1()
 implicit none
 real(8) :: drand1
 real(4) :: r
  call random_number(r)
  drand1=dble(r)
 end function

end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

end module


