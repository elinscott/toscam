 module fast_fouriermod

  use genvar
  use linalg
  use mkl_dfti
  use mkl_trig_transforms
  use MKL_DFT_TYPE


  private
  public cfft_rw2rt,manip_fftrw2rt,cfft_rt2rw,fftw

 contains

!----------------------------!
!----------------------------!
!----------------------------!
!----------------------------!
!----------------------------!

  subroutine cfft_rw2rt(func,n)
    integer    :: i,n,status
    complex(8) :: func(n)
    type(DFTI_DESCRIPTOR), pointer :: Handle
    Status = DftiCreateDescriptor(Handle,DFTI_DOUBLE,DFTI_COMPLEX,1,n)
    Status = DftiCommitDescriptor(Handle)
    Status = DftiComputeForward(Handle,func)
  end subroutine

!----------------------------!
!----------------------------!
!----------------------------!
!----------------------------!
!----------------------------!

  subroutine manip_fftrw2rt(func_in,func_out,nhalf)
    implicit none
    integer :: n,i,nhalf
    complex(8),dimension(2*nhalf) :: func_in
    complex(8),dimension(-nhalf:nhalf) :: func_out,dummy
    n=2*nhalf
    !1) [1,2*L=n]---> [-L,L-1]
    do i=1,n
       dummy(i-nhalf-1)=func_in(i)
    enddo
    !2) g[0,L-1]<--- x[-L,-1]
    do i=-nhalf,-1
       func_out(i+nhalf)=dummy(i)
    enddo
    !3) g[-L,-1]<--- x[0,L-1]
    do i=0,nhalf-1
       func_out(i-nhalf)=dummy(i)
    enddo
  end subroutine 

!----------------------------!
!----------------------------!
!----------------------------!
!----------------------------!
!----------------------------!

  subroutine cfft_rt2rw(func,n)
    integer    :: i,n,status
    complex(8) :: func(n)
    real(8)    :: ex
    type(DFTI_DESCRIPTOR), pointer :: Handle

    Status = DftiCreateDescriptor(Handle,DFTI_DOUBLE,DFTI_COMPLEX,1,n)
    Status = DftiCommitDescriptor(Handle)
    Status = DftiComputeBackward(Handle,func)
    ex=-1.d0
    do i=1,n
       ex=-ex
       func(i)=ex*func(i)
    enddo
  end subroutine

!----------------------------!
!----------------------------!
!----------------------------!
!----------------------------!
!----------------------------!
 
 subroutine fftw(n,x,f,y,func)
 implicit none

 interface
  real(8) function func(x)
   real(8) :: x
  end function
 end interface
 optional :: func                                                                                                                                                                

 integer                        :: n,i,k,tt_type
 integer                        :: ir,ipar(128)
 double precision               :: x(n+1)
 real(8),optional               :: y(n+1)
 real(8)                        :: f(n+1), dpar(3*n/2+1)
 type(dfti_descriptor), pointer :: handle

 if(present(y))    f=y    
 if(present(func)) f=(/( func(x(i)), i=1,size(x))/)       
 
           CALL D_INIT_TRIG_TRANSFORM(n,tt_type,ipar,dpar,ir)
           if (ir.ne.0) goto 99
           CALL D_COMMIT_TRIG_TRANSFORM(f,handle,ipar,dpar,ir)
           if (ir.ne.0) goto 99
           CALL D_FORWARD_TRIG_TRANSFORM(f,handle,ipar,dpar,ir)
           if (ir.ne.0) goto 99
           CALL D_BACKWARD_TRIG_TRANSFORM(f,handle,ipar,dpar,ir)
           if (ir.ne.0) goto 99
           CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)
           if (ir.ne.0) goto 99
           
      return

99    continue
      write(*,*) 'FAILED to compute the solution(s)...'

 return
 end subroutine

!----------------------------!
!----------------------------!
!----------------------------!
!----------------------------!
!----------------------------!

end module
