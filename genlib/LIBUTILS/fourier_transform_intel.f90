module FFTW
  use MKL_DFTI
  use MKL_DFT_TYPE
  implicit none 
  private 
  public cfft_rw2rt,manip_fftrw2rt,cfft_rt2rw,cfft_iw2it,extract
contains
  !+-------------------------------------------------------------------+
  !PROGRAM  : CFFT_RW2RT
  !TYPE     : Subroutine
  !PURPOSE  : Evaluate the FFT of a given function from real frequencies
  ! to real time.
  !COMMENT  : n=number of frequencies should be power of 2: log_2(n)=integer
  ! rescaling and reshaping (1,2*L --> -L:L) should be implemented afterward
  !+-------------------------------------------------------------------+
  subroutine cfft_rw2rt(func,n)
    integer    :: i,n,status
    complex(8) :: func(n)
    type(DFTI_DESCRIPTOR), pointer :: Handle
    Status = DftiCreateDescriptor(Handle,DFTI_DOUBLE,DFTI_COMPLEX,1,n)
    Status = DftiCommitDescriptor(Handle)
    Status = DftiComputeForward(Handle,func)
  end subroutine cfft_rw2rt

  !+-------------------------------------------------------------------+
  !PROGRAM  : MANIP_FFTRW2RT
  !TYPE     : Subroutine
  !PURPOSE  : Implement the manipulations necessary for re-shaping the 
  ! array after w-->t FFT (cfft_rw2rt) 
  !COMMENT  : 
  !+-------------------------------------------------------------------+
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
  end subroutine manip_fftrw2rt


  !+-------------------------------------------------------------------+
  !PROGRAM  : CFFT_RT2RW
  !TYPE     : Subroutine
  !PURPOSE  : Evaluate the FFT of a given function on from real time
  ! to real frequencies. 
  !COMMENT  : Some manipulations are needed after calling of this subroutine
  ! to reshape the arrays. These are implemented in manip_fftrw2rt
  !+-------------------------------------------------------------------+
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
  end subroutine cfft_rt2rw


  !+-------------------------------------------------------------------+
  !PROGRAM  : CFFT_IW2IT
  !TYPE     : Subroutine
  !PURPOSE  : Evaluate the FFT of a given function from Matsubara frequencies
  ! to imaginary time. 
  !"fg"= the GF to FFT. dimension  NMAX (fg(NMAX))
  !"funct"= output function. dimension 0:n (funct(0:N))
  !"beta"= inverse temperature (directly passed)
  !COMMENT  : The routine implements all the necessary manipulations required 
  !by the FFT in this formalism: tail subtraction and reshaping. 
  !+-------------------------------------------------------------------+
  subroutine cfft_iw2it(fg,funct,beta)
    implicit none
    integer :: i,n,nmax,nlimit
    real(8) :: pi,wmax,beta,mues,tau,dtau,At,w
    complex(8) :: tail
    complex(8),dimension(:)   :: fg
    complex(8),dimension(:),allocatable :: dummy
    real(8),dimension(0:)   :: funct
    nmax=size(fg)
    n=size(funct)-1
    allocate(dummy(2*n))
    pi=acos(-1.d0)
    dtau=beta/dble(n)
    wmax=pi/beta*(2.d0*dble(n)-1.d0)
    mues=-real(fg(n))*wmax**2
    do i=1,n
       w=pi/beta*(2.d0*dble(i)-1.d0)
       tail=-(mues+w*(0.d0,1.d0))/(mues**2+w**2)
       if(i >nmax)fg(i)=tail
       dummy(2*i)= fg(i)-tail
       dummy(2*i-1)=(0.d0,0.d0)
    enddo
    call cfft_rw2rt(dummy,2*n)
    dummy=2.d0*dummy/beta
    do i=0,n-1
       tau=dfloat(i)*dtau
       At=-exp((beta-tau)*mues)/(exp(beta*mues)+1.d0)
       funct(i)=real(dummy(i+n+1))-At
    enddo
    funct(n)=1.d0-funct(0)
    deallocate(dummy)
  end subroutine cfft_iw2it

  !+-------------------------------------------------------------------+
  !PROGRAM  : EXTRACT
  !TYPE     : Subroutine
  !PURPOSE  : Sample a given function G(tau) over Nfak < N points.
  !COMMENTS : Incoming function is expected to have N+1 points (tau=0,beta)
  !this is modification with respect to the std extract routine used in
  !HFqmc.
  ! g0 has N+1 points
  ! g00 will have Nfak+1 (tau_fak=0,beta)
  !+-------------------------------------------------------------------+
  subroutine extract(g0,g00)
    real(8),dimension(0:)  :: g0 !0:L
    real(8),dimension(0:) :: g00 !0:Lfak
    integer :: N,Nfak
    integer :: i,ip
    real(8) :: p,mismatch

    N=size(g0)-1
    Nfak=size(g00)-1
    !Fix the end points
    g00(0)=g0(0)
    g00(Nfak)=1.d0-g0(0)
    mismatch=dble(N)/dble(Nfak)
    do i=1,Nfak-1
       p=dble(i)*mismatch
       ip=int(p)
       g00(i)=g0(ip)
    enddo
    return
  end subroutine extract
end module FFTW
