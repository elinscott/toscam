!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  subroutine Lanczos_GPU_get_GS_memory(Niter,sizvec,vecp,GS)

    integer    :: Niter,sizvec,start_diagH,n1,n2
    real(8)    :: vecp(Niter)
#ifndef _complex
    real(8)    :: GS(sizvec)
#else
    complex(8) :: GS(sizvec)
#endif
    interface
#ifndef _complex
      subroutine Lanczos_Real_get_GS_cuda(block,Niter,ioffdiagsz,irankoffsz,sizvec,QUART, &
                                       & diagsz,noffsz,rankoffsz,offdiagsz,vecp,GS,rank)
       integer :: rank,Niter,ioffdiagsz,irankoffsz,sizvec,noffsz(sizvec),rankoffsz(irankoffsz),block
       real(8) :: QUART(sizvec),diagsz(sizvec),vecp(Niter),offdiagsz(ioffdiagsz),GS(sizvec)
      end subroutine
#else
      subroutine Lanczos_get_GS_cuda_complex(block,Niter,ioffdiagsz,irankoffsz,sizvec,QUART, &
                                       & diagsz,noffsz,rankoffsz,offdiagsz,vecp,GS,rank)
       integer    :: rank,Niter,ioffdiagsz,irankoffsz,sizvec,noffsz(sizvec),rankoffsz(irankoffsz),block
       real(8)    :: QUART(sizvec),diagsz(sizvec),vecp(Niter)
       complex(8) :: offdiagsz(ioffdiagsz),GS(sizvec)
      end subroutine
#endif
    end interface

    CALL reset_timer(start_diagH)
    n1=size(offdiagsz);n2=size(rankoffsz)
    if(sizvec==0) stop 'error Hilbert space has 0 dimension'
    write(log_unit,*) '.... start Lanczos on gpu, go get GS ....'
#ifndef _complex
    call Lanczos_Real_get_GS_cuda(cuda_blocksize,Niter,n1,n2,sizvec,QUART_INT_SZ,diagsz,noffsz,rankoffsz,offdiagsz,vecp,GS,rank)
#else
    call Lanczos_get_GS_cuda_complex(cuda_blocksize,Niter,n1,n2,sizvec,QUART_INT_SZ,diagsz,noffsz,rankoffsz,offdiagsz,vecp,GS,rank)
#endif
    write(log_unit,*) '.... fast lanczos gpu converged ....'
    CALL timer_fortran(start_diagH,"# GETTING EIGENVECTORS H TOOK "//c2s(i2c(Niter))//" ITERATIONS AND ")

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

  subroutine Lanczos_dynamic_GPU_memory(Niter,sizvec,diag,subdiag,inputvec)
    implicit none
     integer          :: Niter,sizvec,start_diagH,n1,n2
     real(8)          :: diag(Niter),subdiag(Niter)
#ifndef _complex
     real(8)          :: inputvec(sizvec)
#else
     complex(8)       :: inputvec(sizvec)
#endif
   !---------------------------------------------------------------------------------!
    interface
#ifndef _complex
      subroutine lanczos_real_dynamic_cuda(block,Niter,ioffdiagsz,irankoffsz,sizvec,QUART, &
                         & diagsz,noffsz,rankoffsz,offdiagsz,diag,subdiag,inputvec,rank)
       implicit none
       integer :: rank,Niter,ioffdiagsz,irankoffsz,sizvec,noffsz(sizvec),rankoffsz(irankoffsz),block
       real(8) :: QUART(sizvec),diagsz(sizvec),diag(Niter),subdiag(Niter),offdiagsz(ioffdiagsz),inputvec(sizvec)
      end subroutine
#else
      subroutine lanczos_dynamic_cuda_complex(block,Niter,ioffdiagsz,irankoffsz,sizvec,QUART, &
                         & diagsz,noffsz,rankoffsz,offdiagsz,diag,subdiag,inputvec,rank)
       implicit none
       integer    :: rank,Niter,ioffdiagsz,irankoffsz,sizvec,noffsz(sizvec),rankoffsz(irankoffsz),block
       real(8)    :: QUART(sizvec),diagsz(sizvec),diag(Niter),subdiag(Niter)
       complex(8) :: offdiagsz(ioffdiagsz),inputvec(sizvec)
      end subroutine
#endif
    end interface
   !---------------------------------------------------------------------------------!

    CALL reset_timer(start_diagH)
    if(sizvec==0) stop 'error Hilbert space has 0 dimension'
    write(log_unit,*) '.... start Lanczos on gpu ....'
    n1=size(offdiagsz);n2=size(rankoffsz)
    if(size(QUART_INT_SZ)/=size(inputvec)) stop 'error Lanczos Real GPU : dimensions do not match'
#ifndef _complex
    call lanczos_real_dynamic_cuda(cuda_blocksize,Niter,n1,n2,sizvec,QUART_INT_SZ,diagsz,noffsz,rankoffsz,offdiagsz,diag,subdiag,inputvec,rank)
#else
    call lanczos_dynamic_cuda_complex(cuda_blocksize,Niter,n1,n2,sizvec,QUART_INT_SZ,diagsz,noffsz,rankoffsz,offdiagsz,diag,subdiag,inputvec,rank)
#endif
    write(log_unit,*) '.... fast lanczos gpu converged ....'
    CALL timer_fortran(start_diagH,"# GETTING DYNAMIC TOOK "//c2s(i2c(Niter))//" ITERATIONS AND ")

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

  subroutine Lanczos_GPU_memory(Niter,sizvec,diag,subdiag)
    implicit none
     integer          :: Niter,sizvec,start_diagH,n1,n2
     real(8)          :: diag(Niter),subdiag(Niter)
   !---------------------------------------------------------------------------------!
    interface
#ifndef _complex
      subroutine Lanczos_Real_cuda(block,Niter,ioffdiagsz,irankoffsz,sizvec,QUART, &
                         & diagsz,noffsz,rankoffsz,offdiagsz,diag,subdiag,rank)
       implicit none
       integer :: rank,block,Niter,ioffdiagsz,irankoffsz,sizvec,noffsz(sizvec),rankoffsz(irankoffsz)
       real(8) :: QUART(sizvec),diagsz(sizvec),diag(Niter),subdiag(Niter),offdiagsz(ioffdiagsz)
      end subroutine
#else
      subroutine Lanczos_cuda_complex(block,Niter,ioffdiagsz,irankoffsz,sizvec,QUART, &
                         & diagsz,noffsz,rankoffsz,offdiagsz,diag,subdiag,rank)
       implicit none
       integer    :: rank,block,Niter,ioffdiagsz,irankoffsz,sizvec,noffsz(sizvec),rankoffsz(irankoffsz)
       real(8)    :: QUART(sizvec),diagsz(sizvec),diag(Niter),subdiag(Niter)
       complex(8) :: offdiagsz(ioffdiagsz)
      end subroutine
#endif
    end interface
   !---------------------------------------------------------------------------------!

    CALL reset_timer(start_diagH)
    if(sizvec==0) stop 'error Hilbert space has 0 dimension'
    write(log_unit,*) '.... start Lanczos on gpu ....'
    n1=size(offdiagsz);n2=size(rankoffsz)
#ifndef _complex
    if(messages3) write(*,*) 'GPU : call to Lanczos Real Cuda'
    call Lanczos_Real_cuda(cuda_blocksize,Niter,n1,n2,sizvec,QUART_INT_SZ,diagsz,noffsz,rankoffsz,offdiagsz,diag,subdiag,rank)
#else
    if(messages3) write(*,*) 'GPU : call to Lanczos Aimag Cuda'
    call Lanczos_cuda_complex(cuda_blocksize,Niter,n1,n2,sizvec,QUART_INT_SZ,diagsz,noffsz,rankoffsz,offdiagsz,diag,subdiag,rank)
#endif
    write(log_unit,*) '.... fast lanczos gpu converged ....'
    CALL timer_fortran(start_diagH,"#GETTING TRIDIAG MATRIX TOOK "//c2s(i2c(Niter))//" ITERATIONS AND ")

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

  subroutine Lanczos_GPU(Niter,sizvec,diag,subdiag)
    implicit none
     integer          :: Niter,sizvec,start_diagH,n1,ntot
     real(8)          :: diag(Niter),subdiag(Niter)
   !---------------------------------------------------------------------------------!
    interface
#ifndef _complex
      subroutine lanczos_real_fly_cuda(dimen,pblocksize,norbs,Niter_lanczos,ntot,quart,diag,subdiag, &
            &  Eb,Ec,Vbc,sector_states,sector_ranks,bathnorbs,impnorbs,imporbs,bathorbs,maskEb,maskEc,maskVbc,rank) 
       implicit none
       integer :: rank,dimen,pblocksize,Niter_lanczos,norbs,impnorbs,bathnorbs,ntot,sector_states(ntot),sector_ranks(0:dimen-1)
       integer :: imporbs(impnorbs),bathorbs(bathnorbs)
       real(8) :: quart(ntot),subdiag(Niter_lanczos),diag(Niter_lanczos),Eb(bathnorbs*bathnorbs),Ec(impnorbs*impnorbs),Vbc(impnorbs*bathnorbs)
       logical :: maskEb(bathnorbs*bathnorbs),maskEc(impnorbs*impnorbs),maskVbc(impnorbs*bathnorbs)
      end subroutine
#else
      subroutine lanczos_complex_fly_cuda(dimen,pblocksize,norbs,Niter_lanczos,ntot,quart,diag,subdiag, &
            &  Eb,Ec,Vbc,sector_states,sector_ranks,bathnorbs,impnorbs,imporbs,bathorbs,maskEb,maskEc,maskVbc,rank)
       implicit none
       integer    :: rank,dimen,pblocksize,Niter_lanczos,norbs,impnorbs,bathnorbs,ntot,sector_states(ntot),sector_ranks(0:dimen-1)
       integer    :: imporbs(impnorbs),bathorbs(bathnorbs)
       real(8)    :: quart(ntot),subdiag(Niter_lanczos),diag(Niter_lanczos)
       complex(8) :: Eb(bathnorbs*bathnorbs),Ec(impnorbs*impnorbs),Vbc(impnorbs*bathnorbs)
       logical    :: maskEb(bathnorbs*bathnorbs),maskEc(impnorbs*impnorbs),maskVbc(impnorbs*bathnorbs)
      end subroutine
#endif
    end interface
   !---------------------------------------------------------------------------------!

    CALL reset_timer(start_diagH)
    if(sizvec==0) stop 'error Hilbert space has 0 dimension'
    write(log_unit,*) '.... start Lanczos on gpu ....'
    n1=size(sector_sz%rank)
    ntot=size(QUART_INT_SZ)
    if(size(sector_sz%state)/=ntot) stop 'error in Lanczos_GPU size do not match'
#ifndef _complex
  if(messages3) write(*,*) 'GPU : call to Lanczos Fly Real Cuda'
  call lanczos_real_fly_cuda(n1,cuda_blocksize,sector_sz%norbs,Niter,ntot,QUART_INT_SZ,diag,subdiag, &
              &  AIM2sz%Eb%rc%mat,AIM2sz%Ec%rc%mat,AIM2sz%Vbc%rc%mat,sector_sz%state,sector_sz%rank, &
              & AIM2sz%BATHnorbs,AIM2sz%IMPnorbs,AIM2sz%IMPiorb,AIM2sz%BATHiorb,AIM2sz%Eb%rc%MASK%mat, &
              & AIM2sz%Ec%rc%MASK%mat,AIM2sz%Vbc%rc%MASK%mat,rank)
#else
  if(messages3) write(*,*) 'GPU : call to Lanczos Fly Complex Cuda'
  call lanczos_complex_fly_cuda(n1,cuda_blocksize,sector_sz%norbs,Niter,ntot,QUART_INT_SZ,diag,subdiag, &
              &  AIM2sz%Eb%rc%mat,AIM2sz%Ec%rc%mat,AIM2sz%Vbc%rc%mat,sector_sz%state,sector_sz%rank, &
              & AIM2sz%BATHnorbs,AIM2sz%IMPnorbs,AIM2sz%IMPiorb,AIM2sz%BATHiorb,AIM2sz%Eb%rc%MASK%mat, &
              & AIM2sz%Ec%rc%MASK%mat,AIM2sz%Vbc%rc%MASK%mat,rank)
#endif
    write(log_unit,*) '.... fast lanczos gpu converged ....'
    CALL timer_fortran(start_diagH,"#GETTING TRIDIAG MATRIX TOOK "//c2s(i2c(Niter))//" ITERATIONS AND ")

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

  subroutine Lanczos_dynamic_GPU(Niter,sizvec,diag,subdiag,initvec)
    implicit none
     integer          :: Niter,sizvec,start_diagH,n1,ntot
     real(8)          :: diag(Niter),subdiag(Niter)
#ifndef _complex
     real(8)          :: initvec(sizvec)
#else
     complex(8)       :: initvec(sizvec)
#endif
   !---------------------------------------------------------------------------------!
    interface
#ifndef _complex
      subroutine lanczos_real_fly_dynamic_cuda(dimen,pblocksize,norbs,Niter_lanczos,ntot,quart,diag,subdiag, &
            &  Eb,Ec,Vbc,sector_states,sector_ranks,bathnorbs,impnorbs,imporbs,bathorbs,maskEb,maskEc,maskVbc,initvec,rank)
       implicit none
       integer :: dimen,pblocksize,Niter_lanczos,norbs,impnorbs,bathnorbs,ntot,sector_states(ntot),sector_ranks(0:dimen-1),rank
       integer :: imporbs(impnorbs),bathorbs(bathnorbs)
       real(8) :: quart(ntot),subdiag(Niter_lanczos),diag(Niter_lanczos),Eb(bathnorbs*bathnorbs),Ec(impnorbs*impnorbs),Vbc(impnorbs*bathnorbs)
       real(8) :: initvec(ntot)
       logical :: maskEb(bathnorbs*bathnorbs),maskEc(impnorbs*impnorbs),maskVbc(impnorbs*bathnorbs)
      end subroutine
#else
      subroutine lanczos_complex_fly_dynamic_cuda(dimen,pblocksize,norbs,Niter_lanczos,ntot,quart,diag,subdiag, &
            &  Eb,Ec,Vbc,sector_states,sector_ranks,bathnorbs,impnorbs,imporbs,bathorbs,maskEb,maskEc,maskVbc,initvec,rank)
       implicit none
       integer    :: rank,dimen,pblocksize,Niter_lanczos,norbs,impnorbs,bathnorbs,ntot,sector_states(ntot),sector_ranks(0:dimen-1)
       integer    :: imporbs(impnorbs),bathorbs(bathnorbs)
       real(8)    :: quart(ntot),subdiag(Niter_lanczos),diag(Niter_lanczos)
       complex(8) :: Eb(bathnorbs*bathnorbs),Ec(impnorbs*impnorbs),Vbc(impnorbs*bathnorbs)
       complex(8) :: initvec(ntot)
       logical    :: maskEb(bathnorbs*bathnorbs),maskEc(impnorbs*impnorbs),maskVbc(impnorbs*bathnorbs)
      end subroutine
#endif
    end interface
   !---------------------------------------------------------------------------------!

    CALL reset_timer(start_diagH)
    if(sizvec==0) stop 'error Hilbert space has 0 dimension'
    write(log_unit,*) '.... start Lanczos on gpu ....'
    n1=size(sector_sz%rank)
    ntot=size(QUART_INT_SZ)
    if(size(sector_sz%state)/=ntot) stop 'error in Lanczos_GPU size do not match'
#ifndef _complex
  call lanczos_real_fly_dynamic_cuda(n1,cuda_blocksize,sector_sz%norbs,Niter,ntot,QUART_INT_SZ,diag,subdiag, &
              &  AIM2sz%Eb%rc%mat,AIM2sz%Ec%rc%mat,AIM2sz%Vbc%rc%mat,sector_sz%state,sector_sz%rank, &
              & AIM2sz%BATHnorbs,AIM2sz%IMPnorbs,AIM2sz%IMPiorb,AIM2sz%BATHiorb,AIM2sz%Eb%rc%MASK%mat, &
              & AIM2sz%Ec%rc%MASK%mat,AIM2sz%Vbc%rc%MASK%mat,initvec,rank)
#else
  call lanczos_complex_fly_dynamic_cuda(n1,cuda_blocksize,sector_sz%norbs,Niter,ntot,QUART_INT_SZ,diag,subdiag, &
              &  AIM2sz%Eb%rc%mat,AIM2sz%Ec%rc%mat,AIM2sz%Vbc%rc%mat,sector_sz%state,sector_sz%rank, &
              & AIM2sz%BATHnorbs,AIM2sz%IMPnorbs,AIM2sz%IMPiorb,AIM2sz%BATHiorb,AIM2sz%Eb%rc%MASK%mat, &
              & AIM2sz%Ec%rc%MASK%mat,AIM2sz%Vbc%rc%MASK%mat,initvec,rank)
#endif
    write(log_unit,*) '.... fast lanczos gpu converged ....'
    CALL timer_fortran(start_diagH,"#GETTING TRIDIAG MATRIX TOOK "//c2s(i2c(Niter))//" ITERATIONS AND ")

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

subroutine Lanczos_GPU_get_GS(Niter,sizvec,vecp,GS)
    integer    :: Niter,sizvec,start_diagH,n1,n2,ntot
    real(8)    :: vecp(Niter)
#ifndef _complex
    real(8)    :: GS(sizvec)
#else
    complex(8) :: GS(sizvec)
#endif

   !---------------------------------------------------------------------------------!
    interface
#ifndef _complex
      subroutine lanczos_real_fly_gs_cuda(dimen,pblocksize,norbs,Niter_lanczos,ntot,quart, &
            &  Eb,Ec,Vbc,sector_states,sector_ranks,bathnorbs,impnorbs,imporbs,bathorbs,maskEb,maskEc,maskVbc,vecp,GS,rank)
       implicit none
       integer :: dimen,pblocksize,Niter_lanczos,norbs,impnorbs,bathnorbs,ntot,sector_states(ntot),sector_ranks(0:dimen-1)
       integer :: imporbs(impnorbs),bathorbs(bathnorbs),rank
       real(8) :: quart(ntot),Eb(bathnorbs*bathnorbs),Ec(impnorbs*impnorbs)
       real(8) :: Vbc(impnorbs*bathnorbs),vecp(Niter_lanczos),GS(ntot)
       logical :: maskEb(bathnorbs*bathnorbs),maskEc(impnorbs*impnorbs),maskVbc(impnorbs*bathnorbs)
      end subroutine
#else
      subroutine lanczos_complex_fly_gs_cuda(dimen,pblocksize,norbs,Niter_lanczos,ntot,quart, &
            &  Eb,Ec,Vbc,sector_states,sector_ranks,bathnorbs,impnorbs,imporbs,bathorbs,maskEb,maskEc,maskVbc,vecp,GS,rank)
       implicit none
       integer    :: dimen,pblocksize,Niter_lanczos,norbs,impnorbs,bathnorbs,ntot,sector_states(ntot),sector_ranks(0:dimen-1),rank
       integer    :: imporbs(impnorbs),bathorbs(bathnorbs)
       real(8)    :: quart(ntot)
       complex(8) :: Eb(bathnorbs*bathnorbs),Ec(impnorbs*impnorbs),Vbc(impnorbs*bathnorbs),GS(ntot)
       real(8)    :: vecp(Niter_lanczos)
       logical    :: maskEb(bathnorbs*bathnorbs),maskEc(impnorbs*impnorbs),maskVbc(impnorbs*bathnorbs)
      end subroutine
#endif
    end interface
   !---------------------------------------------------------------------------------!

    CALL reset_timer(start_diagH)
    if(sizvec==0) stop 'error Hilbert space has 0 dimension'
    write(log_unit,*) '.... start Lanczos on gpu ....'
    n1=size(sector_sz%rank)
    ntot=size(QUART_INT_SZ)
    if(size(sector_sz%state)/=ntot) stop 'error in Lanczos_GPU size do not match'
#ifndef _complex
  call lanczos_real_fly_gs_cuda(n1,cuda_blocksize,sector_sz%norbs,Niter,ntot,QUART_INT_SZ, &
              &  AIM2sz%Eb%rc%mat,AIM2sz%Ec%rc%mat,AIM2sz%Vbc%rc%mat,sector_sz%state,sector_sz%rank, &
              & AIM2sz%BATHnorbs,AIM2sz%IMPnorbs,AIM2sz%IMPiorb,AIM2sz%BATHiorb,AIM2sz%Eb%rc%MASK%mat, &
              & AIM2sz%Ec%rc%MASK%mat,AIM2sz%Vbc%rc%MASK%mat,vecp,GS,rank)
#else
  call lanczos_complex_fly_gs_cuda(n1,cuda_blocksize,sector_sz%norbs,Niter,ntot,QUART_INT_SZ, &
              &  AIM2sz%Eb%rc%mat,AIM2sz%Ec%rc%mat,AIM2sz%Vbc%rc%mat,sector_sz%state,sector_sz%rank, &
              & AIM2sz%BATHnorbs,AIM2sz%IMPnorbs,AIM2sz%IMPiorb,AIM2sz%BATHiorb,AIM2sz%Eb%rc%MASK%mat, &
              & AIM2sz%Ec%rc%MASK%mat,AIM2sz%Vbc%rc%MASK%mat,vecp,GS,rank)
#endif

    write(log_unit,*) '.... fast lanczos gpu converged ....'
    CALL timer_fortran(start_diagH,"#GETTING TRIDIAG MATRIX TOOK "//c2s(i2c(Niter))//" ITERATIONS AND ")

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

  SUBROUTINE Hmult_sz_real_cuda(vec_out,vec_in)
    REAL(DBL), INTENT(INOUT) :: vec_out(:)
    REAL(DBL), INTENT(IN)    ::  vec_in(:)
    integer                  ::  n1,n2,n3
    interface
     subroutine hmult_sz_real_cuda_rout(block,ioffsz,irank,ntot,QUART,diagsz,vec_in, &
                                  & vec_out,noffsz,rankoffsz,offdiagsz,rank)
      integer :: irank,ioffsz,ntot,rankoffsz(irank),noffsz(ntot),block,rank
      real(8) :: offdiagsz(ioffsz),diagsz(ntot),QUART(ntot),vec_in(ntot),vec_out(ntot)
     end subroutine
    end interface

if(ON_FLY) stop 'on fly not done yet for Hmult sz cuda'
#ifndef _complex
    n1=size(offdiagsz);n2=size(rankoffsz);n3=size(vec_in)
    CALL hmult_sz_real_cuda_rout(cuda_blocksize,n1,n2,n3,QUART_INT_SZ,diagsz,vec_in,vec_out,noffsz,rankoffsz,offdiagsz,rank)
#else
    STOP 'calling hmult sz real cuda, when compiled with complex flag'
#endif

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

  SUBROUTINE Hmult_sz_complex_cuda(vec_out,vec_in)
    COMPLEX(DBL), INTENT(INOUT) :: vec_out(:)
    COMPLEX(DBL), INTENT(IN)    ::  vec_in(:)
    integer                     ::  n1,n2,n3

    interface
     subroutine hmult_sz_complex_cuda_rout(block,ioffsz,irank,ntot,QUART,diagsz,vec_in, &
                                  & vec_out,noffsz,rankoffsz,offdiagsz,rank)
      integer    :: irank,ioffsz,ntot,rankoffsz(irank),noffsz(ntot),block,rank
      complex(8) :: vec_in(ntot),vec_out(ntot),offdiagsz(ioffsz)
      real(8)    :: QUART(ntot),diagsz(ntot)
     end subroutine
    end interface

if(ON_FLY) stop 'on fly not done yet for Hmult sz cuda'
#ifndef _complex
   STOP 'calling hmult sz cuda complex, when compiled with real flag'
#else
   n1=size(offdiagsz);n2=size(rankoffsz);n3=size(vec_in)
   CALL hmult_sz_complex_cuda_rout(cuda_blocksize,n1,n2,n3,QUART_INT_SZ,diagsz,vec_in,vec_out,noffsz,rankoffsz,offdiagsz,rank)
#endif

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
!**************************************************************************

