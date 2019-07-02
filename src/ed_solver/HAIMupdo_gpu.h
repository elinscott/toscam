!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  subroutine Lanczos_GPU_get_GS_updo(Niter, sizvec, vecp, GS)

     integer    :: Niter, sizvec, start_diagH, n1, n2, n3, n4, n5, n6, n7
     real(kind=DP)    :: vecp(Niter)
#ifndef _complex
     real(kind=DP)    :: GS(sizvec)
#else
     complex(kind=DP) :: GS(sizvec)
#endif
     interface
#ifndef _complex
subroutine lanczos_real_updo_gs_cuda(norbs, block, Niter, ioffdiagup, ioffdiagdn, irankoffup, irankoffdn, sizvec, sizup, sizdn, U, &
                           &  diagup, diagdn, noffup, noffdn, rankoffup, rankoffdn, offdiagup, offdiagdn, UMASK, stateup, statedn &
                           & , iorbup, iorbdn, vecp, gs, rank)
           implicit none
  integer :: block, Niter, sizup, sizdn, ioffdiagup, ioffdiagdn, irankoffup, irankoffdn, sizvec, noffup(sizup), noffdn(sizdn), norbs
         integer :: rank, rankoffup(irankoffup), rankoffdn(irankoffdn), stateup(sizup), statedn(sizdn), iorbup(norbs), iorbdn(norbs)
           real(kind=DP) :: U(norbs*norbs), diagup(sizdn), diagdn(sizup), offdiagup(ioffdiagup), offdiagdn(ioffdiagdn)
           real(kind=DP) :: vecp(Niter), gs(sizvec)
           logical :: UMASK(norbs*norbs)
        end subroutine
#else
      subroutine lanczos_complex_updo_gs_cuda(norbs,block,Niter,ioffdiagup,ioffdiagdn,irankoffup,irankoffdn,sizvec,sizup,sizdn,U, &
                           &  diagup, diagdn, noffup, noffdn, rankoffup, rankoffdn, offdiagup, offdiagdn, UMASK, stateup, statedn &
                           & , iorbup, iorbdn, vecp, gs, rank)
           implicit none
       integer    :: block,Niter,sizup,sizdn,ioffdiagup,ioffdiagdn,irankoffup,irankoffdn,sizvec,noffup(sizup),noffdn(sizdn),norbs
      integer    :: rank, rankoffup(irankoffup), rankoffdn(irankoffdn), stateup(sizup), statedn(sizdn), iorbup(norbs), iorbdn(norbs)
           real(kind=DP)    :: U(norbs*norbs), diagup(sizdn), diagdn(sizup)
           complex(kind=DP) :: offdiagup(ioffdiagup), offdiagdn(ioffdiagdn)
           real(kind=DP)    :: vecp(Niter)
           complex(kind=DP) :: gs(sizvec)
           logical    :: UMASK(norbs*norbs)
        end subroutine
#endif
     end interface

     CALL reset_timer(start_diagH)
     write (log_unit, *) '.... start Lanczos on gpu, go get GS ....'

     n1 = size(U, 1)
     n2 = size(offdiagup)
     n3 = size(offdiagdo)
     n4 = size(rankoffup)
     n5 = size(rankoffdo)
     n6 = size(noffup)
     n7 = size(noffdo)

#ifndef _complex
     call lanczos_real_updo_gs_cuda(n1, cuda_blocksize, Niter, n2, n3, n4, n5, &
                          & dimen, n6, n7, U, diagup, diagdo, noffup, noffdo, rankoffup, rankoffdo, offdiagup, offdiagdo, UMASK, &
                          & sector%up%state, sector%down%state, IMPiorbup, IMPiorbdo, vecp, GS, rank)
#else
     call lanczos_complex_updo_gs_cuda(n1, cuda_blocksize, Niter, n2, n3, n4, n5, &
                           & dimen, n6, n7, U, diagup, diagdo, noffup, noffdo, rankoffup, rankoffdo, offdiagup, offdiagdo, UMASK, &
                           & sector%up%state, sector%down%state, IMPiorbup, IMPiorbdo, vecp, GS, rank)
#endif

     write (log_unit, *) '.... fast lanczos gpu converged ....'
     CALL timer_fortran(start_diagH, "# GETTING EIGENVECTORS H TOOK "//c2s(i2c(Niter))//" ITERATIONS AND ")

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

  subroutine Lanczos_dynamic_GPU_updo(Niter, sizvec, diag, subdiag, inputvec)
     implicit none
     integer          :: Niter, sizvec, start_diagH, n1, n2, n3, n4, n5, n6, n7
     real(kind=DP)          :: diag(Niter), subdiag(Niter)
#ifndef _complex
     real(kind=DP)          :: inputvec(sizvec)
#else
     complex(kind=DP)       :: inputvec(sizvec)
#endif
     !---------------------------------------------------------------------------------!
     interface
#ifndef _complex
      subroutine lanczos_real_updo_dynamic_cuda(norbs,block,Niter,ioffdiagup,ioffdiagdn,irankoffup,irankoffdn,sizvec,sizup,sizdn,U, &
             &  diagup, diagdn, noffup, noffdn, rankoffup, rankoffdn, offdiagup, offdiagdn, diag, subdiag, UMASK, stateup, statedn &
                           & , iorbup, iorbdn, vecinit, rank)
           implicit none
  integer :: block, Niter, sizup, sizdn, ioffdiagup, ioffdiagdn, irankoffup, irankoffdn, sizvec, noffup(sizup), noffdn(sizdn), norbs
         integer :: rank, rankoffup(irankoffup), rankoffdn(irankoffdn), stateup(sizup), statedn(sizdn), iorbup(norbs), iorbdn(norbs)
  real(kind=DP) :: U(norbs*norbs), diagup(sizdn), diagdn(sizup), diag(Niter), subdiag(Niter), offdiagup(ioffdiagup), offdiagdn(ioffdiagdn)
           real(kind=DP) :: vecinit(sizvec)
           logical :: UMASK(norbs*norbs)
        end subroutine
#else
     subroutine lanczos_complex_updo_dynamic_cuda(norbs,block,Niter,ioffdiagup,ioffdiagdn,irankoffup,irankoffdn,sizvec,sizup,sizdn,U, &
             &  diagup, diagdn, noffup, noffdn, rankoffup, rankoffdn, offdiagup, offdiagdn, diag, subdiag, UMASK, stateup, statedn &
                            & , iorbup, iorbdn, vecinit, rank)
           implicit none
       integer    :: rank,block,Niter,sizup,sizdn,ioffdiagup,ioffdiagdn,irankoffup,irankoffdn,sizvec,noffup(sizup),noffdn(sizdn),norbs
           integer    :: rankoffup(irankoffup), rankoffdn(irankoffdn), stateup(sizup), statedn(sizdn), iorbup(norbs), iorbdn(norbs)
           real(kind=DP)    :: U(norbs*norbs), diagup(sizdn), diagdn(sizup), diag(Niter), subdiag(Niter)
           complex(kind=DP) :: offdiagup(ioffdiagup), offdiagdn(ioffdiagdn)
           complex(kind=DP) :: vecinit(sizvec)
           logical    :: UMASK(norbs*norbs)
        end subroutine
#endif
     end interface
     !---------------------------------------------------------------------------------!

     CALL reset_timer(start_diagH)
     if (sizvec == 0) stop 'error Hilbert space has 0 dimension'
     write (log_unit, *) '.... start Lanczos on gpu ....'
     n1 = size(U, 1)
     n2 = size(offdiagup)
     n3 = size(offdiagdo)
     n4 = size(rankoffup)
     n5 = size(rankoffdo)
     n6 = size(noffup)
     n7 = size(noffdo)

#ifndef _complex
     call lanczos_real_updo_dynamic_cuda(n1, cuda_blocksize, Niter, n2, n3, n4, n5, &
             & dimen, n6, n7, U, diagup, diagdo, noffup, noffdo, rankoffup, rankoffdo, offdiagup, offdiagdo, diag, subdiag, UMASK, &
                          & sector%up%state, sector%down%state, IMPiorbup, IMPiorbdo, inputvec, rank)
#else
     call lanczos_complex_updo_dynamic_cuda(n1, cuda_blocksize, Niter, n2, n3, n4, n5, &
             & dimen, n6, n7, U, diagup, diagdo, noffup, noffdo, rankoffup, rankoffdo, offdiagup, offdiagdo, diag, subdiag, UMASK, &
                          & sector%up%state, sector%down%state, IMPiorbup, IMPiorbdo, inputvec, rank)
#endif

     write (log_unit, *) '.... fast lanczos gpu converged ....'
     CALL timer_fortran(start_diagH, "# GETTING DYNAMIC TOOK "//c2s(i2c(Niter))//" ITERATIONS AND ")

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

  subroutine Lanczos_GPU_updo(Niter, sizvec, diag, subdiag)
     implicit none
     integer          :: Niter, sizvec, start_diagH, n1, n2, n3, n4, n5, n6, n7
     real(kind=DP)          :: diag(Niter), subdiag(Niter)
     !---------------------------------------------------------------------------------!
     interface
#ifndef _complex
   subroutine lanczos_real_updo_cuda(norbs, block, Niter, ioffdiagup, ioffdiagdn, irankoffup, irankoffdn, sizvec, sizup, sizdn, U, &
             &  diagup, diagdn, noffup, noffdn, rankoffup, rankoffdn, offdiagup, offdiagdn, diag, subdiag, UMASK, stateup, statedn &
                           & , iorbup, iorbdn, rank)
           implicit none
       integer :: rank,block,Niter,sizup,sizdn,ioffdiagup,ioffdiagdn,irankoffup,irankoffdn,sizvec,noffup(sizup),noffdn(sizdn),norbs
           integer :: rankoffup(irankoffup), rankoffdn(irankoffdn), stateup(sizup), statedn(sizdn), iorbup(norbs), iorbdn(norbs)
  real(kind=DP) :: U(norbs*norbs), diagup(sizdn), diagdn(sizup), diag(Niter), subdiag(Niter), offdiagup(ioffdiagup), offdiagdn(ioffdiagdn)
           logical :: UMASK(norbs*norbs)
        end subroutine
#else
subroutine lanczos_complex_updo_cuda(norbs, block, Niter, ioffdiagup, ioffdiagdn, irankoffup, irankoffdn, sizvec, sizup, sizdn, U, &
             &  diagup, diagdn, noffup, noffdn, rankoffup, rankoffdn, offdiagup, offdiagdn, diag, subdiag, UMASK, stateup, statedn &
                            & , iorbup, iorbdn, rank)
           implicit none
       integer    :: block,Niter,sizup,sizdn,ioffdiagup,ioffdiagdn,irankoffup,irankoffdn,sizvec,noffup(sizup),noffdn(sizdn),norbs
      integer    :: rank, rankoffup(irankoffup), rankoffdn(irankoffdn), stateup(sizup), statedn(sizdn), iorbup(norbs), iorbdn(norbs)
           real(kind=DP)    :: U(norbs*norbs), diagup(sizdn), diagdn(sizup), diag(Niter), subdiag(Niter)
           complex(kind=DP) :: offdiagup(ioffdiagup), offdiagdn(ioffdiagdn)
           logical    :: UMASK(norbs*norbs)
        end subroutine
#endif
     end interface
     !---------------------------------------------------------------------------------!

     if (offdiag_coulomb) then
        write (*, *) 'offdiag coulomb not yet implemented in GPU Lanczos'
        stop
     endif

     CALL reset_timer(start_diagH)
     write (log_unit, *) '.... start Lanczos on gpu ....'
     n1 = size(U, 1)
     n2 = size(offdiagup)
     n3 = size(offdiagdo)
     n4 = size(rankoffup)
     n5 = size(rankoffdo)
     n6 = size(noffup)
     n7 = size(noffdo)

#ifndef _complex
     call lanczos_real_updo_cuda(n1, cuda_blocksize, Niter, n2, n3, n4, n5, &
             & dimen, n6, n7, U, diagup, diagdo, noffup, noffdo, rankoffup, rankoffdo, offdiagup, offdiagdo, diag, subdiag, UMASK, &
                          & sector%up%state, sector%down%state, IMPiorbup, IMPiorbdo, rank)
#else
     call lanczos_complex_updo_cuda(n1, cuda_blocksize, Niter, n2, n3, n4, n5, &
             & dimen, n6, n7, U, diagup, diagdo, noffup, noffdo, rankoffup, rankoffdo, offdiagup, offdiagdo, diag, subdiag, UMASK, &
                           & sector%up%state, sector%down%state, IMPiorbup, IMPiorbdo, rank)
#endif
     write (log_unit, *) '.... fast lanczos gpu converged ....'
     CALL timer_fortran(start_diagH, "#GETTING TRIDIAG MATRIX TOOK "//c2s(i2c(Niter))//" ITERATIONS AND ")
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
