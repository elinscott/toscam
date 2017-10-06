MODULE Lanczos_gpu_mod

  USE tridiag_class
  USE H_class

  IMPLICIT NONE

  REAL(DBL),    PARAMETER, PRIVATE   ::  zero=0.0_DBL,one=1.0_DBL,rerror=epsilon(1.d0)
  LOGICAL,      PARAMETER, PRIVATE   ::  F=.FALSE.,T=.TRUE.

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
!**************************************************************************
!**************************************************************************

  SUBROUTINE Lanczos_get_GS_sector_GPU(lowest)
    implicit none
    TYPE(eigenlist_type),intent(in)     :: lowest
    TYPE(eigen_type),POINTER            :: eigen => NULL()
    INTEGER                             :: jj,sizevec

    if(lowest%neigen==0) return
    sizevec=dimen_H()
 
     write(log_unit,*) '====================================================='
     write(log_unit,*) '  GETTING THE GROUD STATE                            '
     write(log_unit,*) '  N eigenvalues    : ', lowest%neigen
     write(log_unit,*) '  iter convergence : ', lowest%eigen(:)%lanczos_iter
     write(log_unit,*) '  eigenvalues      : ', lowest%eigen(:)%val
     write(log_unit,*) '  SECTOR dimension : ', lowest%eigen(:)%dim_space
     write(log_unit,*) '====================================================='

    !-----------------------------------------------------------------------------------------!
      do jj=1,lowest%neigen
       eigen=>lowest%eigen(jj)
       CALL eigen_allocate_vec(eigen,N=sizevec)
       if(.not.ON_FLY)then
        IF(ASSOCIATED(sector_h%updo))THEN
         call Lanczos_GPU_get_GS_updo(eigen%lanczos_iter,sizevec,eigen%lanczos_vecp(1:eigen%lanczos_iter),eigen%vec%rc)
        ELSE
         call Lanczos_GPU_get_GS_memory(eigen%lanczos_iter,sizevec,eigen%lanczos_vecp(1:eigen%lanczos_iter),eigen%vec%rc)
        ENDIF
       else
        IF(ASSOCIATED(sector_h%updo))THEN
         STOP 'error Lanczos gpu on the fly not done yet'
        ELSE 
         call Lanczos_GPU_get_GS(eigen%lanczos_iter,sizevec,eigen%lanczos_vecp(1:eigen%lanczos_iter),eigen%vec%rc)
        ENDIF
       endif
      enddo
    !-----------------------------------------------------------------------------------------!

  return
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
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE Lanczos_fast_diagonalize_GPU(lowest)

    TYPE(eigenlist_type), INTENT(INOUT) :: lowest
    REAL(DBL), ALLOCATABLE, TARGET      :: VECP(:,:), VALP(:)
    TYPE(eigen_type)                    :: eigen(Neigen)
    INTEGER                             :: Niter,sizevec
    TYPE(rcvector_type)                 :: initvec
    TYPE(tridiag_type)                  :: Lmatrix,subLmatrix
    INTEGER                             :: Neigen_,jj
    REAL(DBL)                           :: coeff

    Niter    =  MIN(dimen_H(),Nitermax)
    Neigen_  =  min(Neigen,dimen_H())
    Neigen_  =  min(Neigen_,Niter)
    sizevec  =  dimen_H()

    ALLOCATE(VECP(Niter,Niter),VALP(Niter))
    CALL new_tridiag(Lmatrix,Niter) 

    if(dimen_H()==0) stop 'error Hilbert space has 0 dimension'

    !----------------------------------------------------------------------------------------------!

         if(.not.ON_FLY)then
          IF(ASSOCIATED(sector_h%updo))THEN
           if(messages3) write(*,*) 'Lanczos GPU : updo basis'
           call Lanczos_GPU_updo(Niter,sizevec,Lmatrix%diag,Lmatrix%subdiag)
          ELSE
           if(messages3) write(*,*) 'Lanczos GPU : memory version'
           call Lanczos_GPU_memory(Niter,sizevec,Lmatrix%diag,Lmatrix%subdiag)
          ENDIF
         else
          IF(ASSOCIATED(sector_h%updo))THEN
           write(*,*) 'on fly for up and dn spin with gpu not yet done'
           stop 'on fly updo not done'
          ELSE
           if(messages3) write(*,*) 'Lanczos GPU on fly'
           call Lanczos_GPU(Niter,sizevec,Lmatrix%diag,Lmatrix%subdiag)
          ENDIF
         endif
  
         CALL   submatrix_tridiag(subLmatrix,Lmatrix,(/1,Niter/))
         CALL diagonalize_tridiag(subLmatrix,VALP(1:Niter),VECP(1:Niter,1:Niter),EIGENVAL_ONLY=F)
         CALL      delete_tridiag(subLmatrix)

         do jj=1,Neigen_
           CALL new_eigen(eigen(jj),VALP(jj),initvec,.true.,RANK=jj,no_vector=.true.)
           eigen(jj)%lanczos_iter=Niter
           eigen(jj)%lanczos_vecp(1:eigen(jj)%lanczos_iter)=VECP(1:Niter,jj)
           eigen(jj)%dim_space=dimen_H()
           CALL add_eigen(eigen(jj),lowest)
         enddo
    !----------------------------------------------------------------------------------------------!

    IF(ALLOCATED(VALP)) DEALLOCATE(VALP);  IF(ALLOCATED(VECP))    DEALLOCATE(VECP) 
    do jj=1,Neigen_
      CALL delete_eigen(eigen(jj))
    enddo
    CALL delete_tridiag(Lmatrix); CALL delete_tridiag(subLmatrix)

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
!**************************************************************************
!**************************************************************************
!**************************************************************************

end module
