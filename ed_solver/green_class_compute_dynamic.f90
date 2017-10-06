MODULE green_class_compute_dynamic

  use green_class

  IMPLICIT NONE

  REAL(DBL),    PARAMETER, PRIVATE   :: zero=0.0_DBL,one=1.0_DBL
  LOGICAL,      PARAMETER, PRIVATE   :: F=.FALSE.,T=.TRUE.

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
!**************************************************************************

  SUBROUTINE compute_dynamic(iph,dyn,freq,Opm,stat,title,normvec,iisector,GS,keldysh_level)

    INTEGER,          INTENT(IN)         :: iph
    COMPLEX(DBL),     INTENT(INOUT)      :: dyn(:)
    TYPE(freq_type),  INTENT(IN)         :: freq
    TYPE(eigen_type), INTENT(IN)         :: Opm
    CHARACTER(LEN=*), INTENT(IN)         :: stat,title
    REAL(DBL)                            :: normvec
    TYPE(eigensectorlist_type)           :: GS
    integer                              :: iisector
    INTEGER,optional                     :: keldysh_level

    if(present(keldysh_level))then
      call compute_dynamic_cpu(iph,dyn,freq,Opm,stat,title,normvec,keldysh_level)
      return
    endif
   
    SELECT CASE (which_lanczos)
     CASE('NORMAL')
       call compute_dynamic_cpu(iph,dyn,freq,Opm,stat,title,normvec)
     CASE('ARPACK')
       call compute_dynamic_cpu(iph,dyn,freq,Opm,stat,title,normvec)
     CASE('GPU') 
       call compute_dynamic_gpu(iph,dyn,freq,Opm,stat,title,normvec)
     CASE('FULL_ED')
      if(.not.FLAG_FULL_ED_GREEN)Then
        call compute_dynamic_cpu(iph,dyn,freq,Opm,stat,title,normvec)
      else
        call compute_dynamic_full_ed(iph,dyn,freq,Opm,stat,title,normvec,iisector,GS)
      endif
     CASE DEFAULT
       stop 'error compute dynamic LANCZOS kind not defined'
    END SELECT

  END SUBROUTINE


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE compute_dynamic_full_ed(iph,dyn,freq,Opm,stat,title,normvec,iisector,GS)
    implicit none
    INTEGER,          INTENT(IN)      ::  iph
    COMPLEX(DBL),     INTENT(INOUT)   ::  dyn(:)
    TYPE(freq_type),  INTENT(IN)      ::  freq
    TYPE(eigen_type), INTENT(IN)      ::  Opm
    CHARACTER(LEN=*), INTENT(IN)      ::  stat,title
    COMPLEX(DBL)                      ::  z
    REAL(DBL)                         ::  fermion_sign,ph_sign,norm_Opmvec2,normvec,norm_
    INTEGER                           ::  iw,iter,Niter
    INTEGER                           ::  ii,dimenvec,iisector
    TYPE(eigensectorlist_type)        ::  GS
    LOGICAL                           ::  conv_one_step

    IF(SIZE(dyn)/=freq%Nw) STOP "ERROR IN compute_dynamic: INCONSISTENT DIMENSIONS!"

   !--------------------------------------------------------------------------------!
   ! COMPUTE DYNAMIC CORRELATIONS OF OPERATORS O,O^+ I.E.  <psi_(N-1)| O^+ |psi_N>  !
   !--------------------------------------------------------------------------------1

    IF(iph==1) ph_sign =  one ; IF(iph==2) ph_sign = -one
    SELECT CASE(stat)
      CASE(FERMIONIC)
        fermion_sign = -one
      CASE(BOSONIC)
                   fermion_sign =  one
        if(iph==1) fermion_sign = -one
    END SELECT
   !==============================================================================!
   if(size(Opm%vec%rc)/=dimen_H())then
      write(*,*) '===================================='
      write(*,*) 'size of sector : ', dimen_H()
      write(*,*) 'Opm dim_space  : ', Opm%dim_space
      write(*,*) 'Opm dim vector : ', Opm%dim_sector
      write(*,*) 'size of Opm    : ', size(Opm%vec%rc)
      write(*,*) '===================================='
      STOP 'compute_dynamic inconsistent dimensions, CRITICAL'
   endif

    dimenvec = dimen_H()
    Niter = dimen_H()
    if(abs(normvec)>1.d-10)then
      norm_Opmvec2 = normvec**2
    else
      norm_Opmvec2 = norm_rcvector(Opm%vec)**2
    endif
    if(norm_Opmvec2<cutoff_dynamic) then
      dyn=0.d0; goto 44   
    endif

   if(rank==0)then
    write(*,*) '-----------------------------------'
    write(*,*) ' Observable    = ', title
    write(*,*) 'fermionic sign = ', fermion_sign
    write(*,*) ' case          = ', stat
    write(*,*) ' prefactor     = ', norm_Opmvec2
    write(*,*) 'part/hole part = ', iph
    write(*,*) ' iisectors     = ', iisector
    write(*,*) '-----------------------------------'
   endif

   !==============================================================================!
   !==============================================================================!
   !==============================================================================!
    dyn=0.d0
    call mpibarrier  
    if(iisector/=0)then
     write(*,*) 'nstates = ', GS%es(iisector)%lowest%neigen
     DO iter=rank+1,GS%es(iisector)%lowest%neigen,size2
      DO iw=1,freq%Nw
       if(size(GS%es(iisector)%lowest%eigen(iter)%vec%rc)/=size(Opm%vec%rc)) stop 'error compute dyn full ED size do not match'
       if(size2==1)then
        !open_mp
         norm_=abs(MPI_DOT_PRODUCT(GS%es(iisector)%lowest%eigen(iter)%vec%rc,Opm%vec%rc))**2
       else
         norm_=abs(DOT_PRODUCT(GS%es(iisector)%lowest%eigen(iter)%vec%rc,Opm%vec%rc))**2
       endif
       dyn(iw) = dyn(iw) - fermion_sign/(freq%vec(iw)-ph_sign*(GS%es(iisector)%lowest%eigen(iter)%val-Opm%val))*norm_
      ENDDO
     ENDDO
     call mpisum(dyn)
    endif
   !==============================================================================!
   !==============================================================================!
   !==============================================================================!

44  continue

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
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE compute_dynamic_cpu(iph,dyn,freq,Opm,stat,title,normvec,keldysh_level)
    implicit none
    INTEGER,          INTENT(IN)         ::  iph
    INTEGER                              ::  Niter_,ortho
    COMPLEX(DBL),     INTENT(INOUT)      ::  dyn(:)
    TYPE(freq_type),  INTENT(IN)         ::  freq
    TYPE(eigen_type), INTENT(IN)         ::  Opm
    CHARACTER(LEN=*), INTENT(IN)         ::  stat,title
    TYPE(rcvector_type)                  ::  lastvec,initvec,tmp
    COMPLEX(8),allocatable               ::  store_vec(:,:)
    TYPE(tridiag_type)                   ::  Lanczos_matrix,tri,subLmatrix
    COMPLEX(DBL)                         ::  z
    REAL(DBL)                            ::  fermion_sign,ph_sign,norm_Opmvec2,normvec
    INTEGER                              ::  iw,iter,Niter
    INTEGER                              ::  i,ii,dimenvec
    LOGICAL                              ::  conv_one_step
    INTEGER,optional                     ::  keldysh_level
    REAL(DBL),ALLOCATABLE                ::  VECP(:,:),VALP(:)
    COMPLEX(DBL),ALLOCATABLE             ::  DD(:,:),full_mat(:,:)
   
    IF(SIZE(dyn)/=freq%Nw) STOP "ERROR IN compute_dynamic: INCONSISTENT DIMENSIONS!"

   !---------------------------------------------------------------------------!
   ! COMPUTE DYNAMIC CORRELATIONS OF OPERATORS O,O^+ I.E.  <0| O(z) * O^+ |0>  !
   !---------------------------------------------------------------------------1

    IF(iph==1) ph_sign =  one ; IF(iph==2) ph_sign = -one

    IF(iph==2.and.present(keldysh_level)) then
     write(*,*) 'ERROR keldysh in ED and should not compute hole part of the G_K function'
     stop
    endif
 
    SELECT CASE(stat)
      CASE(FERMIONIC)
        fermion_sign = -one
      CASE(BOSONIC)
                   fermion_sign =  one
        if(iph==1) fermion_sign = -one
    END SELECT

   !==============================================================================!

   if(.not.USE_TRANSPOSE_TRICK_MPI)then
    if(size(Opm%vec%rc)/=dimen_H())then
      write(*,*) '===================================='
      write(*,*) 'size of sector : ', dimen_H()
      write(*,*) 'Opm dim_space  : ', Opm%dim_space
      write(*,*) 'Opm dim vector : ', Opm%dim_sector
      write(*,*) 'size of Opm    : ', size(Opm%vec%rc)
      write(*,*) '===================================='
      STOP 'compute_dynamic inconsistent dimensions, CRITICAL'
    endif
   endif

    if(.not.USE_TRANSPOSE_TRICK_MPI)then
     dimenvec = dimen_H()
    else
     dimenvec = (  sector_h%updo%up%istatemax(iproc) -   sector_h%updo%up%istatemin(iproc)+1) *  &
                (sector_h%updo%down%istatemax(iproc) - sector_h%updo%down%istatemin(iproc)+1)
    endif

    Niter = MIN(Nitergreenmax,dimen_H())
    CALL new_tridiag(Lanczos_matrix,Niter) 
    CALL new_rcvector(lastvec,dimenvec); CALL new_rcvector(initvec,dimenvec); CALL new_rcvector(tmp,dimenvec)

    if(abs(normvec)>1.d-10)then
      norm_Opmvec2 = normvec**2
    else
      norm_Opmvec2 = norm_rcvector(Opm%vec)**2
    endif

    if(norm_Opmvec2<cutoff_dynamic) then
      dyn=0.d0; goto 44   
    endif

    if(messages3) write(*,*) 'start compute dynamical correlations by Lanczos'

    initvec%rc=Opm%vec%rc/sqrt(norm_Opmvec2)

    if(present(keldysh_level))then
       Niter_=Niter-1
       if(allocated(store_vec)) deallocate(store_vec)
       allocate(store_vec(size(initvec%rc),Niter_))
       initvec%rc=Opm%vec%rc/norm_rcvector(Opm%vec)
       store_vec(:,1)=initvec%rc
    endif

    DO iter=1,Niter
       CALL one_step_Lanczos_fast(iter,initvec%rc,tmp%rc,lastvec%rc,Lanczos_matrix,conv_one_step) 
       if(present(keldysh_level))then
          if(iter<Niter_) then
            store_vec(:,iter+1)=lastvec%rc/Lanczos_matrix%subdiag(iter+1)
            if(keldysh_ortho/=0)then
                                 ortho=keldysh_ortho
             if(keldysh_ortho<0) ortho=iter
             do i=1,ortho
              store_vec(:,iter+1)=store_vec(:,iter+1) - internal_scalprod_(conjg(store_vec(:,iter+1)),store_vec(:,i))*store_vec(:,i)
              store_vec(:,iter+1)=store_vec(:,iter+1)/norme(store_vec(:,iter+1))
             enddo
            endif
          endif
       endif
       if(conv_one_step) then
         write(*,*) 'LANCZOS CONVERGENCE AFTER ONE ITERATION'
          lastvec%rc=initvec%rc; exit
       endif
    ENDDO
   !==============================================================================!

    !--------------------------------------------------------------------------------------------------------!
    ! IF iph=1 THEN WE COMPUTE THE 'PARTICLE PART' (POLES AT w > 0) I.E. <0|O   * [ 1 / (w+E0-H) ] * O^+|0>  !
    ! IF iph=2 THEN WE COMPUTE THE 'HOLE PART'     (POLES AT w < 0) I.E. <0|O^+ * [ 1 / (w-E0+H) ] * O  |0>  !
    !                                                                                                        !
    ! THEN WE COMPUTE THE TOTAL DYNAMICAL CORRELATION                                                        !
    ! I.E. <0|O(w)*O^+(w)|0> = <0|O*[1/(w+E0-H)]*O^+|0> +/- <0|O^+*[1/(w-E0+H)]*O|0> [FERMIONIC/BOSONIC]     !
    !--------------------------------------------------------------------------------------------------------!
 
      CALL new_tridiag(tri,Lanczos_matrix)
      tri%diag = ph_sign * ( tri%diag - Opm%val ) 

      write(log_unit,*) ' Observable     = ', title
      write(log_unit,*) ' fermionic sign = ', fermion_sign
      write(log_unit,*) ' case           = ', stat
      write(log_unit,*) ' prefactor      = ', norm_Opmvec2
      write(log_unit,*) ' Keldysh        = ', present(keldysh_level)

    !-------------------------------------------------------------------------------!
   
    if(.not.present(keldysh_level))then
     if(FLAG_MPI_GREENS==0)then 
       dyn=0.d0
       call mpibarrier  
       DO iw=rank+1,freq%Nw,size2
         dyn(iw) =-fermion_sign*invert_zmtridiag(freq%vec(iw),tri)*norm_Opmvec2 
       ENDDO
       call mpisum(dyn)
     else
       dyn=0.d0
       DO iw=1,freq%Nw
         dyn(iw) =-fermion_sign*invert_zmtridiag(freq%vec(iw),tri)*norm_Opmvec2
       ENDDO
     endif
    else
        if(allocated(VECP)) deallocate(VECP,VALP,DD)
        allocate(VECP(Niter_,Niter_),VALP(Niter_),DD(Niter_,Niter_))
        CALL   submatrix_tridiag(subLmatrix,Lanczos_matrix,(/1,Niter_/))
        CALL diagonalize_tridiag(subLmatrix,VALP(1:Niter_),VECP(1:Niter_,1:Niter_),EIGENVAL_ONLY=F)
        CALL      delete_tridiag(subLmatrix)
        store_vec=matmul(store_vec,VECP)
        DD=0.d0
        do iter=1,Niter_
          DD(iter,iter)= mplx( -VALP(iter) * keldysh_delta )  
        enddo
        initvec%rc=Opm%vec%rc
        if(.true.)then !for testing...
          do iter=1,Niter_
            DD(iter,1) = DD(iter,iter) * internal_scalprod(conjg(store_vec(:,iter)),initvec%rc)
          enddo
          store_vec(:,1)=MATMUL(store_vec,DD(:,1))
        else
          if(allocated(full_mat)) deallocate(full_mat)
          allocate(full_mat(size(initvec%rc),size(initvec%rc)))
          full_mat= MATMUL( store_vec, MATMUL(DD,transpose(conjg(store_vec))) ) 
          store_vec(:,1)=MATMUL(full_mat,initvec%rc)
          deallocate(full_mat)
        endif
        if(size(Opm%vec%rc)/=size(store_vec,1))then
          write(*,*) 'error size not matching in keldysh'
          stop
        endif
#ifndef _complex
        write(*,*) 'ERROR, KELDYSH ONLY IN COMPLEX LANCZOS MODE, RECOMPILE'
        stop
#endif
        Opm%vec%rc=store_vec(:,1)
       if(allocated(store_vec)) deallocate(store_vec)
       if(allocated(VECP))      deallocate(VECP,VALP,DD)
    endif
    !-------------------------------------------------------------------------------!

    CALL delete_tridiag(tri); CALL delete_tridiag(Lanczos_matrix)

44  continue

    CALL delete_rcvector(lastvec); CALL delete_rcvector(initvec); CALL delete_rcvector(tmp)

  contains

  complex(8) function internal_scalprod_(x1,x2)
  implicit none
  complex(8) :: x1(:)
  complex(8) :: x2(:)
  integer    :: i,j,k
  internal_scalprod_=0.d0
  do i=1,size(x1)
   internal_scalprod_=internal_scalprod_+x1(i)*x2(i)
  enddo
  end function

  complex(8) function internal_scalprod(x1,x2)
  implicit none
  complex(8) :: x1(:)
#ifdef _complex
  complex(8) :: x2(:)
#else
  real(8)    :: x2(:)
#endif
  integer    :: i,j,k 
  internal_scalprod=0.d0
  do i=1,size(x1)
   internal_scalprod=internal_scalprod+x1(i)*x2(i)
  enddo
  end function

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
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE compute_dynamic_gpu(iph,dyn,freq,Opm,stat,title,normvec)
    implicit none
    INTEGER,          INTENT(IN)      ::  iph
    COMPLEX(DBL),     INTENT(INOUT)   ::  dyn(:)
    TYPE(freq_type),  INTENT(IN)      ::  freq
    TYPE(eigen_type), INTENT(IN)      ::  Opm
    CHARACTER(LEN=*), INTENT(IN)      ::  stat,title
    TYPE(tridiag_type)                ::  tri
    COMPLEX(DBL)                      ::  z
    REAL(DBL)                         ::  fermion_sign,ph_sign,norm_Opmvec2,normvec
    INTEGER                           ::  iw,iter,Niter
    INTEGER                           ::  ii,sizevec

    IF(iph==1) ph_sign = one ; IF(iph==2) ph_sign = -one

    SELECT CASE(stat)
      CASE(FERMIONIC)
        fermion_sign = -one
      CASE(BOSONIC)
                   fermion_sign =  one
        if(iph==1) fermion_sign = -one
    END SELECT

    Niter = MIN(Nitergreenmax,dimen_H())
    CALL new_tridiag(tri,Niter) 

    sizevec=dimen_H()

    if(abs(normvec)>1.d-10)then
     norm_Opmvec2 = normvec**2
    else
     norm_Opmvec2 = norm_rcvector(Opm%vec)**2
    endif

    if(norm_Opmvec2<cutoff_dynamic) then
      dyn=0.d0; goto 44
    endif

   if(messages3) write(*,*) 'start compute dynamical correlations by Lanczos GPU'
   if(.not.ON_FLY)then
    IF(ASSOCIATED(sector_h%updo))THEN
     call Lanczos_dynamic_GPU_updo(Niter,sizevec,tri%diag,tri%subdiag,Opm%vec%rc)
    else
     call Lanczos_dynamic_GPU_memory(Niter,sizevec,tri%diag,tri%subdiag,Opm%vec%rc) 
    endif
   else
    IF(ASSOCIATED(sector_h%updo))THEN
     stop 'on fly updo GPU not done'
    else
     call Lanczos_dynamic_GPU(Niter,sizevec,tri%diag,tri%subdiag,Opm%vec%rc)
    endif
   endif

    tri%diag=ph_sign * (tri%diag - Opm%val) 
    DO iw=1,freq%Nw
      dyn(iw) = - fermion_sign * invert_zmtridiag(freq%vec(iw),tri) * norm_Opmvec2 
    ENDDO
44  continue

    CALL delete_tridiag(tri)

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
!**************************************************************************
!**************************************************************************


end module
