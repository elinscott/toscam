MODULE ED_ARPACK

  use Lanczos_fast

  IMPLICIT NONE

 !-------------------------------!
 ! DIAGONALIZATION FULL SPETRUM  !
 !-------------------------------!

  REAL(DBL), PARAMETER, PRIVATE  ::  zero=0.0_DBL,one=1.0_DBL

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
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE ARPACK_diago(lowest)
    implicit none
    TYPE(eigenlist_type), INTENT(INOUT) :: lowest
#ifdef _complex
    COMPLEX(DBL),ALLOCATABLE :: VECP(:,:)
#else
    REAL(DBL),ALLOCATABLE    :: VECP(:,:)
#endif
    REAL(DBL), ALLOCATABLE              :: VALP(:)
    INTEGER                             :: start_diagH
    REAL(DBL)                           :: a
    INTEGER                             :: jj,ii,i,j
    REAL(DBL)                           :: coeff
    INTEGER                             :: dimenvec,Neigen_,nconv

    CALL reset_timer(start_diagH)
    if(dimen_H()==0) stop 'error Hilbert space has 0 dimension'
    dimenvec=dimen_H()

#ifdef _complex
    Neigen_=min(dimen_H()-1,  Neigen)
#else
    Neigen_=min(dimen_H()-1,2*Neigen)
#endif

    ALLOCATE(VECP(dimen_H(),Neigen_),VALP(Neigen_))

    !----------------------------------------------------------------------------------------------!
     if(Neigen_>0)then
     write(*,*) '======= RUNNING ARPACK WITH NEIGEN = ', Neigen_ , '========'
#ifdef _complex
     call arpack_eigenvector_sym_matrix_(.true.,'SR',tolerance,dimenvec,.true.,VALP(1:Neigen_),VECP(1:dimenvec,1:Neigen_),Neigen_,nconv,av=Hmultc)
#else
     call arpack_eigenvector_sym_matrix (.true.,'BE',tolerance,dimenvec,.true.,VALP(1:Neigen_),VECP(1:dimenvec,1:Neigen_),Neigen_,nconv,av=Hmultr)
#endif

#ifndef _complex
     if(mod(nconv,2)==0)then
      nconv=nconv/2
     else
      if(nconv/=1)then
       nconv=(nconv+1)/2
      endif
     endif
#endif

     j=1
     do i=2,nconv
      if(.not.is_eigen_in_window(VALP(i),[VALP(1),VALP(1)+dEmax]))then
       j=i
       exit 
      endif
     enddo
     nconv=j

     write(*,*) 'ARPACK converged ---> : ', nconv
     CALL timer_fortran(start_diagH,"# BUILDING OF "//TRIM(ADJUSTL(title_H_()))//" TOOK ")
     if(nconv>0)then
      CALL reset_timer(start_diagH)
      CALL add_eigen(nconv,VECP,VALP,lowest)
      CALL timer_fortran(start_diagH,"# STORING "//TRIM(ADJUSTL(title_H_()))//" TOOK ")
      if(rank==0)then
       write(*,*) '====================================================='
       write(*,*) '  GETTING THE GROUD STATE                            '
       write(*,*) '  N eigenvalues    : ', lowest%neigen
       write(*,*) '  eigenvalues      : ', lowest%eigen(:)%val
       write(*,*) '  SECTOR dimension : ', lowest%eigen(1)%dim_space
       write(*,*) '====================================================='
      endif
     endif
    else
      nconv=0
     endif
    !----------------------------------------------------------------------------------------------!

    IF(ALLOCATED(VALP)) DEALLOCATE(VALP); IF(ALLOCATED(VECP)) DEALLOCATE(VECP)
 
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

  SUBROUTINE ED_diago(lowest)
    implicit none
    TYPE(eigenlist_type), INTENT(INOUT) :: lowest
#ifdef _complex
    COMPLEX(DBL), ALLOCATABLE           :: VECP(:,:)
#else
    REAL(DBL),    ALLOCATABLE           :: VECP(:,:)
#endif
    REAL(DBL),    ALLOCATABLE           :: VALP(:)
    TYPE(rcvector_type)                 :: invec,outvec
    INTEGER                             :: start_diagH
    REAL(DBL)                           :: a
    INTEGER                             :: jj,ii,i,j
    REAL(DBL)                           :: coeff
    INTEGER                             :: dimenvec
    
    CALL reset_timer(start_diagH)
    ALLOCATE(VECP(dimen_H(),dimen_H()),VALP(dimen_H()))
    if(dimen_H()==0) stop 'error Hilbert space has 0 dimension'
    dimenvec=dimen_H()
    CALL new_rcvector(invec,dimenvec); CALL new_rcvector(outvec,dimenvec);
    invec%rc=0.d0


    !----------------------------------------------------------------------------------------------!
    DO ii=1,dimenvec
     ! H |v_i> to build H_ij !
      invec%rc(ii)=1.d0
      CALL Hmult__(dimenvec,outvec%rc,invec%rc)
      VECP(:,ii)=outvec%rc
      invec%rc(ii)=0.d0
    ENDDO
 
    CALL timer_fortran(start_diagH,"# BUILDING OF "//TRIM(ADJUSTL(title_H_()))//" TOOK ")
    CALL reset_timer(start_diagH)
    call eigenvector_matrix(lsize=dimenvec,mat=VECP,vaps=VALP,eigenvec=VECP)
    CALL timer_fortran(start_diagH,"# DIAGONALIZATION OF "//TRIM(ADJUSTL(title_H_()))//" TOOK ")
    CALL reset_timer(start_diagH)

    if(.not.FLAG_FULL_ED_GREEN)then
     j=1
     do i=2,size(VALP)
      write(*,*) 'diff val, dEmax : ', VALP(i)-VALP(1),dEmax,is_eigen_in_window(VALP(i),[VALP(1),VALP(1)+dEmax])
      if(.not.is_eigen_in_window(VALP(i),[VALP(1),VALP(1)+dEmax]))then
       j=i
       exit
      endif
     enddo
    else
     j=dimenvec
    endif

    CALL add_eigen(j,VECP,VALP,lowest)
    !----------------------------------------------------------------------------------------------!

    IF(ALLOCATED(VALP)) DEALLOCATE(VALP); IF(ALLOCATED(VECP)) DEALLOCATE(VECP) 

    CALL delete_rcvector(invec) 
    CALL delete_rcvector(outvec)
    CALL timer_fortran(start_diagH,"# STORING "//TRIM(ADJUSTL(title_H_()))//" TOOK ")

    if(rank==0)then
      write(*,*) '====================================================='
      write(*,*) '  GETTING THE GROUD STATE                            '
      write(*,*) '  N eigenvalues    : ', lowest%neigen
      write(*,*) '  dEmax            : ', dEmax
      write(*,*) '  j                : ', j
      write(*,*) '  eigenvalues      : ', lowest%eigen(:)%val
      write(*,*) '  SECTOR dimension : ', lowest%eigen(1)%dim_space
      write(*,*) '====================================================='
    endif

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
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

END MODULE 
