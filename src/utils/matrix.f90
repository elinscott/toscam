module matrix

   use genvar
   use linalg, only : matrixinverse, same_array, erase_divergence, q_zsum, csign1, &
       zdot, zaxpy__, q_zamax, q_zscal, swap, outerprod, cabs1, zsum, qsign, zaxpyb
   use random
   use tools_algebra
   use sorting
   use common_def,   only : reset_timer, timer_fortran, dump_message, c2s, i2c, &
       create_seg_fault
   ! use fortran_cuda

private
public :: average_matrix
public :: average_vec
public :: bande_mat
public :: bande_matc
public :: decomposevec
public :: diag
public :: diagc
public :: diagi
public :: diagonalize
public :: diagonalize_real
public :: diagr
public :: eigenvector_matrix
public :: eigenvector_matrix_c_
public :: Id
public :: invmat
public :: matmul_x
public :: matmul_x_c 
public :: new_diag
public :: new_Id
public :: offdiag
public :: offdiagc
public :: qr_decomp
public :: rearrange_columns_to_identity
public :: rearrange_columns_to_identity_c
public :: write_array
public :: write_real_array_rank_2
public :: write_cplx_array_rank_2

!--------------------------------------------------------------------------------!
INTERFACE rearrange_columns_to_identity
 MODULE PROCEDURE rearrange_columns_to_identity_r,rearrange_columns_to_identity_c
END INTERFACE
!--------------------------------------------------------------------------------!
! INTERFACE invmat_sym
!  MODULE PROCEDURE invmat_sym_c,invmat_sym_r
! END INTERFACE
!--------------------------------------------------------------------------------!
INTERFACE bande_mat
 MODULE PROCEDURE bande_matr,bande_matc 
END INTERFACE 
!--------------------------------------------------------------------------------!
! INTERFACE getrf__
!  MODULE PROCEDURE getrfc__,getrfcs__,getrfr__,getrfs__
! END INTERFACE
! !--------------------------------------------------------------------------------!
! INTERFACE getri__
!  MODULE PROCEDURE getric__,getrics__,getrir__,getrirs__
! END INTERFACE
! !--------------------------------------------------------------------------------!
! INTERFACE gedi__
!  MODULE PROCEDURE gedic__,gedicq__
! END INTERFACE
! !--------------------------------------------------------------------------------!
! INTERFACE geco__
!  MODULE PROCEDURE gecoc__,gecocs__,gecocq__
! END INTERFACE
!--------------------------------------------------------------------------------!
! INTERFACE svd_wrapper
!     MODULE PROCEDURE svd_wrapper_cs,svd_wrapper_c,svd_wrapper_r,svd_wrapper_rs,svd_wrapper_rq
! END INTERFACE
! !--------------------------------------------------------------------------------!
INTERFACE MATMUL_
    MODULE PROCEDURE MATMULr,MATMULc
END INTERFACE
! !--------------------------------------------------------------------------------!
! INTERFACE tqli
!     MODULE PROCEDURE tqli_,tqli__,tqliq
! END INTERFACE
!--------------------------------------------------------------------------------!
interface invmat
   module procedure invert_lapack_real
   module procedure invert_lapack_complex
end interface
!--------------------------------------------------------------------------------!
! INTERFACE invmat
!     MODULE PROCEDURE invmat_comp,invmat_comps,invmat_real,invmat_reals,invmat_real_quad,invmat_comp_quad
! END INTERFACE
! !--------------------------------------------------------------------------------!
! INTERFACE invmat_jordan
!     MODULE PROCEDURE inverse_,inverse__,inverses_,inverses__
! END INTERFACE
!--------------------------------------------------------------------------------!
! INTERFACE ludcmp
!     MODULE PROCEDURE ludcmpd,ludcmpq
! END INTERFACE
!--------------------------------------------------------------------------------!
! INTERFACE lubksb
!     MODULE PROCEDURE lubksbd,lubksbq
! END INTERFACE
!--------------------------------------------------------------------------------!
! INTERFACE MATMUL_keep_diag
!     MODULE PROCEDURE MATMUL_keep_diag_,MATMUL_keep_diag__
! END INTERFACE
! !--------------------------------------------------------------------------------!
! INTERFACE MATMUL_sum_diag
!     MODULE PROCEDURE MATMUL_sum_diag_,MATMUL_sum_diag__
! END INTERFACE
! !--------------------------------------------------------------------------------!
! INTERFACE symmetrize_mat
!   MODULE PROCEDURE symmetrize_mat_c,symmetrize_mat_r
! END INTERFACE
!--------------------------------------------------------------------------------!
! INTERFACE rescale
!  MODULE PROCEDURE rescale_r,rescale_c,rescale_r2,rescale_c2,rescale_r3,rescale_c3
! END INTERFACE
!--------------------------------------------------------------------------------!
INTERFACE eigenvector_matrix
 MODULE PROCEDURE      eigenvector_matrix_cc,eigenvector_matrix_c,eigenvector_matrix_r,    &
                     & eigenvector_matrix_rr,eigenvector_matrix_cc_,eigenvector_matrix_c_, &
                     & eigenvector_matrix_rc,eigenvector_matrix_rrc
END INTERFACE
!--------------------------------------------------------------------------------!
! INTERFACE eigenvector_matrix_b_
!  MODULE PROCEDURE eigenvector_matrixa_,eigenvector_matrixc_
! END INTERFACE
! !--------------------------------------------------------------------------------!
! INTERFACE eigenvalue_matrix
!  MODULE PROCEDURE eigenvalue_matrix__,eigenvalue_matrix_,eigenvalue_matrixr_
! END INTERFACE
! !--------------------------------------------------------------------------------!
! INTERFACE print_eigenvalue
!  MODULE PROCEDURE print_eigenvalue__
! END INTERFACE
!--------------------------------------------------------------------------------!
INTERFACE MATMUL_x
  MODULE PROCEDURE MATMUL_x_c,MATMUL_x_r
END INTERFACE
!--------------------------------------------------------------------------------!
! INTERFACE get_det_from_zgeco
!  MODULE PROCEDURE get_det_from_zgeco_,get_det_from_zgeco__,get_det_from_zgeco___
! END INTERFACE
!--------------------------------------------------------------------------------!
! INTERFACE INVMG
!  MODULE PROCEDURE D_INVMG,Q_INVMG
! END INTERFACE
!--------------------------------------------------------------------------------!
INTERFACE diag
 MODULE PROCEDURE diagr,diagc,diagi,diagrr,diagr_,diagc_,diagc__,diagr__
END INTERFACE
!--------------------------------------------------------------------------------!
INTERFACE offdiag
 MODULE PROCEDURE offdiagr,offdiagc
END INTERFACE
!--------------------------------------------------------------------------------!
INTERFACE write_array
    MODULE PROCEDURE write_bool_array_rank_3  !  write A(n1,n2,n3) boolean
    MODULE PROCEDURE write_intg_array_rank_3  !  write A(n1,n2,n3) integer
    MODULE PROCEDURE write_real_array_rank_3  !  write A(n1,n2,n3) real
    MODULE PROCEDURE write_cplx_array_rank_3  !  write A(n1,n2,n3) complex
    MODULE PROCEDURE write_bool_array_rank_2  !  write A(n1,n2) boolean
    MODULE PROCEDURE write_intg_array_rank_2  !  write A(n1,n2) integer
    MODULE PROCEDURE write_real_array_rank_2  !  write A(n1,n2) real
    MODULE PROCEDURE write_cplx_array_rank_2  !  write A(n1,n2) complex
    MODULE PROCEDURE write_bool_array_rank_1  !  write A(n1) boolean
    MODULE PROCEDURE write_intg_array_rank_1  !  write A(n1) integer
    MODULE PROCEDURE write_real_array_rank_1  !  write A(n1) real
    MODULE PROCEDURE write_cplx_array_rank_1  !  write A(n1) complex
END INTERFACE
!--------------------------------------------------------------------------------!
INTERFACE diagonalize
    MODULE PROCEDURE diagonalize_comp,diagonalize_real
END INTERFACE
!--------------------------------------------------------------------------------!
! INTERFACE invert_lapack
!     MODULE PROCEDURE invert_lapack_r,invert_lapack_C
! END INTERFACE
!--------------------------------------------------------------------------------!
! INTERFACE invert_openmp
!     MODULE PROCEDURE invert_pivot_double_,invert_pivot_complex_,invert_pivot_single_,invert_pivot_complexs_
! END INTERFACE

!--------------------------------------------------------------------------------!
!  real(8),parameter,private :: rerror=1.d-12,invmat_error=1.d-13,rrerror=1.d-28 
!  real(8),parameter,private :: smallest_pivot=1.d-12,zero=0.d0,one=1.d0
!  integer,parameter,private :: n_cuda_rout=100,n_openmp_rout=30
!--------------------------------------------------------------------------------!

contains

   subroutine invert_lapack_real(n, A)
   
      implicit none
   
      integer,           intent(in   ) :: n
      real(8),           intent(inout) :: A(:,:)
   
      integer, allocatable :: piv(:)
      real(8), allocatable :: WORK(:)
      integer              :: ierr
   
      external :: dgetrf, dgetri
   
      ! Check we have correct dimensions
      if (size(A, dim=1) /= n .or. size(A, dim=2) /= n) then
         write(*,*) "Error in invert_lapack_real: matrix dimensions do not match input"
         stop
      end if

      ! Allocate local variables
      allocate(piv(n), stat = ierr)
      allocate(WORK(n), stat = ierr)
   
      ! Perform inversion
      call dgetrf(n, n, A, n, piv, ierr)
      call dgetri(n, A, n, piv, WORK, n, ierr)
   
      ! Check error flag
      if (ierr /= 0) then
         write(*,*) 'Error in invert_lapack_real: matrix has no inverse'
      endif
   
      ! Deallocate local variables
      deallocate(piv, stat=ierr)
      deallocate(WORK, stat=ierr)
   
   end subroutine
   
   subroutine invert_lapack_complex(n, A)
   
      implicit none
   
      integer,           intent(in   ) :: n
      complex(8),        intent(inout) :: A(:,:)
   
      integer, allocatable    :: piv(:)
      complex(8), allocatable :: WORK(:)
      integer                 :: ierr
   
      external :: zgetrf, zgetri
   
      ! Check we have correct dimensions
      if (size(A, dim=1) /= n .or. size(A, dim=2) /= n) then
         write(*,*) "Error in invert_lapack_complex: matrix dimensions do not match input"
         stop
      end if
   
      ! Allocate local variables
      allocate(piv(n), stat = ierr)
      allocate(WORK(n), stat = ierr)
   
      ! Perform inversion
      call zgetrf(n, n, A, n, piv, ierr)
      call zgetri(n, A, n, piv, WORK, n, ierr)
   
      ! Check error flag
      if (ierr /= 0) then
         write(*,*) 'Error in invert_lapack_complex: matrix has no inverse'
      endif
   
      ! Deallocate local variables
      deallocate(piv, stat=ierr)
      deallocate(WORK, stat=ierr)
   
   end subroutine

! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!          !---------------------!
! 
!   subroutine invmat_sym_c(cmatrix,num)
!     implicit none
!     integer, intent(in) :: num
!     complex(kind=DP), intent(inout) :: cmatrix(num,num) ! inverse on exit
!     integer :: work_length, info, row, col
!     integer, allocatable, dimension(:) :: ipiv
!     complex(kind=DP), allocatable, dimension(:) :: work_array
!     integer :: ierr ! error flag
! 
!     work_length=3*num
!     allocate(ipiv(num),stat=ierr)
!     allocate(work_array(work_length),stat=ierr)
! 
!     call zsytrf('L', num, cmatrix, num, ipiv, work_array, work_length, info)
!     if (info.ne.0) then
!              write(6,'(a,i5)') 'ERROR in &
!             &wrappers_invert_sym_cmatrix: Problem with zsytrf, info=',info
!        stop
!     endif
! 
!     call zsytri('L', num, cmatrix, num, ipiv, work_array, info )
!     if (info.ne.0) then
!            write(6,'(a,i5)') 'ERROR in &
!             &wrappers_invert_sym_cmatrix: Problem with zsytri, info=',info
!        stop
!     endif
! 
!     do row=1,num
!        do col=1,row-1
!           cmatrix(col,row)=cmatrix(row,col)
!        enddo
!     enddo
! 
!     deallocate(work_array,stat=ierr)
!     deallocate(ipiv,stat=ierr)
! 
!   end subroutine
! 
!          !---------------------!
! 
!   subroutine invmat_sym_r(matrix,num)
!     implicit none
!     integer, intent(in) :: num
!     real(8), intent(inout) :: matrix(num,num) ! inverse on exit
!     integer :: work_length, info, row, col
!     integer, allocatable, dimension(:) :: ipiv
!     real(8), allocatable, dimension(:) :: work_array
!     integer :: ierr ! error flag
! 
! 
!     work_length=3*num
!     allocate(ipiv(num),stat=ierr)
!     allocate(work_array(work_length),stat=ierr)
!     call dsytrf('L', num, matrix, num, ipiv, work_array, work_length, info)
!     if (info.ne.0) then
!          write(6,'(a,i5)') 'ERROR in &
!             &wrappers_invert_sym_matrix: Problem with dsytrf, info=',info
!        stop
!     endif
!     call dsytri('L', num, matrix, num, ipiv, work_array, info )
!     if (info.ne.0) then
!             write(6,'(a,i5)') 'ERROR in &
!             &wrappers_invert_sym_matrix: Problem with dsytri, info=',info
!        stop
!     endif
!     do row=1,num
!        do col=1,row-1
!           matrix(col,row)=matrix(row,col)
!        enddo
!     enddo
!     deallocate(work_array,stat=ierr)
!     deallocate(ipiv,stat=ierr)
!   end subroutine 
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!        !-----------------------------!
! 
!  subroutine map_vector_permutation(vec,mapping)
!  implicit none
!  real(8) :: vec(:),temp(size(vec))
!  integer :: mapping(:),i,ii
!     temp=vec;vec=0
!     do i=1,size(vec,1)
!        ii=mapping(i); vec(i)=temp(ii)
!     enddo
!  end subroutine
! 
!        !-----------------------------!
! 
!  subroutine map_matrix_permutation(mat,mapping)
!  implicit none
!  complex(8) :: mat(:,:),temp(size(mat,1),size(mat,2))
!  integer    :: mapping(:),i,j,k,l,ii,jj
!      temp=mat;mat=0
!      do i=1,size(mat,1); do j=1,size(mat,2)
!          ii=mapping(i); jj=mapping(j); mat(i,j)=temp(ii,jj)
!      enddo; enddo
!  end subroutine
! 
!        !-----------------------------!
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
! 
 subroutine average_vec(mat,MASK_AVERAGE_)
 implicit none
 real(8)    :: mat(:)
 integer    :: MASK_AVERAGE_(:),MASK_AVERAGE(size(mat,1))
 integer    :: i,l,m,tot,i_
 real(8)    :: temp

  MASK_AVERAGE=MASK_AVERAGE_

  do i=1,size(mat,1)
    if(abs(MASK_AVERAGE(i))>0)then
     m=abs(MASK_AVERAGE(i));temp=0.d0;tot=0
      do i_=1,size(mat,1)
       l=MASK_AVERAGE(i_)
       if(abs(l)==m) then
        temp=temp+dble(SIGN(1,l))*mat(i_); tot=tot+1
       endif
      enddo
      temp=temp/dble(tot)
      do i_=1,size(mat,1)
       l=MASK_AVERAGE(i_)
       if(abs(l)==m) then
        mat(i_)=dble(SIGN(1,l))*temp
        MASK_AVERAGE(i_)=0
       endif
      enddo
    endif
    if(MASK_AVERAGE_(i)==0) mat(i)=0.d0
  enddo

 return
 end subroutine
! 
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
! 
!  subroutine average_matrix_by_block(mat,M)
!  implicit none
!   complex(8) :: mat(:,:)
!   integer    :: M(:,:)
!   integer    :: i,j,siz1,siz2,siz1s,siz2s
! 
!    siz1 =size(mat,1) ; siz2=size(mat,2) ; siz1s=size(mat,1)/2 ; siz2s=size(mat,2)/2 ;
!    call average_matrix(mat(1:siz1s,1:siz2s),M(1:siz1s,1:siz2s))
!    call average_matrix(mat(siz1s+1:siz1,siz2s+1:siz2),M(siz1s+1:siz1,siz2s+1:siz2))
!    call average_matrix(mat(1:siz1s,siz2s+1:siz2),M(1:siz1s,siz2s+1:siz2))
!    call average_matrix(mat(siz1s+1:siz1,1:siz2s),M(siz1s+1:siz1,1:siz2s))
! 
!  return
!  end subroutine
! 
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
! 
 subroutine average_matrix(mat,MASK_AVERAGE_,offdiag_also)
 implicit none
 complex(8)        :: mat(:,:)
 integer,intent(in):: MASK_AVERAGE_(:,:)
 integer           :: MASK_AVERAGE(size(mat,1),size(mat,2))
 integer           :: i,j,l,m,tot,i_,j_
 complex(8)        :: temp
 logical,optional  :: offdiag_also
 logical           :: only_diag

  only_diag=.false.

  if(present(offdiag_also))then
   only_diag=.not.offdiag_also
  endif 

  MASK_AVERAGE=MASK_AVERAGE_
  
  do i=1,size(mat,1)
   do j=1,size(mat,2)
   if(.not.only_diag.or.j==i)then

    if(abs(MASK_AVERAGE(i,j))>0)then
      m=abs(MASK_AVERAGE(i,j));temp=0.d0;tot=0
      do i_=1,size(mat,1)
       do j_=1,size(mat,2)
        l=MASK_AVERAGE(i_,j_)
        if(abs(l)==m)then
         temp=temp+dble(SIGN(1,l))*mat(i_,j_)
         tot=tot+1
        endif
       enddo
      enddo
      temp=temp/dble(tot)
      do i_=1,size(mat,1)
       do j_=1,size(mat,2)
        l=MASK_AVERAGE(i_,j_)
        if(abs(l)==m)then
          mat(i_,j_)=dble(SIGN(1,l))*temp
          MASK_AVERAGE(i_,j_)=0 
        endif
       enddo
      enddo
    endif

    if(MASK_AVERAGE_(i,j)==0)  mat(i,j)=0.d0

   endif
   enddo
  enddo 

 return
 end subroutine
! 
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
! 
!  subroutine average_diag_by_block(MASK,vec,vec2)
!  implicit none
!  real(8) :: vec(:),vec2(:),vectemp(size(vec))
!  integer :: MASK(:)
!  integer :: siz
! 
!   siz=size(vec)
!   if(mod(siz,2)/=0)then
!    write(*,*) 'average by block, size does not match'
!    stop 'critical'
!   endif
!   vectemp(1:siz/2)=vec(1:siz/2)
!   vectemp(siz/2+1:siz)=vec2(siz/2+1:siz)
!   call average_vec(vectemp,MASK)
!   vec=vectemp
! 
!  return
!  end subroutine
! 
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
!          !---------------------!
! 
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

!     !-----------------------!
! 
!  subroutine getrfs__(mat,piv)
!  implicit none
!  integer       :: nnn,piv(:),INFO
!  real(4)       :: mat(:,:)
!     nnn=size(mat,1)
!     ! ebl: Removing GPU functionality
!     ! if(use_cula_routines)then
!     !  call cula__getrf(nnn,mat,piv)
!     !  return
!     ! endif
!     CALL SGETRF(nnn,nnn,mat,nnn,piv,INFO)
!     IF(info>0)THEN
!       CALL dump_message(TEXT="ERROR IN 'invert[GETRF]': MATRIX IS SINGULAR : ")
!       STOP 'error invert_lapack singulat matrix , getrf'
!     ENDIF
!  end subroutine
! 
!     !-----------------------!
! 
!  subroutine getrfcs__(mat,piv)
!  implicit none
!  integer       :: nnn,piv(:),INFO
!  complex(4)    :: mat(:,:)
!     nnn=size(mat,1)
!     ! ebl: Removing GPU functionality
!     ! if(use_cula_routines)then
!     !   call cula__getrf(nnn,mat,piv)
!     !   return
!     ! endif
!     CALL CGETRF(nnn,nnn,mat,nnn,piv,INFO)
!     IF(info>0)THEN
!       CALL dump_message(TEXT="ERROR IN 'invert[GETRF]': MATRIX IS SINGULAR : ")
!       STOP 'error invert_lapack singulat matrix , getrf'
!     ENDIF
!  end subroutine
! 
!     !-----------------------!
! 
!  subroutine getrfc__(mat,piv)
!  implicit none
!  integer       :: nnn,piv(:),INFO
!  complex(8)    :: mat(:,:)
!     nnn=size(mat,1)
!     ! ebl: Removing GPU functionality
!     ! if(use_cula_routines)then
!     !   call cula__getrf(nnn,mat,piv)
!     !   return
!     ! endif
!     CALL ZGETRF(nnn,nnn,mat,nnn,piv,INFO)
!     IF(info>0)THEN
!       CALL dump_message(TEXT="ERROR IN 'invert[GETRF]': MATRIX IS SINGULAR : ")
!       CALL write_array(mat,'mat')
!       STOP 'error invert_lapack singulat matrix , getrf'
!     ENDIF
!  end subroutine
! 
!     !-----------------------!
! 
!  subroutine getrfr__(mat,piv)
!  implicit none
!  integer    :: nnn,piv(:),INFO
!  real(8)    :: mat(:,:)
!     nnn=size(mat,1)
!     ! ebl: Removing GPU functionality
!     ! if(use_cula_routines)then
!     !  call cula__getrf(nnn,mat,piv)
!     !  return
!     ! endif
!     CALL DGETRF(nnn,nnn,mat,nnn,piv,INFO) 
!     IF(info>0)THEN
!       CALL dump_message(TEXT="ERROR IN 'invert[GETRF]': MATRIX IS SINGULAR : ")
!       CALL write_array(mat,'mat')
!       STOP 'error invert_lapack singulat matrix , getrf'
!     ENDIF
!  end subroutine
! 
!     !-----------------------!
! 
!  subroutine getrirs__(mat,piv)
!  implicit none
!  integer    :: nnn,piv(:),INFO
!  real(4)    :: mat(:,:),WORK(size(mat,1))
!     nnn=size(mat,1)
!     ! ebl: Removing GPU functionality
!     ! if(use_cula_routines)then
!     !  call cula__getri(nnn,mat,piv)
!     !  return
!     ! endif
!     call SGETRI(nnn,mat,nnn,piv,WORK,nnn,INFO)
!     IF(info>0)THEN
!       CALL dump_message(TEXT="ERROR IN 'invert[GETRI]': MATRIX IS SINGULAR : ")
!       STOP 'error invert_lapack singulat matrix , getri'
!     ENDIF
!  end subroutine
! 
!     !-----------------------!
! 
!  subroutine getrir__(mat,piv)
!  implicit none
!  integer    :: nnn,piv(:),INFO
!  real(8)    :: mat(:,:),WORK(size(mat,1))
!     nnn=size(mat,1)
!     ! ebl: Removing GPU functionality
!     ! if(use_cula_routines)then
!     !  call cula__getri(nnn,mat,piv)
!     !  return
!     ! endif
!     call DGETRI(nnn,mat,nnn,piv,WORK,nnn,INFO)
!     IF(info>0)THEN
!       CALL dump_message(TEXT="ERROR IN 'invert[GETRI]': MATRIX IS SINGULAR : ")
!       CALL write_array(mat,'mat')
!       STOP 'error invert_lapack singulat matrix , getri'
!     ENDIF
!  end subroutine
! 
!     !-----------------------!
! 
!  subroutine getrics__(mat,piv)
!  implicit none
!  integer    :: nnn,piv(:),INFO
!  complex(4) :: mat(:,:),WORK(size(mat,1))
!     nnn=size(mat,1)
!     ! ebl: Removing GPU functionality
!     ! if(use_cula_routines)then
!     !  call cula__getri(nnn,mat,piv)
!     !  return
!     ! endif
!     call CGETRI(nnn,mat,nnn,piv,WORK,nnn,INFO)
!     IF(info>0)THEN
!       CALL dump_message(TEXT="ERROR IN 'invert[GETRI]': MATRIX IS SINGULAR : ")
!       STOP 'error invert_lapack singulat matrix , getri'
!     ENDIF
!  end subroutine
! 
!     !-----------------------!
! 
!  subroutine getric__(mat,piv)
!  implicit none
!  integer    :: nnn,piv(:),INFO
!  complex(8) :: mat(:,:),WORK(size(mat,1))
!     nnn=size(mat,1)
!     ! ebl: Removing GPU functionality
!     ! if(use_cula_routines)then
!     !  call cula__getri(nnn,mat,piv)
!     !  return
!     ! endif
!     call ZGETRI(nnn,mat,nnn,piv,WORK,nnn,INFO)
!     IF(info>0)THEN
!       CALL dump_message(TEXT="ERROR IN 'invert[GETRI]': MATRIX IS SINGULAR : ")
!       CALL write_array(mat,'mat')
!       STOP 'error invert_lapack singulat matrix , getri'
!     ENDIF
!  end subroutine
! 
!     !-----------------------!
! 
!  subroutine gecoc__(nnn,mat,piv)
!  implicit none
!  integer      :: nnn
!  real(8)      :: rcond
!  integer      :: piv(nnn)
!  complex(8)   :: mat(nnn,nnn),WORK2(size(mat,1))
!    if(.not.fast_invmat)then
!      call zgeco(mat,nnn,nnn,piv,rcond,WORK2)
!     else
!      call zgeco_(mat,nnn,nnn,piv,rcond,WORK2)
!    endif
!  end subroutine
! 
!   !-----------------------!
! 
!  subroutine gedic__(nnn,mat,piv,deti)
!  implicit none
!  integer      :: nnn
!  real(8)      :: rcond
!  integer      :: piv(nnn)
!  complex(8)   :: mat(nnn,nnn),WORK(2*size(mat,1)),deti(2)
!    if(.not.fast_invmat)then
!      call zgedi(mat,nnn,nnn,piv,deti,WORK,11)
!     else
!      call zgedi_(mat,nnn,nnn,piv,deti,WORK,11)
!    endif
!  end subroutine
! 
!     !-----------------------!
! 
!  subroutine gecocs__(nnn,mat,piv)
!  implicit none
!  integer      :: nnn
!  real(4)      :: rcond
!  integer      :: piv(nnn)
!  complex(4)   :: mat(nnn,nnn),WORK2(size(mat,1))
!    call cgeco(mat,nnn,nnn,piv,rcond,WORK2)
!  end subroutine
! 
!     !-----------------------!
! 
!  subroutine gedics__(nnn,mat,piv,deti)
!  implicit none
!  integer      :: nnn
!  real(4)      :: rcond
!  integer      :: piv(nnn)
!  complex(4)   :: mat(nnn,nnn),WORK(2*size(mat,1)),deti(2)
!    call cgedi(mat,nnn,nnn,piv,deti,WORK,11)
!  end subroutine
! 
!     !-----------------------!
! 
!  subroutine gecocq__(nnn,mat,piv)
!  implicit none
!  integer      :: nnn
!  real(16)     :: rcond
!  integer      :: piv(nnn)
!  complex(16)  :: mat(nnn,nnn),WORK2(size(mat,1))
!    call q_zgeco(mat,nnn,nnn,piv,rcond,WORK2)
!  end subroutine
! 
!     !-----------------------!
! 
!  subroutine gedicq__(nnn,mat,piv,deti)
!  implicit none
!  integer      :: nnn
!  real(16)     :: rcond
!  integer      :: piv(nnn)
!  complex(16)  :: mat(nnn,nnn),WORK(2*size(mat,1)),deti(2)
!      call q_zgedi(mat,nnn,nnn,piv,deti,WORK,11)
!  end subroutine
! 
!     !-----------------------!

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
! 
! !---------!
! 
! subroutine svd_wrapper_cs(A,values,V,M,N,method)
! implicit none
! integer :: M,N
! complex(4) :: A(M,N),A_BACK(M,N),U(M,M),V(N,N)
! real(4)    :: values(N),E(N),work(M)
! integer    :: method,ierr
!   A_BACK=A 
!   if(M<N) then
!    stop 'svd_wrapper, leading dim of matrix should be rows'
!   endif
!   SELECT CASE(method) 
!   CASE(5)
!     call cula_svd(M,N,A,values,U,V)
!   CASE DEFAULT
!     stop 'svd_not_defined'
!   END SELECT
!   A=A_BACK
! end subroutine
! 
! !---------!
! 
! subroutine svd_wrapper_c(A,values,V,M,N,method)
! implicit none
! integer :: M,N
! complex(8) :: A(M,N),A_BACK(M,N),U(M,M),V(N,N)
! real(8)    :: values(N),E(N),work(M)
! integer    :: method,ierr
!   A_BACK=A
!   if(M<N) then
!    stop 'svd_wrapper, leading dim of matrix should be rows'
!   endif
!   SELECT CASE(method)
!   CASE(5)
!     call cula_svd(M,N,A,values,U,V)
!   CASE DEFAULT
!     stop 'svd_not_defined'
!   END SELECT
!   A=A_BACK
! end subroutine
! 
! !---------!
! 
! subroutine svd_wrapper_r(A,values,V,M,N,method)
! implicit none
! integer :: M,N
! real(8) :: A(M,N),A_BACK(M,N),U(M,M),V(N,N)
! real(8) :: values(N),E(N),work(M)
! integer :: method,ierr
!   A_BACK=A
!   if(M<N) then
!    stop 'svd_wrapper, leading dim of matrix should be rows'
!   endif
!   SELECT CASE(method)
!   CASE(1)
!     call SVDCMP(A,M,N,M,N,values,V)
!   CASE(2)
!     call SVD(A,values,.true.,U,.true.,V,ierr)
!   CASE(3)
! #ifndef NOSVD
!     call svdcmp_dp(A,values,V)  
! #else
!     write(*,*) 'no svd stop'
!     stop
! #endif
!   CASE(4)
!     call dsvdc(A,M,M,N,values,E,U,M,V,N,work,11,ierr)
!   CASE(5)
!     call cula_svd(M,N,A,values,U,V)
!   CASE DEFAULT
!     stop 'svd_not_defined'
!   END SELECT
!   A=A_BACK
! end subroutine
! 
! !---------!
! 
! subroutine svd_wrapper_rs(A,values,V,M,N,method)
! implicit none
! integer :: M,N
! real(4) :: A(M,N),A_BACK(M,N),U(M,M),V(N,N)
! real(4) :: values(N),E(N),work(M)
! integer :: method,ierr
!   A_BACK=A
!   if(M<N) then
!    stop 'svd_wrapper, leading dim of matrix should be rows'
!   endif
!   SELECT CASE(method)
!   CASE(2)
!     call SVD(A,values,.true.,U,.true.,V,ierr)
!   CASE(3)
! #ifndef NOSVD
!     call svdcmp_sp(A,values,V)
! #else
!     write(*,*) 'no svd stop'
!     stop
! #endif
!   CASE(5)
!     call cula_svd(M,N,A,values,U,V)
!   CASE DEFAULT
!     stop 'svd_not_defined'
!   END SELECT
!   A=A_BACK
! end subroutine
! 
! !---------!
! 
! subroutine svd_wrapper_rq(A,values,V,M,N,method)
! implicit none
! integer  :: M,N
! real(16) :: A(M,N),A_BACK(M,N),U(M,M),V(N,N)
! real(16) :: values(N),E(N),work(M)
! integer  :: method,ierr
!   A_BACK=A
!   if(M<N) then
!    stop 'svd_wrapper, leading dim of matrix should be rows'
!   endif
!   SELECT CASE(method)
!   CASE(1)
!     call Q_SVDCMP(A,M,N,M,N,values,V)
!   CASE(2)
!     call SVD(A,values,.true.,U,.true.,V,ierr) 
!   CASE DEFAULT
!     stop 'svd_not_defined'
!   END SELECT
!   A=A_BACK
! end subroutine
! 
! !---------!
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! subroutine test_cuda_eigenvectors
! implicit none
! integer :: i
! 
!  write(*,*) 'double precision....'
!  do i=100,800,100
!    call testit(i)
!  enddo
!  write(*,*) 'single precision....'
!  do i=100,800,100
!    call testit_(i)
!  enddo
!  stop 'testing done'
! 
! contains 
! 
! !---------!
! !---------!
! !---------!
! 
! subroutine testit(nnn)
! implicit none
! integer              :: nnn 
! real(8)              :: A(nnn,nnn),B(nnn,nnn),C(nnn,nnn),vec(nnn),vec2(nnn)
! integer              :: i,j,n,jjj,jj
! use_cuda_routines=.false.
! A=0.d0
! do i=1,nnn
!  A(i,i)=drand1()
!  do j=i+1,nnn
!    A(i,j)=drand1()
!    A(j,i)=A(i,j)
!  enddo
! enddo
! call reset_timer(jjj)
! write(*,*) 'start cpu calculations'
! call eigenvector_matrix(lsize=nnn,mat=A,vaps=vec,eigenvec=B)
! write(*,*) 'done'
! call timer_fortran(jjj,'CPU TOOK : ', unit_=6)
! use_cuda_routines=.true.
! call reset_timer(jjj)
! write(*,*) 'start gpu calculations'
! call eigenvector_matrix(lsize=nnn,mat=A,vaps=vec2,eigenvec=C)
! write(*,*) 'done'
! call timer_fortran(jjj,'GPU TOOK : ', unit_=6)
! use_cuda_routines=.false.
! write(*,*) 'maxval gpu eigen  : ', maxval(abs(C))
! write(*,*) 'maxval gpu values : ', maxval(abs(vec2))
! write(*,*) 'checking gpu : '
! do i=1,1
!  write(*,*) 'checking vector : ', i
!  write(*,*) 'A vec - lambda*vec : ', maxval(abs(matmul(A,C(:,i))-vec2(i)*C(:,i)))
!  write(*,*) 'norme vec, eigenvalue : ', norme(C(:,i)),vec2(i)
! enddo
! end subroutine
! 
! !---------!
! !---------!
! !---------!
! 
! subroutine testit_(nnn)
! implicit none
! integer              :: nnn
! real(4)              :: A(nnn,nnn),B(nnn,nnn),C(nnn,nnn),vec(nnn),vec2(nnn)
! integer              :: i,j,n,jjj,jj
! A=0.d0
! do i=1,nnn
!  A(i,i)=drand1()
!  do j=i+1,nnn
!    A(i,j)=drand1()
!    A(j,i)=A(i,j)
!  enddo
! enddo
! call reset_timer(jjj)
! write(*,*) 'start cpu calculations'
! use_cuda_routines=.false.
! call eigenvector_matrix_rr(lsize=nnn,mat=A,vaps=vec,eigenvec=B)
! write(*,*) 'done'
! call timer_fortran(jjj,'CPU TOOK : ', unit_=6)
! use_cuda_routines=.true.
! call reset_timer(jjj)
! write(*,*) 'start gpu calculations'
! call eigenvector_matrix_rr(lsize=nnn,mat=A,vaps=vec2,eigenvec=C)
! write(*,*) 'done'
! call timer_fortran(jjj,'GPU TOOK : ', unit_=6)
! use_cuda_routines=.false.
! write(*,*) 'eig gpu : ', minval(vec2),maxval(vec2)
! write(*,*) 'eig cpu : ', minval(vec),maxval(vec)
! write(*,*) 'checking gpu : '
! do i=1,1
!  write(*,*) 'checking vector : ', i
!  write(*,*) 'A vec - lambda*vec : ', maxval(abs(matmul(A,C(:,i))-vec2(i)*C(:,i)))
!  write(*,*) 'norme eig gpu : ', vec2(i)
!  write(*,*) 'norme eig cpu : ', vec(i)
! enddo
! end subroutine
! 
! !---------!
! !---------!
! !---------!
! 
! end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! subroutine test_invmat_array_matrix_collect
! implicit none
! integer,parameter    :: nnn=20,nf=2000
! complex(8)           :: A(nnn,nnn),B(nnn,nnn,nf),C(nnn,nnn),vec(nf)
! real(8)              :: c1(nnn,nnn,nf),c2(nnn,nnn,nf)
! integer              :: i,j,n,jjj,jj
! 
! do i=1,nnn
! do j=1,nnn
!  A(i,j)=drand1()
! enddo
! enddo
! do i=1,nf
!  vec(i)=drand1()+imi*drand1()
! enddo
! 
! call reset_timer(jjj)
! write(*,*) 'start cpu calculations'
! B=0.
! do i=1,nf
!  C=-A
!  do j=1,nnn
!   C(j,j)=C(j,j)+vec(i)
!  enddo
!  call invmat(n=nnn,mat=C)
!  B(:,:,i)=C
! enddo
! 
! write(*,*) 'CPU, real : ', maxval(abs(real(B(:,:,:))))
! write(*,*) 'CPU, imag : ', maxval(abs(aimag(B(:,:,:))))
! 
! write(*,*) 'done'
! call timer_fortran(jjj,'IT TOOK : ', unit_=6)
! 
! call cuda_array_of_inverse_collect(nnn,nf,A,c1,c2,vec,1)
! 
! write(*,*) 'GPU, real : ', maxval(abs(c1(:,:,:)))
! write(*,*) 'GPU, imag : ', maxval(abs(c2(:,:,:)))
! 
! write(*,*) 'real part comparison : ', maxval(abs(real(B)-c1))
! write(*,*) 'imag part comparison : ', maxval(abs(aimag(B)-c2))
! 
! stop
! 
! end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! subroutine test_invmat_array_matrix
! implicit none
! integer,parameter    :: nnn=20,nf=2000
! complex(8)           :: A(nnn,nnn),B(nnn,nnn),C(nnn,nnn),vec(nf)
! integer              :: i,j,n,jjj,jj
! 
! 
! do i=1,nnn
! do j=1,nnn
!  A(i,j)=drand1()
! enddo
! enddo
! do i=1,nf
!  vec(i)=drand1()
! enddo
! 
! call reset_timer(jjj)
! write(*,*) 'start cpu calculations'
! do jj=1,10
! B=0.
! do i=1,nf
!  C=-A
!  do j=1,nnn
!   C(j,j)=C(j,j)+vec(i)
!  enddo
!  call invmat(n=nnn,mat=C)
!  B=B+C
! enddo
! enddo
! write(*,*) 'done'
! call timer_fortran(jjj,'IT TOOK : ', unit_=6)
! 
!  write(*,*) 'start gpu calculations'
!  call reset_timer(jjj)
!  do i=1,10
!    if(mod(i,10)==0) write(*,*) i
!    if(i==1)then
!    call cuda_array_of_inverse(nnn,nf,A,C,vec,1)
!    elseif(i>1.and.i<100)then
!    call cuda_array_of_inverse(nnn,nf,A,C,vec,0)
!    elseif(i==100)then
!    call cuda_array_of_inverse(nnn,nf,A,C,vec,2)
!   endif
!  enddo
! 
!  call timer_fortran(jjj,'IT TOOK : ', unit_=6)
!  write(*,*) 'done'
! 
!  if(nnn<10)then
!  call write_array( B/dble(nf), " obtained by CPU ", unit = 6, short=.true. )
!  call write_array( C/dble(nf), " obtained by GPU C", unit = 6, short=.true. )
!  endif
! 
!  write(*,*) 'diff = ', maxval(abs(B-C))
! 
!  stop 'done'
! 
! end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! subroutine test_invmat_gpu
! implicit none
! 
! complex(8),allocatable  :: A(:,:),B(:,:)
! complex(8),allocatable  :: AAA(:,:),BBB(:,:),CCC(:,:)
! real(8),allocatable     :: Ar(:,:),Br(:,:)
! integer                 :: nnn,nnnn
! integer                 :: ii,i,j,k,l,m1,n1,q1
! 
! complex(8) :: Ac(8,8),Acinv(8,8),Bc_(8,8),Cc_(8,8),Dc_(8,8)
! real(8)    :: AA(16,24),BB(24,32),CC(16,32),DD(16,32)
! interface 
!  subroutine invert_ge(nn,a,inva)
!    integer :: nnn
!    real(8) :: a(nn,nn),inva(nn,nn)
!  end subroutine
!  subroutine invert_spd(nn,a,inva)
!    integer :: nnn
!    real(8) :: a(nn,nn),inva(nn,nn)
!  end subroutine
! end interface
! 
!  nnn=20
!  allocate(A(nnn,nnn),B(nnn,nnn),Ar(nnn,nnn),Br(nnn,nnn))
! 
!  do i=1,nnn
!   do j=1,nnn
!    Ar(i,j)=drand1() 
!   enddo
!   Ar(i,i)=100.+Ar(i,i)
!  enddo
!  Ar=Ar+TRANSPOSE(Ar)
! 
!  do i=1,nnn
!    write(*,'(300f9.3)') (Ar(i,j),j=1,nnn)
!  enddo
! 
!  Br=0.
!  call invert_ge(nnn,Ar,Br)
!  Ar=MATMUL(Ar,Br)
!  write(*,*) 'GE------'
!  do i=1,nnn
!    write(*,'(300f9.3)') (Ar(i,j),j=1,nnn)
!  enddo
! 
! 
!  do i=1,nnn
!   do j=1,nnn
!    Ar(i,j)=drand1() 
!   enddo
!   Ar(i,i)=100.+Ar(i,i)
!  enddo
!  Ar=Ar+TRANSPOSE(Ar)
!  Br=0.
!  call invert_spd(nnn,Ar,Br)
!  Ar=MATMUL(Ar,Br)
!  write(*,*) 'SPD------'
!  do i=1,nnn
!    write(*,'(300f9.3)') (Ar(i,j),j=1,nnn)
!  enddo
! 
! 
!  write(*,*) 'testing matinv_sym_complex'
!  nnn=8
!  do i=1,nnn
!   do j=1,nnn
!    Ac(i,j)=drand1() +imi*drand1()*j
!   enddo
!   Ac(i,i)=1000.*imi+1000.+Ac(i,i)
!  enddo
!  Ac=Ac+TRANSPOSE(Ac)
!  Acinv=Ac
!  call matinv_sym_complex(nnn,Acinv)
! 
!  Ac=MATMUL(Ac,Acinv)
!  do i=1,nnn
!  write(*,'(200f6.2)') (Ac(i,j),j=1,nnn)
!  enddo
! 
! 
!  write(*,*) 'testing matinv complex diago it -----------------'
!  do ii=1,9
!  write(*,*) '##############################################'
!  write(*,*) 'START NEW ITERATION : ',ii
!  if(ii==1) nnnn=5
!  if(ii==2) nnnn=10
!  if(ii==3) nnnn=300
!  if(ii==4) nnnn=600
!  if(ii==5) nnnn=670
!  if(ii==6) nnnn=1000
!  if(ii==7) nnnn=1200
!  if(ii==8) nnnn=1400
!  if(ii==9) nnnn=1670
!  if(allocated(AAA)) deallocate(AAA,BBB,CCC)
!  allocate(AAA(nnnn,nnnn),BBB(nnnn,nnnn),CCC(nnnn,nnnn))
!  do i=1,nnnn
!   do j=1,nnnn
!    AAA(i,j)=drand1() +imi*drand1()
!   enddo
!  enddo
!  BBB=AAA
!  write(*,*) 'start diag'
!  call diago_cuda_it_c(nnnn,BBB)
!  write(*,*) 'done diag'
!  call matmulcuda_c_cublas(AAA,BBB,CCC,nnnn,nnnn,nnnn)
!  write(*,*) 'max deviation from Id diago_it: ',ii, nnnn,maxval(abs(CCC-Id(nnnn)))
!  do i=1,nnnn
!   do j=1,nnnn
!    AAA(i,j)=drand1() +imi*drand1()
!   enddo
!  enddo
!  BBB=AAA
! 
! !write(*,*) 'start diag MAGMAg'
!  !call matinv_magma_complex(nnnn,BBB)
!   call matinv_magma_complex_cxx(nnnn,BBB)
! 
!  write(*,*) 'done diag'
!  call matmulcuda_c_cublas(AAA,BBB,CCC,nnnn,nnnn,nnnn)
!  write(*,*) 'max deviation from Id magma : ',ii, nnnn,maxval(abs(CCC-Id(nnnn)))
!  enddo
! 
!  write(*,*) 'testing matinv_magma_complex(n,mat)'
!  nnn=8
!  do i=1,nnn
!   do j=1,nnn
!    Ac(i,j)=drand1() +imi*drand1()*j
!   enddo
!  enddo
!  Acinv=Ac
!  call matinv_magma_complex(nnn,Acinv) 
!  Ac=MATMUL(Ac,Acinv)
!  do i=1,nnn
!   write(*,'(200f9.2)') (Ac(i,j),j=1,nnn)
!  enddo
! 
! 
!  write(*,*) '---->testing matinv_magma_double'
!  nnn=20
!  do i=1,nnn
!   do j=1,nnn
!    Ar(i,j)=drand1()
!   enddo
!  enddo
!  Br=Ar
!  call matinv_magma_double(nnn,Br)
!  Ar=MATMUL(Ar,Br)
!  write(*,*) 'GE------'
!  do i=1,nnn
!    write(*,'(300f9.3)') (Ar(i,j),j=1,nnn)
!  enddo
! 
!  write(*,*) 'testing matmulcuda_r_cublas(M,N,Q,dima,dimb,dimc)'
!  nnn=16
!  do i=1,16
!   do j=1,24
!    AA(i,j)=drand1() 
!   enddo
!  enddo
!  do i=1,24
!   do j=1,32
!    BB(i,j)=drand1()
!   enddo
!  enddo
!  m1=16
!  n1=24
!  q1=32
!  call matmulcuda_r_cublas(AA,BB,CC,m1,n1,q1)
!  DD=MATMUL(AA,BB)
!  write(*,*) 'max diff: '
!  write(*,*) maxval(abs(CC-DD)) 
!  write(*,*) 'out of cuda'
!   do i=1,16
!  write(*,'(200f9.2)') (CC(i,j),j=1,32)
!  enddo
!  write(*,*) 'out of cpu'
!  do i=1,16
!  write(*,'(200f9.2)') (DD(i,j),j=1,32)
!  enddo
! 
!  write(*,*) 'testing complex MATMUL on cublas'
! 
!  nnn=8
!  do i=1,nnn
!   do j=1,nnn
!    Ac(i,j)=drand1() +imi*drand1()*j
!   enddo
!  enddo
!  do i=1,nnn
!   do j=1,nnn
!    Bc_(i,j)=drand1() +imi*drand1()*j
!   enddo
!  enddo
!  call matmulcuda_c_cublas(Ac,Bc_,Cc_,nnn,nnn,nnn)
!  Dc_=MATMUL(Ac,Bc_)
!  write(*,*) 'max diff: '
!  write(*,*) maxval(abs(Cc_-Dc_))
!  write(*,*) 'out of cuda'
!   do i=1,8
!  write(*,'(200f9.2)') (Cc_(i,j),j=1,8)
!  enddo
!  write(*,*) 'out of cpu'
!  do i=1,8
!  write(*,'(200f9.2)') (Dc_(i,j),j=1,8)
!  enddo
! 
! 
! deallocate(A,B,Ar,Br)
! stop
! end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! subroutine test_invmat
! implicit none
! complex(16),allocatable :: Aqc(:,:),Bqc(:,:)
! real(16),allocatable    :: Aq(:,:),Bq(:,:)
! complex(8),allocatable  :: A(:,:),B(:,:)
! real(8),allocatable     :: Ar(:,:),Br(:,:)
! real(4),allocatable     :: As(:,:),Bs(:,:)
! complex(4),allocatable  :: Acs(:,:),Bcs(:,:)
! integer                 :: i,j,nnn,jjj,ii
! 
!  do ii=0,6
!  call define_flag(ii)
!  write(*,*) '################################'
!  write(*,*) 'STARTING PROCESS WITH FLAG TEST : ', ii
!  write(*,*) '################################'
!  do nnn=6,6
!   write(*,*) 'MATRIX SIZE : ', nnn
!   if(allocated(A)) deallocate(A,B,Br,Ar,As,Bs,Acs,Bcs,Aqc,Bqc,Aq,Bq)
!   allocate(A(nnn,nnn),B(nnn,nnn),Ar(nnn,nnn),Br(nnn,nnn),As(nnn,nnn),Bs(nnn,nnn),Acs(nnn,nnn),Bcs(nnn,nnn),&
!            & Aqc(nnn,nnn),Bqc(nnn,nnn),Aq(nnn,nnn),Bq(nnn,nnn))
!   call run_test(nnn)
!  enddo
!  enddo
!  stop 
! 
!  contains
! 
!  subroutine define_flag(flag)
!  implicit none
!  integer :: flag
!  fast_invmat=.false.
!  use_cuda_routines=.false.
!  use_openmp_invmat=.false.
!  diag_use_LU_instead_of_pivot=.false.
!  flag_use_invmat_jordan=.false.
!  flag_use_invmat_jordan_real=.false.
!  force_invmat_single_prec=.false.
! 
!  SELECT CASE(flag)
!   CASE(0)
! 
!   CASE(1)
!   fast_invmat=.true.
!   CASE(2)
!   use_cuda_routines=.true.
!   CASE(3)
!   use_openmp_invmat=.true.
!   CASE(4)
!   diag_use_LU_instead_of_pivot=.true.
!   CASE(5)
!   flag_use_invmat_jordan=.true.
!   flag_use_invmat_jordan_real=.true.
!   CASE(6)
!   force_invmat_single_prec=.true.
!  end SELECT
!  end subroutine
! 
!  subroutine run_test(n)
!  implicit none
!  integer :: n
!   write(*,*) '============================'
!   write(*,*) 'FULL MATRIX'
!   call define_mat(n,.false.)
!   call run_1(n,.false.)
!   if(mod(n,2)==0)then
!   write(*,*) 'BLOCK MATRIX'
!   call define_mat(n,.true.)
!   call run_1(n,.true.)
!   endif
!   write(*,*) '============================'
!  end subroutine
! 
! 
!  subroutine run_1(n,block)
!  implicit none
!  logical :: block
!  integer :: n
!  real(4) :: detr
!  complex(4)::det2r
!  real(8) :: det
!  real(16):: detq
!  complex(8) :: det2
!  complex(16) :: det2q
!  integer     :: pdet
! 
!      write(*,*) '--- QUAD ---'
!      call reset_timer(jjj)
!      call invmat(n,Aq,block_matrix=block,det2b=det2q,detb=detq,pdetb=pdet)
!      call timer_fortran(jjj,'IT TOOK : ', unit_=6)
!      write(*,*) maxval(abs(MATMUL(Bq,Aq)))
!      write(*,*) det2q,detq,pdet
!      write(*,*) '----------------'
! 
!      write(*,*) '--- QUAD COMPLEX ---'
!      call reset_timer(jjj)
!      call invmat(n,Aqc,block_matrix=block,det2b=det2q,detb=detq,pdetb=pdet)
!      call timer_fortran(jjj,'IT TOOK : ', unit_=6)
!      write(*,*) maxval(abs(MATMUL(Bqc,Aqc)))
!      write(*,*) det2q,detq,pdet
!      write(*,*) '----------------'
! 
!      write(*,*) '--- DOUBLE COMPLEX ---'
!      call reset_timer(jjj)
!      call invmat(n,A,block_matrix=block,det2b=det2,detb=det,pdetb=pdet)
!      call timer_fortran(jjj,'IT TOOK : ', unit_=6)
!      write(*,*) maxval(abs(MATMUL(B,A)))
!      write(*,*) det2,det,pdet
!      write(*,*) '----------------'
!      
!      write(*,*) '--- COMPLEX ---'
!      call reset_timer(jjj)
!      call invmat(n,Acs,block_matrix=block,det2b=det2r,detb=detr,pdetb=pdet)
!      call timer_fortran(jjj,'IT TOOK : ', unit_=6)
!      write(*,*) maxval(abs(MATMUL(Bcs,Acs)))
!      write(*,*) det2r,detr,pdet
!      write(*,*) '----------------'
!      
!      write(*,*) '--- DOUBLE ---'
!      call reset_timer(jjj)
!      call invmat(n,Ar,block_matrix=block,det2b=det2,detb=det,pdetb=pdet)
!      call timer_fortran(jjj,'IT TOOK : ', unit_=6)
!      write(*,*) maxval(abs(MATMUL(Br,Ar)))
!      write(*,*) det2,det,pdet
!      write(*,*) '----------------'
!      
!      write(*,*) '--- SINGLE ---'
!      call reset_timer(jjj)
!      call invmat(n,As,block_matrix=block,det2b=det2r,detb=detr,pdetb=pdet)
!      call timer_fortran(jjj,'IT TOOK : ', unit_=6)
!      write(*,*) maxval(abs(MATMUL(Bs,As)))
!      write(*,*) det2r,detr,pdet
!      write(*,*) '----------------'
!  end subroutine
! 
!  subroutine define_mat(n,block)
!  implicit none
!  logical :: block
!  integer :: n
!  do i=1,n
!   do j=1,n
!    A(i,j)=drand1()+imi*drand1()*j
!   enddo
!  enddo
!  if(block) then
!   A(1:n/2,n/2+1:n)=0.
!   A(n/2+1:n,1:n/2)=0.
!  endif
!  Ar =A
!  Br =A
!  B  =A
!  As =A
!  Acs=A
!  Bs =A
!  Bcs=A
!  do i=1,n
!   do j=1,n
!    Aq(i,j)=real(A(i,j))
!   enddo
!  enddo
!  do i=1,n
!   do j=1,n
!    Aqc(i,j)=real(A(i,j))
!   enddo
!  enddo
!  if(block) then
!   Aq(1:n/2,n/2+1:n)=0.
!   Aq(n/2+1:n,1:n/2)=0.
!   Aqc(1:n/2,n/2+1:n)=0.
!   Aqc(n/2+1:n,1:n/2)=0.
!  endif
!  Bq =Aq
!  Bqc=Aqc
!  end subroutine
! 
! 
! end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!   SUBROUTINE invert_zmblocktridiag(RES,z,A,B,N)
!     COMPLEX(8), INTENT(INOUT) :: RES(:,:)
!     COMPLEX(8), INTENT(IN)    :: z
!     REAL(8),    INTENT(IN)    :: A(:,:,:),B(:,:,:)
!     COMPLEX(8)                :: det_ratio_i(SIZE(A,2),SIZE(A,3)),det_ratio_ip1(SIZE(A,2),SIZE(A,3))
!     REAL(8)                   :: Id(SIZE(A,2),SIZE(A,3))
!     INTEGER                   :: iter,N
! 
!     IF(SIZE(A,1)/=SIZE(B,1).OR.ANY(SHAPE(A(1,:,:))/=SHAPE(B(1,:,:)))) &
!     STOP "ERROR IN invert_zmblocktridiag: INCONSISTENT INPUT DIMENSIONS!"
!     IF(SIZE(A,2)/=SIZE(A,3))  &
!     STOP "ERROR IN invert_zmblocktridiag: SQUARE BLOCKS EXPECTED!"
!     IF(SIZE(RES,1)/=SIZE(A,2).OR.SIZE(RES,2)/=SIZE(A,3)) &
!     STOP "ERROR IN invert_zmblocktridiag: INCONSISTENT RESOLVANT DIMENSIONS!"
! 
!     CALL new_Id(Id)
! 
!   ! COMPUTE <1/(z-T)> WHERE T IS A BLOCK-TRIDIAGONAL MATRIX AND <> IS THE FIRST BLOCK
! 
!     det_ratio_ip1 = z * Id - A(N,:,:)
!     CALL invert_lapack(det_ratio_ip1)
! 
!     DO iter = N-1,1,-1
!       det_ratio_i   =   z * Id - A(iter,:,:) - MATMUL(B(iter+1,:,:), &
!              & MATMUL(det_ratio_ip1,TRANSPOSE(CONJG( CMPLX(B(iter+1,:,:),0.d0,8) ))))
!       CALL invert_lapack(det_ratio_i)
!       det_ratio_ip1 = det_ratio_i
!     ENDDO
! 
!     RES = det_ratio_i
! 
!   END SUBROUTINE
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! subroutine invert_pivot_complexs_(matrix)
!  implicit none
!    complex(4) :: matrix(:,:)
!    integer    :: size_
!    integer    :: switch(size(matrix,1),2),k,jj,kp1,i,j,l,krow,irow
!    complex(4) :: pivot,temp,dkk
!  
!     size_=size(matrix,1)
! 
!         do  k = 1,size_
! 
!                 jj = k
!                 if (k .ne. size_) then
!                         kp1 = k + 1
!                         pivot = (matrix(k, k))
! !$OMP PARALLEL PRIVATE(i,pivot,temp,jj) 
! !$OMP DO
!                         do i = kp1,size_
!                                 temp = (matrix(i, k))
!                                 if ( abs(pivot) .lt.  abs(temp)) then
!                                         pivot = temp
!                                         jj = i
!                                 endif
!                         enddo
! !$OMP END DO
! !$OMP END PARALLEL
!                 endif
! 
!                 switch(k,1)=k; switch(k,2)=jj
! 
!                 if (jj .ne. k) then
! !$OMP PARALLEL PRIVATE(j,temp)
! !$OMP DO
!                     do  j = 1 ,size_
!                          temp = matrix(jj, j)
!                          matrix(jj, j) = matrix(k, j)
!                          matrix(k, j) = temp
!                     enddo
! !$OMP END DO
! !$OMP END PARALLEL
!                 endif
! 
!                 if(abs(matrix(k,k))>rerror)then
!                  dkk=matrix(k,k)
!                 else
!                  dkk=rerror
!                 endif
! 
! !$OMP PARALLEL PRIVATE(j)
! !$OMP DO
!                 do j = 1,size_
!                    if (j.ne.k) then 
!                       matrix(k,j) = matrix(k,j)/dkk
!                    endif
!                 enddo
! !$OMP END DO
! !$OMP END PARALLEL
! 
!               if(abs(matrix(k,k))>rerror)then
!                 matrix(k, k) = 1.d0 / matrix(k, k)
!               else
!                 matrix(k, k) = 1.d0 / rerror 
!               endif
! 
! !$OMP PARALLEL PRIVATE(i,j)
! !$OMP DO
!                 do  i = 1,size_
!                  if (i.ne.k) then
!                      do  j = 1,size_
!                        if(j.ne.k) matrix(i,j)=matrix(i,j)-matrix(k,j)*matrix(i,k)
!                      enddo
!                  endif
!                 enddo
! !$OMP END DO
! !$OMP END PARALLEL
! 
! 
! 
! !$OMP PARALLEL PRIVATE(i)
! !$OMP DO
!                 do i = 1, size_
!                         if (i.ne.k)  matrix(i,k) = -matrix(i,k) * matrix(k,k)
!                 enddo
! !$OMP END DO
! !$OMP END PARALLEL
! 
!         enddo
! 
!         do  l = 1,size_
!                 k = size_ - l + 1
!                 krow = switch(k, 1)
!                 irow = switch(k, 2)
!                 if (krow.ne.irow) then
! !$OMP PARALLEL PRIVATE(i,temp)
! !$OMP DO
!                         do  i = 1,size_
!                          temp = matrix(i, krow)
!                          matrix(i, krow) = matrix(i, irow)
!                          matrix(i, irow) = temp
!                         enddo
! !$OMP END DO
! !$OMP END PARALLEL
!                 endif
!         enddo
! 
! 
! end subroutine

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

! subroutine invert_pivot_single_(matrix)
!  implicit none
!    real(4) :: matrix(:,:)
!    integer :: size_
!    integer :: switch(size(matrix,1),2),k,jj,kp1,i,j,l,krow,irow
!    real(4) :: pivot,temp,dkk
!  
!     size_=size(matrix,1)
! 
!         do  k = 1,size_
! 
!                 jj = k
!                 if (k .ne. size_) then
!                         kp1 = k + 1
!                         pivot = (matrix(k, k))
! !$OMP PARALLEL PRIVATE(i,pivot,temp,jj) 
! !$OMP DO
!                         do i = kp1,size_
!                                 temp = (matrix(i, k))
!                                 if ( abs(pivot) .lt.  abs(temp)) then
!                                         pivot = temp
!                                         jj = i
!                                 endif
!                         enddo
! !$OMP END DO
! !$OMP END PARALLEL
!                 endif
! 
!                 switch(k,1)=k; switch(k,2)=jj
! 
!                 if (jj .ne. k) then
! !$OMP PARALLEL PRIVATE(j,temp)
! !$OMP DO
!                     do  j = 1 ,size_
!                          temp = matrix(jj, j)
!                          matrix(jj, j) = matrix(k, j)
!                          matrix(k, j) = temp
!                     enddo
! !$OMP END DO
! !$OMP END PARALLEL
!                 endif
! 
!                 if(abs(matrix(k,k))>rerror)then
!                  dkk=matrix(k,k)
!                 else
!                  dkk=rerror
!                 endif
! 
! !$OMP PARALLEL PRIVATE(j)
! !$OMP DO
!                 do j = 1,size_
!                    if (j.ne.k) then 
!                       matrix(k,j) = matrix(k,j)/dkk
!                    endif
!                 enddo
! !$OMP END DO
! !$OMP END PARALLEL
! 
!               if(abs(matrix(k,k))>rerror)then
!                 matrix(k, k) = 1.d0 / matrix(k, k)
!               else
!                 matrix(k, k) = 1.d0 / rerror 
!               endif
! 
! !$OMP PARALLEL PRIVATE(i,j)
! !$OMP DO
!                 do  i = 1,size_
!                  if (i.ne.k) then
!                      do  j = 1,size_
!                        if(j.ne.k) matrix(i,j)=matrix(i,j)-matrix(k,j)*matrix(i,k)
!                      enddo
!                  endif
!                 enddo
! !$OMP END DO
! !$OMP END PARALLEL
! 
! 
! 
! !$OMP PARALLEL PRIVATE(i)
! !$OMP DO
!                 do i = 1, size_
!                         if (i.ne.k)  matrix(i,k) = -matrix(i,k) * matrix(k,k)
!                 enddo
! !$OMP END DO
! !$OMP END PARALLEL
! 
!         enddo
! 
!         do  l = 1,size_
!                 k = size_ - l + 1
!                 krow = switch(k, 1)
!                 irow = switch(k, 2)
!                 if (krow.ne.irow) then
! !$OMP PARALLEL PRIVATE(i,temp)
! !$OMP DO
!                         do  i = 1,size_
!                          temp = matrix(i, krow)
!                          matrix(i, krow) = matrix(i, irow)
!                          matrix(i, irow) = temp
!                         enddo
! !$OMP END DO
! !$OMP END PARALLEL
!                 endif
!         enddo
! 
! 
! end subroutine

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

! subroutine invert_pivot_double_(matrix)
!  implicit none
!    real(8) :: matrix(:,:)
!    integer :: size_
!    integer :: switch(size(matrix,1),2),k,jj,kp1,i,j,l,krow,irow
!    real(8) :: pivot,temp,dkk
!  
!     size_=size(matrix,1)
! 
!         do  k = 1,size_
! 
!                 jj = k
!                 if (k .ne. size_) then
!                         kp1 = k + 1
!                         pivot = (matrix(k, k))
! !$OMP PARALLEL PRIVATE(i,pivot,temp,jj) 
! !$OMP DO
!                         do i = kp1,size_
!                                 temp = (matrix(i, k))
!                                 if ( abs(pivot) .lt.  abs(temp)) then
!                                         pivot = temp
!                                         jj = i
!                                 endif
!                         enddo
! !$OMP END DO
! !$OMP END PARALLEL
!                 endif
! 
!                 switch(k,1)=k; switch(k,2)=jj
! 
!                 if (jj .ne. k) then
! !$OMP PARALLEL PRIVATE(j,temp)
! !$OMP DO
!                     do  j = 1 ,size_
!                          temp = matrix(jj, j)
!                          matrix(jj, j) = matrix(k, j)
!                          matrix(k, j) = temp
!                     enddo
! !$OMP END DO
! !$OMP END PARALLEL
!                 endif
! 
!                 if(abs(matrix(k,k))>rerror)then
!                  dkk=matrix(k,k)
!                 else
!                  dkk=rerror
!                 endif
! 
! !$OMP PARALLEL PRIVATE(j)
! !$OMP DO
!                 do j = 1,size_
!                    if (j.ne.k) then 
!                       matrix(k,j) = matrix(k,j)/dkk
!                    endif
!                 enddo
! !$OMP END DO
! !$OMP END PARALLEL
! 
!               if(abs(matrix(k,k))>rerror)then
!                 matrix(k, k) = 1.d0 / matrix(k, k)
!               else
!                 matrix(k, k) = 1.d0 / rerror 
!               endif
! 
! !$OMP PARALLEL PRIVATE(i,j)
! !$OMP DO
!                 do  i = 1,size_
!                  if (i.ne.k) then
!                      do  j = 1,size_
!                        if(j.ne.k) matrix(i,j)=matrix(i,j)-matrix(k,j)*matrix(i,k)
!                      enddo
!                  endif
!                 enddo
! !$OMP END DO
! !$OMP END PARALLEL
! 
! 
! 
! !$OMP PARALLEL PRIVATE(i)
! !$OMP DO
!                 do i = 1, size_
!                         if (i.ne.k)  matrix(i,k) = -matrix(i,k) * matrix(k,k)
!                 enddo
! !$OMP END DO
! !$OMP END PARALLEL
! 
!         enddo
! 
!         do  l = 1,size_
!                 k = size_ - l + 1
!                 krow = switch(k, 1)
!                 irow = switch(k, 2)
!                 if (krow.ne.irow) then
! !$OMP PARALLEL PRIVATE(i,temp)
! !$OMP DO
!                         do  i = 1,size_
!                          temp = matrix(i, krow)
!                          matrix(i, krow) = matrix(i, irow)
!                          matrix(i, irow) = temp
!                         enddo
! !$OMP END DO
! !$OMP END PARALLEL
!                 endif
!         enddo
! 
! 
! end subroutine

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

! subroutine invert_pivot_complex_(matrix)
!  implicit none
!    complex(8) :: matrix(:,:)
!    integer    :: size_
!    integer    :: switch(size(matrix,1),2),k,jj,kp1,i,j,l,krow,irow
!    complex(8) :: pivot,temp,dkk
!  
!     size_=size(matrix,1)
! 
!         do  k = 1,size_
! 
!                 jj = k
!                 if (k .ne. size_) then
!                         kp1 = k + 1
!                         pivot = (matrix(k, k))
! !$OMP PARALLEL PRIVATE(i,pivot,temp,jj)
! !$OMP DO
!                         do i = kp1,size_
!                                 temp = (matrix(i, k))
!                                 if ( abs(pivot) .lt.  abs(temp)) then
!                                         pivot = temp
!                                         jj = i
!                                 endif
!                         enddo
! !$OMP END DO
! !$OMP END PARALLEL
!                 endif
! 
!                 switch(k, 1) = k
!                 switch(k, 2) = jj
! 
!                 if (jj .ne. k) then
! !$OMP PARALLEL PRIVATE(j,temp)
! !$OMP DO
!                     do  j = 1 ,size_
!                          temp = matrix(jj, j)
!                          matrix(jj, j) = matrix(k, j)
!                          matrix(k, j) = temp
!                     enddo
! !$OMP END DO
! !$OMP END PARALLEL
!                 endif
! 
!                 if(abs(matrix(k,k))>rerror)then
!                  dkk=matrix(k,k)
!                 else
!                  dkk=rerror
!                 endif
! 
! !$OMP PARALLEL PRIVATE(j)
! !$OMP DO
!                 do j = 1,size_
!                    if (j.ne.k) then 
!                       matrix(k,j) = matrix(k,j) / dkk
!                    endif
!                 enddo
! !$OMP END DO
! !$OMP END PARALLEL
! 
!               if(abs(matrix(k,k))>rerror)then
!                 matrix(k, k) = 1.d0 / matrix(k, k)
!               else
!                 matrix(k, k) = 1.d0 / rerror 
!               endif
! 
! !$OMP PARALLEL PRIVATE(i,j)
! !$OMP DO
!                 do  i = 1,size_
!                  if (i.ne.k) then
!                      do  j = 1,size_
!                        if(j.ne.k) matrix(i,j)=matrix(i,j)-matrix(k,j)*matrix(i,k)
!                      enddo
!                  endif
!                 enddo
! !$OMP END DO
! !$OMP END PARALLEL
! 
! 
! 
! !$OMP PARALLEL PRIVATE(i)
! !$OMP DO
!                 do i = 1, size_
!                         if (i.ne.k)  matrix(i, k) = -matrix(i, k) * matrix(k, k)
!                 enddo
! !$OMP END DO
! !$OMP END PARALLEL
! 
!         enddo
! 
!         do  l = 1,size_
!                 k = size_ - l + 1
!                 krow = switch(k, 1)
!                 irow = switch(k, 2)
!                 if (krow.ne.irow) then
! !$OMP PARALLEL PRIVATE(i,temp)
! !$OMP DO
!                         do  i = 1,size_
!                          temp = matrix(i, krow)
!                          matrix(i, krow) = matrix(i, irow)
!                          matrix(i, irow) = temp
!                         enddo
! !$OMP END DO
! !$OMP END PARALLEL
!                 endif
!         enddo
! 
! 
! end subroutine

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

function MATMULr(AA,BB)
implicit none
real(8) :: AA(:,:),BB(:,:),MATMULr(size(AA,1),size(BB,2))
integer :: i,j,k,ii
 i=size(AA,1); j=size(AA,2); k=size(BB,2)
 ii=max(max(i,j),k)
 ! ebl: Removing GPU functionality
 ! if(ii>n_cuda_rout.and.use_cuda_routines)then
 !  CALL matmulcuda_r(AA,BB,MATMULr,i,j,k)
 !  return
 ! endif
 MATMULr=MATMUL(AA,BB)
end function

 !----------------------------------------!

function MATMULc(AA,BB)
implicit none
complex(8) :: AA(:,:),BB(:,:),MATMULc(size(AA,1),size(BB,2))
integer    :: i,j,k,ii
 i=size(AA,1); j=size(AA,2); k=size(BB,2)
 ii=max(max(i,j),k)
 ! ebl: Removing GPU functionality
 ! if(ii>n_cuda_rout.and.use_cuda_routines)then
 !  CALL matmulcuda_c(AA,BB,MATMULc,i,j,k)
 !  return
 ! endif
 MATMULc=MATMUL(AA,BB)
end function

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
! 
! SUBROUTINE tred2(a,d,e,novectors)
!   IMPLICIT NONE
! 
!   REAL(8), DIMENSION(:,:), INTENT(inout) :: a
!   REAL(8), DIMENSION(:), INTENT(out)     :: d,e
!   LOGICAL, OPTIONAL, INTENT(in)          :: novectors
!   INTEGER                                :: i,j,l,n
!   REAL(8)                                :: f,g,h,hh,tscale
!   REAL(8), DIMENSION(SIZE(a,1))          :: gg
!   LOGICAL, SAVE                          :: yesvec=.TRUE.
!   
!   n=SIZE(a,1)
!   
!   IF (PRESENT(novectors)) yesvec = .NOT. novectors
!   DO i=n,2,-1
!      l=i-1
!      h=0.0
!      IF(l>1) THEN
!         tscale=SUM(ABS(a(i,1:l)))
!         IF(tscale==0.0) THEN
!            e(i) = a(i,l)
!         ELSE
!            a(i,1:l) = a(i,1:l)/tscale
!            h=SUM(a(i,1:l)**2)
!            f=a(i,l)
!            g=-SIGN(SQRT(h),f)
!            e(i) = tscale*g
!            h=h-f*g
!            a(i,l)=f-g
!            IF(yesvec) a(1:l,i) = a(i,1:l)/h
!            DO j=1,l
!               e(j)=(DOT_PRODUCT(a(j,1:j),a(i,1:j)) &
!                    + DOT_PRODUCT(a(j+1:l,j),a(i,j+1:l)))/h
!            END DO
!            f=DOT_PRODUCT(e(1:l),a(i,1:l))
!            hh=f/(h+h)
!            e(1:l)=e(1:l)-hh*a(i,1:l)
!            DO j=1,l
!               a(j,1:j)=a(j,1:j) - a(i,j)*e(1:j)-e(j)*a(i,1:j)
!            END DO
!         END IF
!      ELSE
!         e(i)=a(i,l)
!      END IF
!      d(i)=h
!   END DO
!   IF(yesvec) d(1)=0.0
!   e(1)=0.0
!   DO i=1,n
!      IF(yesvec) THEN
!         l=i-1
!         IF(d(i) /=0.0) THEN
!            gg(1:l)=MATMUL(a(i,1:l),a(1:l,1:l))
!            a(1:l,1:l) = a(1:l,1:l) - outerprod(a(1:l,i),gg(1:l))
!         END IF
!         d(i)=a(i,i)
!         a(i,i)=1.0
!         a(i,1:l)=0.0
!         a(1:l,i)=0.0
!      ELSE
!         d(i)=a(i,i)
!      END IF
!   END DO
! END SUBROUTINE 
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! SUBROUTINE tqli__(d,e,z)
!   IMPLICIT NONE
! 
!   REAL(8), DIMENSION(:), INTENT(inout)            :: d,e
!   REAL(8),DIMENSION(:,:), OPTIONAL, INTENT(inout) :: z
!   INTEGER                                         :: i,iter,l,mm,n,ndum
!   REAL(8)                                         :: b,c,dd,f,g,p,r,s
!   REAL(8), DIMENSION(SIZE(e))                     :: ff
! 
!   n=SIZE(d)
!   IF(PRESENT(z)) ndum=n
!   e(:)=EOSHIFT(e(:),1) !convenient to renumber the elements of e
!   DO l=1,n
!      iter=0
!      iterate: DO
!         DO mm=l,n-1
!            dd=ABS(d(mm))+ABS(d(mm+1))
!            IF(ABS(e(mm)) + dd == dd) EXIT
!         END DO
!         IF(mm==l) EXIT iterate
!         IF(iter==100) WRITE(*,*) 'too many iterations in tqli'
!         iter=iter+1
!         g=(d(l+1)-d(l))/(2.d0*e(l))!Form shift
!         r=pythag(g,one)
!         g=d(mm)-d(l)+e(l)/(g+SIGN(r,g)) !this is d_m - k_s.
!         s=one
!         c=one
!         p=zero
!         DO i =mm-1,l,-1
!            f=s*e(i)
!            b=c*e(i)
!            r=pythag(f,g)
!            e(i+1)=r
!            IF(r==zero) THEN
!               d(i+1)=d(i+1)-p
!               e(mm)=zero
!               CYCLE iterate
!            END IF
!            s=f/r
!            c=g/r
!            g=d(i+1)-p
!            r=(d(i)-g)*s+2.d0*c*b
!            p=s*r
!            d(i+1)=g+p
!            g=c*r-b
!            IF(PRESENT(z)) THEN
!               ff(1:n)=z(1:n,i+1)
!               z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
!               z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
!            END IF
!         END DO
!         d(l)=d(l)-p
!         e(l)=g
!         e(mm)=zero
!      END DO iterate
!   END DO
!   
! 
!   
! CONTAINS
! 
!  
!   FUNCTION pythag(a,b)
!     IMPLICIT NONE
!     REAL(8), INTENT(in) :: a,b
!     REAL(8)             :: pythag
!     REAL(8)             :: absa,absb
!     absa = ABS(a)
!     absb = ABS(b)
!     IF(absa > absb) THEN
!        pythag=absa*SQRT(1.d0 + (absb/absa)**2)
!     ELSE
!        IF(absb==0.0) THEN
!           pythag=0.0
!        ELSE
!           pythag=absb*SQRT(1.d0 + (absa/absb)**2)
!        END IF
!     END IF
!   END FUNCTION 
! 
! 
! END SUBROUTINE 
! 
! 
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

!           !-------------------------------------!
! 
!     subroutine get_det_from_zgeco___(deti,det,det2,pdet)
!     implicit none
!     complex(4) :: deti(2),det2
!     real(4)    :: det
!     integer    :: pdet
!      ! convert it to the form 2**pdet * abs(=det) * phase(=det2)
!         pdet=NINT(real(deti(2))/LOG10(2.d0))
!         det2=(deti(1))*(2.d0**(real(deti(2))/LOG10(2.d0)-dble(pdet)))
!         det=abs(det2)
!         if(det>1.d-23)then
!          det2=det2/det
!         else
!          det2=1.
!         endif
!         call rescale(det,pdet)
!     end subroutine
! 
!           !-------------------------------------!
! 
!     subroutine get_det_from_zgeco__(deti,det,det2,pdet)
!     implicit none
!     complex(16) :: deti(2),det2
!     real(16)    :: det
!     integer     :: pdet
!      ! convert it to the form 2**pdet * abs(=det) * phase(=det2)
!         pdet=NINT(real(deti(2))/LOG10(2.d0))
!         det2=(deti(1))*(2.d0**(real(deti(2))/LOG10(2.d0)-dble(pdet)))
!         det=abs(det2)
!         if(det>1.d-23)then
!          det2=det2/det
!         else
!          det2=1.
!         endif
!         call rescale(det,pdet)
!     end subroutine
! 
!           !-------------------------------------!
! 
!     subroutine get_det_from_zgeco_(deti,det,det2,pdet)
!     implicit none
!     complex(8) :: deti(2),det2
!     real(8)    :: det
!     integer    :: pdet
!      ! convert it to the form 2**pdet * abs(=det) * phase(=det2)
!         pdet=NINT(real(deti(2))/LOG10(2.d0))
!         det2=(deti(1))*(2.d0**(real(deti(2))/LOG10(2.d0)-dble(pdet)))
!         det=abs(det2)
!         call erase_divergence(det)
!         call erase_divergence(det2)
!         if(det>1.d-12)then
!          det2=det2/det
!         else
!          det2=1.
!         endif
!         call rescale(det,pdet)
!     end subroutine

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
! 
!       !------------------------------------!
! 
!  real(8) function trace_rp(mat)
!  implicit none
!  complex(8)::mat(:,:)
!  integer::k
!  trace_rp=0.d0
!  do k=1,size(mat,1)
!  trace_rp=trace_rp+real(mat(k,k),8)
!  enddo
!  end function
! 
!  real(8) function trace_ip(mat)
!  implicit none
!  complex(8)::mat(:,:)
!  integer::k
!  trace_ip=0.d0
!  do k=1,size(mat,1)
!  trace_ip=trace_ip+aimag(mat(k,k))
!  enddo
!  end function
! 
      !------------------------------------!

 function diagr(mat)
 implicit none
 real(8)::mat(:,:)
 real(8)::diagr(size(mat(1,:)))
 integer::k
   diagr=(/(mat(k,k),k=1,size(mat(1,:)))/)
 end function

      !------------------------------------!

 function diagr_(mat)
 implicit none
 real(8) :: mat(:,:,:)
 real(8) :: diagr_(size(mat(:,1,1)),size(mat(1,1,:)))
 integer :: k,i
  do i=1,size(mat(1,1,:))
   diagr_(:,i)=(/( mat(k,k,i), k=1,size(mat(:,1,1))  )/)
  enddo
 end function

      !------------------------------------!

 function diagc_(mat,i1,i2)
 implicit none
 complex(8) :: mat(:,:,:)
 complex(8) :: diagc_(size(mat(:,1,1)),size(mat(1,1,:)))
 integer    :: k,i,i1,i2
  do i=1,size(mat(1,1,:))
   diagc_(:,i)=(/(mat(k,k,i),k=1,size(mat(:,1,1)))/)
  enddo
 end function

     !------------------------------------!

 function offdiagr(mat)
 implicit none
 real(8) :: mat(:,:)
 real(8) :: offdiagr(size(mat(:,1)),size(mat(1,:)))
 integer :: i
  offdiagr=mat
  do i=1,size(mat(1,:))
   offdiagr(i,i)=0.
  enddo
 end function

      !------------------------------------!

 function offdiagc(mat)
 implicit none
 complex(8) :: mat(:,:)
 complex(8) :: offdiagc(size(mat(:,1)),size(mat(1,:)))
 integer    :: i
  offdiagc=mat
  do i=1,size(mat(1,:))
   offdiagc(i,i)=0.
  enddo
 end function

      !------------------------------------!

 function diagr__(mat)
 implicit none
 real(8)::mat(:,:,:,:)
 real(8)::diagr__(size(mat(:,1,1,1)),size(mat(1,:,1,1)),size(mat(1,1,1,:)))
 integer::k,i,j
  do j=1,size(mat(:,1,1,1))
   do i=1,size(mat(1,1,1,:))
     diagr__(j,:,i)=(/(mat(j,k,k,i),k=1,size(mat(1,:,1,1)))/)
   enddo
  enddo
 end function

      !------------------------------------!
 
 function diagc__(mat)
 implicit none
 complex(8)::mat(:,:,:,:)
 complex(8)::diagc__(size(mat(:,1,1,1)),size(mat(1,:,1,1)),size(mat(1,1,1,:)))
 integer::k,i,j
  do j=1,size(mat(:,1,1,1))
   do i=1,size(mat(1,1,1,:))
     diagc__(j,:,i)=(/(mat(j,k,k,i),k=1,size(mat(1,:,1,1)))/)
   enddo
  enddo
 end function

      !------------------------------------!

 function diagrr(mat)
 implicit none
 real(4)::mat(:,:)
 real(4)::diagrr(size(mat(1,:)))
 integer::k
   diagrr=(/(mat(k,k),k=1,size(mat(1,:)))/)
 end function

      !------------------------------------!

 function diagc(mat)
 implicit none
 complex(8)::mat(:,:)
 complex(8)::diagc(size(mat(1,:)))
 integer::k
   diagc=(/(mat(k,k),k=1,size(mat(1,:)))/)
 end function

      !------------------------------------!

 function diagi(mat)
 implicit none
 integer::mat(:,:)
 integer::diagi(size(mat(1,:)))
 integer::k
   diagi=(/(mat(k,k),k=1,size(mat(1,:)))/)
 end function
! 
!       !------------------------------------!
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!  pure integer function linear_size_hermitian_matrix_for_vec(siz)
!  implicit none
!  integer,intent(in) :: siz
!  real(8)            :: tt
! 
!    !le vecteur d entree contient la diagonale + les elements superieures
!    tt=(-1.d0+sqrt(1.d0+8.d0*dble(siz)))/2.d0
!    linear_size_hermitian_matrix_for_vec = NINT(tt)
!    tt=tt-NINT(tt)
! 
!    if(abs(tt)>1.d-3.or.tt<0.d0)then
!       !le vecteur d entree ne contient que la diagonale
!       linear_size_hermitian_matrix_for_vec = siz
!    endif
! 
!  return
!  end function
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       !------------------------------------!
!       !------------------------------------!
!       !------------------------------------!
! 
!     subroutine test_invmat_tri
!     implicit none
!     real(8)    :: mat(10,10),mat_backup(10,10)
!     complex(8) :: mat_(10,10),mat_backup_(10,10)
!     integer    :: i,j
!  
!      mat=0.d0
!      do i=1,10
!                 mat(i,i  )=drand1()+0.1
!        if(i<10) mat(i,i+1)=drand1()
!      enddo
!      do i=1,9
!      do j=i+1,10
!       mat(j,i)=mat(i,j)
!      enddo
!      enddo
!      call write_array( mat, 'MATRIX TO INVERT', SHORT=.true., UNIT=6 )
!      mat_backup=mat
!        
!      call invmat_tridiag(10,mat)
!      write(*,*) 'inverse?? : '
!      call write_array( MATMUL(mat,mat_backup), 'MAt * INV', SHORT=.true., UNIT=6) 
! 
!      mat_        = mat_backup
!      mat_backup_ = mat_
!      do i=1,10
!       mat_backup_(i,i)=mat_backup_(i,i)+imi*2.d0
!      enddo
! 
!      call invmat_tridiag_complex(10,mat_,imi*2.d0)
! 
!      call write_array( REAL(MATMUL(mat_backup_,mat_)), 'MAt * INV', SHORT=.true., UNIT=6)
!      call write_array( AIMAG(MATMUL(mat_backup_,mat_)), 'MAt * INV', SHORT=.true., UNIT=6)
! 
!      stop 'done'
!     end subroutine
! 
!       !------------------------------------!
!       !------------------------------------!
!       !------------------------------------!
! 
!     subroutine invmat_tridiag_complex(siz,mat_,iw_)
!     implicit none
!     integer    :: siz
!     real(8)    :: mat(siz,siz),eigenvalues(siz)
!     complex(8) :: iw_,mat_(siz,siz)
!       mat=real(mat_)
!       call eigenvector_tridiag(siz,eigenvalues,mat)
!       where(abs(eigenvalues)<1.d-13) eigenvalues=1.d-13
!       mat_= MATMUL ( mat, MATMUL(bande_mat(siz,1.d0/(eigenvalues+iw_)),transpose(mat)))
!     end subroutine
! 
!       !------------------------------------!
!  
!     subroutine invmat_tridiag(siz,mat)
!     implicit none
!     integer :: siz
!     real(8) :: mat(siz,siz),eigenvalues(siz)
!       call eigenvector_tridiag(siz,eigenvalues,mat)   
!       where(abs(eigenvalues)<1.d-13) eigenvalues=1.d-13
!       mat= MATMUL ( mat, MATMUL(bande_mat(siz,1.d0/eigenvalues),transpose(mat)) )
!     end subroutine
! 
!       !------------------------------------!
! 
!     subroutine eigenvector_tridiag(siz,eigenvalues,mat)
!     implicit none
!     integer :: siz
!     real(8) :: mat(siz,siz),eigenvalues(siz)
!     INTEGER :: INFO,i
!     REAL(8) :: sdiag(siz-1),work(2*siz-2)
!      eigenvalues =      diag(mat)
!      sdiag       =  (/( mat(i,i+1),i=1,siz-1 )/)
!      call DSTEV('V',siz,eigenvalues,sdiag,mat,siz,WORK,INFO)
!     end subroutine
!   
!       !------------------------------------!
! 
!     subroutine eigenvector_matrixa_(sizin,mat,covd)
!     implicit none
!     integer :: sizin
!     real(8) :: cove(sizin),covd(sizin)
!     real(8) :: mat(sizin,sizin)
!       call tred(mat,sizin,sizin,covd,cove)
!       call tqli(covd,cove,sizin,sizin,mat)
!     return
!     end subroutine
! 
!       !------------------------------------!
! 
!     subroutine eigenvector_matrixb_(sizin,mat,covd)
!     implicit none
!     integer  :: sizin
!     real(8)  :: cove(sizin),covd(sizin)
!     real(8)  :: mat(sizin,sizin)
!     real(16) :: matq(sizin,sizin),coveq(sizin),covdq(sizin)
!       matq=mat
!       call tredq(matq,sizin,sizin,covdq,coveq)
!       call tqliq(covdq,coveq,sizin,sizin,matq)
!       covd=covdq
!     return
!     end subroutine
! 
!       !------------------------------------!
! 
!     subroutine eigenvector_matrixc_(sizin,matq,qcovdq)
!     implicit none
!     integer   :: sizin
!     real(16)  :: cove(sizin),qcovdq(sizin)
!     real(16)  :: matq(sizin,sizin),coveq(sizin),covdq(sizin)
!       call tredq(matq,sizin,sizin,qcovdq,coveq)
!       call tqliq(qcovdq,coveq,sizin,sizin,matq)
!     return
!     end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       SUBROUTINE TQLI_(D,E,N,NP,Z)
!       implicit real(8)(a-h,o-z)
!       DIMENSION D(NP),E(NP),Z(NP,NP)
!       DIMENSION WKSP(Np),IWKSP(Np)
!       IF (N.GT.1) THEN
!       DO 11 I=2,N
!       E(I-1)=E(I)
! 11    CONTINUE
!       E(N)=0.d0
!       DO 15 L=1,N
!       ITER=0
! 1     DO 12 M=L,N-1
!       DD=ABS(D(M))+ABS(D(M+1))
!       IF (ABS(E(M))+DD.EQ.DD) GO TO 2
! 12    CONTINUE
!       M=N
! 2     IF(M.NE.L)THEN
!       ITER=ITER+1
!       G=(D(L+1)-D(L))/(2.d0*E(L))
!       R=SQRT(G**2+1.d0)
!       G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
!       S=1.d0
!       C=1.d0
!       P=0.d0
!       DO 14 I=M-1,L,-1
!       F=S*E(I)
!       B=C*E(I)
!       IF(ABS(F).GE.ABS(G))THEN
!       C=G/F
!       R=SQRT(C**2+1.d0)
!       E(I+1)=F*R
!       S=1.d0/R
!       C=C*S
!       ELSE
!       S=F/G
!       R=SQRT(S**2+1.d0)
!       E(I+1)=G*R
!       C=1.d0/R
!       S=S*C
!       ENDIF
!       G=D(I+1)-P
!       R=(D(I)-G)*S+2.d0*C*B
!       P=S*R
!       D(I+1)=G+P
!       G=C*R-B
!       DO 13 K=1,N
!       F=Z(K,I+1)
!       Z(K,I+1)=S*Z(K,I)+C*F
!       Z(K,I)=C*Z(K,I)-S*F
! 13    CONTINUE
! 14    CONTINUE
!       D(L)=D(L)-P
!       E(L)=G
!       E(M)=0.d0
!       GO TO 1
!       ENDIF
! 15    CONTINUE
!       ENDIF
!       call SORT3(N,np,D,Z,WKSP,IWKSP)
!       RETURN
!       END subroutine
! 
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       SUBROUTINE TQLIq(D,E,N,NP,Z)
!       implicit real(16)(a-h,o-z)
!       DIMENSION D(NP),E(NP),Z(NP,NP)
!       DIMENSION WKSP(Np),IWKSP(Np)
!       IF (N.GT.1) THEN
!       DO 11 I=2,N
!       E(I-1)=E(I)
! 11    CONTINUE
!       E(N)=0.d0
!       DO 15 L=1,N
!       ITER=0
! 1     DO 12 M=L,N-1
!       DD=ABS(D(M))+ABS(D(M+1))
!       IF (ABS(E(M))+DD.EQ.DD) GO TO 2
! 12    CONTINUE
!       M=N
! 2     IF(M.NE.L)THEN
!       ITER=ITER+1
!       G=(D(L+1)-D(L))/(2.d0*E(L))
!       R=SQRT(G**2+1.d0)
!       G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
!       S=1.d0
!       C=1.d0
!       P=0.d0
!       DO 14 I=M-1,L,-1
!       F=S*E(I)
!       B=C*E(I)
!       IF(ABS(F).GE.ABS(G))THEN
!       C=G/F
!       R=SQRT(C**2+1.d0)
!       E(I+1)=F*R
!       S=1.d0/R
!       C=C*S
!       ELSE
!       S=F/G
!       R=SQRT(S**2+1.d0)
!       E(I+1)=G*R
!       C=1.d0/R
!       S=S*C
!       ENDIF
!       G=D(I+1)-P
!       R=(D(I)-G)*S+2.d0*C*B
!       P=S*R
!       D(I+1)=G+P
!       G=C*R-B
!       DO 13 K=1,N
!       F=Z(K,I+1)
!       Z(K,I+1)=S*Z(K,I)+C*F
!       Z(K,I)=C*Z(K,I)-S*F
! 13    CONTINUE
! 14    CONTINUE
!       D(L)=D(L)-P
!       E(L)=G
!       E(M)=0.d0
!       GO TO 1
!       ENDIF
! 15    CONTINUE
!       ENDIF
!       call SORT3q(N,np,D,Z,WKSP,IWKSP)
!       RETURN
!       END subroutine
! 
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       ! Nom : TRED
!       SUBROUTINE TRED(A,N,NP,D,E)
!       implicit real(8)(a-h,o-z)
!       DIMENSION A(NP,NP),D(NP),E(NP)
!       IF(N.GT.1)THEN
!       DO 18 I=N,2,-1
!       L=I-1
!       H=0.d0
!       SCALE=0.d0
!       IF(L.GT.1)THEN
!       DO 11 K=1,L
!       SCALE=SCALE+ABS(A(I,K))
! 11    CONTINUE
!       IF(SCALE<epsilonr)THEN
!       E(I)=A(I,L)
!       ELSE
!       DO 12 K=1,L
!       A(I,K)=A(I,K)/SCALE
!       H=H+A(I,K)**2
! 12    CONTINUE
!       F=A(I,L)
!       G=-SIGN(SQRT(H),F)
!       E(I)=SCALE*G
!       H=H-F*G
!       A(I,L)=F-G
!       F=0.d0
!       DO 15 J=1,L
!       A(J,I)=A(I,J)/H
!       G=0.d0
!       DO 13 K=1,J
!       G=G+A(J,K)*A(I,K)
! 13    CONTINUE
!       IF(L.GT.J)THEN
!       DO 14 K=J+1,L
!       G=G+A(K,J)*A(I,K)
! 14    CONTINUE
!       ENDIF
!       E(J)=G/H
!       F=F+E(J)*A(I,J)
! 15    CONTINUE
!       HH=F/(H+H)
!       DO 17 J=1,L
!       F=A(I,J)
!       G=E(J)-HH*F
!       E(J)=G
!       DO 16 K=1,J
!       A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
! 16    CONTINUE
! 17    CONTINUE
!       ENDIF
!       ELSE
!       E(I)=A(I,L)
!       ENDIF
!       D(I)=H
! 18    CONTINUE
!       ENDIF
!       D(1)=0.d0
!       E(1)=0.d0
!       DO 23 I=1,N
!       L=I-1
!       IF(D(I).NE.0.d0)THEN
!       DO 21 J=1,L
!       G=0.d0
!       DO 19 K=1,L
!       G=G+A(I,K)*A(K,J)
! 19    CONTINUE
!       DO 20 K=1,L
!       A(K,J)=A(K,J)-G*A(K,I)
! 20    CONTINUE
! 21    CONTINUE
!       ENDIF
!       D(I)=A(I,I)
!       A(I,I)=1.d0
!       IF(L.GE.1)THEN
!       DO 22 J=1,L
!       A(I,J)=0.d0
!       A(J,I)=0.d0
! 22    CONTINUE
!       ENDIF
! 23    CONTINUE
!       RETURN
!       END subroutine
! 
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       ! Nom : TRED
!       SUBROUTINE TREDq(A,N,NP,D,E)
!       implicit real(16)(a-h,o-z)
!       DIMENSION A(NP,NP),D(NP),E(NP)
!       IF(N.GT.1)THEN
!       DO 18 I=N,2,-1
!       L=I-1
!       H=0.d0
!       SCALE=0.d0
!       IF(L.GT.1)THEN
!       DO 11 K=1,L
!       SCALE=SCALE+ABS(A(I,K))
! 11    CONTINUE
!       IF(SCALE<epsilonq)THEN
!       E(I)=A(I,L)
!       ELSE
!       DO 12 K=1,L
!       A(I,K)=A(I,K)/SCALE
!       H=H+A(I,K)**2
! 12    CONTINUE
!       F=A(I,L)
!       G=-SIGN(SQRT(H),F)
!       E(I)=SCALE*G
!       H=H-F*G
!       A(I,L)=F-G
!       F=0.d0
!       DO 15 J=1,L
!       A(J,I)=A(I,J)/H
!       G=0.d0
!       DO 13 K=1,J
!       G=G+A(J,K)*A(I,K)
! 13    CONTINUE
!       IF(L.GT.J)THEN
!       DO 14 K=J+1,L
!       G=G+A(K,J)*A(I,K)
! 14    CONTINUE
!       ENDIF
!       E(J)=G/H
!       F=F+E(J)*A(I,J)
! 15    CONTINUE
!       HH=F/(H+H)
!       DO 17 J=1,L
!       F=A(I,J)
!       G=E(J)-HH*F
!       E(J)=G
!       DO 16 K=1,J
!       A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
! 16    CONTINUE
! 17    CONTINUE
!       ENDIF
!       ELSE
!       E(I)=A(I,L)
!       ENDIF
!       D(I)=H
! 18    CONTINUE
!       ENDIF
!       D(1)=0.d0
!       E(1)=0.d0
!       DO 23 I=1,N
!       L=I-1
!       IF(D(I).NE.0.d0)THEN
!       DO 21 J=1,L
!       G=0.d0
!       DO 19 K=1,L
!       G=G+A(I,K)*A(K,J)
! 19    CONTINUE
!       DO 20 K=1,L
!       A(K,J)=A(K,J)-G*A(K,I)
! 20    CONTINUE
! 21    CONTINUE
!       ENDIF
!       D(I)=A(I,I)
!       A(I,I)=1.d0
!       IF(L.GE.1)THEN
!       DO 22 J=1,L
!       A(I,J)=0.d0
!       A(J,I)=0.d0
! 22    CONTINUE
!       ENDIF
! 23    CONTINUE
!       RETURN
!       END subroutine
! 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

! Fonction        : Inversion d'une matrice carree
! Notes           : A -> A^(-1)=AI
!                   Provient de Numerical Recipes, p.40.
! Dependances     : LUDCMP,LUBKSB (externes)
! 
!      SUBROUTINE D_INVMG(NP,NM,A)
!      IMPLICIT NONE
!      real(8)                                     :: D          ! Var. locale
!      INTEGER                                     :: I,J        ! var. locale
!      INTEGER,INTENT(IN)                          :: NP,NM
!      INTEGER,DIMENSION (NP)                      :: INDX       ! Var. locale
!      real(8),DIMENSION (NP,NP)                   :: Y          ! Var. locale
!      real(8),DIMENSION (NP,NP),INTENT(INOUT)     :: A
! 
!      DO I=1,NM
!         DO J=1,NM
!            Y(I,J)=0.D0
!         END DO
!         Y(I,I)=1.D0
!      END DO
!      CALL LUDCMP(A,NM,NP,INDX,D)
!      DO J=1,NM
!         CALL LUBKSBd(A,NM,NP,INDX,Y(1,J))
!      END DO
!      A=Y
! 
!      END SUBROUTINE
! 
! 
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

! Nom            : LUDCMP 
! Fonction       : Decomposition gaussienne A = LU
! Notes          : Provient entierement de Numerical Recipes
! 
!    SUBROUTINE LUDCMPd(A,N,NP,INDX,D)
!    IMPLICIT real(8) (A-H,O-Z)
!    PARAMETER (TINY=1.0D-20)
!    DIMENSION A(NP,NP),INDX(N),VV(N)
! 
!           D=1.D0
!           DO 12 I=1,N
!            AAMAX=0.D0
!            DO 11 J=1,N
!               IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
! 11         CONTINUE
!            IF(AAMAX.EQ.0.) AAMX=1.d-22
!            VV(I)=1.D0/AAMAX
! 12           CONTINUE
!             DO 19 J=1,N
!            IF(J.GT.1) THEN
!               DO 14 I=1,J-1
!                  SUM=A(I,J)
!                  IF(I.GT.1) THEN
!                     DO 13 K=1,I-1
!                        SUM=SUM-A(I,K)*A(K,J)
! 13                  CONTINUE
!                     A(I,J)=SUM
!                  ENDIF
! 14            CONTINUE
!            ENDIF
!            AAMAX=0.D0
!            DO 16 I=J,N
!               SUM=A(I,J)
!               IF(J.GT.1) THEN
!                  DO 15 K=1,J-1
!                     SUM=SUM-A(I,K)*A(K,J)
! 15               CONTINUE
!                  A(I,J)=SUM
!               ENDIF
!               DUM=VV(I)*ABS(SUM)
!               IF(DUM.GE.AAMAX) THEN
!                  IMAX=I
!                  AAMAX=DUM
!               ENDIF
! 16         CONTINUE
!            IF(J.NE.IMAX) THEN
!               DO 17 K=1,N
!                  DUM=A(IMAX,K)
!                  A(IMAX,K)=A(J,K)
!                  A(J,K)=DUM
! 17            CONTINUE
!               D=-D
!               VV(IMAX)=VV(J)
!            ENDIF
!            INDX(J)=IMAX
!            IF(J.NE.N) THEN
!              IF(A(J,J).EQ.0.) A(J,J)=TINY
!              DUM=1./A(J,J)
!              DO 18 I=J+1,N
!                 A(I,J)=A(I,J)*DUM
! 18           CONTINUE
!            ENDIF
! 19          CONTINUE
!             IF(A(N,N).EQ.0.) A(N,N)=TINY
!             RETURN
!             END subroutine
! 
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

! Nom            : LUDCMP 
! Fonction       : Decomposition gaussienne A = LU
! Notes          : Provient entierement de Numerical Recipes
! 
!    SUBROUTINE LUDCMPq(A,N,NP,INDX,D)
!    IMPLICIT real(16) (A-H,O-Z)
!    PARAMETER (TINY=1.0D-34)
!    DIMENSION A(NP,NP),INDX(N),VV(N)
!  
!           D=1.D0
!           DO 12 I=1,N
!            AAMAX=0.D0
!            DO 11 J=1,N
!               IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
! 11         CONTINUE
!            IF(AAMAX.EQ.0.) AAMX=1.d-22
!            VV(I)=1.D0/AAMAX
! 12           CONTINUE
!             DO 19 J=1,N
!            IF(J.GT.1) THEN
!               DO 14 I=1,J-1
!                  SUM=A(I,J)
!                  IF(I.GT.1) THEN
!                     DO 13 K=1,I-1
!                        SUM=SUM-A(I,K)*A(K,J)
! 13                  CONTINUE
!                     A(I,J)=SUM
!                  ENDIF
! 14            CONTINUE
!            ENDIF
!            AAMAX=0.D0
!            DO 16 I=J,N
!               SUM=A(I,J)
!               IF(J.GT.1) THEN
!                  DO 15 K=1,J-1
!                     SUM=SUM-A(I,K)*A(K,J)
! 15               CONTINUE
!                  A(I,J)=SUM
!               ENDIF
!               DUM=VV(I)*ABS(SUM)
!               IF(DUM.GE.AAMAX) THEN
!                  IMAX=I
!                  AAMAX=DUM
!               ENDIF
! 16         CONTINUE
!            IF(J.NE.IMAX) THEN
!               DO 17 K=1,N
!                  DUM=A(IMAX,K)
!                  A(IMAX,K)=A(J,K)
!                  A(J,K)=DUM
! 17            CONTINUE
!               D=-D
!               VV(IMAX)=VV(J)
!            ENDIF
!            INDX(J)=IMAX
!            IF(J.NE.N) THEN
!              IF(A(J,J).EQ.0.) A(J,J)=TINY
!              DUM=1./A(J,J)
!              DO 18 I=J+1,N
!                 A(I,J)=A(I,J)*DUM
! 18           CONTINUE
!            ENDIF
! 19          CONTINUE
!             IF(A(N,N).EQ.0.) A(N,N)=TINY
!             RETURN
!             END subroutine
! 
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

! Nom            : LUBKSB  
! Fonction       : Backsubstitution (pour l'inversion)
! Notes          : Provient de Numerical Recipes
! 
!       SUBROUTINE LUBKSBd(A,N,NP,INDX,B)
!       INTEGER N,NP,INDX(N)
!       real(8) A(NP,NP),B(N)
!       INTEGER I,II,J,LL
!       real(8) SUM
!       II=0
!       DO I=1,N
!          LL=INDX(I)
!          SUM=B(LL)
!          B(LL)=B(I)
!          IF(II.NE.0) THEN
!             DO J=II,I-1
!                SUM=SUM-A(I,J)*B(J)
!             END DO
!          ELSE IF(SUM.NE.0) THEN
!             II=I
!          END IF
!          B(I)=SUM
!       END DO
!       DO I=N,1,-1
!          SUM=B(I)
!          DO J=I+1,N
!             SUM=SUM-A(I,J)*B(J)
!          END DO
!          B(I)=SUM/A(I,I)
!       END DO
!       RETURN
!       END subroutine
! 
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

! Nom            : LUBKSB  
! Fonction       : Backsubstitution (pour l'inversion)
! Notes          : Provient de Numerical Recipes
! 
!       SUBROUTINE LUBKSBq(A,N,NP,INDX,B)
!       INTEGER N,NP,INDX(N)
!       real(16) A(NP,NP),B(N)
!       INTEGER I,II,J,LL
!       real(16) SUM
!       II=0
!       DO I=1,N
!          LL=INDX(I)
!          SUM=B(LL)
!          B(LL)=B(I)
!          IF(II.NE.0) THEN
!             DO J=II,I-1
!                SUM=SUM-A(I,J)*B(J)
!             END DO
!          ELSE IF(SUM.NE.0) THEN
!             II=I
!          END IF
!          B(I)=SUM
!       END DO
!       DO I=N,1,-1
!          SUM=B(I)
!          DO J=I+1,N
!             SUM=SUM-A(I,J)*B(J)
!          END DO
!          B(I)=SUM/A(I,I)
!       END DO
!       RETURN
!       END subroutine

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
! 
! !     Fonction:  Single value decomposition
! 
!       SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
!       integer,intent(in)    :: M,N
!       real(8),intent(inout) :: A(MP,NP),W(NP),V(NP,NP)
!       real(8)               :: RV1(N)
!       real(8)               :: S,G,H,F,SCALE,ANORM
!       INTEGER               :: I,J,L,K
! 
!       G=0.0
!       SCALE=0.0; ANORM=0.0
! 
!       DO 25 I=1,N
!         L=I+1
!         RV1(I)=SCALE*G
!         G=0.0
!         S=0.0
!         SCALE=0.0
!         IF (I.LE.M) THEN
!           DO 11 K=I,M
!             SCALE=SCALE+ABS(A(K,I))
! 11        CONTINUE
!           IF (SCALE.NE.0.0) THEN
!             DO 12 K=I,M
!               A(K,I)=A(K,I)/SCALE
!               S=S+A(K,I)*A(K,I)
! 12          CONTINUE
!             F=A(I,I)
!             G=-SIGN(SQRT(S),F)
!             H=F*G-S
!             A(I,I)=F-G
!             IF (I.NE.N) THEN
!               DO 15 J=L,N
!                 S=0.0
!                 DO 13 K=I,M
!                   S=S+A(K,I)*A(K,J)
! 13              CONTINUE
!                 F=S/H
!                 DO 14 K=I,M
!                   A(K,J)=A(K,J)+F*A(K,I)
! 14              CONTINUE
! 15            CONTINUE
!             ENDIF
!             DO 16 K= I,M
!               A(K,I)=SCALE*A(K,I)
! 16          CONTINUE
!           ENDIF
!         ENDIF
!         W(I)=SCALE *G
!         G=0.0
!         S=0.0
!         SCALE=0.0
! 
!         IF ((I<=M).AND.(I/=N)) THEN
!           DO 17 K=L,N
!             SCALE=SCALE+ABS(A(I,K))
! 17        CONTINUE
!           IF (SCALE.NE.0.0) THEN
!             DO 18 K=L,N
!               A(I,K)=A(I,K)/SCALE
!               S=S+A(I,K)*A(I,K)
! 18          CONTINUE
!             F=A(I,L)
!             G=-SIGN(SQRT(S),F)
!             H=F*G-S
!             A(I,L)=F-G
!             DO 19 K=L,N
!               RV1(K)=A(I,K)/H
! 19          CONTINUE
!             IF (I.NE.M) THEN
!               DO 23 J=L,M
!                 S=0.0
!                 DO 21 K=L,N
!                   S=S+A(J,K)*A(I,K)
! 21              CONTINUE
!                 DO 22 K=L,N
!                   A(J,K)=A(J,K)+S*RV1(K)
! 22              CONTINUE
! 23            CONTINUE
!             ENDIF
!             DO 24 K=L,N
!               A(I,K)=SCALE*A(I,K)
! 24          CONTINUE
!           ENDIF
!         ENDIF
!         ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
! 25    CONTINUE
!       DO 32 I=N,1,-1
!         IF (I<N) THEN
!           IF (G/=0.0) THEN
!             DO 26 J=L,N
!               if(I>M) stop 'error pivot, bad shape 1'
!               V(J,I)=(A(I,J)/A(I,L))/G
! 26          CONTINUE
!             DO 29 J=L,N
!               S=0.0
!               DO 27 K=L,N
!                 if(I>M) stop 'error pivot, bad shape 2'
!                 S=S+A(I,K)*V(K,J)
! 27            CONTINUE
!               DO 28 K=L,N
!                 V(K,J)=V(K,J)+S*V(K,I)
! 28            CONTINUE
! 29          CONTINUE
!           ENDIF
!           DO 31 J=L,N
!             V(I,J)=0.0
!             V(J,I)=0.0
! 31        CONTINUE
!         ENDIF
!         V(I,I)=1.0
!         G=RV1(I)
!         L=I
! 32    CONTINUE
!       DO 39 I=N,1,-1
!         L=I+1
!         G=W(I)
!         IF (I.LT.N) THEN
!           DO 33 J=L,N
!             A(I,J)=0.0
! 33        CONTINUE
!         ENDIF
!         IF (G.NE.0.0) THEN
!           G=1.0/G
!           IF (I/=N) THEN
!             DO 36 J=L,N
!               S=0.0
!               DO 34 K=L,M
!                 S=S+A(K,I)*A(K,J)
! 34            CONTINUE
!               F=(S/A(I,I))*G
!               DO 35 K=I,M
!                 A(K,J)=A(K,J)+F*A(K,I)
! 35            CONTINUE
! 36          CONTINUE
!           ENDIF
!           DO 37 J=I,M
!             A(J,I)=A(J,I)*G
! 37        CONTINUE
!         ELSE
!           DO 38 J= I,M
!             A(J,I)=0.0
! 38        CONTINUE
!         ENDIF
!         A(I,I)=A(I,I)+1.0
! 39    CONTINUE
!       DO 49 K=N,1,-1
!         DO 48 ITS=1,30
!           DO 41 L=K,1,-1
!             NM=L-1
!             IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
!             IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
! 41        CONTINUE
! 1         C=0.0
!           S=1.0
!           DO 43 I=L,K
!             F=S*RV1(I)
!             IF ((ABS(F)+ANORM).NE.ANORM) THEN
!               G=W(I)
!               H=SQRT(F*F+G*G)
!               W(I)=H
!               H=1.0/H
!               C= (G*H)
!               S=-(F*H)
!               DO 42 J=1,M
!                 Y=A(J,NM)
!                 Z=A(J,I)
!                 A(J,NM)=(Y*C)+(Z*S)
!                 A(J,I)=-(Y*S)+(Z*C)
! 42            CONTINUE
!             ENDIF
! 43        CONTINUE
! 2         Z=W(K)
!           IF (L.EQ.K) THEN
!             IF (Z.LT.0.0) THEN
!               W(K)=-Z
!               DO 44 J=1,N
!                 V(J,K)=-V(J,K)
! 44            CONTINUE
!             ENDIF
!             GO TO 3
!           ENDIF
!           X=W(L)
!           NM=K-1
!           Y=W(NM)
!           G=RV1(NM)
!           H=RV1(K)
!           F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
!           G=SQRT(F*F+1.0)
!           F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
!           C=1.0
!           S=1.0
!           DO 47 J=L,NM
!             I=J+1
!             G=RV1(I)
!             Y=W(I)
!             H=S*G
!             G=C*G
!             Z=SQRT(F*F+H*H)
!             RV1(J)=Z
!             C=F/Z
!             S=H/Z
!             F= (X*C)+(G*S)
!             G=-(X*S)+(G*C)
!             H=Y*S
!             Y=Y*C
!             DO 45 NM=1,N
!               X=V(NM,J)
!               Z=V(NM,I)
!               V(NM,J)= (X*C)+(Z*S)
!               V(NM,I)=-(X*S)+(Z*C)
! 45          CONTINUE
!             Z=SQRT(F*F+H*H)
!             W(J)=Z
!             IF (Z.NE.0.0) THEN
!               Z=1.0/Z
!               C=F*Z
!               S=H*Z
!             ENDIF
!             F= (C*G)+(S*Y)
!             X=-(S*G)+(C*Y)
!             DO 46 NM=1,M
!               Y=A(NM,J)
!               Z=A(NM,I)
!               A(NM,J)= (Y*C)+(Z*S)
!               A(NM,I)=-(Y*S)+(Z*C)
! 46          CONTINUE
! 47        CONTINUE
!           RV1(L)=0.0
!           RV1(K)=F
!           W(K)=X
! 48      CONTINUE
! 3       CONTINUE
! 49    CONTINUE
! 
!       RETURN
!       END subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!   !*************************************************!
!   !*     resolution avec la decomposition en SVD    !
!   !*************************************************!
! 
!       subroutine svbksb(u,w,v,m,n,mp,np,b,x)
!       implicit REAL(8) (a-h,o-z)
!       dimension u(mp,np),w(np),v(np,np),b(mp),x(np),tmp(n)
! 
!       do 12 j=1,n
!             s=0.d0
!             if(w(j).ne.0.d0)then
!                     do 11 i=1,m
!                           s=s+u(i,j)*b(i)
! 11                  continue
!                     s=s/w(j)
!             endif
!             tmp(j)=s
! 12    continue
! 
!       do 14 j=1,n
!             s=0.d0
!             do 13 jj=1,n
!                   s=s+v(j,jj)*tmp(jj)
! 13          continue
!             x(j)=s
! 14    continue
!       return
!       end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! subroutine fix_degeneracies(n,val,ordertab,err)
! implicit none
! integer :: n,i,j,k,l
! real(8)  :: err,val(n)
! integer :: ordertab(n)
! 
! !---> ordertab(i)=k : state number i degenerate with the k-1 following state
! 
!  ordertab=0; i=0
! 
!  do 
!   i=i+1
!   if(i>n) exit
!   ordertab(i)=1
!   k=1
!   do j=i+1,n
!    if(abs(val(j)-val(i))<err)then
!     k=k+1
!     ordertab(i)=k
!    else
!     exit
!    endif
!   enddo
!   i=i+k-1
!  enddo
!  
!  k=sum(ordertab(1:n))
!  if(k/=n)then
!   write(*,*) 'should be ... states : ', n
!   write(*,*) 'there are ... states : ', k
!   write(*,*) 'ordertab tab : ', ordertab
!   write(*,*) 'values    : ', val
!   stop 'error fix_degeneracies bad check....'
!  endif
! 
! return
! end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! subroutine perturbation_eigenvectors(n,vecin,valin,Hin,Vin,eigenvec,deg, &
!        & use_nelson,V_not_deg,use_lishu,first_order,permut,permut_typ,half_only)
! implicit none
! integer                     :: n,i,j,k,l,m,nnn,iloop
! complex(8)                  :: temp1,temp2,temp3
! real(8)                     :: valin(n)
! real(8)                     :: eigenval(n)
! complex(8)                  :: vecin(n,n),vec_back(n,n),perturb(n,n),unitary(n,n)
! complex(8)                  :: Vin(n,n),Hin(n,n)
! complex(8)                  :: eigenvec(n,n),ak
! integer                     :: ordertab(n)
! integer                     :: nstate,n_in
! logical                     :: cancel_second_term,degenerate
! logical                     :: use_permut
! 
!  !=================================================================!
!     complex(8),optional         :: permut_typ(n,2)
!     logical,optional            :: use_nelson,V_not_deg,deg,use_lishu
!     integer,optional            :: permut(n,2)
!     logical,optional            :: first_order
!     logical,optional            :: half_only
!  !=================================================================!
! 
! 
!  !=======================================================================!
!  ! Routine donne la correction aux valeurs propres et vecteurs propres   !
!  ! vecin : composantes vecteurs propres en ligne                         !
!  !=======================================================================!
! 
!  use_permut=.false.
!  if(present(permut))then
!    use_permut=.true.
!    if(maxval(abs(permut))==0) use_permut=.false.
!  endif
!  if(present(half_only)) then
!   n_in=n/2+1
!  else
!   n_in=n
!  endif
! 
! !---------------------------------------------------!
! call fix_degeneracies(n,valin,ordertab,0.02d0)
!                        degenerate=maxval(ordertab)>1
! if(present(deg))       deg=degenerate
!                        cancel_second_term=.false.
! if(present(V_not_deg)) cancel_second_term=.true.
! 
! eigenvec=0.d0;eigenval=0.;perturb=0.
! !---------------------------------------------------!
! 
! if(degenerate.and.cancel_second_term) call remove_degeneracies
! 
! !---------------------------------------------------!
! if(.not.degenerate.and.present(use_nelson))then
!   call first_order_NELSON
!   return
! endif
! 
! if(.not.use_permut)then
!  call build_V_matrix_full
! else
!  call build_V_matrix_sparse
! endif
! !---------------------------------------------------!
! 
! 
! !************* FIRST ORDER ************!
! 
! do iloop=1,n_in
!  if(ordertab(iloop)==1)then
!      call first_order_not_deg
!   elseif(ordertab(iloop)>1)then
!     if(present(use_lishu))then
!      call first_order_LiShu(iloop,iloop+ordertab(iloop)-1)
!     else
!      call first_order_deg(iloop,iloop+ordertab(iloop)-1)
!     endif
!  endif
! enddo
! 
! if(cancel_second_term) vecin=vec_back
! 
!  
! if(degenerate.or.present(first_order)) then
!  return
! endif
! 
! 
! !*********** SECOND ORDER *************!
! if(messages) write(*,*) 'go for second order perturbation theory...'
! do iloop=1,n
!  if(ordertab(iloop)==0) stop 'error second order pert. theory'
!  if(ordertab(iloop)==1) call second_order_not_deg
! enddo
! !**************************************!
! 
! return
! 
! contains
! 
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
! 
!   subroutine remove_degeneracies
!    call unitary_matrix(n,perturb,unitary)
!    vec_back=vecin
!    do i=1,n
!     if(ordertab(i)>1)then
!      do j=i,i+ordertab(i)-1
!       vecin(j,:)=MATMUL(unitary,vecin(j,:))
!      enddo
!    endif
!    enddo
!   end subroutine
! 
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
! 
!    subroutine build_V_matrix_sparse
!    implicit none
!    integer :: i,j,k,l
!    complex(8) :: tempv
!     do i=1,n
!      do j=1,n
!        perturb(i,j)=0.d0
!         do k=1,n
!          tempv=Vin(k,k)*vecin(j,k)
!          if(abs(permut(k,1))>0) tempv=tempv+permut_typ(k,1)*vecin(j,abs(permut(k,1)))
!          if(abs(permut(k,2))>0) tempv=tempv+permut_typ(k,2)*vecin(j,abs(permut(k,2)))
!          perturb(i,j)=perturb(i,j)+conjg(vecin(i,k))*tempv   
!        enddo
!      enddo
!     enddo
!    return
!    end subroutine
! 
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
! 
!    subroutine build_V_matrix_full
!    implicit none
!    integer :: i,j
!     do i=1,n
!       do j=1,n
!        perturb(i,j)=SCALPROD(vecin(i,:),MATMUL(Vin,vecin(j,:)))
!       enddo
!     enddo
!    end subroutine
! 
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
! 
!  subroutine first_order_LiShu(j1,j2)
!  implicit none
!  integer    :: j1,j2
!  complex(8) :: dtemp,phi1(n,j2-j1+1),phi2(n,n-(j2-j1+1))
!  complex(8) :: diag(n-(j2-j1+1),n-(j2-j1+1))
!  real(8)     :: norm
!  integer    :: i,j,k,r,kk
! 
!  r=j2-j1+1; k=0; diag=0.d0
!  do i=1,n
!   if(i<j1.or.i>j2)then
!    k=k+1
!    norm=valin(j1)-valin(i)
!    if(abs(norm)<1.d-4) then
!      write(*,*) 'error LiShu'
!      norm=1.d-3
!      if(strongstop) stop 'error : strongstop activated therefore stop'
!    endif
!    diag(k,k)=1.d0/norm
!   endif
!  enddo
! 
!   do i=1,j1-1
!      phi2(:,i)=vecin(i,:)
!   enddo
! 
!   if(.not.use_permut)then
!     do i=j1,j2
!        phi1(:,i-j1+1)=vecin(i,:)
!     enddo
!   else
!     do i=j1,j2
!      do k=1,n
!        phi1(k,i-j1+1)=Vin(k,k)*vecin(i,k)
!        if(abs(permut(k,1))>0) &
!            & phi1(k,i-j1+1)=phi1(k,i-j1+1) + vecin(i,abs(permut(k,1)))*permut_typ(k,1)
!        if(abs(permut(k,2))>0) &
!            & phi1(k,i-j1+1)=phi1(k,i-j1+1) + vecin(i,abs(permut(k,2)))*permut_typ(k,2)
!      enddo
!     enddo
!   endif
! 
!   do i=j2+1,n
!      phi2(:,i-j2+j1-1)=vecin(i,:)
!   enddo
! 
!  if(.not.use_permut) then
!      phi1=MATMUL(phi2,MATMUL_x(diag,MATMUL &
!                & (MATMUL(TRANSPOSE(CONJG(phi2)),Vin),phi1),IdL=.true.))
!  else
!      phi1=MATMUL(phi2,MATMUL_x(diag, & 
!                & MATMUL(TRANSPOSE(CONJG(phi2)),phi1),IdL=.true.))
!  endif
! 
!  k=0
!  do i=j1,j2
!   k=k+1
!   eigenvec(i,:)=phi1(:,k)
!  enddo
! 
!  end subroutine
! 
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
! 
! subroutine first_order_NELSON
! complex(8) :: h(n,n),f(n),d(n,n),v(n)
! complex(8) :: ci,dtemp
! integer    :: i,j,ii,k,jjj
! 
!  do jjj=1,n
!    f       =  MATMUL(Vin,vecin(jjj,:))
!    dtemp   =  scalprod(vecin(jjj,:),f)
!    f       =  dtemp*vecin(jjj,:)
!    f       =  f - MATMUL(Vin,vecin(jjj,:))
!    f(jjj)  =  0.d0
!    h       =  Hin
!    do i=1,n
!      h(i,i)=h(i,i)-valin(jjj)
!    enddo
!    h(jjj,:)=0.; h(:,jjj)=0.; h(jjj,jjj)=1.d0
!    call invmat(n,h)
!    eigenvec(jjj,:)=MATMUL(h,f)
!    ci=-(scalprod(vecin(jjj,:),eigenvec(jjj,:)))
!    eigenvec(jjj,:)=eigenvec(jjj,:)+ci*vecin(jjj,:)
!  enddo
! 
! end subroutine
! 
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
! 
! subroutine first_order_not_deg
!    do j=1,n
!     if(j/=iloop)then
!       ak = perturb(j,iloop)/(valin(iloop)-valin(j))
!     else
!       ak = 0.d0  !=1.d0; si pas correction mais vecteur
!     endif
!     if(abs(ak)>1.d-3) eigenvec(iloop,:) = eigenvec(iloop,:) + ak*vecin(j,:)
!   enddo
! end subroutine
! 
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
! 
! subroutine second_order_not_deg
! integer :: nn
! 
!  do j=1,n
!   ak=0.
!   if(j/=iloop)then
!      do nn=1,n
!       if(nn/=iloop)then
!         if(abs(valin(iloop)-valin(j))<1.d-3.and.strongstop) stop 'error 1 second_order_not_deg'
!         if(abs(valin(iloop)-valin(nn))<1.d-3.and.strongstop) stop 'error 2 second_order_not_deg'
!         ak=ak + perturb(j,nn)*perturb(iloop,nn)/(valin(iloop)-valin(nn))/(valin(iloop)-valin(j))
!         ak=ak - perturb(iloop,iloop)*perturb(j,iloop)/((valin(iloop)-valin(j))**2)
!       endif
!      enddo
!   else
!      do nn=1,n
!       if(nn/=iloop)then
!         if(abs(valin(iloop)-valin(nn))<1.d-3) stop 'error 3 second_order_not_deg'
!         ak = ak - 0.5d0 * (abs(perturb(iloop,nn))**2)  /  ((valin(iloop)-valin(nn))**2)
!       endif
!      enddo
!   endif
!   eigenvec(iloop,:)=eigenvec(iloop,:)+ak*vecin(j,:)
!  enddo
! 
! end subroutine
! 
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
! 
! subroutine first_order_deg(j1,j2)
! integer    :: jj,j1,j2,ii,j,k,l,m
! integer    :: q,qm,qj,nn,nnp,Beta
! complex(8) :: Cq,norm,normb
! 
!  do q=j1,j2
! 
!         !---------------------------------------!
!    if(.not.cancel_second_term)then
!     do qm=j1,j2
! 
!      norm=perturb(q,q)-perturb(qm,qm)
! 
!      if(abs(norm)>1.d-3)then
!  
!      if(abs(norm)<1.d-3) then
!        write(*,*) 'trouble, perturbation to hamiltonian has identical diagonal value'
!        write(*,*) 'q,qm : ', q,qm
!        write(*,*) 'Vq,Vqm : ', perturb(q,q),perturb(qm,qm)
!        write(*,*) '=======> perturbation does not break entirely the degenerescence......'
!        if(strongstop) stop
!      endif
! 
!      if(qm/=q)then
!      if(abs(norm)>1.d-3)then
! 
!      Cq=0.
!      do nnp=1,n
!       normb=valin(q)-valin(nnp)
!       if(abs(normb)>1.d-3)then
!         Cq = Cq + perturb(qm,nnp)*perturb(nnp,q)/normb
!       endif
!      enddo
! 
!      eigenvec(q,:) = eigenvec(q,:) + Cq /norm * vecin(qm,:)
! 
!      endif
!      endif
!      endif
!     enddo
!   endif
!         !---------------------------------------!
! 
!     do nnp=1,n
!      normb=valin(q)-valin(nnp)
!      if(abs(normb)>1.d-3)then
!         eigenvec(q,:) = eigenvec(q,:) + vecin(nnp,:)*perturb(nnp,q)/normb
!      endif
!     enddo
! 
!         !---------------------------------------!
! 
!  enddo
! 
! end subroutine
! 
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
!   !-----------------!
! 
! end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! function MATMUL_line(a,b,s2)
! implicit none
! complex(8) :: a(:),b(:,:),MATMUL_line(size(a(:)))
! integer    :: i,j,k,s2,siz
! siz=size(a(:))
! MATMUL_line=0.
!  do j=1,siz
!   do k=1,size(b(:,1))
!    MATMUL_line(j)=MATMUL_line(j)+a(k)*b(k,j)
!   enddo
!  enddo
! return
! end function
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! function MATMUL_keep_diag_(aa,bb)
! implicit none
! complex(8) :: aa(:,:),bb(:,:),MATMUL_keep_diag_(size(aa(:,1)))
! integer    :: i,j,k,l,siz1,siz2
! siz1=size(aa(:,1))
! siz2=size(bb(1,:))
! do i=1,siz1
!  MATMUL_keep_diag_(i)=0.
!  do k=1,siz2
!    MATMUL_keep_diag_(i)=MATMUL_keep_diag_(i) + aa(i,k)*bb(k,i) 
!  enddo
! enddo
! end function
! 
!        !-----------------------------------!
! 
! function MATMUL_keep_diag__(aa,bb)
! implicit none
! real(8)    :: aa(:,:),bb(:,:),MATMUL_keep_diag__(size(aa(:,1)))
! integer    :: i,j,k,l,siz1,siz2
! siz1=size(aa(:,1))
! siz2=size(bb(1,:))
! do i=1,siz1
!  MATMUL_keep_diag__(i)=0.
!  do k=1,siz2
!    MATMUL_keep_diag__(i)=MATMUL_keep_diag__(i) + aa(i,k)*bb(k,i)
!  enddo
! enddo
! end function
! 
!        !-----------------------------------!
! 
! function MATMUL_sum_diag_(aa,bb)
! implicit none
! complex(8) :: aa(:,:),bb(:,:),MATMUL_sum_diag_
! integer    :: i,j,k,l,siz1,siz2
! siz1=size(aa(:,1))
! siz2=size(bb(1,:))
! MATMUL_sum_diag_=0.
! do i=1,siz1
!  do k=1,siz2
!    MATMUL_sum_diag_=MATMUL_sum_diag_+aa(i,k)*bb(k,i)
!  enddo
! enddo
! end function
! 
!        !-----------------------------------!
! 
! function MATMUL_sum_diag__(aa,bb)
! implicit none
! real(8)    :: aa(:,:),bb(:,:),MATMUL_sum_diag__
! integer    :: i,j,k,l,siz1,siz2
! siz1=size(aa(:,1))
! siz2=size(bb(1,:))
! MATMUL_sum_diag__=0.
! do i=1,siz1
!  do k=1,siz2
!    MATMUL_sum_diag__=MATMUL_sum_diag__ + aa(i,k)*bb(k,i)
!  enddo
! enddo
! end function
! 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

function MATMUL_x_c(aa,bb,ind,n,IdL,IdR,byblock,a_by_block)
implicit none
integer                :: ii,i,j,k,s1,s2,smiddle
complex(8),intent(in)  :: aa(:,:),bb(:,:)
complex(8)             :: MATMUL_x_c(size(aa(:,1)),size(bb(1,:)))
integer,optional       :: n,ind(:)
logical,optional       :: IdL,IdR,byblock,a_by_block

if(present(byblock))then
 if(byblock)then
   ii=size(aa(:,1))
   i =size(aa(1,:))
   if(i/=ii) stop 'MATMUL_x by block should be only for square matrices'
   if(mod(ii,2)>0) stop 'MATMUL_x by block should only be for even linear sizes'
   ii=ii/2
   MATMUL_x_c(1    :   ii, ii+1 : 2*ii)= 0.0
   MATMUL_x_c(ii+1 : 2*ii,    1 :   ii)= 0.0
   MATMUL_x_c(   1 :   ii,    1 :   ii)= MATMUL_(aa(    1:   ii,    1 :   ii),bb(    1 :  ii,    1 :  ii))
   MATMUL_x_c(ii+1 : 2*ii, ii+1 : 2*ii)= MATMUL_(aa(ii +1: 2*ii, ii+1 : 2*ii),bb( ii+1 :2*ii, ii+1 :2*ii))
   return
 endif
endif

!ind: matrice aa n a que les lignes dans ind qui different de l identite

s1=size(aa(:,1)); s2=size(bb(1,:)); smiddle=size(aa(1,:))
if(smiddle/=size(bb(:,1))) stop 'error MATMUL_x_c'

if(present(IdL))then
 if(testing)then
  if(abs(aa(1,2))>1.d-3) stop 'MATMUL_x_c bad diagonal matrix IdL'
 endif
 do i=1,size(bb(:,1))
  MATMUL_x_c(i,:)=bb(i,:)*aa(i,i)
 enddo
 return
endif

if(present(IdR))then
 if(testing)then
  if(abs(aa(1,2))>1.d-3) stop 'MATMUL_x_c bad diagonal matrix IdR'
 endif
 do i=1,size(bb(1,:))
  MATMUL_x_c(:,i)=bb(:,i)*aa(i,i)
 enddo
 return
endif

 if(present(ind))then
  MATMUL_x_c=bb
  if(size(ind)/=n) stop 'error MATMUL_x_c bad ind shape'
  do ii=1,n
   i=ind(ii)
   do j=1,s2
    MATMUL_x_c(i,j)=0.
    do k=1,smiddle
     MATMUL_x_c(i,j)=MATMUL_x_c(i,j)+aa(i,k)*bb(k,j)
    enddo
   enddo
  enddo
 return
 endif

 if(present(a_by_block))then
 if(a_by_block)then
 do i=1,s1
  do j=1,s2
   MATMUL_x_c(i,j)=0.d0
   if(i<=s1/2)then
    do k=1,smiddle/2
      MATMUL_x_c(i,j)=MATMUL_x_c(i,j)+aa(i,k)*bb(k,j)
    enddo
   else
    do k=smiddle/2+1,smiddle
      MATMUL_x_c(i,j)=MATMUL_x_c(i,j)+aa(i,k)*bb(k,j)
    enddo
   endif
   enddo
  enddo
 return
 endif
 endif

 do i=1,s1
  do j=1,s2
    MATMUL_x_c(i,j)=0.d0
   do k=1,smiddle
     MATMUL_x_c(i,j)=MATMUL_x_c(i,j)+aa(i,k)*bb(k,j)
   enddo
  enddo
 enddo

return
end function

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

function MATMUL_x_r(aa,bb,ind,n,IdL,IdR,byblock,a_by_block)
implicit none
integer             :: ii,i,j,k,s1,s2,smiddle
real(8),intent(in)  :: aa(:,:),bb(:,:)
real(8)             :: MATMUL_x_r(size(aa(:,1)),size(bb(1,:)))
integer,optional    :: n,ind(:)
logical,optional    :: IdL,IdR,byblock,a_by_block

if(present(byblock))then
 if(byblock)then
   ii=size(aa(:,1))
   i =size(aa(1,:))
   if(i/=ii) stop 'MATMUL_x by block should be only for square matrices'
   if(mod(ii,2)>0) stop 'MATMUL_x by block should only be for even linear sizes'
   ii=ii/2
   MATMUL_x_r(1    :   ii, ii+1 : 2*ii)=0.
   MATMUL_x_r(ii+1 : 2*ii,    1 :   ii)=0.
   MATMUL_x_r(   1 :   ii,    1 :   ii)= MATMUL_(aa(    1:   ii,    1 :   ii),bb( 1 :  ii,    1 :  ii))
   MATMUL_x_r(ii+1 : 2*ii, ii+1 : 2*ii)= MATMUL_(aa(ii +1: 2*ii, ii+1 : 2*ii),bb( ii+1 :2*ii, ii+1 :2*ii))
   return
 endif
endif


!ind: matrice aa n a que les lignes dans ind qui different de l identite

s1=size(aa(:,1)); s2=size(bb(1,:)); smiddle=size(aa(1,:))
if(smiddle/=size(bb(:,1))) stop 'error MATMUL_x_r'

if(present(IdL))then
 if(testing)then
  if(abs(aa(1,2))>1.d-3) stop 'MATMUL_x_r bad diagonal matrix IdL'
 endif
 do i=1,size(bb(:,1))
  MATMUL_x_r(i,:)=bb(i,:)*aa(i,i)
 enddo
 return
endif

if(present(IdR))then
 if(testing)then
  if(abs(aa(1,2))>1.d-3) stop 'MATMUL_x_r bad diagonal matrix IdR'
 endif
 do i=1,size(bb(1,:))
  MATMUL_x_r(:,i)=bb(:,i)*aa(i,i)
 enddo
 return
endif

if(present(ind))then
 MATMUL_x_r=bb
 if(size(ind)/=n) stop 'error MATMUL_x_r bad ind shape'
 do ii=1,n
  i=ind(ii)
  do j=1,s2
   MATMUL_x_r(i,j)=0.
   do k=1,smiddle
    MATMUL_x_r(i,j)=MATMUL_x_r(i,j)+aa(i,k)*bb(k,j)
   enddo
  enddo
 enddo
return
endif

 if(present(a_by_block))then
 if(a_by_block)then
 do i=1,s1
  do j=1,s2
   MATMUL_x_r(i,j)=0.d0
   if(i<=s1/2)then
    do k=1,smiddle/2
      MATMUL_x_r(i,j)=MATMUL_x_r(i,j)+aa(i,k)*bb(k,j)
    enddo
   else
    do k=smiddle/2+1,smiddle
      MATMUL_x_r(i,j)=MATMUL_x_r(i,j)+aa(i,k)*bb(k,j)
    enddo
   endif
   enddo
  enddo
 return
 endif
 endif

do i=1,s1
 do j=1,s2
  MATMUL_x_r(i,j)=0.d0
  do k=1,smiddle
   MATMUL_x_r(i,j)=MATMUL_x_r(i,j)+aa(i,k)*bb(k,j)
  enddo
 enddo
enddo

return
end function

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
! 
!    subroutine inv_mat_small_numb_of_lines(mat,n,nind,ind,vind,test,messages)
!    implicit none
!     integer           :: nind,n,kk,kkk(nind+1),ii,jj
!     complex(8)        :: Diag(n,n),mat(n,n),vind(nind,n)
!     complex(8)        :: small(nind,nind),vectors(nind,nind),smallinv(nind,nind),vn(nind)
!     integer           :: ind(nind),i,j,k,l,m
!     real(8)           :: RWORK(3*nind),dnorm
!     complex(8)        :: WORK(3*nind),W(nind),mattest(n,n),DUMMY(1,1),matinv(n,n)
!     integer(4)        :: INFO
!     logical,optional  :: test,messages
! 
!        !matrix : Id + row vind(1,:) + row vind(2,:) + etc...
! 
!         !--------------------------------------------------------!
!         small=0.
!         if(present(messages))write(*,*) 'define matrix small'
!         l=0
!         do i=1,n
!          if(askinset(i,ind))then
!           l=l+1
!           k=0
!           do j=1,n
!            if(askinset(j,ind))then
!              k=k+1
!              if(i==j) then 
!               small(l,k)=1.d0+vind(l,j) 
!              else
!               small(l,k)=vind(l,j)
!              endif
!            endif
!           enddo
!          endif
!         enddo 
! 
!         if(present(messages))write(*,*) 'call ZGEEV'
!         call ZGEEV('N','V',nind,small,nind,W, DUMMY,1,vectors,nind,WORK,3*nind,RWORK,INFO)
!         if(present(messages))write(*,*) 'done, build large eigenvectors'
! 
!         !--------------------------------------------------------!
!         mat=0.;l=0
!         do i=1,n !eigenvectors
!          if(askinset(i,ind))then
!           l=l+1
!           k=0
!           do j=1,n
!            if(askinset(j,ind))then
!             k=k+1
!              mat(j,i)=vectors(k,l)
!            endif
!           enddo
!          else
!           kkk=ind_cycle(i,nind+1,n)
!           do ii=1,nind
!            do jj=1,nind
!             smallinv(ii,jj)=vind(ii,kkk(jj+1))
!            enddo
!            vn(ii)=-vind(ii,kkk(1))
!           enddo
!           call invmat_comp(nind,smallinv)
!           vn=MATMUL(smallinv,vn) 
!           mat(kkk(1),i)=1.d0
!           do ii=1,nind
!            mat(kkk(ii+1),i)=vn(ii) 
!           enddo
!          endif
!         enddo
!         do i=1,n
!          dnorm=sum(abs(mat(:,i))**2)
!          mat(:,i)=mat(:,i)/dnorm
!          if(present(messages)) write(*,*) 'norme vec i :', dnorm
!         enddo
! 
!         !--------------------------------------------------------!
!         if(present(test)) then    
!         do i=1,n !eigenvectors
!          do j=1,n
!           dnorm=abs(scalprod(mat(:,i),mat(:,j)))
!           if(dnorm>1.d-4.and.i/=j) then
!             write(*,*) 'non-orthogonal basis : ',i,j,dnorm
!           endif
!          enddo
!         enddo
!         endif
!         !--------------------------------------------------------!
! 
!         matinv=mat; call invmat(n,matinv)
!         if(present(messages))write(*,*) 'done, build diagonal matrix'
!         Diag=0.;l=0
!         do i=1,n
!           if(askinset(i,ind))then
!            l=l+1 
!            Diag(i,i)=W(l)
!           else
!            Diag(i,i)=1.d0
!           endif
!         enddo
! 
!        if(present(messages))write(*,*) 'done build unitary matrix'
!        if(present(test))then
!          mattest=Id(n)
!          l=0
!          do i=1,n 
!           if(askinset(i,ind))then
!            l=l+1
!            mattest(i,:)=mattest(i,:)+vind(l,:)
!           endif
!          enddo
!          write(*,*) 'test mult line invmat matrix, max DIFF : '
!          write(*,*) maxval(abs(mattest-MATMUL(mat,MATMUL(Diag,matinv))))
!        endif
! 
!        if(present(messages))write(*,*) 'build inverse diagonal matrix'
!        do i=1,n
!         Diag(i,i)=1.d0/Diag(i,i)
!        enddo
!        if(present(messages))write(*,*) 'final inverse matrix'
!        mat=MATMUL(mat,MATMUL(Diag,matinv))
!        if(present(messages))write(*,*) 'max DIFF : '
!        if(present(messages))write(*,*) maxval(abs(MATMUL(mattest,mat)-Id(n)))
!        if(present(messages))write(*,*) 'done....'
! 
!    return
!    end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! subroutine unitary_matrix(lsize,mat,unitary,vaps)
! implicit none
! integer                   :: lsize,i
! real(8)                   :: RWORK(3*lsize),W(lsize)
! complex(8)                :: WORK(3*lsize)
! complex(8),intent(inout)  :: mat(lsize,lsize)
! complex(8),optional       :: unitary(lsize,lsize)
! complex(8)                :: rrr
! integer(4)                :: INFO
! real(8),optional          :: vaps(lsize)
! 
!   if(testing) call check_hermitian('unitary_matrix, mat not hermitian',mat)
! 
!   if(present(unitary))then
!     if(lsize/=size(mat(:,1))) stop 'bad shape in unitary_matrix'
!     unitary=mat
!     call ZHEEV('V','U',lsize,unitary,lsize,W,WORK,3*lsize,RWORK,INFO)
!     do i=1,lsize
!      rrr=PHASE(unitary(1,i))
!      unitary(:,i)=unitary(:,i)/rrr
!     enddo
!   else
!     if(lsize/=size(mat(:,1)))stop 'bad shape in unitary_matrix'
!     call ZHEEV('V','U',lsize,mat,lsize,W,WORK,3*lsize,RWORK,INFO)
!     do i=1,lsize
!      rrr=PHASE(mat(1,i))
!      mat(:,i)=mat(:,i)/rrr
!     enddo
!   endif
! 
!   if(INFO/=0)then
!     write(*,*) 'BAD unitary_matrix , info = :', INFO
!     write(*,*) 'stop calculations...'
!     stop
!   endif
! 
!   if(present(vaps)) vaps=W
! 
! end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! subroutine print_eigenvalue__(mat)
! implicit none
! complex(8) :: mat(:,:)
! real(8)    :: vap(size(mat,1))
!  call eigenvalue_matrix(size(mat,1),mat,vap)
!  write(*,*) 'eigenvalues are : ' 
!  write(*,'(1000f10.4)') vap
! end subroutine
! 
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! subroutine eigenvalue_matrix__(lsize,mat,vaps)
! implicit none
! integer                   :: lsize,i
! real(8)                   :: RWORK(3*lsize)
! complex(8)                :: WORK(3*lsize)
! complex(8)                :: mat(lsize,lsize),temp(lsize,1),temp2(lsize,1)
! complex(8)                :: rrr
! integer                   :: INFO
! complex(8)                :: vaps(lsize)
!    call ZGEEV('N','N', lsize, mat, lsize, vaps, temp, 1, temp2, 1, WORK, 3*lsize, RWORK, INFO )
! end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! subroutine eigenvalue_matrix_(lsize,mat,vap)
! implicit none
! integer                   :: lsize,i
! real(8)                   :: RWORK(3*lsize)
! complex(8)                :: WORK(3*lsize)
! complex(8)                :: mat(lsize,lsize)
! complex(8)                :: rrr
! integer(4)                :: INFO
! real(8)                   :: vap(lsize),rr
! 
!    if(testing) then
!      rr=maxval(abs( mat - transpose(conjg(mat)) ))
!      if(rr>1.d-4) then
!        write(*,*) 'A-A^\dagger : ', rr
!        stop 'error eigenvalue matrix not sym'
!      endif
!    endif
! 
!    call ZHEEV('N','U',lsize,mat,lsize,vap,WORK,3*lsize,RWORK,INFO)
! 
!    if(INFO/=0)then
!      write(*,*) 'BAD eigenvalue calculation , info = :', INFO
!      write(*,*) 'stop calculations...'
!      stop
!    endif
! 
! end subroutine
! 
! subroutine eigenvalue_matrixr_(lsize,mat,vap)
! implicit none
! integer                   :: lsize,i
! real(8)                   :: RWORK(3*lsize)
! real(8)                   :: WORK(3*lsize)
! real(8)                   :: mat(lsize,lsize)
! real(8)                   :: rrr
! integer(4)                :: INFO
! real(8)                   :: vap(lsize),rr
! 
!    if(testing) then
!      rr=maxval(abs( mat - transpose(mat) ))
!      if(rr>1.d-4) then
!        write(*,*) 'A-A^\T : ', rr
!        stop 'error eigenvalue matrix not sym'
!      endif
!    endif
! 
!    call DSYEV('N','U',lsize,mat,lsize,vap,WORK,3*lsize,RWORK,INFO)
! 
!    if(INFO/=0)then
!      write(*,*) 'BAD eigenvalue calculation , info = :', INFO
!      write(*,*) 'stop calculations...'
!      stop
!    endif
! 
! end subroutine
! 
! 
! 
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

  !--------------------------------------------------!

subroutine eigenvector_matrix_c(lsize,mat,vaps,eigenvec)
implicit none
integer     :: lsize
real(8)     :: RWORK(3*lsize)
complex(8)  :: WORK(3*lsize)
complex(8)  :: mat(lsize,lsize)
complex(8)  :: eigenvec(lsize,lsize)
integer     :: INFO
! integer     :: order(lsize)
real(8)     :: vaps(lsize)

   if(testing)then
     if(maxval(abs( mat - transpose(conjg(mat)) ))>1.d-4)  then
        call write_array(mat,'matrix',unit=6,short=.true.)
        stop 'error : eigenvector_matrix_c : complex eigenvalue, should call eigenvector_matrix_c_'
     endif
   endif

   if(force_invmat_single_prec) then
     call go_for_single_prec  
     return
   endif

   if(.not.same_array(mat,eigenvec)) eigenvec=mat

   ! ebl: Removing GPU functionality
   ! if(lsize>n_cuda_rout.and.lsize<10000.and.use_cula_routines)then
   !   call cula_eigenvector_square_c(lsize,eigenvec,vaps,mat)
   !   call qsort_array(vaps,order)
   !   call qsort_adj_array(eigenvec,order)
   !   return
   ! endif

   call ZHEEV('V','U',lsize,eigenvec,lsize,vaps,WORK,3*lsize,RWORK,INFO)

   if(INFO/=0)then
     write(*,*) 'BAD eigenvalue calculation , info = :', INFO
     write(*,*) 'stop calculations...'
     stop
   endif
 
 contains

 !-----------------------!

 subroutine go_for_single_prec
 complex(4) :: mat_(lsize,lsize),eigenvec_(lsize,lsize)
 real(4)    :: vaps_(lsize)
  mat_=cmplx(mat, kind=4)
  call eigenvector_matrix(lsize,mat_,vaps_,eigenvec_)
  eigenvec=eigenvec_
  vaps=vaps_
 end subroutine

 !-----------------------!

end subroutine

  !--------------------------------------------------!

subroutine rearrange_columns_to_identity_c(lsize,mat,diagdens)
implicit none
integer                   :: lsize
complex(8)                :: diagdenstemp(lsize),diagdens(lsize),mat(lsize,lsize),mat_temp(lsize,lsize)
integer                   :: i
integer                   :: uu(1),uuu(lsize),uuu_not_placed(lsize),j,k

 uuu=0
 uuu_not_placed=0
 k=0

 do i=1,lsize
  uu=maxloc(abs(mat(:,i)))  
  j=uu(1)
  if(uuu(j)==0)then
    uuu(j)=i     
  else
    k=k+1
    uuu_not_placed(k)=i
  endif
 enddo
 
 k=0
 do i=1,lsize
  if(uuu(i)==0)then
    k=k+1
    uuu(i)=uuu_not_placed(k)
  endif
 enddo

 mat_temp=mat
 diagdenstemp=diagdens
 do i=1,lsize
   mat(:,i)=mat_temp(:,uuu(i))
   diagdens(i)=diagdenstemp(uuu(i))
 enddo

return
end subroutine

!--------------------------------------------------!

subroutine rearrange_columns_to_identity_r(lsize,mat,diagdens)
implicit none
integer                   :: lsize
real(8)                   :: diagdenstemp(lsize),diagdens(lsize),mat(lsize,lsize),mat_temp(lsize,lsize)
integer                   :: i
integer                   :: uu(1),uuu(lsize),uuu_not_placed(lsize),j,k

 uuu=0
 uuu_not_placed=0
 k=0

 do i=1,lsize
  uu=maxloc(abs(mat(:,i)))
  j=uu(1)
  if(uuu(j)==0)then
    uuu(j)=i
  else
    k=k+1
    uuu_not_placed(k)=i
  endif
 enddo

 k=0
 do i=1,lsize
  if(uuu(i)==0)then
    k=k+1
    uuu(i)=uuu_not_placed(k)
  endif
 enddo

 mat_temp=mat
 diagdenstemp=diagdens
 do i=1,lsize
   mat(:,i)=mat_temp(:,uuu(i))
   diagdens(i)=diagdenstemp(uuu(i))
 enddo

return
end subroutine

!--------------------------------------------------!

subroutine eigenvector_matrix_c_(lsize,mat,vaps,eigenvec,symmetric)
implicit none
integer                   :: lsize
real(8)                   :: RWORK(3*lsize)
complex(8)                :: WORK(3*lsize)
complex(8)                :: mat(lsize,lsize)
complex(8)                :: mat_temp(lsize,lsize)
complex(8)                :: eigenvec(lsize,lsize),temp(lsize,1)
integer                   :: INFO,n
complex(8)                :: vaps(lsize)
logical,optional          :: symmetric

   if(.not.present(symmetric))Then
     mat_temp=mat
     call ZGEEV('N','V', lsize, mat_temp, lsize, vaps, temp, 1, eigenvec, lsize, WORK, 3*lsize, RWORK, INFO )
   else
     !** SEigensystem diagonalizes a complex symmetric n-by-n matrix.
     !** Input: n, A = n-by-n matrix, complex symmetric
     !** (only the upper triangle of A needs to be filled).
     !** Output: d = vector of eigenvalues, U = transformation matrix
     !** these fulfill diag(d) = U A U^T = U A U^-1 with U U^T = 1.
     n=size(mat,1)
     mat_temp=mat
     call SEigensystem(n,mat_temp,n,vaps,eigenvec,n,1)
     write(*,*) 'CALL S_eigen, check U U^T : '
     write(*,*)  minval(abs( matmul(eigenvec,transpose(eigenvec))))
     write(*,*)  maxval(abs( matmul(eigenvec,transpose(eigenvec))))
     write(*,*) 'eigenvalues : ', vaps
     eigenvec=transpose(eigenvec)
     call write_array( matmul( matmul(eigenvec,mat),transpose(eigenvec)),' check U A U^T ', unit=6,short=.true.)
     call write_array( matmul( matmul(transpose(eigenvec),mat),eigenvec),' check U^T A U ', unit=6,short=.true.)
   endif

end subroutine

  !--------------------------------------------------!

subroutine eigenvector_matrix_cc(lsize,mat,vaps,eigenvec)
implicit none
integer                   :: lsize
real(4)                   :: RWORK(3*lsize),W(lsize)
complex(4)                :: WORK(3*lsize)
complex(4)                :: mat(lsize,lsize)
complex(4)                :: eigenvec(lsize,lsize)
integer                   :: INFO
! integer     :: order(lsize)
real(4)                   :: vaps(lsize)

   if(testing)then
    if(maxval(abs( mat - transpose(conjg(mat)) ))>1.d-4)  then
       call write_array(real(mat,kind=8),'matrix',unit=6,short=.true.)
       stop 'error : eigenvector_matrix_c : complex eigenvalue, should call eigenvector_matrix_c_'
    endif
   endif

   if(.not.same_array(mat,eigenvec))eigenvec=mat

   ! ebl: Removing GPU functionality
   ! if(lsize>n_cuda_rout.and.lsize<10000.and.use_cula_routines)then
   !   call cula_eigenvector_square_cs(lsize,eigenvec,vaps,mat)
   !   call qsort_array(vaps,order)
   !   call qsort_adj_array(eigenvec,order)
   !   return
   ! endif

   call CHEEV('V','U',lsize,eigenvec,lsize,W,WORK,3*lsize,RWORK,INFO)

   if(INFO/=0)then
    write(*,*) 'BAD eigenvalue calculation , info = :', INFO
    write(*,*) 'stop calculations...'
    stop
   endif

end subroutine

  !--------------------------------------------------!

subroutine eigenvector_matrix_cc_(lsize,mat,vaps,eigenvec)
implicit none
integer                   :: lsize
real(4)                   :: RWORK(3*lsize)
complex(4)                :: WORK(3*lsize)
complex(4)                :: mat_temp(lsize,lsize),mat(lsize,lsize),temp(lsize,1)
complex(4)                :: eigenvec(lsize,lsize)
integer                   :: INFO
complex(4)                :: vaps(lsize)
   mat_temp=mat
   call CGEEV('N','V', lsize, mat_temp, lsize, vaps, temp, 1, eigenvec, lsize, WORK, 3*lsize, RWORK, INFO )
end subroutine

  !--------------------------------------------------!

subroutine eigenvector_matrix_r(lsize,mat,vaps,eigenvec)
implicit none
integer                   :: lsize
real(8)                   :: WORK(3*lsize)
real(8)                   :: mat(lsize,lsize)
real(8)                   :: eigenvec(lsize,lsize)
integer                   :: INFO
! integer                   :: order(lsize)
real(8)                   :: vaps(lsize)

   if(lsize<1) stop 'error eigenvector_matrix, 0 dim'

   if(testing)then
      if(maxval(abs( mat - transpose((mat)) ))>1.d-4)write(*,*) 'error eigenvector_matrix_r, matrix should be symmetric'
   endif

   if(force_invmat_single_prec) then
    call go_for_single_prec
    return
   endif

   if(.not.same_array(mat,eigenvec))then
    eigenvec=mat
   endif

   ! ebl: Removing GPU functionality
   ! if(lsize>n_cuda_rout.and.lsize<10000.and.use_cula_routines)then
   !   write(*,*) 'eigenvector_matrix, calling cula'
   !   call cula_eigenvector_square_d(lsize,eigenvec,vaps,mat)
   !   call qsort_array(vaps,order)
   !   call qsort_adj_array(eigenvec,order)
   !   return
   ! endif

   call DSYEV('V','U',lsize,eigenvec,lsize,vaps,WORK,3*lsize,INFO)

   if(INFO/=0)then
    write(*,*) 'BAD eigenvalue calculation , info = :', INFO
    write(*,*) 'stop calculations...'
    stop
   endif

 contains

 !-----------------------!

 subroutine go_for_single_prec
 real(4)    :: mat_(lsize,lsize),eigenvec_(lsize,lsize)
 real(4)    :: vaps_(lsize)
  mat_=real(mat, kind=4)
  call eigenvector_matrix(lsize,mat_,vaps_,eigenvec_)
  eigenvec=eigenvec_
  vaps=vaps_
 end subroutine

 !-----------------------!

end subroutine


  !--------------------------------------------------!

subroutine eigenvector_matrix_rc(lsize,mat,vaps,eigenvec)
implicit none
integer                   :: lsize
real(8)                   :: WORK(3*lsize),RWORK(3*lsize)
real(8)                   :: mat(lsize,lsize)
real(8)                   :: eigenvec(lsize,lsize)
integer                   :: INFO
complex(8)                :: vaps(lsize)
real(8)                   :: temp(1,lsize),vl(lsize),vr(lsize)

   if(.not.same_array(mat,eigenvec))eigenvec=mat
   call DGEEV('N','V', lsize, mat , lsize, vl, vr , temp, 1, eigenvec, lsize, WORK, 3*lsize, RWORK, INFO )
   vaps=vl+imi*vr

end subroutine

  !--------------------------------------------------!

subroutine eigenvector_matrix_rr(lsize,mat,vaps,eigenvec)
implicit none
integer                   :: lsize
real(4)                   :: WORK(3*lsize)
real(4)                   :: mat(lsize,lsize)
real(4)                   :: eigenvec(lsize,lsize)
integer                   :: INFO
! integer                   :: order(lsize)
real(4)                   :: vaps(lsize)

   if(lsize<1) stop 'error eigenvector_matrix, 0 dim'

   if(testing)then
     if(maxval(abs( mat - transpose((mat)) ))>1.d-4) write(*,*) 'error eigenvector_matrix_rr, matrix should be symmetric'
   endif

   if(.not.same_array(mat,eigenvec)) eigenvec=mat

   ! ebl: Removing GPU functionality
   ! if(lsize>n_cuda_rout.and.lsize<10000.and.use_cula_routines)then
   !   call cula_eigenvector_square_r(lsize,eigenvec,vaps,mat)
   !   call qsort_array(vaps,order)
   !   call qsort_adj_array(eigenvec,order)
   !   return 
   ! endif

   call SSYEV('V','U',lsize,eigenvec,lsize,vaps,WORK,3*lsize,INFO)

   if(INFO/=0)then
    write(*,*) 'BAD eigenvalue calculation , info = :', INFO
    write(*,*) 'stop calculations...'
    stop
   endif

end subroutine


  !--------------------------------------------------!

subroutine eigenvector_matrix_rrc(lsize,mat,vaps,eigenvec)
implicit none
integer                   :: lsize
real(4)                   :: WORK(3*lsize),RWORK(3*lsize)
real(4)                   :: mat(lsize,lsize)
real(4)                   :: eigenvec(lsize,lsize)
integer                   :: INFO
complex(4)                :: vaps(lsize)
real(4)                   :: temp(lsize,1),vl(lsize),vr(lsize)
   if(.not.same_array(mat,eigenvec))eigenvec=mat
   call SGEEV('N','V', lsize, mat, lsize, vl, vr, temp, 1, eigenvec, lsize, WORK, 3*lsize, RWORK, INFO )
   vaps=cmplx(vl+imi*vr, kind=4)
end subroutine

  !--------------------------------------------------!

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
! 
! function expMat(lsize,mat,test)
! implicit none
! integer(4)                :: i,j,k,l,m,lsize
! real(8)                   :: RWORK(3*lsize),W(lsize),ddnorm
! complex(8)                :: WORK(3*lsize),BB(lsize,lsize)
! complex(8),intent(in)     :: mat(:,:)
! complex(8)                :: expMat(lsize,lsize),dnorm
! integer(4)                :: INFO
! logical,optional          :: test
! 
!   expMat=mat
! 
!   if(size(mat(1,:))/=lsize) stop 'error size matrix expMat'
! 
!   call ZHEEV('V','U',lsize,expMat,lsize,W,WORK,3*lsize,RWORK,INFO)
! 
!   if(INFO/=0)then
!    write(*,*) 'BAD expMAT , info = :', INFO
!    write(*,*) 'stop calculations...'
!    stop
!   endif   
! 
!   if(present(test))then
!    call check_hermitian('expMat error, mat not hermitian',mat)
!    BB=0.d0 
!    do i=1,lsize
!     BB(i,i)=W(i)
!    enddo
!    if(size(mat(:,1))/=size(mat(1,:))) write(*,*) 'mat not square'
!    if(size(mat(:,1))/=lsize) write(*,*) 'wrong size exp mat'
!    ddnorm=maxval(abs(mat-MATMUL(expMat,MATMUL(BB,TRANSPOSE(CONJG(expMat))))))
!    if(ddnorm>1.d-2) then
!     write(*,*) 'max DIFF : '
!     write(*,*) ddnorm
!     write(*,*) 'max value exp mat : ',maxval(abs(expMat))
!     write(*,*) 'max vlaue mat     : ',maxval(abs(mat))
!     write(*,*) 'INFO : ', INFO
!     stop 'error expMat'
!    endif
!   endif
! 
!   BB=0.
!   do i=1,lsize
!    BB(i,i)=exp(W(i))
!   enddo
!   expMat=MATMUL(expmat,(MATMUL_x(BB,TRANSPOSE(CONJG(expMat)),IdL=.true.)))
! 
! return
! end function
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! subroutine inverse_sym_mat(mat,lsize,test,c)
! implicit none
! integer                   :: i,j,k,l,m,lsize
! real(8)                   :: RWORK(3*lsize),W(lsize),ddnorm
! complex(8)                :: WORK(3*lsize),BB(lsize,lsize)
! complex(8),intent(inout)  :: mat(lsize,lsize)
! complex(8)                :: mat_back(lsize,lsize),dnorm
! integer(4)                :: INFO
! logical,optional          :: test
! integer,optional          :: c
! 
!  if(present(c))then
!   call randomize_matrix(mat,amp=invmat_error,flag=.true.,kk2=c)
!  else
!   call randomize_matrix(mat,amp=invmat_error)
!  endif
! 
!   mat_back=mat
!   call ZHEEV('V','U',lsize,mat_back,lsize,W,WORK,3*lsize,RWORK,INFO)
!   if(INFO/=0)then
!    write(*,*) 'BAD eigenvect. in inverse_mat_sym , info = :', INFO
!    write(*,*) 'stop calculations...'
!    stop
!   endif
! 
!   if(present(test))then
!    write(*,*) 'mat - mat^dagger : ', maxval(abs(mat-conjg(transpose(mat))))
!    write(*,*) 'max el           : ', maxval(abs(mat))
!    BB=0.d0
!    do i=1,lsize
!     BB(i,i)=W(i)
!    enddo
!    if(size(mat(:,1))/=size(mat(1,:))) stop 'inverse_sym_mat : mat not square'
!    ddnorm=maxval(abs(mat-MATMUL(mat_back,MATMUL(BB,TRANSPOSE(CONJG(mat_back))))))
!    if(ddnorm>1.d-2) then
!     write(*,*) 'max DIFF           : ', ddnorm
!     write(*,*) 'max element in mat : ', maxval(abs(mat))
!     stop 'error mat_back: bad decomposition'
!    endif
!   endif
! 
!   BB=0.d0
!   do i=1,lsize
!    if(abs(W(i))>1.d-25) then
!     BB(i,i)=1.d0/W(i)
!    else 
!     if(messages3) write(*,*) 'DANGER In inverse_sym_mat, very large inverse element'
!     BB(i,i)=1.d25
!    endif
!   enddo
! 
!   mat_back=MATMUL(mat_back,MATMUL(BB,TRANSPOSE(CONJG(mat_back))))
!   
!   if(present(test))then
!     ddnorm=maxval(abs(Id(lsize)-MATMUL(mat,mat_back)))
!      if(ddnorm>1.d-2) then
!         mat=MATMUL(mat,mat_back)
!         write(*,*) 'error,inverse_sym_mat, divergence, here we have (mat*mat^-1)-Id :'
!         write(*,*)  ddnorm
!         write(*,*) 'max element :', maxval(abs(mat)),maxval(abs(mat_back)),maxval(abs(BB))
!         do i=1,min(lsize,3)
!          write(*,*) mat(i,:)
!          write(*,*) 
!         enddo
!         stop 'bad inverse...'
!      endif
!   endif
! 
!   mat=mat_back
! 
! return
! end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       subroutine update_inv(n,invmat,v1,v2,imouv,jmouv,flagerror,temdet,temdet2,onlydet,only_i_term)
!       implicit none
!       integer,intent(in)    :: imouv,jmouv,n
!       integer,intent(out)   :: flagerror
!       integer               :: k
!       real(8)               :: error
!       complex(8)            :: invmat(n,n),tem(n,n),mqj(n),v1(n),vect(n)
!       complex(8)            :: v2(n)
!       complex(8)            :: piq(n)
!       real(8)               :: temdet2,rrr
!       complex(8)            :: temdet,garb,garba,garbb,garbf
!       complex(8)            :: tempscal
!       logical,optional      :: onlydet
!       integer,optional      :: only_i_term
! 
!       error=1.d-9
! 
!       flagerror=0
!       if(n==0) stop '0 dimension in update_inv' 
!       !-------------------------------------!
! 
!       if(imouv==0) temdet=sum(invmat(jmouv,:)*v2(:))
!       if(jmouv==0) temdet=sum(invmat(:,imouv)*v1(:))
! 
!       if(jmouv>0.and.imouv>0) then
!         do k=1,n
!          piq(k)=sum(v1(:)*invmat(:,k))
!         enddo
!         tempscal = v1(jmouv)
!         garb     = sum(piq(:)*v2(:))
!         garb     = invmat(jmouv,imouv) * (tempscal-garb)
!         garba    = sum( v1(:) * invmat(:,imouv) )
!         garbb    = sum( invmat(jmouv,:) * v2(:) )
!         temdet   = (garb + garba*garbb)
!       endif
! 
!       rrr=abs(temdet) 
!       if(rrr>1.d-9) then
!        temdet2=rrr
!        temdet=temdet/rrr
!       else
!        temdet2=rrr
!        temdet=1. 
!       endif
! 
!       if(present(onlydet)) return
! 
! 
!       if(present(only_i_term))then
!       !-------------------------------------!
!       if(imouv==0) then
!         mqj(:)=invmat(:,only_i_term)*v2(only_i_term)
!         if(abs(mqj(jmouv))<error) then
!           if(messages) then
!            write(*,*) 'error in update_inv'
!            write(*,*) 'mqj, jmouv : ', mqj(jmouv),jmouv
!           endif
!          flagerror=1
!          return
!         endif
!         vect(:)=invmat(jmouv,:)/mqj(jmouv)
!         do k=1,n
!          invmat(k,:)=invmat(k,:)-mqj(k)*vect(:)
!         enddo
!         invmat(jmouv,:)=vect(:)
!       endif
!       !-------------------------------------!
!       if(jmouv==0)then
!          mqj(:)=v1(only_i_term)*invmat(only_i_term,:)
!          if(abs(mqj(imouv))<error) then
!           if(messages) then
!            write(*,*) 'error in update_inv'
!            write(*,*) 'mqj, imouv : ', mqj(imouv),imouv
!           endif
!           flagerror=1
!           return
!          endif
!         vect(:)=invmat(:,imouv)/mqj(imouv)
!         do k=1,n
!           invmat(k,:)=invmat(k,:)-vect(k)*mqj(:)
!         enddo
!         invmat(:,imouv)=vect(:)
!        endif
!       !-------------------------------------!
! 
! 
!       return
! 
!       endif
! 
!       !-------------------------------------!
!       if(imouv==0) then
!         mqj(:)=0.d0
!         do k=1,n
!          mqj(:)=invmat(:,k)*v2(k)+mqj(:)
!         enddo
!         if(abs(mqj(jmouv))<error) then
!           if(messages) then
!            write(*,*) 'error in update_inv'
!            write(*,*) 'mqj, jmouv : ', mqj(jmouv),jmouv
!           endif
!          flagerror=1
!          return
!         endif
!         vect(:)=invmat(jmouv,:)/mqj(jmouv)
!         do k=1,n
!          invmat(k,:)=invmat(k,:)-mqj(k)*vect(:)
!         enddo
!         invmat(jmouv,:)=vect(:)
!       endif
!       !-------------------------------------!
!       if(jmouv==0)then
!         mqj(:)=0.d0
!         do k=1,n
!          mqj(:)=v1(k)*invmat(k,:)+mqj(:)
!         enddo
!          if(abs(mqj(imouv))<error) then
!           if(messages) then
!            write(*,*) 'error in update_inv'
!            write(*,*) 'mqj, imouv : ', mqj(imouv),imouv
!           endif
!           flagerror=1 
!           return
!          endif
!         vect(:)=invmat(:,imouv)/mqj(imouv)
!         do k=1,n
!           invmat(k,:)=invmat(k,:)-vect(k)*mqj(:)
!         enddo
!         invmat(:,imouv)=vect(:)
!        endif
!       !-------------------------------------!
!        if(imouv>0.and.jmouv>0) then
!          if(abs(piq(imouv))<error)then
!           flagerror=1
!           return
!          endif
!          vect(:)=invmat(:,imouv)/piq(imouv)
!          do k=1,n
!           tem(k,:)=invmat(k,:)-vect(k)* piq(:)
!          enddo
!          tem(:,imouv)=vect(:)
!          mqj(:)=0.d0
!          do k=1,n
!           mqj(:)=tem(:,k)*v2(k)+mqj(:)
!          enddo
!          if(abs(mqj(jmouv))<error) then
!           if(messages) write(*,*) 'error in update_inv'
!           flagerror=1
!           return
!          endif
!          vect(:)=tem(jmouv,:)/mqj(jmouv)
!          do k=1,n
!           invmat(k,:)=tem(k,:)-vect(:)*mqj(k)
!          enddo
!          invmat(jmouv,:)=vect(:)
!        endif
!       !-------------------------------------!
! 
!       return
!       end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!  subroutine addId(mat,n)
!  implicit none
!  integer    :: n
!  complex(8) :: mat(n,n)
!  integer    :: i
!    do i=1,n
!     mat(i,i)=mat(i,i)+1.
!    enddo
!  end subroutine
! 
!    !-------!
! 
 function Id(n)
 implicit none
 integer :: n,i
 real(8)  :: Id(n,n)
   Id=0.
   do i=1,n
    Id(i,i)=1.
   enddo
 end function
! 
!    !-------!
! 
!  function Idc(n)
!  implicit none
!  integer    :: n,i
!  complex(8) :: Idc(n,n)
!    Idc=0.d0
!    do i=1,n
!     Idc(i,i)=1.
!    enddo
!  end function
! 
   !-------!

 function bande_matr(n,bande)
 implicit none
 integer :: n,i
 real(8) :: bande(n),bande_matr(n,n)
  bande_matr=0
  do i=1,n
    bande_matr(i,i)=bande(i)
  enddo
 end function

   !-------!

 function bande_matc(n,bande)
 implicit none
 integer    :: n,i
 complex(8) :: bande(n),bande_matc(n,n)
  bande_matc=0
  do i=1,n
    bande_matc(i,i)=bande(i)
  enddo
 end function

   !-------!

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
! 
!     subroutine decomposemat(matrice,base)
!     implicit none
!     real(8),intent(inout) :: matrice(:,:),base(:,:)
!     integer :: sizemat
!      sizemat=size(matrice(:,1))
!      if(sizemat/=size(matrice(1,:))) &
!          & stop 'error decomposemat : matrice pas carree'
!      matrice=transpose(base)
!      call invmat_real(n=sizemat,mat=matrice)
!     end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!     subroutine symmetrize_mat_r(matrice)
!     implicit none
!     real(8),intent(inout) :: matrice(:,:)
!     integer :: sizemat,i,j
!      sizemat=size(matrice(:,1))
!      if(sizemat/=size(matrice(1,:))) &
!          & stop 'error sym_mat : matrice pas carree'
!      do i=1,sizemat
!       do j=i+1,sizemat
!        matrice(j,i)=matrice(i,j)
!       enddo
!      enddo
!     return
!     end subroutine
! 
!       !--------------------------!
! 
!     subroutine symmetrize_mat_c(matrice)
!     implicit none
!     complex(8),intent(inout) :: matrice(:,:)
!     integer :: sizemat,i,j
!      sizemat=size(matrice(:,1))
!      if(sizemat/=size(matrice(1,:))) &
!          & stop 'error sym_mat : matrice pas carree'
!      do i=1,sizemat
!       do j=i+1,sizemat
!        matrice(j,i)=conjg(matrice(i,j))
!       enddo
!      enddo
!     return
!     end subroutine
! 
! 
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

     subroutine decomposevec(coef,vec,base)
     implicit none
     real(8),intent(inout)   :: vec(:),coef(:)
     real(8)                 :: matrice(size(vec),size(vec)),base(:,:)
     integer                 :: sizemat
      sizemat=size(vec(:))
      if(sizemat/=size(matrice(1,:))) stop 'error decomposevec : dimensions'
      matrice=transpose(base)
      call invmat(sizemat,matrice)
      coef=matmul(matrice,vec)
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
! 
!       subroutine rescale_c(detphi,detphi2,dpdet)
!       implicit none
!        integer,intent(inout)    :: dpdet
!        real(8),intent(inout)     :: detphi
!        complex(8),intent(inout) :: detphi2
!        real(8)                   :: rdet
! 
!         rdet=abs(detphi2)
!         if(rdet>1.d-8)then
!          detphi2=detphi2/rdet
!          detphi=detphi*rdet
!         else
!          detphi2=1.
!          detphi=detphi*rdet 
!         endif         
!         
!         if(detphi==0.) goto 112
!   97    continue
!         if (abs(detphi)<=1.d0) goto 102
!         detphi = detphi/16.d0
!         dpdet = dpdet + 4
!         goto 97
!   102   continue
!         if (abs(detphi)>=(1.d0/16.d0)) goto 112
!         detphi = detphi * 16.d0
!         dpdet = dpdet - 4
!         goto 102
!   112   continue
!       end subroutine
! 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
! 
!       subroutine rescale_c3(detphi,detphi2,dpdet)
!       implicit none
!        integer,intent(inout)     :: dpdet
!        real(4),intent(inout)     :: detphi
!        complex(4),intent(inout)  :: detphi2
!        real(4)                   :: rdet
!         rdet=abs(detphi2)
!         if(rdet>1.d-8)then
!          detphi2=detphi2/rdet
!          detphi=detphi*rdet
!         else
!          detphi2=1.
!          detphi=detphi*rdet
!         endif
!         if(detphi==0.) goto 112
!   97    continue
!         if (abs(detphi)<=1.0) goto 102
!         detphi = detphi/16.0
!         dpdet = dpdet + 4
!         goto 97
!   102   continue
!         if (abs(detphi)>=(1.0/16.0)) goto 112
!         detphi = detphi * 16.0
!         dpdet = dpdet - 4
!         goto 102
!   112   continue
!       end subroutine
! 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
! 
!       subroutine rescale_r3(detphi,dpdet)
!       implicit none
!       integer,intent(inout) :: dpdet
!       real(4),intent(inout)  :: detphi
!         if(detphi==0.) goto 112
!   97    continue
!         if (abs(detphi)<=1.0) goto 102
!         detphi = detphi/16.0
!         dpdet = dpdet + 4
!         goto 97
!   102   continue
!         if (abs(detphi)>=(1.0/16.0)) goto 112
!         detphi = detphi * 16.0
!         dpdet = dpdet - 4
!         goto 102
!   112   continue
!       end subroutine
! 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
! 
!       subroutine rescale_r(detphi,dpdet)
!       implicit none
!       integer,intent(inout) :: dpdet
!       real(8),intent(inout)  :: detphi
!         if(detphi==0.) goto 112
!   97    continue
!         if (abs(detphi)<=1.d0) goto 102
!         detphi = detphi/16.d0
!         dpdet = dpdet + 4
!         goto 97
!   102   continue
!         if (abs(detphi)>=(1.d0/16.d0)) goto 112
!         detphi = detphi * 16.d0
!         dpdet = dpdet - 4
!         goto 102
!   112   continue
!       end subroutine
! 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
! 
!       subroutine rescale_c2(detphi,detphi2,dpdet)
!       implicit none
!        integer,intent(inout)    :: dpdet
!        real(16),intent(inout)   :: detphi
!        complex(16),intent(inout):: detphi2
!        real(16)                 :: rdet
! 
!         rdet=abs(detphi2)
!         if(rdet>1.d-22)then
!          detphi2=detphi2/rdet
!          detphi=detphi*rdet
!         else
!          detphi2=1.d0
!          detphi=detphi*rdet
!         endif
! 
!         if(abs(detphi)<1.d-24) goto 112
!   97    continue
!         if (abs(detphi)<=1.d0) goto 102
!         detphi = detphi/16.d0
!         dpdet = dpdet + 4
!         goto 97
!   102   continue
!         if (abs(detphi)>=(1.d0/16.d0)) goto 112
!         detphi = detphi * 16.d0
!         dpdet = dpdet - 4
!         goto 102
!   112   continue
! 
!       end subroutine
! 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
! 
!       subroutine rescale_r2(detphi,dpdet)
!       implicit none
!       integer,intent(inout) :: dpdet
!       real(16),intent(inout):: detphi
!         if(abs(detphi)<1.d-22) goto 112
!   97    continue
!         if (abs(detphi)<=1.d0) goto 102
!         detphi = detphi/16.d0
!         dpdet = dpdet + 4
!         goto 97
!   102   continue
!         if (abs(detphi)>=(1.d0/16.d0)) goto 112
!         detphi = detphi * 16.d0
!         dpdet = dpdet - 4
!         goto 102
!   112   continue
!       end subroutine
! 
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
! 
!        !---------------------!
! 
!       function invmat_(n,mat)
!       implicit none
!        complex(8) :: mat(n,n)
!        complex(8) :: invmat_(n,n)
!        integer :: n
!         invmat_=mat
!         call invmat(n,invmat_)
!       end function
! 
!        !---------------------!
! 
!       subroutine invmat_comp(n,mat,det2b,detb,pdetb,check_nan,c,block_matrix,diagmat)
!       implicit none
!       integer                   :: pdet,ppdet,n,nnn,i,j
!       real(8)                   :: big,det,ddet,rcond
!       complex(8),intent(inout)  :: mat(:,:)
!       real(8),optional          :: detb
!       integer,optional          :: pdetb
!       complex(8),optional       :: det2b
!       complex(8)                :: det2,ddet2,deti(2),q
!       complex(8)                :: WORK(2*size(mat,1)),WORK2(size(mat,1))
!       integer(4)                :: INFO,INFO2
!       integer                   :: piv(size(mat,1))
!       logical,optional          :: check_nan,block_matrix,diagmat
!       integer,optional          :: c
! 
!      if(present(diagmat))then
!       if(diagmat)then
!        do i=1,size(mat,1) 
!         if(abs(mat(i,i))>1.d-15)then
!          mat(i,i)=1.d0/mat(i,i)
!         else
!          mat(i,i)=1.d15
!         endif
!        enddo
!        return
!       endif
!      endif
! 
!      if(mat_3_3_c(mat,detb,det2b,pdetb)) return
! 
!      if(size(mat,1)/=size(mat,2)) then
!        write(*,*) 'DIMENSION OF MATRIX : ', shape(mat)
!        call create_seg_fault
!        stop 'error invmat_comp try to inverse rectangular matrix'
!      endif
!  
!      if(present(det2b)) det2b=0.; if(present(detb))  detb=0.; if(present(pdetb)) pdetb=0
!      if(present(check_nan)) call erase_divergence(mat)
!      call force_real_inv
! 
!      if(present(c))then
!        call randomize_matrix(mat,amp=invmat_error,flag=.true.,kk2=c)
!      else
!        call randomize_matrix(mat,amp=invmat_error,flag=.true.)
!      endif
! 
!       if(present(block_matrix))then
!        if(block_matrix)then
!          nnn=n/2
!          if(mod(n,2)>0) stop 'error invmat_comp block_matrix, but linear size is odd'
!          if(present(detb))then
!            call invit(1,nnn,det,det2,pdet)
!            call invit(nnn+1,n,ddet,ddet2,ppdet)
!            mat(1:nnn,nnn+1:n)=0.
!            mat(nnn+1:n,1:nnn)=0.
!            pdetb = pdet + ppdet
!            detb  = det  * ddet
!            det2b = det2 * ddet2
!          else
!            call invit(1,nnn)
!            call invit(nnn+1,n)
!            mat(1:nnn,nnn+1:n)=0.
!            mat(nnn+1:n,1:nnn)=0.
!          endif
!          goto 35
!        endif
!       endif
! 
!       call invit(1,n,detb,det2b,pdetb)
! 35    continue
! 
!       if(present(check_nan)) call erase_divergence(mat)
! 
!       return
! 
!       contains
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!    subroutine zgeco_zgedi(i1,i2)
!    implicit none
!    integer :: i1,i2,nnn
!      nnn=i2-i1+1
!      call geco__(nnn,mat(i1:i2,i1:i2),piv(i1:i2))
!      call gedi__(nnn,mat(i1:i2,i1:i2),piv(i1:i2),deti)
!    end subroutine
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!     subroutine invit(i1,i2,det,det2,pdet)
!     implicit none
!     integer                   :: i1,i2
!     real(8),optional          :: det
!     integer,optional          :: pdet
!     complex(8),optional       :: det2
!     integer                   :: nn
! 
!        nn = i2-i1+1
!   
!        if(mat_3_3_c(mat(i1:i2,i1:i2),det,det2,pdet)) return
! 
!        if(force_invmat_single_prec)then
!         call force_sing_prec(mat(i1:i2,i1:i2),det,det2,pdet)
!         return 
!        endif
! 
!        if(present(det).or.flag_force_invmat_lapack)then
!         call zgeco_zgedi(i1,i2)
!         call get_det_from_zgeco(deti,det,det2,pdet)
!        else
!         if(nn>n_openmp_rout.and.use_openmp_invmat)then
!           call invert_openmp(mat(i1:i2,i1:i2)) 
!           return
!         endif
!         ! ebl: Removing GPU functionality
!         ! if(use_cuda_routines.and.nn>n_cuda_rout)then
!         !   call diago_cuda_it_c(nn,mat(i1:i2,i1:i2))
!         !   return
!         ! else
!          if(.not.flag_use_invmat_jordan)then
!            call invmat_comp2(n,mat(i1:i2,i1:i2))
!          else
!            call invmat_jordan(n,mat(i1:i2,i1:i2))
!          endif
!         ! endif
!        endif
! 
!     end subroutine
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!     subroutine force_sing_prec(mat,det,det2,pdet)
!      complex(8)                :: mat(:,:)
!      complex(4)                :: mat_(size(mat,1),size(mat,2))
!      real(8),optional          :: det
!      integer,optional          :: pdet
!      complex(8),optional       :: det2
!      real(4)                   :: det_
!      complex(4)                :: det2_
!      integer                   :: pdet_
!      if(present(det))then 
!       mat_=mat
!       call invmat_comps(size(mat,1),mat_,det2_,det_,pdet_)
!       mat=mat_
!       pdet=pdet_
!       det2=det2_
!       det=det_
!      else
!       mat_=mat
!       call invmat_comps(size(mat_,1),mat_)
!       mat=mat_
!      endif
!     end subroutine
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!     subroutine force_real_inv
!      real(8) :: matr(n,n)
!      if(testing)then
!        if(maxval(abs(real(mat)))>1.d-9)then
!        if(maxval(abs(aimag(mat))/maxval(abs(real(mat))))<1.d-4)then
!         matr=mat
!         if(messages3) &
!      &  write(*,*) 'DANGER: invmat_comp, &
!      &              diagonalize real matrix, too small complex part'
!         if(present(det2b)) then
!         call invmat_real(n,matr,det2b,detb,pdetb)
!         else
!          call invmat_real(n,matr)
!         endif
!         mat=matr
!         return
!        endif
!        endif
!       endif
!     end subroutine
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!       end subroutine

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

!       subroutine invmat_comps(n,mat,det2b,detb,pdetb,check_nan,c,block_matrix,diagmat)
!       implicit none
!       integer                   :: pdet,ppdet,n,nnn,i,j
!       real(4)                   :: big,det,ddet
!       complex(4),intent(inout)  :: mat(:,:)
!       real(4),optional          :: detb
!       integer,optional          :: pdetb
!       complex(4),optional       :: det2b
!       complex(4)                :: det2,ddet2,deti(2),q
!       integer(4)                :: INFO,INFO2
!       integer                   :: piv(size(mat,1))
!       logical,optional          :: check_nan,block_matrix,diagmat
!       integer,optional          :: c
! 
!      if(present(diagmat))then
!       if(diagmat)then
!        do i=1,size(mat,1) 
!         if(abs(mat(i,i))>1.d-15)then
!          mat(i,i)=1.d0/mat(i,i)
!         else
!          mat(i,i)=1.d15
!         endif
!        enddo
!        return
!       endif
!      endif
! 
!      if(size(mat,1)/=size(mat,2)) then
!        write(*,*) 'DIMENSION OF MATRIX : ', shape(mat)
!        call create_seg_fault
!        stop 'error invmat_comp try to inverse rectangular matrix'
!      endif
! 
!      if(mat_3_3_cs(mat,detb,det2b,pdetb)) return
!  
!      if(present(det2b)) det2b=0.; if(present(detb))  detb=0.; if(present(pdetb)) pdetb=0
! 
!       if(present(block_matrix))then
!        if(block_matrix)then
!          nnn=n/2
!          if(mod(n,2)>0) stop 'error invmat_comp block_matrix, but linear size is odd'
!          if(present(detb))then
!            call invit(1,nnn,det,det2,pdet)
!            call invit(nnn+1,n,ddet,ddet2,ppdet)
!            mat(1:nnn,nnn+1:n)=0.
!            mat(nnn+1:n,1:nnn)=0.
!            pdetb = pdet + ppdet
!            detb  = det  * ddet
!            det2b = det2 * ddet2
!          else
!            call invit(1,nnn)
!            call invit(nnn+1,n)
!            mat(1:nnn,nnn+1:n)=0.
!            mat(nnn+1:n,1:nnn)=0.
!          endif
!          goto 35
!        endif
!       endif
! 
!       call invit(1,n,detb,det2b,pdetb)
! 35    continue
! 
!       return
! 
!       contains
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!     subroutine zgeco_zgedi(i1,i2)
!     implicit none
!     integer :: i1,i2,nnn
!       nnn=i2-i1+1
!      call geco__(nnn,mat(i1:i2,i1:i2),piv(i1:i2))
!      call gedi__(nnn,mat(i1:i2,i1:i2),piv(i1:i2),deti)
!     end subroutine
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!     subroutine invit(i1,i2,det,det2,pdet)
!     implicit none
!     integer :: i1,i2
!     real(4),optional          :: det
!     integer,optional          :: pdet
!     complex(4),optional       :: det2
!        if(mat_3_3_cs(mat(i1:i2,i1:i2),det,det2,pdet)) return
!        if(present(det).or.flag_force_invmat_lapack)then
!         call zgeco_zgedi(i1,i2)
!         call get_det_from_zgeco(deti,det,det2,pdet)
!        else
!         if(n>n_openmp_rout.and.use_openmp_invmat)then
!           call invert_openmp(mat(i1:i2,i1:i2)) 
!           return         
!         endif
!           call invmat_jordan(n,mat(i1:i2,i1:i2))
!        endif
!     end subroutine
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!       end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       subroutine invmat_comp_quad(n,mat,det2b,detb,pdetb,check_nan,c,block_matrix)
!       implicit none
!       real(16),optional          :: detb
!       integer,optional           :: pdetb
!       complex(16),optional       :: det2b
!       integer                    :: ppdet,pdet,n,i,j,piv(n)
!       real(16)                   :: big,det,ddet
!       complex(16),intent(inout)  :: mat(n,n)
!       complex(16)                :: ddet2,det2,deti(2)
!       complex(16)                :: q
!       integer                    :: INFO,INFO2
!       integer                    :: nnn
!       logical,optional           :: check_nan,block_matrix
!       integer,optional           :: c
! 
!       if(mat_3_3_qc(mat,detb,det2b,pdetb)) return
! 
!       if(present(det2b))then; det2b=0.;detb=0.;pdetb=0; endif
! 
!       if(present(block_matrix))then
!        if(block_matrix)then
!          nnn=n/2
!          if(mod(n,2)>0) stop 'error invmat_comp block_matrix, but linear size is odd'
!          if(present(detb))then
!            call invit(1,nnn,det,det2,pdet)
!            call invit(nnn+1,n,ddet,ddet2,ppdet)
!            mat(1:nnn,nnn+1:n)=0.
!            mat(nnn+1:n,1:nnn)=0.
!            pdetb = pdet + ppdet
!            detb  = det  * ddet
!            det2b = det2 * ddet2
!          else
!            call invit(1,nnn)
!            call invit(nnn+1,n)
!            mat(1:nnn,nnn+1:n)=0.
!            mat(nnn+1:n,1:nnn)=0.
!          endif
!          goto 35
!        endif
!       endif
! 
!       call invit(1,n,detb,det2b,pdetb)
! 35    continue
! 
!       return
! 
!       contains
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!    subroutine zgeco_zgedi(i1,i2)
!    implicit none
!    integer :: i1,i2,nnn,i,j
!      nnn=i2-i1+1
!      call geco__(nnn,mat(i1:i2,i1:i2),piv(i1:i2))
!      call gedi__(nnn,mat(i1:i2,i1:i2),piv(i1:i2),deti)
!    end subroutine
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!     subroutine invit(i1,i2,det,det2,pdet)
!     implicit none
!     integer :: i1,i2
!     real(16),optional          :: det
!     integer,optional           :: pdet
!     complex(16),optional       :: det2
!        if(mat_3_3_qc(mat(i1:i2,i1:i2),det,det2,pdet)) return
!        call zgeco_zgedi(i1,i2)
!        if(present(det)) call get_det_from_zgeco(deti,det,det2,pdet)
!     end subroutine
! 
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!       end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       subroutine invmat_comp2(n,mat,detb2,detb,pdetb,diagmat)
!       implicit none
!       integer                   :: pdet,n,i,j
!       integer,optional          :: pdetb
!       real(8)                   :: big,det
!       real(8),optional          :: detb
!       complex(8),intent(inout)  :: mat(:,:)
!       complex(8)                :: det2,q
!       complex(8),optional       :: detb2
!       integer(4)                :: INFO,INFO2
!       integer                   :: piv(size(mat,1))
!       logical,optional          :: diagmat
! 
!       if(size(mat,1)/=size(mat,2)) stop 'error invmat_comp2 try invert rectangular matrix'
! 
!       if(present(diagmat))then
!        if(diagmat)then
!         do i=1,size(mat,1)
!          if(abs(mat(i,i))>1.d-13)then
!           mat(i,i)=1.d0/mat(i,i)
!          else
!           mat(i,i)=1.d14
!          endif
!         enddo
!         return
!        endif
!       endif
! 
!       if(mat_3_3_c(mat,detb,detb2,pdetb)) return
!       det=1.d0; det2=1.d0; pdet=0
!       call GETRF__(mat,piv)
! 
!       if(present(detb2))then
!        do i=1,n
!         det2 = det2 * mat(i,i)
!        enddo
!        det=abs(det2)
!        if(det<error) then
!         det2=1.
!        else
!         det2=det2/det
!        endif
!        call rescale(det,pdet)
!       endif
! 
!       call GETRI__(mat,piv)
! 
!       if(present(detb2))then; detb2=det2; detb=det; pdetb=pdet; endif
! 
!       return
!       end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       subroutine invmat_reals(n,mat,det2b,detb,pdetb,check_nan,c,block_matrix,diagmat)
!       implicit none
!       real(4),optional      :: detb
!       integer,optional      :: pdetb
!       complex(4),optional   :: det2b
!       real(4)               :: det,ddet,q
!       complex(4)            :: det2
!       integer               :: pdet,ppdet,nnn
!       real(4)               :: big
!       real(4),intent(inout) :: mat(:,:)
!       real(4)               :: rtemp
!       logical,optional      :: check_nan,diagmat
!       integer,optional      :: c
!       logical,optional      :: block_matrix
!       integer               :: i,j,k,n
! 
!       if(size(mat,1)/=size(mat,2)) stop 'error invmat_real rectangular matrix'
! 
!       if(present(diagmat))then
!        if(diagmat)then
!         do i=1,size(mat,1)
!          if(abs(mat(i,i))>1.d-15)then
!           mat(i,i)=1.d0/mat(i,i)
!          else
!           mat(i,i)=1.d15
!          endif
!         enddo
!         return
!        endif
!       endif
! 
!       if(mat_3_3_rs(mat,detb,det2b,pdetb)) return
! 
!       if(present(check_nan)) call erase_divergence(mat)
!       if(size(mat(:,1))/=n.or.size(mat(1,:))/=n) stop 'error invmat_real matrix not square'
!  
!       if(present(block_matrix))then
!        if(block_matrix)then
!         nnn=n/2
!         if(mod(n,2)>0) stop 'error invmat_reals block_matrix, but linear size is odd'
!         mat(1:nnn,nnn+1:n)=0.d0 ; mat(nnn+1:n,1:nnn)=0.d0
!         if(present(detb))then
!          call invit(mat,1,nnn,det,det2,pdet)
!          call invit(mat,nnn+1,n,ddet,det2,ppdet)
!          pdetb = pdet + ppdet
!          detb  = det  * ddet
!          det2b = 1.
!         else
!          call invit(mat,1,nnn)
!          call invit(mat,nnn+1,n)
!         endif
!         goto 35
!        endif
!       endif
!       call invit(mat,1,n,detb,det2b,pdetb)
! 35    continue
!       call erase_divergence(mat)
! 
!       return
!       
!       
!       contains
!       
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!       subroutine invit(mat_in,i1,i2,det,det2,pdet)
!       implicit none
!       integer                   :: i1,i2,n
!       integer                   :: piv(i2-i1+1)
!       real(4),target            :: mat_in(:,:)
!       real(4),pointer           :: mat(:,:)
!       real(4)                   :: stock(i2-i1+1)
!       real(4)                   :: det_
!       integer                   :: pdet_
!       real(4),optional          :: det
!       integer,optional          :: pdet
!       complex(4),optional       :: det2
!   
!       n=i2-i1+1
!       mat=>mat_in(i1:i2,i1:i2) 
!       if(mat_3_3_rs(mat,det,det2,pdet)) return
! 
!       if(.not.present(det))then
!        if(n>n_openmp_rout.and.use_openmp_invmat)then
!          call invert_openmp(mat)
!          return
!        endif
!        if(flag_use_invmat_jordan_real)then
!          call invmat_jordan(n,mat)
!          return
!        endif
!       endif
!  
!       det_ = 1.d0 ; pdet_ = 0
!       do i=1,n
!         big = abs(mat(i,i))
!         piv(i) = i
!         do j=i,n
!         if(abs(mat(i,j))>big) then
!           big = abs(mat(i,j))
!           piv(i) = j
!         endif
!         enddo
!         big = mat(i,piv(i))
!         !-------------TEST---------------------------------------!
!         if(abs(big)<rerror) then
!           write(*,*) 'invmat, small term : danger in inverting matrix'
!           write(*,*) 'term is : ', big,mat(i,piv(i))
!           write(*,*) 'i,j     : ', i,piv(i)
!           big=smallest_pivot
!           if(strongstop) stop 'error in invmat real'
!           write(*,*) 'anyway, goes on...'
!         endif
!         !-------------TEST---------------------------------------!
!        det_ = abs(det_*big)
!        call rescale(det_,pdet_)
!         mat(i,piv(i)) = 0.d0
!         rtemp = mat(i,i)
!         do j=1,n
!           stock(j) = mat(j,piv(i))/big
!         enddo
!         if(piv(i)/=i) then
!          do j=1,n
!            mat(j,piv(i)) = mat(j,i)
!          enddo
!         endif
!         mat(i,piv(i)) = rtemp
!         do k=1,n
!           do j=1,n
!             mat(j,k) = mat(j,k) - stock(j)*mat(i,k)
!           enddo
!         enddo
!         do j=1,n
!            mat(j,i) = stock(j)
!         enddo
!         do j=1,n
!           mat(i,j) = -mat(i,j)/big
!         enddo
!         mat(i,i)=1/big
!        enddo
!        do i=n,1,-1
!         do j=1,n
!          stock(j) = mat(i,j)
!          mat(i,j) = mat(piv(i),j)
!          mat(piv(i),j) = stock(j)
!         enddo
!        enddo     
!       if(present(det))then; det=det_; det2=1.d0; pdet=pdet_ ; endif
!       end subroutine
!       
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!         
!       
!       end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       subroutine invmat_real(n,mat,det2b,detb,pdetb,check_nan,c,block_matrix,diagmat)
!       implicit none
!       real(8),optional      :: detb
!       integer,optional      :: pdetb
!       complex(8),optional   :: det2b
!       real(8)               :: det,ddet,q
!       complex(8)            :: det2
!       integer               :: pdet,ppdet,nnn
!       real(8)               :: big
!       real(8),intent(inout) :: mat(:,:)
!       real(8)               :: rtemp
!       logical,optional      :: check_nan,diagmat
!       integer,optional      :: c
!       logical,optional      :: block_matrix
!       integer               :: i,j,k,n
! 
!       if(size(mat,1)/=size(mat,2)) stop 'error invmat_real rectangular matrix'
! 
!       if(present(diagmat))then
!        if(diagmat)then
!         do i=1,size(mat,1)
!          if(abs(mat(i,i))>1.d-15)then
!           mat(i,i)=1.d0/mat(i,i)
!          else
!           mat(i,i)=1.d15
!          endif
!         enddo
!         return
!        endif
!       endif
! 
!       if(mat_3_3_r(mat,detb,det2b,pdetb)) return
! 
!       if(present(check_nan)) call erase_divergence(mat)
!       if(size(mat(:,1))/=n.or.size(mat(1,:))/=n) stop 'error invmat_real matrix not square'
!  
!       if(present(c))then      
!        call randomize_matrix(mat,amp=invmat_error,flag=.true.,kk2=c)
!       else
!        call randomize_matrix(mat,amp=invmat_error)
!       endif
! 
!       if(present(block_matrix))then
!        if(block_matrix)then
!         nnn=n/2
!         if(mod(n,2)>0) stop 'error invmat_real block_matrix, but linear size is odd'
!         mat(1:nnn,nnn+1:n)=0.d0 ; mat(nnn+1:n,1:nnn)=0.d0
!         if(present(detb))then
!          call invit(mat,1,nnn,det,det2,pdet)
!          call invit(mat,nnn+1,n,ddet,det2,ppdet)
!          pdetb = pdet + ppdet
!          detb  = det  * ddet 
!          det2b = 1.
!         else
!          call invit(mat,1,nnn)
!          call invit(mat,nnn+1,n)
!         endif
!         goto 35
!        endif
!       endif
!       call invit(mat,1,n,detb,det2b,pdetb)
! 35    continue
! 
!       call erase_divergence(mat)
! 
!       return
!       
!       
!       contains
!       
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!     subroutine force_sing_prec(mat,det,det2,pdet)
!      real(8)                   :: mat(:,:)
!      real(4)                   :: mat_(size(mat,1),size(mat,2))
!      real(8),optional          :: det
!      integer,optional          :: pdet
!      complex(8),optional       :: det2
!      real(4)                   :: det_
!      complex(4)                :: det2_
!      integer                   :: pdet_
!      if(present(det))then
!       mat_=mat
!       call invmat_reals(size(mat,1),mat_,det2_,det_,pdet_)
!       mat=mat_
!       pdet=pdet_
!       det2=det2_
!       det=det_
!      else
!       mat_=mat
!       call invmat_reals(size(mat_,1),mat_)
!       mat=mat_
!      endif
!     end subroutine
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!       
!       subroutine invit(mat_in,i1,i2,det,det2,pdet)
!       implicit none
!       integer                   :: i1,i2,n
!       integer                   :: piv(i2-i1+1)
!       real(8),target            :: mat_in(:,:)
!       real(8),pointer           :: mat(:,:)
!       real(8)                   :: stock(i2-i1+1)
!       real(8)                   :: det_
!       integer                   :: pdet_
!       real(8),optional          :: det
!       integer,optional          :: pdet
!       complex(8),optional       :: det2
! 
!       mat=>mat_in(i1:i2,i1:i2) 
!       n=i2-i1+1
!       if(mat_3_3_r(mat,det,det2,pdet)) return
! 
!       if(force_invmat_single_prec)then
!         call force_sing_prec(mat,det,det2,pdet)
!         return
!       endif
!       if(.not.present(det))then
!        if(n>n_openmp_rout.and.use_openmp_invmat)then
!          call invert_openmp(mat)
!          return
!        endif
!        ! ebl: Removing GPU functionality
!        ! if(use_cuda_routines.and.n>n_cuda_rout)then
!        !   call diago_cuda_it_r(n,mat)
!        !   return
!        ! endif
!        if(flag_use_invmat_jordan_real)then
!          call invmat_jordan(n,mat)
!          return
!        endif
!        if(diag_use_LU_instead_of_pivot)then
!          call D_INVMG(size(mat,1),size(mat,2),mat)
!          return
!        endif
!       endif
!  
!       det_ = 1.d0 ; pdet_ = 0
!       do i=1,n
!         big = abs(mat(i,i))
!         piv(i) = i
!         do j=i,n
!         if(abs(mat(i,j))>big) then
!           big = abs(mat(i,j))
!           piv(i) = j
!         endif
!         enddo
!         big = mat(i,piv(i))
!         !-------------TEST---------------------------------------!
!         if(abs(big)<rerror) then
!           write(*,*) 'invmat, small term : danger in inverting matrix'
!           write(*,*) 'term is : ', big,mat(i,piv(i))
!           write(*,*) 'i,j     : ', i,piv(i)
!           big=smallest_pivot
!           if(strongstop) stop 'error in invmat real'
!           write(*,*) 'anyway, goes on...'
!         endif
!         !-------------TEST---------------------------------------!
!         det_ = abs(det_*big)
!         call rescale(det_,pdet_)
!         mat(i,piv(i)) = 0.d0
!         rtemp = mat(i,i)
! 
!         do j=1,n
!           stock(j) = mat(j,piv(i))/big
!         enddo
!         if(piv(i)/=i) then
!          do j=1,n
!            mat(j,piv(i)) = mat(j,i)
!          enddo
!         endif
!         mat(i,piv(i)) = rtemp
!         do k=1,n
!           do j=1,n
!             mat(j,k) = mat(j,k) - stock(j)*mat(i,k)
!           enddo
!         enddo
!         do j=1,n
!            mat(j,i) = stock(j)
!         enddo
!         do j=1,n
!           mat(i,j) = -mat(i,j)/big
!         enddo
!         mat(i,i)=1/big
!        enddo
!        do i=n,1,-1
!         do j=1,n
!          stock(j) = mat(i,j)
!          mat(i,j) = mat(piv(i),j)
!          mat(piv(i),j) = stock(j)
!         enddo
!        enddo     
!       if(present(det))then; det=det_; det2=1.d0; pdet=pdet_ ; endif
!       end subroutine
!       
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!         
!       
!       end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       subroutine invmat_real_quad(n,mat,det2b,detb,pdetb,check_nan,c,block_matrix,diagmat)
!       implicit none
!       real(16),intent(inout):: mat(:,:)
!       real(16),optional     :: detb
!       integer,optional      :: pdetb
!       complex(16),optional  :: det2b
!       real(16)              :: det,ddet,q
!       complex(16)           :: det2
!       integer               :: pdet,ppdet
!       real(16)              :: big
!       real(16)              :: rtemp
!       logical,optional      :: check_nan,block_matrix,diagmat
!       integer,optional      :: c
!       integer               :: i,j,k,n,nnn
! 
!       if(size(mat,1)/=size(mat,2)) stop 'error invmat quad rectangular matrix'
! 
!       if(present(diagmat))then
!        if(diagmat)then
!         do i=1,size(mat,1)
!          if(abs(mat(i,i))>1.d-25)then
!           mat(i,i)=1.d0/mat(i,i)
!          else
!           mat(i,i)=1.d25
!          endif
!         enddo
!         return
!        endif
!       endif
! 
!       if(mat_3_3_q(mat,detb,det2b,pdetb)) return
! 
!       if(present(block_matrix))then
!        if(block_matrix)then
!         nnn=n/2
!         if(mod(n,2)>0) stop 'error invmat_real_quad block_matrix, but linear size is odd'
!         mat(1:nnn,nnn+1:n)=0.d0 ; mat(nnn+1:n,1:nnn)=0.d0
!         if(present(detb))then
!          call invit(mat,1,nnn,det,det2,pdet)
!          call invit(mat,nnn+1,n,ddet,det2,ppdet)
!          pdetb = pdet + ppdet
!          detb  = det  * ddet
!          det2b = 1.
!         else
!          call invit(mat,1,nnn)
!          call invit(mat,nnn+1,n)
!         endif
!         goto 35
!        endif
!       endif
!       call invit(mat,1,n,detb,det2b,pdetb)
! 35    continue
! 
! 
! 
!    contains
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!       subroutine invit(mat_in,i1,i2,det,det2,pdet)
!       implicit none
!       integer                    :: i1,i2,n
!       integer                    :: piv(i2-i1+1)
!       real(16),target            :: mat_in(:,:)
!       real(16),pointer           :: mat(:,:)
!       real(16)                   :: stock(i2-i1+1)
!       real(16)                   :: det_
!       integer                    :: pdet_
!       real(16),optional          :: det
!       integer,optional           :: pdet
!       complex(16),optional       :: det2
! 
!       mat=>mat_in(i1:i2,i1:i2)
!       n=i2-i1+1
!       if(mat_3_3_q(mat,det,det2,pdet)) return
! 
!       if(.not.present(det))then
!        if(diag_use_LU_instead_of_pivot)then
!         call Q_INVMG(size(mat,1),size(mat,2),mat)
!         return
!        endif
!       endif
! 
!       det_ = 1.d0; pdet_ = 0
!       do i=1,n
!         big = abs(mat(i,i))
!         piv(i) = i
!         do j=i,n
!         if(abs(mat(i,j))>big) then
!           big = abs(mat(i,j))
!           piv(i) = j
!         endif
!         enddo
!         big = mat(i,piv(i))
! 
!         det_ = abs(det_*big)
!         call rescale(det_,pdet_)
! 
!         mat(i,piv(i)) = 0.d0
!         rtemp = mat(i,i)
! 
!         do j=1,n
!          stock(j) = mat(j,piv(i))/big
!         enddo
!         if(piv(i)/=i) then
!          do j=1,n
!            mat(j,piv(i)) = mat(j,i)
!          enddo
!         endif
!         mat(i,piv(i)) = rtemp
!         do k=1,n
!           do j=1,n
!             mat(j,k) = mat(j,k) - stock(j)*mat(i,k)
!           enddo
!         enddo
!         do j=1,n
!            mat(j,i) = stock(j)
!         enddo
!         do j=1,n
!           mat(i,j) = -mat(i,j)/big
!         enddo
!         mat(i,i)=1/big
!        enddo
! 
!        do i=n,1,-1
!         do j=1,n
!          stock(j) = mat(i,j)
!          mat(i,j) = mat(piv(i),j)
!          mat(piv(i),j) = stock(j)
!         enddo
!        enddo
! 
!        if(present(det))then; det=det_; det2=1.d0; pdet=pdet_ ; endif
! 
!       end subroutine
! 
!     !------------------!
!     !------------------!
!     !------------------!
!     !------------------!
! 
!       end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! !regular one : no control
! 
!     subroutine zgedi_(a,lda,n,ipvt,det,work,job)
!      !subroutines of LINPACK
!       implicit none
!       integer    :: lda,n,ipvt(1),job
!       complex(8) :: a(lda,1),det(2),work(1)
!       complex(8) :: t
!       real(8)     :: ten
!       integer    :: i,j,k,kb,kp1,l,nm1
!       complex(8) :: zdum
!       complex(8) :: zdumr,zdumi
! 
!       if (job/10 == 0) go to 70
!          det(1) = (1.0d0,0.0d0)
!          det(2) = (0.0d0,0.0d0)
!          ten = 10.0d0
!          do 50 i = 1, n
!             if (ipvt(i) /= i) det(1) = -det(1)
!             det(1) = a(i,i)*det(1)
!             if (cabs1(det(1)) == 0.0d0) go to 60
!    10       if (cabs1(det(1)) >= 1.0d0) go to 20
!                det(1) = cmplx(ten,0.0d0,kind=8)*det(1)
!                det(2) = det(2) - (1.0d0,0.0d0)
!             go to 10
!    20       continue
!    30       if (cabs1(det(1))<ten) go to 40
!                det(1) = det(1)/cmplx(ten,0.0d0,kind=8)
!                det(2) = det(2) + (1.0d0,0.0d0)
!             go to 30
!    40       continue
!    50    continue
!    60    continue
!    70    continue
!          if (mod(job,10) == 0) go to 150
!          do k = 1, n
!             a(k,k) = (1.0d0,0.0d0)/a(k,k)
!             t = -a(k,k)
!             call zscal(k-1,t,a(1,k),1)
!             kp1 = k + 1
!             if (n .lt. kp1) go to 90
!             do j = kp1, n
!                t = a(k,j)
!                a(k,j) = (0.0d0,0.0d0)
!                call zaxpy(k,t,a(1,k),1,a(1,j),1)
!             enddo 
!    90       continue
!          enddo 
!          nm1 = n - 1
!          if (nm1 .lt. 1) go to 140
!          do 130 kb = 1, nm1
!             k = n - kb
!             kp1 = k + 1
!             do 110 i = kp1, n
!                work(i) = a(i,k)
!                a(i,k) = (0.0d0,0.0d0)
!   110       continue
!             do j = kp1, n
!                t = work(j)
!                call zaxpy(n,t,a(1,j),1,a(1,k),1)
!             enddo 
!             l = ipvt(k)
!             if (l/=k) call zswap(n,a(1,k),1,a(1,l),1)
!   130    continue
!   140    continue
!   150 continue
! 
!       return
!       end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! !optimized one : overflow control
! 
!       subroutine zgedi(a,lda,n,ipvt,det,work,job)
!       !subroutines of LINPACK
!       implicit none
!       integer    :: lda,n,ipvt(1),job
!       complex(8) :: a(lda,1),det(2),work(1)
!       complex(8) :: t
!       real(8)    :: ten
!       integer    :: i,j,k,kb,kp1,l,nm1
! 
!       if (job/10 == 0) go to 70
!          det(1) = (1.0d0,0.0d0)
!          det(2) = (0.0d0,0.0d0)
!          ten = 10.0d0
!          do 50 i = 1, n
!             if (ipvt(i) /= i) det(1) = -det(1)
!             det(1) = a(i,i)*det(1)
!             if (cabs1(det(1)) == 0.0d0) go to 60
!    10       if (cabs1(det(1)) >= 1.0d0) go to 20
!                det(1) = cmplx(ten,0.0d0,kind=8)*det(1)
!                det(2) = det(2) - (1.0d0,0.0d0)
!             go to 10
!    20       continue
!    30       if (cabs1(det(1))<ten) go to 40
!                det(1) = det(1)/cmplx(ten,0.0d0,kind=8)
!                det(2) = det(2) + (1.0d0,0.0d0)
!             go to 30
!    40       continue
!    50    continue
!    60    continue
!    70    continue
!          if (mod(job,10) == 0) go to 150
!          do k = 1, n
!             a(k,k) = (1.0d0,0.0d0)/a(k,k)
!             t = -a(k,k)
!             call zscal(k-1,t,a(1:k-1,k),1)
!             kp1 = k + 1
!             if (n .lt. kp1) go to 90
!             do j = kp1, n
!                t = a(k,j)
!                a(k,j) = (0.0d0,0.0d0)
!                call zaxpyb(k,t,a(:,k),1,a(:,j),1)
!             enddo 
!    90       continue
!          enddo 
!          nm1 = n - 1
!          if (nm1 .lt. 1) go to 140
!          do 130 kb = 1, nm1
!             k = n - kb
!             kp1 = k + 1
!             do 110 i = kp1, n
!                work(i) = a(i,k)
!                a(i,k) = (0.0d0,0.0d0)
!   110       continue
!             do j = kp1, n
!                t = work(j)
!                call zaxpyb(n,t,a(:,j),1,a(:,k),1)
!             enddo 
!             l = ipvt(k)
!             if (l/=k) call swap(n,a(:,k),1,a(:,l),1)
!   130    continue
!   140    continue
!   150 continue
! 
!       return
!       end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! !regular one : no control
! 
!      subroutine zgeco_(a,lda,n,ipvt,rcond,z)
!       implicit none
!       integer    :: lda,n,ipvt(1)
!       complex(8) :: a(lda,1),z(1)
!       real(8)    :: rcond
!       complex(8) :: zdotc,ek,t,wk,wkm
!       real(8)    :: anorm,s,dzasum,sm,ynorm,ss
!       integer    :: info,j,k,kb,kp1,l
!       complex(8) :: zdum,zdum1,zdum2
!       complex(8) :: zdumr,zdumi
! 
!       anorm = 0.0d0
!       do j=1,n
!          anorm=dmax1(anorm,dzasum(n,a(1,j),1))
!       enddo
!       call zgefa(a,lda,n,ipvt,info)
!       ek =CMPLX(1.0d0,0.0d0,kind=8)
!       z  =CMPLX(0.0d0,0.0d0,kind=8)
! 
!       do k=1,n
!          if(cabs1(z(k)).ne.0.0d0) ek=csign1(ek,-z(k))
!          if(cabs1(ek-z(k)).le.cabs1(a(k,k))) go to 30
!             ss=cabs1(ek-z(k))
!             if(abs(ss)<epsilonr) ss=epsilonr
!             s = cabs1(a(k,k))/ss
!             call zdscal(n,s,z,1)
!             ek = s*ek
!    30    continue
!          wk  =  ek - z(k)
!          wkm = -ek - z(k)
!          s   =  cabs1(wk)
!          sm  =  cabs1(wkm)
!          if(cabs1(a(k,k))<rerror) go to 40
!             wk = wk/conjg(a(k,k))
!             wkm = wkm/conjg(a(k,k))
!          go to 50
!    40    continue
!             wk  = (1.0d0,0.0d0)
!             wkm = (1.0d0,0.0d0)
!    50    continue
!          kp1 = k + 1
!          if (kp1 > n) go to 90
!             do 60 j = kp1, n
!                sm   = sm   + cabs1(z(j)+wkm*conjg(a(k,j)))
!                z(j) = z(j) + wk*conjg(a(k,j))
!                s    = s    + cabs1(z(j))
!    60       continue
!             if (s>=sm) go to 80
!                t  = wkm - wk
!                wk = wkm
!                do 70 j = kp1, n
!                   z(j) = z(j) + t*conjg(a(k,j))
!    70          continue
!    80       continue
!    90    continue
!          z(k) = wk
!       enddo
! 
!       if(abs(dzasum(n,z,1))>rerror) then
!        s = 1.d0/dzasum(n,z,1)
!       else
!        s = 1.d0/rerror
!       endif
!       call zdscal(n,s,z,1)
!       do 120 kb = 1, n
!          k = n + 1 - kb
!          if (k .lt. n) z(k) = z(k) + zdotc(n-k,a(k+1,k),1,z(k+1),1)
!          if (cabs1(z(k)) .le. 1.0d0) go to 110
!             s = 1.0d0/cabs1(z(k))
!             call zdscal(n,s,z,1)
!   110    continue
!          l = ipvt(k)
!          t = z(l)
!          z(l) = z(k)
!          z(k) = t
!   120 continue
! 
!       if(abs(dzasum(n,z,1))>rerror) then
!        s = 1.d0/dzasum(n,z,1)
!       else
!        s = 1.d0/rerror
!       endif
! 
!       call zdscal(n,s,z,1)
!       ynorm = 1.0d0
!       do 140 k = 1, n
!          l = ipvt(k)
!          t = z(l)
!          z(l) = z(k)
!          z(k) = t
!          if (k .lt. n) call zaxpy(n-k,t,a(k+1,k),1,z(k+1),1)
!          if (cabs1(z(k)) .le. 1.0d0) go to 130
!             s = 1.0d0/cabs1(z(k))
!             call zdscal(n,s,z,1)
!             ynorm = s*ynorm
!   130    continue
!   140 continue
! 
!       if(abs(dzasum(n,z,1))>rerror) then
!        s = 1.0d0/dzasum(n,z,1)
!       else
!        s = 1.d0/rerror
!       endif
! 
!       call zdscal(n,s,z,1)
!       ynorm = s*ynorm
!       do 160 kb = 1, n
!          k = n + 1 - kb
!          if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 150
!             s = cabs1(a(k,k))/cabs1(z(k))
!             call zdscal(n,s,z,1)
!             ynorm = s*ynorm
!   150    continue
!          if (cabs1(a(k,k)) > rerror) then
!           z(k) = z(k)/a(k,k)
!          else
!           z(k) = (1.0d0,0.0d0)
!          endif
!          t = -z(k)
!          call zaxpy(k-1,t,a(1,k),1,z(1),1)
!   160 continue
! 
!       if(abs(dzasum(n,z,1))>rerror) then
!        s = 1.d0/dzasum(n,z,1)
!       else
!        s = 1.d0/rerror
!       endif
! 
!       call zdscal(n,s,z,1)
!       ynorm = s*ynorm
!       if(abs(anorm) > rerror) then
!        rcond = ynorm/anorm
!       else
!        rcond = 1.d0/rerror
!       endif
!       
!       return
!       end subroutine
! 
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! !optimized one : control on the floating overflow
! 
!       subroutine zgeco(a,lda,n,ipvt,rcond,z)
!       implicit none
!       integer    :: lda,n,ipvt(1)
!       complex(8) :: a(lda,1),z(1)
!       real(8)    :: rcond
!       complex(8) :: ek,t,wk,wkm
!       real(8)    :: anorm,ss
!       real(16)   :: s,sm,ynorm,tt
!       real(8)    :: ssum
!       integer    :: info,j,k,kb,kp1,l
! 
!       anorm = 0.0d0
!       do j=1,n
!          anorm=max(anorm,zsum(n,a(:,j),1))
!       enddo
!       call zgefa(a,lda,n,ipvt,info)
!       ek =CMPLX(1.0d0,0.0d0,kind=8)
!       z  =CMPLX(0.0d0,0.0d0,kind=8)
! 
!       do k=1,n
!          if(cabs1(z(k)).ne.0.0d0) ek=csign1(ek,-z(k))
!          if(cabs1(ek-z(k)).le.cabs1(a(k,k))) go to 30
!             ss=cabs1(ek-z(k))
!             if(abs(ss)<epsilonr) ss=epsilonr
!             s = cabs1(a(k,k))/ss
!             if(abs(s)>MAX_REAL) s=qsign(s)*MAX_REAL
!             z=z*dble(s)
!             ek = s*ek
!    30    continue
!          wk  =  ek - z(k)
!          wkm = -ek - z(k)
!          s   =  cabs1(wk)
!          sm  =  cabs1(wkm)
!          if(cabs1(a(k,k))<rrerror) go to 40
!             tt = wk
!             tt = tt/conjg(a(k,k)) 
!             if(abs(tt)>MAX_REAL) goto 40
!             wk = tt
!             tt = wkm  
!             tt = tt / conjg(a(k,k))
!             if(abs(tt)>MAX_REAL) goto 40
!             wkm = tt
!          go to 50
!    40    continue
!             wk  = (1.0d0,0.0d0)
!             wkm = (1.0d0,0.0d0)
!    50    continue
!          kp1 = k + 1
!          if (kp1 > n) go to 90
!             do 60 j = kp1, n
!                tt   =  cabs1(z(j)+wkm*conjg(a(k,j))) 
!                if(abs(sm)>MAX_REAL.or.abs(tt)>MAX_REAL)then
!                 sm=MAX_REAL
!                 s=MAX_REAL-1
!                 goto 80
!                endif
!                sm   =   sm   + tt
!                z(j) = z(j) + wk*conjg(a(k,j))
!                s    = s    + cabs1(z(j))
!    60       continue
!             if (s>=sm) go to 80
!                t  = wkm - wk
!                wk = wkm
!                do 70 j = kp1, n
!                   z(j) = z(j) + t*conjg(a(k,j))
!    70          continue
!    80       continue
!    90    continue
!          z(k) = wk
!       enddo
! 
!       ssum=zsum(n,z,1)
!       if(abs(ssum)>rrerror) then
!        s = 1.d0/ssum
!       else
!        s = 1.d0/rrerror
!       endif
! 
!       if(abs(s)>MAX_REAL) s=qsign(s)*MAX_REAL
!       z=z*dble(s)
! 
!       do 120 kb = 1, n
!          k = n + 1 - kb
!          if (k .lt. n) z(k) = z(k) + zdot(n-k,a(k+1:n,k),1,z(k+1:n),1)
!          if (cabs1(z(k)) .le. 1.0d0) go to 110
!           if(abs(cabs1(z(k)))>rrerror)then
!             s = 1.0d0/cabs1(z(k))
!           else
!             s = 1.d0/rrerror
!           endif
!           if(abs(s)>MAX_REAL) s=qsign(s)*MAX_REAL
!           z=z*dble(s)
!   110    continue
!          l = ipvt(k)
!          t = z(l)
!          z(l) = z(k)
!          z(k) = t
!   120 continue
! 
!       ssum=zsum(n,z,1)
!       if(abs(ssum)>rrerror) then
!        s = 1.d0/ssum
!       else
!        s = 1.d0/rrerror
!       endif
! 
!       if(abs(s)>MAX_REAL) s=qsign(s)*MAX_REAL
!       z=z*dble(s)
! 
!       ynorm = 1.0d0
!       do 140 k = 1, n
!          l = ipvt(k)
!          t = z(l)
!          z(l) = z(k)
!          z(k) = t
!          if (k .lt. n) call zaxpyb(n-k,t,a(k+1:n,k),1,z(k+1:n),1)
!          if (cabs1(z(k)) .le. 1.0d0) go to 130
!             if(abs(cabs1(z(k)))>rrerror)then
!              s = 1.0d0/cabs1(z(k))
!             else
!              s = 1.d0/rrerror
!             endif
!             if(abs(s)>MAX_REAL) s=qsign(s)*MAX_REAL
!             z=z*dble(s)
!             ynorm = s*ynorm
!   130    continue
!   140 continue
! 
!       ssum=zsum(n,z,1)
!       if(abs(ssum)>rrerror) then
!        s = 1.0d0/ssum
!       else
!        s = 1.d0/rrerror
!       endif
! 
!       if(abs(s)>MAX_REAL) s=qsign(s)*MAX_REAL
!       z=z*dble(s)
! 
!       ynorm = s*ynorm
!       do 160 kb = 1, n
!          k = n + 1 - kb
!          if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 150
!             s = cabs1(a(k,k))/cabs1(z(k))
!             z=z*dble(s)
!             ynorm = s*ynorm
!   150    continue
!          if(abs(a(k,k))>rerror .and. abs(z(k))<MAX_REAL) then
!           z(k) = z(k)/a(k,k)
!          else
!           z(k) = (1.0d0,0.0d0)
!          endif
!          t = -z(k)
!          call zaxpyb(k-1,t,a(1:k-1,k),1,z(1:k-1),1)
!   160 continue
! 
!       ssum=zsum(n,z,1)
!       if(abs(ssum)>rrerror) then
!        s = 1.d0/ssum
!       else
!        s = 1.d0/rrerror
!       endif
! 
!       if(abs(s)>MAX_REAL) s=qsign(s)*MAX_REAL
!       z=z*dble(s)
! 
!       ynorm = s*ynorm
! 
!       if(abs(anorm) > rrerror) then
!        tt=ynorm/anorm
!        if(abs(tt)<MAX_REAL)then
!         rcond = tt
!        else
!         rcond = qsign(ynorm)*sign(1.d0,anorm)*MAX_REAL
!        endif
!       else
!        rcond = MAX_REAL
!       endif
!       
!       return
!       end subroutine

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
! 
!       subroutine zgefa(a,lda,n,ipvt,info)
!       implicit none
!       integer    :: lda,n,ipvt(1),info
!       complex(8) :: a(lda,1)
!       complex(8) :: t
!       integer    :: izamax,j,k,kp1,l,nm1
! 
!       info = 0
!       nm1 = n - 1
!       if (nm1 < 1) go to 70
!       do 60 k = 1, nm1
!          kp1 = k + 1
!          l = izamax(n-k+1,a(k,k),1) + k - 1
!          ipvt(k) = l
!          if (cabs1(a(l,k)) == 0.0d0) go to 40
!             if (l == k) go to 10
!                t = a(l,k)
!                a(l,k) = a(k,k)
!                a(k,k) = t
!    10       continue
!             t = -(1.0d0,0.0d0)/a(k,k)
!             call zscal(n-k,t,a(k+1,k),1)
!             do 30 j = kp1, n
!                t = a(l,j)
!                if (l == k) go to 20
!                   a(l,j) = a(k,j)
!                   a(k,j) = t
!    20          continue
!                call zaxpyb(n-k,t,a(k+1:n,k),1,a(k+1:n,j),1)
!    30       continue
!          go to 50
!    40    continue
!             info = k
!    50    continue
!    60 continue
!    70 continue
!       ipvt(n) = n
!       if (cabs1(a(n,n)) == 0.0d0) info = n
!       return
!       end subroutine

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

subroutine QR_decomp(L,R,matin,Qmat,Rmat)
implicit none
integer               :: L,R,INFO,i,j
complex(8),intent(in) :: matin(L,R)
complex(8),intent(out):: Qmat(L,R),Rmat(R,R)
complex(8)            :: mat2(L,R),TAU(R)
complex(8)            :: WORK(2*R)

 mat2=matin
 call ZGEQRF(L,R,mat2,L, TAU, WORK, 2*R,INFO)
 Rmat=0.
 do i=1,R
  do j=i,R
   Rmat(i,j)=mat2(i,j)
  enddo
 enddo
 call ZUNGQR(L,R,R,mat2,L,TAU,WORK,2*R,INFO)
 Qmat=mat2

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

! Fonction        : Inversion d'une matrice carree
! Notes           : A -> A^(-1)=AI
!                   Provient de Numerical Recipes, p.40.
! Dependances     : LUDCMP,LUBKSB (externes)
! 
!      SUBROUTINE Q_INVMG(NP,NM,A)
!      IMPLICIT NONE
!      REAL(16)                                   :: D     
!      INTEGER                                    :: I,J    
!      INTEGER,INTENT(IN)                         :: NP,NM
!      INTEGER,DIMENSION (NP)                     :: INDX  
!      REAL(16),DIMENSION (NP,NP)                 :: Y   
!      REAL(16),DIMENSION (NP,NP),INTENT(INOUT)   :: A
! 
!      DO I=1,NM
!         DO J=1,NM
!            Y(I,J)=0.D0
!         END DO
!         Y(I,I)=1.D0
!      END DO
!      CALL Q_LUDCMP(A,NM,NP,INDX,D)
!      DO J=1,NM
!         CALL Q_LUBKSB(A,NM,NP,INDX,Y(1,J))
!      END DO
!      A=Y
! 
!      END SUBROUTINE
! 
! 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

! Nom            : LUDCMP 
! Fonction       : Decomposition gaussienne A = LU
! Notes          : Provient entierement de Numerical Recipes
! 
!    SUBROUTINE Q_LUDCMP(A,N,NP,INDX,D)
!    IMPLICIT REAL(16) (A-H,O-Z)
!    PARAMETER (TINY=1.0D-20)
!    DIMENSION A(NP,NP),INDX(N),VV(N)
!  
!           D=1.D0
!           DO 12 I=1,N
!            AAMAX=0.D0
!            DO 11 J=1,N
!               IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
! 11         CONTINUE
!            IF(AAMAX.EQ.0.) AAMX=1.d-22
!            VV(I)=1.D0/AAMAX
! 12           CONTINUE
!             DO 19 J=1,N
!            IF(J.GT.1) THEN
!               DO 14 I=1,J-1
!                  SUM=A(I,J)
!                  IF(I.GT.1) THEN
!                     DO 13 K=1,I-1
!                        SUM=SUM-A(I,K)*A(K,J)
! 13                  CONTINUE
!                     A(I,J)=SUM
!                  ENDIF
! 14            CONTINUE
!            ENDIF
!            AAMAX=0.D0
!            DO 16 I=J,N
!               SUM=A(I,J)
!               IF(J.GT.1) THEN
!                  DO 15 K=1,J-1
!                     SUM=SUM-A(I,K)*A(K,J)
! 15               CONTINUE
!                  A(I,J)=SUM
!               ENDIF
!               DUM=VV(I)*ABS(SUM)
!               IF(DUM.GE.AAMAX) THEN
!                  IMAX=I
!                  AAMAX=DUM
!               ENDIF
! 16         CONTINUE
!            IF(J.NE.IMAX) THEN
!               DO 17 K=1,N
!                  DUM=A(IMAX,K)
!                  A(IMAX,K)=A(J,K)
!                  A(J,K)=DUM
! 17            CONTINUE
!               D=-D
!               VV(IMAX)=VV(J)
!            ENDIF
!            INDX(J)=IMAX
!            IF(J.NE.N) THEN
!              IF(A(J,J).EQ.0.) A(J,J)=TINY
!              DUM=1./A(J,J)
!              DO 18 I=J+1,N
!                 A(I,J)=A(I,J)*DUM
! 18           CONTINUE
!            ENDIF
! 19          CONTINUE
!             IF(A(N,N).EQ.0.) A(N,N)=TINY
! 
!             RETURN
!             END subroutine

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

! Nom            : LUBKSB  
! Fonction       : Backsubstitution (pour l'inversion)
! Notes          : Provient de Numerical Recipes
! 
!       SUBROUTINE Q_LUBKSB(A,N,NP,INDX,B)
!       INTEGER N,NP,INDX(N)
!       REAL(16) A(NP,NP),B(N)
!       INTEGER I,II,J,LL
!       REAL(16) SUM
!       II=0
!       DO I=1,N
!          LL=INDX(I)
!          SUM=B(LL)
!          B(LL)=B(I)
!          IF(II.NE.0) THEN
!             DO J=II,I-1
!                SUM=SUM-A(I,J)*B(J)
!             END DO
!          ELSE IF(SUM.NE.0) THEN
!             II=I
!          END IF
!          B(I)=SUM
!       END DO
!       DO I=N,1,-1
!          SUM=B(I)
!          DO J=I+1,N
!             SUM=SUM-A(I,J)*B(J)
!          END DO
!          B(I)=SUM/A(I,I)
!       END DO
!       RETURN
!       END subroutine

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
! 
! !     Fonction:  Single value decomposition, QUAD precision
! 
!       SUBROUTINE Q_SVDCMP(A,M,N,MP,NP,W,V)
!       integer,intent(in)      :: M,N
!       real(16),intent(inout)  :: A(MP,NP),W(NP),V(NP,NP)
!       real(16)                :: RV1(N)
!       real(16)                :: S,G,H,F,SCALE,ANORM
!       INTEGER                 :: I,J,L,K
! 
!       G=0.0; SCALE=0.0; ANORM=0.0
! 
!       DO 25 I=1,N
!         L=I+1
!         RV1(I)=SCALE*G
!         G=0.0
!         S=0.0
!         SCALE=0.0
!         IF (I.LE.M) THEN
!           DO 11 K=I,M
!             SCALE=SCALE+ABS(A(K,I))
! 11        CONTINUE
!           IF (SCALE.NE.0.0) THEN
!             DO 12 K=I,M
!               A(K,I)=A(K,I)/SCALE
!               S=S+A(K,I)*A(K,I)
! 12          CONTINUE
!             F=A(I,I)
!             G=-SIGN(SQRT(S),F)
!             H=F*G-S
!             A(I,I)=F-G
!             IF (I.NE.N) THEN
!               DO 15 J=L,N
!                 S=0.0
!                 DO 13 K=I,M
!                   S=S+A(K,I)*A(K,J)
! 13              CONTINUE
!                 F=S/H
!                 DO 14 K=I,M
!                   A(K,J)=A(K,J)+F*A(K,I)
! 14              CONTINUE
! 15            CONTINUE
!             ENDIF
!             DO 16 K= I,M
!               A(K,I)=SCALE*A(K,I)
! 16          CONTINUE
!           ENDIF
!         ENDIF
!         W(I)=SCALE *G
!         G=0.0
!         S=0.0
!         SCALE=0.0
! 
!         IF ((I<=M).AND.(I/=N)) THEN
!           DO 17 K=L,N
!             SCALE=SCALE+ABS(A(I,K))
! 17        CONTINUE
!           IF (SCALE.NE.0.0) THEN
!             DO 18 K=L,N
!               A(I,K)=A(I,K)/SCALE
!               S=S+A(I,K)*A(I,K)
! 18          CONTINUE
!             F=A(I,L)
!             G=-SIGN(SQRT(S),F)
!             H=F*G-S
!             A(I,L)=F-G
!             DO 19 K=L,N
!               RV1(K)=A(I,K)/H
! 19          CONTINUE
!             IF (I.NE.M) THEN
!               DO 23 J=L,M
!                 S=0.0
!                 DO 21 K=L,N
!                   S=S+A(J,K)*A(I,K)
! 21              CONTINUE
!                 DO 22 K=L,N
!                   A(J,K)=A(J,K)+S*RV1(K)
! 22              CONTINUE
! 23            CONTINUE
!             ENDIF
!             DO 24 K=L,N
!               A(I,K)=SCALE*A(I,K)
! 24          CONTINUE
!           ENDIF
!         ENDIF
!         ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
! 25    CONTINUE
!       DO 32 I=N,1,-1
!         IF (I<N) THEN
!           IF (G/=0.0) THEN
!             DO 26 J=L,N
!               if(I>M) stop 'error pivot, bad shape 1'
!               V(J,I)=(A(I,J)/A(I,L))/G
! 26          CONTINUE
!             DO 29 J=L,N
!               S=0.0
!               DO 27 K=L,N
!                 if(I>M) stop 'error pivot, bad shape 2'
!                 S=S+A(I,K)*V(K,J)
! 27            CONTINUE
!               DO 28 K=L,N
!                 V(K,J)=V(K,J)+S*V(K,I)
! 28            CONTINUE
! 29          CONTINUE
!           ENDIF
!           DO 31 J=L,N
!             V(I,J)=0.0
!             V(J,I)=0.0
! 31        CONTINUE
!         ENDIF
!         V(I,I)=1.0
!         G=RV1(I)
!         L=I
! 32    CONTINUE
!       DO 39 I=N,1,-1
!         L=I+1
!         G=W(I)
!         IF (I.LT.N) THEN
!           DO 33 J=L,N
!             A(I,J)=0.0
! 33        CONTINUE
!         ENDIF
!         IF (G.NE.0.0) THEN
!           G=1.0/G
!           IF (I/=N) THEN
!             DO 36 J=L,N
!               S=0.0
!               DO 34 K=L,M
!                 S=S+A(K,I)*A(K,J)
! 34            CONTINUE
!               F=(S/A(I,I))*G
!               DO 35 K=I,M
!                 A(K,J)=A(K,J)+F*A(K,I)
! 35            CONTINUE
! 36          CONTINUE
!           ENDIF
!           DO 37 J=I,M
!             A(J,I)=A(J,I)*G
! 37        CONTINUE
!         ELSE
!           DO 38 J= I,M
!             A(J,I)=0.0
! 38        CONTINUE
!         ENDIF
!         A(I,I)=A(I,I)+1.0
! 39    CONTINUE
!       DO 49 K=N,1,-1
!         DO 48 ITS=1,30
!           DO 41 L=K,1,-1
!             NM=L-1
!             IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
!             IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
! 41        CONTINUE
! 1         C=0.0
!           S=1.0
!           DO 43 I=L,K
!             F=S*RV1(I)
!             IF ((ABS(F)+ANORM).NE.ANORM) THEN
!               G=W(I)
!               H=SQRT(F*F+G*G)
!               W(I)=H
!               H=1.0/H
!               C= (G*H)
!               S=-(F*H)
!               DO 42 J=1,M
!                 Y=A(J,NM)
!                 Z=A(J,I)
!                 A(J,NM)=(Y*C)+(Z*S)
!                 A(J,I)=-(Y*S)+(Z*C)
! 42            CONTINUE
!             ENDIF
! 43        CONTINUE
! 2         Z=W(K)
!           IF (L.EQ.K) THEN
!             IF (Z.LT.0.0) THEN
!               W(K)=-Z
!               DO 44 J=1,N
!                 V(J,K)=-V(J,K)
! 44            CONTINUE
!             ENDIF
!             GO TO 3
!           ENDIF
!           X=W(L)
!           NM=K-1
!           Y=W(NM)
!           G=RV1(NM)
!           H=RV1(K)
!           F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
!           G=SQRT(F*F+1.0)
!           F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
!           C=1.0
!           S=1.0
!           DO 47 J=L,NM
!             I=J+1
!             G=RV1(I)
!             Y=W(I)
!             H=S*G
!             G=C*G
!             Z=SQRT(F*F+H*H)
!             RV1(J)=Z
!             C=F/Z
!             S=H/Z
!             F= (X*C)+(G*S)
!             G=-(X*S)+(G*C)
!             H=Y*S
!             Y=Y*C
!             DO 45 NM=1,N
!               X=V(NM,J)
!               Z=V(NM,I)
!               V(NM,J)= (X*C)+(Z*S)
!               V(NM,I)=-(X*S)+(Z*C)
! 45          CONTINUE
!             Z=SQRT(F*F+H*H)
!             W(J)=Z
!             IF (Z.NE.0.0) THEN
!               Z=1.0/Z
!               C=F*Z
!               S=H*Z
!             ENDIF
!             F= (C*G)+(S*Y)
!             X=-(S*G)+(C*Y)
!             DO 46 NM=1,M
!               Y=A(NM,J)
!               Z=A(NM,I)
!               A(NM,J)= (Y*C)+(Z*S)
!               A(NM,I)=-(Y*S)+(Z*C)
! 46          CONTINUE
! 47        CONTINUE
!           RV1(L)=0.0
!           RV1(K)=F
!           W(K)=X
! 48      CONTINUE
! 3       CONTINUE
! 49    CONTINUE
! 
!       RETURN
!       END subroutine
! 
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
! 
!       subroutine q_zgedi(a,lda,n,ipvt,det,work,job)
!      !subroutines of LINPACK
!       implicit none
!       integer     :: lda,n,ipvt(n),job
!       complex(16) :: a(lda,n),det(2),work(n)
!       complex(16) :: t
!       integer     :: i,j,k,kb,kp1,l,nm1
! 
!          if (job/10 == 0) go to 70
!          det(1) = CMPLX(1.0,0.0,16)
!          det(2) = CMPLX(0.0,0.0,16)
!            
!          do 50 i = 1, n
!             if (ipvt(i) /= i) det(1) = -det(1)
!             det(1) = a(i,i)*det(1)
!             if (abs(det(1)) == 0.0d0) go to 60
!    10       if (abs(det(1)) >= 1.0d0) go to 20
!                det(1) = cmplx(10.0,0.0d0,kind=16)*det(1)
!                det(2) = det(2) - cmplx(1.0,0.0,16)
!             go to 10
!    20       continue
!    30       if (abs(dble(det(1)))<10.d0) go to 40
!                det(1) = det(1)  /cmplx(10.0,0.0,16)
!                det(2) = det(2) + cmplx(1.0d0,0.0d0,16)
!             go to 30
!    40       continue
!    50    continue
!    60    continue
!    70    continue
! 
!          if (mod(job,10) == 0) go to 150
! 
!          do k = 1, n
!             a(k,k) = cmplx(1.0,0.0,16)/a(k,k)
!             t = -a(k,k)
!             call q_zscal(k-1,t,a(1:k-1,k),1)
!             kp1 = k + 1
!             if (n .lt. kp1) go to 90
!             do j = kp1, n
!                t = a(k,j)
!                a(k,j) = cmplx(0.0,0.0,16)
!                call zaxpy__(k,t,a(1:k,k),1,a(1:k,j),1)
!             enddo 
!    90       continue
!          enddo 
!          nm1 = n - 1
!          if (nm1 .lt. 1) go to 140
!          do 130 kb = 1, nm1
!             k = n - kb
!             kp1 = k + 1
!             do 110 i = kp1, n
!                work(i) = a(i,k)
!                a(i,k) = cmplx(0.0,0.0,16)
!   110       continue
!             do j = kp1, n
!                t = work(j)
!                call zaxpy__(n,t,a(:,j),1,a(:,k),1)
!             enddo 
!             l = ipvt(k)
!             if (l/=k) call swap(n,a(:,k),1,a(:,l),1)
!   130    continue
!   140    continue
!   150 continue
! 
!       return
!       end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       subroutine q_zgeco(a,lda,n,ipvt,rcond,z)
!       implicit none
!       integer     :: lda,n,ipvt(n)
!       complex(16) :: a(lda,n),z(n)
!       real(16)    :: rcond
!       complex(16) :: ek,t,wk,wkm
!       real(16)    :: anorm,s,sm,ynorm,rrerror
!       integer     :: info,j,k,kb,kp1,l
! 
!       rrerror = 1.d-24
!       anorm = 0.0d0
! 
!       do j=1,n
!          anorm=max(anorm,q_zsum(n,a(:,j),1))
!       enddo
!       call q_zgefa(a,lda,n,ipvt,info)
!       ek =CMPLX(1.0d0,0.0d0,16)
!       z  =CMPLX(0.0d0,0.0d0,16)
! 
!       do k=1,n
!          if(abs(z(k)).ne.0.0d0) ek=csign1(ek,-z(k))
!          if(abs(ek-z(k)).le.abs(a(k,k))) go to 30
!             s = abs(a(k,k))/abs(ek-z(k))
!             z=z*s
!             ek = cmplx(s,0.0d0,kind=16)*ek
!    30    continue
!          wk  =  ek - z(k)
!          wkm = -ek - z(k)
!          s   =  abs(wk)
!          sm  =  abs(wkm)
!          if(abs(a(k,k))<rrerror) go to 40
!             wk = wk/conjg(a(k,k))
!             wkm = wkm/conjg(a(k,k))
!          go to 50
!    40    continue
!             wk  = cmplx(1.0d0,0.0d0,16)
!             wkm = cmplx(1.0d0,0.0d0,16)
!    50    continue
!          kp1 = k + 1
!          if (kp1 > n) go to 90
!             do 60 j = kp1, n
!                sm   = sm   + abs(z(j)+wkm*conjg(a(k,j)))
!                z(j) = z(j) + wk*conjg(a(k,j))
!                s    = s    + abs(z(j))
!    60       continue
!             if (s>=sm) go to 80
!                t  = wkm - wk
!                wk = wkm
!                do 70 j = kp1, n
!                   z(j) = z(j) + t*conjg(a(k,j))
!    70          continue
!    80       continue
!    90    continue
!          z(k) = wk
!       enddo
! 
!       s = 1.d0/q_zsum(n,z,1)
!       z=z*s
!       do 120 kb = 1, n
!          k = n + 1 - kb
!          if (k .lt. n) z(k) = z(k) + zdot(n-k,a(k+1:n,k),1,z(k+1:n),1)
!          if (abs(z(k)) .le. 1.0d0) go to 110
!             s = 1.0d0/abs(z(k))
!             z=z*s
!   110    continue
!          l = ipvt(k)
!          t = z(l)
!          z(l) = z(k)
!          z(k) = t
!   120 continue
! 
!       s = 1.d0/q_zsum(n,z,1)
! 
!       z=z*s
!       ynorm = 1.0d0
!       do 140 k = 1, n
!          l = ipvt(k)
!          t = z(l)
!          z(l) = z(k)
!          z(k) = t
!          if (k .lt. n) call zaxpy__(n-k,t,a(k+1:n,k),1,z(k+1:n),1)
!          if (abs(z(k)) .le. 1.0d0) go to 130
!             s = 1.0d0/abs(z(k))
!             z=z*s
!             ynorm = s*ynorm
!   130    continue
!   140 continue
! 
!       s = 1.d0/q_zsum(n,z,1)
! 
!       z=z*s
!       ynorm = s*ynorm
!       do 160 kb = 1, n
!          k = n + 1 - kb
!          if (abs(z(k)) .le. abs(a(k,k))) go to 150
!             s = abs(a(k,k))/abs(z(k))
!             z=z*s
!             ynorm = s*ynorm
!   150    continue
!          if (abs(a(k,k)) > rrerror) then
!           z(k) = z(k)/a(k,k)
!          else
!           z(k) = CMPLX(1.0d0,0.0d0,16)
!          endif
!          t = -z(k)
!          call zaxpy__(k-1,t,a(1:k-1,k),1,z(1:k-1),1)
!   160 continue
! 
!       s = 1.d0/q_zsum(n,z,1)
! 
!       z=z*s
!       ynorm = s*ynorm
!       if(abs(anorm) > rrerror) then
!        rcond = ynorm/anorm
!       else
!        rcond = 1.d0/rrerror
!       endif
!       
!       return
!       end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       subroutine q_zgefa(a,lda,n,ipvt,info)
!       implicit none
!       integer     :: lda,n,ipvt(n),info
!       complex(16) :: a(lda,n)
!       complex(16) :: t
!       integer     :: j,k,kp1,l,nm1
! 
!       info = 0
!       nm1 = n - 1
!       if (nm1 < 1) go to 70
!       do 60 k = 1, nm1
!          kp1 = k + 1
!          l = q_zamax(n-k+1,a(k,k),1) + k - 1
!          ipvt(k) = l
!          if (abs(a(l,k)) == 0.0d0) go to 40
!             if (l == k) go to 10
!                t = a(l,k)
!                a(l,k) = a(k,k)
!                a(k,k) = t
!    10       continue
!             t = -CMPLX(1.0d0,0.0d0,16)/a(k,k)
!             call q_zscal(n-k,t,a(k+1,k),1)
!             do 30 j = kp1, n
!                t = a(l,j)
!                if (l == k) go to 20
!                   a(l,j) = a(k,j)
!                   a(k,j) = t
!    20          continue
!                call zaxpy__(n-k,t,a(k+1:n,k),1,a(k+1:n,j),1)
!    30       continue
!          go to 50
!    40    continue
!             info = k
!    50    continue
!    60 continue
!    70 continue
!       ipvt(n) = n
!       if (abs(a(n,n)) == 0.0d0) info = n
!       return
!       end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
! 
! SUBROUTINE inverses_(nn,a)
!   IMPLICIT NONE
!   REAL(4), DIMENSION(:,:), INTENT(INOUT) :: a
!   INTEGER, DIMENSION(SIZE(a,1))          :: ipiv,indxr,indxc
!   INTEGER                                :: nn
!   LOGICAL, DIMENSION(SIZE(a,1))          :: lpiv
!   REAL(4)                                :: pivinv
!   REAL(4), DIMENSION(SIZE(a,1))          :: dumc
!   INTEGER, TARGET                        :: irc(2)
!   INTEGER                                :: i,l,n
!   INTEGER, POINTER                       :: irow,icol
!   n=SIZE(a,1)
!   irow => irc(1)
!   icol => irc(2)
!   ipiv=0
!   DO i=1,n
!      lpiv = (ipiv == 0)
!      irc=MAXLOC(ABS(a),outerand(lpiv,lpiv))
!      ipiv(icol)=ipiv(icol)+1
!      IF (ipiv(icol) > 1) STOP 'gaussj:singular matrix (1)'
!      IF (irow /= icol) CALL swap(a(irow,:),a(icol,:))
!      indxr(i)=irow !We are now ready to divide the pivot row by the pivot element,
!                    !located at irow and icol.
!      indxc(i)=icol
!      IF (a(icol,icol) == zero) STOP 'gaussj:singular matrix (2)'
!      pivinv=one/a(icol,icol)
!      a(icol,icol)=CMPLX(one,zero)
!      a(icol,:)=a(icol,:)*pivinv
!      dumc=a(:,icol)
!      a(:,icol)     = CMPLX(zero,zero)
!      a(icol,icol)  = pivinv
!      a(1:icol-1,:) = a(1:icol-1,:) - outerprod(dumc(1:icol-1),a(icol,:))
!      a(icol+1:,:)  = a(icol+1:,:)  - outerprod(dumc(icol+1:),a(icol,:))
!   END DO
!   DO l=n,1,-1
!      CALL swap(a(:,indxr(l)),a(:,indxc(l)))
!   END DO
! CONTAINS
!  FUNCTION outerand(a,b)
!    IMPLICIT NONE
!    LOGICAL, DIMENSION(:), INTENT(IN)   :: a,b
!    LOGICAL, DIMENSION(SIZE(a),SIZE(b)) :: outerand
!    outerand = SPREAD(a,dim=2,ncopies=SIZE(b)).AND.SPREAD(b,dim=1,ncopies=SIZE(a))
!  END FUNCTION
! 
! end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! SUBROUTINE inverses__(nn,a) 
!   IMPLICIT NONE 
!   COMPLEX(4), DIMENSION(:,:), INTENT(INOUT) :: a
!   INTEGER, DIMENSION(SIZE(a,1))             :: ipiv,indxr,indxc 
!   INTEGER                                   :: nn
!   LOGICAL, DIMENSION(SIZE(a,1))             :: lpiv 
!   COMPLEX(4)                                :: pivinv 
!   COMPLEX(4), DIMENSION(SIZE(a,1))          :: dumc 
!   INTEGER, TARGET                           :: irc(2) 
!   INTEGER                                   :: i,l,n
!   INTEGER, POINTER                          :: irow,icol 
!   n=SIZE(a,1)
!   irow => irc(1) 
!   icol => irc(2) 
!   ipiv=0 
!   DO i=1,n 
!      lpiv = (ipiv == 0) 
!      irc=MAXLOC(ABS(a),outerand(lpiv,lpiv))
!      ipiv(icol)=ipiv(icol)+1 
!      IF (ipiv(icol) > 1) STOP 'gaussj:singular matrix (1)'
!      IF (irow /= icol) CALL swap(a(irow,:),a(icol,:)) 
!      indxr(i)=irow !We are now ready to divide the pivot row by the pivot element, 
!                    !located at irow and icol.
!      indxc(i)=icol 
!      IF (a(icol,icol) == zero) STOP 'gaussj:singular matrix (2)'
!      pivinv=one/a(icol,icol) 
!      a(icol,icol)=CMPLX(one,zero)
!      a(icol,:)=a(icol,:)*pivinv 
!      dumc=a(:,icol)
!      a(:,icol)     = CMPLX(zero,zero)
!      a(icol,icol)  = pivinv 
!      a(1:icol-1,:) = a(1:icol-1,:) - outerprod(dumc(1:icol-1),a(icol,:)) 
!      a(icol+1:,:)  = a(icol+1:,:)  - outerprod(dumc(icol+1:),a(icol,:)) 
!   END DO
!   DO l=n,1,-1 
!      CALL swap(a(:,indxr(l)),a(:,indxc(l))) 
!   END DO
! 
! CONTAINS
! 
!  FUNCTION outerand(a,b) 
!    IMPLICIT NONE
!    LOGICAL, DIMENSION(:), INTENT(IN)   :: a,b 
!    LOGICAL, DIMENSION(SIZE(a),SIZE(b)) :: outerand 
!    outerand = SPREAD(a,dim=2,ncopies=SIZE(b)).AND.SPREAD(b,dim=1,ncopies=SIZE(a)) 
!  END FUNCTION
! 
! END SUBROUTINE
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! SUBROUTINE inverse_(nn,a) 
!   IMPLICIT NONE 
!   COMPLEX(8), DIMENSION(:,:), INTENT(INOUT) :: a
!  
!   !--------------------------------------------------------------------------! 
!   ! Linear equation solution by Gauss-Jordan elimination                     !
!   ! a is an N x N input coefficient matrix. On output, a is replaced by its  !
!   ! matrix inverse.                                                          !
!   !--------------------------------------------------------------------------!
!   
!   INTEGER, DIMENSION(SIZE(a,1))      :: ipiv,indxr,indxc 
! 
!   !These arrays are used for bookkeeping on the pivoting. 
! 
!   INTEGER                            :: nn
!   LOGICAL, DIMENSION(SIZE(a,1))      :: lpiv 
!   COMPLEX(8)                         :: pivinv 
!   COMPLEX(8), DIMENSION(SIZE(a,1))   :: dumc 
!   INTEGER, TARGET                    :: irc(2) 
!   INTEGER                            :: i,l,n
!   INTEGER, POINTER                   :: irow,icol 
! 
!   n=SIZE(a,1)
! 
!   irow => irc(1) 
!   icol => irc(2) 
! 
!   ipiv=0 
!   
!   DO i=1,n 
!      !Main loop over columns to be reduced. 
!      lpiv = (ipiv == 0) 
!      !Begin search for a pivot element. 
!      irc=MAXLOC(ABS(a),outerand(lpiv,lpiv))
!      ipiv(icol)=ipiv(icol)+1 
!      IF (ipiv(icol) > 1) STOP 'gaussj:singular matrix (1)'
!      
!      !We now have the pivot element, so we interchange
!      !rows, if needed, to put the pivot element on the diagonal. The columns
!      !are not physically interchanged, only relabeled:
!      !indxc(i),the column of the ith pivot element, is the ith column that is
!      !reduced, while indxr(i) is the row in which that pivot element was
!      !originally located. If indxr(i) = indxc(i) there is an implied column
!      !interchange. With this form of bookkeeping, the inverse matrix will be
!      !scrambled by
!      !columns. 
! 
!      IF (irow /= icol) CALL swap(a(irow,:),a(icol,:)) 
! 
!      indxr(i)=irow !We are now ready to divide the pivot row by the pivot element, 
!                    !located at irow and icol.
!      indxc(i)=icol 
!      
!      IF (a(icol,icol) == zero) STOP 'gaussj:singular matrix (2)'
!      pivinv=one/a(icol,icol) 
!      a(icol,icol)=CMPLX(one,zero)
!      a(icol,:)=a(icol,:)*pivinv 
!      dumc=a(:,icol)
! 
!      !Next, we reduce the rows, except for the pivot one, of course. 
!      a(:,icol)     = CMPLX(zero,zero)
!      a(icol,icol)  = pivinv 
!      a(1:icol-1,:) = a(1:icol-1,:) - outerprod(dumc(1:icol-1),a(icol,:)) 
!      a(icol+1:,:)  = a(icol+1:,:)  - outerprod(dumc(icol+1:),a(icol,:)) 
!   END DO
! 
!   !It only remains to unscramble the solution in view of the column
!   !interchanges. 
!   !We do this by interchanging pairs of columns in the reverse order that the
!   !permutation 
!   !was built up. 
! 
!   DO l=n,1,-1 
!      CALL swap(a(:,indxr(l)),a(:,indxc(l))) 
!   END DO
! 
! CONTAINS
! 
!  FUNCTION outerand(a,b) 
!    IMPLICIT NONE
!    LOGICAL, DIMENSION(:), INTENT(IN)   :: a,b 
!    LOGICAL, DIMENSION(SIZE(a),SIZE(b)) :: outerand 
!    outerand = SPREAD(a,dim=2,ncopies=SIZE(b)).AND.SPREAD(b,dim=1,ncopies=SIZE(a)) 
!  END FUNCTION
! 
! END SUBROUTINE
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! SUBROUTINE inverse__(nn,a) 
!   IMPLICIT NONE 
!   REAL(8), DIMENSION(:,:), INTENT(INOUT) :: a
!   INTEGER, DIMENSION(SIZE(a,1))          :: ipiv,indxr,indxc 
!   INTEGER                                :: nn
!   LOGICAL, DIMENSION(SIZE(a,1))          :: lpiv 
!   REAL(8)                                :: pivinv 
!   REAL(8), DIMENSION(SIZE(a,1))          :: dumc 
!   INTEGER, TARGET                        :: irc(2) 
!   INTEGER                                :: i,l,n
!   INTEGER, POINTER                       :: irow,icol 
! 
!   n=SIZE(a,1)
! 
!   irow => irc(1) 
!   icol => irc(2) 
! 
!   ipiv=0 
!   
!   DO i=1,n 
!      !Main loop over columns to be reduced. 
!      lpiv = (ipiv == 0) 
!      !Begin search for a pivot element. 
!      irc=MAXLOC(ABS(a),outerand(lpiv,lpiv))
!      ipiv(icol)=ipiv(icol)+1 
!      IF (ipiv(icol) > 1) STOP 'gaussj:singular matrix (1)'
!      
!      !We now have the pivot element, so we interchange
!      !rows, if needed, to put the pivot element on the diagonal. The columns
!      !are not physically interchanged, only relabeled:
!      !indxc(i),the column of the ith pivot element, is the ith column that is
!      !reduced, while indxr(i) is the row in which that pivot element was
!      !originally located. If indxr(i) = indxc(i) there is an implied column
!      !interchange. With this form of bookkeeping, the inverse matrix will be
!      !scrambled by
!      !columns. 
! 
!      IF (irow /= icol) CALL swap(a(irow,:),a(icol,:)) 
! 
!      indxr(i)=irow !We are now ready to divide the pivot row by the pivot element, 
!                    !located at irow and icol.
!      indxc(i)=icol 
!      
!      IF (a(icol,icol) == zero) STOP 'gaussj:singular matrix (2)'
!      pivinv=one/a(icol,icol) 
!      a(icol,icol)=CMPLX(one,zero)
!      a(icol,:)=a(icol,:)*pivinv 
!      dumc=a(:,icol)
! 
!      !Next, we reduce the rows, except for the pivot one, of course. 
!      a(:,icol)     = CMPLX(zero,zero)
!      a(icol,icol)  = pivinv 
!      a(1:icol-1,:) = a(1:icol-1,:) - outerprod(dumc(1:icol-1),a(icol,:)) 
!      a(icol+1:,:)  = a(icol+1:,:)  - outerprod(dumc(icol+1:),a(icol,:)) 
!   END DO
! 
!   !It only remains to unscramble the solution in view of the column
!   !interchanges. 
!   !We do this by interchanging pairs of columns in the reverse order that the
!   !permutation 
!   !was built up. 
!   DO l=n,1,-1 
!      CALL swap(a(:,indxr(l)),a(:,indxc(l))) 
!   END DO
! 
! CONTAINS
! 
!  FUNCTION outerand(a,b) 
!    IMPLICIT NONE
!    LOGICAL, DIMENSION(:), INTENT(IN)   :: a,b 
!    LOGICAL, DIMENSION(SIZE(a),SIZE(b)) :: outerand 
!    outerand = SPREAD(a,dim=2,ncopies=SIZE(b)).AND.SPREAD(b,dim=1,ncopies=SIZE(a)) 
!  END FUNCTION
! 
! END SUBROUTINE

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
! 
!   subroutine test_hermicity(Amat)
!   logical    :: not_hermitian
!   COMPLEX(8) :: Amat(:,:)
!      not_hermitian = ANY( abs(Amat-TRANSPOSE(CONJG(Amat)))>1.d-8 )
!      IF(not_hermitian)THEN
!          CALL dump_message(TEXT="ERROR MATRICE ISN T HERMITIC!")
!          CALL write_cplx_array_rank_2(Amat,"Amat = ")
!          STOP 'matrix not hermitian'
!      ENDIF
!   end subroutine
! 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE diagonalize_comp(matrix,VALP,VECP,EIGENVAL_ONLY)
    implicit none 
  !--------------------------------------------------------------------------------------! 
  ! CONVENIENT WRAPPER OF DSYEVD/ZHEEVD (DIVIDE & CONQUER OF SYMMETRIC/HERMITIAN MATRIX) !
  !--------------------------------------------------------------------------------------!

    COMPLEX(8)                         :: matrix(:,:)
    COMPLEX(8),       INTENT(INOUT)    :: VECP(:,:)
    REAL(8),          INTENT(INOUT)    :: VALP(:)
    LOGICAL,  OPTIONAL, INTENT(IN)     :: EIGENVAL_ONLY
    INTEGER                            :: LRWORK
    CHARACTER(LEN=1)                   :: FLAG 
    INTEGER                            :: LWORK,LIWORK,N,info
    complex(8), allocatable            :: work(:)
    real(8), allocatable               :: rwork(:)
    INTEGER, allocatable               :: iwork(:)

    N = size(VALP)

    write(*,*) '......... DIAGONALIZE complex matrix ..........'
    if(any(shape(matrix)-shape(VECP)>0)) then
      write(*,*) 'error diagonalize_comp shape do not match'
      write(*,*) 'shape matrix : ', shape(matrix)
      write(*,*) 'shape VECP   : ', shape(VECP)
      stop 'critical'
    endif
    IF(N/=SIZE(matrix,1).OR.N/=SIZE(matrix,2).OR.N/=SIZE(VECP,1).OR.N/=SIZE(VECP,2)) &
   & STOP "ERROR IN diagonalize_dsyevd: INCONSISTENT DIMENSIONS!"

                        FLAG = 'V' 
    IF(PRESENT(EIGENVAL_ONLY))THEN
      IF(EIGENVAL_ONLY) FLAG = 'N'
    ENDIF

    IF(N<=1)THEN
      LWORK  = 1
      LIWORK = 1
      LRWORK = 1
    ELSE
      SELECT CASE(FLAG)
        CASE('N')
          LWORK  = N + 1
          LIWORK = 1
          LRWORK = N 
        CASE('V')
          LWORK  = 2*N + N*N + 5
          LRWORK = 5*N + 2*N*N + 5
          LIWORK = 3 + 5*N + 5 
      END SELECT
    ENDIF

    ALLOCATE(work(LWORK), rwork(LRWORK), iwork(LIWORK))
    VECP=matrix; 
    write(*,*) '.... call ZHEEVD ....',FLAG
    CALL ZHEEVD(FLAG,'U',N,VECP,N,VALP,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,info)
    DEALLOCATE(work, rwork, iwork)

    IF(info>0)THEN
      SELECT CASE(FLAG)
        CASE('N')
          CALL dump_message(TEXT="ERROR IN diagonalize: FAILED TO CONVERGE")
          CALL dump_message(TEXT=c2s(i2c(info))//" ELEMENTS OF AN INTERMEDIATE TRIDIAGONAL FORM DIDN'T CONVERGE TO 0!")
        CASE('V')
          CALL dump_message(TEXT="ERROR IN diagonalize: FAILED TO COMPUTE AN EIGENVALUE")
          CALL dump_message(TEXT="WHILE WORKING ON SUBMATRIX "//c2s(i2c(info/(N+1)))//" x "//c2s(i2c(MOD(info,N+1))))
      END SELECT
    ENDIF
    IF(info<0) CALL dump_message(TEXT=c2s(i2c(info))//"-TH ARGUMENT HAD AN ILLEGAL VALUE")


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

  SUBROUTINE diagonalize_real(matrix,VALP,VECP,EIGENVAL_ONLY)
    implicit none

    !--------------------------------------------------------------------------------------!
    ! CONVENIENT WRAPPER OF DSYEVD/ZHEEVD (DIVIDE & CONQUER OF SYMMETRIC/HERMITIAN MATRIX) !
    !--------------------------------------------------------------------------------------!

    REAL(8)                         :: matrix(:,:)
    REAL(8),          INTENT(INOUT) :: VECP(:,:)
    REAL(8),          INTENT(INOUT) :: VALP(:)
    LOGICAL,  OPTIONAL, INTENT(IN)  :: EIGENVAL_ONLY
    CHARACTER(LEN=1)                :: FLAG 
    INTEGER                         :: LWORK,LIWORK,N,info
    real(8), allocatable            :: work(:)
    INTEGER, allocatable            :: iwork(:)

    N = SIZE(VALP)
    IF(N/=SIZE(matrix,1).OR.N/=SIZE(matrix,2).OR.N/=SIZE(VECP,1).OR.N/=SIZE(VECP,2)) &
        &  STOP "ERROR IN diagonalize_dsyevd: INCONSISTENT DIMENSIONS!"

    if(any(shape(matrix)-shape(VECP)>0)) then
        write(*,*) 'error diagonalize_comp shape do not match'
        write(*,*) 'shape matrix : ', shape(matrix)
        write(*,*) 'shape VECP   : ', shape(VECP)
        stop 'critical'
    endif

                        FLAG = 'V' 
    IF(PRESENT(EIGENVAL_ONLY))THEN
      IF(EIGENVAL_ONLY) FLAG = 'N'
    ENDIF

    IF(N<=1)THEN
      LWORK  = 1
      LIWORK = 1
    ELSE
      SELECT CASE(FLAG)
        CASE('N')
          LWORK  = 2 * N + 1
          LIWORK = 1
        CASE('V')
          LWORK  = 2 * N**2 + 6 * N + 1
          LIWORK =            5 * N + 3
      END SELECT
    ENDIF

    ALLOCATE(work(LWORK),iwork(LIWORK))
    VECP=matrix
    write(*,*) '.... call DSYEVD ....', FLAG
    CALL DSYEVD(FLAG,'U',N,VECP,N,VALP,WORK,LWORK,IWORK,LIWORK,info)
    DEALLOCATE(work,iwork)

    IF(info>0)THEN
      SELECT CASE(FLAG)
        CASE('N')
          CALL dump_message(TEXT="ERROR IN diagonalize: FAILED TO CONVERGE")
          CALL dump_message(TEXT=c2s(i2c(info))//" ELEMENTS OF AN INTERMEDIATE TRIDIAGONAL FORM DIDN'T CONVERGE TO 0!")
        CASE('V')
          CALL dump_message(TEXT="ERROR IN diagonalize: FAILED TO COMPUTE AN EIGENVALUE")
          CALL dump_message(TEXT="WHILE WORKING ON SUBMATRIX "//c2s(i2c(info/(N+1)))//" x "//c2s(i2c(MOD(info,N+1))))
      END SELECT
    ENDIF
    IF(info<0) CALL dump_message(TEXT=c2s(i2c(info))//"-TH ARGUMENT HAD AN ILLEGAL VALUE")


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
! 
!   SUBROUTINE invert_lapack_c(A)
!     COMPLEX(8), INTENT(INOUT) :: A(:,:)
!     INTEGER                   :: bIPVT(size(A,1))
!     IF(SIZE(A,1)/=SIZE(A,2)) STOP "ERROR IN invert: SQUARE MATRIX REQUIRED!"
!     CALL GETRF__(A,bIPVT)
!     CALL GETRI__(A,bIPVT)
!   RETURN
!   END SUBROUTINE
! 
! 
!       !---------------------------------------------------------!
! 
!   SUBROUTINE invert_lapack_r(A)
!     REAL(8), INTENT(INOUT)    :: A(:,:)
!     INTEGER                   :: bIPVT(size(A,1))
!     IF(SIZE(A,1)/=SIZE(A,2)) STOP "ERROR IN invert: SQUARE MATRIX REQUIRED!"
!     CALL GETRF__(A,bIPVT)
!     CALL GETRI__(A,bIPVT)
!   RETURN
!   END SUBROUTINE
! 
! 
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE new_Id(Id,N)
    REAL(8), INTENT(INOUT)            :: Id(:,:)
    INTEGER,   INTENT(IN),   OPTIONAL :: N
    INTEGER :: i,N_
    N_ = SIZE(Id,1)
    IF(PRESENT(N)) N_ = N
    Id = 0.0_DBL
    DO i=1,N_
      Id(i,i) = 1.0_DBL
    ENDDO
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
! 
!   SUBROUTINE new_cId(cId,N) 
!     INTEGER,      INTENT(IN)    :: N
!     COMPLEX(8), INTENT(INOUT)   :: cId(:,:)
!     INTEGER                     :: i
!     cId = zero
!     DO i=1,N
!       cId(i,i) =  one
!     ENDDO
!   END SUBROUTINE 
! 
!        !--------------!
! 
  SUBROUTINE new_diag(Id,N) 
    LOGICAL, INTENT(INOUT) :: Id(:,:)
    INTEGER, INTENT(IN)    :: N
    INTEGER                :: i
    Id = .false.
    DO i=1,N
     Id(i,i) = .true.
    ENDDO
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

  SUBROUTINE write_cplx_array_rank_3(A,title,SHORT,UNIT)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ WRITE COMPLEX ARRAY A(n1,n2,n3) $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    CHARACTER(LEN=*),               INTENT(IN) :: title
    COMPLEX(8), DIMENSION(:,:,:), INTENT(IN)   :: A
    LOGICAL,              INTENT(IN), OPTIONAL :: SHORT
    INTEGER,              INTENT(IN), OPTIONAL :: UNIT
    INTEGER                                    :: i1,i2,i3,n1,n2,n3,unit_
    LOGICAL                                    :: short_
    CHARACTER(LEN=400)                         :: fmt_A

                      unit_ = log_unit ! DEFAULT: STANDARD OUTPUT
    IF(PRESENT(UNIT)) unit_ = UNIT

    short_ = present(short)

    n1 = SIZE(A,1); n2 = SIZE(A,2); n3 = SIZE(A,3)

    CALL dump_message(UNIT=unit_,TEXT=title)

    IF(short_)THEN
      WRITE(fmt_A,*) '(',n1-1,'(',n3,'(',n2,'(2(a,f10.6),a,x),2x)/),',n3,'(',n2,'(2(a,f10.6),a,x),2x))'
    ELSE
      WRITE(fmt_A,*) '(',n1-1,'(',n3,'(',n2,'(2(a,f20.16),a,x),2x)/),',n3,'(',n2,'(2(a,f20.16),a,x),2x))'
    ENDIF

    WRITE(unit_,fmt_A) ((('(',DBLE(A(i1,i2,i3)),',',AIMAG(A(i1,i2,i3)),')',i2=1,n2),i3=1,n3),i1=1,n1)

    CALL flush(unit_)

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

  SUBROUTINE write_real_array_rank_3(A,title,SHORT,UNIT)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ WRITE REAL ARRAY A(n1,n2,n3) $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    CHARACTER(LEN=*),               INTENT(IN) :: title
    REAL(8),    DIMENSION(:,:,:), INTENT(IN)   :: A
    LOGICAL,              INTENT(IN), OPTIONAL :: SHORT
    INTEGER,              INTENT(IN), OPTIONAL :: UNIT
    INTEGER                                    :: i1,i2,i3,n1,n2,n3,unit_
    LOGICAL                                    :: short_
    CHARACTER(LEN=400)                         :: fmt_A

    unit_ = log_unit ! DEFAULT: STANDARD OUTPUT
    IF(PRESENT(UNIT)) unit_ = UNIT

    short_ = present(short) 

    n1 = SIZE(A,1)
    n2 = SIZE(A,2)
    n3 = SIZE(A,3)

    CALL dump_message(UNIT=unit_,TEXT=title)

    IF(short_)THEN
      WRITE(fmt_A,*) '(',n1-1,'(',n3,'(',n2,'(f10.6,x),2x)/),',n3,'(',n2,'(f10.6,x),2x))'
    ELSE
      WRITE(fmt_A,*) '(',n1-1,'(',n3,'(',n2,'(f20.16,x),2x)/),',n3,'(',n2,'(f20.16,x),2x))'
    ENDIF
    WRITE(unit_,fmt_A) (((A(i1,i2,i3),i2=1,n2),i3=1,n3),i1=1,n1)

    CALL flush(unit_)
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


  SUBROUTINE write_intg_array_rank_3(A,title,UNIT)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ WRITE INTEGER ARRAY A(n1,n2,n3) $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    CHARACTER(LEN=*),               INTENT(IN) :: title
    INTEGER,      DIMENSION(:,:,:), INTENT(IN) :: A
    INTEGER,              INTENT(IN), OPTIONAL :: UNIT
    INTEGER                                    :: i1,i2,i3,n1,n2,n3,unit_
    CHARACTER(LEN=1000)                        :: fmt_A

    unit_ = log_unit ! DEFAULT: STANDARD OUTPUT
    IF(PRESENT(UNIT)) unit_ = UNIT

    n1 = SIZE(A,1)
    n2 = SIZE(A,2)
    n3 = SIZE(A,3)

    CALL dump_message(UNIT=unit_,TEXT=title)

    WRITE(fmt_A,*) '(',n1-1,'(2x,',n3,'(',n2,'(I3,x),2x)/),2x,',n3,'(',n2,'(I3,x),2x))'
    WRITE(unit_,fmt_A) (((A(i1,i2,i3),i2=1,n2),i3=1,n3),i1=1,n1)

    CALL flush(unit_)
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

  SUBROUTINE write_bool_array_rank_3(A,title,UNIT)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ WRITE BOOLEAN ARRAY A(n1,n2,n3) $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    CHARACTER(LEN=*),               INTENT(IN)  :: title
    LOGICAL,      DIMENSION(:,:,:), INTENT(IN)  :: A
    INTEGER,              INTENT(IN), OPTIONAL  :: UNIT
    INTEGER                                     :: i1,i2,i3,n1,n2,n3,unit_
    CHARACTER(LEN=400)                          :: fmt_A

    unit_ = log_unit ! DEFAULT: STANDARD OUTPUT
    IF(PRESENT(UNIT)) unit_ = UNIT

    n1 = SIZE(A,1)
    n2 = SIZE(A,2)
    n3 = SIZE(A,3)

    CALL dump_message(UNIT=unit_,TEXT=title)

    WRITE(fmt_A,*) '(',n1-1,'(2x,',n3,'(',n2,'(L2,x),2x)/),',n3,'(2x,',n2,'(L2,x)))'
    WRITE(unit_,fmt_A) (((A(i1,i2,i3),i2=1,n2),i3=1,n3),i1=1,n1)

    CALL flush(unit_)
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


  SUBROUTINE write_cplx_array_rank_2(A,title,UNIT,SHORT,ULTRASHORT)
    use common_def, only : dump_message

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ WRITE COMPLEX ARRAY A(n1,n2) $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    CHARACTER(LEN=*), INTENT(IN)           :: title
    COMPLEX(8),     INTENT(IN)             :: A(:,:)
    LOGICAL,          INTENT(IN), OPTIONAL :: SHORT,ULTRASHORT
    INTEGER,          INTENT(IN), OPTIONAL :: UNIT
    INTEGER                                :: i1,i2,n1,n2,unit_
    LOGICAL                                :: short_
    CHARACTER(LEN=400)                     :: fmt_A

    unit_ = log_unit ! DEFAULT: STANDARD OUTPUT
    IF(PRESENT(UNIT)) unit_ = UNIT

    short_ = present(short) 

    n1 = SIZE(A,1); n2 = SIZE(A,2)

    CALL dump_message(UNIT=unit_,TEXT=title)

    if(.not.present(ULTRASHORT))then
     IF(short_)THEN
      WRITE(fmt_A,*) '(',n2,'( 2f10.6))'
     ELSE
      WRITE(fmt_A,*) '(',n2,'( 2f20.10))'
     ENDIF
    else
      WRITE(fmt_A,*) '(',n2,'( 2f5.1))'
    endif

    WRITE(unit_,fmt_A,err=10) ((DBLE(A(i1,i2)),AIMAG(A(i1,i2)),i2=1,n2),i1=1,n1)
 10 continue
    CALL flush(unit_)

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

  SUBROUTINE write_real_array_rank_2(A,title,UNIT,SHORT,ULTRASHORT)
    use common_def, only : dump_message

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ WRITE REAL ARRAY A(n1,n2) $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    CHARACTER(LEN=*),             INTENT(IN)           :: title
    REAL(8),    DIMENSION(:,:), INTENT(IN)             :: A
    LOGICAL,                      INTENT(IN), OPTIONAL :: SHORT,ULTRASHORT
    INTEGER,                      INTENT(IN), OPTIONAL :: UNIT
    INTEGER                                            :: i1,i2,n1,n2,unit_
    LOGICAL                                            :: short_
    CHARACTER(LEN=400)                                 :: fmt_A

    unit_ = log_unit ! STANDARD OUTPUT
    IF(PRESENT(UNIT)) unit_ = UNIT

    short_ = present(short)

    n1 = SIZE(A,1)
    n2 = SIZE(A,2)
    if(n1==1.and.n2==1) return

    CALL dump_message(UNIT=unit_,TEXT=title)

   if(.not.present(ULTRASHORT))then
    IF(short_)THEN
      WRITE(fmt_A,*) '(',n1-1,'(',n2,'(f10.6,x)/),',n2,'(f10.6,x))'
    ELSE
      WRITE(fmt_A,*) '(',n1-1,'(',n2,'(f20.16,x)/),',n2,'(f20.16,x))'
    ENDIF
   else
      WRITE(fmt_A,*) '(',n1-1,'(',n2,'(f5.1,x)/),',n2,'(f5.1,x))'
   endif

    WRITE(unit_,fmt_A) ((A(i1,i2),i2=1,n2),i1=1,n1)

    CALL flush(unit_)

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


  SUBROUTINE write_intg_array_rank_2(A,title,UNIT)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ WRITE INTEGER ARRAY A(n1,n2) $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    CHARACTER(LEN=*),  INTENT(IN) :: title
    INTEGER,           INTENT(IN) :: A(:,:)
    INTEGER, OPTIONAL, INTENT(IN) :: UNIT
    INTEGER                       :: i1,i2,n1,n2,unit_
    CHARACTER(LEN=400)            :: fmt_A

    unit_ = log_unit ! STANDARD OUTPUT
    IF(PRESENT(UNIT)) unit_ = UNIT

    n1 = SIZE(A,1)
    n2 = SIZE(A,2)

    CALL dump_message(UNIT=unit_,TEXT=title)

!CEDRIC NOVEMBER 2011
    if(n1>1.and.n2>1)then
       WRITE(fmt_A,*) '(',n1-1,'(2x,',n2,'(I4,x)/),2x,',n2,'(I4,x))'
       WRITE(unit_,fmt_A) ((A(i1,i2),i2=1,n2),i1=1,n1)
    else
       WRITE(unit_,*) ((A(i1,i2),i2=1,n2),i1=1,n1)
    endif

    CALL flush(unit_)

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


  SUBROUTINE write_bool_array_rank_2(A,title,UNIT)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ WRITE BOOLEAN ARRAY A(n1,n2) $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    CHARACTER(LEN=*),        INTENT(IN) :: title
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: A
    INTEGER,       INTENT(IN), OPTIONAL :: UNIT
    INTEGER                             :: i1,i2,n1,n2,unit_
    CHARACTER(LEN=400)                  :: fmt_A

    unit_ = log_unit ! DEFAULT: STANDARD OUTPUT
    IF(PRESENT(UNIT)) unit_ = UNIT

    n1 = SIZE(A,1)
    n2 = SIZE(A,2)

    CALL dump_message(UNIT=unit_,TEXT=title)

!CEDRIC NOVEMBER 2011
    if(n1>1.and.n2>1)then
       WRITE(fmt_A,*) '(',n1-1,'(2x,',n2,'(L2,x)/),2x,',n2,'(L2,x))'
       WRITE(unit_,fmt_A) ((A(i1,i2),i2=1,n2),i1=1,n1)
    else
       WRITE(unit_,*) ((A(i1,i2),i2=1,n2),i1=1,n1)
    endif

    CALL flush(unit_)

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


  SUBROUTINE write_cplx_array_rank_1(A,title,UNIT,SHORT)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ WRITE COMPLEX ARRAY A(n1) $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    CHARACTER(LEN=*),           INTENT(IN)  :: title
    COMPLEX(8), DIMENSION(:), INTENT(IN)    :: A
    LOGICAL,          INTENT(IN), OPTIONAL  :: SHORT
    INTEGER,          INTENT(IN), OPTIONAL  :: UNIT
    INTEGER                                 :: i1,n1,unit_
    LOGICAL                                 :: short_
    CHARACTER(LEN=400)                      :: fmt_A

    unit_ = log_unit ! DEFAULT: STANDARD OUTPUT
    IF(PRESENT(UNIT)) unit_ = UNIT

    short_ = present(short)

    n1 = SIZE(A)
    if(n1==1) return

    CALL dump_message(UNIT=unit_,TEXT=title)

   if(n1>0)then
    IF(short_)THEN
      WRITE(fmt_A,*)  '(',n1,'(2(a,f10.6),a,x))'
    ELSE
      WRITE(fmt_A,*)  '(',n1,'(2(a,f20.16),a,x))'
    ENDIF
    WRITE(unit_,fmt_A) ('(',DBLE(A(i1)),',',AIMAG(A(i1)),')',i1=1,n1)
   endif

    CALL flush(unit_)
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


  SUBROUTINE write_real_array_rank_1(A,title,UNIT,SHORT)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ WRITE REAL ARRAY A(n1) $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    CHARACTER(LEN=*),        INTENT(IN) :: title
    REAL(8), DIMENSION(:), INTENT(IN)   :: A
    LOGICAL,       INTENT(IN), OPTIONAL :: SHORT
    INTEGER,       INTENT(IN), OPTIONAL :: UNIT
    INTEGER                             :: i1,n1,unit_
    LOGICAL                             :: short_
    CHARACTER(LEN=400)                  :: fmt_A

    unit_ = log_unit ! DEFAULT: STANDARD OUTPUT
    IF(PRESENT(UNIT)) unit_ = UNIT

    short_ = present(short)

    n1 = SIZE(A)

    if(n1==1) return

    CALL dump_message(UNIT=unit_,TEXT=title)

    IF(short_)THEN
      WRITE(fmt_A,*) '(',n1,'(f10.6,x))'
    ELSE
      WRITE(fmt_A,*) '(',n1,'(f20.16,x))'
    ENDIF
    WRITE(unit_,fmt_A) (A(i1),i1=1,n1)

    CALL flush(unit_)
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


  SUBROUTINE write_intg_array_rank_1(A,title,UNIT)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ WRITE INTEGER ARRAY A(n1) $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    CHARACTER(LEN=*),      INTENT(IN) :: title
    INTEGER, DIMENSION(:), INTENT(IN) :: A
    INTEGER,     INTENT(IN), OPTIONAL :: UNIT
    INTEGER                           :: i1,n1,unit_
    CHARACTER(LEN=400)                :: fmt_A

    unit_ = log_unit ! DEFAULT: STANDARD OUTPUT
    IF(PRESENT(UNIT)) unit_ = UNIT

    n1 = SIZE(A)
    if(n1==1) return

    CALL dump_message(UNIT=unit_,TEXT=title)

    WRITE(fmt_A,*) '(2x,',n1,'(I2,x))'
    WRITE(unit_,fmt_A) (A(i1),i1=1,n1)

    CALL flush(unit_)
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


  SUBROUTINE write_bool_array_rank_1(A,title,UNIT)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$$ WRITE BOOLEAN ARRAY A(n1) $$
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    CHARACTER(LEN=*),      INTENT(IN) :: title
    LOGICAL, DIMENSION(:), INTENT(IN) :: A
    INTEGER,     INTENT(IN), OPTIONAL :: UNIT
    INTEGER                           :: i1,n1,unit_
    CHARACTER(LEN=400)                :: fmt_A

    unit_ = log_unit ! STANDARD OUTPUT
    IF(PRESENT(UNIT)) unit_ = UNIT

    n1 = SIZE(A)
    if(n1==1) return

    CALL dump_message(UNIT=unit_,TEXT=title)

    WRITE(fmt_A,*) '(2x,',n1,'(L2,x))'
    WRITE(unit_,fmt_A) (A(i1),i1=1,n1)

    CALL flush(unit_)

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
! 
! 
!   logical function mat_3_3_rs(mat,detb,det2b,pdetb)
!   implicit none
!      integer               :: i,j,k,n
!      real(4),optional      :: detb
!      integer,optional      :: pdetb
!      complex(4),optional   :: det2b
!      real(4)               :: q
!      real(4),intent(inout) :: mat(:,:)
! 
!      n=size(mat,1)
!      mat_3_3_rs=n<4
!      if(n>=4) return
! 
!      SELECT CASE(n)
!        
!        CASE(3)
! 
!         mat=matrixinverse(mat,q)
!         if(present(detb))then
!          det2b=1.; pdetb=0; detb=q
!         endif
! 
!        CASE(1)
! 
!         q=mat(1,1)
!         if(abs(q)<epsilonr) q=epsilonr
!         mat=1.d0/q
!         if(present(detb))then
!          det2b=1.; pdetb=0; detb=q
!         endif
! 
!        CASE(2)
! 
!         q=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
!         if(abs(q)<epsilonr) q=epsilonr
!         call swap(mat(1,1),mat(2,2))
!         mat(1,2)=-mat(1,2)
!         mat(2,1)=-mat(2,1)
!         mat=mat/q
!         if(present(detb))then
!          det2b=1.; pdetb=0; detb=q
!         endif
! 
!     END SELECT
! 
!   end function
! 
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! 
!   logical function mat_3_3_r(mat,detb,det2b,pdetb)
!   implicit none
!      integer               :: i,j,k,n
!      real(8),optional      :: detb
!      integer,optional      :: pdetb
!      complex(8),optional   :: det2b
!      real(8)               :: q
!      real(8),intent(inout) :: mat(:,:)
! 
!      n=size(mat,1)
!      mat_3_3_r=n<4
!      if(n>=4) return
! 
!      SELECT CASE(n)
!        
!        CASE(3)
! 
!         mat=matrixinverse(mat,q)
!         if(present(detb))then
!          det2b=1.; pdetb=0; detb=q
!         endif
! 
!        CASE(1)
! 
!         q=mat(1,1)
!         if(abs(q)<epsilonr) q=epsilonr
!         mat=1.d0/q
!         if(present(detb))then
!          det2b=1.; pdetb=0; detb=q
!         endif
! 
!        CASE(2)
! 
!         q=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
!         if(abs(q)<epsilonr) q=epsilonr
!         call swap(mat(1,1),mat(2,2))
!         mat(1,2)=-mat(1,2)
!         mat(2,1)=-mat(2,1)
!         mat=mat/q
!         if(present(detb))then
!          det2b=1.; pdetb=0; detb=q
!         endif
! 
!     END SELECT
! 
!   end function
! 
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! 
!   logical function mat_3_3_cs(mat,detb,detb2,pdetb)
!   implicit none
!     integer                   :: n,i,j
!     integer,optional          :: pdetb
!     real(4),optional          :: detb
!     complex(4),intent(inout)  :: mat(:,:)
!     complex(4),optional       :: detb2
!     complex(4)                :: q
! 
!      n=size(mat,1)
!      mat_3_3_cs=n<4
!      if(n>=4) return
! 
!       SELECT CASE(n) 
! 
!        CASE(3)
! 
!         mat=matrixinverse(mat,q)
!         if(present(detb2))then
!          detb2=q; pdetb=0; detb=1.
!         endif
! 
!        CASE(1)
! 
!         q=mat(1,1)
!         if(abs(q)<epsilonr) q=epsilonr
!         mat=1.d0/q
!         if(present(detb2))then
!          detb2=q; pdetb=0; detb=1.
!         endif
! 
!        CASE(2)
! 
!         q=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
!         if(abs(q)<epsilonr) q=epsilonr
!         call swap(mat(1,1),mat(2,2))
!         mat(1,2)=-mat(1,2)
!         mat(2,1)=-mat(2,1)
!         mat=mat/q
!         if(present(detb2))then
!          detb2=q; pdetb=0; detb=1.
!         endif
! 
!        END SELECT
! 
!   end function
! 
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! 
!   logical function mat_3_3_c(mat,detb,detb2,pdetb)
!   implicit none
!     integer                   :: n,i,j
!     integer,optional          :: pdetb
!     real(8),optional          :: detb
!     complex(8),intent(inout)  :: mat(:,:)
!     complex(8),optional       :: detb2
!     complex(8)                :: q
! 
!      n=size(mat,1)
!      mat_3_3_c=n<4
!      if(n>=4) return
! 
!       SELECT CASE(n) 
! 
!        CASE(3)
! 
!         mat=matrixinverse(mat,q)
!         if(present(detb2))then
!          detb2=q; pdetb=0; detb=1.
!         endif
! 
!        CASE(1)
! 
!         q=mat(1,1)
!         if(abs(q)<epsilonr) q=epsilonr
!         mat=1.d0/q
!         if(present(detb2))then
!          detb2=q; pdetb=0; detb=1.
!         endif
! 
!        CASE(2)
! 
!         q=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
!         if(abs(q)<epsilonr) q=epsilonr
!         call swap(mat(1,1),mat(2,2))
!         mat(1,2)=-mat(1,2)
!         mat(2,1)=-mat(2,1)
!         mat=mat/q
!         if(present(detb2))then
!          detb2=q; pdetb=0; detb=1.
!         endif
! 
!        END SELECT
! 
!   end function
! 
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! 
!   logical function mat_3_3_q(mat,detb,det2b,pdetb)
!   implicit none
!     integer                :: i,j,k,n
!     real(16),optional      :: detb
!     integer,optional       :: pdetb
!     complex(16),optional   :: det2b
!     real(16)               :: q
!     real(16),intent(inout) :: mat(:,:)
! 
!      n=size(mat,1)
!      mat_3_3_q=n<4
!      if(n>=4) return
! 
!     SELECT CASE(n)
! 
!       CASE(3)
! 
!         mat=matrixinverse(mat,q)
!         if(present(detb))then
!          det2b=1.; pdetb=0; detb=q
!         endif
! 
!       CASE(1)
! 
!         q=mat(1,1)
!         if(abs(q)<epsilonr) q=epsilonr
!         mat=1.d0/q
!         if(present(detb))then
!          det2b=1.; pdetb=0; detb=q
!         endif
! 
!        CASE(2)
! 
!         q=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
!         if(abs(q)<epsilonr) q=epsilonr
!         call swap(mat(1,1),mat(2,2))
!         mat(1,2)=-mat(1,2)
!         mat(2,1)=-mat(2,1)
!         mat=mat/q
!         if(present(detb))then
!          det2b=1.; pdetb=0; detb=q
!         endif
! 
!       END SELECT
! 
!   end function
! 
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
! 
!   logical function mat_3_3_qc(mat,detb,det2b,pdetb)
!   implicit none
!   real(16),optional          :: detb
!   integer,optional           :: pdetb
!   complex(16),optional       :: det2b
!   integer                    :: n,i,j
!   complex(16),intent(inout)  :: mat(:,:)
!   complex(16)                :: q
! 
!     n=size(mat,1)
!     mat_3_3_qc=n<4
!     if(n>=4) return
! 
!      SELECT CASE(n)
! 
!        CASE(3)
! 
!         mat=matrixinverse(mat,q)
!         if(present(det2b))then
!          det2b=q; pdetb=0; detb=1.
!         endif
! 
!        CASE(1)
! 
!         q=mat(1,1)
!         if(abs(q)<epsilonr) q=epsilonr
!         mat=1.d0/q
!         if(present(det2b))then
!          det2b=q; pdetb=0; detb=1.
!         endif
! 
!        CASE(2)
! 
!         q=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
!         if(abs(q)<epsilonr) q=epsilonr
!         call swap(mat(1,1),mat(2,2))
!         mat(1,2)=-mat(1,2)
!         mat(2,1)=-mat(2,1)
!         mat=mat/q
!         if(present(det2b))then
!          det2b=q; pdetb=0; detb=1.
!         endif
! 
!       END SELECT
! 
!   end function 

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
! 
!          !------------------------!
!          !------------------------!
!          !------------------------!
!          !------------------------!
! 
!       subroutine test_arpack_vector
!       implicit none
!       integer,parameter    :: nnn=150
!       real(8)              :: A(nnn,nnn),B(nnn,nnn),C(nnn,nnn),vec(nnn),vec2(nnn),v(nnn,nnn)
!       integer              :: i,j,n,jjj,jj,nval,nconv
! 
!       A=0.d0
!       do i=1,nnn
!        A(i,i)=drand1()
!        do j=i+1,nnn
!          A(i,j)=drand1()
!          A(j,i)=A(i,j)
!        enddo
!       enddo
!       
!       call reset_timer(jjj)
!       write(*,*) 'start lapack calculations'
!       call eigenvector_matrix(lsize=nnn,mat=A,vaps=vec,eigenvec=B)
!       write(*,*) 'done'
!       call timer_fortran(jjj,'LAPACK TOOK : ', unit_=6)
! 
!       if(maxval(abs(A-transpose(A)))>1.d-5) stop 'not symm'
!       call reset_timer(jjj)
!       write(*,*) 'start ARPACK calculations'
! 
!       write(*,*) 'please enter nval'
!       read(*,*) nval
! 
!       write(*,*) '---------------------------------------'
!       call arpack_eigenvector_sym_matrix(.false.,'BE',1d-8,nnn,.true.,vec2(1:nval),v(1:nnn,1:nval),nval,nconv,mat_)
!       call timer_fortran(jjj,'ARPACK TOOK : ', unit_=6)
!       write(*,*) '---------------------------------------'
!       write(*,*) 'ARPACK EIGEN BE: ', vec2(1:nconv/2)
!       write(*,*) '---------------------------------------'
!       call arpack_eigenvector_sym_matrix(.false.,'SA',1d-8,nnn,.true.,vec2(1:nval),v(1:nnn,1:nval),nval,nconv,mat_)
!       write(*,*) '---------------------------------------'
!       write(*,*) 'ARPACK EIGEN SA: ', vec2(1:nconv)
!       write(*,*) '---------------------------------------'
!       call timer_fortran(jjj,'ARPACK TOOK : ', unit_=6)
!       write(*,*) '---------------------------------------'
!       write(*,*) 'LAPACK EIGEN : ', vec(1:nval)
!       write(*,*) '---------------------------------------'
!       stop 'done'
! 
!       contains
! 
!        subroutine mat_(n,w,v)
!         integer n
!         Double precision,intent(in)    ::  v(n) 
!         Double precision,intent(inout) ::  w(n)
!          w=MATMUL(A,v)
!        return
!        end subroutine
! 
!       end subroutine
! 
!          !------------------------!
!          !------------------------!
!          !------------------------!
!          !------------------------!
! 
!       subroutine test_arpack_vector_
!       implicit none
!       integer,parameter    :: nnn=150
!       complex(8)           :: A(nnn,nnn),B(nnn,nnn),C(nnn,nnn),v(nnn,nnn)
!       real(8)              :: vec2(nnn),vec(nnn)
!       integer              :: i,j,n,jjj,jj,nval,nconv
! 
!       A=0.d0
!       do i=1,nnn
!        A(i,i)=drand1()
!        do j=i+1,nnn
!          A(i,j)=drand1()+imi*drand1()
!          A(j,i)=conjg(A(i,j))
!        enddo
!       enddo
! 
!       call reset_timer(jjj)
!       write(*,*) 'start lapack calculations'
!       call eigenvector_matrix(lsize=nnn,mat=A,vaps=vec,eigenvec=B)
!       write(*,*) 'done'
!       call timer_fortran(jjj,'LAPACK TOOK : ', unit_=6)
! 
!       call reset_timer(jjj)
!       write(*,*) 'start ARPACK calculations'
!       write(*,*) 'please enter nval'
!       read(*,*) nval
! 
!       write(*,*) '---------------------------------------'
!       call arpack_eigenvector_sym_matrix_(.false.,'SR',1d-8,nnn,.true.,vec2(1:nval),v(1:nnn,1:nval),nval,nconv,mat_)
!       call timer_fortran(jjj,'ARPACK TOOK : ', unit_=6)
!       write(*,*) '---------------------------------------'
!       write(*,*) 'ARPACK EIGEN SR: ', vec2(1:nconv)
!       write(*,*) '---------------------------------------'
!       call arpack_eigenvector_sym_matrix_(.false.,'SM',1d-8,nnn,.true.,vec2(1:nval),v(1:nnn,1:nval),nval,nconv,mat_)
!       write(*,*) '---------------------------------------'
!       write(*,*) 'ARPACK EIGEN SM: ', vec2(1:nconv)
!       write(*,*) '---------------------------------------'
!       call timer_fortran(jjj,'ARPACK TOOK : ', unit_=6)
!       write(*,*) '---------------------------------------'
!       write(*,*) 'LAPACK EIGEN : ', vec(1:nval)
!       write(*,*) '---------------------------------------'
!       stop 'done'
! 
!       contains
! 
!        subroutine mat_(n,w,v)
!         integer n
!         complex(8),intent(in)   ::  v(n)  
!         complex(8),intent(inout)::  w(n)
!          w=MATMUL(A,v)
!        return
!        end subroutine
! 
!       end subroutine
! 
!          !------------------------!
!          !------------------------!
!          !------------------------!
!          !------------------------!
!  
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       subroutine arpack_eigenvector_sym_matrix(verbose,amode,tol,maxn,rvec,values,v,maxncv,nconv_,av,sigma_,mode_)
!       implicit none
! 
!       !---------------------------------------------!
!       !     Solve A*x = lambda*x in regular mode    !
!       !---------------------------------------------!
! 
!       integer              :: maxn, maxnev, maxncv,nconv_
!       Double precision     :: values(maxncv),v(maxn,maxncv), workl(maxncv*(maxncv+8))
!       Double precision     :: workd(3*maxn), d(maxncv,2), resid(maxn), ax(maxn)
!       logical              :: select(maxncv)
!       integer              :: iparam(11), ipntr(11)
!       character            :: bmat*1, which*2
!       integer              :: ido, n, nev, ncv, lworkl, info, ierr, j, nconv, maxitr, mode, ishfts
!       logical              :: rvec
!       Double precision     :: tol, sigma,dnrm2
!       external             :: dnrm2
!       Double precision     :: zero
!       parameter        (zero = 0.0D+0)
!       character(2)         :: amode
!       logical              :: verbose
!       real(8),optional     :: sigma_
!       integer,optional     :: mode_
!      !----------------------------------!
!      ! amode : SM around 0              !
!      ! amode : BE extremas eigenvalues  !
!      !----------------------------------!
!       interface
!        subroutine av(n,w,v)
!         integer          :: n
!         Double precision,intent(in)    :: v(n)
!         Double precision,intent(inout) :: w(n)
!        end subroutine
!       end interface
! 
! #ifdef _ARPACK
! 
!       if(present(sigma_))then
!         sigma=sigma_
!       else
!         sigma=0.d0
!       endif
! 
!       maxnev=maxncv-1
!       n=maxn
!       nev=maxnev
!       ncv=maxncv
! 
! !     %----------------------------------------------------%
! !     | A standard eigenvalue                              |
! !     | problem is solved (BMAT = 'I'). NEV is the number  |
! !     | of eigenvalues to be approximated.  The user can   |
! !     | modify NEV, NCV, WHICH to solve problems of        |
! !     | different sizes, and to get different parts of the |
! !     | spectrum.  However, The following conditions must  |
! !     | be satisfied:                                      |
! !     |                   N <= MAXN,                       | 
! !     |                 NEV <= MAXNEV,                     |
! !     |             NEV + 1 <= NCV <= MAXNCV               | 
! !     %----------------------------------------------------% 
!       if ( n .gt. maxn ) then
!          write(*,*) 'n,maxn : ',n,maxn       
!          print *, ' ERROR with _SDRV1: N is greater than MAXN '
!          stop 'arpack error'
!       else if ( nev .gt. maxnev ) then
!          print *, ' ERROR with _SDRV1: NEV is greater than MAXNEV '
!          stop 'arpack error'
!       else if ( ncv .gt. maxncv ) then
!          print *, ' ERROR with _SDRV1: NCV is greater than MAXNCV '
!          stop 'arpack error'
!       end if
!       bmat = 'I'
!       which = amode
! !     %--------------------------------------------------%
! !     | The work array WORKL is used in DSAUPD as        |
! !     | workspace.  Its dimension LWORKL is set as       |
! !     | illustrated below.  The parameter TOL determines |
! !     | the stopping criterion.  If TOL<=0, machine      |
! !     | precision is used.  The variable IDO is used for |
! !     | reverse communication and is initially set to 0. |
! !     | Setting INFO=0 indicates that a random vector is |
! !     | generated in DSAUPD to start the Arnoldi         |
! !     | iteration.                                       |
! !     %--------------------------------------------------%
!       lworkl = ncv*(ncv+8)
!       info = 0
!       ido = 0
! !     %---------------------------------------------------%
! !     | This program uses exact shifts with respect to    |
! !     | the current Hessenberg matrix (IPARAM(1) = 1).    |
! !     | IPARAM(3) specifies the maximum number of Arnoldi |
! !     | iterations allowed.  Mode 1 of DSAUPD is used     |
! !     | (IPARAM(7) = 1).  All these options may be        |
! !     | changed by the user. For details, see the         |
! !     | documentation in DSAUPD.                          |
! !     %---------------------------------------------------%
!       ishfts = 1
!       maxitr = 300
!       mode   = 1
! 
!       if(present(mode_))then
!        mode=mode_
!        if(mode/=1)then
!           write(*,*) 'ARPACK routine needs to be modified for mode/=1';stop
!        endif
!       else
!        mode=1
!       endif
! 
!       iparam(1) = ishfts 
!       iparam(3) = maxitr 
!       iparam(7) = mode 
! 
!  10   continue
! !        %---------------------------------------------%
! !        | Repeatedly call the routine DSAUPD and take | 
! !        | actions indicated by parameter IDO until    |
! !        | either convergence is indicated or maxitr   |
! !        | has been exceeded.                          |
! !        %---------------------------------------------%
!          call dsaupd ( ido, bmat, n, which, nev, tol, resid, ncv, v, maxn, iparam, ipntr, workd, workl, lworkl, info )
!          if (ido .eq. -1 .or. ido .eq. 1) then
! !           %--------------------------------------%
! !           | Perform matrix vector multiplication |
! !           |              y <--- OP*x             |
! !           %--------------------------------------%
!             call av (n, workd(ipntr(2)), workd(ipntr(1)))
!             go to 10
!          end if 
!       if ( info .lt. 0 ) then
!          print *, ' '
!          print *, ' Error with _saupd, info = ', info
!          print *, ' Check documentation in _saupd '
!          print *, ' '
!       else 
! !        %-------------------------------------------%
! !        | No fatal errors occurred.                 |
! !        | Post-Process using DSEUPD.                |
! !        | Computed eigenvalues may be extracted.    |  
! !        | Eigenvectors may also be computed now if  |
! !        | desired.  (indicated by rvec = .true.)    | 
! !        %-------------------------------------------%
!          call dseupd ( rvec, 'All', select, d, v, maxn, sigma, bmat, n, which, &
!              & nev, tol, resid, ncv, v, maxn, iparam, ipntr, workd, workl, lworkl, ierr )
! !        %----------------------------------------------%
! !        | Eigenvalues are returned in the first column |
! !        | of the two dimensional array D and the       |
! !        | corresponding eigenvectors are returned in   |
! !        | the first NEV columns of the two dimensional |
! !        | array V if requested.  Otherwise, an         |
! !        | orthogonal basis for the invariant subspace  |
! !        | corresponding to the eigenvalues in D is     |
! !        | returned in V.                               |
! !        %----------------------------------------------%
!          if ( ierr .ne. 0) then
!              print *, ' '
!              print *, ' Error with _seupd, info = ', ierr
!              print *, ' Check the documentation of _seupd. '
!              print *, ' '
!          else
!              nconv =  iparam(5)
!              do 20 j=1, nconv
! !               %---------------------------%
! !               | Compute the residual norm |
! !               |   ||  A*x - lambda*x ||   |
! !               | for the NCONV accurately  |
! !               | computed eigenvalues and  |
! !               | eigenvectors.  (iparam(5) |
! !               | indicates how many are    |
! !               | accurate to the requested |
! !               | tolerance)                |
! !               %---------------------------%
!                 call av(n, ax, v(1,j))
!                 call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
!                 d(j,2) = dnrm2(n, ax, 1)
!                 d(j,2) = d(j,2) / abs(d(j,1))
!  20          continue
!             if(verbose)  call dmout(6, nconv, 2, d, maxncv, -6, 'Ritz values and relative residuals')
!          end if
! !        %------------------------------------------%
! !        | Print additional convergence information |
! !        %------------------------------------------%
!          if ( info .eq. 1) then
!            if(verbose)then
!             print *, ' '
!             print *, ' Maximum number of iterations reached.'
!             print *, ' '
!            endif
!          else if ( info .eq. 3) then
!            if(verbose)then
!             print *, ' ' 
!             print *, ' No shifts could be applied during implicit', ' Arnoldi update, try increasing NCV.'
!             print *, ' '
!            endif
!          end if      
!         if(verbose)then
!          print *, ' '
!          print *, ' _SDRV1 '
!          print *, ' ====== '
!          print *, ' '
!          print *, ' Size of the matrix is ', n
!          print *, ' The number of Ritz values requested is ', nev
!          print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
!          print *, ' What portion of the spectrum: ', which
!          print *, ' The number of converged Ritz values is ',   nconv 
!          print *, ' The number of Implicit Arnoldi update',  ' iterations taken is ', iparam(3)
!          print *, ' The number of OP*x is ', iparam(9)
!          print *, ' The convergence criterion is ', tol
!          print *, ' '
!         endif
!       end if
!       nconv_=nconv
!       values=d(:,1)
! 
!  9000 continue
! 
! #endif
! 
!  end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       subroutine arpack_eigenvector_sym_matrix_(verbose,amode,tol,maxn,rvec,values,v,maxncv,nconv_,av,sigma_,mode_)
!       implicit none
!       !---------------------------------------------!
!       !     Solve A*x = lambda*x in regular mode    !
!       !---------------------------------------------!
!       integer           maxn, maxnev, maxncv, nconv_
!       integer           iparam(11), ipntr(14)
!       logical           select(maxncv)
!       Complex(8)        ax(maxn),d(maxncv),v(maxn,maxncv),workd(3*maxn), workev(3*maxncv), resid(maxn), workl(3*maxncv*maxncv+5*maxncv)
!       Double precision  rwork(maxncv), rd(maxncv,3),values(:)
!       character         bmat*1, which*2
!       integer           ido, n, nev, ncv, lworkl, info, j, ierr, nconv, maxitr, ishfts, mode
!       Complex(8)        sigma
!       Double precision  tol
!       logical           rvec
!       Double precision  dznrm2 , dlapy2 
!       external          dznrm2 , dlapy2  
!       character(2)      amode
!       logical           verbose
!      !----------------------------------!
!      ! amode : SM around 0              !
!      ! amode : BE extremas eigenvalues  !
!      !----------------------------------!
! 
!       interface
!        subroutine av(n,w,v)
!         integer                  :: n
!         complex(8),intent(in)    :: v(n)
!         complex(8),intent(inout) :: w(n)
!        end subroutine
!       end interface
! 
!       complex(8),optional :: sigma_
!       integer,optional    :: mode_
! 
! 
! #ifdef _ARPACK
! 
!       if(present(sigma_))then
!         sigma=sigma_
!       else
!         sigma=0.d0
!       endif
! 
!       maxnev=maxncv-1
!       n=maxn
!       nev=maxnev
!       ncv=maxncv
! 
!       if ( n .gt. maxn ) then
!          print *, ' ERROR with _NDRV1: N is greater than MAXN '
!          go to 9000
!       else if ( nev .gt. maxnev ) then
!          print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
!          go to 9000
!       else if ( ncv .gt. maxncv ) then
!          print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
!          go to 9000
!       end if
! 
!       bmat  = 'I'
!       which = amode
!       lworkl  = 3*ncv**2+5*ncv 
!       ido    = 0
!       info   = 0
!       ishfts = 1
!       maxitr = 300
!       mode   = 1
! 
!       if(present(mode_))then
!        mode=mode_
!        if(mode/=1)then
!           write(*,*) 'ARPACK routine needs to be modified for mode/=1';stop
!        endif
!       else
!        mode=1
!       endif
! 
!       iparam(1) = ishfts
!       iparam(3) = maxitr 
!       iparam(7) = mode 
! 
!  10   continue
! 
!       call znaupd  ( ido, bmat, n, which, nev, tol, resid, ncv, v, maxn, iparam, ipntr, workd, workl, lworkl, rwork,info )
!       if (ido .eq. -1 .or. ido .eq. 1) then
!          call av (n, workd(ipntr(2)), workd(ipntr(1)))
!          go to 10
!       end if
!       if ( info .lt. 0 ) then
!          print *, ' '
!          print *, ' Error with _naupd, info = ', info
!          print *, ' Check the documentation of _naupd'
!          print *, ' '
!       else 
!          call zneupd  (rvec, 'A', select, d, v, maxn, sigma, workev, bmat, n, which, &
!             & nev, tol, resid, ncv, v, maxn, iparam, ipntr, workd, workl, lworkl, rwork, ierr)
!          if ( ierr .ne. 0) then
!              print *, ' '
!              print *, ' Error with _neupd, info = ', ierr
!              print *, ' Check the documentation of _neupd. '
!              print *, ' '
!          else
!              nconv = iparam(5)
!              do 20 j=1, nconv
!                 call av(n, ax,v(1,j))
!                 call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
!                 rd(j,1) = dble (d(j))
!                 rd(j,2) = dimag (d(j))
!                 rd(j,3) = dznrm2 (n, ax, 1)
!                 rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
!  20          continue
!              if(verbose) call dmout (6, nconv, 3, rd, maxncv, -6, 'Ritz values (Real, Imag) and relative residuals')
!           end if
!         if(verbose)then
!          if ( info .eq. 1) then
!              print *, ' '
!              print *, ' Maximum number of iterations reached.'
!              print *, ' '
!          else if ( info .eq. 3) then
!              print *, ' ' 
!              print *, ' No shifts could be applied during implicit', ' Arnoldi update, try increasing NCV.'
!              print *, ' '
!          end if      
!          print *, ' '
!          print *, '_NDRV1'
!          print *, '====== '
!          print *, ' '
!          print *, ' Size of the matrix is ', n
!          print *, ' The number of Ritz values requested is ', nev
!          print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
!          print *, ' What portion of the spectrum: ', which
!          print *, ' The number of converged Ritz values is ', nconv 
!          print *, ' The number of Implicit Arnoldi update', ' iterations taken is ', iparam(3)
!          print *, ' The number of OP*x is ', iparam(9)
!          print *, ' The convergence criterion is ', tol
!          print *, ' '
!         endif
!       end if
!       nconv_=nconv
!       values=real(d(:))
! 
!  9000 continue
! 
! #endif
!       end subroutine
! 
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! !**************************************************************************
! 
!       Subroutine GramSchmidt (vect, dimreal, dim, numvects) 
!       IMPLICIT none 
!       INTEGER  :: dim, dimreal, numvects, iii, j 
!       REAL(8)  :: vect (dimreal, * ) 
!       REAL(8)  :: anor, aux (numvects) , temp
! 
!       DO iii = 1, numvects 
!        temp = norme(vect (1:dim,iii))
!        if(temp==0.d0)then
!          write(*,*) 'GramSchmidt error, 0 division'
!          write(*,*) 'temp        : ', temp
!          write(*,*) 'dim,dimreal : ', dim,dimreal
!          write(*,*) 'vect        : ', vect(1:dim,iii)
!          stop 'critical'
!        endif
!        anor = 1.d0 / temp
!        vect(1:dim,iii)=vect(1:dim,iii)*anor
!       enddo 
!                                                                         
!       DO iii = 2, numvects 
!        DO j = 1, iii - 1 
!         aux(j) = - DOT_PRODUCT(vect(1:dim,iii),vect(1:dim,j)) 
!        enddo 
!        DO j = 1, iii - 1 
!         vect(1:dim,iii)=vect(1:dim,iii)+aux(j)*vect(1:dim,j)
!        enddo 
!        temp = norme(vect(1:dim,iii))
!        if(temp==0.d0)then
!          write(*,*) 'number of vectors : ', numvects
!          write(*,*) 'GramSchmidt error, 0 division bis'
!          stop 'critical'
!        endif
!        anor = 1.d0 / temp
!        vect(1:dim,iii)=vect(1:dim,iii)*anor
!       enddo 
! 
!       RETURN 
!       END subroutine
! 
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
!                                                                         
!       SUBROUTINE EHOBKS (A, N, M1, M2, Z, IZ) 
!       DIMENSION A ( * ), Z (IZ, * ) 
!       real(8) A, Z, H, S 
!       IF (N.EQ.1) GOTO 30 
!       DO 25 I = 2, N 
!          L = I - 1 
!          IA = (I * L) / 2 
!          H = A (IA + I) 
!          IF (H.EQ.0.D0) GOTO 25 
! !                                  DERIVES EIGENVECTORS M1 TO M2 OF     
! !                                  THE ORIGINAL MATRIX FROM EIGENVECTORS
! !                                  M1 TO M2 OF THE SYMMETRIC            
! !                                  TRIDIAGONAL MATRIX                   
!          DO 20 J = M1, M2 
!             S = 0.0D0 
!             DO 10 K = 1, L 
!                S = S + A (IA + K) * Z (K, J) 
!    10       END DO 
!             S = S / H 
!             DO 15 K = 1, L 
!                Z (K, J) = Z (K, J) - S * A (IA + K) 
!    15       END DO 
!    20    END DO 
!    25 END DO 
!    30 RETURN 
!       END SUBROUTINE
!                                                                       
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
!                                                                        
!       SUBROUTINE EHOUSS (A, N, D, E, E2) 
!       DIMENSION A ( * ), D (N), E (N), E2 (N) 
!       real(8) A, D, E, E2, ZERO, H, SCALE, F, G, HH 
!       DATA ZERO / 0.0D0 / 
! !                                  FIRST EXECUTABLE STATEMENT           
!       NP1 = N + 1 
!       NN = (N * NP1) / 2 - 1 
!       NBEG = NN + 1 - N 
!       DO 70 II = 1, N 
!          I = NP1 - II 
!          L = I - 1 
!          H = ZERO 
!          SCALE = ZERO 
!          IF (L.LT.1) GOTO 10 
! !                                  SCALE ROW (ALGOL TOL THEN NOT NEEDED)
!          NK = NN 
!          DO 5 K = 1, L 
!             SCALE = SCALE+DABS (A (NK) ) 
!             NK = NK - 1 
!     5    END DO 
!          IF (SCALE.NE.ZERO) GOTO 15 
!    10    E (I) = ZERO 
!          E2 (I) = ZERO 
!          GOTO 65 
!    15    NK = NN 
!          DO 20 K = 1, L 
!             A (NK) = A (NK) / SCALE 
!             H = H + A (NK) * A (NK) 
!             NK = NK - 1 
!    20    END DO 
!          E2 (I) = SCALE * SCALE * H 
!          F = A (NN) 
!          G = - DSIGN (DSQRT (H), F) 
!          E (I) = SCALE * G 
!          H = H - F * G 
!          A (NN) = F - G 
!          IF (L.EQ.1) GOTO 55 
!          F = ZERO 
!          JK1 = 1 
!          DO 40 J = 1, L 
!             G = ZERO 
!             IK = NBEG + 1 
!             JK = JK1 
! !                                  FORM ELEMENT OF A*U                  
!             DO 25 K = 1, J 
!                G = G + A (JK) * A (IK) 
!                JK = JK + 1 
!                IK = IK + 1 
!    25       END DO 
!             JP1 = J + 1 
!             IF (L.LT.JP1) GOTO 35 
!             JK = JK + J - 1 
!             DO 30 K = JP1, L 
!                G = G + A (JK) * A (IK) 
!                JK = JK + K 
!                IK = IK + 1 
!    30       END DO 
! !                                  FORM ELEMENT OF P                    
!    35       E (J) = G / H 
!             F = F + E (J) * A (NBEG + J) 
!             JK1 = JK1 + J 
!    40    END DO 
!          HH = F / (H + H) 
! !                                  FORM REDUCED A                       
!          JK = 1 
!          DO 50 J = 1, L 
!             F = A (NBEG + J) 
!             G = E (J) - HH * F 
!             E (J) = G 
!             DO 45 K = 1, J 
!                A (JK) = A (JK) - F * E (K) - G * A (NBEG + K) 
!                JK = JK + 1 
!    45       END DO 
!    50    END DO 
!    55    DO 60 K = 1, L 
!             A (NBEG + K) = SCALE * A (NBEG + K) 
!    60    END DO 
!    65    D (I) = A (NBEG + I) 
!          A (NBEG + I) = H * SCALE * SCALE 
!          NBEG = NBEG - I + 1 
!          NN = NN - I 
!    70 END DO 
!       RETURN 
!       END SUBROUTINE
!      
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
!                                                                   
!       SUBROUTINE EIGRS (A, N, JOBN, D, Z, IZ, WK, IER) 
!       IMPLICIT none 
!       INTEGER N, JOBN, IZ, IER 
!       real(8) A ( * ), D ( * ), WK ( * ), Z (IZ, * ) 
!       INTEGER IJOB, IR, JR, IJ, JI, NP1 
!       INTEGER JER, NA, ND, IIZ, IBEG, IL, KK, LK, I, J, K, L 
!       real(8) ANORM, ASUM, PI, SUMZ, SUMR, AN, S, TEN, RDELP,   &
!       ZERO, ONE, THOUS                                                  
!       DATA RDELP / 0.222045D-15 / 
!       DATA ZERO, ONE / 0.0D0, 1.0D0 /, TEN / 10.0D0 /, THOUS / 1000.0D0 &
!       /                                                                 
! !                                  INITIALIZE ERROR PARAMETERS          
! !                                  FIRST EXECUTABLE STATEMENT           
!       IER = 0 
!       JER = 0 
!       IF (JOBN.LT.10) GOTO 15 
! !                                  CONVERT TO SYMMETRIC STORAGE MODE    
!       K = 1 
!       JI = N - 1 
!       IJ = 1 
!       DO 10 J = 1, N 
!          DO 5 I = 1, J 
!             A (K) = A (IJ) 
!             IJ = IJ + 1 
!             K = K + 1 
!     5    END DO 
!          IJ = IJ + JI 
!          JI = JI - 1 
!    10 END DO 
!    15 IJOB = MOD (JOBN, 10) 
!       IF (IJOB.GE.0.AND.IJOB.LE.3) GOTO 20 
! !                                  WARNING ERROR - IJOB IS NOT IN THE   
! !                                    RANGE                              
!       IER = 66 
!       IJOB = 1 
!       GOTO 25 
!    20 IF (IJOB.EQ.0) GOTO 35 
!    25 IF (IZ.GE.N) GOTO 30 
! !                                  WARNING ERROR - IZ IS LESS THAN N    
! !                                    EIGENVECTORS CAN NOT BE COMPUTED,  
! !                                    IJOB SET TO ZERO                   
!       IER = 67 
!       IJOB = 0 
!    30 IF (IJOB.EQ.3) GOTO 75 
!    35 NA = (N * (N + 1) ) / 2 
!       IF (IJOB.NE.2) GOTO 45 
!       DO 40 I = 1, NA 
!          WK (I) = A (I) 
!    40 END DO 
! !                                  SAVE INPUT A IF IJOB = 2             
!    45 ND = 1 
!       IF (IJOB.EQ.2) ND = NA + 1 
! !                                  REDUCE A TO SYMMETRIC TRIDIAGONAL    
! !                                    FORM                               
!       CALL EHOUSS (A, N, D, WK (ND), WK (ND) ) 
!       IIZ = 1 
!       IF (IJOB.EQ.0) GOTO 60 
!       IIZ = IZ 
! !                                  SET Z TO THE IDENTITY MATRIX         
!       DO 55 I = 1, N 
!          DO 50 J = 1, N 
!             Z (I, J) = ZERO 
!    50    END DO 
!          Z (I, I) = ONE 
!    55 END DO 
! !                                  COMPUTE EIGENVALUES AND EIGENVECTORS 
!    60 CALL EQRT2S (D, WK (ND), N, Z, IIZ, JER) 
!       IF (IJOB.EQ.0) GOTO 9000 
!       IF (JER.GT.128) GOTO 65 
! !                                  BACK TRANSFORM EIGENVECTORS          
!       CALL EHOBKS (A, N, 1, N, Z, IZ) 
!    65 IF (IJOB.LE.1) GOTO 9000 
! !                                  MOVE INPUT MATRIX BACK TO A          
!       DO 70 I = 1, NA 
!          A (I) = WK (I) 
!    70 END DO 
!       WK (1) = THOUS 
!       IF (JER.NE.0) GOTO 9000 
! !                                  COMPUTE 1 - NORM OF A                
!    75 ANORM = ZERO 
!       IBEG = 1 
!       DO 85 I = 1, N 
!          ASUM = ZERO 
!          IL = IBEG 
!          KK = 1 
!          DO 80 L = 1, N 
!             ASUM = ASUM + DABS (A (IL) ) 
!             IF (L.GE.I) KK = L 
!             IL = IL + KK 
!    80    END DO 
!          ANORM = DMAX1 (ANORM, ASUM) 
!          IBEG = IBEG + I 
!    85 END DO 
!       IF (ANORM.EQ.ZERO) ANORM = ONE 
! !                                  COMPUTE PERFORMANCE INDEX            
!       PI = ZERO 
!       DO 100 I = 1, N 
!          IBEG = 1 
!          S = ZERO 
!          SUMZ = ZERO 
!          DO 95 L = 1, N 
!             LK = IBEG 
!             KK = 1 
!             SUMZ = SUMZ + DABS (Z (L, I) ) 
!             SUMR = - D (I) * Z (L, I) 
!             DO 90 K = 1, N 
!                SUMR = SUMR + A (LK) * Z (K, I) 
!                IF (K.GE.L) KK = K 
!                LK = LK + KK 
!    90       END DO 
!             S = S + DABS (SUMR) 
!             IBEG = IBEG + L 
!    95    END DO 
!          IF (SUMZ.EQ.ZERO) GOTO 100 
!          PI = DMAX1 (PI, S / SUMZ) 
!   100 END DO 
!       AN = N 
!       PI = PI / (ANORM * TEN * AN * RDELP) 
!       WK (1) = PI 
!       IF (JOBN.LT.10) GOTO 9000 
! !                                  CONVERT BACK TO FULL STORAGE MODE    
!       NP1 = N + 1 
!       IJ = (N - 1) * NP1 + 2 
!       K = (N * (NP1) ) / 2 
!       DO 110 JR = 1, N 
!          J = NP1 - JR 
!          DO 105 IR = 1, J 
!             IJ = IJ - 1 
!             A (IJ) = A (K) 
!             K = K - 1 
!   105    END DO 
!          IJ = IJ - JR 
!   110 END DO 
!       JI = 0 
!       K = N - 1 
!       DO 120 I = 1, N 
!          IJ = I - N 
!          DO 115 J = 1, I 
!             IJ = IJ + N 
!             JI = JI + 1 
!             A (IJ) = A (JI) 
!   115    END DO 
!          JI = JI + K 
!          K = K - 1 
!   120 END DO 
!  9000 CONTINUE 
!       IF (IER.NE.0) CALL UERTST (IER, 'EIGRS ') 
!       IF (JER.EQ.0) GOTO 9005 
!       IER = JER 
!       CALL UERTST (IER, 'EIGRS ') 
!  9005 RETURN 
!       END SUBROUTINE
!      
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
!                                                                   
!       SUBROUTINE EQRT2S (D, E, N, Z, IZ, IER) 
!       DIMENSION D ( * ), E ( * ), Z (IZ, * ) 
!       real(8) D, E, Z, B, C, F, G, H, P, R, S, RDELP, ONE, ZERO 
!       DATA RDELP / 0.222045D-15 / 
!       DATA ZERO, ONE / 0.0D0, 1.0D0 / 
! !                                  MOVE THE LAST N-1 ELEMENTS           
! !                                  OF E INTO THE FIRST N-1 LOCATIONS    
! !                                  FIRST EXECUTABLE STATEMENT           
!       IER = 0 
!       K = 0 
!       IF (N.EQ.1) GOTO 9005 
!       DO 5 I = 2, N 
!          E (I - 1) = E (I) 
!     5 END DO 
!       E (N) = ZERO 
!       B = ZERO 
!       F = ZERO 
!       DO 60 L = 1, N 
!          J = 0 
!          H = RDELP * (DABS (D (L) ) + DABS (E (L) ) ) 
!          IF (B.LT.H) B = H 
! !                                  LOOK FOR SMALL SUB-DIAGONAL ELEMENT  
!          DO 10 M = L, N 
!             K = M 
!             IF (DABS (E (K) ) .LE.B) GOTO 15 
!    10    END DO 
!    15    M = K 
!          IF (M.EQ.L) GOTO 55 
!    20    IF (J.EQ.30) GOTO 85 
!          J = J + 1 
!          L1 = L + 1 
!          G = D (L) 
!          P = (D (L1) - G) / (E (L) + E (L) ) 
!          R = DABS (P) 
!          IF (RDELP * DABS (P) .LT.1.0D0) R = DSQRT (P * P + ONE) 
!          D (L) = E (L) / (P + DSIGN (R, P) ) 
!          H = G - D (L) 
!          DO 25 I = L1, N 
!             D (I) = D (I) - H 
!    25    END DO 
!          F = F + H 
! !                                  QL TRANSFORMATION                    
!          P = D (M) 
!          C = ONE 
!          S = ZERO 
!          MM1 = M - 1 
!          MM1PL = MM1 + L 
!          IF (L.GT.MM1) GOTO 50 
!          DO 45 II = L, MM1 
!             I = MM1PL - II 
!             G = C * E (I) 
!             H = C * P 
!             IF (DABS (P) .LT.DABS (E (I) ) ) GOTO 30 
!             C = E (I) / P 
!             R = DSQRT (C * C + ONE) 
!             E (I + 1) = S * P * R 
!             S = C / R 
!             C = ONE / R 
!             GOTO 35 
!    30       C = P / E (I) 
!             R = DSQRT (C * C + ONE) 
!             E (I + 1) = S * E (I) * R 
!             S = ONE / R 
!             C = C * S 
!    35       P = C * D (I) - S * G 
!             D (I + 1) = H + S * (C * G + S * D (I) ) 
!             IF (IZ.LT.N) GOTO 45 
! !                                  FORM VECTOR                          
!             DO 40 K = 1, N 
!                H = Z (K, I + 1) 
!                Z (K, I + 1) = S * Z (K, I) + C * H 
!                Z (K, I) = C * Z (K, I) - S * H 
!    40       END DO 
!    45    END DO 
!    50    E (L) = S * P 
!          D (L) = C * P 
!          IF (DABS (E (L) ) .GT.B) GOTO 20 
!    55    D (L) = D (L) + F 
!    60 END DO 
! !                                  ORDER EIGENVALUES AND EIGENVECTORS   
!       DO 80 I = 1, N 
!          K = I 
!          P = D (I) 
!          IP1 = I + 1 
!          IF (IP1.GT.N) GOTO 70 
!          DO 65 J = IP1, N 
!             IF (D (J) .GE.P) GOTO 65 
!             K = J 
!             P = D (J) 
!    65    END DO 
!    70    IF (K.EQ.I) GOTO 80 
!          D (K) = D (I) 
!          D (I) = P 
!          IF (IZ.LT.N) GOTO 80 
!          DO 75 J = 1, N 
!             P = Z (J, I) 
!             Z (J, I) = Z (J, K) 
!             Z (J, K) = P 
!    75    END DO 
!    80 END DO 
!       GOTO 9005 
!    85 IER = 128 + L 
!  9000 CONTINUE 
!       CALL UERTST (IER, 'EQRT2S') 
!  9005 RETURN 
!       END SUBROUTINE
! 
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
!  
!       SUBROUTINE UERTST (IER, NAME) 
!       INTEGER IER 
!       CHARACTER NAME * ( * ) 
!       INTEGER I, IEQDF, IOUNIT, LEVEL, LEVOLD, NIN, NMTB 
!       CHARACTER IEQ, NAMEQ (6), NAMSET (6), NAMUPK (6) 
!       DATA NAMSET / 'U', 'E', 'R', 'S', 'E', 'T' / 
!       DATA NAMEQ / 6 * ' ' / 
!       DATA LEVEL / 4 / , IEQDF / 0 / , IEQ / '=' / 
! !                                  UNPACK NAME INTO NAMUPK              
! !                                  FIRST EXECUTABLE STATEMENT           
!       CALL USPKD (NAME, 6, NAMUPK, NMTB) 
! !                                  GET OUTPUT UNIT NUMBER               
!       CALL UGETIO (1, NIN, IOUNIT) 
! !                                  CHECK IER                            
!       IF (IER.GT.999) GOTO 25 
!       IF (IER.LT. - 32) GOTO 55 
!       IF (IER.LE.128) GOTO 5 
!       IF (LEVEL.LT.1) GOTO 30 
! !                                  PRINT TERMINAL MESSAGE               
!       IF (IEQDF.EQ.1) WRITE (IOUNIT, 35) IER, NAMEQ, IEQ, NAMUPK 
!       IF (IEQDF.EQ.0) WRITE (IOUNIT, 35) IER, NAMUPK 
!       GOTO 30 
!     5 IF (IER.LE.64) GOTO 10 
!       IF (LEVEL.LT.2) GOTO 30 
! !                                  PRINT WARNING WITH FIX MESSAGE       
!       IF (IEQDF.EQ.1) WRITE (IOUNIT, 40) IER, NAMEQ, IEQ, NAMUPK 
!       IF (IEQDF.EQ.0) WRITE (IOUNIT, 40) IER, NAMUPK 
!       GOTO 30 
!    10 IF (IER.LE.32) GOTO 15 
! !                                  PRINT WARNING MESSAGE                
!       IF (LEVEL.LT.3) GOTO 30 
!       IF (IEQDF.EQ.1) WRITE (IOUNIT, 45) IER, NAMEQ, IEQ, NAMUPK 
!       IF (IEQDF.EQ.0) WRITE (IOUNIT, 45) IER, NAMUPK 
!       GOTO 30 
!    15 CONTINUE 
! !                                  CHECK FOR UERSET CALL                
!       DO 20 I = 1, 6 
!          IF (NAMUPK (I) .NE.NAMSET (I) ) GOTO 25 
!    20 END DO 
!       LEVOLD = LEVEL 
!       LEVEL = IER 
!       IER = LEVOLD 
!       IF (LEVEL.LT.0) LEVEL = 4 
!       IF (LEVEL.GT.4) LEVEL = 4 
!       GOTO 30 
!    25 CONTINUE 
!       IF (LEVEL.LT.4) GOTO 30 
! !                                  PRINT NON-DEFINED MESSAGE            
!       IF (IEQDF.EQ.1) WRITE (IOUNIT, 50) IER, NAMEQ, IEQ, NAMUPK 
!       IF (IEQDF.EQ.0) WRITE (IOUNIT, 50) IER, NAMUPK 
!    30 IEQDF = 0 
!       RETURN 
!    35 FORMAT(19H *** TERMINAL ERROR,10X,7H(IER = ,I3,                   &
!      &       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)                        
!    40 FORMAT(27H *** WARNING WITH FIX ERROR,2X,7H(IER = ,I3,            &
!      &       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)                        
!    45 FORMAT(18H *** WARNING ERROR,11X,7H(IER = ,I3,                    &
!      &       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)                        
!    50 FORMAT(20H *** UNDEFINED ERROR,9X,7H(IER = ,I5,                   &
!      &       20H) FROM IMSL ROUTINE ,6A1,A1,6A1)                        
! !                                                                       
! !                                  SAVE P FOR P = R CASE                
! !                                    P IS THE PAGE NAMUPK               
! !                                    R IS THE ROUTINE NAMUPK            
!    55 IEQDF = 1 
!       DO 60 I = 1, 6 
!    60 NAMEQ (I) = NAMUPK (I) 
!    65 RETURN 
!       END SUBROUTINE
!       
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
!  
!       SUBROUTINE UGETIO (IOPT, NIN, NOUT) 
!       INTEGER IOPT, NIN, NOUT 
!       INTEGER NIND, NOUTD 
!       DATA NIND / 5 /, NOUTD / 6 / 
!       IF (IOPT.EQ.3) GOTO 10 
!       IF (IOPT.EQ.2) GOTO 5 
!       IF (IOPT.NE.1) GOTO 9005 
!       NIN = NIND 
!       NOUT = NOUTD 
!       GOTO 9005 
!     5 NIND = NIN 
!       GOTO 9005 
!    10 NOUTD = NOUT 
!  9005 RETURN 
!       END SUBROUTINE
!                                                                  
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
!  
!       SUBROUTINE USPKD (PACKED, NCHARS, UNPAKD, NCHMTB) 
!       INTEGER NC, NCHARS, NCHMTB 
!       CHARACTER UNPAKD ( * ), IBLANK 
!       CHARACTER ( * ) PACKED 
!       DATA IBLANK / ' ' / 
! !                                  INITIALIZE NCHMTB                    
!       NCHMTB = 0 
! !                                  RETURN IF NCHARS IS LE ZERO          
!       IF (NCHARS.LE.0) RETURN 
! !                                  SET NC=NUMBER OF CHARS TO BE DECODED 
!       NC = MIN0 (129, NCHARS) 
!       READ (PACKED, 150) (UNPAKD (I), I = 1, NC) 
!   150 FORMAT (129A1) 
! !                                  CHECK UNPAKD ARRAY AND SET NCHMTB    
! !                                  BASED ON TRAILING BLANKS FOUND       
!       DO 200 N = 1, NC 
!          NN = NC - N + 1 
!          IF (UNPAKD (NN) .NE.IBLANK) GOTO 210 
!   200 END DO 
!       NN = 0 
!   210 NCHMTB = NN 
!       RETURN 
!       END SUBROUTINE
! 
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! !*************************************************************************!
! 
! ! subroutine for determination of e-vectors and e-values                
! ! from numerical recipes     
!                                            
!       SUBROUTINE JACOBI (A, N, NP, D, V, NROT) 
!       IMPLICIT none 
!       INTEGER NMAX, NP, N 
!       PARAMETER (NMAX = 2500) 
!       REAL(8) A (NP, NP), D (NP), V (NP, NP), B (NMAX), Z (NMAX) 
!       INTEGER IP, IQ, NROT, I, J 
!       REAL(8) sm, tresh, G, H, T, theta, C, s, tau 
!       REAL(8) tmp1, tmp2 
!                                                                         
!       DO 12 IP = 1, N 
!          DO 11 IQ = 1, N 
!             V (IP, IQ) = 0.d0 
!    11    END DO 
!          V (IP, IP) = 1.d0 
!    12 END DO 
!       DO 13 IP = 1, N 
!          B (IP) = A (IP, IP) 
!          D (IP) = B (IP) 
!          Z (IP) = 0.d0 
!    13 END DO 
!       NROT = 0 
!       DO 24 I = 1, 50 
!          SM = 0.d0 
!          DO 15 IP = 1, N - 1 
!             DO 14 IQ = IP + 1, N 
!                SM = SM + ABS (A (IP, IQ) ) 
!    14       END DO 
!    15    END DO 
!          IF (SM.EQ.0.) RETURN 
!          IF (I.LT.4) THEN 
!             TRESH = 0.2d0 * SM / N**2 
!          ELSE 
!             TRESH = 0.d0 
!          ENDIF 
!          DO 22 IP = 1, N - 1 
!             DO 21 IQ = IP + 1, N 
!                G = 100.d0 * ABS (A (IP, IQ) ) 
!                tmp1 = ABS (D (IP) ) + G 
!                tmp2 = ABS (D (IQ) ) + G 
!                IF ( (I.GT.4) .AND. (tmp1.EQ.ABS (D (IP) ) ) .AND. (     &
!                tmp2.EQ.ABS (D (IQ) ) ) ) THEN                           
!                   A (IP, IQ) = 0.d0 
!                ELSEIF (ABS (A (IP, IQ) ) .GT.TRESH) THEN 
!                   H = D (IQ) - D (IP) 
!                   tmp1 = ABS (H) + G 
!                   IF (tmp1.EQ.ABS (H) ) THEN 
!                      T = A (IP, IQ) / H 
!                   ELSE 
!                      THETA = 0.5d0 * H / A (IP, IQ) 
!                      T = 1. / (ABS (THETA) + SQRT (1. + THETA**2) ) 
!                      IF (THETA.LT.0.d0) T = - T 
!                   ENDIF 
!                   C = 1.d0 / SQRT (1.d0 + T**2) 
!                   S = T * C 
!                   TAU = S / (1.d0 + C) 
!                   H = T * A (IP, IQ) 
!                   Z (IP) = Z (IP) - H 
!                   Z (IQ) = Z (IQ) + H 
!                   D (IP) = D (IP) - H 
!                   D (IQ) = D (IQ) + H 
!                   A (IP, IQ) = 0.d0 
!                   DO 16 J = 1, IP - 1 
!                      G = A (J, IP) 
!                      H = A (J, IQ) 
!                      A (J, IP) = G - S * (H + G * TAU) 
!                      A (J, IQ) = H + S * (G - H * TAU) 
!    16             END DO 
!                   DO 17 J = IP + 1, IQ - 1 
!                      G = A (IP, J) 
!                      H = A (J, IQ) 
!                      A (IP, J) = G - S * (H + G * TAU) 
!                      A (J, IQ) = H + S * (G - H * TAU) 
!    17             END DO 
!                   DO 18 J = IQ + 1, N 
!                      G = A (IP, J) 
!                      H = A (IQ, J) 
!                      A (IP, J) = G - S * (H + G * TAU) 
!                      A (IQ, J) = H + S * (G - H * TAU) 
!    18             END DO 
!                   DO 19 J = 1, N 
!                      G = V (J, IP) 
!                      H = V (J, IQ) 
!                      V (J, IP) = G - S * (H + G * TAU) 
!                      V (J, IQ) = H + S * (G - H * TAU) 
!    19             END DO 
!                   NROT = NROT + 1 
!                ENDIF 
!    21       END DO 
!    22    END DO 
!          DO 23 IP = 1, N 
!             B (IP) = B (IP) + Z (IP) 
!             D (IP) = B (IP) 
!             Z (IP) = 0.d0 
!    23    END DO 
!    24 END DO 
!       PAUSE '50 iterations should never happen' 
!       RETURN 
!       END SUBROUTINE
!                                 
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


!************************************************************************
!** SEigensystem diagonalizes a complex symmetric n-by-n matrix.
!** Input: n, A = n-by-n matrix, complex symmetric
!** (only the upper triangle of A needs to be filled).
!** Output: d = vector of eigenvalues, U = transformation matrix
!** these fulfill diag(d) = U A U^T = U A U^-1 with U U^T = 1.
!************************************************************************


      subroutine SEigensystem(n, A,ldA, d, U,ldU, sort)
      implicit none
      integer           :: n, ldA, ldU, sort
      Complex(8)        :: A(ldA,ldA), U(ldU,ldU), d(ldA)
      real(8),parameter :: SYM_EPS=2.D0**(-52)
      integer           :: p, q, j
      Real(8)           :: red, off, thresh
      Complex(8)        :: delta, t, invc, s
      Complex(8)        :: x, y
      Complex(8)        :: ev(2,n)
      integer           :: sweep

      do p = 1, n
        ev(1,p) = 0
        ev(2,p) = A(p,p)
        d(p)    = ev(2,p)
      enddo

      do p = 1, n
        do q = 1, n
          U(q,p) = 0.d0
        enddo
        U(p,p) = 1.d0
      enddo

      red = .04D0/dble(n**4)

      do sweep = 1, 200

        off = 0
        do q = 2, n
          do p = 1, q - 1
            off = off + Sq(A(p,q))
          enddo
        enddo
        if( .not. off .gt. SYM_EPS ) goto 1

        thresh = 0
        if( sweep .lt. 4 ) thresh = off*red

        do q = 2, n
          do p = 1, q - 1
            delta = A(p,q)
            off = Sq(delta)
            if( sweep .gt. 4 .and. off .lt. SYM_EPS*(Sq(ev(2,p)) + Sq(ev(2,q))) ) then
              A(p,q) = 0.
            else if( off .gt. thresh ) then
              x = 0.5D0*(ev(2,p) - ev(2,q))
              y = sqrt(x**2 + delta**2)
              t = x - y
              s = x + y
              if( Sq(t) .lt. Sq(s) ) t = s

              t = delta/t
              delta = delta*t
              ev(1,p) = ev(1,p) + delta
              ev(2,p) = d(p) + ev(1,p)
              ev(1,q) = ev(1,q) - delta
              ev(2,q) = d(q) + ev(1,q)

              invc = sqrt(t**2 + 1.d0)
              s = t/invc
              t = t/(invc + 1.d0)

              do j = 1, p - 1
                x = A(j,p)
                y = A(j,q)
                A(j,p) = x + s*(y - t*x)
                A(j,q) = y - s*(x + t*y)
              enddo

              do j = p + 1, q - 1
                x = A(p,j)
                y = A(j,q)
                A(p,j) = x + s*(y - t*x)
                A(j,q) = y - s*(x + t*y)
              enddo

              do j = q + 1, n
                x = A(p,j)
                y = A(q,j)
                A(p,j) = x + s*(y - t*x)
                A(q,j) = y - s*(x + t*y)
              enddo

              A(p,q) = 0

              do j = 1, n
                x = U(p,j)
                y = U(q,j)
                U(p,j) = x + s*(y - t*x)
                U(q,j) = y - s*(x + t*y)
              enddo
            endif
          enddo
        enddo

        do p = 1, n
          ev(1,p) = 0
          d(p) = ev(2,p)
        enddo

      enddo

      print *, "Bad convergence in SEigensystem"

1     if( sort .eq. 0 ) return

! sort the eigenvalues by their real part

      do p = 1, n - 1
        j = p
        t = d(p)
        do q = p + 1, n
          if( sort*(dble(t) - dble(d(q))) .gt. 0 ) then
            j = q
            t = d(q)
          endif
        enddo

        if( j .ne. p ) then
          d(j) = d(p)
          d(p) = t
          do q = 1, n
            x = U(p,q)
            U(p,q) = U(j,q)
            U(j,q) = x
          enddo
        endif
      enddo

     contains

     Real(8) function Sq(c)
      Complex(8)  :: c
      Sq = dble(c*Conjg(c))
     end function

     end subroutine

!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!*************************************************************************!
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

end module

