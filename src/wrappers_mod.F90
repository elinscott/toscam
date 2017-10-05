! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris, Arash A. Mostofi and Peter D. Haynes
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module wrappers

!CW
#if defined (GPU_SPEEDUP_WRAPPER) || defined (FFTW3GPU)
 use fortran_cuda
#endif
!END CW

  implicit none

  private

  public :: wrappers_ddot
  public :: wrappers_dsygv_lt
  public :: wrappers_dsygv_lt_2
  public :: wrappers_dcopy
  public :: wrappers_invert_sym_matrix
  public :: wrappers_invert_sym_cmatrix
  public :: wrappers_vcos_sin
  public :: wrappers_dscal
  public :: wrappers_daxpy
  public :: wrappers_dgemm
  public :: wrappers_1d_fft
  public :: wrappers_dsyev_lt ! lpl
  public :: wrappers_dgesv
  public :: wrappers_dgelss
!CW
  public :: wrappers_zhegv
!END CW

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wrappers_dgemm(cc, & ! output
       aa, bb, num)               ! input

    !=========================================================================!
    ! This subroutine returns in cc the product C = A x B of square matrices  !
    ! A and B. The multiplication operation is parallelised. All processors   !
    ! hold copies of the whole matrices A and B and also C (on return from    !
    ! this subroutine).                                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  cc (out)     : matrix product (matrix C)                               !
    !  aa  (input)  : A square matrix                                         !
    !  bb  (input)  : B square matrix                                         !
    !  num (input)  : dimension, all matrices are num x num square matrices   !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !  1) All matrices are square but not necessarily symmetric               !
    !                                                                         !
    !  2) The dummy array cc is overwritten during the calculation so it must !
    !     under no circustances point to the same memory as aa or bb          !
    !-------------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 2/3/2004                            !
    !=========================================================================!

    use comms, only: comms_bcast, pub_my_node_id, pub_total_num_nodes
    use constants, only: DP,stdout
    use timer, only: timer_clock
    implicit none

    integer, intent(in) :: num
    real(kind=DP), intent(out):: cc(num, num) ! cks: not the same array as aa or bb!
    real(kind=DP), intent(in) :: aa(num, num)
    real(kind=DP), intent(in) :: bb(num, num)

    ! cks: <<local variables>>
    integer :: col_bb_batch_size
    integer :: col_bb_start
    integer :: col_bb_end
    integer :: local_num

    integer :: node_count
    integer :: col_node_start
    integer :: col_node_end

!CW
#ifdef GPU_SPEEDUP_WRAPPER
 CALL matmulcuda_r(aa,bb,cc,num,num,num)
 return
#endif
!END CW

    call timer_clock('wrappers_dgemm',1)

    ! cks: parallelisation of the dgemm - each node
    ! cks: does a subset of the columns of bb
    col_bb_batch_size =num/pub_total_num_nodes
    col_bb_start      =1 +pub_my_node_id*col_bb_batch_size
    col_bb_end        =col_bb_start +col_bb_batch_size -1
    if (pub_my_node_id .eq. (pub_total_num_nodes-1) ) col_bb_end= num
    local_num         =col_bb_end -col_bb_start +1

    call timer_clock('wrappers_dgemm_call_dgemm',1)

    ! cks: calculate and return aa x bb
    call dgemm('n', 'n', num, local_num, num, 1.0_DP, aa, num, &
         bb(:, col_bb_start: col_bb_end), num, 0.0_DP, &
         cc(:, col_bb_start: col_bb_end), num )

    call timer_clock('wrappers_dgemm_call_dgemm',2)


    call timer_clock('wrappers_dgemm_bcasts', 1)

    ! cks: every node braodcasts the columns of cc it has calculated
    ! cks: to the other nodes
    do node_count =0, pub_total_num_nodes -1

       col_node_start =1 +node_count*col_bb_batch_size
       col_node_end   =col_node_start +col_bb_batch_size -1
       if (node_count .eq. (pub_total_num_nodes-1) ) col_node_end= num
       local_num =col_node_end -col_node_start +1

       call comms_bcast(node_count, cc(:, &
            col_node_start : col_node_end), num*local_num)

    end do

    call timer_clock('wrappers_dgemm_bcasts', 2)

    call timer_clock('wrappers_dgemm',2)

  end subroutine wrappers_dgemm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function wrappers_ddot(n,x,incx,y,incy)

    !===========================================!
    ! ddot wrapper.                             !
    !-------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000. !
    !===========================================!

    use constants, only: DP
    implicit none

    real(kind=DP) :: wrappers_ddot
    integer, intent(in) :: n, incx, incy
    real(kind=DP), intent(in), dimension(:) :: x, y

    real(kind=DP) :: ddot

    wrappers_ddot=ddot(n,x,incx,y,incy)

  end function wrappers_ddot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wrappers_dsygv_lt(eigenvectors,eigenvalues, &
       overlap_square,num)

    !=======================================================!
    ! This subroutine solves the FC=SCe type of generalised !
    ! eigenvalue problem and obtains both eigenvalues and   !
    ! eigenvectors. F is stored as a lower triangle in a    !
    ! full square matrix. The eigenvalues are returned in   !
    ! ascending order.                                      !
    !-------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.             !
    !=======================================================!

    use comms, only: comms_abort, pub_on_root
    use constants, only: DP, stdout
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: num
    real(kind=DP), intent(inout) :: eigenvectors(:,:)
    real(kind=DP), intent(inout) :: overlap_square(:,:)
    real(kind=DP), intent(out) :: eigenvalues(:)

    ! Local Variables
    real(kind=DP), dimension(:), allocatable :: work
    integer :: info
    integer :: lda,ldb
    integer :: ierr ! error flag

    call timer_clock('wrappers_dsygv_lt', 1)
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering wrappers_dsygv_lt'
#endif

    ! ndmh: check arguments
    if (size(eigenvectors)<num) then
       call utils_abort('Error in wrappers_dsygv_lt: eigenvectors array is &
            &too small')
    end if
    if ((size(overlap_square,1)<num).or.(size(overlap_square,2)<num)) then
       call utils_abort('Error in wrappers_dsygv_lt: overlap matrix is &
            &too small')
    end if
    if ((size(eigenvectors,1)<num).or.(size(eigenvectors,2)<num)) then
       call utils_abort('Error in wrappers_dsygv_lt: eigenvectors matrix is &
            &too small')
    end if

    ! ndmh: find leading dimension of arrays
    lda = size(eigenvectors,1)
    ldb = size(overlap_square,1)

    allocate(work(3*num), stat=ierr)
    call utils_alloc_check('wrappers_dsygv_lt','work',ierr)

    ! cks: we want to solve the FC=SCe type of generalised eigenvalue
    !      problem and obtain both eigenvalues and eigenvectors.
    !      F is stored as a lower triangle in a full square matrix.
    !      The eigenvalues are returned in ascending order.
    call dsygv(1,'V','L',num,eigenvectors,lda,overlap_square,ldb,&
         eigenvalues,work,3*num,info)

    if (info.ne.0) then
       if (pub_on_root) write(stdout,'(a,i5)') 'DSYGV in subroutine &
            & wrappers_dsygv_lt returned info=',info
       if (pub_on_root) write(stdout,'(a)') 'ONETEP execution stops'
       call comms_abort
    end if

    deallocate(work,stat=ierr)
    call utils_dealloc_check('wrappers_dsygv_lt','work',ierr)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving wrappers_dsygv_lt'
#endif
    call timer_clock('wrappers_dsygv_lt', 2)


  end subroutine wrappers_dsygv_lt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wrappers_zhegv(eigenvectors,eigenvalues, &
       overlap_square,num)

    !=======================================================!
    ! This subroutine solves the FC=SCe type of generalised !
    ! eigenvalue problem and obtains both eigenvalues and   !
    ! eigenvectors. F is stored as a lower triangle in a    !
    ! full square matrix. The eigenvalues are returned in   !
    ! ascending order.                                      !
    !-------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.             !
    !=======================================================!

    use comms, only: comms_abort, pub_on_root
    use constants, only: DP, stdout
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: num
    complex(kind=DP), intent(inout) :: eigenvectors(:,:)
    complex(kind=DP), intent(inout) :: overlap_square(:,:)
    real(kind=DP), intent(out) :: eigenvalues(:)

    ! Local Variables
    complex(8),allocatable :: work(:)
    real(kind=DP), dimension(:), allocatable :: rwork
    integer :: info
    integer :: lda,ldb
    integer :: ierr ! error flag

    lda = size(eigenvectors,1)
    ldb = size(overlap_square,1)

    allocate(work(2*num-1),rwork(3*num-2), stat=ierr)

    call zhegv(1,'V','L',num,eigenvectors,lda,overlap_square,ldb,&
         eigenvalues,work,2*num-1,rwork,info)

    if(info.ne.0) then
       if (pub_on_root) write(stdout,'(a,i5)') 'DSYGV in subroutine &
            & wrappers_dsygv_lt returned info=',info
       if (pub_on_root) write(stdout,'(a)') 'ONETEP execution stops'
       call comms_abort
    endif

    deallocate(work,rwork,stat=ierr)

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!END CW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wrappers_dsygv_lt_2(eigenvectors, eigenvalues, &
       overlap_square, num)

    !=======================================================!
    ! This subroutine solves the KSC=Cn type of generalised !
    ! eigenvalue problem and obtains both eigenvalues and   !
    ! eigenvectors. F is stored as a lower triangle in a    !
    ! full square matrix. The eigenvalues are returned in   !
    ! ascending order.                                      !
    !-------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.             !
    !=======================================================!

    use comms, only: comms_abort, pub_on_root
    use constants, only: DP, stdout
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: num
    real(kind=DP), intent(inout) :: eigenvectors(:,:)
    real(kind=DP), intent(inout) :: overlap_square(:,:)
    real(kind=DP), intent(out) :: eigenvalues(:)

    ! Local Variables
    real(kind=DP), dimension(:), allocatable :: work
    integer :: info
    integer :: lda,ldb
    integer :: ierr ! error flag


    call timer_clock('wrappers_dsygv_lt_2', 1)
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering wrappers_dsygv_lt_2'
#endif

    ! ndmh: check arguments
    if (size(eigenvectors)<num) then
       call utils_abort('Error in wrappers_dsygv_lt_2: eigenvectors array is &
            &too small')
    end if
    if ((size(overlap_square,1)<num).or.(size(overlap_square,2)<num)) then
       call utils_abort('Error in wrappers_dsygv_lt_2: overlap matrix is &
            &too small')
    end if
    if ((size(eigenvectors,1)<num).or.(size(eigenvectors,2)<num)) then
       call utils_abort('Error in wrappers_dsygv_lt_2: eigenvectors matrix is &
            &too small')
    end if

    ! ndmh: find leading dimension of arrays
    lda = size(eigenvectors,1)
    ldb = size(overlap_square,1)

    allocate(work(3*num),stat=ierr)
    call utils_alloc_check('wrappers_dsygv_lt_2','work',ierr)

    ! cks: we want to solve the KSC=Cn type of generalised eigenvalue
    !      problem and obtain both eigenvalues and eigenvectors.
    !      K is stored as a lower triangle in a full square matrix.
    !      The eigenvalues are returned in ascending order.
    call dsygv(2,'V','L',num,eigenvectors,lda,overlap_square,ldb,&
         eigenvalues,work,3*num,info)

    if (info.ne.0) then
       if (pub_on_root) write(stdout,'(a,i5)') 'DSYGV in subroutine &
            & wrappers_dsygv_lt_2 returned info=',info
       if (pub_on_root) write(stdout,'(a)') 'ONETEP execution stops'
       call comms_abort
    end if

    deallocate(work,stat=ierr)
    call utils_dealloc_check('wrappers_dsygv_lt_2','work',ierr)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving wrappers_dsygv_lt_2'
#endif
    call timer_clock('wrappers_dsygv_lt_2', 2)

  end subroutine wrappers_dsygv_lt_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wrappers_invert_sym_matrix(matrix,num)

    !=======================================================!
    ! This subroutine returns the inverse of a full square  !
    ! symmatric matrix.                                     !
    !-------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 17/4/2001         !
    !=======================================================!

    use comms, only: comms_abort, pub_on_root
    use constants, only: DP, stdout
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    integer, intent(in) :: num
    real(kind=DP), intent(inout) :: matrix(:,:) ! inverse on exit

    ! cks: internal declarations
    integer :: work_length, info, row, col, lda
    integer, allocatable, dimension(:) :: ipiv
    real(kind=DP), allocatable, dimension(:) :: work_array
    integer :: ierr ! error flag

!CW
#ifdef GPU_SPEEDUP_WRAPPER
    call magma_fortran_double_(num,matrix)
    return
#endif
!END CW

    !call timer_clock('wrappers_invert_sym_matrix',1)

    work_length=3*num
    lda = size(matrix,1)

    if (lda<num) then
       if (pub_on_root) write(stdout,'(a)') 'ERROR in &
            &wrappers_invert_sym_matrix: invalid matrix sizes'
       call comms_abort
    endif

    allocate(ipiv(num),stat=ierr)
    call utils_alloc_check('wrappers_invert_sym_matrix','ipiv',ierr)
    allocate(work_array(work_length),stat=ierr)
    call utils_alloc_check('wrappers_invert_sym_matrix','work_array',ierr)

    ! cks: compute the factorization of a real symmetric matrix A using
    !      the Bunch-Kaufman diagonal pivoting method
    call dsytrf('L', num, matrix, lda, ipiv, work_array, work_length, info)
    if (info.ne.0) then
       if (pub_on_root) write(stdout,'(a,i5)') 'ERROR in &
          &wrappers_invert_sym_matrix: Problem with dsytrf, info=',info
           do row=1,num
              write(*,*) (matrix(row,col),col=1,row)
           enddo
!CW
       return 
    !  call comms_abort
!END CW

    endif

    ! cks: compute the inverse of a real symmetric indefinite matrix A using
    ! cks: the factorization A = U*D*U**T or A = L*D*L**T computed by DSYTRF
    call dsytri('L', num, matrix, lda, ipiv, work_array, info )
    if (info.ne.0) then
       if (pub_on_root) write(stdout,'(a,i5)') 'ERROR in &
            &wrappers_invert_sym_matrix: Problem with dsytri, info=',info
       call comms_abort
    endif

    ! cks: fill in the upper triangular part of the symmetric matrix
    do row=1,num
       do col=1,row-1
          matrix(col,row)=matrix(row,col)
       enddo
    enddo

    deallocate(work_array,stat=ierr)
    call utils_dealloc_check('wrappers_invert_sym_matrix','work_array',ierr)
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('wrappers_invert_sym_matrix','ipiv',ierr)

    !call timer_clock('wrappers_invert_sym_matrix',2)

  end subroutine wrappers_invert_sym_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine wrappers_invert_sym_cmatrix(cmatrix,num)

    !=======================================================!
    ! This subroutine returns the inverse of a full square  !
    ! complex indefinite symmetric matrix.                  !
    !-------------------------------------------------------!
    ! Written by David O'Regan in April 2009 based on       !
    ! wrappers_invert_sym_matrix by Chris-Kriton Skylaris   !
    !=======================================================!

    use comms, only: comms_abort, pub_on_root
    use constants, only: DP, stdout
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    integer, intent(in) :: num
    complex(kind=DP), intent(inout) :: cmatrix(num,num) ! inverse on exit

    ! cks: internal declarations
    integer :: work_length, info, row, col
    integer, allocatable, dimension(:) :: ipiv
    complex(kind=DP), allocatable, dimension(:) :: work_array
    integer :: ierr ! error flag

!CW
#ifdef GPU_SPEEDUP_WRAPPER
    call magma_fortran_comp_(num,cmatrix)
    return
#endif
!END CW

    !call timer_clock('wrappers_invert_sym_cmatrix',1)

    work_length=3*num

    allocate(ipiv(num),stat=ierr)
    call utils_alloc_check('wrappers_invert_sym_cmatrix','ipiv',ierr)
    allocate(work_array(work_length),stat=ierr)
    call utils_alloc_check('wrappers_invert_sym_cmatrix','work_array',ierr)

    ! cks: compute the factorization of a complex symmetric matrix A using
    !      the Bunch-Kaufman diagonal pivoting method
    call zsytrf('L', num, cmatrix, num, ipiv, work_array, work_length, info)
    if (info.ne.0) then
       if (pub_on_root) write(stdout,'(a,i5)') 'ERROR in &
            &wrappers_invert_sym_cmatrix: Problem with zsytrf, info=',info
       call comms_abort
    endif

    ! cks: compute the inverse of a complex symmetric indefinite matrix A using
    ! cks: the factorization A = U*D*U**T or A = L*D*L**T computed by CSYTRF
    call zsytri('L', num, cmatrix, num, ipiv, work_array, info )
    if (info.ne.0) then
       if (pub_on_root) write(stdout,'(a,i5)') 'ERROR in &
            &wrappers_invert_sym_cmatrix: Problem with zsytri, info=',info
       call comms_abort
    endif

    ! cks: fill in the upper triangular part of the complex symmetric matrix
    do row=1,num
       do col=1,row-1
          cmatrix(col,row)=cmatrix(row,col)
       enddo
    enddo

    deallocate(work_array,stat=ierr)
    call utils_dealloc_check('wrappers_invert_sym_cmatrix','work_array',ierr)
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('wrappers_invert_sym_cmatrix','ipiv',ierr)

    !call timer_clock('wrappers_invert_sym_cmatrix',2)

  end subroutine wrappers_invert_sym_cmatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wrappers_dcopy(nn,xx,incx,yy,incy)

    !================================================!
    ! dcopy wrapper.                                 !
    !------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/7/2001. !
    !================================================!

    ! cks: BLAS vector copy. yy<--xx
    use constants, only: DP
    implicit none

    integer :: nn, incx, incy
    real(kind=DP), intent(in) :: xx(nn)
    real(kind=DP), intent(out) :: yy(nn)

    call dcopy(nn,xx,incx,yy,incy)

  end subroutine wrappers_dcopy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wrappers_vcos_sin(xx,incx,yy,incy,zz,incz,nn)

    !===========================================!
    ! Wrapper for vector cosine and sine        !
    !-------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2001. !
    !===========================================!

    ! wrapper for BLAS vector cosine and sine
    use constants, only: DP
    implicit none

    integer, intent(in) :: incx, incy, incz, nn
    real(kind=DP), intent(in) :: xx(nn)
    real(kind=DP), intent(out) :: yy(nn), zz(nn)

#ifdef ALPHA
    call vcos_sin(xx,incx,yy,incy,zz,incz,nn)
#else
#ifdef SUN
stop "hack by cks to compile on the sun: vsincos does not exist. STOP!"
!    call vsincos(nn,xx,incx,zz,incz,yy,incy)
#else
    integer i
    complex(KIND=DP) eiarg

    if (incx == 1 .and. incy == 1 .and. incz == 1) then
       do i=1,nn
          eiarg = exp(cmplx(0.0_DP,xx(i),DP))
          yy(i) = real(eiarg,DP)
          zz(i) = aimag(eiarg)
       end do
    else
       do i=1,nn
          eiarg = exp(cmplx(0.0_DP,xx(1+(i-1)*incx),DP))
          yy(1+(i-1)*incy) = real(eiarg,DP)
          zz(1+(i-1)*incz) = aimag(eiarg)
       end do
    end if
#endif
#endif

  end subroutine wrappers_vcos_sin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wrappers_dscal(nn, aa, xx, incx )

    !================================================!
    ! dscal wrapper.                                 !
    !------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 27/7/2001. !
    !================================================!

    ! cks: wrapper for BLAS scale of vector  by scalar
    use constants, only: DP
    implicit none

    integer, intent(in) :: nn, incx
    real(kind=DP), intent(in) :: aa
    real(kind=DP), intent(out) :: xx(nn)

    call dscal(nn, aa, xx, incx)

  end subroutine wrappers_dscal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wrappers_daxpy(nn,aa,xx,incx,yy,incy)

    !================================================!
    ! daxpy wrapper.                                 !
    !------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 22/8/2001. !
    !================================================!

    use constants, only: DP
    implicit none

    integer, intent(in) :: nn, incx, incy
    real(kind=DP), intent(in) :: aa
    real(kind=DP), intent(in) :: xx(nn)
    real(kind=DP), intent(inout) :: yy(nn)

    call daxpy(nn,aa,xx,incx,yy,incy)

  end subroutine wrappers_daxpy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wrappers_1d_fft(dir, isize, in_array, out_array)

    !================================================!
    ! Wrapper for 1-dimensional FFT                  !
    !------------------------------------------------!
    ! Written by David O'Regan in March 2009         !
    ! based on fftbench by A. Mostofi                !
    !================================================!

    use comms, only: pub_on_root
    use constants, only: DP, stdout
    use timer, only: timer_clock

    implicit none

#ifdef FFTW
  integer, parameter :: FFTW_FORWARD = -1
  integer, parameter :: FFTW_BACKWARD = 1
  integer, parameter :: FFTW_REAL_TO_COMPLEX = -1
  integer, parameter :: FFTW_COMPLEX_TO_REAL = 1
  integer, parameter :: FFTW_ESTIMATE = 0
  integer, parameter :: FFTW_MEASURE = 1
  integer, parameter :: FFTW_OUT_OF_PLACE = 0
  integer, parameter :: FFTW_IN_PLACE = 8
  integer, parameter :: FFTW_USE_WISDOM = 16
  integer :: the_plan
#endif

#ifdef FFTW3
  integer, parameter :: FFTW_FORWARD = -1
  integer, parameter :: FFTW_BACKWARD = 1
  !integer, parameter :: FFTW_MEASURE = 0
  integer, parameter :: FFTW_ESTIMATE = 64
  integer :: the_plan
#endif

!CW
#ifdef FFTW3GPU
  integer            :: the_plan
  integer            :: CUFFT_FORWARD, CUFFT_INVERSE
  integer,parameter  :: batch=1
  parameter(CUFFT_FORWARD=-1, CUFFT_INVERSE=1)
#endif
!END CW

  character, intent(in) :: dir
  integer, intent(in) :: isize
  complex(kind=DP), intent(in) :: in_array(isize)
  complex(kind=DP), intent(inout) :: out_array(isize)
  integer :: zwork_len
  complex(kind=DP), allocatable :: zwork(:)

!CW
    write(*,*) 'CALLING WRAPPER 1D FFT'
#ifdef FFTW3GPU
    write(*,*) 'please debug and crosscheck wrt CPU, this is experimental'
    write(*,*) 'please also remove stop point in code to continue'
    stop
#endif
!END CW

    call timer_clock('wrappers_1d_fft', 1)
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering wrappers_1d_fft'
#endif

     ! Copy original into out_array
     out_array = in_array

     ! Initialise FFT routine
     call fourier_1d_init(dir,isize)

     ! Do forwards or backwards FFT
     call fourier_1d_apply(dir,isize,out_array)

     ! Finalise FFT routine
     call fourier_1d_exit

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving wrappers_1d_fft'
#endif
    call timer_clock('wrappers_1d_fft', 2)


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

contains

  subroutine fourier_1d_init(dir,n)

    use utils, only: utils_alloc_check

    implicit none

    character, intent(in) :: dir
    integer, intent(in) :: n

    ! Local variables

    integer :: ierr

    zwork_len = n

       allocate(zwork(zwork_len),stat=ierr)
       call utils_alloc_check('fourier_1d_fft','zwork',ierr)

    ! Platform-dependent initialisation

#ifdef FFTW
    if (dir == 'B' .or. dir == 'b') then
    call fftw_f77_create_plan(the_plan,n, &
         FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
    else
    call fftw_f77_create_plan(the_plan,n, &
         FFTW_FORWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
    endif
#endif

#ifdef FFTW3
    if (dir == 'B' .or. dir == 'b') then
    call dfftw_plan_dft_1d(the_plan,n,zwork,zwork,FFTW_BACKWARD,FFTW_ESTIMATE)
    else
    call dfftw_plan_dft_1d(the_plan,n,zwork,zwork,FFTW_FORWARD,FFTW_ESTIMATE)
    endif
#endif

!CW
#ifdef FFTW3GPU
    call cufftplan1d_Z2Z(the_plan, n, BATCH)
#endif
!END CW

  end subroutine fourier_1d_init

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_1d_apply(dir,n,array)

    use constants, only: DP

    implicit none

    character, intent(in) :: dir
    integer, intent(in) :: n
    complex(kind=DP), intent(inout) :: array(n)

    real(kind=DP) :: scale

    ! Platform-dependent call

#ifdef FFTW
       call fftw_f77_one(the_plan,array,zwork)
    if (dir == 'B' .or. dir == 'b') then
       scale = 1.0_dp / n
       array = array * scale
    end if
#endif

#ifdef FFTW3
       call dfftw_execute_dft(the_plan,array,array)
    if (dir == 'B' .or. dir == 'b') then
       scale = 1.0_DP / n
       array = array * scale
    end if
#endif

!CW
#ifdef FFTW3GPU
    if (dir == 'B' .or. dir == 'b') then
    call cufftexecz2z(the_plan, array, array, CUFFT_INVERSE,n,BATCH)
    else
    call cufftexecz2z(the_plan, array, array, CUFFT_FORWARD,n,BATCH)
    endif
    if (dir == 'B' .or. dir == 'b') then
       scale = 1.0_DP / n
       array = array * scale
    end if
#endif
!END CW

  end subroutine fourier_1d_apply

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine fourier_1d_exit

    use utils, only: utils_dealloc_check

    implicit none

    integer :: ierr

    ! Platform-dependent finalisation

#ifdef FFTW
    call fftw_f77_destroy_plan(the_plan)
#endif

#ifdef FFTW3
    call dfftw_destroy_plan(the_plan)
#endif

!CW
#ifdef FFTW3GPU
    call cufftdestroy(the_plan)
#endif
!END CW

    ! Workspace deallocation

    if (allocated(zwork)) then
       deallocate(zwork,stat=ierr)
       call utils_dealloc_check('fourier_1d_exit','zwork',ierr)
    endif

  end subroutine fourier_1d_exit

  end subroutine wrappers_1d_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wrappers_dsyev_lt(eigenvectors,eigenvalues,num)
  ! lpl: Solves normal symmetric eigenvalue equation Av=av

    use comms, only: comms_abort, pub_on_root
    use constants, only: DP, stdout
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    integer, intent(in) :: num
    real(kind=DP), intent(inout) :: eigenvectors(num,num)
    real(kind=DP), intent(out) :: eigenvalues(num)

    real(kind=DP), dimension(:), allocatable :: work
    integer :: info
    integer :: ierr ! error flag

    call timer_clock('wrappers_dsyev_lt', 1)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering wrappers_dsyev_lt'
#endif

    allocate(work(5*num), stat=ierr)
    call utils_alloc_check('wrappers_dsyev_lt','work',ierr)

    call dsyev('V','U',num,eigenvectors,num,eigenvalues,work,5*num,info)

    if (info.ne.0) then
       if (pub_on_root) write(stdout,'(a,i5)') 'DSYEV in subroutine &
            & wrappers_dsyev_lt returned info=',info
       if (pub_on_root) write(stdout,'(a)') 'ONETEP execution stops'
       call comms_abort
    end if

    deallocate(work,stat=ierr)
    call utils_dealloc_check('wrappers_dsyev_lt','work',ierr)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving wrappers_dsyev_lt'
#endif
    call timer_clock('wrappers_dsyev_lt', 2)

  end subroutine wrappers_dsyev_lt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wrappers_dgesv(xx, & ! output
       aa, bb, num, nrhs)         ! input

    !=========================================================================!
    ! This subroutine solves a linear system AX=B where A is a general matrix.!
    ! B and X are N by NRHS matrices. NRHS is the number of right hand sides. !
    ! Uses LAPACK.                                                            !
    !-------------------------------------------------------------------------!
    ! Written by Quintin Hill on 21/12/2007  for dposv                        !
    ! Modified by Quintin Hill on 30/01/2007 for dgesv                        !
    !=========================================================================!

    use comms, only: comms_abort
    use constants, only: stdout, dp
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    integer, intent(in) :: num
    integer, intent(in) :: nrhs
    real(kind=DP), intent(out):: xx(num,nrhs) !X
    real(kind=DP), intent(in) :: aa(num, num) !A
    real(kind=DP), intent(in) :: bb(num,nrhs) !B
    integer, allocatable :: pivots(:)
    real(kind=DP), allocatable :: aainternal(:,:)
    integer :: ierr
    integer :: info

    allocate(aainternal(num,num),stat=ierr)
    call utils_alloc_check('wrappers_dgesv','aainternal',ierr)
    allocate(pivots(num), stat=ierr)
    call utils_alloc_check('wrappers_dgesv','pivots',ierr)

    aainternal = aa
    xx = bb

    call timer_clock('wrappers_dgesv',1)

        call dgesv(num, nrhs, aainternal, num, pivots, xx, num, info)

        if (info/=0) then
        write(stdout,'(a,i4)') &
            'Error in wrapper_dgesv: &
            &dgesv returned abnormal exit status ', info
       call comms_abort
       end if


    call timer_clock('wrappers_dgesv',2)

    deallocate(aainternal, stat=ierr)
    call utils_dealloc_check('wrappers_dgesv','aainternal',ierr)
    deallocate(pivots, stat=ierr)
    call utils_dealloc_check('wrappers_dgesv','pivots',ierr)


  end subroutine wrappers_dgesv

  subroutine wrappers_dgelss(xx, & ! output
       aa, bb, num, nrhs)         ! input

    !=========================================================================!
    ! This subroutine solves a linear system AX=B where A is a general matrix.!
    ! B and X are N by NRHS matrices. NRHS is the number of right hand sides. !
    ! Uses LAPACK.                                                            !
    !-------------------------------------------------------------------------!
    ! Written by Quintin Hill on 27/01/2009  for dgelss                       !
    !=========================================================================!

    use comms, only: comms_abort
    use constants, only: stdout, dp
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    implicit none

    integer, intent(in) :: num
    integer, intent(in) :: nrhs
    real(kind=DP), intent(out):: xx(num,nrhs) !X
    real(kind=DP), intent(in) :: aa(num, num) !A
    real(kind=DP), intent(in) :: bb(num,nrhs) !B
    real(kind=DP), allocatable :: ss(:)
    real(kind=DP), allocatable :: aainternal(:,:)
    real(kind=DP), allocatable :: work(:)
    real(kind=DP) :: rcond
    integer :: rank
    integer :: ierr
    integer :: info
    integer :: lwork

    lwork = 3*num + max(2*num,nrhs)

    allocate(aainternal(num,num),stat=ierr)
    call utils_alloc_check('wrappers_dgelss','aainternal',ierr)
    allocate(ss(num), stat=ierr)
    call utils_alloc_check('wrappers_dgelss','ss',ierr)
    allocate(work(lwork), stat=ierr)
    call utils_alloc_check('wrappers_dgelss','work',ierr)

    rcond = -1.0_DP
    aainternal = aa
    xx = bb

    call timer_clock('wrappers_dgelss',1)

    call dgelss(num, num, nrhs, aainternal, num, xx, num, ss, rcond, rank, &
         work, lwork, info)

    if (info/=0) then
       write(stdout,'(a,i4)') &
            'Error in wrapper_dgelss: &
            &dgelss returned abnormal exit status ', info
       call comms_abort
    end if

    call timer_clock('wrappers_dgelss',2)

    deallocate(work, stat=ierr)
    call utils_dealloc_check('wrappers_dgelss','work',ierr)
    deallocate(ss, stat=ierr)
    call utils_dealloc_check('wrappers_dgelss','ss',ierr)
    deallocate(aainternal, stat=ierr)
    call utils_dealloc_check('wrappers_dgelss','aainternal',ierr)

  end subroutine wrappers_dgelss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module wrappers
