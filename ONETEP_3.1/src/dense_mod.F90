! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Dense Matrix Algebra Module
!
!   The subroutines in this file were written by
!
!   Nicholas D.M. Hine
!
!   Thomas Young Centre,
!   Department of Materials,
!   Imperial College London,
!   Exhibition Road, London SW7 2AZ
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dense

  use constants, only: DP, stdout

!CW
#ifdef GPU_SPEEDUP_DENSE
 use magma
 use fortran_cuda
#endif

  implicit none

  private

  ! Type definition for the structure of a dense matrix
  type DEM

    ! Variables required for both LAPACK and ScaLAPACK matrices
    logical :: iscmplx
    !logical :: islocal TODO - support for local matrices
    integer :: nrows, mcols
    real(kind=DP),allocatable :: dmtx(:,:)
    complex(kind=DP),allocatable :: zmtx(:,:)

#ifdef SCALAPACK
    ! Variables only required for the BLACS distributed matrices used
    ! by ScaLAPACK
    integer :: blacs_desc(50)
    integer :: blacs_nb
    integer :: blacs_ld
    integer :: blacs_ncol
#endif

  end type DEM

  ! Public Subroutines and types
  public :: DEM
  public :: dense_init
  public :: dense_exit
  public :: dense_create
  public :: dense_destroy
  public :: dense_product
  public :: dense_eigensolve
  public :: dense_invert
  public :: dense_convert
  public :: dense_axpy
  public :: dense_rank1_update
  public :: dense_scale
  public :: dense_copy
  public :: dense_get_element
  public :: dense_put_element
  public :: dense_get_col
  public :: dense_put_col

  ! AXPY routines
  interface dense_axpy
     module procedure dense_axpy_real
     module procedure dense_axpy_complex
  end interface

  ! Rank one update routines
  interface dense_rank1_update
     module procedure dense_rank1_update_real
     module procedure dense_rank1_update_complex
  end interface

  ! Scale and shift routines
  interface dense_scale
     module procedure dense_scale_real
     module procedure dense_scale_complex
  end interface

  ! Conversion routines
  interface dense_convert
     module procedure dense_convert_densetosparse
     module procedure dense_convert_sparsetodense
  end interface

  ! Element operation routines
  interface dense_get_element
     module procedure dense_get_element_real
     module procedure dense_get_element_complex
  end interface
  interface dense_put_element
     module procedure dense_put_element_real
     module procedure dense_put_element_complex
  end interface
  interface dense_get_col
     module procedure dense_get_col_real
     module procedure dense_get_col_complex
  end interface
  interface dense_put_col
     module procedure dense_put_col_real
     module procedure dense_put_col_complex
  end interface

#ifdef SCALAPACK
  ! Module variables describing BLACS context
  integer :: blacs_nprow, blacs_npcol
  integer :: blacs_myrow, blacs_mycol
  integer :: blacs_iam, blacs_nprocs
  integer :: blacs_context
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine initialises the variables in the dense matrix algebra      !
  ! module.                                                                    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  !============================================================================!

  subroutine dense_init

    use comms, only: pub_total_num_nodes

    implicit none

#ifdef SCALAPACK

    ! Choose a suitable value for nprow
    blacs_npcol = pub_total_num_nodes
    do blacs_nprow=int(sqrt(real(pub_total_num_nodes,kind=DP))),1,-1
       if (int(pub_total_num_nodes/blacs_nprow)*blacs_nprow== &
            pub_total_num_nodes) then
          blacs_npcol = pub_total_num_nodes/blacs_nprow
          exit
       end if
    end do

    ! Setup BLACS grid
    call blacs_pinfo(blacs_iam, blacs_nprocs)
    if (blacs_nprocs < 1) then
       call blacs_setup(blacs_iam, blacs_nprow*blacs_npcol)
    end if
    call blacs_get(-1, 0, blacs_context)
    call blacs_gridinit(blacs_context, 'R', blacs_nprow, blacs_npcol)
    call blacs_gridinfo(blacs_context, blacs_nprow, blacs_npcol, blacs_myrow, &
         blacs_mycol)

#endif

  end subroutine dense_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine cleans up variables associated with the dense matrix       !
  ! algebra module.                                                            !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  !============================================================================!

  subroutine dense_exit

    implicit none

#ifdef SCALAPACK
    ! Remove BLACS context
    call blacs_gridexit(blacs_context)
#endif

  end subroutine dense_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine allocates storage for a dense matrix, and, if applicable,  !
  ! sets up a BLACS descriptor for it.                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The dense matrix to create                            !
  !   nrows    (input) : The number of rows in the matrix                      !
  !   mcols    (input) : The number of columns in the matrix                   !
  !   iscmplx  (input) : Whether the matrix is to be complex                   !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  !============================================================================!

  subroutine dense_create(mat,nrows,mcols,iscmplx)

    use comms, only: comms_abort, pub_on_root
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: mat
    integer, intent(in) :: nrows
    integer, intent(in) :: mcols
    logical, intent(in), optional :: iscmplx

    ! Local Variables
    logical :: loc_iscmplx ! Local copy of optional iscmplx argument
    integer :: ierr        ! Error flag
#ifdef SCALAPACK
    integer, external :: numroc
#endif

    ! Set default for optional argument
    loc_iscmplx = .false.
    if (present(iscmplx)) loc_iscmplx = iscmplx
    mat%iscmplx = loc_iscmplx
    mat%nrows = nrows
    mat%mcols = mcols

    ! Check arguments
    if (nrows<0) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_create: nrows < 0'
       call comms_abort
    end if
    if (mcols<0) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_create: mcols < 0'
       call comms_abort
    end if

#ifdef SCALAPACK

    ! Choose a suitable blocking factor
    mat%blacs_nb = min(mat%nrows/blacs_nprow,mat%mcols/blacs_npcol)

    ! Initialise ld and ncol
    mat%blacs_ld = numroc(mat%nrows,mat%blacs_nb,blacs_myrow,0,blacs_nprow)
    mat%blacs_ncol = numroc(mat%mcols,mat%blacs_nb,blacs_mycol,0,blacs_npcol)

    ! Create matrix description
    call descinit(mat%blacs_desc, mat%nrows, mat%mcols, mat%blacs_nb, mat%blacs_nb, &
         0, 0, blacs_context, mat%blacs_ld, ierr)
    ! Check return value ierr
    if (ierr/=0) then
       if (pub_on_root) write(stdout,'(a,i5)') 'Error in dense_create: &
            &BLACS descinit call returned ierr=',ierr
       call comms_abort
    end if

    ! Allocate pointers
    if (mat%iscmplx) then
       allocate(mat%zmtx(mat%blacs_ld,mat%blacs_ncol),stat=ierr)
       call utils_alloc_check('dense_create','mat%zmtx',ierr)
    else
       allocate(mat%dmtx(mat%blacs_ld,mat%blacs_ncol),stat=ierr)
       call utils_alloc_check('dense_create','mat%dmtx',ierr)
    end if

#else

    ! Allocate pointers
    if (mat%iscmplx) then
       allocate(mat%zmtx(mat%nrows,mat%mcols),stat=ierr)
       call utils_alloc_check('dense_create','mat%zmtx',ierr)
    else
       allocate(mat%dmtx(mat%nrows,mat%mcols),stat=ierr)
       call utils_alloc_check('dense_create','mat%dmtx',ierr)
    end if

#endif

  end subroutine dense_create


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine deallocates storage for a dense matrix.                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat     (input) : The dense matrix to deallocate.                        !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  !============================================================================!

  subroutine dense_destroy(mat)

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: mat

    ! Local Variables
    integer :: ierr        ! Error flag

#ifdef SCALAPACK
    ! Deallocate pointers
    if (mat%iscmplx) then
       deallocate(mat%zmtx,stat=ierr)
       call utils_dealloc_check('dense_destroy','mat%zmtx',ierr)
    else
       deallocate(mat%dmtx,stat=ierr)
       call utils_dealloc_check('dense_destroy','mat%dmtx',ierr)
    end if
#else
    ! Deallocate pointers
    if (mat%iscmplx) then
       deallocate(mat%zmtx,stat=ierr)
       call utils_dealloc_check('dense_destroy','mat%zmtx',ierr)
    else
       deallocate(mat%dmtx,stat=ierr)
       call utils_dealloc_check('dense_destroy','mat%dmtx',ierr)
    end if
#endif

  end subroutine dense_destroy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine calculates the product of two dense matrices A and B and   !
  ! stores it in the dense matrix C.                                           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   cmat      (inout) : The dense matrix C                                   !
  !   amat      (in)    : The dense matrix A                                   !
  !   bmat      (in)    : The dense matrix B                                   !
  !   transpose_bmat (in,optional)  : Whether to transpose the B matrix        !
  !   first_k, last_k (in,optional) : Range of k-values to sum over in product !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23/02/2010, based partially on wrappers_dgemm,   !
  ! which was written by Chris-Kriton Skylaris on 2/3/2004.                    !
  !============================================================================!

  subroutine dense_product(cmat,amat,bmat,transpose_amat,transpose_bmat,first_k,last_k)

    use comms, only: comms_abort, comms_bcast, pub_my_node_id, pub_on_root, &
         pub_total_num_nodes
    implicit none

    ! Arguments
    type(DEM), intent(inout) :: cmat
    type(DEM), intent(in) :: amat
    type(DEM), intent(in) :: bmat
    logical, optional, intent(in) :: transpose_amat
    logical, optional, intent(in) :: transpose_bmat
    integer, optional, intent(in) :: first_k, last_k

    ! Local Variables
    integer :: m,n,k ! sizes of the matrices (see DGEMM manpages)
#ifndef SCALAPACK
    integer :: batch_size, local_start, local_end, local_n
    integer :: ldb
    integer :: node
#endif
    integer :: ia,ja,ib,jb
    character :: transa,transb
!CW
#ifdef GPU_SPEEDUP_DENSE
    real(8),allocatable    :: matar(:,:),matbr(:,:),matcr(:,:),tmpr(:,:)
    complex(8),allocatable :: matac(:,:),matbc(:,:),matcc(:,:),tmpc(:,:)
    integer                :: na1,na2,nb1,nb2
    logical                :: transpose_amat_,transpose_bmat_

     na1 = amat%nrows; na2 = amat%mcols; nb1 = bmat%nrows; nb2 = bmat%mcols;

      if (present(first_k).and.present(last_k)) then
        if(first_k/=1) goto 37 ; na2=last_k; nb2=last_k
      endif 

      transpose_amat_=.false.
      transpose_bmat_=.false.
     if(present(transpose_amat))then
      transpose_amat_=transpose_amat
     endif
     if(present(transpose_bmat))then
      transpose_bmat_=transpose_bmat
     endif
     if(pub_my_node_id==0)then
        if(.not.amat%iscmplx)then 
         if(.not.transpose_amat_)then
          allocate(matar(na1,na2))
          matar=amat%dmtx(1:na1,1:na2)
         else
          allocate(matar(na1,na2))
          matar=amat%dmtx(1:na1,1:na2)
          allocate(tmpr(na2,na1)); tmpr=transpose(matar)
          deallocate(matar);allocate(matar(na2,na1));
          matar=tmpr;deallocate(tmpr)
         endif
         if(.not.transpose_bmat_)then
          allocate(matbr(nb1,nb2))
          matbr=bmat%dmtx(1:nb1,1:nb2)
         else
          allocate(matbr(nb1,nb2))
          matbr=bmat%dmtx(1:nb1,1:nb2)
          allocate(tmpr(nb2,nb1)); tmpr=transpose(matbr)
          deallocate(matbr);allocate(matbr(nb2,nb1));
          matbr=tmpr;deallocate(tmpr)
         endif
         allocate(matcr(size(matar,1),size(matbr,2)))
         write(*,*) 'MATMUL GPU REAL'
         write(*,*) 'shape matrix a  = ', amat%nrows,amat%mcols
         write(*,*) 'shape matrix b  = ', bmat%nrows,bmat%mcols
         write(*,*) 'shape matrix c  = ', cmat%nrows,cmat%mcols
         write(*,*) 'shape matrix a  = ', shape(matar)
         write(*,*) 'shape matrix b  = ', shape(matbr)
         write(*,*) 'shape matrix c  = ', shape(matcr)
         call matmulcuda_r(matar,matbr,matcr,size(matar,1),size(matar,2),size(matbr,2))
         write(*,*) 'DONE'
         cmat%dmtx(1:size(matcr,1),1:size(matcr,2))=matcr
         deallocate(matcr)
        else
         if(.not.transpose_amat_)then
          allocate(matac(na1,na2))
          matac=amat%zmtx(1:na1,1:na2)
         else
          allocate(matac(na1,na2))
          matac=amat%zmtx(1:na1,1:na2)
          allocate(tmpc(na2,na1)); tmpc=transpose(matac)
          deallocate(matac);allocate(matac(na2,na1));
          matac=tmpc;deallocate(tmpc)
         endif
         if(.not.transpose_bmat_)then
          allocate(matbc(nb1,nb2))
          matbc=bmat%zmtx(1:nb1,1:nb2)
         else
          allocate(matbc(nb1,nb2))
          matbc=bmat%zmtx(1:nb1,1:nb2)
          allocate(tmpc(nb2,nb1)); tmpc=transpose(matbc)
          deallocate(matbc);allocate(matbc(nb2,nb1));
          matbc=tmpc;deallocate(tmpc)
         endif
         allocate(matcc(size(matac,1),size(matbc,2)))
         write(*,*) 'MATMUL GPU COMPLEX'
         write(*,*) 'shape matrix a  = ', amat%nrows,amat%mcols
         write(*,*) 'shape matrix b  = ', bmat%nrows,bmat%mcols
         write(*,*) 'shape matrix c  = ', cmat%nrows,cmat%mcols
         call matmulcuda_c(matac,matbc,matcc,size(matac,1),size(matac,2),size(matbc,2))
         write(*,*) 'DONE'
         cmat%zmtx(1:size(matcc,1),1:size(matcc,2))=matcc
         deallocate(matcc)
        endif
     endif
     if(cmat%iscmplx)then
      call comms_bcast(0,cmat%zmtx)
     else
      call comms_bcast(0,cmat%dmtx)
     endif
     return
     37 continue
#endif
!END CW

    ! Deal with optional argument for transposing matrices
    transa = 'n'
    if (present(transpose_amat)) then
       if (transpose_amat) then
          transa = 't'
       else
          transa = 'n'
       end if
    end if
    transb = 'n'
    if (present(transpose_bmat)) then
       if (transpose_bmat) then
          transb = 't'
       else
          transb = 'n'
       end if
    end if

    ! Check Matrix Sizes
    m = cmat%nrows
    n = cmat%mcols
    k = amat%mcols
    if (amat%nrows /= m) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_product: &
            &amat%nrows /= cmat%nrows'
       call comms_abort
    end if
    if (transb=='n') then
       if (bmat%mcols /= n) then
          if (pub_on_root) write(stdout,'(a)') 'Error in dense_product: &
               &bmat%mcols /= cmat%mcols'
          call comms_abort
       end if
       if (bmat%nrows /= k) then
          if (pub_on_root) write(stdout,'(a)') 'Error in dense_product: &
               &bmat%nrows /= amat%mcols'
          call comms_abort
       end if
    else ! transb=='t'
       if (bmat%nrows /= n) then
          if (pub_on_root) write(stdout,'(a)') 'Error in dense_product: &
               &bmat%nrows /= cmat%mcols'
          call comms_abort
       end if
       if (bmat%mcols /= k) then
          if (pub_on_root) write(stdout,'(a)') 'Error in dense_product: &
               &bmat%mcols /= amat%mcols'
          call comms_abort
       end if
    end if

    ! Deal with optional arguments for only multiplying certain row ranges
    if (present(first_k).and.present(last_k)) then
       ! Check range is valid
       if ((first_k < 0).or.(last_k > k).or.(last_k < first_k)) then
          if (pub_on_root) write(stdout,'(a)') 'Error in dense_product: &
               &invalid k range specification'
          call comms_abort
       end if
       ! Start from first_k'th column
       ia = 1
       ja = first_k
       if (transb=='n') then
          ib = first_k
          jb = 1
       else
          ib = 1
          jb = first_k
       end if
       ! Set range of k
       k = last_k - first_k + 1
    else
       ! Use whole range
       ia = 1
       ja = 1
       ib = 1
       jb = 1
    end if

    ! Perform matrix product with LAPACK (d/z)gemm or ScaLAPACK p(d/z)gemm
#ifdef SCALAPACK

    if (amat%iscmplx) then
       call PZGEMM(transa, transb, m, n, k, 1.0_DP, &
            amat%zmtx(1,1), ia, ja, amat%blacs_desc, &
            bmat%zmtx(1,1), ib, jb, bmat%blacs_desc, 0.0_DP, &
            cmat%zmtx(1,1), 1, 1, cmat%blacs_desc)
    else
       call PDGEMM(transa, transb, m, n, k, 1.0_DP, &
            amat%dmtx(1,1), ia, ja, amat%blacs_desc, &
            bmat%dmtx(1,1), ib, jb, bmat%blacs_desc, 0.0_DP, &
            cmat%dmtx(1,1), 1, 1, cmat%blacs_desc)
    end if

#else

    ! Semi-parallelised DGEMM call - nodes work out their columns of result
    batch_size = n/pub_total_num_nodes
    local_start = pub_my_node_id*batch_size + 1
    local_end = local_start + batch_size - 1
    if (pub_my_node_id == pub_total_num_nodes-1) local_end = n
    local_n = local_end - local_start + 1

    ! Sort out what local cols/rows of B to select, based on transb
    if (transb=='n') then
       jb = local_start
       ldb = bmat%nrows
    else
       ib = local_start
       ldb = n
    end if

    ! Perform the LAPACK call
    if (amat%iscmplx) then
       call ZGEMM(transa, transb, m, local_n, k, 1.0_DP, amat%zmtx(ia,ja), m, &
            bmat%zmtx(ib,jb), ldb, 0.0_DP, &
            cmat%zmtx(1,local_start), m)
    else
       call DGEMM(transa, transb, m, local_n, k, 1.0_DP, amat%dmtx(ia,ja), m, &
            bmat%dmtx(ib,jb), ldb, 0.0_DP, &
            cmat%dmtx(1,local_start), m)
    end if

    ! Broadcast results from each node to fill whole matrix of C
    do node=0,pub_total_num_nodes-1
       local_start = node*batch_size + 1
       local_end = local_start + batch_size - 1
       if (node == pub_total_num_nodes-1) local_end = n
       local_n = local_end - local_start + 1

       if (amat%iscmplx) then
          call comms_bcast(node,cmat%zmtx(1,local_start),m*local_n)
       else
          call comms_bcast(node,cmat%dmtx(1,local_start),m*local_n)
       end if
    end do

#endif

  end subroutine dense_product


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine finds the solution to the generalised eigenproblem defined !
  ! by two dense matrices amat and bmat, and returns the real eigenvalues and  !
  ! a dense matrix storing the eigenvectors.                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   eigenvals (out)   : Real array storing the eigenvalues of amat.          !
  !   amat (input)      : Dense Symmetric or Hermitian matrix.                 !
  !   bmat (input)      : Dense Symmetric or Hermitian matrix.                 !
  !   ibtype (input)    : Which type of eigenproblem to solve (see LAPACK help)!
  !   eigenvecs (out)   : Dense matrix storing eigenvectors defining solution. !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  ! Fixed to broadcast resulting eigenvectors from root node by Nicholas Hine  !
  ! on 10 March 2010, to fix problem of different eigenvectors on different    !
  ! nodes.                                                                     !
  !============================================================================!

  subroutine dense_eigensolve(n_eig,eigenvals,amat,bmat,ibtype,eigenvecs)

    use comms, only: comms_abort, comms_barrier, comms_bcast, pub_on_root, &
         pub_root_node_id, pub_total_num_nodes
!CW
    use comms, only : pub_my_node_id
!END CW
    use rundat, only: pub_devel_code
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: n_eig
    real(kind=DP), intent(out) :: eigenvals(:)
    type(DEM), intent(inout) :: eigenvecs
    type(DEM), intent(in) :: amat
    type(DEM), intent(in) :: bmat
    integer, intent(in) :: ibtype

    ! Local Variables
#ifdef SCALAPACK
    integer :: nzfound
    real(kind=DP) :: abstol, orfac
    character(len=200) :: dense_devel_code
    integer            :: start_pos, stop_pos, test_pos
#endif
    integer :: nfound
    integer :: ierr, info
    character :: jobz, range

    ! Work arrays and sizes
    integer :: lwork
#ifdef SCALAPACK
    integer :: liwork, lrwork
#endif
    integer, allocatable :: iwork(:), ifail(:)
    real(kind=DP), allocatable :: work(:)
    complex(kind=DP), allocatable :: zwork(:)
#ifdef SCALAPACK
    integer, allocatable :: iclustr(:)
    real(kind=DP), allocatable :: gap(:)
#endif

    ! Start timer
    call timer_clock('dense_eigensolve',1)

    ! Check arguments
    if (n_eig<=0) then
       jobz = 'N'
       range = 'A'
    else
       jobz = 'V'
       if ((n_eig>0).and.(n_eig<amat%nrows)) then
          range = 'I'
       else
          range = 'A'
       end if
    end if

    ! Allocate work arrays
    call internal_allocate_work

#ifdef SCALAPACK
    ! ndmh: call ScaLAPACK routine to solve the generalised eigenproblem
    abstol = 1e-9_DP
    ! ars: set orfac
    orfac = 1e-4_DP
    if (pub_on_root) then
       dense_devel_code=pub_devel_code
       if (len_trim(dense_devel_code)>0) then
          start_pos=index(dense_devel_code,'DENSE:')
          stop_pos=index(dense_devel_code,':DENSE')
          if (stop_pos<=0) stop_pos=len_trim(dense_devel_code) !missing end so go to end of string
          if (start_pos>0) then

             test_pos=index(dense_devel_code,'ORFAC=')
             if (test_pos>start_pos.and.test_pos<stop_pos) then
                test_pos=test_pos+len('ORFAC=')
                read(dense_devel_code(test_pos:test_pos+ &
                     & index(dense_devel_code(test_pos:stop_pos),':')-2),*) orfac
             end if
          end if
       end if
    end if
    call comms_bcast(pub_root_node_id,orfac)

    if (amat%iscmplx) then
       call PZHEGVX(ibtype, jobz, range, 'L', amat%nrows, &
            amat%zmtx(1,1), 1, 1, amat%blacs_desc, &
            bmat%zmtx(1,1), 1, 1, bmat%blacs_desc, &
            0.0_DP, 0.0_DP, 1, n_eig, abstol, nfound, nzfound, &
            eigenvals, orfac, eigenvecs%zmtx(1,1), 1, 1, eigenvecs%blacs_desc, &
            zwork, lwork, work, lrwork, iwork, liwork, ifail, iclustr, gap, &
            info)
    else
       call PDSYGVX(ibtype, jobz, range, 'L', amat%nrows, &
            amat%dmtx(1,1), 1, 1, amat%blacs_desc, &
            bmat%dmtx(1,1), 1, 1, bmat%blacs_desc, &
            0.0_DP, 0.0_DP, 1, n_eig, abstol, nfound, nzfound, &
            eigenvals, orfac, eigenvecs%dmtx(1,1), 1, 1, eigenvecs%blacs_desc, &
            work, lwork, iwork, liwork, ifail, iclustr, gap, info)
    end if
#else
    ! ndmh: call LAPACK routine to solve the generalised eigenproblem
    if (amat%iscmplx) then
!CW
#ifndef GPU_SPEEDUP_DENSE
       call ZHEGVX(ibtype,jobz,range,'L',amat%nrows, &
            amat%zmtx(1,1),amat%nrows,bmat%zmtx(1,1),bmat%nrows,&
            0.0_DP, 0.0_DP, 1, n_eig, -1.0_DP, nfound, eigenvals, &
            eigenvecs%zmtx(1,1),amat%nrows,zwork,lwork,work,iwork,ifail,info)
#else
       if(amat%nrows/=bmat%nrows)then
        write(*,*) 'NOT SUPPOSED TO HAPPEN a/=b rows (double complex)'
        stop
       endif 
       ifail=0;info=0;
#ifndef GPU_SPEEDUP_SINGLEPREC
       if(pub_my_node_id==0) &
   &    call eigenvector_cuda_type_c(ibtype,amat%nrows,amat%zmtx,bmat%zmtx,eigenvals,eigenvecs%zmtx,.false.)
#else
       if(pub_my_node_id==0) &
   &    call eigenvector_cuda_type_c(ibtype,amat%nrows,amat%zmtx,bmat%zmtx,eigenvals,eigenvecs%zmtx,.true.) 
#endif
        call comms_bcast(0,eigenvecs%zmtx);call comms_bcast(0,eigenvals);
#endif
!END CW
       ! ndmh: degenerate eigenvectors may have come out in different linear
       ! ndmh: combinations on different nodes, so synchronise to eigenvectors
       ! ndmh: found on root node
       call comms_bcast(pub_root_node_id,eigenvecs%zmtx, &
            eigenvecs%nrows*eigenvecs%mcols)
    else

!CW
#ifndef GPU_SPEEDUP_DENSE
       call DSYGVX(ibtype,jobz,range,'L',amat%nrows, &
            amat%dmtx(1,1),amat%nrows,bmat%dmtx(1,1),bmat%nrows, &
            0.0_DP, 0.0_DP, 1, n_eig, -1.0_DP, nfound, eigenvals, &
            eigenvecs%dmtx(1,1),amat%nrows,work,lwork,iwork,ifail,info)
#else
       if(amat%nrows/=bmat%nrows)then
        write(*,*) 'NOT SUPPOSED TO HAPPEN a/=b rows (double real)'
        stop
       endif 
       ifail=0;info=0;
#ifndef GPU_SPEEDUP_SINGLEPREC
       if(pub_my_node_id==0) &
   &   call eigenvector_cuda_type_r(ibtype,amat%nrows,amat%dmtx,bmat%dmtx,eigenvals,eigenvecs%dmtx,.false.)
#else
       if(pub_my_node_id==0) &
   &   call eigenvector_cuda_type_r(ibtype,amat%nrows,amat%dmtx,bmat%dmtx,eigenvals,eigenvecs%dmtx,.true.)
#endif
       if(pub_my_node_id==0) write(*,*) 'EIGENVAL GPU : ', eigenvals(1:10)
       call comms_bcast(0,eigenvecs%dmtx);call comms_bcast(0,eigenvals);
#endif
!END CW

       ! ndmh: degenerate eigenvectors may have come out in different linear
       ! ndmh: combinations on different nodes, so synchronise to eigenvectors
       ! ndmh: found on root node
       call comms_bcast(pub_root_node_id,eigenvecs%dmtx, &
            eigenvecs%nrows*eigenvecs%mcols)
    end if
#endif

    ! ndmh: check for errors
    if (info/=0) then
       if (amat%iscmplx) then
          if (pub_on_root) then
             write(stdout,'(a,i5)') '(P)ZHEGVX in subroutine &
                  &dense_eigensolve returned info=',info
             write(stdout,*) 'ifail=',ifail(1:info)
          end if
       else
          if (pub_on_root) then
             write(stdout,'(a,i5)') '(P)DSYGVX in subroutine &
                  &dense_eigensolve returned info=',info
             write(stdout,*) 'ifail=',ifail(1:info)
          end if
       end if
       ! ndmh: continue even if some eigenvectors did not converge
       if ((info<1).or.(info>amat%nrows)) then
          call comms_abort
       end if
    end if

!CW
#ifdef debug
    write(*,*) 'EIGENSOLVE DONE, ID = ', pub_my_node_id
    write(*,*) 'shape matrix a      = ', amat%nrows,amat%mcols
    write(*,*) 'shape matrix b      = ', bmat%nrows,bmat%mcols
#endif
!ENDCW

    ! Deallocate work arrays
    call internal_deallocate_work

    ! Stop timer
    call timer_clock('dense_eigensolve',2)

    return

contains

    ! Allocate diagonalisation work arrays
    subroutine internal_allocate_work

      implicit none

#ifdef SCALAPACK

      ! Local Variables
      integer :: n, nn, nb, np0, mq0, nq0, anb, sqnpc, nps
      integer :: nsytrd_lwopt, nsygst_lwopt
      integer :: liclustr, lgap
      integer, external :: numroc, pjlaenv, iceil

      ! Calculate optimal sizes of work arrays
      n = amat%nrows
      nb = amat%blacs_nb
      nn = max(amat%nrows,nb,2)
      np0 = numroc(nn,nb,0,0,blacs_nprow)
      nq0 = numroc(nn,nb,0,0,blacs_npcol)
      mq0 = numroc(nn,nb,0,0,blacs_npcol)
      if (amat%iscmplx) then
         anb = pjlaenv(blacs_context,3,'PZHETTRD','L',0,0,0,0)
      else
         anb = pjlaenv(blacs_context,3,'PDSYTTRD','L',0,0,0,0)
      end if
      sqnpc = int(sqrt(real(blacs_nprow*blacs_npcol,kind=DP)))
      nps = max(numroc(n,1,0,0,sqnpc),2*anb)
      nsytrd_lwopt = n + 2*(anb+1)*(4*nps+2) + (nps+3)*nps
      nsygst_lwopt = 2*np0*nb + nq0*nb + nb*nb
      lwork = 5*n + max(5*n, np0*nq0 + 2*nb*nb) + iceil(n, blacs_nprow*blacs_npcol)*nn
      lwork = max(lwork, 5*n + nsytrd_lwopt,nsygst_lwopt)
      liwork = 6*max(n,blacs_nprow*blacs_npcol+1,4)
      liclustr = 2*blacs_nprow*blacs_npcol
      lgap = blacs_nprow*blacs_npcol

      ! Allocate work arrays
      if (amat%iscmplx) then
         allocate(zwork(lwork),stat=ierr)
         call utils_alloc_check('dense_eigensolve','zwork',ierr)
      end if
      allocate(work(lwork),stat=ierr)
      call utils_alloc_check('dense_eigensolve','work',ierr)
      allocate(iwork(liwork),stat=ierr)
      call utils_alloc_check('dense_eigensolve','iwork',ierr)
      allocate(iclustr(liclustr),stat=ierr)
      call utils_alloc_check('dense_eigensolve','iclustr',ierr)
      allocate(gap(liclustr),stat=ierr)
      call utils_alloc_check('dense_eigensolve','gap',ierr)
      allocate(ifail(n),stat=ierr)
      call utils_alloc_check('dense_eigensolve','ifail',ierr)

#else

      ! Allocate work arrays
      lwork = 8*amat%nrows
      if (amat%iscmplx) then
         allocate(zwork(lwork),stat=ierr)
         call utils_alloc_check('dense_eigensolve','zwork',ierr)
      end if
      allocate(work(lwork),stat=ierr)
      call utils_alloc_check('dense_eigensolve','work',ierr)
      allocate(iwork(5*amat%nrows),stat=ierr)
      call utils_alloc_check('dense_eigensolve','iwork',ierr)
      allocate(ifail(amat%nrows),stat=ierr)
      call utils_alloc_check('dense_eigensolve','ifail',ierr)


#endif

    end subroutine internal_allocate_work

    ! Deallocate diagonalisation work arrays
    subroutine internal_deallocate_work

      implicit none

#ifdef SCALAPACK

      deallocate(ifail,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','ifail',ierr)
      deallocate(gap,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','gap',ierr)
      deallocate(iclustr,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','iclustr',ierr)
      deallocate(iwork,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','iwork',ierr)
      deallocate(work,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','work',ierr)
      if (amat%iscmplx) then
         deallocate(zwork,stat=ierr)
         call utils_dealloc_check('dense_eigensolve','zwork',ierr)
      end if

#else

      deallocate(ifail,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','ifail',ierr)
      deallocate(iwork,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','iwork',ierr)
      deallocate(work,stat=ierr)
      call utils_dealloc_check('dense_eigensolve','work',ierr)
      if (amat%iscmplx) then
         deallocate(zwork,stat=ierr)
         call utils_dealloc_check('dense_eigensolve','zwork',ierr)
      end if

#endif

    end subroutine internal_deallocate_work

  end subroutine dense_eigensolve


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine finds the inverse of a general dense matrix.               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat (input)      : Dense Symmetric or Hermitian matrix.                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 1 March 2010.                                    !
  !============================================================================!

  subroutine dense_invert(mat)

    use comms, only: comms_abort, comms_barrier, pub_on_root, &
         pub_total_num_nodes
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(DEM), intent(in) :: mat

    ! Local Variables
    integer, allocatable :: ipiv(:)
    real(kind=DP), allocatable :: work(:)
!CW
    complex(kind=DP), allocatable :: cwork(:)
!END CW

#ifdef SCALAPACK
    integer, allocatable :: iwork(:)
    integer :: liwork
#endif
    integer :: lwork
    integer :: ierr
    integer :: info

    ! Check Arguments
    if (mat%nrows /= mat%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'ERROR in dense_invert: &
            &Matrix is not square'
       call comms_abort
    end if

    ! Allocate work arrays
#ifdef SCALAPACK
    ! TODO FIX FIX FIX - uses too much memory
    liwork = 5*mat%nrows
    lwork = mat%nrows**2
    allocate(ipiv(mat%nrows),stat=ierr)
    call utils_alloc_check('dense_invert','ipiv',ierr)
    allocate(iwork(liwork),stat=ierr)
    call utils_alloc_check('dense_invert','iwork',ierr)
    allocate(work(lwork),stat=ierr)
    call utils_alloc_check('dense_invert','work',ierr)
#else
    lwork = 3*mat%nrows
    allocate(ipiv(mat%nrows),stat=ierr)
    call utils_alloc_check('dense_invert','ipiv',ierr)
    allocate(work(lwork),stat=ierr)
    call utils_alloc_check('dense_invert','work',ierr)
#endif

!CW
    allocate(cwork(lwork),stat=ierr)
    call utils_alloc_check('dense_invert','cwork_array',ierr)
!END CW

    ! cks: compute the factorization of a real symmetric matrix A using
    !      the Bunch-Kaufman diagonal pivoting method
#ifdef SCALAPACK
    if (mat%iscmplx) then
       ! ndmh: TODO not done yet
!CW
       write(*,*) 'PZETRF - SCALAPACK'
       call PZGETRF(mat%nrows, mat%nrows, mat%zmtx(1,1), 1, 1, mat%blacs_desc, &
            ipiv, info)
!END CW
    else
       call PDGETRF(mat%nrows, mat%nrows, mat%dmtx(1,1), 1, 1, mat%blacs_desc, &
            ipiv, info)
    end if
#else
    if (mat%iscmplx) then
       ! ndmh: TODO not done yet
!CW
       write(*,*) 'ZGETRF'
       call ZGETRF(mat%nrows, mat%nrows, mat%zmtx(1,1), mat%nrows, ipiv, info)
!END CW
    else
       call DGETRF(mat%nrows, mat%nrows, mat%dmtx(1,1), mat%nrows, ipiv, info)
    end if
#endif

    if (info/=0) then
!CW
       if (pub_on_root) write(stdout,'(a,i5)') 'ERROR in dense_invert: &
            &(P)(Z/D)GETRF returned info=',info
!END CW
!CW
       if (mat%iscmplx) write(*,*) 'matrix min/max : ',minval(abs(mat%zmtx)),maxval(abs(mat%zmtx)) 
!END CW
       call comms_abort
    endif

    ! cks: compute the inverse of a real symmetric indefinite matrix A using
    ! cks: the factorization A = U*D*U**T or A = L*D*L**T computed by DSYTRF
#ifdef SCALAPACK
    if (mat%iscmplx) then
       ! ndmh: TODO not done yet
!CW
       write(*,*) 'PZGETRi'
       call PZGETRI(mat%nrows, mat%zmtx(1,1), 1, 1, mat%blacs_desc, &
            ipiv, cwork, lwork, iwork, liwork, info)
!END CW
    else
       call PDGETRI(mat%nrows, mat%dmtx(1,1), 1, 1, mat%blacs_desc, &
            ipiv, work, lwork, iwork, liwork, info)
    end if
#else
    if (mat%iscmplx) then
       ! ndmh: TODO not done yet
!CW
       write(*,*) 'ZGETRI'
       call ZGETRI(mat%nrows, mat%zmtx(1,1), mat%nrows, ipiv, cwork, lwork, info)
!END CW
    else
       call DGETRI(mat%nrows, mat%dmtx(1,1), mat%nrows, ipiv, work, lwork, info)
    end if
#endif

    if (info/=0) then
!CW
       if (pub_on_root) write(stdout,'(a,i5)') 'ERROR in dense_invert: &
            &(P)(Z/D)GETRI returned info=',info
!END CW
       call comms_abort
    endif

    ! Deallocate work arrays
#ifdef SCALAPACK
    deallocate(work,stat=ierr)
    call utils_dealloc_check('dense_invert','work',ierr)
    deallocate(iwork,stat=ierr)
    call utils_dealloc_check('dense_invert','iwork',ierr)
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('dense_invert','ipiv',ierr)
#else
    deallocate(work,stat=ierr)
    call utils_dealloc_check('dense_invert','work',ierr)
    deallocate(ipiv,stat=ierr)
    call utils_dealloc_check('dense_invert','ipiv',ierr)
#endif

!CW
   deallocate(cwork,stat=ierr)
   call utils_dealloc_check('dense_invert','cwork',ierr)
!END CW

  end subroutine dense_invert


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !============================================================================!
  ! This subroutine converts a SPAM3 matrix to a dense matrix.                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dmat      (input) : The destination dense matrix                         !
  !   smat      (input) : The source sparse matrix                             !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  !============================================================================!

  subroutine dense_convert_sparsetodense(dmat_out,smat_in)

    use comms, only: comms_abort, pub_on_root
    use constants, only: stdout
#ifdef SCALAPACK
    use sparse, only: SPAM3,sparse_spam3toblacs
#else
    use sparse, only: SPAM3,sparse_convert
#endif

    implicit none

    ! Arguments
    type(DEM),intent(inout) :: dmat_out
    type(SPAM3),intent(in) :: smat_in

    ! Check compatability
    if (dmat_out%iscmplx.neqv.smat_in%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &dense_convert_densetosparse: matrices must be both real or both &
            &complex'
       call comms_abort
    end if

    ! Fill dense matrix values from SPAM3 matrix
#ifdef SCALAPACK
    if (smat_in%iscmplx) then
       call sparse_spam3toblacs(dmat_out%zmtx,smat_in,dmat_out%nrows, &
            dmat_out%blacs_ncol,dmat_out%blacs_desc)
    else
       call sparse_spam3toblacs(dmat_out%dmtx,smat_in,dmat_out%blacs_ld, &
            dmat_out%blacs_ncol,dmat_out%blacs_desc)
    end if
#else
    if (smat_in%iscmplx) then
       call sparse_convert(dmat_out%zmtx,smat_in)
    else
       call sparse_convert(dmat_out%dmtx,smat_in)
    end if
#endif

  end subroutine dense_convert_sparsetodense


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine converts a dense matrix to a SPAM3 matrix.                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dmat      (input) : The destination sparse matrix                        !
  !   smat      (input) : The source dense matrix                              !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 23 February 2010.                                !
  !============================================================================!

  subroutine dense_convert_densetosparse(smat_out,dmat_in)

    use comms, only: comms_abort, pub_on_root
    use constants, only: stdout
#ifdef SCALAPACK
    use sparse, only: SPAM3,sparse_blacstospam3
#else
    use sparse, only: SPAM3,sparse_convert
#endif

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: smat_out
    type(DEM),intent(in) :: dmat_in

    ! Check compatability
    if (smat_out%iscmplx.neqv.dmat_in%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &dense_convert_densetosparse: matrices must be both real or both &
            &complex'
       call comms_abort
    end if
    ! Fill SPAM3 matrix values from dense matrix
#ifdef SCALAPACK
    if (dmat_in%iscmplx) then
       call sparse_blacstospam3(smat_out,dmat_in%zmtx,dmat_in%nrows, &
            dmat_in%blacs_ncol,dmat_in%blacs_desc)
    else
       call sparse_blacstospam3(smat_out,dmat_in%dmtx,dmat_in%nrows, &
            dmat_in%blacs_ncol,dmat_in%blacs_desc)
    end if
#else
    if (dmat_in%iscmplx) then
       call sparse_convert(smat_out,dmat_in%zmtx)
    else
       call sparse_convert(smat_out,dmat_in%dmtx)
    end if
#endif

  end subroutine dense_convert_densetosparse


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine copies one dense matrix to another.                        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest      (input) : The destination dense matrix                         !
  !   dest      (input) : The source dense matrix                              !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_copy(dest,src)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: dest      ! The destination dense matrix
    type(DEM), intent(in) :: src          ! The source dense matrix

    ! Check Arguments
    if ((dest%nrows /= dest%nrows).or.(src%mcols /= src%mcols)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_copy: &
            &src and dest matrix sizes do not match'
       call comms_abort
    end if

    ! Trivial axpy of whole data arrays
    if (src%iscmplx) then
       if (dest%iscmplx) then
          dest%zmtx = src%zmtx
       else
          dest%dmtx(:,:) = real(src%zmtx(:,:),kind=DP)
       end if
    else
       if (dest%iscmplx) then
          dest%zmtx(:,:) = cmplx(src%dmtx(:,:),0.0_DP,kind=DP)
       else
          dest%dmtx = src%dmtx
       end if
    end if

  end subroutine dense_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine performs an axpy operation on the matrices:                !
  !   y := y + alpha * x                                                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The dense matrix y                                       !
  !   xmat  (input) : The dense matrix x                                       !
  !   alpha (input) : The real parameter alpha                                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_axpy_real(ymat,xmat,alpha)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: ymat      ! The dense matrix y
    type(DEM), intent(in) :: xmat         ! The dense matrix x
    real(kind=DP), intent(in) :: alpha    ! The parameter alpha

    ! Check Arguments
    if ((ymat%nrows /= xmat%nrows).or.(ymat%mcols /= xmat%mcols)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_axpy_real: &
            &xmat and ymat sizes do not match'
       call comms_abort
    end if

    ! Trivial axpy of whole data arrays
    if (xmat%iscmplx) then
       if (ymat%iscmplx) then
          ymat%zmtx = ymat%zmtx + alpha * xmat%zmtx
       else
          ymat%dmtx(:,:) = ymat%dmtx(:,:) + &
               alpha * real(xmat%zmtx(:,:),kind=DP)
       end if
    else
       if (ymat%iscmplx) then
          ymat%zmtx(:,:) = ymat%zmtx(:,:) + &
               alpha * cmplx(xmat%dmtx(:,:),0.0_DP,kind=DP)
       else
          ymat%dmtx = ymat%dmtx + alpha * xmat%dmtx
       end if
    end if

  end subroutine dense_axpy_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine performs an axpy operation on the matrices:                !
  !   y := y + alpha * x                                                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The dense matrix y                                       !
  !   xmat  (input) : The dense matrix x                                       !
  !   alpha (input) : The complex parameter alpha                              !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_axpy_complex(ymat,xmat,alpha)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: ymat      ! The dense matrix y
    type(DEM), intent(in) :: xmat         ! The dense matrix x
    complex(kind=DP), intent(in) :: alpha ! The parameter alpha

    ! Check Arguments
    if ((ymat%nrows /= xmat%nrows).or.(ymat%mcols /= xmat%mcols)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_axpy_complex: &
            &xmat and ymat sizes do not match'
       call comms_abort
    end if

    ! Trivial axpy of data arrays
    if (xmat%iscmplx) then
       if (ymat%iscmplx) then
          ymat%zmtx = ymat%zmtx + alpha * xmat%zmtx
       else
          ymat%dmtx(:,:) = ymat%dmtx(:,:) + &
               real(alpha * xmat%zmtx(:,:),kind=DP)
       end if
    else
       if (ymat%iscmplx) then
          ymat%zmtx(:,:) = ymat%zmtx(:,:) + &
               alpha * cmplx(xmat%dmtx(:,:),0.0_DP,kind=DP)
       else
          ymat%dmtx = ymat%dmtx + real(alpha * xmat%dmtx,kind=DP)
       end if
    end if

  end subroutine dense_axpy_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine performs an rank 1 update operation on the matrices:       !
  !   A := A + alpha * x * y^T                                                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat  (inout) : The dense matrix A                                       !
  !   xmat  (inout) : The dense matrix x                                       !
  !   ymat  (input) : The dense matrix y                                       !
  !   ix    (input) : The row of matrix x to use as the updating row           !
  !   jy    (input) : The row of matrix y to use as the updating row           !
  !   alpha (input) : The real parameter alpha                                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 6 February 2012.                                 !
  !============================================================================!

  subroutine dense_rank1_update_real(amat,xmat,ymat,ix,iy,alpha)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: amat      ! The dense matrix y
    type(DEM), intent(in) :: xmat         ! The dense matrix x
    type(DEM), intent(in) :: ymat         ! The dense matrix y
    integer, intent(in) :: ix,iy          ! Updating rows
    real(kind=DP), intent(in) :: alpha    ! The parameter alpha

    ! Local Variables
    complex(kind=DP) :: alpha_cmplx

    ! Check Arguments
    if ((ymat%nrows /= xmat%nrows).or.(ymat%mcols /= xmat%mcols)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_rank1_update_real: &
            &xmat and ymat sizes do not match'
       call comms_abort
    end if
    if ((amat%iscmplx.neqv.xmat%iscmplx) .or. &
         (amat%iscmplx.neqv.xmat%iscmplx)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_rank1_update_real: &
            &amat, xmat and ymat must be all real or all complex'
       call comms_abort
    end if

#ifdef SCALAPACK
    if (amat%iscmplx) then
       alpha_cmplx = cmplx(alpha,0.0_DP)
       call PZGERU(amat%nrows, amat%mcols, alpha_cmplx, &
            xmat%zmtx(1,1), 1,ix,xmat%blacs_desc, 1, &
            ymat%zmtx(1,1), 1,iy,ymat%blacs_desc, 1, &
            amat%zmtx(1,1), 1, 1,amat%blacs_desc)
    else
       call PDGER(amat%nrows, amat%mcols, alpha, &
            xmat%dmtx(1,1), 1,ix,xmat%blacs_desc, 1, &
            ymat%dmtx(1,1), 1,iy,ymat%blacs_desc, 1, &
            amat%dmtx(1,1), 1, 1,amat%blacs_desc)
    end if
#else
    if (amat%iscmplx) then
       alpha_cmplx = cmplx(alpha,0.0_DP)
       call ZGERU(amat%nrows, amat%mcols, alpha_cmplx, &
            xmat%zmtx(1,iy), 1, ymat%zmtx(1,iy), 1, &
            amat%zmtx(1,1), amat%nrows)
    else
       call DGER(amat%nrows, amat%mcols, alpha, &
            xmat%dmtx(1,iy), 1, ymat%dmtx(1,iy), 1, &
            amat%dmtx(1,1), amat%nrows)
    end if
#endif

  end subroutine dense_rank1_update_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine performs an rank 1 update operation on the matrices:       !
  !   A := A + alpha * x * y^T                                                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat  (inout) : The dense matrix A                                       !
  !   xmat  (inout) : The dense matrix x                                       !
  !   ymat  (input) : The dense matrix y                                       !
  !   ix    (input) : The row of matrix x to use as the updating row           !
  !   jy    (input) : The row of matrix y to use as the updating row           !
  !   alpha (input) : The complex parameter alpha                              !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 6 February 2012.                                 !
  !============================================================================!

  subroutine dense_rank1_update_complex(amat,xmat,ymat,alpha)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: amat      ! The dense matrix y
    type(DEM), intent(in) :: xmat         ! The dense matrix x
    type(DEM), intent(in) :: ymat         ! The dense matrix y
    complex(kind=DP), intent(in) :: alpha ! The parameter alpha

    ! Check Arguments
    if ((ymat%nrows /= xmat%nrows).or.(ymat%mcols /= xmat%mcols)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_rank1_update_real: &
            &xmat and ymat sizes do not match'
       call comms_abort
    end if
    if ((amat%iscmplx.neqv.xmat%iscmplx) .or. &
         (amat%iscmplx.neqv.xmat%iscmplx)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_rank1_update_real: &
            &amat, xmat and ymat must be all real or all complex'
       call comms_abort
    end if
    if (amat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &dense_rank1_update_complex: cannot update real amat with complex &
            alpha'
       call comms_abort
    end if

    if (pub_on_root) write(stdout,'(a)') 'Error in &
         &dense_rank1_update_complex: routine not implemented yet. Contact &
         &Nicholas Hine if you need it'
    call comms_abort

#ifdef SCALAPACK
    
#else
    
#endif

  end subroutine dense_rank1_update_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues of a matrix.  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (inout) : The matrix to be rescaled and shifted                   !
  !   alpha  (input) : The real scaling parameter                              !
  !   beta   (input) : The optional real shift parameter                       !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, March 2010 based on SPAM3 version.               !
  !============================================================================!

  subroutine dense_scale_real(mat,alpha,beta)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    type(DEM), intent(inout)            :: mat    ! The matrix to be operated on
    real(kind=DP), intent(in)           :: alpha  ! The scaling parameter
    real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

    ! Local Variables
    integer :: ielem

    ! Rescale the elements if required
    if (alpha /= 1.0_DP) then
       if (mat%iscmplx) then
          mat%zmtx = alpha * mat%zmtx
       else
          mat%dmtx = alpha * mat%dmtx
       end if
    end if

    ! Shift the eigenvalues if required
    if (present(beta)) then

       ! Can only shift eigenvalues of square matrix
       if (mat%mcols/=mat%nrows) then
          if (pub_on_root) write(stdout,'(a)') 'Error in dense_scale_real: &
               &cannot shift eigenvalues of non-square matrices'
          call comms_abort
       end if

       ! Add identity matrix scaled by beta to all elements
       if (mat%iscmplx) then
          do ielem=1,mat%nrows
             mat%zmtx(ielem,ielem) = mat%zmtx(ielem,ielem) + &
                  cmplx(beta,0.0_DP,kind=DP)
          end do
       else
          do ielem=1,mat%nrows
             mat%dmtx(ielem,ielem) = mat%dmtx(ielem,ielem) + beta
          end do
       end if
    end if

  end subroutine dense_scale_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues of a matrix.  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (inout) : The matrix to be rescaled and shifted                   !
  !   alpha  (input) : The real scaling parameter                              !
  !   beta   (input) : The optional real shift parameter                       !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, March 2010 based on SPAM3 version.               !
  !============================================================================!

  subroutine dense_scale_complex(mat,alpha,beta)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    type(DEM), intent(inout)            :: mat    ! The matrix to be operated on
    complex(kind=DP), intent(in)        :: alpha  ! The scaling parameter
    real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

    ! Local Variables
    integer :: ielem

    ! This only makes sense if mat is complex...
    if (.not.mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_scale_complex: &
            &mat must be complex.'
       call comms_abort
    end if

    ! Rescale the elements if required
    if (alpha /= 1.0_DP) then
       if (mat%iscmplx) then
          mat%zmtx = alpha * mat%zmtx
       end if
    end if

    ! Shift the eigenvalues if required
    if (present(beta)) then

       ! Can only shift eigenvalues of square matrix
       if (mat%mcols/=mat%nrows) then
          if (pub_on_root) write(stdout,'(a)') 'Error in dense_scale_real: &
               &cannot shift eigenvalues of non-square matrices'
          call comms_abort
       end if

       ! Add identity matrix scaled by beta to all elements
       if (mat%iscmplx) then
          do ielem=1,mat%nrows
             mat%zmtx(ielem,ielem) = mat%zmtx(ielem,ielem) + &
                  cmplx(beta,0.0_DP,kind=DP)
          end do
       end if
    end if

  end subroutine dense_scale_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine retrieves an element from a given dense matrix, where the  !
  ! element can be local to any processor.                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (output)    : The element of the real matrix retrieved                !
  !   mat  (input)   : The real dense matrix                                   !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_get_element_real(el,mat,jrow,icol)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: el          ! The element to insert
    type(DEM), intent(in) :: mat              ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! Check matrix is real
    if (mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_get_element_real: &
            &mat must be real.'
       call comms_abort
    end if

    ! Get element directly from matrix, or call pdelget
#ifdef SCALAPACK
    call pdelget('A',' ',el,mat%dmtx,jrow,icol,mat%blacs_desc)
#else
    el = mat%dmtx(jrow,icol)
#endif

  end subroutine dense_get_element_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine inserts an element into a given dense matrix, where the    !
  ! element can be local to any processor.                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (output)    : The element to insert into the real matrix              !
  !   mat  (input)   : The real dense matrix                                   !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_put_element_real(el,mat,jrow,icol)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: el           ! The element to insert
    type(DEM), intent(inout) :: mat           ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! Check matrix is real
    if (mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_put_element_real: &
            &mat must be real.'
       call comms_abort
    end if

    ! Put element directly in matrix, or call pdelset
#ifdef SCALAPACK
    call pdelset(mat%dmtx,jrow,icol,mat%blacs_desc,el)
#else
    mat%dmtx(jrow,icol) = el
#endif

  end subroutine dense_put_element_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine retrieves an element from a given dense matrix, where the  !
  ! element can be local to any processor.                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (output)    : The element of the complex matrix retrieved             !
  !   mat  (input)   : The complex dense matrix                                !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_get_element_complex(el,mat,jrow,icol)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: el       ! The element to insert
    type(DEM), intent(in) :: mat              ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! Check matrix is complex
    if (.not.mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_get_element_complex: &
            &mat must be complex.'
       call comms_abort
    end if

    ! Get element directly from matrix, or call pdelget
#ifdef SCALAPACK
    call pzelget('A',' ',el,mat%zmtx,jrow,icol,mat%blacs_desc)
#else
    el = mat%zmtx(jrow,icol)
#endif

  end subroutine dense_get_element_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine inserts an element into a given dense matrix, where the    !
  ! element can be local to any processor.                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (output)    : The element to insert into the complex matrix           !
  !   mat  (input)   : The complex dense matrix                                !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_put_element_complex(el,mat,jrow,icol)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    complex(kind=DP), intent(in) :: el        ! The element to insert
    type(DEM), intent(inout) :: mat           ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! Check matrix is complex
    if (.not.mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_put_element_complex: &
            &mat must be complex.'
       call comms_abort
    end if

    ! Put element directly in matrix, or call pzelset
#ifdef SCALAPACK
    call pzelset(mat%zmtx,jrow,icol,mat%blacs_desc,el)
#else
    mat%zmtx(jrow,icol) = el
#endif

  end subroutine dense_put_element_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine retrieves a column from a given dense matrix.              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   col (output)   : The column of the real matrix retrieved                 !
  !   mat  (input)   : The real dense matrix                                   !
  !   icol (input)   : The column index of the column required                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_get_col_real(col,mat,icol)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: col(:)      ! The column to insert
    type(DEM), intent(in) :: mat              ! The matrix
    integer, intent(in) :: icol               ! The index of the column

    ! Local variables
    integer :: jrow

    ! Check matrix is real
    if (mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_get_col_real: &
            &mat must be real.'
       call comms_abort
    end if
    ! Check size of col
    if (size(col)/=mat%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_get_col_real: &
            &incorrect size of input column.'
       call comms_abort
    end if

    ! Get column directly from matrix, or call pdelget
    do jrow=1,mat%nrows
#ifdef SCALAPACK
       call pdelget('A',' ',col(jrow),mat%dmtx,jrow,icol,mat%blacs_desc)
#else
       col(jrow) = mat%dmtx(jrow,icol)
#endif
    end do

  end subroutine dense_get_col_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine inserts a column into a given dense matrix.                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   col (output)   : The column to insert into the real matrix               !
  !   mat  (input)   : The real dense matrix                                   !
  !   icol (input)   : The column index of the column required                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_put_col_real(col,mat,icol)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: col(:)       ! The column to insert
    type(DEM), intent(inout) :: mat           ! The matrix
    integer, intent(in) :: icol               ! The index of the column

    ! Local variables
    integer :: jrow

    ! Check matrix is real
    if (mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_put_col_real: &
            &mat must be real.'
       call comms_abort
    end if
    ! Check size of col
    if (size(col)/=mat%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_put_col_real: &
            &incorrect size of input column.'
       call comms_abort
    end if

    ! Put column directly in matrix, or call pdelset
    do jrow=1,mat%nrows
#ifdef SCALAPACK
       call pdelset(mat%dmtx,jrow,icol,mat%blacs_desc,col(jrow))
#else
       mat%dmtx(jrow,icol) = col(jrow)
#endif
    end do

  end subroutine dense_put_col_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !============================================================================!
  ! This subroutine retrieves an column from a given dense matrix, where the  !
  ! column can be local to any processor.                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   col (output)   : The column of the complex matrix retrieved             !
  !   mat  (input)   : The complex dense matrix                                !
  !   icol (input)   : The column index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 February 2010.                                !
  !============================================================================!

  subroutine dense_get_col_complex(col,mat,icol)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: col(:)   ! The column to insert
    type(DEM), intent(in) :: mat              ! The matrix
    integer, intent(in) :: icol               ! The index of the column

    ! Local variables
    integer :: jrow

    ! Check matrix is complex
    if (.not.mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_get_col_complex: &
            &mat must be complex.'
       call comms_abort
    end if
    ! Check size of col
    if (size(col)/=mat%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_get_col_complex: &
            &incorrect size of input column.'
       call comms_abort
    end if

    ! Get column directly from matrix, or call pdelget
    do jrow=1,mat%nrows
#ifdef SCALAPACK
       call pzelget('A',' ',col(jrow),mat%zmtx,jrow,icol,mat%blacs_desc)
#else
       col(jrow) = mat%zmtx(jrow,icol)
#endif
    end do

  end subroutine dense_get_col_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !============================================================================!
  ! This subroutine inserts an column into a given dense matrix.               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   col  (input)   : The column to insert into the complex matrix            !
  !   mat  (inout)   : The complex dense matrix                                !
  !   jrow (input)   : The column index of the row required                    !
  !   icol (input)   : The column index of the column required                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, 24 January 2011.                                 !
  !============================================================================!

  subroutine dense_put_col_complex(col,mat,icol)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    complex(kind=DP), intent(in) :: col(:)    ! The column to insert
    type(DEM), intent(inout) :: mat           ! The matrix
    integer, intent(in) :: icol               ! The index of the column

    ! Local variables
    integer :: jrow

    ! Check matrix is complex
    if (.not.mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_put_col_complex: &
            &mat must be complex.'
       call comms_abort
    end if
    ! Check size of col
    if (size(col)/=mat%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in dense_put_col_complex: &
            &incorrect size of input column.'
       call comms_abort
    end if

    ! Put column directly in matrix, or call pzelset
    do jrow=1,mat%nrows
#ifdef SCALAPACK
       call pzelset(mat%zmtx,jrow,icol,mat%blacs_desc,col(jrow))
#else
       mat%zmtx(jrow,icol) = col(jrow)
#endif
    end do

  end subroutine dense_put_col_complex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module dense
