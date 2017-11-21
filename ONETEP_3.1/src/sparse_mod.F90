! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-

module sparse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   The module in this file was originally created by Nicholas Hine in May    !
!   2009, based largely on the existing SPAM2 code, written by Peter Haynes   !
!   and modified by Nicholas Hine, between 2004 and 2009.                     !
!                                                                             !
!   TCM Group, Cavendish laboratory, University of Cambridge                  !
!   Madingley Road, Cambridge CB3 0HE, UK                                     !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=============================================================================!
!                                                                             !
! Block sparse matrix module - documentation                                  !
!                                                                             !
!-----------------------------------------------------------------------------!
!                                                                             !
! The SPAM3 type and its associated module are an update to the SPAM2 type on !
! which ONETEP v2.0-2.3 relied, to improve both the performance and the       !
! capabilities of the matrix algebra. Most significantly, SPAM3 allows        !
! matrices of arbitrary (non-square) dimensions, for applications such as     !
! < NGWF | projector > type overlaps, < NGWF | Hubbard U projector >          !
! matrices, and many others. The new code also represents a redesign          !
! of the data structures for even more efficient parallelisation.             !
! Whereas in SPAM2, matrices were always either wholly-sparse or wholly-dense !
! SPAM3 divides the matrix into "segments", which are blocks corresponding to !
! the sections of the columns belonging to a given node associated with rows  !
! of a second node. Segments are assigned sparse or dense depending on the    !
! number of nonzero elements in that segment.                                 !
!                                                                             !
! Sparse matrices in ONETEP are divided column-wise over the different nodes  !
! on which the simulation is running, attempting to share out the NGWFs       !
! equally over the nodes. A given node has a certain set of atoms (and thus   !
! NGWFs) associated with it, and only stores the data corresponding to those  !
! NGWF columns. When representing non-square matrices such as an              !
! NGWF-Projector overlap matrix, the same atom-blocking scheme is retained,   !
! with the sizes of the blocks modified accordingly.                          !
!                                                                             !
! The block sparse matrices in ONETEP arise because the structure of the      !
! sparse matrices containing matrix elements between NGWFs (the localised     !
! functions) depends upon whether the spherical atom-centred regions overlap. !
! The matrix elements are therefore grouped into blocks whose size is         !
! determined by the number of NGWFs on each atom.                             !
!                                                                             !
! Imagine a sparse matrix with the structure the overlap matrix of a          !
! simplified linear butane molecule, with 1 NGWF on each hydrogen atom and 4  !
! on each carbon atom                                                         !
!                                                                             !
!            H3   H6   H9   H12                                               !
!            |    |    |    |                                                 !
!       H1 - C2 - C5 - C8 - C11- H14                                          !
!            |    |    |    |                                                 !
!            H4   H7   H10  H13                                               !
!                                                                             !
! Assume the NGWF radii are small enough that only atoms associated with      !
! neighbouring carbon atoms overlap (for illustrative purposes only).         !
! Imagine simulating this on 4 nodes of a parallel computer, numbered 0,1,2,3,!
! with 4,3,3, and 4 atoms respectively on them. The sparsity pattern of the   !
! resulting matrix would look something like this, with X's representing      !
! nonzero elements of the overlap matrix:                                     !
!                                                                             !
!  |--------------------||------------------------------------------||        !
!  |                Node|| 0        || 1      || 2      || 3        ||        !
!  |--------------------||----------++--------++--------++----------||        !
!  |    |    |     |Atom||1| 2  |3|4|| 5  |6|7|| 8  |910|| 11 121314||        !
!  |    |    |     |Type||H|CCCC|H|H||CCCC|H|H||CCCC|H|H||CCCC|H|H|H||        !
!  |Node|Atom| Type|NGWF||1|1234|1|1||1234|1|1||1234|1|1||1234|1|1|1||        !
!  |--------------------||----------++--------++--------++----------||        !
!  |    | 1  | H   | 1  ||X|XXXX|X|X||OOOO|O|O||OOOO|O|O||OOOO|O|O|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    |    |     | 1  ||X|XXXX|X|X||XXXX|X|X||OOOO|O|O||OOOO|O|O|O||        !
!  |    | 2  | C   | 2  ||X|XXXX|X|X||XXXX|X|X||OOOO|O|O||OOOO|O|O|O||        !
!  | 0  |    |     | 3  ||X|XXXX|X|X||XXXX|X|X||OOOO|O|O||OOOO|O|O|O||        !
!  |    |    |     | 4  ||X|XXXX|X|X||XXXX|X|X||OOOO|O|O||OOOO|O|O|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 3  | H   | 1  ||X|XXXX|X|X||XXXX|X|X||OOOO|O|O||OOOO|O|O|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 4  | H   | 1  ||X|XXXX|X|X||XXXX|X|X||OOOO|O|O||OOOO|O|O|O||        !
!  |--------------------||==========++========++========++==========||        !
!  |    |    |     | 1  ||O|XXXX|X|X||XXXX|X|X||XXXX|X|X||OOOO|O|O|O||        !
!  |    | 5  | C   | 2  ||O|XXXX|X|X||XXXX|X|X||XXXX|X|X||OOOO|O|O|O||        !
!  | 1  |    |     | 3  ||O|XXXX|X|X||XXXX|X|X||XXXX|X|X||OOOO|O|O|O||        !
!  |    |    |     | 4  ||O|XXXX|X|X||XXXX|X|X||XXXX|X|X||OOOO|O|O|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 6  | H   | 1  ||O|XXXX|X|X||XXXX|X|X||XXXX|X|X||OOOO|O|O|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 7  | H   | 1  ||O|XXXX|X|X||XXXX|X|X||XXXX|X|X||OOOO|O|O|O||        !
!  |--------------------||==========++========++========++==========||        !
!  |    |    |     | 1  ||O|OOOO|O|O||XXXX|X|X||XXXX|X|X||XXXX|X|X|O||        !
!  |    | 8  | C   | 2  ||O|OOOO|O|O||XXXX|X|X||XXXX|X|X||XXXX|X|X|O||        !
!  |    |    |     | 3  ||O|OOOO|O|O||XXXX|X|X||XXXX|X|X||XXXX|X|X|O||        !
!  | 2  |    |     | 4  ||O|OOOO|O|O||XXXX|X|X||XXXX|X|X||XXXX|X|X|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 9  | H   | 1  ||O|OOOO|O|O||XXXX|X|X||XXXX|X|X||XXXX|X|X|O||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 10 | H   | 1  ||O|OOOO|O|O||XXXX|X|X||XXXX|X|X||XXXX|X|X|O||        !
!  |--------------------||==========++========++========++==========||        !
!  |    |    |     | 1  ||O|OOOO|O|O||OOOO|O|O||XXXX|X|X||XXXX|X|X|X||        !
!  |    | 11 | C   | 2  ||O|OOOO|O|O||OOOO|O|O||XXXX|X|X||XXXX|X|X|X||        !
!  |    |    |     | 3  ||O|OOOO|O|O||OOOO|O|O||XXXX|X|X||XXXX|X|X|X||        !
!  | 3  |    |     | 4  ||O|OOOO|O|O||OOOO|O|O||XXXX|X|X||XXXX|X|X|X||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 12 | H   | 1  ||O|OOOO|O|O||OOOO|O|O||XXXX|X|X||XXXX|X|X|X||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 13 | H   | 1  ||O|OOOO|O|O||OOOO|O|O||XXXX|X|X||XXXX|X|X|X||        !
!  |    |    |     |    ||----------++--------++--------++----------||        !
!  |    | 14 | H   | 1  ||O|OOOO|O|O||OOOO|O|O||OOOO|O|O||XXXX|X|X|X||        !
!  |--------------------||----------++--------++--------++----------||        !
!  |--------------------||------------------------------------------||        !
!                                                                             !
! We can see that, for example, the columns stored on node 0 have no nonzero  !
! elements in rows associated with nodes 2 and 3, whereas it has some but not !
! all nonzero elements corresponding to rows associated with node 1 and every !
! nonzero element corresponding to rows stored on row 0. Putting this into a  !
! matrix of size Nnodes x Nnodes, where 'D' represents a segment with dense   !
! storage, 'S' represents a segment with sparse storage and 'B' is a blank    !
! segment, we have                                                            !
!                                                                             !
! ++-+-------++                                                               !
! ||N|0|1|2|3||                                                               !
! ++-+-+-+-+-++                                                               !
! ||0|D|S|B|B||                                                               !
! ||1|S|D|D|B||                                                               !
! ||2|B|D|D|S||                                                               !
! ||3|B|B|S|D||                                                               !
! +++--------++                                                               !
!                                                                             !
! For the sparse segments, the segment is split into 'blocks', which are      !
! rectangular sections of the matrix corresponding to the overlap of all the  !
! elements associated with one atom with all the elements associated with     !
! another (or itself). The sizes of the blocks are determined by the number   !
! of elements on each atom, be they NGWFs, projectors, etc.                   !
!                                                                             !
! Matrices arising in ONETEP (in the limit of large system sizes) will contain!
! many blocks of zeroes, and therefore it is more efficient to index atomic   !
! blocks rather than individual matrix elements. However, if all the blocks   !
! within a given segment are nonzero, then the indexing is unnecessary and    !
! this part of the matrix may should be stored as a dense segment. Equally, if!
! there are no blocks within a given segment, there is no need to index them. !
!                                                                             !
!-----------------------------------------------------------------------------!
!                                                                             !
! Whether in sparse or dense format, the elements are stored in               !
! column-indexed format. This is because in Fortran, the first array index    !
! changes fastest as one steps through memory i.e. the order in which         !
! elements are stored for an array a(n,m) is:                                 !
!                                                                             !
! a(1,1) | a(2,1) | a(3,1) | ... | a(n,1) | a(1,2) | a(2,2) | ... | a(n,2) ...!
!                                                                             !
! The first index conventionally refers to the rows, and the second to the    !
! columns, of a matrix i.e. a matrix A may be represented by an array of      !
! dimension 2 where the elements are stored as:                               !
!                                                                             !
!    A    <->  a(i,j)                                                         !
!     ij                                                                      !
!                                                                             !
! Thus stepping through memory corresponds to moving down a column.           !
!                                                                             !
! For efficiency, one wishes to address memory sequentially if possible. So   !
! for a simple operation such as axpy one would code as follows:              !
!                                                                             !
!    A    := A   + alpha B                                                    !
!     ij      ij          ij                                                  !
!                                                                             !
!    do j=1,m                                ! loop over columns              !
!       do i=1,n                             ! loop over rows in column j     !
!          a(i,j) = a(i,j) + alpha * b(i,j)  ! element-wise axpy              !
!       end do                               ! loop over rows in column j     !
!    end do                                  ! loop over columns              !
!                                                                             !
! Matrix-matrix multiplication is more complicated. In Fortran maths libraries!
! matrices are represented by arrays in this column-wise fashion.             !
!                                                                             !
!-----------------------------------------------------------------------------!
!                                                                             !
! Sparse segments are represented using an atom-blocked column-indexed sparse !
! storage method. Two integer arrays and one double precision real array are  !
! used to store the matrices as follows:                                      !
!                                                                             !
! Let the number of block-rows and block-columns (i.e. atoms) be nblk, the    !
! number of segments with any nonzero blocks in be nseg, and the number of    !
! nonzero blocks to be stored be nzb.                                         !
!                                                                             !
!  integer :: seg_idx(0:nnodes)          ! Starting indices for each segment  !
!  integer :: blk_idx((nblk+1)*nseg+nzb) ! an index of nonzero blocks         !
!  integer :: seg_ptr(0:nnodes)          ! Pointers to each segment's data    !
!  integer :: blk_ptr((nblk+1)*nseg+nzb) ! pointers to the data in the blocks !
!  real(kind=DP) :: dmtx(nze)        ! the data (matrix elements) themselves  !
!                                                                             !
! The index array works as follows. The segment index points to the starting  !
! indices of each segment, and then the first nblk+1 entries for each         !
! index segment describe the index itself i.e. they point to a list of        !
! nonzero block-rows in each block-column for that segment.                   !
! This list begins at blk_idx(nblk+2) and there is one entry in this list for !
! each of the nzb nonzero blocks in the matrix. The list entries themselves   !
! are the indexes of the nonzero block-rows in a given block-column. For the  !
! carbon dioxide example above we have nnodes=4 (0-3), and on node 0, we have !
! nblk=4 and nzb=25, with the '0' segment having 16 nonzero blocks and the    !
! '1' segment having 9. Hence, the segment index on node 0 would be:          !
!                                                                             !
! seg_idx(0)  =                                 1                             !
! seg_idx(1)  = seg_idx(0) + nblk+1 + nzb(0) = 22                             !
! seg_idx(2)  = seg_idx(1) + nblk+1 + nzb(1) = 36                             !
! seg_idx(3)  = seg_idx(2) (as 2 is blank)   = 36                             !
! seg_idx(4)  = seg_idx(3) (as 3 is blank)   = 36                             !
!                                                                             !
!                                                                             !
! blk_idx(1)  = nblk+2 = 6      ! start of list of nonzero block-rows in col 1!
! blk_idx(2)  = 10              ! start of list of nonzero block-rows in col 2!
! blk_idx(3)  = 14              ! start of list of nonzero block-rows in col 3!
! blk_idx(4)  = 18              ! start of list of nonzero block-rows in col 4!
! blk_idx(5)  = nblk+nzb+2 = 22 ! end+1 of list of nonzero block-rows in col 4!
! blk_idx(6)  = 1               ! index of first nonzero block-row in col 1   !
! blk_idx(7)  = 2               ! index of second nonzero block-row in col 1  !
! blk_idx(8)  = 3               ! index of third nonzero block-row in col 1   !
! blk_idx(9)  = 4               ! index of fourth nonzero block-row in col 1  !
! blk_idx(10) = 1               ! index of first nonzero block-row in col 2   !
! blk_idx(11) = 2               ! index of second nonzero block-row in col 2  !
! ...                                                                         !
! blk_idx(21) = 4               ! index of fourth nonzero block-row in col 4  !
! blk_idx(22) = 27              ! start of list of nonzero block-rows in col 1!
! blk_idx(23) = 27              ! start of list of nonzero block-rows in col 2!
! blk_idx(24) = 30              ! start of list of nonzero block-rows in col 3!
! blk_idx(25) = 33              ! start of list of nonzero block-rows in col 4!
! blk_idx(26) = 36              ! end+1 of list of nonzero block-rows in col 4!
! blk_idx(27) = 5               ! index of first nonzero block-row in col 2   !
! blk_idx(28) = 6               ! index of second nonzero block-row in col 2  !
! blk_idx(29) = 7               ! index of third nonzero block-row in col 2   !
! blk_idx(30) = 5               ! index of first nonzero block-row in col 3   !
! blk_idx(31) = 6               ! index of second nonzero block-row in col 3  !
! ...                                                                         !
! blk_idx(35) = 7               ! index of third nonzero block-row in col 4   !
!                                                                             !
! Thus one can straightforwardly loop over all nonzero blocks in a given      !
! column of the whole matrix as follows:                                      !
!                                                                             !
! do seg=0,nnodes-1                                                           !
!   seg_start = seg_idx(seg)                                                  !
!   do blk_col=1,nblk                          ! loop over all columns on node!
!      do idx=blk_idx(seg_start+blk_col-1), &  ! loop along list of nonzero   !
!             blk_idx(seg_start+blk_col)-1     ! block-rows in block-column   !
!                                              ! blk_col for this segment     !
!         blk_row = blk_idx(idx)               ! index of block-row           !
!         ...                                                                 !
!                                                                             !
! Thus the number of non-zero block-rows in block-column blk_col in segment   !
! seg is given by:                                                            !
!    blk_idx(seg_idx(seg)+blk_col) - blk_idx(seg_idx(seg)+blk_col-1)          !
! and this is why blk_idx(seg_idx(seg)+nblk) is set to nblk+nzb+2.            !
!                                                                             !
! The matrix elements in each nonzero block are stored sequentially as        !
! matrices in the data array in the same order in memory as they would be if  !
! stored as dense two-dimensional arrays **for each segment** as in section 1 !
! above. The blk_ptr points to the start of each block's matrix elements      !
! in the data array. The first nblk entries of blk_ptr contain additional     !
! pointers to the diagonal blocks. The entry blk_ptr(nblk+1) is unused. There !
! is also an extra pointer at the end to mark the end of the data. For the    !
! above example this works out as follows:                                    !
!                                                                             !
! blk_ptr(1)  = 1               ! pointer to diagonal block (1,1)             !
! blk_ptr(2)  = 9               ! pointer to diagonal block (2,2)             !
! blk_ptr(3)  = 41              ! pointer to diagonal block (3,3)             !
! blk_ptr(4)  = 49              ! pointer to diagonal block (4,4)             !
! blk_ptr(5)  = 0               ! unused                                      !
! blk_ptr(6)  = 1               ! pointer to 1x1 array for block (1,1)        !
! blk_ptr(7)  = 2               ! pointer to 1x4 array for block (2,1)        !
! blk_ptr(8)  = 6               ! pointer to 1x1 array for block (3,1)        !
! blk_ptr(9)  = 7               ! pointer to 1x1 array for block (4,1)        !
! blk_ptr(10) = 8               ! pointer to 4x1 array for block (1,2)        !
! blk_ptr(11) = 12              ! pointer to 4x4 array for block (2,2)        !
! ...                                                                         !
! blk_ptr(21) = 50              ! pointer to end+1 of data for block (4,4)    !
! blk_ptr(22) = 0               ! unused                                      !
! blk_ptr(23) = 0               ! unused                                      !
! blk_ptr(24) = 0               ! unused                                      !
! blk_ptr(25) = 0               ! unused                                      !
! blk_ptr(26) = 50              ! pointer to 4x4 array for block (5,2)        !
! blk_ptr(27) = 56              ! pointer to 4x1 array for block (6,2)        !
! blk_ptr(28) = 60              ! pointer to 4x1 array for block (7,2)        !
! blk_ptr(30) = 64              ! pointer to 1x4 array for block (5,3)        !
! blk_ptr(31) = 68              ! pointer to 1x1 array for block (6,3)        !
! ...                                                                         !
! blk_ptr(35) = 86              ! pointer to end+1 of data array for this node!
!                                                                             !
! The columns of the  sparse matrix on this node thus contain 85 nonzero      !
! elements and this is the length of the dmtx array nze on this node.         !
!                                                                             !
! Note that on nodes beyond node zero, the list of diagonal blocks would      !
! appear later, in the section of the array corresponding to that segment of  !
! that node.                                                                  !
!                                                                             !
!-----------------------------------------------------------------------------!


  use constants, only: DP, LONG, stdout

  implicit none

  private

  ! Public subroutines and functions

  ! sparse_init routines
  public :: sparse_count_ss
  public :: sparse_index_ss
  public :: sparse_count_union
  public :: sparse_index_union

  ! sparse_base routines
  public :: sparse_create
  public :: sparse_destroy
  public :: sparse_index_length
  public :: sparse_generate_index
  public :: sparse_transpose_structure
  public :: sparse_get_element
  public :: sparse_put_element
  public :: sparse_get_block
  public :: sparse_put_block
  public :: sparse_get_col
  public :: sparse_put_col
  public :: sparse_clr_col

  ! sparse_inquiry routines
  public :: sparse_rms_element
  public :: sparse_max_abs_element
  public :: sparse_element_exists
  public :: sparse_fill_fac_denom
  public :: sparse_is_dense
  public :: sparse_num_element
  public :: sparse_node_num_element
  public :: sparse_num_rows
  public :: sparse_num_cols
  public :: sparse_first_elem_on_node
  public :: sparse_first_elem_on_atom
  public :: sparse_num_elems_on_node
  public :: sparse_num_elems_on_atom
  public :: sparse_atom_of_elem
  public :: sparse_node_of_elem
  public :: sparse_memory
  public :: sparse_show_memory_usage

  ! sparse_ops routines
  public :: sparse_copy
  public :: sparse_scale
  public :: sparse_axpy
  public :: sparse_product
  public :: sparse_trace
  public :: sparse_transpose
  public :: sparse_expand
  public :: sparse_extremal_eigenvalue
  public :: sparse_hotelling_init
  public :: sparse_hotelling_invert

  ! sparse_utils routines
  public :: sparse_init_blocking_scheme
  public :: sparse_mod_init
  public :: sparse_exit
  public :: sparse_convert
#ifdef SCALAPACK
  public :: sparse_spam3toblacs
  public :: sparse_blacstospam3
#endif
  public :: sparse_show_matrix
  public :: sparse_show_network
  public :: sparse_write
  public :: sparse_read
  public :: sparse_convert_unsegment_real
  public :: sparse_convert_segment_real

!CW
  public :: my_first_blk,my_last_blk
!END CW

  ! Enquiry routines
  interface sparse_num_element
     module procedure sparse_num_element_lib
     module procedure sparse_num_element_mat
  end interface
  interface sparse_num_elems_on_atom
     module procedure sparse_num_elems_on_atom_mat
     module procedure sparse_num_elems_on_atom_blks
  end interface

  ! Element operation routines
  interface sparse_get_element
     module procedure sparse_get_element_real
     module procedure sparse_get_element_complex
  end interface
  interface sparse_put_element
     module procedure sparse_put_element_real
     module procedure sparse_put_element_complex
  end interface

  ! Block operation routines
  interface sparse_get_block
     module procedure sparse_get_block_real
     module procedure sparse_get_block_complex
  end interface
  interface sparse_put_block
     module procedure sparse_put_block_real
     module procedure sparse_put_block_complex
  end interface

  ! Column operation routines
  interface sparse_get_col
     module procedure sparse_get_col_real
     module procedure sparse_get_col_complex
  end interface
  interface sparse_put_col
     module procedure sparse_put_col_real
     module procedure sparse_put_col_complex
  end interface
  interface sparse_clr_col
     module procedure sparse_clr_col_real
     module procedure sparse_clr_col_complex
  end interface

  ! AXPY routines
  interface sparse_axpy
     module procedure sparse_axpy_real
     module procedure sparse_axpy_complex
  end interface

  ! Scale and shift routines
  interface sparse_scale
     module procedure sparse_scale_real
     module procedure sparse_scale_complex
  end interface

  ! Conversion routines
  interface sparse_convert
     module procedure sparse_spam3tofull_real
     module procedure sparse_spam3tofull_complex
     module procedure sparse_fulltospam3_real
     module procedure sparse_fulltospam3_complex
  end interface
#ifdef SCALAPACK
  interface sparse_spam3toblacs
     module procedure sparse_spam3toblacs_real
     module procedure sparse_spam3toblacs_complex
  end interface
  interface sparse_blacstospam3
     module procedure sparse_blacstospam3_real
     module procedure sparse_blacstospam3_complex
  end interface
#endif

  ! Read/write routines
  interface sparse_write
     module procedure sparse_write_scalar
     module procedure sparse_write_vector
  end interface
  interface sparse_read
     module procedure sparse_read_scalar
     module procedure sparse_read_vector
  end interface

  ! Type definition for the structure of a block sparse matrix
  type STRUC3

    character(len=10) :: structure ! The structure code for this matrix type
    character(len=10) :: transpose_structure ! The structure code for the
                                             ! transpose of this structure
    integer :: nrows          ! The total number of rows in the matrix
    integer :: mcols          ! The total number of cols in the matrix
    integer :: nblk           ! The total number of block-rows (atoms)
    integer(kind=LONG) :: nze ! The total number of nonzero elements (<n*m)
    integer :: nzb            ! The total number of nonzero blocks stored
    integer :: my_mcols       ! The number of cols on this node
    integer :: my_nblks       ! The number of blocks on this node
    integer :: my_nze         ! The number of nonzero elements on this node
    integer :: my_nzb         ! The number of nonzero blocks on this node
    integer :: max_nze        ! Maximum number of nonzero elements on any node
    integer :: max_nzb        ! Maximum number of nonzero blocks on any node
    integer :: max_nblks      ! Maximum number of block-cols on any node
    integer :: col_blks       ! Identifier for column blocking scheme
    integer :: row_blks       ! Identifier for row blocking scheme
    integer,allocatable :: idx_lens(:) ! Lengths of index arrays on all nodes
    integer,allocatable :: idx_seg_lens(:) ! Lengths of the segments of idx
                                       ! array on other nodes corresponding to
                                       ! the segment for this node
    integer,allocatable :: mtx_seg_lens(:) ! Lengths of the segments of data
                                       ! array on other nodes corresponding to
                                       ! the segment for this node
    ! Segment information: (s_type,:) elements store whether each segment is
    ! dense, sparse or blank. (s_idx,:) elements store the starting position of
    ! the index for each segment. (s_ptr,:) elements store the pointer to the
    ! start of the data for each segment in the dmtx/zmtx arrays.
    integer,allocatable :: seg_info(:,:)
    integer,allocatable :: blk_idx(:)   ! Column indexed list of nonzero blocks
    integer,allocatable :: blk_ptr(:)   ! Pointers to nonzero blocks
  end type STRUC3

  ! Type definition for a communicator to send matrix data between MPI processes
  type, public :: COM3
    integer :: lib
    integer :: buflen
    integer :: datlen
    integer :: num_handles
    integer :: idx_hdl       ! Index of index recv msg in handles array
    integer :: info_hdl      ! Index of seginfo recv msg in handles array
    integer :: ptr_hdl       ! Index of pointer recv msg in handles array
    integer :: data_hdl      ! Index of data recv msg in handles array
    integer :: req_send_hdl  ! Index of request send msg in handles array
    integer :: data_send_hdl ! Index of data send msg in handles array
    integer :: nreqs
    integer :: reqdatlen
    logical :: send_buffer_free ! Whether or not d/zmtxsendbuf is in use
    logical :: iscmplx
    logical :: cropped
    integer, allocatable :: seginfobuf(:,:,:)
    integer, allocatable :: ptrbuf(:,:)
    integer, allocatable :: idxbuf(:,:)
    real(kind=DP), allocatable :: dmtxrecvbuf(:,:)
    real(kind=DP), allocatable :: dmtxsendbuf(:,:)
    complex(kind=DP), allocatable :: zmtxrecvbuf(:,:)
    complex(kind=DP), allocatable :: zmtxsendbuf(:,:)
    integer, allocatable :: ptrreqrecvbuf(:,:)
    integer, allocatable :: ptrreqsendbuf(:,:)
    integer, allocatable :: index_reqs(:) ! Receive buffer for index requests
    integer, allocatable :: data_reqs(:)  ! Receive buffer for data requests
    integer, allocatable :: handles(:)    ! MPI message handles for each recv
  end type COM3

  ! Type definition for the matrix elements themselves, with a pointer to
  ! the appropriate structure
  type, public :: SPAM3
    integer :: lib     ! The index for the structure in the library
    logical :: iscmplx ! TRUE if matrix is complex, otherwise FALSE
    real(kind=DP), allocatable :: dmtx(:)  ! If real, the non-zero elements of
                                           ! the matrix, stored in blocks
    complex(kind=DP), allocatable :: zmtx(:)  ! If complex, the non-zero elmnts
                                           ! of the matrix, stored in blocks
    character(len=10) :: structure         ! Character string identifying the
                                           ! sparse matrix structure
  end type SPAM3

  ! Segment info list components
  integer, parameter :: s_type = 1
  integer, parameter :: s_idx = 2
  integer, parameter :: s_ptr = 3

  ! Segment types
  integer, parameter :: SEG_DENSE      =  222
  integer, parameter :: SEG_SPARSE     =  333
  integer, parameter :: SEG_BLANK      = -111

  ! Block size list indices
  integer, public, parameter :: BLKS_NGWF      =  1
  integer, public, parameter :: BLKS_PROJ      =  2
  integer, public, parameter :: BLKS_HUB_PROJ  =  3
  integer, public, parameter :: BLKS_COND      =  4
  integer, public, parameter :: BLKS_JOINT     =  5
  integer, public, parameter :: BLKS_SW        =  6
  integer, public, parameter :: BLKS_CORE      =  7
  integer, public, parameter :: BLKS_AUX       =  8

  ! Element patterns for sparse_expand
  integer, public, parameter :: PATTERN_LOWER     = 444
  integer, public, parameter :: PATTERN_ALTERNATE = 555
  integer, public, parameter :: PATTERN_FULL      = 666

  ! Allow creation of new matrix structures?
  logical, public :: pub_sparse_allow_new_matrices

  ! Tag identifiers and special send values for sparse_product
  integer, parameter :: SEGINFO_TAG  = 200000001
  integer, parameter :: BLKIDX_TAG   = 200000002
  integer, parameter :: BLKPTR_TAG   = 200000003
  integer, parameter :: DATA_TAG     = 200000004
  integer, parameter :: IDX_REQ_TAG  = 500000001
  integer, parameter :: DATA_REQ_TAG = 500000002
  integer, parameter :: PTR_REQ_TAG  = 500000003
  integer, parameter :: index_sent = -444
  integer, parameter :: index_not_sent = -555
  integer, parameter :: index_needed = -666
  integer, parameter :: data_sent = -777
  integer, parameter :: data_not_sent = -888

  ! Library of known sparse matrix structures
  integer, parameter :: max_library = 200 ! Maximum number of indices in library
  integer :: num_library = 0             ! Number of indices in library
  type(STRUC3) :: library(max_library)   ! Library of sparse indices

  ! Record of matrices currently in use
  integer :: nalloc
  integer(kind=LONG) :: global_struc_mem
  integer(kind=LONG) :: local_struc_mem
  integer(kind=LONG) :: global_mat_mem
  integer(kind=LONG) :: local_mat_mem

  ! Parallelisation variables
  integer :: my_first_blk                ! First block-column for this node
  integer :: my_last_blk                 ! Last block-column for this node

  ! Internal arrays for commonly-used 'reverse row/column look-up'
  integer, allocatable :: found_idx(:)
  integer, allocatable :: found_ptr(:)
  logical, allocatable :: found_blk(:)

  ! Local copies of parallel strategy arrays
  integer, allocatable :: first_elem_on_atom(:,:)
  integer, allocatable :: first_elem_on_node(:,:)
  integer, allocatable :: num_elems_on_atom(:,:)
  integer, allocatable :: num_elems_on_node(:,:)
  integer, allocatable :: atom_of_elem(:,:)
  integer, allocatable :: node_of_elem(:,:)

  ! File version number
  real(kind=DP), parameter :: file_version = 1.0_DP

contains

  subroutine sparse_init_blocking_scheme(scheme,num,num_funcs_on_node, &
       num_funcs_on_atom, first_func_on_node, first_func_on_atom, &
       atom_of_func, node_of_func)

    use comms, only: pub_total_num_nodes
    use rundat, only: pub_any_nl_proj, pub_PAW, pub_hubbard, pub_hfxsw, &
         pub_cond_calculate, pub_eels_calculate, pub_use_aux_ngwfs
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: scheme
    integer, intent(in) :: num
    integer, intent(in) :: num_funcs_on_node(0:pub_total_num_nodes-1)
    integer, intent(in) :: num_funcs_on_atom(1:pub_cell%nat)
    integer, intent(in) :: first_func_on_node(0:pub_total_num_nodes)
    integer, intent(in) :: first_func_on_atom(1:pub_cell%nat)
    integer, intent(in) :: atom_of_func(1:num)
    integer, intent(in) :: node_of_func(1:num)

    ! Local Variables
    integer :: ierr
    integer :: num_elem_types
    integer :: max_elems

    ! Allocate arrays for local copies of all blocking information
    if (.not.allocated(first_elem_on_atom)) then

       ! Find number of blocking schemes
       num_elem_types = BLKS_NGWF
       if (pub_any_nl_proj.or.pub_paw) num_elem_types = BLKS_PROJ
       if (pub_hubbard) num_elem_types = BLKS_HUB_PROJ
       if (pub_cond_calculate) num_elem_types = BLKS_JOINT
       if (pub_hfxsw) num_elem_types = BLKS_SW
       if (pub_eels_calculate) num_elem_types = BLKS_CORE
       if (pub_use_aux_ngwfs) num_elem_types = BLKS_AUX

       ! Find max number of elements in any blocking scheme
       max_elems = pub_cell%num_ngwfs
       if (pub_any_nl_proj) max_elems = max(max_elems,pub_cell%num_projectors)
       if (pub_paw) max_elems = max(max_elems,pub_cell%num_pawpws)
       if (pub_hubbard) max_elems = max(max_elems,pub_cell%num_hub_proj)
       if (pub_cond_calculate) max_elems = max(max_elems, &
            pub_cell%num_ngwfs_cond+pub_cell%num_ngwfs)
       if (pub_hfxsw) max_elems = max(max_elems, &
            pub_cell%num_sw)
       if (pub_eels_calculate) max_elems = max(max_elems,pub_cell%num_corewfs)
       if (pub_use_aux_ngwfs) max_elems = max(max_elems,pub_cell%num_ngwfs_aux)

       allocate(first_elem_on_atom(pub_cell%nat,num_elem_types),stat=ierr)
       call utils_alloc_check('sparse_init','first_elem_on_atom',ierr)

       allocate(first_elem_on_node(0:pub_total_num_nodes,num_elem_types), &
            stat=ierr)
       call utils_alloc_check('sparse_init','first_elem_on_node',ierr)

       allocate(num_elems_on_atom(pub_cell%nat,num_elem_types),stat=ierr)
       call utils_alloc_check('sparse_init','num_elems_on_atom',ierr)

       allocate(num_elems_on_node(0:pub_total_num_nodes-1,num_elem_types), &
            stat=ierr)
       call utils_alloc_check('sparse_init','num_elems_on_node',ierr)

       allocate(atom_of_elem(max_elems,num_elem_types),stat=ierr)
       call utils_alloc_check('sparse_init','atom_of_elem',ierr)

       allocate(node_of_elem(max_elems,num_elem_types),stat=ierr)
       call utils_alloc_check('sparse_init','node_of_elem',ierr)
    end if

    first_elem_on_atom(:,scheme) = first_func_on_atom(:)
    first_elem_on_node(:,scheme) = first_func_on_node(:)
    num_elems_on_atom(:,scheme) = num_funcs_on_atom(:)
    num_elems_on_node(:,scheme) = num_funcs_on_node(:)
    atom_of_elem(1:num,scheme) = atom_of_func(:)
    node_of_elem(1:num,scheme) = node_of_func(:)

  end subroutine sparse_init_blocking_scheme


  !==========================================================================!
  ! This routine performs standard initialisation tasks for the module, such !
  ! as allocating temporary arrays and setting up counters and shorthands.   !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine in April 2011 from bits of sparse_init.         !
  !==========================================================================!

  subroutine sparse_mod_init

    use comms, only: pub_my_node_id, pub_total_num_nodes
    use parallel_strategy, only: parallel_strategy_list_overlaps, &
         pub_first_atom_on_node, pub_num_overlaps, pub_overlap_list
    use rundat, only: pub_any_nl_proj,pub_paw,pub_hubbard
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check

    implicit none

    ! Local Variables
    integer :: ierr

    global_mat_mem = 0
    local_mat_mem = 0
    global_struc_mem = 0
    local_struc_mem = 0
    nalloc = 0

    allocate(found_idx(pub_cell%nat),stat=ierr)
    call utils_alloc_check('sparse_init','found_idx',ierr)

    allocate(found_ptr(pub_cell%nat),stat=ierr)
    call utils_alloc_check('sparse_init','found_ptr',ierr)

    allocate(found_blk(pub_cell%nat),stat=ierr)
    call utils_alloc_check('sparse_init','found_blk',ierr)

    ! Initialise logical array of flags to FALSE
    found_blk = .false.

    ! Set shorthand variables
    my_first_blk = pub_first_atom_on_node(pub_my_node_id)
    my_last_blk = pub_first_atom_on_node(pub_my_node_id + 1) - 1

  end subroutine sparse_mod_init


  !============================================================================!
  ! This subroutine performs all the finalisation routines involving the       !
  ! sparse matrices.                                                           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   None                                                                     !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Updated for SPAM3 by Nicholas Hine, May 2009.                              !
  !============================================================================!

  subroutine sparse_exit

    use utils, only: utils_dealloc_check

    implicit none

    ! Local variables
    integer :: ierr       ! Status flag
    integer :: ilibrary   ! Counter over library entries

    ! Loop over library entries
    do ilibrary=num_library,1,-1

       ! Deallocate arrays within matrix structure
       call sparse_struc_dealloc(library(ilibrary))

    end do

    ! Library is empty
    num_library = 0

    ! Deallocate copies of parallel strategy arrays
    deallocate(node_of_elem,stat=ierr)
    call utils_dealloc_check('sparse_exit','node_of_elem',ierr)

    deallocate(atom_of_elem,stat=ierr)
    call utils_dealloc_check('sparse_exit','atom_of_elem',ierr)

    deallocate(num_elems_on_node,stat=ierr)
    call utils_dealloc_check('sparse_exit','num_elems_on_node',ierr)

    deallocate(num_elems_on_atom,stat=ierr)
    call utils_dealloc_check('sparse_exit','num_elems_on_atom',ierr)

    deallocate(first_elem_on_node,stat=ierr)
    call utils_dealloc_check('sparse_exit','first_elem_on_node',ierr)

    deallocate(first_elem_on_atom,stat=ierr)
    call utils_dealloc_check('sparse_exit','first_elem_on_atom',ierr)

    ! Deallocate memory for internal arrays
    deallocate(found_blk,stat=ierr)
    call utils_dealloc_check('sparse_exit','found_blk',ierr)

    deallocate(found_ptr,stat=ierr)
    call utils_dealloc_check('sparse_exit','found_ptr',ierr)

    deallocate(found_idx,stat=ierr)
    call utils_dealloc_check('sparse_exit','found_idx',ierr)

  end subroutine sparse_exit

  !==========================================================================!
  ! This routine returns the numbers of matrix elements and blocks (total,   !
  ! per node and per segment) which are nonzero due to overlap of pairs of   !
  ! spheres, as stored in the pub_overlap_list array.                        !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   nze     (output) : Total number of nonzero elements                    !
  !   nzb     (output) : Total number of nonzero blocks                      !
  !   my_nze  (output) : Number of nonzero elements on this node             !
  !   my_nzb  (output) : Number of nonzero blocks on this node               !
  !   seg_nze (output) : Number of nonzero elements per segment on this node !
  !   seg_nzb (output) : Number of nonzero blocks per segment on this node   !
  !--------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                     !
  ! Revised for distributed data, June 2006.                                 !
  ! Updated for SPAM3 by Nicholas Hine, May 2009.                            !
  !==========================================================================!

  subroutine sparse_count_ss(col_blks,row_blks,nze,nzb,my_nze,my_nzb, &
       seg_nze,seg_nzb)

    use comms, only: comms_reduce, pub_total_num_nodes
    use parallel_strategy, only: pub_overlap_list, pub_num_overlaps

    implicit none

    ! Arguments
    integer, intent(in) :: col_blks ! Number identifying col blk sizes
    integer, intent(in) :: row_blks ! Number identifying row blk sizes
    integer(kind=LONG), intent(out) :: nze ! Number of nonzero elements
    integer, intent(out) :: nzb     ! Number of nonzero blocks (atoms)
    integer, intent(out) :: my_nze  ! Number of nonzero elements on this proc
    integer, intent(out) :: my_nzb  ! Number of nonzero blocks on this proc
    integer, intent(out) :: seg_nze(0:pub_total_num_nodes-1)
    integer, intent(out) :: seg_nzb(0:pub_total_num_nodes-1)

    ! Local variables
    integer :: node        ! Node of each found block
    integer :: iblk, jblk  ! Atom counters
    integer :: iovlap      ! Atomic overlap counter
    integer :: ielems      ! Number of elements on atom iblk
    integer :: ibuf(2)     ! Communications buffer

    ! Zero counters
    my_nze = 0
    my_nzb = 0
    seg_nzb(:) = 0
    seg_nze(:) = 0

    ! Loop over atoms on this node (block-columns)
    do iblk=my_first_blk,my_last_blk

       ! Number of elements on atom iblk
       ielems = num_elems_on_atom(iblk,col_blks)

       ! Skip this atom if there are no elements on it - this is the case
       ! for atoms which are not Hubbard atoms when constructing V and W
       if (ielems==0) pub_num_overlaps(iblk) = 0

       ! Loop over all atoms jblk which overlap atom iblk (block-rows)
       do iovlap=1,pub_num_overlaps(iblk)
          jblk = pub_overlap_list(iovlap,iblk)

          ! Skip this atom if there are no elements on it
          if (num_elems_on_atom(jblk,row_blks)==0) cycle

          ! Count contribution to my_nze and my_nzb totals
          my_nze = my_nze + ielems * num_elems_on_atom(jblk,row_blks)
          my_nzb = my_nzb + 1

          ! Find node to which jblk belongs and record in segment nze and nzb
          node = node_of_elem(first_elem_on_atom(jblk,row_blks),row_blks)

          seg_nze(node) = seg_nze(node) + ielems * &
               num_elems_on_atom(jblk,row_blks)
          seg_nzb(node) = seg_nzb(node) + 1

       end do  ! Loop over atoms jblk which overlap atom iblk (block-rows)

    end do  ! Loop over atoms iblk on this node (block-columns)

    ! Sum up over all nodes
    ibuf(1) = my_nze
    ibuf(2) = my_nzb
    call comms_reduce('SUM',ibuf)
    nze = ibuf(1)
    nzb = ibuf(2)

  end subroutine sparse_count_ss


  !==========================================================================!
  ! This routine generates a block index for a sparse matrix from a list of  !
  ! direct sphere-sphere overlaps and adds it to the sparse library.         !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   nze     (input) : Number of nonzero elements                           !
  !   nzb     (input) : Number of nonzero blocks                             !
  !   my_nze  (inout) : Number of nonzero elements on this node              !
  !   my_nzb  (input) : Number of nonzero blocks on this node                !
  !   seg_nze (inout) : Number of nonzero elements per segment on this node  !
  !   seg_nzb (input) : Number of nonzero blocks per segment on this node    !
  !   name    (input) : Identifying name                                     !
  !--------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                     !
  ! Revised for distributed data, June 2006.                                 !
  ! Updated for SPAM3 by Nicholas Hine, May 2009.                            !
  !==========================================================================!

  subroutine sparse_index_ss(col_blks,row_blks,name,nze,nzb,my_nze,my_nzb, &
       seg_nze,seg_nzb,rlib,transpose_name)

    use comms, only: comms_abort, pub_total_num_nodes, pub_my_node_id, &
         pub_on_root
    use parallel_strategy, only: pub_overlap_list, pub_num_overlaps, &
         pub_first_atom_on_node

    implicit none

    ! Arguments
    integer, intent(in) :: col_blks       ! Number identifying col blk sizes
    integer, intent(in) :: row_blks       ! Number identifying row blk sizes
    character(len=*), intent(in) :: name  ! Identifying name
    integer(kind=LONG), intent(in) :: nze ! Number of nonzero elements
    integer, intent(in) :: nzb            ! Number of nonzero blocks
    integer, intent(inout) :: my_nze      ! Num of nonzero els on this proc
    integer, intent(in) :: my_nzb         ! Num of nonzero blocks on this proc
    integer, intent(inout) :: seg_nze(0:pub_total_num_nodes-1)
    integer, intent(in) :: seg_nzb(0:pub_total_num_nodes-1)
    integer, intent(out), optional :: rlib ! Library entry index return val
    character(len=*), intent(in), optional :: transpose_name ! Name of transpose

    ! Local variables
    integer :: my_nat     ! Number of atoms on this proc
    integer :: iblk,jblk  ! Atomic loop counters
    integer :: loc_iblk   ! Atomic loop counter in terms of local atoms
    integer :: ielems     ! Number of elements associated with atom iblk
    integer :: idx,ptr    ! Index and pointer
    integer :: dptr       ! Pointer to record (does not increment if dense)
    integer :: iovlap     ! Atomic overlap counter
    integer :: seg        ! Segment index
    integer :: seg_type   ! Segment type for this segment
    integer :: seg_start  ! Start position of this segment in the index

    ! Increment library counter
    num_library = num_library + 1
    if (num_library > max_library) &
         call sparse_library_full('sparse_index_ss')

    ! Set identification name
    library(num_library)%structure = name

    ! Set transpose name
    if (.not.present(transpose_name).and.(row_blks==col_blks)) then
       library(num_library)%transpose_structure = name
    else if (present(transpose_name)) then
       library(num_library)%transpose_structure = transpose_name
    else
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_index_ss: &
            &Tranpose structure not specified'
       call comms_abort
    end if

    ! Set row and column block size list identifiers
    library(num_library)%col_blks = col_blks
    library(num_library)%row_blks = row_blks

    ! Set segment types and increase nze accordingly
    call sparse_segments_alloc(library(num_library),my_nze,seg_nze,seg_nzb)

    ! Initialise dimensions of matrix
    call sparse_struc_alloc(library(num_library),nze,nzb,my_nze,my_nzb)

    ! Ensure list of overlaps for each atom is in ascending order
    do iblk=my_first_blk,my_last_blk
       call internal_heapsort(pub_num_overlaps(iblk),pub_overlap_list(:,iblk))
    end do

    ! Initialise counters and pointers
    my_nat = my_last_blk - my_first_blk + 1
    idx = 1
    ptr = 1
    dptr = -1

    ! Loop over segments of the index
    do seg=0,pub_total_num_nodes-1

       ! Set up entries in the seg_info(s_idx,:) and seg_info(s_ptr,:) arrays
       library(num_library)%seg_info(s_ptr,seg) = ptr
       library(num_library)%seg_info(s_idx,seg) = idx
       seg_type = library(num_library)%seg_info(s_type,seg)
       seg_start = idx

       ! Skip indexing of this segment if it contains no nonzero elements
       if (seg_type == SEG_BLANK) cycle

       if (seg/=pub_my_node_id) then
          ! Set unused pointers to zero
          library(num_library)%blk_ptr(idx:idx+my_nat) = 0
       end if

       ! Go past list of column start positions for this segment
       idx = idx + my_nat + 1
       library(num_library)%blk_idx(seg_start) = idx

       ! Loop over atom block-columns iblk on this node for this segment
       loc_iblk = 0
       do iblk=my_first_blk,my_last_blk
          loc_iblk = loc_iblk + 1

          ! Number of elements on atom iblk
          ielems = num_elems_on_atom(iblk,col_blks)

          ! Skip this atom if there are no elements on it - this is the case
          ! for atoms which are not Hubbard atoms when constructing V and W
          if (ielems==0) pub_num_overlaps(iblk) = 0

          ! Loop over overlapping block-rows (atoms)
          do iovlap=1,pub_num_overlaps(iblk)
             jblk = pub_overlap_list(iovlap,iblk)

             if (jblk < pub_first_atom_on_node(seg)) cycle
             if (jblk >= pub_first_atom_on_node(seg+1)) exit

             ! Skip this atom if there are no elements on it
             if (num_elems_on_atom(jblk,row_blks)==0) cycle

             ! Index contains block-row
             library(num_library)%blk_idx(idx) = jblk

             ! Pointer points to start of data
             if ( seg_type == SEG_SPARSE ) then
                dptr = ptr
                ! Increment block pointer by size of block
                ptr = ptr + ielems * num_elems_on_atom(jblk,row_blks)
             else if ( seg_type == SEG_DENSE ) then
                ! Find ptr to first element of this block in segment
                dptr = ptr + (first_elem_on_atom(iblk,col_blks) - &
                     first_elem_on_node(pub_my_node_id,col_blks)) * &
                     num_elems_on_node(seg,row_blks) &
                     + first_elem_on_atom(jblk,row_blks) - &
                     first_elem_on_node(seg,row_blks)
             end if

             library(num_library)%blk_ptr(idx) = dptr
             if (jblk == iblk) &
                  library(num_library)%blk_ptr(seg_start+loc_iblk-1) = dptr

             ! Increment block index by one
             idx = idx + 1

          end do  ! jblk

          ! Store idx in list of column start positions for this segment
          library(num_library)%blk_idx(seg_start+loc_iblk) = idx
          library(num_library)%blk_ptr(seg_start+loc_iblk) = 0

       end do  ! iblk

       ! For dense segments, move ptr on by total number of elements
       if ( seg_type == SEG_DENSE ) then
          ptr = ptr + num_elems_on_node(pub_my_node_id,col_blks) * &
               num_elems_on_node(seg,row_blks)
       end if

       ! Set final pointer to end of data
       library(num_library)%blk_ptr(idx) = ptr

    end do  ! seg

    ! Record end of last segment
    library(num_library)%seg_info(s_idx,pub_total_num_nodes) = idx
    library(num_library)%seg_info(s_ptr,pub_total_num_nodes) = ptr

    ! Record library value of new matrix in rlib
    if (present(rlib)) rlib = num_library

  end subroutine sparse_index_ss


  !==========================================================================!
  ! This routine returns the numbers of matrix elements and blocks (total,   !
  ! per node and per segment) in a sparse matrix whose sparsity pattern is   !
  ! the union of the nonzero elements of two already-defined sparse matrix   !
  ! indices.                                                                 !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   alib    (input)  : Library entry for first predefined matrix           !
  !   blib    (input)  : Library entry for second predefined matrix          !
  !   nze     (output) : Total number of nonzero elements                    !
  !   nzb     (output) : Total number of nonzero blocks                      !
  !   my_nze  (output) : Number of nonzero elements on this node             !
  !   my_nzb  (output) : Number of nonzero blocks on this node               !
  !   seg_nze (output) : Number of nonzero elements per segment on this node !
  !   seg_nzb (output) : Number of nonzero blocks per segment on this node   !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                     !
  !==========================================================================!

  subroutine sparse_count_union(alib,blib,nze,nzb,my_nze,my_nzb, &
       seg_nze,seg_nzb)

    use comms, only: comms_abort, comms_reduce, pub_total_num_nodes, pub_on_root

    implicit none

    ! Arguments
    integer, intent(in) :: alib
    integer, intent(in) :: blib
    integer(kind=LONG), intent(out) :: nze
    integer, intent(out) :: nzb
    integer, intent(out) :: my_nze
    integer, intent(out) :: my_nzb
    integer, intent(out) :: seg_nze(0:pub_total_num_nodes-1)
    integer, intent(out) :: seg_nzb(0:pub_total_num_nodes-1)

    ! Locals
    integer :: row_blks    ! Identifier for column blocking scheme
    integer :: col_blks    ! Identifier for row blocking scheme
    integer :: iblk,jblk   ! Atomic loop counters
    integer :: loc_iblk    ! Atomic loop counter in terms of local atoms
    integer :: iidx        ! Overlap index
    integer :: seg         ! Segment index
    integer :: aseg_start  ! Segment index start for library A
    integer :: bseg_start  ! Segment index start for library B
    integer :: aseg_type   ! Segment type for this segment of A
    integer :: bseg_type   ! Segment type for this segment of B
    integer :: num_ovlaps  ! Number of overlaps for this atom in this segment
    integer :: ielems      ! Number of elements associated with atom iblk
    integer :: ibuf(2)     ! Communications buffer

    ! Zero counters
    my_nze = 0
    my_nzb = 0

    ! Check inputs
    if ((alib<1).or.(alib>num_library)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_count_union: &
            &invalid library entry provided for alib'
       call comms_abort
    end if
    if ((blib<1).or.(blib>num_library)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_count_union: &
            &invalid library entry provided for blib'
       call comms_abort
    end if

    ! Find row and column blocking schemes
    row_blks = library(alib)%row_blks
    col_blks = library(alib)%col_blks

    ! Sanity check
    if (library(blib)%row_blks/=row_blks) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_count_union: &
            &A and B row blocking schemes do not match'
       call comms_abort
    end if
    if (library(blib)%col_blks/=col_blks) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_count_union: &
            &A and B col blocking schemes do not match'
       call comms_abort
    end if

    ! Loop over segments of both structures
    do seg=0,pub_total_num_nodes-1

       seg_nze(seg) = 0
       seg_nzb(seg) = 0
       aseg_start = library(alib)%seg_info(s_idx,seg)
       bseg_start = library(blib)%seg_info(s_idx,seg)
       aseg_type = library(alib)%seg_info(s_type,seg)
       bseg_type = library(blib)%seg_info(s_type,seg)

       ! Loop over atoms on this node
       loc_iblk = 0
       do iblk=my_first_blk,my_last_blk
          loc_iblk = loc_iblk + 1

          ielems = num_elems_on_atom(iblk,col_blks)

          ! Set counter of number of overlaps found to zero
          num_ovlaps = 0
          if (aseg_type /= SEG_BLANK) then
             ! Loop over nonzero blocks in this column in structure A
             do iidx=library(alib)%blk_idx(aseg_start+loc_iblk-1), &
                  library(alib)%blk_idx(aseg_start+loc_iblk)-1
                jblk = library(alib)%blk_idx(iidx)
                found_blk(jblk) = .true.
                num_ovlaps = num_ovlaps + 1
                found_idx(num_ovlaps) = jblk
             end do  ! iidx
          end if

          if (bseg_type /= SEG_BLANK) then
             ! Loop over nonzero blocks in this column in structure B
             do iidx=library(blib)%blk_idx(bseg_start+loc_iblk-1), &
                  library(blib)%blk_idx(bseg_start+loc_iblk)-1
                jblk = library(blib)%blk_idx(iidx)

                ! Check if this block has already been found
                if (.not.found_blk(jblk)) then

                   found_blk(jblk) = .true.
                   num_ovlaps = num_ovlaps + 1
                   found_idx(num_ovlaps) = jblk

                end if

             end do  ! iidx
          end if

          ! Loop over all overlaps kat (block-rows)
          do iidx=1,num_ovlaps
             jblk = found_idx(iidx)

             ! Whole atomic block counts
             my_nze = my_nze + ielems * num_elems_on_atom(jblk,row_blks)
             my_nzb = my_nzb + 1
             seg_nze(seg) = seg_nze(seg) + &
                  ielems * num_elems_on_atom(jblk,row_blks)
             seg_nzb(seg) = seg_nzb(seg) + 1

             ! Reset flags
             found_blk(jblk) = .false.

          end do  ! iidx

       end do  ! iblk

    end do  ! seg

    ! Sum up over all nodes
    ibuf(1) = my_nze
    ibuf(2) = my_nzb
    call comms_reduce('SUM',ibuf)
    nze = ibuf(1)
    nzb = ibuf(2)

  end subroutine sparse_count_union


  !==========================================================================!
  ! This routine generates a block index for a sparse matrix whose sparsity  !
  ! pattern is the union of the nonzero elements of two already-defined      !
  ! sparse matrix indices.                                                   !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   alib    (input) : Library entry for first predefined matrix            !
  !   blib    (input) : Library entry for second predefined matrix           !
  !   nze     (input) : Total number of nonzero elements                     !
  !   nzb     (input) : Total number of nonzero blocks                       !
  !   my_nze  (inout) : Number of nonzero elements on this node              !
  !   my_nzb  (input) : Number of nonzero blocks on this node                !
  !   seg_nze (inout) : Number of nonzero elements per segment on this node  !
  !   seg_nzb (input) : Number of nonzero blocks per segment on this node    !
  !   name    (input) : Identifying name                                     !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                     !
  !==========================================================================!

  subroutine sparse_index_union(alib_in,blib_in,name,nze,nzb,my_nze,my_nzb, &
       seg_nze,seg_nzb,rlib)

    use comms, only: comms_abort, pub_my_node_id, pub_total_num_nodes, &
         pub_on_root

    implicit none

    ! Arguments
    integer, intent(in) :: alib_in
    integer, intent(in) :: blib_in
    character(len=*), intent(in) :: name
    integer(kind=LONG), intent(in) :: nze
    integer, intent(in) :: nzb
    integer, intent(inout) :: my_nze
    integer, intent(in) :: my_nzb
    integer, intent(inout) :: seg_nze(0:pub_total_num_nodes-1)
    integer, intent(in) :: seg_nzb(0:pub_total_num_nodes-1)
    integer, intent(out), optional :: rlib ! Library entry index return val

    ! Locals
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Atomic loop counters
    integer :: loc_iblk      ! Atomic loop counter in terms of local atoms
    integer :: num_ovlaps    ! Number of overlaps for this atom in this segment
    integer :: iidx          ! Overlap index
    integer :: idx,ptr,dptr  ! Index, pointer and dense pointer
    integer :: my_nat        ! Number of atoms on this nod
    integer :: ielems        ! Number of elements associated with atom iblk
    integer :: seg           ! Segment index
    integer :: seg_type      ! Type for this segment of C
    integer :: seg_start     ! Index start position for this segment of C
    integer :: aseg_start    ! Segment index start for library A
    integer :: bseg_start    ! Segment index start for library B
    integer :: aseg_type     ! Segment type for this segment of A
    integer :: bseg_type     ! Segment type for this segment of B
    integer :: ulib          ! Library index of union
    integer :: temp_lib      ! Temporary library index
    integer :: alib          ! Local copy of alib variable
    integer :: blib          ! Local copy of alib variable

    ! Check inputs
    if ((alib_in<1).or.(alib_in>num_library)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_count_union: &
            &invalid library entry provided for alib'
       call comms_abort
    end if
    if ((blib_in<1).or.(blib_in>num_library)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_count_union: &
            &invalid library entry provided for blib'
       call comms_abort
    end if
    alib = alib_in
    blib = blib_in

    ! Find row and column blocking schemes
    row_blks = library(alib)%row_blks
    col_blks = library(alib)%col_blks

    ! Sanity check
    if (library(blib)%row_blks/=row_blks) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_count_union: &
               &A and B row blocking schemes do not match'
          call comms_abort
    end if
    if (library(blib)%col_blks/=col_blks) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_count_union: &
               &A and B col blocking schemes do not match'
          call comms_abort
    end if

    ! Check if the new structure name already exists
    call sparse_search_library(name,ulib)
    temp_lib = 0

    ! If the structure code was not found, create a new one
    if (ulib == 0) then

       ! Increment library counter
       num_library = num_library + 1
       if (num_library > max_library) &
            call sparse_library_full('sparse_index_union')

       ! Set identification name
       library(num_library)%structure = name

       ! Set row and column block size list identifiers
       library(num_library)%col_blks = col_blks
       library(num_library)%row_blks = row_blks

       ulib = num_library

    else ! this structure already exists

       ! Check if ulib is also either alib or blib
       if (ulib==alib) then
          temp_lib = num_library + 1
          call sparse_struc_copy(library(temp_lib),library(alib))
          alib = temp_lib
          if (ulib==blib) blib = temp_lib
       else if (ulib==blib) then
          temp_lib = num_library + 1
          call sparse_struc_copy(library(temp_lib),library(blib))
          blib = temp_lib
       end if

       ! Deallocate previous storage associated with ulib
       call sparse_segments_dealloc(library(ulib))
       call sparse_struc_dealloc(library(ulib))

    end if

    ! Set segment types and increase nze accordingly
    call sparse_segments_alloc(library(ulib),my_nze,seg_nze,seg_nzb)

    ! Initialise dimensions of matrix
    call sparse_struc_alloc(library(ulib),nze,nzb,my_nze,my_nzb)

    ! Initialise counters and pointers
    my_nat = my_last_blk - my_first_blk + 1
    idx = 1
    ptr = 1
    dptr = -1

    ! Loop over segments of the index
    do seg=0,pub_total_num_nodes-1

       ! Set up entries in the seg_info(s_idx,:) and seg_info(s_ptr,:) arrays
       library(ulib)%seg_info(s_ptr,seg) = ptr
       library(ulib)%seg_info(s_idx,seg) = idx
       seg_type = library(ulib)%seg_info(s_type,seg)
       seg_start = idx

       aseg_start = library(alib)%seg_info(s_idx,seg)
       bseg_start = library(blib)%seg_info(s_idx,seg)
       aseg_type = library(alib)%seg_info(s_type,seg)
       bseg_type = library(blib)%seg_info(s_type,seg)

       ! Skip indexing of this segment if it contains no nonzero elements
       if (seg_type == SEG_BLANK) cycle

       if (seg/=pub_my_node_id) then
          ! Set unused pointers to zero
          library(ulib)%blk_ptr(idx:idx+my_nat) = 0
       end if

       ! Go past list of column start positions for this segment
       idx = idx + my_nat + 1
       library(ulib)%blk_idx(seg_start) = idx

       ! Loop over atom block-columns iblk on this node for this segment
       loc_iblk = 0
       do iblk=my_first_blk,my_last_blk
          loc_iblk = loc_iblk + 1

          ! Number of elements on atom iblk
          ielems = num_elems_on_atom(iblk,col_blks)

          ! Set counter of number of overlaps found to zero
          num_ovlaps = 0

          if (aseg_type/=SEG_BLANK) then
             ! Loop over nonzero blocks in this column in structure A
             do iidx=library(alib)%blk_idx(aseg_start+loc_iblk-1), &
                  library(alib)%blk_idx(aseg_start+loc_iblk)-1
                jblk = library(alib)%blk_idx(iidx)
                found_blk(jblk) = .true.
                num_ovlaps = num_ovlaps + 1
                found_idx(num_ovlaps) = jblk
             end do  ! iidx
          end if

          if (bseg_type/=SEG_BLANK) then
             ! Loop over nonzero blocks in this column in structure B
             do iidx=library(blib)%blk_idx(bseg_start+loc_iblk-1), &
                  library(blib)%blk_idx(bseg_start+loc_iblk)-1
                jblk = library(blib)%blk_idx(iidx)

                ! Check if this block has already been found
                if (.not.found_blk(jblk)) then

                   found_blk(jblk) = .true.
                   num_ovlaps = num_ovlaps + 1
                   found_idx(num_ovlaps) = jblk

                end if

             end do  ! iidx
          end if

          ! Sort list of overlaps to ensure it is in ascending order
          call internal_heapsort(num_ovlaps,found_idx)

          do iidx=1,num_ovlaps
             jblk = found_idx(iidx)

             ! Index contains block-row
             library(ulib)%blk_idx(idx) = jblk

             ! Pointer points to start of data
             if ( seg_type == SEG_SPARSE ) then
                dptr = ptr
                ! Increment block pointer by size of block
                ptr = ptr + ielems * num_elems_on_atom(jblk,row_blks)
             else if ( seg_type == SEG_DENSE ) then
                ! Find ptr to first element of this block in segment
                dptr = ptr + (first_elem_on_atom(iblk,col_blks) - &
                     first_elem_on_node(pub_my_node_id,col_blks)) * &
                     num_elems_on_node(seg,row_blks) &
                     + first_elem_on_atom(jblk,row_blks) - &
                     first_elem_on_node(seg,row_blks)
             end if

             library(ulib)%blk_ptr(idx) = dptr
             if (jblk == iblk) &
                  library(ulib)%blk_ptr(seg_start+loc_iblk-1) = dptr

             ! Increment block index by one
             idx = idx + 1

            ! Reset flags
            found_blk(jblk) = .false.

          end do  ! jblk

          ! Store idx in list of column start positions for this segment
          library(ulib)%blk_idx(seg_start+loc_iblk) = idx
          library(ulib)%blk_ptr(seg_start+loc_iblk) = 0

       end do  ! iblk

       ! For dense segments, move ptr on by total number of elements
       if ( seg_type == SEG_DENSE ) then
          ptr = ptr + num_elems_on_node(pub_my_node_id,col_blks) * &
               num_elems_on_node(seg,row_blks)
       end if

       ! Set final pointer to end of data
       library(ulib)%blk_ptr(idx) = ptr

    end do  ! seg

    ! Record end of last segment
    library(ulib)%seg_info(s_idx,pub_total_num_nodes) = idx
    library(ulib)%seg_info(s_ptr,pub_total_num_nodes) = ptr

    if (temp_lib > 0) then
       call sparse_struc_dealloc(library(temp_lib))
    end if

    ! Record library value of new matrix in rlib
    if (present(rlib)) rlib = ulib

  end subroutine sparse_index_union


  subroutine sparse_show_segment_filling

    use comms, only: pub_on_root

    implicit none

    ! Local variables
    type(SPAM3) :: q,h,s,k,ks,ksks

    q%structure = 'Q'
    call sparse_create(q)
    h%structure = 'H'
    call sparse_create(h)
    s%structure = 'S'
    call sparse_create(s)
    k%structure = 'K'
    call sparse_create(k)
    ks%structure = 'KS'
    call sparse_create(ks)
    ksks%structure = 'KSKS'
    call sparse_create(ksks)

    if (pub_on_root) write(stdout,'(a)') 'Q'
    call internal_show_seg_frac(q)
    if (pub_on_root) write(stdout,'(a)') 'H'
    call internal_show_seg_frac(h)
    if (pub_on_root) write(stdout,'(a)') 'S'
    call internal_show_seg_frac(s)
    if (pub_on_root) write(stdout,'(a)') 'K'
    call internal_show_seg_frac(k)
    if (pub_on_root) write(stdout,'(a)') 'KS'
    call internal_show_seg_frac(ks)
    if (pub_on_root) write(stdout,'(a)') 'KSKS'
    call internal_show_seg_frac(ksks)

    call sparse_destroy(ksks)
    call sparse_destroy(ks)
    call sparse_destroy(k)
    call sparse_destroy(s)
    call sparse_destroy(h)
    call sparse_destroy(q)

  contains

    subroutine internal_show_segs(mat)

      use comms, only: pub_on_root,pub_my_node_id,pub_total_num_nodes, &
           comms_barrier,comms_bcast
      use utils, only: utils_flush,utils_alloc_check,utils_dealloc_check

      implicit none

      ! Arguments
      type(SPAM3),intent(in) :: mat

      ! Locals
      integer :: node,seg
      integer :: ierr
      character(4), allocatable :: segchar(:,:),segline(:)
      character(20) :: fmt, tmp

      allocate(segchar(0:pub_total_num_nodes-1,0:pub_total_num_nodes-1), &
           stat=ierr)
      call utils_alloc_check('internal_show_segs','segchar',ierr)
      allocate(segline(0:pub_total_num_nodes-1),stat=ierr)
      call utils_alloc_check('internal_show_segs','segline',ierr)

      do seg=0,pub_total_num_nodes-1
         if (library(mat%lib)%seg_info(s_type,seg)==SEG_SPARSE) then
            segchar(seg,pub_my_node_id) = '+'
         else if (library(mat%lib)%seg_info(s_type,seg)==SEG_DENSE) then
            segchar(seg,pub_my_node_id) = '#'
         else
            segchar(seg,pub_my_node_id) = '·'
         end if
      end do

      do node=0,pub_total_num_nodes-1
         call comms_bcast(node,segchar(:,node),pub_total_num_nodes)
      end do

      if (pub_on_root) then
         write(stdout,'(a)',advance='no') '     '
         do node=0,pub_total_num_nodes-1
            write(stdout,'(i2)',advance='no') node
         end do
         write(stdout,*)
         write(tmp,'(i10)') pub_total_num_nodes
         write(fmt,'(a,a,a)')'(i5,',trim(adjustl(tmp)),'a2)'
         do seg=0,pub_total_num_nodes-1
            segline(:) = segchar(seg,:)
            write(stdout,fmt) seg,segline(:)
         end do
      end if

      call comms_barrier
      call utils_flush
      call comms_barrier

      deallocate(segline,stat=ierr)
      call utils_dealloc_check('internal_show_segs','segline',ierr)
      deallocate(segchar,stat=ierr)
      call utils_dealloc_check('internal_show_segs','segchar',ierr)

    end subroutine internal_show_segs

    subroutine internal_show_seg_frac(mat)

      use comms, only: pub_on_root,pub_my_node_id,pub_total_num_nodes, &
           comms_barrier,comms_bcast
      use utils, only: utils_flush,utils_alloc_check,utils_dealloc_check

      implicit none

      ! Arguments
      type(SPAM3),intent(in) :: mat

      ! Locals
      integer :: node,seg
      integer :: ierr
      real(kind=DP), allocatable :: segfrac(:,:),segline(:)
      character(20) :: fmt, tmp

      allocate(segfrac(0:pub_total_num_nodes-1,0:pub_total_num_nodes-1), &
           stat=ierr)
      call utils_alloc_check('internal_show_seg_frac','segchar',ierr)
      allocate(segline(0:pub_total_num_nodes-1),stat=ierr)
      call utils_alloc_check('internal_show_seg_frac','segline',ierr)

      do seg=0,pub_total_num_nodes-1
         segfrac(seg,pub_my_node_id) = (library(mat%lib)%seg_info(s_ptr,seg+1) - &
              library(mat%lib)%seg_info(s_ptr,seg))/ &
              real(num_elems_on_node(pub_my_node_id,library(mat%lib)%col_blks)*&
              num_elems_on_node(seg,library(mat%lib)%row_blks),kind=DP)
      end do

      do node=0,pub_total_num_nodes-1
         call comms_bcast(node,segfrac(:,node),pub_total_num_nodes)
      end do

      if (pub_on_root) then
         write(stdout,'(a)',advance='no') '     '
         do node=0,pub_total_num_nodes-1
            write(stdout,'(i5)',advance='no') node
         end do
         write(stdout,*)
         write(tmp,'(i10)') pub_total_num_nodes
         write(fmt,'(a,a,a)')'(i5,',trim(adjustl(tmp)),'f5.2)'
         do seg=0,pub_total_num_nodes-1
            segline(:) = segfrac(seg,:)
            write(stdout,fmt) seg,segline(:)
         end do
      end if

      call comms_barrier
      call utils_flush
      call comms_barrier

      deallocate(segline,stat=ierr)
      call utils_dealloc_check('internal_show_seg_frac','segline',ierr)
      deallocate(segfrac,stat=ierr)
      call utils_dealloc_check('internal_show_seg_frac','segchar',ierr)

    end subroutine internal_show_seg_frac

  end subroutine sparse_show_segment_filling

  !==========================================================================!
  ! This subroutine sorts a list of integers (to ensure that blocks are      !
  ! always indexed in ascending order.                                       !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   n     (input) : The length of the list                                 !
  !   list  (inout) : The list to be sorted                                  !
  !--------------------------------------------------------------------------!
  ! Written by Peter Haynes, June 2006.                                      !
  !==========================================================================!

  subroutine internal_heapsort(n,list)

    implicit none

    ! Arguments
    integer, intent(in)    :: n
    integer, intent(inout) :: list(:)

    ! Local variables
    integer :: i,j                 ! Loop variables over list
    integer :: ic,is               ! Counters for heap creation/selection
    integer :: temp                ! Temporary copy of list item

    ! If there is only one item in the list, there's no sorting to be done!
    if (n <= 1) return

    ! Do the heap sort
    ic = n/2+1
    is = n
    do
       if (ic > 1) then ! in heap creation phase
          ic = ic - 1
          temp = list(ic)
       else ! in heap selection phase
          temp = list(is)
          list(is) = list(1)
          is = is - 1
          if (is == 1) then
             list(1) = temp
             exit
          end if
       end if
       ! Sift down temporary copy to correct level in heap
       i = ic
       j = 2 * ic
       do while (j <= is)
          if (j < is) then
             if (list(j) < list(j+1)) j = j + 1
          end if
          if (temp < list(j)) then ! demote temporary copy
             list(i) = list(j)
             i = j
             j = 2 * j
          else ! found level for temporary copy
             j = is + 1
          end if
       end do
       list(i) = temp
    end do

  end subroutine internal_heapsort


  !============================================================================!
  ! This subroutine creates a new sparse matrix by:                            !
  !   (i) looking its name up in the library of known structures, or           !
  !  (ii) copying the structure from another matrix, or                        !
  ! (iii) generating a new structure from the product of two sparse matrices.  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   new_spam (output) : The sparse matrix to be created                      !
  !   spama    (input)  : An optional sparse matrix whose structure is used    !
  !   spamb    (input)  : A second (optional) sparse matrix whose structure    !
  !                       will be used                                         !
  !   iscmplx  (input)  : A flag to indicate the data stored in this matrix is !
  !                       complex                                              !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009, based on SPAM2 sparse_create routine,  !
  ! Written by Peter Haynes, March 2004,                                       !
  ! Revised for distributed data, June 2006,                                   !
  ! Modified for dense matrices by Nicholas Hine, Dec 2007.                    !
  !============================================================================!

  subroutine sparse_create(new_spam,spama,spamb,iscmplx,rlib)

    use comms, only: comms_abort, comms_reduce, comms_send, &
         comms_barrier, comms_free, comms_bcast, comms_probe, comms_wait, &
         pub_on_root, pub_my_node_id, pub_total_num_nodes, pub_null_handle
    use parallel_strategy, only: pub_first_atom_on_node, &
         pub_num_atoms_on_node
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: new_spam      ! Sparse matrix to be created
    type(SPAM3), optional, intent(in) :: spama  ! First (optional) matrix to use
    type(SPAM3), optional, intent(in) :: spamb  ! Other (optional) matrix to use
    logical, optional, intent(in) :: iscmplx    ! Flag for complex data (opt)
    integer, optional, intent(out) :: rlib      ! Optional library index

    ! Local variables
    logical :: loc_iscmplx ! Local copy of optional iscmplx argument
    integer :: ierr        ! Error flag
    integer :: ilib        ! Library entry
    integer(kind=LONG) :: nze ! Number of nonzero elements in new matrix
    integer :: nzb         ! Number of nonzero blocks in new matrix
    integer :: my_nze      ! Number of nonzero els in new matrix on this proc
    integer :: my_nzb      ! Number of nonzero blks in new matrix on this proc
    integer :: alib,blib   ! Library entries for spama and spamb
    integer, allocatable :: num_ovlaps(:) ! Num of overlaps found for each atom
    integer,allocatable :: seg_nzb(:) ! Number of nonzero blocks on each seg
    integer,allocatable :: seg_nze(:) ! Number of nonzero elements on each seg
    type(COM3) :: acom

    if ((.not. present(spama)) .and. (.not. present(spamb))) then

       ! Check name of new matrix against library of known structures

       ! Search library for matching entry
       call sparse_search_library(new_spam%structure,ilib)

       ! Quit if not found
       if (ilib == 0) then
          if (pub_on_root) write(stdout,'(2a)') 'Error in sparse_create: &
               &no library structure found to match ',trim(new_spam%structure)
          call comms_abort
       end if

       ! Copy structure from this library entry
       new_spam%lib = ilib

       ! Set default for optional argument
       loc_iscmplx = .false.
       if (present(iscmplx)) loc_iscmplx = iscmplx

       ! Allocate space for data
       call sparse_data_alloc(new_spam,loc_iscmplx)

    else if (.not. present(spamb)) then

       ! Copy structure from given sparse matrix: spama
       new_spam%lib = spama%lib

       ! Set default for optional argument
       loc_iscmplx = spama%iscmplx
       if (present(iscmplx)) loc_iscmplx = iscmplx

       ! Allocate space for data
       call sparse_data_alloc(new_spam,loc_iscmplx)

    else
       ! Generate a new structure from the product of spama and spamb

       ! Set default for optional argument
       loc_iscmplx = (spama%iscmplx .or. spamb%iscmplx)
       if (present(iscmplx)) loc_iscmplx = iscmplx

       ! Generate the new structure identifier
       write(new_spam%structure,'(a8)') trim(spama%structure)// &
            trim(spamb%structure)
       new_spam%structure = adjustl(new_spam%structure)

       ! Check if this already exists in the library
       call sparse_search_library(new_spam%structure,ilib)

       ! If it does...
       if (ilib > 0) then

          ! Copy the known structure from the library
          new_spam%lib = ilib

          ! Allocate space for data
          call sparse_data_alloc(new_spam,loc_iscmplx)

       else    ! ... otherwise some work must be done!

          if (.not.pub_sparse_allow_new_matrices) then
             !if (pub_on_root) then
             !   write(stdout,'(a)') 'WARNING in sparse_create: &
             !        &Creation of new matrix structures is discouraged'
             !   write(stdout,'(a)') 'outside of the routine &
             !        &sparse_init_rep'
             !   write(stdout,'(6a)') 'Matrix Structures: ', &
             !        trim(spama%structure),',',trim(spamb%structure), ',', &
             !        trim(new_spam%structure)
             !end if
          end if

          ! Find library entries for spama and spamb
          alib = spama%lib ; blib = spamb%lib

          if (library(alib)%col_blks /= library(blib)%row_blks) then
             if (pub_on_root) then
                write(stdout,'(4a)') 'Error in sparse_create: &
                     &column block scheme of matrix structure ', &
                     library(alib)%structure
                write(stdout,'(3a)')'does not match row block scheme of &
                     &matrix structure ', library(blib)%structure
                write(stdout,'(a)')'Matrix product cannot be indexed'
             end if
             call comms_abort
          end if

          ! Increment library counter
          num_library = num_library + 1
          if (num_library > max_library) &
               call sparse_library_full('sparse_create')

          ! Enter structure into library
          library(num_library)%structure = new_spam%structure
          library(num_library)%col_blks = library(blib)%col_blks
          library(num_library)%row_blks = library(alib)%row_blks

          ! Allocate workspace
          allocate(num_ovlaps(my_last_blk-my_first_blk+1),stat=ierr)
          call utils_alloc_check('sparse_create','num_ovlaps',ierr)
          allocate(seg_nzb(0:pub_total_num_nodes-1),stat=ierr)
          call utils_alloc_check('sparse_create','seg_nzb',ierr)
          allocate(seg_nze(0:pub_total_num_nodes-1),stat=ierr)
          call utils_alloc_check('sparse_create','seg_nze',ierr)

          ! Allocate arrays and initialise comms
          call sparse_com_allocate(acom,spama,2, &
               alloc_mtx=.false.,cropped=.false.,seg=.false.)

          ! Count the number of nonzero blocks and elements using the library
          call internal_count_product(nze,nzb,my_nze,my_nzb)

          ! Allocate arrays for the structure segment arrays
          call sparse_segments_alloc(library(num_library),my_nze,seg_nze, &
               seg_nzb)

          ! Allocate arrays for matrix structure
          call sparse_struc_alloc(library(num_library),nze,nzb,my_nze,my_nzb)

          ! Generate the new block sparse index
          call internal_index_product

          ! Allocate arrays and initialise comms
          call sparse_com_deallocate(acom,dealloc_mtx=.false.,cropped=.false.)

          ! Deallocate workspace
          deallocate(seg_nze,stat=ierr)
          call utils_dealloc_check('sparse_create','seg_nze',ierr)
          deallocate(seg_nzb,stat=ierr)
          call utils_dealloc_check('sparse_create','seg_nzb',ierr)
          deallocate(num_ovlaps,stat=ierr)
          call utils_dealloc_check('sparse_create','num_ovlaps',ierr)

          ! Point new_spam to new structure in the library
          new_spam%lib = num_library

          ! Allocate space for data
          call sparse_data_alloc(new_spam,loc_iscmplx)

       end if

    end if

    ! Return library index if requested
    if (present(rlib)) rlib = new_spam%lib

  contains

    !==========================================================================!
    ! This subroutine counts the nonzero elements and blocks in a structure    !
    ! derived from the product of two existing sparse matrices.                !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   nze    (output) : Number of nonzero elements                           !
    !   nzb    (output) : Number of nonzero blocks                             !
    !   my_nze (output) : Number of nonzero elements on this node              !
    !   my_nzb (output) : Number of nonzero blocks on this node                !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine, May 2009, based on the SPAM3 routine:          !
    ! Written by Peter Haynes, March 2004,                                     !
    ! Revised for distributed data, June 2006,                                 !
    ! Modified for dense matrices by Nicholas Hine, Dec 2007.                  !
    !==========================================================================!

    subroutine internal_count_product(nze,nzb,my_nze,my_nzb)

      implicit none

      ! Arguments
      integer(kind=LONG), intent(out) :: nze ! Number of nonzero elements
      integer, intent(out) :: nzb    ! Number of nonzero blocks
      integer, intent(out) :: my_nze ! Number of nonzero elements on this proc
      integer, intent(out) :: my_nzb ! Number of nonzero blocks on this proc

      ! Local variables
      integer :: col_blks             ! Identifier for column blocking scheme
      integer :: row_blks             ! Identifier for row blocking scheme
      integer :: step                 ! Node step loop counter
      integer :: reqnode              ! Node to request index from for next step
      integer :: recvnode             ! Node to receive index from this step
      integer :: my_nat               ! Number of atoms on this node
      integer :: iblk,jblk,kblk       ! Block counters
      integer :: loc_iblk,loc_jblk    ! Block counters local to node
      integer :: iovlap,jovlap,kovlap ! Block overlap counters
      integer :: ielems               ! Number of elements on atom iblk
      integer :: iblk_ovlaps          ! Number of overlaps involving atom iblk
      integer :: max_ovlaps           ! Maximum num overlaps found for any atom
      integer :: seg                  ! Segment index
      integer :: b_seg_type           ! Segment type for this segment of B
      integer :: seg_start            ! Start position of this segment in index
      integer, allocatable :: ovlap_list(:,:) ! Ovlap list for atoms on my proc
      integer, allocatable :: node_nze(:) ! Numbers of nonzero elements on nodes

      ! Zero counters
      num_ovlaps = 0

      ! First obtain a quick over-estimate of the number of overlaps
      ! involving each local atom
      call comms_barrier
      call comms_free

      ! Initialise amat comms
      acom%index_reqs(:) = index_not_sent
      call sparse_init_comms(acom)

      ! Find block sizes from the product matrices
      col_blks = library(num_library)%col_blks
      row_blks = library(num_library)%row_blks

      ! Loop over other processors for whole columns of A
      do step=0,pub_total_num_nodes-1

         ! Receive data for this step
         recvnode = modulo(pub_my_node_id-step+pub_total_num_nodes, &
              pub_total_num_nodes)
         call sparse_get_step_index(acom,recvnode)

         ! Request index and data for next step if required
         reqnode = modulo(pub_my_node_id-step-1+pub_total_num_nodes, &
              pub_total_num_nodes)
         if (reqnode/=pub_my_node_id) then
            ! Send request for index, pointers and data
            call comms_send(reqnode,index_needed,1,tag=IDX_REQ_TAG)
            ! Start asynchronous receives for index and pointers
            call sparse_recv_index(acom,reqnode,1,.false.,async=.true.)
         end if

         ! Find segment type of this part of index of B and cycle if blank
         b_seg_type = library(blib)%seg_info(s_type,recvnode)
         if (b_seg_type == SEG_BLANK) cycle

         ! Loop over atom block-columns iblk of B on this node
         loc_iblk = 0
         do iblk=my_first_blk,my_last_blk
            loc_iblk = loc_iblk + 1

            ! Reset counter for number of overlaps for atom iblk
            iblk_ovlaps = 0

            ! Loop over atoms jblk overlapping atom iblk: block-rows in
            ! block-col iblk of B
            seg_start = library(blib)%seg_info(s_idx,recvnode)

            do iovlap=library(blib)%blk_idx(seg_start+loc_iblk-1), &
                 library(blib)%blk_idx(seg_start+loc_iblk)-1
               jblk = library(blib)%blk_idx(iovlap)

               ! Loop over atoms kblk overlapping atom jblk: block-rows in
               ! block-column jblk of A
               do seg=0,pub_total_num_nodes-1

                  ! If this segment of A is blank, move on to next
                  if (acom%seginfobuf(s_type,seg,2)==SEG_BLANK) cycle

                  loc_jblk = jblk - pub_first_atom_on_node(recvnode) + 1
                  do jovlap=acom%idxbuf(acom%seginfobuf(s_idx,seg,2)+loc_jblk-1,2), &
                       acom%idxbuf(acom%seginfobuf(s_idx,seg,2)+loc_jblk,2)-1
                     kblk = acom%idxbuf(jovlap,2)

                     ! Check whether this sph-sph overlap has already been found
                     if (.not. found_blk(kblk)) then

                        ! Mark block as found
                        found_blk(kblk) = .true.
                        iblk_ovlaps = iblk_ovlaps + 1
                        found_idx(iblk_ovlaps) = kblk

                     end if

                  end do  ! kblk

               end do ! seg

            end do  ! jblk

            ! Keep count of total number of overlaps involving iblk found
            num_ovlaps(loc_iblk) = num_ovlaps(loc_iblk) + iblk_ovlaps

            ! Reset flags for overlaps found
            do iovlap=1,iblk_ovlaps
               kblk = found_idx(iovlap)
               found_blk(kblk) = .false.
            end do

            call sparse_check_send_requests(acom)

         end do  ! Loop over atoms iblk on this node

      end do  ! Loop over nodes

      ! Finalise amat comms
      call sparse_exit_comms(acom,spama)

      ! Set up holding array to store overlap details
      max_ovlaps = maxval(num_ovlaps) + 1
      my_nat = my_last_blk - my_first_blk + 1
      allocate(ovlap_list(max_ovlaps,my_nat),stat=ierr)
      call utils_alloc_check('internal_count_product (sparse_create)', &
           'ovlap_list',ierr)
      ovlap_list = 0

      ! Now count again, avoiding over-counting between processors

      ! Zero counters
      my_nze = 0
      my_nzb = 0
      num_ovlaps = 0
      seg_nzb(:) = 0
      seg_nze(:) = 0
      call comms_barrier

      ! Initialise amat comms
      acom%index_reqs(:) = index_not_sent
      call sparse_init_comms(acom)

      ! Loop over processors
      do step=0,pub_total_num_nodes-1

         ! Receive data for this step
         recvnode = modulo(pub_my_node_id-step+pub_total_num_nodes, &
              pub_total_num_nodes)
         call sparse_get_step_index(acom,recvnode)

         ! Request index and data for next step if required
         reqnode = modulo(pub_my_node_id-step-1+pub_total_num_nodes, &
              pub_total_num_nodes)
         if (reqnode/=pub_my_node_id) then
            ! Send request for index, pointers and data
            call comms_send(reqnode,index_needed,1,tag=IDX_REQ_TAG)
            ! Start asynchronous receives for index and pointers
            call sparse_recv_index(acom,reqnode,1,.false.,async=.true.)
         end if

         ! Find segment type of this part of index of B and cycle if blank
         b_seg_type = library(blib)%seg_info(s_type,recvnode)
         if (b_seg_type == SEG_BLANK) cycle

         ! Loop over atom block-columns iblk on this node
         loc_iblk = 0
         do iblk=my_first_blk,my_last_blk
            loc_iblk = loc_iblk + 1

            ! Reset counter for number of overlaps for atom iblk
            iblk_ovlaps = 0

            ! Number of elements on atom iblk
            ielems = num_elems_on_atom(iblk,col_blks)

            ! Loop over atoms jblk overlapping atom iblk: block-rows in
            ! block-col iblk of B
            seg_start = library(blib)%seg_info(s_idx,recvnode)
            do iovlap=library(blib)%blk_idx(seg_start+loc_iblk-1), &
                 library(blib)%blk_idx(seg_start+loc_iblk)-1
               jblk = library(blib)%blk_idx(iovlap)

               ! Loop over atoms kblk overlapping atom jblk: block-rows in
               ! block-column jblk of A
               do seg=0,pub_total_num_nodes-1

                  ! If this segment of A is blank, move on to next
                  if (acom%seginfobuf(s_type,seg,2)==SEG_BLANK) cycle

                  loc_jblk = jblk - pub_first_atom_on_node(recvnode) + 1
                  do jovlap=acom%idxbuf(acom%seginfobuf(s_idx,seg,2)+loc_jblk-1,2), &
                       acom%idxbuf(acom%seginfobuf(s_idx,seg,2)+loc_jblk,2)-1
                     kblk = acom%idxbuf(jovlap,2)

                     ! Check whether this sph-sph overlap has already been found
                     if (.not. found_blk(kblk)) then

                        ! Mark block as found
                        found_blk(kblk) = .true.
                        iblk_ovlaps = iblk_ovlaps + 1
                        found_idx(iblk_ovlaps) = kblk

                     end if

                  end do  ! kblk

               end do ! seg

            end do  ! jblk

            ! Loop over overlaps found for block iblk (i.e. block-rows kblk in
            ! block-column iblk of the product AB)
            do iovlap=1,iblk_ovlaps
               kblk = found_idx(iovlap)

               ! Find position for kblk in the list
               do jovlap=1,num_ovlaps(loc_iblk)
                  if (ovlap_list(jovlap,loc_iblk) >= kblk) exit
               end do

               ! Check that kblk does not already occur in the list
               if (.not. ovlap_list(jovlap,loc_iblk) == kblk) then

                  ! Increase number of overlaps found for iblk
                  num_ovlaps(loc_iblk) = num_ovlaps(loc_iblk) + 1

                  ! Shuffle up entries for atoms >kblk
                  do kovlap=num_ovlaps(loc_iblk),jovlap+1,-1
                     ovlap_list(kovlap,loc_iblk) = ovlap_list(kovlap-1,loc_iblk)
                  end do

                  ! Insert entry for kblk
                  ovlap_list(jovlap,loc_iblk) = kblk

                  ! Count elements and blocks
                  my_nze = my_nze + ielems * num_elems_on_atom(kblk,row_blks)
                  my_nzb = my_nzb + 1

                  ! Count contributions to this segment of product AB
                  seg = node_of_elem(first_elem_on_atom(kblk,row_blks),row_blks)
                  seg_nzb(seg) = seg_nzb(seg) + 1
                  seg_nze(seg) = seg_nze(seg) + &
                       ielems * num_elems_on_atom(kblk,row_blks)

               end if

               ! Reset flags
               found_blk(kblk) = .false.

            end do

            call sparse_check_send_requests(acom)

         end do  ! Loop over atoms iblk on this node

      end do  ! Loop over nodes

      ! Finalise amat comms
      call sparse_exit_comms(acom,spama)

      ! Allocate array to calculate non-overflowing nze
      allocate(node_nze(0:pub_total_num_nodes-1),stat=ierr)
      call utils_alloc_check('internal_count_product (sparse_create)', &
           'node_nze',ierr)

      ! ndmh: the following calculates nze allowing for it to exceed the
      ! length of a standard integer. Assumes only nze may require kind=LONG,
      ! but that the individual my_nze values can be stored as standard
      ! integers

      ! Get node_nze from each node
      node_nze(:) = 0
      node_nze(pub_my_node_id) = my_nze
      call comms_reduce('SUM',node_nze(0:pub_total_num_nodes-1),&
           pub_total_num_nodes)

      ! ndmh: non-overflowing way to add up nze (defined LONG)
      nze = 0
      do recvnode=0,pub_total_num_nodes-1
         nze = nze + int(node_nze(recvnode),kind=LONG)
      end do

      ! Dellocate array to calculate non-overflowing nze
      deallocate(node_nze,stat=ierr)
      call utils_dealloc_check('internal_count_product (sparse_create)', &
           'node_nze',ierr)

      ! Sum up nzb over all nodes
      acom%idxbuf(1,1) = my_nzb
      call comms_reduce('SUM',acom%idxbuf,1)
      nzb = acom%idxbuf(1,1)

      ! Deallocate holding array to store overlap details
      deallocate(ovlap_list,stat=ierr)
      call utils_dealloc_check('internal_count_product (sparse_create)', &
           'ovlap_list',ierr)

    end subroutine internal_count_product

    !==========================================================================!
    ! This subroutine indexes the nonzero blocks in a structure derived from   !
    ! the product of two existing sparse matrices.                             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   None                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine, May 2009, based on the SPAM3 routine:          !
    ! Written by Peter Haynes, March 2004,                                     !
    ! Revised for distributed data, June 2006,                                 !
    ! Modified for dense matrices by Nicholas Hine, Dec 2007.                  !
    !==========================================================================!

    subroutine internal_index_product

      implicit none

      ! Local variables
      integer :: col_blks             ! Identifier for column blocking scheme
      integer :: row_blks             ! Identifier for row blocking scheme
      integer :: step                 ! Node step loop counter
      integer :: reqnode              ! Node to request index from for next step
      integer :: recvnode             ! Node to receive index from this step
      integer :: my_nat               ! Number of atoms on this node
      integer :: iblk,jblk,kblk       ! Atom block counters
      integer :: loc_iblk,loc_jblk    ! Atom counters local to node
      integer :: iovlap,jovlap,kovlap ! Block overlap counters
      integer :: ielems               ! Number of elements on atom iblk
      integer :: iblk_ovlaps          ! Number of overlaps found for this atom
      integer :: max_ovlaps           ! Maximum num overlaps found for any atom
      integer :: idx                  ! Index counter
      integer :: ptr                  ! Pointer counter
      integer :: dptr                 ! Dense segment pointer
      integer :: seg                  ! Segment index
      integer :: c_seg_type           ! Segment type for this segment
      integer :: seg_start            ! Start position of this segment in index
      integer, allocatable :: ovlap_list(:,:) ! Ovlap list for atoms on my proc
      integer :: handles(2)

      ! Find block sizes from the product matrix
      col_blks = library(num_library)%col_blks
      row_blks = library(num_library)%row_blks

      ! Set up holding array to store overlap details
      max_ovlaps = maxval(num_ovlaps) + 1
      allocate(ovlap_list(max_ovlaps,pub_num_atoms_on_node(pub_my_node_id)),&
           stat=ierr)
      call utils_alloc_check('internal_index_product (sparse_create)', &
           'ovlap_list',ierr)
      ovlap_list = 0
      num_ovlaps = 0
      dptr = -1

      ! Initialise amat comms
      acom%index_reqs(:) = index_not_sent
      call sparse_init_comms(acom)

      ! Loop over processors
      do step=0,pub_total_num_nodes-1

         ! Receive data for this step
         recvnode = modulo(pub_my_node_id-step+pub_total_num_nodes, &
              pub_total_num_nodes)
         call sparse_get_step_index(acom,recvnode)

         ! Request index and data for next step if required
         reqnode = modulo(pub_my_node_id-step-1+pub_total_num_nodes, &
              pub_total_num_nodes)
         if (reqnode/=pub_my_node_id) then
            ! Send request for index, pointers and data
            call comms_send(reqnode,index_needed,1,tag=IDX_REQ_TAG)
            ! Start asynchronous receives for index and pointers
            call sparse_recv_index(acom,reqnode,1,.false.,async=.true.)
         end if

         ! Find segment type of this part of index of B and cycle if blank
         if (library(blib)%seg_info(s_type,recvnode) == SEG_BLANK) cycle

         ! Loop over atom block-columns iblk on this node
         loc_iblk = 0
         do iblk=my_first_blk,my_last_blk
            loc_iblk = loc_iblk + 1

            ! Reset counter of number of overlaps for atom iblk
            iblk_ovlaps = 0

            ! Loop over atoms jblk overlapping atom iblk according to spamb
            seg_start = library(blib)%seg_info(s_idx,recvnode)
            do iovlap=library(blib)%blk_idx(seg_start+loc_iblk-1), &
                 library(blib)%blk_idx(seg_start+loc_iblk)-1
               jblk = library(blib)%blk_idx(iovlap)

               ! Loop over atoms kblk overlapping atom jblk: block-rows in
               ! block-column jblk of spama
               do seg=0,pub_total_num_nodes-1
                  loc_jblk = jblk - pub_first_atom_on_node(recvnode) + 1

                  ! If this segment of A is blank, move on
                  if (acom%seginfobuf(s_type,seg,2)==SEG_BLANK) cycle

                  do jovlap=acom%idxbuf(acom%seginfobuf(s_idx,seg,2)+loc_jblk-1,2), &
                       acom%idxbuf(acom%seginfobuf(s_idx,seg,2)+loc_jblk,2)-1
                     kblk = acom%idxbuf(jovlap,2)

                     ! Check whether this sph-sph overlap has already been found
                     if (.not. found_blk(kblk)) then

                        ! Mark block as found
                        found_blk(kblk) = .true.
                        iblk_ovlaps = iblk_ovlaps + 1
                        found_idx(iblk_ovlaps) = kblk

                     end if

                  end do  ! kblk

               end do  ! seg

            end do  ! jblk

            ! Loop over overlaps found for atom iblk
            do iovlap=1,iblk_ovlaps
               kblk = found_idx(iovlap)

               ! Find position for kblk in the list
               do jovlap=1,num_ovlaps(loc_iblk)
                  if (ovlap_list(jovlap,loc_iblk) >= kblk) exit
               end do

               ! Check that kblk does not already occur in the list
               if (.not. ovlap_list(jovlap,loc_iblk) == kblk) then

                  ! Increase number of overlaps found for iblk
                  num_ovlaps(loc_iblk) = num_ovlaps(loc_iblk) + 1

                  ! Shuffle up entries for atoms >kblk
                  do kovlap=num_ovlaps(loc_iblk),jovlap+1,-1
                     ovlap_list(kovlap,loc_iblk) = ovlap_list(kovlap-1,loc_iblk)
                  end do

                  ! Insert entry for kblk
                  ovlap_list(jovlap,loc_iblk) = kblk

               end if

               ! Reset flag
               found_blk(kblk) = .false.

            end do  ! Overlaps for atom iblk

            call sparse_check_send_requests(acom)

         end do  ! Loop over atoms iblk

      end do  ! Loop over steps

      ! Finalise amat comms
      call sparse_exit_comms(acom,spama)

      ! Initialise counters and pointers
      my_nat = my_last_blk - my_first_blk + 1
      idx = 1
      ptr = 1

      ! Loop over segments of the index of the product AB
      do seg=0,pub_total_num_nodes-1

         ! Set up entries in the seg_info(s_idx,:) and seg_info(s_ptr,:) arrays
         library(num_library)%seg_info(s_ptr,seg) = ptr
         library(num_library)%seg_info(s_idx,seg) = idx
         c_seg_type = library(num_library)%seg_info(s_type,seg)
         seg_start = idx

         ! Skip indexing of this segment if it contains no nonzero elements
         if (c_seg_type == SEG_BLANK) cycle

         if (seg/=pub_my_node_id) then
            ! Set unused pointers to zero
            library(num_library)%blk_ptr(idx:idx+my_nat) = 0
         end if

         ! Go past list of column start positions for this segment
         idx = idx + my_nat + 1
         library(num_library)%blk_idx(seg_start) = idx

         ! Loop over atom block-columns iblk on this node for this segment
         loc_iblk = 0
         do iblk=my_first_blk,my_last_blk
            loc_iblk = loc_iblk + 1

            ! Number of elements on atom iblk
            ielems = num_elems_on_atom(iblk,col_blks)

            ! Loop over overlaps found with this atom
            do iovlap=1,num_ovlaps(loc_iblk)
               jblk = ovlap_list(iovlap,loc_iblk)

               if (jblk < pub_first_atom_on_node(seg)) cycle
               if (jblk >= pub_first_atom_on_node(seg+1)) exit

               ! Index contains block-row
               library(num_library)%blk_idx(idx) = jblk

               ! Pointer points to start of data
               if ( c_seg_type == SEG_SPARSE ) then
                  dptr = ptr
                  ! Increment block pointer by size of block
                  ptr = ptr + ielems * num_elems_on_atom(jblk,row_blks)
               else if ( c_seg_type == SEG_DENSE ) then
                  ! Find ptr to first element of this block in segment
                  dptr = ptr + (first_elem_on_atom(iblk,col_blks) - &
                       first_elem_on_node(pub_my_node_id,col_blks)) * &
                       num_elems_on_node(seg,row_blks) &
                       + first_elem_on_atom(jblk,row_blks) - &
                       first_elem_on_node(seg,row_blks)
               end if

               library(num_library)%blk_ptr(idx) = dptr
               if (jblk == iblk) &
                    library(num_library)%blk_ptr(seg_start+loc_iblk-1) = dptr

               ! Increment block index by one
               idx = idx + 1

            end do  ! jblk

            ! Store idx in list of column start positions for this segment
            library(num_library)%blk_idx(seg_start+loc_iblk) = idx
            library(num_library)%blk_ptr(seg_start+loc_iblk) = 0

         end do  ! iblk

         ! For dense segments, move ptr on by total number of elements
         if ( c_seg_type == SEG_DENSE ) then
            ptr = ptr + num_elems_on_node(pub_my_node_id,col_blks) * &
                 num_elems_on_node(seg,row_blks)
         end if

         ! Set final pointer to end of data
         library(num_library)%blk_ptr(idx) = ptr

      end do  ! seg

      ! Record end of last segment
      library(num_library)%seg_info(s_idx,pub_total_num_nodes) = idx
      library(num_library)%seg_info(s_ptr,pub_total_num_nodes) = ptr

      ! Deallocate holding array to store overlap details
      deallocate(ovlap_list,stat=ierr)
      call utils_dealloc_check('internal_index_product (sparse_create)', &
           'ovlap_list',ierr)

      ! Do not exit until all sends are completed
      call comms_free

      my_nat = pub_total_num_nodes - &
           count(library(num_library)%seg_info(s_type,:)==SEG_BLANK)

      iovlap = (library(num_library)%my_nblks + 1)*my_nat + &
           library(num_library)%my_nzb
      if (iovlap+1 /= idx) then
          write(stdout,'(a,i5)') 'Error in internal_index_product: Consistency &
               & check failed on node ', pub_my_node_id
          write(stdout,'(a,i10,a,i10)') 'iovlap = ',iovlap, &
               ', final value of idx-1 = ',idx-1
          call comms_abort
      end if

    end subroutine internal_index_product

  end subroutine sparse_create


  !============================================================================!
  ! This subroutine destroys a sparse matrix structure, freeing up the memory  !
  ! allocated to internal arrays.                                              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   old_spam (inout) : The sparse matrix to be destroyed                     !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Revised for distributed data, June 2006.                                   !
  ! Updated for SPAM3 by Nicholas Hine, May 2009.                              !
  !============================================================================!

  subroutine sparse_destroy(old_spam)

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: old_spam

    ! Deallocate memory assigned to internal structures
    call sparse_data_dealloc(old_spam)

  end subroutine sparse_destroy


  !============================================================================!
  ! This function returns the index length of a sparse matrix.                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, June 2006.                                        !
  !============================================================================!

  integer function sparse_index_length(mat)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Local variable
    integer :: ilib                  ! Library entry for mat

    ! Very simple...
    ilib = mat%lib
    sparse_index_length = library(ilib)%my_nblks + library(ilib)%my_nzb + 1

  end function sparse_index_length

  !============================================================================!
  ! This subroutine returns the index of a sparse matrix in a column-indexed,  !
  ! non-segmented form suitable for use with integrals_mod routines.           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   idx   (output) : The array for the index                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                        !
  !============================================================================!

  subroutine sparse_generate_index(idx,mat)

    use comms, only: comms_abort, pub_on_root, pub_total_num_nodes

    implicit none

    ! Arguments
    integer, intent(out) :: idx(:)  ! The array for the index
    type(SPAM3), intent(in) :: mat  ! The sparse matrix

    ! Local variables
    integer :: ilib                 ! The library entry for mat
    integer :: ilen                 ! The length of the index
    integer :: seg                  ! Segment index
    integer :: seg_start            ! Start position of seg in index
    integer :: iovlap               ! Block overlap counter
    integer :: iblk                 ! Block-column counter
    integer :: loc_iblk             ! Block-column counter local to this node
    integer :: jblk                 ! Block-row counter
    integer :: iidx                 ! Indices

    ! Basic information about the matrix mat
    ilib = mat%lib
    ilen = library(ilib)%my_nblks + library(ilib)%my_nzb + 1

    ! Check array is sufficiently long
    if (size(idx) < ilen) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_generate_index: idx is not long enough for index.'
       call comms_abort
    end if

    ! Index is sorted into node-node segments but we want the whole index
    ! for a single node's columns combined, so we need to reorganise it

    ! Loop over blocks on this node
    iidx = my_last_blk - my_first_blk + 3
    loc_iblk = 0
    do iblk=my_first_blk,my_last_blk
       loc_iblk = loc_iblk + 1

       ! Store index counter in the initial list of column start positions
       idx(loc_iblk) = iidx

       ! Loop over node-segments of the index
       do seg=0,pub_total_num_nodes-1

          ! Nothing to add if this segment is blank
          if (library(ilib)%seg_info(s_type,seg)==SEG_BLANK) cycle

          ! Loop over overlaps for this atom in this segment
          seg_start = library(ilib)%seg_info(s_idx,seg)
          do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
               library(ilib)%blk_idx(seg_start+loc_iblk)-1

             ! Store this block in the new list
             jblk = library(ilib)%blk_idx(iovlap)
             idx(iidx) = jblk
             iidx = iidx + 1

          end do ! iovlap

       end do  ! seg

    end do  ! iblk
    idx(loc_iblk+1) = iidx

  end subroutine sparse_generate_index

  !============================================================================!
  ! This subroutine returns the structure code which corresponds to the        !
  ! transpose of the sparse matrix passed in. For symmetric structures this    !
  ! is the same as the structure of the matrix itself, but for rectangular     !
  ! matrices it may be different.                                              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be transposed                             !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, March 2011                                       !
  !============================================================================!

  subroutine sparse_transpose_structure(struc_trans,mat)

    implicit none

    ! Arguments
    character(len=10), intent(out) :: struc_trans
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Local variable
    integer :: ilib                  ! Library entry for mat

    ! Get transpose strcture
    ilib = mat%lib
    struc_trans = library(ilib)%transpose_structure

  end subroutine sparse_transpose_structure

  !============================================================================!
  ! Returns the index of the first element on a given node in the row or       !
  ! column blocking scheme of a given matrix.                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   node   (input) : The node for which the first element is required        !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_first_elem_on_node(node,mat,rowcol)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix
    integer, intent(in) :: node        ! The node
    character, intent(in) :: rowcol    ! Whether to use row or col blocks

    if (rowcol == 'R') then
       sparse_first_elem_on_node = &
            first_elem_on_node(node,library(mat%lib)%row_blks)
    else if (rowcol == 'C') then
       sparse_first_elem_on_node = &
            first_elem_on_node(node,library(mat%lib)%col_blks)
    else
       sparse_first_elem_on_node = 0
    end if

  end function sparse_first_elem_on_node

  !============================================================================!
  ! Returns the index of the first element on a given atom in the row or       !
  ! column blocking scheme of a given matrix.                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   atom   (input) : The atom for which the first element is required        !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_first_elem_on_atom(atom,mat,rowcol)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix
    integer, intent(in) :: atom        ! The atom
    character, intent(in) :: rowcol    ! Whether to use row or col blocks

    if (rowcol == 'R') then
       sparse_first_elem_on_atom = &
            first_elem_on_atom(atom,library(mat%lib)%row_blks)
    else if (rowcol == 'C') then
       sparse_first_elem_on_atom = &
            first_elem_on_atom(atom,library(mat%lib)%col_blks)
    else
       sparse_first_elem_on_atom = 0
    end if

  end function sparse_first_elem_on_atom

  !============================================================================!
  ! Returns the number of elements on a given node in the row or column        !
  ! blocking scheme of a given matrix.                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   node   (input) : The node for which the number of elements is required   !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_num_elems_on_node(node,mat,rowcol)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix
    integer, intent(in) :: node        ! The node
    character, intent(in) :: rowcol    ! Whether to use row or col blocks

    if (rowcol == 'R') then
       sparse_num_elems_on_node = &
            num_elems_on_node(node,library(mat%lib)%row_blks)
    else if (rowcol == 'C') then
       sparse_num_elems_on_node = &
            num_elems_on_node(node,library(mat%lib)%col_blks)
    else
       sparse_num_elems_on_node = 0
    end if

  end function sparse_num_elems_on_node

  !============================================================================!
  ! Returns the number of elements on a given atom in the row or column        !
  ! blocking scheme of a given matrix.                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   atom   (input) : The atom for which the number of elements is required   !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_num_elems_on_atom_mat(atom,mat,rowcol)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix
    integer, intent(in) :: atom        ! The atom
    character, intent(in) :: rowcol    ! Whether to use row or col blocks

    if (rowcol == 'R') then
       sparse_num_elems_on_atom_mat = &
            num_elems_on_atom(atom,library(mat%lib)%row_blks)
    else if (rowcol == 'C') then
       sparse_num_elems_on_atom_mat = &
            num_elems_on_atom(atom,library(mat%lib)%col_blks)
    else
       sparse_num_elems_on_atom_mat = 0
    end if

  end function sparse_num_elems_on_atom_mat

  !============================================================================!
  ! Returns the number of elements on a given atom in the row or column        !
  ! blocking scheme of a given matrix.                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   atom   (input) : The atom for which the number of elements is required   !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_num_elems_on_atom_blks(atom,blks)

    implicit none

    ! Arguments
    integer, intent(in) :: atom        ! The atom
    integer, intent(in) :: blks        ! The blocking scheme

    sparse_num_elems_on_atom_blks = num_elems_on_atom(atom,blks)

  end function sparse_num_elems_on_atom_blks

  !============================================================================!
  ! Returns the node index of a given element in the row or column blocking    !
  ! scheme of a given matrix.                                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   elem   (input) : The element index required                              !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_node_of_elem(elem,mat,rowcol)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix
    integer, intent(in) :: elem        ! The element
    character, intent(in) :: rowcol    ! Whether to use row or col blocks

    if (rowcol == 'R') then
       sparse_node_of_elem = &
            node_of_elem(elem,library(mat%lib)%row_blks)
    else if (rowcol == 'C') then
       sparse_node_of_elem = &
            node_of_elem(elem,library(mat%lib)%col_blks)
    else
       sparse_node_of_elem = 0
    end if

  end function sparse_node_of_elem

  !============================================================================!
  ! Returns the atom index of a given element in the row or column blocking    !
  ! scheme of a given matrix.                                                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (input) : The matrix                                              !
  !   elem   (input) : The element index required                              !
  !   rowcol (input) : 'R' to use row blocks, 'C' to use column blocks         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_atom_of_elem(elem,mat,rowcol)

    implicit none

    ! Arguments
    integer, intent(in) :: elem        ! The element
    type(SPAM3), intent(in) :: mat     ! The matrix
    character, intent(in) :: rowcol      ! Whether to use row or col blocks

    if (rowcol == 'R') then
       sparse_atom_of_elem = &
            atom_of_elem(elem,library(mat%lib)%row_blks)
    else if (rowcol == 'C') then
       sparse_atom_of_elem = &
            atom_of_elem(elem,library(mat%lib)%col_blks)
    else
       sparse_atom_of_elem = 0
    end if

  end function sparse_atom_of_elem

  !============================================================================!
  ! Tests whether a specified element is listed in the index of a sparse       !
  ! matrix.                                                                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat  (input)   : The matrix                                              !
  !   irow (input)   : The element index of the row required                   !
  !   jcol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, July 2006.                                        !
  ! Modified for SPAM3 by Nicholas Hine, May 2009.                             !
  !============================================================================!

  logical function sparse_element_exists(mat,irow,jcol)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix
    integer, intent(in) :: irow        ! The element index of the row
    integer, intent(in) :: jcol        ! The element index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: iblk      ! The block-row to which the element belongs
    integer :: jblk      ! The block-col to which the element belongs
    integer :: kblk      ! Block counter
    integer :: iovlap    ! Block overlap counter
    integer :: loc_jblk  ! Block-column counter local to this node
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index

    ! Get library entry
    ilib = mat%lib

    ! Check arguments
#ifdef DEBUG
    if (irow < 1 .or. irow > library(ilib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_element_exists: &
            &invalid row index'
       call comms_abort
    end if
    if (jcol < 1 .or. jcol > library(ilib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_element_exists: &
            &invalid column index'
       call comms_abort
    end if
#endif

    ! Obtain information about the relevant block for this element
    seg = node_of_elem(irow,library(ilib)%row_blks)
    iblk = atom_of_elem(irow,library(ilib)%row_blks)
    jblk = atom_of_elem(jcol,library(ilib)%col_blks)

    ! Check this block-column is local to this processor
    if (jblk < my_first_blk .or. jblk > my_last_blk) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_element_exists: &
            &requested column is not local to this node'
       call comms_abort
    end if

    ! Assume not found
    sparse_element_exists = .false.

    ! Element definitely doesn't exist if segment is blank
    if (library(ilib)%seg_info(s_type,seg) == SEG_BLANK) return

    ! Find local index for the block-column
    loc_jblk = jblk - my_first_blk + 1

    ! Loop over the nonzero blocks in this segment
    seg_start = library(ilib)%seg_info(s_idx,seg)
    do iovlap=library(ilib)%blk_idx(seg_start+loc_jblk-1), &
         library(ilib)%blk_idx(seg_start+loc_jblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk == iblk) then
          sparse_element_exists = .true.
          exit
       end if
    end do

  end function sparse_element_exists


  !============================================================================!
  ! This function returns true the denominator for the filling factor of a     !
  ! sparse matrix, based on the row and column blocking schemes.               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   row_blks   (input)  : The row blocking scheme.                           !
  !   col_blks   (input)  : The col blocking scheme.                           !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in April 2011.                                    !
  !============================================================================!

  real(kind=DP) function sparse_fill_fac_denom(row_blks,col_blks)

    implicit none

    ! Arguments
    integer, intent(in) :: row_blks
    integer, intent(in) :: col_blks

    sparse_fill_fac_denom = 100.0_DP / (&
         real(sum(num_elems_on_node(:,row_blks)),kind=DP) * &
         real(sum(num_elems_on_node(:,col_blks)),kind=DP))

  end function sparse_fill_fac_denom


  !============================================================================!
  ! This function returns true if every element in the matrix is nonzero, or   !
  ! false otherwise.                                                           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in December 2009.                                 !
  !============================================================================!

  logical function sparse_is_dense(mat)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Local Variables
    integer(kind=LONG) :: num_sq

    ! Check if the number of nonzero elements is equal to nrows*mcols
    num_sq = int(library(mat%lib)%nrows,kind=LONG) * &
         int(library(mat%lib)%mcols,kind=LONG)

    ! Return true if so, false otherwise
    if (library(mat%lib)%nze==num_sq) then
       sparse_is_dense = .true.
    else
       sparse_is_dense = .false.
    end if

  end function sparse_is_dense


  !============================================================================!
  ! This function calculates the RMS element value of a sparse matrix.         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, June 2006.                                        !
  !============================================================================!

  real(kind=DP) function sparse_rms_element(mat)

    use comms, only: comms_reduce

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Local variables
    integer :: ilib            ! Library entry for mat
    integer :: iel             ! Element loop counter
    real(kind=DP) :: loc_sum   ! Sum accumulated on this processor
    real(kind=DP) :: glob_sum  ! Sum accumulated across all processors

    ! Get library entry for mat
    ilib = mat%lib

    ! Sum modulus squared of each element on this processor
    loc_sum = 0.0_DP
    if (mat%iscmplx) then
       do iel=1,library(ilib)%my_nze
          loc_sum = loc_sum + real(mat%zmtx(iel)*conjg(mat%zmtx(iel)),kind=DP)
       end do
    else
       do iel=1,library(ilib)%my_nze
          loc_sum = loc_sum + mat%dmtx(iel)*mat%dmtx(iel)
       end do
    end if

    ! Sum up across processors
    glob_sum = loc_sum
    call comms_reduce('SUM',glob_sum)

    ! Calculate RMS value
    sparse_rms_element = sqrt(glob_sum / real(library(ilib)%nze,kind=DP))

  end function sparse_rms_element


  !============================================================================!
  ! This function calculates the maximum absolute element value of a sparse    !
  ! matrix.                                                                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, June 2006.                                        !
  !============================================================================!

  real(kind=DP) function sparse_max_abs_element(mat)

    use comms, only: comms_reduce

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Local variables
    integer :: ilib            ! Library entry for mat
    integer :: iel             ! Element loop counter
    real(kind=DP) :: loc_max   ! Maximum value on this processor
    real(kind=DP) :: glob_max  ! Maximum value across all processors

    ! Get library entry for mat
    ilib = mat%lib

    ! Find maximum absolute value of every element on this processor
    loc_max = 0.0_DP
    if (mat%iscmplx) then
       do iel=1,library(ilib)%my_nze
          loc_max = max(loc_max, &
               real(mat%zmtx(iel) * conjg(mat%zmtx(iel)),kind=DP))
       end do
       loc_max = sqrt(loc_max)
    else
       do iel=1,library(ilib)%my_nze
          loc_max = max(loc_max, abs(mat%dmtx(iel)))
       end do
    end if

    !  across processors
    glob_max = loc_max
    call comms_reduce('MAX',glob_max)

    ! Calculate RMS value
    sparse_max_abs_element = glob_max

  end function sparse_max_abs_element


  !============================================================================!
  ! This function returns the number of non-zero elements in a sparse matrix   !
  ! in the columns stored on this node.                                        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  integer function sparse_node_num_element(mat)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Very simple...
    sparse_node_num_element = library(mat%lib)%my_nze

  end function sparse_node_num_element


  !============================================================================!
  ! This function returns the number of non-zero elements in a sparse matrix   !
  ! specified by its library entry, in a floating-point real.                  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, April 2011.                                      !
  !============================================================================!

  real(kind=DP) function sparse_num_element_lib(ilib)

    implicit none

    ! Arguments
    integer, intent(in) :: ilib   ! The library entry

    ! Very simple...
    sparse_num_element_lib = real(library(ilib)%nze,kind=DP)

  end function sparse_num_element_lib


  !============================================================================!
  ! This function returns the number of non-zero elements in a sparse matrix   !
  ! in a floating-point real (so as to avoid integer overflows).               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, December 2009, based on sparse_num_element,      !
  ! originally written by Peter Haynes, June 2006.                             !
  !============================================================================!

  real(kind=DP) function sparse_num_element_mat(mat)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Very simple...
    sparse_num_element_mat = real(library(mat%lib)%nze,kind=DP)

  end function sparse_num_element_mat


  !============================================================================!
  ! This function returns the number of rows in a sparse matrix by consulting  !
  ! the relevant library entry.                                                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, January 2010.                                    !
  !============================================================================!

  real(kind=DP) function sparse_num_rows(mat)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Very simple...
    sparse_num_rows = library(mat%lib)%nrows

  end function sparse_num_rows


  !============================================================================!
  ! This function returns the number of cols in a sparse matrix by consulting  !
  ! the relevant library entry.                                                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, January 2010.                                    !
  !============================================================================!

  real(kind=DP) function sparse_num_cols(mat)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat   ! The sparse matrix

    ! Very simple...
    sparse_num_cols = library(mat%lib)%mcols

  end function sparse_num_cols


  !============================================================================!
  ! This function returns the number of non-zero elements in a sparse matrix   !
  ! in the columns stored on this node.                                        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat   (input)  : The matrix to be assessed                               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, February 2010.                                   !
  !============================================================================!

  subroutine sparse_memory(local_mem,global_mem,structure)

    use comms, only: comms_abort, pub_on_root, pub_total_num_nodes
    use constants, only: char_size, int_size, logical_size, real_size

    implicit none

    ! Arguments
    integer, intent(out) :: local_mem
    integer, intent(out) :: global_mem
    character(*), intent(in) :: structure   ! Sparse matrix structure code

    ! Local Variables
    integer :: ilib

    call sparse_search_library(structure,ilib)

    ! Quit if not found
    if (ilib == 0) then
       if (pub_on_root) write(stdout,'(2a)') 'Error in sparse_memory: &
            &no library structure found to match ',trim(structure)
       call comms_abort
    end if

    local_mem = library(ilib)%my_nze
    global_mem = library(ilib)%nze

  end subroutine sparse_memory


  !============================================================================!
  ! This function displays the currently allocated memory of all the sparse    !
  ! matrices currently in existence.                                           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, February 2010.                                   !
  !============================================================================!

  subroutine sparse_show_memory_usage

    use comms, only: comms_abort, pub_on_root, pub_total_num_nodes
    use constants, only: char_size, int_size, logical_size, real_size

    implicit none

    if (pub_on_root) then
       write(stdout,'(a)') 'Currently allocated Sparse Matrix Memory'
       write(stdout,'(a,i20)') 'Global           :',global_mat_mem
       write(stdout,'(a,i20)') 'Local (root node):',local_mat_mem
    end if

  end subroutine sparse_show_memory_usage

  !============================================================================!
  ! This subroutine displays the elements of the matrix smat for debugging     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   smat  (inout) : The matrix to print                                      !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                        !
  !============================================================================!

 !subroutine sparse_show_matrix(smat,outunit,show_elems,matlab_format)
!CW
  subroutine sparse_show_matrix(smat,outunit,show_elems,matlab_format,labels)
!END CW

    use comms, only: comms_abort, comms_barrier, comms_bcast, comms_send, &
         comms_recv, pub_on_root, pub_my_node_id, pub_root_node_id, &
         pub_total_num_nodes, comms_free
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_elements_on_node, pub_first_atom_on_node, &
         pub_orig_atom
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_flush

    implicit none

    ! Arguments
    type(SPAM3),intent(in) :: smat
    integer, intent(in), optional :: outunit
    logical, intent(in), optional :: show_elems
    logical, intent(in), optional :: matlab_format

    ! Locals
    integer :: loc_outunit        !
    logical :: loc_show_elems     !
    logical :: loc_matlab_format  !
    integer :: ierr               !
    integer :: seg                !
    integer :: node               !
    integer :: nrows              ! Total number of rows in this block/segment
    integer :: mcols              !
    integer :: datlen             !
    integer :: ptr                !
    integer :: iorig              !
    integer :: row_start, row_end !
    integer :: row_len            !
    integer :: col_start, col_end !
    integer :: col_len            !
    integer :: col_blks           ! Identifier for col blocking scheme
    integer :: row_blks           ! Identifier for row blocking scheme
    integer :: ielem,global_ielem !
    integer :: iatom              !
    integer :: loc_iatom          !
    integer :: max_seg_rows       !
    integer :: predpwidth         ! Pre decimal point width
    integer :: width              ! Width of number for display
    real(kind=DP), allocatable :: dmtxbuf(:)    ! Buffer for real data
    complex(kind=DP), allocatable :: zmtxbuf(:) ! Buffer for complex data
    real(kind=DP), allocatable :: dmtxout(:)    !
    complex(kind=DP), allocatable :: zmtxout(:) !
    character(len=1000) :: fmt                  !
    character(len=12) :: tmp

    character(len=2),allocatable :: row_symbols(:) !
    character(len=2),allocatable :: col_symbols(:) !
    integer,allocatable :: row_elem_on_atom(:)     !
    integer,allocatable :: col_elem_on_atom(:)     !
!CW
    character(len=2),optional :: labels(library(smat%lib)%mcols)
!END CW

    ! Process optional arguments
    if (present(outunit)) then
       loc_outunit = outunit
    else
       loc_outunit = 6
    end if
    if (present(show_elems)) then
       loc_show_elems = show_elems
    else
       loc_show_elems = .false.
    end if
    if (present(matlab_format)) then
       loc_matlab_format = matlab_format
    else
       loc_matlab_format = .false.
    end if

    ! Consistency Check
    if (loc_matlab_format.and.loc_show_elems) then
       ! error
       call comms_abort
    end if

    mcols = library(smat%lib)%mcols
    nrows = library(smat%lib)%nrows
    col_blks = library(smat%lib)%col_blks
    row_blks = library(smat%lib)%row_blks
    max_seg_rows = maxval(num_elems_on_node(:,row_blks))

    ! Allocate communication buffer
    if (smat%iscmplx) then
       allocate(zmtxbuf(mcols*max_seg_rows),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','zmtxbuf',ierr)
       allocate(zmtxout(mcols),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','zmtxout',ierr)
    else
       allocate(dmtxbuf(mcols*max_seg_rows),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','dmtxbuf',ierr)
       allocate(dmtxout(mcols),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','dmtxout',ierr)
    end if

    ! Allocate storage for symbols and elem_on_atom lists
    if (loc_show_elems) then

       allocate(col_symbols(mcols),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','col_symbols',ierr)
       allocate(row_symbols(nrows),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','row_symbols',ierr)
       allocate(col_elem_on_atom(mcols),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','col_elem_on_atom',ierr)
       allocate(row_elem_on_atom(nrows),stat=ierr)
       call utils_alloc_check('sparse_show_matrix','row_elem_on_atom',ierr)

       ! Collect this node's symbol strings and element indices along cols
       do ielem=first_elem_on_node(pub_my_node_id,col_blks), &
            first_elem_on_node(pub_my_node_id+1,col_blks)-1
         iatom = atom_of_elem(ielem,col_blks)
         loc_iatom = iatom - pub_first_atom_on_node(pub_my_node_id) + 1
         col_symbols(ielem) = pub_elements_on_node(loc_iatom)%symbol
         col_elem_on_atom(ielem) = ielem - first_elem_on_atom(iatom,col_blks)+1
       end do

       ! Collect this node's symbol strings and element indices along rows
       do ielem=first_elem_on_node(pub_my_node_id,row_blks), &
            first_elem_on_node(pub_my_node_id+1,row_blks)-1
         iatom = atom_of_elem(ielem,row_blks)
         loc_iatom = iatom - pub_first_atom_on_node(pub_my_node_id) + 1
         row_symbols(ielem) = pub_elements_on_node(loc_iatom)%symbol
         row_elem_on_atom(ielem) = ielem - first_elem_on_atom(iatom,row_blks)+1
       end do

       ! Share symbol strings and element indices
       do node=0,pub_total_num_nodes-1
          col_start = first_elem_on_node(node,col_blks)
          col_end = first_elem_on_node(node+1,col_blks) - 1
          col_len = col_end - col_start + 1
          row_start = first_elem_on_node(node,row_blks)
          row_end = first_elem_on_node(node+1,row_blks) - 1
          row_len = row_end - row_start + 1
          call comms_bcast(node,col_symbols(col_start:col_end),col_len*2)
          call comms_bcast(node,row_symbols(row_start:row_end),row_len*2)
          call comms_bcast(node,col_elem_on_atom(col_start:col_end),col_len)
          call comms_bcast(node,row_elem_on_atom(row_start:row_end),row_len)
       end do
!CW
       if(present(labels)) then
          labels=col_symbols
       endif
!END CW
    end if

    ! Only the root node shows the matrix
    if (pub_on_root.and.loc_show_elems) then
       ! Write the line identifying the atoms and elements on each atom
       write(loc_outunit,"(12x)",advance='no')
       do ielem=1,mcols
          iorig = pub_orig_atom(atom_of_elem(ielem,col_blks))
          if (smat%iscmplx) then
             write(loc_outunit,"(15x,a2,i4,2x,i2,1x)",advance='no') &
                  trim(col_symbols(ielem)),iorig,col_elem_on_atom(ielem)
          else
             write(loc_outunit,"(7x,a2,i4,2x,i2,1x)",advance='no') &
                  trim(col_symbols(ielem)),iorig,col_elem_on_atom(ielem)
          end if
       end do
       write(loc_outunit,*)
    end if

    ! qoh: Find width neccessary to show widest number in matrix
!CW
    if(sparse_max_abs_element(smat)/=0.0_DP)then
!END CW
     predpwidth = ceiling(log10(sparse_max_abs_element(smat)))
     predpwidth = max(predpwidth,1)
!CW
    else
     predpwidth=1
    endif
!END CW 

    if (pub_on_root) then
       ! Write the format string to use for the main data
       write(tmp,'(i10)') mcols
       if (smat%iscmplx) then
!CW
         ! if (loc_matlab_format) then
         !    write(fmt,"(a,a,a)")'(sp,',trim(adjustl(tmp)), &
         !         "(e16.10,e16.10,'i',x),';')"
         ! else if (loc_show_elems) then
         !    write(fmt,"(a,a,a)")'(a2,i4,i4,2x,sp,',trim(adjustl(tmp)), &
         !         "(f12.8,f12.8,'i',x))"
         ! else
         !    write(fmt,"(a,a,a)")'(sp,',trim(adjustl(tmp)), &
         !         "(f20.14,f20.14,'i',x))"
         ! end if
!END CW
!CW
          if (loc_matlab_format) then
!             write(fmt,"(a,a,a)")'(sp,',trim(adjustl(tmp)),  "(e16.10,e16.10,'i',x),';')"
              write(fmt,"(a,a,a)")'(sp,',trim(adjustl(tmp)), "(e5.2,e5.2,'i',x),';')"
          else if (loc_show_elems) then
!             write(fmt,"(a,a,a)")'(a2,i4,i4,2x,sp,',trim(adjustl(tmp)), "(f12.8,f12.8,'i',x))"
              write(fmt,"(a,a,a)")'(a2,i4,i4,2x,sp,',trim(adjustl(tmp)), "(f5.2,f5.2,'i',x))"
          else
!             write(fmt,"(a,a,a)")'(sp,',trim(adjustl(tmp)), "(f20.14,f20.14,'i',x))"
              write(fmt,"(a,a,a)")'(sp,',trim(adjustl(tmp)), "(f7.3,f7.3,'i',x))"
          end if
!END CW
       else
          if (loc_matlab_format) then
             write(fmt,"(a,a,a)")'(',trim(adjustl(tmp)),"f20.12,';')"
          else if (loc_show_elems) then
             write(fmt,"(a,a,a)")'(a2,i4,i4,2x,',trim(adjustl(tmp)),'f18.10)'
          else
             width = 17 + predpwidth
             write(fmt,"(a,a,a,i2,a)")'(',trim(adjustl(tmp)),'f',width,'.14)'
          end if
       end if
    end if

    ! Loop over segments
    if (pub_on_root.and.loc_matlab_format) &
         write(loc_outunit,'(a)',advance='no') '['
    do seg=0,pub_total_num_nodes-1
       call comms_free

       ! Re-create the dense version of this segment if required
       if (smat%iscmplx) then
          call internal_dense_segment_complex(smat,seg,zmtxbuf)
       else
          call internal_dense_segment_real(smat,seg,dmtxbuf)
       end if

       nrows = num_elems_on_node(seg,row_blks)
       datlen = nrows * num_elems_on_node(pub_my_node_id,col_blks)
       ptr = 1
       do node=0,pub_total_num_nodes-1
          datlen = nrows * num_elems_on_node(node,col_blks)
          if (pub_on_root.and.(node/=pub_my_node_id)) then
             if (smat%iscmplx) then
                call comms_recv(node,zmtxbuf(ptr:ptr+datlen-1),datlen)
             else
                call comms_recv(node,dmtxbuf(ptr:ptr+datlen-1),datlen)
             end if
          else if (.not.(pub_on_root).and.(node==pub_my_node_id)) then
             if (smat%iscmplx) then
                call comms_send(pub_root_node_id,zmtxbuf(1:datlen),datlen)
             else
                call comms_send(pub_root_node_id,dmtxbuf(1:datlen),datlen)
             end if
          end if
          ptr = ptr + datlen
       end do

       if (pub_on_root) then

          ! Loop over rows in this segment
          do ielem=1,num_elems_on_node(seg,row_blks)
             global_ielem = ielem + first_elem_on_node(seg,row_blks) - 1
             iorig = pub_orig_atom(atom_of_elem(global_ielem,row_blks))
             if (smat%iscmplx) then
                ptr = ielem
                do node=0,pub_total_num_nodes-1
                   datlen = nrows * num_elems_on_node(node,col_blks)
                   zmtxout(first_elem_on_node(node,col_blks): &
                        first_elem_on_node(node+1,col_blks)-1) = &
                        zmtxbuf(ptr:ptr+datlen-ielem:nrows)
                   ptr = ptr + datlen
                end do
                if (loc_show_elems) then
                   write(loc_outunit,fmt) trim(row_symbols(global_ielem)), &
                        iorig,row_elem_on_atom(global_ielem),zmtxout(1:mcols)
                else
                   write(loc_outunit,fmt) zmtxout(1:mcols)
                end if
             else
                ptr = ielem
                do node=0,pub_total_num_nodes-1
                   datlen = nrows * num_elems_on_node(node,col_blks)
                   dmtxout(first_elem_on_node(node,col_blks): &
                        first_elem_on_node(node+1,col_blks)-1) = &
                        dmtxbuf(ptr:ptr+datlen-ielem:nrows)
                   ptr = ptr + datlen
                end do
                if (loc_show_elems) then
                   write(loc_outunit,fmt) trim(row_symbols(global_ielem)), &
                        iorig,row_elem_on_atom(global_ielem),dmtxout(1:mcols)
                else
                   write(loc_outunit,fmt) dmtxout(1:mcols)
                end if
             endif
          end do

       end if

       call comms_barrier

    end do
    if (pub_on_root.and.loc_matlab_format) &
         write(loc_outunit,'(a)') ']'

    call utils_flush
    call comms_barrier
    call utils_flush

    ! Dellocate storage for symbols and elem_on_atom list
    if (loc_show_elems) then
       deallocate(row_elem_on_atom,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','row_elem_on_atom',ierr)
       deallocate(col_elem_on_atom,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','col_elem_on_atom',ierr)
       deallocate(row_symbols,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','row_symbols',ierr)
       deallocate(col_symbols,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','col_symbols',ierr)
    end if

    ! Deallocate communication buffers
    if (smat%iscmplx) then
       deallocate(zmtxout,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','zmtxout',ierr)
       deallocate(zmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','zmtxbuf',ierr)
    else
       deallocate(dmtxout,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','dmtxout',ierr)
       deallocate(dmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_show_matrix','dmtxbuf',ierr)
    end if

contains

    !========================================================================!
    ! Collates the data of a real-valued segment of a matrix (sparse or      !
    ! dense) so that it can be printed.                                      !
    !------------------------------------------------------------------------!
    ! Written by Nicholas Hine, May 2009.                                    !
    !========================================================================!

    subroutine internal_dense_segment_real(mat,seg,dbuf)

      implicit none

      ! Arguments
      type(SPAM3),intent(in) :: mat
      integer,intent(in) :: seg
      real(kind=DP), intent(out) :: dbuf(:)

      ! Locals
      integer :: ilib          ! Library entry for mat
      integer :: iblk          ! Block-column counter
      integer :: loc_iblk      ! Block-column counter local to this node
      integer :: jblk          ! Block-row counter
      integer :: iovlap        ! Block overlap counter
      integer :: ielem,jelem   ! Element row/col counters
      integer :: ptr           ! Pointer to data in block
      integer :: datlen        ! Length of data in block
      integer :: bufptr        ! Pointer to data in buffer
      integer :: row_blks      ! Identifier for column blocking scheme
      integer :: col_blks      ! Identifier for row blocking scheme
      integer :: seg_type      ! Segment type for this segment
      integer :: seg_start     ! Start position of this segment in the index

      ! Get library entry
      ilib = mat%lib
      row_blks = library(ilib)%row_blks
      col_blks = library(ilib)%col_blks

      seg_start = library(ilib)%seg_info(s_idx,seg)
      seg_type = library(ilib)%seg_info(s_type,seg)

      datlen = num_elems_on_node(pub_my_node_id,col_blks) * &
           num_elems_on_node(seg,row_blks)
      dbuf(1:datlen) = 0.0_DP

      if (seg_type == SEG_DENSE) then

         ptr = library(ilib)%seg_info(s_ptr,seg)
         datlen = library(ilib)%seg_info(s_ptr,seg + 1) - ptr
         dbuf(1:datlen) = mat%dmtx(ptr:ptr+datlen-1)

      else if (seg_type == SEG_SPARSE) then

         ! Loop over atom block-columns iblk on this node
         loc_iblk = 0
         seg_start = library(ilib)%seg_info(s_idx,seg)
         do iblk=my_first_blk,my_last_blk
            loc_iblk = loc_iblk + 1

            ! Loop over atoms jblk overlapping atom iblk: block-rows in
            ! block-col iblk of spam
            do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
                 library(ilib)%blk_idx(seg_start+loc_iblk)-1
               jblk = library(ilib)%blk_idx(iovlap)
               ptr = library(ilib)%blk_ptr(iovlap)

               bufptr = (first_elem_on_atom(iblk,col_blks) - &
                    first_elem_on_node(pub_my_node_id,col_blks))* &
                    num_elems_on_node(seg,row_blks) + &
                    first_elem_on_atom(jblk,row_blks) - &
                    first_elem_on_node(seg,row_blks) + 1
               do ielem=1,num_elems_on_atom(iblk,col_blks)
                  do jelem=1,num_elems_on_atom(jblk,row_blks)
                     dbuf(bufptr) = mat%dmtx(ptr)
                     bufptr = bufptr + 1
                     ptr = ptr + 1
                  end do ! jelem
                  bufptr = bufptr - num_elems_on_atom(jblk,row_blks)
                  bufptr =  bufptr + num_elems_on_node(seg,row_blks)
               end do ! ielem

            end do ! iovlap

         end do ! iblk

      end if ! seg_type

    end subroutine internal_dense_segment_real

    !========================================================================!
    ! Collates the data of a complex-valued segment of a matrix (sparse or   !
    ! dense) so that it can be printed.                                      !
    !------------------------------------------------------------------------!
    ! Written by Nicholas Hine, May 2009.                                    !
    !========================================================================!

    subroutine internal_dense_segment_complex(mat,seg,zbuf)

      implicit none

      ! Arguments
      type(SPAM3),intent(in) :: mat
      integer,intent(in) :: seg
      complex(kind=DP), intent(out) :: zbuf(:)

      ! Locals
      integer :: ilib          ! Library entry for mat
      integer :: iblk          ! Block-column counter
      integer :: loc_iblk      ! Block-column counter local to this node
      integer :: jblk          ! Block-row counter
      integer :: iovlap        ! Block overlap counter
      integer :: ielem,jelem   ! Element row/col counters
      integer :: ptr           ! Pointer to data in block
      integer :: datlen        ! Length of data in block
      integer :: bufptr        ! Pointer to data in buffer
      integer :: row_blks      ! Identifier for column blocking scheme
      integer :: col_blks      ! Identifier for row blocking scheme
      integer :: seg_type      ! Segment type for this segment
      integer :: seg_start     ! Start position of this segment in the index

      ! Get library entry
      ilib = mat%lib
      row_blks = library(ilib)%row_blks
      col_blks = library(ilib)%col_blks

      seg_start = library(ilib)%seg_info(s_idx,seg)
      seg_type = library(ilib)%seg_info(s_type,seg)

      datlen = num_elems_on_node(pub_my_node_id,col_blks) * &
           num_elems_on_node(seg,row_blks)
      zbuf(1:datlen) = (0.0_DP,0.0_DP)

      if (seg_type == SEG_DENSE) then

         ptr = library(ilib)%seg_info(s_ptr,seg)
         datlen = library(ilib)%seg_info(s_ptr,seg + 1) - ptr
         zbuf(1:datlen) = mat%zmtx(ptr:ptr+datlen-1)

      else if (seg_type == SEG_SPARSE) then

         ! Loop over atom block-columns iblk on this node
         loc_iblk = 0
         seg_start = library(ilib)%seg_info(s_idx,seg)
         do iblk=my_first_blk,my_last_blk
            loc_iblk = loc_iblk + 1

            ! Loop over atoms jblk overlapping atom iblk: block-rows in
            ! block-col iblk of spam
            do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
                 library(ilib)%blk_idx(seg_start+loc_iblk)-1
               jblk = library(ilib)%blk_idx(iovlap)
               ptr = library(ilib)%blk_ptr(iovlap)

               bufptr = (first_elem_on_atom(iblk,col_blks) - &
                    first_elem_on_node(pub_my_node_id,col_blks))* &
                    num_elems_on_node(seg,row_blks) + &
                    first_elem_on_atom(jblk,row_blks) - &
                    first_elem_on_node(seg,row_blks) + 1
               do ielem=1,num_elems_on_atom(iblk,col_blks)
                  do jelem=1,num_elems_on_atom(jblk,row_blks)
                     zbuf(bufptr) = mat%zmtx(ptr)
                     bufptr = bufptr + 1
                     ptr = ptr + 1
                  end do ! jelem
                  bufptr = bufptr - num_elems_on_atom(jblk,row_blks)
                  bufptr =  bufptr + num_elems_on_node(seg,row_blks)
               end do ! ielem

            end do ! iovlap

         end do ! iblk

      end if ! seg_type

    end subroutine internal_dense_segment_complex

  end subroutine sparse_show_matrix


  subroutine sparse_show_network(mat,iunit)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id, &
         pub_total_num_nodes, pub_root_node_id, comms_barrier
    use parallel_strategy, only: pub_first_atom_on_node, &
         pub_num_atoms_on_node

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: iunit

    ! Local Variables
    type(COM3) :: com
    integer :: recvnode
    integer :: iblk,loc_iblk
    integer :: jblk,loc_jblk
    integer :: kblk
    integer :: jidx
    integer :: seg, seg_start
    integer :: nrow,ncol
    logical :: found

    ! Allocate arrays and initialise comms
    call sparse_com_allocate(com,mat,1, &
         alloc_mtx=.false.,cropped=.false.,seg=.false.)

    do recvnode=0,pub_total_num_nodes-1

       call comms_barrier

       ! Check if this node needs to send yet
       if ((recvnode==pub_my_node_id).and.(.not.pub_on_root)) then
          call sparse_send_index(com,pub_root_node_id,.false.)
       end if

       ! Root node receives and writes
       if (pub_on_root) then
          call sparse_recv_index(com,recvnode,1,.false.,async=.false.)

          do iblk=pub_first_atom_on_node(recvnode), &
               pub_first_atom_on_node(recvnode+1)-1
             loc_iblk = iblk - pub_first_atom_on_node(recvnode) + 1
             ncol = num_elems_on_atom(iblk,library(mat%lib)%row_blks)

             do seg=0,pub_total_num_nodes-1
                seg_start = com%seginfobuf(s_idx,seg,1)
                if (com%seginfobuf(s_type,seg,1)==SEG_BLANK) then
                   do loc_jblk=1,pub_num_atoms_on_node(seg)
                      write(iunit,'(i4)',advance='no') 0
                   end do
                   cycle
                end if
                do loc_jblk=1,pub_num_atoms_on_node(seg)
                   jblk = loc_jblk + pub_first_atom_on_node(seg) -1
                   found = .false.
                   ! Loop over block-rows kblk of A in column jblk
                   do jidx=com%idxbuf(seg_start+loc_iblk-1,1), &
                        com%idxbuf(seg_start+loc_iblk,1)-1
                      kblk = com%idxbuf(jidx,1)
                      if (kblk==jblk) then
                         found=.true.
                         exit
                      end if
                   end do
                   if (found) then
                      nrow = num_elems_on_atom(jblk,library(mat%lib)%row_blks)
                      write(iunit,'(i4)',advance='no') nrow*ncol
                   else
                      write(iunit,'(i4)',advance='no') 0
                   end if
                end do
             end do
             write(iunit,*)
          end do
       end if
    end do

    ! Allocate arrays and initialise comms
    call comms_barrier
    call sparse_com_deallocate(com,dealloc_mtx=.false.,cropped=.false.)

  end subroutine sparse_show_network

  !============================================================================!
  ! This subroutine inserts an element into a given block sparse matrix local  !
  ! to this processor.                                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (input)     : The element of the real matrix to insert                !
  !   mat (inout)    : The real matrix                                         !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in June 2006.                       !
  !============================================================================!

  subroutine sparse_put_element_real(el,mat,jrow,icol)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: el           ! The element to insert
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-row to which the element belongs
    integer :: jblk      ! The block-col to which the element belongs
    integer :: kblk      ! Block counter
    integer :: nrows     ! Total number of rows in this block/segment
    integer :: jelem0    ! The first element in the requested block-row
    integer :: ielem0    ! The first element in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks
    nrows = 0

    ! Check arguments
    if (mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_put_element_real: real matrices only'
       call comms_abort
    end if
#ifdef DEBUG
    if (jrow < 1 .or. jrow > library(ilib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_put_element_real: invalid row index'
       call comms_abort
    end if
    if (icol < 1 .or. icol > library(ilib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_put_element_real: invalid column index'
       call comms_abort
    end if
#endif
    ! Check this column is local to this processor
    if (icol<first_elem_on_node(pub_my_node_id,col_blks).or.&
         icol>=first_elem_on_node(pub_my_node_id+1,col_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_put_element_real: requested column is not local to &
            &this node'
       call comms_abort
    end if

    ! Obtain information about the relevant block & segment for this element
    seg = node_of_elem(jrow,row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type==SEG_BLANK) return
    seg_start = library(ilib)%seg_info(s_idx,seg)
    jblk = atom_of_elem(jrow,row_blks)
    iblk = atom_of_elem(icol,col_blks)
    jelem0 = first_elem_on_atom(jblk,row_blks)
    ielem0 = first_elem_on_atom(iblk,col_blks)
    if (seg_type==SEG_DENSE) then
       nrows = num_elems_on_node(seg,row_blks)
    else if (seg_type==SEG_SPARSE) then
       nrows = num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - my_first_blk + 1

    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       mat%dmtx(ptr + (icol - ielem0)*nrows + (jrow - jelem0)) = el
    end do

  end subroutine sparse_put_element_real

  !============================================================================!
  ! This subroutine retrieves an element into a given block sparse matrix      !
  ! local to this processor.                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (output)    : The element of the real matrix retrieved                !
  !   mat (inout)    : The real matrix                                         !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in June 2006.                       !
  !============================================================================!

  subroutine sparse_get_element_real(el,mat,jrow,icol)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: el          ! The element to insert
    type(SPAM3), intent(in) :: mat            ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-row to which the element belongs
    integer :: jblk      ! The block-col to which the element belongs
    integer :: kblk      ! Block counter
    integer :: nrows     ! Total number of rows in this block/segment
    integer :: jelem0    ! The first element in the requested block-row
    integer :: ielem0    ! The first element in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! Check arguments
    if (mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_element_real: &
            &real matrices only'
       call comms_abort
    end if
#ifdef DEBUG
    if (jrow < 1 .or. jrow > library(ilib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_element_real: &
            &invalid row index'
       call comms_abort
    end if
    if (icol < 1 .or. icol > library(ilib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_element_real: &
            &invalid column index'
       call comms_abort
    end if
#endif
    ! Check this column is local to this processor
    if (icol<first_elem_on_node(pub_my_node_id,col_blks).or.&
         icol>=first_elem_on_node(pub_my_node_id+1,col_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_element_real: &
          &requested column is not local to this node'
       call comms_abort
    end if

    el = 0.0_DP
    nrows = 0

    ! Obtain information about the relevant block & segment for this element
    seg = node_of_elem(jrow,row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type==SEG_BLANK) return
    seg_start = library(ilib)%seg_info(s_idx,seg)
    jblk = atom_of_elem(jrow,row_blks)
    iblk = atom_of_elem(icol,col_blks)
    jelem0 = first_elem_on_atom(jblk,row_blks)
    ielem0 = first_elem_on_atom(iblk,col_blks)
    if (seg_type==SEG_DENSE) then
       nrows = num_elems_on_node(seg,row_blks)
    else if (seg_type==SEG_SPARSE) then
       nrows = num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - my_first_blk + 1

    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       el = mat%dmtx(ptr + (icol - ielem0)*nrows + (jrow - jelem0))
    end do

  end subroutine sparse_get_element_real

  !============================================================================!
  ! This subroutine returns a block of a given block sparse matrix local to    !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   blk (output)   : The required block of the real matrix                   !
  !   mat (input)    : The real matrix whose element is required               !
  !   jblk (input)   : The atom index of the row required                      !
  !   iblk (input)   : The atom index of the column required                   !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in July 2006.                       !
  !============================================================================!

  subroutine sparse_get_block_real(blk,mat,jblk,iblk)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id
    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: blk(:,:)    ! The required block
    type(SPAM3), intent(in) :: mat            ! The matrix
    integer, intent(in) :: jblk               ! The atom index of the row
    integer, intent(in) :: iblk               ! The atom index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: kblk      ! Block counter
    integer :: ielem     ! Element row counter
    integer :: jelem     ! Element column counter
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    integer :: nrows     ! Total number of rows in this block/segment

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! Check arguments
    if(ilib == 0) then ! jd: Otherwise uninitialized matrices segfault below
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_block_real: &
            &Matrix does not exist'
       call comms_abort
    end if
    if (mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_block_real: &
            &real matrices only'
       call comms_abort
    end if
#ifdef DEBUG
    if (jblk < 1 .or. jblk > library(ilib)%nblk) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_block_real: &
            &invalid row index'
       call comms_abort
    end if
    if (iblk < 1 .or. iblk > library(ilib)%nblk) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_block_real: &
            &invalid column index'
       call comms_abort
    end if
    if (size(blk,1) < num_elems_on_atom(jblk,row_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_block_real: &
            &too few rows in blk'
       call comms_abort
    end if
    if (size(blk,2) < num_elems_on_atom(iblk,col_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_block_real: &
            &too few columns in blk'
       call comms_abort
    end if
#endif

    ! Check this block-column is local to this processor
    if (iblk < my_first_blk .or. iblk > my_last_blk) then
       write(stdout,'(a,i6)') 'Error in sparse_get_block_real: &
            &requested block is not local to node',pub_my_node_id
       call comms_abort
    end if

    ! Set block to zero
    blk = 0.0_DP
    nrows = 0

    ! Find information about this segment
    seg = node_of_elem(first_elem_on_atom(jblk,row_blks),row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type == SEG_BLANK) return

    ! Find number of rows in this block/segment
    if (seg_type == SEG_DENSE) then
       nrows = num_elems_on_node(seg,row_blks)
    else if (seg_type == SEG_SPARSE) then
       nrows = num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - my_first_blk + 1
    seg_start = library(ilib)%seg_info(s_idx,seg)
    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       do ielem=1,num_elems_on_atom(iblk,col_blks)
          do jelem=1,num_elems_on_atom(jblk,row_blks)
             blk(jelem,ielem) = mat%dmtx(ptr+(ielem-1)*nrows+jelem-1)
          end do
       end do
    end do

  end subroutine sparse_get_block_real

  !============================================================================!
  ! This subroutine inserts a block of a given block sparse matrix local to    !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   blk (input)    : The block to insert into the real matrix                !
  !   mat (output)   : The real matrix into which to insert the block          !
  !   jblk (input)   : The atom index of the row to insert                     !
  !   iblk (input)   : The atom index of the column to insert                  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in July 2006.                       !
  !============================================================================!

  subroutine sparse_put_block_real(blk,mat,jblk,iblk)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id
    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: blk(:,:)     ! The block to insert
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: jblk               ! The atom index of the row
    integer, intent(in) :: iblk               ! The atom index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: kblk      ! Block counter
    integer :: ielem     ! Element row counter
    integer :: jelem     ! Element column counter
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    integer :: nrows     ! Total number of rows in this block/segment

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks
    nrows = 0

    ! Check arguments
    if (mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_block_real: &
            &real matrices only'
       call comms_abort
    end if
#ifdef DEBUG
    if (jblk < 1 .or. jblk > library(ilib)%nblk) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_block_real: &
            &invalid row index'
       call comms_abort
    end if
    if (iblk < 1 .or. iblk > library(ilib)%nblk) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_block_real: &
            &invalid column index'
       call comms_abort
    end if
    if (size(blk,1) < num_elems_on_atom(jblk,row_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_block_real: &
            &too few rows in blk'
       call comms_abort
    end if
    if (size(blk,2) < num_elems_on_atom(iblk,col_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_block_real: &
            &too few columns in blk'
       call comms_abort
    end if
#endif

    ! Check this block-column is local to this processor
    if (iblk < my_first_blk .or. iblk > my_last_blk) then
       write(stdout,'(a,i6)') 'Error in sparse_put_block_real: &
            &requested block is not local to node',pub_my_node_id
       call comms_abort
    end if

    ! Find information about this segment
    seg = node_of_elem(first_elem_on_atom(jblk,row_blks),row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type == SEG_BLANK) return

    ! Find number of rows in this block/segment
    if (seg_type == SEG_DENSE) then
       nrows = num_elems_on_node(seg,row_blks)
    else if (seg_type == SEG_SPARSE) then
       nrows = num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - my_first_blk + 1
    seg_start = library(ilib)%seg_info(s_idx,seg)
    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       do ielem=1,num_elems_on_atom(iblk,col_blks)
          do jelem=1,num_elems_on_atom(jblk,row_blks)
             mat%dmtx(ptr+(ielem-1)*nrows+jelem-1) = blk(jelem,ielem)
          end do
       end do
    end do

  end subroutine sparse_put_block_real

  !============================================================================!
  ! This subroutine returns a column of a given block sparse matrix local to   !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   data (output)  : The required column of the real matrix                  !
  !   mat (input)    : The real matrix whose column is required                !
  !   elem (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in May 2004.                        !
  !============================================================================!

  subroutine sparse_get_col_real(data,mat,elem)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id, &
         pub_total_num_nodes

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: data(:)    ! The required column
    type(SPAM3), intent(in) :: mat             ! The matrix
    integer, intent(in) :: elem                ! The elem index of the col

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-column of the requested column
    integer :: ielem0    ! The first elem in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: jblk      ! Block-row counter
    integer :: jelems    ! Number of rows in block-row
    integer :: jelem0    ! First row in block-row
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! Check arguments
    if (mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_col_real: &
            &real matrices only'
       call comms_abort
    end if
    if (size(data) < library(ilib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_col_real: &
            &insufficient space in output array'
       call comms_abort
    end if
    if (elem < 1 .or. elem > library(ilib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_col_real: &
            &invalid column index'
       call comms_abort
    end if
    ! Check this column is local to this processor
    if (elem<first_elem_on_node(pub_my_node_id,col_blks).or.&
       elem>first_elem_on_node(pub_my_node_id+1,col_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_col_real: &
          &requested column is not local to this node'
       call comms_abort
    end if

    ! Loop over node-node segments
    do seg=0,pub_total_num_nodes-1

       ! Segment is dense
       if (library(ilib)%seg_info(s_type,seg) == SEG_DENSE) then

          jelems = num_elems_on_node(seg,row_blks)
          ptr = library(ilib)%seg_info(s_ptr,seg) + &
               (elem - first_elem_on_node(pub_my_node_id,col_blks)) * &
               num_elems_on_node(seg,row_blks)

          data(first_elem_on_node(seg,row_blks):first_elem_on_node(seg+1, &
               row_blks)-1) = mat%dmtx(ptr:ptr+jelems-1)

       ! Segment is sparse
       else if (library(ilib)%seg_info(s_type,seg) == SEG_SPARSE) then

          ! Obtain information about the relevant block-column for this elem
          iblk = atom_of_elem(elem,col_blks)
          ielem0 = first_elem_on_atom(iblk,col_blks)

          ! Find local index for the block-column
          loc_iblk = iblk - my_first_blk + 1
          seg_start = library(ilib)%seg_info(s_idx,seg)
          do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
               library(ilib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(ilib)%blk_idx(iovlap)   ! block-row
             jelems = num_elems_on_atom(jblk,row_blks)
             jelem0 = first_elem_on_atom(jblk,row_blks)
             ptr = library(ilib)%blk_ptr(iovlap) + (elem - ielem0) * jelems
             data(jelem0:jelem0+jelems-1) = mat%dmtx(ptr:ptr+jelems-1)
          end do

       end if

    end do  ! seg

  end subroutine sparse_get_col_real

  !============================================================================!
  ! This subroutine inserts a column of a given block sparse matrix local to   !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   data (input)   : The column of the real matrix to insert                 !
  !   mat  (output)  : The real matrix whose column is required                !
  !   elem (input)   : The element index of the column to insert               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in May 2004.                        !
  !============================================================================!

  subroutine sparse_put_col_real(data,mat,elem)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id, &
         pub_total_num_nodes

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: data(:)      ! The required column
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: elem               ! The elem index of the col

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-column of the requested column
    integer :: ielem0    ! The first elem in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: jblk      ! Block-row counter
    integer :: jelems    ! Number of rows in block-row
    integer :: jelem0    ! First row in block-row
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index

    ! get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! Check arguments
    if (mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_col_real: &
            &real matrices only'
       call comms_abort
    end if
    if (size(data) < library(ilib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_col_real: &
            &insufficient space in output array'
       call comms_abort
    end if
    if (elem < 1 .or. elem > library(ilib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_col_real: &
            &invalid column index'
       call comms_abort
    end if
    ! Check this column is local to this processor
    if (elem<first_elem_on_node(pub_my_node_id,col_blks).or.&
       elem>first_elem_on_node(pub_my_node_id+1,col_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_col_real: &
          &requested column is not local to this node'
       call comms_abort
    end if

    ! Loop over node-node segments
    do seg=0,pub_total_num_nodes-1

       ! Segment is dense
       if (library(ilib)%seg_info(s_type,seg) == SEG_DENSE) then

          jelems = num_elems_on_node(seg,row_blks)
          ptr = library(ilib)%seg_info(s_ptr,seg) + &
               (elem - first_elem_on_node(pub_my_node_id,col_blks)) * &
               num_elems_on_node(seg,row_blks)

          mat%dmtx(ptr:ptr+jelems-1) = data(first_elem_on_node(seg,row_blks): &
               first_elem_on_node(seg+1,row_blks)-1)

       ! Segment is sparse
       else if (library(ilib)%seg_info(s_type,seg) == SEG_SPARSE) then

          ! Obtain information about the relevant block-column for this elem
          iblk = atom_of_elem(elem,col_blks)
          ielem0 = first_elem_on_atom(iblk,col_blks)

          ! Find local index for the block-column
          loc_iblk = iblk - my_first_blk + 1
          seg_start = library(ilib)%seg_info(s_idx,seg)
          do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
               library(ilib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(ilib)%blk_idx(iovlap)   ! block-row
             jelems = num_elems_on_atom(jblk,row_blks)
             jelem0 = first_elem_on_atom(jblk,row_blks)
             ptr = library(ilib)%blk_ptr(iovlap) + (elem - ielem0) * jelems
             mat%dmtx(ptr:ptr+jelems-1) = data(jelem0:jelem0+jelems-1)
          end do

       ! Nothing to do if segment is blank

       end if

    end do  ! seg

  end subroutine sparse_put_col_real

  !============================================================================!
  ! This subroutine clears an array according to the sparsity pattern of a     !
  ! column of a given block sparse matrix local to this processor.             !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   data (output)  : The real array to clear                                 !
  !   mat  (output)  : The real matrix whose sparsity pattern is to be used    !
  !   elem (input)   : The element index of the column to use                  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in May 2004.                        !
  !============================================================================!

  subroutine sparse_clr_col_real(data,mat,elem)

    use comms, only: comms_abort, pub_on_root, pub_total_num_nodes

    implicit none

    ! Arguments
    real(kind=DP), intent(inout) :: data(:)    ! The column to insert
    type(SPAM3), intent(in) :: mat             ! The matrix
    integer, intent(in) :: elem                ! The elem index of the col

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-column of the requested column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: jblk      ! Block-row counter
    integer :: jelems    ! Number of rows in block-row
    integer :: jelem0    ! First row in block-row
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! Check arguments
    if (mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_clr_col_real: &
            &real matrices only'
       call comms_abort
    end if
    if (size(data) < library(ilib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_clr_col_real: &
            &insufficient space in input array'
       call comms_abort
    end if
    if (elem < 1 .or. elem > library(ilib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_clr_col_real: &
            &invalid column index'
       call comms_abort
    end if

    ! Obtain information about the relevant block-column for this elem
    iblk = atom_of_elem(elem,col_blks)

    ! Check this block-column is local to this processor
    if (iblk < my_first_blk .or. iblk > my_last_blk) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_clr_col_real: &
            &requested column is not local to this node'
       call comms_abort
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - my_first_blk + 1

    do seg=0,pub_total_num_nodes-1

       ! Nothing to do if segment is blank
       if (library(ilib)%seg_info(s_type,seg)==SEG_BLANK) cycle

       ! Loop over the nonzero blocks of this segment
       seg_start = library(ilib)%seg_info(s_idx,seg)
       do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
            library(ilib)%blk_idx(seg_start+loc_iblk)-1
          jblk = library(ilib)%blk_idx(iovlap)   ! block-row

          ! Set the nonzero elements of this col to zero
          jelems = num_elems_on_atom(jblk,row_blks)
          jelem0 = first_elem_on_atom(jblk,row_blks)
          data(jelem0:jelem0+jelems-1) = 0.0_DP

       end do  ! iovlap

    end do  ! seg

  end subroutine sparse_clr_col_real


  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3 into a dense     !
  ! full square real array, or collates a dense matrix across nodes to one     !
  ! single full square array                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination full square matrix                     !
  !   src    (input)  : The source block sparse matrix                         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  subroutine sparse_spam3tofull_real(dest,src)

    use comms, only: comms_abort, comms_bcast, pub_on_root, pub_total_num_nodes

    implicit none

    ! Arguments
    real(kind=DP), intent(out) :: dest(:,:)     ! The full square matrix
    type(SPAM3), intent(in) :: src              ! The block sparse matrix

    ! Local variables
    integer :: srclib        ! Library entry for src
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! elem in atom block iblk counter
    integer :: jelemonat     ! elem in atom block jblk counter
    integer :: ielem,jelem   ! elem counters
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in this block/segment
    integer :: node          ! Node loop counter
    integer :: seg           ! Segment index
    integer :: seg_type      ! Start position of this segment in the index
    integer :: seg_start     ! Start position of this segment in the index

    ! Get library entry for src
    srclib = src%lib
    row_blks = library(srclib)%row_blks
    col_blks = library(srclib)%col_blks
    nrows = 0

    ! Check arguments
    if (size(dest,1) /= library(srclib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_spam3tofull_real: &
            &incompatible numbers of rows in arguments.'
       call comms_abort
    end if
    if (size(dest,2) /= library(srclib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_spam3tofull_real: &
            &incompatible numbers of cols in arguments.'
       call comms_abort
    end if
    if (src%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_spam3tofull_real: &
            &real matrices only.'
       call comms_abort
    end if

    ! Zero full square matrix
    dest(:,:) = 0.0_DP

    ! Loop over the segments of src on this node
    do seg=0,pub_total_num_nodes-1

       ! Cycle if this segment is blank
       seg_type = library(srclib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = num_elems_on_node(seg,row_blks)

       ! Loop over block-columns of src on this node
       seg_start = library(srclib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=my_first_blk,my_last_blk
          loc_iblk = loc_iblk + 1
          ielem = first_elem_on_atom(iblk,col_blks)
          ielems = num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of src
          do idx=library(srclib)%blk_idx(seg_start+loc_iblk-1), &
               library(srclib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(srclib)%blk_idx(idx)
             jelem = first_elem_on_atom(jblk,row_blks)
             jelems = num_elems_on_atom(jblk,row_blks)

             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ptr = library(srclib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                do jelemonat=0,jelems-1   ! Element row in block
                   dest(jelem+jelemonat,ielem+ielemonat) = &
                        src%dmtx(ptr+ielemonat*nrows+jelemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of src

       end do  ! Loop over block-columns of src

    end do  ! seg

    ! Broadcast result across all nodes
    do node=0,pub_total_num_nodes-1

       ! Get first column and number of columns
       ielem = first_elem_on_node(node,col_blks)
       ielems = num_elems_on_node(node,col_blks)

       ! Broadcast this part of the matrix
       call comms_bcast(node,dest(1,ielem),ielems*library(srclib)%nrows)

    end do  ! Loop over nodes

  end subroutine sparse_spam3tofull_real

  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3 into a dense     !
  ! full square real array, or collates a dense matrix across nodes to one     !
  ! single full square array                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination full square matrix                     !
  !   src    (input)  : The source block sparse matrix                         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  subroutine sparse_fulltospam3_real(dest,src)

    use comms, only: comms_abort, pub_on_root, pub_total_num_nodes

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: dest       ! The block sparse matrix
    real(kind=DP), intent(in) :: src(:,:)    ! The full square matrix

    ! Local variables
    integer :: destlib       ! Library entry for dest
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! Element in atom block iblk counter
    integer :: jelemonat     ! Element in atom block jblk counter
    integer :: ielem,jelem   ! Element start positions of a block
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in the block/segment
    integer :: seg           ! Segment index
    integer :: seg_type      ! Segment type for this segment
    integer :: seg_start     ! Start position of this segment in the index

    ! Get library entry for src
    destlib = dest%lib
    row_blks = library(destlib)%row_blks
    col_blks = library(destlib)%col_blks
    nrows = 0

    ! Check arguments
    if (size(src,1) /= library(destlib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_fulltospam3_real: &
            &incompatible numbers of rows in arguments.'
       call comms_abort
    end if
    if (size(src,2) /= library(destlib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_fulltospam3_real: &
            &incompatible numbers of cols in arguments.'
       call comms_abort
    end if
    if (dest%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_fulltospam3_real: &
            &real matrices only.'
       call comms_abort
    end if

    ! Loop over the segments of dest on this node
    do seg=0,pub_total_num_nodes-1

       ! Cycle if this segment is blank
       seg_type = library(destlib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = num_elems_on_node(seg,row_blks)

       ! Loop over block-columns of dest on this node
       seg_start = library(destlib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=my_first_blk,my_last_blk
          loc_iblk = loc_iblk + 1
          ielems = num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of dest
          do idx=library(destlib)%blk_idx(seg_start+loc_iblk-1), &
               library(destlib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(destlib)%blk_idx(idx)
             jelems = num_elems_on_atom(jblk,row_blks)
             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ielem = first_elem_on_atom(iblk,col_blks)
             ptr = library(destlib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                jelem = first_elem_on_atom(jblk,row_blks)
                do jelemonat=0,jelems-1   ! Element row in block
                   dest%dmtx(ptr+ielemonat*nrows+jelemonat) = &
                        src(jelem+jelemonat,ielem+ielemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of dest

       end do  ! Loop over block-columns of dest

    end do  ! seg

  end subroutine sparse_fulltospam3_real

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues of a matrix.  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (inout) : The matrix to be rescaled and shifted                   !
  !   alpha  (input) : The real scaling parameter                              !
  !   beta   (input) : The optional real shift parameter                       !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  !============================================================================!

  subroutine sparse_scale_real(mat,alpha,beta)

    use comms, only: comms_abort, pub_my_node_id

    implicit none

    ! Arguments
    type(SPAM3), intent(inout)          :: mat    ! The matrix to be operated on
    real(kind=DP), intent(in)           :: alpha  ! The scaling parameter
    real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: blks      ! Identifier for blocking scheme
    integer :: iblk      ! Block-column loop counter
    integer :: loc_iblk  ! Loop counter for iblk on this node
    integer :: ielems    ! Number of columns in block
    integer :: ptr       ! Pointer to diagonal entries
    integer :: ielem     ! Element in block-column loop counter
    integer :: seg       ! Segment index
    integer :: seg_start ! Segment type for this segment


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

       ! Get library entry for mat
       ilib = mat%lib
       blks = library(ilib)%row_blks
       seg = pub_my_node_id

       ! Can only shift eigenvalues of square matrix
       if (library(ilib)%row_blks/=library(ilib)%col_blks) then
          write(stdout,'(a)') 'Error in sparse_scale: cannot shift &
               &eigenvalues of non-square matrices'
          call comms_abort
       end if

       if (library(ilib)%seg_info(s_type,seg)==SEG_DENSE) then

          ! Set pointer to first diagonal element
          ptr = library(ilib)%seg_info(s_ptr,seg)

          ! Get number of columns in this segment
          ielems = num_elems_on_node(seg,blks)

          ! Add identity matrix scaled by beta to all elements
          if (mat%iscmplx) then
             ! Loop over columns on my node
             do ielem=1,ielems
                mat%zmtx(ptr) = mat%zmtx(ptr) + cmplx(beta,0.0_DP,kind=DP)
                ptr = ptr + ielems + 1
             end do
          else
             ! Loop over columns on my node
             do ielem=1,num_elems_on_node(seg,blks)
                mat%dmtx(ptr) = mat%dmtx(ptr) + beta
                ptr = ptr + ielems + 1
             end do
          end if

       else if (library(ilib)%seg_info(s_type,seg)==SEG_SPARSE) then

          seg_start = library(ilib)%seg_info(s_idx,seg)

          ! Loop over atom-blocks iblk on my node
          loc_iblk = 0
          do iblk=my_first_blk,my_last_blk
             loc_iblk = loc_iblk + 1

             ! Get number of elements on this atom
             ielems = num_elems_on_atom(iblk,blks)

             ! Find pointer to start of diagonal block in this column
             ptr = library(ilib)%blk_ptr(seg_start+loc_iblk-1)

             ! Loop over diagonal elements in this diagonal block and add
             ! identity matrix scaled by beta
             if (mat%iscmplx) then
                do ielem=1,ielems
                   mat%zmtx(ptr) = mat%zmtx(ptr) + cmplx(beta,0.0_DP,kind=DP)
                   ptr = ptr + ielems + 1
                end do
             else
                do ielem=1,ielems
                   mat%dmtx(ptr) = mat%dmtx(ptr) + beta
                   ptr = ptr + ielems + 1
                end do
             end if

          end do  ! Loop over atom-blocks iblk on this node

       endif

    end if


  end subroutine sparse_scale_real


  !============================================================================!
  ! This subroutine performs an axpy operation on the matrices:                !
  !   y := y + alpha * x                                                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The sparse matrix y                                      !
  !   xmat  (input) : The sparse matrix x                                      !
  !   alpha (input) : The real parameter alpha                                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  !============================================================================!

  subroutine sparse_axpy_real(ymat,xmat,alpha)

    use comms, only: comms_abort, pub_total_num_nodes

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: ymat    ! The sparse matrix y
    type(SPAM3), intent(in) :: xmat       ! The sparse matrix x
    real(kind=DP), intent(in) :: alpha    ! The parameter alpha

    ! Local variables
    integer :: xlib,ylib    ! Library pointers for x and y
    integer :: row_blks     ! Identifier for column blocking scheme
    integer :: col_blks     ! Identifier for row blocking scheme
    integer :: iblk,jblk    ! Atom loop counters
    integer :: icol         ! Atom loop counters
    integer :: loc_iblk     ! Atom loop counter on this node
    integer :: idx          ! Index loop counter
    integer :: iel          ! Element loop counter
    integer :: ifound       ! Loop counter over columns found
    integer :: nfound       ! Number of columns found
    integer :: xstart,xend  ! Start and end of block data of xmat
    integer :: ystart,yend  ! Start and end of block data of ymat
    integer :: ystride      ! Stride to get to next column of src
    integer :: xstride      ! Stride to get to next column of dest
    integer :: seg          ! Segment index
    integer :: x_seg_type   ! Segment type for segments of x
    integer :: y_seg_type   ! Segment type for segments of y
    integer :: x_seg_start  ! Start position of segments of x in the index
    integer :: y_seg_start  ! Start position of segments of y in the index


    ! Extract library entries
    xlib = xmat%lib
    ylib = ymat%lib
    row_blks = library(xlib)%row_blks
    col_blks = library(xlib)%col_blks
    xstride = 0; ystride = 0

    ! Check matrices share same blocking scheme
    if ( (library(xlib)%row_blks /= row_blks) .or. &
         (library(xlib)%col_blks /= col_blks) ) then
       write(stdout,'(a)') 'Error in sparse_axpy_real: x and y matrices &
            &have different blocking schemes'
       call comms_abort
    end if

    ! Check for alpha = 0 (no operation required)
    if (alpha == 0.0_DP) then
       return
    end if

    ! Check whether the matrices share the same structure
    if ( xlib == ylib ) then

       ! Trivial axpy of whole data arrays
       if (xmat%iscmplx) then
          if (ymat%iscmplx) then
             ymat%zmtx = ymat%zmtx + alpha * xmat%zmtx
          else
             do iel=1,size(ymat%dmtx)
                ymat%dmtx(iel) = ymat%dmtx(iel) + &
                     alpha * real(xmat%zmtx(iel),kind=DP)
             end do
          end if
       else
          if (ymat%iscmplx) then
             do iel=1,size(ymat%zmtx)
                ymat%zmtx(iel) = ymat%zmtx(iel) + &
                     alpha * cmplx(xmat%dmtx(iel),0.0_DP,kind=DP)
             end do
          else
             ymat%dmtx = ymat%dmtx + alpha * xmat%dmtx
          end if
       end if

    ! Structures differ, so proceed segment-by-segment
    else

       do seg=0,pub_total_num_nodes-1

          x_seg_type = library(xlib)%seg_info(s_type,seg)
          y_seg_type = library(ylib)%seg_info(s_type,seg)

          ! Nothing to add to if this segment in y is blank
          if (y_seg_type==SEG_BLANK) cycle

          ! Nothing to add if this segment in x is blank
          if (x_seg_type == SEG_BLANK) cycle

          ! Find information about these segments
          x_seg_start = library(xlib)%seg_info(s_idx,seg)
          y_seg_start = library(ylib)%seg_info(s_idx,seg)
          if (y_seg_type == SEG_DENSE) ystride = &
               num_elems_on_node(seg,row_blks)
          if (x_seg_type == SEG_DENSE) xstride = &
               num_elems_on_node(seg,row_blks)

          ! Loop over atom blocks (block-columns) on this node
          loc_iblk = 0
          do iblk=my_first_blk,my_last_blk
             loc_iblk = loc_iblk + 1

             ! Reset counters of number of rows found
             nfound = 0

             ! Loop over block-rows in x for this segment
             do idx=library(xlib)%blk_idx(x_seg_start+loc_iblk-1), &
                  library(xlib)%blk_idx(x_seg_start+loc_iblk)-1
                jblk = library(xlib)%blk_idx(idx)

                ! Mark flag for this row and add to list
                nfound = nfound + 1
                found_blk(jblk) = .true.
                found_idx(nfound) = jblk
                found_ptr(jblk) = library(xlib)%blk_ptr(idx)

             end do

             ! Loop over block-rows in y
             do idx=library(ylib)%blk_idx(y_seg_start+loc_iblk-1), &
                  library(ylib)%blk_idx(y_seg_start+loc_iblk)-1
                jblk = library(ylib)%blk_idx(idx)

                ! Check whether this block-row exists in source
                if (found_blk(jblk)) then

                   ! Find pointers to start of dest and src data
                   ystart = library(ylib)%blk_ptr(idx)
                   yend = ystart + num_elems_on_atom(jblk,row_blks) - 1
                   xstart = found_ptr(jblk)
                   xend = xstart + num_elems_on_atom(jblk,row_blks) - 1

                   ! Find stride to next column of dest data
                   if (y_seg_type == SEG_SPARSE) then
                      ystride = num_elems_on_atom(jblk,row_blks)
                   end if

                   ! Find stride to next column of source data
                   if (x_seg_type == SEG_SPARSE) then
                      xstride = num_elems_on_atom(jblk,row_blks)
                   end if

                   ! Copy columns from source to destination
                   do icol=1,num_elems_on_atom(iblk,col_blks)

                      if (xmat%iscmplx) then
                         if (ymat%iscmplx) then
                            ymat%zmtx(ystart:yend) = &
                                 ymat%zmtx(ystart:yend) + &
                                 alpha * xmat%zmtx(xstart:xend)
                         else
                            do iel=ystart,yend
                               ymat%dmtx(iel) = ymat%dmtx(iel) + alpha * &
                                    real(xmat%zmtx(iel-ystart+xstart), &
                                    kind=DP)
                            end do
                         end if
                      else
                         if (ymat%iscmplx) then
                            do iel=ystart,yend
                               ymat%zmtx(iel) = ymat%zmtx(iel) + alpha * &
                                    cmplx(xmat%dmtx(iel-ystart+xstart), &
                                    0.0_DP,kind=DP)
                            end do
                         else
                            ymat%dmtx(ystart:yend) = &
                                 ymat%dmtx(ystart:yend) + &
                                 alpha * xmat%dmtx(xstart:xend)
                         end if
                      end if

                      ystart = ystart + ystride
                      yend = yend + ystride
                      xstart = xstart + xstride
                      xend = xend + xstride

                   end do

                end if

             end do  ! Loop over destination block-rows

             ! Reset flags
             do ifound=1,nfound
                jblk = found_idx(ifound)
                found_blk(jblk) = .false.
             end do

          end do  ! Loop over atom block-columns

       end do  ! seg

    end if


  end subroutine sparse_axpy_real


  !============================================================================!
  ! This subroutine inserts an element into a given block sparse matrix local  !
  ! to this processor.                                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (input)     : The element of the real matrix to insert                !
  !   mat (inout)    : The real matrix                                         !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in June 2006.                       !
  !============================================================================!

  subroutine sparse_put_element_complex(el,mat,jrow,icol)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id
    implicit none

    ! Arguments
    complex(kind=DP), intent(in) :: el        ! The element to insert
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-row to which the element belongs
    integer :: jblk      ! The block-col to which the element belongs
    integer :: kblk      ! Block counter
    integer :: nrows     ! Total number of rows in this block/segment
    integer :: jelem0    ! The first element in the requested block-row
    integer :: ielem0    ! The first element in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks
    nrows = 0

    ! Check arguments
    if (.not.mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_put_element_complex: complex matrices only'
       call comms_abort
    end if
#ifdef DEBUG
    if (jrow < 1 .or. jrow > library(ilib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_put_element_complex: invalid row index'
       call comms_abort
    end if
    if (icol < 1 .or. icol > library(ilib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_put_element_complex: invalid column index'
       call comms_abort
    end if
#endif
    ! Check this column is local to this processor
    if (icol<first_elem_on_node(pub_my_node_id,col_blks).or.&
         icol>=first_elem_on_node(pub_my_node_id+1,col_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_put_element_complex: requested column is not local to &
            &this node'
       call comms_abort
    end if

    ! Obtain information about the relevant block & segment for this element
    seg = node_of_elem(jrow,row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type==SEG_BLANK) return
    seg_start = library(ilib)%seg_info(s_idx,seg)
    jblk = atom_of_elem(jrow,row_blks)
    iblk = atom_of_elem(icol,col_blks)
    jelem0 = first_elem_on_atom(jblk,row_blks)
    ielem0 = first_elem_on_atom(iblk,col_blks)
    if (seg_type==SEG_DENSE) then
       nrows = num_elems_on_node(seg,row_blks)
    else if (seg_type==SEG_SPARSE) then
       nrows = num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - my_first_blk + 1

    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       mat%zmtx(ptr + (icol - ielem0)*nrows + (jrow - jelem0)) = el
    end do

  end subroutine sparse_put_element_complex

  !============================================================================!
  ! This subroutine retrieves an element into a given block sparse matrix      !
  ! local to this processor.                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   el (output)    : The element of the real matrix retrieved                !
  !   mat (inout)    : The real matrix                                         !
  !   jrow (input)   : The element index of the row required                   !
  !   icol (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in June 2006.                       !
  !============================================================================!

  subroutine sparse_get_element_complex(el,mat,jrow,icol)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id

    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: el       ! The element to insert
    type(SPAM3), intent(in) :: mat            ! The matrix
    integer, intent(in) :: jrow               ! The index of the row
    integer, intent(in) :: icol               ! The index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-row to which the element belongs
    integer :: jblk      ! The block-col to which the element belongs
    integer :: kblk      ! Block counter
    integer :: nrows     ! Total number of rows in this block/segment
    integer :: jelem0    ! The first element in the requested block-row
    integer :: ielem0    ! The first element in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! Check arguments
    if (.not.mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_element_complex: &
            &complex matrices only'
       call comms_abort
    end if
#ifdef DEBUG
    if (jrow < 1 .or. jrow > library(ilib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_element_complex: &
            &invalid row index'
       call comms_abort
    end if
    if (icol < 1 .or. icol > library(ilib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_element_complex: &
            &invalid column index'
       call comms_abort
    end if
#endif
    ! Check this column is local to this processor
    if (icol<first_elem_on_node(pub_my_node_id,col_blks).or.&
         icol>=first_elem_on_node(pub_my_node_id+1,col_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_element_complex: &
          &requested column is not local to this node'
       call comms_abort
    end if

    nrows = 0
    el = 0.0_DP

    ! Obtain information about the relevant block & segment for this element
    seg = node_of_elem(jrow,row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type==SEG_BLANK) return
    seg_start = library(ilib)%seg_info(s_idx,seg)
    jblk = atom_of_elem(jrow,row_blks)
    iblk = atom_of_elem(icol,col_blks)
    jelem0 = first_elem_on_atom(jblk,row_blks)
    ielem0 = first_elem_on_atom(iblk,col_blks)
    if (seg_type==SEG_DENSE) then
       nrows = num_elems_on_node(seg,row_blks)
    else if (seg_type==SEG_SPARSE) then
       nrows = num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - my_first_blk + 1

    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       el = mat%zmtx(ptr + (icol - ielem0)*nrows + (jrow - jelem0))
    end do

  end subroutine sparse_get_element_complex

  !============================================================================!
  ! This subroutine returns a block of a given block sparse matrix local to    !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   blk (output)   : The required block of the real matrix                   !
  !   mat (input)    : The real matrix whose element is required               !
  !   jblk (input)   : The atom index of the row required                      !
  !   iblk (input)   : The atom index of the column required                   !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in July 2006.                       !
  !============================================================================!

  subroutine sparse_get_block_complex(blk,mat,jblk,iblk)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id
    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: blk(:,:) ! The required block
    type(SPAM3), intent(in) :: mat            ! The matrix
    integer, intent(in) :: jblk               ! The atom index of the row
    integer, intent(in) :: iblk               ! The atom index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: kblk      ! Block counter
    integer :: ielem     ! Element row counter
    integer :: jelem     ! Element column counter
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    integer :: nrows     ! Total number of rows in this block/segment

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! Check arguments
    if (.not.mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_block_complex: &
            &complex matrices only'
       call comms_abort
    end if
#ifdef DEBUG
    if (jblk < 1 .or. jblk > library(ilib)%nblk) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_block_complex: &
            &invalid row index'
       call comms_abort
    end if
    if (iblk < 1 .or. iblk > library(ilib)%nblk) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_block_complex: &
            &invalid column index'
       call comms_abort
    end if
    if (size(blk,1) < num_elems_on_atom(jblk,row_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_block_complex: &
            &too few rows in blk'
       call comms_abort
    end if
    if (size(blk,2) < num_elems_on_atom(iblk,col_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_block_complex: &
            &too few columns in blk'
       call comms_abort
    end if
#endif

    ! Check this block-column is local to this processor
    if (iblk < my_first_blk .or. iblk > my_last_blk) then
       write(stdout,'(a,i6)') 'Error in sparse_get_block_complex: &
            &requested block is not local to node',pub_my_node_id
       call comms_abort
    end if

    ! Set block to zero
    blk = 0.0_DP
    nrows = 0

    ! Find information about this segment
    seg = node_of_elem(first_elem_on_atom(jblk,row_blks),row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type == SEG_BLANK) return

    ! Find number of rows in this block/segment
    if (seg_type == SEG_DENSE) then
       nrows = num_elems_on_node(seg,row_blks)
    else if (seg_type == SEG_SPARSE) then
       nrows = num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - my_first_blk + 1
    seg_start = library(ilib)%seg_info(s_idx,seg)
    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       do ielem=1,num_elems_on_atom(iblk,col_blks)
          do jelem=1,num_elems_on_atom(jblk,row_blks)
             blk(jelem,ielem) = mat%zmtx(ptr+(ielem-1)*nrows+jelem-1)
          end do
       end do
    end do

  end subroutine sparse_get_block_complex

  !============================================================================!
  ! This subroutine inserts a block of a given block sparse matrix local to    !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   blk (input)    : The block to insert into the real matrix                !
  !   mat (output)   : The real matrix into which to insert the block          !
  !   jblk (input)   : The atom index of the row to insert                     !
  !   iblk (input)   : The atom index of the column to insert                  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in July 2006.                       !
  !============================================================================!

  subroutine sparse_put_block_complex(blk,mat,jblk,iblk)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id
    implicit none

    ! Arguments
    complex(kind=DP), intent(in) :: blk(:,:)  ! The block to insert
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: jblk               ! The atom index of the row
    integer, intent(in) :: iblk               ! The atom index of the column

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: kblk      ! Block counter
    integer :: ielem     ! Element row counter
    integer :: jelem     ! Element column counter
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_type  ! Segment type for this segment
    integer :: seg_start ! Start position of this segment in the index
    integer :: nrows     ! Total number of rows in this block/segment

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks
    nrows = 0

    ! Check arguments
    if (.not.mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_block_complex: &
            &complex matrices only'
       call comms_abort
    end if
#ifdef DEBUG
    if (jblk < 1 .or. jblk > library(ilib)%nblk) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_block_complex: &
            &invalid row index'
       call comms_abort
    end if
    if (iblk < 1 .or. iblk > library(ilib)%nblk) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_block_complex: &
            &invalid column index'
       call comms_abort
    end if
    if (size(blk,1) < num_elems_on_atom(jblk,row_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_block_complex: &
            &too few rows in blk'
       call comms_abort
    end if
    if (size(blk,2) < num_elems_on_atom(iblk,col_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_block_complex: &
            &too few columns in blk'
       call comms_abort
    end if
#endif

    ! Check this block-column is local to this processor
    if (iblk < my_first_blk .or. iblk > my_last_blk) then
       write(stdout,'(a,i6)') 'Error in sparse_put_block_complex: &
            &requested block is not local to node',pub_my_node_id
       call comms_abort
    end if

    ! Find information about this segment
    seg = node_of_elem(first_elem_on_atom(jblk,row_blks),row_blks)
    seg_type = library(ilib)%seg_info(s_type,seg)
    if (seg_type == SEG_BLANK) return

    ! Find number of rows in this block/segment
    if (seg_type == SEG_DENSE) then
       nrows = num_elems_on_node(seg,row_blks)
    else if (seg_type == SEG_SPARSE) then
       nrows = num_elems_on_atom(jblk,row_blks)
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - my_first_blk + 1
    seg_start = library(ilib)%seg_info(s_idx,seg)
    do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
         library(ilib)%blk_idx(seg_start+loc_iblk)-1
       kblk = library(ilib)%blk_idx(iovlap)   ! block-row
       if (kblk < jblk) cycle
       if (kblk > jblk) exit
       ptr = library(ilib)%blk_ptr(iovlap)
       do ielem=1,num_elems_on_atom(iblk,col_blks)
          do jelem=1,num_elems_on_atom(jblk,row_blks)
             mat%zmtx(ptr+(ielem-1)*nrows+jelem-1) = blk(jelem,ielem)
          end do
       end do
    end do

  end subroutine sparse_put_block_complex

  !============================================================================!
  ! This subroutine returns a column of a given block sparse matrix local to   !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   data (output)  : The required column of the real matrix                  !
  !   mat (input)    : The real matrix whose column is required                !
  !   elem (input)   : The element index of the column required                !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in May 2004.                        !
  !============================================================================!

  subroutine sparse_get_col_complex(data,mat,elem)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id, &
         pub_total_num_nodes

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout) :: data(:) ! The required column
    type(SPAM3), intent(in) :: mat             ! The matrix
    integer, intent(in) :: elem                ! The elem index of the col

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-column of the requested column
    integer :: ielem0    ! The first elem in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: jblk      ! Block-row counter
    integer :: jelems    ! Number of rows in block-row
    integer :: jelem0    ! First row in block-row
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! Check arguments
    if (.not.mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_col_complex: &
            &complex matrices only'
       call comms_abort
    end if
    if (size(data) < library(ilib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_col_complex: &
            &insufficient space in output array'
       call comms_abort
    end if
    if (elem < 1 .or. elem > library(ilib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_col_complex: &
            &invalid column index'
       call comms_abort
    end if
    ! Check this column is local to this processor
    if (elem<first_elem_on_node(pub_my_node_id,col_blks).or.&
       elem>first_elem_on_node(pub_my_node_id+1,col_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_get_col_complex: &
          &requested column is not local to this node'
       call comms_abort
    end if

    ! Loop over node-node segments
    do seg=0,pub_total_num_nodes-1

       ! Segment is dense
       if (library(ilib)%seg_info(s_type,seg) == SEG_DENSE) then

          jelems = num_elems_on_node(seg,row_blks)
          ptr = library(ilib)%seg_info(s_ptr,seg) + &
               (elem - first_elem_on_node(pub_my_node_id,col_blks)) * &
               num_elems_on_node(seg,row_blks)

          data(first_elem_on_node(seg,row_blks):first_elem_on_node(seg+1, &
               row_blks)-1) = mat%zmtx(ptr:ptr+jelems-1)

       ! Segment is sparse
       else if (library(ilib)%seg_info(s_type,seg) == SEG_SPARSE) then

          ! Obtain information about the relevant block-column for this elem
          iblk = atom_of_elem(elem,col_blks)
          ielem0 = first_elem_on_atom(iblk,col_blks)

          ! Find local index for the block-column
          loc_iblk = iblk - my_first_blk + 1
          seg_start = library(ilib)%seg_info(s_idx,seg)
          do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
               library(ilib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(ilib)%blk_idx(iovlap)   ! block-row
             jelems = num_elems_on_atom(jblk,row_blks)
             jelem0 = first_elem_on_atom(jblk,row_blks)
             ptr = library(ilib)%blk_ptr(iovlap) + (elem - ielem0) * jelems
             data(jelem0:jelem0+jelems-1) = mat%zmtx(ptr:ptr+jelems-1)
          end do

       end if

    end do  ! seg

  end subroutine sparse_get_col_complex

  !============================================================================!
  ! This subroutine inserts a column of a given block sparse matrix local to   !
  ! this processor.                                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   data (input)   : The column of the real matrix to insert                 !
  !   mat  (output)  : The real matrix whose column is required                !
  !   elem (input)   : The element index of the column to insert               !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in May 2004.                        !
  !============================================================================!

  subroutine sparse_put_col_complex(data,mat,elem)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id, &
         pub_total_num_nodes

    implicit none

    ! Arguments
    complex(kind=DP), intent(in) :: data(:)   ! The required column
    type(SPAM3), intent(inout) :: mat         ! The matrix
    integer, intent(in) :: elem               ! The elem index of the col

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-column of the requested column
    integer :: ielem0    ! The first elem in the requested block-column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: jblk      ! Block-row counter
    integer :: jelems    ! Number of rows in block-row
    integer :: jelem0    ! First row in block-row
    integer :: ptr       ! Pointer to data array
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! Check arguments
    if (.not.mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_col_complex: &
            &complex matrices only'
       call comms_abort
    end if
    if (size(data) < library(ilib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_col_complex: &
            &insufficient space in output array'
       call comms_abort
    end if
    if (elem < 1 .or. elem > library(ilib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_col_complex: &
            &invalid column index'
       call comms_abort
    end if
    ! Check this column is local to this processor
    if (elem<first_elem_on_node(pub_my_node_id,col_blks).or.&
       elem>first_elem_on_node(pub_my_node_id+1,col_blks)) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_put_col_complex: &
          &requested column is not local to this node'
       call comms_abort
    end if

    ! Loop over node-node segments
    do seg=0,pub_total_num_nodes-1

       ! Segment is dense
       if (library(ilib)%seg_info(s_type,seg) == SEG_DENSE) then

          jelems = num_elems_on_node(seg,row_blks)
          ptr = library(ilib)%seg_info(s_ptr,seg) + &
               (elem - first_elem_on_node(pub_my_node_id,col_blks)) * &
               num_elems_on_node(seg,row_blks)

          mat%zmtx(ptr:ptr+jelems-1) = data(first_elem_on_node(seg,row_blks): &
               first_elem_on_node(seg+1,row_blks)-1)

       ! Segment is sparse
       else if (library(ilib)%seg_info(s_type,seg) == SEG_SPARSE) then

          ! Obtain information about the relevant block-column for this elem
          iblk = atom_of_elem(elem,col_blks)
          ielem0 = first_elem_on_atom(iblk,col_blks)

          ! Find local index for the block-column
          loc_iblk = iblk - my_first_blk + 1
          seg_start = library(ilib)%seg_info(s_idx,seg)
          do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
               library(ilib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(ilib)%blk_idx(iovlap)   ! block-row
             jelems = num_elems_on_atom(jblk,row_blks)
             jelem0 = first_elem_on_atom(jblk,row_blks)
             ptr = library(ilib)%blk_ptr(iovlap) + (elem - ielem0) * jelems
             mat%zmtx(ptr:ptr+jelems-1) = data(jelem0:jelem0+jelems-1)
          end do

       ! Nothing to do if segment is blank

       end if

    end do  ! seg

  end subroutine sparse_put_col_complex

  !============================================================================!
  ! This subroutine clears an array according to the sparsity pattern of a     !
  ! column of a given block sparse matrix local to this processor.             !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   data (output)  : The real array to clear                                 !
  !   mat  (output)  : The real matrix whose sparsity pattern is to be used    !
  !   elem (input)   : The element index of the column to use                  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, in May 2004.                        !
  !============================================================================!

  subroutine sparse_clr_col_complex(data,mat,elem)

    use comms, only: comms_abort, pub_on_root, pub_total_num_nodes

    implicit none

    ! Arguments
    complex(kind=DP), intent(inout) :: data(:) ! The column to insert
    type(SPAM3), intent(in) :: mat             ! The matrix
    integer, intent(in) :: elem                ! The elem index of the col

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: row_blks  ! Identifier for column blocking scheme
    integer :: col_blks  ! Identifier for row blocking scheme
    integer :: iblk      ! The block-column of the requested column
    integer :: iovlap    ! Block overlap counter
    integer :: loc_iblk  ! Block-column counter local to this node
    integer :: jblk      ! Block-row counter
    integer :: jelems    ! Number of rows in block-row
    integer :: jelem0    ! First row in block-row
    integer :: seg       ! Segment index
    integer :: seg_start ! Start position of this segment in the index

    ! Get library entry
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! Check arguments
    if (.not.mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_clr_col_complex: &
            &complex matrices only'
       call comms_abort
    end if
    if (size(data) < library(ilib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_clr_col_complex: &
            &insufficient space in input array'
       call comms_abort
    end if
    if (elem < 1 .or. elem > library(ilib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_clr_col_complex: &
            &invalid column index'
       call comms_abort
    end if

    ! Obtain information about the relevant block-column for this elem
    iblk = atom_of_elem(elem,col_blks)

    ! Check this block-column is local to this processor
    if (iblk < my_first_blk .or. iblk > my_last_blk) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_clr_col_complex: &
            &requested column is not local to this node'
       call comms_abort
    end if

    ! Find local index for the block-column
    loc_iblk = iblk - my_first_blk + 1

    do seg=0,pub_total_num_nodes-1

       ! Nothing to do if segment is blank
       if (library(ilib)%seg_info(s_type,seg)==SEG_BLANK) cycle

       ! Loop over the nonzero blocks of this segment
       seg_start = library(ilib)%seg_info(s_idx,seg)
       do iovlap=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
            library(ilib)%blk_idx(seg_start+loc_iblk)-1
          jblk = library(ilib)%blk_idx(iovlap)   ! block-row

          ! Set the nonzero elements of this col to zero
          jelems = num_elems_on_atom(jblk,row_blks)
          jelem0 = first_elem_on_atom(jblk,row_blks)
          data(jelem0:jelem0+jelems-1) = cmplx(0.0_DP,0.0_DP,kind=DP)

       end do  ! iovlap

    end do  ! seg

  end subroutine sparse_clr_col_complex


  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3 into a dense     !
  ! full square real array, or collates a dense matrix across nodes to one     !
  ! single full square array                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination full square matrix                     !
  !   src    (input)  : The source block sparse matrix                         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  subroutine sparse_spam3tofull_complex(dest,src)

    use comms, only: comms_abort, comms_bcast, pub_on_root, &
         pub_total_num_nodes

    implicit none

    ! Arguments
    complex(kind=DP), intent(out) :: dest(:,:)  ! The full square matrix
    type(SPAM3), intent(in) :: src              ! The block sparse matrix

    ! Local variables
    integer :: srclib        ! Library entry for src
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! elem in atom block iblk counter
    integer :: jelemonat     ! elem in atom block jblk counter
    integer :: ielem,jelem   ! elem counters
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in this block/segment
    integer :: node          ! Node loop counter
    integer :: seg           ! Segment index
    integer :: seg_type      ! Start position of this segment in the index
    integer :: seg_start     ! Start position of this segment in the index

    ! Get library entry for src
    srclib = src%lib
    row_blks = library(srclib)%row_blks
    col_blks = library(srclib)%col_blks

    ! Check arguments
    if (size(dest,1) /= library(srclib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_spam3tofull_complex: &
            &incompatible numbers of rows in arguments.'
       call comms_abort
    end if
    if (size(dest,2) /= library(srclib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_spam3tofull_complex: &
            &incompatible numbers of cols in arguments.'
       call comms_abort
    end if
    if (.not.src%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_spam3tofull_complex: &
            &complex matrices only.'
       call comms_abort
    end if

    ! Zero full square matrix
    dest = 0.0_DP
    nrows = 0

    ! Loop over the segments of src on this node
    do seg=0,pub_total_num_nodes-1

       ! Cycle if this segment is blank
       seg_type = library(srclib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = num_elems_on_node(seg,row_blks)

       ! Loop over block-columns of src on this node
       seg_start = library(srclib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=my_first_blk,my_last_blk
          loc_iblk = loc_iblk + 1
          ielem = first_elem_on_atom(iblk,col_blks)
          ielems = num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of src
          do idx=library(srclib)%blk_idx(seg_start+loc_iblk-1), &
               library(srclib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(srclib)%blk_idx(idx)
             jelem = first_elem_on_atom(jblk,row_blks)
             jelems = num_elems_on_atom(jblk,row_blks)

             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ptr = library(srclib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                do jelemonat=0,jelems-1   ! Element row in block
                   dest(jelem+jelemonat,ielem+ielemonat) = &
                        src%zmtx(ptr+ielemonat*nrows+jelemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of src

       end do  ! Loop over block-columns of src

    end do  ! seg

    ! Broadcast result across all nodes
    do node=0,pub_total_num_nodes-1

       ! Get first column and number of columns
       ielem = first_elem_on_node(node,col_blks)
       ielems = num_elems_on_node(node,col_blks)

       ! Broadcast this part of the matrix
       call comms_bcast(node,dest(1,ielem),ielems*library(srclib)%nrows)

    end do  ! Loop over nodes

  end subroutine sparse_spam3tofull_complex

  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3 into a dense     !
  ! full square real array, or collates a dense matrix across nodes to one     !
  ! single full square array                                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination full square matrix                     !
  !   src    (input)  : The source block sparse matrix                         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  subroutine sparse_fulltospam3_complex(dest,src)

    use comms, only: comms_abort, pub_on_root, pub_total_num_nodes

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: dest       ! The block sparse matrix
    complex(kind=DP), intent(in) :: src(:,:) ! The full square matrix

    ! Local variables
    integer :: destlib       ! Library entry for dest
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! Element in atom block iblk counter
    integer :: jelemonat     ! Element in atom block jblk counter
    integer :: ielem,jelem   ! Element start positions of a block
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in the block/segment
    integer :: seg           ! Segment index
    integer :: seg_type      ! Segment type for this segment
    integer :: seg_start     ! Start position of this segment in the index

    ! Get library entry for src
    destlib = dest%lib
    row_blks = library(destlib)%row_blks
    col_blks = library(destlib)%col_blks
    nrows = 0

    ! Check arguments
    if (size(src,1) /= library(destlib)%nrows) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_fulltospam3_complex: &
            &incompatible numbers of rows in arguments.'
       call comms_abort
    end if
    if (size(src,2) /= library(destlib)%mcols) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_fulltospam3_complex: &
            &incompatible numbers of cols in arguments.'
       call comms_abort
    end if
    if (.not.dest%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_fulltospam3_complex: &
            &complex matrices only.'
       call comms_abort
    end if

    ! Loop over the segments of dest on this node
    do seg=0,pub_total_num_nodes-1

       ! Cycle if this segment is blank
       seg_type = library(destlib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = num_elems_on_node(seg,row_blks)

       ! Loop over block-columns of dest on this node
       seg_start = library(destlib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=my_first_blk,my_last_blk
          loc_iblk = loc_iblk + 1
          ielems = num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of dest
          do idx=library(destlib)%blk_idx(seg_start+loc_iblk-1), &
               library(destlib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(destlib)%blk_idx(idx)
             jelems = num_elems_on_atom(jblk,row_blks)
             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ielem = first_elem_on_atom(iblk,col_blks)
             ptr = library(destlib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                jelem = first_elem_on_atom(jblk,row_blks)
                do jelemonat=0,jelems-1   ! Element row in block
                   dest%zmtx(ptr+ielemonat*nrows+jelemonat) = &
                        src(jelem+jelemonat,ielem+ielemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of dest

       end do  ! Loop over block-columns of dest

    end do  ! seg

  end subroutine sparse_fulltospam3_complex

  !============================================================================!
  ! This subroutine scales and optionally shifts the eigenvalues of a matrix.  !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat    (inout) : The matrix to be rescaled and shifted                   !
  !   alpha  (input) : The complex scaling parameter                           !
  !   beta   (input) : The optional real shift parameter                       !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  !============================================================================!

  subroutine sparse_scale_complex(mat,alpha,beta)

    use comms, only: comms_abort, pub_on_root, pub_my_node_id
    implicit none

    ! Arguments
    type(SPAM3), intent(inout)          :: mat    ! The matrix to be operated on
    complex(kind=DP), intent(in)        :: alpha  ! The scaling parameter
    real(kind=DP), optional, intent(in) :: beta   ! The shift parameter

    ! Local variables
    integer :: ilib      ! Library entry for mat
    integer :: blks      ! Identifier for blocking scheme
    integer :: iblk      ! Block-column loop counter
    integer :: loc_iblk  ! Loop counter for iblk on this node
    integer :: ielems    ! Number of columns in block
    integer :: ptr       ! Pointer to diagonal entries
    integer :: ielem     ! Element in block-column loop counte
    integer :: seg       ! Segment index
    integer :: seg_start ! Segment type for this segment


    ! This only makes sense if mat is complex...
    if (.not. mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_scale_complex: &
            &mat must be complex.'
       call comms_abort
    end if

    ! Rescale the elements if required
    if (alpha /= (1.0_DP,0.0_DP)) then
       if (mat%iscmplx) then
          mat%zmtx = alpha * mat%zmtx
       else
          mat%dmtx = real(alpha * mat%dmtx,kind=DP)
       end if
    end if

    ! Shift the eigenvalues if required
    if (present(beta)) then

       ! Get library entry for mat
       ilib = mat%lib
       blks = library(ilib)%row_blks
       seg = pub_my_node_id

       ! Can only shift eigenvalues of square matrix
       if (library(ilib)%row_blks/=library(ilib)%col_blks) then
          write(stdout,'(a)') 'Error in sparse_scale: cannot shift &
               &eigenvalues of non-square matrices'
          call comms_abort
       end if

       if (library(ilib)%seg_info(s_type,seg)==SEG_DENSE) then

          ! Set pointer to first diagonal element
          ptr = library(ilib)%seg_info(s_ptr,seg)

          ! Get number of columns in this segment
          ielems = num_elems_on_node(seg,blks)

          ! Add identity matrix scaled by beta to all elements
          if (mat%iscmplx) then
             ! Loop over columns on my node
             do ielem=1,ielems
                mat%zmtx(ptr) = mat%zmtx(ptr) + cmplx(beta,0.0_DP,kind=DP)
                ptr = ptr + ielems + 1
             end do
          else
             ! Loop over columns on my node
             do ielem=1,num_elems_on_node(seg,blks)
                mat%dmtx(ptr) = mat%dmtx(ptr) + beta
                ptr = ptr + ielems + 1
             end do
          end if

       else if (library(ilib)%seg_info(s_type,seg)==SEG_SPARSE) then

          seg_start = library(ilib)%seg_info(s_idx,seg)

          ! Loop over atom-blocks iblk on my node
          loc_iblk = 0
          do iblk=my_first_blk,my_last_blk
             loc_iblk = loc_iblk + 1

             ! Get number of elements on this atom
             ielems = num_elems_on_atom(iblk,blks)

             ! Find pointer to start of diagonal block in this column
             ptr = library(ilib)%blk_ptr(seg_start+loc_iblk-1)

             ! Loop over diagonal elements in this diagonal block and add
             ! identity matrix scaled by beta
             if (mat%iscmplx) then
                do ielem=1,ielems
                   mat%zmtx(ptr) = mat%zmtx(ptr) + cmplx(beta,0.0_DP,kind=DP)
                   ptr = ptr + ielems + 1
                end do
             else
                do ielem=1,ielems
                   mat%dmtx(ptr) = mat%dmtx(ptr) + beta
                   ptr = ptr + ielems + 1
                end do
             end if

          end do  ! Loop over atom-blocks iblk on this node

       endif

    end if


  end subroutine sparse_scale_complex


  !============================================================================!
  ! This subroutine performs an axpy operation on the matrices:                !
  !   y := y + alpha * x                                                       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   ymat  (inout) : The sparse matrix y                                      !
  !   xmat  (input) : The sparse matrix x                                      !
  !   alpha (input) : The real parameter alpha                                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  !============================================================================!

  subroutine sparse_axpy_complex(ymat,xmat,alpha)

    use comms, only: comms_abort, pub_total_num_nodes

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: ymat    ! The sparse matrix y
    type(SPAM3), intent(in) :: xmat       ! The sparse matrix x
    complex(kind=DP), intent(in) :: alpha ! The parameter alpha

    ! Local variables
    integer :: xlib,ylib    ! Library pointers for x and y
    integer :: row_blks     ! Identifier for column blocking scheme
    integer :: col_blks     ! Identifier for row blocking scheme
    integer :: iblk,jblk    ! Atom loop counters
    integer :: icol         ! Atom loop counters
    integer :: loc_iblk     ! Atom loop counter on this node
    integer :: idx          ! Index loop counter
    integer :: iel          ! Element loop counter
    integer :: ifound       ! Loop counter over columns found
    integer :: nfound       ! Number of columns found
    integer :: xstart,xend  ! Start and end of block data of xmat
    integer :: ystart,yend  ! Start and end of block data of ymat
    integer :: ystride      ! Stride to get to next column of src
    integer :: xstride      ! Stride to get to next column of dest
    integer :: seg          ! Segment index
    integer :: x_seg_type   ! Segment type for segments of x
    integer :: y_seg_type   ! Segment type for segments of y
    integer :: x_seg_start  ! Start position of segments of x in the index
    integer :: y_seg_start  ! Start position of segments of y in the index


    ! Extract library entries
    xlib = xmat%lib
    ylib = ymat%lib
    row_blks = library(xlib)%row_blks
    col_blks = library(xlib)%col_blks
    xstride = 0; ystride = 0

    ! Check matrices share same blocking scheme
    if ( (library(xlib)%row_blks /= row_blks) .or. &
         (library(xlib)%col_blks /= col_blks) ) then
       write(stdout,'(a)') 'Error in sparse_axpy_complex: x and y matrices &
            &have different blocking schemes'
       call comms_abort
    end if

    ! Check for alpha = 0 (no operation required)
    if (alpha == (0.0_DP,0.0_DP)) then
       return
    end if

    ! Check whether the matrices share the same structure
    if (xlib == ylib) then

       ! Trivial axpy of whole data arrays
       if (xmat%iscmplx) then
          if (ymat%iscmplx) then
             ymat%zmtx = ymat%zmtx + alpha * xmat%zmtx
          else
             do iel=1,size(ymat%dmtx)
                ymat%dmtx(iel) = ymat%dmtx(iel) + &
                     alpha * real(xmat%zmtx(iel),kind=DP)
             end do
          end if
       else
          if (ymat%iscmplx) then
             do iel=1,size(ymat%zmtx)
                ymat%zmtx(iel) = ymat%zmtx(iel) + &
                     alpha * cmplx(xmat%dmtx(iel),0.0_DP,kind=DP)
             end do
          else
             ymat%dmtx = ymat%dmtx + real(alpha * xmat%dmtx,kind=DP)
          end if
       end if

    ! Structures differ, so proceed segment-by-segment
    else

       do seg=0,pub_total_num_nodes-1

          x_seg_type = library(xlib)%seg_info(s_type,seg)
          y_seg_type = library(ylib)%seg_info(s_type,seg)

          ! Nothing to add to if this segment in y is blank
          if (y_seg_type==SEG_BLANK) cycle

          ! Nothing to add if this segment in x is blank
          if (x_seg_type == SEG_BLANK) cycle

          ! Find information about these segments
          x_seg_start = library(xlib)%seg_info(s_idx,seg)
          y_seg_start = library(ylib)%seg_info(s_idx,seg)
          if (y_seg_type == SEG_DENSE) ystride = &
               num_elems_on_node(seg,row_blks)
          if (x_seg_type == SEG_DENSE) xstride = &
               num_elems_on_node(seg,row_blks)

          ! Loop over atom blocks (block-columns) on this node
          loc_iblk = 0
          do iblk=my_first_blk,my_last_blk
             loc_iblk = loc_iblk + 1

             ! Reset counters of number of rows found
             nfound = 0

             ! Loop over block-rows in x for this segment
             do idx=library(xlib)%blk_idx(x_seg_start+loc_iblk-1), &
                  library(xlib)%blk_idx(x_seg_start+loc_iblk)-1
                jblk = library(xlib)%blk_idx(idx)

                ! Mark flag for this row and add to list
                nfound = nfound + 1
                found_blk(jblk) = .true.
                found_idx(nfound) = jblk
                found_ptr(jblk) = library(xlib)%blk_ptr(idx)

             end do

             ! Loop over block-rows in y
             do idx=library(ylib)%blk_idx(y_seg_start+loc_iblk-1), &
                  library(ylib)%blk_idx(y_seg_start+loc_iblk)-1
                jblk = library(ylib)%blk_idx(idx)

                ! Check whether this block-row exists in source
                if (found_blk(jblk)) then

                   ! Find pointers to start of dest and src data
                   ystart = library(ylib)%blk_ptr(idx)
                   yend = ystart + num_elems_on_atom(jblk,row_blks) - 1
                   xstart = found_ptr(jblk)
                   xend = xstart + num_elems_on_atom(jblk,row_blks) - 1

                   ! Find stride to next column of dest data
                   if (y_seg_type == SEG_SPARSE) then
                      ystride = num_elems_on_atom(jblk,row_blks)
                   end if

                   ! Find stride to next column of source data
                   if (x_seg_type == SEG_SPARSE) then
                      xstride = num_elems_on_atom(jblk,row_blks)
                   end if

                   ! Copy columns from source to destination
                   do icol=1,num_elems_on_atom(iblk,col_blks)

                      if (xmat%iscmplx) then
                         if (ymat%iscmplx) then
                            ymat%zmtx(ystart:yend) = &
                                 ymat%zmtx(ystart:yend) + &
                                 alpha * xmat%zmtx(xstart:xend)
                         else
                            do iel=ystart,yend
                               ymat%dmtx(iel) = ymat%dmtx(iel) + real( &
                                    alpha * xmat%zmtx(iel-ystart+xstart), &
                                    kind=DP)
                            end do
                         end if
                      else
                         if (ymat%iscmplx) then
                            do iel=ystart,yend
                               ymat%zmtx(iel) = ymat%zmtx(iel) + alpha * &
                                    cmplx(xmat%dmtx(iel-ystart+xstart), &
                                    0.0_DP,kind=DP)
                            end do
                         else
                            ymat%dmtx(ystart:yend) = &
                                 real(ymat%dmtx(ystart:yend) + &
                                 alpha * xmat%dmtx(xstart:xend),kind=DP)
                         end if
                      end if

                      ystart = ystart + ystride
                      yend = yend + ystride
                      xstart = xstart + xstride
                      xend = xend + xstride

                   end do

                end if

             end do  ! Loop over destination block-rows

             ! Reset flags
             do ifound=1,nfound
                jblk = found_idx(ifound)
                found_blk(jblk) = .false.
             end do

          end do  ! Loop over atom block-columns

       end do  ! seg

    end if


  end subroutine sparse_axpy_complex


  !============================================================================!
  ! This subroutine copies the data from one sparse matrix to another, taking  !
  ! different structures into account if necessary.                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest (inout) : The destination sparse matrix                             !
  !   src  (input) : The source sparse matrix                                  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009 based on SPAM2 version.                 !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Minor Modification for dense matrices by Nicholas Hine, Dec 2007.          !
  !============================================================================!

  subroutine sparse_copy(dest,src)

    use comms, only: comms_abort, pub_total_num_nodes

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: dest   ! The destination sparse matrix
    type(SPAM3), intent(in) :: src       ! The source sparse matrix

    ! Local variables
    integer :: destlib        ! Library entry for dest
    integer :: srclib         ! Library entry for src
    integer :: row_blks       ! Identifier for column blocking scheme
    integer :: col_blks       ! Identifier for row blocking scheme
    integer :: iblk,jblk      ! Atom loop counters
    integer :: loc_iblk       ! Atom loop counter on this node
    integer :: idx            ! Index loop counter
    integer :: iel            ! Element loop counter
    integer :: ifound         ! Loop counter over columns found
    integer :: nfound         ! Number of columns found
    integer :: icol           ! Column element loop counter
    integer :: dstart,dend    ! Start and end of block data of dest
    integer :: sstart,send    ! Start and end of block data of src
    integer :: dstride        ! Stride to get to next column of src
    integer :: sstride        ! Stride to get to next column of dest
    integer :: seg            ! Segment index
    integer :: src_seg_type   ! Segment type for segments of src
    integer :: dest_seg_type  ! Segment type for segments of dest
    integer :: src_seg_start  ! Start position of segments of src in the index
    integer :: dest_seg_start ! Start position of segments of dest in the index

    ! Extract library entries
    destlib = dest%lib
    srclib = src%lib
    row_blks = library(srclib)%row_blks
    col_blks = library(srclib)%col_blks

    ! Check matrices share same blocking scheme
    if ( (library(destlib)%row_blks /= row_blks) .or. &
         (library(destlib)%col_blks /= col_blks) ) then
       write(stdout,'(a)') 'Error in sparse_copy: src and dest matrices &
            &have different blocking schemes'
       call comms_abort
    end if

    ! Check whether the matrices share the same structure
    if ( srclib == destlib ) then

       ! Trivial copy of data arrays
       if (src%iscmplx) then
          if (dest%iscmplx) then
             dest%zmtx = src%zmtx
          else
             do iel=1,size(dest%dmtx)
                dest%dmtx(iel) = real(src%zmtx(iel),kind=DP)
             end do
          end if
       else
          if (dest%iscmplx) then
             do iel=1,size(dest%zmtx)
                dest%zmtx(iel) = cmplx(src%dmtx(iel),0.0_DP,kind=DP)
             end do
          else
             dest%dmtx = src%dmtx
          end if
       end if

    else

       do seg=0,pub_total_num_nodes-1

          src_seg_type = library(srclib)%seg_info(s_type,seg)
          dest_seg_type = library(destlib)%seg_info(s_type,seg)
          sstart = library(srclib)%seg_info(s_ptr,seg)
          send = library(srclib)%seg_info(s_ptr,seg+1) - 1
          dstart = library(destlib)%seg_info(s_ptr,seg)
          dend = library(destlib)%seg_info(s_ptr,seg+1) - 1

          ! Nothing to copy if destination segment is blank
          if (dest_seg_type==SEG_BLANK) cycle

          ! Fill dest data with zeros if source segment is blank
          if (src_seg_type == SEG_BLANK) then

             if (dest%iscmplx) then
                dest%zmtx(dstart:dend) = cmplx(0.0_DP,0.0_DP,kind=DP)
             else
                dest%dmtx(dstart:dend) = 0.0_DP
             end if

             cycle

          end if

          src_seg_start = library(srclib)%seg_info(s_idx,seg)
          dest_seg_start = library(destlib)%seg_info(s_idx,seg)

          ! Loop over block-columns on this node
          loc_iblk = 0
          do iblk=my_first_blk,my_last_blk
             loc_iblk = loc_iblk + 1

             ! Reset counters of number of rows found
             nfound = 0

             ! Loop over block-rows in source for this segment
             do idx=library(srclib)%blk_idx(src_seg_start+loc_iblk-1), &
                  library(srclib)%blk_idx(src_seg_start+loc_iblk)-1
                jblk = library(srclib)%blk_idx(idx)

                ! Mark flag for this row and add to list
                nfound = nfound + 1
                found_blk(jblk) = .true.
                found_idx(nfound) = jblk
                found_ptr(jblk) = library(srclib)%blk_ptr(idx)

             end do  ! idx

             ! Loop over block-rows in destination
             do idx=library(destlib)%blk_idx(dest_seg_start+loc_iblk-1), &
                  library(destlib)%blk_idx(dest_seg_start+loc_iblk)-1
                jblk = library(destlib)%blk_idx(idx)

                ! Check whether this block-row exists in source
                if (found_blk(jblk)) then

                   ! Find pointers to start of dest and src data
                   dstart = library(destlib)%blk_ptr(idx)
                   dend = dstart + num_elems_on_atom(jblk,row_blks) - 1
                   sstart = found_ptr(jblk)
                   send = sstart + num_elems_on_atom(jblk,row_blks) - 1

                   ! Find stride to next column of dest data
                   if (dest_seg_type == SEG_DENSE) then
                      dstride = num_elems_on_node(seg,row_blks)
                   else
                      dstride = num_elems_on_atom(jblk,row_blks)
                   end if

                   ! Find stride to next column of source data
                   if (src_seg_type == SEG_DENSE) then
                      sstride = num_elems_on_node(seg,row_blks)
                   else
                      sstride = num_elems_on_atom(jblk,row_blks)
                   end if

                   ! Copy columns from source to destination
                   do icol=1,num_elems_on_atom(iblk,col_blks)

                      if (src%iscmplx) then
                         if (dest%iscmplx) then
                            dest%zmtx(dstart:dend) = src%zmtx(sstart:send)
                         else
                            do iel=dstart,dend
                               dest%dmtx(iel) = &
                                    real(src%zmtx(iel-dstart+sstart),kind=DP)
                            end do
                         end if
                      else
                         if (dest%iscmplx) then
                            do iel=dstart,dend
                               dest%zmtx(iel) = cmplx(src%dmtx(iel-&
                                    dstart+sstart),0.0_DP,kind=DP)
                            end do
                         else
                            dest%dmtx(dstart:dend) = src%dmtx(sstart:send)
                         end if
                      end if

                      dstart = dstart + dstride
                      dend = dend + dstride
                      sstart = sstart + sstride
                      send = send + sstride

                   end do

                else

                   ! Set data for this block to zero
                   dstart = library(destlib)%blk_ptr(idx)
                   dend = dstart + num_elems_on_atom(jblk,row_blks) - 1

                   if (dest_seg_type == SEG_DENSE) then
                      dstride = num_elems_on_node(seg,row_blks)
                   else
                      dstride = num_elems_on_atom(jblk,row_blks)
                   end if

                   do icol=1,num_elems_on_atom(iblk,col_blks)

                      if (dest%iscmplx) then
                         dest%zmtx(dstart:dend) = (0.0_DP,0.0_DP)
                      else
                         dest%dmtx(dstart:dend) = 0.0_DP
                      end if

                      dstart = dstart + dstride
                      dend = dend + dstride

                   end do  ! icol

                end if

             end do  ! Loop over destination block-rows

             ! Reset flags
             do ifound=1,nfound
                jblk = found_idx(ifound)
                found_blk(jblk) = .false.
             end do

          end do  ! Loop over atom block-columns

       end do  ! seg

    end if

  end subroutine sparse_copy


  !============================================================================!
  ! This subroutine performs a matmul operation on the matrices:               !
  !   C := A . B                                                               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   cmat  (inout) : The sparse matrix C                                      !
  !   amat  (input) : The sparse matrix A                                      !
  !   bmat  (input) : The sparse matrix B                                      !
  !----------------------------------------------------------------------------!
  ! SPAM3 version written by Nicholas Hine, June 2009.                         !
  ! Allows mixed sparse-dense matrices divided segment-by-segment.             !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Modified by Nicholas Hine, November 2007 to loop-unroll block matrix       !
  ! multiplication in the inner loop.                                          !
  ! Modified by Nicholas Hine, January 2008 to remove blocking comms.          !
  ! Re-written by Nicholas Hine, November 2008 to improve comms efficiency on  !
  ! very large numbers of nodes.                                               !
  !============================================================================!

  subroutine sparse_product(cmat,amat,bmat)

    use comms, only: comms_abort, comms_send, comms_free, comms_barrier, &
         pub_on_root, pub_my_node_id, pub_total_num_nodes
    use comms, only: comms_allgather, comms_reduce, pub_comms_group_size, &
         pub_group_comm, pub_my_rank_in_group, pub_my_comms_group, &
         pub_num_comms_groups, pub_first_node_in_group, pub_last_node_in_group
    use parallel_strategy, only: pub_first_atom_on_node
    use rundat, only: timings_level
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: cmat   ! The sparse matrix C
    type(SPAM3), intent(in) :: amat      ! The sparse matrix A
    type(SPAM3), intent(in) :: bmat      ! The sparse matrix B

!CW
    integer :: ii,i,j,k
!END CW

    ! Local variables
    integer :: alib,blib,clib   ! Pointers to library entries
    integer :: arow_blks        ! Identifier for A row, C row blocking scheme
    integer :: bcol_blks        ! Identifier for B col, C col blocking scheme
    integer :: acol_blks        ! Identifier for A col, B row blocking scheme
    integer :: iblk,jblk,kblk   ! Atom block loop counters
    integer :: recvnode         ! Node index of recv origin
    integer :: loc_iblk         ! Atom block loop counter local to this node
    integer :: loc_jblk         ! Atom block loop counter local to this node
    integer :: iidx,jidx        ! Indices
    integer :: ierr             ! Status flag
    integer :: step             ! Node step counter
    integer :: nrows            ! Number of rows of A and C in block/segment
    integer :: ncols            ! Number of cols of B and C in block/segment
    integer :: nsum             ! Number of cols of A and rows of B
    integer :: aptr,bptr,cptr   ! Pointer to blocks/segs in matrices A, B and C
    integer :: ielem,jelem      ! Element row/col counters
    integer :: ielem0,jelem0    ! Element row/col block start positions
    integer :: ielem1,jelem1    ! Element row/col block start positions
    integer :: ac_row_seg       ! Segment index of row segment in A and C
    integer :: bc_col_seg       ! Segment index of col segment in B and C
    integer :: bc_col_seg_node  ! Node corresponding to bc_col_seg
    integer :: aseg_start       ! Start position in A index of current segment
    integer :: bseg_start       ! Start position in B index of current segment
    integer :: cseg_start       ! Start position in C index of current segment
    integer :: lda              ! Leading dimension of array storing block of A
    integer :: ldb              ! Leading dimension of array storing block of B
    integer :: ldc              ! Leading dimension of array storing block of C
    integer :: ifound,nfound    ! Index and number of found blocks of C
    logical :: a_dense,a_sparse ! Shorthand flags for type of segment of A
    logical :: b_dense,b_sparse ! Shorthand flags for type of segment of B
    logical :: c_dense,c_sparse ! Shorthand flags for type of segment of C
    integer :: reqnode          ! Node to request from
    type(COM3) :: acom
    integer, allocatable :: bseginfobuf(:,:,:) ! Recv buffer for segment info
    integer, allocatable :: bidxbuf(:,:)   ! Recv buffer for index
    integer, allocatable :: bptrbuf(:,:)   ! Recv buffer for ptr list
    real(kind=DP), allocatable :: bdmtxbuf(:,:)   ! Buffer for real data
    complex(kind=DP), allocatable :: bzmtxbuf(:,:)! Buffer for cplx data
    integer, allocatable :: cseginfobuf(:,:,:) ! Recv buffer for segment info
    integer, allocatable :: cidxbuf(:,:)   ! Recv buffer for index
    integer, allocatable :: cptrbuf(:,:)   ! Recv buffer for ptr list
    real(kind=DP), allocatable :: cdmtxbuf(:,:)   ! Buffer for real data
    complex(kind=DP), allocatable :: czmtxbuf(:,:)! Buffer for cplx data

    ! Extract library entries
    alib = amat%lib
    blib = bmat%lib
    clib = cmat%lib
    arow_blks = library(alib)%row_blks
    bcol_blks = library(blib)%col_blks
    acol_blks = library(alib)%col_blks

    ! Special timer level for measuring individual sparse_product calls
    if (timings_level==9) then
       call timer_clock('sparse_product_' &
            //trim(library(alib)%structure)//'_' &
            //trim(library(blib)%structure)//'_' &
            //trim(library(clib)%structure),1)
    end if
#ifdef ITC_TRACE
    call VTBEGIN(vt_sparse_product,vt_err)
#endif

    ! Check consistency of matrix types and blocking schemes
    if ( (acol_blks /= library(blib)%row_blks).or.&
         (arow_blks /= library(clib)%row_blks).or.&
         (bcol_blks /= library(clib)%col_blks) ) then
       if (pub_on_root) write(stdout,'(a)') &
            'Error in sparse_product: incompatible blocking schemes.'
       call comms_abort
    end if
    if ((amat%iscmplx .neqv. bmat%iscmplx) .or. &
         (cmat%iscmplx .neqv. bmat%iscmplx)) then
       if (pub_on_root) write(stdout,'(a)') &
            'Error in sparse_product: incompatible argument types.'
       call comms_abort
    end if

    ! Allocate arrays and initialise comms
    call sparse_com_allocate(acom,amat,2, &
         alloc_mtx=.true.,cropped=.true.,seg=.false.)

    ! Set up which other nodes will be requesting data
    do recvnode=0,pub_total_num_nodes-1
       if (library(blib)%idx_seg_lens(recvnode)/=0) then
          acom%index_reqs(recvnode) = index_not_sent
          acom%data_reqs(recvnode) = data_not_sent
       end if
    end do

    ! Initialise B, C data sharing
    call internal_init_shares

    ! Initialise amat comms
    call sparse_init_comms(acom)

    ! The product to be formed is:
    !  C   = A   B
    !   ki    kj  ji

    ! Loop over processors
    do step=0,pub_num_comms_groups-1

       ! Receive index, pointers and data for this step if required
       recvnode = modulo(pub_my_node_id-step*pub_comms_group_size+pub_total_num_nodes, &
            pub_total_num_nodes)
       if (any(bseginfobuf(s_type,recvnode,0:pub_comms_group_size-1)/=SEG_BLANK)) then
          call sparse_get_crop_step_data(acom,amat,recvnode, &
               maxval(library(blib)%idx_lens),maxval(library(clib)%idx_lens), &
               bidxbuf,cidxbuf,bseginfobuf,cseginfobuf)
       end if

       ! Request index and data for next step if required
       reqnode = modulo(pub_my_node_id-(step+1)*pub_comms_group_size+pub_total_num_nodes, &
            pub_total_num_nodes)
       if (any(bseginfobuf(s_type,reqnode,0:pub_comms_group_size-1)/=SEG_BLANK).and. &
            (reqnode/=pub_my_node_id)) then

          ! Send request for index and pointers
          call comms_send(reqnode,index_needed,1,tag=IDX_REQ_TAG)

          ! Start asynchronous receives for index and pointers
          call sparse_recv_index(acom,reqnode,1,.false.,async=.true.)
          call sparse_recv_pointers(acom,reqnode,1,.false.,async=.true.)
       end if

       ! Loop over col segments of B and C in this group of nodes
       do bc_col_seg=0,pub_comms_group_size-1

          ! Find the node index in the full list for this bc_col_seg
          bc_col_seg_node = bc_col_seg + pub_first_node_in_group

          ! No work for this recvnode if corresponding segment of B is blank
          if (bseginfobuf(s_type,recvnode,bc_col_seg)==SEG_BLANK) cycle

          ! Find out if corresponding section of B is sparse or dense
          b_sparse = (bseginfobuf(s_type,recvnode,bc_col_seg)==SEG_SPARSE)
          b_dense = (bseginfobuf(s_type,recvnode,bc_col_seg)==SEG_DENSE)

          ! Loop over (row-) segments of A and C
          do ac_row_seg=0,pub_total_num_nodes-1

             ! Skip this segment if A or C are blank
             if ((acom%seginfobuf(s_type,ac_row_seg,2)==SEG_BLANK).or.&
                  (cseginfobuf(s_type,ac_row_seg,bc_col_seg)==SEG_BLANK)) cycle

             ! Find out information about these segments of A and C
             a_sparse = (acom%seginfobuf(s_type,ac_row_seg,2)==SEG_SPARSE)
             a_dense = (acom%seginfobuf(s_type,ac_row_seg,2)==SEG_DENSE)
             c_sparse = (cseginfobuf(s_type,ac_row_seg,bc_col_seg)==SEG_SPARSE)
             c_dense = (cseginfobuf(s_type,ac_row_seg,bc_col_seg)==SEG_DENSE)
             nrows = num_elems_on_node(ac_row_seg,arow_blks) ! rows of A and C
             ncols = num_elems_on_node(bc_col_seg_node,bcol_blks) ! cols of C
             nsum = num_elems_on_node(recvnode,acol_blks) ! cols of A
             aptr = acom%seginfobuf(s_ptr,ac_row_seg,2)
             bptr = bseginfobuf(s_ptr,recvnode,bc_col_seg)
             cptr = cseginfobuf(s_ptr,ac_row_seg,bc_col_seg)

             if (a_dense.and.b_dense.and.c_dense) then
#ifdef ITC_TRACE
                call VTBEGIN(vt_sparse_product_dense,vt_err)
#endif
                ! Multiply segments of A and B with Lapack call
                if (cmat%iscmplx) then
                   call zgemm('N','N',nrows,ncols,nsum,(1.0_DP,0.0_DP), &
                        acom%zmtxrecvbuf(aptr,2),nrows,&
                        bzmtxbuf(bptr,bc_col_seg),nsum,&
                        (1.0_DP,0.0_DP),czmtxbuf(cptr,bc_col_seg),nrows)
                else
                   call dgemm('N','N',nrows,ncols,nsum,1.0_DP, &
                        acom%dmtxrecvbuf(aptr,2),nrows,&
                        bdmtxbuf(bptr,bc_col_seg),nsum,&
                        1.0_DP,cdmtxbuf(cptr,bc_col_seg),nrows)
                end if
#ifdef ITC_TRACE
                call VTEND(vt_sparse_product_dense,vt_err)
#endif
             else if (a_dense.and.b_dense.and.c_sparse) then
#ifdef ITC_TRACE
                call VTBEGIN(vt_sparse_product_picks,vt_err)
#endif
                ! Loop over elements of C only, calculating product explicitly
                cseg_start = cseginfobuf(s_idx,ac_row_seg,bc_col_seg)

                ! Loop over col-blocks of C
                loc_iblk = 0
                do iblk=pub_first_atom_on_node(bc_col_seg_node), &
                     pub_first_atom_on_node(bc_col_seg_node+1)-1
                   loc_iblk = loc_iblk + 1

                   ! Find range of col offsets in this block of C
                   ielem0 = first_elem_on_atom(iblk,bcol_blks) - &
                        first_elem_on_node(bc_col_seg_node,bcol_blks)
                   ielem1 = ielem0 + num_elems_on_atom(iblk,bcol_blks) - 1

                   ! Loop over block-rows jblk in block-column iblk of C
                   do iidx=cidxbuf(cseg_start+loc_iblk-1,bc_col_seg), &
                        cidxbuf(cseg_start+loc_iblk,bc_col_seg)-1
                      jblk = cidxbuf(iidx,bc_col_seg)
                      cptr = cptrbuf(iidx,bc_col_seg)

                      ! Find range of row offsets in this block of C
                      jelem0 = first_elem_on_atom(jblk,arow_blks) - &
                           first_elem_on_node(ac_row_seg,arow_blks)
                      jelem1 = jelem0 + num_elems_on_atom(jblk,arow_blks) - 1

                      ! Loop over elements of this block of C
                      if (cmat%iscmplx) then
                         do ielem=ielem0,ielem1
                            do jelem=jelem0,jelem1
                               czmtxbuf(cptr,bc_col_seg) = &
                                    czmtxbuf(cptr,bc_col_seg) + sum( &
                                    bzmtxbuf(bptr+ielem*nsum:bptr+(ielem+1)*nsum-1:1,bc_col_seg) * &
                                    acom%zmtxrecvbuf(aptr+jelem:aptr+jelem+nrows*(nsum-1):nrows,2))
                               cptr = cptr + 1
                            end do
                         end do
                      else
                        do ielem=ielem0,ielem1
                            do jelem=jelem0,jelem1

                               cdmtxbuf(cptr,bc_col_seg) = cdmtxbuf(cptr,bc_col_seg) + &
                                    sum(bdmtxbuf(bptr+ielem*nsum:bptr+(ielem+1)*nsum-1:1,bc_col_seg) * &
                                    acom%dmtxrecvbuf(aptr+jelem:aptr+jelem+nrows*(nsum-1):nrows,2))
                               cptr = cptr + 1
                            end do
                         end do
                      end if

                   end do  ! iidx
                end do  ! iblk
#ifdef ITC_TRACE
                call VTEND(vt_sparse_product_picks,vt_err)
#endif
             else ! any other combination of sparse and dense segments
#ifdef ITC_TRACE
                call VTBEGIN(vt_sparse_product_sparse,vt_err)
#endif
                ! Set up shorthand variables for these segments
                aseg_start = acom%seginfobuf(s_idx,ac_row_seg,2)
                bseg_start = bseginfobuf(s_idx,recvnode,bc_col_seg)
                cseg_start = cseginfobuf(s_idx,ac_row_seg,bc_col_seg)
                if (a_dense) lda = num_elems_on_node(ac_row_seg,arow_blks)
                if (b_dense) ldb = num_elems_on_node(recvnode,acol_blks)
                if (c_dense) ldc = num_elems_on_node(ac_row_seg,arow_blks)

                ! Loop over block-columns iblk of B and C on this node
                loc_iblk = 0
                do iblk=pub_first_atom_on_node(bc_col_seg_node), &
                     pub_first_atom_on_node(bc_col_seg_node+1)-1
                   loc_iblk = loc_iblk + 1

                   ! Set number of columns in this block-column of C
                   ncols = num_elems_on_atom(iblk,bcol_blks)

                   ! Only proceed if there are nonzero block-rows in this column
                   ! in this segment of B
                   if (bidxbuf(bseg_start+loc_iblk-1,bc_col_seg)== &
                        bidxbuf(bseg_start+loc_iblk,bc_col_seg)) cycle

                   ! Set number of blocks found to zero
                   nfound = 0

                   ! Loop over block-rows kblk of C in column iblk
                   do iidx=cidxbuf(cseg_start+loc_iblk-1,bc_col_seg), &
                        cidxbuf(cseg_start+loc_iblk,bc_col_seg)-1
                      kblk = cidxbuf(iidx,bc_col_seg)

                      ! Mark this block-row as found and note row and pointer
                      found_blk(kblk) = .true.
                      nfound = nfound + 1
                      found_idx(nfound) = kblk
                      found_ptr(kblk) = cptrbuf(iidx,bc_col_seg)

                   end do  ! Loop over block-rows kblk of C in column iblk

                   ! Don't check B index if C has no nze's in this col of
                   ! this A-row segment
                   if (nfound==0) cycle

                   ! Loop over block-rows jblk of B in column iblk
                   do iidx=bidxbuf(bseg_start+loc_iblk-1,bc_col_seg), &
                        bidxbuf(bseg_start+loc_iblk,bc_col_seg)-1
                      jblk = bidxbuf(iidx,bc_col_seg)

                      ! Find number of this atom on recvnode
                      loc_jblk = jblk - pub_first_atom_on_node(recvnode) + 1

                      ! Set number of rows in this block-row of B
                      nsum = num_elems_on_atom(jblk,acol_blks)
                      if (b_sparse) ldb = nsum

                      ! Copy pointer to this block in B
                      bptr = bptrbuf(iidx,bc_col_seg)

                      ! Loop over block-rows kblk of A in column jblk
                      do jidx=acom%idxbuf(aseg_start+loc_jblk-1,2), &
                           acom%idxbuf(aseg_start+loc_jblk,2)-1
                         kblk = acom%idxbuf(jidx,2)
                         ! If this row occurs in C do block multiply
                         if (found_blk(kblk)) then
                            nrows = num_elems_on_atom(kblk,arow_blks)
                            if (a_sparse) lda = nrows
                            if (c_sparse) ldc = nrows
                            cptr = found_ptr(kblk)
                            aptr = acom%ptrbuf(jidx,2)

                            if (amat%iscmplx) then
                               if ((nrows==1.or.nrows==4.or.nrows==9).and. &
                                   (ncols==1.or.ncols==4.or.ncols==9).and. &
                                   (nsum ==1.or.nsum ==4.or.nsum ==9)) then
                                  call sparse_zgemm_149(nrows,ncols,nsum, &
                                       acom%zmtxrecvbuf(aptr,2),lda, &
                                       bzmtxbuf(bptr,bc_col_seg),ldb, &
                                       czmtxbuf(cptr,bc_col_seg),ldc)
                               else
                                  call zgemm('N','N',nrows,ncols,nsum, &
                                       (1.0_DP,0.0_DP), &
                                       acom%zmtxrecvbuf(aptr,2),lda, &
                                       bzmtxbuf(bptr,bc_col_seg),ldb,(1.0_DP,0.0_DP), &
                                       czmtxbuf(cptr,bc_col_seg),ldc)
                               end if
                            else
                               if ((nrows==1.or.nrows==4.or.nrows==9).and. &
                                   (ncols==1.or.ncols==4.or.ncols==9).and. &
                                   (nsum ==1.or.nsum ==4.or.nsum ==9)) then
                                  call sparse_dgemm_149(nrows,ncols,nsum, &
                                       acom%dmtxrecvbuf(aptr,2),lda, &
                                       bdmtxbuf(bptr,bc_col_seg),ldb, &
                                       cdmtxbuf(cptr,bc_col_seg),ldc)
                               else
                                  call dgemm('N','N',nrows,ncols,nsum,1.0_DP, &
                                       acom%dmtxrecvbuf(aptr,2),lda, &
                                       bdmtxbuf(bptr,bc_col_seg),ldb,1.0_DP, &
                                       cdmtxbuf(cptr,bc_col_seg),ldc)
                               end if
                            end if

                         end if

                      end do  ! Column blocks kblk of A in row jblk

                   end do  ! Block-rows jblk of B in column iblk

                   ! Reset found block-rows flags
                   do ifound=1,nfound
                      iidx = found_idx(ifound)
                      found_blk(iidx) = .false.
                   end do

                end do  ! Loop over block-columns iblk of B and C
#ifdef ITC_TRACE
                call VTEND(vt_sparse_product_sparse,vt_err)
#endif
             end if

             ! Check incoming requests after each non-skipped segment
             call sparse_check_send_requests(acom,amat)

          end do  ! ac_row_seg

       end do  ! bc_col_seg

    end do  ! step

    ! Wait for final sends and deallocate communication buffers
    call sparse_exit_comms(acom,amat)
    call sparse_com_deallocate(acom,dealloc_mtx=.true.,cropped=.true.)
    call internal_exit_shares

    call comms_barrier

    if (any(library(clib)%seg_info(s_type,:)==SEG_DENSE)) then
       call sparse_enforce_sparsity(cmat)
    end if

    ! ndmh: special timer level for measuring individual sparse_product calls
    if (timings_level==9) then
       call timer_clock('sparse_product_' &
            //trim(library(alib)%structure)//'_' &
            //trim(library(blib)%structure)//'_' &
            //trim(library(clib)%structure),2)
    end if
#ifdef ITC_TRACE
    call VTEND(vt_sparse_product,vt_err)
#endif

contains

    ! This subroutine allocates temporary arrays required for the
    ! sharing of the B and C matrix data between a group
    ! Written by Nicholas Hine in July 2011
    subroutine internal_init_shares

      implicit none

      integer :: sendnode
      integer :: buflen, datlen
      integer, allocatable :: displs(:), lengths(:)

      ! Allocate buffers for bmat data
      allocate(bseginfobuf(1:3,0:pub_total_num_nodes, &
           0:pub_comms_group_size-1),stat=ierr)
      call utils_alloc_check('sparse_product','bseginfobuf',ierr)
      buflen = maxval(library(blib)%idx_lens)
      allocate(bidxbuf(1:buflen,0:pub_comms_group_size-1),stat=ierr)
      call utils_alloc_check('sparse_product','bidxbuf',ierr)
      allocate(bptrbuf(1:buflen,0:pub_comms_group_size-1),stat=ierr)
      call utils_alloc_check('sparse_product','bptrbuf',ierr)
      datlen = library(blib)%max_nze
      if (bmat%iscmplx) then
         allocate(bzmtxbuf(datlen,0:pub_comms_group_size-1),stat=ierr)
         call utils_alloc_check('sparse_product','bzmtxbuf',ierr)
         allocate(bdmtxbuf(1,1),stat=ierr)
         call utils_alloc_check('sparse_product','bdmtxbuf',ierr)
      else
         allocate(bdmtxbuf(datlen,0:pub_comms_group_size-1),stat=ierr)
         call utils_alloc_check('sparse_product','bdmtxbuf',ierr)
         allocate(bzmtxbuf(1,1),stat=ierr)
         call utils_alloc_check('sparse_product','bzmtxbuf',ierr)
      end if

      ! Allocate buffers for cmat data
      allocate(cseginfobuf(1:3,0:pub_total_num_nodes, &
           0:pub_comms_group_size-1),stat=ierr)
      call utils_alloc_check('sparse_product','cseginfobuf',ierr)
      buflen = maxval(library(clib)%idx_lens)
      allocate(cidxbuf(1:buflen,0:pub_comms_group_size-1),stat=ierr)
      call utils_alloc_check('sparse_product','cidxbuf',ierr)
      allocate(cptrbuf(1:buflen,0:pub_comms_group_size-1),stat=ierr)
      call utils_alloc_check('sparse_product','cptrbuf',ierr)
      datlen = library(clib)%max_nze
      if (cmat%iscmplx) then
         allocate(czmtxbuf(datlen,0:pub_comms_group_size-1),stat=ierr)
         call utils_alloc_check('sparse_product','czmtxbuf',ierr)
         allocate(cdmtxbuf(1,1),stat=ierr)
         call utils_alloc_check('sparse_product','cdmtxbuf',ierr)
      else
         allocate(cdmtxbuf(datlen,0:pub_comms_group_size-1),stat=ierr)
         call utils_alloc_check('sparse_product','cdmtxbuf',ierr)
         allocate(czmtxbuf(1,1),stat=ierr)
         call utils_alloc_check('sparse_product','czmtxbuf',ierr)
      end if

      ! Zero result matrix C
      if (cmat%iscmplx) then
         czmtxbuf = (0.0_DP,0.0_DP)
      else
         cdmtxbuf = 0.0_DP
      end if

      ! Length and displacement buffers for allgatherv
      allocate(displs(0:pub_comms_group_size-1),stat=ierr)
      call utils_alloc_check('sparse_product','displs',ierr)
      allocate(lengths(0:pub_comms_group_size-1),stat=ierr)
      call utils_alloc_check('sparse_product','lengths',ierr)

      ! Share index, pointers and data for bmat over all nodes
      call comms_allgather(bseginfobuf(1:3,0:pub_total_num_nodes,0), &
           library(blib)%seg_info(1:3,0:pub_total_num_nodes), &
           comm=pub_group_comm)
      lengths(0:pub_comms_group_size-1) = library(blib)%idx_lens( &
           pub_first_node_in_group:pub_last_node_in_group)
      do sendnode=0,pub_comms_group_size-1
         displs(sendnode) = sendnode*size(bidxbuf,1)
      end do
      buflen = max(library(blib)%idx_lens(pub_my_node_id),1)
      call comms_allgather(bidxbuf(:,0),library(blib)%blk_idx(1:buflen), &
           length_src=lengths(pub_my_rank_in_group), &
           lengths_dest=lengths,displs_dest=displs,comm=pub_group_comm)
      call comms_allgather(bptrbuf(:,0),library(blib)%blk_ptr(1:buflen), &
           length_src=lengths(pub_my_rank_in_group), &
           lengths_dest=lengths,displs_dest=displs,comm=pub_group_comm)
      lengths(0:pub_comms_group_size-1) = &
           bseginfobuf(s_ptr,pub_total_num_nodes,0:pub_comms_group_size-1) - 1
      do sendnode=0,pub_comms_group_size-1
         if (bmat%iscmplx) then
            displs(sendnode) = sendnode*size(bzmtxbuf,1)
         else
            displs(sendnode) = sendnode*size(bdmtxbuf,1)
         end if
      end do
      datlen = max(library(blib)%my_nze,1)
      if (bmat%iscmplx) then
         call comms_allgather(bzmtxbuf(:,0),bmat%zmtx(1:datlen), &
              length_src=lengths(pub_my_rank_in_group), &
              lengths_dest=lengths,displs_dest=displs,comm=pub_group_comm)
      else
         call comms_allgather(bdmtxbuf(:,0),bmat%dmtx(1:datlen), &
              length_src=lengths(pub_my_rank_in_group), &
              lengths_dest=lengths,displs_dest=displs,comm=pub_group_comm)
      end if

      ! Share index of cmat to all nodes in group (no need to share data
      ! as it is about to be overwritten)
      call comms_allgather(cseginfobuf(1:3,0:pub_total_num_nodes,0), &
           library(clib)%seg_info(1:3,0:pub_total_num_nodes), &
           comm=pub_group_comm)
      lengths(0:pub_comms_group_size-1) = library(clib)%idx_lens( &
           pub_first_node_in_group:pub_last_node_in_group)
      do sendnode=0,pub_comms_group_size-1
         displs(sendnode) = sendnode*size(cidxbuf,1)
      end do
      buflen = max(library(clib)%idx_lens(pub_my_node_id),1)
      call comms_allgather(cidxbuf(:,0),library(clib)%blk_idx(1:buflen), &
           length_src=lengths(pub_my_rank_in_group), &
           lengths_dest=lengths,displs_dest=displs,comm=pub_group_comm)
      call comms_allgather(cptrbuf(:,0),library(clib)%blk_ptr(1:buflen), &
           length_src=lengths(pub_my_rank_in_group), &
           lengths_dest=lengths,displs_dest=displs,comm=pub_group_comm)

      ! Free up memory in comms_mod
      call comms_free

      ! Deallocate length and displacement buffers for allgatherv
      deallocate(displs,stat=ierr)
      call utils_dealloc_check('sparse_product','displs',ierr)
      deallocate(lengths,stat=ierr)
      call utils_dealloc_check('sparse_product','lengths',ierr)

   end subroutine internal_init_shares

    ! This subroutine ensures all the requests and ongoing send operations
    ! have completed then deallocates temporary arrays, and in the case of
    ! comms groups being used, sums the result matrix over the comms group
    ! Written by Nicholas Hine in July 2009
    subroutine internal_exit_shares

      integer :: node
      integer :: datlen

      ! Sum resulting cmat data over all nodes in group
      do node=0,pub_comms_group_size-1
         datlen = cseginfobuf(s_ptr,pub_total_num_nodes,node) - 1
         if (cmat%iscmplx) then
            call comms_reduce('SUM',cmat%zmtx(:), &
                 length=datlen,comm=pub_group_comm,root=node, &
                 z_array_src=czmtxbuf(1:datlen,node))
         else
            call comms_reduce('SUM',cmat%dmtx(:), &
                 length=datlen,comm=pub_group_comm,root=node, &
                 d_array_src=cdmtxbuf(1:datlen,node))
         end if
      end do
      call comms_free

      ! Deallocate communication buffers
      deallocate(cptrbuf,stat=ierr)
      call utils_dealloc_check('sparse_product','cptrbuf',ierr)
      deallocate(cidxbuf,stat=ierr)
      call utils_dealloc_check('sparse_product','cidxbuf',ierr)
      deallocate(cseginfobuf,stat=ierr)
      call utils_dealloc_check('sparse_product','cseginfobuf',ierr)
      deallocate(cdmtxbuf,stat=ierr)
      call utils_dealloc_check('sparse_product','cdmtxbuf',ierr)
      deallocate(czmtxbuf,stat=ierr)
      call utils_dealloc_check('sparse_product','czmtxbuf',ierr)
      deallocate(bptrbuf,stat=ierr)
      call utils_dealloc_check('sparse_product','bptrbuf',ierr)
      deallocate(bidxbuf,stat=ierr)
      call utils_dealloc_check('sparse_product','bidxbuf',ierr)
      deallocate(bseginfobuf,stat=ierr)
      call utils_dealloc_check('sparse_product','bseginfobuf',ierr)
      deallocate(bdmtxbuf,stat=ierr)
      call utils_dealloc_check('sparse_product','bdmtxbuf',ierr)
      deallocate(bzmtxbuf,stat=ierr)
      call utils_dealloc_check('sparse_product','bzmtxbuf',ierr)

    end subroutine internal_exit_shares

    subroutine internal_print_cropped_index

      !========================================================================!
      ! This subroutine prints the cropped indices of A, B and C for debugging !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   None                                                                 !
      !------------------------------------------------------------------------!
      ! Written by Nicholas Hine, Nov 2008.                                    !
      !========================================================================!

       integer :: nprint
       integer, parameter :: pp=10

       ! A
       write(pub_my_node_id+60,'(a,i3,a,i3)')'NODE ',pub_my_node_id, &
            ' from ',recvnode
       write(pub_my_node_id+60,'(a,a)') 'amat%structure=', amat%structure
       nprint = library(alib)%idx_lens(recvnode)-1
       call internal_fmt_int_array('acom%idxrecvbuf(:,1)',nprint,pp, &
            acom%idxbuf(1:nprint,1))
       nprint = library(alib)%idx_lens(recvnode)-1
       call internal_fmt_int_array('acom%idxrecvbuf(:,2)',nprint,pp, &
            acom%idxbuf(1:nprint,2))
       nprint = library(alib)%idx_lens(recvnode)-1
       call internal_fmt_int_array('acom%ptrrecvbuf(:,1)',nprint,pp, &
            acom%ptrbuf(1:nprint,1))
       nprint = library(alib)%idx_lens(recvnode)-1
       call internal_fmt_int_array('acom%ptrrecvbuf(:,2)',nprint,pp, &
            acom%ptrbuf(1:nprint,2))
       ! B
       write(pub_my_node_id+60,'(a,a)') 'bmat%structure=', bmat%structure
       nprint = library(blib)%idx_lens(pub_my_node_id)-1
       call internal_fmt_int_array('library(blib)%blk_idx',nprint,pp, &
            library(blib)%blk_idx(1:nprint))
       call internal_fmt_int_array('library(blib)%blk_idx',nprint,pp, &
            library(blib)%blk_ptr(1:nprint+1))
       ! C
       write(pub_my_node_id+60,'(a,a)') 'cmat%structure=', cmat%structure
       nprint = library(clib)%idx_lens(pub_my_node_id)-1
       call internal_fmt_int_array('library(clib)%blk_idx',nprint,pp, &
            library(clib)%blk_idx(1:nprint))
       call internal_fmt_int_array('library(clib)%blk_idx',nprint,pp, &
            library(clib)%blk_ptr(1:nprint+1))

    end subroutine internal_print_cropped_index

    subroutine internal_fmt_int_array(label,nprint,pp,array)

      !========================================================================!
      ! This subroutine prints an array of integers in a nicely-spaced format  !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      !   None                                                                 !
      !------------------------------------------------------------------------!
      ! Written by Nicholas Hine, Nov 2008.                                    !
      !========================================================================!

       ! Arguments
       character(*),intent(in) :: label
       integer,intent(in) :: nprint
       integer,intent(in) :: pp
       integer,intent(in) :: array(nprint)

       ! Locals
       integer :: nlines,nrem
       character(len=20) :: fmt, tmp1, tmp2, tmp3

       write(pub_my_node_id+60,'(a)') label
       nlines = nprint / pp
       nrem = modulo(nprint,pp)
       write(tmp1,'(i10)') nlines
       write(tmp2,'(i10)') pp
       write(tmp3,'(i10)') nrem
       write(fmt,'(7a)')'(',trim(adjustl(tmp1)),'(',trim(adjustl(tmp2)), &
            'i8/)',trim(adjustl(tmp3)),'i8)'
       write(pub_my_node_id+60,fmt) array(1:nprint)

    end subroutine internal_fmt_int_array

  end subroutine sparse_product

  ! This subroutine crops the received index of A so that only the nonzero
  ! blocks which contribute to C on this node are present, and forms
  ! the pointer request buffer for the required sections of the data of A
  ! Written by Nicholas Hine in July 2009
  subroutine sparse_crop_index(com,recvnode, &
       bbuflen,cbuflen,bidxbuf,cidxbuf,bseginfobuf,cseginfobuf)

    use comms, only: pub_comms_group_size, pub_first_node_in_group, &
         pub_total_num_nodes, pub_my_node_id
    use parallel_strategy, only: pub_first_atom_on_node, &
         pub_num_atoms_on_node
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    integer, intent(in) :: recvnode
    integer, intent(in) :: bbuflen, cbuflen
    integer, intent(in) :: bidxbuf(1:bbuflen,0:pub_comms_group_size-1)
    integer, intent(in) :: cidxbuf(1:cbuflen,0:pub_comms_group_size-1)
    integer, intent(in) :: bseginfobuf(1:3,0:pub_total_num_nodes, &
           0:pub_comms_group_size-1)
    integer, intent(in) :: cseginfobuf(1:3,0:pub_total_num_nodes, &
           0:pub_comms_group_size-1)

    integer :: jidx_crop
    integer :: aseg_start_crop
    integer :: nstride
    logical :: any_b_dense
    logical :: a_dense,a_sparse ! Shorthand flags for type of segment of A
    integer :: nrows, ncols, ifound, nfound
    integer :: ac_row_seg
    integer :: iidx, loc_iblk, kblk
    integer :: aseg_start, bseg_start, cseg_start
    integer :: bc_col_seg, bc_col_seg_node
    integer :: jidx, jblk, loc_jblk


#ifdef ITC_TRACE
    call VTBEGIN(vt_sparse_product_crop,vt_err)
#endif

    ! Reset the cropped index for A
    com%idxbuf(1:library(com%lib)%idx_lens(recvnode),2) = 0
    com%ptrbuf(1:library(com%lib)%idx_lens(recvnode),2) = 0
    com%seginfobuf(:,:,2) = 0

    ! Initialisations
    com%nreqs = 0
    nstride = -1
    com%reqdatlen = 1

    ! Check if any of the row-segments corresponding to recvnode
    ! in this group of col-segments of B is dense
    any_b_dense = any(bseginfobuf(s_type,recvnode, &
         0:pub_comms_group_size-1)==SEG_DENSE)

    ! Loop over (row) segments of index of A and C
    do ac_row_seg=0,pub_total_num_nodes-1

       ! Nothing to do here if A is blank for this segment
       if (com%seginfobuf(s_type,ac_row_seg,1)==SEG_BLANK) cycle
       ! Or if this row of C is blank for all the segments in this group
       if (all(cseginfobuf(s_type,ac_row_seg, &
            0:pub_comms_group_size-1)==SEG_BLANK)) cycle

       ! Find index start position for this segment of A
       aseg_start = com%seginfobuf(s_idx,ac_row_seg,1)

       ! Find out if this segment of A is sparse or dense and set nstride
       a_sparse = (com%seginfobuf(s_type,ac_row_seg,1)==SEG_SPARSE)
       a_dense = (com%seginfobuf(s_type,ac_row_seg,1)==SEG_DENSE)
       if (a_dense) nstride = num_elems_on_node(ac_row_seg, &
            library(com%lib)%row_blks)

       ! Mark where this segment will start in received data buffer
       com%seginfobuf(s_ptr,ac_row_seg,2) = com%reqdatlen

       ! If both A and any B segs are dense, request this whole segment of A
       if (a_dense.and.any_b_dense) then
          ! Find number of rows and columns of this segment
          ncols = num_elems_on_node(recvnode,library(com%lib)%col_blks)
          nrows = num_elems_on_node(ac_row_seg,library(com%lib)%row_blks)
          ! Record request pointer for A data on other node
          com%nreqs = com%nreqs + 1
          com%ptrreqsendbuf(1,com%nreqs) = com%seginfobuf(s_ptr,ac_row_seg,1)
          com%ptrreqsendbuf(2,com%nreqs) = ncols
          com%ptrreqsendbuf(3,com%nreqs) = nrows
          com%ptrreqsendbuf(4,com%nreqs) = nstride

          ! Add this whole segment to the cropped index of A
          ! Loop over atoms on recvnode (cols of segment of A)
          do loc_iblk=1,pub_num_atoms_on_node(recvnode)

             ! Loop over block-rows kblk of A in column iblk
             do jidx=com%idxbuf(aseg_start+loc_iblk-1,1), &
                  com%idxbuf(aseg_start+loc_iblk,1)-1

                ! Add this block to the cropped index, adjusting pointer
                ! to correct offset
                com%idxbuf(jidx,2) = com%idxbuf(jidx,1)
                com%ptrbuf(jidx,2) = com%ptrbuf(jidx,1) &
                     - com%seginfobuf(s_ptr,ac_row_seg,1) + com%reqdatlen
             end do

          end do

          ! Advance the request pointer past this whole segment
          com%reqdatlen = com%reqdatlen + nrows*ncols

       else

       ! None of the B segments which will multiply this A segment in this
       ! group are dense, so we have to go through all the B segments and
       ! mark which blocks need to be requested from the node that holds A.
       do bc_col_seg=0,pub_comms_group_size-1

          if (bseginfobuf(s_type,recvnode,bc_col_seg)==SEG_BLANK) cycle
          if (cseginfobuf(s_type,ac_row_seg,bc_col_seg)==SEG_BLANK) cycle

          bc_col_seg_node = bc_col_seg + pub_first_node_in_group

          cseg_start = cseginfobuf(s_idx,ac_row_seg,bc_col_seg)
          bseg_start = bseginfobuf(s_idx,recvnode,bc_col_seg)

          ! Loop over block columns of B and C
          do loc_iblk=1,pub_num_atoms_on_node(bc_col_seg_node)

             ! Set number of blocks found to zero
             nfound = 0

             ! Loop over block-rows kblk of C in column iblk
             do iidx=cidxbuf(cseg_start+loc_iblk-1,bc_col_seg), &
                  cidxbuf(cseg_start+loc_iblk,bc_col_seg)-1
                kblk = cidxbuf(iidx,bc_col_seg)

                ! Mark this block-row as found and note index and pointer
                found_blk(kblk) = .true.
                nfound = nfound + 1
                found_idx(nfound) = kblk

             end do  ! Loop over block-rows kblk of C in column iblk

             ! Loop over block-rows jblk of B in column iblk
             do iidx=bidxbuf(bseg_start+loc_iblk-1,bc_col_seg), &
                  bidxbuf(bseg_start+loc_iblk,bc_col_seg)-1
                jblk = bidxbuf(iidx,bc_col_seg)

                ! Find local number of this atom
                loc_jblk = jblk - pub_first_atom_on_node(recvnode) + 1

                ! Set number of rows in this block-row
                ncols = num_elems_on_atom(jblk,library(com%lib)%col_blks)

                ! Loop over block-rows kblk of A in column jblk
                do jidx=com%idxbuf(aseg_start+loc_jblk-1,1), &
                     com%idxbuf(aseg_start+loc_jblk,1)-1
                   kblk = com%idxbuf(jidx,1)

                   ! If this row occurs in C, add this block to the
                   ! cropped index of A and request it from recvnode
                   if (found_blk(kblk)) then

                      nrows = num_elems_on_atom(kblk, &
                           library(com%lib)%row_blks)
                      if (a_sparse) nstride = nrows

                      ! If this entry has not yet been recorded in the
                      ! cropped index aidxrecvbuf(:,2)
                      if (com%idxbuf(jidx,2)/=kblk) then
                         com%nreqs = com%nreqs + 1
                         ! Record request pointer for A data on other node
                         com%ptrreqsendbuf(1,com%nreqs) = com%ptrbuf(jidx,1)
                         com%ptrreqsendbuf(2,com%nreqs) = ncols
                         com%ptrreqsendbuf(3,com%nreqs) = nrows
                         com%ptrreqsendbuf(4,com%nreqs) = nstride
                         ! Add this block to the cropped index of A
                         com%idxbuf(jidx,2) = kblk
                         com%ptrbuf(jidx,2) = com%reqdatlen
                         com%reqdatlen = com%reqdatlen + nrows*ncols
                      end if

                   end if

                end do  ! Column blocks kblk of A in row jblk

             end do  ! Block-rows jblk of B in column iblk

             ! Reset found block-rows flags
             do ifound=1,nfound
                iidx = found_idx(ifound)
                found_blk(iidx) = .false.
             end do

          end do  ! Loop over block-columns iblk of B and C

       end do ! Col-Segments bc_col_seg of B,C in this group

       end if

    end do ! Row-Segments ac_row_seg of remote A cols

    ! Request pointer finishes at 1 beyond end of array, so set to end
    com%reqdatlen = com%reqdatlen - 1

    ! Shift the cropped indices of A down
    ! Find the start of the block list for first column of A
    jidx_crop = 1
    aseg_start_crop = 1

    do ac_row_seg=0,pub_total_num_nodes-1

       if (com%seginfobuf(s_type,ac_row_seg,1)==SEG_DENSE) then
          a_dense = .true.
       else
          a_dense = .false.
       end if

       ! Mark segment start in cropped segment info
       com%seginfobuf(s_idx,ac_row_seg,2) = aseg_start_crop

       ! Nothing to do here if A or C are blank for this segment
       if ((com%seginfobuf(s_type,ac_row_seg,1)==SEG_BLANK).or.&
            (all(cseginfobuf(s_type,ac_row_seg,0:pub_comms_group_size-1) &
            ==SEG_BLANK))) then
          com%seginfobuf(s_type,ac_row_seg,2) = SEG_BLANK
          cycle
       end if

       ! Find segment start positions
       aseg_start = com%seginfobuf(s_idx,ac_row_seg,1)

       ! Leave space for column start positions
       jidx_crop = jidx_crop + pub_num_atoms_on_node(recvnode) + 1

       ! Loop over block-columns jblk of A
       do loc_jblk=1,pub_num_atoms_on_node(recvnode)

          ! Mark the column start in the cropped index of A
          com%idxbuf(aseg_start_crop+loc_jblk-1,2) = jidx_crop

          ! Loop over block-rows kblk of A in column jblk
          do jidx=com%idxbuf(aseg_start+loc_jblk-1,1), &
               com%idxbuf(aseg_start+loc_jblk,1)-1
             ! Check if we recorded this block in first pass
             if (com%idxbuf(jidx,2)>0) then
                ! Add it to the cropped index of A
                com%idxbuf(jidx_crop,2) = com%idxbuf(jidx,2)
                com%ptrbuf(jidx_crop,2) = com%ptrbuf(jidx,2)
                jidx_crop = jidx_crop + 1
             end if
          end do

       end do

       ! Mark the last column end in the cropped index of A
       com%idxbuf(aseg_start_crop+pub_num_atoms_on_node(recvnode),2) &
            = jidx_crop

       ! If there are no blocks in this segment of the cropped index, mark
       ! it as blank so it is skipped
       if (aseg_start_crop == jidx_crop) then
          com%seginfobuf(s_type,ac_row_seg,2) = SEG_BLANK
       else
          ! If these segments of A and B are both dense, mark cropped segment
          ! as dense, otherwise sparse
          if (a_dense.and.any_b_dense) then
             com%seginfobuf(s_type,ac_row_seg,2) = SEG_DENSE
          else
             com%seginfobuf(s_type,ac_row_seg,2) = SEG_SPARSE
          end if
          aseg_start_crop = jidx_crop
       end if

    end do

    ! Record end of last segment
    com%seginfobuf(s_idx,pub_total_num_nodes,2) = jidx_crop
    com%seginfobuf(s_ptr,pub_total_num_nodes,2) = com%reqdatlen + 1

#ifdef ITC_TRACE
    call VTEND(vt_sparse_product_crop,vt_err)
#endif

  end subroutine sparse_crop_index

  !============================================================================!
  ! This subroutine receives the index, segment info, pointers and data of A   !
  ! from recvnode, while checking for any new requests arriving                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com  (input) : The sparse matrix communicator of the matrix (COM3 type)  !
  !   mat  (input) : The sparse matrix itself (SPAM3 type)                     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in July 2009, later made standalone.              !
  !============================================================================!

  subroutine sparse_get_crop_step_data(com,mat,recvnode, &
       bbuflen,cbuflen,bidxbuf,cidxbuf,bseginfobuf,cseginfobuf)

    use comms, only: comms_send, comms_irecv, comms_wait, comms_waitany, &
         pub_my_node_id, pub_null_handle, pub_total_num_nodes, &
         pub_comms_group_size

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    type(SPAM3), intent(in) :: mat
    integer, intent(in) :: recvnode
    integer, intent(in) :: bbuflen, cbuflen
    integer, intent(in) :: bidxbuf(1:bbuflen,0:pub_comms_group_size-1)
    integer, intent(in) :: cidxbuf(1:cbuflen,0:pub_comms_group_size-1)
    integer, intent(in) :: bseginfobuf(1:3,0:pub_total_num_nodes, &
           0:pub_comms_group_size-1)
    integer, intent(in) :: cseginfobuf(1:3,0:pub_total_num_nodes, &
           0:pub_comms_group_size-1)

    ! Local Variables
    logical :: got_index, got_reqs

    got_index = .false.
    got_reqs = .false.

    ! Local node
    if (recvnode==pub_my_node_id) then
       ! Call 'comms' routines to copy local data only
       call sparse_recv_index(com,recvnode,2,.false.,.true.)
       call sparse_recv_pointers(com,recvnode,2,.false.,.true.)
       call sparse_recv_data(com,mat,recvnode,2,.false.,.true.)
    ! Remote node
    else
       ! Loop until all required components have been received
       do

         ! Test if any data needs to be sent
         call sparse_check_send_requests(com,mat)

         ! Check if the index receives have been completed
         if ( (com%handles(com%idx_hdl)==pub_null_handle) .and. &
              (com%handles(com%info_hdl)==pub_null_handle) .and. &
              (com%handles(com%ptr_hdl)==pub_null_handle).and. &
              (.not.got_index)) then
             call comms_wait(com%handles(com%idx_hdl))
             call comms_wait(com%handles(com%info_hdl))
             call comms_wait(com%handles(com%ptr_hdl))

             ! Crop the index of A received from recvnode
             ! and form the list of required sections of data
             call sparse_crop_index(com,recvnode, &
                  bbuflen,cbuflen,bidxbuf,cidxbuf,bseginfobuf,cseginfobuf)

             ! If there are no blocks required from recvnode, send 0 to
             ! recvnode and move on to next node (all segs of A will be
             ! blank in the cropped index)
             if (com%nreqs==0) then
                call comms_send(recvnode,0,1,tag=DATA_REQ_TAG)
                exit
             end if

             ! Send number of blocks required to recvnode as request
             call comms_send(recvnode,com%nreqs,1,tag=DATA_REQ_TAG)

             ! Start asynchronous receives for the matrix data
             if (mat%iscmplx) then
                call comms_irecv(recvnode,com%zmtxrecvbuf(:,2),com%reqdatlen, &
                     tag=DATA_TAG, handle=com%handles(com%data_hdl))
             else
                call comms_irecv(recvnode,com%dmtxrecvbuf(:,2),com%reqdatlen, &
                     tag=DATA_TAG, handle=com%handles(com%data_hdl))
             end if

             ! Send the request list to recvnode
             call comms_send(recvnode,com%ptrreqsendbuf,com%nreqs*4, &
                  tag=PTR_REQ_TAG,return_handle=com%handles(com%req_send_hdl), &
                  add_to_stack=.false.)

             got_index = .true.
             cycle

         end if

         ! Check if the request send has been completed
         if (got_index .and. (.not.got_reqs) .and. &
              (com%handles(com%req_send_hdl)==pub_null_handle)) then
             call comms_wait(com%handles(com%req_send_hdl))

             got_reqs = .true.
             cycle

         end if

         ! Check if the data receive has been completed
         if (got_index .and. got_reqs .and. &
              (com%handles(com%data_hdl)==pub_null_handle)) then
            call comms_wait(com%handles(com%data_hdl))
            ! All done, so leave loop
            exit
         end if

         ! This point is only reached if there are no outstanding requests
         ! and we are still waiting for the data
         call comms_waitany(com%num_handles,com%handles)

       end do
    end if

  end subroutine sparse_get_crop_step_data


  !============================================================================!
  ! This subroutine receives the index, segment info, pointers and data of A   !
  ! from recvnode, while checking for any new requests arriving                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com  (input) : The sparse matrix communicator of the matrix (COM3 type)  !
  !   mat  (input) : The sparse matrix itself (SPAM3 type)                     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in July 2009, later made standalone.              !
  !============================================================================!

  subroutine sparse_get_step_data(com,mat,recvnode)

    use comms, only: comms_send, comms_irecv, comms_wait, comms_waitany, &
         pub_my_node_id, pub_null_handle, pub_total_num_nodes
    use parallel_strategy, only: pub_num_atoms_on_node

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    type(SPAM3), intent(in) :: mat
    integer, intent(in) :: recvnode

    ! Local Variables
    integer :: seg_start, iblk, buflen

    ! Local node
    if (recvnode==pub_my_node_id) then
       ! Call 'comms' routines to copy local data only
       call sparse_recv_index(com,recvnode,2,.true.,.true.)
       call sparse_recv_pointers(com,recvnode,2,.true.,.true.)
       call sparse_recv_data(com,mat,recvnode,2,.true.,.true.)

    ! Remote node
    else
       ! Loop until all required components have been received
       do

         ! Test if any data needs to be sent
         call sparse_check_send_requests(com,mat)

         ! Check if the index receives have been completed
         if ( (com%handles(com%idx_hdl)==pub_null_handle) .and. &
              (com%handles(com%info_hdl)==pub_null_handle) .and. &
              (com%handles(com%ptr_hdl)==pub_null_handle).and. &
              (com%handles(com%data_hdl)==pub_null_handle)) then
             call comms_wait(com%handles(com%idx_hdl))
             call comms_wait(com%handles(com%info_hdl))
             call comms_wait(com%handles(com%ptr_hdl))
             call comms_wait(com%handles(com%data_hdl))

             ! Shift to other buffer
             com%seginfobuf(:,0,2) = com%seginfobuf(:,0,1)
             com%idxbuf(:,2) = com%idxbuf(:,1)
             com%ptrbuf(:,2) = com%ptrbuf(:,1)
             if (com%iscmplx) then
                com%zmtxrecvbuf(:,2) = com%zmtxrecvbuf(:,1)
             else
                com%dmtxrecvbuf(:,2) = com%dmtxrecvbuf(:,1)
             end if

            ! All done, so leave loop
            exit
         end if

         ! This point is only reached if there are no outstanding requests
         ! and we are still waiting for the data
         call comms_waitany(com%num_handles,com%handles)

       end do
    end if

    ! Index counts from start position in full list on recvnode
    ! so shift it down by seg_start, so that the nzb's start at
    ! index position pub_num_atoms_on_node(recvnode) + 2
    seg_start = com%idxbuf(1,2)
    do iblk=1,pub_num_atoms_on_node(recvnode)+1
       com%idxbuf(iblk,2) = com%idxbuf(iblk,2) - seg_start + &
            pub_num_atoms_on_node(recvnode) + 2
    end do

    ! Pointers count from first element of segment in full array
    ! so shift them down so that they start from 1 in the received array
    buflen = library(com%lib)%idx_seg_lens(recvnode)
    do iblk=pub_num_atoms_on_node(recvnode)+2,buflen
       com%ptrbuf(iblk,2) = com%ptrbuf(iblk,2) - com%seginfobuf(s_ptr,0,2) + 1
    end do
    if (recvnode == pub_my_node_id) then
       do iblk=1,pub_num_atoms_on_node(recvnode)
          com%ptrbuf(iblk,2) = com%ptrbuf(iblk,2) - com%seginfobuf(s_ptr,0,2) + 1
       end do
    end if

  end subroutine sparse_get_step_data


  !============================================================================!
  ! This subroutine receives the index and segment info of a matrix            !
  ! from recvnode, while checking for any new requests arriving                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com  (input) : The sparse matrix communicator of the matrix (COM3 type)  !
  !   mat  (input) : The sparse matrix itself (SPAM3 type)                     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in July 2009, later made standalone.              !
  !============================================================================!

  subroutine sparse_get_step_index(com,recvnode)

    use comms, only: comms_send, comms_irecv, comms_wait, comms_waitany, &
         pub_my_node_id, pub_null_handle, pub_total_num_nodes
    use parallel_strategy, only: pub_num_atoms_on_node

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    integer, intent(in) :: recvnode

    ! Local node
    if (recvnode==pub_my_node_id) then
       ! Call 'comms' routines to copy local data only
       call sparse_recv_index(com,recvnode,2,.false.,.true.)

    ! Remote node
    else
       ! Loop until all required components have been received
       do

         ! Test if any data needs to be sent
         call sparse_check_send_requests(com)

         ! Check if the index receives have been completed
         if ( (com%handles(com%idx_hdl)==pub_null_handle) .and. &
              (com%handles(com%info_hdl)==pub_null_handle)) then
             call comms_wait(com%handles(com%idx_hdl))
             call comms_wait(com%handles(com%info_hdl))

             ! Shift to other buffer
             com%seginfobuf(:,:,2) = com%seginfobuf(:,:,1)
             com%idxbuf(:,2) = com%idxbuf(:,1)

            ! All done, so leave loop
            exit
         end if

         ! This point is only reached if there are no outstanding requests
         ! and we are still waiting for the data
         call comms_waitany(com%num_handles,com%handles)

       end do
    end if

  end subroutine sparse_get_step_index

  !============================================================================!
  ! This subroutine checks the index and data requests from all nodes and      !
  ! sends off indices, pointers and data as required. Does not exit until      !
  ! either there are no more outstanding requests, or there are outstanding    !
  ! data requests only but the data send buffer is still in use.               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com  (input) : The sparse matrix communicator of the matrix (COM3 type)  !
  !   mat  (input) : The sparse matrix itself (SPAM3 type)                     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in July 2009, later made standalone.              !
  !============================================================================!

  subroutine sparse_check_send_requests(com,mat)

    use comms, only: comms_probe, comms_recv, comms_send, &
         comms_test, comms_wait, pub_my_node_id, &
         pub_total_num_nodes, pub_comms_group_size
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    type(SPAM3), intent(in), optional :: mat

    ! Local Variables
    integer :: sendnode
    integer :: jreq,sendreqs
    integer :: sendptr,origptr
    integer :: nrows,ncols,nstride,icol
    integer :: nodestride
    logical :: can_exit
    logical :: probe_test

    ! Test whether previous data send has completed yet
    ! and ensure that ongoing mpi calls complete
    if (.not.com%send_buffer_free) then
       call comms_test(com%send_buffer_free,com%handles(com%data_send_hdl))
    else
       call comms_probe(probe_test,modulo(pub_my_node_id+1, &
            pub_total_num_nodes))
    end if

    ! Loop over nodes until no more messages are left
    sendnode = pub_my_node_id
    can_exit = .true.
    nodestride = 1
    if (com%cropped) nodestride = pub_comms_group_size
    do

       ! Check if index request has been received
       if (com%index_reqs(sendnode)==index_needed) then

          ! Make sure request has finished being received
          call comms_wait(com%handles(sendnode+1))

          ! Different response depending on whether this routine is
          ! called from sparse_product (hence cropped=true), or
          ! sparse_trace/transpose/expand (hence seg_only=true), or
          ! sparse_create (hence send only index).
          if (com%cropped) then
             ! Send the index and pointers only
             call sparse_send_index(com,sendnode,.false.)
             call sparse_send_pointers(com,sendnode,.false.)
          else if (allocated(com%dmtxrecvbuf)) then
             ! Send segment index, pointers and data
             call sparse_send_index(com,sendnode,.true.)
             call sparse_send_pointers(com,sendnode,.true.)
             call sparse_send_data(mat,sendnode,.true.)
          else
             ! Send segment index only
             call sparse_send_index(com,sendnode,.false.)
          end if

          ! Mark completion of index send to sendnode
          com%index_reqs(sendnode) = index_sent

          ! Need to restart check, since a request came in
          can_exit = .false.

       end if

       ! Check if the send buffer has become free recently
       if ((com%data_reqs(sendnode)>=0).and.(.not.com%send_buffer_free)) then
          call comms_test(com%send_buffer_free,com%handles(com%data_send_hdl))
       end if

       if ((com%data_reqs(sendnode)>=0).and.(com%send_buffer_free)) then
#ifdef ITC_TRACE
          call VTBEGIN(vt_sparse_product_data_reply,vt_err)
#endif
          ! Make sure request has finished being received
          call comms_wait(com%handles(sendnode+pub_total_num_nodes+1))

          ! No data to send, so just mark this node as done and move on
          if (com%data_reqs(sendnode)==0) then
             com%data_reqs(sendnode) = data_sent
          else
             ! Receive the pointer request buffer
             sendreqs = 4*com%data_reqs(sendnode)
             call comms_recv(sendnode,com%ptrreqrecvbuf,sendreqs,tag=PTR_REQ_TAG)
             sendreqs = com%data_reqs(sendnode)
             sendptr = 1

             ! Do not re-fill the send buffer unless the previous data send
             ! is done. The corresponding receive was posted before the
             ! send started so there is no chance of lockup.
             call comms_wait(com%handles(com%data_send_hdl))

             ! Fill a(d/z)mtxsendbuf with the requested data
             do jreq=1,sendreqs
                origptr = com%ptrreqrecvbuf(1,jreq)
                ncols = com%ptrreqrecvbuf(2,jreq)
                nrows = com%ptrreqrecvbuf(3,jreq)
                nstride = com%ptrreqrecvbuf(4,jreq)
                if (mat%iscmplx) then
                   do icol=1,ncols
                      com%zmtxsendbuf(sendptr:sendptr+nrows-1,1) = &
                           mat%zmtx(origptr:origptr+nrows-1)
                      origptr = origptr + nstride
                      sendptr = sendptr + nrows
                   end do
                else
                   do icol=1,ncols
                      com%dmtxsendbuf(sendptr:sendptr+nrows-1,1) = &
                           mat%dmtx(origptr:origptr+nrows-1)
                      origptr = origptr + nstride
                      sendptr = sendptr + nrows
                   end do
                end if
             end do

             ! Set sendptr to the length of the data
             sendptr = sendptr - 1

             ! Send the data requested by sendnode
             ! Do not add these to the send stack - the handles are
             ! managed explicitly to prevent buffer-overwrite
             if (mat%iscmplx) then
                call comms_send(sendnode,com%zmtxsendbuf(:,1),sendptr,tag=DATA_TAG, &
                     return_handle=com%handles(com%data_send_hdl), &
                     add_to_stack=.false.)
             else
                call comms_send(sendnode,com%dmtxsendbuf(:,1),sendptr,tag=DATA_TAG, &
                     return_handle=com%handles(com%data_send_hdl), &
                     add_to_stack=.false.)
             end if

             ! Mark completion of data send to sendnode
             com%data_reqs(sendnode) = data_sent

             ! Prevent send buffer being re-used until this send is finished
             com%send_buffer_free = .false.

          end if

          ! Need to restart check, since a new request may have come in
          can_exit = .false.
#ifdef ITC_TRACE
          call VTEND(vt_sparse_product_data_reply,vt_err)
#endif
       end if

       ! Advance to next node there might be a request from
       ! Only exit if nothing happened this time round
       sendnode = sendnode + nodestride
       if (sendnode>=pub_total_num_nodes) &
            sendnode = sendnode - pub_total_num_nodes
       if (sendnode==pub_my_node_id) then
          if (can_exit) then
             exit
          else
             sendnode = pub_my_node_id
             can_exit = .true.
          end if
       end if

    end do

  end subroutine sparse_check_send_requests


  !============================================================================!
  ! This subroutine initialises the requests and handles required for the      !
  ! communications operations of sparse_product                                !
  ! The communicator must already have had its %index_reqs and %data_reqs      !
  ! set up prior to entry to this routine.                                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com  (input) : The sparse matrix communicator of the matrix (COM3 type)  !
  !   mat  (input) : The sparse matrix itself (SPAM3 type)                     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in July 2009, later made standalone.              !
  !============================================================================!

  subroutine sparse_init_comms(com)

    use comms, only: comms_irecv, pub_my_node_id, pub_total_num_nodes, &
         pub_my_rank_in_group, pub_comms_group_size

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com

    ! Local Variables
    integer :: reqnode, recvnode
    integer :: rank_in_group,group_size
    integer :: first_in_req_group,last_in_req_group

    ! Set up handle indices for various messages
    if (com%cropped) then
       com%num_handles = 2*pub_total_num_nodes + 6
       com%data_send_hdl = com%num_handles - 5   ! send operation for data
       com%info_hdl = com%num_handles - 4 ! recv operation for com%seginfobuf
       com%idx_hdl = com%num_handles - 3  ! recv operation for index
       com%ptr_hdl = com%num_handles - 2  ! recv operation for pointers
       com%data_hdl = com%num_handles - 1 ! recv operation for data
       com%req_send_hdl = com%num_handles ! send operation for requests
    else
       com%num_handles = 2*pub_total_num_nodes + 4
       com%info_hdl = com%num_handles - 3 ! recv operation for com%seginfobuf
       com%idx_hdl = com%num_handles - 2  ! recv operation for index
       com%ptr_hdl = com%num_handles - 1  ! recv operation for pointers
       com%data_hdl = com%num_handles     ! recv operation for data
    end if

    ! Initialise arrays
    com%send_buffer_free = .true.
    com%index_reqs(pub_my_node_id) = index_sent
    com%data_reqs(pub_my_node_id) = data_sent

    if (com%cropped) then
       rank_in_group = pub_my_rank_in_group
       group_size = pub_comms_group_size
    else
       rank_in_group = 0
       group_size = 1
    end if

    ! Start asynchronous request receive for all those nodes which
    ! will be requesting anything from this node, by dint of having the
    ! same rank within their group as this node
    do reqnode=rank_in_group,pub_total_num_nodes-1,group_size
       ! Find range of nodes belonging to the same group as reqnodes
       first_in_req_group = reqnode - rank_in_group
       last_in_req_group = first_in_req_group + group_size - 1
       ! Node never needs to send to itself, so override for pub_my_node_id
       if (reqnode==pub_my_node_id) then
          com%index_reqs(first_in_req_group:last_in_req_group) = index_sent
          if (com%cropped) com%data_reqs(first_in_req_group:last_in_req_group) = data_sent
       end if
       ! Check if any of the group of nodes reqnode is in have nonzero
       ! index_reqs
       if (any(com%index_reqs(first_in_req_group:last_in_req_group)/=index_sent)) then
          ! Record the fact that we are waiting for an index and data request
          ! from this node, and clear the index_reqs and data_reqs for the
          ! other nodes in the same group as reqnode, since they will not
          ! request anything from this node
          com%index_reqs(reqnode) = index_not_sent
          if (com%cropped) com%data_reqs(reqnode) = data_not_sent
          do recvnode=first_in_req_group,last_in_req_group
             if (recvnode/=reqnode) then
                com%index_reqs(recvnode) = index_sent
                if (com%cropped) com%data_reqs(recvnode) = data_sent
             end if
          end do
         ! Start receive for index request
          call comms_irecv(reqnode,com%index_reqs(reqnode),1,tag=IDX_REQ_TAG, &
               handle=com%handles(reqnode+1))
          ! Start receive for data request flag
          if (com%cropped) call comms_irecv(reqnode,com%data_reqs(reqnode),1,tag=DATA_REQ_TAG, &
               handle=com%handles(reqnode+pub_total_num_nodes+1))
       end if
    end do

  end subroutine sparse_init_comms

  !============================================================================!
  ! This subroutine ensures all the requests and ongoing send operations       !
  ! have completed then deallocates temporary arrays, and in the case of       !
  ! comms groups being used, sums the result matrix over the comms group       !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com  (input) : The sparse matrix communicator of the matrix (COM3 type)  !
  !   mat  (input) : The sparse matrix itself (SPAM3 type)                     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine in July 2009, later made standalone.              !
  !============================================================================!

  subroutine sparse_exit_comms(com,mat)

    use comms, only: comms_free, comms_probe, comms_wait, comms_waitany, &
         pub_null_handle, pub_my_node_id, pub_total_num_nodes

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    type(SPAM3), intent(in), optional :: mat

    ! Ensure all requests are dealt with before deallocating
    do
       ! If there is something left to send, check incoming requests
       call sparse_check_send_requests(com,mat)

       ! Check all requests have been dealt with and that all request handles
       ! are null before exiting
       if ((all(com%index_reqs(:) == index_sent)) .and. &
            (all(com%data_reqs(:) == data_sent)) .and. &
            (all(com%handles(1:2*pub_total_num_nodes+1)==pub_null_handle))) exit

       ! Wait until something new happens: check only the requests and the
       ! data send handle (as this may be blocking the final outgoing
       ! sends from starting), as the index and pointer and data receive
       ! handles will already be null.
       call comms_waitany(com%num_handles-5,com%handles)
    end do

    ! Complete all outgoing sends before deallocating
    if (com%cropped) call comms_wait(com%handles(com%data_send_hdl))
    call comms_free

  end subroutine sparse_exit_comms


  !============================================================================!
  ! This function performs a trace on a sparse matrix or the product of two    !
  ! sparse matrices (the whole product is not calculated).                     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   amat  (input) : The sparse matrix A                                      !
  !   bmat  (input) : The (optional) sparse matrix B                           !
  !----------------------------------------------------------------------------!
  ! SPAM3 version written by Nicholas Hine, June 2009.                         !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  !============================================================================!

  real(kind=DP) function sparse_trace(amat,bmat)

    use comms, only: comms_abort, comms_free, comms_reduce, comms_send, &
         pub_on_root, pub_my_node_id, pub_total_num_nodes
    use parallel_strategy, only: pub_first_atom_on_node
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
#ifdef ITC_TRACE
    use vt
#endif

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: amat                ! The sparse matrix A
    type(SPAM3), optional, intent(in) :: bmat      ! The sparse matrix B

    ! Local variables
    integer :: alib,blib       ! Pointers to library entries
    integer :: iblk,jblk,kblk  ! Atom block loop counters
    integer :: recvnode        ! Node index of recv origin
    integer :: reqnode         ! Node index of next recv origin
    integer :: step            ! Node loop counter
    integer :: loc_iblk        ! Atom block loop counter local to this node
    integer :: loc_jblk        ! Atom block loop counter local to this node
    integer :: iidx,jidx       ! Indices
    integer :: ielems,jelems   ! Numbers of element rows/cols in this block
    integer :: ielem           ! Loop counter over elems on atom
    integer :: jelem           ! Loop counter over elems on atom
    integer :: aptr,bptr       ! Pointers to blocks in matrices A and B
    integer :: seg_start       ! Segment type for this segment
    integer :: seg_type        ! Start position of this segment in the index
    integer :: arow_blks       ! Identifier for A row, B col row blocking scheme
    integer :: brow_blks       ! Identifier for B row, A col row blocking scheme
    real(kind=DP) :: local_tr  ! Local trace on this node
    type(COM3) :: acom

#ifdef ITC_TRACE
    call VTBEGIN(vt_sparse_trace,vt_err)
#endif

    call timer_clock('sparse_trace',1)

    ! qoh: Initialise variables
    blib = -1
    brow_blks = -1

    ! Extract library entries
    alib = amat%lib
    if (present(bmat)) blib = bmat%lib

    ! Check consistency of arguments
    if (.not.present(bmat)) then
       if (library(alib)%row_blks/=library(alib)%col_blks) then
          if (pub_on_root) write(stdout,'(a)') 'Error in sparse_trace: &
               &matrix is not square'
          call comms_abort
       end if
       arow_blks = library(alib)%row_blks
    else
       if (amat%iscmplx .neqv. bmat%iscmplx) then
          if (pub_on_root) write(stdout,'(a)') &
               'Error in sparse_trace: incompatible argument types.'
          call comms_abort
       end if
       if (library(alib)%row_blks/=library(blib)%col_blks) then
          if (pub_on_root) write(stdout,'(a)') 'Error in sparse_trace: &
               &incompatible matrix blocking schemes'
          call comms_abort
       end if
       if (library(alib)%col_blks/=library(blib)%row_blks) then
          if (pub_on_root) write(stdout,'(a)') 'Error in sparse_trace: &
               &incompatible matrix blocking schemes'
          call comms_abort
       end if
       arow_blks = library(alib)%row_blks
       brow_blks = library(blib)%row_blks
    end if

    ! Zero trace accumulator for this node
    local_tr = 0.0_DP
    ielems = -1

    ! If only one argument is given, calculate the trace of matrix A
    if (.not. present(bmat)) then

       seg_type = library(alib)%seg_info(s_type,pub_my_node_id)

       if (seg_type/=SEG_BLANK) then

          ! Loop over atom-blocks iblk on this node
          loc_iblk = 0
          seg_start = library(alib)%seg_info(s_idx,pub_my_node_id)
          do iblk=my_first_blk,my_last_blk
             loc_iblk = loc_iblk + 1

             ! Get number of elems in this block
             if (seg_type==SEG_DENSE) then
                ielems = num_elems_on_node(pub_my_node_id,arow_blks)
             else if (seg_type==SEG_SPARSE) then
                ielems = num_elems_on_atom(iblk,arow_blks)
             end if

             ! Find pointer to start of diagonal block in this column
             aptr = library(alib)%blk_ptr(seg_start+loc_iblk-1)

             ! Accumulate trace from this diagonal block
             if (amat%iscmplx) then
                do ielem=1,num_elems_on_atom(iblk,arow_blks)
                   local_tr = local_tr + real(amat%zmtx(aptr),kind=DP)
                   aptr = aptr + ielems + 1
                end do
             else
                do ielem=1,num_elems_on_atom(iblk,arow_blks)
                   local_tr = local_tr + amat%dmtx(aptr)
                   aptr = aptr + ielems + 1
                end do
             end if

          end do  ! Loop over atom-blocks iblk on this node

       end if

    else  ! Calculate the trace of the product of A and B

       ! Allocate communication buffers
       call sparse_com_allocate(acom,amat,2, &
            alloc_mtx=.true.,cropped=.false.,seg=.true.)

       ! Set up which other nodes will be requesting data
       do recvnode=0,pub_total_num_nodes-1
          if (library(alib)%seg_info(s_type,recvnode)/=SEG_BLANK) then
             acom%index_reqs(recvnode) = index_not_sent
          end if
       end do

       ! Initialise amat comms
       call sparse_init_comms(acom)

       ! Loop over processors
       do step=0,pub_total_num_nodes-1

          ! Receive data for this step
          recvnode = modulo(pub_my_node_id-step+pub_total_num_nodes, &
               pub_total_num_nodes)
          if (library(alib)%idx_seg_lens(recvnode)>0) &
               call sparse_get_step_data(acom,amat,recvnode)

          ! Request index and data for next step if required
          reqnode = modulo(pub_my_node_id-step-1+pub_total_num_nodes, &
               pub_total_num_nodes)
          if ((library(alib)%idx_seg_lens(reqnode)/=0).and. &
               (reqnode/=pub_my_node_id)) then
             ! Send request for index, pointers and data
             call comms_send(reqnode,index_needed,1,tag=IDX_REQ_TAG)
             ! Start asynchronous receives for index and pointers
             call sparse_recv_index(acom,reqnode,1,.true.,async=.true.)
             call sparse_recv_pointers(acom,reqnode,1,.true.,async=.true.)
             call sparse_recv_data(acom,amat,reqnode,1,.true.,async=.true.)
          end if
          if (library(alib)%idx_seg_lens(recvnode)==0) cycle

          ! Find type and start position of the corresponding segment in B
          seg_type = library(blib)%seg_info(s_type,recvnode)
          seg_start = library(blib)%seg_info(s_idx,recvnode)

          ! Override for dense-dense trace (skips indexing)
          if ((acom%seginfobuf(s_type,0,2) == SEG_DENSE) .and. &
               (seg_type == SEG_DENSE)) then

             ! Find range of rows on A and cols on B
             ielems = num_elems_on_node(pub_my_node_id,arow_blks)
             ! Find range of cols on A and rows on B
             jelems = num_elems_on_node(recvnode,brow_blks)

             do jelem=0,jelems-1

                ! Find first elem in src needed for dest
                aptr = jelem*ielems + 1
                bptr = library(blib)%seg_info(s_ptr,recvnode) + jelem

                if (amat%iscmplx) then
                   do ielem=0,ielems-1
                      local_tr = local_tr + real(acom%zmtxrecvbuf(aptr,2) * &
                           bmat%zmtx(bptr),kind=DP)
                      aptr = aptr + 1
                      bptr = bptr + jelems
                   end do
                else
                   do ielem=0,ielems-1
                      local_tr = local_tr + acom%dmtxrecvbuf(aptr,2) * bmat%dmtx(bptr)
                      aptr = aptr + 1
                      bptr = bptr + jelems
                   end do
                endif

             end do

          ! one or other segment is sparse, so go through block-by-block
          else if (.not.((seg_type == SEG_BLANK).or. &
               (acom%seginfobuf(s_type,0,2) == SEG_BLANK))) then

             ! Calculate contribution to trace
             ! Loop over block-columns iblk of B on this node
             loc_iblk = 0
             do iblk=my_first_blk,my_last_blk
                loc_iblk = loc_iblk + 1

                ! Number of elems in this block
                if (acom%seginfobuf(s_type,0,2) == SEG_DENSE) then
                   ielems = num_elems_on_node(pub_my_node_id,arow_blks)
                else
                   ielems = num_elems_on_atom(iblk,arow_blks)
                end if

                ! Loop over block-rows jblk in block-column iblk of B
                do iidx=library(blib)%blk_idx(seg_start+loc_iblk-1), &
                     library(blib)%blk_idx(seg_start+loc_iblk)-1
                   jblk = library(blib)%blk_idx(iidx)

                   ! Find local number of this atom
                   loc_jblk = jblk - pub_first_atom_on_node(recvnode) + 1

                   ! Copy pointer to this block in B
                   bptr = library(blib)%blk_ptr(iidx)
                   if (seg_type == SEG_DENSE) then
                      jelems = num_elems_on_node(recvnode,brow_blks)
                   else
                      jelems = num_elems_on_atom(jblk,brow_blks)
                   end if

                   ! Loop over block-rows kblk in block-column jblk of A
                   do jidx=acom%idxbuf(loc_jblk,2),acom%idxbuf(loc_jblk+1,2)-1
                      kblk = acom%idxbuf(jidx,2)

                      ! If this is a diagonal block, it will contribute
                      if (kblk == iblk) then

                         ! Copy pointer to this block in A
                         aptr = acom%ptrbuf(jidx,2)

                         ! Perform block trace
                         if (amat%iscmplx) then
                            do ielem=0,num_elems_on_atom(iblk,arow_blks)-1
                               do jelem=0,num_elems_on_atom(jblk,brow_blks)-1
                                  local_tr = local_tr + real(acom%zmtxrecvbuf( &
                                       aptr+ielem+jelem*ielems,2) * &
                                       bmat%zmtx(bptr+jelem+ielem*jelems), &
                                       kind=DP)
                               end do
                            end do
                         else
                            do ielem=0,num_elems_on_atom(iblk,arow_blks)-1
                               do jelem=0,num_elems_on_atom(jblk,brow_blks)-1
                                  local_tr = local_tr + acom%dmtxrecvbuf( &
                                       aptr+ielem+jelem*ielems,2) * &
                                       bmat%dmtx(bptr+jelem+ielem*jelems)
                               end do
                            end do
                         end if
                      end if  ! Diagonal block

                   end do  ! Block-rows kblk in block-column jblk of A on node

                end do  ! Block-rows jblk in block-column iblk of B

             end do  ! Loop over atom-blocks iblk on this node

          end if

          ! Check incoming requests after each non-skipped segment
          call sparse_check_send_requests(acom,amat)

       end do  ! Loop over processors

       ! Wait for final sends and deallocate communication buffers
       call sparse_exit_comms(acom,amat)
       call sparse_com_deallocate(acom,dealloc_mtx=.true.,cropped=.false.)

    end if

    ! Sum up contributions over all nodes
    call comms_reduce('SUM',local_tr)
    sparse_trace = local_tr

    call timer_clock('sparse_trace',2)
#ifdef ITC_TRACE
    call VTEND(vt_sparse_trace,vt_err)
#endif

  end function sparse_trace


  !============================================================================!
  ! This subroutine forms the transpose of a sparse matrix.                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest  (output) : The resulting transpose                                 !
  !   src   (input)  : The matrix to be transposed                             !
  !----------------------------------------------------------------------------!
  ! SPAM3 version written by Nicholas Hine, June 2009.                         !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  !============================================================================!

  subroutine sparse_transpose(dest,src)

    use comms, only: comms_abort, comms_barrier, comms_free, comms_send, &
         pub_my_node_id, pub_on_root, pub_total_num_nodes
    use parallel_strategy, only: pub_first_atom_on_node
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: dest    ! The resulting transpose matrix
    type(SPAM3), intent(in)    :: src     ! The matrix to be transposed

    ! Local variables
    integer :: destlib         ! Library entry for dest
    integer :: srclib          ! Library entry for src
    integer :: drow_scol_blks  ! Identifier for D row, S col blocking scheme
    integer :: srow_dcol_blks  ! Identifier for S row, D col blocking scheme
    integer :: reqnode         ! Node index of send destination
    integer :: recvnode        ! Node index of recv origin
    integer :: step            ! Node loop counter
    integer :: iblk,jblk,kblk  ! Block loop counters
    integer :: loc_iblk        ! Loop counter for blocks on node
    integer :: loc_jblk        ! Loop counter for blocks on node
    integer :: iidx,jidx       ! Indices
    integer :: ielems,jelems   ! Numbers of elements in blocks/segments
    integer :: ielem,jelem     ! Element counter in block/segment
    integer :: seg_start       ! Start position of this segment in the index
    integer :: srcptr          ! Pointer to data in src
    integer :: destptr         ! Pointer to start of data in dest
    logical :: src_dense       ! Whether this source matrix segment is dense
    logical :: src_sparse      ! Whether this source matrix segment is sparse
    logical :: dest_dense      ! Whether this dest matrix segment is dense
    logical :: dest_sparse     ! Whether this dest matrix segment is sparse
    type(COM3) :: srccom

    ! Get library entries for matrices
    destlib = dest%lib
    srclib = src%lib
    drow_scol_blks = library(srclib)%col_blks
    srow_dcol_blks = library(srclib)%row_blks
    ielems = 0; jelems = 0

    ! Check arguments
    if (drow_scol_blks /= library(destlib)%row_blks) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_transpose: &
            &source cols and dest rows have different blocking schemes.'
       call comms_abort
    end if
    if (srow_dcol_blks /= library(destlib)%col_blks) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_transpose: &
            &source rows and dest cols have different blocking schemes.'
       call comms_abort
    end if
    if (dest%iscmplx .neqv. src%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_transpose: &
            &mixed real/complex argument types.'
       call comms_abort
    end if

    ! Allocate communication buffers
    call sparse_com_allocate(srccom,src,2, &
         alloc_mtx=.true.,cropped=.false.,seg=.true.)

    ! Set up which other nodes will be requesting data
    do recvnode=0,pub_total_num_nodes-1
       if (library(srclib)%seg_info(s_type,recvnode)/=SEG_BLANK) then
          srccom%index_reqs(recvnode) = index_not_sent
       end if
    end do

    ! Initialise amat comms
    call sparse_init_comms(srccom)

    ! Zero dest matrix
    if (dest%iscmplx) then
       dest%zmtx = (0.0_DP,0.0_DP)
    else
       dest%dmtx = 0.0_DP
    end if

    ! Loop over nodes
    do step=0,pub_total_num_nodes-1

       ! Receive data for this step
       recvnode = modulo(pub_my_node_id-step+pub_total_num_nodes, &
            pub_total_num_nodes)
       if (library(srclib)%idx_seg_lens(recvnode)>0) &
            call sparse_get_step_data(srccom,src,recvnode)

       ! Request index and data for next step if required
       reqnode = modulo(pub_my_node_id-step-1+pub_total_num_nodes, &
            pub_total_num_nodes)
       if ((library(srclib)%idx_seg_lens(reqnode)/=0).and. &
            (reqnode/=pub_my_node_id)) then
          ! Send request for index, pointers and data
          call comms_send(reqnode,index_needed,1,tag=IDX_REQ_TAG)
          ! Start asynchronous receives for index and pointers
          call sparse_recv_index(srccom,reqnode,1,.true.,async=.true.)
          call sparse_recv_pointers(srccom,reqnode,1,.true.,async=.true.)
          call sparse_recv_data(srccom,src,reqnode,1,.true.,async=.true.)
       end if
       if (library(srclib)%idx_seg_lens(recvnode)==0) cycle

       dest_dense = (library(destlib)%seg_info(s_type,recvnode)==SEG_DENSE)
       dest_sparse = (library(destlib)%seg_info(s_type,recvnode)==SEG_SPARSE)
       src_dense = (srccom%seginfobuf(s_type,0,2)==SEG_DENSE)
       src_sparse = (srccom%seginfobuf(s_type,0,2)==SEG_SPARSE)

       ! Override for dense-dense transpose (skips indexing)
       if (src_dense.and.dest_dense) then

          ! Find range of cols on src and rows on dest
          jelems = num_elems_on_node(recvnode,drow_scol_blks)
          ! Find range of rows on src and cols on dest
          ielems = num_elems_on_node(pub_my_node_id,srow_dcol_blks)

          do jelem=0,jelems-1

             ! Find first elem in src needed for dest
             srcptr = jelem*ielems + 1
             destptr = library(destlib)%seg_info(s_ptr,recvnode) + jelem

             ! Copy columns of node to rows of my_node
             if (src%iscmplx) then
                do ielem=0,ielems-1
                   dest%zmtx(destptr) = conjg(srccom%zmtxrecvbuf(srcptr,2))
                   srcptr = srcptr + 1
                   destptr = destptr + jelems
                end do
             else
                do ielem=0,ielems-1
                   dest%dmtx(destptr) = srccom%dmtxrecvbuf(srcptr,2)
                   srcptr = srcptr + 1
                   destptr = destptr + jelems
                end do
             end if

          end do

          call sparse_enforce_sparsity(dest)

       ! one or other segment is sparse, so go through block-by-block
       else

          if (dest_dense) jelems = &
               num_elems_on_node(recvnode,drow_scol_blks)
          if (src_dense) ielems = &
               num_elems_on_node(pub_my_node_id,srow_dcol_blks)

          ! Loop over block-columns iblk of dest on this node
          loc_iblk = 0
          seg_start = library(destlib)%seg_info(s_idx,recvnode)
          do iblk=my_first_blk,my_last_blk
             loc_iblk = loc_iblk + 1

             ! Number of elems in this block
             if (src_sparse) ielems = num_elems_on_atom(iblk,srow_dcol_blks)

             ! Loop over block-rows jblk in block-column iblk of dest
             do iidx=library(destlib)%blk_idx(seg_start+loc_iblk-1), &
                     library(destlib)%blk_idx(seg_start+loc_iblk)-1
                jblk = library(destlib)%blk_idx(iidx)

                ! Find local number of this atom on recvnode
                loc_jblk = jblk - pub_first_atom_on_node(recvnode) + 1

                ! Find pointer to this block in dest, and stride length
                destptr = library(destlib)%blk_ptr(iidx)
                if (dest_sparse) jelems = num_elems_on_atom(jblk,drow_scol_blks)

                ! Now find this block in src
                do jidx=srccom%idxbuf(loc_jblk,2),srccom%idxbuf(loc_jblk+1,2)-1
                   kblk = srccom%idxbuf(jidx,2)
                   if (kblk == iblk) then

                      ! Find pointer to first element in block and stride length
                      srcptr = srccom%ptrbuf(jidx,2)

                      ! Copy elements while transposing block
                      if (src%iscmplx) then
                         do jelem=0,num_elems_on_atom(jblk,drow_scol_blks)-1
                            do ielem=0,num_elems_on_atom(iblk,srow_dcol_blks)-1
                               dest%zmtx(destptr+ielem*jelems+jelem) = &
                                    conjg(srccom%zmtxrecvbuf(srcptr+jelem*ielems+ielem,2))
                            end do
                         end do
                      else
                         do jelem=0,num_elems_on_atom(jblk,drow_scol_blks)-1
                            do ielem=0,num_elems_on_atom(iblk,srow_dcol_blks)-1
                               dest%dmtx(destptr+ielem*jelems+jelem) = &
                                    srccom%dmtxrecvbuf(srcptr+jelem*ielems+ielem,2)
                            end do
                         end do
                      end if

                      exit
                   end if

                end do  ! Loop finding block in src

             end do  ! Loop over block-rows in col iblk of dest

          end do  ! Loop over block-colks iblk of dest

       end if

       ! Check incoming requests after each non-skipped segment
       call sparse_check_send_requests(srccom,src)

    end do  ! step

    ! Finish comms and deallocate communication buffers
    call sparse_exit_comms(srccom,src)
    call sparse_com_deallocate(srccom,dealloc_mtx=.true.,cropped=.false.)

    call comms_free
    call comms_barrier

  end subroutine sparse_transpose

  !============================================================================!
  ! This subroutine expands a set of elements related by symmetry to fill a    !
  ! square sparse matrix. Choice of elements depends on variable 'pattern'     !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat  (inout)  : The matrix to be expanded                                !
  !   pattern (in)  : PATTERN_LOWER means lower triangle elements filled in    !
  !                   from upper triangle ones                                 !
  !                   PATTERN_ALTERNATE means alternating elements along a row !
  !                   or down a column come from upper and lower triangles     !
  !   sym  (in)     : Whether the matrix is symmetric (sym=.false. means       !
  !                   signs of elements are flipped on transposition).         !
  !----------------------------------------------------------------------------!
  ! SPAM3 version written by Nicholas Hine, May 2009 based on SPAM2 version.   !
  !----------------------------------------------------------------------------!
  ! SPAM2 version written by Peter Haynes, March 2004.                         !
  ! Revised for distributed data, June 2006.                                   !
  ! Dense matrix version by Nicholas Hine, Dec 2007.                           !
  ! Facility for anti-symmetric (-Hermitian) matrices added PDH May 2008.      !
  !============================================================================!

  subroutine sparse_expand(mat,pattern,sym)

    use comms, only: comms_abort, comms_barrier, comms_free, comms_send, &
         pub_my_node_id, pub_on_root, pub_total_num_nodes
    use parallel_strategy, only: pub_first_atom_on_node
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat    ! The matrix
    integer,intent(in) :: pattern
    logical,optional,intent(in) :: sym

    ! Local variables
    integer :: ilib            ! Library entry for matrix
    integer :: col_blks        ! Identifier for column blocking scheme
    integer :: row_blks        ! Identifier for row blocking scheme
    integer :: reqnode         ! Node index of next recv origin
    integer :: recvnode        ! Node index of recv origin
    integer :: step            ! Node loop counter
    integer :: iblk,jblk,kblk  ! Block loop counters
    integer :: loc_iblk        ! Loop counter for blocks on node
    integer :: loc_jblk        ! Loop counter for blocks on node
    integer :: iidx,jidx       ! Indices
    integer :: ielems,jelems   ! Numbers of elements in blocks/segments
    integer :: ielem,jelem     ! Element counter in block/segment
    integer :: ielem0,jelem0   ! Element row/col counter to top left of block
    integer :: seg_type        ! Segment index
    integer :: seg_start       ! Start position of this segment in the index
    integer :: srcptr          ! Pointer to source data
    integer :: destptr         ! Pointer to dest data
    real(kind=DP) :: sgn       ! Sign to multiply by while transposing block
    type(COM3) :: matcom

    ! Deal with optional argument
    sgn = 1.0_DP
    if (present(sym)) then
       if (.not. sym) sgn = -1.0_DP
    end if

    ! Get library entries for matrices
    ilib = mat%lib
    row_blks = library(ilib)%row_blks
    col_blks = library(ilib)%col_blks

    ! Check this matrix can be expanded
    if (row_blks /= col_blks) then
       if (pub_on_root) write(stdout,'(a)') 'Error in sparse_expand: &
            &Non-square matrices can not be expanded'
       call comms_abort
    end if

    ! Allocate communication buffers
    call sparse_com_allocate(matcom,mat,2, &
         alloc_mtx=.true.,cropped=.false.,seg=.true.)

    ! Set up which other nodes will be requesting data
    do recvnode=0,pub_total_num_nodes-1
       if (library(ilib)%seg_info(s_type,recvnode)/=SEG_BLANK) then
          matcom%index_reqs(recvnode) = index_not_sent
       end if
    end do

    ! Initialise amat comms
    call sparse_init_comms(matcom)

    ! Loop over nodes
    do step=0,pub_total_num_nodes-1

       ! Receive data for this step
       recvnode = modulo(pub_my_node_id-step+pub_total_num_nodes, &
            pub_total_num_nodes)
       if (library(ilib)%idx_seg_lens(recvnode)>0) &
            call sparse_get_step_data(matcom,mat,recvnode)

       ! Request index and data for next step if required
       reqnode = modulo(pub_my_node_id-step-1+pub_total_num_nodes, &
            pub_total_num_nodes)
       if ((library(ilib)%idx_seg_lens(reqnode)/=0).and. &
            (reqnode/=pub_my_node_id)) then

          ! Send request for index, pointers and data
          call comms_send(reqnode,index_needed,1,tag=IDX_REQ_TAG)

          ! Start asynchronous receives for index and pointers
          call sparse_recv_index(matcom,reqnode,1,.true.,async=.true.)
          call sparse_recv_pointers(matcom,reqnode,1,.true.,async=.true.)
          call sparse_recv_data(matcom,mat,reqnode,1,.true.,async=.true.)
       end if
       if (library(ilib)%idx_seg_lens(recvnode)==0) cycle

       seg_type = library(ilib)%seg_info(s_type,recvnode)

       ! Override for dense-dense expand (skips indexing)
       if ((matcom%seginfobuf(s_type,0,2) == SEG_DENSE) .and. &
            (seg_type == SEG_DENSE)) then

          ! Find range of cols in this segment
          jelems = num_elems_on_node(recvnode,col_blks)
          jelem0 = first_elem_on_node(recvnode,col_blks)
          ! Find range of rows in this segment
          ielems = num_elems_on_node(pub_my_node_id,row_blks)
          ielem0 = first_elem_on_node(pub_my_node_id,row_blks)

          do jelem=0,jelems-1

             ! Find first elem in src needed for dest
             srcptr = jelem*ielems + 1
             destptr = library(ilib)%seg_info(s_ptr,recvnode) + jelem

             ! Copy columns of node to rows of my_node
             if (mat%iscmplx) then
                do ielem=0,ielems-1
                   if (pattern_test()) mat%zmtx(destptr) = &
                        sgn*matcom%zmtxrecvbuf(srcptr,2)
                   srcptr = srcptr + 1
                   destptr = destptr + jelems
                end do
             else
                do ielem=0,ielems-1
                   if (pattern_test()) mat%dmtx(destptr) = &
                        sgn*matcom%dmtxrecvbuf(srcptr,2)
                   srcptr = srcptr + 1
                   destptr = destptr + jelems
                end do
             end if

          end do

       ! One or other segment is sparse, so go through block-by-block
       else if (.not.((seg_type == SEG_BLANK).or. &
            (matcom%seginfobuf(s_type,0,2) == SEG_BLANK))) then

          ! Loop over block-columns iblk of dest on this node
          loc_iblk = 0
          seg_start = library(ilib)%seg_info(s_idx,recvnode)
          do iblk=my_first_blk,my_last_blk
             loc_iblk = loc_iblk + 1

             ! Number of elems in this block/segment, and origin column
             if (matcom%seginfobuf(s_type,0,2) == SEG_DENSE) then
                ielems = num_elems_on_node(pub_my_node_id,row_blks)
             else
                ielems = num_elems_on_atom(iblk,row_blks)
             end if
             ielem0 = first_elem_on_atom(iblk,col_blks)

             ! Loop over block-rows jblk in block-column iblk of dest
             do iidx=library(ilib)%blk_idx(seg_start+loc_iblk-1), &
                     library(ilib)%blk_idx(seg_start+loc_iblk)-1
                jblk = library(ilib)%blk_idx(iidx)
                jelem0 = first_elem_on_atom(jblk,row_blks)

                ! Find local number of this atom on recvnode
                loc_jblk = jblk - pub_first_atom_on_node(recvnode) + 1

                ! Find pointer to this block in dest, and stride length
                destptr = library(ilib)%blk_ptr(iidx)
                if (seg_type == SEG_DENSE) then
                   jelems = num_elems_on_node(recvnode,col_blks)
                else
                   jelems = num_elems_on_atom(jblk,col_blks)
                end if

                ! Now find this block in src
                do jidx=matcom%idxbuf(loc_jblk,2),matcom%idxbuf(loc_jblk+1,2)-1
                   kblk = matcom%idxbuf(jidx,2)
                   if (kblk == iblk) then

                      ! Find pointer to first element in block and stride length
                      srcptr = matcom%ptrbuf(jidx,2)

                      ! Copy elements while expanding block
                      if (mat%iscmplx) then
                         do jelem=0,num_elems_on_atom(jblk,col_blks)-1
                            do ielem=0,num_elems_on_atom(iblk,row_blks)-1
                               if (pattern_test()) mat%zmtx(destptr+ &
                                    ielem*jelems+jelem) = sgn*conjg( &
                                    matcom%zmtxrecvbuf(srcptr+jelem*ielems+ielem,2))
                            end do
                         end do
                      else
                         do jelem=0,num_elems_on_atom(jblk,col_blks)-1
                            do ielem=0,num_elems_on_atom(iblk,row_blks)-1
                               if (pattern_test()) mat%dmtx(destptr+ &
                                    ielem*jelems+jelem) = sgn*matcom%dmtxrecvbuf( &
                                    srcptr+jelem*ielems+ielem,2)
                            end do
                         end do
                      end if

                      exit
                   end if

                end do  ! Loop finding block in src

             end do  ! Loop over block-rows in col iblk of dest

          end do  ! Loop over block-colks iblk of dest

       end if

       ! Check incoming requests after each non-skipped segment
       call sparse_check_send_requests(matcom,mat)

    end do  ! step

    ! Finish comms and deallocate communication buffers
    call sparse_exit_comms(matcom,mat)
    call sparse_com_deallocate(matcom,dealloc_mtx=.true.,cropped=.false.)

    call comms_free
    call comms_barrier

contains

    logical function pattern_test()

       ! Locals
       integer :: row,col   ! row,col of src (also col,row of dest)

       pattern_test = .false.

       row = ielem0 + ielem
       col = jelem0 + jelem

       ! Never include diagonal elements
       if (row == col) return

       ! Check for i>j for lower triangle expansion
       if (pattern == PATTERN_LOWER) then
          if (row > col) pattern_test = .true.
       end if

       if (pattern == PATTERN_ALTERNATE) then
          if (((col<row).and.(mod(col+row,2)==1)).or. &
              ((col>row).and.(mod(col+row,2)==0))) pattern_test = .true.
       end if

       return

    end function pattern_test

  end subroutine sparse_expand

  !============================================================================!
  ! This subroutine obtains the maximum eigenvalue of a given sparse matrix by !
  ! an iterative conjugate gradients procedure.                                !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat (input)    : The matrix whose eigenvalue is desired                  !
  !   met (input)    : The metric to use                                       !
  !   eval (output)  : The eigenvalue estimate                                 !
  !   tol (input)    : Optional tolerance for eigenvalue estimate              !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Revised for distributed data, June 2006.                                   !
  !============================================================================!

  subroutine sparse_extremal_eigenvalue(mat,met,eval,tol)

    use comms, only: comms_abort, comms_bcast, comms_reduce, pub_on_root, &
         pub_my_node_id, pub_root_node_id, pub_total_num_nodes
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat     ! The matrix whose eigenvalue is desired
    type(SPAM3), intent(in) :: met     ! The metric to be used
    real(kind=DP), intent(out) :: eval ! The eigenvalue
    real(kind=DP), optional, intent(in) :: tol  ! Optional tolerance

    ! Local variables
    integer, parameter :: maxiter = 100  ! Maximum number of iterations
    integer, parameter :: cgreset = 5    ! Conjugate gradients reset frequency
    integer :: n                         ! Vector length
    integer :: iter                      ! Iteration counter
    integer :: conv                      ! Converged iteration counter
    real(kind=DP) :: local_tol           ! Local copy of tolerance
    real(kind=DP) :: norm                ! Normalisation
    real(kind=DP) :: normfac             ! Renormalisation factor
    real(kind=DP) :: gdotg               ! Norm of gradient
    real(kind=DP) :: gamma               ! Conjugate gradients gamma
    real(kind=DP) :: xdotd,ddotd,ddoth,ddoty    ! Dot products for line search
    real(kind=DP) :: a,b,c,disc,q,step,oldeval  ! Variables for line search
    real(kind=DP), allocatable :: xvec(:)  ! Trial eigenvector
    real(kind=DP), allocatable :: yvec(:)  ! Product of matrix and search dirn
    real(kind=DP), allocatable :: grad(:)  ! Gradient
    real(kind=DP), allocatable :: dirn(:)  ! Search direction
    real(kind=DP), allocatable :: matx(:)  ! Product of matrix and trial vector
    real(kind=DP), allocatable :: metx(:)  ! Product of metric and trial vector

    ! External function
    real(kind=DP) :: ddot


    ! This routine will only function for real square matrices
    if (mat%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_extremal_eigenvalue: real matrices only.'
       call comms_abort
    end if
    if (library(mat%lib)%row_blks/=library(mat%lib)%col_blks) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_extremal_eigenvalue: square matrices only.'
       call comms_abort
    end if

    ! Deal with optional argument
    if (present(tol)) then
       local_tol = tol
    else
       local_tol = 0.01_DP
    end if

    ! Find vector length
    n = library(mat%lib)%nrows

    ! Allocate workspace
    call internal_alloc_work

    ! Initialise eigenvector to random guess on one node and broadcast
    if (pub_on_root) call random_number(xvec)
    call comms_bcast(pub_root_node_id,xvec)

    ! Normalise
    call internal_matvec(metx,met,xvec)
    norm = abs(ddot(n,xvec,1,metx,1))
    normfac = 1.0_DP / sqrt(norm)
    xvec = xvec * normfac
    metx = metx * normfac

    ! Initialise CG
    gdotg = 1.0_DP
    gamma = 0.0_DP
    dirn = 0.0_DP

    ! Obtain initial estimate of eigenvalue and matx := mat . x
    call internal_matvec(matx,mat,xvec)
    eval = ddot(n,metx,1,matx,1)

    ! ndmh: suggested by Victor Milman to avoid possible WIN32 parallel desync
    call comms_bcast(pub_root_node_id,eval)

    ! Start of conjugate gradients loop
    conv = 0
    do iter=1,maxiter

       ! Obtain gradient grad
       grad = matx - eval * xvec

       ! Obtain conjugate direction dirn
       gamma = 1.0_DP / gdotg
       call internal_matvec(yvec,met,grad)
       gdotg = ddot(n,grad,1,yvec,1)
       if (mod(iter,cgreset) == 0) then
          gamma = 0.0_DP
       else
          gamma = gamma * gdotg
       end if
       dirn = -grad + gamma * dirn

       ! Orthogonalise to xvec
       xdotd = ddot(n,metx,1,dirn,1)
       dirn = dirn - xdotd * xvec

       ! Obtain yvec := mat . dirn
       call internal_matvec(yvec,mat,dirn)

       ! Obtain metx := met .dirn
       call internal_matvec(metx,met,dirn)

       ! Obtain line search coefficients
       ddotd = ddot(n,metx,1,dirn,1)
       ddoth = ddot(n,metx,1,matx,1)
       ddoty = ddot(n,metx,1,yvec,1)

       a = ddotd * ddoth
       call comms_bcast(pub_root_node_id,a)
       if (abs(a) < epsilon(1.0_DP)) exit

       b = eval * ddotd - ddoty
       call comms_bcast(pub_root_node_id,b)
       if (abs(b) < epsilon(1.0_DP)) exit

       c = -ddoth
       disc = b*b - 4.0_DP*a*c
       call comms_bcast(pub_root_node_id,disc)
       if (disc < 0.0_DP) exit

       q = -0.5_DP * (b + sign(sqrt(disc),b))
       if (b < 0.0_DP) then
          step = q / a
       else
          step = c / q
       end if

       ! Update vector xvec and matx
       xvec = xvec + step * dirn
       matx = matx + step * yvec

       ! Renormalise
       call internal_matvec(metx,met,xvec)
       norm = abs(ddot(n,xvec,1,metx,1))
       normfac = 1.0_DP / sqrt(norm)
       xvec = xvec * normfac
       matx = matx * normfac
       metx = metx * normfac

       ! Re-estimate eigenvalue
       oldeval = eval
       eval = ddot(n,metx,1,matx,1)
       ! ndmh: prevent WIN32 desync as above
       call comms_bcast(pub_root_node_id,eval)

       ! Check convergence (to better than tol)
       if (abs(eval - oldeval) < local_tol) then
          conv = conv + 1
       else
          conv = 0
       end if
       ! ndmh: prevent WIN32 desync as above
       call comms_bcast(pub_root_node_id,conv)
       if (conv > cgreset .or. abs(gdotg) < epsilon(1.0_DP)) exit

    end do

    ! Deallocate workspace
    call internal_dealloc_work


  contains

    !==========================================================================!
    ! This subroutine performs a matrix-vector product between a block-sparse  !
    ! matrix and a single dense vector.                                        !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   yvec (output) : the result vector                                      !
    !   amat (input)  : the block sparse matrix                                !
    !   xvec (input)  : the original vector                                    !
    !--------------------------------------------------------------------------!
    ! Written by Peter Haynes, March 2004.                                     !
    ! Revised for distributed data, June 2006.                                 !
    ! Modified for SPAM3 by Nicholas Hine, June 2009.                          !
    !==========================================================================!

    subroutine internal_matvec(yvec,amat,xvec)

      implicit none

      ! Arguments
      real(kind=DP), intent(out) :: yvec(:)   ! result vector
      type(SPAM3), intent(in) :: amat         ! block sparse matrix
      real(kind=DP), intent(in) :: xvec(:)    ! input vector

      ! Local variables
      integer :: alib        ! Library entry for amat
      integer :: blks        ! Identifier for row blocking scheme
      integer :: iblk        ! Block-column counter
      integer :: loc_iblk    ! Block-column counter local to this node
      integer :: jblk        ! Block-row counter
      integer :: idx         ! Index
      integer :: ielems      ! Number of elems in block-column
      integer :: jelems      ! Number of elems in block-row
      integer :: ielem       ! Row elem counter for input vector
      integer :: jelem       ! Row elem counter for result vector
      integer :: seg         ! Segment index
      integer :: seg_start   ! Start position of this segment in the index

      ! Get library entry
      alib = amat%lib
      blks = library(alib)%row_blks

      ! Set result vector to zero
      yvec(:) = 0.0_DP

      ! The product to be formed is:
      !  y  = A   x
      !   j    ji  i

      ! Loop over segments of the matrix on this node
      do seg=0,pub_total_num_nodes-1

         if (library(alib)%seg_info(s_type,seg)==SEG_DENSE) then

            ! Do matrix-vector multiplication
            ielems = num_elems_on_node(pub_my_node_id,blks)
            ielem = first_elem_on_node(pub_my_node_id,blks)
            jelems = num_elems_on_node(seg,blks)
            jelem = first_elem_on_node(seg,blks)
            call dgemv('N',jelems,ielems,1.0_DP, &
                 amat%dmtx(library(alib)%seg_info(s_ptr,seg)), &
                 jelems,xvec(ielem),1,1.0_DP,yvec(jelem),1)

         else if (library(alib)%seg_info(s_type,seg)==SEG_SPARSE) then

            ! Loop over block-columns of matrix and input vector on this node
            seg_start = library(alib)%seg_info(s_idx,seg)
            loc_iblk = 0
            do iblk=my_first_blk,my_last_blk
               loc_iblk = loc_iblk + 1

               ! Number of elems and first elem in this block-column
               ielems = num_elems_on_atom(iblk,blks)
               ielem = first_elem_on_atom(iblk,blks)

               ! Loop over block-rows jblk in block-col iblk of matrix
               do idx=library(alib)%blk_idx(seg_start+loc_iblk-1), &
                    library(alib)%blk_idx(seg_start+loc_iblk)-1
                  jblk = library(alib)%blk_idx(idx)

                  ! Number of elems and first elem in this block-row
                  jelems = num_elems_on_atom(jblk,blks)
                  jelem = first_elem_on_atom(jblk,blks)

                  ! Do block multiplication
                  call dgemv('N',jelems,ielems,1.0_DP, &
                       amat%dmtx(library(alib)%blk_ptr(idx)), &
                       jelems,xvec(ielem),1,1.0_DP,yvec(jelem),1)

               end do  ! Loop over block-rows jblk in block-col iblk of matrix

            end do  ! Loop over block-columns of matrix and input vector

         end if  ! seg_info(s_type,seg)

      end do  ! seg

      ! Sum contributions from all nodes
      call comms_reduce('SUM',yvec)

    end subroutine internal_matvec


    !==========================================================================!
    ! This subroutine allocates workspace for the parent subroutine.           !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   None                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Peter Haynes, March 2004.                                     !
    ! Revised for distributed data, June 2006.                                 !
    !==========================================================================!

    subroutine internal_alloc_work

      implicit none

      ! Local variables
      integer :: ierr    ! Error flag

      ! Allocate workspace
      allocate(xvec(n),stat=ierr)
      call utils_alloc_check('internal_alloc_work &
           &(sparse_extremal_eigenvalue)','xvec',ierr)

      allocate(yvec(n),stat=ierr)
      call utils_alloc_check('internal_alloc_work &
           &(sparse_extremal_eigenvalue)','yvec',ierr)

      allocate(grad(n),stat=ierr)
      call utils_alloc_check('internal_alloc_work &
           &(sparse_extremal_eigenvalue)','grad',ierr)

      allocate(dirn(n),stat=ierr)
      call utils_alloc_check('internal_alloc_work &
           &(sparse_extremal_eigenvalue)','dirn',ierr)

      allocate(matx(n),stat=ierr)
      call utils_alloc_check('internal_alloc_work &
           &(sparse_extremal_eigenvalue)','matx',ierr)

      allocate(metx(n),stat=ierr)
      call utils_alloc_check('internal_alloc_work &
           &(sparse_extremal_eigenvalue)','metx',ierr)

    end subroutine internal_alloc_work


    !==========================================================================!
    ! This subroutine deallocates workspace for the parent subroutine.         !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    !   None                                                                   !
    !--------------------------------------------------------------------------!
    ! Written by Peter Haynes, March 2004.                                     !
    ! Revised for distributed data, June 2006.                                 !
    !==========================================================================!

    subroutine internal_dealloc_work

      implicit none

      ! Local variables
      integer :: ierr    ! Error flag

      ! Deallocate workspace
      deallocate(metx,stat=ierr)
      call utils_dealloc_check('internal_dealloc_work &
           &(sparse_extremal_eigenvalue)','metx',ierr)

      deallocate(matx,stat=ierr)
      call utils_dealloc_check('internal_dealloc_work &
           &(sparse_extremal_eigenvalue)','matx',ierr)

      deallocate(dirn,stat=ierr)
      call utils_dealloc_check('internal_dealloc_work &
           &(sparse_extremal_eigenvalue)','dirn',ierr)

      deallocate(grad,stat=ierr)
      call utils_dealloc_check('internal_dealloc_work &
           &(sparse_extremal_eigenvalue)','grad',ierr)

      deallocate(yvec,stat=ierr)
      call utils_dealloc_check('internal_dealloc_work &
           &(sparse_extremal_eigenvalue)','yvec',ierr)

      deallocate(xvec,stat=ierr)
      call utils_dealloc_check('internal_dealloc_work &
           &(sparse_extremal_eigenvalue)','xvec',ierr)

    end subroutine internal_dealloc_work

  end subroutine sparse_extremal_eigenvalue


  !============================================================================!
  ! This subroutine enforces the sparsity pattern of a given sparse matrix, in !
  ! the cases where some or all of its segments are declared 'dense'.          !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat  (inout) : The sparse matrix                                         !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  subroutine sparse_enforce_sparsity(mat)

    use comms, only: comms_abort, pub_my_node_id, pub_total_num_nodes
    use parallel_strategy, only: pub_num_atoms_on_node
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat   ! The sparse matrix

    ! Local variables
    integer :: lib             ! Library pointers for x and y
    integer :: row_blks        ! Identifier for column blocking scheme
    integer :: col_blks        ! Identifier for row blocking scheme
    integer :: seg             ! Segment index
    integer :: seg_start       ! Start position of this segment in the index
    integer :: iblk            ! Block-column counter
    integer :: loc_iblk        ! Block-column counter local to this node
    integer :: jblk            ! Block-row counter
    integer :: iidx            ! Index, Pointer
    integer :: ptr,bptr        ! Matrix and block pointers
    integer :: col_nzb         ! Number of nonzero blocks in this column
    integer :: ielems,jelems   ! Numbers of elements in each block
    integer :: ielem           ! ELement in block loop counter
    integer :: nrows           ! Total number of rows in this segment
    integer :: ierr            ! Error flag

    real(kind=DP), allocatable :: dblk_col(:)
    complex(kind=DP), allocatable :: zblk_col(:)

    ! Get library entries
    lib = mat%lib
    row_blks = library(lib)%row_blks
    col_blks = library(lib)%col_blks

    ! Allocate temporary storage for block column
    if (mat%iscmplx) then
       allocate(zblk_col(maxval(num_elems_on_node(:,row_blks)) * &
            maxval(num_elems_on_atom(:,col_blks))),stat=ierr)
       call utils_alloc_check('sparse_enforce_sparsity','zblk_col',ierr)
       zblk_col(:) = (0.0_DP,0.0_DP)
    else
       allocate(dblk_col(maxval(num_elems_on_node(:,row_blks)) * &
            maxval(num_elems_on_atom(:,col_blks))),stat=ierr)
       call utils_alloc_check('sparse_enforce_sparsity','dblk_col',ierr)
       dblk_col(:) = 0.0_DP
    end if

    ! Loop over segments
    do seg=0,pub_total_num_nodes-1

       ! Only anything to do if this is a dense segment
       if (library(lib)%seg_info(s_type,seg)/=SEG_DENSE) cycle

       seg_start = library(lib)%seg_info(s_idx,seg)
       nrows = num_elems_on_node(seg,row_blks)

       ! Loop over block columns on this node
       loc_iblk = 0
       do iblk=my_first_blk,my_last_blk
          loc_iblk = loc_iblk + 1

          ielems = num_elems_on_atom(iblk,col_blks)
          if (ielems==0) cycle

          ! Find number of nonzero blocks in this column of this segment
          col_nzb = library(lib)%blk_idx(seg_start+loc_iblk) - &
               library(lib)%blk_idx(seg_start+loc_iblk-1)
          ! Move on if all the blocks are nonzero (nothing to do)
          if (col_nzb == pub_num_atoms_on_node(seg)) cycle

          ! Otherwise, go through and copy the blocks in the index to blk_col
          do iidx=library(lib)%blk_idx(seg_start+loc_iblk-1), &
               library(lib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(lib)%blk_idx(iidx)
             jelems = num_elems_on_atom(jblk,row_blks)
             ptr = library(lib)%blk_ptr(iidx)
             bptr = first_elem_on_atom(jblk,row_blks) - &
                  first_elem_on_node(seg,row_blks) + 1

             if (mat%iscmplx) then
                do ielem=1,ielems
                   zblk_col(bptr:bptr+jelems-1) = mat%zmtx(ptr:ptr+jelems-1)
                   ptr = ptr + nrows
                   bptr = bptr + nrows
                end do
             else
                do ielem=1,ielems
                   dblk_col(bptr:bptr+jelems-1) = mat%dmtx(ptr:ptr+jelems-1)
                   ptr = ptr + nrows
                   bptr = bptr + nrows
                end do
             end if

          end do

          ! Copy blk_col back to mtx
          ptr = library(lib)%seg_info(s_ptr,seg) + nrows * &
               (first_elem_on_atom(iblk,col_blks) - &
               first_elem_on_node(pub_my_node_id,col_blks))
          bptr = 1
          if (mat%iscmplx) then
             do ielem=1,ielems
                mat%zmtx(ptr:ptr+nrows-1) = zblk_col(bptr:bptr+nrows-1)
                ptr = ptr + nrows
                bptr = bptr + nrows
             end do
          else
             do ielem=1,ielems
                mat%dmtx(ptr:ptr+nrows-1) = dblk_col(bptr:bptr+nrows-1)
                ptr = ptr + nrows
                bptr = bptr + nrows
             end do
          end if

          ! Reset nonzero blocks of blk_col to zero
          do iidx=library(lib)%blk_idx(seg_start+loc_iblk-1), &
               library(lib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(lib)%blk_idx(iidx)
             jelems = num_elems_on_atom(jblk,row_blks)
             bptr = first_elem_on_atom(jblk,row_blks) - &
                  first_elem_on_node(seg,row_blks) + 1

             if (mat%iscmplx) then
                do ielem=1,ielems
                   zblk_col(bptr:bptr+jelems-1) = (0.0_DP,0.0_DP)
                   bptr = bptr + nrows
                end do
             else
                do ielem=1,ielems
                   dblk_col(bptr:bptr+jelems-1) = 0.0_DP
                   bptr = bptr + nrows
                end do
             end if
          end do

       end do  ! iblk

    end do  ! seg

    if (mat%iscmplx) then
       deallocate(zblk_col,stat=ierr)
       call utils_dealloc_check('sparse_enforce_sparsity','zblk_col',ierr)
    else
       deallocate(dblk_col,stat=ierr)
       call utils_dealloc_check('sparse_enforce_sparsity','dblk_col',ierr)
    end if

  end subroutine sparse_enforce_sparsity


  !==========================================================================!
  ! This subroutine initialises the inverse of a sparse matrix so that it is !
  ! ready to be passed to sparse_hotelling_invert and inverted using the     !
  ! Hotelling Algorithm.                                                     !
  ! Based on the theory in T. Ozaki, Phys. Rev. B., vol 64, page 195110.     !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   inv_mat  (input) : The inverse to initialise                           !
  !   mat      (input) : The matrix to be inverted                           !
  !--------------------------------------------------------------------------!
  ! Originally written by Chris-Kriton Skylaris on 31/7/2003.                !
  ! Modified to support SPAM 2 by Peter Haynes on 21/7/2004.                 !
  ! Adapted for SPAM3 matrices by Nicholas Hine, June 2009.                  !
  ! Moved to sparse_mod by Nicholas Hine, September 2009.                    !
  !==========================================================================!

  subroutine sparse_hotelling_init(inv_mat,mat)

    use comms, only: pub_on_root, comms_reduce, pub_my_node_id
    use constants, only: DP, stdout
    use simulation_cell, only : pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    type(SPAM3), intent(inout) :: inv_mat
    type(SPAM3), intent(in)  :: mat

    ! cks: local variables
    real(kind=DP) :: sigma, sum_col
    real(kind=DP), allocatable, dimension(:) :: drow
    complex(kind=DP), allocatable, dimension(:) :: zrow
    integer :: nn, row, col
    integer :: ierr
    integer :: ilib, col_blks


    ! Allocate workspace
    ilib = mat%lib
    nn = library(ilib)%nrows
    col_blks = library(ilib)%col_blks
    if (library(ilib)%row_blks /= col_blks) then
       ! Error
    end if
    if (mat%iscmplx) then
       allocate(zrow(nn),stat=ierr)
       call utils_alloc_check('sparse_hotelling_init','zrow',ierr)
       zrow = 0.0_DP
    else
       allocate(drow(nn),stat=ierr)
       call utils_alloc_check('sparse_hotelling_init','drow',ierr)
       drow = 0.0_DP
    end if

    sigma = 0.0_DP

    ! Loop over all cols of S on this processor
    do col=first_elem_on_node(pub_my_node_id,col_blks), &
         first_elem_on_node(pub_my_node_id+1,col_blks)-1

       if (mat%iscmplx) then

          ! Get column of S
          call sparse_get_col(zrow,mat,col)

          ! Sum column of S
          sum_col = 0.0_DP
          do row=1,nn
             sum_col = sum_col + abs(zrow(row))
          end do
          sigma = max(sigma,sum_col)

          ! Clear column of S
          call sparse_clr_col(zrow,mat,col)

       else

          ! Get column of S
          call sparse_get_col(drow,mat,col)

          ! Sum column of S
          sum_col = 0.0_DP
          do row=1,nn
             sum_col = sum_col + abs(drow(row))
          end do
          sigma = max(sigma,sum_col)

          ! Clear column of S
          call sparse_clr_col(drow,mat,col)

       end if

    end do

    call comms_reduce('MAX',sigma)

    ! Deallocate workspace
    if (mat%iscmplx) then
       deallocate(zrow,stat=ierr)
       call utils_dealloc_check('sparse_hotelling_init','zrow',ierr)
    else
       deallocate(drow,stat=ierr)
       call utils_dealloc_check('sparse_hotelling_init','drow',ierr)
    end if

    ! Copy mat into inv_mat
    call sparse_copy(inv_mat,mat)

    if (sigma > epsilon(1.0_DP)) then
       sigma = 1.0_DP / (sigma*sigma)
    else
       sigma = 0.001_DP
       if (pub_on_root) write(stdout,'(a)') &
            'WARNING in sparse_hotelling_init: zero overlap matrix'
    end if

    ! cks: scale by sigma
    call sparse_scale(inv_mat,sigma)


  end subroutine sparse_hotelling_init

  !==========================================================================!
  ! This subroutine applies successive quadratically convergent Hotelling    !
  ! iterations to improve an approximate inverse overlap matrix, based on    !
  ! the theory in the paper by T. Ozaki, Phys. Rev. B., vol 64, page 195110. !
  ! The iterations continue (unless they exceed the num_iter input           !
  ! parameter) until convergence to machine precision is reached - but       !
  ! taking into account limitations arising from the truncation of the       !
  ! inverse overlap matrix.                                                  !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   inv_mat  (input) : The inverse to calculate                            !
  !   mat      (input) : The matrix to be inverted                           !
  !--------------------------------------------------------------------------!
  ! Originally written by Chris-Kriton Skylaris on 31/7/2003                 !
  ! Modified to support SPAM 2 by Peter Haynes on 21/7/2004.                 !
  ! Modified to test for convergence by Chris-Kriton Skylaris on 4/10/2004.  !
  ! Minor modifications for parallel SPAM 2 by Peter Haynes                  !
  ! Adapted to use SPAM3 matrices by Nicholas Hine, June 2009.               !
  ! Moved to sparse_mod by Nicholas Hine, September 2009.                    !
  !==========================================================================!

  subroutine sparse_hotelling_invert(inv_mat,mat,show_output, &
       max_resid_converged, num_iter, final_max_resid)

    use comms, only: pub_on_root,comms_abort
    use constants, only: DP, stdout
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: inv_mat
    type(SPAM3), intent(in) :: mat
    real(kind=DP), intent(in) :: max_resid_converged
    logical, intent(in) :: show_output
    integer, intent(in) :: num_iter
    real(kind=DP), intent(out), optional :: final_max_resid

    ! Local Variables
    integer :: iter
    logical :: quit_early  ! pdh: flag
    real(kind=DP) :: max_resid  ! maximum value of residual
    real(kind=DP) :: frob_norm
    real(kind=DP) :: previous_frob_norm
    type(SPAM3) :: mim,inv_tmp

    ! Start timer
    call timer_clock('sparse_hotelling_invert',1)

    ! Initialisations
    previous_frob_norm = huge(1.0_DP)
    quit_early = .false.

    ! Allocate workspace
    call sparse_create(inv_tmp,inv_mat)
    call sparse_create(mim,mat,inv_mat)

    ! Start of Hotelling iteration loop
    hotel_sparse_loop: do iter=1,num_iter

       ! MIM := M*(M^-1)
       call sparse_product(mim,mat,inv_mat)

       ! cks: MIM := I-M*(M^-1)
       call sparse_scale(mim,-1.0_DP,1.0_DP)

       ! cks: Calculate Frobenius norm of MIM -->0
       frob_norm = sparse_rms_element(mim) * &
            sqrt(sparse_num_element(mim))

       ! cks: Maximum element of residual I-M*M_n^-1
       max_resid = sparse_max_abs_element(mim)

       if (show_output.and.pub_on_root) then
          write(stdout,'(t12,i5,tr5,e16.8,tr3,e16.8)') iter, frob_norm, &
               max_resid
       end if

       ! cks: Test for convergence due to machine precision
       ! cks: or inverse overlap truncation
       if (frob_norm >= previous_frob_norm .or. &
            max_resid <= max_resid_converged) then
          quit_early = .true.
          exit
       else
          previous_frob_norm = frob_norm
       endif

       ! cks: MIM := I + MIM := 2I -M*(M^-1)
       call sparse_scale(mim,1.0_DP,1.0_DP)

       ! (M^-1) := (M^-1)*(2*I-M*(M^-1))
       call sparse_copy(inv_tmp,inv_mat)
       call sparse_product(inv_mat,inv_tmp,mim)

    end do hotel_sparse_loop

    if (present(final_max_resid)) final_max_resid = max_resid

    if (pub_on_root .and. .not. quit_early) write(stdout,*) &
         'WARNING: max Hotelling iterations exceeded. Last norm=',frob_norm

    ! Deallocate workspace
    call sparse_destroy(mim)
    call sparse_destroy(inv_tmp)

    ! Stop timer
    call timer_clock('sparse_hotelling_invert',2)

  end subroutine sparse_hotelling_invert

! These routines only need to exist if the code is being compiled with ScaLAPACK
#ifdef SCALAPACK

  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3 into a dense     !
  ! array in the BLACS distributed format suitable for use with ScaLAPACK      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination BLACS matrix (must be real)            !
  !   src    (input)  : The source SPAM3 matrix (must be real)                 !
  !   ld     (input)  : Leading dimension of the array dest                    !
  !   nc     (input)  : Number of columns of the array dest                    !
  !   desc   (input)  : BLACS descriptor for the array dest                    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, October 2009.                                    !
  !============================================================================!

  subroutine sparse_spam3toblacs_real(dest,src,ld,nc,desc)

    use comms, only: comms_abort, comms_bcast, pub_my_node_id, pub_on_root, &
         pub_total_num_nodes
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    real(kind=DP),intent(out) :: dest(ld,nc)
    type(SPAM3),intent(in) :: src
    integer,intent(in) :: ld,nc
    integer,intent(in) :: desc(50)

    ! Local Variables
    integer :: ierr
    integer :: srclib        ! Library entry for src
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! elem in atom block iblk counter
    integer :: jelemonat     ! elem in atom block jblk counter
    integer :: ielem,jelem   ! elem counters
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in this block/segment
    integer :: node          ! Node loop counter
    integer :: seg           ! Segment index
    integer :: seg_type      ! Start position of this segment in the index
    integer :: seg_start     ! Start position of this segment in the index
    real(kind=DP), allocatable :: buffer(:,:),loc_buffer(:,:)

    ! Get library entry for src
    srclib = src%lib
    row_blks = library(srclib)%row_blks
    col_blks = library(srclib)%col_blks

    ! Check arguments
    if (src%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_spam3toblacs_real: real matrices only.'
       call comms_abort
    end if

    ! Allocate buffers
    nrows = library(srclib)%nrows
    allocate(buffer(nrows,maxval(num_elems_on_node(:,col_blks))),stat=ierr)
    call utils_alloc_check('sparse_spam3toblacs_real','buffer',ierr)
    allocate(loc_buffer(nrows,num_elems_on_node(pub_my_node_id,col_blks)), &
         stat=ierr)
    call utils_alloc_check('sparse_spam3toblacs_real','loc_buffer',ierr)

    ! Zero destination matrix and buffer
    dest(:,:) = 0.0_DP
    buffer(:,:) = 0.0_DP
    loc_buffer(:,:) = 0.0_DP

    ! Loop over the segments of src on this node
    do seg=0,pub_total_num_nodes-1

       ! Cycle if this segment is blank
       seg_type = library(srclib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = num_elems_on_node(seg,row_blks)

       ! Loop over block-columns of src on this node
       seg_start = library(srclib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=my_first_blk,my_last_blk
          loc_iblk = loc_iblk + 1
          ielem = first_elem_on_atom(iblk,col_blks) - &
               first_elem_on_node(pub_my_node_id,col_blks) + 1
          ielems = num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of src
          do idx=library(srclib)%blk_idx(seg_start+loc_iblk-1), &
               library(srclib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(srclib)%blk_idx(idx)
             jelem = first_elem_on_atom(jblk,row_blks)
             jelems = num_elems_on_atom(jblk,row_blks)

             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ptr = library(srclib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                do jelemonat=0,jelems-1   ! Element row in block
                   loc_buffer(jelem+jelemonat,ielem+ielemonat) = &
                        src%dmtx(ptr+ielemonat*nrows+jelemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of src

       end do  ! Loop over block-columns of src

    end do  ! seg

    ! Broadcast results across all nodes
    nrows = library(srclib)%nrows
    do node=0,pub_total_num_nodes-1

       ! Copy local buffer to buffer if node is local node
       if (node==pub_my_node_id) then
          ielems = num_elems_on_node(pub_my_node_id,col_blks)
          buffer(1:nrows,1:ielems) = loc_buffer(1:nrows,1:ielems)
       end if

       ! Broadcast this part of the matrix
       call comms_bcast(node,buffer(1,1), &
            num_elems_on_node(node,col_blks)*library(srclib)%nrows)

       ! Loop over rows of dest, setting values
       do ielem=first_elem_on_node(node,col_blks), &
            first_elem_on_node(node+1,col_blks)-1
          ielems = ielem - first_elem_on_node(node,col_blks) + 1
          do jelem=1,nrows
             call pdelset(dest,jelem,ielem,desc,buffer(jelem,ielems))
          end do
       end do

    end do  ! Loop over nodes

    ! Dellocate buffers
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('sparse_spam3toblacs_real','buffer',ierr)
    deallocate(loc_buffer,stat=ierr)
    call utils_dealloc_check('sparse_spam3toblacs_real','loc_buffer',ierr)

  end subroutine sparse_spam3toblacs_real

  !============================================================================!
  ! This subroutine converts a block sparse matrix in a SPAM3 into a dense     !
  ! array in the BLACS distributed format suitable for use with ScaLAPACK      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination BLACS matrix (must be complex)         !
  !   src    (input)  : The source SPAM3 matrix (must be complex)              !
  !   ld     (input)  : Leading dimension of the array dest                    !
  !   nc     (input)  : Number of columns of the array dest                    !
  !   desc   (input)  : BLACS descriptor for the array dest                    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, October 2009.                                    !
  !============================================================================!

  subroutine sparse_spam3toblacs_complex(dest,src,ld,nc,desc)

    use comms, only: comms_abort, comms_bcast, pub_my_node_id, pub_on_root, &
         pub_total_num_nodes
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    complex(kind=DP),intent(out) :: dest(ld,nc)
    type(SPAM3),intent(in) :: src
    integer,intent(in) :: ld,nc
    integer,intent(in) :: desc(50)

    ! Local Variables
    integer :: ierr
    integer :: srclib        ! Library entry for src
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! elem in atom block iblk counter
    integer :: jelemonat     ! elem in atom block jblk counter
    integer :: ielem,jelem   ! elem counters
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in this block/segment
    integer :: node          ! Node loop counter
    integer :: seg           ! Segment index
    integer :: seg_type      ! Start position of this segment in the index
    integer :: seg_start     ! Start position of this segment in the index
    complex(kind=DP), allocatable :: buffer(:,:),loc_buffer(:,:)

    ! Get library entry for src
    srclib = src%lib
    row_blks = library(srclib)%row_blks
    col_blks = library(srclib)%col_blks

    ! Check arguments
    if (.not.src%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_spam3toblacs_complex: complex matrices only.'
       call comms_abort
    end if

    ! Allocate buffers
    nrows = library(srclib)%nrows
    allocate(buffer(nrows,maxval(num_elems_on_node(:,col_blks))),stat=ierr)
    call utils_alloc_check('sparse_spam3toblacs_complex','buffer',ierr)
    allocate(loc_buffer(nrows,num_elems_on_node(pub_my_node_id,col_blks)), &
         stat=ierr)
    call utils_alloc_check('sparse_spam3toblacs_complex','loc_buffer',ierr)

    ! Zero destination matrix and buffer
    dest(:,:) = cmplx(0.0_DP,0.0_DP)
    buffer(:,:) = cmplx(0.0_DP,0.0_DP)
    loc_buffer(:,:) = cmplx(0.0_DP,0.0_DP)

    ! Loop over the segments of src on this node
    do seg=0,pub_total_num_nodes-1

       ! Cycle if this segment is blank
       seg_type = library(srclib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = num_elems_on_node(seg,row_blks)

       ! Loop over block-columns of src on this node
       seg_start = library(srclib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=my_first_blk,my_last_blk
          loc_iblk = loc_iblk + 1
          ielem = first_elem_on_atom(iblk,col_blks) - &
               first_elem_on_node(pub_my_node_id,col_blks) + 1
          ielems = num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of src
          do idx=library(srclib)%blk_idx(seg_start+loc_iblk-1), &
               library(srclib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(srclib)%blk_idx(idx)
             jelem = first_elem_on_atom(jblk,row_blks)
             jelems = num_elems_on_atom(jblk,row_blks)

             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ptr = library(srclib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                do jelemonat=0,jelems-1   ! Element row in block
                   loc_buffer(jelem+jelemonat,ielem+ielemonat) = &
                        src%zmtx(ptr+ielemonat*nrows+jelemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of src

       end do  ! Loop over block-columns of src

    end do  ! seg

    ! Broadcast results across all nodes
    nrows = library(srclib)%nrows
    do node=0,pub_total_num_nodes-1

       ! Copy local buffer to buffer if node is local node
       if (node==pub_my_node_id) then
          ielems = num_elems_on_node(pub_my_node_id,col_blks)
          buffer(1:nrows,1:ielems) = loc_buffer(1:nrows,1:ielems)
       end if

       ! Broadcast this part of the matrix
       call comms_bcast(node,buffer(1,1), &
            num_elems_on_node(node,col_blks)*library(srclib)%nrows)

       ! Loop over rows of dest, setting values
       do ielem=first_elem_on_node(node,col_blks), &
            first_elem_on_node(node+1,col_blks)-1
          ielems = ielem - first_elem_on_node(node,col_blks) + 1
          do jelem=1,nrows
             call pzelset(dest,jelem,ielem,desc,buffer(jelem,ielems))
          end do
       end do

    end do  ! Loop over nodes

    ! Dellocate buffers
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('sparse_spam3toblacs_complex','buffer',ierr)
    deallocate(loc_buffer,stat=ierr)
    call utils_dealloc_check('sparse_spam3toblacs_complex','loc_buffer',ierr)

  end subroutine sparse_spam3toblacs_complex

  !============================================================================!
  ! This subroutine converts a dense array in the BLACS distributed format     !
  ! suitable for use with ScaLAPACK into a SPAM3 block sparse matrix           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination SPAM3 matrix (must be real)            !
  !   src    (input)  : The source BLACS dense matrix (real)                   !
  !   ld     (input)  : Leading dimension of the array dest                    !
  !   nc     (input)  : Number of columns of the array dest                    !
  !   desc   (input)  : BLACS descriptor for the array dest                    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, February 2010.                                   !
  !============================================================================!

  subroutine sparse_blacstospam3_real(dest,src,ld,nc,desc)

    use comms, only: comms_abort, comms_reduce, pub_my_node_id, pub_on_root, &
         pub_total_num_nodes
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: dest
    real(kind=DP),intent(in) :: src(ld,nc)
    integer,intent(in) :: ld,nc
    integer,intent(in) :: desc(50)

    ! Local Variables
    integer :: ierr
    integer :: destlib       ! Library entry for dest
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! elem in atom block iblk counter
    integer :: jelemonat     ! elem in atom block jblk counter
    integer :: ielem,jelem   ! elem counters
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in this block/segment
    integer :: node          ! Node loop counter
    integer :: seg           ! Segment index
    integer :: seg_type      ! Start position of this segment in the index
    integer :: seg_start     ! Start position of this segment in the index
    real(kind=DP), allocatable :: buffer(:,:),loc_buffer(:,:)

    ! Get library entry for src
    destlib = dest%lib
    row_blks = library(destlib)%row_blks
    col_blks = library(destlib)%col_blks
    nrows = library(destlib)%nrows

    ! Check arguments
    if (dest%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_blacstospam3_real: real matrices only.'
       call comms_abort
    end if

    ! Allocate buffers
    allocate(buffer(nrows,maxval(num_elems_on_node(:,col_blks))),stat=ierr)
    call utils_alloc_check('sparse_blacstospam3_real','buffer',ierr)
    allocate(loc_buffer(nrows,num_elems_on_node(pub_my_node_id,col_blks)), &
         stat=ierr)
    call utils_alloc_check('sparse_blacstospam3_real','loc_buffer',ierr)

    ! Broadcast results across all nodes
    nrows = library(destlib)%nrows
    do node=0,pub_total_num_nodes-1

       ! Zero buffer
       buffer(:,:) = 0.0_DP

       ! Loop over rows of src, getting values
       do ielem=first_elem_on_node(node,col_blks), &
            first_elem_on_node(node+1,col_blks)-1
          ielems = ielem - first_elem_on_node(node,col_blks) + 1
          do jelem=1,nrows
             call pdelget('O','O',buffer(jelem,ielems),src,jelem,ielem,desc)
          end do
       end do

       ! Sum this part of the matrix over nodes
       call comms_reduce('SUM',buffer, &
            num_elems_on_node(node,col_blks)*library(destlib)%nrows)

       ! Copy buffer to local buffer if node is local node
       if (node==pub_my_node_id) then
          ielems = num_elems_on_node(pub_my_node_id,col_blks)
          loc_buffer(1:nrows,1:ielems) = buffer(1:nrows,1:ielems)
       end if

    end do  ! Loop over nodes

    ! Loop over the segments of dest on this node
    do seg=0,pub_total_num_nodes-1

       ! Cycle if this segment is blank
       seg_type = library(destlib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = num_elems_on_node(seg,row_blks)

       ! Loop over block-columns of dest on this node
       seg_start = library(destlib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=my_first_blk,my_last_blk
          loc_iblk = loc_iblk + 1
          ielem = first_elem_on_atom(iblk,col_blks) - &
               first_elem_on_node(pub_my_node_id,col_blks) + 1
          ielems = num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of dest
          do idx=library(destlib)%blk_idx(seg_start+loc_iblk-1), &
               library(destlib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(destlib)%blk_idx(idx)
             jelem = first_elem_on_atom(jblk,row_blks)
             jelems = num_elems_on_atom(jblk,row_blks)

             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ptr = library(destlib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                do jelemonat=0,jelems-1   ! Element row in block
                   dest%dmtx(ptr+ielemonat*nrows+jelemonat) = &
                        loc_buffer(jelem+jelemonat,ielem+ielemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of src

       end do  ! Loop over block-columns of src

    end do  ! seg

    ! Dellocate buffers
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('sparse_blacstospam3_real','buffer',ierr)
    deallocate(loc_buffer,stat=ierr)
    call utils_dealloc_check('sparse_blacstospam3_real','loc_buffer',ierr)

  end subroutine sparse_blacstospam3_real

  !============================================================================!
  ! This subroutine converts a dense array in the BLACS distributed format     !
  ! suitable for use with ScaLAPACK into a SPAM3 block sparse matrix           !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   dest   (output) : The destination SPAM3 matrix (must be complex)         !
  !   src    (input)  : The source BLACS dense matrix (complex)                !
  !   ld     (input)  : Leading dimension of the array dest                    !
  !   nc     (input)  : Number of columns of the array dest                    !
  !   desc   (input)  : BLACS descriptor for the array dest                    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, February 2010.                                   !
  !============================================================================!

  subroutine sparse_blacstospam3_complex(dest,src,ld,nc,desc)

    use comms, only: comms_abort, comms_reduce, pub_my_node_id, pub_on_root, &
         pub_total_num_nodes
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: dest
    complex(kind=DP),intent(in) :: src(ld,nc)
    integer,intent(in) :: ld,nc
    integer,intent(in) :: desc(50)

    ! Local Variables
    integer :: ierr
    integer :: destlib       ! Library entry for dest
    integer :: row_blks      ! Identifier for column blocking scheme
    integer :: col_blks      ! Identifier for row blocking scheme
    integer :: iblk,jblk     ! Block counters
    integer :: loc_iblk      ! Local counter for iblk
    integer :: idx           ! Index
    integer :: ptr           ! Pointer
    integer :: ielemonat     ! elem in atom block iblk counter
    integer :: jelemonat     ! elem in atom block jblk counter
    integer :: ielem,jelem   ! elem counters
    integer :: ielems,jelems ! Numbers of element rows/cols in this block
    integer :: nrows         ! Total number of rows in this block/segment
    integer :: node          ! Node loop counter
    integer :: seg           ! Segment index
    integer :: seg_type      ! Start position of this segment in the index
    integer :: seg_start     ! Start position of this segment in the index
    complex(kind=DP), allocatable :: buffer(:,:),loc_buffer(:,:)

    ! Get library entry for src
    destlib = dest%lib
    row_blks = library(destlib)%row_blks
    col_blks = library(destlib)%col_blks
    nrows = library(destlib)%nrows

    ! Check arguments
    if (.not.dest%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_blacstospam3_real: complex matrices only.'
       call comms_abort
    end if

    ! Allocate buffers
    allocate(buffer(nrows,maxval(num_elems_on_node(:,col_blks))),stat=ierr)
    call utils_alloc_check('sparse_blacstospam3_complex','buffer',ierr)
    allocate(loc_buffer(nrows,num_elems_on_node(pub_my_node_id,col_blks)), &
         stat=ierr)
    call utils_alloc_check('sparse_blacstospam3_complex','loc_buffer',ierr)

    ! Broadcast results across all nodes
    nrows = library(destlib)%nrows
    do node=0,pub_total_num_nodes-1

       ! Zero buffer
       buffer(:,:) = 0.0_DP

       ! Loop over rows of src, getting values
       do ielem=first_elem_on_node(node,col_blks), &
            first_elem_on_node(node+1,col_blks)-1
          ielems = ielem - first_elem_on_node(node,col_blks) + 1
          do jelem=1,nrows
             call pzelget('O','O',buffer(jelem,ielems),src,jelem,ielem,desc)
          end do
       end do

       ! Sum this part of the matrix over nodes
       call comms_reduce('SUM',buffer, &
            num_elems_on_node(node,col_blks)*library(destlib)%nrows)

       ! Copy buffer to local buffer if node is local node
       if (node==pub_my_node_id) then
          ielems = num_elems_on_node(pub_my_node_id,col_blks)
          loc_buffer(1:nrows,1:ielems) = buffer(1:nrows,1:ielems)
       end if

    end do  ! Loop over nodes

    ! Loop over the segments of dest on this node
    do seg=0,pub_total_num_nodes-1

       ! Cycle if this segment is blank
       seg_type = library(destlib)%seg_info(s_type,seg)
       if (seg_type == SEG_BLANK) cycle
       if (seg_type == SEG_DENSE) nrows = num_elems_on_node(seg,row_blks)

       ! Loop over block-columns of dest on this node
       seg_start = library(destlib)%seg_info(s_idx,seg)
       loc_iblk = 0
       do iblk=my_first_blk,my_last_blk
          loc_iblk = loc_iblk + 1
          ielem = first_elem_on_atom(iblk,col_blks) - &
               first_elem_on_node(pub_my_node_id,col_blks) + 1
          ielems = num_elems_on_atom(iblk,col_blks)

          ! Loop over block-rows jblk in column iblk of dest
          do idx=library(destlib)%blk_idx(seg_start+loc_iblk-1), &
               library(destlib)%blk_idx(seg_start+loc_iblk)-1
             jblk = library(destlib)%blk_idx(idx)
             jelem = first_elem_on_atom(jblk,row_blks)
             jelems = num_elems_on_atom(jblk,row_blks)

             if (seg_type == SEG_SPARSE) nrows = jelems

             ! Copy block
             ptr = library(destlib)%blk_ptr(idx)
             do ielemonat=0,ielems-1      ! Element col in block
                do jelemonat=0,jelems-1   ! Element row in block
                   dest%zmtx(ptr+ielemonat*nrows+jelemonat) = &
                        loc_buffer(jelem+jelemonat,ielem+ielemonat)
                end do
             end do

          end do  ! Loop over block-rows in column iblk of src

       end do  ! Loop over block-columns of src

    end do  ! seg

    ! Dellocate buffers
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check('sparse_blacstospam3_complex','buffer',ierr)
    deallocate(loc_buffer,stat=ierr)
    call utils_dealloc_check('sparse_blacstospam3_complex','loc_buffer',ierr)

  end subroutine sparse_blacstospam3_complex

#endif

  !============================================================================!
  ! This subroutine is a wrapper to write a single sparse matrix to file by    !
  ! calling sparse_write_vector with a single-element array.                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The sparse matrix to be written                       !
  !   filename (input) : The filename to use                                   !
  !   unit     (input) : The (optional) unit number to use if file is open     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009.                                       !
  !============================================================================!

  subroutine sparse_write_scalar(mat,filename,unit)

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat            ! Matrix to write
    character(len=*), intent(in) :: filename  ! Filename to use
    integer, optional, intent(in) :: unit     ! I/O unit

    if (present(unit)) then
       call sparse_write_vector((/mat/),filename,unit)
    else
       call sparse_write_vector((/mat/),filename)
    end if

  end subroutine sparse_write_scalar


  !============================================================================!
  ! This subroutine writes a vector SPAM3 sparse matrix to a file.             !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The sparse matrix to be written                       !
  !   filename (input) : The filename to use                                   !
  !   unit     (input) : The (optional) unit number to use if file is open     !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, January 2007.                                     !
  ! Modification for dense matrices by Nicholas Hine, Dec 2007.                !
  ! Significantly modified for SPAM3 by Nicholas Hine, June 2009.              !
  !============================================================================!

  subroutine sparse_write_vector(mat,filename,unit)

    use comms, only: comms_abort, comms_barrier, comms_bcast, comms_send, &
         comms_recv, comms_reduce, pub_my_node_id, pub_on_root, &
         pub_root_node_id, pub_total_num_nodes
    use parallel_strategy, only: pub_orig_atom, pub_first_atom_on_node
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: mat(:)         ! Matrix to write
    character(len=*), intent(in) :: filename  ! Filename to use
    integer, optional, intent(in) :: unit     ! I/O unit

    ! Local variables
    integer :: ncomps               ! Number of components
    integer :: icomp                ! Loop counter for components
    integer :: blks                 ! Identifier for blocking scheme
    integer :: ierr                 ! Error flag
    integer :: ilib                 ! Library entry for mat
    integer :: iunit                ! Unit number for writing file
    integer :: iblk                 ! Atom loop counter
    integer :: iat_orig             ! Atom iblk in original order
    integer :: jblk                 ! Second atom
    integer :: jat_orig             ! Atom jblk in original order
    integer :: loc_iblk             ! Local atom counter for iat
    integer :: iidx                 ! Index loop counter
    integer :: jidx                 ! Index info counter
    integer :: idxlen               ! Index info length
    integer :: node                 ! Processor loop counter
    integer :: datlen               ! Buffer length for data
    integer :: blk_start,blk_end    ! First and last blocks on node
    integer :: ielems               ! Number of elements on atom iat
    integer, allocatable :: block_sizes(:)  ! Block sizes
    integer, allocatable :: idxbuf(:)       ! Comm buffer for index
    integer, allocatable :: idxinfo(:)      ! Index info for writing
    character(len=11) :: tag        ! Field tag
    real(kind=DP), allocatable :: dmtxbuf(:)    ! Comm buffer for real data
    complex(kind=DP), allocatable :: zmtxbuf(:) ! Comm buffer for complex data

    ! Obtain number of components
    ncomps = size(mat)

    ! Obtain library entry for mat
    ilib = mat(1)%lib
    blks = library(ilib)%row_blks

    ! Check that the matrix is square
    if (blks/=library(ilib)%col_blks) then
       write(stdout,'(a)') 'Error in sparse_write_vector: &
            &row and column blocking schemes do not match'
       call comms_abort
    end if

    ! Preliminaries on root node only
    if (pub_on_root) then

       if (present(unit)) then
          iunit = unit
       else
          iunit = utils_unit()
          open(unit=iunit,file=trim(filename),form='unformatted', &
               action='write',iostat=ierr,status='replace')
          if (ierr /= 0) then
             write(stdout,'(3a,i6)') 'Error in sparse_write_vector: &
                  &opening file "',trim(filename),'" failed with code ',ierr
             call comms_abort
          end if
       end if

       ! Version - version*100 so 150 means version 1.5
       tag = 'VERSION'
       write(iunit) tag//'I'
       write(iunit) 1
       write(iunit) nint(file_version*100)

       ! Number of atoms/blocks
       tag = 'NUM_ATOMS'
       write(iunit) tag//'I'
       write(iunit) 1
       write(iunit) library(ilib)%nblk

       ! Number of NGWFs
       tag = 'NUM_NGWFS'
       write(iunit) tag//'I'
       write(iunit) 1
       write(iunit) library(ilib)%nrows

       ! Data type: real or complex
       tag = 'TYPE'
       write(iunit) tag//'I'
       write(iunit) 1
       if (mat(1)%iscmplx) then
          write(iunit) -1
       else
          write(iunit) 0
       end if

       ! Number of components
       tag = 'NUM_COMPS'
       write(iunit) tag//'I'
       write(iunit) 1
       write(iunit) ncomps

       ! List of block sizes (in original order)
       tag = 'BLOCK_SIZES'
       write(iunit) tag//'I'
       write(iunit) library(ilib)%nblk
       allocate(block_sizes(library(ilib)%nblk),stat=ierr)
       call utils_alloc_check('sparse_write_vector','block_sizes',ierr)
       do iblk=1,library(ilib)%nblk
          iat_orig = pub_orig_atom(iblk)
          block_sizes(iat_orig) = num_elems_on_atom(iblk,blks)
       end do
       write(iunit) block_sizes
       deallocate(block_sizes,stat=ierr)
       call utils_dealloc_check('sparse_write_vector','block_sizes',ierr)

       ! End of root node preliminaries
    else
       iunit = -1 ! qoh: Initialise to avoid compiler warning
    end if

    ! Allocate buffers to receive data communicated from nodes
    if (mat(1)%iscmplx) then
       allocate(zmtxbuf(library(ilib)%max_nze*ncomps),stat=ierr)
       call utils_alloc_check('sparse_write_vector','zmtxbuf',ierr)
    else
       allocate(dmtxbuf(library(ilib)%max_nze*ncomps),stat=ierr)
       call utils_alloc_check('sparse_write_vector','dmtxbuf',ierr)
    end if

    ! Global preliminaries
    idxlen = sparse_index_length(mat(1)) + 1
    call comms_reduce('MAX',idxlen)
    allocate(idxinfo(idxlen*3),stat=ierr)
    call utils_alloc_check('sparse_write_vector','idxinfo',ierr)
    allocate(idxbuf(idxlen),stat=ierr)
    call utils_alloc_check('sparse_write_vector','idxbuf',ierr)

    ! Get the index for this node in non-segmented format
    call sparse_generate_index(idxbuf,mat(1))

    ! Read segmented data into unsegmented buffers
    call internal_unsegment

    call comms_barrier

    ! Loop over processors and communicate with root node
    do node=0,pub_total_num_nodes-1
       blk_start = pub_first_atom_on_node(node)
       blk_end = pub_first_atom_on_node(node+1) - 1

       ! Send index of mat stored on node to root node
       if ((pub_my_node_id == node).and.(.not.pub_on_root)) then
          idxbuf(idxlen) = datlen
          call comms_send(pub_root_node_id,idxbuf,idxlen)
       end if
       if (node /= pub_root_node_id .and. pub_on_root) then
          call comms_recv(node,idxbuf,idxlen)
       end if

       ! Send data of mat stored on node to root node
       if (pub_my_node_id == node) then

          if (.not.pub_on_root) then
             if (mat(1)%iscmplx) then
                call comms_send(pub_root_node_id,zmtxbuf,datlen*ncomps)
             else
                call comms_send(pub_root_node_id,dmtxbuf,datlen*ncomps)
             end if
          end if

       end if

       ! Receive data of mat on root node from other node
       if (pub_on_root .and. node /= pub_root_node_id) then
          datlen = idxbuf(idxlen)
          if (mat(1)%iscmplx) then
             call comms_recv(node,zmtxbuf,datlen*ncomps)
          else
             call comms_recv(node,dmtxbuf,datlen*ncomps)
          end if
       end if

       if (pub_on_root) then
          ! Now write out indexing information for this data
          jidx = 0
          loc_iblk = 0
          do iblk=blk_start,blk_end
             loc_iblk = loc_iblk + 1
             iat_orig = pub_orig_atom(iblk)
             ielems = num_elems_on_atom(iblk,blks)

             ! Loop over block-rows in column iat
             do iidx=idxbuf(loc_iblk),idxbuf(loc_iblk+1)-1
                jblk = idxbuf(iidx)
                jat_orig = pub_orig_atom(jblk)

                ! Generate index info: row, column and number of elements
                idxinfo(jidx+1) = jat_orig
                idxinfo(jidx+2) = iat_orig
                idxinfo(jidx+3) = ielems * num_elems_on_atom(jblk,blks)
                jidx = jidx + 3
             end do
          end do
          tag = 'INDEX'
          write(iunit) tag//'I'
          write(iunit) jidx
          write(iunit) idxinfo(1:jidx)

          ! Now the data itself
          tag = 'DATA'
          if (mat(1)%iscmplx) then
             write(iunit) tag//'Z'
             write(iunit) datlen*ncomps
             write(iunit) zmtxbuf(1:datlen*ncomps)
          else
             write(iunit) tag//'D'
             write(iunit) datlen*ncomps
             write(iunit) dmtxbuf(1:datlen*ncomps)
          end if

       end if

       ! Loop over nodes
    end do

    if (pub_on_root) then
       ! Write a terminating tag
       tag = 'ENDFILE'
       write(iunit) tag//'I'
       write(iunit) 1
       write(iunit) 0

       ! Close file for output
       close(unit=iunit,iostat=ierr)
       if (ierr /= 0) then
          write(stdout,'(3a,i6)') 'Error in sparse_write_vector: &
               &closing file "',trim(filename),'" failed with code ',ierr
          call comms_abort
       end if
    end if

    ! Re-sync nodes
    call comms_barrier

    ! Deallocate workspace
    deallocate(idxbuf,stat=ierr)
    call utils_dealloc_check('sparse_write_vector','idxbuf',ierr)
    deallocate(idxinfo,stat=ierr)
    call utils_dealloc_check('sparse_write_vector','idxinfo',ierr)

    if (mat(1)%iscmplx) then
       deallocate(zmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_write_vector','zmtxbuf',ierr)
    else
       deallocate(dmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_write_vector','dmtxbuf',ierr)
    end if

  contains

    subroutine internal_unsegment

      implicit none

      ! Locals
      integer :: iidx,ptr
      integer :: iblk,loc_iblk
      integer :: jblk
      integer :: jelems
      integer :: ielem
      integer :: max_elems_on_atom
      real (kind=DP), allocatable :: dblk(:,:)
      complex (kind=DP), allocatable :: zblk(:,:)

      ! Allocate temporary array to store block
      max_elems_on_atom = maxval(num_elems_on_atom(:,blks))
      if (mat(1)%iscmplx) then
         allocate(zblk(max_elems_on_atom,max_elems_on_atom),stat=ierr)
         call utils_alloc_check('internal_unsegment (sparse_write_vector)', &
              'zblk',ierr)
      else
         allocate(dblk(max_elems_on_atom,max_elems_on_atom),stat=ierr)
         call utils_alloc_check('internal_unsegment (sparse_write_vector)', &
              'dblk',ierr)
      end if

      ptr = 1

      ! Loop over components of the matrix array
      do icomp=1,ncomps

         loc_iblk = 0
         ! Loop over block-cols on this node
         do iblk=my_first_blk,my_last_blk
            loc_iblk = loc_iblk + 1

            ! Loop over nonzero block-rows in this block-column
            do iidx=idxbuf(loc_iblk),idxbuf(loc_iblk+1)-1
               jblk = idxbuf(iidx)
               jelems = num_elems_on_atom(jblk,blks)

               ! Get this block and write it into buffer
               if (mat(1)%iscmplx) then
                  call sparse_get_block(zblk,mat(icomp),jblk,iblk)
                  do ielem=1,num_elems_on_atom(iblk,blks)
                     zmtxbuf(ptr:ptr+jelems-1) = zblk(1:jelems,ielem)
                     ptr = ptr + jelems
                  end do
               else
                  call sparse_get_block(dblk,mat(icomp),jblk,iblk)
                  do ielem=1,num_elems_on_atom(iblk,blks)
                     dmtxbuf(ptr:ptr+jelems-1) = dblk(1:jelems,ielem)
                     ptr = ptr + jelems
                  end do
               end if

            end do

         end do

      end do

      ! Record total length of data in each matrix
      datlen = (ptr - 1) / ncomps

      ! Deallocate temporary array for block
      if (mat(1)%iscmplx) then
         deallocate(zblk,stat=ierr)
         call utils_dealloc_check('internal_unsegment (sparse_write_vector)', &
              'zblk',ierr)
      else
         deallocate(dblk,stat=ierr)
         call utils_dealloc_check('internal_unsegment (sparse_write_vector)', &
              'dblk',ierr)
      end if

    end subroutine internal_unsegment

  end subroutine sparse_write_vector

  !============================================================================!
  ! This subroutine is a wrapper to read a single sparse matrix from a file by !
  ! calling sparse_read_vector with a single-element array.                    !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input) : The sparse matrix to be read                          !
  !   filename (input) : The filename to use                                   !
  !   unit     (input) : The (optional) unit number to use if file is open     !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, June 2009                                        !
  !============================================================================!

  subroutine sparse_read_scalar(mat,filename,unit)

    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat         ! Matrix to write
    character(len=*), intent(in) :: filename  ! Filename to use
    integer, optional, intent(in) :: unit     ! I/O unit

    ! Locals
    integer :: ierr
    type(SPAM3),allocatable :: mat_copy(:)

    ! Allocate and create a temporary SPAM3 array of length 1 to read into
    allocate(mat_copy(1),stat=ierr)
    call utils_alloc_check('sparse_read_scalar','mat_copy',ierr)
    call sparse_create(mat_copy(1),mat)

    ! Call vector version of sparse_read
    if (present(unit)) then
       call sparse_read_vector(mat_copy,filename,unit)
    else
       call sparse_read_vector(mat_copy,filename)
    end if

    ! Copy temporary array version into argument
    call sparse_copy(mat,mat_copy(1))

    ! Destroy and deallocate temporary SPAM3 array
    call sparse_destroy(mat_copy(1))
    deallocate(mat_copy,stat=ierr)
    call utils_dealloc_check('sparse_read_scalar','mat_copy',ierr)

  end subroutine sparse_read_scalar


  !============================================================================!
  ! This subroutine reads a vector of sparse matrices from a file.             !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat     (output) : The sparse matrix to be read                          !
  !   filename (input) : The filename                                          !
  !   unit     (input) : The (optional) unit number to use if file is open     !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, February 2007.                                    !
  ! Modification for SPAM3 by Nicholas Hine, June 2009                         !
  !============================================================================!

  subroutine sparse_read_vector(mat,filename,unit)
    use comms, only: comms_abort, comms_barrier, comms_bcast, comms_send, &
         comms_recv, pub_my_node_id, pub_on_root, pub_root_node_id, &
         pub_total_num_nodes
    use parallel_strategy, only: pub_orig_atom, pub_distr_atom, &
         pub_node_of_atom
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat(:)      ! Matrix to read
    character(len=*), intent(in) :: filename  ! Filename to use
    integer, optional, intent(in) :: unit     ! I/O unit

    ! Local variables
    integer :: ncomps               ! Number of components
    integer :: icomp                ! Loop counter for components
    integer :: ierr                 ! Error flag
    integer :: ilib                 ! Library entry for mat
    integer :: blks                 ! Identifier for blocking scheme
    integer :: iunit                ! Unit number for writing file
    integer :: idum                 ! Dummy variable
    integer :: file_nat             ! Number of atoms from file
    integer :: file_num             ! Number of NGWFs from file
    integer :: iblk                 ! Block loop counter
    integer :: iat_orig             ! Atom iblk in original order
    integer :: jblk                 ! Second block
    integer :: jat_orig             ! Atom jblk in original order
    integer :: iidx                 ! Index loop counter
    integer :: ielem,jelems         ! Element loop counter
    integer :: iptr                 ! Pointer
    integer :: idxlen               ! Index info length
    integer :: node                 ! Processor loop counter
    integer :: idxptr               ! Index pointer
    integer :: datptr               ! Data pointer
    integer :: idxptr0              ! Index pointer start for a node
    integer :: datptr0              ! Data pointer start for a node
    integer :: comp0                ! Offset between components in data
    integer :: blksize              ! Block size
    integer :: datlen               ! Buffer length for data
    integer :: size_idxinfo         ! Workaround for apparent bug in PGI F90
    integer, allocatable :: block_sizes(:)  ! Block sizes
    integer, allocatable :: ibuf(:)         ! Integer buffer
    real(kind=DP), allocatable :: dbuf(:)   ! Double precision real buffer
    complex(kind=DP), allocatable :: zbuf(:)! Complex buffer
    integer, allocatable :: idxbuf(:)       ! Comm buffer for index
    integer, allocatable :: idxinfo(:)      ! Index info for reading
    character(len=12) :: tag        ! Field tag
    character(len=12) :: cbuf       ! Character buffer
    real(kind=DP), allocatable :: dmtxbuf(:)    ! Comm buffer for real data
    complex(kind=DP), allocatable :: zmtxbuf(:) ! Comm buffer for complex data
    integer :: max_elems_on_atom                ! Maximum number of elements
    real (kind=DP), allocatable :: dblk(:,:)    ! Buffer for real block
    complex (kind=DP), allocatable :: zblk(:,:) ! Buffer for complex block

!CW
    integer :: ncomps_backup
!END CW

    ! Obtain number of components
    ncomps = size(mat)
    icomp = 1

!CW
    ncomps_backup=ncomps
!END CW

    ! Obtain library entry for mat
    ilib = mat(1)%lib
    blks = library(ilib)%row_blks

    ! Check that the matrix is square
    if (blks/=library(ilib)%col_blks) then
       write(stdout,'(a)') 'Error in sparse_read_vector: &
            &row and column blocking schemes do not match'
       call comms_abort
    end if

    ! Preliminaries on root node only
    if (pub_on_root) then

       if (present(unit)) then
          iunit = unit
       else
          iunit = utils_unit()
          open(unit=iunit,file=trim(filename),form='unformatted', &
               action='read',iostat=ierr,status='old')
          if (ierr /= 0) then
             write(stdout,'(3a,i6)') 'Error in sparse_read_vector: &
                  &opening file "',trim(filename),'" failed with code ',ierr
             call comms_abort
          end if
       end if

       ! Initialise some variables which should be set within the file
       file_nat = -1
       file_num = -1
    else
       iunit = -1 ! qoh: Initialise to avoid compiler warning
    end if

    ! Global preliminaries
    ! Allocate temporary array to store block
    max_elems_on_atom = maxval(num_elems_on_atom(:,blks))
    if (mat(1)%iscmplx) then
       allocate(zblk(max_elems_on_atom,max_elems_on_atom),stat=ierr)
       call utils_alloc_check('sparse_read_vector', &
            'zblk',ierr)
    else
       allocate(dblk(max_elems_on_atom,max_elems_on_atom),stat=ierr)
       call utils_alloc_check('sparse_read_vector', &
            'dblk',ierr)
    end if
    allocate(ibuf(2*pub_total_num_nodes),stat=ierr)
    call utils_alloc_check('sparse_read_vector','ibuf',ierr)

    do icomp=1,ncomps
       if (mat(1)%iscmplx) then
          mat(icomp)%zmtx = (0.0_DP,0.0_DP)
       else
          mat(icomp)%dmtx = 0.0_DP
       end if
    end do

    ! This is what the root node does: reads from file and sends data to
    ! nodes as appropriate
    if (pub_on_root) then

       size_idxinfo = 0 ! qoh: Initialise to prevent compiler warning
       ! Loop until end of file reached
       do

          ! Read tag and work out what to do with record
          read(iunit) tag

          select case (trim(tag(1:11)))

             ! File version information
          case ('VERSION')
             read(iunit) idum
             if (idum /= 1) then
                write(stdout,'(3a)') 'Error in sparse_read_vector: &
                     &malformed VERSION record in file "',trim(filename),'"'
                call comms_abort
             end if
             read(iunit) idum
             if (idum > nint(file_version*100)) then
                write(stdout,'(a,f5.2,3a,f5.2,a)') &
                     'WARNING in sparse_read_vector: version', &
                     real(idum,kind=DP)*0.01_DP,' of file "',trim(filename), &
                     '" is more recent than implemented version ', &
                     file_version,': continuing'
             end if

             ! Number of atoms
          case ('NUM_ATOMS')
             read(iunit) idum
             if (idum /= 1) then
                write(stdout,'(3a)') 'Error in sparse_read_vector: &
                     &malformed NUM_ATOMS record in file "',trim(filename),'"'
                call comms_abort
             end if
             read(iunit) file_nat
             if (file_nat /= library(ilib)%nblk) then
                write(stdout,'(3a,i6,a)') 'Error in sparse_read_vector: &
                     &atom number mismatch - file "',trim(filename), &
                     '" specifies',file_nat,' atoms'
                call comms_abort
             end if

             ! Number of NGWFs
          case ('NUM_NGWFS')
             read(iunit) idum
             if (idum /= 1) then
                write(stdout,'(3a)') 'Error in sparse_read_vector: &
                     &malformed NUM_NGWFS record in file "',trim(filename),'"'
                call comms_abort
             end if
             read(iunit) file_num
             if (file_num /= library(ilib)%nrows) then
                write(stdout,'(3a,i6,a)') 'Error in sparse_read_vector: &
                     &NGWF number mismatch - file "',trim(filename), &
                     '" specifies',file_num,' NGWFs'
                call comms_abort
             end if

             ! Data type of matrix
          case ('TYPE')
             read(iunit) idum
             if (idum /= 1) then
                write(stdout,'(3a)') 'Error in sparse_read_vector: &
                     &malformed TYPE record in file "',trim(filename),'"'
                call comms_abort
             end if
             read(iunit) idum
             if (idum == -1 .and. (.not. mat(1)%iscmplx)) then
                write(stdout,'(3a,i6,a)') 'Error in sparse_read_vector: &
                     &data type mismatch - file "',trim(filename), &
                     '" specifies complex data but expecting real data'
                call comms_abort
             else if (idum == 0 .and. mat(1)%iscmplx) then
                write(stdout,'(3a,i6,a)') 'Error in sparse_read_vector: &
                     &data type mismatch - file "',trim(filename), &
                     '" specifies real data but expecting complex data'
                call comms_abort
             end if

             ! Number of components
          case ('NUM_COMPS')
             read(iunit) idum
             if (idum /= 1) then
                write(stdout,'(3a)') 'Error in sparse_read_vector: &
                     &malformed NUM_COMPS record in file "',trim(filename),'"'
                call comms_abort
             end if
             read(iunit) idum
             if (idum /= ncomps) then
!CW
                   write(*,*) 'WARNING : components do not match, expecting : ', ncomps
                   write(*,*) '          we now copy the first component to the second one'
                   write(*,*) '          ncomps set to : ', idum
                   ncomps=idum
                   goto 44
             endif
!END CW
             if (idum /= ncomps) then
                write(stdout,'(3a,i3,a,i3)') 'Error in sparse_read_vector: &
                     &component number mismatch - file "',trim(filename), &
                     '" specifies',idum,' components but expecting',ncomps
                call comms_abort
             end if
!CW
             44 continue
!END CW
             ! Block sizes (number of NGWFs on each atom)
          case ('BLOCK_SIZES')
             read(iunit) idum
             if (file_nat < 0) file_nat = idum
             if (file_nat /= idum) then
                write(stdout,'(3a)') 'Error in sparse_read_vector: &
                     &inconsistent atom number data in file "',trim(filename), &
                     '"'
                call comms_abort
             end if
             if (file_nat /= library(ilib)%nblk) then
                write(stdout,'(3a,i6,a)') 'Error in sparse_read_vector: &
                     &atom number mismatch - file "',trim(filename), &
                     '" specifies',file_nat,' atoms'
                call comms_abort
             end if
             allocate(block_sizes(library(ilib)%nblk),stat=ierr)
             call utils_alloc_check('sparse_read_vector','block_sizes',ierr)
             read(iunit) block_sizes
             do iblk=1,library(ilib)%nblk
                iat_orig = pub_orig_atom(iblk)
                if (block_sizes(iat_orig) /= num_elems_on_atom(iblk,blks)) then
                   write(stdout,'(3a,i6,a)') 'Error in sparse_read_vector: &
                        &incompatible structure in file "',trim(filename),'"'
                   call comms_abort
                end if
             end do
             deallocate(block_sizes,stat=ierr)
             call utils_dealloc_check('sparse_read_vector','block_sizes',ierr)

             ! Block data
          case ('INDEX')

             ! Do we have all relevant info?
             if (file_nat < 0) then
                write(stdout,'(3a)') 'WARNING in sparse_read_vector: &
                     &file "',trim(filename), &
                     '" does not specify number of atoms: continuing'
                file_nat = library(ilib)%nblk
             end if
             if (file_num < 0) then
                write(stdout,'(3a)') 'WARNING in sparse_read_vector: &
                     &file "',trim(filename), &
                     '" does not specify number of NGWFs: continuing'
                file_num = library(ilib)%nrows
             end if

             ! Read in index information for data
             read(iunit) idxlen
             if (allocated(idxinfo)) then
                if (size_idxinfo < idxlen) then
                   deallocate(idxinfo,stat=ierr)
                   call utils_dealloc_check('sparse_read_vector','idxinfo',ierr)
                end if
             end if
             if (.not. allocated(idxinfo)) then
                allocate(idxinfo(idxlen),stat=ierr)
                call utils_alloc_check('sparse_read_vector','idxinfo',ierr)
                size_idxinfo = idxlen
             end if
             read(iunit) idxinfo(1:idxlen)

             ! Read in data which follows in next record
             read(iunit) cbuf
             if (cbuf(1:4) /= 'DATA') then
                write(stdout,'(3a)') 'Error in sparse_read_vector: &
                     &no DATA record following INDEX record in file "', &
                     trim(filename),'"'
                call comms_abort
             end if
             if ((cbuf(12:12) == 'D' .and. mat(1)%iscmplx) .or. &
                  (cbuf(12:12) == 'Z' .and. (.not. mat(1)%iscmplx))) then
                write(stdout,'(3a)') 'Error in sparse_read_vector: &
                     &incompatible data type in DATA block in file "', &
                     trim(filename),'"'
                call comms_abort
             end if
             read(iunit) datlen
             if (mod(datlen,ncomps) /= 0) then
                write(stdout,'(3a)') 'Error in sparse_read_vector: &
                     &length of DATA record in file "',trim(filename), &
                     '" is incompatible with number of components'
                call comms_abort
             end if
             datlen = datlen / ncomps
             if (mat(1)%iscmplx) then
                if (allocated(zbuf)) then
                   if (size(zbuf) < datlen*ncomps) then
                      deallocate(zbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_vector','zbuf',ierr)
                   end if
                end if
                if (.not. allocated(zbuf)) then
                   allocate(zbuf(datlen*ncomps),stat=ierr)
                   call utils_alloc_check('sparse_read_vector','zbuf',ierr)
                end if
                read(iunit) zbuf(1:datlen*ncomps)
             else
                if (allocated(dbuf)) then
                   if (size(dbuf) < datlen*ncomps) then
                      deallocate(dbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_vector','dbuf',ierr)
                   end if
                end if
                if (.not. allocated(dbuf)) then
                   allocate(dbuf(datlen*ncomps),stat=ierr)
                   call utils_alloc_check('sparse_read_vector','dbuf',ierr)
                end if
                read(iunit) dbuf(1:datlen*ncomps)
             end if

             ! Work out how much data must be sent to each node
             ibuf(1:2*pub_total_num_nodes) = 0
             do iidx=1,idxlen,3
                iat_orig = idxinfo(iidx+1)
                node = pub_node_of_atom(iat_orig)
                ibuf(2*node+1) = ibuf(2*node+1) + 3
                ibuf(2*node+2) = ibuf(2*node+2) + idxinfo(iidx+2)
             end do
!CW
             call comms_bcast(pub_root_node_id,ncomps)
!END CW
             call comms_bcast(pub_root_node_id,ibuf,2*pub_total_num_nodes)
             call comms_barrier

             ! Allocate buffers for communication
             if (allocated(idxbuf)) then
                if (size(idxbuf) < idxlen) then
                   deallocate(idxbuf,stat=ierr)
                   call utils_dealloc_check('sparse_read_vector','idxbuf',ierr)
                end if
             end if
             if (.not. allocated(idxbuf)) then
                allocate(idxbuf(idxlen),stat=ierr)
                call utils_alloc_check('sparse_read_vector','idxbuf',ierr)
             end if
             if (mat(1)%iscmplx) then
                if (allocated(zmtxbuf)) then
                   if (size(zmtxbuf) < datlen*ncomps) then
                      deallocate(zmtxbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_vector', &
                           'zmtxbuf',ierr)
                   end if
                end if
                if (.not. allocated(zmtxbuf)) then
                   allocate(zmtxbuf(datlen*ncomps),stat=ierr)
                   call utils_alloc_check('sparse_read_vector','zmtxbuf',ierr)
                end if
             else
                if (allocated(dmtxbuf)) then
                   if (size(dmtxbuf) < datlen*ncomps) then
                      deallocate(dmtxbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_vector', &
                           'dmtxbuf',ierr)
                   end if
                end if
                if (.not. allocated(dmtxbuf)) then
                   allocate(dmtxbuf(datlen*ncomps),stat=ierr)
                   call utils_alloc_check('sparse_read_vector','dmtxbuf',ierr)
                end if
             end if

             ! Deposit data which is stored on root node (no comms needed)
             iptr = 1
             do iidx=1,idxlen,3
                iat_orig = idxinfo(iidx+1)
                blksize = idxinfo(iidx+2)
                if (pub_node_of_atom(iat_orig) == pub_root_node_id) then
                   iblk = pub_distr_atom(iat_orig)
                   jat_orig = idxinfo(iidx)
                   jblk = pub_distr_atom(jat_orig)
                   jelems = num_elems_on_atom(jblk,blks)
                   comp0 = 0
                   do icomp=1,ncomps
                      if (mat(1)%iscmplx) then
                         do ielem=1,num_elems_on_atom(iblk,blks)
                            zblk(1:jelems,ielem) = zbuf(comp0+iptr+&
                                 (ielem-1)*jelems:comp0+iptr+ielem*jelems-1)
                         end do
                         call sparse_put_block(zblk,mat(icomp),jblk,iblk)
                      else
                         do ielem=1,num_elems_on_atom(iblk,blks)
                            dblk(1:jelems,ielem) = dbuf(comp0+iptr+&
                                 (ielem-1)*jelems:comp0+iptr+ielem*jelems-1)
                         end do
                         call sparse_put_block(dblk,mat(icomp),jblk,iblk)
                      end if
                      comp0 = comp0 + datlen
                   end do
                end if
                iptr = iptr + blksize
             end do

             ! Now deal with remaining nodes (comms required)
             idxptr = 1
             datptr = 1
             do node=1,pub_total_num_nodes-1
                if (ibuf(2*node+1) == 0) cycle
                iptr = 1
                idxptr0 = idxptr
                datptr0 = datptr
                do iidx=1,idxlen,3
                   iat_orig = idxinfo(iidx+1)
                   blksize = idxinfo(iidx+2)
                   if (node == pub_node_of_atom(iat_orig)) then
                      idxbuf(idxptr) = idxinfo(iidx)
                      idxbuf(idxptr+1) = iat_orig
                      idxbuf(idxptr+2) = blksize
                      comp0 = 0
                      do icomp=1,ncomps
                         if (mat(1)%iscmplx) then
                            zmtxbuf(datptr:datptr+blksize-1) = &
                                 zbuf(comp0+iptr:comp0+iptr+blksize-1)
                         else
                            dmtxbuf(datptr:datptr+blksize-1) = &
                                 dbuf(comp0+iptr:comp0+iptr+blksize-1)
                         end if
                         comp0 = comp0 + datlen
                         datptr = datptr + blksize
                      end do
                      idxptr = idxptr + 3
                   end if
                   iptr = iptr + blksize
                end do
                call comms_send(node,idxbuf(idxptr0),idxptr-idxptr0)
                if (mat(1)%iscmplx) then
                   call comms_send(node,zmtxbuf(datptr0),datptr-datptr0)
                else
                   call comms_send(node,dmtxbuf(datptr0),datptr-datptr0)
                end if
             end do

             ! End of file
          case ('ENDFILE')
!CW
             call comms_bcast(pub_root_node_id,ncomps)
!END CW
             ibuf(1:2*pub_total_num_nodes) = -1
             call comms_bcast(pub_root_node_id,ibuf,2*pub_total_num_nodes)
             exit

             ! Unknown tag
          case default
             if (tag(12:12) == 'I' .or. tag(12:12) == 'D' .or. &
                  tag(12:12) == 'Z') then
                write(stdout,'(5a)') 'WARNING in sparse_read_vector: &
                     &unknown record tag ',trim(tag(1:11)),' found in file "', &
                     trim(filename),'": continuing'
             else
                write(stdout,'(5a)') 'Error in sparse_read_vector: &
                     &malformed record tag ',trim(tag(1:11)), &
                     ' found in file "',trim(filename),'"'
                call comms_abort
             end if
             read(iunit) idum
             if (idum <= 0) then
                write(stdout,'(5a)') 'Error in sparse_read_vector: &
                     &malformed record length ',trim(tag(1:11)), &
                     ' found in file "',trim(filename),'"'
                call comms_abort
             end if
             if (tag(12:12) == 'I') then
                if (allocated(ibuf)) then
                   if (size(ibuf) < idum) then
                      deallocate(ibuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_vector','ibuf',ierr)
                   end if
                end if
                if (.not. allocated(ibuf)) then
                   allocate(ibuf(idum),stat=ierr)
                   call utils_alloc_check('sparse_read_vector','ibuf',ierr)
                end if
                read(iunit) ibuf(1:idum)
             else if (tag(12:12) == 'D') then
                if (allocated(dbuf)) then
                   if (size(dbuf) < idum) then
                      deallocate(dbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_vector','dbuf',ierr)
                   end if
                end if
                if (.not. allocated(dbuf)) then
                   allocate(dbuf(idum),stat=ierr)
                   call utils_alloc_check('sparse_read_vector','dbuf',ierr)
                end if
                read(iunit) dbuf(1:idum)
             else if (tag(12:12) == 'Z') then
                if (allocated(zbuf)) then
                   if (size(zbuf) < idum) then
                      deallocate(zbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_vector','zbuf',ierr)
                   end if
                end if
                if (.not. allocated(zbuf)) then
                   allocate(zbuf(idum),stat=ierr)
                   call utils_alloc_check('sparse_read_vector','zbuf',ierr)
                end if
                read(iunit) zbuf(1:idum)
             end if
          end select
       end do

       ! Other nodes wait for information from root node
    else
       ! Loop until end-of-file signal is received
       do
!CW
          call comms_bcast(pub_root_node_id,ncomps)
!END CW
          call comms_bcast(pub_root_node_id,ibuf,2*pub_total_num_nodes)
          idxlen = ibuf(2*pub_my_node_id+1)
          datlen = ibuf(2*pub_my_node_id+2)
          if (idxlen == -1) exit  ! signal that end of file reached
          call comms_barrier
          if (idxlen > 0) then
             if (allocated(idxbuf)) then
                if (size(idxbuf) < idxlen) then
                   deallocate(idxbuf,stat=ierr)
                   call utils_dealloc_check('sparse_read_vector','idxbuf',ierr)
                end if
             end if
             if (.not. allocated(idxbuf)) then
                allocate(idxbuf(idxlen),stat=ierr)
                call utils_alloc_check('sparse_read_vector','idxbuf',ierr)
             end if
             call comms_recv(pub_root_node_id,idxbuf,idxlen)
             if (mat(1)%iscmplx) then
                if (allocated(zmtxbuf)) then
                   if (size(zmtxbuf) < datlen*ncomps) then
                      deallocate(zmtxbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_vector', &
                           'zmtxbuf',ierr)
                   end if
                end if
                if (.not. allocated(zmtxbuf)) then
                   allocate(zmtxbuf(datlen*ncomps),stat=ierr)
                   call utils_alloc_check('sparse_read_vector','zmtxbuf',ierr)
                end if
                call comms_recv(pub_root_node_id,zmtxbuf,datlen*ncomps)
             else
                if (allocated(dmtxbuf)) then
                   if (size(dmtxbuf) < datlen*ncomps) then
                      deallocate(dmtxbuf,stat=ierr)
                      call utils_dealloc_check('sparse_read_vector', &
                           'dmtxbuf',ierr)
                   end if
                end if
                if (.not. allocated(dmtxbuf)) then
                   allocate(dmtxbuf(datlen*ncomps),stat=ierr)
                   call utils_alloc_check('sparse_read_vector','dmtxbuf',ierr)
                end if
                call comms_recv(pub_root_node_id,dmtxbuf,datlen*ncomps)
             end if

             iptr = 1
             do iidx=1,idxlen,3
                jat_orig = idxbuf(iidx)
                jblk = pub_distr_atom(jat_orig)
                iat_orig = idxbuf(iidx+1)
                iblk = pub_distr_atom(iat_orig)
                blksize = idxbuf(iidx+2)
                jelems = num_elems_on_atom(jblk,blks)
                do icomp=1,ncomps
                   if (mat(1)%iscmplx) then
                      do ielem=1,num_elems_on_atom(iblk,blks)
                         zblk(1:jelems,ielem) = &
                              zmtxbuf(iptr+(ielem-1)*jelems:iptr+ielem*jelems-1)
                      end do
                      call sparse_put_block(zblk,mat(icomp),jblk,iblk)
                   else
                      do ielem=1,num_elems_on_atom(iblk,blks)
                         dblk(1:jelems,ielem) = &
                              dmtxbuf(iptr+(ielem-1)*jelems:iptr+ielem*jelems-1)
                      end do
                      call sparse_put_block(dblk,mat(icomp),jblk,iblk)
                   end if
                   iptr = iptr + blksize
                end do

             end do

          end if

       end do
    end if

    ! Close file for input
    if (pub_on_root .and. (.not. present(unit))) then
       close(unit=iunit,iostat=ierr)
       if (ierr /= 0) then
          write(stdout,'(3a,i6)') 'Error in sparse_read_vector: &
               &closing file "',trim(filename),'" failed with code ',ierr
          call comms_abort
       end if
    end if

    ! Re-sync before continuing
    call comms_barrier

    ! Global finalisation
    if (allocated(zmtxbuf)) then
       deallocate(zmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_read_vector','zmtxbuf',ierr)
    end if
    if (allocated(dmtxbuf)) then
       deallocate(dmtxbuf,stat=ierr)
       call utils_dealloc_check('sparse_read_vector','dmtxbuf',ierr)
    end if
    if (allocated(idxbuf)) then
       deallocate(idxbuf,stat=ierr)
       call utils_dealloc_check('sparse_read_vector','idxbuf',ierr)
    end if
    if (allocated(idxinfo)) then
       deallocate(idxinfo,stat=ierr)
       call utils_dealloc_check('sparse_read_vector','idxinfo',ierr)
    end if
    if (allocated(zbuf)) then
       deallocate(zbuf,stat=ierr)
       call utils_dealloc_check('sparse_read_vector','zbuf',ierr)
    end if
    if (allocated(dbuf)) then
       deallocate(dbuf,stat=ierr)
       call utils_dealloc_check('sparse_read_vector','dbuf',ierr)
    end if
    if (mat(1)%iscmplx) then
       deallocate(zblk,stat=ierr)
       call utils_dealloc_check('sparse_read_vector','zblk',ierr)
    else
       deallocate(dblk,stat=ierr)
       call utils_dealloc_check('sparse_read_vector','dblk',ierr)
    end if
    deallocate(ibuf,stat=ierr)
    call utils_dealloc_check('sparse_read_vector','ibuf',ierr)
!CW
    if(ncomps_backup/=ncomps)then
     do icomp=ncomps+1,ncomps_backup
      write(*,*)' COPY DATA FROM COMPONENT 1 TO : ', icomp
      call sparse_copy(mat(icomp),mat(1))
     enddo
    endif
!END CW
  end subroutine sparse_read_vector


  !============================================================================!
  ! This subroutine allocates and initialises the segment type array for a     !
  ! SPAM3 hierarchical-sparsity matrix                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   str      (inout) : The matrix structure to be allocated                  !
  !   my_nze   (inout) : Number of nonzero elements on node                    !
  !   seg_nzb  (in)    : Number of nonzero blocks in each segment on node      !
  !   seg_nze  (in)    : Number of nonzero elements in each segment on node    !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009                                         !
  !============================================================================!

  subroutine sparse_segments_alloc(str,my_nze,seg_nze,seg_nzb)

    use comms, only: comms_alltoall, pub_my_node_id, pub_total_num_nodes
    use parallel_strategy, only: pub_num_atoms_on_node
    use rundat, only: pub_dense_threshold
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(STRUC3), intent(inout) :: str ! The matrix structure
    integer, intent(inout) :: my_nze   ! Number of nonzero elements on this proc
    integer, intent(inout) :: seg_nze(0:pub_total_num_nodes-1)
    integer, intent(in) :: seg_nzb(0:pub_total_num_nodes-1)

    ! Local variables
    integer :: ierr                    ! Error flag
    integer :: node                    ! Node counter
    real(kind=DP) :: seg_density       ! Density of the current segment
    integer :: seg_blks                ! Total number of blocks in this segment

    ! Allocate the segment info array
    allocate(str%seg_info(3,0:pub_total_num_nodes),stat=ierr)
    call utils_alloc_check('sparse_segments_alloc','str%seg_info('// &
         trim(str%structure)//')',ierr)

    ! Decide on which segments are dense, sparse and blank
    do node=0,pub_total_num_nodes-1

       seg_blks = pub_num_atoms_on_node(node) * &
            pub_num_atoms_on_node(pub_my_node_id)

       ! Calculate the density of nonzero blocks in this segment on this node
       seg_density = real(seg_nzb(node),kind=DP)/real(seg_blks,kind=DP)

       ! Determine what type of segment this should be
       if (seg_density >= pub_dense_threshold) then

          ! Segment is dense
          str%seg_info(s_type,node) = SEG_DENSE

          ! Remove the number of counted nonzero elements from total
          ! and add on the total number in the segment
          my_nze = my_nze - seg_nze(node)
          seg_nze(node) = num_elems_on_node(pub_my_node_id,str%col_blks) * &
               num_elems_on_node(node,str%row_blks)
          my_nze = my_nze + seg_nze(node)

       else if (seg_density > 0.0_DP) then

          ! Segment sparse-indexed
          str%seg_info(s_type,node) = SEG_SPARSE

       else

          ! Segment contains no elements so can always be skipped
          str%seg_info(s_type,node) = SEG_BLANK

       end if

    end do

    ! Zero final seg_type (unused)
    str%seg_info(s_type,pub_total_num_nodes) = 0

    ! Allocate list of lengths of local node's segments of index on other nodes
    allocate(str%idx_seg_lens(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('sparse_segments_alloc','str%idx_seg_lens('// &
         trim(str%structure)//')',ierr)

    ! Allocate list of lengths of local node's segments of data on other nodes
    allocate(str%mtx_seg_lens(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('sparse_segments_alloc','str%mtx_seg_lens('// &
         trim(str%structure)//')',ierr)

    ! Find lengths of other node's segments on local nodes
    do node=0,pub_total_num_nodes-1
       if (str%seg_info(s_type,node)==SEG_BLANK) then
          str%idx_seg_lens(node) = 0
       else
          str%idx_seg_lens(node) = pub_num_atoms_on_node(pub_my_node_id) &
               + 1 + seg_nzb(node)
       end if
       str%mtx_seg_lens(node) = seg_nze(node)
    end do

    ! Share with other nodes (only need to retain other nodes' segment lengths
    ! corresponding to this node)
    call comms_alltoall(str%idx_seg_lens,1)
    call comms_alltoall(str%mtx_seg_lens,1)

  end subroutine sparse_segments_alloc

  !============================================================================!
  ! This subroutine allocates the internal arrays of a block sparse matrix     !
  ! structure.                                                                 !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   str     (inout) : The matrix structure to be allocated                   !
  !   nze     (input) : Number of nonzero elements                             !
  !   nzb     (input) : Number of nonzero blocks                               !
  !   my_nze  (input) : Number of nonzero elements on this processor           !
  !   my_nzb  (input) : Number of nonzero blocks on this processor             !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Revised for distributed data, June 2006.                                   !
  ! Revised for SPAM3 by Nicholas Hine, May 2009.                              !
  !============================================================================!

  subroutine sparse_struc_alloc(str,nze,nzb,my_nze,my_nzb)

    use comms, only: comms_abort, comms_reduce, &
         pub_my_node_id, pub_total_num_nodes
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(STRUC3), intent(inout) :: str ! The matrix structure
    integer(kind=LONG), intent(in) :: nze ! Number of nonzero elements
    integer, intent(in) :: nzb         ! Number of nonzero blocks
    integer, intent(in) :: my_nze      ! Number of nonzero elements on this proc
    integer, intent(in) :: my_nzb      ! Number of nonzero blocks on this proc

    ! Local variables
    integer :: ierr                    ! Error flag
    integer :: nindices                ! Number of non-blank segment indices
    integer :: my_idx_len              ! Index length on this node
    integer, allocatable :: ibuf(:)    ! Buffer for communication

    ! Store sizes
    str%nrows = sum(num_elems_on_node(:,str%row_blks))
    str%mcols = sum(num_elems_on_node(:,str%col_blks))
    str%nze = nze
    str%nblk = pub_cell%nat
    str%nzb = nzb
    str%my_nze = my_nze
    str%my_nzb = my_nzb
    str%my_nblks = my_last_blk - my_first_blk + 1

    nindices = pub_total_num_nodes - count(str%seg_info(s_type,:)==SEG_BLANK)

    my_idx_len = (str%my_nblks + 1)*nindices + str%my_nzb

    ! Ensure arrays are allocated to at least size 1
    my_idx_len = max(my_idx_len,1)

    ! Allocate list of index lengths
    allocate(str%idx_lens(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','str%idx_lens('// &
         trim(str%structure)//')',ierr)

    ! Allocate block index
    allocate(str%blk_idx(my_idx_len),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','str%blk_idx('// &
         trim(str%structure)//')',ierr)

    ! Allocate block pointers
    allocate(str%blk_ptr(my_idx_len+1),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','str%blk_ptr('// &
         trim(str%structure)//')',ierr)

    ! Allocate buffer for communication
    allocate(ibuf(0:max(2,pub_total_num_nodes-1)),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','ibuf',ierr)

    ! Find maximum sizes
    ibuf = 0
    ibuf(0) = my_nze ; ibuf(1) = str%my_nblks ; ibuf(2) = my_nzb
    call comms_reduce('MAX',ibuf(0:2))
    str%max_nze = ibuf(0)
    str%max_nblks = ibuf(1)
    str%max_nzb = ibuf(2)

    ! Set up list of index lengths
    ibuf = 0
    ibuf(pub_my_node_id) = my_idx_len
    call comms_reduce('SUM',ibuf)
    str%idx_lens(0:pub_total_num_nodes-1) = ibuf(0:pub_total_num_nodes-1)

    ! Deallocate buffer for communication
    deallocate(ibuf,stat=ierr)
    call utils_dealloc_check('sparse_struc_alloc','ibuf',ierr)

  end subroutine sparse_struc_alloc

  !============================================================================!
  ! This subroutine copies the internal arrays of a block sparse matrix        !
  ! structure from one structure to another.                                   !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   str_dest (inout) : The matrix structure to be copied into                !
  !   str_src  (input) : Number of nonzero elements                            !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, April 2011.                                      !
  !============================================================================!

  subroutine sparse_struc_copy(str_dest,str_src)

    use comms, only: comms_abort, comms_reduce, &
         pub_my_node_id, pub_total_num_nodes
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(STRUC3), intent(inout) :: str_dest ! The matrix structure
    type(STRUC3), intent(in) :: str_src  ! The matrix structure

    ! Local variables
    integer :: ierr                    ! Error flag
    integer :: my_idx_len              ! Index length on this node

    ! Store basic info
    str_dest%structure = str_src%structure
    str_dest%nrows = str_src%nrows
    str_dest%mcols = str_src%mcols
    str_dest%nze = str_src%nze
    str_dest%nblk = str_src%nblk
    str_dest%nzb = str_src%nzb
    str_dest%my_nze = str_src%my_nze
    str_dest%my_nzb = str_src%my_nzb
    str_dest%my_nblks = str_src%my_nblks

    ! Allocate the segment info array
    allocate(str_dest%seg_info(3,0:pub_total_num_nodes),stat=ierr)
    call utils_alloc_check('sparse_struc_copy','str%seg_info('// &
         trim(str_dest%structure)//')',ierr)
    str_dest%seg_info(:,:) = str_src%seg_info(:,:)

    ! Allocate list of index lengths
    allocate(str_dest%idx_lens(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','str%idx_lens('// &
         trim(str_dest%structure)//')',ierr)
    str_dest%idx_lens = str_src%idx_lens

    ! Allocate block index
    my_idx_len = size(str_src%blk_idx)
    allocate(str_dest%blk_idx(my_idx_len),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','str%blk_idx('// &
         trim(str_dest%structure)//')',ierr)
    str_dest%blk_idx = str_src%blk_idx

    ! Allocate block pointers
    allocate(str_dest%blk_ptr(my_idx_len+1),stat=ierr)
    call utils_alloc_check('sparse_struc_alloc','str%blk_ptr('// &
         trim(str_dest%structure)//')',ierr)
    str_dest%blk_ptr = str_src%blk_ptr

    str_dest%max_nze = str_src%max_nze
    str_dest%max_nblks = str_src%max_nblks
    str_dest%max_nzb = str_src%max_nzb

  end subroutine sparse_struc_copy

  !============================================================================!
  ! This subroutine allocates the data arrays of a block sparse matrix.        !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat     (inout)           : The matrix structure to be allocated         !
  !   iscmplx (optional, input) : TRUE if data is complex, otherwise FALSE     !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Revised for distributed data, June 2006.                                   !
  !============================================================================!

  subroutine sparse_data_alloc(mat,iscmplx)

    use constants, only: real_size, cmplx_size
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat        ! The matrix structure
    logical, optional, intent(in) :: iscmplx ! Flag to indicate complex data

    ! Local variables
    integer :: ierr          ! Error flag
    logical :: loc_iscmplx   ! Local copy of iscmplx input flag
    integer :: my_nze        ! Length of data arrays for this matrix

    ! Make copy of optional input argument
    if (present(iscmplx)) then
       loc_iscmplx = iscmplx
    else
       loc_iscmplx = .false.
    end if

    ! Find length of data arrays from library
    my_nze = max(library(mat%lib)%my_nze,1)

    ! Allocate space for matrix elements
    if (loc_iscmplx) then

       allocate(mat%zmtx(my_nze),stat=ierr)
       call utils_alloc_check('sparse_data_alloc','mat%zmtx('// &
            trim(library(mat%lib)%structure)//')',ierr)
       mat%zmtx = (0.0_DP,0.0_DP)

    else

       allocate(mat%dmtx(my_nze),stat=ierr)
       call utils_alloc_check('sparse_data_alloc','mat%dmtx('// &
            trim(library(mat%lib)%structure)//')',ierr)
       mat%dmtx = 0.0_DP

    end if

    ! Fill in remaining details
    mat%iscmplx = loc_iscmplx
    mat%structure = library(mat%lib)%structure

    ! Record matrix memory usage
    if (mat%iscmplx) then
       local_mat_mem = local_mat_mem + int(library(mat%lib)%my_nze,kind=LONG)* &
            cmplx_size
       global_mat_mem = global_mat_mem + library(mat%lib)%nze*cmplx_size
    else
       local_mat_mem = local_mat_mem + int(library(mat%lib)%my_nze,kind=LONG)* &
            real_size
       global_mat_mem = global_mat_mem + library(mat%lib)%nze*real_size
    end if

  end subroutine sparse_data_alloc

  !============================================================================!
  ! This subroutine deallocates a block sparse matrix structure.               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   str     (inout) : The matrix structure to be deallocated                 !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                        !
  !============================================================================!

  subroutine sparse_segments_dealloc(str)

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(STRUC3), intent(inout) :: str        ! The matrix structure

    ! Local variables
    integer :: ierr          ! Error flag

    ! Deallocate list of lengths of local node's segment of data on other nodes
    if (allocated(str%mtx_seg_lens)) then
       deallocate(str%mtx_seg_lens,stat=ierr)
       call utils_dealloc_check('sparse_segments_dealloc','str%mtx_seg_lens(' &
            //trim(str%structure)//')',ierr)
    end if

    ! Deallocate list of local node's segment of index lengths on other nodes
    if (allocated(str%idx_seg_lens)) then
       deallocate(str%idx_seg_lens,stat=ierr)
       call utils_dealloc_check('sparse_segments_dealloc','str%idx_seg_lens(' &
            //trim(str%structure)//')',ierr)
    end if

    ! Deallocate segment information
    if (allocated(str%seg_info)) then
       deallocate(str%seg_info,stat=ierr)
       call utils_dealloc_check('sparse_segments_dealloc','str%seg_info(' &
            //trim(str%structure)//')',ierr)
    end if

  end subroutine sparse_segments_dealloc


  !============================================================================!
  ! This subroutine deallocates a block sparse matrix structure.               !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   str     (inout) : The matrix structure to be deallocated                 !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Revised for distributed data, June 2006.                                   !
  ! Adapted for SPAM3 by Nicholas Hine, May 2009.                              !
  !============================================================================!

  subroutine sparse_struc_dealloc(str)

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(STRUC3), intent(inout) :: str        ! The matrix structure

    ! Local variables
    integer :: ierr          ! Error flag

    ! Deallocate block pointers
    if (allocated(str%blk_ptr)) then
       deallocate(str%blk_ptr,stat=ierr)
       call utils_dealloc_check('sparse_struc_dealloc','str%blk_ptr('// &
            trim(str%structure)//')',ierr)
    end if

    ! Deallocate block index
    if (allocated(str%blk_idx)) then
       deallocate(str%blk_idx,stat=ierr)
       call utils_dealloc_check('sparse_struc_dealloc','str%blk_idx('// &
            trim(str%structure)//')',ierr)
    end if

    ! Deallocate list of index lengths
    if (allocated(str%idx_lens)) then
       deallocate(str%idx_lens,stat=ierr)
       call utils_dealloc_check('sparse_struc_dealloc','str%idx_lens('// &
            trim(str%structure)//')',ierr)
    end if

    ! Deallocate segment related arrays
    call sparse_segments_dealloc(str)

  end subroutine sparse_struc_dealloc


  !============================================================================!
  ! This subroutine deallocates the data arrays of a block sparse matrix.      !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat     (inout) : The matrix structure to be deallocated                 !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                       !
  ! Revised for distributed data, June 2006.                                   !
  !============================================================================!

  subroutine sparse_data_dealloc(mat)

    use constants, only: real_size, cmplx_size
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat        ! The matrix structure

    ! Local variables
    integer :: ierr          ! Error flag

    ! Deallocate data arrays
    if (allocated(mat%dmtx)) then
       deallocate(mat%dmtx,stat=ierr)
       call utils_dealloc_check('sparse_data_dealloc','mat%dmtx('// &
            trim(library(mat%lib)%structure)//')',ierr)
    end if

    if (allocated(mat%zmtx)) then
       deallocate(mat%zmtx,stat=ierr)
       call utils_dealloc_check('sparse_data_dealloc','mat%zmtx('// &
            trim(library(mat%lib)%structure)//')',ierr)
    end if

    ! Record matrix memory usage
    ! Record matrix memory usage
    if (mat%iscmplx) then
       local_mat_mem = local_mat_mem - int(library(mat%lib)%my_nze,kind=LONG)* &
            cmplx_size
       global_mat_mem = global_mat_mem - library(mat%lib)%nze*cmplx_size
    else
       local_mat_mem = local_mat_mem - int(library(mat%lib)%my_nze,kind=LONG)* &
            real_size
       global_mat_mem = global_mat_mem - library(mat%lib)%nze*real_size
    end if

  end subroutine sparse_data_dealloc


  !==========================================================================!
  ! This subroutine searches the library for an entry which matches the      !
  ! given structure identification string.                                   !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   name (inout) : Structure to search for in library                      !
  !   idx (output) : Entry number in library if found, otherwise 0           !
  !--------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2004.                                     !
  ! Moved to private module subroutine, Nicholas Hine, Dec 2007              !
  !==========================================================================!

  subroutine sparse_search_library(name,idx)

    implicit none

    ! Arguments
    character(len=*), intent(in) :: name     ! Structure identifier to search
    integer, intent(out) :: idx              ! Index in library if found

    ! Local variables
    integer :: ilib    ! Counter over library entries

    ! Set index to zero
    idx = 0

    ! Loop over library entries
    do ilib=1,num_library

       ! Compare structure identification strings
       if (trim(library(ilib)%structure) == trim(name)) then

          ! This library entry matches
          idx = ilib

          exit
       end if

    end do

  end subroutine sparse_search_library


  !============================================================================!
  ! This routine lists the library in the event that it fills up.              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   name   (input) : Identifying routine name                                !
  !----------------------------------------------------------------------------!
  ! Written by Peter Haynes, March 2006.                                       !
  !============================================================================!

  subroutine sparse_library_full(name)

    use comms, only: comms_abort, pub_on_root

    implicit none

    ! Argument
    character(len=*), intent(in) :: name  ! Identifying name

    ! Local variables
    integer :: ilib      ! Library entry

    ! Report error and library information
    if (pub_on_root) then
       write(stdout,'(3a)') 'Error in ',trim(name),': library full'
       write(stdout,'(a,i3)') '  Library size: ',max_library
       write(stdout,'(a,i3)') 'Current number: ',num_library
       write(stdout,'(a)') '  Library contents:'
       do ilib=1,max_library
          write(stdout,'(4x,i3,2a)') ilib,': ', &
               adjustl(library(ilib)%structure)
       end do
    end if

    ! And quit
    call comms_abort

  end subroutine sparse_library_full


  !============================================================================!
  ! This routine allocates a COM3 communicator. Takes multiple options for     !
  ! allocating different arrays, depending on which routine called it.         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com   (inout) : COM3 type with communicator arrays.                      !
  !   mat   (input) : SPAM3 matrix.                                            !
  !   nbuf  (input) : number of buffers to allocate.                           !
  !   alloc_mtx (in): whether to allocate the data parts of the communicator.  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, February 2012.                                   !
  !============================================================================!

  subroutine sparse_com_allocate(com,mat,nbuf,alloc_mtx,cropped,seg)

    use comms, only: pub_total_num_nodes, pub_null_handle
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    type(SPAM3), intent(in) :: mat
    integer, intent(in) :: nbuf
    logical, intent(in) :: alloc_mtx
    logical, intent(in) :: cropped
    logical, intent(in) :: seg

    ! Local Variables
    integer :: ierr

    ! Set up basic info
    com%lib = mat%lib
    com%iscmplx = mat%iscmplx
    com%cropped = cropped
    if (seg) then
       com%buflen = maxval(library(mat%lib)%idx_seg_lens)
       com%datlen = maxval(library(mat%lib)%mtx_seg_lens)
    else
       com%buflen = maxval(library(mat%lib)%idx_lens)
       com%datlen = library(mat%lib)%max_nze
    end if

    ! Allocate communication buffers
    allocate(com%seginfobuf(3,0:pub_total_num_nodes,1:nbuf),stat=ierr)
    call utils_alloc_check('sparse_com_allocate','com%seginfobuf',ierr)
    allocate(com%idxbuf(1:com%buflen,1:nbuf),stat=ierr)
    call utils_alloc_check('sparse_com_allocate','com%idxbuf',ierr)
    allocate(com%ptrbuf(1:com%buflen,1:nbuf),stat=ierr)
    call utils_alloc_check('sparse_com_allocate','com%ptrbuf',ierr)
    if (alloc_mtx) then
       if (com%iscmplx) then
          allocate(com%zmtxrecvbuf(com%datlen,1:nbuf),stat=ierr)
          call utils_alloc_check('sparse_com_allocate','com%zmtxrecvbuf',ierr)
          allocate(com%dmtxrecvbuf(1,2:2),stat=ierr)
          call utils_alloc_check('sparse_com_allocate','com%dmtxrecvbuf',ierr)
       else
          allocate(com%dmtxrecvbuf(com%datlen,1:nbuf),stat=ierr)
          call utils_alloc_check('sparse_com_allocate','com%dmtxrecvbuf',ierr)
          allocate(com%zmtxrecvbuf(1,2:2),stat=ierr)
          call utils_alloc_check('sparse_com_allocate','com%zmtxrecvbuf',ierr)
       end if
    end if

    ! Allocate comms buffers for cropped send/recv
    if (cropped) then
       allocate(com%ptrreqrecvbuf(4,1:com%buflen),stat=ierr)
       call utils_alloc_check('sparse_com_allocate','com%ptrreqrecvbuf',ierr)
       allocate(com%ptrreqsendbuf(4,1:com%buflen),stat=ierr)
       call utils_alloc_check('sparse_com_allocate','com%ptrreqsendbuf',ierr)
       if (alloc_mtx) then
          if (com%iscmplx) then
             allocate(com%zmtxsendbuf(com%datlen,1:nbuf),stat=ierr)
             call utils_alloc_check('sparse_com_allocate','com%zmtxsendbuf',ierr)
             allocate(com%dmtxsendbuf(1,2:2),stat=ierr)
             call utils_alloc_check('sparse_com_allocate','com%dmtxsendbuf',ierr)
          else
             allocate(com%dmtxsendbuf(com%datlen,1:nbuf),stat=ierr)
             call utils_alloc_check('sparse_com_allocate','com%dmtxsendbuf',ierr)
             allocate(com%zmtxsendbuf(1,2:2),stat=ierr)
             call utils_alloc_check('sparse_com_allocate','com%zmtxsendbuf',ierr)
          end if
       end if
    end if

    ! Allocate request and handle buffers
    com%num_handles = 2*pub_total_num_nodes + 6
    allocate(com%index_reqs(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('sparse_com_allocate','com%index_reqs',ierr)
    allocate(com%data_reqs(0:pub_total_num_nodes-1),stat=ierr)
    call utils_alloc_check('sparse_com_allocate','com%data_reqs',ierr)
    allocate(com%handles(1:com%num_handles),stat=ierr)
    call utils_alloc_check('sparse_com_allocate','com%handles',ierr)
    com%index_reqs(:) = index_sent
    com%data_reqs(:) = data_sent
    com%handles(:) = pub_null_handle

  end subroutine sparse_com_allocate


  !============================================================================!
  ! This routine deallocates a COM3 communicator.                              !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   com   (inout) : COM3 type with communicator arrays.                      !
  !   dealloc_mtx (in): whether to deallocate the data parts.                  !
  !----------------------------------------------------------------------------!
  ! Written by Nicholas Hine, February 2012.                                   !
  !============================================================================!

  subroutine sparse_com_deallocate(com,dealloc_mtx,cropped)

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    logical, intent(in) :: dealloc_mtx
    logical, intent(in) :: cropped

    ! Local Variables
    integer :: ierr

    ! Deallocate request and handle buffers
    deallocate(com%handles,stat=ierr)
    call utils_dealloc_check('sparse_com_deallocate','com%handles',ierr)
    deallocate(com%data_reqs,stat=ierr)
    call utils_dealloc_check('sparse_com_deallocate','com%data_reqs',ierr)
    deallocate(com%index_reqs,stat=ierr)
    call utils_dealloc_check('sparse_com_deallocate','com%index_reqs',ierr)

    ! Allocate comms buffers for cropped send/recv
    if (cropped) then
       if (dealloc_mtx) then
          deallocate(com%dmtxsendbuf,stat=ierr)
          call utils_dealloc_check('sparse_com_deallocate','com%dmtxsendbuf',ierr)
          deallocate(com%zmtxsendbuf,stat=ierr)
          call utils_dealloc_check('sparse_com_deallocate','com%zmtxsendbuf',ierr)
       end if
       deallocate(com%ptrreqsendbuf,stat=ierr)
       call utils_dealloc_check('sparse_com_deallocate','com%ptrreqsendbuf',ierr)
       deallocate(com%ptrreqrecvbuf,stat=ierr)
       call utils_dealloc_check('sparse_com_deallocate','com%ptrreqrecvbuf',ierr)
    end if

    ! Deallocate communication buffers
    if (dealloc_mtx) then
       deallocate(com%zmtxrecvbuf,stat=ierr)
       call utils_dealloc_check('sparse_com_deallocate','com%zmtxrecvbuf',ierr)
       deallocate(com%dmtxrecvbuf,stat=ierr)
       call utils_dealloc_check('sparse_com_deallocate','com%dmtxrecvbuf',ierr)
    end if

    deallocate(com%ptrbuf,stat=ierr)
    call utils_dealloc_check('sparse_com_deallocate','com%ptrbuf',ierr)
    deallocate(com%idxbuf,stat=ierr)
    call utils_dealloc_check('sparse_com_deallocate','com%idxbuf',ierr)
    deallocate(com%seginfobuf,stat=ierr)
    call utils_dealloc_check('sparse_com_deallocate','com%seginfobuf',ierr)

  end subroutine sparse_com_deallocate


  !==========================================================================!
  ! This subroutine sends either a segment of the index or the whole of the  !
  ! index of a given sparse matrix structure to another node                 !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   lib                                                                    !
  !   destnode                                                               !
  !   segonly                                                                !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                      !
  !==========================================================================!

  subroutine sparse_send_index(com,destnode,segonly,idx_handle,info_handle)

    use comms, only: comms_send, pub_my_node_id, pub_total_num_nodes, &
         pub_null_handle

    implicit none

    ! Arguments
    type(COM3), intent(in) :: com
    integer, intent(in) :: destnode
    logical, intent(in) :: segonly
    integer, intent(inout), optional :: idx_handle
    integer, intent(inout), optional :: info_handle

    ! Locals
    integer :: buflen, seg_start

    if (pub_my_node_id==destnode) return

    if (segonly) then
       seg_start = library(com%lib)%seg_info(s_idx,destnode)
       buflen = library(com%lib)%seg_info(s_idx,destnode+1) - seg_start
       if (buflen < 1) return
       if (present(idx_handle).and.present(info_handle)) then
          call comms_send(destnode,library(com%lib)%seg_info(:,destnode), &
               3,tag=SEGINFO_TAG,return_handle=info_handle)
          call comms_send(destnode,library(com%lib)%blk_idx(seg_start), &
               buflen,tag=BLKIDX_TAG,return_handle=idx_handle)
       else
          call comms_send(destnode,library(com%lib)%seg_info(:,destnode), &
               3,tag=SEGINFO_TAG)
          call comms_send(destnode,library(com%lib)%blk_idx(seg_start), &
               buflen,tag=BLKIDX_TAG)
       end if
    else
       buflen = library(com%lib)%idx_lens(pub_my_node_id)
       if (buflen < 1) then
          if (present(idx_handle).and.present(info_handle)) then
             idx_handle = pub_null_handle
             info_handle = pub_null_handle
          end if
          return
       end if
       if (present(idx_handle).and.present(info_handle)) then
          call comms_send(destnode,library(com%lib)%seg_info(:,:), &
               3*pub_total_num_nodes+3,tag=SEGINFO_TAG, &
               return_handle=info_handle)
          call comms_send(destnode,library(com%lib)%blk_idx, &
               buflen,tag=BLKIDX_TAG,return_handle=idx_handle)
       else ! synchronous receives - no need for handles
          call comms_send(destnode,library(com%lib)%seg_info(:,:), &
               3*pub_total_num_nodes+3,tag=SEGINFO_TAG)
          call comms_send(destnode,library(com%lib)%blk_idx, &
               buflen,tag=BLKIDX_TAG)
       end if
    end if

  end subroutine sparse_send_index


  !==========================================================================!
  ! This subroutine receives either the whole index or a segment of the      !
  ! index of a given sparse matrix structure from another node               !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   lib                                                                    !
  !   srcnode                                                                !
  !   segonly                                                                !
  !   idxbuf                                                                 !
  !   seginfobuf                                                             !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                      !
  !==========================================================================!

  subroutine sparse_recv_index(com,srcnode,ibuf,segonly,async)

    use comms, only: comms_recv, comms_irecv, pub_my_node_id, &
         pub_total_num_nodes
    use parallel_strategy, only: pub_num_atoms_on_node

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com      ! Sparse Matrix Communicator object
    integer, intent(in) :: srcnode        ! Source node from which to receive
    integer, intent(in) :: ibuf           ! Which buffer to use
    logical, intent(in) :: segonly        ! Segment only or whole matrix
    logical, intent(in) :: async          ! Whether to send asynchronously

    ! Locals
    integer :: buflen      ! Buffer Length
    integer :: seg_start   ! Segment start position

    ! If seg==.true., we are receiving the segment of the index only
    if (segonly) then

       buflen = library(com%lib)%idx_seg_lens(srcnode)
       if (buflen == 0) then
          com%seginfobuf(s_type,0,ibuf) = SEG_BLANK
          return
       end if
       if (pub_my_node_id /= srcnode) then
          if (async) then
             call comms_irecv(srcnode,com%seginfobuf(:,0,ibuf),3, &
                  tag=SEGINFO_TAG,handle=com%handles(com%info_hdl))
             call comms_irecv(srcnode,com%idxbuf(:,ibuf),buflen, &
                  tag=BLKIDX_TAG,handle=com%handles(com%idx_hdl))
          else
             call comms_recv(srcnode,com%seginfobuf(:,0,ibuf),3,tag=SEGINFO_TAG)
             call comms_recv(srcnode,com%idxbuf(:,ibuf),buflen,tag=BLKIDX_TAG)
          end if
       else
          seg_start = library(com%lib)%seg_info(s_idx,srcnode)
          com%seginfobuf(:,0,ibuf) = library(com%lib)%seg_info(:,srcnode)
          com%idxbuf(1:buflen,ibuf) = &
               library(com%lib)%blk_idx(seg_start:seg_start+buflen-1)
       end if

    else ! else we are receiving the whole index

       buflen = library(com%lib)%idx_lens(srcnode)
       if (buflen == 0) then
          com%seginfobuf(s_type,:,ibuf) = SEG_BLANK
          return
       end if
       if (pub_my_node_id /= srcnode) then
          if (async) then
             call comms_irecv(srcnode,com%seginfobuf(1,0,ibuf), &
                  3*pub_total_num_nodes+3,tag=SEGINFO_TAG, &
                  handle=com%handles(com%info_hdl))
             call comms_irecv(srcnode,com%idxbuf(1,ibuf),buflen, &
                  tag=BLKIDX_TAG,handle=com%handles(com%idx_hdl))
          else
             call comms_recv(srcnode,com%seginfobuf(1,0,ibuf), &
                  3*pub_total_num_nodes+3,tag=SEGINFO_TAG)
             call comms_recv(srcnode,com%idxbuf(1,ibuf),buflen,tag=BLKIDX_TAG)
          end if
       else
          com%seginfobuf(:,:,ibuf) = library(com%lib)%seg_info(:,:)
          com%idxbuf(1:buflen,ibuf) = library(com%lib)%blk_idx(1:buflen)
       end if

    end if

  end subroutine sparse_recv_index


  !==========================================================================!
  ! This subroutine sends either a segment of the pointers or the whole set  !
  ! of pointers for a given sparse matrix structure to another node          !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   lib                                                                    !
  !   destnode                                                               !
  !   segonly                                                                !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                      !
  !==========================================================================!

  subroutine sparse_send_pointers(com,destnode,segonly,ptr_handle)

    use comms, only: comms_send, pub_my_node_id

    implicit none

    ! Arguments
    type(COM3), intent(in) :: com     ! Matrix Communicator
    integer, intent(in) :: destnode   ! Destination node to which to send
    logical, intent(in) :: segonly    ! Segment only or whole matrix
    integer, intent(inout), optional :: ptr_handle

    ! Locals
    integer :: buflen, seg_start

    if (pub_my_node_id == destnode) return

    if (segonly) then
       seg_start = library(com%lib)%seg_info(s_idx,destnode)
       buflen = library(com%lib)%seg_info(s_idx,destnode+1) - seg_start
       call comms_send(destnode,library(com%lib)%blk_ptr(seg_start), &
            buflen,tag=BLKPTR_TAG)
    else
       buflen = library(com%lib)%idx_lens(pub_my_node_id)
       if (present(ptr_handle)) then
          call comms_send(destnode,library(com%lib)%blk_ptr,buflen, &
               tag=BLKPTR_TAG,return_handle=ptr_handle)
       else
          call comms_send(destnode,library(com%lib)%blk_ptr,buflen, &
               tag=BLKPTR_TAG)
       end if
    end if

  end subroutine sparse_send_pointers


  !==========================================================================!
  ! This subroutine receives either the whole set of pointers or a segment   !
  ! of the pointers for a given sparse matrix structure from another node.   !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   lib                                                                    !
  !   srcnode                                                                !
  !   segonly                                                                !
  !   segptrstart                                                            !
  !   ptrbuf                                                                 !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                      !
  !==========================================================================!

  subroutine sparse_recv_pointers(com,srcnode,ibuf,segonly,async)

    use comms, only: comms_recv, comms_irecv, pub_my_node_id
    use parallel_strategy, only: pub_num_atoms_on_node

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com      ! Sparse Matrix Communicator object
    integer, intent(in) :: srcnode        ! Source node from which to receive
    integer, intent(in) :: ibuf           ! Which buffer to use
    logical, intent(in) :: segonly        ! Segment only or whole matrix
    logical, intent(in) :: async          ! Whether to send asynchronously

    ! Locals
    integer :: buflen
    integer :: seg_start

    ! If segonly==.true. we are receiving this node's segment of the index only
    if (segonly) then
       buflen = library(com%lib)%idx_seg_lens(srcnode)
       seg_start = library(com%lib)%seg_info(s_idx,srcnode)
    else ! else we are receiving the whole set of pointers
       buflen = library(com%lib)%idx_lens(srcnode)
       seg_start = 1
    end if

    if (pub_my_node_id /= srcnode) then
       if (async) then
          call comms_irecv(srcnode,com%ptrbuf(:,ibuf),buflen,tag=BLKPTR_TAG, &
               handle=com%handles(com%ptr_hdl))
       else
          call comms_recv(srcnode,com%ptrbuf(:,ibuf),buflen,tag=BLKPTR_TAG)
       end if
    else
       com%ptrbuf(1:buflen,ibuf) = &
            library(com%lib)%blk_ptr(seg_start:seg_start+buflen-1)
    end if

  end subroutine sparse_recv_pointers


  !==========================================================================!
  ! This subroutine sends either the whole of the matrix data of a given     !
  ! node or the segment of it corresponding to the local node's rows         !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   mat                                                                    !
  !   destnode                                                               !
  !   segonly                                                                !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                      !
  !==========================================================================!

  subroutine sparse_send_data(mat,destnode,segonly)

    use comms, only: comms_send, pub_my_node_id, pub_total_num_nodes

    implicit none

    ! Arguments
    type(SPAM3),intent(in) :: mat
    integer, intent(in) :: destnode
    logical, intent(in) :: segonly

    ! Locals
    integer :: datlen
    integer :: seg_start

    if (pub_my_node_id == destnode) return

    ! If seg==.true., we are sending the segment of the data only
    if (segonly) then
       seg_start = library(mat%lib)%seg_info(s_ptr,destnode)
       datlen = library(mat%lib)%seg_info(s_ptr,destnode+1) - seg_start
    else
       seg_start = 1
       datlen = library(mat%lib)%seg_info(s_ptr,pub_total_num_nodes) - 1
    end if

    if (mat%iscmplx) then
       call comms_send(destnode,mat%zmtx(seg_start:seg_start+datlen-1), &
            datlen,tag=DATA_TAG)
    else
       call comms_send(destnode,mat%dmtx(seg_start:seg_start+datlen-1), &
            datlen,tag=DATA_TAG)
    end if

  end subroutine sparse_send_data


  !==========================================================================!
  ! This subroutine receives either the whole set of data or a segment       !
  ! of the data for a given sparse matrix structure from another node.       !
  !--------------------------------------------------------------------------!
  ! Arguments:                                                               !
  !   mat                                                                    !
  !   srcnode                                                                !
  !   segonly                                                                !
  !   dbuf                                                                   !
  !   zbuf                                                                   !
  !   totlen                                                                 !
  !--------------------------------------------------------------------------!
  ! Written by Nicholas Hine, May 2009.                                      !
  !==========================================================================!

  subroutine sparse_recv_data(com,mat,srcnode,ibuf,segonly,async)

    use comms, only: comms_recv, comms_irecv, pub_my_node_id, &
         pub_total_num_nodes

    implicit none

    ! Arguments
    type(COM3), intent(inout) :: com
    type(SPAM3),intent(in) :: mat
    integer, intent(in) :: srcnode
    integer, intent(in) :: ibuf
    logical, intent(in) :: segonly
    logical, intent(in) :: async          ! Whether to send asynchronously


    ! Locals
    integer :: seg_start
    integer :: datlen

    if (segonly) then
       datlen = library(mat%lib)%mtx_seg_lens(srcnode)
    else
       datlen = com%seginfobuf(s_ptr,pub_total_num_nodes,2)-1
    end if

    if (pub_my_node_id /= srcnode) then
       if (mat%iscmplx) then
          if (async) then
             call comms_irecv(srcnode,com%zmtxrecvbuf(:,ibuf),datlen, &
                  tag=DATA_TAG,handle=com%handles(com%data_hdl))
          else
             call comms_recv(srcnode,com%zmtxrecvbuf(:,ibuf),datlen, &
                  tag=DATA_TAG)
          end if
       else
          if (async) then
             call comms_irecv(srcnode,com%dmtxrecvbuf(:,ibuf),datlen, &
                  tag=DATA_TAG,handle=com%handles(com%data_hdl))
          else
             call comms_recv(srcnode,com%dmtxrecvbuf(:,ibuf),datlen, &
                  tag=DATA_TAG)
          end if
       end if
    else
       if (segonly) then
          seg_start = library(mat%lib)%seg_info(s_ptr,srcnode)
       else
          seg_start = 1
       end if
       if (mat%iscmplx) then
          com%zmtxrecvbuf(1:datlen,ibuf) = &
               mat%zmtx(seg_start:seg_start+datlen-1)
       else
          com%dmtxrecvbuf(1:datlen,ibuf) = &
               mat%dmtx(seg_start:seg_start+datlen-1)
       end if
    end if

  end subroutine sparse_recv_data


  !============================================================================!
  ! This subroutine converts a SPAM3 matrix into a column-indexed,             !
  ! non-segmented form                                                         !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (input)  : The sparse matrix to be converted                    !
  !   idxinfo  (output) : The block-index in a non-segmented form              !
  !   dmtx     (output) : The actual non-zero matrix elements                  !
  !----------------------------------------------------------------------------!
  ! Written by Simon M.-M. Dubois (January 2010)                               !
  !    * largely inspired by the sparse_write_vector subroutine written        !
  !      by Nicholas Hine                                                      !
  !============================================================================!

  subroutine sparse_convert_unsegment_real(mat,idxlen,idxinfo,dmtxlen, &
       dmtx,sizeflag)

    use comms, only: comms_abort, comms_barrier, comms_bcast, comms_send, &
         comms_recv, comms_reduce, pub_my_node_id, pub_on_root, &
         pub_root_node_id, pub_total_num_nodes
    use parallel_strategy, only: pub_orig_atom, pub_first_atom_on_node
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit

    implicit none

    ! Arguments
    type(SPAM3), intent(in)               :: mat(:)
    integer, intent(inout)                :: dmtxlen
    integer, intent(inout)                :: idxlen
    integer, intent(out)                  :: idxinfo(:)
    real(kind=DP), intent(out)            :: dmtx(:)
    logical                               :: sizeflag

    ! Local variables
    integer :: ncomps               ! Number of components
    integer :: icomp                ! Loop counter for components
    integer :: blks                 ! Identifier for blocking scheme
    integer :: my_nblks              ! Number of blocks stored on this node
    integer :: my_nzb               ! Number of non-zero blocks on this node
    integer :: ierr                 ! Error flag
    integer :: ilib                 ! Library entry for mat
    integer :: iblk                 ! Atom loop counter
    integer :: iat_orig             ! Atom iblk in original order
    integer :: jblk                 ! Second atom
    integer :: jat_orig             ! Atom jblk in original order
    integer :: loc_iblk             ! Local atom counter for iat
    integer :: iidx                 ! Index loop counter
    integer :: jidx                 ! Index info counter
    integer :: ielems               ! Number of elements on atom iat
    integer, allocatable :: idxbuf(:)       ! Comm buffer for index


    ! Check arguments
    if (mat(1)%iscmplx) then
       if (pub_on_root) write(stdout,'(a)') 'Error in &
            &sparse_convert_unsegment_real: real matrices only'
       call comms_abort
    end if

    ! Obtain number of components
    ncomps = size(mat)

    ! Obtain library entry for mat
    ilib = mat(1)%lib
    blks = library(ilib)%row_blks
    my_nblks = library(ilib)%my_nblks
    my_nzb = library(ilib)%my_nzb

    ! Check that the matrix is square
    if (blks/=library(ilib)%col_blks) then
       write(stdout,'(a)') 'Error in sparse_convert_unsegment_real: &
            &row and column blocking schemes do not match'
       call comms_abort
    end if

    ! Determine the space needed to store dmtx and idxinfo
    dmtxlen = library(ilib)%my_nze*ncomps
    idxlen  = my_nzb*3 + 2
    if (sizeflag) return

    ! Check that dmtx and idxinfo have the convenient dimension
    if (dmtxlen .ne. size(dmtx)) then
       write(stdout,'(a)') 'Error in sparse_convert_unsegment_real: &
            &storage space for data incorrectly allocated'
       call comms_abort
    end if
    if (idxlen .ne. size(idxinfo)) then
       write(stdout,'(a)') 'Error in sparse_convert_unsegment_real: &
            &storage space for index incorrectly allocated'
       call comms_abort
    end if

    ! Get the index for this node in non-segmented format
    allocate(idxbuf(my_nblks+my_nzb+2),stat=ierr)
    call utils_alloc_check('sparse_convert_unsegment_real','idxbuf',ierr)
    call sparse_generate_index(idxbuf,mat(1))

    ! Read segmented data into unsegmented buffers
    call internal_unsegment

    ! Write indexing information for the unsegment data into idxinfo :
    ! idxinfo(1)           > is set to idxlen
    ! idxinfo(2)           > is set to dmtxlen
    ! idxinfo(n*3+3:n*3+5) > contains the information regarding the nth
    !                        block of mat written in a column-indexed
    !                        non-segmented form

    idxinfo(1)=idxlen
    idxinfo(2)=dmtxlen


    jidx = 2
    loc_iblk = 0
    do iblk=my_first_blk,my_last_blk
       loc_iblk = loc_iblk + 1
       iat_orig = pub_orig_atom(iblk)
       ielems = num_elems_on_atom(iblk,blks)

       ! Loop over block-rows in column iat
       do iidx=idxbuf(loc_iblk),idxbuf(loc_iblk+1)-1
          jblk = idxbuf(iidx)
          jat_orig = pub_orig_atom(jblk)

          ! Generate index info: row, column and number of elements
          idxinfo(jidx+1) = jat_orig
          idxinfo(jidx+2) = iat_orig
          idxinfo(jidx+3) = ielems * num_elems_on_atom(jblk,blks)

          jidx = jidx + 3
       end do
    end do


    ! Re-sync nodes
    call comms_barrier

    ! Deallocate workspace
    deallocate(idxbuf,stat=ierr)
    call utils_dealloc_check('sparse_convert_unsegment_real','idxbuf',ierr)

  contains

    subroutine internal_unsegment

      implicit none

      ! Locals
      integer :: iidx,ptr
      integer :: iblk,loc_iblk
      integer :: jblk
      integer :: jelems
      integer :: ielem
      integer :: max_elems_on_atom
      real (kind=DP), allocatable :: dblk(:,:)

      ! Allocate temporary array to store block
      max_elems_on_atom = maxval(num_elems_on_atom(:,blks))
      allocate(dblk(max_elems_on_atom,max_elems_on_atom),stat=ierr)
      call utils_alloc_check('internal_unsegment  &
           &(sparse_convert_unsegment_real)', 'dblk',ierr)

      ptr = 1

      ! Loop over components of the matrix array
      do icomp=1,ncomps

         loc_iblk = 0

         ! Loop over block-cols on this node
         do iblk=my_first_blk,my_last_blk
            loc_iblk = loc_iblk + 1

            ! Loop over nonzero block-rows in this block-column
            do iidx=idxbuf(loc_iblk),idxbuf(loc_iblk+1)-1
               jblk = idxbuf(iidx)
               jelems = num_elems_on_atom(jblk,blks)

               ! Get this block and write it into buffer
               call sparse_get_block(dblk,mat(icomp),jblk,iblk)
               do ielem=1,num_elems_on_atom(iblk,blks)
                  dmtx(ptr:ptr+jelems-1) = dblk(1:jelems,ielem)
                  ptr = ptr + jelems
               end do

            end do

         end do

      end do

      ! Deallocate temporary array for block
      deallocate(dblk,stat=ierr)
      call utils_dealloc_check('internal_unsegment  &
           &(sparse_convert_unsegment)','dblk',ierr)

    end subroutine internal_unsegment

  end subroutine sparse_convert_unsegment_real


  !============================================================================!
  ! This subroutine converts a block sparse matrix stored in a column-indexed, !
  ! non-segmented form into a SPAM3                                            !
  !----------------------------------------------------------------------------!
  ! Arguments:                                                                 !
  !   mat      (output)  : The sparse matrix to be converted                   !
  !   idxinfo  (input)   : The block-index in a non-segmented form             !
  !   dmtx     (input)   : The actual non-zero matrix elements                 !
  !----------------------------------------------------------------------------!
  ! Written by Simon M.-M. Dubois (January 2010)                               !
  !    * largely inspired by the sparse_read_vector subroutine written         !
  !      by Nicholas Hine                                                      !
  !============================================================================!

  subroutine sparse_convert_segment_real(mat,idxinfo,dmtx)

    use comms, only: comms_alltoall, comms_abort, comms_barrier, &
         comms_bcast, comms_send, comms_recv, pub_my_node_id, &
         pub_on_root, pub_root_node_id, pub_total_num_nodes, &
         comms_irecv, comms_probe, comms_test, comms_waitany
    use parallel_strategy, only: pub_orig_atom, pub_distr_atom, &
         pub_node_of_atom
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_flush

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: mat(:)        ! Matrix to fill
    integer, intent(in)        :: idxinfo(:)
    real(kind=DP), intent(in)  :: dmtx(:)

    ! Local variables
    integer :: ncomps               ! Number of components
    integer :: icomp                ! Loop counter for components
    integer :: ierr                 ! Error flag
    integer :: ilib                 ! Library entry for mat
    integer :: blks                 ! Identifier for blocking scheme
    integer :: iblk                 ! Block loop counter
    integer :: iat_orig             ! Atom iblk in original order
    integer :: jblk                 ! Second block
    integer :: jat_orig             ! Atom jblk in original order
    integer :: iidx                 ! Index loop counter
    integer :: ielem,jelems         ! Element loop counter
    integer :: iptr                 ! Pointer
    integer :: idxlen               ! Index info length
    integer :: node_s               ! Processor loop counters
                                    ! (sender node in p2p commications)
    integer :: node_r               ! Processor loop counters
                                    ! (receiver node in p2p commications)
    integer :: idxptr               ! Index pointer
    integer :: datptr               ! Data pointer
    integer :: comp0                ! Offset between components in data
    integer :: blksize              ! Block size
    integer :: datlen               ! Buffer length for data
    integer :: max_elems_on_atom                ! Maximum number of elements
    integer :: send_num(2)
    integer :: receive_num(2)
    integer :: send_ptr(2)
    integer :: rcv_ptr(2)
    integer  :: ibuf_s(2*pub_total_num_nodes)           ! Integer buffers
    integer  :: ibuf_r(2*pub_total_num_nodes)           ! Integer buffers
    integer, allocatable :: idxbuf_s(:)         ! Comm buffers for index
    integer, allocatable :: idxbuf_r(:)         ! Comm buffers for index
    real(kind=DP), allocatable :: dmtxbuf_s(:)  ! Comm buffer for real data
    real(kind=DP), allocatable :: dmtxbuf_r(:)  ! Comm buffer for real data
    real (kind=DP), allocatable :: dblk(:,:)    ! Buffer for real block

    ! Obtain number of components
    ncomps = size(mat)
    icomp = 1

    ! Obtain library entry for mat
    ilib = mat(1)%lib
    blks = library(ilib)%row_blks

    ! Check that the matrix is square
    if (blks/=library(ilib)%col_blks) then
       write(stdout,'(a)') 'Error in sparse_convert_segment_real: &
            &row and column blocking schemes do not match'
       call comms_abort
    end if

    ! Global preliminaries
    ! Allocate temporary array to store block
    max_elems_on_atom = maxval(num_elems_on_atom(:,blks))
    allocate(dblk(max_elems_on_atom,max_elems_on_atom),stat=ierr)
    call utils_alloc_check('sparse_convert_segment_real', &
         'dblk',ierr)

    do icomp=1,ncomps
       mat(icomp)%dmtx = 0.0_DP
    end do

    idxlen = idxinfo(1)
    datlen = idxinfo(2)/ncomps

    ! Work out how much data must be sent to other nodes
    ibuf_r(:) = 0
    ibuf_s(:) = 0
    do iidx=3,idxlen,3
       iat_orig = idxinfo(iidx+1)
       node_r = pub_node_of_atom(iat_orig)
       if (node_r .ne. pub_my_node_id) then
          ibuf_s(2*node_r+1) =  ibuf_s(2*node_r+1) + 3
          ibuf_s(2*node_r+2) =  ibuf_s(2*node_r+2) +  idxinfo(iidx+2)
       endif
    enddo

    ibuf_r(:) = ibuf_s(:)             ! data to be sent to other nodes
    call comms_alltoall(ibuf_r,2)     ! data to be recieved by other nodes

    receive_num = 0
    send_num = 0
    do node_r=0,pub_total_num_nodes-1
       send_num(:) = send_num(:) + ibuf_s(2*node_r+1:2*node_r+2)
       receive_num(:) = receive_num(:) + ibuf_r(2*node_r+1:2*node_r+2)
    enddo

    ! Allocate buffers for data to be send
    if (send_num(1).ne.0) then
       allocate(idxbuf_s(send_num(1)),stat=ierr)
       call utils_alloc_check('sparse_convert_segment_real', &
            'idxbuf_s',ierr)
       allocate(dmtxbuf_s(send_num(2)*ncomps),stat=ierr)
       call utils_alloc_check('sparse_convert_segment_real', &
            'dmtxbuf_s',ierr)
       idxbuf_s = 0
       dmtxbuf_s = 0.0d0
    endif

    ! Allocate buffers for data to be received
    if (receive_num(1).ne.0) then
       allocate(idxbuf_r(receive_num(1)),stat=ierr)
       call utils_alloc_check('sparse_convert_segment_real', &
            'idxbuf_r',ierr)
       allocate(dmtxbuf_r(receive_num(2)*ncomps),stat=ierr)
       call utils_alloc_check('sparse_convert_segment_real', &
            'dmtxbuf_r',ierr)
       idxbuf_r = 0
       dmtxbuf_r = 0.0d0
    endif


    ! If required, fill the buffer arrays for outgoing communications
    if (send_num(1).ne.0) then

       idxptr = 1
       datptr = 1
       do node_r=0,pub_total_num_nodes-1

          if (ibuf_s(2*node_r+1) == 0) then
             cycle
          endif

          iptr = 1
          do iidx=3,idxlen,3
             iat_orig = idxinfo(iidx+1)
             blksize = idxinfo(iidx+2)
             if (node_r == pub_node_of_atom(iat_orig)) then
                idxbuf_s(idxptr) = idxinfo(iidx)
                idxbuf_s(idxptr+1) = iat_orig
                idxbuf_s(idxptr+2) = blksize
                comp0 = 0
                do icomp=1,ncomps
                   dmtxbuf_s(datptr:datptr+blksize-1) = &
                        dmtx(comp0+iptr:comp0+iptr+blksize-1)
                   comp0 = comp0 + datlen
                   datptr = datptr + blksize
                enddo
                idxptr = idxptr + 3
             endif
             iptr = iptr + blksize
          enddo

       enddo
    endif

    ! Synchronous loop over all nodes to send/receive data
    rcv_ptr(:) = 1
    do node_s=0,pub_total_num_nodes-1

       if (pub_my_node_id == node_s) then
          send_ptr(:) = 1
          do node_r=0,pub_total_num_nodes-1
             if (ibuf_s(2*node_r+1) == 0) then
                cycle
             endif
             call comms_send(node_r,idxbuf_s(send_ptr(1):send_ptr(1)+ibuf_s(2*node_r+1)-1), &
                        length=ibuf_s(2*node_r+1),tag=1)
             call comms_send(node_r,dmtxbuf_s(send_ptr(2):send_ptr(2)+ibuf_s(2*node_r+2)*ncomps-1), &
                        length=ibuf_s(2*node_r+2)*ncomps,tag=2)
             send_ptr(1) = send_ptr(1) + ibuf_s(2*node_r+1)
             send_ptr(2) = send_ptr(2) + ibuf_s(2*node_r+2)*ncomps
          enddo
       else
          if (ibuf_r(2*node_s+1) .ne. 0) then
             call comms_recv(node_s,idxbuf_r(rcv_ptr(1):rcv_ptr(1)+ibuf_r(2*node_s+1)-1), &
                        length=ibuf_r(2*node_s+1),tag=1)
             call comms_recv(node_s,dmtxbuf_r(rcv_ptr(2):rcv_ptr(2)+ibuf_r(2*node_s+2)*ncomps-1), &
                        length=ibuf_r(2*node_s+2)*ncomps,tag=2)
             rcv_ptr(1) = rcv_ptr(1) + ibuf_r(2*node_s+1)
             rcv_ptr(2) = rcv_ptr(2) + ibuf_r(2*node_s+2)*ncomps
          endif
       endif

       call comms_barrier
    enddo

    ! Deposit data that are stored locally
    iptr = 1
    do iidx=3,idxlen,3
       iat_orig = idxinfo(iidx+1)
       blksize = idxinfo(iidx+2)
       if (pub_node_of_atom(iat_orig) == pub_my_node_id) then

          iblk = pub_distr_atom(iat_orig)
          jat_orig = idxinfo(iidx)
          jblk = pub_distr_atom(jat_orig)
          jelems = num_elems_on_atom(jblk,blks)
          comp0 = 0
          do icomp=1,ncomps
             do ielem=1,num_elems_on_atom(iblk,blks)
                dblk(1:jelems,ielem) = dmtx(comp0+iptr+&
                     (ielem-1)*jelems:comp0+iptr+ielem*jelems-1)
             end do
             call sparse_put_block(dblk,mat(icomp),jblk,iblk)
             comp0 = comp0 + datlen
          end do
       end if
       iptr = iptr + blksize
    end do

    ! If required, empty the buffer arrays for incoming communications
    if (receive_num(1).ne.0) then

       iptr = 1
       do iidx=1,receive_num(1),3
          jat_orig = idxbuf_r(iidx)
          jblk = pub_distr_atom(jat_orig)
          iat_orig = idxbuf_r(iidx+1)
          iblk = pub_distr_atom(iat_orig)
          blksize = idxbuf_r(iidx+2)
          jelems = num_elems_on_atom(jblk,blks)
          do icomp=1,ncomps
             do ielem=1,num_elems_on_atom(iblk,blks)
                dblk(1:jelems,ielem) = &
                     dmtxbuf_r(iptr+(ielem-1)*jelems:iptr+ielem*jelems-1)
             end do
             call sparse_put_block(dblk,mat(icomp),jblk,iblk)
             iptr = iptr + blksize
          end do
       end do

    endif


    ! Deallocate communication buffers
    if (allocated(dmtxbuf_s)) then
       deallocate(dmtxbuf_s,stat=ierr)
       call utils_dealloc_check('sparse_convert_segment_real', &
            'dmtxbuf_s',ierr)
    end if
    if (allocated(dmtxbuf_r)) then
       deallocate(dmtxbuf_r,stat=ierr)
       call utils_dealloc_check('sparse_convert_segment_real', &
            'dmtxbuf_r',ierr)
    end if
    if (allocated(idxbuf_s)) then
       deallocate(idxbuf_s,stat=ierr)
       call utils_dealloc_check('sparse_convert_segment_real', &
            'idxbuf_s',ierr)
    end if
    if (allocated(idxbuf_r)) then
       deallocate(idxbuf_r,stat=ierr)
       call utils_dealloc_check('sparse_convert_segment_real', &
            'idxbuf_r',ierr)
    end if

    deallocate(dblk,stat=ierr)
    call utils_dealloc_check('sparse_convert_segment_real','dblk',ierr)

    ! Re-sync before continuing
    call comms_barrier

  end subroutine sparse_convert_segment_real


end module sparse

!==========================================================================!
! These subroutines override d/zgemm when nrows, ncols, nsum are all in    !
! the set {1,4,9}                                                          !
! They perform the operation c = c + a.b.                                  !
!--------------------------------------------------------------------------!
! Arguments:                                                               !
!   ncol,nrow : Dimensions of output matrix c,                             !
!                also number of rows of b and cols of a respectively       !
!   nsum      : Number of cols of a and rows of b                          !
!   a (input) : The block matrix a                                         !
!   a (input) : The block matrix b                                         !
!   c (inout) : The block matrix c                                         !
!--------------------------------------------------------------------------!
! This table shows the 27 types of block multiply that might be happening: !
! type nrow ncol nsum  matrix multiply                                     !
!    1    1    1    1  (1,1)*(1,1)->(1,1)                                  !
!    2    1    1    4  (1,4)*(4,1)->(1,1)                                  !
!    3    1    1    9  (1,9)*(9,1)->(1,1)                                  !
!    4    1    4    1  (1,1)*(1,4)->(1,4)                                  !
!    5    1    4    4  (1,4)*(4,4)->(1,4)                                  !
!    6    1    4    9  (1,9)*(9,4)->(1,4)                                  !
!    7    1    9    1  (1,1)*(1,9)->(1,9)                                  !
!    8    1    9    4  (1,4)*(4,9)->(1,9)                                  !
!    9    1    9    9  (1,9)*(9,9)->(1,9)                                  !
!   10    4    1    1  (4,1)*(1,1)->(4,1)                                  !
!   11    4    1    4  (4,4)*(4,1)->(4,1)                                  !
!   12    4    1    9  (4,9)*(9,1)->(4,1)                                  !
!   13    4    4    1  (4,1)*(1,4)->(4,4)                                  !
!   14    4    4    4  (4,4)*(4,4)->(4,4)                                  !
!   15    4    4    9  (4,9)*(9,4)->(4,4)                                  !
!   16    4    9    1  (4,1)*(1,9)->(4,9)                                  !
!   17    4    9    4  (4,4)*(4,9)->(4,9)                                  !
!   18    4    9    9  (4,9)*(9,9)->(4,9)                                  !
!   19    9    1    1  (9,1)*(1,1)->(9,1)                                  !
!   20    9    1    4  (9,4)*(4,1)->(9,1)                                  !
!   21    9    1    9  (9,9)*(9,1)->(9,1)                                  !
!   22    9    4    1  (9,1)*(1,4)->(9,4)                                  !
!   23    9    4    4  (9,4)*(4,4)->(9,4)                                  !
!   24    9    4    9  (9,9)*(9,4)->(9,4)                                  !
!   25    9    9    1  (9,1)*(1,9)->(9,9)                                  !
!   26    9    9    4  (9,4)*(4,9)->(9,9)                                  !
!   27    9    9    9  (9,9)*(9,9)->(9,9)                                  !
!--------------------------------------------------------------------------!
! Written by Nicholas Hine, November 2008                                  !
!==========================================================================!
subroutine sparse_dgemm_149(nrow,ncol,nsum,a,lda,b,ldb,c,ldc)

  use constants, only: DP

  implicit none

  ! Arguments
  integer,intent(in) :: nrow,ncol,nsum
  integer,intent(in) :: lda,ldb,ldc
  real(kind=DP),intent(in) :: a(lda,9),b(ldb,9)
  real(kind=DP),intent(inout) :: c(ldc,9)

  ! Locals
  integer :: i,j

  if      (nrow==1.and.ncol==1.and.nsum==1) then ! 1
     c(1,1)=c(1,1)+a(1,1)*b(1,1)
  else if (nrow==1.and.ncol==1.and.nsum==4) then ! 2
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
  else if (nrow==1.and.ncol==1.and.nsum==9) then ! 3
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
  else if (nrow==1.and.ncol==4.and.nsum==1) then ! 4
     c(1,1:4)=c(1,1:4)+a(1,1)*b(1,1:4)
  else if (nrow==1.and.ncol==4.and.nsum==4) then ! 5
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
     c(1,2)=c(1,2)+sum(a(1,1:4)*b(1:4,2))
     c(1,3)=c(1,3)+sum(a(1,1:4)*b(1:4,3))
     c(1,4)=c(1,4)+sum(a(1,1:4)*b(1:4,4))
  else if (nrow==1.and.ncol==4.and.nsum==9) then ! 6
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
     c(1,2)=c(1,2)+sum(a(1,1:9)*b(1:9,2))
     c(1,3)=c(1,3)+sum(a(1,1:9)*b(1:9,3))
     c(1,4)=c(1,4)+sum(a(1,1:9)*b(1:9,4))
  else if (nrow==1.and.ncol==9.and.nsum==1) then ! 7
     c(1,1:9)=c(1,1:9)+a(1,1)*b(1,1:9)
  else if (nrow==4.and.ncol==1.and.nsum==1) then ! 10
     c(1:4,1)=c(1:4,1)+a(1:4,1)*b(1,1)
  else if (nrow==4.and.ncol==1.and.nsum==4) then ! 11
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
     c(2,1)=c(2,1)+sum(a(2,1:4)*b(1:4,1))
     c(3,1)=c(3,1)+sum(a(3,1:4)*b(1:4,1))
     c(4,1)=c(4,1)+sum(a(4,1:4)*b(1:4,1))
  else if (nrow==4.and.ncol==1.and.nsum==9) then ! 12
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
     c(2,1)=c(2,1)+sum(a(2,1:9)*b(1:9,1))
     c(3,1)=c(3,1)+sum(a(3,1:9)*b(1:9,1))
     c(4,1)=c(4,1)+sum(a(4,1:9)*b(1:9,1))
  else if (nrow==4.and.ncol==4.and.nsum==1) then ! 13
     c(1:4,1)=c(1:4,1)+a(1:4,1)*b(1,1)
     c(1:4,2)=c(1:4,2)+a(1:4,1)*b(1,2)
     c(1:4,3)=c(1:4,3)+a(1:4,1)*b(1,3)
     c(1:4,4)=c(1:4,4)+a(1:4,1)*b(1,4)
  else if (nrow==4.and.ncol==4.and.nsum==4) then ! 14
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
     c(2,1)=c(2,1)+sum(a(2,1:4)*b(1:4,1))
     c(3,1)=c(3,1)+sum(a(3,1:4)*b(1:4,1))
     c(4,1)=c(4,1)+sum(a(4,1:4)*b(1:4,1))
     c(1,2)=c(1,2)+sum(a(1,1:4)*b(1:4,2))
     c(2,2)=c(2,2)+sum(a(2,1:4)*b(1:4,2))
     c(3,2)=c(3,2)+sum(a(3,1:4)*b(1:4,2))
     c(4,2)=c(4,2)+sum(a(4,1:4)*b(1:4,2))
     c(1,3)=c(1,3)+sum(a(1,1:4)*b(1:4,3))
     c(2,3)=c(2,3)+sum(a(2,1:4)*b(1:4,3))
     c(3,3)=c(3,3)+sum(a(3,1:4)*b(1:4,3))
     c(4,3)=c(4,3)+sum(a(4,1:4)*b(1:4,3))
     c(1,4)=c(1,4)+sum(a(1,1:4)*b(1:4,4))
     c(2,4)=c(2,4)+sum(a(2,1:4)*b(1:4,4))
     c(3,4)=c(3,4)+sum(a(3,1:4)*b(1:4,4))
     c(4,4)=c(4,4)+sum(a(4,1:4)*b(1:4,4))
  else if (nrow==4.and.ncol==4.and.nsum==9) then ! 15
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
     c(2,1)=c(2,1)+sum(a(2,1:9)*b(1:9,1))
     c(3,1)=c(3,1)+sum(a(3,1:9)*b(1:9,1))
     c(4,1)=c(4,1)+sum(a(4,1:9)*b(1:9,1))
     c(1,2)=c(1,2)+sum(a(1,1:9)*b(1:9,2))
     c(2,2)=c(2,2)+sum(a(2,1:9)*b(1:9,2))
     c(3,2)=c(3,2)+sum(a(3,1:9)*b(1:9,2))
     c(4,2)=c(4,2)+sum(a(4,1:9)*b(1:9,2))
     c(1,3)=c(1,3)+sum(a(1,1:9)*b(1:9,3))
     c(2,3)=c(2,3)+sum(a(2,1:9)*b(1:9,3))
     c(3,3)=c(3,3)+sum(a(3,1:9)*b(1:9,3))
     c(4,3)=c(4,3)+sum(a(4,1:9)*b(1:9,3))
     c(1,4)=c(1,4)+sum(a(1,1:9)*b(1:9,4))
     c(2,4)=c(2,4)+sum(a(2,1:9)*b(1:9,4))
     c(3,4)=c(3,4)+sum(a(3,1:9)*b(1:9,4))
     c(4,4)=c(4,4)+sum(a(4,1:9)*b(1:9,4))
  else if (nrow==9.and.ncol==1.and.nsum==1) then ! 19
     c(1:9,1)=c(1:9,1)+a(1:9,1)*b(1,1)
  else !  any of the rest
     do i=1,ncol
        do j=1,nrow
           c(j,i) = c(j,i) + sum(a(j,1:nsum)*b(1:nsum,i))
        end do
     end do
  end if

end subroutine sparse_dgemm_149

subroutine sparse_zgemm_149(nrow,ncol,nsum,a,lda,b,ldb,c,ldc)

  use constants, only: DP

  implicit none

  ! Arguments
  integer,intent(in) :: nrow,ncol,nsum
  integer,intent(in) :: lda,ldb,ldc
  complex(kind=DP),intent(in) :: a(lda,9),b(ldb,9)
  complex(kind=DP),intent(inout) :: c(ldc,9)

  ! Locals
  integer :: i,j

  if      (nrow==1.and.ncol==1.and.nsum==1) then ! 1
     c(1,1)=c(1,1)+a(1,1)*b(1,1)
  else if (nrow==1.and.ncol==1.and.nsum==4) then ! 2
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
  else if (nrow==1.and.ncol==1.and.nsum==9) then ! 3
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
  else if (nrow==1.and.ncol==4.and.nsum==1) then ! 4
     c(1,1:4)=c(1,1:4)+a(1,1)*b(1,1:4)
  else if (nrow==1.and.ncol==4.and.nsum==4) then ! 5
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
     c(1,2)=c(1,2)+sum(a(1,1:4)*b(1:4,2))
     c(1,3)=c(1,3)+sum(a(1,1:4)*b(1:4,3))
     c(1,4)=c(1,4)+sum(a(1,1:4)*b(1:4,4))
  else if (nrow==1.and.ncol==4.and.nsum==9) then ! 6
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
     c(1,2)=c(1,2)+sum(a(1,1:9)*b(1:9,2))
     c(1,3)=c(1,3)+sum(a(1,1:9)*b(1:9,3))
     c(1,4)=c(1,4)+sum(a(1,1:9)*b(1:9,4))
  else if (nrow==1.and.ncol==9.and.nsum==1) then ! 7
     c(1,1:9)=c(1,1:9)+a(1,1)*b(1,1:9)
  else if (nrow==4.and.ncol==1.and.nsum==1) then ! 10
     c(1:4,1)=c(1:4,1)+a(1:4,1)*b(1,1)
  else if (nrow==4.and.ncol==1.and.nsum==4) then ! 11
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
     c(2,1)=c(2,1)+sum(a(2,1:4)*b(1:4,1))
     c(3,1)=c(3,1)+sum(a(3,1:4)*b(1:4,1))
     c(4,1)=c(4,1)+sum(a(4,1:4)*b(1:4,1))
  else if (nrow==4.and.ncol==1.and.nsum==9) then ! 12
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
     c(2,1)=c(2,1)+sum(a(2,1:9)*b(1:9,1))
     c(3,1)=c(3,1)+sum(a(3,1:9)*b(1:9,1))
     c(4,1)=c(4,1)+sum(a(4,1:9)*b(1:9,1))
  else if (nrow==4.and.ncol==4.and.nsum==1) then ! 13
     c(1:4,1)=c(1:4,1)+a(1:4,1)*b(1,1)
     c(1:4,2)=c(1:4,2)+a(1:4,1)*b(1,2)
     c(1:4,3)=c(1:4,3)+a(1:4,1)*b(1,3)
     c(1:4,4)=c(1:4,4)+a(1:4,1)*b(1,4)
  else if (nrow==4.and.ncol==4.and.nsum==4) then ! 14
     c(1,1)=c(1,1)+sum(a(1,1:4)*b(1:4,1))
     c(2,1)=c(2,1)+sum(a(2,1:4)*b(1:4,1))
     c(3,1)=c(3,1)+sum(a(3,1:4)*b(1:4,1))
     c(4,1)=c(4,1)+sum(a(4,1:4)*b(1:4,1))
     c(1,2)=c(1,2)+sum(a(1,1:4)*b(1:4,2))
     c(2,2)=c(2,2)+sum(a(2,1:4)*b(1:4,2))
     c(3,2)=c(3,2)+sum(a(3,1:4)*b(1:4,2))
     c(4,2)=c(4,2)+sum(a(4,1:4)*b(1:4,2))
     c(1,3)=c(1,3)+sum(a(1,1:4)*b(1:4,3))
     c(2,3)=c(2,3)+sum(a(2,1:4)*b(1:4,3))
     c(3,3)=c(3,3)+sum(a(3,1:4)*b(1:4,3))
     c(4,3)=c(4,3)+sum(a(4,1:4)*b(1:4,3))
     c(1,4)=c(1,4)+sum(a(1,1:4)*b(1:4,4))
     c(2,4)=c(2,4)+sum(a(2,1:4)*b(1:4,4))
     c(3,4)=c(3,4)+sum(a(3,1:4)*b(1:4,4))
     c(4,4)=c(4,4)+sum(a(4,1:4)*b(1:4,4))
  else if (nrow==4.and.ncol==4.and.nsum==9) then ! 15
     c(1,1)=c(1,1)+sum(a(1,1:9)*b(1:9,1))
     c(2,1)=c(2,1)+sum(a(2,1:9)*b(1:9,1))
     c(3,1)=c(3,1)+sum(a(3,1:9)*b(1:9,1))
     c(4,1)=c(4,1)+sum(a(4,1:9)*b(1:9,1))
     c(1,2)=c(1,2)+sum(a(1,1:9)*b(1:9,2))
     c(2,2)=c(2,2)+sum(a(2,1:9)*b(1:9,2))
     c(3,2)=c(3,2)+sum(a(3,1:9)*b(1:9,2))
     c(4,2)=c(4,2)+sum(a(4,1:9)*b(1:9,2))
     c(1,3)=c(1,3)+sum(a(1,1:9)*b(1:9,3))
     c(2,3)=c(2,3)+sum(a(2,1:9)*b(1:9,3))
     c(3,3)=c(3,3)+sum(a(3,1:9)*b(1:9,3))
     c(4,3)=c(4,3)+sum(a(4,1:9)*b(1:9,3))
     c(1,4)=c(1,4)+sum(a(1,1:9)*b(1:9,4))
     c(2,4)=c(2,4)+sum(a(2,1:9)*b(1:9,4))
     c(3,4)=c(3,4)+sum(a(3,1:9)*b(1:9,4))
     c(4,4)=c(4,4)+sum(a(4,1:9)*b(1:9,4))
  else if (nrow==9.and.ncol==1.and.nsum==1) then ! 19
     c(1:9,1)=c(1:9,1)+a(1:9,1)*b(1,1)
  else !  any of the rest
     do i=1,ncol
        do j=1,nrow
           c(j,i) = c(j,i) + sum(a(j,1:nsum)*b(1:nsum,i))
        end do
     end do
  end if

end subroutine sparse_zgemm_149

