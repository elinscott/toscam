! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!   The subroutines in this file were written by                              !
!                                                                             !
!   Chris-Kriton Skylaris, Arash A. Mostofi, Peter D. Haynes and              !
!   Nicholas D.M. Hine                                                        !
!                                                                             !
!   TCM Group, Cavendish laboratory                                           !
!   Madingley Road                                                            !
!   Cambridge CB3 0HE                                                         !
!   UK                                                                        !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module basis

  use constants, only: DP
  use geometry, only: point

  implicit none

  private

  ! This structure holds all the information for the localised fft grid
  ! (contained in the "tight box") that will be used for each NGWF.
  ! It also holds all the information necessary for extracting a NGWF
  ! function from the ppd-in-standard-grid representation and placing it in
  ! the localised fft grid and vice-versa.
  type FUNCTION_TIGHT_BOX

     ! the number of zero padding points in each lattice vector direction
     integer :: pad1,pad2,pad3

     ! the very first and very last ppd number in each lattice vector
     ! direction that contains NGWF values.
     integer :: start_ppds1,start_ppds2,start_ppds3
     integer :: finish_ppds1,finish_ppds2,finish_ppds3

     ! the very first and the very last NGWF value point on the
     ! very first and very last ppd with function values in each lattice
     ! vector direction.
     integer :: start_pts1,start_pts2,start_pts3
     integer :: finish_pts1,finish_pts2,finish_pts3

     ! the number of tight points in the localised tight box grid in each
     ! lattice vector direction.
     integer :: tight_pts1, tight_pts2, tight_pts3

  end type FUNCTION_TIGHT_BOX


  type SPHERE

     ! the cartesian coordinates of the atom to which sphere belongs in a.u.
     type(POINT)       :: centre

     ! the radius of the sphere in atomic units
     ! cks: For NGWFs with halos, note that this is equal to the
     ! cks: NGWF radius (without the halo) and therefore it is different
     ! cks: from element%radius which includes the halo
     real(kind=DP) :: radius

     ! the number of ppds belonging to the sphere
     integer :: n_ppds_sphere

     ! ppd_list(1,:) contains a list of the indices (in the global ppd
     ! counting scheme) of the ppds that belong to the sphere
     ! ppd_list(2,:) contains a list of integers from -13 to 13 for each
     ! ppd belonging to the sphere. This number shows if the contribution
     ! to the sphere comes from the simulation cell or from a periodic image
     ! of the sphere in one of the neighbouring simulation cells
     integer, allocatable, dimension(:,:) :: ppd_list

     ! This is an offset showing the starting index (for the function inside
     ! the sphere) in the large array where all the functions are stored.
     integer :: offset

  end type SPHERE

  public :: SPHERE, FUNCTION_TIGHT_BOX

  ! ndmh: Location routines for finding points/vectors etc in cells & boxes
  public :: basis_function_origin_wrt_tb
  public :: basis_func_centre_wrt_fftbox
  public :: basis_find_function_wrt_box
  public :: basis_ket_start_wrt_fftbox
  public :: basis_location_fb_wrt_box
  public :: basis_location_func_wrt_cell
  public :: basis_point_wrt_box
  public :: basis_start_of_box_wrt_cell
  public :: basis_box_origin_to_atom

  ! ndmh: Sphere initialisation routines
  public :: basis_initialise_sphere
  public :: basis_sphere_deallocate
  public :: basis_ppd_location
  public :: basis_copy_sphere

  ! ndmh: Function-in-box manipulation routines
  public :: basis_copy_function_to_box
  public :: basis_add_function_to_box
  public :: basis_extract_function_from_box
  public :: basis_multiply_function_by_box
  public :: basis_clean_function
  public :: basis_put_tightbox_in_fftbox
  public :: basis_copy_tightbox_to_fftbox
  public :: basis_phase_on_fftbox_recip

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_function_origin_wrt_tb(tb_orig1, tb_orig2, tb_orig3,  &
       orig_wrt_cell, tight_box)

    !=======================================================================!
    ! This subroutine returns the coordinates of the atom centre of an NGWF !
    ! from the origin of its tightbox expressed in terms of                 !
    ! (in general non-integer) numbers of standard grid points              !
    ! along each lattice vector direction.                                  !
    !-----------------------------------------------------------------------!
    !                                                                       !
    !-----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 12/4/2004                         !
    !=======================================================================!

    use constants, only: DP, two_pi
    use geometry, only : POINT, operator(.DOT.)
    use simulation_cell, only : pub_cell

    implicit none

    real(kind =DP), intent(out) :: tb_orig1
    real(kind =DP), intent(out) :: tb_orig2
    real(kind =DP), intent(out) :: tb_orig3
    type(POINT), intent(in) :: orig_wrt_cell
    type(FUNCTION_TIGHT_BOX), intent(in) :: tight_box


    ! cks: << local variables>>
    integer :: start1
    integer :: start2
    integer :: start3
    real(kind=DP) :: ff
    real(kind=DP) :: inv_two_pi

    inv_two_pi = 1.0_DP / TWO_PI

    ! cks: origin of tightbox wrt cell in terms of integer numbers of grid
    !      points

    start1 = (tight_box%start_ppds1-1)*pub_cell%n_pt1 + tight_box%start_pts1 - 1
    start2 = (tight_box%start_ppds2-1)*pub_cell%n_pt2 + tight_box%start_pts2 - 1
    start3 = (tight_box%start_ppds3-1)*pub_cell%n_pt3 + tight_box%start_pts3 - 1

    ! cks: origin of function wrt tightbox in terms of fractional numbers of
    !      grid points
    ff = (orig_wrt_cell .DOT. pub_cell%b1)
    tb_orig1 = ff*inv_two_pi*real(pub_cell%total_pt1, kind=DP) - &
         real(start1, kind=DP)

    ff = (orig_wrt_cell .DOT. pub_cell%b2)
    tb_orig2 = ff*inv_two_pi*real(pub_cell%total_pt2, kind=DP) - &
         real(start2, kind=DP)

    ff = (orig_wrt_cell .DOT. pub_cell%b3)
    tb_orig3 = ff*inv_two_pi*real(pub_cell%total_pt3, kind=DP) - &
         real(start3, kind=DP)

  end subroutine basis_function_origin_wrt_tb


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function basis_func_centre_wrt_fftbox(ngwf_centre, &
       ngwf_start1, ngwf_start2, ngwf_start3, &
       ngwf_cell_start1, ngwf_cell_start2, ngwf_cell_start3)

    !=============================================================!
    ! Returns a vector from the origin of an FFTbox to the centre !
    ! of a function.                                              !
    !-------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 2/5/2003.               !
    ! Modified by Quintin Hill on 10/09/2008 to use pub_cell      !
    !=============================================================!

    use geometry, only: point, operator(+), operator(*)
    use simulation_cell, only: pub_cell

    implicit none

    type(POINT) :: basis_func_centre_wrt_fftbox

    ! Arguments
    type(POINT), intent(in) :: ngwf_centre
    integer, intent(in) :: ngwf_start1, ngwf_start2, ngwf_start3
    integer, intent(in) :: ngwf_cell_start1, ngwf_cell_start2, ngwf_cell_start3

    ! cks: internal
    type(POINT) :: buffvec

    buffvec = ngwf_centre + &
         real((ngwf_start1 -ngwf_cell_start1), kind =DP)*&
         &pub_cell%d1*pub_cell%a1_unit + &
         real((ngwf_start2 -ngwf_cell_start2), kind =DP)*&
         &pub_cell%d2*pub_cell%a2_unit + &
         real((ngwf_start3 -ngwf_cell_start3), kind =DP)*&
         &pub_cell%d3*pub_cell%a3_unit

    basis_func_centre_wrt_fftbox = buffvec

  end function basis_func_centre_wrt_fftbox


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_find_function_wrt_box( &
       row_start1, row_start2, row_start3, &           ! output
       box_start1, box_start2, box_start3, tight_box)  ! input

    !============================================================!
    ! This subroutine returns the starting gridpoint, in each    !
    ! lattice vector direction, of the current function with     !
    ! respect to the start of the FFTbox.                        !
    !------------------------------------------------------------!
    ! Written by Arash A. Mostofi, September 2002.               !
    ! Modified by Chris-Kriton Skylaris on 11/08/2007 to work in !
    ! the case where the FFT box coincides with the simulation   !
    ! cell.                                                      !
    ! Moved to basis_mod by Nicholas Hine, 29/06/2009.           !
    !============================================================!

    use simulation_cell, only : pub_cell, pub_fftbox

    implicit none

    type(FUNCTION_TIGHT_BOX), intent(in) :: tight_box
    integer, intent(out) :: row_start1,row_start2,row_start3
    integer, intent(in)  :: box_start1,box_start2,box_start3


    if (pub_fftbox%coin1) then
       row_start1 = 1
    else
       row_start1 = (tight_box%start_ppds1 - 1)*pub_cell%n_pt1 &
            + tight_box%start_pts1 - box_start1 + 1
       if (row_start1 .lt. 1) then
          row_start1 = row_start1 + pub_cell%total_pt1
       elseif (row_start1 .gt. pub_cell%total_pt1) then
          row_start1 = row_start1 - pub_cell%total_pt1
       endif
    endif

    if (pub_fftbox%coin2) then
       row_start2 = 1
    else
       row_start2 = (tight_box%start_ppds2 - 1)*pub_cell%n_pt2 &
            + tight_box%start_pts2 - box_start2 + 1
       if (row_start2 .lt. 1) then
          row_start2 = row_start2 + pub_cell%total_pt2
       elseif (row_start2 .gt. pub_cell%total_pt2) then
          row_start2 = row_start2 - pub_cell%total_pt2
       endif
    endif

    if (pub_fftbox%coin3) then
       row_start3 = 1
    else
       row_start3 = (tight_box%start_ppds3 - 1)*pub_cell%n_pt3 &
            + tight_box%start_pts3 - box_start3 + 1
       if (row_start3 .lt. 1) then
          row_start3 = row_start3 + pub_cell%total_pt3
       elseif (row_start3 .gt. pub_cell%total_pt3) then
          row_start3 = row_start3 - pub_cell%total_pt3
       endif
    endif

  end subroutine basis_find_function_wrt_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_ket_start_wrt_fftbox(row_start1,row_start2,row_start3, &
       n1, n2, n3)

    !===============================================================!
    ! Returns the starting position of the tightbox of the          !
    ! 'ket' function with respect to the beginning of the FFT box.  !
    !---------------------------------------------------------------!
    ! This subroutine also works when the length of the FFT box     !
    ! coincides with that of the simulation cell along one of more  !
    ! lattice vector directions.                                    !
    !---------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 6/4/2007, based on the    !
    ! basis_location_fa_wrt_box subroutine by Arash Mostofi.        !
    !===============================================================!

    use simulation_cell, only: pub_fftbox
    implicit none

    ! Arguments
    integer, intent(out) :: row_start1
    integer, intent(out) :: row_start2
    integer, intent(out) :: row_start3
    integer, intent(in)  :: n1,n2,n3


    if (pub_fftbox%coin1) then
       row_start1 =1
    else
       row_start1 =int(n1/3) + 1
    endif


    if (pub_fftbox%coin2) then
       row_start2 =1
    else
       row_start2 = int(n2/3) + 1
    endif


    if (pub_fftbox%coin3) then
       row_start3 =1
    else
       row_start3 = int(n3/3) + 1
    endif


  end subroutine basis_ket_start_wrt_fftbox


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_location_fb_wrt_box( &
       col_start1,col_start2,col_start3, &
       row_start1,row_start2,row_start3, &
       row_cell_start1,row_cell_start2,row_cell_start3, &
       col_cell_start1,col_cell_start2,col_cell_start3, &
       cell_n_pt1,cell_n_pt2,cell_n_pt3)

    !===================================================!
    ! Returns position of 'col' function in the FFTbox. !
    !---------------------------------------------------!
    ! Written by Arash A. Mostofi, September 2002.      !
    !===================================================!

    implicit none

    ! Arguments
    integer, intent(out) :: col_start1,col_start2,col_start3
    integer, intent(in)  :: row_start1,row_start2,row_start3
    integer, intent(in)  :: row_cell_start1,row_cell_start2,row_cell_start3
    integer, intent(in)  :: col_cell_start1,col_cell_start2,col_cell_start3
    integer, intent(in)  :: cell_n_pt1,cell_n_pt2,cell_n_pt3

    ! for 1st dimension
    col_start1 = col_cell_start1 - row_cell_start1 + row_start1
    if (col_start1 < 1) then
       col_start1 = col_start1 + cell_n_pt1
    else if (col_start1 > cell_n_pt1) then
       col_start1 = col_start1 - cell_n_pt1
    end if

    ! for 2nd dimension
    col_start2 = col_cell_start2 - row_cell_start2 + row_start2
    if (col_start2 < 1) then
       col_start2 = col_start2 + cell_n_pt2
    else if (col_start2 > cell_n_pt2) then
       col_start2 = col_start2 - cell_n_pt2
    end if

    ! for 3rd dimension
    col_start3 = col_cell_start3 - row_cell_start3 + row_start3
    if (col_start3 < 1) then
       col_start3 = col_start3 + cell_n_pt3
    else if (col_start3 > cell_n_pt3) then
       col_start3 = col_start3 - cell_n_pt3
    end if

  end subroutine basis_location_fb_wrt_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_location_func_wrt_cell( &
       start1,start2,start3,tight_box)

    !==========================================================!
    ! Returns the position of a function wrt simulation cell.  !
    !----------------------------------------------------------!
    ! Written by Arash A. Mostofi, September 2002.             !
    !==========================================================!

    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    type(FUNCTION_TIGHT_BOX), intent(in) :: tight_box
    integer, intent(out)               :: start1
    integer, intent(out)               :: start2
    integer, intent(out)               :: start3

    start1 = (tight_box%start_ppds1 - 1)*pub_cell%n_pt1 + tight_box%start_pts1
    start2 = (tight_box%start_ppds2 - 1)*pub_cell%n_pt2 + tight_box%start_pts2
    start3 = (tight_box%start_ppds3 - 1)*pub_cell%n_pt3 + tight_box%start_pts3

  end subroutine basis_location_func_wrt_cell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_point_wrt_box(dg1, dg2, dg3, &
       box_centre_point, box_pts1, box_pts2, box_pts3)

    !===================================================================!
    ! This subroutine returns the distance along each lattice vector    !
    ! direction, in gridspacings, from the origin of the FFT box to     !
    ! a given point. The point is positioned as close to the middle of  !
    ! the box as possible                                               !
    !-------------------------------------------------------------------!
    ! Written by Arash A. Mostofi, September 2002, originally named     !
    ! pseudo_centre_of_proj_wrt_box.                                    !
    ! Updated by Arash A. Mostofi on  10/07/2003 for non-orthogonal     !
    ! lattice geometries.                                               !
    ! Modified by Chris-Kriton Skylaris on 11/08/2007 to work in the    !
    ! case where the FFTbox coincides with the simulation cell.         !
    ! Moved to basis_mod by Nicholas Hine, 29/06/2009                   !
    !===================================================================!

    use constants, only: PI
    use geometry, only: POINT, OPERATOR(.dot.)
    use simulation_cell, only : pub_cell, pub_fftbox

    implicit none

    type(POINT), intent(in)     :: box_centre_point
    integer, intent(in)         :: box_pts1,box_pts2,box_pts3
    real(kind=DP), intent(out)  :: dg1,dg2,dg3 !{D}istance in {G}rid-spacings

    ! internal declarations
    real(kind=DP), parameter :: inv_two_pi = 0.5_DP/PI
    real(kind=DP) :: dd,xx
    real(kind=DP) :: tp1,tp2,tp3

    tp1 = real(pub_cell%total_pt1,kind=DP)
    tp2 = real(pub_cell%total_pt2,kind=DP)
    tp3 = real(pub_cell%total_pt3,kind=DP)

    if (pub_fftbox%coin1) then
       dg1=tp1*inv_two_pi*(box_centre_point.DOT.pub_cell%b1)
    else
       dd=real(box_pts1/2-1,kind=DP)
       xx=MOD(tp1*inv_two_pi*(box_centre_point.DOT.pub_cell%b1),1.0_DP)
       dg1 = dd + xx
    endif

    if (pub_fftbox%coin2) then
       dg2=tp2*inv_two_pi*(box_centre_point.DOT.pub_cell%b2)
    else
       dd=real(box_pts2/2-1,kind=DP)
       xx=MOD(tp2*inv_two_pi*(box_centre_point.DOT.pub_cell%b2),1.0_DP)
       dg2 = dd + xx
    endif

    if (pub_fftbox%coin3) then
       dg3=tp3*inv_two_pi*(box_centre_point.DOT.pub_cell%b3)
    else
       dd=real(box_pts3/2-1,kind=DP)
       xx=MOD(tp3*inv_two_pi*(box_centre_point.DOT.pub_cell%b3),1.0_DP)
       dg3 = dd + xx
    endif

  end subroutine basis_point_wrt_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_start_of_box_wrt_cell( &
       box_start1, box_start2, box_start3, &         ! output
       box_centre_point, dg1, dg2, dg3)        ! input

    !=================================================================!
    ! This subroutine returns the starting gridpoint in each lattice  !
    ! vector direction of the FFTbox in the simulation cell.          !
    !-----------------------------------------------------------------!
    ! Written by Arash A. Mostofi, September 2002.                    !
    ! Updated by Arash A. Mostofi on  10/07/2003 for non-orthogonal   !
    ! lattice geometries.                                             !
    ! Moved to basis_mod by Nicholas Hine, 29/06/2009                 !
    !=================================================================!

    use constants, only: PI
    use geometry, only: POINT, OPERATOR(.dot.)
    use simulation_cell, only : pub_cell

    implicit none

    type(POINT), intent(in) :: box_centre_point
    integer, intent(out) :: box_start1,box_start2,box_start3
    ! !{D}istance in {G}rid-spacings of box centre point from origin
    ! of box in each *lattice vector direction*. Not necessarily integer values
    real(kind=DP), intent(in) :: dg1,dg2,dg3

    ! internal declarations
    real(kind=DP), parameter :: inv_two_pi = 0.5_DP/PI
    real(kind=DP) :: aa
    real(kind=DP) :: tp1,tp2,tp3

    tp1=real(pub_cell%total_pt1,kind=DP)
    tp2=real(pub_cell%total_pt2,kind=DP)
    tp3=real(pub_cell%total_pt3,kind=DP)

    aa=tp1*inv_two_pi*(box_centre_point.DOT.pub_cell%b1) - dg1
    box_start1 = NINT(aa) + 1

    aa=tp2*inv_two_pi*(box_centre_point.DOT.pub_cell%b2) - dg2
    box_start2 = NINT(aa) + 1

    aa=tp3*inv_two_pi*(box_centre_point.DOT.pub_cell%b3) - dg3
    box_start3 = NINT(aa) + 1

  end subroutine basis_start_of_box_wrt_cell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_box_origin_to_atom(box_to_atom,box_origin,atom_origin, &
       cell_start1,cell_start2,cell_start3,da1,da2,da3)

    use geometry, only: POINT, OPERATOR(.dot.), OPERATOR(+), OPERATOR(-), &
         OPERATOR(*), geometry_magnitude
    use simulation_cell, only: pub_cell

    ! Arguments
    type(POINT),intent(out) :: box_origin
    type(POINT),intent(out) :: box_to_atom
    type(POINT),intent(in) :: atom_origin
    integer,intent(in) :: cell_start1,cell_start2,cell_start3
    type(POINT),intent(in) :: da1, da2, da3

    ! Find vector to origin of box
    box_origin = (cell_start1-1)*da1 + (cell_start2-1)*da2 + &
         (cell_start3-1)*da3

    ! ndmh: vector from origin of box to centre of atom
    box_to_atom = atom_origin - box_origin

    ! ndmh: check components of this vector in terms of lattice vectors
    ! ndmh: are all positive. If not, the box has been looped back into the
    ! ndmh: cell, so shift the origin back accordingly
    if ((box_to_atom.DOT.pub_cell%b1) < 0.0_DP) then
       box_origin = box_origin - pub_cell%a1
       box_to_atom = box_to_atom + pub_cell%a1
    end if
    if ((box_to_atom.DOT.pub_cell%b2) < 0.0_DP) then
       box_origin = box_origin - pub_cell%a2
       box_to_atom = box_to_atom + pub_cell%a2
    end if
    if ((box_to_atom.DOT.pub_cell%b3) < 0.0_DP) then
       box_origin = box_origin - pub_cell%a3
       box_to_atom = box_to_atom + pub_cell%a3
    end if

  end subroutine basis_box_origin_to_atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_phase_on_fftbox_recip(f_out, f_in, &
       n1,n2,n3,ld1,ld2,len1,len2,len3)

    !====================================================================!
    ! Application of any desired phase factor to a quantity (f_in)       !
    ! represented on the reciprocal space coarse grid pair box.          !
    ! The phase factor is characterised by the three quantities          !
    ! len1, len2 and len3, which give the phase vector in terms of       !
    ! numbers of gridpoints (not necessarily integer) along each lattice !
    ! direction. n1, n2 and n3 can be even or odd.                       !
    !--------------------------------------------------------------------!
    ! Written by Arash A. Mostofi, August 2002.                          !
    ! Modified January 2003 - got rid of make_even.                      !
    ! Modified April 2003 - complex-to-complex FFTs throughout.          !
    ! Modified April 2008 and August 2008 by Nick Hine for speed         !
    !====================================================================!

    use constants, only: DP, two_pi, cmplx_i
    use timer, only: timer_clock

    implicit none

    ! Arguments
    integer, intent(in)   :: n1,n2,n3,ld1,ld2
    real(kind=DP), intent(in)  :: len1,len2,len3
    complex(kind=DP), intent(out)  :: f_out(ld1,ld2,n3)
    complex(kind=DP), intent(in)  :: f_in(ld1,ld2,n3)

    ! <<< local variables >>>
    complex(kind=DP) :: phase1,phase2,phase3
    complex(kind=DP) :: phase_neg1,phase_neg2,phase_neg3
    complex(kind=DP) :: z1(n1)
    complex(kind=DP) :: z2(n2)
    complex(kind=DP) :: z3(n3)
    complex(kind=DP) :: z2z3
    integer  :: i1,i2,i3
    real(kind=DP) :: r_n1, r_n2, r_n3

    ! Start timer

    r_n1 = real(n1,kind=DP)
    r_n2 = real(n2,kind=DP)
    r_n3 = real(n3,kind=DP)

    ! the elementary phase factors for the three directions in space
    phase1 = exp(two_pi * cmplx_i * len1 / r_n1)
    phase2 = exp(two_pi * cmplx_i * len2 / r_n2)
    phase3 = exp(two_pi * cmplx_i * len3 / r_n3)

    ! for the negative reciprocal lattice points:
    ! e^(-i*2*pi*len2) and e^(-i*2*pi*len3)
    phase_neg1 = exp(-two_pi * cmplx_i * len1)
    phase_neg2 = exp(-two_pi * cmplx_i * len2)
    phase_neg3 = exp(-two_pi * cmplx_i * len3)

    ! initialise phase factor arrays
    z1(1) = cmplx(1.0_DP,0.0_DP,kind=DP)
    do i1=2,n1/2+1
       z1(i1) = z1(i1-1)*phase1
    end do
    z1(n1/2+2) = z1(n1/2+1)*phase_neg1*phase1
    do i1=n1/2+3,n1
       z1(i1) = z1(i1-1)*phase1
    end do

    z2(1) = cmplx(1.0_DP,0.0_DP,kind=DP)
    do i2=2,n2/2+1
       z2(i2) = z2(i2-1)*phase2
    end do
    z2(n2/2+2) = z2(n2/2+1)*phase_neg2*phase2
    do i2=n2/2+3,n2
       z2(i2) = z2(i2-1)*phase2
    end do

    z3(1) = cmplx(1.0_DP,0.0_DP,kind=DP)
    do i3=2,n3/2+1
       z3(i3) = z3(i3-1)*phase3
    end do
    z3(n3/2+2) = z3(n3/2+1)*phase_neg3*phase3
    do i3=n3/2+3,n3
       z3(i3) = z3(i3-1)*phase3
    end do

    ! Now copy over the f_in array while multiplying by phase
    do i3=1,n3

       do i2=1,n2

          z2z3 = z2(i2) * z3(i3)

          f_out(1:n1,i2,i3) = z1(1:n1)*z2z3*f_in(1:n1,i2,i3)

       end do

    end do

    ! Stop timer

  end subroutine basis_phase_on_fftbox_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_copy_sphere(current_sphere,old_sphere,offset)

    !=======================================================================!
    ! This subroutine initialises sphere_a by copying its values from       !
    ! sphere_b                                                              !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! sphere_a       (output)                                               !
    ! sphere_b       (input)                                                !
    ! offset         (input)
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine, 16/09/2008                                  !
    !=======================================================================!

    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(SPHERE), intent(out) :: current_sphere ! sphere to be initialised
    type(SPHERE), intent(in) :: old_sphere      ! sphere to copy from
    integer, intent(in) :: offset               ! current offset of NGWF data

    ! Local Variables
    integer :: ierr             ! allocation error flag
    integer :: ippd             ! loop counter

    current_sphere%centre%x = old_sphere%centre%x
    current_sphere%centre%y = old_sphere%centre%y
    current_sphere%centre%z = old_sphere%centre%z
    current_sphere%radius = old_sphere%radius
    current_sphere%n_ppds_sphere = old_sphere%n_ppds_sphere
    current_sphere%offset = offset

    allocate(current_sphere%ppd_list(2,current_sphere%n_ppds_sphere),stat=ierr)
    call utils_alloc_check('basis_copy_sphere','current_sphere%ppd_list',ierr)

    do ippd=1,current_sphere%n_ppds_sphere
       current_sphere%ppd_list(:,ippd) = old_sphere%ppd_list(:,ippd)
    end do

  end subroutine basis_copy_sphere


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_initialise_sphere(current_sphere,centre,radius,&
       current_offset,ppd_list,ppd_loc)

    !=======================================================================!
    ! This subroutine initialises a single sphere for one function on       !
    ! one atom.                                                             !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! centre         (input)                                                !
    ! radius         (input)                                                !
    ! cell           (input)                                                !
    ! current_offset (input)                                                !
    ! current_sphere (output)                                               !
    !-----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.                             !
    ! Updated to work when the FFT box coincides with the simulation cell   !
    ! by Chris-Kriton Skylaris on 131/08/2007.                              !
    ! Modified by Nicholas Hine on 24/06/2009 to combine ppd_loc            !
    ! with ppd_list for MPI efficiency.                                     !
    ! Modified to use function basis_ppd_location rather than ppd_location  !
    ! array by Nicholas Hine, October 2009.                                 !
    !=======================================================================!

    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout, TWO_PI
    use geometry, only: point, operator(+), operator(*), geometry_distance, &
         operator(.DOT.), magnitude
    use simulation_cell, only: pub_cell, pub_fftbox
    use utils, only: utils_alloc_check,utils_dealloc_check

    implicit none

    ! Arguments
    type(SPHERE), intent(out) :: current_sphere ! sphere to be initialised
    type(POINT), intent(in) :: centre           ! centre of sphere
    real(kind=DP), intent(in) :: radius         ! radius of sphere
    integer, intent(in) :: current_offset       ! psinc-counting offset for spere (local proc)
    integer, intent(inout) :: ppd_list(pub_cell%n_ppds) ! temporary ppd index list
    integer, intent(inout) :: ppd_loc(pub_cell%n_ppds)  ! temporary ppd loc list

    ! cks: internal variables
    integer :: ierr             ! allocation error flag
    integer :: a3_step          ! ppd counter along lattice vector a3
    integer :: a2_step          ! ppd counter along lattice vector a2
    integer :: a1_step          ! ppd counter along lattice vector a1
    integer :: a1_neighbour     ! neighbour simcell counter along a1
    integer :: a2_neighbour     ! neighbour simcell counter along a2
    integer :: a3_neighbour     ! neighbour simcell counter along a3
    integer :: ppd              ! counter of all ppds in sim cell
    integer :: one_or_zero      ! ppd belongs to sphere flag
    type(POINT) :: periodic_centre  ! centre od sphere in neighbour of simc ell
    logical :: ppd_near_sphere  ! result of quick test to select ppds near sphere

    real(kind=DP) :: a1_centre  ! ppd start along a1 (fractional coords)
    real(kind=DP) :: a2_centre  ! ppd start along a2 (fractional coords)
    real(kind=DP) :: a3_centre  ! ppd start along a3 (fractional coords)
    real(kind=DP) :: a1_radius  ! sphere radius expressed in ppds along a1
    real(kind=DP) :: a2_radius  ! sphere radius expressed in ppds along a2
    real(kind=DP) :: a3_radius  ! sphere radius expressed in ppds along a3
    integer :: a1_start         ! ppd start along a1 (periodic)
    integer :: a2_start         ! ppd start along a2 (periodic)
    integer :: a3_start         ! ppd start along a3 (periodic)
    integer :: a1_end           ! ppd finish along a1 (periodic)
    integer :: a2_end           ! ppd finish along a2 (periodic)
    integer :: a3_end           ! ppd finish along a3 (periodic)
    integer :: a1_ppds, a2_ppds, a3_ppds ! shorthand variables
    integer :: ppd_count        ! counter for ppd_list and ppd_loc
    integer :: ppd_check        ! counter for avoiding double inclusion
    integer :: ppd_lowest
    logical :: ppd_already_found

    ! cks: initialise radius of sphere
    current_sphere%radius = radius

    ! cks: initialise the centre of the sphere
    current_sphere%centre = centre

    ! cks: initialise the offset
    current_sphere%offset = current_offset

    a1_ppds = pub_cell%n_ppds_a1
    a2_ppds = pub_cell%n_ppds_a2
    a3_ppds = pub_cell%n_ppds_a3

    ! find centre in ppd coordinates along a1
    a1_centre = (centre .dot. pub_cell%b1) * real(a1_ppds,kind=DP) / TWO_PI
    ! find radius projected onto a1 axis in terms of a1 ppd lengths
    a1_radius = radius * magnitude(pub_cell%b1) * real(a1_ppds,kind=DP) / TWO_PI
    ! find range of ppds to check along a1
    a1_start = floor(a1_centre - a1_radius + 1_DP)
    a1_end =   floor(a1_centre + a1_radius + 1_DP)
    ! ensure we are not checking last slab twice
    if ((a1_radius>=1_DP).and.(modulo(a1_end + a1_ppds - 1, a1_ppds) == &
         modulo(a1_start + a1_ppds - 1, a1_ppds))) a1_end = a1_end - 1

    ! find centre in ppd coordinates along a2
    a2_centre = (centre .dot. pub_cell%b2) * real(a2_ppds,kind=DP) / TWO_PI
    a2_radius = radius * magnitude(pub_cell%b2) * real(a2_ppds,kind=DP) / TWO_PI
    a2_start = floor(a2_centre - a2_radius + 1_DP)
    a2_end =   floor(a2_centre + a2_radius + 1_DP)
    if ((a2_radius>=1_DP).and.(modulo(a2_end + a2_ppds - 1, a2_ppds) == &
         modulo(a2_start + a2_ppds - 1, a2_ppds))) a2_end = a2_end - 1

    ! find centre in ppd coordinates along a3
    a3_centre = (centre .dot. pub_cell%b3) * real(a3_ppds,kind=DP) / TWO_PI
    a3_radius = radius * magnitude(pub_cell%b3) * real(a3_ppds,kind=DP) / TWO_PI
    a3_start = floor(a3_centre - a3_radius + 1_DP)
    a3_end =   floor(a3_centre + a3_radius + 1_DP)
    if ((a3_radius>=1_DP).and.(modulo(a3_end + a3_ppds - 1, a3_ppds) == &
         modulo(a3_start + a3_ppds - 1, a3_ppds))) a3_end = a3_end - 1

    ppd_count = 0
    ! cks: loop over all ppds in simulation cell first
    do a3_step=a3_start,a3_end
       do a2_step=a2_start,a2_end
          do a1_step=a1_start,a1_end

             ppd = (modulo(a3_step + a3_ppds - 1, a3_ppds)*a2_ppds &
                  + modulo(a2_step + a2_ppds - 1, a2_ppds))*a1_ppds &
                  + modulo(a1_step + a1_ppds - 1, a1_ppds) + 1

             do a1_neighbour = -1,1
                do a2_neighbour = -1,1
                   do a3_neighbour = -1,1

                      ! this is the centre of the NGWF_sphere in one of the
                      ! neighbouring simulation cells
                      periodic_centre = current_sphere%centre + &
                           real(a1_neighbour,kind=DP)*pub_cell%a1 + &
                           real(a2_neighbour,kind=DP)*pub_cell%a2 + &
                           real(a3_neighbour,kind=DP)*pub_cell%a3

                      ! consider if a ppd belongs to the current NGWF
                      ! function (or its periodic image in one of the
                      ! surrounding simulation cells) only if the distance of
                      ! one of its vertices from the NGWF centre is
                      ! smaller than the sum of the NGWF radius and the
                      ! length of the three ppd edges. This test selects, for
                      ! all the spheres, a number of ppds that scales linearly
                      ! with system size.
                      ppd_near_sphere = (geometry_distance(periodic_centre, &
                           basis_ppd_location(ppd)) <= current_sphere%radius + &
                           real(pub_cell%n_pt1,kind=DP)*pub_cell%d1 + &
                           real(pub_cell%n_pt2,kind=DP)*pub_cell%d2 + &
                           real(pub_cell%n_pt3,kind=DP)*pub_cell%d3)

                      if (ppd_near_sphere) then

                         one_or_zero = basis_ppd_belongs_to_sphere( &
                              periodic_centre,basis_ppd_location(ppd), &
                              current_sphere%radius)

                         if (one_or_zero == 1) then

                            ppd_already_found = .false.
                            do ppd_check = 1,ppd_count
                               if (ppd_list(ppd_check)==ppd) then
                                  ppd_already_found = .true.
                                  exit
                               end if
                            end do
                            if (ppd_already_found .and. &
                                 ((.not.pub_fftbox%coin1) .and. &
                                 (.not.pub_fftbox%coin2) .and. &
                                 (.not.pub_fftbox%coin3))) then
                               if (pub_on_root) then
                                  write(stdout,'(a)') &
                                       'Error in basis_initialise_sphere:'
                                  write(stdout,'(a,i6)')'  a NGWF sphere and &
                                       &its periodic image have values on ppd',&
                                       ppd
                               end if
                               call comms_abort
                            end if

                            if (ppd_already_found) then
                               ppd_loc(ppd_check) = (9*a1_neighbour+ &
                                    3*a2_neighbour+a3_neighbour) * &
                                    one_or_zero - 14*(1-one_or_zero)
                            end if
                            if (.not.ppd_already_found) then
                               ppd_count = ppd_count + 1
                               ppd_list(ppd_count) = ppd
                               ppd_loc(ppd_count) = (9*a1_neighbour+ &
                                    3*a2_neighbour+a3_neighbour) * &
                                    one_or_zero - 14*(1-one_or_zero)
                            end if

                         end if

                      end if

                   end do
                end do
             end do

          end do
       end do
    end do

    ! cks: now determine the other quantities relevant to the sphere
    current_sphere%n_ppds_sphere = ppd_count

    allocate(current_sphere%ppd_list(2,current_sphere%n_ppds_sphere),stat=ierr)
    call utils_alloc_check('basis_initialise_sphere', &
         'current_sphere%ppd_list',ierr)

    ! ndmh: sort the ppds in ascending order
    do ppd_check=1,ppd_count
       ppd_lowest = 1
       do ppd=1,ppd_count
          if (ppd_list(ppd)<ppd_list(ppd_lowest)) ppd_lowest=ppd
       end do
       current_sphere%ppd_list(1,ppd_check) = ppd_list(ppd_lowest)
       current_sphere%ppd_list(2,ppd_check) = ppd_loc(ppd_lowest)
       ppd_list(ppd_lowest) = pub_cell%n_ppds + 1
    end do

  end subroutine basis_initialise_sphere


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  integer function basis_ppd_belongs_to_sphere(centre,ppd_start, &
       NGWF_radius)

    !===================================================================!
    ! This function returns the value 1 if the current ppd belongs to   !
    ! the current sphere and zero otherwise.                            !
    !-------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.                         !
    ! Updated by Chris-Kriton Skylaris on 13/08/2007.                   !
    !===================================================================!

    use constants, only: DP, TWO_PI
    use geometry, only: point, operator(.DOT.), geometry_magnitude, &
         operator(-), local_displacement, operator(+), geometry_distance
    use simulation_cell, only :pub_cell

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: NGWF_radius
    TYPE(POINT), intent(in) :: centre,ppd_start

    ! cks: internal variable declarations
    real(kind=DP) :: first_plane_a1_a2, second_plane_a1_a2
    real(kind=DP) :: first_plane_a2_a3, second_plane_a2_a3
    real(kind=DP) :: first_plane_a3_a1, second_plane_a3_a1
    integer :: more_than_two, loc_1, loc_2, loc_3
    real(kind=DP) :: local_1, local_2, local_3, length
    TYPE(POINT) :: current_displacement, current_point

    ! ======================= SPHERE CENTRE INSIDE PPD? ======================

    basis_ppd_belongs_to_sphere = 0

    ! cks: find the distance of the sphere from the a1_a2 planes defined by
    !      boundary points of the current ppd
    first_plane_a1_a2 = ((ppd_start - centre) .DOT. pub_cell%b3) / &
         geometry_MAGNITUDE(pub_cell%b3)

    second_plane_a1_a2 = (((ppd_start - centre) .DOT. pub_cell%b3) + &
         TWO_PI*(real(pub_cell%n_pt3-1,kind=DP)*pub_cell%d3) / &
         geometry_MAGNITUDE(pub_cell%a3)) / geometry_MAGNITUDE(pub_cell%b3)

    ! cks: find the distance of the sphere from the a2_a3 planes defined
    !      by boundary points of the current ppd
    first_plane_a2_a3 = ((ppd_start - centre) .DOT. pub_cell%b1) / &
         geometry_MAGNITUDE(pub_cell%b1)

    second_plane_a2_a3 = (((ppd_start - centre) .DOT. pub_cell%b1) + &
         TWO_PI*(real(pub_cell%n_pt1-1,kind=DP)*pub_cell%d1) / &
         geometry_MAGNITUDE(pub_cell%a1)) / geometry_MAGNITUDE(pub_cell%b1)

    ! cks: find the distance of the sphere from the a3_a1 planes defined
    !      by boundary points of the current ppd
    first_plane_a3_a1 = ((ppd_start - centre) .DOT. pub_cell%b2) / &
         geometry_MAGNITUDE(pub_cell%b2)

    second_plane_a3_a1 = (((ppd_start - centre) .DOT. pub_cell%b2) + &
         TWO_PI*(real(pub_cell%n_pt2-1,kind=DP)*pub_cell%d2) / &
         geometry_MAGNITUDE(pub_cell%a2)) / geometry_MAGNITUDE(pub_cell%b2)

    ! cks: check if ppd contains the centre of sphere and therefore belongs
    !      to the sphere
    more_than_two = 0
    if (first_plane_a1_a2*second_plane_a1_a2 <= 0.0_DP) &
         more_than_two = more_than_two + 1
    if (first_plane_a2_a3*second_plane_a2_a3 <= 0.0_DP) &
         more_than_two = more_than_two + 1
    if (first_plane_a3_a1*second_plane_a3_a1 <= 0.0_DP) &
         more_than_two = more_than_two + 1

    if (more_than_two == 3) then
       basis_ppd_belongs_to_sphere = 1
       return
    end if


    ! =================== END SPHERE CENTRE INSIDE PPD? =======================


    !$$$$$$$$$$$$ CHECK IF PPD BORDER POINTS BELONG TO SPHERE $$$$$$$$$$$$$$$$$

    ! cks: reaching this point means that the centre of the sphere is
    !      outside the ppd

    ! cks: Loop over all the boundary grid points of the current ppd and
    !      examine if any of them belongs to the current sphere.
    do loc_3=0, pub_cell%n_pt3-1
       local_3=real(loc_3,kind=DP)*pub_cell%d3
       do loc_2=0,pub_cell%n_pt2-1
          local_2=real(loc_2,kind=DP)*pub_cell%d2

          current_displacement = local_displacement(pub_cell%a1_unit, &
               pub_cell%a2_unit,pub_cell%a3_unit,0.0_DP,local_2,local_3)

          current_point = ppd_start + current_displacement
          length = geometry_distance(current_point,centre)

          ! cks: test if the points at the plane at the beginning of the a1
          !      direction belong
          if (length <= NGWF_radius) then
             basis_ppd_belongs_to_sphere = 1
             return
          end if

          current_displacement = local_displacement(pub_cell%a1_unit, &
               pub_cell%a2_unit,pub_cell%a3_unit, &
               real(pub_cell%n_pt1-1,kind=DP)*pub_cell%d1, &
               local_2,local_3)

          current_point = ppd_start + current_displacement
          length = geometry_distance(current_point,centre)

          ! cks: test if the points at the plane at the end of the a1 direction
          !      belong
          if (length <= NGWF_radius) then
             basis_ppd_belongs_to_sphere = 1
             return
          end if

          if (loc_3 == 0 .or. loc_3 == pub_cell%n_pt3-1 .or. loc_2 == 0 .or. &
               loc_2 == pub_cell%n_pt2-1) then

             do loc_1=1, pub_cell%n_pt1-2
                local_1=real(loc_1,kind=DP)*pub_cell%d1

                current_displacement = local_displacement(pub_cell%a1_unit,&
                     pub_cell%a2_unit, pub_cell%a3_unit,local_1,local_2,local_3)

                current_point = ppd_start + current_displacement
                length = geometry_distance(current_point,centre)

                ! cks: test if the points at the planes at the end and at the
                !      beginning of the a2 and a3 directions belong
                if (length <= NGWF_radius) then
                   basis_ppd_belongs_to_sphere = 1
                   return
                end if

             end do

          end if
       end do
    end do

    !$$$$$$$$$$$$ END CHECK IF PPD BORDER POINTS BELONG TO SPHERE $$$$$$$$$$$$$

  end function basis_ppd_belongs_to_sphere


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_sphere_deallocate(current_sphere)

    !===================================================!
    ! Deallocate pointers in the sphere object.         !
    !---------------------------------------------------!
    ! Written by Victor Milman on 03/11/2006.           !
    ! Modified by Nicholas Hine on 24/06/2009 to        !
    ! combine ppd_loc with ppd_list for MPI efficiency. !
    !===================================================!

    use utils, only : utils_dealloc_check
    implicit none
    type(SPHERE), intent(inout) :: current_sphere

    integer :: ierr

    ! deallocate "sphere" pointers
    if (allocated(current_sphere%ppd_list)) then
       deallocate(current_sphere%ppd_list,stat=ierr)
       call utils_dealloc_check('basis_sphere_deallocate', &
            'current_sphere%ppd_list',ierr)
    end if

  end subroutine basis_sphere_deallocate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  type(POINT) function basis_ppd_location(ppd)

    !====================================================================!
    ! This function returns the cartesian coordinates of the origin of a !
    ! given ppd.                                                         !
    !--------------------------------------------------------------------!
    ! Written by Nicholas Hine on 19/10/2009 to replace the previous     !
    ! array of ppd locations as that was O(N) memory and not O(1/Nproc). !
    !====================================================================!

    use geometry, only: operator(*), operator(+)
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    integer, intent(in) :: ppd

    ! Local variables
    integer :: d1,d2,d3
    real(kind=DP) :: rd1, rd2, rd3

    ! ndmh: find integer coordinates of this ppd along each lattice direction
    d1 = modulo(ppd-1,pub_cell%n_ppds_a1)
    d2 = modulo((ppd-d1)/pub_cell%n_ppds_a1, pub_cell%n_ppds_a2)
    d3 = (ppd-d1-d2*pub_cell%n_ppds_a1)/(pub_cell%n_ppds_a2*pub_cell%n_ppds_a1)

    ! ndmh: find distance along each lattice direction of this ppd
    rd1 = real(d1,kind=DP)*real(pub_cell%n_pt1,kind=DP)*pub_cell%d1
    rd2 = real(d2,kind=DP)*real(pub_cell%n_pt2,kind=DP)*pub_cell%d2
    rd3 = real(d3,kind=DP)*real(pub_cell%n_pt3,kind=DP)*pub_cell%d3

    ! ndmh: find cartesian coordinates of this ppd in terms of unit vectors
    ! ndmh: along each lattice vector (avoiding geometry_mod function calls)
    basis_ppd_location%x = rd1*pub_cell%a1_unit%x + &
                           rd2*pub_cell%a2_unit%x + &
                           rd3*pub_cell%a3_unit%x
    basis_ppd_location%y = rd1*pub_cell%a1_unit%y + &
                           rd2*pub_cell%a2_unit%y + &
                           rd3*pub_cell%a3_unit%y
    basis_ppd_location%z = rd1*pub_cell%a1_unit%z + &
                           rd2*pub_cell%a2_unit%z + &
                           rd3*pub_cell%a3_unit%z

  end function basis_ppd_location


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_copy_function_to_box(fa_box,box_n1,box_n2,box_n3, & ! in-out
       offset1,offset2,offset3,fa_tightbox,fa_on_grid,fa_sphere)       ! input

    !=================================================================!
    ! This subroutine copies an NGWF stored in its ppd representation !
    ! to a three-dimensional array with the size of the function's    !
    ! tightbox.                                                       !
    !-----------------------------------------------------------------!
    ! Copied from basis_copy_function_to_fftbox by Nicholas Hine on   !
    ! 23/07/2009 to replace old double-loop version of same. See      !
    ! comments on above.                                              !
    ! Amalgamated all the copy_function_to_*box routines into one     !
    ! multi-purpose routine, Nicholas Hine, 28/02/2011.               !
    !=================================================================!

    use simulation_cell, only : pub_cell
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(FUNCTION_TIGHT_BOX), intent(in) :: fa_tightbox ! function tightbox
    integer, intent(in) :: box_n1,box_n2,box_n3
    integer, intent(in) :: offset1,offset2,offset3
    real(kind=DP), intent(out) :: fa_box(box_n1,box_n2,box_n3) ! box data
    real(kind=DP), intent(in) :: fa_on_grid(:) ! function data to copy from
    type(SPHERE), intent(in) :: fa_sphere ! sphere of function's PPDs

    ! Local Variables
    integer :: fa_ppd      ! PPD counter within sphere
    integer :: fa_point    ! point counter within fa_on_grid
    integer :: a1_pos      ! absolute 1-coordinate of current PPD
    integer :: a2_pos      ! absolute 2-coordinate of current PPD
    integer :: a3_pos      ! absolute 3-coordinate of current PPD
    integer :: start1      ! start of points within box along direction 1 of current PPD
    integer :: start2      ! start of points within box along direction 2 of current PPD
    integer :: start3      ! start of points within box along direction 3 of current PPD
    integer :: finish1     ! end of points within box along direction 1 of current PPD
    integer :: finish2     ! end of points within box along direction 2 of current PPD
    integer :: finish3     ! end of points within box along direction 3 of current PPD
    integer :: box_offset1 ! offset to express 1-start in box
    integer :: box_offset2 ! offset to express 2-start in box
    integer :: box_offset3 ! offset to express 3-start in box
    integer :: box_pt1     ! current point in box along direction-1
    integer :: box_pt2     ! current point in box along direction-2
    integer :: box_pt3     ! current point in box along direction-3
    integer :: box_start1  ! start point of box points along direction-1
    integer :: box_end1    ! end point of box points along direction-1
    integer :: point2      ! point counter wrt current PPD along direction-2
    integer :: point3      ! point counter wrt current PPD along direction-3

    call timer_clock('basis_copy_function_to_box',1)

    fa_box(:,:,:) = 0.0_DP

    if (pub_cell%n_pts == 1) then
       ! cks: case of only one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               pub_cell%n_ppds_a1, pub_cell%n_ppds_a2, pub_cell%n_ppds_a3)

          box_pt1 = a1_pos - fa_tightbox%start_ppds1 + offset1
          box_pt2 = a2_pos - fa_tightbox%start_ppds2 + offset2
          box_pt3 = a3_pos - fa_tightbox%start_ppds3 + offset3

          fa_point = fa_sphere%offset + fa_ppd -1

          ! ndmh: copy ngwf value to box
          fa_box(box_pt1, box_pt2, box_pt3) = fa_on_grid(fa_point)
       enddo

    else

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               pub_cell%n_ppds_a1, pub_cell%n_ppds_a2, pub_cell%n_ppds_a3)

          call basis_lims_1d_in_ppd_in_tight(start1, finish1, box_offset1, &
               a1_pos, fa_tightbox%start_ppds1, fa_tightbox%finish_ppds1, &
               fa_tightbox%start_pts1, fa_tightbox%finish_pts1, pub_cell%n_pt1)

          call basis_lims_1d_in_ppd_in_tight(start2, finish2, box_offset2, &
               a2_pos, fa_tightbox%start_ppds2, fa_tightbox%finish_ppds2, &
               fa_tightbox%start_pts2, fa_tightbox%finish_pts2, pub_cell%n_pt2)

          call basis_lims_1d_in_ppd_in_tight(start3, finish3, box_offset3, &
               a3_pos, fa_tightbox%start_ppds3, fa_tightbox%finish_ppds3, &
               fa_tightbox%start_pts3, fa_tightbox%finish_pts3, pub_cell%n_pt3)

          ! ndmh: find the first point in the current ppd that needs to be copied
          fa_point = fa_sphere%offset + (fa_ppd-1)*pub_cell%n_pts &
               + (start3-1)*pub_cell%n_pt2*pub_cell%n_pt1 &
               + (start2-1)*pub_cell%n_pt1 + (start1-1)

          ! ndmh: find range of points in '1' direction
          box_start1 = box_offset1 + offset1
          box_end1 = box_start1 + (finish1 - start1)

          do point3=start3,finish3
             ! ndmh: coordinates of points with respect to the box,
             !      i.e. the first point inside the box in a certain
             !      dimension starts from 1, etc.
             box_pt3 = point3 - start3 + box_offset3 + offset3
             do point2=start2,finish2

                box_pt2 = point2 - start2 + box_offset2 + offset2

                ! ndmh: copy '1' direction line of function values to box
                fa_box(box_start1:box_end1,box_pt2,box_pt3) = &
                     fa_on_grid(fa_point:fa_point+(finish1-start1))

                ! ndmh: skip to next line in '2' direction
                fa_point = fa_point + pub_cell%n_pt1

             end do

             ! ndmh: skip to next line in '3' direction
             fa_point = fa_point + (pub_cell%n_pt2-finish2+start2-1)* &
                  pub_cell%n_pt1
          end do

       end do

    end if

    call timer_clock('basis_copy_function_to_box',2)

  end subroutine basis_copy_function_to_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_add_function_to_box(fa_box, &                  ! in-out
       box_n1, box_n2, box_n3, fa_start1, fa_start2, fa_start3, & ! input
       fa_tightbox, fa_on_grid, fa_sphere, factor)                ! input

    !==================================================================!
    ! This subroutine adds an NGWF multiplied by some factor to a      !
    ! specified position with a three-dimensional array with the size  !
    ! of the fftbox, reading the ngwf directly from its ppd            !
    ! representation.                                                  !
    !------------------------------------------------------------------!
    ! Adapted from code originally written by Chris-Kriton Skylaris    !
    ! on 11/6/2001 for the ONES code, for extracting ppds from a box,  !
    ! called "basis_extract_ppds_from_pairbox"                         !
    ! Renamed by Arash A. Mostofi in August 2003 and removed the       !
    ! "make_pair_even" increment from the FFTbox dimensions.           !
    ! Modified for speed by Chris-Kriton Skylaris on 08/09/2005.       !
    ! Modified to work faster when the ppds contain only one point     !
    ! by Chris-Kriton Skylaris on 31/12/2006.                          !
    ! Modified to work in reverse to deposit ppds to the fftbox, and   !
    ! loop-unrolled for speed by Nicholas Hine on 09/05/2008.          !
    ! Modified by Nicholas Hine on 24/06/2009 to combine ppd_loc       !
    ! with ppd_list for MPI efficiency.                                !
    ! Renamed by Nicholas Hine on 28/02/2011 and adapted to work     !
    ! box of any size rather than just fftbox.                       !
    !==================================================================!

    use simulation_cell, only : pub_cell
    use timer, only: timer_clock

    implicit none

    ! Arguments
    integer, intent(in) :: box_n1,box_n2,box_n3
    real(kind=DP), intent(inout) :: fa_box(box_n1,box_n2,box_n3) ! box data
    integer, intent(in) :: fa_start1 ! 1-beginning of NGWF tightbox wrt box
    integer, intent(in) :: fa_start2 ! 2-beginning of NGWF tightbox wrt box
    integer, intent(in) :: fa_start3 ! 3-beginning of NGWF tightbox wrt box
    type(FUNCTION_TIGHT_BOX), intent(in) :: fa_tightbox ! NGWF tightbox
    real(kind=DP), intent(in), dimension(:) :: fa_on_grid ! NGWF to add
    type(SPHERE), intent(in) :: fa_sphere ! sphere of NGWFs PPDs
    real(kind=DP), intent(in) :: factor ! ndmh: factor to multiply NGWF by

    ! Local Variables
    integer :: fa_ppd      ! PPD counter within sphere
    integer :: fa_point    ! point counter within fa_on_grid
    integer :: a1_pos      ! absolute 1-coordinate of current PPD
    integer :: a2_pos      ! absolute 2-coordinate of current PPD
    integer :: a3_pos      ! absolute 3-coordinate of current PPD
    integer :: start1      ! start of points within tightbox along direction 1 of current PPD
    integer :: start2      ! start of points within tightbox along direction 2 of current PPD
    integer :: start3      ! start of points within tightbox along direction 3 of current PPD
    integer :: finish1     ! end of points within tightbox along direction 1 of current PPD
    integer :: finish2     ! end of points within tightbox along direction 2 of current PPD
    integer :: finish3     ! end of points within tightbox along direction 3 of current PPD
    integer :: box_offset1 ! offset to express 1-start of function in box
    integer :: box_offset2 ! offset to express 2-start of function in box
    integer :: box_offset3 ! offset to express 3-start of function in box
    integer :: box_pt1     ! point counter wrt box along direction-1
    integer :: box_pt2     ! point counter wrt box along direction-2
    integer :: box_pt3     ! point counter wrt box along direction-3
    integer :: box_start1  ! start point of box points along direction-1
    integer :: box_end1    ! end point of box points along direction-1
    integer :: point2      ! point counter wrt current PPD along direction-2
    integer :: point3      ! point counter wrt current PPD along direction-3

    call timer_clock('basis_add_function_to_box',1)

    if (pub_cell%n_pts == 1) then
       ! cks: case of only one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               pub_cell%n_ppds_a1, pub_cell%n_ppds_a2, pub_cell%n_ppds_a3)

          box_pt1 = a1_pos - fa_tightbox%start_ppds1 + fa_start1
          box_pt2 = a2_pos - fa_tightbox%start_ppds2 + fa_start2
          box_pt3 = a3_pos - fa_tightbox%start_ppds3 + fa_start3

          fa_point = fa_sphere%offset + fa_ppd - 1

          ! ndmh: add factor*(ngwf value) to box
          fa_box(box_pt1, box_pt2, box_pt3) = &
               fa_box(box_pt1, box_pt2, box_pt3) + &
               fa_on_grid(fa_point)*factor
       enddo

    else
       ! cks: case for more than one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               pub_cell%n_ppds_a1, pub_cell%n_ppds_a2, pub_cell%n_ppds_a3)

          call basis_lims_1d_in_ppd_in_tight(start1, finish1, box_offset1, &
               a1_pos, fa_tightbox%start_ppds1, fa_tightbox%finish_ppds1, &
               fa_tightbox%start_pts1, fa_tightbox%finish_pts1, pub_cell%n_pt1)

          call basis_lims_1d_in_ppd_in_tight(start2, finish2, box_offset2, &
               a2_pos, fa_tightbox%start_ppds2, fa_tightbox%finish_ppds2, &
               fa_tightbox%start_pts2, fa_tightbox%finish_pts2, pub_cell%n_pt2)

          call basis_lims_1d_in_ppd_in_tight(start3, finish3, box_offset3, &
               a3_pos, fa_tightbox%start_ppds3, fa_tightbox%finish_ppds3, &
               fa_tightbox%start_pts3, fa_tightbox%finish_pts3, pub_cell%n_pt3)

          ! cks: loop over the points of the current ppd that also belong to
          !      the box

          fa_point = fa_sphere%offset + (fa_ppd-1)*pub_cell%n_pts - 1 &
               + (start3-1)*pub_cell%n_pt2*pub_cell%n_pt1 &
               + (start2-1)*pub_cell%n_pt1 + start1

          box_start1 = box_offset1 + fa_start1
          box_end1 = box_start1 + finish1 - start1

          do point3=start3,finish3
             ! cks: coordinates of points with respect to the box,
             !      i.e. the first point inside the box in a certain
             !      dimension starts from 1, etc.
             box_pt3 = point3 - start3 + box_offset3 + fa_start3
             do point2=start2,finish2

                box_pt2 = point2 - start2 + box_offset2 + fa_start2

                ! ndmh: removed inner loop for speed
                ! ndmh: add factor*(ngwf value) to box
                fa_box(box_start1:box_end1,box_pt2,box_pt3) = &
                     fa_box(box_start1:box_end1,box_pt2,box_pt3) + &
                     fa_on_grid(fa_point:fa_point+(finish1-start1))*factor

                fa_point = fa_point + pub_cell%n_pt1

             end do

             fa_point = fa_point + &
                  (pub_cell%n_pt2-finish2+start2-1)*pub_cell%n_pt1

          end do

       end do

    endif

    call timer_clock('basis_add_function_to_box',2)

  end subroutine basis_add_function_to_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_extract_function_from_box(fa_on_grid, &          ! in-out
       box_n1, box_n2, box_n3, fa_box, fa_sphere, fa_tightbox, &    ! input
       fa_start1, fa_start2, fa_start3, ngwf_offset)                ! input

    !================================================================!
    ! This subroutines extracts PPDs from an FFTbox.                 !
    ! The PPDs which belong to the NGWF "fa" are extracted           !
    ! and stored in the order dictated by the sphere of "fa".        !
    !----------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris on 11/6/2001 for   !
    ! the ONES code and called "basis_extract_ppds_from_pairbox"     !
    ! Renamed by Arash A. Mostofi in August 2003 and removed the     !
    ! "make_pair_even" increment from the FFTbox dimensions.         !
    ! Modified for speed by Chris-Kriton Skylaris on 08/09/2005.     !
    ! Modified to work faster when the ppds contain only one point   !
    ! by Chris-Kriton Skylaris on 31/12/2006.                        !
    ! Modified by Nicholas Hine on 24/06/2009 to combine ppd_loc     !
    ! with ppd_list for MPI efficiency.                              !
    ! Renamed by Nicholas Hine on 28/02/2011 and adapted to work     !
    ! box of any size rather than just fftbox.                       !
    !================================================================!

    use simulation_cell, only : pub_cell
    use timer, only: timer_clock

    implicit none

    ! Arguments
    integer, intent(in) :: box_n1,box_n2,box_n3
    type(SPHERE), intent(in) :: fa_sphere ! sphere of NGWFs PPDs we want
    type(FUNCTION_TIGHT_BOX), intent(in) :: fa_tightbox ! tightbox "
    integer, intent(in) :: fa_start1 ! 1-beginning of extract region wrt box
    integer, intent(in) :: fa_start2 ! 2-beginning of extract region wrt box
    integer, intent(in) :: fa_start3 ! 3-beginning of extract region wrt box
    integer, intent(in) :: ngwf_offset ! beginning of NGWF in PPD rep
    real(kind=DP), intent(in) :: fa_box(box_n1,box_n2,box_n3) ! extraction box
    real(kind=DP), intent(inout), dimension(:) :: fa_on_grid ! extracted NGWF

    ! Local Variables
    integer :: fa_ppd     ! PPD counter within sphere
    integer :: fa_point   ! point counter within fa_on_grid
    integer :: a1_pos     ! absolute 1-coordinate of current PPD
    integer :: a2_pos     ! absolute 2-coordinate of current PPD
    integer :: a3_pos     ! absolute 3-coordinate of current PPD
    integer :: start1     ! start of points within tightbox along direction 1 of current PPD
    integer :: start2     ! start of points within tightbox along direction 2 of current PPD
    integer :: start3     ! start of points within tightbox along direction 3 of current PPD
    integer :: finish1    ! end of points within tightbox along direction 1 of current PPD
    integer :: finish2    ! end of points within tightbox along direction 2 of current PPD
    integer :: finish3    ! end of points within tightbox along direction 3 of current PPD
    integer :: box_offset1 ! offset to express 1-start of function in box
    integer :: box_offset2 ! offset to express 2-start of function in box
    integer :: box_offset3 ! offset to express 3-start of function in box
    integer :: box_pt1    ! point counter wrt box along direction-2
    integer :: box_pt2    ! point counter wrt box along direction-2
    integer :: box_pt3    ! point counter wrt box along direction-3
    integer :: box_start1 ! start point in box along direction-1
    integer :: box_end1   ! end point in box along direction-1
    integer :: point2     ! point counter wrt current PPD along direction-2
    integer :: point3     ! point counter wrt current PPD along direction-3

    call timer_clock('basis_extract_function_from_box', 1)

    ! cks: initialise
    fa_on_grid(ngwf_offset: &
         ngwf_offset - 1 + fa_sphere%n_ppds_sphere*pub_cell%n_pts) = 0.0_DP

    if (pub_cell%n_pts == 1) then
       ! cks: case of only one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               pub_cell%n_ppds_a1, pub_cell%n_ppds_a2, pub_cell%n_ppds_a3)

          box_pt1 = a1_pos - fa_tightbox%start_ppds1 + fa_start1
          box_pt2 = a2_pos - fa_tightbox%start_ppds2 + fa_start2
          box_pt3 = a3_pos - fa_tightbox%start_ppds3 + fa_start3

          fa_point = ngwf_offset + fa_ppd -1

          fa_on_grid(fa_point) = fa_box(box_pt1, box_pt2, box_pt3)

       enddo

    else
       ! cks: case for more than one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               pub_cell%n_ppds_a1, pub_cell%n_ppds_a2, pub_cell%n_ppds_a3)

          call basis_lims_1d_in_ppd_in_tight(start1, finish1, box_offset1, &
               a1_pos, fa_tightbox%start_ppds1, fa_tightbox%finish_ppds1, &
               fa_tightbox%start_pts1, fa_tightbox%finish_pts1, pub_cell%n_pt1)

          call basis_lims_1d_in_ppd_in_tight(start2, finish2, box_offset2, &
               a2_pos, fa_tightbox%start_ppds2, fa_tightbox%finish_ppds2, &
               fa_tightbox%start_pts2, fa_tightbox%finish_pts2, pub_cell%n_pt2)

          call basis_lims_1d_in_ppd_in_tight(start3, finish3, box_offset3, &
               a3_pos, fa_tightbox%start_ppds3, fa_tightbox%finish_ppds3, &
               fa_tightbox%start_pts3, fa_tightbox%finish_pts3, pub_cell%n_pt3)

          fa_point = ngwf_offset + (fa_ppd-1)*pub_cell%n_pts - 1 &
               + (start3-1)*pub_cell%n_pt2*pub_cell%n_pt1 &
               + (start2-1)*pub_cell%n_pt1 + start1

          box_start1 = box_offset1 + fa_start1
          box_end1 = box_start1 + finish1 - start1

          do point3=start3,finish3
             ! cks: coordinates of points with respect to the FFT box,
             !      i.e. the first point inside the FFT box in a certain
             !      dimension starts from 1, etc.
             box_pt3 = point3 - start3 + box_offset3 + fa_start3
             do point2=start2,finish2

                box_pt2 = point2 - start2 + box_offset2 + fa_start2

                fa_on_grid(fa_point:fa_point+(finish1-start1)) = &
                     fa_box(box_start1:box_end1,box_pt2, box_pt3)

                fa_point = fa_point + pub_cell%n_pt1

             end do

             fa_point = fa_point + &
                  (pub_cell%n_pt2-finish2+start2-1)*pub_cell%n_pt1

          end do

       end do

    endif

    call timer_clock('basis_extract_function_from_box', 2)

  end subroutine basis_extract_function_from_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_multiply_function_by_box(fa_on_grid, &           ! in-out
       box_n1, box_n2, box_n3, fa_box, fa_sphere, fa_tightbox, &    ! input
       fa_start1, fa_start2, fa_start3, ngwf_offset)                ! input

    !================================================================!
    ! This subroutines multiplies the values of a function in PPDs   !
    ! by a function held in a box.                                   !
    ! The PPDs which belong to the NGWF "fa" are multiplied          !
    ! and stored in the order dictated by the sphere of "fa".        !
    !----------------------------------------------------------------!
    ! Adapted by Nicholas Hine from basis_extract_function_from_box  !
    ! on 28/02/2011.                                                 !
    !================================================================!

    use simulation_cell, only : pub_cell
    use timer, only: timer_clock

    implicit none

    ! Arguments
    integer, intent(in) :: box_n1,box_n2,box_n3
    type(SPHERE), intent(in) :: fa_sphere ! sphere of NGWFs PPDs we want
    type(FUNCTION_TIGHT_BOX), intent(in) :: fa_tightbox ! tightbox "
    integer, intent(in) :: fa_start1 ! 1-beginning of region wrt box
    integer, intent(in) :: fa_start2 ! 2-beginning of region wrt box
    integer, intent(in) :: fa_start3 ! 3-beginning of region wrt box
    integer, intent(in) :: ngwf_offset ! beginning of NGWF in PPD rep
    real(kind=DP), intent(in) :: fa_box(box_n1,box_n2,box_n3) ! box
    real(kind=DP), intent(inout), dimension(:) :: fa_on_grid ! multiplied NGWF

    ! Local Variables
    integer :: fa_ppd     ! PPD counter within sphere
    integer :: fa_point   ! point counter within fa_on_grid
    integer :: a1_pos     ! absolute 1-coordinate of current PPD
    integer :: a2_pos     ! absolute 2-coordinate of current PPD
    integer :: a3_pos     ! absolute 3-coordinate of current PPD
    integer :: start1     ! start of points within tightbox along direction 1 of current PPD
    integer :: start2     ! start of points within tightbox along direction 2 of current PPD
    integer :: start3     ! start of points within tightbox along direction 3 of current PPD
    integer :: finish1    ! end of points within tightbox along direction 1 of current PPD
    integer :: finish2    ! end of points within tightbox along direction 2 of current PPD
    integer :: finish3    ! end of points within tightbox along direction 3 of current PPD
    integer :: box_offset1 ! offset to express 1-start of function in box
    integer :: box_offset2 ! offset to express 2-start of function in box
    integer :: box_offset3 ! offset to express 3-start of function in box
    integer :: box_pt1    ! point counter wrt box along direction-2
    integer :: box_pt2    ! point counter wrt box along direction-2
    integer :: box_pt3    ! point counter wrt box along direction-3
    integer :: box_start1 ! start point in box along direction-1
    integer :: box_end1   ! end point in box along direction-1
    integer :: point2     ! point counter wrt current PPD along direction-2
    integer :: point3     ! point counter wrt current PPD along direction-3

    call timer_clock('basis_multiply_function_by_box', 1)

    if (pub_cell%n_pts == 1) then
       ! cks: case of only one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               pub_cell%n_ppds_a1, pub_cell%n_ppds_a2, pub_cell%n_ppds_a3)

          box_pt1 = a1_pos - fa_tightbox%start_ppds1 + fa_start1
          box_pt2 = a2_pos - fa_tightbox%start_ppds2 + fa_start2
          box_pt3 = a3_pos - fa_tightbox%start_ppds3 + fa_start3

          fa_point = ngwf_offset + fa_ppd -1

          fa_on_grid(fa_point) = fa_on_grid(fa_point) * &
               fa_box(box_pt1, box_pt2, box_pt3)
       enddo

    else
       ! cks: case for more than one grid point per ppd

       ! cks: loop over the ppds that belong to the sphere fa
       do fa_ppd=1,fa_sphere%n_ppds_sphere

          call basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, &
               fa_sphere%ppd_list(1,fa_ppd), fa_sphere%ppd_list(2,fa_ppd), &
               pub_cell%n_ppds_a1, pub_cell%n_ppds_a2, pub_cell%n_ppds_a3)

          call basis_lims_1d_in_ppd_in_tight(start1, finish1, box_offset1, &
               a1_pos, fa_tightbox%start_ppds1, fa_tightbox%finish_ppds1, &
               fa_tightbox%start_pts1, fa_tightbox%finish_pts1, pub_cell%n_pt1)

          call basis_lims_1d_in_ppd_in_tight(start2, finish2, box_offset2, &
               a2_pos, fa_tightbox%start_ppds2, fa_tightbox%finish_ppds2, &
               fa_tightbox%start_pts2, fa_tightbox%finish_pts2, pub_cell%n_pt2)

          call basis_lims_1d_in_ppd_in_tight(start3, finish3, box_offset3, &
               a3_pos, fa_tightbox%start_ppds3, fa_tightbox%finish_ppds3, &
               fa_tightbox%start_pts3, fa_tightbox%finish_pts3, pub_cell%n_pt3)

          fa_point = ngwf_offset + (fa_ppd-1)*pub_cell%n_pts - 1 &
               + (start3-1)*pub_cell%n_pt2*pub_cell%n_pt1 &
               + (start2-1)*pub_cell%n_pt1 + start1

          box_start1 = box_offset1 + fa_start1
          box_end1 = box_start1 + finish1 - start1

          do point3=start3,finish3
             ! cks: coordinates of points with respect to the FFT box,
             !      i.e. the first point inside the FFT box in a certain
             !      dimension starts from 1, etc.
             box_pt3 = point3 - start3 + box_offset3 + fa_start3
             do point2=start2,finish2

                box_pt2 = point2 - start2 + box_offset2 + fa_start2

                fa_on_grid(fa_point:fa_point+(finish1-start1)) = &
                     fa_on_grid(fa_point:fa_point+(finish1-start1)) * &
                     fa_box(box_start1:box_end1,box_pt2,box_pt3)

                fa_point = fa_point + pub_cell%n_pt1

             end do

             fa_point = fa_point + &
                  (pub_cell%n_pt2-finish2+start2-1)*pub_cell%n_pt1

          end do

       end do

    endif

    call timer_clock('basis_multiply_function_by_box', 2)

  end subroutine basis_multiply_function_by_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_clean_function(functions_on_grid,  & ! inout
       function_sphere, n_func_ppds)  ! in

    !===============================================================!
    ! This subroutine cleans (shaves) an NGWF in ppd representation !
    ! so that it is zero outside its sphere.                        !
    !---------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 5/7/2001.                 !
    ! Modified to use pub_cell by Quintin Hill on 15/10/2008.       !
    ! Tidied up and moved to basis_mod by Nick Hine on 12/07/2009   !
    !===============================================================!

    use constants, only: DP
    use geometry, only: operator(*), point, local_displacement, operator(+), &
         geometry_distance
    use simulation_cell, only: pub_cell, pub_fftbox

    implicit none

    ! Arguments
    integer, intent(in) :: n_func_ppds
    real(kind=DP), intent(inout) :: functions_on_grid(n_func_ppds*pub_cell%n_pts)
    type(SPHERE), intent(in) :: function_sphere

    ! Local Variables
    integer :: point_counter
    integer :: ppd
    integer :: loc
    integer :: loc_1, loc_2, loc_3
    integer :: a1_neighbour, a2_neighbour, a3_neighbour
    integer :: a1_upper,a2_upper,a3_upper
    integer :: a1_lower,a2_lower,a3_lower
    integer :: a1_cell,a2_cell,a3_cell
    integer :: ppd_count
    logical :: in_region
    real(kind=DP) :: local_1, local_2, local_3
    type(POINT) :: current_point,current_displacement,periodic_centre,ppd_origin

    ! ndmh: to check for ppds included multiple times, when cell=fftbox
    if (pub_fftbox%coin1) then
       a1_lower = -1
       a1_upper = 1
    else
       a1_lower = 0
       a1_upper = 0
    end if
    if (pub_fftbox%coin2) then
       a2_lower = -1
       a2_upper = 1
    else
       a2_lower = 0
       a2_upper = 0
    end if
    if (pub_fftbox%coin3) then
       a3_lower = -1
       a3_upper = 1
    else
       a3_lower = 0
       a3_upper = 0
    end if

    point_counter = function_sphere%offset-1

    ! cks: loop over the ppds that belong to the sphere
    do ppd_count=1,function_sphere%n_ppds_sphere

       ppd = function_sphere%ppd_list(1,ppd_count)
       ppd_origin = basis_ppd_location(ppd)

       loc = function_sphere%ppd_list(2,ppd_count)

       ! cks: loop over points of this ppd
       do loc_3=0,pub_cell%n_pt3 -1
          local_3 = real(loc_3, DP)*pub_cell%d3

          do loc_2=0,pub_cell%n_pt2-1
             local_2 = real(loc_2, DP)*pub_cell%d2

             do loc_1=0,pub_cell%n_pt1-1
                local_1 = real(loc_1, DP)*pub_cell%d1

                current_displacement = local_displacement(pub_cell%a1_unit,&
                     pub_cell%a2_unit,pub_cell%a3_unit,local_1,local_2,local_3)

                current_point = ppd_origin + current_displacement

                point_counter = point_counter + 1

                in_region = .false.

                ! ndmh: loop over possible origins of this NGWF from different cells
                ! ndmh: NB: a1_lower = a1_upper = 0 if pub_fftbox%coin1=.false.
                do a1_cell=a1_lower,a1_upper
                   do a2_cell=a2_lower,a2_upper
                      do a3_cell=a3_lower,a3_upper

                         a1_neighbour=nint(real(loc,kind=DP)/9.0_DP)
                         a2_neighbour=nint(real(loc-9*a1_neighbour,DP)/3.0_DP)
                         a3_neighbour=loc-9*a1_neighbour-3*a2_neighbour

                         ! ndmh: In cases where FFTbox==cell, check over nearest
                         ! ndmh: neighbours regardless of loc value
                         if (pub_fftbox%coin1) a1_neighbour = a1_cell
                         if (pub_fftbox%coin2) a2_neighbour = a2_cell
                         if (pub_fftbox%coin3) a3_neighbour = a3_cell

                         periodic_centre = function_sphere%centre &
                              + real(a1_neighbour,DP)*pub_cell%a1 &
                              + real(a2_neighbour,DP)*pub_cell%a2 &
                              + real(a3_neighbour,DP)*pub_cell%a3

                         if (geometry_distance(periodic_centre,current_point) <= &
                              function_sphere%radius) then
                            in_region = .true.
                         end if

                      enddo
                   enddo
                enddo

                if (.not.in_region) functions_on_grid(point_counter) = 0.0_DP

             enddo
          enddo
       enddo

    enddo

  end subroutine basis_clean_function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  integer function basis_count_psincs(function_sphere)

    !===============================================================!
    ! This subroutine counts the number of unique psinc values in   !
    ! this function sphere.                                         !
    !---------------------------------------------------------------!
    ! Written by Nicholas Hine on 11/11/2010.                       !
    !===============================================================!

    use constants, only: DP
    use geometry, only: operator(*), point, local_displacement, operator(+), &
         geometry_distance
    use simulation_cell, only: pub_cell, pub_fftbox

    implicit none

    ! Arguments
    type(SPHERE), intent(in) :: function_sphere

    ! Local Variables
    integer :: ppd
    integer :: loc
    integer :: loc_1, loc_2, loc_3
    integer :: a1_neighbour, a2_neighbour, a3_neighbour
    integer :: a1_upper,a2_upper,a3_upper
    integer :: a1_lower,a2_lower,a3_lower
    integer :: a1_cell,a2_cell,a3_cell
    integer :: ppd_count
    logical :: in_region
    real(kind=DP) :: local_1, local_2, local_3
    type(POINT) :: current_point,current_displacement,periodic_centre,ppd_origin

    ! ndmh: to check for ppds included multiple times, when cell=fftbox
    if (pub_fftbox%coin1) then
       a1_lower = -1
       a1_upper = 1
    else
       a1_lower = 0
       a1_upper = 0
    end if
    if (pub_fftbox%coin2) then
       a2_lower = -1
       a2_upper = 1
    else
       a2_lower = 0
       a2_upper = 0
    end if
    if (pub_fftbox%coin3) then
       a3_lower = -1
       a3_upper = 1
    else
       a3_lower = 0
       a3_upper = 0
    end if

    basis_count_psincs = 0

    ! cks: loop over the ppds that belong to the sphere
    do ppd_count=1,function_sphere%n_ppds_sphere

       ppd = function_sphere%ppd_list(1,ppd_count)
       ppd_origin = basis_ppd_location(ppd)

       loc = function_sphere%ppd_list(2,ppd_count)

       ! cks: loop over points of this ppd
       do loc_3=0,pub_cell%n_pt3 -1
          local_3 = real(loc_3, DP)*pub_cell%d3

          do loc_2=0,pub_cell%n_pt2-1
             local_2 = real(loc_2, DP)*pub_cell%d2

             do loc_1=0,pub_cell%n_pt1-1
                local_1 = real(loc_1, DP)*pub_cell%d1

                current_displacement = local_displacement(pub_cell%a1_unit,&
                     pub_cell%a2_unit,pub_cell%a3_unit,local_1,local_2,local_3)

                current_point = ppd_origin + current_displacement

                in_region = .false.

                ! ndmh: loop over possible origins of this NGWF from different cells
                ! ndmh: NB: a1_lower = a1_upper = 0 if pub_fftbox%coin1=.false.
                do a1_cell=a1_lower,a1_upper
                   do a2_cell=a2_lower,a2_upper
                      do a3_cell=a3_lower,a3_upper

                         a1_neighbour=nint(real(loc,kind=DP)/9.0_DP)
                         a2_neighbour=nint(real(loc-9*a1_neighbour,DP)/3.0_DP)
                         a3_neighbour=loc-9*a1_neighbour-3*a2_neighbour
                         if (pub_fftbox%coin1) a1_neighbour = a1_cell
                         if (pub_fftbox%coin2) a2_neighbour = a2_cell
                         if (pub_fftbox%coin3) a3_neighbour = a3_cell

                         periodic_centre = function_sphere%centre &
                              + real(a1_neighbour,DP)*pub_cell%a1 &
                              + real(a2_neighbour,DP)*pub_cell%a2 &
                              + real(a3_neighbour,DP)*pub_cell%a3

                         if (geometry_distance(periodic_centre,current_point) <= &
                              function_sphere%radius) then
                            in_region = .true.
                         end if

                      enddo
                   enddo
                enddo

                if (in_region) basis_count_psincs = basis_count_psincs + 1

             enddo
          enddo
       enddo

    enddo

  end function basis_count_psincs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_put_tightbox_in_fftbox(fftbox_inout, &   ! input/output
       start1, start2, start3, &                            ! input
       tightbox_in, npts1, npts2, npts3, factor)            ! input

    !================================================!
    ! Deposits tightbox data into an FFTbox.         !
    !------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/2/2004. !
    !================================================!

    use simulation_cell, only: pub_fftbox
    use timer, only: timer_clock

    implicit none

    ! Arguments
    integer, intent(in)        :: start1, start2, start3
    integer, intent(in)        :: npts1, npts2, npts3
    real(kind =DP), intent(in)    :: factor
    ! pdh: changed to avoid unnecessary copy onto stack
    real(kind =DP), intent(in)    :: tightbox_in(:,:,:)
    real(kind =DP), intent(inout) :: fftbox_inout(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2, pub_fftbox%total_pt3)

    ! cks: << local variables>>
    integer :: end1
    integer :: end2
    integer :: end3

    call timer_clock('basis_put_tightbox_in_fftbox',1)

    end1 = start1 + npts1 - 1
    end2 = start2 + npts2 - 1
    end3 = start3 + npts3 - 1

    ! cks: deposit tightbox (multiplied by factor) to the appropriate
    ! cks: position in fftbox_inout
    fftbox_inout(start1: end1, start2: end2, start3: end3) = &
         fftbox_inout(start1: end1, start2: end2, start3: end3) &
         +factor *tightbox_in(1: npts1, 1: npts2, 1: npts3)

    call timer_clock('basis_put_tightbox_in_fftbox',2)

  end subroutine basis_put_tightbox_in_fftbox


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_copy_tightbox_to_fftbox(box_out, &
       start1, start2, start3, &
       data_in, npts1, npts2, npts3)

    !===========================================!
    ! Copies tightbox into FFTbox.              !
    !-------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2001. !
    !===========================================!

    use simulation_cell, only: pub_fftbox
    use timer, only: timer_clock

    implicit none

    ! Arguments
    integer, intent(in)        :: start1,start2,start3
    integer, intent(in)        :: npts1,npts2,npts3
    real(kind=DP), intent(in)  :: data_in(:,:,:)
    real(kind=DP), intent(out) :: box_out(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2, pub_fftbox%total_pt3)

    call timer_clock('basis_copy_tightbox_to_fftbox',1)

    ! cks: initialise pair_data_row
    box_out = 0.0_DP

    ! cks: copy all the elements of data_in to the appropriate
    !      position in box_out
    box_out(start1:start1+npts1-1,start2:start2+npts2-1, &
         start3:start3+npts3-1) = data_in(1:npts1, 1:npts2, 1:npts3)

    call timer_clock('basis_copy_tightbox_to_fftbox',2)

  end subroutine basis_copy_tightbox_to_fftbox


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_lims_1d_in_ppd_in_tight(start, finish, box_offset, &
       pos, start_ppds, finish_ppds, start_pts, finish_pts, n_pt)

    !==================================================================!
    ! This subroutine works in 1d and given a tightbox and a ppd       !
    ! inside it returns the starting and ending grid points of the ppd !
    ! within the tightbox and the number of points of the tightbox     !
    ! before the starting grid point of the current ppd.               !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2001, improved in 2007.      !
    !==================================================================!

    ! Arguments
    implicit none
    integer, intent(out) :: start      ! first point (of ppd) in tightbox wrt start of ppd
    integer, intent(out) :: finish     ! last point (of ppd) in tightbox wrt start of ppd
    integer, intent(out) :: box_offset ! number of points in tightbox belonging to previous ppds

    integer, intent(in) :: pos         ! position of current ppd
    integer, intent(in) :: start_ppds  ! first ppd of tightbox
    integer, intent(in) :: finish_ppds ! last ppd of tightbox
    integer, intent(in) :: start_pts   ! first point in tightbox of current ppd
    integer, intent(in) :: finish_pts  ! last point in tightbox of current ppd
    integer, intent(in) :: n_pt        ! number of points in current ppd

    if (pos == start_ppds .and. pos == finish_ppds) then
       ! cks: the box contains only 1 ppd in this dimension
       start = start_pts
       finish = finish_pts
       box_offset = 0
    else if (pos == start_ppds) then
       ! cks: this is the first ppd in the box in this dimension
       start = start_pts
       finish = n_pt
       box_offset = 0
    else if (pos == finish_ppds) then
       ! cks: this is the last ppd in the box in this dimension
       start = 1
       finish = finish_pts
       box_offset = -start_pts + 1 + (finish_ppds-start_ppds) * n_pt
    else
       ! cks: this is neither the first, nor the last ppd in the fftbox in
       !      this dimension
       start = 1
       finish = n_pt
       box_offset = -start_pts + 1 + (pos-start_ppds) * n_pt
    end if

  end subroutine basis_lims_1d_in_ppd_in_tight


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine basis_find_ppd_in_neighbour(a1_pos, a2_pos, a3_pos, & ! output
       global_count, global_loc, n_ppds_a1, n_ppds_a2, n_ppds_a3)  ! input

    !==========================================================!
    ! This subroutine returns the integer PPD-grid coordinates !
    ! of the current PPD along each lattice vector direction.  !
    ! These coordinates are "absolute" in the sense that they  !
    ! belong to PPDs which may not even be inside the          !
    ! simulation cell so that they always belong to an NGWF    !
    ! sphere whose centre is inside the simulation cell.       !
    !----------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2001.                !
    ! Modified by Chris-Kriton Skylaris on 21/06/2007 to work  !
    ! also when simulation cell and FFT box coincide in one or !
    ! more dimensions.                                         !
    ! Modified by Nicholas Hine in June 2008 for optimal speed !
    !==========================================================!

    use simulation_cell, only: pub_fftbox
    implicit none

    ! Arguments
    integer, intent(out) :: a1_pos, a2_pos, a3_pos
    integer, intent(in) :: global_count, global_loc, &
         n_ppds_a1, n_ppds_a2, n_ppds_a3

    ! cks: internal variable declarations
    integer ::  map_a1a2, a1_neighbour, a2_neighbour, a3_neighbour
    integer, parameter :: neighbours(-13:13,1:3) = &
         reshape((/-1, -1, -1, -1, -1, -1, -1, -1, -1, &
         0,  0,  0,  0,  0,  0,  0,  0,  0, &
         +1, +1, +1, +1, +1, +1, +1, +1, +1, &
         -1, -1, -1,  0,  0,  0, +1, +1, +1, &
         -1, -1, -1,  0,  0,  0, +1, +1, +1, &
         -1, -1, -1,  0,  0,  0, +1, +1, +1, &
         -1,  0, +1, -1,  0, +1, -1,  0, +1, &
         -1,  0, +1, -1,  0, +1, -1,  0, +1, &
         -1,  0, +1, -1,  0, +1, -1,  0, +1 /),(/27,3/))

    ! These are the coordinates of the current ppd in the simulation cell
    a3_pos = (global_count-1) / (n_ppds_a1*n_ppds_a2) + 1
    map_a1a2 = global_count - (n_ppds_a1*n_ppds_a2) * (a3_pos - 1)

    a2_pos = (map_a1a2 - 1) / n_ppds_a1 + 1

    a1_pos = map_a1a2 - n_ppds_a1 * (a2_pos-1)

    ! If the current ppd is not due to a function in the simulation cell
    ! but due to a function in a neighbouring simulation cell
    ! (periodic image), find the coordinates of the ppd that corresponds
    ! to the actual function inside the simulation cell and therefore
    ! is a ppd outside the simulation cell.

    ! ndmh: use lookup table to avoid floating-point arithmetic since this
    ! ndmh: function is called millions of times per iteration
    a1_neighbour = neighbours(global_loc, 1)
    a2_neighbour = neighbours(global_loc, 2)
    a3_neighbour = neighbours(global_loc, 3)

    ! cks: these are now the positions of ppds that can be even outside the
    ! cks: simulation cell, but they all belong to a single function whose
    ! cks: centre is INSIDE the simulation cell.
    if (.not.pub_fftbox%coin3) a3_pos = a3_pos - a3_neighbour * n_ppds_a3
    if (.not.pub_fftbox%coin2) a2_pos = a2_pos - a2_neighbour * n_ppds_a2
    if (.not.pub_fftbox%coin1) a1_pos = a1_pos - a1_neighbour * n_ppds_a1

  end subroutine basis_find_ppd_in_neighbour


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module basis
