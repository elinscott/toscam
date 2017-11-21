! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!========================================================================!
!                          PBC corrections module                        !
!------------------------------------------------------------------------!
! Implements the Martyna-Tuckerman correction to Hartree and Vloc terms, !
! according to [1] Martyna, Tuckerman, 1999 J. Chem. Phys. 110, 2810).   !
! Written by Jacek Dziedzic under the supervision of Chris Skylaris      !
! in January 2010.                                                       !
!========================================================================!
! The Martyna-Tuckerman screening function only needs to be computed     !
! once. This is done in pbc_corr_initialise and the result is stored     !
! in pbc_corr_screening_term. Initialization is performed automatically  !
! on first use of pbc_corr_hartree or pbc_corr_vloc. Clean-up is         !
! realized by pbc_corr_exit.                                             !
!========================================================================! 

module pbc_corrections

  use constants, only: DP
  use utils, only: utils_trace_in, utils_trace_out

  implicit none

  private

  ! Public subroutines
  public :: pbc_corr_check_cell_size
  public :: pbc_corr_hartree
  public :: pbc_corr_initialise
  public :: pbc_corr_exit
  public :: pbc_corr_vloc

  ! Private variables
  logical :: pbc_corr_is_initialised = .false.
  real(kind=DP), allocatable :: pbc_corr_screening_term(:,:,:)


contains

  subroutine pbc_corr_initialise(grid)

    !==========================================================================!
    ! Allocates memory and precalculates the Martyna-Tuckerman screening       !
    ! function. There is no need to call this routine, because it is invoked   !
    ! automatically upon first use of pbc_corr_hartree of pbc_corr_vloc,       !
    ! yet if need arises to do the initialization earlier, it may be called.   !
    !--------------------------------------------------------------------------!
    ! Arguments: None                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, January 2010                                  !
    !==========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_node_id, comms_abort, comms_reduce
    use constants, only: DP, PI, stdout
    use fourier, only: fourier_apply_cell_forward
    use geometry, only: magnitude, POINT
    use rundat, only: pub_mt_cutoff
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_erf

    implicit none

    ! Arguments
    type(GRID_INFO), intent(inout) :: grid

    ! Local variables
    integer :: ierr                                ! Error flag
    real(kind=DP), parameter :: fourpi = 4.0_DP*PI ! 4 pi

    real(kind=DP) :: alpha               ! Martyna-Tuckermann cutoff
    real(kind=DP) :: alphasq             ! alpha**2
    real(kind=DP) :: a1,a2,a3            ! cell vectors
    integer :: islab12, i1,i2,i3,islab23 ! indices on the fine grid
    real(kind=DP) :: r                   ! dist. from (0,0,0) to fine grid point
    real(kind=DP) :: gvec(3)             ! G vector
    real(kind=DP) :: int_of_erf_alpha_r_over_r  ! integral of erf(alpha*r)/r over cell 
    complex(kind=DP), allocatable :: screening_term_cplx(:,:,:)  ! workspace
    real(kind=DP), allocatable :: erf_alpha_r_over_r(:,:,:)      ! workspace
    real(kind=DP), allocatable :: one_over_r(:,:,:)              ! workspace

    !------------------------------------------------------------------------

    if(pbc_corr_is_initialised) return

    ! jd: Start timer
    call timer_clock('pbc_corr_initialise',1)
    call utils_trace_in('pbc_corr_initialise')

    ! jd: Compute alpha from alpha*L given in the input file
    a1 = magnitude(pub_cell%a1)
    a2 = magnitude(pub_cell%a2)
    a3 = magnitude(pub_cell%a3)
    alpha = pub_mt_cutoff / min(a1,a2,a3)
    alphasq = alpha * alpha

    ! jd: Allocate persistent array
    allocate(pbc_corr_screening_term(grid%ld3,&
         grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('pbc_corr_initialise','pbc_corr_screening_term',ierr)

    ! jd: Allocate work arrays
    allocate(screening_term_cplx(grid%ld3,&
         grid%ld2,grid%max_slabs23),stat=ierr)
    call utils_alloc_check('pbc_corr_initialise','screening_term_cplx',ierr)
    allocate(erf_alpha_r_over_r(grid%ld1, &
         grid%ld2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('pbc_corr_initialise','erf_alpha_r_over_r',ierr)
    allocate(one_over_r(grid%ld1, &
         grid%ld2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('pbc_corr_initialise','one_over_r',ierr)

    ! jd: Will fill these arrays in a moment, now zero what's between pt and ld
    erf_alpha_r_over_r = 0.0_DP
    one_over_r = 0.0_DP

    ! jd: Compute erf(alpha*r)/r and 1/r on the fine grid
    do islab12 = 1, grid%num_my_slabs12
       i3 = grid%first_slab12(pub_my_node_id) + islab12 - 1

       do i2 = 1,grid%n2
          do i1 = 1,grid%n1

             r = fine_grid_point_to_min_img_loc(i1,i2,i3)

             ! jd: Calculate erf(alpha*r)/r, considering r==0 separately
             if (r < 1E-10_DP) then
                ! jd: lim_{r->0} erf(alpha*r)/r = 2*alpha/sqrt(pi)
                erf_alpha_r_over_r(i1,i2,islab12) = 2.0_DP*alpha/sqrt(PI) * &
                     grid%weight
             else
                erf_alpha_r_over_r(i1,i2,islab12) = utils_erf(alpha*r)/r * &
                     grid%weight
             end if

          end do ! over i1
       end do ! over i2
    end do ! over 12-slabs

    ! jd: Calculate integral of erf(alpha/r)/r over the cell
    int_of_erf_alpha_r_over_r = sum(erf_alpha_r_over_r)
    call comms_reduce('SUM',int_of_erf_alpha_r_over_r)

    ! jd: Fourier transform erf(alpha*r)/r to reciprocal space
    call fourier_apply_cell_forward(erf_alpha_r_over_r,screening_term_cplx, &
         grid)

    ! jd: Strip the zero imaginary part to save memory
    pbc_corr_screening_term = real(screening_term_cplx,kind=DP)

    ! jd: Take care of G=0 term
    if (grid%first_slab23(pub_my_node_id)==1) then

       ! jd: Eq. B3 in [1]. This is *much* more accurate than trying to
       !     integrate 1/r over all cell because it avoids the singularity.
       !     Although an analytical expression can be obtained for the
       !     integral of 1/r over all cell, the approach below is still prefera-
       !     ble, because there are no complications due to minimum images
       !     inherent in the analytical approach (and the integral over all
       !     cell is not strictly equal to 8 * integral over lower corner either
       !     as lines of x=0, y=0 and z=0 are 8-ple counted).
       pbc_corr_screening_term(1,1,1) = int_of_erf_alpha_r_over_r + &
            PI/(alpha*alpha)

    end if

    ! Loop over reciprocal space grid points
    do islab23=1,grid%num_slabs23                ! along b1
       do i2=1,grid%n2                           ! along b2
          do i3=1,grid%n3                        ! along b3

             call cell_grid_recip_pt(gvec,islab23 + &
                  grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)

             ! jd: Subtract \tilde{phi} (Coulomb,long) (eqs. 2.9 and B2 in [1])
             pbc_corr_screening_term(i3,i2,islab23) = &
                  pbc_corr_screening_term(i3,i2,islab23) - &
                  fourpi * grid%coulomb_recip(i3,i2,islab23) * &
                  exp(-0.25_DP*sum(gvec(1:3)**2)/alphasq)
                  
          end do
       end do
    end do

    ! jd: Clean up everything but the screening term array, which must persist
    deallocate(one_over_r,stat=ierr)
    call utils_dealloc_check('pbc_corr_initialise','one_over_r',ierr)
    deallocate(erf_alpha_r_over_r,stat=ierr)
    call utils_dealloc_check('pbc_corr_initialise','erf_alpha_r_over_r',ierr)
    deallocate(screening_term_cplx,stat=ierr)
    call utils_dealloc_check('pbc_corr_initialise','screening_term_cplx',ierr)

    pbc_corr_is_initialised = .true.

    ! jd: Stop timer
    call utils_trace_out('pbc_corr_initialise')
    call timer_clock('pbc_corr_initialise',2)

  contains

    real(kind=DP) function fine_grid_point_to_min_img_loc(i1,i2,i3) 
      !========================================================================!
      ! Returns the distance from (0,0,0) to a point on the fine grid point,   !
      ! given the indices of the point. Minimum image convention is applied.   !
      !------------------------------------------------------------------------!
      ! Arguments:                                                             !
      ! i1, i2, i3 (input): indices of the grid point.                         !
      ! Returns:                                                               !
      ! Distance from (0,0,0) to minimum image of specified fine grid point.   !
      ! Caveats:                                                               !
      ! Last argument is i3, _not_ islab12                                     !
      !------------------------------------------------------------------------!
      ! Written by Jacek Dziedzic, January 2010                                !
      !========================================================================!

      use geometry, only: POINT, operator(+), operator(*), magnitude
      use simulation_cell, only: pub_cell

      implicit none

      ! Arguments: 
      integer, intent(in) :: i1,i2,i3 ! Indices of the grid point

      ! Local variables
      type(POINT) :: rvec0, rvec
      real(kind=DP) :: r, rmin
      integer :: ii1, ii2, ii3
      !------------------------------------------------------------------------

      ! jd: Calculate physical location of grid point in simulation cell
      rvec0 = &
           real((i1-1),kind=DP) * grid%da1 + &
           real((i2-1),kind=DP) * grid%da2 + &
           real((i3-1),kind=DP) * grid%da3

      ! jd: Check eight candidate images to see which one is the minimum one,
      !     including this one. Since we're measuring distance from (0,0,0)
      !     there is no need to examine 27 cells, just 8.
      rmin = huge(1.0_DP)
      do ii1 = -1, 0
         do ii2 = -1, 0
           do ii3 = -1, 0
              rvec = rvec0 + &
                   real(ii1,kind=DP) * pub_cell%a1 + &
                   real(ii2,kind=DP) * pub_cell%a2 + &
                   real(ii3,kind=DP) * pub_cell%a3
              r = magnitude(rvec)
              if (r < rmin) rmin=r
           end do
         end do
      end do

      fine_grid_point_to_min_img_loc = rmin

    end function fine_grid_point_to_min_img_loc
    ! --------------------------------------------------

  end subroutine pbc_corr_initialise
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine pbc_corr_exit
    !==========================================================================!
    ! Cleans up the persistent array allocated in pbc_corr_initialise.         !
    ! If pbc_corr_initialise was not called, does nothing.                     !
    !--------------------------------------------------------------------------!
    ! Arguments: None                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, January 2010                                  !
    !==========================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Local variables
    integer :: ierr                                ! Error flag
    !------------------------------------------------------------------------

    if(.not. pbc_corr_is_initialised) return
    call utils_trace_in('pbc_corr_exit')

    deallocate(pbc_corr_screening_term,stat=ierr)
    call utils_dealloc_check('pbc_corr_exit','pbc_corr_screening_term',ierr)

    ! jd: Needed for TS search, where we'll have to re-initialise in a moment
    pbc_corr_is_initialised = .false.

    call utils_trace_out('pbc_corr_exit')

  end subroutine pbc_corr_exit
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine pbc_corr_hartree(zwork,grid)
    !==========================================================================!
    ! Calculates the Hartree energy, applying the Martyna-Tuckerman technique. !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! zwork (input):  Electronic density on reciprocal grid.                   !
    ! zwork (output): Corrected Hartree potential on reciprocal grid.          !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, January 2010                                  !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use constants, only: DP, PI
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(GRID_INFO), intent(inout) :: grid
    complex(kind=DP), intent(inout) :: zwork( &
         grid%ld3, grid%ld2, grid%max_slabs23) 

    ! Local variables
    real(kind=DP), parameter :: fourpi = 4.0_DP*PI ! 4 pi
    !------------------------------------------------------------------------

    ! jd: Start timer
    call timer_clock('pbc_corr_hartree',1)
    call utils_trace_in('pbc_corr_hartree')

    ! jd: Initialise
    if(.not. pbc_corr_is_initialised) call pbc_corr_initialise(grid)

    ! jd: Apply the correction: V_H = n * ( 4pi/g^2 + screening term(g) )
    zwork = zwork * (grid%coulomb_recip(:,:,:) * fourpi + &
         pbc_corr_screening_term)

    ! jd: Stop timer
    call utils_trace_out('pbc_corr_hartree')
    call timer_clock('pbc_corr_hartree',2)

  end subroutine pbc_corr_hartree
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine pbc_corr_vloc(i3,i2,islab23,ion_charge,v_loc_value,grid)
    !==========================================================================!
    ! Applies the Martyna-Tuckerman correction to the local pseudopotential    !
    ! on a reciprocal grid point.                                              !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! i3, i2, islab23 (input): Indices of the reciprocal grid point.           !
    ! ion_charge (input): Charge of the species due to which the v_loc is.     !
    ! v_loc_value (input): Uncorrected value of the v_loc.                     !
    ! v_lov_value (output): Corrected value of the v_loc.                      !
    ! Caveats:                                                                 !
    ! The v_loc_value that is passed in is expected to be the one _before_     !
    ! scaling by 1/weight_fine.                                                !
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, January 2010                                  !
    !==========================================================================!

    use cell_grid, only: GRID_INFO

    implicit none

    ! Arguments
    type(GRID_INFO), intent(inout) :: grid
    integer, intent(in) :: i3,i2,islab23
    real(kind=DP), intent(in) :: ion_charge
    real(kind=DP), intent(inout) :: v_loc_value
    !------------------------------------------------------------------------

    call utils_trace_in('pbc_corr_vloc')

    ! jd: Initialise
    if(.not. pbc_corr_is_initialised) call pbc_corr_initialise(grid)

    ! jd: Apply the correction
    v_loc_value = v_loc_value - ion_charge * &
         pbc_corr_screening_term(i3,i2,islab23)

    call utils_trace_out('pbc_corr_vloc')

  end subroutine pbc_corr_vloc
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine pbc_corr_check_cell_size(elements)
    !==========================================================================!
    ! Checks that the simulation cell is at least twice as large, along every  !
    ! axis, as the localized charge. If not, a warning is produced.            !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! elements(input): The array describing the ions.                          ! 
    !--------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, February 2010                                 !
    !==========================================================================!

    use constants, only: DP, stdout
    use geometry, only: POINT, operator(.DOT.), magnitude
    use ion, only: ELEMENT
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! Local variables
    type(POINT) :: unit_vectors(3), cell_vectors(3)
    real(kind=DP) :: this_position_plus_ngwf(3)
    real(kind=DP) :: this_position_minus_ngwf(3)
    real(kind=DP) :: min_pos(3), max_pos(3)
    integer :: i, axis

    ! jd: Useful constants
    unit_vectors(1) =  pub_cell%a1_unit
    unit_vectors(2) =  pub_cell%a2_unit
    unit_vectors(3) =  pub_cell%a3_unit
    cell_vectors(1) =  pub_cell%a1
    cell_vectors(2) =  pub_cell%a2
    cell_vectors(3) =  pub_cell%a3

    ! jd: Project the positions of all ions onto the cell vectors, taking into
    !     account the NGWF around every ion. The NGWF localization sphere is
    !     unaffected by the projection. Then find the maximum and minimum.
    do axis=1, 3
       min_pos(axis)=huge(1.0_DP)
       max_pos(axis)=0.0_DP
    end do

    do i=1, pub_cell%nat
       do axis=1, 3
          this_position_plus_ngwf(axis) = &
               (elements(i)%centre .DOT. unit_vectors(axis))+elements(i)%radius
          this_position_minus_ngwf(axis) = &
               (elements(i)%centre .DOT. unit_vectors(axis))-elements(i)%radius
          if(this_position_plus_ngwf(axis) > max_pos(axis)) then 
             max_pos(axis) = this_position_plus_ngwf(axis)
          end if
          if(this_position_minus_ngwf(axis) < min_pos(axis)) then 
             min_pos(axis) = this_position_minus_ngwf(axis)
          end if
       end do
    end do

    ! jd: Check if the cell is at least twice as large along each axis
    do axis=1, 3
       if(max_pos(axis)-min_pos(axis) > &
            magnitude(cell_vectors(axis)) * 0.5_DP) then
          write(stdout,'(a)') 'WARNING: When using pbc_correction&
               &_cutoff, the simulation cell should be at least TWICE as large &
               &as the extent of the localized charge, along every axis.'
          write(stdout,'(a,i1,a,f10.5,a,f10.5,a)') &
               '         This is not the case for axis ', axis, &
               ' where the charge density extends across ', &
               max_pos(axis)-min_pos(axis),' Bohr and the cell size is ', &
               magnitude(cell_vectors(axis)),' Bohr.'
          write(stdout,'(a,f10.5,a)') '         Increase the length of the &
               &cell to at least ', 2.0_DP * (max_pos(axis)-min_pos(axis)), &
               ' or expect the results to be unreliable.'
       end if
    end do

  end subroutine pbc_corr_check_cell_size
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module pbc_corrections

