! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Nicholas D.M. Hine
!
!   based in many cases on existing versions by Chris-Kriton Skylaris,
!   Arash A. Mostofi, Nicholas D.M. Hine and Peter D. Haynes
!
!   Several extra routines have been added by Laura E. Ratcliff
!
!   Thomas Young Centre
!   Imperial College London
!   Exhibition Road
!
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module projectors

  use constants, only: DP, PI, stdout
  use geometry, only: POINT

  implicit none

  private

  ! Derived type to store details of a set of projectors
  type PROJECTOR_SET

     ! ndmh: number of projector species
     integer :: n_proj_species

     ! ndmh: definition of radial reciprocal space grid for each species
     real(kind=DP), allocatable :: gmax(:)
     integer, allocatable :: n_rad_pts(:)

     ! ndmh: projectors on radial reciprocal space grid
     integer, allocatable :: num_shells(:)
     integer, allocatable :: ang_mom(:,:)
     real(kind=DP), allocatable :: rad_proj_recip(:,:,:)

     ! ndmh: centres of the projectors of each atom
     type(POINT), allocatable, dimension(:) :: proj_centre

     ! ndmh: max radius of all the projectors on each atom
     real(kind=DP), allocatable, dimension(:) :: proj_max_radius

     ! ndmh: species of each atom
     integer, allocatable, dimension(:) :: proj_species

     ! ndmh: number of projectors on atoms of each species
     integer, allocatable, dimension(:) :: species_num_proj

     ! ndmh: first entry for this species in fftbox_proj_recip
     integer, allocatable, dimension(:) :: species_first_proj

     ! ndmh: reciprocal space FFTbox containing each type of projector
     complex(kind=DP), allocatable, dimension(:,:,:,:) :: fftbox_proj_recip

     ! ndmh: real-space projectors in PPDs (optional)
     real(kind=DP), allocatable :: projs_on_grid(:)

     ! ndmh: whether to normalise the projectors of this set
     logical :: normalise

  end type PROJECTOR_SET

  public :: PROJECTOR_SET

  public :: projectors_allocate_set
  public :: projectors_deallocate_set
  public :: projectors_init_fftbox_recip
  public :: projectors_exit_fftbox_recip
  public :: projector_in_box_recip
  public :: projector_in_box_recip_gradg
  public :: projectors_create_real
  public :: projectors_func_ovlp_box
  public :: projectors_func_grad_ovlp_box
  public :: projectors_func_pos_ovlp_box
  public :: projectors_gradient_batch
  public :: projectors_commutator

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_allocate_set(proj_set,max_shells,max_npts)

    !==============================================================!
    ! This subroutine allocates storage for a set of projectors    !
    !--------------------------------------------------------------!
    ! Written by Nicholas Hine on 19/05/2010.                      !
    !==============================================================!

    use simulation_cell, only : pub_cell, pub_fftbox
    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(PROJECTOR_SET), intent(inout) :: proj_set
    integer, intent(in) :: max_shells, max_npts

    ! Local Variables
    integer :: ierr      ! error flag

    ! Allocate internal arrays in proj_set relating to projector species
    allocate(proj_set%species_num_proj(proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%species_num_proj',ierr)
    allocate(proj_set%species_first_proj(proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%species_first_proj',ierr)
    allocate(proj_set%gmax(proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%gmax',ierr)
    allocate(proj_set%n_rad_pts(proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%n_rad_pts',ierr)
    allocate(proj_set%num_shells(proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%num_shells',ierr)
    allocate(proj_set%ang_mom(max_shells,proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%ang_mom',ierr)
    allocate(proj_set%rad_proj_recip(max_npts,max_shells, &
         proj_set%n_proj_species),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%rad_proj_recip',ierr)

    ! Allocate arrays in proj_set relating to actual projectors
    allocate(proj_set%proj_centre(pub_cell%nat),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%proj_centre',ierr)
    allocate(proj_set%proj_max_radius(pub_cell%nat),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%proj_radius',ierr)
    allocate(proj_set%proj_species(pub_cell%nat),stat=ierr)
    call utils_alloc_check('projectors_allocate_set', &
         'proj_set%proj_species',ierr)

    ! Default is no normalisation
    proj_set%normalise = .false.

  end subroutine projectors_allocate_set


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_deallocate_set(proj_set)

    !==============================================================!
    ! This subroutine deallocates storage for a set of projectors  !
    !--------------------------------------------------------------!
    ! Written by Nicholas Hine on 19/05/2010.                      !
    !==============================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(PROJECTOR_SET), intent(inout) :: proj_set

    ! Local Variables
    integer :: ierr      ! error flag

    ! Deallocate storage of evaluated projectors in real/recip space
    if (allocated(proj_set%projs_on_grid)) then
       deallocate(proj_set%projs_on_grid,stat=ierr)
       call utils_dealloc_check('projectors_deallocate_set', &
            'proj_set%projs_on_grid',ierr)
    end if
    if (allocated(proj_set%fftbox_proj_recip)) then
       deallocate(proj_set%fftbox_proj_recip,stat=ierr)
       call utils_dealloc_check('projectors_deallocate_set', &
            'proj_set%fftbox_proj_recip',ierr)
    end if

    ! Deallocate arrays in proj_set relating to actual projectors
    deallocate(proj_set%proj_species,stat=ierr)
    call utils_dealloc_check('projectors_deallocate_set', &
         'proj_set%proj_species',ierr)
    deallocate(proj_set%proj_max_radius,stat=ierr)
    call utils_dealloc_check('projectors_deallocate_set', &
         'proj_set%proj_radius',ierr)
    deallocate(proj_set%proj_centre,stat=ierr)
    call utils_dealloc_check('projectors_deallocate_set', &
         'proj_set%proj_centre',ierr)

    ! Deallocate internal arrays in proj_set relating to projector species
    deallocate(proj_set%rad_proj_recip,stat=ierr)
    call utils_dealloc_check('projectors_allocate_set', &
         'proj_set%rad_proj_recip',ierr)
    deallocate(proj_set%ang_mom,stat=ierr)
    call utils_dealloc_check('projectors_allocate_set', &
         'proj_set%ang_mom',ierr)
    deallocate(proj_set%num_shells,stat=ierr)
    call utils_dealloc_check('projectors_allocate_set', &
         'proj_set%num_shells',ierr)
    deallocate(proj_set%n_rad_pts,stat=ierr)
    call utils_dealloc_check('projectors_allocate_set', &
         'proj_set%n_rad_pts',ierr)
    deallocate(proj_set%gmax,stat=ierr)
    call utils_dealloc_check('projectors_allocate_set', &
         'proj_set%gmax',ierr)
    deallocate(proj_set%species_first_proj,stat=ierr)
    call utils_dealloc_check('projectors_deallocate_set', &
         'proj_set%species_first_proj',ierr)
    deallocate(proj_set%species_num_proj,stat=ierr)
    call utils_dealloc_check('projectors_deallocate_set', &
         'proj_set%species_num_proj',ierr)

  end subroutine projectors_deallocate_set


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_init_fftbox_recip(proj_set,kpt,delta,cart,swap_rc)

    !==============================================================!
    ! This subroutine initialises a set of projectors in FFTboxes  !
    !--------------------------------------------------------------!
    ! Written by Nicholas Hine on 12/12/2011.                      !
    !==============================================================!

    use simulation_cell, only: pub_fftbox
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check

    implicit none

    ! Arguments
    type(PROJECTOR_SET), intent(inout) :: proj_set
    type(POINT), intent(in), optional :: kpt
    logical, intent(in), optional :: swap_rc
    real(kind=dp), intent(in), optional :: delta
    integer, intent(in), optional :: cart

    ! Local Variables
    integer :: proj_count, proj_count_on_species
    integer :: shell, am, azim
    integer :: isp
    integer :: ierr      ! error flag
    type(POINT) :: kpt_loc
    logical :: loc_swap_rc, loc_grad
    real(kind=DP) :: proj_norm

    call timer_clock('projectors_init_fftbox_recip',1)

    ! ndmh: check for optional arguments
    if (present(kpt)) then
       kpt_loc = kpt
    else
       kpt_loc%x = 0.0_DP ; kpt_loc%y = 0.0_DP ; kpt_loc%z = 0.0_DP
    end if
    if (present(delta).and.present(cart)) then
       if ((cart<1).or.(cart>3)) then
          call utils_abort('Error in projectors_init_fftbox_recip: Invalid &
               &Cartesian direction supplied')
       end if
       loc_grad = .true.
    else
       loc_grad = .false.
    end if

    if (present(swap_rc)) then
       loc_swap_rc = swap_rc
    else
       loc_swap_rc = .false.
    end if

    ! Count total projectors in this set (including m)
    proj_count = 0
    do isp=1,proj_set%n_proj_species
       proj_count = proj_count + proj_set%species_num_proj(isp)
    end do

    ! Allocate projector FFTboxes
    allocate(proj_set%fftbox_proj_recip(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3,proj_count),stat=ierr)
    call utils_alloc_check('projectors_init_fftbox_recip', &
         'proj_set%fftbox_proj_recip',ierr)

    ! cks: now initialise projectors in fftbox reciprocal representation
    proj_count = 0
    do isp=1,proj_set%n_proj_species

       proj_count_on_species = 0
       do shell=1,proj_set%num_shells(isp)

          am = proj_set%ang_mom(shell,isp)
          do azim=-am,am

             proj_count = proj_count + 1
             proj_count_on_species = proj_count_on_species + 1

             ! Put projector in reciprocal space box
             ! pdh: change to avoid copy in/out of fftbox_proj_recip subarray
             ! ndmh: change to allow use of proj_box_rec by other modules
             if (.not.loc_grad) then
                call projector_in_box_recip( &
                     proj_set%fftbox_proj_recip(:,:,:,proj_count), &
                     proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                     shell,isp), am, azim, pub_fftbox%total_pt1, &
                     pub_fftbox%total_pt2, pub_fftbox%total_pt3, &
                     proj_set%n_rad_pts(isp), proj_set%gmax(isp), kpt_loc)
             else
                call projector_in_box_recip_gradg( &
                     proj_set%fftbox_proj_recip(:,:,:,proj_count), &
                     proj_set%rad_proj_recip(1:proj_set%n_rad_pts(isp), &
                     shell,isp), am, azim, pub_fftbox%total_pt1, &
                     pub_fftbox%total_pt2, pub_fftbox%total_pt3, &
                     proj_set%n_rad_pts(isp), proj_set%gmax(isp), delta, cart)
             end if

             ! ndmh: normalise the projector if required
             ! ndmh: (only used for Hubbard hydrogenic projectors)
             if (proj_set%normalise) then
                proj_norm = sqrt(sum(abs(proj_set%fftbox_proj_recip(:,:,:, &
                     proj_count))**2)/real(pub_fftbox%total_pt1 * &
                     pub_fftbox%total_pt2*pub_fftbox%total_pt3,kind=DP) &
                     / pub_fftbox%weight)
                proj_set%fftbox_proj_recip(:,:,:,proj_count) = &
                     proj_set%fftbox_proj_recip(:,:,:,proj_count) / proj_norm                     
             end if

             ! ndmh: if we passed in swap_rc=T, then swap round real and
             ! ndmh: imaginary parts of projector
             if (loc_swap_rc) then
                proj_set%fftbox_proj_recip(:,:,:,proj_count) = &
                     conjg(proj_set%fftbox_proj_recip(:,:,:,proj_count)) &
                     * cmplx(0.0_DP,1.0_DP)
             end if

             ! Go to next species if we have set up all the projectors on
             ! this atom (even if there are more according to num_shells)
             if (proj_count_on_species>=proj_set%species_num_proj(isp)) exit

          end do
          if (proj_count_on_species>=proj_set%species_num_proj(isp)) exit

       end do

    end do

    call timer_clock('projectors_init_fftbox_recip',2)

  end subroutine projectors_init_fftbox_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_exit_fftbox_recip(proj_set)

    !==============================================================!
    ! This subroutine deallocates a set of projectors in FFTboxes  !
    !--------------------------------------------------------------!
    ! Written by Nicholas Hine on 12/12/2011.                      !
    !==============================================================!

    use simulation_cell, only: pub_fftbox
    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(PROJECTOR_SET), intent(inout) :: proj_set

    ! Local Variables
    integer :: ierr      ! error flag

    ! Allocate projector FFTboxes
    deallocate(proj_set%fftbox_proj_recip,stat=ierr)
    call utils_dealloc_check('projectors_init_fftbox_recip', &
         'proj_set%fftbox_proj_recip',ierr)

  end subroutine projectors_exit_fftbox_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projector_in_box_recip(fftbox_recip, &
       radial_projector,ang_mom,azimuthal_ang_mom,n1,n2,n3, &
       n_rad_pts, g_max, kpt, fine)

    !==============================================================!
    ! This subroutine returns a single projector in an FFTbox      !
    ! in reciprocal space.                                         !
    !--------------------------------------------------------------!
    ! Written by Arash A. Mostofi in 2001 based on a subroutine    !
    ! developed earlier by Chris-Kriton Skylaris.                  !
    ! Modified by Arash A. Mostofi in April 2003 to work with      !
    ! complex-to-complex FFTs.                                     !
    ! Modified by Chris-Kriton Skylaris on 25/1/2004 to work with  !
    ! the pseudo_species type (in the calling subroutines).        !
    ! Extension to f-type projectors by Chris-Kriton Skylaris      !
    ! on 17/08/2007.                                               !
    ! Modified to use pub_cell by Quintin Hill on 15/10/2008.      !
    ! Modified to use an output fftbox rather than only depositing !
    ! to the p_species array by  Nicholas Hine on 15/06/2009.      !
    ! Modified to use sw_real_sph_harm by Quintin Hill on          !
    ! 18/09/2009.                                                  !
    ! Moved to projectors_mod by Nicholas Hine on 19/05/2010.      !
    !==============================================================!

    use comms, only: pub_on_root, comms_abort
    use geometry, only: POINT, OPERATOR(*), OPERATOR(+), geometry_magnitude
    use services, only: services_1d_interpolation
    use simulation_cell, only: pub_cell
    use spherical_wave, only: sw_real_sph_harm

    implicit none

    ! Arguments
    integer, intent(in) :: n1,n2,n3
    complex(kind=DP), intent(out) :: fftbox_recip(n1,n2,n3)
    integer, intent(in) :: ang_mom,azimuthal_ang_mom
    integer, intent(in) :: n_rad_pts
    real(kind =DP), intent(in) :: g_max
    real(kind=DP), intent(in)  :: radial_projector(n_rad_pts)
    type(POINT), intent(in) :: kpt
    character(len=*), optional, intent(in) :: fine

    ! Local Variables
    integer     :: row1,row2,row3
    integer     :: K3,K2,K1
    type(POINT) :: g_vector3,g_vector2,g_vector1
    type(POINT) :: pairvec3,pairvec2,pairvec1
    real(kind=DP) :: g_length
    complex(kind=DP) :: phase_fac
    real(kind=DP) :: projector_value,harmonic_factor
    real(kind=DP), parameter :: four_pi=4.0_DP*PI

    pairvec3=( real(pub_cell%total_pt3,DP)/real(n3,DP) ) * pub_cell%b3
    pairvec2=( real(pub_cell%total_pt2,DP)/real(n2,DP) ) * pub_cell%b2
    pairvec1=( real(pub_cell%total_pt1,DP)/real(n1,DP) ) * pub_cell%b1

    if (present(fine)) then
       pairvec3=2.0_DP*pairvec3
       pairvec2=2.0_DP*pairvec2
       pairvec1=2.0_DP*pairvec1
    endif

    ! calculate the phase factor (-i)^l
    select case (ang_mom)
    case (0)
       phase_fac = cmplx(1.0_DP,0.0_DP)
    case (1)
       phase_fac = cmplx(0.0_DP,-1.0_DP)
    case (2)
       phase_fac = cmplx(-1.0_DP,0.0_DP)
    case (3)
       phase_fac = cmplx(0.0_DP,1.0_DP)
    case default
       if (pub_on_root) then
          write(stdout,*) 'ERROR in projector_in_box_recip: Angular momentum &
               &too high'
          write(stdout,*) 'Only projectors up to f-type symmetry &
               &are currently supported.'
       end if
       call comms_abort
       ! qoh: Initialise to avoid compiler warning
       phase_fac = cmplx(0.0_DP,0.0_DP)
    end select

    ! loop over the points of the reciprocal FFT box
    do row3=0,n3-1

       k3=row3
       if ( row3.gt.(n3/2) ) k3=row3-n3

       g_vector3 = real(k3,kind=DP)*pairvec3 + kpt

       do row2=0,n2-1

          k2=row2
          if ( row2.gt.(n2/2) ) k2=row2-n2

          g_vector2 = real(k2,kind=DP)*pairvec2 + g_vector3

          do row1=0,n1-1

             k1=row1
             if ( row1.gt.(n1/2) ) k1=row1-n1

             g_vector1 = real(k1,kind=DP)*pairvec1 + g_vector2

             g_length=geometry_magnitude(g_vector1)

             ! cks: interpolation of radial projector
             projector_value=services_1d_interpolation(  &
                  radial_projector(:), n_rad_pts, &
                  (g_length/g_max)*real(n_rad_pts-1, kind=DP), ang_mom)

             ! multiply value of projector at each g-point by the appropriate
             ! phase factor and the appropriate real spherical harmonic factor
             harmonic_factor = four_pi * sw_real_sph_harm(g_vector1%x, &
                  g_vector1%y, g_vector1%z, g_length, ang_mom, azimuthal_ang_mom)

             fftbox_recip(row1+1,row2+1,row3+1) &
                  = projector_value * harmonic_factor * phase_fac

          enddo
       enddo
    enddo

  end subroutine projector_in_box_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projector_in_box_recip_gradg(fftbox_recip, &
       radial_projector,ang_mom,azimuthal_ang_mom,n1,n2,n3, &
       n_rad_pts, g_max, fin_diff_shift, cart, fine)

    !==============================================================!
    ! This subroutine calculates the gradient of the projectors    !
    ! with respect to G in reciprocal space.                       !
    !--------------------------------------------------------------!
    ! Written by Laura Ratcliff in March 2011 based on             !
    ! projector_in_box_recip.                                      !
    !==============================================================!

    use comms, only: pub_on_root, comms_abort
    use geometry, only: POINT, OPERATOR(*), OPERATOR(+), geometry_magnitude
    use services, only: services_1d_interpolation
    use simulation_cell, only: pub_cell
    use spherical_wave, only: sw_real_sph_harm, sw_grad_real_sph_harm !lr408

    implicit none

    ! Arguments
    integer, intent(in) :: n1,n2,n3
    complex(kind=DP), intent(out) :: fftbox_recip(n1,n2,n3)
    integer, intent(in) :: ang_mom,azimuthal_ang_mom
    integer, intent(in) :: n_rad_pts
    real(kind=DP), intent(in) :: g_max
    real(kind=DP), intent(in)  :: radial_projector(n_rad_pts)
    real(kind=dp), intent(in) :: fin_diff_shift
    integer, intent(in) :: cart
    character(len=*), optional, intent(in) :: fine

    ! Local Variables
    integer     :: row1,row2,row3
    integer     :: K3,K2,K1
    type(POINT) :: g_vector3,g_vector2,g_vector1
    type(POINT) :: pairvec3,pairvec2,pairvec1
    real(kind=DP) :: g_length
    complex(kind=DP) :: phase_fac
    real(kind=DP) :: projector_value,harmonic_factor
    real(kind=DP), parameter :: four_pi=4.0_DP*PI

    real(kind=DP) :: grad_sph_harm(3), grad_proj(3)
    real(kind=DP) :: g_length2, projector_value2, grad_projector_value

    pairvec3=( real(pub_cell%total_pt3,DP)/real(n3,DP) ) * pub_cell%b3
    pairvec2=( real(pub_cell%total_pt2,DP)/real(n2,DP) ) * pub_cell%b2
    pairvec1=( real(pub_cell%total_pt1,DP)/real(n1,DP) ) * pub_cell%b1

    if (present(fine)) then
       pairvec3=2.0_DP*pairvec3
       pairvec2=2.0_DP*pairvec2
       pairvec1=2.0_DP*pairvec1
    endif

    ! calculate the phase factor (-i)^l
    select case (ang_mom)
    case (0)
       phase_fac = cmplx(1.0_DP,0.0_DP)
    case (1)
       phase_fac = cmplx(0.0_DP,-1.0_DP)
    case (2)
       phase_fac = cmplx(-1.0_DP,0.0_DP)
    case (3)
       phase_fac = cmplx(0.0_DP,1.0_DP)
    case default
       if (pub_on_root) write(stdout,*) 'Angular momentum too high too high in', &
            ' projector_in_box_recip_gradg. Only projectors up to f-type symmetry ',&
            ' are currently supported.'
       call comms_abort
       ! qoh: Initialise to avoid compiler warning
       phase_fac = cmplx(0.0_DP,0.0_DP)
    end select

    ! loop over the points of the reciprocal FFT box
    do row3=0,n3-1

       k3=row3
       if ( row3.gt.(n3/2) ) k3=row3-n3

       g_vector3 = real(k3,kind=DP)*pairvec3 !+ kpt

       do row2=0,n2-1

          k2=row2
          if ( row2.gt.(n2/2) ) k2=row2-n2

          g_vector2 = real(k2,kind=DP)*pairvec2 + g_vector3

          do row1=0,n1-1

             k1=row1
             if ( row1.gt.(n1/2) ) k1=row1-n1

             g_vector1 = real(k1,kind=DP)*pairvec1 + g_vector2

             g_length=geometry_magnitude(g_vector1)

             ! cks: interpolation of radial projector
             projector_value=services_1d_interpolation(  &
                  radial_projector(:), n_rad_pts, &
                  (g_length/g_max)*real(n_rad_pts-1, kind=DP), ang_mom)


             ! lr408: interpolate projector at a second point for finite diff
             g_length2 = g_length + fin_diff_shift

             projector_value2=services_1d_interpolation(  &
                  radial_projector(:), n_rad_pts, &
                  (g_length2/g_max)*real(n_rad_pts-1, kind=DP), ang_mom)

             grad_projector_value = (projector_value2 - projector_value) / fin_diff_shift

             harmonic_factor = four_pi * sw_real_sph_harm(g_vector1%x, &
                  g_vector1%y, g_vector1%z, g_length, ang_mom, azimuthal_ang_mom)


             if (g_length > 1d-10) then
                grad_proj(1) = grad_projector_value * harmonic_factor &
                     * g_vector1%x / g_length
                grad_proj(2) = grad_projector_value * harmonic_factor &
                     * g_vector1%y / g_length
                grad_proj(3) = grad_projector_value * harmonic_factor &
                     * g_vector1%z / g_length
             else
                if (ang_mom == 1) then  !special value for L=1: lim_{q->0} g(q)/q
                   projector_value = grad_projector_value
                   grad_proj(:) = 0.0_dp                
                else
                   projector_value = 0.0_dp
                   grad_proj(:) = 0.0_dp
                end if
             end if

             ! multiply value of projector at each g-point by the appropriate
             ! phase factor (n.b. real part and imaginary part stored as
             ! separate consecutive elements in the x-direction of the array),
             ! and the appropriate real spherical harmonic factor.
             call sw_grad_real_sph_harm(grad_sph_harm, g_vector1%x, &
                  g_vector1%y, g_vector1%z, g_length, ang_mom, azimuthal_ang_mom)

             ! grad_proj * sph_harm
             fftbox_recip(row1+1,row2+1,row3+1) &
                  = grad_proj(cart) * phase_fac

             ! + proj * grad_sph_harm
             fftbox_recip(row1+1,row2+1,row3+1) &
                  = fftbox_recip(row1+1,row2+1,row3+1) + four_pi &
                  * projector_value * grad_sph_harm(cart) * phase_fac

          enddo
       enddo
    enddo

  end subroutine projector_in_box_recip_gradg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_create_real(proj_basis,proj_set,projs_on_grid)

    !===============================================================!
    ! This subroutine allocates space for and calculates a set of   !
    ! projectors in real space on PPDs                              !
    !---------------------------------------------------------------!
    ! Written by Nicholas Hine in March 2011.                       !
    !===============================================================!

    use basis, only: basis_phase_on_fftbox_recip, basis_point_wrt_box, &
         basis_start_of_box_wrt_cell, basis_extract_function_from_box, &
         basis_find_function_wrt_box
    use comms, only: comms_abort, comms_bcast, comms_reduce, &
         pub_on_root, pub_my_node_id
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_brappd_ketfftbox
    use parallel_strategy, only: pub_first_atom_on_node
    use simulation_cell, only: pub_cell, pub_fftbox
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(PROJECTOR_SET), intent(inout) :: proj_set
    real(kind=DP), intent(out), optional :: projs_on_grid(proj_basis%n_ppds* &
         pub_cell%n_pts)

    ! Local Variables
    integer :: local_proj
    integer :: global_proj
    integer :: atom_of_proj
    integer :: loc_atom_of_proj
    integer :: species_number
    integer :: proj_count, proj_count_atom
    integer :: ierr             ! Error flag
    integer :: n1,n2,n3,ld1,ld2
    integer :: proj_start(1:3)    ! Start of projector function in FFTbox
    integer :: fftbox_start(1:3)  ! Start of FFTbox in cell
    ! {P}rojector {C}entre wrt fft {B}ox in terms of number of {G}ridpoints
    ! in each lattice vector direction. Not necessarily integer values.
    real(kind=DP) :: pcbg1,pcbg2,pcbg3
    real(kind=DP), allocatable :: projector_fftbox(:,:,:)
    complex(kind=DP), allocatable :: proj_complex(:,:,:)   ! Workspace

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &projectors_create_real'
#endif

    ! Start timer
    call timer_clock('projectors_create_real',1)

    ! Initialise shorthand variables
    ld1 = pub_fftbox%total_ld1
    ld2 = pub_fftbox%total_ld2
    n1 = pub_fftbox%total_pt1
    n2 = pub_fftbox%total_pt2
    n3 = pub_fftbox%total_pt3

    ! Allocate internal storage in projector set for ppds of projectors
    ! (only if we did not provide an array to hold the result)
    if (.not.present(projs_on_grid)) then
       allocate(proj_set%projs_on_grid(proj_basis%n_ppds*pub_cell%n_pts), &
            stat=ierr)
       call utils_alloc_check('projectors_create_real', &
            'proj_set%projs_on_grid',ierr)
    end if

    ! Allocate workspace
    allocate(proj_complex(ld1, ld2, n3), stat=ierr)
    call utils_alloc_check('projectors_create_real','proj_complex',ierr)
    allocate(projector_fftbox(ld1,ld2,n3),stat=ierr)
    call utils_alloc_check('projectors_create_real','projector_fftbox',ierr)

    do local_proj=1,proj_basis%node_num

       ! Find information about this projector
       global_proj = local_proj + proj_basis%first_on_node(pub_my_node_id) - 1
       atom_of_proj = proj_basis%atom_of_func(global_proj)
       loc_atom_of_proj = atom_of_proj - &
            pub_first_atom_on_node(pub_my_node_id) + 1
       species_number = proj_set%proj_species(atom_of_proj)
       proj_count_atom = global_proj - &
            proj_basis%first_on_atom(atom_of_proj) + 1
       proj_count = proj_set%species_first_proj(species_number) + &
            proj_count_atom - 1

       ! Centre of projector wrt fftbox in terms of grid spacings
       call basis_point_wrt_box(pcbg1, pcbg2, pcbg3, &  ! output
            proj_set%proj_centre(atom_of_proj), n1, n2, n3)

       ! Start of fftbox wrt cell in terms of grid-point number
      call basis_start_of_box_wrt_cell(fftbox_start(1), &
           fftbox_start(2),fftbox_start(3), &
           proj_set%proj_centre(atom_of_proj),pcbg1,pcbg2,pcbg3)

       ! Find position of tightbox start wrt fftbox
       call basis_find_function_wrt_box(proj_start(1), proj_start(2), &
            proj_start(3), fftbox_start(1),fftbox_start(2),fftbox_start(3), &
            proj_basis%tight_boxes(local_proj))

       ! Apply phase factor e^-(ik.proj_centre_wrt_box) to projector
       call basis_phase_on_fftbox_recip(proj_complex, &
            proj_set%fftbox_proj_recip(:,:,:,proj_count),n1,n2,n3,&
            ld1,ld2,-pcbg1,-pcbg2,-pcbg3)

       ! g=0 element must be real
       proj_complex(1,1,1) &
            = cmplx(real(proj_complex(1,1,1),kind=DP),0.0_DP,kind=DP)

       ! Fourier transform to real space:
       call fourier_apply_box('Coarse','Backward',proj_complex)

       ! Put real part into projector_box
       projector_fftbox(1:n1,1:n2,1:n3) &
            = real(proj_complex(1:n1,1:n2,1:n3),kind=DP)/pub_cell%weight

       ! Extract projector to ppds
       if (present(projs_on_grid)) then
          call basis_extract_function_from_box(projs_on_grid(:), &
               ld1,ld2,n3,projector_fftbox,proj_basis%spheres(local_proj), &
               proj_basis%tight_boxes(local_proj),proj_start(1),proj_start(2), &
               proj_start(3),proj_basis%spheres(local_proj)%offset)
       else
          call basis_extract_function_from_box(proj_set%projs_on_grid(:), &
               ld1,ld2,n3,projector_fftbox,proj_basis%spheres(local_proj), &
               proj_basis%tight_boxes(local_proj),proj_start(1),proj_start(2), &
               proj_start(3),proj_basis%spheres(local_proj)%offset)
       end if

    end do

    ! Deallocate workspace
    deallocate(projector_fftbox,stat=ierr)
    call utils_dealloc_check('projectors_create_real','projector_fftbox',ierr)
    deallocate(proj_complex,stat=ierr)
    call utils_dealloc_check('projectors_create_real','proj_complex',ierr)

    ! Stop timer
    call timer_clock('projectors_create_real',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &projectors_create_real'
#endif

  end subroutine projectors_create_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_func_ovlp_box(sp_overlap, &
       ngwfs_on_grid,ngwf_basis,proj_basis,proj_set,kshift,delta,cart,swap_rc)

    !===============================================================!
    ! This subroutine calculates the ngwf-projector overlap matrix  !
    ! and stores it in a SPAM3 matrix.                              !
    !---------------------------------------------------------------!
    ! Written by Nicholas Hine in June 2009.                        !
    ! Some code from reused from the old version of                 !
    ! pseudopotentials_sp_ovlp_box written by Arash A. Mostofi in   !
    ! January 2004 and modified by Chris-Kriton Skylaris, Arash     !
    ! Mostofi and Nicholas Hine in 2004-2008.                       !
    ! Moved to projectors_mod by Nicholas Hine on 19/05/2010.       !
    ! Optional k-vector added as input for calculating overlap at a !
    ! different k-point by Laura Ratcliff 17/2/11.                  !
    !===============================================================!

    use basis, only: basis_phase_on_fftbox_recip, basis_point_wrt_box, &
         basis_start_of_box_wrt_cell
    use comms, only: comms_abort, comms_bcast, comms_reduce, &
         pub_on_root, pub_my_node_id
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_brappd_ketfftbox
    use parallel_strategy, only: pub_first_atom_on_node
    use rundat, only: locpot_int_batch_size
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: sp_overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(PROJECTOR_SET), intent(inout) :: proj_set
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(POINT), optional, intent(in) :: kshift
    real(kind=DP), optional, intent(in) :: delta
    integer, optional, intent(in) :: cart
    logical, optional, intent(in) :: swap_rc

    ! Local Variables
    integer :: n_batches
    integer :: batch_size
    integer :: batch_count
    integer :: local_start
    integer :: local_end
    integer :: local_proj
    integer :: global_proj
    integer :: batch_proj
    integer :: atom_of_proj
    integer :: loc_atom_of_proj
    integer :: species_number
    integer :: proj_count, proj_count_atom
    integer :: ierr             ! Error flag
    integer :: idx_len          ! Length of sparse index
    integer :: n1,n2,n3,ld1,ld2
    ! {P}rojector {C}entre wrt fft {B}ox in terms of number of {G}ridpoints
    ! in each lattice vector direction. Not necessarily integer values.
    real(kind=DP) :: pcbg1,pcbg2,pcbg3
    real(kind=DP), allocatable :: projector_fftbox(:,:,:,:)
    integer, allocatable :: fftbox_start(:,:)              ! Ket start positions
    integer, allocatable :: sp_overlap_idx(:)              ! Sparse index
    real(kind=DP), allocatable :: bra_on_grid_buffer(:)    ! Workspace
    real(kind=DP), allocatable :: ket_on_bragrid_buffer(:) ! Workspace
    complex(kind=DP), allocatable :: proj_complex(:,:,:)   ! Workspace
    logical :: loc_swap_rc ! lr408

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &projectors_func_ovlp_box'
#endif

    ! Start timer
    call timer_clock('projectors_func_ovlp_box',1)

    if (present(swap_rc)) then
       loc_swap_rc = swap_rc
    else
       loc_swap_rc = .false.
    end if

    ! Initialise shorthand variables
    ld1 = pub_fftbox%total_ld1
    ld2 = pub_fftbox%total_ld2
    n1 = pub_fftbox%total_pt1
    n2 = pub_fftbox%total_pt2
    n3 = pub_fftbox%total_pt3

    ! Obtain index of sp_overlap
    idx_len = sparse_index_length(sp_overlap)
    allocate(sp_overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('projectors_func_ovlp_box','sp_overlap_idx',ierr)
    call sparse_generate_index(sp_overlap_idx,sp_overlap)

    ! Set batch size
    batch_size = locpot_int_batch_size

    ! Allocate workspace
    allocate(proj_complex(ld1, ld2, n3), stat=ierr)
    call utils_alloc_check('projectors_func_ovlp_box', &
         'proj_complex',ierr)
    allocate(projector_fftbox(ld1,ld2,n3,batch_size),stat=ierr)
    call utils_alloc_check('projectors_func_ovlp_box', &
         'projector_fftbox',ierr)
    allocate(fftbox_start(3,batch_size),stat=ierr)
    call utils_alloc_check('projectors_func_ovlp_box', &
         'fftbox_start',ierr)
    allocate(bra_on_grid_buffer(ngwf_basis%func_on_grid_buffer_size),stat=ierr)
    call utils_alloc_check('projectors_func_ovlp_box', &
         'bra_on_grid_buffer',ierr)
    allocate(ket_on_bragrid_buffer( &
         ngwf_basis%max_n_ppds_sphere*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('projectors_func_ovlp_box', &
         'ket_on_bragrid_buffer',ierr)

    call projectors_init_fftbox_recip(proj_set,kshift,delta,cart)

    ! Number of row-steps per row-block
    n_batches = proj_basis%max_on_node / batch_size
    if (mod(proj_basis%max_on_node,batch_size) > 0) n_batches = n_batches + 1

    local_start = 1
    do batch_count=1,n_batches
       local_end = min(local_start+batch_size-1,proj_basis%node_num)
       do local_proj=local_start,local_end

          batch_proj = local_proj - local_start + 1

          ! Find information about this projector
          global_proj = local_proj + proj_basis%first_on_node(pub_my_node_id) &
               - 1
          atom_of_proj = proj_basis%atom_of_func(global_proj)
          loc_atom_of_proj = atom_of_proj - &
               pub_first_atom_on_node(pub_my_node_id) + 1
          species_number = proj_set%proj_species(atom_of_proj)
          proj_count_atom = global_proj - &
               proj_basis%first_on_atom(atom_of_proj) + 1
          proj_count = proj_set%species_first_proj(species_number) + &
               proj_count_atom - 1

          ! Centre of projector wrt fftbox in terms of grid spacings
          call basis_point_wrt_box(pcbg1, pcbg2, pcbg3, &  ! output
               proj_set%proj_centre(atom_of_proj), n1, n2, n3)

          ! Start of fftbox wrt cell in terms of grid-point number
          call basis_start_of_box_wrt_cell(fftbox_start(1,batch_proj), &
               fftbox_start(2,batch_proj),fftbox_start(3,batch_proj), &
               proj_set%proj_centre(atom_of_proj),pcbg1,pcbg2,pcbg3)

          ! Apply phase factor e^-(ik.proj_centre_wrt_box) to projector
          call basis_phase_on_fftbox_recip(proj_complex, &
               proj_set%fftbox_proj_recip(:,:,:,proj_count),n1,n2,n3,&
               ld1,ld2,-pcbg1,-pcbg2,-pcbg3)

          ! g=0 element must be real
          proj_complex(1,1,1) &
               = cmplx(real(proj_complex(1,1,1),kind=DP),0.0_DP,kind=DP)

          ! Fourier transform to real space:
          call fourier_apply_box('Coarse','Backward',proj_complex)

          ! Put real part into projector_box
          projector_fftbox(1:n1,1:n2,1:n3,batch_proj) &
               = real(proj_complex(1:n1,1:n2,1:n3),kind=DP)/pub_cell%weight

          if (loc_swap_rc) then
             ! Put imag part into projector_box
             projector_fftbox(1:n1,1:n2,1:n3,batch_proj) &
                  = aimag(proj_complex(1:n1,1:n2,1:n3))/pub_cell%weight
          end if

       end do

       ! Calculate overlap integrals
       call integrals_brappd_ketfftbox(sp_overlap,  &                   ! inout
            ngwfs_on_grid, ngwf_basis, &                                ! input
            projector_fftbox, fftbox_start, batch_size,  &              ! input
            local_start, local_end, idx_len, sp_overlap_idx, 'FULL', &  ! input
            bra_on_grid_buffer,ket_on_bragrid_buffer)                   ! work

       local_start = local_start + batch_size

    end do

    call projectors_exit_fftbox_recip(proj_set)

    ! Deallocate workspace
    deallocate(ket_on_bragrid_buffer,stat=ierr)
    call utils_dealloc_check('projectors_func_ovlp_box', &
         'ket_on_bragrid_buffer',ierr)
    deallocate(bra_on_grid_buffer,stat=ierr)
    call utils_dealloc_check('projectors_func_ovlp_box', &
         'bra_on_grid_buffer',ierr)
    deallocate(fftbox_start,stat=ierr)
    call utils_dealloc_check('projectors_func_ovlp_box', &
         'fftbox_start',ierr)
    deallocate(projector_fftbox,stat=ierr)
    call utils_dealloc_check('projectors_func_ovlp_box', &
         'projector_fftbox',ierr)
    deallocate(proj_complex,stat=ierr)
    call utils_dealloc_check('projectors_func_ovlp_box', &
         'proj_complex',ierr)
    deallocate(sp_overlap_idx,stat=ierr)
    call utils_dealloc_check('projectors_func_ovlp_box', &
         'sp_overlap_idx',ierr)

    ! Stop timer
    call timer_clock('projectors_func_ovlp_box',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &projectors_func_ovlp_box'
#endif

  end subroutine projectors_func_ovlp_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_func_grad_ovlp_box(siGp_overlap, &
       ngwfs_on_grid,ngwf_basis,proj_basis,proj_set,cart)

    !========================================================================!
    ! This subroutine calculates the overlap matrix between the NGWFs and    !
    ! the first derivatives of the projectors in each of the Cartesian       !
    ! directions, and stores it in a SPAM3 matrix.                           !
    !------------------------------------------------------------------------!
    ! Arguments                                                              !
    ! 1)  siGp_overlap  : output : <ngwf|iG*proj> overlap matrix             !
    ! 2)  ngwfs_on_grid : input  : ngwf values on the grid                   !
    ! 3)  ngwfs_basis   : input  : function basis type for the NGWFs         !
    !------------------------------------------------------------------------!
    ! Written by Nicholas Hine in July 2009.                                 !
    ! Some code reused from the old version of                               !
    ! pseudopotentials_sp_ovlp_box written by Arash A. Mostofi in            !
    ! January 2004 and modified by Chris-Kriton Skylaris, Arash              !
    ! Mostofi and Nicholas Hine in 2004-2008.                                !
    ! Moved to projectors_mod by Nicholas Hine on 19/05/2010.                !
    !========================================================================!

    use basis, only: basis_phase_on_fftbox_recip, basis_point_wrt_box, &
         basis_start_of_box_wrt_cell
    use comms, only: comms_abort, comms_bcast, comms_reduce, &
         pub_on_root, pub_my_node_id
    use constants, only: cmplx_i
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(*), operator(+)
    use integrals, only: integrals_brappd_ketfftbox
    use parallel_strategy, only: pub_first_atom_on_node
    use rundat, only: locpot_int_batch_size
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: siGp_overlap
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(PROJECTOR_SET), intent(inout) :: proj_set
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    integer, intent(in) :: cart

    ! Local Variables
    integer :: n_batches
    integer :: batch_size
    integer :: batch_count
    integer :: local_start
    integer :: local_end
    integer :: local_proj
    integer :: global_proj
    integer :: batch_proj
    integer :: atom_of_proj
    integer :: loc_atom_of_proj
    integer :: species_number
    integer :: proj_count, proj_count_atom
    integer :: ierr             ! Error flag
    integer :: idx_len          ! Length of sparse index
    integer :: n1,n2,n3,ld1,ld2
    ! {P}rojector {C}entre wrt fft {B}ox in terms of number of {G}ridpoints
    ! in each lattice vector direction. Not necessarily integer values.
    real(kind=DP) :: pcbg1,pcbg2,pcbg3
    real(kind=DP), allocatable :: projector_fftbox(:,:,:,:)
    integer, allocatable :: fftbox_start(:,:)              ! Ket start positions
    integer, allocatable :: sp_overlap_idx(:)              ! Sparse index
    real(kind=DP), allocatable :: bra_on_grid_buffer(:)    ! Workspace
    real(kind=DP), allocatable :: ket_on_bragrid_buffer(:) ! Workspace
    complex(kind=DP), allocatable :: proj_complex(:,:,:)   ! Workspace

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &projectors_func_grad_ovlp_box'
#endif

    ! Start timer
    call timer_clock('projectors_func_grad_ovlp_box',1)

    ! Initialise shorthand variables
    ld1 = pub_fftbox%total_ld1
    ld2 = pub_fftbox%total_ld2
    n1 = pub_fftbox%total_pt1
    n2 = pub_fftbox%total_pt2
    n3 = pub_fftbox%total_pt3

    ! Obtain index of siGp_overlap
    idx_len = sparse_index_length(siGp_overlap)
    allocate(sp_overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('projectors_func_grad_ovlp_box','sp_overlap_idx',ierr)
    call sparse_generate_index(sp_overlap_idx,siGp_overlap)

    ! Set batch size
    batch_size = locpot_int_batch_size

    ! Allocate workspace
    allocate(proj_complex(ld1, ld2, n3), stat=ierr)
    call utils_alloc_check('projectors_func_grad_ovlp_box', &
         'proj_complex',ierr)
    allocate(projector_fftbox(ld1,ld2,n3,batch_size),stat=ierr)
    call utils_alloc_check('projectors_func_grad_ovlp_box', &
         'projector_fftbox',ierr)
    allocate(fftbox_start(3,batch_size),stat=ierr)
    call utils_alloc_check('projectors_func_grad_ovlp_box', &
         'fftbox_start',ierr)
    allocate(bra_on_grid_buffer(ngwf_basis%func_on_grid_buffer_size),stat=ierr)
    call utils_alloc_check('projectors_func_grad_ovlp_box', &
         'bra_on_grid_buffer',ierr)
    allocate(ket_on_bragrid_buffer( &
         ngwf_basis%max_n_ppds_sphere*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('projectors_func_grad_ovlp_box', &
         'ket_on_bragrid_buffer',ierr)

    call projectors_init_fftbox_recip(proj_set)

    ! Number of row-steps per row-block
    n_batches = proj_basis%max_on_node / batch_size
    if (mod(proj_basis%max_on_node,batch_size) > 0) n_batches = n_batches + 1

    local_start = 1
    do batch_count=1,n_batches
       local_end = min(local_start+batch_size-1,proj_basis%node_num)
       do local_proj=local_start,local_end

          batch_proj = local_proj - local_start + 1

          ! Find information about this projector
          global_proj = local_proj + proj_basis%first_on_node(pub_my_node_id) - 1
          atom_of_proj = proj_basis%atom_of_func(global_proj)
          loc_atom_of_proj = atom_of_proj - &
               pub_first_atom_on_node(pub_my_node_id) + 1
          species_number = proj_set%proj_species(atom_of_proj)
          proj_count_atom = global_proj - &
               proj_basis%first_on_atom(atom_of_proj) + 1
          proj_count = proj_set%species_first_proj(species_number) + &
               proj_count_atom - 1

          ! Centre of projector wrt fftbox in terms of grid spacings
          call basis_point_wrt_box(pcbg1, pcbg2, pcbg3, &  ! output
               proj_set%proj_centre(atom_of_proj), n1, n2, n3)

          ! Start of fftbox wrt cell in terms of grid-point number
          call basis_start_of_box_wrt_cell(fftbox_start(1,batch_proj), &
               fftbox_start(2,batch_proj), fftbox_start(3,batch_proj), &
               proj_set%proj_centre(atom_of_proj), pcbg1, pcbg2, &
               pcbg3)

          ! Apply phase factor e^-(ik.proj_centre_wrt_box) to projector
          call basis_phase_on_fftbox_recip(proj_complex, &
               proj_set%fftbox_proj_recip(:,:,:,proj_count),n1,n2,n3, &
               ld1,ld2,-pcbg1,-pcbg2,-pcbg3)

          ! Multiply projector (workspace) by iG and store as proj_complex
          proj_complex(:,:,:) = &
               proj_complex(:,:,:)*pub_fftbox%recip_grid(cart,:,:,:)*cmplx_i

          ! g=0 element must be real
          proj_complex(1,1,1) = cmplx( &
               real(proj_complex(1,1,1),kind=DP), 0.0_DP, kind=DP )

          ! Fourier transform to real space:
          ! ndmh: explicitly in-place transform
          call fourier_apply_box('Coarse','Backward',proj_complex)

          ! Put real part into projector_box
          projector_fftbox(1:n1,1:n2,1:n3,batch_proj) &
               = real(proj_complex(1:n1,1:n2,1:n3),kind=DP)/pub_cell%weight

       end do

       ! Calculate overlap integrals
       call integrals_brappd_ketfftbox(siGp_overlap, &           ! inout
            ngwfs_on_grid, ngwf_basis, &                         ! input
            projector_fftbox(:,:,:,:), fftbox_start, &           ! input
            batch_size, local_start, local_end, idx_len, &       ! input
            sp_overlap_idx, 'FULL', bra_on_grid_buffer, &        ! input
            ket_on_bragrid_buffer)

       local_start = local_start + batch_size

    end do

    call projectors_exit_fftbox_recip(proj_set)

    ! Deallocate workspace
    deallocate(ket_on_bragrid_buffer,stat=ierr)
    call utils_dealloc_check('projectors_func_grad_ovlp_box', &
         'ket_on_bragrid_buffer',ierr)
    deallocate(bra_on_grid_buffer,stat=ierr)
    call utils_dealloc_check('projectors_func_grad_ovlp_box', &
         'bra_on_grid_buffer',ierr)
    deallocate(fftbox_start,stat=ierr)
    call utils_dealloc_check('projectors_func_grad_ovlp_box', &
         'fftbox_start',ierr)
    deallocate(projector_fftbox,stat=ierr)
    call utils_dealloc_check('projectors_func_grad_ovlp_box', &
         'projector_fftbox',ierr)
    deallocate(proj_complex,stat=ierr)
    call utils_dealloc_check('projectors_func_grad_ovlp_box', &
         'proj_complex',ierr)
    deallocate(sp_overlap_idx,stat=ierr)
    call utils_dealloc_check('projectors_func_grad_ovlp_box', &
         'sp_overlap_idx',ierr)

    ! Stop timer
    call timer_clock('projectors_func_grad_ovlp_box',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &projectors_func_grad_ovlp_box'
#endif

  end subroutine projectors_func_grad_ovlp_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_func_pos_ovlp_box(sp_overlap, &
       ngwfs_on_grid,ngwf_basis,proj_basis,proj_set)

    !===============================================================!
    ! This subroutine calculates the overlap matrix between a set   !
    ! of NGWFs, and a set of projectors, with the position operator !
    ! in between the two, and stores the result in a SPAM3 matrix.  !
    !---------------------------------------------------------------!
    ! Written by Laura Ratcliff in November 2011, based on          !
    ! the routine projectors_func_ovlp_box by Nicholas Hine.        !
    !===============================================================!

    use basis, only: basis_phase_on_fftbox_recip, basis_point_wrt_box, &
         basis_start_of_box_wrt_cell
    use comms, only: comms_abort, comms_bcast, comms_reduce, &
         pub_on_root, pub_my_node_id
    use constants, only: cmplx_i
    use fourier, only: fourier_apply_box, fourier_filter, fourier_interpolate
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(*), operator(+)
    use integrals, only: integrals_brappd_ketfftbox
    use parallel_strategy, only: pub_first_atom_on_node
    use rundat, only: locpot_int_batch_size
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(SPAM3),intent(inout) :: sp_overlap(3)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(PROJECTOR_SET), intent(inout) :: proj_set
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)

    ! Local Variables
    type(POINT) :: r1,r2,r3 ! lr408
    type(POINT) :: a1,a2,a3 ! lr408
    integer :: i1,i2,i3,count

    integer :: n_batches
    integer :: batch_size
    integer :: batch_count
    integer :: local_start
    integer :: local_end
    integer :: local_proj
    integer :: global_proj
    integer :: batch_proj, old_batch_proj
    integer :: atom_of_proj
    integer :: loc_atom_of_proj
    integer :: species_number
    integer :: proj_count, proj_count_atom
    integer :: ierr             ! Error flag
    integer :: idx_len          ! Length of sparse index
    integer :: n1,n2,n3,ld1,ld2
    integer :: ld1_dbl,ld2_dbl ! lr408
    ! {P}rojector {C}entre wrt fft {B}ox in terms of number of {G}ridpoints
    ! in each lattice vector direction. Not necessarily integer values.
    integer:: xyz
    real(kind=DP) :: pcbg1,pcbg2,pcbg3
    real(kind=DP), allocatable :: projector_fftbox(:,:,:,:)
    integer, allocatable :: fftbox_start(:,:)              ! Ket start positions
    integer, allocatable :: sp_overlap_idx(:)              ! Sparse index
    real(kind=DP), allocatable :: bra_on_grid_buffer(:)    ! Workspace
    real(kind=DP), allocatable :: ket_on_bragrid_buffer(:) ! Workspace
    complex(kind=DP), allocatable :: proj_complex(:,:,:)   ! Workspace
    real(kind=DP), allocatable :: old_proj(:,:,:)   ! Workspace 
    real(kind=DP), allocatable :: new_proj(:,:,:)   ! Workspace 
    real(kind=DP), dimension(:,:,:,:), allocatable :: r_op
    real(kind=DP), dimension(:,:,:), allocatable :: fftbox1_dbl,fftbox2_dbl

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &projectors_func_pos_ovlp_box'
#endif

    ! Start timer
    call timer_clock('projectors_func_pos_ovlp_box',1)

    ! Initialise shorthand variables
    ld1 = pub_fftbox%total_ld1
    ld2 = pub_fftbox%total_ld2
    n1 = pub_fftbox%total_pt1
    n2 = pub_fftbox%total_pt2
    n3 = pub_fftbox%total_pt3

    ld1_dbl = pub_fftbox%total_ld1_dbl ! lr408
    ld2_dbl = pub_fftbox%total_ld2_dbl ! lr408

    ! Obtain index of sp_overlap
    idx_len = sparse_index_length(sp_overlap(1))
    allocate(sp_overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box','sp_overlap_idx',ierr)
    call sparse_generate_index(sp_overlap_idx,sp_overlap(1))

    ! Set batch size
    batch_size = locpot_int_batch_size

    ! Allocate workspace
    allocate(proj_complex(ld1, ld2, n3), stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box', &
         'proj_complex',ierr)
    allocate(projector_fftbox(ld1,ld2,n3,batch_size),stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box', &
         'projector_fftbox',ierr)
    allocate(fftbox_start(3,batch_size),stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box', &
         'fftbox_start',ierr)
    allocate(bra_on_grid_buffer(ngwf_basis%func_on_grid_buffer_size),stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box', &
         'bra_on_grid_buffer',ierr)
    allocate(ket_on_bragrid_buffer( &
         ngwf_basis%max_n_ppds_sphere*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box', &
         'ket_on_bragrid_buffer',ierr)
    allocate(fftbox1_dbl(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box','fftbox1_dbl',ierr)
    allocate(fftbox2_dbl(ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box','fftbox2_dbl',ierr)
    allocate(r_op(3,ld1_dbl,ld2_dbl,2*n3), stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box','r_op',ierr)
    allocate(new_proj(ld1, ld2, n3), stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box', &
         'new_proj',ierr)
    allocate(old_proj(ld1, ld2, n3), stat=ierr)
    call utils_alloc_check('projectors_func_pos_ovlp_box', &
         'old_proj',ierr)

    call projectors_init_fftbox_recip(proj_set)

    ! calculate vectors between fine grid points
    a1 = (1.0_DP / pub_fftbox%total_pt1_dbl) * pub_fftbox%a1
    a2 = (1.0_DP / pub_fftbox%total_pt2_dbl) * pub_fftbox%a2
    a3 = (1.0_DP / pub_fftbox%total_pt3_dbl) * pub_fftbox%a3

    ! ddor: construct position operator to power `order' in FFTbox
    r_op = 0.0_DP
    do i3=1,2*n3
       r3 = real(i3-1,kind=DP) * a3
       do i2=1,ld2_dbl
          r2 = r3 + real(i2-1,kind=DP) * a2
          do i1=1,ld1_dbl
             r1 = r2 + real(i1-1,kind=DP) * a1
             r_op(1,i1,i2,i3) = (r1%X)!**order
             r_op(2,i1,i2,i3) = (r1%Y)!**order
             r_op(3,i1,i2,i3) = (r1%Z)!**order
          enddo
       enddo
    enddo

    ! Number of row-steps per row-block
    n_batches = proj_basis%max_on_node / batch_size
    if (mod(proj_basis%max_on_node,batch_size) > 0) n_batches = n_batches + 1

    ! loop over x,y,z 
    do xyz=1,3

       local_start = 1
       do batch_count=1,n_batches
          local_end = min(local_start+batch_size-1,proj_basis%node_num)

          count =1
          do local_proj=local_start,local_end
   
             batch_proj = local_proj - local_start + 1
   
             ! Find information about this projector
             global_proj = local_proj + proj_basis%first_on_node(pub_my_node_id) &
                  - 1
             atom_of_proj = proj_basis%atom_of_func(global_proj)
             loc_atom_of_proj = atom_of_proj - &
                  pub_first_atom_on_node(pub_my_node_id) + 1
             species_number = proj_set%proj_species(atom_of_proj)
             proj_count_atom = global_proj - &
                  proj_basis%first_on_atom(atom_of_proj) + 1
                proj_count = proj_set%species_first_proj(species_number) + &
                  proj_count_atom - 1
   
             ! Centre of projector wrt fftbox in terms of grid spacings
             call basis_point_wrt_box(pcbg1, pcbg2, pcbg3, &  ! output
                  proj_set%proj_centre(atom_of_proj), n1, n2, n3)
   
             ! Start of fftbox wrt cell in terms of grid-point number
             call basis_start_of_box_wrt_cell(fftbox_start(1,batch_proj), &
                  fftbox_start(2,batch_proj),fftbox_start(3,batch_proj), &
                  proj_set%proj_centre(atom_of_proj),pcbg1,pcbg2,pcbg3)

             ! Apply phase factor e^-(ik.proj_centre_wrt_box) to projector
             call basis_phase_on_fftbox_recip(proj_complex, &
                  proj_set%fftbox_proj_recip(:,:,:,proj_count),n1,n2,n3,&
                  ld1,ld2,-pcbg1,-pcbg2,-pcbg3)
   
             ! g=0 element must be real
             proj_complex(1,1,1) &
                  = cmplx(real(proj_complex(1,1,1),kind=DP),0.0_DP,kind=DP)
   
             ! Fourier transform to real space:
             call fourier_apply_box('Coarse','Backward',proj_complex)

             ! really not so sure about this bit
             ! only want real bit
             if (count == 2 .or. local_proj == local_end) then ! for odd cases
                new_proj = real(proj_complex,dp)

                ! interpolate fftboxes to fine grid
                call fourier_interpolate(old_proj,new_proj,&
                     fftbox1_dbl,fftbox2_dbl)

                ! apply FFTbox position operator and increment batch_count
                fftbox1_dbl = r_op(xyz,:,:,:) * fftbox1_dbl
                fftbox2_dbl = r_op(xyz,:,:,:) * fftbox2_dbl

                ! filter fftboxes to standard grid
                call fourier_filter(fftbox1_dbl,fftbox2_dbl,old_proj,&
                     new_proj)

             else
                old_proj = real(proj_complex,dp)
                old_batch_proj = batch_proj
             endif

             ! Put real part into projector_box              
             if (count == 2) then

                projector_fftbox(1:n1,1:n2,1:n3,old_batch_proj) &
                     = old_proj(1:n1,1:n2,1:n3)/pub_cell%weight

                projector_fftbox(1:n1,1:n2,1:n3,batch_proj) &
                     = new_proj(1:n1,1:n2,1:n3)/pub_cell%weight

               count = 1
             else if (local_proj == local_end) then

                projector_fftbox(1:n1,1:n2,1:n3,batch_proj) &
                     = new_proj(1:n1,1:n2,1:n3)/pub_cell%weight
             else
                count = count + 1   
             end if

          end do

      
          ! Calculate overlap integrals
          call integrals_brappd_ketfftbox(sp_overlap(xyz),  &              ! inout
               ngwfs_on_grid, ngwf_basis, &                                ! input
               projector_fftbox, fftbox_start, batch_size,  &              ! input
               local_start, local_end, idx_len, sp_overlap_idx, 'FULL', &  ! input
               bra_on_grid_buffer,ket_on_bragrid_buffer)                   ! work

          local_start = local_start + batch_size

       end do

    end do !xyz

    call projectors_exit_fftbox_recip(proj_set)

    deallocate(old_proj,stat=ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box', &
         'old_proj',ierr)
    deallocate(new_proj,stat=ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box', &
         'new_proj',ierr)
    deallocate(r_op,stat =ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box','r_op',ierr)
    deallocate(fftbox2_dbl,stat =ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box','fftbox2_dbl',ierr)
    deallocate(fftbox1_dbl,stat =ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box','fftbox1_dbl',ierr)

    ! Deallocate workspace
    deallocate(ket_on_bragrid_buffer,stat=ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box', &
         'ket_on_bragrid_buffer',ierr)
    deallocate(bra_on_grid_buffer,stat=ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box', &
         'bra_on_grid_buffer',ierr)
    deallocate(fftbox_start,stat=ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box', &
         'fftbox_start',ierr)
    deallocate(projector_fftbox,stat=ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box', &
         'projector_fftbox',ierr)
    deallocate(proj_complex,stat=ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box', &
         'proj_complex',ierr)
    deallocate(sp_overlap_idx,stat=ierr)
    call utils_dealloc_check('projectors_func_pos_ovlp_box', &
         'sp_overlap_idx',ierr)

    ! Stop timer
    call timer_clock('projectors_func_pos_ovlp_box',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &projectors_func_pos_ovlp_box'
#endif

  end subroutine projectors_func_pos_ovlp_box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_gradient_batch(fftbox_batch, &
       mmat,nmat,local_start,local_end,batch_size, &
       fa_start1,fa_start2,fa_start3,ngwf_basis,proj_basis,coeff_mat, &
       ps_overlap,proj_set)

    !====================================================================!
    ! This subroutine calculates the non-local potential contribution    !
    ! to the gradient of the energy with respect to the NGWF expansion   !
    ! coefficients in the psinc basis.                                   !
    !--------------------------------------------------------------------!
    ! Original version written by Arash A. Mostofi in January 2004.      !
    ! Modified for speed by Chris-Kriton Skylaris on 31/1/2004.          !
    ! Bug related to initialisation of overlap_test_fa located and       !
    ! fixed by Chris-Kriton Skylaris on 2/5/2004.                        !
    ! Modified to use parallel SPAM by Peter Haynes, July 2006.          !
    ! Modified to use batch system by Nicholas Hine, July 2008.          !
    ! Modified to use SPAM3 matrices by Nicholas Hine, July 2009.        !
    !--------------------------------------------------------------------!
    ! Re-written by Nicholas Hine in October 2009. Outer loop is now     !
    ! over NGWFs fa, inner loop split into two parts: first find overlaps!
    ! for each fa and and store the coefficients of RK and RSK for each  !
    ! NGWF-projector overlap, plus phase arrays e^iG.R for each overlap. !
    ! Next, loop over the overlaps found, and calculate a structure      !
    ! factor for each projector type, multiply this by the projector in  !
    ! reciprocal space and FFT back to real space to get the gradient.   !
    !====================================================================!

    use comms, only: pub_on_root, pub_my_node_id
    use constants, only: two_pi, cmplx_i
    use fourier, only: fourier_apply_box
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, geometry_magnitude, OPERATOR(.dot.), &
         OPERATOR(+), OPERATOR(*)
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_first_atom_on_node
    use simulation_cell, only : pub_fftbox, pub_cell
    use sparse, only: SPAM3, sparse_index_length, sparse_generate_index, &
         sparse_get_element
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer, intent(in) :: local_start,local_end,batch_size
    integer, intent(in) :: mmat,nmat
    real(kind=DP), intent(inout) :: fftbox_batch(pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3, pub_cell%num_spins, &
         nmat, batch_size)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    integer, intent(in) :: fa_start1,fa_start2,fa_start3
    type(SPAM3), intent(in) :: coeff_mat(pub_cell%num_spins,nmat)
    type(SPAM3), intent(in) :: ps_overlap
    type(PROJECTOR_SET), intent(inout) :: proj_set

    ! Local Variables
    integer :: n1,n2,n3,ld1,ld2
    integer :: fa, local_fa
    integer :: loc_fa_atom, loc_fa_atom_prev
    integer :: ovlp_atom
    integer :: species
    integer :: proj_count
    integer :: global_proj
    integer :: ierr
    integer :: idx_len
    integer :: fa_atom_overlaps
    integer :: ovlp_list_len
    integer :: iproj, iovlp, novlp
    integer :: i1,i2,i3
    integer :: is,idx,imat
    integer :: batch_count
    complex(kind=DP) :: phase23
    integer, allocatable :: ps_overlap_idx(:)
    integer, allocatable :: ovlp_species(:)
    real(kind=DP), allocatable :: ovlp_proj_coeffs(:,:,:)
    complex(kind=DP), allocatable :: out_complex(:,:,:)
    complex(kind=DP), allocatable :: struct_fac(:,:,:)
    complex(kind=DP), allocatable :: phases1(:,:),phases2(:,:),phases3(:,:)

    ! cks: time it
    call timer_clock('projectors_gradient_batch', 1)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering projectors_gradient_batch'
#endif

    ! Initialisations
    n1 = pub_fftbox%total_pt1
    n2 = pub_fftbox%total_pt2
    n3 = pub_fftbox%total_pt3
    ld1 = pub_fftbox%total_ld1
    ld2 = pub_fftbox%total_ld2
    loc_fa_atom_prev = -9999

    ! Memory allocation
    allocate(out_complex(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','out_complex',ierr)
    allocate(struct_fac(ld1, ld2, n3),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','struct_fac',ierr)

    ! Obtain index of ps_overlap
    idx_len = sparse_index_length(ps_overlap)
    allocate(ps_overlap_idx(idx_len),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','ps_overlap_idx',ierr)
    call sparse_generate_index(ps_overlap_idx,ps_overlap)

    ! Find maximum number of projectors with which any fa in batch overlaps
    ovlp_list_len = 0
    do local_fa=local_start,local_end
       fa = local_fa + ngwf_basis%first_on_node(pub_my_node_id) - 1
       loc_fa_atom = ngwf_basis%atom_of_func(fa) - &
            pub_first_atom_on_node(pub_my_node_id) + 1
       fa_atom_overlaps = ps_overlap_idx(loc_fa_atom+1) - &
            ps_overlap_idx(loc_fa_atom) + 1
       ovlp_list_len = max(ovlp_list_len,fa_atom_overlaps)
    end do

    ! Allocate lists to store overlaps for a given fa
    allocate(ovlp_species(ovlp_list_len),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','ovlp_species',ierr)
    allocate(ovlp_proj_coeffs(ovlp_list_len,proj_basis%max_on_atom,mmat:nmat), &
         stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','ovlp_proj_coeffs',ierr)
    allocate(phases1(1:n1,ovlp_list_len),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','phases1',ierr)
    allocate(phases2(1:n2,ovlp_list_len),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','phases2',ierr)
    allocate(phases3(1:n3,ovlp_list_len),stat=ierr)
    call utils_alloc_check('projectors_gradient_batch','phases3',ierr)

    call projectors_init_fftbox_recip(proj_set)

    ! Loop over spins
    do is=1,pub_cell%num_spins

       ! Loop over NGWFs in batch
       batch_count = 0
       do local_fa=local_start,local_end
          batch_count = batch_count + 1
          fa = local_fa + ngwf_basis%first_on_node(pub_my_node_id) - 1
          loc_fa_atom = ngwf_basis%atom_of_func(fa) - &
               pub_first_atom_on_node(pub_my_node_id) + 1

          ! Reset number of overlaps and projectors for this fa
          novlp = 0

          ! Loop over atoms with projectors overlapping fa
          do idx=ps_overlap_idx(loc_fa_atom),ps_overlap_idx(loc_fa_atom+1)-1
             ovlp_atom = ps_overlap_idx(idx)
             novlp = novlp + 1

             ! If fa is on a new atom, we need to update the entries in the
             ! phases arrays for this overlap
             if (loc_fa_atom/=loc_fa_atom_prev) call internal_phases

             ! Store information about the species of this atom
             ovlp_species(novlp) = proj_set%proj_species(ovlp_atom)

             ! Loop over all the projectors on this atom storing coefficients
             do iproj=1,proj_basis%num_on_atom(ovlp_atom)
                global_proj = proj_basis%first_on_atom(ovlp_atom)+iproj-1

                ! Store \sum_bj 4*Kab*D_ij*<pj|fb> in ovlp_proj_coeffs(:,i,3)
                ! and \sum_bj 4*Kac*Scb*D_ij*<pj|fb> in ovlp_proj_coeffs(:,i,4)
                do imat=mmat,nmat
                   call sparse_get_element(ovlp_proj_coeffs(novlp,iproj,imat), &
                        coeff_mat(is,imat),global_proj,fa)
                end do

             end do

          end do

          ! Loop over the types of gradient fftbox we are accumulating
          do imat=mmat,nmat ! 3 is ham_fftbox_batch, 4 is tc_ham_fftbox_batch

             ! Zero reciprocal space nonlocal gradient box
             out_complex(:,:,:) = cmplx(0.0_DP,0.0_DP)

             ! Loop over atom species
             do species=1,proj_set%n_proj_species

                ! Cycle if fa does not overlap any atoms of this species
                if (.not.any(ovlp_species(1:novlp)==species)) cycle

                ! Loop over all the projector types for this species
                do iproj=1,proj_set%species_num_proj(species)

                   ! Find entry in global projector fftboxes array
                   proj_count = proj_set%species_first_proj(species) + &
                        iproj - 1

                   ! Reset structure factor to zero for new projector type
                   struct_fac(:,:,:) = cmplx(0.0_DP,0.0_DP)

                   ! Loop over slabs of reciprocal space points of fftbox
                   do i3=1,n3

                      ! Loop over atoms with projectors overlapping fa
                      do iovlp=1,novlp

                         ! Cycle if overlapping atom is not of the right species
                         if (.not.ovlp_species(iovlp)==species) cycle

                         ! Loop over recip points in this slab of fftbox to
                         ! calculate structure factor for this projector type
                         do i2=1,n2
                            phase23 = phases2(i2,iovlp) * phases3(i3,iovlp)
                            struct_fac(1:n1,i2,i3) = struct_fac(1:n1,i2,i3) + &
                                 phases1(1:n1,iovlp) * phase23 * &
                                 ovlp_proj_coeffs(iovlp,iproj,imat)
                         end do

                      end do

                      ! Multiply the structure factor for this projector type by
                      ! the projector in reciprocal space, and add it to the
                      ! current reciprocal space gradient
                      out_complex(:,:,i3) = out_complex(:,:,i3) + &
                           struct_fac(:,:,i3) * &
                           proj_set%fftbox_proj_recip(:,:,i3,proj_count)

                   end do

                end do  ! iproj

             end do  ! species_number

             ! Fourier transform gradient contribution to real space:
             call fourier_apply_box('Coarse', 'Backward', out_complex(:,:,:))

             ! Put real part of this into relevant gradient fftbox
             fftbox_batch(1:n1,1:n2,1:n3,is,imat,batch_count) = &
                  fftbox_batch(1:n1,1:n2,1:n3,is,imat,batch_count) + &
                  real(out_complex(1:n1,1:n2,1:n3),kind=DP)

          end do

          ! Store previous local index of fa atom, so that we can tell if
          ! phases arrays need regenerating
          loc_fa_atom_prev = loc_fa_atom

       end do ! loop over NGWFs in batch

    end do ! loop over spins

    call projectors_exit_fftbox_recip(proj_set)

    ! Deallocate Memory
    deallocate(phases3,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','phases3', &
         ierr)
    deallocate(phases2,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','phases2', &
         ierr)
    deallocate(phases1,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','phases1', &
         ierr)
    deallocate(ovlp_proj_coeffs,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','ovlp_proj_coeffs', &
         ierr)
    deallocate(ovlp_species,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','ovlp_species', &
         ierr)
    deallocate(ps_overlap_idx,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','ps_overlap_idx',&
         ierr)
    deallocate(struct_fac,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','struct_fac',&
         ierr)
    deallocate(out_complex,stat=ierr)
    call utils_dealloc_check('projectors_gradient_batch','out_complex',&
         ierr)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Leaving projectors_gradient_batch'
#endif

    ! cks: time it
    call timer_clock('projectors_gradient_batch', 2)

    return

contains

    ! This subroutine calculates the phase shifts for each reciprocal-space
    ! grid point in the FFTbox centred on NGWF fa, given a projector at
    ! position proj_centre.
    ! Written by Nicholas Hine in October 2009 based on bits of
    ! basis_phase_on_fftbox_recip.
    subroutine internal_phases

      implicit none

      real(kind=DP),parameter :: inv_two_pi = 0.5_DP / PI
      real(kind=DP) :: a1,a2,a3
      real(kind=DP) :: box_start1,box_start2,box_start3
      real(kind=DP) :: proj_radius
      complex(kind=DP) :: phase_inc, phase_neg
      type(POINT) :: proj_centre

      proj_centre = proj_set%proj_centre(ovlp_atom)
      proj_radius = proj_set%proj_max_radius(ovlp_atom)

      ! Start of fa fftbox wrt sim cell in terms of numbers of grid points
      call internal_origin_nl_fa(box_start1,box_start2, &
           box_start3,fa_start1,fa_start2,fa_start3, &
           ngwf_basis%tight_boxes(local_fa),proj_radius, &
           proj_centre,ngwf_basis%spheres(local_fa))

      ! Calculate (box_start - proj_centre) in terms of number of grid
      ! points in each lattice direction.
      a1 = (proj_centre .DOT. (pub_cell%b1))
      box_start1 = box_start1 - a1 * inv_two_pi * &
           real(pub_cell%total_pt1,kind=DP)

      a2 = (proj_centre .DOT. (pub_cell%b2))
      box_start2 = box_start2 - a2 * inv_two_pi * &
           real(pub_cell%total_pt2,kind=DP)

      a3 = (proj_centre .DOT. (pub_cell%b3))
      box_start3 = box_start3 - a3 * inv_two_pi * &
           real(pub_cell%total_pt3,kind=DP)

      ! Phases e^(iG1.R) along '1' direction
      phase_inc = exp(two_pi*cmplx_i/real(n1,kind=DP) * box_start1)
      phase_neg = exp(-two_pi*cmplx_i * box_start1)
      phases1(1,novlp) = cmplx(1.0_DP,0.0_DP,kind=DP)
      do i1=2,n1/2+2
         phases1(i1,novlp) = phases1(i1-1,novlp)*phase_inc
      end do
      i1=n1/2+2; phases1(i1,novlp) = phases1(i1,novlp)*phase_neg
      do i1=n1/2+3,n1
         phases1(i1,novlp) = phases1(i1-1,novlp)*phase_inc
      end do

      ! Phases e^(iG2.R) along '2' direction
      phase_inc = exp(two_pi*cmplx_i/real(n2,kind=DP) * box_start2)
      phase_neg = exp(-two_pi*cmplx_i * box_start2)
      phases2(1,novlp) = cmplx(1.0_DP,0.0_DP,kind=DP)
      do i2=2,n2/2+2
         phases2(i2,novlp) = phases2(i2-1,novlp)*phase_inc
      end do
      i2=n2/2+2; phases2(i2,novlp) = phases2(i2,novlp)*phase_neg
      do i2=n2/2+3,n2
         phases2(i2,novlp) = phases2(i2-1,novlp)*phase_inc
      end do

      ! Phases e^(iG3.R) along '3' direction
      phase_inc = exp(two_pi*cmplx_i/real(n3,kind=DP) * box_start3)
      phase_neg = exp(-two_pi*cmplx_i * box_start3)
      phases3(1,novlp) = cmplx(1.0_DP,0.0_DP,kind=DP)
      do i3=2,n3/2+2
         phases3(i3,novlp) = phases3(i3-1,novlp)*phase_inc
      end do
      i3=n3/2+2; phases3(i3,novlp) = phases3(i3,novlp)*phase_neg
      do i3=n3/2+3,n3
         phases3(i3,novlp) = phases3(i3-1,novlp)*phase_inc
      end do

    end subroutine internal_phases

    subroutine internal_origin_nl_fa( &
         box_start1, box_start2, box_start3, &               ! output
         fa_start1, fa_start2, fa_start3, fa_fft_box, &      ! input
         proj_radius, proj_centre, ngwf_sphere)              ! input

      !============================================================!
      ! Returns the distance, in terms of gridpoints (real values) !
      ! in each lattice direction, from the simulation cell origin !
      ! to the nearest vertex of the fftbox.                       !
      !------------------------------------------------------------!
      ! Written by Arash A. Mostofi, August 2002, using some parts !
      ! from older subroutines written by Chris-Kriton Skylaris.   !
      ! Modified by Arash A. Mostofi on July 2003 and on           !
      ! January 2004.                                              !
      !============================================================!

      use basis, only: FUNCTION_TIGHT_BOX, SPHERE
      use geometry, only: POINT, geometry_magnitude, &
           OPERATOR(.dot.), OPERATOR(+), OPERATOR(*)
      use simulation_cell, only : pub_cell

      implicit none

      type(FUNCTION_TIGHT_BOX), intent(in) :: fa_fft_box
      type(SPHERE), intent(in)           :: ngwf_sphere
      type(POINT), intent(in)            :: proj_centre
      integer, intent(in)                :: fa_start1,fa_start2,fa_start3
      real(kind=DP), intent(in)          :: proj_radius
      real(kind=DP), intent(out)         :: box_start1,box_start2,box_start3

      real(kind=DP),parameter :: inv_two_pi = 0.5_DP / PI

      ! aam: internal declarations
      integer       :: s1,s2,s3
      real(kind=DP) :: sum_radii
      real(kind=DP) :: dist1,dist2,dist3
      type(POINT)   :: ngwf_proj

      ! aam: find start of box with respect to origin of
      !      simulation cell in terms of numbers of gridpoints.
      s1 = (fa_fft_box%start_ppds1 - 1)*pub_cell%n_pt1 &
           + fa_fft_box%start_pts1 - fa_start1

      s2 = (fa_fft_box%start_ppds2 - 1)*pub_cell%n_pt2 &
           + fa_fft_box%start_pts2 - fa_start2

      s3 = (fa_fft_box%start_ppds3 - 1)*pub_cell%n_pt3 &
           + fa_fft_box%start_pts3 - fa_start3

      box_start1 = real(s1,kind=DP)
      box_start2 = real(s2,kind=DP)
      box_start3 = real(s3,kind=DP)

      ! cks: sum of projector and nwgf radii
      sum_radii = proj_radius + ngwf_sphere%radius

      ! cks: vector from centre of projector to centre of ngwf
      ngwf_proj = ngwf_sphere%centre+(-1.0_DP)*proj_centre

      ! cks: a1_unit, a2_unit and a3_unit components of ngwf_proj
      dist1 = (ngwf_proj.DOT.pub_cell%b1) * &
           inv_two_pi*geometry_magnitude(pub_cell%a1)
      dist2 = (ngwf_proj.DOT.pub_cell%b2) * &
           inv_two_pi*geometry_magnitude(pub_cell%a2)
      dist3 = (ngwf_proj.DOT.pub_cell%b3) * &
           inv_two_pi*geometry_magnitude(pub_cell%a3)

      ! cks: if necessary, move the origin of the ngwf-fftbox by a lattice
      ! cks: vector in each direction. This need arises when the projector
      ! cks: and the ngwf overlap not inside the simulation cell but through
      ! cks: its faces due to the periodic boundary conditions. Then, the
      ! cks: displacement of the ngwf-fftbox places it in the position of
      ! cks: one of its periodic images so that the origin of the projector
      ! cks: w.r.t the origin of the box will have the correct relationship
      ! cks: with the origin of the ngwf w.r.t. the box.
      if ( abs(dist1) > sum_radii ) box_start1 = box_start1 &
           - SIGN(1.0_DP,dist1)*real(pub_cell%total_pt1,kind=DP)

      if ( abs(dist2) > sum_radii ) box_start2 = box_start2 &
           - SIGN(1.0_DP,dist2)*real(pub_cell%total_pt2,kind=DP)

      if ( abs(dist3) > sum_radii ) box_start3 = box_start3 &
           - SIGN(1.0_DP,dist3)*real(pub_cell%total_pt3,kind=DP)

    end subroutine internal_origin_nl_fa

  end subroutine projectors_gradient_batch


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine projectors_commutator(nonlocpot_com, proj_basis, &
       ngwf_basis, ngwfs_on_grid, sp_overlap, proj_set, delta, dij)

    !==================================================================!
    ! This subroutine calculates the commutator between the nonlocal   !
    ! potential and the position operator for the 3 Cartesian          !
    ! directions.                                                      !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff February 2011.                         !
    ! Modified by Nicholas Hine in April 2011 to not use NGWF_REP and  !
    ! to simplify and cleanup existing code, and for changes to        !
    ! projector initialisation routines.                               !
    ! Moved to projectors_mod by Nicholas Hine in December 2011.       !
    ! Modified by Laura Ratcliff Dec 2011 to fix bug involving complex !
    ! matrices and phase factors.                                      !
    !==================================================================!

    use comms, only: pub_on_root
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product, &
         sparse_transpose, sparse_axpy, sparse_scale, sparse_copy, &
         sparse_transpose_structure
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3), intent(inout) :: nonlocpot_com(3)
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(SPAM3), intent(in) :: sp_overlap
    type(PROJECTOR_SET), intent(inout) :: proj_set
    real(kind=DP), intent(in) :: delta
    type(SPAM3), intent(in) :: dij

    ! Local Variables
    type(SPAM3) :: sp_overlap_dij, ps_overlap, dij_ps_overlap
    type(SPAM3) :: sDGp_overlap, DGps_overlap, sp_overlap_cmplx
    type(SPAM3) :: ps_overlap_rc, sDGp_overlap_rc
    type(SPAM3) :: nonlocpot_com2
    integer :: cart

    ! NB not allowing for sp_overlap or ps_overlap to be at diff k-points
    call sparse_create(nonlocpot_com2,nonlocpot_com(1),iscmplx=.true.)
    call sparse_create(sp_overlap_dij,sp_overlap,iscmplx=.true.)
    call sparse_create(sDGp_overlap,sp_overlap,iscmplx=.true.)
    call sparse_create(sDGp_overlap_rc,sp_overlap,iscmplx=.false.)

    ! ndmh: Create ps_overlap (modifying structure code if required)
    call sparse_transpose_structure(ps_overlap%structure,sp_overlap)
    ps_overlap_rc%structure = ps_overlap%structure
    dij_ps_overlap%structure = ps_overlap%structure
    DGps_overlap%structure = ps_overlap%structure

    call sparse_create(dij_ps_overlap,iscmplx=.true.)
    call sparse_create(DGps_overlap,iscmplx=.true.)
    call sparse_create(sp_overlap_cmplx,sp_overlap,iscmplx=.true.)
    call sparse_create(ps_overlap,iscmplx=.true.)
    call sparse_create(ps_overlap_rc,iscmplx=.false.)

    ! Create the transpose of the NGWF-projector overlap matrix
    call sparse_transpose(ps_overlap_rc,sp_overlap)

    call sparse_copy(ps_overlap,ps_overlap_rc)
    call sparse_copy(sp_overlap_cmplx,sp_overlap)

    call sparse_destroy(ps_overlap_rc)

    call sparse_product(sp_overlap_dij,sp_overlap_cmplx,dij)
    call sparse_product(dij_ps_overlap,dij,ps_overlap)

    call sparse_destroy(ps_overlap)
    call sparse_destroy(sp_overlap_cmplx)

    ! Loop over Cartesian directions
    do cart=1,3

       ! Calculate <phi|Del_G(proj(G))> overlap matrix
       ! real part
       call projectors_func_ovlp_box(sDGp_overlap_rc, &
            ngwfs_on_grid, ngwf_basis, proj_basis, proj_set, delta=delta, &
            cart=cart, swap_rc=.false.)

       ! ndmh: add real part to complex sp_overlap matrix
       call sparse_copy(sDGp_overlap,sDGp_overlap_rc)

       ! imag part
       call projectors_func_ovlp_box(sDGp_overlap_rc, &
            ngwfs_on_grid, ngwf_basis, proj_basis, proj_set, delta=delta, &
            cart=cart, swap_rc=.true.)

       ! ndmh: add imaginary part to complex sp_overlap matrix
       call sparse_axpy(sDGp_overlap,sDGp_overlap_rc,(0.0_DP,1.0_DP))
 
       ! Transpose <phi|Del_G(proj(G))> to get <Del_G(proj(G))|phi> for ngwf_basis
       call sparse_transpose(DGps_overlap,sDGp_overlap)

       ! Calculate the matrix \sum_i (<NGWF_a|Del_G(Proj(G))_i>D_ij<Proj_j|NGWF_b>)
       call sparse_product(nonlocpot_com2,sDGp_overlap,dij_ps_overlap)

       ! Calculate the matrix \sum_i (<NGWF_a|Proj_i>D_ij<Del_G(Proj(G))_i|NGWF_b>)
       call sparse_product(nonlocpot_com(cart),sp_overlap_dij,DGps_overlap)

       ! Sum the two
       call sparse_axpy(nonlocpot_com(cart),nonlocpot_com2,(1.0_DP,0.0_DP))

       ! Multiply by -i
       call sparse_scale(nonlocpot_com(cart),(0.0_DP,-1.0_DP))

    end do

    call sparse_destroy(DGps_overlap)
    call sparse_destroy(dij_ps_overlap)
    call sparse_destroy(sDGp_overlap_rc)
    call sparse_destroy(sDGp_overlap)
    call sparse_destroy(sp_overlap_dij)
    call sparse_destroy(nonlocpot_com2)

  end subroutine projectors_commutator


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module projectors
