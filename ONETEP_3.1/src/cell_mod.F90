! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by Chris-Kriton Skylaris
!
!   August 2000
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module simulation_cell

  use geometry, only: point
  use constants, only: DP

  implicit none

  private

  type CELL_INFO
     ! cks: lattice vectors
     type(POINT)      :: a1, a2, a3
     ! cks: lattice vectors normalised to unit length
     type(POINT)      :: a1_unit, a2_unit, a3_unit
     ! cks: reciprocal lattice vectors
     type(POINT)      :: b1, b2, b3
     ! cks: regular grid spacing in atomic units
     real(kind=DP)   :: d1, d2, d3
     ! cks: weight for integrals on the standard grid
     real(kind=DP)   :: weight
     ! cks: plane-wave kinetic energy cutoff (in Eh)
     real(kind=DP)   :: ke_cutoff
     ! cks: total number of ppds in the simulation cell
     integer          :: n_ppds
     ! cks: total number of ppds in each lattice vector direction
     integer          :: n_ppds_a1,n_ppds_a2,n_ppds_a3
     ! cks: number of points per ppd on the standard grid
     integer          :: n_pts
     ! cks: number of points per ppd in each lattice vector dir on standard grid
     integer          :: n_pt1,n_pt2,n_pt3
     ! cks: total number of points in each lattice direction on the standard grid
     integer          :: total_pt1,total_pt2,total_pt3
     ! cks: the total number of atoms in the cell
     integer          :: nat
     ! ndmh: the total number of Hubbard atoms in the cell
     integer          :: nat_hub
     ! ars: the total number of classical atoms in the cell
     integer          :: nat_classical
     ! cks: the total number of NGWFs in the cell
     integer          :: num_ngwfs
     ! ndmh: the total number of conduction NGWFs in the cell
     integer          :: num_ngwfs_cond
     ! ndmh: the total number of auxiliary NGWFs in the cell
     integer          :: num_ngwfs_aux
     ! cks: the total number of (non local pseudopotential) projectors in the cell
     integer          :: num_projectors
     ! ndmh: the total number of PAW partial waves in the cell
     integer          :: num_pawpws
     ! ndmh: the number of hubbard projectors in the cell
     integer          :: num_hub_proj
     ! lr408: the total number of core wavefunctions in the cell ?!
     integer          :: num_corewfs
     ! qoh: the total number of distinct pseudopotential species in the cell
     integer          :: num_pspecies
     ! qoh: the total number of distinct atomic species in the cell
     integer          :: num_species
     ! ndmh: the total number of hubbard species in the cell
     integer          :: num_hub_species
     ! ndmh: the total number of hubbard species in the cell
     integer          :: num_cdft_species
     ! qoh: the total number of SW basis functions in the cell
     integer          :: num_sw
     ! pdh: number of spins
     integer          :: num_spins
     ! ndmh: spin degeneracy factor
     real(kind=DP)    :: spin_fac
  end type CELL_INFO


  type FFTBOX_INFO
     ! lattice vectors
     type(POINT)      :: a1,a2,a3
     ! lattice vectors normalised to unit length
     type(POINT)      :: a1_unit,a2_unit,a3_unit
     ! reciprocal lattice vectors
     type(POINT)      :: b1,b2,b3
     ! regular grid spacing in atomic units
     real(kind=DP)    :: d1,d2,d3
     ! weight for integrals on the standard grid
     real(kind=DP)    :: weight
     ! total number of points in each lattice direction on the standard grid
     integer          :: total_pt1,total_pt2,total_pt3
     ! the total number of points in each lattice direction on the double grid
     integer          :: total_pt1_dbl,total_pt2_dbl,total_pt3_dbl
     ! pdh: additions to improve FFT efficiency
     integer          :: total_ld1,total_ld2
     integer          :: total_ld1_dbl,total_ld2_dbl
     ! reciprocal lattice vector grid
     real(kind=DP), pointer, dimension(:,:,:,:) :: recip_grid
     !
     ! Note: Allocatable arrays are not valid structure components in Fortran
     !       90/95 but this feature was introduced in Fortran 2003.  Pointers
     !       are therefore required instead.
     !
     ! reciprocal lattice vector fine grid
     real(kind=DP), pointer, dimension(:,:,:,:) :: recip_grid_dbl
     ! cks: true if FFT box length coincides with sim cell along a1
     logical          :: coin1
     ! cks: true if FFT box length coincides with sim cell along a2
     logical          :: coin2
     ! cks: true if FFT box length coincides with sim cell along a3
     logical          :: coin3

  end type FFTBOX_INFO

  type castep_cell                                         ! CASTEP cell definition
     real(kind=dp), dimension(3,3) :: real_lattice
     real(kind=dp), dimension(3,3) :: recip_lattice
     real(kind=dp) :: volume
     real(kind=dp), dimension(:,:,:), pointer :: ionic_positions
     integer :: num_ions
     integer :: num_species
     integer :: max_ions_in_species
     integer, dimension(:), pointer :: num_ions_in_species
     real(kind=dp), dimension(:), pointer :: ionic_charge  ! pa: changed from integer to real
     real(kind=dp), dimension(:), pointer :: species_mass
     character(len=2), dimension(:), pointer:: species_symbol
  end type castep_cell


  type castep_model                                        ! CASTEP model definition
     type(castep_cell) :: cell
     type(castep_cell) :: orig_cell
     type(castep_cell) :: ref_cell
     type(castep_cell) :: reac_cell ! Reactant geometry for transition state runs
     type(castep_cell) :: prod_cell ! Product geometry for transition state runs
     type(castep_cell) :: intm_cell ! Intermediate geometry for transition state runs
     real(kind=dp) :: total_energy
     real(kind=dp), dimension(:,:,:), pointer :: forces
     integer :: bfgs_iteration
     integer :: lbfgs_iteration
     real(kind=dp), dimension(:,:), pointer :: bfgs_inv_hessian
     real(kind=dp), dimension(:,:), pointer :: lbfgs_position_updates
     real(kind=dp), dimension(:,:), pointer :: lbfgs_gradient_updates
     integer :: lbfgs_num_updates
     logical :: bfgs_optimisation
     logical :: lbfgs_optimisation
     logical :: found_ground_state
     real(kind=dp), dimension(6) :: stress   ! stress tensor
     real(kind=dp), dimension(3,3) :: strain ! strain tensor
     ! aam: in presence of constraints we keep a copy of the uncontrained forces
     real(kind=dp), dimension(:,:,:), pointer :: orig_forces
  end type castep_model



  public :: CELL_INFO
  public :: FFTBOX_INFO
  public :: castep_cell
  public :: castep_model
  public :: simulation_cell_initialise
  public :: simulation_cell_add_padding
  public :: simulation_cell_fftbox_init
  public :: simulation_cell_fftbox_exit
  public :: castep_cell_copy
  public :: castep_cell_alloc
  public :: castep_cell_dealloc
  public :: castep_cell_nullify
  public :: castep_cell_cart_to_frac
  public :: castep_cell_frac_to_cart
  public :: castep_model_alloc
  public :: castep_model_dealloc
  public :: castep_model_nullify
  public :: castep_model_cell_changed
  public :: castep_cell_cart_lattice_to_abc
  public :: castep_cell_2_elements
  public :: copy_to_castep_cell
  public :: minimum_image_distance

  ! cks: information for the fft_box stored in a CELL_INFO structure
  ! cks: Now pub_fft_box is a public variable available to any subroutine
  ! cks: that uses the simulation_cell module
  type(FFTBOX_INFO), public :: pub_fftbox

  ! cks: same for public (simulation) cell structure
  type(CELL_INFO), public :: pub_cell

  ! ndmh: In cutoff coulomb calculations, pub_padded_cell is a CELL_INFO
  ! ndmh: structure for the padded cell
  type(CELL_INFO), public :: pub_padded_cell

  ! ndmh: Definitions for universal tightbox, for FFTs in tbs
  integer, public :: pub_maxtight_pts1
  integer, public :: pub_maxtight_pts2
  integer, public :: pub_maxtight_pts3

  ! ndmh: sizes of augmentation function box in each direction
  integer, public :: pub_aug_box_n1
  integer, public :: pub_aug_box_n2
  integer, public :: pub_aug_box_n3

  ! ndmh: Reciprocal grid vectors and magnitudes for unversal tightboxes
  real(kind=DP), allocatable, public :: pub_tb_recip_grid(:,:,:,:)
  real(kind=DP), allocatable, public :: pub_tb_recip_grid_fine(:,:,:,:)

  ! ndmh: Cutoff Coulomb extra padding required
  real(kind=DP), parameter, public :: cutoff_coulomb_tol = 10.0_DP

  !--------------------------------------------------------------------------!
  ! Overload subroutines
  !--------------------------------------------------------------------------!

  interface castep_cell_frac_to_cart
     module procedure castep_cell_frac_to_cart_cell
     module procedure castep_cell_frac_to_cart_vector
  end interface


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simulation_cell_fftbox_init(n1 ,n2, n3)

    !============================================================!
    ! This subroutine initialises the components of an           !
    ! FFTBOX_INFO  type variable which groups together           !
    ! information about the FFTbox.                              !
    !------------------------------------------------------------!
    ! WARNING: This subroutine should be called only AFTER       !
    !          the pub_cell variable has been initialised!       !
    !------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 13/5/2003.             !
    ! Reciprocal grid array added by Peter Haynes on 15/7/2004   !
    ! Reciprocal grid fine array added by Quintin Hill on        !
    ! 24/04/2008.                                                !
    !============================================================!

    use comms, only: comms_abort, pub_on_root
    use constants, only: stdout
    use geometry, only: operator(*)
    use rundat, only: pub_usehfx, pub_dbl_grid_scale
    use utils, only : utils_alloc_check
    implicit none

    ! Arguments
    integer, intent(in) :: n1, n2, n3

    ! Local variables
    integer :: ierr
    integer :: i1,i2,i3
    integer :: k1,k2,k3
    integer :: n1half,n2half,n3half
    real(kind=DP) :: b1(3),b2(3),b3(3)
    real(kind=DP) :: g3(3),g23(3),g(3)
    real(kind=DP) :: gsq
    real(kind=DP) :: coulomb_cutoff

    if (mod(n1,2) == 0 .or. mod(n2,2) == 0 .or. mod(n3,2) == 0) then
       if (pub_on_root) then
          write(stdout,'(a)') 'Error in simulation_cell_fftbox_init: &
               &at least one dimension of the FFT box is even'
          write(stdout,'(3(a,i8))') ', n1 = ',n1,', n2 = ',n2,', n3 = ',n3
       end if
       call comms_abort
    end if

    pub_fftbox%total_pt1 = n1
    pub_fftbox%total_pt2 = n2
    pub_fftbox%total_pt3 = n3

    ! cks: set flags if FFT box coincides with sim cell along any lattice vector
    pub_fftbox%coin1 =.false.
    if (pub_fftbox%total_pt1 == pub_cell%total_pt1) pub_fftbox%coin1 =.true.
    pub_fftbox%coin2 =.false.
    if (pub_fftbox%total_pt2 == pub_cell%total_pt2) pub_fftbox%coin2 =.true.
    pub_fftbox%coin3 =.false.
    if (pub_fftbox%total_pt3 == pub_cell%total_pt3) pub_fftbox%coin3 =.true.

    if (pub_dbl_grid_scale > 1.0_DP) then
       pub_fftbox%total_pt1_dbl = 2*n1
       pub_fftbox%total_pt2_dbl = 2*n2
       pub_fftbox%total_pt3_dbl = 2*n3
    else
       pub_fftbox%total_pt1_dbl = n1
       pub_fftbox%total_pt2_dbl = n2
       pub_fftbox%total_pt3_dbl = n3
    end if

    pub_fftbox%a1 =( real(n1, kind=DP)/real(pub_cell%total_pt1, kind=DP) ) &
         * pub_cell%a1
    pub_fftbox%a2 =( real(n2, kind=DP)/real(pub_cell%total_pt2, kind=DP) ) &
         * pub_cell%a2
    pub_fftbox%a3 =( real(n3, kind=DP)/real(pub_cell%total_pt3, kind=DP) ) &
         * pub_cell%a3

    pub_fftbox%a1_unit =pub_cell%a1_unit
    pub_fftbox%a2_unit =pub_cell%a2_unit
    pub_fftbox%a3_unit =pub_cell%a3_unit

    pub_fftbox%b1 =( real(pub_cell%total_pt1, kind=DP)/real(n1, kind=DP) ) &
         * pub_cell%b1
    pub_fftbox%b2 =( real(pub_cell%total_pt2, kind=DP)/real(n2, kind=DP) ) &
         * pub_cell%b2
    pub_fftbox%b3 =( real(pub_cell%total_pt3, kind=DP)/real(n3, kind=DP) ) &
         * pub_cell%b3

    pub_fftbox%d1 =pub_cell%d1
    pub_fftbox%d2 =pub_cell%d2
    pub_fftbox%d3 =pub_cell%d3

    pub_fftbox%weight =pub_cell%weight

    ! pdh: improve efficiency of FFTs - assumed complex-to-complex here
    ! aam: only ever do C2C FFTs on the FFT box, therefore don't require
    !      the extra two rows in the first dimension of either the coarse
    !      or fine FFT boxes.
    ! pdh: And FFTw will not work if ld1 /= n1 and ld2 /= n2!
    !    pub_fftbox%total_ld1=pub_fftbox%total_pt1+2
    pub_fftbox%total_ld1=pub_fftbox%total_pt1
    pub_fftbox%total_ld2=pub_fftbox%total_pt2
    !    pub_fftbox%total_ld1_dbl=pub_fftbox%total_pt1_dbl+2
    pub_fftbox%total_ld1_dbl=pub_fftbox%total_pt1_dbl
    pub_fftbox%total_ld2_dbl=pub_fftbox%total_pt2_dbl

    ! pdh: add reciprocal grid
    allocate(pub_fftbox%recip_grid(6,pub_fftbox%total_ld1, &
         pub_fftbox%total_ld2,pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('simulation_cell_fftbox_init','pub_fftbox%recip_grid', ierr)

    ! local copies of FFT box reciprocal lattice vectors
    b1(1) = pub_fftbox%b1%x ; b1(2) = pub_fftbox%b1%y ; b1(3) = pub_fftbox%b1%z
    b2(1) = pub_fftbox%b2%x ; b2(2) = pub_fftbox%b2%y ; b2(3) = pub_fftbox%b2%z
    b3(1) = pub_fftbox%b3%x ; b3(2) = pub_fftbox%b3%y ; b3(3) = pub_fftbox%b3%z

    coulomb_cutoff = 0.48_DP*min(&
         (real(pub_fftbox%total_pt1,kind=DP)*pub_fftbox%d1),&
         (real(pub_fftbox%total_pt2,kind=DP)*pub_fftbox%d2),&
         (real(pub_fftbox%total_pt3,kind=DP)*pub_fftbox%d3))

    ! loop over FFT box reciprocal grid
    pub_fftbox%recip_grid = 0.0_DP
    n1half = pub_fftbox%total_pt1/2+1
    n2half = pub_fftbox%total_pt2/2+1
    n3half = pub_fftbox%total_pt3/2+1
    do i3=1,pub_fftbox%total_pt3
       if (i3 > n3half) then
          k3 = i3 - pub_fftbox%total_pt3 - 1
       else
          k3 = i3 - 1
       end if
       g3 = k3 * b3
       do i2=1,pub_fftbox%total_pt2
          if (i2 > n2half) then
             k2 = i2 - pub_fftbox%total_pt2 - 1
          else
             k2 = i2 - 1
          end if
          g23 = g3 + k2 * b2
          do i1=1,pub_fftbox%total_pt1
             if (i1 > n1half) then
                k1 = i1 - pub_fftbox%total_pt1 - 1
             else
                k1 = i1 - 1
             end if
             g = g23 + k1 * b1
             pub_fftbox%recip_grid(1:3,i1,i2,i3) = g
             gsq = g(1)*g(1) + g(2)*g(2) + g(3)*g(3)
             pub_fftbox%recip_grid(4,i1,i2,i3) = sqrt(gsq)
             pub_fftbox%recip_grid(5,i1,i2,i3) = 0.5_DP * gsq
             ! qoh: Treat case of G=0 specially
             if (i1 == 1 .and. i2 == 1 .and. i3 == 1) then
                pub_fftbox%recip_grid(6,1,1,1)=0.5_DP*coulomb_cutoff**2
             else
                pub_fftbox%recip_grid(6,i1,i2,i3) &
                     = (1.0_DP - cos(sqrt(gsq)*coulomb_cutoff)) / gsq
             end if
          end do
       end do
    end do

    !qoh: Initialise pub_fftbox%recip_grid_dbl
    !qoh: pub_fftbox%recip_grid_dbl(6,:,:,:) is the cutoff coulomb energy
    !qoh: pub_fftbox%recip_grid_dbl(1:5,:,:,:) are reserved for future use and
    !qoh: will correspond to pub_fftbox%recip_grid(1:5,:,:,:)

    if (pub_usehfx) then

       allocate(pub_fftbox%recip_grid_dbl(6:6,pub_fftbox%total_ld1_dbl,&
            pub_fftbox%total_ld2_dbl,pub_fftbox%total_pt3_dbl),stat=ierr)
       call utils_alloc_check('simulation_cell_fftbox_init',&
            'pub_fftbox%recip_grid_dbl', ierr)
       pub_fftbox%recip_grid_dbl = 0.0_DP

       n1half = pub_fftbox%total_pt1_dbl/2+1
       n2half = pub_fftbox%total_pt2_dbl/2+1
       n3half = pub_fftbox%total_pt3_dbl/2+1
       coulomb_cutoff = 0.48_DP*min(&
            &(real(pub_fftbox%total_pt1,kind=DP)*pub_fftbox%d1),&
            &(real(pub_fftbox%total_pt2,kind=DP)*pub_fftbox%d2),&
            &(real(pub_fftbox%total_pt3,kind=DP)*pub_fftbox%d3))

       do i3=1,pub_fftbox%total_pt3_dbl
          if (i3 > n3half) then
             k3 = i3 - pub_fftbox%total_pt3_dbl - 1
          else
             k3 = i3 - 1
          end if
          g3 = k3 * b3
          do i2=1,pub_fftbox%total_pt2_dbl
             if (i2 > n2half) then
                k2 = i2 - pub_fftbox%total_pt2_dbl - 1
             else
                k2 = i2 - 1
             end if
             g23 = g3 + k2 * b2
             do i1=1,pub_fftbox%total_pt1_dbl
                if (i1 > n1half) then
                   k1 = i1 - pub_fftbox%total_pt1_dbl - 1
                else
                   k1 = i1 - 1
                end if
                g = g23 + k1 * b1
                !qoh: Not needed yet
                !pub_fftbox%recip_grid_dbl(1:3,i1,i2,i3) = g
                gsq = g(1)*g(1) + g(2)*g(2) + g(3)*g(3)
                !qoh: Not needed yet
                !pub_fftbox%recip_grid_dbl(4,i1,i2,i3) = sqrt(gsq)
                !pub_fftbox%recip_grid_dbl(5,i1,i2,i3) = 0.5_DP * gsq
                ! qoh: Treat case of G=0 specially
                if (i1 == 1 .and. i2 == 1 .and. i3 == 1) then
                   pub_fftbox%recip_grid_dbl(6,1,1,1)=0.5_DP*coulomb_cutoff**2
                else
                   pub_fftbox%recip_grid_dbl(6,i1,i2,i3) &
                        = (1.0_DP - cos(sqrt(gsq)*coulomb_cutoff)) / gsq
                end if
             end do
          end do
       end do
    else
       !qoh: nullify() to prevent the dreaded undefined pointer status.
       nullify(pub_fftbox%recip_grid_dbl)
    end if

  end subroutine simulation_cell_fftbox_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simulation_cell_fftbox_exit

    !============================================================!
    ! This subroutine deallocates pointers in fftbox arrays.     !
    !------------------------------------------------------------!
    ! Written by Victor Milman on 03/11/2006.                    !
    ! Modified by Quintin Hill on 24/04/2008.                    !
    !============================================================!
    use utils, only : utils_dealloc_check

    implicit none

    integer :: ierr

    if (associated(pub_fftbox%recip_grid)) then
       deallocate(pub_fftbox%recip_grid,stat=ierr)
       call utils_dealloc_check('simulation_cell_fftbox_exit (cell_mod.F90)', &
            'pub_fftbox%recip_grid',ierr)
    end if

    !qoh: Deallocate pub_fftbox%recip_grid_dbl
    if (associated(pub_fftbox%recip_grid_dbl)) then
       deallocate(pub_fftbox%recip_grid_dbl,stat=ierr)
       call utils_dealloc_check('simulation_cell_fftbox_exit (cell_mod.F90)', &
             'pub_fftbox%recip_grid_dbl',ierr)
    end if

  end subroutine simulation_cell_fftbox_exit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine simulation_cell_add_padding(d1,d2,d3,n_pt1,n_pt2,n_pt3, &
       a1,a2,a3,a1_pad,a2_pad,a3_pad)

    !================================================================!
    ! This subroutine checks that the padded cell is commensurate    !
    ! with the original cell, or sets up the padding if it is left   !
    ! to the default.                                                !
    !----------------------------------------------------------------!
    ! Written by Nicholas Hine in September 2010.                    !
    !================================================================!

    use geometry, only: POINT, operator(+), operator(*), geometry_magnitude
    use rundat, only: pub_coulomb_radius, pub_coulomb_length
    use utils, only: utils_abort

    implicit none

    ! Arguments
    real(kind=DP), intent(in) :: d1,d2,d3
    integer, intent(in) :: n_pt1, n_pt2, n_pt3
    type(POINT), intent(in) :: a1, a2, a3
    type(POINT), intent(inout) :: a1_pad, a2_pad, a3_pad

    ! Local Variables
    type(POINT) :: da1, da2, da3
    real(kind=DP) :: len, len_pad
    real(kind=DP) :: cutoff

    cutoff = sqrt(pub_coulomb_radius**2 + pub_coulomb_length**2)

    ! ndmh: check if the default is still set, if so, override with a suitable
    ! ndmh: padded cell
    if ((a1_pad%x==-1.0_DP).and.(a1_pad%y== 0.0_DP).and.(a1_pad%z== 0.0_DP).and. &
        (a2_pad%x== 0.0_DP).and.(a2_pad%y==-1.0_DP).and.(a2_pad%z== 0.0_DP).and. &
        (a3_pad%x== 0.0_DP).and.(a3_pad%y== 0.0_DP).and.(a3_pad%z==-1.0_DP)) then

       ! Copy the original cell
       a1_pad = a1
       a2_pad = a2
       a3_pad = a3

       ! Increase a1 by one ppd at a time till we are large enough
       len = geometry_magnitude(a1)
       da1 = (real(n_pt1,kind=DP)/(len/d1)) * a1
       do
           len_pad = geometry_magnitude(a1_pad)
           if (len_pad >= len + cutoff + cutoff_coulomb_tol) exit
           a1_pad = a1_pad + da1
       end do

       ! Now increase a2 by one ppd at a time till we are large enough
       len = geometry_magnitude(a2)
       da2 = (real(n_pt2,kind=DP)/(len/d2)) * a2
       do
           len_pad = geometry_magnitude(a2_pad)
           if (len_pad >= len + cutoff + cutoff_coulomb_tol) exit
           a2_pad = a2_pad + da2
       end do

       ! Now increase a3 by one ppd at a time till we are large enough
       len = geometry_magnitude(a3)
       da3 = (real(n_pt3,kind=DP)/(len/d3)) * a3
       do
           len_pad = geometry_magnitude(a3_pad)
           if (len_pad >= len + cutoff + cutoff_coulomb_tol) exit
           a3_pad = a3_pad + da3
       end do

    end if

    ! ndmh: Check for badly set-up padded cells
    if (geometry_magnitude(a1_pad)<geometry_magnitude(a1)) then
        call utils_abort('ERROR in simulation_cell_add_padding: Padded cell &
             &is smaller than original cell along a1')
    end if
    if (geometry_magnitude(a2_pad)<geometry_magnitude(a2)) then
        call utils_abort('ERROR in simulation_cell_add_padding: Padded cell &
             &is smaller than original cell along a2')
    end if
    if (geometry_magnitude(a3_pad)<geometry_magnitude(a3)) then
        call utils_abort('ERROR in simulation_cell_add_padding: Padded cell &
             &is smaller than original cell along a3')
    end if

  end subroutine simulation_cell_add_padding


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine simulation_cell_initialise(cell,elements,num,num_cond,num_aux, &
       nat,class_nat,a1,a2,a3,n_pt1,n_pt2,n_pt3,d1,d2,d3)

    !============================================================!
    ! This subroutine initialises the components of a CELL_INFO  !
    ! type which groups together information about the           !
    ! simulation cell.                                           !
    !------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 21/6/2001.             !
    ! Changed to accept cell as an argument to allow setting up  !
    ! of padded cell for cutoff coulomb by Nick Hine 01/04/2008  !
    !============================================================!

    use constants, only: DP, PI
    use geometry, only: point, unit_vector, operator(.cross.), operator(*),&
        operator(.dot.), geometry_magnitude
    use ion, only: element

    implicit none

    type(CELL_INFO), intent(out) :: cell
    integer, intent(in) :: num, num_cond, num_aux, nat, class_nat
    type(ELEMENT), intent(in) :: elements(nat)
    type(POINT), intent(in) :: a1, a2, a3
    real(kind=DP), intent(in) :: d1, d2, d3
    integer, intent(in) :: n_pt1, n_pt2, n_pt3

    ! cks: internal variables
    integer :: atom, row
    logical :: new_species

    ! cks: initilise grid spacing
    cell%d1 = d1
    cell%d2 = d2
    cell%d3 = d3

    ! cks: initialise primitive lattice vectors
    cell%a1 = a1
    cell%a2 = a2
    cell%a3 = a3

    ! cks: initialise number of points in ppd in each lattice vector direction
    cell%n_pt1 = n_pt1
    cell%n_pt2 = n_pt2
    cell%n_pt3 = n_pt3


    ! cks: initialise number of atoms
    cell%nat = nat

    ! ars : set pub_cell%class_nat
    pub_cell%nat_classical = class_nat

    ! cks: initialise number of NGWFs
    cell%num_ngwfs = num

    ! ndmh: initialise number of Conduction NGWFs
    cell%num_ngwfs_cond = num_cond

    ! ndmh: initialise number of Auxiliary NGWFs
    cell%num_ngwfs_aux = num_aux

    ! cks: initialise primitive lattice vectors with unit length
    cell%a1_unit = UNIT_VECTOR(cell%a1)
    cell%a2_unit = UNIT_VECTOR(cell%a2)
    cell%a3_unit = UNIT_VECTOR(cell%a3)

    ! cks: initialise reciprocal lattice vectors
    cell%b1=(2.0_DP*PI/( cell%a1.DOT.(cell%a2.CROSS.cell%a3) ) ) * &
         (cell%a2.CROSS.cell%a3)
    cell%b2=(2.0_DP*PI/( cell%a1.DOT.(cell%a2.CROSS.cell%a3) ) ) * &
         (cell%a3.CROSS.cell%a1)
    cell%b3=(2.0_DP*PI/( cell%a1.DOT.(cell%a2.CROSS.cell%a3) ) ) * &
         (cell%a1.CROSS.cell%a2)

    ! cks: initialise grid point weights
    cell%weight=cell%d1*cell%d2*cell%d3 * &
         ( (cell%a1_unit.CROSS.cell%a2_unit).DOT.cell%a3_unit )


    ! cks: initialise number of ppds in each lattice vector direction
    cell%n_ppds_a1=nint(geometry_MAGNITUDE(cell%a1) / &
         ( real(cell%n_pt1,DP)*cell%d1) )
    cell%n_ppds_a2=nint(geometry_MAGNITUDE(cell%a2) / &
         ( real(cell%n_pt2,DP)*cell%d2) )
    cell%n_ppds_a3=nint(geometry_MAGNITUDE(cell%a3) / &
         ( real(cell%n_pt3,DP)*cell%d3) )
    cell%n_ppds=cell%n_ppds_a1*cell%n_ppds_a2*cell%n_ppds_a3

    ! cks: initialise total number of points per ppd
    cell%n_pts = n_pt1 * n_pt2 * n_pt3

    ! cks: initialise total number of points in simulation cell in each
    !      lattice vector direction
    cell%total_pt1 = cell%n_ppds_a1 * cell%n_pt1
    cell%total_pt2 = cell%n_ppds_a2 * cell%n_pt2
    cell%total_pt3 = cell%n_ppds_a3 * cell%n_pt3

    ! cks: initialise total number of different atomic pseudopotential species
    ! cks, 25/1/2004: A distint "species" is determined by the
    ! cks: pseudopotential name rather than the atomic number
    cell%num_pspecies = 1
    do atom=2,cell%nat

       new_species = .true.
       do row=1,atom-1
          if (elements(atom)%pseudo_name == elements(row)%pseudo_name) &
               new_species = .false.
       end do

       if (new_species) cell%num_pspecies = cell%num_pspecies + 1

    end do

    ! qoh: initialise total number of different atomic species
    ! qoh: (as listed in % block species).

    cell%num_species = 1
    do atom=2,cell%nat

       new_species = .true.
       do row=1,atom-1
          if (elements(atom)%species_number == elements(row)%species_number) &
               new_species = .false.
       end do

       if (new_species) cell%num_species = cell%num_species + 1

    end do

  end subroutine simulation_cell_initialise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_nullify(cell)
    !=========================================================================!
    ! Nullify pointers in cell                                                !
    !-------------------------------------------------------------------------!
    ! Arash A Mostofi, 2004                                                   !
    !=========================================================================!

    implicit none

    type(castep_cell), intent(out) :: cell

    ! Nullify pointers
    nullify(cell%ionic_positions)
    nullify(cell%num_ions_in_species)
    nullify(cell%ionic_charge)
    nullify(cell%species_mass)
    nullify(cell%species_symbol)

    return
  end subroutine castep_cell_nullify

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_copy(in_cell,out_cell)
    !=========================================================================!
    ! Copy CASTEP cell data                                                   !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi 2004                                         !
    !=========================================================================!
    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout
    use utils, only : utils_alloc_check
    implicit none

    type(castep_cell), intent(inout) :: out_cell
    type(castep_cell), intent(in) :: in_cell

    ! <<< local variables >>>
    integer :: i,j,ierr

    ! Check that 'in_cell' pointers are associated - they should be! - abort if not.
    if (.not.associated(in_cell%ionic_positions)) then
       if (pub_on_root) write(stdout,'(a)') &
            'Error in ONETEP: in_cell%ionic_positions not associated in castep_cell_copy'
       call comms_abort
    end if
    if (.not.associated(in_cell%num_ions_in_species)) then
       if (pub_on_root) write(stdout,'(a)') &
            'Error in ONETEP: in_cell%ionic_positions not associated in castep_cell_copy'
       call comms_abort
    endif
    if (.not.associated(in_cell%ionic_charge)) then
       if (pub_on_root) write(stdout,'(a)') &
            'Error in ONETEP: in_cell%ionic_positions not associated in castep_cell_copy'
       call comms_abort
    endif
    if (.not.associated(in_cell%species_mass)) then
       if (pub_on_root) write(stdout,'(a)') &
            'Error in ONETEP: in_cell%species_mass not associated in castep_cell_copy'
       call comms_abort
    endif
    if (.not.associated(in_cell%species_symbol)) then
       if (pub_on_root) write(stdout,'(a)') &
            'Error in ONETEP: in_cell%species_symbol not associated in castep_cell_copy'
       call comms_abort
    endif

    ! Real and reciprocal space lattices
    do i=1,3
       do j=1,3
          out_cell%real_lattice(i,j)  = in_cell%real_lattice(i,j)
          out_cell%recip_lattice(i,j) = in_cell%recip_lattice(i,j)
       enddo
    enddo

    ! Cell volume
    out_cell%volume = in_cell%volume

    ! Check that 'out_cell' pointers are associated - they should be! - allocate if not.
    if (.not.associated(out_cell%ionic_positions)) then
       allocate(out_cell%ionic_positions(1:3,1:1,1:pub_cell%nat),stat=ierr)
       call utils_alloc_check('castep_cell_copy','cell%ionic_positions', ierr)
    endif
    if (.not.associated(out_cell%num_ions_in_species)) then
       allocate(out_cell%num_ions_in_species(1:pub_cell%nat),stat=ierr)
       call utils_alloc_check('castep_cell_copy','cell%num_ions_in_species', ierr)
    endif
    if (.not.associated(out_cell%ionic_charge)) then
       allocate(out_cell%ionic_charge(1:pub_cell%nat),stat=ierr)
       call utils_alloc_check('castep_cell_copy','cell%ionic_charge', ierr)
    endif
    if (.not.associated(out_cell%species_mass)) then
       allocate(out_cell%species_mass(1:pub_cell%nat),stat=ierr)
       call utils_alloc_check('castep_cell_copy','cell%species_mass', ierr)
    endif
    if (.not.associated(out_cell%species_symbol)) then
       allocate(out_cell%species_symbol(1:pub_cell%nat),stat=ierr)
       call utils_alloc_check('castep_cell_copy','cell%species_symbol', ierr)
    endif

    ! Ion data
    out_cell%num_ions            = in_cell%num_ions
    out_cell%num_species         = in_cell%num_species
    out_cell%max_ions_in_species = in_cell%max_ions_in_species
    do i=1,in_cell%num_ions
       out_cell%ionic_positions(1,1,i) = in_cell%ionic_positions(1,1,i)
       out_cell%ionic_positions(2,1,i) = in_cell%ionic_positions(2,1,i)
       out_cell%ionic_positions(3,1,i) = in_cell%ionic_positions(3,1,i)
       out_cell%num_ions_in_species(i) = in_cell%num_ions_in_species(i)
       out_cell%ionic_charge(i)        = in_cell%ionic_charge(i)
       out_cell%species_mass(i)        = in_cell%species_mass(i)
       out_cell%species_symbol(i)      = in_cell%species_symbol(i)
    end do


    return
  end subroutine castep_cell_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_frac_to_cart_vector(current_cell,frac,cart)
    !=========================================================================!
    ! Obtain absolute Cartesian co-ordinates from fractional ones             !
    !-------------------------------------------------------------------------!
    ! Arash A Mostofi, 2004                                                   !
    !=========================================================================!

    implicit none

    type(castep_cell), intent(in)            :: current_cell
    real(kind=dp), dimension(3), intent(in)  :: frac
    real(kind=dp), dimension(3), intent(out) :: cart

    ! <<< local variables >>>
    integer :: ii

    do ii=1,3
       cart(ii) = current_cell%real_lattice(1,ii)*frac(1) &
            + current_cell%real_lattice(2,ii)*frac(2) &
            + current_cell%real_lattice(3,ii)*frac(3)
    end do

    return
  end subroutine castep_cell_frac_to_cart_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_frac_to_cart_cell(current_cell,cart)
    !==========================================================================!
    ! Obtain absolute Cartesian co-ordinates from fractional ones              !
    ! This one acts on the whole array of coordinates as stored in cell object !
    !--------------------------------------------------------------------------!
    ! Victor Milman 21 Apr 2006                                                !
    !==========================================================================!

    implicit none

    type(castep_cell), intent(in)                :: current_cell
    real(kind=dp), dimension(:,:,:), intent(out) :: cart

    ! <<< local variables >>>
    integer :: ispec,iatom

    do ispec = 1,current_cell%num_species
       do iatom = 1,current_cell%num_ions_in_species(ispec)
          call castep_cell_frac_to_cart(current_cell, &
                current_cell%ionic_positions(:,iatom,ispec),cart(:,iatom,ispec))
       enddo
    enddo

    return
  end subroutine castep_cell_frac_to_cart_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_cart_to_frac(current_cell,cart,frac)
    !=========================================================================!
    ! Obtain fractional coordinates from Cartesian ones                       !
    !-------------------------------------------------------------------------!
    ! Arash A Mostofi, 2004                                                   !
    !=========================================================================!

    use constants, only: two_pi
    implicit none

    type(castep_cell), intent(in)           :: current_cell
    real(kind=dp), dimension(3), intent(in) :: cart
    real(kind=dp), dimension(3), intent(out):: frac

    ! <<< local variables >>>
    integer :: ii

    do ii=1,3
       frac(ii)=current_cell%recip_lattice(ii,1)*cart(1) &
            + current_cell%recip_lattice(ii,2)*cart(2) &
            + current_cell%recip_lattice(ii,3)*cart(3)
    end do

    frac(:) = frac(:) / two_pi

    return
  end subroutine castep_cell_cart_to_frac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_model_cell_changed(mdl)
    !=========================================================================!
    ! Sets castep_model data when ionic positions change                      !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, 2004                                        !
    !=========================================================================!

    implicit none

    type(castep_model), intent(inout) :: mdl

    mdl%found_ground_state = .false.
!!$    mdl%found_forces = .false.
!!$    mdl%bfgs_optimisation = .false.

    return
  end subroutine castep_model_cell_changed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_dealloc(cell)
    !=========================================================================!
    ! Dellocate CASTEP cell data                                              !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi 2004                                         !
    !=========================================================================!
    use comms, only: pub_on_root, comms_abort
    use utils, only: utils_dealloc_check

    implicit none

    type(castep_cell), intent(inout) :: cell

    ! <<< local variables >>>
    integer :: ierr

    deallocate(cell%species_symbol,stat=ierr)
    call utils_dealloc_check('castep_cell_copy','cell%species_symbol', ierr)
    deallocate(cell%species_mass,stat=ierr)
    call utils_dealloc_check('castep_cell_copy','cell%species_mass', ierr)
    deallocate(cell%ionic_positions,stat=ierr)
    call utils_dealloc_check('castep_cell_copy','cell%ionic_positions', ierr)
    deallocate(cell%num_ions_in_species,stat=ierr)
    call utils_dealloc_check('castep_cell_copy','cell%num_ions_in_species', ierr)
    deallocate(cell%ionic_charge,stat=ierr)
    call utils_dealloc_check('castep_cell_copy','cell%ionic_charge', ierr)


    return
  end subroutine castep_cell_dealloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_alloc(cell)
    !=========================================================================!
    ! Allocate CASTEP cell data                                               !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi 2004                                         !
    !=========================================================================!
    use comms, only: pub_on_root, comms_abort
    use constants, only: stdout
    use utils, only: utils_alloc_check

    implicit none

    type(castep_cell), intent(inout) :: cell

    ! <<< local variables >>>
    integer :: ierr

    if  ((associated(cell%ionic_positions)) .or. &
         (associated(cell%num_ions_in_species)) .or. &
         (associated(cell%ionic_charge)) .or. &
         (associated(cell%species_mass)) .or. &
         (associated(cell%species_symbol))) then
       if (pub_on_root) write (stdout,'(a)') 'ERROR in ONETEP: cell already &
            &allocated'
       call comms_abort
    end if

    ! For simplicity each atom has its own species
    allocate(cell%ionic_positions(1:3,1:1, 1: &
         (pub_cell%nat +pub_cell%nat_classical)),stat=ierr)
    call utils_alloc_check('castep_cell_copy','cell%ionic_positions', ierr)
    allocate(cell%num_ions_in_species(1: &
         (pub_cell%nat +pub_cell%nat_classical)),stat=ierr)
    call utils_alloc_check('castep_cell_copy','cell%num_ions_in_species', ierr)
    allocate(cell%ionic_charge(1: &
         (pub_cell%nat +pub_cell%nat_classical)),stat=ierr)
    call utils_alloc_check('castep_cell_copy','cell%ionic_charge', ierr)
    allocate(cell%species_mass(1: &
         (pub_cell%nat+pub_cell%nat_classical)),stat=ierr)
    call utils_alloc_check('castep_cell_copy','cell%species_mass', ierr)
    allocate(cell%species_symbol(1: &
         pub_cell%nat+pub_cell%nat_classical),stat=ierr)
    call utils_alloc_check('castep_cell_copy','cell%species_symbol', ierr)

    return
  end subroutine castep_cell_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_model_alloc(mdl,ndim,constrain_ions,do_lbfgs)
    !=========================================================================!
    ! Allocate CASTEP model data                                              !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi 2004                                         !
    !=========================================================================!

    use constants, only: stdout
    use comms, only: pub_on_root, comms_abort
    use utils, only : utils_alloc_check
    use rundat, only : geom_lbfgs_block_length

    implicit none

    type(castep_model), intent(inout) :: mdl
    integer, intent(in)               :: ndim
    logical, intent(in)               :: constrain_ions
    logical, intent(in)               :: do_lbfgs

    ! <<< local variables >>>
    integer :: length
    integer :: ierr

    ! Allocate mdl%cell
    call castep_cell_alloc(mdl%cell)

    ! Allocate mdl%orig_cell
    call castep_cell_alloc(mdl%orig_cell)

    ! Allocate mdl%ref_cell
    call castep_cell_alloc(mdl%ref_cell)

    ! Allocate mdl%reac_cell
    call castep_cell_alloc(mdl%reac_cell)

    ! Allocate mdl%prod_cell
    call castep_cell_alloc(mdl%prod_cell)

    ! Allocate mdl%intm_cell
    call castep_cell_alloc(mdl%intm_cell)

    ! Allocate mdl%forces
    ! For simplicity each atom has its own species
    allocate(mdl%forces(1:3,1:1,1:pub_cell%nat),stat=ierr)
    call utils_alloc_check('castep_model_alloc','mdl%forces', ierr)

    ! Allocate mdl%orig_forces if constrain_ions is TRUE
    ! For simplicity each atom has its own species
    if (constrain_ions) then
       allocate(mdl%orig_forces(1:3,1:1,1:pub_cell%nat),stat=ierr)
       call utils_alloc_check('castep_model_alloc','mdl%orig_forces', ierr)
    endif

!    if(present(method_string)) then
!       if((method_string(1).eq.'l').or.(method_string(1).eq.'L')) then
!          lbfgs=.true.
!       else if((method_string(1).eq.'b').or.(method_string(1).eq.'B')) then
!          lbfgs=.false.
!       else
!          write(stdout,*) "castep_model_alloc : method_string argument unrecognised"
!          call comms_abort
!       end if
!    else
!       lbfgs=.false.
!    end if

    ! Allocate mdl%bfgs_inv_Hessian
    if (((.not.do_lbfgs).or.(mdl%bfgs_optimisation)).and.pub_on_root) then
       allocate(mdl%bfgs_inv_Hessian(1:ndim,1:ndim),stat=ierr)
       call utils_alloc_check('castep_model_alloc','mdl%bfgs_inv_Hessian', ierr)
    end if

    if ((do_lbfgs.or.mdl%lbfgs_optimisation).and.pub_on_root) then
!       if(geom_lbfgs_max_updates.eq.0) then
!          length=geom_lbfgs_block_length+mdl%lbfgs_num_updates
!       else
!          length=geom_lbfgs_block_length
!       end if
       length=max(mdl%lbfgs_num_updates,geom_lbfgs_block_length) ! might be too large, but this will be sorted out after...
       allocate(mdl%lbfgs_position_updates(1:ndim,1:length),stat=ierr)
       call utils_alloc_check('castep_model_alloc','mdl%lbfgs_position_updates', ierr)
       allocate(mdl%lbfgs_gradient_updates(1:ndim,1:length),stat=ierr)
       call utils_alloc_check('castep_model_alloc','mdl%lbfgs_gradient_updates', ierr)
    end if

    return
  end subroutine castep_model_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_model_dealloc(mdl,do_lbfgs)
    !=========================================================================!
    ! Deallocate CASTEP model data                                            !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi 2004                                         !
    !=========================================================================!

    use comms, only : pub_on_root, comms_abort
    use utils, only : utils_dealloc_check

    implicit none

    type(castep_model), intent(inout) :: mdl
    logical :: do_lbfgs

    ! <<< local variables >>>
    integer :: ierr

    ! Deallocate mdl%bfgs_inv_Hessian
    if ((.not.do_lbfgs).and.pub_on_root) then
       deallocate(mdl%bfgs_inv_Hessian,stat=ierr)
       call utils_dealloc_check('castep_model_dealloc','mdl%bfgs_inv_Hessian', ierr)
    end if

    if (do_lbfgs.and.pub_on_root) then
       deallocate(mdl%lbfgs_position_updates,stat=ierr)
       call utils_dealloc_check('castep_model_dealloc','mdl%lbfgs_position_updates', ierr)
       deallocate(mdl%lbfgs_gradient_updates,stat=ierr)
       call utils_dealloc_check('castep_model_dealloc','mdl%lbfgs_gradient_updates', ierr)
    end if

    ! Dellocate mdl%orig_forces if constrain_ions is TRUE
    if (associated(mdl%orig_forces)) then
       deallocate(mdl%orig_forces,stat=ierr)
       call utils_dealloc_check('castep_model_dealloc','mdl%orig_forces', ierr)
    endif

    ! Dellocate mdl%forces
    deallocate(mdl%forces,stat=ierr)
    call utils_dealloc_check('castep_model_dealloc','mdl%forces', ierr)

    ! Deallocate mdl%intm_cell
    call castep_cell_dealloc(mdl%intm_cell)

    ! Deallocate mdl%prod_cell
    call castep_cell_dealloc(mdl%prod_cell)

    ! Deallocate mdl%reac_cell
    call castep_cell_dealloc(mdl%reac_cell)

    ! Deallocate mdl%ref_cell
    call castep_cell_dealloc(mdl%ref_cell)

    ! Deallocate mdl%orig_cell
    call castep_cell_dealloc(mdl%orig_cell)

    ! Deallocate mdl%cell
    call castep_cell_dealloc(mdl%cell)

    return
  end subroutine castep_model_dealloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_model_nullify(model)
    !=========================================================================!
    ! Nullify pointers in model                                               !
    !-------------------------------------------------------------------------!
    ! Arash A Mostofi, 2005                                                   !
    !=========================================================================!

    implicit none

    type(castep_model), intent(inout) :: model

    ! Nullify pointers
    nullify(model%forces)
    nullify(model%orig_forces)
    nullify(model%bfgs_inv_Hessian)
    nullify(model%lbfgs_position_updates)
    nullify(model%lbfgs_gradient_updates)
    call castep_cell_nullify(model%cell)
    call castep_cell_nullify(model%orig_cell)
    call castep_cell_nullify(model%ref_cell)
    call castep_cell_nullify(model%reac_cell)
    call castep_cell_nullify(model%prod_cell)
    call castep_cell_nullify(model%intm_cell)

    return
  end subroutine castep_model_nullify

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_cart_lattice_to_abc(cell,a,b,c,alpha,beta,gamma)

    !--------------------------------------------------------------------------!
    ! Purpose:                                                                 !
    ! This routine generates the lattice vectors of `cell' and                 !
    ! stores them in a,b,c,alpha,beta,gamma                                    !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! cell(in): cell from which real_lattice will be used                      !
    ! a(out): a lattice parameter                                              !
    ! b(out): b lattice parameter                                              !
    ! c(out): c lattice parameter                                              !
    ! alpha(out): angle between b and c                                        !
    ! beta(out) : angle between c and a                                        !
    ! gamma(out): angle between a and b                                        !
    !--------------------------------------------------------------------------!
    ! Parent module variables used: unit_cell                                  !
    !--------------------------------------------------------------------------!
    ! Modules used: none                                                       !
    !--------------------------------------------------------------------------!
    ! Key internal variables: none                                             !
    !--------------------------------------------------------------------------!
    ! Necessary conditions: none                                               !
    !--------------------------------------------------------------------------!
    ! Transferred from CASTEP by Victor Milman 20/04/06                        !
    !--------------------------------------------------------------------------!

    use constants, only : pi

    implicit none

    type(castep_cell), intent(in)::cell    ! The cell
    real(kind=dp), intent(out)::a,b,c,alpha,beta,gamma  ! The lattice parameters

    ! Calculate a
    a=sqrt(cell%real_lattice(1,1)**2+ &
         cell%real_lattice(1,2)**2+ &
         cell%real_lattice(1,3)**2)

    ! Calculate b
    b=sqrt(cell%real_lattice(2,1)**2+ &
         cell%real_lattice(2,2)**2+ &
         cell%real_lattice(2,3)**2)

    ! Calculate c
    c=sqrt(cell%real_lattice(3,1)**2+ &
         cell%real_lattice(3,2)**2+ &
         cell%real_lattice(3,3)**2)

    ! Calculate cos(alpha), temp stored in alpha
    alpha=(cell%real_lattice(2,1)*cell%real_lattice(3,1)+ &
         cell%real_lattice(2,2)*cell%real_lattice(3,2)+ &
         cell%real_lattice(2,3)*cell%real_lattice(3,3))/(b*c)

    ! Take acos and tranfrom into degrees
    alpha=180.0_dp*acos(alpha)/pi

    ! Calculate cos(beta), temp stored in beta
    beta =(cell%real_lattice(3,1)*cell%real_lattice(1,1)+ &
         cell%real_lattice(3,2)*cell%real_lattice(1,2)+ &
         cell%real_lattice(3,3)*cell%real_lattice(1,3))/(c*a)

    ! Take acos and tranfrom into degrees
    beta =180.0_dp*acos(beta)/pi

    ! Calculate cos(gamma), temp stored in gamma
    gamma=(cell%real_lattice(1,1)*cell%real_lattice(2,1)+ &
         cell%real_lattice(1,2)*cell%real_lattice(2,2)+ &
         cell%real_lattice(1,3)*cell%real_lattice(2,3))/(a*b)

    ! Take acos and tranfrom into degrees
    gamma=180.0_dp*acos(gamma)/pi

    return

  end subroutine castep_cell_cart_lattice_to_abc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine castep_cell_2_elements(cell,elements)
    !=========================================================================!
    ! Update elements with new ionic positions from castep cell               !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi, 2004                                        !
    ! AAM: Modified July 2005 to work with negative fractional co-ordinates   !
    !=========================================================================!

    use ion, only : element
    implicit none

    type(castep_cell), intent(in) :: cell
    type(element), intent(inout)  :: elements(1:pub_cell%nat)

    ! <<< local variable >>>
    integer :: atom,k
    real(kind=dp) :: frac(1:3), cart(1:3)

    do atom=1,pub_cell%nat
       frac(1:3) = cell%ionic_positions(1:3,1,atom)
       do k=1,3
          ! rationalise fractional positions to the interval [0,1]
          if ( frac(k).lt.0.0_dp ) then
             frac(k) = frac(k) + 1.0_dp
          else if ( frac(k).gt.1.0_dp ) then
             frac(k) = frac(k) - 1.0_dp
          else
             continue
          endif
       enddo

       call castep_cell_frac_to_cart(cell,frac,cart)

       elements(atom)%centre%x = cart(1)
       elements(atom)%centre%y = cart(2)
       elements(atom)%centre%z = cart(3)
    enddo

    return
  end subroutine castep_cell_2_elements

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine copy_to_castep_cell(current_cell,elements, classical_elements)
    !=========================================================================!
    ! Create and fill castep cell from pub_cell data structures               !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman 27/03/06                                       !
    ! Modified by Chris-Kriton Skylaris to inlcude classical info, 2010/08/30.!
    !=========================================================================!

    use comms, only : pub_on_root,comms_abort
    use constants, only: DP, periodic_table_name, two_pi, &
         periodic_table_mass, stdout
    use geometry, only: operator(.dot.)
    use ion, only: element
    implicit none

    type(castep_cell), intent(inout) :: current_cell
    type(element), intent(in)        :: elements(1:pub_cell%nat)
    type(element), intent(in), optional :: &
         classical_elements(1:pub_cell%nat_classical)

    ! local variables
    integer :: ion_i, ion_k
    logical :: found_ion

    !***********************************************************************!
    !         Copy data from ONETEP structures into CASTEP structures       !
    !***********************************************************************!

    ! Nullify pointers in current_cell
    call castep_cell_nullify(current_cell)

    ! Allocate current_cell
    call castep_cell_alloc(current_cell)

    ! Real space lattice
    current_cell%real_lattice(1,1) = pub_cell%a1%x
    current_cell%real_lattice(1,2) = pub_cell%a1%y
    current_cell%real_lattice(1,3) = pub_cell%a1%z
    current_cell%real_lattice(2,1) = pub_cell%a2%x
    current_cell%real_lattice(2,2) = pub_cell%a2%y
    current_cell%real_lattice(2,3) = pub_cell%a2%z
    current_cell%real_lattice(3,1) = pub_cell%a3%x
    current_cell%real_lattice(3,2) = pub_cell%a3%y
    current_cell%real_lattice(3,3) = pub_cell%a3%z

    ! Reciprocal lattice
    current_cell%recip_lattice(1,1) = pub_cell%b1%x
    current_cell%recip_lattice(1,2) = pub_cell%b1%y
    current_cell%recip_lattice(1,3) = pub_cell%b1%z
    current_cell%recip_lattice(2,1) = pub_cell%b2%x
    current_cell%recip_lattice(2,2) = pub_cell%b2%y
    current_cell%recip_lattice(2,3) = pub_cell%b2%z
    current_cell%recip_lattice(3,1) = pub_cell%b3%x
    current_cell%recip_lattice(3,2) = pub_cell%b3%y
    current_cell%recip_lattice(3,3) = pub_cell%b3%z

    ! Cell volume
    current_cell%volume = &
         pub_cell%a1%x*(pub_cell%a2%y*pub_cell%a3%z-pub_cell%a2%z*pub_cell%a3%y) + &
         pub_cell%a1%y*(pub_cell%a2%z*pub_cell%a3%x-pub_cell%a2%x*pub_cell%a3%z) + &
         pub_cell%a1%z*(pub_cell%a2%x*pub_cell%a3%y-pub_cell%a2%y*pub_cell%a3%x)

    ! Ion data plus classical ion data
    current_cell%num_ions = pub_cell%nat +pub_cell%nat_classical
    current_cell%num_species = pub_cell%nat +pub_cell%nat_classical
    current_cell%max_ions_in_species = 1 ! Give each ion its own species for now

    ! ONETEP ionic positions are in absolute Cartesian coordinates (in bohr)
    ! CASTEP uses fractional co-ordinates
    do ion_i=1,pub_cell%nat
       current_cell%ionic_positions(1,1,ion_i) = &
            (pub_cell%b1.dot.elements(ion_i)%centre) / two_pi
       current_cell%ionic_positions(2,1,ion_i) = &
            (pub_cell%b2.dot.elements(ion_i)%centre) / two_pi
       current_cell%ionic_positions(3,1,ion_i) = &
            (pub_cell%b3.dot.elements(ion_i)%centre) / two_pi
       current_cell%num_ions_in_species(ion_i) = 1
       current_cell%ionic_charge(ion_i)        = elements(ion_i)%ion_charge
       current_cell%species_symbol(ion_i)      = & ! Transform species_symbol
            species_symbol(elements(ion_i)%symbol) ! characters to standard form
    enddo

    ! cks: do the same for "classical" atoms
    ! ONETEP ionic positions are in absolute Cartesian coordinates (in bohr)
    ! CASTEP uses fractional co-ordinates
    if (present(classical_elements)) then
       do ion_i =pub_cell%nat +1, pub_cell%nat +pub_cell%nat_classical
          current_cell%ionic_positions(1,1,ion_i) = &
               (pub_cell%b1.dot.classical_elements(ion_i-pub_cell%nat)%centre)/&
               two_pi
          current_cell%ionic_positions(2,1,ion_i) = &
               (pub_cell%b2.dot.classical_elements(ion_i-pub_cell%nat)%centre)/&
               two_pi
          current_cell%ionic_positions(3,1,ion_i) = &
               (pub_cell%b3.dot.classical_elements(ion_i-pub_cell%nat)%centre)/&
               two_pi
          current_cell%num_ions_in_species(ion_i) = 1
          current_cell%ionic_charge(ion_i)        = &
               classical_elements(ion_i -pub_cell%nat)%ion_charge
          current_cell%species_symbol(ion_i)      = & ! Transform species_symbol
               species_symbol(classical_elements(ion_i-pub_cell%nat)%symbol)
                                                   ! characters to standard form
       enddo
    endif

    ! Ionic masses
    do ion_i=1,pub_cell%nat
       found_ion =.false.
       periodic_table: do ion_k=1,size(periodic_table_mass)
          if ( current_cell%species_symbol(ion_i).eq.periodic_table_name(ion_k) ) then
             current_cell%species_mass(ion_i) = periodic_table_mass(ion_k)
             found_ion = .true.
             exit periodic_table
          endif
       enddo periodic_table
       if (.not.found_ion) then
          if (pub_on_root) write(stdout,'(3a)') &
               'Error in ONETEP: species_symbol ',current_cell%species_symbol(ion_i),&
               ' in geometry_optimise not recognised'
          call comms_abort
       endif
    end do
    !*********************************************************************!
    !                       Cell data copy completed                      !
    !*********************************************************************!

  end subroutine copy_to_castep_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function species_symbol(string)
    !=========================================================================!
    ! Transform the species symbol character to a standard form (first letter !
    ! uppercase, second lowercase)                                            !
    !-------------------------------------------------------------------------!
    ! Written by Arash A Mostofi 2004                                         !
    !=========================================================================!
    use comms, only : pub_on_root,comms_abort
    use constants, only: stdout
    implicit none

    character(*), intent(in) :: string
    character(2)             :: species_symbol

    ! <<< local variables >>>
    integer :: iA,iaa,iZ,izz,ishift,ic,ln

    iaa = ichar('a')
    izz = ichar('z')
    iA = ichar('A')
    iZ = ichar('Z')
    ishift = iA-iaa

    ln = len(string)
    if (ln.lt.1 .or. ln.gt.2) then
       if (pub_on_root) write(stdout,'(3a)') &
            'Error in ONETEP: incorrect length of element symbol:',string,'.'
       call comms_abort
    endif

    species_symbol(1:ln) = string(1:ln)
    if (ln.lt.2) species_symbol(2:2)=' '

    ! Make first character uppercase
    ic = ichar(species_symbol(1:1))
    if ((ic.ge.iaa).and.(ic.le.izz)) species_symbol(1:1) = char(ishift+ic)

    ! Make second character lowercase
    if (species_symbol(2:2).ne.' ') then
       ic = ichar(species_symbol(2:2))
       if ((ic.ge.iA).and.(ic.le.iZ)) species_symbol(2:2) = char(ic-ishift)
    endif

    return
  end function species_symbol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function minimum_image_distance(r1,r2)
    !=========================================================================!
    ! Returns the distance between two points r1, r2 in the minimum image     !
    ! convention.                                                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   r1, r2 (input): the coordinates of the points in question.            !
    ! Returns:                                                                !
    !   |r1 - r2| in the minimum image convention.                            !
    ! Caveats, preconditions:                                                 !
    !   Both r1 and r2 are assumed to lie within the simulation cell.         !
    !   While this routine will sometimes work if the above is not satisified,!
    !   this should not be relied on.                                         !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 05/2010                                      !
    !=========================================================================!

    use geometry, only: POINT, operator(+), operator(-), operator(*), magnitude

    ! jd: Arguments
    type(POINT), intent(in) :: r1, r2

    ! jd: Local variables
    type(POINT) :: rvec0, rvec
    integer :: i1, i2, i3
    real(kind=DP) :: r, rmin

    !------------------------------------------------------------------------

    ! jd: Calculate the original displacement between the two
    rvec0 = r2-r1

    ! jd: Check 27 candidate displacement vectors to see which one is the
    !     shortest one, including the original one. This is the simplest way
    !     in non-orthorombic cells. Note that this fails if r1 or r2 lie
    !     sufficiently far away from the cell
    rmin = huge(1.0_DP)
    do i1 = -1, 1
       do i2 = -1, 1
          do i3 = -1, 1
             rvec = rvec0 + &
                  real(i1,kind=DP) * pub_cell%a1 + &
                  real(i2,kind=DP) * pub_cell%a2 + &
                  real(i3,kind=DP) * pub_cell%a3
             r = magnitude(rvec)
             if (r < rmin) rmin=r
          end do
       end do
    end do

    minimum_image_distance = rmin

  end function minimum_image_distance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module simulation_cell
