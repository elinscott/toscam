module xc_component

    use augmentation, only: augmentation_density_on_grid, &
         augmentation_screen_dij
    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP
    use cutoff_coulomb, only: cutoff_coulomb_hartree
    use enthalpy, only: enthalpy_terms
    use function_basis, only: FUNC_BASIS
    use hartree, only: hartree_on_grid, hartree_via_multigrid
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_coulomb_cutoff, pub_nlcc, pub_paw, &
         pub_multigrid_hartree, pub_aug_den_dim, pub_aug, pub_usp, &
         pub_nhat_in_xc, pub_external_pressure
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check
!    use xc, only: xc_energy_potential

    implicit none

    private

    public :: xc_energy_potential_comp
    public :: xc_init_comp
  
    real(kind=DP), parameter :: dentol = 1.0e-15_DP
    real(kind=DP), parameter :: THIRD = 1.0_DP / 3.0_DP
    real(kind=DP), parameter :: FTHRD = 4.0_DP / 3.0_DP
  
    ! Parsed version of public "xc_functional" variable
  character(len=80) :: functional

  ! Fraction of Hartree-Fock exchange
  real(kind=DP), save, public :: pub_hfxfraction

  ! Whether gradient corrections are required
  logical, save, public :: pub_xc_gradient_corrected

  ! Private internal variables used by the module
#ifdef LIBXC
  logical :: use_libxc
  logical :: libxc_gc
  type(xc_f90_pointer_t) :: libxc_func(2)
  type(xc_f90_pointer_t) :: libxc_info(2)
#endif

    ! Arguments
    ! Local Variables
    integer :: ierr   ! error flag
    integer :: is
    real(kind=DP), allocatable, dimension(:,:,:,:,:) :: nhat_den_grad


! ====================== END XC ENERGY =======================

contains

  subroutine xc_energy_potential_comp(density_fine_1,density_fine_2,xc_energy, &
      grid, nhat_dim,nhat_den_grad,tddftxc)

    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the exchange-correlation energy and the exchange-correlation potential,   !
    ! according to the functional defined by the input file.                    !
    !---------------------------------------------------------------------------!
    ! Arguments:                                                                !
    !    density_fine (in) : Input density on the fine grid                     !
    !    grid (in)         : GRID_INFO grid definition                          !
    !    xc_energy (out)   : Exchange-correlation energy                        !
    !    xc_pot_fine (out) : Exchange-correlation potential on grid             !
    !    nhat_dim (in)     : Number of components of compensation density (if   !
    !                        present - ignored if not present)                  !
    !    nhat_den_grad (in): Compensation density and Compensation density      !
    !                        gradient (if required).                            !
    !    tddftxc (in)      : Switch to override exchange-correlation functional !
    !                        with TDDFT exchange-correlation functional         !
    !---------------------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004                         !
    ! Addition of hybrid functionals and modification to use new subroutines by !
    ! Quintin Hill 11/03/2009.                                                  !
    !===========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_abort, pub_on_root
    use constants, only: stdout, DP, PI, VERBOSE
    use rundat, only: pub_output_detail
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check
!    use xc_funcs

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid    ! Grid definition
    real(kind=DP), intent(inout) :: density_fine_1(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
    real(kind=DP), intent(inout) :: density_fine_2(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
    real(kind=DP), intent(out) :: xc_energy
    integer, intent(in) :: nhat_dim
    real(kind=DP), optional, intent(in) :: nhat_den_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins,0:nhat_dim)
    logical, optional, intent(in) :: tddftxc

    ! Local Variables
#ifdef WIN32
  real(kind=DP), allocatable, target :: density_grad(:,:,:,:,:)
  real(kind=DP), allocatable, target :: density_aux(:,:,:,:)
  complex(kind=DP), allocatable, target :: recip_work(:,:,:,:)
#else
  real(kind=DP), allocatable :: density_grad(:,:,:,:,:)
  real(kind=DP), allocatable :: density_aux(:,:,:,:)
  complex(kind=DP), allocatable :: recip_work(:,:,:,:)
#endif
    logical          :: loc_tddftxc ! ddor: A local copy of tddftxc

    ! Start timer
    call timer_clock('xc_energy_potential_comp',1)

    ! ddor: Set default for loc_tddftxc
    loc_tddftxc = .false.
    if (present(tddftxc)) loc_tddftxc = tddftxc

    ! Check input density
    if (pub_output_detail == VERBOSE) call internal_check_density_comp

    ! Allocate workspace if required
#ifdef LIBXC
    if (pub_xc_gradient_corrected.or.use_libxc) then
        call internal_allocate_workspace_comp
    endif
#else
    if (pub_xc_gradient_corrected) call internal_allocate_workspace_comp
#endif

    ! Call appropriate routine according to functional type
    select case (functional)
    case ('CAPZ')
       !----------------------------------------------------------------------!
       ! Local density approximation - from Ceperley-Alder Monte Carlo data,  !
       !  and Gell-Mann-Brueckner expansion, parameterised by Perdew & Zunger !
       !----------------------------------------------------------------------!
       call xc_lda_comp(density_fine_1,density_fine_2,xc_energy,grid,xc_capz_c_point_comp)
    case ('VWN')
       !----------------------------------------------------------------------!
       ! Local density approximation - Vosko, Wilk & Nusair                   !
       !----------------------------------------------------------------------!
       call xc_lda_comp(density_fine_1,density_fine_2,xc_energy,grid,xc_vwn_c_point_comp)
    case ('PW92')
       !----------------------------------------------------------------------!
       ! Local density approximation - Perdew & Wang 1992 functional          !
       !----------------------------------------------------------------------!
       call xc_lda_comp(density_fine_1,density_fine_2,xc_energy,grid,xc_pw92_c_point_comp)
    case default
       !----------------------------------------------------------------------!
       ! Unrecognised functional type                                         !
       !----------------------------------------------------------------------!
       if (pub_on_root) write(stdout,'(2a)') 'WARNING: unknown exchange-&
            &correlation functional: ', trim(functional)
       xc_energy = 0.0_DP
    end select

    ! Free workspace if required
    if (pub_xc_gradient_corrected) call internal_deallocate_workspace_comp

    ! Stop timer
    call timer_clock('xc_energy_potential_comp',2)

  contains


    subroutine internal_allocate_workspace_comp

      !=========================================================================!
      ! This subroutine allocates workspace for the gradients of the density if !
      ! required.                                                               !
      !-------------------------------------------------------------------------!
      ! Originally written by Peter Haynes in February 2004.                    !
      ! Modified by Quintin Hill in November 2008 to use utils_alloc_check.     !
      ! Modified by Quintin Hill in March 2009 to use pub_cell%num_spins.       !
      ! Made subroutine-internal by Nicholas Hine, October 2010.                !
      !=========================================================================!

      implicit none

      ! Local Variables
      integer :: ierr    ! Error flag

      ! Allocate workspace for calculating gradients of the density
      allocate(density_grad(grid%ld1,grid%ld2,grid%max_slabs12,3, &
           pub_cell%num_spins),stat=ierr)
      call utils_alloc_check('internal_allocate_workspace comp (xc)', &
           'density_grad',ierr)

      ! Allocate workspace for the auxiliary density array
      allocate(density_aux(grid%ld1,grid%ld2,grid%max_slabs12, &
           pub_cell%num_spins),stat=ierr)
      call utils_alloc_check('internal_allocate_workspace comp (xc)', &
           'density_aux',ierr)

      ! Allocate reciprocal space workspace
      allocate(recip_work(grid%ld3,grid%ld2,grid%max_slabs23,3),stat=ierr)
      call utils_alloc_check('internal_allocate_workspace comp (xc)', &
           'recip_work',ierr)

    end subroutine internal_allocate_workspace_comp

    subroutine internal_deallocate_workspace_comp

      !========================================================================!
      ! This subroutine frees workspace for the gradients of the density.      !
      !------------------------------------------------------------------------!
      ! Originally written by Peter Haynes in February 2004.                   !
      ! Modified by Quintin Hill in November 2008 to use utils_dealloc_check.  !
      ! Made subroutine-internal by Nicholas Hine, October 2010.               !
      !========================================================================!

      implicit none

      ! Local Variables
      integer :: ierr    ! Error flag

      ! Deallocate reciprocal space workspace
      if (allocated(recip_work)) then
         deallocate(recip_work,stat=ierr)
         call utils_dealloc_check('internal_deallocate_workspace (xc)', &
              'recip_work', ierr)
      end if

      ! Deallocate workspace for calculating gradients of the density
      if (allocated(density_aux)) then
         deallocate(density_aux,stat=ierr)
         call utils_dealloc_check('internal_deallocate_workspace (xc)', &
              'density_aux', ierr)
      end if

      ! Deallocate workspace for calculating gradients of the density
      if (allocated(density_grad)) then
         deallocate(density_grad,stat=ierr)
         call utils_dealloc_check('internal_deallocate_workspace (xc)', &
              'density_grad', ierr)
      end if

    end subroutine internal_deallocate_workspace_comp

    subroutine internal_check_density_comp

      !------------------------------------------------------------------------!
      ! This subroutine prints a warning if the charge density of either spin  !
      ! component is negative at any point.                                    !
      !------------------------------------------------------------------------!

      use comms, only: comms_reduce, pub_on_root
      implicit none

      !------------------------
      ! Local variables:
      !------------------------

      integer :: i1,i2,islab12 ! Loop counters over real-space fine grid
      integer :: is            ! Loop counter over spin
      real(kind=DP) :: denmin  ! Minimum values of density

      spins: do is=1, pub_cell%num_spins

         denmin = density_fine_1(1,1,1,is)

         a3: do islab12=1,grid%num_my_slabs12
            a2: do i2=1,grid%n2
               a1: do i1=1,grid%n1

                  denmin = min(denmin,density_fine_1(i1,i2,islab12,is))

               end do a1
            end do a2
         end do a3

         call comms_reduce('MIN',denmin)

         if (pub_cell%num_spins == 1) then
            if (denmin <= -0.0001_DP .and. pub_on_root) &
                 write(stdout,'(a,f7.4)') 'WARNING: minimum value of charge &
                 &density: ',denmin
         else
            if (denmin <= -0.0001_DP .and. pub_on_root) &
                 write(stdout,'(a,i1,a,f7.4)') 'WARNING: minimum value of charge&
                 & density for spin ',is,': ',denmin
         end if

      end do spins

    end subroutine internal_check_density_comp

  end subroutine xc_energy_potential_comp
  
  subroutine xc_init_comp(use_tddftxc)

    !=========================================================================!
    ! This subroutine initialises the internal variables in the module.       !
    !-------------------------------------------------------------------------!
    ! Originally written by Peter Haynes in February 2004.                    !
    ! Modified by Quintin Hill to remove unnecessary internal uppercase       !
    ! subroutine in November 2008.                                            !
    ! Modified by Quintin Hill to cope with additional functionals and remove !
    ! mention of nspins.                                                      !
    ! Modified by Quintin Hill on 06/04/2009 to remove unnecessary            !
    ! initialisations, workspace allocation is already done elsewhere.        !
    ! Modified for TDDFT by David D. O'Regan in May 2009.                     !
    ! Modified to be called at start only by Nicholas Hine in October 2010.   !
    !=========================================================================!

    use rundat, only: xc_functional, tddft_xc_functional, &
         pub_libxc_x_func_id, pub_libxc_c_func_id, pub_aug
    use simulation_cell, only: pub_cell
    use utils, only: utils_abort
#ifdef LIBXC
    use xc_f90_types_m
    use xc_f90_lib_m
#endif

    implicit none

    ! Arguments
    ! ddor: Switch between DFT xc functional and TDDFT xc functional
    logical, intent(in), optional :: use_tddftxc

    ! Local Variables
    character(len=6), parameter :: gcfunctionals(14) &
         = (/'BLYP  ','PBE   ','PW91  ','REVPBE','RPBE  ','XLYP  ','WC    ', &
         'B1LYP ','B1PW91', 'B3LYP ','B3PW91','X3LYP ','PBE0  ','PBEX  '/)
#ifdef LIBXC
    integer :: libxc_family
    integer :: libxc_polarized
#endif

    ! Make a local copy of the xc_functional parameter from the input file
    ! ddor: loc_tddftxc switches between functionals defined in
    !       xc_functional and tddft_xc_functional
    if (present(use_tddftxc)) then
       if (use_tddftxc) then
          functional = tddft_xc_functional
       else
          functional = xc_functional
       end if
    else
       functional = xc_functional
    endif

    ! Apply defaults to functional type
    if (functional == 'LDA') functional = 'CAPZ'
    if (functional == 'GGA') functional = 'RPBE'

    ! ndmh: initialisations for LIBXC 'functional'
    if (functional == 'LIBXC') then

#ifdef LIBXC
       use_libxc = .true.

       ! ndmh: determined polarization
       if (pub_cell%num_spins==1) then
          libxc_polarized = XC_UNPOLARIZED
       else
          libxc_polarized = XC_POLARIZED
       end if

       pub_xc_gradient_corrected = .false.

       if (pub_libxc_x_func_id > 0) then
          call utils_abort('Error in xc_init_comp: libxc not supported')
!          call xc_f90_func_init(libxc_func(1), libxc_info(1), &
!               pub_libxc_x_func_id, libxc_polarized)
!          libxc_family = xc_f90_info_family(libxc_info(1))
!          libxc_gc = (libxc_family==XC_FAMILY_GGA) .or. &
!               (libxc_family==XC_FAMILY_HYB_GGA)
       end if
       if (pub_libxc_c_func_id > 0) then
          call utils_abort('Error in xc_init_comp: libxc not supported')
!          call xc_f90_func_init(libxc_func(2), libxc_info(2), &
!               pub_libxc_c_func_id, libxc_polarized)
!          libxc_family = xc_f90_info_family(libxc_info(2))
!          if (libxc_gc.neqv.((libxc_family==XC_FAMILY_GGA) .or. &
!               (libxc_family==XC_FAMILY_HYB_GGA))) then
!             call utils_abort('Error in xc_init: Mixing of LIBXC functionals &
!                  &with and without gradient corrections is not supported')
!          end if
       end if

       ! ndmh: check if this functional is gradient corrected
       pub_xc_gradient_corrected = libxc_gc

       if (libxc_family==XC_FAMILY_MGGA) then
          call utils_abort('Error in xc_init: Meta-GGA functionals &
               &not supported')
       end if
#else
       call utils_abort('Error in xc_init: Functional "LIBXC" &
            &specified, but LIBXC not present in build')
#endif

    else

#ifdef LIBXC
       use_libxc = .false.
#endif

       ! Determine whether functional requires gradients to be calculated
       pub_xc_gradient_corrected = (any(gcfunctionals == functional))

    end if

    if (pub_xc_gradient_corrected.and.pub_aug) then
       call utils_abort('Error in xc_init: Augmentation charges with &
            &gradient-corrected functionals are not yet supported')
    end if

  end subroutine xc_init_comp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lda_comp(density_fine_1,density_fine_2,xc_energy,grid,c_functional)

    !==========================================================================!
    ! This subroutine calculates the exchange-correlation energy and potential !
    ! according to the local density approximation (LDA), given the electronic !
    ! density on the fine grid. The correlation functional to be used is given !
    ! as an argument.                                                          !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in March 2009 based on existing xc_capz          !
    ! subroutine written by Peter Haynes in February 2004.                     !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP, PI
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    ! Grid definition
    type(GRID_INFO), intent(in) :: grid

    ! Input total density on the fine grid
    real(kind=DP), intent(in) :: density_fine_1(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)

    ! Input partial density on the fine grid
    real(kind=DP), intent(in) :: density_fine_2(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)

    ! Exchange-correlation energy
    real(kind=DP), intent(out) :: xc_energy

    interface
       subroutine c_functional(den1,den2,c_energy)
         use constants, only: DP, PI
         real(kind=DP), intent(in)  :: den1
         real(kind=DP), intent(in)  :: den2
         real(kind=DP), intent(out) :: c_energy
       end subroutine c_functional
    end interface

    ! Local variables:
    integer :: i1,i2,islab12   ! Fine grid loop counters
    real(kind=DP) :: ex        ! Exchange energy in Ha
    real(kind=DP) :: ec        ! Correlation energy in Ha
    real(kind=DP) :: den1      ! Total Charge density
    real(kind=DP) :: den2      ! Partial Charge density
    real(kind=DP) :: x_energy  ! Exchange energy at point
    real(kind=DP) :: c_energy  ! Correlation energy at point

    ex = 0.0_DP ; ec = 0.0_DP

    if (pub_cell%num_spins == 1) then

       ! qoh: Spin-unpolarised case:
       a3: do islab12=1,grid%num_my_slabs12
          a2: do i2=1,grid%n2
             a1: do i1=1,grid%n1

                den1 = density_fine_1(i1,i2,islab12,1)
                den2 = density_fine_2(i1,i2,islab12,1)

                call xc_lda_x_point_comp(den1,den2,x_energy)
                call c_functional(den1,den2,c_energy)

                ex = ex + x_energy
                ec = ec + c_energy

             end do a1
          end do a2
       end do a3

    end if

    ! Global reductions
    call comms_reduce('SUM',ex)
    call comms_reduce('SUM',ec)

    xc_energy = (ex + ec) * grid%weight

  end subroutine xc_lda_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lda_x_point_comp(den1, den2, x_energy)

    !==========================================================================!
    ! This subroutine calculates the local density appoximation (LDA) exchange !
    ! energy and potential at a point given the density at that point.         !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in March 2009 based on existing xc_capz          !
    ! subroutine written by Peter Haynes in February 2004.                     !
    !==========================================================================!

    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den1      ! Total density
    real(kind=DP), intent(in)  :: den2      ! Partial Density
    real(kind=DP), intent(out) :: x_energy ! Exchange energy

    real(kind=DP) :: rs
    real(kind=DP) :: epsx
    real(kind=DP), parameter :: xf = -0.458165293283142893475554_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    if (den1 > 0.0_DP) then   ! Positive charge density

       rs = (TONFPI / den1)**THIRD   ! Wigner-Seitz radius
       epsx = xf / rs
       x_energy = den2 * epsx
    else
       x_energy = 0.0_DP
    end if

  end subroutine xc_lda_x_point_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_capz_c_point_comp(den1,den2,c_energy)

    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the exchange-correlation energy and potential, according to the local     !
    ! density approximation (LDA), from the Monte Carlo data by Ceperley and    !
    ! Alder, the Gell-Mann & Brueckner expansion, and the parametrisation by    !
    ! Perdew & Zunger.                                                          !
    !                                                                           !
    ! References:                                                               !
    !  D M Ceperley & B J Alder, Phys. Rev. Lett. 45, 566 (1980)                !
    !  M Gell-Mann & K Brueckner, Phys. Rev. 106, 364 (1957)                    !
    !     see also N W Ashcroft & N D Mermin 'Solid State Physics',             !
    !              International Edition, p. 336                                !
    !  J P Perdew & A Zunger, Phys. Rev. B 23, 5048 (1981)                      !
    !---------------------------------------------------------------------------!
    ! Originally written by Peter Haynes, February 2004.                        !
    ! Split into point contribution by Quintin Hill on 09/03/2009.              !
    !===========================================================================!


    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den1
    real(kind=DP), intent(in)  :: den2
    real(kind=DP), intent(out) :: c_energy

    real(kind=DP) :: rs           ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs       ! Square root of rs
    real(kind=DP) :: logrs        ! Natural logarithm of rs
    real(kind=DP) :: epsc         ! Correlation energy per particle

    ! Constants
    real(kind=DP), parameter :: SSXTH = 7.0_DP / 6.0_DP
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    ! Constants for spin-unpolarised case:
    real(kind=DP), parameter :: g = -0.1423_DP
    real(kind=DP), parameter :: b1 = 1.0529_DP
    real(kind=DP), parameter :: b2 = 0.3334_DP
    real(kind=DP), parameter :: aa = 0.0311_DP
    real(kind=DP), parameter :: bb = -0.048_DP
    real(kind=DP), parameter :: cc = 0.0020_DP
    real(kind=DP), parameter :: dd = -0.0116_DP

    if (den1 > 0.0_DP) then   ! Positive charge density

       rs = (TONFPI / den1)**THIRD   ! Wigner-Seitz radius

       if (rs >= 1.0_DP) then  ! Low charge density

          sqrtrs = sqrt(rs)
          epsc = g / (1.0_DP + b1 * sqrtrs + b2 * rs)
          c_energy = den2 * epsc

       else  ! High charge density

          logrs = log(rs)
          epsc = aa * logrs + bb + &
               (cc * logrs + dd) * rs
          c_energy = den2 * epsc

       end if

    else   ! Negative charge density
       c_energy = 0.0_DP

    end if

  end subroutine xc_capz_c_point_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vwn_c_point_comp(den1,den2,c_energy)
    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the correlation energy and potential, according to the local density      !
    ! approximation (LDA), from the parameterisation by Vosko, Wilk and Nusair. !
    !                                                                           !
    ! Spin-unpolarized version.                                                 !
    !                                                                           !
    ! Reference:                                                                !
    !  [1] S H Vosko, L Wilk & M Nusair, Can. J. Phys. 58, 1200 (1980)          !
    !  [2] S H Vosko and L Wilk, Phys. Rev. B 22, 3812 (1980)                   !
    !                                                                           !
    !---------------------------------------------------------------------------!
    ! Written in 2011/05/31 by Jacek Dziedzic, basing on a template by Peter    !
    ! Haynes, modified by Quintin Hill.                                         !
    !===========================================================================!
    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den1
    real(kind=DP), intent(in)  :: den2
    real(kind=DP), intent(out) :: c_energy

    real(kind=DP) :: rs ! Wigner-Seitz radius corresponding to den

    ! Constants
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)
    real(kind=DP), parameter :: minus_fourpi_over_nine = -4.0_DP/9.0_DP * PI
    real(kind=DP), parameter :: MIN_DENS = 5.0D-13

    ! Local variables
    real(kind=DP) :: c_eps     ! Correlation energy density
    real(kind=DP) :: drs_drho  ! Derivative of rs wrt electronic density

    ! -------------------------------------------------------------------------

    ! We're only concerned about positive densities.
    ! Furthermore, to get perfect agreement with LIBXC, we ignore densities
    ! below MIN_DENS just like they do
    if (den1 > MIN_DENS) then

       rs = (TONFPI / den1)**THIRD ! the Wigner-Seitz radius

       ! Calculate energy and potential (eq. (4.4) in [1] and its derivative)
       call xc_vwn_eps_c_helper_comp(c_eps,rs,1)

       ! Energy density -> energy
       c_energy = c_eps * den2

    else
       c_energy = 0.0_DP
    end if

  end subroutine xc_vwn_c_point_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_vwn_eps_c_helper_comp(c_eps, rs, paramset)
    !===========================================================================!
    ! Given the Wigner-Seitz radius, rs, calculates the correlation energy      !
    ! density, according to expression (4.4) in [1]. The first derivative       !
    ! wrt rs is also calculated.                                                !
    !---------------------------------------------------------------------------!
    ! Arguments:                                                                !
    !   c_eps    (out): Returned energy density.                                !
    !   c_pot    (out): Returned potential.                                     !
    !   rs       (in):  Wigner-Seitz radius (a simple function of density)      !
    !   paramset (in):  Selects the VWN parameter set:                          !
    !                   1 - paramagnetic case                                   !
    !                   2 - ferromagnetic case                                  !
    !                   3 - calculation of alpha_c (cf. text below (9) in [2].  !
    !---------------------------------------------------------------------------!
    ! Reference:                                                                !
    !  [1] S H Vosko, L Wilk & M Nusair, Can. J. Phys. 58, 1200 (1980)          !
    !  [2] S H Vosko and L Wilk, Phys. Rev. B 22, 3812 (1980)                   !
    !---------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 01/06/2011.                                    !
    !===========================================================================!
    use constants, only: DP, PI
    use utils, only: utils_abort
    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: c_eps
    real(kind=DP), intent(in)  :: rs
    integer                    :: paramset

    ! jd: Parameters of the method, following [1], [2]
    !     Note A(:) are divided by 2 because of Ry->Ha conversion.
    !     This also applies to A(3), which becomes 1/(6 pi^2) rather than
    !     1/(3 pi^2)
    real(kind=DP), parameter :: As(3) = &
         (/ 0.0310907_DP, 0.01554535_DP, -0.01688686394038963_DP /)
    real(kind=DP), parameter :: bs(3) = (/ 3.72744_DP, 7.06042_DP, 1.13107_DP /)
    real(kind=DP), parameter :: cs(3) = (/ 12.9352_DP, 18.0578_DP, 13.0045_DP /)
    real(kind=DP), parameter :: x0s(3) = &
         (/ -0.10498_DP, -0.325_DP, -0.0047584_DP /)
    real(kind=DP), parameter :: Qs(3) = &
         (/ 6.15199082_DP, 4.73092691_DP, 7.123108918_DP /) ! sqrt(4c-b2)
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    ! jd: Internal variables
    real(kind=DP) :: x                 ! jd: =rs^(1/2), cf. [1], below [4.3]
    real(kind=DP) :: X_of_x, X_of_x0   ! jd:            cf. [1], below [4.4]
    real(kind=DP) :: den               ! jd: the density corresponding to rs
    real(kind=DP) :: ln1, ln2, arctan1 ! jd: Temporaries
    real(kind=DP) :: A,b,c,x0,Q        ! jd: Parameters for current paramset

    ! -------------------------------------------------------------------------

    if(paramset < 1 .or. paramset > 3) then
       call utils_abort("Bad paramset in xc_vwn_eps_c_helper_comp")
    endif

    A = As(paramset)
    b = bs(paramset)
    c = cs(paramset)
    x0 = x0s(paramset)
    Q = Qs(paramset)

    x = sqrt(rs) ! rs>0 is guaranteed by the caller
    den = TONFPI / (rs*rs*rs)
    X_of_x = x*x + b*x + c
    X_of_x0 = x0*x0 + b*x0 + c

    ln1 = log(x*x/X_of_x)
    ln2 = log((x-x0)*(x-x0)/X_of_x)
    arctan1 = atan(Q/(2.0_DP*x+b))

    ! Calculate energy density according to (4.4)
    c_eps = A* &
         (ln1 + 2.0_DP*b/Q * arctan1 - b*x0/X_of_x0 * &
         (ln2 + 2.0_DP*(b+2.0_DP*x0)/Q * arctan1))

  end subroutine xc_vwn_eps_c_helper_comp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_pw92_c_point_comp(den1,den2,c_energy)

    !===========================================================================!
    ! Given the electronic density on the fine grid, this subroutine calculates !
    ! the exchange-correlation energy and potential, according to the local     !
    ! density approximation (LDA), from the Monte Carlo data by Ceperley and    !
    ! Alder, the Gell-Mann & Brueckner expansion, and the parametrisation by    !
    ! Perdew & Wang.                                                            !
    !                                                                           !
    ! References:                                                               !
    !  D M Ceperley & B J Alder, Phys. Rev. Lett. 45, 566 (1980)                !
    !  M Gell-Mann & K Brueckner, Phys. Rev. 106, 364 (1957)                    !
    !     see also N W Ashcroft & N D Mermin 'Solid State Physics',             !
    !              International Edition, p. 336                                !
    !  J P Perdew and Y Wang, Phys. Rev. B 45, 13244 (1992)                     !
    !---------------------------------------------------------------------------!
    ! Originally written by Nicholas Hine, March 2011.                          !
    !===========================================================================!


    use constants, only: DP, PI
    implicit none

    real(kind=DP), intent(in)  :: den1
    real(kind=DP), intent(in)  :: den2
    real(kind=DP), intent(out) :: c_energy

    real(kind=DP) :: rs           ! Wigner-Seitz radius
    real(kind=DP) :: sqrtrs       ! Square root of rs
    real(kind=DP) :: ecunif,eurs

    ! Constants
    real(kind=DP), parameter :: TONFPI = 3.0_DP / (4.0_DP * PI)

    ! Constants for spin-unpolarised case:

    tol: if (den1 > dentol) then   ! Positive charge density
       rs = (TONFPI / den1)**THIRD   ! Wigner-Seitz radius
       sqrtrs = sqrt(rs)
       call xc_pw91_eq10_comp(0.0310907_DP,0.21370_DP,7.5957_DP, &
            3.5876_DP,1.6382_DP,0.49294_DP,rs,sqrtrs,ecunif,eurs)
       c_energy = den2 * ecunif

    else
       c_energy = 0.0_DP
    end if tol

  end subroutine xc_pw92_c_point_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine xc_pw91_eq10_comp(A,alpha1,beta1,beta2,beta3,beta4,rs,sqrtrs,G,dGdrs)

    !---------------------------------------------------------!
    ! This subroutine evaluates Eq. 10 of:                    !
    !  J P Perdew & Y Wang, Phys. Rev. B 45, 13244 (1992)     !
    ! for G(rs,A,alpha1,beta1,beta2,beta3,beta4,p=1)          !
    ! and its derivative with respect to rs.                  !
    !---------------------------------------------------------!

    implicit none

    !------------------------
    ! INPUT/OUTPUT variables:
    !------------------------

    real(kind=DP), intent(in)  :: A,alpha1,beta1,beta2,beta3,beta4,rs,sqrtrs
    real(kind=DP), intent(out) :: G,dGdrs

    !------------------------
    ! Local variables:
    !------------------------

    real(kind=DP) :: q0,rs32,rs2,q1,q2,q3

    q0 = -2.0_DP*A*(1.0_DP+alpha1*rs)
    rs32 = sqrtrs*rs
    rs2 = rs*rs
    q1 = 2.0_DP*A*(beta1*sqrtrs+beta2*rs+beta3*rs32+beta4*rs2)
    q2 = log(1.0_DP+1.0_DP/q1)
    G = q0*q2
    q3 = A*(beta1/sqrtrs+2.0_DP*beta2+3.0_DP*beta3*sqrtrs+4.0_DP*beta4*rs)
    dGdrs = -2.0_DP*A*alpha1*q2-q0*q3/(q1*q1+q1)

  end subroutine xc_pw91_eq10_comp


end module xc_component 
