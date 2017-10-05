! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!==============================================================================!
!                                                                              !
!  ONETEP exchange-correlation module: xc_mod.F90                              !
!                                                                              !
!  The subroutines in this file were written by Peter Haynes, February 2004    !
!                                                                              !
!  Theory of Condensed Matter, Cavendish Laboratory, University of Cambridge   !
!  Madingley Road, Cambridge CB3 0HE, UK                                       !
!                                                                              !
!  Addition of hybrid functionals, BLYP and XLYP, and major restructuring of   !
!  module by Quintin Hill in February/March 2009 under the supervision of      !
!  Chris-Kriton Skylaris.                                                      !
!  Note for hybrid functionals the Hartree-Fock exchange contribution to the   !
!  energy and potential is calculated and applied outside this subroutine.     !
!                                                                              !
!  Modified for TDDFT by David D. O'Regan in May 2009.                         !
!  Added support for libxc functionals, Nicholas Hine, July 2010.              !
!  Added support for radial xc calculations, Nicholas Hine, September 2010.    !
!------------------------------------------------------------------------------!
!                                                                              !
! The following exchange-correlation functionals (specified by the input file  !
! keyword XC_FUNCTIONAL) are available:                                        !
!                                                                              !
!  Local density approximation functionals:                                    !
!   CAPZ   - L(S)DA based on Ceperley-Alder Monte Carlo data parametrised by   !
!            Perdew & Zunger                                                   !
!   VWN    - L(S)DA based on parametrisation by Vosko, Wilk & Nusair           !
!   PW92   - L(S)DA based on Ceperley-Alder Monte Carlo data parametrised by   !
!            Perdew & Wang                                                     !
!                                                                              !
!  Generalised gradient approximation functionals:                             !
!   BLYP   - GG(S)A based on Becke 88 + LYP (Lee, Yang, Parr)                  !
!   PBE    - GG(S)A based on 1996 functional of Perdew, Burke Ernzerhof        !
!   REVPBE - GG(S)A based on 1998 revised PBE by Zhang & Yang                  !
!   RPBE   - GG(S)A based on 1999 RPBE by Hammer, Hansen & Norskov             !
!   PW91   - GG(S)A based on 1991 functional of Perdew & Wang                  !
!   XLYP   - GG(S)A based on 2004 functional of Xu and Goddard                 !
!   WC     - GG(S)A based on 2006 functional of Wu and Cohen                   !
!                                                                              !
!  Hybrid functionals:                                                         !
!   B1LYP  - hybrid based on 1997 functional of Adamo and Barone               !
!   B1PW91 - hybrid based on 1997 functional of Adamo and Barone               !
!   B3LYP  - hybrid based on 1993 functional of Becke using LYP instead of PW91!
!   B3PW91 - hybrid based on 1993 functional of Becke                          !
!   PBE0   - hybrid based on 1999 functional of Adamo, Cossi and Barone        !
!   X3LYP  - hybrid based on 2004 functional of Xu and Goddard                 !
!                                                                              !
!  Libxc functionals:                                                          !
!   LIBXC - specify actual functional id with the libxc_x/cfunc_id parameters  !
!                                                                              !
!  Non-functionals:                                                            !
!   HF     - only Hartree-Fock exchange (for solving Hartree-Fock equations)   !
!   NONE   - no exchange-correlation (for solving Hartree equations)           !
!                                                                              !
!   In addition, the following defaults apply:                                 !
!                                                                              !
!   LDA -> CAPZ                                                                !
!   GGA -> RPBE                                                                !
!                                                                              !
!==============================================================================!

module xc

  use constants, only: DP
#ifdef LIBXC
   use xc_f90_types_m
#endif

  implicit none

  private

  ! Public subroutines
  public :: xc_energy_potential
  public :: xc_hfxinit
  public :: xc_init
  public :: xc_exit
  public :: xc_radial
  !public :: xc_test

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

contains

  subroutine xc_energy_potential(density_fine,xc_energy,xc_pot_fine,grid, &
       nhat_dim,nhat_den_grad,tddftxc)

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
    use xc_funcs

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid    ! Grid definition
    real(kind=DP), intent(inout) :: density_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
    real(kind=DP), intent(out) :: xc_energy
    real(kind=DP), intent(out) :: xc_pot_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
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
    call timer_clock('xc_energy_potential',1)

    ! ddor: Set default for loc_tddftxc
    loc_tddftxc = .false.
    if (present(tddftxc)) loc_tddftxc = tddftxc

    ! Check input density
    if (pub_output_detail == VERBOSE) call internal_check_density

    ! Allocate workspace if required
#ifdef LIBXC
    if (pub_xc_gradient_corrected.or.use_libxc) call internal_allocate_workspace
#else
    if (pub_xc_gradient_corrected) call internal_allocate_workspace
#endif

    ! Calculate gradients of the density if needed
    ! ndmh: supplying compensation density and its gradient if necessary
    if (pub_xc_gradient_corrected) then
       if (present(nhat_den_grad)) then
          call xc_gradients(density_fine,density_grad,recip_work,grid, &
               nhat_dim,nhat_den_grad)
       else
          call xc_gradients(density_fine,density_grad,recip_work,grid,0)
       end if
    end if

    ! Zero potential
    xc_pot_fine = 0.0_DP

    ! Call appropriate routine according to functional type
    select case (functional)
    case ('CAPZ')
       !----------------------------------------------------------------------!
       ! Local density approximation - from Ceperley-Alder Monte Carlo data,  !
       !  and Gell-Mann-Brueckner expansion, parameterised by Perdew & Zunger !
       !----------------------------------------------------------------------!
       call xc_lda(density_fine,xc_energy, xc_pot_fine, grid, xc_capz_c_point, &
            xc_capz_c_point_sp)
    case ('VWN')
       !----------------------------------------------------------------------!
       ! Local density approximation - Vosko, Wilk & Nusair                   !
       !----------------------------------------------------------------------!
       call xc_lda(density_fine,xc_energy, xc_pot_fine, grid, xc_vwn_c_point, &
            xc_vwn_c_point_sp)
    case ('PW92')
       !----------------------------------------------------------------------!
       ! Local density approximation - Perdew & Wang 1992 functional          !
       !----------------------------------------------------------------------!
       call xc_lda(density_fine,xc_energy, xc_pot_fine, grid, xc_pw92_c_point, &
            xc_pw92_c_point_sp)
    case ('BLYP')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: Becke88 + LYP            !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy, xc_pot_fine, xc_b88_x_point,xc_b88_x_point_sp, &
            xc_lyp_c_point, xc_lyp_c_point_sp)
    case ('PBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1996 version of          !
       ! Perdew, Burke & Ernzerhof                                            !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy,xc_pot_fine,xc_pbe_x_point,&
            xc_pbe_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp)
    case ('PBEX')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1996 version of          !
       ! Perdew, Burke & Ernzerhof                                            !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy,xc_pot_fine,xc_pbe_x_point,&
            xc_pbe_x_point_sp,xc_none_c_point,xc_none_c_point_sp)
    case ('REVPBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1998 revised version of  !
       ! PBE due to Zhang & Yang                                              !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy,xc_pot_fine,xc_revpbe_x_point,&
            xc_revpbe_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp)
    case ('RPBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1999 revised version of  !
       ! PBE due to Hammer, Hansen & Norskov                                  !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy,xc_pot_fine,xc_rpbe_x_point,&
            xc_rpbe_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp)
    case ('PW91')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1991 version of          !
       ! Perdew and Wang                                                      !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy,xc_pot_fine,xc_pw91_x_point, &
            xc_pw91_x_point_sp,xc_pw91_c_point,xc_pw91_c_point_sp)
    case ('XLYP')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: X + LYP                  !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy, xc_pot_fine, xc_x_x_point, &
            xc_x_x_point_sp, xc_lyp_c_point, xc_lyp_c_point_sp)
    case ('WC')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 2006 functional of Wu &  !
       ! Cohen                                                                !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy, xc_pot_fine, xc_wc_x_point, &
            xc_wc_x_point_sp, xc_pbe_c_point, xc_pbe_c_point_sp)
    case ('B1LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional B1LYP:                                             !
       ! 1/4 HF_X + 3/4 B88_X - LDA_X) + LYP_C                                !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy, xc_pot_fine, xc_b1_x_point, &
            xc_b1_x_point_sp, xc_lyp_c_point, xc_lyp_c_point_sp)
    case ('B1PW91')
       !----------------------------------------------------------------------!
       ! Hybrid functional B1PW91:                                             !
       ! 1/4 HF_X + 3/4 B88_X + PW91_C                                        !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy, xc_pot_fine, xc_b1_x_point, &
            xc_b1_x_point_sp, xc_pw91_c_point, xc_pw91_c_point_sp)
    case ('B3LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional B3LYP:                                             !
       ! LDA + a(HF_X - LDA_X) + b (B88_X - LDA_X) + c(LYP_C - VWN_C)         !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy, xc_pot_fine, xc_b3_x_point, &
            xc_b3_x_point_sp, xc_b3lyp_c_point, xc_b3lyp_c_point_sp)
    case ('B3PW91')
       !----------------------------------------------------------------------!
       ! Hybrid functional B3PW91:                                            !
       ! LDA + a(HF_X - LDA_X) + b (B88_X - LDA_X) + c(PW91_C - VWN_C)        !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy, xc_pot_fine, xc_b3_x_point, &
            xc_b3_x_point_sp, xc_b3pw91_c_point, xc_b3pw91_c_point_sp)
    case ('PBE0')
       !----------------------------------------------------------------------!
       ! Hybrid functional PBE0: PBE + 1/4 (HF_X - PBE_X)                     !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy, xc_pot_fine, xc_pbe0_x_point,&
            xc_pbe0_x_point_sp, xc_pbe_c_point, xc_pbe_c_point_sp)
    case ('X3LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional X3LYP:                                             !
       ! LDA + a(HF_X - LDA_X) + b ((d B88_X + e PW91_X - LDA_X)              !
       !     + c(LYP_C - VWN_C)                                               !
       !----------------------------------------------------------------------!
       call xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
            xc_energy, xc_pot_fine, xc_x3_x_point, &
            xc_x3_x_point_sp, xc_x3lyp_c_point, xc_x3lyp_c_point_sp)
    case ('NONE','HF')
       !----------------------------------------------------------------------!
       ! No exchange-correlation functional (Hartree approximation) or only   !
       ! Hartree-Fock exchange (Hartree-Fock approximation).                  !
       !----------------------------------------------------------------------!
       xc_energy = 0.0_DP
       xc_pot_fine = 0.0_DP
    case ('LIBXC')
       !----------------------------------------------------------------------!
       ! Functional from LIBXC library - type specified by pub_libxc_func_id  !
       !----------------------------------------------------------------------!
#ifdef LIBXC
       if (pub_xc_gradient_corrected) then
          call xc_libxc(density_fine,recip_work,xc_energy,xc_pot_fine,grid, &
               density_grad,density_aux)
       else
          call xc_libxc(density_fine,recip_work,xc_energy,xc_pot_fine,grid)
       end if
#else
       call utils_abort('Error in xc_energy_potential: Functional "LIBXC" &
            &specified, but LIBXC not present in build')
#endif
    case default
       !----------------------------------------------------------------------!
       ! Unrecognised functional type                                         !
       !----------------------------------------------------------------------!
       if (pub_on_root) write(stdout,'(2a)') 'WARNING: unknown exchange-&
            &correlation functional: ', trim(functional)
       xc_energy = 0.0_DP
       xc_pot_fine = 0.0_DP
    end select

    ! Free workspace if required
    if (pub_xc_gradient_corrected) call internal_deallocate_workspace

    ! Stop timer
    call timer_clock('xc_energy_potential',2)

  contains

    subroutine internal_allocate_workspace

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
      call utils_alloc_check('internal_allocate_workspace (xc)', &
           'density_grad',ierr)

      ! Allocate workspace for the auxiliary density array
      allocate(density_aux(grid%ld1,grid%ld2,grid%max_slabs12, &
           pub_cell%num_spins),stat=ierr)
      call utils_alloc_check('internal_allocate_workspace (xc)', &
           'density_aux',ierr)

      ! Allocate reciprocal space workspace
      allocate(recip_work(grid%ld3,grid%ld2,grid%max_slabs23,3),stat=ierr)
      call utils_alloc_check('internal_allocate_workspace (xc)', &
           'recip_work',ierr)

    end subroutine internal_allocate_workspace

    subroutine internal_deallocate_workspace

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

    end subroutine internal_deallocate_workspace

    subroutine internal_check_density

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

         denmin = density_fine(1,1,1,is)

         a3: do islab12=1,grid%num_my_slabs12
            a2: do i2=1,grid%n2
               a1: do i1=1,grid%n1

                  denmin = min(denmin,density_fine(i1,i2,islab12,is))

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

    end subroutine internal_check_density

  end subroutine xc_energy_potential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_init(use_tddftxc)

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
          call xc_f90_func_init(libxc_func(1), libxc_info(1), &
               pub_libxc_x_func_id, libxc_polarized)
          libxc_family = xc_f90_info_family(libxc_info(1))
          libxc_gc = (libxc_family==XC_FAMILY_GGA) .or. &
               (libxc_family==XC_FAMILY_HYB_GGA)
       end if
       if (pub_libxc_c_func_id > 0) then
          call xc_f90_func_init(libxc_func(2), libxc_info(2), &
               pub_libxc_c_func_id, libxc_polarized)
          libxc_family = xc_f90_info_family(libxc_info(2))
          if (libxc_gc.neqv.((libxc_family==XC_FAMILY_GGA) .or. &
               (libxc_family==XC_FAMILY_HYB_GGA))) then
             call utils_abort('Error in xc_init: Mixing of LIBXC functionals &
                  &with and without gradient corrections is not supported')
          end if
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

  end subroutine xc_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_exit

    !=========================================================================!
    ! This subroutine deallocates the internal variables in the module.       !
    !-------------------------------------------------------------------------!
    ! Originally written by Nicholas Hine in October 2010.                    !
    !=========================================================================!

    use rundat, only: xc_functional, tddft_xc_functional, &
         pub_libxc_x_func_id, pub_libxc_c_func_id
    use simulation_cell, only: pub_cell
    use utils, only: utils_abort
#ifdef LIBXC
    use xc_f90_types_m
    use xc_f90_lib_m
#endif

    implicit none

#ifdef LIBXC
    ! ndmh: de-initialisations for LIBXC 'functional'
    if (use_libxc) then
       if (pub_libxc_x_func_id > 0) then
          call xc_f90_func_end(libxc_func(1))
       end if
       if (pub_libxc_c_func_id > 0) then
          call xc_f90_func_end(libxc_func(2))
       end if
    end if
#endif

  end subroutine xc_exit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine xc_hfxinit

    !===========================================================!
    ! This subroutine sets pub_usehfx and pub_hfxfraction.      !
    !-----------------------------------------------------------!
    ! Written by Quintin Hill in October 2008.                  !
    ! Quintin Hill added hybrid functionals on 11/03/2009.      !
    !===========================================================!

    use rundat, only: xc_functional, pub_usehfx, pub_tightbox_fft_fine, &
         pub_tightbox_fft_coarse

    implicit none

    select case (xc_functional)

    case ('HF')
       pub_usehfx = .true.
       pub_hfxfraction = 1.0_DP
    case ('B1LYP','B1PW91','PBE0')
       pub_usehfx = .true.
       pub_hfxfraction = 0.25_DP
    case ('B3LYP','B3PW91')
       pub_usehfx = .true.
       pub_hfxfraction = 0.2_DP
    case ('X3LYP')
       pub_usehfx = .true.
       pub_hfxfraction = 0.218_DP
    case default
       pub_usehfx = .false.
       pub_hfxfraction = 0.0_DP
    end select

    ! ars/jd: allow toggling pub_tightbox_fft_coarse ON, but not OFF
    pub_tightbox_fft_coarse = pub_usehfx .or. pub_tightbox_fft_coarse

  end subroutine xc_hfxinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if 0
  subroutine xc_test(grid)

    !=========================================================================!
    ! This subroutine tests routines associated with the module.              !
    !-------------------------------------------------------------------------!
    ! Originally written by Peter Haynes in February 2004.                    !
    ! Modified by Quintin Hill in November 2008 to use utils_(de)alloc_check. !
    ! Modified by Quintin Hill in March 2009 to use new subroutines and add   !
    ! test of BLYP.                                                           !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce, pub_on_root, pub_my_node_id
    use constants, only: PI, DP, stdout
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check
    use xc_funcs

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid

    ! Local variables
    integer :: ierr                                     ! Error flag
    integer :: i1,i2,i3                                 ! Grid loop counters
    integer :: islab12                                  ! Slab loop counter
    integer :: m1,m2,m3                                 ! Grid midpoints
    integer :: npts                                     ! Number of points
    logical :: grad_alloc                               ! Flag for workspace
    real(kind=DP), parameter :: threshold = 0.005_DP    ! Threshold
    real(kind=DP) :: xc_energy                          ! XC energy
    real(kind=DP) :: xc_potint                          ! XC potential integral
    real(kind=DP) :: denint                             ! Density integral
    real(kind=DP) :: normfac                            ! Normalising factor
    real(kind=DP) :: r1,r2,r3                           ! Distances
    real(kind=DP) :: r3sq,r23sq,rsq                     ! Squared distances
    real(kind=DP) :: gamma                              ! Parameter for test
    real(kind=DP) :: grad(3)                            ! Gradient
    real(kind=DP) :: mod_grad                           ! Modulus of gradient
    real(kind=DP) :: ana_grad                           ! Analytic gradient
    real(kind=DP) :: num_grad                           ! Numerical gradient
    real(kind=DP) :: err,errsq,rmserr                   ! Errors in gradient
    real(kind=DP), allocatable :: density_fine(:,:,:,:) ! Density array
    real(kind=DP), allocatable :: xc_pot_fine(:,:,:,:)  ! Potential array

    if (pub_on_root) then
       write(stdout,'(/a/)') '======================= &
            &Exchange-correlation module test ======================='
       write(stdout,'(a/)') '  N.B. Orthorhombic cells only!'
    end if

    ! Initialise workspace for gradients if necessary
    grad_alloc = allocated(density_grad)
    if (.not. grad_alloc) call xc_workspace(grid)

    ! Allocate workspace
    allocate(density_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('xc_test', 'density_fine', ierr)
    allocate(xc_pot_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('xc_test', 'xc_pot_fine', ierr)

    ! Calculate grid midpoints
    m1 = grid%n1 / 2
    m2 = grid%n2 / 2
    m3 = grid%n3 / 2

    ! Calculate optimum exponent
    gamma = 1.0_DP / sqrt(pi * min(grid%n1, &
         grid%n2, grid%n3))

    ! Create a Gaussian density to test the routines
    normfac = pub_cell%nat * 2 * (0.5_DP*gamma/pi)**1.5_DP
    denint = 0.0_DP
    do islab12=1,grid%num_my_slabs12
       i3=grid%first_slab12(pub_my_node_id) + islab12 - 1
       r3 = (i3-m3) * 0.5_DP * pub_cell%d3
       r3sq = r3*r3
       do i2=1,grid%n2
          r2 = (i2-m2) * 0.5_DP * pub_cell%d2
          r23sq = r2*r2 + r3sq
          do i1=1,grid%n1
             r1 = (i1-m1) * 0.5_DP * pub_cell%d1
             rsq = r1*r1 + r23sq
             density_fine(i1,i2,islab12,1) = normfac * &
                  exp(-0.5_DP * gamma * rsq)
             denint = denint + density_fine(i1,i2,i3,1)
          end do
       end do
    end do
    call comms_reduce('SUM',denint)
    denint = denint * grid%weight

    ! Calculate gradients of the density
    call xc_gradients(density_fine,grid)

    ! Compare gradients against analytic result
    npts = 0
    errsq = 0.0_DP
    do islab12=1,grid%num_my_slabs12
       i3=grid%first_slab12(pub_my_node_id) + islab12 - 1
       r3 = (i3-m3) * 0.5_DP * pub_cell%d3
       r3sq = r3*r3
       do i2=1,grid%n2
          r2 = (i2-m2) * 0.5_DP * pub_cell%d2
          r23sq = r2*r2 + r3sq
          do i1=1,grid%n1
             r1 = (i1-m1) * 0.5_DP * pub_cell%d1
             rsq = r1*r1 + r23sq
             if (gamma * rsq < threshold) then
                grad = density_grad(i1,i2,islab12,:,1)
                mod_grad = sqrt(grad(1)*grad(1) + grad(2)*grad(2) + &
                     grad(3)*grad(3))
                num_grad = mod_grad / density_fine(i1,i2,islab12,1)
                ana_grad = gamma * sqrt(rsq)
                err = num_grad - ana_grad
                errsq = errsq + err*err
                npts = npts + 1
             end if
          end do
       end do
    end do
    call comms_reduce('SUM',errsq)
    call comms_reduce('SUM',npts)
    rmserr = sqrt(errsq / npts)
    if (pub_on_root) then
       write(stdout,'(a)') '  Gradient calculation test:'
       write(stdout,'(a,f12.7)') '    Gaussian density variance  : ',&
            1.0_DP / gamma
       write(stdout,'(a,f12.7)') '    Gaussian density integral  : ',denint
       write(stdout,'(a,e12.4)') '    Error in density integral  : ', &
            abs(denint-2.0_DP*pub_cell%nat)
       write(stdout,'(a,e12.4)') '    RMS error in gradient      : ',rmserr
    end if

    ! Call appropriate routines for each functional type available:

    ! Ceperley-Alder LDA:
    call xc_lda(density_fine, xc_energy, xc_pot_fine, grid, xc_capz_c_point, &
         xc_capz_c_point_sp, xc_potint)
    if (pub_on_root) then
       write(stdout,'(/a)') '  Ceperley-Alder LDA:'
       write(stdout,'(a,f16.8)') '    Exchange-correlation energy: ',xc_energy
       write(stdout,'(a,f16.8)') '    Potential-density integral : ',xc_potint
    end if

    ! Vosko-Wilk-Nusair LDA:
    call xc_lda(density_fine, xc_energy, xc_pot_fine, grid, xc_vwn_c_point, &
         xc_vwn_c_point_sp, xc_potint)
    if (pub_on_root) then
       write(stdout,'(/a)') '  Vosko-Wilk-Nusair LDA:'
       write(stdout,'(a,f16.8)') '    Exchange-correlation energy: ',xc_energy
       write(stdout,'(a,f16.8)') '    Potential-density integral : ',xc_potint
    end if

    ! Perdew-Burke-Ernzerhof GGA:
    call xc_gc(density_fine, xc_energy, xc_pot_fine, grid, xc_pbe_x_point,&
         xc_pbe_x_point_sp, xc_pbe_c_point, xc_pbe_c_point_sp, xc_potint)
    if (pub_on_root) then
       write(stdout,'(/a)') '  Perdew-Burke-Ernzerhof GGA:'
       write(stdout,'(a,f16.8)') '    Exchange-correlation energy: ',xc_energy
       write(stdout,'(a,f16.8)') '    Potential-density integral : ',xc_potint
    end if

    ! Perdew-Wang 91 GGA:
    call xc_gc(density_fine, xc_energy, xc_pot_fine, grid, xc_pw91_x_point,&
         xc_pw91_x_point_sp, xc_pw91_c_point, xc_pw91_c_point_sp,xc_potint)
    if (pub_on_root) then
       write(stdout,'(/a)') '  Perdew-Wang ''91 GGA:'
       write(stdout,'(a,f16.8)') '    Exchange-correlation energy: ',xc_energy
       write(stdout,'(a,f16.8)') '    Potential-density integral : ',xc_potint
    end if

    ! BLYP GGA:
    call xc_gc(density_fine, xc_energy, xc_pot_fine, grid, xc_b88_x_point, &
         xc_b88_x_point_sp, xc_lyp_c_point, xc_lyp_c_point_sp,xc_potint)
    if (pub_on_root) then
       write(stdout,'(/a)') '  BLYP GGA:'
       write(stdout,'(a,f16.8)') '    Exchange-correlation energy: ',xc_energy
       write(stdout,'(a,f16.8)') '    Potential-density integral : ',xc_potint
    end if

    ! Deallocate workspace
    deallocate(xc_pot_fine,stat=ierr)
    call utils_dealloc_check('xc_test', 'xc_pot_fine', ierr)
    deallocate(density_fine,stat=ierr)
    call utils_alloc_check('xc_test', 'density_fine', ierr)

    ! Free workspace if required
    if (.not. grad_alloc) call xc_free

    if (pub_on_root) write(stdout,'(/a/)') '=======================&
         &========================================================='

  end subroutine xc_test
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_radial(npts,npts_max,ns,den_rad,grad_rad,vxc_rad,dvxc_dn_rad, &
       exc_rad)

    !==========================================================================!
    ! This subroutine calculates the exchange-correlation energy and potential !
    ! of a radial density distribution, for atoms and PAW calculations.        !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in June 2010.                                   !
    !==========================================================================!

    use comms, only: pub_on_root
    use constants, only: max_spins, stdout
    use utils, only: utils_abort
    use xc_funcs

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    integer, intent(in) :: npts_max
    integer, intent(in) :: ns
    real(kind=DP), intent(in) :: den_rad(npts_max,ns)
    real(kind=DP), intent(in) :: grad_rad(npts_max,ns)
    real(kind=DP), intent(out) :: vxc_rad(npts_max,ns)
    real(kind=DP), intent(out), optional :: dvxc_dn_rad(npts_max,ns)
    real(kind=DP), intent(out), optional :: exc_rad(npts_max,ns)

    ! Local Variables
    integer :: ipt
    real(kind=DP) :: x_energy, x_pot(max_spins)
    real(kind=DP) :: c_energy, c_pot(max_spins)

    select case (functional)
    case ('CAPZ')
       !----------------------------------------------------------------------!
       ! Local density approximation - from Ceperley-Alder Monte Carlo data,  !
       !  and Gell-Mann-Brueckner expansion, parameterised by Perdew & Zunger !
       !----------------------------------------------------------------------!
       call internal_xc_radial_lda(xc_capz_c_point,xc_capz_c_point_sp)
    case ('VWN')
       !----------------------------------------------------------------------!
       ! Local density approximation - Vosko, Wilk & Nusair                   !
       !----------------------------------------------------------------------!
       call internal_xc_radial_lda(xc_vwn_c_point,xc_vwn_c_point_sp)
    case ('PW92')
       !----------------------------------------------------------------------!
       ! Local density approximation - Perdew-Wang 1992 functional            !
       !----------------------------------------------------------------------!
       call internal_xc_radial_lda(xc_pw92_c_point,xc_pw92_c_point_sp)
    case ('BLYP')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: Becke88 + LYP            !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_b88_x_point,xc_b88_x_point_sp, &
            xc_lyp_c_point, xc_lyp_c_point_sp)
    case ('PBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1996 version of          !
       ! Perdew, Burke & Ernzerhof                                            !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_pbe_x_point,xc_pbe_x_point_sp, &
            xc_pbe_c_point,xc_pbe_c_point_sp)
    case ('REVPBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1998 revised version of  !
       ! PBE due to Zhang & Yang                                              !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_revpbe_x_point,xc_revpbe_x_point_sp, &
            xc_pbe_c_point,xc_pbe_c_point_sp)
    case ('RPBE')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1999 revised version of  !
       ! PBE due to Hammer, Hansen & Norskov                                  !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_rpbe_x_point,&
            xc_rpbe_x_point_sp,xc_pbe_c_point,xc_pbe_c_point_sp)
    case ('PW91')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 1991 version of          !
       ! Perdew and Wang                                                      !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_pw91_x_point, &
            xc_pw91_x_point_sp,xc_pw91_c_point,xc_pw91_c_point_sp)
    case ('XLYP')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: X + LYP                  !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_x_x_point, &
            xc_x_x_point_sp, xc_lyp_c_point, xc_lyp_c_point_sp)
    case ('WC')
       !----------------------------------------------------------------------!
       ! Generalised gradient density approximation: 2006 functional of Wu &  !
       ! Cohen                                                                !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_wc_x_point, &
            xc_wc_x_point_sp, xc_pbe_c_point, xc_pbe_c_point_sp)
    case ('B1LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional B1LYP:                                             !
       ! 1/4 HF_X + 3/4 B88_X - LDA_X) + LYP_C                                !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_b1_x_point, &
            xc_b1_x_point_sp, xc_lyp_c_point, xc_lyp_c_point_sp)
    case ('B1PW91')
       !----------------------------------------------------------------------!
       ! Hybrid functional B1PW91:                                             !
       ! 1/4 HF_X + 3/4 B88_X + PW91_C                                        !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_b1_x_point, &
            xc_b1_x_point_sp, xc_pw91_c_point, xc_pw91_c_point_sp)
    case ('B3LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional B3LYP:                                             !
       ! LDA + a(HF_X - LDA_X) + b (B88_X - LDA_X) + c(LYP_C - VWN_C)         !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_b3_x_point, &
            xc_b3_x_point_sp, xc_b3lyp_c_point, xc_b3lyp_c_point_sp)
    case ('B3PW91')
       !----------------------------------------------------------------------!
       ! Hybrid functional B3PW91:                                            !
       ! LDA + a(HF_X - LDA_X) + b (B88_X - LDA_X) + c(PW91_C - VWN_C)        !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_b3_x_point, &
            xc_b3_x_point_sp, xc_b3pw91_c_point, xc_b3pw91_c_point_sp)
    case ('PBE0')
       !----------------------------------------------------------------------!
       ! Hybrid functional PBE0: PBE + 1/4 (HF_X - PBE_X)                     !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_pbe0_x_point,&
            xc_pbe0_x_point_sp, xc_pbe_c_point, xc_pbe_c_point_sp)
    case ('X3LYP')
       !----------------------------------------------------------------------!
       ! Hybrid functional X3LYP:                                             !
       ! LDA + a(HF_X - LDA_X) + b ((d B88_X + e PW91_X - LDA_X)              !
       !     + c(LYP_C - VWN_C)                                               !
       !----------------------------------------------------------------------!
       call internal_xc_radial_gc(xc_x3_x_point, &
            xc_x3_x_point_sp, xc_x3lyp_c_point, xc_x3lyp_c_point_sp)
    case ('LIBXC')
       !----------------------------------------------------------------------!
       ! Functional from LIBXC library - type specified by pub_libxc_func_id  !
       !----------------------------------------------------------------------!
#ifdef LIBXC
       call internal_xc_radial_libxc
#else
       call utils_abort('Error in xc_energy_potential: Functional "LIBXC" &
            &specified, but LIBXC not present in build')
#endif
    case ('NONE','HF')
       !----------------------------------------------------------------------!
       ! No exchange-correlation functional (Hartree approximation) or only   !
       ! Hartree-Fock exchange (Hartree-Fock approximation).                  !
       !----------------------------------------------------------------------!
       vxc_rad = 0.0_DP
       if (present(exc_rad)) exc_rad = 0.0_DP
    case default
       !----------------------------------------------------------------------!
       ! Unrecognised functional type                                         !
       !----------------------------------------------------------------------!
       if (pub_on_root) write(stdout,'(2a)') 'WARNING: unknown exchange-&
            &correlation functional: ', trim(functional)
    end select

contains

    subroutine internal_xc_radial_lda(c_functional,c_functional_sp)

      interface
         subroutine c_functional(den,c_energy,c_pot)
           use constants, only: DP, PI
           real(kind=DP), intent(in)  :: den
           real(kind=DP), intent(out) :: c_energy
           real(kind=DP), intent(out) :: c_pot
         end subroutine c_functional
      end interface

      interface
         subroutine c_functional_sp(den,den1,den2,c_energy,c_pot1,c_pot2)
           use constants, only: DP, PI
           real(kind=DP), intent(in)  :: den
           real(kind=DP), intent(in)  :: den1
           real(kind=DP), intent(in)  :: den2
           real(kind=DP), intent(out) :: c_energy
           real(kind=DP), intent(out) :: c_pot1
           real(kind=DP), intent(out) :: c_pot2
         end subroutine c_functional_sp
      end interface

      if (ns==1) then
         do ipt=1,npts
            call xc_lda_x_point(den_rad(ipt,1),x_energy,x_pot(1))
            call c_functional(den_rad(ipt,1),c_energy,c_pot(1))
            vxc_rad(ipt,1) = x_pot(1) + c_pot(1)
            if (present(exc_rad)) exc_rad(ipt,1) = x_energy + c_energy
         end do
      else
         do ipt=1,npts
            call xc_lsda_x_point(den_rad(ipt,1),den_rad(ipt,2),x_energy, &
                 x_pot(1),x_pot(2))
            call c_functional_sp(sum(den_rad(ipt,:)),den_rad(ipt,1), &
                 den_rad(ipt,2),c_energy,c_pot(1),c_pot(2))
            vxc_rad(ipt,:) = x_pot(:) + c_pot(:)
            if (present(exc_rad)) exc_rad(ipt,1) = x_energy + c_energy
         end do
      end if

    end subroutine internal_xc_radial_lda

    subroutine internal_xc_radial_gc(x_functional,x_functional_sp, &
         c_functional,c_functional_sp)

      interface
         subroutine x_functional(den,mgd,x_energy,x_pot,x_dfdmgd)
           use constants, only: DP, PI
           real(kind=DP), intent(in)  :: den
           real(kind=DP), intent(in)  :: mgd
           real(kind=DP), intent(out) :: x_energy
           real(kind=DP), intent(out) :: x_pot
           real(kind=DP), intent(out) :: x_dfdmgd
         end subroutine x_functional
      end interface

      interface
         subroutine x_functional_sp(den1,den2,mgd1,mgd2,x_energy,x_pot1,x_pot2,&
              x_dfdmgd1,x_dfdmgd2)
           use constants, only: DP, PI
           real(kind=DP), intent(in)  :: den1
           real(kind=DP), intent(in)  :: den2
           real(kind=DP), intent(in)  :: mgd1
           real(kind=DP), intent(in)  :: mgd2
           real(kind=DP), intent(out) :: x_energy
           real(kind=DP), intent(out) :: x_pot1
           real(kind=DP), intent(out) :: x_pot2
           real(kind=DP), intent(out) :: x_dfdmgd1
           real(kind=DP), intent(out) :: x_dfdmgd2
         end subroutine x_functional_sp
      end interface

      interface
         subroutine c_functional(den,mgd,c_energy,c_pot,c_dfdmgd)
           use constants, only: DP, PI
           real(kind=DP), intent(in)  :: den
           real(kind=DP), intent(in)  :: mgd
           real(kind=DP), intent(out) :: c_energy
           real(kind=DP), intent(out) :: c_pot
           real(kind=DP), intent(out) :: c_dfdmgd
         end subroutine c_functional
      end interface

      interface
         subroutine c_functional_sp(den,den1,den2,mgd,mgd1,mgd2,c_energy,&
              c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)
           use constants, only: DP, PI
           real(kind=DP), intent(in)  :: den
           real(kind=DP), intent(in)  :: den1
           real(kind=DP), intent(in)  :: den2
           real(kind=DP), intent(in)  :: mgd
           real(kind=DP), intent(in)  :: mgd1
           real(kind=DP), intent(in)  :: mgd2
           real(kind=DP), intent(out) :: c_energy
           real(kind=DP), intent(out) :: c_pot1
           real(kind=DP), intent(out) :: c_pot2
           real(kind=DP), intent(out) :: c_dfdmgd
           real(kind=DP), intent(out) :: c_dfdmgd1
           real(kind=DP), intent(out) :: c_dfdmgd2
           real(kind=DP), intent(out) :: c_dfdmgd12
         end subroutine c_functional_sp
      end interface

      ! Local Variables
      real(kind=DP) :: mgd1, mgd2, mgd
      real(kind=DP) :: x_dfdmgd(2)
      real(kind=DP) :: c_dfdmgd(0:3)

      if (ns==1) then
         do ipt=1,npts
            mgd = abs(grad_rad(ipt,1))
            call x_functional(den_rad(ipt,1),mgd, x_energy, &
                 x_pot(1), x_dfdmgd(1))
            call c_functional(den_rad(ipt,1),mgd, c_energy, &
                 c_pot(1), c_dfdmgd(1))
            vxc_rad(ipt,1) = x_pot(1) + c_pot(1)
            if (present(dvxc_dn_rad)) dvxc_dn_rad(ipt,1) = &
                 x_dfdmgd(1) + c_dfdmgd(1)
            if (present(exc_rad)) exc_rad(ipt,1) = x_energy + c_energy
         end do
      else
         do ipt=1,npts
            mgd1 = abs(grad_rad(ipt,1))
            mgd1 = abs(grad_rad(ipt,1))
            mgd = abs(grad_rad(ipt,1)+grad_rad(ipt,2))
            call x_functional_sp(den_rad(ipt,1),den_rad(ipt,2),mgd1, &
                 mgd2,x_energy,x_pot(1),x_pot(2),x_dfdmgd(1), &
                 x_dfdmgd(2))
            call c_functional_sp(sum(den_rad(ipt,:)),den_rad(ipt,1), &
                 den_rad(ipt,2),mgd,mgd1,mgd2,c_energy,c_pot(1),c_pot(2), &
                 c_dfdmgd(0),c_dfdmgd(1),c_dfdmgd(2),c_dfdmgd(3))
            vxc_rad(ipt,:) = x_pot(:) + c_pot(:)
            if (present(dvxc_dn_rad)) then
               dvxc_dn_rad(ipt,1) = x_dfdmgd(1) + c_dfdmgd(0) + c_dfdmgd(1)
               dvxc_dn_rad(ipt,2) = x_dfdmgd(2) + c_dfdmgd(0) + c_dfdmgd(2)
               ! TODO: deal with c_dfdmgd(3)
               if (c_dfdmgd(3)>0.0_DP) call utils_abort('Error in xc_radial: &
                    &nonzero c_dfdmgd(3): Functional not yet supported by &
                    &xc_radial')
            end if
            if (present(exc_rad)) exc_rad(ipt,1) = x_energy + c_energy
         end do
      end if

    end subroutine internal_xc_radial_gc

#ifdef LIBXC
    subroutine internal_xc_radial_libxc

      use rundat, only: pub_libxc_x_func_id, pub_libxc_c_func_id

      ! Local Variables
      real(kind=xc_f90_kind) :: den(max_spins) ! Density
      real(kind=xc_f90_kind) :: grad(3,max_spins)  ! Density gradient
      real(kind=xc_f90_kind) :: sigma(3) ! Modulus of density gradient
      real(kind=xc_f90_kind) :: ex                  ! Accumulated xc energy
      real(kind=xc_f90_kind) :: ec                  ! Accumulated xc energy
      real(kind=xc_f90_kind) :: lx_pot(max_spins)   ! xc potential at point
      real(kind=xc_f90_kind) :: lc_pot(max_spins)   ! xc potential at point
      real(kind=xc_f90_kind) :: x_dfdsigma(3)       ! Derivative df_{xc}/dsigma
      real(kind=xc_f90_kind) :: c_dfdsigma(3)       ! Derivative df_{xc}/dsigma
      real(kind=xc_f90_kind), parameter :: two = real(2,kind=xc_f90_kind)
      integer :: ierr

      do ipt=1,npts

         den(1:ns) = den_rad(ipt,1:ns)

         ! ndmh: Calculate energy and potential at this point
         if (.not.libxc_gc) then  ! Non-Gradient Corrected

            ! ndmh: Exchange part
            if (pub_libxc_x_func_id > 0) then
               call xc_f90_lda_exc_vxc(libxc_func(1), 1, den(1), &
                    ex, lx_pot(1))
            else
               ex = 0.0_DP; lx_pot(:) = 0.0_DP
            end if

            ! ndmh: Correlation part
            if (pub_libxc_c_func_id > 0) then
               call xc_f90_lda_exc_vxc(libxc_func(2), 1, den(1), &
                    ec, lc_pot(1))
            else
               ec = 0.0_DP; lc_pot(:) = 0.0_DP
            end if

         else  ! Gradient Corrected

            ! Calculate reduced gradients
            if (ns==2) then
               sigma(1) = abs(grad_rad(ipt,1))**2
               sigma(2) = grad_rad(ipt,1)*grad_rad(ipt,2)
               sigma(3) = abs(grad_rad(ipt,2))**2
            else
               sigma(1) = abs(grad_rad(ipt,1))**2
            end if

            ! ndmh: Exchange part
            if (pub_libxc_x_func_id > 0) then
               call xc_f90_gga_exc_vxc(libxc_func(1), 1, den(1), sigma(1), &
                    ex, x_pot(1), x_dfdsigma(1))
            else
               ex = 0.0_DP; x_pot(:) = 0.0_DP; x_dfdsigma(:) = 0.0_DP
            end if

            ! ndmh: Correlation part
            if (pub_libxc_c_func_id > 0) then
               call xc_f90_gga_exc_vxc(libxc_func(2), 1, den(1), sigma(1), &
                    ec, c_pot(1), c_dfdsigma(1))
            else
               ec = 0.0_DP; c_pot(:) = 0.0_DP; c_dfdsigma(:) = 0.0_DP
            end if

            if (present(dvxc_dn_rad)) then
               dvxc_dn_rad(ipt,1) = (x_dfdsigma(1) + c_dfdsigma(1)) &
                    * abs(grad_rad(ipt,1)) / 2.0_DP
               if (ns==2) dvxc_dn_rad(ipt,2) = (x_dfdsigma(3) + c_dfdsigma(3)) &
                    * abs(grad_rad(ipt,2)) / 2.0_DP
               ! TODO: deal with c_dfdsigma(2)
               if (c_dfdsigma(2)>0.0_DP) call utils_abort('Error in xc_radial: &
                    &nonzero c_dfdsigma(2): Functional not yet supported by &
                    xc_radial')
            end if

         end if

      end do

    end subroutine internal_xc_radial_libxc
#endif

  end subroutine xc_radial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_lda(density_fine,xc_energy,xc_pot_fine,grid,c_functional,&
       c_functional_sp,xc_potint)

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
    use xc_funcs, only: xc_lda_x_point, xc_lsda_x_point

    implicit none

    ! Arguments
    ! Grid definition
    type(GRID_INFO), intent(in) :: grid

    ! Input density on the fine grid
    real(kind=DP), intent(in) :: density_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)

    ! Exchange-correlation energy
    real(kind=DP), intent(out) :: xc_energy

    ! Exchange-correlation potential
    real(kind=DP), intent(out) :: xc_pot_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)

    ! Exchange-correlation potential-density integral
    real(kind=DP), optional, intent(out) :: xc_potint

    interface
       subroutine c_functional(den,c_energy,c_pot)
         use constants, only: DP, PI
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(out) :: c_energy
         real(kind=DP), intent(out) :: c_pot
       end subroutine c_functional
    end interface

    interface
       subroutine c_functional_sp(den,den1,den2,c_energy,c_pot1,c_pot2)
         use constants, only: DP, PI
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(in)  :: den1
         real(kind=DP), intent(in)  :: den2
         real(kind=DP), intent(out) :: c_energy
         real(kind=DP), intent(out) :: c_pot1
         real(kind=DP), intent(out) :: c_pot2
       end subroutine c_functional_sp
    end interface

    ! Local variables:
    integer :: i1,i2,islab12   ! Fine grid loop counters
    real(kind=DP) :: ex        ! Exchange energy in Ha
    real(kind=DP) :: ec        ! Correlation energy in Ha
    real(kind=DP) :: potint    ! Potential-density integral
    real(kind=DP) :: den       ! Charge density
    real(kind=DP) :: den1      ! Charge densities for spin 1
    real(kind=DP) :: den2      ! Charge densities for spin 2
    real(kind=DP) :: x_pot     ! Exchange potential
    real(kind=DP) :: x_pot1    ! Exchange potential spin 1
    real(kind=DP) :: x_pot2    ! Exchange potential spin2
    real(kind=DP) :: x_energy  ! Exchange energy at point
    real(kind=DP) :: c_energy  ! Correlation energy at point
    real(kind=DP) :: c_pot     ! Correlation potential
    real(kind=DP) :: c_pot1    ! Correlation potential spin 1
    real(kind=DP) :: c_pot2    ! Correlation potential spin 2

    ex = 0.0_DP ; ec = 0.0_DP
    potint = 0.0_DP

    xc_pot_fine = 0.0_DP ! jd: To clear margins between n1 and ld1 etc.

    if (pub_cell%num_spins == 1) then

       ! qoh: Spin-unpolarised case:
       a3: do islab12=1,grid%num_my_slabs12
          a2: do i2=1,grid%n2
             a1: do i1=1,grid%n1

                den = density_fine(i1,i2,islab12,1)

                call xc_lda_x_point(den,x_energy,x_pot)
                call c_functional(den,c_energy,c_pot)

                ex = ex + x_energy
                ec = ec + c_energy

                xc_pot_fine(i1,i2,islab12,1) = x_pot + c_pot
                potint = potint + den * xc_pot_fine(i1,i2,islab12,1)

             end do a1
          end do a2
       end do a3

    else

       ! qoh: Spin-polarised case:
       a3_sp: do islab12=1,grid%num_my_slabs12
          a2_sp: do i2=1,grid%n2
             a1_sp: do i1=1,grid%n1

                den1 = density_fine(i1,i2,islab12,1)
                den2 = density_fine(i1,i2,islab12,2)
                den = den1 + den2

                call xc_lsda_x_point(den1,den2,x_energy,x_pot1,x_pot2)
                call c_functional_sp(den,den1,den2,c_energy,c_pot1,c_pot2)

                ex = ex + x_energy
                ec = ec + c_energy

                xc_pot_fine(i1,i2,islab12,1) = x_pot1 + c_pot1
                xc_pot_fine(i1,i2,islab12,2) = x_pot2 + c_pot2
                potint = potint + den1 * xc_pot_fine(i1,i2,islab12,1) + &
                     den2 * xc_pot_fine(i1,i2,islab12,2)

             end do a1_sp
          end do a2_sp
       end do a3_sp

    end if

    ! Global reductions
    call comms_reduce('SUM',ex)
    call comms_reduce('SUM',ec)
    call comms_reduce('SUM',potint)

    xc_energy = (ex + ec) * grid%weight

    if (present(xc_potint)) xc_potint = potint * grid%weight

  end subroutine xc_lda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_gc(density_fine,density_grad,density_aux,recip_work,grid, &
       xc_energy,xc_pot_fine,x_functional,x_functional_sp, &
       c_functional,c_functional_sp,xc_potint)

    !==========================================================================!
    ! This subroutine calculates the exchange-correlation energy and potential !
    ! for a gradient corrected functional (including hybrids), given the       !
    ! electronic density on the fine grid.  The exchange and correlation       !
    ! functionals to be used are given as arguments.                           !
    !--------------------------------------------------------------------------!
    ! Arguments                                                                !
    !   grid (in) : Grid definition                                            !
    !   density_fine (in) : Input density on the fine grid                     !
    !   xc_energy (out) : Exchange-correlation energy                          !
    !   xc_pot_fine (out) : Exchange-correlation potential                     !
    !   xc_potint (out) : Exchange-correlation potential-density integral      !
    !   x_functional,x_functional_sp (in) : X subroutines given as arguements  !
    !   c_functional,c_functional_sp (in) : C subroutines given as arguements  !
    !   density_grad (inout) : Density gradient on the fine grid               !
    !   density_aux (inout) : Auxiliary density on the fine grid               !
    !   recip_work (inout) : Fine Grid Recip-space workspace                   !
    !--------------------------------------------------------------------------!
    ! Written by Quintin Hill in February / March 2009 based on existing       !
    ! xc_pbe subroutine written by Peter Haynes in February 2004.              !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP, PI
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in) :: density_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
    real(kind=DP), intent(out) :: xc_energy
    real(kind=DP), intent(out) :: xc_pot_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
    real(kind=DP), optional, intent(out) :: xc_potint
    real(kind=DP),intent(inout) :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_cell%num_spins)
    real(kind=DP),intent(inout) :: density_aux(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
    complex(kind=DP),intent(inout) :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)

    interface
       subroutine x_functional(den,mgd,x_energy,x_pot,x_dfdmgd)
         use constants, only: DP, PI
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(in)  :: mgd
         real(kind=DP), intent(out) :: x_energy
         real(kind=DP), intent(out) :: x_pot
         real(kind=DP), intent(out) :: x_dfdmgd
       end subroutine x_functional
    end interface

    interface
       subroutine x_functional_sp(den1,den2,mgd1,mgd2,x_energy,x_pot1,x_pot2,&
            x_dfdmgd1,x_dfdmgd2)
         use constants, only: DP, PI
         real(kind=DP), intent(in)  :: den1
         real(kind=DP), intent(in)  :: den2
         real(kind=DP), intent(in)  :: mgd1
         real(kind=DP), intent(in)  :: mgd2
         real(kind=DP), intent(out) :: x_energy
         real(kind=DP), intent(out) :: x_pot1
         real(kind=DP), intent(out) :: x_pot2
         real(kind=DP), intent(out) :: x_dfdmgd1
         real(kind=DP), intent(out) :: x_dfdmgd2
       end subroutine x_functional_sp
    end interface

    interface
       subroutine c_functional(den,mgd,c_energy,c_pot,c_dfdmgd)
         use constants, only: DP, PI
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(in)  :: mgd
         real(kind=DP), intent(out) :: c_energy
         real(kind=DP), intent(out) :: c_pot
         real(kind=DP), intent(out) :: c_dfdmgd
       end subroutine c_functional
    end interface

    interface
       subroutine c_functional_sp(den,den1,den2,mgd,mgd1,mgd2,c_energy,&
            c_pot1,c_pot2,c_dfdmgd,c_dfdmgd1,c_dfdmgd2,c_dfdmgd12)
         use constants, only: DP, PI
         real(kind=DP), intent(in)  :: den
         real(kind=DP), intent(in)  :: den1
         real(kind=DP), intent(in)  :: den2
         real(kind=DP), intent(in)  :: mgd
         real(kind=DP), intent(in)  :: mgd1
         real(kind=DP), intent(in)  :: mgd2
         real(kind=DP), intent(out) :: c_energy
         real(kind=DP), intent(out) :: c_pot1
         real(kind=DP), intent(out) :: c_pot2
         real(kind=DP), intent(out) :: c_dfdmgd
         real(kind=DP), intent(out) :: c_dfdmgd1
         real(kind=DP), intent(out) :: c_dfdmgd2
         real(kind=DP), intent(out) :: c_dfdmgd12
       end subroutine c_functional_sp
    end interface

    real(kind=DP) :: den      ! Density
    real(kind=DP) :: den1     ! Density spin 1
    real(kind=DP) :: den2     ! Density spin 2
    real(kind=DP) :: mgd      ! Modulus of density gradient
    real(kind=DP) :: mgd1     ! Modulus of density gradient for spin 1
    real(kind=DP) :: mgd2     ! Modulus of density gradient for spin 2
    real(kind=DP) :: grad(3)  ! Density gradient
    real(kind=DP) :: grad1(3) ! Density gradient for spin 1
    real(kind=DP) :: grad2(3) ! Density gradient for spin 2
    real(kind=DP) :: potint   ! Potential integral
    !qoh: Exchange variables
    real(kind=DP) :: ex        ! Accumulated exchange energy
    real(kind=DP) :: x_energy  ! Exchange energy at point
    real(kind=DP) :: x_pot     ! Exchange potential at point
    real(kind=DP) :: x_pot1    ! Exchange potential at point for spin 1
    real(kind=DP) :: x_pot2    ! Exchange potential at point for spin 2
    real(kind=DP) :: x_dfdmgd! Derivative df_{x}/d|grad n|
    real(kind=DP) :: x_dfdmgd1! Derivative df_{x}/d|grad n| spin 1
    real(kind=DP) :: x_dfdmgd2! Derivative df_{x}/d|grad n| spin 2
    !qoh: Correlation variables
    real(kind=DP) :: ec        ! Accumulated correlation energy
    real(kind=DP) :: c_energy  ! Correlation energy at point
    real(kind=DP) :: c_pot     ! Correlation potential at point
    real(kind=DP) :: c_pot1    ! Correlation potential at point for spin 1
    real(kind=DP) :: c_pot2    ! Correlation potential at point for spin 2
    real(kind=DP) :: c_dfdmgd  ! Derivative df_{c}/d|grad n|
    real(kind=DP) :: c_dfdmgd1 ! Derivative df_{c}/d|grad n| spin 1
    real(kind=DP) :: c_dfdmgd2 ! Derivative df_{c}/d|grad n| spin 2
    real(kind=DP) :: c_dfdmgd12 ! df_{c}/d(grad n_1.grad n_2) component
    real(kind=DP) :: gd1_dot_gd2 ! (grad n_1).(grad n_2)
    integer :: islab12 ! Loop counter
    integer :: i1 ! Loop counter
    integer :: i2 ! Loop counter
    integer :: is ! Spin counter

    ex = 0.0_DP
    ec = 0.0_DP

    xc_pot_fine = 0.0_DP ! jd: To clear margins between n1 and ld1 etc.

    spin: if ( pub_cell%num_spins == 1 ) then
       !qoh: Spin unpolarised case
       a3: do islab12=1,grid%num_my_slabs12
          a2: do i2=1,grid%n2
             a1: do i1=1,grid%n1

                den = density_fine(i1,i2,islab12,1)
                grad = density_grad(i1,i2,islab12,:,1)
                mgd = sqrt(grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3))

                ! qoh: Calculate energy and potential at this point
                call x_functional(den,mgd,x_energy, x_pot, x_dfdmgd)
                call c_functional(den,mgd,c_energy, c_pot, c_dfdmgd)

                ! qoh: Store this point's energy and potential
                ex = x_energy + ex
                ec = c_energy + ec
                xc_pot_fine(i1,i2,islab12,1) = x_pot + c_pot

                ! qoh: Make sure we don't divide by zero
                if (mgd > 0.0_DP) then
                   density_grad(i1,i2,islab12,:,1) = &
                        -(c_dfdmgd + x_dfdmgd) * grad(:) / mgd
                else
                   density_grad(i1,i2,islab12,:,1) = 0.0_DP
                end if

             end do a1
          end do a2
       end do a3

    else
       !qoh: Spin polarised case
       a3_sp: do islab12=1,grid%num_my_slabs12
          a2_sp: do i2=1,grid%n2
             a1_sp: do i1=1,grid%n1

                den1 = density_fine(i1,i2,islab12,1)
                den2 = density_fine(i1,i2,islab12,2)
                den = den1 + den2

                grad1 = density_grad(i1,i2,islab12,:,1)
                grad2 = density_grad(i1,i2,islab12,:,2)
                grad = grad1 + grad2

                mgd1 = sqrt(grad1(1)*grad1(1) + grad1(2)*grad1(2) + &
                     grad1(3)*grad1(3))
                mgd2 = sqrt(grad2(1)*grad2(1) + grad2(2)*grad2(2) + &
                     grad2(3)*grad2(3))
                mgd = sqrt(grad(1)*grad(1)+grad(2)*grad(2)+grad(3)*grad(3))
                gd1_dot_gd2 = grad1(1)*grad2(1) + grad1(2)*grad2(2) + &
                     grad1(3)*grad2(3)
                ! qoh: Calculate energy and potential at this point
                call x_functional_sp(den1,den2,mgd1,mgd2,&
                     x_energy, x_pot1, x_pot2, x_dfdmgd1, x_dfdmgd2)
                call c_functional_sp(den,den1,den2,mgd,mgd1,mgd2,&
                     c_energy, c_pot1, c_pot2, &
                     c_dfdmgd, c_dfdmgd1, c_dfdmgd2, c_dfdmgd12)

                ! qoh: Store this point's energy and potential
                ex = x_energy + ex
                ec = c_energy + ec
                xc_pot_fine(i1,i2,islab12,1) = x_pot1 + c_pot1
                xc_pot_fine(i1,i2,islab12,2) = x_pot2 + c_pot2

                ! qoh: Make sure we don't divide by zero
                if (mgd1 > 0.0_DP) then
                   density_grad(i1,i2,islab12,:,1) = &
                        -(x_dfdmgd1+c_dfdmgd1) * grad1(:) / mgd1
                else
                   density_grad(i1,i2,islab12,:,1) = 0.0_DP
                end if
                if (mgd2 > 0.0_DP) then
                   density_grad(i1,i2,islab12,:,2) = &
                        -(x_dfdmgd2+c_dfdmgd2) * grad2(:) / mgd2
                else
                   density_grad(i1,i2,islab12,:,2) = 0.0_DP
                end if
                if (abs(gd1_dot_gd2)>0.0_DP) then
                   density_grad(i1,i2,islab12,:,1) = &
                        density_grad(i1,i2,islab12,:,1) &
                        -c_dfdmgd12 * grad2(:) / gd1_dot_gd2
                   density_grad(i1,i2,islab12,:,2) = &
                        density_grad(i1,i2,islab12,:,2) &
                        -c_dfdmgd12 * grad1(:) / gd1_dot_gd2
                end if
                if (mgd > 0.0_DP) then
                   density_grad(i1,i2,islab12,:,1) = &
                        density_grad(i1,i2,islab12,:,1) &
                        -c_dfdmgd * grad(:) / mgd
                   density_grad(i1,i2,islab12,:,2) = &
                        density_grad(i1,i2,islab12,:,2) &
                        -c_dfdmgd * grad(:) / mgd
                end if

                !if (den > 0.01_DP) then
                !   write(26,'(3i6,2f20.12,6f14.10,19f20.12)') &
                !        i1,i2,islab12,den1,den2, &
                !        grad1(1:3),grad2(1:3),x_energy,x_pot1,x_pot2, &
                !        c_energy,c_pot1,c_pot2,x_dfdmgd1,0.0_DP,x_dfdmgd2, &
                !        c_dfdmgd,c_dfdmgd1,c_dfdmgd12,c_dfdmgd2,&
                !        density_grad(i1,i2,islab12,1:3,1:2)
                !end if

             end do a1_sp
          end do a2_sp
       end do a3_sp
    end if spin

    ! Calculate divergence of gradient dependent part and add to potential
    call xc_divergence(density_grad,density_aux,recip_work,grid)
    xc_pot_fine = xc_pot_fine + density_aux

    ! Global reductions
    call comms_reduce('SUM',ex)
    call comms_reduce('SUM',ec)
    xc_energy = (ex+ec) * grid%weight

    ! qoh: Calculate potential if necessary

    if (present(xc_potint)) then
       potint = 0.0_DP
       do is=1,pub_cell%num_spins
          do islab12=1,grid%num_my_slabs12
             do i2=1,grid%n2
                do i1=1,grid%n1
                   potint = potint + xc_pot_fine(i1,i2,islab12,is) * &
                        density_fine(i1,i2,islab12,is)
                end do
             end do
          end do
       end do
       call comms_reduce('SUM',potint)
       xc_potint = potint * grid%weight
    end if

  end subroutine xc_gc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef LIBXC
  subroutine xc_libxc(density_fine,recip_work,xc_energy,xc_pot_fine,grid, &
       density_grad,density_aux,xc_potint)

    !==========================================================================!
    ! This subroutine calculates the exchange-correlation energy and potential !
    ! for a gradient corrected functional (including hybrids), given the       !
    ! electronic density on the fine grid.  The exchange and correlation       !
    ! functionals to be used are given as arguments.                           !
    !--------------------------------------------------------------------------!
    ! Arguments                                                                !
    !   grid (in) : Grid definition                                            !
    !   density_fine (in) : Input density on the fine grid                     !
    !   xc_energy (out) : Exchange-correlation energy                          !
    !   xc_pot_fine (out) : Exchange-correlation potential                     !
    !   xc_potint (out) : Exchange-correlation potential-density integral      !
    !   density_grad (inout) : Density gradient on the fine grid               !
    !   density_aux (inout) : Auxiliary density on the fine grid               !
    !   recip_work (inout) : Fine Grid Recip-space workspace                   !
    !--------------------------------------------------------------------------!
    ! Written by Nicholas Hine in July 2010, based partly on xc_gc by Quintin  !
    ! Hill.                                                                    !
    !==========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: comms_reduce
    use constants, only: DP, PI, max_spins
    use rundat, only: pub_libxc_x_func_id, pub_libxc_c_func_id
    use simulation_cell, only: pub_cell
#ifdef LIBXC
    use xc_f90_types_m
    use xc_f90_lib_m
#endif

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(in) :: density_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
    complex(kind=DP),intent(inout) :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)
    real(kind=DP), intent(out) :: xc_energy
    real(kind=DP), intent(out) :: xc_pot_fine(grid%ld1, grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
    real(kind=DP), optional, intent(inout) :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_cell%num_spins)
    real(kind=DP), optional, intent(inout) :: density_aux(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
    real(kind=DP), optional, intent(out) :: xc_potint

    ! Local Variables
    real(kind=DP) :: potint   ! Potential integral
    real(kind=xc_f90_kind) :: den(max_spins) ! Density
    real(kind=xc_f90_kind) :: grad(3,max_spins)  ! Density gradient
    real(kind=xc_f90_kind) :: sigma(3) ! Modulus of density gradient
    real(kind=xc_f90_kind) :: ex                  ! Accumulated xc energy
    real(kind=xc_f90_kind) :: ec                  ! Accumulated xc energy
    real(kind=xc_f90_kind) :: x_pot(max_spins)    ! xc potential at point
    real(kind=xc_f90_kind) :: c_pot(max_spins)    ! xc potential at point
    real(kind=xc_f90_kind) :: x_dfdsigma(3)       ! Derivative df_{xc}/dsigma
    real(kind=xc_f90_kind) :: c_dfdsigma(3)       ! Derivative df_{xc}/dsigma
    real(kind=xc_f90_kind), parameter :: two = real(2,kind=xc_f90_kind)
    integer :: islab12 ! Loop counter
    integer :: i1 ! Loop counter
    integer :: i2 ! Loop counter
    integer :: is ! Spin counter
    integer :: ns ! Spin counter

    ns = pub_cell%num_spins

    xc_energy = 0.0_DP
    xc_pot_fine = 0.0_DP ! jd: To clear margins between n1 and ld1 etc.

    a3: do islab12=1,grid%num_my_slabs12
       a2: do i2=1,grid%n2
          a1: do i1=1,grid%n1

             den(1:ns) = density_fine(i1,i2,islab12,1:ns)
             do is=1,ns
                if (den(is)<0.0_DP) den(is) = 0.0_DP
             end do

             ! ndmh: Calculate energy and potential at this point
             if (.not.libxc_gc) then  ! Non-Gradient Corrected

                ! ndmh: Exchange part
                if (pub_libxc_x_func_id > 0) then
                   call xc_f90_lda_exc_vxc(libxc_func(1), 1, den(1), &
                        ex, x_pot(1))
                else
                   ex = 0.0_DP; x_pot(:) = 0.0_DP
                end if

                ! ndmh: Correlation part
                if (pub_libxc_c_func_id > 0) then
                   call xc_f90_lda_exc_vxc(libxc_func(2), 1, den(1), &
                        ec, c_pot(1))

                else
                   ec = 0.0_DP; c_pot(:) = 0.0_DP
                end if

             else  ! Gradient Corrected

                grad(1:3,1:ns) = density_grad(i1,i2,islab12,1:3,1:ns)

                ! Calculate reduced gradients
                if (ns==2) then
                   sigma(1) = grad(1,1)*grad(1,1) + &
                        grad(2,1)*grad(2,1)+grad(3,1)*grad(3,1)
                   sigma(2) = grad(1,1)*grad(1,2) + &
                        grad(2,1)*grad(2,2)+grad(3,1)*grad(3,2)
                   sigma(3) = grad(1,2)*grad(1,2) + &
                        grad(2,2)*grad(2,2)+grad(3,2)*grad(3,2)
                else
                   sigma(1) = grad(1,1)*grad(1,1) + &
                        grad(2,1)*grad(2,1)+grad(3,1)*grad(3,1)
                end if

                ! ndmh: Exchange part
                if (pub_libxc_x_func_id > 0) then
                   call xc_f90_gga_exc_vxc(libxc_func(1), 1, den(1), sigma(1), &
                        ex, x_pot(1), x_dfdsigma(1))
                else
                   ex = 0.0_DP; x_pot(:) = 0.0_DP; x_dfdsigma(:) = 0.0_DP
                end if

                ! ndmh: Correlation part
                if (pub_libxc_c_func_id > 0) then
                   call xc_f90_gga_exc_vxc(libxc_func(2), 1, den(1), sigma(1), &
                        ec, c_pot(1), c_dfdsigma(1))
                else
                   ec = 0.0_DP; c_pot(:) = 0.0_DP; c_dfdsigma(:) = 0.0_DP
                end if

                if (ns==2) then
                   !x_dfdsigma(1) = x_dfdsigma(1)*2.0_DP*sqrt(sigma(1)) &
                   !     + x_dfdsigma(2)*sqrt(sigma(3))
                   !x_dfdsigma(3) = x_dfdsigma(3)*2.0_DP*sqrt(sigma(3)) &
                   !     + x_dfdsigma(2)*sqrt(sigma(1))
                   !c_dfdsigma(1) = c_dfdsigma(1)*2.0_DP*sqrt(sigma(1)) &
                   !     + c_dfdsigma(2)*sqrt(sigma(3))
                   !c_dfdsigma(3) = c_dfdsigma(3)*2.0_DP*sqrt(sigma(3)) &
                   !     + c_dfdsigma(2)*sqrt(sigma(1))

                   if (sigma(1) > 0.0_DP) then
                      density_grad(i1,i2,islab12,:,1) = &
                           -two*(x_dfdsigma(1)+c_dfdsigma(1)) * grad(:,1) &
                           -(x_dfdsigma(2)+c_dfdsigma(2)) * grad(:,2)
                   else
                      density_grad(i1,i2,islab12,:,1) = 0.0_DP
                   end if
                   if (sigma(3) > 0.0_DP) then
                      density_grad(i1,i2,islab12,:,2) = &
                           -two*(x_dfdsigma(3)+c_dfdsigma(3)) * grad(:,2) &
                           -(x_dfdsigma(2)+c_dfdsigma(2)) * grad(:,1)
                   else
                      density_grad(i1,i2,islab12,:,2) = 0.0_DP
                   end if

                else
                   density_grad(i1,i2,islab12,:,1) = &
                        -two*(x_dfdsigma(1)+c_dfdsigma(1)) * grad(:,1)
                end if

                !if (sum(den(1:ns))>0.01_DP) then
                !   write(25,'(3i6,2f20.12,6f14.10,19f20.12)') &
                !       i1,i2,islab12,den(1:ns),grad(1:3,1:ns), &
                !       ex*sum(den(1:ns)),x_pot(1:ns), &
                !       ec*sum(den(1:ns)),c_pot(1:ns), &
                !       x_dfdsigma(1:3),0.0_DP,c_dfdsigma(1:3), &
                !       density_grad(i1,i2,islab12,:,:)
                !end if

             end if

             ! ndmh: Store this point's energy and potential
             xc_energy = xc_energy + (ex + ec)*sum(den(1:ns))
             !c_energy = c_energy + ec*sum(den(1:ns))
             !x_energy = x_energy + ex*sum(den(1:ns))
             xc_pot_fine(i1,i2,islab12,1:ns) = x_pot(1:ns) + c_pot(1:ns)

          end do a1
       end do a2
    end do a3

    if (libxc_gc) then
       ! Calculate divergence of gradient dependent part and add to potential
       call xc_divergence(density_grad,density_aux,recip_work,grid)
       xc_pot_fine = xc_pot_fine + density_aux
    end if

    ! Global reductions
    call comms_reduce('SUM',xc_energy)
    xc_energy = xc_energy * grid%weight

    ! qoh: Calculate potential if necessary
    if (present(xc_potint)) then
       potint = 0.0_DP
       do is=1,pub_cell%num_spins
          do islab12=1,grid%num_my_slabs12
             do i2=1,grid%n2
                do i1=1,grid%n1
                   potint = potint + xc_pot_fine(i1,i2,islab12,is) * &
                        density_fine(i1,i2,islab12,is)
                end do
             end do
          end do
       end do
       call comms_reduce('SUM',potint)
       xc_potint = potint * grid%weight
    end if

  end subroutine xc_libxc
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_gradients(density_fine,density_grad,recip_work,grid, &
       nhat_dim,nhat_den_grad)

    !========================================================================!
    ! This subroutine calculates the gradients of the density in preparation !
    ! for a gradient-corrected functional.                                   !
    !------------------------------------------------------------------------!
    ! Arguments                                                              !
    !   grid (in) : Grid definition                                          !
    !   density_fine (in) : Input density on the fine grid                   !
    !   density_grad (inout) : Density gradient on the fine grid             !
    !   recip_work (inout) : Fine grid recip-space workspace                 !
    !   nhat_den_grad (in) : Compensation density (and gradient) on fine grid!
    !------------------------------------------------------------------------!
    ! Originally written by Peter Haynes in February 2004.                   !
    ! Modified to avoid module-level arrays by Nicholas Hine in October 2010. !
    !========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_node_id
    use constants, only: PI, cmplx_i
    use fourier, only: fourier_apply_cell_forward, fourier_apply_cell_backward
    use rundat, only: pub_nhat_in_xc
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    type(GRID_INFO),intent(in) :: grid
#ifdef WIN32
    real(kind=DP), intent(inout),target :: density_fine(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_cell%num_spins)
    complex(kind=DP),intent(inout),target :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)
    real(kind=DP),intent(inout),target :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_cell%num_spins)
#else
    real(kind=DP), intent(inout) :: density_fine(grid%ld1,&
         grid%ld2,grid%max_slabs12,pub_cell%num_spins)
    complex(kind=DP),intent(inout) :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)
    real(kind=DP),intent(inout) :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_cell%num_spins)
#endif
    integer, intent(in) :: nhat_dim
    real(kind=DP), optional, intent(in) :: nhat_den_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins,0:nhat_dim)

    ! Local Variables
    integer :: is                      ! Loop variable over spin component
    integer :: i3,i2,islab23           ! Loop variables over fine grids
    integer :: k                       ! Loop variables for Cartesian components
    complex(kind=DP) :: dengtmp        ! Temporary copies of reciprocal density
    real(kind=DP) :: gvec(3)
    ! Introduce pointers to avoid copies to stack on Win32
#ifdef WIN32
    real(kind=DP), pointer :: dens_real(:,:,:)
    complex(kind=DP), pointer :: dens_cmpl(:,:,:)
#endif

    spin: do is=1, pub_cell%num_spins

       ! ndmh: subtract off the compensation density (if present) before
       ! ndmh: taking gradient: analytic gradient will be added afterwards
       if (present(nhat_den_grad).and.pub_nhat_in_xc) then
          density_fine(:,:,:,is) = density_fine(:,:,:,is) - &
               nhat_den_grad(:,:,:,is,0)
       end if

       ! Fourier transform to reciprocal space
#ifdef WIN32
       dens_real => density_fine(:,:,:,is)
       dens_cmpl => recip_work(:,:,:,1)
       call fourier_apply_cell_forward(dens_real,dens_cmpl,grid)
#else
       call fourier_apply_cell_forward(density_fine(:,:,:,is), &
            recip_work(:,:,:,1),grid)
#endif

       ! Apply grad in reciprocal space (multiply by iG)
       b1: do islab23=1,grid%num_slabs23
          b2: do i2=1,grid%n2
             b3: do i3=1,grid%n3

                ! Apply grad
                dengtmp = recip_work(i3,i2,islab23,1)

                call cell_grid_recip_pt(gvec,islab23 + &
                     grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)

                cart_comps1: do k=1,3
                   recip_work(i3,i2,islab23,k) = dengtmp * cmplx_i * gvec(k)
                end do cart_comps1

             end do b3
          end do b2
       end do b1

       ! Transform gradient back to real space
       cart_comps2: do k=1,3
#ifdef WIN32
          dens_real => density_grad(:,:,:,k,is)
          dens_cmpl => recip_work(:,:,:,k)
          call fourier_apply_cell_backward(dens_real,dens_cmpl,grid)
#else
          call fourier_apply_cell_backward(density_grad(:,:,:,k,is), &
               recip_work(:,:,:,k),grid)
#endif
       end do cart_comps2

       ! ndmh: add back on the compensation density itself and add on
       ! ndmh: the compensation density gradient (if present)
       if (present(nhat_den_grad).and.pub_nhat_in_xc) then
          density_fine(:,:,:,is) = density_fine(:,:,:,is) + &
               nhat_den_grad(:,:,:,is,0)
          do k=1,3
             density_grad(:,:,:,k,is) = density_grad(:,:,:,k,is) + &
                  nhat_den_grad(:,:,:,is,k)
          end do
       end if

    end do spin

  end subroutine xc_gradients

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xc_divergence(density_grad,density_aux,recip_work,grid)

    !=========================================================================!
    ! This subroutine calculates the divergence of a vector field (stored in  !
    ! the module array density_grad) which is returned in the array           !
    ! density_aux.                                                            !
    !-------------------------------------------------------------------------!
    ! Originally written by Peter Haynes in February 2004.                    !
    ! Modified to avoid module-level arrays by Nicholas Hine in October 2010. !
    !=========================================================================!

    use cell_grid, only: GRID_INFO, cell_grid_recip_pt
    use comms, only: pub_my_node_id
    use constants, only: DP, cmplx_i
    use fourier, only: fourier_apply_cell_forward, fourier_apply_cell_backward
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
#ifdef WIN32
    real(kind=DP),intent(inout),target :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_cell%num_spins)
    complex(kind=DP),intent(inout),target :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)
    real(kind=DP),intent(inout),target :: density_aux(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
#else
    real(kind=DP),intent(inout) :: density_grad(grid%ld1,grid%ld2, &
         grid%max_slabs12,3,pub_cell%num_spins)
    complex(kind=DP),intent(inout) :: recip_work(grid%ld3,grid%ld2, &
         grid%max_slabs23,3)
    real(kind=DP),intent(inout) :: density_aux(grid%ld1,grid%ld2, &
         grid%max_slabs12,pub_cell%num_spins)
#endif

    ! Local Variables
    integer :: is                      ! Loop variable over spin component
    integer :: i2,i3,islab23           ! Loop variables over fine grids
    integer :: k                       ! Loop variables for Cartesian components
    complex(kind=DP) :: div
    real(kind=DP) :: gvec(3)
#ifdef WIN32
    real(kind=DP), pointer :: dens_real(:,:,:)
    complex(kind=DP), pointer :: dens_cmpl(:,:,:)
#endif

    spin: do is=1, pub_cell%num_spins

       ! Loop over Cartesian components of the vector field
       cart_comps1: do k=1,3

          ! Fourier transform to reciprocal space
#ifdef WIN32
          dens_real => density_grad(:,:,:,k,is)
          dens_cmpl => recip_work(:,:,:,k)
          call fourier_apply_cell_forward(dens_real,dens_cmpl,grid)
#else
          call fourier_apply_cell_forward(density_grad(:,:,:,k,is), &
               recip_work(:,:,:,k),grid)
#endif
       end do cart_comps1

       ! Apply divergence in reciprocal space (multiply by iG)

       b1: do islab23=1,grid%num_slabs23
          b2: do i2=1,grid%n2
             b3: do i3=1,grid%n3

                ! Calculate divergence
                div = 0.0_DP

                call cell_grid_recip_pt(gvec,islab23 + &
                     grid%first_slab23(pub_my_node_id) - 1,i2,i3,grid)

                ! Loop over Cartesian components of the vector field
                cart_comps2: do k=1,3

                   div = div + recip_work(i3,i2,islab23,k) * cmplx_i * gvec(k)

                end do cart_comps2

                recip_work(i3,i2,islab23,1) = div

             end do b3
          end do b2
       end do b1

       ! Transform divergence back to real space
#ifdef WIN32
       dens_real => density_aux(:,:,:,is)
       dens_cmpl => recip_work(:,:,:,1)
       call fourier_apply_cell_backward(dens_real,dens_cmpl,grid)
#else
       call fourier_apply_cell_backward(density_aux(:,:,:,is),&
            recip_work(:,:,:,1),grid)
#endif
    end do spin

  end subroutine xc_divergence

end module xc
