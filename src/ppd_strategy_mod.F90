! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                        PPD strategy module                     !
!----------------------------------------------------------------!
! Written by Chris-Kriton Skylaris, 6/11/2004                    !
!================================================================!




module ppd_strategy

  implicit none

  private


  public :: ppd_strategy_determine
  public :: ppd_strategy_fftbox_spec
  public :: ppd_strategy_check_and_print

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ppd_strategy_determine(d1, d2, d3, n_ppd_pt1, n_ppd_pt2, n_ppd_pt3, & !output
       a1, a2, a3, elements)                                                       !input
    !====================================================================!
    ! This subroutine determines the PSINC grid spacing and the number   !
    ! of points per PPD.                                                 !
    !--------------------------------------------------------------------!
    ! NOTE. The grid spacing is determined according to the following    !
    ! algorithm: First it is calculated according to the given kinetic   !
    ! energy cutoff based on the condition of equal volumes between the  !
    ! ONETEP cube of plane wave vectors and the sphere of wave vectors   !
    ! of conventional plane wave codes. Then the grid spacing in each    !
    ! lattice direction is made finer until it corresponds to an odd     !
    ! number of points with as many prime factors (for efficient FFTs)   !
    ! as possible. If a nonzero value is given in the input for any      !
    ! lattice vector direction, this value is used instead and overrides !
    ! what is selected by the previous procedure. The number of points   !
    ! per ppd is detrmined so that an integer number of ppds fits in     !
    ! each lattice vector direction and again there, the value obtained  !
    ! can be overriden by a non-zero value specified in the input.       !
    !--------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 06/11/2004.                    !
    ! Modified by Chris-Kriton Skylaris on 07/08/2005 to add even psinc  !
    ! grid capability.                                                   !
    !====================================================================!

    use comms, only: comms_abort, comms_barrier, pub_on_root
    use constants, only: DP, VERBOSE, stdout, PI
    use geometry, only: point, operator(.cross.), operator(.dot.),operator(*), &
        geometry_magnitude
    use fourier, only: fourier_size
    use ion, only : ELEMENT
    use rundat, only: cutoff_energy, pub_output_detail, ppd_npoints, &
        psinc_spacing, pub_even_psinc_grid, pub_odd_psinc_grid, &
        pub_cond_calculate
    implicit none


    real(kind=DP), intent(out) :: d1 ! Standard grid spacing in direction 1
    real(kind=DP), intent(out) :: d2 ! Standard grid spacing in direction 2
    real(kind=DP), intent(out) :: d3 ! Standard grid spacing in direction 3
    integer, intent(out) :: n_ppd_pt1 ! Number of points per ppd in direction 1
    integer, intent(out) :: n_ppd_pt2 ! Number of points per ppd in direction 2
    integer, intent(out) :: n_ppd_pt3 ! Number of points per ppd in direction 3
    type(POINT), intent(in) :: a1  ! Lattice vector 1
    type(POINT), intent(in) :: a2  ! Lattice vector 2
    type(POINT), intent(in) :: a3  ! Lattice vector 3
    type(ELEMENT), intent(in) :: elements(:)

    ! cks: <<local variables>>
    real(kind=DP) :: ddd    ! Grid-spacing corresponding to KE cutoff
    real(kind=DP) :: min_fftbox_width
    integer :: n_lat_pt1     ! Number of grid points along a1
    integer :: n_lat_pt2     ! Number of grid points along a2
    integer :: n_lat_pt3     ! Number of grid points along a3
    type(POINT) :: b1
    type(POINT) :: b2
    type(POINT) :: b3

    ! cks: make output info more clear
    if ( pub_on_root .and. pub_output_detail == VERBOSE) write(stdout,*)'   '


    b1=(2.0_DP*PI/( a1.DOT.(a2.CROSS.a3) ) ) * (a2.CROSS.a3)
    b2=(2.0_DP*PI/( a1.DOT.(a2.CROSS.a3) ) ) * (a3.CROSS.a1)
    b3=(2.0_DP*PI/( a1.DOT.(a2.CROSS.a3) ) ) * (a1.CROSS.a2)


    ! cks: grid spacing in 1D corresponding to KE energy cutoff
    if (cutoff_energy > 0.0_DP) then
       ddd = ( (6.0 *PI*PI)**(1.0_DP/3.0_DP) ) &
            /sqrt(2.0_DP*cutoff_energy)
    else
       if (pub_on_root) write(stdout,*)'ERROR: Kinetic energy cutoff is',cutoff_energy,'Eh'
       if (pub_on_root) write(stdout,*)'in ppd_strategy_determine. ONETEP stops.'
       call comms_barrier
       call comms_abort
    endif




    ! cks: minimum required width of FFT-box along any direction
    min_fftbox_width = 6.0_DP*maxval(elements(:)%radius)
    ! ndmh: also check radius of conduction NGWFs if required
    if (pub_cond_calculate) then
       min_fftbox_width = &
            max(min_fftbox_width,6.0_DP*maxval(elements(:)%radius_cond))
    end if



    ! cks: ------------- DETERMINE d1 and n_ppd_pt1 -----------------------------
    call internal_find_dpsinc_in_latvec(d1, n_lat_pt1, &
         a1, b1, psinc_spacing(1), "a1")

    call internal_find_npts_per_ppd_1d(n_ppd_pt1, &
         n_lat_pt1)
    if (n_ppd_pt1 == 0) then
       if (pub_on_root) &
            write(stdout,*)'Failed to initialise n_ppd_pt1 in ppd_strategy_determine.'
       if (pub_on_root) write(stdout,*)'ONETEP stops.'
       call comms_barrier
       call comms_abort
    end if
    if (ppd_npoints(1) > 0) then
       n_ppd_pt1 =ppd_npoints(1)
       if (pub_on_root .and. pub_output_detail == VERBOSE) &
            write(stdout,'(a,i2,a)') '  Number of PPD points along a1 set to',&
            n_ppd_pt1,' from input file.'
    endif
    ! cks: ---------- END DETERMINE d1 and n_ppd_pt1 ---------------------------



    ! cks: ------------- DETERMINE d2 and n_ppd_pt2 ----------------------------
    call internal_find_dpsinc_in_latvec(d2, n_lat_pt2, &
         a2, b2, psinc_spacing(2), "a2")

    call internal_find_npts_per_ppd_1d(n_ppd_pt2, &
         n_lat_pt2)
    if (n_ppd_pt2 == 0) then
       if (pub_on_root) &
            write(stdout,*)'Failed to initialise n_ppd_pt2 in ppd_strategy_determine.'
       if (pub_on_root) write(stdout,*)'ONETEP stops.'
       call comms_barrier
       call comms_abort
    end if
    if (ppd_npoints(2) > 0) then
       n_ppd_pt2 =ppd_npoints(2)
       if (pub_on_root .and. pub_output_detail == VERBOSE) &
            write(stdout,'(a,i2,a)') '  Number of PPD points along a2 set to',&
            n_ppd_pt2,' from input file.'
    endif
    ! cks: ---------- END DETERMINE d2 and n_ppd_pt2 ---------------------------



    ! cks: ------------- DETERMINE d3 and n_ppd_pt3 ----------------------------
    call internal_find_dpsinc_in_latvec(d3, n_lat_pt3, &
         a3, b3, psinc_spacing(3), "a3")

    call internal_find_npts_per_ppd_1d(n_ppd_pt3, &
         n_lat_pt3)
    if (n_ppd_pt3 == 0) then
       if (pub_on_root) &
            write(stdout,*)'Failed to initialise n_ppd_pt3 in ppd_strategy_determine.'
       if (pub_on_root) write(stdout,*)'ONETEP stops.'
       call comms_barrier
       call comms_abort
    end if
    if (ppd_npoints(3) > 0) then
       n_ppd_pt3 =ppd_npoints(3)
       if (pub_on_root .and. pub_output_detail == VERBOSE) &
            write(stdout,'(a,i2,a)') '  Number of PPD points along a3 set to',&
            n_ppd_pt3,' from input file.'
    endif
    ! cks: ---------- END DETERMINE d3 and n_ppd_pt3 ---------------------------



  contains


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_find_dpsinc_in_latvec(dxx, n_lat_pt, &
         a_lat, b_recip, psinc_spacing, lat_text)
      ! cks: Fixed bugs in the manual specification of psinc-grid, 23/04/2010

      implicit none

      real(kind=DP), intent(out) :: dxx
      integer,        intent(out) :: n_lat_pt
      real(kind=DP), intent(in)  :: psinc_spacing
      type(POINT), intent(in)     :: a_lat
      type(POINT), intent(in)     :: b_recip
      character(len=*), intent(in):: lat_text

      ! cks: << local variables >>
      real(kind=DP) :: len_a
      real(kind=DP) :: plane_dist
      logical        :: force_odd


      len_a =geometry_magnitude(a_lat)
      n_lat_pt =nint(len_a/ddd)


      ! cks: distance between consecutive planes defined by other two lattice vectors
      plane_dist =( a_lat .DOT. b_recip)/geometry_magnitude(b_recip)
      ! cks: decide whether FFT box should coincide with simulation cell
      ! cks: and force odd number of points if true
      force_odd =.false.
      if ( plane_dist < (min_fftbox_width+2.0_DP)) force_odd =.true.
      
      if (force_odd.and.pub_even_psinc_grid) then

            if (pub_on_root) write(stdout,'(2a)')'FFT box needs to coincide &
                 &with simulation cell along ', lat_text
            if (pub_on_root) write(stdout,'(a)')'but even_psinc_grid is &
                 &specified. These are incompatible.'

            call comms_abort
      end if

      if ( (psinc_spacing > 0.0_DP) ) then
         ! cks: apply psinc_spacing settings from input file

         dxx =psinc_spacing
         if (pub_on_root .and. pub_output_detail == VERBOSE) &
              write(stdout,'(3a,f16.12,a)') &
              '  PSINC grid spacing along ',lat_text,' set to ',dxx,'a0 from input file.'

         n_lat_pt =nint(len_a/dxx)

         if ( ( abs(mod(real(n_lat_pt, kind=DP), 2.0_DP)) < tiny(1.0_DP) ) .and. force_odd) then

            if (pub_on_root) write(stdout,'(2a)')'FFT box needs to coincide &
                 &with simulation cell along ', lat_text
            if (pub_on_root) write(stdout,'(a)')'Adjust value of psinc_spacing to &
                 &produce an odd number of psinc points.'

            call comms_abort
         endif

         ! ndmh: enforce chosen psinc grid spacing
         dxx =len_a/real(n_lat_pt, kind=DP)

      else

         ! cks: determine optimum number of psinc-grid points
         if (pub_odd_psinc_grid .or. force_odd) then
            call fourier_size(n_lat_pt, 1)
         else if (pub_even_psinc_grid) then
            call fourier_size(n_lat_pt, 0)
         else
            call fourier_size(n_lat_pt)
         endif

         ! cks: optimum psinc grid spacing
         dxx =len_a/real(n_lat_pt, kind=DP)

      endif


    end subroutine internal_find_dpsinc_in_latvec


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_find_npts_per_ppd_1d(n_pt_in_ppd, &  ! output
         n_pt_in_lat)                                        ! input

      implicit none
      !====================================================================!
      ! This subroutine subdivides the number of psinc grid points along   !
      ! a lattice vector into an odd-integer number of ppds. The odd number!
      ! of grid points in the ppd is returned.                             !
      !--------------------------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 10/11/2004.                    !
      ! Modified by Chris-Kriton Skylaris on 07/08/2005 to add even psinc  !
      ! grid capability.                                                   !
      !====================================================================!
      integer, intent(out) :: n_pt_in_ppd
      integer, intent(in)  :: n_pt_in_lat


      ! cks: <<local variables>>
      integer :: row
      integer, parameter :: num_fac =11
      ! cks: PPD sizes in order of preference
      integer, parameter :: primes(num_fac) = (/ 5, 6, 7, 8, 9, 10, 11, 12, 4, 3, 13 /)


      n_pt_in_ppd =0
      ppd_npt_finder: do row= 1, num_fac

         if ( mod( n_pt_in_lat, primes(row)) == 0 ) then
            n_pt_in_ppd =primes(row)
            exit ppd_npt_finder
         endif

      enddo ppd_npt_finder



    end subroutine internal_find_npts_per_ppd_1d

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  end subroutine ppd_strategy_determine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ppd_strategy_fftbox_spec(n1, n2, n3, &  ! output
       elements)                                     ! input
    !==================================================================!
    ! This subroutine returns the number of grid points in each        !
    ! lattice vector direction for the FFTbox.                         !
    !------------------------------------------------------------------!
    ! WARNING: At present, the number of points returned by this       !
    !          subroutine is correct only when orthorhombic simulation !
    !          cells are used. For other types of simulation cells,    !
    !          where either of the alpha or beta or gamma angles are   !
    !          not 90 degrees, this subroutine will underestimate      !
    !          the size of the FFTbox, and the only way to set its     !
    !          size right is by the user specifying it in the          !
    !          input file.                                             !
    !------------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris in 2001 as           !
    ! basis_pair_box_specification to work with a "pair-box".          !
    ! Modified by Arash Mostofi in April 2003 to work with             !
    ! a "triple-box".                                                  !
    ! Modified by Chris-Kriton Skylaris on 18/11/2003 so that it       !
    ! runs in parallel.                                                !
    ! Modified by Peter D. Haynes in June 2004 so that it accepts      !
    ! larger, user-defined FFTbox dimensions.                          !
    ! Modified by Chris-Kriton Skylaris on 06/11/2004 in order to      !
    ! work as part of ppd_strategy.                                    !
    ! Modified by Chris-Kriton Skylaris on 11/02/2008 in order to      !
    ! work with non-orthogonal lattice vectors.                        !
    !==================================================================!

    use comms, only: pub_on_root, comms_abort
    use constants, only: DP, stdout, VERBOSE
    use fourier, only: fourier_size
    use geometry, only: geometry_magnitude, operator(.DOT.)
    use ion, only: ELEMENT
    use rundat, only: fftbox_pref, pub_ngwf_halo, pub_output_detail, &
         pub_cond_calculate
    use simulation_cell, only: pub_cell
    implicit none

    ! Arguments
    integer, intent(out)        :: n1, n2, n3
    type(ELEMENT), intent(in)   :: elements(:) ! elements for all atoms

    ! Local variables
    integer :: same_n ! common value for pairs of n1, n2, n3
    real(kind=DP) :: max_rad  ! maximum ngwf radius
    integer :: n1_temp ! temporary storage for n1
    integer :: n2_temp ! temporary storage for n2
    integer :: n3_temp ! temporary storage for n3
    integer :: odd_cell_pt1 ! highest odd number <= pub_cell%total_pt1
    integer :: odd_cell_pt2 ! highest odd number <= pub_cell%total_pt2
    integer :: odd_cell_pt3 ! highest odd number <= pub_cell%total_pt3
    real(kind=DP)  :: proj1 ! cosine of angle between a1 and a2 x a3
    real(kind=DP)  :: proj2 ! cosine of angle between a2 and a3 x a1
    real(kind=DP)  :: proj3 ! cosine of angle between a3 and a2 x a1


    ! cks: Determine the maximum number of grid points in each lattice vector direction.
    max_rad= maxval(elements(:)%radius)

    ! ndmh: Also account for possibility of conduction NGWFs being bigger
    ! ndmh: than standard ones
    if (pub_cond_calculate) then
       max_rad = max(max_rad,maxval(elements(:)%radius_cond))
    end if

    call internal_check_rad_vs_grid(max_rad, pub_cell%d1)
    call internal_check_rad_vs_grid(max_rad, pub_cell%d2)
    call internal_check_rad_vs_grid(max_rad, pub_cell%d3)

    proj1 =( pub_cell%a1_unit .DOT. pub_cell%b1)&
         /geometry_magnitude(pub_cell%b1)
    proj2 =( pub_cell%a2_unit .DOT. pub_cell%b2)&
         /geometry_magnitude(pub_cell%b2)
    proj3 =( pub_cell%a3_unit .DOT. pub_cell%b3)&
         /geometry_magnitude(pub_cell%b3)

    n1 =3*( int(2.0_DP*real(max_rad, kind=DP)/(abs(proj1)*pub_cell%d1)) +1)
    n2 =3*( int(2.0_DP*real(max_rad, kind=DP)/(abs(proj2)*pub_cell%d2)) +1)
    n3 =3*( int(2.0_DP*real(max_rad, kind=DP)/(abs(proj3)*pub_cell%d3)) +1)

    ! cks: incease points to account for halo region if there is one
    if (pub_ngwf_halo > 0.0_DP) then
       n1=n1 +int(6.0_DP*pub_ngwf_halo/pub_cell%d1)
       n2=n2 +int(6.0_DP*pub_ngwf_halo/pub_cell%d2)
       n3=n3 +int(6.0_DP*pub_ngwf_halo/pub_cell%d3)
    endif


    ! ensures that all dimensions are odd
    n1 = n1 + 1 - mod(n1,2)
    n2 = n2 + 1 - mod(n2,2)
    n3 = n3 + 1 - mod(n3,2)

    ! ndmh: find odd number equal to or below cell dimensions
    odd_cell_pt1 = pub_cell%total_pt1 - 1 + mod(pub_cell%total_pt1,2)
    odd_cell_pt2 = pub_cell%total_pt2 - 1 + mod(pub_cell%total_pt2,2)
    odd_cell_pt3 = pub_cell%total_pt3 - 1 + mod(pub_cell%total_pt3,2)

    ! cks: Increase dimensions until Fourier transform efficiency
    ! cks: provided it does not exceed the simulation cell
    ! ndmh: added protection against cell size becoming even
    n1_temp =n1
    call fourier_size(n1, 1)
    if (n1 > odd_cell_pt1) n1 = n1_temp
    n2_temp =n2
    call fourier_size(n2, 1)
    if (n2 > odd_cell_pt2) n2 = n2_temp
    n3_temp =n3
    call fourier_size(n3, 1)
    if (n3 > odd_cell_pt3) n3 = n3_temp

    ! cks: if grid spacing is the same along two directions make sure FFT box
    ! cks: points are the same
    ! ndmh: added protection against dimensions becoming even
    ! jd: added missing abs that ndmh had spotted
    if ( abs(pub_cell%d1 -pub_cell%d2) < epsilon(1.0_DP) ) then
       same_n =max(n1, n2)
       if (same_n <= odd_cell_pt1) n1 = same_n
       if (same_n <= odd_cell_pt2) n2 = same_n
    endif
    if ( abs(pub_cell%d1 -pub_cell%d3) < epsilon(1.0_DP) ) then
       same_n =max(n1, n3)
       if (same_n <= odd_cell_pt1) n1 = same_n
       if (same_n <= odd_cell_pt3) n3 = same_n
    endif
    if ( abs(pub_cell%d3 -pub_cell%d2) < epsilon(1.0_DP) ) then
       same_n =max(n3, n2)
       if (same_n <= odd_cell_pt3) n3 = same_n
       if (same_n <= odd_cell_pt2) n2 = same_n
    endif


    ! cks: If the FFT box is calculated to be larger than the
    ! cks: simulation cell along a lattice vector, set it
    ! cks: equal to the simulation cell along that vector.
    if (n1 > pub_cell%total_pt1) n1 =pub_cell%total_pt1
    if (n2 > pub_cell%total_pt2) n2 =pub_cell%total_pt2
    if (n3 > pub_cell%total_pt3) n3 =pub_cell%total_pt3


    ! pdh: increase FFT box dimensions to specified preference
    if (pub_on_root .and. pub_output_detail == VERBOSE .and. &
         (fftbox_pref(1) > n1 .or. fftbox_pref(2) > n2 .or. &
         fftbox_pref(3) > n3)) &
         write(stdout,'(/a)') 'WARNING in ppd_strategy_fftbox_spec:'
    if (fftbox_pref(1) > n1) then
       if (pub_on_root .and. pub_output_detail == VERBOSE) &
            write(stdout,'(2(a,i3))') '  FFT-box dimension 1 increased from ',n1, &
            ' to preferred size ',fftbox_pref(1)
       n1 = fftbox_pref(1)
    end if
    if (fftbox_pref(2) > n2) then
       if (pub_on_root.and. pub_output_detail == VERBOSE) &
            write(stdout,'(2(a,i3))') '  FFT-box dimension 2 increased from ',n2, &
            ' to preferred size ',fftbox_pref(2)
       n2 = fftbox_pref(2)
    end if
    if (fftbox_pref(3) > n3) then
       if (pub_on_root.and. pub_output_detail == VERBOSE) &
            write(stdout,'(2(a,i3))') '  FFT-box dimension 3 increased from ',n3, &
            ' to preferred size ',fftbox_pref(3)
       n3 = fftbox_pref(3)
    end if


    ! cks: declare that FFT box has been specified
    if (pub_on_root .and. pub_output_detail == VERBOSE) &
         write(stdout,'(a)') '... FFT-box specified'


  contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_check_rad_vs_grid(max_rad_xx, d_xx)

      implicit none
      real(kind=DP), intent(in) :: max_rad_xx  ! maximum ngwf radius
      real(kind=DP), intent(in) :: d_xx        ! psinc grid spacing


      if (max_rad_xx < d_xx) then
         write(stdout, *)"ERROR: Maximum NGWF radius is"&
              &" smaller than psinc grid spacing."
         call comms_abort
      endif


    end subroutine internal_check_rad_vs_grid

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  end subroutine ppd_strategy_fftbox_spec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ppd_strategy_check_and_print
    !=============================================================!
    ! This subroutine prints information on the PSINC grids used  !
    ! by ONETEP. It also checks the sizes of the FFTbox and PPDs  !
    ! and stops if they are not sensible.                         !
    !-------------------------------------------------------------!
    ! Originally written by Chris-Kriton Skylaris in 2001 for the !
    ! ONES code.                                                  !
    ! Modified by Chris-Kriton Skylaris on 6/11/2004 to become    !
    ! part of the ppd_strategy module.                            !
    !=============================================================!
    use comms, only : pub_on_root
    use constants, only: DP, stdout, PI, VERBOSE, HARTREE_IN_eVs
    use rundat, only: pub_output_detail
    use simulation_cell, only: pub_cell, pub_fftbox

    implicit none

    ! cks: internal declarations
    real(kind=DP) :: ke1    ! Kinetic energy cutoff along direction a1
    real(kind=DP) :: ke2    ! Kinetic energy cutoff along direction a2
    real(kind=DP) :: ke3    ! Kinetic energy cutoff along direction a3
    real(kind=DP) :: ke_fac ! Grid space to kinetic energy conversion factor


    ! cks: print PSINC-grids related info
    if (pub_on_root) then

       write(stdout,'(a)') ''
       write(stdout,'(a)') '============================== PSINC grid sizes ================================'
       write(stdout,'(3(a,i4))') '                      Simulation cell: ',&
            pub_cell%total_pt1 ,' x',pub_cell%total_pt2 ,' x',pub_cell%total_pt3

       write(stdout,'(3(a,i4))') '                              FFT-box: ',&
            pub_fftbox%total_pt1 ,' x',pub_fftbox%total_pt2 ,' x',pub_fftbox%total_pt3

       write(stdout,'(3(a,i4))') '                                  PPD: ',&
            pub_cell%n_pt1 ,' x',pub_cell%n_pt2 ,' x',pub_cell%n_pt3

       if (pub_output_detail == VERBOSE ) then
          ! cks: kinetic energy in Eh corresponding to each grid space
          ke_fac =( (6.0_DP*(PI**2))**(2.0_DP/3.0_DP) ) / 2.0_DP
          ke1 = ke_fac/(pub_cell%d1 **2)
          ke2 = ke_fac/(pub_cell%d2 **2)
          ke3 = ke_fac/(pub_cell%d3 **2)
          write(stdout,'(a,f16.12,a, f10.5,a,f10.5,a)') 'Grid space d1=',pub_cell%d1,'a0 (KE cutoff=',ke1, &
               'Eh =',ke1*HARTREE_IN_eVs,'eV)'
          write(stdout,'(a,f16.12,a, f10.5,a,f10.5,a)') 'Grid space d2=',pub_cell%d2,'a0 (KE cutoff=',ke2, &
               'Eh =',ke2*HARTREE_IN_eVs,'eV)'
          write(stdout,'(a,f16.12,a, f10.5,a,f10.5,a)') 'Grid space d3=',pub_cell%d3,'a0 (KE cutoff=',ke3, &
               'Eh =',ke3*HARTREE_IN_eVs,'eV)'
       endif

       write(stdout,'(a)') '================================================================================'

    end if


    ! cks: make sure the size of the fftbox does not exceed the simulation cell
    call internal_check_fftbox_size


    ! cks: make sure that an integer number of grid points and ppds fit
    ! cks: along each lattice vector.
    call internal_check_ppd_shape( &
         pub_cell%a1, pub_cell%d1, pub_cell%n_pt1, 1)
    call internal_check_ppd_shape( &
         pub_cell%a2, pub_cell%d2, pub_cell%n_pt2, 2)
    call internal_check_ppd_shape( &
         pub_cell%a3, pub_cell%d3, pub_cell%n_pt3, 3)



  contains


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_check_fftbox_size
      !================================================================!
      ! This subroutine checks that the size of the fftbox does not    !
      ! exceed the size of the simulation cell.                        !
      !----------------------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 6/11/2004.                 !
      !================================================================!

      use comms, only: comms_abort, comms_barrier, pub_on_root
      implicit none



      if ( pub_fftbox%total_pt1 > pub_cell%total_pt1 .or. &
           pub_fftbox%total_pt2 > pub_cell%total_pt2 .or. &
           pub_fftbox%total_pt3 > pub_cell%total_pt3) then

         if (pub_on_root) then
            write(stdout,'(a/a)') 'Error in ppd_strategy_check_and_print:', &
                 '  one or more of the dimensions of the FFT-box is larger &
                 & than the simulation cell dimensions.'
         end if
         call comms_barrier
         call comms_abort
      end if


      ! cks: make sure that the FFT box has an odd number of points
      ! cks: along each lattice vector direction
     if ( (mod(pub_fftbox%total_pt1, 2) == 0) .or. &
          (mod(pub_fftbox%total_pt2, 2) == 0) .or. &
          (mod(pub_fftbox%total_pt3, 2) == 0) ) then

         if (pub_on_root) then
            write(stdout,'(a/a)') 'Error in ppd_strategy_check_and_print:', &
                 '  the number of psinc grid points along one or more of the &
                 & lattice vectors of the FFT-box is even.'
         end if
         call comms_barrier
         call comms_abort
      end if

    end subroutine internal_check_fftbox_size

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine internal_check_ppd_shape(lat_vec, dx, ppd_npt, &  ! input
         direction_counter)                                      ! input
      !================================================================!
      ! This subroutine checks in a lattice vector direction that the  !
      ! similation cell grid contains an integer number of grid points !
      ! and an integer number of ppds.                                 !
      !----------------------------------------------------------------!
      ! Written by Chris-Kriton Skylaris on 21/6/2001.                 !
      !================================================================!
      use comms, only : pub_on_root, comms_abort, comms_barrier
      use geometry, only: point, geometry_magnitude
      use constants, only: DP, stdout
      implicit none

      type(POINT), intent(in) :: lat_vec
      real(kind=DP), intent(in) :: dx
      integer, intent(in) :: ppd_npt, direction_counter

      ! cks: internal declarations
      real(kind=DP) :: total_npt_real
      integer :: total_npt

      if (ppd_npt < 0) then
         if (pub_on_root) write(stdout,'(a/a,i1,a,i3,a)') &
              'Error in ppd_strategy_check_and_print:', &
              '  the number of points per ppd in direction ',direction_counter, &
              ' is ',ppd_npt, ' < 0'
         call comms_barrier
         call comms_abort
      end if

      total_npt_real = geometry_magnitude(lat_vec) / dx
      total_npt = nint(total_npt_real)

      ! cks: make sure that an essentially integral number of grid points fits
      ! cks: in the current edge of the simulation cell
      if (abs(total_npt_real - real(total_npt,kind=DP)) > 1.0e-5_DP) then
         if (pub_on_root) write(stdout,'(a/a,i1,a,f7.2,a)') &
              'Error in ppd_strategy_check_and_print:', &
              '  the number of points in direction ',direction_counter, &
              ' of the simulation cell (',total_npt_real, &
              ') is not an integer - adjust grid spacing'
         call comms_barrier
         call comms_abort
      end if


      ! cks: make sure that an integer number of ppds fits in the current edge of
      !      the simulation cell.
      if (mod(total_npt, ppd_npt) /= 0) then
         if (pub_on_root) write(stdout,'(a/a,i1,a,f7.2,a)') &
              'Error in ppd_strategy_check_and_print:', &
              '  the number of ppds in direction ',direction_counter, &
              ' of the simulation cell (', &
              real(total_npt,kind=DP) / real(ppd_npt,kind=DP), &
              ') is not an integer - adjust points in ppd'
         call comms_barrier
         call comms_abort
      end if

    end subroutine internal_check_ppd_shape

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  end subroutine ppd_strategy_check_and_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module ppd_strategy




