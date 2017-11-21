! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   The subroutines in this file were written by
!
!   Chris-Kriton Skylaris and Arash A. Mostofi
!
!   TCM Group, Cavendish laboratory
!   Madingley Road
!   Cambridge CB3 0HE
!   UK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



module services

  use constants, only: DP

  implicit none

  private

  public :: services_symm_fd_sec_der
  public :: services_symm_fd_sec_der_KAX
  public :: services_equally_spaced_numbers
  public :: services_locate_interp
  public :: services_linear_interpolation
  public :: services_1d_interpolation
  public :: services_polak_cg_coeff
  public :: services_line_search_parabola
  public :: services_parabolic_step
  public :: services_polynomial_step
  public :: services_cubic_fit_minimum
  public :: services_cubic_minimum
  public :: services_cubic_fit_maximum
  public :: services_flush
  public :: services_sbessj
  public :: services_radial_transform
  public :: services_regular_transform
  public :: services_radial_integral
  public :: services_radial_integral_rmax
  public :: services_radial_derivative
  public :: services_regular_integral
  public :: services_analytic_limit
  public :: services_write_xyz
  public :: services_print_num_species
  public :: services_rationalise_coords
  public :: services_maxboltzdist
  public :: services_rms_fit

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine services_write_xyz(elements, xyz_root, title_line)

    !==================================================================!
    ! This subroutine outputs the cartesian coordinates of the atoms   !
    ! in a .xyz format file in Angstroms.                              !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/03/2007.                  !
    !==================================================================!

    use comms, only: comms_abort, pub_on_root
    use constants, only: DP, stdout, stderr, ANGSTROM
    use ion, only: ELEMENT
    use rundat, only: pub_rootname
    use simulation_cell, only: pub_cell
    use utils, only: utils_unit

    implicit none

    type(ELEMENT), intent(in) :: elements(:)   ! elements for all atoms on all nodes
    character(len=*), intent(in) :: xyz_root   ! root name of file to output
    character(len=*), intent(in) :: title_line ! title line to output in file

    ! Local Variables
    integer :: xyz_unit                ! unit for output file
    integer :: ierr                    ! error flag
    integer :: row                     ! atom counter
    character(len=128) :: xyz_filename ! name of file to output
    character(len=10)  :: cbuf         ! character buffer

#ifdef DEBUG
    if (pub_on_root) &
         write(stdout,'(a)') 'DEBUG: Entering services_write_xyz'
#endif

    if (pub_on_root) then

       write(xyz_filename,'(2a)') trim(xyz_root), ".xyz"

       write(stdout,'(3a)',advance='no') &
       ' Writing xyz coordinates to file "',trim(xyz_filename),'" ...'


       ! cks: find free unit to open for the write operation
       xyz_unit =utils_unit()

       ! cks: open file
       open(unit=xyz_unit, form="formatted", file= trim(xyz_filename), &
            action="write", iostat=ierr, position='append')
       if (ierr /= 0) then
          write(stderr,'(3a,i6)') 'Error in services_write_xyz: &
               &opening file "',trim(xyz_filename),'" failed with code ',ierr
          call comms_abort
       end if

       ! cks: write the number of atoms
       write(cbuf,'(i9)') pub_cell%nat
       write(xyz_unit,'(a)',iostat=ierr) adjustl(cbuf)
       if (ierr /= 0) then
          write(stderr,'(a)') 'Error in services_write_xyz: &
               &writing pub_cell%nat failed'
          call comms_abort
       end if

       ! cks: write the title line
       write(xyz_unit,'(a)',iostat=ierr) adjustl(title_line)
       if (ierr /= 0) then
          write(stderr,'(a)') 'Error in services_write_xyz: &
               &writing title_line failed'
          call comms_abort
       end if

       ! cks: write the atomic coordinates
       do row =1, pub_cell%nat
          write(xyz_unit,'(a,3f12.6)',iostat=ierr)adjustl(elements(row)%symbol), &
               elements(row)%centre%x/ANGSTROM, &
               elements(row)%centre%y/ANGSTROM, &
               elements(row)%centre%z/ANGSTROM
          if (ierr /= 0) then
             write(stderr,'(a)') 'Error in services_write_xyz: &
                  &writing elements(row)%centre%x/ANGSTROM failed'
             call comms_abort
          end if
       end do

       ! close file
       close(unit=xyz_unit, iostat=ierr)
       if (ierr /= 0) then
          write(stderr,'(3a,i6)') 'Error in services_write_xyz: &
               &closing file "',trim(xyz_filename),'" failed with code ',ierr
          call comms_abort
       end if

       write(stdout,'(a)') ' done'

    endif


#ifdef DEBUG
    if (pub_on_root) &
         write(stdout,'(a)') 'DEBUG: Leaving services_write_xyz'
#endif

    return
  end subroutine services_write_xyz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_print_num_species(elements)

    !=========================================================-!
    ! This subroutine prints out the number of atoms and NGWFs !
    ! for each species.                                        !
    !----------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/03/2007.          !
    ! Conduction NGWFs added by Nicholas Hine on 12/10/2011.   !
    !==========================================================!

    use constants, only: stdout
    use ion, only: ELEMENT
    use rundat, only: pub_cond_calculate, pub_paw
    use simulation_cell, only: pub_cell

    implicit none

    ! Arguments
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)

    ! Local Variables
    integer :: n_atoms
    integer :: n_ngwfs
    integer :: n_ngwfs_cond
    integer :: n_projs
    integer :: species_counter
    integer :: row
    character(len=2) :: el_symbol

    write(stdout,"(a)")"------------------------- Atom counting information&
         & ---------------------------"

    if (pub_cond_calculate) then
       write(stdout,"(a)")"Symbol    Natoms  Nvalngwfs Ncondngwfs  Nprojs"
    else
       write(stdout,"(a)")"Symbol    Natoms    Nngwfs    Nprojs"
    end if

    do species_counter=1,pub_cell%num_pspecies

       n_atoms = 0
       n_ngwfs = 0
       n_ngwfs_cond = 0
       n_projs = 0
       do row=1,pub_cell%nat

          if (elements(row)%pspecies_number == species_counter) then
             n_atoms = n_atoms + 1
             n_ngwfs = n_ngwfs + elements(row)%nfunctions
             n_ngwfs_cond = n_ngwfs_cond + elements(row)%nfunctions_cond
             if (pub_paw) then
                n_projs = n_projs + elements(row)%npawpws
             else
                n_projs = n_projs + elements(row)%nprojectors
             end if
             el_symbol = elements(row)%symbol
          endif
       enddo

       if (pub_cond_calculate) then
          write(stdout,"(a4,4i10)")el_symbol,n_atoms,n_ngwfs,n_ngwfs_cond,n_projs
       else
          write(stdout,"(a4,3i10)")el_symbol,n_atoms,n_ngwfs,n_projs
       end if

    enddo

    if (pub_paw) then
       n_projs = pub_cell%num_pawpws
    else
       n_projs = pub_cell%num_projectors
    end if
    if (pub_cond_calculate) then
       write(stdout,"(a)")".......   ......    ......    ......    ......"
       write(stdout,"(a7,i7,3i10)")"Totals:",pub_cell%nat, pub_cell%num_ngwfs, &
            pub_cell%num_ngwfs_cond, n_projs
       write(stdout,"(a)")"---------------------------------------&
            &----------------------------------------"
    else
       write(stdout,"(a)")".......   ......    ......    ......"
       write(stdout,"(a7,i7,2i10)")"Totals:",pub_cell%nat, pub_cell%num_ngwfs, &
            n_projs
       write(stdout,"(a)")"---------------------------------------&
            &----------------------------------------"
    end if
    write(stdout,*)

  end subroutine services_print_num_species


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function services_sbessj(l,x)

    !===============================================!
    ! SPHERICAL BESSEL function OF THE FIRST KIND.  !
    !-----------------------------------------------!
    ! Written by Peter D. Haynes in early 2004.     !
    !===============================================!

    ! Arguments
    integer, intent(in) :: l
    real(kind=DP), intent(in) :: x

    ! Local variables
    integer :: j
    real(kind=DP), parameter :: third = 1.0_DP / 3.0_DP
    real(kind=DP), parameter :: ftnth = 1.0_DP / 14.0_DP
    real(kind=DP) :: x2, sb0, sb1, by, bym, byp, ux
    real(kind=DP) :: sbessj

    x2 = 0.5_DP*x*x
    if (abs(x) > 0.001_DP) then
       sb0 = sin(x)/x
    else
       sb0 = 1.0_DP - third*x2*(1.0_DP - 0.1_DP*x2)
    end if
    if (l == 0) then
       sbessj = sb0
    else
       if (abs(x) > 0.001_DP) then
          sb1 = (sb0 - cos(x)) / x
       else
          sb1 = third*x*(1.0_DP - (0.2_DP*x2)*(1.0_DP - ftnth*x2))
       end if
       if (l == 1) then
          sbessj = sb1
       else if (x == 0.0_DP) then
          sbessj = 0.0_DP
       else
          by = sb1
          bym = sb0
          ux = 1.0_DP / x
          do j=1,l-1
             byp = real(2*J+1,DP)*ux*by - bym
             bym = by
             by = byp
          end do
          sbessj = by
       end if
    end if

    services_sbessj =sbessj

  end function services_sbessj


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! cks: written by Chris-Kriton Skylaris in 2000
  function services_symm_fd_sec_der(val,delta_space,order)

    use comms, only: comms_abort
    use constants, only: DP, stdout
    implicit none

    real(kind=DP) :: services_symm_fd_sec_der

    ! Arguments
    integer, intent(in) :: order
    real(kind=DP), intent(in) :: val(-order/2:order/2), delta_space

    ! cks: internal variables
    real(kind=DP) :: sd


    if (order.eq.2) then
       sd= val(-1) -2.0_DP*val(0) +val(1)
    elseif (order.eq.4) then
       sd= -val(-2)/12.0_DP + 4.0_DP*val(-1)/3.0_DP  -5.0_DP*val(0)/2.0_DP &
            +4.0_DP*val(1)/3.0_DP -val(2)/12.0_DP
    elseif (order.eq.6) then
       sd= ( val(-3) + val(3) )/90.0_DP -3.0_DP*(val(-2) + val(2) )/20.0_DP &
            +3.0_DP*(val(-1)+val(1))/2.0_DP -49.0_DP*val(0)/18.0_DP
    elseif (order.eq.8) then
       sd= -(val(-4)+val(4))/560.0_DP +8.0_DP*(val(-3)+val(3))/315.0_DP &
            -(val(-2)+val(2))/5.0_DP +8.0_DP*(val(-1)+val(1))/5.0_DP &
            -205.0_DP*val(0)/72.0_DP
    elseif (order.eq.10) then
       sd= (val(-5)+val(5))/3150.0_DP -5.0_DP*(val(-4)+val(4))/1008.0_DP &
            +5.0_DP*(val(-3)+val(3))/126.0_DP -5.0_DP*(val(-2)+val(2))/21.0_DP &
            +5.0_DP*(val(-1)+val(1))/3.0_DP -5269.0_DP*val(0)/1800.0_DP
    elseif (order.eq.12) then
       sd= -(val(-6)+val(6))/16632.0_DP + 2.0_DP*(val(-5)+val(5))/1925.0_DP &
            -(val(-4)+val(4))/112.0_DP +10.0_DP*(val(-3)+val(3))/189.0_DP &
            -15.0_DP*(val(-2)+val(2))/56.0_DP + 12.0_DP*(val(-1)+val(1))/7.0_DP &
            -5369.0_DP*val(0)/1800.0_DP
    elseif (order.eq.16) then
       sd= -(val(-8)+val(8))*2.428127428127428127428127428127428E-6_DP &
            +(val(-7)+val(7))*5.074290788576502862217147931433646E-5_DP &
            -(val(-6)+val(6))*5.180005180005180005180005180005181E-4_DP &
            +(val(-5)+val(5))*3.480963480963480963480963480987078E-3_DP &
            -(val(-4)+val(4))*1.767676767676767676767676767673481E-2_DP &
            +(val(-3)+val(3))*7.542087542087542087542087542076068E-2_DP &
            -(val(-2)+val(2))*0.311111111111111111111111111111111_DP  &
            +(val(-1)+val(1))*1.777777777777777777777777777777777_DP   &
            -val(0)*3.05484410430839005820287374945639_DP
    elseif (order.eq.20) then
       sd= -(val(-10)+val(10))*1.082508822446902942258979410682197E-7_DP &
            +(val(-9)+val(9))*2.672861289992352943849331878227647E-6_DP &
            -(val(-8)+val(8))*3.213698066639243109831345125462773E-5_DP &
            +(val(-7)+val(7))*2.518489913447896641173952098321847E-4_DP &
            -(val(-6)+val(6))*1.456876456876456876456876456876457E-3_DP &
            +(val(-5)+val(5))*6.713286713286713286713286713332222E-3_DP &
            -(val(-4)+val(4))*2.622377622377622377622377622392322E-2_DP &
            +(val(-3)+val(3))*9.324009324009324009324009324009325E-2_DP &
            -(val(-2)+val(2))*0.3409090909090909090909090909090909_DP   &
            +(val(-1)+val(1))*1.8181818181818181818181818181818181_DP   &
            -val(0)*3.09953546233308141622756510748108
    elseif (order.eq.24) then
       sd= -(val(-12)+val(12))*5.136127090629715478281907141780610E-9_DP &
            +(val(-11)+val(11))*1.466979770679032784540683560495355E-7_DP &
            -(val(-10)+val(10))*2.041302350899874119688361174429285E-6_DP &
            +(val(-9)+val(9))*1.848092663366141178318680898655856E-5_DP &
            -(val(-8)+val(8))*1.227970945463205525125029768992618E-4_DP &
            +(val(-7)+val(7))*6.415521674256747233306277976777757E-4_DP &
            -(val(-6)+val(6))*2.765208647561589252051467324785477E-3_DP &
            +(val(-5)+val(5))*1.023917259211376858435681965093730E-2_DP &
            -(val(-4)+val(4))*3.399725274725274725274725274737918E-2_DP &
            +(val(-3)+val(3))*0.107448107448107448107448107448107_DP &
            -(val(-2)+val(2))*0.362637362637362637362637362637_DP &
            +(val(-1)+val(1))*1.846153846153846153846153846153_DP &
            -val(0)*3.1299532768418050158602556492717
    elseif (order.eq.28) then
       sd = -(val(-14)+val(14))*2.543605797264240046387230203548493E-10_DP &
            +(val(-13)+val(13))*8.259945926263993712765159382884103E-9_DP &
            -(val(-12)+val(12))*1.308685182692451503866229939725700E-7_DP &
            +(val(-11)+val(11))*1.349784386777007832086822284940225E-6_DP &
            -(val(-10)+val(10))*1.020774442500112173015659352986046E-5_DP &
            +(val(-9)+val(9))*6.049033733333998062315018388065456E-5_DP &
            -(val(-8)+val(8))*2.934726522187822497420020639834882E-4_DP &
            +(val(-7)+val(7))*1.204692403277100313809734420083823E-3_DP &
            -(val(-6)+val(6))*4.304265565875473796285376112133183E-3_DP &
            +(val(-5)+val(5))*1.377364981080151358789129686981574E-2_DP &
            -(val(-4)+val(4))*4.089052287581699346405228758192856E-2_DP &
            +(val(-3)+val(3))*0.118954248366013083531987019838242_DP &
            -(val(-2)+val(2))*0.379166666666666666666666666666666_DP &
            +(val(-1)+val(1))*1.866666666666666666666666666666666_DP &
            -val(0)*3.15199167800108529601965668779366_DP
    else
       write (stdout,'(a,i6)')  &
            "Finite difference second derivative for accuracy order=",order
       write (stdout,'(a)')  "not available. Program execution stops"
       call comms_abort
       sd = 0.0_DP ! qoh: Prevent compiler warning
    endif


    services_symm_fd_sec_der=sd/(delta_space**2)

  end function services_symm_fd_sec_der


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! cks: written by Chris-Kriton Skylaris, 9/11/2001
  function services_symm_fd_sec_der_KAX(val,delta_space,order)

    use comms, only: comms_abort
    use constants, only: DP, stdout, PI

    implicit none

    real(kind=DP) :: services_symm_fd_sec_der_KAX

    ! Arguments
    integer, intent(in) :: order
    real(kind=DP), intent(in) :: val(-order/2:order/2), delta_space

    ! cks: internal variables
    real(kind=DP) :: sd
    real(kind=DP), parameter :: psq=PI**2

    if (order.eq.12) then
       sd=  (val(-6)+val(6))*(1.0_DP/600.0_DP-psq/4096.0_DP) &
            + (val(-5)+val(5))*( 3.0_DP*psq/1024.0_DP - 31.0_DP/1575.0_DP  )  &
            +(val(-4)+val(4))*( 2647.0_DP/25200.0_DP-33.0_DP*psq/2048.0_DP ) &
            +(val(-3)+val(3))*( 55.0_DP*psq/1024.0_DP - 103.0_DP/315.0_DP ) &
            +(val(-2)+val(2))*( 493.0_DP/840.0_DP - 495.0_DP*psq/4096.0_DP ) &
            +(val(-1)+val(1))*( 26.0_DP/75.0_DP + 99.0_DP*psq/512.0_DP ) &
            -val(0)*(2497.0_DP/1800.0_DP+231.0_DP*psq/1024.0_DP )
    else
       write (stdout,'(a,i6)')  &
            "Finite difference second derivative for accuracy order=",order
       write (stdout,'(a)')  "not available. Program execution stops"
       call comms_abort
       sd = 0.0_DP ! qoh: Prevent compiler warning
    endif


    services_symm_fd_sec_der_KAX=sd/(delta_space**2)

  end function services_symm_fd_sec_der_KAX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_equally_spaced_numbers(numbers,min_number,max_number, &
       num_numbers)

    use comms, only: comms_abort
    use constants, only: DP, stdout

    implicit none

    ! Arguments
    integer, intent(in) :: num_numbers
    real(kind=DP), intent(out) :: numbers(num_numbers)
    real(kind=DP), intent(in) :: min_number, max_number

    ! Local Variables
    integer :: row
    real(kind=DP) :: delta_number

    if (num_numbers<=1) then
       write (stdout,'(a)') 'num_numbers has unpermitted value in &
            &services_equally_spaced_numbers.'
       call comms_abort
    end if

    delta_number = (max_number-min_number)/real(num_numbers-1,kind=DP)

    do row=1,num_numbers
       numbers(row) = min_number + real(row-1,kind=DP)*delta_number
    end do

  end subroutine services_equally_spaced_numbers


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  integer function services_locate_interp(length,points,num_points)

    use constants, only: DP, stdout
    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: num_points
    real(kind=DP), intent(in) :: length
    real(kind=DP), intent(in) :: points(num_points)

    ! Local Variables
    integer :: n_low,n_high,n_mid

    ! CKS: FIRST USE BISECTION METHOD TO FIND THE POINTS BETWEEN WHICH LENGTH LIES
    n_low=0
    n_high=num_points+1
    do
       if (n_high-n_low==1) exit

       n_mid=(n_low+n_high)/2
       if ( length>=points(n_mid) ) then
          n_low=n_mid
       else
          n_high=n_mid
       endif

    enddo

    ! CKS: TAKE CARE OF THE CASE WHERE LENGTH COINCIDES WITH THE FIRST OR THE
    ! CKS: LAST POINT OR STOP IF IT IS OUT OF BOUNDS.
    if (length==points(1)) then
       services_locate_interp=1
    else if (length==points(num_points)) then
       services_locate_interp=num_points-1
    else
       services_locate_interp=n_low
       if ((n_low==0) .or. (n_low==num_points)) then
          call utils_abort('Error in services_locate_interp: length out of &
               &interpolating points range')
       endif

    endif


  end function services_locate_interp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function services_linear_interpolation(x_i,y1,y2,x1,x2)

    ! CKS : 1D REAL LINEAR INTERPOLATION BETWEEN TWO GIVEN POINTS

    use constants, only: DP, stdout
    use utils, only: utils_abort

    implicit none


    ! Arguments
    real(kind=DP), intent(in) :: x_i,y1,y2,x1,x2

    if ((x2-x1)==0.0_DP) then
       call utils_abort('Error in services_linear_interpolation: Attempted &
            &division by zero')
    end if

    services_linear_interpolation = ( y1*(x2-x_i) + y2*(x_i-x1) )/ (x2-x1)


  end function services_linear_interpolation



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function services_1d_interpolation(values, num_values, xx, l_mom)

    !===============================================================!
    ! Interpolation in one dimension suitable for pseudopotentials. !
    !---------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 10/9/2001.                !
    ! Modified by Chris-Kriton Skylaris on 29/3/2004.               !
    ! Cleanup by Nicholas Hine, May 2011.                           !
    !===============================================================!

    use comms, only: pub_on_root, comms_abort
    use constants, only: DP, stdout

    implicit none

    ! Arguments
    real(kind=DP) :: services_1d_interpolation
    integer, intent(in) ::  num_values
    integer, intent(in) ::  l_mom
    real(kind=DP), intent(in) :: values(num_values)
    real(kind=DP), intent(in) :: xx

    ! Local Variables
    integer :: start_index
    real(kind=DP) :: x2
    real(kind=DP) :: off, f0, f1, f2, f3, f4, t0, t1, t2, t3, t4

    ! initialisations
    start_index = int(xx) + 1
    x2 = real(start_index-1,DP)
    off = xx - x2

    if (start_index .lt. 1) then

       write(stdout,*)'start_index=',start_index
       write(stdout,*)'in services_1d_interpolation. Error, ONETEP stops now'
       call comms_abort
       services_1d_interpolation = 0.0_DP !qoh: Prevent compiler warning

    elseif (start_index .eq. 1) then

       ! cks: take into account whether function is even or odd
       f0 = values(3)*( (-1.0_DP)**l_mom)
       f1 = values(2)*( (-1.0_DP)**l_mom)
       f2 = values(1)
       f3 = values(2)
       f4 = values(3)

       ! Quartic interpolation
       t0 = f2
       t1 = 2.0_dp*f0-16.0_dp*f1+16.0_dp*f3-2.0_dp*f4
       t2 = -f0+16.0_dp*f1-30.0_dp*f2+16.0_dp*f3-f4
       t3 = -2.0_dp*f0+4.0_dp*f1-4.0_dp*f3+2.0_dp*f4
       t4 = f0-4.0_dp*f1+6.0_dp*f2-4.0_dp*f3+f4

       services_1d_interpolation = t0+off*(t1+off*(t2+off*(t3+off*t4)))/24.0_dp

    elseif (start_index.gt.(num_values-2)) then

       services_1d_interpolation = 0.0_DP

    else

       f1 = values(start_index-1)
       f2 = values(start_index)
       f3 = values(start_index+1)
       f4 = values(start_index+2)

       t0 = f2
       t1 = ((6.0_dp*f3)-(2.0_dp*f1)-(3.0_dp*f2)-f4)/6.0_dp
       t2 = (f1+f3-(2.0_dp*f2))/2.0_dp
       t3 = (f4-f1+(3.0_dp*(f2-f3)))/6.0_dp

       services_1d_interpolation = t0+off*(t1+off*(t2+off*t3))

    endif


  end function services_1d_interpolation



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  function services_polak_cg_coeff(prev_cov_dir,&
       cov_grad,contra_grad,prev_contra_grad, vec_size)

    !===============================================================!
    ! This subroutine calculates the conjugate gradient coefficient !
    ! according to the Polak formula:                               !
    !  b_{r+1}=g_{r+1}*(g_{r+1}-g_r)/( p_r * (g_{r+1}-g_r) )        !
    ! taking into account the covariant and contravariant           !
    ! character of the quantities involved.                         !
    !---------------------------------------------------------------!
    !                                                               !
    !                                                               !
    !---------------------------------------------------------------!
    ! Originaly written by Chris-Kriton Skylaris on 25/6/2001       !
    ! for the ONES code.                                            !
    ! Slightly modified by Arash A. Mostofi in March 2003.          !
    ! Parallelised by Chris-Kriton Skylaris on 30/12/2003.          !
    !===============================================================!

    use comms, only: pub_on_root, comms_reduce
    use constants, only: DP, stdout, NORMAL
    use rundat, only : pub_output_detail
    implicit none

    real(kind=DP)             :: services_polak_cg_coeff
    integer, intent(in)       :: vec_size
    real(kind=DP), intent(in) :: prev_cov_dir(vec_size)
    real(kind=DP), intent(in) :: cov_grad(vec_size)
    real(kind=DP), intent(in) :: contra_grad(vec_size)
    real(kind=DP), intent(in) :: prev_contra_grad(vec_size)

    real(kind=DP) :: denominator,eps


    eps =epsilon(1.0_DP)

    ! cks: calculate my node contribution
    denominator =sum( prev_cov_dir(1: vec_size) &
         *( contra_grad(1: vec_size) -prev_contra_grad(1: vec_size) )   )
    ! cks: add up contributions from all nodes
    call comms_reduce('SUM', denominator)


    if ( abs(denominator) .gt. eps ) then
       ! cks: contribution from my node
       services_polak_cg_coeff = &
            sum(cov_grad(1: vec_size) &
            *(contra_grad(1: vec_size) -prev_contra_grad(1: vec_size) ) ) / denominator
       ! cks: sum of contributions from all nodes
       call comms_reduce('SUM', services_polak_cg_coeff)
    else
       if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
            write(stdout,'(a/a)') &
            'WARNING: zero denominator in services_polak_cg_coeff', &
            '         Setting to zero'
       services_polak_cg_coeff=0.0_DP
    endif


    if ( abs(services_polak_cg_coeff).gt.(2.0_DP) ) then

       if (pub_on_root .and. pub_output_detail >= NORMAL) &
            write(stdout,'(a,f11.5,a)') &
            'WARNING: services_polak_cg_coeff too large &
            &(',services_polak_cg_coeff,'). Setting to zero'
       services_polak_cg_coeff=0.0_DP
    endif


  end function services_polak_cg_coeff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! LINE SEARCH SERVICES !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_cubic_fit_minimum( &
       min_length, min_value, success, &         ! output
       val1, val2, val3, grad1, coor2, coor3, &  ! input
       max_step)

    !==================================================================!
    ! Do line search by fitting cubic polynomial to f0, f1, f2 and g0. !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 4/12/2001 based on code      !
    ! written earlier by Peter D. Haynes.                              !
    ! Modified by Chris-Kriton Skylaris on 16/7/2003.                  !
    ! Return value to indicate success added by Nick Hine on 04/10/2008!
    !==================================================================!

    use comms, only : pub_on_root
    use constants, only: DP, stdout, NORMAL
    use rundat, only: pub_output_detail

    implicit none

    ! Arguments
    real(kind=DP), intent(out):: min_length, min_value
    logical, intent(out) :: success
    real(kind=DP), intent(in) :: val1, val2, val3, grad1, coor2, coor3
    real(kind=DP), intent(in) :: max_step
    ! coor1 is zero

    ! Local Variables
    real(kind=DP) :: aa, bb, xx, yy, disc, aon3b

    ! cks: initialise
    min_value =0.0_DP ; xx = coor2 ; yy = coor3

    aa = ((yy*yy*yy*val2 - xx*xx*xx*val3)/(yy-xx) &
         - (xx*xx+xx*yy+yy*yy)*val1 - xx*yy*(xx+yy)*grad1) / (xx*xx*yy*yy)

    bb = ((yy*yy*val2 - xx*xx*val3)/(xx-yy) + (xx+yy)*val1 &
         + xx*yy*grad1) / (xx*xx*yy*yy)

    disc=-1.0_DP

    if (abs(bb*yy/aa) > epsilon(1.0_DP)) then        ! avoid div by zero
       aon3b = aa / (3.0_DP * bb)
       disc = aon3b * aon3b - grad1 / (3.0_DP * bb)  ! discriminant

       if (disc >= 0.0_DP) then

          ! cks: set optimal_step to cubic minimum
          min_length = -aon3b + sign(sqrt(disc), bb)

          ! cks: Value of cubic polynomial at minimum
          min_value = val1 + min_length * (grad1 + min_length * &
               (aa + min_length * bb))

       end if
    end if

    ! cks: see if the cubic fit was successful and choose safe length if not.
    if (disc < 0.0_DP ) then
       if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
            write(stdout,'(a)') 'WARNING: Cubic fit was unsuccessful.'
       min_length=0.15_DP
       success = .false.
    ! cks: protection from crazy line search coefficients
    ! ndmh: now uses max_step input parameter
    else if (min_length > max_step) then
       if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
            write(stdout,'(a/a,f11.5,a)') &
            'WARNING in services_cubic_fit_minimum:','  cubic step (', &
            min_length,') too large - setting to safe value'
       min_length = 0.10_DP
       success = .false.
    else if (min_length < 0.0_DP) then
       if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
            write(stdout,'(a/a,f11.5,a)') &
            'WARNING in services_cubic_fit_minimum:','  cubic step (', &
            min_length,') less than zero - setting to safe value'
       min_length = 0.10_DP
       success = .false.
    else
       success = .true.
    endif

  end subroutine services_cubic_fit_minimum


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_line_search_parabola( &
       min_length, min_value, success, &                       ! output
       grad_at_zero, value_at_zero, value_at_trial_length, &   ! input
       trial_length,max_step)                                  ! input

    !==================================================================!
    ! Do line search by fitting a parabola.                            !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/6/2001.                   !
    ! Improved by Chris-Kriton Skylaris on 20/6/2003 and on 26/4/2004. !
    ! Return value to indicate success added by Nick Hine on 04/10/2008!
    !==================================================================!

    use comms, only : pub_on_root
    use constants, only: DP, stdout, NORMAL
    use rundat, only: pub_output_detail
    implicit none

    real(kind=DP), intent(out) :: min_length, min_value
    logical, intent(out) :: success
    real(kind=DP), intent(in) :: trial_length, grad_at_zero, value_at_zero &
         , value_at_trial_length,max_step

    ! cks: internal declarations
    real(kind =DP) :: linear_term
    real(kind =DP) :: quadratic_term
    real(kind =DP) :: pscale
    real(kind =DP) :: q_diff


    ! cks: determine pscale so that precision is maintained in polynomial fitting
    pscale =abs(value_at_trial_length -value_at_zero)
    if ( pscale > tiny(1.0_DP)) then
       pscale =1.0_DP/pscale
       if (value_at_trial_length > value_at_zero) then
          q_diff =1.0_DP
       else
          q_diff =-1.0_DP
       endif
    else
       pscale =100000.0_DP
       q_diff =pscale*value_at_trial_length -pscale*value_at_zero
    endif


    linear_term =  pscale*grad_at_zero  ! cks: b

    quadratic_term = q_diff -linear_term*trial_length


    if ( abs(  quadratic_term/ (pscale*value_at_zero)  ) > epsilon(1.0_DP)) then

       min_length = -0.5_DP *linear_term /quadratic_term ! -b/2a
       min_length = min_length *(trial_length**2)


       min_value = value_at_zero &
            -(0.25_DP * linear_term*linear_term/ (pscale*quadratic_term) ) &
            *(trial_length**2)  ! c -b^2/4a
       success = .true.
    else
       min_length = 0.1_DP
       success = .false.
    endif

    ! cks: protection from crazy line search coefficients
    ! ndmh: changed to use keyword for max line search step
    if (abs(min_length) > max_step) then
       if (pub_on_root .and. pub_output_detail >= NORMAL) &
            write(stdout,'(a/a,f11.5,a)') &
            'WARNING in services_line_search_parabola:', &
            '  quadratic step (', min_length, &
            ') too large - setting to safe value'
            write(stdout,'(a,f11.5,a)') &
            'value at zero =', value_at_zero
            write(stdout,'(a,f11.5,a)') &
            'linear_term=', linear_term
            write(stdout,'(a,f11.5,a)') &
            'quadratic_term=', quadratic_term
            write(stdout,'(a,f11.5,a)') &
            'trial_length=', trial_length
            write(stdout,'(a,f11.5,a)') &
            'pscale=', pscale
       min_length = 0.15_DP
       success = .false.
    end if


  end subroutine services_line_search_parabola


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_fit_parabola_ptsgrad(aa,bb,cc &
       ,grad_at_zero, value_at_zero, value_at_point, point)

    !=========================================================!
    ! Fits a parabola using its value and gradient at point 0 !
    ! and its value at point x.                               !
    !---------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 25/6/2001           !
    !=========================================================!

    use constants, only: DP
    implicit none

    real(kind=DP), intent(out) :: aa, bb, cc
    real(kind=DP), intent(in) :: grad_at_zero, value_at_zero &
         ,value_at_point, point

    ! cks: aa= ( f(0+x)-f(0)-f'(0)*x ) / ( x^2 )
    aa=(value_at_point - (value_at_zero + grad_at_zero*point) ) / ( point**2 )

    ! cks: bb= f'(0)
    bb=grad_at_zero

    ! cks: cc=f(0)
    cc=value_at_zero

  end subroutine services_fit_parabola_ptsgrad



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  function services_parabolic_step(Qinit, Q1, Q2, L1, L2)

    !==========================================================!
    ! Fits a parabola given value of function at three points. !
    !----------------------------------------------------------!
    ! Written by Arash A. Mostofi, December 2002.              !
    !==========================================================!

    use comms, only : pub_on_root
    use constants, only: DP, stdout
    implicit none

    real(kind=DP) :: services_parabolic_step
    real(kind=DP), intent(in) :: Qinit  ! initial value of function
    real(kind=DP), intent(in) :: Q1,Q2  ! value at trial steps 1 and 2
    real(kind=DP), intent(in) :: L1,L2  ! trial steps 1 and 2

    real(kind=DP) :: ll,aa,bb,cc,d1,d2,eps
    logical, parameter :: info=.false.

    eps = epsilon(1.0_DP)
    services_parabolic_step=10.0_DP
    cc=Qinit ; ll=L2/L1
    d1=Q1-Qinit ; d2 = Q2-Qinit

    bb = d1*ll - d2/ll ; aa = (d2 - d1*ll)/L2

    if (info .and. pub_on_root) then
       write(stdout,'(a,f22.15)') 'LINE SEARCH =====> quadratic coefficient =', aa
       write(stdout,'(a,f22.15)') 'LINE SEARCH =====> linear    coefficient =', bb
    endif

    if (abs(aa).gt.eps) then
       services_parabolic_step = -bb/(2*aa)
    else
       if (pub_on_root) write(stdout,'(a)') 'WARNING: Quadratic fit unsuccessful'
       services_parabolic_step = 0.5_DP
    endif

    !    if (abs(services_parabolic_step).gt.(1.0_DP)) then
    !       print*,'services_parabolic_step=' &
    !            ,services_parabolic_step
    !       print*,'setting services_parabolic_step equal to 0.07'
    !       services_parabolic_step=0.07_DP
    !    endif

  end function services_parabolic_step


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_polynomial_step(poly,poly_order,poly_num,Yi,Xi,Xj)

    !==========================================================!
    ! Fits a series of "poly_num" polynomial given values of   !
    ! "poly_num" functions at "poly_order" points.             !
    !----------------------------------------------------------!
    ! Written by Simon M.-M. Dubois, June 2011.                !
    !==========================================================!

    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout
    implicit none

    integer, intent(in) :: poly_order                    ! polynomial order + 1
    integer, intent(in) :: poly_num                      ! number of polynomials
    real(kind=DP), intent(in) :: Yi(poly_order,poly_num) ! set of ordinates
    real(kind=DP), intent(in) :: Xi(poly_order)          ! set of abcissa
    real(kind=DP), intent(in) :: Xj                      ! interpolated abcissa
    real(kind=DP), intent(out) :: poly(poly_num)

    real(kind=DP) :: Vdm_mat(poly_order,poly_order)
    real(kind=DP) :: Vdm_vec(poly_order)
    !real(kind=DP) :: Yj(poly_num,1)
    integer :: ipiv(poly_order)
    integer :: info, io, jo


    ! Compose Vandermonde matrix
    do io = 1, poly_order
       do jo = 1, poly_order
          Vdm_mat(jo,io) = Xi(jo)**(poly_order-io)
       enddo
       Vdm_vec(io) = Xj**(poly_order-io)
    enddo 

    ! Solve Vandermonde system of equations (i.e. compute polynomial coeffs)
    call dgesv(poly_order,poly_num,Vdm_mat,poly_order,ipiv,Yi,poly_order,info)
    if (info .ne. 0 ) then
       write(stdout,'(a,i4,a)') 'Error in services_polynomial_step : &
           &computation of the polynomial coefficients with &
           &lapack_dgesv failed with error, ', info, ' !'
       call comms_abort
    endif
 
    poly(:) = 0.0_dp
    do io = 1, poly_num
       do jo = 1, poly_order
          poly(io) = Yi(jo,io)*Vdm_vec(jo)
       enddo
    enddo
    !! Compute the interpolated ordinates
    !call dgemm('T','N',poly_num,1,poly_order,1.0_dp,Yi,poly_order,&
    !       Vdm_vec,1,0.0_dp,Yj,poly_num)

    !poly(:) = Yj(:,1) 

  end subroutine services_polynomial_step


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function services_cubic_minimum(aa,bb,cc,dd)


    !==================================================================!
    ! Returns the x value at which the cubic polynomial                !
    ! f(x)=aa*x^3+bb*x^2+cc*x+dd has a minimium. If the polynomial has !
    ! no minimum the function stops.                                   !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 19/4/2001.                   !
    !==================================================================!


    use comms, only : pub_on_root
    use constants, only: DP, stdout
    implicit none

    real(kind=DP) :: services_cubic_minimum
    real(kind=DP), intent(in) :: aa, bb, cc, dd

    ! cks: internal declarations
    real(kind=DP) :: x1, x2, fx1, fx2, discriminant

    discriminant=4.0_DP*(bb**2-3.0_DP*aa*cc)
    if (discriminant.lt.(0.0_DP)) then
       if (pub_on_root) write(stdout,'(a)')&
            'negative discriminant in services_cubic_minimum'
       !       print*,'ONES execution stops'
       !       stop
       if (pub_on_root) write(stdout,'(a)')&
            'RETURNING services_cubic_minimum=0.4_DP'
       services_cubic_minimum=0.4_DP
       return
    endif
    discriminant=sqrt(discriminant)

    if (aa.ne.(0.0_DP)) then
       x1=(-2.0_DP*bb+discriminant)/(6.0_DP*aa)
       x2=(-2.0_DP*bb-discriminant)/(6.0_DP*aa)
    else if (bb.ne.(0.0_DP)) then
       x1=-cc/(2.0_DP*bb)
       x2=x1
    else
       x1=0.0_DP
       x2=x1
    endif

    fx1=aa*x1**3+bb*x1**2+cc*x1+dd
    fx2=aa*x2**3+bb*x2**2+cc*x2+dd

    if (fx1.lt.fx2) then
       services_cubic_minimum=x1
    else
       services_cubic_minimum=x2
    endif

  end function services_cubic_minimum


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!gibo-start


!  subroutine services_cubic_fit_maximum( &
!       min_length, min_value, success, &         ! output
!       val1, val2, val3, grad1, coor2, coor3, &  ! input
!       max_step)

  subroutine services_cubic_fit_maximum( &
       max_length, max_value, success, &         ! output
       y0, y1, y2, g0, x1, x2, &  ! input
       max_step)

    !==================================================================!
    ! Do line search by fitting cubic polynomial to y0, y1, y2 and g0. !
    !          y(x) = ax^3 + bx^2 + cx + d                             !
    !      dy(x)/dx = 3ax^2 + 2b*x + c                                 !
    ! d^2y(x)/da^2x = 6ax + 2b                                         !
    !                                                                  !
    ! Mind that for computational convenience the x-frame is shifted   !
    ! to have x0= 0, yielding                                          !
    ! *** y(x0) = d                                                    !
    ! *** dy(x0)/dx = 3a(x0)^2 +2b(x0) +c = c                          !
    !------------------------------------------------------------------!
    ! Adapted from services_cubic_fit_minimum                          !
    ! by Gilberto Teobaldi on 9/12/11                                  !
    !==================================================================!

    use comms, only : pub_on_root
    use constants, only: DP, stdout, NORMAL
    use rundat, only: pub_output_detail

    implicit none

    ! Arguments
    !real(kind=DP), intent(out):: min_length, min_value
    real(kind=DP), intent(out):: max_length, max_value
    logical, intent(out) :: success
    !real(kind=DP), intent(in) :: val1, val2, val3, grad1, coor2, coor3
    real(kind=DP), intent(in) :: y0, y1, y2, g0, x1, x2
    real(kind=DP), intent(in) :: max_step
    ! coor1 is zero

    ! Local Variables
    !real(kind=DP) :: aa, bb, xx, yy, disc, aon3b
    real(kind=DP) :: a, b, c, d
    !real(kind=DP) :: x0, x1, x2
    real(kind=DP) :: x1_2, x1_3, x2_2
    real(kind=DP) :: x_stat_1, x_stat_2
    real(kind=DP) :: second_deriv_1, second_deriv_2
    !real(kind=DP) :: y0, y1, y2
    !real(kind=DP) :: g0
    real(kind=DP) :: den, num
    real(kind=DP) :: disc
    !real(kind=DP) :: max_value, max_length

  !initialise a negative discriminant for dy/dx= 3ax^2 +2bx + c = 0
  !disc = -1.0_DP
  disc = -1.0

  !====STEP-1:
  !fit cubic polynomial (y = a*x^3 + b*x^2 + c*x + d)
  !to (available) y(x0), y(x1), y(x2) and g0= [dy/dx]@x0 
  !MIND: from computational convenience the x-frame is shifted
  !to have x0=0., therefore y(x0)=d

  x1_2 = x1*x1
  x1_3 = x1_2*x1

  x2_2 = x2*x2

  !den = x2_2*x1 - 1.0_DP
  den = x2_2*x1 - 1.0

  num = y2*x1_3 - y1 + y0 + g0*x1 -g0*x2*x1_3 - y0*x1_3

  NUM_STAB: if (abs(den) > epsilon(1.0_DP)) then  ! avoid division by zero
  !NUM_STAB: if (abs(den) > 0.) then  ! avoid division by zero

    a = (y1 - y0 - g0*x1 - num/den) / x1_2
    b = num /(x1_2*den)
    c = g0
    d = y0

    !====STEP-2 
    ! Calculate the x-coordinate of the stationary points
    ! dy/dx= 3ax^2 +2bx + c = 0
    ! and check sign of 2nd derivative (6ax+2b) at stat.points
    ! [we are after a maximu i.e. d^2y/dx^2 < 0.

    ! discriminant of dy/dx= 3ax^2 +2bx + c = 0
    disc = 4.0_DP*b*b -12.0_DP*a*c
    !disc = 4.0*b*b -12.0*a*c

    if (disc >= 0.0_DP) then
    !if (disc >= 0.) then

        x_stat_1 = (-2.0_DP*b - disc)/( 6.0_DP*a)
        !x_stat_1 = (-2.0*b - disc)/( 6.0*a)
  
        second_deriv_1 = 6.0_DP*a*x_stat_1 +2.0_DP*b
        !second_deriv_1 = 6.0*a*x_stat_1 +2.0*b
  
        x_stat_2 = (-2.0_DP*b + disc)/( 6.0_DP*a)
        !x_stat_2 = (-2.0*b + disc)/( 6.0*a)
  
        second_deriv_2 = 6.0_DP*a*x_stat_2 +2.0_DP*b
        !second_deriv_2 = 6.0*a*x_stat_2 +2.0*b
  
        if (second_deriv_1 < 0.0_DP) then
        !if (second_deriv_1 < 0.0) then
  
           ! set optimum step (from x0) to reach cubic maximum
           max_length = x_stat_1
        
           ! value of cubic polynomial at cubic maximum
           max_value = a*x_stat_1*x_stat_1*x_stat_1 + b*x_stat_1*x_stat_1 &
                     + c*x_stat_1 + d
  
        elseif (second_deriv_2 > 0.0_DP) then
        !elseif (second_deriv_2 > 0.0) then
  
           ! set optimum step (from x0) to reach cubic maximum
           max_length = x_stat_2
        
           ! value of cubic polynomial at cubic maximum
           max_value = a*x_stat_2*x_stat_2*x_stat_2 + b*x_stat_2*x_stat_2 &
                     + c*x_stat_2 + d
        else
  
         if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
               write(stdout,'(a/a)') &
               'WARNING in services_cubic_fit_maximum:',&
               &' none of the (two) stationary point is a maximum!!!'
  
        endif

    endif

  endif NUM_STAB

  ! cks: see if the cubic fit was successful and choose safe length if not.
  if (disc < 0.0_DP ) then
     if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
          write(stdout,'(a)') 'WARNING: Cubic fit was unsuccessful.'
     !min_length=0.15_DP
     max_length=0.15_DP
     success = .false.
  ! cks: protection from crazy line search coefficients
  ! ndmh: now uses max_step input parameter
  !else if (min_length > max_step) then
  else if (max_length > max_step) then
     if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
          write(stdout,'(a/a,f11.5,a)') &
          'WARNING in services_cubic_fit_minimum:','  cubic step (', &
          !min_length,') too large - setting to safe value'
          max_length,') too large - setting to safe value'
     !min_length = 0.10_DP
     max_length = 0.10_DP
     success = .false.
  !else if (min_length < 0.0_DP) then
  else if (max_length < 0.0_DP) then
     if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
          write(stdout,'(a/a,f11.5,a)') &
          'WARNING in services_cubic_fit_minimum:','  cubic step (', &
          !min_length,') less than zero - setting to safe value'
          max_length,') less than zero - setting to safe value'
     !min_length = 0.10_DP
     max_length = 0.10_DP
     success = .false.
  else
     success = .true.
  endif

  end subroutine services_cubic_fit_maximum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gibo-end



  subroutine services_flush(unit)

    !==================================================================!
    ! Flushes output (stdout by default)                               !
    !------------------------------------------------------------------!
    ! Written by Peter Haynes 12 November 2004                         !
    !==================================================================!

    use comms, only : pub_on_root
    use constants, only: stdout

    implicit none

    ! Argument
    integer, optional, intent(in) :: unit

    ! Local variables
    integer :: unit_local     ! Local copy of argument
#ifdef SGI
    integer :: ierr           ! Error flag
#endif
#ifdef SUN
    integer :: ierr           ! Error flag
    integer :: flush
#endif
#ifndef INTFLUSH
    external flush
#endif

    if (present(unit)) then
       unit_local = unit
    else
       unit_local = stdout
    end if

    if (pub_on_root) then

#ifdef SUN
       ierr = flush(unit_local)
       if (ierr /= 0) write(stdout,'(a,i6)') 'WARNING in services_flush: &
            &flush failed with code ',ierr

#else
#ifdef SGI
       ! SGI introduced a second argument starting from the v7.4
       ! of their runtime environment. This form is still acceptable
       ! by the older environments - the return value is ignored in those
       ! older cases, and the code exits on error.
       ! This second argument is also allowed under UNICOS as an
       ! optional argument.
       call flush(unit_local,ierr)
       if (ierr /= 0) write(stdout,'(a,i6)') 'WARNING in services_flush: &
            &flush failed with code ',ierr

#else
!CW
#ifndef FLUSH_EXT_NAME
          call flush(unit_local)
#else
          call flush_(unit_local)
#endif
!END CW
#endif
#endif

    end if

  end subroutine services_flush


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_radial_transform(l,power,nrpts,rlog,rab,nqpts,qmax, &
       rfunc,qfunc)

    !=========================================================================!
    ! Transforms a radial real space function on a logarithmic grid to a      !
    ! regular reciprocal space grid                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  nrpts  (input)  : number of real space points                          !
    !  rlog   (input)  : the logarithmic grid                                 !
    !  rab    (input)  : the "spacings" of the logarithmic grid               !
    !  nqpts  (input)  : number of points in the regular reciprocal grid      !
    !  qmax   (input)  : the maximum q-vector for the transform               !
    !  rfunc  (input)  : the real space radial function                       !
    !  qfunc  (output) : the reciprocal space function                        !
    !-------------------------------------------------------------------------!
    ! Adapted for ONETEP by Nicholas Hine from a similar routine in CASTEP,   !
    ! which was originally written by Chris J. Pickard on 08/06/2000.         !
    !=========================================================================!

    use constants, only: stdout
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer,intent(in) :: nrpts
    integer,intent(in) :: l
    integer,intent(in) :: power
    real(kind=dp),intent(in) :: rlog(nrpts)
    real(kind=dp),intent(in) :: rab(nrpts)
    integer,intent(in) :: nqpts
    real(kind=dp),intent(in) :: qmax
    real(kind=dp),intent(in) :: rfunc(nrpts)
    real(kind=dp),intent(out) :: qfunc(nqpts)

    ! Local Variables
    integer :: ierr        ! Error flag
    integer :: nq,nr       ! Real and recip point counters
    real(kind=dp) :: q,qr        ! q-vector and product q.r
    real(kind=dp),allocatable :: rwgt(:)     ! Integration weights
    real(kind=dp),allocatable :: lookup(:,:) ! Bessel function lookup table

    ! Allocate temporary variables
    allocate(rwgt(nrpts),stat=ierr)
    call utils_alloc_check('services_radial_transform','rwgt',ierr)
    allocate(lookup(nqpts,nrpts),stat=ierr)
    call utils_alloc_check('services_radial_transform','lookup',ierr)

    ! Check arguments
    if (modulo(nrpts,2)==0) call utils_abort('Error in &
         &services_radial_transform: number of grid points must be odd')

    ! Integration weights
    rwgt(1) = rab(1)*4.0_dp/12.0_dp*rlog(1)**power
    do nr=2,nrpts-1,2
       rwgt(nr)=rab(nr)*16.0_dp/12.0_dp*rlog(nr)**power
    end do
    do nr=3,nrpts-2,2
       rwgt(nr)=rab(nr)*8.0_dp/12.0_dp*rlog(nr)**power
    end do
    rwgt(nrpts) = rab(nr)*4.0_dp/12.0_dp*rlog(nrpts)**power

    ! Loop over the radial reciprocal space grid
    do nq=1,nqpts
       q = real(nq-1,dp)*qmax/real(nqpts-1,dp)
       lookup(nq,1) = services_sbessj(l,0.0_dp)*rwgt(1)
       do nr=2,nrpts
          qr = q*rlog(nr)
          lookup(nq,nr) = services_sbessj(l,qr)*rwgt(nr)
       end do
    end do

    ! Do the matrix-vector multiplication
    call dgemv('N',nqpts,nrpts,1.0_dp,lookup(1,1),nqpts,rfunc,1,0.0_dp, &
         qfunc,1)

    ! Deallocate temporary variables
    deallocate(lookup,stat=ierr)
    call utils_dealloc_check('services_radial_transform','lookup',ierr)
    deallocate(rwgt,stat=ierr)
    call utils_dealloc_check('services_radial_transform','rwgt',ierr)

  end subroutine services_radial_transform


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_regular_transform(l,power,nrpts,rmax,nqpts,qmax, &
       rfunc,qfunc)

    !=========================================================================!
    ! Transforms a radial real space function on a logarithmic grid to a      !
    ! regular reciprocal space grid                                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  nrpts  (input)  : number of real space points                          !
    !  rmax   (input)  : the maximum r-point of the regular real grid         !
    !  nqpts  (input)  : number of points in the regular reciprocal grid      !
    !  qmax   (input)  : the maximum q-vector for the transform               !
    !  rfunc  (input)  : the real space radial function                       !
    !  qfunc  (output) : the reciprocal space function                        !
    !-------------------------------------------------------------------------!
    ! Adapted for ONETEP by Nicholas Hine from a similar routine in CASTEP,   !
    ! which was originally written by Chris J. Pickard on 08/06/2000.         !
    !=========================================================================!

    use constants, only: stdout
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    integer,intent(in) :: l
    integer,intent(in) :: power
    integer,intent(in) :: nrpts
    real(kind=dp),intent(in) :: rmax
    integer,intent(in) :: nqpts
    real(kind=dp),intent(in) :: qmax
    real(kind=dp),intent(in) :: rfunc(nrpts)
    real(kind=dp),intent(out) :: qfunc(nqpts)

    ! Local Variables
    integer :: iq,ir         ! Real and recip point counters
    real(kind=dp) :: q,r     ! q-vector and r-point
    real(kind=dp) :: dr
    real(kind=dp), allocatable :: work(:)
    
    allocate(work(nrpts))
    dr = rmax/real(nrpts-1,dp)

    ! Check arguments
    if (modulo(nrpts,2)==0) call utils_abort('Error in &
         &services_radial_transform: number of grid points must be odd')

    ! Loop over the radial reciprocal space grid
    do iq=1,nqpts
       q = real(iq-1,dp)*qmax/real(nqpts-1,dp)
       do ir=1,nrpts
          r = real(ir-1,dp)*rmax/real(nrpts-1,dp)
          work(ir) = rfunc(ir)*r**power*services_sbessj(l,q*r)
       end do
       qfunc(iq) = work(1) + work(nrpts)
       do ir=2,nrpts-1,2
          qfunc(iq) = qfunc(iq) + 4.0_DP*work(ir) + 2.0_DP*work(ir+1)
       end do
       qfunc(iq) = qfunc(iq) * dr / 3.0_DP
    end do

    deallocate(work)

  end subroutine services_regular_transform


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function services_radial_integral(npts,rab,f,inter)

    !=========================================================================!
    ! Simpson's rule integrator for a function on a radial mesh.              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  npts (input)   : number of points in radial mesh                       !
    !  rab (input)    : "spacings" of the radial grid                         !
    !  func (input)   : function to be integrated                             !
    !  inter (output) : intermediate results (optional)                       !
    !-------------------------------------------------------------------------!
    ! Adapted for ONETEP by Nicholas Hine from a similar routine in CASTEP,   !
    ! which was originally written by Chris J. Pickard on 09/06/2000.         !
    !=========================================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP), intent(in) :: rab(:)
    real(kind=DP), intent(in) :: f(:)
    real(kind=DP), optional, intent(out) :: inter(:)

    ! Local Variables
    real(kind=DP), parameter :: c1o12 = 1.0_DP/12.0_DP
    real(kind=DP), parameter :: c5o12 = 5.0_DP/12.0_DP
    real(kind=DP), parameter :: c8o12 = 8.0_DP/12.0_DP
    real(kind=DP) :: fm1,f0,fp1
    integer :: ipt

    ! Check arguments
    if (modulo(npts,2)==0) call utils_abort('Error in &
         &services_radial_integral: number of grid points must be odd')

    ! Initialise
    services_radial_integral = 0.0_DP
    if (present(inter)) inter(1) = 0.0_DP

    ! Loop over alternate grid points, adding up integral and calculating
    ! intermediate values if required
    do ipt=2,npts-1,2

       fm1 = f(ipt-1)*rab(ipt-1)
       f0  = f(ipt  )*rab(ipt  )
       fp1 = f(ipt+1)*rab(ipt+1)

       services_radial_integral = services_radial_integral &
            + c5o12*fm1 + c8o12*f0 - c1o12*fp1
       if (present(inter)) inter(ipt) = services_radial_integral

       services_radial_integral = services_radial_integral &
            - c1o12*fm1 + c8o12*f0 + c5o12*fp1
       if (present(inter)) inter(ipt+1) = services_radial_integral

    end do

  end function services_radial_integral


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function services_radial_integral_rmax(npts,rab,r,rmax,f,inter)

    !=========================================================================!
    ! Wrapper for services_radial_integral to integrate up to a fixed rmax    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  npts (input)   : number of points in radial grid                       !
    !  rab (input)    : "spacings" of the radial grid                         !
    !  r   (input)    : radial grid position                                  !
    !  rmax (input)   : upper limit of integral                               !
    !  func (input)   : function to be integrated                             !
    !  inter (output) : array to hold intermediate results                    !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 02/06/2010.                                 !
    !=========================================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP), intent(in) :: rab(:)
    real(kind=DP), intent(in) :: r(:)
    real(kind=DP), intent(in) :: rmax
    real(kind=DP), intent(in) :: f(:)
    real(kind=DP), intent(out) :: inter(:)

    ! Local Variables
    integer :: ir
    integer :: npts_int
    real(kind=DP) :: int_full

    ! Calculate the interpolation point
    ir = services_locate_interp(rmax,r,npts)
    ! Find endpoint: min of a) position of rmax + 4 pts, or b) end of array
    npts_int = min(npts,ir+4)
    npts_int = npts_int + modulo(npts_int,2) - 1
    ! Calculate the integral up to npts_int
    int_full = services_radial_integral(npts_int,rab,f,inter)
    ! Do the interpolation
    services_radial_integral_rmax = services_linear_interpolation(rmax, &
         inter(ir),inter(ir+1),r(ir),r(ir+1))

  end function services_radial_integral_rmax


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real(kind=DP) function services_regular_integral(npts,dr,f)

    !=========================================================================!
    ! Integrates a function on a regular grid.                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  npts (input)   : number of points in grid                              !
    !  dr (input)    :  spacing of the grid                                   !
    !  r   (input)    : radial grid position                                  !
    !  func (input)   : function to be integrated                             !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 09/09/2010.                                 !
    !=========================================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP), intent(in) :: dr
    real(kind=DP), intent(in) :: f(:)

    ! Local Variables
    real(kind=DP) :: a(4)
    real(kind=DP) :: intf

    a(1)=17.0_DP/48.0_DP
    a(2)=59.0_DP/48.0_DP
    a(3)=43.0_DP/48.0_DP
    a(4)=49.0_DP/48.0_DP

    intf = a(1)*(f(1)+f(npts)) + &
         a(2)*(f(2)+f(npts-1)) + &
         a(3)*(f(3)+f(npts-2)) + &
         a(4)*(f(4)+f(npts-3))
    intf = intf + sum(f(5:npts-4))
    services_regular_integral = intf*dr

  end function services_regular_integral


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_radial_derivative(grad,func,npts,xmax)

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP),intent(in) :: xmax
    real(kind=DP),intent(in) :: func(:)
    real(kind=DP),intent(out) :: grad(:)

    ! Local Variables
    integer :: ipt
    real(kind=DP) :: t1,t2,t3

    ! Gradient is zero at origin by symmetry
    grad(1) = 0.0_DP
    t1 = 0.0_DP; t2 = 0.0_DP; t3 = 0.0_DP

    ! Cubic interpolation
    do ipt=2,npts-2
       t1 = (6.0_DP*func(ipt+1) - 2.0_DP*func(ipt-1) - 3.0_DP*func(ipt) &
            - 1.0_DP*func(ipt+2)) / 6.0_DP
       t2 = (1.0_DP*func(ipt-1) + 1.0_DP*func(ipt+1) - 2.0_DP*func(ipt)) &
            / 2.0_DP
       t3 = (1.0_DP*func(ipt+2) - 1.0_DP*func(ipt-1) + 3.0_DP*func(ipt) &
            - 3.0_DP*func(ipt+1)) / 6.0_DP
       grad(ipt) = t1
    end do

    ! Last two points
    grad(npts-1) = t1 + 2.0_DP*t2 + 3.0_DP*t3
    grad(npts) = t1 + 4.0_DP*t2 + 12.0_DP*t3

    ! Normalise for dr
    grad(:) = grad(:)*real(npts-1,kind=DP)/xmax

  end subroutine services_radial_derivative


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine services_analytic_limit(npts,rad,func)

    !=========================================================================!
    ! Finds the r->0 limit of a radial function through a spline fit to the   !
    ! first few points on the radial grid.                                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  npts (input)   : number of points in radial grid                       !
    !  rad (input)    : radial grid position                                  !
    !  func (input)   : function to be integrated                             !
    !-------------------------------------------------------------------------!
    ! Written by Nicholas Hine on 20/07/2010.                                 !
    !=========================================================================!

    use utils, only: utils_abort

    implicit none

    ! Arguments
    integer, intent(in) :: npts
    real(kind=DP), intent(in) :: rad(:)
    real(kind=DP), intent(inout) :: func(:)

    if ((size(rad)<npts).or.(size(func)<npts)) then
       call utils_abort('Error in services_analytic_limit: mismatching array &
            & sizes')
    end if

    if (npts>3) then
       func(1) = func(4)+3*(func(2)-func(3))
    else
       func(1) = func(2)
    end if

  end subroutine services_analytic_limit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine services_rationalise_coords(nr,r_abs)

    !=============================================================!
    ! This subroutine rationalises cartesian coordinates by       !
    ! translating them back into the home simulation cell.        !
    !-------------------------------------------------------------!
    ! Originally, md_update_coordinates subroutine written by     !
    ! Arash A. Mostofi in 2006. Modified and displaced by         ! 
    ! Simon M.-M. Dubois in Oct. 2010.                            !
    !=============================================================!

    use constants,       only: dp,two_pi,stdout
    use ion,             only: element
    use simulation_cell, only: pub_cell
    use utils,           only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments 
    integer, intent(in) :: nr
    real(kind=dp), intent(inout) :: r_abs(3,nr)

    ! Local variables
    integer :: j,iat
    integer :: rshift
    real(kind=dp) :: r_frac(3,nr)

    ! convert Cartesian coordinates to fractional and rationalise
    do iat=1,nr
       r_frac(1,iat) = pub_cell%b1%x*r_abs(1,iat) + pub_cell%b1%y*r_abs(2,iat) + pub_cell%b1%z*r_abs(3,iat)
       r_frac(2,iat) = pub_cell%b2%x*r_abs(1,iat) + pub_cell%b2%y*r_abs(2,iat) + pub_cell%b2%z*r_abs(3,iat)
       r_frac(3,iat) = pub_cell%b3%x*r_abs(1,iat) + pub_cell%b3%y*r_abs(2,iat) + pub_cell%b3%z*r_abs(3,iat)
       r_frac(:,iat) = r_frac(:,iat) / two_pi
       do j=1,3
          if (r_frac(j,iat).lt.0.0_dp .or. r_frac(j,iat).ge.1.0_dp) then 
             rshift = floor(r_frac(j,iat))
             r_frac(j,iat) = r_frac(j,iat) - real(rshift)
          endif
       enddo
    enddo

    ! convert rationalised fractional coordinates back to Cartesian
    do iat=1,nr
       r_abs(1,iat) = pub_cell%a1%x*r_frac(1,iat) + pub_cell%a2%x*r_frac(2,iat) + pub_cell%a3%x*r_frac(3,iat)
       r_abs(2,iat) = pub_cell%a1%y*r_frac(1,iat) + pub_cell%a2%y*r_frac(2,iat) + pub_cell%a3%y*r_frac(3,iat)
       r_abs(3,iat) = pub_cell%a1%z*r_frac(1,iat) + pub_cell%a2%z*r_frac(2,iat) + pub_cell%a3%z*r_frac(3,iat)
    enddo

    return

  end subroutine services_rationalise_coords


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function services_maxboltzdist() 

    !=============================================================!
    ! This function generates a pseudo-random number normally     !
    ! distributed (zero expectation, unit variance)               !
    !-------------------------------------------------------------!
    ! Originally written by Arash A. Mostofi as md_maxboltzdist   !
    ! Displaced by Simon M.-M. Dubois in July 2011.               !
    !=============================================================!

      use constants, only: stdout

      implicit none

      real(kind=dp) :: services_maxboltzdist
      real(kind=dp) :: ran(2),x(2),w,f

      do 
         call random_number(ran)
         x(:) = 2.0_dp * ran(:) - 1.0_dp
         w = x(1)*x(1) + x(2)*x(2)
         if (w.lt.1.0_dp) exit
      enddo

      f = sqrt((-2.0_dp*log(w))/w)
      services_maxboltzdist = x(1) * f

      return

    end function services_maxboltzdist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine services_rms_fit(alpha,nr,ncol,nrow,r_in,r_out)

    !========================================================================!
    ! Computes the coefficients alpha such that the root mean square         !
    ! deviation between Rout and R = \Sum_{i=1}^{nr} (alpha_i * Rin_i)       !
    ! is minimized.                                                          !
    !                                                                        !
    ! Arguments :                                                            !
    ! alpha	(output) : computed coefficients                             !   
    ! nr        (input)  : the number of elements in the basis set           !
    ! nrow      (input)  : the first dimension of Rin and Rout               !
    ! ncol      (input)  : the second dimension of Rin and Rout              !
    ! Rin       (input)  : basis set [dim(nrow,ncol,nr)]                     !
    ! Rout      (input)  : reference vector [dim(nrow,ncol)]                 !
    !------------------------------------------------------------------------!
    ! Written by Simon M.-M. Dubois 06/04/2011                               !
    !========================================================================!

    use comms, only : pub_on_root, comms_abort
    use constants, only: DP, stdout

    implicit none

    ! Argument
    integer, intent(in)          :: nr     ! The number of input elements
    integer, intent(in)          :: ncol   !  
    integer, intent(in)          :: nrow                ! # of row in each dataset
    real(kind=DP), intent(in)    :: r_in(nrow,ncol,nr)  ! datasets
    real(kind=DP), intent(in)    :: r_out(nrow,ncol)    ! 
    real(kind=DP), intent(out)   :: alpha(nr)

    ! Internal variables
    integer                   :: istep, jstep
    integer                   :: ic, ir, il
    integer                   :: length
    real(kind=DP)             :: dr(nrow*ncol,nr)
    real(kind=DP)             :: drdr_ij(nr-1,nr-1)
    real(kind=DP)             :: drdr_in(nr-1)
    integer                   :: ipiv(nr-1)
    integer                   :: ierr

    !======================================================================!

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
       'DEBUG: Entering services_rms_fit'
#endif

    ! Compose variations
    dr(:,:) = 0.0d0
    do istep = 1, nr-1
       do ic = 1, ncol
          do ir = 1, nrow
             dr((ic-1)*ncol+ir,istep) = r_in(ir,ic,istep) - r_in(ir,ic,istep+1)
          enddo
       enddo
    enddo
    do ic = 1, ncol
       do ir = 1, nrow
          dr((ic-1)*ncol+ir,nr) = r_out(ir,ic) - r_in(ir,ic,1)
       enddo
    enddo
    length = nrow*ncol

    ! Compose the linear matrix
    drdr_ij(:,:) = 0.0d0
    drdr_in(:) = 0.0d0
    do istep = 1, nr-1
       do jstep = 1, nr-1
          do il = 1, length
             drdr_ij(istep,jstep) = drdr_ij(istep,jstep) + &
                        dr(il,istep)*dr(il,jstep)
          enddo
       enddo
    enddo
    do istep = 1, nr-1
       do il = 1, length
          drdr_in(istep) = drdr_in(istep) + &
                     dr(il,istep)*dr(il,nr)
       enddo
    enddo

    ! Solve the system of linear equations     
    call dgesv(nr-1,1,drdr_ij,nr-1,ipiv,drdr_in,nr-1,ierr)
    if (ierr .ne. 0 ) then
       write(stdout,'(a,i4,a)') 'Error in services_rms_fit : &
           &computation of the extrapolation coefficients with &
           &lapack_dgesv failed with error, ', ierr, ' !'
    endif

    alpha(1) = 1 + drdr_in(1)
    do istep = 2, nr-1
       alpha(istep) = drdr_in(istep) - drdr_in(istep-1)
    enddo
    alpha(nr) = -drdr_in(nr-1)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') &
       'DEBUG: Leave services_rms_fit'
#endif

    return

  end subroutine services_rms_fit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module services
