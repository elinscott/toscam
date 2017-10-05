! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*- !
!================================================================!
!                                                                !
!                         Phonons Module                         !
!                                                                !
! This module calculates the Gamma-point vibrational frequencies !
! of the cell using the finite displacement method.              !
!----------------------------------------------------------------!
! Written by Fabiano Corsetti (May 2011)                         !
!================================================================!

module phonon

  use constants, only: DP

  implicit none

  private

  ! Conversion constants
  real(kind=DP), parameter :: electron_mass_u =5.4857990943E-4_DP
  real(kind=DP), parameter :: au2THz = 1.0_DP/2.41888468d-5
  real(kind=DP), parameter :: au2inv_cm = au2THz/2.99792458d-2
  ! Boltzmann constant in Eh/K
  real(kind=DP), parameter :: k_B = 3.1668115744561575d-06

  public :: phonon_main

contains

  subroutine phonon_main(total_energy,forces,elements)

    !==================================================================!
    ! This subroutine calculates and outputs the phonon frequencies of !
    ! a system using a finite-displacement method.                     !
    !------------------------------------------------------------------!
    ! Arguments:                                                       !
    !  total_energy (inout)  : Total energy of the system.             !
    !  forces (inout)        : Total force on each ion.                !
    !  elements (inout)      : List of ions stored in ELEMENT type.    !
    !------------------------------------------------------------------!
    ! Written by Fabiano Corsetti in June 2011.                        !
    ! Reformatted by Nicholas Hine in October 2011.                    !
    !==================================================================!

    use dense, only: DEM, dense_create, dense_destroy, dense_put_element, &
         dense_get_element, dense_eigensolve, dense_init, dense_exit
    use constants, only: DP,TWO_PI,stdout,periodic_table_mass
    use energy_and_force, only: energy_and_force_calculate
    use ion, only: element
    use comms, only: comms_abort, comms_bcast,pub_on_root, pub_root_node_id
    use simulation_cell, only: pub_cell
    use utils, only: utils_unit, utils_alloc_check, utils_dealloc_check
    use rundat, only: pub_rootname, write_denskern, write_tightbox_ngwfs, &
         read_denskern, read_tightbox_ngwfs, pub_phonon_disp, &
         pub_phonon_fmax, pub_phonon_farming_task, pub_have_phonon_disp_list, &
         pub_num_disp, pub_phonon_disp_list, pub_phonon_tmin, pub_phonon_tmax, &
         pub_phonon_deltat, pub_phonon_min_freq, print_qc, &
         pub_geom_reuse_dk_ngwfs, maxit_ngwf_cg

    implicit none

    ! Arguments
    type(ELEMENT), intent(inout) :: elements(pub_cell%nat)
    real(kind=DP), intent(inout) :: total_energy
    real(kind=DP), intent(inout) :: forces(1:3,pub_cell%nat)

    ! Local Variables
    type(DEM) :: dynamical_dem
    type(DEM) :: dynamical_symm_dem
    type(DEM) :: eigenvecs_dem
    type(DEM) :: unit_mat_dem
    character(len=1) :: dir_char
    character(len=1) :: cart_char
    character(len=256) :: output_file
    character(len=256) :: disp_number
    logical :: do_disp
    logical :: converged
    logical :: print_warning
    integer :: i, j, k, l, T_num_iter
    integer :: output_unit
    integer :: ierr
    real(kind=DP) :: Fmax
    real(kind=DP) :: eql_pos
    real(kind=DP) :: el
    real(kind=DP) :: el_T
    real(kind=DP) :: el_symm
    real(kind=DP) :: freq_conv
    real(kind=DP) :: zero_point_E, T, F, S, U, C_v, beta, sqrt_freq, exp_freq
    real(kind=DP), allocatable, dimension(:)   :: freqs
    real(kind=DP), allocatable, dimension(:,:) :: force_consts

    ! Display Banner
    if (pub_on_root) then
       write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
            &-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
       write(stdout,'(a)') 'Starting phonon calculation...'
       write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
            &-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
       write(stdout,'()')
    end if

    ! Check inputs
    if ((pub_phonon_farming_task<0) .or. (pub_phonon_farming_task>3)) then
       if (pub_on_root) then
          write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-&
               &+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
          write(stdout,'(a)') 'ERROR: invalid value of phonon_farming_task'
          write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-&
               &+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
       end if
       call comms_abort
    end if

    if (pub_have_phonon_disp_list) then
       if ((any(pub_phonon_disp_list(1:pub_num_disp)<0)) .or. &
            (any(pub_phonon_disp_list(1:pub_num_disp)>3*pub_cell%nat))) then
          if (pub_on_root) then
             write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
                  &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
             write(stdout,'(a)') 'ERROR: invalid displacement number in &
                  &phonon_disp_list'
             write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
                  &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
          end if
          call comms_abort
       end if
    end if

    !==========================================================================!
    ! STAGE 1: calculate forces for unperturbed system and check that they are !
    !          close to zero                                                   !
    !==========================================================================!

    if ((pub_phonon_farming_task==0) .or. (pub_phonon_farming_task==1)) then

       if (.not. write_denskern) then
          write_denskern=.true.
          if (pub_on_root) write(stdout,'(a)') 'Phonon calculation: &
               &write_denskern parameter over-ridden to TRUE'
       end if
       if (.not. write_tightbox_ngwfs .and. (maxit_ngwf_cg>0)) then
          write_tightbox_ngwfs=.true.
          if (pub_on_root) write(stdout,'(a)') 'Phonon calculation: &
               &write_tightbox_ngwfs parameter over-ridden to TRUE'
       end if

       call energy_and_force_calculate(total_energy,forces,elements, &
            return_converged=converged)

       ! Check initial calculation converged
       !if (.not.converged) then
       !   if (pub_on_root) then
       !      write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
       !           &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
       !      write(stdout,'(a)') 'ERROR in phonon_main: Calculation of initial &
       !           &configuration did not converge.'
       !      write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
       !           &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
       !   end if
       !   call comms_abort
       !end if

       ! Calculate |F|max for system
       Fmax=0.0_DP
       do i=1,pub_cell%nat
          Fmax=max(Fmax,dot_product(forces(1:3,i),forces(1:3,i)))
       end do
       Fmax=sqrt(Fmax)

       ! If |F|max is too large, abort
       if (Fmax>=pub_phonon_fmax) then
          if (pub_on_root) then
             write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
                  &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
             write(stdout,'(a)') 'ERROR: Forces for starting configuration are &
                  &too large for a meaningful phonon'
             write(stdout,'(a)') '       calculation. Please perform a &
                  &geometry optimization and try again.'
             write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
                  &-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
          end if
          call comms_abort
       end if

    end if

    !==================================================================!
    ! STAGE 2: displace atoms in turn and save force constants to file !
    !==================================================================!

    if ((pub_phonon_farming_task==0) .or. (pub_phonon_farming_task==2)) then

       if ((.not. read_denskern).and.pub_geom_reuse_dk_ngwfs) then
          read_denskern=.true.
          if (pub_on_root) write(stdout,'(a)') 'Phonon calculation: &
               &read_denskern parameter over-ridden to TRUE'
       end if
       if ((.not. read_tightbox_ngwfs) .and. pub_geom_reuse_dk_ngwfs .and. &
           (maxit_ngwf_cg>0)) then
          read_tightbox_ngwfs=.true.
          if (pub_on_root) write(stdout,'(a)') 'Phonon calculation: &
               &read_tightbox_ngwfs parameter over-ridden to TRUE'
       end if
       if (write_denskern) then
          write_denskern=.false.
          if (pub_on_root) write(stdout,'(a)') 'Phonon calculation: &
               &write_denskern parameter over-ridden to FALSE'
       end if
       if (write_tightbox_ngwfs) then
          write_tightbox_ngwfs=.false.
          if (pub_on_root) write(stdout,'(a)') 'Phonon calculation: &
               &write_tightbox_ngwfs parameter over-ridden to FALSE'
       end if

       ! Allocate array for holding a row of the force constants matrix
       if (pub_on_root) then
          allocate(force_consts(1:3,1:pub_cell%nat),stat=ierr)
          call utils_alloc_check('phonon','force_consts',ierr)
       end if

       ! Loop over atoms in unit cell
       do i=1,pub_cell%nat

          ! Loop over x/y/z directions
          do j=1,3

             if ((pub_phonon_farming_task==0) .or. &
                  (.not. pub_have_phonon_disp_list)) then
                do_disp=.true.
             else if (any(pub_phonon_disp_list(1:pub_num_disp)==(i-1)*3+j)) then
                do_disp=.true.
             else
                do_disp=.false.
             end if

             if (do_disp) then

                if (j==1) then
                   cart_char='x'
                else if (j==2) then
                   cart_char='y'
                else if (j==3) then
                   cart_char='z'
                end if

                if (pub_on_root) then
                   write(disp_number,*) (i-1)*3+j
                   write(output_file,*) trim(pub_rootname)//'.force_consts_'//&
                        trim(adjustl(disp_number))
                   output_unit=utils_unit()
                   open(unit=output_unit,form='unformatted',&
                        file=trim(output_file),action='write')
                end if

                ! Loop over +ve/-ve directions
                do k=0,1

                   if (k==0) then
                      dir_char='+'
                   else if (k==1) then
                      dir_char='-'
                   end if

                   if (pub_on_root) then
                      write(stdout,'()')
                      write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-&
                           &+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
                      write(stdout,'(a,i4,a,i4,a,a1,a,a1,a)') 'Displacing &
                           &atom ', i, ' of ', pub_cell%nat, ' in the ', &
                           dir_char, 've ', cart_char, '-direction'
                      write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-&
                           &+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
                      write(stdout,'()')
                   end if

                   ! displace atomic position in required direction
                   if (j==1) then
                      eql_pos=elements(i)%centre%x
                      elements(i)%centre%x=eql_pos+((-1.0_DP)**k)*&
                           pub_phonon_disp
                   else if (j==2) then
                      eql_pos=elements(i)%centre%y
                      elements(i)%centre%y=eql_pos+((-1.0_DP)**k)*&
                           pub_phonon_disp
                   else if (j==3) then
                      eql_pos=elements(i)%centre%z
                      elements(i)%centre%z=eql_pos+((-1.0_DP)**k)*&
                           pub_phonon_disp
                   end if

                   ! Calculate forces for this displacement
                   call energy_and_force_calculate(total_energy,forces,elements)

                   ! Check last calculation converged
                   !if (.not.converged) then
                   !   if (pub_on_root) then
                   !      write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-&
                   !           &+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
                   !      write(stdout,'(a)') 'ERROR in phonon_main: &
                   !           &Calculation of last configuration did not &
                   !           &converge.'
                   !      write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+&
                   !           &+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
                   !   end if
                   !   call comms_abort
                   !end if

                   ! Restore equilibrium position of displaced atom
                   if (j==1) then
                      elements(i)%centre%x=eql_pos
                   else if (j==2) then
                      elements(i)%centre%y=eql_pos
                   else if (j==3) then
                      elements(i)%centre%z=eql_pos
                   end if

                   ! Calculate force constants and write to file
                   if (pub_on_root) then
                      if (k==0) then
                         force_consts=forces
                      else
                         force_consts=-0.5_DP*(force_consts-forces)/&
                              pub_phonon_disp
                         write(output_unit) force_consts
                      end if
                   end if

                end do

                if (pub_on_root) close(output_unit)

             end if

          end do

       end do

       if (pub_on_root) then
          deallocate(force_consts,stat=ierr)
          call utils_dealloc_check('phonon','force_consts',ierr)
       end if

    end if

    !=================================================================!
    ! STAGE 3: construct dynamical matrix and find phonon frequencies !
    !=================================================================!

    if ((pub_phonon_farming_task==0) .or. (pub_phonon_farming_task==3)) then

       call dense_init

       ! allocate full 3Nx3N dynamical matrix
       call dense_create(dynamical_symm_dem,3*pub_cell%nat,3*pub_cell%nat,&
            .false.)
       call dense_create(dynamical_dem,3*pub_cell%nat,3*pub_cell%nat,.false.)

       ! read force constants back in from file and calculate dynamical matrix
       ! elements
       allocate(force_consts(1:3,1:pub_cell%nat),stat=ierr)
       call utils_alloc_check('phonon','force_consts',ierr)
       do i=1,pub_cell%nat
          do j=1,3
             if (pub_on_root) then
                write(disp_number,*) (i-1)*3+j
                write(output_file,*) trim(pub_rootname)//'.force_consts_'//&
                     trim(adjustl(disp_number))
                output_unit=utils_unit()
                open(unit=output_unit,form='unformatted',&
                     file=trim(output_file),action='read',status='old')
                read(output_unit) force_consts
             end if
             call comms_bcast(pub_root_node_id,force_consts)
             do k=1,pub_cell%nat
                do l=1,3
                   el=force_consts(l,k)/sqrt(periodic_table_mass(elements(i)%&
                        atomic_number)*periodic_table_mass(elements(k)%&
                        atomic_number)/(electron_mass_u**2))
                   call dense_put_element(el,dynamical_dem,3*(k-1)+l,3*(i-1)+j)
                end do
             end do
             if (pub_on_root) close(output_unit)
          end do
       end do
       deallocate(force_consts,stat=ierr)
       call utils_dealloc_check('phonon','force_consts',ierr)

       !symmetrize dynamical matrix
       do i=1,3*pub_cell%nat
          call dense_get_element(el,dynamical_dem,i,i)
          call dense_put_element(el,dynamical_symm_dem,i,i)
       end do

       do i=1,3*pub_cell%nat
          do j=i+1,3*pub_cell%nat
             call dense_get_element(el,dynamical_dem,i,j)
             call dense_get_element(el_T,dynamical_dem,j,i)
             el_symm=(el+el_T)*0.5_DP
             call dense_put_element(el_symm,dynamical_symm_dem,i,j)
             call dense_put_element(el_symm,dynamical_symm_dem,j,i)
          end do
       end do

       call dense_destroy(dynamical_dem)

       ! diagonalize dynamical matrix to find phonon frequencies
       call dense_create(eigenvecs_dem,3*pub_cell%nat,3*pub_cell%nat,.false.)
       allocate(freqs(1:3*pub_cell%nat),stat=ierr)
       call utils_alloc_check('phonon','freqs',ierr)
       call dense_create(unit_mat_dem,3*pub_cell%nat,3*pub_cell%nat,.false.)

       do i=1,3*pub_cell%nat
          do j=1,3*pub_cell%nat
             if (i==j) then
                call dense_put_element(1.0_DP,unit_mat_dem,j,i)
             else
                call dense_put_element(0.0_DP,unit_mat_dem,j,i)
             end if
          end do
       end do

       call dense_eigensolve(3*pub_cell%nat,freqs,dynamical_symm_dem,&
            unit_mat_dem,1,eigenvecs_dem)

       call dense_destroy(unit_mat_dem)

       ! print out phonon frequencies
       if (pub_on_root) then
          write(stdout,'()')
          write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-&
               &+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
          write(stdout,'(a)') 'Phonon frequencies               (cm^-1)        &
               &         (THz)'
          write(stdout,'()')
          do i=1,3*pub_cell%nat
             freq_conv=sign(sqrt(abs(freqs(i)))/TWO_PI,freqs(i))
             write(stdout,'(12x,i5,a,2(6x,f16.11))') i, ":", &
                  freq_conv*au2inv_cm, freq_conv*au2THz
          end do
          write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-&
               &+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'

          ! ndmh: print QC test data
          if (print_qc) then
             i=1
             freq_conv=sign(sqrt(abs(freqs(i)))/TWO_PI,freqs(i))
             write(stdout,'(a,f24.12)') &
                  '<QC>               [phonon_freq_1]: ',freq_conv
             i=3*pub_cell%nat
             freq_conv=sign(sqrt(abs(freqs(i)))/TWO_PI,freqs(i))
             write(stdout,'(a,f24.12)') &
                  '<QC>              [phonon_freq_3N]: ',freq_conv
          end if

       end if

       ! calculate and print out zero-point energy
       if (pub_on_root) then
          zero_point_E=0.0_DP
          print_warning=.true.
          do i=1,3*pub_cell%nat
             if (freqs(i)<0.0_DP) then
                if (print_warning) then
                   print_warning=.false.
                   write(stdout,'(a)') 'WARNING: Imaginary phonon frequencies &
                        &excluded from the computation of the zero-point energy'
                end if
             else
                zero_point_E=zero_point_E+sqrt(freqs(i))
             end if
          end do
          zero_point_E=0.5_DP*zero_point_E
          write(stdout,'()')
          write(stdout,'(a,1x,f16.11,1x,a)') "Zero-point energy =", &
               zero_point_E, "Eh"

          ! ndmh: print QC test data
          if (print_qc) write(stdout,'(a,f24.12)') &
                  '<QC>           [zero_point_energy]: ',zero_point_E

       end if

       ! calculate and print out F, S, U, and C_v for specified range of T
       if (pub_on_root) then
          T_num_iter=nint((pub_phonon_tmax-pub_phonon_tmin)/pub_phonon_deltat)
          T=pub_phonon_tmin
          print_warning=.false.
          write(stdout,'()')
          write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-&
               &+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-'
          write(stdout,'(a)') '                             Helmholtz        &
               &                            Internal              Specific'
          write(stdout,'(a)') ' Temperature (K)      free energy (Eh)        &
               &Entropy (Eh/K)           energy (Eh)           heat (Eh/K)'
          write(stdout,'()')
          do i=0,T_num_iter
             if (T==0.0_DP) then
                write(stdout,'(f16.11,4(6x,f16.11))') 0.0_DP, zero_point_E, &
                     0.0_DP, zero_point_E, 0.0_DP
             else
                beta=1.0_DP/T
                F=0.0_DP
                S=0.0_DP
                C_v=0.0_DP
                do j=1,3*pub_cell%nat
                   sqrt_freq=sign(sqrt(abs(freqs(j))),freqs(j))
                   if (sqrt_freq/TWO_PI<pub_phonon_min_freq) then
                      print_warning=.true.
                   else
                      exp_freq=exp(-beta*sqrt_freq)
                      F=F+log(1.0_DP-exp_freq)
                      S=S+sqrt_freq*exp_freq/(1.0_DP-exp_freq)
                      C_v=C_v+sqrt_freq**2*exp_freq/((1.0_DP-exp_freq)**2)
                   end if
                end do
                F=T*F
                S=S-F
                F=F+zero_point_E
                U=F+S
                S=k_B*S/T
                C_v=k_B*C_v/(T**2)
                write(stdout,'(f16.11,4(6x,f16.11))') T/k_B, F, S, U, C_v
             end if
             T=T+pub_phonon_deltat
          end do

          ! ndmh: print QC test data
          if (print_qc) then
             T=T-pub_phonon_deltat
             write(stdout,'(a,f24.12)') &
                '<QC>                 [temperature]: ',T/k_B
             write(stdout,'(a,f24.12)') &
                '<QC>               [helmholtz_f_e]: ',F
             write(stdout,'(a,f24.12)') &
                '<QC>                     [entropy]: ',S
             write(stdout,'(a,f24.12)') &
                '<QC>                  [internal_e]: ',U
             write(stdout,'(a,f24.12)') &
                '<QC>               [specific_heat]: ',C_v
          end if

          write(stdout,'(a)') '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-&
               &+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-'
          if (print_warning) then
             write(stdout,'(a)') 'WARNING: Low phonon frequencies &
                  &excluded from the computation of thermodynamic quantities'
          end if
       end if

       deallocate(freqs,stat=ierr)
       call utils_dealloc_check('phonon','freqs',ierr)
       call dense_destroy(eigenvecs_dem)
       call dense_destroy(dynamical_symm_dem)
       call dense_exit

    end if

  end subroutine phonon_main

end module phonon
