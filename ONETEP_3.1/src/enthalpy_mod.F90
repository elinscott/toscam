! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!               Electronic enthalpy module                       !
!----------------------------------------------------------------!
! Calculates an electronic enthalpy term according to the        !
! electronic isosurface method of Cococcioni et al.              !
!                                                                !
! Written by Niccolo Corsini in July 2011.                       !
!================================================================!

module enthalpy

  use constants, only: DP, VERBOSE

  implicit none

  private

  public :: enthalpy_terms
  public :: enthalpy_volume_test

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine enthalpy_terms(pot_fine, enthalpy_energy, density_fine,grid)

    !=========================================================================!
    ! Calculates the volume of the system by integrating the electronic       !
    ! density up to an isosurface defined by a cutoff value. A smeared step   !
    ! function is used for this purpose with parameter alpha defining the     !
    ! smear. Pressure is defined as an input parameter and we can thus        !
    ! calculate the PV contribution both to the energy and to the potential   !
    ! by taking the functional derivative with respect to the density.        !
    !-------------------------------------------------------------------------!
    ! Written by Niccolo Corsini in July 2011.                                !
    ! Minor modifications by Nicholas Hine, October 2011.                     !
    !=========================================================================!

    use cell_grid, only: GRID_INFO
    use comms, only: pub_on_root
    use constants, only: DN, UP, PI, stdout
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_external_pressure,pub_smoothing_factor, &
         pub_isosurface_cutoff, pub_output_detail
    use simulation_cell, only: pub_cell
    use utils, only: utils_abort,utils_alloc_check, utils_dealloc_check, &
         utils_erfc

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(out)   :: pot_fine(grid%ld1,grid%ld2, &
         grid%max_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(inout) :: density_fine(grid%ld1, &
         grid%ld2, grid%max_slabs12, pub_cell%num_spins)
    real(kind=DP), intent(out), optional :: enthalpy_energy
   
    ! Local variables
    integer :: ierr,i1,i2,i3  
    real(kind=DP) :: volume
    real(kind=DP) :: prefac
    real(kind=DP), allocatable, dimension(:,:,:) :: step_func

    ! nc: allocate workspace
    allocate(step_func(grid%ld1,grid%ld2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('enthalpy_terms','step_func',ierr)

    ! ndmh: combine both spin channels in density_fine(:,:,:,1)
    if (pub_cell%num_spins==2) then
       density_fine(:,:,:,1) = density_fine(:,:,:,1) + density_fine(:,:,:,2)
    end if

    ! nc: calculate step function of density with step at pub_isosurface_cutoff
    do i3=1,grid%num_my_slabs12
       do i2=1,grid%n2
          do i1=1,grid%n1
             step_func(i1,i2,i3) = 0.5_DP * utils_erfc( &
                  sqrt((pub_smoothing_factor**2) / &
                  (2.0_DP*(pub_isosurface_cutoff**2)))* & 
                  (pub_isosurface_cutoff-density_fine(i1,i2,i3,1)))
          end do
       end do
    end do
   
    ! nc: calculate 'volume' of density above pub_isosurface_cutoff
    volume = integrals_product_on_grid(grid,step_func)

    ! nc: calculate potential associated with this enthalpy
    prefac = pub_smoothing_factor*pub_external_pressure / &
         (pub_isosurface_cutoff*sqrt(2.0_DP*PI))
    pot_fine = prefac * exp((-(pub_smoothing_factor**2)* & 
         ((density_fine-pub_isosurface_cutoff)**2)) / &
         (2.0_DP*(pub_isosurface_cutoff**2)))

    ! nc: calculate enthalpy associated with this volume
    enthalpy_energy = pub_external_pressure * volume

    ! nc: report results
    if ((pub_output_detail>=VERBOSE).and.(pub_on_root)) then
       write(stdout,'(a,f20.12)') '  Vol. Enthalpy: ',enthalpy_energy
    end if

    ! nc: deallocate workspace
    deallocate(step_func,stat=ierr)
    call utils_dealloc_check('enthalpy_terms','step_func',ierr)

    ! ndmh: undo combination of spin channels
    if (pub_cell%num_spins==2) then
       density_fine(:,:,:,1) = density_fine(:,:,:,1) - density_fine(:,:,:,2)
       pot_fine(:,:,:,2) = pot_fine(:,:,:,1)
    end if
  
  end subroutine enthalpy_terms


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine enthalpy_volume_test(density_fine,grid)

    use cell_grid, only: GRID_INFO
    use constants, only: DP, DN, UP, PI, stdout
    use integrals, only: integrals_product_on_grid 
    use rundat, only:pub_devel_code,pub_external_pressure, &
         pub_smoothing_factor,pub_isosurface_cutoff
    use simulation_cell, only: pub_cell
    use utils, only: utils_abort, utils_alloc_check, utils_dealloc_check, & 
         utils_erfc

    implicit none

    ! Arguments
    type(GRID_INFO), intent(in) :: grid
    real(kind=DP), intent(inout) :: density_fine(grid%ld1, &
         grid%ld2, grid%max_slabs12, pub_cell%num_spins)

    ! Local variables
    integer :: ierr,i1,i2,i3
    real(kind=DP) :: volume_sharp_step
    real(kind=DP) :: volume_smooth_step
    real(kind=DP), allocatable, dimension(:,:,:) :: step_func

    ! Allocate workspace
    allocate(step_func(grid%ld1,grid%ld2,grid%max_slabs12),stat=ierr)
    call utils_alloc_check('enthalpy_volume_test','step_func',ierr)

    do i3=1,grid%num_my_slabs12
       do i2=1,grid%n2
          do i1=1,grid%n1
            step_func(i1,i2,i3)=0.5_DP*utils_erfc( &
                 sqrt((pub_smoothing_factor**2) / &
                 (2.0_DP*(pub_isosurface_cutoff**2)))* &
                 (pub_isosurface_cutoff-density_fine(i1,i2,i3,1)))
          end do
       end do
    end do
   
    volume_smooth_step = integrals_product_on_grid(grid,step_func)

    volume_sharp_step = integrals_product_on_grid(grid,density_fine)

    write(stdout,*) ' The smoothed volume is', volume_smooth_step, &
         '  The sharp volume is', volume_sharp_step

    ! Deallocate workspace
    deallocate(step_func,stat=ierr)
    call utils_dealloc_check('enthalpy_volume_test','step_func',ierr)
  
  end subroutine enthalpy_volume_test


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module enthalpy

