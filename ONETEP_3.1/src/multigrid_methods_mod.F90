!=============================================================================!
!                          MULTIGRID_METHODS MODULE                           !
!=============================================================================!
! This module defines the onetep interface to the multigrid solver and        !
! subroutines to calculate electrostatic quantities that are related to the   !
! solution of the Poisson equation.                                           !
!                                                                             !
! For more info on multigrid methods see the following useful references      !
!                                                                             !
! [1] 'A Multigrid Tutorial' Briggs, Henson and McCormick.                    !
!                                                                             !
! [2] TL Beck 'Real-space mesh techniques in density-functional theory'       !
!     Rev Mod Phys 72 pg 1041.                                                !
!                                                                             !
! Portions of this module use techniques or formulas described in             !
! [3] Scherlis, Fattebert, Gygi, Cococcioni and Marzari, J. Chem. Phys. 124,  !
!     074103 (2006).                                                          !
!-----------------------------------------------------------------------------!
! Written for castep by Hatem H Helal, starting in 10/2007                    !
! Adapted for onetep by Jacek Dziedzic, 04-05/2010.                           !
! Made parallel-ready by Jacek Dziedzic, 04-05/2010.                          !
!                                                                             !
! Implementation of defect correction solver by HH Helal 06/2010              !
! Described in some detail within the following references:                   !
! [4] U Trottenberg, C Oosterlee and A Schuller, "Multigrid", 2001.           !
! [5] S Schaffer, Math of Comp, 43, 89-115 (1984).                            !
!-----------------------------------------------------------------------------!

module multigrid_methods

  use constants, only: DP
  use finite_differences, only: FD_GRID_INFO
  use multigrid, only: dummy

  implicit none

  private

  ! jd: Initialisation flags
  logical :: multigrid_initialised = .false.

  ! jd: Debug flag, set to true by DEVEL_CODE MG: ... :MG
  logical :: multigrid_debug = .false.

  !---------------------------------------------------------------------------!
  !                   P u b l i c   D a t a s t r u c t u r e s               !
  !---------------------------------------------------------------------------!

  ! jd: Geometry of the multigrid
  !     Mostly used privately, but is_solvation needs it to calculate FD
  !     gradients on the MG grid.
  type(FD_GRID_INFO), public, save :: mg

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!

  public  ::  multigrid_initialise
  public  ::  multigrid_calculate_hartree
  public  ::  multigrid_defect_corr_solver
  public  ::  multigrid_electrostatic_energy
  public  ::  multigrid_prepare_bound_cond

contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine multigrid_initialise
    !=========================================================================!
    ! Initialises the multigrid.                                              !
    ! Calculates the size of the multigrid, sanity-checks that the            !
    ! distribution makes sense, examines DEVEL_CODE.                          !
    ! There is no need to call this subroutine, it will be automatically      !
    ! called the first time you use any of multigrid_*() routines, but you can!
    ! call it explicitly, any number of times, if so desired.                 !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, April 2010.                                  !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use constants, only: stdout
    use comms, only: pub_on_root, pub_root_node_id, comms_bcast
    use finite_differences, only: finite_difference_initialise, &
         finite_difference_set_geometry
    use rundat, only: pub_is_multigrid_nlevels, pub_devel_code
    use simulation_cell, only: pub_cell
    use utils, only: utils_assert, utils_flush, utils_trace_in, utils_trace_out

    implicit none

    ! jd: Local variables
    integer            :: multigrid_granularity
    character(len=200) :: mg_devel_code
    integer            :: start_pos, stop_pos, test_pos

    ! ------------------------------------------------------------------------

    if(multigrid_initialised) return

    call utils_trace_in('multigrid_initialise')

    multigrid_granularity = 2**(pub_is_multigrid_nlevels+1)

    ! jd: Initialise the MG grid
    call finite_difference_initialise
    mg = finite_difference_set_geometry(pub_fine_grid,multigrid_granularity)

    if(pub_on_root) then
       write(stdout,'(a)',advance='no') 'ONETEP fine grid is '
       write(stdout,'(i0,a,i0,a,i0)',advance='no') &
            mg%pt1f,' x ',mg%pt2f,' x ',mg%pt3f
       write(stdout,'(a,F0.4,a,F0.4,a,F0.4,a)') ' gridpoints, ', &
            mg%pt1f*mg%d1f, ' x ',mg%pt2f*mg%d2f,' x ', mg%pt3f*mg%d3f,' bohr.'
       write(stdout,'(a)',advance='no') 'FD multigrid is     '
       write(stdout,'(i0,a,i0,a,i0)',advance='no') &
            mg%pq1f,' x ',mg%pq2f,' x ',mg%pq3f
       write(stdout,'(a,F0.4,a,F0.4,a,F0.4,a/)') ' gridpoints, ', &
            mg%pq1f*mg%d1f, ' x ',mg%pq2f*mg%d2f,' x ', mg%pq3f*mg%d3f,' bohr.'

       if(mg%pt1f-mg%pq1f>10 .or. mg%pt2f-mg%pq2f>10 .or. &
            mg%pt3f-mg%pq3f>10) then
          write(stdout,'(a)') 'WARNING: A considerable portion of the fine grid&
               & will be hidden from the'
          write(stdout,'(a)') 'multigrid calculation. If you''re sure that&
               & there is no charge density in the '
          write(stdout,'(a)') 'margin that constitues the above difference,&
               & you''re OK. If not, try to increase '

          if(mg%pt1f-mg%pq1f>10) then
             write(stdout,'(a,i0,a,i0,a)') 'the 1st dimension, ',mg%pt1f, &
                  ', to ', mg%pq1f+multigrid_granularity+1,','
          end if
          if(mg%pt2f-mg%pq2f>10) then
             write(stdout,'(a,i0,a,i0,a)') 'the 2nd dimension, ',mg%pt2f, &
                  ', to ', mg%pq2f+multigrid_granularity+1,','
          end if
          if(mg%pt3f-mg%pq3f>10) then
             write(stdout,'(a,i0,a,i0,a)') 'the 3rd dimension, ',mg%pt3f, &
                  ', to ', mg%pq3f+multigrid_granularity+1,','
          end if
          write(stdout,'(a)') 'and re-run.'
       end if

       ! jd: Check if in devel-code low-level debug mode
       mg_devel_code=pub_devel_code
       if (len_trim(mg_devel_code)>0) then
          start_pos=index(mg_devel_code,'MG:')
          stop_pos=index(mg_devel_code,':MG')
          if (stop_pos<=0) stop_pos=len_trim(mg_devel_code)
          if (start_pos>0) then
             test_pos=index(mg_devel_code,'DEBUG_DUMP')
             if (test_pos>start_pos.and.test_pos<stop_pos) then
                multigrid_debug = .true.
                write(stdout,'(a)') &
                     'NB: Multigrid debugging dumps will be produced.'
             end if
          end if
       end if

    end if

    call comms_bcast(pub_root_node_id,multigrid_debug)

    call utils_flush

    multigrid_initialised = .true.

    call utils_trace_out('multigrid_initialise')

  end subroutine multigrid_initialise
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine parallel_multigrid_solver(pot, residual, rho, eps_half, &
       bound, tolerance_override)
    !=========================================================================!
    ! For a given electronic density (rho) and dielectric functional (eps)    !
    ! this subroutine will compute the electrostatic potential (pot) which is !
    ! the solution to the Poisson equation:                                   !
    !                                                                         !
    !        -  \nabla \cdot [ eps \nabla (pot) ] = 4*pi*rho                  !
    !                                                                         !
    ! As the name implies a multigrid solver is used to find the potential    !
    ! and the parameter pub_is_discretization_order will select the level of  !
    ! accuracy for the resulting discretized equations.  The high order defect!
    ! method is used to obtain the higher order accuracy result.  This method !
    ! is summarized in the references [4] and [5].                            !
    !                                                                         !
    ! This routine takes care of the distributed vs. serial nastiness -- it   !
    ! takes care of gathering the distributed arrays ld1f x ld2f x max_slabs12!
    ! onto a local grid, suitably sized for the solver (pq1f x pq2f x pq3f)   !
    ! and then scatters the solution to the distributed arrays.               !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   pot,               (output),     the potential we seek                !
    !                                                                         !
    !   residual,          (output)      useful for error analysis            !
    !                                                                         !
    !   rho,               (input),      the input electron density           !
    !                                                                         !
    !   eps_half,          (input),      the dielectric functional            !
    !                      optional      interpolated to 'half' grids         !
    !                      Defaults to 1.0 everywhere if omitted              !
    !                                                                         !
    !   bound,             (input)       initial guess for the potential      !
    !                      optional      and fixed boundary condition         !
    !                      Defaults to 0.0 everywhere if omitted              !
    !                                                                         !
    !   tolerance_override (input)       tolerance for stop criterion in the  !
    !                      optional      mg solver.                           !
    !                      Defaults to pub_is_multigrid_error_tol if omitted. !
    !                                                                         !
    ! All the arrays are assumed to be in the distributed slab representation.!
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal, 19/10/2007                                    !
    ! Adapted for onetep by Jacek Dziedzic, 04-05/2010                        !
    ! Made parallel-ready by Jacek Dziedzic, 04-05/2010                       !
    ! Updated for defect correction by HH Helal, 06/2010                      !
    ! Adapted for fully parallel operation by Jacek Dziedzic, 03/2011.        !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: pub_on_root, pub_my_node_id, comms_barrier
    use constants, only: pi, stdout, VERBOSE
    use geometry, only: operator(.DOT.)
    use integrals, only: integrals_product_on_grid
    !use multigrid, only: multigrid_solver
    use rundat, only: pub_is_multigrid_error_tol, pub_is_multigrid_max_iters, &
         pub_is_multigrid_nlevels, pub_output_detail, &
         pub_is_discretization_order
    use simulation_cell, only: pub_cell
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_assert, utils_alloc_check, &
         utils_dealloc_check, utils_unit, utils_sanity_check, utils_trace_in, &
         utils_trace_out

    implicit none

    ! --------------------------------------------------------------------------
    ! Output arguments
    real(kind=DP), dimension(mg%ld1f, mg%ld2f, mg%max_slabs12), intent(out) &
         :: pot
    real(kind=DP), dimension(mg%ld1f, mg%ld2f, mg%max_slabs12), intent(out) &
         :: residual

    ! Input arguments
    real(kind=DP), dimension(mg%ld1f, mg%ld2f, mg%max_slabs12), &
         intent(in)                                                 :: rho
    real(kind=DP), dimension(mg%ld1f, mg%ld2f, mg%max_slabs12, 3), &
         intent(in), optional, target                               :: eps_half
    real(kind=DP), dimension(mg%ld1f, mg%ld2f, mg%max_slabs12), &
         intent(in), optional                                       :: bound
    real(kind=DP), intent(in), optional :: tolerance_override

    ! --------------------------------------------------------------------------
    ! Internal variables
    real(kind=DP), dimension(:,:,:,:), pointer :: eps_half_loc

    real(kind=DP) :: res_norm   ! jd: Norm of the residual for the defect-corr'n
    integer :: ierr             ! jd: Error flag
    real(kind=DP) :: tolerance  ! jd: Requested multigrid tolerance
    character(len=*), parameter :: myself = 'parallel_multigrid_solver'

    integer       :: mgunit     ! jd: Output unit for mg_output_file
    character(len=*), parameter :: mg_output_file = 'multigrid_iteration.data'

    !------------------------------------------------------------------------

    call utils_trace_in(myself)
    call comms_barrier()    ! jd: Synchronize to get accurate timings for solver
    call timer_clock(myself,1)

    ! jd: Break, if multigrid not yet initialised. It's too late to initialise
    !     now, since the initialisation takes care of proper array dimensioning.
    call utils_assert(multigrid_initialised, &
         'Must call multigrid_initialise() first.')

    ! jd: Sanity checks on input
    call utils_sanity_check(rho,'rho in '//myself)
    if(present(eps_half)) then
       call utils_sanity_check(eps_half(:,:,:,1),'eps_half/x in '//myself)
       call utils_sanity_check(eps_half(:,:,:,2),'eps_half/y in '//myself)
       call utils_sanity_check(eps_half(:,:,:,3),'eps_half/z in '//myself)
    end if
    if(present(bound)) call utils_sanity_check(bound,'bound in '//myself)

    ! hhh: We have to be careful here as we implicitly assume our cell is
    !      orthorombic for using the multigrid solver
    ! jd:  Explicitly check if we are safe
    call utils_assert(&
         (pub_cell%a1_unit .dot. pub_cell%a2_unit) == 0.0_DP .and. &
         (pub_cell%a1_unit .dot. pub_cell%a3_unit) == 0.0_DP .and. &
         (pub_cell%a2_unit .dot. pub_cell%a3_unit) == 0.0_DP, &
         'The multigrid solver only supports orthorombic cells, sorry.')

    ! jd: If boundary conditions are specified in 'bound', use them,
    !     else use zero BC. Since bound is intent(in) and the solver expects
    !     the BCs in the pot array, let's use the pot array.
    !     First zero pot to take care of padding.
    pot = 0D0
    if(present(bound)) then
       pot(1:mg%pq1f,1:mg%pq2f,1:mg%num_slabs12_pq) = &
            bound(1:mg%pq1f,1:mg%pq2f,1:mg%num_slabs12_pq)
    end if

    ! jd: If permittivity specified in 'eps_half', use it,
    !     else use unity everywhere, but must allocate the array.
    if(present(eps_half)) then
       eps_half_loc => eps_half
    else
       allocate(eps_half_loc(mg%ld1f,mg%ld2f,mg%max_slabs12,3),stat=ierr)
       call utils_alloc_check(myself,'eps_half_loc',ierr)
       eps_half_loc = 1.0_DP
    end if

    ! jd: Allow the number of iterations to be overridden, so that crude
    !     solutions are possible from within the defect_corr_solver
    if(present(tolerance_override)) then
       tolerance = tolerance_override
    else
       tolerance = pub_is_multigrid_error_tol
    end if

    ! jd: Take care of padding in residual
    residual = 0D0

    ! jd: Open the multigrid details file
    mgunit=-1
    if(pub_on_root) then
       mgunit = utils_unit()
       open(mgunit, file=mg_output_file, position='append', err=10)
    end if

    ! **********************************************************************
    ! **********************************************************************
    ! jd: Call the parallel solver
    ! **********************************************************************
    ! **********************************************************************
    call utils_abort('Error in parallel_multigrid_solver: Multigrid not &
         &present in this version.')
    !call multigrid_solver( &
    !     eps_half_loc(1:mg%pq1f,1:mg%pq2f,1:mg%num_slabs12_pq,:), &
    !     pot(1:mg%pq1f,1:mg%pq2f,1:mg%num_slabs12_pq), &
    !     rho(1:mg%pq1f,1:mg%pq2f,1:mg%num_slabs12_pq)*(4.0*PI), &
    !     residual(1:mg%pq1f,1:mg%pq2f,1:mg%num_slabs12_pq), &
    !     mg%pq1f,mg%pq2f,mg%pq3f, &
    !     pub_fine_grid%first_slab12(pub_my_node_id), &
    !     min(mg%pq3f,pub_fine_grid%last_slab12(pub_my_node_id)), &
    !     mg%d1f,mg%d2f,mg%d3f,tolerance, &
    !     mgunit,pub_is_multigrid_max_iters, &
    !     pub_is_multigrid_nlevels)
    ! **********************************************************************
    ! **********************************************************************

    ! jd: Close the multigrid details file
    if(pub_on_root) then
       close(mgunit, err=20)
    end if

    ! jd: Calculate the residual norm
    call utils_sanity_check(residual,'residual in '//myself)
    res_norm = sqrt(integrals_product_on_grid(pub_fine_grid,&
         residual,residual,mg%pq1f,mg%pq2f,mg%num_slabs12_pq))

    if(.not. present(eps_half)) then
       deallocate(eps_half_loc,stat=ierr)
       call utils_dealloc_check(myself,'eps_half_loc',ierr)
    end if

    call utils_trace_out(myself)
    call timer_clock(myself,2)

    return

    ! jd: I/O error handling
10  call utils_abort('Error during creation of file: '//mg_output_file)
20  call utils_abort('Error during closing of file: '//mg_output_file)

  end subroutine parallel_multigrid_solver
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine multigrid_defect_corr_solver(pot, res, & ! out
       rho, bound, eps_full, eps_half)  ! in
    !=========================================================================!
    ! For a given electronic density (rho) and dielectric functional (eps)    !
    ! this subroutine will compute the electrostatic potential (pot) which is !
    ! the solution to the Poisson equation:                                   !
    !                                                                         !
    !        -  \nabla \cdot [ eps \nabla (pot) ] = 4*pi*rho                  !
    !                                                                         !
    ! This subroutine uses the defect correction method with multigrid to     !
    ! obtain a solution to the Poisson equation which has high order accuracy !
    ! determined by the finite difference operators used.                     !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   pot,       intent=out,             the potential we seek              !
    !   res,       intent=out,             useful for error analysis          !
    !              optional
    !   rho,       intent=in,              the density                        !
    !                                                                         !
    !   bound,     intent=in,              the initial guess for the potential!
    !              optional                and its fixed boundary conditions  !
    !              If omitted, will default to 0.0 everywhere.                !
    !   eps_full,  intent=in,              the dielectric functional on       !
    !              optional                'full' grid points                 !
    !              If omitted, will default to 1.0 everywhere.                !
    !   eps_half,  intent=in,              the dielectric functional          !
    !              optional                interpolated to 'half' grids       !
    !              If omitted, will default to 1.0 everywhere.                !
    !-------------------------------------------------------------------------!
    ! Defect correction written for castep by HH Helal 29/04/2010.            !
    ! Mended for onetep by HH Helal 7/6/2010.                                 !
    ! Adapted to work in parallel by Jacek Dziedzic 11/06/2010.               !
    !=========================================================================!

    use constants, only: pi, stdout, VERBOSE
    use cell_grid, only: pub_fine_grid
    use comms, only: pub_on_root
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_is_multigrid_error_tol, &
         pub_is_multigrid_defect_error_t, pub_is_multigrid_max_iters, &
         pub_output_detail, pub_is_discretization_order
    use timer, only: timer_clock
    use utils, only: utils_abort, utils_assert, utils_flush, &
         utils_alloc_check, utils_dealloc_check, utils_sanity_check, &
         utils_trace_in, utils_trace_out
    use visual, only: visual_scalarfield

    implicit none

    ! Output variables
    real(kind=DP), dimension(mg%ld1f,mg%ld2f,mg%max_slabs12), &
         intent(out)           :: pot
    real(kind=DP), dimension(mg%ld1f,mg%ld2f,mg%max_slabs12), &
         intent(out), optional :: res

    ! Input variables
    real(kind=DP), dimension(mg%ld1f,mg%ld2f,mg%max_slabs12), &
         intent(in)            :: rho
    real(kind=DP), dimension(mg%ld1f,mg%ld2f,mg%max_slabs12), &
         intent(in), optional  :: bound
    real(kind=DP), dimension(mg%ld1f,mg%ld2f,mg%max_slabs12), &
         intent(in), optional  :: eps_full
    real(kind=DP), dimension(mg%ld1f,mg%ld2f,mg%max_slabs12,3), &
         intent(in), optional  :: eps_half

    ! Internal variables
    integer       :: ierr ! jd: Error flag
    real(kind=DP) :: res_norm, defect_norm, error_norm
    character(len=*), parameter :: myself = 'multigrid_defect_corr_solver'
    character(len=*), parameter :: explanation = 'This could mean that the&
         & multigrid is too coarse to represent the quickly changing&
         & quantities (such as the permittivity).'

    ! Defect correction method variables
    real(kind=DP), allocatable :: defect(:,:,:) ! hhh: curr. defect (residual)
    real(kind=DP), allocatable :: error(:,:,:)  ! hhh: curr. error approx.
    real(kind=DP), allocatable :: eps_full_1(:,:,:) ! jd: Eps of unity
    real(kind=DP), allocatable :: res_local(:,:,:)  ! jd: Workspace
    integer                    :: counter
    logical                    :: converged

    !------------------------------------------------------------------------

    call timer_clock(myself,1)
    call utils_trace_in(myself)

    ! jd: Break, if multigrid not yet initialised. It's too late to initialise
    !     now, since the initialisation takes care of proper array dimensioning.
    call utils_assert(multigrid_initialised, &
         'Must call multigrid_initialise() first.')

    ! jd: Sanity checks
    call utils_sanity_check(rho,'rho in '//myself)
    if(present(bound)) call utils_sanity_check(bound, 'bound in '//myself)
    if(present(eps_full)) call utils_sanity_check(eps_full, &
         'eps_full in '//myself)
    if(present(eps_half)) then
       call utils_sanity_check(eps_half(:,:,:,1), 'eps_half/x in '//myself)
       call utils_sanity_check(eps_half(:,:,:,2), 'eps_half/y in '//myself)
       call utils_sanity_check(eps_half(:,:,:,3), 'eps_half/z in '//myself)
    end if

    if(multigrid_debug) then
       call visual_scalarfield(rho,pub_fine_grid, &
            'Total density','_rho_in_mg-dcsolver')
       if(present(eps_half)) then
          call visual_scalarfield(eps_half(:,:,:,1),pub_fine_grid, &
               'x-shifted permittivity', '_epshalfx_in_mg-dcsolver')
          call visual_scalarfield(eps_half(:,:,:,2),pub_fine_grid, &
               'y-shifted permittivity', '_epshalfy_in_mg-dcsolver')
          call visual_scalarfield(eps_half(:,:,:,3),pub_fine_grid, &
               'z-shifted permittivity', '_epshalfz_in_mg-dcsolver')
       end if
       if(present(eps_full)) then
          call visual_scalarfield(eps_full,pub_fine_grid, &
               'Permittivity','_epsfull_in_mg-dcsolver')
       end if
       if(present(bound)) then
          call visual_scalarfield(bound,pub_fine_grid, &
               'Boundary conditions','_bound_in_mg-dcsolver')
       end if
    end if

    if(pub_on_root .and. pub_output_detail == VERBOSE) then
       write(stdout,'(a)',advance='no') 'MG Solving with 2nd order solver... '
       call utils_flush(stdout,.true.)
    end if

    ! jd: If res is not specified, will need to allocate this workspace
    if(.not. present(res)) then
       allocate(res_local(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
       call utils_alloc_check(myself,'res_local',ierr)
       res_local = 0.0_DP ! jd: Takes care of max-num, ld-pt padding
    end if

    ! jd: Solve the equation with the 2nd order solver
    if(present(res)) then
       call parallel_multigrid_solver(pot,res,rho,eps_half,bound)
    else
       ! NB: This residual is ignored if the discretization order is >2
       call parallel_multigrid_solver(pot,res_local,rho,eps_half,bound)
    end if
    if(pub_on_root .and. pub_output_detail == VERBOSE) then
       write(stdout,'(a)') 'done'
       call utils_flush(stdout,.true.)
    end if

    ! jd: If higher order was requested, employ defect correction
    if (pub_is_discretization_order > 2) then

       ! Allocate memory for local arrays needed for defect correction
       allocate(defect(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
       call utils_alloc_check(myself,'defect',ierr)
       allocate(error(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
       call utils_alloc_check(myself,'error',ierr)

       ! jd: If eps_full is not specified, assume unity everywhere
       if(.not. present(eps_full)) then
          allocate(eps_full_1(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
          call utils_alloc_check(myself,'eps_full_1', &
               ierr)
          eps_full_1 = 1.0_DP
       end if

       ! hhh: Error initialized to zero -this is fine for possibly arbitrary BCs
       !      set earlier by optional bound argument in multigrid_poisson_solver
       error(:,:,:) = 0.0_DP

       ! hhh: We now have all the inputs ready to start the defect correction
       !      iteration
       converged = .false.
       counter = 0

       if (pub_on_root .and. pub_output_detail == VERBOSE ) then
          write(stdout,'(a,i0)') 'MG Starting defect correction with &
               &a discretization order of: ', pub_is_discretization_order
          write(stdout,'(a)') 'MG Iteration ***   |error|   ***  |defect|  ***&
               & |error eq''n residual|  <-- D I'
       end if

       ! jd: Defect correction loop
       do while (.not. converged)

         counter = counter + 1

         ! hhh: Compute the defect using high order finite difference
         !      Level of accuracy is selected through the input param
         !      discretization_order
         !      defect = source + \nabla \cdot [\epsilon \grad potmg]
         if(present(eps_full)) then
            call compute_current_defect(defect, defect_norm, &
                 pub_is_discretization_order, rho*4.0_DP*pi, eps_full, pot)
         else ! jd: Assume unity everywhere
            call compute_current_defect(defect, defect_norm, &
                 pub_is_discretization_order, rho*4.0_DP*pi, eps_full_1, pot)
         end if

         call utils_sanity_check(defect,'defect in '//myself)

         ! hhh: Now that we have the high order defect we very approximately
         !      solve the error equation using second order multigrid method
         !      with a low convergence threshold
         !      - \nabla \cdot [epsilon \grad error] = defect
         defect = defect/(4.0_DP*pi)
         if(present(res)) then
            call parallel_multigrid_solver(error,res,defect,&
                 eps_half,tolerance_override = pub_is_multigrid_defect_error_t)
         else
            call parallel_multigrid_solver(error,res_local,defect,&
                 eps_half,tolerance_override = pub_is_multigrid_defect_error_t)
         end if
         defect = defect*(4.0_DP*pi)

         ! Correct the current approximation to get the new approximation
         pot(:,:,:) = pot(:,:,:) + error(:,:,:)

         ! hhh: Compute error_norm and res_norm
         error_norm = sqrt(integrals_product_on_grid(pub_fine_grid, &
              error,error,mg%pq1f,mg%pq2f,mg%num_slabs12_pq))
         res_norm = sqrt(integrals_product_on_grid(pub_fine_grid, &
              res_local,res_local,mg%pq1f,mg%pq2f,mg%num_slabs12_pq))

         if (pub_on_root .and. pub_output_detail == VERBOSE) then
            write(stdout,'(a,i9,a,e13.7,a,e12.6,a,e15.9,a)') 'MG ',counter, &
                 '    ', error_norm, '   ',defect_norm, '       ',res_norm, &
                 '     <-- D I'
            call utils_flush(stdout,.true.)
         end if

         call utils_assert(error_norm < 1D6, 'The error norm is absurdly &
              &high. '//trim(explanation),error_norm)
         call utils_assert(error_norm < 1D6, 'The defect norm is absurdly &
              &high. '//trim(explanation),defect_norm)
         call utils_assert(error_norm < 1D6, 'The residual norm is absurdly &
              &high. '//trim(explanation),res_norm)

         call utils_assert(counter < pub_is_multigrid_max_iters, &
              'The multigrid defect-correction solver failed to converge within&
              & the prescribed number of iterations. '//trim(explanation))

         ! hhh: Evaluate stop criteria
         if (error_norm < pub_is_multigrid_error_tol) converged=.true.

       end do

       if (pub_on_root .and. pub_output_detail == VERBOSE) then
         write(stdout,'(a)') 'MG Defect correct solution obtained.'
         write(stdout,'(a,e12.4)') 'MG  - final defect norm: ',defect_norm
         write(stdout,'(a,e12.4)') 'MG  - final error norm : ',error_norm
         write(stdout,'(a,i0)') 'MG  - defect iterations:   ',counter
       end if

       ! jd: Clean up
       if(.not. present(res) .and. pub_is_discretization_order /= 2) then
          deallocate(res_local,stat=ierr)
          call utils_dealloc_check(myself,'res_local', ierr)
       end if
       if(.not. present(eps_full)) then
          deallocate(eps_full_1,stat=ierr)
          call utils_dealloc_check(myself,'eps_full_1', ierr)
       end if
       deallocate(error,stat=ierr)
       call utils_dealloc_check(myself,'error',ierr)
       deallocate(defect,stat=ierr)
       call utils_dealloc_check(myself,'defect',ierr)

    end if ! if defect correction used

    if(multigrid_debug) then
       call visual_scalarfield(pot,pub_fine_grid, &
            'Total potential','_pot_in_mg-dcsolver')
       if(present(res)) then
          call visual_scalarfield(res,pub_fine_grid, &
               'Error equation residual','_res_in_mg-dcsolver')
       end if
       call utils_assert(.false.,'Intentionally terminating ONETEP after MG &
            &debug dumps have been produced.')
    end if

    call utils_trace_out(myself)
    call timer_clock(myself,2)

  end subroutine multigrid_defect_corr_solver
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function multigrid_electrostatic_energy(phi,eps)
    !=========================================================================!
    ! Calculates the electrostatic energy due to a potential phi in a dielec- !
    ! tric with permittivity eps, via (4) in [3], i.e.                        !
    ! E_es = 1/8pi \int eps(r) * (grad phi(r))^2 dr.                          !
    ! This formula is only valid when integrating over all space, in a finite !
    ! box this is only an approximation, rather accurate for neutral molecules!
    ! and quite crude for charged species.                                    !
    ! NB: This value is never used in the calculation, it is only displayed   !
    !     to the user, so that the discretization error can be assessed.      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   phi,               intent=in,  the potential                          !
    !   eps,               intent=in,  the dielectric permittivity            !
    !   ... both in the usual distributed representation.                     !
    ! Returns:                                                                !
    !   calculated electrostatic energy.                                      !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 13/12/2010.                                  !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use constants, only: DP, PI
    use finite_differences, only: finite_difference_mod_grad_sqd
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_is_discretization_order
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in)           :: phi(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional :: eps(mg%ld1f,mg%ld2f,mg%max_slabs12)

    ! jd: Local variables
    real(kind=DP), allocatable :: grad_phi_sqd(:,:,:) ! jd: (\nabla \phi)^2
    real(kind=DP) :: integral                         ! jd: see below
    integer :: ierr                                   ! jd: Error flag

    !------------------------------------------------------------------------

    call utils_trace_in('multigrid_electrostatic_energy')

    ! jd: Break, if multigrid not yet initialised. It's too late to initialise
    !     now, since the initialisation takes care of proper array dimensioning.
    call utils_assert(multigrid_initialised, &
         'Must call multigrid_initialise() first.')

    allocate(grad_phi_sqd(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
    call utils_alloc_check('multigrid_electrostatic_energy','grad_phi_sqd',ierr)

    ! jd: Compute (\nabla \phi)^2
    call finite_difference_mod_grad_sqd(grad_phi_sqd, phi, &
         pub_is_discretization_order,mg)

    ! jd: Compute \int (eps * (\nabla \phi)^2) dr
    if(present(eps)) then
       integral = integrals_product_on_grid(pub_fine_grid, &
            grad_phi_sqd,eps,mg%pq1f, mg%pq2f, mg%num_slabs12_pq)
    else
       integral = integrals_product_on_grid(pub_fine_grid, &
            grad_phi_sqd, m1 = mg%pq1f, m2 = mg%pq2f, m3 = mg%num_slabs12_pq)
    end if

    ! jd: Return 1/8pi times the above integral, cf (4) in [3].
    multigrid_electrostatic_energy = 1.0_DP/(8.0_DP * pi) * integral

    deallocate(grad_phi_sqd,stat=ierr)
    call utils_dealloc_check('multigrid_electrostatic_energy','grad_phi_sqd', &
         ierr)

    call utils_trace_out('multigrid_electrostatic_energy')

  end function multigrid_electrostatic_energy
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   real(kind=DP) function multigrid_finite_box_correction(phi,eps)
    !=========================================================================!
    ! Calculates the correction to the electrostatic energy returned by       !
    ! multigrid_electrostatic_energy() resulting from the finiteness of the   !
    ! box, i.e. the term -1/8pi * surfint eps(r) phi(r) grad_phi(r) dS.       !
    ! NB: This value is never used in the calculation, it is only displayed   !
    !     to the user, so that the discretization error can be assessed.      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   phi,               intent=in,  the potential                          !
    !   eps,               intent=in,  the dielectric permittivity            !
    !   ... both in the usual distributed representation.                     !
    ! Returns:                                                                !
    !   calculated correction.                                                !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 13/12/2010.                                  !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: pub_my_node_id, pub_total_num_nodes
    use constants, only: DP, PI
    use finite_differences, only: finite_difference_gradient
    use integrals, only: integrals_product_on_grid
    use rundat, only: pub_is_discretization_order
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_trace_in, utils_trace_out

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(in)           :: phi(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional :: eps(mg%ld1f,mg%ld2f,mg%max_slabs12)

    ! jd: Local variables
    real(kind=DP), allocatable :: grad_phi(:,:,:,:)       ! jd: grad phi
    real(kind=DP), allocatable :: phi_gradphi_ds(:,:,:)   ! jd: phi*grad_phi*dS
    real(kind=DP) :: integral                             ! jd: see below
    real(kind=DP) :: gradphi_ds                           ! jd: grad_phi * dS
    real(kind=DP) :: surf_weight_x, surf_weight_y, surf_weight_z ! jd: weights
    integer :: i1, i2, islab12                            ! jd: indices
    integer :: ierr                                       ! jd: Error flag

    !------------------------------------------------------------------------

    call utils_trace_in('multigrid_finite_box_correction')

    ! jd: Break, if multigrid not yet initialised. It's too late to initialise
    !     now, since the initialisation takes care of proper array dimensioning.
    call utils_assert(multigrid_initialised, &
         'Must call multigrid_initialise() first.')

    ! jd: Allocate work arrays
    allocate(grad_phi(mg%ld1f,mg%ld2f,mg%max_slabs12,3),stat=ierr)
    call utils_alloc_check('multigrid_finite_box_correction','grad_phi',ierr)
    allocate(phi_gradphi_ds(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
    call utils_alloc_check('multigrid_finite_box_correction', &
         'phi_gradphi_ds',ierr)

    ! jd: Compute grad \phi
    call finite_difference_gradient(grad_phi,phi,pub_is_discretization_order,mg)

    ! jd: Weights for dS -- the grid is orthorombic, but not necessarily cubic
    surf_weight_x = pub_fine_grid%da2%y*pub_fine_grid%da3%z
    surf_weight_y = pub_fine_grid%da1%x*pub_fine_grid%da3%z
    surf_weight_z = pub_fine_grid%da1%x*pub_fine_grid%da2%y

    ! jd: Compute phi * grad_phi * dS on the surfaces of the box and
    !     a lot of zeros inside the box, will allow easy integration
    !     with the volume integrator.
    do islab12=1, mg%num_slabs12_pq
       do i2=1, mg%pq2f
          do i1=1, mg%pq1f
             gradphi_ds = 0D0
             if(i1 == 1) gradphi_ds = &
                  gradphi_ds - grad_phi(i1,i2,islab12,1) * surf_weight_x
             if(i1 == mg%pq1f) gradphi_ds = &
                  gradphi_ds + grad_phi(i1,i2,islab12,1) * surf_weight_x
             if(i2 == 1) gradphi_ds = &
                  gradphi_ds - grad_phi(i1,i2,islab12,2) * surf_weight_y
             if(i2 == mg%pq2f) gradphi_ds = &
                  gradphi_ds + grad_phi(i1,i2,islab12,2) * surf_weight_y
             if(islab12 == 1 .and. pub_my_node_id == 0) gradphi_ds = &
                  gradphi_ds - grad_phi(i1,i2,islab12,3) * surf_weight_z
             if(islab12 == mg%num_slabs12_pq .and. &
                  pub_my_node_id == pub_total_num_nodes-1) gradphi_ds = &
                  gradphi_ds + grad_phi(i1,i2,islab12,3) * surf_weight_z

             phi_gradphi_ds(i1,i2,islab12) = phi(i1,i2,islab12) * gradphi_ds
          end do
       end do
    end do

    ! jd: Compute \int (eps * gphi * (grad \phi) dS)
    if(present(eps)) then
       integral = integrals_product_on_grid(pub_fine_grid, &
            phi_gradphi_ds,eps,mg%pq1f, mg%pq2f, mg%num_slabs12_pq)
    else
       integral = integrals_product_on_grid(pub_fine_grid, &
            phi_gradphi_ds, m1 = mg%pq1f, m2 = mg%pq2f, m3 = mg%num_slabs12_pq)
    end if

    ! jd: Take -1/8pi into account, also remove the dV weight that
    !     integrals_product_on_grid() adds, dS weights have already been used.
    multigrid_finite_box_correction = &
         -1.0_DP/(8.0_DP * pi) * integral / pub_fine_grid%weight

    ! jd: Clean up
    deallocate(phi_gradphi_ds,stat=ierr)
    call utils_dealloc_check('multigrid_finite_box_correction', &
         'phi_gradphi_ds', ierr)
    deallocate(grad_phi,stat=ierr)
    call utils_dealloc_check('multigrid_finite_box_correction','grad_phi_sqd', &
         ierr)

    call utils_trace_out('multigrid_finite_box_correction')

  end function multigrid_finite_box_correction
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine multigrid_prepare_bound_cond(bound,rho_in,uniform_eps)
    !=========================================================================!
    ! Calculates the Coulombic boundary conditions for a potential due to a   !
    ! density rho, on the surface points of the cell. Points in the bulk are  !
    ! zeroed, even though the MG solver never uses them. Coarse-graining of   !
    ! charge is used to keep the O(L^5) scaling in check.                     !
    ! In solvent an approximation is used in which it is assumed that the     !
    ! whole cell is filled with uniform dielectric, with a permittivity of    !
    ! uniform_eps. This will become inaccurate for solvents with low permitti-!
    ! vities. If uniform_eps is omitted, 1.0 is assumed.                      !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   bound,       intent=out,  calculated Coulombic potential on the surf. !
    !   rho_in,      intent=in,   charge density generating the BCs           !
    !   ... both in the usual distributed representation.                     !
    !   uniform_eps, intent=in,   (optional) the uniform diel. permittivity   !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 13/12/2010.                                  !
    ! Revised by Jacek Dziedzic, 15/04/2011 to allow blocks spanning core     !
    !                                       boundaries.                       !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid, &
         cell_grid_extract_box, cell_grid_deposit_box
    use comms, only: pub_my_node_id, pub_total_num_nodes, &
         comms_bcast, comms_barrier, comms_reduce, pub_on_root, pub_root_node_id
    use constants, only: VERBOSE, stdout
    use geometry, only: POINT, magnitude, OPERATOR(-), OPERATOR(*), OPERATOR(+)
    use rundat, only: pub_is_bc_coarseness, pub_is_bc_surface_coarseness, &
         pub_is_bc_threshold, pub_output_detail, pub_charge
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_trace_in, utils_trace_out, utils_assert, utils_abort
    use visual, only: visual_scalarfield

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out)  :: bound(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in)   :: rho_in(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional :: uniform_eps

    ! jd: Internal variables
    real(kind=DP) :: inv_uniform_eps   ! 1/uniform_eps
    integer :: i1, i2, i3, islab12     ! jd: Indices
    integer :: j_packed, n_packed, n_packed_max ! jd: Refer to the packed arrays

    ! jd: Packed coarse-grained charges
    real(kind=DP), allocatable, dimension(:)      :: packed_q
    ! jd: ... and their locations
    real(kind=DP), allocatable, dimension(:,:)    :: packed_loc
    ! jd: Stores a one-block-thick layer of rho_in, extracted from other nodes
    real(kind=DP), allocatable, dimension(:,:,:)  :: rho
    ! jd: Buffer for cell_grid_extract_box, cell_grid_deposit_box
    real(kind=DP), allocatable, dimension(:,:,:)  :: buffer
    ! jd: Value of charge per gridpoint below which all charges are ignored
    real(kind=DP) :: q_threshold
    real(kind=DP) :: q_i
    real(kind=DP) :: q_factor
    real(kind=DP) :: v_here
    type(POINT) :: r_i, r_j, r_of_block, r_0 ! jd: Various vectors
    integer :: n_packed_first, n_packed_last ! jd: CG'ed charges on this node
    integer :: my_portion       ! jd: # of CG'ed charges on this node
    integer :: typical_portion  ! jd: # of CG'ed charges on every node but last
    integer :: tail             ! jd: Difference between the two, on last node
    real(kind=DP) :: d          ! jd: Distance
    real(kind=DP) :: q_sum      ! jd: Charge a in CG'ed block
    real(kind=DP) :: q_in_block ! jd: Sum of the above
    integer :: n_in_block       ! jd: Over that many elements
    integer :: bsize, bheight   ! jd: Desired and actual block height
    integer :: nb1, nb2, nb3    ! jd: Number of blocks along x, y, z
    integer :: b1, b2, b3       ! jd: Block indices
    integer :: i1min, i1max, i2min, i2max, i3min, i3max ! jd: Bounds
    real(kind=DP) :: q ! jd: Current coarse-grained charge point

    integer :: face    ! jd: Face of the cell under consideration (1-6).
    real(kind=DP), allocatable, dimension(:,:) :: face_pot ! jd: Calculated q'ty
    integer :: i1_pq, i1_ld, i2_pq, i2_ld, i3_pq ! jd: Bounds of generic indices
    type(POINT) :: i1_da, i2_da, i3_da ! Grid widths along each generic index
    ! jd: These refer to Cartesian (non-generic) representation
    integer :: width_x, width_y, width_z, start_x, start_y, start_z, ld_x, ld_y

    ! jd: The following refer to surface CG'ing (stage 3).
    integer :: prev1, prev2, next1, next2
    integer :: offset1, offset2
    real(kind=DP) :: t, u
    real(kind=DP) :: width1, width2
    real(kind=DP) :: val_prev_t_prev_u, val_prev_t_next_u, val_next_t_prev_u, &
         val_next_t_next_u

    ! jd: Varia
    integer :: ierr ! jd: Error flag
    character(len=*), parameter :: myself = 'multigrid_prepare_bound_cond'
    logical :: ringing_warning_issued = .false.

    !------------------------------------------------------------------------

    call timer_clock(myself,1)
    call utils_trace_in(myself)

    ! jd: Break, if multigrid not yet initialised. It's too late to initialise
    !     now, since the initialisation takes care of proper array dimensioning.
    call utils_assert(multigrid_initialised, &
         'Must call multigrid_initialise() first.')

    if(multigrid_debug) then
       call visual_scalarfield(rho_in,pub_fine_grid, &
             'Charge density','_rhoin_in_prepare_bound_cond')
    end if

    bound = 0.0_DP

    ! jd: pub_is_bc_coarseness of 0 means zero BC's
    if(pub_is_bc_coarseness == 0) then
      call utils_trace_out(myself)
      call timer_clock(myself,2)
      return
    end if

    if(present(uniform_eps)) then
       inv_uniform_eps = 1.0_DP/uniform_eps
    else
       inv_uniform_eps = 1.0_DP
    end if

    ! jd: Determine the number of blocks along each axis, globally
    bsize = pub_is_bc_coarseness
    nb3 = mg%pq3f/bsize
    nb2 = mg%pq2f/bsize
    nb1 = mg%pq1f/bsize

    ! jd: Count partial blocks as well
    if(mod(mg%pq3f,bsize) /= 0) nb3 = nb3 + 1
    if(mod(mg%pq2f,bsize) /= 0) nb2 = nb2 + 1
    if(mod(mg%pq1f,bsize) /= 0) nb1 = nb1 + 1

    ! jd: Allocate temporaries
    n_packed_max = nb1*nb2*nb3
    allocate(packed_q(n_packed_max),stat=ierr)
    call utils_alloc_check(myself,'packed_q',ierr)
    allocate(packed_loc(3,n_packed_max),stat=ierr)
    call utils_alloc_check(myself,'packed_loc',ierr)
    allocate(rho(mg%ld1f,mg%ld2f,bsize),stat=ierr)
    call utils_alloc_check(myself,'rho',ierr)
    allocate(buffer(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'buffer_fine',ierr)

    ! ---------------------------------------------------------------------
    ! STAGE 1
    ! jd: Go over rho_in in layers, one-block thick.
    !     Get the layer to root node, regardless of who owns the slabs.
    !     On the root node, go over all blocks in the layer, examine all
    !     entries in the block that are above a tiny threshold.
    !     Each entry in packed_q represents the total charge in the block.
    !     Each entry in packed_loc is the location of the centre-of-charge.
    !     A value of -1 in packed_loc signifies end of the array,
    !     beyond which the entries in packed_q are undefined
    !     This stage is very cheap (compared to stage 2), so we can get away
    !     with doing this on the root node.
    ! ---------------------------------------------------------------------

    packed_loc = 1D33     ! jd: Garbage
    ! jd: All grid points with |charge| lower than that will be ignored.
    q_threshold = pub_is_bc_threshold
    n_packed = 0
    q_sum = 0.0_DP

    ! jd: Go over layers, each one-block thick
    do b3 = 1, nb3

       ! jd: Height of the block in this layer
       !     (=bsize except in last global layer)
       if(b3 < nb3) then
          bheight = bsize
       else
          bheight = mod(mg%pq3f,bsize)
       end if
       ! jd: Extract bheight slabs, mg%pq1f x mg%pq2f wide into rho,
       !     regardless of who owns them
       rho = 0D0
       call cell_grid_extract_box(rho, buffer, rho_in, pub_fine_grid, &
            mg%ld1f, mg%ld2f, bheight, mg%ld1f, mg%ld2f, &
            1, 1, (b3-1)*bsize+1, &
            pub_on_root, .false.)

       ! jd: Go over blocks in the current layer, packing them
       if(pub_on_root) then

         i3min = (b3-1)*bsize+1
         i3max = b3*bsize
         if(i3max > mg%pq3f) i3max = mg%pq3f

         do b2=1, nb2
            i2min = (b2-1)*bsize+1
            i2max = b2*bsize
            if(i2max > mg%pq2f) i2max = mg%pq2f

            do b1=1, nb1
               i1min = (b1-1)*bsize+1
               i1max = b1*bsize
               if(i1max > mg%pq1f) i1max = mg%pq1f

               ! jd: Centre of the block
               r_0 = 0.5_DP * ( &
                    real((i1min-1),kind=DP) * pub_fine_grid%da1 + &
                    real((i2min-1),kind=DP) * pub_fine_grid%da2 + &
                    real((i3min-1),kind=DP) * pub_fine_grid%da3 + &
                    real((i1max-1),kind=DP) * pub_fine_grid%da1 + &
                    real((i2max-1),kind=DP) * pub_fine_grid%da2 + &
                    real((i3max-1),kind=DP) * pub_fine_grid%da3 )

               n_in_block = 0
               q_in_block = 0.0_DP
               r_of_block%X = 0.0_DP
               r_of_block%Y = 0.0_DP
               r_of_block%Z = 0.0_DP

               ! jd: Go over points in current block
               do islab12=1, bheight
                  i3 = i3min + islab12 - 1
                  do i2=i2min, i2max
                     do i1=i1min, i1max

                        call utils_assert(i1<=mg%pq1f .and. i1>0, &
                             'Bound check [1]')
                        call utils_assert(i2<=mg%pq2f .and. i2>0, &
                             'Bound check [2]')
                        call utils_assert(i3<=mg%pq3f .and. i3>0, &
                             'Bound check [3]')
                        q_i = rho(i1,i2,islab12) * pub_fine_grid%weight

                        if(abs(q_i) > q_threshold) then
                           r_i = &
                                real((i1-1),kind=DP) * pub_fine_grid%da1 + &
                                real((i2-1),kind=DP) * pub_fine_grid%da2 + &
                                real((i3-1),kind=DP) * pub_fine_grid%da3

                           q_in_block = q_in_block + q_i
                           r_of_block = r_of_block + q_i* (r_i-r_0)
                           n_in_block = n_in_block + 1
                        end if

                     end do  ! }
                  end do     ! } over gridpoints in block
               end do        ! }

               ! jd: Store packed block
               if(n_in_block /= 0) then
                  n_packed = n_packed + 1
                  packed_q(n_packed) = q_in_block
                  q_sum = q_sum + q_in_block
                  r_of_block = 1.0_DP/q_in_block * r_of_block + r_0
                  packed_loc(1,n_packed) = r_of_block%X
                  packed_loc(2,n_packed) = r_of_block%Y
                  packed_loc(3,n_packed) = r_of_block%Z
               end if

            end do ! }
         end do    ! } over blocks
      end if ! if on root
    end do  ! over layers

    if(pub_on_root .and. pub_output_detail == VERBOSE) then
       write(stdout,'(a,f0.7,a,i0,a)') 'MG boundary conditions: total CG''ed &
            &q: ', q_sum,' over ',n_packed,' coarse points'
    end if

    if(pub_on_root) then
       call utils_assert(abs(q_sum - real(nint(q_sum),kind=DP)) < 0.05_DP, &
            'The coarse-grained charge is suspiciously far from an integer valu&
            &e. First check that your charge is correctly localized within the &
            &FD grid. Allow up to 10 bohr of extra vacuum to account for the ri&
            &nging, especially if fine_grid_scale is different from 2.0. If you&
            & get this error during a properties calculation, make sure that &
            &the calculation is well-converged first. If none of the above appl&
            &y, reduce is_bc_threshold (to, say, 1E-10) to make charge coarse-&
            &graining more accurate. If this doesn''t help, this &
            &error might indicate that you need to decrease ''is_bc_coarseness&
            &'' to improve accuracy. The obtained coarse-grained charge was', &
            q_sum)

       if(n_packed == 0) then
          write(stdout,'(/a)') &
               'WARNING: Charge density is zero (or extremely low) everywhere!'
          return ! jd: Leaving bound(:)=0D0, which is fine
       end if
    end if

    ! jd: Because we killed off all charge below the threshold, we've left out
    !     the ringing. Better adjust the CGed charge by multiplying it by
    !     q_desired / q_cged (a factor in the order of 1.00001) so that the
    !     the total charge agrees.
    q_factor = real(nint(q_sum),kind=DP) / q_sum
    if(q_factor /= 0.0_DP) then ! jd: But not when the molecule is neutral
       packed_q(1:n_packed) = packed_q(1:n_packed) * q_factor
    end if

    ! jd: These are no longer needed
    deallocate(buffer,stat=ierr)
    call utils_dealloc_check(myself,'buffer',ierr)
    deallocate(rho,stat=ierr)
    call utils_dealloc_check(myself,'rho',ierr)

    ! jd: Broadcast the packed arrays to everyone
    call comms_bcast(pub_root_node_id,packed_q)
    call comms_bcast(pub_root_node_id,packed_loc)
    call comms_bcast(pub_root_node_id,n_packed)

    ! ---------------------------------------------------------------------
    ! STAGE 2
    ! jd: Go over packed array and evaluate potential on coarse surface points.
    !     Then fill in the missing points by bilinear interpolation.
    !     Each processor deals with *all* the points on the boundary, but only
    !     a subset of the packed point charges. It doesn't make sense to use the
    !     usual slab distribution, because then the first and last nodes have
    !     much more work -- they deal with top and bottom surfaces. Also the
    !     surface coarse-graining would become more involved (as surface blocks
    !     would cross processor boundaries).
    ! ---------------------------------------------------------------------

    ! jd: Determine how big our portion of the packed array is
    !     When n_packed is not divisible by n_nodes, give the remainder to
    !     the last node -- we have n_packed in the order of 10^3-10^6, so
    !     the extra tail is negligible.
    my_portion = n_packed / pub_total_num_nodes
    typical_portion = my_portion
    tail = mod(n_packed,pub_total_num_nodes)
    if(pub_my_node_id == pub_total_num_nodes-1) my_portion = my_portion + tail
    n_packed_first = pub_my_node_id * typical_portion + 1
    n_packed_last = n_packed_first + my_portion - 1

    bsize = pub_is_bc_surface_coarseness

    ! --------------------------------------------------
    ! --- Go over the 6 faces, fill them and deposit ---
    ! --------------------------------------------------
    ! jd: The face is covered by indices i1, i2.

    allocate(buffer(pub_fine_grid%ld1,pub_fine_grid%ld2, &
         pub_fine_grid%max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'face_pot',ierr)

    do face = 1, 6
       ! jd: Silence compiler warning
       i1_ld = 0; i1_pq = 0; i1_ld = 0; i2_ld = 0; i2_pq = 0; i3_pq = 0
       ! jd: Top or bottom face (XY): i1=x, i2=y, i3=z
       if(face == 1 .or. face == 2) then
          i1_ld = pub_fine_grid%ld1
          i1_pq = mg%pq1f
          i1_da = pub_fine_grid%da1
          i2_ld = pub_fine_grid%ld2
          i2_pq = mg%pq2f
          i2_da = pub_fine_grid%da2
          i3_pq = mg%pq3f
          i3_da = pub_fine_grid%da3
          ld_x = i1_ld
          ld_y = i2_ld
          width_x = i1_pq
          width_y = i2_pq
          width_z = 1
          start_x = 1
          start_y = 1
          if(face == 1 ) then
             start_z = 1           ! top
          else
             start_z = i3_pq       ! bottom
          end if
       end if
       ! jd: Left or right face (YZ): i1=y, i2=z, i3=x
       if(face == 3 .or. face == 4) then
          i1_ld = pub_fine_grid%ld2
          i1_pq = mg%pq2f
          i1_da = pub_fine_grid%da2
          i2_ld = pub_fine_grid%ld3
          i2_pq = mg%pq3f
          i2_da = pub_fine_grid%da3
          i3_pq = mg%pq1f
          i3_da = pub_fine_grid%da1
          ld_x = 1
          ld_y = i1_ld
          width_x = 1
          width_y = i1_pq
          width_z = i2_pq
          start_y = 1
          start_z = 1
          if(face == 3 ) then
             start_x = 1           ! left
          else
             start_x = i3_pq       ! right
          end if
       end if
       ! jd: Front or back (XZ): i1=x, i2=z, i3=y
       if(face == 5 .or. face == 6) then
          i1_ld = pub_fine_grid%ld1
          i1_pq = mg%pq1f
          i1_da = pub_fine_grid%da1
          i2_ld = pub_fine_grid%ld3
          i2_pq = mg%pq3f
          i2_da = pub_fine_grid%da3
          i3_pq = mg%pq2f
          i3_da = pub_fine_grid%da2
          ld_x = i1_ld
          ld_y = 1
          width_x = i1_pq
          width_y = 1
          width_z = i2_pq
          start_x = 1
          start_z = 1
          if(face == 5 ) then
             start_y = 1           ! front
          else
             start_y = i3_pq       ! back
          end if
       end if

       ! jd: Allocate and zero storage for this cell face
       allocate(face_pot(i1_ld,i2_ld),stat=ierr)
       call utils_alloc_check(myself,'face_pot',ierr)
       face_pot = 0D0

       ! jd: Low-index-value or high-index-value faces?
       if(face == 1 .or. face == 3 .or. face == 5) then
          i3 = 1     ! jd: Face location: top, left or front
       else
          i3 = i3_pq ! jd: Face location: bottom, right or back
       end if

       ! jd: Loop over points on current face.
       !     Everything's been cleared earlier, so no need to worry about
       !     padding margins.
       do i2=1, i2_pq
          do i1=1, i1_pq
             ! jd: If coarse-graining surfaces, skip inner points
             if(bsize /=1 .and. &
                  ((mod(i1,bsize) /= 1 .and. i1 /= i1_pq) .or. &
                   (mod(i2,bsize) /= 1 .and. i2 /= i2_pq))) cycle

             r_i = &
                  real((i1-1),kind=DP) * i1_da + &
                  real((i2-1),kind=DP) * i2_da + &
                  real((i3-1),kind=DP) * i3_da

             v_here = 0.0_DP

             ! jd: Go over all my points in the packed array
             do j_packed = n_packed_first, n_packed_last
                r_j%X = packed_loc(1,j_packed)
                r_j%Y = packed_loc(2,j_packed)
                r_j%Z = packed_loc(3,j_packed)

                q = packed_q(j_packed)
                d = magnitude(r_i-r_j)

                ! jd: Packed charge points very close to the face indicate
                !     something's seriously wrong.
                if(abs(d) < 1D-2) then
                   if(.not. ringing_warning_issued) then
                      write(stdout,'(a)') 'WARNING: Non-zero charge at, or very&
                           & close to, one of the simulation cell faces.'
                      write(stdout,'(a,f0.5,a,f0.5,a,f0.5,a)') '     &
                           &    Charge at:     ',r_j%X,',',r_j%Y,',',r_j%Z,'.'
                      write(stdout,'(a,f0.5,a,f0.5,a,f0.5,a)') '     &
                           &    Point on face: ',r_i%X,',',r_i%Y,',',r_i%Z,'.'
                      write(stdout,'(a,e16.5)') &
                           'WARNING: Charge:                 ',q
                      write(stdout,'(a,e16.5)') &
                           'WARNING: Distance from box face: ',d
                   end if
                   if(abs(q) > 1D-7) then
                      call utils_abort('Non-zero charge density&
                           & detected on or very close to cell boundaries.&
                           & Check your cell bounds and charge localization.')
                   else
                      ! jd: Unless it's just the ringing, then ignore it
                      if(.not. ringing_warning_issued) then
                         write(stdout,'(a)') 'WARNING: Density ringing detected on the&
                              & face of the simulation cell.'
                         write(stdout,'(a)') '         This may lead to&
                              & inaccuracies in the charge coarse-graining&
                              & procedure.'
                         write(stdout,'(a)') '         Consider making your &
                              &simulation cell slightly larger.'
                         ringing_warning_issued = .true.
                      end if
                      ! jd: Avoid the singularity
                      q = 0D0
                      d = 1D0
                   end if
                end if
                v_here = v_here + q / d
             end do
             v_here = v_here * inv_uniform_eps
             face_pot(i1,i2) = v_here

          end do  ! } over points on this face
       end do     ! }

       ! ---------------------------------------------------------------------
       ! STAGE 3
       ! jd: Now fill in the missing points on the surfaces by bilinear
       !     interpolation
       ! ---------------------------------------------------------------------

       if(bsize /= 1) then ! jd: No need to do this if all points filled already
          ! jd: Loop over points on current face.
          do i2=1, i2_pq
             do i1=1, i1_pq

                ! jd: Ignore points that *are* on the coarsened surface
                if(  (i1==i1_pq .or. mod(i1,bsize) == 1) .and. &
                     (i2==i2_pq .or. mod(i2,bsize) == 1)) then
                     if(face_pot(i1,i2) == 0D0) then
                        call utils_abort("Internal error [1] in "//myself)
                     end if
                     cycle
                end if

                r_i = &
                     real((i1-1),kind=DP) * i1_da + &
                     real((i2-1),kind=DP) * i2_da + &
                     real((i3-1),kind=DP) * i3_da

                call utils_assert(face_pot(i1,i2) == 0D0, &
                     'Internal error [2] in '//myself)

                ! jd: Our offset within a CG surface block (0..bsize-1)
                offset1 = mod(i1-1,bsize)
                offset2 = mod(i2-1,bsize)

                ! jd: Determine what the coarse-grained neighbours are
                prev1 = i1 - offset1
                prev2 = i2 - offset2
                next1 = i1 + (bsize-offset1)
                next2 = i2 + (bsize-offset2)
                if(prev1 < 1) prev1 = 1
                if(prev2 < 1) prev2 = 1
                if(next1 > i1_pq) next1 = i1_pq
                if(next2 > i2_pq) next2 = i2_pq

                ! jd: Determine the size of the bilinear interpolation block
                !     This is usually bsize, except close to pq boundaries
                width1=real(next1-prev1,kind=DP)
                width2=real(next2-prev2,kind=DP)

                ! jd: Take care of the corner case where prev==next.
                !     This can happen e.g. when i1==pq1, and prev is i1,
                !     but simultaneously next is i1, because it cannot lie
                !     outside ot the boundary. We then want 't' or 'u' to
                !     be zero, not NaN.
                if(width1==0.0_DP) width1=1.0_DP
                if(width2==0.0_DP) width2=1.0_DP

                ! jd: Pick right neighbours
                val_prev_t_prev_u = face_pot(prev1,prev2)
                val_prev_t_next_u = face_pot(prev1,next2)
                val_next_t_prev_u = face_pot(next1,prev2)
                val_next_t_next_u = face_pot(next1,next2)
                t=real(offset1,kind=DP)/width1
                u=real(offset2,kind=DP)/width2

                ! jd: Apply bilinear interpolation
                v_here = &
                     (1.0_DP-t)*(1.0_DP-u)*val_prev_t_prev_u + &
                     (1.0_DP-t)*u*val_prev_t_next_u + &
                     t*(1.0_DP-u)*val_next_t_prev_u + &
                     t*u*val_next_t_next_u

                face_pot(i1,i2) = v_here

             end do  ! } over points
          end do     ! } on the face
       end if ! jd: If surface coarse-graining

       ! jd: Careful not to double-count edges
       face_pot(1,:) = face_pot(1,:)/2.0
       face_pot(:,1) = face_pot(:,1)/2.0
       face_pot(i1_pq,:) = face_pot(i1_pq,:)/2.0
       face_pot(:,i2_pq) = face_pot(:,i2_pq)/2.0

       ! jd: Ditto for corners
       face_pot(1,1) = face_pot(1,1) * 4.0/3.0
       face_pot(i1_pq,1) = face_pot(i1_pq,1) * 4.0/3.0
       face_pot(1,i2_pq) = face_pot(1,i2_pq) * 4.0/3.0
       face_pot(i1_pq,i2_pq) = face_pot(i1_pq,i2_pq) * 4.0/3.0

       ! jd: Reduce contributions to this face over all nodes to root
       call comms_reduce('SUM',face_pot,root=pub_root_node_id)

       ! jd: Deposit this from root to the slab distrib'n on respective nodes
       call cell_grid_deposit_box(bound, face_pot, buffer, pub_fine_grid, &
            width_x, width_y, width_z, ld_x, ld_y, &
            start_x, start_y, start_z, pub_on_root, .false.)

       deallocate(face_pot,stat=ierr)
       call utils_dealloc_check(myself,'face_pot',ierr)

    end do ! over faces

    deallocate(buffer,stat=ierr)
    call utils_dealloc_check(myself,'buffer',ierr)

    ! jd: Sync to get accurate timings
    call comms_barrier

    ! jd: Clean up
    deallocate(packed_loc,stat=ierr)
    call utils_dealloc_check(myself,'packed_loc',ierr)
    deallocate(packed_q,stat=ierr)
    call utils_dealloc_check(myself,'packed_q',ierr)

    call utils_trace_out(myself)
    call timer_clock(myself,2)

  end subroutine multigrid_prepare_bound_cond
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  real(kind=DP) function multigrid_calculate_hartree(phi, rho, uniform_eps, &
       eps_full, eps_half)
    !=========================================================================!
    ! This is a back-end subroutine for the calculation of the Hartree energy !
    ! and Hartree potential using a multigrid solver. There are three,        !
    ! mutually-exclusive cases when this is in order:                         !
    !   a) An open BC calculation, without smeared ions, in vacuum.           !
    !   b) An open BC calculation, with smeared ions, in vacuum.              !
    !   c) An open BC calculation, with smeared ions, in solvent.             !
    ! In a) the Hartree energy and potential are only due to electrons, and   !
    !       rho is the electronic density.                                    !
    ! In b) and c) the Hartree energy and potential are due to the whole      !
    !       molecule, i.e. electrons and smeared ions, and rho is the total   !
    !       density.                                                          !
    ! In c) the dielectric permittivity must be supplied through eps.         !
    ! However, this is assumed to be a scalar value, we thus assume that in   !
    ! the solvated case the whole cell is filled with bulk-permittivity die-  !
    ! lectric.
    !                                                                         !
    ! Note that because this subroutine cannot know whether the supplied      !
    ! density is electronic or total, it must be spin-unaware, i.e. if the    !
    ! system is spin-polarized, the electronic density must be summed across  !
    ! spins prior to calling this subroutine and the second component of the  !
    ! potential must be initialized from phi after the call.                  !
    !                                                                         !
    ! Because the multigrid only works with the fine grid, there is no need   !
    ! to pass a grid argument to this subroutine.                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   phi,               intent=out, the calculated Hartree potential       !
    !   rho,               intent=in,  the input density, electronic or total !
    !   uniform_eps, (opt) intent=in,  the dielectric permittivity (scalar)   !
    !                                  to use for approximate BCs             !
    !   eps_full, (opt)    intent=in,  the dielectric permittivity (for all   !
    !                                  points)                                !
    !   eps_half, (opt)    intent=in,  the dielectric permittivity (for all   !
    !                                  points, shifted by half a grid point)  !
    ! All of the optional arguments default to 1, or arrays of 1, if omitted. !
    ! Returns:                                                                !
    !   calculated Hartree potential due to rho.                              !
    !-------------------------------------------------------------------------!
    ! Written by Jacek Dziedzic, 10/12/2010.                                  !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: pub_on_root
    use constants, only: stdout
    use integrals, only: integrals_product_on_grid
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_assert, &
         utils_trace_in, utils_trace_out

    ! jd: Arguments
    ! jd: Electronic or total density in real space.
    real(kind=DP), intent(in)  :: rho(mg%ld1f,mg%ld2f,mg%max_slabs12)
    ! jd: Calculated potential due to rho
    real(kind=DP), intent(out) :: phi(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional :: uniform_eps
    real(kind=DP), intent(in), optional &
         :: eps_full(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in), optional &
         :: eps_half(mg%ld1f,mg%ld2f,mg%max_slabs12,3)

    ! jd: Boundary conditions
    real(kind=DP), allocatable :: bound(:,:,:)

    real(kind=DP) :: hartree_from_phi_formula
    real(kind=DP) :: finite_box_correction

    integer :: ierr ! jd: Error flag

    !------------------------------------------------------------------------

    call utils_trace_in('multigrid_calculate_hartree')

    ! jd: Break, if multigrid not yet initialised. It's too late to initialise
    !     now, since the initialisation takes care of proper array dimensioning.
    call utils_assert(multigrid_initialised, &
         'Must call multigrid_initialise() first.')

    ! jd: Allocate the array that will hold the BCs for the MG calculation
    allocate(bound(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
    call utils_alloc_check('multigrid_calculate_hartree','bound',ierr)

    ! jd: Determine the correct boundary conditions
    call multigrid_prepare_bound_cond(bound,rho,uniform_eps)

    ! jd: Use the multigrid solver to get the Hartree potential
    call multigrid_defect_corr_solver(pot = phi, & ! output
         rho = rho, bound = bound, &               ! input
         eps_full = eps_full, eps_half = eps_half) ! input

    ! jd: Calculate E_Hartree, the returned value
    multigrid_calculate_hartree = &
         0.5_DP * integrals_product_on_grid(pub_fine_grid, &
         phi, rho, mg%pq1f, mg%pq2f, mg%num_slabs12_pq)

    if(pub_on_root) write(stdout,'(a,f14.6)') &
         "MG Hartree energy from (\int phi*rho dr):           ", &
         multigrid_calculate_hartree

    ! jd: Calculate 1/8pi \int eps (grad phi)^2, the other estimate
    hartree_from_phi_formula = multigrid_electrostatic_energy(phi,eps_full)

    ! jd: Calculate \surfint phi * eps * grad_phi dS, the correction
    finite_box_correction = multigrid_finite_box_correction(phi,eps_full)

    if(pub_on_root) then
       write(stdout,'(a,f14.6)') &
            "MG Hartree energy from (\int eps (grad phi)^2 dr):  ", &
            hartree_from_phi_formula
       write(stdout,'(a,f14.6,a,f10.5,a)') &
            "MG Discrepancy due to finite box and discretization:", &
            hartree_from_phi_formula - multigrid_calculate_hartree, ' (', &
            100.0_DP*(hartree_from_phi_formula - multigrid_calculate_hartree) /&
            multigrid_calculate_hartree,'%)'
       write(stdout,'(a,f14.6)') &
            "MG Correction to account for finite box:            ", &
            finite_box_correction
       hartree_from_phi_formula = hartree_from_phi_formula+finite_box_correction
       write(stdout,'(a,f14.6,a,f10.5,a)') &
            "MG Remaining discrepancy, due to discretization:    ", &
            hartree_from_phi_formula - multigrid_calculate_hartree, ' (', &
            100.0_DP*(hartree_from_phi_formula - multigrid_calculate_hartree) /&
            multigrid_calculate_hartree,'%)'

    end if

    ! jd: Clean up
    deallocate(bound,stat=ierr)
    call utils_dealloc_check('multigrid_calculate_hartree','bound',ierr)

    call utils_trace_out('multigrid_calculate_hartree')

  end function multigrid_calculate_hartree
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine compute_current_defect(defect,defect_norm,order,source,eps,v_i)
    !=========================================================================!
    ! This subroutine will compute the current defect for the Poisson equation!
    ! in a dielectric medium.  The defect is defined as the amount that the   !
    ! current approximation fails to satisfy the equation or in other words:  !
    !                                                                         !
    !    defect = source + \nabla \cdot [eps \grad v_i]                       !
    !                                                                         !
    ! The discretization of the operators is selected by the                  !
    ! discretization_order parameter.  This subroutine is intended for use    !
    ! within the defect correction method which will produce a higher         !
    ! accuracy solution than is possible with normal multigrid iteration.     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !    defect      (output): defined as above                               !
    !    defect_norm (output): discrete L2 norm of the defect function        !
    !    order       (input): finite difference order to use                  !
    !    source      (input): generally the 4\pi\rho or RHS of the Poisson eq !
    !    eps         (input): the diel. functional on the full grid points    !
    !    v_i         (input): the current approximation of the potential      !
    !                                                                         !
    ! Arrays are assumed to be in the distributed slab representation.        !
    !-------------------------------------------------------------------------!
    ! Written by Hatem H Helal 30/04/2010                                     !
    ! Mended for onetep by HH Helal 7/6/2010                                  !
    ! Adapted for parallel operation by Jacek Dziedzic, 10-11/06/2010.        !
    !=========================================================================!

    use cell_grid, only: pub_fine_grid
    use finite_differences, only: finite_difference_gradient, &
         finite_difference_laplacian
    use integrals, only: integrals_product_on_grid
    use utils, only: utils_alloc_check, utils_dealloc_check, &
         utils_sanity_check, utils_trace_in, utils_trace_out
    use rundat, only: pub_is_bulk_permittivity

    implicit none

    ! jd: Arguments
    real(kind=DP), intent(out) :: defect(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(out) :: defect_norm

    integer,       intent(in)  :: order
    real(kind=DP), intent(in)  :: source(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in)  :: eps(mg%ld1f,mg%ld2f,mg%max_slabs12)
    real(kind=DP), intent(in)  :: v_i(mg%ld1f,mg%ld2f,mg%max_slabs12)

    ! jd: Internal variables
    real(kind=DP), allocatable     :: lap_v(:,:,:)         ! nabla^2 v
    real(kind=DP), allocatable     :: grad_eps(:,:,:,:)    ! nabla eps
    real(kind=DP), allocatable     :: grad_v(:,:,:,:)      ! nabla v
    integer                        :: ierr                 ! jd: Error flag
    character(len=*), parameter    :: myself = 'compute_current_defect'

    !------------------------------------------------------------------------

    call utils_trace_in(myself)

    call utils_sanity_check(source,'source in '//myself)
    call utils_sanity_check(eps,'eps in '//myself)
    call utils_sanity_check(v_i,'v_i in'//myself)

    ! Allocate memory for local variables
    allocate(lap_v(mg%ld1f,mg%ld2f,mg%max_slabs12),stat=ierr)
    call utils_alloc_check(myself,'lap_v',ierr)
    allocate(grad_eps(mg%ld1f,mg%ld2f,mg%max_slabs12,3),stat=ierr)
    call utils_alloc_check(myself,'grad_eps',ierr)
    allocate(grad_v(mg%ld1f,mg%ld2f,mg%max_slabs12,3),stat=ierr)
    call utils_alloc_check(myself,'grad_v',ierr)

    ! Calculate the gradient and Laplacian of the potential
    call finite_difference_gradient(grad_v,v_i,order,mg)
    call finite_difference_laplacian(lap_v,v_i,order,mg)
    ! jd: Calculate the gradient of epsilon
    call finite_difference_gradient(grad_eps,eps,order,mg)

    call utils_sanity_check(lap_v,'lap_v in '//myself)
    call utils_sanity_check(grad_v(:,:,:,1),'grad_v/x in '//myself)
    call utils_sanity_check(grad_v(:,:,:,2),'grad_v/y in '//myself)
    call utils_sanity_check(grad_v(:,:,:,3),'grad_v/z in '//myself)
    call utils_sanity_check(grad_eps(:,:,:,1),'grad_eps/x in '//myself)
    call utils_sanity_check(grad_eps(:,:,:,2),'grad_eps/y in '//myself)
    call utils_sanity_check(grad_eps(:,:,:,3),'grad_eps/z in '//myself)

    ! Now have everything we need to compute the current defect and defect_norm
    defect(:,:,:) = source(:,:,:) + eps(:,:,:) * lap_v(:,:,:) + &
         grad_eps(:,:,:,1) * grad_v(:,:,:,1) + &
         grad_eps(:,:,:,2) * grad_v(:,:,:,2) + &
         grad_eps(:,:,:,3) * grad_v(:,:,:,3)

    defect_norm = sqrt(integrals_product_on_grid(pub_fine_grid, &
         defect, defect, mg%pq1f, mg%pq2f, mg%num_slabs12_pq))

    ! jd: Clean up
    deallocate(grad_v,stat=ierr)
    call utils_dealloc_check(myself,'grad_v',ierr)
    deallocate(grad_eps,stat=ierr)
    call utils_dealloc_check(myself,'grad_eps',ierr)
    deallocate(lap_v,stat=ierr)
    call utils_dealloc_check(myself,'lap_v',ierr)

    call utils_trace_out(myself)

  end subroutine compute_current_defect
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module multigrid_methods

