! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!=============================================================================!
!              G E O M E T R Y _ D E L O C  _ A L G O R                       !
!=============================================================================!
!                                                                             !
! $Id: geometry_deloc_algor.F90,v 1.9 2009/04/22 15:11:53 cks22 Exp $         !
!-----------------------------------------------------------------------------!
!                                                                             !
!         Algorithms for delocalized coordinates optimization                 !
!                                                                             !
!=============================================================================!
!                                                                             !
!                                                                             !
!-----------------------------------------------------------------------------!
! This module contains utilities required for geometry optimization using     !
! delocalized internal coordinates. This technology is based on the DMol3     !
! implementation which is described in:                                       !
! J. Andzelm, R. D. King-Smith, G. Fitzgerald,                                !
! Chem. Phys. Lett. 335 (2001) 321-326                                        !
!                                                                             !
!  This module can locate minima and transition states                        !
!  using delocalized internal coordinates that are generated                  !
!  automatically or with the help of MDF file.                                !
!  It also handles fixed distance, bond angle and dihedral                    !
!  angle constraints and Cartesian coordinates.                               !
!  It can deal with disconnected fragments.                                   !
!  Failure in these routines sets ip_error to 1 and that will                 !
!  make optimization in Cartesians only, input permiting.                     !
!-----------------------------------------------------------------------------!
! Written from  "Module Specification for New Plane Wave Code",  M. Segall,   !
!    P. Lindan, M. Probert, C. Pickard, P. Hasnip, S. Clark and M. Payne      !
!                           Copyright 1999/2000                               !
!-----------------------------------------------------------------------------!
! Written by Victor Milman, v1.0, 20/11/2002                                  !
!-----------------------------------------------------------------------------!
!                                                                             !
! $Log: geometry_deloc_algor.F90,v $
! Revision 1.9  2009/04/22 15:11:53  cks22
!
! cks_2009_04_22: several minor bug fixes and the inclusion of "classical" point charges to the external potential.
!
! Revision 1.8  2009/02/05 03:20:59  cks22
!
! cks_2009_02_05: Hartree-Fock exchange (slow algorithm, using the FFT box equal to the simulation cell). Also a bug fix by Nick Hine in restart_ngwfs_tightbox_input to follow the F95 standard and fix wrong reading of NGWFs with pg compilers.
!
! Revision 1.19  2009/02/02 17:53:01  cks
!
! cks_2009_02_02: Version 2.2.16
!
! Revision 1.7  2008/12/18 12:55:35  ndmh3
!
! ndmh_2008_12_18: Checked in changes provided by Victor Milman
! (predominantly involving changing '*' to 'stdout' throughout - hence
! the large number of files involved. Also fixed various output
! strings, and grd output for very large files, and added some changes
! from Quintin to nho_mod.F90
!
! Revision 1.3  2008/12/01 14:08:43  vmilman
! Merge with the updated ODG version 2.2.14
!
! Revision 1.6  2008/07/31 17:48:23  cks22
!
! cks_2008_07_31: (1) Trivial changes in spherical_mod and (2) incorporation of empirical dispersion forces and energy in order to be able to use dispersion contributions in geometry optimisation.
!
! Revision 1.12  2008/07/31 17:41:16  cks
!
! cks_2008_07_31: VdW forces and geopt.
!
! Revision 1.5  2007/08/19 08:54:53  cks22
!
! cks_2007_08_19: Extension to use simulation cells smaller than the 6 NGWF radii rule, along
! one or more lattice vector directions. The minimum simulation cell is now defined by the twice
! the maximum NGWF radius.
!
! Revision 1.3  2007/08/17 17:17:14  cks
!
! cks_2007_08_17: Merged flex and current versions.
!
! Revision 1.4  2007/01/31 14:10:45  pdh1001
! Victor's changes to TS search and delocalised internals, plus Windows fix.
!
! Revision 1.3  2006/05/17 08:41:21  cks22
!
! cks_2006_05_17: Some cosmetic changes in properties_mod and latest minor
! changes by Victor.
!
! Revision 1.2  2006/04/21 15:42:09  cks22
!
! cks_2004_04_21: Added Victor's transition state search files.
! Also added ngwf_halo, a new experimentlal feature which can take into account
! small but non-zero matrix elements.
!
! Revision 1.1  2006/02/10 16:45:06  pdh1001
! Merged with source from Victor Milman (mainly to incorporate delocalized internals) and some minor bug fixes.
!
! Revision 1.4  2003/10/31 14:59:20  mijp2
! added constraints bug-fix in deloc-internals by Victor
!
! Revision 1.1  2003/04/30 22:17:19  mijp2
! Initial revision
!                                                                             !
!                                                                             !
!=============================================================================!

module geometry_deloc_algor

  use geometry_deloc_utils, only: di_dp, di_pi, map_fixed_coords, &
       bmat_disconn_fragm, coords_cart, trans, elat, gradient_cart, ntrans, &
       number_atoms, di_on_root, di_iprint, di_stdout, length_conv, symflag, &
       ncons, num_fixed_coords, iupdat, num_rigid_body_atoms, np_p, ip_error, &
       num_fixed_constr_rigid_body, ncycle, xc_redressed, eold, lp_mdf, &
       lp_mdf_del, bohr_di, ec, ldokpt, l_scan_pes, hess_cartesian, lp_del, &
       list_fixed_coords, lp_idb, list_rigid_body_atoms, ncons_uns, period, &
       ldokpt, neqatm, atsymb, atomic_numbers, tsflag, maxdiis, lp_discf, &
       global_data_exists, lp_del_test, di_energy_label, use_deloc_int_output,&
       deloc_utils_write_di_data, deloc_utils_io_abort, &
       deloc_utils_unpack_atom_indices, deloc_utils_pack_atom_indices, &
       deloc_utils_read_di_data, deloc_utils_rationalize_coords, &
       deloc_utils_save_structure, deloc_utils_read_opt_mode, &
       deloc_utils_geom_converged, deloc_utils_read_fix_cartesians, &
       deloc_utils_read_hessian, deloc_utils_read_fix_atoms, &
       deloc_utils_read_np_p, deloc_utils_get_alt_bond_list, &
       deloc_utils_write_hessian, deloc_utils_read_di_constraints
 
  implicit none                                 !Impose strong typing

  private                                       !Everything is private ...

  !---------------------------------------------------------------------------!
  !                       P u b l i c   R o u t i n e s                       !
  !---------------------------------------------------------------------------!
  public :: geom_algor_DELOC                    ! del_cor
  !---------------------------------------------------------------------------!
  !                        P u b l i c   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!

  !---------------------------------------------------------------------------!
  !                      P r i v a t e   R o u t i n e s                      !
  !---------------------------------------------------------------------------!

  private :: di_algor_allocate
  private :: di_algor_deallocate
  private :: di_algor_optimize                  ! del_cor2
  private :: di_algor_get_RB_deloc_coords       ! get_del3
  private :: di_algor_get_init_deloc_coords     ! get_del1
  private :: di_algor_add_constrain_to_fixed    ! cpyicc
  private :: di_algor_get_final_deloc_coords    ! get_del2
  private :: di_algor_remove_0_weight_coords    ! tidyupp
  private :: di_algor_remove_0_grad_coords      ! symnicp
  private :: di_algor_get_weight_primitives     ! getwghtp
  private :: di_algor_make_bond_distances       ! bondl2
  private :: di_algor_make_bond_list            ! fbond
  private :: di_algor_combine_bond_lists        ! fbondup
  private :: di_algor_print_bonds               ! fbondpri
  private :: di_algor_make_B_matrix             ! makebpr
  private :: di_algor_make_B_matrix_old         ! makebmd
  private :: di_algor_initialize_H_internal     ! defhedp
  private :: di_algor_invert_B_matrix           ! binvtd
  private :: di_algor_svd_matrix_invert         ! invmatd
  private :: di_algor_U_transform_B_matrix      ! tranbmp
  private :: di_algor_transf_Hess_to_intern     ! hssindp
  private :: di_algor_make_bonds_B_matrix       ! makebp
  private :: di_algor_make_angles_B_matrix      ! makeap
  private :: di_algor_make_torsions_B_matrix    ! makedp
  private :: di_algor_make_list_primitives      ! mklipe
  private :: di_algor_make_F_primitives         ! domig
  private :: di_algor_make_F_angles             ! findang
  private :: di_algor_make_F_torsions           ! finddih
  private :: di_algor_make_F_disc_fragm         ! sepfrag
  private :: di_algor_diagonalize_F_matrix      ! makef_d
  private :: di_algor_get_step                  ! get_step
  private :: di_algor_get_step_unsat_constr     ! get_step_u
  private :: di_algor_write_constr_di_data      ! ficcwr
  private :: di_algor_init_constraints          ! const_p
  private :: di_algor_find_internal_coord       ! find_int
  private :: di_algor_stop_fixed_constraints    ! confix
  private :: di_algor_tidy_constraints          ! tidycon
  private :: di_algor_schmidt_constraints       ! schmidt_p
  private :: di_algor_print_constraint          ! pr_consat
  private :: di_algor_order_unsat_constr        ! orderuns
  private :: di_algor_back_trans                ! getcard
  private :: di_algor_symmetrize_cartesian      ! symmveclattop
  private :: di_algor_back_trans_converged      ! cnvback_p
  private :: di_algor_get_cart_from_DI          ! newcartd
  private :: di_algor_cart_hessian_from_DI      ! hintbkp
  private :: di_algor_make_U_matrix             ! makeu
  private :: di_algor_check_vector_length       ! chkdp
  private :: di_algor_update_hessian            ! updhesp
  private :: di_algor_make_F_fixed_cart         ! cartfrz
  private :: di_algor_expand_hessian_1          ! exphes1
  private :: di_algor_expand_hessian_2          ! exphes2
  private :: di_algor_expand_hessian_3          ! exphes3
  private :: di_algor_expand_hessian_4          ! exphes4
  private :: di_algor_fixed_atoms_symmetry      ! symfix
  private :: di_algor_constr_from_rigid_body    ! fixatcon
  private :: di_algor_is_atom_in_list           ! klistcar
  private :: di_algor_eigenvector_following     ! optefp
  private :: di_algor_ef_check_H_scale_unsat    ! chkskald
  private :: di_algor_ef_check_H_scale          ! chkskalp
  private :: di_algor_ef_RFOstep                ! formdp
  private :: di_algor_ef_RFOstep_constr         ! conformdp
  private :: di_algor_update_constraints        ! upcons
  private :: di_algor_constr_to_primitives      ! addcons
  private :: di_algor_num_negative_hess_EVs     ! findnegp
  private :: di_algor_cleanup_hessian           ! chkhesp
  private :: di_algor_adjust_hess_EVs_range     ! chkmagp
  private :: di_algor_find_mode_to_follow       ! modeoverlap
  private :: di_algor_add_periodic_box          ! p_box
  private :: di_algor_print_matrix              ! prntmat
  private :: di_algor_print_constraints         ! wrconsti
  private :: di_algor_print_cartesians          ! prcart
  private :: di_algor_string_length             ! lstr
  private :: di_algor_s2                        ! s2         
  private :: di_algor_arc1                      ! arc1
  private :: di_algor_vector_product            ! cross
  private :: di_algor_norm_vector_product       ! normal
  private :: di_algor_norm_vector_difference    ! vektor
  private :: di_algor_distance                  ! dista
  private :: di_algor_unfold_atom_coords        ! cordn
  private :: di_algor_fake_gradients            ! rdgrad0
  private :: di_algor_symmetrize_hessian        ! symhes
  !---------------------------------------------------------------------------!
  !                      P r i v a t e   V a r i a b l e s                    !
  !---------------------------------------------------------------------------!
  integer,save :: ihess !  Hessian status flag
                   ! 0 - have "exact" or initial Cartesian Hessian
                   !     use as is for Cartesian; transform if internals
                   ! 1 - have Hessian from previous step   need to update
                   !-1 - set up default diagonal Hessian

  integer,save           :: num_primitive_bonds ! number of unique bonds
  integer,save           :: num_primitive_angles ! number of unique bends
  integer,save           :: num_primitive_torsions ! number of unique torsions
  integer,save           :: nunsat,nunsat0
  real(kind=di_dp),save :: eold_u

  integer,parameter :: max_number_of_connections = 100 
  integer,save :: max_number_of_bonds 
  integer,save :: nvib ! 3N-3
  integer,save :: intcor ! total number of primitives, including Lrt and Constraints
                    ! if disconnected fragments
                    !   intcor = intcor + num_disconn_fragm 
                    ! if fixed coordinates
                    !   intcor = intcor + num_fixed_coords
  integer,save :: ndeg ! number degrees of freedom, typically 3N-3 (=nvib)
                  ! but could be 3N with constraints
  integer,save :: ns   ! size of the active space
                  ! = ndeg - num_fixed_coords (in schmidt)
                  ! = ndeg - zero grad (in symnic)
  integer,save :: num_disconn_fragm ! number of diconnected fragments = ired - 3
  integer,save :: ncons_1  ! 1 or ncons
  !-------------------------- memory pointers
  integer,save :: num_fixed_coords_2 ! 1 or num_fixed_coords
  integer,save :: ncons_2 ! ncons_1 + num_fixed_constr_rigid_body
  integer,save :: num_disconn_fragm_2 ! 1 or num_disconn_fragm 
  integer,save :: maxdiis_2 ! 1 or maxdiis
  !------------------------------ thresholds
  real (kind=di_dp),parameter :: tp_vdw=1.0_di_dp ! scale initial atomic vdw radii 
  real (kind=di_dp),save      :: tp_benlin=0.01_di_dp !  eliminate linear bends
  real (kind=di_dp),parameter :: tp_smlden = 0.0001_di_dp ! small denominator
  real (kind=di_dp),parameter :: tp_prprim = 1.0e-12_di_dp ! projection of primitives 
  real (kind=di_dp),parameter :: tp_symgrd = 1.0e-9_di_dp ! small internal gradients
  real (kind=di_dp),parameter :: tp_weight = 1.0e-9_di_dp ! weights of coeff matrix
  real (kind=di_dp),save      :: tp_eigvz = 1.0e-3_di_dp ! zero eigenvalue
  real (kind=di_dp),parameter :: tp_schmidt = 1.0e-3_di_dp ! zero in Schmidt routine to select vectors
  real (kind=di_dp),parameter :: tp_schmidt_a = 1.0e-6_di_dp ! zero in Schmidt routine to activate constraints
  real (kind=di_dp),parameter :: tp_backz = 5.0e-9_di_dp ! zero for backtransform
  real (kind=di_dp),parameter :: tp_f_bond = 0.5_di_dp ! starting force constant for bond
  real (kind=di_dp),parameter :: tp_f_angle= 0.2_di_dp ! starting force constant for angle
  real (kind=di_dp),parameter :: tp_f_tors = 0.1_di_dp ! starting force constant for torsions
  real (kind=di_dp),parameter :: tp_f_crt = 0.1_di_dp ! starting force constant for Crt
  real (kind=di_dp),parameter :: tp_back_car = 1.0_di_dp ! maximum Cartesian change in back trns
  real (kind=di_dp),parameter :: tp_bad_hes = -2.0_di_dp ! bad hessian element causes redoing delocalized
  real (kind=di_dp),parameter :: tp_con_sat = 5.0e-4_di_dp ! threshold to recognize satisfied constraints
  real (kind=di_dp),parameter :: thrbnd=0.60_di_dp ! distance ratio for bonding
  real (kind=di_dp),parameter :: toang = 180.0_di_dp/di_pi
  !------------------------------ parameters
  real (kind=di_dp),save      :: dmax = 0.3_di_dp    ! maximum displacement
 
  real (kind=di_dp),parameter :: one=1.0_di_dp
  real (kind=di_dp),parameter :: zero=0.0_di_dp

  !------------------------------
  integer,dimension(:,:),allocatable,save                  :: icc  ! atoms involved in constraint
                                                          ! IC1-IC2           distance constraint
                                                          ! IC1-IC2-IC3       angle constraint
                                                          ! IC1-IC2-IC3-IC4   dihedral constraint
  real(kind=di_dp),dimension(:,:),allocatable,save         :: rcon !  constraint value (in atomic units)

  real(kind=di_dp),dimension(:),allocatable,save           :: coords_DI   ! internal coordinates
  real(kind=di_dp),dimension(:),allocatable,save           :: gradient_DI ! internal gradient
  real(kind=di_dp),dimension(:),allocatable,save           :: disp_DI     ! previous displacement vector (step)
  real(kind=di_dp),dimension(:),allocatable,save           :: gradient_DI_old ! previous internal gradient
  real(kind=di_dp),dimension(:,:),allocatable,save         :: hess_DI     ! Hessian matrix in Internal coordinates
                                                                 ! the size is Ns x Ns ; Ns is the active space < 3*NATOMS
  real (kind=di_dp),dimension(:),allocatable,save          :: vmdel  ! mode to follow (in delocalized coordinates)
  real (kind=di_dp),dimension(:),allocatable,save          :: savtor ! saved torsions

contains


!------------------------------------------------------------------------------
  subroutine geom_algor_DELOC(verify,iGeoCycle,status)
    !=========================================================================!
    ! The top level driver for geometry optimization in delocalized           !
    ! coordinates                                                             !
    !                                                                         !
    !  geom_algor_DELOC itself is simply the "main wrapper" which             !
    !  reads pertinent information from the input files as to                 !
    !  the size of the current system and the optimization                    !
    !  options requested and, based on this information,                      !
    !  allocates the necessary memory. geom_algor_DELOC                       !
    !  completes job input and is responsible for all file I/O.               !
    !  It then calls di_algor_optimize which is the main driver.              !
    !                                                                         !
    !  This routine can be used to just verify correctness                    !
    !  of the OPT data. If verify is true, program will read and analyze      !
    !  input but will not run any optimization.                               !
    !  If verify is .false. verification will be repeated and optimize        !
    !  routine will be called.                                                !
    !                                                                         !
    !                         [del_cor in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) verify (in) - only check the input if .true.                         !
    ! 2) ideocycle (inout) - step number                                      !
    ! 3) status (out)                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:  none                                     !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) ncycle                                                               !
    ! 2) ip_error                                                             !
    ! 3) max_number_of_bonds                                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !  coords_cart and gradient_cart arrays are already filled in             !
    !  (and symmetrized)                                                      !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(inout)     :: iGeoCycle
    logical,intent(in)        :: verify
    integer,intent(out)       :: status
    !-------------------------------------------------------------------------!
    real(kind=di_dp), external :: ddot
    real (kind=di_dp), parameter :: small=1.0e-7_di_dp
    integer :: nat3,ndim,itor,nile,ierr
    real (kind=di_dp) :: rmsg
    !-------------------------------------------------------------------------!

    ip_error = 0
    ncycle = iGeoCycle
    max_number_of_bonds = number_atoms * max_number_of_connections
    ncons_1 = ncons
    ncons_2 = max(ncons_1,1)

    if (l_scan_pes .and. verify) go to 990

    ! if verify (or restart) and periodic and delocalized
    ! this routine does
    ! 1) make delocalized internals
    ! 2) check that backtransformation works

    if (verify) ncycle = 0
    if (l_scan_pes .and. igeocycle==0) ncycle=0

    nat3 = 3*number_atoms
    ndim = nat3

    ! on second pass number of constraints may be larger
    if (ncycle>0) ncons_2 = ncons_2 + num_fixed_constr_rigid_body

    maxdiis_2 = max(maxdiis,1)
    num_fixed_coords_2 = max(num_fixed_coords,1)

    !--- printout

    if (di_on_root .and. di_iprint>1) then
       write(di_stdout,1014)
       if (num_fixed_coords/=0) write(di_stdout,1015)
       if (ncons_1/=0) write(di_stdout,1016)
       if (tsflag) then
          write(di_stdout,1020)
       else
          write(di_stdout,1030)
       endif
    endif


    ! get np_p for redres
    if (global_data_exists) then

       call deloc_utils_read_np_p(np_p,ierr)

       if (verify) then
          write(di_stdout,*) ' Message: No need to construct delocalized internals'
          go to 999
       endif !ncycle=0

    endif ! exst

    call di_algor_allocate(ndim)
 
    if (di_on_root .and. use_deloc_int_output) write(di_stdout,1040) ncycle

    if (ncycle<=0) then

       ! fresh start; verify only
       ihess = -1
       itor = 1

    elseif (ncycle==1) then

       ! accurate gradients available; finish delocalized and do Hessian
       ! make first step
       ihess = -1
       itor = 2

    elseif (ncycle>1) then

       ! continuation; only update Hessian 
       ihess = 1
       itor = 2

    endif


    if (period) then

       ! redres  to [0,1] and not [-.5,.5] in order to define constraints
       ! redres  to [0,1]  cell - at that level internal constarints are defined
       ! fresh start; verify only option
       ! redress atoms into unit cell [0,1] and develop new np_p
       if (ncycle==0) then
          call deloc_utils_rationalize_coords(1)

       elseif (ncycle>0) then
          ! just redres using old Ip_p
          call deloc_utils_rationalize_coords(2)
       endif

    endif

    if (di_on_root .and. use_deloc_int_output) then
       write(di_stdout,'(/,30X,a)') 'Input Coordinates ('//di_energy_label//')'
       call di_algor_print_cartesians
    endif

    if (.not.period) then
       ! make a big box around molecule
       call di_algor_add_periodic_box
    endif

    ! ----- Get Gradients (already there!)

    if (ncycle<=0.and.(.not.L_SCAN_PES)) then

       ! dummy gradient
       call di_algor_fake_gradients

    else

       ! ----- Check gradients ; if very large steepest descent
       !                         if very small exit

       rmsg = ddot(nat3,gradient_cart(1,1),1,gradient_cart(1,1),1)
       rmsg = sqrt(rmsg/float(nat3))

       !  if rmsg is very small then skip geometry optimization

       if (rmsg<small) then
          if (di_on_root) then
             write(di_stdout,'(a, e10.1)') ' Message: Gradients are very small. RMS cartesian gradient: ', rmsg
             write(di_stdout,'(a)') ' Message: Geometry is converged'
          endif
          go to 999
       endif
    endif


    ! ----- Get Hessian from HESSIAN file; only if we have accurate gradients
    if (ncycle==1) then
       ! no Cartesian Hessian found; ihess becomes -1; otherwise ihess = 0
       call deloc_utils_read_hessian(nat3,hess_cartesian,ihess)
    endif

    ! SYMMETRY section: symmetrize coords_cart and gradient_cart, define nvib
    ! In CASTEP we expect forces and coordinates to be symmetrized already.
    nvib = nat3 - 3

    ! zero disp_DI
    disp_DI = zero

    ! the first pass - get all the input data that is not collected yet
    if (ncycle<=0) then
       Eold_u = zero

       ! INTERNAL CONSTRAINTS
       if (ncons_1>0) then
          ncons = ncons_1
          call deloc_utils_read_DI_constraints(ncons,icc,rcon)

          ! WARNING
          if (ntrans>1 .and. di_on_root) then
             write(di_stdout,*) ' Warning: System symmetry may not be correct if constraints are present'
             write(di_stdout,*) ' Warning: Optimization may abort'
          endif

          ! print satisfied constraints
          call di_algor_init_constraints

       else

          nunsat0 = 0
          nunsat = 0

       endif

       ! FIX  CONSTRAINT section; read only once
       if (num_fixed_coords>0) then

          ! WARNING
          if (ntrans>1 .and. di_on_root) then
             write(di_stdout,*) ' Warning: System symmetry may not be correct if atoms are fixed'
             write(di_stdout,*) ' Warning: Optimization may abort'
          endif

          call deloc_utils_read_fix_cartesians

          ! symmetry check
          call di_algor_fixed_atoms_symmetry

          ! check if all atoms of constraints are fixed; stop if yes
          call di_algor_stop_fixed_constraints

       endif

       ! FIX ATOMS (with internals all fixed as constraints)
       if (num_rigid_body_atoms/=0) then
          call deloc_utils_read_fix_atoms
       endif

    endif  ! done only at the beginning
 
    nile = ncycle

    ! restart Hessian if change in the unsatisfied constraints
    if (nunsat/=nunsat0) then
       ihess = -1
       nunsat0 = nunsat
       nile = 1
    endif

    if (di_iprint>4 .and. di_on_root) then
       write(di_stdout,*) 'geom_algor_DELOC:'
       write(di_stdout,*) 'ncycle,ihess,itor ',ncycle,ihess,itor
       write(di_stdout,*) 'ncons,nunsat,ntrans,mfxcon',ncons,nunsat,ntrans,num_fixed_constr_rigid_body
    endif

    ! THIS IS IT - TIME TO WORK
    call di_algor_optimize (ndim,itor,nile)

    if (ip_error==1) then
       ! exit continue with cartesian optimization, if allowed
       if (di_on_root) write(di_stdout,*) ' Warning: switching to optimization with Cartesians'

       ! redress atoms back into the original frame (from cell [0,1]
       if (period) call deloc_utils_rationalize_coords(-1)

       ! update OPT
       if (ncons_2/=1) then
          if (.not.period .and. num_fixed_constr_rigid_body==0) then
             ! call upd_opt(natoms,ndim,icc,rcon) 
          else
             call deloc_utils_io_abort('Error: Constraints may not be valid for optimization with Cartesians')
          endif
       endif

    endif  ! ip_error==1

    global_data_exists = .true.

999 continue

    call di_algor_deallocate

    ! all done
    if (lp_del_test) then
       call deloc_utils_io_abort('End of testing Optimizer')
    endif

    if (ip_error/=0) then
       ncycle = -2
       if (igeocycle==0) igeocycle=1

       ! This next line might be OK for DMol; in CASTEP we have to switch back to Cartesians
       ! on a higher level, based on the status returned by geom_algor_DELOC

       ! ioptc = 0

    endif ! ip_error

    ! return the error code
    status = ip_error

990 continue

 1014 format(//,'** GEOMETRY OPTIMIZATION IN DELOCALIZED COORDINATES **')
 1015 format(/,'** with  FIXED COORDINATES **')
 1016 format(/,'** with  INTERNAL CONSTRAINTS **')
 1020 format(/,'   Searching for a Transition State')
 1030 format(/,'   Searching for a Minimum')
 1040 format(/,'   Optimization Cycle: ',I3)

    return

  end subroutine geom_algor_DELOC



!------------------------------------------------------------------------------

  subroutine di_algor_allocate(ndim)
    !=========================================================================!
    ! Allocate a few arrays that are used throughout DI optimization          !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndim (in) -  number of degrees of freedom                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) ncons_2                                                              !
    ! 2) maxdiis_2                                                            !
    ! 3) num_fixed_coords_2                                                   !
    ! 4) ndeg                                                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: allocate memory for -                 !
    ! 1) coords_DI                                                            !
    ! 2) gradient_DI                                                          !
    ! 3) disp_DI                                                              !
    ! 4) gradient_DI_old                                                      !
    ! 5) hess_DI                                                              !
    ! 6) icc                                                                  !
    ! 7) rcon                                                                 !
    ! 8) map_fixed_coords                                                     !
    ! 9) vmdel                                                                !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in) :: ndim
    !-------------------------------------------------------------------------!
    integer :: ierr
    !-------------------------------------------------------------------------!
    allocate(coords_DI(1:ndim),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating coords_DI in di_algor_allocate')
    coords_DI = zero

    allocate(gradient_DI(1:ndim),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating gradient_DI in di_algor_allocate')
    gradient_DI = zero

    allocate(disp_DI(1:ndim),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating disp_DI in di_algor_allocate')
    disp_DI = zero

    allocate(gradient_DI_old(1:ndim),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating gradient_DI_old in di_algor_allocate')
    gradient_DI_old = zero

    allocate(hess_DI(1:ndim,1:ndim),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating hess_DI in di_algor_allocate')
    hess_DI = zero

    allocate(icc(1:ncons_2,1:21),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating icc in di_algor_allocate')
    icc = 0

    allocate(rcon(1:ncons_2,1:2),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating rcon in di_algor_allocate')
    rcon = zero

    allocate(map_fixed_coords(1:num_fixed_coords_2),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating map_fixed_coords in di_algor_allocate')
    map_fixed_coords = 0

    allocate(vmdel(1:ndeg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating vmdel in di_algor_allocate')
    vmdel = zero

    return
  end subroutine di_algor_allocate

!------------------------------------------------------------------------------

  subroutine di_algor_deallocate
    !=========================================================================!
    ! Dellocate arrays that are used throughout DI optimization               !
    !-------------------------------------------------------------------------!
    ! Arguments: none                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: deallocate memory for -               !
    ! 1) coords_DI                                                            !
    ! 2) gradient_DI                                                          !
    ! 3) disp_DI                                                              !
    ! 4) gradient_DI_old                                                      !
    ! 5) hess_DI                                                              !
    ! 6) icc                                                                  !
    ! 7) rcon                                                                 !
    ! 8) map_fixed_coords                                                     !
    ! 9) vmdel                                                                !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none

    integer :: ierr

    if (allocated(coords_DI)) then
       deallocate(coords_DI,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating coords_DI in di_algor_deallocate')
    endif

    if (allocated(gradient_DI)) then
       deallocate(gradient_DI,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating gradient_DI in di_algor_deallocate')
    endif

    if (allocated(disp_DI)) then
       deallocate(disp_DI,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating disp_DI in di_algor_deallocate')
    endif

    if (allocated(gradient_DI_old)) then
       deallocate(gradient_DI_old,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating gradient_DI_old in di_algor_deallocate')
    endif

    if (allocated(hess_DI)) then
       deallocate(hess_DI,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hess_DI in di_algor_deallocate')
    endif

    if (allocated(icc)) then
       deallocate(icc,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icc in di_algor_deallocate')
    endif

    if (allocated(rcon)) then
       deallocate(rcon,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating rcon in di_algor_deallocate')
    endif

    if (allocated(map_fixed_coords)) then
       deallocate(map_fixed_coords,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating map_fixed_coords in di_algor_deallocate')
    endif

    if (allocated(vmdel)) then
       deallocate(vmdel,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating vmdel in di_algor_deallocate')
    endif

    if (allocated(savtor)) then
       deallocate(savtor,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating savtor in di_algor_deallocate')
    endif

    return
  end subroutine di_algor_deallocate

!------------------------------------------------------------------------------

  subroutine di_algor_optimize (ndim,itor,nile)
    !=========================================================================!
    ! The main driver for geometry optimization in delocalized coordinates    !
    !                                                                         !
    !                         [del_cor2 in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndim (in) -  number of degrees of freedom                            !
    ! 2) itor (out) - integer flag for saving/checking values of primitive    !
    !             torsions (this is done due to possible convergence          !
    !             problems in iterative generation of new Cartesians)         !
    ! 3) nile(out) - internal description of iteration                        !
    !               0 - there are no accurate gradients;                      !
    !                   skip Symmetry and Hessian part;                       !
    !                   store intcor,ndeg                                     !
    !               1 - we have accurate gradients;                           !
    !                   do Symmetry and Hessian part                          !
    !              -1 - redo delocalized (after di_algor_back_trans failed)   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:  essentially all of them                  !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: essentially all of them               !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                  :: ndim
    integer,intent(out)                                 :: itor
    integer,intent(inout)                               :: nile
    !-------------------------------------------------------------------------!

    real(kind=di_dp),parameter :: thrsh=1.0e-8_di_dp
    integer :: ndiis,neg,negreq
    integer :: jumps,jumph,mode,ierr,nat3
    integer :: nfixcons_2
    real(kind=di_dp),dimension(:),allocatable   :: vmode ! Actual mode (vector) followed on previous cycle

    real(kind=di_dp),dimension(:,:),allocatable :: coords_cart0
    real(kind=di_dp),dimension(:,:),allocatable     :: ut
    integer,dimension(:,:),allocatable              :: lidelp
    real(kind=di_dp),dimension(:),allocatable       :: xprim
    integer,dimension(:),allocatable                :: ktyp
    integer,dimension(:,:),allocatable              :: ictyp !  constraint type
                                                         !  1 - fixed distance
                                                         !  2 - fixed bond angle
                                                         !  3 - fixed dihedral angle
    real(kind=di_dp),dimension(:,:),allocatable     :: hpad
    real(kind=di_dp),dimension(:,:),allocatable     :: xstr
    real(kind=di_dp),dimension(:,:),allocatable     :: gstr
    logical :: cnvgd,modeok

    !-------------------------------------------------------------------------!
    allocate(vmode(1:ndim),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating vmode in di_algor_optimize')

    eold = zero
    ndiis = 0

    ! make initial delocalized coordinates and save (PCHK file or whatever
    ! else deloc_utils_write_di would do)

    jumps = 0
    jumph = 0

    ! gxf: read vmode if available
    modeok = .false.
    if (ncycle<=1) then

       call deloc_utils_read_opt_mode(vmode,modeok)

       if (modeok) mode = -1
    endif

    ! a special entry point if di_algor_back_trans fails jumps here only once
 101 continue

    ! after failure of di_algor_back_trans start Update fresh and do Hi = Unew Hcar Unew
    ! make sure itor = 1 if nile==-1
    if (nile==-1) then
       itor=1
       ! bring in dihedrals  if missing
       tp_eigvz = 0.01_di_dp
    endif

    ! beginning (or redoing): generate delocalized internals
    if (nile<=0) then
       if (di_on_root) write(di_stdout,'(a)') 'Message: Generating delocalized internals'

       if (num_fixed_constr_rigid_body/=0) then
          ! fresh start but redefine constraints
          call di_algor_get_RB_deloc_coords(ndim,itor,nile)
       else

          ! fresh start; verify and do a preliminary set of delocalized
          call di_algor_get_init_deloc_coords(ndim,icc,rcon,itor,nile)
       endif

       if (ip_error==1) then
          ! exit continue with cartesian optimization, if allowed
          if (di_on_root) write(di_stdout,*) ' Warning: switching to optimization with Cartesians'
          go to 999 
       endif

       ! jump over if L_SCAN_PES
       if (L_SCAN_PES .and. nile==0) then
          nile = 1
          itor = 2
          go to 103
       endif

       ! return only if this is a genuine begining and not the jump from
       ! di_algor_get_init_deloc_coords

       if (di_on_root) then
          write(di_stdout,'(a)')'Message: Generation of delocalized internals is successful'
          write(di_stdout,'(a)')' '
       endif

       if (nile==0) go to 999  ! end of verify

       ihess = -1
    endif  ! nile

 103 continue

    nat3 = 3*number_atoms
    nfixcons_2 = num_fixed_coords_2 + ncons_2

    ! symmetrize HESSIAN
    if (ncycle==1 .and. ihess==0 .and. ntrans>1) then
       if (di_iprint>2 .and. di_on_root) write(di_stdout,'(a)') 'Symmetrizing Hessian matrix'
       call di_algor_symmetrize_hessian(hess_cartesian,thrsh)
    endif


    ! limited read of PCHK if that was first iteration (Ncycle =1)
    if (allocated(ictyp)) then
       deallocate(ictyp,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ictyp in di_algor_optimize')
    endif
    allocate(ictyp(1:nfixcons_2,1:2),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating ictyp in di_algor_optimize')

    if (allocated(hpad)) then
       deallocate(hpad,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hpad in di_algor_optimize')
    endif
    allocate(hpad(1:ns+ncons,1:ncons_2),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating hpad in di_algor_optimize')

    if (allocated(xstr)) then
       deallocate(xstr,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating xstr in di_algor_optimize')
    endif
    allocate(xstr(1:3*number_atoms,1:MaxDiis_2),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating xstr in di_algor_optimize')

    if (allocated(gstr)) then
       deallocate(gstr,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating gstr in di_algor_optimize')
    endif
    allocate(gstr(1:3*number_atoms,1:MaxDiis_2),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating gstr in di_algor_optimize')

    if (allocated(ut)) then
       deallocate(ut,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ut in di_algor_optimize')
    endif
    allocate(ut(1:intcor,1:ndeg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating ut in di_algor_optimize')

    if (allocated(lidelp)) then
       deallocate(lidelp,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating lidelp in di_algor_optimize')
    endif
    allocate(lidelp(1:4,1:intcor),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating lidelp in di_algor_optimize')

    if (allocated(xprim)) then
       deallocate(xprim,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating xprim in di_algor_optimize')
    endif
    allocate(xprim(1:intcor),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating xprim in di_algor_optimize')

    if (allocated(savtor)) then
       deallocate(savtor,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating savtor in di_algor_optimize')
    endif
    allocate(savtor(1:intcor),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating savtor in di_algor_optimize')

    if (allocated(ktyp)) then
       deallocate(ktyp,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ktyp in di_algor_optimize')
    endif
    allocate(ktyp(1:intcor),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating ktyp in di_algor_optimize')

    call deloc_utils_read_di_data (np_p,ut,lidelp,xprim,savtor,ktyp,map_fixed_coords,bmat_disconn_fragm, &
     &     ictyp,icc,rcon,disp_DI,gradient_DI_old,hess_DI,hpad,vmdel,Xstr,Gstr, &
     &     list_fixed_coords,list_rigid_body_atoms,num_disconn_fragm,nile,ierr)

    if (ierr/=0)  then
       if (di_on_root) write(di_stdout,*) ' Warning: problem with reading deloc_utils_read_di_data     '
       ip_error = 1
       go to 999
    endif

    eold = eold_u
     
    ! redefine new delocalized coords_DI, gradient_DI

    call di_algor_get_final_deloc_coords (ut,lidelp,xprim,ktyp,ictyp, &
     &             hess_cartesian,vmode,modeok,cnvgd,itor,nile)

    if (cnvgd) go to 999

    if (ip_error==1) then
       ! exit continue with cartesian optimization, if allowed
       if (di_on_root) write(di_stdout,*) ' Warning: switching to optimization with Cartesians'
       go to 999
    endif

    !  make a new step disp_DI; add this step to make a new coords_DI
    !  save gradient_DI into gradient_DI_old
    !  update hess_DI
    if (nunsat/=0) then
       ! unsatisfied
       call di_algor_get_step_unsat_constr(ns+nunsat,mode,neg,negreq,hpad,ierr)
    else
       ! satisfied
       call di_algor_get_step(ns,mode,ndiis,neg,negreq,ierr)
    endif
         
    if (ierr/=0)  then

       ! try to redo once more delocalized internals;
       ! switch to cartesians on the second failure
       jumph = jumph  + 1
       if (jumph<=1) then
          ! do di_algor_get_init_deloc_coords again
          nile= -1
          if (di_on_root) then
             write(di_stdout,*) ' Warning:  Hessian has bad eigenvalues'
             write(di_stdout,*) ' Message: New Coordinates will be constructed'
          endif
          go to 101
       endif

       if (di_on_root) write(di_stdout,*) ' Warning: problem with di_algor_get_step '
       ip_error = 1

       !exit if ncons
       if (ncons/=0) go to 125

       go to 999

    else

       nile = ncycle

    endif
 
    !  back transformation
    if (.not.allocated(coords_cart0)) then
       allocate(coords_cart0(1:3,1:number_atoms),stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating coords_cart0 in di_algor_optimize')
    endif
    coords_cart0 = coords_cart

    call di_algor_back_trans(ndim,ut(1:intcor,1:ndeg),lidelp,xprim,ktyp,coords_cart0,hess_cartesian,ierr)

    if (ierr/=0)  then
       if (di_on_root) write(di_stdout,*) ' Warning: problem with di_algor_back_trans '
       ! regenerate initial coords_cart
       coords_cart = coords_cart0

       ! try to redo once more delocalized internals; 
       ! switch to cartesians on the second failure
       jumps = jumps  + 1
       if (jumps<=1) then
          ! do di_algor_get_init_deloc_coords again
          nile= -1
          if (di_on_root) then
             write(di_stdout,*) ' Warning: Backtransformation failed'
             write(di_stdout,*) ' Message: New Coordinates will be constructed'
          endif
          go to 101
       endif

       ip_error = 1

       ! exit if ncons
       if (ncons/=0) go to 125

       go to 999
    else

       nile = ncycle

    endif

    ! check convergence of cartesian gradients and displacement
    call deloc_utils_geom_converged(coords_cart0,ns,neg,negreq,gradient_DI,cnvgd)

    ! write constraints, fixed atoms
    if (ncons/=0) then

       call di_algor_print_constraints(lidelp,ktyp,ictyp)
       !if # of unsatisfied changed jump up; redo delocalized

       if (nunsat/=nunsat0) then
          nile = 0
          itor=1

          if (di_on_root) write(di_stdout,'(a)') 'Message: Generating delocalized internals'

          if (num_fixed_constr_rigid_body/=0) then
             call di_algor_get_RB_deloc_coords(ndim,itor,nile)
          else
             call di_algor_get_init_deloc_coords(ndim,icc,rcon,itor,nile)
          endif

          if (di_on_root) write(di_stdout,'(/,a,/)') 'Message: Generation of delocalized internals is successful'

          if (ip_error==1) then
             ! exit continue with cartesian optimization, if allowed
             if (di_on_root) write(di_stdout,*) ' Warning: switching to optimization with Cartesians'
             go to 999 
          endif
          ! redress atoms back into the original frame (from cell [-0.5,0.5]
          if (period) call deloc_utils_rationalize_coords(-1)
 
          ! write INCOOR (CASTEP does not do anything here)
          call deloc_utils_save_structure

          eold_u = ec                  ! save current energy
          go to 999
       endif
    endif  ! ncons=|= 0

    eold = ec                  ! save current energy
    eold_u = ec

    ! write PCHK
    call deloc_utils_write_di_data (np_p,ut(1:intcor,1:ndeg),lidelp(:,1:intcor),xprim(1:intcor), &
     &     savtor(1:intcor),ktyp(1:intcor),map_fixed_coords,bmat_disconn_fragm, &
     &     ictyp,icc,rcon,disp_DI,gradient_DI,hess_DI,list_fixed_coords,list_rigid_body_atoms, &
     &     num_disconn_fragm,nile,ierr,hpad_loc=hpad,vmode_loc=vmdel,xstr_loc=xstr,gstr_loc=gstr)

    ! redress atoms back into the original frame (from cell [-0.5,0.5]
    if (period) call deloc_utils_rationalize_coords(-1)
 
    ! write INCOOR(CASTEP does not do anything here)
    call deloc_utils_save_structure

    ! write HESS independently of PCHK
    call deloc_utils_write_hessian(nat3,hess_cartesian)

    go to 999

 125 continue

    ! if UNSATISFIED or Symmetry stop
    if (Ncons/=0) then
       if (di_on_root) write(di_stdout,'(/,a)') ' Message: Program failed likely due to imposed constraints'
       if (Nunsat/=0) then
          call deloc_utils_io_abort('Target and actual value of constraints may differ too much; change target')
       endif

       if (ntrans>1) then
          call deloc_utils_io_abort('Symmetry may not be correct if constraints are present; change symmetry to c1(p1)')
       endif

       call deloc_utils_io_abort ('Error termination in di_algor_optimize')
    endif

 999 continue

    ! clean up
    if (allocated(ictyp)) then
       deallocate(ictyp,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ictyp in di_algor_optimize')
    endif
    if (allocated(hpad)) then
       deallocate(hpad,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hpad in di_algor_optimize')
    endif
    if (allocated(xstr)) then
       deallocate(xstr,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating xstr in di_algor_optimize')
    endif
    if (allocated(gstr)) then
       deallocate(gstr,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating gstr in di_algor_optimize')
    endif
    if (allocated(ut)) then
       deallocate(ut,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ut in di_algor_optimize')
    endif
    if (allocated(lidelp)) then
       deallocate(lidelp,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating lidelp in di_algor_optimize')
    endif
    if (allocated(xprim)) then
       deallocate(xprim,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating xprim in di_algor_optimize')
    endif
    if (allocated(ktyp)) then
       deallocate(ktyp,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ktyp in di_algor_optimize')
    endif
    if (allocated(vmode)) then
       deallocate(vmode,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating vmode in di_algor_optimize')
    endif
    if (allocated(coords_cart0)) then
       deallocate(coords_cart0,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating coords_cart0 in di_algor_optimize')
    endif

    return
  end subroutine di_algor_optimize

!---------------------------------------------------------------

  subroutine di_algor_get_RB_deloc_coords(ndim,itor,nile)
    !=========================================================================!
    ! The same as di_algor_get_init_deloc_coords: make B matrix and           !
    ! generate deloaclized coordinates. Also modify the set of constraints -  !
    ! this time take out rigid body constraints from ncons (RB=rigid body)    !
    !                                                                         !
    !                         [get_del3 in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndim (in) -  number of degrees of freedom                            !
    ! 2) itor (out) - integer flag for saving/checking values of primitive    !
    !             torsions (this is done due to possible convergence          !
    !             problems in iterative generation of new Cartesians)         !
    ! 3) nile(in) - internal description of iteration                         !
    !               0 - there are no accurate gradients;                      !
    !                   skip Symmetry and Hessian part;                       !
    !                   store intcor,ndeg                                     !
    !               1 - we have accurate gradients;                           !
    !                   do Symmetry and Hessian part                          !
    !              -1 - redo delocalized (after di_algor_back_trans failed)   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) ncons_2                                                              !
    ! 2) ncons                                                                !
    ! 3) num_fixed_constr_rigid_body                                          !
    ! 4) icc                                                                  !
    ! 5) rcon                                                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) ncons                                                                !
    ! 2) ncons_2                                                              !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                    :: ndim
    integer,intent(out)                                   :: itor
    integer,intent(in)                                    :: nile
    !-------------------------------------------------------------------------!
    integer :: i,j,ierr
    integer,dimension(:,:),allocatable          :: icc_c  
    real(kind=di_dp),dimension(:,:),allocatable :: rcon_c
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(icc,1)<ncons_2 .or. size(icc,2)/=21) &
    &  call deloc_utils_io_abort('Error in di_algor_get_RB_deloc_coords: ICC array is too small')
    if (size(rcon,1)<ncons_2 .or. size(rcon,2)/=2) &
    &  call deloc_utils_io_abort('Error in di_algor_get_RB_deloc_coords: RCON array is too small')
#endif

    ! redefine ncons_2
    ncons_2 = ncons_2 - num_fixed_constr_rigid_body
    ncons   = ncons   - num_fixed_constr_rigid_body

    allocate(icc_c(1:ncons_2,1:21),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating icc_c in di_algor_get_RB_deloc_coords')
    icc_c = 0
    allocate(rcon_c(1:ncons_2,1:2),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating icc_c in di_algor_get_RB_deloc_coords')
    rcon_c = zero

    itor = 1

    ! copy
    do i=1,ncons
       do j=1,21
          icc_c(i,j) = icc(i,j)
       enddo
       rcon_c(i,1) = rcon(i,1)
       rcon_c(i,2) = rcon(i,2)
    enddo
    ! redefine ncons_2

    call di_algor_get_init_deloc_coords(ndim,icc_c,rcon_c,itor,nile)

    deallocate(icc_c,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icc_c in di_algor_get_RB_deloc_coords')
    deallocate(rcon_c,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icc_c in di_algor_get_RB_deloc_coords')

    return
  end subroutine di_algor_get_RB_deloc_coords
!------------------------------------------------------------------------------

  subroutine di_algor_get_init_deloc_coords(ndim,icc,rcon,itor,nile)
    !=========================================================================!
    ! Constructs B-matrix for primitive internals                             !
    ! Makes the final set of delocalized internals.                           !
    !                                                                         !
    !                         [get_del1 in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndim (in) -  number of degrees of freedom                            !
    ! 2) icc (inout)  -  which primitive to constrain                         !
    ! 3) rcon (inout) -  values of the constraints                            !
    ! 4) itor (in) - integer flag for saving/checking values of primitive     !
    !             torsions (this is done due to possible convergence          !
    !             problems in iterative generation of new Cartesians)         !
    ! 5) nile(in) - internal description of iteration                         !
    !               0 - there are no accurate gradients;                      !
    !                   skip Symmetry and Hessian part;                       !
    !                   store intcor,ndeg                                     !
    !               1 - we have accurate gradients;                           !
    !                   do Symmetry and Hessian part                          !
    !              -1 - redo delocalized (after di_algor_back_trans failed)   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_iprint                                                            !
    ! 2) di_on_root                                                           !
    ! 3) di_stdout                                                            !
    ! 4) ncons_2                                                              !
    ! 5) max_number_of_bonds                                                  !
    ! 6) lp_del                                                               !
    ! 7) lp_mdf                                                               !
    ! 8) lp_mdf_del                                                           !
    ! 9) max_number_of_connections                                            !
    !10) number_atoms                                                         !
    !11) num_primitive_bonds                                                  !
    !12) num_primitive_angles                                                 !
    !13) num_primitive_torsions                                               !
    !14) nvib                                                                 !
    !15) num_fixed_coords                                                     !
    !16) num_rigid_body_atoms                                                 !
    !17) lp_discf                                                             !
    !18) num_disconn_fragm_2                                                  !
    !19) num_disconn_fragm                                                    !
    !20) num_fixed_constr_rigid_body                                          !
    !21) num_fixed_coords_2                                                   !
    !22) list_fixed_coords                                                    !
    !23) list_rigid_body_atoms                                                !
    !24) map_fixed_coords                                                     !
    !25) np_p                                                                 !
    !26) coords_di                                                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) hess_DI                                                              !
    ! 2) gradient_DI                                                          !
    ! 3) ip_error                                                             !
    ! 4) intcor                                                               !
    ! 5) savtor                                                               !
    ! 6) bmat_disconn_fragm                                                   !
    ! 7) ns                                                                   !
    ! 8) ndeg                                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                      :: ndim
    integer,dimension(:,:),intent(inout)                    :: icc
    real(kind=di_dp),dimension(:,:),intent(inout)           :: rcon
    integer,intent(in)                                      :: itor
    integer,intent(in)                                      :: nile
    !-------------------------------------------------------------------------!
    real (kind=di_dp),dimension(:),allocatable      :: bndlen
    real (kind=di_dp),dimension(:),allocatable      :: vscal
    integer,dimension (:,:),allocatable             :: latomb
    integer,dimension (:,:,:),allocatable           :: latomc
    integer,dimension (:),allocatable               :: latomcn
    integer,dimension (:),allocatable               :: ktyp
    integer,dimension (:),allocatable               :: ktyp_new
    integer,dimension (:,:),allocatable             :: jctyp_new
    integer,dimension (:,:),allocatable             :: lidelp
    integer,dimension (:,:),allocatable             :: lidelp_new
    real (kind=di_dp),dimension(:,:),allocatable    :: bmat
    real (kind=di_dp),dimension(:),allocatable      :: xprim
    real (kind=di_dp),dimension(:),allocatable      :: xprim_new
    real (kind=di_dp),dimension(:,:),allocatable    :: h
    real (kind=di_dp),dimension(:,:),allocatable    :: hh
    real (kind=di_dp),dimension(:),allocatable      :: eigv
    integer,dimension(:,:),allocatable              :: indexb
    real(kind=di_dp),dimension(:),allocatable       :: savtor_new 
    real(kind=di_dp),dimension(:),allocatable       :: rcon_new 
    real(kind=di_dp),dimension(:,:),allocatable     :: ut_new 
    integer,dimension(:),allocatable                :: iz
    integer,dimension(:,:),allocatable              :: icc_c
    integer,dimension(:),allocatable                :: iab
    integer,dimension(:,:),allocatable              :: icbond
    real(kind=di_dp),dimension(:,:),allocatable     :: rcon_c 

    logical ldone

    real(kind=di_dp) :: scalea,scalep

    integer :: nat3,ierr,i,j,nbonds,natsnum,it
    integer :: listb,lad,nic0,intcor0,ltotal,nic
    integer :: ired,intcor1,intcor2,nic1,ncon
    integer :: nfixcons_2,nfixcons
    integer :: nvb,nfi1
    integer :: ncatom,ncons0
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(icc,1)<ncons_2 .or. size(icc,2)/=21) &
    &  call deloc_utils_io_abort('Error in di_algor_get_init_deloc_coords: ICC array is too small')
    if (size(rcon,1)<ncons_2 .or. size(rcon,2)/=2) &
    &  call deloc_utils_io_abort('Error in di_algor_get_init_deloc_coords: RCON array is too small')
#endif

    nat3 = 3*number_atoms
    ncon = ncons

    if (allocated(gradient_DI)) then
       deallocate(gradient_DI,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating gradient_DI in di_algor_get_init_deloc_coords')
    endif
    allocate(gradient_DI(1:ndim),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating gradient_DI in di_algor_get_init_deloc_coords')
    gradient_DI = zero

    if (allocated(hess_DI)) then
       deallocate(hess_DI,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hess_DI in di_algor_get_init_deloc_coords')
    endif
    allocate(hess_DI(1:ndim,1:ndim),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating hess_DI in di_algor_get_init_deloc_coords')
    hess_DI = zero


    if (di_iprint>4 .and. di_on_root) write(di_stdout,1000)

    ! only for genuine delocalized internals

    !--- FIND UNIQUE BONDS

    !--- MDF READ

    nbonds = 0 ! bonds found in MDF file
    listb = 0
    ldone = .false.
      
    !-------- if mdf read only
    allocate(icbond(1:max_number_of_bonds,1:2),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating icbond in di_algor_get_init_deloc_coords')
    icbond = 0

    if (.not.lp_del) then

       ! Bonds are derived from atomic distances AND from some other source of input
       ! (DMol can use MDF file, for example)
       call deloc_utils_get_alt_bond_list(icbond,nbonds,ierr)
 
       if (ierr==0)  then
          if (lp_mdf) ldone = .true.
          if (lp_mdf_del) ldone = .false.
       else
          if (.not.lp_mdf_del) go to 990
       endif  ! error

    endif ! MDF read 

    if (.not.ldone .and. (lp_mdf_del.or.lp_del)) then
       ! do bonds; add them to the list of alternative bonds, if any

       ! knowing atomic symbols get atomic numbers
       !       call getatno(natoms,atsymb,ian)

       ! get bonds VWaals
       scalea = 1.5_di_dp
       scalep = tp_vdw

       ! scalep is a global scaling
       ! natsnum is the number of atoms, with individual scaling  lvscal
       allocate(vscal(1:number_atoms),stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating vscal in di_algor_get_init_deloc_coords')
       allocate(bndlen(1:number_atoms),stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating bndlen in di_algor_get_init_deloc_coords')
       allocate(latomb(1:max_number_of_bonds,1:5),stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating latomb in di_algor_get_init_deloc_coords')
       allocate(latomc(1:max_number_of_connections,1:4,1:number_atoms),stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating latomc in di_algor_get_init_deloc_coords')
       allocate(latomcn(1:number_atoms),stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating latomcn in di_algor_get_init_deloc_coords')

       vscal = one
       natsnum = 0 
       it=0
       latomb = 0
       latomc = 0
       latomcn = 0

 120   continue

       call di_algor_make_bond_distances(bndlen,scalep,natsnum,vscal)

       if (di_iprint>1 .and. di_on_root) write(di_stdout,'(/,a)') ' Automatic generation of unique bonds:'

       ! find unique bonds and connectivity
       call di_algor_make_bond_list(bndlen,listb,latomb,latomc,latomcn,natsnum,vscal,scalea)

       ! Check natsnum looking for atoms with no connections.
       ! increase scaling of van der waals radii for those atoms and do di_algor_make_bond_distances
       ! repeat this rescue three times before aborting.
      
       it=it+1
       if (natsnum/=0.and. it<3) then
          if (di_on_root .and. di_iprint>1) write(di_stdout,*) 'scale atom vdw  radii ',scalea,natsnum
          scalea = scalea*1.3_di_dp
          go to 120
       endif

    endif  ! lp_del

    ! Add alternative bonds to the list (check that there are no duplicates)
    if (nbonds/=0)  then

       call di_algor_combine_bond_lists(listb,nbonds,lad,latomb,icbond)

       ! print final bonds; calculate final numbers of the primitives
       if (lad/=0) then
          if (di_on_root .and. di_iprint>=1) write(di_stdout,'(a,/)') ' Final bonds including alternative data:'

          call di_algor_print_bonds(listb,latomb,latomc,latomcn)

       endif ! lad

    endif


    ! calculate ltotal (maximum number of primitives)
       
    ! lidelp is not used in this call
    nic0 = 0
    call di_algor_make_list_primitives(nic0,intcor0,listb,latomb,latomc,latomcn,icc)

    ltotal = num_primitive_bonds + num_primitive_angles + num_primitive_torsions
         
    !==   end===================unique bonds only


    ! increase ltotal with number of constraints
    ltotal = ltotal + ncons
    nic0 = ltotal

    allocate(lidelp(1:4,1:ltotal),stat=ierr)  ! z2(1)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating lidelp in di_algor_get_init_deloc_coords')

    !  get intcor0 < nic ; since we eliminate linear bends, dihedrals

    !--- FINAL LIST of PRIMITIVES

    !  check if ncons are present in the list lidelp;
    !  if not add them to the list and correct intcor0
    !  find which primitives correspond to constraints
    call di_algor_make_list_primitives(nic0,intcor0, &
     &          listb,latomb,latomc,latomcn,icc,lidelp=lidelp)

    !    B-MATRIX CONSTRUCTION FOR PRIMITIVE INTERNALS

    !  form the B-Matrix

    nic = intcor0
    intcor = nic

    allocate(bmat(1:12,1:nic),stat=ierr) ! zb(ib)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating bmat in di_algor_get_init_deloc_coords')
    bmat = zero
    allocate(indexb(1:12,1:nic),stat=ierr) ! zb(ibx)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating indexb in di_algor_get_init_deloc_coords')
    indexb = 0
    allocate(xprim(1:nic),stat=ierr) ! z3(ixt)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating xprim in di_algor_get_init_deloc_coords')
    xprim = zero
    if (allocated(savtor)) then
       deallocate(savtor,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating savtor in di_algor_get_init_deloc_coords')
    endif
    allocate(savtor(1:nic),stat=ierr) ! z3(isav)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating savtor in di_algor_get_init_deloc_coords')
    savtor = zero
    allocate(ktyp(1:nic),stat=ierr) ! z3(iktyp)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating ktyp in di_algor_get_init_deloc_coords')
    ktyp = 0

    call di_algor_make_B_matrix(itor,ktyp,bmat,indexb,xprim,lidelp,ierr)

    if (ierr>0) then
       if (di_on_root) write(di_stdout,*) ' Warning: Error in di_algor_make_B_matrix '
       go to 990
    endif

    !    MAKE THE DELOCALIZED INTERNAL COORDINATES
    !    form B(t)*B and diagonalize

    allocate(h(1:nat3,1:nat3),stat=ierr)  ! z3(ih)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating h in di_algor_get_init_deloc_coords')
    h = zero
    allocate(hh(1:nat3,1:nat3),stat=ierr) ! z3(ihh)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating hh in di_algor_get_init_deloc_coords')
    hh = zero
    allocate(eigv(1:nat3),stat=ierr)  ! z3(ie)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating eigv in di_algor_get_init_deloc_coords')
    eigv = zero

    call di_algor_make_F_primitives(nat3,nic,ired,h,hh,eigv,bmat,indexb,lidelp,xprim,ktyp,icc,ierr)

    if (ierr>0) then
       if (di_on_root) write(di_stdout,*) ' Warning: Error in di_algor_make_F_primitives '
       go to 990
    endif

    ! here we have already a final intcor which is le.intcor0
    intcor1 = intcor
    intcor2 = intcor
        
    !    DISCONNECTED FRAGMENTS or CONSTRAINTS
    ! if there are constraints; in that case intcor = intcor+constraints
    ! nvib should be 3n-3
    nic1 = nat3 - ired

    if (di_iprint>4 .and. di_on_root) &
     & write(di_stdout,*) 'nic1,ired, nvib,intcor,nic,ncon ,nfix,mfixat', &
     &    nic1,ired,nvib,intcor,nic,ncon,num_fixed_coords,num_rigid_body_atoms

    num_disconn_fragm = 0
    if (nic1/=nvib .and. lp_discf) then

       ! intcor changes maximum by adding  ired -3
       ! solid
       num_disconn_fragm = ired - 3
       num_disconn_fragm_2 = max(num_disconn_fragm,1)
       if (allocated(bmat_disconn_fragm)) then
          if (size(bmat_disconn_fragm,2)/=num_disconn_fragm_2) then
             deallocate(bmat_disconn_fragm,stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bmat_disconn_fragm in di_algor_get_init_deloc_coords')
             allocate(bmat_disconn_fragm(1:nat3,1:num_disconn_fragm_2),stat=ierr)
             if (ierr/=0) call deloc_utils_io_abort('Error in allocating bmat_disconn_fragm in di_algor_get_init_deloc_coords')
             bmat_disconn_fragm = 0.0_di_dp
          endif
       else
          allocate(bmat_disconn_fragm(1:nat3,1:num_disconn_fragm_2),stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in allocating bmat_disconn_fragm in di_algor_get_init_deloc_coords')
          bmat_disconn_fragm = 0.0_di_dp
       endif

       intcor = intcor + num_disconn_fragm 

    elseif (nic1/=nvib .and. (.not. lp_discf)) then

       ! could not recover
       if (di_on_root) write(di_stdout,1800) nvib,nic1
 1800 format(/,2x,'** WARNING ** Wrong Number of Internal Coordinates', &
     &       ' in <di_algor_make_F_primitives>', &
     &       /,14x,'Should be: ',i3,'  Found: ',i3)

       go to 990
    endif
    
    num_disconn_fragm_2 = max(num_disconn_fragm,1)
         
    ndeg = nvib

    ! cartesian constraints
    if (num_fixed_coords/=0) then
       intcor2 = intcor
       intcor = intcor + num_fixed_coords
       ndeg = nvib + 3
    endif

    ! Fixed atoms
    num_fixed_constr_rigid_body = 0
    if (num_rigid_body_atoms/=0) then

       ! identify Atom_Fixed constraints in the list  lidelp
       ! there are num_fixed_constr_rigid_body constraints amongst intcor coordinates
       ! num_fixed_constr_rigid_body may be close to intcor

       call di_algor_constr_from_rigid_body(intcor1,lidelp)

    endif

    ! sum FIXED and CONSTRAINTS
    nfixcons_2 = num_fixed_coords_2+ncons_2 + num_fixed_constr_rigid_body
    nfixcons = num_fixed_coords+ncons+num_fixed_constr_rigid_body

    ! make ut matrix , save lidel into klist,
    ! savtor,ktyp, ixt into intcor dimensions of z4
    ! jctyp, ircon serves both FIX and CONSTRAINT

    allocate(xprim_new(1:intcor),stat=ierr) ! z4(jxt)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating xprim_new in di_algor_get_init_deloc_coords')
    xprim_new = zero

    allocate(savtor_new(1:intcor),stat=ierr) ! z4(jsav)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating savtor_new in di_algor_get_init_deloc_coords')
    savtor_new = zero

    allocate(ktyp_new(1:intcor),stat=ierr) ! z4(jktyp)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating ktyp_new in di_algor_get_init_deloc_coords')
    ktyp_new = 0

    allocate(jctyp_new(1:nfixcons_2,1:2),stat=ierr) ! z4(jctyp)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating jctyp_new in di_algor_get_init_deloc_coords')
    jctyp_new = 0

    allocate(lidelp_new(1:4,1:intcor),stat=ierr)  ! z4(jklist)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating lidelp_new in di_algor_get_init_deloc_coords')
    lidelp_new = 0

    allocate(rcon_new(1:nfixcons_2),stat=ierr) ! z4(ircon)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating rcon_new in di_algor_get_init_deloc_coords')
    rcon_new = zero

    allocate(ut_new(1:intcor,1:nat3 + nfixcons_2),stat=ierr) ! z4(jut)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating ut_new in di_algor_get_init_deloc_coords')
    ut_new = zero

    do j=1,intcor1
       savtor_new(j) = savtor(j)
       xprim_new(j)  = xprim(j)
       ktyp_new(j)   = ktyp(j)
       do i=1,4
          lidelp_new(i,j) = lidelp(i,j)
       enddo
    enddo

    !----- DISCONNECTED FRAGMENTS

    if (nic1/=nvib .and. lp_discf) then
       call di_algor_make_F_disc_fragm(nat3,intcor,intcor1,ired,h,hh,eigv, &
     &              lidelp_new,xprim_new,ktyp_new,ierr)
       if (ierr/=0) then
          if (di_on_root) write(di_stdout,*) ' Warning: Error in di_algor_make_F_disc_fragm '
          go to 990
       endif

    endif ! disconnected fragments

    !=============CONSTRAINED INTERNALS+++

    if (ncons>0)  then
       ! just update  z4(jctyp), z4(ircon)
       call di_algor_update_constraints(jctyp_new,rcon_new)
    endif

    ! FIXED ATOMS
    if (num_rigid_body_atoms/=0) then
       ! store num_fixed_constr_rigid_body Atom_Fixed constraints from the list  z2(lidel)
       ! in icc  array
       ! update jctyp, ircon  and ncons
       call di_algor_constr_from_rigid_body(intcor1,lidelp_new,jctyp_new,rcon_new,xprim_new)
    endif

    !---- CARTESIAN CONSTRAINTS


    if (num_fixed_coords/=0) then
       call di_algor_make_F_fixed_cart(nat3,intcor,intcor2,ired,h,hh,eigv, &
     &            jctyp_new,rcon_new,xprim_new,ktyp_new,savtor_new,ierr)
       if (ierr/=0) then
          if (di_on_root) write(di_stdout,*) ' Warning: Error in di_algor_make_F_fixed_cart '
          go to 990
       endif
    endif  ! fixed coordinates

    ! redefine ndeg
    nvb = nat3 - ired
    if (nvb/=ndeg) then
        ndeg = nvb
        if (di_on_root) write(di_stdout,*) ' Message: Number degrees of freedom was redefined ',ndeg
    endif

    !---- MAKE UT

    call di_algor_make_U_matrix(ired,ut_new(1:intcor,1:ndeg),h,eigv,bmat,indexb)

    ! schmidt if num_fixed_coords
    if (nfixcons>0) then
       allocate(iz(1:nfixcons_2),stat=ierr)  ! z5(isc4)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating iz in di_algor_get_init_deloc_coords')
       iz = 0
       allocate(icc_c(1:nfixcons_2,1:21),stat=ierr) ! z5(isc1)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating icc_c in di_algor_get_init_deloc_coords')
       icc_c = 0
       allocate(rcon_c(1:nfixcons_2,1:2),stat=ierr) ! z5(isc2)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating rcon_c in di_algor_get_init_deloc_coords')
       rcon_c = zero
       allocate(iab(1:nfixcons_2),stat=ierr) ! z5(isc5)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating iab in di_algor_get_init_deloc_coords')
       iab = 0

       ! copy icc21
       call di_algor_add_constrain_to_fixed(icc(1:,21),jctyp_new(1:,1),iz)
       
    endif

    ns = ndeg
    !=============FIXED CARTESIANS or CONSTRAINED INTERNALS+++
     
    nfi1 = nfixcons
    if (nfixcons>0)  then

       ncatom = ncons-num_fixed_constr_rigid_body
       ncons0 = ncons-num_fixed_constr_rigid_body
       call di_algor_schmidt_constraints(ndeg,intcor,ncatom,nfixcons, & 
     &               jctyp_new,rcon_new,iz,ns,ut_new,iab,ierr)

       ! ncatom a new value for num_fixed_constr_rigid_body = num_fixed_constr_rigid_body - ncatom
       num_fixed_constr_rigid_body = num_fixed_constr_rigid_body - ncatom

       nfi1 = nfi1 - nfixcons 
       if (ierr/=0) go to 990
    endif   ! nfixcons

    if (nfi1/=0) then
       call di_algor_tidy_constraints(ncons-nfi1,ncons0,jctyp_new,iab)
    endif

    ! ns is the optimization space:   ns + num_fixed_coords + ncons = ndeg  
    if (di_iprint>4 .and. di_on_root) then
       write(di_stdout,'(a)')'Values of Primitive Internals'
       do i=1,intcor
          write(di_stdout,'(i5,2f20.14)')i,xprim_new(i)
       enddo
    endif

    ierr = 0
    if (num_fixed_constr_rigid_body/=0) then
       ! you have new constraints ; update icc and rcon and ncons_2
       ! however copy first into z5 since you do not have yet memory for icc
       ! 
       ! update icc , rcon and  call deloc_utils_write_di_data
       call di_algor_write_constr_di_data(ncons_2+num_fixed_constr_rigid_body,ncons0,ut_new(1:intcor,1:ndeg),lidelp_new, &
                 xprim_new,savtor_new,ktyp_new,jctyp_new,nile,ierr)
    else

       ! write data on PCHK (or save in memory)
       call deloc_utils_write_di_data (np_p,ut_new(1:intcor,1:ndeg),lidelp_new(1:4,1:intcor), &
     &           xprim_new(1:intcor),savtor_new(1:intcor),ktyp_new(1:intcor), &
     &           map_fixed_coords,bmat_disconn_fragm,jctyp_new,icc,rcon, &
     &           coords_DI,gradient_DI,hess_DI,list_fixed_coords,list_rigid_body_atoms,num_disconn_fragm,nile,ierr)
    endif

    if (ierr/=0) then
       if (di_on_root) write(di_stdout,'(a)') ' Warning: Error in saving intermediate delocalized data'
       go to 990
    endif

    if (di_iprint>1 .and. di_on_root) then
       write(di_stdout,1400)
       write(di_stdout,1500) ns,ndeg, intcor
    endif

    nic = intcor        ! set NIC to actual number of primitives used

 999 continue

    ! remove memory for primitive finding 
    if (allocated(lidelp)) then
       deallocate(lidelp,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating lidelp in di_algor_get_init_deloc_coords')
    endif
    if (allocated(icbond)) then
       deallocate(icbond,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icbond in di_algor_get_init_deloc_coords')
    endif
    if (allocated(vscal)) then
       deallocate(vscal,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating vscal in di_algor_get_init_deloc_coords')
    endif
    if (allocated(bndlen)) then
       deallocate(bndlen,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bndlen in di_algor_get_init_deloc_coords')
    endif
    if (allocated(latomb)) then
       deallocate(latomb,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating latomb in di_algor_get_init_deloc_coords')
    endif
       if (allocated(latomc)) then
       deallocate(latomc,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating latomc in di_algor_get_init_deloc_coords')
    endif
    if (allocated(latomcn)) then
       deallocate(latomcn,stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating latomcn in di_algor_get_init_deloc_coords')
    endif
    if (allocated(bmat)) then
       deallocate(bmat,stat=ierr) ! zb(ib)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bmat in di_algor_get_init_deloc_coords')
    endif
    if (allocated(indexb)) then
       deallocate(indexb,stat=ierr) ! zb(ibx)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating indexb in di_algor_get_init_deloc_coords')
    endif
    if (allocated(xprim)) then
       deallocate(xprim,stat=ierr) ! z3(ixt)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating xprim in di_algor_get_init_deloc_coords')
    endif
    if (allocated(ktyp)) then
       deallocate(ktyp,stat=ierr) ! z3(iktyp)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ktyp in di_algor_get_init_deloc_coords')
    endif
    if (allocated(h)) then
       deallocate(h,stat=ierr)  ! z3(ih)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating h in di_algor_get_init_deloc_coords')
    endif
    if (allocated(hh)) then
       deallocate(hh,stat=ierr) ! z3(ihh)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hh in di_algor_get_init_deloc_coords')
    endif
    if (allocated(eigv)) then
       deallocate(eigv,stat=ierr)  ! z3(ie)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating eigv in di_algor_get_init_deloc_coords')
    endif
    if (allocated(xprim_new)) then
       deallocate(xprim_new,stat=ierr) ! z4(jxt)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating xprim_new in di_algor_get_init_deloc_coords')
    endif
    if (allocated(savtor_new)) then
       deallocate(savtor_new,stat=ierr) ! z4(jsav)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating savtor_new in di_algor_get_init_deloc_coords')
    endif
    if (allocated(ktyp_new)) then
       deallocate(ktyp_new,stat=ierr) ! z4(jktyp)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ktyp_new in di_algor_get_init_deloc_coords')
    endif
    if (allocated(jctyp_new)) then
       deallocate(jctyp_new,stat=ierr) ! z4(jctyp)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating jctyp_new in di_algor_get_init_deloc_coords')
    endif
    if (allocated(lidelp_new)) then
       deallocate(lidelp_new,stat=ierr)  ! z4(jklist)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating lidelp_new in di_algor_get_init_deloc_coords')
    endif
    if (allocated(rcon_new)) then
       deallocate(rcon_new,stat=ierr) ! z4(ircon)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating rcon_new in di_algor_get_init_deloc_coords')
    endif
    if (allocated(ut_new)) then
       deallocate(ut_new,stat=ierr) ! z4(jut)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ut_new in di_algor_get_init_deloc_coords')
    endif
    if (allocated(iz)) then
       deallocate(iz,stat=ierr)  ! z5(isc4)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating iz in di_algor_get_init_deloc_coords')
    endif
    if (allocated(icc_c)) then
       deallocate(icc_c,stat=ierr) ! z5(isc1)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icc_c in di_algor_get_init_deloc_coords')
    endif
    if (allocated(rcon_c)) then
       deallocate(rcon_c,stat=ierr) ! z5(isc2)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating rcon_c in di_algor_get_init_deloc_coords')
    endif
    if (allocated(iab)) then
       deallocate(iab,stat=ierr) ! z5(isc5)
       if (ierr/=0) call deloc_utils_io_abort('Error in deallocating iab in di_algor_get_init_deloc_coords')
    endif

       
    return

! error exit
 990  continue
 
    ip_error = 1
    go to 999

 1000 format(/,' Attempting to Generate Delocalized Internal Coordinates')
 1400 format(/,'*** OPTIMIZATION USES DELOCALIZED INTERNALS ***',/)
 1500 format(' The preliminary size of active space is : ',i7,/, &
     &         ' There are :      ',i7, '  degrees of freedom',/, &
     &         ' There are likely:',i7, '  primitive internals',/)

  end subroutine di_algor_get_init_deloc_coords
!------------------------------------------------------------------------------

  subroutine di_algor_add_constrain_to_fixed (icc21,jctyp,iz)
    !=========================================================================!
    ! Combines internal constraints and fixed rigid body atom constraints     !
    !                                                                         !
    !                         [cpyicc in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) icc21 (in) -  which primitive to constrain                           !
    ! 2) jctyp (in) -  which coordinates to fix                               !
    ! 3) iz (out) - combined list                                             !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) ncons                                                                !
    ! 2) num_fixed_constr_rigid_body                                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,dimension(:),intent(in)         :: icc21
    integer,dimension(:),intent(in)         :: jctyp
    integer,dimension(:),intent(out)        :: iz
    !-------------------------------------------------------------------------!
    integer :: i,ii
    !-------------------------------------------------------------------------!

    ii = ncons - num_fixed_constr_rigid_body
    do i = 1,ii
       ! this keeps a position of constraint
       iz(i) = icc21(i)
    enddo

    do i = ii+1,ncons 
       ! this keeps position of fixed internal
       iz(i) = jctyp(i)
    enddo

    return
  end subroutine di_algor_add_constrain_to_fixed

!------------------------------------------------------------------------------


  subroutine di_algor_get_final_deloc_coords(ut,lidelp,xprim,ktyp,ictyp, &
  &             hess,vmode,modeok,cnvgd,itor,nile)
    !=========================================================================!
    ! Constructs B-matrix for primitive internals                             !
    ! Makes the final set of delocalized internals.                           !
    !                                                                         !
    !                         [get_del2 in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ut (inout) -  transformation matrix                                  !
    !                  i.e. which linear combination of primitive internals   !
    !                  make up each compound (natural) internal coordinate    !
    ! 2) lidelp (inout) -  list of atoms/cell indices of primitive internal   !
    ! 3) xprim (inout)  -  values of primitive coordinates                    !
    ! 4) ktyp (inout) -  array indicating each primitive internal type.       !
    !                currently these are                                      !
    !                1 - stretch                                              !
    !                3 - bend                                                 !
    !                5 - torsion                                              !
    ! 5) ictyp (inout)  -  which primitive to fix                             !
    ! 6) hess (inout) -  Hessian in cartesians - we might have read one in    !
    !                    It is transformed to internals to set hess_DI        !
    !                    HESS is modified if the gradient term is added in    !
    !                    di_algor_transf_Hess_to_intern (lp_idb=TRUE)         !
    ! 7) vmode (in)  -  mode to follow                                        !
    ! 8) modeok (in)  -  is there a mode to follow?                           !
    ! 9) cnvgd (out)  -  can we skip geometry optimization? Returns TRUE      !
    !                    if initial gradients are very small                  !
    !10) itor (in) - integer flag for saving/checking values of primitive     !
    !             torsions (this is done due to possible convergence          !
    !             problems in iterative generation of new Cartesians)         !
    !11) nile(in) - internal description of iteration                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_iprint                                                            !
    ! 2) di_on_root                                                           !
    ! 3) di_stdout                                                            !
    ! 4) intcor                                                               !
    ! 5) number_atoms                                                         !
    ! 6) ndeg                                                                 !
    ! 7) gradient_cart                                                        !
    ! 8) num_fixed_coords                                                     !
    ! 9) ns                                                                   !
    !10) ncons                                                                !
    !11) ncons_2                                                              !
    !12) ihess                                                                !
    !13) num_fixed_coords_2                                                   !
    !14) icc                                                                  !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) coords_DI                                                            !
    ! 2) gradient_DI                                                          !
    ! 3) vmdel                                                                !
    ! 4) ihess                                                                !
    ! 5) icc                                                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none

    real (kind=di_dp),dimension(:,:),intent(inout)             :: ut
    integer,dimension(:,:),intent(inout)                       :: lidelp
    real (kind=di_dp),dimension(1:intcor),intent(inout)        :: xprim
    integer,dimension(:),intent(inout)                         :: ktyp
    integer,dimension(:,:),intent(inout)                       :: ictyp
    real (kind=di_dp),dimension(:,:),intent(inout)             :: hess
    real (kind=di_dp),dimension(:),intent(in)                  :: vmode
    logical,intent(in)                                         :: modeok
    logical,intent(out)                                        :: cnvgd
    integer,intent(in)                                         :: itor
    integer,intent(in)                                         :: nile
    !-------------------------------------------------------------------------!

    real(kind=di_dp), external :: ddot

    real (kind=di_dp), parameter :: small=1.0e-7_di_dp
    integer :: i,nic,nat3,ierr,istat,ns1
    real (kind=di_dp) :: rmsg

    real (kind=di_dp),dimension(:,:),allocatable             :: bmat
    real (kind=di_dp),dimension(:,:),allocatable             :: bmbt
    real (kind=di_dp),dimension(:,:),allocatable             :: bnew 
    real (kind=di_dp),dimension(:,:),allocatable             :: binv
    real (kind=di_dp),dimension(:),allocatable               :: coeff
    integer,dimension(:,:),allocatable                       :: indexb
    !-------------------------------------------------------------------------!
    nat3 = 3*number_atoms
#ifdef debug
    if (size(ut,1)<intcor .or. size(ut,2)<ndeg) &
    &  call deloc_utils_io_abort('Error in di_algor_get_final_deloc_coords: UT array is too small')
    if (size(lidelp,1)/=4 .or. size(lidelp,2)<intcor) &
    &  call deloc_utils_io_abort('Error in di_algor_get_final_deloc_coords: LIDELP array is too small')
    if (size(xprim)<intcor) &
    &  call deloc_utils_io_abort('Error in di_algor_get_final_deloc_coords: XPRIM array is too small')
    if (size(ktyp)<intcor) &
    &  call deloc_utils_io_abort('Error in di_algor_get_final_deloc_coords: KTYP array is too small')
    if (size(icc,1)<ncons_2 .or. size(icc,2)/=21) &
    &  call deloc_utils_io_abort('Error in di_algor_get_final_deloc_coords: ICC array is too small')
    if (size(vmode,1)<number_atoms) &
    &  call deloc_utils_io_abort('Error in di_algor_get_final_deloc_coords: ICC array is too small')
    if (size(hess,1)<nat3 .or. size(hess,2)<nat3) &
    &  call deloc_utils_io_abort('Error in di_algor_get_final_deloc_coords: HESS array is too small')
    if (size(ictyp,1)<ncons_2+num_fixed_coords_2 .or. size(ictyp,2)/=2) &
    &  call deloc_utils_io_abort('Error in di_algor_get_final_deloc_coords: ICTYP array is too small')
#endif

    ! allocate scratch space
    allocate(binv(1:3*number_atoms,1:ndeg),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating binv in di_algor_get_final_deloc_coords')
    binv = zero
    allocate(bnew(1:3*number_atoms,1:ndeg),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating bnew in di_algor_get_final_deloc_coords')
    bnew = zero
    allocate(bmat(1:12,1:intcor),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating bmat in di_algor_get_final_deloc_coords')
    bmat = zero
    allocate(bmbt(1:ndeg,1:ndeg),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating bmbt in di_algor_get_final_deloc_coords')
    bmbt = zero
    allocate(indexb(1:12,1:intcor),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating indexb in di_algor_get_final_deloc_coords')
    indexb = 0

    !    B-MATRIX CONSTRUCTION FOR PRIMITIVE INTERNALS
    !
    !  form the B-Matrix
    cnvgd = .false.
    nic = intcor

    !--- PERIODIC B 

    call di_algor_make_B_matrix(itor,ktyp,bmat,indexb,xprim,lidelp,ierr)

    if (ierr>0) then
       if (di_on_root) write(di_stdout,*) ' Warning: Error in di_algor_make_B_matrix '
       ip_error = 1
       go to 999
    endif

    !--- transform B(d,x) = U(d,x) * Bp(q,x)

    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,'(/,6x,a)')'UT Matrix'
       call di_algor_print_matrix(ndeg,intcor,ndeg, ut)
       write(di_stdout,'(a)')'END of UT Matrix'
    endif

    call di_algor_U_transform_B_matrix(ut(1:intcor,1:ndeg),bnew,bmat,indexb)

    !  now invert the new B-Matrix
    call di_algor_invert_B_matrix(bnew,bmbt,binv,ierr)


    if (ierr/=0) then
       if (di_on_root) write(di_stdout,1200)
       ip_error = 1
       go to 999
    endif

    !    TRANSFORMATION SECTION - CARTESIANS TO INTERNALS
    !
  
    do i=1,ndeg
       !  coordinates
       coords_DI(i)   = ddot(intcor,ut(1,i),1,xprim(1),1)

       !  gradients
       gradient_DI(i) = ddot(nat3,binv(1,i),1,gradient_cart(1,1),1)
    end do

    ! transform vmode to delocalized internals if requested
    if (modeok) then
       call dgemv('t',nat3,ndeg,one,bnew,nat3,vmode,1,zero,vmdel,1)
    endif

    ! printout
    if (di_iprint>4 .and. di_on_root) then
       write(di_stdout,'(a)')'Values of Primitive Internals'
       do i=1,intcor
          write(di_stdout,'(i5,2f20.14)')i,xprim(i)   
       enddo
       write(di_stdout,'(a)')'END Values of Primitive Internals'
       
       write(di_stdout,'(a)')'Values of Delocalized Internals and gradients'
       do i=1,ndeg
          write(di_stdout,'(i5,2f20.14)')i,coords_DI(i),gradient_DI(i)   
       enddo
       write(di_stdout,'(a)')'END Values of Delocalized Internals'
    endif ! (di_iprint>4)

    ! return if Ncycle >1
    ! with the exception if Nile = -1 ; means di_algor_back_trans failed 
    !                                   delocalized were redone
    !                                   a new U is done
    if (Nile>1) go to 999

    ! if there is no constraints(?) redefine ndeg if gradients are small
    if (num_fixed_coords==0) then
                   ! to be redefined in first iteration 
       ns1 = ns              ! size of optimization space

       call di_algor_remove_0_grad_coords(ut,binv)
       if (ns1/=ns) then
          if (di_iprint>1 .and. di_on_root) write(di_stdout,*) ' Active space was refined in di_algor_remove_0_grad_coords',ns   
       endif

       if (di_iprint>4 .and. di_on_root) write(di_stdout,*) 'j: symnic intcor,nvib,ndeg ',intcor,nvib,ndeg
    endif  ! symmetry part

    ! if gradients very small exit
    if (ns>0) then
       rmsg = ddot(ns,gradient_DI(1),1,gradient_DI(1),1)
       rmsg = sqrt(rmsg/float(ns))
    else
       rmsg=zero
    endif

    !  if rmsg is very small then skip geometry optimization

    cnvgd = rmsg<small
    if (cnvgd) then
       if (di_on_root) then
          write(di_stdout,'(a, e10.1)') ' Gradients are very small. RMS cartesian gradient: ', rmsg
          write(di_stdout,'(a)') ' Geometry is converged'
       endif
       go to 999
    endif

    !  transform the original B-Matrix to the new coordinates
    !
    !  initialize primitive weights

    allocate(coeff(1:intcor),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating coeff in di_algor_get_final_deloc_coords')
    coeff = one

    !  find the weight of each primitive in the final
    !  non-redundant optimization space

    if (di_iprint>4 .and. di_on_root) write(di_stdout,*) 'j: getwght intcor,ns ',intcor,ns

    call di_algor_get_weight_primitives(num_fixed_coords+ncons,ut,coeff)

    !  transform the Hessian if requested

    ! no update hess_DI first time , read cartesian and make hess_DI
    !   ihess = 0 restart
    !   ihess = 1 continue
    !   ihess = -1 fresh start from diagonal

    if (ihess==0 .or. ihess==1) then

       if (di_on_root) write(di_stdout,*) ' Hessian in Cartesians is converted to Internal Coordinates'
        
       call di_algor_transf_Hess_to_intern(hess,ktyp,lidelp,coeff,ut,binv)

    elseif (ihess==-1) then

       !  set up default Hessian
       if (di_on_root) write(di_stdout,'(/,a)') ' Startup Hessian in Internal Coordinates'

       call di_algor_initialize_H_internal(ktyp,ut)

       ihess = 0

    endif  ! ihess

    !    TIDY UP PRIMITIVE SPACE
    call di_algor_remove_0_weight_coords(coeff,ktyp,lidelp,xprim,ictyp,icc(1:ncons_2,18:21),ut)
       
    if (di_iprint>1 .and. di_on_root) write(di_stdout,1400)
    if (di_on_root) write(di_stdout,1500) ns,ndeg, intcor

999 continue

! free scratch memory
    if (allocated(coeff)) then
       deallocate(coeff,stat=istat)
       if (istat/=0) call deloc_utils_io_abort('Error in deallocating coeff in di_algor_get_final_deloc_coords')
    endif
    deallocate(binv,stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in deallocating binv in di_algor_get_final_deloc_coords')
    deallocate(bnew,stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in deallocating bnew in di_algor_get_final_deloc_coords')
    deallocate(bmat,stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in deallocating bmat in di_algor_get_final_deloc_coords')
    deallocate(bmbt,stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in deallocating bmbt in di_algor_get_final_deloc_coords')
    deallocate(indexb,stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in deallocating indexb in di_algor_get_final_deloc_coords')

    return

 1200 format(/,2x,' Warning: Unable to Invert B-Matrix')
 1400 format(/,'*** OPTIMIZATION USES DELOCALIZED INTERNALS ***',/)
 1500 format(/,' The size of active space is : ',i5,/, &
     &         ' There are :',i7, '  degrees of freedom',/, &
     &         ' There are :',i7, '  primitive internals',/)

  end subroutine di_algor_get_final_deloc_coords

!------------------------------------------------------------------------------

  subroutine di_algor_remove_0_weight_coords(coeff,ktyp,lidelp,xint,ictyp,icc,ut) 
    !=========================================================================!
    ! tidyup delocalized ccordinates by considering symmetry and weights      !
    ! Remove all internals from the primitive space with                      !
    ! weights less than threshold                                             !
    !                                                                         !
    !                         [tidyupp in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  ndeg    -  size of active space                                        !
    !  intcor  -  on input current size of primitive space                    !
    !             on output new size of space                                 !
    !  threshold   -  weight threshold for eliminating primitive              !
    !  coeff   -  weighting of primitives in active space                     !
    !  ktyp    -  array indicating each primitive internal type.              !
    !             currently these are                                         !
    !                1 - stretch                                              !
    !                3 - bend                                                 !
    !                4 - out-of-plane bend                                    !
    !                5 - torsion                                              !
    !  lidelp  -  list of atoms/cell indices of primitive internal            !
    !  savtor  -  array for storing primitive torsions                        !
    !             (possible sign changes near limiting values)                !
    !  ictyp   -  which primitive to fix                                      !
    !  ncons   -  >0 primitive internal constraints                           !
    !  icc     -  which primitive to constrain                                !
    !  ut      -  transformation matrix                                       !
    !              i.e. which linear combination of primitive internals       !
    !                   make up each compound (natural) internal coordinate   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_iprint                                                            !
    ! 2) di_on_root                                                           !
    ! 3) di_stdout                                                            !
    ! 4) intcor                                                               !
    ! 5) ndeg                                                                 !
    ! 6) ncons                                                                !
    ! 7) ncons_2                                                              !
    ! 8) tp_weight                                                            !
    ! 9) num_fixed_coords                                                     !
    !10) num_fixed_coords_2                                                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) intcor                                                               !
    ! 2) num_primitive_bonds                                                  !
    ! 3) num_primitive_angles                                                 !
    ! 4) num_primitive_torsions                                               !
    ! 5) savtor                                                               !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none

    real (kind=di_dp), dimension(:), intent(inout)        :: coeff
    integer, dimension(:),intent(inout)                   :: ktyp
    integer, dimension(:,:), intent(inout)                :: lidelp
    real (kind=di_dp), dimension(:), intent(inout)        :: xint
    integer, dimension(:,:),intent(inout)                 :: ictyp
    integer, dimension(:,:),intent(inout)                 :: icc
    real (kind=di_dp), dimension(:,:), intent(inout)      :: ut
    !-------------------------------------------------------------------------!

    integer :: i,ii,j,it,ni,kc,i0,j0,k0,ia,jb,ld,ierr
    integer :: itra,itrb,itrc,ila,ilb,ilc,ika,ikb,ikc
    integer :: intcor0,intcor00
    real (kind=di_dp) :: threshold
    integer, dimension(1:19) :: ityp
    real (kind=di_dp), dimension(:,:),allocatable   :: utt
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(ut,1)<intcor .or. size(ut,2)<ndeg) &
    &  call deloc_utils_io_abort('Error in di_algor_remove_0_weight_coords: UT array is too small')
    if (size(coeff,1)<intcor) call deloc_utils_io_abort('Error in di_algor_remove_0_weight_coords: COEFF array is too small')
    if (size(ktyp,1)<intcor) call deloc_utils_io_abort('Error in di_algor_remove_0_weight_coords: KTYP array is too small')
    if (size(xint,1)<intcor) call deloc_utils_io_abort('Error in di_algor_remove_0_weight_coords: XINT array is too small')
    if (size(icc,1)<ncons_2 .or. size(icc,2)/=4) &
    &  call deloc_utils_io_abort('Error in di_algor_remove_0_weight_coords: ICC array is too small')
    if (size(lidelp,1)/=4 .or. size(lidelp,2)<intcor) &
    &  call deloc_utils_io_abort('Error in di_algor_remove_0_weight_coords: LIDELP array is too small')
#endif

    allocate(utt(1:ndeg,1:intcor),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating utt in di_algor_remove_0_weight_coords')
    utt = zero

    threshold = tp_weight

    do i=1,19
       ityp(i)=0
    enddo

    !  first check if we actually have anything to eliminate

    do i=1,intcor
       if (coeff(i)<threshold) go to 20
    end do
      
    !  at this point every primitive should be kept
    !  so why are we still here?

    go to 999

    !  there is something to discard
 20 continue

    !  initialize
    intcor0=intcor

    ii = 0

    ! now check weights for elimination

    do i=1,intcor

       ! do not eliminate supplementary coodinates
       if (coeff(i)>=threshold .or. ktyp(i)==6) then

          !  keep primitive coordinate

          ii = ii+1
          if (i/=ii) then
             ktyp(ii) = ktyp(i)

             lidelp(1,ii) = lidelp(1,i)
             lidelp(2,ii) = lidelp(2,i)
             lidelp(3,ii) = lidelp(3,i)
             lidelp(4,ii) = lidelp(4,i)
             coeff(ii) = coeff(i)
             savtor(ii) = savtor(i)
             xint(ii) = xint(i)
          endif

          !  take care with UT as leading dimension is changing
          !  copy into UTT array as a transpose

          do j=1,ndeg
             utt(j,ii) = ut(i,j)
          end do

          ! fixed section
          if (num_fixed_coords+ncons>0) then
             do j=1,num_fixed_coords+ncons
               if (ictyp(j,1)==i) ictyp(j,2) = ii
             enddo
          endif

          ! constrained section
          if (ncons>0) then
             do j=1,ncons
                if (icc(j,1)==i) icc(j,2) = ii
             enddo
          endif

       else

          !  eliminate it
          !  (keep tabs for print out)

          it = ktyp(i)
          ityp(it) = ityp(it) + 1

       endif

    end do ! i=1,intcor
     
    !  reset number of primitives

    intcor = ii

    !  restore revised internal coordinates

    do i=1,ndeg
       do j=1,intcor
          ut(j,i) = utt(i,j)
       end do
    end do

    if (di_iprint>=2 .and. di_on_root) then
       write(di_stdout,1000)
       if (ityp(1)>0) write(di_stdout,1100) ityp(1)
       if (ityp(2)>0) write(di_stdout,1150) ityp(2)
       if (ityp(3)>0) write(di_stdout,1200) ityp(3)
       if (ityp(4)>0) write(di_stdout,1300) ityp(4)
       if (ityp(5)>0) write(di_stdout,1400) ityp(5)
       if (ityp(6)>0) write(di_stdout,1420) ityp(6)
       if (ityp(7)>0) write(di_stdout,1440) ityp(7)
    endif

    if (di_iprint>=1 .and. di_on_root) write(di_stdout,1500) intcor

    if (intcor/=intcor0) then

       intcor00=intcor
       if (num_fixed_coords>0) intcor00=intcor-num_fixed_coords

       num_primitive_bonds=0
       num_primitive_angles=0
       num_primitive_torsions=0

       ! reset numbers of bonds/angles/torsions
       do ni=1,intcor00
          call deloc_utils_unpack_atom_indices(lidelp(1,ni), ika, ikb, ikc, kc)
          call deloc_utils_unpack_atom_indices(lidelp(2,ni),  i0,  j0,  k0, ia)
          call deloc_utils_unpack_atom_indices(lidelp(3,ni),itra,itrb,itrc, jb)
          call deloc_utils_unpack_atom_indices(lidelp(4,ni), ila, ilb, ilc, ld)

          if (ktyp(ni)==1) then
             num_primitive_bonds=num_primitive_bonds+1
             jb=0
             ld=0
          endif
          
          if (ktyp(ni)==3) then
             num_primitive_angles=num_primitive_angles+1
             ld=0
          endif

          if (ktyp(ni)==5) num_primitive_torsions=num_primitive_torsions+1

          if (di_iprint>4 .and. di_on_root) then
             write(di_stdout,'(i4,2x,i4,2x,5i5)') ni,ktyp(ni),kc,ia,jb,ld
          endif
       enddo

       if (di_iprint>5 .and. di_on_root) then
          write(di_stdout,'(/,6x,a)')'UT Matrix'
          call di_algor_print_matrix(ndeg,intcor,ndeg, ut)
          write(di_stdout,'(a)')'END of UT Matrix'
       endif

       if (num_fixed_coords+ncons>0) then
           do i=1,num_fixed_coords+ncons
              ictyp(i,1) = ictyp(i,2)
              if (di_iprint>4 .and. di_on_root) write(di_stdout,*) 'nfixcons ',i,ictyp(i,1),ictyp(i,2)
           enddo
       endif

       if (ncons>0) then
          do i=1,ncons
             icc(i,1) = icc(i,2)
             if (di_iprint>4 .and. di_on_root) write(di_stdout,*) ' cons ',i,icc(i,1),icc(i,2)
          enddo
       endif

    endif ! change intcor

 999 continue

    deallocate(utt,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating utt in di_algor_remove_0_weight_coords')

    return

 1000 format(/,' Eliminating Redundant Primitive Internals from Space')
 1100 format('  Removed ',i6,' Stretches')
 1150 format('  Removed ',i6,' Logarithmic Stretches')
 1200 format('  Removed ',i6,' Bends')
 1300 format('  Removed ',i6,' Out-of-Plane Bends')
 1400 format('  Removed ',i6,' Torsions')
 1420 format('  Removed ',i6,' Supplementary Coordinates')
 1440 format('  Removed ',i6,' Linear Perpendicular Bends')
 1500 format(' There are now ',i6,' Primitive Internals')

  end subroutine di_algor_remove_0_weight_coords

!------------------------------------------------------------------------------

  subroutine di_algor_remove_0_grad_coords(ut,binv)
    !=========================================================================!
    ! Remove symmetry-redundant composite internals (actually, simply         !
    ! remove coordinates with very small gradients - the old comments seem    !
    ! to assume that would happen onyl because of symmetry).                  !
    !                                                                         !
    !                         [symnicp in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ut (inout)   - transformation matrix                                 !
    !             (vectors of linear coefficients relating                    !
    !              primitive and natural internals)                           !
    ! 2) binv (inout) - inverse b-matrix                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) number_atoms                                                         !
    ! 2) di_iprint                                                            !
    ! 3) di_on_root                                                           !
    ! 4) di_stdout                                                            !
    ! 5) intcor                                                               !
    ! 6) ns                                                                   !
    ! 7) ndeg                                                                 !
    ! 8) ncons                                                                !
    ! 9) nunsat                                                               !
    !10) gradient_DI                                                          !
    !11) ihess (if hessian needs to be transformed we need to eliminate       !
    !           symmetry-redundant columns of binv)                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) ndeg                                                                 !
    ! 2) ns                                                                   !
    ! 3) gradient_DI                                                          !
    ! 4) coords_DI                                                            !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    real (kind=di_dp),dimension(:,:),intent(inout)        :: ut
    real (kind=di_dp),dimension(:,:),intent(inout)        :: binv
    !-------------------------------------------------------------------------!

    real (kind=di_dp) :: threshold,val
    integer :: i,iii,it,icstartu,icstarts,nat3
    character (len=10) :: string
    !-------------------------------------------------------------------------!
    nat3 = 3*number_atoms
    !  eliminate all natural internals with gradient below threshold

    threshold = tp_symgrd

    !  threshold   -  threshold for elimination of natural internals
    !             (if gradient below threshold - eliminate)

    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,1200)
       call di_algor_print_matrix(ns,intcor,ns,ut)
       write(di_stdout,'(a,a)')'END of Final Set of Delocalized Internal Coordinates'
    endif

    !  unsatisfied  constraints start
    icstartu = ndeg - ncons

    !  satisfied  constraints start
    icstarts = ndeg - ncons + nunsat

    it = 0
    do i=1,ndeg
       val = abs(gradient_DI(i))

       if (val<=threshold) then
          if (di_iprint>5 .and. di_on_root) write(di_stdout,'(a,i5)')'Eliminating Internal ',i

          ! if the internal i belongs to constraints
          if (i>icstartu) then

             if (i<=icstarts) then
                write (string,'(i8)') i - icstartu

                !  1. stop the program if this is not
                call deloc_utils_io_abort('Error in di_algor_remove_0_grad_coords: constraint'//string//'will break symmetry')

             elseif (i>icstarts) then

                !  2. drop this constraint if this is satisfied
                if (di_on_root) write(di_stdout,*) ' This constraint is satisfied by symmetry'
                it = it + 1

                ! shift all the data (coords, gradients, hessian, etc.)
                call dcopy(intcor,ut(1,i),1,ut(1,it),1)
                if (ihess==0) call dcopy(nat3,binv(1,i),1,binv(1,it),1)
                coords_DI(it) = coords_DI(i)
                gradient_DI(it) = gradient_DI(i)
             endif
          endif
 
       else

          ! shift all the data (coords, gradients, hessian, etc.)
          it = it + 1
          call dcopy(intcor,ut(1,i),1,ut(1,it),1)
          if (ihess==0) call dcopy(nat3,binv(1,i),1,binv(1,it),1)
          coords_DI(it) = coords_DI(i)
          gradient_DI(it) = gradient_DI(i)

       endif

    enddo

    iii = ndeg-it
    if (di_iprint>=2 .and. iii/=0 .and. di_on_root) write(di_stdout,1000) ndeg-it

    ns   = ns   -iii
    ndeg = ndeg -iii

    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,1200)
       call di_algor_print_matrix(ns,intcor,ns,ut)
       write(di_stdout,'(a,a)')'END of Final Set of Delocalized Internal Coordinates'
    endif

    return

 1000 format(i6,'   Coordinates Fixed Due to Symmetry')
 1200 format(/,'   Final Set of Delocalized Internal Coordinates')

  end subroutine di_algor_remove_0_grad_coords

!------------------------------------------------------------------------------
  subroutine di_algor_get_weight_primitives(ncon,ut,coeff)
    !=========================================================================!
    ! Determine the weight of each primitive in the final                     !
    ! symmetry non-redundant optimization space                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ncon (in)   -  number of constraints (fixed primitives)              !
    ! 2) ut (in)     -  transformation matrix                                 !
    !                   i.e. which linear combination of primitive internals  !
    !                   make up each compound (natural) internal coordinate   !
    ! 3) coeff (in)  -  on exit contains primitive weights                    !
    !                                                                         !
    !                         [getwghtp in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) tp_weight                                                            !
    ! 2) di_iprint                                                            !
    ! 3) di_on_root                                                           !
    ! 4) di_stdout                                                            !
    ! 5) intcor                                                               !
    ! 6) ns                                                                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer, intent(in)                                        :: ncon
    real (kind=di_dp),dimension(:,:), intent(in)               :: ut
    real (kind=di_dp),dimension(:), intent(out)                :: coeff
    !-------------------------------------------------------------------------!
    real (kind=di_dp) :: threshold,val
    integer :: i,j
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(coeff)<intcor) call deloc_utils_io_abort('Error in di_algor_get_weight_primitives: COEFF array is too small')
    if (size(ut,1)<intcor .and. size(ut,2)<ns+ncon) &
    &   call deloc_utils_io_abort('Error in di_algor_get_weight_primitives: UT array is too small')
#endif

    threshold = tp_weight

    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,*) 'intcor, ndeg, ncon ',intcor, ns, ncon
       call di_algor_print_matrix(ns+ncon,intcor,ns+ncon,ut)
       write(di_stdout,'(a,a)')'END of Final Set of Delocalized Internal Coordinates'
    endif

    do i=1,intcor
       val = zero
       do j=1,ns+ncon
          val = val + ut(i,j)*ut(i,j)
       end do
       if (abs(val)<threshold) val = zero
       coeff(i) = val
    end do

    !  make sure any fixed primitives get a weighting or else
    !  the Cartesian back-transformation will fail

    do i=ns+1,ns+ncon
       do j=1,intcor
          if (coeff(j)==zero.and.abs(ut(j,i))>threshold) coeff(j) = one
       end do
    end do

    if (di_iprint>4 .and. di_on_root) then
       write(di_stdout,1000)
       do i=1,intcor
          write(di_stdout,1100) i,coeff(i)
       enddo
       write(di_stdout,'(/,a)')'END of Primitive Weights'
    endif

    return

 1000 format(/,'  Primitive Weights in Final Non-Redundant Optimization Space',/)
 1100 format('  Weight for primitive ',i5,' is  ',f10.6)

  end subroutine di_algor_get_weight_primitives

  !---------------------------------------------------------------------------!
  !                P R I V A T E   R O U T I N E S                            !
  !---------------------------------------------------------------------------!
  subroutine di_algor_make_bond_distances (bndlena,scale,natsnum,vscal)
    !=========================================================================!
    ! Sets bond distances for every atom based on vdw radii and scaling       !
    ! MATERIALS STUDIO COVALENT RADII , Dec.2000                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) bndlena(out) - bond distances for every atom                         !
    ! 2) scale(in) - scaling factor  of vdw (default is 1.0)   (tp_vdw)       !
    ! 3) natsnum(in) - how many atoms with individual scaling                 !
    ! 4) vscal(in) - individual scaling factors                               !
    !                                                                         !
    !                         [fbond in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) number_atoms                                                         !
    ! 2) di_iprint                                                            !
    ! 3) di_on_root                                                           !
    ! 4) di_stdout                                                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !  bndlen - list of van der Waals radii as used by Materials Studio       !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none

    real (kind=di_dp),dimension(:),intent(out)        :: bndlena
    real (kind=di_dp), intent(in)                     :: scale
    integer, intent(in)                               :: natsnum
    real (kind=di_dp),dimension(:),intent(in)         :: vscal
    !-------------------------------------------------------------------------!

    real (kind=di_dp),dimension(1:103) :: bndlen(103)      ! covalent radii in Angstrom

    !                    H-H  He-He Li-Li Be-Be  B-B   C-C   N-N   O-O
    data bndlen(1:8)   /0.37, 0.32, 1.31, 0.91, 0.71, 0.77, 0.74, 0.74/

    !                    F-F  Ne-Ne Na-Na Mg-Mg Al-Al Si-Si  P-P   S-S
    data bndlen(9:16)  /0.72, 0.69, 1.66, 1.36, 1.11, 1.17, 1.10, 1.04/

    !                   Cl-Cl Ar-Ar  K-K  Ca-Ca Sc-Sc Ti-Ti  V-V  Cr-Cr
    data bndlen(17:24) /0.99, 0.97, 2.06, 1.66, 1.46, 1.26, 1.21, 1.26/

    !                   Mn-Mn Fe-Fe Co-Co Ni-Ni Cu-Cu Zn-Zn Ga-Ga Ge-Ge
    data bndlen(25:32) /1.26, 1.26, 1.21, 1.21, 1.21, 1.21, 1.16, 1.22/

    !                   As-As Se-Se Br-Br Kr-Kr Rb-Rb Sr-Sr  Y-Y  Zr-Zr
    data bndlen(33:40) /1.21, 1.17, 1.14, 1.10, 2.21, 1.86, 1.66, 1.41/

    !                   Nb-Nb Mo-Mo Tc-Tc Ru-Ru Rh-Rh Pd-Pd Ag-Ag Cd-Cd
    data bndlen(41:48) /1.31, 1.31, 1.21, 1.16, 1.21, 1.26, 1.46, 1.41/

    !                   In-In Sn-Sn Sb-Sb Te-Te  I-I  Xe-Xe Cs-Cs Ba-Ba
    data bndlen(49:56) /1.41, 1.40, 1.41, 1.37, 1.33, 1.30, 2.46, 2.01/

    !                   La-La Ce-Ce Pr-Pr Nd-Nd Pm-Pm Sm-Sm Eu-Eu Gd-Gd
    data bndlen(57:64) /1.81, 1.71, 1.71, 1.71, 1.71, 1.71, 1.71, 1.66/

    !                   Tb-Tb Dy-Dy Ho-Ho Er-Er Tm-Tm Yb-Yb Lu-Lu Hf-Hf
    data bndlen(65:72) /1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.61, 1.41/

    !                   Ta-Ta  W-W  Re-Re Os-Os Ir-Ir Pt-Pt Au-Au Hg-Hg
    data bndlen(73:80) /1.31, 1.21, 1.21, 1.16, 1.21, 1.21, 1.21, 1.36/

    !                   Tl-Tl Pb-Pb Bi-Bi Po-Po At-At Rn-Rn Fr-Fr Ra-Ra
    data bndlen(81:88) /1.76, 1.66, 1.46, 1.76, 1.45, 1.45, 2.82, 2.01/

    !                   Ac-Ac Th-Th Pa-Pa  U-U   Np-Np Pu-Pu Am-Am  Cm-Cm
    data bndlen(89:96) /1.81, 1.66, 1.66, 1.61,  1.61, 1.61, 1.61,  1.74/
    
    !                    Bk-Bk Cf-Cf Es-Es Fm-Fm Md-Md No-No Lr-Lr   
    data bndlen(97:103) /1.70, 1.86, 1.86, 1.86, 1.86, 1.86, 1.86/

    real (kind=di_dp),parameter :: bond_threshold=1.15_di_dp
    real (kind=di_dp) :: thrb

    integer :: i,j
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(bndlena)<number_atoms) call deloc_utils_io_abort('Error in di_algor_make_bond_distances: bndlena array is too small')
    if (size(vscal)<number_atoms) call deloc_utils_io_abort('Error in di_algor_make_bond_distances: vscal array is too small')
#endif

    thrb = bond_threshold * scale * bohr_di

    if (di_on_root .and. natsnum/=0 .and. di_iprint>4) then
       do j=1,number_atoms 
          write(di_stdout,'(a,i3,a,f10.4)') 'Size of atom # ',j,'  is scaled by ',vscal(j)
       enddo
    endif

    do i = 1,number_atoms
       bndlena(i) = bndlen(atomic_numbers(i))*thrb   ! vdW radius times bonding threshold
                                                     ! (contains global scaling factor, SCALE)
       bndlena(i) = bndlena(i)*vscal(i)              ! Add atom-specific scaling factor
    enddo

    return

  end subroutine di_algor_make_bond_distances

!------------------------------------------------------------------------------

  subroutine di_algor_make_bond_list(bndlena,listb,latomb,latomc,latomcn,natsnum,vscal,scale)
    !=========================================================================!
    ! Find unique bonds and connectivity                                      !
    !                                                                         !
    !                         [fbond in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) bndlena (in) : list of bonding distances (vdW) per atom              !
    ! 2) listb (out) : number of unique bonds found in the structure          !
    ! 3) latomb (out) : list of atoms involved in each bond (IDs and cells)   !
    ! 4) latomc (out) : list of all connections                               !
    ! 5) latomcn (out) : number of connections on each atom                   !
    ! 6) natsnum (out) : number of atoms without bonds                        !
    ! 7) vscal (inout) : scaling factors for vdW radii - used in bond         !
    !                    calculations. the values for atoms with no bonds     !
    !                    are reset in this routine to SCALE (see below)       !
    ! 8) scale (in) : scaling factor to use for atoms with no bonds; it is    !
    !                 assigned to the appropriate VSCAL element and is used   !
    !                 in subsequent calls to di_algor_make_bond_distances     !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) max_number_of_connections                                            !
    ! 2) di_iprint                                                            !
    ! 3) di_on_root                                                           !
    ! 4) di_stdout                                                            !
    ! 5) number_atoms                                                         !
    ! 6) coords_cart                                                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    real (kind=di_dp),dimension(:),intent(in)             :: bndlena
    integer,intent(out)                                   :: listb
    integer,dimension(:,:),intent(inout)                  :: latomb
    integer,dimension(:,:,:),intent(out)                  :: latomc
    integer,dimension(:),intent(out)                      :: latomcn
    integer,intent(out)                                   :: natsnum
    real (kind=di_dp),dimension(:),intent(inout)          :: vscal
    real (kind=di_dp),intent(in)                          :: scale
    !-------------------------------------------------------------------------!
    integer,dimension (:,:),allocatable          :: icnnct
    integer,dimension (:),allocatable            :: icnnct2
    integer,dimension (:),allocatable            :: icnnm
    integer,dimension (:,:,:),allocatable        :: icntrns
    integer,dimension (:,:),allocatable          :: icntrns2

    real (kind=di_dp),dimension(1:3) :: tt
    logical :: skip
    real (kind=di_dp) :: bondl1,dist1s,dist2s

    integer :: i,ii,j,jj,ia,jb,m1,m2,m3,inmj,itra,itrb,itrc,ita,itb,itc,icn
    integer :: iii,jb2,jb1,list,ierr
    character (len=80) :: string
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(bndlena)<number_atoms) call deloc_utils_io_abort('Error in di_algor_make_bond_list: Array BNDLENA is too small')
    if (size(vscal)<number_atoms) call deloc_utils_io_abort('Error in di_algor_make_bond_list: Array VSCAL is too small')
    if (size(latomb,1)<max_number_of_bonds .and. size(latomb,2)<5) &
       &  call deloc_utils_io_abort('Error in di_algor_make_bond_list : Array LATOMB is too small')
    if (size(latomcn,1)<number_atoms) &
       &  call deloc_utils_io_abort('Error in di_algor_make_bond_list : Array LATOMCN is too small')
    if (size(latomb,1)<max_number_of_connections .and. size(latomb,2)<4 .and. size(latomc,3)<number_atoms) &
       &  call deloc_utils_io_abort('Error in di_algor_make_bond_list : Array LATOMC is too small')
#endif

    ! Initialize workspace
    allocate(icnnct(1:max_number_of_connections,1:number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating icnnct in di_algor_make_bond_list')
    allocate(icnnct2(1:max_number_of_connections),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating icnnct2 in di_algor_make_bond_list')
    allocate(icnnm(1:number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating icnnm in di_algor_get_init_deloc_coords')
    allocate(icntrns(1:3,1:max_number_of_connections,1:number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating icntrns in di_algor_get_init_deloc_coords')
    allocate(icntrns2(1:3,1:max_number_of_connections),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating icntrns2 in di_algor_get_init_deloc_coords')
    icnnct = 0
    icnnct2 = 0
    icnnm = 0
    icntrns = 0
    icntrns2 = 0

    ! itr=0 path
    do i=1,number_atoms
       ia = i
       do j=1,number_atoms   
          jb = j

          if (i/=j) then
             bondl1 = bndlena(i)+bndlena(j)

             call di_algor_distance(coords_cart(1,i),coords_cart(1,j), dist1s)

             if (di_on_root .and. di_iprint>6) write(di_stdout,*) ' ========atoms =============',i,j,bondl1,dist1s

             if (dist1s<=bondl1) then
                if (icnnm(ia)<max_number_of_connections) then
                   icnnm(ia)=icnnm(ia)+1
                   icnnct(icnnm(ia),ia) = jb
 
                   icntrns(1,icnnm(ia),ia) = 0   
                   icntrns(2,icnnm(ia),ia) = 0   
                   icntrns(3,icnnm(ia),ia) = 0   

                   if (di_iprint>6) then
                      write(di_stdout,'(a,2i5,2i5,4f10.3)') 'Connectcluster',ia,jb,icnnm(ia),icnnct(icnnm(ia),ia),dist1s
                   endif
                else
                   write(string,'(a,i5,a,i5)') ' Error: atom ',i,' has more bonds then',max_number_of_connections
                   call deloc_utils_io_abort('Error in di_algor_make_bond_list -'//string)
                endif
             endif
          endif
       enddo
    enddo


    ! first cell
    m1=1
    m2=1
    m3=1

    do i=1,number_atoms
       inmj=icnnm(i)
       ia = i

       do j=1,number_atoms   
          jb = j
          bondl1 = bndlena(i)+bndlena(j)

          if (di_iprint>6 .and. di_on_root) &
           &  write(di_stdout,*) ' ========atoms =============',i,j,bondl1*length_conv

          do itra=-m1,m1
             do itrb=-m2,m2
                do itrc=-m3,m3
                !  get interatomic distance I-J

                   if (itra==0 .and. itrb==0 .and. itrc==0 .and. ia==jb) cycle

                   call di_algor_unfold_atom_coords(coords_cart(1,j),itra,itrb,itrc,tt)

                   call di_algor_distance(coords_cart(1,i),tt,dist2s)

                   if (di_iprint>6 .and. di_on_root) write(di_stdout,'(a,5i5,f20.14)')'Dist',i,j, &
                     &  itra,itrb,itrc,dist2s*length_conv
                   ! new bond is dist2s?
                   if (dist2s<=bondl1) then

                      ! check if such bond to atom ww  was already stored in icnnct

                      skip =.false.
                      do jj= 1,icnnm(ia)
                         ita = icntrns(1,jj,ia)
                         itb = icntrns(2,jj,ia)
                         itc = icntrns(3,jj,ia)
                         icn = icnnct(jj,ia)
      
                         if (jb==icn .and. ita==itra .and. itb==itrb .and. itc==itrc) skip=.true.
          
                         if (di_iprint>6 .and. di_on_root) write(di_stdout,'(2i3,6i4)') ita,itb,itc,itra,itrb,itrc
                      enddo

                      if (.not.skip) then
                         if (icnnm(ia)<max_number_of_connections) then
                            icnnm(ia)=icnnm(ia)+1
                            icnnct(icnnm(ia),ia) = jb
                            icntrns(1,icnnm(ia),ia) = itra
                            icntrns(2,icnnm(ia),ia) = itrb
                            icntrns(3,icnnm(ia),ia) = itrc

                            if (di_iprint>6 .and. di_on_root) then
                               write(di_stdout,'(a,i5,2i5,3i5,f10.3)') 'Connect*******',ia,icnnm(ia), &
                                 &  icnnct(icnnm(ia),ia),itra,itrb,itrc,dist2s*length_conv
                            endif ! di_iprint>4
                         else
                            write(string,'(a,i5,a,i5)') ' atom ',i,' has more bonds then',max_number_of_connections
                            call deloc_utils_io_abort('Error in di_algor_make_bond_list -'//string)
                         endif

                      endif ! not.skip

                   endif ! (bondl1>dist2s) 

                enddo ! itrc
             enddo ! itrb
          enddo ! itra

          if (di_iprint>6 .and. di_on_root) write(di_stdout,*) 'iab ',ia,jb,icnnm(ia),inmj
       enddo ! J

       inmj = icnnm(ia)+1

       ! possibly truncate additional images
       ! pruning time

    enddo ! I

    ! what is the max max_number_of_connections
    ii = maxval(icnnm)
    if (di_iprint>4 .and. di_on_root) write(di_stdout,*) 'Make unique bonds ',max_number_of_connections,ii

    ! print connectivities

        
    do ia=1,number_atoms
       if (di_iprint>6 .and. di_on_root) write(di_stdout,'(/a,3f8.3)') 'atom ', &
            & coords_cart(1,ia)*length_conv,coords_cart(2,ia)*length_conv, &
            & coords_cart(3,ia)*length_conv
       iii=0
       do jj=1,icnnm(ia)
          itra =  icntrns(1,jj,ia) 
          itrb =  icntrns(2,jj,ia) 
          itrc =  icntrns(3,jj,ia) 
          jb = icnnct(jj,ia)

          skip=.false.
          if ((itra+itrb+itrc)>0) then
             ! a possibility of unnecesary bond
             !  (ia,0,0,0) and (jb,itra,itrb,itrc)

             !  do we have another (jb2,ita,itb,itc) so itra-ita=0...?
             do i=1,jj-1
                ita = icntrns(1,i,ia) 
                itb = icntrns(2,i,ia) 
                itc = icntrns(3,i,ia) 
                jb2 = icnnct(i,ia)
                if (jb==2 .and. (itra==ita .and. itrb==itb .and. itrc==itc)) skip=.true.
             enddo

             if (ia>jb) then
             !  a possibility that ia,jb was already found
                do i=1,icnnm(jb)
                   ita =  icntrns(1,i,jb) 
                   itb =  icntrns(2,i,jb) 
                   itc =  icntrns(3,i,jb) 
                   jb2 = icnnct(i,jb)
                   if (ia==jb2 .and. (itra==-ita .and. itrb==-itb .and. itrc==-itc)) skip=.true.
                enddo
             endif

          endif ! (itra+itrb+itrc)>0
       
          if (di_iprint>6 .and. di_on_root) write(di_stdout,*) 'itrabc ', icnnct(jj,ia),itra,itrb,itrc,skip

          if (.not.skip) then
             iii=iii+1
             icnnct2(iii) = icnnct(jj,ia)
             do i=1,3
                icntrns2(i,iii)=icntrns(i,jj,ia)
             enddo
          endif

       enddo ! jj=1,icnnm(ia)

       ! redefine icnnm
       icnnm(ia)=iii
       do jj=1,icnnm(ia)
          icnnct(jj,ia) = icnnct2(jj)
          do i=1,3
             icntrns(i,jj,ia)=icntrns2(i,jj)
          enddo
       enddo

    enddo ! ia=1,number_atoms

    do ia=1,number_atoms
       if (di_iprint>6 .and. di_on_root) write(di_stdout,'(/a,i3,3f8.3)') 'atom ',ia, &
            &  coords_cart(1,ia)*length_conv,coords_cart(2,ia)*length_conv, &
            &  coords_cart(3,ia)*length_conv

       do jj=1,icnnm(ia)
          itra =  icntrns(1,jj,ia) 
          itrb =  icntrns(2,jj,ia) 
          itrc =  icntrns(3,jj,ia) 
          jb = icnnct(jj,ia)

          call di_algor_unfold_atom_coords(coords_cart(1,jb),itra,itrb,itrc,tt)

          call di_algor_distance(coords_cart(1,ia),tt,dist2s)

          if (di_iprint>6 .and. di_on_root) then
             tt=tt*length_conv
  
             write(di_stdout,'(i4,3i3,3f8.3,f10.3)') icnnct(jj,ia),itra,itrb,itrc, &
         & tt(1),tt(2),tt(3),dist2s*length_conv
          endif

       enddo ! jj=1,icnnm(ia)

    enddo ! ia=1,number_atoms

    ! the final list
    list = 0
    do ia=1,number_atoms
       if (di_iprint>6 .and. di_on_root) write(di_stdout,'(/a,i3,3f8.3)') 'atom ', ia, &
            & coords_cart(1,ia)*length_conv,coords_cart(2,ia)*length_conv, &
            & coords_cart(3,ia)*length_conv

       do jj=1,icnnm(ia)
          itra =  icntrns(1,jj,ia) 
          itrb =  icntrns(2,jj,ia) 
          itrc =  icntrns(3,jj,ia) 
          jb   =  icnnct(jj,ia)

          ! could it be that we have exactly the same bond already stored??
          skip=.false.
          do j=1,list

             ita =  latomb(j,3) 
             itb =  latomb(j,4) 
             itc =  latomb(j,5) 
             jb1 = latomb(j,1)
             jb2 = latomb(j,2)

             if (ia==jb1 .and. jb==jb2) then
                if (itra==ita .and. itrb==itb .and. itrc==itc) skip=.true.
             endif

             if (ia==jb2 .and. jb==jb1) then
                 if (itra==-ita .and. itrb==-itb .and. itrc==-itc) skip=.true.
             endif
          enddo ! j=1,list
   
          if (.not.skip) then

             list = list +1
             if (list>max_number_of_bonds) then
                write(string,*) ' Error: too many bonds',list
                call deloc_utils_io_abort('Error in di_algor_make_bond_list -'//string)
             endif

             latomb(list,1)=ia
             latomb(list,2)=jb
             latomb(list,3)=itra
             latomb(list,4)=itrb
             latomb(list,5)=itrc

             if (di_iprint>4 .and. di_on_root) write(di_stdout,'(a,5i5)') 'latomb ',ia,jb,itra,itrb,itrc
          endif
             
       enddo ! jj=1,icnnm(ia)

    enddo ! ia=1,number_atoms
    
    listb = list
     
    ! print out the list of bonds   
    call di_algor_print_bonds(listb,latomb,latomc,latomcn)
       
    ! analyze latomcn - look for atoms with no bonds, reset their vdW radius scaling
    natsnum=0
    do ia=1,number_atoms
       if (latomcn(ia)==0) then
          natsnum=natsnum+1
          vscal(ia)= scale
        endif
    enddo

    ! clean up
    deallocate(icnnct,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icnnct in di_algor_make_bond_list')
    deallocate(icnnct2,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icnnct2 in di_algor_make_bond_list')
    deallocate(icnnm,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icnnm in di_algor_make_bond_list')
    deallocate(icntrns,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icntrns in di_algor_make_bond_list')
    deallocate(icntrns2,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icntrns2 in di_algor_make_bond_list')

    return
  end subroutine di_algor_make_bond_list 

!------------------------------------------------------------------------------
  subroutine di_algor_combine_bond_lists (listb,nbonds,lad,latomb,icbond)
    !=========================================================================!
    ! Combine bond lists in latomb (from current interatomic distances) and   !
    ! in icbond (comes from, for example, .mdf file)                          !
    !                                                                         !
    !                         [fbondup in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) listb (inout) : number of unique bonds                               !
    ! 2) nbonds (in) : number of bonds in the additional list                 !
    ! 3) lad (out) : number of added bonds                                    !
    ! 4) latomb (inout) : list of atoms involved in each bond (IDs and cells) !
    ! 4) icbond (in) : alternative list of bonds (packed IDs and cells)       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) max_number_of_bonds                                                  !
    ! 2) di_iprint                                                            !
    ! 3) di_on_root                                                           !
    ! 4) di_stdout                                                            !
    ! 5) np_p                                                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(inout)                        :: listb
    integer,intent(in)                           :: nbonds
    integer,intent(out)                          :: lad
    integer,dimension(:,:),intent(inout)         :: latomb
    integer,dimension(:,:),intent(in)            :: icbond
    !-------------------------------------------------------------------------!
    logical :: skip
    integer :: list,list0,n,ia,itra,itrb,itrc,jb,ia1,ia2,ia3,ib1,ib2,ib3,j
    integer :: ita,itb,itc,jb1,jb2,ll
    character (len=80) :: string
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(latomb,1)<max_number_of_bonds .and. size(latomb,2)<5) &
       &  call deloc_utils_io_abort('Error in di_algor_combine_bond_lists : Array LATOMB is too small')
    if (size(icbond,1)<nbonds .and. size(icbond,2)<2) &
       &  call deloc_utils_io_abort('Error in di_algor_combine_bond_lists : Array ICBOND is too small')
#endif

    list = listb
    list0 = list

    ! get Nbond

    do n=1,nbonds

       ia = icbond(n,1)
       call deloc_utils_unpack_atom_indices(icbond(n,2), itra, itrb, itrc, jb)

       ! correct itra,itrb,itrc in case deloc_utils_rationalize_coords was called
       ia1 = np_p(1,ia)
       ia2 = np_p(2,ia)
       ia3 = np_p(3,ia)
       ib1 = np_p(1,jb) + itra
       ib2 = np_p(2,jb) + itrb
       ib3 = np_p(3,jb) + itrc
       itra = ib1 - ia1
       itrb = ib2 - ia2
       itrc = ib3 - ia3

       ! a pair  ia, jb, itra, itrb, itrc: is this anywhere in the set of latomb?
       skip = .false.

       do j=1,list

          ita = latomb(j,3)
          itb = latomb(j,4)
          itc = latomb(j,5)
          jb1 = latomb(j,1)
          jb2 = latomb(j,2)
          if (ia==jb1 .and. jb==jb2) then
             if (itra==ita .and. itrb==itb .and. itrc==itc) skip=.true.
          endif
          if (ia==jb2 .and. jb==jb1) then
             if (itra==-ita .and. itrb==-itb .and. itrc==-itc) skip=.true.
          endif
       enddo

       if (.not.skip) then

          list = list +1
          if (list>max_number_of_bonds) then
             write(string,*) ' too many bonds',num_primitive_bonds
             call deloc_utils_io_abort('Error in di_algor_combine_bond_lists - '//string)
          endif

          if (ia>=jb) then
             latomb(list,1)=jb
             latomb(list,2)=ia
             latomb(list,3)=-itra
             latomb(list,4)=-itrb
             latomb(list,5)=-itrc
          elseif (ia<jb) then
             latomb(list,1)=ia
             latomb(list,2)=jb
             latomb(list,3)=itra
             latomb(list,4)=itrb
             latomb(list,5)=itrc
          endif

          if (di_iprint>=8 .and. di_on_root) write(di_stdout,'(a,5i5)') 'latomb= ',ia,jb,itra,itrb,itrc

       end if
    enddo ! nbonds

    listb = list

    ll = list-list0
    if (ll/=0) then
       if (di_iprint>=3 .and. di_on_root) write(di_stdout,'(/,i5,a,/)') ll, ' additional bonds found in MDF file'
       do j = list0+1,list
          ita =  latomb(j,3)
          itb =  latomb(j,4)
          itc =  latomb(j,5)
          ia  = latomb(j,1)
          jb  = latomb(j,2)
          if (di_iprint>2 .and. di_on_root) write(di_stdout,'(a,i5,3x,i5,1x,3i5)') 'latomb + ',ia,jb,ita,itb,itc
       enddo
    else
       if (di_iprint>=3 .and. di_on_root) write(di_stdout,'(/,a)') ' No additional bonds found in MDF file'
    endif

    lad = ll

    return
  end subroutine di_algor_combine_bond_lists
!------------------------------------------------------------------------------

  subroutine di_algor_print_bonds(listb,latomb,latomc,latomcn)
    !=========================================================================!
    ! Print number of unique bonds.                                           !
    ! Set latomc, latomcn  mapping                                            !
    !                                                                         !
    !                         [fbondpri in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) listb (in) : number of unique bonds                                  !
    ! 2) latomb (in) : list of atoms involved in each bond (IDs and cells)    !
    ! 3) latomc (out) : list of all connections                               !
    ! 4) latomcn (out) : number of connections on each atom                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) max_number_of_bonds                                                  !
    ! 2) di_iprint                                                            !
    ! 3) di_on_root                                                           !
    ! 4) di_stdout                                                            !
    ! 5) max_number_of_connections                                            !
    ! 6) number_atoms                                                         !
    ! 7) coords_cart                                                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none

    integer,intent(in)                                :: listb
    integer,dimension(:,:),intent(in)                 :: latomb
    integer,dimension(:),intent(out)                  :: latomcn
    integer,dimension(:,:,:),intent(out)              :: latomc
    !-------------------------------------------------------------------------!

    real (kind=di_dp),dimension(1:3) :: tt
    real (kind=di_dp) :: dist2s
    logical :: skip
    integer :: ia,kk,itra,itrb,itrc,jb,jj
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(latomb,1)<max_number_of_bonds .and. size(latomb,2)<5) &
       &  call deloc_utils_io_abort('Error in di_algor_print_bonds : Array LATOMB is too small')
    if (size(latomcn,1)<number_atoms) &
       &  call deloc_utils_io_abort('Error in di_algor_print_bonds : Array LATOMCN is too small')
    if (size(latomb,1)<max_number_of_connections .and. size(latomb,2)<4 .and. size(latomc,3)<number_atoms) &
       &  call deloc_utils_io_abort('Error in di_algor_print_bonds : Array LATOMC is too small')
#endif

    if (di_iprint>=3 .and. di_on_root) then
       write(di_stdout,'(a,i5,/)') 'Total number of unique bonds is ',listb
       write(di_stdout,'(/,a,/)') ' atom #   l1 l2 l3       x       y       z       distance '
    endif

    ! the final print
    do ia=1,number_atoms
       if (di_iprint>=3 .and. di_on_root) write(di_stdout,'(/a,i4,a,3f8.3,/,a,/)') 'atom ', ia, '  with xyz: ', &
                &       coords_cart(1,ia)*length_conv, coords_cart(2,ia)*length_conv, &
                &       coords_cart(3,ia)*length_conv, 'is connected to: '


       kk=0
       do jj=1,listb

          skip = .true.
          if (latomb(jj,2)==ia) then

             itra = -latomb(jj,3)
             itrb = -latomb(jj,4)
             itrc = -latomb(jj,5)
             jb = latomb(jj,1)

             call di_algor_unfold_atom_coords(coords_cart(1,jb),itra,itrb,itrc,tt)

             skip=.false.

          elseif (latomb(jj,1)==ia) then
             itra =  latomb(jj,3)
             itrb =  latomb(jj,4)
             itrc =  latomb(jj,5)
             jb = latomb(jj,2)

             call di_algor_unfold_atom_coords(coords_cart(1,jb),itra,itrb,itrc,tt)

             skip=.false.
          endif

          if (.not.skip) then
             call di_algor_distance(coords_cart(1,ia),tt,dist2s)
             tt = tt*length_conv

             if (di_iprint>=3 .and. di_on_root) write(di_stdout,'(i4,4x,3i3,3x,3f8.3,f10.3)') jb,itra,itrb,itrc, &
                               &      tt(1),tt(2),tt(3),dist2s*length_conv

             kk = kk+1
             latomc(kk,1,ia) = jb
             latomc(kk,2,ia) = itra
             latomc(kk,3,ia) = itrb
             latomc(kk,4,ia) = itrc

             if (latomb(jj,1)==ia .and. latomb(jj,2)==ia) then
                ! the same index do this atom again
                itra = -itra
                itrb = -itrb
                itrc = -itrc
                call di_algor_unfold_atom_coords(coords_cart(1,ia),itra,itrb,itrc,tt)
                call di_algor_distance(coords_cart(1,ia),tt,dist2s)

                tt = tt*length_conv

                if (di_iprint>=3 .and. di_on_root) write(di_stdout,'(i4,4x,3i3,3x,3f8.3,f10.3)') jb,itra,itrb,itrc, &
                                  &      tt(1),tt(2),tt(3),dist2s*length_conv

                kk = kk+1
                latomc(kk,1,ia) = jb
                latomc(kk,2,ia) = itra
                latomc(kk,3,ia) = itrb
                latomc(kk,4,ia) = itrc
             endif  ! the same atoms
          endif ! not skip
       enddo

       latomcn(ia) = kk
    enddo

    return
  end subroutine di_algor_print_bonds
!------------------------------------------------------------------------------


  subroutine di_algor_make_B_matrix(itor,ktyp,bmat,indexb,xprim,lidelp,iexit)
    !=========================================================================!
    ! Make B matrix for bonds, angles, dihedrals, frozen atoms, and           !
    ! disconnected fragments.                                                 !
    !                                                                         !
    !                         [makebpr in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) itor (in) : flag defining saving of torsions                         !
    ! 2) ktyp (inout) : type of primitives (0,1,3,5,6)                        !
    ! 3) bmat (out) : full B matrix                                           !
    ! 4) indexb (out) : B matrix index                                        !
    ! 5) xprim (out) : values of primitive coordinates                        !
    ! 6) lidelp (in) : list of atoms constituting internal coordinates        !
    ! 7) iexit (out)                                                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) intcor                                                               !
    ! 2) di_iprint                                                            !
    ! 3) di_on_root                                                           !
    ! 4) di_stdout                                                            !
    ! 5) num_primitive_bonds                                                  !
    ! 6) num_primitive_angles                                                 !
    ! 7) num_primitive_torsions                                               !
    ! 8) number_atoms                                                         !
    ! 9) num_disconn_fragm                                                    !
    !10) bmat_disconn_fragm                                                   !
    !11) num_fixed_coords                                                     !
    !12) map_fixed_coords                                                     !
    !13) coords_cart                                                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                             :: itor
    integer,dimension(:),intent(inout)                             :: ktyp
    real (kind=di_dp),dimension(:,:),intent(out)                   :: bmat
    integer,dimension(:,:),intent(out)                             :: indexb
    real (kind=di_dp),dimension(:),intent(out)                     :: xprim
    integer,dimension(:,:),intent(in)                              :: lidelp
    integer,intent(out)                                            :: iexit
    !-------------------------------------------------------------------------!
    integer :: i,int1,j,k,ii
    character (len=80) :: string
    real (kind=di_dp) :: s1
    !-------------------------------------------------------------------------!

    iexit=0
          
    ! bonds
    if (num_primitive_bonds<=intcor) then
       call di_algor_make_bonds_B_matrix(bmat,indexb,xprim,lidelp,1,num_primitive_bonds)
       do i=1,num_primitive_bonds
          ktyp(i)=1
       enddo
    else
       write(string,*) ' too many bonds',num_primitive_bonds
       call deloc_utils_io_abort('Error in di_algor_make_B_matrix - '//string)
    endif

    ! angles
    if (num_primitive_bonds+num_primitive_angles<=intcor) then
       call di_algor_make_angles_B_matrix(bmat,indexb,xprim,lidelp,num_primitive_bonds+1,num_primitive_angles+num_primitive_bonds, &
     &            iexit)
       do i=num_primitive_bonds+1,num_primitive_angles+num_primitive_bonds
          ktyp(i)=3
       enddo
    else
       write(string,*) ' too many bonds and angles',num_primitive_bonds,num_primitive_angles
       call deloc_utils_io_abort('Error in di_algor_make_B_matrix - '//string)
    endif

    ! int1 not intcor in case we have constraints
    int1 = num_primitive_angles+num_primitive_bonds+num_primitive_torsions
         
    ! dihedrals
    call di_algor_make_torsions_B_matrix(itor,bmat,indexb,xprim,lidelp,num_primitive_bonds+num_primitive_angles+1,int1,iexit)
    do i=num_primitive_bonds+num_primitive_angles+1,int1     
       ktyp(i)=5
    enddo


    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,*) 'di_algor_make_B_matrix:'
       write(di_stdout,*) 'nuna,nunb,nund,intcor,nfix,itor,nmfree ', num_primitive_angles,num_primitive_bonds, &
       &             num_primitive_torsions,intcor,num_fixed_coords,itor,num_disconn_fragm
    endif

    ! first time
    if (itor==1) return

    !----------DISCONNECTED FRAGMENTS
    if (num_disconn_fragm/=0 .and. itor>1)  then

       ! fill xprim
       do k=1,num_disconn_fragm
          ii=0
          s1 = zero
          do j=1,number_atoms
             do i=1,3
                ii=ii+1
                s1 = s1 + coords_cart(i,j)*bmat_disconn_fragm(ii,k)
             enddo
          enddo
          xprim(k+int1) = s1
       enddo

       int1 = int1 + num_disconn_fragm
    endif  ! num_disconn_fragm
       

    !---------FROZEN coordinates
    if (num_fixed_coords/=0 .and. itor>1) then

    ! fill xprim
       do k=1,num_fixed_coords
          ii=0
          do j=1,number_atoms
             do i=1,3
                ii=ii+1
                if (map_fixed_coords(k)==ii) xprim(k+int1) = coords_cart(i,j)
             enddo
          enddo
       enddo

    endif  ! num_fixed_coords

    return
  end subroutine di_algor_make_B_matrix
!------------------------------------------------------------------------------

  subroutine di_algor_make_B_matrix_old(ktyp,klist,coeff,itor_loc,ierr,b,savtor_loc,coords_loc)
    !=========================================================================!
    ! Calculates internal coordinate values and the B-Matrix                  !
    ! over the primitive internal coordinates.                                !
    !                                                                         !
    !  NOTE: This routine is only called from di_algor_transf_Hess_to_intern  !
    !        (ex hssindp) when the gradient contribution to hessian is        !
    !        evaluated. This functionality is currently unsupported, and thus !
    !        this routine is not being properly used or tested.               !
    !                                                                         !
    !        The same problem as in hssindp: we don't have the "right" klist  !
    !        data array (has to be constructed from the packed LIDELP one)    !
    !                                                                         !
    ! For an easy to read description of the B-matrix see                     !
    ! S.Califano, Vibrational States, Wiley, London 1976                      !
    !                                                                         !
    ! This routine recognizes several different primitive                     !
    ! internal coordinates:                                                   !
    !  1 stre    stretch: a bond distance between 2 atoms                     !
    !  2 logr    logarithmic bond distance (log(r))                           !
    !            (1/r supposedly gives better convergence in certain cases)   !
    !  3 bend    a bond angle involving 3 atoms (a planar bend)               !
    !  4 outp    out-of-plane bend involving 4 atoms                          ! 
    !  5 tors    standard torsion involving 4 atoms                           !
    ! 11 xcar   cartesian X coordinate                                        !
    ! 12 ycar   cartesian Y coordinate                                        !
    ! 13 zcar   cartesian Z coordinate                                        !
    !                                                                         !
    !                        [makebmd in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ktyp (in) -  array indicating each primitive internal type           !   
    !               currently these are                                       !
    !                1 - stretch                                              !
    !                2 - logarithmic stretch                                  !
    !                3 - bend                                                 !
    !                4 - out-of-plane bend                                    !
    !                5 - torsion                                              !
    !                6 - linear coplanar bend (currently not used)            !   
    !                7 - linear perpendicular bend (currently not used)       ! 
    !             8-10 - (currently not used)                                 !
    !               11 - cartesian X coordinate                               !
    !               12 - cartesian Y coordinate                               !
    !               13 - cartesian Z coordinate                               !
    !               14 - rigid X translation coordinate                       !
    !               15 - rigid Y translation coordinate                       !
    !               16 - rigid Z translation coordinate                       ! 
    ! 2) klist (in)  -  list of atoms involved in each primitive              !
    ! 3) coeff (in)  -  weighting of primitive in non-redundant               ! 
    !             natural internal coordinate space                           !
    !             Note: if the weighting of any given primitive is            !
    !                   zero, it can be ignored in B-Matrix construction      ! 
    ! 4) itor_loc (in) -  integer flag for saving/checking values of primitive!    
    !             torsions (this is done due to possible convergence          !
    !             problems in iterative generation of new Cartesians)         !  
    !              0 - no action                                              !
    !              1 - save initial primitive torsions in SavTOR              !
    !              2 - check current torsions against initial                 !
    !                  (may need to change sign of angles > PI)               !
    ! 5) ierr (out)   -  error flag   0 - success                             !
    !                              -1 - something went wrong                  !
    ! 6) B (out,optional)      -  B-Matrix                                    !
    ! 7) savtor_loc (out,optional) -  array for storing primitive torsions    !
    !                           (possible sign changes near limiting values)  !
    ! 8) coords_loc (out,optional) -  internal coordinates                    !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) number_atoms                                                         !
    ! 2) intcor                                                               !
    ! 3) di_on_root                                                           !
    ! 4) di_stdout                                                            !
    ! 5) di_iprint                                                            !
    ! 6) coords_cart                                                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none

    integer,dimension(1:intcor),intent(in)                       :: ktyp
    integer,dimension(1:4,1:intcor),intent(in)                   :: klist
    real (kind=di_dp),dimension(1:intcor),intent(in)             :: coeff
    integer,intent(in)                                           :: itor_loc
    integer,intent(out)                                          :: ierr
    real (kind=di_dp),dimension(:,:),intent(out),optional        :: b
    real (kind=di_dp),dimension(:),intent(out),optional          :: coords_loc
    real (kind=di_dp),dimension(:),intent(out),optional          :: savtor_loc
    !-------------------------------------------------------------------------!
    logical :: make_coords,make_bmatrix,save_torsions

    real(kind=di_dp), external :: ddot,dnrm2
    real (kind=di_dp),dimension(1:3) :: t,u,v,w,x,y,z,ss,tt,yy,uu,vv,ww,zz
    real (kind=di_dp),dimension(1:12) :: uv
    equivalence (uv(1),uu(1)),(uv(4),vv(1)),(uv(7),ww(1)),(uv(10),zz(1))
    real (kind=di_dp),dimension(1:3,1:4) :: xat
    real (kind=di_dp),parameter :: pt5=0.5_di_dp
    real (kind=di_dp),parameter :: two=2.0_di_dp
    real (kind=di_dp),parameter :: small=1.0e-6_di_dp
    integer,dimension(1:19) :: limtyp
    data limtyp / 2,2,3,4,4,4,4,3,3,3,1,1,1,0,0,0,1,1,1 /
    integer i,ntyp,k,j,l,kk,nat3
    real (kind=di_dp) :: co,si,sir1,sir2,qq,r1,r2,r3,wnorm,xnorm,cp,sj,cq,cr,d,st2,st3,cx
    real (kind=di_dp) :: rinv,sjr3,s
    character (len=80) :: string
    !-------------------------------------------------------------------------!
    ierr = -1
    nat3 = 3*number_atoms

    ! set "make" switches
    save_torsions = present(savtor_loc)
    make_bmatrix = present(b)
    make_coords = present(coords_loc)
    if (make_bmatrix) b = zero

    !  Loop over all primitive internal coordinates

    do i=1,intcor

       !  ignore this primitive if its weighting in the
       !  symmetry non-redundant space is zero

       if (coeff(i)==zero) cycle

       !  initialize for this coordinate

       ntyp = limtyp(ktyp(i))

       if (di_iprint>6 .and. di_on_root) then
          write(di_stdout,*) ' Calculating B-Matrix Column:',i

          if (ktyp(i)==3) then
             write(di_stdout,"(i5,2i5,'(',3i5,')',i5,'(',3i5,')')") i,klist(1,i),klist(2,i),klist(3,i)
          endif ! bend

          do k=1,ntyp
             write(di_stdout,'(i5,3f20.14)') klist(k,i),(coords_cart(j,klist(k,i)),j=1,3)
          enddo ! k
       endif

       do k = 1,ntyp
          do j = 1,3
             xat(j,k)=coords_cart(j,klist(k,i))
          enddo ! j
       enddo ! k

       if (ktyp(i)==1) then

          ! ------STRETCH------

          call di_algor_norm_vector_difference(uu,qq,xat(1,1),xat(1,2))
          vv(1) = -uu(1)
          vv(2) = -uu(2)
          vv(3) = -uu(3)

       else if (ktyp(i)==2) then

          ! ------LOG OF STRETCH------

          call di_algor_norm_vector_difference(uu,qq,xat(1,1),xat(1,2))
          rinv = one/qq
          qq = log(qq)
          uu(1) = uu(1)*rinv
          uu(2) = uu(2)*rinv
          uu(3) = uu(3)*rinv
          vv(1) = -uu(1)
          vv(2) = -uu(2)
          vv(3) = -uu(3)

       else if (ktyp(i)==3) then

          ! ------BEND------

          call di_algor_norm_vector_difference(u,r1,xat(1,1),xat(1,3))
          call di_algor_norm_vector_difference(v,r2,xat(1,2),xat(1,3))

          co = ddot(3,u(1),1,v(1),1)
          si = di_algor_s2(co)
          sir1 = si*r1
          sir2 = si*r2

          if (abs(sir1)<small .or. abs(sir2)<small) then
             if (di_iprint>3 .and. di_on_root) then
                write(di_stdout,2000) i
                write(di_stdout,2200)
                write(di_stdout,2100) i
                write(di_stdout,'(3f20.14)')((xat(j,k),j=1,3),k=1,4)
                write(di_stdout,'(/,a,/)')'Try starting the optimization again.'
             endif
             write(string,2000) i
             call deloc_utils_io_abort('Error in di_algor_make_B_matrix_old:'//string)
             return
          endif

          do l=1,3
             uu(l) = (co*u(l) - v(l))/sir1
             vv(l) = (co*v(l) - u(l))/sir2
             ww(l) = - (uu(l) + vv(l))
          end do

          If (Abs(co)>one) co = SIGN(one,co)  
          qq = acos(co)

       else if (ktyp(i)==4) then

          ! ------OUT OF PLANE BEND------

          call di_algor_norm_vector_difference(u,r1,xat(1,2),xat(1,1))
          call di_algor_norm_vector_difference(v,r2,xat(1,3),xat(1,1))
          call di_algor_norm_vector_difference(w,r3,xat(1,4),xat(1,1))
          co = ddot(3,v(1),1,w(1),1)
          si = di_algor_s2(co)

          if (si<small) then
             if (di_iprint>3 .and. di_on_root) then
                write(di_stdout,2000) i
                write(di_stdout,2300)
                write(di_stdout,2100) i
                write(di_stdout,'(3f20.14)')((xat(j,k),j=1,3),k=1,4)
                write(di_stdout,'(/,a,/)')'Try starting the optimization again.'
             endif
             write(string,2000) i
             call deloc_utils_io_abort('Error in di_algor_make_B_matrix_old:'//string)
             return
          endif

          call di_algor_norm_vector_product(v,w,z)
          cp = ddot(3,u(1),1,z(1),1)
          sj = di_algor_s2(cp)

          if (sj<small) then
             if (di_iprint>3 .and. di_on_root) then
                write(di_stdout,2000) i
                write(di_stdout,2400)
                write(di_stdout,2100) i
             endif
             write(string,2000) i
             call deloc_utils_io_abort('Error in di_algor_make_B_matrix_old:'//string)
             return
          endif

          cq = ddot(3,w(1),1,u(1),1)
          cr = ddot(3,v(1),1,u(1),1)
          d = sj*si*si
          st2 = (co*cq-cr)/(r2*d)
          st3 = (co*cr-cq)/(r3*d)

          do l=1,3
             ww(l) = z(l)*st2
             zz(l) = z(l)*st3
          end do

          call di_algor_norm_vector_product(z,u,x)
          call di_algor_norm_vector_product(u,x,z)

          do l=1,3
             vv(l) = z(l)/r1
             uu(l) = - (vv(l) + ww(l) + zz(l))
          end do

          cx = -one
          if (cp<zero) cx = one
          if (abs(sj)>one) sj = sign(one,sj)  

          qq = -cx*acos(sj)
 
       else if (ktyp(i)==5) then

          ! ------TORSION------

          call di_algor_norm_vector_difference(u,r1,xat(1,1),xat(1,2))
          call di_algor_norm_vector_difference(v,r2,xat(1,3),xat(1,2))
          call di_algor_norm_vector_difference(w,r3,xat(1,3),xat(1,4))
          co = ddot(3,u(1),1,v(1),1)
          cp = ddot(3,v(1),1,w(1),1)
          si = di_algor_s2(co)
          sj = di_algor_s2(cp)
          sir1 = si*r1
          sjr3 = sj*r3

          if (abs(sir1)<small .or. abs(sjr3)<small) then
             if (di_iprint>3 .and. di_on_root) then
                write(di_stdout,2000) i
                write(di_stdout,2500)
                write(di_stdout,2100) i
             endif
             return
          endif

          call di_algor_norm_vector_product(u,v,z)
          call di_algor_norm_vector_product(w,v,x)
          do l=1,3
             uu(l) = z(l)/sir1
             zz(l) = x(l)/sjr3
             vv(l) = (r1*co/r2 - one)*uu(l) - (r3*cp/r2)*zz(l)
             ww(l) = - (uu(l) + vv(l) + zz(l))
          end do

          co = ddot(3,z(1),1,x(1),1)
          call di_algor_vector_product(z,x,u)
          si = dnrm2(3,u,1)
          cp = ddot(3,u(1),1,v(1),1)

          s = di_algor_arc1(-co,si)
          if (cp<zero) s = -s
          qq = -s

          !    ** saving/checking of primitive torsions **

          if (save_torsions .and. itor_loc==1) then
             savtor_loc(i) = qq
          else if (save_torsions .and. itor_loc==2) then
             if (abs(qq)>one .and. sign(one,qq)/=sign(one,savtor_loc(i))) then
                if (di_iprint>3 .and. di_on_root) then
                   write(di_stdout,*) ' Torsional sign change  primitive torsion',i
                   write(di_stdout,*) '   previous value: ',savtor_loc(i)
                   write(di_stdout,*) ' calculated value: ',qq
                endif
                d = di_pi - abs(qq)
                qq = sign(di_pi+d,savtor_loc(i))
                if (di_iprint>3) write(di_stdout,*) ' reassigned value: ',qq
             endif
          endif

       else if (ktyp(i)==6) then

          ! ------LINEAR COPLANAR BEND------

          call di_algor_norm_vector_difference(u,r1,xat(1,1),xat(1,3))
          call di_algor_norm_vector_difference(v,r2,xat(1,4),xat(1,3))
          call di_algor_norm_vector_difference(x,r3,xat(1,2),xat(1,3))
          co = ddot(3,v(1),1,u(1),1)
          cp = ddot(3,x(1),1,v(1),1)

          if (abs(co)>one) co = sign(one,co)  
          if (abs(cp)>one) cp = sign(one,cp)  

          qq = di_pi - acos(co) - acos(cp)

          call di_algor_norm_vector_product(v,u,w)
          call di_algor_norm_vector_product(u,w,z)
          call di_algor_norm_vector_product(w,v,y)
          call di_algor_norm_vector_product(x,v,w)
          call di_algor_norm_vector_product(w,x,u)
          call di_algor_norm_vector_product(v,w,t)

          !  internal coordinate +Ve if atom b moves towards atom d

          do l=1,3
             uu(l) = z(l)/r1
             vv(l) = u(l)/r3
             zz(l) = (y(l)+t(l))/r2
             ww(l) = - (uu(l) + vv(l) + zz(l))
          end do

       else if (ktyp(i)==7) then

          ! ------LINEAR PERPENDICULAR BEND------

          call di_algor_norm_vector_difference(u,r1,xat(1,1),xat(1,3))
          call di_algor_norm_vector_difference(v,r2,xat(1,4),xat(1,3))
          call di_algor_norm_vector_difference(z,r3,xat(1,2),xat(1,3))
          call di_algor_vector_product(v,u,w)
          wnorm = dnrm2(3,w,1)

          if (wnorm<small) then
             if (di_iprint>3 .and. di_on_root) then
                write(di_stdout,2000) i
                write(di_stdout,2600)
                write(di_stdout,2100) i
             endif
             return
          endif

          wnorm = one/wnorm
          call dscal(3,wnorm,w(1),1)
          call di_algor_vector_product(z,v,x)
          xnorm = dnrm2(3,x,1)

          if (xnorm<small) then
             if (di_iprint>3 .and. di_on_root) then
                write(di_stdout,2000) i
                write(di_stdout,2600)
                write(di_stdout,2100) i
             endif
             return
          endif

          xnorm = one/xnorm
          call dscal(3,xnorm,x(1),1)
          co = ddot(3,u(1),1,x(1),1)
          cp = ddot(3,z(1),1,w(1),1)

          if (abs(co)>one) co = sign(one,co)  
          if (abs(cp)>one) cp = sign(one,cp)  

          qq = pt5*(di_pi - acos(co) - acos(cp))

          si = two*di_algor_s2(co)
          sj = two*di_algor_s2(cp)

          if (si<small .or. sj<small) then
             if (di_iprint>3 .and. di_on_root) then
                write(di_stdout,2000) i
                write(di_stdout,2700)
                write(di_stdout,2100) i
             endif
             return
          endif

          call di_algor_vector_product(z,w,y)
          call di_algor_vector_product(w,y,t)
          call di_algor_vector_product(u,t,y)
          call di_algor_vector_product(v,t,yy)
          call di_algor_vector_product(u,x,t)
          call di_algor_vector_product(x,t,ss)
          call di_algor_vector_product(z,ss,t)
          call di_algor_vector_product(v,ss,tt)
          do l=1,3
             uu(l) = x(l)/(r1*si) - yy(l)*wnorm/(r1*sj)
             vv(l) = w(l)/(r3*sj) + tt(l)*xnorm/(r3*si)
             zz(l) = y(l)*wnorm/(r2*sj) - t(l)*xnorm/(r2*si)
          end do
          call di_algor_vector_product(uu,u,t)
          call di_algor_vector_product(u,t,uu)
          call di_algor_vector_product(vv,z,t)
          call di_algor_vector_product(z,t,vv)
          call di_algor_vector_product(zz,v,t)
          call di_algor_vector_product(v,t,zz)
          do l=1,3
             ww(l) = - (uu(l) + vv(l) + zz(l))
          end do

       else if (ktyp(i)==11) then

          ! ------CARTESIAN X COORDINATE------

          qq = xat(1,1)
          uu(1)=one
          uu(2)=zero
          uu(3)=zero

       else if (ktyp(i)==12) then

          ! ------CARTESIAN Y COORDINATE------

          qq = xat(2,1)
          uu(1)=zero
          uu(2)=one
          uu(3)=zero
 
       else if (ktyp(i)==13) then

          ! ------CARTESIAN Z COORDINATE------

          qq = xat(3,1)
          uu(1)=zero
          uu(2)=zero
          uu(3)=one

       else if (ktyp(i)==14) then

          ! ------RIGID X TRANSLATION COORDINATE------

          qq = zero
          do k=1,number_atoms
             qq = qq+coords_cart(1,k)
          enddo ! k
          if (make_coords) coords_loc(i) = qq
          if (make_bmatrix) then
             do k=1,number_atoms
                b(2*(k-1)+1,i)=one
             enddo ! k
          endif ! make_bmatrix

       else if (ktyp(i)==15) then

          ! ------RIGID Y TRANSLATION COORDINATE------

          qq = zero
          do k=1,number_atoms
             qq = qq+coords_cart(2,k)
          enddo ! k
          if (make_coords) coords_loc(i) = qq
          if (make_bmatrix) then
             do k=1,number_atoms
                b(2*(k-1)+2,i)=one
             enddo ! k
          endif ! make_bmatrix

       else if (ktyp(i)==16) then

          ! ------RIGID Z TRANSLATION COORDINATE------

          qq = zero
          do k=1,number_atoms
             qq = qq+coords_cart(3,k)
          enddo ! k
          if (make_coords) coords_loc(i) = qq
          if (make_bmatrix) then
             do k=1,number_atoms
                b(2*(k-1)+3,i)=one
             enddo ! k
          endif ! make_bmatrix

       else

          !  Unknown Internal coordinate type

          if (di_on_root) then
             write(di_stdout,2000) i
             write(di_stdout,2800) ktyp(i), (klist(j,i),j=1,4)
          end if
          return

       endif

       if ((ktyp(i)>=14) .and. (ktyp(i)<=19)) cycle ! B matrix already done

       !  Set up column of B-Matrix for this coordinate

       if (make_coords) coords_loc(i) = qq

       if (make_bmatrix) then
          do k=1,ntyp
             kk = 3*(k-1)
             l = 3*(klist(k,i)-1)
             b(l+1,i) = b(l+1,i) + uv(kk+1)
             b(l+2,i) = b(l+2,i) + uv(kk+2)
             b(l+3,i) = b(l+3,i) + uv(kk+3)
          end do
       endif
    end do  ! i=1,intcor

    !  end of loop over internal coordinates

    if  (make_bmatrix .and. di_iprint>6 .and. di_on_root) then
       write(di_stdout,1000)
       call di_algor_print_matrix(intcor,nat3,intcor,b)
    endif

    ierr = 0

    return

 1000 format(/,6x,'B-Matrix for Primitive Internals',/)
 2000 format(/,2x,'***ERROR*** B-Matrix Construction  Internal Coordinate ',i3)
 2100 format(5x,'Internal Coordinate ',i3,' became Ill-Conditioned')
 2200 format(5x,'Bending coordinate is virtually Linear')
 2300 format(5x,'Atoms in base plane for out-of-plane bend are virtually linear')
 2400 format(5x,'Out-of-plane bend ill-defined')
 2500 format(5x,'Three or more atoms for torsion are virtually linear')
 2600 format(5x,'Linear perpendicular bend is degenerate')
 2700 format(5x,'Linear perpendicular bend ill-defined')
 2800 format(5x,'Unknown internal coordinate type',i5,1x,4i5)

  end subroutine di_algor_make_B_matrix_old
!------------------------------------------------------------------------------

  subroutine di_algor_initialize_H_internal(ktyp,ut)
    !=========================================================================!
    ! Calculate Hessian by forming UT(t)*H(prim)*UT                           !
    ! where H(prim) is a "Hessian" diagonal matrix in the                     !
    ! primitive space with appropriate force constants                        !
    ! given for each primitive type                                           !
    !                                                                         !
    !                        [defhedp in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ktyp (in) -  array indicating each primitive internal type           !
    !                 currently these are                                     !
    !                1 - stretch                                              !
    !                3 - bend                                                 !
    !                4 - out-of-plane bend                                    !
    !                5 - torsion                                              !
    ! 2) ut  (in)  -  set of active delocalized internal coordinates          !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) ndeg                                                                 !
    ! 2) intcor                                                               !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) hess_DI                                                              !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    real (kind=di_dp),dimension(:,:),intent(in)       :: ut
    integer,dimension(:),intent(in)                   :: ktyp
    !-------------------------------------------------------------------------!
    real(kind=di_dp),dimension(1:13) :: cnvrt
    real(kind=di_dp) :: val
    integer :: i,j,k
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(ut,1)<intcor .or. size(ut,2)<ns) call deloc_utils_io_abort('di_algor_initialize_H_internal: Array UT is too small')
    if (size(ktyp,1)<intcor) call deloc_utils_io_abort('di_algor_initialize_H_internal: Array KTYP is too small')
#endif

    cnvrt = zero
    cnvrt(1) = tp_f_bond
    cnvrt(3) = tp_f_angle
    cnvrt(5) = tp_f_tors
    cnvrt(6) = tp_f_crt 

    !  now form the Hessian
    !  this is simply UT(t)*UT scaled

    do i=1,ns
       do j=1,i
          val = zero
          do k=1,intcor
             if (ktyp(k)>0 .and. ktyp(k)<13) val = val + ut(k,i)*ut(k,j)*cnvrt(ktyp(k))
          end do
          hess_DI(i,j) = val
          hess_DI(j,i) = val
       end do
    end do

    return
  end subroutine di_algor_initialize_H_internal
!------------------------------------------------------------------------------

  subroutine di_algor_invert_B_matrix(b,bmbt,binv,ierr)
    !=========================================================================!
    ! "Invert" the B-matrix                                                   !
    ! i.e. set up  B**(-1) = Bt*(B*Bt)**(-1)                                  !
    ! (the m-matrix usual in this formalism is taken as a unit matrix)        !
    !                                                                         !
    !                        [binvtd in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) b (in) -  symmetrized B-matrix                                       !
    ! 2) bmbt (out) -  contains (B*Bt), on exit it has (B*Bt)**(-1)           !
    ! 3) binv (out) -  inverse of B-matrix (is used as scratch in             !
    !                  di_algor_svd_matrix_invert)                            !
    ! 4) ierr (out) -  error flag   0 - success                               !
    !                              -1 - could not invert B-matrix             !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) ndeg                                                                 !
    ! 2) number_atoms                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none

    real (kind=di_dp),dimension(:,:),intent(in)                   :: b
    real (kind=di_dp),dimension(:,:),intent(out)                  :: bmbt
    real (kind=di_dp),dimension(:,:),intent(out)                  :: binv
    integer,intent(out)                                           :: ierr
    !-------------------------------------------------------------------------!
    real(kind=di_dp), external :: ddot

    real(kind=di_dp),parameter :: thrsh=1.0e-5_di_dp
    real(kind=di_dp) :: sum
    integer :: i,j,k,nat3
    !-------------------------------------------------------------------------!
    nat3 = 3*number_atoms
#ifdef debug
    if (size(b,1)<nat3 .or. size(b,2)<ndeg) call deloc_utils_io_abort('Error in di_algor_invert_B_matrix: Array B is too small')
    if (size(binv,1)<nat3 .or. size(binv,2)<ndeg) &
       & call deloc_utils_io_abort('Error in di_algor_invert_B_matrix: Array BINV is too small')
    if (size(bmbt,1)<ndeg .or. size(bmbt,2)<ndeg) &
       & call deloc_utils_io_abort('Error in di_algor_invert_B_matrix: Array BMBT is too small')
#endif

    ierr = -1

    !  set up (B*Bt)

    do i=1,ndeg
       do j=1,i
          bmbt(i,j) = ddot(nat3,b(1,i),1,b(1,j),1)
          bmbt(j,i) = bmbt(i,j)
       end do
    end do

    ! now invert it

    if (ndeg==1) then
       if (abs(bmbt(1,1))<thrsh) go to 999
       bmbt(1,1) = one/bmbt(1,1)
       ierr = 0
    else

       call di_algor_svd_matrix_invert(bmbt,ndeg,ierr)
       if (ierr/=0) go to 999
    endif

    !  set up B**(-1) = Bt*(B*Bt)**(-1)

    do i=1,nat3
       do j=1,ndeg
          sum = zero
          do k=1,ndeg
             sum = sum + b(i,k)*bmbt(k,j)
          end do
          binv(i,j) = sum
       end do
    end do
    
 999 continue

    return
  end subroutine di_algor_invert_B_matrix
!------------------------------------------------------------------------------

  subroutine di_algor_svd_matrix_invert(a,n,ierr)
    !=========================================================================!
    ! Invert A matrix using Singular value decomposition algorithm            !
    !                                                                         !
    !                        [invmatd in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) A (inout)  -  input matrix (returns BB(T)-1  ==A-1)                  !
    ! 2) N (in)     -  dimension of A                                         !
    ! 3) IErr (out) -  error flag                                             !
    !             0 - matrix successfully inverted                            !
    !            -1 - SVD decomposition failed                                !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                 :: n
    real (kind=di_dp),dimension(1:n,1:n),intent(inout) :: a
    integer,intent(out)                                :: ierr

    integer :: lw,i,k,j,istat
    real (kind=di_dp),dimension(:),allocatable         :: w,s
    real (kind=di_dp),dimension(:,:),allocatable       :: u,v
    !-------------------------------------------------------------------------!
    ! Relative Threshold of Zero
    real(kind=di_dp), parameter :: svdthr = 1.0e-14_di_dp
    real(kind=di_dp) :: ss
    !-------------------------------------------------------------------------!

    !  Singular value decomposition:

    !       A = U * S * transpose(V)
    !

    lw = 6*n
    allocate(w(1:lw),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating W in di_algor_svd_matrix_invert')
    w = zero
    allocate(s(1:n),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating S in di_algor_svd_matrix_invert')
    s = zero
    allocate(u(1:n,1:n),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating V in di_algor_svd_matrix_invert')
    u = zero
    allocate(v(1:n,1:n),stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in allocating V in di_algor_svd_matrix_invert')
    v = zero

    call dgesvd('A','A', n,n, a,n, s,u,n, v,n, w,lw, ierr)
    if (ierr/=0) go to 999

    !  u*u(T) = 1;   v*v(T) = 1

    !   calculate A-1
    !      A^(-1) = V*S^(-1)*U(T)

    do i=1,n
       if (s(i)<s(1)*svdthr ) then
          if (di_iprint>4 .and. di_on_root) write(di_stdout,*) 'small diagonal element of bb(t) ',i,s(i)
          s(i)=zero
       endif
    enddo

    !  find inverse of A  by columns
    do k = 1,n
       do j = 1,n

          ! 1/s * U(T)
          if (s(j)/=zero) then
              w(j) = u(k,j)/s(j)
          else
              w(j) = zero
          endif
       end do

       ! V* A(J,K)
       do j=1,n
          ss = zero
          do i = 1,n
             ss = ss + v(i,j)*w(i)
          end do
          a(j,k) = ss
       end do

    end do  ! k column

 999 continue

    deallocate (w,stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in deallocating W in di_algor_svd_matrix_invert')
    deallocate (s,stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in deallocating S in di_algor_svd_matrix_invert')
    deallocate (u,stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in deallocating U in di_algor_svd_matrix_invert')
    deallocate (v,stat=istat)
    if (istat/=0) call deloc_utils_io_abort('Error in deallocating V in di_algor_svd_matrix_invert')

    return
  end subroutine di_algor_svd_matrix_invert
!------------------------------------------------------------------------------

  subroutine di_algor_U_transform_B_matrix(ut,vm,b,indexb)
    !=========================================================================!
    ! Transformation with U for B, X, G, making  B-1                          !
    !                                                                         !
    !                        [tranbmp in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ut (in) : transformation matrix, i.e. which linear combination of    !
    !              primitive internals make up each compound (natural)        !
    !              internal coordinate                                        !
    ! 2) vm(out) : transformed B matrix                                       !
    ! 3) b (in) : B matrix to be transformed                                  !
    ! 4) indexb (in) :  indexing of B matrix                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) map_fixed_coords                                                     !
    ! 2) num_primitive_bonds                                                  !
    ! 3) num_primitive_angles                                                 !
    ! 4) num_primitive_torsions                                               !
    ! 5) num_disconn_fragm                                                    !
    ! 6) bmat_disconn_fragm                                                   !
    ! 7) num_fixed_coords                                                     !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) hess_DI  (the main result of the routine)                            !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none

    real (kind=di_dp),dimension(:,:),intent(in)               :: ut
    real (kind=di_dp),dimension(:,:),intent(out)              :: vm
    real (kind=di_dp),dimension (:,:),intent(in)              :: b
    integer,dimension (:,:),intent(in)                        :: indexb
    !-------------------------------------------------------------------------!
    real (kind=di_dp),dimension(:),allocatable :: v

    integer :: nat3,j,in,k,i,intc1,intc2,intc3,inn,ierr
    real (kind=di_dp) :: uval
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(ut,1)<intcor .or. size(ut,2)<ndeg) &
      & call deloc_utils_io_abort('Error in di_algor_U_transform_B_matrix: Array UT is too small')
#endif
    nat3 = 3*number_atoms
    allocate(v(1:nat3),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating v in di_algor_U_transform_B_matrix')

    v = zero
    vm = zero

    do j=1,ndeg

       ! primitives
       intc1 = num_primitive_bonds+num_primitive_angles+num_primitive_torsions

       do in =1,intc1
          ! unroll b into v
          do k=1,nat3
             v(k) = zero
          enddo
          do i = 1,12
             if (indexb(i,in)/=0) v(indexb(i,in)) = v(indexb(i,in)) + b(i,in)
          enddo
          uval = ut(in,j)
          do k=1,nat3
             vm(k,j) = vm(k,j) + uval*v(k)   
          enddo
       enddo

       ! free   fragments
       intc2 = intc1 + num_disconn_fragm
       if (num_disconn_fragm/=0) then
          do in = intc1+1,intc2
             inn = in - intc1
             uval = ut(in,j)
             do k=1,nat3
                vm(k,j) = vm(k,j) + uval*bmat_disconn_fragm(k,inn)
             enddo
          enddo
       endif

       ! frozen atoms
       if (num_fixed_coords/=0) then
          intc3 = intc2 + num_fixed_coords   
          do in = intc2+1,intc3
             !  unroll b into v
             do k=1,nat3
                v(k) = zero
             enddo
             inn = in - intc2
             v(map_fixed_coords(inn)) = one
             uval = ut(in,j)
             do k=1,nat3
                vm(k,j) = vm(k,j) + uval*v(k)
             enddo
          enddo
       endif
    end do 

    deallocate(v,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating v in di_algor_U_transform_B_matrix')

    return
  end subroutine di_algor_U_transform_B_matrix
!------------------------------------------------------------------------------

  subroutine di_algor_transf_Hess_to_intern(hess,ktyp,lidelp,coeff,ut,binv)    
    !=========================================================================!
    ! Transformation of cartesian Hessian to Hessian over the set             !
    ! of internal coordinates according to:                                   !
    ! hess_DI = B**(-1)t * (HESS - gradient_DI*dB/dcart) * B**(-1)            !
    !                                                                         !
    !  NOTE: Current version does not treat lp_idb=.true. case correctly.     !
    !        In other words, the gradient term should not be used for now.    !
    !                                                                         !
    !                        [hssindp in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) hess(inout) : hessian in cartesian coordinates.                      !
    !            gradient_DI*dB/cart term is evaluated numerically and added. !
    ! 2) ktyp(inout) : array indicating each primitive internal type          !
    !                1 - stretch                                              !
    !                3 - bend                                                 !
    !                5 - torsion                                              !
    ! 3) lidelp (inout) : list of atoms involved in each primitive internal   !
    ! 4) coeff (in) : weighting of primitive in non-redundant                 !
    !                 delocalized internal coordinate space                   !
    ! 5) ut (in) : transformation matrix, i.e. which linear combination of    !
    !              primitive internals make up each compound (natural)        !
    !              internal coordinate                                        !
    ! 6) binv (in) : current inverse B-Matrix                                 !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) number_atoms                                                         !
    ! 5) intcor                                                               !
    ! 6) ndeg                                                                 !
    ! 7) coords_cart                                                          !
    ! 8) lp_idb                                                               !
    ! 9) ns                                                                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) hess_DI  (the main result of the routine)                            !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !  BT      -  will hold B-Matrix in primitive internals                   !
    !  BP      -  will hold transformed B-Matrix forward step                 !
    !  BM      -  will hold transformed B-Matrix backward step                !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none

    real (kind=di_dp),dimension(:,:),intent(inout)                :: hess
    integer,dimension(:),intent(inout)                            :: ktyp
    integer,dimension(:,:),intent(inout)                          :: lidelp
    real (kind=di_dp),dimension(:),intent(in)                     :: coeff
    real (kind=di_dp),dimension(:,:),intent(in)                   :: ut
    real (kind=di_dp),dimension(:,:),intent(in)                   :: binv
    !-------------------------------------------------------------------------!
    real(kind=di_dp),external :: ddot

    real(kind=di_dp),parameter :: delta=0.0005_di_dp
    integer :: nat3,i,j,k,ierr,ic,itor_loc
    real(kind=di_dp) :: gdbdx
    real (kind=di_dp),dimension(:,:),allocatable :: bt
    real (kind=di_dp),dimension(:,:),allocatable :: bp
    real (kind=di_dp),dimension(:,:),allocatable :: bm
    !-------------------------------------------------------------------------!
    nat3 = 3*number_atoms

    allocate(bt(1:nat3,1:intcor),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating bt in di_algor_transf_Hess_to_intern')
    allocate(bp(1:nat3,1:ndeg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating bp in di_algor_transf_Hess_to_intern')
    allocate(bm(1:nat3,1:ndeg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating bm in di_algor_transf_Hess_to_intern')
    bt=zero
    bp=zero
    bm=zero


    if (di_iprint>2 .and. di_on_root) then
       write(di_stdout,1000)
       if (lp_idb) write(di_stdout,1100)
    endif


    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,*) ' Input Hessian in cartesian coordinates is:',nat3
       call di_algor_print_matrix(nat3,nat3,nat3,hess)
    endif

    ! This branch needs actual writing - it will not work the way
    ! the code is. LIDELP is the wrong array to be passed to di_algor_make_B_matrix_old.

    lp_idb = .false.

    if (lp_idb) then
       if (di_iprint>1 .and. di_on_root) write(di_stdout,1200)

       do k = 1,number_atoms

          do i = 1,3

             !  evaluate the gradient_DI*dB/cart term and incorporate
             !  into HESS - this will be done numerically

             !  forward step

             coords_cart(i,k) = coords_cart(i,k) + delta

             itor_loc = 0
             call di_algor_make_B_matrix_old(ktyp,lidelp,coeff,itor_loc,ierr,b=bt)

             if (ierr/=0) go to 95

             !  transform the original B-Matrix to the new coordinates

             ! call di_algor_U_transform_B_matrix1(ut,bt,bp)
             bp = matmul(bt,ut)

             !  backward step

             coords_cart(i,k) = coords_cart(i,k) - delta - delta

             call di_algor_make_B_matrix_old(ktyp,lidelp,coeff,itor_loc,ierr,b=bt)

             if (ierr/=0) go to 95

             !  transform the original B-Matrix to the new coordinates

             ! call di_algor_U_transform_B_matrix1(ut,bt,bm)
             bm = matmul(bt,ut)


             !  restore coordinate

             coords_cart(i,k) = coords_cart(i,k) + delta

             !  form the derivative contribution

             do j=1,nat3
                gdbdx = zero
                do ic=1,ndeg
                   gdbdx = gdbdx + gradient_DI(ic)*(bp(j,ic) - bm(j,ic))
                end do
                gdbdx = gdbdx/(delta + delta)

                !  add contribution to cartesian Hessian prior
                !  to transformation

                hess(i,j) = hess(i,j) - gdbdx
             end do ! do j=1,nat3

          end do ! do i = 1,3

       end do ! do k = 1,number_atoms

    endif   ! lp_idb

    !  Now do the transformation
    !  Transform columns

    do i=1,nat3
       do j=1,ndeg
          bp(i,j) = ddot(nat3,hess(1,i),1,binv(1,j),1)
       end do
    end do

    !  Transform rows

    do i=1,ndeg
       do j=1,i
          bm(i,j) = ddot(nat3,binv(1,i),1,bp(1,j),1)
       end do
    end do

    do i=1,ns
       do j=1,i
          hess_DI(i,j) = bm(i,j)
          hess_DI(j,i) = hess_DI(i,j)
       end do
    end do

    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,*) ' Hessian in internal coordinates is:',ns
       call di_algor_print_matrix(ns,ns,ns,hess_DI)
    endif

    deallocate(bt,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bt in di_algor_transf_Hess_to_intern')
    deallocate(bp,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bp in di_algor_transf_Hess_to_intern')
    deallocate(bm,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bm in di_algor_transf_Hess_to_intern')

    return

 95 continue
    call deloc_utils_io_abort('Error in di_algor_transf_Hess_to_intern: Error Transforming Hessian to Internal Coordinates')

 1000 format(/,' Transforming Hessian to Internal Coordinates')
 1100 format(' Hessian Transformation does not Include Derivative of B-matrix')
 1200 format(' Hessian Transformation Includes Derivative of B-matrix')

  end subroutine di_algor_transf_Hess_to_intern
!------------------------------------------------------------------------------


  subroutine di_algor_make_bonds_B_matrix(bmat,indexb,xprim,lidelp,nun1,nun2)
    !=========================================================================!
    ! Make B matrix for angles primitives                                     !
    !                                                                         !
    !                        [makebp in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) bmat(inout) : B-matrix                                               !
    ! 2) indexb(inout) : indices of the B-matrix                              !
    ! 3) xprim(inout) : primitive internals                                   !
    ! 4) lidelp(in) : atom IDs and lattice translations for internals         !
    ! 5) nun1(in) : first index (after bonds and angles)                      !
    ! 6) nun2(in) : last index (total number of primitives)                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) coords_cart                                                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none

    real (kind=di_dp),dimension(:,:),intent(inout)     :: bmat       
    integer,dimension(:,:),intent(inout)               :: indexb
    real (kind=di_dp),dimension(:),intent(inout)       :: xprim
    integer,dimension(:,:),intent(in)                  :: lidelp
    integer,intent(in)                                 :: nun1
    integer,intent(in)                                 :: nun2
    !-------------------------------------------------------------------------!

    real (kind=di_dp),dimension(1:3) :: tt,uu,vv
    real (kind=di_dp) :: dist2s,qq
    integer i,ni,i0,j0,k0,ia,jb,iaa,ibb,itra,itrb,itrc
    !-------------------------------------------------------------------------!

    !  STRETCH: a bond distance between 2 atoms

    ! loop over all unique atoms

    do ni=nun1,nun2

       ! zero bmat
       do i=1,12
          bmat(i,ni) = zero
       enddo
       do i=7,12
          indexb(i,ni) = 0    
       enddo

       call deloc_utils_unpack_atom_indices(lidelp(1,ni),  i0,  j0,  k0, ia)
       call deloc_utils_unpack_atom_indices(lidelp(2,ni),itra,itrb,itrc, jb)

       call di_algor_unfold_atom_coords(coords_cart(1,jb),itra,itrb,itrc,tt)

       call di_algor_distance(coords_cart(1,ia),tt,dist2s)

       call di_algor_norm_vector_difference(uu,qq,coords_cart(1,ia),tt)
       vv(1) = -uu(1)
       vv(2) = -uu(2)
       vv(3) = -uu(3)

       if (di_iprint>5 .and. di_on_root) then
          tt = tt*length_conv
          write(di_stdout,'(i4,4x,3i3,3x,3f8.3,f10.3)') jb,itra,itrb,itrc,tt(1),tt(2),tt(3), &
          & dist2s*length_conv
       endif

       !  Set up column of B-Matrix for this coordinate

       xprim(ni) = qq

       if (di_iprint>5)  then
          write(di_stdout,'(a,i4,2i5,2f20.14)')'BB=',ni,ia,jb,qq
       endif

       iaa = 3*(ia-1)
       ibb = 3*(jb-1)

       bmat(1,ni) = bmat(1,ni) + uu(1)
       bmat(2,ni) = bmat(2,ni) + uu(2)
       bmat(3,ni) = bmat(3,ni) + uu(3)
       bmat(4,ni) = bmat(4,ni) + vv(1)
       bmat(5,ni) = bmat(5,ni) + vv(2)
       bmat(6,ni) = bmat(6,ni) + vv(3)

       indexb(1,ni) = iaa+1
       indexb(2,ni) = iaa+2
       indexb(3,ni) = iaa+3
       indexb(4,ni) = ibb+1
       indexb(5,ni) = ibb+2
       indexb(6,ni) = ibb+3

    enddo

    return
  end subroutine di_algor_make_bonds_B_matrix

!------------------------------------------------------------------------------


  subroutine di_algor_make_angles_B_matrix(bmat,indexb,xprim,lidelp,nun1,nun2,iexit)
    !=========================================================================!
    ! Make B matrix for angles primitives                                     !
    !                                                                         !
    !                        [makeap in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) bmat(inout) : B-matrix                                               !
    ! 2) indexb(inout) : indices of the B-matrix                              !
    ! 3) xprim(inout) : primitive internals                                   !
    ! 4) lidelp(in) : atom IDs and lattice translations for internals         !
    ! 5) nun1(in) : first index (after bonds and angles)                      !
    ! 6) nun2(in) : last index (total number of primitives)                   !
    ! 7) iexit(out)                                                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) tp_smlden                                                            !
    ! 5) coords_cart                                                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none

    real (kind=di_dp),dimension(:,:),intent(inout)     :: bmat       
    integer,dimension(:,:),intent(inout)               :: indexb
    real (kind=di_dp),dimension(:),intent(inout)       :: xprim
    integer,dimension(:,:),intent(in)                  :: lidelp
    integer,intent(in)                                 :: nun1
    integer,intent(in)                                 :: nun2
    integer,intent(out)                                :: iexit
    !-------------------------------------------------------------------------!

    real (kind=di_dp),external :: ddot

    real (kind=di_dp),dimension(1:3) :: tt,uu,vv,ww,ss,u,v
    real (kind=di_dp) :: r1,r2,co,si,sir1,sir2,qq,qang
    integer ni,i,jb,jc,ija,ijb,ijc,ika,ikb,ikc,i0,j0,k0,iaa,ibb,icc,ia,l
    !-------------------------------------------------------------------------!

    !  ANGLES

    ! loop over all unique atoms

    do ni=nun1,nun2

       ! zero bmat
       do i=1,12
          bmat(i,ni) = zero
       enddo
       do i=10,12
          indexb(i,ni) = 0 
       enddo

       call deloc_utils_unpack_atom_indices(lidelp(1,ni), ija, ijb, ijc, jb)
       call deloc_utils_unpack_atom_indices(lidelp(2,ni),  i0,  j0,  k0, ia)
       call deloc_utils_unpack_atom_indices(lidelp(3,ni), ika, ikb, ikc, jc)

       call di_algor_unfold_atom_coords(coords_cart(1,jb),ija,ijb,ijc,tt)
       call di_algor_unfold_atom_coords(coords_cart(1,jc),ika,ikb,ikc,ss)

       ! a triangle    ( jb, ia,  jc)
       !         (ija, ijb, ijc, 0,0,0, ika, ikb, ikc)
       if (di_iprint>7 .and. di_on_root)  then
          write(di_stdout,'(a,4i3)') '==',ni,jb,ia,jc
          write(di_stdout,'(a,6i3)') 'ij ',ija,ijb,ijc,ika,ikb,ikc
       endif

       call di_algor_norm_vector_difference(u,r1,tt, coords_cart(1,ia))
       call di_algor_norm_vector_difference(v,r2,ss, coords_cart(1,ia))

       co = ddot(3,u(1),1,v(1),1)
       si = di_algor_s2(co)
       sir1 = si*r1
       sir2 = si*r2


       if (abs(sir1)<tp_smlden .or. abs(sir2)<tp_smlden) then
          if (di_on_root) write(di_stdout,'(a,3i3,a)') 'Angle ',jb,ia,jc,' was skipped'
          iexit=1
          return
       endif

       do l=1,3
          uu(l) = (co*u(l) - v(l))/sir1
          vv(l) = (co*v(l) - u(l))/sir2
          ww(l) = - (uu(l) + vv(l))
       enddo

       If (Abs(co)>one) co = SIGN(one,co)
       qq = acos(co)

       qang = qq*toang

       !  Set up column of B-Matrix for this coordinate

       xprim(ni) = qq

       if (di_iprint>5)  then
         write(di_stdout,'(a,i4,3i5,2f20.14)') 'AA=',ni,jb,ia,jc,qq,qang
       endif

       iaa = 3*(ia-1)
       ibb = 3*(jb-1)
       icc = 3*(jc-1)

       bmat(1,ni) = bmat(1,ni) + ww(1)
       bmat(2,ni) = bmat(2,ni) + ww(2)
       bmat(3,ni) = bmat(3,ni) + ww(3)
       bmat(4,ni) = bmat(4,ni) + uu(1)
       bmat(5,ni) = bmat(5,ni) + uu(2)
       bmat(6,ni) = bmat(6,ni) + uu(3)
       bmat(7,ni) = bmat(7,ni) + vv(1)
       bmat(8,ni) = bmat(8,ni) + vv(2)
       bmat(9,ni) = bmat(9,ni) + vv(3)
       indexb(1,ni) = iaa+1
       indexb(2,ni) = iaa+2
       indexb(3,ni) = iaa+3
       indexb(4,ni) = ibb+1
       indexb(5,ni) = ibb+2
       indexb(6,ni) = ibb+3
       indexb(7,ni) = icc+1
       indexb(8,ni) = icc+2
       indexb(9,ni) = icc+3

    enddo

    return
  end subroutine di_algor_make_angles_B_matrix


!------------------------------------------------------------------------------

  subroutine di_algor_make_torsions_B_matrix(itor,bmat,indexb,xprim,lidelp,nun1,nun2,iexit)
    !=========================================================================!
    ! Make B matrix for dihedral primitives                                   !
    !                                                                         !
    !                        [makedp in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) itor (in) : integer flag for saving/checking values of primitive     !
    !             torsions (this is done due to possible convergence          !
    !             problems in iterative generation of new Cartesians)         !
    !              0 - no action                                              !
    !              1 - save initial primitive torsions in savtor              !
    !              2 - check current torsions against initial                 !
    !                  (may need to change sign of angles > PI)               !
    ! 2) bmat(inout) : B-matrix                                               !
    ! 3) indexb(out) : indices of the B-matrix                                !
    ! 4) xprim(inout) : primitive internals                                   !
    ! 5) lidelp(in) : atom IDs and lattice translations for internals         !
    ! 6) nun1(in) : first index (after bonds and angles)                      !
    ! 7) nun2(in) : last index (total number of primitives)                   !
    ! 8) iexit(out)                                                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) tp_smlden                                                            !
    ! 5) coords_cart                                                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) savtor                                                               !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none

    integer,intent(in)                                 :: itor
    real (kind=di_dp),dimension(:,:),intent(inout)     :: bmat       
    integer,dimension(:,:),intent(inout)               :: indexb
    real (kind=di_dp),dimension(:),intent(inout)       :: xprim    
    integer,dimension(:,:),intent(in)                  :: lidelp
    integer,intent(in)                                 :: nun1
    integer,intent(in)                                 :: nun2
    integer,intent(out)                                :: iexit
    !-------------------------------------------------------------------------!
    real (kind=di_dp),external :: ddot,dnrm2

    real (kind=di_dp),dimension(1:3) :: tt,uu,vv,ww,zz,ss,xx,u,v,w,z,x
    integer :: ichd,i,ni,ika,ikb,ikc,kc,itra,itrb,itrc,jb,ila,ilb,ilc,ld,ia,i0,j0,k0,l
    integer :: iaa,ibb,icc,idd
    real (kind=di_dp) :: co,cp,si,sir1,sjr3,r1,r2,r3,sj,s,qq,d,qang
    !-------------------------------------------------------------------------!


    !  DIHEDRAL

    ! loop over all unique bonds 
    ichd = 0

    do ni=nun1,nun2

       ! zero bmat
       do i=1,12
          bmat(i,ni) = zero
       enddo

       call deloc_utils_unpack_atom_indices(lidelp(1,ni), ika, ikb, ikc, kc)
       call deloc_utils_unpack_atom_indices(lidelp(2,ni),  i0,  j0,  k0, ia)
       call deloc_utils_unpack_atom_indices(lidelp(3,ni),itra,itrb,itrc, jb)
       call deloc_utils_unpack_atom_indices(lidelp(4,ni), ila, ilb, ilc, ld)

       call di_algor_unfold_atom_coords(coords_cart(1,jb),itra,itrb,itrc,tt)

       call di_algor_unfold_atom_coords(coords_cart(1,kc),ika,ikb,ikc,ss)
 
       call di_algor_unfold_atom_coords(coords_cart(1,ld),ila+itra,ilb+itrb,ilc+itrc,xx)

       ! a rectangular ( kc, ia,  jb, ld)
       !         (ika, ikb, ikc, 0,0,0, itra, itrb, itrc, ila, ilb,ilc)
       ! coordinates     ss  coords_cart(1,ia) tt  xx

       if (di_iprint>6 .and. di_on_root)  then
          write(di_stdout,'(a,4i3)') 'j ',kc,ia,jb,ld
          write(di_stdout,'(a,9i3)') 'ij ',ika,ikb,ikc,itra,itrb,itrc,ila+itra,ilb+itrb,ilc+itrc
       endif

       ! ------TORSION------

       call di_algor_norm_vector_difference(u,r1,ss, coords_cart(1,ia))
       call di_algor_norm_vector_difference(v,r2,tt, coords_cart(1,ia))
       call di_algor_norm_vector_difference(w,r3,tt, xx)
       co = ddot(3,u(1),1,v(1),1)
       cp = ddot(3,v(1),1,w(1),1)
       si = di_algor_s2(co)
       sj = di_algor_s2(cp)
       sir1 = si*r1
       sjr3 = sj*r3

       if (abs(sir1)<tp_smlden .or. abs(sjr3)<tp_smlden) then
          iexit=1
          return
       endif

       call di_algor_norm_vector_product(u,v,z)
       call di_algor_norm_vector_product(w,v,x)

       do l=1,3
          uu(l) = z(l)/sir1
          zz(l) = x(l)/sjr3
          vv(l) = (r1*co/r2 - one)*uu(l) - (r3*cp/r2)*zz(l)
          ww(l) = - (uu(l) + vv(l) + zz(l))
       enddo     

       co = ddot(3,z(1),1,x(1),1)
       call di_algor_vector_product(z,x,u)
       si = dnrm2(3,u,1)
       cp = ddot(3,u(1),1,v(1),1)

       s = di_algor_arc1(-co,si)
       if (cp<zero) s = -s
       qq = -s

       !    ** saving/checking of primitive torsions **

       if (itor==1) then
          savtor(ni) = qq
       else if (itor==2) then
          if (abs(qq)>one .and. sign(one,qq)/=sign(one,savtor(ni))) then
             ichd = ichd + 1
             if (di_iprint>2 .and. di_on_root) then
                write(di_stdout,*) ' Torsional sign change  primitive torsion',ni
                write(di_stdout,*) '   previous value: ',savtor(ni)
                write(di_stdout,*) ' calculated value: ',qq
             endif
             d = di_pi - abs(qq)
             qq = sign(di_pi+d,savtor(ni))
             if (di_iprint>1 .and. di_on_root) write(di_stdout,*) ' reassigned value: ',qq
          endif
       endif

       qang = qq*toang

       !  Set up column of B-Matrix for this coordinate

       xprim(ni) = qq

       if (di_iprint>5 .and. itor/=2 .and. di_on_root)  then
         write(di_stdout,'(a,i4,4i5,2f20.14)')'DD=',ni,kc,ia,jb,ld,qq, qang
       endif

       iaa = 3*(kc-1)
       ibb = 3*(ia-1)
       icc = 3*(jb-1)
       idd = 3*(ld-1)

       bmat(1,ni) = bmat(1,ni) + uu(1)
       bmat(2,ni) = bmat(2,ni) + uu(2)
       bmat(3,ni) = bmat(3,ni) + uu(3)
       bmat(4,ni) = bmat(4,ni) + vv(1)
       bmat(5,ni) = bmat(5,ni) + vv(2)
       bmat(6,ni) = bmat(6,ni) + vv(3)
       bmat(7,ni) = bmat(7,ni) + ww(1)
       bmat(8,ni) = bmat(8,ni) + ww(2)
       bmat(9,ni) = bmat(9,ni) + ww(3)
       bmat(10,ni) = bmat(10,ni) + zz(1)
       bmat(11,ni) = bmat(11,ni) + zz(2)
       bmat(12,ni) = bmat(12,ni) + zz(3)
       indexb(1,ni) = iaa+1
       indexb(2,ni) = iaa+2
       indexb(3,ni) = iaa+3
       indexb(4,ni) = ibb+1
       indexb(5,ni) = ibb+2
       indexb(6,ni) = ibb+3
       indexb(7,ni) = icc+1
       indexb(8,ni) = icc+2
       indexb(9,ni) = icc+3
       indexb(10,ni) = idd+1
       indexb(11,ni) = idd+2
       indexb(12,ni) = idd+3

    enddo

    if (di_iprint>1 .and. ichd/=0 .and. di_on_root) then
       write(di_stdout,'(a,i4,a)') ' Torsional sign changed for ',ichd, ' internals'
    endif

    return
  end subroutine di_algor_make_torsions_B_matrix
!------------------------------------------------------------------------------




  subroutine di_algor_make_list_primitives(nic,intcor0,listb,latomb,latomc,latomcn,icc,lidelp)
    !=========================================================================!
    ! Make the final list of bonds, angles, dihedrals out of                  !
    ! di_algor_make_bond_list output.                                         !
    ! Correct in case we have internal constraints.                           !
    !                                                                         !
    !                        [mklipe in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) nic (in) : max dimension for primitive set                           !
    ! 2) intcor0(out) : number of primitives                                  !
    ! 3) listb(in) : list of bonds                                            !
    ! 4) latomb(in) : pointers to bonds                                       !
    ! 5) latomc(in) : pointers to unique atoms                                !
    ! 6) latomcn(in) : unique connections per atom                            !
    ! 7) icc(inout) : list of constraints                                     !
    ! 8) lidelp(out,optional) :atom IDs and lattice translations for internals!
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) ncons                                                                !
    ! 5) tp_benlin                                                            !
    ! 6) coords_cart                                                          !
    ! 7) num_primitive_bonds                                                  !
    ! 8) number_atoms                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) num_primitive_angles                                                 !
    ! 2) num_primitive_torsions                                               !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! first_pass - TRUE  --> only calculate number of primitives              !
    !              FALSE --> save LIDELP as well                              !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none

    integer,intent(in)                                 :: nic
    integer,intent(out)                                :: intcor0
    integer,intent(in)                                 :: listb
    integer,dimension(:,:),intent(in)                  :: latomb
    integer,dimension(:,:,:),intent(in)                :: latomc
    integer,dimension(:),intent(in)                    :: latomcn
    integer,dimension(:,:),intent(inout)               :: icc
    integer,dimension(:,:),intent(out),optional        :: lidelp
    !-------------------------------------------------------------------------!

    real(kind=di_dp), external :: ddot

    logical :: linear
    real(kind=di_dp) :: tt(1:3),ss(1:3),xx(1:3),u(1:3),v(1:3),w(1:3),z(1:3),x(1:3)
    real(kind=di_dp) :: tt1(1:3),ss1(1:3),uu(1:3),vv(1:3),zz(1:3),almostlinear
    real(kind=di_dp) :: planar,co,si,sir1,r1,r2,r3,qq,cp,sj,sjr3,qq1,qq2,s
    integer :: j,ni,jj,ia,itra,itrb,itrc,jb,ibad,iad1,jlist,kk,ika,ikb,ikc,ll,jc
    integer :: inj,ija,ijb,nicmax,ijc,igo,inj1,klist1,jj1,jb1,ija1,ijb1
    integer :: ijc1,kk1_skip,klist2,kk1,llist,kc,ikc1,jc1,ika1,ikb1
    integer :: ld,ila,ilb,ilc,l
    character (len=80) :: string
    logical :: first_pass
    !-------------------------------------------------------------------------!
    ! if the indexing is not requested, it has to be the first pass
    first_pass = .not. present(lidelp)

#ifdef debug
    if (.not.first_pass) then
       if (size(lidelp,1)<4 .or. size(lidelp,2)<nic) &
        & call deloc_utils_io_abort('Error in di_algor_make_list_primitives: Wrong dimensions of LIDELP')
    endif
    if (size(latomb,1)<listb .and. size(latomb,2)/=5) &
      & call deloc_utils_io_abort('Error in di_algor_make_list_primitives: Wrong dimensions of LATOMB')
    if (size(latomcn,1)<number_atoms) &
      & call deloc_utils_io_abort('Error in di_algor_make_list_primitives: Wrong dimensions of LATOMCN')
    if (size(latomc,2)<4 .and. size(latomc,3)<number_atoms) &
      & call deloc_utils_io_abort('Error in di_algor_make_list_primitives: Wrong dimensions of LATOMC')
#endif

    almostlinear = tp_benlin
    planar = di_pi
! 21/Apr/2006 aperlov: nicmax is used only if .not. first_pass, but to satisfy ftncheck: 
    nicmax=-1

    if (.not. first_pass) then

       if (ncons/=0 .and. di_iprint>=4 .and. di_on_root) write(di_stdout,'(/,a)') 'CONSTRAINT INTERNAL COORDINATES'

       ! zero lidelp
       do j=1,4
          do ni=1,nic   
             lidelp(j,ni) = 0
          enddo
       enddo
           
       nicmax = nic
    endif ! .not. first_pass

    ! BOND
    if (first_pass) then

       ! nothing to do here
       ni = listb

    else

       ni = 0  
       do jj=1,listb

          ia   =  latomb(jj,1)  ! atom in the 0 cell
          jb    = latomb(jj,2)  ! atom in the (itra itrb itrc) cell
          itra =  latomb(jj,3)
          itrb =  latomb(jj,4)
          itrc =  latomb(jj,5)

          ni = ni +1
          if (ni>nicmax) then
             write (string,'(i5)') ni
             call deloc_utils_io_abort('Error in di_algor_make_list_primitives - too many primitives '//string)
          endif

          ! store distance

          call deloc_utils_pack_atom_indices(lidelp(1,ni),   0,   0,   0, ia)
          call deloc_utils_pack_atom_indices(lidelp(2,ni),itra,itrb,itrc, jb)

       enddo ! jj

    endif ! .not. first_pass

    !  check if ncons are present in the list lidelp
    !  if not add them to the list and correct num_primitive_bonds     
    if (ncons/=0 .and. .not. first_pass) then
       call di_algor_constr_to_primitives(nic,lidelp,1,ni,icc,almostlinear,1)
    endif
    num_primitive_bonds=ni

    !  ANGLES

    ! loop over all unique atoms

    ibad=0
    iad1 = 0

    do ia=1,number_atoms

       jlist = latomcn(ia)

       inj=0
       do jj=1,jlist

          jb  =  latomc(jj,1,ia)
          ija =  latomc(jj,2,ia)
          ijb =  latomc(jj,3,ia)
          ijc =  latomc(jj,4,ia)

          call di_algor_unfold_atom_coords(coords_cart(1,jb),ija,ijb,ijc,tt)
          call di_algor_norm_vector_difference(u,r1,tt, coords_cart(1,ia))

          do kk = 1,jj-1 

             jc =  latomc(kk,1,ia)
             ika =  latomc(kk,2,ia)
             ikb =  latomc(kk,3,ia)
             ikc =  latomc(kk,4,ia)

             call di_algor_unfold_atom_coords(coords_cart(1,jc),ika,ikb,ikc,ss)
             call di_algor_norm_vector_difference(v,r2,ss, coords_cart(1,ia))

             ! a triangle    ( jb, ia,  jc)
             !         (ija, ijb, ijc, 0,0,0, ika, ikb, ikc)
             if (di_iprint>7 .and. di_on_root)  then
                inj=inj+1
                write(di_stdout,'(a,i3,5x,3i3)') '---- j ',inj, jb,ia,jc
                write(di_stdout,'(a,6i3)') 'ij ',ija,ijb,ijc,ika,ikb,ikc
                write(di_stdout,'(a,3f14.3)') 'xb ',(coords_cart(j,jb),j=1,3)
                write(di_stdout,'(a,3f14.3)') 'xa ',(coords_cart(j,ia),j=1,3)
                write(di_stdout,'(a,3f14.3)') 'xc ',(coords_cart(j,jc),j=1,3)
             endif

             co = ddot(3,u(1),1,v(1),1)
             si = di_algor_s2(co)
             sir1 = si*r1

             if (abs(co)>one) co = sign(one,co)
             qq = acos(abs(co))
 
             linear = .false.
             if (abs(qq)<almostlinear) then
                ibad=ibad+1
                linear = .true.
                ! this is valid only if ia has two bonds ia-jb and ia-jc
                if (jlist>2) linear=.false.

                if (di_iprint>7 .and. di_on_root) then
                   write(di_stdout,2200)
                   write(di_stdout,'(a,3i3,a,f10.3)') 'Angle ',jb,ia,jc,' is skipped', qq*toang
                endif
             else  

                ni = ni +1
                if (.not. first_pass) then

                   if (ni>nicmax) then
                      write (string,'(i5)') ni
                      call deloc_utils_io_abort('Error in di_algor_make_list_primitives - too many primitives '//string)
                   endif

                   ! store angle

                   call deloc_utils_pack_atom_indices(lidelp(1,ni), ija, ijb, ijc, jb)
                   call deloc_utils_pack_atom_indices(lidelp(2,ni),   0,   0,   0, ia)
                   call deloc_utils_pack_atom_indices(lidelp(3,ni), ika, ikb, ikc, jc)
                endif ! .not. first_pass
             endif ! not linear

             if (linear) then
                ! how about including angles of ia with the neighbors of jb,jc
                ! this is valid only if ia has two bonds ia-jb and ia-jc

                ! r1b distance  ia-jb
                ! r2c distance  ia-jc

                do igo=1,2
                   inj1 = 0
                   klist1 = latomcn(jb)
                   do jj1=1,klist1+1

                      if (jj1>1) then
                         jb1  =  latomc(jj1-1,1,jb)
                         ija1 =  latomc(jj1-1,2,jb) + ija
                         ijb1 =  latomc(jj1-1,3,jb) + ijb
                         ijc1 =  latomc(jj1-1,4,jb) + ijc
                      else
                         ! just jb
                         ! admit as few angles as possible 
                         ! if (jb) - (ia) - (neighbour of jc) is successful skip the others
                         ! to maintain symmetry do that only if jj1=1
                         jb1 = jb
                         ija1 =  ija
                         ijb1 =  ijb
                         ijc1 =  ijc
                      endif  !jj1>1

                      if (jb1/=ia) then

                         call di_algor_unfold_atom_coords(coords_cart(1,jb1),ija1,ijb1,ijc1,tt1)
                         ! r11 distance jb-jb1
                         ! r1  distance ia-jb1
                         call di_algor_norm_vector_difference(u,r1,tt1, coords_cart(1,ia))

                         ! if r1 > r1b+r11 skip, if r1=0 skip

                         kk1_skip = 0
                         klist2 = latomcn(jc)
                         do kk1 = 1,klist2+1

                            ! no for jb-ia-jc angle
                            if (kk1_skip==1 .or. jj1==1) go to 202

                            if (kk1>1) then
                               ! neighbours of js
                               jc1  =  latomc(kk1-1,1,jc)
                               ika1 =  latomc(kk1-1,2,jc) + ika
                               ikb1 =  latomc(kk1-1,3,jc) + ikb
                               ikc1 =  latomc(kk1-1,4,jc) + ikc
                            else
                               ! just jc
                               ! admit as few angles as possible 
                               ! if (neighbour of jb) - (ia) - (jc) is successful skip the others
                               ! to maintain symmetry do that only if kk1=1
                               jc1  = jc
                               ika1 = ika
                               ikb1 = ikb
                               ikc1 = ikc
                            endif  !kk1>1

                            if (jc1/=ia) then

                               call di_algor_unfold_atom_coords(coords_cart(1,jc1),ika1,ikb1,ikc1,ss1)
                               ! r22 distance jc-jc1
                               ! r2  distance ia-jb1
                               call di_algor_norm_vector_difference(v,r2,ss1, coords_cart(1,ia))

                               if (di_iprint>7 .and. di_on_root)  then
                                  inj1 = inj1+1
                                  write(di_stdout,'(a,i3,5x,3i3)') '++++ j ',inj1, jb1,ia,jc1
                                  write(di_stdout,'(a,6i3)') 'ij ',ija1,ijb1,ijc1,ika1,ikb1,ikc1
                               endif

                               co = ddot(3,u(1),1,v(1),1)
                               si = di_algor_s2(co)
                               sir1 = si*r1

                               if (abs(co) > one) co = sign( one,co)
                               qq = acos(abs(co))

                               if (abs(qq)<almostlinear) then
                                  if (di_iprint>7 .and. di_on_root) then
                                     write(di_stdout,2200)
                                     write(di_stdout,'(a,3i3,a,f10.3)') 'Angle ',jb1,ia,jc1,' is skipped', qq*toang
                                  endif

                               else  ! qq<almostlinear

                                  ni = ni +1
                                  if (.not. first_pass) then
                                     iad1 = iad1 +1
                                     if (ni>nicmax) then
                                        write (string,'(i5)') ni
                                        call deloc_utils_io_abort &
                                       & ('Error in di_algor_make_list_primitives - too many primitives '//string)
                                     endif

                                     if (kk1==1) kk1_skip=1
                                     ! store angle

                                     call deloc_utils_pack_atom_indices(lidelp(1,ni), ija1, ijb1, ijc1, jb1)
                                     call deloc_utils_pack_atom_indices(lidelp(2,ni),   0,   0,   0, ia)
                                     call deloc_utils_pack_atom_indices(lidelp(3,ni), ika1, ikb1, ikc1, jc1)
                                  endif ! .not. first_pass

                               endif ! qq<almostlinear
                               
                            endif !jc1/=ia

 202                        continue

                         enddo  !kk1

                      endif ! jb1/=ia

                   enddo !jj1

                   ! exchange
                   jb1 = jc
                   jc = jb
                   jb = jb1

                enddo ! igo
 
             endif ! linear

          enddo  !kk
       enddo !jj
    enddo !ia

    !  check if ncons are present in the list lidelp
    !  if not add them to the list and correct num_primitive_angles     
    if (ncons/=0 .and. .not. first_pass) then
       call di_algor_constr_to_primitives(nic,lidelp,num_primitive_bonds+1,ni,icc,almostlinear,3)
    endif ! ncons

    num_primitive_angles = ni-num_primitive_bonds

    if (.not. first_pass) then

       if (di_iprint>=4 .and. di_on_root) then
          write (di_stdout,'(/,a,i5)') 'Maximum number of primitive angles:       ',num_primitive_angles+ibad
          if (ibad/=0) write(di_stdout,2202) ibad 
          if (iad1/=0) write(di_stdout,2203) iad1
       endif
    endif


 2200 format(/,5x,'Bending coordinate is virtually Linear')
 2202 format(i5,'  Linear Bend angles were skipped')
 2203 format(i5,'  Additional bends replacing Linear Bends')

    !  DIHEDRAL
    ! loop over all unique bonds 

    ibad=0

    do jj=1,listb

       ! atom ia has itr 0,0,0
       ! atom jb may have itr different from 0

       ia   =  latomb(jj,1)
       jb   =  latomb(jj,2)
       itra =  latomb(jj,3)
       itrb =  latomb(jj,4)
       itrc =  latomb(jj,5)

       call di_algor_unfold_atom_coords(coords_cart(1,jb),itra,itrb,itrc,tt)

       ! loop over all the neighbours of ia 
       jlist = latomcn(ia)

       inj=0
       do kk=1,jlist

          kc =  latomc(kk,1,ia)
          ika =  latomc(kk,2,ia)
          ikb =  latomc(kk,3,ia)
          ikc =  latomc(kk,4,ia)

          if (kc/=jb .or. itra/=ika .or. itrb/=ikb .or. itrc/=ikc) then

             call di_algor_unfold_atom_coords(coords_cart(1,kc),ika,ikb,ikc,ss)

             ! loop over all the neighbours of jb 
             llist = latomcn(jb)

             do ll=1,llist

                ld =  latomc(ll,1,jb)
                ila =  latomc(ll,2,jb)
                ilb =  latomc(ll,3,jb)
                ilc =  latomc(ll,4,jb)

                if (ld/=ia .or. (ila+itra)/=0 .or. (ilb+itrb)/=0 .or. (ilc+itrc)/=0) then
                   ! lc may not be at zero

                   call di_algor_unfold_atom_coords(coords_cart(1,ld),ila+itra,ilb+itrb,ilc+itrc,xx)

                   ! a rectangular ( kc, ia,  jb, ld)
                   !         (ika, ikb, ikc, 0,0,0, itra, itrb, itrc, ila, ilb,ilc)
                   ! coordinates     ss  coords_cart(1,ia) tt  xx

                   if (di_iprint>7 .and. di_on_root)  then
                      inj=inj+1
                      write(di_stdout,'(a,i3,5x,4i3)') '----j ',inj, kc,ia,jb,ld
                      write(di_stdout,'(a,9i3)') 'ij ',ika,ikb,ikc,itra,itrb,itrc,ila+itra,ilb+itrb,ilc+itrc
                   endif  ! print


                   ! ------TORSION------

                   call di_algor_norm_vector_difference(u,r1,ss, coords_cart(1,ia))
                   call di_algor_norm_vector_difference(v,r2,tt, coords_cart(1,ia))
                   call di_algor_norm_vector_difference(w,r3,tt, xx)
                   co = ddot(3,u(1),1,v(1),1)
                   cp = ddot(3,v(1),1,w(1),1)
                   si = di_algor_s2(co)
                   sj = di_algor_s2(cp)
                   sir1 = si*r1
                   sjr3 = sj*r3

                   ! check for planar bends

                   !       --compute value of first planar bend (kc,ia,jb)
                   if (abs(co)> one ) co = sign( one, co)
                   qq1 = acos(abs(co))

                   !       --compute value of second planar bend (ia,jb,ld)
                   if (abs(cp)> one ) cp = sign( one, cp)
                   qq2 = acos(abs(cp))

                   if ( (abs(qq1)<almostlinear) .or. (abs(qq2)<almostlinear) )  then
                      ibad=ibad+1
                      if (di_iprint>7 .and. di_on_root) then
                         write(di_stdout,2500)
                         write(di_stdout,'(a,4i3,a,2f10.3)') 'Dihedral ',kc,ia,jb,ld,' is skipped',qq1*toang,qq2*toang
                      endif  ! print

                   else

                      ! check for dihedral=180, eliminate

                      call di_algor_norm_vector_product(u,v,z)
                      call di_algor_norm_vector_product(w,v,x)
                      do l=1,3
                         uu(l) = z(l)/sir1
                         zz(l) = x(l)/sjr3
                         vv(l) = (r1*co/r2 - one)*uu(l) - (r3*cp/r2)*zz(l)
                      enddo

                      co = ddot(3,z(1),1,x(1),1)
                      call di_algor_vector_product(z,x,u)
                      si = sqrt(ddot(3,u(1),1,u(1),1))
                      cp = ddot(3,u(1),1,v(1),1)

                      s = di_algor_arc1(-co,si)
                      if (cp<zero) s = -s
                      qq = -s

                      ! check for dihedral=0, eliminate
    
                      ni = ni +1
                      if (.not. first_pass) then
                         if (ni>nicmax) then
                            write (string,'(i5)') ni
                            call deloc_utils_io_abort('Error in di_algor_make_list_primitives - too many primitives '//string)
                         endif

                         ! store this dihedral

                         call deloc_utils_pack_atom_indices(lidelp(1,ni), ika, ikb, ikc, kc)
                         call deloc_utils_pack_atom_indices(lidelp(2,ni),   0,   0,   0, ia)
                         call deloc_utils_pack_atom_indices(lidelp(3,ni),itra,itrb,itrc, jb)
                         call deloc_utils_pack_atom_indices(lidelp(4,ni), ila, ilb, ilc, ld)

                      endif ! .not. first_pass

                   endif   ! planar bonds

                endif ! ld/=ia

             enddo  ! ll

          endif ! kc/=jb

       enddo  ! kk

    enddo  ! jj

    !  check if ncons are present in the list lidelp
    !  if not add them to the list and correct num_primitive_angles     
    if (ncons/=0 .and. .not. first_pass) then
        call di_algor_constr_to_primitives(nic,lidelp,num_primitive_bonds+num_primitive_angles+1,ni,icc,almostlinear,5)
    endif ! ncons

    num_primitive_torsions = ni-num_primitive_bonds-num_primitive_angles
    intcor0 = ni

    if (.not. first_pass) then
       if (di_iprint>=4 .and. di_on_root) then
          write(di_stdout,'(/,a,i5)') 'Maximum number of primitive dihedrals:    ',num_primitive_torsions+ibad
          if (ibad/=0) write(di_stdout,2502) ibad 
       endif

       ! final printout
       if (di_on_root) then
          write(di_stdout,'(/,a,i7)') 'Total number of primitive bonds:        ',num_primitive_bonds
          write(di_stdout,'(a,i7)')   'Total number of primitive angles:       ',num_primitive_angles
          write(di_stdout,'(a,i7,/)') 'Total number of primitive dihedrals:    ',num_primitive_torsions
          write(di_stdout,'(a,i7,/)') 'Total number of primitive internals:    ',ni
       endif
    endif ! .not. first_pass

 2500 format(5x,'Three or more atoms for torsion are virtually linear')
 2502 format(i5,'  Dihedral angles with linear bends were skipped')
    return
  end subroutine di_algor_make_list_primitives

!------------------------------------------------------------------------------

  subroutine di_algor_make_F_primitives(nat3,nic,ired,h,hh,eigv,b,indexb,lidelp,xprim,ktyp,icc,ierr)
    !=========================================================================!
    ! Make F matrix for bonds, angles, torsions  [F=B(T)*B]                   !
    !                                                                         !
    !                        [domig in DMol]                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) nat3 (in) : dimension (3*number_atoms)                               !
    ! 2) nic (in) : dimension (3*number_atoms)                                !
    ! 3) ired(out) : number of redundant coordinates                          !
    ! 4) h(out) : F matrix  eigenvectors                                      !
    ! 5) hh(out) : F matrix  itself                                           !
    ! 6) eigv(out) : F matrix eigenvalues                                     !
    ! 7) b(inout) : B matrix                                                  !
    ! 8) indexb(inout) : index of B matrix                                    !
    ! 9) lidelp(inout) : atom IDs and lattice translations for internals      !
    !10) xprim(inout) : internal coordinates                                  !
    !11) ktyp(inout) : types of internals                                     !
    !12) icc(inout) : list of constraints                                     !
    !13) ierr(out)                                                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) ncons_2                                                              !
    ! 5) ncons                                                                !
    ! 6) num_primitive_bonds                                                  !
    ! 7) num_primitive_angles                                                 !
    ! 8) num_primitive_torsions                                               !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) num_primitive_angles                                                 !
    ! 2) num_primitive_torsions                                               !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none

    integer,intent(in)                                    :: nat3
    integer,intent(in)                                    :: nic
    integer,intent(out)                                   :: ired
    real(kind=di_dp),dimension(1:nat3,1:nat3),intent(out) :: h
    real(kind=di_dp),dimension(1:nat3,1:nat3),intent(out) :: hh
    real(kind=di_dp),dimension(1:nat3),intent(out)        :: eigv
    real(kind=di_dp),dimension(1:12,1:nic),intent(inout)  :: b
    integer,dimension(1:12,1:nic),intent(inout)           :: indexb
    integer,dimension(:,:),intent(inout)                  :: lidelp
    real(kind=di_dp),dimension(:),intent(inout)           :: xprim
    integer,dimension(1:intcor),intent(inout)             :: ktyp
    integer,dimension(1:ncons_2,1:21),intent(inout)       :: icc
    integer,intent(out)                                   :: ierr
    !-------------------------------------------------------------------------!
    integer :: nunb0,nuna0,ni,i,j,k,ii,jj,intnew,intold,nundn
    !-------------------------------------------------------------------------!
    ierr = 0

    nunb0 = num_primitive_bonds
    nuna0 = num_primitive_angles

      
    ! form F matrix for 1) bonds, 
    !                   2) bonds + "non zero" angles
    !                   3) bonds + "non zero" angles + "non zero" torsions
    ! decide if the primitive is "non zero" by projecting its b vector
    ! from zero vectors
    !
    ! make sure that every constraint from icc list ends up as active primitive

    ni = num_primitive_bonds
    ! bonds
    h = zero

    do i = 1,6
       do j = 1,6
          do k = 1,ni
             ii = indexb(i,k)
             jj = indexb(j,k)
             if (ii*jj/=0) h(ii,jj) = h(ii,jj) + b(i,k)*b(j,k)
          enddo
       enddo
    enddo
    ! save h
    hh = h

    if (di_iprint>=2 .and. di_on_root) write(di_stdout,'(/,a,3i7)') 'Diagonalize F:  bonds',nat3,ni
    call di_algor_diagonalize_F_matrix(nat3,h,eigv,ired,ierr)
    if (ierr/=0)  then
       if (di_on_root) write(di_stdout,*) ' problem in di_algor_make_F_primitives (bonds)'
       go to 999
    endif

    ! second pass
    ! scan angles , find overlap with zero H and add to the B matrix
    intold = num_primitive_bonds
    num_primitive_angles=0
    num_primitive_torsions=0
    call di_algor_make_F_angles(nat3,ired,h,b,indexb,ktyp,xprim, &
     &           lidelp,num_primitive_bonds,nunb0+nuna0,intnew,nundn,icc)
    ni = intnew
    num_primitive_angles = nundn

    if (intnew/=intold) then
       ! bends (angles) 
       h = zero

       do i = 1,9
          do j = 1,9
             do k = 1,ni
                ii = indexb(i,k)
                jj = indexb(j,k)
                if (ii*jj/=0) h(ii,jj) = h(ii,jj) + b(i,k)*b(j,k)
             enddo
          enddo
       enddo
       ! save h
       hh = h

       if (di_iprint>=2 .and. di_on_root) write(di_stdout,'(/,a,3i7)') 'Diagonalize F:  bends',nat3,num_primitive_angles
       call di_algor_diagonalize_F_matrix(nat3,h,eigv,ired,ierr )
       if (ierr/=0)  then
          if (di_on_root) write(di_stdout,*) ' problem in di_algor_make_F_primitives (bends)'
          go to 999
       endif
    endif ! intnew/=intold

    ! third pass
    ! scan dihedrals (torsions), find overlap with zero H and add to the B matrix
    intold = ni   

    call di_algor_make_F_torsions(nat3,ired,h,b,indexb,ktyp,xprim, &
     &            lidelp,num_primitive_bonds+num_primitive_angles,nunb0+nuna0, &
     &            intnew,nundn,icc)
    ni = intnew
    num_primitive_torsions = nundn

    if (intnew/=intold) then
       ! torsions
       h = zero


       do i = 1,12
          do j = 1,12
             do k = 1,ni
                ii = indexb(i,k)
                jj = indexb(j,k)
                if (ii*jj/=0) h(ii,jj) = h(ii,jj) + b(i,k)*b(j,k)
             enddo
          enddo
       enddo
       ! save h
       hh = h
          
       if (di_iprint>=2 .and. di_on_root) write(di_stdout,'(/,a,3i7)') 'Diagonalize F:  torsions',nat3,num_primitive_torsions
       call di_algor_diagonalize_F_matrix(nat3,h,eigv,ired,ierr )
       if (ierr/=0)  then
          if (di_on_root) write(di_stdout,*) ' problem in di_algor_make_F_primitives (torsions)'
          go to 999
       endif
    endif ! intnew/=intold

    ! redefine intcor
    if (intcor/=intnew) intcor = intnew
 
    ! redefine icc if constraints
    do j=1,ncons
       icc(j,18) = icc(j,19)
       if (di_iprint>=4 .and. di_on_root) write(di_stdout,*) 'icc is redefined ',j,icc(j,18)
    enddo

 999 continue 

    return
  end subroutine di_algor_make_F_primitives

!------------------------------------------------------

  subroutine di_algor_make_F_angles(nat3,ired,h,b,indexb,ktyp,xprim, &
     &             lidelp,nunba,intcor0,intnew,nundn,icc)
    !=========================================================================!
    ! Prepare B matrix for bond angles                                        !
    !                                                                         !
    !                        [findang in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) nat3 (in) : dimension (3*number_atoms)                               !
    ! 2) ired(in) : number of redundant coordinates                           !
    ! 3) h(in) : F matrix  eigenvectors                                       !
    ! 4) b(inout) : B matrix                                                  !
    ! 5) indexb(in) : index of B matrix                                       !
    ! 6) ktyp(inout) : types of internals                                     !
    ! 7) xprim(inout) : internal coordinates                                  !
    ! 8) lidelp(inout) : atom IDs and lattice translations for internals      !
    ! 9) nunba(in) : start point in the list of internals                     !
    !10) intcor0(in) : current number of internals                            !
    !11) intnew(out) : new number of internals                                !
    !12) nundn(out) : number of newly found internals                         !
    !13) icc(inout) : list of constraints                                     !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) ncons_2                                                              !
    ! 5) ncons                                                                !
    ! 6) tp_prprim                                                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                   :: nat3
    integer,intent(in)                                   :: ired
    real(kind=di_dp),dimension(1:nat3,1:nat3),intent(in) :: h
    real(kind=di_dp),dimension(:,:),intent(inout)        :: b
    integer,dimension(:,:),intent(inout)                 :: indexb
    integer,dimension(:),intent(inout)                   :: ktyp
    real(kind=di_dp),dimension(:),intent(inout)          :: xprim
    integer,dimension(:,:),intent(inout)                 :: lidelp
    integer,intent(in)                                   :: nunba
    integer,intent(in)                                   :: intcor0
    integer,intent(out)                                  :: intnew
    integer,intent(out)                                  :: nundn
    integer,dimension(1:ncons_2,1:21),intent(inout)      :: icc
    !-------------------------------------------------------------------------!
    real(kind=di_dp),dimension(:),allocatable       :: v
    real(kind=di_dp) :: thrsh,ss,ss1,sw
    integer :: ni,in,i,j,k,iz,ija,ijb,ijc,jb,i0,j0,k0
    integer :: ika,ikb,ikc,jc,ia,ierr
    logical :: skip
    !-------------------------------------------------------------------------!

    thrsh = tp_prprim
    allocate(v(1:nat3),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating v in di_algor_make_F_angles')
    v = zero

    ! loop over all b vectors which are bends  and
    ! find if they have overlap with the first ired of zeros

    ! (it is NOT assumed that translations are still in)

    if (di_iprint>4 .and. di_on_root) then
       write(di_stdout,*) '===di_algor_make_F_angles ===',nunba,intcor0,ired,thrsh
    endif

    ni = nunba     
      
    do in = nunba+1, intcor0

       skip =.true.
        
       ! go over zeros
       ss1=zero

       ! unroll b into v
       do k=1,nat3
          v(k) = zero
       enddo
       do i = 1,9
          v(indexb(i,in)) = v(indexb(i,in)) + b(i,in)
       enddo

       ! go over zeros
       do iz = 1,ired

          ss = zero
          do k=1,nat3
             ss =  ss + v(k)*h(k,iz)    
          enddo

          ! if not zero, check if this overlap with the rest
          if (abs(ss)>thrsh) skip=.false.
          if (abs(ss)>ss1) then
             ss1= abs(ss)
          endif
       enddo   ! iz


       ! make sure constraints are always included
       ! but their corresponding primitive number may be redefined
       do j=1,ncons
          if (icc(j,17)==3 .and. icc(j,18)==in) then
             skip =.false.
             icc(j,19) = ni+1   
          endif
       enddo

       if (.not.skip) then

          if (di_iprint>4 .and. di_on_root) then

             call deloc_utils_unpack_atom_indices(lidelp(1,in), ija, ijb, ijc, jb)
             call deloc_utils_unpack_atom_indices(lidelp(2,in),  i0,  j0,  k0, ia)
             call deloc_utils_unpack_atom_indices(lidelp(3,in), ika, ikb, ikc, jc)

             sw = xprim(in)*toang
             write(di_stdout,1122) jb,'(',ija,ijb,ijc,')', &
     &                ia,'(',i0,i0,i0,')', &
     &                jc,'(',ika,ikb,ikc,')',sw,'=',ss1
          endif

          ! save lidelp
          ni=ni+1
          do i=1,3
             lidelp(i,ni) = lidelp(i,in)
          enddo
          lidelp(4,ni)=0

          ! save xprim
          xprim(ni) = xprim(in)

          ! save ktyp
          ktyp(ni) = ktyp(in)

          ! save b
          do k = 1,9
             b(k,ni) = b(k,in)
             indexb(k,ni) = indexb(k,in)
          enddo
       endif  ! .not.skip

    enddo
    intnew = ni
    nundn =  ni - nunba

    deallocate(v,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating v in di_algor_make_F_angles')

 1122   format (3(i3,a1,3i2,a1),f8.2,a1,e10.2)

    return
  end subroutine di_algor_make_F_angles



!------------------------------------------------------

  subroutine  di_algor_make_F_torsions(nat3,ired,h,b,indexb,ktyp,xprim, &
     &             lidelp,nunba,nunba0,intnew,nundn,icc)
    !=========================================================================!
    ! Prepare B matrix for torsions (dihedrals)                               !
    !                                                                         !
    !                        [finddih in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) nat3 (in) : dimension (3*number_atoms)                               !
    ! 2) ired(in) : number of redundant coordinates                           !
    ! 3) h(in) : F matrix  eigenvectors                                       !
    ! 4) b(inout) : B matrix                                                  !
    ! 5) indexb(in) : index of B matrix                                       !
    ! 6) ktyp(inout) : types of internals                                     !
    ! 7) xprim(inout) : internal coordinates                                  !
    ! 8) lidelp(inout) : atom IDs and lattice translations for internals      !
    ! 9) nunba(in) : start point in the list of internals                     !
    !10) nunba0(in) : start point in the list of internals                    !
    !11) intnew(out) : new number of internals                                !
    !12) nundn(out) : number of newly found internals                         !
    !13) icc(inout) : list of constraints                                     !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) ncons_2                                                              !
    ! 5) ncons                                                                !
    ! 6) tp_prprim                                                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                   :: nat3
    integer,intent(in)                                   :: ired
    real(kind=di_dp),dimension(1:nat3,1:nat3),intent(in) :: h
    real(kind=di_dp),dimension(:,:),intent(inout)        :: b
    integer,dimension(:,:),intent(inout)                 :: indexb
    integer,dimension(1:intcor),intent(inout)            :: ktyp
    real(kind=di_dp),dimension(:),intent(inout)          :: xprim
    integer,dimension(:,:),intent(inout)                 :: lidelp
    integer,intent(in)                                   :: nunba
    integer,intent(in)                                   :: nunba0
    integer,intent(out)                                  :: intnew
    integer,intent(out)                                  :: nundn
    integer,dimension(1:ncons_2,1:21),intent(inout)      :: icc
    !-------------------------------------------------------------------------!
    real(kind=di_dp),dimension(:),allocatable       :: v
    logical skip
    real(kind=di_dp) :: thrsh,ss1,ss,sw
    integer :: ni,in,i,k,iz,ika,ikb,ikc,kc,j,i0,j0,k0,itra,itrb,itrc,jb
    integer :: ila,ilb,ilc,ld,ia,ierr
    !-------------------------------------------------------------------------!

    thrsh = tp_prprim

    allocate(v(1:nat3),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating v in di_algor_make_F_torsions')
    v = zero

    ! loop over all B-matrix vectors which are dihedrals and
    ! find if they have overlap with the first ired of zeros

    ! (it is NOT assumed that translations are still in)
      
    if (di_iprint>4 .and. di_on_root) then
       write(di_stdout,*) '===di_algor_make_F_torsions ===',nunba,nunba0,intcor,ired,thrsh
    endif

    ni = nunba     
      
    do in = nunba0+1, intcor

       ! unroll B into V
       do k=1,nat3
          v(k) = zero
       enddo
       do i = 1,12
          v(indexb(i,in)) = v(indexb(i,in)) + b(i,in)
       enddo

       ss1 = zero
       skip =.true. 
       ! go over zeros
       do iz = 1,ired

          ss = zero
          do k=1,nat3
             ss=ss+v(k)*h(k,iz)    
          enddo

          ! if not zero, check if this overlap with the rest
          if (abs(ss)>thrsh) skip=.false.
          if (abs(ss)>ss1) ss1 = abs(ss)
       enddo

       ! make sure constraints are always included
       ! but their corresponding primitive number may be redefined
       do j=1,ncons
          if (icc(j,17)==5 .and. icc(j,18)==in) then
             skip =.false.
             icc(j,19) = ni+1   
          endif
       enddo

       if (.not.skip) then

          if (di_iprint>4 .and. di_on_root) then

             call deloc_utils_unpack_atom_indices(lidelp(1,in), ika, ikb, ikc, kc)
             call deloc_utils_unpack_atom_indices(lidelp(2,in),  i0,  j0,  k0, ia)
             call deloc_utils_unpack_atom_indices(lidelp(3,in),itra,itrb,itrc, jb)
             call deloc_utils_unpack_atom_indices(lidelp(4,in), ila, ilb, ilc, ld)

             sw = xprim(in)*toang
             write(di_stdout,1122) kc,'(',ika,ikb,ikc,')', &
     &                ia,'(',i0,i0,i0,')', &
     &                jb,'(',itra,itrb,itrc,')', &
     &                ld,'(',ila, ilb, ilc,')',sw,'=',ss1
          endif

          ni = ni+1

          if (ni/=in) then

             ! save lidelp
             do i=1,4
                lidelp(i,ni) = lidelp(i,in)
             enddo

             ! save savtor
             savtor(ni) = savtor(in)

             ! save xprim
             xprim(ni) = xprim(in)

             ! save ktyp
             ktyp(ni) = ktyp(in)

             ! save b
             do k = 1,12
                b(k,ni) = b(k,in)
                indexb(k,ni) = indexb(k,in)
             enddo
          endif ! ni/=in

       endif  ! .not.skip

    enddo

    intnew = ni
    nundn =  ni - nunba

    deallocate(v,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating v in di_algor_make_F_torsions')

 1122   format (4(i3,a1,3i2,a1),f8.2,a1,e10.2)

    return
  end subroutine di_algor_make_F_torsions


!------------------------------------------------------

  subroutine di_algor_make_F_disc_fragm(nat3,nic,intcor,ired,h,hh,eigv,lidelp,xint,ktyp,status)
    !=========================================================================!
    ! Make F matrix for disconnected fragments                                !
    !                                                                         !
    !                        [sepfrag in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) nat3 (in) : dimension (3*number_atoms)                               !
    ! 2) nic(in) : max dimension of DI space                                  !
    ! 3) intcor(inout) : actual number of DI variables                        !
    ! 4) ired(out) : number of redundant coordinates                          !
    ! 5) h(out) : F matrix  eigenvectors                                      !
    ! 6) hh(inout) : F matrix itself                                          !
    ! 7) eigv(out) : eigenvalues of F matrix                                  !
    ! 8) lidelp(out) : list of atoms and cell translations for DI             !
    ! 9) ktyp(out) : types of internals                                       !
    !10) xint(out) : internal coordinates                                     !
    !11) ierr(out)                                                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) tp_eigvz                                                             !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) bmat_disconn_fragm                                                   !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                      :: nat3
    integer,intent(in)                                      :: nic
    integer,intent(inout)                                   :: intcor
    integer,intent(inout)                                   :: ired
    real(kind=di_dp),dimension(1:nat3,1:nat3),intent(inout) :: h 
    real(kind=di_dp),dimension(1:nat3,1:nat3),intent(inout) :: hh
    real(kind=di_dp),dimension(1:nat3),intent(out)          :: eigv
    integer,dimension(:,:),intent(inout)                    :: lidelp
    integer,dimension(:),intent(inout)                      :: ktyp
    real(kind=di_dp),dimension(:),intent(inout)             :: xint
    integer,intent(out)                                     :: status
    !-------------------------------------------------------------------------!
    real(kind=di_dp),external :: ddot
    real(kind=di_dp) :: thrsh,snorm,s1,cf,ss
    integer :: nel,i,j,k,kk,ii,nsr,intnew,ierr
    real(kind=di_dp),dimension(:,:),allocatable  :: vm
    real(kind=di_dp),dimension(:),allocatable    :: v
    !-------------------------------------------------------------------------!
    status = 0

    allocate(vm(1:nat3,1:nat3),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating vm in di_algor_make_F_disc_fragm')
    vm = zero
    allocate(v(1:nat3),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating v in di_algor_make_F_disc_fragm')
    v = zero

    thrsh = tp_schmidt_a

    ! check for translations 
    ! project the x, y, z vectors on the zeros


    ! eliminate 3 translations (for solid)
    nel = 3

    ! translations
    do i=1,3 
       do j=1,nat3 
          do k=1,nat3
             kk= mod(k,3)
             v(k) = zero
             ! x
             if (i==1 .and. kk==1) v(k) = h(k,j)
             ! y
             if (i==2 .and. kk==2) v(k) = h(k,j)
             ! z
             if (i==3 .and. kk==0) v(k) = h(k,j)
          enddo
          do k=1,nat3
             vm(k,i) = vm(k,i) + v(k)*h(k,j)
          enddo
       enddo
    enddo

    ! normalize translations
    do i=1,nel  
       snorm = ddot(nat3,vm(1,i),1,vm(1,i),1)
       if (snorm<thrsh) then
          if (di_on_root) write (di_stdout,*) 'problem in sch_tr',snorm,i
       else
          snorm = one/sqrt(snorm)
          call dscal(nat3,snorm,vm(1,i),1)
       endif
    enddo

    ! Schmidt-Orthogonalize the translations 
    ! and rotations if nel=6

    ! do the zeros

    ii=nel+1

    do i=1,ired

       vm(:,ii) = h(:,i)

       do j=ii-1,1,-1

          cf = ddot(nat3,vm(1,j),1,vm(1,ii),1)
          do k=1,nat3   
             vm(k,ii) = vm(k,ii) - cf*vm(k,j)
          enddo
       enddo

       snorm = ddot(nat3,vm(1,ii),1,vm(1,ii),1)

       if (di_iprint>4 .and. di_on_root) write(di_stdout,*) 'snorm sep ',i,ii,ired,snorm

       if (snorm>=thrsh .and. ii<=ired) then
          snorm = one/sqrt(snorm)
          call dscal(nat3,snorm,vm(1,ii),1)
          ii = ii+1
       endif

    enddo

    ! number of active rotations is ii-4
    nsr = ii - nel - 1

    ! redefine intcor
    intnew = intcor
    intcor = intcor + nsr 
    num_disconn_fragm = nsr
    num_disconn_fragm_2 = max(num_disconn_fragm,1)

    if (allocated(bmat_disconn_fragm)) then
       if (size(bmat_disconn_fragm,2)/=num_disconn_fragm_2) then
          deallocate(bmat_disconn_fragm,stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bmat_disconn_fragm in di_algor_get_init_deloc_coords')
          allocate(bmat_disconn_fragm(1:nat3,1:num_disconn_fragm_2),stat=ierr)
          if (ierr/=0) call deloc_utils_io_abort('Error in allocating bmat_disconn_fragm in di_algor_get_init_deloc_coords')
          bmat_disconn_fragm = 0.0_di_dp
       endif
    else
       allocate(bmat_disconn_fragm(1:nat3,1:num_disconn_fragm_2),stat=ierr)
       if (ierr/=0) call deloc_utils_io_abort('Error in allocating bmat_disconn_fragm in di_algor_get_init_deloc_coords')
       bmat_disconn_fragm = 0.0_di_dp
    endif

    if (intcor>nic) then
       if (di_on_root) write(di_stdout,*) 'ni + nfrc is too large ',intcor,nic
       status = 1
       go to 999
    endif


    ! fill ktyp with 6  ,lidelp with 0
    do k=1,nsr  
       ktyp(k+intnew) = 6
       call deloc_utils_pack_atom_indices(lidelp(1,k+intnew),   0,   0,   0,  0 )
       call deloc_utils_pack_atom_indices(lidelp(2,k+intnew),   0,   0,   0,  0 )
       call deloc_utils_pack_atom_indices(lidelp(3,k+intnew),   0,   0,   0,  0 )
       call deloc_utils_pack_atom_indices(lidelp(4,k+intnew),   0,   0,   0,  0 )
       ! make xint
       ii=0
       s1 = zero
       do j=1,number_atoms
          do i=1,nel
             ii=ii+1
             s1 = s1 + coords_cart(i,j)*vm(ii,k+3)
          enddo
       enddo
       xint(k+intnew) = s1

    enddo

    ! expand b
    do k=1,nsr
       do i=1,nat3
          bmat_disconn_fragm(i,k) =  vm(i,k+nel)
       enddo
    enddo


    if (di_on_root) then
       write(di_stdout,'(/,a,/)') 'Found disconnected fragments'
       if (di_iprint>1) write(di_stdout,'(a,i5)') ' Number of primitive internals increased by ',num_disconn_fragm

       if (di_iprint>4) write(di_stdout,*) 'intcor,nic, nsr,nat3 ',intcor,nic, nsr,nat3
    endif

    ! make a new F and diagonalize
    ! restore h
    h = hh

    do i=1,nat3
       do j=1,i
          ss=zero
          do k=1,nsr
             ss=ss + bmat_disconn_fragm(i,k)*bmat_disconn_fragm(j,k)
          enddo
          h(i,j) = h(i,j) + ss
          h(j,i) = h(i,j)
       enddo
    enddo

    ! save h
    hh = h

    if (di_iprint>1 .and. di_on_root) write(di_stdout,'(/,a,2i5)') &
    &    'Diagonalize F: disconnected fragments ',nat3,num_disconn_fragm

    call di_algor_diagonalize_F_matrix(nat3,h,eigv,ired,status)

    if (status/=0)  then
       if (di_on_root) write(di_stdout,*) ' problem in di_algor_make_F_disc_fragm'
       go to 999
    endif

 999 continue

    deallocate(vm,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating vm in di_algor_make_F_disc_fragm')
    deallocate(v,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating vm in di_algor_make_F_disc_fragm')

    return

  end subroutine di_algor_make_F_disc_fragm


!------------------------------------------------------

  subroutine di_algor_diagonalize_F_matrix(n,h,e,ired,status)
    !=========================================================================!
    ! Make F matrix = B(T) * B                                                !
    !                                                                         !
    !                        [makef_d in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) n (in) : dimension (3*number_atoms)                                  !
    ! 2) h(inout) : F matrix (contains eigenvecotrs on exit)                  !
    ! 3) e(out) : ordered eigenvalues of the F matrix                         !
    ! 4) ired(out) : number of redundant coordinates                          !
    ! 5) ierr(out)                                                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) tp_eigvz                                                             !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                 :: n
    real(kind=di_dp),dimension(1:n,1:n),intent(inout)  :: h
    real(kind=di_dp),dimension(1:n),intent(out)        :: e
    integer,intent(out)                                :: ired
    integer,intent(out)                                :: status
    !-------------------------------------------------------------------------!
    integer :: info,i,nwrk,ierr
    !integer :: lwrk
    real(kind=di_dp),dimension(:),allocatable :: v
    !integer,dimension(:),allocatable          :: lv
    !-------------------------------------------------------------------------!
    status = 0

    ! Allocate workspace (as per LAPACK 2.0)
    !i = 1
    !do while (2**i < n)
    !   i = i + 1
    !enddo
    ! dimension for DSYEVD in LAPACK 2
    !nwrk = 1 + 5*n + 2*n*i + 3*n*n
    !lwrk = 3 + 5*n

    ! nwrk in LAPACK 3.0 is smaller:
    !  nwrk = 1 + 6*n + 2*n*n

    ! nwrk for DSYEV
    nwrk = n*(n+2)
    allocate(v(1:nwrk),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating v in di_algor_diagonalize_F_matrix')
    v = zero
    !allocate(lv(1:lwrk),stat=ierr)
    !if (ierr/=0) call deloc_utils_io_abort('Error in allocating lv in di_algor_diagonalize_F_matrix')
    !lv = 0

    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,'(/,6x,a)')'B(t)*B Matrix'
       call di_algor_print_matrix(n,n,n,h)
       write(di_stdout,'(a)')'END of B(t)*B Matrix'
    endif

    ! find eigenvectors and eigenvalues of the F matrix
!    call dsyevd('V','U',n,h,n,e,v,nwrk,lv,lwrk,info)
    
!    if (info/=0) then
!       if (di_on_root) write(di_stdout,*) 'Error in dsyevd; try old dsyev'

    ! removed DSYEVD call since apparently it is not supported on CRAY
    call dsyev('V','U',n,h,n,e,v,nwrk,info)
!    endif

    ! clear up memory
    deallocate(v,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating v in di_algor_diagonalize_F_matrix')
    !deallocate(lv,stat=ierr)
    !if (ierr/=0) call deloc_utils_io_abort('Error in deallocating lv in di_algor_diagonalize_F_matrix')

    if (info/=0) then
       write(di_stdout,1200)
       go to 999
    endif

    ! Print out eigenvalues and eigenvectors
    if (di_iprint>1 .and. di_on_root) then
       write(di_stdout,'(/,3x,a)')'Eigenvalues of B(t)*B:'
       write(di_stdout,'(1x,6f12.6)') (e(i),i=1,n)
       write(di_stdout,'(/)')
    endif

    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,'(/,3x,a)')'Eigenvectors of B(t)*B:'
       call di_algor_print_matrix(n,n,n,h)
       write(di_stdout,'(a)') 'END of Eigenvectors of B(t)*B'
    endif

    !  Remove redundancies from coordinate set and check that
    !  we have a full set of 3*NATOMS-3 internal coordinates

    do i=1,n   
       if (abs(e(i))>tp_eigvz) go to 21
    end do

    !  should never get here

    status = 1
    if (di_on_root) write(di_stdout,1600)
    go to 999

 21 continue

    ired = i-1

999 continue
    return

 1200 format(/,2x,' Warning: Unable to Diagonalize B(t)*B in di_algor_diagonalize_F_matrix')
 1600 format(/,2x,' Warning: B(t)*B Matrix has all zero eigenvalues in di_algor_diagonalize_F_matrix!')

  end subroutine di_algor_diagonalize_F_matrix

!---------------------------------------------------------------

  subroutine di_algor_get_step(ndeg,mode,ndiis,neg,negreq,status)
    !=========================================================================!
    ! Find the step for geometry optimization in delocalized internals:       !
    ! the version with unsatisfied internal constraints                       !
    !                                                                         !
    ! GDIIS route is currently commented out                                  !
    !                                                                         !
    !  References                                                             !
    !  ----------                                                             !
    !                                                                         !
    !  "An Algorithm for the Location of Transition States"                   !
    !   J.Baker  J.Comp.Chem.  7 (1986) 385                                   !
    !                                                                         !
    !  "Geometry Optimization by Direct Inversion in the Iterative Subspace"  !
    !   P.Csaszar and P.Pulay  J.Mol.Struct.  114 (1984) 31                   !
    !                                                                         !
    !  "Geometry Optimization in Cartesian Coordinates:                       !
    !   The End of the Z-Matrix?"                                             !
    !   J.Baker and W.J.Hehre  J.Comp.Chem.  12 (1991) 606                    !
    !                                                                         !
    !                        [get_step in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndeq (in) : number of internal degrees of freedom (ns+nunsat)        !
    ! 2) mode(inout) : Hessian mode being followed during a TS search         !
    !               if mode=0 mode following is switched off and              !
    !               maximization will occur along the lowest mode             !
    ! 3) ndiis(inout) : current size of GDIIS subspace                        !
    ! 4) neg(out) : number of negative eigenvalues                            !
    ! 5) negreq(out) : required number of negative eigenvalues                !
    ! 6) status(out)                                                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) maxdiis                                                              !
    ! 5) disp_DI                                                              !
    ! 6) coords_DI                                                            !
    ! 7) gradient_DI                                                          !
    ! 8) gradient_DI_old                                                      !
    ! 9) tsflag                                                               !
    !10) tp_bad_hes                                                           !
    !11) hess_DI                                                              !
    !12) ihess     : Hessian status flag                                      !
    !               0 - have "exact" or initial Hessian   use as is           !
    !               1 - have Hessian from previous step   need to update      !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) hess_DI                                                              !
    ! 2) disp_DI                                                              !
    ! 3) coords_DI                                                            !
    ! 4) gradient_DI                                                          !
    ! 5) gradient_DI_old                                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                        :: ndeg
    integer,intent(inout)                                     :: mode
    integer,intent(inout)                                     :: ndiis
    integer,intent(out)                                       :: neg
    integer,intent(out)                                       :: negreq
    integer,intent(out)                                       :: status
    !-------------------------------------------------------------------------!
    real(kind=di_dp),dimension(:,:),allocatable       :: u
    real(kind=di_dp),dimension(:),allocatable         :: eigv
    real(kind=di_dp),dimension(:),allocatable         :: work1
    real(kind=di_dp) :: eigmin,eigmax,scale,toldiis,rmsg,rmsd
    real(kind=di_dp), external :: ddot
    integer :: iunsat,i,j,nc,ierr
    real(kind=di_dp), parameter :: small=1.0e-8_di_dp
    logical :: diisflag
    !-------------------------------------------------------------------------!

    status = 0
    if (di_iprint>4 .and. di_on_root) then
       write(di_stdout,*) 'di_algor_get_step:'
       write(di_stdout,*) 'ndeg, ncycle, di_iprint, ihess ',ndeg, ncycle, di_iprint, ihess
    endif

    !  Check options and set up defaults

    eigmin = 0.0001_di_dp
    eigmax = 25.0_di_dp
    scale = one
    toldiis = 0.1_di_dp

    iunsat = 0
    negreq = iunsat
    if (tsflag) negreq = iunsat + 1

    allocate(work1(1:(ndeg+2)*ndeg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating work1 in di_algor_get_step')
    allocate(u(1:ndeg,1:ndeg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating u in di_algor_get_step')
    allocate(eigv(1:ndeg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating eigv in di_algor_get_step')

    !  START THE OPTIMIZATION PROPER
    !  First update the Hessian if necessary

    if (ihess==1) then
       call di_algor_update_hessian(ndeg-iunsat,disp_DI,gradient_DI,gradient_DI_old,hess_DI)
    endif

 10 continue

    do i = 1,ndeg
      do j = 1,ndeg
        u(j,i) = hess_DI(j,i)
      end do
    end do
    
    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,1000)
       call di_algor_print_matrix(ndeg,ndeg,ndeg,u)
    endif

    !  Diagonalize

    call dsyev('V','L',ndeg,u,ndeg,eigv,work1,(ndeg+2)*ndeg,status)
    if (status/=0) then
       if (di_on_root) write(di_stdout,1010)
       go to 999
    endif

    ! if Eigenvalues are very bad; = very small;
    ! restore the old Hessian, and make a new set of delocalized internals
    if (eigv(1)<tp_bad_hes .and. (.not.TSflag)) then
       status = 2
       ! print eigenvalues
       if (di_on_root) then
          write(di_stdout,1120)
          write(di_stdout,1130) (eigv(i),i=1,ndeg)
       endif
       go to 999
    endif

    !  check for and remove eigenvectors corresponding to zero
    !  eigenvalues and possible symmetry-redundant modes

    call di_algor_cleanup_hessian(ndeg,iunsat,u,eigv,gradient_DI,small,nc,status)
    if (status/=0) go to 999

    !  on exit, NC is the number of Hessian modes that will be
    !  used to form the next step.

    !  check eigenvalue magnitudes are between allowed values

    call di_algor_adjust_hess_EVs_range(nc,eigmin,eigmax,eigv)

    !  find the number of negative Hessian eigenvalues

    call di_algor_num_negative_hess_EVs(nc,eigv,neg)

    !  form the RMS gradient

    rmsg = ddot(ndeg,gradient_DI(1),1,gradient_DI(1),1)
    rmsg = sqrt(rmsg/float(nc))

    !  calculate the next step
    !  this can be done using either GDIIS or the
    !  standard EF algorithm

    !  see if GDIIS can be used (minimization only)

    diisflag = maxdiis>1 .and. rmsg<toldiis .and. neg==0 .and. .not.tsflag

    if (diisflag) then
       
!       call gdiis(nc,ndeg,maxdiis,ndiis,u,eigv,coords_DI,gradient_cart,xstr,gstr,d)
       if (ndiis<=1) diisflag = .false.

    else

       call di_algor_eigenvector_following(nc,ndeg,iunsat,neg,negreq,mode,u,eigv,gradient_DI,disp_DI,status)

    endif


    if (status==-1) then
       call di_algor_ef_check_H_scale(ndeg-iunsat,eigv,scale)
       if (scale==one) then
          go to 95
       else
          go to 10
       endif
    endif
    if (status>0) go to 999

    !  we have a new step in D
    !  check the stepsize

    call di_algor_check_vector_length(ndeg,disp_DI)

    !    Standard Printout

    if (di_iprint>2) then

       !  form the RMS displacement

       rmsd = ddot(ndeg,disp_DI(1),1,disp_DI(1),1)
       rmsd = sqrt(rmsd/float(nc))

       if (di_on_root) write(di_stdout,2000) ncycle,ec,rmsg,rmsd
    endif


    gradient_DI_old = gradient_DI

    coords_DI = coords_DI + disp_DI

    go to 999

 95 continue
    status = 3

 999 continue

    deallocate(work1,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating work1 in di_algor_get_step')
    deallocate(u,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating u in di_algor_get_step')
    deallocate(eigv,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating eigv in di_algor_get_step')

    return

 1000 format(/,'   Hessian Matrix')
 1010 format(/,2X,' Warning: ',/,' Unable to Diagonalize Hessian Matrix')
 2000 format(/,' ** Cycle ',I3,'  Energy ',F16.9,'   RMSG ',F8.5, &
     &         '   RMSD ',F8.5,' **')
 1120 format('  Hessian Eigenvalues:')
 1130 format(1X,6F12.6)

  end subroutine di_algor_get_step

!------------------------------------------------------------------------------

  subroutine di_algor_get_step_unsat_constr(ndeg,mode,neg,negreq,hpad,status)
    !=========================================================================!
    ! Find the step for geometry optimization in delocalized internals:       !
    ! the version with unsatisfied internal constraints                       !
    !                                                                         !
    !  References                                                             !
    !  ----------                                                             !
    !                                                                         ! 
    !  "An Algorithm for the Location of Transition States"                   !
    !   J.Baker  J.Comp.Chem.  7 (1986) 385                                   !
    !                                                                         !
    !  "Geometry Optimization by Direct Inversion in the Iterative Subspace"  !
    !   P.Csaszar and P.Pulay  J.Mol.Struct.  114 (1984) 31                   !
    !                                                                         !
    !  "Geometry Optimization in Cartesian Coordinates:                       !
    !   The End of the Z-Matrix?"                                             !
    !   J.Baker and W.J.Hehre  J.Comp.Chem.  12 (1991) 606                    !
    !                                                                         !
    !                        [get_step_u in DMol]                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndeq (in) : number of internal degrees of freedom (ns+nunsat)        !
    ! 2) mode(inout) : Hessian mode being followed during a TS search         !
    !               if mode=0 mode following is switched off and              !
    !               maximization will occur along the lowest mode             !
    ! 3) neg(out) : number of negative eigenvalues                            !
    ! 4) negreq(out) : required number of negative eigenvalues                !
    ! 5) hpad(out) : padded strip of hamiltonian                              !
    ! 6) status(out)                                                          !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) nunsat                                                               !
    ! 5) disp_DI                                                              !
    ! 6) coords_DI                                                            !
    ! 7) gradient_DI                                                          !
    ! 8) gradient_DI_old                                                      !
    ! 9) tsflag                                                               !
    !10) tp_bad_hes                                                           !
    !11) hess_DI                                                              !
    !12) ihess :     Hessian status flag                                      !
    !               0 - have "exact" or initial Hessian   use as is           !
    !               1 - have Hessian from previous step   need to update      !
    !13) icc                                                                  !
    !14) rcon                                                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) hess_DI                                                              !
    ! 2) disp_DI                                                              !
    ! 3) coords_DI                                                            !
    ! 4) gradient_DI                                                          !
    ! 5) gradient_DI_old                                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                     :: ndeg
    integer,intent(inout)                                  :: mode
    integer,intent(out)                                    :: neg
    integer,intent(out)                                    :: negreq
    real(kind=di_dp),dimension(:,:),intent(out)            :: hpad
    integer,intent(out)                                    :: status
    !-------------------------------------------------------------------------!
    real(kind=di_dp),external :: ddot

    real(kind=di_dp),parameter :: small=1.0e-8_di_dp
    real(kind=di_dp) :: eigmin,eigmax,scale,rmsd,rmsg,cx,cn,gx,dx
    real(kind=di_dp),dimension(:),allocatable             :: disp_DI_local
    real(kind=di_dp),dimension(:),allocatable             :: coords_DI_local
    real(kind=di_dp),dimension(:),allocatable             :: gradient_DI_local
    real(kind=di_dp),dimension(:),allocatable             :: gradient_DI_old_local
    real(kind=di_dp),dimension(:),allocatable             :: eigvals
    real(kind=di_dp),dimension(:,:),allocatable           :: h_local
    real(kind=di_dp),dimension(:,:),allocatable           :: hnew
    real(kind=di_dp),dimension(:),allocatable             :: work1
    real(kind=di_dp),dimension(1:100) :: xrcon

    real(kind=di_dp),dimension(1:100,1:3),save :: vlam_dgx
    real(kind=di_dp),save :: vdad,vgrd

    integer :: nsat,j,jj,ndg,iunsat,i,it,nc,ierr
    !-------------------------------------------------------------------------!

    if (di_iprint>4 .and. di_on_root) then
       write(di_stdout,*) 'di_algor_get_step_unsat_constr :'
       write(di_stdout,*) 'ndeg, ncycle, di_iprint, nunsat, ihess ',ndeg,ncycle,di_iprint,nunsat,ihess
    endif

    status = 0

    !  Check options and set up defaults

    eigmin = 0.0001_di_dp
    eigmax = 25.0_di_dp
    scale = one

    iunsat = nunsat
    negreq = iunsat
    ndg = ndeg + iunsat

    if (tsflag) negreq = iunsat + 1

    !  Allocate memory for intermediate storage

    allocate(disp_DI_local(1:ndg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating disp_DI_local in di_algor_get_step_unsat_constr')
    allocate(coords_DI_local(1:ndg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating coords_DI_local in di_algor_get_step_unsat_constr')
    allocate(gradient_DI_local(1:ndg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating gradient_DI_local in di_algor_get_step_unsat_constr')
    allocate(gradient_DI_old_local(1:ndg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating gradient_DI_old_local in di_algor_get_step_unsat_constr')
    allocate(eigvals(1:ndg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating eigvals in di_algor_get_step_unsat_constr')
    allocate(h_local(1:ndg,1:ndg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating h_local in di_algor_get_step_unsat_constr')
    allocate(hnew(1:ndg,1:ndg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating hnew in di_algor_get_step_unsat_constr')
    allocate(work1(1:(ndg+2)*ndg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating work1 in di_algor_get_step_unsat_constr')
    disp_DI_local = zero
    coords_DI_local = zero
    gradient_DI_local = zero
    gradient_DI_old_local = zero
    h_local = zero
    hnew = zero
    eigvals = zero

    do i = 1,ndeg
       disp_DI_local(i)    = disp_DI(i)
       coords_DI_local(i) = coords_DI(i)
       gradient_DI_local(i)   = gradient_DI(i)
       gradient_DI_old_local(i) = gradient_DI_old(i)
       do j = 1,ndeg
          h_local(i,j) = hess_DI(i,j)
       enddo
    enddo

    ! keep HESS as (NDEG-iunsat)**2 for update and change to cartesian
    ! if (iunsat/=0) expand Hess to (NDEG+iunsat)**2

    if (iunsat/=0) then
       nsat = 1
       if (ihess==0) then
          ! augment initial Hessian HESS  with diagonal forces of 
          ! constraints to make Ndeg-iunsat Hessian
          call di_algor_expand_hessian_1(ndeg-iunsat,ndeg,iunsat,xrcon,rcon(nsat,2:2), &
     &                 hess_DI,h_local,hpad,icc(nsat,17:17),gradient_DI_local)

         
       elseif (ihess==1) then
          ! pad Hess with previous Hess_uns strip
          ! extend gradient_DI,XINT,D
          call di_algor_expand_hessian_2(ndeg-iunsat,ndeg,iunsat,xrcon,rcon(nsat,2:2), &
     &                 hess_DI,h_local,hpad,coords_DI_local,gradient_DI_local)
          ! unsatisfied pieces 
          do j=1,iunsat
             jj = ndeg - iunsat + j
             disp_DI_local(jj)  = vlam_dgx(j,1)
             gradient_DI_old_local(jj) = vlam_dgx(j,2)
          enddo

          disp_DI_local(ndeg+iunsat)  = vdad
          gradient_DI_old_local(ndeg+iunsat) = vgrd
       endif

    elseif (iunsat==0) then

       do i = 1,ndeg
          do j = 1,ndeg
             h_local(i,j) = hess_DI(i,j)
          enddo
       enddo

    endif

    !  First update the Hessian if necessary

    if (ihess==1) call di_algor_update_hessian(ndeg,disp_DI_local,gradient_DI_local,gradient_DI_old_local,h_local)

    ! expand H (ndeg,ndeg) to  H(ndg,ndg)
    if (iunsat/=0) then

       call di_algor_expand_hessian_3(ndeg,ndg,ndeg-iunsat,iunsat,hpad,h_local,hnew) 

       ! copy back hnew into iu
       h_local = hnew

    elseif (iunsat==0) then
       do i = 1,ndeg
          do j = 1,ndeg
             hess_DI(i,j) = h_local(i,j)
          enddo
       enddo
    endif

 10 continue

    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,1000)
       call di_algor_print_matrix(ndg,ndg,ndg,h_local)
    endif

    !  Diagonalize

    call dsyev('v','l',ndg,h_local,ndg,eigvals,work1,(ndg+2)*ndg,status)
    if (status/=0) then
       if (di_on_root) write(di_stdout,1010)
       go to 999
    endif

    ! if Eigenvalues are very bad; = very small;
    ! restore the old Hessian, and make a new set of delocalized internals
    if (iunsat==0) then
       if (eigvals(1)<tp_bad_hes .and. (.not.TSflag)) then
          status= 2
          ! print eigenvalues
          if (di_on_root) then
             write(di_stdout,1120)
             write(di_stdout,1130) (eigvals(i),i=1,ndeg)
          endif
          go to 999
       endif
    endif

    !  check for and remove eigenvectors corresponding to zero
    !  eigenvalues and possible symmetry-redundant modes

    call di_algor_cleanup_hessian(ndg,iunsat,h_local,eigvals,gradient_DI_local,small,nc,status)
    if (status/=0) go to 999

    !  on exit, NC is the number of Hessian modes that will be
    !  used to form the next step.
    !
    !  check eigenvalue magnitudes are between allowed values

    call di_algor_adjust_hess_EVs_range(nc,eigmin,eigmax,eigvals)

    !  find the number of negative Hessian eigenvalues

    call di_algor_num_negative_hess_EVs(nc,eigvals,neg)

    !  form the RMS gradient

    rmsg = ddot(ndg,gradient_DI_local(1),1,gradient_DI_local(1),1)
    rmsg = sqrt(rmsg/float(nc))

    !  calculate the next step
    !  this can be done using either GDIIS or the
    !  standard EF algorithm

    !  modify iunsat gradients by adding constraint lamda
    if (iunsat/=0) then
       it = ndeg-iunsat
       do j=1,iunsat
          gradient_DI_local(it + j) = gradient_DI_local(it + j) + xrcon(j)
       enddo
    endif

    call di_algor_eigenvector_following(nc,ndg,iunsat,neg,negreq,mode,h_local,eigvals,gradient_DI_local,disp_DI_local,status)


    if (status==-1) then
       ! scale hess, hpad, and expand like di_algor_expand_hessian_3
       call di_algor_ef_check_H_scale_unsat(ndeg,ndeg-iunsat,iunsat,hpad,eigvals,scale)
       if (scale==one) go to 95

       ! expand 
       call di_algor_expand_hessian_4(ndeg,ndg,ndeg-iunsat,iunsat,hpad,h_local)

       if (scale/=one) go to 10
    endif
    if (status>0) go to 999

    !  we have a new step in D
    !  check the stepsize

    call di_algor_check_vector_length(ndeg,disp_DI_local)

    ! find xrcon (a new step for Lagrangian)
    if (iunsat/=0) then

       if (di_iprint>2 .and. di_on_root) write(di_stdout,1070)
 1070  format(/,21X,' Lagrange Multipliers for Constraints',/, &
     &       '       Constraint  Current Value    Gradient  Displacement   New Value')

       do j=1,iunsat
          it = ndeg+j
          gx = gradient_DI_local(it)
          dx = disp_DI_local(it)

          cx = coords_DI_local(it)
          cn = cx + dx
          if (di_iprint>2 .and. di_on_root) then
             write(di_stdout,1080) j,cx,gx,dx,cn
          endif

          ! a new step for Lagrangian
          xrcon(j) = cn
       enddo
 1080  format(10x,i3,7x,f11.6,4x,f9.6,4x,f9.6,2x,f11.6)

    endif ! Iunsat

    !    Standard Printout

    if (di_iprint>2) then

       !  form the RMS displacement

       rmsd = ddot(ndg,disp_DI_local(1),1,disp_DI_local(1),1)
       rmsd = sqrt(rmsd/float(nc))

       if (di_on_root) write(di_stdout,2000) ncycle,ec,rmsg,rmsd
    endif

    ! restore D
    do i = 1,ndeg
       disp_DI(i) = disp_DI_local(i) 
       ! restore gradient_cart
       gradient_DI(i) = gradient_DI_local(i)
       gradient_DI_old(i) = gradient_DI(i)
       coords_DI(i) = coords_DI(i) + disp_DI(i)
    enddo

    ! save lambda - d,gradient_DI_old
    ! save things that are < ndeg-iunsat
    if (iunsat/=0) then
       do j=1,iunsat
          jj = ndeg-iunsat+j   
          vlam_dgx(j,1) = disp_DI(jj)
          vlam_dgx(j,2) = gradient_DI_old(jj) 
          vlam_dgx(j,3) = coords_DI(jj) 
       enddo
       vdad = disp_DI_local(ndeg+iunsat)
       vgrd = gradient_DI_local(ndeg+iunsat)
    endif
    go to 999

 95 continue
    status = 3

 999 continue

    deallocate(disp_DI_local,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating disp_DI_local in di_algor_get_step_unsat_constr')
    deallocate(coords_DI_local,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating coords_DI_local in di_algor_get_step_unsat_constr')
    deallocate(gradient_DI_local,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating gradient_DI_local in di_algor_get_step_unsat_constr')
    deallocate(gradient_DI_old_local,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating gradient_DI_old_local in di_algor_get_step_unsat_constr')
    deallocate(h_local,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating h_local in di_algor_get_step_unsat_constr')
    deallocate(hnew,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hnew in di_algor_get_step_unsat_constr')
    deallocate(eigvals,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating eigvals in di_algor_get_step_unsat_constr')
    deallocate(work1,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating work1 in di_algor_get_step_unsat_constr')

    return

 1000 format(/,'   Hessian Matrix')
 1010 format(/,2X,'Error: ',/,' Unable to Diagonalize Hessian Matrix')
 2000 format(/,' ** Cycle ',I3,'  Energy ',F16.9,'   RMSG ',F8.5,'   RMSD ',F8.5,' **')
 1120 format('  Hessian Eigenvalues:')
 1130 format(1X,6F12.6)

  end subroutine di_algor_get_step_unsat_constr

!------------------------------------------------------

  subroutine di_algor_write_constr_di_data (ncons22,ncons0,ut,lidelp,xprim,savtor,ktyp, &
    &                ictyp,nile,status)
    !=========================================================================!
    ! If num_fixed_constr_rigid_body != 0, update icc and rcon                !
    ! and call deloc_utils_write_di_data                                      !
    !                                                                         !
    !                        [ficcwr in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ncons22 (in) : dimensions of icc and rcon                            !
    ! 2) ncons0 (in) : number of constraints                                  !
    ! 3) ut (in) : U matrix                                                   !
    ! 4) lidelp (in) : list of atom IDs and cell translations                 !
    ! 5) xprim (in) : primitive internals                                     !
    ! 6) savtor (in) : saved torsions                                         !
    ! 7) ktyp (in) : constraint type                                          !
    ! 8) ictyp (in) : list of constraints                                     !
    ! 9) nile(in) : internal description of iteration                         !
    !10) ierr(out)                                                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) coords_DI                                                            !
    ! 5) gradient_DI                                                          !
    ! 6) ncons                                                                !
    ! 7) num_fixed_constr_rigid_body                                          !
    ! 8) np_p                                                                 !
    ! 9) map_fixed_coords                                                     !
    !10) bmat_disconn_fragm                                                   !
    !11) hess_DI                                                              !
    !12) list_fixed_coords                                                    !
    !13) num_disconn_fragm                                                    !
    !14) icc                                                                  !
    !15) rcon                                                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                      :: ncons22
    integer,intent(in)                                      :: ncons0
    real(kind=di_dp),dimension(:,:),intent(in)              :: ut
    integer,dimension(:,:),intent(in)                       :: lidelp
    real(kind=di_dp),dimension(:),intent(in)                :: xprim
    real(kind=di_dp),dimension(:),intent(in)                :: savtor
    integer,dimension(:),intent(in)                         :: ktyp
    integer,dimension(:,:),intent(in)                       :: ictyp
    integer,intent(in)                                      :: nile
    integer,intent(out)                                     :: status
    !-------------------------------------------------------------------------!
    integer :: inn,i,j,nc,in,ika,ikb,ikc,kc,ia1,ia2,ia3,ia,ib1,ib2,ib3,jb,ierr
    integer :: ila,ilb,ilc,ld
    integer,dimension(:,:),allocatable            :: icc_c
    real(kind=di_dp),dimension(:,:),allocatable   :: rcon_c
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(ut,1)<intcor .and. size(ut,2)<ndeg) &
       &  call deloc_utils_io_abort('Error in di_algor_write_constr_di_data : Array UT is too small')
    if (size(lidelp,1)<4 .and. size(lidelp,2)<intcor) &
       &  call deloc_utils_io_abort('Error in di_algor_write_constr_di_data : Array LIDELP is too small')
    if (size(xprim,1)<intcor) &
       &  call deloc_utils_io_abort('Error in di_algor_write_constr_di_data : Array XPRIM is too small')
    if (size(savtor,1)<intcor) &
       &  call deloc_utils_io_abort('Error in di_algor_write_constr_di_data : Array SAVTOR is too small')
    if (size(ktyp,1)<intcor) &
       &  call deloc_utils_io_abort('Error in di_algor_write_constr_di_data : Array KTYP is too small')
#endif
    status = 0
    allocate(icc_c(1:ncons22,1:21),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating icc_c in di_algor_write_constr_di_data')
    icc_c = 0
    allocate(rcon_c(1:ncons22,1:2),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating rcon_c in di_algor_write_constr_di_data')
    rcon_c = zero

    ! loop over all num_fixed_constr_rigid_body
    inn = ncons - num_fixed_constr_rigid_body

    ! just copy
    if (inn>0) then
       do i = 1,inn
          do j=1,21
            icc_c(i,j) = icc(i,j)
          enddo
          do j=1,2
            rcon_c(i,j) = rcon(i,j)
          enddo
       enddo
    else
       inn = 0
    endif

    do i = 1,num_fixed_constr_rigid_body
       nc = inn + i
         
       in = ictyp(nc,1)
       call deloc_utils_unpack_atom_indices(lidelp(1,in), ika, ikb, ikc, kc)
       call deloc_utils_unpack_atom_indices(lidelp(2,in), ia1, ia2, ia3, ia)
       call deloc_utils_unpack_atom_indices(lidelp(3,in), ib1, ib2, ib3, jb)
       call deloc_utils_unpack_atom_indices(lidelp(4,in), ila, ilb, ilc, ld)


       ! update icc
       icc_c(nc,1) = kc
       icc_c(nc,2) = ika
       icc_c(nc,3) = ikb
       icc_c(nc,4) = ikc
       icc_c(nc,5) = ia 
       icc_c(nc,6) = ia1
       icc_c(nc,7) = ia2
       icc_c(nc,8) = ia3

       if (ktyp(in)==1) then
          jb  = 0
          ib1 = 0
          ib2 = 0
          ib3 = 0
       endif

       if (ktyp(in)==1 .or. ktyp(in)==3) then
          ld  = 0
          ila = 0
          ilb = 0
          ilc = 0
       endif

       icc_c(nc,9)  = jb 
       icc_c(nc,10) = ib1
       icc_c(nc,11) = ib2
       icc_c(nc,12) = ib3
       icc_c(nc,13) = ld 
       icc_c(nc,14) = ila
       icc_c(nc,15) = ilb
       icc_c(nc,16) = ilc
       icc_c(nc,17) = ktyp(in)
       icc_c(nc,18) = in
       icc_c(nc,19) = in
       ! satisfied
       icc_c(nc,20) = 1
       ! this is not right if num_fixed_constr_rigid_body changes , but does not hurt now!!!!!!!!!!
       icc_c(nc,21) = ncons0 + i

       ! what is the  rcon?
       rcon_c(nc,1) = xprim(in)
       rcon_c(ncons,2) = zero
    enddo
  
    ! save data
    call deloc_utils_write_di_data (np_p,ut(1:intcor,1:ndeg),lidelp(:,1:intcor),xprim(1:intcor), &
     &          savtor(1:intcor),ktyp(1:intcor),map_fixed_coords,bmat_disconn_fragm,ictyp, &
     &          icc_c,rcon_c,coords_DI,gradient_DI,hess_DI,list_fixed_coords,list_rigid_body_atoms,num_disconn_fragm,nile,status)

    if (di_iprint>4 .and. di_on_root) then
       do i = 1,ncons
          write(di_stdout,'(a,21i3)') 'icc_c ',(icc_c(i,j),j=1,21)
          write(di_stdout,*) 'rcon ', rcon_c(i,1),rcon_c(i,2)
       enddo
    endif

    deallocate(icc_c,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icc_c in ficwr')
    deallocate(rcon_c,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating rcon_c in ficwr')

    return

  end subroutine di_algor_write_constr_di_data
!------------------------------------------------------


  subroutine di_algor_init_constraints
    !=========================================================================!
    ! Print values of constraint bonds, angles, dihedrals.                    !
    ! memory partition:                                                       !
    !  icc (1-16) 4*4  Info up to torsions                                    !
    !  icc 17          ityp                                                   !
    !  icc 18,19       constraint maping to primitive                         !
    !  icc 20          0 unsatisfied, 1 satisfied                             !
    !  icc 21          initial constraint number                              !
    !  rcon(,1)        value of internal                                      !
    !  rcon(,2)        value of gradient                                      !
    !                                                                         !
    !                        [const_p in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments: none                                                         !
    ! 1) ncons (in) : number of constraints                                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) ncons_uns                                                            !
    ! 5) tp_con_sat                                                           !
    ! 6) tp_benlin                                                            !
    ! 7) ncons                                                                !
    ! 8) icc                                                                  !
    ! 9) rcon                                                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) tp_benlin                                                            !
    ! 2) icc                                                                  !
    ! 3) rcon                                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    !-------------------------------------------------------------------------!
    integer :: j,ityp,ia,ita,itb,itc,jb,isa,isb,isc,iwa,iwb,iwc,jc,jd,ixa,ixb,ixc,ibad
    real(kind=di_dp) :: almostlinear,qq,almost0,diff,d11,vali,rc,scale
    !-------------------------------------------------------------------------!

    nunsat = 0
    almostlinear = tp_benlin
    
    if (di_on_root) then
       write(di_stdout,'(/,a)') '  INTERNAL CONSTRAINTS'
       write(di_stdout,'(a,/)') ' #  Type     Target    Actual           Definition'
    endif

    do j=1,ncons

       ityp = icc(j,17)

       ia = icc(j,1)
       ita = icc(j,2)
       itb = icc(j,3)
       itc = icc(j,4)
    
       jb = icc(j,5)
       isa = icc(j,6)
       isb = icc(j,7)
       isc = icc(j,8)

       jc = icc(j,9)
       iwa = icc(j,10)
       iwb = icc(j,11)
       iwc = icc(j,12)

       jd = icc(j,13)
       ixa = icc(j,14)
       ixb = icc(j,15)
       ixc = icc(j,16)

       select case (ityp)

       case (1)
          ! bonds

          isa = isa - ita
          isb = isb - itb
          isc = isc - itc
          ita = 0
          itb = 0
          itc = 0
          icc(j,2) = ita
          icc(j,3) = itb
          icc(j,4) = itc
          icc(j,6) = isa
          icc(j,7) = isb
          icc(j,8) = isc
          call di_algor_find_internal_coord(qq,1, &
     &                  ia,ita,itb,itc, &
     &                  jb,isa,isb,isc, &
     &                  jc,iwa,iwb,iwc, &
     &                  jd,ixa,ixb,ixc, &
     &                  almostlinear, ibad )

          vali = qq*length_conv   ! convert to angstrom
          if (j>Ncons_uns) rcon(j,1) = qq
          rc = rcon(j,1)*length_conv   ! convert to angstrom

          scale = one
          call di_algor_print_constraint(j, 1, &
     &                   icc(j,1),ita,itb,itc, &
     &                   icc(j,5),isa,isb,isc, &
     &                   icc(j,9),iwa,iwb,iwc, &
     &                   icc(j,13),ixa,ixb,ixc, &
     &                   vali,rc, tp_con_sat*scale)


       case (3)
          ! bends
  
          ita = ita - isa
          itb = itb - isb
          itc = itc - isc
          iwa = iwa - isa
          iwb = iwb - isb
          iwc = iwc - isc
          isa = 0
          isb = 0
          isc = 0
          icc(j,2) = ita
          icc(j,3) = itb
          icc(j,4) = itc
          icc(j,6) = isa
          icc(j,7) = isb
          icc(j,8) = isc
          icc(j,10) = iwa
          icc(j,11) = iwb
          icc(j,12) = iwc
          ! a triangle    (ia, jb, jc)
          call di_algor_find_internal_coord(qq,  3, &
     &                  ia,ita,itb,itc, &
     &                  jb,isa,isb,isc, &
     &                  jc,iwa,iwb,iwc, &
     &                  jd,ixa,ixb,ixc, &
     &                  almostlinear, ibad )

          vali = qq*toang
          if (j>Ncons_uns) rcon(j,1) = qq
          rc = rcon(j,1)*toang

          scale = 10.0_di_dp
          call di_algor_print_constraint(j, 3, &
     &                   icc(j,1),ita,itb,itc, &
     &                   icc(j,5),isa,isb,isc, &
     &                   icc(j,9),iwa,iwb,iwc, &
     &                   icc(j,13),ixa,ixb,ixc, &
     &                   vali,rc, tp_con_sat*scale)

          if (ibad>=1) then
             ! I have an angle close to 0 or 180 within almostlinear 
             d11 = min(abs(qq-di_pi),qq)
             ! can I rescue run by changing almostlinear (to 0.01 deg)
             almost0 = 0.00016_di_dp

             if (d11>almost0) then
                almostlinear = almost0
                tp_benlin = almostlinear

                if (di_on_root) write (di_stdout,'(a,i4,a)') &
     &            'Warning: Constraint # ',j,' is close to linear bend; Optimization may be unstable'
             else
                ! failure
                if (di_on_root) then
                   write(di_stdout,*) ' Error: Constraint # ',j, ' is a linear bend'
                   write(di_stdout,*)' Currently you can impose such constraint in two ways:'
                   write(di_stdout,*)' 1)  use penalty method and optimization in Cartesians ,or'
                   write(di_stdout,*)' 2)  modify structure and constraint by at least 0.01 deg'
                   write(di_stdout,*) '     and use delocalized coordinates without symmetry'
                endif
                call deloc_utils_io_abort &
                & ('Error in di_algor_init_constraints (found linear bend constraint): Please modify input and submit again')
             endif
          endif

       case (5) 
          ! torsions
          ita = ita - isa
          itb = itb - isb
          itc = itc - isc
          iwa = iwa - isa
          iwb = iwb - isb
          iwc = iwc - isc
          ixa = ixa - isa - iwa
          ixb = ixb - isb - iwb
          ixc = ixc - isc - iwc
          isa = 0
          isb = 0
          isc = 0
          icc(j,2) = ita
          icc(j,3) = itb
          icc(j,4) = itc
          icc(j,6) = isa
          icc(j,7) = isb
          icc(j,8) = isc
          icc(j,10) = iwa
          icc(j,11) = iwb
          icc(j,12) = iwc
          icc(j,14) = ixa
          icc(j,15) = ixb
          icc(j,16) = ixc
          ! a rectangular ( ia, jb, jc, jd)
          call di_algor_find_internal_coord(qq,  5, &
     &                  ia,ita,itb,itc, &
     &                  jb,isa,isb,isc, &
     &                  jc,iwa,iwb,iwc, &
     &                  jd,ixa,ixb,ixc, &
     &                  almostlinear, ibad )

          vali = qq*toang
          if (j>Ncons_uns) rcon(j,1) = qq
          rc = rcon(j,1)*toang

          scale = 10.0_di_dp
          call di_algor_print_constraint(j, 5, &
     &                   icc(j,1),ita,itb,itc, &
     &                   icc(j,5),isa,isb,isc, &
     &                   icc(j,9),iwa,iwb,iwc, &
     &                   icc(j,13),ixa,ixb,ixc, &
     &                   vali,rc, tp_con_sat*scale)

          if (ibad==1)  then
             if (di_on_root) then
                write(di_stdout,*) ' Error: Constraint # ',j, ' has three or more atoms that are virtually linear'
                write(di_stdout,*) ' Error: Please eliminate this constraint and submit again'
             endif
             call deloc_utils_io_abort('Error in di_algor_init_constraints: found constraint with three nearly linear atoms')
          endif
       case default
          scale = -1.0_di_dp ! this is to satisfy ftncheck. anyway we are exiting....
          call deloc_utils_io_abort('Error in di_algor_init_constraints: unknown type of internal constraint')
       end select
       
       if (abs(vali-rc)>tp_con_sat*scale) then
          icc(j,20)=0
          nunsat = nunsat + 1 
          ! WARNING
          ! 0.1 A and 10 deg max change
          diff = abs(qq-rcon(j,1))
          if (diff>0.3_di_dp .and. di_on_root) then
             write(di_stdout,'(/,a)') ' Warning: Target value differs too much from Actual'
             write(di_stdout,*) ' Warning: Optimization may abort'
          endif
       else
          icc(j,20)=1
          rcon(j,1) = qq   
       endif

       ! register type and gradient
       icc(j,21) = j
       rcon(j,2) = qq - rcon(j,1)
 
    enddo ! ncons 
    nunsat0 = nunsat

    ! order constraints to put unsatisfied at the front
    call di_algor_order_unsat_constr(ncons)

    if (nunsat>100) call deloc_utils_io_abort(' Error in di_algor_init_constraints: too many unsatisfied constraints ')

    return

  end subroutine di_algor_init_constraints

!------------------------------------------------------

  subroutine  di_algor_find_internal_coord(qq,ityp,ia,ita,itb,itc,jb,isa,isb,isc, &
     &                         jc,iwa,iwb,iwc,jd,ixa,ixb,ixc,almostlinear,ibad )
    !=========================================================================!
    ! Find internal coordinate qq,                                            !
    ! ityp = 1 (bond), = 3 (angle), = 5 (torsion),                            !
    ! from cartesian coordinates coords_cart knowing atoms numbers ia,jb,jc,jd!
    ! and unit cell indices.                                                  !
    !                                                                         !
    !                        [confix in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) qq (out) : value of the internal coordinate                          !
    ! 2) ityp (in) : type of the constraint                                   !
    ! 3) ia (in) : atom ID                                                    !
    ! 4) ita (in) : cell translation along 'a', atom ia                       !
    ! 5) itb (in) : cell translation along 'b', atom ia                       !
    ! 6) itc (in) : cell translation along 'c', atom ia                       !
    ! 7) jb (in) : atom ID                                                    !
    ! 8) isa (in) : cell translation along 'a', atom jb                       !
    ! 9) isb (in) : cell translation along 'b', atom jb                       !
    !10) isc (in) : cell translation along 'c', atom jb                       !
    !11) jc (in) : atom ID                                                    !
    !12) iwa (in) : cell translation along 'a', atom jc                       !
    !13) iwb (in) : cell translation along 'b', atom jc                       !
    !14) iwc (in) : cell translation along 'c', atom jc                       !
    !15) jd (in) : atom ID                                                    !
    !16) ixa (in) : cell translation along 'a', atom jd                       !
    !17) ixb (in) : cell translation along 'b', atom jd                       !
    !18) ixc (in) : cell translation along 'c', atom jd                       !
    !19) almostlinear (in) : threshold for linear constraints                 !
    !20) ibad (out) : 0 if all is OK, 1 if linear bend/torsion is found       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) number_atoms                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    real(kind=di_dp),intent(out)                    :: qq
    integer,intent(in)                              :: ityp
    integer,intent(in)                              :: ia
    integer,intent(in)                              :: ita
    integer,intent(in)                              :: itb
    integer,intent(in)                              :: itc
    integer,intent(in)                              :: jb
    integer,intent(in)                              :: isa
    integer,intent(in)                              :: isb
    integer,intent(in)                              :: isc
    integer,intent(in)                              :: jc
    integer,intent(in)                              :: iwa
    integer,intent(in)                              :: iwb
    integer,intent(in)                              :: iwc
    integer,intent(in)                              :: jd
    integer,intent(in)                              :: ixa
    integer,intent(in)                              :: ixb
    integer,intent(in)                              :: ixc
    real(kind=di_dp),intent(in)                     :: almostlinear
    integer,intent(out)                             :: ibad
    !-------------------------------------------------------------------------!
    real(kind=di_dp),external :: ddot
    real(kind=di_dp),dimension(1:3) :: tt,ss,ww,xx,u,v,w,z,x
    real(kind=di_dp) :: co,si,cp,sj,qq1,qq2,s,r1,r2,r3
    !-------------------------------------------------------------------------!
    ibad = 0
    select case (ityp)
    
    case (1)

       ! bond
       call di_algor_unfold_atom_coords(coords_cart(1,ia),ita,itb,itc,tt)
       call di_algor_unfold_atom_coords(coords_cart(1,jb),isa,isb,isc,ss)
       call di_algor_distance(ss,tt,qq)
          
    case (3)

       ! angle
       call di_algor_unfold_atom_coords(coords_cart(1,ia),ita,itb,itc,tt)
       call di_algor_unfold_atom_coords(coords_cart(1,jb),isa,isb,isc,ss)
       call di_algor_unfold_atom_coords(coords_cart(1,jc),iwa,iwb,iwc,ww)

       call di_algor_norm_vector_difference(u,r1,tt,ss)
       call di_algor_norm_vector_difference(v,r2,ww,ss)
       co = ddot(3,u(1),1,v(1),1)
       si = di_algor_s2(co)

       if (abs(co)>one) co = sign(one,co)
       qq = acos(abs(co))
       if (abs(qq)<almostlinear) ibad = 1 
       qq = acos(co)

    case (5)
        ! torsion

        call di_algor_unfold_atom_coords(coords_cart(1,ia),ita,itb,itc,tt)
        call di_algor_unfold_atom_coords(coords_cart(1,jb),isa,isb,isc,ss)
        call di_algor_unfold_atom_coords(coords_cart(1,jc),iwa,iwb,iwc,ww)
        call di_algor_unfold_atom_coords(coords_cart(1,jd),ixa+iwa,ixb+iwb,ixc+iwc,xx)

        call di_algor_norm_vector_difference(u,r1,tt,ss)
        call di_algor_norm_vector_difference(v,r2,ww,ss)
        call di_algor_norm_vector_difference(w,r3,ww, xx)
        co = ddot(3,u(1),1,v(1),1)
        cp = ddot(3,v(1),1,w(1),1)
        si = di_algor_s2(co)
        sj = di_algor_s2(cp)

        !       --compute value of first planar bend (kc,ia,jb)
        if (abs(co)>one) co = sign(one,co)
        qq1 = acos(abs(co))
        !       --compute value of second planar bend (ia,jb,ld)
        if (abs(cp)>one) cp = sign(one,cp)
        qq2 = acos(abs(cp))

        if (abs(qq1)<almostlinear .or. abs(qq2)<almostlinear) ibad = 1 

        call di_algor_norm_vector_product(u,v,z)
        call di_algor_norm_vector_product(w,v,x)
        co = ddot(3,z(1),1,x(1),1)
        call di_algor_vector_product(z,x,u)
        si = sqrt(ddot(3,u(1),1,u(1),1))
        cp = ddot(3,u(1),1,v(1),1)

        s = di_algor_arc1(-co,si)
        if (cp<zero) s = -s
        qq = -s

     case default

        call deloc_utils_io_abort('Error in di_algor_find_internal_coord: unknown type of internal constraint')

     end select
          
     return

  end subroutine di_algor_find_internal_coord

!------------------------------------------------------

  subroutine di_algor_stop_fixed_constraints
    !=========================================================================!
    ! If all coordinates of the atom are fixed and                            !
    ! all the atoms of a particular internal are fixed, then                  !
    ! stop optimization.                                                      !
    !                                                                         !
    ! This routine is called on first entry to optimizer,                     !
    ! so will not capture FIX ATOMS (with internals all fixed).               !
    !                                                                         !
    ! The algorithm is to compare icc(1,5,9,13) - those are atom indices -    !
    ! with list_fixed_coords numbers. If all icc are fixed atoms, stop here.  !
    !                                                                         !
    !                        [confix in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:   none                                                       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) number_atoms                                                         !
    ! 5) ncons                                                                !
    ! 6) icc                                                                  !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    !-------------------------------------------------------------------------!
    integer :: i,j,ierr,kc,ia,jb,ld,iss,isave,istop
    integer, allocatable, dimension(:) :: iatm 
    !-------------------------------------------------------------------------!

    allocate(iatm(number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating iatm in di_algor_stop_fixed_constraints')
    iatm = 0

    ! make list of fixed atoms
    do i=1,number_atoms
       iatm(i)=0
       if ((list_fixed_coords(1,i)+list_fixed_coords(2,i)+list_fixed_coords(3,i))>0) iatm(i)=i
    enddo
        
    if (di_on_root) write(di_stdout,'(/)')

    istop = 0
    do j=1,ncons

       kc = icc(j,1)
       ia = icc(j,5)
       jb = icc(j,9)
       ld = icc(j,13)

       call di_algor_is_atom_in_list(kc,number_atoms,iatm,isave)
       iss = isave

       if (iss==0) call di_algor_is_atom_in_list(ia,number_atoms,iatm,isave)
       iss = iss + isave

       if (iss==0 .and. jb/=0) call di_algor_is_atom_in_list(jb,number_atoms,iatm,isave)
       iss = iss + isave

       if (iss==0 .and. ld/=0) call di_algor_is_atom_in_list(ld,number_atoms,iatm,isave)
       iss = iss + isave

       if (iss==0) then
          ! this constraint is build out of fixed atoms
          istop = 1
          if (di_on_root) write (di_stdout,'(a,i5,a)') ' Constraint # ',j,' has all atoms fixed'
       endif
    enddo

    deallocate(iatm,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating iatm in di_algor_stop_fixed_constraints')

    if (istop==1) call deloc_utils_io_abort('Error in input: Please modify constraints and/or fixed atoms')

    return
  end subroutine di_algor_stop_fixed_constraints
!------------------------------------------------------

  subroutine di_algor_tidy_constraints(newc,ncons0,ictyp,iab)
    !=========================================================================!
    ! If nfixcons changed that means some constraints are eliminated.         !
    ! Notice only constraints may be eliminated and not fixed cartesians.     !
    ! Rearrange icc, rcon,ictyp, redefine ncons                               !
    !                                                                         !
    !                        [tidycon in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) newc (in) : dimension for constraints arrays                         !
    ! 2) ncons0 (in) : number of "genuine"  constraints                       !
    ! 3) ictyp (inout) : list of which primitives to fix                      !
    ! 4) iab (in) : active internals left after Schmidt-Orthogonalization     !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) num_fixed_coords                                                     !
    ! 5) nunsat                                                               !
    ! 6) ncons                                                                !
    ! 7) icc                                                                  !
    ! 8) rcon                                                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) ncons                                                                !
    ! 2) icc                                                                  !
    ! 3) rcon                                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                   :: newc
    integer,intent(in)                                   :: ncons0
    integer,dimension(:,:),intent(inout)                 :: ictyp
    integer,dimension(:),intent(in)                      :: iab
    !-------------------------------------------------------------------------!
    integer,dimension(:,:),allocatable            :: icc_c
    real(kind=di_dp),dimension(:,:),allocatable   :: rcon_c
    integer,dimension(:),allocatable              :: ictyp_c
    integer :: ncons1,ii,i,ncold,isave,j,ierr,nfixcons_2
    !-------------------------------------------------------------------------!
    nfixcons_2 = size(iab)
    allocate(icc_c(1:nfixcons_2,1:21),stat=ierr)  
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating icc_c in di_algor_tidy_constraints')
    allocate(rcon_c(1:nfixcons_2,1:2),stat=ierr)    
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating rcon_c in di_algor_tidy_constraints')
    allocate(ictyp_c(1:nfixcons_2),stat=ierr)   
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating ictyp_c in di_algor_tidy_constraints')
    icc_c = 0
    rcon_c = zero
    ictyp_c = 0

    ncons1 = 0
    ! loop over the old ncons
    ii=0
    do i = 1,ncons
       ! the old constraint maps into primitive:
       ncold = ictyp(i,1)
       ! is this ncold present amongst iab or not?
       call di_algor_is_atom_in_list(ncold,newc,iab,isave)
       if (isave==0) then
          ! found update _c arrays
          ii=ii+1
          ictyp_c(ii) = ictyp(i,1)

          ! access only if genuine constraints and not fixed atoms constraints
          if (i<=ncons0) then
             ncons1 = ii
             do j=1,21
                icc_c(ii,j) = icc(i,j)
             enddo
                
             do j=1,2
                rcon_c(ii,j) = rcon(i,j)
             enddo
          endif
       endif
    enddo

    if (ii/=newc) call deloc_utils_io_abort('Internal error in di_algor_tidy_constraints - mismatch in number of constraints')

    ! redefine ncons
    ncons = newc

    ! reset arrays
    do i=1,ncons
       ictyp(i,1) = ictyp_c(i)
       ictyp(i,2) = ictyp_c(i)
    enddo

    if (ncons1/=0) then
       do i=1,ncons1
          do j=1,21
             icc(i,j) = icc_c(i,j)
          enddo
          do j=1,2
             rcon(i,j) = rcon_c(i,j)
          enddo
       enddo
    endif

    deallocate(icc_c,stat=ierr)  
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icc_c in di_algor_tidy_constraints')
    deallocate(rcon_c,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating rcon_c in di_algor_tidy_constraints')
    deallocate(ictyp_c,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating ictyp_c in di_algor_tidy_constraints')

    return
  end subroutine di_algor_tidy_constraints


  subroutine di_algor_schmidt_constraints(ndeg,intcor,ncatom,ncon,icon,rcon,icc21, &
     &                 ns,ut,icob,status)
    !=========================================================================!
    !  Modify the non-redundant set of natural internals                      !
    !  to account for any imposed constraints                                 !
    !                                                                         !
    !  What is done is to set up constraint vectors (with, e.g.,              !
    !  unit components corresponding to the primitives to be fixed)           !
    !  and Schmidt-Orthogonalize the current set of active                    !
    !  internals to these vectors                                             !
    !                                                                         !
    !                        storage for constraints:                         !
    ! [1 - nunsat] [ nunsat+1 - ncatom] [ncatom+1 - ncons]   [ncons+1- ncon]  !
    ! nonsatisfied     satisfied  num_fixed_constr_rigid_body num_fixed_coords!
    !                                                                         !
    !  Algorithm stops if constraints num_fixed_coords or nunsat              !
    !  are not acceptable.                                                    !
    !  If satisfied are skipped that's ok , just print message,               !
    !  and redefine pointers                                                  !
    !                                                                         !
    !                        [schmidt_p in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndeg (in) : number of active internals (number of degrees of freedom)!
    ! 2) intcor (in) : number of primitive internals                          !
    ! 3) ncatom (inout) : end of normal satisfied                             !
    ! 4) ncon (inout) : number of constrained (fixed) primitives              !
    !                   on exit actual number of constraints                  !
    ! 5) icon (in) : list of which primitives to fix                          !
    ! 6) rcon (in) : constraint coefficients for possible compound constraint !
    ! 7) icc21 (in) : original ordering of constraints                        !
    ! 8) ns (out) : active internals left after Schmidt-Orthogonalization     !
    ! 9) ut (inout) : on entry contains current set of natural internals      !
    !                 on exit contains Schmidt-Orthogonalized set             !
    !10) icob (out) : final arrangement of constraints                        !
    !11) ierr (out)                                                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    ! 3) di_iprint                                                            !
    ! 4) num_fixed_coords                                                     !
    ! 5) nunsat                                                               !
    ! 6) ncons                                                               !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                        :: ndeg
    integer,intent(in)                                        :: intcor
    integer,intent(inout)                                     :: ncatom
    integer,intent(inout)                                     :: ncon
    integer,dimension(:,:),intent(in)                         :: icon
    real(kind=di_dp),dimension(:),intent(in)                  :: rcon
    integer,dimension(:),intent(in)                           :: icc21
    integer,intent(out)                                       :: ns
    real(kind=di_dp),dimension(:,:),intent(inout)             :: ut
    integer,dimension(:),intent(out)                          :: icob
    integer,intent(out)                                       :: status
    !-------------------------------------------------------------------------!
    real(kind=di_dp),dimension(:,:),allocatable :: vm
    integer,dimension(:),allocatable :: icoa
    integer,dimension(:),allocatable :: icoc
    real(kind=di_dp),external :: ddot
    !real(kind=di_dp),parameter :: thcomp=1.0e-7_di_dp
    real(kind=di_dp) :: thrsh,thrsha,rval,proj,snorm,cf
    integer :: ister,isterun,istfix,istatom,i,ii,j,k,ncons_loc,is,ij,ncndm,ierr
    !-------------------------------------------------------------------------!
    ncndm = size(icob)
    allocate(vm(1:intcor,1:ncon+ndeg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating vm in di_algor_schmidt_constraints')
    allocate(icoa(1:ncndm),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating icoa in di_algor_schmidt_constraints')
    allocate(icoc(1:ncndm),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating icoc in di_algor_schmidt_constraints')
    vm = zero
    icoa = 0
    icoc = 0

    thrsh = tp_schmidt
    thrsha = tp_schmidt_a
    ister=0    
    isterun=0    
    istfix=0
    istatom = 0

    status = -1

    if (di_iprint>2 .and. di_on_root) write(di_stdout,1000)

    if (di_iprint>4 .and. di_on_root) write(di_stdout,*) &
     & ' schmidt: ndeg,intcor,nfix,ncons,ncatom,ncon,ncndm,nunsat ' &
     & , ndeg, intcor, num_fixed_coords, ncons, ncatom, ncon, ncndm ,nunsat

    !  project the constraints onto the active space

    do i=1,ncon
       ii = icon(i,2)
       rval = abs(rcon(i))
       if (rval==zero) rval = one
       if (di_iprint>2 .and. di_on_root)  write(di_stdout,*) 'constraints ',i,ii,rval

       do j=1,ndeg
          proj = ut(ii,j)*rval
          do k=1,intcor
             vm(k,i) = vm(k,i) + proj*ut(k,j)
          enddo
       enddo
    enddo

    !  assign to ncons_loc the total number of constraints
    !  (including compound constraints)

    ncons_loc = ncon

    !  normalize all the constraint vectors
    !  (at the same time eliminate zero constraint vectors)

    ii = 0
    do i=1,ncons_loc
       ii = ii+1
       icoa(ii)=icon(i,2)
       ! reuse icob
       icob(ii)= i
       snorm = ddot(intcor,vm(1,i),1,vm(1,i),1)

       if (snorm<thrsha) then

          !   [1 - nunsat] [ nunsat+1 - ncatom] [ncatom+1 - ncons] [ncons+1- ncon]
          !   nonsatisfied     satisfied   num_fixed_constr_rigid_body   num_fixed_coords 
          !   BAD if:
          !     nonsatisfied fail
          !     satisfied fail
          !     frozen atoms fail
          !   OK if
          !     fixed atoms (num_fixed_constr_rigid_body) fail if there are no num_fixed_coords 

          if (di_on_root) write(di_stdout,1100) icc21(i)
          ! stoping for bad unsatisfied
          if (i<=nunsat) isterun = isterun+1
          ! stoping for fixed coordinates zeros
          if (i>ncons) istfix = istfix+1
          ! identify  num_fixed_constr_rigid_body
          if (i>ncatom) istatom = istatom + 1

          ii = ii-1
       else
          snorm = one/sqrt(snorm)
          call dscal(intcor,snorm,vm(1,i),1)
          if (ii/=i) call dcopy(intcor,vm(1,i),1,vm(1,ii),1)
       endif
    enddo

    ! stoping if ii/=ncons_loc
    if (ii/=ncons_loc) ister=1

    !  number of constraints surviving is II

    if (ii/=ncons_loc) then
       if (di_iprint>2 .and. di_on_root) write(di_stdout,1350) ncons_loc-ii
    endif

    ncons_loc = ii


    if (di_iprint>5 .and. di_on_root) then
       ! projected constraints 
       write(di_stdout,1200)
       call di_algor_print_matrix(ncons_loc,intcor,ncons_loc,vm)
    endif

    icoc(1) = icoa(1)

    !  Constraints may not be independent
    !  see if any can be eliminated by Schmidt-Orthogonalization

    ii = 2
    do i=2,ncons_loc
       loop60 : do j=ii-1,1,-1

          cf = ddot(intcor,vm(1,j),1,vm(1,i),1)
          do k=1,intcor
             vm(k,i) = vm(k,i) - cf*vm(k,j)
          enddo

          snorm = ddot(intcor,vm(1,i),1,vm(1,i),1)

          if (snorm<thrsha) then
             if (di_on_root) write(di_stdout,1150) icc21(icob(i))
             ! stoping for bad unsatisfied
             if (icob(i)<=nunsat) isterun = isterun+1
             ! stoping for fixed coordinates  dependency
             if (icob(i)>ncons) istfix = istfix+1
             ! identify  num_fixed_constr_rigid_body
             if (icob(i)>ncatom) istatom = istatom + 1
             ister=2
             exit loop60
          endif

       enddo loop60

       !  normalize current constraint vector

       snorm = one/sqrt(snorm)
       call dscal(intcor,snorm,vm(1,i),1)
       if (ii/=i) call dcopy(intcor,vm(1,i),1,vm(1,ii),1)
       icoc(ii) = icoa(i)
       ii = ii+1

    enddo

    !  number of constraints surviving is II-1

    if (ii-1/=ncons_loc) then
       if (di_iprint>2 .and. di_on_root) write(di_stdout,1300) ncons_loc-ii+1
       ncons_loc = ii-1
    endif

    ! stoping if ii/=ncons_loc
    if (isterun>0) call deloc_utils_io_abort('Error: Please correct input constraints and submit again')

    ! stoping for problems with fixed coordinates
    if (istfix>0) call deloc_utils_io_abort(' Error: Please correct input fixed  coordinates and submit again ')

    ! stoping if istatom and num_fixed_coords
    if (istatom>0 .and. num_fixed_coords>0) call deloc_utils_io_abort(' Error: Please correct input fixed atoms and submit again ')

    if (ister>0) then
       if (di_on_root) write(di_stdout,'(/,a)') ' Message: Input constraints were modified; program will continue'
    endif

    ! how many fixed atom internals
    ncatom = istatom

    ! restore icob
    do i=1,ncons_loc
       icob(i)= icoc(i)
    enddo       

    if (di_iprint>5 .and. di_on_root) then
       ! final set of constraints 
       write(di_stdout,1400)
       call di_algor_print_matrix(ncons_loc,intcor,ncons_loc,vm)
    endif

    !  Now Schmidt-Orthogonalize the active internals
    !  to the remaining constraint vectors

    ii = ncons_loc+1
    do i=1,ndeg
       call dcopy(intcor,ut(1,i),1,vm(1,ii),1)

       do j=ii-1,1,-1

          cf = ddot(intcor,vm(1,j),1,vm(1,ii),1)
          do k=1,intcor
             vm(k,ii) = vm(k,ii) - cf*vm(k,j)
          enddo

          snorm = ddot(intcor,vm(1,ii),1,vm(1,ii),1)

       enddo

       !  normalize current internal coordinate

       if (di_iprint>4 .and. di_on_root) write(di_stdout,*) 'snorm ',i,ii,ndeg,snorm

       if (snorm>thrsh .and. ii<=ndeg) then
          snorm = one/sqrt(snorm)
          call dscal(intcor,snorm,vm(1,ii),1)
          ii = ii+1
       endif

    enddo

    !  copy orthogonalized active internals back into UT

    ns = ii-ncons_loc-1

    if (di_iprint>=2 .and. di_on_root) write(di_stdout,1600) ns

    do i=1,ns
       call dcopy(intcor,vm(1,i+ncons_loc),1,ut(1,i),1)
    enddo

    !  ------- WARNING   NEED TO REVISIT THIS SECTION ------
    !  now restore constraint vectors
    !  (these are needed in order to span the full space
    !   for the back-transformation)

    is = 0
    do i=1,ncons_loc
       ii = i + ns

       if (ii/=is) ut(:,ii) = zero
       is = ii

       ij = icob(i)
       ut(ij,ii) = one  
    enddo

    ncon = ncons_loc         ! final number of constraints

    !  normalize the final constraint vectors
    !  and copy them into UT

    do i=1,ncon
       ii = i+ns
       snorm = ddot(intcor,ut(1,ii),1,ut(1,ii),1)
       snorm = one/sqrt(snorm)
       call dscal(intcor,snorm,ut(1,ii),1)
    enddo

    !  possibly see what we've got

    if (di_on_root) then
       if (di_iprint>1) write(di_stdout,1700) ns,ncon
       if (di_iprint>5) then
          write(di_stdout,*) 'ncons,ns,ndeg,intcor', ncons_loc,ns,ndeg,intcor
          call di_algor_print_matrix(ns+ncon,intcor,ns+ncon,ut)
       endif
    endif

    !  set all small components to zero

    do i=1,ndeg
       do j=1,intcor
          if (abs(ut(j,i))<thrsh) ut(j,i) = zero
       enddo
    enddo

    !  check for error

    if (ns+ncon/=ndeg) then
       if (di_on_root) write(di_stdout,1800)
       status = 1
    else
       status = 0
    endif

    deallocate(vm,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating vm in di_algor_schmidt_constraints')
    deallocate(icoa,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icoa in di_algor_schmidt_constraints')
    deallocate(icoc,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating icoc in di_algor_schmidt_constraints')

    return

 1000 format(/,' Imposing Constraints by Schmidt Orthogonalization')
 1100 format(' Constraint ',i5,' is not active ')
 1150 format(' Constraint ',i5,' is not independent')
 1200 format(' Normalized Projected Constraint Vectors')
 1300 format(' Imposed Constraints are not Independent.  Eliminating ',i5,' Constraints')
 1350 format(' Imposed Constraints are not Active.  Eliminating ',i5,' Constraints')
 1400 format(' Final Set of Constraint Vectors')
 1600 format(' Number of Internal Coordinates Left in Active Space is ',i7)
 1700 format(/,' Schmidt-Orthogonalized Set of ',i7,' Active and ',i5,' Constraint Vectors')
 1800 format(/,2x,' Warning: Optimization Space has Incorrect Dimension',/, &
     &            '   after Schmidt-Orthogonalization of Constraints')

  end subroutine di_algor_schmidt_constraints

!------------------------------------------------------
  subroutine  di_algor_print_constraint(j,itype, &
     &   icc1,ita,itb,itc,icc5,isa,isb,isc,icc9,iwa,iwb,iwc, &
     &   icc13,ixa,ixb,ixc,vali,rc,thrsh)
    !=========================================================================!
    ! Print given constraint as specifed by its type, value, and atoms        !
    ! involved in the constraint.                                             !
    !                                                                         !
    !                        [pr_consat in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) j (in) : number of the constraint                                    !
    ! 2) itype (in) : type of the constraint                                  !
    ! 3) icc1 (in) : atom ID                                                  !
    ! 4) ita (in) : cell translation along 'a', atom icc1                     !
    ! 5) itb (in) : cell translation along 'b', atom icc1                     !
    ! 6) itc (in) : cell translation along 'c', atom icc1                     !
    ! 7) icc5 (in) : atom ID                                                  !
    ! 8) isa (in) : cell translation along 'a', atom icc5                     !
    ! 9) isb (in) : cell translation along 'b', atom icc5                     !
    !10) isc (in) : cell translation along 'c', atom icc5                     !
    !11) icc9 (in) : atom ID                                                  !
    !12) iwa (in) : cell translation along 'a', atom icc9                     !
    !13) iwb (in) : cell translation along 'b', atom icc9                     !
    !14) iwc (in) : cell translation along 'c', atom icc9                     !
    !15) icc13 (in) : atom ID                                                 !
    !16) ixa (in) : cell translation along 'a', atom icc13                    !
    !17) ixb (in) : cell translation along 'b', atom icc13                    !
    !18) ixc (in) : cell translation along 'c', atom icc13                    !
    !19) vali (in) : current value of the constrained internal                !
    !20) rc (in) : prescribed value of the constrained internal               !
    !21) thrsh (in) : threshold for deciding whether constrain is satisfied   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root                                                           !
    ! 2) di_stdout                                                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)            :: j
    integer,intent(in)            :: itype
    integer,intent(in)            :: icc1
    integer,intent(in)            :: ita
    integer,intent(in)            :: itb
    integer,intent(in)            :: itc
    integer,intent(in)            :: icc5
    integer,intent(in)            :: isa
    integer,intent(in)            :: isb
    integer,intent(in)            :: isc
    integer,intent(in)            :: icc9
    integer,intent(in)            :: iwa
    integer,intent(in)            :: iwb
    integer,intent(in)            :: iwc
    integer,intent(in)            :: icc13
    integer,intent(in)            :: ixa
    integer,intent(in)            :: ixb
    integer,intent(in)            :: ixc
    real(kind=di_dp),intent(in)   :: vali
    real(kind=di_dp),intent(in)   :: rc
    real(kind=di_dp),intent(in)   :: thrsh
    !-------------------------------------------------------------------------!
    if (.not. di_on_root) return

    if (abs(vali-rc)>thrsh) then
       ! nonsatisfied
        
       if (itype==5) then
          write(di_stdout,'(i3,a,f9.2,f10.2,1x,4(i3,a,3i2,a),a)') j,' Torsion', &
     &   rc, vali, &
     &   icc1,' (',ita,itb,itc,')', &
     &   icc5,' (',isa,isb,isc,')', &
     &   icc9,' (',iwa,iwb,iwc,')', &
     &   icc13,' (',ixa,ixb,ixc,')' 
       endif

       if (itype==3) then
          write(di_stdout,'(i3,a,f9.2,f10.2,1x,3(i3,a,3i2,a),a)') j,' Angle  ',  &
     &   rc, vali, &
     &   icc1,' (',ita,itb,itc,')', &
     &   icc5,' (',isa,isb,isc,')', &
     &   icc9,' (',iwa,iwb,iwc,')'
        endif

       if (itype==1) then
          write(di_stdout,'(i3,a,f9.3,f10.3,1x,2(i3,a,3i2,a),a)') j,' Bond   ', &
     &   rc, vali, &
     &   icc1,' (',ita,itb,itc,')', &
     &   icc5,' (',isa,isb,isc,')'
        endif

     else 
       ! satisfied

       if (itype==5) then
          write(di_stdout,'(i3,a,f9.2,a,1x,4(i3,a,3i2,a),a)') j,' Torsion', &
     &   vali, ' Satisfied', &
     &   icc1,' (',ita,itb,itc,')', &
     &   icc5,' (',isa,isb,isc,')', &
     &   icc9,' (',iwa,iwb,iwc,')', &
     &   icc13,' (',ixa,ixb,ixc,')'
       endif

       if (itype==3) then
          write(di_stdout,'(i3,a,f9.2,a,1x,3(i3,a,3i2,a),a)') j,' Angle  ', &
     &   vali, ' Satisfied', &
     &   icc1,' (',ita,itb,itc,')', &
     &   icc5,' (',isa,isb,isc,')', &
     &   icc9,' (',iwa,iwb,iwc,')'
       endif

       if (itype==1) then
          write(di_stdout,'(i3,a,f9.3,a,1x,2(i3,a,3i2,a),a)') j,' Bond   ', &
     &   vali, ' Satisfied', &
     &   icc1,' (',ita,itb,itc,')', &
     &   icc5,' (',isa,isb,isc,')'
       endif

    endif
    return

  end subroutine di_algor_print_constraint

!------------------------------------------------------
  subroutine di_algor_order_unsat_constr(ncons)
    !=========================================================================!
    ! Order unsatisfied constraints                                           !
    ! icc(20) has 0 (unsatisfied);  and 1 (satisfied)                         !
    ! icc is ordered ; first are unsatisfied, last satisfied                  !
    ! order icc  so 0 are at the end                                          !
    !                                                                         !
    !                        [orderuns in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ncons (in) : number of constraints                                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) coords_DI                                                            !
    ! 2) di_iprint                                                            !
    ! 3) di_on_root                                                           !
    ! 4) di_stdout                                                            !
    ! 5) disp_DI                                                              !
    ! 6) coords_cart                                                          !
    ! 7) ndeg                                                                 !
    ! 8) number_atoms                                                         !
    ! 9) icc                                                                  !
    !10) rcon                                                                 !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) icc                                                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                            :: ncons
    !-------------------------------------------------------------------------!
    integer,dimension(1:21) :: is
    integer :: j,i,i1,l
    real(kind=di_dp) :: r1,r2
    !-------------------------------------------------------------------------!

    if (ncons==1) return

    j=0
    do i1=1,ncons
       j=j+1
       if (icc(j,20)==1) then
          is(1:21) = icc (j,1:21)
          r1    = rcon(j,1)
          r2    = rcon(j,2)

          ! shift everything and put 0 at the end
          do l=j,ncons-1
             icc (l,1:21) = icc (l+1,1:21)
             rcon(l,1) = rcon(l+1,1)
             rcon(l,2) = rcon(l+1,2)
          enddo

          icc (ncons,1:21) = is(1:21)
          rcon(ncons,1) = r1
          rcon(ncons,2) = r2
          if (icc(j,20)==1) j=j-1
       endif
    enddo

    if (di_iprint>4 .and. di_on_root) then
       do i = 1,ncons
          write(di_stdout,'(a,21i3)') 'icc ',(icc(i,j),j=1,21)
          write(di_stdout,*) 'rcon ', rcon(i,1),rcon(i,2)
       enddo

    endif

    return
  end subroutine di_algor_order_unsat_constr


!--------------------------------------------------------------------
  subroutine di_algor_back_trans(ndim,ut,lidelp,xprim,ktyp,coords_cart0,hess,status)
    !=========================================================================!
    !  Conversion of Internal Coordinates back into Cartesians                !
    !  This is done iteratively according to:                                 !
    !                                                                         !
    ! coords_cart(K) = coords_cart(K-1) + B**(-1)*(coords_DI - coords_DI(K-1))!
    !                                                                         !
    !                    [getcard in DMol]                                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndim (in) : dimension (3*number_atoms)                               !
    ! 2) ut (in) : transformation matrix (which linear combination of         !
    !              primitive internals make up each compound delocalized      !
    !              coordinate                                                 !
    ! 3) lidelp (in) : list of atoms and cells involved in each primitive     !
    ! 4) xprim (out) : values of primitive internals                          !
    ! 5) ktyp (inout) : integer array containing internal coordinate type     !
    ! 6) coords_cart0 (in) : old cartesian coordinates                        !
    ! 7) hess (out) : Hessian in Cartesian                                    !
    ! 8) ierr (out)                                                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) coords_DI                                                            !
    ! 2) di_iprint                                                            !
    ! 3) di_on_root                                                           !
    ! 4) di_stdout                                                            !
    ! 5) disp_DI                                                              !
    ! 6) coords_cart                                                          !
    ! 7) ndeg                                                                 !
    ! 8) number_atoms                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) coords_cart                                                          !
    ! 2) coords_DI                                                            !
    ! 3) disp_DI                                                              !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                           :: ndim
    real(kind=di_dp),dimension(:,:),intent(in)                   :: ut
    integer,dimension(:,:),intent(in)                            :: lidelp
    real(kind=di_dp),dimension(:)                                :: xprim
    integer,dimension(:),intent(inout)                           :: ktyp
    real(kind=di_dp),dimension(:,:),intent(in)                   :: coords_cart0
    real(kind=di_dp),dimension(1:ndim,1:ndim),intent(out)        :: hess
    integer,intent(out)                                          :: status
    !-------------------------------------------------------------------------!
    real(kind=di_dp),external :: ddot
    logical :: cnvgd
  !  real(kind=di_dp),parameter :: small=1.0e-8_di_dp
    integer,parameter :: maxit=60
    integer,parameter :: itrdmax=3
    real(kind=di_dp),parameter :: onet=0.1_di_dp
    integer :: nat3,nic,it,itor,idmax,ifinish,ijump,i,j,k,ierr
    real(kind=di_dp) :: scald,scaldi,dmax0,dmax1,scald01,dci,disc,disi,disi0,xcmax
    real(kind=di_dp),dimension(:,:),allocatable :: coords_cart_store ! 1:3,1:natoms
    real(kind=di_dp),dimension(:,:),allocatable :: disp_cart ! 1:3,1:natoms
    real(kind=di_dp),dimension(:),allocatable   :: xdint_store ! 1:ndeg
    real(kind=di_dp),dimension(:,:),allocatable :: bmat ! 1:12,1:intcor
    real(kind=di_dp),dimension(:,:),allocatable :: bnew ! 1:3*natoms,1:ndeg
    real(kind=di_dp),dimension(:,:),allocatable :: bmbt ! 1:ndeg,1:ndeg
    real(kind=di_dp),dimension(:,:),allocatable :: binv ! 1:3*natoms,1:ndeg
    real(kind=di_dp),dimension(:),allocatable   :: eigv ! 1:ndeg
    real(kind=di_dp),dimension(:),allocatable   :: coords_cart_change ! 1:3*natoms
    integer,dimension(:,:),allocatable          :: indexb ! 1:12,1:intcor
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(ut,1)<intcor .or. size(ut,2)<ndeg) &
      & call deloc_utils_io_abort('Error in di_algor_back_trans: Array UT is too small')
    if (size(lidelp,1)<4 .or. size(lidelp,2)<intcor) &
      & call deloc_utils_io_abort('Error in di_algor_back_trans: Array LIDELP is too small')
    if (size(xprim,1)<intcor) &
      & call deloc_utils_io_abort('Error in di_algor_back_trans: Array XPRIM is too small')
    if (size(ktyp,1)<intcor) &
      & call deloc_utils_io_abort('Error in di_algor_back_trans: Array KTYP is too small')
#endif

    disi0=-1 ! to satisfy FTNCHECK: I've checked it can not be used uninitialised 21/Apr/2006 aperlov
    status = 0
    nat3 = 3*number_atoms
    nic = intcor

    ! allocate memory for internal needs
    allocate(coords_cart_store(1:3,1:number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating coords_cart_store in di_algor_back_trans')
    coords_cart_store = zero
    allocate(disp_cart(1:3,1:number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating disp_cart in di_algor_back_trans')
    disp_cart = zero
    allocate(xdint_store(1:ndeg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating xdint_store in di_algor_back_trans')
    xdint_store = zero
    allocate(bmat(1:12,1:nic),stat=ierr) ! zz(ib)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating bmat in di_algor_back_trans')
    bmat = zero
    allocate(indexb(1:12,1:nic),stat=ierr) ! zz(ibx)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating indexb in di_algor_back_trans')
    indexb = 0
    allocate(bnew(1:nat3,1:ndeg),stat=ierr) ! zz(ig)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating bnew in di_algor_back_trans')
    bnew = zero
    allocate(bmbt(1:ndeg,1:ndeg),stat=ierr) ! zz(isc1)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating bmbt in di_algor_back_trans')
    bmbt = zero
    allocate(binv(1:nat3,1:ndeg),stat=ierr) ! zz(ibv)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating binv in di_algor_back_trans')
    binv = zero
    allocate(eigv(1:ndeg),stat=ierr) ! zz(ieig)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating eigv in di_algor_back_trans')
    eigv = zero
    allocate(coords_cart_change(1:nat3),stat=ierr) ! zz(isc6)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating coords_cart_change in di_algor_back_trans')
    coords_cart_change = zero

    it = 0
    itor = 2
    Scald = one   
    Scaldi = one
    dmax0 = zero 
    idmax = 0
    ifinish=0
    ijump=0

    if (di_iprint>1 .and. di_on_root) write(di_stdout,1000)

    ! save original coords_cart
    coords_cart_store = coords_cart

    ! save original coords_DI-disp_DI(i)
    do i=1,ndeg
       xdint_store(i)=coords_DI(i)-disp_DI(i)
    enddo
    go to 100

 50 continue

    ! start again with the initial coords_cart and a small scald
    scald = onet
    ifinish = 1

    if (di_iprint>4 .and. di_on_root) write(di_stdout,*) 'scale factor for back transformation ',scald

    ! restore original coords_cart
    coords_cart = coords_cart_store

    ! modify coords_DI
    do i=1,ndeg
       coords_DI(i) = xdint_store(i) + scald*disp_DI(i)
    enddo

    dmax0 = zero 
    idmax = 0

 100 continue
    it = it+1

    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,'(/,a)')'Values of Cartesian coordinates1'
       do i=1,number_atoms
          write(di_stdout,'(i5,3f20.14)') i,(coords_cart(j,i),j=1,3)
       enddo
    endif

    !  form the B-matrix
    !
    !  using coords_cart cartesian coodinates form B matrix (bmat) and 
    !  internal primitives in xprim

    call di_algor_make_B_matrix(itor,ktyp,bmat,indexb,xprim,lidelp,status)

    if (status/=0) then
       if (di_on_root) write(di_stdout,*) ' Warning: error in di_algor_make_B_matrix '
       ip_error = 1
       go to 999
    endif

    !  transform the original B-Matrix to the new coordinates
    call di_algor_U_transform_B_matrix(ut,bnew,bmat,indexb)

    !  new way to invert the new B-Matrix ; find (b*bt(1-)*b in z(ib)

    call di_algor_invert_B_matrix(bnew,bmbt,binv,status)
    if (status/=0) then
       if (di_on_root) write(di_stdout,1200)
       ip_error = 11
       go to 999
    endif

    !  get the internal coordinates from the primitives

    do i=1,ndeg
       eigv(i) = ddot(intcor,ut(1,i),1,xprim(1),1)
    enddo

    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,'(a)')'Values of Primitive Internals'
       do i=1,intcor
          write(di_stdout,'(i5,f20.14,i5)') i,xprim(i),ktyp(i)
       enddo
       write(di_stdout,'(a/)')'END Values of Primitive Internals'
       write(di_stdout,'(a)')'Values of Delocalized Internals'
       do i=1,ndeg
          write(di_stdout,'(i5,f20.14)') i,eigv(i)
       enddo
       write(di_stdout,'(a)')'END Values of Delocalized Internals'
    endif ! (di_iprint>5)

    !  check for convergence

    call di_algor_back_trans_converged(it,ndeg,eigv,cnvgd,dmax1)

    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,*) 'cnvgd,scald,ifinish,scald,idmax,dmax0,dmax1 ', cnvgd,scald,ifinish,scald,idmax,dmax0,dmax1
    endif

    if (cnvgd) then
       scald01 = scald + onet 
       if (scald01>=one .or. ifinish==2) go to 95

       ! save the current success ; increase scald and go to 100
       ! save original coords_cart
       coords_cart_store = coords_cart

       if (di_iprint>4 .and. di_on_root) write(di_stdout,*) 'Back transformation successful ',scald
       scald = scald + onet 
       idmax = 0
       dmax0 = zero
       ifinish = 1

       ! modify coords_DI
       do i=1,ndeg
          coords_DI(i) = xdint_store(i) + scald*disp_DI(i)
       enddo
       go to 100
    endif  ! cnvgd

    ! dmax increased or not; are we progressing here or not
    if (dmax1>=dmax0) idmax = idmax + 1 
    dmax0 = dmax1

    ! too many failures; reset scaling 
    if (idmax>itrdmax) then

       ! if this is a beginning go to 50
       if (ifinish==0) go to 50

       ! otherwise go back to the latest success
       ! restore original coords_cart
       coords_cart = coords_cart_store
       scald = scald - onet 

       ! modify coords_DI
       do i=1,ndeg
          coords_DI(i) = xdint_store(i) + scald*disp_DI(i)
       enddo
       ifinish=2
  
       if (di_iprint>4 .and. di_on_root) write(di_stdout,*) 'back transformation done ', scald
       if (scald>zero) then
          go to 100
       else
          if (di_on_root) write(di_stdout,1200)
          status = 10
          go to 999
       endif

    endif ! idmax>itrdmax

    if (it>maxit .and. ifinish==0) then
       if (di_on_root) write(di_stdout,1200)
       status = 20
       go to 999
    endif

    !  generate new Cartesian coordinates in coords_cart

    disp_cart = coords_cart
    call di_algor_get_cart_from_DI(ndeg,binv,eigv)

    ! save displacement as NEW-OLD cartesian coordinates
    disp_cart = coords_cart - disp_cart

    ! SYMMETRY
    if (ntrans>1) then

       ! solid
       if (ldokpt) then
          if (di_iprint>2 .and. di_on_root) write(di_stdout,'(a)') 'Symmetrizing coordinates '
    
          ! put previous coordinates into coords_cart, symmetrize and add displacements
          coords_cart = coords_cart - disp_cart
          call di_algor_symmetrize_cartesian(disp_cart)
          coords_cart = coords_cart + disp_cart

       endif

       ! molecule
       if (.not.period) then
          if (di_iprint>2) write(di_stdout,'(a)')'Symmetrizing coordinates '

          !vmilman 
          !commented out (molecular case)

          !call symvec(ntrans,neqatm,trans,small)
       endif
    endif ! ntrans>1


    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,'(/,a)')'Values of Cartesian coordinates'
       do i=1,number_atoms
          write(di_stdout,'(i5,3f20.14)') i,(coords_cart(j,i),j=1,3)
       enddo ! i
    endif

 
    !  go back and get new internals

    go to 100

 95 continue

    if (scald/=one .and. di_iprint>4 .and. di_on_root) write(di_stdout,'(/,a,f10.6)') 'Final scaling factor used: ',scald

    ! redefine disp_DI
    if (scald/=one) then
       do i=1,ndeg
          disp_DI(i) = coords_DI(i) - xdint_store(i) 
       enddo
    endif


    if (ijump<=2) then
       ! check if the changes in Cartesians are too large
       ! if so decrease disp_DI by half and redo search again
       k=0
       do i=1,3
          do j=1,number_atoms
             k=k+1
             coords_cart_change(k) = coords_cart(i,j)-coords_cart0(i,j)
          enddo
       enddo
       !  check validity of new Cartesians by comparing internal
       !  and Cartesian displacements

       disi = sqrt(ddot(ndeg,disp_DI(1),1,disp_DI(1),1))
       disc = sqrt(ddot(nat3,coords_cart_change(1),1,coords_cart_change(1),1))
       dci=100.0_di_dp

       if (di_iprint>1 .and. di_on_root) then
          write(di_stdout,'(a,f10.6)') ' Norm of Displacement of Delocalized Coordinates: ',disI
          write(di_stdout,'(a,f10.6)') ' Norm of Displacement of Cartesian Coordinates:   ',disC
       endif
    
       ! maximum Cartesian change
       xcmax = zero
       do i=1,nat3
          if (abs(coords_cart_change(i))>xcmax) xcmax = abs(coords_cart_change(i))
       enddo

       if (ijump==0) disi0 = disi

       if (disc>dci*disi0 .or. xcmax>tp_back_car) then
          coords_cart = coords_cart0
          scaldi = scaldi*0.5_di_dp 
 
          if (di_iprint>4 .and. di_on_root) write(di_stdout,'(/,a,f10.6)') '     Scaling   disp_DI ',scaldi
          do i=1,ndeg
             disp_DI(i) = scaldi*disp_DI(i)
          enddo
          ijump=ijump+1
          scald = one
          ifinish = 0
          dmax0 = zero  
          idmax = 0

          ! restore original coords_cart
          coords_cart = coords_cart_store

          ! modify coords_DI
          do i=1,ndeg
             coords_DI(i) = xdint_store(i) + scald*disp_DI(i)
          enddo
          go to 100
       endif

    else

       if (di_on_root) write(di_stdout,1200)
       status = 30
       go to 999

    endif ! ijump=2

    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,'(/,a)')'Values of Cartesian coordinates'
       do i=1,number_atoms
          write(di_stdout,'(i5,3f20.14)') i,(coords_cart(j,i),j=1,3)
       enddo ! i
    endif

    call di_algor_cart_hessian_from_DI(ndeg,ns,bnew,hess)
    if (di_iprint>4 .and. di_on_root) write(di_stdout,*) 'j: hess_DI transformed into hess'

 999 continue

    if (status/=0) write(di_stdout,*) 'status: ',status
    ! free memory
    deallocate(coords_cart_store,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating coords_cart_store in di_algor_back_trans')
    deallocate(disp_cart,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating disp_cart in di_algor_back_trans')
    deallocate(xdint_store,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating xdint_store in di_algor_back_trans')
    deallocate(bmat,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bmat in di_algor_back_trans')
    deallocate(indexb,stat=ierr) 
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating indexb in di_algor_back_trans')
    deallocate(bnew,stat=ierr) 
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bnew in di_algor_back_trans')
    deallocate(bmbt,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bmbt in di_algor_back_trans')
    deallocate(binv,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating binv in di_algor_back_trans')
    deallocate(eigv,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating eigv in di_algor_back_trans')
    deallocate(coords_cart_change,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating coords_cart_change in di_algor_back_trans')

    return

 1000 format(/,' New Cartesian Coordinates Obtained by Inverse Iteration')
 1200 format(/,2x,' Warning: Exceeded allowed number of iterative cycles in GetCART'/2x, &
     & ' an internal coordinate may have approached 180 degrees,',/, &
     & '   new set of delocalized internals will be generated')

  end subroutine di_algor_back_trans

!--------------------------------------------------------------------

  subroutine di_algor_symmetrize_cartesian(carts)
    !=========================================================================!
    ! This applies the symmetry operations in locally held 'trans' array      !
    ! to cartesian coordinates.                                               !
    !                                                                         !
    !                   [symmveclattop in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) carts(inout) array to symmetrize                                     !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) number_atoms                                                         !
    ! 2) ntrans                                                               !
    ! 3) trans                                                                !
    ! 4) neqatm                                                               !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:  none                                 !
    !-------------------------------------------------------------------------!
    ! Modules used:  none                                                     !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! 1) work : store the unsymmetrised forces                                !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    implicit none
    real(kind=di_dp),dimension(:,:),intent(inout)       :: carts
    !-------------------------------------------------------------------------!
    integer       :: iatom,nsym,ierr
    real(kind=di_dp),dimension(:,:),allocatable :: work
    !-------------------------------------------------------------------------!

    if (ntrans==0) return

    ! Allocate work array
    allocate(work(1:3,1:number_atoms),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating work in di_algor_symmetrize_cartesian')

    ! ** Store the unsymmetrised coordinates in work

    work = carts

    ! ** Zero the coordinates

    carts = zero

    ! ** Apply the symmetrised average

    do iatom = 1,number_atoms

       do nsym = 1,ntrans

          carts(:,iatom) = carts(:,iatom) + &
                  matmul(transpose(trans(:,:,nsym)),work(:,neqatm(iatom,nsym)))          

       end do
    end do

    ! ** Scale by the number of symmetry operations

    carts = carts/ntrans

    ! Deallocate work array
    deallocate(work,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating work in di_algor_symmetrize_cartesian')

    return
  end subroutine di_algor_symmetrize_cartesian

!--------------------------------------------------------
  subroutine di_algor_back_trans_converged(iteration,nic,qq,cnvgd,disp_max)
    !=========================================================================!
    !  Checks for convergence during iterative back-transformation            !
    !  of natural internal coordinates to Cartesians.                         !
    !  Convergence is attained when each internal coordinate                  !
    !  generated from the current estimate of the corresponding               !
    !  Cartesian coordinates differs by less than tp_backz from its           !
    !  actual (known) value.                                                  !
    !                                                                         !
    !                         [cnvback_p in DMol]                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) iteration (in) : current cycle number                                !
    ! 2) nic (in) : number of internal coordinates                            !
    ! 3) qq (in) : internals as generated from current estimate of Cartesians !
    ! 4) cnvgd (out) : logical flag indicating convergence                    !
    ! 5) disp_max (out) : maximum deviation                                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) coords_DI                                                            !
    ! 2) di_iprint                                                            !
    ! 3) di_on_root                                                           !
    ! 4) di_stdout                                                            !
    ! 5) tp_backz                                                             !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/03/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                  :: iteration
    integer,intent(in)                                  :: nic
    real(kind=di_dp),dimension(1:nic),intent(in)        :: qq
    logical,intent(out)                                 :: cnvgd
    real(kind=di_dp),intent(out)                        :: disp_max
    !-------------------------------------------------------------------------!
    integer :: i
    real(kind=di_dp) :: dx
    !-------------------------------------------------------------------------!
    if (di_iprint>5 .and. di_on_root) write(di_stdout,1000)

    disp_max = abs(coords_DI(1)-qq(1))

    do i=1,nic
       dx = coords_DI(i)-qq(i)
       if (abs(dx)>disp_max) disp_max = abs(dx)
       if (di_iprint>5 .and. di_on_root) write(di_stdout,1100) i,coords_DI(i),qq(i),dx
    enddo

    if (di_iprint>3 .and. di_on_root) write(di_stdout,1200) iteration,disp_max

    !  Converged?

    cnvgd = disp_max<tp_backz

    return

 1000 format(/,5x,'Iterative generation of Cartesian Coordinates',/, &
     &       '     Internal   True Value     Estimate     Difference')
 1100 format(5x,i5,2x,3(2x,f12.8))
 1200 format(5x,'Cycle: ',i3,'  Maximum deviation: ',f12.8)

  end subroutine di_algor_back_trans_converged
!--------------------------------------------------------------------
  subroutine di_algor_get_cart_from_DI(nic,binv,qq)
    !=========================================================================!
    ! Gets new Cartesians from current estimate Cartesians and internals      !
    !                                                                         !
    !                       [newcartd in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) nic (in) : dimension of internal space                               !
    ! 2) binv (in) : inverse of the B matrix                                  !
    ! 2) qq (in) : new internals                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) coords_DI                                                            !
    ! 2) number_atoms                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) coords_cart                                                          !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none

    integer,intent(in)                                    :: nic
    real(kind=di_dp),dimension(:,:),intent(in)            :: binv
    real(kind=di_dp),dimension(1:nic),intent(in)          :: qq
    !-------------------------------------------------------------------------!
    integer :: i,j,k,ind
    real(kind=di_dp) :: qdif
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(binv,1)<3*number_atoms .or. size(binv,2)<nic) &
      & call deloc_utils_io_abort('Error in di_algor_get_cart_from_DI: Array BINV is too small')
#endif
    do i=1,nic
       qdif = coords_DI(i) - qq(i)
       ind = 1
       do j=1,number_atoms
          do k=1,3
             coords_cart(k,j) = coords_cart(k,j) + binv(ind,i)*qdif
             ind = ind + 1
          enddo
       enddo
    enddo

    return
  end subroutine di_algor_get_cart_from_DI

!--------------------------------------------------------------------
  subroutine di_algor_cart_hessian_from_DI(ndeg,ns,bm,hess_local)
    !=========================================================================!
    ! Transform Hessian in internal coordinates to cartesians                 !
    !                                                                         !
    !                       [hintbkp in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndeg (in) : dimension of H matrices                                  !
    ! 2) ns (in) : number of primitives                                       !
    ! 2) bm (in) : B matrix for the transformation (3*number_atoms,ndeg)      !
    ! 4) hess_local (out) : transformed matrix (3*number_atoms,3*number_atoms)!
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) number_atoms                                                         !
    ! 2) di_on_root,di_iprint,di_stdout                                       !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                         :: ndeg
    integer,intent(in)                                         :: ns
    real(kind=di_dp),dimension(:,:),intent(in)                    :: bm
    real(kind=di_dp),dimension(:,:),intent(out)                   :: hess_local
    !-------------------------------------------------------------------------!
    integer :: nat3,i,j,k,ierr
    real(kind=di_dp) :: ss
    real(kind=di_dp),dimension(:,:),allocatable :: bp
    !-------------------------------------------------------------------------!
    nat3 = 3*number_atoms

#ifdef debug
    if (size(bm,1)<nat3 .or. size(bm,2)<ndeg) call deloc_utils_io_abort('Error in hinbkp: array BM is too small') 
    if (size(hess_local,1)<nat3 .or. size(hess_local,2)<nat3) call deloc_utils_io_abort('Error in hinbkp: array HESS is too small') 
#endif

    allocate(bp(1:nat3,1:ndeg),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating bp in di_algor_cart_hessian_from_DI')

    hess_local = zero 

    ! expand hint matrix 
    do i=1,ns   
       do j=1,ns    
          hess_local(i,j) = hess_DI(i,j)
       enddo
    enddo

   ! transform back hess_DI to hess
    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,*) ' h in internal coordinates is:',ndeg,ns
       call di_algor_print_matrix(nat3,nat3,nat3,hess_local)
    endif

    ! make Bp=Hint*B(T)
    do i=1,nat3
       do j=1,ndeg
          ss=zero
          do k=1,ndeg
             ss=ss+hess_local(j,k)*bm(i,k)
          enddo
          bp(i,j)=ss
       enddo
    enddo

    ! make Hcart = B * Hint * B(T)
    do i=1,nat3
       do j=1,i
          ss=zero
          do k=1,ndeg
             ss=ss+bm(i,k)*bp(j,k)
          enddo
          hess_local(i,j)=ss
          hess_local(j,i)=hess_local(i,j)
       enddo
    enddo

    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,*) ' h in cartesian coordinates is:',nat3
       call di_algor_print_matrix(nat3,nat3,nat3,hess_local)
    endif

    deallocate(bp,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating bp in di_algor_cart_hessian_from_DI')

    return
  end subroutine di_algor_cart_hessian_from_DI
!------------------------------------------------------
  subroutine di_algor_make_U_matrix(ired,ut,h,eigv,b,indexb)
    !=========================================================================!
    ! Generate U matrix                                                       !
    !                                                                         !
    !                           [makeu in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ired (in) : number of reduced dimensions                             !
    ! 2) ut (out) : U  matrix                                                 !
    ! 3) h (in) : hessian eigenvectors                                        !
    ! 4) eigv (in) : hessian eigenvalues                                      !
    ! 5) b (in) : B matrix                                                    !
    ! 6) indexb (in) : index of B matrix elements                             !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) tp_eigvz                                                             !
    ! 2) di_on_root,di_iprint,di_stdout                                       !
    ! 3) num_primitive_bonds                                                  !
    ! 4) num_primitive_angles                                                 !
    ! 5) num_primitive_torsions                                               !
    ! 6) num_disconn_fragm                                                    !
    ! 7) bmat_disconn_fragm                                                   !
    ! 8) num_fixed_coords                                                     !
    ! 9) map_fixed_coords                                                     !
    !10) ndeg                                                                 !
    !11) intcor                                                               !
    !12) number_atoms                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                    :: ired
    real(kind=di_dp),dimension(:,:),intent(out)           :: ut
    real(kind=di_dp),dimension(:,:),intent(in)            :: h 
    real(kind=di_dp),dimension(:),intent(in)              :: eigv
    real(kind=di_dp),dimension(:,:),intent(in)            :: b
    integer,dimension(:,:),intent(in)                     :: indexb
    !-------------------------------------------------------------------------!
    integer :: i,j,k,in,intc1,intc2,inn,intc3,ierr,nat3
    real(kind=di_dp) :: thrsh,ss,s
    real(kind=di_dp),dimension(:),allocatable             :: v
    real(kind=di_dp),dimension(:,:),allocatable           :: vm
    !-------------------------------------------------------------------------!
    nat3 = 3 * number_atoms
#ifdef debug
    if (size(h,1)<nat3 .or. size(h,2)<nat3) &
      & call deloc_utils_io_abort('Error in di_algor_make_U_matrix: Array H is too small')
    if (size(eigv,1)<nat3) &
      & call deloc_utils_io_abort('Error in di_algor_make_U_matrix: Array RIGV is too small')
#endif

    allocate(v(1:nat3),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating v in di_algor_make_U_matrix')
    v = zero

    allocate(vm(1:nat3,1:nat3),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating vm in di_algor_make_U_matrix')
    vm=zero

    thrsh = tp_eigvz

    do i=1,ndeg
       if (abs(eigv(i+ired))>thrsh) then
          do j=1,nat3
             vm(j,i) = one/sqrt(eigv(i+ired))*h(j,i+ired)
          enddo
       endif
    enddo

    do j=1,ndeg

       ! primitives
       intc1 = num_primitive_bonds+num_primitive_angles+num_primitive_torsions
       do in = 1,intc1
          ! unroll b into v

          v = zero

          do i = 1,12
             if (indexb(i,in)/=0) v(indexb(i,in)) = v(indexb(i,in)) + b(i,in)
          enddo

          ss=zero   
          do k=1,nat3
             ss=ss+ v(k)*vm(k,j)
          enddo
          ut(in,j) = ss
       enddo

       intc2 = intc1 + num_disconn_fragm
       ! free   fragments
       if (num_disconn_fragm/=0) then
          do in = intc1+1,intc2  
             ss=zero
             inn = in - intc1
             do k=1,nat3
                ss=ss + bmat_disconn_fragm(k,inn)*vm(k,j)
             enddo
             ut(in,j) = ss
          enddo
       endif

       ! frozen atoms
       if (num_fixed_coords/=0) then
          intc3 = intc2 + num_fixed_coords 
          do in = intc2+1,intc3
             !  unroll b into v
             v = zero

             inn = in - intc2
             v(map_fixed_coords(inn)) = one

             ss=zero
             do k=1,nat3
                ss=ss+ v(k)*vm(k,j)
             enddo
             ut(in,j) = ss
          enddo
       endif


    enddo ! ndeg

    ! U is normalized
    do i=1,ndeg
       s=zero
       do k=1,intcor
          s=s + ut(k,i)*ut(k,i)
       enddo
       ss = s-one
       if (abs(ss)>1.0e-6_di_dp .and. di_on_root)  write(di_stdout,*) 'diag u ',i,s
    enddo


    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,'(/,6x,a)')'UT Matrix'
       call di_algor_print_matrix(ndeg,intcor,ndeg, ut)
       write(di_stdout,'(a)')'END of UT Matrix'
    endif

    deallocate(v,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating v in di_algor_make_U_matrix')
    deallocate(vm,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating v in di_algor_make_U_matrix')

    return
  end subroutine di_algor_make_U_matrix

!------------------------------------------------------------------------------
  subroutine di_algor_check_vector_length(n,d)
    !=========================================================================!
    ! checks the length of vector d                                           !
    ! if |d| > dmax, scales d down accordingly                                !
    !                                                                         !
    !                              [chkdp in DMol]                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) n (in) : dimension of D vector                                       !
    ! 2) d (inout) : displacement vector                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) dmax                                                                 !
    ! 2) di_on_root,di_iprint,di_stdout                                       !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                            :: n
    real(kind=di_dp),dimension(1:n)               :: d
    !-------------------------------------------------------------------------!
    real(kind=di_dp) :: dd,skal
    real(kind=di_dp),external :: ddot
    !-------------------------------------------------------------------------!

    dd = sqrt(ddot(n,d(1),1,d(1),1))

    if (dd>dmax) then
       skal = dmax/dd
       if (di_iprint>2 .and. di_on_root) write(di_stdout,1000) skal
       call dscal(n,skal,d(1),1)
       dd = dmax
    endif

    if (di_iprint>2 .and. di_on_root) write(di_stdout,1010) dd

    return

 1000 format(' calculated step too large.  step scaled by ',f9.6)
 1010 format(' step taken.  stepsize is ',f9.6)

  end subroutine di_algor_check_vector_length

!------------------------------------------------------------------------------

  subroutine di_algor_update_hessian(n,d,f,oldf,hess_local)
    !=========================================================================!
    ! This routine is responsible for updating the Hessian                    !
    ! during a geometry optimization                                          !
    !                                                                         !
    !  Several updating procedures are available:                             !
    !   (a) The Murtagh-Sargent update  [IUPDAT=4]                            !
    !       A simple rank-one update formula which allows                     !
    !       the Hessian eigenvalue structure to change                        !
    !       (see Comp.J. 13 (1970) 185)                                       !
    !   (b) The Powell update   [IUPDAT=1]                                    !
    !       This is a flexible, general rank-two update which allows          !
    !       the Hessian eigenvalue structure to change.                       !
    !       (see Math.Prog. 1 (1971) 26)                                      !
    !   (c) A Composite Powell/Murtagh-Sargent update  [IUPDAT=5]             ! 
    !       Suggested by Bofill as an improved update in practice             !
    !       (see J.Comp.Chem. 15 (1994) 1)                                    !
    !   (d) The BFGS update  [IUPDAT=2]                                       !
    !       This update is more likely to retain positive definiteness.       !
    !       Default update for minimization                                   !
    !   (e) BFGS update with positive definite check [IUPDAT=3]               !
    !                  (skip update if threatened)                            !
    !                                                                         !
    !                        [updhesp in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) n (in) : dimension of Hessian matrices                               !
    ! 2) d (in) : geometry displacement on previous cycle                     !
    ! 3) f (in) : current gradient                                            !
    ! 4) oldf (in) : previous gradient                                        !
    ! 5) hess_local (inout) : hessian                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) iupdat                                                               !
    ! 2) di_on_root,di_iprint,di_stdout                                       !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                  :: n
    real(kind=di_dp),dimension(1:n),intent(in)          :: d
    real(kind=di_dp),dimension(1:n),intent(in)          :: f
    real(kind=di_dp),dimension(1:n),intent(in)          :: oldf
    real(kind=di_dp),dimension(:,:),intent(inout)       :: hess_local
    !-------------------------------------------------------------------------!
    real(kind=di_dp),dimension(:),allocatable :: v1,v2
    real(kind=di_dp), external :: ddot
    real(kind=di_dp) :: dt,dd,dtdd,zp,zm,tmp,tt,tmm
    real(kind=di_dp),parameter :: tollzero=1.0e-8_di_dp
    integer :: i,j,ierr
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(hess_local,1)<n .or. size(hess_local,2)<n) &
      & call deloc_utils_io_abort('Hessian dimension in di_algor_update_hessian is too small')
#endif

    !  Skip update if requested
    if (iupdat==0) then
       if (di_iprint>1 .and. di_on_root) write(di_stdout,1000)
       return
    endif

    !  allocate scratch arrays
    allocate(v1(1:n),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating v1 in di_algor_update_hessian')
    allocate(v2(1:n),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating v2 in di_algor_update_hessian')
    v1 = zero
    v2 = zero

    !  Prepare for update
    !  Form and save Hess*D in V1

    do j=1,n
       tmp = d(j)
       do i=1,n
          v1(i) = v1(i) + hess_local(i,j)*tmp
       enddo
    enddo

    if (iupdat==4) then

       !  (a) Murtagh-Sargent update

       if (di_iprint>1 .and. di_on_root) write(di_stdout,1010)

       do i=1,n
          v2(i) = f(i) - oldf(i) - v1(i)
       enddo

       dt = ddot(n,v2(1),1,d(1),1)

       if (abs(dt)<tollzero) then
          if (di_iprint>1 .and. di_on_root) write(di_stdout,1015)
          go to 999
       endif

       do i=1,n
          do j=1,i
             hess_local(i,j) = hess_local(i,j) + v2(i)*v2(j)/dt
             hess_local(j,i) = hess_local(i,j)
          enddo
       enddo

    elseif (iupdat==1) then

       !  (b) Powell update

       if (di_iprint>1 .and. di_on_root) write(di_stdout,1020)

       do i=1,n
          v2(i) = f(i) - oldf(i) - v1(i)
       enddo

       dd = ddot(n,d(1),1,d(1),1)

       if (dd<tollzero) then
          if (di_iprint>1 .and. di_on_root) write(di_stdout,1015)
          go to 999
       endif

       dt = ddot(n,v2(1),1,d(1),1)
       dtdd = dt/dd

       do i=1,n
          do j=1,i
             tmp = v2(i)*d(j) + d(i)*v2(j) - d(i)*dtdd*d(j)
             hess_local(i,j) = hess_local(i,j) + tmp/dd
             hess_local(j,i) = hess_local(i,j)
          enddo
       enddo

    elseif (iupdat==5) then

       !  (c) Powell/Murtagh-Sargent update

       if (di_iprint>1 .and. di_on_root) write(di_stdout,1030)

       do i=1,n
          v2(i) = f(i) - oldf(i) - v1(i)
       end do

       dd = ddot(n,d(1),1,d(1),1)
       dt = ddot(n,v2(1),1,d(1),1)

       if (dd<tollzero .or. abs(dt)<tollzero) then
          if (di_iprint>1 .and. di_on_root) write(di_stdout,1015)
          go to 999
       endif

       tt = ddot(n,v2(1),1,v2(1),1)
       dtdd = dt/dd

       !  find mixing factor

       zp = one - (dt*dt)/(dd*tt)
       zm = one - zp

       if (di_iprint>1 .and. di_on_root) write(di_stdout,1031) zp,zm

       do i=1,n
          do j=1,i
             tmp = v2(i)*d(j) + d(i)*v2(j) - d(i)*dtdd*d(j)
             tmm = v2(i)*v2(j)
             hess_local(i,j) = hess_local(i,j) + zp*tmp/dd + zm*tmm/dt
             hess_local(j,i) = hess_local(i,j)
          enddo
       enddo

    else

       !  (d) BFGS update

       if (di_iprint>1 .and. di_on_root) write(di_stdout,1040)

       do i=1,n
          v2(i) = f(i) - oldf(i)
       end do

       dd = ddot(n,v2(1),1,d(1),1)

       !  If DD is negative, retention of positive definiteness is not
       !  guaranteed. Print a warning and skip update if requested

       if (dd<tollzero) then
          if (di_iprint>1 .and. di_on_root) write(di_stdout,1041)

          ! skip update for IUPDAT==3
          if (iupdat==3) then
             if (di_iprint>1 .and. di_on_root) write(di_stdout,1042)
             go to 999
          endif
       endif

       dt = ddot(n,d(1),1,v1(1),1)

       if (abs(dt)<tollzero .or. abs(dd)<tollzero) then
          if (di_iprint>1 .and. di_on_root) write(di_stdout,1015)
          go to 999
       endif

       do i=1,n
          do j=1,i
             tmp = (v2(i)*v2(j))/dd - (v1(i)*v1(j))/dt
             hess_local(i,j) = hess_local(i,j) + tmp
             hess_local(j,i) = hess_local(i,j)
          enddo
       enddo

    endif

 999 continue

    deallocate(v1,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating v1 in di_algor_update_hessian')
    deallocate(v2,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating v2 in di_algor_update_hessian')

    return

 1000 format(' Hessian Update Skipped by Request')
 1010 format(' Hessian Updated using Murtagh-Sargent Update')
 1015 format('**WARNING** Small Dot Product',/,' Hessian Update Skipped this cycle')
 1020 format(' Hessian Updated using Powell Update')
 1030 format(' Hessian Updated using Powell/Murtagh-Sargent Update')
 1031 format('  Mixing factors: ',F12.6,' Powell',/18X,F12.6,' Murtagh-Sargent')
 1040 format(' Hessian Updated using BFGS Update')
 1041 format('**WARNING** Hereditary positive definiteness endangered')
 1042 format(' Hessian Update Skipped this cycle')

  end subroutine di_algor_update_hessian
!------------------------------------------------------------------------------

  subroutine di_algor_make_F_fixed_cart(nat3,nic,intcor,ired,h,hh,eigv,ictyp,rcon, &
     &               xint,ktyp,savtor_loc,ierr)
    !=========================================================================!
    ! make F matrix with fixed cartesian coordinates                          !
    !                                                                         !
    !                       [cartfrz in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) nat3 (in) : dimension of H matrices (3*number_atoms)                 !
    ! 2) nic  (in) : upper limit for the number of primitives                 !
    ! 3) intcor (in) : number of primitives                                   !
    ! 4) ired (out) : number of reduced dimensions                            !
    ! 5) h (inout) : F matrix  eigenvectors                                   !
    ! 6) hh (inout) : F matrix itself                                         !
    ! 7) eigv (out) : F matrix eigenvalues                                    !
    ! 8) ictyp (inout) : which primitive internal is fixed                    !
    ! 9) rcon (out) : the fixed value (=/= 0)                                 !
    !10) ktyp (out) : type of internal (bond, angle, etc.)                    !
    !11) xint (inout) : primitive coordinates                                 !
    !12) savtor_loc (out) : saved torsions                                    !
    !13) ierr (out)                                                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) num_fixed_coords                                                     !
    ! 2) number_atoms                                                         !
    ! 3) list_fixed_coords                                                    !
    ! 4) di_on_root,di_iprint,di_stdout                                       !
    ! 5) ncycle                                                               !
    ! 6) coords_cart                                                                   !
    ! 7) ncons                                                                !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) map_fixed_coords                                                     !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: ndg>=ndeg, ndg>=ndeg+iunsat                       !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                         :: nat3
    integer,intent(in)                                         :: nic
    integer,intent(in)                                         :: intcor
    integer,intent(out)                                        :: ired
    real(kind=di_dp),dimension(1:nat3,1:nat3),intent(inout)    :: h
    real(kind=di_dp),dimension(1:nat3,1:nat3),intent(inout)    :: hh
    real(kind=di_dp),dimension(1:nat3),intent(out)             :: eigv
    integer,dimension(:,:),intent(inout)                       :: ictyp
    real(kind=di_dp),dimension(:),intent(out)                  :: rcon
    integer,dimension(:),intent(out)                           :: ktyp
    real(kind=di_dp),dimension(:),intent(inout)                :: xint
    real(kind=di_dp),dimension(:),intent(out)                  :: savtor_loc
    integer,intent(out)                                        :: ierr
    !-------------------------------------------------------------------------!
    real(kind=di_dp),dimension(:,:),allocatable  :: vm
    integer :: k,kk,ii,i,ni,j,ind
    real(kind=di_dp) :: s,ss
    character (len=80) :: string
    !-------------------------------------------------------------------------!
    allocate(vm(1:nat3,1:nat3),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating vm in di_algor_make_F_fixed_cart')
    vm = zero

    ! fixed coordinates section ; augment b
    ! form F again, diagonalize

    ! fill ictyp, rcon
    do k=1,num_fixed_coords     
       kk = k + ncons
       ictyp(kk,2) =  k+intcor
       ictyp(kk,1) = ictyp(kk,2)
       rcon(kk) = one
    enddo

    ! map list_fixed_coords into map_fixed_coords
    ii= 0
    ind = 1
    do j=1,number_atoms
       do i=1,3
          if (list_fixed_coords(i,j)==1) then
             ii=ii+1
             map_fixed_coords(ii) = ind
          endif
          ind = ind + 1
       enddo
    enddo
    if (di_iprint>4 .and. di_on_root) write(di_stdout,'(/,a)') '  FROZEN COORDINATES '

    ! fill xint
    do k=1,num_fixed_coords
       ii=0
       do j=1,number_atoms
          do i=1,3
             ii=ii+1
             if (map_fixed_coords(k)==ii) s = coords_cart(i,j)
          enddo
       enddo
       xint(k+intcor) = s

       if (di_iprint>4 .and. di_on_root) write(di_stdout,'(a,i3,f10.6)')  'coordinate: ',map_fixed_coords(k),s
    enddo

    if (ncycle<=0.and.di_iprint>1) write(di_stdout,'(/,a,i5,/)') ' Number of primitive internals increased by ',num_fixed_coords

    do k=1,num_fixed_coords
       ! zero b
       do i=1,nat3
          vm(i,k) =  zero 
       enddo
       ! expand b
       vm(map_fixed_coords(k),k) = one   
    enddo

    ! redefine intcor
    ni = intcor + num_fixed_coords

    if (ni>nic) then
       write (string,'(3i6)') intcor,ni,nic
       call deloc_utils_io_abort(' Error in di_algor_make_F_fixed_cart: ni + num_fixed_coords is too large ' // string)
    endif

    do k=1,num_fixed_coords
       savtor_loc(k+intcor) = zero
       ktyp(k+intcor) = 0
    enddo
 
    ! restore h
    h = hh

    do i=1,nat3
       do j=1,i
          ss=zero
          do k=1,num_fixed_coords
             ss=ss + vm(i,k)*vm(j,k)
          enddo
          h(i,j) = h(i,j) + ss
          h(j,i) = h(i,j)
       enddo
    enddo

    deallocate(vm,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating vm in di_algor_make_F_fixed_cart')

    ! save h
    hh = h

    if (di_iprint>1 .and. di_on_root) write(di_stdout,'(/,a,2i5)') 'Diagonalize F: fixed coordinates ',nat3,num_fixed_coords
    call di_algor_diagonalize_F_matrix (nat3,h,eigv,ired,ierr )
    if (ierr/=0 .and. di_on_root) write(di_stdout,*) ' problem in di_algor_make_F_fixed_cart'
    
    return
  end subroutine di_algor_make_F_fixed_cart

!------------------------------------------------------------------------------

  subroutine di_algor_expand_hessian_1 (ndeg,ndg,iunsat,xrcon,grcon,hold,hnew,hpad,ltypus,gnew)
    !=========================================================================!
    ! this is valid only for the first entry to di_algor_update_hessian       !
    !  copy hold into hnew                                                    !
    ! pad hnew with hpad                                                      !
    ! add derivatives of lagrangians                                          !
    ! zero for lagrangian value, xrcon                                        !
    !  ndg = ndeg + iunsat                                                    !
    ! fill xnew and gnew with xrcon, grcon                                    !
    !                                                                         !
    !                          [exphes1 in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndeg (in) : dimension of HOLD                                        !
    ! 2) ndg  (in) : 1st dimension of HPAD, dimension of HNEW                 !
    ! 3) iunsat (in) : number of unsatisfied constraints (2nd dim of HPAD)    !
    ! 4) xrcon (out) : constrained coordinates (set to zero)                  !
    ! 5) grcon (in) : gradients on constrained coordinates                    !
    ! 6) hold (in) : old hamiltonian, ndeg x ndeg                             !
    ! 7) hnew (out) : new combined hamiltonian, ndg x ndg                     !
    ! 8) hpad (out) : padded hamiltonian, ndg x iunsat                        !
    ! 9) ltypus (in) : types of constraints (bonds, angles, etc.)             !
    !10) gnew (out) : new gradients (copied from grcon)                       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: ndg>=ndeg, ndg>=ndeg+iunsat                       !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                     :: ndeg
    integer,intent(in)                                     :: ndg
    integer,intent(in)                                     :: iunsat
    real(kind=di_dp),dimension(:),intent(out)              :: xrcon
    real(kind=di_dp),dimension(:)                          :: grcon
    real(kind=di_dp),dimension(1:ndeg,1:ndeg)              :: hold
    real(kind=di_dp),dimension(1:ndg,1:ndg),intent(out)    :: hnew
    real(kind=di_dp),dimension(1:ndg,1:iunsat)             :: hpad
    integer,dimension(1:iunsat)                            :: ltypus
    real(kind=di_dp),dimension(:)                          :: gnew
    !-------------------------------------------------------------------------!
    integer :: i,j,ndegi
    real(kind=di_dp) :: cnvrt
    !-------------------------------------------------------------------------!
    !   check the dimensions
    ndegi = ndeg + iunsat
#ifdef debug
    if (ndeg>ndg) call deloc_utils_io_abort('Error in di_algor_expand_hessian_1: Dimension NDG is less than NDEG')
    if (ndegi>ndg) call deloc_utils_io_abort('Error in di_algor_expand_hessian_1: Dimension NDG is less than NDEG+IUNSAT')
#endif

    ! zero Hnew
    hnew = zero

    do i=1,ndeg
       do j=1,ndeg
          hnew(i,j) = hold(i,j)
       enddo
    enddo

    ! fill  hpad
    hpad = zero

    do i=1,iunsat

       select case (ltypus(i))
       case (1)
          cnvrt = tp_f_bond
       case (3)
          cnvrt = tp_f_angle
       case (5)
          cnvrt = tp_f_tors
       case (6) 
          cnvrt = tp_f_crt 
       case default
          cnvrt=-1 ! anyway we are dying....
          call deloc_utils_io_abort('Error in di_algor_expand_hessian_1: unknown type of internal constraint')
       end select

       Hnew(ndeg+i,ndeg+i) = cnvrt
       Hpad(ndeg+i,i) = cnvrt
    enddo
    
    do i=1,iunsat
       ! initialize Lagrangian coordinate
       xrcon(i) = zero

       ! add Lagrangian gradients
       gnew(ndegi+i) = grcon(i)
    enddo

    return

  end subroutine di_algor_expand_hessian_1

!---------------------------------------------------------------------------
  subroutine di_algor_expand_hessian_2 (ndeg,ndg,iunsat,xrcon,grcon,hold,hnew,hpad,xnew,gnew )
    !=========================================================================!
    !  this is valid for the next (ihess=1) entry to di_algor_update_hessian  !
    !  copy hold into hnew                                                    !
    ! pad hnew with hpad                                                      !
    ! add derivatives of lagrangians                                          !
    !  ndg = ndeg + iunsat                                                    !
    ! fill xnew and gnew with xrcon, grcon                                    !
    !                                                                         !
    !                          [exphes2 in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndeg (in) : dimension of HOLD                                        !
    ! 2) ndg  (in) : 1st dimension of HPAD, dimension of HNEW                 !
    ! 3) iunsat (in) : number of unsatisfied constraints (2nd dim of HPAD)    !
    ! 4) xrcon (in) : constrained coordinates                                 !
    ! 5) grcon (in) : gradients on constrained coordinates                    !
    ! 6) hold (in) : old hamiltonian, ndeg x ndeg                             !
    ! 7) hnew (out) : new combined hamiltonian, ndg x ndg                     !
    ! 8) hpad (in) : padded hamiltonian, ndg x iunsat                         !
    ! 9) xnew (out) : new coordinates                                         !
    !10) gnew (out) : new gradients                                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: ndg>=ndeg, ndg>=ndeg+iunsat                       !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                     :: ndeg
    integer,intent(in)                                     :: ndg
    integer,intent(in)                                     :: iunsat
    real(kind=di_dp),dimension(:),intent(in)               :: xrcon
    real(kind=di_dp),dimension(:),intent(in)               :: grcon
    real(kind=di_dp),dimension(1:ndeg,1:ndeg),intent(in)   :: hold
    real(kind=di_dp),dimension(1:ndg,1:ndg),intent(out)    :: hnew
    real(kind=di_dp),dimension(1:ndg,1:iunsat),intent(in)  :: hpad
    real(kind=di_dp),dimension(:),intent(out)              :: xnew
    real(kind=di_dp),dimension(:),intent(out)              :: gnew
    !-------------------------------------------------------------------------!
    integer :: i,j,ndegi
    !-------------------------------------------------------------------------!
    !   check the dimensions
    ndegi = ndeg + iunsat
#ifdef debug
    if (ndeg>ndg) call deloc_utils_io_abort('Error in di_algor_expand_hessian_2: Dimension NDG is less than NDEG')
    if (ndegi>ndg) call deloc_utils_io_abort('Error in di_algor_expand_hessian_2: Dimension NDG is less than NDEG+IUNSAT')
#endif
    ! zero Hnew
    hnew = zero

    do i=1,ndeg
       do j=1,ndeg
          hnew(i,j) = hold(i,j)
       enddo
    enddo

    ! fill  with hpad
    do i=1,iunsat
       do j=1,ndegi 
          hnew(j,ndeg+i) = hpad(j,i)
          hnew(ndeg+i,j) = hpad(j,i)
       enddo
    enddo

    ! add Lagrangian gradients and coordinates
    do i=1,iunsat
       gnew(ndegi+i) = grcon(i) ! gradients
       xnew(ndegi+i) = xrcon(i) ! coordinates
    enddo

    return
  end subroutine di_algor_expand_hessian_2

!---------------------------------------------------------------------------
  subroutine di_algor_expand_hessian_3 (ndeg,ndg,ns1,iunsat,hpad,hold,hnew)
    !=========================================================================!
    !  this is valid for the next (ihess=1) entry to di_algor_update_hessian  !
    !  copy hold into hnew                                                    !
    !                                                                         !
    !                          [exphes3 in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndeg (in) : 1st dimension of HPAD                                    !
    ! 2) ndg  (in) : dimension of HNEW                                        !
    ! 3) ns1 (in) : actual dimension of the space (hess_DI)                   !
    ! 4) iunsat (in) : number of unsatisfied constraints (2nd dim of HPAD)    !
    ! 5) hpad (in) : padded hamiltonian, ndeg x iunsat                        !
    ! 6) hold (in) : old hamiltonian, ndeg x ndeg                             !
    ! 7) hnew (out) : new combined hamiltonian, ndg x ndg                     !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) hess_DI                                                              !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: ndg<ns1, ndeg<ns1                                 !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                     :: ndeg
    integer,intent(in)                                     :: ndg
    integer,intent(in)                                     :: ns1
    integer,intent(in)                                     :: iunsat
    real(kind=di_dp),dimension(1:ndeg,1:iunsat),intent(out):: hpad
    real(kind=di_dp),dimension(1:ndeg,1:ndeg),intent(in)   :: hold
    real(kind=di_dp),dimension(1:ndg,1:ndg),intent(out)    :: hnew
    !-------------------------------------------------------------------------!
    integer :: i,j,ni,nii,iq
    !-------------------------------------------------------------------------!
#ifdef debug
    !   check the dimensions
    if (ns>ndeg) call deloc_utils_io_abort('Error in di_algor_expand_hessian_3: Dimension NDEG is less than NS')
    if (ndeg>ndg) call deloc_utils_io_abort('Error in di_algor_expand_hessian_3: Dimension NDG is less than NDEG')
    if (ns1>size(hess_DI,1)) call deloc_utils_io_abort('Error in di_algor_expand_hessian_4: ns1>size(hess_DI,1)')
#endif

    ! zero hnew
    hnew = zero

    do i=1,ndeg
       do j=1,ndeg
          hnew(i,j) = hold(i,j)
       enddo
    enddo

    ni = ndg - iunsat
    nii = ni - iunsat
    do iq = 1,iunsat
       hnew(nii+iq, ni+iq) = one
       hnew(ni+iq, nii+iq) = one
    enddo

    ! save hpdat
    do i=1,iunsat
       hpad(:,i) = hold(:,ns1+i) 
    enddo

    ! save hess
    do i=1,ns1
       do j=1,ns1
          hess_DI(i,j) = hold(i,j)
       enddo
    enddo

    return
  end subroutine di_algor_expand_hessian_3

!---------------------------------------------------------------------------
  subroutine di_algor_expand_hessian_4(ndeg,ndg,ns1,iunsat,hpad,hnew)
    !=========================================================================!
    !  this is valid for the copying of scalled Hessian from HESS AND HPAD    !
    !  into a new Hessian that will be augmented with derivatives of          !
    !  constraints and diagonalized                                           !
    !                                                                         !
    !                          [exphes4 in DMol]                              !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ndeg (in) : 1st dimension of HPAD                                    !
    ! 2) ndg  (in) : dimension of HNEW                                        !
    ! 3) ns1 (in) : actual dimension of the space (hess_DI)                   !
    ! 4) iunsat (in) : number of unsatisfied constraints (2nd dim of HPAD)    !
    ! 5) hpad (in) : padded hamiltonian, ndeg x iunsat                        !
    ! 6) hnew (out) : new combined hamiltonian, ndg x ndg                     !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) hess_DI                                                              !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: ndg<ns1, ndeg<ns1                                 !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                     :: ndeg
    integer,intent(in)                                     :: ndg
    integer,intent(in)                                     :: ns1
    integer,intent(in)                                     :: iunsat
    real(kind=di_dp),dimension(1:ndeg,1:iunsat),intent(in) :: hpad
    real(kind=di_dp),dimension(1:ndg,1:ndg),intent(out)    :: hnew
    !-------------------------------------------------------------------------!
    integer :: i,j,ni,nii,iq
    !-------------------------------------------------------------------------!

#ifdef debug
    !   check the dimensions
    if (ns1>ndg) call deloc_utils_io_abort('Error in di_algor_expand_hessian_4: Dimension NDG is less than NS1')
    if (ndeg>ndg) call deloc_utils_io_abort('Error in di_algor_expand_hessian_4: Dimension NDG is less than NDEG')
    if (ns1>size(hess_DI,1)) call deloc_utils_io_abort('Error in di_algor_expand_hessian_4: ns1>size(hess_DI,1)')
#endif

    ! zero Hnew
    hnew = zero

    do i=1,ns1   
       do j=1,ns1   
          hnew(i,j) = hess_DI(i,j)
       enddo
    enddo

    ! fill  with hpdat
    do i=1,iunsat
       do j=1,ndeg 
          hnew(j,ns1+i) = hpad(j,i)
          hnew(ns1+i,j) = hpad(j,i)
       enddo
    enddo

    ! add lagrangian gradients
    ni = ndg - iunsat
    nii = ni - iunsat

    !vmilman Let's hope these indices never get negative ...

    do iq = 1,iunsat
       hnew(nii+iq, ni+iq) = one
       hnew(ni+iq, nii+iq) = one
    enddo

    return
  end subroutine di_algor_expand_hessian_4

!---------------------------------------------------------------------------

  subroutine di_algor_fixed_atoms_symmetry
    !=========================================================================!
    !  Checks whether fixed atom constraints will preserve the                !
    !  molecular symmetry.  If not, prints out the additional                 !
    !  constraints required and Exits                                         !
    !                                                                         !
    !                          [symfix in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments: none                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) ntrans                                                               !
    ! 2) number_atoms                                                         !
    ! 3) list_fixed_coords                                                    !
    ! 4) neqatm                                                               !
    ! 5) di_on_root,di_iprint,di_stdout                                       !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified: none                                  !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    !-------------------------------------------------------------------------!
    character(len=1),dimension(1:3) :: cfix
    data cfix/'X','Y','Z'/
    integer :: i,j,it,iop,ierr
    !-------------------------------------------------------------------------!

    if (ntrans==1) return     ! there is no symmetry to conserve
    ierr = 0

    do i=1,number_atoms
       do j=1,3
          if (list_fixed_coords(j,i)==1) then

             !  we have a fixed coordinate
             !  check for equivalent atoms
             !  these must also be fixed for symmetry
             !  to be maintained

             do iop=2,ntrans
                it = neqatm(i,iop)
                if (it/=i .and. list_fixed_coords(j,it)/=1) then
                   ierr = -1
                   if (di_on_root) write(di_stdout,1000) it,cfix(j)
                endif
             enddo

          endif
       enddo
    enddo

    if (ierr/=0) then
       if (di_on_root) then
          write (di_stdout,'(a)')'Error: Only fixing certain atomic coordinates will break symmetry'
          write (di_stdout,'(a)')'Error: Please turn off symmetry and start the calculation again'
       endif
       call deloc_utils_io_abort('Error: Only fixing certain atomic coordinates will break symmetry')
    endif

    return

 1000 format(' To Preserve Symmetry  Atom ',i3,2x,a1,' Coordinate Should be Fixed')

  end subroutine di_algor_fixed_atoms_symmetry
!---------------------------------------------------------------------------

  subroutine di_algor_constr_from_rigid_body(num_primitives,lidelp,jctyp,xrcon,xprim_loc)
    !=========================================================================!
    !  Identify num_fixed_constr_rigid_body constraints belonging to          !
    !  num_rigid_body_atoms atoms from the list list_rigid_body_atoms.        !
    !  Constraints are identified from the list of all internals, lidelp.     !
    !                                                                         !
    !                        [fixatcon in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) num_primitives (in) - number of primitives                           !
    ! 2) lidelp (in) - list of internals                                      !
    ! 3) jctyp (inout,optional) - list of internals constrained               !
    ! 4) xrcon (inout,optional) - value of internal constrained               !
    ! 5) xprim_loc (in,optional) - values of primitive internals              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) list_rigid_body_atoms                                                !
    ! 2) num_rigid_body_atoms                                                 !
    ! 3) num_primitive_bonds                                                  !
    ! 4) num_primitive_angles                                                 !
    ! 4) ncons                                                                !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) num_fixed_constr_rigid_body                                          !
    ! 2) ncons                                                                !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! mode - determines whether only num_fixed_constr_rigid_body is found (0) !
    !        or test if the constraints found are on the list of original     !
    !        constraints (1), "mode" is determined by the presence of all     !
    !        three optional arguments                                         !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                    :: num_primitives
    integer,dimension(:,:),intent(in)                     :: lidelp
    integer,dimension(:,:),intent(inout),optional         :: jctyp
    real(kind=di_dp),dimension(:),intent(inout),optional  :: xrcon
    real(kind=di_dp),dimension(:),intent(in),optional     :: xprim_loc
    !-------------------------------------------------------------------------!
    integer :: intn,intn0,in,ika,ikb,ikc,kc,ia1,ia2,ia3,ia,ib1,ib2,ib3,jb
    integer :: ila,ilb,ilc,ld,iss,isave,i
    integer :: mode
    !-------------------------------------------------------------------------!

    mode = 0 
    if (present(jctyp) .and. present(xrcon) .and. present(xprim_loc)) mode = 1

    ! loop over all primitives
    intn = 0
    intn0 = 0
    do in = 1,num_primitives

       call deloc_utils_unpack_atom_indices(lidelp(1,in), ika, ikb, ikc, kc)
       call deloc_utils_unpack_atom_indices(lidelp(2,in), ia1, ia2, ia3, ia)
       call deloc_utils_unpack_atom_indices(lidelp(3,in), ib1, ib2, ib3, jb)
       call deloc_utils_unpack_atom_indices(lidelp(4,in), ila, ilb, ilc, ld)

       ! is the atoms kc & ia & jb & ld on a list of atoms list_rigid_body_atoms

       call di_algor_is_atom_in_list(kc,num_rigid_body_atoms,list_rigid_body_atoms,isave)
       iss = isave
       if (iss==0) call di_algor_is_atom_in_list(ia,num_rigid_body_atoms,list_rigid_body_atoms,isave)
       iss = iss + isave
       if (iss==0 .and. in>num_primitive_bonds) call di_algor_is_atom_in_list(jb,num_rigid_body_atoms,list_rigid_body_atoms,isave)
       iss = iss + isave
       if (iss==0 .and. in>(num_primitive_angles+num_primitive_bonds)) call di_algor_is_atom_in_list(ld,num_rigid_body_atoms, &
                                                                        & list_rigid_body_atoms,isave)
       iss = iss + isave

       ! survive; this internal will be optimized iss/=0

       if (iss==0) then
          ! constrained: this internal will be constrained iss=0
          intn0 = intn0 + 1

          if (mode==1)  then
             ! check if this in constrained internal is already present on the jctyp list
             isave = 1
             if (ncons/=0) call di_algor_is_atom_in_list(in,ncons,jctyp(1:ncons,1),isave)

             if (isave/=0) then
                intn  = intn + 1

                ! update ncons
                ncons = ncons + 1

                jctyp(ncons,1) = in
                jctyp(ncons,2) = in
                xrcon(ncons) = xprim_loc(in)
             endif
          endif ! mode
       endif ! iss
    enddo ! in = 1,num_primitives
  
    num_fixed_constr_rigid_body = intn
    if (mode==0) num_fixed_constr_rigid_body = intn0

    if (di_iprint>4 .and. mode>0 .and. di_on_root) then
       write(di_stdout,*) 'di_algor_constr_from_rigid_body: ncons,mfxcon ',ncons,num_fixed_constr_rigid_body
       do i=1,ncons 
          write(di_stdout,*) 'jct ',i,jctyp(i,1)
       enddo
    endif

    return
  end subroutine di_algor_constr_from_rigid_body
!------------------------------------------------------------------------------

  subroutine di_algor_is_atom_in_list(klist_atom,mixcar,mixci,isave)
    !=========================================================================!
    !  Search among the mixcar atoms in mixci looking for klist_atom          !
    !                                                                         !
    !                          [klistcar in DMol]                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) klist_atom (in) - atom ID to look for                                !
    ! 2) mixcar (in) - number of elements of mixci to search through          !
    ! 3) mixci (in) - list of atom IDs                                        !
    ! 4) isave (out) - 0 if atom ID is found, 1 otherwise                     !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                     :: klist_atom
    integer,intent(in)                     :: mixcar
    integer,dimension(:),intent(in)        :: mixci
    integer,intent(out)                    :: isave
    !-------------------------------------------------------------------------!
    integer :: i
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(mixci)<mixcar) call deloc_utils_io_abort('Error in di_algor_is_atom_in_list: mixcar is too big')
#endif

    isave = 0
    do i=1,mixcar
       if (mixci(i)==klist_atom) return
    enddo

    isave = 1

    return
    
  end subroutine di_algor_is_atom_in_list

!------------------------------------------------------------------------------

  subroutine di_algor_eigenvector_following(nc,n,ncon,neg,negreq,mode,u,eigval,f,d,status)
    !=========================================================================!
    !  Driving Routine for Eigenvector (mode) Following algorithm             !
    !                                                                         !
    !  Based on:                                                              !
    !   "On Finding Transition States"                                        !
    !    Cerjan & Miller, J.Chem.Phys. 75 (1981) 2800                         !
    !                                                                         !
    !   "Walking on Potential Energy Surfaces"                                !
    !    Simons, Jorgenson, Taylor & Ozment, J.Phys.Chem. 87 (1983) 2745      !
    !                                                                         !
    !   "Searching for Stationary Points on Surfaces"                         !
    !    Banerjee, Adams, Simons & Shepard, J.Phys.Chem. 89 (1985) 52         ! 
    !                                                                         !
    !  For full details see:                                                  !
    !                                                                         !
    !   "An Algorithm for the Location of Transition States"                  !
    !    J.Baker, J.Comp.Chem. 7 (1986) 385                                   !
    !                                                                         !
    !                                                                         !
    !                            [optefp in DMol]                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) nc (in) - number of modes used to make the step                      !
    ! 2) n (in) - number of primitives (dimension of the problem)             !
    ! 3) ncon (in) - number of constraints                                    !
    ! 4) neg (in) - dnumber of negative eigenvalues of Hessian                !
    ! 5) negreq (in) - desired number of negative eigenvalues                 !
    ! 6) mode (inout) - mode being followed (zero if mode following is off)   !
    ! 7) u (in) - Hessian eigenvectors                                        !
    ! 8) eigval (inout) - Hessian eigenvalues                                 !
    ! 9) f (in) - gradient                                                    !
    !10) d (out) - new displacement vector (a step)                           !
    !11) ierr (out)                                                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root,di_iprint,di_stdout                                       !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: fx:gradient in the basis of Hessian eigenvectors!
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                           :: nc
    integer,intent(in)                           :: n
    integer,intent(in)                           :: ncon
    integer,intent(in)                           :: neg
    integer,intent(in)                           :: negreq
    integer,intent(inout)                        :: mode
    real(kind=di_dp),dimension(:,:),intent(in)   :: u
    real(kind=di_dp),dimension(:),intent(inout)  :: eigval
    real(kind=di_dp),dimension(:),intent(in)     :: f
    real(kind=di_dp),dimension(:),intent(out)    :: d
    integer,intent(out)                          :: status
    !-------------------------------------------------------------------------!
    integer :: i,ierr
    real(kind=di_dp), external :: ddot
    real(kind=di_dp),dimension(:),allocatable   :: fx
    !-------------------------------------------------------------------------!
    status = 0

    allocate(fx(1:nc),stat=ierr) 
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating fx in di_algor_eigenvector_following')

    !  Use the selected NC eigenvectors of the Hessian (stored in U)
    !  to form the next step

    !  First form the FX vector
    !  (the components of F along the local Hessian modes)

    do i=1,nc
       fx(i) = ddot(n,u(1,i),1,f(1),1)
    enddo

    if (ncon==0) then

       !  Normal Unconstrained Optimization

       !  Take the P-RFO step for a TS search
       !  Take the simple RFO step for a minimum search

       if (di_iprint>1 .and. di_on_root) then
          if (neg/=negreq) write(di_stdout,1000)
          if (negreq==0) write(di_stdout,1010)
          if (negreq==1) write(di_stdout,1020)
       endif

       call di_algor_ef_RFOstep(nc,n,negreq,mode,u,eigval,fx,d,status)

    else

       !  Constrained optimization in Cartesian coordinates

       !  Take the P-RFO step

       if (di_iprint>1 .and. di_on_root) then
          if (neg/=negreq) write(di_stdout,1000)
          if (negreq==ncon) write(di_stdout,1030)
          if (negreq==ncon+1) write(di_stdout,1020)
       endif

       call di_algor_ef_RFOstep_constr(nc,n,ncon,negreq,mode,u,eigval,fx,d,status)

    endif

    deallocate(fx,stat=ierr) 
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating fx in di_algor_eigenvector_following')

    return

 1000 format(' Warning: Hessian does not have the Desired Local Structure')
 1010 format(/,' Minimum Search - Taking Simple RFO Step')
 1020 format(/,' Transition State Search - Taking P-RFO Step')
 1030 format(/,' Minimum Search - Taking P-RFO Step')

  end subroutine di_algor_eigenvector_following
!------------------------------------------------------------------------------
  subroutine di_algor_ef_check_H_scale_unsat(n,n1,iunsat,hpad,eigval,scale)
    !=========================================================================!
    !  Checks whether it is likely to be profitable to scale                  !
    !  the Hessian in an attempt to calculate a suitable                      !
    !  shift parameter for EF during a constrained optimization               !
    !  (The original value for the shift parameter was unacceptable)          !
    !                                                                         !
    !  The assumption here is that the lowest (positive) Hessian              !
    !  eigenvalue scales linearly whilst the (initially too high)             !
    !  shift parameter scales at a slower rate; thus it may be                !
    !  possible to find a suitable shift parameter (i.e. one lower            !
    !  than the Hessian eigenvalue) by scaling the Hessian.                   !
    !                                                                         !
    !                            [chkskald in DMol]                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) n (in) - dimension of Hessian                                        !
    ! 2) n1 (in) - auxiliary dimension                                        !
    ! 3) iunsat (in) - number of unsatisfied constraints                      !
    ! 4) hpad (inout) - padded strip of a Hessian                             !
    ! 5) eigval (in) - EigVal(1)=lowest "positive" Hessian Eigenvalue         !
    !                  Eigval(2)  -  previously calculated shift parameter    !
    !                  * NOTE *  These values set in CONFormD                 !
    ! 6) scale (out) - contains Hessian scale factor on exit                  !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) hess_DI                                                              !
    ! 2) di_on_root,di_iprint,di_stdout                                       !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                      :: n
    integer,intent(in)                                      :: n1
    integer,intent(in)                                      :: iunsat
    real(kind=di_dp),dimension(1:n,1:iunsat),intent(inout)  :: hpad
    real(kind=di_dp),dimension(1:2),intent(in)              :: eigval
    real(kind=di_dp),intent(inout)                          :: scale

!------------------------------------------------------------------------------
    real(kind=di_dp),parameter :: eigmin=0.001_di_dp
    real(kind=di_dp),parameter :: scalmax=10.0_di_dp
    real(kind=di_dp) :: sc
!------------------------------------------------------------------------------

    if (scale/=one) then

       !  already made an attempt at scaling which failed
       !  restore original Hessian and exit

       call dscal(n1*n1,one/scale,hess_DI(1,1),1)
       call dscal(n*iunsat,one/scale,hpad(1,1),1)
       scale = one

    else

       !  calculate likely scale factor

       if (eigval(1)<eigmin .or. eigval(2)<zero) go to 95

       scale = eigval(2)/eigval(1)
       scale = scale+scale
       if (scale<=scalmax) then
          if (di_on_root) write(di_stdout,1000) scale
          call dscal(n1*n1,scale,hess_DI(1,1),1)
          call dscal(n*iunsat,scale,hpad(1,1),1)
       else
          scale = one
          go to 95
       endif

    endif

    return

 95 continue

    sc = 0.5_di_dp
    if (eigval(1)/=zero) sc = eigval(2)/eigval(1)
    if (di_on_root) write(di_stdout,1010) sc*2.0_di_dp
    if (di_iprint>4 .and. di_on_root) write(di_stdout,'(a,3f9.6)') &
     & ' scaling hessian : e1,e2,scale ',eigval(1),eigval(2),scale

    return

 1000 format(' ** Scaling Hessian by ',F9.6,' **',/)
 1010 format(' ** It is not possible to scale  Hessian **',f9.6,/)

  end subroutine di_algor_ef_check_H_scale_unsat

!------------------------------------------------------------------------------

  subroutine di_algor_ef_check_H_scale(n,eigval,scale)
    !=========================================================================!
    !  Checks whether it is likely to be profitable to scale                  !
    !  the Hessian in an attempt to calculate a suitable                      !
    !  shift parameter for EF during a constrained optimization               !
    !  (The original value for the shift parameter was unacceptable)          !
    !                                                                         !
    !  The assumption here is that the lowest (positive) Hessian              !
    !  eigenvalue scales linearly whilst the (initially too high)             !
    !  shift parameter scales at a slower rate; thus it may be                !
    !  possible to find a suitable shift parameter (i.e. one lower            !
    !  than the Hessian eigenvalue) by scaling the Hessian.                   !
    !                                                                         !
    !                            [chkskalp in DMol]                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) n (in) - dimension of Hessian                                        !
    ! 2) eigval (in) - EigVal(1)=lowest "positive" Hessian Eigenvalue         !
    !                  Eigval(2)  -  previously calculated shift parameter    !
    !                  * NOTE *  These values set in CONFormD                 !
    ! 3) scale (out) - contains Hessian scale factor on exit                  !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) hess_DI                                                              !
    ! 2) di_on_root,di_iprint,di_stdout                                       !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                         :: n
    real(kind=di_dp),dimension(1:2),intent(in)                 :: eigval
    real(kind=di_dp),intent(inout)                             :: scale
    !-------------------------------------------------------------------------!
    real(kind=di_dp),parameter :: eigmin=0.001_di_dp
    real(kind=di_dp),parameter :: scalmax=5.0_di_dp
    real(kind=di_dp) :: sc
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(hess_DI,1)<n .or. size(hess_DI,2)<n) call deloc_utils_io_abort('Error in di_algor_ef_check_H_scale: N is too small')
#endif

    if (scale/=one) then

       !  already made an attempt at scaling which failed
       !  restore original hessian and exit

       call dscal(n*n,one/scale,hess_DI(1,1),1)
       scale = one

    else

       !  calculate likely scale factor

       if (eigval(1)<eigmin .or. eigval(2)<zero) go to 95

       scale = eigval(2)/eigval(1)
       scale = scale+scale

       if (scale<=scalmax) then
          if (di_on_root) write(di_stdout,1000) scale
          call dscal(n*n,scale,hess_DI(1,1),1)
       else
          scale = one
          go to 95
       endif

    endif

    return

 95 continue

    sc = 0.5_di_dp
    if (eigval(1)/=zero) sc = eigval(2)/eigval(1)
    if (di_on_root) write(di_stdout,1010) sc*2.0_di_dp
    if (di_iprint>4 .and. di_on_root) write(di_stdout,'(a,3f9.6)') &
     & ' scaling hessian : e1,e2,scale ',eigval(1),eigval(2),scale

    return

 1000 format(' ** Scaling Hessian by ',F9.6,' **',/)
 1010 format(' ** It is not possible to scale  Hessian **',f9.6,/)

  end subroutine di_algor_ef_check_H_scale
!------------------------------------------------------------------------------
  subroutine di_algor_ef_RFOstep(nc,n,negreq,mode,u,eigval,fx,d,ierr)
    !=========================================================================!
    ! This routine works out the RFO step (Rational Function Approximation)   !
    ! MIN Search:  Forms a step by simple RFO) that attempts to               !
    !              minimize along all Hessian modes                           !
    !  TS Search:  Forms a step by P-RFO (partitioned RFO) that attempts to   !
    !              maximize along the direction of a chosen Hessian mode      !
    !              (stored in VMODE) and minimize along all other modes       !
    !                                                                         !
    !                            [formdp in DMol]                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) nc (in) - number of modes used to make the step                      !
    ! 2) n (in) - number of primitives (dimension of the problem)             !
    ! 3) negreq (in) - desired number of negative eigenvalues                 !
    ! 4) mode (inout) - mode being followed                                   !
    ! 5) u (in) - Hessian eigenvectors                                        !
    ! 6) eigval (inout) - Hessian eigenvalues                                 !
    ! 7) fx (in) - gradient in Hessian eigenvector basis                      !
    ! 8) d (out) - new displacement vector (a step)                           !
    ! 9) ierr (out)                                                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root,di_iprint,di_stdout                                       !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                          :: nc
    integer,intent(in)                                          :: n
    integer,intent(in)                                          :: negreq
    integer,intent(inout)                                       :: mode
    real(kind=di_dp),dimension(:,:),intent(in)                  :: u
    real(kind=di_dp),dimension(1:nc),intent(inout)              :: eigval
    real(kind=di_dp),dimension(1:nc),intent(in)                 :: fx
    real(kind=di_dp),dimension(1:n),intent(out)                 :: d
    integer,intent(out)                                         :: ierr
    !-------------------------------------------------------------------------!
    real(kind=di_dp) :: lamda,lamda0,xleft,xright,temp
    real(kind=di_dp),parameter :: toll=1.0e-8_di_dp
    real(kind=di_dp),parameter :: step=0.05_di_dp
    !real(kind=di_dp),parameter :: big=1.0e+3_di_dp
    integer,parameter :: maxit=999
    integer :: numit,it,nmode,jt,i,j
    !-------------------------------------------------------------------------!
    ierr = 0

    numit = 0
    it = 0

    if (negreq==0) go to 10

    !  (a)  Maximize along selected Hessian mode

    if (mode/=0) then

       call di_algor_find_mode_to_follow(n,mode,nmode,u)

       !  On return from di_algor_find_mode_to_follow, Nmode is the new mode along which
       !  the energy is to be maximized

       if (di_iprint>1 .and. nmode/=mode .and. di_on_root) write (di_stdout,1000) mode,nmode
       mode = nmode
       if (di_iprint>1 .and. di_on_root) write(di_stdout,1010) mode
       it = mode

       !  If the mode now being followed is the lowest mode
       !  switch off mode following

       if (mode==1) then
          mode = 0
          if (di_iprint>1 .and. di_on_root) write(di_stdout,1020)
       endif

    else ! mode/=0

       if (di_iprint>1 .and. di_on_root) write(di_stdout,1030)
       it = 1

    endif ! mode/=0

    !  Calculate the Shift

    lamda0 = eigval(it) + sqrt( eigval(it)**2 + 4.0_di_dp*fx(it)*fx(it) )
    lamda0 = 0.5_di_dp*lamda0

    if (di_iprint>2 .and. di_on_root) write(di_stdout,1040) lamda0
    if (nc==1) go to 40

    !  (b)  Minimize along all other Hessian modes

 10 continue

    jt = 1+it
    if (jt>2) jt = 1

    if (di_iprint>1 .and. di_on_root) then
       if (negreq==0) write(di_stdout,1050)
       if (negreq==1) write(di_stdout,1060)
    endif

    !  Solve Iteratively for Lamda
    !  Initial guess for Lamda is zero EXCEPT Lamda < EigVal(JT)

    lamda = zero
    if (eigval(jt)<zero) then
       lamda  = eigval(jt)-step
    endif

 15 continue

    temp = zero
    do i=1,nc
       if (i/=it) then
          temp = temp + ( fx(i)*fx(i) ) / ( lamda - eigval(i) )
       endif
    enddo

    if (temp<lamda) then
       xleft  = temp
       xright = lamda
    else
       if (temp>eigval(jt)) then
          lamda = 0.5_di_dp*(lamda + eigval(jt))
          goto 15
       endif
       xleft  = lamda
       xright = temp
    endif

    lamda = 0.5_di_dp*(xleft+xright)

 20 continue

    numit = numit+1
    temp = zero
    do i=1,nc
       if (i/=it) then
          temp = temp + ( fx(i)*fx(i) ) / ( lamda - eigval(i) )
       endif
    enddo

    if (di_iprint>3 .and. di_on_root) write(di_stdout,1111) lamda,temp

    !  Check for Convergence of Lamda

    if (abs(lamda-temp)<toll) go to 30

    !  Check for Maximum Iterations Exceeded

    if (numit>maxit) go to 90

    !  Cautious Bracketing Scheme (C.Koelmel, 1/9/96)

    if (temp<lamda) then
       if (temp>xleft ) xleft = temp
       xright = lamda
    else
       xleft  = lamda
       if (temp<xright) xright = temp
    endif

    lamda = 0.5_di_dp*(xleft+xright)

    go to 20

    !  At this point we should have an acceptable value for Lamda
    !  Make final check

 30 continue

    if (di_iprint>2 .and. di_on_root) write(di_stdout,1040) lamda
    if (lamda>eigval(jt)) go to 91
    if (eigval(jt)>zero .and. lamda>zero) go to 91

    !  Everything OK
    !  Calculate the new Step

 40 continue

    d = zero
      
    do i=1,nc
       temp = fx(i)/( lamda-eigval(i) )
       if (i==it) temp = fx(i)/( lamda0-eigval(i) )
       do j=1,n
          d(j) = d(j) + temp*u(j,i)
       enddo
    enddo

    return

    !    ** ERROR SECTION **

 90 continue
    if (di_on_root) write(di_stdout,1070)
    ierr = 1
    return

 91 continue
    if (di_on_root) write(di_stdout,1080)
    ierr = 1
    return

 1000 format('**WARNING** Mode Switching: Was Following mode ',I3,' Now Following mode ',I3)
 1010 format(' Searching for Lamda that Maximizes along mode ',I3)
 1020 format(' Mode Following Switched Off')
 1030 format(' Searching for Lamda that Maximizes Along the Lowest mode')
 1040 format(' Value Taken    Lamda = ',F12.8)
 1050 format(' Searching for Lamda that Minimizes Along All modes')
 1060 format(' Searching for Lamda that Minimizes Along All other modes')
 1070 format(//,' Warning: ',/,' ** unable to determine Lamda in FormD **',/)
 1080 format(//,' Warning: ',/,' ** Error in determining Lamda in FormD **',/)
 1111 format(' in iterative cycle:  Lamda = ',f12.8,' TEMP = ',f12.8)

  end subroutine di_algor_ef_RFOstep
!------------------------------------------------------------------------------
  subroutine di_algor_ef_RFOstep_constr(nc,n,ncons,negreq,mode,u,eigval,fx,d,ierr )
    !=========================================================================!
    ! This routine works out the RFO step (Rational Function Approximation)   !
    ! MIN Search:  Forms a step by P-RFO (partitioned RFO) that attempts to   !
    !              maximize along the NCons constraint modes and minimize     !
    !              along all other Hessian modes                              !
    !  TS Search:  Forms a step by P-RFO that attempts to maximize            !
    !              along the NCons constraint modes plus one other            !
    !              chosen mode (stored in VMDEL) and minimize                 !
    !              along all other Hessian modes                              !
    !                                                                         !
    !                         [conformdp in DMol]                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) nc (in) - number of modes used to make the step                      !
    ! 2) n (in) - number of primitives (dimension of the problem)             !
    ! 3) ncons (in) - number of internal constraints                          !
    ! 4) negreq (in) - desired number of negative eigenvalues                 !
    ! 5) mode (inout) - mode being followed                                   !
    ! 6) fx (in) - gradient in Hessian eigenvector basis                      !
    ! 7) u (in) - Hessian eigenvectors                                        !
    ! 8) eigval (inout) - Hessian eigenvalues                                 !
    ! 9) d (out) - new displacement vector (a step)                           !
    !10) ierr (out)                                                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root,di_iprint,di_stdout                                       !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                          :: nc
    integer,intent(in)                                          :: n
    integer,intent(in)                                          :: ncons
    integer,intent(in)                                          :: negreq
    integer,intent(inout)                                       :: mode
    real(kind=di_dp),dimension(:,:),intent(in)                  :: u
    real(kind=di_dp),dimension(1:nc),intent(inout)              :: eigval
    real(kind=di_dp),dimension(1:nc),intent(in)                 :: fx
    real(kind=di_dp),dimension(1:n),intent(out)                 :: d
    integer,intent(out)                                         :: ierr
    !-------------------------------------------------------------------------!
    real(kind=di_dp) :: lamda,lamda1,lamda2,temp
    real(kind=di_dp),parameter :: toll=1.0e-8_di_dp
    real(kind=di_dp),parameter :: step=0.05_di_dp
    real(kind=di_dp),parameter :: big=1.0e+3_di_dp
    integer,parameter :: maxit=999
    integer :: numit,it,i,j,jt,nmode
    !-------------------------------------------------------------------------!
    ierr = -1

! Both lambdas are used in a proper way. this initialisation is just to calm down ftncheck: (21/Apr/2006 aperlov)

    lamda1 = -big 
    lamda2 = -big 

    numit = 0
    it = ncons
    d = zero

    if (ncons==0) go to 50

    ! (a)  Maximize along NegReq Hessian modes

    if (mode/=0) then

       call di_algor_find_mode_to_follow(n,mode,nmode,u)

       ! On return from di_algor_find_mode_to_follow, Nmode is the additional mode
       ! along which the energy is to be maximized

       if (di_iprint>1 .and. nmode/=mode .and. di_on_root) write(di_stdout,1000) mode,nmode
       mode = nmode
       if (di_iprint>1 .and. di_on_root) write(di_stdout,1010) mode
       it = mode

       ! If the mode now being followed is now amongst the NCons+1
       ! lowest modes, switch off mode following

       if (mode<=negreq) then
          mode = 0
          if (di_iprint>1 .and. di_on_root) write(di_stdout,1020)
       endif

    else ! mode /=0

       if (di_iprint>1 .and. di_on_root) then
          if (negreq==ncons) write(di_stdout,1030)
          if (negreq>ncons)  write(di_stdout,1035)
       endif
       if (negreq>ncons) it = negreq

    endif ! mode /=0

    if (negreq==1) then

       ! There is only one constraint and we are minimizing

       lamda = eigval(it) + sqrt( eigval(it)**2 + 4.0_di_dp*fx(it)*fx(it) )
       lamda = 0.5_di_dp*lamda

       if (di_iprint>2 .and. di_on_root) write(di_stdout,1040) lamda

       if (eigval(it+1)>zero .and. lamda>eigval(it+1)) then

          ! Set up for potential Hessian scaling in subroutine ChkSKAL
          eigval(1) = eigval(it+1)
          eigval(2) = lamda
          go to 91
       endif

    else ! negreq==1

       !  We are maximizing along the lowest NCons modes
       !  and possibly mode IT
       !
       !  Solve iteratively for Lamda
       !  Initial guess for Lamda is zero EXCEPT Lamda > EigVal(IT)

       lamda = zero
       if (eigval(it)>zero) lamda = eigval(it) + step

 10    continue
       numit = numit+1
       temp = zero

       do i=1,ncons
          temp = temp + ( fx(i)*fx(i) )/( lamda-eigval(i) )
       enddo

       if (it>ncons) temp = temp + ( fx(it)*fx(it) )/( lamda-eigval(it) )
       if (di_iprint>3 .and. di_on_root) write(di_stdout,1111) lamda,temp

       ! Check for Convergence of Lamda

       if (abs(lamda-temp)<toll) go to 20

       !  Check for Maximum Iterations Exceeded

       if (numit>maxit) go to 90

       lamda = temp
       go to 10

 20    continue

       ! At this point we should have an acceptable value for Lamda
       ! Make final check

       if (di_iprint>2 .and. di_on_root) write(di_stdout,1040) lamda
       if (lamda<zero) go to 91
       if (eigval(it+1)>zero .and. lamda>eigval(it+1)) then

       ! Set up for potential Hessian scaling in subroutine ChkSKAL
          eigval(1) = eigval(it+1)
          eigval(2) = lamda
          go to 91
       endif

    endif ! negreq==1

    ! Calculate the Maximization Step

    do i=1,ncons
       temp = fx(i)/( lamda-eigval(i) )
       do j=1,n
          d(j) = d(j) + temp*u(j,i)
       enddo
    enddo

    if (it>ncons) then
       temp = fx(it)/( lamda-eigval(it) )
       do j=1,n
          d(j) = d(j) + temp*u(j,it)
       enddo
    endif

    numit = 0

 50 continue

    jt = 1+it
    if (jt>ncons+2) jt = ncons+1

    ! (b)  Minimize along all other Hessian modes

    if (di_iprint>1 .and. di_on_root) then
       if (ncons==0) write(di_stdout,1050)
       if (ncons>0)  write(di_stdout,1060)
    endif

    !  Solve Iteratively for Lamda
    !  Initial guess for Lamda is zero EXCEPT Lamda < EigVal(JT)

    lamda = zero
    if (eigval(jt)<zero) then
       lamda  = eigval(jt) - step
       lamda1 = eigval(jt)
       lamda2 = -big
    endif

 60 continue
    numit = numit+1
    temp = zero
    do i=ncons+1,nc
       if (i/=it) temp = temp + ( fx(i)*fx(i) )/( lamda-eigval(i) )
    enddo
    if (di_iprint>3 .and. di_on_root) write(di_stdout,1111) lamda,temp

    ! Check for Convergence of Lamda

    if (abs(lamda-temp)<toll) go to 70

    ! Check for Maximum Iterations Exceeded

    if (numit>maxit) go to 90

    if (eigval(jt)<zero) then

       ! (ii) Cautious Bracketing Scheme

       if (temp<lamda) lamda1 = lamda
       if (temp>lamda) lamda2 = lamda
       if (lamda2>-big)  lamda = 0.5_di_dp*(lamda1+lamda2)
       if (lamda2==-big) lamda = lamda-step

       go to 60

    else

       ! (i)  Simple Iterative Scheme

       lamda = temp
       go to 60

    endif

    ! At this point we should have an acceptable value for Lamda
    ! Make final check

 70 continue

    if (di_iprint>2 .and. di_on_root) write(di_stdout,1040) lamda
    if (lamda>eigval(jt)) go to 91
    if (eigval(jt)>zero .and. lamda>zero) go to 91

    !  Calculate the Minimization Step

    do i=ncons+1,nc
       if (i/=it) then
          temp = fx(i)/( lamda-eigval(i) )
          do j=1,n
             d(j) = d(j) + temp*u(j,i)
          enddo
       endif
    enddo

    ierr = 0
    return

    !    ** ERROR SECTION **

 90 continue
    if (di_on_root) write (di_stdout,1070)
    return

 91 continue
    if (di_on_root) write(di_stdout,1080)
    return

 1000 format('**WARNING** Mode Switching: Was Following mode ',I3,' Now Following mode ',I3)
 1010 format(' Searching for Lamda that Maximizes along mode ',I3)
 1020 format(' Mode Following Switched Off')
 1030 format(' Searching for Lamda that Maximizes Along the Constraint modes Only')
 1035 format(' Searching for Lamda that Maximizes Along the Lowest Non-Constraint mode')
 1040 format(' Value Taken    Lamda = ',F12.8)
 1050 format(' Searching for Lamda that Minimizes Along All modes')
 1060 format(' Searching for Lamda that Minimizes Along All other modes')
 1070 format(/,' Warning: Unable to determine Lamda IN CONFormD ',/)
 1080 format(/,' Warning: Error in determining Lamda IN CONFormD ',/,' Scaling of the Hessian will be attempted',/)

 1111 format(' in iterative cycle:  Lamda = ',f12.8,' TEMP = ',f12.8)

  end subroutine di_algor_ef_RFOstep_constr

!------------------------------------------------------------------------------
  subroutine di_algor_update_constraints (jctyp,xcon)
    !=========================================================================!
    ! update jctyp, xcon with icc, rcon                                       !
    ! This implementation simply copies them over                             !
    !                                                                         !
    !                         [upcons in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) icc (in) - index of constraints                                      !
    ! 2) rcon (in) - values of constraints                                    !
    ! 3) jctyp (out) - index of constraints                                   !
    ! 4) xcon (out) - values of constraints                                   !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) ncons                                                                !
    ! 2) icc                                                                  !
    ! 3) rcon                                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,dimension(:,:),intent(inout)                    :: jctyp
    real(kind=di_dp),dimension(:),intent(inout)             :: xcon
    !-------------------------------------------------------------------------!
    integer :: i
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(jctyp,1)<ncons .or. size(jctyp,2)<2) &
    & call deloc_utils_io_abort('Error in di_algor_update_constraints: Array JCTYP is too small')
    if (size(icc,1)<ncons .or. size(icc,2)<2) &
    & call deloc_utils_io_abort('Error in di_algor_update_constraints: Array ICC is too small')
    if (size(rcon,1)<ncons .or. size(rcon,2)<1) &
    & call deloc_utils_io_abort('Error in di_algor_update_constraints: Array RCON is too small')
    if (size(xcon,1)<ncons) &
    & call deloc_utils_io_abort('Error in di_algor_update_constraints: Array XCON is too small')
#endif
    ! this version does not do much, does it?
    do i=1,ncons
       jctyp(i,1) = icc(i,18)
       jctyp(i,2) = icc(i,19)
       xcon(i) = rcon(i,1)
    enddo
      
    return
  end subroutine di_algor_update_constraints
!------------------------------------------------------------------------------

  subroutine di_algor_constr_to_primitives(nic,lidelp,ni1,num_prim_coord,icc,almostlinear,ipath)
    !=========================================================================!
    ! Make the final list of bonds, angles, dihedrals                         !
    ! by considering internal constraints:                                    !
    ! 1) check if they are on the list lidelp                                 !
    ! 2) if not add to the list                                               !
    !                                                                         !
    !                     [addcons in DMol]                                   !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) nic (in) - maximum dimension of the primitive set                    !
    ! 2) lidelp (inout) - store cell indices and atom numbers                 !
    ! 3) ni1 (in) - start checking the list of primitives from this position  !
    ! 4) num_prim_coord (inout) - number of primitives                        !
    ! 5) icc (inout) - index of constraints                                   !
    ! 6) almostlinear (in) - threshold to determine linear angles             !
    ! 7) ipath (in) - type of the internal (1,3,5)                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root,di_iprint,di_stdout                                       !
    ! 2) ncons                                                                !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none

    integer,intent(in)                                       :: nic
    integer,dimension(:,:),intent(inout)                     :: lidelp
    integer,intent(in)                                       :: ni1
    integer,intent(inout)                                    :: num_prim_coord
    integer,dimension(:,:),intent(inout)                     :: icc
    real(kind=di_dp),intent(in)                              :: almostlinear
    integer,intent(in)                                       :: ipath
    !-------------------------------------------------------------------------!
    logical :: skip
    integer :: j,list,itra,itrb,itrc,jb,i0,j0,k0,ityp,jaf,ita,itb,itc,jbf,isa,isb,isc
    integer :: jcf,iwa,iwb,iwc,jdf,ixa,ixb,ixc,ifound,ni,ia,ika,ikb,ikc,jc,kc
    integer :: ija,ijb,ijc,isk,ibad,ila,ilb,ilc,ld,jj
    real(kind=di_dp) :: qq
    character(len=10) :: string
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(icc,1)<ncons .or. size(icc,2)<21) &
    & call deloc_utils_io_abort('Error in di_algor_constr_to_primitives: Array ICC is too small')
    if (size(lidelp,1)<4 .or. size(lidelp,2)<num_prim_coord) &
    & call deloc_utils_io_abort('Error in di_algor_constr_to_primitives: Array LIDELP is too small')
#endif

    list = num_prim_coord

    do j=1,ncons

       skip = .false.
       ityp = icc(j,17)

       jaf = icc(j,1)
       ita = icc(j,2)
       itb = icc(j,3)
       itc = icc(j,4)

       jbf = icc(j,5)
       isa = icc(j,6)
       isb = icc(j,7)
       isc = icc(j,8)

       jcf = icc(j,9)
       iwa = icc(j,10)
       iwb = icc(j,11)
       iwc = icc(j,12)

       jdf = icc(j,13)
       ixa = icc(j,14)
       ixb = icc(j,15)
       ixc = icc(j,16)

       if (ityp==1 .and. ityp==ipath) then

          ! bonds
          ifound = 0
          do ni = ni1,num_prim_coord
             call deloc_utils_unpack_atom_indices(lidelp(1,ni),  i0,  j0,  k0, ia)
             call deloc_utils_unpack_atom_indices(lidelp(2,ni),itra,itrb,itrc, jb)

             if (ia==jaf .and. jb==jbf) then
                if (itra==isa .and. itrb==isb .and. itrc==isc) skip=.true.
             endif
             if (ia==jbf .and. jb==jaf) then
                if (itra==-isa .and. itrb==-isb .and. itrc==-isc) skip=.true.
             endif

             if (skip) then
                ifound = ni
                exit
             endif
          enddo ! num_prim_coord

          if (.not.skip) then
             list = list +1
             ifound = list
          
             if (list>nic) then
                write (string,'(i8)') list
                call deloc_utils_io_abort('Error: too many bonds' // string)
             endif

             call deloc_utils_pack_atom_indices(lidelp(1,list),   0,   0,   0, jaf)
             call deloc_utils_pack_atom_indices(lidelp(2,list), isa, isb, isc, jbf)
          endif ! not skip

          jj = ifound
          if (ifound==0) jj = list
          icc(j,18) = jj     
          icc(j,19) = jj    

          if (di_iprint>=1 .and. di_on_root) then
             if (skip) then
                write(di_stdout,'(2(a,i5))') ' constraint # ',j,' is a primitive # ',jj
             else
                write(di_stdout,'(2(a,i5))') ' constraint # ',j,' is a new primitive # ',jj
             endif
          endif

       endif ! bonds

       if (ityp==3 .and. ityp==ipath) then

          ! angles
          ifound = 0
          do ni = ni1,num_prim_coord

             call deloc_utils_unpack_atom_indices(lidelp(1,ni), ija, ijb, ijc, jb)
             call deloc_utils_unpack_atom_indices(lidelp(2,ni),  i0,  j0,  k0, ia)
             call deloc_utils_unpack_atom_indices(lidelp(3,ni), ika, ikb, ikc, jc)

             if (ia==jbf) then

                isk = 0
                if (jb==jaf .and. ita==ija .and. itb==ijb .and. itc==ijc) isk = isk+1  
                if (jc==jcf .and. ika==iwa .and. ikb==iwb .and. ikc==iwc) isk = isk+1
                if (jc==jaf .and. ita==ika .and. itb==ikb .and. itc==ikc) isk = isk+1  
                if (jb==jcf .and. iwa==ija .and. iwb==ijb .and. iwc==ijc) isk = isk+1  

                if (isk>=2) then
                   skip =.true.
                   ifound = ni
                   exit
                endif
             endif ! ia=jbf
          enddo ! ni

          if (.not.skip) then
             list = list +1
             ifound = list
             if (list>nic) then
                write (string,'(i8)') list
                call deloc_utils_io_abort('Error: too many angles' // string)
             endif

             ! is this  a linear bend
             call di_algor_find_internal_coord(qq,  3, &
     &                     jaf,ita,itb,itc, &
     &                     jbf,isa,isb,isc, &
     &                     jcf,iwa,iwb,iwc, &
     &                     jdf,ixa,ixb,ixc, &
     &                     almostlinear, ibad)
             if (ibad==1) then
                write (string,'(i8)') j
                call deloc_utils_io_abort('Error: Constraint # '// string // ' is a linear bend')
             endif

             call deloc_utils_pack_atom_indices(lidelp(1,list), ita, itb, itc, jaf)
             call deloc_utils_pack_atom_indices(lidelp(2,list),   0,   0,   0, jbf)
             call deloc_utils_pack_atom_indices(lidelp(3,list), iwa, iwb, iwc, jcf)

          endif ! skip

          jj = ifound
          if (ifound==0) jj = list
          icc(j,18) = jj     
          icc(j,19) = jj    
          if (di_iprint>=1 .and. di_on_root) then
             if (skip) then
                write(di_stdout,'(2(a,i5))') ' constraint # ',j,' is a primitive # ',jj
             else
                write(di_stdout,'(2(a,i5))') ' constraint # ',j,' is a new primitive # ',jj
             endif
          endif

       endif ! angles

       if (ityp==5 .and. ityp==ipath) then

          ! torsions
          ifound = 0
          do ni = ni1,num_prim_coord

             call deloc_utils_unpack_atom_indices(lidelp(1,ni), ika, ikb, ikc, kc)
             call deloc_utils_unpack_atom_indices(lidelp(2,ni),  i0,  j0,  k0, ia)
             call deloc_utils_unpack_atom_indices(lidelp(3,ni),itra,itrb,itrc, jb)
             call deloc_utils_unpack_atom_indices(lidelp(4,ni), ila, ilb, ilc, ld)


             if (ia==jbf .and. jb==jcf .and. itra==iwa .and. itrb==iwb .and. itrc==iwc) then
                isk = 0
                if (kc==jaf .and. ita==ika .and. itb==ikb .and. itc==ikc) isk = isk+1  
                if (ld==jdf .and. ixa==ila .and. ixb==ilb .and. ixc==ilc) isk = isk+1

                if (isk>=2) then
                   skip =.true.
                   ifound = ni
                   exit
                endif
             endif ! ia=jbf
          enddo ! ni

          if (.not.skip) then
             list = list + 1
             ifound = list
             if (list>nic) then
                write (string,'(i8)') list
                call deloc_utils_io_abort('Error: too many torsions' // string)
             endif

             ! is this  a linear torsion
             call di_algor_find_internal_coord(qq,  5, &
     &                     jaf,ita,itb,itc, &
     &                     jbf,isa,isb,isc, &
     &                     jcf,iwa,iwb,iwc, &
     &                     jdf,ixa,ixb,ixc, &
     &                     almostlinear, ibad )
             if (ibad==1) then
                write (string,'(i8)') j
                call deloc_utils_io_abort('Error: Constraint # '// string // ' has three or more atoms that are virtually linear')
             endif

             call deloc_utils_pack_atom_indices(lidelp(1,list), ita, itb, itc, jaf)
             call deloc_utils_pack_atom_indices(lidelp(2,list),   0,   0,   0, jbf)
             call deloc_utils_pack_atom_indices(lidelp(3,list), iwa, iwb, iwc, jcf)
             call deloc_utils_pack_atom_indices(lidelp(4,list), ixa, ixb, ixc, jdf)

          endif ! skip

          jj = ifound
          if (ifound==0) jj = list
          icc(j,18) = jj     
          icc(j,19) = jj    

          if (di_iprint>=1 .and. di_on_root) then
             if (skip) then
                write(di_stdout,'(2(a,i5))') ' constraint # ',j,' is a primitive # ',jj
             else
                write(di_stdout,'(2(a,i5))') ' constraint # ',j,' is a new primitive # ',jj
             endif
          endif

       endif ! torsions
    
    enddo ! j,ncons

    num_prim_coord = list
    return
  end subroutine di_algor_constr_to_primitives

!------------------------------------------------------------------------------
  subroutine di_algor_num_negative_hess_EVs(n,eigval,num_neg)
    !=========================================================================!
    ! Find the number of negative Hessian eigenvalues                         !
    !                                                                         !
    !                     [findnegp in DMol]                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) n (in) - dimension of eigenvalues array                              !
    ! 2) eigval (in) - eigenvalues                                            !
    ! 3) num_neg (out) - number of negative eigenvalues                       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: assumed that eigenvalues are ordered              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                          :: n
    real(kind=di_dp),dimension(1:n),intent(in)  :: eigval
    integer,intent(out)                         :: num_neg
    !-------------------------------------------------------------------------!
    integer :: i
    !-------------------------------------------------------------------------!
    do i=1,n
       if (eigval(i)>zero) exit
    enddo

    num_neg = i-1

    return
  end subroutine di_algor_num_negative_hess_EVs
!------------------------------------------------------------------------------
  subroutine di_algor_cleanup_hessian(n,ncons,u,eigval,g,thrsh,nc,ierr)
    !=========================================================================!
    !  Removes eigenvectors with "zero" eigenvalues from a set of             !
    !  N Hessian modes. Additionally, if symmetry is switched on,             !
    !  removes all symmetry-breaking modes including, if applicable,          !
    !  constraint modes                                                       !
    !                                                                         !
    !                         [chkhesp in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) n (in) - number of modes to check (i.e. total dimension of problem)  !
    ! 2) ncons (in) - initial number of constraints                           !
    ! 3) u (inout) - eigenvectors stored as columns (n x n)                   !
    ! 4) eigval (inout) - eigenvalues (assumed to be ordered)                 !
    ! 5) g (in) - current gradient vector                                     !
    !             (if symflag is t all modes whose overlap with the gradient  !
    !             is less than thrsh will be rejected)                        !
    ! 6) thrsh (in) - eigenvalues below thrsh will be rejected as "zero"      !
    ! 7) nc (out)   - number of modes remaining                               !
    ! 8) ierr (out)                                                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root,di_iprint,di_stdout                                       !
    ! 2) symflag                                                              !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: assumed that eigenvalues are ordered              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                        :: n
    integer,intent(in)                                        :: ncons
    real(kind=di_dp),dimension(1:n,1:n),intent(inout)         :: u
    real(kind=di_dp),dimension(1:n),intent(inout)             :: eigval
    real(kind=di_dp),dimension(1:n),intent(in)                :: g
    real(kind=di_dp),intent(in)                               :: thrsh
    integer,intent(out)                                       :: nc
    integer,intent(out)                                       :: ierr
    !-------------------------------------------------------------------------!
    real(kind=di_dp),parameter :: ltol=1.0e-8_di_dp
    real(kind=di_dp) :: gdv
    real(kind=di_dp),external :: ddot
    integer :: nat3,imin,i,imax,it,j
    !-------------------------------------------------------------------------!
    ierr = 0

    nat3 = n-ncons

    ! first remove zeroes

    imin = 0
    do i=1,n
       if (abs(eigval(i))>thrsh) then
          if (eigval(i)<zero) then
             imin = i
          else
             go to 20
          endif
       endif
    enddo

    ! should never get here except for di- and maybe tri-atomics
    if (n>9) then
       if (di_on_root) write(di_stdout,1000)
       ierr =1
    endif

 20 continue
    imax = i-1

    if (imin>=imax) then

       ! there are no zero eigenvalues
       nc = n

    else

       ! All eigenvalues between IMIN+1 and IMAX to be deleted
       do i=1,n-imax
          eigval(imin+i) = eigval(imax+i)
          call dcopy(n,u(1,imax+i),1,u(1,imin+i),1)
       enddo

       nc = n - (imax-imin)

    endif

    ! "zero" eigenvalues removed
    ! now check for symmetry

    if (symflag) then

       ! systematically remove all modes that break symmetry
       ! (whose overlap with the gradient is less than thrsh)

       it = 0
 40    continue
       it = it+1
       if (it>nc) go to 60

       ! form scalar product with gradient

       gdv = ddot(nat3,g(1),1,u(1,it),1)
       if (abs(gdv)<thrsh) then

          ! mode has no overlap with gradient - remove
          nc = nc-1

          ! if we are doing a constrained optimization and the eigenvalue
          ! is negative, check if a constraint mode is being removed
          ! **  WARNING  **  THIS CHECK IS NOT BULLET-PROOF

          if (ncons>0 .and. eigval(it)<zero) then

             ! check if the mode has any constraint component
             do i=nat3+1,n
                if (abs(u(i,it))>ltol) &
                   call deloc_utils_io_abort('Symmetry breaking of constrained mode is not allowed')
             enddo

          endif

          ! shift all eigenvalues and vectors down

          do j=it,nc
             eigval(j) = eigval(j+1)
             call dcopy(n,u(1,j+1),1,u(1,j),1)
          enddo

          it = it-1

       endif ! abs(gdv)<thrsh

       go to 40

 60    continue

    endif ! symflag

    ! NC is the total number of modes remaining
    if (di_iprint>1 .and. di_on_root) then
       write(di_stdout,1010) nc
       write(di_stdout,1020)
       write(di_stdout,1030) (eigval(i),i=1,nc)
    endif
    if (di_iprint>5 .and. di_on_root) then
       write(di_stdout,1050)
       call di_algor_print_matrix(nc,n,nc,u)
    endif

    return

 1000 format(/,2X,' Warning: ',/,' Hessian Appears to have all zero or negative eigenvalues')
 1010 format(/,I3,' Hessian modes will be used to form the next step')
 1020 format('  Hessian Eigenvalues:')
 1030 format(1X,6F12.6)
 1050 format(/,'  Hessian Eigenvectors:')

  end subroutine di_algor_cleanup_hessian

!------------------------------------------------------------------------------

  subroutine di_algor_adjust_hess_EVs_range(n,eigmin,eigmax,eigval)
    !=========================================================================!
    ! Checks absolute magnitudes of a set of N eigenvalues.                   !
    ! All eigenvalues below EigMIN in magnitude are increased                 !
    ! to this value while all those above EigMAX are reduced                  !
    ! to this value.                                                          !
    !                                                                         !
    !                       [chkmagp in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !  1) n (in) : number of eigenvalues to search                            !
    !  2) eigmin (in) : lower boundary for eigenvalues                        !
    !  3) eigmax (in) : upper boundary for eigenvalues                        !
    !  4) eigval (inout) : Hessian eigenvalues                                !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: It is assumed that the eigenvalues are ordered.   !
    !                       (eigmin<eigmax, eigmin, eigmax>0)                 !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 28/02/2003                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                    :: n
    real(kind=di_dp),intent(in)                           :: eigmin
    real(kind=di_dp),intent(in)                           :: eigmax
    real(kind=di_dp),dimension(1:n),intent(inout)         :: eigval ! ordered
    !-------------------------------------------------------------------------!
    integer :: i
    logical :: finished
    !-------------------------------------------------------------------------!
#ifdef debug
    if (eigmin>eigmax) call deloc_utils_io_abort('Error in di_algor_adjust_hess_EVs_range: eigmin>eigmax')
    if (eigmin<zero) call deloc_utils_io_abort('Error in di_algor_adjust_hess_EVs_range: negative eigmin')
    if (eigmax<zero) call deloc_utils_io_abort('Error in di_algor_adjust_hess_EVs_range: negative eigmax')
#endif

    finished = .true.
    ! First check if magnitudes are below eigmin

    do i=1,n
       if (eigval(i)>eigmin) then
          finished = .false.
          exit
       endif
       if (abs(eigval(i))<eigmin) then
          if (eigval(i)<zero) eigval(i) = -eigmin
          if (eigval(i)>zero) eigval(i) =  eigmin
          if (di_on_root) write (di_stdout,1000) i,eigval(i)
       endif
    enddo

    ! No need to check against eigmax now
    if (finished) go to 999

    ! Now check if magnitudes exceed eigmax

    do i=1,n
       if (abs(eigval(i))<eigmax) then
          finished = .false.
          exit
       endif
       if (eigval(i)<zero) eigval(i) = -eigmax
       if (eigval(i)>zero) eigval(i) =  eigmax
       if (di_on_root) write(di_stdout,1010) i,eigval(i)
    enddo 

    if (finished) go to 999

    do i=n,1,-1
       if (eigval(i)<eigmax) return
       eigval(i) = eigmax
       if (di_on_root) write(di_stdout,1010) i,eigval(i)
    enddo

999 continue
    return

 1000 format('**WARNING** Magnitude of eigenvalue ',I3,' too small. Replaced by ',F12.6)
 1010 format('**WARNING** Magnitude of eigenvalue ',I3,' too large. Replaced by ',F12.6)

  end subroutine di_algor_adjust_hess_EVs_range
!------------------------------------------------------------------------------

  subroutine di_algor_find_mode_to_follow(ncoord,mode,nmode,u)
    !=========================================================================!
    ! This routine determines which Hessian mode to follow during             !
    ! a Transition State search.                                              ! 
    !  On the first cycle this has been set by the user.                      !
    !  Subsequently it is determined by criterion of maximum                  !
    !  overlap with the mode followed on the previous cycle.                  !
    !  If the user supplies an explicit vector rather than a mode #,          !
    !  then we flag this with mode=-1.                                        !
    ! The result is saved in the module variable VMDEL                        !
    !                                                                         !
    !                         [modeoverlap in DMol]                           !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) ncoord (in) : number of coordinates (1st dimension of U)             !
    ! 2) mode (in)  : 1st step, can tell which mode to use                    !
    ! 3) nmode (out) : mode which is followed (==mode if ncycle==1)           !
    ! 4) u (in) : matrix ncoord x ncoord, contains all the modes              !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) ncycle                                                               !
    ! 2) di_on_root,di_iprint,di_stdout                                       !
    ! 3) vmdel                                                                !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) vmdel                                                                !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                                        :: ncoord
    integer,intent(in)                                        :: mode
    integer,intent(out)                                       :: nmode
    real(kind=di_dp),dimension(:,:),intent(in)                :: u
    !-------------------------------------------------------------------------!
    integer :: it,i
    real(kind=di_dp) :: tovlp,ovlp
    real(kind=di_dp),external :: ddot
    !-------------------------------------------------------------------------!
#ifdef debug
    if (size(u,1)<ncoord .or. size(u,2)<ncoord) &
            & call deloc_utils_io_abort('Error in di_algor_find_mode_to_follow: U array is too small')
    if (ncoord>size(vmdel)) call deloc_utils_io_abort('Error in di_algor_find_mode_to_follow: NCOORD is too big')
#endif

    if (ncycle==1 .and. mode/=-1) then

       it = mode
       if (di_iprint>1 .and. di_on_root) write(di_stdout,1000) mode

    else

       it = 1
       tovlp = ddot(ncoord,u(1,1),1,vmdel(1),1)
       tovlp = abs(tovlp)

       do i=2,ncoord
          ovlp = ddot(ncoord,u(1,i),1,vmdel(1),1)
          ovlp = abs(ovlp)
          if (ovlp>tovlp) then
             tovlp = ovlp
             it = i
          endif
       enddo

       if (di_iprint>1 .and. di_on_root) write(di_stdout,1010) tovlp

    endif

    ! Store the mode to be followed in VMDEL
    do i=1,ncoord
       vmdel(i) = u(i,it)
    enddo

    nmode = it

    return

 1000 format(' Hessian Mode Following Switched On',/,' Following mode: ',I3)
 1010 format(' Overlap of current mode with previous mode is: ',F12.6)

  end subroutine di_algor_find_mode_to_follow

!------------------------------------------------------------------------------
  subroutine di_algor_add_periodic_box
    !=========================================================================!
    ! Make a funny, big box around molecule with molecule inside              !
    !                                                                         !
    !                           [p_box in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments: none                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) coords_cart                                                          !
    ! 2) lattice_dimension                                                    !
    ! 3) number_atoms                                                         !
    ! 4) elat                                                                 !
    ! 5) elatp                                                                !
    ! 6) ipvlt                                                                !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) np_p                                                                 !
    ! 2) elat                                                                 !
    ! 3) xc_redressed                                                         !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    !-------------------------------------------------------------------------!
    integer :: i
    real (kind=di_dp) :: xmin,ymin,zmin,xmax,ymax,zmax,xl,yl,zl,box
    !-------------------------------------------------------------------------!
    ! find edges of molecule  
    xmin = 10000.0_di_dp 
    ymin = 10000.0_di_dp 
    zmin = 10000.0_di_dp 
    xmax = zero 
    ymax = zero 
    zmax = zero 

    do i=1,number_atoms
       xmin = min(xmin,coords_cart(1,i))
       ymin = min(ymin,coords_cart(2,i))
       zmin = min(zmin,coords_cart(3,i))
       xmax = max(xmax,coords_cart(1,i))
       ymax = max(ymax,coords_cart(2,i))
       zmax = max(zmax,coords_cart(3,i))
    enddo
     
    ! zero np_p
    np_p = 0

    xl = xmax - xmin
    yl = ymax - ymin
    zl = zmax - zmin

    ! make funny elat
    elat = zero
    box = 40.0_di_dp
    xc_redressed = .false.
    elat(1,1) = xl + box
    elat(2,2) = yl + box
    elat(3,3) = zl + box

    return
  end subroutine di_algor_add_periodic_box


!------------------------------------------------------------------------------
  subroutine di_algor_print_matrix(n,nrow,ncol,a)
    !=========================================================================!
    ! Prints out N columns of an NRow * NCol matrix A                         !
    ! Format - six columns on a page.                                         !
    !                                                                         !
    !                         [prntmat in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) n (in) : number of columns to print                                  !
    ! 2) nrow (in) : 1st dimension of A                                       !
    ! 3) ncol (in) : 2nd dimension of A                                       !
    ! 4) a (in) : nrow x ncol matrix to be printed                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  has to be used on the root node                  !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    integer, intent(in)                                   :: nrow
    integer, intent(in)                                   :: ncol
    integer, intent(in)                                   :: n
    real (kind=di_dp),dimension(:,:),intent(in)           :: a
    !-------------------------------------------------------------------------!
    integer, parameter :: maxcol=6
    integer :: np,nt,nq,imin,imax,i,j,k,nleft
    !-------------------------------------------------------------------------!
#ifdef debug
    if (nrow>size(a,1)) call deloc_utils_io_abort('Error in di_algor_print_matrix - nrow is too big')
    if (ncol>size(a,2)) call deloc_utils_io_abort('Error in di_algor_print_matrix - ncol is too big')
#endif

    np = min(n,ncol)  ! can't print more than NCol

    nt = np/maxcol

    do i=1,nt
       imin = (i-1)*maxcol + 1
       imax = i*maxcol
       write(di_stdout,1000)
       do j=1,nrow
          write(di_stdout,1100) (a(j,k),k=imin,imax)
       end do
    end do

    nq = nt*maxcol
    nleft = np - nq
    if (nleft==0) return

    write(di_stdout,1000)
    do j=1,nrow
       write(di_stdout,1100) (a(j,k),k=nq+1,np)
    end do

    return

 1000 format(/)
 1100 format(1x,6f12.6)

  end subroutine di_algor_print_matrix

!------------------------------------------------------------------------------

  subroutine di_algor_print_constraints(lidelp,ktyp,ictyp)
    !=========================================================================!
    ! print values of constraint bonds, angles, dihedrals                     !
    !                                                                         !
    !                           [wrconsti in DMol]                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) lidelp (in) : store primitive indices                                !
    ! 2) ktyp (in) : type of primitives (0,1,3,5,6)                           !
    ! 3) ictyp (in) : pointer to constraints                                  !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) num_fixed_constr_rigid_body                                          !
    ! 2) di_on_root,di_iprint,di_stdout                                       !
    ! 3) tp_con_sat                                                           !
    ! 4) nunsat0                                                              !
    ! 5) icc                                                                  !
    ! 6) rcon                                                                 !
    ! 7) ncons                                                                !
    ! 8) tp_benlin                                                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables modified:                                       !
    ! 1) nunsat                                                               !
    ! 2) icc                                                                  !
    ! 3) rcon                                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  none                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    integer,dimension(:,:)                           :: lidelp
    integer,dimension(:)                             :: ktyp
    integer,dimension(:,:)                           :: ictyp
    !-------------------------------------------------------------------------!
    real(kind=di_dp) :: qq,vali,scale,rc
    integer :: j,ni,ita,itb,itc,ia,isa,isb,isc,jb,jc,iwa,iwb,iwc,jd,ixa,ixb,ixc,ibad
    integer :: jj,ityp
    !-------------------------------------------------------------------------!
    ! new nunsat
    nunsat = 0

    if (di_on_root) then
       write(di_stdout,'(/,a)') '  INTERNAL CONSTRAINTS'
       write(di_stdout,'(a,/)') ' #  Type     Target    Actual           Definition'
    endif

    ! bonds
    do j = 1,ncons

       ! where is ncons on a list of lidelp       
       ni = ictyp(j,2)

       ! what type is that:   1,3,5 ?
       ityp = ktyp(ni)

       select case (ityp)
       case (1) 
          ! bond
          call deloc_utils_unpack_atom_indices(lidelp(1,ni),ita,itb,itc,ia)
          call deloc_utils_unpack_atom_indices(lidelp(2,ni),isa,isb,isc,jb)
          scale = 1.0_di_dp
          call di_algor_find_internal_coord( qq,  1, &
     &                   ia,ita,itb,itc, &
     &                   jb,isa,isb,isc, &
     &                   jc,iwa,iwb,iwc, &
     &                   jd,ixa,ixb,ixc, &
     &                   tp_benlin, ibad )

          vali = qq*length_conv   ! convert to angstrom

          jj = j
          rc = rcon(jj,1)*length_conv   ! convert to angstrom

          call di_algor_print_constraint(j, 1, &
     &                   ia,ita,itb,itc, &
     &                   jb,isa,isb,isc, &
     &                   jc,iwa,iwb,iwc, &
     &                   jd,ixa,ixb,ixc, &
     &                   vali,rc, tp_con_sat)

       case (3) 
          ! angle
          call deloc_utils_unpack_atom_indices(lidelp(1,ni), ita, itb, itc, ia)
          call deloc_utils_unpack_atom_indices(lidelp(2,ni), isa, isb, isc, jb)
          call deloc_utils_unpack_atom_indices(lidelp(3,ni), iwa, iwb, iwc, jc)

          call di_algor_find_internal_coord( qq,  3, &
     &                   ia,ita,itb,itc, &
     &                   jb,isa,isb,isc, &
     &                   jc,iwa,iwb,iwc, &
     &                   jd,ixa,ixb,ixc, &
     &                   tp_benlin, ibad )

          vali = qq*toang
          jj = j
          rc = rcon(jj,1)*toang

          scale = 10.0_di_dp
          call di_algor_print_constraint(j, 3, &
     &                   ia,ita,itb,itc, &
     &                   jb,isa,isb,isc, &
     &                   jc,iwa,iwb,iwc, &
     &                   jd,ixa,ixb,ixc, &
     &                   vali,rc, tp_con_sat*scale)

       case (5)
          ! torsion
          call deloc_utils_unpack_atom_indices(lidelp(1,ni), ita, itb, itc, ia)
          call deloc_utils_unpack_atom_indices(lidelp(2,ni), isa, isb, isc, jb)
          call deloc_utils_unpack_atom_indices(lidelp(3,ni), iwa, iwb, iwc, jc)
          call deloc_utils_unpack_atom_indices(lidelp(4,ni), ixa, ixb, ixc, jd)


          call di_algor_find_internal_coord( qq,  5, &
     &                   ia,ita,itb,itc, &
     &                   jb,isa,isb,isc, &
     &                   jc,iwa,iwb,iwc, &
     &                   jd,ixa,ixb,ixc, &
     &                   tp_benlin, ibad )

          vali = qq*toang
          jj = j
          rc = rcon(jj,1)*toang

          scale = 10.0_di_dp
          call di_algor_print_constraint(j, 5, &
     &                   ia,ita,itb,itc, &
     &                   jb,isa,isb,isc, &
     &                   jc,iwa,iwb,iwc, &
     &                   jd,ixa,ixb,ixc, &
     &                   vali,rc, tp_con_sat*scale)

       case default
          scale = -1.0_di_dp ! this is to satisfy ftncheck. anyway we are exiting....
          jj=-1
          call deloc_utils_io_abort('Error in di_algor_print_constraints:  unknown type of internal constraint')

       end select

       if (abs(vali-rc)>tp_con_sat*scale) then
          icc(jj,20)=0
          nunsat = nunsat + 1
       else
          icc(jj,20)=1
       endif

       ! register type and gradient
       icc(jj,21) = jj
       rcon(jj,2) = qq - rcon(j,1)

    enddo  ! ni     

    ! order constraints to put unsatisfied at the end
    if (nunsat0/=nunsat .and. nunsat/=0) then
       call di_algor_order_unsat_constr(ncons-num_fixed_constr_rigid_body)
    endif

    if (di_on_root) write(di_stdout,'(/)')
    if (nunsat/=0 .and. di_iprint>4 .and. di_on_root) then
       write(di_stdout,'(/,a,/,a)') &
  &  ' Constrained Optimization is done in Delocalized Internals using Lagrange Multiplier Algorithm'
       write(di_stdout,'(/)')
! end unsatisfied
    endif

    return
  end subroutine di_algor_print_constraints 
!------------------------------------------------------------------------------

  subroutine di_algor_print_cartesians
    !=========================================================================!
    ! Print out Cartesian coordinates (in specified units)                    !
    !                                                                         !
    !                       [prcart in DMol]                                  !
    !-------------------------------------------------------------------------!
    ! Arguments: none                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) number_atoms                                                         !
    ! 2) coords_cart                                                          !
    ! 3) atsymb                                                               !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  has to be used on the root node                  !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    character(len=8) :: symb
    character(len=1) :: symb1
    character(len=2) :: symb2
    !-------------------------------------------------------------------------!
    integer :: iatm,lensym
    real(kind=di_dp) :: cx,cy,cz
    !-------------------------------------------------------------------------!

    write(di_stdout,1000)

    do iatm=1,number_atoms
       cx = coords_cart(1,iatm)*length_conv
       cy = coords_cart(2,iatm)*length_conv
       cz = coords_cart(3,iatm)*length_conv

       symb = adjustl(atsymb(iatm))

       lensym = di_algor_string_length(symb,len(symb))
       if (lensym==1) then

          symb1 = symb(1:1)
          write(di_stdout,1011) iatm,symb1,cx,cy,cz

       elseif (lensym==2) then

          symb2 = symb(1:2)
          write(di_stdout,1012) iatm,symb2,cx,cy,cz

       endif ! lensym

    enddo

    write(di_stdout,2000)
    return

 1000 format(10x,70('-'),/,15x,'ATOM',15x,'X',19x,'Y',19x,'Z')
 1011 format(10x,i3,4x,a,1x,3(2x,f18.6))
 1012 format(10x,i3,4x,a,3(2x,f18.6))
 2000 format(10x,70('-'),/)

  end subroutine di_algor_print_cartesians
!------------------------------------------------------------------------------


  integer function di_algor_string_length(str,l)
    !=========================================================================!
    ! find rightmost non-blank in string str of length l                      !
    ! and return its position                                                 !
    !                                                                         !
    !                         [lstr in DMol]                                  !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) str (in): a string of length l                                       !
    ! 2) l (in)  : length of the string                                       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:  none                                     !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  none                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    integer,intent(in)                :: l  
    character(len=l),intent(in)       :: str
    !-------------------------------------------------------------------------!
    integer :: i
    !-------------------------------------------------------------------------!

    di_algor_string_length=0
    i = min(len(str),l) + 1

 10 continue
    i = i-1
    if (i==0) return
    if (str(i:i)==' ' .or. str(i:i)==char(0)) go to 10
    di_algor_string_length=i

    return 
  end function di_algor_string_length
!------------------------------------------------------------------------------

  function di_algor_s2(x)
    !=========================================================================!
    ! di_algor_s2(x) = sqrt ( 1 - x*x )                                       !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) x (in)                                                               !
    ! 2) y (in)                                                               !
    !                                                                         !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) di_on_root,di_iprint,di_stdout                                       !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  none                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    real (kind=di_dp), intent(in)            :: x
    real (kind=di_dp)                        :: di_algor_s2
    !-------------------------------------------------------------------------!

    if (abs(x)>1.00001_di_dp) then
       di_algor_s2 = zero
       if (di_on_root .and. di_iprint>2) write(di_stdout,*) ' WARNING : Argument of function di_algor_s2 is ',x,' > 1 '
    elseif (abs(x)>=one) then
       di_algor_s2 = zero
    else
       di_algor_s2 = sqrt(one-x*x)
    endif

    return
  end function di_algor_s2
!------------------------------------------------------------------------------

  function di_algor_arc1 (x,y)
    !=========================================================================!
    ! computes di_algor_arc1 = arctan(y/x)                                    !
    !                                                                         !
    !                           [arc1 in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) x (in)                                                               !
    ! 2) y (in)                                                               !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) zero                                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  none                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    real (kind=di_dp), intent(in)            :: x,y
    real (kind=di_dp)                        :: di_algor_arc1
    !-------------------------------------------------------------------------!

    if (abs(x)<1.0e-11_di_dp) then
       di_algor_arc1 = di_pi/2.0_di_dp
    else
       di_algor_arc1 = atan(y/x)
       if (x<zero) di_algor_arc1 = di_algor_arc1 + di_pi
    end if

    return
  end function di_algor_arc1

!------------------------------------------------------------------------------
  subroutine di_algor_vector_product(a,b,c)
    !=========================================================================!
    ! vector c = vector product of vectors a and b                            !
    !                                                                         !
    !                           [cross in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) a (in) : Cartesian vector                                            !
    ! 2) b (in) : Cartesian vector                                            !
    ! 3) c (out) : vector product of a and b                                  !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  none                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    real(kind=di_dp),dimension(1:3),intent(in)  :: a
    real(kind=di_dp),dimension(1:3),intent(in)  :: b
    real(kind=di_dp),dimension(1:3),intent(out) :: c
    !-------------------------------------------------------------------------!

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

    return
  end subroutine di_algor_vector_product

!------------------------------------------------------------------------------

  subroutine di_algor_norm_vector_product(u,v,w)
    !=========================================================================!
    ! compute normalized vector cross product w = u x v                       !
    !                                                                         !
    !                        [normal in DMol]                                 !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) u (in) : Cartesian vector                                            !
    ! 2) v (in) : Cartesian vector                                            !
    ! 3) w (out) : normalized vector cross product                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  none                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none

    real(kind=di_dp),dimension(1:3),intent(in)  :: u
    real(kind=di_dp),dimension(1:3),intent(in)  :: v
    real(kind=di_dp),dimension(1:3),intent(out) :: w
    !-------------------------------------------------------------------------!
    real(kind=di_dp), external :: dnrm2
    real(kind=di_dp),parameter :: tol = 1.0e-08_di_dp
    real(kind=di_dp) :: x
    !-------------------------------------------------------------------------!

    call di_algor_vector_product(u,v,w)

    x = dnrm2(3,w,1)

    if (x>tol) then
       x = one/x
       call dscal(3,x,w(1),1)
    else
       w = zero
    endif

    return

  end subroutine di_algor_norm_vector_product

!------------------------------------------------------------------------------

  subroutine di_algor_norm_vector_difference(u,r,c1,c2)
    !=========================================================================!
    ! compute the difference u between vectors c1 and c2, u = c1-c2           !
    !     (all vectors of dimension 3) and return the norm of u in r          !
    ! The difference vector is normalized on exit                             !
    !                                                                         !
    !                          [vektor in DMol]                               !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) u (out) : unit vector along c1-c2                                    !
    ! 2) r (out) : norm of the c1-c2 difference vector                        !
    ! 3) c1 (in) : Cartesian vector                                           !
    ! 4) c2 (in) : Cartesian vector                                           !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  none                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none

    real (kind=di_dp),dimension(1:3),intent(out)  :: u
    real (kind=di_dp),intent(out)                 :: r       
    real (kind=di_dp),dimension(1:3),intent(in)   :: c1
    real (kind=di_dp),dimension(1:3),intent(in)   :: c2
    !-------------------------------------------------------------------------!
    
    real(kind=di_dp), external :: dnrm2
    real(kind=di_dp),parameter :: tol = 1.0e-08_di_dp
    real(kind=di_dp) :: x
    !-------------------------------------------------------------------------!

    u=c1-c2
    r=dnrm2(3,u,1)

    if (r>tol) then
       x = one/r
       call dscal(3,x,u(1),1)
    else
       r = zero
       u = zero
    endif

    return
  end subroutine di_algor_norm_vector_difference
!------------------------------------------------------------------------------

  subroutine di_algor_distance(v,w,dist)
    !=========================================================================!
    ! computes the difference between vectors V and W, returns its length     !
    ! in dist                                                                 !
    !                                                                         !
    !                              [dista in DMol]                            !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) v (in) : Cartesian vector                                            !
    ! 2) w (in) : Cartesian vector                                            !
    ! 3) dist (out) : Norm of V-W vector                                      !
    !-------------------------------------------------------------------------!
    ! Parent module variables used: none                                      !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  none                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none

    real (kind=di_dp),dimension(1:3),intent(in) :: v
    real (kind=di_dp),dimension(1:3),intent(in) :: w
    real (kind=di_dp),intent(out)               :: dist
    !-------------------------------------------------------------------------!

    real (kind=di_dp),parameter :: eps_dist = 1.0e-14_di_dp
    real (kind=di_dp) :: sum
    !-------------------------------------------------------------------------!

    sum = (v(1)-w(1))*(v(1)-w(1)) + (v(2)-w(2))*(v(2)-w(2)) + (v(3)-w(3))*(v(3)-w(3))

    dist = zero
    if (sum>eps_dist)  dist = sqrt(sum)

    return
  end subroutine di_algor_distance

!------------------------------------------------------------------------------
  subroutine di_algor_unfold_atom_coords(coord,itra,itrb,itrc,tt)
    !=========================================================================!
    ! vector tt = coord + cell *itrabc                                        !
    ! where itrabc describes the cell index                                   !
    !                                                                         !
    !                             [cordn in DMol]                             !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    ! 1) coord (in): Cartesian vector in the zeroth cell                      !
    ! 2) itra (in) : Cell translation along a axis                            !
    ! 3) itrb (in) : Cell translation along b axis                            !
    ! 4) itrc (in) : Cell translation along c axis                            !
    ! 5) tt (out)  : Cartesian vector in the (itra,itrb,itrc) cell            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) elat                                                                 !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:  none                                           !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:  none                                             !
    !-------------------------------------------------------------------------!
    ! Written by Victor Milman, v0.1, 06/12/2002                              !
    !=========================================================================!
    implicit none
    real (kind=di_dp),dimension(1:3),intent(in)   :: coord
    integer,intent(in)                            :: itra,itrb,itrc
    real (kind=di_dp),dimension(1:3),intent(out)  :: tt
    !-------------------------------------------------------------------------!
    real (kind=di_dp)  :: tra,trb,trc
    !-------------------------------------------------------------------------!

    tra = real(itra,kind=di_dp)
    trb = real(itrb,kind=di_dp)
    trc = real(itrc,kind=di_dp)

    tt(1) = coord(1) + tra*elat(1,1)+trb*elat(1,2)+trc*elat(1,3)
    tt(2) = coord(2) + tra*elat(2,1)+trb*elat(2,2)+trc*elat(2,3)
    tt(3) = coord(3) + tra*elat(3,1)+trb*elat(3,2)+trc*elat(3,3)
    return
  end subroutine di_algor_unfold_atom_coords


  subroutine di_algor_fake_gradients
    !=========================================================================!
    ! Initialize Cartesian gradients with 0.1 at the beginning of the run     !
    !                                                                         !
    !                    [rdgrad0 in DMol]                                    !
    !-------------------------------------------------------------------------!
    ! Arguments:  None                                                        !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) number_atoms                                                         !
    ! 2) gradient_cart                                                        !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables: none                                            !
    !-------------------------------------------------------------------------!
    ! Necessary conditions: none                                              !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    implicit none

    !-------------------------------------------------------------------------!
    integer :: i
    !-------------------------------------------------------------------------!

    ! dummy gradient of 0.1 on all but one atom
    do i=1,number_atoms-1
       gradient_cart(1,i)=0.1_di_dp
       gradient_cart(2,i)=0.1_di_dp
       gradient_cart(3,i)=0.1_di_dp
    enddo

    gradient_cart(1,number_atoms)=-(number_atoms-1)*0.1_di_dp
    gradient_cart(2,number_atoms)=-(number_atoms-1)*0.1_di_dp
    gradient_cart(3,number_atoms)=-(number_atoms-1)*0.1_di_dp

    return

  end subroutine di_algor_fake_gradients

!------------------------------------------------------------------------------
  subroutine di_algor_symmetrize_hessian(hold,thrsh,hnew)
    implicit none
    !=========================================================================!
    ! This applies the symmetry operations in locally held 'trans' array      !
    ! to hessian matrix. Additionally zeroes all elements below thrsh         !
    !                                                                         !
    !                         [symhes in DMol]                                !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !                                                                         !
    ! 1) hold (inout): Entry: unsymmetrised hessian, exit: symmetrised hessian!
    !                  (not modified if hnew argument is present)             !
    ! 2) thrsh (in) : threshold below which elements of H will be set to zero !
    ! 3) hnew (out,optional) : symmetrised hessian                            !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    ! 1) number_atoms                                                         !
    ! 2) ntrans                                                               !
    ! 3) trans                                                                !
    ! 4) neqatm                                                               !
    !-------------------------------------------------------------------------!
    ! Modules used: none                                                      !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    ! 1) R and U matrices (keep symmetry matrices and R*H product             !
    ! 2) hcopy matrix: stores the result that is later copied to hold or hnew !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    ! 1) size(hold,1) == size(hold,2) == 3*number_atoms                       !
    ! 2) size(hnew,1) == size(hnew,2) == 3*number_atoms                       !
    !-------------------------------------------------------------------------!
    ! Victor Milman, V0.1, 4th March 2003                                     !
    !=========================================================================!
    
    real(kind=di_dp),dimension(:,:),intent(inout)                 :: hold
    real(kind=di_dp),intent(in)                                   :: thrsh
    real(kind=di_dp),dimension(:,:),intent(out),optional          :: hnew
    !-------------------------------------------------------------------------!
    real(kind=di_dp),dimension(:,:),allocatable :: r
    real(kind=di_dp),dimension(:,:),allocatable :: u
    real(kind=di_dp),dimension(:,:),allocatable :: hcopy
    integer :: iop,iatm,jatm,ii,jj,i,j,k,l,kk,ll,ierr,nat3
    real(kind=di_dp) :: val
    !-------------------------------------------------------------------------!
    nat3 = 3*number_atoms

#ifdef debug
    if (size(hold,1)/=nat3 .or. size(hold,2)/=nat3) &
      &  call deloc_utils_io_abort('Error in di_algor_symmetrize_hessian: dimension of HOLD is wrong')
    if (present(hnew)) then
       if (size(hnew,1)/=nat3 .or. size(hnew,2)/=nat3) &
          &call deloc_utils_io_abort('Error in di_algor_symmetrize_hessian: dimension of HNEW is wrong')
    endif
#endif

    ! allocate memory
    allocate(hcopy(1:nat3,1:nat3),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating hcopy in di_algor_symmetrize_hessian')
    allocate(r(1:nat3,1:nat3),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating r in di_algor_symmetrize_hessian')
    allocate(u(1:nat3,1:nat3),stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in allocating u in di_algor_symmetrize_hessian')

    hcopy = hold
    if (ntrans==1) go to 75

    do iop=2,ntrans

       !  Construct the R matrix for each symmetry operation

       r = zero

       do iatm=1,number_atoms
          ii = (iatm-1)*3
          jatm = neqatm(iatm,iop)
          jj = (jatm-1)*3

          do k=1,3
             kk = ii+k
             do l=1,3
                ll = jj+l
                r(kk,ll) = trans(l,k,iop)
             enddo
          enddo
       enddo

       !  Form R * HESS

       u = zero

       do j=1,nat3
          do k=1,nat3
             val = hold(k,j)
             do i=1,nat3
                u(i,j) = u(i,j) + r(i,k)*val
             enddo
          enddo
       enddo

       !  Form R * HESS * R(t)

       do j=1,nat3
          do k=1,nat3
             val = r(j,k)
             do i=1,nat3
                hcopy(i,j) = hcopy(i,j) + u(i,k)*val
             enddo
          enddo
       enddo

    enddo ! iop

    !  scale the final matrix

    val = one/float(ntrans)
    call dscal(nat3*nat3,val,hcopy(1,1),1)

 75 continue

    !  zero out elements below thrsh

    do i=1,nat3
       do j=1,i
          if (abs(hcopy(i,j))<thrsh) hcopy(i,j) = zero
          hcopy(j,i) = hcopy(i,j)
       enddo
    enddo

    if (present(hnew)) then
       hnew = hcopy
    else
       hold = hcopy
    endif

    ! Free memory
    deallocate(r,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating r in di_algor_symmetrize_hessian')
    deallocate(u,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating u in di_algor_symmetrize_hessian')
    deallocate(hcopy,stat=ierr)
    if (ierr/=0) call deloc_utils_io_abort('Error in deallocating hcopy in di_algor_symmetrize_hessian')

    return
  end subroutine di_algor_symmetrize_hessian
!------------------------------------------------------------------------------

end module geometry_deloc_algor
