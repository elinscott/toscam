! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                          NGWFs  module                         !
!                                                                !
! This module initialises and manipulates NGWFs.                 !
!----------------------------------------------------------------!
! Originally written by Chris-Kriton Skylaris in May 2000.       !
! Subsequent modifications and parallelisation by                !
! Chris-Kriton Skylaris.                                         !
! Automatic NGWF initialisation added on 31/03/2006 by           !
! Chris-Kriton Skylaris.                                         !
! Reorganisation and cleanup by Nicholas Hine in September 2010  !
!================================================================!


module ngwfs

  use constants, only: DP
  use ngwf_data, only: GTO_SET

  implicit none
  private

  ! cks: public functions and subroutines
  public :: ngwfs_initialise
  public :: ngwfs_merge_sets

  ! ndmh: Define type for NGWF radial functions
  type RADIAL_NGWF_TYPE

     ! ndmh: number of shells
     integer :: nshells

     ! ndmh: total number of functions
     integer :: nfuncs

     ! ndmh: number of points of radial grid
     integer :: npts

     ! ndmh: angular momentum of each shell
     integer, pointer, dimension(:) :: angmom

     ! ndmh: cutoff radius for each shell
     real(kind =DP), pointer, dimension(:) :: rc

     ! ndmh: radial grid positions
     real(kind =DP), pointer, dimension(:) :: rad

     ! ndmh: radial parts of NGWF on real space radial grid
     real(kind =DP), pointer, dimension(:,:) :: func_real

  end type RADIAL_NGWF_TYPE


  ! cks: GTO_SET for each supported element of the periodic table
  type(GTO_SET), allocatable, dimension(:) :: gbasis


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_initialise(ngwfs_on_grid,ngwf_basis,elements,read_cond)

    !=======================================================================!
    ! This subroutine initialises the NGWFs of my_node_id ppd-storage form. !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! ngwfs_on_grid (output)  : NGWF data in ppds on psinc grid             !
    ! ngwf_basis (input)      : function basis describing NGWFs             !
    ! elements (input)        : Array of ELEMENT type with data about ions  !
    ! read_cond (input)       : Flag to denote Conduction NGWFs are needed  !
    !-----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000 and revised on 30/5/2001     !
    ! Modified by Chris-Kriton Skylaris on 21/8/2003 so that it works with  !
    ! the parallel version (ONETEP).                                        !
    ! Modified by Chris-Kriton Skylaris on 2/1/2004 so that it can          !
    ! initialise NGWFs from restart info from file.                         !
    ! Modified by Chris-Kriton Skylaris on 15/03/2005 so that it reads      !
    ! fireball sets only from the root proc.                                !
    ! Modified by Alvaro Ruiz Serrano on March 2009 so that it can          !
    ! initialise the NGWFs from its spherical waves expansion saved in a    !
    ! restart file                                                          !
    !-----------------------------------------------------------------------!
    ! Reorganised and mostly re-written by Nicholas Hine to use the new     !
    ! RADIAL_NGWF_TYPE in September 2010.                                   !
    ! Also now warns about unfilled angular momentum shells, September 2010.!
    ! Reorganised again in October 2011 to allow for NGWF radii to be       !
    ! overridden by choices in the pseudoatomic solver, by Nicholas Hine    !
    !=======================================================================!

    use comms, only: comms_bcast, pub_my_node_id, pub_on_root, pub_root_node_id
    use constants, only: DP, stdout, NORMAL
    use function_basis, only: FUNC_BASIS
    use ion, only: element
    use parallel_strategy, only: pub_distr_atom, pub_orig_atom
    use restart, only: restart_ngwfs_tightbox_input, &
         restart_sph_waves_input, restart_ngwfs_tightbox_compose, &
         retrieve_tightbox_ngwfs
    use rundat, only: read_tightbox_ngwfs, pub_read_sw_ngwfs, &
         cond_read_tightbox_ngwfs, pub_output_detail, pub_devel_code
!CW
    use rundat, only: pub_dmft_plot_real_space,pub_dmft_optics
!END CW
    use simulation_cell, only : pub_cell, pub_fftbox
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: ngwf_basis
    real(kind=DP), intent(out) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    ! lr408: Flag to read from conduction NGWF tightbox file if needed
    logical, optional, intent(in) :: read_cond

    ! Local Variables
    character(80) :: ngwf_set
    character(10) :: fmt, tmp
    integer :: iat, orig_iat
    integer :: ingwf, loc_ingwf
    integer :: shell,em      ! shell counter for atomic set
    integer :: isp
    integer :: ierr
    integer :: ingwf_on_atom, ingwf_in_species
    logical :: loc_cond
    logical :: initialise_new
    logical :: shell_filled
    logical :: store_density
    real(kind=DP) :: radius
    type(RADIAL_NGWF_TYPE), allocatable :: radial_ngwfs(:)
!CW
    logical :: check,check_nabla
    INQUIRE(file='store_ham1'  ,EXIST=check)
    INQUIRE(file='store_nabla1',EXIST=check_nabla)
    if(check .and. .not. pub_dmft_plot_real_space .and. .not. (pub_dmft_optics.and..not.check_nabla) ) return
!END CW

    ! ndmh: check for optional arguments
    if (present(read_cond)) then
       loc_cond = read_cond
    else
       loc_cond = .false.
    end if

    ! cks: allocate and initialise the gbasis array, in case it is needed
    if (pub_on_root) then
       call ngwfs_sto3g_create
    endif

    ! ndmh: Create storage for the radial functions for each l
    allocate(radial_ngwfs(pub_cell%num_species),stat=ierr)
    call utils_alloc_check('ngwfs_initialise','radial_ngwfs',ierr)

    ! ndmh: Evaluate or load the radial functions
    do isp=1,pub_cell%num_species

       ! ndmh: Find an example of this species
       do orig_iat=1,pub_cell%nat
          if (elements(orig_iat)%species_number==isp) exit
       end do
       iat = pub_distr_atom(orig_iat)
       radial_ngwfs(isp)%nfuncs = ngwf_basis%num_on_atom(iat)

       ! ndmh: get NGWF file name / solver string for this species
       if (ngwf_basis%name=='ngwfs') then
          ngwf_set = elements(orig_iat)%ngwf_set
       else if (ngwf_basis%name=='ngwfs_cond') then
          ngwf_set = elements(orig_iat)%cond_ngwf_set
       else if (ngwf_basis%name=='ngwfs_aux') then
          ngwf_set = elements(orig_iat)%aux_ngwf_set
       end if

       ! ndmh: get the default sphere radius for this species
       radius = elements(orig_iat)%radius
       if (ngwf_basis%name=='ngwfs_cond') then
          radius = elements(orig_iat)%radius_cond
       end if

       ! ndmh: allocate storage for NGWF radial functions
       call ngwfs_radial_create(radial_ngwfs(isp))

       ! ndmh: now choose how to fill the radial_ngwfs arrays for this species
       if (pub_on_root) then

          if (index(ngwf_set,'AUTO')>0) then

             ! ndmh: STO-3G orbital
             call ngwfs_generate_sto3g(radial_ngwfs(isp),elements(orig_iat))

          else if (index(ngwf_set,'SOLVE')>0) then

             ! ndmh: pseudoatomic solver orbital
             store_density = .false.
             if ((ngwf_basis%name=='ngwfs')) store_density=.true.
             call ngwfs_solve_atom(radial_ngwfs(isp),elements(orig_iat), &
                  ngwf_set,radius,ngwf_basis%num_on_atom(iat),store_density)

          else

             ! ndmh: any other string is assumed to be the name of a fireball
             call ngwfs_read_fireball(radial_ngwfs(isp),elements(orig_iat))

          end if

       end if

       ! ndmh: share the radial function data to all nodes
       call comms_bcast(pub_root_node_id,radial_ngwfs(isp)%nshells)
       call comms_bcast(pub_root_node_id,radial_ngwfs(isp)%npts)
       call comms_bcast(pub_root_node_id,radial_ngwfs(isp)%angmom)
       call comms_bcast(pub_root_node_id,radial_ngwfs(isp)%rc)
       call comms_bcast(pub_root_node_id,radial_ngwfs(isp)%rad)
       call comms_bcast(pub_root_node_id,radial_ngwfs(isp)%func_real)

       ! ndmh: test that this species has no shells of NGWFs where only
       ! ndmh: part of the shell is included in the NGWF set
       ingwf_in_species = 0
       shell_filled = .false.
       do shell=1,radial_ngwfs(isp)%nshells
          shell_filled = .false.
          do em=-radial_ngwfs(isp)%angmom(shell),radial_ngwfs(isp)%angmom(shell)
             ingwf_in_species = ingwf_in_species + 1
             if (em==radial_ngwfs(isp)%angmom(shell)) shell_filled = .true.
             if (ingwf_in_species.eq.ngwf_basis%num_on_atom(iat)) exit
          enddo
          if (ingwf_in_species.eq.ngwf_basis%num_on_atom(iat)) exit
       enddo

       ! ndmh: for final NGWF on atom, check the whole current shell is included
       ! ndmh: in the set of NGWFs. If not, print warning.
       if (.not.shell_filled.and.pub_on_root) then
          write(stdout,*)
          write(stdout,'(a,i2,3a)') 'WARNING in ngwfs_initialise: Setting &
               &radial NGWF functions for species ',isp,' (', &
               trim(elements(orig_iat)%species_id),'):'
          write(stdout,'(a,2(i2,a))') 'Last NGWF on atom has m = ',em, &
               ' and does not complete set of m-values for l =', &
               radial_ngwfs(isp)%angmom(shell),'.'
          write(tmp,'(i10)') radial_ngwfs(isp)%nshells
          write(fmt,'(a,a,a)') '(a,',trim(adjustl(tmp)),'i5)'
          write(stdout,'(a)') 'This breaks symmetry and may cause unstable &
               &NGWF optimisation and poor '
          write(stdout,'(a,2(i3,a))') 'convergence. Suggest increasing number &
               &of NGWFs per atom from ',ngwf_basis%num_on_atom(iat),' to ', &
               ngwf_basis%num_on_atom(iat)-em+radial_ngwfs(isp)%angmom(shell), &
               '.'
          write(stdout,fmt) 'NB: All NGWF radial function angular momenta: ', &
               radial_ngwfs(isp)%angmom(1:radial_ngwfs(isp)%nshells)
          write(stdout,*)
       end if

    end do

    ! ndmh: Check for sphere radius overrides
    do loc_ingwf=1,ngwf_basis%num_on_node(pub_my_node_id)

       ! ndmh: find out global NGWF index, atom index, atom index in input-file
       ! ndmh: order, species number and NGWF index on this atom
       ingwf = loc_ingwf + ngwf_basis%first_on_node(pub_my_node_id) - 1
       iat = ngwf_basis%atom_of_func(ingwf)
       orig_iat = pub_orig_atom(iat)
       isp = elements(orig_iat)%species_number
       ingwf_on_atom = ingwf - ngwf_basis%first_on_atom(iat) + 1

       ! cks: determine the shell number and the z-angular momentum component
       ! cks: of the current function
       ingwf_in_species = 0
       do shell=1,radial_ngwfs(isp)%nshells
          do em=-radial_ngwfs(isp)%angmom(shell),radial_ngwfs(isp)%angmom(shell)
             ingwf_in_species = ingwf_in_species + 1
             if (ingwf_in_species.eq.ingwf_on_atom) exit
          enddo
          if (ingwf_in_species.eq.ingwf_on_atom) exit
       enddo

       ! ndmh: override sphere radius if function is smaller than sphere radius
       if (ngwf_basis%spheres(loc_ingwf)%radius > &
            radial_ngwfs(isp)%rc(shell)) then
          ngwf_basis%spheres(loc_ingwf)%radius = radial_ngwfs(isp)%rc(shell)
       end if

    end do 

    ! ndmh: load NGWFs in from file if required
    initialise_new = .false.
    if (loc_cond .and. cond_read_tightbox_ngwfs) then
       ! cks: initialise NGWFs by reading values stored on disk in:
       ! cks: universal tightbox representation
       call restart_ngwfs_tightbox_input(ngwfs_on_grid, &    ! output
            ngwf_basis, elements, 'tightbox_'//ngwf_basis%name)
    else if (.not. loc_cond .and. read_tightbox_ngwfs &
       .and. .not. retrieve_tightbox_ngwfs) then
       ! cks: initialise NGWFs by reading values stored on disk in
       ! cks: universal tightbox representation
       call restart_ngwfs_tightbox_input(ngwfs_on_grid, &    ! output
            ngwf_basis, elements, 'tightbox_'//ngwf_basis%name)
    else if (.not. loc_cond .and. retrieve_tightbox_ngwfs) then
       ! smmd: initialise NGWFs by mixing NGWFs at previous steps
       call restart_ngwfs_tightbox_compose(ngwfs_on_grid, &    ! output
            ngwf_basis, elements)
    else if (pub_read_sw_ngwfs) then
       ! ars: initialise NGWFs by reading values stored on disk in
       ! ars: spherical waves representation
       call restart_sph_waves_input(ngwfs_on_grid, &     ! output
            ngwf_basis, 'sw_'//ngwf_basis%name)          ! input
    else
       initialise_new = .true.
    end if
    ! ndmh: override: never load in auxiliary NGWF basis
    if (ngwf_basis%name=='ngwfs_aux') then
       initialise_new = .true.
    end if

    ! ndmh: NGWFs were not loaded in, so create new ones
    if (initialise_new) then

       if (pub_on_root.and.(pub_output_detail>=NORMAL)) then
          if (ngwf_basis%name=='ngwfs') then
             write(stdout,'(a)',advance='no') 'NGWF initialisation ...'
          else if (ngwf_basis%name=='ngwfs_cond') then
             write(stdout,'(a)',advance='no') 'Conduction NGWF &
                  &initialisation ...'
          else if (ngwf_basis%name=='ngwfs_aux') then
             write(stdout,'(a)',advance='no') 'Auxiliary NGWF &
                  &initialisation ...'
          end if
       end if

       ! ndmh: Expand radial NGWFs to psinc grid
       if (index(pub_devel_code,'INIT_NGWFS_RECIP')>0) then
          call ngwfs_initialise_from_recip(ngwfs_on_grid,ngwf_basis, &
               radial_ngwfs,elements)
       else
          call ngwfs_initialise_from_radial(ngwfs_on_grid,ngwf_basis, &
               radial_ngwfs,elements)
       end if

       if (pub_on_root.and.(pub_output_detail>=NORMAL)) &
            write(stdout,'(a/)') '... done'

    end if

    ! ndmh: Deallocate the arrays inside the radial functions type
    do isp=pub_cell%num_species,1,-1
       call ngwfs_radial_destroy(radial_ngwfs(isp))
    end do

    ! ndmh: Deallocate storage for the radial functions
    deallocate(radial_ngwfs,stat=ierr)
    call utils_dealloc_check('ngwfs_initialise','radial_ngwfs',ierr)

    ! cks: deallocate the gbasis array
    if (pub_on_root) then
       call ngwfs_sto3g_destroy
    endif

  end subroutine ngwfs_initialise


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_initialise_from_radial(ngwfs_on_grid,ngwf_basis, &
      radial_ngwfs,elements)

    !=======================================================================!
    ! This subroutine initialises the NGWFs of my_node_id ppd-storage form, !
    ! according to the radial_ngwfs defined in the parent routine.          !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! ngwfs_on_grid (output)  : NGWF data in ppds on psinc grid             !
    ! ngwf_basis (input)      : function basis describing NGWFs             !
    ! elements (input)        : Array of ELEMENT type with data about ions  !
    ! radial_ngwfs (input)    : Array of RADIAL_NGWF type with radial part  !
    !                           f(r) for each species                       !
    !-----------------------------------------------------------------------!
    ! Spun off from parent routine by Nicholas Hine, October 2011.          !
    !=======================================================================!

    use basis, only: basis_clean_function, basis_ppd_location
    use comms, only: pub_my_node_id
    use constants, only: DP, stdout, NORMAL
    use geometry, only: POINT, operator(*), operator(+), local_displacement
    use function_basis, only: FUNC_BASIS
    use ion, only: element
    use parallel_strategy, only: pub_orig_atom
    use simulation_cell, only : pub_cell, pub_fftbox
    use spherical_wave, only: sw_init, sw_exit

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: ngwf_basis
    real(kind=DP), intent(out) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    type(RADIAL_NGWF_TYPE), intent(in) :: radial_ngwfs(pub_cell%num_species)

    ! Local Variables
    integer :: point_counter
    integer :: ppd,ppd_loc,loc_1,loc_2,loc_3
    integer :: a1_neighbour,a2_neighbour,a3_neighbour
    integer :: a1_upper,a2_upper,a3_upper
    integer :: a1_lower,a2_lower,a3_lower
    integer :: a1_cell,a2_cell,a3_cell
    integer :: ingwf, loc_ingwf
    integer :: ingwf_on_atom, ingwf_in_species
    integer :: iat, orig_iat
    integer :: shell,em      ! shell counter for atomic set
    integer :: ppd_count     ! ndmh: sphere ppd list counter
    integer :: isp
    integer :: lmax
    logical :: init_sw
    real(kind=DP) :: local_1, local_2, local_3
    type(POINT) :: current_point,current_displacement,periodic_centre

    ! ndmh: to check for ppds included multiple times, when cell=fftbox
    if (pub_fftbox%coin1) then
       a1_lower = -1; a1_upper = 1
    else
       a1_lower =  0; a1_upper = 0
    end if
    if (pub_fftbox%coin2) then
       a2_lower = -1; a2_upper = 1
    else
       a2_lower =  0; a2_upper = 0
    end if
    if (pub_fftbox%coin3) then
       a3_lower = -1; a3_upper = 1
    else
       a3_lower =  0; a3_upper = 0
    end if

    ! ndmh: for high-l
    init_sw = .false.
    lmax = 0
    do isp=1,pub_cell%num_species
       lmax = max(lmax,maxval(radial_ngwfs(isp)%angmom(:)))
    end do
    if (lmax>=5) init_sw = .true.
    if (init_sw) call sw_init(lmax,lmax*2+1)

    ! ndmh: Loop over NGWFs on this node
    ngwfs_on_grid(:) = 0.0_DP
    do loc_ingwf=1,ngwf_basis%num_on_node(pub_my_node_id)

       ! ndmh: find out global NGWF index, atom index, atom index in input-file
       ! ndmh: order, species number and NGWF index on this atom
       ingwf = loc_ingwf + ngwf_basis%first_on_node(pub_my_node_id) - 1
       iat = ngwf_basis%atom_of_func(ingwf)
       orig_iat = pub_orig_atom(iat)
       isp = elements(orig_iat)%species_number
       ingwf_on_atom = ingwf - ngwf_basis%first_on_atom(iat) + 1

       ! cks: determine the shell number and the z-angular momentum component
       ! cks: of the current function
       ingwf_in_species = 0
       do shell=1,radial_ngwfs(isp)%nshells
          do em=-radial_ngwfs(isp)%angmom(shell),radial_ngwfs(isp)%angmom(shell)
             ingwf_in_species = ingwf_in_species + 1
             if (ingwf_in_species.eq.ingwf_on_atom) exit
          enddo
          if (ingwf_in_species.eq.ingwf_on_atom) exit
       enddo

       ! ndmh: loop over all ppds in the NGWF sphere
       do ppd_count=1,ngwf_basis%spheres(loc_ingwf)%n_ppds_sphere

          ppd = ngwf_basis%spheres(loc_ingwf)%ppd_list(1,ppd_count)
          ppd_loc = ngwf_basis%spheres(loc_ingwf)%ppd_list(2,ppd_count)

          do a1_cell=a1_lower,a1_upper
             do a2_cell=a2_lower,a2_upper
                do a3_cell=a3_lower,a3_upper

                   ! ndmh: find which copy of the supercell this NGWF comes from
                   a1_neighbour=nint(real(ppd_loc,kind=DP)/9.0_DP)
                   a2_neighbour=nint(real(ppd_loc-9*a1_neighbour,DP)/3.0_DP)
                   a3_neighbour=ppd_loc-9*a1_neighbour-3*a2_neighbour
                   if (pub_fftbox%coin1) a1_neighbour = a1_cell
                   if (pub_fftbox%coin2) a2_neighbour = a2_cell
                   if (pub_fftbox%coin3) a3_neighbour = a3_cell

                   periodic_centre=ngwf_basis%spheres(loc_ingwf)%centre &
                        + real(a1_neighbour,DP)*pub_cell%a1 &
                        + real(a2_neighbour,DP)*pub_cell%a2 &
                        + real(a3_neighbour,DP)*pub_cell%a3

                   point_counter = ngwf_basis%spheres(loc_ingwf)%offset + &
                        pub_cell%n_pts*(ppd_count-1)

                   ! cks: loop over actual grid points in the current ppd
                   do loc_3=0,pub_cell%n_pt3-1
                      local_3=real(loc_3,DP)*pub_cell%d3

                      do loc_2=0,pub_cell%n_pt2-1
                         local_2=real(loc_2,DP)*pub_cell%d2

                         do loc_1= 0, pub_cell%n_pt1-1
                            local_1=real(loc_1,DP)*pub_cell%d1

                            current_displacement = local_displacement( &
                                 pub_cell%a1_unit,pub_cell%a2_unit, &
                                 pub_cell%a3_unit,local_1,local_2,local_3)

                            current_point = basis_ppd_location(ppd) + &
                                 current_displacement

                            ngwfs_on_grid(point_counter) = &
                                 ngwfs_on_grid(point_counter) + &
                                 ngwfs_at_point(shell,em,radial_ngwfs(isp), &
                                 periodic_centre,current_point, &
                                 ngwf_basis%spheres(loc_ingwf)%radius)

                            point_counter = point_counter + 1

                         enddo  ! loc_1
                      enddo  ! loc_2
                   enddo  ! loc_3

                enddo  ! a3_cell
             enddo  ! a2_cell
          enddo  ! a1_cell

       enddo  ! ppd_count

    end do  ! loc_ingwf

    if (init_sw) call sw_exit

  end subroutine ngwfs_initialise_from_radial


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_initialise_from_recip(ngwfs_on_grid,ngwf_basis, &
      radial_ngwfs,elements)

    !=======================================================================!
    ! This subroutine initialises the NGWFs of my_node_id ppd-storage form, !
    ! according to the radial_ngwfs defined in the parent routine, by first !
    ! transforming them to reciprocal space as projectors, building them in !
    ! the reciprocal space FFTbox, and then extracting them to the ppds.    !
    !-----------------------------------------------------------------------!
    ! Arguments:                                                            !
    ! ngwfs_on_grid (output)  : NGWF data in ppds on psinc grid             !
    ! ngwf_basis (input)      : function basis describing NGWFs             !
    ! elements (input)        : Array of ELEMENT type with data about ions  !
    ! radial_ngwfs (input)    : Array of RADIAL_NGWF type with radial part  !
    !                           f(r) for each species                       !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine, November 2011.                              !
    !=======================================================================!

    use basis, only: basis_clean_function
    use comms, only: pub_my_node_id
    use constants, only: DP, stdout, NORMAL
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT
    use ion, only: element
    use parallel_strategy, only: pub_orig_atom
    use projectors, only: PROJECTOR_SET, projectors_allocate_set, &
         projectors_deallocate_set, projectors_create_real, &
         projectors_init_fftbox_recip, projectors_exit_fftbox_recip
    use services, only: services_regular_transform
    use simulation_cell, only : pub_cell, pub_fftbox
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(inout) :: ngwf_basis
    real(kind=DP), intent(out) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    type(RADIAL_NGWF_TYPE), intent(in) :: radial_ngwfs(pub_cell%num_species)

    ! Local Variables
    type(PROJECTOR_SET) :: ngwf_projectors
    integer :: ingwf
    integer :: iat, orig_iat
    integer :: shell      ! shell counter for atomic set
    integer :: isp
    integer :: proj_count
    integer, parameter :: n_recip_pts = 2001
    real(kind=DP), parameter :: g_max = 50.0_DP

    ! ndmh: allocate storage for ngwf projectors
    ngwf_projectors%n_proj_species = pub_cell%num_species
    call projectors_allocate_set(ngwf_projectors, &
         maxval(radial_ngwfs(:)%nshells),n_recip_pts)

    ! ndmh: set species_num_proj and species_first_proj values
    proj_count = 1
    do isp=1,ngwf_projectors%n_proj_species
       ! Get an example atom of this species
       do iat=1,pub_cell%nat
          orig_iat = pub_orig_atom(iat)
          if (elements(orig_iat)%species_number==isp) exit
       end do
       ! Set number of projectors of each species
       ngwf_projectors%species_num_proj(isp) = ngwf_basis%num_on_atom(iat)
       ngwf_projectors%species_first_proj(isp) = proj_count
       ngwf_projectors%gmax(isp) = g_max
       ngwf_projectors%n_rad_pts(isp) = n_recip_pts
       ngwf_projectors%num_shells(isp) = radial_ngwfs(isp)%nshells
       ngwf_projectors%ang_mom(:,isp) = 0
       ngwf_projectors%rad_proj_recip(:,:,isp) = 0.0_DP
       do shell=1,radial_ngwfs(isp)%nshells
          ngwf_projectors%ang_mom(shell,isp) = radial_ngwfs(isp)%angmom(shell)
          call services_regular_transform(radial_ngwfs(isp)%angmom(shell),2, &
               radial_ngwfs(isp)%npts,maxval(radial_ngwfs(isp)%rad), &
               n_recip_pts,g_max,radial_ngwfs(isp)%func_real(:,shell), &
               ngwf_projectors%rad_proj_recip(1:n_recip_pts,shell,isp))
       end do
       proj_count = proj_count + ngwf_projectors%species_num_proj(isp)

    end do

    ! ndmh: copy in projector centres and radii
    do iat=1,pub_cell%nat
       orig_iat = pub_orig_atom(iat)
       isp = elements(orig_iat)%species_number
       ngwf_projectors%proj_centre(iat) = elements(orig_iat)%centre
       ngwf_projectors%proj_max_radius(iat) = &
            maxval(ngwf_basis%spheres(:)%radius)
       ngwf_projectors%proj_species(iat) = isp
    end do

    ! ndmh: initialise projectors in reciprocal-space fftbox representation
    call projectors_init_fftbox_recip(ngwf_projectors)

    ! ndmh: create projectors in recip space & extract from FFTbox to ppds
    call projectors_create_real(ngwf_basis,ngwf_projectors,ngwfs_on_grid)

    ! ndmh: zero regions outside sphere boundary but inside nonzero ppds
    do ingwf=1,ngwf_basis%node_num
       call basis_clean_function(ngwfs_on_grid,ngwf_basis%spheres(ingwf), &
            ngwf_basis%n_ppds)
    end do

    call projectors_exit_fftbox_recip(ngwf_projectors)
    call projectors_deallocate_set(ngwf_projectors)

  end subroutine ngwfs_initialise_from_recip


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_merge_sets(ngwfs_on_grid,ngwf_basis,ngwfs_on_grid_src1, &
       ngwf_basis_src1,ngwfs_on_grid_src2,ngwf_basis_src2)

    !====================================================================!
    ! This subroutine merges two sets of NGWFs together to form a joint  !
    ! set.                                                               !
    !====================================================================!
    ! Arguments:                                                         !
    !   ngwfs_on_grid (output) : Merged NGWFs in ppd storage             !
    !   ngwf_basis (input)     : Spheres etc describing merged set       !
    !   ngwfs_on_grid_src1 (in): NGWF set 1 to merge (in ppds)           !
    !   ngwf_basis_src1 (in)   : Spheres etc describing NGWF set 1       !
    !   ngwfs_on_grid_src2 (in): NGWF set 2 to merge (in ppds)           !
    !   ngwf_basis_src2 (in)   : Spheres etc describing NGWF set 2       !
    !====================================================================!
    ! Written by Nicholas Hine in April 2011.                            !
    !====================================================================!

    use comms, only: pub_my_node_id
    use function_basis, only: FUNC_BASIS
    use parallel_strategy, only: pub_num_atoms_on_node, pub_first_atom_on_node
    use simulation_cell, only: pub_cell, pub_fftbox
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(out) :: ngwfs_on_grid(ngwf_basis%n_ppds* &
         pub_cell%n_pts)
    type(FUNC_BASIS), intent(in) :: ngwf_basis_src1
    real(kind=DP), intent(in) :: ngwfs_on_grid_src1(ngwf_basis_src1%n_ppds* &
         pub_cell%n_pts)
    type(FUNC_BASIS), intent(in) :: ngwf_basis_src2
    real(kind=DP), intent(in) :: ngwfs_on_grid_src2(ngwf_basis_src2%n_ppds* &
         pub_cell%n_pts)

    ! Local Variables
    integer :: ifunc_src1, ifunc_src2
    integer :: iat, loc_iat
    integer :: count, count_src1, count_src2
    integer :: start, start_src, finish, finish_src
    character(len=10) :: iat_str

    ! ndmh: loop over all atoms on node and copy the spheres of each one
    count = 0
    count_src1 = 0
    count_src2 = 0
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       iat = pub_first_atom_on_node(pub_my_node_id) + loc_iat - 1

       if ((ngwf_basis_src1%num_on_atom(iat) + &
            ngwf_basis_src2%num_on_atom(iat)) &
            /=ngwf_basis%num_on_atom(iat)) then
          write(iat_str,'(i10)') iat
          iat_str = adjustl(trim(iat_str))
          call utils_abort('Error in ngwfs_merge_sets: &
               &Mismatching NGWF counts on atom '//iat_str)
       end if

       do ifunc_src1=1,ngwf_basis_src1%num_on_atom(iat)
          count_src1 = count_src1 + 1
          count = count + 1

          ! Find start and finish in dest array
          start = ngwf_basis%spheres(count)%offset
          finish = start - 1 + &
               ngwf_basis%spheres(count)%n_ppds_sphere*pub_cell%n_pts

          ! Find start and finish in src array
          start_src = ngwf_basis_src1%spheres(count_src1)%offset
          finish_src = start_src - 1 + &
               ngwf_basis_src1%spheres(count_src1)%n_ppds_sphere*pub_cell%n_pts

          ! Copy NGWF from src1 to dest array
          ngwfs_on_grid(start:finish) = ngwfs_on_grid_src1(start_src:finish_src)
       end do

       do ifunc_src2=1,ngwf_basis_src2%num_on_atom(iat)
          count_src2 = count_src2 + 1
          count = count + 1

          ! Find start and finish in dest array
          start = ngwf_basis%spheres(count)%offset
          finish = start - 1 + &
               ngwf_basis%spheres(count)%n_ppds_sphere*pub_cell%n_pts

          ! Find start and finish in src array
          start_src = ngwf_basis_src2%spheres(count_src2)%offset
          finish_src = start_src - 1 + &
               ngwf_basis_src2%spheres(count_src2)%n_ppds_sphere*pub_cell%n_pts

          ! Copy NGWF from src2 to dest array
          ngwfs_on_grid(start:finish) = ngwfs_on_grid_src2(start_src:finish_src)

       end do

    end do

    if (count/=ngwf_basis%node_num) then
       call utils_abort('Error in ngwfs_merge_sets: Wrong total number of &
            &NGWFs copied into merged set')
    end if

  end subroutine ngwfs_merge_sets


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function ngwfs_at_point(shell,em,radial_ngwf,periodic_centre, &
       current_point,ngwf_radius)

    !=================================================================!
    ! This function returns the value of an ngwf at a given point in  !
    ! real space given a radial representation of the ngwf.           !
    !-----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.                       !
    ! f-functions added by Chris-Kriton Skylaris on 9/5/2007.         !
    ! Modified by Nicholas Hine on 03/09/2010 to use RADIAL_NGWF_TYPE.!
    !=================================================================!

    use comms, only: comms_abort, pub_my_node_id
    use constants, only: stdout
    use geometry,  only: point, operator(-), geometry_distance
    use services, only: services_linear_interpolation, services_locate_interp
    use spherical_wave, only: sw_real_sph_harm
    use utils, only: utils_isnan

    implicit none

    ! Arguments
    real(kind=DP) ::  ngwfs_at_point
    integer, intent(in) :: shell
    integer, intent(in) :: em
    type(RADIAL_NGWF_TYPE),intent(in) :: radial_ngwf
    real(kind=DP), intent(in) :: ngwf_radius
    type(POINT), intent(in) :: current_point
    type(POINT), intent(in) :: periodic_centre

    ! Local variables
    real(kind=DP) :: local_distance, radial_value, angular_value
    integer :: interp_index
    type(point) :: r_vec

    ngwfs_at_point = 0.0_DP ! qoh: initialise to prevent compiler warning

    ! cks: here centre can be the centre of a ngwf function inside the
    ! cks: simulation cell but it can also be the centre of a periodic
    ! cks: image of a ngwf function which happens to have values in a
    ! cks: ppd of the simulation cell.

    local_distance = geometry_distance(periodic_centre,current_point)
    if (local_distance <= ngwf_radius) then

       ! cks: find the array index of the lowest point for the intepolation,
       ! cks: provided the point does not fall out of the region of validity
       ! cks: of the fireball ngwf function.in the opposite case set
       ! cks: radial_value=0.0_DP
       if (local_distance<radial_ngwf%rad(radial_ngwf%npts)) then

          interp_index = services_locate_interp(local_distance, &
               radial_ngwf%rad,radial_ngwf%npts)

          radial_value = services_linear_interpolation(local_distance, &
               radial_ngwf%func_real(interp_index,shell), &
               radial_ngwf%func_real(interp_index+1,shell), &
               radial_ngwf%rad(interp_index), &
               radial_ngwf%rad(interp_index+1))
       else
          radial_value = 0.0_DP
       endif

       r_vec = current_point - periodic_centre
       angular_value = sw_real_sph_harm(r_vec%x,r_vec%y,r_vec%z, &
            local_distance,radial_ngwf%angmom(shell),em)

       ngwfs_at_point = radial_value*angular_value

    endif

  end function ngwfs_at_point


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_read_fireball(radial_ngwf,current_element)

    !==================================================================!
    ! This function reads from a file a set of NGWFs for a single atom.!
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris in 2000.                        !
    ! Improved by Peter D. Haynes in 2004.                             !
    ! Modified by Nicholas Hine on 03/09/2010 to use RADIAL_NGWF_TYPE. !
    !==================================================================!

    use comms, only : comms_abort
    use constants, only: DP, stdout
    use ion, only: ELEMENT
    use utils, only: utils_unit, utils_open_unit_check, utils_close_unit_check

    implicit none

    ! Arguments
    type(RADIAL_NGWF_TYPE), intent(inout) :: radial_ngwf
    type(ELEMENT), intent(in) :: current_element

    ! CKS: INTERNAL VARIABLES
    character(len=80) :: shell_name,line
    character(len=3) :: shell_string, mom_string
    character(len=64) :: ngwf_set
    character(len=80) :: current_job
    integer :: row, shell, num_shells
    integer :: ierr
    integer :: iunit

    iunit = utils_unit()
    ngwf_set = current_element%ngwf_set

    open(iunit,file=trim(ngwf_set),status='old',position='rewind', &
         iostat=ierr)
    call utils_open_unit_check('ngwfs_read_fireball',ngwf_set,ierr)

    ! cks : read angular momenta
    current_job = 'seeking marker "ANGULAR_MOMENTA"'
    do
       read(iunit,'(a80)',err=100,end=200) line
       if (index(line,'ANGULAR_MOMENTA') > 0) exit
    end do

    current_job = 'reading block "ANGULAR_MOMENTA"'
    num_shells = 0
    do
       read(iunit,'(a80)',err=100,end=200) line
       if (index(line,'END_ANGULAR_MOMENTA') > 0) exit

       num_shells = num_shells + 1
       if (num_shells <= radial_ngwf%nshells) then
          read(line,*) radial_ngwf%angmom(num_shells)
       end if
    end do

    radial_ngwf%nshells = num_shells

    ! cks : read radial positions
    rewind(iunit,iostat=ierr,err=400)
400 if (ierr /= 0) then
       write(stdout,'(3a,i6)') &
            'Error in ngwfs_read_set: rewinding file "', &
            trim(ngwf_set),'" failed with code ',ierr
       call comms_abort
    end if

    current_job = 'seeking marker "RADIAL_POSITIONS"'
    do
       read(iunit,'(a80)',err=100,end=200) line
       if (index(line,'RADIAL_POSITIONS') > 0) exit
    end do

    current_job = 'reading block "RADIAL_POSITIONS"'
    do row=1,radial_ngwf%npts
       read(iunit,'(a80)',err=100,end=200) line
       if (index(line,'END RADIAL_POSITIONS') > 0) exit
       if (index(line,'END_RADIAL_POSITIONS') > 0) exit
       read(line,*,err=100,end=200) radial_ngwf%rad(row)
    end do
    radial_ngwf%npts = row - 1
    radial_ngwf%rc(:) = max(radial_ngwf%rad(radial_ngwf%npts), &
         current_element%radius)

    ! cks : read radial functions
    do shell=1,num_shells

       write(shell_string,'(i3)') shell
       write(mom_string,'(i3)') radial_ngwf%angmom(shell)
       write(shell_name,'(a80)') 'SHELL_'//trim(adjustl(shell_string))// &
            '_ANGMOM_'//adjustl(mom_string)
       shell_name = adjustl(shell_name)

       rewind(iunit,iostat=ierr,err=500)
500    if (ierr /= 0) then
          write(stdout,'(3a,i6)') &
               'Error in ngwfs_read_set: rewinding file "', &
               trim(ngwf_set),'" failed with code ',ierr
          call comms_abort
       end if

       write(current_job,'(a80)') 'seeking marker "'//trim(shell_name)//'"'
       current_job = adjustl(current_job)
       do
          read(iunit,'(a80)',err=100,end=200) line
          if (trim(adjustl(line)) == shell_name) exit
       end do

       write(current_job,'(a80)') 'reading block "'//trim(shell_name)//'"'
       current_job = adjustl(current_job)
       do row=1,radial_ngwf%npts
          read(iunit,*,err=100,end=200) radial_ngwf%func_real(row,shell)
       end do

    end do

    close(iunit,iostat=ierr)
    call utils_close_unit_check('ngwfs_read_fireball',ngwf_set,ierr)

    return

100 write(stdout,'(4a)') 'Error in ngwfs_read_fireball: reading file "', &
         trim(ngwf_set),'" failed while ',trim(current_job)
    call comms_abort

200 write(stdout,'(4a)') 'Error in ngwfs_read_fireball: file "', &
         trim(ngwf_set),'" ended unexpectedly while ',trim(current_job)
    call comms_abort

  end subroutine ngwfs_read_fireball


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_generate_sto3g(radial_ngwf,current_element)

    !=========================================================!
    ! This subroutine initialises a radial fireball set for a !
    ! single atom to a set of STO-3G contracted Gaussian (CG) !
    ! functions taking into account also that first CG should !
    ! correspond to the first valence atomic function of the  !
    ! atom.                                                   !
    !---------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 31/03/2006.         !
    ! Modified by Nicholas Hine on 03/09/2010 to use          !
    ! the new RADIAL_NGWF_TYPE.                               !
    !=========================================================!

    use comms, only : comms_abort
    use constants, only: DP, stdout
    use ion, only: element
    use ngwf_data, only: pub_num_gtatoms
    use services, only: services_equally_spaced_numbers

    implicit none

    ! Arguments
    type(RADIAL_NGWF_TYPE), intent(inout) :: radial_ngwf
    type(ELEMENT), intent(in) :: current_element

    ! Local Variables
    real(kind=DP), parameter :: max_rad = 10.0_DP ! max extent of Gaussians
    real(kind=DP)  :: rval   ! distance in radial grid
    real(kind=DP)  :: rvsq   ! squared distance in radial grid
    integer :: shell          ! shell counter
    integer :: start_shell    ! starting STO-3G valence shell to use
    real(kind=DP) :: shell_charge   ! electron charge of core shells
    real(kind=DP) :: ion_charge     ! ionic charge of currrent atom
    integer :: fbl_shells     ! number of fireball shells constructed from STO-3G
    integer :: atnum          ! atomic number of current atom
    integer :: row            ! radial grid point counter
    integer :: num_shells

    atnum = current_element%atomic_number
    ion_charge = current_element%ion_charge

    ! cks: stop if set not available
    if (atnum > pub_num_gtatoms) then
       write(stdout,*)'Error: Automatic NGWF initialisation for element', &
            'with atomic number:',atnum,' not supported!'
       call comms_abort
    endif

    ! cks: find starting shell of STO-3G basis from which to start fireball set
    ! cks: (depends on ion charge of pseudopotential)
    start_shell = -1
    shell_charge = 0
    find_start_shell_loop: do shell=1,gbasis(atnum)%nshells
       shell_charge = shell_charge + 2*(2*gbasis(atnum)%angmom(shell) + 1)
       if (shell_charge >= ( atnum-ion_charge ) ) then
          if (shell_charge == ( atnum-ion_charge ) ) start_shell = shell + 1
          if (shell_charge > ( atnum-ion_charge ) ) start_shell = shell
          exit find_start_shell_loop
       endif
    enddo find_start_shell_loop

    ! ndmh: total number of shells in STO3G set
    fbl_shells = gbasis(atnum)%nshells - start_shell + 1
    if (fbl_shells > radial_ngwf%nshells) then
       write(stdout,'(a,i3,a,i3,a)') 'Error in internal_generate_sto3g_set: &
            &number of shells in sto3g set for atomic number ',atnum, &
            ' exceeds maximum (',radial_ngwf%nshells,')'
       call comms_abort
    end if

    ! cks: now initialise the radial grid point positions
    call services_equally_spaced_numbers(radial_ngwf%rad(:), &
         0.0_DP,max_rad,radial_ngwf%npts)

    num_shells = 0
    do shell=start_shell,gbasis(atnum)%nshells

       num_shells = num_shells + 1

       ! cks: initialise angular momentum for current shell
       radial_ngwf%angmom(num_shells) = gbasis(atnum)%angmom(shell)
       radial_ngwf%rc(num_shells) = radial_ngwf%rad(radial_ngwf%npts)

       ! cks: generate current radial contracted Gaussian
       do row=1,radial_ngwf%npts

          rval = radial_ngwf%rad(row)
          rvsq = rval**2

          ! cks: contract the Gaussian
          radial_ngwf%func_real(row,num_shells) = &
               gbasis(atnum)%coco(shell, 1) * &
               exp(-gbasis(atnum)%expo(shell,1)*rvsq) &
               + gbasis(atnum)%coco(shell, 2) * &
               exp(-gbasis(atnum)%expo(shell,2)*rvsq) &
               + gbasis(atnum)%coco(shell, 3) * &
               exp(-gbasis(atnum)%expo(shell,3)*rvsq)

          ! cks: multiply contraction with r^ang_mom
          if (radial_ngwf%angmom(num_shells) > 0) then
             radial_ngwf%func_real(row,num_shells) = &
                  radial_ngwf%func_real(row,num_shells) * &
                  rval**radial_ngwf%angmom(num_shells)
          endif

       enddo
    enddo

    radial_ngwf%nshells = num_shells

  end subroutine ngwfs_generate_sto3g


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_solve_atom(radial_ngwf,current_element,elem_string, &
       target_radius,target_nfuncs,store_density)

    !==================================================================!
    ! This routine solves the Schrodinger equation for a single ion,   !
    ! and creates the corresponding starting radial functions for the  !
    ! NGWFs.                                                           !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 03/09/2010.                          !
    !==================================================================!

    use atom, only: ATOM_TYPE, BASIS_TYPE, atom_create, atom_destroy, &
         atom_create_basis, atom_destroy_basis, atom_solve, &
         atom_write_orbitals, atom_get_lmax, atom_split_orbital, &
         atom_polarise_orbital
    use comms, only : comms_abort
    use constants, only: DP, stdout
    use density, only: density_radial_store
    use ion, only: ELEMENT
    use rundat, only: cutoff_energy, pub_rootname, &
         pub_initial_dens_realspace, pub_cond_calculate
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(RADIAL_NGWF_TYPE), intent(inout) :: radial_ngwf
    type(ELEMENT), intent(in) :: current_element
    character(80), intent(in) :: elem_string
    real(kind=DP), intent(in) :: target_radius
    integer, intent(in) :: target_nfuncs
    logical, intent(in) :: store_density

    ! Local Variables
    type(ATOM_TYPE) :: tatom
    type(BASIS_TYPE) :: tbasis
    character(len=120) :: report(400)
    character(80) :: config
    character(80) :: filename
    character(20) :: tmp
    real(kind=DP) :: cwidth, cwidth_l(0:5), cscale
    real(kind=DP) :: rmax, rmax_l(0:5)
    integer :: nreport
    integer :: isplitnorm, nsplitnorm
    integer :: iorb, jorb, oorb
    integer :: ishell, ireport
    integer :: ipol
    integer :: pos
    integer :: lmax
    logical :: fix_derivs

    ! Parse the string following "SOLVE"
    pos = index(elem_string,'conf=')
    if (pos>0) then
       config = elem_string(pos+5:)
    else
       config = ''
    end if

    ! Read the confining potential width
    cwidth_l(:) = 3.0_DP
    pos = index(elem_string,'w=')
    if (pos>0) then
       tmp = elem_string(pos+2:)
       read(tmp,*) cwidth
       cwidth_l(:) = cwidth
    end if
    pos = index(elem_string,'ws=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) cwidth_l(0)
    end if
    pos = index(elem_string,'wp=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) cwidth_l(1)
    end if
    pos = index(elem_string,'wd=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) cwidth_l(2)
    end if
    pos = index(elem_string,'wf=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) cwidth_l(3)
    end if
    pos = index(elem_string,'wg=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) cwidth_l(4)
    end if
    pos = index(elem_string,'wh=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) cwidth_l(5)
    end if

    ! Read the cutoff radius/radii
    rmax_l(:) = target_radius
    pos = index(elem_string,'R=')
    if (pos>0) then
       tmp = elem_string(pos+2:)
       read(tmp,*) rmax
       rmax_l(:) = rmax
    end if
    pos = index(elem_string,'Rs=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) rmax_l(0)
    end if
    pos = index(elem_string,'Rp=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) rmax_l(1)
    end if
    pos = index(elem_string,'Rd=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) rmax_l(2)
    end if
    pos = index(elem_string,'Rf=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) rmax_l(3)
    end if
    pos = index(elem_string,'Rg=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) rmax_l(4)
    end if
    pos = index(elem_string,'Rh=')
    if (pos>0) then
       tmp = elem_string(pos+3:)
       read(tmp,*) rmax_l(5)
    end if

    ! Read the confining potential height
    pos = index(elem_string,'S=')
    if (pos>0) then
       tmp = elem_string(pos+2:)
       read(tmp,*) cscale
    else
       cscale = 100.0_DP
    end if

    ! Parse the string following "SOLVE"
    pos = index(elem_string,'FD')
    if (pos>0) then
       fix_derivs = .true.
    else
       fix_derivs = .false.
    end if

    ! Set up the basis and the atom
    call atom_get_lmax(lmax,current_element,target_nfuncs,config)
    call atom_create_basis(tbasis,rmax_l,lmax,6.0_DP*cutoff_energy,fix_derivs)
    call atom_create(tatom,tbasis,current_element,target_nfuncs, &
         config,cwidth_l,cscale)

    ! Solve the atom
    call atom_solve(tatom,tbasis,report,nreport)

    ! Print report
    do ireport=1,nreport
       write(stdout,'(a)') trim(report(ireport))
    end do

    ! Copy the atom solution into the radial NGWF
    if (tbasis%npts > size(radial_ngwf%rad)) then
       call utils_abort('Error in ngwfs_solve_atom: Atom grid is larger than &
            &RADIAL_NGWF grid')
    end if
    radial_ngwf%npts = tbasis%npts
    radial_ngwf%rad(1:radial_ngwf%npts) = tbasis%rad(1:tbasis%npts)
    radial_ngwf%nshells = 0

    ! Find orbitals which need polarising, and generate polarisation orbitals
    jorb = tatom%norbs - sum(tatom%polarise) + 1
    do iorb=1,tatom%norbs - sum(tatom%polarise)
       ! If splitnorm == -1, ignore this orbital
       if (tatom%splitnorm(1,iorb)==-1.0_DP) cycle

       ! Generate the polarisation orbital
       do ipol=1,tatom%polarise(iorb)
          if (tatom%orb_ang_mom(iorb)+ipol/=tatom%orb_ang_mom(jorb)) then
             write(stdout,'(a,i3,a)') &
                 'Error polarising orbital ',iorb,': angular momenta do not &
                 &match'
          end if
          if (ipol==1) then
             oorb = iorb
          else
             oorb = jorb - 1
          end if
          write(stdout,'(a,i3,a,i2,a,i2,a)') &
              'Polarising orbital',oorb,' to generate l=', &
               tatom%orb_ang_mom(oorb) + 1,' function (orbital',jorb,')'
          call atom_polarise_orbital(tatom%psi_r(1:tbasis%npts,jorb), &
               tatom%psi_r(1:tbasis%npts,oorb),tatom%eigs(iorb), &
               tatom%orb_ang_mom(oorb),tatom,tbasis)
          jorb = jorb + 1
       end do
    end do

    ! Find orbitals which need splitting, and transfer them to radial NGWF type
    ishell = 0
    do iorb=1,tatom%norbs
       ! If splitnorm == -1, ignore this orbital
       if (tatom%splitnorm(1,iorb)==-1.0_DP) cycle

       ! Find number of target functions to split out of this orbital
       nsplitnorm = count(tatom%splitnorm(:,iorb)>0.0_DP) + 1
       ! Copy the orbital into the radial NGWF nsplitnorm times over
       do isplitnorm=1,nsplitnorm
          ishell = ishell + 1
          radial_ngwf%nshells = radial_ngwf%nshells + 1
          radial_ngwf%angmom(ishell) = tatom%orb_ang_mom(iorb)
          radial_ngwf%rc(ishell) = tbasis%rmax_l(tatom%orb_ang_mom(iorb))
          radial_ngwf%func_real(1:radial_ngwf%npts,ishell) = &
               tatom%psi_r(1:tbasis%npts,iorb)
          if (isplitnorm>1) write(stdout,'(a,i3,a,f12.9)') &
              'Splitting orbital ',iorb,', splitnorm=', &
               tatom%splitnorm(isplitnorm-1,iorb)
       end do
       ! Now do the splitting
       call atom_split_orbital( &
            radial_ngwf%func_real(:,ishell-nsplitnorm+1:ishell), &
            tatom%psi_r(1:tbasis%npts,iorb),nsplitnorm, &
            tatom%splitnorm(:,iorb),tatom%orb_ang_mom(iorb), &
            tatom%work,tbasis)

    end do

    ! Copy the atom density into the radial density storage
    if (pub_initial_dens_realspace.and.store_density) then
       call density_radial_store(current_element%species_number,tbasis%npts, &
            tbasis%rad(1:tbasis%npts),tatom%den(1:tbasis%npts), &
            tatom%comp_den(1:tbasis%npts))
    end if

    ! Write the orbitals out, if required
    if (.true.) then
       write(tmp,'(i3.3)') current_element%species_number
       if (any(tatom%splitnorm(:,:)>0.0_DP).or.(any(tatom%polarise(:)>0))) then
          filename = trim(pub_rootname)//'_atom_orbs_'//tmp
          call atom_write_orbitals(tatom,tbasis,filename)
       end if
       filename = trim(pub_rootname)//'_initial_rad_ngwf_'//tmp
       call ngwfs_radial_write(radial_ngwf,filename)
    end if

    ! Clean up the atom and basis
    call atom_destroy(tatom)
    call atom_destroy_basis(tbasis)

  end subroutine ngwfs_solve_atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_radial_write(radial_ngwf,filename)

    !=====================================================================!
    ! This subroutine allocates storage for a RADIAL_NGWF_TYPE structure. !
    !---------------------------------------------------------------------!
    ! Written by Nicholas Hine on 03/09/2010.                             !
    !=====================================================================!

    use utils, only: utils_open_unit_check, utils_close_unit_check, utils_unit

    implicit none

    ! Arguments
    type(RADIAL_NGWF_TYPE), intent(inout) :: radial_ngwf
    character(len=80), intent(in) :: filename

    ! Local Variables
    integer :: ierr
    integer :: iunit
    integer :: ipt
    integer :: ishell

    ! ndmh: open the file
    iunit = utils_unit()
    open(iunit,file=trim(filename),iostat=ierr)
    call utils_open_unit_check('ngwfs_radial_write',trim(filename),ierr)

    ! ndmh: write the orbitals
    do ipt=1,radial_ngwf%npts
       write(iunit,'(f22.15)',advance='no') radial_ngwf%rad(ipt)
       do ishell=1,radial_ngwf%nshells
          write(iunit,'(f22.15)',advance='no') radial_ngwf%func_real(ipt,ishell)
       end do
       write(iunit,*)
    end do

    ! ndmh: close the file
    close(iunit,iostat=ierr)
    call utils_close_unit_check('ngwfs_radial_write',trim(filename),ierr)

  end subroutine ngwfs_radial_write


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_radial_create(radial_ngwf)

    !=====================================================================!
    ! This subroutine allocates storage for a RADIAL_NGWF_TYPE structure. !
    !---------------------------------------------------------------------!
    ! Written by Nicholas Hine on 03/09/2010.                             !
    !=====================================================================!

    use utils, only: utils_alloc_check

    implicit none

    ! Arguments
    type(RADIAL_NGWF_TYPE), intent(inout) :: radial_ngwf

    ! Local Variables
    integer :: ierr

    ! ndmh: Default maximum sizes (not necessarily all used)
    radial_ngwf%nshells = 100
    radial_ngwf%npts = 4001

    allocate(radial_ngwf%angmom(radial_ngwf%nshells),stat=ierr)
    call utils_alloc_check('ngwfs_radial_create','radial_ngwf%angmom',ierr)
    allocate(radial_ngwf%rc(radial_ngwf%nshells),stat=ierr)
    call utils_alloc_check('ngwfs_radial_create','radial_ngwf%rc',ierr)
    allocate(radial_ngwf%rad(radial_ngwf%npts),stat=ierr)
    call utils_alloc_check('ngwfs_radial_create','radial_ngwf%rad',ierr)
    allocate(radial_ngwf%func_real(radial_ngwf%npts,radial_ngwf%nshells),stat=ierr)
    call utils_alloc_check('ngwfs_radial_create','radial_ngwf%func_real',ierr)

    ! ndmh: Default maximum size for a fireball
    radial_ngwf%npts = 2001
    radial_ngwf%angmom = 0
    radial_ngwf%rc = 0.0_DP
    radial_ngwf%rad = 0.0_DP
    radial_ngwf%func_real = 0.0_DP

  end subroutine ngwfs_radial_create


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_radial_destroy(radial_ngwf)

    !=======================================================================!
    ! This subroutine deallocates storage for a RADIAL_NGWF_TYPE structure. !
    !-----------------------------------------------------------------------!
    ! Written by Nicholas Hine on 03/09/2010.                               !
    !=======================================================================!

    use utils, only: utils_dealloc_check

    implicit none

    ! Arguments
    type(RADIAL_NGWF_TYPE), intent(inout) :: radial_ngwf

    ! Local Variables
    integer :: ierr

    deallocate(radial_ngwf%func_real,stat=ierr)
    call utils_dealloc_check('ngwfs_radial_destroy','radial_ngwf%func_real',ierr)
    deallocate(radial_ngwf%rad,stat=ierr)
    call utils_dealloc_check('ngwfs_radial_destroy','radial_ngwf%rad',ierr)
    deallocate(radial_ngwf%rc,stat=ierr)
    call utils_dealloc_check('ngwfs_radial_destroy','radial_ngwf%rc',ierr)
    deallocate(radial_ngwf%angmom,stat=ierr)
    call utils_dealloc_check('ngwfs_radial_destroy','radial_ngwf%angmom',ierr)

  end subroutine ngwfs_radial_destroy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_normalise(ngwf_on_grid, &   ! INPUT-OUTPUT
       ngwf_basis)                             ! INPUT

    !=========================================================!
    ! Normalise a set of NGWFs to unity.                      !
    !---------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 23/10/2001.         !
    ! Modified by Quintin Hill to use pub_cell on 15/10/2008. !
    ! Modified to use function basis type by Nicholas Hine on !
    ! 14/12/2009.                                             !
    !=========================================================!

    use comms, only: comms_abort
    use constants, only: DP,stdout
    use function_basis, only: FUNC_BASIS
    use simulation_cell, only: pub_cell

    implicit none

    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(inout) :: ngwf_on_grid( &
         ngwf_basis%n_ppds*pub_cell%n_pts)

    ! cks: internal declarations
    integer :: row, start, finish
    real(kind=DP) :: norm_constant

    do row=1,ngwf_basis%num
       start=ngwf_basis%spheres(row)%offset
       finish=start+(ngwf_basis%spheres(row)%n_ppds_sphere)*pub_cell%n_pts -1

       norm_constant=sum(ngwf_on_grid(start:finish)*ngwf_on_grid(start:finish) )
       norm_constant=sqrt(norm_constant*pub_cell%weight)

       if (norm_constant.gt.0) then
          ngwf_on_grid(start:finish)=ngwf_on_grid(start:finish)/norm_constant
       else
          write (stdout,*) 'norm_constant=',norm_constant,'in ngwfs_normalise'
          write (stdout,*) 'for function ',row,'  ONETEP execution stops'
          call comms_abort
       endif

    enddo

  end subroutine ngwfs_normalise


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_sto3g_create

    !===================================================================!
    ! Allocate appropriate memory and initialise it to hold Gaussian    !
    ! basis set parameters (all STO-3G and only the polarisation        !
    ! functions of 6-31G*) for each element up to atomic number         !
    ! pub_num_gtatoms
    !-------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 31/03/2006.                   !
    !===================================================================!

    use constants, only: DP
    use ngwf_data, only: pub_num_gtatoms, &
       ngwf_data_cocos_and_expos01_20, &
       ngwf_data_cocos_and_expos21_40, &
       ngwf_data_cocos_and_expos41_60, &
       ngwf_data_cocos_and_expos61_80, &
       ngwf_data_cocos_and_expos81_103
    use utils, only: utils_alloc_check
    implicit none

    ! cks: <<local variables>>
    integer :: atom     ! atom counter
    integer :: n_shells ! number of shells
    integer :: n_prim   ! number of primitives in given shell
    integer :: ierr     ! allocation error flag


    ! cks: allocate gbasis array for all supported elements
    allocate(gbasis(pub_num_gtatoms), stat=ierr)
    call utils_alloc_check('ngwfs_sto3g_create','gbasis',ierr)

    ! cks: Number of shells per element, including one polarisation shell
    ! cks: This data was typed by Shyong Chen
    gbasis(1)%nshells=2
    gbasis(2)%nshells=2
    gbasis(3)%nshells=4
    gbasis(4)%nshells=4
    gbasis(5)%nshells=4
    gbasis(6)%nshells=4
    gbasis(7)%nshells=4
    gbasis(8)%nshells=4
    gbasis(9)%nshells=4
    gbasis(10)%nshells=4
    gbasis(11)%nshells=6
    gbasis(12)%nshells=6
    gbasis(13)%nshells=6
    gbasis(14)%nshells=6
    gbasis(15)%nshells=6
    gbasis(16)%nshells=6
    gbasis(17)%nshells=6
    gbasis(18)%nshells=6
    gbasis(19)%nshells=8
    gbasis(20)%nshells=8
    gbasis(21)%nshells=9
    gbasis(22)%nshells=9
    gbasis(23)%nshells=9
    gbasis(24)%nshells=9
    gbasis(25)%nshells=9
    gbasis(26)%nshells=9
    gbasis(27)%nshells=9
    gbasis(28)%nshells=9
    gbasis(29)%nshells=9
    gbasis(30)%nshells=9
    gbasis(31)%nshells=9
    gbasis(32)%nshells=9
    gbasis(33)%nshells=9
    gbasis(34)%nshells=9
    gbasis(35)%nshells=9
    gbasis(36)%nshells=9
    gbasis(37)%nshells=11
    gbasis(38)%nshells=11
    gbasis(39)%nshells=12
    gbasis(40)%nshells=12
    gbasis(41)%nshells=12
    gbasis(42)%nshells=12
    gbasis(43)%nshells=12
    gbasis(44)%nshells=12
    gbasis(45)%nshells=12
    gbasis(46)%nshells=12
    gbasis(47)%nshells=12
    gbasis(48)%nshells=12
    gbasis(49)%nshells=12
    gbasis(50)%nshells=12
    gbasis(51)%nshells=12
    gbasis(52)%nshells=12
    gbasis(53)%nshells=12
    gbasis(54)%nshells=13
    gbasis(55)%nshells=14
    gbasis(56)%nshells=14
    gbasis(57)%nshells=15
    gbasis(58)%nshells=15
    gbasis(59)%nshells=15
    gbasis(60)%nshells=15
    gbasis(61)%nshells=15
    gbasis(62)%nshells=15
    gbasis(63)%nshells=15
    gbasis(64)%nshells=15
    gbasis(65)%nshells=15
    gbasis(66)%nshells=15
    gbasis(67)%nshells=15
    gbasis(68)%nshells=15
    gbasis(69)%nshells=15
    gbasis(70)%nshells=15
    gbasis(71)%nshells=15
    gbasis(72)%nshells=15
    gbasis(73)%nshells=15
    gbasis(74)%nshells=15
    gbasis(75)%nshells=15
    gbasis(76)%nshells=15
    gbasis(77)%nshells=15
    gbasis(78)%nshells=15
    gbasis(79)%nshells=15
    gbasis(80)%nshells=15
    gbasis(81)%nshells=15
    gbasis(82)%nshells=15
    gbasis(83)%nshells=15
    gbasis(84)%nshells=15
    gbasis(85)%nshells=15
    gbasis(86)%nshells=16
    gbasis(87)%nshells=16
    gbasis(88)%nshells=16
    gbasis(89)%nshells=16
    gbasis(90)%nshells=17
    gbasis(91)%nshells=18
    gbasis(92)%nshells=18
    gbasis(93)%nshells=18
    gbasis(94)%nshells=18
    gbasis(95)%nshells=17
    gbasis(96)%nshells=17
    gbasis(97)%nshells=17
    gbasis(98)%nshells=17
    gbasis(99)%nshells=17
    gbasis(100)%nshells=17
    gbasis(101)%nshells=17
    gbasis(102)%nshells=17
    gbasis(103)%nshells=17


    ! cks: allocate memory for the components of each shell
    do atom=1, pub_num_gtatoms

       n_shells = gbasis(atom)%nshells
       n_prim =3

       ! cks: allocate angular momentum memory
       allocate(gbasis(atom)%angmom(n_shells), stat=ierr)
       call utils_alloc_check('ngwfs_sto3g_create','gbasis%angmom',ierr)
       ! cks: initialise
       gbasis(atom)%angmom =-1


       ! cks: allocate exponents memory
       allocate(gbasis(atom)%expo(n_shells, n_prim), stat=ierr)
       call utils_alloc_check('ngwfs_sto3g_create','gbasis%expo',ierr)
       ! cks: initialise
       gbasis(atom)%expo =0.0_DP

       ! cks: allocate contraction coefficients memory
       allocate(gbasis(atom)%coco(n_shells, n_prim), stat=ierr)
       call utils_alloc_check('ngwfs_sto3g_create','gbasis%coco',ierr)
       ! cks: initialise
       gbasis(atom)%coco =0.0_DP

    enddo



    ! cks: load cocos and expos in the allocated memory
    call ngwf_data_cocos_and_expos01_20(gbasis(:))
    call ngwf_data_cocos_and_expos21_40(gbasis(:))
    call ngwf_data_cocos_and_expos41_60(gbasis(:))
    call ngwf_data_cocos_and_expos61_80(gbasis(:))
    call ngwf_data_cocos_and_expos81_103(gbasis(:))


  end subroutine ngwfs_sto3g_create


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_sto3g_destroy

    !=======================================================!
    ! Deallocate all memory of gbasis array.                !
    !---------------------------------------- --------------!
    ! Written by Chris-Kriton Skylaris on 31/03/2006.       !
    !=======================================================!

    use ngwf_data, only: pub_num_gtatoms
    use utils, only: utils_dealloc_check
    implicit none

    ! cks: <<local variables>>
    integer :: atom     ! atom counter
    integer :: ierr     ! allocation error flag

    ! cks: deallocate memory for the components of each shell
    do atom=1, pub_num_gtatoms

       ! cks: deallocate angular momentum memory
       deallocate(gbasis(atom)%angmom, stat=ierr)
       call utils_dealloc_check('ngwfs_sto3g_destroy','gbasis%angmom',ierr)

       ! cks: deallocate exponents memory
       deallocate(gbasis(atom)%expo, stat=ierr)
       call utils_dealloc_check('ngwfs_sto3g_destroy','gbasis%expo',ierr)

       ! cks: deallocate contraction coefficients memory
       deallocate(gbasis(atom)%coco, stat=ierr)
       call utils_dealloc_check('ngwfs_sto3g_destroy','gbasis%coco',ierr)

    enddo

    ! cks: deallocate the gbasis array
    deallocate(gbasis, stat=ierr)
    call utils_dealloc_check('ngwfs_sto3g_destroy','gbasis',ierr)


  end subroutine ngwfs_sto3g_destroy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module ngwfs



