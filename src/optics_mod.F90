! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!=================================================================!
!                                                                 !
!                          Optics Module                          !
!                                                                 !
! This module calculates optical absorption spectra.              !
!-----------------------------------------------------------------!
! Relocated from properties and slightly re-organised by Laura    !
! Ratcliff in November 2011.                                      !
!=================================================================!

module optics

  use constants, only: DP

  implicit none

  private

  public :: optics_calculate_spectra

!CW
  public :: optics_grad_mat_els
!END CW

contains

  subroutine optics_calculate_spectra(cur_spin,eigs_dens,eigen_en,&
       rep,ngwf_basis,proj_basis,num_opt_states)

    !==================================================================!
    ! This subroutine calculates optical absorption spectra using      !
    ! either position or momentum matrix elements (including the       !
    ! commutator between the non-local potential and the position      !
    ! operator).  Both the matrix elements and the imaginary component !
    ! of the dielectric function are written to file.                  !
    !------------------------------------------------------------------!
    ! Moved from properties_spectra, which was written by Laura        !
    ! Ratcliff in December 2010 and updated to include the output of   !
    ! the imaginary component of the dielectric function by Laura      !
    ! Ratcliff in October 2011.                                        ! 
    !==================================================================!

    use comms, only: pub_on_root, comms_barrier
    use constants, only: stdout, HARTREE_IN_EVS, UP, DN, PI
    use dense, only: DEM, dense_get_element, dense_create, dense_destroy
    use function_basis, only: FUNC_BASIS
    use geometry, only: operator(.CROSS.), operator(.DOT.)
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use rundat, only: pub_rootname, pub_calc_mom_mat_els, pub_any_nl_proj, &
         pub_spectra_print_mat_els, pub_opt_smear
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy
    use utils, only:  utils_unit, utils_open_unit_check, utils_close_unit_check,&
         utils_assert, utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(NGWF_REP), intent(in) :: rep
    type(DEM), intent(in) :: eigs_dens
    real(kind=DP), intent(in) :: eigen_en(:) ! hamiltonian eigenvalues
    integer, intent(in) :: cur_spin
    ! lr408: integer defining the total number of states to include for optical spectra
    integer, intent(in) :: num_opt_states

    ! Local variables
    type(DEM) :: opt_mat_elements(3) ! Optical matrix elements of eigenvectors
    type(SPAM3)   :: pos_elements(3)
    real(kind=DP) :: cell_volume
    real(kind=DP) :: opt_mat_el
    real(kind=DP), allocatable, dimension(:) :: trans_energies, trans_weights, &
         trans_weights_av
    complex(kind=DP) :: opt_mat_el_cmplx
    integer :: xyz, homo, i, f, count, num_energies
    integer :: basis_num
    integer :: output_unit, io_status
    integer :: max ! number of states to output for matrix elements
    integer :: num_val, num_cond, ierr
    character(len=256) :: output_file  ! output file name buffer
    character(len=6) :: file_type

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &optics_calculate_spectra'
#endif

    ! Create dense matrices to store optical matrix elements between eigenstates
    do xyz=1,3
       call dense_create(opt_mat_elements(xyz),ngwf_basis%num,ngwf_basis%num, &
            iscmplx=pub_calc_mom_mat_els)
    end do

    ! lr408: At some point, actually make use of this to scale matrix elements
    ! lr408: within ONETEP
    cell_volume = abs((pub_cell%a1 .CROSS. pub_cell%a2) .DOT. pub_cell%a3)
    call utils_assert(cell_volume /= 0.0_DP,'Cell volume is zero')

    ! lr408: Create sparse matrices
    if (pub_calc_mom_mat_els) then
       if (pub_any_nl_proj) then
          do xyz=1,3
             call sparse_create(pos_elements(xyz),rep%nonlocpot(1), &
                  iscmplx=.true.)
          end do
       else ! non-local commutator not needed - but still complex matrix elems
          do xyz=1,3
             call sparse_create(pos_elements(xyz),rep%overlap,iscmplx=.true.)
          end do
       end if
    else ! position matrix elements only (all real)
       do xyz=1,3
          call sparse_create(pos_elements(xyz),rep%overlap)
       end do
    end if

    ! lr408: Calculate position or momentum matrix elements
    if (.not. pub_calc_mom_mat_els) then
       call optics_pos_mat_els(rep%ngwfs_on_grid, ngwf_basis, &
            rep%ngwfs_on_grid, ngwf_basis, proj_basis, rep%overlap, &
            rep%sp_overlap, pos_elements)
    else
       call optics_grad_mat_els(pos_elements,rep,ngwf_basis,proj_basis)
    end if

    basis_num = ngwf_basis%num

    homo = rep%n_occ(cur_spin)
    if (homo==0) homo=1

    do xyz=1,3
       call optics_mat_els(opt_mat_elements(xyz),ngwf_basis%num, &
            eigs_dens,pos_elements(xyz),pub_calc_mom_mat_els)
    end do

    ! write out matrix elements for comparison with CASTEP

    ! Determine prefix from NGWF_REP postfix
    if (rep%postfix=='') then
       file_type='_val'
    else if (rep%postfix=='j') then
       file_type='_joint'
    end if

    ! lr408: changed so matrix elements are printed and spectra are calculated
    ! lr408: for all valence + number of optimised conduction states
    max = num_opt_states 
   
    num_val = homo
    num_cond = max - homo

    ! only print out matrix elements if required
    if (pub_on_root .and. pub_spectra_print_mat_els) then

       write(stdout,'(a)') ''

       write(stdout,'(a)')&
            '================ Writing optical matrix elements &
            &================'
       if (pub_cell%num_spins == 1) then
          write(output_file,*)trim(pub_rootname)//trim(file_type)//'_OPT_MAT_ELS.txt'
       else if (cur_spin == UP) then
          write(output_file,*)trim(pub_rootname)//trim(file_type)//'_OPT_MAT_ELS_UP.txt'
       else if (cur_spin == DN) then
          write(output_file,*)trim(pub_rootname)//trim(file_type)//'_OPT_MAT_ELS_DN.txt'
       end if

       output_file = adjustl(output_file)

       ! cks: print output warning
       write(stdout,'(3a)') 'Writing "', trim(output_file),'" ...  done'

       write(stdout,'(a)')'================================&
            &================================================'
       write(stdout,'(a)') ''

       ! cks: get a unit number that is free
       output_unit = utils_unit()

       open(unit=output_unit, form="formatted" ,file=trim(output_file), &
            action="write",iostat=io_status)
       call utils_open_unit_check('optics_calculate_spectra','output_file', &
            io_status)

       ! cks: write first line
       write(output_unit,'(a)')'# initial | Energy (eV) | final  | Energy (eV) |&
           & Matrix el.  | Trans. E (eV) '
       !write(output_unit,'(a,2x,I2)') '# Spin component',cur_spin

       write(output_unit,'(I5,2x,I5)') max,max,1
       write(output_unit,'(F24.12)') cell_volume

    end if

    ! calculate imaginary part of dielectric function and joint DOS

    if (pub_opt_smear > 0) then

       num_energies = num_val * num_cond
 
       allocate(trans_energies(num_energies),stat=ierr)
       call utils_alloc_check('optics_calculate_spectra','trans_energies',ierr)
 
       allocate(trans_weights(num_energies),stat=ierr)
       call utils_alloc_check('optics_calculate_spectra','trans_weights',ierr)

       allocate(trans_weights_av(num_energies),stat=ierr)
       call utils_alloc_check('optics_calculate_spectra','trans_weights_av',ierr)

       trans_weights_av = 0

       count = 1

       ! only between valence and optimized conduction states
       do i=1,homo
          do f=homo+1,max
             trans_energies(count) = eigen_en(f) - eigen_en(i)
             count = count+1
          end do
       end do

       ! first do JDOS - strictly speaking the scissor operator is meaningless here
       ! but we might want to compare JDOS with weighted JDOS (i.e. imag diel fn)
       call optics_dos_gp(num_energies,trans_energies,file_type,cur_spin)

    end if

    ! eae32: dense_get_element calls across all nodes
    do xyz=1,3
       
       if (pub_on_root .and. pub_spectra_print_mat_els) then
          write(output_unit,'(a,I1)') '# Cartesian component ',xyz
       end if

       count = 1

       do i=1,max ! initial state
          do f=1,max ! final state
             
             if (.not.pub_calc_mom_mat_els) then
                call dense_get_element(opt_mat_el,opt_mat_elements(xyz),i,f)
                if (pub_on_root .and. pub_spectra_print_mat_els) then
                   write(output_unit,'(2x,2(I5,4x,F12.6,2x),2(F24.12,2x))') &
                        i,eigen_en(i)*HARTREE_IN_EVS,&
                        f,eigen_en(f)*HARTREE_IN_EVS,&
                        opt_mat_el,&
                        (eigen_en(f)-eigen_en(i))*HARTREE_IN_EVS
                end if

                ! get matrix elements for this direction and generate imag diel fn
                if (pub_opt_smear > 0) then
                   ! only between valence and (optimised) conduction states
                   if (i <= homo .and. f > homo) then
                      trans_weights(count) = opt_mat_el
                      count = count + 1
                   end if

                end if

             else ! momentum matrix elements

                ! print imag mat els as well
                call dense_get_element(opt_mat_el_cmplx, &
                     opt_mat_elements(xyz),i,f)
                if (pub_on_root .and. pub_spectra_print_mat_els) then
                   if (i /= f) then
                      ! eae32: set matrix elements to zero for like states
                      write(output_unit,'(2x,2(I5,4x,F12.6,2x),4(F24.12,2x))') &
                           i,eigen_en(i)*HARTREE_IN_EVS,&
                           f,eigen_en(f)*HARTREE_IN_EVS,&
                           abs(opt_mat_el_cmplx/(eigen_en(f)-eigen_en(i))),&
                           (eigen_en(f)-eigen_en(i))*HARTREE_IN_EVS,&
                           opt_mat_el_cmplx/(eigen_en(f)-eigen_en(i))
                   else
                      write(output_unit,'(2x,2(I5,4x,F12.6,2x),4(F24.12,2x))') &
                           i,eigen_en(i)*HARTREE_IN_EVS,&
                           f,eigen_en(f)*HARTREE_IN_EVS,&
                           0.0,(eigen_en(f)-eigen_en(i))*HARTREE_IN_EVS,&
                           0.0,0.0
                   end if
                end if

                ! get matrix elements for this direction and generate imag diel fn
                if (pub_opt_smear > 0) then

                   ! only between valence and conduction states - assume for now want up to
                   ! max cond states                   
                   if (i <= homo .and. f > homo) then
                      trans_weights(count) = abs(opt_mat_el_cmplx/(eigen_en(f)-eigen_en(i)))
                      count = count + 1
                   end if

                end if

             end if ! position / momentum matrix elements             
             
          end do
       end do

       trans_weights = trans_weights**2
       trans_weights = trans_weights * (8.0_dp * PI**2 / cell_volume) * HARTREE_IN_EVS

       trans_weights_av = trans_weights_av + trans_weights

       call optics_dos_gp(num_energies,trans_energies,file_type,cur_spin,trans_weights,pos=xyz)

    end do ! loop over xyz
    
    if (pub_opt_smear > 0) then
       ! generate imag diel fn for average of all 3 directions
       trans_weights_av = trans_weights_av / 3.0_dp 

       call optics_dos_gp(num_energies,trans_energies,file_type,cur_spin,&
            trans_weights_av,pos=5)   
    end if

    if (pub_opt_smear > 0) then    

       deallocate(trans_energies,stat=ierr)
       call utils_dealloc_check('optics_calculate_spectra','trans_energies',ierr)
 
       deallocate(trans_weights,stat=ierr)
       call utils_dealloc_check('optics_calculate_spectra','trans_weights',ierr)

       deallocate(trans_weights_av,stat=ierr)
       call utils_dealloc_check('optics_calculate_spectra','trans_weights_av',ierr)

    end if

    if (pub_on_root .and. pub_spectra_print_mat_els) then
       
       close(unit=output_unit,iostat=io_status)
       call utils_close_unit_check('optics_calculate_spectra','output_unit', &
            io_status)

    end if

    call comms_barrier

    do xyz=3,1,-1
       call sparse_destroy(pos_elements(xyz))
    end do


    do xyz=3,1,-1
       call dense_destroy(opt_mat_elements(xyz))
    end do


#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving optics_calculate_spectra'
#endif

  end subroutine optics_calculate_spectra


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine optics_mat_els(opt_mat_elements,basis_num, &
       eigs_dens,pos_elements,calc_grad)

    !==================================================================!
    ! This subroutine takes position or momentum matrix elements       !
    ! between NGWFs and multiplies them by eigenvector coefficients    !
    ! to give matrix elements between states for use in generating     !
    ! optical absorption spectra.                                      !
    !------------------------------------------------------------------!
    ! Arguments                                                        !
    !   opt_mat_elements (output) :                                    !
    !   basis_num (input)         :                                    !
    !   eigs_dens (input)         :                                    !
    !   pos_elements (input)      :                                    !
    !   calc_grad (input)         :                                    !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff in December 2010.                      !
    ! Cleaned up, parallelised and made to use DEM and dense_product   !
    ! by Nicholas Hine, April 2011.                                    !
    ! Moved from properties_opt_mat_els to optics_mat_els Nov 2011.    !
    !==================================================================!

    use comms, only: pub_on_root, comms_abort, comms_barrier
    use constants, only: stdout, cmplx_i
    use dense, only: DEM, dense_axpy, dense_scale, &
         dense_create, dense_destroy, dense_convert, dense_product
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_axpy, &
         sparse_scale

    implicit none

    ! Arguments
    type(DEM), intent(inout) :: opt_mat_elements
    type(DEM), intent(in) :: eigs_dens
    type(SPAM3), intent(in) :: pos_elements
    integer, intent(in) :: basis_num
    logical, intent(in) :: calc_grad

    ! Local variables
    type(SPAM3) :: pos_elements_real
    type(SPAM3) :: pos_elements_imag
    type(SPAM3) :: pos_elements_swap

    type(DEM) :: pos_dens_real
    type(DEM) :: pos_dens_imag
    type(DEM) :: eigs_pos_prod
    type(DEM) :: result_rc

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering optics_mat_els'
#endif

    ! ndmh: Set opt_mat_elements to zero
    call dense_scale(opt_mat_elements,0.0_DP)

    ! ndmh: Create real/imag parts of pos_dens
    if (calc_grad) then
       call dense_create(pos_dens_real,basis_num,basis_num)
       call dense_create(pos_dens_imag,basis_num,basis_num)

       ! eae32: create and initialise to sparse matrix to contain the 
       !        real and imaginary part of the complex sparse matrix pos_elements
       call sparse_create(pos_elements_real,pos_elements,iscmplx=.false.)
       call sparse_create(pos_elements_imag,pos_elements,iscmplx=.false.)
       call sparse_create(pos_elements_swap,pos_elements,iscmplx=.true.)
       call sparse_axpy(pos_elements_swap,pos_elements,-cmplx_i)

       call sparse_scale(pos_elements_real,0.0_DP)
       call sparse_scale(pos_elements_imag,0.0_DP)
       call dense_scale(pos_dens_real,0.0_DP)
       call dense_scale(pos_dens_imag,0.0_DP)

       ! eae32: no complex version of spam3toblacs for SCALAPACK
       !        --> instead copy contents of complex sparse pos_elements
       !            to real sparse matrices pos_elements_real and 
       !            pos_elements_imag
       !        --> then convert separately to dense pos_dens_real
       !            and pos_dens_imag
       call sparse_axpy(pos_elements_real,pos_elements,1.0_DP)
       ! lr408: no sparse_axpy option for ymat(real), xmat(cmplx), alpha(cmplx)
       ! lr408: so added pos_elements_swap matrix
       call sparse_axpy(pos_elements_imag,pos_elements_swap,1.0_DP)
       call sparse_destroy(pos_elements_swap)

       call dense_convert(pos_dens_real,pos_elements_real)
       call dense_convert(pos_dens_imag,pos_elements_imag)
       call sparse_destroy(pos_elements_real)
       call sparse_destroy(pos_elements_imag)

    else ! not calc_grad so everything real
       call dense_create(pos_dens_real,basis_num,basis_num)
       call dense_convert(pos_dens_real,pos_elements)
    end if

    ! ndmh: Create temporary matrices
    call dense_create(eigs_pos_prod,basis_num,basis_num)
    call dense_create(result_rc,basis_num,basis_num)

    ! ndmh: Find real part of R_ij = M_i^a Re(<phi_a|r|phi_b>) (M*)^b_j
    call dense_product(eigs_pos_prod,pos_dens_real,eigs_dens, &
         transpose_amat=.true.)
    call dense_product(result_rc,eigs_pos_prod,eigs_dens, &
         transpose_amat=.true.)
    ! ndmh: Add to R_ij matrix
    call dense_axpy(opt_mat_elements,result_rc,1.0_DP)

    if (calc_grad) then
       ! ndmh: Find imag part of R_ij = M_i^a Im(<phi_a|r|phi_b>) (M*)^b_j
       call dense_product(eigs_pos_prod,pos_dens_imag,eigs_dens, &
            transpose_amat=.true.)
       call dense_product(result_rc,eigs_dens,eigs_pos_prod, &
            transpose_amat=.true.)
       ! ndmh: Add to R_ij matrix
       call dense_axpy(opt_mat_elements,result_rc,cmplx_i)
    end if

    ! eae32: Synchronise nodes before writing optical matrix elements to file
    call comms_barrier

    ! ndmh: Clean up temporary matrices
    call dense_destroy(result_rc)
    call dense_destroy(eigs_pos_prod)
    if (calc_grad) then
       call dense_destroy(pos_dens_imag)
    end if
    call dense_destroy(pos_dens_real)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving optics_mat_els'
#endif

  end subroutine optics_mat_els


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine optics_pos_mat_els(bras_on_grid, bra_basis, &
            kets_on_grid, ket_basis, proj_basis, overlap, sp_overlap, &
            r_elements)

    !==================================================================!
    ! This subroutine calculates position matrix elements for use in   !
    ! optical absorption spectra calculations.  This is only expected  !
    ! to be useful for molecules where no NGWFs overlap with copies    !
    ! from neigbouring cells.                                          !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff in December 2010, based on             !
    ! properties_polarisation.                                         !
    ! Modified for PAW by Nicholas Hine in December 2011.              !
    !==================================================================!

    use augmentation, only: augmentation_pos
    use basis, only: basis_location_func_wrt_cell, basis_ket_start_wrt_fftbox
    use comms, only: pub_on_root, pub_my_node_id
    use constants, only: stdout
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(*), operator(+)
    use integrals, only: integrals_pos
    use ion, only: element
    use rundat, only: pub_aug
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: sparse_create, sparse_get_element, SPAM3, &
         sparse_destroy, sparse_put_element, sparse_element_exists, &
         sparse_expand, PATTERN_ALTERNATE
    use timer, only: timer_clock

    implicit none

    ! Arguments
    type(SPAM3), intent(in) :: overlap
    type(SPAM3), intent(inout):: r_elements(3)
    type(FUNC_BASIS), intent(in) :: bra_basis
    type(FUNC_BASIS), intent(in) :: ket_basis
    real(kind=DP), intent(in) :: bras_on_grid(bra_basis%n_ppds*pub_cell%n_pts)
    real(kind=DP), intent(in) :: kets_on_grid(ket_basis%n_ppds*pub_cell%n_pts)
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(SPAM3), intent(in) :: sp_overlap

    ! Local Variables
    type(POINT)   :: a1,a2,a3
    type(POINT)   :: r
    real(kind=DP) :: R_fft(3),r_el,o_el
    integer       :: cs1,cs2,cs3
    integer       :: bs1,bs2,bs3
    integer       :: xyz,jngwf,loc_ingwf,ingwf

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering optics_pos_mat_els'
#endif

    ! Start timer
    call timer_clock('optics_pos_mat_els',1)

    ! calculate matrix elements < phi_a | r_op - R_fft | phi_b >
    call integrals_pos(r_elements,&
         bras_on_grid,bra_basis,kets_on_grid,ket_basis,1)

    ! Add R_fft * S_ab to matrix elements
    ! NOT IMPLEMENTED YET: accounting for atoms moving across cell boundary
    ! ndmh: NB this involves an O(N^2) check - could use the
    ! ndmh: index of S to calculate which elements need shifting
    ! ndmh: but this is a one-time calculation and the prefactor is tiny

    ! Calculate vectors between grid points
    a1 = (1.0_DP / pub_cell%total_pt1) * pub_cell%a1
    a2 = (1.0_DP / pub_cell%total_pt2) * pub_cell%a2
    a3 = (1.0_DP / pub_cell%total_pt3) * pub_cell%a3

    ! Determine where tightbox begins wrt fftbox
    call basis_ket_start_wrt_fftbox(bs1,bs2,bs3,pub_fftbox%total_pt1,&
         pub_fftbox%total_pt2,pub_fftbox%total_pt3)

    ! Loop over all ngwfs on this node
    do loc_ingwf=1,ket_basis%num_on_node(pub_my_node_id)
       ingwf = loc_ingwf + ket_basis%first_on_node(pub_my_node_id) - 1

       ! Find start of tightbox of ingwf
       call basis_location_func_wrt_cell(cs1,cs2,cs3, &
            ket_basis%tight_boxes(loc_ingwf))

       ! Find vector to origin of FFTbox
       r = real(cs1-bs1,kind=DP) * a1 &
            + real(cs2-bs2,kind=DP) * a2 &
            + real(cs3-bs3,kind=DP) * a3

       R_fft(1) = r%X
       R_fft(2) = r%Y
       R_fft(3) = r%Z

       ! Loop over all row ngwfs
       do jngwf=1,bra_basis%num

          ! Test if these NGWFs overlap
          if (.not.sparse_element_exists(overlap,jngwf,ingwf)) cycle

          ! Extract overlap element
          call sparse_get_element(o_el,overlap,jngwf,ingwf)

          ! Extract element from r_elements and shift by R_fft*o_el
          do xyz=1,3

             call sparse_get_element(r_el,r_elements(xyz),jngwf,ingwf)
             r_el = R_fft(xyz) * o_el + r_el
             call sparse_put_element(r_el,r_elements(xyz),jngwf,ingwf)
          enddo

       enddo
    enddo

    ! ndmh: in PAW/USP formalism, add contribution from sphere part
    if (pub_aug) call augmentation_pos(r_elements,proj_basis,sp_overlap)

    ! Stop timer
    call timer_clock('optics_pos_mat_els',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving optics_pos_mat_els'
#endif


  end subroutine optics_pos_mat_els


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine optics_grad_mat_els(pos_elements,rep,ngwf_basis,proj_basis)

    !==================================================================!
    ! This subroutine calculates momentum matrix elements (including   !
    ! the commutator between the non-local potential and the position  !
    ! operator) between NGWFs for optical absorption spectra           !
    ! calculations.                                                    !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff in November 2011, based on             !
    ! properties_spectra, which was written by Laura Ratcliff in       ! 
    ! December 2010.                                                   !
    ! Modified for PAW by Nicholas Hine in December 2011.              !
    !==================================================================!

    use augmentation, only: augmentation_grad,aug_nonlocal_commutator_mat
    use comms, only: pub_on_root, comms_barrier
    use constants, only: stdout
    use dense, only: DEM, dense_get_element
    use function_basis, only: FUNC_BASIS
    use geometry, only: operator(.CROSS.), operator(.DOT.)
    use integrals, only: integrals_grad
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use pseudopotentials, only: pseudo_nonlocal_commutator_mat, &
         pseudo_nonlocal_com_mat_fd
    use rundat, only: pub_calc_mom_mat_els, pub_any_nl_proj, &
         pub_spec_nonloc_fin_diff, pub_calc_nonloc_comm, &
         pub_spec_cont_deriv, pub_spectra_print_mat_els, pub_opt_smear, &
         pub_paw, pub_aug
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_scale, &
         sparse_axpy, sparse_copy
    use timer, only: timer_clock
    use utils, only: utils_abort

    implicit none

    ! Arguments
    type(SPAM3), intent(inout)   :: pos_elements(3)
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(NGWF_REP), intent(in) :: rep
    ! Local variables
    type(SPAM3)   :: pos_elements_rc(3)
    type(SPAM3)   :: nonloc_com_elements(3)
    real(kind=DP) :: delta
    integer :: xyz

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering optics_grad_mat_els'
#endif

    ! Start timer
    call timer_clock('optics_grad_mat_els',1)

    delta = pub_spec_nonloc_fin_diff

    ! lr408: Use relationship q.<f|r|i> = q.{<f|p|i>/(im*w_fi) + <f|[V_NL,r]|i>/(hbar*w_fi)}
    ! lr408: to calculate momentum matrix elements rather than position

    ! lr408: Calculate grad mat els then commutator
    ! lr408: (NB both need dividing by w_fi and we are in atomic units, so can
    ! lr408: ignore constants)
    do xyz=1,3
       call sparse_create(pos_elements_rc(xyz),pos_elements(xyz),iscmplx=.false.)
    end do

    call integrals_grad(pos_elements_rc, rep%ngwfs_on_grid, ngwf_basis, &
         rep%ngwfs_on_grid, ngwf_basis)

    ! ndmh: for PAW/USP, augment the grad operator matrix
    if (pub_aug) call augmentation_grad(pos_elements_rc,rep%sp_overlap)

    if (pub_any_nl_proj.or.pub_aug) then

       if (pub_calc_nonloc_comm) then

          do xyz=1,3
             call sparse_create(nonloc_com_elements(xyz),rep%nonlocpot(1),iscmplx=.true.)
          end do

          if (.not. pub_spec_cont_deriv) then ! i.e. finite diff

             if (pub_any_nl_proj) then
                call pseudo_nonlocal_com_mat_fd(nonloc_com_elements, &
                     proj_basis, ngwf_basis, rep%ngwfs_on_grid, &
                     rep%sp_overlap, rep%nonlocpot(1), delta)
             else if (pub_aug) then
                call utils_abort('Error in optics_grad_mat_els: &
                     &aug_nonlocal_com_mat_fd does not exist')
             end if

             do xyz=1,3
                call sparse_scale(pos_elements_rc(xyz),-1.0_dp)
                call sparse_copy(pos_elements(xyz),pos_elements_rc(xyz))
                call sparse_axpy(pos_elements(xyz),nonloc_com_elements(xyz),(1.0_DP,0.0_DP))
             end do

          else ! continuous derivative

             if (pub_any_nl_proj) then
                call pseudo_nonlocal_commutator_mat(nonloc_com_elements, &
                     proj_basis, ngwf_basis, rep%ngwfs_on_grid, &
                     rep%sp_overlap, delta)
             else if (pub_aug) then
                !call aug_nonlocal_commutator_mat(nonloc_com_elements, &
                !     proj_basis, elements, dijhat, rho_ij, ngwf_basis, &
                !     rep%ngwfs_on_grid, rep%sp_overlap, delta)
             end if

             do xyz=1,3
                call sparse_scale(pos_elements_rc(xyz),-1.0_dp)
                call sparse_copy(pos_elements(xyz),pos_elements_rc(xyz))
                call sparse_axpy(pos_elements(xyz),nonloc_com_elements(xyz),(1.0_DP,0.0_DP))
             end do

          end if ! finite/cont diff

          do xyz=3,1,-1
             call sparse_destroy(nonloc_com_elements(xyz))
          end do

       else ! don't calculate the commutator

          do xyz=1,3
             call sparse_copy(pos_elements(xyz),pos_elements_rc(xyz))
          end do

       end if ! method of calculating non-local commutator

    else ! no non-local projectors so commutator not needed anyway

       do xyz=1,3
          call sparse_copy(pos_elements(xyz),pos_elements_rc(xyz))
       end do

    end if ! non-local commutator present or not

    do xyz=3,1,-1
       call sparse_destroy(pos_elements_rc(xyz))
    end do

    ! Stop timer
    call timer_clock('optics_grad_mat_els',2)


#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving optics_grad_mat_els'
#endif

  end subroutine optics_grad_mat_els


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! NB spin loop outside of this routine
  subroutine optics_dos_gp(num_energies,trans_energies,ham_type, &
       cur_spin,trans_weights,pos)

    !==================================================================!
    ! This subroutine prepares and outputs in simple txt format a      !
    ! joint density of states (JDOS) between valence and conduction    !
    ! states plot which has been generated with the Gaussian smearing  !
    ! method using gamma point only energies. The imaginary component  !
    ! of the dielectric function is also calculated for the x, y and z !
    ! directions as well as an average of the three.  This subroutine  !
    ! is based strongly on properties_dos_gp.                          !
    !------------------------------------------------------------------!
    ! Written by Laura Ratcliff in October 2011 based on the routine   !
    ! properties_dos_gp by Chris-Kriton Skylaris.                      !
    !==================================================================!

    use comms, only: comms_abort, comms_barrier, pub_on_root
    use constants, only: UP, DN, stdout, HARTREE_IN_EVS, PI
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_rootname, pub_opt_smear, pub_scissor
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_open_unit_check, utils_close_unit_check

    implicit none

    real(kind=DP), intent(in) :: trans_energies(:)
    integer, intent(in) :: num_energies, cur_spin
    character(len=6), intent(in) :: ham_type
    real(kind=DP), optional, intent(in) :: trans_weights(:)
    integer, optional, intent(in) :: pos

    ! Local Variables
    real(kind =DP) :: alpha     ! Gaussian exponent
    real(kind =DP) :: delta_e   ! energy increment
    real(kind =DP) :: e_point   ! energy point
    real(kind =DP) :: histo_val ! histogram value
    real(kind =DP) :: dist_sq   ! squared distance
    real(kind =DP) :: gnorm     ! Gaussian normalisation factor
    real(kind =DP), allocatable, dimension(:) :: energies ! buffer for energies in eV
    real(kind =DP) :: en_offset

    integer, parameter :: histonum=2001 ! number of histogram points
    integer :: row  ! row counting index
    integer :: output_unit ! fortran output unit number
    integer :: ierr ! memory allocation error flag
    integer :: io_status ! file access error flag
    integer :: orbital_index ! orbital counting index

    character(len=256) :: output_file  ! output file name buffer
    character(len=1) :: pos_string

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering optics_dos_gp'
#endif

    if (pub_on_root) then

       allocate(energies(num_energies),stat=ierr)
       call utils_alloc_check('optics_dos_gp','energies',ierr)

       ! cks: store eigen energies in buffer and convert to eV
       energies = trans_energies * HARTREE_IN_EVS
       energies = energies + pub_scissor

       ! smearing Gaussian exponent
       ! pub_opt_smear is the half-width of the Gaussian
       alpha = 0.5_dp /((pub_opt_smear*HARTREE_IN_EVS)**2)

       ! normalisation factor for Gaussians
       gnorm = sqrt(alpha/PI)
 
       en_offset = pub_opt_smear * 10.0_dp

       ! lr408: Write to different output files
       if (present(pos)) then
          write(pos_string,'(I1)') pos
       end if

       write(stdout,'(a)')' '

       write(output_file,*)trim(pub_rootname)//trim(ham_type)//'_'
       if ((pub_cell%num_spins==2).and.(cur_spin == UP)) &
            write(output_file,*) trim(output_file)//'UP'//'_'
       if ((pub_cell%num_spins==2).and.(cur_spin == DN)) &
            write(output_file,*) trim(output_file)//'DN'//'_'
       if (.not. present(pos)) then
          write(output_file,*) trim(output_file)//'JDOS.txt'
       else if (pos /= 5) then
          write(output_file,*) trim(output_file)//'OPT_SPEC'//pos_string//'.txt'
       else
          write(output_file,*) trim(output_file)//'OPT_SPEC_AV.txt'
       end if

       if (.not. present(pos)) then
          if (pub_cell%num_spins == 1) then
             write(stdout,'(a)')&
                  '================== Joint density of States (JDOS) calculation &
                  &=================='
          elseif (cur_spin == UP) then
             write(stdout,'(a)')&
                  '============ Joint density of States (JDOS) calculation &
                  &for UP spin ============='
          elseif (cur_spin == DN) then
             write(stdout,'(a)')&
                  '=========== Joint density of States (JDOS) calculation &
                  &for DOWN spin ============'
          endif
       else
          if (pub_cell%num_spins == 1) then
             write(stdout,'(a)')&
                  '================ Imaginary component of the dielectric function &
                  &================'
          elseif (cur_spin == UP) then
             write(stdout,'(a)')&
                  '========== Imaginary component of the dielectric function &
                  &for UP spin ===='
          elseif (cur_spin == DN) then
             write(stdout,'(a)')&
                  '========== Imaginary component of the dielectric function &
                  &for DOWN spin ===='
          endif
       end if

       output_file = adjustl(output_file)

       ! cks: print output warning
       write(stdout,'(3a)',advance ='no') &
            'Writing "', trim(output_file),'" ...'

       ! cks: get a unit number that is free
       output_unit = utils_unit()

       open(unit=output_unit, form="formatted" ,file=trim(output_file), &
            action="write",iostat=io_status)
       call utils_open_unit_check('optics_dos_gp','output_file', &
            io_status)

       ! cks: write first line
       if (.not. present(pos)) then
          write(output_unit,'(a)',err =100)'#  Energy (eV) |  JDOS'
       else
          write(output_unit,'(a)',err =100)'#  Energy (eV) |  Imag. diel. fn.'
       end if

       delta_e =(2.0_DP*en_offset +maxval(energies(:))&
            -minval(energies(:)))/real(histonum-1, kind=DP)

       e_point = minval(energies(:)) -en_offset

       ! cks: Loop over DOS histogram points
       do row=1, histonum

          histo_val = 0.0_DP
          ! cks: accumulate value at current histogram point
          do orbital_index=1,num_energies
             dist_sq =(e_point -energies(orbital_index))**2
             ! lr408: If appropriate, multiply by weighting factor
             if (present(trans_weights)) then
                histo_val = histo_val +gnorm*exp(-alpha *dist_sq) &
                     *trans_weights(orbital_index)
             else
                histo_val = histo_val +gnorm*exp(-alpha *dist_sq)
             end if
          end do

          ! cks: write to file current histo-point
          write(output_unit,'(f14.8,f16.10)') e_point, histo_val

          ! cks: update next histogram coordinate
          e_point = e_point +delta_e
       end do


       close(unit=output_unit,iostat=io_status)
       call utils_close_unit_check('optics_dos_gp','output_unit', &
            io_status)

       ! cks: notify of end of output
       write(stdout,*)' done'

       write(stdout,'(a)')'================================&
            &================================================'

       deallocate(energies,stat=ierr)

       call utils_dealloc_check('optics_dos_gp','energies',ierr)

    endif
    call comms_barrier

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving optics_dos_gp'
#endif

    return

100 write(stdout,*)'Problem writing to file in optics_dos_gp. ONETEP stops'
    call comms_abort

  end subroutine optics_dos_gp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module optics
