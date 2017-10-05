! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                  Energy Eigenstate module                      !
!                                                                !
! This module analyses the energy eigenstates obtained by        !
! diagonalising the hamiltonian in the NGWF basis, after a       !
! converged ONETEP calculation has completed.                    !
!----------------------------------------------------------------!
! Module created from parts of properties_mod.F90 by Nicholas    !
! Hine on 16/11/2011.                                            !
! The original routines were written by Chris-Kriton Skylaris,   !
! Peter Haynes, Nicholas Hine, and others, 2002-2011.            !
!================================================================!

module eigenstates

  use constants, only: DP

  implicit none

  private

  public :: eigenstates_calculate
  public :: eigenstates_plot_mo

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CW
!proj_basis used as optional argument now
  subroutine eigenstates_calculate(denskern, ham, rep, &
       ngwf_basis, hub_proj_basis, elements, ham_type, &
       num_opt_states,proj_basis,eigen_en_copy)
!END CW

    !======================================================================!
    ! This subroutine generates and diagonalises a Hamiltonian and         !
    ! performs various kinds of analysis with the resulting eigenvectors   !
    ! and eigenenergies.                                                   !
    !----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 10/05/2006 as part of the routine!
    ! properties_calculate.                                                !
    ! Modified by Nicholas Hine in July 2009 to use function_basis and     !
    ! SPAM3 routines.                                                      !
    ! Modified by Nicholas Hine in October 2009 to make use of ScaLAPACK   !
    ! for diagonalisation of large systems.                                !
    ! Modified by Nicholas Hine in February 2010 to remove #ifdef's by     !
    ! moving code relating to eigensolving to dense_mod.                   !
    ! Modified by Nicholas Hine in October 2010 to use NGWF_REP.           !
    ! Modified by Laura Ratcliff in October 2010 for conduction            !
    ! calculations.                                                        !
    ! Simplified by Nicholas Hine in April 2011 to reduce branching of     !
    ! code in spectra calculations.                                        !
    ! Moved to eigenstates_mod by Nicholas Hine in November 2011.          !
    !======================================================================!

    use augmentation, only: augmentation_overlap
    use cell_grid, only: pub_fine_grid
    use comms, only: pub_on_root, comms_abort
    use constants, only: max_spins, stdout, UP, DN, HARTREE_IN_EVS
    use dense, only: DEM, dense_create, dense_destroy, dense_eigensolve, &
         dense_product, dense_get_col, dense_convert, dense_get_element
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use optics, only: optics_calculate_spectra
    use rundat, only: pub_homo_dens_plot, pub_lumo_dens_plot, pub_homo_plot, &
         pub_lumo_plot, pub_dos_smear, pub_ldos_smear, pub_ldos_ngroups, &
         pub_num_eigenvalues, pub_aug, pub_cond_calculate, &
         pub_spectra_calculate, pub_calc_mom_mat_els, pub_hubbard
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_axpy, &
         sparse_destroy, sparse_copy, sparse_product, &
         sparse_scale, sparse_transpose, sparse_transpose_structure
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in) :: rep
    type(NGWF_HAM), intent(inout) :: ham
    type(SPAM3), intent(inout) :: denskern(pub_cell%num_spins)
    type(FUNC_BASIS), intent(in) :: ngwf_basis

!CW
    !type(FUNC_BASIS), intent(in) :: proj_basis
     type(FUNC_BASIS), optional   :: proj_basis
!END CW
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(ELEMENT), intent(in) :: elements(:) ! elements of all atoms in  input file order
    ! lr408: Optional string defining which type of Hamiltonian this is for
    character(len=*), intent(in) :: ham_type
    ! lr408: Optional integer defining the total number of states to include for optical spectra
    integer, optional, intent(in) :: num_opt_states

    ! Local Variables
    type(SPAM3)  :: mo_kern(1)   ! density kernel of one molecular orbital
    type(SPAM3)  :: hub_proj     ! Hubbard projection in SPAM3 format
    type(SPAM3)  :: aug_overlap  ! Augmentation part of overlap matrix
    type(SPAM3)  :: cross_olap_trans ! temporary cross overlap transpose
    type(SPAM3)  :: kt,ukt       ! Density kernel projected into current rep
    type(DEM) :: eigs_dens       ! dense matrix for eigenvectors
    type(DEM) :: buffer_dens     ! dense matrix for hamiltonian and density kernel
    type(DEM) :: overlap_dens    ! dense matrix for overlap
    type(DEM) :: hub_proj_dens   ! dense matrix for Hubbard projection
    real(kind=DP), allocatable, dimension(:,:) :: eigen_en     ! hamiltonian eigenvalues
!CW
    real(kind=DP),optional                     :: eigen_en_copy(:,:)
!END CW
    real(kind=DP), allocatable, dimension(:,:) :: eigen_occ    ! denskern occupancies
    real(kind=DP), allocatable, dimension(:) :: mo_coeff    ! coefficients of MO's
    integer :: num                 ! number of NGWFs to diagonalise over (may be joint)
    integer :: ierr                ! memory allocation error flag
    integer :: is                  ! spin loop counter
    integer :: orbital_index       ! orbital counting index
    integer :: homo_max            ! maximum index of occupied orbital to plot
    integer :: homo_min            ! minimum index of occupied orbital to plot
    integer :: lumo_max            ! maximum index of virtual orbital to plot
    integer :: lumo_min            ! minimum index of virtual orbital to plot
    character(len=256) :: output_file  ! file names
    character(len=256) :: title_line   ! title line in output file

    ! Start timer
    call timer_clock("eigenstates_calculate",1)

    ! lr408: Print diagonalisation headers
    if (pub_cond_calculate) then
       if (pub_on_root) then
          if (ham_type == 'valence') then
             write(stdout,'(/a)') ''
             write(stdout,'(/a)') '=========================== &
                  &Valence diagonalisation &
                  &============================'
          else if (ham_type == 'proj') then
             write(stdout,'(/a)') ''
             write(stdout,'(/a)') '==================== &
                  &Projected conduction diagonalisation &
                  &======================'
          else if (ham_type == 'cond') then
             write(stdout,'(/a)') ''
             write(stdout,'(/a)') '=================== &
                  &Unprojected conduction diagonalisation &
                  &====================='
          else if (ham_type == 'joint') then
             write(stdout,'(/a)') ''
             write(stdout,'(/a)') '=============== &
                  &Joint valence and conduction diagonalisation &
                  &=================='
          end if
       end if
    end if

    ! ndmh: Allocate storage for eigenstates
    num = ngwf_basis%num
    allocate(eigen_en(num,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_calculate','eigen_en',ierr)
    allocate(eigen_occ(num,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('eigenstates_calculate','eigen_occ',ierr)
    allocate(mo_coeff(num),stat=ierr)
    call utils_alloc_check('eigenstates_calculate','mo_coeff',ierr)

    ! Create dense matrices for Hamiltonian and overlap
    call dense_create(buffer_dens,num,num)
    call dense_create(overlap_dens,num,num)

    ! No need for eigenvectors unless they are needed for plotting
    if ((pub_homo_dens_plot >=0 .or. pub_lumo_dens_plot >= 0 .or. &
         pub_homo_plot >= 0 .or. pub_lumo_plot >=0) &
         .or. ((pub_ldos_ngroups > 0).and.(pub_ldos_smear > 0.0_DP))) then
       call dense_create(eigs_dens,num,num)
    else
       call dense_create(eigs_dens,num,num)
    end if

    ! ndmh: Create matrix for <chi_c|phi_b> K^ba <phi_a|chi_d>
    if ((ham_type=='valence').or.(ham_type=='cond').or. &
         (ham_type == 'proj')) then
       call sparse_create(kt,denskern(1),rep%overlap)
       call sparse_create(ukt,rep%overlap,kt)
    else if (ham_type=='joint') then
       call sparse_create(kt,denskern(1),rep%cross_overlap)
       call sparse_transpose_structure(cross_olap_trans%structure, &
            rep%cross_overlap)
       call sparse_create(cross_olap_trans)
       call sparse_transpose(cross_olap_trans,rep%cross_overlap)
       call sparse_create(ukt,cross_olap_trans,kt)
    end if

    ! ddor: Create matrices for DFT+U projection
    if (pub_hubbard) then
       call dense_create(hub_proj_dens,num,num)
       call sparse_create(hub_proj,rep%hub_overlap,rep%hub_overlap_t)
    endif

    ! ndmh: Create matrix to hold MO kernel
    if ((pub_aug.and.(pub_homo_plot >= 0 .or. pub_lumo_plot >= 0)).or. &
         (pub_homo_dens_plot >= 0 .or. pub_lumo_dens_plot >= 0)) then
       call sparse_create(mo_kern(1),denskern(1))
    end if

    ! ndmh: Create matrix for aug part of overlap
    if (pub_aug.and.(pub_homo_plot >= 0 .or. pub_lumo_plot >= 0)) then
       call sparse_create(aug_overlap,rep%overlap)
       call augmentation_overlap(aug_overlap,rep%sp_overlap)
    end if

    ! pdh: loop over spins
    do is=1,pub_cell%num_spins

       ! ndmh: diagonalise hamiltonian for eigenvalues
       call dense_convert(overlap_dens,rep%overlap)
       ! lr408: Choose correct Hamiltonian for conduction calculation
       if (ham_type=='cond') then
          call dense_convert(buffer_dens,ham%cond_non_proj_ham(is))
       else
          call dense_convert(buffer_dens,ham%ham(is))
       end if
       call dense_eigensolve(num,eigen_en(:,is),buffer_dens, &
            overlap_dens,1,eigs_dens)

       ! ndmh: now calculate occupancy of each state by projection onto
       ! ndmh: valence density matrix for this spin
       if ((ham_type=='valence').or.(ham_type=='cond').or. &
            (ham_type == 'proj')) then
          call sparse_product(kt,denskern(is),rep%overlap)
          call sparse_product(ukt,rep%overlap,kt)
       else if (ham_type=='joint') then
          call sparse_product(kt,denskern(is),rep%cross_overlap)
          call sparse_product(ukt,cross_olap_trans,kt)
       end if
       call dense_convert(overlap_dens,ukt)
       call dense_product(buffer_dens,eigs_dens,overlap_dens, &
            transpose_amat=.true.)
       call dense_product(overlap_dens,buffer_dens,eigs_dens)
       do orbital_index=1,ngwf_basis%num
          call dense_get_element(eigen_occ(orbital_index,is),overlap_dens, &
               orbital_index,orbital_index)
       end do

       ! ndmh: find range of states to plot
       if (pub_cond_calculate .and. ham_type == 'proj') then
          lumo_min = 1
          lumo_max = 1  + pub_lumo_plot
       else
          lumo_min = rep%n_occ(is) + 1
          lumo_max = rep%n_occ(is) + 1 + pub_lumo_plot
       end if
       if (lumo_max > ngwf_basis%num) lumo_max = ngwf_basis%num
       homo_max = rep%n_occ(is)
       homo_min = rep%n_occ(is) - pub_homo_plot
       if (homo_min <= 0) homo_min = 1

       ! cks: Create canonical orbital plot files
       call internal_plot_orbitals

       ! ndmh: Calculate and output LDOS for this spin
       if ((pub_ldos_smear > 0.0_DP).and.(pub_ldos_ngroups>0)) then

          ! ndmh: Standard LDOS, projected onto atoms
          call dense_convert(overlap_dens,rep%overlap)
          call dense_product(buffer_dens,overlap_dens,eigs_dens)
          call eigenstates_ldos_gp(eigen_en(:,is),eigs_dens,buffer_dens, &
               is,elements,ngwf_basis,.false.,ham_type)

          ! ddor: Compute projection of LDOS onto Hubbard manifold
          if (pub_hubbard) then
             call sparse_product(hub_proj,&
                  rep%hub_overlap,rep%hub_overlap_t)
             call dense_convert(hub_proj_dens,hub_proj)
             call dense_product(buffer_dens,hub_proj_dens,eigs_dens)
             call eigenstates_ldos_gp(eigen_en(:,is),eigs_dens,buffer_dens, &
                  is,elements,ngwf_basis,.true.,ham_type)
          endif

       end if

       ! lr408: Calculate matrix elements for spectra - only in conduction
       ! lr408: calculation and for valence and joint basis sets only
       if (pub_spectra_calculate .and. pub_cond_calculate) then
!CW
       if(.not.present(proj_basis))then
        write(*,*) ' proj basis argument not present, stop'
        stop
       endif 
!END CW 
         if (ham_type=='valence'.or.ham_type=='joint') then
             call optics_calculate_spectra(is,eigs_dens,eigen_en(:,is), &
                  rep,ngwf_basis,proj_basis,num_opt_states)
          end if
       end if

    end do

    ! cks: output gamma-point Gaussian-smeared density of states
    ! cks: and CASTEP format .bands file
    if (pub_dos_smear > 0.0_DP) then
       call eigenstates_dos_gp(eigen_en, rep%n_occ, ngwf_basis,ham_type)
       call eigenstates_write_bands(ngwf_basis, eigen_en, rep%n_occ,ham_type)
    end if

    ! cks: print out orbital energies and unpurified denskern occupancies
    if (pub_num_eigenvalues > 0 .or. pub_cond_calculate) then
       call eigenstates_print_ens_occs(eigen_en, eigen_occ, rep%n_occ, &
            ngwf_basis%num)
    end if

    ! ndmh: Destroy matrix for aug part of overlap
    if (pub_aug.and.(pub_homo_plot >= 0 .or. pub_lumo_plot >= 0)) then
       call sparse_destroy(aug_overlap)
    end if

    ! ndmh: Destroy matrix for MO kernel
    if ((pub_aug.and.(pub_homo_plot >= 0 .or. pub_lumo_plot >= 0)).or. &
         (pub_homo_dens_plot >= 0 .or. pub_lumo_dens_plot >= 0)) then
       call sparse_destroy(mo_kern(1))
    end if

    ! ndmh: Destroy matrices for kernel projection
    call sparse_destroy(ukt)
    if (ham_type=='joint') call sparse_destroy(cross_olap_trans)
    call sparse_destroy(kt)

    ! Destroy dense matrices for Hamiltonian, denskern and overlap
    call dense_destroy(eigs_dens)
    if (pub_hubbard) then
       call sparse_destroy(hub_proj)
       call dense_destroy(hub_proj_dens)
    endif
    call dense_destroy(overlap_dens)
    call dense_destroy(buffer_dens)

!CW
    eigen_en_copy=eigen_en
!END CW

    deallocate(mo_coeff,stat=ierr)
    call utils_dealloc_check('eigenstates_calculate','mo_coeff',ierr)
    deallocate(eigen_occ,stat=ierr)
    call utils_dealloc_check('eigenstates_calculate','eigen_occ',ierr)
    deallocate(eigen_en, stat=ierr)
    call utils_dealloc_check('eigenstates_calculate','eigen_en',ierr)

    ! cks: stop timer
    call timer_clock("eigenstates_calculate",2)

    contains

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine internal_plot_orbitals

        if ((pub_homo_plot >= 0 .or. pub_lumo_plot >= 0)) then

           if (pub_on_root) then
              write(stdout,'(/a)') '===================== &
                   &Writing canonical molecular orbitals &
                   &====================='
              if (pub_cell%num_spins == 2) write(stdout,'(a,i1)') &
                   '    Spin ',is
              write(stdout,'(a)')'| Orbital  Energy (eV)  Norm     |'
           end if

           ! ndmh: output homo and below for plotting
           do orbital_index=homo_min,homo_max
              call internal_orb_title_and_filename('HOMO',rep%n_occ(is))

              ! ndmh: if augmentation is active, we need to calculate the
              ! ndmh: aug part of the norm, so we need the MO 'kernel'
              if (pub_aug) then
                 ! ndmh: create dense 'kernel' for MO from selected eigenvector
                 call dense_product(buffer_dens,eigs_dens,eigs_dens, &
                      transpose_amat=.false.,transpose_bmat=.true., &
                      first_k=orbital_index,last_k=orbital_index)
                 ! ndmh: convert to SPAM3
                 call dense_convert(mo_kern(1),buffer_dens)
              end if

              call dense_get_col(mo_coeff,eigs_dens,orbital_index)
              call eigenstates_plot_mo(mo_coeff, rep%ngwfs_on_grid, &
                   ngwf_basis, elements, aug_overlap, mo_kern(1), &
                   title_line, output_file)
           end do

           ! ndmh: output lumo and above for plotting
           do orbital_index=lumo_min,lumo_max
              if (ham_type == 'proj') then
                 call internal_orb_title_and_filename('LUMO',1)
              else
                 call internal_orb_title_and_filename('LUMO',rep%n_occ(is)+1)
              end if

              ! ndmh: if augmentation is active, we need to calculate the
              ! ndmh: aug part of the norm, so we need the MO 'kernel'
              if (pub_aug) then
                 ! ndmh: create dense 'kernel' for MO from selected eigenvector
                 call dense_product(buffer_dens,eigs_dens,eigs_dens, &
                      transpose_amat=.false.,transpose_bmat=.true., &
                      first_k=orbital_index,last_k=orbital_index)
                 ! ndmh: convert to SPAM3
                 call dense_convert(mo_kern(1),buffer_dens)
              end if

              call dense_get_col(mo_coeff,eigs_dens,orbital_index)
              call eigenstates_plot_mo(mo_coeff, rep%ngwfs_on_grid, &
                   ngwf_basis, elements, aug_overlap, mo_kern(1), &
                   title_line, output_file)
           end do

           if (pub_on_root) write(stdout,'(a)')&
                '========================================================&
                &========================'
        end if

        ! ndmh: find range of states to plot
        if (pub_cond_calculate .and. ham_type=='proj') then
           lumo_min = 1
           lumo_max = 1 + pub_lumo_dens_plot
        else
           lumo_min = rep%n_occ(is) + 1
           lumo_max = rep%n_occ(is) + 1 + pub_lumo_dens_plot
        end if
        if (lumo_max > ngwf_basis%num) lumo_max = ngwf_basis%num
        homo_max = rep%n_occ(is)
        homo_min = rep%n_occ(is) - pub_homo_dens_plot
        if (homo_min <= 0) homo_min = 1

        if ((pub_homo_dens_plot >= 0 .or. pub_lumo_dens_plot >= 0)) then

           if (pub_on_root) then
              write(stdout,'(/a)') '============== &
                   &Writing densities of canonical molecular orbitals &
                   &==============='
              if (pub_cell%num_spins == 2) write(stdout,'(a,i1)') &
                   '    Spin ',is
              write(stdout,'(a)')'| Orbital  Energy (eV)  Norm     |'
           end if

           ! cks: output homo and below for plotting
           do orbital_index=homo_min,homo_max
              ! ndmh: create dense 'kernel' for MO from selected eigenvector
              call dense_product(buffer_dens,eigs_dens,eigs_dens, &
                   transpose_amat=.false.,transpose_bmat=.true., &
                   first_k=orbital_index,last_k=orbital_index)
              ! ndmh: convert to SPAM3
              call dense_convert(mo_kern(1),buffer_dens)
              ! ndmh: plot orbital density scalar field
              call internal_orb_title_and_filename('HOMO_density',rep%n_occ(is))
              call eigenstates_plot_squared_mo( &
                   mo_kern(1), rep%ngwf_overlap, rep%ngwfs_on_grid, &
                   ngwf_basis, elements, title_line, output_file)
           end do

           ! cks: output lumo and above for plotting
           if (lumo_max > ngwf_basis%num) lumo_max = ngwf_basis%num
           do orbital_index=lumo_min,lumo_max
              ! ndmh: create dense 'kernel' for MO from selected eigenvector
              call dense_product(buffer_dens,eigs_dens,eigs_dens, &
                   transpose_amat=.false.,transpose_bmat=.true., &
                   first_k=orbital_index,last_k=orbital_index)
              ! ndmh: convert to SPAM3
              call dense_convert(mo_kern(1),buffer_dens)
              ! ndmh: plot orbital density scalar field
              if (ham_type == 'proj') then
                 call internal_orb_title_and_filename('LUMO_density',1)
              else
                 call internal_orb_title_and_filename('LUMO_density',rep%n_occ(is)+1)
              end if
              call eigenstates_plot_squared_mo( &
                   mo_kern(1), rep%ngwf_overlap, rep%ngwfs_on_grid, &
                   ngwf_basis, elements, title_line, output_file)
           end do

           if (pub_on_root) write(stdout,'(a)')&
                '========================================================&
                &========================'
        end if

      end subroutine internal_plot_orbitals

      subroutine internal_orb_title_and_filename(id,origin)

        ! Arguments
        character(*),intent(in) :: id
        integer,intent(in) :: origin

        ! Local Variables
        character(len=50)  :: txt_buffer   ! text buffer
        character(len=50)  :: orb_buffer   ! text buffer
        character(len=256) :: title_buffer ! text buffer

        write(title_buffer,'(a30,i4)')'Density of canonical orbital:', &
             orbital_index

        if (orbital_index == origin) then
           if (pub_cell%num_spins == 1) then
              write(output_file,*) '_',id
           else
              if (is == UP) then
                 write(output_file,*) '_',id,'_UP'
              else
                 write(output_file,*) '_',id,'_DN'
              end if
           end if
           if (pub_cond_calculate) write(output_file,*) &
                '_'//trim(ham_type)//trim(adjustl(output_file))
           if (pub_on_root) write(stdout,'(a8,f14.8,a)',advance ='no') &
                adjustl('| '//id//'  '),&
                eigen_en(orbital_index,is)*HARTREE_IN_EVS,' '
           write(title_line,*)trim(adjustl(title_buffer)),' (',id,')'
        else
           if ((orbital_index-origin)>0) then
              write(txt_buffer,*) orbital_index-origin
              write(orb_buffer,*) '+',trim(adjustl(txt_buffer))
           else
              write(orb_buffer,*) orbital_index-origin
           end if
           if (pub_cell%num_spins == 1) then
              write(output_file,*)'_',id,trim(adjustl(orb_buffer))
           else
              if (is == UP) then
                 write(output_file,*)'_',id, &
                      trim(adjustl(orb_buffer)),'_UP'
              else
                 write(output_file,*)'_',id, &
                      trim(adjustl(orb_buffer)),'_DN'
              end if
           end if
           if (pub_cond_calculate) write(output_file,*) &
                '_'//trim(ham_type)//trim(adjustl(output_file))
           write(title_line,*)trim(adjustl(title_buffer)),' (',id, &
                trim(adjustl(orb_buffer)),')'
           write(txt_buffer,*) '| ',id,trim(adjustl(orb_buffer))
           if (pub_on_root) write(stdout,'(a9,f13.8,a)',advance ='no')&
                adjustl(txt_buffer),eigen_en(orbital_index,is)* &
                HARTREE_IN_EVS,' '
        end if

      end subroutine internal_orb_title_and_filename

  end subroutine eigenstates_calculate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_plot_mo(mo_coeff, &
     ngwfs_on_grid, ngwf_basis, elements, aug_overlap, mo_kern, &
     file_header, scalar_name)

    !==================================================================!
    ! This subroutine prepares and outputs a plotfile for the value of !
    ! a canonical molecular orbital.                                   !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 13/01/2011.                          !
    !==================================================================!

    use cell_grid, only: pub_std_grid
    use comms, only: pub_on_root
    use constants, only: stdout
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_product_on_grid
    use ion, only: ELEMENT
    use rundat, only: pub_cube_format, pub_dx_format, pub_grd_format, pub_aug
    use sparse, only: SPAM3, sparse_trace
    use utils, only: utils_alloc_check, utils_dealloc_check
    use visual, only: visual_scalarfield

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: mo_coeff(ngwf_basis%num) ! NGWF coeffs of orb.
    type(ELEMENT), intent(in) :: elements(:) ! all elements in input file order
    real(kind=DP), intent(in) :: ngwfs_on_grid(:) ! ngwfs on this proc
    type(SPAM3), intent(in) :: aug_overlap    ! sp overlaps
    type(SPAM3), intent(in) :: mo_kern       ! for aug norm
    character(len=*), intent(in) :: file_header ! part of header line of output file
    character(len=*), intent(in) :: scalar_name ! part of file-name to output (scalafield name)

    ! Local Variables
    integer :: ierr
    real(kind=DP) :: norm, aug_norm
    real(kind=DP),allocatable :: orbital_std(:,:,:)

    ! ndmh: allocate storage for orbital
    allocate(orbital_std(pub_std_grid%ld1,pub_std_grid%ld2, &
         pub_std_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('eigenstates_plot_mo','orbital_std',ierr)

    ! lr408: set to zero here instead to allow joint val-cond basis calculations
    orbital_std = 0.0_DP

    ! ndmh: calculate the M.O. on the standard grid
    call eigenstates_calculate_mo(orbital_std, mo_coeff, ngwfs_on_grid, &
         ngwf_basis)

    ! ndmh: get norm of M.O.
    norm = integrals_product_on_grid(pub_std_grid,orbital_std,orbital_std)

    ! ndmh: augment norm with augmentation part in PAW/USP formalism
    if (pub_aug) then
       aug_norm = sparse_trace(mo_kern,aug_overlap)
       norm = norm + aug_norm
    end if

    if (pub_on_root) write(stdout,'(f9.6,a)',advance='no') norm, ' | '

    ! ndmh: plot this orbital on standard grid
    call visual_scalarfield(orbital_std(:,:,:), pub_std_grid, &
         file_header, scalar_name, elements, 1.0_DP)  ! input

    ! ndmh: print end of line if no message was written
    if ((.not.pub_grd_format).and.(.not.pub_cube_format).and. &
         (.not.pub_dx_format).and.pub_on_root) write(stdout,*)

    ! ndmh: deallocate storage for orbital
    deallocate(orbital_std,stat=ierr)
    call utils_dealloc_check('eigenstates_plot_mo','orbital_std',ierr)

  end subroutine eigenstates_plot_mo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_calculate_mo(orbital_std, mo_coeff, ngwfs_on_grid, &
       ngwf_basis)

    !==================================================================!
    ! This subroutine calculates a canonical molecular orbital from    !
    ! coefficients of each NGWF, returned by a diagonalisation of the  !
    ! Hamiltonian.                                                     !
    !------------------------------------------------------------------!
    ! Written by Nicholas Hine on 13/01/2011.                          !
    !==================================================================!

    use basis, only: basis_copy_function_to_box, &
         basis_location_func_wrt_cell
    use cell_grid, only: pub_std_grid, cell_grid_deposit_box
    use comms, only: pub_my_node_id
    use function_basis, only: FUNC_BASIS
    use utils, only: utils_alloc_check, utils_dealloc_check

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP),intent(in) :: mo_coeff(ngwf_basis%num) ! NGWF coeffs of orb.
    real(kind=DP), intent(in) :: ngwfs_on_grid(:)        ! NGWFs on this proc
    real(kind=DP), intent(inout) :: orbital_std(pub_std_grid%ld1, &
         pub_std_grid%ld2,pub_std_grid%max_slabs12)

    ! Local Variables
    integer :: ierr
    integer :: ingwf, loc_ingwf
    integer :: ngwf_cell_start1, ngwf_cell_start2, ngwf_cell_start3
    logical :: i_have_box
    integer :: n1, n2, n3
    real(kind=DP),allocatable :: ngwf_box_std(:,:,:)
    real(kind=DP),allocatable :: buffer_std(:,:,:)

    ! ndmh: initialisations
    n1 = ngwf_basis%maxtight_pts1
    n2 = ngwf_basis%maxtight_pts2
    n3 = ngwf_basis%maxtight_pts3

    ! ndmh: allocate storage for buffer and NGWF tightbox
    allocate(buffer_std(n1,n2,pub_std_grid%max_slabs12),stat=ierr)
    call utils_alloc_check('eigenstates_calculate_mo','buffer_std',ierr)
    allocate(ngwf_box_std(n1,n2,n3),stat=ierr)
    call utils_alloc_check('eigenstates_calculate_mo','ngwf_box_std',ierr)
    ngwf_box_std = 0.0_DP
    buffer_std = 0.0_DP

    ! ndmh: Loop over NGWFs on each node
    do loc_ingwf=1,ngwf_basis%max_on_node
       ingwf = loc_ingwf + ngwf_basis%first_on_node(pub_my_node_id) - 1

       if (loc_ingwf<=ngwf_basis%num_on_node(pub_my_node_id)) then

          ! ndmh: copy NGWF to tightbox
          call basis_copy_function_to_box(ngwf_box_std,n1,n2,n3, &
               1,1,1,ngwf_basis%tight_boxes(loc_ingwf), &
               ngwfs_on_grid,ngwf_basis%spheres(loc_ingwf))
          ngwf_box_std = ngwf_box_std * mo_coeff(ingwf)

          ! ndmh: find start of tightbox of loc_ingwf in simulation cell
          call basis_location_func_wrt_cell( &
               ngwf_cell_start1, ngwf_cell_start2, ngwf_cell_start3, &
               ngwf_basis%tight_boxes(loc_ingwf))

          i_have_box = .true.
       else
          ngwf_cell_start1 = -1234
          ngwf_cell_start2 = -1234
          ngwf_cell_start3 = -1234
          i_have_box = .false.
       end if

       ! ndmh: deposit NGWF scaled by orbital coefficient to standard grid
       call cell_grid_deposit_box(orbital_std,ngwf_box_std,buffer_std, &
            pub_std_grid,n1,n2,n3,n1,n2,ngwf_cell_start1,ngwf_cell_start2, &
            ngwf_cell_start3,i_have_box,.false.)
    end do

    ! ndmh: deallocate universal tightbox array
    deallocate(ngwf_box_std,stat=ierr)
    call utils_dealloc_check('eigenstates_calculate_mo','ngwf_box_std',ierr)
    deallocate(buffer_std,stat=ierr)
    call utils_dealloc_check('eigenstates_calculate_mo','buffer_std',ierr)

  end subroutine eigenstates_calculate_mo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_plot_squared_mo( mo_kern, &
     overlap, ngwfs_on_grid, ngwf_basis, elements, &
     file_header, scalar_name)

    !==================================================================!
    ! This subroutine prepares and outputs a plotfile for the square   !
    ! of a selected canonical molecular orbital.                       !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 10/05/2006.                  !
    ! Modified by Nicholas Hine on 04/03/2010 to take in mo_kern       !
    ! already formed, rather than create it from the orbitals.         !
    ! Trivially modified by Jacek Dziedzic on 14/05/2010 to delegate   !
    ! the unit conversion to visual_scalarfield.                       !
    !==================================================================!

    use cell_grid, only: pub_fine_grid
    use comms, only: comms_abort, pub_my_node_id, pub_on_root
    use constants, only: DP, ANGSTROM, stdout
    use density, only: density_on_grid
    use function_basis, only: FUNC_BASIS
    use integrals, only: integrals_trace_on_grid
    use ion, only: ELEMENT
    use rundat, only: pub_cube_format, pub_dx_format, pub_grd_format
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3
    use utils, only: utils_alloc_check, utils_dealloc_check
    use visual, only: visual_scalarfield

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(SPAM3), intent(in) :: mo_kern(1) ! denskern of single mo as SPAM3
    type(SPAM3), intent(in)    :: overlap    ! overlap as SPAM3
    type(ELEMENT), intent(in) :: elements(:) ! all elements in input file order
    real(kind =DP), intent(in) :: ngwfs_on_grid(:) ! ngwfs on this proc
    character(len =*), intent(in) :: file_header ! part of header line of output file
    character(len =*), intent(in) :: scalar_name ! part of file-name to output (scalafield name)

    ! Local Variables
    real(kind=DP), allocatable, dimension(:,:,:,:) :: mo_density_fine
    real(kind=DP) :: norm
    integer :: ierr        ! memory allocation error flag
    integer :: nspins

    ! cks: allocate memory for orbital charge density
    allocate(mo_density_fine(pub_fine_grid%ld1, &
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, 1),stat=ierr)
    call utils_alloc_check('eigenstates_plot_squared_mo','mo_density_fine',ierr)

    ! cks: generate MO "density"
    ! pdh: nasty hack
    nspins = pub_cell%num_spins
    pub_cell%num_spins = 1
    call density_on_grid(mo_density_fine, pub_fine_grid, &          ! output
         mo_kern(1), overlap, ngwfs_on_grid, ngwf_basis)            ! input
    pub_cell%num_spins = nspins

    ! ndmh: get sum of M.O. density
    norm = integrals_trace_on_grid(mo_density_fine(:,:,:,1),pub_fine_grid)
    if (pub_on_root) write(stdout,'(f9.6,a)',advance='no') norm, ' | '

    ! cks: output MO density to plotfile formats
    ! cks: output orbital density in Angstrom^-3 rather than in Bohr^-3
    ! jd:  ... but leave the unit conversion to visual_scalarfield
    call visual_scalarfield( &
         mo_density_fine(:,:,:,1), pub_fine_grid, file_header, &
         scalar_name, elements, ANGSTROM**3)  ! input

    ! ndmh: print end of line if no message was written
    if ((.not.pub_grd_format).and.(.not.pub_cube_format).and. &
         (.not.pub_dx_format).and.pub_on_root) write(stdout,*)

    ! cks: deallocate memory for orbital charge density
    deallocate(mo_density_fine,stat=ierr)
    call utils_dealloc_check('eigenstates_plot_squared_mo', &
         'mo_density_fine',ierr)

  end subroutine eigenstates_plot_squared_mo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_dos_gp(eigen_en, n_occ, ngwf_basis, ham_type)

    !==================================================================!
    ! This subroutine prepares and outputs in simple txt format a      !
    ! density of states (DOS) plot which has been generated with the   !
    ! Gaussian smearing method using gamma point only energies. The    !
    ! DOS is calculated as a histogram of the sum of all Gaussian-     !
    ! broadened energies. Energies in the output are in eV and so is   !
    ! the half-width of the smearing Gaussians which is set by the     !
    ! input parameter dos_smear.                                       !
    !------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 10/05/2006.                  !
    ! Moved to eigenstates_mod by Nicholas Hine in November 2011.      !
    !==================================================================!

    use comms, only: comms_abort, comms_barrier, pub_on_root
    use constants, only: DP, UP, DN, stdout, HARTREE_IN_EVS, PI
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_rootname, pub_dos_smear, pub_cond_calculate
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_open_unit_check, utils_close_unit_check

    implicit none

    real(kind=DP), intent(in) :: eigen_en(:,:) ! molecular orbital energies
    integer, intent(in) :: n_occ(:)            ! number of occupied orbitals
    type(FUNC_BASIS), intent(in) :: ngwf_basis ! Function basis for NGWFs
    character(len=*), intent(in) :: ham_type

    ! Local Variables
    real(kind =DP) :: alpha     ! Gaussian exponent
    real(kind =DP) :: fermi_e   ! Fermi energy
    real(kind =DP) :: delta_e   ! energy increment
    real(kind =DP) :: e_point   ! energy point
    real(kind =DP) :: histo_val ! histogram value
    real(kind =DP) :: dist_sq   ! squared distance
    real(kind =DP) :: gnorm     ! Gaussian normalisation factor
    real(kind =DP), allocatable, dimension(:,:) :: energies ! buffer for energies in eV
    real(kind =DP), parameter :: en_offset =3.0_DP ! left and right offset for energy scale
    integer, parameter :: histonum=2001 ! number of histogram points
    integer :: row  ! row counting index
    integer :: is   ! spin loop counter
    integer :: output_unit ! fortran output unit number
    integer :: ierr ! memory allocation error flag
    integer :: io_status ! file access error flag
    integer :: orbital_index ! orbital counting index
    character(len=256) :: output_file  ! output file name buffer


    if (pub_on_root) then


       allocate(energies(ngwf_basis%num,pub_cell%num_spins),stat=ierr)
       call utils_alloc_check('eigenstates_dos_gp','energies',ierr)

       ! cks: store eigen energies in buffer and convert to eV
       energies = eigen_en
       energies = energies*HARTREE_IN_EVS


       fermi_e = 0.0_DP
       do is=1,pub_cell%num_spins
          ! cks: assume Fermi energy as average of HOMO and LUMO and spins
          if (n_occ(is)<ngwf_basis%num.and.n_occ(is)>0) then
             fermi_e = fermi_e + 0.5_DP*(energies(n_occ(is)+1,is)+ &
                  energies(n_occ(is),is))
          ! ndmh: unless we have all orbitals occupied
          else if (n_occ(is) > 0) then
             fermi_e = fermi_e + energies(n_occ(is),is) + tiny(1.0_DP)
          else
             fermi_e = fermi_e + energies(1,is) + tiny(1.0_DP)
          end if
       end do
       fermi_e = fermi_e / pub_cell%num_spins

       ! cks: set fermi energy as the zero of energy
       energies = energies - fermi_e

       ! cks: smearing Gaussian exponent
       ! cks: pub_dos_smear is the half-width of the Gaussian
       alpha =log(2.0_DP) /((pub_dos_smear*HARTREE_IN_EVS)**2)

       ! cks: normalisation factor for Gaussians
       gnorm = sqrt(alpha/PI)


       ! cks: generate DOS for each spin
       spin_loop: do is=1,pub_cell%num_spins

          if (pub_on_root) write(stdout,'(a)')' '

          write(output_file,*)trim(pub_rootname)//'_'
          if (pub_cond_calculate) &
               write(output_file,*)trim(output_file)//trim(ham_type)//'_'
          if ((pub_cell%num_spins==2).and.(is == UP)) &
               write(output_file,*) trim(output_file)//'UP'//'_'
          if ((pub_cell%num_spins==2).and.(is == DN)) &
               write(output_file,*) trim(output_file)//'DN'//'_'
          write(output_file,*) trim(output_file)//'DOS.txt'

          if (pub_cell%num_spins == 1) then
             write(stdout,'(a)')&
                  '===================== Density of States (DOS) calculation &
                  &======================'
          elseif (is == UP) then
             write(stdout,'(a)')&
                  '=============== Density of States (DOS) calculation &
                  &for UP spin ================'
          elseif (is == DN) then
             write(stdout,'(a)')&
                  '============== Density of States (DOS) calculation &
                  &for DOWN spin ==============='
          endif

          output_file = adjustl(output_file)

          ! cks: print output warning
          write(stdout,'(3a)',advance ='no') &
               'Writing "', trim(output_file),'" ...'

          ! cks: get a unit number that is free
          output_unit = utils_unit()

          open(unit=output_unit, form="formatted" ,file=trim(output_file), &
               action="write",iostat=io_status)
          call utils_open_unit_check('eigenstates_dos_gp','output_file', &
               io_status)

          ! cks: write first line
          write(output_unit,'(a)',err =100)'#  Energy (eV) |  DOS (states/eV)'

          delta_e =(2.0_DP*en_offset + maxval(energies(ngwf_basis%num,:)) &
               - minval(energies(1,:)))/real(histonum-1, kind=DP)

          e_point = minval(energies(1,:)) - en_offset

          ! cks: Loop over DOS histogram points
          do row=1, histonum

             histo_val = 0.0_DP
             ! cks: accumulate value at current histogram point
             do orbital_index=1,ngwf_basis%num
                dist_sq =(e_point -energies(orbital_index, is))**2
!CW
               if(abs(alpha*dist_sq)<500.d0)then
!END CW

                histo_val = histo_val +gnorm*exp(-alpha *dist_sq)
!CW
               endif
!END CW
             end do

             ! cks: write to file current histo-point
             write(output_unit,'(f14.8,f16.10)',err=100) e_point, histo_val

             ! cks: update next histogram coordinate
             e_point = e_point +delta_e
          end do

          close(unit=output_unit,iostat=io_status)
          call utils_close_unit_check('eigenstates_dos_gp','output_unit', &
               io_status)

          ! cks: notify of end of output
          write(stdout,*)' done'

          write(stdout,'(a)')'================================&
               &================================================'

       enddo spin_loop



       deallocate(energies,stat=ierr)
       call utils_dealloc_check('eigenstates_dos_gp','energies',ierr)

    endif
    call comms_barrier

    return


100 write(stdout,*)'Problem writing to file in eigenstates_dos_gp. ONETEP stops'
    call comms_abort

  end subroutine eigenstates_dos_gp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_ldos_gp(eigen_en, eigenvecs_dens, s_eigenvecs_dens, &
       is, elements, ngwf_basis, hubbard_switch, ham_type)

    !===================================================================!
    ! This subroutine prepares and outputs in simple txt format a       !
    ! Local Density of States (LDOS) plot which has been generated      !
    ! with the Gaussian smearing method, using gamma point only         !
    ! energies. The groups over which to sum the LDOS are defined by    !
    ! pub_ldos_groups, which is set through the species_ldos_groups     !
    ! block in the input file. LDOS is calculated as a histogram of the !
    ! sum of all Gaussian-broadened energies multiplied by a coefficient!
    ! expressing the contribution of the NGWFs on each atom in the LDOS !
    ! group to that orbital. Energies in the output are in eV and so is !
    ! the half-width of the smearing Gaussians which is set by the      !
    ! input parameter ldos_smear.                                       !
    !-------------------------------------------------------------------!
    ! Written by Nicholas Hine on 15th February 2010, based on parts of !
    ! eigenstates_dos_gp which was written by Chris-Kriton Skylaris on  !
    ! 10/05/2006.                                                       !
    ! Modified to include Hubbard subspaces by David O'Regan, 10/2/2011.!
    ! Moved to eigenstates_mod by Nicholas Hine in November 2011.       !
    !===================================================================!

    use comms, only: comms_abort, comms_barrier, pub_my_node_id, &
         pub_on_root
    use constants, only: DP, UP, DN, stdout, HARTREE_IN_EVS, PI
    use dense, only: DEM, dense_get_col
    use function_basis, only: FUNC_BASIS
    use ion, only: ELEMENT
    use parallel_strategy, only: pub_orig_atom
    use rundat, only: pub_rootname, pub_ldos_smear, pub_ldos_ngroups, &
         pub_ldos_group_nsp, pub_ldos_groups, pub_cond_calculate
    use simulation_cell, only: pub_cell
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_open_unit_check, utils_close_unit_check

    implicit none

    type(FUNC_BASIS), intent(in) :: ngwf_basis   ! Function basis for NGWFs
    real(kind=DP), intent(in) :: eigen_en(ngwf_basis%num)
    type(DEM), intent(in) :: eigenvecs_dens
    type(DEM), intent(in) :: s_eigenvecs_dens
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    integer, intent(in) :: is                    ! Current spin
    logical, intent(in) :: hubbard_switch        ! ddor: DFT+U projected LDOS
    character(len=*), intent(in) :: ham_type

    ! Local Variables
    real(kind=DP) :: alpha     ! Gaussian exponent
    real(kind=DP) :: delta_e   ! energy increment
    real(kind=DP) :: e_point   ! energy point
    real(kind=DP) :: gnorm     ! Gaussian normalisation factor
    real(kind=DP) :: orb_coeff ! Coefficient of orbital for this NGWF
    real(kind=DP), allocatable, dimension(:) :: energies ! buffer for energies in eV
    real(kind=DP), allocatable, dimension(:,:) :: histo_val ! histogram values
    real(kind=DP), parameter :: en_offset =3.0_DP ! left and right offset for energy scale
    real(kind=DP), allocatable :: eig_col(:), s_eig_col(:)
    integer, parameter :: histonum=2001 ! number of histogram points
    integer :: row  ! row counting index
    integer :: output_unit ! fortran output unit number
    integer :: ierr ! memory allocation error flag
    integer :: igroup
    integer :: iat
    integer :: isp
    integer :: ingwf
    integer :: io_status ! file access error flag
    integer :: orbital_index ! orbital counting index
    logical :: iat_in_group
    character(len=256) :: output_file  ! output file name buffer
    character(len=4) :: iat_orig_id

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering eigenstates_ldos_gp'
#endif

    ! ndmh: Allocate temporary arrays
    allocate(energies(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('eigenstates_ldos_gp','energies',ierr)
    allocate(histo_val(histonum,0:pub_ldos_ngroups+1),stat=ierr)
    call utils_alloc_check('eigenstates_ldos_gp','histo_val',ierr)
    allocate(eig_col(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('eigenstates_ldos_gp','eig_col',ierr)
    allocate(s_eig_col(ngwf_basis%num),stat=ierr)
    call utils_alloc_check('eigenstates_ldos_gp','s_eig_col',ierr)

    ! cks: store eigen energies in buffer and convert to eV
    energies = eigen_en
    energies = energies*HARTREE_IN_EVS
    ! qoh: Initialise to avoid compiler warning
    output_unit = -1

    if (pub_on_root) then

       write(stdout,'(/a)')&
            '================== Local Density of States (LDOS) calculation &
            &=================='

       ! ndmh: find name for output files
       ! ddor: Include possibility of DFT+U projected LDOS
       ! ndmh: Include possibility of various COND task ham_types
       write(output_file,*) trim(pub_rootname)
       if (pub_cond_calculate) write(output_file,*) &
            trim(output_file)//'_'//trim(ham_type)
       if (hubbard_switch) write(output_file,*) &
            trim(output_file)//'_'//'HUB'
       if ((pub_cell%num_spins==2).and.(is == UP)) &
            write(output_file,*) trim(output_file)//'_'//'UP'
       if ((pub_cell%num_spins==2).and.(is == DN)) &
            write(output_file,*) trim(output_file)//'_'//'DN'
       write(output_file,*) trim(output_file)//'_LDOS.txt'
 
       output_file = adjustl(output_file)

       ! cks: print output warning
       write(stdout,'(3a)',advance ='no') 'Writing "', trim(output_file),'" ...'

       ! cks: get a unit number that is free
       output_unit = utils_unit()

       open(unit=output_unit, form="formatted" ,file=trim(output_file), &
            action="write",iostat=io_status)
       call utils_open_unit_check('eigenstates_ldos_gp','output_file', &
            io_status)

       ! ndmh: write first line
       write(output_unit,'(a)',advance='no',err =100)'#  Energy (eV) |'
       do igroup=1,pub_ldos_ngroups
          write(output_unit,'(a)',advance='no',err=100) '  LDOS ('
          do isp=1,pub_ldos_group_nsp(igroup)
             write(output_unit,'(a)',advance='no',err=100) &
                  pub_ldos_groups(isp,igroup)
             if (isp<pub_ldos_group_nsp(igroup)) &
                  write(output_unit,'(a)',advance='no',err=100) ', '
          end do
          write(output_unit,'(a)',advance='no',err=100) ') |'
       end do
       write(output_unit,'(a)',err=100) '  Total DOS    |'

    end if

    ! ndmh: initialisations
    histo_val(:,:) = 0.0_DP
    delta_e = (2.0_DP*en_offset + energies(ngwf_basis%num) - energies(1)) / &
         real(histonum-1, kind=DP)

    ! ndmh: smearing Gaussian exponent
    ! ndmh: pub_ldos_smear is the half-width in eVs
    alpha = log(2.0_DP) / ((pub_ldos_smear*HARTREE_IN_EVS)**2)
    gnorm = sqrt(alpha/PI)

    ! ndmh: first calculate the total DOS, as last column of histogram
    ! ddor: However this is not needed for DFT+U
    if (.not. hubbard_switch) then
       do orbital_index=1,ngwf_basis%num

          e_point = energies(1) - en_offset
          do row=1,histonum
!CW        
            if( abs(alpha * (e_point - energies(orbital_index))**2) <  500.d0)then
!END CW

             histo_val(row,pub_ldos_ngroups+1) = &
                  histo_val(row,pub_ldos_ngroups+1) &
                  + gnorm*exp(-alpha*(e_point - energies(orbital_index))**2)
!CW 
           endif
!END CW
             e_point = e_point + delta_e
          end do

       end do
    endif

    ! ndmh: accumulate contributions to LDOS for each orbital
    do orbital_index=1,ngwf_basis%num

       ! ndmh: Loop over LDOS histogram points to set up exp(-alpha*(e-ei)^2)
       ! ndmh: array for each point, where ei is the current eigenvalue
       e_point = energies(1) - en_offset
       do row=1,histonum
!CW        
            if( abs(alpha * (e_point - energies(orbital_index))**2) <  500.d0)then
!END CW
          histo_val(row,0) = &
               gnorm*exp(-alpha*(e_point - energies(orbital_index))**2)
!CW 
           endif
!END CW
          e_point = e_point + delta_e
       end do

       call dense_get_col(eig_col,eigenvecs_dens,orbital_index)
       call dense_get_col(s_eig_col,s_eigenvecs_dens,orbital_index)

       ! ndmh: loop over atoms to add up contributions to the local DOS
       do iat=1,pub_cell%nat
          iat_orig_id = elements(pub_orig_atom(iat))%species_id

          ! ndmh: go to next atom if this atom does not occur in any LDOS groups
          iat_in_group = .false.
          do igroup=1,pub_ldos_ngroups
             if (any(pub_ldos_groups(1:pub_ldos_group_nsp(igroup),igroup) &
                  ==iat_orig_id)) iat_in_group=.true.
          end do
          ! ddor: Full summation for DFT+U total LDOS
          if ((.not.iat_in_group) .and. (.not. hubbard_switch)) cycle

          ! ndmh: loop over NGWFs on this atom
          do ingwf=ngwf_basis%first_on_atom(iat), &
               ngwf_basis%first_on_atom(iat)+ngwf_basis%num_on_atom(iat)-1

             ! ndmh: evaluate orbital coefficient for NGWF a
             ! ndmh: given by M(dagger)_an * sum_b(S_ab M_nb)
             orb_coeff = eig_col(ingwf) * s_eig_col(ingwf)

             ! ndmh: add contribution of this eigenvalue to any of the LDOS
             ! ndmh: groups in which this atom is listed
             ! ddor: iat_in_group check for DFT+U
             if (iat_in_group) then
                do igroup=1,pub_ldos_ngroups

                   if (any(pub_ldos_groups(1:pub_ldos_group_nsp(igroup),igroup) &
                        ==iat_orig_id)) histo_val(:,igroup) = &
                        histo_val(:,igroup) + orb_coeff * histo_val(:,0)

                end do  ! igroup
             endif

             ! ddor: Total DFT+U LDOS
             if (hubbard_switch) then
                histo_val(:,pub_ldos_ngroups+1) = histo_val(:,pub_ldos_ngroups+1) + &
                     orb_coeff * histo_val(:,0)
             endif

          end do  ! ingwf

       end do  ! iat

    end do  ! orbital_index

    if (pub_on_root) then

       ! ndmh: loop over histogram points writing output file
       e_point = energies(1) - en_offset
       do row=1,histonum

          ! ndmh: write to file current histo-point for each group, and total
          write(output_unit,'(f14.8,f16.10)',advance='no',err=100) e_point
          do igroup=1,pub_ldos_ngroups
             write(output_unit,'(f16.10)',advance='no',err=100) &
                  histo_val(row,igroup)
          end do
          write(output_unit,'(f16.10)') histo_val(row,pub_ldos_ngroups+1)

          e_point = e_point + delta_e

       end do

       ! ndmh: close output file
       close(unit=output_unit,iostat=io_status)
       call utils_close_unit_check('eigenstates_ldos_gp','output_unit', &
            io_status)

       ! ndmh: finish writing message to stdout
       write(stdout,*)' done'
       write(stdout,'(a)')'================================&
            &================================================'

    endif

    ! ndmh: deallocate temporary arrays
    deallocate(s_eig_col,stat=ierr)
    call utils_dealloc_check('eigenstates_ldos_gp','s_eig_col',ierr)
    deallocate(eig_col,stat=ierr)
    call utils_dealloc_check('eigenstates_ldos_gp','eig_col',ierr)
    deallocate(histo_val,stat=ierr)
    call utils_dealloc_check('eigenstates_ldos_gp','histo_val',ierr)
    deallocate(energies,stat=ierr)
    call utils_dealloc_check('eigenstates_ldos_gp','energies',ierr)

    call comms_barrier

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving eigenstates_ldos_gp'
#endif

    return


100 write(stdout,*)'Problem writing to file in eigenstates_ldos_gp. ONETEP stops'
    call comms_abort

  end subroutine eigenstates_ldos_gp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_write_bands(ngwf_basis, eigen_en, n_occ, ham_type)

    !============================================================!
    ! This subroutine writes to a file the orbitals in CASTEP    !
    ! .bands format.                                             !
    !------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 18/02/2007.            !
    ! Bug in Fermi level choice for fully-spin polarised systems !
    ! fixed by Nicholas Hine on 08/06/2010. Also cleaned up the  !
    ! code somewhat.                                             !
    ! Moved to eigenstates_mod by Nicholas Hine in November 2011.!
    !============================================================!

    use comms, only: pub_on_root, comms_abort, comms_barrier
    use constants, only: DP, max_spins, stdout, stderr
    use function_basis, only: FUNC_BASIS
    use rundat, only: pub_rootname
    use simulation_cell, only: pub_cell
    use utils, only: utils_abort, utils_unit, utils_open_unit_check, &
         utils_close_unit_check

    implicit none

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind =DP), intent(in) :: eigen_en(:,:) ! band eigen-energies
    integer, intent(in)        :: n_occ(:)      ! number of occupied bands
    ! lr408: Optional argument to set different filenames for different diagonalisations in
    ! lr408: a conduction calculation
    character(*), optional, intent(in) :: ham_type

    ! Local Variables
    integer :: output_unit  ! unit name of output file
    integer :: ierr         ! unit open/close error flag
    integer :: is           ! spin counter
    integer :: neigen       ! eigenvalues counter
    character(len=256) :: filename  ! output file name buffer
    real(kind =DP)     :: efermi(max_spins) ! "Fermi level" for each spin

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &eigenstates_write_bands'
#endif

    if (pub_on_root) then

       ! cks: Initialise "Fermi levels" to HOMO energies
       ! ndmh: modified to check for invalid n_occ
       do is=1,pub_cell%num_spins
          if ((n_occ(is) > 0).and.(n_occ(is) <= ngwf_basis%num)) then
             efermi(is) = eigen_en(n_occ(is),is)
          else if (n_occ(is) == 0) then
             efermi(is) = eigen_en(1,is)
          else
             call utils_abort('Error in eigenstates_write_bands: n_occ(is) < 0 &
                  &or n_occ(is) > N_NGWFs')
          end if
       end do

       ! cks: name of output file
       if (present(ham_type)) then
          if (ham_type=='valence') then
             write(filename,*)trim(pub_rootname)//'.val_bands'
          else if (ham_type=='proj') then
             write(filename,*)trim(pub_rootname)//'.proj_bands'
          else if (ham_type=='cond') then
             write(filename,*)trim(pub_rootname)//'.cond_bands'
          else if (ham_type=='joint') then
             write(filename,*)trim(pub_rootname)//'.joint_bands'
          end if
       else
          write(filename,*)trim(pub_rootname)//'.bands'
       end if

       filename = adjustl(filename)

       write(stdout,*)''
       write(stdout,'(3a)', advance='no')&
            'Writing bands to file "',trim(adjustl(filename)),'" ...'

       ! cks: find free unit to open
       output_unit =utils_unit()

       ! cks: open output unit
       open(unit=output_unit, file=trim(filename), form='formatted', &
            action='write', iostat=ierr, status='replace')
       call utils_open_unit_check('eigenstates_write_bands','output_unit',ierr)

       ! Write the number of k-points, spin components and eigenvalues to file
       write(output_unit,'(a,i4)') 'Number of k-points        ', 1
       write(output_unit,'(a,i4)') 'Number of spin components ', &
            pub_cell%num_spins

       if (pub_cell%num_spins == 1) then
          write(output_unit,'(a,g10.4)') 'Number of electrons       ', &
               2*n_occ(1)
       else
          write(output_unit,'(a,2g10.4)') 'Number of electrons       ', &
               n_occ(1), n_occ(2)
       end if

       if (pub_cell%num_spins == 1) then
          write(output_unit,'(a,i6)') 'Number of eigenvalues     ', &
               ngwf_basis%num
       else
          write(output_unit,'(a,2i6)') 'Number of eigenvalues     ', &
               ngwf_basis%num, ngwf_basis%num
       end if

       if( pub_cell%num_spins == 1 ) then
          write(output_unit,'(a,f12.6)') 'Fermi energy (in atomic units) ', &
               efermi(1)
       else
          write(output_unit,'(a,2f12.6)') 'Fermi energy (in atomic units) ', &
               efermi(1:pub_cell%num_spins)
       end if

       write(output_unit,'(a)') 'Unit cell vectors'
       write(output_unit,'(3f12.6)') pub_cell%a1%x, pub_cell%a1%y, pub_cell%a1%z
       write(output_unit,'(3f12.6)') pub_cell%a2%x, pub_cell%a2%y, pub_cell%a2%z
       write(output_unit,'(3f12.6)') pub_cell%a3%x, pub_cell%a3%y, pub_cell%a3%z

       ! cks: write the actual band energies for each spin component
       write(output_unit,'(a,i4,4f12.8)') 'K-point ', &
            1, 0.0_DP, 0.0_DP, 0.0_DP, 1.0_DP
       do is=1,pub_cell%num_spins
          write(output_unit,'(a,i1)') 'Spin component ',is
          do neigen=1,ngwf_basis%num
             write(output_unit,'(f14.8)') eigen_en(neigen,is)
          end do
       end do

       ! cks: close the output unit (file)
       close(unit=output_unit,iostat=ierr)
       call utils_close_unit_check('eigenstates_write_bands','output_unit',ierr)
       write(stdout,'(a)') ' done'

    endif
    call comms_barrier

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &eigenstates_write_bands'
#endif

    return

  end subroutine eigenstates_write_bands


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine eigenstates_print_ens_occs(eigen_en, eigen_occ, n_occ, neig)

    !=================================================================!
    ! This subroutine prints the eigenvalues of the Hamiltonian and   !
    ! the density kernel and prints maximum and minimum values as     !
    ! well as a range of them around the band gap.                    !
    !-----------------------------------------------------------------!
    ! NOTE:                                                           !
    !   The occupancies printed by this subroutine are the projection !
    !   of the orbitals |phi_c> M^c_n onto the density operator       !
    !   given by |phi_a>K^ab<phi_b|.                                  !
    !-----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 21/05/2006 by modifying     !
    ! internal_find_eigen_info which was written by Chris-Kriton      !
    ! Skylaris and Peter Haynes in 2004-2005.                         !
    ! Modified format strings, and enabled printing of wider range    !
    ! of unoccupied states if required, by Nicholas Hine on 19/10/09. !
    ! Moved to eigenstates_mod by Nicholas Hine in November 2011.     !
    !=================================================================!

    use comms, only: pub_on_root, comms_barrier
    use constants, only: DP, stdout
    use rundat, only: pub_num_eigenvalues, print_qc
    use simulation_cell, only: pub_cell

    implicit none

    real(kind =DP), intent(in) :: eigen_en(:,:)  ! orbital energies from Hamiltonian
    real(kind =DP), intent(in) :: eigen_occ(:,:) ! orbital occupancies from density kernel
    integer, intent(in)        :: n_occ(:)       ! number of occupied orbitals of each spin
    integer, intent(in)        :: neig           ! number of NGWFs

    ! Local variables
    integer :: is         ! Spin counter
    integer :: ieig       ! Loop counter
    integer :: num_print  ! Number of eigenvalues to actually print


#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering &
         &eigenstates_print_ens_occs'
#endif

    ! Diagonalisation is not parallelised, so only work on root
    if (pub_on_root) then

       ! pdh: loop over spins
       do is=1,pub_cell%num_spins

          ! Print out the required results
          if (pub_cell%num_spins == 1) then
             write(stdout,'(/a/)') '=================== &
                  &Orbital energy and occupancy information ==================='
          else
             write(stdout,'(/a,i1,a/)') '============== &
                  &Orbital energy and occupancy information for spin ',is, &
                  ' ============='
          end if
          write(stdout,'(a,i8)') '      Total number of orbitals:', neig
          write(stdout,'(a,i8)') '   Number of occupied orbitals:', n_occ(is)
          write(stdout,'(a,f18.9)') '                 Occupancy sum:', &
               sum(eigen_occ(:,is))

          ! pdh: check for no states occupied as well as all
          if (n_occ(is)<neig .and. n_occ(is) > 0) then
             write(stdout,'(a,f18.9,a)') '                 HOMO-LUMO gap:', &
                  eigen_en(n_occ(is)+1,is) - eigen_en(n_occ(is),is),' Eh'
             write(stdout,'(a,f18.9,a)') '                 Mid-gap level:', &
                  0.5_DP*(eigen_en(n_occ(is)+1,is) + eigen_en(n_occ(is),is)),' Eh'
          else if (n_occ(is)==0) then
             write(stdout,'(a,a)') '                 HOMO-LUMO gap:', &
                  '  Unavailable - no orbitals occupied'
             write(stdout,'(a,a)') '                 Mid-gap level:', &
                  '  Unavailable - no orbitals occupied'
          else
             write(stdout,'(a,a)') '                 HOMO-LUMO gap:', &
                  '  Unavailable - all orbitals occupied'
             write(stdout,'(a,a)') '                 Mid-gap level:', &
                  '  Unavailable - all orbitals occupied'
          end if
          write(stdout,'(/a)') '                        &
               &Orbital | Energy (Eh) | Occupancy'

          ! pdh: don't print lowest state if no states are occupied
          if (n_occ(is) > 0) write(stdout,'(a,i6,f16.9,f12.7)') &
               '                       ',1,eigen_en(1,is), eigen_occ(1,is)

          ! pdh: avoid repeating lowest state
          num_print = max(n_occ(is)-pub_num_eigenvalues,1)+1

          if (num_print > 2) write(stdout,'(a)') &
               '                        .......   ...........   .........'

          do ieig=num_print,n_occ(is)
             write(stdout,'(a,i8,f16.9,f12.7)')'                     ', &
                  ieig, eigen_en(ieig,is), eigen_occ(ieig,is)
          end do

          ! pdh: don't print gap if no states are occupied
          if (n_occ(is) > 0) write(stdout,'(a)') '                        &
               &.......   --- gap ---   .........'

          ! Don't print more virtual orbitals than are available
          num_print = min(n_occ(is)+pub_num_eigenvalues,neig-1)

          do ieig=n_occ(is)+1,num_print
             write(stdout,'(a,i6,f16.9,f12.7)') '                       ', &
                  ieig, eigen_en(ieig,is), eigen_occ(ieig,is)
          end do
          if (num_print < neig-1) write(stdout,'(a)') '                        &
               &.......   ...........   .........'

          ! pdh: don't print highest state if all states are occupied
          if (n_occ(is) < neig) write(stdout,'(a,i6,f16.9,f12.7)') &
               '                       ',neig, eigen_en(neig,is), &
               eigen_occ(neig,is)

          if (print_qc) then
             write(stdout,*)
             write(stdout,'(a,f20.12)') '<QC>     [energy_eigenvalue_1]:', &
                  eigen_en(1,is)
             write(stdout,'(a,f20.12)') '<QC>     [energy_eigenvalue_n]:', &
                  eigen_en(neig,is)
             write(stdout,'(a,f14.06)') '<QC>             [occupancy_1]:', &
                  eigen_occ(1,is)
             write(stdout,'(a,f14.06)') '<QC>             [occupancy_n]:', &
                  eigen_occ(neig,is)
          end if

          write(stdout,'(/a)') '=========================================&
               &======================================='

       end do

    end if

    call comms_barrier

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving &
         &eigenstates_print_ens_occs'
#endif

  end subroutine eigenstates_print_ens_occs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module eigenstates

