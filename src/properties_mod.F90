! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!                    Property analysis module                    !
!                                                                !
! This module performs property analysis using canonical         !
! molecular orbitals obtained from converged ONETEP calculations.!
!----------------------------------------------------------------!
! Written by Chris-Kriton Skylaris on 10/05/2006                 !
! Population analysis added by Peter Haynes, July 2006           !
! NGWF analysis added by Mark Robinson, January 2008             !
! General re-writing and updating by Nicholas Hine, 2009-2011.   !
!================================================================!

module properties

  use constants, only: DP

  implicit none

  private

  public :: properties_calculate
  public :: properties_polarisation

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine properties_calculate(denskern, ham, &
       rep, ngwf_basis, proj_basis, hub_proj_basis, hub, &
       elements, localpseudo_fine, core_density_fine, lhxc_fine, &
       properties_only, ham_type, num_opt_states)

    !======================================================================!
    ! This subroutine generates and diagonalises a Hamiltonian and         !
    ! performs various kinds of analysis with the resulting eigenvectors   !
    ! and eigenenergies.                                                   !
    !----------------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 10/05/2006.                      !
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
    !======================================================================!

!CW
    use comms, only: pub_on_root, comms_abort, comms_exit, comms_barrier
    use rundat, only: pub_dmft_points,pub_dmft_sc,pub_dmft_nkpoints
    use restart, only: restart_kernel_write,restart_kernel_read
    use kernel, only : kernel_init_core_ham
!END CW

    use bandstructure, only: bandstructure_calculate
    use cell_grid, only: pub_fine_grid
    use comms, only: pub_on_root
    use constants, only: max_spins, stdout, UP, DN, HARTREE_IN_EVS, &
         paw_en_size
    use eigenstates, only: eigenstates_calculate
    use function_basis, only: FUNC_BASIS, &
         function_basis_sum_ppd_funcs ! lpl: NBO NGWF analysis
    use hamiltonian, only: hamiltonian_dens_dep_matrices, &
         hamiltonian_build_matrix
    use hubbard_build, only: HUBBARD_MODEL
    use ion, only: ELEMENT
    use ngwf_representation, only: NGWF_REP, NGWF_HAM
    use npa, only: npa_main
    use rundat, only: pub_homo_dens_plot, pub_lumo_dens_plot, pub_homo_plot, &
         pub_lumo_plot, pub_dos_smear, pub_ldos_smear, pub_ldos_ngroups, &
         pub_num_eigenvalues, pub_write_density_plot, pub_grd_format, &
         pub_cube_format, pub_dx_format, pub_popn_calculate, &
         pub_ngwf_analysis, pub_nnho, pub_do_bandstructure, &
         pub_polarisation_calculate, pub_spread_calculate, &
         pub_aug, pub_cond_calculate, pub_etrans_calculate, &
         pub_plot_nbo, pub_write_nbo, pub_nbo_pnao_analysis
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_scale
    use timer, only: timer_clock
    use transport, only: etrans_calculate
    use utils, only: utils_alloc_check, utils_dealloc_check
!CW
    use kernel, only : kernel_purify
    use hubbard_build, only: hubbard_dmft_interface
    use sparse,only : sparse_trace
    use rundat,only : pub_dmft_kernel,pub_dmft_purify_sc
!END CW

    implicit none

    ! Arguments
    type(NGWF_REP), intent(in) :: rep
    type(NGWF_HAM), intent(inout) :: ham   ! NGWF Hamiltonian type - changed to argument by lr408
    type(SPAM3), intent(inout) :: denskern(pub_cell%num_spins) ! denskern matrix
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(FUNC_BASIS), intent(in) :: hub_proj_basis
    type(HUBBARD_MODEL), intent(inout) :: hub
    type(ELEMENT), intent(in) :: elements(:) ! elements of all atoms in  input file order
    real(kind=DP), intent(in) :: localpseudo_fine(:,:,:) ! local pseudo on this proc slab
    real(kind=DP), intent(in) :: core_density_fine(:,:,:) ! core density on this proc slab
    real(kind=DP), intent(inout) :: lhxc_fine(pub_fine_grid%ld1,&
         pub_fine_grid%ld2, pub_fine_grid%max_slabs12, pub_cell%num_spins)
    logical, intent(in) :: properties_only ! flag determining whether Hamiltonian needs calculating
    ! lr408: Optional string defining which type of Hamiltonian this is for
    character(len=*), intent(in) :: ham_type
    ! lr408: Optional integer defining the total number of states to include for optical spectra
    integer, optional, intent(in) :: num_opt_states

    ! Local Variables
    real(kind=DP), dimension(:), allocatable :: orth_on_grid   ! on-site orthogonal NGWFs
    real(kind=DP) :: total_energy      ! total energy
    real(kind=DP) :: lhxc_energy       ! pseudo+hartree+xc energy
    real(kind=DP) :: hubbard_energy    ! Hubbard energy
    real(kind=DP) :: paw_sphere_energies(paw_en_size)
    integer :: ierr                    ! memory allocation error flag
    type(SPAM3) :: pnao_tr(1)          ! NBO NGWF analysis
!CW
    logical :: check,dmft_sc_back
    real(kind=DP) :: eigen_en(ngwf_basis%num,pub_cell%num_spins)
    integer :: is
    TYPE(SPAM3) :: purkern(pub_cell%num_spins)
!END CW

!CW
   if(pub_dmft_points>0)then
     INQUIRE(file='store_ham1',EXIST=check)
     if(check) goto 77
   endif
!END CW

    ! cks: start timer
    call timer_clock("properties_calculate",1)

    ! cks: === BUILD HAMILTONIAN & WRITE DENSITY/POTENTIAL ===================

    ! ndmh: If this is a stand-alone properties calculation, we need to
    ! ndmh: calculate the Hamiltonian. Otherwise, it will already have been
    ! ndmh: calculated and this can be skipped.
    ! ndmh: 25/01/2011: removed this option to skip recalculation, since there
    ! ndmh: were too many circumstances where the hamiltonian was inconsistent
    ! ndmh: with the current NGWFs or density kernel and thus wrong.
    if ((properties_only.or.(.true.)).and.(.not. pub_cond_calculate)) then
       ! cks: scale to take into account alpha plus beta electrons

!CW
!added lines for DFT+DMFT
   dmft_sc_back=pub_dmft_sc
   pub_dmft_sc=.false.
   if(pub_dmft_purify_sc.and.pub_dmft_kernel/=0)then
        do is=1,pub_cell%num_spins
         call sparse_create(purkern(is),denskern(is))
        enddo
        call kernel_purify(purkern,denskern,rep%overlap,rep%inv_overlap,rep%n_occ)
        if (pub_cell%num_spins == 1) call sparse_scale(purkern(1),2.0_DP)
        call hamiltonian_dens_dep_matrices(ham, lhxc_fine, total_energy, &
            lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
            ngwf_basis, hub_proj_basis, hub, purkern, 0.0_DP, elements, &
            pub_fine_grid, localpseudo_fine, core_density_fine, &
            ham_update=.true., lhxc_fixed=.false.,h_xc_write=.true.)
        if (pub_cell%num_spins == 1) call sparse_scale(purkern(1),0.5_DP)
        if(pub_dmft_points>0)then
           open(unit=5050,file='store_tot_energy')
           write(5050,*) total_energy,lhxc_energy,hubbard_energy
           close(5050)
        endif
          call hamiltonian_build_matrix(ham, rep)
       do is=1,pub_cell%num_spins
          call sparse_destroy(purkern(is))
       enddo
   else
!END CW
       if (pub_cell%num_spins == 1) call sparse_scale(denskern(1),2.0_DP)
       ! ndmh: calculate density dependent energies and matrices
       call hamiltonian_dens_dep_matrices(ham, lhxc_fine, total_energy, &
            lhxc_energy, hubbard_energy, paw_sphere_energies, rep, &
            ngwf_basis, hub_proj_basis, hub, denskern, 0.0_DP, elements, &
            pub_fine_grid, localpseudo_fine, core_density_fine, &
            ham_update=.true., lhxc_fixed=.false.,h_xc_write=.true.)
!CW
        if(pub_dmft_points>0)then
           open(unit=5050,file='store_tot_energy')
           write(5050,*) total_energy,lhxc_energy,hubbard_energy
           close(5050)
        endif
!END CW
       ! ndmh: build Hamiltonian matrix from its component matrices
       call hamiltonian_build_matrix(ham, rep)

       ! cks: undo scaling
       if (pub_cell%num_spins == 1) call sparse_scale(denskern(1),0.5_DP)
    end if
!CW
   endif
   pub_dmft_sc=dmft_sc_back
!END CW

    if (pub_write_density_plot .and. &
         (pub_cube_format .or. pub_grd_format .or. pub_dx_format)) then
       ! cks: write potential and density plot files
       call internal_plot_dens_and_pot
    endif

    ! cks: END BUILD HAMILTONIAN & WRITE DENSITY/POTENTIAL ===================

    ! ndmh: population analysis, NGWF analysis and polarisation analysis
    ! ndmh: all only work for valence Hamiltonian
    if (ham_type=='valence') then

       ! pdh: POPULATION ANALYSIS
       if (pub_popn_calculate) call &
            properties_popn_analysis(rep%overlap,denskern,rep%inv_overlap, &
            elements,ngwf_basis,rep%n_occ)
       ! pdh: END POPULATION ANALYSIS

       ! mr: NGWF ANALYSIS
       if ((pub_nnho) .and. (pub_ngwf_analysis .or. pub_spread_calculate)) then
          ! orthogonalise NNHO before characterising
          allocate(orth_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
          call utils_alloc_check('properties_calculate','orth_on_grid',ierr)
          call ngwfs_orthogonalise(orth_on_grid,rep%ngwfs_on_grid,ngwf_basis, &
               rep%overlap)
          ! mr: print s/p/d character of NGWFs
          if (pub_ngwf_analysis) then
             call properties_ngwfs_char(orth_on_grid,ngwf_basis,elements,3)
          endif
          ! ddor: Print NGWF spreads
          if (pub_spread_calculate) then
             call properties_spread(orth_on_grid,ngwf_basis,elements, &
                  rep%inv_overlap,rep%overlap,.false.)
          endif
          deallocate(orth_on_grid,stat=ierr)
          call utils_dealloc_check('properties_calculate','orth_on_grid',ierr)
       else
          ! mr: print s/p/d character of NGWFs
          if (pub_ngwf_analysis) then
             call properties_ngwfs_char(rep%ngwfs_on_grid,ngwf_basis,elements,3)
          endif
          ! ddor: Print NGWF spreads
          if (pub_spread_calculate) then
             call properties_spread(rep%ngwfs_on_grid,ngwf_basis,elements, &
                  rep%inv_overlap,rep%overlap,.true.)
          endif
       end if
       ! mr: END NGWF ANALYSIS


       ! mr: POLARISATION
       ! ndmh: flag to activate polarisation calculations
       if (pub_polarisation_calculate) then
          call properties_polarisation(rep,ngwf_basis,proj_basis,elements, &
               denskern)
       end if
       ! mr: END POLARISATION

    end if

    ! pdh: Only diagonalise Hamiltonian if required
!CW
    if ((pub_num_eigenvalues > 0) .or. (pub_dos_smear > 0.0_DP) .or. &
         ((pub_ldos_smear > 0).and.(pub_ldos_ngroups > 0)) .or. &
         (pub_homo_dens_plot >= 0) .or. (pub_lumo_dens_plot >= 0) .or. &
         (pub_homo_plot >= 0) .or. (pub_lumo_plot >= 0) .or. &
         pub_cond_calculate .or. pub_dmft_points>0 ) then
         call eigenstates_calculate(denskern, ham, &
                       rep, ngwf_basis, hub_proj_basis, elements, ham_type, &
                       num_opt_states,proj_basis=proj_basis,eigen_en_copy=eigen_en)
!         call eigenstates_calculate(denskern, ham, &
!                       rep, ngwf_basis, proj_basis, hub_proj_basis, elements, ham_type, &
!                       num_opt_states,eigen_en)
!END CW
    end if


    ! pdh: BANDSTRUCTURE
!CW    
    ! if (pub_do_bandstructure.and.(ham_type/='cond')) &
      if ( (pub_do_bandstructure.and.(ham_type/='cond')) .or.  abs(pub_dmft_nkpoints)>1 ) &
    !     call bandstructure_calculate(ham,rep,ngwf_basis,proj_basis, &
    !     elements,ham_type)
          call bandstructure_calculate(ham,rep,ngwf_basis,proj_basis, &
          elements,ham_type,hub=hub,hub_proj_basis=hub_proj_basis)
    ! pdh: BANDSTRUCTURE
!END CW


!CW
!##################################################################!
       if (pub_dmft_points > 0) then
          77 continue
          write(*,*) 'hubbard dmft interface'
          call hubbard_dmft_interface(eigen_en,rep%n_occ, &
               ham%ham,rep%overlap,rep%inv_overlap,ngwf_basis,hub_proj_basis, &
               hub,rep,elements,denskern,proj_basis=proj_basis)
          write(*,*) 'hubbard dmft interface done'
          call properties_popn_analysis(rep%overlap,denskern,rep%inv_overlap,elements,ngwf_basis,rep%n_occ,nopurify=.true.)
          call comms_barrier
          call comms_abort
       endif
!##################################################################!
!END CW

    ! lpl: NBO
    if(pub_write_nbo) then
       if(pub_nbo_pnao_analysis) then
            pnao_tr(1)%structure = 'D'
            call sparse_create(pnao_tr(1))
            call npa_main(denskern,rep,ham,elements,ngwf_basis,pnao_tr(1))

            allocate(orth_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts),stat=ierr)
            call utils_alloc_check('properties_calculate','orth_on_grid',ierr)

            orth_on_grid = 0.0_DP
            call function_basis_sum_ppd_funcs(orth_on_grid,ngwf_basis, &
                 pnao_tr,1,1,pnao_tr(1),rep%ngwfs_on_grid,ngwf_basis)

            ! lpl: NGWF analysis
            if(pub_on_root) write(stdout,'(a)') &
                 '============ PNAO s/p/d/f Character ============'
            call properties_ngwfs_char(orth_on_grid,ngwf_basis,elements,3)
            if(pub_on_root) write(stdout,'(a)') &
                 '================================================'
            deallocate(orth_on_grid,stat=ierr)
            call utils_dealloc_check('properties_calculate','orth_on_grid',ierr)

            call sparse_destroy(pnao_tr(1))
       else
            call npa_main(denskern,rep,ham,elements,ngwf_basis)
       end if
    end if
    if(pub_plot_nbo) call internal_plot_nbo
    ! lpl: NBO

    ! smmd: TRANSMISSION FUNCTION 
    if (pub_etrans_calculate .and. (ham_type/='cond')) then
      call etrans_calculate(ham,rep,ngwf_basis,elements)
    end if
    ! smmd: TRANSMISSION FUNCTION 

    ! cks: stop timer
    call timer_clock("properties_calculate",2)

    contains

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine internal_plot_dens_and_pot
      ! Trivially modified by Jacek Dziedzic on 14/05/2010 to delegate   !
      ! the unit conversion to visual_scalarfield.                       !
      ! ndmh: added PAW compensation density, and changed variable names !
      ! to match other parts of the code 29/09/2010.                     !

        use augmentation, only: augmentation_density_on_grid, &
             aug_projector_denskern
        use cell_grid, only: pub_fine_grid
        use constants, only: DP, ANGSTROM
        use density, only: density_on_grid
        use paw, only: paw_sphere_density_on_grid
        use rundat, only: pub_aug_den_dim
        use sparse, only: SPAM3, sparse_create, sparse_destroy
        use visual, only: visual_scalarfield

        implicit none

        ! Local Variables
        real(kind=DP), allocatable :: density_fine(:,:,:,:)! density
        real(kind=DP), allocatable :: density_plot(:,:,:,:)! density
        real(kind=DP), allocatable :: nhat_den_grad(:,:,:,:,:)! compensation den
        type(SPAM3), allocatable :: rhoij(:)
        integer :: ierr          ! memory allocation error flag
        integer :: is

        if (pub_on_root) then
           write(stdout,'(/a)') '======================== &
                &Writing density and potential &
                &========================='
        end if

        ! cks:########## DENSITY OUTPUT ##########

        ! cks: Allocate charge density slabs
        allocate(density_fine(pub_fine_grid%ld1,pub_fine_grid%ld2, &
             pub_fine_grid%max_slabs12,pub_cell%num_spins),stat=ierr)
        call utils_alloc_check('internal_plot_dens_and_pot &
             &(properties_calculate)','density_fine',ierr)

        ! ndmh: Scale denskern for spin degeneracy
        if (pub_cell%num_spins == 1) then
           call sparse_scale(denskern(1),pub_cell%spin_fac)
        end if

        ! cks: Calculate data-parallelised charge density
        call density_on_grid(density_fine,pub_fine_grid, &
             denskern,rep%ngwf_overlap,rep%ngwfs_on_grid,ngwf_basis)

        ! ndmh: Calculate compensation density on fine grid
        if (pub_aug) then ! PAW version

           allocate(nhat_den_grad(pub_fine_grid%ld1,pub_fine_grid%ld2, &
                pub_fine_grid%max_slabs12,pub_cell%num_spins, &
                0:pub_aug_den_dim),stat=ierr)
           call utils_alloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','nhat_den_grad',ierr)
           allocate(density_plot(pub_fine_grid%ld1,pub_fine_grid%ld2, &
                pub_fine_grid%max_slabs12,pub_cell%num_spins),stat=ierr)
           call utils_alloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','density_plot',ierr)

           ! CREATE PSEUDO DENSITY (INCLUDING NHAT) AND OUTPUT

           nhat_den_grad = 0.0_DP
           call augmentation_density_on_grid(nhat_den_grad,pub_fine_grid, &
                denskern,rep%sp_overlap)

           density_plot = density_fine + nhat_den_grad(:,:,:,:,0)

           ! ndmh: spin polarisation: calculate total density and spin density
           if (pub_cell%num_spins == 2) then
              density_plot(:,:,:,1) = density_plot(:,:,:,UP) + &
                   density_plot(:,:,:,DN)
              density_plot(:,:,:,2) = density_plot(:,:,:,UP) - &
                   2.0_DP * density_plot(:,:,:,DN)
           end if

           ! cks: output density in plot format file
           ! vm: output density in Angstrom rather than in Bohr
           ! jd: ... but leave the unit conversion to visual_scalarfield
           call visual_scalarfield( &
                density_plot(:,:,:,1), pub_fine_grid, &
                'Electronic pseudo density (in e/ang^3) for:', '_ps_density', &
                elements, ANGSTROM**3)
           if (pub_cell%num_spins == 2) call visual_scalarfield( &
                density_plot(:,:,:,2), pub_fine_grid, &
                'Electronic pseudo spin density (in e/ang^3) for:', &
                '_ps_spindensity', elements, ANGSTROM**3)


           ! CREATE AE DENSITY AND OUTPUT

           allocate(rhoij(pub_cell%num_spins),stat=ierr)
           call utils_alloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','rhoij',ierr)
           do is=1,pub_cell%num_spins
              rhoij(is)%structure = 'E'
              call sparse_create(rhoij(is))
           end do
           call aug_projector_denskern(rhoij,denskern,rep%sp_overlap)

           nhat_den_grad(:,:,:,:,0) = 0.0_DP
           call paw_sphere_density_on_grid(nhat_den_grad(:,:,:,:,0), &
                pub_fine_grid,rhoij,1.0_DP,-1.0_DP)
           density_plot = density_fine + nhat_den_grad(:,:,:,:,0)

           ! ndmh: spin polarisation: calculate total density and spin density
           if (pub_cell%num_spins == 2) then
              density_plot(:,:,:,1) = density_plot(:,:,:,UP) + &
                   density_plot(:,:,:,DN)
              density_plot(:,:,:,2) = density_plot(:,:,:,UP) - &
                   2.0_DP * density_plot(:,:,:,DN)
           end if

           ! cks: output density in plot format file
           ! vm: output density in Angstrom rather than in Bohr
           ! jd: ... but leave the unit conversion to visual_scalarfield
           call visual_scalarfield( &
                density_plot(:,:,:,1), pub_fine_grid, &
                'All-electron density (in e/ang^3) for:', '_ae_density', &
                elements, ANGSTROM**3)
           if (pub_cell%num_spins == 2) call visual_scalarfield( &
                density_plot(:,:,:,2), pub_fine_grid, &
                'All-electron spin density (in e/ang^3) for:', &
                '_ae_spindensity', elements, ANGSTROM**3)

           do is=pub_cell%num_spins,1,-1
              call sparse_destroy(rhoij(is))
           end do
           deallocate(rhoij,stat=ierr)
           call utils_dealloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','rhoij',ierr)
           deallocate(density_plot,stat=ierr)
           call utils_dealloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','density_plot',ierr)
           deallocate(nhat_den_grad,stat=ierr)
           call utils_dealloc_check('internal_plot_dens_and_pot &
                &(properties_calculate)','nhat_den_grad',ierr)

        else ! non-PAW version

           ! ndmh: spin polarisation: calculate total density and spin density
           if (pub_cell%num_spins == 2) then
              density_fine(:,:,:,1) = density_fine(:,:,:,UP) + &
                   density_fine(:,:,:,DN)
              density_fine(:,:,:,2) = density_fine(:,:,:,UP) - &
                   2.0_DP * density_fine(:,:,:,DN)
           end if

           ! cks: output density in plot format file
           ! vm: output density in Angstrom rather than in Bohr
           ! jd: ... but leave the unit conversion to visual_scalarfield
           call visual_scalarfield( &
                density_fine(:,:,:,1), pub_fine_grid, &
                'Electronic density (in e/ang^3) for:', '_density', elements, &
                ANGSTROM**3)
           if (pub_cell%num_spins == 2) call visual_scalarfield( &
                density_fine(:,:,:,2), pub_fine_grid, &
                'Electronic spin density (in e/ang^3) for:', '_spindensity', &
                elements, ANGSTROM**3)

        end if

        ! ndmh: restore normalisation of denskern
        if (pub_cell%num_spins == 1) then
           call sparse_scale(denskern(1),1.0_DP/pub_cell%spin_fac)
        end if

        ! cks: Deallocate charge density slabs
        deallocate(density_fine,stat=ierr)
        call utils_dealloc_check('internal_plot_dens_and_pot &
             &(properties_calculate)','density_fine',ierr)

        ! cks:####### END DENSITY OUTPUT ##########



        ! cks:########## POTENTIAL OUTPUT ##########

        !cks: output potential in plot format file
        ! vm: output potential in eV rather than in Hartree
        !     also use "-" to do ESP: potential seen by a positive probe
        ! jd: ... but leave the unit conversion to visual_scalarfield

        ! pdh: only one potential even if spin-polarised
!        if (pub_cell%num_spins == 1) then
        call visual_scalarfield( &
             lhxc_fine(:,:,:,1), pub_fine_grid, &
             'Local potential (Ion+Hartree+XC) for:', &
             '_potential', elements, -HARTREE_IN_EVS)
!        else
!           call visual_scalarfield( &
!                lhxc_fine(:,:,:,1), pub_fine_grid, &
!                'Local up-spin potential (Ion+Hartree+XC) for:', &
!                '_potential_UP', elements, -HARTREE_IN_EVS)
!           call visual_scalarfield( &
!                lhxc_fine(:,:,:,2), pub_fine_grid, &
!                'Local down-spin potential (Ion+Hartree+XC) for:', &
!                '_potential_DN', elements, -HARTREE_IN_EVS)
!        end if

        ! cks:###### END POTENTIAL OUTPUT ##########

        if (pub_on_root) write(stdout,'(a)')'================================&
             &================================================'


      end subroutine internal_plot_dens_and_pot

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine internal_plot_nbo

        ! lpl: Reads GENNBO output files to plot relevant orbitals
        !      Uses properties_plot_mo
        use augmentation, only: augmentation_overlap
        use comms, only: comms_abort, comms_barrier, comms_bcast, &
             pub_root_node_id
        use dense, only: DEM, dense_create, dense_destroy, dense_get_col
        use eigenstates, only: eigenstates_plot_mo
        use rundat, only: pub_nbo_list_plotnbo, pub_nbo_plotorbtype, &
             pub_rootname
        use npa, only: npa_read_nbodata
        implicit none
 
        ! Local Variables
        type(SPAM3) :: mo_kern(max_spins), aug_overlap
        type(DEM) :: nbongwf_tr
        real(kind=DP), allocatable :: col_buffer(:)
        integer :: it, num_ao, icol
        character(len=64) :: output_tag, output_file, title, nbo_trfile
        logical :: pre_ao

#ifdef DEBUG
        if(pub_on_root) write(stdout,'(a)') 'Entering internal_plot_nbo'
#endif

        ! lpl: Initialize matrices and arrays
        call dense_create(nbongwf_tr,ngwf_basis%num,ngwf_basis%num)
  
        allocate(col_buffer(ngwf_basis%num),stat=ierr)
        call utils_alloc_check('properties_plot_nbo','col_buffer',ierr)

        ! ndmh: NBO 'kernel' and aug overlap part for dealing with augmentation
        if (pub_aug) then
           call sparse_create(mo_kern(1),denskern(1))
           call sparse_create(aug_overlap,rep%overlap)
           call augmentation_overlap(aug_overlap,rep%sp_overlap)
        end if

        ! lpl: Determine type of orbital and appropriate NBO input filename
        pub_nbo_plotorbtype = adjustl(pub_nbo_plotorbtype)
        if(pub_on_root) then
           select case(trim(pub_nbo_plotorbtype))
              case ('PNAO')
                 pre_ao = .true.
                 nbo_trfile = 'NONE'
                 output_tag = '_pnao'
              case ('NAO')
                 pre_ao = .false.
                 nbo_trfile = 'NONE'
                 output_tag = '_nao'
              case ('PNHO')
                 pre_ao = .true.
                 nbo_trfile = '.35'
                 output_tag = '_pnho'
              case ('NHO')
                 pre_ao = .false.
                 nbo_trfile = '.35'
                 output_tag = '_nho'
              case ('PNBO')
                 pre_ao = .true.
                 nbo_trfile = '.37'
                 output_tag = '_pnbo'
              case ('NBO')
                 pre_ao = .false.
                 nbo_trfile = '.37'
                 output_tag = '_nbo'
              case ('NLMO')
                 pre_ao = .false.
                 nbo_trfile = '.39'
                 output_tag = '_nlmo'
              case default
                 write(stdout,'(a)') 'ERROR: Invalid nbo_plot_orbtype'
                 ierr = 1
           end select

           nbo_trfile = adjustl(nbo_trfile)
           if(nbo_trfile /= 'NONE') then
              write(nbo_trfile,'(a)') trim(pub_rootname)//'_nao'//trim(nbo_trfile)
              nbo_trfile = adjustl(nbo_trfile)
           end if

           output_tag = adjustl(output_tag)

        end if
        if(ierr /= 0) call comms_abort

        call comms_bcast(pub_root_node_id,pre_ao)
        call comms_bcast(pub_root_node_id,nbo_trfile)
        call comms_bcast(pub_root_node_id,output_tag)

        ! lpl: Obtain full NGWF to NBO's orbital transformation
        call npa_read_nbodata(nbongwf_tr,ngwf_basis,trim(nbo_trfile), &
             pre_ao,num_ao)

#ifdef DEBUG
        if(pub_on_root) then
           write(stdout,'(a)') 'DEBUG info: NBOs read for plotting'
           write(stdout,'(a)') '----------------------------------'
           do it=1,pub_nbo_list_plotnbo(0)
              write(stdout,'(1x,I5,1x,I5)') it,pub_nbo_list_plotnbo(it)
           end do
           write(stdout,'(a)') '----------------------------------'
        end if
#endif
        ! lpl: Plot orbitals
        do it=1,pub_nbo_list_plotnbo(0)
           write(title,'(I5)') pub_nbo_list_plotnbo(it)

           output_file = trim(output_tag)//trim(adjustl(title))
           title = trim(pub_nbo_plotorbtype)//' # '//trim(adjustl(title))

           if(pub_nbo_list_plotnbo(it) <= num_ao) then
              ! lpl: Get NBO in run-time NGWF basis
              call dense_get_col(col_buffer,nbongwf_tr, &
                   pub_nbo_list_plotnbo(it))
              call comms_barrier

              call eigenstates_plot_mo(col_buffer(:),rep%ngwfs_on_grid, &
                   ngwf_basis,elements,aug_overlap,mo_kern(1), &
                   trim(adjustl(title)),trim(adjustl(output_file)))
           else
              write(stdout,'(a,I5)') 'ERROR: # orbs > num_ao. Skipping...',it
           end if

           call comms_barrier

        end do

        if (pub_aug) then
           call sparse_destroy(aug_overlap)
           call sparse_destroy(mo_kern(1))
        end if

        deallocate(col_buffer,stat=ierr)
        call utils_dealloc_check('properties_plot_nbo','col_buffer',ierr)

        call dense_destroy(nbongwf_tr)

#ifdef DEBUG
        if(pub_on_root) write(stdout,'(a)') 'Exiting internal_plot_nbo'
#endif

      end subroutine internal_plot_nbo


  end subroutine properties_calculate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CW
! added nopurify argument, overlap->overlap_in
  subroutine properties_popn_analysis(overlap_in,denskern,inv_overlap_in,elements, &
       ngwf_basis,n_occ,nopurify)
!ENDCW
    !======================================================================!
    ! This subroutine performs Mulliken population analysis using the      !
    ! NGWFs as the set of local orbitals.                                  !
    !----------------------------------------------------------------------!
    ! Written by Peter Haynes on 31/07/2006.                               !
    !======================================================================!
!CW added comms_barrier
    use comms, only: comms_bcast, pub_on_root, pub_my_node_id, &
         pub_root_node_id, comms_reduce,comms_barrier
!END CW
    use constants, only: DP, max_spins, stdout, UP, DN, TWO_PI
    use function_basis, only: FUNC_BASIS
    use geometry, only: point, operator(*), operator(.dot.), operator(+),&
         magnitude
    use ion, only: element
    use kernel, only: kernel_purify
    use parallel_strategy, only: pub_num_atoms_on_node, pub_num_overlaps, &
         pub_first_atom_on_node, pub_orig_atom, pub_distr_atom, &
         pub_overlap_list, pub_max_overlaps, parallel_strategy_list_overlaps
    use rundat, only : pub_popn_bond_cutoff, print_qc
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_get_block, &
         sparse_product
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_heapsort
!CW
    use sparse, only: sparse_copy
!END CW
    implicit none

    ! Arguments:
!CW
    type(SPAM3), intent(in) :: overlap_in                      ! Overlap matrix
    type(SPAM3), intent(inout) :: denskern(pub_cell%num_spins) ! Density kernel
    type(SPAM3), intent(in) :: inv_overlap_in                  ! Inverse overlap
!END CW
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)        ! Elements
    type(FUNC_BASIS), intent(in) :: ngwf_basis                 ! NGWF Basis type
    integer, intent(in) :: n_occ(max_spins)                    ! Occupied bands

    ! Local variables:
    type(SPAM3) :: ks                              ! K.S product
    type(SPAM3), allocatable :: purkern(:)         ! Purified kernel K
    type(POINT) :: bond                            ! Bond
    integer :: ierr                                ! Error flag
    integer :: is                                  ! Spin counter
    integer :: loc_iat                             ! Local atom counter
    integer :: iat,jat                             ! Atom counters
    integer :: iat_orig                            ! Atom iat in original order
    integer :: ingwf,jngwf                         ! NGWF counters
    integer :: ingwfs                              ! Number of NGWFs on iat
    integer :: iel                                 ! Element counter
    integer :: max_ngwfs_on_atom                   ! Max NGWFs on atom
    integer :: iovlap                              ! Overlap counter
    integer :: num_bonds                           ! The number of bonds
    integer :: nbond                               ! Number of bonds for output
    integer :: nb                                  ! Running index of bonds
    integer, allocatable :: bond_index(:)          ! Stores the index for
                                                   ! ordering output arrays
    real(kind=DP) :: block_tr                      ! Block trace
    real(kind=DP) :: popn                          ! Population
    real(kind=DP) :: charge,spin                   ! Atomic charge and spin
    real(kind=DP) :: bondvec(3)                    ! Bond vector
    real(kind=DP) :: frac(3,2)                     ! Fractional coordinates
    real(kind=DP), allocatable :: q_atom(:,:)      ! Atomic charges
    real(kind=DP), allocatable :: q_bond(:,:,:)    ! Bond charges
    real(kind=DP), allocatable :: over_block(:,:)  ! Block of overlap matrix
    real(kind=DP), allocatable :: kern_block(:,:)  ! Block of kernel
    real(kind=DP), allocatable :: bond_len(:)      ! Stores the length of each bond (for output)
    real(kind=DP), allocatable :: bond_pop(:)      ! Stores the population of each bond (for output)
    real(kind=DP), allocatable :: bond_spin(:)     ! Stores the spin of each bond (for output)
    character(len=26), allocatable :: bond_name(:) ! Stores the label of each bond (for output)
!CW
    logical,optional    :: nopurify
    type(SPAM3)         :: overlap,inv_overlap 
    logical,allocatable :: orb_in_plane(:,:)
    integer             :: in1,in2,ufii,nplanes
    logical             :: isorb
    integer             :: jp

    call sparse_create(     overlap ,     overlap_in, iscmplx=.false.)
    call sparse_create( inv_overlap , inv_overlap_in, iscmplx=.false.)
    call sparse_copy(    overlap,     overlap_in)
    call sparse_copy(inv_overlap, inv_overlap_in) 
    call check_sparse_in_file(    overlap,    'store_overlap')
    call check_sparse_in_file(inv_overlap,'store_inv_overlap')
  
   nplanes=0; ufii=2010

   inquire(file='group_of_atoms_label_sorted',exist=isorb)
   if(isorb)then
      open(unit=ufii,file='group_of_atoms_label_sorted',form='unformatted')
        read(ufii) nplanes
        read(ufii) in1,in2
        if(allocated(orb_in_plane)) deallocate(orb_in_plane); allocate(orb_in_plane(in1,in2))
        read(ufii) orb_in_plane
      close(ufii)
   endif

!ENDCW
    ! Allocate space for purified kernel (according to overlap sparsity)
    allocate(purkern(pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('properties_popn_analysis','purkern',ierr)
    do is=1,pub_cell%num_spins
       call sparse_create(purkern(is),overlap)
    end do
    ! Calculate purified kernel
!CW
   if(.not.present(nopurify)) then
!ENDCW
     call kernel_purify(purkern,denskern,overlap,inv_overlap,n_occ )
!CW
   else
     do is=1,pub_cell%num_spins
       call sparse_copy(purkern(is),denskern(is)) 
     enddo
   endif
!ENDCW

    !--------------------
    ! Atomic populations
    !--------------------

    ! Space for blocks
    max_ngwfs_on_atom = maxval(ngwf_basis%num_on_atom)
    allocate(over_block(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
    call utils_alloc_check('properties_popn_analysis','over_block',ierr)
    allocate(kern_block(max_ngwfs_on_atom,max_ngwfs_on_atom),stat=ierr)
    call utils_alloc_check('properties_popn_analysis','kern_block',ierr)

    ! Generate block diagonal matrix for K.S product
    ks%structure = 'D'
    call sparse_create(ks)

    ! Allocate space for atomic populations
    allocate(q_atom(pub_cell%nat,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('properties_popn_analysis','q_atom',ierr)
    q_atom = 0.0_DP

    ! Loop over spins
    do is=1,pub_cell%num_spins

       ! Calculate product of density kernel and overlap
       call sparse_product(ks,purkern(is),overlap)

       iat = pub_first_atom_on_node(pub_my_node_id)
       do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)

          ! Calculate atomic populations on this node
          call sparse_get_block(over_block,ks,iat,iat)

          block_tr = 0.0_DP
          do ingwf=1,ngwf_basis%num_on_atom(iat)
             block_tr = block_tr + over_block(ingwf,ingwf)
          end do
          q_atom(pub_orig_atom(iat),is) = block_tr
          iat = iat + 1
       end do

    end do

    ! Sum up over all procs
    call comms_reduce('SUM',q_atom)

    ! Write out results
    if (pub_on_root) then

       ! ndmh: print results for first three atoms as a QC test
       if (print_qc) then
          if (pub_cell%num_spins == 1) then
             do iat=1,min(3,pub_cell%nat)
                popn = 2.0_DP * q_atom(iat,1)
                charge = elements(iat)%ion_charge - popn
                write(stdout,'(a,i1,a,f10.6)') &
                     '<QC>   [atom_',iat,'_population]: ',popn
                write(stdout,'(a,i1,a,f10.6)') &
                     '<QC>       [atom_',iat,'_charge]: ',charge
             end do
          else
             do iat=1,min(3,pub_cell%nat)
                popn = q_atom(iat,UP) + q_atom(iat,DN)
                charge = elements(iat)%ion_charge - popn
                spin = 0.5_DP * (q_atom(iat,UP) - q_atom(iat,DN))
                write(stdout,'(a,i1,a,f10.6)') &
                     '<QC>   [atom_',iat,'_population]: ',popn
                write(stdout,'(a,i1,a,f10.6)') &
                     '<QC>       [atom_',iat,'_charge]: ',charge
                write(stdout,'(a,i1,a,f10.6)') &
                     '<QC>         [atom_',iat,'_spin]: ',spin
             end do
          end if
       end if

       write(stdout,'(/a)') '    Atomic Populations'
       write(stdout,'(a)')  '    ------------------'
       if (pub_cell%num_spins == 1) then
          write(stdout,'(a)') 'Species  Ion    Total   Charge (e)'
          write(stdout,'(a)') '=================================='
       else
          write(stdout,'(a)') 'Species  Ion    Total   Charge (e)    &
               &Spin (hbar)'
          write(stdout,'(a)') '======================================&
               &==========='
       end if

       if (pub_cell%num_spins == 1) then
          do iat=1,pub_cell%nat
             popn = 2.0_DP * q_atom(iat,1)
             charge = elements(iat)%ion_charge - popn
             write(stdout,'(2x,a2,1x,i5,2x,f10.6,2x,f9.6)') &
                  elements(iat)%symbol, iat, popn, charge
          end do
          write(stdout,'(a)') '=================================='
       else
          do iat=1,pub_cell%nat
             popn = q_atom(iat,UP) + q_atom(iat,DN)
             charge = elements(iat)%ion_charge - popn
             spin = 0.5_DP * (q_atom(iat,UP) - q_atom(iat,DN))
             write(stdout,'(2x,a2,1x,i5,2x,f10.6,2x,f9.6,4x,f6.2)') &
                  elements(iat)%symbol, iat, popn, charge, spin
          end do
          write(stdout,'(a)') '========================================&
               &========='
       end if
!CW
   if(nplanes>0.and.isorb)Then
      do jp=1,nplanes
       write(stdout,'(a)') '=================================='
       write(stdout,'(a)') '=================================='
       write(stdout,'(a,i4)') 'plane : ', jp
       if (pub_cell%num_spins == 1) then
          do iat=1,pub_cell%nat
           if(.not.orb_in_plane(jp,iat)) cycle
             popn = 2.0_DP * q_atom(iat,1)
             charge = elements(iat)%ion_charge - popn
             write(stdout,'(2x,a2,1x,i5,2x,f10.6,2x,f9.6)') &
                  elements(iat)%symbol, iat, popn, charge
          end do
          write(stdout,'(a)') '=================================='
       else
          do iat=1,pub_cell%nat
             if(.not.orb_in_plane(jp,iat)) cycle
             popn = q_atom(iat,UP) + q_atom(iat,DN)
             charge = elements(iat)%ion_charge - popn
             spin = 0.5_DP * (q_atom(iat,UP) - q_atom(iat,DN))
             write(stdout,'(2x,a2,1x,i5,2x,f10.6,2x,f9.6,4x,f6.2)') &
                  elements(iat)%symbol, iat, popn, charge, spin
          end do
          write(stdout,'(a)') '========================================&
               &========='
       end if
      enddo
    endif
!END CW
    end if

    ! Deallocate workspace
    deallocate(q_atom,stat=ierr)
    call utils_dealloc_check('properties_popn_analysis','q_atom',ierr)
    call sparse_destroy(ks)

    !------------------
    ! Bond populations
    !------------------

    ! Obtain list of overlaps
    call parallel_strategy_list_overlaps(elements,'Fixed','Fixed', &
         0.5_DP*pub_popn_bond_cutoff)

    ! Allocate space for output arrays (find estimate for number of bonds first)
    num_bonds = 0
    if (pub_on_root) then
       do iat = 1, pub_cell%nat
          num_bonds = num_bonds + pub_num_overlaps(iat)
       end do
    end if
    call comms_bcast(pub_root_node_id,num_bonds)

    allocate(bond_pop(num_bonds),stat=ierr)
    call utils_alloc_check('properties_popn_analysis','bond_pop',ierr)
    if (pub_cell%num_spins == 2) then
        allocate(bond_spin(num_bonds),stat=ierr)
        call utils_alloc_check('properties_popn_analysis','bond_spin',ierr)
    end if
    allocate(bond_len(num_bonds),stat=ierr)
    call utils_alloc_check('properties_popn_analysis','bond_len',ierr)
    allocate(bond_name(num_bonds),stat=ierr)
    call utils_alloc_check('properties_popn_analysis','bond_name',ierr)
    allocate(bond_index(num_bonds),stat=ierr)
    call utils_alloc_check('properties_popn_analysis','bond_index',ierr)

    ! Allocate space for atomic populations
    allocate(q_bond(pub_max_overlaps,pub_cell%nat,pub_cell%num_spins),stat=ierr)
    call utils_alloc_check('properties_popn_analysis','q_bond',ierr)
    q_bond = 0.0_DP

    ! Loop over atoms on this processor
    iat = pub_first_atom_on_node(pub_my_node_id)
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)
       ingwfs = ngwf_basis%num_on_atom(iat)
       iat_orig = pub_orig_atom(iat)
       do iovlap=1,pub_num_overlaps(iat_orig)
          jat = pub_distr_atom(pub_overlap_list(iovlap,iat_orig))
          call sparse_get_block(over_block,overlap,jat,iat)
          do is=1,pub_cell%num_spins
             call sparse_get_block(kern_block,purkern(is),jat,iat)
             block_tr = 0.0_DP
             do ingwf=1,ingwfs
                do jngwf=1,ngwf_basis%num_on_atom(jat)
                   block_tr = block_tr + over_block(jngwf,ingwf) * &
                        kern_block(jngwf,ingwf)
                end do
             end do
             q_bond(iovlap,iat_orig,is) = block_tr
          end do
       end do
       iat = iat + 1
    end do

    ! Sum up over all procs
    call comms_reduce('SUM',q_bond)

    ! Factor of two in definition
    q_bond = q_bond * 2.0_DP

    nbond = 0

    ! Write out results
    if (pub_on_root) then
       do iat=1,pub_cell%nat
          frac(1,1) = (elements(iat)%centre .dot. pub_cell%b1) / TWO_PI
          frac(2,1) = (elements(iat)%centre .dot. pub_cell%b2) / TWO_PI
          frac(3,1) = (elements(iat)%centre .dot. pub_cell%b3) / TWO_PI
          do iovlap=1,pub_num_overlaps(iat)
             jat = pub_overlap_list(iovlap,iat)
             if (jat <= iat) cycle
             frac(1,2) = (elements(jat)%centre .dot. pub_cell%b1) / TWO_PI
             frac(2,2) = (elements(jat)%centre .dot. pub_cell%b2) / TWO_PI
             frac(3,2) = (elements(jat)%centre .dot. pub_cell%b3) / TWO_PI
             ! Bond length
             bondvec = abs(frac(:,1)-frac(:,2))
             bondvec = modulo(bondvec,1.0_DP)
             do iel=1,3
                if (bondvec(iel) > 0.5_DP) bondvec(iel) = bondvec(iel) - 1.0_DP
             end do
             bond = bondvec(1) * pub_cell%a1 + bondvec(2) * pub_cell%a2 + &
                  bondvec(3) * pub_cell%a3
             nbond = nbond + 1
             bond_len(nbond) = magnitude(bond)

             write (bond_name(nbond),'(1x,a2,i5,a4,a2,i5)') &
                   elements(iat)%symbol, iat, ' -- ', elements(jat)%symbol, jat

             if (pub_cell%num_spins == 1) then
                bond_pop(nbond)  = 2.0_DP * q_bond(iovlap,iat,1)
             else if(pub_cell%num_spins == 2) then
                bond_pop(nbond)  = q_bond(iovlap,iat,UP) + q_bond(iovlap,iat,DN)
                bond_spin(nbond) = 0.5_DP * (q_bond(iovlap,iat,UP) - q_bond(iovlap,iat,DN))
             end if
          end do
       end do

       call utils_heapsort(nbond,bond_len,bond_index)

       if (pub_cell%num_spins == 1) then
          write(stdout,'(/a)') '        Bond         Population      &
               &Length (bohr)'
          write(stdout,'(a)') '======================================&
               &============'
          do nb = 1,nbond
              write(stdout,730)  bond_name(bond_index(nb)),bond_pop(bond_index(nb)),&
                                 bond_len(nb)
          end do
          write(stdout,'(a/)') '========================================&
               &=========='

       else
          write(stdout,'(/a)') '        Bond         Population      &
               &Spin       Length (bohr)'
          write(stdout,'(a)') '======================================&
               &======================='
          do nb = 1,nbond
              write(stdout,740)  bond_name(bond_index(nb)),bond_pop(bond_index(nb)),&
                                 bond_spin(bond_index(nb)),bond_len(nb)
          end do
          write(stdout,'(a/)') '========================================&
               &====================='
       end if

730    format(a23,t23,f9.2,t35,f10.5)
740    format(a23,t23,f9.2,t34,f8.2,t45,f10.5)

    end if

    ! Deallocate output arrays
    deallocate(bond_pop,stat=ierr)
    call utils_dealloc_check('properties_popn_analysis','bond_pop',ierr)
    if (pub_cell%num_spins == 2) then
       deallocate(bond_spin,stat=ierr)
       call utils_dealloc_check('properties_popn_analysis','bond_spin',ierr)
    end if
    deallocate(bond_len,stat=ierr)
    call utils_dealloc_check('properties_popn_analysis','bond_len',ierr)
    deallocate(bond_name,stat=ierr)
    call utils_dealloc_check('properties_popn_analysis','bond_name',ierr)
    deallocate(bond_index,stat=ierr)
    call utils_dealloc_check('properties_popn_analysis','bond_index',ierr)

    ! Destroy workspace
    deallocate(q_bond,stat=ierr)
    call utils_dealloc_check('properties_popn_analysis','q_bond',ierr)
    deallocate(kern_block,stat=ierr)
    call utils_dealloc_check('properties_popn_analysis','kern_block',ierr)
    deallocate(over_block,stat=ierr)
    call utils_dealloc_check('properties_popn_analysis','over_block',ierr)
    do is=pub_cell%num_spins,1,-1
       call sparse_destroy(purkern(is))
    end do
    deallocate(purkern,stat=ierr)
    call utils_dealloc_check('properties_popn_analysis','purkern',ierr)

!CW
    call sparse_destroy(    overlap)
    call sparse_destroy(inv_overlap)

    if(allocated(orb_in_plane)) deallocate(orb_in_plane)
     
    contains
      !---------------------------------------!
     subroutine check_sparse_in_file(mat,filename)
     implicit none
     type(spam3)    :: mat
     character*(*)  :: filename
     logical        :: check
       call comms_barrier
       INQUIRE(file=filename,EXIST=check)
       call comms_barrier
       if(.not.check)then
        if(pub_my_node_id==0) write(*,*) 'popn analysis : no store file present'
       else
        if(pub_my_node_id==0) write(*,*) 'popn analysis : store file present=',trim(adjustl(filename))
        call sparse_read_(mat,filename=filename)
       endif
       call comms_barrier
     end subroutine
      !---------------------------------------!
     subroutine sparse_read_(mat,filename)
     use sparse,   only: sparse_convert,sparse_num_rows,sparse_num_cols
     implicit none
     type(spam3)                                :: mat
     real(kind=DP), allocatable, dimension(:,:) :: mat_square
     integer                                    :: n1,n2
     character*(*)                              :: filename
        call comms_barrier
        n1=NINT(sparse_num_rows(mat))
        n2=NINT(sparse_num_cols(mat))
        if(n1==0.and.n2==0)then
          write(*,*) 'SPARSE READ, return, 0-shape matrix'
          return
        endif
        open(unit=20001,file=filename,form='unformatted')
        allocate(mat_square(n1,n2),stat=ierr)
        read(20001) mat_square
        call sparse_convert(mat,mat_square)
        deallocate(mat_square,stat=ierr)
        close(20001)
     end subroutine
      !---------------------------------------!
!ENDCW
  end subroutine properties_popn_analysis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine ngwfs_orthogonalise(orth_on_grid,ngwfs_on_grid,ngwf_basis,overlap)

    !==========================================================================!
    ! Symmetrically orthogonalise NGWFs on each atom using the Lowdin          !
    ! transformation.                                                          !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! ngwfs_on_grid (input)    : NGWFs in ppd representation on this node      !
    ! ngwf_basis (input)       : The function basis for the NGWFs              !
    ! orth_on_grid (output)    : Orthogonalised NGWFs in ppd representation    !
    !--------------------------------------------------------------------------!
    ! Written by Mark Robinson, January 2008                                   !
    ! Cleaned up and made to use new code by Nicholas Hine, 08/06/2010         !
    ! Moved to ngwfs_mod by Nicholas Hine.                                     !
    !==========================================================================!

    use comms, only: pub_my_node_id, pub_on_root
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS, function_basis_sum_ppd_funcs
    use parallel_strategy, only: pub_num_atoms_on_node, pub_orig_atom, &
         pub_first_atom_on_node
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_get_block, &
         sparse_put_block
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check
    use wrappers, only: wrappers_dsygv_lt

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(out)    :: orth_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts,1)
    real(kind=DP), intent(in)     :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(SPAM3), intent(in) :: overlap

    ! Local variables
    type(SPAM3) :: dmat(1)
    integer :: num_ngwfs
    integer :: loc_iat, iat
    integer :: i
    integer :: ierr

    ! matrices and arrays
    real(kind=DP), allocatable, dimension(:,:) :: overlap_block,transform,d,temp
    real(kind=DP), allocatable, dimension(:)   :: evals

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering ngwfs_orthogonalise'
#endif

    ! Start timer
    call timer_clock('ngwfs_orthogonalise',1)

    ! Allocate temporary arrays
    num_ngwfs = ngwf_basis%max_on_atom
    allocate(overlap_block(num_ngwfs,num_ngwfs),stat=ierr)
    call utils_alloc_check('ngwfs_orthogonalise','overlap_block',ierr)
    allocate(transform(num_ngwfs,num_ngwfs),stat=ierr)
    call utils_alloc_check('ngwfs_orthogonalise','transform',ierr)
    allocate(d(num_ngwfs,num_ngwfs),stat=ierr)
    call utils_alloc_check('ngwfs_orthogonalise','d',ierr)
    allocate(temp(num_ngwfs,num_ngwfs),stat=ierr)
    call utils_alloc_check('ngwfs_orthogonalise','temp',ierr)
    allocate(evals(num_ngwfs),stat=ierr)
    call utils_alloc_check('ngwfs_orthogonalise','evals',ierr)

    dmat(1)%structure = 'D'
    call sparse_create(dmat(1))

    ! Loop over atoms on this node
    do loc_iat=1,pub_num_atoms_on_node(pub_my_node_id)

       iat = loc_iat + pub_first_atom_on_node(pub_my_node_id) - 1

       ! Find number of NGWFs on atom
       num_ngwfs = ngwf_basis%num_on_atom(iat)

       ! Only need to orthogonalise if more than one NGWF
       if (num_ngwfs > 1) then

          ! Get lower triangle of overlap matrix of NGWFs on atom
          overlap_block = 0
          call sparse_get_block(overlap_block,overlap,iat,iat)

          ! Solve the eigenvalue problem
          transform = overlap_block
          temp = 0
          do i=1,num_ngwfs
             temp(i,i) = 1.0_DP
          enddo
          call wrappers_dsygv_lt(transform,evals,temp,num_ngwfs)

          ! Construct d(-1/2) from evals
          d = 0
          do i=1,num_ngwfs
             d(i,i) = 1 / sqrt(evals(i))
          enddo

          ! Transform d^(-1/2) to S^(-1/2)
          ! X=d*transform'
          call dgemm('N','T',num_ngwfs,num_ngwfs,num_ngwfs,1.0_DP, &
                     d,num_ngwfs,transform,num_ngwfs,&
                     0.0_DP,temp,num_ngwfs)
          ! S^-1/2=transform*X
          call dgemm('N','N',num_ngwfs,num_ngwfs,num_ngwfs,1.0_DP, &
                     transform,num_ngwfs,temp,num_ngwfs,&
                     0.0_DP,d,num_ngwfs)

          ! Put resulting matrix block into dmat
          call sparse_put_block(d,dmat(1),iat,iat)

       else
          ! Copy block straight from overlap matrix
          call sparse_get_block(d,overlap,iat,iat)
          call sparse_put_block(d,dmat(1),iat,iat)
       endif

    enddo

    orth_on_grid = 0
    call function_basis_sum_ppd_funcs(orth_on_grid,ngwf_basis, &
         dmat,1,1,dmat(1),ngwfs_on_grid,ngwf_basis)

    ! Deallocate temporary arrays
    call sparse_destroy(dmat(1))
    deallocate(evals,stat=ierr)
    call utils_dealloc_check('ngwfs_orthogonalise','evals',ierr)
    deallocate(temp,stat=ierr)
    call utils_dealloc_check('ngwfs_orthogonalise','temp',ierr)
    deallocate(d,stat=ierr)
    call utils_dealloc_check('ngwfs_orthogonalise','d',ierr)
    deallocate(transform,stat=ierr)
    call utils_dealloc_check('ngwfs_orthogonalise','transform',ierr)
    deallocate(overlap_block,stat=ierr)
    call utils_dealloc_check('ngwfs_orthogonalise','overlap_block',ierr)

    ! Stop timer
    call timer_clock('ngwfs_orthogonalise',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving ngwfs_orthogonalise'
#endif

  end subroutine ngwfs_orthogonalise


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine properties_ngwfs_char(ngwfs_on_grid,ngwf_basis,elements,max_l, &
       nmax,output)

    !==========================================================================!
    ! Prints out an analysis of the s/p/d/f character of the NGWFs             !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! ngwfs_on_grid (input)    : NGWFs in ppd representation on this node      !
    ! ngwf_basis (input)       : The function basis for the NGWFs              !
    ! elements (input)         : Element array for atoms                       !
    ! lmax (input)             : Largest angular momentum SW to include        !
    ! nmax (input)             : Size of spherical wave basis set              !
    ! output (input)           : Should we write out SW coefficients to file   !
    !--------------------------------------------------------------------------!
    ! Written by Mark Robinson, November 2007                                  !
    !==========================================================================!

    use basis, only: basis_extract_function_from_box, &
         basis_location_func_wrt_cell
    use comms, only: pub_on_root, pub_root_node_id, pub_total_num_nodes, &
         pub_my_node_id, comms_recv, comms_send, comms_barrier, comms_reduce
    use constants, only: DP, stdout, two_pi
    use function_basis, only: FUNC_BASIS
    use geometry, only: point, operator(.DOT.), operator(*), operator(+)
    use ion, only: element
    use parallel_strategy, only: pub_orig_atom, pub_distr_atom
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit, &
         utils_open_unit_check, utils_close_unit_check
    use simulation_cell, only: pub_cell, pub_fftbox
    use spherical_wave, only: sw_bessel_zeros,sw_init, sw_recp_generate
    use timer, only: timer_clock

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    integer, intent(in)       :: max_l
    integer, optional, intent(in) :: nmax
    logical, optional, intent(in) :: output

    ! default to not printing out SW coefficients
    logical :: sw2file                           ! if true write SW coeffs to file
    integer :: max_n

    ! counters
    integer       :: atom,orig_atom,loc_ngwf,ngwf
    integer       :: l,m,n,node        ! loop counters
    integer       :: ngwf_on_atom      ! more counters

    ! variables
    type(POINT)   :: fbcentre,centre,frac,disp
    real(kind=DP) :: radius,tot,norm
    integer       :: ierr,offset,npoints,lastp
    integer       :: sx,sy,sz

    ! spherical bessel function zeros array
    real(kind=DP), allocatable :: qnl(:)

    ! SW FFTbox and tightbox
    real(kind=DP), allocatable :: sw_fftbox(:,:,:)
    real(kind=DP), allocatable :: sw_on_grid(:)
    complex(kind=DP), allocatable :: sw_work(:,:,:)

    ! summation arrays
    real(kind=DP), allocatable :: my_spd_sum(:,:)
    real(kind=DP), allocatable :: snd_spd_sum(:,:)
    real(kind=DP), allocatable :: spd_sum(:,:)

    ! file data (when printing out SW coefficients)
    integer           :: fileunit
    character(len=16) :: filename

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Entering properties_ngwfs_char'
#endif

    ! start timer
    call timer_clock('properties_ngwfs_char',1)

    ! check optional arguments
    if (present (nmax)) then
       max_n = nmax
    else
       ! default is number of psincs across largest NGWF
       max_n = ceiling(maxval(ngwf_basis%spheres(:)%radius) / &
                       min(pub_cell%d1,pub_cell%d2,pub_cell%d3))
       call comms_reduce('MAX',max_n)
    endif

    sw2file = .false.
    if (present (output)) sw2file = output

    ! allocate qnl array
    allocate(qnl(max_n),stat=ierr)
    call utils_alloc_check('properties_ngwfs_char','qnl',ierr)

    ! init SW module if not already - with max kr
    radius = maxval(ngwf_basis%spheres(:)%radius)**2 * 2 &
         * maxval(pub_fftbox%recip_grid(5,:,:,:))
    call comms_reduce('MAX',radius)
    call sw_init(max_l,max_n,radius)

    ! calculate shift required to centre SW in FFTbox on a grid point
    fbcentre = real(pub_fftbox%total_pt1-1,kind=DP)/2/pub_fftbox%total_pt1 * pub_fftbox%a1 + &
               real(pub_fftbox%total_pt2-1,kind=DP)/2/pub_fftbox%total_pt2 * pub_fftbox%a2 + &
               real(pub_fftbox%total_pt3-1,kind=DP)/2/pub_fftbox%total_pt3 * pub_fftbox%a3

    ! obtain free unit and number of ngwfs prior to this node
    ! (if writing SW coefficients to file)
    fileunit = utils_unit ()

    ! allocate workspace
    allocate(my_spd_sum(0:max_l+2,ngwf_basis%node_num),stat=ierr)
    call utils_alloc_check('properties_ngwfs_char','my_spd_sum',ierr)
    allocate(sw_fftbox(pub_fftbox%total_ld1,pub_fftbox%total_ld2,&
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('properties_ngwfs_char','sw_fftbox',ierr)
    allocate(sw_work(pub_fftbox%total_ld1,pub_fftbox%total_ld2,&
         pub_fftbox%total_pt3),stat=ierr)
    call utils_alloc_check('properties_ngwfs_char','sw_work',ierr)
    ! allocate sw_on_grid array
    allocate(sw_on_grid(ngwf_basis%max_n_ppds_sphere*pub_cell%n_pts),stat=ierr)
    call utils_alloc_check('properties_ngwfs_char','sw_on_grid',ierr)

    ! zero summation array
    my_spd_sum = 0

    ! loop over ngwfs on this node
    loc_ngwf = 0
    do ngwf=ngwf_basis%first_on_node(pub_my_node_id), &
         ngwf_basis%first_on_node(pub_my_node_id+1)-1
       loc_ngwf = loc_ngwf + 1

       orig_atom = pub_orig_atom(ngwf_basis%atom_of_func(ngwf))
       ngwf_on_atom = ngwf - &
            ngwf_basis%first_on_atom(ngwf_basis%atom_of_func(ngwf)) + 1

       ! open file to write SW coefficients to if requested
       if (sw2file) then
          write(filename,'(a5,i5.5,a6)') 'atom_',orig_atom, &
               'ngwf_',ngwf_on_atom, '.swdat'
          open(unit=fileunit,file=trim(filename),action='write',iostat=ierr)
          call utils_open_unit_check('properties_ngwfs_char','fileunit',ierr)
       endif

       ! NGWF center and location in large array
       centre = ngwf_basis%spheres(loc_ngwf)%centre
       radius = ngwf_basis%spheres(loc_ngwf)%radius
       offset = ngwf_basis%spheres(loc_ngwf)%offset
       npoints = ngwf_basis%spheres(loc_ngwf)%n_ppds_sphere * pub_cell%n_pts
       lastp = offset + npoints -1

       ! calculate product of NGWF with itself
       my_spd_sum(max_l+1,loc_ngwf) = sum(ngwfs_on_grid(offset:lastp)*&
                                      ngwfs_on_grid(offset:lastp))

       ! calculate fractional coordinates of ngwf centre
       frac%X = (centre .DOT. pub_cell%b1) / two_pi * real(pub_cell%total_pt1,kind=DP)
       frac%Y = (centre .DOT. pub_cell%b2) / two_pi * real(pub_cell%total_pt2,kind=DP)
       frac%Z = (centre .DOT. pub_cell%b3) / two_pi * real(pub_cell%total_pt3,kind=DP)

       ! calculate displacement of centre from nearest grid point
       disp = mod(frac%X,1.0_DP) / pub_fftbox%total_pt1 * pub_fftbox%a1 + &
              mod(frac%Y,1.0_DP) / pub_fftbox%total_pt2 * pub_fftbox%a2 + &
              mod(frac%Z,1.0_DP) / pub_fftbox%total_pt3 * pub_fftbox%a3

       ! calculate first grid point of tightbox wrt cell
       call basis_location_func_wrt_cell(sx,sy,sz,ngwf_basis%tight_boxes(loc_ngwf))

       ! calculate location of tightbox wrt FFTbox and displacement of centre wrt
       ! FFTbox depending whether FFTbox side coincides with sim cell
       if(pub_fftbox%coin1) then
          sx = 1
          disp%X = centre%X
       else
          sx = (pub_fftbox%total_pt1-1)/2 - floor(frac%X) + sx
          disp%X = fbcentre%X + disp%X
       endif
       if(pub_fftbox%coin2) then
          sy = 1
          disp%Y = centre%Y
       else
          sy = (pub_fftbox%total_pt2-1)/2 - floor(frac%Y) + sy
          disp%Y = fbcentre%Y + disp%Y
       endif
       if(pub_fftbox%coin3) then
          sz = 1
          disp%Z = centre%Z
       else
          sz = (pub_fftbox%total_pt3-1)/2 - floor(frac%Z) + sz
          disp%Z = fbcentre%Z + disp%Z
       endif

       ! loop over angular momentum l
       do l=0,max_l

          ! generate bessel function roots for current NGWF radius and l
          qnl(:) = sw_bessel_zeros(:,l) / radius

          ! loop over azimuthal angular momentum m
          do m=-l,+l

             ! write l,m to file if requested
             if (sw2file) write(fileunit,'(/2i4,1x)',advance='no') l, m

             ! loop over n
             do n=1,max_n

                ! calculate SW in an FFTbox
                call sw_recp_generate(sw_fftbox,sw_work,l,m,qnl(n),radius,disp)

                ! extract SW from FFTbox into ppds
                call basis_extract_function_from_box(sw_on_grid, &
                     pub_fftbox%total_ld1,pub_fftbox%total_ld2, &
                     pub_fftbox%total_pt3,sw_fftbox, &
                     ngwf_basis%spheres(loc_ngwf), &
                     ngwf_basis%tight_boxes(loc_ngwf),sx,sy,sz,1)

                ! perform integrals
                tot  = sum(sw_on_grid(1:npoints) * ngwfs_on_grid(offset:lastp))
                norm = sum(sw_on_grid(1:npoints) * sw_on_grid(1:npoints))

                ! sum the square of the integral after normalising spherical wave
                my_spd_sum(l,loc_ngwf) = my_spd_sum(l,loc_ngwf) + tot**2/norm

                ! write SW coefficient to file if requested
                if (sw2file) write(fileunit,'(f10.5)',advance='no') &
                     tot*sqrt(pub_cell%weight/norm)

             enddo  ! end n loop
          enddo     ! end m loop
       enddo        ! end l loop

       ! calculate total and spilling
       my_spd_sum(max_l+2,loc_ngwf) = sum(my_spd_sum(:max_l,loc_ngwf))
       my_spd_sum(max_l+1,loc_ngwf) = my_spd_sum(max_l+1,loc_ngwf) - &
            my_spd_sum(max_l+2,loc_ngwf)

       if (sw2file) then
          close(unit=fileunit,iostat=ierr)
          call utils_open_unit_check('properties_ngwfs_char','fileunit',ierr)
       end if

    enddo  ! end ngwf loop

    ! weight integrals
    my_spd_sum(:,:) = my_spd_sum(:,:) * pub_cell%weight

    ! deallocate sw_on_grid array
    deallocate(sw_on_grid,stat=ierr)
    call utils_dealloc_check('properties_ngwfs_char','sw_on_grid',ierr)
    ! deallocate workspace
    deallocate(sw_work,stat=ierr)
    call utils_dealloc_check('properties_ngwfs_char','sw_work',ierr)
    deallocate(sw_fftbox,stat=ierr)
    call utils_dealloc_check('properties_ngwfs_char','sw_fftbox',ierr)

    ! deallocate qnl array
    deallocate(qnl,stat=ierr)
    call utils_dealloc_check('properties_ngwfs_char','qnl',ierr)

    ! allocate workspace
    allocate(snd_spd_sum(0:max_l+2,ngwf_basis%max_on_node),stat=ierr)
    call utils_alloc_check('properties_ngwfs_char','snd_spd_sum',ierr)
    allocate(spd_sum(0:max_l+2,ngwf_basis%num),stat=ierr)
    call utils_alloc_check('properties_ngwfs_char','spd_sum',ierr)

    ! gather my_spd_sums to root node
    ! this is effectively an MPI_GATHERV
    if (pub_on_root) then
       offset = 1
       do node=0,pub_total_num_nodes-1
          if (node .eq. pub_root_node_id) then
             spd_sum(:,offset:offset+ngwf_basis%node_num-1) = &
                my_spd_sum(:,:ngwf_basis%node_num)
          else
             call comms_recv(node,snd_spd_sum)
             spd_sum(:,offset:offset+ngwf_basis%num_on_node(node)-1) = &
                snd_spd_sum(:,:ngwf_basis%num_on_node(node))
          endif
          offset = offset + ngwf_basis%num_on_node(node)
       enddo
    else
       call comms_send(pub_root_node_id,my_spd_sum)
    endif
    call comms_barrier

    ! write out table
    if (pub_on_root) then
       write(stdout,'(/a,i2.2,a)') '************************ NGWF s/p/d/f Character Analysis &
                             &*********** (n=',max_n,') ****'
       write(stdout,'(a1,3x,2x,1x,a6,1x,a6,4x,a1,2x,4a9,4x,a9,3x,a1)') &
          '|','Atom','NGWF',':','%s','%p','%d','%f','spilling','|'
       do atom=1,pub_cell%nat
          ngwf = ngwf_basis%first_on_atom(pub_distr_atom(atom)) - 1
          do ngwf_on_atom=1,ngwf_basis%num_on_atom(pub_distr_atom(atom))
             ngwf = ngwf + 1           ! original NGWF counting
             if (ngwf_on_atom .eq. 1) then
                write(stdout,'(a1,3x,a2,1x,i6,1x,i6,4x,a1,2x,4f9.3,4x,e9.2,3x,a1)') &
                   '|',elements(atom)%symbol,atom,ngwf_on_atom,':',&
                   100.0_DP*spd_sum(:3,ngwf)/spd_sum(max_l+2,ngwf),&
                   spd_sum(max_l+1,ngwf),'|'
             else
                write(stdout,'(a1,3x,2x,1x,6x,1x,i6,4x,a1,2x,4f9.3,4x,e9.2,3x,a1)') &
                   '|',ngwf_on_atom,':',&
                   100.0_DP*spd_sum(:3,ngwf)/spd_sum(max_l+2,ngwf),&
                   spd_sum(max_l+1,ngwf),'|'
             endif
          enddo
       enddo
       write(stdout,'(a/)') '****************************************&
                             &****************************************'
    endif

    ! deallocate workspace
    deallocate(snd_spd_sum,stat=ierr)
    call utils_dealloc_check('properties_ngwfs_char','snd_spd_sum',ierr)
    deallocate(spd_sum,stat=ierr)
    call utils_dealloc_check('properties_ngwfs_char','spd_sum',ierr)
    deallocate(my_spd_sum,stat=ierr)
    call utils_dealloc_check('properties_ngwfs_char','my_spd_sum',ierr)

    ! stop timer
    call timer_clock('properties_ngwfs_char',2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: Leaving properties_ngwfs_char'
#endif

  end subroutine properties_ngwfs_char


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine properties_polarisation_simcell(ngwfs_on_grid,ngwf_basis,&
       elements,denskern,overlap)

    !==========================================================================!
    ! Calculates dipole moment of a molecule by integrating over the           !
    ! simulation cell (cannot handle atoms moving across cell boundary)        !
    !                                                                          !
    ! For large systems this will be slower than properties_polarisation and   !
    ! will also need a larger stack size                                       !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! ngwfs_on_grid (input)    : NGWFs in ppd representation on this node      !
    ! ngwf_basis (input)       : The function basis of NGWFs                   !
    ! elements (input)         : Element array for atoms                       !
    ! denskern (input)         : Density kernel                                !
    ! overlap (input)          : Overlap matrix                                !
    !--------------------------------------------------------------------------!
    ! Written by Mark Robinson, July 2008                                      !
    ! Modified 22/01/2011 by Nicholas Hine to use cell_grid_real_pt routine to !
    ! save on storage.                                                         !
    !==========================================================================!

    use cell_grid, only: pub_fine_grid, cell_grid_real_pt
    use comms, only: comms_reduce, pub_on_root
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use ion, only: element
    use density, only: density_on_grid
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3
    use timer, only: timer_clock

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in) :: ngwfs_on_grid(ngwf_basis%n_ppds*pub_cell%n_pts)
    type(ELEMENT), intent(in) :: elements(pub_cell%nat)
    type(SPAM3), intent(in)   :: denskern(pub_cell%num_spins)
    type(SPAM3), intent(in)   :: overlap

    ! Local Variables
    real(kind=DP) :: density_fine(pub_fine_grid%ld1,pub_fine_grid%ld2,&
         pub_fine_grid%max_slabs12,pub_cell%num_spins)
    real(kind=DP) :: d(3),q(3),t(3),rpt(3)
    integer       :: islab12,i2,i1,is

    ! Start timer
    call timer_clock('properties_polarisation_simcell',1)

    ! calculate dipole moment for each coordinate
    d = 0.0_DP ! ddor: Set dipole moment initially to zero
    ! calculate density on fine grid
    call density_on_grid(density_fine,pub_fine_grid, &
         denskern,overlap,ngwfs_on_grid,ngwf_basis)

    ! loop over fine grid points and over spins
    do islab12=1,pub_fine_grid%num_my_slabs12
       do i2=1,pub_fine_grid%n2
          do i1=1,pub_fine_grid%n1
             do is=1,pub_cell%num_spins
                 ! integrate density*position over simulation cell
                 call cell_grid_real_pt(rpt,i1,i2,islab12,pub_fine_grid)
                 d(1:3) = d(1:3) + density_fine(i1,i2,islab12,is)*rpt(1:3)
             end do
          end do
       enddo
    enddo

    ! sum contributions from all nodes
    call comms_reduce('SUM',d)

    ! weight integral
    d = d * pub_fine_grid%weight

    ! multiplication factor depending on number of spins
    if (pub_cell%num_spins == 1) d(:) = 2.0_DP * d(:)

    ! calculate dipole moment from ions
    q(1) = sum(elements(:)%ion_charge * elements(:)%centre%X)
    q(2) = sum(elements(:)%ion_charge * elements(:)%centre%Y)
    q(3) = sum(elements(:)%ion_charge * elements(:)%centre%Z)

    ! calculate total dipole moment
    d = -d
    t = q + d

    ! print pretty table
    if(pub_on_root) then
       write(stdout,'(/a)') '========================== Dipole Moment &
                            &Calculation (SimCell) ================='
       write(stdout,'(a34,4x,a)') 'Type:','Single Molecule'
       write(stdout,*)
       write(stdout,'(a34,4x,a4,1x,f24.12)') 'Electronic dipole moment &
            &(e.bohr):','dx =',d(1)
       write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',d(2)
       write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',d(3)
       write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
            sqrt(sum(d(:)**2))
       write(stdout,*)
       write(stdout,'(a34,4x,a4,1x,f24.12)') 'Ionic dipole moment (e.bohr):', &
            'dx =',q(1)
       write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',q(2)
       write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',q(3)
       write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
            sqrt(sum(q(:)**2))
       write(stdout,*)
       write(stdout,'(a34,4x,a4,1x,f24.12)') 'Total dipole moment (e.bohr):', &
            'dx =',t(1)
       write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',t(2)
       write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',t(3)
       write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
            sqrt(sum(t(:)**2))
       write(stdout,*)
       write(stdout,'(a)') '==================================================&
                           &=============================='
    endif

    ! Stop timer
    call timer_clock('properties_polarisation_simcell',2)

  end subroutine properties_polarisation_simcell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine properties_polarisation(rep,ngwf_basis,proj_basis,elements, &
       denskern,axis,dipole)

    !==========================================================================!
    ! Calculates dipole moment of a molecule by computing matrix elements      !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! rep (input)        : NGWF representation for valence                     !
    ! ngwf_basis (input) : Function basis type describing the NGWFs            !
    ! elements (input)   : Element array for atoms                             !
    ! denskern (input)   : Density kernel                                      !
    !--------------------------------------------------------------------------!
    ! Written by Mark Robinson, April 2008                                     !
    ! Modified by David O'Regan in March 2009 so that it can be                !
    ! called for a single direction only.                                      !
    ! Extended for PAW sphere terms by Nicholas Hine in November 2011.         !
    !==========================================================================!

    use augmentation, only: aug_projector_denskern, augmentation_pos
    use basis, only : basis_location_func_wrt_cell, basis_ket_start_wrt_fftbox
    use comms, only: pub_on_root, pub_my_node_id
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(*), operator(+)
    use integrals, only: integrals_pos
    use ion, only: element
    use md_ions, only: md_ionic_polarisation, md_elec_polarisation, &
         md_total_polarisation
    use ngwf_representation, only: NGWF_REP
    use parallel_strategy, only: pub_orig_atom
    use paw, only: paw_position_operator, paw_projector_overlap
    use rundat, only: print_qc, pub_aug, task
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: sparse_create, sparse_get_element, SPAM3, sparse_trace, &
         sparse_destroy, sparse_put_element, sparse_element_exists, &
         sparse_expand, PATTERN_ALTERNATE
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check, utils_unit
!CW code of David O'regan
    use rundat, only : pub_write_polarisation_plot
    use constants, only: ANGSTROM, UP, DN, max_spins
    use cell_grid, only: pub_fine_grid
    use sparse, only: sparse_product,sparse_transpose,sparse_axpy,sparse_scale
    use density, only: density_on_grid
    use visual, only: visual_scalarfield
!ENDCW
    ! Arguments
    type(NGWF_REP), intent(in) :: rep
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    type(FUNC_BASIS), intent(in) :: proj_basis
    type(ELEMENT), intent(in)            :: elements(pub_cell%nat)
    type(SPAM3), intent(in)              :: denskern(pub_cell%num_spins)
    ! ddor: When axis is present, only the electronic contribution to the dipole
    !       moment in direction axis is calculated, and this is output to dipole.
    integer, optional, intent(in) :: axis
    real(kind=DP), optional, intent(out) :: dipole

    ! Local Variables
    type(SPAM3)   :: r_elements(3)
    type(POINT)   :: a1,a2,a3
    type(POINT)   :: r
    real(kind=DP) :: d(3),q(3),t(3)
    real(kind=DP) :: R_fft(3),r_el,o_el
    integer       :: cs1,cs2,cs3
    integer       :: bs1,bs2,bs3
    integer       :: xyz,axmin,axmax,is
    integer       :: jngwf,loc_ingwf,ingwf
!CW for piece of code David
    integer       :: ierr ! memory allocation error flag
    character(len=1) :: dir_string
    real(kind=DP), allocatable :: polarisation_density_fine(:,:,:,:)
    type(SPAM3) :: kr_elements, srk_elements, krs_elements(max_spins)
!END CW

    ! Start timer
    call timer_clock('properties_polarisation',1)

    ! ddor: Allocate only one matrix if axis is specified
    if (present(axis)) then
       axmin = axis
       axmax = axis
    else
       axmin = 1
       axmax = 3
    end if

    ! Allocate sparse matrices
    do xyz=axmin,axmax
       call sparse_create(r_elements(xyz),rep%overlap)
    end do

    ! calculate matrix elements < phi_a | r_op - R_fft | phi_b >
    call integrals_pos(r_elements,&
         rep%ngwfs_on_grid,ngwf_basis,rep%ngwfs_on_grid,ngwf_basis,1,axis)

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
    do loc_ingwf=1,ngwf_basis%num_on_node(pub_my_node_id)
       ingwf = loc_ingwf + ngwf_basis%first_on_node(pub_my_node_id) - 1

       ! Find start of tightbox of ingwf
       call basis_location_func_wrt_cell(cs1,cs2,cs3, &
            ngwf_basis%tight_boxes(loc_ingwf))

       ! Find vector to origin of FFTbox
       r = real(cs1-bs1,kind=DP) * a1 &
            + real(cs2-bs2,kind=DP) * a2 &
            + real(cs3-bs3,kind=DP) * a3
       R_fft(1) = r%X
       R_fft(2) = r%Y
       R_fft(3) = r%Z

       ! Loop over all row ngwfs
       do jngwf=1,ngwf_basis%num

          ! Test if these NGWFs overlap
          if (.not.sparse_element_exists(rep%overlap,jngwf,ingwf)) cycle

          ! Extract overlap element
          call sparse_get_element(o_el,rep%overlap,jngwf,ingwf)

          ! Extract element from r_elements and shift by R_fft*o_el
          ! ddor: get elements for one direction only if axis is specified
          do xyz=axmin,axmax
             call sparse_get_element(r_el,r_elements(xyz),jngwf,ingwf)
             r_el = R_fft(xyz) * o_el + r_el
             call sparse_put_element(r_el,r_elements(xyz),jngwf,ingwf)
          end do

       end do
    end do

    ! ndmh: in PAW/USP formalism, add contribution from sphere part
    if (pub_aug) call augmentation_pos(r_elements,proj_basis,rep%sp_overlap, &
         axis=axis)

    ! Calculate dipole moment from electrons by Tr(KR)
    d = 0.0_DP

    ! ddor: Calculate for one direction only if axis is specified
    axiscase: if (present(axis)) then

       do is=1,pub_cell%num_spins
          d(axis) = d(axis) - sparse_trace(denskern(is),r_elements(axis))
       end do

       ! multiplication factor depending on number of spins
       if (pub_cell%num_spins == 1) d(:) = 2.0_DP * d(:)

       ! calculate dipole moment from ions
       q(1) = sum(elements(:)%ion_charge * elements(:)%centre%X)
       q(2) = sum(elements(:)%ion_charge * elements(:)%centre%Y)
       q(3) = sum(elements(:)%ion_charge * elements(:)%centre%Z)

       ! calculate total dipole moment
       t = q + d

       ! print pretty table
       if (pub_on_root) then
          write(stdout,'(/a)') '========================== Dipole Moment &
                               &Calculation ==========================='
          write(stdout,'(a34,4x,a)') 'Type:','Single Molecule'
          write(stdout,'(a34,4x,a2,1x,i1,1x,a3,1x,f24.18)') 'Electronic dipole moment (e.bohr):',&
                               &'d(',axis,') =',d(axis)
          write(stdout,'(a34,4x,a2,1x,i1,1x,a3,1x,f24.18)') 'Ionic dipole moment (e.bohr):',&
                               &'q(',axis,') =',q(axis)
          write(stdout,'(a34,4x,a2,1x,i1,1x,a3,1x,f24.18)') 'Total dipole moment (e.bohr):',&
                               &'t(',axis,') =',t(axis)
          write(stdout,'(a)') '==================================================&
                           &=============================='
       end if

       dipole = d(axis)

       ! deallocate sparse matrices
       call sparse_destroy(r_elements(axis))

    else

       do xyz=1,3
!CW piece of code of David
          ! ddor: Write out the electronic "polarisation density" defined as
          !     : -1.0 * phi_a(r) K^ab <r>^i_bc S^cd phi_d(r) to .cube file
          !     : N.B. only the change in this quantity between different
          !     : densities is of any interest, and NOT the quantity itself.
          if (pub_write_polarisation_plot) then

             ! ddor: Allocate polarisation density slabs
             allocate(polarisation_density_fine(&
                  pub_fine_grid%ld1,pub_fine_grid%ld2, &
                  pub_fine_grid%max_slabs12,pub_cell%num_spins),stat=ierr)
             call utils_alloc_check('properties_polarisation',&
                  'polarisation_density_fine',ierr)

             ! ddor: Create matrices with sparsity KSK for symmetrisation
             call sparse_create(kr_elements,denskern(1),r_elements(xyz))
             do is=1,pub_cell%num_spins
                call sparse_create(krs_elements(is),kr_elements,denskern(1))
             enddo
             call sparse_create(srk_elements,krs_elements(1))

             ! ddor: Loop over spins
             do is=1,pub_cell%num_spins

                ! ddor: Compute the dipole moment
                call sparse_product(kr_elements,denskern(is),r_elements(xyz))
                d(xyz) = d(xyz) - sparse_trace(kr_elements)

                ! ddor: Compute the symmetrised matrix elements of K<r>S^^
                call sparse_product(krs_elements(is),kr_elements,&
                     rep%inv_overlap)
                call sparse_transpose(srk_elements,krs_elements(is))
                call sparse_axpy(krs_elements(is),srk_elements,1.0_DP)
                if (pub_cell%num_spins .ne. 1) & ! Scale for spin if necessary
                     call sparse_scale(krs_elements(is),0.5_DP)

             enddo

             ! ddor: Calculate data-parallelised polarisation density
             !CW modfication
             call density_on_grid(polarisation_density_fine,pub_fine_grid, &
                  krs_elements,rep%ngwf_overlap,rep%ngwfs_on_grid,ngwf_basis)
                  !, &
                  !rep%ngwfs_on_grid,ngwf_basis)
             !END CW

             ! ddor: spin polarisation: calculate total polarisation density 
             !     : and polarisation spin density
             if (pub_cell%num_spins == 2) then
                polarisation_density_fine(:,:,:,1) = &
                     polarisation_density_fine(:,:,:,UP) + &
                     polarisation_density_fine(:,:,:,DN)
                polarisation_density_fine(:,:,:,2) = &
                     polarisation_density_fine(:,:,:,UP) - &
                     2.0_DP * polarisation_density_fine(:,:,:,DN)
             end if

             ! ddor: Direction string for file name
             if (xyz==1) dir_string='x'
             if (xyz==2) dir_string='y'
             if (xyz==3) dir_string='z'

             ! cks: output density in plot format file
             ! vm: output density in Angstrom rather than in Bohr
             ! jd: ... but leave the unit conversion to visual_scalarfield
             ! ddor: The minus sign for electrons occurs here
             call visual_scalarfield( &
                  polarisation_density_fine(:,:,:,1), pub_fine_grid, &
                  'Electronic polarisation density (in e/ang^2) for:', &
                  '_polarisation_density_'//trim(dir_string), &
                  elements, -1.0_DP * ANGSTROM**2)
             if (pub_cell%num_spins == 2) call visual_scalarfield( &
                  polarisation_density_fine(:,:,:,2), pub_fine_grid, &
                  'Electronic polarisation spin density (in e/ang^2) for:', &
                  '_polarisation_spindensity_'//trim(dir_string), &
                  elements, -1.0_DP * ANGSTROM**2)

             !ddor: Destroy matrices
             call sparse_destroy(srk_elements) 
             do is=pub_cell%num_spins,1,-1
                call sparse_destroy(krs_elements(is))
             enddo
             call sparse_destroy(kr_elements) 

             ! ddor: Deallocate polarisation density slabs
             deallocate(polarisation_density_fine,stat=ierr)
             call utils_dealloc_check('properties_polarisation',&
                  'polarisation_density_fine',ierr)

          else
!END CW
           do is=1,pub_cell%num_spins
              d(xyz) = d(xyz) - sparse_trace(denskern(is),r_elements(xyz))
           end do
!CW 
          endif
!END CW
       end do

       ! Multiplication factor depending on number of spins
       if (pub_cell%num_spins == 1) d(:) = 2.0_DP * d(:)

       ! Calculate dipole moment from ions
       q(1) = sum(elements(:)%ion_charge * elements(:)%centre%X)
       q(2) = sum(elements(:)%ion_charge * elements(:)%centre%Y)
       q(3) = sum(elements(:)%ion_charge * elements(:)%centre%Z)

       ! Calculate total dipole moment
       t = q + d

       ! Print pretty table
       if(pub_on_root) then
          write(stdout,'(/a)') '========================== Dipole Moment &
                               &Calculation ==========================='
          write(stdout,'(a34,4x,a)') 'Type:','Single Molecule'
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Electronic dipole moment &
               &(e.bohr):','dx =',d(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',d(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',d(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(d(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Ionic dipole moment (e.bohr):', &
               'dx =',q(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',q(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',q(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(q(:)**2))
          write(stdout,*)
          write(stdout,'(a34,4x,a4,1x,f24.12)') 'Total dipole moment (e.bohr):', &
               'dx =',t(1)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dy =',t(2)
          write(stdout,'(34x,4x,a4,1x,f24.12)') 'dz =',t(3)
          write(stdout,'(a34,4x,4x,1x,f24.12)') 'Magnitude (e.bohr):', &
               sqrt(sum(t(:)**2))
          write(stdout,*)
          if (print_qc) then
             write(stdout,'(a,f24.12)') &
                  '<QC>        [elec_polarisation_d1]: ',d(1)
             write(stdout,'(a,f24.12)') &
                  '<QC>        [elec_polarisation_d2]: ',d(2)
             write(stdout,'(a,f24.12)') &
                  '<QC>        [elec_polarisation_d3]: ',d(3)
             write(stdout,'(a,f24.12)') &
                  '<QC>       [ionic_polarisation_d1]: ',q(1)
             write(stdout,'(a,f24.12)') &
                  '<QC>       [ionic_polarisation_d2]: ',q(2)
             write(stdout,'(a,f24.12)') &
                  '<QC>       [ionic_polarisation_d3]: ',q(3)
             write(stdout,'(a,f24.12)') &
                  '<QC>       [total_polarisation_d1]: ',t(1)
             write(stdout,'(a,f24.12)') &
                  '<QC>       [total_polarisation_d2]: ',t(2)
             write(stdout,'(a,f24.12)') &
                  '<QC>       [total_polarisation_d3]: ',t(3)
             write(stdout,*)
          end if
          write(stdout,'(a)') '===============================================&
                              &================================='

          ! smmd: print data in md file
          if (task == 'MOLECULARDYNAMICS') then
             md_elec_polarisation(:) = d(:)
             md_ionic_polarisation(:) = q(:)
             md_total_polarisation(:) = t(:)
          end if

       endif

       ! Deallocate sparse matrices
       do xyz=1,3
          call sparse_destroy(r_elements(xyz))
       end do

    end if axiscase

    ! Stop timer
    call timer_clock('properties_polarisation',2)

  end subroutine properties_polarisation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine properties_spread(ngwfs_on_grid,ngwf_basis, &
       elements,inv_overlap,overlap,nonorthogonal)

    !==========================================================================!
    ! Calculates the spread (second central moment) of the NGWFs.              !
    !--------------------------------------------------------------------------!
    ! Arguments:                                                               !
    ! ngwfs_on_grid (input)    : NGWFs in ppd representation on this node      !
    ! ngwf_basis (input)       : Function basis type describing the NGWFs      !
    ! elements (input)         : Element array for atoms                       !
    ! inv_overlap (input)      : Inverse overlap matrix                        !
    ! overlap (input)          : Overlap matrix                                !
    !--------------------------------------------------------------------------!
    ! Written by David O'Regan in September 2009,                              !
    ! based on properties_polarisation by M. Robinson.                         !
    ! Automatic arrays removed by Nicholas Hine, December 2009.                !
    !==========================================================================!

    use basis, only : basis_location_func_wrt_cell, basis_ket_start_wrt_fftbox
    use comms, only: pub_on_root, pub_my_node_id, comms_reduce, comms_barrier
    use constants, only: DP, stdout
    use function_basis, only: FUNC_BASIS
    use geometry, only: POINT, operator(*), operator(+)
    use integrals, only: integrals_brappd_ketppd, integrals_pos
    use ion, only: element
    use parallel_strategy, only : pub_distr_atom
    use simulation_cell, only: pub_cell, pub_fftbox
    use sparse, only: sparse_create, sparse_get_element, SPAM3, sparse_trace, &
         sparse_destroy, sparse_put_element, sparse_product, &
         sparse_axpy, sparse_expand, sparse_element_exists
    use timer, only: timer_clock
    use utils, only: utils_alloc_check, utils_dealloc_check

    ! Arguments
    type(FUNC_BASIS), intent(in) :: ngwf_basis
    real(kind=DP), intent(in)            :: ngwfs_on_grid( &
         ngwf_basis%n_ppds*pub_cell%n_pts)
    type(ELEMENT), intent(in)            :: elements(pub_cell%nat)
    type(SPAM3), intent(in)              :: inv_overlap
    type(SPAM3), intent(in)              :: overlap
    ! ddor: Whether or not the input NGWFs are nonorthogonal
    logical, intent(in)                  :: nonorthogonal

    ! Local Variables
    type(SPAM3)   :: r_elements(3)
    type(SPAM3)   :: r2_elements(3)
    type(SPAM3)   :: rs(3), rsr(3), thespread(3)
    type(SPAM3)   :: closetounity
    type(POINT)   :: a1,a2,a3
    type(POINT)   :: r
    real(kind=DP) :: R_fft(3),r_el(3),r2_el(3),o_el
    integer       :: cs1,cs2,cs3
    integer       :: bs1,bs2,bs3
    integer       :: xyz
    integer       :: loc_ingwf, ingwf, jngwf
    integer       :: atom, ngwf_on_atom
    integer       :: ierr
    real(kind=DP), allocatable :: covariant(:,:)
    real(kind=DP), allocatable :: invariant(:,:)

    ! Start timer
    call timer_clock('properties_spread',1)

    ! ndmh: allocate temporary arrays
    allocate(covariant(ngwf_basis%num,4),stat=ierr)
    call utils_alloc_check('properties_spread','covariant',ierr)
    allocate(invariant(ngwf_basis%num,4),stat=ierr)
    call utils_alloc_check('properties_spread','invariant',ierr)

    do xyz=1,3
       call sparse_create(r_elements(xyz),overlap)
       call sparse_create(r2_elements(xyz),overlap)
       if (nonorthogonal) then
          call sparse_create(rs(xyz),r_elements(xyz),inv_overlap)
          call sparse_create(rsr(xyz),rs(xyz),r_elements(xyz))
          call sparse_create(thespread(xyz),rsr(xyz),inv_overlap)
       else
          call sparse_create(thespread(xyz),overlap,overlap)
       endif
    enddo

    if (.not. nonorthogonal) then
       call sparse_create(closetounity,overlap)
       call integrals_brappd_ketppd(closetounity, &  !input-output
            ngwfs_on_grid, ngwf_basis, &
            ngwfs_on_grid, ngwf_basis )
    endif

    covariant = 0.0_DP
    invariant = 0.0_DP

    ! Calculate vectors between grid points
    a1 = (1.0_DP / pub_cell%total_pt1) * pub_cell%a1
    a2 = (1.0_DP / pub_cell%total_pt2) * pub_cell%a2
    a3 = (1.0_DP / pub_cell%total_pt3) * pub_cell%a3

    ! Determine where tightbox begins wrt fftbox
    call basis_ket_start_wrt_fftbox(bs1,bs2,bs3,pub_fftbox%total_pt1,&
         pub_fftbox%total_pt2,pub_fftbox%total_pt3)

    ! calculate matrix elements < phi_a | r_op - R_fft | phi_b >
    ! ndmh: use integrals_pos routine
    call integrals_pos(r_elements,&
         ngwfs_on_grid,ngwf_basis,ngwfs_on_grid,ngwf_basis,1)

    ! calculate matrix elements < phi_a | ( r_op - R_fft )^2 | phi_b >
    call integrals_pos(r2_elements,&
         ngwfs_on_grid,ngwf_basis,ngwfs_on_grid,ngwf_basis,2)

    ! Add R_fft * S_ab to matrix elements, taking into account atoms moving
    ! across cell boundary (NOT IMPLEMENTED YET)
    do loc_ingwf=1,ngwf_basis%num_on_node(pub_my_node_id)
       call basis_location_func_wrt_cell(cs1,cs2,cs3, &
            ngwf_basis%tight_boxes(loc_ingwf))
       r = real(cs1-bs1,kind=DP) * a1 &
            + real(cs2-bs2,kind=DP) * a2 &
            + real(cs3-bs3,kind=DP) * a3
       R_fft(1) = r%X
       R_fft(2) = r%Y
       R_fft(3) = r%Z

       ! ddor: add shifts to matrix elements of position
       !       operator and position operator squared
       ingwf = loc_ingwf + ngwf_basis%first_on_node(pub_my_node_id) - 1
       do jngwf=ingwf,ngwf_basis%num

          ! Test if these NGWFs overlap
          if (.not.sparse_element_exists(overlap,jngwf,ingwf)) cycle

          if (nonorthogonal) then
             call sparse_get_element(o_el,overlap,jngwf,ingwf)
          else
             call sparse_get_element(o_el,closetounity,jngwf,ingwf)
          endif
          do xyz=1,3
             call sparse_get_element(r_el(xyz),r_elements(xyz),jngwf,ingwf)
             call sparse_get_element(r2_el(xyz),r2_elements(xyz),jngwf,ingwf)
             !ddor: <a| r | b>  =  <a| (r - R) | b> +  R <a | b>
             r_el(xyz) = (R_fft(xyz) * o_el) + r_el(xyz)
             !ddor: <a| r** 2 |b>  =  <a| (r - R)**2 | b> +  2 <a| R.r |b> -
             !                        <a| R**2 |b>
             !                     =  <a| (r - R)**2 | b> +  2 R.<a| r |b> -
             !                        R**2 <a | b>
             r2_el(xyz) = r2_el(xyz) + &
                  (2.0_DP * R_fft(xyz) * r_el(xyz)) - &
                  (R_fft(xyz) * R_fft(xyz) * o_el )
             call sparse_put_element(r2_el(xyz),r2_elements(xyz),jngwf,ingwf)
             call sparse_put_element(r_el(xyz),r_elements(xyz),jngwf,ingwf)
          enddo
       enddo
    enddo

    call comms_barrier

    do xyz=1,3

       if (nonorthogonal) then

          ! ddor: <a| r |b> S^bc
          call sparse_product(rs(xyz),r_elements(xyz),inv_overlap)
          ! ddor: <a| r |b> S^bc <c| r | d>
          call sparse_product(rsr(xyz),rs(xyz),r_elements(xyz))
          ! ddor: The negative of the covariant NGWF spread matrix
          !       <a| r |b> S^bc <c| r | d> - <a| r^2 |d>
          call sparse_axpy(rsr(xyz),r2_elements(xyz),-1.0_DP)

          ! ddor: The negative of the invariant NGWF spread matrix
          !       ( <a| r |b> S^bc <c| r | d> - <a| r^2 |d> ) S^de
          call sparse_product(thespread(xyz),rsr(xyz),inv_overlap)

       else
          ! ddor: Case of orthogonalised NGWFs

          ! ddor: <a| r |c> <c| r |b>
          call sparse_product(thespread(xyz),r_elements(xyz),r_elements(xyz))
          ! ddor: The negative of the NGWF spread matrix
          !       <a| r |c> <c| r | b> - <a| r^2 |b>
          call sparse_axpy(thespread(xyz),r2_elements(xyz),-1.0_DP)

       endif

       ! ddor: Collect diagonal elements of Spread matrices
       do loc_ingwf=1,ngwf_basis%num_on_node(pub_my_node_id)
          ingwf = loc_ingwf + ngwf_basis%first_on_node(pub_my_node_id) - 1
          if (nonorthogonal) then
             call sparse_get_element(covariant(ingwf,xyz),&
                  rsr(xyz),ingwf,ingwf)
             covariant(ingwf,xyz) = -1.0_DP * covariant(ingwf,xyz)
          endif
          call sparse_get_element(invariant(ingwf,xyz),&
               thespread(xyz),ingwf,ingwf)
          invariant(ingwf,xyz) = -1.0_DP * invariant(ingwf,xyz)
       enddo

    enddo

    ! ddor: Sum of spreads in each direction
    !       Could also calculate spread anisotropy and output Wannier centres
    do loc_ingwf=1,ngwf_basis%num_on_node(pub_my_node_id)
       ingwf = loc_ingwf + ngwf_basis%first_on_node(pub_my_node_id) - 1
       if (nonorthogonal) then
          covariant(ingwf,4) = SUM(covariant(ingwf,1:3))
       endif
       invariant(ingwf,4) = SUM(invariant(ingwf,1:3))
    enddo

    call comms_reduce('SUM',invariant)

    ! ddor: If we need to print out both covariant and invariant spreads
    if (nonorthogonal) then

       call comms_reduce('SUM',covariant)

       ! Print table of spreads, adapted from properties_ngwfs_char by mr
       if (pub_on_root) then
          write(stdout,'(/a)') '****************************************&
               &****************************************'
          write(stdout,'(a1,2x,a4,2x,a4,1x,8x,a1,1x,a25,1x,a25,4x,a1)') &
               &'|','Atom','NGWF',':','Covariant Spread (bohr^2)',&
               &'Invariant Spread (bohr^2)','|'
          do atom=1,pub_cell%nat
             ingwf = ngwf_basis%first_on_atom(pub_distr_atom(atom))
             do ngwf_on_atom=1,ngwf_basis%num_on_atom(pub_distr_atom(atom))
                if (ngwf_on_atom==1) then
                   write(stdout,'(a1,2x,a2,2x,i6,1x,i6,2x,a1,2x,f24.10,2x,&
                        &f24.10,4x,a1)') &
                        &'|',elements(atom)%symbol,atom,ngwf_on_atom,':',&
                        covariant(ingwf,4),invariant(ingwf,4),'|'
                else
                   write(stdout,'(a1,2x,2x,2x,6x,1x,i6,2x,a1,2x,f24.10,2x,&
                        &f24.10,4x,a1)') &
                        &'|',ngwf_on_atom,':',&
                        covariant(ingwf,4),invariant(ingwf,4),'|'
                endif
                ingwf = ingwf + 1           ! original NGWF counting
             enddo
          enddo

          write(stdout,'(a/)') '****************************************&
               &****************************************'
       endif

    else

       ! Print table of spreads, adapted from properties_ngwfs_char by mr
       if (pub_on_root) then
          write(stdout,'(/a)') '****************************************&
               &****************************************'
          write(stdout,'(a1,2x,a4,2x,a4,1x,8x,a1,1x,a25,4x,a1)') &
               &'|','Atom','NGWF',':','Spread (bohr^2)','|'
          do atom=1,pub_cell%nat
             ingwf = ngwf_basis%first_on_atom(pub_distr_atom(atom))
             do ngwf_on_atom=1,ngwf_basis%num_on_atom(pub_distr_atom(atom))
                if (ngwf_on_atom==1) then
                   write(stdout,'(a1,2x,a2,2x,i6,1x,i6,2x,a1,2x,&
                        &f24.10,4x,a1)') &
                        &'|',elements(atom)%symbol,atom,ngwf_on_atom,':',&
                        invariant(ingwf,4),'|'
                else
                   write(stdout,'(a1,2x,2x,2x,6x,1x,i6,2x,a1,2x,&
                        &f24.10,4x,a1)') &
                        &'|',ngwf_on_atom,':',&
                        invariant(ingwf,4),'|'
                endif
                ingwf = ingwf + 1           ! original NGWF counting
             enddo
          enddo

          write(stdout,'(a/)') '****************************************&
               &****************************************'
       endif

    endif

    ! Deallocate sparse matrices
    if (.not. nonorthogonal) then
       call sparse_destroy(closetounity)
    endif

    do xyz=1,3
       call sparse_destroy(thespread(xyz))
       if (nonorthogonal) then
          call sparse_destroy(rsr(xyz))
          call sparse_destroy(rs(xyz))
       endif
       call sparse_destroy(r2_elements(xyz))
       call sparse_destroy(r_elements(xyz))
    enddo

    ! ndmh: deallocate temporary arrays
    deallocate(covariant,stat=ierr)
    call utils_dealloc_check('properties_spread','covariant',ierr)
    deallocate(invariant,stat=ierr)
    call utils_dealloc_check('properties_spread','invariant',ierr)

    ! Stop timer
    call timer_clock('properties_spread',2)

  end subroutine properties_spread


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module properties

