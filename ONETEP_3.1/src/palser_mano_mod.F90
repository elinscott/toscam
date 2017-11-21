! -*- mode: F90 ; mode: font-lock ; column-number-mode: true -*-
!================================================================!
!                                                                !
!         Palser-Manolopoulos kernel optimisation module         !
!                                                                !
! This module implements variants of the method by               !
! Palser and Manolopoulos                                        !
!         Phys. Rev. B. 58(19), 12704 (1998)                     !
! for optimising the density kernel using purification-based     ! 
! algorithms.                                                    !
!----------------------------------------------------------------!
! Written by Chris-Kriton Skylaris on 24/07/2006                 !
!================================================================!

module palser_mano
  implicit none

  private

  public :: palser_mano_kernel_optimise

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine palser_mano_kernel_optimise(k_new, &
       ham, olap, s_inv, n_occ, num_iter)
    !================================================================!
    ! This subroutine implements the "Canonical purification" method !
    ! of Palser and Manolopoulos                                     !
    !         Phys. Rev. B. 58(19), 12704 (1998)                     !
    ! for obtaining the density kernel which optimises the band      !
    ! structure energy for a fixed Hamiltonian. When no truncation   !
    ! is applied to the density kernel, the results obtained are     !
    ! equivalent to the density kernel that would be obtained        !    
    ! by diagonalising the Hamiltonian and building a kernel from    !
    ! the eigenvectors obtained.                                     !
    !----------------------------------------------------------------!
    ! Note: The method is suitable to sparse matrices with truncation!
    !       applied. The iterations are stopped and maximum          !  
    !       convergence is assumed when either of the following      !
    !       conditions is satisfied:                                 !
    !       1) The stable fixed point value c_n is outside [0,1].    !
    !       2) The band structure energy ceases to be monotonically  !
    !          decreasing.                                           !
    !       3) The absolute change in band structure energy per atom !
    !          becomes less than the delta_e_thresh parameter.       !  
    !----------------------------------------------------------------!
    ! Written by Chris-Kriton Skylaris on 24/07/2006                 !
    !================================================================!
    use comms, only: pub_on_root, comms_abort, comms_barrier
    use constants, only: DP, stdout, verbose
    use rundat, only: pub_output_detail
    use services, only: services_flush
    use simulation_cell, only: pub_cell
    use sparse, only: SPAM3, sparse_create, sparse_destroy, sparse_product, &
         sparse_axpy, sparse_scale, sparse_trace, sparse_copy, &
         sparse_extremal_eigenvalue, sparse_num_rows
    use timer, only: timer_clock
    implicit none

    type(SPAM3), intent(inout)    :: k_new   ! optimised kernel
    type(SPAM3), intent(in)       :: ham     ! hamiltonian
    type(SPAM3), intent(in)       :: olap    ! overlap
    type(SPAM3), intent(in)       :: s_inv   ! inverse overlap
    integer, intent(in)           :: n_occ   ! number of up/down electrons
    integer, optional, intent(in) :: num_iter! maximum number of PM iterations 

    ! cks: << local variables>>
    type(SPAM3) :: k_old            ! previous iteration kernel
    type(SPAM3) :: ks               ! current K.S
    type(SPAM3) :: ksk              ! current K.S.K
    type(SPAM3) :: ksks             ! K.S.K.S structure, for calculating efermi
    type(SPAM3) :: sks              ! current S.K.S
    type(SPAM3) :: sinv_ham         ! inverse overlap times hamiltonian
    real(kind=DP) :: min_en         ! minimum hamiltonian eigenvalue
    real(kind=DP) :: max_en         ! maximum hamiltonian eigenvalue
    real(kind=DP) :: average_en     ! everage hamiltonian eigenvalue
    real(kind=DP) :: lambda         ! PM initialisation parameter
    real(kind=DP) :: i_fac          ! PM purification formula factor
    real(kind=DP) :: h_fac          ! PM purification formula factor
    real(kind=DP) :: c_numer        ! numerator of c_n
    real(kind=DP) :: c_denom        ! denominator of c_n
    real(kind=DP) :: c_n            ! unstable fixed point of PM formula
    real(kind=DP) :: scale_fac      ! PM formula scaling factor
    real(kind=DP) :: band_energy    ! current band structure energy
    real(kind=DP) :: old_band_energy! previous band structure energy 
    real(kind=DP),parameter :: delta_e_thresh =1.0E-12_DP ! energy/atom threshold
    real(kind=DP),parameter :: eval_prec = 0.0001_DP ! eigenvalue bound precision
    real(kind=DP) :: lbound         ! minimum kernel eigenvalue
    real(kind=DP) :: ubound         ! maximum kernel eigenvalue
    real(kind=DP) :: mid_occ        ! closest kernel eigenvalue to 0.5
    integer        :: iter          ! iteration counter
    integer        :: loc_num_iter  ! local maximum number of PM iterations 
    logical        :: quitearly     ! loop exit flag
#ifdef DEBUG
    real(kind=DP) :: debug_trace    ! matrix trace for debugging purposes
#endif

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') '  '
    if (pub_on_root) write(stdout,'(a)') &
         'DEBUG: Entering palser_mano_kernel_optimise'
#endif

    ! pdh: nothing to be done if no electrons of this spin
    if (n_occ < 1) then
       k_new%dmtx = 0.0_DP
       return
    end if

    ! cks: start timer
    call timer_clock("palser_mano_kernel_optimise", 1)

    ! cks: initialisations
    call sparse_create(k_old, k_new)
    call sparse_create(ks, k_new, olap)
    call sparse_create(ksk, ks, k_new)
    call sparse_create(sks, olap, ks)
    call sparse_create(sinv_ham, s_inv, ham)
    old_band_energy =huge(1.0_DP)
    quitearly = .false.
    if (present(num_iter)) then
       loc_num_iter = num_iter
    else
       loc_num_iter = 1
    end if


    ! cks: --------- PRINT PRELIMINARY INFO -----------------------
    if ( (pub_output_detail == VERBOSE) .and. pub_on_root) then
       write(stdout,'(a)')' '
       write(stdout,'(a)')'======= &
            &Optimisation of K by Palser-Manolopoulos &
            &canonical purification ======== '
       write(stdout,'(a)')'         Iteration  |  Band structure energy&
            & |    c_n   '
    end if
    ! cks: -----END PRINT PRELIMINARY INFO -----------------------


    ! cks: ======== INITIALISE DENSKERN =====================

    ! cks: S^-1.H
    call sparse_product(sinv_ham, s_inv, ham)

    ! cks: maximum orbital energy
    call sparse_extremal_eigenvalue(sinv_ham, olap, max_en, delta_e_thresh)
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a,f24.14)') &
         'DEBUG: Max Hamiltonian eigenvalue:', max_en
#endif

    ! cks: -S^-1.H
    call sparse_scale(sinv_ham, -1.0_DP)
    ! cks: minimum orbital energy
    call sparse_extremal_eigenvalue(sinv_ham, olap, min_en, delta_e_thresh)
    min_en =-min_en
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a,f24.14)') &
         'DEBUG: Min Hamiltonian eigenvalue:', min_en
#endif

    ! cks: truncate the eigenenergy bounds precision to eval_prec digits as
    ! cks: sparse_extremal_eigenvalue results can become compiler-dependent 
    ! cks: if more digits are used
    max_en =max_en -mod(max_en, eval_prec) + eval_prec
    min_en =min_en -mod(min_en, eval_prec) - eval_prec
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a,f24.14)') &
         'DEBUG: Truncated max Ham eigenvalue:', max_en
#endif
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a,f24.14)') &
         'DEBUG: Truncated min Ham eigenvalue:', min_en
#endif

    ! cks: S^-1.H
    call sparse_scale(sinv_ham, -1.0_DP)

    ! cks: average orbital energy
    average_en = sparse_trace(s_inv, ham) / sparse_num_rows(ham)
#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a,f24.14)') &
         'DEBUG: Average orbital energy:', average_en
#endif

    lambda =min( real(n_occ, kind=DP)/(max_en -average_en), & 
         (sparse_num_rows(ham) -real(n_occ,kind=DP))/(average_en -min_en)  )

    ! cks: (lamda*mu + N_e)/N
    i_fac =(lambda*average_en +n_occ)/sparse_num_rows(ham)

    ! cks: -lambda/N
    h_fac = -lambda/sparse_num_rows(ham)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: lambda:', lambda
    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG:  i_fac:', i_fac
    if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG:  h_fac:', h_fac
#endif

    ! cks: S^-1.H.S^-1
    call sparse_product(k_old, sinv_ham, s_inv )

    ! cks: -lambda/N*S^-1.H.S^-1
    call sparse_scale(k_old, h_fac)

    ! cks: -lambda/N*S^-1.H.S^-1 + (lamda*mu + N_e)/N*S^-1
    call sparse_axpy(k_old, s_inv, i_fac)
    ! cks: ==== END INITIALISE DENSKERN =====================

    !########################################################
    ! cks: ######### PALSER-MANO ITERATIONS #################
    pm_loop: do iter=1,loc_num_iter

       band_energy =sparse_trace(k_old, ham)

#ifdef DEBUG
       if (iter == 1) then
          debug_trace =sparse_trace(s_inv, ham)
          if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[Sinv*H]:',&
               debug_trace 
          debug_trace =sparse_trace(olap, ham)
          if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[S*H]:',&
               debug_trace
          debug_trace =sparse_trace(ham, ham)
          if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[H*H]:',&
               debug_trace
          debug_trace =sparse_trace(olap, olap)
          if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[S*S]:',&
               debug_trace
          debug_trace =sparse_trace(k_old, k_old)
          if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[K*K]:',&
               debug_trace
          debug_trace =sparse_trace(k_old, olap)
          if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[K*S]:',&
               debug_trace
          debug_trace =sparse_trace(s_inv, s_inv)
          if (pub_on_root) write(stdout,'(a,f24.14)') 'DEBUG: tr[Sinv*Sinv]:',&
               debug_trace

       endif
#endif


       ! cks: K.S
       call sparse_product(ks, k_old, olap)

       ! cks: K.S.K
       call sparse_product(ksk, ks, k_old)

       ! cks: S.K.S
       call sparse_product(sks, olap, ks)

       ! cks: -S.K.S
       call sparse_scale(sks, -1.0_DP)

       ! cks: (S-S.K.S)
       call sparse_axpy(sks,olap, 1.0_DP)

       ! cks: tr[K.(S -S.K.S)]
       c_denom =sparse_trace(k_old, sks)
       if ( abs(c_denom) < tiny(1.0_DP)) then
          if (pub_on_root) write(stdout,*)&
               'Error in palser_mano_kernel_optimise: c_denom is zero'
          call comms_barrier
          call comms_abort
       endif

       ! cks: tr[K.S.K.(S -S.K.S)]
       c_numer =sparse_trace(ksk, sks)

       c_n =c_numer/c_denom

       ! cks: print values for current K_old
       if ( (pub_output_detail == VERBOSE) .and. pub_on_root) then
          write(stdout,'(t12,i5,tr5,f22.14,tr3,f10.7)') iter, band_energy, &
               c_n
       end if

       ! cks: Check criteria to quit based on machine precision or 
       ! cks: density kernel truncation and quit if any are satisfied
       if (band_energy > old_band_energy) quitearly =.true.
       if ( (abs(band_energy -old_band_energy)/pub_cell%nat ) < delta_e_thresh)&
            quitearly =.true.
       if (c_n > 1.0_DP) quitearly =.true.
       if (c_n < 0.0_DP) quitearly =.true.

       if (quitearly) then
          ! ndmh: copy k_old to k_new if we quit on first iteration
          if (iter==1) call sparse_copy(k_new,k_old)
          exit pm_loop
       end if


       ! cks: (1+cn)*I - K.S 
       call sparse_scale(ks, -1.0_DP, 1.0_DP+c_n)

       ! cks: (1+cn)*K.S.K - K.S.K.S.K in k_new
       call sparse_product(k_new, ks, ksk)


       ! cks: use different formula to obtain K_new depending on the value of c_n
       if ( c_n >= 0.5_DP) then 
          scale_fac =1.0_DP/c_n
          call sparse_scale(k_new, scale_fac )
       else
          scale_fac = 1.0_DP -2.0_DP*c_n
          call sparse_axpy(k_new, k_old, scale_fac)
          scale_fac = 1.0_DP/(1.0_DP -c_n)
          call sparse_scale(k_new, scale_fac )
       endif

       ! cks: store new density kernel and band energy
       call sparse_copy(k_old, k_new)
       old_band_energy =band_energy

    enddo pm_loop
    ! cks: ###### END PALSER-MANO ITERATIONS #################
    !#########################################################

    ! cks: print warning if convergence criteria not satisfied
    if (pub_on_root .and. .not. quitearly) write(stdout,*) &
         'WARNING: Maximum Palser-Manolopoulos iterations exceeded.'
         
    ! ndmh: Check eigenvalue bounds of ks

    ! ndmh: create temporary matrix for KSKS structures
    call sparse_create(ksks,ks,ks)
    ! ndmh: find upper bound of occupancies
    call sparse_product(ks, k_new, olap)
    call sparse_extremal_eigenvalue(ks,olap,ubound)
    ! ndmh: find lower bound of occupancies
    call sparse_scale(ks,-1.0_DP,1.0_DP)
    call sparse_extremal_eigenvalue(ks,olap,lbound)
    lbound = 1.0_DP - lbound
    ! ndmh: find closest eigenvalue to 0.5
    call sparse_product(ks, k_new, olap)
    call sparse_scale(ks,-1.0_DP,0.5_DP)
    call sparse_product(ksks,ks,ks)
    call sparse_scale(ksks,-1.0_DP,0.0_DP)
    call sparse_extremal_eigenvalue(ksks,olap,mid_occ)
    mid_occ = sqrt(abs(mid_occ)) + 0.5_DP
    
    ! ndmh: warn if eigenvalue suggests degeneracy at eF
    if (pub_on_root.and.(mid_occ>0.15_DP).and.(mid_occ<0.85_DP)) then
       write(stdout,'(a/a)') 'WARNING in palser_mano_kernel_optimise: &
            &Possible degeneracy at fermi level'
       ! ndmh: if the calculation was reckoned to be converged,
       ! ndmh: suggest diagonalisation instead.
       if (quitearly) write(stdout,'(a)') 'Diagonalisation &
            &(maxit_palser_mano=-1) may be required.'
    end if 

    if (pub_on_root.and.(pub_output_detail == VERBOSE)) then
       write(stdout,'(a,f20.14,a,f20.14,a)') &
            'Final occupancy bounds: [',lbound,':',ubound,']'
       write(stdout,'(a,f20.14)') &
            'Closest occupancy to 0.5: ',mid_occ
       write(stdout,'(a)')'===================================&
            &============================================='
    end if
         
    call sparse_destroy(ksks)

    ! cks: deallocate SPAM3 memory
    call sparse_destroy(sinv_ham)
    call sparse_destroy(sks)
    call sparse_destroy(ksk)
    call sparse_destroy(ks)
    call sparse_destroy(k_old)

    ! cks: stop timer
    call timer_clock("palser_mano_kernel_optimise", 2)

#ifdef DEBUG
    if (pub_on_root) write(stdout,'(a)') 'DEBUG: &
         &Leaving palser_mano_kernel_optimise'
#endif

    ! Flush output
    call services_flush

    return

  end subroutine palser_mano_kernel_optimise


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module palser_mano
