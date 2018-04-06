   subroutine reduced_density_matrix_toolbox()

      ! COMPUTE DENSITY MATRIX & ENTANGLEMENT ENTROPY

      use rcmatrix_class,      only: new_rcmatrix, write_raw_rcmatrix
      use globalvar_ed_solver, only: beta_ed, nvecout
      use density_matrix,      only: analyze_density_matrix, compute_density_matrix
      use common_def,          only: dump_message, stats_func
      use overlap_module,      only: compute_overlap

      implicit none

      CALL dump_message(TEXT = "##############################")
      CALL dump_message(TEXT = "### COMPUTE DENSITY MATRIX ###")
      CALL dump_message(TEXT = "##############################")

      CALL new_rcmatrix(dmat, AIM%impurity%nstates)
      write(*, *) ' ==== compute density matrix ==== '
      CALL compute_density_matrix(dmat, AIM, beta_ED, GS)
      CALL write_raw_rcmatrix(dmat, FILEOUT = './ED_out/density_matrix.dat')
      write(*, *) ' ==== analyze density matrix ==== '
      CALL analyze_density_matrix(dmat%rc, AIM%impurity%iorb, NAMBU = SUPER)

      !-------------------------------------------!
      ! COMPUTE WEIGHT OF GS ON REFERENCE VECTORS !
      !-------------------------------------------!

      if(nvecout > 0)then
         write(*, *) ' ==== compute overlap ==== '
         CALL compute_overlap(AIM, beta_ED, GS, vec_list)
      endif

      if(dump_stat) CALL stats_func()

   end subroutine

   subroutine matsubara_super_green()

      use apply_c,               only: apply_cdo, apply_cup, apply_n_cdo, &
           apply_n_cup, cdo_sector, cup_sector
      use green_class_computeab, only: compute_greenab
      use green_class_computeaa, only: compute_greenaa
      use globalvar_ed_solver,   only: beta_ed, flag_ncup, FLAG_GUP_IS_GDN, &
           para_state
      use genvar,                only: matsubara, retarded
      use green_class,           only: copy_green, pad_green
      use matrix,                only: write_array

      implicit none

      integer :: ii, jj, iii, jjj

      write(*, *) ' ==== compute matsubara green ==== '

      !---------!
      ! SPIN UP !
      !---------!

      write(*, *) '--- COMPUTE PARTICLE/HOLE UP SPIN PART ---'
      CALL compute_greenAA(G(1), AIM, beta_ED, GS, Cup_sector, apply_Cup, &
           COMPUTE_DYN = .not.only_dens)
      CALL pad_green(G(2), G(1))

      if(FLAG_NCUP)then
         write(*, *) '--- COMPUTE PARTICLE/HOLE NDN_C_UP SPIN PART ---'
         CALL compute_greenAA(GN(1), AIM, beta_ED, GS, Cup_sector, &
              apply_N_Cup, COMPUTE_DYN = .not.only_dens)
         CALL pad_green(GN(2), GN(1))
      endif

      !-----------!
      ! SPIN DOWN !
      !-----------!

      if((.not.para_state .or. AIM%BATH%SUPER) .and. .not.FLAG_GUP_IS_GDN)then
         write(*, *) '--- COMPUTE DOWN PART AS WELL ---'
         CALL compute_greenAA(G(2), AIM, beta_ED, GS, Cdo_sector, apply_Cdo, &
              COMPUTE_DYN = .not.only_dens)
         if(FLAG_NCUP)then
            write(*, *) 'compute correlated hopping spin down'
            CALL compute_greenAA(GN(2), AIM, beta_ED, GS, Cdo_sector, &
                 apply_N_Cdo, COMPUTE_DYN = .not.only_dens)
         endif
      else
         write(*, *) '---ASSUME PARAMAGNETIC STATE---'
         CALL copy_green(G(2), G(1))
         if(FLAG_NCUP)then
            CALL copy_green(GN(2), GN(1))
         endif
      endif

      if(FLAG_NCUP)then
         write(*, *) 'compute correlated hopping non-hermitian sandwitch'
         CALL compute_greenAB(GNF(1), GNF(2), GN(1), G(1), AIM, beta_ED, GS, &
              Cup_sector, Cup_sector, apply_N_Cup, apply_Cup, COMPUTE_DYN = &
              .not.only_dens)
         write(*, *) 'done COR HOPPINGS'
         write(*, *) 'maxval GNF 1', maxval(abs(GNF(1)%correl(2, 1)%fctn))
         write(*, *) 'maxval GNF 2', maxval(abs(GNF(2)%correl(2, 1)%fctn))
         write(*, *) 'maxval GN  1', maxval(abs(GN(1)%correl(2, 1)%fctn))
         write(*, *) 'maxval GN  2', maxval(abs(GN(2)%correl(2, 1)%fctn))
         write(*, *) 'maxval G   1', maxval(abs(G(1)%correl(2, 1)%fctn))
         write(*, *) 'maxval G   2', maxval(abs(G(2)%correl(2, 1)%fctn))
         CALL copy_green(GNF(2), GNF(1))
         CALL Nambu_green(GNAMBUN, G = GNF)
         if((.not.para_state .or. AIM%BATH%SUPER) .and. .not.FLAG_GUP_IS_GDN) &
              then
            write(*, *) 'C UP cor hopping only for para state and SZ basis'
            stop
         endif
      endif

      if(.not.only_dens)then
         IF(AIM%BATH%SUPER .AND. SUPER)THEN
            write(*, *) 'compute Green Matsubara'
            CALL compute_greenAB(GF(1), GF(2), G(1), G(2), AIM, beta_ED, GS, &
                 Cup_sector, Cdo_sector, apply_Cup, apply_Cdo, COMPUTE_DYN = &
                 .not.only_dens)
            call write_array( GF(1)%correl(1, 1)%fctn(:, :, 4), "GF(1) --", &
                 unit = 98, short = .true.)
            call write_array( GF(1)%correl(2, 2)%fctn(:, :, 4), "GF(2) + + ", &
                 unit = 98, short = .true.)
            CALL Nambu_green(GNAMBU, G = G, GF = GF)
         ELSE
            CALL Nambu_green(GNAMBU, G = G)
         ENDIF
      endif

      if(.not.only_dens) CALL compute_self_energy(SNAMBU, GNAMBU%correl(2, 1), &
           AIM, retarded = .false.)
      write(*, *) '..... done .....'

   end subroutine

   subroutine compute_retarded_functions()

      ! COMPUTE RETARDED GREEN S FUNCTION MATSUBARA + RETARDED BOSONIC
      ! CORRELATIONS

      use apply_c, only: apply_cdo, apply_cup, cdo_sector, cup_sector
      use green_class_computeab, only: compute_greenab
      use green_class_computeaa, only: compute_greenaa
      use globalvar_ed_solver,   only: beta_ed, flag_ncup, FLAG_GUP_IS_GDN, &
           para_state
      use green_class,           only: copy_green, pad_green
      use genvar,                only: retarded

      implicit none

      write(*, *) ' ==== compute retarted green  ==== '

      !---------!
      ! SPIN UP !
      !---------!

      write(*, *) '--- COMPUTE PARTICLE/HOLE UP SPIN PART ---'
      CALL compute_greenAA(Gret(1), AIM, beta_ED, GS, Cup_sector, apply_Cup)
      CALL pad_green(Gret(2), Gret(1))

      !------------!
      ! SPIN DOWN  !
      !------------!

      if((.not.para_state .or. AIM%BATH%SUPER) .and. .not.FLAG_GUP_IS_GDN)then
         write(*, *) '--- COMPUTE DOWN PART AS WELL ---'
         CALL compute_greenAA(Gret(2), AIM, beta_ED, GS, Cdo_sector, apply_Cdo)
      else
         write(*, *) '---ASSUME PARAMAGNETIC STATE---'
         call copy_green(Gret(2), Gret(1))
      endif

      IF(AIM%BATH%SUPER .AND. SUPER)THEN
         write(*, *) 'COMPUTE GF RET'
         CALL compute_greenAB(GFret(1), GFret(2), Gret(1), Gret(2), AIM, &
              beta_ED, GS, Cup_sector, Cdo_sector, apply_Cup, apply_Cdo)
         CALL Nambu_green(GNAMBUret, G = Gret, GF = GFret)
      ELSE
         CALL Nambu_green(GNAMBUret, G = Gret)
      ENDIF

      CALL compute_self_energy(SNAMBUret, GNAMBUret%correl(2, 1), AIM, &
           retarded = .true.)

   end subroutine

   subroutine bosonic_equal_time()

      ! COMPUTE BOSONIC CORRELATIONS: EQUAL-TIME + MATSUBARA IF NECESSARY

      use apply_ns, only: apply_n, apply_s, apply_sz, n_sector, s_sector, &
           sz_sector
      use green_class_computeaa, only: compute_greenaa
      use apply_p,               only: apply_p3, apply_p4, p_sector
      use globalvar_ed_solver,   only: always_compute_static_obs, beta_ed
      use common_def,            only: stats_func

      implicit none

      write(*, *) ' ==== compute Sz ==== '
      if(always_compute_static_obs .or. compute_Sz) CALL compute_greenAA(Sz, &
           AIM, beta_ED, GS, Sz_sector, apply_Sz, COMPUTE_DYN = .false.)
      if(dump_stat) call stats_func("APRES Sz")

      write(*, *) ' ==== compute Spm ==== '
      if(always_compute_static_obs .or. compute_Spm) CALL compute_greenAA(Spm, &
           AIM, beta_ED, GS, S_sector, apply_S, COMPUTE_DYN = .false.)
      if(dump_stat) call stats_func("APRES Spm")

      write(*, *) ' ==== compute N ==== '
      if(always_compute_static_obs .or. compute_N) CALL compute_greenAA(N, &
           AIM, beta_ED, GS, N_sector, apply_N, COMPUTE_DYN = .false.)
      if(dump_stat) call stats_func("APRES N")

      write(*, *) ' ==== P3 ==== '
      if(compute_P3) CALL compute_greenAA(P3, AIM, beta_ED, GS, P_sector, &
           apply_P3, COMPUTE_DYN = .false.)
      if(dump_stat) call stats_func("APRES P3")

      write(*, *) ' ==== P4 ==== '
      if(compute_P4) CALL compute_greenAA(P4, AIM, beta_ED, GS, P_sector, &
           apply_P4, COMPUTE_DYN = .false.)
      if(dump_stat) call stats_func("APRES P4")

   end subroutine

   subroutine bosonic_dynamic()

      ! COMPUTE DYNAMICAL BOSONIC CORRELATIONS IF NECESSARY

      use apply_ns, only: apply_n, apply_s, apply_sz, n_sector, s_sector, &
           sz_sector
      use green_class_computeaa, only: compute_greenaa
      use apply_p,               only: apply_p3, apply_p4, p_sector
      use genvar,                only: iproc, matsubara, quarter
      use globalvar_ed_solver,   only: beta_ed
      use green_class,           only: write_green
      use correl_class,          only: dump_chi0, write_correl

      implicit none




      IF(ANY(Sz%compute))THEN
         write(*, *)  ' ==== compute Sz ==== '
         CALL compute_greenAA(Sz,   AIM, beta_ED, GS, Sz_sector, apply_Sz)
         write(*, *) ' ==== compute Szret ==== '
         CALL compute_greenAA(Szret, AIM, beta_ED, GS, Sz_sector, apply_Sz)
         IF(iproc == 1 .and. write_cor) CALL write_green(Sz)
         IF(iproc == 1 .and. write_cor) CALL write_green(Szret)
         write(*, *) 'dump Chi_0 on matsubara axis'
         call dump_chi0(G(1)%correl(2, 1), beta)
      ENDIF

      IF(ANY(Spm%compute))THEN
         write(*, *) ' ==== compute Spm ==== '
         CALL compute_greenAA(Spm,   AIM, beta_ED, GS, S_sector, apply_S)
         write(*, *) ' ==== compute Spmret ==== '
         CALL compute_greenAA(Spmret, AIM, beta_ED, GS, S_sector, apply_S)
         IF(iproc == 1 .and. write_cor) CALL write_green(Spm)
         IF(iproc == 1 .and. write_cor) CALL write_green(Spmret)
      ENDIF

      IF(ANY(N%compute))THEN
         write(*, *) ' ==== compute N ==== '
         CALL compute_greenAA(N,   AIM, beta_ED, GS, N_sector, apply_N)
         write(*, *) ' ==== compute Nret ==== '
         CALL compute_greenAA(Nret, AIM, beta_ED, GS, N_sector, apply_N)
         IF(iproc == 1 .and. write_cor) CALL write_green(N)
         IF(iproc == 1 .and. write_cor) CALL write_green(Nret)
      ENDIF

      IF(ANY(P3%compute))THEN
         write(*, *) ' ==== compute P3 ==== '
         CALL compute_greenAA(P3,   AIM, beta_ED, GS, P_sector, apply_P3)
         write(*, *) ' ==== compute P3ret ==== '
         CALL compute_greenAA(P3ret, AIM, beta_ED, GS, P_sector, apply_P3)
         IF(iproc == 1 .and. write_cor) CALL write_green(P3)
         IF(iproc == 1 .and. write_cor) CALL write_green(P3ret)
         CHI%fctn = - ( P3%correl(1, 1)%fctn + P3%correl(2, 2)%fctn - &
              P3%correl(1, 2)%fctn - P3%correl(2, 1)%fctn ) * quarter**2
         CHIret%fctn = - ( P3ret%correl(1, 1)%fctn + P3ret%correl(2, 2)%fctn - &
              P3ret%correl(1, 2)%fctn - P3ret%correl(2, 1)%fctn ) * quarter**2
         IF(iproc == 1 .and. write_cor) CALL write_correl(CHI)
         IF(iproc == 1 .and. write_cor) CALL write_correl(CHIret)
      ENDIF

      IF(ANY(P4%compute))THEN
         write(*, *) ' ==== compute P4 ==== '
         CALL compute_greenAA(P4,   AIM, beta_ED, GS, P_sector, apply_P4)
         write(*, *) ' ==== compute P4ret ==== '
         CALL compute_greenAA(P4ret, AIM, beta_ED, GS, P_sector, apply_P4)
         IF(iproc == 1 .and. write_cor) CALL write_green(P4)
         IF(iproc == 1 .and. write_cor) CALL write_green(P4ret)
      ENDIF

   end subroutine

