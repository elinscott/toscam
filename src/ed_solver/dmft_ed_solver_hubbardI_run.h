   subroutine hub1_run()

      use linalg,              only: minloci
      use genvar,              only: imi, pi
      use globalvar_ed_solver, only: beta_ed, fmos, jjmatrix, &
           slater_coulomb_r, Jhund_Slater_type, uumatrix
      use mesh,                only: build1dmesh
      use matrix,              only: invmat, write_array

      implicit none

      integer                 :: Nc, nw_r, nw_m, nmoments, i, j, k, l, ii, jj
      logical                 :: check
      real(8)                 :: target, target1, target2
      real(8), allocatable    :: mout(:), nout(:)
      complex(8), allocatable :: frequ__(:), Tail(:,:,:), matsubara(:), &
                                 gimp(:,:,:), T1(:,:), T2(:,:)
      real(8), allocatable    :: umn(:,:), ujmn(:,:), Ur(:,:,:,:)
      integer                 :: nmu_adapt
      integer                 :: imu, mu_min
      real(8)                 :: mu_range
      real(8), allocatable    :: mu_cor(:), denshub(:)
      real(8)                 :: mmu_

      Nc       = size(Eimp, 1)/2
      nw_r     = size(sigw, 3)
      nw_m     = size(g_out, 3)
      nmoments = 2
      mmu_     = mmu

      do ii = 1, Nc
         do jj = 1, Nc
            if(ii /= jj .and. abs(UUmatrix(ii, jj)) > 1.d-6)then
               if(.not.Jhund_Slater_type)then
                  UUmatrix(ii, jj) = UUmatrix(ii, jj) + 2.5d0*JJmatrix(ii, jj)
               else
                  UUmatrix(ii, jj) = UUmatrix(ii, jj) + 2.0d0*JJmatrix(ii, jj)
               endif
            endif
         enddo
      enddo

      allocate(matsubara(nw_m), gimp(2*Nc, 2*Nc, nw_m))
      allocate(Tail(nmoments, size(Eimp, 1), size(Eimp, 2)), T1(2*Nc, 2*Nc), &
           T2(2*Nc, 2*Nc))
      allocate(umn(Nc, Nc), ujmn(Nc, Nc), Ur(Nc, Nc, Nc, Nc))
      allocate(mout(Nc), nout(2*Nc))

      umn = 0.
      ujmn = 0.

      if(.not.flag_slater_int)then
         write(*, *) 'HUBBARD I SOLVER IMPLEMENTED FOR SLATER INTERACTION'
         stop
      else
         write(*, *) 'HUBBARD I SOLVER SLATER INTERACTION'
         do i = 1, Nc
            do j = 1, Nc
               do ii = 1, Nc
                  do jj = 1, Nc
                     Ur(i, j, ii, jj) = Slater_Coulomb_r( i, j, jj, ii )
                  enddo
               enddo
            enddo
         enddo
         do i = 1, Nc
            do j = 1, Nc
               umn(i, j)  = Ur(i, j, i, j)
            enddo
         enddo
         do i = 1, Nc
            do j = 1, Nc
               if(i == j) ujmn(i, j) = 0.0
               if(i /= j) ujmn(i, j) = Ur(i, j, i, j) - Ur(i, j, j, i)
            enddo
         enddo
         call write_array(  umn, " U_{ mn} ", unit = 6, short = .true. )
         call write_array( ujmn, " U_{Jmn} ", unit = 6, short = .true. )
      endif

      matsubara(:) = (/( imi*pi/beta_ED*dble(2*i-1), i = 1, nw_m)/)

      if(allocated(frequ__)) deallocate(frequ__)
      allocate(frequ__(nw_r))
      frequ__ = bath%hybridret%freq%vec

      write(*, *) ' FMOS WITH [x] REAL  FREQUENCIES : ', nw_r
      write(*, *) ' FMOS WITH [x] MATSU FREQUENCIES : ', nw_m
      write(*, *) ' INVERSE TEMPERATURE             : ', beta_ED
      write(*, *) ' CHEMICAL POTENTIAL              : ', mmu_

      mu_range = 0.00001d0
      nmu_adapt = 1
      target = 0.0
      target1 = 0.
      target2 = 0.
      inquire(file = 'ed.hub1.target', exist = check)
      if(check)then
         open(unit = 10015, file = 'ed.hub1.target')
         read(10015, *) target1
         read(10015, *) target2
         target = target1 + target2
         read(10015, *) nmu_adapt
         read(10015, *) mu_range
         if(mod(nmu_adapt, 2) == 0) nmu_adapt = nmu_adapt + 1
         write(*, *) ' ======= TARGET IS ======== ', target
         write(*, *) ' ======= NMU_ADAPT ======== ', nmu_adapt
         write(*, *) ' ======= MU_RANGE  ======== ', mu_range
         close(10015)
      endif
      allocate(mu_cor(nmu_adapt), denshub(nmu_adapt))
      if(nmu_adapt > 1)then
         call build1Dmesh(mu_cor, nmu_adapt, -mu_range, mu_range)
      else
         mu_cor = 0.0
      endif

      open(unit = 101010, file = 'filling_mu_hub1')
      do imu = 1, nmu_adapt
         write(*, *) 'MU CORRECTION : ', mu_cor(imu)
         write(*, *) 'MU POTENTIAL  : ', mmu_
         do i = 1, 2*Nc
            Eimp(i, i) = Eimp(i, i)-mmu_-mu_cor(imu)
         enddo
         call gf_HI_fullU( GF = g_out, Tail = Tail, e0f = Eimp, Ur = Ur, umn = &
              umn, ujmn = ujmn, zmsb = matsubara, nlm = Nc, Iwmax = nw_m, nmom &
              = 2, ns = 2, atocc = nout, atmag = mout, temp = 1.d0/beta_ED, &
              verbosity = 0 )
         call sigma_atomic_fullU(GF = g_out, Sigma_mat = self_out, e0f = Eimp, &
              zmsb = matsubara, nlm = Nc, Iwmax = nw_m, ns = 2)
         do i = 1, 2*Nc
            Eimp(i, i) = Eimp(i, i) + mmu_ + mu_cor(imu)
         enddo
         gimp = 0.0
         do i = 1, nw_m
            gimp(:, :, i) = g_out(:, :, i)
            call invmat(2*Nc, gimp(:, :, i))
            gimp(:, :, i) = gimp(:, :, i) - hybrid_in(:, :, i)
            call invmat(2*Nc, gimp(:, :, i))
         enddo
         do i = 1, size(gimp, 1)
            rdens(i) =  2.0*sum( real (  gimp(i, i, 1:nw_m-10) ) )/beta_ED+0.5d0
            !rdens(i) = 2.0*sum( real ( g_out(i, i, 1:nw_m-10) )
            ! )/beta_ED+0.5d0 !BUG TO TEST
         enddo

         denshub(imu) = sum(rdens)
         write(*, '(a, 100f10.3)') ' occ. HUB           : ', sum(nout)
         write(*, '(a, 100f10.3)') ' occ. atom / target : ', denshub(imu), &
              target
         write(101010, *) mmu_ + mu_cor(imu), denshub(imu), sum(nout)
      enddo
      close(101010)

      if(nmu_adapt > 1)then
         mu_min = minloci(abs( denshub - target ))
         mmu_   = mmu_ + mu_cor( mu_min )
         write(*, *) ' =========== DIFF FROM TARGET ============ ', &
              abs( denshub(mu_min) - target )
         do i = 1, 2*Nc
            Eimp(i, i) = Eimp(i, i)-mmu_
         enddo
         call gf_HI_fullU( GF = g_out, Tail = Tail, e0f = Eimp, Ur = Ur, umn = &
              umn, ujmn = ujmn, zmsb = matsubara, nlm = Nc, Iwmax = nw_m, nmom &
              = 2, ns = 2, atocc = nout, atmag = mout, temp = 1.d0/beta_ED, &
              verbosity = 0)
         call sigma_atomic_fullU(GF = g_out, Sigma_mat = self_out, e0f = Eimp, &
              zmsb = matsubara, nlm = Nc, Iwmax = nw_m, ns = 2)
         do i = 1, 2*Nc
            Eimp(i, i) = Eimp(i, i) + mmu_
         enddo
         do i = 1, nw_m
            gimp(:, :, i) = g_out(:, :, i)
            call invmat(2*Nc, gimp(:, :, i))
            gimp(:, :, i) = gimp(:, :, i)-hybrid_in(:, :, i)
            call invmat(2*Nc, gimp(:, :, i))
         enddo
         do i = 1, size(gimp, 1)
            rdens(i) = 2.0*sum( real ( gimp(i, i, 1:nw_m-10)) )/beta_ED+0.5d0
         enddo
      endif

      g_out = gimp

      write(*, *) ' DENSITY HUB1 : ', nout
      write(*, *) ' DENSITY IMP  : ', rdens
      write(*, *) ' DENSITY TOT  : ', sum(rdens)

      do i = 1, 2*Nc
         Eimp(i, i) = Eimp(i, i)-mmu_
      enddo
      call gf_HI_fullU( GF = gw, Tail = Tail, e0f = Eimp, Ur = Ur, umn = umn, &
           ujmn = ujmn, zmsb = frequ__, nlm = Nc, Iwmax = nw_r, nmom = 2, ns = &
           2, atocc = nout, atmag = mout, temp = 1.d0/beta_ED, verbosity = 0)
      call sigma_atomic_fullU(GF = gw, Sigma_mat = sigw, e0f = Eimp, zmsb = &
           frequ__, nlm = Nc, Iwmax = nw_r, ns = 2)
      do i = 1, 2*Nc
         Eimp(i, i) = Eimp(i, i) + mmu_
      enddo

      do i = 1, 2*Nc
         !call plotarray(aimag(matsubara), real(g_out(i, i, :)), 'g_out_re_'
         ! // "_band_" // trim(adjustl(toString(i))) // "_")
         !call plotarray(aimag(matsubara), real(self_out(i, i, :)),
         ! 'self_out_re_' // "_band_" // trim(adjustl(toString(i))) // "_")
         !call plotarray(real(frequ__), real(gw(i, i, :)), 'gw_re_' //
         ! "_band_" // trim(adjustl(toString(i))) // "_")
         !call plotarray(real(frequ__), real(sigw(i, i, :)), 'sigw_re_' //
         ! "_band_" // trim(adjustl(toString(i))) // "_")
      enddo
      do i = 1, 2*Nc
         !call plotarray(aimag(matsubara), aimag(g_out(i, i, :)), 'g_out_im_'
         ! // "_band_" // trim(adjustl(toString(i))) // "_")
         !call plotarray(aimag(matsubara), aimag(self_out(i, i, :)),
         ! 'self_out_im_' // "_band_" // trim(adjustl(toString(i))) // "_")
         !call plotarray(real(frequ__), aimag(gw(i, i, :)), 'gw_im_' //
         ! "_band_" // trim(adjustl(toString(i))) // "_")
         !call plotarray(real(frequ__), aimag(sigw(i, i, :)), 'sigw_im_' //
         ! "_band_" // trim(adjustl(toString(i))) // "_")
      enddo

   end subroutine
