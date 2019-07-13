   subroutine fmos_solver()

      ! matsubara frequencies : matsubara
      ! real frequ            : frequ
      ! temperature           : beta_ED

      use genvar, only: imi, pi
      use globalvar_ed_solver, only: beta_ed, fmos, fmos_fluc, fmos_iter, &
         fmos_mix, jjmatrix, JHUND_SLATER_TYPE, uumatrix
      use matrix, only: diag, id, invmat, write_array

      implicit none

      integer                 :: Nc, Nb, s1, s2, iter, i, j, k, l, m, kp, mp, &
                                 ip, jp, nw_r, nw_m, Ntot, lp, ii, jj
      real(kind=DP), allocatable    :: Eb(:, :, :), Ei(:, :, :), Ei_(:, :, :)
      complex(kind=DP), allocatable :: Vb(:, :, :), Vbath(:, :, :)
      real(kind=DP), allocatable    :: densav(:, :, :), ueff(:, :), Umma(:, :), &
                                 Ummb(:, :, :), Ummc(:, :, :), UUeff(:, :, :, :)
      complex(kind=DP), allocatable :: SIAM_T(:, :, :)
      complex(kind=DP), allocatable :: I1a(:, :, :), I1b(:, :, :, :), I2a(:, :, :), &
                                 I2b(:, :, :, :), I1c(:, :, :, :), I2c(:, :, :, :)
      complex(kind=DP), allocatable :: Gdyn_m(:, :, :, :), Gdyn_r(:, :, :, :)
      complex(kind=DP), allocatable :: G_ff_m(:, :, :, :), G_ff_r(:, :, :, :), &
                                 S_ff_m_old(:, :, :, :), S_ff_m(:, :, :, :), &
                                 S_ff_r(:, :, :, :)
      complex(kind=DP), allocatable :: Dmasig(:, :, :), Dmbsig(:, :, :, :), &
                                 Dmcsig(:, :, :, :), DmasigTild(:, :, :), &
                                 DmbsigTild(:, :, :, :), DmcsigTild(:, :, :, :)
      real(kind=DP), allocatable    :: matsubara(:), UU(:, :, :, :)
      logical                 :: check_previous
      real(kind=DP)                 :: temp, ttt
      complex(kind=DP), allocatable :: Dmsig(:, :, :), frequ__(:)
      complex(kind=DP), allocatable :: AA(:, :, :), BB(:, :, :, :), CC(:, :, :, :), tempvec(:)
      integer                 :: mstep
      logical, parameter      :: use_def_feng = .true. !BUG
      logical                 :: fmos_fluc_

      Nc = size(Eimp, 1)/2
      Nb = AIM%Nb
      Ntot = Nc + Nb
      Eimp = flip_matrix_Eimp(Eimp, .true.)
      nw_r = size(bath%hybrid%freq%vec(:))
      nw_m = size(g_out, 3)

      if (allocated(frequ__)) deallocate (frequ__)
      allocate (frequ__(nw_r))

      frequ__ = bath%hybridret%freq%vec

      write (*, *) 'FMOS WITH [x] REAL  FREQUENCIES : ', nw_r
      write (*, *) 'FMOS WITH [x] MATSU FREQUENCIES : ', nw_m
      write (*, *) 'TOTAL NUMBER OF SITES           : ', Ntot
      write (*, *) 'INVERSE TEMPERATURE             : ', beta_ED
      write (*, *) 'CHEMICAL POTENTIAL              : ', mmu

      if (allocated(Eb)) deallocate (Eb, Vb, Ei, Ei_, matsubara, UU)
      allocate (Eb(-1:1, Nb, Nb), Vb(-1:1, Nc, Nb), Vbath(-1:1, Nc, Nb/Nc), &
                Ei_(-1:1, Nc, Nc), Ei(-1:1, Nc, Nc), matsubara(nw_m), UU(-1:1, &
                                                                         -1:1, Nc, Nc))
      matsubara(:) = (/(pi/beta_ED*dble(2*i - 1), i=1, nw_m)/)
      write (*, *) 'FIRST MATSUBARA FREQU           : ', matsubara(1)

      if (allocated(G_ff_m)) deallocate (G_ff_r, G_ff_m, S_ff_m, S_ff_r)
      allocate (G_ff_m(-1:1, Nc, Nc, nw_m), G_ff_r(-1:1, Nc, Nc, nw_r), &
                S_ff_m(-1:1, Nc, Nc, nw_m), S_ff_r(-1:1, Nc, Nc, nw_r))
      if (allocated(densav)) deallocate (densav, ueff, Umma, Ummb, Ummc, UUeff, &
                                         Gdyn_m, Gdyn_r)
      allocate (densav(-1:1, Nc + Nb, Nc + Nb), ueff(-1:1, Nc), UUeff(-1:1, &
                                                                      -1:1, Nc, Nc))
      allocate (S_ff_m_old(-1:1, Nc, Nc, nw_m), Gdyn_m(-1:1, Nc + Nb, Nc + Nb, &
                                                       nw_m), Gdyn_r(-1:1, Nc + Nb, Nc + Nb, nw_r))
      allocate (Umma(-1:1, Nc), Ummb(-1:1, Nc, Nc), Ummc(-1:1, Nc, Nc))
      if (allocated(I1a)) deallocate (I1a, I1b, I2a, I2b, I1c, I2c)
      allocate (I1a(-1:1, Nc, nw_m), I1b(-1:1, Nc, Nc, nw_m), I2a(-1:1, Nc, &
                                                                nw_m), I2b(-1:1, Nc, Nc, nw_m), I1c(-1:1, Nc, Nc, nw_m), I2c(-1:1, &
                                                                                                                      Nc, Nc, nw_m))
      allocate (Dmasig(-1:1, Nc, nw_m), Dmbsig(-1:1, Nc, Nc, nw_m), &
                Dmcsig(-1:1, Nc, Nc, nw_m), DmasigTild(-1:1, Nc, nw_m), &
                DmbsigTild(-1:1, Nc, Nc, nw_m), DmcsigTild(-1:1, Nc, Nc, nw_m))
      allocate (Dmsig(-1:1, Nc, nw_m))
      allocate (AA(-1:1, Nc, nw_m), BB(-1:1, Nc, Nc, nw_m), CC(-1:1, Nc, Nc, &
                                                               nw_m), tempvec(nw_m))

      do ii = 1, Nc
         do jj = 1, Nc
            if (ii /= jj .and. abs(UUmatrix(ii, jj)) > 1.d-6) then
               if (.not. Jhund_Slater_type) then
                  UUmatrix(ii, jj) = UUmatrix(ii, jj) + 2.5d0*JJmatrix(ii, jj)
               else
                  UUmatrix(ii, jj) = UUmatrix(ii, jj) + 2.0d0*JJmatrix(ii, jj)
               endif
            endif
         enddo
      enddo

      UU = 0.
      do s1 = -1, 1, 2
         do s2 = -1, 1, 2
            do i = 1, Nc
               do j = 1, Nc
                  if (s1 /= s2 .and. i == j) UU(s1, s2, i, j) = UUmatrix(i, j)
                  if (s1 == s2 .and. i /= j) UU(s1, s2, i, j) = (UUmatrix(i, j) &
                                                                 - 2.0*JJmatrix(i, j))
                  if (s2 /= s1 .and. i /= j) UU(s1, s2, i, j) = (UUmatrix(i, j) &
                                                                 - JJmatrix(i, j))
               enddo
            enddo
         enddo
      enddo

      call write_array(UUmatrix(:, :), 'U matrix', short=.true., unit=6)
      call write_array(JJmatrix(:, :), 'J matrix', short=.true., unit=6)
      call write_array(UU(1, 1, :, :), 'COULOMB + + ', short=.true., unit=6)
      call write_array(UU(1, -1, :, :), 'COULOMB + -', short=.true., unit=6)
      call write_array(UU(-1, 1, :, :), 'COULOMB - + ', short=.true., unit= &
                       6)
      call write_array(UU(-1, -1, :, :), 'COULOMB --', short=.true., unit=6)

      Eb(-1, :, :) = BATH%Eb(2)%rc%mat
      Eb(1, :, :) = BATH%Eb(1)%rc%mat
      !BUG
      Eb(-1, :, :) = Eb(-1, :, :) + Id(Nc)*mmu
      Eb(1, :, :) = Eb(1, :, :) + Id(Nc)*mmu
      !END BUG
      Vb(-1, :, :) = transpose(BATH%Vbc(2)%rc%mat)
      Vb(1, :, :) = transpose(BATH%Vbc(1)%rc%mat)
      Ei_(1, :, :) = Eimp(1:Nc, 1:Nc)
      Ei_(-1, :, :) = Eimp(Nc + 1:2*Nc, Nc + 1:2*Nc)
      S_ff_m = 0.
      S_ff_r = 0.

      mstep = Nb/Nc
      do s1 = -1, 1, 2
         do m = 1, Nc
            k = 0
            do l = (m - 1)*mstep + 1, m*mstep
               k = k + 1
               Vbath(s1, m, k) = Vb(s1, m, l)
            enddo
         enddo
      enddo

      if (allocated(SIAM_T)) deallocate (SIAM_T)
      allocate (SIAM_T(-1:1, Nb + Nc, Nb + Nc))

      write (*, *) 'FMOS SOLVER WITH [x]/[y] BATH/IMP SITES : ', Nb, Nc
      do s1 = -1, 1, 2
         call write_array(Eb(s1, :, :), 'Epsilon matrix', unit=6, ultrashort &
                          =.true.)
         call write_array(real(Vb(s1, :, :)), 'V matrix', unit=6, short= &
                          .true.)
         call write_array(real(Vbath(s1, :, :)), 'V matrix packed', unit=6, &
                          short=.true.)
         call write_array(Ei_(s1, :, :), 'Eimp matrix', unit=6, ultrashort= &
                          .true.)
      enddo

      !------------------!
      do iter = 1, fmos_iter
         !------------------!

         if (iter < fmos_iter/2) then
            fmos_fluc_ = .false.
         else
            fmos_fluc_ = fmos_fluc
         endif

         write (*, *) 'FMOS ITERATION : ', iter

         !----------------!
         ! static offdiag !
         !----------------!

         !BUG
         if (.false.) then
            Vb = 0.10
            Ei_ = 0.0
            Eb = 0.0
            S_ff_m = 0.0
            mmu = 0.0
         endif
         !END BUG

         SIAM_T = 0.d0
         do s1 = -1, 1, 2
            SIAM_T(s1, 1:Nc, Nc + 1:Nc + Nb) = Vb(s1, 1:Nc, 1:Nb)
            SIAM_T(s1, :, :) = SIAM_T(s1, :, :) + transpose(conjg(SIAM_T(s1, &
                                                                         :, :)))
            SIAM_T(s1, 1:Nc, 1:Nc) = Ei_(s1, 1:Nc, 1:Nc)
            SIAM_T(s1, Nc + 1:Nc + Nb, Nc + 1:Nc + Nb) = Eb(s1, 1:Nb, 1:Nb)
            if (maxval(abs(SIAM_T(s1, :, :) - transpose(conjg(SIAM_T(s1, :, &
                                                                     :))))) > 1.d-4) then
               write (*, *) 'ERROR FMOS SIAM_T not HERMITIC'
               stop
            endif
            do i = 1, nw_m
               Gdyn_m(s1, 1:Ntot, 1:Ntot, i) = (matsubara(i)*imi + &
                                                mmu)*Id(Ntot) - SIAM_T(s1, 1:Ntot, 1:Ntot)
               Gdyn_m(s1, 1:Nc, 1:Nc, i) = Gdyn_m(s1, 1:Nc, 1:Nc, i) - &
                                           S_ff_m(s1, 1:Nc, 1:Nc, i)
               call invmat(Ntot, Gdyn_m(s1, :, :, i))
            enddo
         enddo

         if (iter == 1) then
            do s1 = -1, 1, 2
               do m = 1, Nc
                  !call plotarray(matsubara, real(Gdyn_m(s1, m, m, :)),
                  ! 'FMOS_GREEN0_RE_ITER_' // trim(adjustl(toString(iter))) //
                  ! "_" // trim(adjustl(toString(m))) // "_")
                  !call plotarray(matsubara, aimag(Gdyn_m(s1, m, m, :)),
                  ! 'FMOS_GREEN0_IM_ITER_' // trim(adjustl(toString(iter))) //
                  ! "_" // trim(adjustl(toString(m))) // "_")
               enddo
            enddo
         endif

         if (iter == 1) then
            inquire (file='ed.fmos', exist=check_previous)
            if (check_previous) then
               write (*, *) 'READING FROM G FROM EARLIER CALCULATION'
               open (unit=67998, file='ed.fmos', form='unformatted')
               read (67998) Gdyn_m
               close (67998)
            endif
         endif

         !---------------!
         ! static values !
         !---------------!

         densav = 0.0
         do s1 = -1, 1, 2
            do i = 1, Nc + Nb
               do j = 1, Nc + Nb
                  if (.false.) then
                     densav(s1, i, j) = 1.d0/beta_ED*sum(real(Gdyn_m(s1, i, j, &
                                                                     :)))*2.d0
                     if (i == j) densav(s1, i, j) = densav(s1, i, j) + 0.5d0
                  else
                     call Fourier(nw_m, Gdyn_m(s1, i, j, 1:nw_m), &
                                  matsubara(1:nw_m), densav(s1, i, j), beta_ED, &
                                  tailmax=3.d0, green=i == j)
                  endif
               enddo
               if (rank == 0 .and. s1 == 1) write (*, *) 'FMOS DENSITY :', i, &
                  s1, densav(s1, i, i)
               if (rank == 0 .and. i <= Nc .and. s1 == 1) write (*, *) 'MAX &
                    &HOPPING VAL : ', maxval(abs(densav(s1, i, Nc + 1:Ntot)))
               if (rank == 0 .and. i <= Nc .and. s1 == 1) write (*, *) 'MAX &
                    &HOPPING LOC : ', maxloc(abs(densav(s1, i, Nc + 1:Ntot)))

            enddo
         enddo
         write (*, *) 'TOTAL DENSITY : ', sum(diag(densav(-1, 1:Nc, 1:Nc)) + &
                                              diag(densav(+1, 1:Nc, 1:Nc)))

         if (.not. fmos_fluc_) then
            Ei = 0.0
            do s1 = -1, 1, 2
               do m = 1, Nc
                  temp = 0.d0
                  do l = 1, Nc
                     if (m /= l) then
                        temp = temp + UU(s1, s1, l, m)*densav(s1, l, &
                                                              l)*densav(s1, m, m) + UU(s1, -s1, l, &
                                                                                       m)*densav(-s1, l, l)*densav(s1, m, m)
                     endif
                  enddo
                  Ei(s1, m, m) = Ei_(s1, m, m) + (1.d0 - densav(-s1, m, m))*temp
               enddo
            enddo
         else
            do s1 = -1, 1, 2
               Ei(s1, :, :) = Ei_(s1, :, :)
            enddo
         endif

         !---------------------------!
         ! Ueff, UUeff, Umma, Ummb, Ummc !
         !---------------------------!

         ueff = 0.
         do s1 = -1, 1, 2
            do m = 1, Nc
               ueff(s1, m) = UU(1, -1, m, m)
               do l = 1, Nc
                  if (l /= m) ueff(s1, m) = ueff(s1, m) + UU(s1, s1, l, &
                                                             m)*densav(s1, l, l) + UU(s1, -s1, l, m)*densav(-s1, l, &
                                                                                                            l)
               enddo
            enddo
         enddo

         UUeff = 0.0
         do s1 = -1, 1, 2
            do l = 1, Nc
               do m = 1, Nc
                  if (l /= m) then
                     UUeff(s1, s1, l, m) = UU(s1, s1, l, m) + UU(1, -1, m, &
                                                                 m)*densav(-s1, m, m) + UU(s1, -s1, l, m)*densav(-s1, &
                                                                                                                 l, l)
                     do lp = 1, Nc
                        if (lp /= m .and. lp /= l) then
                           UUeff(s1, s1, l, m) = UUeff(s1, s1, l, m) + UU(s1, &
                                                                          s1, lp, m)*densav(s1, lp, lp) + UU(s1, -s1, &
                                                                                                          lp, m)*densav(-s1, lp, lp)
                        endif
                     enddo
                     UUeff(s1, -s1, l, m) = UU(s1, -s1, l, m) + UU(1, -1, m, &
                                                                   m)*densav(-s1, m, m) + UU(s1, s1, l, m)*densav(s1, &
                                                                                                                  l, l)
                     do lp = 1, Nc
                        if (lp /= m .and. lp /= l) then
                           UUeff(s1, -s1, l, m) = UUeff(s1, -s1, l, m) + &
                                                  UU(s1, s1, lp, m)*densav(s1, lp, lp) + UU(s1, &
                                                                                            -s1, lp, m)*densav(-s1, lp, lp)
                        endif
                     enddo
                  endif
               enddo
            enddo
         enddo

         write (*, *) 'Diag ueff + 1 : ', ueff(1, :)
         write (*, *) 'Diag ueff -1 : ', ueff(-1, :)
         call write_array(UUeff(1, 1, :, :), 'Ueff + + ', unit=6, short= &
                          .true.)
         call write_array(UUeff(1, -1, :, :), 'Ueff + -', unit=6, short= &
                          .true.)

         if (.not. fmos_fluc_) then
            Umma = 0.
            do s1 = -1, 1, 2
               do m = 1, Nc
                  Umma(s1, m) = UU(1, -1, m, m)
                  do l = 1, Nc
                     Umma(s1, m) = Umma(s1, m) + 2.d0*(UU(s1, s1, l, &
                                                          m)*densav(s1, l, l) + UU(s1, -s1, l, m)*densav(-s1, &
                                                                                                         l, l))
                  enddo
               enddo
            enddo
         else
            Umma = 0.
            do s1 = -1, 1, 2
               do m = 1, Nc
                  Umma(s1, m) = ueff(s1, m)
                  do l = 1, Nc
                     if (l /= m) then
                        Umma(s1, m) = Umma(s1, m) + 2.d0*(UU(s1, s1, l, &
                                                             m)*densav(s1, l, l) + UU(s1, -s1, l, &
                                                                                      m)*densav(-s1, l, l) + 2.d0*UU(1, -1, l, &
                                                                                              l)*densav(s1, l, l)*densav(-s1, l, l))
                     endif
                  enddo
               enddo
            enddo
         endif

         Ummb = 0.
         do s1 = -1, 1, 2
            do l = 1, Nc
               do m = 1, Nc
                  if (l /= m) then
                     Ummb(s1, l, m) = UUeff(s1, s1, l, m)
                     Ummb(s1, l, m) = Ummb(s1, l, m) + 2.d0*(UU(1, -1, m, &
                                                                m)*densav(-s1, m, m) + UU(s1, -s1, l, m)*densav(-s1, &
                                                                                             l, l) + UU(1, -1, m, m)*densav(s1, m, &
                                                                              m)*densav(-s1, m, m) + UU(1, -1, l, l)*densav(s1, l, &
                                                                                                               l)*densav(-s1, l, l))
                     do lp = 1, Nc
                        if (lp /= l .and. lp /= m) then
                           Ummb(s1, l, m) = Ummb(s1, l, m) + 2.d0*(UU(s1, s1, &
                                                                      lp, m)*densav(s1, lp, lp) + UU(s1, -s1, lp, &
                                                                                       m)*densav(-s1, lp, lp) + 2.d0*UU(1, -1, lp, &
                                                                                         lp)*densav(s1, lp, lp)*densav(-s1, lp, lp))
                        endif
                     enddo
                  endif
               enddo
            enddo
         enddo

         Ummc = 0.
         do s1 = -1, 1, 2
            do l = 1, Nc
               do m = 1, Nc
                  if (l /= m) then
                     Ummc(s1, l, m) = UUeff(s1, -s1, l, m)
                     Ummc(s1, l, m) = Ummc(s1, l, m) + 2.d0*(UU(1, -1, m, &
                                                                m)*densav(-s1, m, m) + UU(s1, s1, l, m)*densav(s1, &
                                                                              l, l) + UU(1, -1, m, m)*densav(s1, m, m)*densav(-s1, &
                                                                              m, m) + UU(1, -1, l, l)*densav(s1, l, l)*densav(-s1, &
                                                                                                                              l, l))
                     do lp = 1, Nc
                        if (lp /= l .and. lp /= m) then
                           Ummc(s1, l, m) = Ummc(s1, l, m) + 2.d0*(UU(s1, s1, &
                                                                      lp, m)*densav(s1, lp, lp) + UU(s1, -s1, lp, &
                                                                                       m)*densav(-s1, lp, lp) + 2.d0*UU(1, -1, lp, &
                                                                                         lp)*densav(s1, lp, lp)*densav(-s1, lp, lp))
                        endif
                     enddo
                  endif
               enddo
            enddo
         enddo

         !-------------------------------------------------------------!
         ! Dmsig, Dmasig, Dmbsig, Dmcsig, DmasigTild, DmbsigTild, DmcsigTild !
         !-------------------------------------------------------------!

         do s1 = -1, 1, 2
            do m = 1, Nc
               Dmsig(s1, m, :) = 0.
               do k = 1, Nb
                  Dmsig(s1, m, :) = Dmsig(s1, m, :) + conjg(Vb(s1, m, k))* &
                                    Vb(s1, m, k)/(matsubara(:)*imi + mmu - Eb(s1, k, k))
               enddo
            enddo
         enddo

         if (iter == 1) then
            do s1 = -1, 1, 2
               do m = 1, Nc
                  !call plotarray(matsubara, real(1.d0/(matsubara(:)*imi +
                  ! mmu-Ei_(s1, m, m)-Dmsig(s1, m, :))), 'FMOS_WEISS_RE_ITER_'
                  ! // trim(adjustl(toString(iter))) // "_" //
                  ! trim(adjustl(toString(m))) // "_")
                  !call plotarray(matsubara, aimag(1.d0/(matsubara(:)*imi +
                  ! mmu-Ei_(s1, m, m)-Dmsig(s1, m, :))), 'FMOS_WEISS_IM_ITER_'
                  ! // trim(adjustl(toString(iter))) // "_" //
                  ! trim(adjustl(toString(m))) // "_")
                  !call plotarray(matsubara, real(Dmsig(s1, m, :)),
                  ! 'FMOS_HYBRID_RE_ITER_' // trim(adjustl(toString(iter))) //
                  ! "_" // trim(adjustl(toString(m))) // "_")
                  !call plotarray(matsubara, aimag(Dmsig(s1, m, :)),
                  ! 'FMOS_HYBRID_IM_ITER_' // trim(adjustl(toString(iter))) //
                  ! "_" // trim(adjustl(toString(m))) // "_")
               enddo
            enddo
         endif

         do s1 = -1, 1, 2
            do m = 1, Nc
               Dmasig(s1, m, :) = 0.
               do k = 1, Nb
                  Dmasig(s1, m, :) = Dmasig(s1, m, :) + conjg(Vb(-s1, m, k))* &
                                     Vb(-s1, m, k)/(matsubara(:)*imi + mmu + Ei(-s1, m, &
                                                                                m) - Ei(s1, m, m) - Eb(-s1, k, k))
               enddo
            enddo
         enddo

         do s1 = -1, 1, 2
            do m = 1, Nc
               do l = 1, Nc
                  Dmbsig(s1, l, m, :) = 0.
                  do k = 1, Nb
                     Dmbsig(s1, l, m, :) = Dmbsig(s1, l, m, :) + conjg(Vb(s1, &
                                                                          l, k))*Vb(s1, l, k)/(matsubara(:)*imi + mmu + &
                                                                                         Ei(s1, l, l) - Ei(s1, m, m) - Eb(s1, k, k))
                  enddo
               enddo
            enddo
         enddo

         do s1 = -1, 1, 2
            do m = 1, Nc
               do l = 1, Nc
                  Dmcsig(s1, l, m, :) = 0.
                  do k = 1, Nb
                     Dmcsig(s1, l, m, :) = Dmcsig(s1, l, m, :) + conjg(Vb(-s1, &
                                                                          l, k))*Vb(-s1, l, k)/(matsubara(:)*imi + mmu + &
                                                                                       Ei(-s1, l, l) - Ei(s1, m, m) - Eb(-s1, k, k))
                  enddo
               enddo
            enddo
         enddo

         do s1 = -1, 1, 2
            do m = 1, Nc
               DmasigTild(s1, m, :) = 0.
               do k = 1, Nb
                  DmasigTild(s1, m, :) = DmasigTild(s1, m, :) + conjg(Vb(-s1, &
                                                                         m, k))*Vb(-s1, m, k)/(matsubara(:)*imi + mmu + &
                                                                           Eb(-s1, k, k) - Ei(-s1, m, m) - Ei(s1, m, m) - Umma(s1, &
                                                                                                                                 m))
               enddo
            enddo
         enddo

         do s1 = -1, 1, 2
            do m = 1, Nc
               do l = 1, Nc
                  DmbsigTild(s1, l, m, :) = 0.
                  do k = 1, Nb
                     DmbsigTild(s1, l, m, :) = DmbsigTild(s1, l, m, :) + &
                                               conjg(Vb(s1, l, k))*Vb(s1, l, k)/( &
                                               matsubara(:)*imi + mmu + Eb(s1, k, k) - Ei(s1, l, l) &
                                               - Ei(s1, m, m) - Ummb(s1, l, m))
                  enddo
               enddo
            enddo
         enddo

         do s1 = -1, 1, 2
            do m = 1, Nc
               do l = 1, Nc
                  DmcsigTild(s1, l, m, :) = 0.
                  do k = 1, Nb
                     DmcsigTild(s1, l, m, :) = DmcsigTild(s1, l, m, :) + &
                                               conjg(Vb(-s1, l, k))*Vb(-s1, l, k)/( &
                                               matsubara(:)*imi + mmu + Eb(-s1, k, k) - Ei(-s1, l, &
                                                                                           l) - Ei(s1, m, m) - Ummc(s1, l, m))
                  enddo
               enddo
            enddo
         enddo

         !---------------------------!
         ! I1a, I2a, I1b, I2b, I1c, I2c   !
         !---------------------------!

         do s1 = -1, 1, 2
            do m = 1, Nc
               I1a(s1, m, :) = 0.
               do k = 1, Nb
                  I1a(s1, m, :) = I1a(s1, m, :) + conjg(Vb(-s1, m, k))* &
                                  densav(-s1, m, Nc + k)/(matsubara(:)*imi + mmu + &
                                                          Ei(-s1, m, m) - Ei(s1, m, m) - Eb(-s1, k, k))
                  I1a(s1, m, :) = I1a(s1, m, :) - Vb(-s1, m, k)*densav(-s1, &
                                                                       Nc + k, m)/(matsubara(:)*imi + mmu + Eb(-s1, k, &
                                                                                    k) - Ei(-s1, m, m) - Ei(s1, m, m) - Umma(s1, m))
               enddo
            enddo
         enddo

         do s1 = -1, 1, 2
            do m = 1, Nc
               I2a(s1, m, :) = 0.
               do k = 1, Nb
                  do kp = 1, Nb
                     I2a(s1, m, :) = I2a(s1, m, :) - conjg(Vb(-s1, m, k))* &
                                     Vb(-s1, m, k)*densav(-s1, Nc + kp, Nc + k)/( &
                                     matsubara(:)*imi + mmu + Ei(-s1, m, m) - Ei(s1, m, &
                                                                                 m) - Eb(-s1, k, k))
                     I2a(s1, m, :) = I2a(s1, m, :) - conjg(Vb(-s1, m, k))* &
                                     Vb(-s1, m, k)*densav(-s1, Nc + k, Nc + kp)/( &
                                     matsubara(:)*imi + mmu + Eb(-s1, k, k) - Ei(-s1, m, &
                                                                                 m) - Ei(s1, m, m) - Umma(s1, m))
                  enddo
               enddo
            enddo
         enddo

         do s1 = -1, 1, 2
            do m = 1, Nc
               do l = 1, Nc
                  I1b(s1, m, l, :) = 0.
                  do k = 1, Nb
                     I1b(s1, m, l, :) = I1b(s1, m, l, :) + conjg(Vb(s1, l, k)) &
                                        *densav(s1, l, Nc + k)/(matsubara(:)*imi + mmu + &
                                                                Ei(s1, l, l) - Ei(s1, m, m) - Eb(s1, k, k))
                     I1b(s1, m, l, :) = I1b(s1, m, l, :) - Vb(s1, l, k)* &
                                        densav(s1, Nc + k, l)/(matsubara(:)*imi + mmu + &
                                                               Eb(s1, k, k) - Ei(s1, l, l) - Ei(s1, m, m) - &
                                                               Ummb(s1, l, m))
                  enddo
               enddo
            enddo
         enddo

         do s1 = -1, 1, 2
            do m = 1, Nc
               do l = 1, Nc
                  I2b(s1, m, l, :) = 0.
                  do k = 1, Nb
                     do kp = 1, Nb
                        I2b(s1, m, l, :) = I2b(s1, m, l, :) - conjg(Vb(s1, l, &
                                                                       k))*Vb(s1, l, k)*densav(s1, Nc + kp, Nc + k) &
                                           /(matsubara(:)*imi + mmu + Ei(s1, l, l) - Ei(s1, &
                                                                                        m, m) - Eb(s1, k, k))
                        I2b(s1, m, l, :) = I2b(s1, m, l, :) - conjg(Vb(s1, l, &
                                                                       k))*Vb(s1, l, k)*densav(s1, Nc + k, Nc + kp) &
                                           /(matsubara(:)*imi + mmu + Eb(s1, k, k) - Ei(s1, &
                                                                                        m, m) - Ei(s1, l, l) - Ummb(s1, l, m))
                     enddo
                  enddo
               enddo
            enddo
         enddo

         do s1 = -1, 1, 2
            do m = 1, Nc
               do l = 1, Nc
                  I1c(s1, m, l, :) = 0.

                  do k = 1, Nb
                     I1c(s1, m, l, :) = I1c(s1, m, l, :) + conjg(Vb(-s1, l, &
                                                                    k))*densav(-s1, l, Nc + k)/(matsubara(:)*imi + &
                                                                                 mmu + Ei(-s1, l, l) - Ei(s1, m, m) - Eb(-s1, k, k))
                     if (.not. use_def_feng) then
                        I1c(s1, m, l, :) = I1c(s1, m, l, :) - Vb(-s1, l, k)* &
                                           densav(-s1, Nc + k, l)/(matsubara(:)*imi + mmu &
                                                                   + Eb(-s1, k, k) - Ei(-s1, l, l) - Ei(s1, m, &
                                                                                                        m) - Ummc(s1, l, m))
                     endif
                  enddo
                  if (use_def_feng) then
                     do k = 1, Nb/Nc
                        i = (m - 1)*mstep
                        j = (l - 1)*mstep
                        I1c(s1, m, l, :) = I1c(s1, m, l, :) - Vbath(-s1, m, k) &
                                           *densav(-s1, Nc + k + j, l)/(matsubara(:)*imi &
                                                                        + mmu + Eb(-s1, k + j, k + j) - Ei(-s1, l, &
                                                                                                 l) - Ei(s1, m, m) - Ummc(s1, l, m))
                     enddo
                  endif

               enddo
            enddo
         enddo

         do s1 = -1, 1, 2
            do m = 1, Nc
               do l = 1, Nc
                  I2c(s1, m, l, :) = 0.
                  do k = 1, Nb
                     do kp = 1, Nb
                        I2c(s1, m, l, :) = I2c(s1, m, l, :) - conjg(Vb(-s1, l, &
                                                                       k))*Vb(-s1, l, k)*densav(-s1, Nc + kp, Nc + &
                                                                                           k)/(matsubara(:)*imi + mmu + Ei(-s1, l, &
                                                                                                  l) - Ei(s1, m, m) - Eb(-s1, k, k))
                        I2c(s1, m, l, :) = I2c(s1, m, l, :) - conjg(Vb(-s1, l, &
                                                                       k))*Vb(-s1, l, k)*densav(-s1, Nc + k, Nc + &
                                                                                          kp)/(matsubara(:)*imi + mmu + Eb(-s1, k, &
                                                                                 k) - Ei(-s1, l, l) - Ei(s1, m, m) - Ummc(s1, l, m))
                     enddo
                  enddo
               enddo
            enddo
         enddo

         !----------!
         ! AA, BB, CC !
         !----------!

         if (fmos_fluc_) then
            do s1 = -1, 1, 2
               do m = 1, Nc
                  AA(s1, m, :) = 0.
                  do l = 1, Nc
                     if (l /= m) then
                        AA(s1, m, :) = AA(s1, m, :) + UU(s1, s1, l, &
                                                         m)*densav(s1, l, l) + UU(s1, -s1, l, &
                                                                                  m)*densav(-s1, l, l) + 2.0*UU(1, -1, m, &
                                                                                               m)*densav(s1, l, l)*densav(-s1, l, l)
                        AA(s1, m, :) = AA(s1, m, :) + densav(-s1, m, m)*( &
                                       Dmsig(-s1, l, :) + Dmsig(s1, l, :))
                     endif
                  enddo
                  AA(s1, m, :) = -AA(s1, m, :) + matsubara(:)*imi + mmu - Ei(s1, &
                                                                             m, m) - ueff(s1, m) - Dmsig(s1, m, :) - Dmasig(s1, m, &
                                                                                                           :) - DmasigTild(s1, m, :)
                  AA(s1, m, :) = ueff(s1, m)/AA(s1, m, :)
               enddo
            enddo
         else
            do s1 = -1, 1, 2
               do m = 1, Nc
                  AA(s1, m, :) = 0.
                  AA(s1, m, :) = matsubara(:)*imi + mmu - Ei(s1, m, m) - ueff(s1, &
                                                                       m) - Dmsig(s1, m, :) - Dmasig(s1, m, :) - DmasigTild(s1, m, &
                                                                                                                                  :)
                  AA(s1, m, :) = ueff(s1, m)/AA(s1, m, :)
               enddo
            enddo
         endif

         do s1 = -1, 1, 2
            do m = 1, Nc
               do l = 1, Nc
                  BB(s1, m, l, :) = 0.
                  if (l /= m) then
                     do lp = 1, Nc
                        if (lp /= l) then
                           BB(s1, m, l, :) = BB(s1, m, l, :) + (UU(s1, s1, &
                                                                   lp, m)*densav(s1, lp, lp) + UU(s1, -s1, lp, &
                                                                                        m)*densav(-s1, lp, lp) + 2.0*UU(1, -1, lp, &
                                                                                         lp)*densav(s1, lp, lp)*densav(-s1, lp, lp))
                           BB(s1, m, l, :) = BB(s1, m, l, :) + (densav(-s1, &
                                                                       lp, lp) + densav(s1, lp, lp))*(Dmsig(s1, m, &
                                                                                                            :))
                        endif
                     enddo
                     BB(s1, m, l, :) = -BB(s1, m, l, :) + matsubara(:)*imi + &
                                       mmu - Ei(s1, m, m) - UUeff(s1, s1, l, m) - (UU(1, &
                                                                                      -1, m, m)*densav(-s1, m, m) + UU(-s1, s1, l, &
                                                                              m)*densav(-s1, l, l) + UU(1, -1, m, m)*densav(s1, m, &
                                                                              m)*densav(-s1, m, m) + UU(1, -1, l, l)*densav(s1, l, &
                                                                             l)*densav(-s1, l, l)) - (Dmsig(s1, m, :) - Dmbsig(s1, &
                                                                             l, m, :) - DmbsigTild(s1, l, m, :)) - (densav(-s1, m, &
                                                                                           m) + densav(-s1, l, l))*(Dmsig(s1, m, :))
                     BB(s1, m, l, :) = UUeff(s1, s1, l, m)/BB(s1, m, l, :)
                  endif
               enddo
            enddo
         enddo

         do s1 = -1, 1, 2
            do m = 1, Nc
               do l = 1, Nc
                  CC(s1, m, l, :) = 0.
                  if (l /= m) then
                     do lp = 1, Nc
                        if (lp /= l) then
                           CC(s1, m, l, :) = CC(s1, m, l, :) + (UU(s1, s1, &
                                                                   lp, m)*densav(s1, lp, lp) + UU(s1, -s1, lp, &
                                                                                        m)*densav(-s1, lp, lp) + 2.0*UU(1, -1, lp, &
                                                                                         lp)*densav(s1, lp, lp)*densav(-s1, lp, lp))
                           CC(s1, m, l, :) = CC(s1, m, l, :) + (densav(-s1, &
                                                                       lp, lp) + densav(s1, lp, lp))*(Dmsig(s1, m, &
                                                                                                            :))
                        endif
                     enddo
                     CC(s1, m, l, :) = -CC(s1, m, l, :) + matsubara(:)*imi + &
                                       mmu - Ei(s1, m, m) - UUeff(-s1, s1, l, m) - (UU(1, -1, &
                                                                                       m, m)*densav(-s1, m, m) + UU(s1, s1, l, &
                                                                               m)*densav(s1, l, l) + UU(1, -1, m, m)*densav(s1, m, &
                                                                              m)*densav(-s1, m, m) + UU(1, -1, l, l)*densav(s1, l, &
                                                                             l)*densav(-s1, l, l)) - (Dmsig(s1, m, :) - Dmcsig(s1, &
                                                                             l, m, :) - DmcsigTild(s1, l, m, :)) - (densav(-s1, m, &
                                                                                            m) + densav(s1, l, l))*(Dmsig(s1, m, :))
                     CC(s1, m, l, :) = UUeff(-s1, s1, l, m)/CC(s1, m, l, :)
                  endif
               enddo
            enddo
         enddo

         if (.not. fmos_fluc_) then
            BB = 0.d0
            CC = 0.d0
         endif

         do s1 = -1, 1, 2
            Gdyn_m(s1, 1:Nc, 1:Nc, :) = 0.0
            do m = 1, Nc
               do l = 1, Nc
                  if (l /= m) then
                     Gdyn_m(s1, m, m, :) = Gdyn_m(s1, m, m, :) + BB(s1, m, l, &
                                                                    :)*(Dmsig(s1, m, :)*I1b(s1, m, l, :) + I2b(s1, m, l, &
                                                                              :)) + CC(s1, m, l, :)*(Dmsig(s1, m, :)*I1c(s1, m, l, &
                                                                                                              :) + I2c(s1, m, l, :))
                  endif
               enddo
               Gdyn_m(s1, m, m, :) = -Gdyn_m(s1, m, m, :) + matsubara(:)*imi &
                                     + mmu - Ei(s1, m, m) - Dmsig(s1, m, :) - AA(s1, m, &
                                                                                 :)*(Dmsig(s1, m, :)*I1a(s1, m, :) + I2a(s1, m, :))

               tempvec = 0.0
               do l = 1, Nc
                  if (l /= m) then
                     tempvec(:) = tempvec(:) + BB(s1, m, l, :)*(densav(s1, &
                                                                       l, l) + I1b(s1, m, l, :)) + CC(s1, m, l, :)* &
                                  (densav(-s1, l, l) + I1c(s1, m, l, :))
                  endif
               enddo

               if (s1 == -1) then
                  !call plotarray( matsubara, real(BB(s1, m, 1, :)),
                  ! 'FMOS_BB_RE_' // trim(adjustl(toString(iter))) // "_" //
                  ! trim(adjustl(toString(m))) // "_")
                  !call plotarray( matsubara, aimag(BB(s1, m, 1, :)),
                  ! 'FMOS_BB_IM_' // trim(adjustl(toString(iter))) // "_" //
                  ! trim(adjustl(toString(m))) // "_")
                  !call plotarray( matsubara, real(tempvec),
                  ! 'FMOS_TEMPVEC_RE_' // trim(adjustl(toString(iter))) // "_"
                  ! // trim(adjustl(toString(m))) // "_")
                  !call plotarray( matsubara, aimag(tempvec),
                  ! 'FMOS_TEMPVEC_IM_' // trim(adjustl(toString(iter))) // "_"
                  ! // trim(adjustl(toString(m))) // "_")
               endif
               Gdyn_m(s1, m, m, :) = (1.d0/Gdyn_m(s1, m, m, :))*(1.d0 + &
                                                                 AA(s1, m, :)*(densav(-s1, m, m) + I1a(s1, m, :)) + &
                                                                 tempvec(:))
            enddo
         enddo

         S_ff_m_old = S_ff_m
         do s1 = -1, 1, 2
            S_ff_m(s1, :, :, :) = 0.0
            do m = 1, Nc
               !call plotarray(matsubara, real(Gdyn_m(s1, m, m, :)),
               ! 'FMOS_GREEN_RE_ITER_' // trim(adjustl(toString(iter))) // "_"
               ! // trim(adjustl(toString(m))) // "_")
               !call plotarray(matsubara, aimag(Gdyn_m(s1, m, m, :)),
               ! 'FMOS_GREEN_IM_ITER_' // trim(adjustl(toString(iter))) // "_"
               ! // trim(adjustl(toString(m))) // "_")
               S_ff_m(s1, m, m, :) = (matsubara(:)*imi + mmu) - Ei_(s1, m, m) &
                                     - Dmsig(s1, m, :) - 1.d0/Gdyn_m(s1, m, m, :)
               !call plotarray(matsubara, real(S_ff_m(s1, m, m, :)),
               ! 'FMOS_SELF_RE_ITER_' // trim(adjustl(toString(iter))) // "_"
               ! // trim(adjustl(toString(m))) // "_")
               !call plotarray(matsubara, aimag(S_ff_m(s1, m, m, :)),
               ! 'FMOS_SELF_IM_ITER_' // trim(adjustl(toString(iter))) // "_"
               ! // trim(adjustl(toString(m))) // "_")
               do i = 1, nw_m
                  if (aimag(S_ff_m(s1, m, m, i)) > 0.) S_ff_m(s1, m, m, i) = &
                     CMPLX(real(S_ff_m(s1, m, m, i)), 0.d0)
               enddo
            enddo
         enddo

         if (iter < fmos_iter/2) then
            ttt = min(1.d0, 2*fmos_mix)
         else
            ttt = fmos_mix
         endif

         S_ff_m = ttt*S_ff_m + (1.0 - ttt)*S_ff_m_old

         !------------------!
      end do
      !------------------!

      write (*, *) 'FMOS DONE, COLLECT RESULTS'

      rdens(1:Nc) = diag(densav(+1, :, :))
      rdens(Nc + 1:2*Nc) = diag(densav(-1, :, :))

      self_out = 0.
      g_out = 0.
      do i = 1, nw_m
         do m = 1, Nc
            self_out(m, m, i) = S_ff_m(1, m, m, i)
            g_out(m, m, i) = Gdyn_m(1, m, m, i)
         enddo
         do m = Nc + 1, 2*Nc
            self_out(m, m, i) = S_ff_m(-1, m - Nc, m - Nc, i)
            g_out(m, m, i) = Gdyn_m(-1, m - Nc, m - Nc, i)
         enddo
      enddo

      if (rank == 0) then
         open (unit=67998, file='ed.fmos', form='unformatted')
         write (67998) Gdyn_m
         close (67998)
      endif

      !#####################################
      !#####################################
      !#####################################
      !  REAL AXIS CALCULATIONS
      !#####################################
      !#####################################

      if (allocated(Dmasig)) deallocate (Dmasig, Dmbsig, Dmcsig, DmasigTild, &
                                         DmbsigTild, DmcsigTild, Dmsig, AA, BB, CC, tempvec)
      if (allocated(I1a)) deallocate (I1a, I1b, I2a, I2b, I1c, I2c)

      allocate (I1a(-1:1, Nc, nw_r), I1b(-1:1, Nc, Nc, nw_r), I2a(-1:1, Nc, &
                                                                nw_r), I2b(-1:1, Nc, Nc, nw_r), I1c(-1:1, Nc, Nc, nw_r), I2c(-1:1, &
                                                                                                                      Nc, Nc, nw_r))
      allocate (Dmasig(-1:1, Nc, nw_r), Dmbsig(-1:1, Nc, Nc, nw_r))
      allocate (Dmcsig(-1:1, Nc, Nc, nw_r), DmasigTild(-1:1, Nc, nw_r), &
                DmbsigTild(-1:1, Nc, Nc, nw_r), DmcsigTild(-1:1, Nc, Nc, nw_r))
      allocate (Dmsig(-1:1, Nc, nw_r))
      allocate (AA(-1:1, Nc, nw_r), BB(-1:1, Nc, Nc, nw_r), CC(-1:1, Nc, Nc, &
                                                               nw_r), tempvec(nw_r))

      iter = 1

      write (*, *) 'FMOS REAL FREQUENCIES : ', frequ__(1:2)

      !----------------!
      ! static offdiag !
      !----------------!

      SIAM_T = 0.d0
      do s1 = -1, 1, 2
         SIAM_T(s1, 1:Nc, Nc + 1:Nc + Nb) = Vb(s1, 1:Nc, 1:Nb)
         SIAM_T(s1, :, :) = SIAM_T(s1, :, :) + transpose(conjg(SIAM_T(s1, :, &
                                                                      :)))
         SIAM_T(s1, 1:Nc, 1:Nc) = Ei_(s1, :, :)
         SIAM_T(s1, Nc + 1:Nc + Nb, Nc + 1:Nc + Nb) = Eb(s1, :, :)
         if (maxval(abs(SIAM_T(s1, :, :) - transpose(conjg(SIAM_T(s1, :, :) &
                                                           )))) > 1.d-4) then
            write (*, *) 'ERROR FMOS SIAM_T not HERMITIC'
            stop
         endif
         do i = 1, nw_r
            Gdyn_r(s1, :, :, i) = (frequ__(i) + mmu)*Id(Ntot) - SIAM_T(s1, :, &
                                                                       :)
            Gdyn_r(s1, 1:Nc, 1:Nc, i) = Gdyn_r(s1, 1:Nc, 1:Nc, i) - S_ff_r(s1, &
                                                                           1:Nc, 1:Nc, i)
            call invmat(Ntot, Gdyn_r(s1, :, :, i))
         enddo
      enddo

      if (iter == 1) then
         do s1 = -1, 1, 2
            do m = 1, Nc
               !call plotarray(real(frequ__), real(Gdyn_r(s1, m, m, :)),
               ! 'FMOSR_GREEN0_RE_ITER_' // trim(adjustl(toString(iter))) //
               ! "_" // trim(adjustl(toString(m))) // "_")
               !call plotarray(real(frequ__), aimag(Gdyn_r(s1, m, m, :)),
               ! 'FMOSR_GREEN0_IM_ITER_' // trim(adjustl(toString(iter))) //
               ! "_" // trim(adjustl(toString(m))) // "_")
            enddo
         enddo
      endif

      !-------------------------------------------------------------!
      ! Dmsig, Dmasig, Dmbsig, Dmcsig, DmasigTild, DmbsigTild, DmcsigTild !
      !-------------------------------------------------------------!

      do s1 = -1, 1, 2
         do m = 1, Nc
            Dmsig(s1, m, :) = 0.
            do k = 1, Nb
               Dmsig(s1, m, :) = Dmsig(s1, m, :) + conjg(Vb(s1, m, k))* &
                                 Vb(s1, m, k)/(frequ__(:) + mmu - Eb(s1, k, k))
            enddo
         enddo
      enddo

      do s1 = -1, 1, 2
         do m = 1, Nc
            !call plotarray(real(frequ__), real(1.d0/(frequ__(:) + mmu-Ei_(s1,
            ! m, m)-Dmsig(s1, m, :))), 'FMOSR_WEISS_RE_ITER_' //
            ! trim(adjustl(toString(iter))) // "_" //
            ! trim(adjustl(toString(m))) // "_")
            !call plotarray(real(frequ__), aimag(1.d0/(frequ__(:) +
            ! mmu-Ei_(s1, m, m)-Dmsig(s1, m, :))), 'FMOSR_WEISS_IM_ITER_' //
            ! trim(adjustl(toString(iter))) // "_" //
            ! trim(adjustl(toString(m))) // "_")
            !call plotarray(real(frequ__), real(Dmsig(s1, m, :)),
            ! 'FMOSR_HYBRID_RE_ITER_' // trim(adjustl(toString(iter))) // "_"
            ! // trim(adjustl(toString(m))) // "_")
            !call plotarray(real(frequ__), aimag(Dmsig(s1, m, :)),
            ! 'FMOSR_HYBRID_IM_ITER_' // trim(adjustl(toString(iter))) // "_"
            ! // trim(adjustl(toString(m))) // "_")
         enddo
      enddo

      do s1 = -1, 1, 2
         do m = 1, Nc
            Dmasig(s1, m, :) = 0.
            do k = 1, Nb
               Dmasig(s1, m, :) = Dmasig(s1, m, :) + conjg(Vb(-s1, m, k))* &
                                  Vb(-s1, m, k)/(frequ__(:) + mmu + Ei(-s1, m, m) - Ei(s1, &
                                                                                       m, m) - Eb(-s1, k, k))
            enddo
         enddo
      enddo

      do s1 = -1, 1, 2
         do m = 1, Nc
            do l = 1, Nc
               Dmbsig(s1, l, m, :) = 0.
               do k = 1, Nb
                  Dmbsig(s1, l, m, :) = Dmbsig(s1, l, m, :) + conjg(Vb(s1, l, &
                                                                       k))*Vb(s1, l, k)/(frequ__(:) + mmu + Ei(s1, l, &
                                                                                                   l) - Ei(s1, m, m) - Eb(s1, k, k))
               enddo
            enddo
         enddo
      enddo

      do s1 = -1, 1, 2
         do m = 1, Nc
            do l = 1, Nc
               Dmcsig(s1, l, m, :) = 0.
               do k = 1, Nb
                  Dmcsig(s1, l, m, :) = Dmcsig(s1, l, m, :) + conjg(Vb(-s1, l, &
                                                                       k))*Vb(-s1, l, k)/(frequ__(:) + mmu + Ei(-s1, l, &
                                                                                                  l) - Ei(s1, m, m) - Eb(-s1, k, k))
               enddo
            enddo
         enddo
      enddo

      do s1 = -1, 1, 2
         do m = 1, Nc
            DmasigTild(s1, m, :) = 0.
            do k = 1, Nb
               DmasigTild(s1, m, :) = DmasigTild(s1, m, :) + conjg(Vb(-s1, m, &
                                                                      k))*Vb(-s1, m, k)/(frequ__(:) + mmu + Eb(-s1, k, k) - &
                                                                                         Ei(-s1, m, m) - Ei(s1, m, m) - Umma(s1, m))
            enddo
         enddo
      enddo

      do s1 = -1, 1, 2
         do m = 1, Nc
            do l = 1, Nc
               DmbsigTild(s1, l, m, :) = 0.
               do k = 1, Nb
                  DmbsigTild(s1, l, m, :) = DmbsigTild(s1, l, m, :) + &
                                            conjg(Vb(s1, l, k))*Vb(s1, l, k)/(frequ__(:) + mmu &
                                                                           + Eb(s1, k, k) - Ei(s1, l, l) - Ei(s1, m, m) - Ummb(s1, &
                                                                                                                              l, m))
               enddo
            enddo
         enddo
      enddo

      do s1 = -1, 1, 2
         do m = 1, Nc
            do l = 1, Nc
               DmcsigTild(s1, l, m, :) = 0.
               do k = 1, Nb
                  DmcsigTild(s1, l, m, :) = DmcsigTild(s1, l, m, :) + &
                                            conjg(Vb(-s1, l, k))*Vb(-s1, l, k)/(frequ__(:) + &
                                                                              mmu + Eb(-s1, k, k) - Ei(-s1, l, l) - Ei(s1, m, m) - &
                                                                                Ummc(s1, l, m))
               enddo
            enddo
         enddo
      enddo

      !---------------------------!
      ! I1a, I2a, I1b, I2b, I1c, I2c   !
      !---------------------------!

      do s1 = -1, 1, 2
         do m = 1, Nc
            I1a(s1, m, :) = 0.
            do k = 1, Nb
               I1a(s1, m, :) = I1a(s1, m, :) + conjg(Vb(-s1, m, k))* &
                               densav(-s1, m, Nc + k)/(frequ__(:) + mmu + Ei(-s1, m, &
                                                                             m) - Ei(s1, m, m) - Eb(-s1, k, k))
               I1a(s1, m, :) = I1a(s1, m, :) - Vb(-s1, m, k)*densav(-s1, Nc &
                                                                    + k, m)/(frequ__(:) + mmu + Eb(-s1, k, k) - Ei(-s1, m, &
                                                                                                    m) - Ei(s1, m, m) - Umma(s1, m))
            enddo
         enddo
      enddo

      do s1 = -1, 1, 2
         do m = 1, Nc
            I2a(s1, m, :) = 0.
            do k = 1, Nb
               do kp = 1, Nb
                  I2a(s1, m, :) = I2a(s1, m, :) - conjg(Vb(-s1, m, k))* &
                                  Vb(-s1, m, k)*densav(-s1, Nc + kp, Nc + k)/( &
                                  frequ__(:) + mmu + Ei(-s1, m, m) - Ei(s1, m, m) - Eb(-s1, &
                                                                                       k, k))
                  I2a(s1, m, :) = I2a(s1, m, :) - conjg(Vb(-s1, m, k))* &
                                  Vb(-s1, m, k)*densav(-s1, Nc + k, Nc + kp)/( &
                                  frequ__(:) + mmu + Eb(-s1, k, k) - Ei(-s1, m, m) - Ei(s1, &
                                                                                        m, m) - Umma(s1, m))
               enddo
            enddo
         enddo
      enddo

      do s1 = -1, 1, 2
         do m = 1, Nc
            do l = 1, Nc
               I1b(s1, m, l, :) = 0.
               do k = 1, Nb
                  I1b(s1, m, l, :) = I1b(s1, m, l, :) + conjg(Vb(s1, l, k))* &
                                     densav(s1, l, Nc + k)/(frequ__(:) + mmu + Ei(s1, l, &
                                                                                  l) - Ei(s1, m, m) - Eb(s1, k, k))
                  I1b(s1, m, l, :) = I1b(s1, m, l, :) - Vb(s1, l, k)*densav( &
                                     s1, Nc + k, l)/(frequ__(:) + mmu + Eb(s1, k, &
                                                                           k) - Ei(s1, l, l) - Ei(s1, m, m) - Ummb(s1, l, m))
               enddo
            enddo
         enddo
      enddo

      do s1 = -1, 1, 2
         do m = 1, Nc
            do l = 1, Nc
               I2b(s1, m, l, :) = 0.
               do k = 1, Nb
                  do kp = 1, Nb
                     I2b(s1, m, l, :) = I2b(s1, m, l, :) - conjg(Vb(s1, l, k)) &
                                        *Vb(s1, l, k)*densav(s1, Nc + kp, Nc + k)/( &
                                        frequ__(:) + mmu + Ei(s1, l, l) - Ei(s1, m, m) - Eb(s1, &
                                                                                            k, k))
                     I2b(s1, m, l, :) = I2b(s1, m, l, :) - conjg(Vb(s1, l, k)) &
                                        *Vb(s1, l, k)*densav(s1, Nc + k, Nc + kp)/( &
                                        frequ__(:) + mmu + Eb(s1, k, k) - Ei(s1, m, m) - Ei(s1, &
                                                                                            l, l) - Ummb(s1, l, m))
                  enddo
               enddo
            enddo
         enddo
      enddo

      do s1 = -1, 1, 2
         do m = 1, Nc
            do l = 1, Nc
               I1c(s1, m, l, :) = 0.

               do k = 1, Nb
                  I1c(s1, m, l, :) = I1c(s1, m, l, :) + conjg(Vb(-s1, l, k))* &
                                     densav(-s1, l, Nc + k)/(frequ__(:) + mmu + Ei(-s1, &
                                                                                   l, l) - Ei(s1, m, m) - Eb(-s1, k, k))
                  if (.not. use_def_feng) then
                     I1c(s1, m, l, :) = I1c(s1, m, l, :) - Vb(-s1, l, k)* &
                                        densav(-s1, Nc + k, l)/(frequ__(:) + mmu + &
                                                                Eb(-s1, k, k) - Ei(-s1, l, l) - Ei(s1, m, m) - Ummc(s1, l, &
                                                                                                                    m))
                  endif
               enddo

               if (use_def_feng) then
                  do k = 1, Nb/Nc
                     i = (m - 1)*mstep
                     j = (l - 1)*mstep
                     I1c(s1, m, l, :) = I1c(s1, m, l, :) - Vbath(-s1, m, k)* &
                                        densav(-s1, Nc + k + j, l)/(frequ__(:) + mmu + &
                                                                    Eb(-s1, k + j, k + j) - Ei(-s1, l, l) - Ei(s1, m, &
                                                                                                               m) - Ummc(s1, l, m))
                  enddo
               endif

            enddo
         enddo
      enddo

      do s1 = -1, 1, 2
         do m = 1, Nc
            do l = 1, Nc
               I2c(s1, m, l, :) = 0.
               do k = 1, Nb
                  do kp = 1, Nb
                     I2c(s1, m, l, :) = I2c(s1, m, l, :) - conjg(Vb(-s1, l, &
                                                                    k))*Vb(-s1, l, k)*densav(-s1, Nc + kp, Nc + k)/ &
                                        (frequ__(:) + mmu + Ei(-s1, l, l) - Ei(s1, m, &
                                                                               m) - Eb(-s1, k, k))
                     I2c(s1, m, l, :) = I2c(s1, m, l, :) - conjg(Vb(-s1, l, &
                                                                    k))*Vb(-s1, l, k)*densav(-s1, Nc + k, Nc + kp)/ &
                                        (frequ__(:) + mmu + Eb(-s1, k, k) - Ei(-s1, l, &
                                                                               l) - Ei(s1, m, m) - Ummc(s1, l, m))
                  enddo
               enddo
            enddo
         enddo
      enddo

      !----------!
      ! AA, BB, CC !
      !----------!

      if (fmos_fluc) then
         do s1 = -1, 1, 2
            do m = 1, Nc
               AA(s1, m, :) = 0.
               do l = 1, Nc
                  if (l /= m) then
                     AA(s1, m, :) = AA(s1, m, :) + UU(s1, s1, l, m)*densav(s1, &
                                                                           l, l) + UU(s1, -s1, l, m)*densav(-s1, l, l) + &
                                    2.0*UU(1, -1, m, m)*densav(s1, l, l)*densav(-s1, l, &
                                                                                l)
                     AA(s1, m, :) = AA(s1, m, :) + densav(-s1, m, m)*( &
                                    Dmsig(-s1, l, :) + Dmsig(s1, l, :))
                  endif
               enddo
               AA(s1, m, :) = -AA(s1, m, :) + frequ__(:) + mmu - Ei(s1, m, &
                                                                    m) - ueff(s1, m) - Dmsig(s1, m, :) - Dmasig(s1, m, &
                                                                                                           :) - DmasigTild(s1, m, :)
               AA(s1, m, :) = ueff(s1, m)/AA(s1, m, :)
            enddo
         enddo
      else
         do s1 = -1, 1, 2
            do m = 1, Nc
               AA(s1, m, :) = 0.
               AA(s1, m, :) = frequ__(:) + mmu - Ei(s1, m, m) - ueff(s1, m) - &
                              Dmsig(s1, m, :) - Dmasig(s1, m, :) - DmasigTild(s1, m, :)
               AA(s1, m, :) = ueff(s1, m)/AA(s1, m, :)
            enddo
         enddo
      endif

      do s1 = -1, 1, 2
         do m = 1, Nc
            do l = 1, Nc
               BB(s1, m, l, :) = 0.0
               if (l /= m) then
                  do lp = 1, Nc
                     if (lp /= l) then
                        BB(s1, m, l, :) = BB(s1, m, l, :) + (UU(s1, s1, lp, &
                                                                m)*densav(s1, lp, lp) + UU(s1, -s1, lp, &
                                                                                        m)*densav(-s1, lp, lp) + 2.0*UU(1, -1, lp, &
                                                                                         lp)*densav(s1, lp, lp)*densav(-s1, lp, lp))
                        BB(s1, m, l, :) = BB(s1, m, l, :) + (densav(-s1, lp, &
                                                                    lp) + densav(s1, lp, lp))*(Dmsig(s1, m, :))
                     endif
                  enddo
                  BB(s1, m, l, :) = -BB(s1, m, l, :) + frequ__(:) + mmu - Ei(s1, &
                                                                             m, m) - UUeff(s1, s1, l, m) - (UU(1, -1, m, &
                                                                           m)*densav(-s1, m, m) + UU(-s1, s1, l, m)*densav(-s1, l, &
                                                                           l) + UU(1, -1, m, m)*densav(s1, m, m)*densav(-s1, m, m) &
                                                                             + UU(1, -1, l, l)*densav(s1, l, l)*densav(-s1, l, l)) &
                                    - (Dmsig(s1, m, :) - Dmbsig(s1, l, m, :) - DmbsigTild(s1, &
                                                                             l, m, :)) - (densav(-s1, m, m) + densav(-s1, l, l))*( &
                                    Dmsig(s1, m, :))
                  BB(s1, m, l, :) = UUeff(s1, s1, l, m)/BB(s1, m, l, :)
               endif
            enddo
         enddo
      enddo

      do s1 = -1, 1, 2
         do m = 1, Nc
            do l = 1, Nc
               CC(s1, m, l, :) = 0.
               if (l /= m) then
                  do lp = 1, Nc
                     if (lp /= l) then
                        CC(s1, m, l, :) = CC(s1, m, l, :) + (UU(s1, s1, lp, &
                                                                m)*densav(s1, lp, lp) + UU(s1, -s1, lp, &
                                                                                        m)*densav(-s1, lp, lp) + 2.0*UU(1, -1, lp, &
                                                                                         lp)*densav(s1, lp, lp)*densav(-s1, lp, lp))
                        CC(s1, m, l, :) = CC(s1, m, l, :) + (densav(-s1, lp, &
                                                                    lp) + densav(s1, lp, lp))*(Dmsig(s1, m, :))
                     endif
                  enddo
                  CC(s1, m, l, :) = -CC(s1, m, l, :) + frequ__(:) + mmu - Ei(s1, &
                                                                             m, m) - UUeff(-s1, s1, l, m) - (UU(1, -1, m, &
                                                                             m)*densav(-s1, m, m) + UU(s1, s1, l, m)*densav(s1, l, &
                                                                           l) + UU(1, -1, m, m)*densav(s1, m, m)*densav(-s1, m, m) &
                                                                             + UU(1, -1, l, l)*densav(s1, l, l)*densav(-s1, l, l)) &
                                    - (Dmsig(s1, m, :) - Dmcsig(s1, l, m, :) - DmcsigTild(s1, &
                                                                              l, m, :)) - (densav(-s1, m, m) + densav(s1, l, l))*( &
                                    Dmsig(s1, m, :))
                  CC(s1, m, l, :) = UUeff(-s1, s1, l, m)/CC(s1, m, l, :)
               endif
            enddo
         enddo
      enddo

      if (.not. fmos_fluc) then
         BB = 0.d0
         CC = 0.d0
      endif

      do s1 = -1, 1, 2
         Gdyn_r(s1, 1:Nc, 1:Nc, :) = 0.0
         do m = 1, Nc
            do l = 1, Nc
               if (l /= m) then
                  Gdyn_r(s1, m, m, :) = Gdyn_r(s1, m, m, :) + BB(s1, m, l, &
                                                                 :)*(Dmsig(s1, m, :)*I1b(s1, m, l, :) + I2b(s1, m, l, &
                                                                           :)) + CC(s1, m, l, :)*(Dmsig(s1, m, :)*I1c(s1, m, l, :) &
                                                                                                                 + I2c(s1, m, l, :))
               endif
            enddo
            Gdyn_r(s1, m, m, :) = -Gdyn_r(s1, m, m, :) + frequ__(:) + &
                                  mmu - Ei(s1, m, m) - Dmsig(s1, m, :) - AA(s1, m, :)*(Dmsig(s1, &
                                                                                             m, :)*I1a(s1, m, :) + I2a(s1, m, :))
            tempvec = 0.0
            do l = 1, Nc
               if (l /= m) then
                  tempvec(:) = tempvec(:) + BB(s1, m, l, :)*(densav(s1, l, l) &
                                                             + I1b(s1, m, l, :)) + CC(s1, m, l, :)*(densav(-s1, l, &
                                                                                                           l) + I1c(s1, m, l, :))
               endif
            enddo
            Gdyn_r(s1, m, m, :) = 1.d0/Gdyn_r(s1, m, m, :)*(1.d0 + AA(s1, &
                                                                      m, :)*(densav(-s1, m, m) + I1a(s1, m, :)) + tempvec(:))
         enddo
      enddo

      do s1 = -1, 1, 2
         S_ff_r(s1, :, :, :) = 0.0
         do m = 1, Nc
            !call plotarray(real(frequ__), real(Gdyn_r(s1, m, m, :)),
            ! 'FMOSR_GREEN_RE_ITER_' // trim(adjustl(toString(iter))) // "_"
            ! // trim(adjustl(toString(m))) // "_")
            !call plotarray(real(frequ__), aimag(Gdyn_r(s1, m, m, :)),
            ! 'FMOSR_GREEN_IM_ITER_' // trim(adjustl(toString(iter))) // "_"
            ! // trim(adjustl(toString(m))) // "_")
            S_ff_r(s1, m, m, :) = frequ__(:) + mmu - Ei_(s1, m, m) - Dmsig(s1, &
                                                                           m, :) - 1.d0/Gdyn_r(s1, m, m, :)
            !call plotarray(real(frequ__), real(S_ff_r(s1, m, m, :)),
            ! 'FMOSR_SELF_RE_ITER_' // trim(adjustl(toString(iter))) // "_" //
            ! trim(adjustl(toString(m))) // "_")
            !call plotarray(real(frequ__), aimag(S_ff_r(s1, m, m, :)),
            ! 'FMOSR_SELF_IM_ITER_' // trim(adjustl(toString(iter))) // "_" //
            ! trim(adjustl(toString(m))) // "_")
            do i = 1, nw_r
               if (aimag(S_ff_r(s1, m, m, i)) > 0.) S_ff_r(s1, m, m, i) = &
                  CMPLX(real(S_ff_r(s1, m, m, i)), 0.d0)
            enddo
         enddo
      enddo

      sigw = 0.
      gw = 0.
      do i = 1, nw_r
         do m = 1, Nc
            sigw(m, m, i) = S_ff_r(1, m, m, i)
            gw(m, m, i) = Gdyn_r(1, m, m, i)
         enddo
         do m = Nc + 1, 2*Nc
            sigw(m, m, i) = S_ff_r(-1, m - Nc, m - Nc, i)
            gw(m, m, i) = Gdyn_r(-1, m - Nc, m - Nc, i)
         enddo
      enddo

      return
   end subroutine

   subroutine add_substract_tail(nw, Giom, iom, ahh, ddd, add, tail)

      use genvar, only: imi

      implicit none

      integer    :: nw
      complex(kind=DP) :: Giom(1:nw), temp(1:nw)
      real(kind=DP)    :: iom(1:nw), ah, ahh, drsign, ddd, tail
      integer    :: t, n, i, j
      logical    :: add

      if (add) then
         drsign = 1.d0
      else
         drsign = -1.d0
      endif
      temp = 0.d0
      do n = 1, nw
         if (n /= 0) then
            if (abs(iom(n)) > tail) then
               temp(n) = 0.
            else
               temp(n) = Giom(n) + ddd*drsign/(imi*iom(n) + ahh)
            endif
         endif
      enddo

      Giom = temp

   end subroutine

   subroutine get_tail_parameters(nw, Giom, iom, ahh, ddd, green)

      implicit none

      integer    :: nw, kk
      logical    :: green
      complex(kind=DP) :: Giom(1:nw)
      real(kind=DP)    :: iom(1:nw), ddd, ahh

      ddd = get_ddd(nw, iom, Giom, green)
      ahh = get_ah(nw, iom, Giom)
      ahh = ahh/ddd
      return
   end subroutine

   elemental real(kind=DP) function step_func_loc(rr)

      implicit none

      real(kind=DP), intent(in) :: rr

      if (rr > 0.d0) then
         step_func_loc = 1.d0
      else
         step_func_loc = 0.d0
      endif

   end function

   subroutine Fourier(nw, Giom, iom, Gtau_, beta, tailmax, green)

      use linalg, only: dexpc, MPLX

      implicit none

      integer, parameter:: Ntau = 1
      integer           :: nw, kk
      complex(kind=DP)        :: Giom(1:nw), Giomback(1:nw)
      real(kind=DP)           :: iom(1:nw), tau(Ntau), ah, ahh, mm
      real(kind=DP)           :: Gtau_
      real(kind=DP)           :: beta, df, temp
      complex(kind=DP)        :: csum
      real(kind=DP)           :: dww, bb, ddd, ahh2, ddd2, g1, gn, tailmax
      integer           :: t, n, i, j, itail
      real(16)          :: t1, t2, pp1, pp2
      logical           :: green

      tau(1) = -0.0001d0

      Giomback = Giom

      call get_tail_parameters(nw, Giom, iom, ahh, ddd, green)
      call add_substract_tail(nw, Giom, iom, ahh, ddd, add=.false., tail= &
                              tailmax)

      do t = 1, Ntau
         csum = 0.d0
         do n = 1, nw
            csum = csum + MPLX(-iom(n)*tau(t))*Giom(n)
         enddo
         Gtau_ = csum/beta*2.d0
      enddo

      do t = 1, Ntau
         if (ahh > 0.d0) then
            pp1 = DEXPc(ahh*(tau(t) - beta))
            pp1 = pp1/(DEXPc(-beta*ahh) + 1.d0)
            pp2 = DEXPc(ahh*(tau(t)))
            pp2 = pp2/(DEXPc(-beta*ahh) + 1.d0)
         else
            pp1 = DEXPc(ahh*tau(t))
            pp1 = pp1/(DEXPc(beta*ahh) + 1.d0)
            pp2 = DEXPc(ahh*(tau(t) + beta))
            pp2 = pp2/(DEXPc(beta*ahh) + 1.d0)
         endif
         Gtau_ = Gtau_ - ddd*(step_func_loc(tau(t))*pp1 - &
                              step_func_loc(-tau(t))*pp2)
      enddo

      Giom = Giomback

      return
   end subroutine

   real(kind=DP) function get_ddd(nw, iom, Giom, green_)

      implicit none

      integer    :: nw
      logical    :: green_
      real(kind=DP)    :: iom(1:nw), ddd
      complex(kind=DP) :: Giom(1:nw)

      if (green_) then
         get_ddd = 1.d0
      else
         get_ddd = -Aimag(Giom(nw))*iom(nw)
      endif
   end function

   real(kind=DP) function get_ah(nw, iom, Giom)

      implicit none

      integer    :: nw
      real(kind=DP)    :: iom(1:nw), ddd
      complex(kind=DP) :: Giom(1:nw)
      integer    :: i, kk, p, ii, k

      get_ah = Real(Giom(nw))*(iom(nw)**2)
   end function
