
!======================!
!======================!

program dmft_check_convergence
!---------------!
   ! use DMFT_SOLVER_ED
   !use dmftmod
   use genvar
   use init_and_close_my_sim
   ! use fortran_cuda
   use matrix
   use StringManip
   !use DLPLOT
   !use plotlib
   !use plot2d
!---------------!
   implicit none

   integer                :: kk_, i, j, k, l, channels, n_frequ, n_frequ_long
   real(8), allocatable    :: COEF(:, :), rotmat(:, :), temp(:), vvv(:, :), ddens(:)
   complex(8), allocatable :: frequ_(:), green_lat(:, :, :), green_matsu(:, :, :)
   logical                :: check
   real(8)                :: mmu, vv, aa_matsu, aa
   complex(8)             :: ttt
   complex(8), allocatable :: tt(:, :)
   integer                :: paramagnet, fac, mmax
   logical                :: checkfile
   real(8)                :: beta

   call initialize_my_simulation
   testing = .false.
   fast_invmat = .true.
   use_cuda_routines = .false.
   use_cula_routines = .false.
   force_invmat_single_prec = .false.
   use_openmp_invmat = .false.
   diag_use_LU_instead_of_pivot = .false.
   flag_use_invmat_jordan = .true.
   flag_use_invmat_jordan_real = .true.
   enable_mpi_dot = .false.

   write (*, *) 'is it a paramagnet ? (=1 for true)'
   read (*, *) paramagnet

   fac = 4
   if (paramagnet == 1) fac = fac/2
   mmax = 2
   if (paramagnet == 1) mmax = 1

   do kk_ = 1, mmax
      open (unit=10012, file="green_onetep_spin_check_convergence_ren"//TRIM(ADJUSTL(toString(kk_))), form='unformatted')
      read (10012) n_frequ, channels
     if(.not.allocated(temp)) allocate(temp(fac*channels),vvv(channels,2),rotmat(channels,2),COEF(channels,2),frequ_(n_frequ),green_lat(n_frequ,channels,2))
      write (*, *) 'Reading input file, nfrequ/channels : ', n_frequ, channels
      read (10012) vvv(:, kk_)
      write (*, *) 'vvv : ', vvv(:, kk_)
      read (10012) mmu
      write (*, *) 'mu  : ', mmu
      read (10012) COEF(:, kk_)
      write (*, *) 'COEF: ', COEF(:, kk_)
      read (10012) rotmat(:, kk_)
      write (*, *) 'Rot : ', rotmat(:, kk_)
      do i = 1, n_frequ
         write (*, *) 'reading frequ :', i
         read (10012) frequ_(i), green_lat(i, :, kk_)
         green_lat(i, :, kk_) = green_lat(i, :, kk_)*(rotmat(:, kk_)**2)
         write (*, *) 'done...'
      enddo
      close (10012)
   enddo

   INQUIRE (file='green_output_matsu', exist=checkfile)
   if (checkfile) then
      open (unit=1000, file='green_output_matsu')
      j = 0
      do
         read (1000, *, end=212)
         j = j + 1
      enddo
212   continue; n_frequ_long = j; rewind (1000); 
      allocate (green_matsu(n_frequ_long, channels, 2))
      do i = 1, n_frequ_long
         read (1000, *) vv, (temp(k), k=1, fac*channels)
         do j = 1, mmax
         do k = 1, channels
            green_matsu(i, k, j) = temp(channels*(j - 1)*2 + 2*k - 1) + imi*temp(channels*(j - 1)*2 + 2*k)
            green_matsu(i, k, j) = green_matsu(i, k, j)*(rotmat(k, j)**2)
         enddo
         enddo
      enddo
      close (1000)
   endif

   INQUIRE (file='green_output', exist=checkfile)
   if (checkfile) then
      open (unit=1000, file='green_output')
      read (1000, *)
      j = 0
      do
         read (1000, *, end=213)
         j = j + 1
      enddo
213   continue; n_frequ_long = j; rewind (1000); read (1000, *)
      allocate (green_matsu(n_frequ_long, channels, 2))
      do i = 1, n_frequ_long
         read (1000, *) vv, (temp(k), k=1, fac*channels)
         do j = 1, mmax
         do k = 1, channels
            green_matsu(i, k, j) = temp(channels*(j - 1)*2 + 2*k - 1) + imi*temp(channels*(j - 1)*2 + 2*k)
            green_matsu(i, k, j) = green_matsu(i, k, j)*(rotmat(k, j)**2)
         enddo
         enddo
      enddo
      close (1000)
   endif

   if (.not. allocated(green_matsu)) then
      allocate (green_matsu(n_frequ_long, channels, 2))
      green_matsu = green_lat
      write (*, *) '=========== NOT OUTPUT FILE FROM SOLVER, USES LATTICE GF INSTEAD ==========='
   endif

   beta = dacos(-1.d0)/aimag(frequ_(1))

   write (*, *)
   write (*, *) 'COEF up       : ', COEF(:, 1)
   write (*, *) 'Chem          : ', mmu
   write (*, *) 'beta          : ', beta

   open (unit=100, file='rho_fermi_level_from_gimp')
   open (unit=101, file='rho_fermi_level_from_gproj')
   open (unit=102, file='rho_n_tot_from_gimp')

   if (allocated(tt)) deallocate (tt, ddens)
   allocate (tt(channels, n_frequ), ddens(channels))

   do kk_ = 1, mmax

      do k = 1, channels
         open (1010, file='compare_imp_'//TRIM(ADJUSTL(toString(k)))//"_"//TRIM(ADJUSTL(toString(kk_))))
         open (1011, file='compare_lat_'//TRIM(ADJUSTL(toString(k)))//"_"//TRIM(ADJUSTL(toString(kk_))))
         do i = 1, n_frequ
            tt(k, i) = green_matsu(i, k, kk_)
            ttt = green_lat(i, k, kk_)
            if (i == 1 .or. i == 2) then
               write (100, *) real(-aimag(tt(k, i))/dacos(-1.d0)*COEF(k, kk_))
               write (101, *) real(-aimag(ttt)/dacos(-1.d0)*COEF(k, kk_))
            endif
            write (1010, *) aimag(frequ_(i)), real(tt(k, i)), aimag(tt(k, i))
            write (1011, *) aimag(frequ_(i)), real(ttt), aimag(ttt)
         enddo
         close (1010)
         close (1011)
      enddo

      do k = 1, channels
         ddens(k) = get_dens(real(green_matsu(:, k, kk_)), aimag(green_matsu(:, k, kk_)), diag=.true.)*COEF(k, kk_)
         write (102, *) ddens(k)
      enddo
      write (*, *) 'TOTAL DENSITY spin [x] : ', kk_, sum(ddens)

   enddo

   close (100)
   close (101)
   close (102)

   call finalize_my_simulation

   write (*, *) 'CHECK_CONVERGENCE_FINISHED'

contains

!------------------------!
!------------------------!
!------------------------!
!------------------------!
!------------------------!

   real(8) function get_dens(GlocRe, GlocIm, diag)
      use linalg, only: mplx
      implicit none
      integer    :: k1, k2, i, j, k, jj
      complex(8) :: dens
      real(8)    :: frequ, pi, tau0, omega, df, GlocRe(:), GlocIm(:)
      complex(8) :: temp
      logical    :: diag
      real(8)    :: alpha

      pi = dacos(-1.d0)
      jj = size(GlocIm) - 7
      frequ = pi/beta*(2.d0*dble(jj) - 1)
      alpha = -GlocIm(jj)*frequ
      write (*, *) 'GET DENS ALPHA COEF : ', alpha

      tau0 = 1.d-8; dens = 0.d0
      do j = 1, n_frequ_long
         omega = pi/beta*(2.d0*dble(j) - 1.d0)
         if (abs(omega) < 1.d-9) omega = 1.d-9
         temp = CMPLX(GlocRe(j), GlocIm(j), kind=8)
         if (diag) then
            dens = dens + 2.d0/beta*(MPLX(omega*tau0)*(temp + imi/omega*alpha))
         else
            dens = dens + 2.d0/beta*(MPLX(omega*tau0)*(temp))
         endif
      enddo
      if (diag) then
         dens = dens + 0.5*alpha
      endif
      get_dens = dens

      return
   end function

!------------------------!
!------------------------!
!------------------------!
!------------------------!

end program

!======================!
!======================!

