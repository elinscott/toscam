!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!

module local_vars

   use genvar, only: DP
   use mpi
   use mesh, only: findin1dmesh

   implicit none

   real(kind=DP), parameter      ::  me_in_ev = 510998.91d0
   real(kind=DP), parameter      ::  small_scale_derivative = 0.005d0, hartree_in_eV = 27.211396132d0, eV_in_hartree = 0.036749308824d0
   real(kind=DP), parameter      ::  pi = 3.14159265358979323846264338327950288419716939937510d0
   real(kind=DP), parameter      ::  one_angstrom_in_bohr_radius = 1.889725989d0
   real(kind=DP), parameter      ::  hbar_divided_electron_mass_in_meter2_per_second = 0.000115767622d0
   real(kind=DP), parameter      ::  hbar_divided_by_second_in_ev = 6.58211899131*1.d-16
   real(kind=DP), parameter      ::  hbar_divided_by_me_in_Angstrom_square_divided_by_second = 1.15767635048*1.d16
   real(kind=DP), parameter      ::  one_over_Angstrom_in_one_over_cm = 100000000.d0
   complex(kind=DP), parameter   ::  imi = CMPLX(0.d0, 1.0d0, 8)
   real(kind=DP), parameter      ::  MAX_EXP = 700.d0
   real(kind=DP), parameter      ::  MIN_EXP = -700.d0
   integer                ::  rank, size2, ierr
   real(kind=DP), parameter      ::  Angst_in_ev = 1968.d0
   real(kind=DP), parameter      ::  Angst_in_cm = 1.d-8
   real(kind=DP), parameter      ::  hbar_divided_e2_in_ohm = 25812.d0/2.d0/pi
   real(kind=DP), parameter      ::  cm_moins_1_in_ev = 0.000123980262342235
   real(kind=DP), parameter      ::  Kelvin_in_ev = 8.617343d-5
   real(kind=DP), parameter      ::  ev_in_Kelvin = 11604.505008098

contains

   subroutine mpisum(mat)
      implicit none
      real(kind=DP)  :: mat(:, :, :, :), matm(size(mat(:, 1, 1, 1)), size(mat(1, :, 1, 1)), size(mat(1, 1, :, 1)), size(mat(1, 1, 1, :)))
      integer  :: i, j, k, s(4), p
      if (size2 == 1) then
         return
      endif
      s = shape(mat)
      matm = 0.
      p = s(1)*s(2)*s(3)*s(4)
      call MPI_ALLREDUCE(mat, matm, p, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      mat = matm
      return
   end subroutine

   subroutine mpisum_(mat)
      implicit none
      real(kind=DP)  :: mat(:, :), matm(size(mat(:, 1)), size(mat(1, :)))
      integer :: i, j, k, s1, s2
      if (size2 == 1) then
         return
      endif
      s1 = size(mat(:, 1)); s2 = size(mat(1, :)); matm = 0.
      call MPI_ALLREDUCE(mat, matm, s1*s2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      mat = matm
      return
   end subroutine

   INTEGER Function StrInt2(ch)
      character*(*) ch
      integer i, j, ilen
      j = 0; 
      ilen = LEN_TRIM(ch)
      do i = 1, ilen !+1 !BUG corrected
         if (i .eq. ilen) then
            j = j + IACHAR(ch(i:i)) - 48
         else
            j = j + INT((IACHAR(ch(i:i)) - 48)*10**(ilen - i))
         end if
      end do
      StrInt2 = j
   END Function

   subroutine initString(STRING)
      character*(*) :: STRING
      integer :: i
      do i = 1, len(STRING)
         STRING(i:i) = " "
      enddo
   end subroutine

   character*30 function toString(i)
      integer :: i
      call initString(toString)
      write (toString, '(i8)') i
      toString = TRIM(ADJUSTL(toString))
   end function

   elemental real(kind=DP) function fermi_dirac_(rr, mu, T)
      implicit none
      real(kind=DP), intent(in) :: rr, T, mu
      real(kind=DP)            :: aa
      aa = (rr - mu)/T
      if (aa < -100.d0) then
         fermi_dirac_ = 1.d0
         return
      endif
      if (aa > 100.d0) then
         fermi_dirac_ = 0.d0
         return
      endif
      fermi_dirac_ = 1.d0/(1.d0 + DEXPc(aa))
   end function

   elemental real(kind=DP) function derivative_fermi_dirac_(rr, mu, T)
      implicit none
      real(kind=DP), intent(in) :: rr, T, mu
      real(kind=DP)            :: aa
      aa = (rr - mu)/T/2.d0
      if (abs(aa) > 200.d0) then
         derivative_fermi_dirac_ = 0.d0
         return
      endif
      derivative_fermi_dirac_ = -1.d0/T/(DEXPc(-aa) + DEXPc(aa))**2
   end function

   pure real(kind=DP) function DEXPc(rr)
      implicit none
      real(kind=DP), intent(in) :: rr
      if (rr < MAX_EXP) then
         if (rr < MIN_EXP) then
            DEXPc = 0.d0
         else
            DEXPc = EXP(rr)
         endif
      else
         DEXPc = EXP(MAX_EXP)
      endif
   end function

   subroutine initialize_my_simulation
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, size2, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   end subroutine

   subroutine finalize_my_simulation
      call MPI_FINALIZE(ierr)
   end subroutine

end module

!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!
!================================================!

program optical_cond

   use local_vars
   ! use fortran_cuda

   implicit none

   integer                :: imin, imax, directions, ttt, k1, k2, nn, i1, i2, nfrequ, ii, ien, i, j, k, l, is, mpi_rank, mpi_size
   logical                :: check, usexy
   character(50)          :: filename
   complex(kind=DP)             :: energy
   real(kind=DP)                :: units_for_optics
   real(kind=DP)                :: frequ_temp, tttemp(4), a1, a2, a3, a4, Volume
   real(kind=DP)                :: aa(3), bb(3), cc(3), t1(3), t2(3), t3(3), F_factor, dstep, TT, rtemp, chem_shift, fermi_e, mmu
   real(kind=DP), allocatable    :: green(:, :, :), greeny(:, :, :), vertex__(:, :, :)
   real(kind=DP), allocatable    :: fermiv(:), derfermiv(:), frequ(:), Omega(:), optics(:, :), opticstot(:, :)
   character(300)         :: value(100)
   integer, allocatable    :: i_plus_j(:, :)
   integer                :: num, len__, status__, kp, nkpoints, nktemp
   real(kind=DP)                :: normk, nktempd

   !-------------------------------!
   call initialize_my_simulation
   !-------------------------------!

!###########################################################!
!###########################################################!
!###########################################################!
!###########################################################!
!###########################################################!

   call system(" rm scratchdir 2> /dev/null ")

   filename = trim(adjustl(set_name(1, 1, 1)))
   call check_file
   open (unit=30, file=trim(adjustl(filename)), form='unformatted')
   read (30) nkpoints, normk, directions, Volume, TT, nfrequ, mpi_size, nn

   do is = 1, 2
   do kp = 1, nkpoints

      write (*, *) '-----------------------------------'
      write (*, *) ' doing SPIN : ', is
      write (*, *) ' kpoint     : ', kp
      write (*, *) '-----------------------------------'

      filename = trim(adjustl(set_name(is, kp, 1)))
      call check_file
      open (unit=30, file=trim(adjustl(filename)), form='unformatted')
      read (30) nktemp, normk, directions, Volume, TT, nfrequ, mpi_size, nn

      TT = abs(TT)
      write (*, *) 'how many files should be present : ', mpi_size
      write (*, *) 'size of matrices                 : ', nn
      write (*, *) 'total number of frequencies      : ', nfrequ
      write (*, *) 'for temperature                  : ', TT
      write (*, *) 'And unitcell Volume (Ang^3)      : ', Volume
      close (30)

      usexy = .false.

      imin = 100000
      imax = -100000
      do mpi_rank = 1, mpi_size
         filename = set_name(is, kp, mpi_rank)
         call check_file
         write (*, *) 'reading file : ', filename
         open (unit=30, file=trim(adjustl(filename)), form='unformatted')
         do
            read (30, end=1024) nktemp, nktempd, directions, rtemp, rtemp, ii, ii, ii, ien
            if (directions /= 1) then
               usexy = .true.
            endif
            if (ien < imin) imin = ien
            if (ien > imax) imax = ien
         enddo
1024     continue
         write (*, *) 'IMIN,IMAX IN FILE : ', imin, imax, mpi_rank
         close (30)
      enddo

      write (*, *) 'SELECTION OF FREQUENCIES PRESENT IN FILES : ', imin, imax

      if (.not. allocated(green)) then
         allocate (green(imin:imax, nn, nn))
         allocate (frequ(imin:imax), fermiv(imin:imax), derfermiv(imin:imax))
         green = 0.d0
         frequ = 0.d0
         derfermiv = 0.d0
         fermiv = 0.d0
         if (usexy) then
            allocate (greeny(imin:imax, nn, nn))
            greeny = 0.d0
         endif
      endif

      !##################################!
      do mpi_rank = 1, mpi_size
         filename = set_name(is, kp, mpi_rank)
         call check_file
         write (*, *) 'reading file : ', filename
         open (unit=30, file=trim(adjustl(filename)), form='unformatted')
         !---------------------------------!
         do
         read (30, end=24) nktemp, nktempd, directions, rtemp, rtemp, ii, ii, ii, ien, energy, fermi_e, chem_shift, green(ien, :, :)
            if (usexy) then
        read (30, end=24) nktemp, nktempd, directions, rtemp, rtemp, ii, ii, ii, ien, energy, fermi_e, chem_shift, greeny(ien, :, :)
            endif
            frequ(ien) = real(energy) - (fermi_e + chem_shift)
            mmu = fermi_e + chem_shift
         enddo
24       continue
         !---------------------------------!
         close (30)
      enddo
      !##################################!

      if (.not. allocated(Omega)) then
         j = 1
         do i = imin, imax
            if (frequ(i) > 0.d0) then
               j = j + 1
            endif
         enddo
         allocate (Omega(j), optics(j, 4), opticstot(j, 4))
         optics = 0.d0
         opticstot = 0.d0
         Omega = 0.d0
         j = 1
         do i = imin, imax
            if (frequ(i) > 0.d0) then
               j = j + 1
               Omega(j) = frequ(i)
            endif
         enddo
         allocate (i_plus_j(1:size(Omega), imin:imax))
      endif

      write (*, *) 'reading data from files done, now move on to calculations'

      i = (imax + imin)/2
      dstep = abs(frequ(i) - frequ(i - 1))

      write (*, *) ' frequency dstep : ', dstep

      i_plus_j = -1
      do i = 1, size(Omega)
         do j = imin, imax
            !w_i + w_j
            frequ_temp = Omega(i) + frequ(j)
            if (frequ_temp < frequ(imax) .and. frequ_temp > frequ(imin)) then
               call findin1Dmesh(frequ(imin:imax), imax - imin + 1, frequ_temp, i_plus_j(i, j))
               if (i_plus_j(i, j) > 0) then
                  i_plus_j(i, j) = i_plus_j(i, j) + imin
                  if (abs(Omega(i) + frequ(j) - frequ(i_plus_j(i, j))) > dstep) then
                     i_plus_j(i, j) = -1
                     write (*, *) 'FREQUENCY NOT FOUND, SHOULD NOT HAPPEN'
                     write (*, *) 'freq min - max      : ', frequ(imin), frequ(imax)
                     write (*, *) 'Omega_i and frequ_j : ', Omega(i), frequ(j)
                     write (*, *) 'i+j                 : ', Omega(i) + frequ(j)
                     write (*, *) 'dstep               : ', dstep
                  endif
               else
                  i_plus_j(i, j) = -1
               endif
            endif
         enddo
      enddo

      write (*, *) 'frequencies .........'
      do i = imin, imax
         write (*, '(a,i5,200f15.4)') '    frequ : ', i, frequ(i)
         fermiv(i) = fermi_dirac_(frequ(i), 0.d0, TT)
         derfermiv(i) = derivative_fermi_dirac_(frequ(i), 0.d0, TT)
         write (*, '(a,i5,2000f15.4)') 'fermiv   : ', i, fermiv(i)
         write (*, '(a,i5,2000f15.4)') 'derFermi : ', i, derfermiv(i)
      enddo
      write (*, *) '.....................'

      if (allocated(vertex__)) deallocate (vertex__)
      allocate (vertex__(1:size(Omega), imin:imax, 4))
      do i = 1, size(Omega)
         do j = imin, imax
            do ttt = 1, 4
               vertex__(i, j, ttt) = vertex_(frequ(j), Omega(i), ttt)
            enddo
         enddo
      enddo

      write (*, *) 'START OPTICS'
      !----------------------------------------------!
      !----------------------------------------------!
      !----------------------------------------------!

      if (kp == 1) opticstot = 0.d0

      optics = 0.d0

      do i = rank + 1, size(Omega), size2
         if (Omega(i) >= -1.d-7) then
            write (*, *) 'DOING FREQU OPTICS : ', i, size(Omega)
            !################ INTEGRAL #################!
            a1 = 0.d0; a2 = 0.d0; a3 = 0.d0; a4 = 0.d0
!$OMP PARALLEL PRIVATE(j,F_factor,rtemp,k1,k2,ttt),REDUCTION(+:a1),REDUCTION(+:a2),REDUCTION(+:a3),REDUCTION(+:a4)
!$OMP DO
            do j = imin, imax
               if (abs(Omega(i)) < small_scale_derivative) then
                  F_factor = derfermiv(j)
               else
                  if (i_plus_j(i, j) <= imax .and. i_plus_j(i, j) >= imin) then
                     F_factor = (fermiv(i_plus_j(i, j)) - fermiv(j))/Omega(i)
                  else
                     F_factor = 0.d0
                  endif
               endif
               if (i_plus_j(i, j) >= imin .and. i_plus_j(i, j) <= imax .and. abs(F_factor) > 1.d-4) then
                  rtemp = 0.d0
                  if (.not. usexy) then
                     do k1 = 1, nn
                        do k2 = 1, nn
                           rtemp = rtemp + green(i_plus_j(i, j), k1, k2)*green(j, k2, k1)
                        enddo
                     enddo
                  else
                     do k1 = 1, nn
                        do k2 = 1, nn
                           rtemp = rtemp + green(i_plus_j(i, j), k1, k2)*greeny(j, k2, k1)
                        enddo
                     enddo
                  endif
                  rtemp = rtemp*F_factor
                  a1 = a1 + vertex__(i, j, 1)*rtemp
                  a2 = a2 + vertex__(i, j, 2)*rtemp
                  a3 = a3 + vertex__(i, j, 3)*rtemp
                  a4 = a4 + vertex__(i, j, 4)*rtemp
               endif
            enddo
!$OMP END DO
!$OMP END PARALLEL
            tttemp(1) = a1
            tttemp(2) = a2
            tttemp(3) = a3
            tttemp(4) = a4
            optics(i, :) = tttemp(:)
            write (*, *) 'OPTICS   : ', Omega(i), optics(i, 1)
            !################ INTEGRAL #################!
         endif
      enddo

      call mpisum_(optics)
      optics = optics*dstep

      opticstot = opticstot + optics*normk

      units_for_optics = pi/Volume
      units_for_optics = units_for_optics/(hartree_in_eV**2)
      units_for_optics = units_for_optics/hbar_divided_e2_in_ohm
      units_for_optics = units_for_optics*(hbar_divided_by_second_in_ev**2)
      units_for_optics = units_for_optics*(hbar_divided_by_me_in_Angstrom_square_divided_by_second**2)
      units_for_optics = units_for_optics*(one_over_Angstrom_in_one_over_cm)
      units_for_optics = units_for_optics*(one_angstrom_in_bohr_radius**2)

      !----------------------------------------------!
      !----------------------------------------------!
      !----------------------------------------------!

      !---------------------!
      if (kp == nkpoints) then

         if (rank == 0 .and. is == 1) then
            do ttt = 1, 4
               open (unit=1001, file='optical_conductivity_pm_vertex_'//trim(adjustl(toString(ttt))))
               write (1001, *) '#volume cell is (in Angstrom^3) : ', Volume
               write (1001, *) '#units factor was               : ', units_for_optics
               do i = 1, size(Omega)
                  write (1001, *) Omega(i)*hartree_in_eV, opticstot(i, ttt)*units_for_optics*2.d0
               enddo
               close (1001)
            enddo
         endif

         if (rank == 0) then
            do ttt = 1, 4
           open (unit=1001, file='optical_conductivity_spin'//trim(adjustl(toString(is)))//'_vertex_'//trim(adjustl(toString(ttt))))
               write (1001, *) '#volume cell is (in Angstrom^3) : ', Volume
               write (1001, *) '#units factor was               : ', units_for_optics
               do i = 1, size(Omega)
                  write (1001, *) Omega(i)*hartree_in_eV, opticstot(i, ttt)*units_for_optics
               enddo
               close (1001)
            enddo
         endif

         if (rank == 0 .and. is == 1) then
            do ttt = 1, 4
               open (unit=1001, file='optical_conductivity_pm_vertex_nounits_'//trim(adjustl(toString(ttt))))
               write (1001, *) '#volume cell is (in Angstrom^3) : ', Volume
               write (1001, *) '#units factor was               : ', units_for_optics
               do i = 1, size(Omega)
                  write (1001, *) Omega(i)*hartree_in_eV, opticstot(i, ttt)*2.d0
               enddo
               close (1001)
            enddo
         endif

         if (rank == 0) then
            do ttt = 1, 4
   open (unit=1001, file='optical_conductivity_spin'//trim(adjustl(toString(is)))//'_vertex_nounits_'//trim(adjustl(toString(ttt))))
               write (1001, *) '#volume cell is (in Angstrom^3) : ', Volume
               write (1001, *) '#units factor was               : ', units_for_optics
               do i = 1, size(Omega)
                  write (1001, *) Omega(i)*hartree_in_eV, opticstot(i, ttt)
               enddo
               close (1001)
            enddo
         endif

      endif
      !---------------------!

   enddo
   enddo

!###########################################################!
!###########################################################!
!###########################################################!
!###########################################################!
!###########################################################!

   !-------------------------------!
   call finalize_my_simulation
   !-------------------------------!

contains

!----------------------!
!----------------------!
!----------------------!
!----------------------!
!----------------------!
!----------------------!
!----------------------!
!----------------------!
!----------------------!
!----------------------!

   subroutine check_file
      INQUIRE (file=trim(adjustl(filename)), exist=check)
      if (.not. check) then
         write (*, *) 'optics file is not present, filename = ', filename
         stop
      endif
   end subroutine

   character(40) function set_name(is, k, mpi_rank)
      implicit none
      integer :: is, mpi_rank, k
 set_name='optics_data_spin_'//TRIM(ADJUSTL(toString(is)))//'_k_'//TRIM(ADJUSTL(toString(k)))//'_rank_'//TRIM(ADJUSTL(toString(mpi_rank)))
   end function

   real(kind=DP) function vertex_(w, omega, ttt)
      implicit none
      real(kind=DP) :: w, omega
      integer :: ttt
      SELECT CASE (ttt)
      CASE (1)
         vertex_ = 1.d0
      CASE (2)
         vertex_ = (w + omega/2.d0)**2
      CASE (3)
         vertex_ = (w + omega/2.d0)
      CASE (4)
         vertex_ = 1.d0
      CASE DEFAULT
         stop 'error vertex_ in polarization bubble not defined'
      END SELECT
   end function

!----------------------!
!----------------------!
!----------------------!
!----------------------!
!----------------------!
!----------------------!
!----------------------!
!----------------------!
!----------------------!
!----------------------!

end program
