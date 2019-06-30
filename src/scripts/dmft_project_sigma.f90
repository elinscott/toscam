program projectorsigma

   use strings, only: replace_in_string, string, assignment(=)
   use matrix
   use StringManip, only: toString
   use StringManip, only: StrInt2

   implicit none

   integer                :: i, j, k, l, ien, pub_dmft_points, channels, atom, channelsb, atomb
   real(8)                :: mmu, ttemp
   integer                :: num, len__, status__, paramagnetic
   character(200)         :: value(2000)
   character(200)         :: filename
   complex(8)             :: frequ
   complex(8), allocatable :: green_2(:, :), green_(:, :, :, :), green_temp(:, :), frequ_(:)
   type(string)           :: cc_
   logical, allocatable    :: mask_proj(:)
   integer                :: nproj
   integer               :: k_, l_, iii, kk, ll, nnn, kkk, lll
   logical                :: check, no_dimer
   integer                :: ndimer, ii, atom_number

   num = COMMAND_ARGUMENT_COUNT()
   if (num > size(value)) then
      write (*, *) 'ERROR size value too small'
      stop
   endif
   do i = 1, num
      call GET_COMMAND_ARGUMENT(i, value(i), len__, status__)
      write (*, *) 'MY ARGUMENTS', value(i)
   enddo

   if (num /= 4) then
      no_dimer = .true.
      ndimer = 1
   else
      ndimer = 2
      no_dimer = .false.
      INQUIRE (file='mask_dimer', EXIST=check)
      if (.not. check) then
         write (*, *) 'doing single site DMFT'
         stop
      endif
   endif

   INQUIRE (file='compound_tot_orb', EXIST=check)
   if (.not. check) then
      write (*, *) 'ERROR compounds_tot file does not exist'
      stop
   endif
   open (unit=444, file='compound_tot_orb')

   nnn = 0
   call get_atom_number(atom_number, value(1))
   write (*, *) 'THIS IS ATOM : ', atom_number
   do ii = 1, atom_number
      read (444, *, end=32) nnn
   enddo
32 continue

   write (*, *) 'COMPOUNDS WITH [x] total elements : ', nnn
   close (444)

   channels = nnn; channelsb = nnn

   do i = 1, num
      open (unit=10001, file=TRIM(ADJUSTL(value(i))), form='unformatted')
      if (.not. allocated(green_temp)) allocate (green_temp(channels, channelsb))
      j = 0
      do
         read (10001, end=65) green_temp
         j = j + 1
      enddo
65    continue
      rewind (10001)
      write (*, *) 'projecting SIGMA, [x] frequencies : ', j
      write (*, *) 'name file : ', value(1)
      if (j == 0) then
         write (*, *) 'ERROR projecting sigma 0 frequ'
         stop
      endif
      pub_dmft_points = j
      if (.not. allocated(green_)) allocate (green_(pub_dmft_points, channels, channelsb, num))
      j = 0
      do
         read (10001, end=66) green_temp
         j = j + 1
         green_(j, :, :, i) = green_temp
      enddo
66    continue
      write (*, *) 'number of frequencies : ', j
      write (*, *) 'file/tot : ', i, num
      close (10001)
   enddo

   do i = 1, num
      call system("rm "//TRIM(ADJUSTL(value(i))))
   enddo

   filename = TRIM(ADJUSTL(value(1)))
   write (*, *) 'TARGET COMBINED PROJECTED GREEN FUNCTION FILE = ', TRIM(ADJUSTL(filename))
   write (*, *) 'pub_dmft_points   : ', pub_dmft_points
   write (*, *) 'channels          : ', channels
   write (*, *) 'channelsb         : ', channelsb

   open(unit=1010,file='mask_projections')
   write(*,*) 'reading  mask'
   allocate(mask_proj(channels*ndimer))

   do j=1,atom_number
      read(1010,*)  (mask_proj(i),i=1,channels*ndimer)
   enddo

   write(*,*)  (mask_proj(i),i=1,channels*ndimer)
   nproj=0
   do i=1,channels*ndimer
      if(mask_proj(i)) nproj=nproj+1
   enddo
   nproj=nproj/ndimer
   write(*,*) 'THERE ARE [x] KEPT FUNCTIONS : ', nproj
   close(1010)

   allocate (green_2(nproj*ndimer, nproj*ndimer))

   open (unit=2000, file=TRIM(ADJUSTL(filename)), form='unformatted')
   do j = 1, pub_dmft_points
      ! green_func1_1_1 green_func2_1_1 green_func_1_1_1_dimer green_func_2_1_1_dimer
      green_2 = 0.; 
      !===============================!
      !===============================!
      !===============================!
      do iii = 1, ndimer**2
         !-------------------------------!
         select case (iii)
         case (1)
            kk = 0; ll = 0; kkk = 0; lll = 0
         case (2)
            kk = nproj; ll = nproj; kkk = channels; lll = channels
         case (3)
            kk = 0; ll = nproj; kkk = 0; lll = channels
         case (4)
            kk = nproj; ll = 0; kkk = channels; lll = 0
         end select
         k_ = 0
         do k = 1, channels
            if (mask_proj(k + kkk)) then
               k_ = k_ + 1
               l_ = 0
               do l = 1, channels
                  if (mask_proj(l + lll)) then
                     l_ = l_ + 1
                     green_2(k_ + kk, l_ + ll) = green_(j, k, l, iii)
                  endif
               enddo
            endif
         enddo
         !-------------------------------!
      enddo
      !===============================!
      !===============================!
      !===============================!

      write (*, *) '--------------------------------------------------------'
      call write_array(green_2, ' SIGMA PROJECTED ', short=.true., unit=6)
      write (*, *) '--------------------------------------------------------'
      write (2000) green_2
   enddo
   close (2000)

contains

!===================================================================================================!

   subroutine get_atom_number(newatom, filename)
      implicit none
      integer       :: newatom, j, i
      character*(*) :: filename
      do i = LEN(filename), 1, -1
         if (filename(i:i) == "t") then
            exit
         endif
      enddo
      i = i + 1
      do j = i, 2000
         if (filename(j:j) == "_") then
            exit
         endif
      enddo
      j = j - 1
      write (*, *) 'ATOM NUMBER : ', filename(i:j)
      newatom = StrInt2(TRIM(ADJUSTL(filename(i:j))))
   end subroutine
!===================================================================================================!

end program
