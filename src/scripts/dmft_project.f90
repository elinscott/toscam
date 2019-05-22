program projector
   use strings, only: replace_in_string, string, assignment(=)
   use matrix
   use StringManip, only: toString
   implicit none
   integer                :: i, j, k, l, ien, pub_dmft_points, channels, atom, channelsb, atomb
   real(8)                :: mmu, ttemp
   integer                :: num, len__, status__, paramagnetic
   character(200)         :: value(2000)
   character(200)         :: filename
   complex(8)             :: frequ
   complex(8), allocatable :: green_2(:, :), green_(:, :, :, :), green_temp(:, :), frequ_(:)
   real(8), allocatable    :: occupation_matrix(:, :, :, :)
   real(8), allocatable    :: occup_2(:, :, :)
   type(string)           :: cc_
   logical, allocatable    :: mask_proj(:)
   integer                :: nproj
   integer                :: k_, l_, iii, kk, ll
   logical                :: check, no_dimer
   integer                :: ndimer

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
      no_dimer = .false.
      ndimer = 2
      INQUIRE (file='mask_dimer', EXIST=check)
      if (.not. check) then
         write (*, *) 'doing single site DMFT'
         stop
      endif
   endif

   do i = 1, num
      open (unit=10001, file=TRIM(ADJUSTL(value(i))), form='unformatted')
      read (10001, end=66) ien, pub_dmft_points, ttemp, channels, channelsb, atom, atomb
      rewind (10001)
      write (*, *) 'ien,npoints,ttemp,channels,atom : ', ien, pub_dmft_points, ttemp, channels, channelsb, atom, atomb
      if (.not. allocated(green_)) allocate (frequ_(pub_dmft_points), green_(pub_dmft_points, channels, channelsb, num), &
                                        & green_temp(channels, channelsb), occupation_matrix(channels, channels, num, 2))
      green_(:, :, :, i) = 0.d0
      j = 0
      do
         read (10001, end=66) ien, pub_dmft_points, ttemp, channels, channelsb, atom, atomb
         read (10001, end=66) occupation_matrix(:, :, i, 1), occupation_matrix(:, :, i, 2)
         read (10001, end=66) mmu, frequ, green_temp
         green_(ien, :, :, i) = green_(ien, :, :, i) + green_temp
         frequ_(ien) = frequ
      enddo
66    continue
      write (*, *) 'number of frequencies in file : ', j
      write (*, *) 'file/tot                      : ', i, num
      close (10001)
   enddo
   write (*, *) 'TEMP = ', ttemp

   do i = 1, num
      call system("rm "//TRIM(ADJUSTL(value(i))))
   enddo

   filename = TRIM(ADJUSTL(value(1)))
   write (*, *) 'TARGET COMBINED PROJECTED GREEN FUNCTION FILE = ', TRIM(ADJUSTL(filename))
   write (*, *) 'writing file      : ', mmu
   write (*, *) 'pub_dmft_points   : ', pub_dmft_points
   write (*, *) 'ttemp             : ', ttemp
   write (*, *) 'channels          : ', channels
   write (*, *) 'channelsb         : ', channelsb
   write (*, *) 'atom              : ', atom
   write (*, *) 'atomb             : ', atomb

   open (unit=1010, file='mask_projections')
   write (*, *) 'reading  mask'
   allocate (mask_proj(channels))
   read (1010, *) (mask_proj(i), i=1, channels)
   write (*, *) (mask_proj(i), i=1, channels)
   nproj = 0
   do i = 1, channels
      if (mask_proj(i)) nproj = nproj + 1
   enddo
   write (*, *) 'THERE ARE [x] KEPT FUNCTIONS : ', nproj
   close (1010)

   allocate (green_2(nproj*ndimer, nproj*ndimer), occup_2(nproj*ndimer, nproj*ndimer, 2))

   open (unit=2000, file=TRIM(ADJUSTL(filename)), form='unformatted')
   do j = 1, pub_dmft_points
      ! green_func1_1_1 green_func2_1_1 green_func_1_1_1_dimer green_func_2_1_1_dimer
      green_2 = 0.; occup_2 = 0.

      !===============================!
      !===============================!
      !===============================!
      do iii = 1, ndimer**2
         !-------------------------------!
         select case (iii)
         case (1)
            kk = 0; ll = 0
         case (2)
            kk = nproj; ll = nproj; 
         case (3)
            kk = 0; ll = nproj
         case (4)
            kk = nproj; ll = 0
         end select
         k_ = 0
         do k = 1, channels
            if (mask_proj(k)) then
               k_ = k_ + 1
               write (*, *) 'k,k_=', k_
               l_ = 0
               do l = 1, channels
                  if (mask_proj(l)) then
                     l_ = l_ + 1
                     write (*, *) 'l,l_=', l, l_
                     if (iii == 1 .or. iii == 2) occup_2(k_ + kk, l_ + ll, :) = occupation_matrix(k, l, iii, :)
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

      write (*, *) 'MASK : ', mask_proj(:), nproj, channels
      write (*, *) '--------------------------------------------------------'
      call write_array(occupation_matrix(:, :, 1, 1), 'OCCUPATION UP NOT PROJECTED(1)', short=.true., unit=6)
      call write_array(occupation_matrix(:, :, 2, 1), 'OCCUPATION UP NOT PROJECTED(2)', short=.true., unit=6)
      write (*, *) '--------------------------------------------------------'
      call write_array(occup_2(:, :, 1), ' OCCUPATION UP PROJECTED ', short=.true., unit=6)
      write (*, *) '--------------------------------------------------------'
      call write_array(green_2, ' GREEN PROJECTED ', short=.true., unit=6)
      write (*, *) '--------------------------------------------------------'
      write (2000) j, pub_dmft_points, ttemp, ndimer*nproj, ndimer*nproj, atom, atom
      write (2000) occup_2(:, :, 1), occup_2(:, :, 2)
      write (2000) mmu, frequ_(j), green_2
   enddo
   close (2000)

end program
