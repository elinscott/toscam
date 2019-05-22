program collect
   use strings, only: replace_in_string, string, assignment(=)
   implicit none
   integer                :: i, j, ien, pub_dmft_points, channels, channelsb
   real(8)                :: mmu, ttemp
   integer                :: num, len__, status__, LL, paramagnetic
   character(200)         :: value(2000)
   character(200)         :: filename
   complex(8)             :: frequ
   complex(8), allocatable :: sigma_embed_(:, :, :), sigma_embed_temp(:, :), frequ_(:)
   type(string)           :: cc_

   num = COMMAND_ARGUMENT_COUNT()
   if (num > size(value)) then
      write (*, *) 'ERROR size value too small'
      stop
   endif
   do i = 1, num
      call GET_COMMAND_ARGUMENT(i, value(i), len__, status__)
      write (*, *) 'MY ARGUMENTS', value(i)
   enddo

   if (num == 1) then
      call system(" mv "//TRIM(ADJUSTL(value(1)))//" `echo '"//TRIM(ADJUSTL(value(1)))//"' | sed 's/_rank1//g'`  ")
      stop
   endif

   do i = 1, num
      open (unit=10001, file=TRIM(ADJUSTL(value(i))), form='unformatted')
      read (10001, end=66) ien, pub_dmft_points, ttemp, channels, channelsb
      rewind (10001)
      write (*, *) 'ien,npoints,ttemp,channels: ', ien, pub_dmft_points, ttemp, channels, channelsb
      if (.not. allocated(sigma_embed_)) then
         allocate (frequ_(pub_dmft_points), sigma_embed_(pub_dmft_points, channels, channelsb), &
                   & sigma_embed_temp(channels, channelsb))
         sigma_embed_ = 0.d0
      endif
      write (*, *) 'start loop'
      j = 0
      do
         read (10001, end=66) ien, pub_dmft_points, ttemp, channels, channelsb
         read (10001, end=66) mmu, frequ, sigma_embed_temp
         j = j + 1
         sigma_embed_(ien, :, :) = sigma_embed_(ien, :, :) + sigma_embed_temp
         frequ_(ien) = frequ
      enddo
66    continue
      write (*, *) 'number of frequencies : ', j
      write (*, *) 'file/tot : ', i, num
      close (10001)
   enddo

   do i = 1, num
      call system(" rm "//TRIM(ADJUSTL(value(i))))
   enddo

   write (*, *) 'TEMP = ', ttemp
   write (*, *) 'TARGET FILE'
   filename = TRIM(ADJUSTL(value(1)))
   do i = LEN(filename), 1, -1
      if (filename(i:i) == "k") then
         exit
      else
         filename(i:i) = " "
      endif
   enddo
   cc_ = TRIM(ADJUSTL(filename))
   call replace_in_string(cc_, '_rank', '', 'last')
   filename = cc_

   write (*, *) 'TARGET GREEN FUNCTION FILE = ', TRIM(ADJUSTL(filename))
   write (*, *) 'writing file      : ', mmu
   write (*, *) 'pub_dmft_points   : ', pub_dmft_points
   write (*, *) 'ttemp             : ', ttemp
   write (*, *) 'channels          : ', channels
   write (*, *) 'channelsb         : ', channelsb
   call system("rm "//TRIM(ADJUSTL(filename)))
   write (*, *) 'opening file'
   open (unit=2000, file=TRIM(ADJUSTL(filename)), form='unformatted')
   write (*, *) 'start loop'
   do j = 1, pub_dmft_points
      write (2000) j, pub_dmft_points, ttemp, channels, channels
      write (2000) mmu, frequ_(j), sigma_embed_(j, :, :)
   enddo
   close (2000)

   stop
end program
