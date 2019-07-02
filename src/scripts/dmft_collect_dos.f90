program collect
   use genvar, only: DP
   use strings, only: replace_in_string, string, assignment(=)
   implicit none
   integer                :: i, j, jj, ien, pub_dmft_points, channels, atom, channelsb, atomb
   real(kind=DP)                :: mmu
   integer                :: num, len__, status__, LL, paramagnetic
   character(200)         :: value(2000)
   character(200)         :: filename
   real(kind=DP)                :: frequ
   real(kind=DP), allocatable    :: green_(:, :), green_temp(:), frequ_(:)
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

   do i = 1, num
      open (unit=10001, file=TRIM(ADJUSTL(value(i))))
      read (10001, *, end=66) ien, channels, pub_dmft_points
      rewind (10001)
      write (*, *) 'ien,channels,npoints : ', ien, channels, pub_dmft_points
      if (.not. allocated(green_)) then
         allocate (frequ_(pub_dmft_points), green_(pub_dmft_points, channels), green_temp(channels))
         green_ = 0.d0
         frequ_ = -123456789.d0
      endif
      j = 0
      do
         read (10001, *, end=66) ien, channels, pub_dmft_points, frequ, (green_temp(jj), jj=1, channels)
         j = j + 1
         green_(ien, :) = green_(ien, :) + green_temp
         frequ_(ien) = frequ
      enddo
66    continue
      write (*, *) 'number of frequencies : ', j
      write (*, *) 'file/tot              : ', i, num
      close (10001)
   enddo

   do i = 1, num
      call system(" rm "//TRIM(ADJUSTL(value(i))))
   enddo

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

   call system("rm "//TRIM(ADJUSTL(filename)))

   write (*, *) 'opening file'
   open (unit=2000, file=TRIM(ADJUSTL(filename)))
   do j = 1, pub_dmft_points
      if (abs(frequ_(j) - (-123456789.d0)) > 1.d-5) then
         write (2000, '(i5,i5,i5,300f14.4)') j, channels, pub_dmft_points, frequ_(j), (green_(j, jj), jj=1, channels)
      endif
   enddo
   close (2000)

   stop
end program

