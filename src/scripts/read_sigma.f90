program readsigma
   integer :: chan, i
   complex(kind=DP), allocatable :: mat1(:, :), mat2(:, :)
   write (*, *) 'ENTER CHANNELS'
   read (*, *) chan
   allocate (mat1(chan, chan), mat2(chan, chan))
   open (unit=10, file='v1', form='unformatted')
   open (unit=11, file='v2', form='unformatted')
   i = 0
   do
      i = i + 1
      read (10, end=60) mat1
      read (11, end=60) mat2
      if (maxval(abs(mat1 - mat2)) > 1.d-8) then
         write (*, *) 'DIFF RE-IM: ', i, maxval(abs(real(mat1 - mat2))), maxval(abs(aimag(mat1 - mat2)))
      endif
   enddo
60 continue
   close (10); close (11)
   write (*, *) 'total frequ : ', i - 1
end program
