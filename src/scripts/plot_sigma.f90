program plotsigma
   implicit none
   integer                :: chan, i, j
   complex(kind=DP), allocatable :: mat1(:, :)

   write (*, *) 'ENTER CHANNELS'
   read (*, *) chan

   allocate (mat1(chan, chan))

   open (unit=10, file='full', form='unformatted')
   open (unit=11, file='diag')

   i = 0
   do
      i = i + 1
      write (*, *) 'frequ : ', i
      read (10, end=60) mat1
      write (11, '(i6,20f14.7)') i, (real(mat1(j, j)), aimag(mat1(j, j)), j=1, chan)
   enddo
60 continue
   close (10); close (11)
   write (*, *) 'total frequ : ', i - 1
end program
