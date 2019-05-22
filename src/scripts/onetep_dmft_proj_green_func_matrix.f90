program projgreen_software
   use linalg
   use StringManip, only: StrInt2
   use matrix, only: diag
   use mesh
   implicit none
   integer                   :: shapeit(4), i, j, k, l, m, i1, i2, i3
   real(8)                   :: ii, mat(3, 3)
   complex(8), allocatable    :: projgreen(:, :, :, :)
   real(8), allocatable       :: positions(:, :)
   logical                   :: check
   character(300)            :: value(100)
   character(22), allocatable :: label(:)
   integer                   :: nat, num, len__, status__

   num = COMMAND_ARGUMENT_COUNT()
   if (num /= 1) then
      write (*, *) 'wrong number of arguments'
      write (*, *) 'NUM       : ', num
      write (*, *) 'ARGUMENTS : ', value(1:num)
      write (*, *) 'execting 1 argument (name of file containg proj green)'
      stop
   endif
   do i = 1, num
      call GET_COMMAND_ARGUMENT(i, value(i), len__, status__)
   enddo
   inquire (file=trim(adjustl(value(1))), exist=check)
   if (.not. check) then
      write (*, *) 'file does not exist'
      stop
   endif

   open (unit=4, file=trim(adjustl(value(1))), form='unformatted')

   read (4) shapeit
   allocate (projgreen(shapeit(1), shapeit(2), shapeit(3), shapeit(4)))
   read (4) projgreen
   read (4) nat
   allocate (positions(nat, 3))
   read (4) (mat(j, 1), j=1, 3)
   read (4) (mat(j, 2), j=1, 3)
   read (4) (mat(j, 3), j=1, 3)
   do j = 1, nat
      read (4) positions(j, 1:3)
      write (*, '(a,i5,3f10.3)') 'position in unitcell : ', j, (positions(j, i), i=1, 3)
   enddo
   close (4)

   do i = 1, shapeit(3)
      do j = 1, nat
  write (300 + i, *) manhattan_distance(positions(1, :), positions(j, :)), real(projgreen(1, j, i, i)), aimag(projgreen(1, j, i, i))
         write (400 + i, *) manhattan_distance(positions(1, :), positions(j, :)), abs(real(projgreen(1, j, i, i)))
         write (500 + i, *) manhattan_distance(positions(1, :), positions(j, :)), abs(aimag(projgreen(1, j, i, i)))
         write (600 + i, *) manhattan_distance(positions(1, :), positions(j, :)), abs(projgreen(1, j, i, i))
         write (700 + i, *) manhattan_distance(positions(1, :), positions(j, :)), sum(diag(abs(aimag(projgreen(1, j, :, :)))))
      enddo
   enddo

   do ii = 1, size(positions, 1)
   do i = 1, shapeit(3)
      do j = 1, nat
      write(3000+i,*) manhattan_distance(positions(ii,:),positions(j,:)),real(projgreen(ii,j,i,i)),aimag(projgreen(ii,j,i,i))
         write (4000 + i, *) manhattan_distance(positions(ii, :), positions(j, :)), abs(real(projgreen(ii, j, i, i)))
         write (5000 + i, *) manhattan_distance(positions(ii, :), positions(j, :)), abs(aimag(projgreen(ii, j, i, i)))
         write (6000 + i, *) manhattan_distance(positions(ii, :), positions(j, :)), abs(projgreen(ii, j, i, i))
         write (7000 + i, *) manhattan_distance(positions(ii, :), positions(j, :)), sum(diag(abs(aimag(projgreen(ii, j, :, :)))))
      enddo
   enddo
   enddo

contains

   real(8) function manhattan_distance(v1, v2)
      implicit none
      integer :: i1_, i2_, i3_
      real(8) :: v1(3), v2(3), dd, temp, v1_(3)
      dd = 1.d20
      do i1_ = -1, 1
      do i2_ = -1, 1
      do i3_ = -1, 1
         v1_ = dble(i1_)*mat(:, 1) + dble(i2_)*mat(:, 2) + dble(i3_)*mat(:, 3)
         temp = norm(v1 + v1_ - v2)
         if (temp < dd) dd = temp
      enddo
      enddo
      enddo
      manhattan_distance = dd
   end function

end program
