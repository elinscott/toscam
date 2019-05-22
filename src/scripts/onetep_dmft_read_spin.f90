program correctspin

   implicit none

   real(8) :: mat(5, 5, 3), matb(5, 5, 3), a1, a2, a3, sz, mag
   integer :: i, j

   call system(" cat logfile-p1 |grep -a -A 17 'SPIN CORRELATIONS' |tail -15 > spin_temp ")

   open (unit=10, file='spin_temp')

   do
      read (10, *, end=67) i, j, a1, a2, a3
      write (*, '(a,2i3,3f10.3)') 'i,j,S+S-/2,S-S+/2,Sz', i, j, a1, a2, a3
      mat(i, j, 1) = a1
      if (abs(a2) > 1.d-5) then
         mat(i, j, 2) = a2
      else
         mat(i, j, 2) = a1
      endif
      mat(i, j, 3) = a3
      mat(j, i, :) = mat(i, j, :)
   enddo
67 continue
   close (10)

   call system(" cat logfile-p1 | grep -a '# <Stot>' | awk '{ print $4 }' > spin_temp ")
   open (unit=10, file='spin_temp')
   read (10, *) sz
   close (10)

   call system("cat logfile-p1 |grep -a  'mag' | awk '{print $4+$5+$6+$7+$8}' > spin_temp")
   open (unit=10, file='spin_temp')
   read (10, *) mag
   close (10)

   write (*, *) 'TOTAL MAGNETIZATION (NUP-NDN) : ', mag
   write (*, *) 'TOTAL SPIN S                  : ', sz
   write (*, *) 'TOTAL SPIN S=sqrt(S^2)        : ', sqrt(abs(sum(mat) - sz**2))

   matb = mat
   do i = 1, 5
      do j = 1, 5
         if (i /= j) mat(i, j, :) = 0.
      enddo
   enddo
   write (*, *) 'TOTAL SPIN S=sqrt(sum S_i^2 ) : ', sqrt(sum(mat))
   mat = matb

   mat(:, :, 1:2) = 0.
   write (*, *) 'TOTAL SPIN Sz=sqrt(Sz^2) : ', sqrt(sum(mat))

   do i = 1, 5
      do j = 1, 5
         if (i /= j) mat(i, j, :) = 0.
      enddo
   enddo
   write (*, *) 'TOTAL SPIN Sz=sqrt(Sum Sz_i^2) : ', sqrt(sum(mat))

   mat = sqrt(abs(mat))

   write (*, *) 'TOTAL SPIN Sz=Sum sqrt(Sz_i^2) : ', sum(mat)

   call system("rm spin_temp")

end program

