program correctn
   use matrix
   use genvar, only: DP
   implicit none

   real(kind=DP) :: nn, mat(5, 5, 3), matb(5, 5, 3), a1, a2, a3
   integer :: i, j

   call system(" cat logfile-p1 |grep -a -A 17 'DENSITY CORRELATIONS' |tail -15 > n_var_temp ")

   open (unit=10, file='n_var_temp')

   do
      read (10, *, end=67) i, j, a1, a2, a3
      write (*, '(a,2i3,3f10.3)') 'i,j, <n[i]*n[j]>/  <n[i]>*<n[j]> / correlator', i, j, a1, a2, a3
      mat(i, j, 1) = a1
      mat(i, j, 2) = a2
      mat(i, j, 3) = a3
      mat(j, i, :) = mat(i, j, :)
   enddo
67 continue
   close (10)

   nn = sum(sqrt(diag(mat(:, :, 2))))
   write (*, *) 'TOTAL N           : ', nn
   write (*, *) 'Variance of N_tot : ', sqrt(sum(mat(:, :, 1)) - nn**2)
   write (*, *) 'from local corr   : ', sqrt(sum(mat(:, :, 3)))
   call system("rm n_var_temp")

end program

