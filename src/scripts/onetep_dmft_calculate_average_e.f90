program averageenergy
   implicit none

   real(8) :: energies(1000), zpart, ee, tt
   integer :: i, j, countit

   call system("grep -a 'E0 = ' logfile-p1 | wc -l > values")
   call system("grep -a 'E0 = ' logfile-p1  | awk '{print $4}'>> values")
   call system(" cat *.dat |grep temp | awk '{print $3}' > temp ")

   open (unit=30, file='temp')
   read (30, *) tt
   close (30)
   write (*, *) 'temperature is : ', tt

   open (unit=30, file='values')
   read (30, *) countit
   do i = 1, countit
      read (30, *) energies(i)
   enddo

   close (30)
   write (*, *) 'ENERGIES ARE : ', energies(1:countit)

   zpart = 0.
   ee = 0.

   do i = 1, countit
      zpart = zpart + exp(-(energies(i) - minval(energies(1:countit)))/tt)
      ee = ee + exp(-(energies(i) - minval(energies(1:countit)))/tt)*energies(i)
   enddo
   ee = ee/zpart
   write (*, *) 'partition function is : ', zpart
   write (*, *) 'energy is             : ', ee

end program
