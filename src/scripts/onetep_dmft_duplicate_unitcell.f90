program duplicate

   use genvar, only: DP

   implicit none
   integer       :: i, j, k, l, m
   integer       :: ii(3), iii(3), k1, k2, k3
   real(kind=DP)       :: pos_shift, shift(3), a(3, 3), aa(3, 3), v(3), coord(2000, 3)
   character(30) :: filename, atom(2000), temp, species

   open (unit=1010, file='duplicate.input')
   write (*, *) 'please enter file name'
   read (1010, *) filename
   open (unit=10, file=filename)
   write (*, *) 'new unitcell bounds x (min-max integer): '
   read (1010, *) ii(1), iii(1)
   write (*, *) 'new unitcell bounds y (min-max integer): '
   read (1010, *) ii(2), iii(2)
   write (*, *) 'new unitcell bounds z (min-max integer): '
   read (1010, *) ii(3), iii(3)
   write (*, *) 'atom type to duplicate (string) (ALL for all)'
   read (1010, *) species
   write (*, *) 'shift for positive carthesian coordinates (one real number) : '
   read (1010, *) pos_shift
   write (*, *) 'first 3 lines of the file should be the original unitcell vector'
   close (1010)

   do i = 1, 3
      read (10, *) (a(i, j), j=1, 3)
      write (*, *) 'original unitcell vector : ', i
      write (*, *) a(i, :)
   enddo

   do i = 1, 3
      aa(i, :) = dble(iii(i) - ii(i) + 1)*a(i, :)
   enddo

   k = 0
   do
      read (10, *, end=66) temp, v(1), v(2), v(3)
      k = k + 1
      atom(k) = temp
      write (*, *) 'atom  :  ', atom(k)
      write (*, *) 'coord :  ', v
      coord(k, :) = v
   enddo
66 continue
   write (*, *) 'HOW MANY ATOMS IN UNITCELL : ', k

   open (unit=20, file='unitcell_output')
   do i = 1, 3
      write (20, '(3f14.5)') (aa(i, j), j=1, 3)
   enddo
   do k1 = ii(1), iii(1)
   do k2 = ii(2), iii(2)
   do k3 = ii(3), iii(3)
      write (*, *) 'cell : ', k1, k2, k3
      shift = k1*a(1, :) + k2*a(2, :) + k3*a(3, :)
      do i = 1, k
         write (*, *) 'atom : ', i
         if (species == 'ALL' .or. atom(i) == species) then
            write (20, '(a4,3f16.5)') TRIM(ADJUSTL(atom(i))), (coord(i, j) + shift(j) + pos_shift, j=1, 3)
         endif
      enddo
   enddo
   enddo
   enddo
   close (20)

!example file
!    4.8376620217     0.0000000000    -3.0950037256
!    0.0000000000     4.5169982501     0.0000000000
!    0.0000000000     0.0000000000     5.3749979002
! V       1.1707142090    -2.1455741688    -3.3021149042
! V      -1.1707142090     0.1129249563     0.6146159541
! V      -1.1707142090     2.1455741688     3.3021149042
! V       1.1707142090    -0.1129249563    -0.6146159541
! O       0.4837662023     1.3099294925    -1.9219997428
! O      -0.4837662023    -0.9485696325    -0.7654992073
! O      -0.4837662023    -1.3099294925     1.9219997428
! O       0.4837662023     0.9485696325     0.7654992073
! O       1.8866881886    -0.8582296675    -2.3358010117
! O      -1.8866881886     1.4002694575    -0.3516979384
! O      -1.8866881886     0.8582296675     2.3358010117
! O       1.8866881886    -1.4002694575     0.3516979384

   close (10)
end program
