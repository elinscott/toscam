program local_dimers
   use genvar, only: DP
   use linalg
   use sorting
   use random
   implicit none

   integer           :: i, j, k, l, m
   integer           :: k1, k2, k3

   real(kind=DP)           :: ladder(1000, 0:1000), v(3), v_(3), coord(1000, 3), coord_disorder(1000, 3), proj(1000, 1000)
   character(30)     :: ref, filename, temp, species
   real(kind=DP)           :: www(3), mat(3, 3), ex(3), ey(3), ez(3), normv(1000), dist, cutoff, vvv(3), rutile(3)
   logical           :: same_line(1000, 1000)
   integer           :: jj, i1_, i2_, i3_, i1, i2, kk, order(1000, 0:1000), nneigh(1000, 0:1000), counti(1000)
   real(kind=DP)           :: v1_(3), v12(3), disorder_amp, peierls_move, shift_amp(3)
   logical, parameter :: testingloc = .false.

   disorder_amp = 0.d0
   peierls_move = 0.d0
   shift_amp = 0.d0

   if (.not. testingloc) open (unit=100, file='dimer.input')
   write (*, *) 'please enter file name'
   if (.not. testingloc) read (100, *) filename
   write (*, *) 'center for local coordinates (symbol, example V,Cu,O...)'
   if (.not. testingloc) read (100, *) ref
   if (testingloc) ref = 'V'
   write (*, *) 'your local rotations are at : ', trim(ref)
   write (*, *) 'cutoff (0.05): '
   if (.not. testingloc) read (100, *) cutoff
   write (*, *) 'please enter rutile c axis (3 real numbers) : '
   if (.not. testingloc) read (100, *) rutile
   write (*, *) 'please enter disorder amplitude (linear): '
   if (.not. testingloc) read (100, *) disorder_amp
   write (*, *) 'please enter peierls move (linear): '
   if (.not. testingloc) read (100, *) peierls_move
   write (*, *) 'please enter shift amplitude (linear): '
   if (.not. testingloc) read (100, *) shift_amp
   if (.not. testingloc) close (100)

   if (testingloc) then
      rutile(1) = 5.743
      rutile(2) = 0.
      rutile(3) = 0.
      cutoff = 0.2
      filename = 'test'
   endif

   open (unit=10, file=filename)
   read (10, *) (mat(jj, 1), jj=1, 3)
   read (10, *) (mat(jj, 2), jj=1, 3)
   read (10, *) (mat(jj, 3), jj=1, 3)

   write (*, *) '======================================'
   write (*, '(a,3f10.3)') 'VECTOR A : ', (mat(jj, 1), jj=1, 3)
   write (*, '(a,3f10.3)') 'VECTOR B : ', (mat(jj, 2), jj=1, 3)
   write (*, '(a,3f10.3)') 'VECTOR C : ', (mat(jj, 3), jj=1, 3)
   write (*, *) '======================================'

   open (unit=1010, file='disorder_out')
   open (unit=1111, file='lattice_used_for_dimers')
   open (unit=1112, file='lattice_used_for_disorder')
   j = 0
   do
      read (10, *, end=66) temp, v(1), v(2), v(3)
      v_(1) = random_float_from_gaussian()
      v_(2) = random_float_from_gaussian()
      v_(3) = random_float_from_gaussian()
      v_ = v_*disorder_amp
      v_ = v + v_
      write (1010, '(a5,3f15.6)') temp, v_(1), v_(2), v_(3)
      if (TRIM(temp) == TRIM(ref)) then
         j = j + 1
         coord_disorder(j, :) = v_
         coord(j, :) = v
         write (*, '(a,a5,3f15.6)') 'atom ref: ', TRIM(temp), coord(j, :)
         write (1111, '(a5,3f15.6)') TRIM(temp), coord(j, :)
         write (1112, '(a5,3f15.6)') TRIM(temp), coord_disorder(j, :)
      endif
   enddo
66 continue
   close (10)
   close (1010)
   close (1111)
   close (1112)

   write (*, *) 'They are [x] atoms to form dimers/pairs : ', j
   if (j > size(counti)) then
      write (*, *) 'arrays on stack are too small'
      stop
   endif

   same_line = .false.
   nneigh = 0
   counti = 0
   ladder = 0.

   if (norm(vecprod(mat(:, 1), rutile))/norm(rutile) < cutoff) mat(:, 1) = 0.
   if (norm(vecprod(mat(:, 2), rutile))/norm(rutile) < cutoff) mat(:, 2) = 0.
   if (norm(vecprod(mat(:, 3), rutile))/norm(rutile) < cutoff) mat(:, 3) = 0.

   do i = 1, j
   if (.not. ANY(same_line(1:j, i))) then
   do k = 1, j
      if (i /= k) then

         v1_ = 0.d0
         proj(i, k) = norm(vecprod(coord(i, :) + v1_ - coord(k, :), rutile))/norm(rutile)  !/ norm(coord(i,:)+v1_ -coord(k,:))
         if (abs(proj(i, k)) < cutoff) goto 21

         do i1_ = -1, 1
         do i2_ = -1, 1
         do i3_ = -1, 1
            v1_ = dble(i1_)*mat(:, 1) + dble(i2_)*mat(:, 2) + dble(i3_)*mat(:, 3)
            proj(i, k) = norm(vecprod(coord(i, :) + v1_ - coord(k, :), rutile))/norm(rutile)  !/ norm(coord(i,:)+v1_ -coord(k,:))
            if (abs(proj(i, k)) < cutoff) goto 21
         enddo
         enddo
         enddo

21       continue

         if (abs(proj(i, k)) < cutoff) then
            same_line(i, k) = .true.
            counti(i) = counti(i) + 1
            nneigh(i, counti(i)) = k
            nneigh(i, 0) = i
            ladder(i, 0) = 0.
            ladder(i, counti(i)) = scalprod(coord(k, :) - coord(i, :) - v1_, rutile)/norm(rutile)
            write (*, '(a,2i4,3f10.3)') 'vecprod/projection on rutile axis : ', i, k, proj(i, k), ladder(i, counti(i))
         endif

      endif
   enddo
   endif
   enddo

30 continue
   k = 0
   do i = 1, j
      if (counti(i) > 0) then
         k = k + 1
         kk = counti(i)
         call qsort_array(ladder(i, 0:kk), order(i, 0:kk))
         call qsort_adj_array(nneigh(i, 0:kk), order(i, 0:kk))
         write (*, '(a,i4,200i4)') 'line   : ', i, nneigh(i, 0:kk)
         write (*, '(a,i4,200f8.4)') 'ladder : ', i, ladder(i, 0:kk)
     if(manhattan_distance(coord(nneigh(i,0),:),coord(nneigh(i,1),:))>manhattan_distance(coord(nneigh(i,1),:),coord(nneigh(i,2),:)))then
            write (*, *) 'WARNING : DIMER PAIRS ARE SCRAMBLED,PLEASE USE ANOTHER INPUT FILES FOR THE COORDINATES'
            write (*, *) 'YOU CAN USE A SHIFT TO DEFINE THE BOUNDARY SUCH THAT DIMER PAIRS FIT IN THE UNITCELL'
            write (*, *) 'I swap them for you in the meanwhile...'
            ladder(i, 0) = 1.d20
            goto 30
         endif
      endif
   enddo

   write (*, *) 'build dimers out of the 1D lines of atoms along rutile axis'
   write (*, *) 'There are [x] 1D lines of dimers : ', k

   open (unit=1616, file='position_shifted')
   open (unit=10, file=filename)
   read (10, *)
   read (10, *)
   read (10, *)
   do
      read (10, *, end=73) temp, v
      write (1616, '(a,3f15.6)') trim(temp), v + shift_amp
   enddo
73 continue
   close (1616)
   close (10)

   open (unit=1414, file='mask_dimer')
   open (unit=1515, file='mask_dimer_coord')
   open (unit=1516, file='mask_dimer_coord_disorder')
   open (unit=1517, file='mask_dimer_coord_peierls')
   open (unit=1518, file='mask_dimer_coord_peierls_and_disorder')

   do i = 1, j
      if (counti(i) > 0) then
         kk = counti(i)
         nneigh(i, 1:kk + 1) = nneigh(i, 0:kk)
         kk = kk + 1
         write (*, *) 'line sites : ', nneigh(i, 1:kk)
         do l = 1, kk/2
            write (*, *) 'dimers are : ', nneigh(i, 2*l - 1), nneigh(i, 2*l)
            i1 = nneigh(i, 2*l - 1)
            i2 = nneigh(i, 2*l)

            if (l > 1 .and. l < kk/2) then
            write (*, *) 'V distance to other dimer (above): ', manhattan_distance(coord(i2, :), coord(nneigh(i, 2*(l + 1) - 1), :))
               write (*, *) 'V distance to other dimer (below): ', manhattan_distance(coord(i1, :), coord(nneigh(i, 2*(l - 1)), :))
            endif
            write (*, *) 'DIMER BOND LENGTH IS      (here) : ', manhattan_distance(coord(i2, :), coord(i1, :))

            write (1414, *) i1, i2
            write (1515, *) coord(i1, :)
            write (1515, *) coord(i2, :)
            write (1516, *) coord_disorder(i1, :)
            write (1516, *) coord_disorder(i2, :)
            write (1414, *) i2, i1
            write (1515, *) coord(i2, :)
            write (1515, *) coord(i1, :)
            write (1516, *) coord_disorder(i2, :)
            write (1516, *) coord_disorder(i1, :)
            v12 = coord(i2, :) - coord(i1, :)
            coord(i1, :) = coord(i1, :) + v12*peierls_move
            coord(i2, :) = coord(i2, :) - v12*peierls_move
            coord_disorder(i1, :) = coord_disorder(i1, :) + v12*peierls_move
            coord_disorder(i2, :) = coord_disorder(i2, :) - v12*peierls_move
            write (1518, *) coord_disorder(i1, :)
            write (1518, *) coord_disorder(i2, :)
            write (1518, *) coord_disorder(i2, :)
            write (1518, *) coord_disorder(i1, :)
            write (1517, *) coord(i1, :)
            write (1517, *) coord(i2, :)
            write (1517, *) coord(i2, :)
            write (1517, *) coord(i1, :)
         enddo
      endif
   enddo

   close (1414)
   close (1515)
   close (1516)
   close (1517)
   close (1518)

   open (unit=10, file=filename)
   read (10, *)
   read (10, *)
   read (10, *)

   open (unit=1415, file='position_peierls')
   open (unit=1416, file='position_peierls_and_disorder')

   j = 0
   do
      read (10, *, end=68) temp, v(1), v(2), v(3)
      if (TRIM(temp) == TRIM(ref)) then
         j = j + 1
         write (1415, '(a5,3f14.5)') TRIM(temp), coord(j, :)
         write (1416, '(a5,3f14.5)') TRIM(temp), coord_disorder(j, :)
      else
         write (1415, '(a5,3f14.5)') TRIM(temp), v(1:3)
         write (1416, '(a5,3f14.5)') TRIM(temp), v(1:3)
      endif
   enddo
68 continue
   close (10)

   close (1415)
   close (1416)

contains

   real(kind=DP) function manhattan_distance(v1, v2)
      implicit none
      integer :: i1_, i2_, i3_
      real(kind=DP) :: v1(3), v2(3), dd, temp, v1_(3)
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

