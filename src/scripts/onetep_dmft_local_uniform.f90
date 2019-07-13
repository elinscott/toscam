program local_uniform_mask
   use genvar, only: DP
   use linalg
   implicit none

   integer                    :: jj, i, j, k, l, m
   integer                    :: k1, k2, k3
   integer, allocatable :: uniform(:, :), label_in_big(:), label_in_hub(:)
   character(30), allocatable :: atom(:), atom_small(:)
   real(kind=DP), allocatable :: normv(:), coord(:, :), coord_small(:, :)
   real(kind=DP)                    :: v(3)
   character(30)              :: ref, filename, filename2, filename3, temp, species
   integer                    :: natom, natom_small, kk1, kk2, kk3
   real(kind=DP)                    :: coord_tmp(3), www(3), ex(3), ey(3), ez(3), bx(3), by(3), bz(3), dist, cutoff, vvv(3)
   real(kind=DP)                    :: shift(3), max_dist, max_dist2, d1, d2, d3
   logical                    :: check

   open (unit=1010, file='uniform.input')
   write (*, *) 'please enter file name for large cell (onetep lattice)'
   read (1010, *) filename
   write (*, *) 'how many atoms in the lattice?'
   read (1010, *) natom
   write (*, *) 'please enter file name for small cell (original small unitcell)'
   read (1010, *) filename2
   write (*, *) 'how many atoms in the small unitcell?'
   read (1010, *) natom_small
   write (*, *) 'center for local coordinates (symbol, example V,Cu,O...)'
   read (1010, *) ref
   write (*, *) ' shift '
   read (1010, *) shift
   close (1010)

   allocate (atom(natom), uniform(natom, 2), label_in_hub(natom))
   allocate (coord(natom, 3))
   allocate (coord_small(natom_small, 3))
   allocate (atom_small(natom_small))
   allocate (label_in_big(natom_small))

   jj = 0
   label_in_hub = 0
   open (unit=10, file=filename)
   read (10, *) bx
   read (10, *) by
   read (10, *) bz
   do k = 1, natom
      read (10, *) temp, v(1), v(2), v(3)
      atom(k) = temp; coord(k, :) = v
      write (*, '(a,a,3f10.3)') 'atom,coord : ', TRIM(atom(k)), coord(k, :)
      if (trim(temp) == trim(ref)) then
         jj = jj + 1
         label_in_hub(k) = jj
      endif
   enddo
   close (10)

   open (unit=10, file=filename2)
   read (10, *) ex
   read (10, *) ey
   read (10, *) ez
   do k = 1, natom_small
      read (10, *) temp, v(1), v(2), v(3)
      atom_small(k) = temp; coord_small(k, :) = v
      write (*, '(a,a,3f10.3)') 'atom,coord : ', TRIM(atom_small(k)), coord_small(k, :)
   enddo
   close (10)

   kk1 = NINT(norme(bx)/norme(ex))
   kk2 = NINT(norme(by)/norme(ey))
   kk3 = NINT(norme(bz)/norme(ez))

   write (*, *) ' large unitcell is [x,y,z] x small unitcell : ', kk1, kk2, kk3

   write (*, *) 'checking that small unitcell is contained in the large one'
   label_in_big = 0

   do i = 1, natom_small
      do j = 1, natom
         d1 = abs(coord_small(i, 1) - coord(j, 1))
         if (d1 <= 1.d-2) then
            d2 = abs(coord_small(i, 2) - coord(j, 2))
            if (d2 <= 1.d-2) then
               d3 = abs(coord_small(i, 3) - coord(j, 3))
               if (d3 <= 1.d-2) then
                  label_in_big(i) = j
                  exit
               endif
            endif
         endif
      enddo
      if (label_in_big(i) == 0) then
         write (*, *) 'site not found in big lattice : ', i
         write (*, *) 'coord small lattice : ', coord_small(i, :)
         stop
      endif
   enddo

   write (*, *) 'checks done, now build output files'

   open (unit=1221, file='mask_uniform')
   open (unit=1310, file='mask_uniform_coord')

   do i = 1, natom
      jj = 0
      do j = 1, natom_small
      do k1 = -kk1 - 3, kk1 + 3
      do k2 = -kk2 - 3, kk2 + 3
      do k3 = -kk3 - 3, kk3 + 3
         v = dble(k1)*ex + dble(k2)*ey + dble(k3)*ez
         v = v + coord_small(j, :)
         if (norme(v - coord(i, :)) < 1.d-3) then
            jj = j
            goto 89
         endif
      enddo
      enddo
      enddo
      enddo

89    continue

      if (jj == 0) then
         write (*, *) 'site of large lattice not found in unitcell'
         stop
      endif

      if (trim(atom(i)) == trim(ref)) then
         write (*, *) 'SITE [x] of big lattice is obtained from [j] : ', i, jj
         uniform(i, 1) = i
         uniform(i, 2) = label_in_big(jj)
         write (1310, *) coord(uniform(i, 1), :) + shift
         write (1310, *) coord(uniform(i, 2), :) + shift
         write (1221, *) label_in_hub(uniform(i, 1:2))
      endif

   enddo

   close (1221)
   close (1310)

end program
