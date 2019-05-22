program local_rotations
   use genvar, only: pi
   use linalg
   use geometry2, only: plane_exp_point_dist_3d_, angle_rad_3d_
   implicit none

   integer                    :: i, j, k, l, m, jj
   integer                    :: ii(3), k1, k2, k3, k_ez, ninplane, noutofplane
   logical, allocatable :: neigh(:, :)
   character(30), allocatable :: atom(:), atom_small(:)
   integer, allocatable :: nneigh(:), label_in_big(:), his_neigh(:, :), inplane(:, :), outofplane(:, :)
   real(8), allocatable :: normv(:), coord(:, :), coord_small(:, :), mat(:, :, :), mat_small(:, :, :)
   real(8)                    :: flip, pos_shift, shift(3), a(3, 3), aa(3, 3), v(3)
   character(30)              :: ref, filename, filename2, filename3, temp, species
   integer                    :: natom, maxneigh, natom_small, k4, itot, dir_reference
   real(8)                    :: coord_tmp(3), www(3), ex(3), ey(3), ez(3), dist, cutoff, vvv(3), rutile(3)
   real(8)                    :: min_dist, max_dist, max_dist2, d1, d2, d3
   logical                    :: check, keep_perp_to_rutile, all_planes
   logical                    :: checkaxisfile
   logical, parameter :: testingloc = .false.

   min_dist = 2.d-2
   dir_reference = 3
   inquire (file='axis_file', exist=checkaxisfile)

   if (.not. testingloc) open (unit=1010, file='rotation.input')
   call system("rm ./large_cell")
   filename = 'large_cell'
   write (*, *) 'please enter file name for small cell (original lattice)'
   if (.not. testingloc) read (1010, *) filename2
   write (*, *) 'how many atoms in the small unitcell?'
   if (.not. testingloc) read (1010, *) natom_small
   write (*, *) 'center for local coordinates (symbol, example V,Cu,O...)'
   if (.not. testingloc) read (1010, *) ref
   write (*, *) 'enter cutoff parameters to determine which Oxygen are apicals (typically 0.03'
   if (.not. testingloc) read (1010, *) cutoff
   write (*, *) 'please enter rutile c axis (3 real numbers) : '
   if (.not. testingloc) read (1010, *) rutile
   write (*, *) 'please enter max dist for bonds (in Angstrom, 2.2 for example) : '
   if (.not. testingloc) read (1010, *) max_dist
   write (*, *) 'please enter file name for small cell (modified coordinates, optional output)'
   if (.not. testingloc) read (1010, *) filename3
   write (*, *) 'please enter flip value (+1 or -1) to flip y axis'
   if (.not. testingloc) read (1010, *) flip
   write (*, *) 'please enter if you want to discard sites in the octahedra which are perpendicular to rutile axis'
   if (.not. testingloc) read (1010, *) keep_perp_to_rutile
   write (*, *) 'include all 4-sites to define the 4-site plane, or only triplet of points?'
   if (.not. testingloc) read (1010, *) all_planes
   if (.not. testingloc) close (1010)

   if (maxloci(abs(rutile)) == dir_reference) then
      write (*, *) 'WARNING - YOUR RUTILE AXIS IS ALONG REFERENCE DIRECTION'
      write (*, *) 'THE REFERENCE DIR IS USED TO DETERMINE THE DIRECTION OF THE SECOND AXIS'
      dir_reference = minloci(abs(rutile) + (/(1.d-8*dble(jj), jj=1, 3)/))
      write (*, *) 'NOW REFERENCE DIR IS (x,y,z) = ', dir_reference
      write (*, *) 'your rutile axis = ', rutile
   endif

   if (testingloc) then
      all_planes = .false.
      keep_perp_to_rutile = .false.
      rutile(1) = 5.743
      rutile(2) = 0.
      rutile(3) = 0.
      max_dist = 2.2
      ref = 'V'
      flip = 1.d0
      natom_small = 24
      cutoff = 0.03
      filename2 = 'test.small.cell'
      filename3 = 'not.present'
   endif

   call duplicate_unitcell

   max_dist2 = max_dist**2
   if (max_dist > max_dist2) then
      write (*, *) 'sorry, but the code is only working for max_dist>1, please update the code'
      stop
   endif

   allocate (neigh(natom, natom))
   allocate (atom(natom), mat(natom, 3, 3))
   allocate (nneigh(natom), his_neigh(natom, 30), inplane(natom, 200), outofplane(natom, 200))
   allocate (normv(natom), coord(natom, 3))
   allocate (label_in_big(natom_small), coord_small(natom_small, 3), atom_small(natom_small), mat_small(natom_small, 3, 3))

!AXIS Z ALONG APICAL
!AXIS X ALONG RUTILE AXIS
!AXIS Y given by X and Z

   his_neigh = 0
   neigh = .false.
   open (unit=10, file=filename)
   read (10, *)
   read (10, *)
   read (10, *)
   do k = 1, natom
      read (10, *) temp, v(1), v(2), v(3)
      atom(k) = temp; coord(k, :) = v
      write (*, '(a,a,3f10.3)') 'atom,coord : ', TRIM(atom(k)), coord(k, :)
   enddo
   close (10)
   !call system("rm large_cell")

   open (unit=10, file=filename2)
   read (10, *)
   read (10, *)
   read (10, *)
   do k = 1, natom_small
      read (10, *) temp, v(1), v(2), v(3)
      atom_small(k) = temp; coord_small(k, :) = v
      write (*, '(a,a,3f10.3)') 'atom,coord : ', TRIM(atom_small(k)), coord_small(k, :)
   enddo
   close (10)

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

   write (*, *) 'computing bonds'
   do i = 1, natom
      do j = i + 1, natom
         d1 = (coord(i, 1) - coord(j, 1))**2.
         if (d1 <= max_dist2) then
            d2 = (coord(i, 2) - coord(j, 2))**2.
            if (d1 + d2 <= max_dist2) then
               d3 = (coord(i, 3) - coord(j, 3))**2.
               if (d1 + d2 + d3 <= max_dist2) then
                  dist = sqrt(d1 + d2 + d3)
                  if (dist < max_dist) then
                     neigh(i, j) = .true.
                     neigh(j, i) = .true.
                  endif
                  if (dist < min_dist) then
                     write (*, *) 'ERROR - minimum distance reached (atoms at same site??) , distance =  ', dist
                     write (*, *) 'min dist    : ', min_dist
                     write (*, *) 'coord atoms : ', coord(i, :)
                     write (*, *) 'coord atoms : ', coord(j, :)
                     stop
                  endif
               endif
            endif
         endif
      enddo
   enddo

   nneigh = 0
   maxneigh = 0

   do i = 1, natom
      j = count(neigh(i, 1:natom))
      if (TRIM(atom(i)) == TRIM(ref)) then
         nneigh(i) = j
      endif
      k = 0
      if (nneigh(i) > 0) then
      do j = 1, natom
         if (neigh(i, j)) then
            k = k + 1
            his_neigh(i, k) = j
         endif
      enddo
      if (k /= nneigh(i)) then
         write (*, *) 'counting issues....'
         write (*, *) 'k is = ', k
         write (*, *) 'number of neighbors is = ', nneigh(i)
         stop
      endif
      endif
   enddo

   do i = 1, natom
      if (nneigh(i) > 0) then
         write (*, *) 'atom i has [x] neighbors : ', i, nneigh(i)
         do j = 1, nneigh(i)
            write (*, '(a,i4,3f10.3)') 'his neighbor is =', his_neigh(i, j), coord(his_neigh(i, j), :)
         enddo
      endif
   enddo

   inplane = 0
   outofplane = 0
   mat = 0.d0

!===========================================================================!
!===========================================================================!

   do i = 1, natom
      if (nneigh(i) > 0) then

         write (*, *) '-------- ATOM --------', i
         k = nneigh(i)
         if (k /= 6) goto 78

         if (checkaxisfile) then
            write (*, *) 'AXIS ARE NOT CHOSEN ACCORDING TO RUTILE DIRECTION'
            write (*, *) 'BUT RATHER GIVEN AS CONSTANTS FROM THE FILE AXIS_FILE'
            write (*, *) 'SKIP IN-PLANE AND OUT-OF-PLANE CALCULATIONS'
            goto 3031
         endif

         if (.not. all_planes) then
            j = 0
            do k1 = 1, k
            do k2 = k1 + 1, k
            do k3 = k2 + 1, k
call plane_exp_point_dist_3d_(coord(his_neigh(i, k1), :), coord(his_neigh(i, k2), :), coord(his_neigh(i, k3), :), coord(i, :), dist)
               if (dist < cutoff) then
                  write (*, '(a,3i4,2f10.3)') 'k1,k2,k3,dist,cutoff : ', k1, k2, k3, dist, cutoff
                  if (.not. ANY((k1 - inplane(i, :)) == 0)) then
                     j = j + 1
                     inplane(i, j) = k1
                  endif
                  if (.not. ANY((k2 - inplane(i, :)) == 0)) then
                     j = j + 1
                     inplane(i, j) = k2
                  endif
                  if (.not. ANY((k3 - inplane(i, :)) == 0)) then
                     j = j + 1
                     inplane(i, j) = k3
                  endif
               endif

            enddo
            enddo
            enddo
         else
            j = 0
            do k1 = 1, k
            do k2 = k1 + 1, k
            do k3 = k2 + 1, k
            do k4 = k3 + 1, k

call plane_exp_point_dist_3d_(coord(his_neigh(i, k1), :), coord(his_neigh(i, k2), :), coord(his_neigh(i, k3), :), coord(i, :), dist)
               if (dist < cutoff) then
call plane_exp_point_dist_3d_(coord(his_neigh(i, k1), :), coord(his_neigh(i, k2), :), coord(his_neigh(i, k4), :), coord(i, :), dist)
                  if (dist < cutoff) then
call plane_exp_point_dist_3d_(coord(his_neigh(i, k1), :), coord(his_neigh(i, k3), :), coord(his_neigh(i, k4), :), coord(i, :), dist)
                     if (dist < cutoff) then
call plane_exp_point_dist_3d_(coord(his_neigh(i, k2), :), coord(his_neigh(i, k3), :), coord(his_neigh(i, k4), :), coord(i, :), dist)
                        if (dist < cutoff) then

                           if (keep_perp_to_rutile) then
                              itot = 0
                              if (abs(scalprod(coord(his_neigh(i, k1), :) - coord(i, :), rutile))/norme(rutile) < cutoff) then
                                 itot = itot + 1
                              endif
                              if (abs(scalprod(coord(his_neigh(i, k2), :) - coord(i, :), rutile))/norme(rutile) < cutoff) then
                                 itot = itot + 1
                              endif
                              if (abs(scalprod(coord(his_neigh(i, k3), :) - coord(i, :), rutile))/norme(rutile) < cutoff) then
                                 itot = itot + 1
                              endif
                              if (abs(scalprod(coord(his_neigh(i, k4), :) - coord(i, :), rutile))/norme(rutile) < cutoff) then
                                 itot = itot + 1
                              endif
                              if (itot < 2) goto 214
                              if (itot > 2) then
                                 write (*, *) 'ERROR : build axis along octahedra routine'
                                 write (*, *) 'ATOM HAS MORE THAN 2 APICAL, TOTAL =', itot
                                 write (*, *) 'CRITICAL / WILL NOW STOP'
                                 write (*, *) 'RUTILE AXIS : ', rutile
                                 write (*, *) 'site : ', coord(i, :)
                                 do j = 1, size(his_neigh, 2)
                                    if (his_neigh(i, j) /= 0) write (*, *) 'neighbors : ', coord(his_neigh(i, j), :)
                                 enddo
                                 stop
                              endif
                           endif

                           write (*, '(a,3i4,2f10.3)') 'k1,k2,k3,dist,cutoff : ', k1, k2, k3, dist, cutoff
                           if (.not. ANY((k1 - inplane(i, :)) == 0)) then
                              j = j + 1
                              inplane(i, j) = k1
                           endif
                           if (.not. ANY((k2 - inplane(i, :)) == 0)) then
                              j = j + 1
                              inplane(i, j) = k2
                           endif
                           if (.not. ANY((k3 - inplane(i, :)) == 0)) then
                              j = j + 1
                              inplane(i, j) = k3
                           endif
                           if (.not. ANY((k4 - inplane(i, :)) == 0)) then
                              j = j + 1
                              inplane(i, j) = k4
                           endif
                           if (keep_perp_to_rutile .and. j == 4) goto 215

                        endif
                     endif
                  endif
               endif

214         enddo
            enddo
            enddo
            enddo
215         continue

         endif

         write (*, *) 'TOTAL IDENTIFIED IN-plane : ', j
         ninplane = j

         if (ninplane /= 4) then
            write (*, *) 'this routine was written for octahedra, please update it....'
            write (*, *) 'the atoms are : '
            do k1 = 1, k
               write (*, '(200f12.4)') coord(his_neigh(i, k1), :)
            enddo
            write (*, '(a,200f12.4)') 'the site is : ', coord(i, :)
            stop
         endif

         write (*, '(a,i5,a,10i4)') 'IN-plane neighbors label of atom [x] : ', i, ' are ', inplane(i, 1:ninplane)

         j = 0
         do k = 1, nneigh(i)
            if (.not. ANY(inplane(i, :) - k == 0)) then
               j = j + 1
               outofplane(i, j) = k
            endif
         enddo
         noutofplane = j

         write (*, *) 'NUMBER OF OUT OF PLANE : ', noutofplane
         write (*, *) 'OUT OF PLANE SITE      : ', outofplane(i, 1:noutofplane)

         do k = 1, ninplane
            inplane(i, k) = his_neigh(i, inplane(i, k))
         enddo
         do k = 1, noutofplane
            outofplane(i, k) = his_neigh(i, outofplane(i, k))
         enddo
         write (*, '(a,i5,a,10i4)') 'IN-plane oxygen neighbors of atom [x]    : ', i, ' are ', inplane(i, 1:ninplane)
         write (*, '(a,i5,a,10i4)') 'OUTOF-plane oxygen neighbors of atom [x] : ', i, ' are ', outofplane(i, 1:noutofplane)

         write (*, *) 'ANGLE IN-PLANE-TO-VO WITH VO-TO-OUT-OF_PLANE : '
         do k1 = 1, ninplane
            do k2 = 1, noutofplane
               write (*, *) 'inplane atom [x] and out of plane [y] : ', k1, k2
              write (*, *) 'angle = ', (180.d0/pi)*angle_rad_3d_(coord(inplane(i, k1), :), coord(i, :), coord(outofplane(i, k2), :))
            enddo
         enddo

         call choose_my_axis

3031     continue

         if (.not. checkaxisfile) then
            call set_my_axis
         else
            call set_given_axis_compute_all
         endif

         mat(i, :, 1) = ex
         mat(i, :, 2) = ey
         mat(i, :, 3) = ez

         if (scalprod(vecprod(ex, ey), ez) < 0.0) then
            write (*, *) 'DANGER/ERROR : your new axis are not properly oriented'
            write (*, *) ' ex x ey     = ', vecprod(ex, ey)
            write (*, *) '(ex x ey).ez = ', scalprod(vecprod(ex, ey), ez)
            stop
         endif

         write (*, *) 'Mat * Mat^T  (minval) : ', minval(matmul(mat(i, :, :), transpose(mat(i, :, :))))
         write (*, *) 'Mat * Mat^T  (maxval) : ', maxval(matmul(mat(i, :, :), transpose(mat(i, :, :))))
         write (*, *) '==========================='

      endif
78 enddo

!===========================================================================!
!===========================================================================!

   mat_small = 0.d0
   do i = 1, natom_small
      do j = 1, natom
         d1 = abs(coord_small(i, 1) - coord(j, 1))
         if (d1 <= 1.d-4) then
            d2 = abs(coord_small(i, 2) - coord(j, 2))
            if (d2 <= 1.d-4) then
               d3 = abs(coord_small(i, 3) - coord(j, 3))
               if (d3 <= 1.d-4) then
                  mat_small(i, :, :) = mat(j, :, :)
                  exit
               endif
            endif
         endif
      enddo
   enddo

   INQUIRE (file=filename3, EXIST=check)
   if (check) then
      open (unit=2143, file=filename3)
      open (unit=2144, file='mask_local_rotations_coord_of_input2', form='unformatted')
   endif

   k = 0
   open (unit=1414, file='mask_local_rotations', form='unformatted')
   do i = 1, natom_small
      if (nneigh(label_in_big(i)) > 0) then
         if (maxval(abs(mat_small(i, :, :))) < 1.d-3) then
            write (*, *) 'most probably a bug and an index was missed or wronged'
            write (*, *) 'atom [x] and total [x] : ', i, natom_small
            write (*, *) 'coord atom             : ', coord_small(i, :)
            write (*, *) 'mat_small              : ', mat_small(i, :, :)
            write (*, *) 'site in big lattice    : ', label_in_big(i)
            write (*, *) 'number of neighbors    : ', nneigh(label_in_big(i))
            write (*, *) 'atom type              : ', TRIM(atom(i))
            stop
         endif
         k = k + 1

         if (check) then
            read (2143, *, end=98) temp, coord_tmp(1), coord_tmp(2), coord_tmp(3)
            write (2144) i, coord_tmp(:), mat_small(i, :, :)
         endif

         write (1414) i, coord_small(i, :), mat_small(i, :, :)

      endif
   enddo
   close (1414)

   if (check) then
      close (2143)
      close (2144)
   endif

   write (*, *) '-----> THERE ARE [x] local rotations : ', k

   if (.false.) then
98    continue
      write (*, *) ' error, not enough entries in file :  ', trim(filename3)
      stop
   endif

contains
!----------------!
!----------------!
!----------------!
!----------------!
!----------------!

   subroutine set_given_axis_compute_all
      implicit none
      real(8) :: vectors(6, 3), vref1(3), vref2(3), vref3(3)
      integer :: ijk, aa1, aa2, aa3
      real(8) :: projs(3, 6)

      do ijk = 1, 6
         vectors(ijk, :) = -coord(i, :) + coord(his_neigh(i, ijk), :)
      enddo
      do ijk = 1, 6
         vectors(ijk, :) = vectors(ijk, :)/norm(vectors(ijk, :))
      enddo

      inquire (file='axis_file', exist=checkaxisfile)
      if (.not. checkaxisfile) then
         write (*, *) 'ERROR : axis file supposed to be present'
         stop
      endif
      open (unit=90111, file='axis_file')
      read (90111, *) (vref1(ijk), ijk=1, 3)
      read (90111, *) (vref2(ijk), ijk=1, 3)
      read (90111, *) (vref3(ijk), ijk=1, 3)
      vref1 = vref1/norm(vref1)
      vref2 = vref2/norm(vref2)
      vref3 = vref3/norm(vref3)
      do ijk = 1, 6
         projs(1, ijk) = scalprod(vref1, vectors(ijk, :))
         projs(2, ijk) = scalprod(vref2, vectors(ijk, :))
         projs(3, ijk) = scalprod(vref3, vectors(ijk, :))
      enddo
      aa1 = maxloci(projs(1, :))
      aa2 = maxloci(projs(2, :))
      aa3 = maxloci(projs(3, :))
      close (90111)

      write (*, *) '==========================='
      write (*, *) 'LOCAL AXIS SETUP'

      ex = vectors(aa1, :)
      ex = ex/norm(ex)
      write (*, '(a,3f13.4)') 'AXE X : ', ex
      ez = vectors(aa3, :)
      ez = ez/norm(ez)
      write (*, '(a,3f13.4)') 'AXE Z : ', ez

      write (*, *) 'correction due to small non-orthogonality (disorder, distortion, other...)'
      ez = ez - scalprod(ez, ex)*ex
      ez = ez/norm(ez)
      write (*, '(a,3f13.4)') 'AXE Z : ', ez

      ey = vecprod(ez, ex)
      ey = ey/norm(ey)
      write (*, '(a,3f13.4)') 'AXE Y : ', ey
      write (*, *) 'SCALPROD AXEY - supposed axis : ', scalprod(ey, vectors(aa2, :))
      write (*, *) '==========================='

   end subroutine

!----------------!
!----------------!
!----------------!
!----------------!
!----------------!
!----------------!
!----------------!
!----------------!
!----------------!
!----------------!

   subroutine set_my_axis
      write (*, *) '==========================='
      write (*, *) 'LOCAL AXIS SETUP'

      ex = rutile
      !ex=-coord(i,:)+coord(outofplane(i,k2),:)
      ex = ex/norm(ex)
      write (*, '(a,3f13.4)') 'AXE X : ', ex

      ez = -coord(i, :) + coord(inplane(i, k_ez), :)
      ez = ez/norm(ez)
      write (*, '(a,3f13.4)') 'AXE Z : ', ez

      write (*, *) 'correction due to small non-orthogonality (disorder, distortion, other...)'
      ez = ez - scalprod(ez, ex)*ex
      ez = ez/norm(ez)
      write (*, '(a,3f13.4)') 'AXE Z : ', ez

      ey = vecprod(ez, ex)
      ey = ey/norm(ey)
      write (*, '(a,3f13.4)') 'AXE Y : ', ey

      write (*, *) '==========================='
   end subroutine

!----------------!
!----------------!
!----------------!
!----------------!
!----------------!

   subroutine choose_my_axis
      implicit none
      real(8)  ::  tt, in_plane_angles(ninplane)
      real(8)  ::  www1(3), www2(3)
      integer  ::  both_directions(2)

      write (*, *) '---------------------------------'
      write (*, *) 'local geometry informations : '

      normv = 0.d0
      do k = 1, noutofplane
  normv(k) = scalprod(-coord(i, :) + coord(outofplane(i, k), :), rutile)/norm(coord(i, :) - coord(outofplane(i, k), :))/norm(rutile)
         write (*, *) ' k,normv : ', k, normv(k)
      enddo
      k2 = maxloci(normv(1:noutofplane))

      normv = 0.d0
      do k = 1, ninplane
        normv(k) = scalprod(-coord(i, :) + coord(inplane(i, k), :), rutile)/norm(rutile)/norm(coord(i, :) - coord(inplane(i, k), :))
         write (*, *) ' k,normv : ', k, normv(k)
      enddo
      k3 = maxloci(normv(1:ninplane))

     write(*,'(a,f8.3)')    'MY APICAL OXYGEN AXIS             : ', scalprod(-coord(i,:)+coord(outofplane(i,k2),:),rutile)/norm(coord(i,:)-coord(outofplane(i,k2),:))/norm(rutile)
     write(*,'(a,f8.3)')    'MY INPLANE OXYGEN AXIS            : ', scalprod(-coord(i,:)+coord(inplane(i,k3),:),rutile)/norm(rutile)/norm(coord(i,:)-coord(inplane(i,k3),:))

      do k = 1, ninplane
         write (*, '(a,2i5,f8.3)') 'DISTANCE FROM V TO IN PLANE     : ', i, k, norm(coord(i, :) - coord(inplane(i, k), :))
      enddo
      do k = 1, noutofplane
         write (*, '(a,2i5,f8.3)') 'DISTANCE FROM V TO OUTOFPLANE   : ', i, k, norm(coord(i, :) - coord(outofplane(i, k), :))
      enddo

      vvv = coord(outofplane(i, 1), :) - coord(outofplane(i, 2), :)

      write (*, '(a,3f10.3)') 'OCTAHEDRON APICAL-APICAL DIRECTION : ', vvv
      write (*, '(a,3f10.3)') 'RUTILE DIRECTION                   : ', rutile
      write (*, *) 'ANGLE WITH RUTILE AXIS             : ', dacos(scalprod(vvv, rutile)/norm(vvv)/norm(rutile))*180.d0/pi

      do k1 = 1, ninplane
         www = coord(inplane(i, k1), :) - coord(i, :)
         write (*, *) 'APICAL-APICAL ANGLE WITH IN-PLANE  : ', dacos(scalprod(vvv, www)/norm(vvv)/norm(www))*180.d0/pi
      enddo

      do k1 = 1, ninplane
         www = coord(inplane(i, k1), :) - coord(i, :)
         in_plane_angles(k1) = dacos(scalprod(rutile, www)/norm(rutile)/norm(www))*180.d0/pi
         write (*, *) 'IN-PLANE WITH RUTILE               : ', in_plane_angles(k1)
      enddo

      both_directions(1) = minloci(abs(in_plane_angles - 90.d0))
      tt = in_plane_angles(both_directions(1))
      in_plane_angles(both_directions(1)) = 900000.
      both_directions(2) = minloci(abs(in_plane_angles - 90.d0))
      in_plane_angles(both_directions(1)) = tt

      write (*, *) 'ANGLES PERP to RUTILES ARE : ', in_plane_angles(both_directions(1:2))
      www1 = coord(inplane(i, both_directions(1)), :) - coord(i, :)
      www2 = coord(inplane(i, both_directions(2)), :) - coord(i, :)
      write (*, *) 'direction 1 is : ', www1
      write (*, *) 'direction 2 is : ', www2

      if (flip*www1(dir_reference) > flip*www2(dir_reference)) then
         k_ez = both_directions(1)
      else
         k_ez = both_directions(2)
      endif

      write (*, *) '---------------------------------'

      return
   end subroutine

!----------------!
!----------------!
!----------------!
!----------------!
!----------------!
!----------------!
!----------------!
!----------------!
!----------------!
!----------------!

   subroutine duplicate_unitcell
      open (unit=1010, file='duplicate.input')
      write (1010, *) trim(filename2)
      write (1010, *) - 1, 1
      write (1010, *) - 1, 1
      write (1010, *) - 1, 1
      write (1010, *) 'ALL'
      write (1010, *) 0.d0
      close (1010)
      natom = natom_small*3*3*3
      call system(" onetep.dmft.duplicate.unitcell ")
      call system("rm duplicate.input")
      call system(" mv unitcell_output "//trim(filename))
   end subroutine

!----------------!
!----------------!
!----------------!
!----------------!

end program
