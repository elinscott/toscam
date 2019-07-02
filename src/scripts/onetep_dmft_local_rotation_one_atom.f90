program one_atom
   use genvar, only: dp
   use linalg
   implicit none
   real(dp) :: v1(3), v2(3), v3(3), site1(3), site2a(3), site2b(3), site3a(3), site3b(3), mat(3, 3), coord(3)
   integer  :: j, i, jj, kk, kkk, jjj, lll_

   open (unit=1414, file='mask_local_rotations', form='unformatted')

   write (*, *) 'how many atom types'
   read (*, *) kk

   do jj = 1, kk

      mat = 0; coord = 0.; 
      write (*, *) 'ENTER HOW MANY OCCURENCES'
      read (*, *) kkk
      write (*, *) 'ENTER ROTATION LABEL (=0 if no rotation needed)'
      read (*, *) lll_
      if (lll_ == 0) goto 21

      write (*, *) 'ENTER COORDINATES OF ATOMIC SITE'
      read (*, *) (site1(j), j=1, 3)
      write (*, *) 'ENTER COORDINATES OF SITE FOR X AXIS, A'
      read (*, *) (site2a(j), j=1, 3)
      write (*, *) 'ENTER COORDINATES OF SITE FOR X AXIS, B'
      read (*, *) (site2b(j), j=1, 3)
      write (*, *) 'ENTER COORDINATES OF SITE FOR Y AXIS, A'
      read (*, *) (site3a(j), j=1, 3)
      write (*, *) 'ENTER COORDINATES OF SITE FOR Y AXIS, B'
      read (*, *) (site3b(j), j=1, 3)

      v1 = site2b - site2a
      v2 = site3b - site3a

      v1 = v1/norme(v1)
      v2 = v2/norme(v2)

      write (*, *) 'V2 is                               : ', v2
      v2 = v2 - scalprod(v2, v1)*v1
      v2 = v2/norme(v2)
      write (*, *) 'V2 orthogonalized (gram-schmidt) is : ', v2

      v3 = vecprod(v1, v2)

      mat(:, 1) = v1
      mat(:, 2) = v2
      mat(:, 3) = v3
      coord = site1

      write (*, *) 'MAXVAL O * O^T : ', maxval(abs(matmul(mat, transpose(mat))))
      write (*, *) 'MINVAL O * O^T : ', minval(abs(matmul(mat, transpose(mat))))

21    continue
      do jjj = 1, kkk
         write (1414) lll_, coord(:), mat(:, :)
      enddo

   enddo

   close (1414)

end program
