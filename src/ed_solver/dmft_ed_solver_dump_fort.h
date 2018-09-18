
   subroutine dump_to_fort()

      use mpirout, only: mpi_stop

      implicit none

      integer :: j, i

      if(rank == 0)then

         do j = 1, size(G(1)%correl(1, 2)%fctn(1, 1, :))
            !normal
            write(1933, *) aimag(G(1)%correl(1, 2)%freq%vec(j)), &
                 real(G(1)%correl(1, 2)%fctn(1, 1, j)), aimag(G(1)%correl(1, &
                 2)%fctn(1, 1, j))
            write(3933, *) aimag(GF(1)%correl(1, 1)%freq%vec(j)), &
                 real(GF(1)%correl(1, 1)%fctn(1, 1, j)), aimag(GF(1)%correl(1, &
                 1)%fctn(1, 1, j))
            write(4933, *) aimag(GF(1)%correl(2, 2)%freq%vec(j)), &
                 real(GF(1)%correl(2, 2)%fctn(1, 1, j)), aimag(GF(1)%correl(2, &
                 2)%fctn(1, 1, j))
            write(9001, *) aimag(GNAMBU%correl(1, 2)%freq%vec(j)), &
                 real(GNAMBU%correl(1, 2)%fctn(1, 1, j)), &
                 aimag(GNAMBU%correl(1, 2)%fctn(1, 1, j))
            write(5001, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(1, 1, &
                 j)), aimag(SNAMBU%fctn(1, 1, j))
         enddo


         do j = 1, size(G(1)%correl(1, 2)%fctn(1, 1, :))
            write(1934, *) aimag(G(1)%correl(1, 2)%freq%vec(j)), &
                 real(G(1)%correl(1, 2)%fctn(1, 2, j)), aimag(G(1)%correl(1, &
                 2)%fctn(1, 2, j))
            write(1935, *) aimag(G(1)%correl(1, 2)%freq%vec(j)), &
                 real(G(1)%correl(1, 2)%fctn(2, 1, j)), aimag(G(1)%correl(1, &
                 2)%fctn(2, 1, j))
            write(1936, *) aimag(G(1)%correl(1, 2)%freq%vec(j)), &
                 real(G(1)%correl(1, 2)%fctn(2, 2, j)), aimag(G(1)%correl(1, &
                 2)%fctn(2, 2, j))
            write(1937, *) aimag(G(1)%correl(1, 2)%freq%vec(j)), &
                 real(G(1)%correl(1, 2)%fctn(1, 3, j)), aimag(G(1)%correl(1, &
                 2)%fctn(1, 3, j))
            write(1938, *) aimag(G(2)%correl(1, 2)%freq%vec(j)), &
                 real(G(2)%correl(1, 2)%fctn(1, 1, j)), aimag(G(2)%correl(1, &
                 2)%fctn(1, 1, j))

            write(2933, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                 real(G(1)%correl(2, 1)%fctn(1, 1, j)), aimag(G(1)%correl(2, &
                 1)%fctn(1, 1, j))
            write(2934, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                 real(G(1)%correl(2, 1)%fctn(1, 2, j)), aimag(G(1)%correl(2, &
                 1)%fctn(1, 2, j))
            write(2935, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                 real(G(1)%correl(2, 1)%fctn(2, 1, j)), aimag(G(1)%correl(2, &
                 1)%fctn(2, 1, j))
            write(2936, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                 real(G(1)%correl(2, 1)%fctn(2, 2, j)), aimag(G(1)%correl(2, &
                 1)%fctn(2, 2, j))

            !anomalous

            if(associated(GF(1)%correl(1, 1)%freq%vec)) then
               write(3934, *) aimag(GF(1)%correl(1, 1)%freq%vec(j)), &
                    real(GF(1)%correl(1, 1)%fctn(1, 2, j)), &
                    aimag(GF(1)%correl(1, 1)%fctn(1, 2, j))
               write(3935, *) aimag(GF(1)%correl(1, 1)%freq%vec(j)), &
                    real(GF(1)%correl(1, 1)%fctn(2, 1, j)), &
                    aimag(GF(1)%correl(1, 1)%fctn(2, 1, j))
               write(3936, *) aimag(GF(1)%correl(1, 1)%freq%vec(j)), &
                    real(GF(1)%correl(1, 1)%fctn(2, 2, j)), &
                    aimag(GF(1)%correl(1, 1)%fctn(2, 2, j))
               write(3937, *) aimag(GF(1)%correl(1, 1)%freq%vec(j)), &
                    real(GF(1)%correl(1, 1)%fctn(1, 3, j)), &
                    aimag(GF(1)%correl(1, 1)%fctn(1, 3, j))
               write(4934, *) aimag(GF(1)%correl(2, 2)%freq%vec(j)), &
                    real(GF(1)%correl(2, 2)%fctn(1, 2, j)), &
                    aimag(GF(1)%correl(2, 2)%fctn(1, 2, j))
               write(4935, *) aimag(GF(1)%correl(2, 2)%freq%vec(j)), &
                    real(GF(1)%correl(2, 2)%fctn(2, 1, j)), &
                    aimag(GF(1)%correl(2, 2)%fctn(2, 1, j))
               write(4936, *) aimag(GF(1)%correl(2, 2)%freq%vec(j)), &
                    real(GF(1)%correl(2, 2)%fctn(2, 2, j)), &
                    aimag(GF(1)%correl(2, 2)%fctn(2, 2, j))
            endif

         enddo


         !diagonal

         do i = 1, size(G(1)%correl(1, 2)%fctn(1, :, 1))
            do j = 1, size(G(1)%correl(1, 2)%fctn(1, 1, :))
               write(3000 + i, *) aimag(G(1)%correl(1, 2)%freq%vec(j)), &
                    real(G(1)%correl(1, 2)%fctn(i, i, j)), &
                    aimag(G(1)%correl(1, 2)%fctn(i, i, j))
            enddo
         enddo

         !gnambu normal

         do j = 1, size(GNAMBU%correl(1, 2)%fctn(1, 1, :))
            write(9002, *) aimag(GNAMBU%correl(1, 2)%freq%vec(j)), &
                 real(GNAMBU%correl(1, 2)%fctn(1, 2, j)), &
                 aimag(GNAMBU%correl(1, 2)%fctn(1, 2, j))
            write(9003, *) aimag(GNAMBU%correl(1, 2)%freq%vec(j)), &
                 real(GNAMBU%correl(1, 2)%fctn(1, 3, j)), &
                 aimag(GNAMBU%correl(1, 2)%fctn(1, 3, j))
         enddo

         !self energy normal

         do j = 1, size(SNAMBU%fctn(1, 1, :))
            write(5002, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(1, 2, &
                 j)), aimag(SNAMBU%fctn(1, 2, j))
            write(5003, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(2, 1, &
                 j)), aimag(SNAMBU%fctn(2, 1, j))
            write(5004, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(2, 2, &
                 j)), aimag(SNAMBU%fctn(2, 2, j))
            write(5005, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(2, 3, &
                 j)), aimag(SNAMBU%fctn(2, 3, j))
         enddo

         !self energy anomalous

         if(size(SNAMBU%fctn(1, :, 1)) > 5)then
            do j = 1, size(SNAMBU%fctn(1, 1, :))
               write(6001, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(1, &
                    5, j)), aimag(SNAMBU%fctn(1, 5, j))
               write(6005, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(5, &
                    1, j)), aimag(SNAMBU%fctn(5, 1, j))
               write(6002, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(1, &
                    6, j)), aimag(SNAMBU%fctn(1, 6, j))
               write(6003, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(2, &
                    5, j)), aimag(SNAMBU%fctn(2, 5, j))
               write(6004, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(2, &
                    6, j)), aimag(SNAMBU%fctn(2, 6, j))
            enddo
         endif

         !self energy anomalous

         if(size(SNAMBU%fctn(1, :, 1)) > 5)then
            do j = 1, size(SNAMBU%fctn(1, 1, :))
               write(6011, *) aimag(SNAMBU%freq%vec(j)), real(self_out(1, 5, &
                    j)), aimag(self_out(1, 5, j))
               write(6012, *) aimag(SNAMBU%freq%vec(j)), real(self_out(1, 6, &
                    j)), aimag(self_out(1, 6, j))
               write(6013, *) aimag(SNAMBU%freq%vec(j)), real(self_out(2, 5, &
                    j)), aimag(self_out(2, 5, j))
               write(6014, *) aimag(SNAMBU%freq%vec(j)), real(self_out(2, 6, &
                    j)), aimag(self_out(2, 6, j))
               write(6015, *) aimag(SNAMBU%freq%vec(j)), real(self_out(6, 1, &
                    j)), aimag(self_out(6, 1, j))
               write(6016, *) aimag(SNAMBU%freq%vec(j)), real(self_out(6, 2, &
                    j)), aimag(self_out(6, 2, j))
            enddo
         endif

         !**** done ****!

      endif

      call mpi_stop('done')

   end subroutine

   subroutine dump_to_fortb()

      use mpirout, only: mpi_stop

      implicit none

      integer :: j, i

      if(rank == 0)then

         do j = 1, size(G(1)%correl(2, 1)%fctn(1, 1, :))
            !normal
            write(1933, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                 real(G(1)%correl(2, 1)%fctn(1, 1, j)), aimag(G(1)%correl(2, &
                 1)%fctn(1, 1, j))
            write(2933, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                 real(G(1)%correl(2, 1)%fctn(1, 1, j)), aimag(G(1)%correl(2, &
                 1)%fctn(1, 1, j))
            write(9001, *) aimag(GNAMBU%correl(2, 1)%freq%vec(j)), &
                 real(GNAMBU%correl(2, 1)%fctn(1, 1, j)), &
                 aimag(GNAMBU%correl(2, 1)%fctn(1, 1, j))
            write(5001, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(1, 1, &
                 j)), aimag(SNAMBU%fctn(1, 1, j))
            write(6007, *) aimag(SNAMBU%freq%vec(j)), real(self_out(1, 1, j)), &
                 aimag(self_out(1, 1, j))
         enddo

         do j = 1, size(G(1)%correl(2, 1)%fctn(1, 1, :))
            write(1934, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                 real(G(1)%correl(2, 1)%fctn(1, 2, j)), aimag(G(1)%correl(2, &
                 1)%fctn(1, 2, j))
            write(1935, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                 real(G(1)%correl(2, 1)%fctn(2, 1, j)), aimag(G(1)%correl(2, &
                 1)%fctn(2, 1, j))
            write(1936, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                 real(G(1)%correl(2, 1)%fctn(2, 2, j)), aimag(G(1)%correl(2, &
                 1)%fctn(2, 2, j))
            write(1937, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                 real(G(1)%correl(2, 1)%fctn(1, 3, j)), aimag(G(1)%correl(2, &
                 1)%fctn(1, 3, j))
            write(1938, *) aimag(G(2)%correl(2, 1)%freq%vec(j)), &
                 real(G(2)%correl(2, 1)%fctn(1, 1, j)), aimag(G(2)%correl(2, &
                 1)%fctn(1, 1, j))
            write(2934, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                 real(G(1)%correl(2, 1)%fctn(1, 2, j)), aimag(G(1)%correl(2, &
                 1)%fctn(1, 2, j))
            write(2935, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                 real(G(1)%correl(2, 1)%fctn(2, 1, j)), aimag(G(1)%correl(2, &
                 1)%fctn(2, 1, j))
            write(2936, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                 real(G(1)%correl(2, 1)%fctn(2, 2, j)), aimag(G(1)%correl(2, &
                 1)%fctn(2, 2, j))

            !anomalous

            if(associated(GF(1)%correl(1, 1)%freq%vec)) then
               write(3933, *) aimag(GF(1)%correl(1, 1)%freq%vec(j)), &
                    real(GF(1)%correl(1, 1)%fctn(1, 1, j)), &
                    aimag(GF(1)%correl(1, 1)%fctn(1, 1, j))
               write(4933, *) aimag(GF(1)%correl(2, 2)%freq%vec(j)), &
                    real(GF(1)%correl(2, 2)%fctn(1, 1, j)), &
                    aimag(GF(1)%correl(2, 2)%fctn(1, 1, j))
               write(3934, *) aimag(GF(1)%correl(1, 1)%freq%vec(j)), &
                    real(GF(1)%correl(1, 1)%fctn(1, 2, j)), &
                    aimag(GF(1)%correl(1, 1)%fctn(1, 2, j))
               write(3935, *) aimag(GF(1)%correl(1, 1)%freq%vec(j)), &
                    real(GF(1)%correl(1, 1)%fctn(2, 1, j)), &
                    aimag(GF(1)%correl(1, 1)%fctn(2, 1, j))
               write(3936, *) aimag(GF(1)%correl(1, 1)%freq%vec(j)), &
                    real(GF(1)%correl(1, 1)%fctn(2, 2, j)), &
                    aimag(GF(1)%correl(1, 1)%fctn(2, 2, j))
               write(3937, *) aimag(GF(1)%correl(1, 1)%freq%vec(j)), &
                    real(GF(1)%correl(1, 1)%fctn(1, 3, j)), &
                    aimag(GF(1)%correl(1, 1)%fctn(1, 3, j))
               write(4934, *) aimag(GF(1)%correl(2, 2)%freq%vec(j)), &
                    real(GF(1)%correl(2, 2)%fctn(1, 2, j)), &
                    aimag(GF(1)%correl(2, 2)%fctn(1, 2, j))
               write(4935, *) aimag(GF(1)%correl(2, 2)%freq%vec(j)), &
                    real(GF(1)%correl(2, 2)%fctn(2, 1, j)), &
                    aimag(GF(1)%correl(2, 2)%fctn(2, 1, j))
               write(4936, *) aimag(GF(1)%correl(2, 2)%freq%vec(j)), &
                    real(GF(1)%correl(2, 2)%fctn(2, 2, j)), &
                    aimag(GF(1)%correl(2, 2)%fctn(2, 2, j))
            endif
         enddo

         !diagonal
         do i = 1, size(G(1)%correl(2, 1)%fctn(1, :, 1))
            do j = 1, size(G(1)%correl(2, 1)%fctn(1, 1, :))
               write(3000 + i, *) aimag(G(1)%correl(2, 1)%freq%vec(j)), &
                    real(G(1)%correl(2, 1)%fctn(i, i, j)), &
                    aimag(G(1)%correl(2, 1)%fctn(i, i, j))
            enddo
         enddo

         !gnambu normal
         do j = 1, size(GNAMBU%correl(2, 1)%fctn(1, 1, :))
            write(9002, *) aimag(GNAMBU%correl(2, 1)%freq%vec(j)), &
                 real(GNAMBU%correl(2, 1)%fctn(1, 2, j)), &
                 aimag(GNAMBU%correl(2, 1)%fctn(1, 2, j))
            write(9003, *) aimag(GNAMBU%correl(2, 1)%freq%vec(j)), &
                 real(GNAMBU%correl(2, 1)%fctn(1, 3, j)), &
                 aimag(GNAMBU%correl(2, 1)%fctn(1, 3, j))
         enddo

         !self energy normal

         do j = 1, size(SNAMBU%fctn(1, 1, :))
            write(5002, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(1, 2, &
                 j)), aimag(SNAMBU%fctn(1, 2, j))
            write(5003, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(2, 1, &
                 j)), aimag(SNAMBU%fctn(2, 1, j))
            write(5004, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(2, 2, &
                 j)), aimag(SNAMBU%fctn(2, 2, j))
            write(5005, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(2, 3, &
                 j)), aimag(SNAMBU%fctn(2, 3, j))
         enddo

         !self energy anomalous

         if(size(SNAMBU%fctn(1, :, 1)) > 5)then
            do j = 1, size(SNAMBU%fctn(1, 1, :))
               write(6001, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(1, &
                    5, j)), aimag(SNAMBU%fctn(1, 5, j))
               write(6002, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(1, &
                    6, j)), aimag(SNAMBU%fctn(1, 6, j))
               write(6003, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(2, &
                    5, j)), aimag(SNAMBU%fctn(2, 5, j))
               write(6004, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(2, &
                    6, j)), aimag(SNAMBU%fctn(2, 6, j))
               write(6005, *) aimag(SNAMBU%freq%vec(j)), real(SNAMBU%fctn(6, &
                    1, j)), aimag(SNAMBU%fctn(6, 1, j))
            enddo
         endif

         !self energy anomalous

         if(size(SNAMBU%fctn(1, :, 1)) > 5)then
            do j = 1, size(SNAMBU%fctn(1, 1, :))
               write(6008, *) aimag(SNAMBU%freq%vec(j)), real(self_out(1, 2, &
                    j)), aimag(self_out(1, 2, j))
               write(6009, *) aimag(SNAMBU%freq%vec(j)), real(self_out(1, 3, &
                    j)), aimag(self_out(1, 3, j))
               write(6010, *) aimag(SNAMBU%freq%vec(j)), real(self_out(2, 6, &
                    j)), aimag(self_out(2, 6, j))
               write(6011, *) aimag(SNAMBU%freq%vec(j)), real(self_out(1, 5, &
                    j)), aimag(self_out(1, 5, j))
               write(6012, *) aimag(SNAMBU%freq%vec(j)), real(self_out(1, 6, &
                    j)), aimag(self_out(1, 6, j))
               write(6013, *) aimag(SNAMBU%freq%vec(j)), real(self_out(2, 5, &
                    j)), aimag(self_out(2, 5, j))
               write(6014, *) aimag(SNAMBU%freq%vec(j)), real(self_out(6, 1, &
                    j)), aimag(self_out(6, 1, j))
            enddo
         endif

      endif

      call mpi_stop('done')

   end subroutine
