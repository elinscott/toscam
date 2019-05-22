module quantum_algebra

   implicit none

   private
   !-----------------------------------------------------------------!
   TYPE fermion
      integer    :: nstate !2**state
      integer    :: ifull, k1, kk1, k2
      complex(8) :: coef(40000)
      integer    :: state(40000)
   END TYPE
   !-----------------------------------------------------------------!

   type(fermion), private :: BCS_STATE

   INTERFACE ASSIGNMENT(=)
      MODULE PROCEDURE equal_fermion
   END INTERFACE

   public :: fermion

contains

   subroutine c_destroy(f1, f2, aa, m, sigma)

      ! C_m_sigma * |.....1(m)..... > /= 0

      implicit none

      type(fermion), intent(inout) :: f1, f2
      complex(8) :: aa
      integer    :: k, m, mm, ket, sigma, isum, tot, i

      if (sigma == 1) k = m
      if (sigma == -1) k = m + f1%k1
      tot = f1%ifull

      do i = 1, tot
         ket = f1%state(i)
         if (ket > 0 .and. ibits(ket, k - 1, 1) == 1) then
            f2%ifull = f2%ifull + 1
            f2%state(f2%ifull) = ibclr(ket, k - 1)
            !signe fermionique:
            isum = 0
            do mm = 1, k - 1
               if (ibits(ket, mm - 1, 1) == 1) isum = isum + 1
            enddo
            f2%coef(f2%ifull) = f1%coef(i)*aa*((-1)**isum)
         endif
      enddo

      return
   end subroutine

   subroutine c_create(f1, f2, aa, m, sigma)

      ! C_m_sigma * |.....1(m)..... > /= 0

      implicit none

      type(fermion), intent(inout) :: f1, f2
      complex(8) :: aa
      integer    :: k, m, mm, ket, sigma, isum, tot, i

      if (sigma == 1) k = m
      if (sigma == -1) k = m + f1%k1

      tot = f1%ifull

      do i = 1, tot
         ket = f1%state(i)
         if (ket > 0 .and. ibits(ket, k - 1, 1) == 0) then
            f2%ifull = f2%ifull + 1
            f2%state(f2%ifull) = ibset(ket, k - 1)
            !signe fermionique:
            isum = 0
            do mm = 1, k - 1
               if (ibits(ket, mm - 1, 1) == 1) isum = isum + 1
            enddo
            f2%coef(f2%ifull) = f1%coef(i)*aa*((-1)**isum)
         endif
      enddo

      return
   end subroutine

   subroutine c_cicj(f1, f2, aa, m2, sigma2, m, sigma)

      implicit none

      type(fermion), intent(inout) :: f1, f2
      complex(8) :: aa
      integer    :: k, m, mm, ket, sigma, isum, tot, i
      integer    :: m2, sigma2, k2

      if (sigma == 1) k = m
      if (sigma == -1) k = m + f1%k1
      if (sigma2 == 1) k2 = m2
      if (sigma2 == -1) k2 = m2 + f1%k1

      tot = f1%ifull

      do i = 1, tot
         ket = f1%state(i)
         if (ket > 0 .and. ibits(ket, k - 1, 1) == 1) then
            ket = ibclr(ket, k - 1)
            isum = 0
            do mm = 1, k - 1
               if (ibits(ket, mm - 1, 1) == 1) isum = isum + 1
            enddo
            if (ibits(ket, k2 - 1, 1) == 0) then
               f2%ifull = f2%ifull + 1
               ket = ibclr(ket, k2 - 1)
               f2%state(f2%ifull) = ket
               do mm = 1, k2 - 1
                  if (ibits(ket, mm - 1, 1) == 1) isum = isum + 1
               enddo
               f2%coef(f2%ifull) = f1%coef(i)*aa*((-1)**isum)
            endif
         endif
      enddo

      return
   end subroutine

   subroutine apply_BCS_op(vv, f1)

      implicit none

      type(fermion), intent(inout) :: f1
      type(fermion) :: f2
      complex(8)    :: vv(:)
      integer       :: m

      call reset_fermion(f2, f1%nstate)
      do m = f1%k1 + 1, f1%nstate
         call c_destroy(f1, f2, vv(m), m - f1%k1, -1)
      enddo
      do m = 1, f1%k1
         call c_create(f1, f2, vv(m), m, 1)
      enddo
      call sort_Fock_space(f1, f2)
      return
   end subroutine

   subroutine mat_op_part_trou(vv, f1)

      implicit none

      type(fermion), intent(inout) :: f1
      type(fermion) :: f2
      complex(8)    :: vv(:, :, :, :)
      integer       :: m, i, j, siz1, sigma, sigma2, m2

      call reset_fermion(f2, f1%nstate)
      siz1 = size(vv(:, 1, 1, 1))
      do i = 1, siz1
         do j = 1, siz1
            do sigma = -1, 1, 2
               do sigma2 = -1, 1, 2
                  if (abs(vv(i, j, sigma, sigma2)) > 1.d-8) call c_cicj(f1, f2, &
                                                                        vv(i, j, sigma, sigma2), i, sigma, j, sigma2)
               enddo
            enddo
         enddo
      enddo
      call sort_Fock_space(f1, f2)
      return
   end subroutine

   subroutine apply_vec(vv, f1)

      implicit none

      type(fermion), intent(inout) :: f1
      type(fermion) :: f2
      complex(8)    :: vv(:, :)
      integer       :: m, i, j, siz1, sigma

      siz1 = size(vv(:, 1))
      call reset_fermion(f1, f1%nstate)
      call set_to_part_vacuum(f1)

      call reset_fermion(f2, f1%nstate)
      do i = 1, siz1
         do sigma = -1, 1, 2
            call c_create(f1, f2, vv(i, sigma), i, sigma)
         enddo
      enddo
      call sort_Fock_space(f1, f2)

      return
   end subroutine

   subroutine apply_vec_left(vv, f1)

      implicit none

      type(fermion), intent(inout) :: f1
      type(fermion) :: f2
      complex(8)    :: vv(:, :)
      integer       :: m, i, j, siz1, sigma

      siz1 = size(vv(:, 1))
      call reset_fermion(f1, f1%nstate)
      call set_to_part_vacuum(f1)

      call reset_fermion(f2, f1%nstate)
      do i = 1, siz1
         do sigma = -1, 1, 2
            call c_create(f1, f2, vv(i, sigma), i, sigma)
         enddo
      enddo
      call sort_Fock_space(f1, f2)

      return
   end subroutine

   subroutine scal_prod(vec, mat)

      implicit none

      type(fermion) :: f1
      type(fermion) :: f2
      complex(8)    :: mat(:, :, :, :), vec(:, :)
      integer       :: m, i, j, sigma, sigma2

      call reset_fermion(f1, f1%nstate)

      call set_to_part_vacuum(f1)
      call apply_vec(vec, f1)
      call mat_op_part_trou(mat, f1)

      call apply_vec_left(vec, f1)

      return
   end subroutine

   subroutine sort_Fock_space(f1, f2)

      implicit none

      type(fermion), intent(inout) :: f1, f2
      integer :: i, j, k, l, ini, ini2

      do i = 1, f2%ifull
         if (f2%state(i) >= 0) then
            ini = f2%state(i)
            do j = i + 1, f2%ifull
               if (f2%state(j) >= 0) then
                  ini2 = f2%state(j)
                  if (ini2 == ini) then
                     f2%coef(i) = f2%coef(i) + f2%coef(j)
                     f2%coef(j) = 0.d0
                     f2%state(j) = -1
                  endif
               endif
            enddo
         endif
      enddo

      call reset_fermion(f1, f1%nstate)
      do i = 1, f2%ifull
         if (f2%state(i) >= 0) then
            f1%ifull = f1%ifull + 1
            f1%state(f1%ifull) = f2%state(i)
            f1%coef(f1%ifull) = f2%coef(i)
         endif
      enddo

      return
   end subroutine

   subroutine equal_fermion(f1, f2)

      implicit none

      type(fermion), intent(inout) :: f1
      type(fermion), intent(in)    :: f2
      integer :: imax

      imax = max(f1%ifull, f2%ifull)
      f1%k1 = f2%k1
      f1%kk1 = f2%kk1
      f1%k2 = f2%k2
      f1%coef(1:imax) = f2%coef(1:imax)
      f1%state(1:imax) = f2%state(1:imax)
      f1%nstate = f2%nstate
      f1%ifull = f2%ifull
   end subroutine

   subroutine reset_fermion(fermion_in, kk)

      implicit none

      type(fermion), intent(inout) :: fermion_in
      integer :: kk, k1, kk1, k2

      fermion_in%k1 = kk/2
      fermion_in%kk1 = kk/2 + 1
      fermion_in%k2 = kk - 1
      fermion_in%coef = 1.
      fermion_in%state = -1
      fermion_in%nstate = kk
      fermion_in%ifull = 0
   end subroutine

   subroutine set_to_vacuum_BCS(fermion_in)

      implicit none

      type(fermion), intent(inout) :: fermion_in
      integer :: i, k1, k2

      k1 = fermion_in%k1
      k2 = fermion_in%nstate
      !up: 1....nstate/2,
      !dn: ns/2 + 1......nstate
      fermion_in%ifull = 1
      fermion_in%coef(1) = 1.
      fermion_in%state = -1
      fermion_in%state(1) = 0
      do i = k1 + 1, k2
         fermion_in%state(1) = ibset(fermion_in%state(1), i - 1)
      enddo

   end subroutine

   subroutine set_to_part_vacuum(fermion_in)

      implicit none

      type(fermion), intent(inout) :: fermion_in
      integer :: i

      fermion_in%coef = 0.
      fermion_in%ifull = 1
      fermion_in%coef(1) = 1.
      fermion_in%state = -1
      fermion_in%state(1) = 0
      return
   end subroutine

   subroutine check_vec_vacuum(vv, fermion_in)

      implicit none

      type(fermion), intent(inout) :: fermion_in
      integer    :: i, j, k, k1, k2
      complex(8) :: vv(:)

      if (size(vv) /= fermion_in%nstate) stop 'error bad vacuum size'
      k = 0
      do j = 1, fermion_in%nstate
         if (ABS(vv(j)) > 1.d-14) k = k + 1
      enddo
      if (k == 0) then
         write (*, *) 'error bad vacuum'
         stop
      endif
      return
   end subroutine

   subroutine check_normalization(rr)

      implicit none

      complex(8) :: rr
      integer    :: i

      if (ABS(rr) < 1.d-20) then
         write (*, *) 'error null coeff, dnormi : ', ABS(rr)
         stop 'error check_normalization'
      endif
      return
   end subroutine

   subroutine check_spectrum(i, W)

      use genvar, only: rank

      implicit none

      integer :: i, l
      real(8) :: W(:)

      if (W(i) > 1.d-1) then
         if (rank == 0) then
            write (27, *) 'eigenvalue system problem routine Mqcoef : positive &
                 &eigenvalue'
            do l = 1, size(W)
               write (27, *) l, W(l)
            enddo
         endif
      endif
      return
   end subroutine

   subroutine BCS_operator(siz, B, W, icoef, coef)

      use genvar, only: messages

      implicit none

      integer    :: siz ! = q1**2
      integer    :: kk2, i, l
      complex(8) :: dnormi
      real(8)    :: W(:)
      complex(8) :: B(:, :)
      integer    :: icoef(siz, 2)
      complex(8) :: coef(siz)
      integer    :: ivecout(siz, 2)
      complex(8) :: vecout(siz)

      dnormi = 100000000000000000.
      call reset_fermion(BCS_STATE, size(W))
      call set_to_vacuum_BCS(BCS_STATE)

      do i = BCS_STATE%k1, BCS_STATE%nstate
         call check_vec_vacuum(B(i, :), BCS_STATE)
      enddo

      do i = 1, BCS_STATE%k1
         call check_spectrum(i, W)
         call d_k(BCS_STATE, i, B(:, i), vecout, ivecout, dnormi)
      enddo

      call check_normalization(dnormi)
      if (messages) write (40, *) 'dnormi : ', ABS(dnormi)
      coef = vecout/dnormi
      icoef = ivecout

      return
   end subroutine

   subroutine d_k(f1, iter, vec, vecout, ivecout, dnormi)

      use linalg, only: int_to_array

      implicit none

      type(fermion), intent(inout) :: f1
      complex(8) :: vec(:), vecout(:), dnormi
      integer    :: k, i, j, l, m, ini, isum, mm, kk, kkk, iii, iiii(2)
      integer    :: ivecout(:, :), iter

      call apply_BCS_op(vec, f1)

      if (iter /= f1%k1) return

      ! OUTPUT :
      kk = 1
      vecout = 0.
      ivecout = 0

      do i = 1, f1%ifull

         l = sum(int_to_array(f1%nstate, f1%state(i)))

         if (l == 0 .and. iter == f1%k1) then
            dnormi = f1%coef(i)
         endif
         if (mod(l, 2) == 1 .and. iter == f1%k1) then
            write (*, *) 'error subroutine MQcoef'
            write (*, *) 'sum  = ', l
            write (*, *) 'f1%state(i): ', f1%state(i)
            write (*, *) 'array : ', int_to_array(f1%nstate, f1%state(i))
            write (*, *)
            stop
         endif

         if (l == 2) then
            kkk = 1
            do iii = 1, f1%nstate
               if (f1%state(i) > 0 .and. ibits(f1%state(i), iii - 1, 1) == 1) then
                  iiii(kkk) = iii
                  kkk = kkk + 1
               endif
            enddo
            ivecout(kk, 1) = iiii(1)
            ivecout(kk, 2) = iiii(2)
            vecout(kk) = f1%coef(i)
            kk = kk + 1
         endif

      enddo

      if (iter == f1%k1 .and. kk /= f1%k1**2 + 1) then
         write (*, *) 'error routine Mqcoef, kk = /.. : ', kk - 1, f1%k1**2
         write (*, *) 'iter : ', iter
         stop
      endif

   end subroutine

end module
