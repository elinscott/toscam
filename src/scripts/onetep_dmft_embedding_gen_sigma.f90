!****************************!
!****************************!
!****************************!
!****************************!
!****************************!
!****************************!

module variableembed

   use namelistmod, only: namelist_set, namelist_init, putel_in_namelist, &
                               & look_for_namelist_in_file, look_for_command_line_argument

   integer  ::  remove_unconnected_sites
   real(8)  ::  cutoff_remove_sites

contains

   subroutine define_vars
      type(namelist_set) :: nm
      call namelist_init(nm, 200, name_of_namelist='gensigma')
 call putel_in_namelist(nm,remove_unconnected_sites,'remove_unconnected_sites',1,'do you want to remove disconnected sites? (enter 0 or 1)')
 call putel_in_namelist(nm,cutoff_remove_sites,'cutoff_remove_sites',0.0001d0,'enter cutoff on Hamiltonian and overlap (<cutoff is disconnected) (enter a real number)')
      call look_for_namelist_in_file(nm, './var.txt')
      call look_for_command_line_argument(nm)
   end subroutine

end module

!****************************!
!****************************!
!****************************!
!****************************!

program embed_self
   use variableembed
   implicit none
   integer                 ::  i, j, ii, jj, s, t, kkk1, kkk2, jjj, k
   complex(8), allocatable  :: selfb(:, :, :), self(:, :, :), energy(:), tt(:, :), ttt(:, :), GG(:, :)
   real(8), allocatable     ::  ovH1(:, :), ovH0(:, :), H0(:, :), H1(:, :)
   integer, allocatable     ::  index_device(:), orb_num_device(:)
   real(8)                 ::  mmu_lead
   integer                 ::  ien, pub_dmft_points, size_mat0
   real(8)                 ::  pub_dmft_temp, mmu
   integer                 :: totsites, iii, totorb
   real(8)                 :: test_matmul1(10), test_matmul2(10)
   real(8), allocatable     :: S0b(:, :), S1b(:, :), H0b(:, :), H1b(:, :)
   integer, allocatable     :: ib(:), mapatom(:), maporb(:)
   real(8), allocatable     :: muarray(:)

   call define_vars

   open (unit=1001, file='embedding_potentials', form='unformatted')
   open (unit=1002, file='_H0_', form='unformatted')
   open (unit=1003, file='_H1_', form='unformatted')
   open (unit=1004, file='sigma_embedding', form='unformatted')
   open (unit=1005, file='connections_D_device_sites')

!    sigma_output_embed : [input]
!                ien , pub_dmft_points , pub_dmft_temp , size_mat0 , size_mat0
!                fermi_e(is)+pub_dmft_chem_shift,energy,ttt0(1:size_mat0,1:size_mat0)
!
!    _H0_ (_H1_) : [input]
!                write(ufi) kkk
!                write(ufi) tmpH
!                write(ufi) ovH
!
!    sigma_embedding : [output]
!            first line only : size_mat0 (size(vec))
!                  vec(:)
!                self(:,:)
!
!    connections_D_(L,R)_device_sites
!            first line : number entry
!              i1 (first entry)
!              i2 (second entry) etcetera...
!
!    reading : connections_D_(L,R)_device_sites
!
   read (1005, *) i
   allocate (index_device(i), orb_num_device(i))
   do i = 1, size(index_device)
      read (1005, *) index_device(i), orb_num_device(i)
   enddo

!read _H0_
   read (1002) kkk1
   allocate (H0(kkk1, kkk1), ovH0(kkk1, kkk1))
   read (1002) H0
   read (1002) ovH0

   write (*, *) '=========================='
   write (*, *) 'H0   : ', H0(1:2, 1:2)
   write (*, *) '=========================='
   write (*, *) 'ovH0 : ', ovH0(1:2, 1:2)
   write (*, *) '=========================='

!read _H1_
   read (1003) kkk2
   if (kkk1 /= kkk2) then
      write (*, *) 'ERROR : dimensions dont match in _H0_ and _H1_'
      stop
   endif
   allocate (H1(kkk2, kkk2), ovH1(kkk2, kkk2))
   read (1003) H1
   read (1003) ovH1

   if (sum(orb_num_device) /= kkk2) then
      write (*, *) 'ERROR : sum of number of orbitals for device sites does'
      write (*, *) 'NOT match the size of the matrix H0 and H1'
      stop
   endif

   write (*, *) '=========================='
   write (*, *) 'lines'
   write (*, *) 'H0   max  : ', maxval(abs(H1))
   write (*, *) 'ovH1 max  : ', maxval(abs(ovH1))
   write (*, *) '=========================='

!sigma_output_embed
   read (1001) ien, pub_dmft_points, pub_dmft_temp, size_mat0, size_mat0
   allocate (self(size_mat0, size_mat0, pub_dmft_points), energy(pub_dmft_points))
   allocate (tt(size_mat0, size_mat0))
   allocate (ttt(size_mat0, size_mat0))
   allocate (GG(size_mat0, size_mat0))

   read (1001) mmu_lead, energy(1), self(:, :, 1)
   do i = 2, pub_dmft_points
      read (1001) ien, pub_dmft_points, pub_dmft_temp, size_mat0, size_mat0
      read (1001) mmu_lead, energy(i), self(:, :, i)
   enddo

   if (remove_unconnected_sites == 1) then
      totsites = 0; totorb = 0
      i = 0
      do iii = 1, size(orb_num_device)
         if (iii > 1) i = i + orb_num_device(iii - 1)
         j = orb_num_device(iii)
         if (maxval(abs(H1(i + 1:i + j, :))) < cutoff_remove_sites .and. maxval(abs(H1(:, i + 1:i + j))) < cutoff_remove_sites&
 & .and. maxval(abs(ovH1(i + 1:i + j, :))) < cutoff_remove_sites .and. maxval(abs(ovH1(:, i + 1:i + j))) < cutoff_remove_sites) then
            write (*, *) ' ATOM NOT CONNECTED (position in the device list): ', iii
         else
            totorb = totorb + orb_num_device(iii)
            totsites = totsites + 1
         endif
      enddo

      write (*, *) 'THERE ARE IN TOTAL [x] ATOMS : ', totsites
      write (*, *) 'AND IN TOTAL [x] ORBITALS    : ', totorb
      i = totorb
      test_matmul1(1) = sum(matmul(H1, H0))
      test_matmul1(2) = sum(matmul(ovH1, H0))
      test_matmul1(3) = sum(matmul(H1, ovH0))
      test_matmul1(4) = sum(matmul(ovH1, ovH0))
      allocate (ib(totsites), mapatom(totsites))
      allocate (H0b(i, i), S0b(i, i), H1b(i, i), S1b(i, i), maporb(i), selfb(i, i, pub_dmft_points))

      j = 0; i = 0; maporb = 0; mapatom = 0; jjj = 0
      do iii = 1, size(orb_num_device)
         k = orb_num_device(iii)
         if (maxval(abs(H1(i + 1:i + k, :))) < cutoff_remove_sites .and. maxval(abs(H1(:, i + 1:i + k))) < cutoff_remove_sites &
  &.and. maxval(abs(ovH1(i + 1:i + k, :))) < cutoff_remove_sites .and. maxval(abs(ovH1(:, i + 1:i + k))) < cutoff_remove_sites) then
         else
            do k = 1, orb_num_device(iii)
               maporb(j + k) = i + k
            enddo
            jjj = jjj + 1
            mapatom(jjj) = iii
            j = j + orb_num_device(iii)
         endif
         i = i + orb_num_device(iii)
      enddo

      do i = 1, totorb
      do j = 1, totorb
         H1b(i, j) = H1(maporb(i), maporb(j))
         H0b(i, j) = H0(maporb(i), maporb(j))
         S0b(i, j) = ovH0(maporb(i), maporb(j))
         S1b(i, j) = ovH1(maporb(i), maporb(j))
         selfb(i, j, :) = self(maporb(i), maporb(j), :)
      enddo
      enddo

      do i = 1, totsites
         ib(i) = index_device(mapatom(i))
      enddo

      deallocate (H0, H1, ovH0, ovH1, index_device, self)
      i = totorb
      allocate (H0(i, i), ovH0(i, i), ovH1(i, i), H1(i, i), index_device(totsites), self(i, i, pub_dmft_points))
      H0 = H0b; H1 = H1b; ovH0 = S0b; ovH1 = S1b; self = selfb
      index_device = ib
      test_matmul2(1) = sum(matmul(H1, H0))
      test_matmul2(2) = sum(matmul(ovH1, H0))
      test_matmul2(3) = sum(matmul(H1, ovH0))
      test_matmul2(4) = sum(matmul(ovH1, ovH0))
      do i = 1, 4
         write (*, *) 'TESTING SITE REMOVAL (should be the same) : ', test_matmul1(i), test_matmul2(i)
      enddo
      write (*, *) 'SIZE OF MATRIX CHANGED FROM [x] TO [y] : ', size_mat0, totorb
      size_mat0 = totorb
   endif

   write (*, *) 'Write in file sigma_embedding, size matrix : ', size_mat0

   write (1004) size_mat0, size(index_device)

   write (*, *) 'Write index_device mask array, shape : ', shape(index_device)
   write (*, *) 'shape H0 H1 : ', shape(H0), shape(H1)
   write (*, *) 'shape ovH0 ovH1 : ', shape(ovH0), shape(ovH1)
   write (*, *) 'shape self : ', shape(self)
   write (*, *) 'energy 1 : ', energy(1)

   allocate (muarray(pub_dmft_points))
   muarray = mmu_lead

   write (1004) index_device, H0, H1, ovH0, ovH1, self, muarray, energy

!do i=1,pub_dmft_points
! write(1004) mmu_lead,self(:,:,i),energy(i)
! if(i==pub_dmft_points-1)then
!  write(*,*) 'SELF OBTAINED AT w=oo, w=', energy(i)
!  write(*,*) 'SELF = ', maxval(abs(self(:,:,i)))
!  write(*,*) 'chem pot = ', mmu_lead
! endif
!enddo

   close (1001)
   close (1002)
   close (1003)
   close (1004)
   close (1005)

end program

!****************************!
!****************************!
!****************************!
!****************************!
!****************************!
!****************************!
!****************************!
!****************************!
!****************************!
!****************************!

