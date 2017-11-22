program opticaltensor

 use linalg
 use StringManip,           only: StrInt2, toString
 use matrix,                only: diag
 use mesh
 !use plotlib
 use init_and_close_my_sim, only: initialize_my_simulation,finalize_my_simulation

implicit none

integer             :: i,ijk,i1,i2,ii,nspin
real(8)             :: alpha_iso,alpha_non_iso,alpha2(3,3),alpha(3,3),tra,tra2
real(8),allocatable :: sigma(:,:,:),frequ(:)
integer             :: k1,k2,nfrequ

call system("ls 111 && echo 1 > scratch_file ")
call system("ls 211 && echo 2 > scratch_file ")

open(unit=100,file='scratch_file')
read(100,*) nspin
write(*,*) 'THERE ARE [x] SPIN SPECIFIES : ', nspin

do ii=1,nspin
 do k1=1,3
  do k2=1,3
     write(*,*) 'opening file name : ',trim(adjustl(toString(ii)))//trim(adjustl(toString(k1)))//trim(adjustl(toString(k2))) 
     open(unit=101,file=trim(adjustl(toString(ii)))//trim(adjustl(toString(k1)))//trim(adjustl(toString(k2))) )
     read(101,*)
     read(101,*)
     ijk=0
     do 
       read(101,*,end=145)
       ijk=ijk+1
     enddo
 145 continue
     rewind(101)
     read(101,*)
     read(101,*)
     nfrequ=ijk
     if(.not.allocated(sigma))then
       allocate(sigma(3,3,nfrequ),frequ(nfrequ))
     endif
     write(*,*) ' there are [x] frequencies : ', ijk
     do i=1,nfrequ
      read(101,*) frequ(i),sigma(k1,k2,i)
     enddo
    close(101)
  enddo
 enddo
enddo

close(100)

 call system("rm scratch_file")

 open(unit=1000,file='optic_isotropic')
 open(unit=1001,file='optic_anisotropic')

 do i=1,nfrequ
       alpha  =        sigma(:,:,i)
       alpha2 = matmul(sigma(:,:,i),sigma(:,:,i))
       tra =sum(diag(alpha ))
       tra2=sum(diag(alpha2))
       alpha_iso = tra/3.d0
       if( 3.d0/2.d0*tra2 - 1.d0/2.d0 * ( tra )**2.d0<0.d0) then
         write(*,*) 'error negative optical tensor'
         stop
       endif
       alpha_non_iso = sqrt( 3.d0/2.d0*tra2 - 1.d0/2.d0 * ( tra )**2.d0 )
       write(1000,*) frequ(i),alpha_iso
       write(1001,*) frequ(i),alpha_non_iso
 enddo

 close(1000)
 close(1001)

end program


