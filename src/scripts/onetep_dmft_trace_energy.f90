program energies
implicit none
real(8) :: tmp,energy(1000),beta
integer :: i,j,k


open(unit=30,file='energies.txt')

read(30,*) beta
j=0
do 
 read(30,*,end=2121) tmp
 j=j+1
 energy(j)=tmp
 write(*,*) 'ENERGY : ', energy(j)
enddo
2121 continue

close(30)

write(*,*) 'Boltzman average <E> : ', sum( energy(1:j) * exp(-beta*energy(1:j)) ) / sum(exp(-beta*energy(1:j)))

end program
