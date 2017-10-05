! -------------------------------------------------------------------
!
! This is a program for generating Coulomb potentials in reciprocal 
! space in the CASTEP foram.
! 
! Written by Chris-Kriton Skylaris on 20/4/2001
! 
!  TCM, Cavendish Laboratory
!  Madingley Road
!  Cambridge CB3 0HE
!  UK
! 
!----------------------------------------------------------------------


program coulomb_potential


  implicit none

  INTEGER, PARAMETER :: DP  = KIND(1.0d0)

  REAL(kind=DP), PARAMETER :: PI=3.141592653589793238462643383279502884197_DP

  real(kind=DP), parameter :: HARTREE_IN_EVS=27.2116529_DP ! cks: The value of one Hartree in terms of electron-volts

  REAL(kind=DP), PARAMETER :: ANGSTROM=1.889726313_DP ! THE VALUE OF ONE ANGSTROM IN TERMS OF BOHRS


  integer :: ios1, ios2, ionic_charge, line

  character(len=128) :: file_name

  real(kind=DP) :: point, step, four_pi, factor




  ! cks :: set the type of element and its charge
  file_name="H.realpot"
  ionic_charge=1

  open(4,file=file_name,iostat=ios1)

  if (ios1.ne.0) print*,'error opening file',file_name,' in services_c2_form_density, iostat=',ios1


  write(4,'(a)') "ION_CHARGE"
  write(4,'(i3)') ionic_charge 
  write(4,'(a)') "END_ION_CHARGE"

  write(4,'(a)') " "

  write(4,'(a)') "NONLOCAL_MOM_KB_RAD"
  write(4,'(a)') "END_NONLOCAL_MOM_KB_RAD"

  write(4,'(a)') " "
  write(4,'(a)') " "

  write(4,'(a)') "V_LOCAL"


  ! cks : it will create 2001 frequency components corresponding to 2001 
  ! cks : equally spaced frequencies from 0 to 100 reciprocal angstroms.
  point=0.0_DP
  step=100.0_DP/2000_DP 

   
  four_pi=4.0_DP*PI
  factor=four_pi*HARTREE_IN_EVS/ANGSTROM

  do line=1,667

     if (line.eq.1) then
        write(4,'(3F21.10)') 0.0_DP, -factor/(point+step)**2, -factor/(point+2.0_DP*step)**2
     else
        write(4,'(3F21.10)') -factor/point**2, -factor/(point+step)**2, -factor/(point+2.0_DP*step)**2
     endif
     
     point=point+3.0_DP*step

  enddo




  write(4,'(a)') "END_V_LOCAL"


  close(4,iostat=ios2)

  if (ios2.ne.0) print*,'error closing file',file_name,'in services_c2_form_density, iostat=',ios2

 



end program coulomb_potential

