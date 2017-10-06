program summary_script
use linalg
use StringManip , only : StrInt2
use mesh
use init_and_close_my_sim   , only : initialize_my_simulation,finalize_my_simulation

implicit none

integer                   :: i,j,k,l,kk
real(8)                   :: beta,dd1,dd2
real(8),allocatable       :: Rhomat(:),Zmat(:),Nmat(:),xt(:),err(:)
integer,allocatable       :: histoRho(:),histoZ(:),histoN(:)
integer                   :: nn,ll,nbin
character(300)            :: value(100)
character(22),allocatable :: label(:)
integer                   :: norb,niter_dmft,num,len__,status__

 !-------------------------------!
   call initialize_my_simulation
 !-------------------------------!

 num=COMMAND_ARGUMENT_COUNT()

 if(num/=3)then
  write(*,*) 'wrong number of arguments'
  write(*,*) 'NUM       : ', num
  write(*,*) 'ARGUMENTS : ', value(1:num)
  write(*,*) 'execting 3 arguments (niter_dmft and nbin and number of orbitals) '
  stop
 endif
 do i=1,num
   call GET_COMMAND_ARGUMENT(i, value(i), len__, status__)
 enddo
 niter_dmft=StrInt2(value(1))
 nbin=StrInt2(value(2))
 norb=StrInt2(value(3))

 write(*,*) 'DMFT ITER = ', niter_dmft
 write(*,*) 'nbin      = ', nbin
 write(*,*) 'norb      = ', norb

 call system(" onetep.dmft.collectZ_in_file "// trim(adjustl(toString(niter_dmft)))//" "// trim(adjustl(toString(norb))) )

   open(unit=1010,file="histogram_Z")
   open(unit=1011,file="histogram_Nr")
   open(unit=1012,file="histogram_Rho")

   ll=0
   do 
    read(1010,*,end=85) 
    ll=ll+1
   enddo
   85 continue
   write(*,*) 'THERE ARE [x] ELEMENTS IN FILE : ', ll
   rewind(1010)
   allocate(Rhomat(ll),Zmat(ll),Nmat(ll),histoRho(nbin),histoZ(nbin),histoN(nbin),xt(nbin),label(nbin),err(nbin))
   do i=1,ll
    read(1010,*) Zmat(i)
    read(1011,*) Nmat(i)
    read(1012,*) Rhomat(i)
   enddo
   close(1010)
   close(1011)
   close(1012)

   open(unit=2323,file='summary_run')

   write(2323,*) " CORRELATED ATOMS TOTAL Z  : " ,  sum(Rhomat)  /  sum( Rhomat/Zmat ) 
   write(2323,*) " CORRELATED ATOMS TOTAL N  : " ,  sum(Nmat)    /  dble(ll)  
   write(2323,*) " CORRELATED ATOMS TOTAL Nf : " ,  sum(Rhomat)  /  dble(ll)

   dd1=minval(Rhomat)-1.d0
   dd2=maxval(Rhomat)+1.d0
   call bin_function(Rhomat,histoRho,nbin,dd1=dd1,dd2=dd2)
   call build1Dmesh(xt,nbin,dd1,dd2)
   do i=1,nbin
    label(i)=trim(adjustl(toString(xt(i))))
   enddo
   err=0.
#ifdef _plot
   call simplehist_(real(histoRho),nbin,0.,0.,label,real(err),'histoRho','x','y')
#endif
   open(unit=2424,file='histRho_graph')
   do i=1,nbin
    write(2323,*) 'Rho binning : ', label(i), histoRho(i)
    write(2424,*) xt(i),histoRho(i)
   enddo
   close(2424)

   dd1=0.d0
   dd2=1.d0
   call bin_function(Zmat,histoZ,nbin,dd1=dd1,dd2=dd2)
   call build1Dmesh(xt,nbin,dd1,dd2)
   do i=1,nbin
    label(i)=trim(adjustl(toString(xt(i))))
   enddo
   err=0.
#ifdef _plot
   call simplehist_(real(histoZ),nbin,0.,0.,label,real(err),'histoZ','x','y')
#endif
   open(unit=2424,file='histZ_graph')
   do i=1,nbin
    write(2323,*) 'Z binning : ', label(i), histoZ(i)
    write(2424,*) xt(i),histoZ(i)
   enddo
   close(2424)

   dd1=0.d0
   dd2=7.d0
   call bin_function(Nmat,histoN,nbin,dd1=dd1,dd2=dd2)
   call build1Dmesh(xt,nbin,dd1,dd2)    
   do i=1,nbin
    label(i)=trim(adjustl(toString(xt(i))))
   enddo
   err=0.
#ifdef _plot
   call simplehist_(real(histoN),nbin,0.,0.,label,real(err),'histoN','x','y')
   call simplehist_(real(histoN),nbin,0.,0.,label,real(err),'histoN','x','y')
#endif
   open(unit=2424,file='histN_graph')
   do i=1,nbin
    write(2323,*) 'N binning : ', label(i), histoN(i)
    write(2424,*) xt(i),histoN(i)
   enddo
   close(2424)
   close(2323)

 !-------------------------------!
   call finalize_my_simulation
 !-------------------------------!

end program
