program dens_script
use linalg
use StringManip , only : StrInt2

implicit none

integer                  :: i,j,k,l,kk
character*100            :: filename,filename2
real(8)                  :: beta,df
real*8,allocatable       :: mat(:,:),sigmat(:,:),dens(:),rho(:),scattering(:),ZZ(:),ZZ_b(:),slope(:)
integer                  :: num,len__,status__,ll
character(300)           :: value(100)


  num=COMMAND_ARGUMENT_COUNT()
  write(*,*) 'NUMBER OF ARGUMENTS IN DENSITY SCRIPT : ', num
  do i=1,num
   call GET_COMMAND_ARGUMENT(i, value(i), len__, status__)
   write(*,*) 'ARGUMENT : ',i,TRIM(ADJUSTL(value(i)))
  enddo

  ll=StrInt2(value(3))
  
  write(*,*) ' ll = ' , ll

  if(num/=3.and.num/=4)then
   write(*,*) 'SCRIPT DENSITY IS EXPECTING EITHER 3 OR 4 ARGUMENTS'
   write(*,*) 'arg1,arg2 : green_output and sigma_output file names'
   write(*,*) ' ll       : number of orbitals, do not forget spin, so 2 orbitals with spin --> ll=4'
   write(*,*) ' dens     : [optional] density at the Fermi level obtained by script check_convergence'
   write(*,*) '            dens includes the projected overlap matrix, so dens=rho*Smat'
   write(*,*) '            the file name is : rho_fermi_level_from_gimp_spin'
   stop
  endif
   
  filename  =trim(adjustl(value(1)))
  filename2 =trim(adjustl(value(2)))
  
  write(*,*) 'name of files:'
  write(*,*) 'filename 1                  : ', filename
  write(*,*) 'equal to green_output       : ', trim(adjustl(filename))=='green_output'
  write(*,*) 'filename 2                  : ', filename2
  
  open(unit=1000,file=filename)
  open(unit=1001,file=filename2)
  
  l=0
  if( trim(adjustl(filename))=='green_output'  ) read(1000,*) !header
  do
   read(1000,*,end=66)
   l=l+1
  enddo
  66 continue
  
  write(*,*) 'allocating array, shapes : ', ll,l
  allocate(dens(ll),rho(ll),ZZ(ll),ZZ_b(ll),scattering(ll),slope(ll))
  allocate(sigmat(l,2*ll+1),mat(l,2*ll+1))
  
  write(*,*) 'reading green func'
  rewind(1000)
  if( trim(adjustl(filename))=='green_output'  ) read(1000,*)
  do i=1,l
    read(1000,*) mat(i,1),(mat(i,1+2*kk-1),mat(i,1+2*kk),kk=1,ll)
  enddo

  if(num==4)then
     write(*,*) 'READING DENSITY RHO IN FILE : ', trim(adjustl(value(4)))
     open(unit=5555,file=trim(adjustl(value(4))))
     do i=1,ll 
       read(5555,*) mat(1,1+2*i)
       read(5555,*) mat(2,1+2*i)
       mat(1,1+2*i)=mat(1,1+2*i)*(-dacos(-1.d0))
       mat(2,1+2*i)=mat(2,1+2*i)*(-dacos(-1.d0))
     enddo
     close(5555)
  endif

 
  write(*,*) 'reading sigma'
  if( trim(adjustl(filename))=='green_output'  ) read(1001,*)
  do i=1,l
   read(1001,*) sigmat(i,1),(sigmat(i,1+2*kk-1),sigmat(i,1+2*kk),kk=1,ll)
  enddo
  
  beta=dacos(-1.d0)/mat(1,1)
  
  write(*,*) 'temperature is       : ', 1.d0/beta
  write(*,*) 'beta is              : ',      beta
  
  write(*,*) 'SLOPE IMAGINARY PART : '
  do i=1,ll
   slope(i)=sigmat(1,1+2*i)/sigmat(1,1)
  enddo
  ZZ_b=max(1.d-7,1.d0/(1.d0-slope))
  where(ZZ_b>1.d0) ZZ_b=1.d-7
  write(*,*) 'SIGMA SLOPES : ', slope
  write(*,'(a,200f8.3)') 'Z obtained with p1-(0,0) : ', ZZ_b
  
  do i=1,ll
   slope(i)=(sigmat(2,1+2*i)-sigmat(1,1+2*i))/(sigmat(2,1)-sigmat(1,1))
  enddo
  do i=1,ll
   scattering(i)=sigmat(1,1+2*i) - slope(i)*sigmat(1,1)
  enddo
  ZZ=max(1.d-7,1.d0/(1.d0-slope))
  write(*,*) 'SIGMA SLOPES : ', slope
  
  where(ZZ>1.d0) ZZ=1.d-7
  write(*,'(a,200f8.3)') 'Z obtained with p3- : ', ZZ
  write(*,'(a,200f8.3)') 'Im Sigma            : ', scattering
  
  do i=1,ll
   slope(i)=(mat(2,1+2*i)-mat(1,1+2*i))/(mat(2,1)-mat(1,1))
  enddo
  do i=1,ll
   rho(i)=mat(1,1+2*i) - slope(i)*mat(1,1)
  enddo
   rho=-rho/dacos(-1.d0)
  
  write(*,'(a,200f8.3)')  'rho extrapolated to (0,0) : ', rho
  write(*,'(a,200f14.3)') 'C=rho/Z                   : ', rho/ZZ
  write(*,'(a,200f14.3)') 'C=rho/Z_b                 : ', rho/ZZ_b
  
  write(*,*) 'TOTAL RHO : ', sum(rho)
  write(*,*) 'TOTAL Z   : ', max(1.d-7,sum(rho)/sum(rho/ZZ))
  write(*,*) 'TOTAL Z_b : ', max(1.d-7,sum(rho)/sum(rho/ZZ_b))
  
  
  write(*,*) 'DENSITIES (with spin): '
  dens=0.d0
  do i=1,ll
     dens(i)=0.5+2.d0/beta*sum(mat(:,1+2*i-1))
   !Correction 1/omega^2
     df      =  mat(size(mat,1),1) * mat(size(mat,1),1+2*i-1) /dacos(-1.d0)
     dens(i) =  dens(i) + df
     write(*,*) 'w^2 correction : ', i,df 
  enddo
  write(*,'(a,200f8.5)') 'DENS          = ', 2.d0*    dens
  write(*,'(a,200f8.5)') 'TOTAL density = ', 2.d0*sum(dens)
  
  dens=0.
  do i=1,ll
   dens(i)=get_dens(mat(:,1),mat(:,1+2*i-1),mat(:,1+2*i),diag=.true.)
  enddo
  write(*,'(a,200f8.5)') 'DENS          = ', 2.d0*    dens
  write(*,'(a,200f8.5)') 'TOTAL density = ', 2.d0*sum(dens)
  
  close(1000)
  close(1001)



contains



!------------------------!
!------------------------!
!------------------------!
!------------------------!
!------------------------!

real(8) function get_dens(frequ,GlocRe,GlocIm,diag)
   use genvar, only: imi
   use linalg, only: mplx
   implicit none
   integer    :: k1,k2,i,j,k,nmatsu_frequ
   complex(8) :: dens
   real(8)    :: tau0,omega,df,frequ(:),GlocRe(:),GlocIm(:)
   complex(8) :: temp
   logical    :: diag
   nmatsu_frequ=size(frequ)
   tau0 = 1.d-8; dens = 0.d0
   do j=1,nmatsu_frequ
     omega   = frequ(j)
     if(abs(omega)<1.d-9) omega=1.d-9
     temp    = CMPLX(GlocRe(j),GlocIm(j),kind=8)
     if(diag)then
      dens = dens + 2.d0/beta * ( MPLX( omega * tau0 ) * (temp  + imi / omega )  )
     else
      dens = dens + 2.d0/beta * ( MPLX( omega * tau0 ) * (temp)  )
     endif
   enddo
   if(diag)then
    df      = GlocRe(nmatsu_frequ)*frequ(nmatsu_frequ)/dacos(-1.d0)
    dens    = dens + 0.5 + df
   endif
   get_dens = dens
return
end function

!------------------------!
!------------------------!
!------------------------!
!------------------------!
!------------------------!

end program


