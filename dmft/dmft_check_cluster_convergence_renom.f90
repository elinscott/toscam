

!======================!
!======================!

program dmft_check_convergence
!---------------!
  use DMFT_SOLVER_ED
  !use dmftmod
  use genvar
  use init_and_close_my_sim
  use fortran_cuda
  use matrix
  use StringManip
  !use DLPLOT
  !use plotlib
  !use plot2d
!---------------!
implicit none

integer                :: fac,paramagnet
real(8),allocatable    :: temp(:),vvv(:,:)
integer                :: n_frequ_long,i1,i2,kk_,i,j,k,l,channels,n_frequ
real(8),allocatable    :: rotmat(:,:,:),Smat(:,:,:)
complex(8),allocatable :: frequ_(:),green_lat(:,:,:,:),green_matsu(:,:,:,:)
logical                :: check,check2,checkfile
real(8)                :: mmu,vv,mmax
complex(8)             :: tt,ttt
complex(8),allocatable :: t1(:),t2(:),t3(:),dummy(:,:),dummy2(:,:),t_trace1(:,:),t_trace2(:,:)
real(8),allocatable    :: dens(:,:)
real(8)                :: beta,totdens

   call initialize_my_simulation
   testing=.false.
   fast_invmat=.true.
   use_cuda_routines=.false.
   use_cula_routines=.false.
   force_invmat_single_prec=.false.
   use_openmp_invmat=.false.
   diag_use_LU_instead_of_pivot=.false.
   flag_use_invmat_jordan=.true.
   flag_use_invmat_jordan_real=.true.
   enable_mpi_dot=.false.

   write(*,*) 'lattice green function'

   do kk_=1,2
     if(kk_==2)then
      call system(" ls green_onetep_spin_full_check_convergence_ren2 || cp green_onetep_spin_full_check_convergence_ren1 green_onetep_spin_full_check_convergence_ren2 ")
     endif
     open(unit=10012,file="green_onetep_spin_full_check_convergence_ren"//TRIM(ADJUSTL(toString(kk_))),form='unformatted')
     read(10012) n_frequ,channels
     write(*,*) 'Nfrequ channels : ', n_frequ,channels
     if(.not.allocated(Smat)) allocate(Smat(channels,channels,2),rotmat(channels,channels,2),frequ_(n_frequ),green_lat(n_frequ,channels,channels,2))
     read(10012) mmu,Smat(:,:,kk_),rotmat(:,:,kk_) 
     do i=1,n_frequ
       read(10012) frequ_(i),green_lat(i,:,:,kk_)
       green_lat(i,:,:,kk_)=rotit(green_lat(i,:,:,kk_),kk_)
     enddo
     close(10012)
   enddo

   write(*,*) 'read impurity green function'


   inquire(file='green_output_matsu_full_long_1',exist=check)
   !===========================================================!
   if(check)then

     open(unit=1000,file='green_output_matsu_full_long_1',form='unformatted')
     i=0
     do
      read(1000,end=111)
      i=i+1
     enddo
     111 continue
     close(1000)
     write(*,*) 'nmatsu_long is : ', i
     n_frequ_long=i
  
     allocate(green_matsu(n_frequ_long,channels,channels,2))
     do j=1,2
     if(j==2)then
      call system(" ls green_output_matsu_full_long_2 || cp green_output_matsu_full_long_1 green_output_matsu_full_long_2 ")
     endif
     open(unit=1000,file='green_output_matsu_full_long_'//TRIM(ADJUSTL(toString(j))),form='unformatted')
     do i=1,n_frequ_long
      write(*,*) 'reading frequ : ', i,n_frequ_long
      read(1000)  green_matsu(i,:,:,j)
      green_matsu(i,:,:,j)=rotit(green_matsu(i,:,:,j),j)
     enddo
     close(1000)
     enddo
     write(*,*) 'done'

   !===========================================================!

   else
   !===========================================================!

     inquire(file='green_output_matsu_full1',exist=check2)
     if(check2)then
         allocate(green_matsu(n_frequ,channels,channels,2))
         n_frequ_long=n_frequ
         do j=1,2
         if(j==2)then
          call system(" ls green_output_matsu_full2 || cp green_output_matsu_full1 green_output_matsu_full2 ")
         endif
         open(unit=1000,file='green_output_matsu_full'//TRIM(ADJUSTL(toString(j))),form='unformatted')
         do i=1,n_frequ
           read(1000)  green_matsu(i,:,:,j)  
           green_matsu(i,:,:,j)=rotit(green_matsu(i,:,:,j),j)
         enddo
         close(1000)
         enddo
         write(*,*) 'done'
     else

           INQUIRE(file='green_output',exist=checkfile)
           if(checkfile)then
             write(*,*) 'is it a paramagnet ? (=1 for true)'
             read(*,*) paramagnet
             fac=4; if(paramagnet==1) fac=fac/2; mmax=2; if(paramagnet==1) mmax=1
             if(.not.allocated(temp)) allocate(temp(fac*channels),vvv(channels,2))
             open(unit=1000,file='green_output')
             read(1000,*)
             j=0
             do
              read(1000,*,end=213)
              j=j+1
             enddo
             213 continue; n_frequ_long=j; rewind(1000) ; read(1000,*)
             allocate(green_matsu(n_frequ_long,channels,channels,2))
             do i=1,n_frequ_long
              read(1000,*) vv, (temp(k),k=1,fac*channels)
              do j=1,mmax
               green_matsu(i,:,:,j)=0.d0
              do k=1,channels
               green_matsu(i,k,k,j)=temp( channels*(j-1)*2 + 2*k-1 ) + imi * temp( channels*(j-1)*2 + 2*k )
              enddo
               green_matsu(i,:,:,j)=rotit(green_matsu(i,:,:,j),j)
              enddo
               if(paramagnet==1) green_matsu(i,:,:,2)=green_matsu(i,:,:,1)
             enddo
             close(1000)

           else

             allocate(green_matsu(n_frequ,channels,channels,2))
             n_frequ_long=n_frequ
             write(*,*) 'FILE green_output_matsu_full not present (probably running single site DMFT version)'
             write(*,*) 'green_imp=green_lat, in order to carry on and compute density'
             green_matsu=green_lat

           endif

     endif
   endif
   !===========================================================!


   if(.not.allocated(dens))     allocate(dens(channels,channels))
   if(.not.allocated(dummy))    allocate(t1(channels),t2(channels),t3(channels),dummy(channels,channels),dummy2(channels,channels))
   if(.not.allocated(t_trace1)) allocate(t_trace1(channels,channels),t_trace2(channels,channels))

     beta=dacos(-1.d0)/aimag(frequ_(1))
     write(*,*) 'Chem           :  ', mmu
     write(*,*) 'beta           :  ', beta

   do kk_=1,2 
    do k=1,channels
     open(1010+100*k+1000*kk_,file='compare_imp_'//TRIM(ADJUSTL(toString(k)))//"_"//TRIM(ADJUSTL(toString(kk_))))
     open(1011+100*k+1000*kk_,file='compare_lat_'//TRIM(ADJUSTL(toString(k)))//"_"//TRIM(ADJUSTL(toString(kk_))))
     open(1012+100*k+1000*kk_,file='diff_gm1_'//TRIM(ADJUSTL(toString(k)))//"_"//TRIM(ADJUSTL(toString(kk_))))
    enddo
   enddo

   open(unit=100,file='rho_fermi_level_from_gimp')
   open(unit=101,file='rho_fermi_level_from_gproj')
   open(unit=102,file='rho_n_tot_from_gimp')
   open(unit=103,file='rho_matrix_from_gimp')

   totdens=0.d0

   !===============================================================!
   !===============================================================!
   !===============================================================!
   do kk_=1,2

    do i=1,n_frequ
        t1=diag(green_matsu(i,:,:,kk_))
        t2=diag(green_lat(i,:,:,kk_))
        t3=t2-t1
        t_trace1=matmul(green_matsu(i,:,:,kk_),Smat(:,:,kk_))
        t_trace2=matmul(  green_lat(i,:,:,kk_),Smat(:,:,kk_))
        do k=1,channels 
           write(1010+100*k+1000*kk_,*) aimag(frequ_(i)),real(t1(k)),aimag(t1(k))
           write(1011+100*k+1000*kk_,*) aimag(frequ_(i)),real(t2(k)),aimag(t2(k))
           write(1012+100*k+1000*kk_,*) aimag(frequ_(i)),real(t3(k)),aimag(t3(k))
           if(i==1.or.i==2)then
              write(100,*)   -aimag(t_trace1(k,k) ) / dacos(-1.d0)
              write(101,*)   -aimag(t_trace2(k,k) ) / dacos(-1.d0)
           endif
        enddo
    enddo

     do i1=1,channels
      do i2=1,channels
        if(i1==i2) then
          dens(i1,i2)=get_dens(real(green_matsu(:,i1,i2,kk_)),aimag(green_matsu(:,i1,i2,kk_)),diag=.true.)
        else
          dens(i1,i2)=get_dens(real(green_matsu(:,i1,i2,kk_)),aimag(green_matsu(:,i1,i2,kk_)),diag=.false.)
        endif
      enddo
     enddo 

     dens=matmul(dens,Smat(:,:,kk_))

     do i1=1,channels
       write(102,*)  dens(i1,i1)
     enddo

     write(*,*) 'TOTAL DENSITY spin [x] : ', kk_, sum(diag(dens))
     totdens=totdens+sum(diag(dens))

     do i1=1,channels
      do i2=1,channels
       write(103,'(a,2i4,f16.10)') 'i,j=',i1,i2,dens(i1,i2)
      enddo
     enddo

   enddo
   !===============================================================!
   !===============================================================!
   !===============================================================!

   close(100)
   close(101)
   close(102)
   close(103)

   write(*,*) 'TOTAL DENSITY (UP+DN) : ', totdens

   do kk_=1,2
    do k=1,channels
     close(1010+100*k+1000*kk_)
     close(1011+100*k+1000*kk_)
     close(1012+100*k+1000*kk_)
    enddo
   enddo

  call finalize_my_simulation
  write(*,*) 'CHECK_CONVERGENCE_FINISHED'

contains

!------------------------!
!------------------------!
!------------------------!
!------------------------!
!------------------------!

function rotit(mat,j)
implicit none
complex(8) :: mat(:,:),rotit(size(mat,1),size(mat,2))
integer    :: j
  rotit=MATMUL(MATMUL(transpose(rotmat(:,:,j)),mat),rotmat(:,:,j))
end function

!------------------------!
!------------------------!
!------------------------!
!------------------------!
!------------------------!

real(8) function get_dens(GlocRe,GlocIm,diag)
implicit none
integer    :: k1,k2,i,j,k,jj
complex(8) :: dens
real(8)    :: frequ,pi,tau0,omega,df,GlocRe(:),GlocIm(:)
complex(8) :: temp
logical    :: diag
real(8)    :: alpha

   pi=dacos(-1.d0)
   jj=size(GlocIm)-7
   frequ=pi/beta*(2.d0*dble(jj)-1)
   alpha=-GlocIm(jj)*frequ
   write(*,*) 'GET DENS ALPHA COEF : ', alpha

   tau0 = 1.d-8; dens = 0.d0
   do j=1,n_frequ_long
     omega   = pi/beta*(2.d0*dble(j)-1.d0)
     if(abs(omega)<1.d-9) omega=1.d-9
     temp    = CMPLX(GlocRe(j),GlocIm(j),kind=8)
     if(diag)then
      dens = dens + 2.d0/beta * ( MPLX( omega * tau0 ) * (temp  + imi / omega *alpha )  )
     else
      dens = dens + 2.d0/beta * ( MPLX( omega * tau0 ) * (temp)  )
     endif
   enddo
   if(diag)then
    dens    = dens + 0.5 * alpha
   endif
   get_dens = dens

return
end function

!------------------------!
!------------------------!
!------------------------!
!------------------------!

end program

!======================!
!======================!



