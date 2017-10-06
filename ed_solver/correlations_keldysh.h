   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!

  subroutine keldysh
  implicit none
  integer                                :: lesser_bigger,ispin,ii,jj,i,j,Nc,isector,ieigen,isector2,ieigen2
  integer                                :: pmi,up1,do1,up2,do2,sz1,sz2
  real(8)                                :: time_mesh(keldysh_n),n_versus_t(keldysh_n,2)
  real(8),allocatable                    :: n_v_t_matrix(:,:,:,:)
  TYPE(eigensectorlist_type),allocatable :: phi(:),psi(:,:,:),tmp(:)
  TYPE(eigensectorlist_type)             :: tmp0
  LOGICAL,allocatable                    :: ORBMASKvec(:)
  TYPE(eigensector_type)                 :: Apm_es
  TYPE(sector_type)                      :: Asec
  complex(8),allocatable                 :: green_lesser(:,:,:,:,:),green_bigger(:,:,:,:,:),green_ret(:,:,:,:,:),green_ad(:,:,:,:,:)
  REAL(DBL)                              :: Zpart1,boltz,E01,Zpart2,E02,e1,e2
  COMPLEX(DBL)                           :: cscal
  INTEGER                                :: mmmax,ntot

  if(.not.do_keldysh) return

  write(*,*) 'STARTING KELDYSH AT FINITE TEMPERATURE WITH [x] ORBITALS : ', G(1)%N 
  Nc=G(1)%N

  if(allocated(ORBMASKvec)) deallocate(ORBMASKvec)
  allocate(ORBMASKvec(Nc))

  call build1Dmesh(time_mesh,size(time_mesh),keldysh_t0,keldysh_tmax)
  keldysh_delta=(keldysh_tmax-keldysh_t0)/dble(keldysh_n+1)

  write(*,*) ' ==== compute keldysh ==== '
  if(allocated(phi)) deallocate(phi)
  allocate(phi(keldysh_n))
  if(allocated(psi)) deallocate(psi,tmp)
  allocate(psi(keldysh_n,Nc,2),tmp(keldysh_n))

  if(allocated(green_lesser)) deallocate(green_lesser,green_bigger,green_ret,green_ad)
  allocate(green_lesser(keldysh_n,keldysh_n,Nc,Nc,2))
  allocate(green_bigger(keldysh_n,keldysh_n,Nc,Nc,2))
  allocate(green_ret(keldysh_n,keldysh_n,Nc,Nc,2))
  allocate(green_ad(keldysh_n,keldysh_n,Nc,Nc,2))

 !level 0 - evolution |phi(t)>=exp(-iHt)|0>

  call copy_eigensectorlist(phi(1),GS)
  do_quench=do_quench_

  if(keldysh_pert_ground_sector) then
     call solve_AIM(GS,AIM)
     call copy_eigensectorlist(phi(1),GS)
  else
     call solve_AIM(GS,AIM,skip_filter=.true.)
     do isector2=1,phi(1)%nsector
      do ieigen2=1,phi(1)%es(isector2)%lowest%neigen
         do isector=1,GS%nsector
          do ieigen=1,GS%es(isector)%lowest%neigen
            if(.not.AIM%bath%SUPER)then
              up1=    GS%es(isector )%sector%updo%up%npart
              do1=    GS%es(isector )%sector%updo%down%npart
              up2=phi(1)%es(isector2)%sector%updo%up%npart
              do2=phi(1)%es(isector2)%sector%updo%down%npart
              if(up1==up2.and.do1==do2)then
                goto 33
              endif
            else
              sz1 =     GS%es(isector )%sector%sz%npart
              sz2 = phi(1)%es(isector2)%sector%sz%npart
              if(sz1==sz2)then
                goto 33
              endif
            endif
          enddo
         enddo
         goto 989
         33 continue
         phi(1)%es(isector2)%lowest%eigen(ieigen2)%vec%rc=GS%es(isector)%lowest%eigen(ieigen)%vec%rc
         write(*,*) 'NEW EIGENVECTOR FOR SECTOR/EIGENVALUE/VAL : ', isector2,ieigen2,phi(1)%es(isector2)%lowest%eigen(ieigen2)%val
         989 continue
      enddo
     enddo
  endif


  do_quench=0

  do i=1,keldysh_n-1
   write(*,*) 'KELDYSH ITERATION', i
   call copy_eigensectorlist(phi(i+1),phi(i)) 
   CALL compute_greenAA(GKret(1),AIM,beta_ED,phi(i),N_sector,apply_Id,COMPUTE_DYN=.true.,keldysh_level=1,GS_out=phi(i+1))
  enddo

 !level 1 - apply c_{spin,orbital}|phi(t)>

  if(do_keldysh_gbigger)then
   mmmax=2
  else
   mmmax=1
  endif

  green_lesser=0.d0
  green_bigger=0.d0

  do lesser_bigger=1,mmmax

  if(lesser_bigger==1) pmi=2
  if(lesser_bigger==2) pmi=1

  do ispin=1,2
  do jj=1,Nc
  do ii=1,keldysh_n

  if(jj==quench_orb.or.quench_orb==-1)then

       call copy_eigensectorlist(tmp0,phi(ii))
       ORBMASKvec=F
       ORBMASKvec(jj)=T
       call delete_eigensectorlist(psi(ii,jj,ispin))
       DO isector=1,tmp0%nsector
         if(ispin==1)then
          CALL Cup_sector(Asec,pm(pmi),tmp0%es(isector)%sector)
         else
          CALL Cdo_sector(Asec,pm(pmi),tmp0%es(isector)%sector)
         endif
          CALL new_eigensector(Apm_es,Asec)
          CALL delete_eigenlist(Apm_es%lowest)
          DO ieigen=1,tmp0%es(isector)%lowest%neigen
           if(ispin==1)then
             CALL apply_Cup(Apm_es,pm(pmi),ORBMASKvec,tmp0%es(isector),ieigen)
           else
             CALL apply_Cdo(Apm_es,pm(pmi),ORBMASKvec,tmp0%es(isector),ieigen)
           endif
          ENDDO
          CALL add_eigensector(Apm_es,psi(ii,jj,ispin)) 
       ENDDO

  endif
 
  enddo
  enddo
  enddo

 !level 2 - evolve  c_{spin,orbital}|phi(t)> back to t=0

  keldysh_delta=-keldysh_delta

  do ispin=1,2
  do jj=1,Nc
  write(*,*) 'EVOLVING ORBITAL/SPIN : ',jj,ispin
  if(jj==quench_orb.or.quench_orb==-1)then

   do ii=1,keldysh_n
    call copy_eigensectorlist(tmp(1),psi(ii,jj,ispin))
    do i=1,ii-1
      write(*,*) 'KELDYSH ITERATION', i
      call copy_eigensectorlist(tmp(i+1),tmp(i))    
      CALL compute_greenAA(GKret(ispin),AIM,beta_ED,tmp(i),N_sector,apply_Id,COMPUTE_DYN=.true., &
                            & keldysh_level=jj,GS_out=tmp(i+1))
    enddo
    call copy_eigensectorlist(psi(ii,jj,ispin),tmp(ii))
   enddo

  endif
  enddo
  enddo

  do ispin=1,2  
  do i=1,keldysh_n
  do j=i,keldysh_n
  do ii=1,Nc
  do jj=1,Nc

  if((jj==quench_orb.and.ii==quench_orb).or.quench_orb==-1)then

   Zpart1 = partition(beta_ED,psi(i,ii,ispin))
   E01    = GSenergy(psi(i,ii,ispin))
   Zpart2 = partition(beta_ED,psi(j,jj,ispin))
   E02    = GSenergy(psi(j,jj,ispin))

   if(abs(Zpart1-Zpart2)>1.d-5)then
     write(*,*) 'ERROR partition functions not consistent !!'
     write(*,*) 'Zpart1,Zpart2 : ', Zpart1,Zpart2
     stop
   endif
   if(abs(E01-E02)>1.d-5)then
     write(*,*) 'ERROR partition functions not consistent !!'
     write(*,*) 'E01,E02 : ', E01,E02
     stop
   endif

   ntot=0
   do isector=1,psi(i,ii,ispin)%nsector
     ntot=ntot+psi(i,ii,ispin)%es(isector)%lowest%neigen
   enddo

   do isector=1,psi(i,ii,ispin)%nsector
    do ieigen=1,psi(i,ii,ispin)%es(isector)%lowest%neigen  

     do isector2=1,psi(j,jj,ispin)%nsector
      do ieigen2=1,psi(j,jj,ispin)%es(isector2)%lowest%neigen
       e1=psi(i,ii,ispin)%es(isector)%lowest%eigen(ieigen)%val
       e2=psi(j,jj,ispin)%es(isector2)%lowest%eigen(ieigen2)%val
       if(abs(e1-e2)<1.d-13)then
         if(.not.AIM%bath%SUPER)then
           up1=psi(i,ii,ispin)%es(isector)%sector%updo%up%npart
           do1=psi(i,ii,ispin)%es(isector)%sector%updo%down%npart
           up2=psi(j,jj,ispin)%es(isector2)%sector%updo%up%npart
           do2=psi(j,jj,ispin)%es(isector2)%sector%updo%down%npart
           if(up1==up2.and.do1==do2)then
             goto 31
           endif
         else
            sz1 = psi(i,ii,ispin)%es(isector)%sector%sz%npart
            sz2 = psi(j,jj,ispin)%es(isector2)%sector%sz%npart
            if(sz1==sz2)then
              goto 31 
            endif
         endif
       endif
      enddo 
     enddo

     goto 990
     31 continue
     boltz = DEXPc(-beta_ED*(psi(i,ii,ispin)%es(isector)%lowest%eigen(ieigen)%val-E01)) / Zpart1
     if(quench_cancel_statistics) boltz=1.d0 / dble(ntot)
     if(lesser_bigger==1)then
       cscal = scalprod( psi(j,jj,ispin)%es(isector2)%lowest%eigen(ieigen2)%vec%rc , psi(i,ii,ispin)%es(isector)%lowest%eigen(ieigen)%vec%rc )
       green_lesser(i,j,ii,jj,ispin)=green_lesser(i,j,ii,jj,ispin) + imi*boltz*cscal  
       if(i/=j) green_lesser(j,i,ii,jj,ispin) = -conjg(green_lesser(i,j,ii,jj,ispin))
     else
       cscal = scalprod( psi(i,ii,ispin)%es(isector)%lowest%eigen(ieigen)%vec%rc , psi(j,jj,ispin)%es(isector2)%lowest%eigen(ieigen2)%vec%rc )
       green_bigger(i,j,ii,jj,ispin)=green_bigger(i,j,ii,jj,ispin) - imi*boltz*cscal
       if(i/=j) green_bigger(j,i,ii,jj,ispin) = -conjg(green_bigger(i,j,ii,jj,ispin))
     endif
     990 continue

    enddo
   enddo 

  endif

  enddo
  enddo
  enddo
  enddo
  enddo

  enddo !lesser-bigger

  if(maxval(abs( green_lesser(:,:,1,1,1) - (-conjg(transpose(green_lesser(:,:,1,1,1)))) ))>1.d-5)then
   write(*,*) 'ERROR G^< not proper time symmetry'
   write(*,*) 'diff = : ', maxval(abs( green_lesser(:,:,1,1,1) - (-conjg(transpose(green_lesser(:,:,1,1,1)))) ))
  endif
  if(maxval(abs( green_bigger(:,:,1,1,1) - (-conjg(transpose(green_bigger(:,:,1,1,1)))) ))>1.d-5)then
   write(*,*) 'ERROR G^> not proper time symmetry'
   write(*,*) 'diff = ', maxval(abs( green_bigger(:,:,1,1,1) - (-conjg(transpose(green_bigger(:,:,1,1,1)))) ))
  endif

  do ispin=1,2
  do ii=1,Nc
  do jj=1,Nc
  green_ret(:,:,ii,jj,ispin)= +(green_bigger(:,:,ii,jj,ispin)-green_lesser(:,:,ii,jj,ispin))
  green_ad(:,:,ii,jj,ispin) = -(green_bigger(:,:,ii,jj,ispin)-green_lesser(:,:,ii,jj,ispin))
  do i=1,keldysh_n
  do j=1,keldysh_n
   if(i<j) green_ret(i,j,ii,jj,ispin)=0.d0
   if(i>j)  green_ad(i,j,ii,jj,ispin)=0.d0 
  enddo
  enddo 
  enddo
  enddo
  enddo


  if(allocated(n_v_t_matrix)) deallocate(n_v_t_matrix)
  allocate(n_v_t_matrix(keldysh_n,Nc,Nc,2))

  do i=1,keldysh_n
   do ii=1,Nc
    do jj=1,Nc
     do ispin=1,2
      n_v_t_matrix(i,ii,jj,ispin) = -imi*green_lesser(i,i,ii,jj,ispin) 
     enddo
    enddo
   enddo
  enddo

  do ispin=1,2 
   open(unit=918181,file='n_time_matrix_'//trim(adjustl(toString(ispin))),form='unformatted')
   write(918181) keldysh_n,Nc
   do i=1,keldysh_n
     write(918181) n_v_t_matrix(i,:,:,ispin) 
   enddo    
   close(918181)
  enddo

  do jj=1,Nc
  do ispin=1,2 
  call plotarray(time_mesh,real(green_lesser(1,:,jj,jj,ispin)),'Real_glesser_0T_orb_'//trim(adjustl(toString(jj)))//'_'//trim(adjustl(toString(ispin))))
  call plotarray(time_mesh,aimag(green_lesser(1,:,jj,jj,ispin)),'Im_glesser_0T_orb_'//trim(adjustl(toString(jj)))//'_'//trim(adjustl(toString(ispin))))
  call plotarray(time_mesh,real(green_bigger(:,1,jj,jj,ispin)),'Real_gbigger_T0_orb_'//trim(adjustl(toString(jj)))//'_'//trim(adjustl(toString(ispin))))
  call plotarray(time_mesh,aimag(green_bigger(:,1,jj,jj,ispin)),'Im_gbigger_T0_orb_'//trim(adjustl(toString(jj)))//'_'//trim(adjustl(toString(ispin))))
  call plotarray(time_mesh,real(green_lesser(:,1,jj,jj,ispin)),'Real_glesser_T0_orb_'//trim(adjustl(toString(jj)))//'_'//trim(adjustl(toString(ispin))))
  call plotarray(time_mesh,aimag(green_lesser(:,1,jj,jj,ispin)),'Im_glesser_T0_orb_'//trim(adjustl(toString(jj)))//'_'//trim(adjustl(toString(ispin))))
  call plotarray(time_mesh,real(green_bigger(1,:,jj,jj,ispin)),'Real_gbigger_0T_orb_'//trim(adjustl(toString(jj)))//'_'//trim(adjustl(toString(ispin))))
  call plotarray(time_mesh,aimag(green_bigger(1,:,jj,jj,ispin)),'Im_gbigger_0T_orb_'//trim(adjustl(toString(jj)))//'_'//trim(adjustl(toString(ispin))))
  call plotarray(time_mesh, real(green_ret(:,1,jj,jj,ispin)),'Real_gret_T0_orb_'//trim(adjustl(toString(jj)))//'_'//trim(adjustl(toString(ispin))))
  call plotarray(time_mesh,aimag(green_ret(:,1,jj,jj,ispin)),  'Im_gret_T0_orb_'//trim(adjustl(toString(jj)))//'_'//trim(adjustl(toString(ispin))))
  call plotarray(time_mesh, real(green_ret(1,:,jj,jj,ispin)),'Real_gret_0T_orb_'//trim(adjustl(toString(jj)))//'_'//trim(adjustl(toString(ispin))))
  call plotarray(time_mesh,aimag(green_ret(1,:,jj,jj,ispin)),  'Im_gret_0T_orb_'//trim(adjustl(toString(jj)))//'_'//trim(adjustl(toString(ispin))))
  do i=1,keldysh_n
   n_versus_t(i,ispin)=-imi*green_lesser(i,i,jj,jj,ispin)
  enddo
  call plotarray(time_mesh,n_versus_t(:,ispin),'n_versus_t_orb_'//trim(adjustl(toString(jj)))//'_'//trim(adjustl(toString(ispin))))
  enddo
  call plotarray(time_mesh,(n_versus_t(:,1)-n_versus_t(:,2))/2.d0,'mag_versus_t_orb_'//trim(adjustl(toString(jj))))
  enddo

  open(unit=918181,file='Glesser',form='unformatted')
  write(918181) shape(green_lesser)
  write(918181) green_lesser
  close(918181)

  open(unit=918181,file='Gbigger',form='unformatted')
  write(918181) shape(green_bigger)
  write(918181) green_bigger
  close(918181)

  open(unit=918181,file='Gret',form='unformatted')
  write(918181) shape(green_ret)
  write(918181) green_ret
  close(918181)

  open(unit=918181,file='Gad',form='unformatted')
  write(918181) shape(green_ad)
  write(918181) green_ad
  close(918181)

  if(allocated(green_lesser)) deallocate(green_lesser)
  if(allocated(green_bigger)) deallocate(green_bigger)
  if(allocated(green_ret))    deallocate(green_ret)
  if(allocated(green_ad))     deallocate(green_ad)

  call delete_eigensectorlist(tmp0) 
  do i=1,keldysh_n
     call delete_eigensectorlist(tmp(i))
  enddo 
  do ispin=1,2
   do i=1,keldysh_n
    do ii=1,Nc
     call delete_eigensectorlist(psi(i,ii,ispin))
    enddo
   enddo
  enddo 

 return
 end subroutine


   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!
   !------------------------!
