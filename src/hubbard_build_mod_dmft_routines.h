!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine read_gbooth_partial_kernel
implicit none

 if(gbooth)then
   if(ien==1.and..not.silent) write(*,*) 'checking if file [x] is present : ', filename
   call comms_barrier
   INQUIRE(file='store_gbooth_partial_kernel1',EXIST=partial_dens_present)
   call comms_barrier
   if(.not.check)then
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'file is not present. Not calculating subsystem energies.'
!    call sparse_write_(mat,filename=filename)
!    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'done...'
   else
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'file is present, read from this file'
    call sparse_read_(denskern_partial(1),filename='store_gbooth_partial_kernel1')
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'done...'
   endif
   call comms_barrier
   if(num_spins==2) stop 'UHF not coded up'
   if(pub_dmft_paramagnetic .and. pub_cell%num_spins==2) then
       call sparse_copy(denskern_partial(2),denskern_partial(1))
   endif
!   !BUG CORRECTED - EDOARDO
!   call system("rm *.dkn_dmft")
!   call restart_kernel_write(denskern_tmp,write_dmft=.true.)
 endif
! call gbooth_kernel_read_
end subroutine


subroutine read_gbooth_write_kernel
implicit none
logical :: is_present
 if(gbooth)then
   call check_sparse_in_file(denskern_tmp(1),'store_gbooth_kernel1')
   if(num_spins==2)then
     INQUIRE(file='store_gbooth_kernel2',EXIST=is_present)
     if(is_present)then
       call check_sparse_in_file(denskern_tmp(2),'store_gbooth_kernel2')
     endif
   endif
   if(pub_dmft_paramagnetic .and. pub_cell%num_spins==2) then
       call sparse_copy(denskern_tmp(2),denskern_tmp(1))
   endif
   !BUG CORRECTED - EDOARDO
   call system("rm *.dkn_dmft")
   call restart_kernel_write(denskern_tmp,write_dmft=.true.)
 endif
 call gbooth_kernel_read_
end subroutine

   subroutine gbooth_kernel_read_
   use sparse,   only: sparse_convert,sparse_num_rows,sparse_num_cols
      implicit none
        real(kind=DP), allocatable, dimension(:,:) :: mat_square
        integer        :: n1,n2
          call comms_barrier
          n1=NINT(sparse_num_rows(denskern(1)))
          n2=NINT(sparse_num_cols(denskern(1)))
          allocate(mat_square(n1,n2),stat=ierr)
          call sparse_convert(mat_square,denskern(1))
          write(*,*) '# denskern printed from onetep'
          do i = 1,n1
              do j = 1,n1
                  write(*,'(F9.4)',advance='no') mat_square(i,j)
              enddo
              write(*,*)
          enddo
          deallocate(mat_square,stat=ierr)
       end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_parallel_scheme

    !faster way to get the energy when running mpi...
    split=(present(dmft_energy_cor).and.pub_total_num_nodes>1).or.(pub_dmft_split.and.pub_total_num_nodes>1)
    silent=present(dmft_energy_cor)
    if(split) then
       if(pub_my_node_id==0)then
        write(*,*) 'SPLIT IS NOW ACTIVE ...... '
        write(*,*) 'SPLITTING THE FREQUENCY SUMMATION OVER THE NODES'
       endif
    endif
                        breaknfscom=.false.
    if(pub_dmft_splitk) breaknfscom=.true.

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine finalize_kernel
implicit none

    if( pub_dmft_temp>0.0_DP .and. (.not. wont_compute_residues.or.mpi_rank==1) ) then

     if(pub_dmft_kernel==-1.and..not.present(dmft_energy_cor).and..not.present(dmft_kernel))then
       call kernel_purify(pur_denskern_tmp,denskern_tmp,overlap,inv_overlap,n_occ)
       do is=1,pub_cell%num_spins
          write(*,*) 'Trace PURE-DensKern(is) Overlap : ', sparse_trace(pur_denskern_tmp(is),overlap)
          write(*,*) 'Trace PURE-DensKern(is) HAM(is) : ', sparse_trace(pur_denskern_tmp(is),ham(is))
          call sparse_copy(denskern_tmp(is),pur_denskern_tmp(is))
       enddo
      endif

     if( (pub_dmft_kernel/=0.or.present(dmft_kernel)).and..not.present(dmft_energy_cor) )then
          if(mpi_rank==1.and.pub_my_node_id==0) write(*,*) '-!KERNEL MIXING!-'
          if(.not.present(dmft_kernel)) write(*,*) 'READING KERNEL FROM FILE : ',pub_rootname
          if(.not.present(dmft_kernel)) call restart_kernel_read(denskern)

          !BUG September 2014, to impose paramagnetic calculation
          if(pub_cell%num_spins==2.and.pub_dmft_paramagnetic)then
             call sparse_copy(denskern(2),denskern(1))
          endif
          !END BUG

          if(present(dmft_kernel).or.pub_dmft_kernel/=2)then
           if(mpi_rank==1.and.pub_my_node_id==0) write(*,*) 'LINEAR MIXING'
           do is=1,pub_cell%num_spins
              call sparse_scale(denskern_tmp(is),pub_dmft_kernel_mix)
              if(pub_cell%num_spins==1.and.present(dmft_kernel))then
                 !when DMFT is called in subroutine Hamiltonian_build_dens, the kernel has a factor x2
                 call sparse_axpy(denskern_tmp(is),denskern(is),0.5*(1.d0-pub_dmft_kernel_mix))
                 call sparse_scale(denskern_tmp(is),2.d0)
              else
                 call sparse_axpy(denskern_tmp(is),denskern(is),1.d0-pub_dmft_kernel_mix)
              endif
           enddo
          endif

!bug was v1
!#ifndef debugMIXv1
!          sizdiis=min(iterdiis,pub_kernel_diis_max)
!#else
           sizdiis=pub_kernel_diis_max
!#endif
          if(.not.present(dmft_kernel).and.pub_dmft_kernel==2)then
            pub_kernel_diis_c_in  = 1.0 - pub_dmft_kernel_mix
            pub_kernel_diis_c_out =       pub_dmft_kernel_mix
            if(mpi_rank==1.and.pub_my_node_id==0) write(*,*) 'PULLAY MIXING'
             ientry = 1 ; shifted  = .false. ; iterdiis = 1 ;
             call kernel_diis_sparse_create(dkn_in, dkn_out, next_dkn_in, residues, denskern)
             inquire(file='diis_kernel',exist=checkdiis)
             if(checkdiis)then
                open(file='diis_kernel',unit=301020,form='unformatted')
                read(301020) ientry,shifted,iterdiis
                close(301020)
                do is=1,pub_cell%num_spins
                 do i=1,sizdiis
                   call sparse_read_(dkn_in(is,i)  ,'diis_dkn_in_'//trim(adjustl(toString(i)))//'_'//trim(adjustl(toString(is))))  
                   call sparse_read_(dkn_out(is,i) ,'diis_dkn_out_'//trim(adjustl(toString(i)))//'_'//trim(adjustl(toString(is))))  
                   call sparse_read_(residues(is,i),'diis_res_'//trim(adjustl(toString(i)))//'_'//trim(adjustl(toString(is))))  
                 enddo
                enddo
             else
                  call kernel_diis_init(residues, next_dkn_in, dkn_out, dkn_in, denskern)
             endif
             if(mpi_rank==1.and.pub_my_node_id==0) then
              write(*,*) '-------------BEFORE ITER--------------------'
              write(*,*) 'IENTRY,SHIFTED,ITER', ientry,shifted,iterdiis
              write(*,*) '--------------------------------------------'
             endif
             do is=1,pub_cell%num_spins
               call sparse_copy(dkn_out(is,ientry),denskern_tmp(is))
             enddo
             if(mpi_rank==1.and.pub_my_node_id==0) write(*,*) 'residue inout'
             call kernel_diis_residue_inout(residues(:,ientry),dkn_out(:,ientry),dkn_in(:,ientry))
             if(mpi_rank==1.and.pub_my_node_id==0) write(*,*) 'diis mixt'
             call kernel_diis_mix(next_dkn_in,dkn_in,dkn_out,residues,overlap,iterdiis,ientry)
             if(mpi_rank==1.and.pub_my_node_id==0) write(*,*) 'diis shift'
             call kernel_diis_shift(shifted,dkn_in,dkn_out,residues,iterdiis)
             if(mpi_rank==1.and.pub_my_node_id==0) write(*,*) 'diis find ientry'
             call kernel_diis_find_ientry(ientry,iterdiis,shifted)
             iterdiis=iterdiis+1
             do is=1,pub_cell%num_spins
                call sparse_copy(dkn_in(is,ientry),next_dkn_in(is))
                call sparse_copy(denskern_tmp(is),next_dkn_in(is))
             end do

             if(mpi_rank==1.and.pub_my_node_id==0) then
                write(*,*) '-------------AFTER ITER---------------------'
                write(*,*) 'IENTRY,SHIFTED,ITER', ientry,shifted,iterdiis
                write(*,*) '--------------------------------------------'
                open(file='_diis_kernel',unit=301020,form='unformatted')
                write(301020) ientry,shifted,iterdiis
                close(301020)
             endif

             if(mpi_rank==1)then
              if(mpi_rank==1.and.pub_my_node_id==0) write(*,*) 'writing output files'
              do is=1, pub_cell%num_spins
               do i=1,sizdiis
                call sparse_write_(dkn_in(is,i),'_diis_dkn_in_'//trim(adjustl(toString(i)))//'_'//trim(adjustl(toString(is))))  
                call sparse_write_(dkn_out(is,i),'_diis_dkn_out_'//trim(adjustl(toString(i)))//'_'//trim(adjustl(toString(is))))  
                call sparse_write_(residues(is,i),'_diis_res_'//trim(adjustl(toString(i)))//'_'//trim(adjustl(toString(is))))  
               enddo
              enddo
             endif

             if(mpi_rank==1.and.pub_my_node_id==0) write(*,*) 'destroy files'
             call kernel_diis_sparse_destroy(dkn_in, dkn_out, next_dkn_in, residues)

          endif

          if(mpi_rank==1.and.pub_my_node_id==0) write(*,*) 'dumping output DMFT kernel'

          if(mpi_rank==1.and..not.present(dmft_kernel)) call restart_kernel_write(denskern_tmp,write_dmft=.true.)
          if(.not.present(dmft_kernel)) write(*,*) 'WARNING/ERROR : now ERASES DFT kernel, and replace it with DMFT KERNEL'
          if(.not.present(dmft_kernel))then
              do is=1,pub_cell%num_spins
                  call sparse_copy(denskern(is),denskern_tmp(is))
              enddo
          else
              do is=1,pub_cell%num_spins
                  call sparse_copy(dmft_kernel(is),denskern(is))
                  call sparse_copy(denskern(is),denskern_tmp(is))
              enddo
          endif

       endif

    endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrapper_comm_reduce_energy_and_kernel
       if(split)then
         write(*,*) 'SYNCING KERNEL'
         call comms_reduce('SUM',dmft_Energy(is))
        if(nkpoints==1)then
         if(.not.allocated(matrixDensKernel))then
           write(*,*) 'ERROR matrixDensKernel is not allocated after k-sum loop';stop
         endif
         call comms_reduce('SUM',matrixDensKernel)
        else
         if(.not.allocated(matrixDensKernelc))then
           write(*,*) 'ERROR matrixDensKernelc is not allocated after k-sum loop';stop
         endif
         call comms_reduce('SUM',matrixDensKernelc)
        endif
       endif
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrapper_add_self_and_invert_green
implicit none

          if(verboseall)write(*,*) 'compute green-func'
          if(.not.split)then
             if(nkpoints==1)then
               call sparse_axpy(inv_greenf,self_energy,cmplx(-1.0,0.0,kind=DP))
             else
               call sparse_axpy(inv_greenf,selfK      ,cmplx(-1.0,0.0,kind=DP))
             endif
             call invert_green_function
          else

           if(pub_dmft_lin_scaling)then
             if(ien/=pub_dmft_points-3.and.ien/=pub_dmft_points-4.and.ien/=pub_dmft_points-5.and.ien/=pub_dmft_points-6)then
                call internal_invert_linear_scaling(greenf,energy,overlap_cplx,ham(is),self_energy,update=.false.,realorcomplexH=.true.)
             else
                 if(ien==pub_dmft_points-4)then
                   greenfbackup=MATMUL(split_Sm,MATMUL(-split_H,split_Sm))
                 elseif(ien==pub_dmft_points-3.or.ien==pub_dmft_points-6)then
                   greenfbackup=split_Sm
                 elseif(ien==pub_dmft_points-5)then
                   greenfbackup=MATMUL(split_Sm,MATMUL(-sigmabackup,split_Sm))
                 endif
             endif
           else
             if(ien/=pub_dmft_points-3.and.ien/=pub_dmft_points-4.and.ien/=pub_dmft_points-5.and.ien/=pub_dmft_points-6)then
                 greenfbackup = greenfbackup - sigmabackup
                 if(nkpoints==1)then
                    call wrappers_invert_sym_cmatrix(greenfbackup,ngwf_basis%num)
                 else
                    call invert_gen_cmat(ngwf_basis%num,greenfbackup)
                 endif
             else
                 if(ien==pub_dmft_points-4)then
                   greenfbackup=MATMUL(split_Sm,MATMUL(-split_H,split_Sm))
                 elseif(ien==pub_dmft_points-3.or.ien==pub_dmft_points-6)then
                   greenfbackup=split_Sm
                 elseif(ien==pub_dmft_points-5)then
                   greenfbackup=MATMUL(split_Sm,MATMUL(-sigmabackup,split_Sm))
                 endif
             endif
           endif

          endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wrapper_build_self_energy
implicit none

          if(.not.split)then
             kien=ien
             call self_energy_header_write_screen
             call construct_self_energy
             call upfold_self_energy
          else
             do kkien = 1, kstep
               if(pub_my_node_id==0) write(*,*) 'BUILDING SELF ENERGY : ', kkien
               kien = iien - kstart + (kkien-1)
               call construct_self_energy
               call upfold_self_energy
               if(nkpoints==1)then
                call sparse_convert(tmp1,self_energy)
               else
                call sparse_convert(tmp1,selfK      )
               endif
               if(pub_my_node_id==kkien-1)then
                 sigmabackup = tmp1
               endif
             enddo
          endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine build_inv_green_matrix

          !=========================================!
          ! Create G^-1(E) = E S - H - (Self - Edc) !
          !=========================================!

          if(verboseall)write(*,*) 'Building GREEN'
          if(verboseall)write(*,*) 'copy [1]'

          if(.not.split) call sparse_copy(inv_greenf,overlap_cplx)

          if(ien/=pub_dmft_points-3.and.ien/=pub_dmft_points-6)then

             if(verboseall)write(*,*) 'copy [2]'
             if(.not.split)then
               call sparse_scale(inv_greenf,energy)
               if(verboseall)write(*,*) 'copy [3]'
               if(ien/=pub_dmft_points-5) then
                  if(nkpoints==1)then
                   call sparse_axpy(inv_greenf,ham(is),cmplx(-1.0,0.0,kind=DP))
                  else
                   call sparse_axpy(inv_greenf,hamK   ,cmplx(-1.0,0.0,kind=DP))
                  endif
                   if(check_embed_h.and.ien/=pub_dmft_points-4)then
                     call sparse_axpy(inv_greenf,ham_embed(is),cmplx(-1.d0,0.d0,kind=DP))
                   endif
               endif
             else
                if(ien/=pub_dmft_points-5.and.ien/=pub_dmft_points-4) greenfbackup=energy*split_S-split_H
                if(check_embed_h)then
                  write(*,*) 'ERROR embedding and split - not compatible - here we would add ham_embed(is) to greenfbackup' ; stop
                endif
             endif

          endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_kpoints
  implicit none

    jp=abs(pub_dmft_nkpoints)
    if(jp>1)then

      nkpoints=0

      if(pub_dmft_nkpoints<0)then
       do ii1=-(jp-1)/2,(jp-1)/2
        do ii2=-(jp-1)/2,(jp-1)/2
         do ii3=-(jp-1)/2,(jp-1)/2
          if(  ii3>0  .or.  (ii3==0  .and.  ( (ii2>0).or.(ii2==0.and.ii1>=0)  ) ) )then
           nkpoints=nkpoints+1
          endif
         enddo
        enddo
       enddo
       if(allocated(norm_kpoints)) deallocate(norm_kpoints,kpt)
       allocate(norm_kpoints(nkpoints),kpt(nkpoints,3))
       j=0
       do ii1=-(jp-1)/2,(jp-1)/2
        do ii2=-(jp-1)/2,(jp-1)/2
         do ii3=-(jp-1)/2,(jp-1)/2
          if(  ii3>0  .or.  (ii3==0  .and.  ( (ii2>0).or.(ii2==0.and.ii1>=0)  ) ) )then
           j=j+1
           t1_=dble(ii1) /  dble(jp)
           t2_=dble(ii2) /  dble(jp)
           t3_=dble(ii3) /  dble(jp)
           kpt(j,1)=t1_*pub_cell%b1%x+t2_*pub_cell%b2%x+t3_*pub_cell%b3%x
           kpt(j,2)=t1_*pub_cell%b1%y+t2_*pub_cell%b2%y+t3_*pub_cell%b3%y
           kpt(j,3)=t1_*pub_cell%b1%z+t2_*pub_cell%b2%z+t3_*pub_cell%b3%z
          endif
         enddo
        enddo
       enddo
       do i=1,nkpoints
         norm=2.d0
         if(sum(abs(kpt(i,:)**2))<1.d-5) norm=norm/2.0 !Gamma Point
         norm_kpoints(i)=norm
       enddo
      else
       nkpoints=jp**3
       if(allocated(norm_kpoints)) deallocate(norm_kpoints,kpt)
       allocate(norm_kpoints(nkpoints),kpt(nkpoints,3))
       j=0
       do ii1=1,jp
        do ii2=1,jp
         do ii3=1,jp
          j=j+1
          t1_=dble(ii1-1)/dble(2*jp-1)
          t2_=dble(ii2-1)/dble(2*jp-1)
          t3_=dble(ii3-1)/dble(2*jp-1)
          kpt(j,1)=t1_*pub_cell%b1%x+t2_*pub_cell%b2%x+t3_*pub_cell%b3%x
          kpt(j,2)=t1_*pub_cell%b1%y+t2_*pub_cell%b2%y+t3_*pub_cell%b3%y
          kpt(j,3)=t1_*pub_cell%b1%z+t2_*pub_cell%b2%z+t3_*pub_cell%b3%z
         enddo
        enddo
       enddo
       do i=1,nkpoints
         norm=8.0
         do j=1,3
            if(abs(kpt(i,j))<1.d-5) norm=norm/2.0 !we only keep mirrors...
         enddo
         norm_kpoints(i)=norm
       enddo
      endif

      lastkpoint=nkpoints

      if(pub_dmft_kpoints_sym.and.(pub_dmft_kernel==0.or.pub_dmft_kpoints_kernel_gamma))then
      !WARNING: CUBIC SYM ON KPOINTS NOT COMPATIBLE WITH KERNEL
       do i=1,nkpoints
        do j=1,nkpoints
         if(j/=i)then
          do ii1=-1,1,2
          do ii2=-1,1,2
          do ii3=-1,1,2
            kpt(j,1)=dble(ii1)*kpt(j,1) ! mirror
            kpt(j,2)=dble(ii2)*kpt(j,2) ! mirror
            kpt(j,3)=dble(ii3)*kpt(j,3) ! mirror
            if(norm_vector(kpt(i,:)-kpt(j,:))<1.d-5)then
              kpt(j,:)=1.d15
              norm_kpoints(i)=norm_kpoints(i)+norm_kpoints(j)
              norm_kpoints(j)=0.d0
            endif
            kpt(j,1)=dble(ii1)*kpt(j,1) ! mirror
            kpt(j,2)=dble(ii2)*kpt(j,2) ! mirror
            kpt(j,3)=dble(ii3)*kpt(j,3) ! mirror
          enddo
          enddo
          enddo
         endif
        enddo
       enddo
       lastkpoint=0
       do i=1,nkpoints
        if(abs(norm_kpoints(i))>1.d-3) lastkpoint=i
       enddo
      endif


      totkpoint=sum(norm_kpoints)
      norm_kpoints=norm_kpoints/totkpoint
      if(pub_my_node_id==0) then
         write(*,*) 'TOTAL NUMBER OF K POINTS                         : ', totkpoint
         write(*,*) 'TOTAL COMPUTED K POINTS (orthorombic symmetries) : ', nkpoints
         do i=1,nkpoints
          if(abs(norm_kpoints(i))>1.d-4)then
           write(*,*) 'KP - TOT             : ', i,nkpoints
           write(*,*) 'WG - TOT             : ', norm_kpoints(i),totkpoint
           write(*,'(a,3f10.5)') 'KP        : ', kpt(i,:)
          endif
         enddo
      endif
    else
       nkpoints=1;lastkpoint=1
       if(allocated(norm_kpoints)) deallocate(norm_kpoints,kpt)
       allocate(norm_kpoints(nkpoints),kpt(nkpoints,3))
       norm_kpoints=1.d0
       kpt=0.d0
    endif

    if(nkpoints==1)then
      if(pub_dmft_splitk) then
        write(*,*) 'ERROR : splitk because ONLY ONE K POINT' 
        stop
      endif
    endif

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine embedding_setup_planes
   implicit none
   integer :: nn,ufii,i,j,k,ii,jj,kk,iii,jjj,lll
   inquire(file='group_of_atoms_label',exist=planes_check)
   nplanes=0
   if(.not.planes_check) return
   ufii=free_unit('group_of_atoms_label')
   nn=lines_in_file(ufii)-1
   read(ufii,*) nplanes

   if(allocated(orb_in_plane)) deallocate(orb_in_plane); allocate(orb_in_plane(nplanes,ngwf_basis%num))
   if(allocated(atom_in_plane)) deallocate(atom_in_plane); allocate(atom_in_plane(nplanes,pub_cell%nat))

   orb_in_plane=.false.; atom_in_plane=.false.

   !--------------------!
   do i=1,nn
    read(ufii,*) k,j
    atom_in_plane(j,k)=.true. ! atom_in_plane is in ORIGINAL atom notations - for popn !
    jj=pub_distr_atom(k)
    jjj=0
    do ii=1,jj-1
     jjj=jjj+ngwf_basis%num_on_atom(ii)
    enddo
    do ii=1,ngwf_basis%num_on_atom(jj)
     orb_in_plane(j,jjj+ii)=.true.
    enddo
   enddo
   close(ufii)

   if(mpi_rank==1.and.pub_my_node_id==0)then
    open(unit=ufii,file='group_of_atoms_label_sorted',form='unformatted')
    write(ufii) nplanes
    write(ufii) size(atom_in_plane,1),size(atom_in_plane,2)
    write(ufii) atom_in_plane  
    close(ufii)     
   endif
   !--------------------!

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine embedding_prepare_tables
   implicit none
   integer             :: i,j,k,l,ufi,ufi1,kk,s1,s2,t1,t2,kkk,i1,j1,i0,j0,ii,ii0,jj0,ii1,jj1
   integer             :: tot_li,isp
   logical             :: check_embedl,check_embedr
   real(8),allocatable :: tmpH(:,:),ovH(:,:)
   real(8)             :: h_element

   if(pub_my_node_id==0.and.mpi_rank==1) write(*,*) 'EMBEDDING TABLE'

   check_embed=.false.
   inquire(file='connections_D_L',exist=check_embedl)
   inquire(file='connections_D_R',exist=check_embedr)
   check_embed=check_embedl.or.check_embedr
   if(.not.check_embed) return

   if(split)then
     write(*,*) 'ERROR - split NOT compatible with embedding - the reason is that the embedding'
     write(*,*) '        requires all atoms to be on the same node'
     write(*,*) '        otherwise you will need to code a sync routine to access the elements'
     write(*,*) '        otherwise you will need to code a sync routine to access the elements'
     STOP
   endif  
   if(nkpoints>1) then
     write(*,*) 'ERROR embedding AND k-points not yet possible at the moment'
     STOP
   endif
 
   if(check_embedl) ufi =free_unit('connections_D_L')
   if(check_embedr) ufi =free_unit('connections_D_R')
   if(check_embedl) ufi1=free_unit('connections_D_L_lead_sites')
   if(check_embedr) ufi1=free_unit('connections_D_R_lead_sites')

   if(check_embedl.and.check_embedr)then
     write(*,*) 'EMBEDDING ERROR, NOT SUPPOSED TO HAVE BOTH RIGHT AND LEFT FILES'
     stop
   endif

   if(allocated(embed_index)) deallocate(embed_index,embed_size,embed_pos,embed_pos_all,&
                               & maskH0_embed,maskH1_embed,maskH0_embed_atoms,maskH1_embed_atoms)

   kk=lines_in_file(ufi1)
   allocate(embed_index(kk),embed_size(kk),embed_pos(kk),embed_pos_all(pub_cell%nat))
   do i=1,kk
     read(ufi1,*) s2
     embed_index(i) = pub_distr_atom(s2)
   enddo

   kkk=0
   do i=1,kk
      embed_size(i) = ngwf_basis%num_on_atom( embed_index(i) )
      kkk           = kkk + embed_size(i)
      ! s2 in original positions (as in dat file) notations
      ! reminder: pub_orig_atom(iat) gives the original position, it iat is the NGWFs index
      ! reminder: pub_distr_atom is the opposite
   enddo

   do i=1,kk
    embed_pos(i)=0
    do j=1,embed_index(i)-1
       embed_pos(i)=embed_pos(i)+ngwf_basis%num_on_atom( j )
    enddo
   enddo
   do ii=1,pub_cell%nat
    embed_pos_all(ii)=0
    i=pub_distr_atom(ii)
    do j=1,i-1
      embed_pos_all(ii) = embed_pos_all(ii) + ngwf_basis%num_on_atom( j )
    enddo
   enddo

   allocate(maskH0_embed(kkk,kkk,2),maskH1_embed(kkk,kkk,2))
   allocate(maskH0_embed_atoms(kkk,kkk,2),maskH1_embed_atoms(kkk,kkk,2))

   do i=1,kk
    do j=1,kk
     i0=0 ; if(i>1) i0=sum(embed_size(1:i-1))
     j0=0 ; if(j>1) j0=sum(embed_size(1:j-1))
     do i1=1,embed_size(i)
      do j1=1,embed_size(j)
        maskH0_embed( i0+i1,j0+j1,1 )      = embed_pos(i) + i1
        maskH0_embed( i0+i1,j0+j1,2 )      = embed_pos(j) + j1
        maskH0_embed_atoms(i0+i1,j0+j1,1)  = embed_index(i)
        maskH0_embed_atoms(i0+i1,j0+j1,2)  = embed_index(j)
      enddo
     enddo
    enddo
   enddo

   tot_li=lines_in_file(ufi)
   do ii=1,tot_li

     read(ufi,*) s1,s2,t1,t2

     if(s1==0.or.s2==0.or.t1==0.or.t2==0)then
       write(*,*) 'BUILDING EMBEDDING TABLE / ERROR'
       write(*,*) 'INDICES ARE NAUGHT'
       write(*,*) 's1,s2,t1,t2 : ', s1,s2,t1,t2
       stop
     endif

     ii0=0;if(s1>1) ii0=sum(embed_size(1:s1-1))
     jj0=0;if(s2>1) jj0=sum(embed_size(1:s2-1))

     if(embed_size(s1)/=ngwf_basis%num_on_atom( pub_distr_atom(t1)))then
       write(*,*) 'BUILD EMBEDDING TABLE ERROR : size mismatch1'
       write(*,*) 's1,s2,t1,t2                                  : ', s1,s2,t1,t2
       write(*,*) 'S2 atom index                                : ', pub_orig_atom(embed_index(s2))
       write(*,*) 'embed_size(s1)                               : ', embed_size(s1)
       write(*,*) 'ngwf_basis%num_on_atom( pub_distr_atom(t1) ) : ', ngwf_basis%num_on_atom( pub_distr_atom(t1) ) 
       stop
     endif
     if(embed_size(s2)/=ngwf_basis%num_on_atom( pub_distr_atom(t2)))then
       write(*,*) 'BUILD EMBEDDING TABLE ERROR : size mismatch2'
       write(*,*) 's1,s2,t1,t2                                  : ', s1,s2,t1,t2
       write(*,*) 'S2 atom index                                : ', pub_orig_atom(embed_index(s2))
       write(*,*) 'embed_size(s2)                               : ', embed_size(s2)
       write(*,*) 'NUM_NG(S2 atom index)                        : ', ngwf_basis%num_on_atom( pub_distr_atom( pub_orig_atom(embed_index(s2))) )
       write(*,*) 'ngwf_basis%num_on_atom( pub_distr_atom(t2) ) : ', ngwf_basis%num_on_atom( pub_distr_atom(t2) ) 
       stop
     endif

     do i1=1,ngwf_basis%num_on_atom( pub_distr_atom(t1) )
      do j1=1,ngwf_basis%num_on_atom( pub_distr_atom(t2)  )
        maskH1_embed(ii0+i1,jj0+j1,1 )       =  embed_pos_all(t1) + i1
        maskH1_embed(ii0+i1,jj0+j1,2 )       =  embed_pos_all(t2) + j1
        maskH1_embed_atoms(ii0+i1,jj0+j1,1 ) = pub_distr_atom(t1)
        maskH1_embed_atoms(ii0+i1,jj0+j1,2 ) = pub_distr_atom(t2)
      enddo
     enddo
   enddo

   close(ufi)
   close(ufi1)

   allocate(tmpH(kkk,kkk),ovH(kkk,kkk))

   do isp=1,pub_cell%num_spins
   if(mpi_rank==1.and.pub_my_node_id==0) open(unit=ufi,file='_H0_'//trim(adjustl(toString(isp))),form='unformatted')
   tmpH=0.;ovH=0.
   do i=1,kkk
   do j=1,kkk
     i1=maskH0_embed(i,j,1) 
     j1=maskH0_embed(i,j,2)
     ii1=maskH0_embed_atoms(i,j,1)
     jj1=maskH0_embed_atoms(i,j,2)
     if(pub_node_of_atom(pub_orig_atom(ii1))==pub_my_node_id.and. &
      & pub_node_of_atom(pub_orig_atom(jj1))==pub_my_node_id) &
      &  call sparse_get_element(h_element,ham(isp),i1,j1)
     tmpH(i,j) = h_element
     if(pub_node_of_atom(pub_orig_atom(ii1))==pub_my_node_id.and. &
      & pub_node_of_atom(pub_orig_atom(jj1))==pub_my_node_id) & 
      & call sparse_get_element(h_element,overlap,i1,j1)
     ovH(i,j) = h_element
   enddo
   enddo
   call comms_reduce('SUM', tmpH)
   call comms_reduce('SUM', ovH)

   if(mpi_rank==1.and.pub_my_node_id==0) then
     write(ufi) kkk
     write(ufi) tmpH
     write(ufi) ovH
   endif
   if(mpi_rank==1.and.pub_my_node_id==0) close(ufi)
   if(mpi_rank==1.and.pub_my_node_id==0) open(unit=ufi,file='_H1_'//trim(adjustl(toString(isp))),form='unformatted')

   tmpH=0;ovH=0
   do i=1,kkk
   do j=1,kkk
     i1=maskH1_embed(i,j,1)
     j1=maskH1_embed(i,j,2)
     ii1=maskH1_embed_atoms(i,j,1)
     jj1=maskH1_embed_atoms(i,j,2)
    if(pub_node_of_atom(pub_orig_atom(ii1))==pub_my_node_id.and. &
     & pub_node_of_atom(pub_orig_atom(jj1))==pub_my_node_id) & 
     &call sparse_get_element(h_element,ham(isp),i1,j1)
     tmpH(i,j) = h_element
     if(pub_node_of_atom(pub_orig_atom(ii1))==pub_my_node_id.and. &
      & pub_node_of_atom(pub_orig_atom(jj1))==pub_my_node_id) & 
      & call sparse_get_element(h_element,overlap,i1,j1)
     ovH(i,j) = h_element
   enddo
   enddo
   call comms_reduce('SUM', tmpH)
   call comms_reduce('SUM', ovH)

   if(mpi_rank==1.and.pub_my_node_id==0) then
     write(ufi) kkk
     write(ufi) tmpH
     write(ufi) ovH
   endif
   if(mpi_rank==1.and.pub_my_node_id==0) close(ufi)
   enddo

   deallocate(tmpH,ovH)
 
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine embedding_write_sigma
   implicit none
   complex(8),allocatable :: ttt0(:,:)
   integer                :: size_mat0,size_mat1,i,j,s,t,ii,jj
   complex(8)             :: sig_element

   if(.not.check_embed) return

   if(pub_my_node_id==0.and.mpi_rank==1) write(*,*) 'WRITING EMBEDDING POTENTIALS'

   if(allocated(ttt0)) deallocate(ttt0)
   allocate(ttt0(size(maskH0_embed,1),size(maskH0_embed,1)))
   size_mat0=size(maskH0_embed,1)
   size_mat1=size(maskH1_embed,1)
   if(size_mat0/=size_mat1)then
     write(*,*) 'ERROR : embedding, size_mat do not match'
     stop
   endif

   ttt0=0.0
   do i=1,size_mat0
    do j=1,size_mat0
      s  = maskH0_embed(i,j,1)
      t  = maskH0_embed(i,j,2)
      ii = maskH0_embed_atoms(i,j,1)  
      jj = maskH0_embed_atoms(i,j,2) 
      if(pub_node_of_atom(pub_orig_atom(ii))==pub_my_node_id.and. &
       & pub_node_of_atom(pub_orig_atom(jj))==pub_my_node_id) &
       & call sparse_get_element(sig_element,self_energy,s,t)
      ttt0(i,j) = sig_element 
    enddo
   enddo
   call comms_reduce('SUM', ttt0)

   if(pub_my_node_id/=0) then
     deallocate(ttt0)
     return
   endif

   filename='embedding_potentials_'//TRIM(ADJUSTL(toString(is)))
   if(mpi_size>1)then
     filename=TRIM(ADJUSTL(filename))//"_rank"//TRIM(ADJUSTL(toString(mpi_rank)))
   endif

   if(pub_my_node_id==0)then
     if(ien==ffirst)then
      open(unit=1001,file=filename,form='unformatted',STATUS='REPLACE')
     else
      open(unit=1001,file=filename,form='unformatted',position='append',status='old')
     endif
     write(1001)  ien , pub_dmft_points , pub_dmft_temp , size_mat0 , size_mat0
     write(1001)  fermi_e(is)+pub_dmft_chem_shift,energy,ttt0(1:size_mat0,1:size_mat0)
     close(1001)
   endif

   deallocate(ttt0)

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine embedding_read_potential
   use wrappers, only : wrappers_invert_sym_matrix
   implicit none
   integer                :: i,j,s,t
   complex(8)             :: sig_element,tmp,ww_
   real(8)                :: mmuL,mmuR
   integer                :: kk,kkk,ll,lll,jj,iii,jjj,ii1,jj1
   integer                :: l1,l2
   complex(8)             :: www_
   character(30)          :: filename
   complex(8),allocatable :: mat_square(:,:)
   logical                :: checkl,checkr,buildr,buildl

   buildr=.false.
   buildl=.false.

   if(pub_dmft_temp > 0.0_DP)then
    filename='sigma_embedding'
   else
    filename='sigma_embedding_real'
   endif

   inquire(file=trim(adjustl(filename))//'L'//trim(adjustl(toString(is))),exist=checkl)
   inquire(file=trim(adjustl(filename))//'R'//trim(adjustl(toString(is))),exist=checkr)
   check_embed_h=checkr.or.checkl
   if(.not.check_embed_h) return

   if(pub_my_node_id==0.and.mpi_rank==1) write(*,*) 'READING EMBEDDING POTENTIALS, L/R PRESENT ? ',checkl,checkr
   if( (checkr.and..not.checkl) .or. (checkl.and..not.checkr) )then
     write(*,*) 'WARNING : could find file sigma_embedding_(real)L but NOT sigma_embedding_(real)R'
   endif

   if(check_embed)then
     write(*,*) ' ERROR : not supposed to have L-R LEAD FILES'
     write(*,*) '         AND sigma_embedding_(real) file (the latter is for the device only)'
     stop
   endif

   call sparse_scale(ham_embed(is),0.0d0)

   fopenEL=20010+2001*is
   fopenER=10400+3001*is

  if(.not.checkl) goto 1013

   if(ien==ffirst.or.reopen_sigma_file_each_step.or.force_reopen_file) then
     if(mpi_rank==1.and.pub_my_node_id==0) write(*,*) 'opening file : ', trim(adjustl(filename))//'L'//trim(adjustl(toString(is))) 
     open(unit=fopenEL,file=trim(adjustl(filename))//'L'//trim(adjustl(toString(is))),form='unformatted')
     read(fopenEL) size_mat0L,size_indexL
     if(mpi_rank==1.and.pub_my_node_id==0) write(*,*) 'SIZE L : ', size_mat0L
     if(.not.allocated(ttt0L))then
       write(*,*) 'ALLOCATING EMBEDDING ARRAYS LEFT'
       buildl=.true.
       allocate(H0L_(size_mat0L,size_mat0L))
       allocate(H1L_(size_mat0L,size_mat0L))
       allocate(S0L_(size_mat0L,size_mat0L))
       allocate(S1L_(size_mat0L,size_mat0L))
#ifndef savememory
       allocate(SSL_(size_mat0L,size_mat0L,pub_dmft_points))
#else
       allocate(SSL_(size_mat0L,size_mat0L))
#endif
       allocate(ttL_(size_mat0L,size_mat0L))
       allocate(tttL_(size_mat0L,size_mat0L))
#ifdef GPU_SPEEDUP
       allocate( ttLc_(size_mat0L,size_mat0L))
       allocate(tttLc_(size_mat0L,size_mat0L))
#endif
       allocate(GGL_(size_mat0L,size_mat0L))
       allocate(GGLback_(size_mat0L,size_mat0L))
       allocate(ttt0L(size_mat0L,size_mat0L))  
       allocate(maskTLatom(size_mat0L,size_mat0L,2))
       allocate(maskTL(size_mat0L,size_mat0L,2))
       allocate(mask_vectorL(size_indexL))
       allocate(mmuLar(pub_dmft_points),wwLar(pub_dmft_points))
       read(fopenEL) mask_vectorL,H0L_,H1L_,S0L_,S1L_,SSL_,mmuLar,wwLar
     endif
#ifdef savememory
    read(fopenEL) mask_vectorL,H0L_,H1L_,S0L_,S1L_
#endif
   endif

   if(.not.allocated(H0L_).or..not.allocated(SSL_))then
     write(*,*) 'CRITICAL - EMBEDDING - ARRAYS NOT ALLOCATED'
     stop
   endif

   if(verboseall) then
      write(*,*) 'ALLOCATED : ', allocated(mask_vectorL)
      write(*,*) 'size      : ', size(mask_vectorL)
      write(*,*) 'size H0H1 : ', shape(H0L_),shape(H1L_),shape(S0L_),shape(S1L_)
      write(*,*) 'ien       : ', ien
      write(*,*) 'LOGIC     : ', reopen_sigma_file_each_step,force_reopen_file 
   endif

    mmuL=mmuLar(ien)
    ww_=wwLar(ien)

#ifdef savememory
   if(.not.(reopen_sigma_file_each_step.or.force_reopen_file))then
       read(fopenEL) mmuL,SSL_,ww_
    else
     do kk=1,ien
       read(fopenEL) mmuL,SSL_,ww_
     enddo
   endif
#endif

   if(ien==llast.or.reopen_sigma_file_each_step.or.force_reopen_file)then
      close(fopenEL)
   endif

   if(verboseall) write(*,*) 'DYSON EQUATION'
   www_= (ww_-mmuL) + (fermi_e(is)-mmuL) + fermi_e(is)
   if(pub_my_node_id==0)  write(*,*) 'EMBEDDING SHIFT MU L : ', fermi_e(is)-mmuL
   ttt0L  = 0.0
   ttL_   =           H1L_  - www_*          S1L_
   tttL_  = transpose(H1L_) - www_*transpose(S1L_)
#ifdef GPU_SPEEDUP
   ttLc_=ttL_;tttLc_=tttL_;
#endif
   do kk=1,pub_dmft_embed_iter
      if(mpi_rank==1.and.pub_my_node_id==0.and..not.silent) write(*,*) 'ITER EMBED',kk,real(ttt0L(1,1)),aimag(ttt0L(1,1))

#ifndef savememory
      GGL_ = www_*S0L_ - H0L_ - ttt0L - SSL_(:,:,ien)
#else
      GGL_ = www_*S0L_ - H0L_ - ttt0L - SSL_
#endif

#ifdef debug
     GGL_ = (GGL_ + transpose(GGL_))/2.0
     ttt0L=GGL_
#endif

     i=size(GGL_,1); call invert_matc(i,GGL_)

#ifdef debug
     if(pub_my_node_id==0) write(*,*) 'TESTING INVERSION : ', &
            & maxval(abs(matmul(GGL_,ttt0L))),minval(abs(matmul(GGL_,ttt0L)))
#endif

#ifdef GPU_SPEEDUP
      i=size(ttL_,1)
      if(i/=size(GGL_,1))then; write(*,*) 'gpu speedup dim dont match';stop; endif
      if(.not.(kk==1.or.abs(pub_dmft_embed_mix-1.d0)<1.d-4)) GGLback_=ttt0L
#ifdef GPU_SPEEDUP_SINGPREC
      call matmulcuda_c_singleprec_bypass(ttLc_, GGL_ ,ttt0L,i,i,i); GGL_=ttt0L
      call matmulcuda_c_singleprec_bypass(GGL_ ,tttLc_,ttt0L,i,i,i) 
#else
      call matmulcuda_c(ttLc_, GGL_ ,ttt0L,i,i,i); GGL_=ttt0L
      call matmulcuda_c(GGL_ ,tttLc_,ttt0L,i,i,i)
#endif
      if(.not.(kk==1.or.abs(pub_dmft_embed_mix-1.d0)<1.d-4)) ttt0L=pub_dmft_embed_mix*ttt0L+(1.d0-pub_dmft_embed_mix)*GGLback_
#else
     if(kk==1.or.abs(pub_dmft_embed_mix-1.d0)<1.d-4)then
      ttt0L = matmul ( matmul( ttL_, GGL_ ) , tttL_ )
     else
      ttt0L = pub_dmft_embed_mix*matmul ( matmul( ttL_, GGL_ ) , tttL_ )+(1.d0-pub_dmft_embed_mix)*ttt0L
     endif
#endif

     if(pub_my_node_id==0.and.(ien==pub_dmft_points-1.or.ien==pub_dmft_points))then
      write(*,*) 'w=oo, frequ    = ', www_
      write(*,*) 'MAX RE SIGMA_L = ', maxval(abs(real(ttt0L)))
      write(*,*) 'MAX IM SIGMA_L = ', maxval(abs(aimag(ttt0L)))
     endif
#ifdef debug
     if(pub_my_node_id==0) then
       write(*,*) 'DYSON ITERATION L (maxval) : ', kk,maxval(abs(ttt0L))
       write(*,*) 'MAXVAL SELF                : ', maxval(abs(SSL_))
     endif
#endif
   enddo

#ifdef debug
   if(verboseall) write(*,*) 'build mask L'
   i=0
   do j=1,size(mask_vectorL)
    i=i+ngwf_basis%num_on_atom( pub_distr_atom(mask_vectorL(j))  )
   enddo
   if(verboseall) write(*,*) 'TOTAL SIZE OF ARRAYS : ', i
   if( i /= size(maskTL,1)) then
    write(*,*) 'ERROR EMBEDDING'
    write(*,*) 'cumulated dimensions of atoms obtained from device indices : ', i
    write(*,*) 'shape of mask array obtained from file                     : ', shape(maskTL)
    write(*,*) 'mask_vector : ', mask_vectorL
    stop
   endif  
#endif

  if(buildl)then
   maskTL=0
   kk=size_indexL
   do i=1,kk
    do j=1,kk
      ll=0
      do jj=1,pub_distr_atom(mask_vectorL(i))-1
        ll=ll+ngwf_basis%num_on_atom( jj )
      enddo
      lll=0
      do jj=1,pub_distr_atom(mask_vectorL(j))-1
        lll=lll+ngwf_basis%num_on_atom( jj )
      enddo   
      iii=0
      do jj=1,i-1
        iii=iii+ngwf_basis%num_on_atom( pub_distr_atom(mask_vectorL(jj)) )
      enddo
      jjj=0
      do jj=1,j-1
        jjj=jjj+ngwf_basis%num_on_atom( pub_distr_atom(mask_vectorL(jj)) )
      enddo
      do l1=1,ngwf_basis%num_on_atom( pub_distr_atom(mask_vectorL( i  )) )
       do l2=1,ngwf_basis%num_on_atom( pub_distr_atom(mask_vectorL( j  )) )
         maskTL( iii+l1 , jjj+l2 , 1 ) = ll  + l1
         maskTL( iii+l1 , jjj+l2 , 2 ) = lll + l2
         maskTLatom( iii+l1 , jjj+l2 , 1 ) = pub_distr_atom(mask_vectorL( i  ))
         maskTLatom( iii+l1 , jjj+l2 , 2 ) = pub_distr_atom(mask_vectorL( j  ))
       enddo
      enddo
    enddo 
   enddo
  endif

   do i=1,size_mat0L
    do j=1,size_mat0L
      s = maskTL(i,j,1)
      t = maskTL(i,j,2)
      ii1=maskTLatom(i,j,1)
      jj1=maskTLatom(i,j,2)
      if(pub_node_of_atom(pub_orig_atom(ii1))==pub_my_node_id.and..not. &
       & pub_node_of_atom(pub_orig_atom(jj1))==pub_my_node_id) then
       write(*,*) 'ERROR : one atom one my node, but not the other one'
       stop
      endif
      if(.not.pub_node_of_atom(pub_orig_atom(ii1))==pub_my_node_id.and.&
            & pub_node_of_atom(pub_orig_atom(jj1))==pub_my_node_id) then
       write(*,*) 'ERROR : one atom one my node, but not the other one'
       stop
      endif
      if(pub_node_of_atom(pub_orig_atom(ii1))==pub_my_node_id.and. &
       & pub_node_of_atom(pub_orig_atom(jj1))==pub_my_node_id) then
        sig_element=ttt0L(i,j)
        call sparse_put_element(sig_element,ham_embed(is),s,t)
      endif
    enddo
   enddo

  1013 continue
  if(.not.checkr) return
 
  if(ien==ffirst.or.reopen_sigma_file_each_step.or.force_reopen_file) then
     open(unit=fopenER,file=trim(adjustl(filename))//'R'//trim(adjustl(toString(is))),form='unformatted')
     read(fopenER) size_mat0R,size_indexR
     if(mpi_rank==1.and.pub_my_node_id==0) write(*,*) 'SIZE R / size index : ', size_mat0R,size_indexR
     if(.not.allocated(ttt0R))then
       buildr=.true.
       write(*,*) 'ALLOCATING EMBEDDING ARRAYS RIGHT'
       allocate(H0R_(size_mat0R,size_mat0R))
       allocate(H1R_(size_mat0R,size_mat0R))
       allocate(S0R_(size_mat0R,size_mat0R))
       allocate(S1R_(size_mat0R,size_mat0R))
#ifndef savememory
       allocate(SSR_(size_mat0R,size_mat0R,pub_dmft_points))
#else
       allocate(SSR_(size_mat0R,size_mat0R))
#endif
       allocate(ttR_(size_mat0R,size_mat0R))
       allocate(tttR_(size_mat0R,size_mat0R))
#ifdef GPU_SPEEDUP
       allocate( ttRc_(size_mat0R,size_mat0R))
       allocate(tttRc_(size_mat0R,size_mat0R))
#endif
       allocate(GGR_(size_mat0R,size_mat0R))
       allocate(GGRback_(size_mat0R,size_mat0R))
       allocate(ttt0R(size_mat0R,size_mat0R))
       allocate(maskTR(size_mat0R,size_mat0R,2))
       allocate(maskTRatom(size_mat0R,size_mat0R,2))
       allocate(mask_vectorR(size_indexR))
       allocate(mmuRar(pub_dmft_points),wwRar(pub_dmft_points))
       read(fopenER) mask_vectorR,H0R_,H1R_,S0R_,S1R_,SSR_,mmuRar,wwRar
     endif
#ifdef savememory
      read(fopenER) mask_vectorR,H0R_,H1R_,S0R_,S1R_
#endif
   endif

   if(.not.allocated(H0R_).or..not.allocated(SSR_))then
     write(*,*) 'CRITICAL - EMBEDDING - ARRAYS NOT ALLOCATED'
     stop
   endif

   if(verboseall) then
      write(*,*) 'ALLOCATED : ', allocated(mask_vectorR)
      write(*,*) 'size      : ', size(mask_vectorR)
      write(*,*) 'size H0H1 : ', shape(H0R_),shape(H1R_),shape(S0R_),shape(S1R_)
      write(*,*) 'ien       : ', ien
      write(*,*) 'LOGIC     : ', reopen_sigma_file_each_step,force_reopen_file
   endif

    mmuR=mmuRar(ien) 
    ww_=wwRar(ien)

#ifdef savememory
   if(.not.(reopen_sigma_file_each_step.or.force_reopen_file))then
       read(fopenER) mmuR,SSR_,ww_
    else
     do kk=1,ien
       read(fopenER) mmuR,SSR_,ww_
     enddo
   endif
#endif

   if(ien==llast.or.reopen_sigma_file_each_step.or.force_reopen_file)then
      close(fopenER)
   endif

#ifdef debug
   write(*,*) 'ASSYMETRY H0 : ', maxval(abs(H0R_-transpose(H0R_)))
   write(*,*) 'ASSYMETRY S0 : ', maxval(abs(S0R_-transpose(S0R_)))
   write(*,*) 'ASSYMETRY SS : ', maxval(abs(SSR_(:,:,ien)-transpose(SSR_(:,:,ien))))
   write(*,*) 'SYMETRY   H0 : ', maxval(abs(H0R_+transpose(H0R_)))
   write(*,*) 'SYMETRY   S0 : ', maxval(abs(S0R_+transpose(S0R_)))
   write(*,*) 'SYMETRY   SS : ', maxval(abs(SSR_(:,:,ien)+transpose(SSR_(:,:,ien))))
#endif

   if(verboseall) write(*,*) 'DYSON EQUATION'
   www_= (ww_-mmuR) + (fermi_e(is)-mmuR) + fermi_e(is)
   if(pub_my_node_id==0)  write(*,*) 'EMBEDDING SHIFT MU R : ', fermi_e(is)-mmuR
   ttt0R  =  0.0
   ttR_   =           H1R_  - www_*          S1R_
   tttR_  = transpose(H1R_) - www_*transpose(S1R_)
#ifdef GPU_SPEEDUP
   ttRc_=ttR_;tttRc_=tttR_;
#endif
   do kk=1,pub_dmft_embed_iter
      if(mpi_rank==1.and.pub_my_node_id==0.and..not.silent) write(*,*) 'ITER EMBED',kk,real(ttt0R(1,1)),aimag(ttt0R(1,1))
#ifdef savememory
      GGR_ = www_*S0R_ - H0R_ - ttt0R - SSR_
#else
      GGR_ = www_*S0R_ - H0R_ - ttt0R - SSR_(:,:,ien)
#endif
#ifdef debug
      GGR_ = (GGR_+transpose(GGR_))/2.d0
#endif

      i=size(GGR_,1); call invert_matc(i,GGR_)

#ifdef GPU_SPEEDUP
      i=size(ttR_,1)
      if(i/=size(GGR_,1))then; write(*,*) 'gpu speedup dim dont match R';stop; endif
      if(.not.(kk==1.or.abs(pub_dmft_embed_mix-1.d0)<1.d-4)) GGRback_=ttt0R
#ifdef GPU_SPEEDUP_SINGPREC
      call matmulcuda_c_singleprec_bypass(ttRc_, GGR_ ,ttt0R,i,i,i); GGR_=ttt0R
      call matmulcuda_c_singleprec_bypass(GGR_ ,tttRc_,ttt0R,i,i,i)
#else
      call matmulcuda_c(ttRc_, GGR_ ,ttt0R,i,i,i); GGR_=ttt0R
      call matmulcuda_c(GGR_ ,tttRc_,ttt0R,i,i,i)
#endif
      if(.not.(kk==1.or.abs(pub_dmft_embed_mix-1.d0)<1.d-4)) ttt0R=pub_dmft_embed_mix*ttt0R+(1.d0-pub_dmft_embed_mix)*GGRback_
#else
     if(kk==1.or.abs(pub_dmft_embed_mix-1.d0)<1.d-4)then
      ttt0R = matmul ( matmul( ttR_, GGR_ ) , tttR_ )
     else
      ttt0R = pub_dmft_embed_mix*matmul ( matmul( ttR_, GGR_ ) , tttR_ )+(1.d0-pub_dmft_embed_mix)*ttt0R
     endif
#endif

     if(pub_my_node_id==0.and.(ien==pub_dmft_points-1.or.ien==pub_dmft_points))then
       write(*,*) 'w=oo, frequ    = ', www_
       write(*,*) 'MAX RE SIGMA_R = ', maxval(abs(real(ttt0R)))
       write(*,*) 'MAX IM SIGMA_R = ', maxval(abs(aimag(ttt0R)))
     endif
#ifdef debug
     if(pub_my_node_id==0) then
       write(*,*) 'DYSON ITERATION R (max) : ', kk,maxval(abs(ttt0R))
       write(*,*) '  G (max)               : ', maxval(abs(GGR_))
     endif
#endif
   enddo

#ifdef debug
   if(verboseall) write(*,*) 'build mask R'
   i=0; do j=1,size(mask_vectorR)
   i=i+ngwf_basis%num_on_atom( pub_distr_atom(mask_vectorR(j))  )
   enddo
   if(verboseall)  write(*,*) 'TOTAL SIZE OF ARRAYS : ', i
   if( i /= size(maskTR,1)) then
    write(*,*) 'ERROR EMBEDDING'
    write(*,*) 'cumulated dimensions of atoms obtained from device indices : ', i
    write(*,*) 'shape of mask array obtained from file                     : ', shape(maskTR)
    write(*,*) 'mask_vector : ', mask_vectorR
    stop
   endif 
#endif


!erase ham by boundary
  if(buildr)then
   kk  = size_indexR
   kkk = size_indexL
   do i=1,kk
    do j=1,kkk
      if(mask_vectorR(i)==mask_vectorL(j))then
        write(*,*) 'ERROR : same site cannot be both in Right and Left leads'
        write(*,*) 'site : ', i
        write(*,*) 'site : ', j
        stop
      endif
      ll=0
      do jj=1,pub_distr_atom(mask_vectorR(i))-1
        ll=ll+ngwf_basis%num_on_atom( jj )
      enddo
      lll=0
      do jj=1,pub_distr_atom(mask_vectorL(j))-1
        lll=lll+ngwf_basis%num_on_atom( jj )
      enddo
      do l1=1,ngwf_basis%num_on_atom( pub_distr_atom(mask_vectorR( i  )) )
       do l2=1,ngwf_basis%num_on_atom( pub_distr_atom(mask_vectorL( j  )) )
         s=ll +l1  
         t=lll+l2  
         call sparse_put_element(0.d0,ham(is),s,t)
         call sparse_put_element(0.d0,ham(is),t,s)
       enddo
      enddo
    enddo
   enddo
  endif
!end erase ham boundary

  if(buildr)then
   if(verboseall)  write(*,*) 'build mask R '
   maskTR=0
   kk=size_indexR
   do i=1,kk
    do j=1,kk
      ll=0
      do jj=1,pub_distr_atom(mask_vectorR(i))-1
        ll=ll+ngwf_basis%num_on_atom( jj )
      enddo
      lll=0
      do jj=1,pub_distr_atom(mask_vectorR(j))-1
        lll=lll+ngwf_basis%num_on_atom( jj )
      enddo
      iii=0
      do jj=1,i-1
        iii=iii+ngwf_basis%num_on_atom( pub_distr_atom(mask_vectorR(jj)) )
      enddo
      jjj=0
      do jj=1,j-1
        jjj=jjj+ngwf_basis%num_on_atom( pub_distr_atom(mask_vectorR(jj)) )
      enddo
      do l1=1,ngwf_basis%num_on_atom( pub_distr_atom(mask_vectorR( i  )) )
       do l2=1,ngwf_basis%num_on_atom( pub_distr_atom(mask_vectorR( j  )) )
         maskTR( iii+l1 , jjj+l2 , 1 ) = ll  + l1
         maskTR( iii+l1 , jjj+l2 , 2 ) = lll + l2
         maskTRatom( iii+l1 , jjj+l2 , 1 ) = pub_distr_atom(mask_vectorR( i  ))
         maskTRatom( iii+l1 , jjj+l2 , 2 ) = pub_distr_atom(mask_vectorR( j  ))
       enddo
      enddo
    enddo
   enddo
  endif

   do i=1,size_mat0R
    do j=1,size_mat0R
      s = maskTR(i,j,1)
      t = maskTR(i,j,2)
      ii1=maskTRatom(i,j,1)
      jj1=maskTRatom(i,j,2)
      if(pub_node_of_atom(pub_orig_atom(ii1))==pub_my_node_id.and..not. &
       & pub_node_of_atom(pub_orig_atom(jj1))==pub_my_node_id) then
       write(*,*) 'ERROR : one atom one my node, but not the other one'
       stop
      endif
      if(.not.pub_node_of_atom(pub_orig_atom(ii1))==pub_my_node_id.and.&
            & pub_node_of_atom(pub_orig_atom(jj1))==pub_my_node_id) then
       write(*,*) 'ERROR : one atom one my node, but not the other one'
       stop
      endif
      if(pub_node_of_atom(pub_orig_atom(ii1))==pub_my_node_id.and. &
       & pub_node_of_atom(pub_orig_atom(jj1))==pub_my_node_id) then
        sig_element=ttt0R(i,j)
        tmp=0.d0
        call sparse_get_element(        tmp,ham_embed(is),s,t)
#ifdef debug
        if(abs(tmp)>1.d-5)then
          write(*,*) 'SITE 1 : ', pub_orig_atom(ii1) 
          write(*,*) 'SITE 2 : ', pub_orig_atom(jj1)
          write(*,*) 'connected to the L and R leads...'
          write(*,*) 'amplitude : ', tmp
          write(*,*) 'should not happen....'
          stop
        endif
#endif
        sig_element=sig_element+tmp
        call sparse_put_element(sig_element,ham_embed(is),s,t)
      endif
    enddo
   enddo

#ifdef debug
  allocate( mat_square(ngwf_basis%num,ngwf_basis%num) )
  call sparse_convert(mat_square, ham_embed(is) )
  if(maxval(abs( mat_square - transpose(mat_square) )) >1.d-4)then
   write(*,*) 'ham_embed not symmetric'
   write(*,*) 'max diff : ', maxval(abs( mat_square - transpose(mat_square) ))
   stop
  endif
  deallocate( mat_square )
#endif

   return
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine fix_fermi_and_total_energies
   use wrappers, only : wrappers_dsygv_lt
#ifdef GPU_SPEEDUP
  use fortran_cuda, only : eigenvector_gen_cuda_r
#endif
   implicit none
   logical :: check_chem
   integer :: ii

    inquire(file='store_tot_energy',exist=check)
    if(check)then
      open(unit=5050,file='store_tot_energy')
      read(5050,*) total_energy,lhxc_energy, hubbard_energy
      write(*,*) ' tot energy : ', total_energy
      write(*,*) ' lhxc E     : ', lhxc_energy
      write(*,*) ' hub E      : ', hubbard_energy
      close(5050)
    endif

    if(present(dmft_energy_cor)) total_energy=0.d0

    dmft_energy_cor_=0.d0

    ii=ngwf_basis%num

    if(present(dmft_energy_cor))then
      if(nkpoints>1)then
       write(*,*) 'ERROR : DMFT ENERGY COR and k-points not yet possible at the moment'
       STOP
      endif
      do is=1,pub_cell%num_spins
       if(.not.full_prec_dmft_energy.and..not.(use_gpu_onlyme.and.use_gpu_eigenvectors))then
         call eigenvec_fullmat(ii,ham(is),overlap,matrixAA(is,:,:),vectorAA(is,:))
       else
         call sparse_convert(matrixAA(is,:,:),ham(is))
         call sparse_convert(matrixOverlap,overlap)
         if(.not.(use_gpu_onlyme.and.use_gpu_eigenvectors))then
          call wrappers_dsygv_lt(eigenvectors=matrixAA(is,:,:),eigenvalues=vectorAA(is,:),overlap_square=matrixOverlap,num=ii)
         else
#ifdef GPU_SPEEDUP
          call eigenvector_gen_cuda_r(ii,matrixAA(is,:,:),matrixOverlap,vectorAA(is,:),tmpr,.false.)
#else
          write(*,*) 'error gpu 1' ; stop
#endif
          matrixAA(is,:,:)=tmpr
         endif
       endif
       E_sum_eigenvalues_DFT(is)=sum(vectorAA(is,1:n_occ(is)))
       if(pub_my_node_id==0)write(*,*) 'SUM EIGENVALUES OF H : ', E_sum_eigenvalues_DFT(is)
       if(pub_my_node_id==0)write(*,*) 'NOCC                 : ', n_occ(is)
       eigen_en(:,is)=vectorAA(is,:)
       fermi_e(is)=0.5_DP*(eigen_en(n_occ(is),is)+eigen_en(n_occ(is)+1,is))
      enddo
      if(pub_cell%num_spins==1)then
        !paramagnetic
        if(size(fermi_e)   ==2)    fermi_e(2)=fermi_e(1)
        if(size(eigen_en,2)==2) eigen_en(:,2)=eigen_en(:,1)
      endif
    else
                                 E_sum_eigenvalues_DFT(1) = sum(eigen_en(1:n_occ(1),1))
       if(pub_cell%num_spins==2) E_sum_eigenvalues_DFT(2) = sum(eigen_en(1:n_occ(2),2))
       do is=1,pub_cell%num_spins
         fermi_e(is) = 0.5_DP*(eigen_en(n_occ(is),is)+eigen_en(n_occ(is)+1,is))
       enddo
    endif

    if(pub_my_node_id==0.and..not.silent)then
      write(*,*) '=========================================================='
      write(*,*) 'FERMI E IS IN THE Mid-gap, shift it with pub_dmft_chem_pot'
      write(*,*) 'EIGEN_EN SHAPE      : ', shape(eigen_en)
      write(*,*) 'MPI RANK            : ', mpi_rank
      write(*,*) 'MPI SIZE            : ', mpi_size
      write(*,*) 'number of states    : ', n_occ(:)
      do is=1,num_spins
        write(*,*) 'eigen E last state  : ', eigen_en(n_occ(is)  ,is)
        write(*,*) 'eigen E last +1     : ', eigen_en(n_occ(is)+1,is)
      enddo
      write(*,*) 'sum occ eigenvalues : ', E_sum_eigenvalues_DFT
    endif
 
    inquire(file='chem.potential.nmu.iter1',exist=check_chem)
    if(check_chem)then
      write(*,*) 'FILE WITH CHEM POT SPIN UP PRESENT READ IT'
      open(unit=81119,file='chem.potential.nmu.iter1')
      read(81119,*) fermi_e(1)
      close(81119)
      if(pub_cell%num_spins==2)then
         inquire(file='chem.potential.nmu.iter2',exist=check_chem) 
         if(check_chem)then
           write(*,*) 'FILE WITH CHEM POT SPIN DOWN PRESENT READ IT'
           open(unit=81119,file='chem.potential.nmu.iter2')
           read(81119,*) fermi_e(2)
           close(81119)
         else
           fermi_e(2)=fermi_e(1)
         endif
      endif 
      write(*,*) 'OBTAINED FERMI : ', fermi_e(1:pub_cell%num_spins)
      write(*,*) 'DONE'
    endif

    if(pub_my_node_id==0)write(*,*) 'FERMI E (UP/DN)     : ', fermi_e
    if(pub_my_node_id==0)write(*,*) '=========================================================='

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                         !----------------!

   subroutine fix_chemical_potential_loop
   implicit none
          target_N = sum(n_occ(1:num_spins)) + pub_dmft_doping
          if(pub_my_node_id==0) write(*,*) 'target N : ', target_N
          nmu_step = nmu_step+1
          if(nmu_step>1)then
           if(abs( vectorNmu(nmu_step-1) - target_N )< pub_dmft_mu_diff_max)then
             vectorNmu(nmu-1)=vectorNmu(nmu_step-1)
             vectorMu(nmu-1)=vectorMu(nmu_step-1)
             nmu_step=nmu 
             chem_shift=0.0 
           endif
          endif
          chem_left  = minval(fermi_e(1:num_spins)) - chem_left_shift
          chem_right = maxval(fermi_e(1:num_spins)) + chem_right_shift
          if(pub_my_node_id==0)write(*,*) 'chem pot brackets   : ', chem_left,chem_right
          vectorNmu(nmu_step)=0.d0
          do is=1,num_spins
             fermi_e(is)=fermi_e(is)+chem_shift
             if(fermi_e(is)<chem_left )  fermi_e(is)=chem_left
             if(fermi_e(is)>chem_right)  fermi_e(is)=chem_right
          enddo
          vectorMu(nmu_step)=sum(fermi_e(1:num_spins))/dble(num_spins)
   end subroutine
           
                         !----------------!

   subroutine fix_chemical_potential_loop_spin
   implicit none
          target_N = n_occ(is)+pub_dmft_doping/2.d0
          if(pub_my_node_id==0) write(*,*) 'target N : ', target_N
          nmu_step = nmu_step+1
          if(nmu_step>1)then
           if(abs( vectorNmu(nmu_step-1) - target_N )< pub_dmft_mu_diff_max)then
             vectorNmu(nmu-1)=vectorNmu(nmu_step-1)
              vectorMu(nmu-1)=vectorMu(nmu_step-1)
             nmu_step=nmu
             chem_shift=0.0
           endif
          endif
          chem_left  = fermi_e(is) - chem_left_shift
          chem_right = fermi_e(is) + chem_right_shift
          if(pub_my_node_id==0)write(*,*) 'chem pot brackets   : ', chem_left,chem_right
          vectorNmu(nmu_step)=0.d0
          fermi_e(is)=fermi_e(is)+chem_shift
          if(fermi_e(is)<chem_left )  fermi_e(is)=chem_left
          if(fermi_e(is)>chem_right)  fermi_e(is)=chem_right
          vectorMu(nmu_step)=fermi_e(is)
   end subroutine

                         !----------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   pure integer function hub_atom_hub_projs(hat)
   implicit none
       integer            :: theatom
       integer,intent(in) :: hat
       theatom            = pub_distr_atom(hub%orig(hat))
       hub_atom_hub_projs = hub_proj_basis%num_on_atom(theatom)
   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine define_my_frequencies
 implicit none
 integer :: iistart,iistep

    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'START TO BUILD FREQUENCIES'
    if(allocated(myfrequencies)) deallocate(myfrequencies)
    allocate(myfrequencies(pub_dmft_points))
    myfrequencies=0 ; myfrequencies_numb=0
    INQUIRE(file="mask_frequ"//trim(adjustl(toString(mpi_rank))),EXIST=check_profile)
    if(check_profile.and..not.pub_dmft_lin_scaling.and..not.split.and..not.pub_dmft_splitk)then
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'OPENING FREQUENCIES PROFILE'
       open(unit=1716,file="mask_frequ"//trim(adjustl(toString(mpi_rank))))
       myfrequencies_numb=0
       read(1716,*) profile_ratio
       if(abs(profile_ratio)<1.d-3)then
       do
          read(1716,*,end=896) iien
          myfrequencies_numb=myfrequencies_numb+1
          myfrequencies(myfrequencies_numb)=iien
       enddo
896    continue
       else
        if(mpi_size/=2)then
          write(*,*) 'first line in mask_frequ is only for 2 gpus'
          write(*,*) 'it gives the ratio between gpu1/gpu2 number of frequs'
          write(*,*) 'please put 0.00 in the first line of the mask_frequ files'
          stop
        endif
        if(mpi_rank==1)then
          myfrequencies_numb=                  NINT(profile_ratio*dble(pub_dmft_points))
          myfrequencies(1:myfrequencies_numb)=(/( i, i=1,myfrequencies_numb )/)
        else
          myfrequencies_numb=pub_dmft_points - NINT(profile_ratio*dble(pub_dmft_points))
          myfrequencies(1:myfrequencies_numb)=(/( i, i=pub_dmft_points-myfrequencies_numb+1,pub_dmft_points )/)
        endif
       endif
       close(1716)
    else
       is=1
       if(.not.pub_dmft_splitk)then
         iistart=(mpi_rank-1); iistep=mpi_size
       else
         iistart=0; iistep=1
       endif
       do ien=iistart+1,pub_dmft_points,iistep
         call define_energy
         if( .not. restrict_window .or. check_embed .or. &
               &        ( real(energy) - (fermi_e(is)+pub_dmft_chem_shift) > pub_dmft_dos_min   &
               &  .and.   real(energy) - (fermi_e(is)+pub_dmft_chem_shift) < pub_dmft_dos_max ) ) then
              myfrequencies_numb=myfrequencies_numb+1
              myfrequencies(myfrequencies_numb)=ien
         endif
       enddo
    endif

    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'DONE, DEFINE FIRST/LAST FREQUENCIES'

    ffirst = myfrequencies(1)
    llast  = myfrequencies(myfrequencies_numb)

   if(.not.silent.and.pub_my_node_id==0)then 
    write(*,*) '=========================================================='
    write(*,*) 'FIRST FREQU   :  ', ffirst
    write(*,*) 'LAST  FREQU   :  ', llast
    write(*,*) 'MY RANK       :  ', mpi_rank
    write(*,*) 'my frequs     :  ', myfrequencies(1:myfrequencies_numb)
    write(*,*) 'total number  :  ', myfrequencies_numb
    write(*,*) '=========================================================='
   endif
 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine check_if_we_can_catch_hamiltonian_in_a_file
 implicit none

    INQUIRE(file='store_eigen',EXIST=check)

    call comms_barrier
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'read-write H'
    call read_write_hamiltonian
    if(.not.silent.and.pub_my_node_id==0)write(*,*) 'dump some info'
    call dump_some_info_on_h_atom

 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine catch_input_arguments_and_set_module_variables
 implicit none
 logical :: checkit
 integer :: iiki

   !BUG : or a good reason to have them?
    inv_projectors=.false.
#ifdef debug
    silent=.false.
    inv_projectors=.true.
#endif
    if(inv_projectors)then
     pub_dmft_invert_overlap=.true.
    endif


    if(pub_dmft_kernel==0.and.present(dmft_energy_cor))then
      write(*,*) 'computing the dmft energy but not the dmft kernel'
      write(*,*) 'inconsistent flags... stop'
      stop
    endif

    nmu=pub_dmft_nmu_loop
    if(nmu==0.and.(pub_dmft_kernel/=0.or.present(dmft_energy_cor)))then
      nmu=1
    endif

    if(pub_dmft_lin_scaling)then
     nmu=0; pub_dmft_kernel=0 
#ifndef LINEAR_SCALING_INV
     write(*,*) 'ERROR : if flag pub_dmft_lin_scaling set to T, then please compile with preprocessing flag LINEAR_SCALING_INV'
     stop
#endif
    endif

    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'sparse create overlap inv-overlap'
    if(pub_dmft_kernel/=0)then
     do iiki=1,pub_cell%num_spins 
      call sparse_create(     denskern_tmp(iiki) , denskern(iiki), iscmplx=.false. )
      call sparse_scale(      denskern_tmp(iiki) , 0.d0                            )
      call sparse_create( pur_denskern_tmp(iiki) , denskern(iiki), iscmplx=.false. )
     enddo
    endif
    call sparse_create(             overlap ,     overlap_input, iscmplx=.false. )
    call sparse_create(         inv_overlap , inv_overlap_input, iscmplx=.false. )
    call sparse_create( hub_overlap_matrix  ,   rep%hub_overlap, iscmplx=.false. )
    call sparse_create( hub_overlap_tmatrix , rep%hub_overlap_t, iscmplx=.false. )
    do iiki=1,pub_cell%num_spins
      call sparse_create( ham_embed(iiki) , ham(iiki) , iscmplx=.true. )
    enddo

    if(nkpoints>1)then
      call sparse_create(  ovK  , overlap    , iscmplx=.true. )
      call sparse_create(  hamK , ham(1)     , iscmplx=.true. )
      call sparse_create( invK  , inv_overlap, iscmplx=.true. )
      if(pub_dmft_kernel/=0)then
       call sparse_create(denskern_tmpc,denskern(1),iscmplx=.true.)
      endif
    endif

    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'sparse copy input data'
    inquire(file='store_overlap',exist=checkit)
    if(.not.checkit)then 
    call sparse_copy(               overlap ,     overlap_input  )
    endif
    inquire(file='store_inv_overlap',exist=checkit)
    if(.not.checkit)then 
    call sparse_copy(           inv_overlap , inv_overlap_input  )
    endif
    inquire(file='store_hub_overlap_matrix',exist=checkit)
    if(.not.checkit)then 
    call sparse_copy(   hub_overlap_matrix  ,   rep%hub_overlap  )
    endif 
    inquire(file='store_hub_overlap_tmatrix',exist=checkit)
    if(.not.checkit)then 
    call sparse_copy(   hub_overlap_tmatrix , rep%hub_overlap_t  )
    endif

                              num_spins=pub_cell%num_spins
    if(pub_dmft_paramagnetic) num_spins=1
    if(.not.silent.and.pub_my_node_id==0) write(*,*) '=====NUM SPINS=====',num_spins

    inquire(file='gpu_max',exist=check)
    if(check)then
     if(pub_my_node_id==0) write(*,*) 'gpu_max file present'
     open(unit=50000,file='gpu_max')
     read(50000,*) gpu_max
     close(50000)
    else
     gpu_max=8
     write(*,*) 'WARNING gpu_max file not present, assuming there are 8 GPUS'
    endif

    if(pub_dmft_points<=nspecial_frequ+1)then
      write(*,*) 'please use at least [x] dmft points (pub_dmft_points)  : ', nspecial_frequ+2
      write(*,*) 'The last [x] frequencies are used for testing purposes : ', nspecial_frequ
      stop
    endif

    cluster=.false.
    INQUIRE(file="mask_dimer",EXIST=check)
    call comms_barrier
    cluster=check
    if(cluster.and.pub_total_num_nodes/=1)then
      if(pub_my_node_id==0) write(*,*) 'WARNING, dimer mode was removed because of mpi use'
      cluster=.false.
      check=.false.
    endif
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'FILE MASK_DIMER EXISTS?',check

 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(8) function realfrequ(ien)
integer :: ien
 realfrequ = en_start + real(ien-1,kind=DP) * en_step
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(8) function smearing(ien)
integer :: ien
real(8) :: ww
   ww       = (realfrequ(ien)-pub_dmft_smear_shift)
   smearing =  pub_dmft_smear  
   smearing =  smearing + pub_dmft_smear_eta / ( DEXPc(-(ww-pub_dmft_smear_w)/pub_dmft_smear_T) + 1.d0 ) 
   smearing =  smearing + pub_dmft_smear_eta / ( DEXPc( (ww+pub_dmft_smear_w)/pub_dmft_smear_T) + 1.d0 )
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine define_energy_zero
implicit none
integer :: u(1),ii
     if (pub_dmft_temp > 0.0_DP) then
         energyzero=1
     else
         u=minloc(abs(   (/(   en_start+real(ii-1,kind=DP)*en_step, ii=1,pub_dmft_points-nspecial_frequ  )/)  ))
         energyzero=u(1) 
     endif
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine define_energy

          if (pub_dmft_temp > 0.0_DP) then
             !-----------------------!
             ! Matsubara frequencies !
             !-----------------------!
             energy=cmplx( (fermi_e(is)+pub_dmft_chem_shift) , en_start+real(ien-1,kind=DP)*en_step,kind=DP)
             if(ien==pub_dmft_points    &
          & .or.ien==pub_dmft_points-1) &
          &                             energy=cmplx( fermi_e(is)+pub_dmft_chem_shift, pub_dmft_free_green_frequ,kind=DP)
             if(ien==pub_dmft_points-2) energy=cmplx( fermi_e(is)+pub_dmft_chem_shift, 0.d0                     ,kind=DP)
             if(ien==pub_dmft_points-3) energy=cmplx( 1.d0                           , 0.d0                     ,kind=DP)
             if(ien==pub_dmft_points-6) energy=cmplx( 1.d0                           , 0.d0                     ,kind=DP)
             if(ien==pub_dmft_points-4) energy=cmplx( 0.d0                           , 0.d0                     ,kind=DP)
             if(ien==pub_dmft_points-5) energy=cmplx( 0.d0                           , 0.d0                     ,kind=DP)
          else
             !------------------!
             ! Real frequencies !
             !------------------!
             energy = cmplx( (fermi_e(is)+pub_dmft_chem_shift) + realfrequ(ien)     , smearing(ien) , kind=DP)
             if(ien==pub_dmft_points    &
          & .or.ien==pub_dmft_points-1) energy=cmplx(  pub_dmft_free_green_frequ, smearing(ien)      ,kind=DP)
             if(ien==pub_dmft_points-2) energy=cmplx(  fermi_e(is)+pub_dmft_chem_shift  , 0.d0       ,kind=DP)
             if(ien==pub_dmft_points-3) energy=cmplx(  1.d0                         , 0.d0           ,kind=DP)
             if(ien==pub_dmft_points-6) energy=cmplx(  1.d0                         , 0.d0           ,kind=DP)
             if(ien==pub_dmft_points-4) energy=cmplx(  0.d0                         , 0.d0           ,kind=DP)
             if(ien==pub_dmft_points-5) energy=cmplx(  0.d0                         , 0.d0           ,kind=DP)
             if(ien==1)                 energy=cmplx( -pub_dmft_free_green_frequ, smearing(ien)      ,kind=DP)
          endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine build_projected_green

!===PROJECTED GREEN===!

         if(pub_my_node_id==0) write(*,*) 'DISTRIBUTE GF AND DUMP IN FILES' 

         if(.not.split) call comms_barrier

!PART 1

        if(split)then
         i=pub_hub_nat
        else
         i=pub_num_hub_atoms_on_node(pub_my_node_id)
        endif

        if(.not.allocated(site_greenf_buffer)) allocate(site_greenf_buffer(i,i,12,12))
        if(size(site_greenf_buffer,1)/=i)then
          deallocate(site_greenf_buffer); allocate(site_greenf_buffer(i,i,12,12))
        endif

        site_greenf_buffer(:,:,:,:)=0.0_DP

        if(ien==1.and..not.silent) write(*,*) 'Loop over Hubbard atoms on my node'

        call extract_green_function_elements

!PART 2
        if(ien==1.and..not.silent) write(*,*) 'dump buffer into file'

        call write_down_green_function

        if(ien==1.and..not.silent) write(*,*) 'deallocating green'

        deallocate(site_greenf_buffer,stat=ierr)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine upfold_self_energy(icall)
use sparse, only: sparse_convert,sparse_num_cols,sparse_num_rows
implicit none
integer                :: nn,nnn
complex(8),allocatable :: mat_square(:,:)
integer,optional       :: icall

           nn=ngwf_basis%num
          nnn=NINT(sparse_num_rows(hub_overlap_tmatrix_cplx))

          if(verboseall) write(*,*) 'build up self energy in NGWFs'
          if(verboseall) write(*,*) 'Find NGWF representation of self-energy operator'

          if(nkpoints>1)then

           if(same_self_for_all)then
              if(.not.present(icall))then
                write(*,*) 'ERROR icall should be present'
                stop
              endif
              if(use_gpu_some_only.and.mpi_size>1.and..not.breaknfscom)then
                allocate(mat_square(nnn,nnn))
                if( mpi_rank/=1)then ; mat_square=0.0 ; goto 6163; endif
              endif
           endif

       !--------------------------------------------------------!
       ! Sigma_imp -> invOk x O_bar x Sigma_imp x O_bar x invOk !
       !--------------------------------------------------------!

#ifndef GPU_SPEEDUP_MATMUL 
             call sparse_product( selfK1 , O_bar           , invOk  ) ! K1= O_bar x invOk
#ifdef debug
             write(*,*) 'SHOWING O_bar x O_k'
             call sparse_show_matrix(selfK1,show_elems=.true.) 
#endif
             call sparse_product( selfK2 , w_self_energy_v , selfK1 ) ! K2= Sigma_imp x O_bar x invOk
             call sparse_product( selfK1 , invOk           , O_bar  ) ! K1= invOk x O_bar
             call sparse_product( selfK3 , selfK1          , selfK2 ) ! K3= K1 x K2
#else
           if(use_gpu_onlyme)then
             call matmul_in_dense_format( selfK1 , O_bar           , invOk  ) ! K1= O_bar x invOk
#ifdef debug
             write(*,*) 'SHOWING O_bar x O_k'
             call sparse_show_matrix(selfK1,show_elems=.true.) 
#endif
             call matmul_in_dense_format( selfK2 , w_self_energy_v , selfK1 ) ! K2= Sigma_imp x O_bar x invOk
             call matmul_in_dense_format( selfK1 , invOk           , O_bar  ) ! K1= invOk x O_bar
             call matmul_in_dense_format( selfK3 , selfK1          , selfK2 ) ! K3= K1 x K2
           else
             call sparse_product( selfK1 , O_bar           , invOk  ) ! K1= O_bar x invOk
#ifdef debug
             write(*,*) 'SHOWING O_bar x O_k'
             call sparse_show_matrix(selfK1,show_elems=.true.) 
#endif
             call sparse_product( selfK2 , w_self_energy_v , selfK1 ) ! K2= Sigma_imp x O_bar x invOk
             call sparse_product( selfK1 , invOk           , O_bar  ) ! K1= invOk x O_bar
             call sparse_product( selfK3 , selfK1          , selfK2 ) ! K3= K1 x K2
           endif
#endif

           6163 continue

           if(same_self_for_all)then
             if(mpi_size>1.and.use_gpu_some_only.and..not.breaknfscom)then
                if(mpi_rank==1)then; call sparse_convert(mat_square,selfK3); endif;
                call sync_nodes_by_shared_nfs_matc(mat_square,'upfold_self_imp'//trim(adjustl(tostring(icall))),mpi_rank,mpi_size,nnn)
                call sparse_convert(selfK3,mat_square); deallocate(mat_square)
             endif
           endif

#ifndef GPU_SPEEDUP_MATMUL
          call sparse_product(self_energy_vK,hub_overlap_matrix_cplx,selfK3)
          call sparse_product(selfK,self_energy_vK,hub_overlap_tmatrix_cplx)
#else
          if(use_gpu_onlyme)then
           call matmul_in_dense_format(self_energy_vK,hub_overlap_matrix_cplx,selfK3)
           call matmul_in_dense_format(selfK,self_energy_vK,hub_overlap_tmatrix_cplx)
          else
           call sparse_product(self_energy_vK,hub_overlap_matrix_cplx,selfK3)
           call sparse_product(selfK,self_energy_vK,hub_overlap_tmatrix_cplx)
          endif
#endif

#ifdef debug
           write(*,*) 'WARNING : checking : selfK'
           call check_sparsec_(selfK)
           write(*,*) 'WARNING : checking : self_energy'
           call check_sparsec_(self_energy)
#endif

           return

         endif     


! GAMMA ONLY VERSION
#ifndef GPU_SPEEDUP_MATMUL
          call sparse_product(self_energy_v,hub_overlap_matrix_cplx,w_self_energy_v)
#else
          if(use_gpu_onlyme)then
           call matmul_in_dense_format(self_energy_v,hub_overlap_matrix_cplx,w_self_energy_v)
          else
                 call sparse_product(self_energy_v,hub_overlap_matrix_cplx,w_self_energy_v)
          endif
#endif

          if(verboseall)write(*,*) 'running transfo of Sigma step2'
#ifndef GPU_SPEEDUP_MATMUL
          call sparse_product(self_energy,self_energy_v,hub_overlap_tmatrix_cplx)
#else
          if(use_gpu_onlyme)then
           call matmul_in_dense_format(self_energy,self_energy_v,hub_overlap_tmatrix_cplx)
          else
                 call sparse_product(self_energy,self_energy_v,hub_overlap_tmatrix_cplx)
          endif
#endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine construct_self_energy

!======PROJECTED SELF=====!

!PART 1
            i=pub_num_hub_atoms_on_node(pub_my_node_id)

            if(.not.allocated(site_self_energy_buffer)) allocate(site_self_energy_buffer(i,i,12,12),stat=ierr)
            if(size(site_self_energy_buffer,1)/=i)then
             deallocate(site_self_energy_buffer); allocate(site_self_energy_buffer(i,i,12,12),stat=ierr)
            endif

            if(.not.silent.and.pub_my_node_id==0) write(*,*) '===== BUILD SELF ENERGY FROM EXTENDED ATOMS ====='
            site_self_energy_buffer=0.0_DP
            call build_self_energy_from_extended_atoms
            if(.not.silent.and.pub_my_node_id==0) write(*,*) 'node done with building self energy : ', pub_my_node_id
            call comms_barrier

!PART 2
           !site_self_energy_buffer is defined from previous section
            if(verboseall.and.pub_my_node_id==0) write(*,*) '====== DISTRIBUTE SELF ENERGY TO LOCAL ATOMS ======'
            call distribute_self_energy_to_atomic_components

            deallocate(site_self_energy_buffer,stat=ierr)
            call comms_barrier

           if(verboseall) then
            write(*,*) '===== SHOW THE SELF ENERGY MATRIX IF NECESSARY ====='
            call show_the_atomic_self
           endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine renorm_Z(i)
 implicit none
 integer :: k,l,m,i
       do k=1,size(Z_density_fine,1)
        do l=1,size(Z_density_fine,2)
         do m=1,size(Z_density_fine,3)
           Z_density_fine(k,l,m,i)=1.d0/(1.d0-Z_density_fine(k,l,m,i))
           if(Z_density_fine(k,l,m,i)>1.0) Z_density_fine(k,l,m,i)=0.d0
           if(Z_density_fine(k,l,m,i)<0.0) Z_density_fine(k,l,m,i)=0.d0
         enddo
        enddo
       enddo
 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine one_minus_Z(i)
 implicit none
 integer :: k,l,m,i
       do k=1,size(Z_density_fine,1)
        do l=1,size(Z_density_fine,2)
         do m=1,size(Z_density_fine,3)
           Z_density_fine(k,l,m,i)=1.d0-Z_density_fine(k,l,m,i)
         enddo
        enddo
       enddo
 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 logical function dimer_connected(hat,hatb)
 implicit none
 integer :: hat,hatb
  dimer_connected=dimer_table(hat,2)==hatb
  if(hubard_on_my_node(hat)==0.or.hubard_on_my_node(hatb)==0)then
    write(*,*) 'connected atoms are not on my node in dimer_connected test'
    stop
  endif   
 end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 integer function hubard_on_my_node(j)
 implicit none
 integer :: j,i,k,hub_atom,hat
 hubard_on_my_node=0
 do hub_atom = 1, pub_num_hub_atoms_on_node(pub_my_node_id)
   hat = pub_hub_atoms_on_node(hub_atom,pub_my_node_id);
   if(hat==j)then
    hubard_on_my_node=hub_atom
    return
   endif
 enddo
 end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plotting_real_space_quantities(build_or_plot)

 use cell_grid,    only: pub_fine_grid
 use constants,    only: DP
 use density,      only: density_on_grid
 use visual,       only: visual_scalarfield

 implicit none

 integer     :: u(1),i,ii,ierr,build_or_plot,k,l,m
 type(SPAM3) :: tmp_spam_complex

   if(.not. pub_dmft_plot_real_space) return 
   if(      nkpoints>1              ) then
      write(*,*) 'ERROR skipping real_space_quantities plots due to k points'
      write(*,*) ' if you want to plot real space quantities, you will have to average over k'
      return
   endif

   if(split)then
     write(*,*) 'ERROR plotting real space quantities but this is not compatible with split'
     stop
   endif

   write(*,*) 'in routine plotting_real_space_quantities, build matrix or plot ? ', build_or_plot

  !do not use energyzero here...check below...!
   u  = minloc (  abs (  (/( en_start+real(i-1,kind=DP)*en_step,i=1,pub_dmft_points )/) - pub_dmft_e_plot ) )
   ii = u(1)

   if(build_or_plot==1)then

   if(ien==1.and.pub_dmft_integrate_green)then
      call sparse_scale(green_tot_spam(is),0.0d0)
   endif
  
   if(  ((ien==1.or.ien==2.or.(pub_dmft_integrate_green.and.ien<=pub_dmft_points-nspecial_frequ)).and.pub_dmft_temp>0.0_DP) .or. &
      &  (pub_dmft_temp<0.0_DP.and.(ien==ii.or.ien==ii+1.or.(pub_dmft_integrate_green.and.ien<=ii.and.ien>1)) )  ) then

     real_space_exit = ien==2                                 .and. .not.pub_dmft_integrate_green .and. pub_dmft_temp>0.0_DP
     real_space_exit = real_space_exit .or. &
                    & (ien==pub_dmft_points-nspecial_frequ                                        .and. pub_dmft_temp>0.0_DP)
     real_space_exit = real_space_exit .or. &
                    & (ien==ii+1                              .and.                                     pub_dmft_temp<0.0_DP) 

     if(.not.allocated(sigma_density_fine))then
          allocate(sigma_density_fine(pub_fine_grid%ld1,pub_fine_grid%ld2, &
          pub_fine_grid%max_slabs12,pub_cell%num_spins),stat=ierr)
     endif
     if(.not.allocated(green_density_fine))then
          allocate(green_density_fine(pub_fine_grid%ld1,pub_fine_grid%ld2, &
          pub_fine_grid%max_slabs12,pub_cell%num_spins),stat=ierr)
     endif
     if(.not.allocated(green_density_fine_proj))then
          allocate(green_density_fine_proj(pub_fine_grid%ld1,pub_fine_grid%ld2, &
          pub_fine_grid%max_slabs12,pub_cell%num_spins),stat=ierr)
     endif
     if(.not.allocated(Z_density_fine))then
          allocate(Z_density_fine(pub_fine_grid%ld1,pub_fine_grid%ld2, &
          pub_fine_grid%max_slabs12,pub_cell%num_spins),stat=ierr)
     endif
     if(.not.allocated(green_tot_fine))then
          allocate(green_tot_fine(pub_fine_grid%ld1,pub_fine_grid%ld2, &
          pub_fine_grid%max_slabs12,pub_cell%num_spins),stat=ierr)
     endif


     if(pub_dmft_integrate_green)then
       call sparse_create(tmp_spam_complex,greenf,iscmplx=.true.)
       call sparse_copy(tmp_spam_complex,greenf)
       if(pub_dmft_temp>0.d0)then
        if(ien<pub_dmft_points-nspecial_frequ)then
          call sparse_axpy(green_tot_spam(is),tmp_spam_complex,2.0*abs(pub_dmft_temp)) 
        endif
        if(ien==pub_dmft_points-nspecial_frequ)then
          call sparse_axpy(green_tot_spam(is),tmp_spam_complex,2.0*abs(pub_dmft_temp)) 
          call sparse_scale(tmp_spam_complex,cmplx(0.0,0.5*(en_start+real(ien-1,kind=DP)*en_step),kind=DP))
          call sparse_axpy(green_tot_spam(is),tmp_spam_complex,1.0d0)
          call sparse_copy(tmp_spam_complex,green_tot_spam(is))
          write(*,*) 'OCCUPATION SPIN : ', is
          call matmul_in_dense_format(tmp_spam_complex,tmp_spam_complex,overlap_cplx,show_trace=.true.)
          if(is==2)then
            call sparse_copy(tmp_spam_complex,green_tot_spam(1))
            call sparse_axpy(tmp_spam_complex,green_tot_spam(2),-1.0d0) 
            call sparse_scale(tmp_spam_complex,cmplx(0.5,0.0,kind=DP))
            write(*,*) 'INTEGRATED MAGNETISM'
            call matmul_in_dense_format(tmp_spam_complex,tmp_spam_complex,overlap_cplx,show_trace=.true.)
          endif
        endif
       endif
       if(pub_dmft_temp<0.d0)then
        if(ien<=ii.and.ien>1)then
          call sparse_scale(tmp_spam_complex,cmplx(0.0,en_step/dacos(-1.d0),kind=DP))
          call sparse_axpy(green_tot_spam(is),tmp_spam_complex,1.0d0) 
        endif
        if(ien==ii)then
          call sparse_copy(tmp_spam_complex,green_tot_spam(is))
          write(*,*) 'OCCUPATION SPIN : ', is
          call matmul_in_dense_format(tmp_spam_complex,tmp_spam_complex,overlap_cplx,show_trace=.true.)
          if(is==2)then
            call sparse_copy(tmp_spam_complex,green_tot_spam(1))
            call sparse_axpy(tmp_spam_complex,green_tot_spam(2),-1.0d0) 
            call sparse_scale(tmp_spam_complex,cmplx(0.5,0.0,kind=DP))
            write(*,*) 'INTEGRATED MAGNETISM'
            call matmul_in_dense_format(tmp_spam_complex,tmp_spam_complex,overlap_cplx,show_trace=.true.)
          endif
        endif
       endif
       call sparse_destroy(tmp_spam_complex)
     endif

     if( (ien==1.and.pub_dmft_temp>0.d0) .or. (ien==ii.and.pub_dmft_temp<0.d0) )then
       !ImGreen
        call sparse_create(tmp_spam_complex,greenf,iscmplx=.true.)
        call sparse_copy(tmp_spam_complex,greenf)
        call sparse_scale(tmp_spam_complex,cmplx(0.0,1.0/dacos(-1.d0),kind=DP))
        call sparse_copy(green_density_spam(is),tmp_spam_complex)
        call sparse_destroy(tmp_spam_complex)
        if(inv_projectors)then
          !ImGreen_proj
           call sparse_create(tmp_spam_complex,greenf_v,hub_inv_overlap_tmatrix_cplx,iscmplx=.true.)
           call sparse_product(greenf_v,hub_inv_overlap_matrix_cplx,w_greenf_v)
           call sparse_product(tmp_spam_complex,greenf_v,hub_inv_overlap_tmatrix_cplx)
           call sparse_scale(tmp_spam_complex,cmplx(0.0,1.0/dacos(-1.d0),kind=DP))
           call sparse_copy(green_density_spam_proj(is),tmp_spam_complex)
           call sparse_destroy(tmp_spam_complex)
          !ImSelf
           call sparse_create(tmp_spam_complex,self_energy_v,hub_inv_overlap_tmatrix_cplx,iscmplx=.true.)
           call sparse_product(self_energy_v,hub_inv_overlap_matrix_cplx,w_self_energy_v)
           call sparse_product(tmp_spam_complex,self_energy_v,hub_inv_overlap_tmatrix_cplx)
           if(pub_dmft_temp>0.d0)then
             call sparse_scale(tmp_spam_complex,cmplx(0.0,-1.0,kind=DP))
           endif
           call sparse_copy(sigma_density_spam(is),tmp_spam_complex)
           call sparse_destroy(tmp_spam_complex)
        endif
     endif

     if( (ien==2.and.pub_dmft_temp > 0.d0) .or. (ien==ii+1.and.pub_dmft_temp<0.d0) )then
      if(inv_projectors)then
        call sparse_create(tmp_spam_complex,self_energy_v,hub_inv_overlap_tmatrix_cplx,iscmplx=.true.)
        call sparse_product(self_energy_v,hub_inv_overlap_matrix_cplx,w_self_energy_v)
        call sparse_product(tmp_spam_complex,self_energy_v,hub_inv_overlap_tmatrix_cplx)
        if(pub_dmft_temp>0.d0)then
          call sparse_scale(tmp_spam_complex,cmplx(0.0,-1.0,kind=DP))
        endif
        call sparse_copy(Z_density_spam(is),tmp_spam_complex)
        call sparse_destroy(tmp_spam_complex)
      endif
     endif

    endif

   else

     if(allocated(sigma_density_fine).and.allocated(green_density_fine).and. &
      & allocated(Z_density_fine)    .and.allocated(green_density_fine_proj) .and. allocated(green_tot_fine) )then
        write(*,*) 'building real-space grid for density plots' 
        if(pub_dmft_paramagnetic.and.pub_cell%num_spins/=1)then
          call sparse_copy( green_density_spam_proj(2) , green_density_spam_proj(1) )
          call sparse_copy( green_density_spam(2)      , green_density_spam(1)      )
          call sparse_copy( sigma_density_spam(2)      , sigma_density_spam(1)      )
          call sparse_copy(     Z_density_spam(2)      ,     Z_density_spam(1)      )
          call sparse_copy(     green_tot_spam(2)      ,     green_tot_spam(1)      )
        endif
        write(*,*) '-------------------------------build density-------------------------------------'
        write(*,*) 'go for density plot, green_density'
        if(.not.pub_dmft_plot_real_space_sigma) call density_on_grid(green_density_fine,pub_fine_grid,green_density_spam,overlap,rep%ngwfs_on_grid,ngwf_basis)
        write(*,*) 'go for density plot, green_density_proj'
        if(.not.pub_dmft_plot_real_space_sigma) call density_on_grid(green_density_fine_proj,pub_fine_grid,green_density_spam_proj,overlap,rep%ngwfs_on_grid,ngwf_basis)
        write(*,*) 'go for density plot, green_tot'
        if(.not.pub_dmft_plot_real_space_sigma.and.pub_dmft_integrate_green) &
                                              & call density_on_grid(green_tot_fine    ,pub_fine_grid,green_tot_spam    ,overlap,rep%ngwfs_on_grid,ngwf_basis)
        write(*,*) 'go for density plot, sigma_density'
                                                call density_on_grid(sigma_density_fine,pub_fine_grid,sigma_density_spam,overlap,rep%ngwfs_on_grid,ngwf_basis)
        write(*,*) 'go for density plot, Z_density'
        if(.not.pub_dmft_plot_real_space_sigma) call density_on_grid(Z_density_fine    ,pub_fine_grid,    Z_density_spam,overlap,rep%ngwfs_on_grid,ngwf_basis)
        write(*,*) '--------------------------done, now go for plots--------------------------------'
        if(pub_my_node_id==0)then
         open(unit=2121,file='density_cut_G'//trim(adjustl(toString(1))))
         write(2121,*) maxval(green_density_fine(:,:,1,1)),minval(green_density_fine(:,:,1,1))
         do i=1,size(green_density_fine,1)
          do j=1,size(green_density_fine,2)
           write(2121,'(i5,i5,f10.5)') i,j,green_density_fine(i,j,1,1)
          enddo
         enddo
         close(2121)
         open(unit=2121,file='density_cut_Sigma'//trim(adjustl(toString(1))))
         write(2121,*) maxval(sigma_density_fine(:,:,1,1)),minval(sigma_density_fine(:,:,1,1))
         do i=1,size(sigma_density_fine,1)
          do j=1,size(sigma_density_fine,2)
           write(2121,'(i5,i5,f10.5)') i,j,sigma_density_fine(i,j,1,1)
          enddo
         enddo
         close(2121)
         open(unit=2121,file='density_cut_Z'//trim(adjustl(toString(1))))
         write(2121,*) maxval(Z_density_fine(:,:,1,1)),minval(Z_density_fine(:,:,1,1))
         do i=1,size(Z_density_fine,1)
          do j=1,size(Z_density_fine,2)
           write(2121,'(i5,i5,f10.5)') i,j,Z_density_fine(i,j,1,1)
          enddo
         enddo
         close(2121)
        endif
        do i=1,pub_cell%num_spins
         if(.not.pub_dmft_plot_real_space_sigma)then
                           !--------------------------------------!
           write(*,*) 'display scalarfield in cube format, go for green_proj'
           call visual_scalarfield(green_density_fine_proj(:,:,:,i), pub_fine_grid,  &
          'density_proj_fermi_level_dmft', 'density_proj_fermi_level_dmft'//trim(adjustl(toString(i))), elements, 1.0_DP)
           write(*,*) 'display scalarfield in cube format, go for green'
           call visual_scalarfield(green_density_fine(:,:,:,i), pub_fine_grid,  &
          'density_fermi_level_dmft_for_spin', 'density_fermi_level_dmft_for_spin'//trim(adjustl(toString(i))), elements, 1.0_DP)
           write(*,*) 'display scalarfield in cube format, go for green tot'
           if(pub_dmft_integrate_green) &
       &   call visual_scalarfield(green_tot_fine(:,:,:,i), pub_fine_grid,  &
       &  'density_dmft_for_spin', 'density_dmft_for_spin'//trim(adjustl(toString(i))), elements, 1.0_DP)
         endif
           write(*,*) 'display scalarfield in cube format, go for sigma'
           call visual_scalarfield(sigma_density_fine(:,:,:,i), pub_fine_grid,  &
          'scattering_fermi_level_dmft', 'scattering_fermi_level_dmft'//trim(adjustl(toString(i))), elements,HARTREE_IN_EVS)
         if(.not.pub_dmft_plot_real_space_sigma)then
           call visual_scalarfield(abs(sigma_density_fine(:,:,:,i)), pub_fine_grid,  &
          'abs_scattering_fermi_level_dmft', 'abs_scattering_fermi_level_dmft'//trim(adjustl(toString(i))), elements,HARTREE_IN_EVS)
                            !--------------------------------------!
           Z_density_fine(:,:,:,i)=(Z_density_fine(:,:,:,i)-sigma_density_fine(:,:,:,i))/en_step
                            !--------------------------------------!
           write(*,*) 'display scalarfield in cube format, go for slope'
           call visual_scalarfield(Z_density_fine(:,:,:,i), pub_fine_grid,  &
           'Slope_fermi_level_dmft', 'Slope_fermi_level_dmft'//trim(adjustl(toString(i))),elements, 1.d0)
           write(*,*) 'display scalarfield in cube format, go for abs slope'
           call visual_scalarfield(abs(Z_density_fine(:,:,:,i)), pub_fine_grid,  &
           'Abs_Slope_fermi_level_dmft', 'Abs_Slope_fermi_level_dmft'//trim(adjustl(toString(i))),elements, 1.d0)
                            !--------------------------------------!
           call renorm_Z(i)
                            !--------------------------------------!
           call visual_scalarfield(Z_density_fine(:,:,:,i), pub_fine_grid,  &
           'Z_fermi_level_dmft', 'Z_fermi_level_dmft'//trim(adjustl(toString(i))),elements, 1.d0)
                           !--------------------------------------!
           if(pub_dmft_temp>0.d0)then
             Z_density_fine(:,:,:,i)=sigma_density_fine(:,:,:,i)/en_start
             call visual_scalarfield(Z_density_fine(:,:,:,i), pub_fine_grid,  &
             'Slope_first_only_fermi_level_dmft_for_spin', 'Slope_first_only_fermi_level_dmft_for_spin'//trim(adjustl(toString(i))),elements, 1.d0)
             call renorm_Z(i)
             call visual_scalarfield(Z_density_fine(:,:,:,i), pub_fine_grid,  &
             'Z_first_only_fermi_level_dmft_for_spin', 'Z_first_only_fermi_level_dmft_for_spin'//trim(adjustl(toString(i))),elements, 1.d0)
             call one_minus_Z(i)
             call visual_scalarfield(Z_density_fine(:,:,:,i), pub_fine_grid,  &
             '1mZ_first_only_fermi_level_dmft_for_spin', '1mZ_first_only_fermi_level_dmft_for_spin'//trim(adjustl(toString(i))),elements, 1.d0)
           endif
          endif
                           !--------------------------------------!

        enddo

        if(.not.pub_dmft_plot_real_space_sigma.and..not.pub_dmft_paramagnetic.and.pub_cell%num_spins/=1)then
           write(*,*) 'display scalarfield in cube format, go for magnetic density at Fermi level'
           call visual_scalarfield((green_density_fine_proj(:,:,:,2)-green_density_fine_proj(:,:,:,1))/2.d0, pub_fine_grid,  &
          'mag_density_proj_fermi_level_dmft', 'mag_density_proj_fermi_level_dmft', elements, 1.0_DP)
           write(*,*) 'display scalarfield in cube format, go for proj magnetic density at Fermi level'
           call visual_scalarfield((green_density_fine(:,:,:,2)-green_density_fine(:,:,:,1))/2.d0, pub_fine_grid,  &
          'mag_density_fermi_level_dmft_for_spin', 'mag_density_fermi_level_dmft_for_spin', elements, 1.0_DP)
           write(*,*) 'display scalarfield in cube format, go for magnetic density'
           if(pub_dmft_integrate_green) &
      &    call visual_scalarfield((green_tot_fine(:,:,:,2)-green_tot_fine(:,:,:,1))/2.d0, pub_fine_grid,  &
      &   'mag_density_dmft_for_spin', 'mag_density_dmft_for_spin', elements, 1.0_DP)
        endif

     endif

   endif

return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_projected_overlap_matrix_Sproj
implicit none
integer                  :: i,k1,k2,j,l1,l2,hub_atom
complex(DP), allocatable :: matout(:,:)

          if(.not.silent.and.pub_my_node_id==0) write(*,*) 'START TO COMPUTE PROJECTED OVERLAP MATRIX' 
          if(.not.inv_projectors)Then
           write(*,*) 'ERROR you need inv_projectors=.true. to compute the projected overlap matrix'
           return
          endif
          if(silent) return

          if(pub_total_num_nodes/=1.or.is/=1.or.nmu_step/=1.or.ikpoint/=1) return

          energy=0.d0; ien=1; kien=1;  ! FOR MATMUL ON GPU BELOW - ENFORCE DOUBLE PREC !
 
#ifndef GPU_SPEEDUP_MATMUL
          if(pub_my_node_id==0) write(*,*) 'Compute < \psi^m | S | phi^\alpha >'
          call sparse_product(w_inv_greenf,hub_inv_overlap_tmatrix_cplx,overlap_cplx)
#else
          if(use_gpu_onlyme)then
          if(pub_my_node_id==0) write(*,*) 'Compute < \psi^m | S | phi^\alpha >'
          call matmul_in_dense_format(w_inv_greenf,hub_inv_overlap_tmatrix_cplx,overlap_cplx)
          else
                call sparse_product(w_inv_greenf,hub_inv_overlap_tmatrix_cplx,overlap_cplx)
          endif
#endif
          
          if(pub_my_node_id==0) write(*,*) 'Compute < \psi^m | S | psi_m_prime >'
#ifndef GPU_SPEEDUP_MATMUL
          call sparse_product(w_S_v,w_inv_greenf,hub_inv_overlap_matrix_cplx)
#else
          if(use_gpu_onlyme)then
          call matmul_in_dense_format(w_S_v,w_inv_greenf,hub_inv_overlap_matrix_cplx)
          else
                  call sparse_product(w_S_v,w_inv_greenf,hub_inv_overlap_matrix_cplx)
          endif
#endif
          if(pub_my_node_id==0) write(*,*) 'show overlap matrix'

          do hub_atom = 1, pub_num_hub_atoms_on_node(pub_my_node_id)
           i =pub_hub_atoms_on_node(hub_atom,pub_my_node_id);
           k1=hub_atom_hub_projs(i)
           k2=hub_atom_hub_projs(i)
           allocate(matout(k1,k2)) 
           if(i==1)then
            call show_hubbard_atom(w_S_v,i,matout,showdiag=i==1,fullmat=overlap_hubbard(:,:,is))
           else
            call show_hubbard_atom(w_S_v,i,matout,showdiag=i==1)
           endif
#ifdef debug
            if(maxval(abs(matout))<1.d-8)then
              write(*,*) 'ERROR - HUBBARD PROJECTION MATRIX IS NAUGHT !!! '
              call write_array(real(matout), ' U^+ S U  ' )
              stop
            endif
#endif 
           if(pub_my_node_id==0)then
            write(*,*) '======================================================='
            write(*,*) 'Re(S) matrix of ATOM [x] / total [y] : ' , i, pub_hub_nat 
            do l1=1,k1
             write(*,'(200f10.3)') (real(matout(l1,l2)),l2=1,k2) 
            enddo
           endif
           deallocate(matout)
          enddo

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine project_green_function(w_greenf_v,greenf,nosplit)
implicit none
logical,optional :: nosplit
TYPE(SPAM3)      :: w_greenf_v,greenf 
#ifdef debug
TYPE(SPAM3)      :: test
#endif

          if(.not.present(nosplit))then
           if(split)then
            if(pub_my_node_id==0) write(*,*) 'BUILD PROJECTED GF'
#ifdef debug
             if(size(split_hubt,1)/=size(split_wgv,1).or.size(split_hub,2)/=size(split_wgv,2))then
               write(*,*) 'ERROR sizes dont match split hubt';stop
             endif
#endif
             split_wgv = MATMUL( MATMUL( split_hubt,greenfbackup),split_hub ); return
           endif
          endif

          if(verboseall)write(*,*) 'Compute < \psi^m | \hat{G} | phi^\alpha >'
#ifndef GPU_SPEEDUP_MATMUL
                  call sparse_product(w_greenf,hub_overlap_tmatrix_cplx,greenf)
#else
          if(use_gpu_onlyme)then
          if(verboseall)write(*,*) 'go for full GPU'
          call matmul_in_dense_format(w_greenf,hub_overlap_tmatrix_cplx,greenf)
          else
                call sparse_product(w_greenf,hub_overlap_tmatrix_cplx,greenf)
#ifdef debug
                write(*,*) 'CHECKING SPARSE PRODUCT'
                call sparse_create(test,w_greenf,iscmplx=.true.);
                call matmul_in_dense_format(test,hub_overlap_tmatrix_cplx,greenf) 
                call sparse_max_diff(test,w_greenf)
                call sparse_destroy(test)
#endif 
          endif
#endif

          if(verboseall)write(*,*) 'Compute < \psi^m | \hat{G} | psi_m_prime >'

#ifndef GPU_SPEEDUP_MATMUL
                  call sparse_product(w_greenf_v,w_greenf,hub_overlap_matrix_cplx)
#else
          if(use_gpu_onlyme)then
                if(verboseall)write(*,*) 'go for full GPU'
                call matmul_in_dense_format(w_greenf_v,w_greenf,hub_overlap_matrix_cplx)
          else
                call sparse_product(w_greenf_v,w_greenf,hub_overlap_matrix_cplx)
#ifdef debug
                write(*,*) 'CHECKING SPARSE PRODUCT'
                call sparse_create(test,w_greenf_v,iscmplx=.true.);
                call matmul_in_dense_format(test,w_greenf,hub_overlap_matrix_cplx)
                call sparse_max_diff(test,w_greenf_v)
                call sparse_destroy(test)
#endif
          endif
#endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine project_hamiltonian
implicit none

          if(.not.inv_projectors)then
           write(*,*) 'ERROR you need to set inv_projectors=.true. to compute the projected hamiltonian'
           return
          endif

          if(verboseall)write(*,*) '<<<<<< START HAMILTONIAN PROJECTION >>>>>>>'
#ifndef GPU_SPEEDUP_MATMUL
          if(verboseall)write(*,*) 'Compute < \psi^m | E-\hat{H} | phi^\alpha >'
          call sparse_product(w_inv_greenf,hub_inv_overlap_tmatrix_cplx,inv_greenf)
#else
          if(use_gpu_onlyme)then
          if(verboseall)write(*,*) 'Compute < \psi^m | E-\hat{H} | phi^\alpha >'
                call matmul_in_dense_format(w_inv_greenf,hub_inv_overlap_tmatrix_cplx,inv_greenf)
          else
                call sparse_product(w_inv_greenf,hub_inv_overlap_tmatrix_cplx,inv_greenf)
          endif
#endif

          if(verboseall)write(*,*) 'Compute < \psi^m | E-\hat{H} | psi_m_prime >'
#ifndef GPU_SPEEDUP_MATMUL
          call sparse_product(w_inv_greenf_v,w_inv_greenf,hub_inv_overlap_matrix_cplx)
#else
          if(use_gpu_onlyme)then
                call matmul_in_dense_format(w_inv_greenf_v,w_inv_greenf,hub_inv_overlap_matrix_cplx)
          else
                call sparse_product(w_inv_greenf_v,w_inv_greenf,hub_inv_overlap_matrix_cplx)
          endif
#endif
          if(verboseall)write(*,*) 'done'
          if(verboseall)write(*,*) '<<<<<< END HAMILTONIAN PROJECTION >>>>>>>'
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_sparse_and_dimers
implicit none
integer :: ii

    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [1]'
    greenf%structure = 'K'
    call sparse_create(greenf,iscmplx=.true.)
    call sparse_create(inv_greenf,ham(1),iscmplx=.true.)
    call sparse_create(inv_greenf_backup,ham(1),iscmplx=.true.)
    call sparse_create(greenf_backup,greenf,iscmplx=.true.)
    call sparse_create(overlap_cplx,overlap,iscmplx=.true.)
    call sparse_create(hub_overlap_matrix_cplx, hub_overlap_matrix, iscmplx=.true.)
    call sparse_create(hub_overlap_tmatrix_cplx,hub_overlap_tmatrix,iscmplx=.true.)
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [2]'
    call sparse_create(w_greenf,hub_overlap_tmatrix_cplx,greenf,iscmplx=.true.) ! WGr

    select case(cluster_scheme)
      case(1)
      if(.not.silent.and.pub_my_node_id==0) write(*,*) 'using cluster scheme 1'
      w_greenf_v%structure = 'G' ! block-diagonal 
      case(2)
      if(.not.silent.and.pub_my_node_id==0) write(*,*) 'using cluster scheme 2'
      w_greenf_v%structure = 'FULL'
      case(3)
      if(.not.silent.and.pub_my_node_id==0) write(*,*) 'using cluster scheme 3'
      call sparse_create(w_greenf_v,w_greenf,hub_overlap_matrix_cplx,iscmplx=.true.)
    end select
    if(cluster_scheme/=3) then
      if(.not.silent.and.pub_my_node_id==0) write(*,*) 'allocating w_greenf_v'
      call sparse_create(w_greenf_v,iscmplx=.true.) ! WGrV->'G'
    endif

    call build_my_dimers
 
    if(nkpoints>1)then
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'building arrays for K points'
       Ok%structure='FULL'
       call sparse_create(       Ok,iscmplx=.true.)
       call sparse_create( invOk,Ok,iscmplx=.true.) 
       call sparse_create(selfK1,Ok,iscmplx=.true.)
       call sparse_create(selfK2,Ok,iscmplx=.true.)
       call sparse_create(selfK3,Ok,iscmplx=.true.)
       call sparse_create( O_bar,Ok,iscmplx=.true.)
       call sparse_create(self_energy_vK,hub_overlap_matrix_cplx,selfK2,iscmplx=.true.)
       call sparse_create(selfK,greenf,iscmplx=.true.)
    endif

    if(inv_projectors)then
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'building arrays for inv. projectors'
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [3]'
       call sparse_create(hub_inv_overlap_matrix,inv_overlap,hub_overlap_matrix)
       call sparse_create(hub_inv_overlap_matrix_cplx,hub_inv_overlap_matrix,iscmplx=.true.)
       call sparse_create(hub_inv_overlap_tmatrix,hub_overlap_tmatrix,inv_overlap)
       call sparse_create(hub_inv_overlap_tmatrix_cplx,hub_inv_overlap_tmatrix,iscmplx=.true.)
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [4]'
       call sparse_create(w_inv_greenf,hub_inv_overlap_tmatrix_cplx,inv_greenf,iscmplx=.true.)
       select case(cluster_scheme)
        case(1)
          w_inv_greenf_v%structure = 'G' ! Blank to block-diagonal 'G' structure
        case(2)
          w_inv_greenf_v%structure = 'FULL'
        case(3)
          call sparse_create(w_inv_greenf_v,w_greenf,hub_overlap_matrix_cplx,iscmplx=.true.)
       end select
       if(.not.cluster_scheme==3) call sparse_create(w_inv_greenf_v,iscmplx=.true.)
    endif

    select case(cluster_scheme)
    case(1)
      w_self_energy_v%structure = 'G' ! Blank to block-diagonal 'G' structure
    case(2)
      w_self_energy_v%structure = 'FULL'
    case(3)
      call sparse_create(w_self_energy_v,w_greenf,hub_overlap_matrix_cplx,iscmplx=.true.)
    end select
    if(.not.cluster_scheme==3) call sparse_create(w_self_energy_v,iscmplx=.true.)

    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [5]'
    if(.not.silent.and.pub_my_node_id==0) write(*,*) '[5.1]'
    call sparse_create(greenf_v,hub_overlap_matrix_cplx,w_greenf_v,iscmplx=.true.)
    call sparse_create(self_energy_v,hub_overlap_matrix_cplx,w_self_energy_v,iscmplx=.true.)
    if(.not.silent.and.pub_my_node_id==0) write(*,*) '[5.2]'
    call sparse_create(self_energy,self_energy_v,hub_overlap_tmatrix_cplx,iscmplx=.true.)

    do i=1,pub_cell%num_spins
      if(nkpoints==1)then
       call sparse_create(self_infinity(i),ham(i))
      else
       call sparse_create(self_infinity(i),ham(i),iscmplx=.true.)
      endif
    enddo

    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [6]'

    if(.not.silent.and.pub_my_node_id==0) write(*,*)' -------------- PRELIMINARIES COMPLETE ----------------- '
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [7]'

    call sparse_copy(overlap_cplx,overlap)
    call sparse_copy(hub_overlap_matrix_cplx,hub_overlap_matrix)
    call sparse_copy(hub_overlap_tmatrix_cplx,hub_overlap_tmatrix)



    if(inv_projectors)then
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [8]'
#ifndef GPU_SPEEDUP_MATMUL_R
       call sparse_product(hub_inv_overlap_matrix,inv_overlap,hub_overlap_matrix)
#else
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [8.1]'
       if(use_gpu_onlyme)then
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [8.2]'
       call matmul_in_dense_format_r(hub_inv_overlap_matrix,inv_overlap,hub_overlap_matrix)
       else
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [8.3]'
               call sparse_product(hub_inv_overlap_matrix,inv_overlap,hub_overlap_matrix)
       endif
#endif
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [8.4]'
       call sparse_copy(hub_inv_overlap_matrix_cplx,hub_inv_overlap_matrix)
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [9]'

#ifndef GPU_SPEEDUP_MATMUL_R
       call sparse_product(hub_inv_overlap_tmatrix,hub_overlap_tmatrix,inv_overlap)
#else
       if(use_gpu_onlyme)then
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [9.1]'
       call matmul_in_dense_format_r(hub_inv_overlap_tmatrix,hub_overlap_tmatrix,inv_overlap)
       else
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [9.2]'
               call sparse_product(hub_inv_overlap_tmatrix,hub_overlap_tmatrix,inv_overlap)
       endif
#endif

       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [10]'

       call sparse_copy(hub_inv_overlap_tmatrix_cplx,hub_inv_overlap_tmatrix)
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'initialization sparse step [11]'
    endif



    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'allocating real-space density SPAM matrix'
    do i=1,pub_cell%num_spins
      call sparse_create( green_density_spam_proj(i), greenf                ,iscmplx=.false.)
      call sparse_create( green_density_spam(i)     , greenf                ,iscmplx=.false.)
      call sparse_create(     green_tot_spam(i)     , greenf                ,iscmplx=.false.)
      call sparse_create( sigma_density_spam(i)     , self_energy           ,iscmplx=.false.)
      call sparse_create(     Z_density_spam(i)     , sigma_density_spam(i) ,iscmplx=.false.)
    enddo

    call sparse_create(w_S_v,w_greenf_v,iscmplx=.true.)

    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'init sparse done'

    ii=ngwf_basis%num

    allocate(vectorAA(pub_cell%num_spins,ii))
#ifdef debug
    allocate(vectorAAback(pub_cell%num_spins,ii))
#endif
    if(nkpoints==1)then
      allocate(matrixAA(pub_cell%num_spins,ii,ii),matrixContributionTail(ii,ii),matrixOverlap(ii,ii),matrixH(pub_cell%num_spins,ii,ii))
    else
      allocate(matrixAAc(pub_cell%num_spins,ii,ii),matrixContributionTailc(ii,ii),matrixOverlapc(ii,ii),matrixHc(pub_cell%num_spins,ii,ii))
    endif

    if(nkpoints==1)then
     allocate(matrixDensKernel(ii,ii)) 
    else
     allocate(matrixDensKernelc(ii,ii))
    endif

    allocate(matrixAAt(pub_cell%num_spins,ii,ii),tmpc(ii,ii),matrixSigStat(pub_cell%num_spins,ii,ii), &
      & tmpr(ii,ii),tmpr2(ii,ii),tmp1(ii,ii),tmp2(ii,ii),greenfbackup(ii,ii), &
      & matrixAsig(pub_cell%num_spins,ii,ii),matrixGan(ii,ii),sigmabackup(ii,ii))

   if(nkpoints>1)then
    allocate(tmpc2(ii,ii))
   endif
#ifdef debug
   if(nkpoints==1)then
    allocate(tmpc2(ii,ii),tmpc3(ii,ii))
   else
    allocate(tmpc3(ii,ii))
   endif
#endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine self_energy_header_write_screen
         if(silent) return
         if(ien==1)then
          write(*,*) 'done...'
          write(*,*) '--------------------------------'
          write(*,*) '----- EXTRACT SELF-ENERGY ------'
          write(*,*) 'THE ATOMIC SELF ENERGY IS IN THE CUBIC (or real) HARMONIC basis with : '
          if(channels==5)then
           write(*,*) ' ELEMENT 1 IS m=-2 :dxy    '
           write(*,*) ' ELEMENT 2 IS m=-1 :dzy    '
           write(*,*) ' ELEMENT 3 IS m= 0 :d3z2-r2'
           write(*,*) ' ELEMENT 4 IS m= 1 :dxz    '
           write(*,*) ' ELEMENT 5 IS m= 2 :dx2-y2 '
          elseif(channels==7)then
           write(*,*) ' ELEMENT 1 IS m=-3,fy(3x^2-y^2) '
           write(*,*) ' ELEMENT 2 IS m=-2,fxyz         '
           write(*,*) ' ELEMENT 3 IS m=-1,fyz^2        '
           write(*,*) ' ELEMENT 4 IS m=0, fz^3         '
           write(*,*) ' ELEMENT 5 IS m=1, fxz^2        '
           write(*,*) ' ELEMENT 6 IS m=2, fz(x2-y2)    '
           write(*,*) ' ELEMENT 7 IS m=3, fx(x^2-3y^2) '
          endif
          write(*,*) '--------------------------------'
         endif
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_routine_variables
implicit none
logical :: ham_is_present
integer :: jij

      same_self_for_all=.false.
      force_reopen_file=.false.

      pub_hub_nproj = sum( (/( hub_atom_hub_projs(jij),jij=1,pub_hub_nat )/) )

      if(allocated(overlap_hubbard)) deallocate(overlap_hubbard)
      allocate(overlap_hubbard(pub_hub_nproj,pub_hub_nproj,2))

                  cluster_scheme=1
      if(cluster) cluster_scheme=2

      show_matrices=.false.

      mpi_size=1;mpi_rank=1
      num=COMMAND_ARGUMENT_COUNT()
      if(.not.silent.and.pub_my_node_id==0) write(*,*) 'NUMBER OF ARGUMENTS IN ONETEP : ', num
      do i=1,num
        call GET_COMMAND_ARGUMENT(i, value(i), len__, status__)
        if(.not.silent.and.pub_my_node_id==0) write(*,*) 'ARGUMENT : ',i,TRIM(ADJUSTL(value(i)))
      enddo
      if(num>=2)then
       mpi_size=StrInt2(value(3))
       mpi_rank=StrInt2(value(2))
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'RUNNING ONETEP AS MULTI_NODES : '
       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'MY RANK / SIZE                : ', mpi_rank,mpi_size
      endif
      if(num>=4)then
       pub_dmft_points = StrInt2(value(4))
       sign_temp       = StrInt2(value(5))
       pub_dmft_emin   = StringToReal(value(6))
       pub_dmft_emax   = StringToReal(value(7))
       if(sign_temp/=0) pub_dmft_temp=-abs(pub_dmft_temp)
       write(*,*) 'GOTCHA, DMFT POINTS IS        : ', pub_dmft_points
       write(*,*) '        DMFT TEMP IS          : ', pub_dmft_temp
       write(*,*) '        DMFT EMIN             : ', pub_dmft_emin
       write(*,*) '        DMFT EMAX             : ', pub_dmft_emax
      endif
      
      restrict_window = ((pub_dmft_dos_min>pub_dmft_emin).or.(pub_dmft_dos_max<pub_dmft_emax)).and.(pub_dmft_temp<0.0_DP)
      reopen_sigma_file_each_step =  mpi_size>1 .or. restrict_window .or. split

      use_gpu=(pub_total_num_nodes<=gpu_max).and.(mpi_size<=gpu_max)
      if(mpi_size==1.and.pub_total_num_nodes>1.and.gpu_max==1.and..not.pub_dmft_splitk)then
         write(*,*) 'WARNING switching off split due to the presence of a GPU'
         use_gpu=.true.;split=.false.;
      endif
      use_gpu_onlyme=use_gpu
      if(mpi_size/=1.and.pub_total_num_nodes==1.and.gpu_max>=1.and..not.pub_dmft_splitk)then
        use_gpu_onlyme=mpi_rank<=gpu_max 
      endif
      use_gpu_some_only=mpi_size/=1.and.pub_total_num_nodes==1.and.gpu_max>=1.and.gpu_max<mpi_size.and..not.pub_dmft_splitk

      if(pub_dmft_splitk) then
        if(use_gpu_onlyme) write(*,*) 'WARNING turning off GPU, because of splitk' 
        use_gpu=.false.; use_gpu_some_only=.false.; use_gpu_onlyme=.false. 
      endif 

      if(pub_my_node_id==0) write(*,*) 'THERE ARE [x] NODES RUNNING AS A BATCH : ', pub_total_num_nodes

      if(.not.silent.and.pub_my_node_id==0)write(*,*) 'MPI_RANK/NODE_IDE, USING GPU ? ', mpi_rank, pub_my_node_id, use_gpu
      if(.not.silent.and.pub_my_node_id==0)write(*,*) 'GPU MAX?                       ', gpu_max

      if(.not.silent.and.pub_my_node_id==0)write(*,*) 'am I using a GPU               ', use_gpu_onlyme
      if(.not.silent.and.pub_my_node_id==0)write(*,*) 'are some of us only using a GPU', use_gpu_some_only

#ifdef GPU_SPEEDUP
      if(use_gpu_onlyme)then
        call init_gpu_device
        if(pub_total_num_nodes==1)then
          call choose_gpu(mpi_rank-1)
        else
         if(pub_my_node_id<gpu_max)then
          call choose_gpu(pub_my_node_id)
         endif
        endif
        call cublas_init()
      endif
#endif

      do hat=1,pub_hub_nat
            if(allocated(h_atoms(hat)%occupancy)) deallocate(h_atoms(hat)%occupancy)
            channels = hub_atom_hub_projs(hat)
            allocate(h_atoms(hat)%occupancy(1:channels,1:channels,2))
            h_atoms(hat)%occupancy=0.d0
      enddo

      if(.not.silent.and.pub_my_node_id==0) write(*,*) 'build occupancy matrices'

      do hub_atom = 1, pub_num_hub_atoms_on_node(pub_my_node_id)
            hat      = pub_hub_atoms_on_node(hub_atom,pub_my_node_id)
            channels = hub_atom_hub_projs(hat)
            if(pub_node_of_atom(hub%orig(hat))/=pub_my_node_id)then
              write(*,*) 'ERROR in getting h_atom, node id do not match'
              write(*,*) ' hat                             : ', hat
              write(*,*) ' orig atom                       : ', hub%orig(hat)
              write(*,*) ' hub_atom                        : ', hub_atom
              write(*,*) ' pub_node_of_atom(hub%orig(hat)) : ', pub_node_of_atom(hub%orig(hat))
              write(*,*) ' pub node id                     : ', pub_my_node_id
              stop
            endif
            h_atoms(hat)%occupancy=0.d0
            h_atoms(hat)%occupancy(1:channels,1:channels,1)=h_atoms_occupancy(hub_proj_basis,hub,hat,1)
           if(num_spins==2)then
             h_atoms(hat)%occupancy(1:channels,1:channels,2)=h_atoms_occupancy(hub_proj_basis,hub,hat,2)
           endif
      enddo
      do hat=1,pub_hub_nat
            call mpisum(h_atoms(hat)%occupancy)
      enddo

       if(gbooth)then
         call restart_kernel_read(denskern_tmp, read_xc=.true.)
         call check_sparse_in_file(denskern_tmp(1),'store_xc1')
         if(num_spins==2)then
          INQUIRE(file='store_xc2',EXIST=ham_is_present)
          if(ham_is_present)then
            write(*,*) 'FILE store_xc2 is present'
            call check_sparse_in_file(denskern_tmp(2),'store_xc2')
          else
            write(*,*) 'FILE store_xc2 is not present, using store_xc1'
            call sparse_copy(denskern_tmp(2),denskern_tmp(1))
          endif
         endif

         call restart_kernel_read(denskern_tmp, read_hartree=.true.)
         call check_sparse_in_file(denskern_tmp(1),'store_hartree1')
         if(num_spins==2)then
          INQUIRE(file='store_hartree2',EXIST=ham_is_present)
          if(ham_is_present)then
            write(*,*) 'FILE store_hartree2 is present'
            call check_sparse_in_file(denskern_tmp(2),'store_hartree2')
          else
            write(*,*) 'FILE store_hartree2 is not present, using store_hartree1'
            call sparse_copy(denskern_tmp(2),denskern_tmp(1))
          endif
         endif
       endif

       if(num_spins==2)then
        INQUIRE(file='store_ham1',EXIST=ham_is_present)
        if(ham_is_present)then
          INQUIRE(file='store_ham2',EXIST=ham_is_present)
          if(.not.ham_is_present)then
            write(*,*) 'WARNING: FILE store_ham1 is present, but not store_ham2'
            write(*,*) 'copying ham_2 to ham_1'
            call check_sparse_in_file(ham(1),'store_ham1')
            call sparse_copy(ham(2),ham(1))
            call check_sparse_in_file(ham(2),'store_ham2')
          endif 
        endif
       endif

       call check_sparse_in_file(ham(1),'store_ham1')
      if(num_spins==2)then
       call check_sparse_in_file(ham(2),'store_ham2')
      endif

      call check_sparse_in_file(overlap,'store_overlap')
      call check_sparse_in_file(inv_overlap,'store_inv_overlap')
      call check_sparse_in_file(hub_overlap_matrix,'store_hub_overlap_matrix')
      call check_sparse_in_file(hub_overlap_tmatrix,'store_hub_overlap_tmatrix')

      if(.not.silent.and.pub_my_node_id==0) write(*,*) 'checking if file are presents is done'

      if(pub_dmft_invert_overlap)then
         write(*,*) 'INVERTING OVERLAP MATRIX'
         call invert_mat_sparser(overlap,inv_overlap,1,greenfunc=.false.)
      endif

      if(pub_dmft_norm_proj/=0) then
       if(abs(pub_dmft_nkpoints)>1)then
        write(*,*) 'ERROR combining renormalized Greens Function and k-points, not compatible for now'
        stop
       endif
       if(pub_dmft_norm_proj>0)then
        call normalize_projectors(hub_overlap_matrix,  proj_on_line_or_col=2, mat_t=hub_overlap_tmatrix)
       else
        call normalize_projectors(hub_overlap_matrix,  proj_on_line_or_col=2)
        call normalize_projectors(hub_overlap_tmatrix, proj_on_line_or_col=1)
       endif
      endif

      temp_times_pi = ABS(pub_dmft_temp)*PI
      if(pub_dmft_temp > 0.0_DP)then
         en_range = REAL(2*pub_dmft_points,kind=DP)*temp_times_pi
         en_start = temp_times_pi
         en_step  = 2.0_DP*temp_times_pi
      else
         en_range = pub_dmft_emax-pub_dmft_emin
         en_start = pub_dmft_emin
         en_step  = en_range/REAL(pub_dmft_points-1,kind=DP)
      endif

      if(.not.silent.and.pub_my_node_id==0)write(*,*) 'MY DMFT TEMPERATURE : ', pub_dmft_temp
      if(.not.silent.and.pub_my_node_id==0)write(*,*) 'MY FIRST FREQUENCY  : ', en_start
  
      en_step_over_pi = -1.0_DP * en_step / PI

      call define_energy_zero

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_scratch_dir
implicit none
logical :: checkscratch
real(8) :: lock

    if(pub_dmft_local_scratch)then
      lock=1.0
      inquire(file='scratchdir',exist=checkscratch)
      if(.not.checkscratch)then
       write(*,*) 'FILE scratchdir is NOT PRESENT'
       stop
      endif
      open(unit=81213,file='scratchdir')
      read(81213,'(a)') scratchdir
      close(81213)
      if(mpi_rank==1.and.pub_my_node_id==0)then
        write(*,*) 'THE SCRATCH DIRECTORY IS : ', trim(adjustl(scratchdir))
         ressys=system(" rm               "//trim(adjustl(scratchdir))//"/*  > /dev/null 2>&1 ")
         ressys=system(" cp sigma_output* "//trim(adjustl(scratchdir))//"    > /dev/null 2>&1 ")
      endif
      if(pub_my_node_id==0) call sync_nodes_by_shared_nfs(lock,'copy_sigma_output_scratch_dir',mpi_rank,mpi_size,longdelay=.true.)
      call comms_barrier
    endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine distribute_self_energy_to_atomic_components
implicit none

            do hub_atom = 1, pub_num_hub_atoms_on_node(pub_my_node_id)
            do hub_atomb= 1, pub_num_hub_atoms_on_node(pub_my_node_id)

             call set_loop_variable(hub_atom,hub_atomb)

             if(.not.cluster.and.hub_atom/=hub_atomb) goto 46
             if(cluster)then
              if( connection_table(hat,hatb) )then
              else
                goto 46
              endif
             endif

            !---rotate Self Func-------!
             if(pub_dmft_rotate_green) then
               if(channels/=channelsb) then
                 write(*,*) 'hub_atom, hub_atomb : ', hub_atom,hub_atomb
                 write(*,*) 'hat,hatb            : ', hat, hatb
                 write(*,*) 'table connection    : ', connection_table(hat,hatb)
                 write(*,*) 'channels,channelsb  : ', channels,channelsb
                 write(*,*) 'channels,channelsb are different, not compatible with local rotation'
                 stop
               endif
               call rotate_back_to_carth_basis_sigma(                                 &
                          & h_atoms(hat )%occupancy(1:channels,1:channels,is),        &
                          & h_atoms(hatb)%occupancy(1:channelsb,1:channelsb,is),      &
                          & site_self_energy_buffer(hub_atom,hub_atomb,1:channels,1:channelsb),channels, &
                          & rot_vec_angles(hub_atom,1:3,1:3),rot_vec_angles(hub_atomb,1:3,1:3),hat,hatb)
             endif
            !----------------------------!

             row_counter = 0
             row_proj_loop1: do row_proj = hub_proj_basis%first_on_atom(theatom), &
                                           hub_proj_basis%first_on_atom(theatom) + &
                                           hub_proj_basis%num_on_atom(theatom) - 1
                row_counter = row_counter + 1 ; col_counter = 0
                col_proj_loop1: do col_proj = hub_proj_basis%first_on_atom(theatomb), &
                                              hub_proj_basis%first_on_atom(theatomb) + &
                                              hub_proj_basis%num_on_atom(theatomb) - 1
                   col_counter = col_counter + 1
                   self_energy_element = site_self_energy_buffer(hub_atom,hub_atomb,row_counter,col_counter)
                   call sparse_put_element(self_energy_element, w_self_energy_v,row_proj,col_proj)
                enddo col_proj_loop1
             enddo row_proj_loop1

             if ((row_counter .ne. channels) .or. (col_counter .ne. channelsb) ) then
                write(stdout,'(a)') 'Error in hubbard_dmft_interface: row_counter or col_counter is not amounting to channels'
                write(stdout, *) 'channels=',channels,channelsb, ', row_counter=',row_counter,', col_counter=',col_counter
                call comms_abort
             endif

    46    enddo
          enddo

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine build_self_energy_from_extended_atoms
implicit none
            hubatoms1 : do hub_atom = 1, pub_num_hub_atoms_on_node(pub_my_node_id)
            hubatomsb : do hub_atomb= 1, pub_num_hub_atoms_on_node(pub_my_node_id)

             call set_loop_variable(hub_atom,hub_atomb)

             if(.not.cluster.and.hub_atom/=hub_atomb) goto 66
             if(     cluster)then
              if( (hat==hatb .and. dimer_table(hat,1)/=0) .or. dimer_connected(hat,hatb) ) then
              else
               goto 66
              endif
             endif

             call dump_info1

             if (channels .ne. (2*h_species(sp)%hub_ang_mom+1)) then
                write(stdout,'(a)') 'Error in hubbard_dmft_interface: hub_atom_hub_projs and hub_ang_mom mismatch'
                call comms_abort
             endif
             if(channels>size(site_self_energy_buffer,3).or.channelsb>size(site_self_energy_buffer,4).or.&
        &       hub_atom>size(site_self_energy_buffer,1).or.hub_atomb>size(site_self_energy_buffer,2))then
                write(*,*) 'ERROR onetep, buffer is too small, O boy I didnt think it could happen, my bad'
               stop
             endif
             if(pub_my_node_id==0.and.verboseall) write(*,*) 'reading self energy'
             call read_self_energy
             if(pub_my_node_id==0.and.verboseall) write(*,*) 'done'

    66      enddo hubatomsb
            enddo hubatoms1
            if(pub_my_node_id==0.and..not.silent)write(*,*) 'building self energy finished, my node / ien : ' ,pub_my_node_id,ien
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine show_hubbard_atom(mat,hub_atom,matout,showdiag,fullmat)
      use sparse,  only : sparse_convert
      implicit none
      integer                                       :: nn,hat,hub_atom,channels,theatom
      integer                                       :: i1,i2,j1,j2,i
      type(SPAM3),intent(in)                        :: mat 
      complex(kind=DP), allocatable, dimension(:,:) :: mat_square
      complex(kind=DP)                              :: matout(:,:)
      complex(8),optional                           :: fullmat(:,:)
      logical                                       :: showdiag
      real(8)                                       :: totdiag

      hat=pub_hub_atoms_on_node(hub_atom,pub_my_node_id)
      channels=hub_atom_hub_projs(hat)

      allocate(mat_square(pub_hub_nproj,pub_hub_nproj),stat=ierr)
      call sparse_convert(mat_square,mat)
   
     if(showdiag)then 
      write(*,*) 'total number of hubbard projections : ',pub_hub_nproj
      totdiag=0.0; do i=1,pub_hub_nproj; totdiag=totdiag+real(mat_square(i,i));enddo;
      write(*,*) 'real part of diagonal elements of hubbard projectors :'
      write(*,'(300f7.3)') (real(mat_square(i,i)),i=1,pub_hub_nproj)
      write(*,*) 'SUM IS : ' , totdiag
     endif

      if(size(matout,1)/=channels.or.size(matout,2)/=channels)then
       write(*,*) 'something is wrong in show_hubbard_atom'
       write(*,*) 'shape matout : ', shape(matout)
       write(*,*) 'channels : ', channels 
       stop
      endif

      theatom=pub_distr_atom(hub%orig(hat))
      i1=   hub_proj_basis%first_on_atom(theatom)
      i2=i1+hub_proj_basis%num_on_atom(theatom) - 1
      write(*,*) 'first and last projections on atom : ', i1,i2

      matout=mat_square(i1:i2,i1:i2)
      if(present(fullmat)) then
        if(any(shape(fullmat)-shape(mat_square)>0))then
         write(*,*) 'show hubbard proj, matrix shapes do not match'
         write(*,*) 'shape full mat : ', shape(fullmat)
         write(*,*) 'should be : ', shape(mat_square)
         stop
        endif
        fullmat=mat_square
      endif
      deallocate(mat_square)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine show_the_atomic_self
         if(verboseall)write(*,*) 'now wait other cpus'
         call comms_barrier
         if(show_matrices.and.ien==4)then
            if(pub_my_node_id==0) write(*,*) '-------showing selfenergy in Atom basis--------'
            call sparse_scale(w_self_energy_v,CMPLX(1000.,0,kind=8))
            call sparse_show_matrix(w_self_energy_v,show_elems=.true.)
            call sparse_scale(w_self_energy_v,CMPLX(1.d0/1000.d0,0,kind=8))
            if(pub_my_node_id==0) write(*,*) '-----------------------------------------------'
         endif
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine show_the_self_matrix
  call comms_barrier
        if(.false..and.show_matrices.and.ien==9)then
          if(pub_my_node_id==0) write(*,*) '-------showing selfenergy in NGWFs--------'
          call sparse_show_matrix(self_energy,show_elems=.true.)
          if(pub_my_node_id==0) write(*,*) '------------------------------------------'
        endif
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine invert_green_function
implicit none
          !---------------------------------!
           if(ien==pub_dmft_points-4.or.ien==pub_dmft_points-5)then
             call sparse_copy(inv_greenf_backup,inv_greenf)
             call sparse_copy(inv_greenf,overlap_cplx) 
           endif
          !---------------------------------!
           if(pub_dmft_lin_scaling) then
             if(verboseall)write(*,*) '================================'
             if(verboseall)write(*,*) 'INVERTING ON CPU WITH LINEAR SCALING SCHEME'
             if(nkpoints==1)then
               if(ien/=pub_dmft_points-3.and.ien/=pub_dmft_points-4.and.ien/=pub_dmft_points-5)then
                 call internal_invert_linear_scaling(greenf,energy,overlap_cplx,ham(is),self_energy,update=.false.,realorcomplexH=.true.)
               else
                 call sparse_copy(greenf,inv_overlap)
               endif
             else
               if(ien/=pub_dmft_points-3.and.ien/=pub_dmft_points-4.and.ien/=pub_dmft_points-5)then
                 call internal_invert_linear_scaling(greenf,energy,overlap_cplx,hamK,selfK,update=.false.,realorcomplexH=.false.)
               else
                 call sparse_copy(greenf,invK) 
               endif
             endif
             if(pub_dmft_temp < 0.0_DP) then
                call show_dos(greenf,overlap_cplx,energy)
             endif
           else
          !---------------------------------!
#ifdef GPU_SPEEDUP
           if(verboseall)write(*,*) '================================'
           if(verboseall)write(*,*) 'INVERTING FULL MATRIX ON GPU'
           if(use_gpu_onlyme)then
             call internal_invert_full_gpu(greenf,inv_greenf,overlap_cplx,energy,greenfunc=.true.)
           else
             call internal_invert_full(greenf,inv_greenf,greenfunc=.true.)
             if(pub_dmft_temp < 0.0_DP) then
                call show_dos(greenf,overlap_cplx,energy)
             endif
           endif
           if(verboseall)write(*,*) '================================'
#else
           if(verboseall)write(*,*) '================================'
           if(verboseall)write(*,*) 'INVERTING FULL MATRIX ON CPU (no open-mp)',is,ien,pub_dmft_points
           call internal_invert_full(greenf,inv_greenf,greenfunc=.true.)
           if(verboseall)write(*,*) '================================'
           if(pub_dmft_temp < 0.0_DP) then
              if(verboseall)write(*,*) 'showing DOS'
              call show_dos(greenf,overlap_cplx,energy)
              if(verboseall)write(*,*) 'done'
           endif
#endif
          endif
          !---------------------------------!

          if(ien==pub_dmft_points-4.or.ien==pub_dmft_points-5)then
             call matmul_in_dense_format(greenf_backup,greenf,inv_greenf_backup,aba=.true.)
             call sparse_copy(greenf,greenf_backup) ! in greenf we have now (S^-1*Self*S^-1) or (S^-1*H*S^-1)
          endif

          if(nkpoints>1)then
#ifdef debug
             call check_sparsec_(greenf)
             write(*,*) 'checking overlap_cplx'
             call check_sparsec_(overlap_cplx)
             write(*,*) 'checking ovK'
             call check_sparsec_(ovK)
             write(*,*) 'checking overlap'
             call check_sparser_(overlap)
             write(*,*) 'checking self energy'
             call check_sparsec_(self_energy)
             write(*,*) 'checking selfK'
             call check_sparsec_(selfK)
#endif
          endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine extract_green_function_elements
implicit none
integer :: kkklll,kkkk
integer :: numberHatom


         if(split)then
           numberHatom=pub_hub_nat
         else
           numberHatom=pub_num_hub_atoms_on_node(pub_my_node_id)
         endif

         hubatoms2 : do hub_atom = 1, numberHatom
         hubatoms2b: do hub_atomb= 1, numberHatom

             if(split)then
              call set_loop_variable_split(hub_atom,hub_atomb)
             else
              call set_loop_variable(hub_atom,hub_atomb)
             endif
             if(ien==energyzero.and.pub_my_node_id==0.and.channels==channelsb) goto 24
             if(.not.cluster.and.hub_atom/=hub_atomb) goto 36
             if(cluster)then
              if( connection_table(hat,hatb) .or. (pub_dmft_plot_all_proj.and.hub_atom==hub_atomb) )then
              else
               goto 36
              endif
             endif
          24 continue !bypass connection table to dump full proj GF!

             if( channels .ne. (2 * h_species(sp)%hub_ang_mom + 1) )then
                write(stdout,'(a)') 'Error in hubbard_dmft_interface: hub_atom_hub_projs and hub_ang_mom mismatch'
                call comms_abort
             endif
             if(channels>size(site_greenf_buffer,3).or.channelsb>size(site_greenf_buffer,4).or.&
        &       hub_atom>size(site_greenf_buffer,1).or.hub_atomb>size(site_greenf_buffer,2))then
               write(*,*) 'ERROR onetep, buffer is too small, O boy I didnt think it could happen, my bad'
               stop
             endif

             row_counter = 0
             row_proj_loop2: do row_proj = hub_proj_basis%first_on_atom(theatom), &
                                           hub_proj_basis%first_on_atom(theatom) + &
                                           hub_proj_basis%num_on_atom(theatom) - 1
                row_counter = row_counter + 1
                col_counter = 0
                col_proj_loop2: do col_proj = &
                                           hub_proj_basis%first_on_atom(theatomb), &
                                           hub_proj_basis%first_on_atom(theatomb) + &
                                           hub_proj_basis%num_on_atom(theatomb) - 1
                   col_counter = col_counter + 1
                   if(.not.split)then
                     call sparse_get_element(greenf_element,w_greenf_v,row_proj,col_proj)
                   else
                     greenf_element=split_wgv(row_proj,col_proj)
                   endif
                   site_greenf_buffer(hub_atom,hub_atomb,row_counter,col_counter) = greenf_element
                enddo col_proj_loop2
             enddo row_proj_loop2

            !---rotate Green's Func------!
             if(pub_dmft_rotate_green) then
               if(channels/=channelsb) then
                 write(*,*) 'channels,channelsb are different, not compatible with local rotation'
                 stop
               endif

              if(split)then
                call rotate_to_local_basis_green(                                                     &
                          & h_atoms(hat )%occupancy(1:channels,1:channels,is),                       &
                          & h_atoms(hatb)%occupancy(1:channelsb,1:channelsb,is),                     &
                          & site_greenf_buffer(hub_atom,hub_atomb,1:channels,1:channelsb),channels,  &
                          & rot_vec_angles_split(hub_atom,1:3,1:3),rot_vec_angles_split(hub_atomb,1:3,1:3),hat,hatb,pub_dmft_rotate_green) 
              else
                call rotate_to_local_basis_green(                                                     &
                          & h_atoms(hat )%occupancy(1:channels,1:channels,is),                       &
                          & h_atoms(hatb)%occupancy(1:channelsb,1:channelsb,is),                     &
                          & site_greenf_buffer(hub_atom,hub_atomb,1:channels,1:channelsb),channels,  &
                          & rot_vec_angles(hub_atom,1:3,1:3),rot_vec_angles(hub_atomb,1:3,1:3),hat,hatb,pub_dmft_rotate_green) 
              endif
             endif
            !----------------------------!
           
             if ((row_counter .ne. channels) .or. (col_counter .ne. channelsb) ) then
                write(stdout,'(a)') 'Error in hubbard_dmft_interface: row_counter or col_counter is not amounting to channels'
                write(stdout, *) 'channels=',channels,channelsb, ', row_counter=',row_counter,', col_counter=',col_counter
                call comms_abort
             endif

            !---dump local orbital DOS---!
            if(pub_dmft_temp<0.0_DP) then
             if(hub_atom==hub_atomb)then
             ! BUG CORRECTED JULY 8th 2014 
             if(     cluster ) then
              if(connection_table(hat,hatb))then
                goto 361
              endif
             endif
             ! END BUG
             if(.not.cluster)then
               361 continue
               if(.not.split)then
                 if(pub_my_node_id==0)then
                  if((ien==ffirst.and.one_k_one_file).or.(.not.one_k_one_file.and.ien==ffirst.and.ikpoint==1))then
                      open(unit=40000+hat+10003*is,file='LOCAL_DOS_atom_'//TRIM(ADJUSTL(toString(hat)))// &
                           & '_spin'//TRIM(ADJUSTL(toString(is)))// &
                           & '_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='replace')
                  else
                      open(unit=40000+hat+10003*is,file='LOCAL_DOS_atom_'//TRIM(ADJUSTL(toString(hat)))// &
                           & '_spin'//TRIM(ADJUSTL(toString(is)))// &
                           & '_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='old',position='append')
                  endif
                  mytrace=0.d0
                  do kk=1,channels
                    mytrace=mytrace+AIMAG(site_greenf_buffer(hub_atom,hub_atomb,kk,kk)*en_step_over_pi*norm_kpoints(ikpoint))
                  end do
                  write(40000+hat+10003*is,'(3i5,300f14.6)') ien,channels+1,pub_dmft_points,(real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132, &
                         & (AIMAG(site_greenf_buffer(hub_atom,hub_atomb,kk,kk)*en_step_over_pi*norm_kpoints(ikpoint)),kk=1,channels),mytrace
                  close(40000+hat+10003*is)
                 endif
               else
                 do kkkk=1,pub_total_num_nodes
                 if(kkkk==pub_my_node_id+1)then
                  if(ien==ffirst.and.kkkk==1.and. ( (.not.pub_dmft_splitk.and.(ikpoint==1.or.one_k_one_file))  .or.  (pub_dmft_splitk.and.ikpoint==splitk_first)  )   )then
                   open(unit=40000+hat+10003*is,file='LOCAL_DOS_atom_'//TRIM(ADJUSTL(toString(hat)))// &
                                                             & '_spin'//TRIM(ADJUSTL(toString(is)))// &                                  
                                                             & '_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),&
                                                             &  STATUS='REPLACE')
                  else
                   open(unit=40000+hat+10003*is,file='LOCAL_DOS_atom_'//TRIM(ADJUSTL(toString(hat)))// &
                                                             & '_spin'//TRIM(ADJUSTL(toString(is)))// &                                  
                                                             & '_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),&
                                                             &  position='append',STATUS='OLD' )
                  endif

                  mytrace=0.d0
                  do kk=1,channels
                    mytrace=mytrace+AIMAG(site_greenf_buffer(hub_atom,hub_atomb,kk,kk)*en_step_over_pi*norm_kpoints(ikpoint))
                  end do
                  write(40000+hat+10003*is,'(3i5,300f14.6)') ien,channels+1,pub_dmft_points,(real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132, &
                         & (AIMAG(site_greenf_buffer(hub_atom,hub_atomb,kk,kk)*en_step_over_pi*norm_kpoints(ikpoint)),kk=1,channels),mytrace
                  close(40000+hat+10003*is)

                 endif
                 call comms_barrier
                 enddo

               endif

             endif
             endif
            endif
            !----------------------------!

     36  enddo hubatoms2b
         enddo hubatoms2

         if(pub_my_node_id==0.and.ien==energyzero)then
           if(pub_dmft_temp>0.0_DP)then
             open(unit=4,file='projected_rotated_green_function_zero_energy_matsubara',form='unformatted')
           else
             open(unit=4,file='projected_rotated_green_function_zero_energy_real',form='unformatted')
           endif
           write(4) shape(site_greenf_buffer) 
           write(4) site_greenf_buffer
           write(4) pub_hub_nat
           write(4) pub_cell%a1%x*0.529177249, pub_cell%a1%y*0.529177249, pub_cell%a1%z*0.529177249
           write(4) pub_cell%a2%x*0.529177249, pub_cell%a2%y*0.529177249, pub_cell%a2%z*0.529177249
           write(4) pub_cell%a3%x*0.529177249, pub_cell%a3%y*0.529177249, pub_cell%a3%z*0.529177249
           do kkklll=1,pub_hub_nat
            write(4) coordinates_atom(kkklll)
           enddo 
           close(4)
         endif

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_down_green_function
implicit none
integer :: numberHatom

         numberHatom=pub_num_hub_atoms_on_node(pub_my_node_id)
        if(split)then
         numberHatom=pub_hub_nat 
        endif

        do hub_atom  = 1,numberHatom 
        do hub_atomb = 1,numberHatom
             if(split)then
              call set_loop_variable_split(hub_atom,hub_atomb)
             else
              call set_loop_variable(hub_atom,hub_atomb)
             endif
             if(.not.cluster.and.hub_atom/=hub_atomb) goto 37
             if(     cluster)then
              if( (hat==hatb .and. dimer_table(hat,1)/=0) .or. dimer_connected(hat,hatb) ) then
              else
               goto 37
              endif
             endif
            if(verboseall)then
             write(*,*) 'calling green info'
             write(*,*) 'hub_atom,hub_atomb   : ', hub_atom,hub_atomb
             write(*,*) 'hat,hatb             : ', hat,hatb
             write(*,*) 'DIAG sitegreenbuffer : ', (REAL(site_greenf_buffer(hub_atom,hub_atomb,i,i)),i=1,5)
            endif
            if(split)then
             call hubbard_greenf_info(hub_atom,hub_atomb,is,energy,ien,theatom,theatomb,sp,spb,hat,hatb, &
                                    & totn_atom_merge_split(hub_atom),totn_atom_merge_split(hub_atomb))
            else
             call hubbard_greenf_info(hub_atom,hub_atomb,is,energy,ien,theatom,theatomb,sp,spb,hat,hatb, & 
                                    & totn_atom_merge(hub_atom),totn_atom_merge(hub_atomb))
            endif
            if(verboseall)write(*,*) '.........done..........'
  37    enddo
        enddo

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine show_the_matrices
         if(show_matrices.and.ien==9)then
            if(pub_my_node_id==0) write(*,*) '-------showing greenf in Atom basis--------'
            call sparse_show_matrix(w_greenf_v,show_elems=.true.)
            if(pub_my_node_id==0) write(*,*) '-------------------------------------------'
         endif
         if(inv_projectors)then
          if(show_matrices.and.ien==9)then
            if(pub_my_node_id==0) write(*,*) '-----showing inv greenf in Atom basis------'
            call sparse_show_matrix(w_inv_greenf_v,show_elems=.true.)
            if(pub_my_node_id==0) write(*,*) '-------------------------------------------'
          endif
         endif
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_self_energy
implicit none

             if(ien==1.and..not.silent) write(*,*) 'reads self energy'
             !#######################################################################!
             !READ IN SELF-ENERGY FILE HERE
             !Self energy in orthogonalized basis (Self energy is diagonal)
             if(pub_dmft_ignore_type)then
              filename='sigma_output'//TRIM(ADJUSTL(toString(hat)))// &
                                & "_1_"//TRIM(ADJUSTL( toString(is)))
             else
              filename='sigma_output'//TRIM(ADJUSTL(toString(hat)))// &
                                & "_"//TRIM(ADJUSTL( toString(sp)))// &
                                & "_"//TRIM(ADJUSTL( toString(is)))
             endif
             if(pub_dmft_local_scratch)then
              filename=trim(adjustl(scratchdir))//trim(adjustl(filename))
             endif

             if( hat/=hatb ) filename=TRIM(ADJUSTL(filename))//'_dimer' 
             if(ien==1.and..not.silent) write(*,*) 'checking self-energy file : ', filename
             INQUIRE(file=filename,EXIST=check)

             !No self energy for extrema: to extract in DMFT the non-orthogonal factor for the tail
             if(check)then
                fopen=40010+theatom+2000*theatomb+110*hat+10101*is+2503*sp
                if(ien==1.and..not.silent) write(*,*) 'SELF ENERGY FILE IS PRESENT'
                if(ien==ffirst.or.reopen_sigma_file_each_step.or.force_reopen_file) then
                 if(ien==1.and..not.silent) write(*,*) 'opening file : ', filename
                 if(ien==1.and..not.silent) write(*,*) 'with unit    : ', fopen
                 open(unit=fopen,file=filename,form='UNFORMATTED')
                endif

                if(allocated(self_temp)) then
                 if(size(self_temp,1)/=totn_atom_merge(hub_atom).or.size(self_temp,2)/=totn_atom_merge(hub_atomb))then
                   deallocate(self_temp)
                 endif
                endif
                if(.not.allocated(self_temp)) allocate(self_temp(totn_atom_merge(hub_atom),totn_atom_merge(hub_atomb)))
                if(verboseall.and.pub_my_node_id==0)then
                  write(*,*) 'READING SIGMA WITH DIMENSIONS : ', shape(self_temp)
                  write(*,*) 'hub_atom and hub_atomb        : ', hub_atom,hub_atomb
                  write(*,*) 'tot atom merge                : ', totn_atom_merge(hub_atom),totn_atom_merge(hub_atomb)
                endif

                if(.not.(reopen_sigma_file_each_step.or.force_reopen_file))then
                  read(fopen) self_temp
                else
                 do kk=1,kien
                  read(fopen) self_temp
                 enddo
                endif

               !distribute self_temp in the blocks of site_self_energy_buffer :
                call distribute(self_temp,hub_atom,hub_atomb,totn_atom_merge(hub_atom),totn_atom_merge(hub_atomb))

                if(pub_dmft_temp .gt. 0.0_DP.and. &
              & (          ien==pub_dmft_points.or.ien==pub_dmft_points-2.or.ien==pub_dmft_points-3.or.ien==pub_dmft_points-4)) site_self_energy_buffer=0.d0
                if(pub_dmft_temp .lt. 0.0_DP.and. &
              & (ien==1.or.ien==pub_dmft_points.or.ien==pub_dmft_points-2.or.ien==pub_dmft_points-3.or.ien==pub_dmft_points-4)) site_self_energy_buffer=0.d0
                if(ien==llast.or.reopen_sigma_file_each_step.or.force_reopen_file)then
                  if(ien==1.and..not.silent) write(*,*) 'closing file : ', fopen
                  close(fopen)
                endif
             else
                if(ien==1.and..not.silent) write(*,*) 'NO SELF ENERGY FILE'
             endif
             !########################################################################!
             if(ien==1.and..not.silent) write(*,*) 'upfolds self-energy to NGWFs basis'

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 function coordinates_atom(jj)
 implicit none
 integer :: jj,j,hat,theatom
 real(8) :: coordinates_atom(3)
   !--------------------------!
   ! jj is in "hat" notations !
   !--------------------------!
   j = hub%orig(jj)
   coordinates_atom(1)=elements(j)%centre%x*0.529177249
   coordinates_atom(2)=elements(j)%centre%y*0.529177249
   coordinates_atom(3)=elements(j)%centre%z*0.529177249
 end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function atom_number_of_coordinates(coord)
implicit none
 integer :: i,j
 real(8) :: coord(3),v(3)
  j=0
  do i=1,pub_hub_nat
   v=coordinates_atom(i) 
   if( ALL ( abs(coord-v)<1.d-3 ) )then
     j=i;exit
   endif 
  enddo
  atom_number_of_coordinates=j
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function atom_number_of_coordinates_shifted(coord)
implicit none
integer :: i,j,k,kk,i1,i2,i3
real(8) :: coord(3),v(3)
real(8) :: real_lattice(3,3)

  real_lattice(1,1) = pub_cell%a1%x*0.529177249
  real_lattice(1,2) = pub_cell%a1%y*0.529177249
  real_lattice(1,3) = pub_cell%a1%z*0.529177249
  real_lattice(2,1) = pub_cell%a2%x*0.529177249
  real_lattice(2,2) = pub_cell%a2%y*0.529177249
  real_lattice(2,3) = pub_cell%a2%z*0.529177249
  real_lattice(3,1) = pub_cell%a3%x*0.529177249
  real_lattice(3,2) = pub_cell%a3%y*0.529177249
  real_lattice(3,3) = pub_cell%a3%z*0.529177249

  i1=INT(scalprod(coord,real_lattice(1,:))/norm_vector(real_lattice(1,:))**2)
  i2=INT(scalprod(coord,real_lattice(2,:))/norm_vector(real_lattice(2,:))**2)
  i3=INT(scalprod(coord,real_lattice(3,:))/norm_vector(real_lattice(3,:))**2)

  kk=0
  do i=-i1-3,-i1+3
  do j=-i2-3,-i2+3
  do k=-i3-3,-i3+3
    v=  dble(i)*real_lattice(1,:)
    v=v+dble(j)*real_lattice(2,:)
    v=v+dble(k)*real_lattice(3,:)
    kk=atom_number_of_coordinates(coord+v)
    if(kk/=0) then
     atom_number_of_coordinates_shifted=kk
     return
    endif
  enddo
  enddo
  enddo

  if(kk==0)then
     write(*,*) 'shifting coordinates, range : ', i1,i2,i3
     write(*,*) 'could not find coord in lattice, coord is : ', coord
     write(*,*) 'lattice vectors'
     do i=1,3
       write(*,'(i5,10f15.4)') i,real_lattice(i,:)
     enddo
     write(*,*) ' lattice is : '
     do i=1,pub_hub_nat; write(*,'(i5,10f15.4)') i,coordinates_atom(i); enddo;
     write(*,*) 'scalar product coord-a1 ', scalprod(coord,real_lattice(1,:))/norm_vector(real_lattice(1,:))**2
     write(*,*) 'scalar product coord-a2 ', scalprod(coord,real_lattice(2,:))/norm_vector(real_lattice(2,:))**2
     write(*,*) 'scalar product coord-a3 ', scalprod(coord,real_lattice(3,:))/norm_vector(real_lattice(3,:))**2
     stop
  endif

return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine correct_dimer_file
implicit none
logical :: check_coord
integer :: i,k1,k2,alpha1,alpha2,sss
real(8) :: coord1(3),coord2(3)

    INQUIRE(file='mask_dimer_coord',EXIST=check_coord)
    if(check_coord)then
     if(mpi_rank==1.and.pub_my_node_id==0) then
        ressys=system(" mv mask_dimer mask_dimer_original_coordinates ")
        ressys=system(" mv mask_dimer_coord mask_dimer_coord_original_coordinates ")
       open(unit=3031,file='mask_dimer')
        do i=1,pub_hub_nat
         if(dimer_table(i,1)>0) write(3031,*) dimer_table(i,1),dimer_table(i,2)
        enddo
       close(3031)
     endif
    endif

    INQUIRE(file='mask_uniform_coord',EXIST=check_coord)
    if(check_coord)then
     if(mpi_rank==1.and.pub_my_node_id==0) then
       open(unit=3031,file='mask_uniform_coord')
       open(unit=3032,file='mask_uniform_corrected')
       open(unit=3033,file='mask_uniform')
       do i=1,pub_hub_nat
         read(3031,*)   coord1
         read(3031,*)   coord2
         read(3033,*)   alpha1,alpha2
         if(pub_dmft_in_bohr)then
           coord1=coord1*0.529177249 
           coord2=coord2*0.529177249 
         endif
         sss = sign(1, alpha2 )
         write(3032,*)     atom_number_of_coordinates_shifted(coord1),&
                     & sss*atom_number_of_coordinates_shifted(coord2)
       enddo
       close(3031)
       close(3032)
       close(3033)
        ressys=system(" mv mask_uniform       mask_uniform_original       ")
        ressys=system(" mv mask_uniform_coord mask_uniform_coord_original ")
        ressys=system(" cp mask_uniform_corrected mask_uniform ")
     endif
    endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine define_cluster
implicit none
logical :: allsites,check,check_coord,check_merge_file
integer :: ii,jj,ijk_
integer :: ll,lll,chan,lll_,llll_
real(8) :: v_temp(3),dtemp,coord(3),mat_temp(3,3),v1(3),v2(3)

     allsites=.true.

     if(allocated(rot_vec_angles))       deallocate(rot_vec_angles)
     if(allocated(rot_vec_angles_split)) deallocate(rot_vec_angles_split)

     allocate(rot_vec_angles(pub_hub_nat,3,3))
     allocate(rot_vec_angles_split(pub_hub_nat,3,3))

     write(*,*) 'THERE ARE [X] HUBBARD ATOMS ON MY NODE [Y] : ', pub_num_hub_atoms_on_node(pub_my_node_id),pub_my_node_id

     rot_vec_angles=0.0d0; rot_vec_angles_split=0.0d0
  
     INQUIRE(file='mask_local_rotations',EXIST=check)
  
     call comms_barrier
     if(check)then
      open(unit=236,file='mask_local_rotations',form='unformatted')
      !----------------------------------------------------------------------------------!
      do ijk_=1,pub_hub_nat

       read(236,end=15) lll,coord,mat_temp
       if(pub_dmft_in_bohr)then
         coord=coord*0.529177249
       endif
       if(lll==0)then
         allsites=.false.
         goto 331
       endif

       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'Rotation for hubbard atom : ', ijk_

       j=atom_number_of_coordinates_shifted(coord)
       if(j/=ijk_)then
         write(*,*) 'WARNING, a shift of coordinate was done by onetep!!!!'
         write(*,*) 'might give problems with rotations when in mpi'
       endif 

       if(j<1.or.j>pub_hub_nat)then
         write(*,*) 'hubbard atom found is wrong, j = ', j
         write(*,*) 'total hubbard atoms is : ', pub_hub_nat
         stop
       endif

       rot_vec_angles_split(j,:,:)=mat_temp

       j=hubard_on_my_node(j) !now "j" is in hub_atom notation

       if(j/=0)then
          rot_vec_angles(j,:,:)=mat_temp

          chan=hub_atom_hub_projs(pub_hub_atoms_on_node(j,pub_my_node_id))  !hub_projs take 1..nat argument, so hat notation

          if(.not.silent.and.pub_my_node_id==0)then
             write(*,'(a,i5)') 'ATOM [x] has rotation : ',lll
             do lll_=1,3
               write(*,'(3f10.3)') (rot_vec_angles(j,lll_,llll_),llll_=1,3)
             enddo
             write(*,*) 'j,pub_my_node_id                    : ', j,pub_my_node_id
             write(*,*) 'checking rotations for [x] channels : ', chan
             write(*,*) 'TEST, B*B^T (maxval): ', maxval(abs(   matmul(cmp_rot(rot_vec_angles(j,:,:),(chan-1)/2,chan,j), &
                                                           & transpose(cmp_rot(rot_vec_angles(j,:,:),(chan-1)/2,chan,j) )) ))
             write(*,*) 'TEST, B*B^T (minval): ', minval(abs(   matmul(cmp_rot(rot_vec_angles(j,:,:),(chan-1)/2,chan,j), &
                                                           & transpose(cmp_rot(rot_vec_angles(j,:,:),(chan-1)/2,chan,j) )) ))
          endif


       endif

 331  enddo
      !----------------------------------------------------------------------------------!
      if(.false.)then
        15 continue
        write(*,*) 'FILE MASK_LOCAL_ROTATIONS HAS NOT ENOUGH ENTRIES'
        stop 
      endif 
      do lll_=1,pub_num_hub_atoms_on_node(pub_my_node_id)
       if(maxval(abs(rot_vec_angles(lll_,:,:)))<1.d-4.and.allsites)then
         write(*,*) ' WARNING WARNING WARNING '
         write(*,*) ' WARNING WARNING WARNING rotation matrix was not defined for site : ', lll_
         v_temp=coordinates_atom( pub_hub_atoms_on_node(lll_,pub_my_node_id) )
         write(*,*) ' coordinates of this site                                         : ', v_temp
       endif
      enddo
      close(236)
     endif

     if(cluster)then
       if(allocated(dimer_table)) deallocate(dimer_table,merge_table,merge_neigh)
       allocate(dimer_table(pub_hub_nat,2),merge_table(pub_hub_nat,nmerge),merge_neigh(pub_hub_nat))
       if(.not.silent.and.pub_my_node_id==0)then
          write(*,*) 'NUMBER OF PROJECTIONS ON ATOMS'
          do j=1,pub_hub_nat
             write(*,'(a,2i7)') 'ATOM, NPROJ : ', j,hub_atom_hub_projs(j)
          enddo
          write(*,*) 'READING MASK_DIMER'
       endif

       open(unit=111,file='mask_dimer')
       dimer_table=0
       do i=1,pub_hub_nat
         read(111,*,end=898) j,jj 
         if(j>0)then
           dimer_table(j,1)=j
           dimer_table(j,2)=jj
           if(.not.silent.and.pub_my_node_id==0) write(*,'(a,200i4)') 'DIMER TABLE : ', (dimer_table(j,jj),jj=1,2)
         endif
         if(jj>0.and.j>0.and.dimer_table(j,1)>0.and.dimer_table(j,2)>0)then
          dtemp=norm_vector(coordinates_atom(dimer_table(j,1))-coordinates_atom(dimer_table(j,2)))
          if(.not.silent.and.pub_my_node_id==0) write(*,*) 'DIMER BOND LENGTH     : ', dtemp
         endif
       enddo
  898  continue
       close(111)

       INQUIRE(file='mask_dimer_coord',EXIST=check_coord)
       if(check_coord)then
        dimer_table=0
        open(unit=3031,file='mask_dimer_coord')
        do i=1,pub_hub_nat
         read(3031,*) v1
         read(3031,*) v2
         if(pub_dmft_in_bohr)then
          v1=v1*0.529177249 
          v2=v2*0.529177249 
         endif
                        j=atom_number_of_coordinates_shifted(v1)
         dimer_table(j,1)=j
         dimer_table(j,2)=atom_number_of_coordinates_shifted(v2)
        enddo
        close(3031)
       endif

       INQUIRE(FILE='mask_merge',EXIST=check_merge_file)
       merge_table=0
       if(check_merge_file) open(unit=111,file='mask_merge')
       merge_neigh=0
       do i=1,pub_hub_nat
         if(check_merge_file) then
            read(111,*,end=899) (merge_table(i,j),j=1,nmerge)
         else
            merge_table(i,:)=0
            merge_table(i,1)=i
         endif
         if(.not.silent.and.pub_my_node_id==0) write(*,'(a,200i4)') 'MERGE TABLE : ', (merge_table(i,j),j=1,nmerge)
         merge_neigh(i)=0
         do j=2,nmerge
          if(merge_table(i,j)==0)then
            merge_neigh(i)=j-1
            exit
          else
            merge_neigh(i)=j
          endif
         enddo
       enddo
899    continue
       if(check_merge_file) close(111)

       if(.not.silent.and.pub_my_node_id==0) write(*,*) 'build connection table'
       if(allocated(connection_table)) deallocate(connection_table)
       allocate(connection_table(pub_hub_nat,pub_hub_nat))
       connection_table=.false.
       do i=1,pub_hub_nat
         if(dimer_table(i,1)/=0)then
         do i2=1,2
         do k=1,merge_neigh(i)
          l=merge_table(i,k)
          jj=dimer_table(i,i2)
          if(jj/=0)then 
           do k2=1,merge_neigh(jj)
            l2=merge_table(jj,k2)
            connection_table(l,l2)=.true.
            connection_table(l2,l)=.true.
            if(hubard_on_my_node(l)==0.or.hubard_on_my_node(l2)==0)then
              write(*,*) 'cluster DMFT, connection table not compatible with distribution of atoms on mpi nodes, node_id : ', pub_my_node_id
              stop
            endif
           enddo
          endif
         enddo
         enddo
         endif
       enddo
      !connection_table=.true. !TESTING-DEBUG
       if(.not.silent.and.pub_my_node_id==0)then
       write(*,*) 'TABLE'
       write(*,'(100i0)') (l,l=1,pub_hub_nat)
       do i=1,pub_hub_nat
         write(*,*) i,(connection_table(i,l),l=1,pub_hub_nat)
       enddo
       endif
      else
        if(.not.silent.and.pub_my_node_id==0) write(*,*) 'MASK_DIMER NOT PRESENT,GO FOR SINGLE SITE DMFT'
      endif

      if(.not.silent.and.pub_my_node_id==0)then
      write(*,*) '===================================='
      write(*,*) '===================================='
      write(*,*) '===================================='
      write(*,*) 'HUBBARD ATOMS, POSITIONS:'
      open(unit=1121214,file='positions_hubbard_atoms')
      do ii=1,pub_hub_nat
        if( hub%orig(ii) /= pub_orig_atom(pub_distr_atom(hub%orig(ii))) )then
          write(*,*) 'WARNING WARNING WARNING :'
          write(*,*) '           you might want to check whats going on in coordinates_atom function'
          write(*,*) '           position definition and access to element(j) array ill defined '
          write(*,*) 'this is hub%orig(ii)                                     : ', hub%orig(ii)
          write(*,*) 'this is  pub_orig_atom(pub_distr_atom(hub%orig(ii)))     : ', pub_orig_atom(pub_distr_atom(hub%orig(ii)))
        endif
        write(*,*) 'hubbard atom : ', ii,pub_hub_nat
        v_temp=coordinates_atom(ii)
        write(1121214,'(i5,i5,10f15.4)') ii,hub%species_number(ii),v_temp
      enddo
      close(1121214)
      write(*,*) '===================================='
      write(*,*) '===================================='
      write(*,*) '===================================='
      endif

return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine build_my_dimers
implicit none
integer :: uu(3)

write(*,*) 'building dimers'

if(pub_my_node_id==0)then
 !--------------!
 ! print dimers !
 !--------------!
 if(cluster)then
  if(.not.silent.and.pub_my_node_id==0) write(*,*) 'dump cluster files'

  open(unit=333,file='compound_sites')
  open(unit=334,file='compound_dimers')
  open(unit=335,file='compound_distances')
  open(unit=336,file='compound_blocks')

  if(.not.split)then
   do hub_atom = 1, pub_num_hub_atoms_on_node(pub_my_node_id)
    hat = pub_hub_atoms_on_node(hub_atom,pub_my_node_id);
    do hub_atomb= 1, pub_num_hub_atoms_on_node(pub_my_node_id)
     hatb = pub_hub_atoms_on_node(hub_atomb,pub_my_node_id);
     if( dimer_table(hat,2)==hatb ) then
       write(336,*) hub_atom_hub_projs(hat),hub_atom_hub_projs(hatb)
     endif
    enddo
   enddo
  else
  do hub_atom = 1, pub_hub_nat
    hat = hub_atom
    do hub_atomb= 1, pub_hub_nat
     hatb = hub_atom
     if( dimer_table(hat,2)==hatb ) then
       write(336,*) hub_atom_hub_projs(hat),hub_atom_hub_projs(hatb)
     endif
    enddo
   enddo
  endif

  do hub_atom = 1, pub_num_hub_atoms_on_node(pub_my_node_id)
  hat   = pub_hub_atoms_on_node(hub_atom,pub_my_node_id); 
  write(333,'(i3,10f11.4)')  hub_atom, coordinates_atom(hat) 
  do hub_atomb= 1, pub_num_hub_atoms_on_node(pub_my_node_id)
  hatb  = pub_hub_atoms_on_node(hub_atomb,pub_my_node_id); 
  dist=distance(hat,hatb,uu); 
  if(hub_atom/=hub_atomb) write(335,'(a,i4,i4,f10.4)') 'ATOMS : ',hub_atom,hub_atomb,dist
  if( dimer_table(hat,2)==hatb ) then
     write(334,'(a,i4,10f11.4)')' A= : ' ,hub_atom, coordinates_atom(hat)
     write(334,'(a,i4,10f11.4)')' B= : ', hub_atomb,coordinates_atom(hatb)
  endif
  enddo
  enddo
  close(333)
  close(334)
  close(335)
  close(336)
  if(.not.silent.and.pub_my_node_id==0) write(*,*) 'done'
 endif
 !..........group atoms..............!
 if(.not.silent.and.pub_my_node_id==0) write(*,*) 'dump compound_tot_orb file'
 if(.not.split)then
   open(unit=337,file='compound_tot_orb')
   do hub_atom= 1, pub_num_hub_atoms_on_node(pub_my_node_id)
    write(337,*) totn_atom_merge(hub_atom)
   enddo
   close(337)
 else
   open(unit=337,file='compound_tot_orb')
   do hub_atom= 1, pub_hub_nat
    write(337,*) totn_atom_merge_split(hub_atom)
   enddo
   close(337)
 endif
 if(.not.silent.and.pub_my_node_id==0) write(*,*) 'done'
 !...................................!
endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine read_write_hamiltonian
   implicit none
   real(8),allocatable :: temp_alloc(:,:)
   integer             :: iboundary

   if(.not.check)then
     if(allocated(eigen_en)) deallocate(eigen_en,n_occ)
     allocate(eigen_en(size(eigen_en_input,1),size(eigen_en_input,2)),n_occ(size(n_occ_input,1)))
     eigen_en=eigen_en_input
     n_occ=n_occ_input
     if(.not.silent.and.pub_my_node_id==0) write(*,*) 'file is not present, write current data'
     if(pub_my_node_id==0) open(unit=20001,file='store_eigen',form='unformatted')
     if(.not.silent.and.pub_my_node_id==0) write(*,*) 'dumping file store_eigen,shapes of array are: ',size(eigen_en_input,1),size(eigen_en_input,2),size(n_occ_input,1)
     if(pub_my_node_id==0) write(20001) size(eigen_en_input,1),size(eigen_en_input,2),size(n_occ_input,1)
     if(pub_my_node_id==0) write(20001) n_occ,eigen_en
     if(pub_my_node_id==0) close(20001)
     call comms_barrier
     do hub_atom = 1, pub_num_hub_atoms_on_node(pub_my_node_id)
           hat = pub_hub_atoms_on_node(hub_atom,pub_my_node_id)
           channels=hub_atom_hub_projs(hat)
           open(unit=20001,file='store_occupancy'//trim(adjustl(toString(hat))),form='unformatted')
           write(20001) channels,num_spins
             write(20001) h_atoms(hat)%occupancy(1:channels,1:channels,1)
           if(num_spins==2)then
             write(20001) h_atoms(hat)%occupancy(1:channels,1:channels,2)
           endif
           close(20001)
     enddo
     if(.not.silent.and.pub_my_node_id==0) write(*,*) '...done...',pub_my_node_id
     call comms_barrier
   else
     if(.not.silent.and.pub_my_node_id==0) write(*,*) 'file is present, read from this file'
     open(unit=20001,file='store_eigen',form='unformatted')
     if(.not.silent.and.pub_my_node_id==0) write(*,*) 'reading file store_eigen'
     read(20001) n1_eigen,n2_eigen,n3_eigen
     if(allocated(eigen_en)) deallocate(eigen_en,n_occ)
     allocate(eigen_en(n1_eigen,n2_eigen),n_occ(n3_eigen))
     read(20001) n_occ,eigen_en

     if(size(eigen_en,2)==1.and.num_spins==2)then
      write(*,*) 'WARNING something wrong with num_spins, you do DMFT with num_spins=2'
      write(*,*) 'using store_files from a onetep simulations done with num_spins=1'
      write(*,*) 'For now, I adapt the size of the arrays for you...'
      if(allocated(temp_alloc)) deallocate(temp_alloc)
      allocate(temp_alloc(size(eigen_en,1),size(eigen_en,2)))
      temp_alloc=eigen_en
      deallocate(eigen_en)
      allocate(eigen_en(size(temp_alloc,1),2))
      eigen_en(:,1)=temp_alloc(:,1)
      eigen_en(:,2)=temp_alloc(:,1)
      deallocate(temp_alloc)
     endif

    if(.not.silent.and.pub_my_node_id==0)then
      write(*,*) 'done...'
      write(*,*) 'n_occ 1-2: ',n_occ
      write(*,*)
      write(*,*) 'eigen_en',eigen_en
     endif

     close(20001)


     if(split)then
       iboundary=pub_hub_nat
     else
       iboundary=pub_num_hub_atoms_on_node(pub_my_node_id)
     endif

     do hub_atom =1,iboundary
           if(split)then
            hat=hub_atom 
           else
            hat=pub_hub_atoms_on_node(hub_atom,pub_my_node_id) 
           endif
           channels=hub_atom_hub_projs(hat)
           open(unit=20001,file='store_occupancy'//trim(adjustl(toString(hat))),form='unformatted')
           read(20001) channels,ii
           read(20001) h_atoms(hat)%occupancy(1:channels,1:channels,1)
           if(ii/=num_spins)then
             write(*,*) 'not the right number of spin species, should be : ', ii
             write(*,*) 'check variable dmft_paramagnetic or number of spin species in onetep input file'
           endif
           if(num_spins==2)then
             if(ii==2)then
              if(.not.silent.and.pub_my_node_id==0) write(*,*) 'spin down too'
              read(20001) h_atoms(hat)%occupancy(1:channels,1:channels,2)
             else
              h_atoms(hat)%occupancy(1:channels,1:channels,2)=h_atoms(hat)%occupancy(1:channels,1:channels,1)
             endif
           endif
           close(20001)
     enddo
   endif

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine dump_some_info_on_h_atom
   implicit none
   real(8) :: temp_dens(100,100) 

   if(silent) return

   if(pub_my_node_id==0)then
    write(*,*) '========================================'
    write(*,*) '========================================'
    write(*,*) 'SPECTRUM'
    write(*,*) 'n_occ              : ', n_occ(:)
    write(*,*) 'shape eigen        : ', shape(eigen_en)
     write(*,*) 'eigen up          : ', eigen_en(n_occ(1),1)
    if(num_spins==2)then
     write(*,*) 'eigen dn          : ', eigen_en(n_occ(2),2)
    endif
    write(*,*) 'pub_dmft_smear     : ', pub_dmft_smear
    write(*,*) 'pub_dmft_smear_T   : ', pub_dmft_smear_T
    write(*,*) 'pub_dmft_smear_w   : ', pub_dmft_smear_w
    write(*,*) 'pub_dmft_smear_eta : ', pub_dmft_smear_eta
    open(unit=41,file='density_spin_up_non_rotated')
    open(unit=42,file='density_spin_up_rotated')
    do hub_atom = 1, pub_num_hub_atoms_on_node(pub_my_node_id)
     hat=pub_hub_atoms_on_node(hub_atom,pub_my_node_id)
     channels=hub_atom_hub_projs(hat)
     if(channels>size(temp_dens))then
      write(*,*) 'ERROR dmft module, temp_dens too small'
      stop
     endif
     write(*,*) 'ATOM [x] with SPIN [x] occupancy matrix : ',hat,1
      do i=1,channels
        write(*,'(200f10.5)') ( h_atoms(hat)%occupancy(i,j,1),j=1,channels )
      enddo
      write(41,'(200f8.4)') ( h_atoms(hat)%occupancy(i,i,1),i=1,channels ) , sum( (/( h_atoms(hat)%occupancy(i,i,1),i=1,channels )/) )
      call diago_occupancy(h_atoms(hat)%occupancy(1:channels,1:channels,1),channels,hat,rot_vec_angle=rot_vec_angles(hub_atom,1:3,1:3),tmp_dens=temp_dens)
      write(42,'(200f8.4)') ( temp_dens(i,i),i=1,channels ) , sum( (/( temp_dens(i,i),i=1,channels )/) )
    enddo
    close(41) 
    close(42)
    write(*,*) '========================================'
    write(*,*) '========================================'
    if(num_spins==2)then
     open(unit=41,file='density_spin_dn_non_rotated')
     open(unit=42,file='density_spin_dn_rotated')
     do hub_atom = 1, pub_num_hub_atoms_on_node(pub_my_node_id)
      hat=pub_hub_atoms_on_node(hub_atom,pub_my_node_id)
      channels=hub_atom_hub_projs(hat)
      write(*,*) 'ATOM [x] with SPIN [x] occupancy matrix : ',hat,2
      do i=1,channels
        write(*,'(200f10.5)') ( h_atoms(hat)%occupancy(i,j,2),j=1,channels )
      enddo
      write(41,'(200f8.4)') ( h_atoms(hat)%occupancy(i,i,2),i=1,channels ) , sum( (/( h_atoms(hat)%occupancy(i,i,2),i=1,channels )/) )
      call diago_occupancy(h_atoms(hat)%occupancy(1:channels,1:channels,2),channels,hat,rot_vec_angle=rot_vec_angles(hub_atom,1:3,1:3),tmp_dens=temp_dens)
      write(42,'(200f8.4)') ( temp_dens(i,i),i=1,channels ) , sum( (/( temp_dens(i,i),i=1,channels )/) )
     enddo
     close(41)
     close(42)
    endif
    write(*,*) '========================================'
    write(*,*) '========================================'
   endif

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine clean_dmft_energy
#ifdef GPU_SPEEDUP
  use fortran_cuda, only : matmulcuda_r
#endif
    implicit none
    integer    :: ii,i,j,m,jj,num_denskern,np
    integer    :: iiistart,iiistep,iistart,iistep
    real(8)    :: totN,totE,beta,chem,lambda
    real(8)    :: ccback,cc,slope_mu,edc_correction
    real(8)    :: ttscreen,mymu,dtempd,dtempdd
    complex(8) :: dtempddc,ccc,ccbackc

    write(*,*) 'FINALIZING DMFT ENERGY AND KERNEL'

    ii=ngwf_basis%num
    beta=1.d0/pub_dmft_temp
    chem=fermi_e(is)+pub_dmft_chem_shift

    if(nkpoints==1)then
     matrixDensKernel  = matrixDensKernel  * pub_dmft_temp 
    else
     matrixDensKernelc = matrixDensKernelc * pub_dmft_temp 
    endif

     !--------------------------!
    if(mpi_size/=1.and..not.breaknfscom)then

#ifndef debug
      call sync_nodes_by_shared_nfs(dmft_Energy(is),'summing_dmft_energy'//trim(adjustl(tostring(is))),mpi_rank,mpi_size)
      if(nkpoints==1)then
        call sync_nodes_by_shared_nfs_mat(matrixDensKernel  ,'sum_matrixDensKernel' ,mpi_rank,mpi_size,size(matrixDensKernel,1))
      else
        call sync_nodes_by_shared_nfs_matc(matrixDensKernelc,'sum_matrixDensKernelc',mpi_rank,mpi_size,size(matrixDensKernelc,1))
      endif
#else
      open(file='dens_kern_tmp'//trim(adjustl(tostring(mpi_rank))),unit=981777+mpi_rank,status="replace",position='rewind',form='unformatted') 
      if(nkpoints==1)then
        write(981777+mpi_rank) matrixDensKernel,dmft_Energy(is) 
        write(981777+mpi_rank) matrixDensKernel,dmft_Energy(is) !to flush the file...
      else
        write(981777+mpi_rank) matrixDensKernelc,dmft_Energy(is) 
        write(981777+mpi_rank) matrixDensKernelc,dmft_Energy(is) !to flush the file...
      endif
        close(981777+mpi_rank)
#ifdef DMFT_SAFE_NFS
       ressys=system("sleep 1")
#endif
      do
#ifdef DMFT_SAFE_NFS
          ressys=system("sleep 0.5")
#else
          ressys=system("sleep 0.1")
#endif
          ressys=system(" echo 'I am done'  > onetep_dmft_confirmation_dmft_density_kernel_"//trim(adjustl(tostring(mpi_rank))))

         if(mpi_rank==1)  ressys=system("echo `ls -l onetep_dmft_confirmation_dmft_density_kernel_*  2>&1  | grep -v 'No such' | wc -l` > onetep_dmft_confirmation_dmft_density_count")

         open(unit=10,file='onetep_dmft_confirmation_dmft_density_count')
         jj=0
         read(10,*,end=67) jj
         67 continue
         close(10)
         write(*,*) 'N PROCESS ARE DONE = ',jj
         if(jj==mpi_size) then
           write(*,*) 'ALL  ARE DONE'
           exit
         endif
      enddo
#ifndef DMFT_SAFE_NFS
       ressys=system("sleep 1")
#else
       ressys=system("sleep 2")
#endif
      if(nkpoints==1)then
       matrixDensKernel=0.d0
      else
       matrixDensKernelc=0.d0
      endif
      dmft_Energy(is)=0.d0
      do i=1,mpi_size
         open(file='dens_kern_tmp'//trim(adjustl(tostring(i))),unit=981777+i,action='read',position='rewind',form='unformatted')
         if(nkpoints==1)then
          read(981777+i) tmpr,cc
          matrixDensKernel  = matrixDensKernel  + tmpr 
         else
          read(981777+i) tmpc,cc
          matrixDensKernelc = matrixDensKernelc + tmpc 
         endif
         close(unit=981777+i)
         dmft_Energy(is)  = dmft_Energy(is)  + cc
      enddo
#endif
    endif

     !--------------------------!
   if(ikpoint<=nkpoints)then
    if(abs(norm_kpoints(ikpoint))>1.d-4)then
      !first sync frequencies, then add tail for this k point
      if(nkpoints==1)then
        matrixDensKernel  = matrixDensKernel  + matrixContributionTail
      else
        matrixDensKernelc = matrixDensKernelc + matrixContributionTailc 
      endif
    endif
   endif

    if(pub_my_node_id==0)then
     write(*,*) '==============================================='
     write(*,*) 'DMFT DENSITY KERNEL'
     write(*,*) 'temperature : ', pub_dmft_temp
    endif

    totN=0.d0

    if( optimisation_parallel )Then
     iistart=pub_my_node_id;iistep=pub_total_num_nodes
    else
     iistart=0;iistep=1
    endif

    if(nkpoints==1)then
     do i=iistart+1,ii,iistep
      do m=1,ii
       totN = totN+matrixDensKernel(i,m)*matrixOverlap(m,i) 
      enddo
     enddo
    else
     do i=iistart+1,ii,iistep
      do m=1,ii
       totN = totN+matrixDensKernelc(i,m)*matrixOverlapc(m,i) 
      enddo
     enddo
    endif

    if( optimisation_parallel .and. pub_total_num_nodes>1) call comms_reduce( 'SUM',totN  )


#ifdef debug
    if(nkpoints>1)then
     write(*,*) 'MINVAL-MAXVAL MATRIXAsig : ', minval(abs(matrixAsig)),maxval(abs(matrixAsig))
     write(*,*) 'MINVAL-MAXVAL MATRIXAAC  : ', minval(abs(matrixAAc)),maxval(abs(matrixAAc))
     write(*,*) 'MINVAL-MAXVAL DENSKERNEL : ', minval(abs(matrixDensKernelc)),maxval(abs(matrixDensKernelc))
     write(*,*) 'MINVAL-MAXVAL DENSK TAIL : ', minval(abs(matrixContributionTailc)),maxval(abs(matrixContributionTailc))
     write(*,*) 'MINVAL-MAXVAL overlapc   : ', minval(abs(matrixOverlapc)),maxval(abs(matrixOverlapc))
    endif
#endif

    if(pub_my_node_id==0) write(*,*) 'TOTAL CHARGE OBTAINED FROM DMFT    : ', totN
    if(pub_my_node_id==0) write(*,*) 'TARGET N                           : ', target_N
    if(pub_my_node_id==0) write(*,*) 'K POINT                            : ', ikpoint
    if(ikpoint<=nkpoints)then
     if(pub_my_node_id==0) write(*,'(a,3f10.4)') 'K PT                    : ', kpt(ikpoint,:)
    else
     if(pub_my_node_id==0) write(*,'(a,3f10.4)') 'K PT OUT OF RANGE'
    endif
    if(pub_my_node_id==0) write(*,*) 'TOTAL K POINTS                     : ', nkpoints

    if(ikpoint<=nkpoints)then
     contrib_N=totN*norm_kpoints(ikpoint)
    else
     contrib_N=0.d0
    endif

    if(pub_dmft_splitk)then
     dmft_splitk_iter=dmft_splitk_iter+1
     if(pub_my_node_id==0) write(*,*) 'this is sync number : ', dmft_splitk_iter
     if(pub_my_node_id==0) call sync_nodes_by_shared_nfs(contrib_N,'summing_up_contribN'&
                           & //trim(adjustl(tostring(dmft_splitk_iter))),mpi_rank,mpi_size)
     call comms_barrier
     call comms_bcast(0,contrib_N)
     if(pub_my_node_id==0) write(*,*) 'TOTAL CHARGE THIS BATCH K POINTS   : ', contrib_N
    endif

    vectorNmu(nmu_step) = vectorNmu(nmu_step) + contrib_N

    if(pub_my_node_id==0) write(*,*) 'MU,N                               : ', vectorMu(nmu_step),vectorNmu(nmu_step)

#ifdef debug
    totN=0.d0
    if(nkpoints==1)then
      do i=1,ii
       do m=1,ii
        totN = totN + matrixContributionTail(i,m)*matrixOverlap(m,i)
       enddo
      enddo
    else
      do i=1,ii
       do m=1,ii
        totN = totN + real( matrixContributionTailc(i,m)*matrixOverlapc(m,i) )
       enddo
      enddo
    endif
    if(pub_my_node_id==0) write(*,*) 'TOTAL CHARGE OBTAINED BY TAIL      : ', totN
#endif

    if(nmu_step < nmu) goto 66778

    if(pub_my_node_id==0) write(*,*) '<HU> FROM DYNAMICAL EFFECT         : ', dmft_Energy(is)

   ! --> 1/2 Tr ( Sigma_an * G_an ) !
    if(.not.use_gpu.and..not.use_gpu_some_only)then
     cc=0.d0;ccc=0.d0
     if(breaknfscom)then
       iiistart=0;iiistep=1
     else
       iiistart=mpi_rank-1;iiistep=mpi_size
     endif

     if( optimisation_parallel )Then
     if( (.not.breaknfscom.and.mpi_size==1) .or. (breaknfscom) )then
       if(iiistart>0)then; write(*,*) 'ERROR parallelization conflict';stop; endif
       iiistart=pub_my_node_id;iiistep=pub_total_num_nodes
     endif
     endif

     do m=iiistart+1,ii,iiistep
      lambda = vectorAA(is,m) - chem
      dtempd = 1.d0 / ( -lambda ) * ( 1.d0 - DEXPc(-beta*lambda) ) / ( 1.d0 + DEXPc(-beta*lambda) )
      if(nkpoints==1)then
        do i=1,ii
         dtempdd = dtempd * matrixAA(is,i,m) 
         do j=1,ii
           cc  = cc  + matrixAsig(is,i,j)*matrixAA(is,j,m)*dtempdd
         enddo
        enddo
      else
        do i=1,ii
         dtempddc = dtempd * conjg(matrixAAc(is,i,m))
         do j=1,ii
           ccc = ccc + matrixAsig(is,i,j)*matrixAAc(is,j,m)*dtempddc
         enddo
        enddo
      endif
     enddo

     if(mpi_size/=1.and..not.breaknfscom) then
       write(*,*) 'I AM DONE RANK : ', mpi_rank
       if(nkpoints==1)then
          call sync_nodes_by_shared_nfs(cc,'summing_up_cc' ,mpi_rank,mpi_size)
       else
        call sync_nodes_by_shared_nfsc(ccc,'summing_up_ccc',mpi_rank,mpi_size)
       endif
     endif

     if( optimisation_parallel )Then
     if( (.not.breaknfscom.and.mpi_size==1) .or. (breaknfscom) )then
      if(pub_total_num_nodes>1)then
       if(nkpoints==1)then
         call comms_reduce( 'SUM',cc  )
       else
         call comms_reduce( 'SUM',ccc )
       endif
      endif
     endif
     endif

    else
     write(*,*) 'compute 1/2 Tr ( Sigma_an * G_an ) on GPU'
     cc=0.0;ccc=0.d0

     if(mpi_rank==1.or.breaknfscom)then
      do m=1,ii
       lambda = vectorAA(is,m) - chem
       dtempd = 1.d0 / ( -lambda ) * ( 1.d0 - DEXPc(-beta*lambda) ) / ( 1.d0 + DEXPc(-beta*lambda) )
       if(nkpoints==1)then
        tmpr(m,:)=matrixAA(is,:,m)*dtempd
       else
        tmpc(m,:)=matrixAAc(is,:,m)*dtempd
       endif
      enddo
      write(*,*) 'MATMUL CUDA ', ii
#ifdef GPU_SPEEDUP
      if(nkpoints==1)then
        call matmulcuda_r(       matrixAA(is,:,:), tmpr,tmpr2,ii,ii,ii) ! tmpr2=matrixAA(is,j,m)*tmpr(m,i)
      else
        call matmulcuda_c(conjg(matrixAAc(is,:,:)),tmpc,tmpc2,ii,ii,ii)
      endif
#else
      write(*,*) 'not supposed to get here';  stop
#endif
      write(*,*) 'done cuda'
      if(nkpoints==1)then
         cc=0.0
         do i=1,ii
          do j=1,ii
           cc = cc + matrixAsig(is,i,j)*tmpr2(j,i)
          enddo
         enddo
         write(*,*) 'CC FROM CUDA : ', cc
      else
         ccc=0.0
         do i=1,ii
          do j=1,ii
            ccc = ccc + matrixAsig(is,i,j)*tmpc2(i,j)
          enddo
         enddo

#ifdef debug
        if(pub_my_node_id==0) then
          write(*,*) 'max matrixAsig   : ', maxval(abs(matrixAsig(is,:,:))) 
          write(*,*) 'assym matrixAsig : ', maxval(abs( matrixAsig(is,:,:) - transpose( matrixAsig(is,:,:)) )) 
          write(*,*) 'Asig-Asig^+      : ', maxval(abs( matrixAsig(is,:,:) - transpose(conjg( matrixAsig(is,:,:) )) ))
          write(*,*) 'max tmpc2        : ', maxval(abs(tmpc2)) 
          write(*,*) 'assym tmpc2      : ', maxval(abs( tmpc2 - transpose( tmpc2 ) )) 
          write(*,*) 'tmpc2-tmpc2^+    : ', maxval(abs( tmpc2 - transpose(conjg( tmpc2) ) ))
          write(*,*) 'max matrixAAc    : ', maxval(abs(matrixAAc)) 
          write(*,*) 'max tmpc         : ', maxval(abs(tmpc)) 
        endif
#endif
        if(pub_my_node_id==0) write(*,*) 'CCC FROM CUDA : ', ccc
      endif
     endif
     if(mpi_size/=1.and..not.breaknfscom) then
       write(*,*) 'I AM DONE RANK : ', mpi_rank
       if(nkpoints==1)then
          call sync_nodes_by_shared_nfs(cc  ,'summing_up_cc',mpi_rank,mpi_size)
       else
         call sync_nodes_by_shared_nfsc(ccc,'summing_up_ccc',mpi_rank,mpi_size)
       endif
     endif
     write(*,*) 'done all'
    endif

#ifdef debug
    write(*,*) 'testing'
    ccback=cc;cc=0.d0;ccbackc=ccc;ccc=0.d0

    do m=1,ii
     lambda = vectorAA(is,m) - chem
     dtempd = 1.d0 / ( -lambda ) * ( 1.d0 - DEXPc(-beta*lambda) ) / ( 1.d0 + DEXPc(-beta*lambda) )
     if(nkpoints==1)then
       do i=1,ii
        dtempdd = dtempd * matrixAA(is,i,m)
        do j=1,ii
         cc = cc + matrixAsig(is,i,j)*matrixAA(is,j,m)*dtempdd
        enddo
       enddo
     else
       do i=1,ii
        dtempddc = dtempd * conjg(matrixAAc(is,i,m))
        do j=1,ii
         ccc = ccc + matrixAsig(is,i,j)*matrixAAc(is,j,m)*dtempddc
        enddo
       enddo
     endif
    enddo

    if(nkpoints==1)then
     if(abs(ccback-cc)>1.d-5)then
       write(*,*) 'ERROR                : MPI OBTAINED BY NFS-SHARED IS NOT CONSISTENT'
       write(*,*) 'DIFF OBTAINED        : ', abs(ccback-cc)  
       write(*,*) 'MAX - MIN matrixAsig : ', minval(abs(matrixAsig)),maxval(abs(matrixAsig))
       write(*,*) 'MAX - MIN matrixAA   : ', minval(abs(matrixAA(is,:,:))),maxval(abs(matrixAA(is,:,:)))
     endif
    else
     if(abs(ccbackc-ccc)>1.d-5)then
       write(*,*) 'ERROR : MPI OBTAINED BY NFS-SHARED IS NOT CONSISTENT complex'
       write(*,*) 'DIFF OBTAINED         : ', abs(ccbackc-ccc)
       write(*,*) 'MAX - MIN matrixAsig  : ', minval(abs(matrixAsig)),maxval(abs(matrixAsig))
       write(*,*) 'MAX - MIN matrixAAc   : ', minval(abs(matrixAAc(is,:,:))),maxval(abs(matrixAAc(is,:,:)))
     endif
    endif
#endif

    if(nkpoints==1)then
      if(pub_my_node_id==0) write(*,*) '<HU> Sigma_an * G_an              : ', 0.25*cc
    else
      if(pub_my_node_id==0) write(*,*) '<HU> Sigma_an * G_an              : ', 0.25*ccc
    endif
    if(nkpoints==1)then
      dmft_Energy(is) = dmft_Energy(is) + 0.25d0 * cc 
    else
      dmft_Energy(is) = dmft_Energy(is) + 0.25d0 * ccc
    endif

    if(pub_my_node_id==0) write(*,*) ' starting totE calculation , ii = ', ii

    totE=0.d0

    if( optimisation_parallel )Then
     iistart=pub_my_node_id;iistep=pub_total_num_nodes
    else
     iistart=0;iistep=1
    endif
    if(nkpoints==1)then
      do i=iistart+1,ii,iistep
       do m=1,ii
         totE = totE +      matrixContributionTail(i,m)   * real(matrixSigStat(is,m,i))*0.5   ! now it is correct ...
       enddo
      enddo
    else
      do i=iistart+1,ii,iistep
       do m=1,ii
         totE = totE + real(matrixContributionTailc(i,m)) * real(matrixSigStat(is,m,i))*0.5   ! now it is correct ...
       enddo
      enddo
    endif
    if(pub_total_num_nodes>1.and.optimisation_parallel) call comms_reduce( 'SUM',totE  )


    if(pub_my_node_id==0) write(*,*) '<HU> SIG(oo)                      : ', totE

    dmft_Energy(is) = dmft_Energy(is) + totE 

    if(pub_my_node_id==0) write(*,*) '<HU> TOTAL                        : ', dmft_Energy(is) 

    totE=0.d0
    if( optimisation_parallel )Then
     iistart=pub_my_node_id;iistep=pub_total_num_nodes
    else
     iistart=0;iistep=1
    endif
    if(nkpoints==1)then
      do i=iistart+1,ii,iistep
       do m=1,ii
        totE = totE + matrixDensKernel(i,m)*matrixH(is,m,i)
       enddo
      enddo
    else
      do i=iistart+1,ii,iistep
       do m=1,ii
        totE = totE + matrixDensKernelc(i,m) * matrixHc(is,m,i) 
       enddo
      enddo
    endif
    if(pub_total_num_nodes>1.and.optimisation_parallel) call comms_reduce( 'SUM',totE  )

    if(pub_my_node_id==0) write(*,*) 'ENERGY FROM Tr( G_dmft*H_DFT )    : ', totE

    ! IF NO ENERGY FROM DENS KERNEL CORRECTION DUE TO CORRELATIONS
    if(present(dmft_energy_cor))then 
       if(.not.pub_dmft_KS_shift)then
         E_sum_eigenvalues_DFT=0.d0; totE=0.d0; 
       endif
       total_energy=0.d0
    endif

    dmft_Energy(is) = dmft_Energy(is) + totE - E_sum_eigenvalues_DFT(is)  + total_energy/2.0 

    edc_correction=0.d0
    inquire(file='edc_total',exist=check)
    if(check)then
      open(unit=40000,file='edc_total')
      read(40000,*) edc_correction
      close(40000)
    endif
    if(pub_my_node_id==0)write(*,*) 'EDC CORRECTION (spin up/dn)', edc_correction

    if(is==1)then
     if(pub_my_node_id==0)then
       write(*,*) 'TOTAL E_LDA ENERGY /2                = ', total_energy/2.0
       write(*,*) 'TOTAL E_KS ENERGIES           SPIN 1 = ', E_sum_eigenvalues_DFT(1)
       write(*,*) 'DMFT TOTAL ENERGY             SPIN 1 = ', dmft_Energy(1)
       write(*,*) 'DMFT TOTAL ENERGY   (EDC COR) SPIN 1 = ', dmft_Energy(1)-edc_correction/2.d0
     endif
    elseif(is==2)then
     if(pub_my_node_id==0)then
        write(*,*) 'TOTAL E_LDA ENERGY /2                = ', total_energy/2.0
        write(*,*) 'TOTAL E_KS ENERGIES           SPIN 2 = ', E_sum_eigenvalues_DFT(2)
        write(*,*) 'DMFT TOTAL ENERGY             SPIN 2 = ', dmft_Energy(2)
        write(*,*) 'DMFT TOTAL ENERGY   (EDC COR) SPIN 2 = ', dmft_Energy(2)-edc_correction/2.d0
     endif
    endif

    66778 continue !SKIPPING ENERGY PART IF NMU_STEP/=NMU

    if(pub_my_node_id==0)then
        write(*,*) 'CHEM POTENTIAL : ', fermi_e(is)+pub_dmft_chem_shift
        write(*,*) 'TOTAL DENSITY  : ', vectorNmu(nmu_step)
        write(*,*) '==============================================='
     endif

#ifdef debug
   if(nkpoints==1)then
     num_denskern = NINT(sparse_num_rows(denskern(1)))
     if(num_denskern/=size(matrixDensKernel,1))then
       write(*,*) 'ERROR : dim do not match denskernel in dmft!'
       write(*,*) 'shape denskern ?      : ', num_denskern
       write(*,*) 'shape matrix dmft ?   : ', shape(matrixDensKernel)
       stop
     endif
   endif
#endif


    if(ikpoint>=lastkpoint.and.(is==num_spins.or.pub_dmft_impose_chem_spin))then
     if(abs( vectorNmu(nmu_step) - target_N )< pub_dmft_mu_diff_max .and. pub_dmft_skip_energy)then
        if(pub_my_node_id==0)then
         write(*,*) 'WARNING : TARGET REACHED - SKIPPING OTHER MU STEPS'
         write(*,*) ' target   : ', target_N
         write(*,*) ' obtained : ', vectorNmu(nmu_step) 
        endif
        vectorNmu(nmu)=vectorNmu(nmu_step)
        vectorMu(nmu)=vectorMu(nmu_step)
        nmu_step=nmu
        chem_shift=0.d0
     endif
    endif

    if(nmu_step==nmu.and.pub_dmft_kernel/=0) then
        if(ikpoint<=nkpoints)then
         contrib_energy=(dmft_Energy(is)-edc_correction/2.d0)*norm_kpoints(ikpoint)
        else
         contrib_energy=0.d0
        endif
        if(pub_dmft_splitk)then
           dmft_splitk_iter=dmft_splitk_iter+1
           if(pub_my_node_id==0)call sync_nodes_by_shared_nfs(contrib_energy,'summing_contrib_energy'//trim(adjustl(tostring(dmft_splitk_iter))),mpi_rank,mpi_size)
           call comms_barrier
        endif
        dmft_energy_cor_(is)=dmft_energy_cor_(is)+contrib_energy
        if(pub_my_node_id==0) write(*,*) 'last MU step reached, preparing kernel calculations' 
        if(nkpoints==1)then
          call sparse_convert(denskern_tmp(is),matrixDensKernel)
        else
          if(pub_dmft_nkpoints>0.and..not.pub_dmft_kpoints_kernel_gamma)then
             write(*,*) 'ERROR : if you use symmetries (other than k  inversion)'
             write(*,*) '        the density kernel will be wrong, only fine to compute projected quantities and total charge' 
             stop
          else
           if(pub_dmft_kpoints_kernel_gamma)Then
               if(ikpoint<=nkpoints)then
                if(norm_vector(kpt(ikpoint,:))>1.d-4) matrixDensKernelc=0.d0 !Not Gamma point
               endif
               if(pub_dmft_splitk)then
                  dmft_splitk_iter=dmft_splitk_iter+1
                  if(pub_my_node_id==0) call sync_nodes_by_shared_nfs_matc(matrixDensKernelc,'sum_matrixDensKernelGc'//trim(adjustl(tostring(dmft_splitk_iter))),mpi_rank,mpi_size,size(matrixDensKernelc,1))
                  call comms_barrier
                  call comms_bcast(0,matrixDensKernelc)
               endif
               call sparse_convert(denskern_tmpc,matrixDensKernelc)
               call sparse_axpy(denskern_tmp(is),denskern_tmpc,1.d0)
           else
              if(ikpoint<=nkpoints)then
               if(    abs(norm_kpoints(ikpoint)*totkpoint-1.d0)<1.d-4)then
                 matrixDensKernelc=(matrixDensKernelc                             )*norm_kpoints(ikpoint)
               elseif(abs(norm_kpoints(ikpoint)*totkpoint-2.d0)<1.d-4)then
                 matrixDensKernelc=(matrixDensKernelc+transpose(matrixDensKernelc))*norm_kpoints(ikpoint)/2.d0
               elseif(    norm_kpoints(ikpoint)*totkpoint > 2.001 .or. norm_kpoints(ikpoint)*totkpoint < 0.99  ) then
                 write(*,*) ' ERROR : k points - not supposed to happen stop ' 
                 write(*,*) ' norm kpoint is : ', norm_kpoints(ikpoint)*totkpoint 
                 stop
               endif
              endif
              if(pub_dmft_splitk)then
                 dmft_splitk_iter=dmft_splitk_iter+1
                 if(pub_my_node_id==0) call sync_nodes_by_shared_nfs_matc(matrixDensKernelc,'sum_matrixDensKernelGc'//trim(adjustl(tostring(dmft_splitk_iter))),mpi_rank,mpi_size,size(matrixDensKernelc,1))
                 call comms_barrier
                 call comms_bcast(0,matrixDensKernelc)
              endif
              call sparse_convert(denskern_tmpc,matrixDensKernelc)
              call sparse_axpy(denskern_tmp(is),denskern_tmpc,1.d0)
           endif
          endif
        endif

        if(ikpoint>=lastkpoint)then
         ttscreen=sparse_trace(denskern_tmp(is),overlap)
         if(nkpoints>1)then
          if(pub_my_node_id==0)write(*,*) 'TOTAL CHARGE ESTIMATED WITH rho_loc * S : ', ttscreen
          if(pub_my_node_id==0)write(*,*) 'TOTAL CHARGE with Sum_k rho_k S_k       : ', vectorNmu(nmu_step)
         endif
         if(pub_my_node_id==0)write(*,*) 'Trace DensKern(is) Overlap : ', ttscreen
         ttscreen=sparse_trace(denskern_tmp(is),ham(is))
         if(pub_my_node_id==0)write(*,*) 'Trace DensKern(is) HAM(is) : ', ttscreen
         if(pub_cell%num_spins==1)then
            dmft_energy_cor_(1)=dmft_energy_cor_(1)*2.d0
            dmft_energy_cor_(2)=0.d0
         endif
         if(pub_dmft_paramagnetic .and. pub_cell%num_spins==2) then
           if(pub_my_node_id==0)write(*,*) 'density_kernel copied from spin 1 to spin 2 (paramagnetic dmft)' 
           call sparse_copy(denskern_tmp(2),denskern_tmp(1))
           dmft_energy_cor_(1)=dmft_energy_cor_(1)*2.d0
           dmft_energy_cor_(2)=0.d0
         endif
        endif
    endif

    if(nmu_step<nmu.and.(is==num_spins.or.pub_dmft_impose_chem_spin).and.ikpoint>=lastkpoint) then 
      if(pub_my_node_id==0) write(*,*)  'nmu < nmu_step, computing derivatives...'

#ifdef debug
      if(nmu_step>1)then
        if(abs(vectorMu(nmu_step)-vectorMu(nmu_step-1))<1.d-8)then
          slope_mu =  0.0
        else
          slope_mu = (vectorNmu(nmu_step)-vectorNmu(nmu_step-1)) / (vectorMu(nmu_step) - vectorMu(nmu_step-1)) 
        endif
        if(pub_my_node_id==0) write(*,*)  'slope is : ', slope_mu
        if(abs(slope_mu)<50.0)then
         chem_shift = sign(1.d0,slope_mu) * ( target_N - vectorNmu(nmu_step) ) * 0.01
        else
         chem_shift = slope_mu_mixing * ( target_N - vectorNmu(nmu_step) )  / slope_mu
         if(chem_shift> slope_mu_max_step) chem_shift= slope_mu_max_step
         if(chem_shift<-slope_mu_max_step) chem_shift=-slope_mu_max_step
        endif
      endif

      if(nmu_step==1)then
       if(vectorNmu(nmu_step)>target_N)then
        chem_shift=-0.003
       else
        chem_shift= 0.003
       endif 
      endif
#else
      np=min(nmu_step,abs(pub_dmft_mu_order))  !  if min of (nmu_step,2) recover method if debug mode
      call shoot_next_mu(np,vectorMu(nmu_step-np+1:nmu_step),vectorNmu(nmu_step-np+1:nmu_step),slope_mu_max_step,slope_mu_mixing,target_N,reshoot=pub_dmft_mu_order<0,chem_shift=chem_shift)
#endif

      if(pub_my_node_id==0)then
        write(*,*) 'CHEM_SHIFT       : ', chem_shift
        write(*,*) 'NMU_STEP         : ', nmu_step
        write(*,*) 'TARGET_N         : ', target_N
        write(*,*) 'N_(mustep)       : ', vectorNmu(nmu_step) 
        write(*,*) 'chem(mustep)     : ', vectorMu(nmu_step)
        if(nmu_step>1)then
           write(*,*) 'N_(mustep-1)     : ', vectorNmu(nmu_step-1)      
           write(*,*) 'chem(mustep-1)   : ', vectorMu(nmu_step-1)
        endif
      endif

    endif

    if(nmu>1.and.ikpoint>=lastkpoint)then
      if(pub_my_node_id==0)then
       write(*,*) 'CHEM POT / CHARGE OBTAINED UP TO NOW ....'
       do i=1,nmu_step
         write(*,*) vectorMu(i),vectorNmu(i)
       enddo
       write(*,*) '.................'       
      endif 
      if(pub_my_node_id==0.and.mpi_rank==1)then
        open(unit=81119,file='chem.potential.nmu.iter'//trim(adjustl(toString(is))))
        write(81119,*) vectorMu(nmu_step)
        close(81119)
      endif
    endif

    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine cleanitall
    integer :: iiki

     if(split)then
      if(allocated(split_hub)) deallocate(split_hub,split_hubt,split_wgv,split_S,split_Sm,split_H)
     endif

#ifdef LINEAR_SCALING_INV
     if(pub_dmft_lin_scaling)then
      call internal_invert_linear_scaling(greenf,energy,overlap_cplx,ham(is),self_energy,update=.false.,realorcomplexH=.true.,finalize=.true.)
      if(allocated(lin_GSq))then
        deallocate(lin_ovSq,lin_HSq,lin_sigSq,lin_GSq,lin_vecs,lin_vec,lin_HrSq)
      endif
     endif
#endif

    if(nkpoints>1)then
      call sparse_destroy(    Ok          )
      call sparse_destroy( invOk          )
      call sparse_destroy( ovK            )
      call sparse_destroy( hamK           )
      call sparse_destroy( invK           )
      call sparse_destroy( selfK          )
      call sparse_destroy( selfK1         )
      call sparse_destroy( selfK2         )
      call sparse_destroy( selfK3         )
      call sparse_destroy( O_bar          )
      call sparse_destroy( self_energy_vK )
    endif

    if(allocated(orb_in_plane))  deallocate(orb_in_plane)
    if(allocated(atom_in_plane)) deallocate(atom_in_plane)

    if(allocated(mmuRar)) deallocate(mmuRar,wwRar) 
    if(allocated(mmuLar)) deallocate(mmuLar,wwLar)       
    if(allocated(ttt0R))  deallocate(ttt0R,maskTR,maskTRatom,mask_vectorR)
    if(allocated(H0R_))   deallocate(H0R_,H1R_,S0R_,S1R_,SSR_)
    if(allocated(ttR_))   deallocate(ttR_,tttR_,GGR_,GGRback_)
    if(allocated(ttt0L))  deallocate(ttt0L,maskTL,maskTLatom,mask_vectorL)
    if(allocated(H0L_))   deallocate(H0L_,H1L_,S0L_,S1L_,SSL_)
    if(allocated(ttL_))   deallocate(ttL_,tttL_,GGL_,GGLback_)
    if(allocated(ttLc_))  deallocate(ttLc_,tttLc_,ttRc_,tttRc_)

    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'FINISHED DMFT LOOP, MY ID : ', pub_my_node_id
    call comms_barrier
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'cleaning sparse matrices'

#ifdef debug
    if(allocated(vectorAAback)) deallocate(vectorAAback)
#endif
   
    if(allocated(vectorAA)) then
           deallocate(vectorAA,matrixAAt,tmpc,matrixSigStat, &
         & tmp1,tmpr,tmpr2,tmp2,greenfbackup,sigmabackup,matrixAsig,matrixGan)
        if(nkpoints==1)then
          deallocate(matrixDensKernel,  matrixAA ,matrixOverlap, matrixH ,matrixContributionTail  )
        else
          deallocate(matrixDensKernelc, matrixAAc,matrixOverlapc,matrixHc,matrixContributionTailc )
        endif  
    endif

    if(allocated(tmpc2))   deallocate(tmpc2)
#ifdef debug
    if(allocated(tmpc3))   deallocate(tmpc3)
#endif

    if(allocated(sigma_density_fine))      deallocate(sigma_density_fine,stat=ierr)
    if(allocated(green_density_fine))      deallocate(green_density_fine) 
    if(allocated(green_tot_fine))          deallocate(green_tot_fine)
    if(allocated(Z_density_fine))          deallocate(Z_density_fine)  
    if(allocated(green_density_fine_proj)) deallocate(green_density_fine_proj)
    do iiki=1,pub_cell%num_spins
      call sparse_destroy( ham_embed(iiki) )
    enddo
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'destroy 1'
    do i=1,pub_cell%num_spins
     if(pub_dmft_kernel/=0)then
       call sparse_destroy(     denskern_tmp(i) )
       call sparse_destroy( pur_denskern_tmp(i) )
       if(nkpoints>1)then
         if(i==1)then
           call sparse_destroy( denskern_tmpc )
         endif
       endif
     endif
     call sparse_destroy( green_density_spam_proj(i) )
     call sparse_destroy( green_density_spam(i)      )
     call sparse_destroy(     green_tot_spam(i)      )
     call sparse_destroy( sigma_density_spam(i)      )
     call sparse_destroy(     Z_density_spam(i)      )
    enddo
    if(pub_dmft_temp<0.0_DP.and.pub_dmft_optics)then
      write(*,*) 'destroy optics'
      if(allocated(temp1_real))    deallocate(temp1_real)
      if(allocated(temp2_real))    deallocate(temp2_real)
      if(allocated(temp2_complex)) deallocate(temp2_complex)
      if(allocated(nabla_square))  deallocate(nabla_square)
      do i=1,3
        call sparse_destroy(nabla(i))                       
      enddo
    endif
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'destroy 2'

    do i=1,pub_cell%num_spins
     call sparse_destroy(self_infinity(i))
    enddo

    call sparse_destroy(self_energy)
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'destroy 2'
    call sparse_destroy(self_energy_v)
    call sparse_destroy(w_self_energy_v)
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'destroy 3'

    if( inv_projectors )Then
      call sparse_destroy(w_inv_greenf_v)
      call sparse_destroy(w_inv_greenf)
      call sparse_destroy(hub_inv_overlap_tmatrix_cplx)
      if(.not.silent.and.pub_my_node_id==0) write(*,*) 'destroy 4'
      call sparse_destroy(hub_inv_overlap_tmatrix)
      call sparse_destroy(hub_inv_overlap_matrix_cplx)
      call sparse_destroy(hub_inv_overlap_matrix)
    endif 

    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'destroy 5'
    call sparse_destroy(w_greenf_v)
    call sparse_destroy(w_greenf)
    call sparse_destroy(hub_overlap_tmatrix_cplx)
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'destroy 6'
    call sparse_destroy(hub_overlap_matrix_cplx)
    call sparse_destroy(overlap_cplx)
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'destroy 7'
    call sparse_destroy(inv_greenf)
    call sparse_destroy(inv_greenf_backup)
    call sparse_destroy(greenf_backup)
    call sparse_destroy(greenf)
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'FINISHED CLEANING, MY ID : ', pub_my_node_id
    call comms_barrier
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'DUMPING CONFIRMATION IN FILE'
    open(unit=2020,file='onetep_confirmation_'//trim(adjustl(toString(mpi_rank))),position='append')
    write(2020,*) 1111
    close(2020)
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'done'

    if(pub_my_node_id==0.and.mpi_size>1.and.mpi_rank==1)then
#ifdef DMFT_SAFE_NFS
        ressys=system("sleep 4")
#endif
#ifdef debug
        ressys=system("rm ./onetep_dmft_confirmation_dmft_density_kernel* onetep_dmft_confirmation_dmft_density_count > /dev/null 2>&1 ")
#endif
    endif

    if (pub_on_root) write(stdout,'(a)')'================================&
         &================================================'

#ifdef GPU_SPEEDUP
     if(use_gpu_onlyme)then
      if(.not.silent) write(*,*) 'SHUTTING DOWN CUBLAS , node/mpirank : ', pub_my_node_id, mpi_rank
      call cublas_shutdown()
     endif
#endif
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eigenvec_fullmat(ii,a,b,vecout,val)
    use dense, only : DEM, dense_eigensolve, dense_create, dense_destroy, dense_convert,dense_get_element
    use sparse,   only: sparse_convert
    integer      :: ii,i,j
    TYPE(SPAM3)  :: a,b
    real(8)      :: vecout(:,:),val(:)
    type(DEM)    :: aa,bb,vvec
      call dense_create(aa,ii,ii)
      call dense_create(bb,ii,ii)
      call dense_create(vvec,ii,ii)
      call dense_convert(aa,a)
      call dense_convert(bb,b)
      call dense_eigensolve(ii,val,aa,bb,1,vvec)
      do i=1,ngwf_basis%num
       do j=1,ngwf_basis%num
         call dense_get_element(vecout(i,j),vvec,i,j)
       enddo
      enddo
      call dense_destroy(aa)
      call dense_destroy(bb)
      call dense_destroy(vvec)
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eigenvec_fullmatc(ii,a,b,vecout,val)
    use dense, only : DEM, dense_eigensolve, dense_create, dense_destroy, dense_convert,dense_get_element
    use sparse,   only: sparse_convert
    integer      :: ii,i,j
    TYPE(SPAM3)  :: a,b
    complex(8)   :: vecout(:,:)
    real(8)      :: val(:)
    type(DEM)    :: aa,bb,vvec
      call dense_create(aa,ii,ii,iscmplx=.true.)
      call dense_create(bb,ii,ii,iscmplx=.true.)
      call dense_create(vvec,ii,ii,iscmplx=.true.)
      call dense_convert(aa,a)
      call dense_convert(bb,b)
      call dense_eigensolve(ii,val,aa,bb,1,vvec)
      do i=1,ngwf_basis%num
       do j=1,ngwf_basis%num
         call dense_get_element(vecout(i,j),vvec,i,j)
       enddo
      enddo
      call dense_destroy(aa)
      call dense_destroy(bb)
      call dense_destroy(vvec)
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine internal_invert_full_real(mat_inv,mat,greenfunc)
      use dense, only: DEM, dense_create, dense_destroy, dense_convert,dense_invert
      implicit none
      type(SPAM3),intent(in)    :: mat
      type(SPAM3),intent(inout) :: mat_inv
      type(DEM)                 :: mat_dens
      logical                   :: greenfunc
      call dense_create(mat_dens,ngwf_basis%num,ngwf_basis%num,iscmplx=.false.)
      call dense_convert(mat_dens,mat)
      call dense_invert(mat_dens)
      call dense_convert(mat_inv,mat_dens)
      call dense_destroy(mat_dens)
      if(greenfunc) greenfbackupflag=.false.

    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_invert_full(mat_inv,mat,greenfunc)
      use dense, only: DEM, dense_create, dense_destroy, dense_convert, &
           dense_invert
      use sparse, only: sparse_num_cols,sparse_num_rows
      implicit none
      type(SPAM3),intent(in)    :: mat
      type(SPAM3),intent(inout) :: mat_inv
      type(DEM)                 :: mat_dens
      integer                   :: nrows
      logical                   :: greenfunc

      if(pub_total_num_nodes==1)then
        call internal_invert_full_slow(mat_inv,mat,greenfunc) !back to LAPACK
        return
      endif 
#ifdef INVERT_N_CUBE
        call internal_invert_full_slow(mat_inv,mat,greenfunc) !back to LAPACK
        return
#endif

      nrows = NINT(sparse_num_rows(mat))

      if(nrows/=sparse_num_cols(mat))then
        write(*,*) 'ERROR in internal_invert_full (dense matrix inversion)'
        write(*,*) ' matrix is rectangular ! '
        stop
      endif

      call dense_create(mat_dens,nrows,nrows,iscmplx=.true.)

      if(verboseall) write(*,*) 'dense convert'
      call dense_convert(mat_dens,mat)
      if(verboseall) write(*,*) 'done'

      if(verboseall) write(*,*) 'dense invert'
      call dense_invert(mat_dens)
      if(verboseall) write(*,*) 'done'

      if(greenfunc) greenfbackupflag=.false.

      call dense_convert(mat_inv,mat_dens)

      call dense_destroy(mat_dens)

    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE invert_gen_cmat(n,A)
        integer        ::  n,piv(n),INFO
        complex(8)     ::  A(:,:)
        complex(8)     ::  WORK(n)
        CALL ZGETRF(n,n,A,n,piv,INFO)
        CALL ZGETRI(n,A,n,piv,WORK,n,INFO)
        if(INFO/=0)then
          write(*,*) 'ERROR : matrix has no inverse'
        endif
      END SUBROUTINE

      SUBROUTINE invert_gen_rmat(n,A)
        integer     ::  n,piv(n),INFO
        real(8)     ::  A(:,:)
        real(8)     ::  WORK(n)
        CALL DGETRF(n,n,A,n,piv,INFO)
        CALL DGETRI(n,A,n,piv,WORK,n,INFO)
        if(INFO/=0)then
          write(*,*) 'ERROR : matrix has no inverse'
        endif
      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_invert_linear_scaling(G,energy,ov,H,sig,update,realorcomplexH,finalize)
#ifdef LINEAR_SCALING_INV
    use inversion_with_window, only : invert_window
#endif
    use sparse,   only : sparse_convert
    use wrappers, only : wrappers_invert_sym_cmatrix
    implicit none
    type(SPAM3),intent(in)    ::  ov,H,sig
    type(SPAM3),intent(inout) ::  G
    integer                   ::  nnn
    complex(8)                ::  energy,ww
    real(8)                   ::  chem
    logical                   ::  update,realorcomplexH
    integer                   ::  lin_inv_nval,lin_inv_ns
    real(8)                   ::  lin_inv_win
    real(8)                   ::  lin_inv_tail,lin_inv_centre,cutoffsig,cutoffH
    integer                   ::  nmpi,nconv
    logical,optional          ::  finalize

! COMMENT : 1) this will NOT do the special point correctly, especially Simp Eimp will not be obtained
!              the high frequency limit will be fine though, so Simp obtained from the tail will be correct
!              so eimp_from_onetep and simp_from_onetep need to be turned off in the dmft package
! COMMENT : 2) this wont work for nmu>1 or the density kernel, use nmu=1 AND dmft_kernel=0
! COMMENT : 3) the chemical potential should hence be the one obtained by the DFT - not from the chem.potential file

#ifdef LINEAR_SCALING_INV
      if(.not.pub_dmft_lin_scaling)then
        write(*,*) 'ERROR not supposed to end up in the linear scaling inversion of the green function subroutine'
        stop
      endif
      lin_inv_ns=20
      lin_inv_nval=min(pub_dmft_nval,ngwf_basis%num-1)
      lin_inv_win=pub_dmft_win
      nnn=ngwf_basis%num

      if(.not.present(finalize))then
         if(.not.allocated(lin_GSq))then
           allocate(lin_HrSq(nnn,nnn),lin_ovSq(nnn,nnn),lin_HSq(nnn,nnn),lin_sigSq(nnn,nnn),lin_GSq(nnn,nnn),lin_vecs(nnn,lin_inv_nval),lin_vec(lin_inv_nval))
         endif
         if(update)then
           call sparse_convert(lin_ovSq,ov)
           if(.not.realorcomplexH)then
             call sparse_convert(lin_HSq,H)  ! Hk
           else
             call sparse_convert(lin_HrSq,H) ! H is real
             lin_HSq=lin_HrSq
           endif
         endif
         chem=fermi_e(is)+pub_dmft_chem_shift
         lin_inv_centre=chem 
         lin_inv_tail=pub_dmft_scaling_tail
         ww=energy-chem
         nmpi=pub_dmft_scaling_nmpi
         cutoffH=pub_dmft_scaling_cutoff_h
         cutoffsig=pub_dmft_scaling_cutoff
         if(split)then
             lin_sigSq=sigmabackup
         else
             call sparse_convert(lin_sigSq,sig)
         endif
         if(update)then
             lin_HSq   = lin_HSq   - real(matrixSigStat(is,:,:))
         else
             lin_sigSq = lin_sigSq - real(matrixSigStat(is,:,:))
         endif
      endif

      if(mpi_rank/=1.and.update.and..not.breaknfscom)then
         if(.not.present(finalize)) then 
          lin_vecs=0.d0 ; lin_vec=0.d0; lin_rconv_back=0.0
         endif
         goto 3139
      endif

      if(mpi_size==1.or.pub_total_num_nodes>1.or.breaknfscom) nmpi=1

      call invert_window(  pub_my_node_id,pub_total_num_nodes,lin_HSq,ww,lin_GSq,             &
                       &   lin_inv_nval,nconv,lin_inv_ns,lin_inv_win,lin_inv_centre,.false.,  &
                       &   update,lin_sigSq,cutoffsig,                                        &
                       &   pub_dmft_scaling_meth,lin_ovSq,lin_inv_tail,lin_inv_centre,        &
                       &   mpi_rank,cutoffH,nmpi,pub_dmft_scaling_tol,                        &
                       &   pub_dmft_scaling_maxspace,pub_dmft_scaling_iter,lin_nconv_back,    &
                       &   lin_vecs,lin_vec,finalize)
      lin_rconv_back=dble(lin_nconv_back)
      3139 continue
      if(mpi_size>1.and..not.present(finalize).and.update.and..not.breaknfscom)then
        call sync_nodes_by_shared_nfs_recc(lin_vecs,'nfs_linear_inversionVECS',mpi_rank,mpi_size,nnn,lin_inv_nval,longdelay=.true.)
        call sync_nodes_by_shared_nfs_vec(lin_vec,'nfs_linear_inversionVEC',mpi_rank,mpi_size,lin_inv_nval)
        call sync_nodes_by_shared_nfs(lin_rconv_back,'nfs_linear_inversionSCAL',mpi_rank,mpi_size,longdelay=.true.)
        lin_nconv_back=NINT(lin_rconv_back)
      endif

      if(.not.present(finalize)) then
       if(.not.split)then
         call sparse_convert(G,lin_GSq)
       endif
       greenfbackupflag=.true.
       greenfbackup=lin_GSq
      endif

#endif
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_invert_full_slow(mat_inv,mat,greenfunc)
      use sparse,   only: sparse_convert
      use wrappers, only: wrappers_invert_sym_cmatrix 
      implicit none
      type(SPAM3),intent(in)                        :: mat
      type(SPAM3),intent(inout)                     :: mat_inv
      complex(kind=DP), allocatable, dimension(:,:) :: mat_square      
      logical                                       :: greenfunc

      if(verboseall) write(*,*) 'ALLOCATE MATRICES, basis : ', ngwf_basis%num

      allocate(mat_square(ngwf_basis%num, ngwf_basis%num),stat=ierr)
      call utils_alloc_check('internal_invert_full_slow','mat_square',ierr)

      call sparse_convert(mat_square,mat)
      
      if(verboseall) write(*,*) 'DIAGO NOW'
      if(nkpoints==1)then
       call wrappers_invert_sym_cmatrix(mat_square, ngwf_basis%num)
      else
       call invert_gen_cmat(ngwf_basis%num,mat_square)
      endif

      if(verboseall) write(*,*) 'done'

      call sparse_convert(mat_inv,mat_square)

      if(greenfunc)then
       greenfbackupflag=.true.
       greenfbackup=mat_square
      endif

      deallocate(mat_square,stat=ierr)
      call utils_dealloc_check('internal_invert_full_slow','mat_square',ierr)

    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_optics_calculation
  use sparse,    only:  sparse_convert
  implicit none
  integer :: i
  logical :: check_store_nabla
  real(8) :: norm_optic
    if(pub_dmft_optics)then
     write(*,*) '--initialize optical conductivity calculations--'
     if(allocated(nabla_square)) deallocate(nabla_square)
     allocate(nabla_square(ngwf_basis%num,ngwf_basis%num,3))
     if(allocated(temp1_real))    deallocate(temp1_real)
     if(allocated(temp2_real))    deallocate(temp2_real)
     if(allocated(temp2_complex)) deallocate(temp2_complex)
     if(pub_dmft_optics)then
       allocate(temp1_real(ngwf_basis%num,ngwf_basis%num))
       allocate(temp2_real(ngwf_basis%num,ngwf_basis%num))
       allocate(temp2_complex(ngwf_basis%num,ngwf_basis%num))
     endif
     do i=1,3
        call sparse_create(nabla(i),overlap,iscmplx=.false.) !BUG corrected
       !call sparse_create(nabla(i),greenf, iscmplx=.false.) 
     enddo
     write(*,*) 'build <phi_i | d/dx phi_j>'
     inquire(file='store_nabla1',exist=check_store_nabla)
    !---------------------------------------------------------------!
    if(check_store_nabla)then
      open(unit=40000,file='store_nabla1',form='unformatted')
      read(40000) nabla_square
      close(40000)
    else
     ! BUG (?) :
     ! call integrals_grad( nabla, rep%ngwfs_on_grid , ngwf_basis, rep%ngwfs_on_grid, ngwf_basis )
       if(.not.present(proj_basis))then
         write(*,*) 'ERROR : COMPUTING OPTICS BUT NLPP ARGUMENT PROJ_BASIS'
         write(*,*) '        IS NOT PRESENT, STOP'
         stop
       endif
       call optics_grad_mat_els(nabla, rep, ngwf_basis, proj_basis)
      do i=1,3
        call sparse_convert( nabla_square(:,:,i), nabla(i) )
      enddo
      if(pub_my_node_id==0 .and. mpi_rank==1)then
        open(unit=40000,file='store_nabla1',form='unformatted')
        write(40000) nabla_square
        close(40000)
        write(*,*) 'store_nabla1 file written'
      endif
    endif
    if(abs(pub_dmft_optics_x1)>1.d-5.or.abs(pub_dmft_optics_y1)>1.d-5.or.abs(pub_dmft_optics_z1)>1.d-5)then
      norm_optic=sqrt(pub_dmft_optics_x1**2+pub_dmft_optics_y1**2+pub_dmft_optics_z1**2)
      pub_dmft_optics_x1=pub_dmft_optics_x1/norm_optic
      pub_dmft_optics_y1=pub_dmft_optics_y1/norm_optic
      pub_dmft_optics_z1=pub_dmft_optics_z1/norm_optic
      nabla_square(:,:,1)=pub_dmft_optics_x1*nabla_square(:,:,1)+pub_dmft_optics_y1*nabla_square(:,:,2) +pub_dmft_optics_z1*nabla_square(:,:,3) 
      if(abs(pub_dmft_optics_x1**2+pub_dmft_optics_y1**2+pub_dmft_optics_z1**2-1.d0)>1.d-4)then
        write(*,*) 'ERROR pub_dmft_optics_x1,y1,z1 wrong, should be a unitary vector'
        stop
      endif
      if(pub_dmft_optics_i1/=1.or.pub_dmft_optics_i2/=1)then
        write(*,*) 'ERROR optical conductivity as defined along a vector only works for current correlator along same direction'
        write(*,*) ' please use pub_dmft_optics_i1=pub_dmft_optics_i2=1'
        stop
      endif
    endif
    !---------------------------------------------------------------!
     write(*,'(a,200f8.4)') 'done, testing : nabla_square(1:2,1:2) is ', nabla_square(1:2,1:2,1)
    endif
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine reset_kernel_and_energy
      if(nkpoints==1)then
       matrixDensKernel=0.d0
      else
       matrixDensKernelc=0.d0
      endif
      dmft_Energy(is)=0.d0
 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine load_k_ham_and_overlap
 use sparse, only: sparse_convert,sparse_num_cols,sparse_num_rows
 implicit none
 character(20)                         :: filek
 logical                               :: check_k
 integer                               :: iproj,ingwf,i1,i2,kk
 logical                               :: verbose
#ifdef debug
 complex(8),allocatable,dimension(:,:) :: mat1,mat2,mat3,mat4,mat5,mat6
#endif

         if(split)then
          if(.not.allocated(split_hub))then
            allocate( split_hub(NINT(sparse_num_rows(hub_overlap_matrix_cplx)), NINT(sparse_num_cols(hub_overlap_matrix_cplx))))
            allocate(split_hubt(NINT(sparse_num_rows(hub_overlap_tmatrix_cplx)),NINT(sparse_num_cols(hub_overlap_tmatrix_cplx))))
            allocate(split_wgv(NINT(sparse_num_rows(w_greenf_v)),NINT(sparse_num_cols(w_greenf_v))))
            allocate(split_S(NINT(sparse_num_rows(overlap_cplx)),NINT(sparse_num_cols(overlap_cplx))))
            allocate(split_Sm(NINT(sparse_num_rows(overlap_cplx)),NINT(sparse_num_cols(overlap_cplx))))
            allocate(split_H(NINT(sparse_num_rows(overlap_cplx)),NINT(sparse_num_cols(overlap_cplx))))
          endif
         endif
          
         if(nkpoints==1) then
           if(pub_dmft_lin_scaling) then
             call build_self_infinity
             call internal_invert_linear_scaling(greenf,energy,overlap_cplx,ham(is),self_energy,update=.true.,realorcomplexH=.true.)
           endif
           kp_file_shift=0 
           if(split)then
             call sparse_convert(split_hubt,hub_overlap_tmatrix_cplx)
             call sparse_convert(split_hub , hub_overlap_matrix_cplx)
             allocate(split_Sr(NINT(sparse_num_rows(overlap_cplx)),NINT(sparse_num_cols(overlap_cplx))))
             call sparse_convert(split_Sr,inv_overlap);split_Sm=split_Sr; 
             call sparse_convert(split_Sr,overlap);    split_S =split_Sr;
             deallocate(split_Sr)
             allocate(split_Hr(NINT(sparse_num_rows(ham(is))),NINT(sparse_num_cols(ham(is)))))
             call sparse_convert(split_Hr,ham(is));    split_H=split_Hr; 
             deallocate(split_Hr)
           endif
           return
         endif

         verbose=pub_my_node_id==0.and.mpi_rank==1

         if(.not.pub_dmft_splitk)then
           if(one_k_one_file)then
             kp_file_shift=mpi_size*(ikpoint-1)
           else
             kp_file_shift=0
           endif
         else
             kp_file_shift=0
         endif

         filek="store_hamk_"//trim(adjustl(toString( ikpoint )))//"_"//trim(adjustl(toString( is )))
         inquire(file=trim(adjustl(filek)),exist=check_k)
         if(.not.check_k)then
           write(*,*) ' ERROR : doing calculations with [x] k points : ', nkpoints
           write(*,*) ' The store file is not present : ', trim(Adjustl(filek))
           stop
         endif
         if(verbose) write(*,*) 'loading H_k from store file'
         call sparse_readc_(hamK,trim(adjustl(filek)))
                                      !-------------!
         filek="store_ovk_"//trim(adjustl(toString( ikpoint )))//"_"//trim(adjustl(toString( is )))
         inquire(file=trim(adjustl(filek)),exist=check_k)
         if(.not.check_k)then
           write(*,*) ' ERROR : doing calculations with [x] k points : ', nkpoints
           write(*,*) ' The store file is not present : ', trim(Adjustl(filek))
           stop
         endif

         if(verbose) write(*,*) 'loading Overlap_k from store file'
         call sparse_readc_(ovK,trim(adjustl(filek)))
         call sparse_copy(overlap_cplx,ovK)
         call invert_mat_sparsec(ovK,invK,1,greenfunc=.false.)

         if(pub_my_node_id==0.and.mpi_rank==1) write(*,*) 'Now update the projectors'

         filek="store_hubk_"//trim(adjustl(toString( ikpoint )))
         inquire(file=trim(adjustl(filek)),exist=check_k)
         if(.not.check_k)then
           write(*,*) ' ERROR : doing calculations with [x] k points : ', nkpoints
           write(*,*) ' The store file is not present : ', trim(Adjustl(filek))
           stop
         endif
         call sparse_readc_(hub_overlap_matrix_cplx,trim(adjustl(filek)))

         filek="store_hubTk_"//trim(adjustl(toString( ikpoint )))
         inquire(file=trim(adjustl(filek)),exist=check_k)
         if(.not.check_k)then
           write(*,*) ' ERROR : doing calculations with [x] k points : ', nkpoints
           write(*,*) ' The store file is not present : ', trim(Adjustl(filek))
           stop
         endif
         call sparse_readc_(hub_overlap_tmatrix_cplx,trim(adjustl(filek)))

         if(verbose) write(*,*) 'Computing P x {S_k^-1} x P'
         ! take care : Obar and invOk need to be in the unrotated basis
         !    Pk Ok Obar^-1 Simp Obar^-1 Ok Pk^\dagger
         !      !                          !   
         !       rotation happens here     and here .........
         !   so Ok should be in the unrotated basis
         ! - rotation is NOT done in project_green_function (but in build_green_function)
         ! - Obar was rotated when collected in files green_output, but is rotated back from sigma_output

         if(verboseall)then
           write(*,*) '============================================='
           write(*,*) 'CHECKING CONSTRUCTION OF MATRIX OVK  : '
           call check_sparsec_(ovK)
           write(*,*) 'CHECKING CONSTRUCTION OF MATRIX INVK  : '
           call check_sparsec_(invK)
           write(*,*) 'CHECKING CONSTRUCTION OF MATRIX tmatrix cplx : '
           call check_sparsec_(hub_overlap_tmatrix_cplx)
           write(*,*) 'CHECKING CONSTRUCTION OF MATRIX matrix cplx : '
           call check_sparsec_(hub_overlap_matrix_cplx)
         endif

         call project_green_function(Ok,invK,nosplit=.true.)

         if(verboseall)then
           write(*,*) '============================================='
           write(*,*) 'CHECKING CONSTRUCTION OF MATRIX Ok : '
           call check_sparsec_(Ok)
         endif

         if(verbose) write(*,*) 'Computing inverse of projection'
         call invert_mat_sparsec(Ok,invOk,2,greenfunc=.false.)
         if(verbose) write(*,*) 'done'  

         if(verboseall)then
           write(*,*) 'CHECKING CONSTRUCTION OF MATRIX INVOK : '
           call check_sparsec_(invOk)
         endif

         force_reopen_file=.true.
         kien=pub_dmft_points-6; ien=kien; call construct_self_energy; call sparse_copy(O_bar,w_self_energy_v)
         force_reopen_file=.false.
         if(verboseall)then
           write(*,*) 'CHECKING CONSTRUCTION OF MATRIX OBAR : '
           call check_sparsec_(O_bar)
         endif

         if(pub_dmft_lin_scaling) then
            call build_self_infinity
            call internal_invert_linear_scaling(greenf,energy,overlap_cplx,hamK,selfK,update=.true.,realorcomplexH=.false.)
         endif

         call update_inverse_projectors

         if(split)then
            call sparse_convert(  split_hubt  ,  hub_overlap_tmatrix_cplx  )
            call sparse_convert(  split_hub   ,  hub_overlap_matrix_cplx   )
            call sparse_convert(  split_Sm    ,  invK                      )
            call sparse_convert(  split_S     ,  ovK                       )
            call sparse_convert(  split_H     ,  hamK                      )
         endif

#ifdef debug
           write(*,*) 'DEBUG MODE - RUNS SOME TESTS  MATRIX INVERSIONS ...'
           iproj  = NINT(sparse_num_rows(hub_overlap_tmatrix_cplx))
           ingwf  = NINT(sparse_num_cols(hub_overlap_tmatrix_cplx))
           allocate(mat1(iproj,ingwf),mat2(ingwf,iproj),mat5(ingwf,ingwf),mat6(iproj,iproj))
           call sparse_convert(mat1,hub_overlap_tmatrix_cplx)
           call sparse_convert(mat2,hub_overlap_matrix_cplx)
           call sparse_convert(mat5,overlap_cplx)
           allocate(mat3(size(mat1,1),size(mat2,2)))
           allocate(mat4(size(mat2,1),size(mat1,2)))
           mat3=MATMUL(mat1,mat2)
           mat4=MATMUL(mat2,mat1)
           
          !mat1=Hub^T
          !mat2=Hub
          !proj space
           write(300+ikpoint,*) 'Hub^T Hub'
           do i=1,size(mat3,1)
              write(300+ikpoint,*) 'i -',mat3(i,i)
           enddo
          !NGWFs space
           write(200+ikpoint,*) 'Hub Hub^T'
           do i=1,size(mat4,1)
              write(200+ikpoint,*) 'i -',mat4(i,i)
           enddo

           call sparse_convert(mat1,hub_overlap_tmatrix_cplx)
           call sparse_convert(mat4,overlap_cplx)
           call sparse_convert(mat5,invK)
           write(*,*) 'MAT6'
           call sparse_convert(mat6,invOk)

           write(800+ikpoint,*) 'MAX-MIN S_k*S_k^-1 : ', &
             & minval(abs(MATMUL(mat4,mat5))),maxval(abs(MATMUL(mat4,mat5)))

           mat1=MATMUL(mat1,mat5)
           mat3=MATMUL(mat1,mat2)
           write(800+ikpoint,*) 'Hub^T S^-1 Hub = Id_proj'
           do i=1,size(mat3,1)
              write(800+ikpoint,*) 'i -',mat3(i,i)
           enddo

           mat3=MATMUL(mat6,mat3)
           write(900+ikpoint,*) 'Hub^T S^-1 Hub * invOk'
           do i=1,size(mat3,1)
              write(900+ikpoint,*) 'i -',mat3(i,i)
           enddo
           deallocate(mat1,mat2,mat3,mat4,mat5,mat6)

#endif

 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine update_inverse_projectors
 implicit none

       if(.not.inv_projectors) return

#ifndef GPU_SPEEDUP_MATMUL_R
                 call sparse_product(hub_inv_overlap_matrix_cplx,invK,hub_overlap_matrix_cplx)
#else
       if(use_gpu_onlyme)then
       call matmul_in_dense_format_r(hub_inv_overlap_matrix_cplx,invK,hub_overlap_matrix_cplx)
       else
                 call sparse_product(hub_inv_overlap_matrix_cplx,invK,hub_overlap_matrix_cplx)
       endif
#endif


#ifndef GPU_SPEEDUP_MATMUL_R
                 call sparse_product(hub_inv_overlap_tmatrix_cplx,hub_overlap_tmatrix_cplx,invK)
#else
       if(use_gpu_onlyme)then
       call matmul_in_dense_format_r(hub_inv_overlap_tmatrix_cplx,hub_overlap_tmatrix_cplx,invK)
       else
                 call sparse_product(hub_inv_overlap_tmatrix_cplx,hub_overlap_tmatrix_cplx,invK)
       endif
#endif

 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine build_self_infinity
  implicit none

      force_reopen_file=.true.

      kien=pub_dmft_points-1; ien=kien
      call construct_self_energy
      same_self_for_all=.true.; call upfold_self_energy(1); same_self_for_all=.false.

      kien=pub_dmft_points-1; ien=kien
      call embedding_read_potential

      if(nkpoints==1)then
        call sparse_convert(matrixSigStat(is,:,:),self_energy )
      else
        call sparse_convert(matrixSigStat(is,:,:),selfK )
      endif
      if(check_embed_h) then
        call sparse_convert(tmp1,ham_embed(is))
        matrixSigStat(is,:,:)=matrixSigStat(is,:,:)+tmp1
      endif

      force_reopen_file=.false.

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine prepare_for_dmft_density_kernel
  use wrappers, only : wrappers_dsygv_lt,wrappers_zhegv
  use sparse,   only : sparse_convert,sparse_copy,sparse_axpy
#ifdef GPU_SPEEDUP
  use fortran_cuda, only : eigenvector_gen_cuda_r
#endif
  implicit none
  integer       :: ii,i,j,m,iistart,iistep
  real(8)       :: chem,beta,occup,mintail,dd,ddd,tmpdiff
  type(SPAM3)   :: scratch
  complex(8)    :: dddc

      if(pub_my_node_id==0.and.mpi_size>1.and.mpi_rank==1) then
#ifdef DMFT_SAFE_NFS
          ressys=system("sleep 4")
#endif
#ifdef debug
          ressys=system("rm ./onetep_dmft_confirmation_dmft_density_kernel* onetep_dmft_confirmation_dmft_density_count > /dev/null 2>&1 ")
#endif
      endif

      if(pub_my_node_id==0) write(*,*) 'START PREPARING DMFT DENSITY KERNEL'

      beta=1.d0/pub_dmft_temp
      chem=fermi_e(is)+pub_dmft_chem_shift
      ii=ngwf_basis%num

      if(nkpoints==1)then
       if(nmu_step>1) goto 3030
      else
       !...nothing...!
      endif

      force_reopen_file=.true.

      kien=pub_dmft_points-1; ien=kien
      call construct_self_energy
      same_self_for_all=.true.; call upfold_self_energy(1); same_self_for_all=.false.

      !EMBEDDING POTENTIALS 
      kien=pub_dmft_points-1; ien=kien
      call embedding_read_potential
      !END EMBEDDING

      if(present(dmft_self)) call sparse_copy(dmft_self(is),self_energy)

      if(nkpoints==1)then
       call sparse_convert(matrixSigStat(is,:,:),self_energy )
       call sparse_copy(self_infinity(is),self_energy)
      else
       call sparse_convert(matrixSigStat(is,:,:),selfK )
       call sparse_copy(self_infinity(is),selfK)
       call sparse_realpartc_(self_infinity(is))
      endif
      if(nkpoints==1)then
       call sparse_axpy(self_infinity(is),ham(is),1.d0)
      else
       call sparse_axpy(self_infinity(is),hamK   ,1.d0)
      endif

      !EMBEDDING POTENTIALS
      if(check_embed_h) then
        call sparse_axpy(self_infinity(is),ham_embed(is),1.d0)
        call sparse_convert(tmp1,ham_embed(is))
        matrixSigStat(is,:,:)=matrixSigStat(is,:,:)+tmp1
      endif
      !END EMBEDDING
      if(verboseall.and.pub_my_node_id==0) write(*,*) ' MAX SELF ENERGY STATIC (max,min): ', maxval(real(matrixSigStat(is,:,:))),minval(real(matrixSigStat(is,:,:)))


      if(present(dmft_z))then
         kien=1; ien=kien
         call sparse_create(scratch,overlap,iscmplx=.false.)
         call construct_self_energy
         call upfold_self_energy
         call sparse_scale(scratch,0.d0)
        ! -Z
         call sparse_axpy(scratch,self_energy,CMPLX(0.0,-1.0,8))
         call sparse_scale(scratch,-1.d0/(pi/beta))
        ! S-Z
         call sparse_axpy(scratch,overlap,1.d0)
         call internal_invert_full_real(dmft_z(is),scratch,greenfunc=.false.)
        ! S (S-Z)^-1 
         call sparse_product(scratch,overlap,dmft_z(is)) 
         call sparse_copy(dmft_z(is),scratch) 
         call sparse_destroy(scratch)
      endif

      kien=pub_dmft_points-nspecial_frequ-1; ien=kien
      call construct_self_energy
      same_self_for_all=.true.; call upfold_self_energy(2); same_self_for_all=.false.

      if(nkpoints==1)then
       call sparse_convert( matrixAsig(is,:,:) ,self_energy )
      else
       call sparse_convert( matrixAsig(is,:,:) ,selfK       )
      endif

      !EMBEDDING POTENTIALS
      if(check_embed_h)then
        call embedding_read_potential
        call sparse_convert(tmp2,ham_embed(is))
        matrixAsig(is,:,:)=matrixAsig(is,:,:)+tmp2
      endif
      !END EMBEDDING

      force_reopen_file=.false.
 
      if(pub_my_node_id==0) write(*,*) 'building SELF ENERGY TAIL'

      kkk1=1;kkk2=1;mintail=0.0
      do i=1,ii
       do j=1,ii
          if(abs(aimag(matrixAsig(is,i,j)))>mintail.and.i==j) then
                kkk1=i
                kkk2=j
                mintail=abs(aimag(matrixAsig(is,i,j)))
          endif
          matrixAsig(is,i,j) = -aimag(matrixAsig(is,i,j)) * ( en_start+real(ien-1,kind=DP)*en_step )
       enddo
      enddo
      if(verboseall.and.pub_my_node_id==0) write(*,*) 'MAX SELF ENERGY : ',mintail,kkk1,kkk2

      if(nkpoints==1)then
        call sparse_convert( matrixH(is,:,:), ham(is) )  
        call sparse_convert( matrixOverlap  , overlap )
      else
        if(pub_my_node_id==0.and.mpi_rank==1) write(*,*) 'loading hamK'
        call sparse_convert( matrixHc(is,:,:), hamK   )
        if(pub_my_node_id==0.and.mpi_rank==1) write(*,*) 'loading ovK'
        call sparse_convert( matrixOverlapc  , ovK    )
      endif

      !H v = lambda S v
      !S = L L^T
      !L^-1 H L^-T  x (L^T v) = lambda (L^T v)
      !   C = L^-1 H L^-T  y=L^T v
      ! C y = lambda y
      ! v   = L^-T y 
      ! v1^T . v2 = y1^T L^-1 L^-1^T y2

       if(pub_my_node_id==0) then
           write(*,*) '-START TO DIAGONALIZE HAMILTONIAN-'
       endif
       if(         .not.  full_prec_dmft_energy &
          &  .and. .not. (use_gpu_onlyme.and.use_gpu_eigenvectors) &
          &  .and. .not. (use_gpu_eigenvectors.and.use_gpu_some_only) )then
           if(nkpoints==1)then
             call eigenvec_fullmat(ii,self_infinity(is),overlap,matrixAA(is,:,:),vectorAA(is,:))
           else
             call eigenvec_fullmatc(ii,self_infinity(is),overlap_cplx,matrixAAc(is,:,:),vectorAA(is,:))
           endif
       else
           if( mpi_rank/=1 .and. use_gpu_some_only .and..not.breaknfscom) then
              if(nkpoints==1)then; matrixAA(is,:,:)=0.d0;else;matrixAAc(is,:,:)=0.d0;endif; 
              vectorAA(is,:)=0.d0; 
              goto 19871
           endif

           if(nkpoints==1)then
               matrixAA(is,:,:) =matrixH (is,:,:)+real(matrixSigStat(is,:,:))
           else
               matrixAAc(is,:,:)=matrixHc(is,:,:)+real(matrixSigStat(is,:,:))
           endif

           if(pub_my_node_id==0) then
                write(*,*) '-USE DSYGV-SIZE=',ii
           endif 
           if(.not.(use_gpu_onlyme.and.use_gpu_eigenvectors))then
               if(nkpoints==1)then
                 call wrappers_dsygv_lt(eigenvectors=matrixAA(is,:,:),eigenvalues=vectorAA(is,:),overlap_square=matrixOverlap,num=ii)
               else
                 call wrappers_zhegv(eigenvectors=matrixAAc(is,:,:),eigenvalues=vectorAA(is,:),overlap_square=matrixOverlapc,num=ii)
               endif

           else
#ifdef GPU_SPEEDUP
            if(nkpoints==1)then
              call eigenvector_gen_cuda_r(   ii,matrixAA(is,:,:),matrixOverlap  ,vectorAA(is,:),tmpr,.false.)
            else
              call eigenvector_cuda_type_c(1,ii,matrixAAc(is,:,:),matrixOverlapc,vectorAA(is,:),tmpc,.false.)
#ifdef debug
              write(*,*) 'CHECKING EIGENVECTOR - VALUES CUDA'
              write(*,*) '--> running on CPU'
              call sparse_convert( matrixOverlapc  , ovK    )
              vectorAAback(is,:)=vectorAA(is,:)
              call wrappers_zhegv(eigenvectors=matrixAAc(is,:,:),eigenvalues=vectorAA(is,:),overlap_square=matrixOverlapc,num=ii)
              write(*,*) 'WARNING : CPU MIN-MAX AB VALUE EIGENVECTORS : ', minval(abs(matrixAAc(is,:,:))),maxval(abs(matrixAAc(is,:,:)))
              write(*,*) 'WARNING : GPU MIN-MAX AB VALUE EIGENVECTORS : ', minval(abs(tmpc)),maxval(abs(tmpc))
              write(*,*) 'WARNING : CPU MIN-MAX RE VALUE EIGENVECTORS : ', minval(real(matrixAAc(is,:,:))),maxval(real(matrixAAc(is,:,:)))
              write(*,*) 'WARNING : GPU MIN-MAX RE VALUE EIGENVECTORS : ', minval(real(tmpc)),maxval(real(tmpc))
              write(*,*) 'WARNING : CPU MIN-MAX IM VALUE EIGENVECTORS : ', minval(aimag(matrixAAc(is,:,:))),maxval(aimag(matrixAAc(is,:,:)))
              write(*,*) 'WARNING : GPU MIN-MAX IM VALUE EIGENVECTORS : ', minval(aimag(tmpc)),maxval(aimag(tmpc))
              write(*,*) 'WARNING :  continue with CPU eigenvectors ... safer...'
              tmpc=matrixAAc(is,:,:) 
              tmpdiff=maxval(abs( vectorAAback(is,:)-vectorAA(is,:) ) )
              write(*,*) 'DIFFERENCE IN EIGENVALUES : ', tmpdiff
              if(tmpdiff>1.d-5)then
                  write(*,*) 'ERROR GPU EIGENVALUES WRONG, DIFF FROM CPU = ', tmpdiff
              endif
#endif
            endif
#else
            write(*,*) 'error gpu 2';stop
#endif
            if(nkpoints==1)then
              matrixAA (is,:,:)=tmpr
            else
              matrixAAc(is,:,:)=tmpc
            endif

           endif
           19871 continue
           if(use_gpu_some_only.and..not.breaknfscom)then
            if(mpi_size>1) then
             write(*,*) 'I AM DONE RANK : ', mpi_rank
              if(nkpoints==1)then
                call sync_nodes_by_shared_nfs_mat(  matrixAA(is,:,:),'sumup_AA'    ,mpi_rank,mpi_size,ii)
              else
                call sync_nodes_by_shared_nfs_matc(matrixAAc(is,:,:),'sumup_AAc'   ,mpi_rank,mpi_size,ii)
              endif
                call sync_nodes_by_shared_nfs_vec(vectorAA(is,:),'sumup_vectorAA',mpi_rank,mpi_size,ii)
            endif
           endif
       endif
#ifdef debug
       if(pub_my_node_id==0) then
              write(*,*) 'TEST EIGENVALUES : ', vectorAA(is,1:10)
           if(nkpoints==1)then
              write(*,*) 'TEST EIGENVECTORS: ', matrixAA(is,1,1:10)
              write(*,*) 'TRANSPOSE        : ', matrixAA(is,1:10,1)
              write(*,*) 'MATRIX H         : ', matrixH(is,1,1:10)
              write(*,*) 'MIN-MAX EIGEN    : ', minval(vectorAA(is,:)),maxval(vectorAA(is,:))
              write(*,*) 'MIN-MAX MATRIX   : ', minval(matrixAA(is,:,:)),maxval(matrixAA(is,:,:))
              write(*,*) 'MIN-MAX COLUMN 1 : ', minval(matrixAA(is,:,1)),maxval(matrixAA(is,:,1))
              write(*,*) 'MIN-MAX LINE   1 : ', minval(matrixAA(is,1,:)),maxval(matrixAA(is,1,:))
           else
              write(*,*) 'TEST EIGENVECTORS: ', matrixAAc(is,1,1:10)
              write(*,*) 'TRANSPOSE        : ', matrixAAc(is,1:10,1)
              write(*,*) 'MATRIX H         : ', matrixHc(is,1,1:10)
              write(*,*) 'MIN-MAX EIGEN    : ', minval(vectorAA(is,:)),maxval(vectorAA(is,:))
              write(*,*) 'MIN-MAX MATRIX   : ', minval(abs(matrixAAc(is,:,:))),maxval(abs(matrixAAc(is,:,:)))
              write(*,*) 'MIN-MAX COLUMN 1 : ', minval(abs(matrixAAc(is,:,1))),maxval(abs(matrixAAc(is,:,1)))
              write(*,*) 'MIN-MAX LINE   1 : ', minval(abs(matrixAAc(is,1,:))),maxval(abs(matrixAAc(is,1,:)))
           endif
           write(*,*) '-END OF DIAGONALIZATION HAMILTONIAN-'
       endif
#endif

       if(nkpoints==1)then
         call sparse_convert( matrixOverlap  , overlap )
       else
#ifdef debug
         write(*,*) 'WARNING : copying ovK to overlap_cplx - it should not affect the result, just to make sure'
         write(*,*) '          that overlap_cplx was not affected up to this point'
         call sparse_copy(overlap_cplx,ovK)
#endif
         call sparse_convert( matrixOverlapc , overlap_cplx )
       endif

     !   A = H_LDA + Sig(oo)
     !  (H + Sig) x A_m = lambda_m S x A_m
     !          G_ab(iw)  =       A_m(a) A_m(b) /     (  iw+mu  - lambda_m  )
     !  Sum_iw (G_ab(iw)) = Sum_m A_m(a) A_m(b) / ( 1 + exp ( beta*(lambda_m-mu) ) )

      3030 continue

      if(pub_my_node_id==0) write(*,*) 'BUILDING TAIL OF GF'


     ! GF tail traced analytically 

      if(nkpoints==1)then
        matrixContributionTail=0.d0
      else
       matrixContributionTailc=0.d0
      endif

      if(.not.breaknfscom)then
        iistart=(mpi_rank-1); iistep=mpi_size
      else
        iistart=0; iistep=1
      endif

      if( optimisation_parallel )Then
      if( (.not.breaknfscom.and.mpi_size==1) .or. (breaknfscom) )then
        if(iistart>0)then; write(*,*) 'ERROR parallelization conflict';stop; endif
        iistart=pub_my_node_id; iistep=pub_total_num_nodes
      endif
      endif

      do m=iistart+1,ii,iistep
       dd=1.d0/(1.d0+dexpc(beta*(vectorAA(is,m)-chem)))
       do i=1,ii

       if(nkpoints==1)then 
          ddd=matrixAA(is,i,m)*dd
          do j=i,ii
            ! Sum_iw (G_ab(iw)) = Sum_m A_m(a) A_m(b) / ( 1 + exp ( beta*(lambda_m-mu) ) ) !
              matrixContributionTail(i,j) = matrixContributionTail(i,j) + ddd*matrixAA(is,j,m)
          enddo
       else
          dddc=matrixAAc(is,i,m)*dd 
          do j=i,ii
            ! Sum_iw (G_ab(iw)) = Sum_m A_m(a) A_m(b)* / ( 1 + exp ( beta*(lambda_m-mu) ) ) !
              matrixContributionTailc(i,j) = matrixContributionTailc(i,j) + dddc*conjg(matrixAAc(is,j,m)) 
          enddo
       endif

       enddo
      enddo
 
     if( optimisation_parallel )Then
     if( (.not.breaknfscom.and.mpi_size==1) .or. (breaknfscom) )then
     if(pub_total_num_nodes>1)then
      if(nkpoints==1)then
       call comms_reduce( 'SUM',matrixContributionTail  )
      else
       call comms_reduce( 'SUM',matrixContributionTailc )
      endif
     endif
     endif
     endif

     if(nkpoints==1)then
      do i=1,ii
       do j=i+1,ii
        matrixContributionTail(j,i) = matrixContributionTail(i,j)
       enddo
      enddo
     else
      do i=1,ii
       do j=i+1,ii
        matrixContributionTailc(j,i) = conjg(matrixContributionTailc(i,j))
       enddo
      enddo
     endif

     if(nkpoints==1)then
      if(mpi_size>1.and..not.breaknfscom) then
        write(*,*) 'I AM DONE RANK : ', mpi_rank
        call sync_nodes_by_shared_nfs_mat(matrixContributionTail,'sumup_tail',mpi_rank,mpi_size,ii)
      endif
     else
      if(mpi_size>1.and..not.breaknfscom) then
        write(*,*) 'I AM DONE RANK : ', mpi_rank
        call sync_nodes_by_shared_nfs_matc(matrixContributionTailc,'sumup_tailc',mpi_rank,mpi_size,ii)
      endif
     endif


#ifdef debug
  if(nkpoints==1)then
      tmpr=matrixContributionTail
      matrixContributionTail=0.d0
      do m=1,ii
       dd=1.d0/(1.d0+dexpc(beta*(vectorAA(is,m)-chem)))
       do i=1,ii
        ddd=matrixAA(is,i,m)*dd
        do j=1,ii
          ! Sum_iw (G_ab(iw)) = Sum_m A_m(a) A_m(b) / ( 1 + exp ( beta*(lambda_m-mu) ) ) !
            matrixContributionTail(i,j) = matrixContributionTail(i,j) + ddd*matrixAA(is,j,m)
        enddo
       enddo
      enddo
      if(maxval(abs(matrixContributionTail-tmpr))>1.d-5)then
        write(*,*) 'ERROR something wrong with nfs-share sync between nodes'
        write(*,*) ' diff = ', maxval(abs(matrixContributionTail-tmpr)) 
      endif
  else
      tmpc=matrixContributionTailc
      matrixContributionTailc=0.d0
      do m=1,ii
       dd=1.d0/(1.d0+dexpc(beta*(vectorAA(is,m)-chem)))
       do i=1,ii
        dddc=matrixAAc(is,i,m)*dd
        do j=1,ii
          ! Sum_iw (G_ab(iw)) = Sum_m A_m(a) A_m(b) / ( 1 + exp ( beta*(lambda_m-mu) ) ) !
            matrixContributionTailc(i,j) = matrixContributionTailc(i,j) + dddc*conjg(matrixAAc(is,j,m))
        enddo
       enddo
      enddo
      if(maxval(abs(matrixContributionTailc-tmpc))>1.d-5)then
        write(*,*) 'ERROR something wrong with nfs-share sync between nodes - complex case'
        write(*,*) 'contribution Tail'
        write(*,*) 'diff= ', maxval(abs(matrixContributionTailc-tmpc))
      endif
  endif
#endif

      if(pub_my_node_id==0) write(*,*) 'SETTING UP FERMI OCCUPATIONS OF THE TAIL'

      occup=0.d0
      do m=1,ii
         occup = occup +  1.d0 / (1.d0 + dexpc( beta*(vectorAA(is,m) - chem) ) )
      enddo

      if(.not.silent.and.pub_my_node_id==0) write(*,*) 'TOTAL OCCUPATION (ASSUMED ORTHOGONAL BASIS) : ', occup

      if(pub_my_node_id==0) write(*,*) 'LEAVING THE PREPARATION FOR DMFT KERNEL'

     if(nkpoints==1)then
       matrixAAt(is,:,:) =       transpose(matrixAA(is,:,:))
     else
       matrixAAt(is,:,:) = conjg(transpose(matrixAAc(is,:,:)))
     endif

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compute_dmft_energy
  use sparse,   only : sparse_convert
#ifdef GPU_SPEEDUP
  use fortran_cuda, only : matmulcuda_c,matmulcuda_c_singleprec_bypass
#endif
  implicit none
  complex(8)           :: HU
  integer              :: i,j,k,l
  integer              :: ii,m,INFO
  real(8)              :: chem,beta
  complex(8)           :: cc
  complex(8),parameter :: imi=CMPLX(0.d0,1.d0,KIND=8)
  logical              :: check
  complex(8)           :: tempd,tempdd

      write(*,*) '-------------------------------------------------'
      write(*,*) 'START COMPUTING FREQUENCY CONTRIB. TO DMFT ENERGY'
      write(*,*) '-------------------------------------------------'

      beta=1.d0/pub_dmft_temp
      chem=fermi_e(is)+pub_dmft_chem_shift
      ii=ngwf_basis%num

      if(.not.split)Then
         if(nkpoints==1)then
          call sparse_convert( tmp1, self_energy )
         else
          call sparse_convert( tmp1, selfK       )
         endif
      else
          tmp1=sigmabackup
      endif 

      if(check_embed_h) then
        call sparse_convert(tmp2,ham_embed(is))
        tmp1=tmp1+tmp2
      endif

      if(.not.greenfbackupflag.and..not.split) then
          call sparse_convert( tmp2,   greenf    ) 
      else
          tmp2=greenfbackup
      endif

#ifdef debug
      if(nkpoints>1)then
         write(*,*) 'COMPUTING G(-iw)'
         tmpc=energy*matrixOverlapc-matrixHc(is,:,:)-tmp1
         call invert_matc(size(tmpc,1),tmpc) 
         write(*,*) ' iw + mu : ', energy
         write(*,*) '-iw + mu : ', conjg(energy)
         write(*,*) 'MAX DIFFERENCE BETWEEN G(iw) AND G(iw) : '
         write(*,*) maxval(abs( tmpc - tmp2 )) 
         tmpc=energy*matrixOverlapc-matrixHc(is,:,:)
         call invert_matc(size(tmpc,1),tmpc)
         tmpc2=tmpc
         tmpc=conjg(energy)*matrixOverlapc-matrixHc(is,:,:)
         call invert_matc(size(tmpc,1),tmpc)
         write(*,*) 'MAX DIFFERENCE BETWEEN G(iw) AND G(-iw)^+ : '
         write(*,*) maxval(abs( transpose(conjg(tmpc)) - tmpc2 ))
         tmpc=conjg(energy)*conjg(matrixOverlapc)-conjg(matrixHc(is,:,:))
         call invert_matc(size(tmpc,1),tmpc)
         write(*,*) 'MAX DIFFERENCE BETWEEN G(iw) AND G(-k,-iw)* : '
         write(*,*) maxval(abs( conjg(tmpc) - tmpc2 ))
      endif
#endif

      if(pub_my_node_id==0.and.mpi_rank==1) write(123456781,*) aimag(energy),real(tmp2(kkk1,kkk2)),aimag(tmp2(kkk1,kkk2))     ! G

      write(*,*) 'COMPUTE G_ANALYTIC TERM'
      if(.not.use_gpu_onlyme)then
        write(*,*) 'G tail on CPU'
        if(nkpoints==1)then
         do m=1,ii
           tempd=1.d0/( energy - vectorAA(is,m))
           tmpc(:,m)=matrixAA(is,:,m)*tempd
         enddo
        else
         do m=1,ii
           tempd=1.d0/( energy - vectorAA(is,m))
           tmpc(:,m)=matrixAAc(is,:,m)*tempd 
         enddo
        endif
        matrixGan = MATMUL( tmpc , matrixAAt(is,:,:) ) 
#ifdef debug
        tmpc2=tmp2
#endif
        tmp2=tmp2-matrixGan
      else
        write(*,*) 'Ganalytic on GPU'
        if(nkpoints==1)then
         do m=1,ii
           tempd=1.d0/( energy - vectorAA(is,m))
           tmpc(:,m)=matrixAA(is,:,m)*tempd
         enddo
        else
         do m=1,ii
           tempd=1.d0/( energy - vectorAA(is,m))
           tmpc(:,m)=matrixAAc(is,:,m)*tempd 
         enddo
        endif
#ifdef GPU_SPEEDUP
        if( doubleprec(energy,ien) ) then
                             call matmulcuda_c(tmpc,matrixAAt(is,:,:),matrixGan,ii,ii,ii)
        else
           call matmulcuda_c_singleprec_bypass(tmpc,matrixAAt(is,:,:),matrixGan,ii,ii,ii)
        endif
#else
        write(*,*) 'not supposed to get here';  stop
#endif
        write(*,*) 'done'
#ifdef debug
        tmpc2=tmp2
#endif
        tmp2=tmp2-matrixGan
       endif

#ifdef debug
        write(*,*) 'checking'
        tmpc3=tmp2
        tmpc=matrixGan
        tmp2=tmpc2
        matrixGan=0.d0
        do m=1,ii
         tempd=1.d0/( energy - vectorAA(is,m))
         if(nkpoints==1)then
          do i=1,ii
           tempdd= tempd * matrixAA(is,i,m)
           do j=1,ii
            cc = matrixAA(is,j,m) * tempdd
            matrixGan(i,j) = matrixGan(i,j) + cc
            tmp2(i,j)      = tmp2(i,j)      - cc
           enddo
          enddo
         else
          do i=1,ii
           tempdd = tempd * matrixAAc(is,i,m)
           do j=1,ii
            cc = conjg(matrixAAc(is,j,m)) * tempdd
            matrixGan(i,j) = matrixGan(i,j) + cc
            tmp2(i,j)      = tmp2(i,j)      - cc
           enddo
          enddo
         endif
        enddo

        if(maxval(abs( tmp2 - tmpc3 )) >1.d-5)then
         write(*,*) 'ERROR TMP2 in GPU OPTIMIZATION'
         write(*,*) ' DIFF       : ', maxval(abs( tmp2 - tmpc3 ))
         write(*,*) ' doubleprec : ', doubleprec(energy,ien)
        endif 
        if(maxval(abs( matrixGan - tmpc ))>1.d-5)then
         write(*,*) 'ERROR GPU optimization - matrixGan'
         write(*,*) 'DIFF = ',maxval(abs( matrixGan - tmpc ))
        else
         write(*,*) 'DIFF GAN IS : ', maxval(abs( matrixGan - tmpc )) 
        endif
#endif

      write(*,*) 'SECOND PART'

      if(pub_my_node_id==0.and.mpi_rank==1) write(123456782,*) aimag(energy),real(tmp2(kkk1,kkk2)),aimag(tmp2(kkk1,kkk2))     ! G_num

      HU=0.d0

     !--------------------!
     !only compute Energy once the chemical potential is converged
     if(nmu_step==nmu)then
     ! Sigma * G_num !
      do i=1,ii
       do m=1,ii
        !for positive and negative matsubara frequencies at the same time
         HU = HU + 2.d0*real( tmp1(i,m)*tmp2(m,i) )                                     ! Tr ( Sigma * G_num )
       enddo
      enddo

      tmp1=tmp1-real(matrixSigStat(is,:,:))                                             ! Sigma-Sig(oo)
      if(pub_my_node_id==0.and.mpi_rank==1)write(123456783,*) aimag(energy),real(tmp1(kkk1,kkk2)),aimag(tmp1(kkk1,kkk2))     !

      tmp1=tmp1-matrixAsig(is,:,:)/(imi*myfrequ(energy))                                ! Sigma_num
      if(pub_my_node_id==0.and.mpi_rank==1)write(123456784,*) aimag(energy),real(tmp1(kkk1,kkk2)),aimag(tmp1(kkk1,kkk2))     ! Sigma_num

    ! Sigma_num * G_an
      do i=1,ii
       do m=1,ii
       !for positive and negative matsubara frequencies at the same time
        HU = HU + 2.d0*real( tmp1(i,m)*matrixGan(m,i) )                                 ! Tr ( Sig_num * G_an )
       enddo
      enddo
     endif
     !--------------------!

      dmft_Energy(is) = dmft_Energy(is) + HU*0.5*pub_dmft_temp

     ! conjg transpose tmp2 corresponds to negative matsubara frequencies !

     if(nkpoints==1)then
       matrixDensKernel  = matrixDensKernel  + tmp2 + conjg(transpose(tmp2))
     else
       matrixDensKernelc = matrixDensKernelc + tmp2 + conjg(transpose(tmp2))
     endif

      write(*,*) '-------------------------------------------------'
      write(*,*) '                  DONE                           ' 
      write(*,*) '-------------------------------------------------'

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dump_data_optics(ien,is,greenf,energy,fermi_e,chem_shift)
  use sparse,    only:  sparse_convert
#ifdef GPU_SPEEDUP
  use fortran_cuda, only : matmul_my_cuda,matmulcuda_r_cublas
#endif
  implicit none
  integer                                       :: ien,is,i,j,k,nn,sizbuffer,kkkk
  type(SPAM3),intent(in)                        :: greenf
  complex(8)                                    :: energy
  complex(8),parameter                          :: imi=CMPLX(0.d0,1.d0,KIND=8)
  real(8)                                       :: fermi_e,chem_shift

  if(pub_my_node_id==0.and..not.silent) write(*,*) 'computing optics? ' ,pub_dmft_temp<0.0_DP 
  if(pub_dmft_temp>0.0_DP.or..not.pub_dmft_optics) return


  !-----------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------!
  if( abs( real(energy) - (fermi_e+chem_shift) ) < pub_dmft_optics_window .and. ien<pub_dmft_points-nspecial_frequ+1 ) then

  if(any(abs(norm_kpoints)<1.d-3))then
   write(*,*) 'ERROR : optics wont work if some k points are excluded - some files would be missing critical stop'
   stop
  endif
  nn=ngwf_basis%num
  write(*,*) ' for frequ / energy : ', ien , energy
  write(*,*) 'convert green SPAM to dense'
  if(.not.split)then
   call sparse_convert(temp2_complex,greenf)
  else
   temp2_complex=greenfbackup
  endif
  temp1_real = (temp2_complex-conjg(transpose(temp2_complex))) / (2.d0*pi*imi) 
 !temp1_real =  aimag(temp2_complex)/pi !BUG : faster for symmetric GF (Gamma=0 only)
  write(*,*) ' G*NABMA_x '
#ifdef GPU_SPEEDUP
  if(use_gpu_onlyme)then
   if(pub_my_node_id==0.or.split)then
    CALL matmulcuda_r(temp1_real,nabla_square(:,:,pub_dmft_optics_i1),temp2_real,nn,nn,nn)
   endif
   if(pub_total_num_nodes>1.and..not.split) call comms_bcast(0,temp2_real)
  else
    temp2_real=MATMUL(temp1_real,nabla_square(:,:,pub_dmft_optics_i1))
  endif
#else
    temp2_real=MATMUL(temp1_real,nabla_square(:,:,pub_dmft_optics_i1))
#endif

  if(.not.split)then
  if(pub_my_node_id==0)then
  if((ien==ffirst.and.one_k_one_file).or.(.not.one_k_one_file.and.ien==ffirst.and.ikpoint==1))then
    open(unit=30030+mpi_rank*pub_total_num_nodes+10200*is,file='optics_data_spin_'//TRIM(ADJUSTL(toString(is)))//&
        & '_k_'//TRIM(ADJUSTL(toString(ikpoint)))//'_rank_'//TRIM(ADJUSTL(toString(mpi_rank))),form='unformatted') 
  else
    open(unit=30030+mpi_rank*pub_total_num_nodes+10200*is,file='optics_data_spin_'//TRIM(ADJUSTL(toString(is)))//&
        & '_k_'//TRIM(ADJUSTL(toString(ikpoint)))//'_rank_'//TRIM(ADJUSTL(toString(mpi_rank))),form='unformatted',position='append') 
  endif
  endif
  endif

  !======================================!

   if(.not.split)then
    if(pub_my_node_id==0)then
      if(pub_dmft_optics_i1==pub_dmft_optics_i2)then
        write(30030+mpi_rank*pub_total_num_nodes+10200*is) nkpoints,norm_kpoints(ikpoint),1,volume_cell_angstrom(),pub_dmft_temp,pub_dmft_points,mpi_size,nn,ien,energy,fermi_e,chem_shift,temp2_real
      else
        write(30030+mpi_rank*pub_total_num_nodes+10200*is) nkpoints,norm_kpoints(ikpoint),2,volume_cell_angstrom(),pub_dmft_temp,pub_dmft_points,mpi_size,nn,ien,energy,fermi_e,chem_shift,temp2_real
      endif
    endif
   else
     do kkkk=1,pub_total_num_nodes
      if(kkkk==pub_my_node_id+1)then
         if(ien==ffirst.and.kkkk==1.and.  (  (.not.pub_dmft_splitk.and.(ikpoint==1.or.one_k_one_file))  .or.  (pub_dmft_splitk.and.ikpoint==splitk_first)  )    )then
           open(unit=30030+mpi_rank*pub_total_num_nodes+10200*is,file='optics_data_spin_'//TRIM(ADJUSTL(toString(is)))//'_k_'//TRIM(ADJUSTL(toString(ikpoint)))//'_rank_'//TRIM(ADJUSTL(toString(mpi_rank))),form='unformatted',status='replace') 
         else
           open(unit=30030+mpi_rank*pub_total_num_nodes+10200*is,file='optics_data_spin_'//TRIM(ADJUSTL(toString(is)))//'_k_'//TRIM(ADJUSTL(toString(ikpoint)))//'_rank_'//TRIM(ADJUSTL(toString(mpi_rank))),form='unformatted',position='append',status='old') 
         endif
         if(pub_dmft_optics_i1==pub_dmft_optics_i2)then
          write(30030+mpi_rank*pub_total_num_nodes+10200*is) nkpoints,norm_kpoints(ikpoint),1,volume_cell_angstrom(),pub_dmft_temp,pub_dmft_points,mpi_size,nn,ien,energy,fermi_e,chem_shift,temp2_real
         else
          write(30030+mpi_rank*pub_total_num_nodes+10200*is) nkpoints,norm_kpoints(ikpoint),2,volume_cell_angstrom(),pub_dmft_temp,pub_dmft_points,mpi_size,nn,ien,energy,fermi_e,chem_shift,temp2_real
         endif
         close(30030+mpi_rank*pub_total_num_nodes+10200*is)
                    
      endif
      call comms_barrier
    enddo
   endif
  !======================================!
  !======================================!

  if(pub_dmft_optics_i1/=pub_dmft_optics_i2)then
#ifdef GPU_SPEEDUP
  if(use_gpu_onlyme)then
   if(pub_my_node_id==0.or.split)then
     CALL matmulcuda_r(temp1_real,nabla_square(:,:,pub_dmft_optics_i2),temp2_real,nn,nn,nn)
   endif
   if(pub_total_num_nodes>1.and..not.split) call comms_bcast(0,temp2_real)
  else
    temp2_real=MATMUL(temp1_real,nabla_square(:,:,pub_dmft_optics_i2))
  endif
#else
    temp2_real=MATMUL(temp1_real,nabla_square(:,:,pub_dmft_optics_i2))
#endif
  endif
  !======================================!

  if(pub_dmft_optics_i1/=pub_dmft_optics_i2)then
   if(.not.split)then
    if(pub_my_node_id==0)then
        write(30030+mpi_rank*pub_total_num_nodes+10200*is) nkpoints,norm_kpoints(ikpoint),2,volume_cell_angstrom(),pub_dmft_temp,pub_dmft_points,mpi_size,nn,ien,energy,fermi_e,chem_shift,temp2_real
    endif
   else
     do kkkk=1,pub_total_num_nodes
      if(kkkk==pub_my_node_id+1)then
         if(ien==ffirst.and.kkkk==1.and.  (   (.not.pub_dmft_splitk.and.(ikpoint==1.or.one_k_one_file))  .or.  (pub_dmft_splitk.and.ikpoint==splitk_first)  )   )then
           open(unit=30030+mpi_rank*pub_total_num_nodes+10200*is,file='optics_data_spin_'//TRIM(ADJUSTL(toString(is)))//'_k_'//TRIM(ADJUSTL(toString(ikpoint)))//'_rank_'//TRIM(ADJUSTL(toString(mpi_rank))),form='unformatted',status='replace') 
         else
           open(unit=30030+mpi_rank*pub_total_num_nodes+10200*is,file='optics_data_spin_'//TRIM(ADJUSTL(toString(is)))//'_k_'//TRIM(ADJUSTL(toString(ikpoint)))//'_rank_'//TRIM(ADJUSTL(toString(mpi_rank))),form='unformatted',position='append',status='old') 
         endif
          write(30030+mpi_rank*pub_total_num_nodes+10200*is) nkpoints,norm_kpoints(ikpoint),2,volume_cell_angstrom(),pub_dmft_temp,pub_dmft_points,mpi_size,nn,ien,energy,fermi_e,chem_shift,temp2_real 
         close(30030+mpi_rank*pub_total_num_nodes+10200*is)
      endif
      call comms_barrier
    enddo
   endif
  endif
  !======================================!

  if(.not.split)then
  if(pub_my_node_id==0)then
    close(30030+mpi_rank*pub_total_num_nodes+10200*is)
  endif
  endif

  endif
  !-----------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------!

  return
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine show_dos(greenf,overlap,energy)
  use sparse,   only: sparse_convert
#ifdef GPU_SPEEDUP
      use fortran_cuda, only : matmulcuda_c,matmulcuda_c_cublas
#endif
      implicit none
      integer                                       :: i,j,k,nn,kkkk,iplot
      type(SPAM3),intent(in)                        :: overlap,greenf
      complex(8)                                    :: energy
      complex(kind=DP), allocatable, dimension(:,:) :: mat_square,mat_square2,mat_square3
      real(kind=DP),allocatable                     :: totatom(:)
      real(kind=DP)                                 :: tot
      logical                                       :: computes_full_product
      real(8),allocatable                           :: totplanes(:)

      if(verboseall) write(*,*) 'show dos allocate'     
      allocate(mat_square2(ngwf_basis%num, ngwf_basis%num),stat=ierr)
      call utils_alloc_check('show_dos','mat_square2',ierr)
      allocate(mat_square3(ngwf_basis%num, ngwf_basis%num),stat=ierr)
      call utils_alloc_check('show_dos','mat_square3',ierr)
      if(.not.split)then
        call sparse_convert(mat_square2, greenf)
      else
        mat_square2=greenfbackup
      endif
      call sparse_convert(mat_square3,overlap)
      nn=ngwf_basis%num
      if(verboseall) write(*,*) 'done'
      allocate(totatom(nn))

!-------------------------!
!-------------------------!
!-------------------------!
!BUG Cedric Feb 2015
!if(pub_my_node_id==0)then

computes_full_product=.false.

!#############################!
if(computes_full_product)then
allocate(mat_square(ngwf_basis%num, ngwf_basis%num),stat=ierr)
call utils_alloc_check('show_dos','mat_square',ierr)

#ifdef GPU_SPEEDUP
   if(use_gpu_onlyme)then
    if(pub_my_node_id==0.or.split)then
    ! A = Trace(G * S)
      i=size(mat_square2,1); j=size(mat_square2,2); k=size(mat_square3,2)
      write(*,*) 'call matmul cuda'
!1)
      call matmulcuda_c(mat_square2,mat_square3,mat_square,i,j,k)  ! Using Magma, is fastest 
!2)
!     call matmulcuda_c_cublas(mat_square2,mat_square3,mat_square,i,j,k) ! Using Cublas,  
!3)
    ! call matmul_my_cuda_c(mat_square2,mat_square3,mat_square,i,j,k)
      if(verboseall) write(*,*) 'done'
     endif
     if(pub_total_num_nodes>1.and..not.split) call comms_bcast(0,mat_square)
   else
      mat_square=MATMUL(mat_square2,mat_square3)
   endif
#else
      if(verboseall) write(*,*) 'call CPU matmul'
      mat_square=MATMUL(mat_square2,mat_square3)
      if(verboseall) write(*,*) 'done'
#endif
      tot=0.0
      do i=1,nn
        totatom(i) = -aimag(mat_square(i,i))/pi
        tot        =  tot + totatom(i)
      enddo
      deallocate(mat_square,stat=ierr)
      call utils_dealloc_check('show_dos','mat_square',ierr)
else
     tot=0.0
     do i=1,nn
        totatom(i)  = - aimag(sum(mat_square2(i,:)*mat_square3(:,i))) / pi
        tot         =   tot + totatom(i)
     enddo
endif
!#############################!

     if(nkpoints>1)then
        totatom=totatom * norm_kpoints(ikpoint)
        tot    =tot     * norm_kpoints(ikpoint)
     endif


    if(.not.split)then
     if(pub_my_node_id==0)then
     !in eV
       if((ien==ffirst.and.one_k_one_file).or.(.not.one_k_one_file.and.ien==ffirst.and.ikpoint==1))then 
        open(unit=23,file='FULL_DOS_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='replace')
       else
        open(unit=23,file='FULL_DOS_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='old',position='append')
       endif
       write(23,'(i5,i3,i5,300f14.6)') ien,1,pub_dmft_points,(real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132,tot 
       close(23) 
     endif
    else
      do kkkk=1,pub_total_num_nodes
        if(kkkk==pub_my_node_id+1)then
          if(ien==ffirst.and.kkkk==1.and.  ( (.not.pub_dmft_splitk.and.(ikpoint==1.or.one_k_one_file))  .or.  (pub_dmft_splitk.and.ikpoint==splitk_first)  )    )then
           open(unit=23,file='FULL_DOS_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='replace')
          else
           open(unit=23,file='FULL_DOS_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),position='append',status='old')
          endif
          write(23,'(i5,i3,i5,300f14.6)') ien,1,pub_dmft_points,(real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132,tot
          close(23)
        endif
        call comms_barrier
      enddo
    endif

    if(nn<1000.and.plot_all_orbitals)then

    do iplot=1,nn
     if(.not.split)then
     if(pub_my_node_id==0)then
     !in eV
       if((ien==ffirst.and.one_k_one_file).or.(.not.one_k_one_file.and.ien==ffirst.and.ikpoint==1))then
        open(unit=23,file='DOS_ORB_'//TRIM(ADJUSTL(toString(iplot)))//'_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='replace')
        if(iplot==1) call plot_label_in_file(greenf)
       else
        open(unit=23,file='DOS_ORB_'//TRIM(ADJUSTL(toString(iplot)))//'_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='old',position='append')
       endif
       write(23,'(i5,i3,i5,300f14.6)') ien,1,pub_dmft_points,(real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132,totatom(iplot)
       close(23)
     endif
    else
      do kkkk=1,pub_total_num_nodes
        if(kkkk==pub_my_node_id+1)then
          if(ien==ffirst.and.kkkk==1.and.  ( (.not.pub_dmft_splitk.and.(ikpoint==1.or.one_k_one_file))  .or.  (pub_dmft_splitk.and.ikpoint==splitk_first)  )    )then
           open(unit=23,file='DOS_ORB_'//TRIM(ADJUSTL(toString(iplot)))//'_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='replace')
           if(iplot==1.and.pub_my_node_id==0) call plot_label_in_file(greenf)
          else
           open(unit=23,file='DOS_ORB_'//TRIM(ADJUSTL(toString(iplot)))//'_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),position='append',status='old')
          endif
          write(23,'(i5,i3,i5,300f14.6)') ien,1,pub_dmft_points,(real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132,totatom(iplot)
          close(23)
        endif
        call comms_barrier
      enddo
    endif
    enddo
    endif

    !----------------------------------------------------------------------!
     if(nplanes/=0)then
       allocate(totplanes(nplanes))
       totplanes=0.0
       do j=1,nplanes
        do i=1,nn
         if(orb_in_plane(j,i)) totplanes(j)=totplanes(j)+totatom(i)
        enddo
       enddo


     do i=1,nplanes 
      if(.not.split)then
       if(pub_my_node_id==0)then
         if((ien==ffirst.and.one_k_one_file).or.(.not.one_k_one_file.and.ien==ffirst.and.ikpoint==1))then
          open(unit=390000+i,file='FULL_DOS_PLANE_'//TRIM(ADJUSTL(toString(i)))&
           & //'_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='replace')
         else
          open(unit=390000+i,file='FULL_DOS_PLANE_'//TRIM(ADJUSTL(toString(i)))&
           & //'_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='old',position='append')
         endif
         write(390000+i,'(i5,i3,i5,300f14.6)') ien,1,pub_dmft_points, &
              & (real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132,totplanes(i)
         close(390000+i)
        endif
       else
         do kkkk=1,pub_total_num_nodes
          if(kkkk==pub_my_node_id+1)then
            if(ien==ffirst.and.kkkk==1.and.  ( (.not.pub_dmft_splitk.and.(ikpoint==1.or.one_k_one_file))  .or.  (pub_dmft_splitk.and.ikpoint==splitk_first)  )   )then
              open(unit=390000,file='FULL_DOS_PLANE_'//TRIM(ADJUSTL(toString(i)))&
                        & //'_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='replace')
            else
              open(unit=390000,file='FULL_DOS_PLANE_'//TRIM(ADJUSTL(toString(i)))&
                        & //'_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),position='append',status='old')
            endif
            write(390000,'(i5,i3,i5,300f14.6)') ien,1,pub_dmft_points, &
                        & (real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132,totplanes(i)
            close(390000)
          endif
          call comms_barrier
        enddo

       endif 
      enddo

       deallocate(totplanes)
     endif
    !----------------------------------------------------------------------!

!endif
!BUG Cedric Feb 2015

!-------------------------!
!-------------------------!
!-------------------------!

      if(verboseall) write(*,*) 'clean/deallocate'
      deallocate(mat_square2,stat=ierr)
      call utils_dealloc_check('show_dos','mat_square2',ierr)
      deallocate(mat_square3,stat=ierr)
      call utils_dealloc_check('show_dos','mat_square3',ierr)
      if(verboseall) write(*,*) 'done, return from show dos'
      deallocate(totatom)

    end subroutine show_dos


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     subroutine plot_label_in_file(greenf)
     implicit none
     TYPE(SPAM3),intent(in) :: greenf
     character(2),allocatable :: labels(:)
     integer :: iplot,nn,jj,iplot_,totorb,totorblast
      nn=NINT(sparse_num_rows(greenf))
      open(unit=99098,file='dos_labels')
      if(allocated(labels)) deallocate(labels)
      allocate(labels(nn))
      call sparse_show_matrix(greenf,outunit=99099,show_elems=.true.,labels=labels)
      do iplot=1,nn
       iplot_=0;totorb=0;totorblast=0
       do jj=1,pub_cell%nat
        totorb=totorb+ngwf_basis%num_on_atom(jj)
        if(iplot>totorblast.and.iplot<=totorb) then
          iplot_=jj
          exit 
        endif
        totorblast=totorb
       enddo
       iplot_=pub_orig_atom(iplot_)
       write(99098,'(i5,3f11.5,a4)') iplot,elements(iplot_)%centre%x*0.529177249, &
                                         & elements(iplot_)%centre%y*0.529177249, &
                                         & elements(iplot_)%centre%z*0.529177249,labels(iplot)
      enddo
      deallocate(labels)
      close(99098)
     end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_invert_full_gpu_real(mat_inv,mat,greenfunc)
#ifdef GPU_SPEEDUP
      use fortran_cuda, only : magma_fortran_double_singleprec_bypass 
#endif
      use sparse,only: sparse_convert,sparse_num_cols,sparse_num_rows

      implicit none
      integer :: nn
      type(SPAM3),intent(in)                     :: mat
      type(SPAM3),intent(inout)                  :: mat_inv
      real(8),allocatable, dimension(:,:)        :: mat_square
      logical                                    :: greenfunc

      if(.not.greenfunc)then
         write(*,*) 'internal_invert_full_gpu_real : only for green function ' ; stop
      endif

      allocate(mat_square(ngwf_basis%num, ngwf_basis%num))
      call sparse_convert(mat_square,mat)
      nn=NINT(sparse_num_rows(mat))

#ifdef GPU_SPEEDUP
  if(pub_my_node_id==0.or.split)then
     call magma_fortran_double_singleprec_bypass(nn,mat_square)     ! SinglePrec
  endif
  if(pub_total_num_nodes>1.and..not.split) call comms_bcast(0,mat_square)
#else
   write(*,*) 'not supposed to use GPU right? problem then...'
   stop
#endif
      call sparse_convert(mat_inv,mat_square)
      deallocate(mat_square)

    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine internal_invert_full_gpu(mat_inv,mat,overlap,energy,greenfunc)
#ifdef GPU_SPEEDUP
    use fortran_cuda, only : diago_cuda_it_c,matinv_magma_complex,magma_fortran_comp_,matinv_ge_complex
#endif
      use sparse,   only: sparse_convert
      implicit none
      integer                                       :: nn,kkkk,iplot
      type(SPAM3),intent(in)                        :: overlap
      complex(8)                                    :: energy
      real(kind=DP)                                 :: tot
      logical                                       :: computes_full_product
      type(SPAM3),intent(in)                        :: mat
      type(SPAM3),intent(inout)                     :: mat_inv
      complex(kind=DP), allocatable, dimension(:,:) :: mat_square,overlap_square
      real(kind=DP),allocatable                     :: totatom(:)
      real(8),allocatable                           :: totplanes(:)
      logical                                       :: greenfunc

#ifndef GPU_SPEEDUP
      write(*,*) 'not supposed to enter this branch of code : DMFT+onetep module'
      stop
#endif

      if(.not.greenfunc)then
         write(*,*) 'internal_invert_full_gpu : only for green function ' ; stop
      endif

      allocate(mat_square(ngwf_basis%num, ngwf_basis%num),stat=ierr)
      call utils_alloc_check('internal_invert_full_gpu','mat_square',ierr)
      call sparse_convert(mat_square,mat)
      nn=NINT(sparse_num_rows(mat))
      allocate(totatom(nn))

      if(verboseall) write(*,*) 'SIZE OF MATRIX TO DIAGONALIZE : ', nn
      if(verboseall) write(*,*) 'start diago NOW'
#ifdef GPU_SPEEDUP

! - Double prec for tail, where Green func is small -

if(pub_my_node_id==0.or.split)then
if( doubleprec(energy,ien) ) then
#ifdef GPU_HANGS
!     call diago_cuda_it_c(nn,mat_square)                          ! WORKS , 5.5-6sec
      call matinv_ge_complex(nn,mat_square)                        ! WORKS , 2.5-3.5 sec
#else
      call magma_fortran_comp_(nn,mat_square)                      ! WORKS , 2-3 sec (best), but more memory demanding
#endif
else
      call magma_fortran_comp_singleprec_bypass(nn,mat_square)     ! SinglePrec     
endif
endif
if(pub_total_num_nodes>1.and..not.split) call comms_bcast(0,mat_square)
#endif
      if(verboseall.and.ien==1) write(*,*) 'done'

!DOS
 if(pub_dmft_temp < 0.0_DP) then
     allocate(overlap_square(ngwf_basis%num, ngwf_basis%num),stat=ierr)
     call sparse_convert(overlap_square,overlap)
     !BUG CW FEB
     !if(pub_my_node_id==0)then
      tot=0.0
      do i=1,nn
        totatom(i) = -aimag(sum(mat_square(i,:)*overlap_square(:,i)))/pi
        tot        =  tot + totatom(i)
      enddo

      if(nkpoints>1)then
        totatom=totatom * norm_kpoints(ikpoint)
        tot    =tot     * norm_kpoints(ikpoint)
      endif

      if(.not.split)then
       if(pub_my_node_id==0)then
        !in eV
        if((ien==ffirst.and.one_k_one_file).or.(.not.one_k_one_file.and.ien==ffirst.and.ikpoint==1))then
         open(unit=22,file='FULL_DOS_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='replace')
        else
         open(unit=22,file='FULL_DOS_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='old',position='append')
        endif
        write(22,'(i5,i3,i5,300f14.6)') ien,1,pub_dmft_points,(real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132,tot
        close(22)
       endif
      else
       do kkkk=1,pub_total_num_nodes
        if(kkkk==pub_my_node_id+1)then
          if(ien==ffirst.and.kkkk==1.and.  ( (.not.pub_dmft_splitk.and.(ikpoint==1.or.one_k_one_file))  .or.  (pub_dmft_splitk.and.ikpoint==splitk_first)  )   )then
            open(unit=22,file='FULL_DOS_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='replace')
          else
            open(unit=22,file='FULL_DOS_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),position='append',status='old')
          endif
          write(22,'(i5,i3,i5,300f14.6)') ien,1,pub_dmft_points,(real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132,tot
          close(22)
        endif
        call comms_barrier
       enddo
      endif

! BUG CW FEB
!     endif

     deallocate(overlap_square,stat=ierr)

!plot all orbitals
    if(nn<1000.and.plot_all_orbitals)then
    do iplot=1,nn
      if(.not.split)then
       if(pub_my_node_id==0)then
        !in eV
        if((ien==ffirst.and.one_k_one_file).or.(.not.one_k_one_file.and.ien==ffirst.and.ikpoint==1))then
         open(unit=22,file='DOS_ORB_'//TRIM(ADJUSTL(toString(iplot)))//'_'//&
       & TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='replace')
         if(iplot==1) call plot_label_in_file(greenf)
        else
         open(unit=22,file='DOS_ORB_'//TRIM(ADJUSTL(toString(iplot)))//'_'//&
       & TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='old',position='append')
        endif
        write(22,'(i5,i3,i5,300f14.6)') ien,1,pub_dmft_points,(real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132,totatom(iplot)
        close(22)
       endif
      else
       do kkkk=1,pub_total_num_nodes
        if(kkkk==pub_my_node_id+1)then
          if(ien==ffirst.and.kkkk==1.and.  ( (.not.pub_dmft_splitk.and.(ikpoint==1.or.one_k_one_file))  .or.  (pub_dmft_splitk.and.ikpoint==splitk_first)  )   )then
            open(unit=22,file='DOS_ORB_'//TRIM(ADJUSTL(toString(iplot)))//'_'//&
          & TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='replace')
            if(iplot==1.and.pub_my_node_id==0) call plot_label_in_file(greenf)
          else
            open(unit=22,file='DOS_ORB_'//TRIM(ADJUSTL(toString(iplot)))//'_'//&
          & TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),position='append',status='old')
          endif
          write(22,'(i5,i3,i5,300f14.6)') ien,1,pub_dmft_points,(real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132,totatom(iplot)
          close(22)
        endif
        call comms_barrier
       enddo
      endif
    enddo
    endif


   !----------------------------------------------------------------------!
     if(nplanes/=0)then
       allocate(totplanes(nplanes))
       totplanes=0.0

       do j=1,nplanes
        do i=1,nn
         if(orb_in_plane(j,i)) totplanes(j)=totplanes(j)+totatom(i)
        enddo
       enddo

       if(.not.split)then
        do i=1,nplanes
         if(pub_my_node_id==0)then
          if((ien==ffirst.and.one_k_one_file).or.(.not.one_k_one_file.and.ien==ffirst.and.ikpoint==1))then
          open(unit=390000+i,file='FULL_DOS_PLANE_'//TRIM(ADJUSTL(toString(i)))&
             & //'_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='replace')
          else
          open(unit=390000+i,file='FULL_DOS_PLANE_'//TRIM(ADJUSTL(toString(i)))&
             & //'_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='old',position='append')
          endif 
          write(390000+i,'(i5,i3,i5,300f14.6)') ien,1,pub_dmft_points, &
               & (real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132,totplanes(i)
          close(390000+i)
         endif
        enddo
       else
        do i=1,nplanes
        do kkkk=1,pub_total_num_nodes
          if(kkkk==pub_my_node_id+1)then
            if(ien==ffirst.and.kkkk==1.and.   ( (.not.pub_dmft_splitk.and.(ikpoint==1.or.one_k_one_file))  .or.  (pub_dmft_splitk.and.ikpoint==splitk_first)  )       )then
             open(unit=390000,file='FULL_DOS_PLANE_'//TRIM(ADJUSTL(toString(i)))&
                        & //'_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),status='replace') 
            else
             open(unit=390000,file='FULL_DOS_PLANE_'//TRIM(ADJUSTL(toString(i)))&
                        & //'_'//TRIM(ADJUSTL(toString(is)))//'_rank'//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift))),position='append',status='old') 
            endif
                        write(390000,'(i5,i3,i5,300f14.6)') ien,1,pub_dmft_points, &
                        & (real(energy-(fermi_e(is)+pub_dmft_chem_shift)))*27.211396132,totplanes(i)
                        close(390000)
          endif
          call comms_barrier
        enddo
        enddo
       endif

       deallocate(totplanes)
     endif
    !----------------------------------------------------------------------!

 endif
!END DOS

      greenfbackupflag=.true.
      greenfbackup=mat_square

      call sparse_convert(mat_inv,mat_square)
      deallocate(mat_square,stat=ierr)
      call utils_dealloc_check('internal_invert_full_gpu', 'mat_square',ierr)
      deallocate(totatom)

    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine matmul_in_dense_format_real(matc,mata,matb)
#ifdef GPU_SPEEDUP
      use fortran_cuda, only : matmulcuda_c,matmulcuda_c_cublas,matmulcuda_r_singleprec_bypass
#endif
      use sparse, only: sparse_convert,sparse_num_cols,sparse_num_rows
      implicit none
      integer                                 :: i,nn,n1,n2,n3
      type(SPAM3),intent(in)                  :: mata,matb
      type(SPAM3),intent(inout)               :: matc
      real(8), allocatable, dimension(:,:)    :: matdensa,matdensb,matdensc

      n1 = NINT(sparse_num_rows(mata))  ! row A
      n2 = NINT(sparse_num_cols(mata))  ! sum
      n3 = NINT(sparse_num_cols(matb))  ! col B

#ifdef debugMIXv2
      if(NINT(sparse_num_rows(matb))/=n2)then
        write(*,*) 'dim do not match matmul'
        write(*,*) 'rows b, col a : ',NINT(sparse_num_rows(matb)),n2
        stop
      endif
      if(NINT(sparse_num_rows(matc))/=n1)then
        write(*,*) 'dim do not match matmul'
        write(*,*) 'rows c, rows a : ',NINT(sparse_num_rows(matc)),n1
        stop
      endif
      if(NINT(sparse_num_cols(matc))/=n3)then
        write(*,*) 'dim do not match matmul'
        write(*,*) 'cols c, cols b : ',NINT(sparse_num_cols(matc)),n3
        stop
      endif
#endif

#ifdef GPU_SPEEDUP
 if(use_gpu_onlyme)then
  allocate(matdensa(n1,n2),matdensb(n2,n3),matdensc(n1,n3),stat=ierr)
  call sparse_convert(matdensa,mata)
  call sparse_convert(matdensb,matb)
  if(pub_my_node_id==0.or.split)then
      if(use_gpu_onlyme)then
        call matmulcuda_r_singleprec_bypass(matdensa,matdensb,matdensc,n1,n2,n3)
      else
        matdensc=matmul(matdensa,matdensb)
      endif
  endif
  if(pub_total_num_nodes>1.and..not.split) call comms_bcast(0,matdensc)
  call sparse_convert(matc,matdensc)
  deallocate(matdensa,matdensb,matdensc,stat=ierr)
 else
      call sparse_product(matc,mata,matb)
 endif
#else
      call sparse_product(matc,mata,matb)
#endif

    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine matmul_in_dense_format(matc,mata,matb,aba,show_trace)
  ! C=A*B or C=A*B*A if flag ABA is present !
#ifdef GPU_SPEEDUP
      use fortran_cuda, only : matmulcuda_c,matmulcuda_c_cublas,matmulcuda_c_singleprec_bypass
#endif
      use sparse,   only: sparse_convert,sparse_num_cols,sparse_num_rows
      implicit none
      integer                                 :: i,nn,n1,n2,n3
      type(SPAM3),intent(in)                  :: mata,matb
      type(SPAM3),intent(inout)               :: matc
      complex(8), allocatable, dimension(:,:) :: matdensa,matdensb,matdensc
      complex(8)                              :: traceit
      logical,optional                        :: aba,show_trace
      character(2),allocatable                :: labels(:)

      n1 = NINT(sparse_num_rows(mata))  ! row A
      n2 = NINT(sparse_num_cols(mata))  ! sum 
      n3 = NINT(sparse_num_cols(matb))  ! col B

      if(ien==1.and..not.silent) write(*,*) 'allocate for GPUm n1n2n3 : ', n1,n2,n3
      allocate(matdensa(n1,n2),matdensb(n2,n3),matdensc(n1,n3),stat=ierr)
      call utils_alloc_check('matmul_gpu','matdens',ierr)

      ! Convert matrix to full square
      if(ien==1.and.verboseall) write(*,*) 'convert sparse matrices to dense complex'
      call sparse_convert(matdensa,mata)
      if(ien==1.and.verboseall) write(*,*) 'convert b'
      call sparse_convert(matdensb,matb)

      if(verboseall)write(*,*) 'go for gpu matmul'

      if(verboseall.and.ien==1.and..not.silent)then
        write(*,*) '------------------------------------'
        write(*,*) 'row A     : ', NINT(sparse_num_rows(mata))
        write(*,*) 'col A     : ', NINT(sparse_num_cols(mata))
        write(*,*) 'row B     : ', NINT(sparse_num_rows(matb))
        write(*,*) 'col B     : ', NINT(sparse_num_cols(matb))
        write(*,*) 'row C     : ', NINT(sparse_num_rows(matc))
        write(*,*) 'col C     : ', NINT(sparse_num_cols(matc))
        write(*,*) 'shape A   : ', shape(matdensa)
        write(*,*) 'shape B   : ', shape(matdensb)
        write(*,*) 'shape C   : ', shape(matdensc)
        write(*,*) 'n1n2n3    : ', n1,n2,n3
        write(*,*) 'allocated : ', allocated(matdensa),allocated(matdensb),allocated(matdensc)
        write(*,*) '------------------------------------'
      endif

#ifdef GPU_SPEEDUP
if(use_gpu_onlyme)then
 if(pub_my_node_id==0.or.split)then
    if( doubleprec(energy,ien) ) then
       if(verboseall) write(*,*) 'matmucuda_c'
       CALL matmulcuda_c(matdensa,matdensb,matdensc,n1,n2,n3)           ! fastest
     ! call matmulcuda_c_cublas(matdensa,matdensb,matdensc,n1,n2,n3)    ! can be a factor 2 slower
    else
       if(verboseall.and..not.silent) write(*,*) 'matmucuda_c_singleprec_bypass'
       call matmulcuda_c_singleprec_bypass(matdensa,matdensb,matdensc,n1,n2,n3)
    endif
    if(present(aba))then
       if(n1/=n2.or.n1/=n3)then
        write(*,*) 'ERROR onetep dmft module, trying produt A*B*A with rectangular matrices'
        stop
       endif
       CALL matmulcuda_c(matdensc,matdensa,matdensb,n1,n2,n3)  
       matdensc=matdensb    
    endif
  endif
  if(pub_total_num_nodes>1.and..not.split) call comms_bcast(0,matdensc)
else
    matdensc=matmul(matdensa,matdensb)
    if(present(aba))then
      matdensc=matmul(matdensc,matdensa)
    endif
endif
#else
    matdensc=matmul(matdensa,matdensb) 
    if(present(aba))then
     matdensc=matmul(matdensc,matdensa)
    endif
#endif

      if(present(show_trace))then
      if(show_trace)then
        if(allocated(labels)) deallocate(labels)
        allocate(labels(size(matdensc,1)))
        call sparse_show_matrix(matc,outunit=99099,show_elems=.true.,labels=labels)
        if(NINT(sparse_num_cols(matc))/=size(matdensc,1))then
          write(*,*) 'ERROR, something is not adding up...'
          write(*,*) 'number of labels : ', NINT(sparse_num_cols(matc))
          write(*,*) 'size of matrix   : ', size(matdensc,1)
          stop
        endif
        traceit=0.d0
        if(pub_my_node_id==0.and.mpi_rank==1)then
            write(*,*) 'TOTAL NUMBER OF ELEMENTS : ', size(matdensc,1)
            do i=1,size(matdensc,1)
              traceit=traceit+matdensc(i,i)
              write(*,'(a,i5,4f20.12)') labels(i),i,matdensc(i,i)
            enddo 
            write(*,'(a,4f20.12)') 'TRACE IS : ', traceit
        endif
      endif
      endif

      if(verboseall.and.ien==1.and..not.silent) write(*,*) 'convert back to sparse format'
      ! Convert matrix back to SPAM3
      call sparse_convert(matc,matdensc)
      if(verboseall.and.ien==1.and..not.silent) write(*,*) 'deallocate'
      deallocate(matdensa,matdensb,matdensc,stat=ierr)
      call utils_dealloc_check('matmul gpu','matdens',ierr)
      if(verboseall.and.ien==1.and..not.silent) write(*,*) 'done'

    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine matmul_in_dense_format_r(matc,mata,matb,aba,show_trace)
  ! C=A*B
#ifdef GPU_SPEEDUP
      use fortran_cuda, only : matmul_my_cuda,matmulcuda_r_cublas
#endif
      use sparse,   only: sparse_convert,sparse_num_cols,sparse_num_rows
      implicit none
      integer                              :: i,nn,n1,n2,n3
      type(SPAM3),intent(in)               :: mata,matb
      type(SPAM3),intent(inout)            :: matc
      real(8), allocatable, dimension(:,:) :: matdensa,matdensb,matdensc
      real(8)                              :: traceit
      logical,optional                     :: aba,show_trace
      character(2),allocatable             :: labels(:)

      n1 = NINT(sparse_num_rows(mata))   !row A
      n2 = NINT(sparse_num_cols(mata))   !sum
      n3 = NINT(sparse_num_cols(matb))   !col B

#ifdef debugMIXv2
      if(NINT(sparse_num_rows(matb))/=n2)then
        write(*,*) 'dim do not match matmul'
        write(*,*) 'rows b, col a : ',NINT(sparse_num_rows(matb)),n2
        stop
      endif
      if(NINT(sparse_num_rows(matc))/=n1)then
        write(*,*) 'dim do not match matmul'
        write(*,*) 'rows c, rows a : ',NINT(sparse_num_rows(matc)),n1
        stop
      endif
      if(NINT(sparse_num_cols(matc))/=n3)then
        write(*,*) 'dim do not match matmul'
        write(*,*) 'cols c, cols b : ',NINT(sparse_num_cols(matc)),n3
        stop
      endif
#endif

      if(verboseall)write(*,*) 'allocate for GPUm n1n2n3 : ', n1,n2,n3
      allocate(matdensa(n1,n2),matdensb(n2,n3),matdensc(n1,n3),stat=ierr)
      call utils_alloc_check('matmul_gpu','matdens',ierr)

      if(verboseall)write(*,*) 'convert sparse matrices to dense real'
      call sparse_convert(matdensa,mata)
      if(verboseall)write(*,*) 'convert b'
      call sparse_convert(matdensb,matb)

      if(verboseall)write(*,*) 'go for gpu matmul'

      if(verboseall.and.ien==1.and..not.silent)then
        write(*,*) '------------------------------------'
        write(*,*) 'row A     : ', NINT(sparse_num_rows(mata))
        write(*,*) 'col A     : ', NINT(sparse_num_cols(mata))
        write(*,*) 'row B     : ', NINT(sparse_num_rows(matb))
        write(*,*) 'col B     : ', NINT(sparse_num_cols(matb))
        write(*,*) 'row C     : ', NINT(sparse_num_rows(matc))
        write(*,*) 'col C     : ', NINT(sparse_num_cols(matc))
        write(*,*) 'shape A   : ', shape(matdensa)
        write(*,*) 'shape B   : ', shape(matdensb)
        write(*,*) 'shape C   : ', shape(matdensc)
        write(*,*) 'n1n2n3    : ', n1,n2,n3
        write(*,*) 'allocated : ', allocated(matdensa),allocated(matdensb),allocated(matdensc)
        write(*,*) '------------------------------------'
       endif

#ifdef debugMIXv3
      if(pub_my_node_id==0)then
        write(*,*)
        write(*,'(a,40f15.3)') 'some elements of matrix a : ', matdensa(1,1:min(20,size(matdensa,2)))
        write(*,*) 
        write(*,'(a,40f15.3)') 'some elements of matrix b : ', matdensa(1,1:min(20,size(matdensa,2)))
      endif
#endif

#ifdef GPU_SPEEDUP
if(use_gpu_onlyme)then
if(pub_my_node_id==0.or.split)then
!1)
               CALL matmulcuda_r(matdensa,matdensb,matdensc,n1,n2,n3)
!2)
  !     call matmulcuda_r_cublas(matdensa,matdensb,matdensc,n1,n2,n3)
        if(verboseall)write(*,*) 'done/matmulcuda_r'
!3)
!       call matmul_my_cuda_r(matdensa,matdensb,matdensc,n1,n2,n3)
        if(present(aba))then
            if(n1/=n2.or.n1/=n3)then
             write(*,*) 'ERROR onetep dmft module, trying produt A*B*A with rectangular matrices'
             stop
            endif
            CALL matmulcuda_r_cublas(matdensc,matdensa,matdensb,n1,n2,n3)           
            matdensc=matdensb
        endif
endif
if(pub_total_num_nodes>1.and..not.split) call comms_bcast(0,matdensc)
else
        matdensc=matmul(matdensa,matdensb)
        if(present(aba))then
            matdensc=matmul(matdensc,matdensa)
        endif
endif
#else
        matdensc=matmul(matdensa,matdensb)   
        if(present(aba))then
            matdensc=matmul(matdensc,matdensa)
        endif
#endif

     if(present(show_trace))then
      if(show_trace)then
        if(allocated(labels)) deallocate(labels)
        allocate(labels(size(matdensc,1))) 
        if(NINT(sparse_num_cols(matc))/=size(matdensc,1))then
         write(*,*) 'ERROR, something is not adding up...'
         write(*,*) 'number of labels : ', NINT(sparse_num_cols(matc))
         write(*,*) 'size of matrix   : ', size(matdensc,1)
         stop
        endif
        call sparse_show_matrix(matc,outunit=99099,show_elems=.true.,labels=labels) 
        traceit=0.d0
        if(pub_my_node_id==0.and.mpi_rank==1)then
            write(*,*) 'TOTAL NUMBER OF ELEMENTS : ', size(matdensc,1)
            do i=1,size(matdensc,1)
              traceit=traceit+matdensc(i,i)
              write(*,'(a,i5,4f20.12)') labels(i), i, matdensc(i,i)
            enddo
            write(*,'(a,4f20.12)') 'TRACE IS : ', traceit
        endif
      endif
      endif

      if(verboseall)write(*,*) 'convert back to sparse format'
      call sparse_convert(matc,matdensc)
      if(verboseall)write(*,*) 'deallocate'
      deallocate(matdensa,matdensb,matdensc,stat=ierr)
      call utils_dealloc_check('matmul gpu','matdens',ierr)
      if(verboseall)write(*,*) 'done'

    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine hubbard_greenf_info(hub_atom,hub_atomb,is,energy,ien,theatom,theatomb,species,speciesb,hat,hatb,n_,m_)
      use comms, only: comms_abort
      use hubbard_init, only:  h_species
      use simulation_cell, only: pub_cell
      implicit none
      integer,intent(in)           :: n_,m_
      integer, intent(in)          :: is,hub_atom,hub_atomb
      integer,intent(in)           :: ien,theatom,theatomb,hatb,hat
      integer                      :: channels,channelsb
      complex(kind=DP), intent(in) :: energy
      complex(kind=DP)             :: ttt(n_,m_)  
      complex(8),allocatable       :: ttt_buffer(:,:,:)
      integer,allocatable          :: fre_buffer(:),at_buffer(:)
      complex(8),allocatable       :: ener_buffer(:)
      real(8),allocatable          :: ener_bufferR(:),ener_bufferI(:)
      integer                      :: species,speciesb
      integer                      :: jj,kkkk,sizbuffer
      integer                      :: kk,stdout 

               if(verboseall) write(*,*) '--------------------------------------------------'
               if(verboseall) write(*,*) 'STARTING HUBBARD GREEN INFO'

               if(split)then
                 allocate(ttt_buffer(pub_total_num_nodes,n_,m_)); allocate(fre_buffer(pub_total_num_nodes),at_buffer(pub_total_num_nodes));
                 allocate(ener_buffer(pub_total_num_nodes)); 
                 allocate(ener_bufferR(pub_total_num_nodes)); allocate(ener_bufferI(pub_total_num_nodes));
                 ener_buffer=0;ener_bufferR=0.0;ener_bufferI=0.0;fre_buffer=0;ttt_buffer=0;at_buffer=0
               endif

               if(split)then
                 channels  = totn_atom_merge_split(hub_atom)
                 channelsb = totn_atom_merge_split(hub_atomb)
               else
                 channels  = totn_atom_merge(hub_atom)
                 channelsb = totn_atom_merge(hub_atomb)
               endif

               if(verboseall) write(*,*) 'merging blocks,channels : ',channels,channelsb
               if(verboseall) write(*,*) 'n_,m_  : ',n_,m_

               if(channels/=size(ttt,1).or.channelsb/=size(ttt,2))then
                 write(*,*) 'ERROR shape of ttt is not matching channels'
                 write(*,*) 'channels,channelsb : ',channels,channelsb
                 write(*,*) 'shape ttt          : ', shape(ttt)  
                 stop
               endif

               if(split)then
                ttt       = merged_green_block_split(hub_atom,hub_atomb,channels,channelsb)
               else
                ttt       = merged_green_block(hub_atom,hub_atomb,channels,channelsb) 
               endif
               if(verboseall) write(*,*) 'done,frequ : ',ien
#ifdef debug
               if(pub_my_node_id==0.or..not.split) call check_not_opened(1001)
#endif
               stdout=20000+pub_my_node_id+mpi_rank*5000
               filename='green_output'//TRIM(ADJUSTL(   toString(hat)))// &
                                 & "_"//TRIM(ADJUSTL(toString(species)))// &
                                 & "_"//TRIM(ADJUSTL(     toString(is)))
               if(hub_atom/=hub_atomb)then
                  if(verboseall)write(*,*) 'dimer with atoms : ', hub_atom,hub_atomb
                  filename=TRIM(ADJUSTL(filename))//"_dimer"
               endif

               if(mpi_size>1.or.nkpoints>1)then
                  filename=TRIM(ADJUSTL(filename))//"_rank"//TRIM(ADJUSTL(toString(mpi_rank+kp_file_shift)))
               endif

               if(verboseall) write(*,*) 'opening filename : ', trim(adjustl(filename))

              if(ien==ffirst.and.  ( (.not.pub_dmft_splitk.and.(ikpoint==1.or.one_k_one_file))  .or.  (pub_dmft_splitk.and.ikpoint==splitk_first)  )  )then
               if(pub_my_node_id==0.or..not.split) then
                 if(verboseall)then
                  write(*,*) 'RANK 0 OPENING GREEN_OUTPUT FILE : ', trim(adjustl(filename))
                  write(*,*) 'CHEM POTENTIAL                   : ', fermi_e(is) 
                  write(*,*) 'FREQUENCY                        : ', energy
                  write(*,*) 'FREQUENCY LABEL  / FIRST         : ', ien, ffirst
                  write(*,*) 'KPOINT                           : ', ikpoint 
                 endif
                 open(unit=1001,file=filename,form='unformatted',STATUS='REPLACE') 
               endif
              else
#ifdef debug
               if(pub_my_node_id==0.or..not.split) call check_files_exists(filename)
#endif
               if(pub_my_node_id==0.or..not.split) then
                open(unit=1001,file=filename,form='unformatted',position='append',status='old')
                if(verboseall)then
                 write(*,*) 'APPENDING INTO FILE : ', trim(adjustl(filename))
                 write(*,*) 'FREQUENCY NUMBER    : ', ien
                 write(*,*) 'KPOINT              : ', ikpoint
                endif
               endif
              endif

              if(nkpoints>1)then
                ttt = ttt * norm_kpoints(ikpoint)
              endif        

              if(split)then
               ttt_buffer(pub_my_node_id+1,:,:)=ttt
               fre_buffer(pub_my_node_id+1)=ien
               at_buffer(pub_my_node_id+1)=hat
               ener_bufferR(pub_my_node_id+1)=real(energy)
               ener_bufferI(pub_my_node_id+1)=aimag(energy)
               call comms_reduce('SUM',ttt_buffer(:,:,:))
               call comms_reduce('SUM',fre_buffer)
               call comms_reduce('SUM',at_buffer)
               call comms_reduce('SUM',ener_bufferR)
               call comms_reduce('SUM',ener_bufferI)
               ener_buffer=ener_bufferR+CMPLX(0.d0,1.d0,8)*ener_bufferI
               if(ANY( abs(at_buffer - at_buffer(1)) > 0 ))then
                 write(*,*) 'ERROR split and WE ARE NOT doing the same atom!!!!'
                 write(*,*) 'atoms on node : ', at_buffer
                 stop
               endif
              endif

              if(split)then
                sizbuffer=pub_total_num_nodes
              else
                sizbuffer=1
              endif

              !---------------------!
              do kkkk=1,sizbuffer

              if(split)then
                if(pub_my_node_id==0) write(1001) fre_buffer(kkkk),pub_dmft_points,pub_dmft_temp,channels,channelsb,hat,hatb
              else
                                      write(1001) ien             ,pub_dmft_points,pub_dmft_temp,channels,channelsb,hat,hatb
              endif

               if(num_spins==1)then !paramagnetic
                if(split)then
                 if(pub_my_node_id==0)then
                      write(1001) merged_occupancy_block_split(hub_atom,1,totn_atom_merge_split(hub_atom)), &
                                & merged_occupancy_block_split(hub_atom,1,totn_atom_merge_split(hub_atom))
                 endif
                else
                      write(1001) merged_occupancy_block(hub_atom,1,totn_atom_merge(hub_atom)), &
                                & merged_occupancy_block(hub_atom,1,totn_atom_merge(hub_atom))
                endif
               elseif(num_spins==2)then !spin up and dn
                 if(split)then
                  if(pub_my_node_id==0)then
                     write(1001) merged_occupancy_block_split(hub_atom,1,totn_atom_merge_split(hub_atom)), &
                               & merged_occupancy_block_split(hub_atom,2,totn_atom_merge_split(hub_atom))
                  endif
                 else
                     write(1001) merged_occupancy_block(hub_atom,1,totn_atom_merge(hub_atom)), &
                               & merged_occupancy_block(hub_atom,2,totn_atom_merge(hub_atom))
                 endif
               else
                 write(*,*) 'spin should be 1 or 2,critical'
                 stop
               endif

               if(ien==1) then
                  write(stdout,'(/a)')'########################################'
                  write(stdout,'(a,i6,a,a)') 'DFT+U Greens function of atom ',hat,' of Hubbard species ', h_species(species)%hub_species
                  write(stdout,*) ' hub_atom,hub_atomb            : ',hub_atom,hub_atomb
                  write(stdout,*) ' hat and theatom               : ',hat,theatom
                  write(stdout,*) ' hat                           : ',hat
                  write(stdout,*) ' species                       : ',species
                  write(stdout,'(a)')'########################################'
               endif

               if(split)then
                if(pub_my_node_id==0)then
#ifdef debug
                  write(*,*) 'WRITING TO FILE : ', (fermi_e(is)+pub_dmft_chem_shift),ener_buffer(kkkk)
                  write(*,*) 'WRITING TO FILE : ',  ttt_buffer(kkkk,1:channels,1:channelsb)
#endif
                  write(1001) (fermi_e(is)+pub_dmft_chem_shift),ener_buffer(kkkk),ttt_buffer(kkkk,1:channels,1:channelsb)
                endif
                ttt=ttt_buffer(kkkk,1:channels,1:channelsb)
               else
                  write(1001) (fermi_e(is)+pub_dmft_chem_shift),energy           ,ttt(1:channels,1:channelsb) 
               endif

              if(.not.split)then
               if(ien==1) write(stdout,'(a,i6,a,f15.8,a,f15.8)') 'Greens function of spin ',is,' and energy ',REAL(energy),' + i ',AIMAG(energy)
              else
               if(ien==1) write(stdout,'(a,i6,a,f15.8,a,f15.8)') 'Greens function of spin ',is,' and energy ',REAL(ener_buffer(kkkk)),' + i ',AIMAG(ener_buffer(kkkk))
              endif
               if(ien==1) write(stdout,*) '----------------------------------------'
               if(ien==1) write(stdout,'(a)') 'REAL PART'
               if(ien==1) write(stdout,*) 'channels,channelsb = ',channels,channelsb
               do jj=1,channels
                if(ien==1) write(stdout,'(300f15.8)')(REAL(ttt(jj,kk)),kk=1,channelsb)
               enddo
               if(ien==1) write(stdout,'(a)') 'IMAGINARY PART'
               do jj=1,channels
                if(ien==1) write(stdout,'(300f15.8)')(AIMAG(ttt(jj,kk)),kk=1,channelsb)
               enddo

               if(ien==1) write(stdout,*) '----------occupation matrix-------------'
               if(split)then
                ttt=merged_occupancy_block_split(hub_atom,1,totn_atom_merge_split(hub_atom))
               else
                ttt=merged_occupancy_block(hub_atom,1,totn_atom_merge(hub_atom))
               endif
               do jj=1,channels
                 if(ien==1) write(stdout,'(300f15.8)')(real(ttt(jj,kk)),kk=1,channelsb)
               enddo
               if(ien==1) write(stdout,'(a)')'########################################'

              enddo
              !---------------------!

               if(pub_my_node_id==0.or..not.split) close(1001)   
               if(verboseall)write(*,*) '--------------------------------------------------'
               if(split) deallocate(ttt_buffer,ener_buffer,fre_buffer,ener_bufferR,ener_bufferI,at_buffer)

    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine check_sparse_in_file(mat,filename)
   implicit none
   type(spam3)    :: mat
   character*(*)  :: filename
   logical        :: check
   if(ien==1.and..not.silent) write(*,*) 'checking if file [x] is present : ', filename
   call comms_barrier
   INQUIRE(file=filename,EXIST=check)
   call comms_barrier
   if(.not.check)then
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'file is not present, write current data'
    call sparse_write_(mat,filename=filename)
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'done...'
   else
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'file is present, read from this file'
    call sparse_read_(mat,filename=filename)
    if(.not.silent.and.pub_my_node_id==0) write(*,*) 'done...'
   endif
   call comms_barrier
   if(ien==1.and..not.silent) write(*,*) 'check is done'
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine sparse_write_(mat,filename)
   use sparse,   only: sparse_convert,sparse_num_rows,sparse_num_cols
   implicit none
   integer        :: n1,n2
   type(spam3)    :: mat
   real(kind=DP), allocatable, dimension(:,:) :: mat_square
   character*(*)  :: filename
    call comms_barrier
    n1=NINT(sparse_num_rows(mat))
    n2=NINT(sparse_num_cols(mat))
    if(n1==0.and.n2==0)then
      write(*,*) 'SPARSE WRITE, return, 0-shape matrix'
      return
    endif
    if(pub_my_node_id==0.and.mpi_rank==1) open(unit=20001,file=filename,form='unformatted')
      if(verboseall.and.pub_my_node_id==0) write(*,*) 'writing sparse matrix with geometry : ', n1,n2
      allocate(mat_square(n1,n2),stat=ierr)
      call sparse_convert(mat_square,mat)
#ifdef debugMIXv2
     if(maxval(abs(mat_square))<1.d-8)then
        write(*,*) 'ERROR/WARNING : WRITING 0 MATRIX INTO FILE : ', trim(adjustl(filename))
     endif
     write(*,*) 'MIN-MAX VALUES OF MATRIX WRITTEN : ', minval(abs(mat_square)),maxval(abs(mat_square))
#endif
       if(pub_my_node_id==0.and.mpi_rank==1)write(20001) mat_square
      deallocate(mat_square,stat=ierr)
      if(verboseall.and.pub_my_node_id==0) write(*,*) 'done, return'
     if(pub_my_node_id==0.and.mpi_rank==1)close(20001)
   end subroutine
             !-------------!
   subroutine sparse_read_(mat,filename)
   use sparse,   only: sparse_convert,sparse_num_rows,sparse_num_cols
   implicit none
   type(spam3)    :: mat
   real(kind=DP), allocatable, dimension(:,:) :: mat_square
   integer        :: n1,n2
   character*(*)  :: filename
      call comms_barrier
      n1=NINT(sparse_num_rows(mat))
      n2=NINT(sparse_num_cols(mat))
      if(n1==0.and.n2==0)then
        write(*,*) 'SPARSE READ, return, 0-shape matrix'
        return
      endif
      open(unit=20001,file=filename,form='unformatted')
      if(verboseall.and.pub_my_node_id==0) write(*,*) 'reading sparse matrix with geometry : ', n1,n2
      allocate(mat_square(n1,n2),stat=ierr)
      read(20001) mat_square
#ifdef debugMIXv2
     if(maxval(abs(mat_square))<1.d-8)then
        write(*,*) 'WARNING READING 0 MATRIX INTO FILE : ', trim(adjustl(filename))
     endif
     write(*,*) 'MIN-MAX VALUES OF MATRIX READ : ', minval(abs(mat_square)),maxval(abs(mat_square))
#endif
      call sparse_convert(mat,mat_square)
      deallocate(mat_square,stat=ierr)
      if(verboseall.and.pub_my_node_id==0) write(*,*) 'done, return'
    close(20001)
   end subroutine
             !-------------!
   subroutine sparse_readc_(mat,filename)
   use sparse,   only: sparse_convert,sparse_num_rows,sparse_num_cols
   implicit none
   type(spam3)    :: mat
   complex(kind=DP), allocatable, dimension(:,:) :: mat_square
   integer        :: n1,n2
   character*(*)  :: filename
      call comms_barrier
      n1=NINT(sparse_num_rows(mat))
      n2=NINT(sparse_num_cols(mat))
      if(n1==0.and.n2==0)then
        write(*,*) 'SPARSE READ, return, 0-shape matrix'
        return
      endif
      open(unit=20001,file=filename,form='unformatted')
      if(verboseall.and.pub_my_node_id==0) write(*,*) 'reading sparse matrix with geometry : ', n1,n2
      allocate(mat_square(n1,n2),stat=ierr)
      read(20001) mat_square
#ifdef debugMIXv2
     if(maxval(abs(mat_square))<1.d-8)then
        write(*,*) 'WARNING READING 0 MATRIX INTO FILE : ', trim(adjustl(filename))
     endif
     write(*,*) 'MIN-MAX VALUES OF MATRIX READ : ', minval(abs(mat_square)),maxval(abs(mat_square))
#endif
      call sparse_convert(mat,mat_square)
      deallocate(mat_square,stat=ierr)
      if(verboseall.and.pub_my_node_id==0) write(*,*) 'done, return'
    close(20001)
   end subroutine
            !-------------!
   subroutine sparse_realpartc_(mat)
   use sparse,   only: sparse_convert,sparse_num_rows,sparse_num_cols
   implicit none
   type(spam3)                                   :: mat
   complex(kind=DP), allocatable, dimension(:,:) :: mat_square
   integer                                       :: n1,n2
      call comms_barrier
      n1=NINT(sparse_num_rows(mat))
      n2=NINT(sparse_num_cols(mat))
      allocate(mat_square(n1,n2),stat=ierr)
      call sparse_convert(mat_square,mat)
      mat_square=real(mat_square)
      call sparse_convert(mat,mat_square)
      deallocate(mat_square)
   end subroutine
            !-------------!
   subroutine check_sparsec_(mat)
   use sparse,   only: sparse_convert,sparse_num_rows,sparse_num_cols
   implicit none
   type(spam3)    :: mat
   complex(kind=DP), allocatable, dimension(:,:) :: mat_square
   integer        :: n1,n2
      call comms_barrier
      n1=NINT(sparse_num_rows(mat))
      n2=NINT(sparse_num_cols(mat))
      if(n1==0.and.n2==0)then
        write(*,*) 'WARNING check sparse matrix, 0-shape matrix'
        return
      endif
      allocate(mat_square(n1,n2),stat=ierr)
      call sparse_convert(mat_square,mat)
     if(maxval(abs(mat_square))<1.d-8)then
        write(*,*) 'WARNING check sparse MATRIX - empty matrix'
     endif
     write(*,*) 'WARNING : MIN-MAX VALUES OF MATRIX IN CHECK SPARSE abs  : ', minval(abs(mat_square)),maxval(abs(mat_square))
     write(*,*) 'WARNING : MIN-MAX VALUES OF MATRIX IN CHECK SPARSE imag : ', minval(abs(aimag(mat_square))),maxval(abs(aimag(mat_square)))
     write(*,*) 'WARNING : MIN-MAX VALUES OF MATRIX IN CHECK SPARSE real : ', minval(abs(real(mat_square))),maxval(abs(real(mat_square)))
     write(*,*) 'WARNING : SUM VALUES                                    : ', sum(mat_square)
     write(*,*) 'WARNING : SUM VALUES (xline number)                     : ', sum( (/( sum(mat_square(i,:)*dble(i)),i=1,size(mat_square,1) )/)   )
     write(*,*) 'WARNING : SUM VALUES (xcol number)                      : ', sum( (/( sum(mat_square(:,i)*dble(i)),i=1,size(mat_square,1) )/)   )

     deallocate(mat_square,stat=ierr)
   end subroutine

           !-------------!
   subroutine check_sparser_(mat)
   use sparse,   only: sparse_convert,sparse_num_rows,sparse_num_cols
   implicit none
   type(spam3)    :: mat
   real(kind=DP), allocatable, dimension(:,:) :: mat_square
   integer        :: n1,n2
      call comms_barrier
      n1=NINT(sparse_num_rows(mat))
      n2=NINT(sparse_num_cols(mat))
      if(n1==0.and.n2==0)then
        write(*,*) 'WARNING check sparse matrix, 0-shape matrix'
        return
      endif
      allocate(mat_square(n1,n2),stat=ierr)
      call sparse_convert(mat_square,mat)
     if(maxval(abs(mat_square))<1.d-8)then
        write(*,*) 'WARNING check sparse MATRIX - empty matrix'
     endif
     write(*,*) 'MIN-MAX VALUES OF MATRIX IN CHECK SPARSE : ', minval(abs(mat_square)),maxval(abs(mat_square))
     write(*,*) 'SUM VALUES                               : ', sum(mat_square)
     deallocate(mat_square,stat=ierr)
   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine normalize_projectors(mat,proj_on_line_or_col,mat_t)
   use sparse,   only : sparse_convert,sparse_num_rows,sparse_num_cols
   use wrappers, only : wrappers_invert_sym_matrix
#ifdef GPU_SPEEDUP
  use fortran_cuda, only : eigenvector_cuda_r
#endif
   implicit none
   type(spam3)                                :: mat
   type(spam3),optional                       :: mat_t
   real(kind=DP), allocatable, dimension(:,:) :: mat_square,Sinv
   real(kind=DP), allocatable, dimension(:,:) :: L_LT
   integer                                    :: nn1,nn2,ierr,n1,n2,proj_on_line_or_col,i,j,k
   real(kind=DP)                              :: normhub
   logical                                    :: use_cholesky,u_dagger_u,flip,check
   real(8),allocatable                        :: vaps(:),WORK(:),temp(:,:)

      use_cholesky=abs(pub_dmft_norm_proj)==2.or.abs(pub_dmft_norm_proj)==3.or.abs(pub_dmft_norm_proj)==4
      u_dagger_u=abs(pub_dmft_norm_proj)==3

      call comms_barrier
      n1=NINT(sparse_num_rows(mat))
      n2=NINT(sparse_num_cols(mat))
      if(n1==0.and.n2==0)then
         write(*,*) 'Onetep, DMFT module : PROJECTORS ARE VOID'
         stop
      endif
      allocate(mat_square(n1,n2),stat=ierr)
      call sparse_convert(mat_square,mat)

      write(*,*) 'USE CHOLESKY? ', use_cholesky

      if(use_cholesky)then
         write(*,*) 'USING ORTHONORMALIZATION OF PROJECTORS'
         if((n1<n2.and.u_dagger_u).or.(n1>n2.and..not.u_dagger_u))then
          nn1=n1
          nn2=n2
          flip=.false.
         else
          nn1=n2
          nn2=n1
          flip=.true.
         endif

         if(allocated(L_LT)) deallocate(L_LT); allocate(L_LT(nn1,nn1))
         if(use_gpu_onlyme)then
#ifdef GPU_SPEEDUP
          if(.not.flip)then
           i=size(mat_square,1); j=size(mat_square,2); k=size(L_LT,2)
           if(pub_my_node_id==0)then
            CALL matmulcuda_r(mat_square,transpose(mat_square),L_LT,i,j,k)
           endif
          else
           i=size(mat_square,2); j=size(mat_square,1); k=size(L_LT,2)
           if(pub_my_node_id==0.or.split)then
            CALL matmulcuda_r(transpose(mat_square),mat_square,L_LT,i,j,k)
           endif
          endif
          if(pub_total_num_nodes>1.and..not.split) call comms_bcast(0,L_LT)
#else
          write(*,*) 'error, use GPU but no compiled with libraries'
#endif
         else
          if(.not.flip)then
           L_LT=matmul(mat_square,transpose(mat_square))
          else
           L_LT=matmul(transpose(mat_square),mat_square)
          endif
         endif
         write(*,*) 'O symmetric?', maxval(abs((L_LT-transpose(L_LT))/2.d0))
         L_LT=(L_LT+transpose(L_LT))/2.d0
         if(allocated(WORK)) deallocate(WORK,vaps,temp)
         allocate(WORK(3*nn1),vaps(nn1),temp(nn1,nn1))
         temp=L_LT
         do i=1,size(temp,1);temp(i,i)=0.d0;enddo
         if(maxval(abs(temp))>1.d-3)then
            write(*,*) 'L_LT is not diagonal'
            write(*,*) 'max off-diagonal : ', maxval(abs(temp))
            write(*,*) 'max     diagonal : ', maxval(abs(diag(L_LT)))
            if(use_gpu_eigenvectors.and.use_gpu_onlyme)then
             if(pub_my_node_id==0.or.split)then
#ifdef GPU_SPEEDUP
              call eigenvector_cuda_r(nn1,L_LT,vaps,temp,.false.)
#else
              write(*,*) 'error gpu 3' ; stop
#endif
#ifdef debug
              call DSYEV('V','U',nn1,L_LT,nn1,vaps,WORK,3*nn1,ierr)
              if(maxval(abs(temp-L_LT))>1.d-5)then
                write(*,*) 'ERROR IN EIGENVECTORS LLT!'
                stop
              endif
#endif
              L_LT=temp
             endif
             if(pub_total_num_nodes>1.and..not.split) call comms_bcast(0,L_LT)
             if(pub_total_num_nodes>1.and..not.split) call comms_bcast(0,vaps)
            else
             call DSYEV('V','U',nn1,L_LT,nn1,vaps,WORK,3*nn1,ierr)
             write(*,*) 'LAPACK DIAGONALIZATION INFO : ', ierr
            endif
         else
            write(*,*) 'L_LT is diagonal'
            vaps=diag(L_LT)
         endif
         vaps=1.d0/sqrt(abs(vaps))
         temp=0.d0;do i=1,nn1;temp(i,i)=vaps(i);enddo;
         if(use_gpu_onlyme)then
#ifdef GPU_SPEEDUP
          if(pub_my_node_id==0.or.split)then
           i=size(temp,1); j=size(temp,2); k=size(temp,2)
           CALL matmulcuda_r(temp,transpose(L_LT),temp,i,j,k)
           i=size(L_LT,1); j=size(L_LT,2); k=size(L_LT,2)
           CALL matmulcuda_r(L_LT,temp,L_LT,i,j,k)
           if(.not.flip)then 
            i=size(L_LT,1); j=size(L_LT,2); k=size(mat_square,2)
            CALL matmulcuda_r(L_LT,mat_square,mat_square,i,j,k)
           else
            i=size(mat_square,1); j=size(mat_square,2); k=size(mat_square,2)
            CALL matmulcuda_r(mat_square,L_LT,mat_square,i,j,k)
           endif
          endif
          if(pub_total_num_nodes>1.and..not.split) call comms_bcast(0,mat_square) 
#else
          write(*,*) 'error, use GPU but no compiled with libraries'
#endif
         else
           L_LT=matmul( matmul(L_LT,  temp) , transpose(L_LT) )
           if(.not.flip)then
            mat_square=matmul(L_LT,mat_square)
           else
            mat_square=matmul(mat_square,L_LT)
           endif
         endif

        if(abs(pub_dmft_norm_proj)==4)then
          if(allocated(Sinv)) deallocate(Sinv)
           allocate(Sinv(ngwf_basis%num, ngwf_basis%num))
           call sparse_convert(Sinv,overlap)
           if(allocated(WORK)) deallocate(WORK,vaps,temp)
           allocate(WORK(3*ngwf_basis%num),vaps(ngwf_basis%num),temp(ngwf_basis%num,ngwf_basis%num))
           inquire(file='store_sinv',exist=check)
           if(.not.check)then
             call DSYEV('V','U',ngwf_basis%num,Sinv,ngwf_basis%num,vaps,WORK,3*ngwf_basis%num,ierr)
             if(mpi_rank==1.and.pub_my_node_id==0) then
               open(unit=33981,file='store_sinv',form='unformatted')
               write(33981) Sinv,vaps
               close(33981)
             endif
           else
             open(unit=33981,file='store_sinv',form='unformatted')
             read(33981)  Sinv,vaps
             close(33981)
           endif
           vaps=sqrt(abs(vaps))
           temp=0.d0;do i=1,ngwf_basis%num;temp(i,i)=vaps(i);enddo;
           if(use_gpu_onlyme)then
#ifdef GPU_SPEEDUP
            if(pub_my_node_id==0.or.split)then
             i=size(temp,1); j=size(temp,2); k=size(temp,2)
             CALL matmulcuda_r(temp,transpose(Sinv),temp,i,j,k)
             i=size(Sinv,1); j=size(Sinv,2); k=size(temp,2)
             CALL matmulcuda_r(Sinv,temp,temp,i,j,k)
             if(size(mat_square,1)==ngwf_basis%num)then
                i=size(temp,1); j=size(temp,2); k=size(mat_square,2)
                CALL matmulcuda_r(temp,mat_square,mat_square,i,j,k)
             elseif(size(mat_square,2)==ngwf_basis%num)then
                i=size(mat_square,1); j=size(mat_square,2); k=size(mat_square,2)
                CALL matmulcuda_r(mat_square,temp,mat_square,i,j,k)
             elseif(.true.)then
                write(*,*) 'oups this is not supposed to happen...'
                write(*,*) 'shape mat_square : ', shape(mat_square)
                write(*,*) 'shape Sinv       : ', shape(Sinv)
                stop
             endif
            endif
            if(pub_total_num_nodes>1.and..not.split) call comms_bcast(0,mat_square)
#else
          write(*,*) 'error, use GPU but no compiled with libraries'
#endif
           else
             temp=matmul(matmul(Sinv,temp),transpose(Sinv))
             if(size(mat_square,1)==ngwf_basis%num)then
                mat_square=matmul(temp,mat_square)
             elseif(size(mat_square,2)==ngwf_basis%num)then
                mat_square=matmul(mat_square,temp)
             elseif(.true.)then
                write(*,*) 'oups this is not supposed to happen...'
                write(*,*) 'shape mat_square : ', shape(mat_square)
                write(*,*) 'shape Sinv       : ', shape(Sinv)
                stop
             endif
           endif
         endif

     endif

     if(.not.use_cholesky)then
      if(n1>n2)then
        if(.not.proj_on_line_or_col==2)then
          write(*,*) 'Onetep, DMFT module : more projectors than the NGWFs?'
          stop
        endif
        write(*,*) 'SCANNING HUBOVERLAP MATRIX, COLUMNS'
        write(*,*) 'WARNING: please check the normalization of the projectors...'
        do i=1,n2
          normhub=sum(abs(mat_square(:,i))**2)
          mat_square(:,i) = mat_square(:,i) / sqrt(normhub)
        enddo
      else
        if(.not.proj_on_line_or_col==1)then
          write(*,*) 'Onetep, DMFT module : more projectors than the NGWFs?'
          stop
        endif
        write(*,*) 'SCANNING HUBOVERLAP MATRIX, LINES'
        do i=1,n1
          normhub=sum(abs(mat_square(i,:))**2)
          mat_square(i,:) = mat_square(i,:) / sqrt(normhub)
        enddo
      endif
     endif

     if(present(mat_t))then
      call sparse_convert(mat_t,transpose(mat_square))
     endif

     call sparse_convert(mat,mat_square)
     deallocate(mat_square,stat=ierr)

     write(*,*) 'projectors renormalized, return'

   end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(8) function volume_cell_angstrom()
  implicit none
  real(8)  :: real_lattice(3,3),v(3)
  real_lattice(1,1) = pub_cell%a1%x*0.529177249
  real_lattice(1,2) = pub_cell%a1%y*0.529177249
  real_lattice(1,3) = pub_cell%a1%z*0.529177249
  real_lattice(2,1) = pub_cell%a2%x*0.529177249
  real_lattice(2,2) = pub_cell%a2%y*0.529177249
  real_lattice(2,3) = pub_cell%a2%z*0.529177249
  real_lattice(3,1) = pub_cell%a3%x*0.529177249
  real_lattice(3,2) = pub_cell%a3%y*0.529177249
  real_lattice(3,3) = pub_cell%a3%z*0.529177249
  v = vecprod( real_lattice(1,:), real_lattice(2,:) )
  volume_cell_angstrom = scalprod( v, real_lattice(3,:) )
  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(8) function distance(hat,hatb,uu)
  implicit none
  integer  :: hat,hatb,i,j,k
  real(8)  :: real_lattice(3,3),dist(-1:1,-1:1,-1:1),v(3)
  integer  :: uu(3)
  real_lattice(1,1) = pub_cell%a1%x*0.529177249  
  real_lattice(1,2) = pub_cell%a1%y*0.529177249
  real_lattice(1,3) = pub_cell%a1%z*0.529177249
  real_lattice(2,1) = pub_cell%a2%x*0.529177249
  real_lattice(2,2) = pub_cell%a2%y*0.529177249
  real_lattice(2,3) = pub_cell%a2%z*0.529177249
  real_lattice(3,1) = pub_cell%a3%x*0.529177249
  real_lattice(3,2) = pub_cell%a3%y*0.529177249
  real_lattice(3,3) = pub_cell%a3%z*0.529177249
  do i=-1,1
  do j=-1,1
  do k=-1,1
  v=  dble(i)*real_lattice(1,:) 
  v=v+dble(j)*real_lattice(2,:) 
  v=v+dble(k)*real_lattice(3,:)
  dist(i,j,k)=norm_vector(coordinates_atom(hatb)+v-coordinates_atom(hat))
  enddo
  enddo
  enddo 
  distance=minval(dist(:,:,:))
  uu=minloc(dist(:,:,:))
  uu=uu-2
  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real(8) function myfrequ(energy)
      implicit none
      complex(8) :: energy
       if (pub_dmft_temp .gt. 0.0_DP) then
         myfrequ=abs(aimag(energy))
       else
         myfrequ=abs(real(energy))
       endif
      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !--------------!

      pure integer function totn_atom_merge(hub_atom)
      implicit none
      integer,intent(in) :: hub_atom
      integer :: hat,totn,j
      integer :: vv(0:nmerge)
      if(cluster)then
       vv=totn_atom_merge_detail(hub_atom)
       totn_atom_merge=sum(vv(1:))
      else
       hat = pub_hub_atoms_on_node(hub_atom,pub_my_node_id)
       totn_atom_merge=hub_atom_hub_projs(hat)
      endif
      end function

            !--------------!

      pure integer function totn_atom_merge_split(hub_atom)
      implicit none
      integer,intent(in) :: hub_atom
      integer            :: hat,totn,j
      integer            :: vv(0:nmerge)
      if(cluster)then
       vv=totn_atom_merge_detail_split(hub_atom)
       totn_atom_merge_split=sum(vv(1:))
      else
       hat = hub_atom
       totn_atom_merge_split=hub_atom_hub_projs(hat)
      endif
      end function

            !--------------!

      pure function totn_atom_merge_detail(hub_atom)
      implicit none
      integer,intent(in) :: hub_atom
      integer            :: kk,hat,hat0,totn,j,totn_atom_merge_detail(0:nmerge)
      totn_atom_merge_detail=0
      totn_atom_merge_detail(0:0)=0
      hat0 = pub_hub_atoms_on_node(hub_atom,pub_my_node_id)
      do kk=1,merge_neigh(hat0)
         j=merge_table(hat0,kk)
         if(j>0)then
           totn_atom_merge_detail(kk)=hub_atom_hub_projs(j)
         endif
      enddo

      end function

          !--------------!

      pure function totn_atom_merge_detail_split(hub_atom)
      implicit none
      integer,intent(in) :: hub_atom
      integer            :: kk,hat,hat0,totn,j,totn_atom_merge_detail_split(0:nmerge)
      totn_atom_merge_detail_split=0
      totn_atom_merge_detail_split(0:0)=0
      hat0 = hub_atom
      do kk=1,merge_neigh(hat0)
         j=merge_table(hat0,kk)
         if(j>0)then
           totn_atom_merge_detail_split(kk)=hub_atom_hub_projs(j)
         endif
      enddo

      end function

            !--------------!

      function merged_green_block(hub_atom,hub_atomb,n,m)
      implicit none
      integer    :: n,m,hat,hatb
      integer    :: s1,s2,kk,kkb,hub_atom,hub_atomb,totn,totnb,j,jb,v(0:nmerge),vb(0:nmerge)
      complex(8) :: merged_green_block(n,m) 

      if(.not.cluster)then
        j =totn_atom_merge(hub_atom) 
        jb=totn_atom_merge(hub_atomb)
        merged_green_block(1:j,1:jb) = site_greenf_buffer(hub_atom,hub_atomb,1:j,1:jb) 
        return
      endif

      if(totn_atom_merge(hub_atom)/=totn_atom_merge(hub_atomb))then
       write(*,*) 'ERROR  sorry extended merging scheme not yet done for rectangular subblocks'
       stop
      endif
      merged_green_block=0.d0
      v =totn_atom_merge_detail(hub_atom)
      vb=totn_atom_merge_detail(hub_atomb)
      hat = pub_hub_atoms_on_node(hub_atom,pub_my_node_id)
      hatb= pub_hub_atoms_on_node(hub_atomb,pub_my_node_id)
      do kk =1,merge_neigh(hat)
      do kkb=1,merge_neigh(hatb)
         j  =merge_table(hat,kk)
         jb =merge_table(hatb,kkb)
         s1=sum( v (0:kk-1 ) )
         s2=sum( vb(0:kkb-1) )
        if(verboseall)then
         write(*,*) 'kk (1-nneigh) : ',kk, merge_neigh(hat)
         write(*,*) 'kkb           : ',kkb,merge_neigh(hatb)
         write(*,*) 's1-s2         : ',s1,s2
         write(*,*) 'vkk,vkk2      : ',v(kk),vb(kkb)
         write(*,*) 'j,jb          : ',j,jb
        endif
         merged_green_block(s1+1:s1+v(kk),s2+1:s2+vb(kkb)) = site_greenf_buffer(hubard_on_my_node(j),hubard_on_my_node(jb),1:v(kk),1:vb(kkb)) 
      enddo
      enddo

      return
      end function

            !--------------!

      function merged_green_block_split(hub_atom,hub_atomb,n,m)
      implicit none
      integer    :: n,m,hat,hatb
      integer    :: s1,s2,kk,kkb,hub_atom,hub_atomb,totn,totnb,j,jb,v(0:nmerge),vb(0:nmerge)
      complex(8) :: merged_green_block_split(n,m)

      if(.not.cluster)then
        j =totn_atom_merge_split(hub_atom)
        jb=totn_atom_merge_split(hub_atomb)
        merged_green_block_split(1:j,1:jb) = site_greenf_buffer(hub_atom,hub_atomb,1:j,1:jb)
        return
      endif
 
      if(totn_atom_merge_split(hub_atom)/=totn_atom_merge_split(hub_atomb))then
       write(*,*) 'ERROR  sorry extended merging scheme not yet done for rectangular subblocks'
       stop
      endif
      merged_green_block_split=0.d0
      v =totn_atom_merge_detail_split(hub_atom)
      vb=totn_atom_merge_detail_split(hub_atomb)
      hat = hub_atom
      hatb= hub_atomb
      do kk =1,merge_neigh(hat)
      do kkb=1,merge_neigh(hatb)
         j  =merge_table(hat,kk)
         jb =merge_table(hatb,kkb)
         s1=sum( v (0:kk-1 ) )
         s2=sum( vb(0:kkb-1) )
        if(verboseall)then
         write(*,*) 'kk (1-nneigh) : ',kk, merge_neigh(hat)
         write(*,*) 'kkb           : ',kkb,merge_neigh(hatb)
         write(*,*) 's1-s2         : ',s1,s2
         write(*,*) 'vkk,vkk2      : ',v(kk),vb(kkb)
         write(*,*) 'j,jb          : ',j,jb
        endif
         merged_green_block_split(s1+1:s1+v(kk),s2+1:s2+vb(kkb)) = site_greenf_buffer(j,jb,1:v(kk),1:vb(kkb))
      enddo
      enddo

      return
      end function

            !--------------!
            !--------------!

      subroutine distribute(big,hub_atom,hub_atomb,n_,m_)
      implicit none
      integer     ::  n_,m_
      integer     ::  s1,s2,kk,kkb,hub_atom,hub_atomb,totn,totnb,j,jb,v(0:nmerge),vb(0:nmerge)
      complex(8)  ::  big(n_,m_)  

      if(.not.cluster)then
        j =totn_atom_merge(hub_atom)
        jb=totn_atom_merge(hub_atomb)
        site_self_energy_buffer(hub_atom,hub_atomb,1:j,1:jb) = big(1:j,1:jb)
        return
      endif

      if(totn_atom_merge(hub_atom)/=totn_atom_merge(hub_atomb))then
        write(*,*) 'ERROR  sorry extended merging scheme not yet done for rectangular subblocks'
        stop
      endif

      v    = totn_atom_merge_detail(hub_atom)
      vb   = totn_atom_merge_detail(hub_atomb)
      hat  = pub_hub_atoms_on_node(hub_atom,pub_my_node_id)
      hatb = pub_hub_atoms_on_node(hub_atomb,pub_my_node_id)

      do kk =1,merge_neigh(hat)
       do kkb=1,merge_neigh(hatb)
         j  = merge_table(hat,kk)
         jb = merge_table(hatb,kkb)
         s1 = sum( v (0:kk-1 ) )
         s2 = sum( vb(0:kkb-1) )
         site_self_energy_buffer(hubard_on_my_node(j),hubard_on_my_node(jb),1:v(kk),1:vb(kkb)) = big(s1+1:s1+v(kk),s2+1:s2+vb(kkb)) 
       enddo
      enddo

      end subroutine

            !--------------!

      function merged_occupancy_block(hub_atom,spin,n_)
      implicit none
      integer    :: n_
      integer    :: s1,kk,hub_atom,totn,j,v(0:nmerge),spin
      real(8)    :: merged_occupancy_block(n_,n_)

      if(.not.cluster)then
        j   = totn_atom_merge(hub_atom)
        hat = pub_hub_atoms_on_node(hub_atom,pub_my_node_id)
        merged_occupancy_block(1:j,1:j) = diago_occupancy_func(h_atoms(hat)%occupancy(1:j,1:j,spin),j, &
                                          & pub_dmft_rotate_green,hat,rot_vec_angles(hub_atom,1:3,1:3))
        return
      endif

      merged_occupancy_block=0.d0
      v =totn_atom_merge_detail(hub_atom)
      hat=pub_hub_atoms_on_node(hub_atom,pub_my_node_id)
      do kk =1,merge_neigh(hat)
         j = merge_table(hat,kk)
         s1=sum( v(0:kk-1) )
         merged_occupancy_block(s1+1:s1+v(kk),s1+1:s1+v(kk)) = &
       & diago_occupancy_func(h_atoms(j)%occupancy(1:v(kk),1:v(kk),spin),v(kk),pub_dmft_rotate_green,hat,rot_vec_angles(hubard_on_my_node(j),1:3,1:3))
      enddo

      return
      end function

         !--------------!
         !--------------!

      function merged_occupancy_block_split(hub_atom,spin,n_)
      implicit none
      integer    :: n_
      integer    :: s1,kk,hub_atom,totn,j,v(0:nmerge),spin
      real(8)    :: merged_occupancy_block_split(n_,n_)

      if(.not.cluster)then
        j   = totn_atom_merge_split(hub_atom)
        hat = hub_atom
        merged_occupancy_block_split(1:j,1:j) = diago_occupancy_func(h_atoms(hat)%occupancy(1:j,1:j,spin),j, &
                                          & pub_dmft_rotate_green,hat,rot_vec_angles_split(hub_atom,1:3,1:3))
        return
      endif

      merged_occupancy_block_split=0.d0
      v =totn_atom_merge_detail_split(hub_atom)
      hat=hub_atom
      do kk =1,merge_neigh(hat)
         j = merge_table(hat,kk)
         s1=sum( v(0:kk-1) )
         merged_occupancy_block_split(s1+1:s1+v(kk),s1+1:s1+v(kk)) = &
       & diago_occupancy_func(h_atoms(j)%occupancy(1:v(kk),1:v(kk),spin),v(kk),pub_dmft_rotate_green,hat,rot_vec_angles_split(j,1:3,1:3))
      enddo

      return
      end function

            !--------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine set_loop_variable(hub_atom,hub_atomb)
      implicit none
      integer :: hub_atom,hub_atomb
             hat       =  pub_hub_atoms_on_node(hub_atom,pub_my_node_id)
             theatom   =  pub_distr_atom( hub%orig(hat)  )
             channels  =  hub_atom_hub_projs(hat)
             sp        =  hub%species_number(hat)
             hatb      =  pub_hub_atoms_on_node(hub_atomb,pub_my_node_id)
             theatomb  =  pub_distr_atom( hub%orig(hatb) )
             channelsb =  hub_atom_hub_projs(hatb)
             spb       =  hub%species_number(hatb)
      end subroutine

      subroutine set_loop_variable_split(hub_atom,hub_atomb)
      implicit none
      integer :: hub_atom,hub_atomb
             hat       =  hub_atom
             theatom   =  pub_distr_atom( hub%orig(hat)  )
             channels  =  hub_atom_hub_projs(hat)
             sp        =  hub%species_number(hat)
             hatb      =  hub_atomb
             theatomb  =  pub_distr_atom( hub%orig(hatb) )
             channelsb =  hub_atom_hub_projs(hatb)
             spb       =  hub%species_number(hatb)
      end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine dump_info1
      implicit none
           if(silent) return

           if(ien==1)then
                write(*,*) '----------------------------------------------------'
                write(*,*) ' RANK [x] is doing ATOM [x] : ', pub_my_node_id , hat
                write(*,*) ' HUBBARD ATOM LABEL         : ', pub_distr_atom( hub%orig(hat) )
                write(*,*) ' size of h_atoms            : ', size(h_atoms(:))                
                write(*,*) ' hub%orig(hat)              : ', hub%orig(hat)
                write(*,*) '----> the atom, hat, sp, is : ', theatom,hat,sp,is
                if(pub_my_node_id==0)then
                write(*,'(a,3f11.4)') '----> position   : ', coordinates_atom(hat)
                endif
                write(*,*) 'first atom , num atom       : ', pub_first_atom_on_node(pub_my_node_id),pub_num_atoms_on_node(pub_my_node_id)
                write(*,*) 'pub orig atom hat           : ', pub_orig_atom(hat)
                write(*,*) '----------------------------------------------------'
           endif

           if(verboseall)then
                write(*,*) '=========================='
                write(*,*) 'number of hubbard atom species : ', sp
                write(*,*) '=========================='
           endif

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 logical function doubleprec(energy,ien)
 implicit none
 complex(8) :: energy
 integer    :: ien
 real(8)    :: cchhem

 cchhem=fermi_e(is)+pub_dmft_chem_shift 

 doubleprec = pub_dmft_temp>0. .and. (myfrequ(energy)<pub_dmft_cutoff_small .or. &
                                    & myfrequ(energy)>pub_dmft_cutoff_tail) 
 doubleprec = doubleprec .or. ien>=pub_dmft_points-nspecial_frequ+1 .or. &
                           & (ien==1.and.pub_dmft_temp<0.0_DP)

 doubleprec = doubleprec .or. (pub_dmft_temp<0. .and.   &
                            & (  abs(real(energy)-cchhem)<pub_dmft_cutoff_small .or. &
                            &    abs(real(energy)-cchhem)>pub_dmft_cutoff_tail ) )
 end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine invert_matc(i,mat)
 implicit none
 integer    :: i
 complex(8) :: mat(:,:)

      if(i/=size(mat,1).or.i/=size(mat,2))then
       write(*,*) 'ERROR : invert_matc - dim do not match'
       stop
      endif

#ifdef GPU_SPEEDUP
      if(use_gpu_onlyme)then
#ifdef GPU_SPEEDUP_SINGPREC
       call magma_fortran_comp_singleprec_bypass(i,mat)
#else
#ifdef GPU_HANGS
       call matinv_ge_complex(i,mat)
#else
       call magma_fortran_comp_(i,mat)
#endif
#endif
     else
       call invert_gen_cmat(i,mat)
     endif
#else
#ifdef debug
       if(maxval(abs(mat-transpose(mat)))>1.d-5)then
         write(*,*) 'ERROR invert_matc : matrix NOT SYMMETRIC'
       endif
#endif
       call invert_gen_cmat(i,mat)
#endif

 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine invert_mat_sparsec(mat,mat_inv,iter,greenfunc)
 use sparse,   only : sparse_convert,sparse_num_rows,sparse_num_cols
 implicit none
 integer                                        ::  nn,nnn,iter
 TYPE(SPAM3)                                    ::  mat
 TYPE(SPAM3)                                    ::  mat_inv
 complex(kind=DP),allocatable,dimension(:,:)    ::  mat_square
 logical                                        ::  greenfunc

#ifndef GPU_SPEEDUP
#ifndef INVERT_N_CUBE
      if(pub_total_num_nodes>1)Then
       call internal_invert_full(mat_inv,mat,greenfunc)
       return
      endif
#endif
#endif

#ifdef GPU_SPEEDUP
#ifndef INVERT_N_CUBE
      if(.not.(use_gpu.or.use_gpu_some_only))then
      if(pub_total_num_nodes>1)Then
        call internal_invert_full(mat_inv,mat,greenfunc)
        return
      endif
      endif
#endif
#endif

       nn=NINT(sparse_num_rows(mat))
      nnn=NINT(sparse_num_rows(mat_inv))

      if(nn/=nnn)then
        write(*,*) 'ERROR invert_mat_sparsec dim dont match'
        stop
      endif

      allocate(mat_square(nn,nn))

      if( mpi_rank/=1 .and. use_gpu_some_only .and. .not. breaknfscom)then
         mat_square=0.d0 ; goto 3131
      endif
   
      call sparse_convert(mat_square,mat)

#ifdef GPU_SPEEDUP
      if(use_gpu_onlyme)then
#ifdef GPU_SPEEDUP_SINGPREC
        call magma_fortran_comp_singleprec_bypass(nn,mat_square)
#else
#ifdef GPU_HANGS
        call matinv_ge_complex(nn,mat_square)
#else
        call magma_fortran_comp_(nn,mat_square)
#endif
#endif
     else
        call invert_gen_cmat(nn,mat_square)
     endif
#else
        call invert_gen_cmat(nn,mat_square)
#endif

      3131 continue

      if(mpi_size>1.and.use_gpu_some_only.and. .not. breaknfscom)then
        call sync_nodes_by_shared_nfs_matc(mat_square,'invert_matc'//adjustl(trim(tostring(iter))),mpi_rank,mpi_size,nn)
      endif

      if(greenfunc)then
       greenfbackupflag=.true.
       greenfbackup=mat_square
      endif

      call sparse_convert(mat_inv,mat_square)
      deallocate(mat_square)

 end subroutine

          !----------------------------------!
          !----------------------------------!
 
 subroutine sparse_max_diff( t1, t2 )
 use sparse,   only : sparse_convert,sparse_num_rows,sparse_num_cols
 implicit none
 TYPE(SPAM3)            :: t1,t2
 complex(8),allocatable :: t1_(:,:),t2_(:,:)
 integer                :: n1,n2
   n1=NINT(sparse_num_rows(t1))
   n2=NINT(sparse_num_cols(t2))
   allocate(t1_(n1,n2),t2_(n1,n2))
   call sparse_convert(t1_,t1)
   call sparse_convert(t2_,t2)
   write(*,*) 'MAX DIFFERENCE BETWEEN SPARSE IS : ', maxval(abs(t1_-t2_))
   deallocate(t1_,t2_)
 end subroutine

          !----------------------------------!
          !----------------------------------!

 subroutine invert_mat_sparser(mat,mat_inv,iter,greenfunc)
 use sparse,   only : sparse_convert,sparse_num_rows,sparse_num_cols
 implicit none
 integer                                        ::  nn,nnn,iter
 TYPE(SPAM3)                                    ::  mat
 TYPE(SPAM3)                                    ::  mat_inv
 real(kind=DP),allocatable,dimension(:,:)       ::  mat_square
 logical                                        ::  greenfunc

 if(greenfunc)then
  write(*,*) 'ERROR invert_mat_sparser not supposed to compte GF' ; stop
 endif


#ifndef GPU_SPEEDUP
#ifndef INVERT_N_CUBE
      if(pub_total_num_nodes>1)Then
       call internal_invert_full_real(mat,mat_inv,greenfunc)
       return
      endif
#endif
#endif

#ifdef GPU_SPEEDUP
#ifndef INVERT_N_CUBE
      if(.not.(use_gpu.or.use_gpu_some_only))then
      if(pub_total_num_nodes>1)Then
        call internal_invert_full_real(mat,mat_inv,greenfunc)
        return
      endif
      endif
#endif
#endif

      nn=NINT(sparse_num_rows(mat))
      nnn=NINT(sparse_num_rows(mat_inv))

      if(nn/=nnn)then
        write(*,*) 'ERROR invert_mat_sparsec dim dont match'
        stop
      endif

      allocate(mat_square(nn,nn))

      if( mpi_rank/=1 .and. use_gpu_some_only .and. .not. breaknfscom)then
         mat_square=0.d0 ; goto 3131
      endif

      call sparse_convert(mat_square,mat)

#ifdef GPU_SPEEDUP
      if(use_gpu_onlyme)then
#ifdef GPU_SPEEDUP_SINGPREC
        call magma_fortran_double_singleprec_bypass(nn,mat_square)
#else
        call magma_fortran_double_(nn,mat_square)
#endif
     else
        call invert_gen_rmat(nn,mat_square)
     endif
#else
        call invert_gen_rmat(nn,mat_square)
#endif

      3131 continue

      if(mpi_size>1.and.use_gpu_some_only.and. .not. breaknfscom)then
        call sync_nodes_by_shared_nfs_mat(mat_square,'invert_matr'//adjustl(trim(tostring(iter))),mpi_rank,mpi_size,nn)
      endif

      call sparse_convert(mat_inv,mat_square)
      deallocate(mat_square)

 end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
