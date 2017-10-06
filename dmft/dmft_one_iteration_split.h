!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!

 subroutine initialize_my_simulation_
 use genvar
 implicit none
   if(.not.no_mpi)then
    if(reboot_mpi)then
     call MPI_INIT(ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,size2,ierr)
     call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    endif
   else
    ierr=0
    size2=1
    rank=0
   endif
   myid                    = rank
   numprocs                = size2
   nproc                   = size2
   iproc                   = rank  + 1
   messages                = .true.
   messages2               = .true.
   messages3               = .false.
   messages4               = .false.
   MPIseparate             = .false.
   testing                 = .false.
   strongstop              = .false.
   enable_mpi_dot          = .false.
 end subroutine

!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!

 subroutine loop_over_atoms
 implicit none 
 integer :: i,j
      if(.not.uniform_sigma)then 
        do i=rank+1,ww,size2 
          j=get_atom_number(TRIM(ADJUSTL(list_of_greens(i,1))))
          call one_atom_one_dmft_iter(i,j)
        enddo
        call mpibarrier
     else
        if(rank==0) then
          j=get_atom_number(TRIM(ADJUSTL(list_of_greens(1,1))))
          call one_atom_one_dmft_iter(1,j)         
        endif
        call mpibarrier
     endif 
 end subroutine

!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!

 subroutine copy_files

    if(uniform_sigma)then
       write(*,*) 'copying sigmas...WW=',ww
       if(ww>=2.and.rank==0) then
          call copy_sigmas
       endif
     else
       if(check.and.rank==0)then
          write(*,*) 'copying sigma pattern to other sites'
          call copy_sigmas_pattern
       endif
     endif
 
     if(rank==0) then
      call system("dmft_dimer_to_orbitals_script.out > onetep_projecting_back"//TRIM(ADJUSTL(toString(iter_dmft))))
     endif

 end subroutine

!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!

 subroutine build_list_of_files
 implicit none
 integer :: i,j

     INQUIRE(file='mask_uniform',EXIST=check)
     check=check.and..not.uniform_sigma

     !-------------------------------------------------------------!
     if(check) then
        open(unit=10001,file='mask_uniform')
        i=0
        do
          read(10001,*,end=4545)
          i=i+1
        enddo
        4545 continue
        rewind(10001)
        if(allocated(list_uniform)) deallocate(list_uniform)
        allocate(list_uniform(i,2)) 
        do i=1,size(list_uniform,1)
          read(10001,*) (list_uniform(i,j),j=1,2) 
          write(*,*) 'READING LIST UNIFORM [i],list1,list2: ', i,list_uniform(i,1:2)
        enddo 
        close(10001) 
        ww=0
        do i=1,size(list_uniform,1)
         if(abs(list_uniform(i,1))==abs(list_uniform(i,2))) ww=ww+1
        enddo
        write(*,*) 'THERE ARE [x] SITES IN UNIT-CELL : ',ww

        flag_onetep_producing_only_up_spin=onetep_only_up
        if(.not.allocated(list_of_greens)) allocate(list_of_greens(ww,2)) 

        i=0
        do j=1,size(list_uniform,1)
         if(abs(list_uniform(j,1))==abs(list_uniform(j,2)))then
          i=i+1
                                                         list_of_greens(i,1)='green_output'//trim(adjustl(toString( abs(list_uniform(j,1)) )))//'_1_1'
             if(.not.flag_onetep_producing_only_up_spin) list_of_greens(i,2)='green_output'//trim(adjustl(toString( abs(list_uniform(j,1)) )))//'_1_2'

          if(rank==0)                                             write(*,*) 'green func atom : ',i,TRIM(ADJUSTL(list_of_greens(i,1)))
          if(rank==0.and..not.flag_onetep_producing_only_up_spin) write(*,*) 'green func atom : ',i,TRIM(ADJUSTL(list_of_greens(i,2)))
         endif
        enddo
     !-------------------------------------------------------------!
     else
     !-------------------------------------------------------------!
         write(*,*) 'reading GREEN'
         open(unit=10001,file='list_of_greens') 
         read(10001,*) ww 
                 flag_onetep_producing_only_up_spin=mod(ww,2)/=0.or.onetep_only_up
         if(.not.flag_onetep_producing_only_up_spin) ww=ww/2
         if(.not.allocated(list_of_greens)) allocate(list_of_greens(ww,2)) 
         do i=1,ww
          read(10001,*) list_of_greens(i,1)
          if(.not.flag_onetep_producing_only_up_spin) read(10001,*) list_of_greens(i,2)
          if(rank==0) write(*,*) 'green func atom : ', i, TRIM(ADJUSTL(list_of_greens(i,1)))
          if(rank==0.and..not.flag_onetep_producing_only_up_spin) write(*,*) 'green func atom : ',i,TRIM(ADJUSTL(list_of_greens(i,2)))
         enddo
         close(10001)
     endif
     !-------------------------------------------------------------!

     if(paramagnetic/=1.and.flag_onetep_producing_only_up_spin)then
      write(*,*) 'oups..onetep only produces up spin and we dont want the paramagnetic calculations'
      do i=1,ww
        cc_=trim(adjustl(list_of_greens(i,1)))
        call replace_in_string(cc_,'_1','_2','last')
        list_of_greens(i,2)=cc_
        call system(" cp "//ADJUSTL(TRIM(list_of_greens(i,1)))//" "//ADJUSTL(TRIM(list_of_greens(i,2)))) 
      enddo
     endif 
 
 end subroutine


!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!

 subroutine copy_sigmas_pattern
 implicit none
 integer :: kkk_
 integer :: j

   if(paramagnetic==1)then
     kkk_=1
   else
     kkk_=2
   endif

   do j=1,size(list_uniform,1)
      if(abs(list_uniform(j,1))/=abs(list_uniform(j,2)))then

       if(list_uniform(j,2)>0)then
         call system(" cp  "//'sigma_output'//trim(adjustl(toString( abs(list_uniform(j,2)) )))//'_1_1'//"  "//&
                            & 'sigma_output'//trim(adjustl(toString( abs(list_uniform(j,1)) )))//'_1_1')
        if(kkk_==2)then
         call system(" cp  "//'sigma_output'//trim(adjustl(toString( abs(list_uniform(j,2)) )))//'_1_2'//"  "//&
                            & 'sigma_output'//trim(adjustl(toString( abs(list_uniform(j,1)) )))//'_1_2')
        endif
       endif

       if(list_uniform(j,2)<0)then ! if negative entry in mask_uniform, flip the spin !
         call system(" cp  "//'sigma_output'//trim(adjustl(toString( abs(list_uniform(j,2)) )))//'_1_2'//"  "//&
                            & 'sigma_output'//trim(adjustl(toString( abs(list_uniform(j,1)) )))//'_1_1')
        if(kkk_==2)then
         call system(" cp  "//'sigma_output'//trim(adjustl(toString( abs(list_uniform(j,2)) )))//'_1_1'//"  "//&
                            & 'sigma_output'//trim(adjustl(toString( abs(list_uniform(j,1)) )))//'_1_2')
        endif
       endif

      endif
   enddo

   do j=1,size(list_uniform,1)
      if(abs(list_uniform(j,1))/=abs(list_uniform(j,2)))then

       if(list_uniform(j,2)>0)then
         call system(" cp  "//'real_sigma_output'//trim(adjustl(toString( abs(list_uniform(j,2)) )))//'_1_1'//"  "//&
                            & 'real_sigma_output'//trim(adjustl(toString( abs(list_uniform(j,1)) )))//'_1_1'//&
                            &"   > /dev/null 2>&1 ")
        if(kkk_==2)then
         call system(" cp  "//'real_sigma_output'//trim(adjustl(toString( abs(list_uniform(j,2)) )))//'_1_2'//"  "//&
                            & 'real_sigma_output'//trim(adjustl(toString( abs(list_uniform(j,1)) )))//'_1_2'//&
                            &"   > /dev/null 2>&1 ")
        endif
       endif

       if(list_uniform(j,2)<0)then ! if negative entry in mask_uniform, flip the spin !
         call system(" cp  "//'real_sigma_output'//trim(adjustl(toString( abs(list_uniform(j,2)) )))//'_1_2'//"  "//&
                            & 'real_sigma_output'//trim(adjustl(toString( abs(list_uniform(j,1)) )))//'_1_1'//&
                            &"   > /dev/null 2>&1 ")
        if(kkk_==2)then
         call system(" cp  "//'real_sigma_output'//trim(adjustl(toString( abs(list_uniform(j,2)) )))//'_1_1'//"  "//&
                            & 'real_sigma_output'//trim(adjustl(toString( abs(list_uniform(j,1)) )))//'_1_2'//&
                            &"   > /dev/null 2>&1 ")
        endif
       endif

      endif
   enddo

 end subroutine

!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!

 subroutine copy_sigmas
 implicit none
 integer :: kkk_

   if(paramagnetic==1)then
     kkk_=1
   else
     kkk_=2
   endif

   write(*,*) 'copy sigmas, kkk_=', kkk_

   do j=1,kkk_
    filename_sigma_source=trim(adjustl(list_of_greens(1,j)))
    cc_=trim(Adjustl(filename_sigma_source)); call replace_in_string(cc_,'green','sigma','first'); filename_sigma_source=cc_
    do i=2,ww
     filename_sigma=trim(adjustl(list_of_greens(i,j)))
     cc_=trim(adjustl(filename_sigma)); call replace_in_string(cc_,'green','sigma','first'); filename_sigma=cc_
     write(*,*) 'UNIFORM SIGMA,copying : ',trim(adjustl(filename_sigma_source)),' TO : ', trim(adjustl(filename_sigma))
     call system("cp      "//trim(adjustl(filename_sigma_source))//"      "//trim(adjustl(filename_sigma)))
     call system("cp real_"//trim(adjustl(filename_sigma_source))//" real_"//trim(adjustl(filename_sigma)))
    enddo
   enddo

 end subroutine

      !-------------!
      !-------------!

 subroutine one_atom_one_dmft_iter(i,jj)
 implicit none
 integer        :: i,jj
 character(400) :: command_line

  if(paramagnetic==1)then
    write(*,*) 'imposing paramagnet'
    command_line="onetep_dmft_one_iter_script_pm.out "//TRIM(ADJUSTL(toString(jj)))//" "//TRIM(ADJUSTL(list_of_greens(i,1)))//" "//TRIM(ADJUSTL(toString(iter_dmft)))
    command_line=trim(adjustl(command_line))//" "//trim(adjustl(input_temp_dir))//" > onetep_dmft_part_solver_"//TRIM(ADJUSTL(toString(iter_dmft)))//"_atom_"//TRIM(ADJUSTL(toString(jj)))
    write(*,*) 'calling command line : ', TRIM(ADJUSTL(command_line))
    call system(command_line)
  else
    write(*,*) 'no paramagnetic constraint'
    write(*,*) 'green func files : ', TRIM(ADJUSTL(list_of_greens(i,1))), TRIM(ADJUSTL(list_of_greens(i,2)))
    command_line="onetep_dmft_one_iter_script_polarized.out "//TRIM(ADJUSTL(toString(jj)))//" "//TRIM(ADJUSTL(list_of_greens(i,1)))//" "//TRIM(ADJUSTL(list_of_greens(i,2)))&
                &//" "//TRIM(ADJUSTL(toString(iter_dmft)))
    command_line=trim(adjustl(command_line))//" "//trim(adjustl(input_temp_dir))//" > onetep_dmft_part_solver_"//TRIM(ADJUSTL(toString(iter_dmft)))//"_atom_"//TRIM(ADJUSTL(toString(jj)))
    write(*,*) 'calling command line : ', TRIM(ADJUSTL(command_line))
    call system(command_line)
  endif

  if(flag_donot_keep_all_files)then
     call system("rm -r atom"//trim(adjustl(toString(jj)))//"/*/AGR")
     call system("rm -f atom"//trim(adjustl(toString(jj)))//"/*/pgplot.ps")
     call system("rm -f atom"//trim(adjustl(toString(jj)))//"/*/chiloc_vertex*")
  endif

  if(iter_dmft==niter_dmft)then
    call system("check_if_symlink atom"//trim(adjustl(toString(jj))))
  endif

 end subroutine

!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!

integer function get_atom_number(filename)
implicit none
integer       :: j,i
character*(*) :: filename
  do i=LEN(filename),1,-1
   if(filename(i:i)=="t")then
     exit
   endif
  enddo
  i=i+1
  do j=i,2000
   if(filename(j:j)=="_")then
     exit
   endif
  enddo
  j=j-1
  write(*,*) 'ATOM NUMBER : ', filename(i:j)
  get_atom_number=StrInt2(TRIM(ADJUSTL(filename(i:j))))
end function

!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!
!------------------------------!
