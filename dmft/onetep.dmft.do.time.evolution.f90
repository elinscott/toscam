program keldysh
use linalg
use StringManip , only : StrInt2,toString
implicit none
integer                :: i,j,k,l,ispin,ntime,Nc
real(8),allocatable    :: n_v_t_matrix(:,:,:,:)
complex(8),allocatable :: sigma(:,:,:,:)
integer                :: hat,species,cpus
character(60)          :: mpiline

  write(*,*) 'CHECK THAT REAL_SPACE_PLOT IN RUN.DAT IS SWITCHED ON'
  write(*,*) 'TO SPEED UP CALCULATIONS YOU CAN ALSO SET pub_dmft_plot_real_space_sigma : T '
  write(*,*)
  write(*,*) 'HOW MANY CPUS'
  read(*,*)   cpus
  call system(" get_mpi_env ")
  open(unit=100,file='dir_onetep_mpi')
  read(100,'(a)') mpiline
  close(100)
  write(*,*) 'COMMAND IS : '
  write(*,*) mpiline
  write(*,*) 'PLEASE ENTER HUBARD ATOM NUMBER'
  read(*,*) hat
  write(*,*) 'PLEASE ENTER SPECIES'
  read(*,*) species

  do ispin=1,2
   open(unit=91,file='n_time_matrix_'//trim(adjustl(toString(ispin))),form='unformatted')
   read(91) ntime,Nc
   if(.not.allocated(n_v_t_matrix)) allocate(n_v_t_matrix(ntime,Nc,Nc,2))
   if(.not.allocated(sigma)) allocate(sigma(ntime,Nc,Nc,2))
   do i=1,ntime
     read(91) n_v_t_matrix(i,:,:,ispin)
     sigma(i,:,:,ispin)=-imi*n_v_t_matrix(i,:,:,ispin)
   enddo
   close(91)
  enddo

  call system(" mkdir MOVIE ")

  do i=1,ntime

   do ispin=1,2
    open(unit=91,file='sigma_output_'//trim(adjustl(toString(hat)))//'_'//trim(adjustl(toString(species)))//'_'//trim(adjustl(toString(ispin))),form='unformatted')
    do j=1,10
     write(91) sigma(i,:,:,ispin)
    enddo
    close(91)
   enddo

   write(*,*) 'RUNNING ONETEP TIME FRAME : ', i
   write(*,*) trim(adjustl(mpiline))//" -np "//trim(adjustl(tostring(cpus)))//" ed.biggie.cpu run "
   call system(trim(adjustl(mpiline))//" -np "//trim(adjustl(tostring(cpus)))//" ed.biggie.cpu run " )
   call system(" mv runscattering_fermi_level_dmft*.cube MOVIE/frame"//trim(adjustl(tostring(i)))//".cube")

  enddo

end program

