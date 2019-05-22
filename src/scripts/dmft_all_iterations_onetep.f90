!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

 program onetep_normal_mode
    use StringManip, only: toString
    use mpi
    implicit none
    integer          :: jj, k, i, rank, size2, ierr
    integer          :: nproc_onetep
    integer          :: iter_dmft
    character(2000)  :: exec_onetep
    character(2000)  :: CASE_ONETEP
    logical          :: compute_dos
    real(8)          :: ed_frequ_min
    real(8)          :: ed_frequ_max
    integer          :: ed_real_frequ_last
    character(2000)  :: command_line
    character(200)   :: myhost

    open (unit=55, file='temp_onetep_part', form='unformatted')
    read (55) nproc_onetep, iter_dmft, exec_onetep, &
           & CASE_ONETEP, compute_dos, ed_frequ_min, ed_frequ_max, ed_real_frequ_last
    close (55)

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size2, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    write (*, *) 'STARTING MAIN LOOP ON CPUS'
    call system(" which "//TRIM(ADJUSTL(exec_onetep))//" && echo 'onetep mpi version exists' ")
        call system( " which "//TRIM(ADJUSTL(exec_onetep))//".serial || echo 'ERROR please compile onetep in serial mode and name it with .serial extension'")

    do i = 1, nproc_onetep
       if (rank == i - 1) then
          command_line = TRIM(ADJUSTL(exec_onetep))//".serial "//TRIM(ADJUSTL(CASE_ONETEP))//" "//&
                     & TRIM(ADJUSTL(toString(i)))//" "//TRIM(ADJUSTL(toString(nproc_onetep)))
          if (compute_dos) then
             command_line = TRIM(ADJUSTL(command_line))//" "//TRIM(ADJUSTL(toString(ed_real_frequ_last)))//" 1 "//&
                        & TRIM(ADJUSTL(toString(ed_frequ_min)))//"  "//TRIM(ADJUSTL(toString(ed_frequ_max)))//"  "
          endif
          command_line=TRIM(ADJUSTL(command_line))//" > onetep_output_iter"//TRIM(ADJUSTL(toString(iter_dmft)))//"_rank__"//TRIM(ADJUSTL(toString(i)))
          write (*, *) ' command line for proc [x], total : ', i, nproc_onetep
          write (*, *) TRIM(ADJUSTL(command_line))
          call system(command_line)
          write (*, *) ' ... done ... rank : ', rank
       endif
    enddo

    write (*, *) 'RANK [x] done : ', rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (rank == 0) call system(" rm temp_onetep_part ")
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    write (*, *) 'done finalize'
    call MPI_FINALIZE(ierr)

    write (*, *) 'EXITING ONETEP PART'
    stop

 end program

!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!
!==============================================================!

