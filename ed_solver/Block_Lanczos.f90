MODULE Block_Lanczos

  use ED_ARPACK

  IMPLICIT NONE


  REAL(DBL),    PARAMETER, PRIVATE                 :: zero=0.0_DBL,one=1.0_DBL,two=2.0_DBL,three=3.0_DBL,four=4.0_DBL
  LOGICAL,      PARAMETER, PRIVATE                 :: F=.FALSE.,T=.TRUE.

 ! DIAGONALIZATION STRATEGIES USING LANCZOS ALGORITHM TO GET *SEVERAL* EIGENVALUES

  contains 


!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************


  SUBROUTINE Block_Lanczos_diagonalize(lowest)

    TYPE(eigenlist_type), INTENT(INOUT) :: lowest    
    TYPE(eigen_type)                    :: eigen
    INTEGER                             :: Nconverged,Niter,ivp,start_diagH
    INTEGER                             :: local_istatemin(nproc),local_istatemax(nproc),local_chunk(nproc)
    REAL(DBL), ALLOCATABLE              :: vec_in(:),vec_out(:)
    INTEGER                             :: LISTOR,LRSTOR
    INTEGER                             :: ISTOR_QUERY(17)
    REAL(DBL)                           :: RSTOR_QUERY(5)
    INTEGER,   ALLOCATABLE              :: ISTOR(:)
    REAL(DBL), ALLOCATABLE              :: RSTOR(:),U(:,:),V(:,:),VALP(:,:),VECP(:,:)
    REAL(DBL)                           :: SIGMA,minv_
    INTEGER                             :: NNEIG,LFLAG,NVOPU,Neigen_

    CALL reset_timer(start_diagH)
    CALL split(dimen_H(),local_istatemin,local_istatemax,local_chunk)

    Neigen_=MIN(Neigen,dimen_H())

    write(*,*) '-------------Block Lanczos-------------'
    write(*,*) ' Block_size        : ', Block_size
    write(*,*) ' Nitermax          : ', Nitermax
    write(*,*) ' Neigen_           : ', Neigen_
    write(*,*) ' locistatemin max  : ', local_istatemin,local_istatemax
    write(*,*) ' rank,local_chunk  : ', rank,local_chunk(iproc)
    write(*,*) ' dimen_H           : ', dimen_H()
    write(*,*) '---------------------------------------'

    ISTOR_QUERY      =  0
    ISTOR_QUERY(1)   =  local_chunk(iproc)
    ISTOR_QUERY(2)   =  local_chunk(iproc)
    ISTOR_QUERY(3)   =  Neigen_
    ISTOR_QUERY(4)   =  MIN(Neigen_+10,Neigen_*3) ! >=MIN(Neigen_+10,Neigen_*2) RECOMMENDED
    ISTOR_QUERY(5)   =  Block_size
    ISTOR_QUERY(6)   =  min(Nitermax,dimen_H())
    ISTOR_QUERY(7)   =  Block_size          ! We will provide 'Block size' input vectors
    ISTOR_QUERY(8)   =  0                   ! no input eigenpairs
    ISTOR_QUERY(9)   =  0                   ! standard eigenvalue problem
    ISTOR_QUERY(10)  =  0                   ! ignored in standard ev case
    ISTOR_QUERY(11)  =  0                   ! ignored in standard ev case
    ISTOR_QUERY(12)  =  0                   ! verbose
    ISTOR_QUERY(13)  =  all_log_unit(1)     ! all output has to be redirected to the same file!
    if(.not.no_mpi)then
      ISTOR_QUERY(14)  =  MPI_COMM_WORLD      ! MPI communicator 
    else
      ISTOR_QUERY(14)  =  0
    endif
    ISTOR_QUERY(15)  =  0                   ! let BLZDRD determine the workspace
    
    RSTOR_QUERY(1)   = zero       ! ignored in standard ev case
    RSTOR_QUERY(2)   = zero       ! ignored in standard ev case
    RSTOR_QUERY(3)   = tolerance

    if(Block_size>0)then
     LISTOR=0
     LRSTOR=0                     ! let BLZDRD determine the workspace
    else
     LISTOR=200000
     LRSTOR=200000     
    endif
   
    ISTOR_QUERY(15) = LISTOR
    RSTOR_QUERY( 4) = LRSTOR
    
    SIGMA = zero   ! ignored in standard ev case
    NNEIG = 0      ! ignored in standard ev case
    LFLAG = 0      ! FIRST CALL
    NVOPU = 0      ! FIRST CALL
    
    write(*,*) '.... Block Lanczos allocate main arrays ..... '
    ALLOCATE(VALP(ISTOR_QUERY(4),2),VECP(ISTOR_QUERY(2),ISTOR_QUERY(4)))
 
    !*************************************************************************! 

        if(allocated(ISTOR)) deallocate(ISTOR)
        if(allocated(RSTOR)) deallocate(RSTOR)
        if(allocated(U))     deallocate(U)
        if(allocated(V))     deallocate(V)

        ALLOCATE(ISTOR(17+ISTOR_QUERY(15)))
 ALLOCATE(RSTOR(5+NINT(RSTOR_QUERY(4))))

       if(Block_size>0)then
         ALLOCATE(U(ISTOR_QUERY(2),Block_size),V(ISTOR_QUERY(2),Block_size))
       else
         ALLOCATE(U(ISTOR_QUERY(2),Neigen_),V(ISTOR_QUERY(2),Neigen_))
       endif  

       ISTOR(1:17) = ISTOR_QUERY
       RSTOR(1:5)  = RSTOR_QUERY
       write(*,*) 'CALL BLOCK LANCZOS LIBRARY'
#ifdef _blocklanczos
       CALL BLZDRD(ISTOR,RSTOR,SIGMA,NNEIG,U,V,LFLAG,NVOPU,VALP,VECP)
#else
       write(*,*) 'error recompile with -D_blocklanczos'
       stop
#endif
       ISTOR_QUERY = ISTOR(1:17)
       RSTOR_QUERY = RSTOR(1:5)
       write(*,*) 'after initialization, got flags : '
       write(*,*)  'LFLAG ----> ',LFLAG
       write(*,*)  'NVOPU ----> ',NVOPU

    !*************************************************************************!

    IF(Block_size==0)THEN
        write(*,*) ' WE LET BLZDRD DECIDE THE BLOCK SIZE '
        Block_size = NVOPU
        CALL dump_message(TEXT="# OPTIMAL BLOCK SIZE ADOPTED BY BLZDRD = "//c2s(i2c(Block_size)))
        write(*,*) ' block size adopted : ', Block_size
    ENDIF

    ISTOR_QUERY(5) = Block_size
    ISTOR_QUERY(7) = Block_size  ! NUMBER OF INPUT VECTORS THAT WE SPECIFY

    IF(LFLAG<0) then
       write(*,*) ' QUERY GAVE ERROR FLAG : ', LFLAG
       write(*,*) ' Block_size            : ', Block_size
       write(*,*) ' Neigen_                : ', Neigen_
       write(*,*) ' local chunk           : ', local_chunk(iproc)
       write(*,*) ' ISTOR QUERY           : ', ISTOR_QUERY
       write(*,*) ' RSTOR_QUERY           : ', RSTOR_QUERY(1:4)
       STOP "ERROR IN Block_Lanczos_diagonalize_real: BLZDRD WORKSPACE QUERY FAILED!"
    ENDIF
    
    write(*,*) '-----------------------------------------------'
    write(*,*) '         ISTOR_QUERY   : '
    write(*,*)           ISTOR_QUERY
    write(*,*) '----->   SIZE OF ISTOR : ', 17+ISTOR_QUERY(15)
    write(*,*) '----->   SIZE OF RSTOR : ', 5+NINT(RSTOR_QUERY(4))
    write(*,*) '-----------------------------------------------'

    if(allocated(ISTOR)) deallocate(ISTOR)
    if(allocated(RSTOR)) deallocate(RSTOR)
    if(allocated(U))     deallocate(U)
    if(allocated(V))     deallocate(V)

    ALLOCATE(ISTOR(17+ISTOR_QUERY(15)))
    ISTOR(1:17) = ISTOR_QUERY
    ALLOCATE(RSTOR(5+NINT(RSTOR_QUERY(4))))
    RSTOR(1:5)  = RSTOR_QUERY

    ALLOCATE(U(ISTOR(2),Block_size),V(ISTOR(2),Block_size))
    ALLOCATE(vec_in(dimen_H()),vec_out(dimen_H()))

    LFLAG = 0      ! FIRST CALL
    NVOPU = 0      ! FIRST CALL

    DO

     write(*,*) ' ---- iteration lanczos ---- '

     !**************************************************************!
#ifdef _blocklanczos
      CALL BLZDRD(ISTOR,RSTOR,SIGMA,NNEIG,U,V,LFLAG,NVOPU,VALP,VECP) 
#else
       write(*,*) 'error recompile with -D_blocklanczos'
       stop
#endif
     !**************************************************************!
      
     Nconverged  =  ISTOR(67)
     Niter       =  ISTOR(19)

     write(*,*) ' ------------------------------------------------------ '
     write(*,*) ' Nconverged eigenpairs / Niter ' , LFLAG, Nconverged, Niter
     write(*,*) ' ------------------------------------------------------ '
 
     IF(LFLAG<0) STOP "ERROR IN Block_Lanczos_diagonalize_real: BLZDRD EXITED WITH LFLAG<0!"

     SELECT CASE(LFLAG)

     CASE(0)

          EXIT

     CASE(1)    ! BLDRD REQUIRES V(:,1:NVOPU) = H * U(:,1:NVOPU)

       !------------------------------------------------------------------------!
       DO ivp=1,NVOPU
         vec_in = zero

         ! Hmult REQUIRES FULL VECTOR => GATHER ALL PARTS FROM LOCAL U(:,ivp)
         if(size2>1) then
           CALL MPI_ALLGATHERV(U(:,ivp),local_chunk(iproc),MPI_DOUBLE_PRECISION,          &
           vec_in,local_chunk,local_istatemin-1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
         else
           vec_in=U(:,ivp)
         endif

         ! APPLY |vec_out> = H |vec_in>
         CALL Hmult__(size(vec_in),vec_out,vec_in)

         ! BLZDRD REQUIRES LOCAL CHUNK ONLY => DISTRIBUTE vec_out OVER LOCAL V(:,ivp)
        V(:,ivp) = vec_out(local_istatemin(iproc):local_istatemax(iproc))
       ENDDO
       !------------------------------------------------------------------------!

     CASE(4)  ! BLDRD REQUIRES THE STARTING VECTORS IN V(:,1:NVOPU) 

       CALL randomize_mat(V)

     CASE DEFAULT

       STOP "ERROR IN Block_Lanczos_diagonalize_real: BLZDRD EXITED WITH UNEXPECTED LFLAG!"

     END SELECT

    ENDDO

    ! NOW WE HAVE THE DESIRED EIGENPAIRS THANKS TO BLZPACK IT JUST REMAINS TO STORE THEM IN EIGENLIST

    Nconverged = ISTOR(67)  ! NUMBER OF CONVERGED EIGENPAIRS (MAY BE MORE THAN Neigen_ REQUIRED!)
    Niter      = ISTOR(19)  ! NUMBER OF BLOCK LANCZOS ITERATIONS PERFORMED

    write(*,*) '===================================================================='    
    write(*,*) '  Nconverged eigenpairs  : ', rank,Nconverged, VALP(1:Nconverged,1)
    write(*,*) '  Niter                  : ', rank,Niter
    write(*,*) '===================================================================='

    IF(Nconverged<=0) then
      write(*,*) "ERROR IN Block_Lanczos_diagonalize_real: NO EIGENPAIR IS CONVERGED AFTER "//c2s(i2c(Niter))//" ITERATIONS!" 
      write(*,*) 'ISTOR: ' 
      write(*,*)  ISTOR
      stop
    endif


    minv_=minval(VALP(1:Nconverged,1))
    DO ivp=1,Nconverged
     if(abs(VALP(ivp,1)-minv_)<=dEmax)then
      CALL new_eigen(eigen,dimen_H())
      eigen%val       = VALP(ivp,1)
      eigen%converged = T
      eigen%rank      = ivp
      if(size2>1) then
        CALL MPI_ALLGATHERV(VECP(:,ivp),local_chunk(iproc),MPI_DOUBLE_PRECISION,&
        eigen%vec%rc,local_chunk,local_istatemin-1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      else
        eigen%vec%rc=VECP(:,ivp)
      endif
      CALL add_eigen(eigen,lowest) ! APPEND NEW EIGENPAIR TO LIST
      CALL delete_eigen(eigen)
     endif
    ENDDO



    IF(Nconverged>Neigen_) &
    CALL dump_message(TEXT="# WARNING IN Block_Lanczos_diagonalize_real: "//c2s(i2c(Nconverged))//&
    " EIGENPAIRS CONVERGED > "//c2s(i2c(Neigen_))//" REQUESTED")

    IF(Nconverged<Neigen_)THEN
      CALL dump_message(TEXT="# WARNING IN Block_Lanczos_diagonalize_real: "//c2s(i2c(Neigen_))//" EIGENPAIRS REQUESTED "//&
      "BUT ONLY "//c2s(i2c(Nconverged))//" ARE CONVERGED AFTER"//c2s(i2c(Niter))//" ITERATIONS")
      WRITE(log_unit,'(a,f12.10,a)') "# THE FOLLOWING EIGENVALUES (RESIDUALS) ARE DISCARDED [tolerance = ",tolerance,"] :"
      DO ivp=Nconverged+1,Neigen_
        WRITE(log_unit,'(a,E20.12,a,E9.1,a)') "E"//c2s(i2c(ivp))//" = ",VALP(ivp,1)," (",VALP(ivp,2),")"   
      ENDDO
    ENDIF

    CALL delete_eigen(eigen)
    IF(ALLOCATED(vec_in))  DEALLOCATE(vec_in)
 IF(ALLOCATED(vec_out)) DEALLOCATE(vec_out)

    DEALLOCATE(ISTOR,RSTOR,U,V,VALP,VECP)

    CALL timer_fortran(start_diagH,"# DIAGONALIZATION OF "//TRIM(ADJUSTL(title_H_()))//" TOOK "//c2s(i2c(Niter))//" ITERATIONS AND ")

  END SUBROUTINE


!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************


END MODULE 

