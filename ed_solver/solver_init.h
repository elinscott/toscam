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


  SUBROUTINE init_solver(FILEIN,AIM,OLDGSFILE) 
    use namelistmod, only: putel_in_namelist
   
    !------------------------! 
    ! READ SOLVER PARAMETERS !
    !------------------------!

    TYPE(AIM_type),   INTENT(IN)           :: AIM
    CHARACTER(LEN=*), INTENT(IN)           :: FILEIN
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: OLDGSFILE
    INTEGER                                :: isec,sz2,npart,nup,ndo,isector,ff,nsec_,k,i,j,fillup
    TYPE(eigensector_type)                 :: es ! for clarity
    TYPE(sector_type)                      :: sec 
    type(namelist_set)                     :: nm
    INTEGER,ALLOCATABLE                    :: list0(:)
    logical                                :: discard_full_state

    discard_full_state=.false.

!===============================================================================================================!

    write(log_unit,*) 'initialize SOLVER, number of sites is : ', AIM%bath%Nb+AIM%Nc
    fillup=AIM%bath%Nb+AIM%Nc

    if(nsec>0)then
      nsec_ = nsec 
      if(AIM%BATH%SUPER)then
          if(nsec0+nsec_-1==fillup) nsec_=nsec_-1
      endif
        if(.not.AIM%BATH%SUPER.and..not.SCAN_FULL_NUP_NDN)then
          if(nsec0+nsec_-1==fillup) nsec_=nsec_-1
      endif
    else  
      if(AIM%BATH%SUPER)then
         nsec_  =  AIM%bath%Nb+AIM%Nc-1
         nsec_  =  nsec_-nsec0+1
      else 
         if(.not.SCAN_FULL_NUP_NDN)then
          nsec_  =  AIM%bath%Nb+AIM%Nc-1
         else
          nsec_  =  AIM%bath%Nb+AIM%Nc
         endif
         nsec_  =  nsec_-nsec0+1
      endif
    endif

    call delete_eigensectorlist(sector2diagH)
    call namelist_init(nm,400,name_of_namelist='ED_SOLVER partial input related to solver')
 
    IF(AIM%bath%SUPER)THEN ! READ LIST OF Sz-SECTORS
      write(log_unit,*) '------> GO FOR SUPERCONDUCTING PHASE : '
      !##############################################################################################!
      ALLOCATE(list_sectors(nsec_),list0(nsec_))
      list0=(/(isec,isec=nsec0,nsec0+nsec_-1)/)
      call putel_in_namelist(nm,list_sectors,'list_sectors',list0,'list of Sz sectors')
      call look_for_namelist_in_file(nm,FILEIN)
      call look_for_command_line_argument(nm)
      write(log_unit,*) ' LIST OF SECTOR TO CONSIDER ARE : ', list0
      !##############################################################################################!
      CALL delete_eigensector(es)
      DO isec=1,nsec_
      ! NUMBER OF NAMBU FERMIONS = Ns + Sz*2
        npart = AIM%Ns + list_sectors(isec) 
        write(log_unit,*) '------------------------------------------------------------'
        write(log_unit,*) 'build up new Sz sector, isec, nsec : ' , isec,nsec_
        write(log_unit,*) 'list of sectors                    : ' , list_sectors(isec)
        write(log_unit,*) 'Npart                              : ' , npart
        write(log_unit,*) '------------------------------------------------------------'
        CALL delete_sector(sec)
 ALLOCATE(sec%sz)
        CALL new_fermion_sector(sec%sz,npart,AIM%norbs,SZ=T) 
        CALL new_eigensector(es,sec) 
        CALL add_eigensector(es,sector2diagH) 
      ENDDO
                 !--------------------------------!
    ELSE         ! READ LIST OF (nup,ndo)-SECTORS !
                 !--------------------------------!

      write(log_unit,*) '------------------------------'
      write(log_unit,*) '---- GO FOR NORMAL PHASE -----'
      write(log_unit,*) '------------------------------'

      if(nsec0<0) then
       write(*,*) ' nup/ndn sectors and you want a negative number of particles?'
       write(*,*) ' sec0 is      : ', nsec0
       write(*,*) ' nsec to scan : ', nsec_
       stop 'critical'
      endif

      if(.not.SCAN_FULL_NUP_NDN)then
       !##############################################################################################!
       ALLOCATE(list_sectors(nsec_*2),list0(nsec_*2))
       list0 = (/([isec,isec],isec=nsec0,nsec0+nsec_-1)/)
       call putel_in_namelist(nm,list_sectors,'list_sectors',list0,'list of sectors nup,ndo')
       call look_for_namelist_in_file(nm,FILEIN)
       call look_for_command_line_argument(nm)
       write(log_unit,*) ' LIST OF SECTOR TO CONSIDER ARE : ', list0
       if(size2>1.and..not.no_mpi) then
          write(*,*) 'FINISHED TO READ LIST_SECTOR, RANK = ', rank
          call mpibarrier
       endif
       !##############################################################################################!
      else
       !##############################################################################################!
       k=2*(nsec_**2)
 
       if(nsec0+nsec_-1==fillup.and.discard_full_state) k=k-2
       ALLOCATE(list_sectors(k),list0(k))
       k=0
       do i=nsec0,nsec0+nsec_-1
        do j=nsec0,nsec0+nsec_-1
         if((i/=fillup.or.j/=fillup).or..not.discard_full_state)then
           k=k+1
 list0(k)=i
 k=k+1
 list0(k)=j
         endif
        enddo
       enddo
       if(k/=size(list0)) then 
         write(*,*) 'k                    : ', k
         write(*,*) 'size list of sectors : ', size(list0)
         stop 'error in solver.f90 in ed solver, size of list0 not ok, critical'
       endif
       call putel_in_namelist(nm,list_sectors,'list_sectors',list0, 'list of sectors nup,ndo')
       call look_for_namelist_in_file(nm,FILEIN)
       call look_for_command_line_argument(nm)
       write(log_unit,*) ' LIST OF SECTOR TO CONSIDER ARE : ', list0
       if(size2>1.and..not.no_mpi) then
          write(*,*) 'FINISHED TO READ LIST_SECTOR, RANK = ', rank
          call mpibarrier
       endif
      !##############################################################################################!
     endif

      CALL delete_eigensector(es)
      DO isec=1,size(list_sectors)/2
        nup = list_sectors(2*isec-1)
        ndo = list_sectors(2*isec) 
        CALL delete_sector(sec)
        IF(associated(sec%updo)) deallocate(sec%updo)
        ALLOCATE(sec%updo)
        CALL new_fermion_sector2(sec%updo,nup,ndo,AIM%Ns,NAMBU=F)
        CALL new_eigensector(es,sec)
        CALL add_eigensector(es,sector2diagH) 
      ENDDO
    ENDIF
   
!===============================================================================================================!
 
    ! COMPLETE LIST OF SECTORS WITH SECTOR OF INPUT GS ... IF NECESSARY
    
    start_from_old_gs = F
    IF(PRESENT(OLDGSFILE))THEN
      start_from_old_gs = T
      !----------------------! 
      ! READ GS IN OLDGSFILE !
      !----------------------!
      CALL read_raw_eigensectorlist(oldGS,OLDGSFILE) 
      ! AT THIS POINT WE HAVE THE LIST OF EIGENSECTORS WITH DUMMY EMPTY SECTOR
      ! THAT WE TRY TO RECONCILE NOW
      DO isector=1,oldGS%nsector
        CALL copy_eigensector(es,oldGS%es(isector))
        ! TEST CONSISTENCY
        IF(ASSOCIATED(es%sector%updo))THEN
          IF(es%sector%updo%up%norbs/=AIM%Ns.OR.es%sector%updo%down%norbs/=AIM%Ns) STOP "ERROR IN init_solver: INCONSISTENT INPUT SECTOR!"
        ELSE IF(ASSOCIATED(es%sector%sz))THEN
          IF(es%sector%sz%norbs/=AIM%norbs) STOP "ERROR IN init_solver: INCONSISTENT INPUT SECTOR!"
        ENDIF
        ! THE INPUT GS ISNT IN sector2diagH: WE NEED TO CREATE & APPEND THIS SECTOR TO sector2diagH!
        IF(rank_sector_in_list(es%sector,sector2diagH)==0)&
        CALL add_eigensector(es,sector2diagH) 
        ! NOW WE ARE SURE THE RELEVANT SECTOR IS INDEED IN sector2diagH
      ENDDO
    ENDIF

    if(allocated(list_sectors)) DEALLOCATE(list_sectors,STAT=istati)
    if(allocated(list0))        DEALLOCATE(list0,STAT=istati)
    CALL delete_eigensector(es) 
 CALL delete_sector(sec)

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
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

