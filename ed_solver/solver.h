  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!

  subroutine get_GS
  implicit none
  TYPE(eigenlist_type), POINTER :: lowest
  INTEGER                       :: rank_,size2_

     if(FLAG_MPI_GREENS>0)then
      rank_=rank
 size2_=size2
     else
      rank_=0
 size2_=1
     endif

     !======================================================================!
     DO isector=rank_+1,GS%nsector,size2_
      if(GS%es(isector)%lowest%neigen>0)then
        call get_info_sector

        if(iterdmft==1 .or. restarted .or. .not. track_sectors .or. &
    & (uup/=-1000.and.(uup<=nupmax+iwindow.and.ddn<=ndnmax+iwindow.and.uup>=nupmin-iwindow.and.ddn>=ndnmin-iwindow) ) .or. &
    & (ssz/=-1000.and.(ssz<= szmax+iwindow.and.ssz>= szmin-iwindow)     )  )then

        if(.not.(USE_TRANSPOSE_TRICK_MPI.and.NOT_COMMENSURATE))then

     !----------------------------------------------------!
        IF(Block_size==0 .or. uup==0 .or. ddn==0 .or. uup==1 .or.  ddn==1  .or.  (uup==itot .and. ddn==itot) )THEN

         lowest => GS%es(isector)%lowest

         call dump_info(6)
         call dump_info(log_unit)
         CALL new_H(AIM,GS%es(isector)%sector)

         SELECT CASE (which_lanczos)
          CASE('NORMAL')
            write(*,*) '---start CPU Lanczos---'
            call Lanczos_get_GS_sector(lowest)
          ! ebl: removing GPU functionality
          ! CASE('GPU')
          !   write(*,*) '---start GPU Lanczos---'
          !   call Lanczos_get_GS_sector_GPU(lowest)
          CASE DEFAULT
            stop 'error solver not defined in get GS'
          END SELECT

         CALL delete_H()

        ELSE

         if(USE_TRANSPOSE_TRICK_MPI) stop 'error block lanczos and lanczos vector split among nodes'

#ifdef _complex
#else
          IF(Block_size<0)THEN
           write(log_unit,*) 'Cullum library not installed - to be added in the future'
           stop
          ELSE
           write(log_unit,*) 'Eigenvectors already computed by Block Lanczos during the first Lanczos pass' 
          ENDIF
#endif
        ENDIF
     !----------------------------------------------------!
        endif
        endif
       endif
      ENDDO
     !======================================================================!


     if(FLAG_MPI_GREENS>0) call sync_eigen_vectors

  end subroutine

  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!

  subroutine get_eigen_energies
  implicit none
  REAL(8) :: dE
  INTEGER :: rank_,size2_
  logical :: goforit,sector_bound_file_present

     INQUIRE(file='ed.sector_bound_file',EXIST=sector_bound_file_present)

     if(sector_bound_file_present)then
        write(*,*) 'WARNING SECTOR BOUND FILE PRESENT - RANK', rank
        open(unit=13910,file='ed.sector_bound_file')
        read(13910,*)        nupmin,nupmax
        read(13910,*)        ndnmin,ndnmax  
        read(13910,*,end=31) szmin,szmax
        31 continue
        close(13910)
        write(*,*) 'SECTOR BOUND FILE PRESENT - nup min max', nupmin,nupmax
        write(*,*) 'SECTOR BOUND FILE PRESENT - ndn min max', ndnmin,ndnmax
     endif

     if(FLAG_MPI_GREENS>0)then
       rank_=rank
 size2_=size2
     else
       rank_=0
 size2_=1
     endif

     !======================================================================!
      DO isector=rank_+1,GS%nsector,size2_

        call get_info_sector

        if(GS%es(isector)%lowest%neigen>0) GS%es(isector)%lowest%eigen(:)%val=huge_

        goforit = iterdmft==1 .or. restarted .or. .not. track_sectors 
        if(sector_bound_file_present) goforit = .false.

        goforit = goforit .or. (uup/=-1000.and.(uup<=nupmax+iwindow.and.ddn<=ndnmax+iwindow.and.uup>=nupmin-iwindow.and.ddn>=ndnmin-iwindow))
        goforit = goforit .or. (ssz/=-1000.and.(ssz<=szmax +iwindow.and.ssz>=szmin -iwindow) ) 


        if(goforit)then

        if(.not.(USE_TRANSPOSE_TRICK_MPI.and.NOT_COMMENSURATE))then 

        call dump_info(6)

        call dump_info(log_unit)

        CALL new_H(AIM,GS%es(isector)%sector)

     !----------------------------------------------------!
        IF(Block_size==0 .or. uup==0 .or. ddn==0 .or. uup==1 .or. ddn==1  .or.  (uup==itot .and. ddn==itot) )THEN

         SELECT CASE (which_lanczos)
          CASE('NORMAL')
            write(*,*) '---start CPU Lanczos---'
            CALL Lanczos_fast_diagonalize(GS%es(isector)%lowest)
          ! ebl: removing GPU functionality
          ! CASE('GPU')
          !   write(*,*) '---start GPU Lanczos---'
          !   CALL Lanczos_fast_diagonalize_GPU(GS%es(isector)%lowest) 
          CASE('FULL_ED')
            write(*,*) '---start CPU/GPU ED---'
            CALL ED_diago(GS%es(isector)%lowest)
          CASE('ARPACK')
            CALL ARPACK_diago(GS%es(isector)%lowest)
          CASE DEFAULT
            stop 'error solver not defined in get eigen'
          END SELECT

        ELSE

         if(USE_TRANSPOSE_TRICK_MPI) stop 'error block lanczos and lanczos vector split among nodes, stop - critical '

#ifdef _complex
            write(log_unit,*) 'Cullum library not installed - to be added in the future'
            stop 
#else
          IF(Block_size<0)THEN
           write(log_unit,*) 'Cullum library not installed - to be added in the future'
           stop
          ELSE
            CALL Block_Lanczos_diagonalize(GS%es(isector)%lowest)
          ENDIF
#endif
        ENDIF
     !----------------------------------------------------!

       if(GS%es(isector)%lowest%neigen>0) then
        if(.not.FLAG_FULL_ED_GREEN)then
           write(log_unit,*) '--------------------------------'
           write(log_unit,*) 'FILTERING EIGENVALUES ON THE FLY'
           dE=minval(GS%es(isector)%lowest%eigen(:)%val)
           write(log_unit,*) 'NupNdn File : ', sector_bound_file_present
           write(log_unit,*) 'nupmin-max  : ', nupmin,nupmax
           write(log_unit,*) 'ndnmin-max  : ', ndnmin,ndnmax
           write(log_unit,*) 'iwindow     : ', iwindow
           write(log_unit,*) 'nup ndn     : ', uup,ddn
           write(log_unit,*) 'WINDOW IS   : ', dE, dE+dEmax
           write(log_unit,*) 'SECTOR IS   : ', isector
           write(log_unit,*) 'N eigen     : ', GS%es(isector)%lowest%neigen 
           write(log_unit,*) 'Energies    : ', GS%es(isector)%lowest%eigen(:)%val
           CALL filter_eigen(GS%es(isector)%lowest,[dE,dE+dEmax])
           write(log_unit,*) '--------------------------------'
        endif
       endif

       CALL delete_H()

       endif
       endif

      ENDDO
     !======================================================================!

     if(FLAG_MPI_GREENS>0) call sync_eigen_energies

     if(present(COMPUTE_ALL_CORREL))then
     if(COMPUTE_ALL_CORREL)then 
     if(FLAG_DUMP_INFO_FOR_GAMMA_VERTEX) then 
       call four_leg_vertex_matrices_routine(AIM,GS)
     endif
     endif
     endif

  end subroutine

  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!

  subroutine sync_eigen_vectors
  implicit none
    integer              :: neigen,ii,iii,jj,kk
    integer              :: dim_space
#ifdef _complex
  COMPLEX(8),ALLOCATABLE :: VECP(:,:)
#else
    REAL(8), ALLOCATABLE :: VECP(:,:)
#endif

   call mpibarrier
   DO iii=1,GS%nsector
    ii=mod(iii-1,size2)
    if(ii==rank)then
 neigen=GS%es(iii)%lowest%neigen
 endif
    call mpibcast(neigen,iii=ii)
      !--------------------------------------------------------!
       if(neigen>0)then
        if(ii==rank)then
 dim_space=GS%es(iii)%lowest%eigen(1)%dim_space
 endif
        call mpibcast(dim_space,iii=ii)
        if(allocated(VECP)) deallocate(VECP)
 allocate(VECP(dim_space,neigen))
        if(ii==rank)then
 do kk=1,neigen
 VECP(:,kk)=GS%es(iii)%lowest%eigen(kk)%vec%rc
 enddo
 endif
        call mpibcast(VECP,iii=ii)
        if(rank/=ii)then
         do jj=1,neigen
           CALL new_rcvector(GS%es(iii)%lowest%eigen(jj)%vec,dim_space)
           GS%es(iii)%lowest%eigen(jj)%vec%rc=VECP(:,jj)
         enddo
        endif
       endif
      !--------------------------------------------------------!
   enddo
   call mpibarrier
 
  end subroutine

  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!

  subroutine sync_eigen_energies
  implicit none
    integer              :: neigen,ii,iii,jj,dim_space,kk
    integer, allocatable :: dim_lanc(:)
    REAL(8), ALLOCATABLE :: coef_lanc(:,:)
    REAL(8), ALLOCATABLE :: VAL(:)
    TYPE(eigen_type)     :: eigen
    TYPE(rcvector_type)  :: tempvec

   call mpibarrier 
   DO iii=1,GS%nsector
    ii=mod(iii-1,size2)
    if(ii==rank)then
 neigen=GS%es(iii)%lowest%neigen
 endif 
    call mpibcast(neigen,iii=ii)
      !--------------------------------------------------------!
       if(neigen>0)then
        if(allocated(dim_lanc)) deallocate(dim_lanc)
 allocate(dim_lanc(neigen))
        if(ii==rank)then
 dim_lanc=GS%es(iii)%lowest%eigen(1:neigen)%lanczos_iter
 endif
        call mpibcast(dim_lanc,iii=ii)
        if(ii==rank)then
 dim_space=GS%es(iii)%lowest%eigen(1)%dim_sector
 endif
        call mpibcast(dim_space,iii=ii)
        if(allocated(VAL)) deallocate(VAL)
 allocate(VAL(neigen))
        if(ii==rank)then
 VAL=GS%es(iii)%lowest%eigen(1:neigen)%val
 endif 
        call mpibcast(VAL,iii=ii)

        if(allocated(coef_lanc)) deallocate(coef_lanc)
 allocate(coef_lanc(maxval(dim_lanc(:)),neigen))
        if(ii==rank)then
 do kk=1,neigen
 coef_lanc(1:dim_lanc(kk),kk)=GS%es(iii)%lowest%eigen(kk)%lanczos_vecp(1:dim_lanc(kk)) 
 enddo
 endif
        call mpibcast(coef_lanc,iii=ii)

        if(ii/=rank)then
         do jj=1,neigen
          CALL new_eigen(eigen,VAL(jj),tempvec,.true.,RANK=jj,no_vector=.true.)
          eigen%lanczos_iter=dim_lanc(jj)
          eigen%lanczos_vecp(1:eigen%lanczos_iter)=coef_lanc(1:dim_lanc(jj),jj)
          eigen%rdist=0.
          eigen%dim_space=dim_space
          eigen%dim_sector=dim_space
          CALL add_eigen(eigen,GS%es(iii)%lowest)
 CALL delete_eigen(eigen)
          write(*,*) '  sync eigenvalues from Lanzcos, FLAG_MPI_GREENS '
          write(*,*) '  N eigenvalues    : ', GS%es(iii)%lowest%neigen
          write(*,*) '  eigenvalues      : ', GS%es(iii)%lowest%eigen(jj)%val
          write(*,*) '  SECTOR dimension : ', GS%es(iii)%lowest%eigen(jj)%dim_space
         enddo
        endif

       endif
      !--------------------------------------------------------!
   enddo
   call mpibarrier
 
  return
  end subroutine

  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!

  subroutine get_new_bounds
  implicit none
  logical :: sector_bound_file_present

    szmin= 1000
nupmin= 1000
ndnmin= 1000
    szmax=-1000
nupmax=-1000
ndnmax=-1000
    DO isector=1,GS%nsector
      if(GS%es(isector)%lowest%neigen>0)then
         if(AIM%bath%SUPER)then
           szmin=min(szmin,GS%es(isector)%sector%sz%npart)
           szmax=max(szmax,GS%es(isector)%sector%sz%npart)
         else
           nupmin=min(nupmin,GS%es(isector)%sector%updo%up%npart)
           nupmax=max(nupmax,GS%es(isector)%sector%updo%up%npart)
           ndnmin=min(ndnmin,GS%es(isector)%sector%updo%down%npart)
           ndnmax=max(ndnmax,GS%es(isector)%sector%updo%down%npart)
         endif
      endif 
    END DO

     INQUIRE(file='ed.sector_bound_file',EXIST=sector_bound_file_present)
     if(rank==0.and.sector_bound_file_present)then
#ifndef NO_SYS_CALL
        call system("rm ed.sector_bound_file")
#else
        call remove_filef("ed.sector_bound_file")
#endif
        open(unit=13910,file='ed.sector_bound_file')
        write(13910,*) nupmin,nupmax
        write(13910,*) ndnmin,ndnmax
        write(13910,*)  szmin, szmax
        close(13910)
     endif

    write(log_unit,*) '========================================================='
    write(log_unit,*) ' NEW BOUNDS (Min,Max) sectors to scan for next iteration '
    write(log_unit,*) '  nup bounds : ', nupmin,nupmax
    write(log_unit,*) '  ndn bounds : ', ndnmin,ndnmax
    write(log_unit,*) '  Sz  bounds : ', szmin, szmax 
    write(log_unit,*) '========================================================='

  end subroutine

  !--------------------!
  !--------------------!

#ifdef NO_SYS_CALL
    subroutine remove_filef(filename1)
    implicit none
    character*(*)  :: filename1
    integer        :: unit_,ios
    logical        :: is_it_opened,checkfile
      write(*,*) 'ERASING FILE : ', trim(adjustl(filename1))
      inquire(file=trim(adjustl(filename1)),exist=checkfile)
      if(.not.checkfile) return
      unit_=20
      do
       unit_=unit_+1
       INQUIRE(unit=unit_,OPENED=is_it_opened,iostat=ios)
       if(.not.is_it_opened.and.ios==0)exit
      enddo
      open(unit=unit_,file=trim(adjustl(filename1)))
      close(unit_,status='delete')
    end subroutine
#endif

  !--------------------!
  !--------------------!
  !--------------------!

  subroutine get_info_sector
       if(AIM%bath%SUPER)then
         uup   =  -1000
         ddn   =  -1000
         ssz   =   GS%es(isector)%sector%sz%npart
         if(USE_TRANSPOSE_TRICK_MPI) &
            & stop 'error sz basis and lanczos vector split among nodes, not implemented - critical'
       else
         ssz   =  -1000
         nstd  =   GS%es(isector)%sector%updo%up%nstates
         nstu  =   GS%es(isector)%sector%updo%down%nstates
         uup   =   GS%es(isector)%sector%updo%up%npart
         ddn   =   GS%es(isector)%sector%updo%down%npart
         NOT_COMMENSURATE=not_commensurate_sector(GS,isector)
       endif
  end subroutine

  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!

  subroutine dump_info(unit)
  integer :: unit
       if(AIM%bath%SUPER)then
         write(unit,*) '====================================================='
         write(unit,*) 'diagonalize sector Sz, # sites    : ', ssz,isector,GS%nsector,itot
         write(unit,*) 'label of sector                   : ', isector
         write(unit,*) 'number of mpi nodes               : ', size2
         write(unit,*) '====================================================='
       else
         write(unit,*) '====================================================='
         write(unit,*) 'diagonalize sector up/dn, # sites : ', uup, ddn, itot
         write(unit,*) 'number of up and dn states        : ', nstu,nstd
         write(unit,*) 'label of sector                   : ', isector
         write(unit,*) 'number of mpi nodes               : ', size2
         write(unit,*) 'chunk up spin                     : ', GS%es(isector)%sector%updo%up%chunk
         write(unit,*) 'chunk dn spin                     : ', GS%es(isector)%sector%updo%down%chunk
         write(unit,*) '====================================================='
       endif
  end subroutine 

  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!


  subroutine dumpmess1
      CALL dump_message(TEXT="#")
      CALL dump_message(TEXT="###############################")
      CALL dump_message(TEXT="### READ GS OF THE ANDERSON ###")
      CALL dump_message(TEXT="###   IMPURITY      MODEL   ###")
      CALL dump_message(TEXT="###   FROM  PREVIOUS  RUN   ###")
      CALL dump_message(TEXT="###############################")
      CALL dump_message(TEXT="#")
  end subroutine

  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!

  subroutine dumpmess2
      CALL dump_message(TEXT="#")
      CALL dump_message(TEXT="#############################")
      CALL dump_message(TEXT="### BEGIN DIAGONALIZATION ###")
      CALL dump_message(TEXT="### OF    THE    ANDERSON ###")
      CALL dump_message(TEXT="###   IMPURITY    MODEL   ###")
      CALL dump_message(TEXT="#############################")
      CALL dump_message(TEXT="#")
  end subroutine

  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
  !--------------------!
