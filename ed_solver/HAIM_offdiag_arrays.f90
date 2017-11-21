MODULE HAIM2_class

  USE liste3
  USE AIM2_class

  IMPLICIT NONE

  REAL(DBL), PARAMETER, PRIVATE      :: zero=0.0_DBL
  LOGICAL,   PARAMETER, PRIVATE      :: F=.FALSE.,T=.TRUE.

  !--------------------------------------!  
  ! QUADRATIC AIM HAMILTONIAN IN SECTOR  !
  !--------------------------------------!

  TYPE(AIM2_type),SAVE  :: AIM2up,AIM2do,AIM2sz

  !----------------!
  ! DIAGONAL PART  !
  !----------------!

  REAL(DBL), ALLOCATABLE, TARGET :: diagup(:),diagdo(:),diagsz(:) 

  !--------------------!
  ! OFF-DIAGONAL PART  !
  !--------------------!

  !------------------------------------------------------------------------------------------------------------------------------!
  ! applying Hod on |istate> generates noff(istate) states whose ranks are stored in rankoff(nnnft+1 <= k <= shift+noff(istate)) !
  ! where shift = SUM(noff(1:istate-1)) while the matrix elements are stored in H%offdiag(k) = <rankoff(k)|Hod|istate>           !
  !------------------------------------------------------------------------------------------------------------------------------!

  INTEGER(8), TARGET                :: noffdiagup=0,noffdiagdo=0,noffdiagsz=0
#ifdef _complex
  COMPLEX(DBL), ALLOCATABLE, TARGET ::  offdiagup(:),offdiagdo(:),offdiagsz(:)
#else
  REAL(DBL),    ALLOCATABLE, TARGET ::  offdiagup(:),offdiagdo(:),offdiagsz(:)
#endif
  INTEGER,      ALLOCATABLE, TARGET ::  rankoffup(:),rankoffdo(:),rankoffsz(:)
  INTEGER,      ALLOCATABLE, TARGET ::     noffup(:),   noffdo(:),   noffsz(:)

CONTAINS

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE new_HAIM2(AIM,sector,SPIN) 
    TYPE(AIM_type),            INTENT(IN)         :: AIM
    TYPE(fermion_sector_type), INTENT(IN), TARGET :: sector
    INTEGER, OPTIONAL,         INTENT(IN)         :: SPIN
    INTEGER                                       :: ii 

    IF(PRESENT(SPIN))THEN
      SELECT CASE(SPIN)
       CASE(1)
         !-------------------------!
         ! SECTOR DEPENDANT ARRAYS !
         !-------------------------!
         IF(ALLOCATED(diagup)) DEALLOCATE(diagup);IF(ALLOCATED(noffup)) DEALLOCATE(noffup)
         if(.not.USE_TRANSPOSE_TRICK_MPI)then
          ALLOCATE(noffup(sector%chunk(iproc)),diagup(sector%chunk(iproc)))
         else
          ALLOCATE(noffup(sum(sector%chunk)),diagup(sum(sector%chunk)))
         endif
         !------------------------------!
         ! HAMILTONIAN DEPENDANT ARRAYS !
         !------------------------------!
         CALL new_AIM2(AIM2up,AIM,SPIN=1)
         CALL tab_HAIM2(AIM2up,sector,NO_STORAGE=T)
         IF(ALLOCATED(rankoffup)) DEALLOCATE(rankoffup); IF(ALLOCATED(offdiagup)) DEALLOCATE(offdiagup)
         ALLOCATE(rankoffup(noffdiagup),offdiagup(noffdiagup))
       CASE(2)
         !-------------------------!
         ! SECTOR DEPENDANT ARRAYS !
         !-------------------------!
         IF(ALLOCATED(diagdo)) DEALLOCATE(diagdo); IF(ALLOCATED(noffdo)) DEALLOCATE(noffdo)
         if(.not.USE_TRANSPOSE_TRICK_MPI)then
          ALLOCATE(noffdo(sector%chunk(iproc)),diagdo(sector%chunk(iproc)))
         else
          ALLOCATE(noffdo(sum(sector%chunk)),diagdo(sum(sector%chunk)))
         endif
         !------------------------------!
         ! HAMILTONIAN DEPENDANT ARRAYS !
         !------------------------------!
         CALL new_AIM2(AIM2do,AIM,SPIN=2)
         CALL tab_HAIM2(AIM2do,sector,NO_STORAGE=T)
         IF(ALLOCATED(rankoffdo)) DEALLOCATE(rankoffdo); IF(ALLOCATED(offdiagdo)) DEALLOCATE(offdiagdo)
         ALLOCATE(rankoffdo(noffdiagdo),offdiagdo(noffdiagdo))
      END SELECT

    ELSE 

      ! SECTOR DEPENDANT ARRAYS
      IF(ALLOCATED(noffsz)) DEALLOCATE(noffsz); IF(ALLOCATED(diagsz)) DEALLOCATE(diagsz)
      if(.not.ON_FLY) ALLOCATE(noffsz(sector%chunk(iproc)),diagsz(sector%chunk(iproc)))
 
      ! HAMILTONIAN DEPENDANT ARRAYS
      CALL new_AIM2(AIM2sz,AIM)
      if(.not.ON_FLY) CALL tab_HAIM2(AIM2sz,sector,NO_STORAGE=T)
      IF(ALLOCATED(rankoffsz)) DEALLOCATE(rankoffsz); IF(ALLOCATED(offdiagsz)) DEALLOCATE(offdiagsz)
      if(.not.ON_FLY) ALLOCATE(rankoffsz(noffdiagsz),offdiagsz(noffdiagsz))

    ENDIF


  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE delete_HAIM2(SPIN)
    INTEGER, OPTIONAL, INTENT(IN) :: SPIN
    IF(PRESENT(SPIN))THEN
      SELECT CASE(SPIN)
       CASE(1)
        CALL delete_AIM2(AIM2up)
        IF(ALLOCATED   (diagup)) DEALLOCATE   (diagup)
        IF(ALLOCATED(offdiagup)) DEALLOCATE(offdiagup)
        IF(ALLOCATED(rankoffup)) DEALLOCATE(rankoffup)
        IF(ALLOCATED   (noffup)) DEALLOCATE   (noffup)
       CASE(2)
        CALL delete_AIM2(AIM2do)
        IF(ALLOCATED   (diagdo)) DEALLOCATE   (diagdo)
        IF(ALLOCATED(offdiagdo)) DEALLOCATE(offdiagdo)
        IF(ALLOCATED(rankoffdo)) DEALLOCATE(rankoffdo)
        IF(ALLOCATED   (noffdo)) DEALLOCATE   (noffdo)
      END SELECT
    ELSE
      CALL delete_AIM2(AIM2sz)
      IF(ALLOCATED   (diagsz)) DEALLOCATE   (diagsz)
      IF(ALLOCATED(offdiagsz)) DEALLOCATE(offdiagsz)
      IF(ALLOCATED(rankoffsz)) DEALLOCATE(rankoffsz)
      IF(ALLOCATED   (noffsz)) DEALLOCATE   (noffsz)
    ENDIF
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE tab_HAIM2(AIM2,sector,NO_STORAGE)

    use openmpmod, only: omp_get_num_threads, omp_get_thread_num
    use lockmod,   only: maxt

  !----------------------------------------------------------------------!
  ! TABULATE THE QUADRATIC PART OF THE AIM HAMILTONIAN IN A GIVEN SECTOR ! 
  !----------------------------------------------------------------------!

    TYPE(fermion_sector_type), INTENT(IN) :: sector
    LOGICAL, OPTIONAL,         INTENT(IN) :: NO_STORAGE
    TYPE(AIM2_type),           INTENT(IN) :: AIM2
    LOGICAL                               :: compute_noffdiag_only ! true if we only compute noffdiag
    INTEGER                               :: istate,jstate,istateloc,istatemin,istatemax,iorb,jorb
    INTEGER                               :: IMPnorbs,BATHnorbs,i
    INTEGER,      POINTER                 :: IMPiorb(:) => NULL(),BATHiorb(:) => NULL()
    REAL(DBL),    POINTER                 ::    diag(:) => NULL()
#ifdef _complex
    COMPLEX(DBL), POINTER                 :: offdiag(:) => NULL()
#else
    REAL(DBL),    POINTER                 :: offdiag(:) => NULL()
#endif
    INTEGER,      POINTER                 :: rankoff(:) => NULL(),noff(:) => NULL()
    INTEGER(8),   POINTER                 :: noffdiag => NULL()
    TYPE(fermion_ket_type)                :: ket_in,ket_out 
    TYPE(masked_matrix_type)              :: Ec,Eb,Vbc
    INTEGER                               :: start_tabH
    INTEGER                               :: istatemin0,istatemax0,TID,istatemin_,j
    INTEGER,ALLOCATABLE                   :: imin_(:),imax_(:)
    INTEGER(8)                            :: noffdiagi,check_sum,noffdiag_save
    LOGICAL                               :: go_for_omp


    CALL reset_timer(start_tabH)

    if(verboseall) write(*,*) 'USING OPENMP - RANK : ', OPEN_MP,rank
    if(OPEN_MP)then
    !$OMP PARALLEL 
    if(OMP_GET_NUM_THREADS()/=MAXT)Then
      write(*,*) 'ERROR obtained number of threads is different from MAXT'
      write(*,*) 'MAXT              = ', MAXT
      write(*,*) 'NUMBER OF THREADS = ', OMP_GET_NUM_THREADS()
      stop
    endif
    !$OMP END PARALLEL 
    if(verboseall) write(*,*) 'ALLOCATING IMIN/IMAX WITH SIZE - RANK : ', MAXT, RANK
    allocate(imin_(MAXT),imax_(MAXT))
    endif


    IMPiorb   => AIM2%IMPiorb 
    IMPnorbs  =  AIM2%IMPnorbs
    BATHiorb  => AIM2%BATHiorb 
    BATHnorbs =  AIM2%BATHnorbs

    CALL new_masked_matrix( Ec,  AIM2%Ec  )
    CALL new_masked_matrix( Eb,  AIM2%Eb  )
    CALL new_masked_matrix( Vbc, AIM2%Vbc )

                            compute_noffdiag_only = F
    IF(PRESENT(NO_STORAGE)) compute_noffdiag_only = NO_STORAGE

    SELECT CASE(AIM2%spin)
      CASE(0)
       noffdiag => noffdiagsz
      CASE(1)
       noffdiag => noffdiagup
      CASE(2)
       noffdiag => noffdiagdo
    END SELECT

    SELECT CASE(AIM2%spin)
       CASE(0)
         noff  =>  noffsz
       CASE(1)
         noff  =>  noffup
       CASE(2)
         noff  =>  noffdo
    END SELECT

    IF(.NOT.compute_noffdiag_only)THEN 
      SELECT CASE(AIM2%spin)
       CASE(0)
         diag    =>    diagsz
         offdiag => offdiagsz
         rankoff => rankoffsz
       CASE(1)
         diag    =>    diagup
         offdiag => offdiagup
         rankoff => rankoffup
       CASE(2)
         diag    =>    diagdo
         offdiag => offdiagdo
         rankoff => rankoffdo
      END SELECT
      diag = zero; offdiag = zero; rankoff = 0
    ENDIF
   
    if(size(noff)==0)then
      write(*,*) 'WARNING NOFF IS EMPTY -empty chunck '
      noffdiag=0; check_sum=0; goto 35
    endif
 
    !=================================================================================!
    !=================================================================================!
    !=================================================================================!
    !=================================================================================!

    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!

    if(compute_noffdiag_only) then
     noff = 0
    else
     check_sum=long_sum(noff)
    endif

   if(.not.USE_TRANSPOSE_TRICK_MPI)then
    istatemin = sector%istatemin(iproc)
    istatemax = sector%istatemax(iproc)
   else
    istatemin = minval(sector%istatemin)
    istatemax = maxval(sector%istatemax)
   endif

#ifndef OPENMP_MPI
   if(OPEN_MP) call openmp_split_array(sector%dimen,imin_,imax_)
   if(iproc/=1)then
     write(*,*) 'ERROR : using openmp and MPI at the same time'
     write(*,*) 'either set FLAG_MPI to >0, or recompile with -DOPENMP_MPI'
     stop
   endif
#else
   if(USE_TRANSPOSE_TRICK_MPI)Then
     write(*,*) 'USE_TRANSPOSE_TRICK_MPI not compatible with OPENMP+MPI';stop
   endif
   if(OPEN_MP) then
       if(sector%istatemax(iproc)>0)then
          if(size(imax_)/=MAXT.or.size(imin_)/=MAXT)Then
            write(*,*) 'imin or imax has wrong dim'
            write(*,*) 'imin_ : ', shape(imin_)
            write(*,*) 'imax_ : ', shape(imax_)
            write(*,*) 'MAXT  : ', MAXT
            stop
          endif
          if(verboseall) write(*,*) 'SPLITTING HILBERT SPACE ON OMP, DIM = ', sector%istatemax(iproc)-sector%istatemin(iproc)+1
          call openmp_split_array(sector%istatemax(iproc)-sector%istatemin(iproc)+1,imin_,imax_)
          where(imin_/=0) imin_=imin_+sector%istatemin(iproc)-1
          where(imax_/=0) imax_=imax_+sector%istatemin(iproc)-1
       else
          imin_=0;imax_=0
       endif
   endif
#endif

              go_for_omp=.false.
  if(OPEN_MP) go_for_omp=minval(imax_)>0.and.minval(imin_)>0.and.ALL(sector%istatemax>0)

  !$OMP PARALLEL IF(go_for_omp) PRIVATE(TID,istatemin_,istatemin0,istatemax0,istate,istateloc,jstate,iorb,jorb,noffdiag_save,ket_in,ket_out,noffdiagi)

   if(OPEN_MP.and.go_for_omp)then
    TID        = OMP_GET_THREAD_NUM()+1
    istatemin0 = imin_(TID) ; istatemax0 = imax_(TID)  
    istatemin_ = 1 
#ifdef OPENMP_MPI
    istatemin_ = istatemin
#endif
    if(.NOT.compute_noffdiag_only)then
#ifndef OPENMP_MPI
     if(istatemin0>1)then
#else
     if(istatemin0>istatemin_)then
#endif
#ifdef OPENMP_MPI
        noffdiagi = long_sum(noff(1:istatemin0-istatemin_)) 
#else
        noffdiagi = long_sum(noff(1:istatemin0-1))
#endif
       if(noffdiagi<0) stop 'error in starting openmp, noffdiagi <0'
     else
        noffdiagi = 0
     endif
    else
        noffdiagi = 0
    endif
   else
     istatemin0=istatemin ; istatemax0=istatemax ; istatemin_=istatemin ; noffdiagi = 0  
   endif

   if(istatemax0==0) goto 1992

    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!

    if(noffdiagi<0) stop 'error tabulate, negative noffdiagi'

    DO istate=istatemin0,istatemax0 

      istateloc = istate - istatemin_ + 1

     !CALL check_openmp

      CALL new_ket(ket_in,sector%state(istate),sector%norbs)

      noffdiag_save = noffdiagi 

      IF(.NOT.compute_noffdiag_only)THEN
       ! IMPURITY QUADRATIC
       DO iorb=1,IMPnorbs
          IF(is_occupied(IMPiorb(iorb),ket_in))  diag(istateloc) = diag(istateloc) + Ec%rc%mat(iorb,iorb)
       ENDDO
       ! BATH
       DO iorb=1,BATHnorbs
          IF(is_occupied(BATHiorb(iorb),ket_in)) diag(istateloc) = diag(istateloc) + Eb%rc%mat(iorb,iorb)
       ENDDO
       if(do_quench==1)then
        if(AIM2%spin==0)then !Sz
         DO iorb=1,IMPnorbs
          IF(is_occupied(IMPiorb(iorb),ket_in))  diag(istateloc) = diag(istateloc) + quench_mag
         ENDDO
        endif
        if(AIM2%spin==1)then !up
         DO iorb=1,IMPnorbs
          IF(is_occupied(IMPiorb(iorb),ket_in))  diag(istateloc) = diag(istateloc) + quench_mag
         ENDDO
        endif
        if(AIM2%spin==2)then !do
         DO iorb=1,IMPnorbs
          IF(is_occupied(IMPiorb(iorb),ket_in))  diag(istateloc) = diag(istateloc) - quench_mag
         ENDDO
        endif
       endif
      ENDIF
     
      ! IMPURITY
      DO iorb=1,IMPnorbs; DO jorb=1,IMPnorbs; IF(iorb/=jorb.AND.Ec%rc%MASK%mat(iorb,jorb))THEN

       if(abs(Ec%rc%mat(iorb,jorb))>cutoff_hamilt_param)then
        CALL hop(ket_out,IMPiorb(iorb),IMPiorb(jorb),ket_in)
        IF(.NOT.ket_out%is_nil)THEN
         jstate = sector%rank(ket_out%state)
         IF(USE_CC.or.jstate<istate)THEN
           noffdiagi = noffdiagi + 1
           IF(.NOT.compute_noffdiag_only)THEN
             rankoff(noffdiagi) = jstate
             offdiag(noffdiagi) = Ec%rc%mat(iorb,jorb) * ket_out%fermion_sign  
           ENDIF
         ENDIF
        endif

       ENDIF
      ENDIF; ENDDO; ENDDO



      ! BATH
      DO iorb=1,BATHnorbs; DO jorb=1,BATHnorbs; IF(iorb/=jorb.AND.Eb%rc%MASK%mat(iorb,jorb))THEN
       if(abs(Eb%rc%mat(iorb,jorb))>cutoff_hamilt_param)then

       CALL hop(ket_out,BATHiorb(iorb),BATHiorb(jorb),ket_in)

       IF(.NOT.ket_out%is_nil)THEN
         jstate = sector%rank(ket_out%state)
         IF(USE_CC.or.jstate<istate)THEN
           noffdiagi = noffdiagi + 1
           IF(.NOT.compute_noffdiag_only)THEN
             rankoff(noffdiagi) = jstate
             offdiag(noffdiagi) = Eb%rc%mat(iorb,jorb) * ket_out%fermion_sign  
           ENDIF
         ENDIF
       ENDIF

       endif
      ENDIF; ENDDO; ENDDO


      ! HYBRIDIZATION

      DO iorb=1,BATHnorbs; DO jorb=1,IMPnorbs; IF(Vbc%rc%MASK%mat(iorb,jorb))THEN

       !------------------!
       ! IMPURITY -> BATH !
       !------------------!
      if(abs(Vbc%rc%mat(iorb,jorb))>cutoff_hamilt_param)then
       CALL hop(ket_out,BATHiorb(iorb),IMPiorb(jorb),ket_in)

       IF(.NOT.ket_out%is_nil)THEN
         jstate = sector%rank(ket_out%state)
         IF(USE_CC.or.jstate<istate)THEN
           noffdiagi = noffdiagi + 1
           IF(.NOT.compute_noffdiag_only)THEN
             rankoff(noffdiagi) = jstate
             offdiag(noffdiagi) = Vbc%rc%mat(iorb,jorb) * ket_out%fermion_sign                     
           ENDIF
         ENDIF
       ENDIF

       !------------------!
       ! BATH -> IMPURITY !
       !------------------!
       CALL hop(ket_out,IMPiorb(jorb),BATHiorb(iorb),ket_in)

       IF(.NOT.ket_out%is_nil)THEN
         jstate = sector%rank(ket_out%state)
         IF(USE_CC.or.jstate<istate)THEN
           noffdiagi = noffdiagi + 1
           IF(.NOT.compute_noffdiag_only)THEN
             rankoff(noffdiagi) = jstate
             offdiag(noffdiagi) = conj(Vbc%rc%mat(iorb,jorb)) * ket_out%fermion_sign                     
           ENDIF
         ENDIF
       ENDIF
      endif

      ENDIF; ENDDO; ENDDO

     ! number of off-diagonal elements generated by H |istate>  
      if(istateloc<0.or.noffdiagi<0.or.noffdiag_save<0) stop 'error in tabulate HAIM, negative noffdiag' 

      if(istateloc>size(noff))then
       write(*,*) 'state outside array noff ! ' 
       write(*,*) 'stateloc              : ', istateloc
       write(*,*) 'rank                  : ', rank
       write(*,*) 'size(noff)            : ', size(noff) 
       write(*,*) 'istatemin0,istatemax0 : ', istatemin0,istatemax0
       write(*,*) 'istatemin,istatemax   : ', istatemin,istatemax
       stop
      endif

      noff(istateloc) = noffdiagi - noffdiag_save 

    ENDDO 

    1992 continue

   !$OMP END PARALLEL

    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!
    !------------------------------------------------------------------------!

   noffdiag=long_sum(noff)

   35 continue

   write(log_unit,*) 'IN SECTOR THERE ARE ... OFF DIAG ELEMENTS : ', noffdiag

   if(.not.compute_noffdiag_only .and. check_sum /= noffdiag   ) stop 'error check_sum on noff not valid'
   if( noffdiag<0 ) then
     write(*,*) ' number of OFF_DIAG terms of the hamiltonian matrix are : ', noffdiag
     write(*,*) ' number of diagonal terms                               : ', size(noff) 
     do i=1,size(noff),100000
      noffdiag=long_sum(noff(1:min(size(noff),i)))
      write(*,*) 'partial sum of off-diag terms : ', i,noffdiag,long_sum(noff(1:min(size(noff),i))) 
     enddo
     stop 'error noffdiag <0, critical '
   endif


    !=================================================================================!
    !=================================================================================!
    !=================================================================================!
    !=================================================================================!

    ! CLEAN UP
    CALL delete_masked_matrix(Ec); CALL delete_masked_matrix(Eb); CALL delete_masked_matrix(Vbc); 
    CALL timer_fortran(start_tabH,"# TABULATION OF THE HAMILTONIAN IN SECTOR "//TRIM(sector%title)//" TOOK") 

    if(allocated(imin_)) deallocate(imin_,imax_)


  contains

    !---------------!

   subroutine check_openmp
      if(istate>size(sector%state))then
       write(*,*) 'running open mp ?     : ', OPEN_MP
       write(*,*) 'size of hilbert space : ', size(sector%state),sector%dimen
       write(*,*) 'how many threads      : ', MAXT
       write(*,*) 'min_                  : ', imin_
       write(*,*) 'max_                  : ', imax_
       write(*,*) 'istatemin0,istatemax0 : ', istatemin0,istatemax0
       write(*,*) 'noffdiag start        : ', noffdiagi
       stop
      endif
   end subroutine

    !---------------!

  END SUBROUTINE



!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

END MODULE 
