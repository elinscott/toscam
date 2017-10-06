MODULE HAIMupdo_class

  use HAIMsz_class

  IMPLICIT NONE

  !-------------------------------------------------!
  ! FULL HAMILTONIAN IN FERMION SECTOR (nup,ndo)    !
  !-------------------------------------------------!

  REAL(DBL), PARAMETER      , PRIVATE :: zero=0.0_DBL
  LOGICAL,  PARAMETER       , PRIVATE :: F=.FALSE.,T=.TRUE.
  INTEGER                   , PRIVATE :: iupmin,iupmax,idomin,idomax,strideup,stridedo,dimen
  INTEGER                   , PRIVATE :: Nc,Nb
  INTEGER,   ALLOCATABLE    , PRIVATE :: IMPiorbup(:),IMPiorbdo(:)
  REAL(DBL), ALLOCATABLE    , PRIVATE :: U(:,:)
  LOGICAL,   ALLOCATABLE    , PRIVATE :: UMASK(:,:)
  LOGICAL                   , PRIVATE :: offdiag_coulomb

  TYPE(masked_real_matrix_type),SAVE  :: QUART_INT_UPDO
  TYPE(fermion_sector2_type), POINTER :: sector => NULL()

  INTERFACE HAIMupdo_mult
    MODULE PROCEDURE HAIMupdo_multc,HAIMupdo_multr
  END INTERFACE 

  INTERFACE HAIMupdo_mult_split
    MODULE PROCEDURE HAIMupdo_multr_split,HAIMupdo_multc_split
  END INTERFACE


CONTAINS

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

#include "HAIMupdo_gpu.h"

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE new_HAIMupdo(AIM,sector_in) 
    TYPE(AIM_type),             INTENT(IN)         :: AIM
    TYPE(fermion_sector2_type), INTENT(IN), TARGET :: sector_in
    LOGICAL                                        :: is_diag(AIM%Nc,AIM%Nc)
    INTEGER                                        :: spin,site,ii

    sector => sector_in
 
    !-----------------------------------! 
    ! QUADRATIC PART OF THE HAMILTONIAN !
    !-----------------------------------!

    CALL new_HAIM2(AIM,sector%up,  SPIN=1)
    CALL new_HAIM2(AIM,sector%down,SPIN=2)

    !--------------! 
    ! QUARTIC PART !
    !--------------!
 
    CALL new_masked_real_matrix(QUART_INT_UPDO,AIM%impurity%U)
    QUART_INT_UPDO%MASK%mat = F
    WHERE(QUART_INT_UPDO%mat/=zero)
      QUART_INT_UPDO%MASK%mat = T
    END WHERE

    !------------------------------------------------! 
    ! TABULATE THE QUADRATIC PART OF THE HAMILTONIAN !
    !------------------------------------------------!

    CALL tab_HAIM2(AIM2up,sector%up)
    CALL tab_HAIM2(AIM2do,sector%down)

    dimen    = sector%dimen
    iupmin   = sector%up%istatemin(iproc); idomin   = sector%down%istatemin(iproc)
    iupmax   = sector%up%istatemax(iproc); idomax   = sector%down%istatemax(iproc)
    strideup = sector%strideup;            stridedo = sector%stridedown  
    
    if(messages)then 
      write(log_unit,*) ' ============================================================ ' 
      write(log_unit,*) ' on cpu ', iproc, ' chunk of states is : ', iupmin,iupmax,idomin,idomax 
      write(log_unit,*) '                          with stride  : ', strideup
      write(log_unit,*) ' ============================================================ '
    endif

    if(iupmax<0 .or. iupmin<0)then
      write(*,*) 'ERROR oups something is wrong ! iupmax negative'
      write(*,*) '      iupmax/iupmin for my rank : ', iupmax,iupmin,rank
      stop
    endif

    Nc  =    AIM%Nc
    Nb  = 2**AIM%Nb

    IF(ALLOCATED(IMPiorbup)) DEALLOCATE(IMPiorbup)
    IF(ALLOCATED(IMPiorbdo)) DEALLOCATE(IMPiorbdo)
    ALLOCATE(IMPiorbup(Nc),IMPiorbdo(Nc))
    IMPiorbup = AIM2up%IMPiorb
    IMPiorbdo = AIM2do%IMPiorb

    IF(ALLOCATED(UMASK)) DEALLOCATE(UMASK)
    ALLOCATE(UMASK(Nc,Nc))
    UMASK     = QUART_INT_UPDO%MASK%mat
    IF(ALLOCATED(U)) DEALLOCATE(U)
    ALLOCATE(U(Nc,Nc))
    U         = QUART_INT_UPDO%mat
    CALL new_diag(is_diag,Nc)
    offdiag_coulomb = ANY(UMASK.AND.(.NOT.is_diag))

  END SUBROUTINE 

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE delete_HAIMupdo
    NULLIFY(sector)
    CALL delete_HAIM2(SPIN=1)
    CALL delete_HAIM2(SPIN=2)
    CALL delete_masked_real_matrix(QUART_INT_UPDO)
    IF(ALLOCATED(IMPiorbup)) DEALLOCATE(IMPiorbup)
    IF(ALLOCATED(IMPiorbdo)) DEALLOCATE(IMPiorbdo)
    IF(ALLOCATED(U))         DEALLOCATE(U)
    IF(ALLOCATED(UMASK))     DEALLOCATE(UMASK)
  END SUBROUTINE 

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE HAIMupdo_multc(vec_out,vec_in)
  implicit none

  !-------------------------------------------------------!
  ! WE COMPUTE THE RELEVANT CHUNK OF vec_out = H * vec_in !
  !-------------------------------------------------------!

    COMPLEX(DBL),ALLOCATABLE    :: vec_tot_out(:)
    COMPLEX(DBL)                :: vec_out(:)
    COMPLEX(DBL)                :: vec_in(:)
    COMPLEX(DBL)                :: hoffdiagup,choffdiagup,hoffdiagdo,choffdiagdo
    INTEGER                     :: site1,site2,n1,n2,j
    INTEGER                     :: istate,   istatemin,   istatemax
    INTEGER                     :: jstate,   jstatemin,   jstatemax
    INTEGER                     :: istateloc,istateminloc,istatemaxloc
    INTEGER                     :: iup,ido,jup,jdo,iuploc,idoloc,irank
    INTEGER                     :: noff,I,sec,TID,iupmin0,min_,max_
    INTEGER(4),allocatable      :: imin_(:),imax_(:)
    INTEGER(8)                  :: nhoffdiagup,nhoffdiagdo
    LOGICAL                     :: bn(2),bn1(2),bn2(2),go_for_omp
    INTEGER                     :: norbs_,jj
    TYPE(fermion_ket_type)      :: up_in,do_in,up_out,do_out
    INTEGER                     :: iupb,idob,istate_,jstate_
    INTEGER                     :: up_in_,do_in_,i1,i2,do_out_,up_out_,fermion_sign_
    integer                     :: fermion_signb_,site3,site4,n3,n4,is1,is2,up_outb_,do_outb_
    logical                     :: bn3(2),bn4(2)
    integer                     :: up_outbb_, do_outbb_, up_outbbb_, do_outbbb_, fermion_signbb_, fermion_signbbb_
    integer                     :: IMPstate,BATHup,BATHdo,iupimp,idoimp,iupimp_out,idoimp_out,AIMup,AIMdo
    logical                     :: emptychunks

    if(dimen==1)then
     vec_out(1)=1
     return
    endif

    if(verboseall) write(*,*) 'USING OPENMP - RANK : ', OPEN_MP,rank

    if(OPEN_MP)then
    if(verboseall) write(*,*) 'ALLOCATING IMIN/IMAX WITH SIZE - RANK : ', MAXT, RANK
    if(verboseall) write(*,*) 'IUPMIN IUPMAX                         : ', iupmin,iupmax
    allocate(imin_(MAXT),imax_(MAXT))
    endif

    sec = sector%up%dimen
  
    if(sec==0) write(*,*) 'danger in hmult up do dim up is 0'
    if(iupmax<iupmin) write(*,*) 'WARNING  iupmax<iupmin'

    emptychunks=any(sector%chunk==0)
 
    if(HILBERT_SPACE_SPLITED_AMONG_NODES.and.USE_CC) then
      if(.not.emptychunks) then
         vec_out( 1+sum(sector%chunk(1:iproc-1)) : sum(sector%chunk(1:iproc)) ) = 0.d0
      else
         vec_out=0.d0
      endif
     else
      vec_out=0.d0
    endif

    if(sector%chunk(iproc)==0) then
     write(*,*) 'EMPTY SECTOR - SIZE VEC_OUT : ', size(vec_out)
     goto 35
    endif

    min_=iupmin; max_=iupmax

#ifndef OPENMP_MPI
   if(OPEN_MP) then
     call openmp_split_array(sec,imin_,imax_)
   endif
#else
   if(OPEN_MP) then
     call openmp_split_array(iupmax-iupmin+1,imin_,imax_)
     where(imin_/=0) imin_=imin_+iupmin-1
     where(imax_/=0) imax_=imax_+iupmin-1
     if(verboseall)then
     do i=1,MAXT
     write(*,*) 'CHUNK : ', imin_(i),imax_(i)
     enddo
     endif
   endif
#endif

                go_for_omp=.false.
    if(OPEN_MP) go_for_omp=minval(imax_)>0.and.minval(sector%chunk)>0..and.minval(imin_)>0

    !$OMP PARALLEL IF(go_for_omp) PRIVATE( IMPstate,BATHup,BATHdo,iupimp,idoimp,iupimp_out,idoimp_out,AIMup,AIMdo,up_outbb_, do_outbb_, up_outbbb_, do_outbbb_, fermion_signbb_, fermion_signbbb_,fermion_signb_,site3,site4,n3,n4,is1,is2,up_outb_,do_outb_,bn3,bn4,n1,n2,iupb,idob,iupmin0,iup,iupmin,iupmax,ido,istate,TID,site1,bn,bn1,bn2,site2,iuploc,idoloc,istatemin,istatemax, nhoffdiagup,nhoffdiagdo,noff,irank,hoffdiagup,choffdiagup,hoffdiagdo,choffdiagdo,jdo,jup,jstate,up_in_,do_in_,i1,i2,jj,do_out_,up_out_,fermion_sign_,istate_,jstate_) SHARED(vec_out,vec_in,sec,diagup,diagdo)

    if(OPEN_MP.and.go_for_omp) then
      TID=OMP_GET_THREAD_NUM()+1; iupmin=imin_(TID); iupmax=imax_(TID); 
    else
      iupmin=min_;iupmax=max_
    endif

    if(OPEN_MP.and.go_for_omp)then
    if(OMP_GET_NUM_THREADS()/=MAXT)Then
     write(*,*) 'ERROR obtained number of threads is different from MAXT'
     write(*,*) 'MAXT = ', MAXT
     write(*,*) 'NUMBER OF THREADS = ', OMP_GET_NUM_THREADS()
     stop
    endif
    endif

    if(iupmax==0) goto 1193

    if(do_quench/=2)then
     DO iup=iupmin,iupmax 
       DO site1=1,Nc
        IF(UMASK(site1,site1).AND.BTEST(sector%up%state(iup),IMPiorbup(site1)-1))THEN
          DO ido=idomin,idomax 
            IF(BTEST(sector%down%state(ido),IMPiorbdo(site1)-1))THEN 
              istate          = iup + (ido-1) * sec  
              vec_out(istate) = vec_out(istate) + U(site1,site1) * vec_in(istate)
            ENDIF
          ENDDO
        ENDIF
       ENDDO
     ENDDO
    else
     DO iup=iupmin,iupmax
       DO site1=1,Nc
        IF(UMASK(site1,site1).AND.BTEST(sector%up%state(iup),IMPiorbup(site1)-1))THEN
          DO ido=idomin,idomax
            IF(BTEST(sector%down%state(ido),IMPiorbdo(site1)-1))THEN
              istate          = iup + (ido-1) * sec
              vec_out(istate) = vec_out(istate) + quench_U * vec_in(istate)
            ENDIF
          ENDDO
        ENDIF
       ENDDO
     ENDDO
    endif

    IF(offdiag_coulomb)THEN
     DO iup=iupmin,iupmax
       DO ido=idomin,idomax 
         DO site1=1,Nc-1
           bn1 = (/BTEST(sector%up%state(iup),IMPiorbup(site1)-1),BTEST(sector%down%state(ido),IMPiorbdo(site1)-1)/) 
           n1  = COUNT(bn1)
           IF(n1/=0)THEN
             DO site2=site1+1,Nc
              IF(UMASK(site1,site2))THEN
                bn2 = (/BTEST(sector%up%state(iup),IMPiorbup(site2)-1),BTEST(sector%down%state(ido),IMPiorbdo(site2)-1)/)
                n2  = COUNT(bn2)
                IF(n2/=0)THEN
                  istate = iup + (ido-1) * sec  
                  vec_out(istate) = vec_out(istate) + U(site1,site2) * n1 * n2 * vec_in(istate)

                  if(Jhund_Slater_type)then
                   if(bn1(1).and.bn2(1))then
                     vec_out(istate) = vec_out(istate) - JJmatrix(site1,site2) * vec_in(istate)
                   endif
                   if(bn1(2).and.bn2(2))then
                     vec_out(istate) = vec_out(istate) - JJmatrix(site1,site2) * vec_in(istate)
                   endif
                  endif

                ENDIF
              ENDIF
             ENDDO
           ENDIF
         ENDDO
       ENDDO
     ENDDO
    ENDIF

    if(OPEN_MP.and.go_for_omp)then
#ifndef OPENMP_MPI
     iupmin0=1
#else
     iupmin0=min_
#endif
    else
     iupmin0=iupmin
    endif

    DO iup=iupmin,iupmax
      iuploc = iup - iupmin0 + 1
      DO ido=idomin,idomax
       idoloc          =  ido - idomin + 1
       istate          =  iup + (ido-1) * sec  
       vec_out(istate) =  vec_out(istate) + ( diagup(iuploc) + diagdo(idoloc) ) * vec_in(istate)
      ENDDO
    ENDDO

    IF(OPEN_MP)then
#ifndef OPENMP_MPI
     if(iupmin>1)then
      nhoffdiagup = long_sum(noffup(1:iupmin-1))
#else
     if(iupmin>min_)then
      nhoffdiagup = long_sum(noffup(min_:iupmin-1))
#endif
     else
      nhoffdiagup = 0
     endif
    else
     nhoffdiagup = 0
    endif

    DO iup=iupmin,iupmax 
      istatemin  = iup +(idomin-1)*sec 
      istatemax  = iup +(idomax-1)*sec 
      noff       = noffup(iup-iupmin0+1)
      DO irank=1,noff
       hoffdiagup  = offdiagup(nhoffdiagup+irank)
       choffdiagup = conjg(hoffdiagup)
       jup         = rankoffup(nhoffdiagup+irank)
       jstate      = jup +(idomin-1)*sec  
       DO istate=istatemin,istatemax,stridedo
         if(.not.USE_CC) &
&        vec_out(jstate) = vec_out(jstate) +  hoffdiagup * vec_in(istate)
         vec_out(istate) = vec_out(istate) + choffdiagup * vec_in(jstate)
         jstate = jstate + stridedo
       ENDDO
      ENDDO
      nhoffdiagup = nhoffdiagup + noff
    ENDDO

    nhoffdiagdo = 0
    DO ido=idomin,idomax ! span the local chunk of the ndo-sector
      istatemin = iupmin +(ido-1)*sec  
      istatemax = iupmax +(ido-1)*sec  
      noff      = noffdo(ido-idomin+1)
      DO irank=1,noff
       hoffdiagdo  = offdiagdo(nhoffdiagdo+irank)
       choffdiagdo = conjg(hoffdiagdo)
       jdo         = rankoffdo(nhoffdiagdo+irank)
       jstate      = iupmin +(jdo-1)*sec 
       DO istate=istatemin,istatemax,strideup
         if(.not.USE_CC) &
&        vec_out(jstate) = vec_out(jstate) +  hoffdiagdo * vec_in(istate)
         vec_out(istate) = vec_out(istate) + choffdiagdo * vec_in(jstate)
         jstate = jstate + strideup
       ENDDO
      ENDDO
      nhoffdiagdo = nhoffdiagdo + noff 
    ENDDO


      if(flag_slater_int.and..not.use_precomputed_slater_matrix)then
        DO iup=iupmin,iupmax
          DO ido=idomin,idomax
              istate_  = iup + (ido-1) * sec
              up_in_   = sector%up%state(iup)
              do_in_   = sector%down%state(ido)

              DO site4=1,Nc
              bn4 = (/BTEST(up_in_,IMPiorbup(site4)-1),BTEST(do_in_,IMPiorbdo(site4)-1)/)
              do is2=1,2
                if(.not.bn4(is2)) cycle
                if(is2==1)then
                 up_out_ = IBCLR(up_in_,IMPiorbup(site4)-1)
                 do_out_ = do_in_
                 fermion_sign_=1
                 i2=IMPiorbup(site4)
                 DO jj=1,i2-1; IF(BTEST(up_in_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                else
                 do_out_ = IBCLR(do_in_,IMPiorbdo(site4)-1)
                 up_out_ = up_in_
                 fermion_sign_=1
                 i2=IMPiorbdo(site4)
                 DO jj=1,i2-1; IF(BTEST(do_in_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                endif

                DO site3=1,Nc
                 bn3 = (/BTEST(up_out_,IMPiorbup(site3)-1),BTEST(do_out_,IMPiorbdo(site3)-1)/)
                 do is1=1,2
                   if(.not.bn3(is1)) cycle
                   if(is1==1)then
                    up_outb_ = IBCLR(up_out_,IMPiorbup(site3)-1)
                    do_outb_ = do_out_
                   else
                    do_outb_ = IBCLR(do_out_,IMPiorbdo(site3)-1)
                    up_outb_ = up_out_
                   endif

                   DO site2=1,Nc
                      bn2 = (/BTEST(up_outb_,IMPiorbup(site2)-1),BTEST(do_outb_,IMPiorbdo(site2)-1)/)
                      if(bn2(is1)) cycle
                      if(is1==1)then
                       up_outbb_ = IBSET(up_outb_,IMPiorbup(site2)-1)
                       do_outbb_ = do_outb_
                       fermion_signbb_=fermion_sign_
                       i1=min(IMPiorbup(site2),IMPiorbup(site3)); i2=max(IMPiorbup(site2),IMPiorbup(site3))
                       DO jj=i1+1,i2-1; IF(BTEST(up_outb_,jj-1)) fermion_signbb_=-fermion_signbb_; ENDDO
                      else
                       up_outbb_ = up_outb_
                       do_outbb_ = IBSET(do_outb_,IMPiorbdo(site2)-1)
                       fermion_signbb_=fermion_sign_
                       i1=min(IMPiorbdo(site2),IMPiorbdo(site3)); i2=max(IMPiorbdo(site2),IMPiorbdo(site3))
                       DO jj=i1+1,i2-1; IF(BTEST(do_outb_,jj-1)) fermion_signbb_=-fermion_signbb_; ENDDO
                      endif

                      DO site1=1,Nc
                        bn1 = (/BTEST(up_outbb_,IMPiorbup(site1)-1),BTEST(do_outbb_,IMPiorbdo(site1)-1)/)
                        if(bn1(is2)) cycle
                        if(is2==1)then
                         up_outbbb_ = IBSET(up_outbb_,IMPiorbup(site1)-1)
                         do_outbbb_ = do_outbb_
                         fermion_signbbb_=fermion_signbb_
                         i2=IMPiorbup(site1)
                         DO jj=1,i2-1; IF(BTEST(up_outbb_,jj-1)) fermion_signbbb_=-fermion_signbbb_; ENDDO
                        else
                         up_outbbb_ = up_outbb_
                         do_outbbb_ = IBSET(do_outbb_,IMPiorbdo(site1)-1)
                         fermion_signbbb_=fermion_signbb_
                         i2=IMPiorbdo(site1)
                         DO jj=1,i2-1; IF(BTEST(do_outbb_,jj-1)) fermion_signbbb_=-fermion_signbbb_; ENDDO
                        endif

                        jstate_ = sector%up%rank(up_outbbb_) + (sector%down%rank(do_outbbb_)-1)*sec
                        vec_out(istate_) = vec_out(istate_) + vec_in(jstate_)*Slater_Coulomb_c(site1,site2,site3,site4)*dble(fermion_signbbb_)

                       enddo !site1

                   enddo !site2

                 enddo !is1
                enddo !site3

                enddo !is2
                ENDDO !site4

          ENDDO
        ENDDO

      endif

     if(flag_slater_int.and.use_precomputed_slater_matrix)then
         DO iup=iupmin,iupmax
          DO ido=idomin,idomax
              istate_  = iup + (ido-1) * sec
              up_in_   = sector%up%state(iup)
              do_in_   = sector%down%state(ido)
              BATHup   =  MOD(up_in_,Nb)
              BATHdo   =  MOD(do_in_,Nb)
              iupimp   = (up_in_ - BATHup)/Nb
              idoimp   = (do_in_ - BATHdo)/Nb
              do site1=1,UCC(iupimp,idoimp,1,3)
                 iupimp_out = UCC(iupimp,idoimp,site1,1)
                 idoimp_out = UCC(iupimp,idoimp,site1,2)
                 AIMup      = BATHup + iupimp_out
                 AIMdo      = BATHdo + idoimp_out
                 jstate_    = sector%up%rank(AIMup) + (sector%down%rank(AIMdo)-1)*sec
                 vec_out(istate_) = vec_out(istate_) + vec_in(jstate_)*UCCc(iupimp,idoimp,site1)
              enddo
          ENDDO
         ENDDO 
      endif

      IF(abs(Jhund)>1.d-5.and.offdiag_coulomb)THEN
       !this term is -2 x Jhund x Si Sj (spin operators)
        DO iup=iupmin,iupmax
          DO ido=idomin,idomax
            istate_ = iup + (ido-1) * sec
            DO site1=1,Nc
              bn1 = (/BTEST(sector%up%state(iup),IMPiorbup(site1)-1),BTEST(sector%down%state(ido),IMPiorbdo(site1)-1)/)
              n1  =  COUNT(bn1)
              IF(n1==1)THEN
                DO site2=site1+1,Nc
                 if(UMASK(site1,site2))then
                   bn2= (/BTEST(sector%up%state(iup),IMPiorbup(site2)-1),BTEST(sector%down%state(ido),IMPiorbdo(site2)-1)/)
                   n2 = COUNT(bn2)
                   IF(n2==1)THEN
                     !-------------------------------------------------------!
                      if(bn1(1).and.bn2(2))then
                        up_in_   = sector%up%state(iup)
                        up_out_  = IBCLR(up_in_,IMPiorbup(site1)-1)
                        do_in_   = sector%down%state(ido)
                        do_out_  = IBCLR(do_in_,IMPiorbdo(site2)-1)
                        do_out_  = IBSET(do_out_,IMPiorbdo(site1)-1)
                        up_out_  = IBSET(up_out_,IMPiorbup(site2)-1)
                        i1=min(IMPiorbup(site1),IMPiorbup(site2)); i2=max(IMPiorbup(site1),IMPiorbup(site2))
                        fermion_sign_ = 1; DO jj=i1+1,i2-1; IF(BTEST(up_out_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                        i1=min(IMPiorbdo(site2),IMPiorbdo(site1)); i2=max(IMPiorbdo(site2),IMPiorbdo(site1))
                        DO jj=i1+1,i2-1; IF(BTEST(do_out_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                        iupb = sector%up%rank(up_out_)
                        idob = sector%down%rank(do_out_)
                        jstate_=iupb +(idob-1)*sec
                        vec_out(istate_) = vec_out(istate_) + vec_in(jstate_)*JJmatrix(site1,site2)*dble(fermion_sign_)
                      endif

                      if(bn1(2).and.bn2(1))then !up-dn
                        up_in_   = sector%up%state(iup)
                        up_out_  = IBCLR(up_in_,IMPiorbup(site2)-1)
                        do_in_   = sector%down%state(ido)
                        do_out_  = IBCLR(do_in_,IMPiorbdo(site1)-1)
                        do_out_  = IBSET(do_out_,IMPiorbdo(site2)-1)
                        up_out_  = IBSET(up_out_,IMPiorbup(site1)-1)
                        i1=min(IMPiorbup(site2),IMPiorbup(site1)); i2=max(IMPiorbup(site2),IMPiorbup(site1))
                        fermion_sign_ = 1; DO jj=i1+1,i2-1; IF(BTEST(up_out_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                        i1=min(IMPiorbdo(site1),IMPiorbdo(site2)); i2=max(IMPiorbdo(site1),IMPiorbdo(site2))
                        DO jj=i1+1,i2-1; IF(BTEST(do_out_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                        iupb = sector%up%rank(up_out_)
                        idob = sector%down%rank(do_out_)
                        jstate_=iupb +(idob-1)*sec
                        vec_out(istate_) = vec_out(istate_) + vec_in(jstate_)*JJmatrix(site1,site2)*dble(fermion_sign_)
                      endif
                      if(.not.Jhund_Slater_type)then
                       if((bn1(1).and.bn2(1)).or.(bn1(2).and.bn2(2)))then !up-dn
                      !Si^z Sj^z : up up
                         vec_out(istate_) = vec_out(istate_) - 0.50d0 * JJmatrix(site1,site2) * vec_in(istate_)
                       else 
                         vec_out(istate_) = vec_out(istate_) + 0.50d0 * JJmatrix(site1,site2) * vec_in(istate_)
                       endif
                      endif
                     !-------------------------------------------------------!
                   ENDIF

                  endif
                ENDDO
              ENDIF
             ENDDO

          ENDDO
        ENDDO
      ENDIF

      IF(abs(Jhund)>1.d-5.and.offdiag_coulomb)THEN
        if(Jhund<0.)then
         write(*,*) 'Weird, hund coupling is negative'
         stop
        endif
       !this term is Jhund x double_hopping(up+dn) from site1 to site2
        DO iup=iupmin,iupmax
          DO ido=idomin,idomax
            istate_ = iup + (ido-1) * sec
            DO site1=1,Nc
              bn1 = (/BTEST(sector%up%state(iup),IMPiorbup(site1)-1),BTEST(sector%down%state(ido),IMPiorbdo(site1)-1)/)
              n1  =  COUNT(bn1)
              IF(n1==2)THEN
                DO site2=1,Nc
                 if(site1/=site2.and.UMASK(site1,site2))then
                   bn2= (/BTEST(sector%up%state(iup),IMPiorbup(site2)-1),BTEST(sector%down%state(ido),IMPiorbdo(site2)-1)/)
                   n2 = COUNT(bn2)
                   IF(n2==0)THEN
                   !-------------------------------------------------------!
                        up_in_   = sector%up%state(iup)
                        up_out_  = IBCLR(up_in_,IMPiorbup(site1)-1)
                        do_in_   = sector%down%state(ido)
                        do_out_  = IBCLR(do_in_,IMPiorbdo(site1)-1)
                        do_out_  = IBSET(do_out_,IMPiorbdo(site2)-1)
                        up_out_  = IBSET(up_out_,IMPiorbup(site2)-1)
                        i1=min(IMPiorbup(site1),IMPiorbup(site2)); i2=max(IMPiorbup(site1),IMPiorbup(site2))
                        fermion_sign_ = 1; DO jj=i1+1,i2-1; IF(BTEST(up_out_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                        i1=min(IMPiorbdo(site1),IMPiorbdo(site2)); i2=max(IMPiorbdo(site1),IMPiorbdo(site2))
                        DO jj=i1+1,i2-1; IF(BTEST(do_out_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                        iupb = sector%up%rank(up_out_)
                        idob = sector%down%rank(do_out_)
                        jstate_=iupb +(idob-1)*sec
                        vec_out(istate_) = vec_out(istate_) + vec_in(jstate_)*JJmatrix(site1,site2)*dble(fermion_sign_)
                   !-------------------------------------------------------!
                   ENDIF
                  endif
                ENDDO
              ENDIF
             ENDDO
          ENDDO
        ENDDO
      ENDIF

     1193 continue

   !$OMP END PARALLEL

35  continue
    
    if(HILBERT_SPACE_SPLITED_AMONG_NODES) then
     if(FLAG_MPI_GREENS>0) then
      write(*,*) 'you try to split the hilbert space amongst nodes, but also'
      write(*,*) 'at the same time you are parallelizing the computation of green functions'
      stop 'critical'
     endif
     if(USE_CC.and..not.emptychunks)then
#ifdef debug
      if(size(vec_out)/=dimen)then
       write(*,*) 'ERROR - oups vec_out has wrong dimension'
       write(*,*) ' dim vec_out : ', size(vec_out)
       write(*,*) ' dim sector  : ', dimen
       stop
      endif
      if(sum(sector%chunk)/=dimen)then
       write(*,*) 'ERROR MPI ED : sum of chuncks different from dimension of sector'
       write(*,*) ' chuncks  : ', sector%chunk(1:j-1)
       write(*,*) 'dimension : ', dimen
       stop
      endif
#endif
      call mpibcast(vec_out,sector%chunk,[(sum(sector%chunk(1:j-1)),j=1,size2)]   )
     else
      if(allocated(vec_tot_out)) deallocate(vec_tot_out)
      allocate(vec_tot_out(size(vec_out)))
      vec_tot_out=0.
      call MPI_ALLREDUCE(vec_out,vec_tot_out,dimen,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      vec_out=vec_tot_out
      if(allocated(vec_tot_out)) deallocate(vec_tot_out)
     endif
    endif

  if(allocated(imin_)) deallocate(imin_,imax_)

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
!*********************************************
!*********************************************

  SUBROUTINE HAIMupdo_multr(vec_out,vec_in)

  !-------------------------------------------------------!
  ! WE COMPUTE THE RELEVANT CHUNK OF vec_out = H * vec_in !
  !-------------------------------------------------------!
    IMPLICIT NONE
    REAL(DBL),ALLOCATABLE       :: vec_tot_out(:)
    REAL(DBL)                   :: vec_in(:)
    REAL(DBL)                   :: vec_out(:)
    REAL(DBL)                   :: hoffdiagup,choffdiagup,hoffdiagdo,choffdiagdo 
    INTEGER                     :: site1,site2,n1,n2,j
    INTEGER                     :: istate,   istatemin,   istatemax
    INTEGER                     :: jstate,   jstatemin,   jstatemax
    INTEGER                     :: istateloc,istateminloc,istatemaxloc
    INTEGER                     :: iup,ido,jup,jdo,iuploc,idoloc,irank
    INTEGER                     :: noff,I,sec,TID,iupmin0,min_,max_
    INTEGER(4),allocatable      :: imin_(:),imax_(:)
    INTEGER(8)                  :: nhoffdiagup,nhoffdiagdo
    LOGICAL                     :: bn(2),bn1(2),bn2(2)
    INTEGER                     :: norbs_,jj
    TYPE(fermion_ket_type)      :: up_in,do_in,up_out,do_out
    INTEGER                     :: iupb,idob,istate_,jstate_
    INTEGER                     :: up_in_,do_in_,i1,i2,do_out_,up_out_,fermion_sign_
    integer                     :: fermion_signb_,site3,site4,n3,n4,is1,is2,up_outb_,do_outb_
    logical                     :: bn3(2),bn4(2),go_for_omp
    integer                     :: up_outbb_, do_outbb_, up_outbbb_, do_outbbb_, fermion_signbb_, fermion_signbbb_
    integer                     :: IMPstate,BATHup,BATHdo,iupimp,idoimp,iupimp_out,idoimp_out,AIMup,AIMdo
    logical                     :: emptychunks

    if(verboseall) write(*,*) 'USING OPENMP - RANK : ', OPEN_MP,rank

    if(dimen==1)then
     vec_out(1)=1
     return
    endif

    if(OPEN_MP)then
    if(verboseall)then
     write(*,*) 'ALLOCATING IMIN/IMAX WITH SIZE - RANK : ', MAXT, RANK
     write(*,*) 'IUPMIN IUPMAX                         : ', iupmin,iupmax
     endif
     allocate(imin_(MAXT),imax_(MAXT))
    endif

    sec = sector%up%dimen

    if(sec==0) write(*,*) 'danger in hmult up do dim up is 0'
    if(iupmax<iupmin) write(*,*) 'WARNING  iupmax<iupmin'

    
    emptychunks=any(sector%chunk==0)

    if(HILBERT_SPACE_SPLITED_AMONG_NODES.and.USE_CC) then
      if(.not.emptychunks)then
         vec_out( 1+sum(sector%chunk(1:iproc-1)) : sum(sector%chunk(1:iproc)) ) = 0.d0
      else
         vec_out=0.d0
      endif
     else
      vec_out=0.d0
    endif

    if(sector%chunk(iproc)==0) then
       write(*,*) 'EMPTY SECTOR - SIZE VEC_OUT : ', size(vec_out)
       goto 35
    endif

    min_=iupmin; max_=iupmax

#ifndef OPENMP_MPI
   if(OPEN_MP) then
     call openmp_split_array(sec,imin_,imax_) 
   endif
#else
   if(OPEN_MP) then
     call openmp_split_array(iupmax-iupmin+1,imin_,imax_)
     where(imin_/=0) imin_=imin_+iupmin-1
     where(imax_/=0) imax_=imax_+iupmin-1
     if(verboseall)then
     do i=1,MAXT
     write(*,*) 'CHUNK : ', imin_(i),imax_(i)
     enddo
     endif
   endif
#endif

                 go_for_omp=.false.
     if(OPEN_MP) go_for_omp=minval(imax_)>0.and.minval(sector%chunk)>0.and.minval(imin_)>0

    !$OMP PARALLEL  IF(go_for_omp)  PRIVATE( IMPstate,BATHup,BATHdo,iupimp,idoimp,iupimp_out,idoimp_out,AIMup,AIMdo,up_outbb_, do_outbb_, up_outbbb_, do_outbbb_, fermion_signbb_, fermion_signbbb_,fermion_signb_,site3,site4,n3,n4,is1,is2,up_outb_,do_outb_,bn3,bn4,n1,n2,iupb,idob,iupmin0,iup,iupmin,iupmax,ido,istate,TID,site1,bn,bn1,bn2,site2,iuploc,idoloc,istatemin,istatemax,nhoffdiagup,nhoffdiagdo,noff,irank,hoffdiagup,choffdiagup,hoffdiagdo,choffdiagdo,jdo,jup,jstate,up_in_,do_in_,i1,i2,jj,do_out_,up_out_,fermion_sign_,istate_,jstate_) SHARED(vec_out,vec_in,sec,diagup,diagdo)
   if(OPEN_MP.and.go_for_omp) then
     TID=OMP_GET_THREAD_NUM()+1; iupmin=imin_(TID); iupmax=imax_(TID); 
   else
     iupmin=min_; iupmax=max_
   endif

   if(OPEN_MP.and.go_for_omp)then
   if(OMP_GET_NUM_THREADS()/=MAXT)Then
     write(*,*) 'ERROR obtained number of threads is different from MAXT'
     write(*,*) 'MAXT = ', MAXT
     write(*,*) 'NUMBER OF THREADS = ', OMP_GET_NUM_THREADS()
     stop
   endif
   endif

    if(iupmax==0) goto 1192

    DO iup=iupmin,iupmax 
      DO site1=1,Nc
      !call check_openmp 
       IF(UMASK(site1,site1).AND.BTEST(sector%up%state(iup),IMPiorbup(site1)-1))THEN
         DO ido=idomin,idomax 
           IF(BTEST(sector%down%state(ido),IMPiorbdo(site1)-1))THEN
             istate          = iup + (ido-1) * sec  
             vec_out(istate) = vec_out(istate) + U(site1,site1) * vec_in(istate)
           ENDIF
         ENDDO
       ENDIF
      ENDDO
    ENDDO

    IF(offdiag_coulomb)THEN
     DO iup=iupmin,iupmax
       DO ido=idomin,idomax
         DO site1=1,Nc-1
           bn1 = (/BTEST(sector%up%state(iup),IMPiorbup(site1)-1),BTEST(sector%down%state(ido),IMPiorbdo(site1)-1)/) 
           n1  = COUNT(bn1)
           IF(n1/=0)THEN
             DO site2=site1+1,Nc
              IF(UMASK(site1,site2))THEN
                bn2 = (/BTEST(sector%up%state(iup),IMPiorbup(site2)-1),BTEST(sector%down%state(ido),IMPiorbdo(site2)-1)/)
                n2  = COUNT(bn2)
                IF(n2/=0)THEN
                  istate = iup + (ido-1) * sec  
                  vec_out(istate) = vec_out(istate) + U(site1,site2) * n1 * n2 * vec_in(istate)

                  if(Jhund_Slater_type)then
                   if(bn1(1).and.bn2(1))then
                     vec_out(istate) = vec_out(istate) - JJmatrix(site1,site2) * vec_in(istate)
                   endif
                   if(bn1(2).and.bn2(2))then
                     vec_out(istate) = vec_out(istate) - JJmatrix(site1,site2) * vec_in(istate)
                   endif
                  endif

                ENDIF
              ENDIF
             ENDDO
           ENDIF
         ENDDO
       ENDDO
     ENDDO
    ENDIF

    if(OPEN_MP.and.go_for_omp)then
#ifndef OPENMP_MPI
     iupmin0=1
#else
     iupmin0=min_
#endif
    else
     iupmin0=iupmin
    endif

    DO iup=iupmin,iupmax 
      iuploc = iup - iupmin0 + 1
      DO ido=idomin,idomax 
       idoloc          =  ido - idomin + 1
       istate          =  iup + (ido-1) * sec  
       vec_out(istate) =  vec_out(istate) + ( diagup(iuploc) + diagdo(idoloc) ) * vec_in(istate)
      ENDDO
    ENDDO

    IF(OPEN_MP)then
#ifndef OPENMP_MPI
     if(iupmin>1)then
      nhoffdiagup = long_sum(noffup(1:iupmin-1))
#else
     if(iupmin>min_)then
      nhoffdiagup = long_sum(noffup(min_:iupmin-1))
#endif
     else
      nhoffdiagup = 0
     endif
    else
     nhoffdiagup = 0
    endif

    DO iup=iupmin,iupmax 
      istatemin = iup + (idomin-1) * sec  
      istatemax = iup + (idomax-1) * sec 
      noff      = noffup(iup-iupmin0+1)
      DO irank=1,noff
       hoffdiagup  = offdiagup(nhoffdiagup+irank)
       jup         = rankoffup(nhoffdiagup+irank)
       jstate      = jup + (idomin-1) * sec  
       DO istate=istatemin,istatemax,stridedo
         if(.not.USE_CC) &
&        vec_out(jstate) = vec_out(jstate) +  hoffdiagup * vec_in(istate)
         vec_out(istate) = vec_out(istate) +  hoffdiagup * vec_in(jstate)
         jstate          = jstate          + stridedo
       ENDDO
      ENDDO
      nhoffdiagup = nhoffdiagup + noff
    ENDDO

    nhoffdiagdo = 0
    DO ido=idomin,idomax 
      istatemin = iupmin + (ido-1) * sec 
      istatemax = iupmax + (ido-1) * sec  
      noff      = noffdo(ido-idomin+1)
      DO irank=1,noff
       hoffdiagdo  = offdiagdo(nhoffdiagdo+irank)
       jdo         = rankoffdo(nhoffdiagdo+irank)
       jstate      = iupmin + (jdo-1) * sec   
       DO istate=istatemin,istatemax,strideup
         if(.not.USE_CC) &
&        vec_out(jstate) = vec_out(jstate) +  hoffdiagdo * vec_in(istate)
         vec_out(istate) = vec_out(istate) +  hoffdiagdo * vec_in(jstate)
         jstate = jstate + strideup
       ENDDO
      ENDDO
      nhoffdiagdo = nhoffdiagdo + noff
    ENDDO




      if(flag_slater_int.and..not.use_precomputed_slater_matrix)then
        DO iup=iupmin,iupmax
          DO ido=idomin,idomax
              istate_  = iup + (ido-1) * sec
              up_in_   = sector%up%state(iup)
              do_in_   = sector%down%state(ido)

              DO site4=1,Nc
              bn4 = (/BTEST(up_in_,IMPiorbup(site4)-1),BTEST(do_in_,IMPiorbdo(site4)-1)/)
              do is2=1,2
                if(.not.bn4(is2)) cycle
                if(is2==1)then
                 up_out_ = IBCLR(up_in_,IMPiorbup(site4)-1)
                 do_out_ = do_in_
                 fermion_sign_=1
                 i2=IMPiorbup(site4)
                 DO jj=1,i2-1; IF(BTEST(up_in_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                else
                 do_out_ = IBCLR(do_in_,IMPiorbdo(site4)-1)
                 up_out_ = up_in_
                 fermion_sign_=1
                 i2=IMPiorbdo(site4)
                 DO jj=1,i2-1; IF(BTEST(do_in_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                endif

                DO site3=1,Nc
                 bn3 = (/BTEST(up_out_,IMPiorbup(site3)-1),BTEST(do_out_,IMPiorbdo(site3)-1)/)
                 do is1=1,2
                   if(.not.bn3(is1)) cycle
                   if(is1==1)then
                    up_outb_ = IBCLR(up_out_,IMPiorbup(site3)-1)
                    do_outb_ = do_out_
                   else
                    do_outb_ = IBCLR(do_out_,IMPiorbdo(site3)-1)
                    up_outb_ = up_out_
                   endif

                   DO site2=1,Nc
                      bn2 = (/BTEST(up_outb_,IMPiorbup(site2)-1),BTEST(do_outb_,IMPiorbdo(site2)-1)/)
                      if(bn2(is1)) cycle
                      if(is1==1)then
                       up_outbb_ = IBSET(up_outb_,IMPiorbup(site2)-1)
                       do_outbb_ = do_outb_
                       fermion_signbb_=fermion_sign_
                       i1=min(IMPiorbup(site2),IMPiorbup(site3)); i2=max(IMPiorbup(site2),IMPiorbup(site3))
                       DO jj=i1+1,i2-1; IF(BTEST(up_outb_,jj-1)) fermion_signbb_=-fermion_signbb_; ENDDO
                      else
                       up_outbb_ = up_outb_
                       do_outbb_ = IBSET(do_outb_,IMPiorbdo(site2)-1)
                       fermion_signbb_=fermion_sign_
                       i1=min(IMPiorbdo(site2),IMPiorbdo(site3)); i2=max(IMPiorbdo(site2),IMPiorbdo(site3))
                       DO jj=i1+1,i2-1; IF(BTEST(do_outb_,jj-1)) fermion_signbb_=-fermion_signbb_; ENDDO
                      endif

                      DO site1=1,Nc
                        bn1 = (/BTEST(up_outbb_,IMPiorbup(site1)-1),BTEST(do_outbb_,IMPiorbdo(site1)-1)/)
                        if(bn1(is2)) cycle
                        if(is2==1)then
                         up_outbbb_ = IBSET(up_outbb_,IMPiorbup(site1)-1)
                         do_outbbb_ = do_outbb_
                         fermion_signbbb_=fermion_signbb_
                         i2=IMPiorbup(site1)
                         DO jj=1,i2-1; IF(BTEST(up_outbb_,jj-1)) fermion_signbbb_=-fermion_signbbb_; ENDDO
                        else
                         up_outbbb_ = up_outbb_
                         do_outbbb_ = IBSET(do_outbb_,IMPiorbdo(site1)-1)
                         fermion_signbbb_=fermion_signbb_
                         i2=IMPiorbdo(site1)
                         DO jj=1,i2-1; IF(BTEST(do_outbb_,jj-1)) fermion_signbbb_=-fermion_signbbb_; ENDDO
                        endif

                        jstate_ = sector%up%rank(up_outbbb_) + (sector%down%rank(do_outbbb_)-1)*sec
                        vec_out(istate_) = vec_out(istate_) + vec_in(jstate_)*Slater_Coulomb_r(site1,site2,site3,site4)*dble(fermion_signbbb_)

                       enddo !site1

                   enddo !site2

                 enddo !is1
                enddo !site3

                enddo !is2
                ENDDO !site4

          ENDDO
        ENDDO

      endif

     if(flag_slater_int.and.use_precomputed_slater_matrix)then
         DO iup=iupmin,iupmax
          DO ido=idomin,idomax
              istate_  = iup + (ido-1) * sec
              up_in_   = sector%up%state(iup)
              do_in_   = sector%down%state(ido)
              BATHup   =  MOD(up_in_,Nb)
              BATHdo   =  MOD(do_in_,Nb)
              iupimp   = (up_in_ - BATHup)/Nb
              idoimp   = (do_in_ - BATHdo)/Nb
              do site1=1,UCC(iupimp,idoimp,1,3)
                 iupimp_out = UCC(iupimp,idoimp,site1,1)
                 idoimp_out = UCC(iupimp,idoimp,site1,2)
                 AIMup      = BATHup +  iupimp_out
                 AIMdo      = BATHdo +  idoimp_out
                 jstate_    = sector%up%rank(AIMup) + (sector%down%rank(AIMdo)-1)*sec
                 vec_out(istate_) = vec_out(istate_) + vec_in(jstate_)*UCCr(iupimp,idoimp,site1)
              enddo
          ENDDO
         ENDDO
      endif


      IF(abs(Jhund)>1.d-5.and.offdiag_coulomb)THEN
       !this term is -2 x Jhund x Si Sj (spin operators)
        DO iup=iupmin,iupmax
          DO ido=idomin,idomax
            istate_ = iup + (ido-1) * sec
            DO site1=1,Nc
              bn1 = (/BTEST(sector%up%state(iup),IMPiorbup(site1)-1),BTEST(sector%down%state(ido),IMPiorbdo(site1)-1)/)
              n1  =  COUNT(bn1)
              IF(n1==1)THEN
                DO site2=site1+1,Nc
                 if(UMASK(site1,site2))then
                   bn2= (/BTEST(sector%up%state(iup),IMPiorbup(site2)-1),BTEST(sector%down%state(ido),IMPiorbdo(site2)-1)/)
                   n2 = COUNT(bn2)
                   IF(n2==1)THEN
                     !-------------------------------------------------------!
                      if(bn1(1).and.bn2(2))then
                        up_in_   = sector%up%state(iup)
                        up_out_  = IBCLR(up_in_,IMPiorbup(site1)-1)
                        do_in_   = sector%down%state(ido)
                        do_out_  = IBCLR(do_in_,IMPiorbdo(site2)-1)
                        do_out_  = IBSET(do_out_,IMPiorbdo(site1)-1)
                        up_out_  = IBSET(up_out_,IMPiorbup(site2)-1)
                        i1=min(IMPiorbup(site1),IMPiorbup(site2)); i2=max(IMPiorbup(site1),IMPiorbup(site2))
                        fermion_sign_ = 1; DO jj=i1+1,i2-1; IF(BTEST(up_out_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                        i1=min(IMPiorbdo(site2),IMPiorbdo(site1)); i2=max(IMPiorbdo(site2),IMPiorbdo(site1))
                        DO jj=i1+1,i2-1; IF(BTEST(do_out_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                        iupb = sector%up%rank(up_out_)
                        idob = sector%down%rank(do_out_)
                        jstate_=iupb +(idob-1)*sec
                       !vec_out(jstate_) = vec_out(jstate_) + vec_in(istate_)*JJmatrix(site1,site2)*dble(fermion_sign_)
                        vec_out(istate_) = vec_out(istate_) + vec_in(jstate_)*JJmatrix(site1,site2)*dble(fermion_sign_) !C.C.
                      endif

                      if(bn1(2).and.bn2(1))then !up-dn
                        up_in_   = sector%up%state(iup)
                        up_out_  = IBCLR(up_in_,IMPiorbup(site2)-1)
                        do_in_   = sector%down%state(ido)
                        do_out_  = IBCLR(do_in_,IMPiorbdo(site1)-1)
                        do_out_  = IBSET(do_out_,IMPiorbdo(site2)-1)
                        up_out_  = IBSET(up_out_,IMPiorbup(site1)-1)
                        i1=min(IMPiorbup(site2),IMPiorbup(site1)); i2=max(IMPiorbup(site2),IMPiorbup(site1))
                        fermion_sign_ = 1; DO jj=i1+1,i2-1; IF(BTEST(up_out_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                        i1=min(IMPiorbdo(site1),IMPiorbdo(site2)); i2=max(IMPiorbdo(site1),IMPiorbdo(site2))
                        DO jj=i1+1,i2-1; IF(BTEST(do_out_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                        iupb = sector%up%rank(up_out_)
                        idob = sector%down%rank(do_out_)
                        jstate_=iupb +(idob-1)*sec
                       !vec_out(jstate_) = vec_out(jstate_) + vec_in(istate_)*JJmatrix(site1,site2)*dble(fermion_sign_)
                        vec_out(istate_) = vec_out(istate_) + vec_in(jstate_)*JJmatrix(site1,site2)*dble(fermion_sign_) !C.C.
                      endif
                      if(.not.Jhund_Slater_type)then
                       if((bn1(1).and.bn2(1)).or.(bn1(2).and.bn2(2)))then !up-dn
                      !Si^z Sj^z : up up
                         vec_out(istate_) = vec_out(istate_) - 0.50d0 * JJmatrix(site1,site2) * vec_in(istate_)
                       else 
                         vec_out(istate_) = vec_out(istate_) + 0.50d0 * JJmatrix(site1,site2) * vec_in(istate_)
                       endif
                      endif
                     !-------------------------------------------------------!
                   ENDIF

                  endif
                ENDDO
              ENDIF
             ENDDO

          ENDDO
        ENDDO
      ENDIF

      IF(abs(Jhund)>1.d-5.and.offdiag_coulomb)THEN
        if(Jhund<0.)then
         write(*,*) 'Weird, hund coupling is negative'
         stop
        endif
       !this term is Jhund x double_hopping(up+dn) from site1 to site2
        DO iup=iupmin,iupmax
          DO ido=idomin,idomax
            istate_ = iup + (ido-1) * sec
            DO site1=1,Nc
              bn1 = (/BTEST(sector%up%state(iup),IMPiorbup(site1)-1),BTEST(sector%down%state(ido),IMPiorbdo(site1)-1)/)
              n1  =  COUNT(bn1)
              IF(n1==2)THEN
                DO site2=1,Nc
                 if(site1/=site2.and.UMASK(site1,site2))then
                   bn2= (/BTEST(sector%up%state(iup),IMPiorbup(site2)-1),BTEST(sector%down%state(ido),IMPiorbdo(site2)-1)/)
                   n2 = COUNT(bn2)
                   IF(n2==0)THEN
                   !-------------------------------------------------------!
                        up_in_   = sector%up%state(iup)
                        up_out_  = IBCLR(up_in_,IMPiorbup(site1)-1)
                        do_in_   = sector%down%state(ido)
                        do_out_  = IBCLR(do_in_,IMPiorbdo(site1)-1)
                        do_out_  = IBSET(do_out_,IMPiorbdo(site2)-1)
                        up_out_  = IBSET(up_out_,IMPiorbup(site2)-1)
                        i1=min(IMPiorbup(site1),IMPiorbup(site2)); i2=max(IMPiorbup(site1),IMPiorbup(site2))
                        fermion_sign_ = 1; DO jj=i1+1,i2-1; IF(BTEST(up_out_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                        i1=min(IMPiorbdo(site1),IMPiorbdo(site2)); i2=max(IMPiorbdo(site1),IMPiorbdo(site2))
                        DO jj=i1+1,i2-1; IF(BTEST(do_out_,jj-1)) fermion_sign_=-fermion_sign_; ENDDO
                        iupb = sector%up%rank(up_out_)
                        idob = sector%down%rank(do_out_)
                        jstate_=iupb +(idob-1)*sec
                       !vec_out(jstate_) = vec_out(jstate_) + vec_in(istate_)*JJmatrix(site1,site2)*dble(fermion_sign_)
                        vec_out(istate_) = vec_out(istate_) + vec_in(jstate_)*JJmatrix(site1,site2)*dble(fermion_sign_)
                   !-------------------------------------------------------!
                   ENDIF
                  endif
                ENDDO
              ENDIF
             ENDDO
          ENDDO
        ENDDO
      ENDIF

     1192 continue

!$OMP END PARALLEL

35  continue

    if(HILBERT_SPACE_SPLITED_AMONG_NODES) then
     if(FLAG_MPI_GREENS>0) then
      write(*,*) 'you try to split the hilbert space amongst nodes, but also'
      write(*,*) 'at the same time you are parallelizing the computation of green functions'
      stop 'critical'
     endif
     if(USE_CC.and..not.emptychunks)then
#ifdef debug
      if(size(vec_out)/=dimen)then
       write(*,*) 'ERROR - oups vec_out has wrong dimension'
       write(*,*) ' dim vec_out : ', size(vec_out)
       write(*,*) ' dim sector  : ', dimen
       stop
      endif
      if(sum(sector%chunk)/=dimen)then
       write(*,*) 'ERROR MPI ED : sum of chuncks different from dimension of sector'
       write(*,*) ' chuncks  : ', sector%chunk(1:j-1)
       write(*,*) 'dimension : ', dimen
       stop
      endif
#endif
      call mpibcast(vec_out,sector%chunk,(/(sum(sector%chunk(1:j-1)),j=1,size2)/)   )
     else
      if(allocated(vec_tot_out)) deallocate(vec_tot_out)
      allocate(vec_tot_out(size(vec_out)))
      vec_tot_out=0.
      call MPI_ALLREDUCE(vec_out,vec_tot_out,dimen,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      vec_out=vec_tot_out
      if(allocated(vec_tot_out)) deallocate(vec_tot_out)
     endif
    endif

    if(allocated(imin_)) deallocate(imin_,imax_)

   contains

    !-----------------------------------------------!
    !-----------------------------------------------!
    !-----------------------------------------------!

    subroutine check_openmp
      if(iup>size(sector%up%state))then
       write(*,*) 'running open mp ?     : ', OPEN_MP
       write(*,*) 'size of hilbert space : ', size(sector%up%state),sector%up%dimen
       write(*,*) 'how many threads      : ', MAXT
       write(*,*) 'min_                  : ', min_
       write(*,*) 'max_                  : ', max_
       write(*,*) 'istatemin0,istatemax0 : ', iupmin,iupmax
       stop
      endif
    end subroutine

    !-----------------------------------------------!
    !-----------------------------------------------!
    !-----------------------------------------------!

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
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE HAIMupdo_multr_split(vec_out,vec_in)
  implicit none

  !-------------------------------------------------------!
  ! WE COMPUTE THE RELEVANT CHUNK OF vec_out = H * vec_in !
  !-------------------------------------------------------!

   REAL(DBL)                   :: vec_in(:)
   REAL(DBL)                   :: vec_out(:)
   REAL(DBL)                   :: hoffdiagup,choffdiagup
   INTEGER                     :: site1,site2,n1,n2,j,nup,ndn,range
   INTEGER                     :: istate,   istatemin,   istatemax
   INTEGER                     :: jstate,   jstatemin,   jstatemax
   INTEGER                     :: istateloc,istateminloc,istatemaxloc
   INTEGER                     :: iup,ido,jup,jdo,iuploc,idoloc,irank
   INTEGER                     :: noff,I,sec,imin_(MAXT),imax_(MAXT),TID,iupmin0,min_,max_
   INTEGER(8)                  :: nhoffdiagup,nhoffdiagdo
   LOGICAL                     :: bn(2),bn1(2),bn2(2)

    if(dimen==1)then
     vec_out(1)=1
     return
    endif

   IF(abs(Jhund)>1.d-5)THEN
     write(*,*) 'Jhund not implemented yet for HAIM real split case'
     stop
   endif
 
   vec_out = 0.d0

   if(iupmin/=1) stop 'error use transpose trick for nup ndn basis, the striping of the up spin should be 1'

   DO iup=1,iupmax 
    DO site1=1,Nc
     IF(UMASK(site1,site1).AND.BTEST(sector%up%state(iup),IMPiorbup(site1)-1))THEN
      DO ido=idomin,idomax 
        IF(BTEST(sector%down%state(ido),IMPiorbdo(site1)-1))THEN
          istate          = iup + (ido-idomin) * iupmax  
          vec_out(istate) = vec_out(istate) + U(site1,site1) * vec_in(istate)
        ENDIF
      ENDDO
     ENDIF
    ENDDO
   ENDDO

   IF(offdiag_coulomb)THEN
    DO iup=1,iupmax
       DO ido=idomin,idomax
         DO site1=1,Nc-1
           bn = (/BTEST(sector%up%state(iup),IMPiorbup(site1)-1),BTEST(sector%down%state(ido),IMPiorbdo(site1)-1)/) 
           n1 = COUNT(bn)
           IF(n1/=0)THEN
             DO site2=site1+1,Nc
              IF(UMASK(site1,site2))THEN
                bn = (/BTEST(sector%up%state(iup),IMPiorbup(site2)-1),BTEST(sector%down%state(ido),IMPiorbdo(site2)-1)/)
                n2 = COUNT(bn)
                IF(n2/=0)THEN
                  istate = iup + (ido-idomin) * iupmax 
                  vec_out(istate) = vec_out(istate) + U(site1,site2) * n1 * n2 * vec_in(istate)
                ENDIF
              ENDIF
             ENDDO
           ENDIF
         ENDDO
       ENDDO
     ENDDO
    ENDIF

    DO iup=1,iupmax 
      iuploc = iup
      DO ido=idomin,idomax 
       istate          =  iup + (ido-idomin) * iupmax 
       vec_out(istate) =  vec_out(istate) + ( diagup(iup) + diagdo(ido) ) * vec_in(istate)
      ENDDO
    ENDDO

       nup         =  iupmax
       ndn         =  idomax-idomin+1
       range       =  (ndn-1)*(nup)
       nhoffdiagup =  0

    DO iup=1,nup
         istatemin = iup ; 
         istatemax = iup + range 
         noff      = noffup(iup)
      DO irank=1,noff
       hoffdiagup  = offdiagup(nhoffdiagup+irank)
       jup         = rankoffup(nhoffdiagup+irank)
       jstate      = jup 
       DO istate=istatemin,istatemax,nup
         vec_out(istate) = vec_out(istate) +  hoffdiagup * vec_in(jstate)
         jstate          = jstate          +  nup
       ENDDO
      ENDDO
      nhoffdiagup = nhoffdiagup + noff
    ENDDO
 
    call distributed_memory_transpose_mat(vec_in,  nup , ndn , ndn*size2 ) 
    call distributed_memory_transpose_mat(vec_out, nup , ndn , ndn*size2 )

    nup         =  (idomax-idomin+1)*size2
    ndn         =   iupmax/size2
    range       =  (nup)*(ndn-1)
    nhoffdiagup =  0

    DO iup=1,nup
         istatemin = iup
         istatemax = iup + range
         noff      = noffdo(iup)
      DO irank=1,noff
       hoffdiagup  = offdiagdo(nhoffdiagup+irank)
       jup         = rankoffdo(nhoffdiagup+irank)
       jstate      = jup 
       DO istate=istatemin,istatemax,nup
         vec_out(istate) = vec_out(istate) +  hoffdiagup * vec_in(jstate)
         jstate          = jstate          +  nup
       ENDDO
      ENDDO
      nhoffdiagup = nhoffdiagup + noff
    ENDDO

    call distributed_memory_transpose_mat(vec_in,  nup , ndn , ndn*size2 ) 
    call distributed_memory_transpose_mat(vec_out, nup , ndn , ndn*size2 )

  return
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
!*********************************************
!*********************************************
!*********************************************
!*********************************************

  SUBROUTINE HAIMupdo_multc_split(vec_out,vec_in)

  !-------------------------------------------------------!
  ! WE COMPUTE THE RELEVANT CHUNK OF vec_out = H * vec_in !
  !-------------------------------------------------------!

  COMPLEX(DBL)                :: vec_in(:)
  COMPLEX(DBL)                :: vec_out(:)
  COMPLEX(DBL)                :: hoffdiagup,choffdiagup
  INTEGER                     :: site1,site2,n1,n2,j,nup,ndn,range
  INTEGER                     :: istate,   istatemin,   istatemax
  INTEGER                     :: jstate,   jstatemin,   jstatemax
  INTEGER                     :: istateloc,istateminloc,istatemaxloc
  INTEGER                     :: iup,ido,jup,jdo,iuploc,idoloc,irank
  INTEGER                     :: noff,I,sec,imin_(MAXT),imax_(MAXT),TID,iupmin0,min_,max_
  INTEGER(8)                  :: nhoffdiagup,nhoffdiagdo
  LOGICAL                     :: bn(2)
   

    if(dimen==1)then
     vec_out(1)=1
     return
    endif

   IF(abs(Jhund)>1.d-5)THEN
     write(*,*) 'Jhund not implemented yet for HAIM complex case'
     stop
   endif
 
   vec_out = 0.d0

   if(iupmin/=1) stop 'error use transpose trick for nup ndn basis, the striping of the up spin should be 1'

   DO iup=1,iupmax 
    DO site1=1,Nc
     IF(UMASK(site1,site1).AND.BTEST(sector%up%state(iup),IMPiorbup(site1)-1))THEN
      DO ido=idomin,idomax 
        IF(BTEST(sector%down%state(ido),IMPiorbdo(site1)-1))THEN
          istate          = iup + (ido-idomin) * iupmax  
          vec_out(istate) = vec_out(istate) + U(site1,site1) * vec_in(istate)
        ENDIF
      ENDDO
     ENDIF
    ENDDO
   ENDDO

   IF(offdiag_coulomb)THEN
    DO iup=1,iupmax
       DO ido=idomin,idomax
         DO site1=1,Nc-1
           bn = (/BTEST(sector%up%state(iup),IMPiorbup(site1)-1),BTEST(sector%down%state(ido),IMPiorbdo(site1)-1)/) 
           n1 = COUNT(bn)
           IF(n1/=0)THEN
             DO site2=site1+1,Nc
              IF(UMASK(site1,site2))THEN
                bn = (/BTEST(sector%up%state(iup),IMPiorbup(site2)-1),BTEST(sector%down%state(ido),IMPiorbdo(site2)-1)/)
                n2 = COUNT(bn)
                IF(n2/=0)THEN
                  istate = iup + (ido-idomin) * iupmax 
                  vec_out(istate) = vec_out(istate) + U(site1,site2) * n1 * n2 * vec_in(istate)
                ENDIF
              ENDIF
             ENDDO
           ENDIF
         ENDDO
       ENDDO
     ENDDO
    ENDIF

    DO iup=1,iupmax 
      iuploc = iup
      DO ido=idomin,idomax 
       istate          =  iup + (ido-idomin) * iupmax 
       vec_out(istate) =  vec_out(istate) + ( diagup(iup) + diagdo(ido) ) * vec_in(istate)
      ENDDO
    ENDDO

       nup         =  iupmax
       ndn         =  idomax-idomin+1
       range       =  (ndn-1)*(nup)
       nhoffdiagup =  0

    DO iup=1,nup
         istatemin = iup ; 
         istatemax = iup + range 
         noff      = noffup(iup)
      DO irank=1,noff
       hoffdiagup  = offdiagup(nhoffdiagup+irank)
       hoffdiagup  = conjg(hoffdiagup)
       jup         = rankoffup(nhoffdiagup+irank)
       jstate      = jup 
       DO istate=istatemin,istatemax,nup
         vec_out(istate) = vec_out(istate) +  hoffdiagup * vec_in(jstate)
         jstate          = jstate          +  nup
       ENDDO
      ENDDO
      nhoffdiagup = nhoffdiagup + noff
    ENDDO

   !-------------------------------------!
   ! xx transpose of the upxdn matrix xx !
   !-------------------------------------!
 
    call distributed_memory_transpose_mat(vec_in,  nup , ndn , ndn*size2 ) 
    call distributed_memory_transpose_mat(vec_out, nup , ndn , ndn*size2 )

    nup         =  (idomax-idomin+1)*size2
    ndn         =   iupmax/size2
    range       =  (nup)*(ndn-1)
    nhoffdiagup =  0

    DO iup=1,nup
         istatemin = iup
         istatemax = iup + range
         noff      = noffdo(iup)
      DO irank=1,noff
       hoffdiagup  = offdiagdo(nhoffdiagup+irank)
       hoffdiagup  = conjg(hoffdiagup)
       jup         = rankoffdo(nhoffdiagup+irank)
       jstate      = jup 
       DO istate=istatemin,istatemax,nup
         vec_out(istate) = vec_out(istate) +  hoffdiagup * vec_in(jstate)
         jstate          = jstate          +  nup
       ENDDO
      ENDDO
      nhoffdiagup = nhoffdiagup + noff
    ENDDO

    call distributed_memory_transpose_mat(vec_in,  nup , ndn , ndn*size2 ) 
    call distributed_memory_transpose_mat(vec_out, nup , ndn , ndn*size2 )

   !-------------------------------------!
   ! xx transpose of the upxdn matrix xx !
   !-------------------------------------!

  return
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


END MODULE 
