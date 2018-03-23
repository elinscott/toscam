MODULE correlations

  use green_class_computeAB 

  IMPLICIT NONE

  !---------!
  ! PUBLIC  !
  !---------!

  TYPE(correl_type),SAVE  ::  SNAMBU  
  TYPE(correl_type),SAVE  ::  SNAMBUret
  TYPE(green_type) ,SAVE  ::  GNAMBU,GNAMBUN
  TYPE(green_type) ,SAVE  ::  GNAMBUret
  TYPE(green_type) ,SAVE  ::  G(2),GN(2),GKret(2)
  TYPE(green_type) ,SAVE  ::  GF(2),GNF(2)
  TYPE(green_type) ,SAVE  ::  Spm

  !----------!
  ! PRIVATE  !
  !----------!

  TYPE(green_type),  PRIVATE,SAVE           :: Gret(2),GFret(2)
  TYPE(green_type),  PRIVATE,SAVE           :: Sz, N
  TYPE(green_type),  PRIVATE,SAVE           :: Szret,Spmret,Nret
  TYPE(green_type),  PRIVATE,SAVE           :: P3,   P4
  TYPE(green_type),  PRIVATE,SAVE           :: P3ret,P4ret
  TYPE(correl_type), PRIVATE,SAVE           :: CHI,CHIret

  LOGICAL,   PRIVATE                        :: SUPER
  REAL(DBL), PRIVATE                        :: Nc,beta,width
  INTEGER, ALLOCATABLE, PRIVATE             :: triplets(:,:),quadruplets(:,:)
  TYPE(readable_vec_type), POINTER, PRIVATE :: vec_list(:) => NULL()
  REAL(DBL),    PARAMETER, PRIVATE          :: zero=0.0_DBL,one=1.0_DBL,two=2.0_DBL,three=3.0_DBL,four=4.0_DBL
  LOGICAL,      PARAMETER, PRIVATE          :: F=.FALSE.,T=.TRUE.
  LOGICAL, PRIVATE                          :: compute_Sz,compute_Spm,compute_N,compute_P3,compute_P4  
  LOGICAL, PRIVATE                          :: write_cor=.true.,dump_stat=.false.


 CONTAINS

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************

#include "correlations_write.h"
#include "correlations_init.h"

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

  SUBROUTINE compute_correlations(AIM,GS,retard,COMPUTE_ALL_CORREL,only_dens)

    TYPE(AIM_type),             INTENT(IN) :: AIM
    TYPE(eigensectorlist_type), INTENT(IN) :: GS
    LOGICAL                                :: retard
    LOGICAL, OPTIONAL,          INTENT(IN) :: COMPUTE_ALL_CORREL
    INTEGER                                :: ii,Nc,start_correl,bndsup(2),bndsdo(2)
    TYPE(rcmatrix_type)                    :: dmat
    LOGICAL                                :: do_all_correl,only_dens
    LOGICAL, SAVE                          :: first_call = T
    TYPE(rcvector_type)                    :: vec1,vec2 

                                    do_all_correl = F
    IF(PRESENT(COMPUTE_ALL_CORREL)) do_all_correl = COMPUTE_ALL_CORREL

    CALL reset_timer(start_correl)

    Nc     = AIM%Nc
    bndsup = (/1,Nc/)
    bndsdo = bndsup + Nc

    CALL dump_message(TEXT="#")
    CALL dump_message(TEXT="#########################################") 
    CALL dump_message(TEXT="### BEGIN COMPUTATION OF CORRELATIONS ###")
    CALL dump_message(TEXT="#########################################")   
    CALL dump_message(TEXT="#")

    IF(first_call)THEN
       if(.not.ALL_FIRST_CALL) first_call = F
       CALL init_apply_C (AIM)
       CALL init_apply_NS(AIM)
       CALL init_apply_P(triplets,quadruplets,AIM)
    ENDIF

    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!
     write(log_unit,*) ' ==== compute super green ==== '
     call matsubara_super_green

     if(retard.and..not.only_dens)then
       write(log_unit,*) ' ==== compute retarded functions ==== '
       call compute_retarded_functions
     endif

     IF(do_all_correl.and..not.only_dens)THEN
       write(log_unit,*) ' ==== compute dynamic bosonic functions ==== '
       call bosonic_dynamic
       write(log_unit,*) ' ==== compute reduced density matrix ==== '
       call reduced_density_matrix_toolbox
     ELSE
       write(log_unit,*) ' ==== compute equal time bosonic ==== '
       if(.not.only_dens) call bosonic_equal_time
     ENDIF

     CALL write_static(AIM,GS)
    !-----------------------------------------------------------------!
    !-----------------------------------------------------------------!

    CALL timer_fortran(start_correl,"### COMPUTING CORRELATIONS TOOK ###")

  contains

#include "correlations_dyn_stat.h"

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

  SUBROUTINE compute_self_energy(S,G,AIM,retarded)
  IMPLICIT NONE
  TYPE(correl_type), INTENT(INOUT) :: S
  TYPE(correl_type), INTENT(IN)    :: G
  TYPE(AIM_type),    INTENT(IN)    :: AIM
  TYPE(masked_matrix_type)         :: E
  COMPLEX(DBL)                     :: Gm1(G%N,G%N,G%Nw)
  COMPLEX(DBL)                     :: bath1(G%N,G%N,G%Nw)
  INTEGER                          :: iw,nnn
  logical                          :: retarded,block

    CALL bath2hybrid(AIM%bath,S%freq%title)
    CALL Nambu_Ec(E,AIM%impurity%Ec)

    if(retarded)then
     if(S%Nw/=size(AIM%bath%hybridret%fctn,3))then
       write(*,*) 'S%Nw: ', S%Nw
       write(*,*) 'number of frequencies in hybrid : ', size(AIM%bath%hybrid%fctn,3)
       stop 'error compute_self_energy ED solver'
     endif
    else
     if(S%Nw/=size(AIM%bath%hybrid%fctn,3))then
       write(*,*) 'S%Nw: ', S%Nw
       write(*,*) 'number of frequencies in hybrid : ', size(AIM%bath%hybrid%fctn,3)
       stop 'error compute_self_energy ED solver'
     endif
    endif

    nnn    =   size(Gm1,1)
    block  =   maxval(abs(AIM%BATH%Pb%rc%mat))      < cutoff_rvb .and. &
            &  maxval(abs(AIM%BATH%PVbc(1)%rc%mat)) < cutoff_rvb .and. &
            &  maxval(abs(AIM%BATH%PVbc(2)%rc%mat)) < cutoff_rvb

     Gm1=G%fctn
    if(.not.retarded)then
     bath1=AIM%bath%hybrid%fctn
    else
     bath1=AIM%bath%hybridret%fctn
    endif
   
    if(.not.retarded)then
     if(abs(average_G)>=3) call average_correlations(AIM%bath%hybrid,bath1,average_G>=0,MASK_AVERAGE)
    else
     if(abs(average_G)>=3) call average_correlations(AIM%bath%hybridret,bath1,average_G>=0,MASK_AVERAGE)
    endif

     if(abs(average_G)>=2) call average_correlations(G,Gm1,average_G>=0,MASK_AVERAGE)

    DO iw=1,S%Nw
      CALL invmat(size(Gm1,1),Gm1(:,:,iw),block_matrix=.not.AIM%BATH%SUPER.or.block)
    ENDDO

    DO iw=1,S%Nw
      if(.not.retarded)then
        S%fctn(:,:,iw) = SNAMBU%freq%vec(iw)*Id(nnn) -(E%rc%mat + bath1(:,:,iw) + Gm1(:,:,iw))
      else
        if(aimag(SNAMBUret%freq%vec(iw))<0.00001)then
         write(*,*) 'ED, i*delta too small in compute self energy : ', aimag(SNAMBUret%freq%vec(iw))
         stop
        endif
        S%fctn(:,:,iw) = SNAMBUret%freq%vec(iw)*Id(nnn) -(E%rc%mat + bath1(:,:,iw) + Gm1(:,:,iw))
      endif
    ENDDO

   CALL delete_masked_matrix(E)

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

  SUBROUTINE Nambu_green(GNAMBU,G,GF)

    TYPE(green_type), INTENT(INOUT)        :: GNAMBU
    TYPE(green_type), INTENT(IN), OPTIONAL :: G(2)
    TYPE(green_type), INTENT(IN), OPTIONAL :: GF(2)
    INTEGER                                :: Nc,bndsup(2),bndsdo(2),iw,iorb,jorb,ipm,mipm
    TYPE(mask_type)                        :: MASK
    COMPLEX(DBL)                           :: swapvec(GNAMBU%Nw)

#ifdef _complex
    COMPLEX(DBL)                           :: dswapv
#else
    REAL(DBL)                              :: dswapv
#endif

    Nc     = G(1)%N
    bndsup = (/1,Nc/)
    bndsdo = bndsup + Nc

    IF(PRESENT(G))THEN
     
      !-----------------------------------------------------------------------------! 
      ! CAST GREEN'S FUNCTIONS INTO NAMBU FORM I.E.                                 !
      !               ( G1[-,+]  F1[-,-] )                                          !
      ! GNAMBU[-,+] = ( F2[+,+]  G2[+,-] )                                          !  
      !                                                                             !
      ! WITH                                                                        !
      ! G1[-,+](i,j) = < c[i,up](z) * C[j,up] >                                     !
      ! G2[+,-](i,j) = < C[i,do](z) * c[j,do] >                                     !
      ! F1[-,-](i,j) = < c[i,up](z) * c[j,do] >                                     ! 
      ! F2[+,+](i,j) = < C[i,do](z) * C[j,up] > = TRANSPOSE(F1[-,-]) IF H IS REAL   ! 
      ! AND                                                                         ! 
      !               ( G1[+,-]  F1[+,+] )                                          !
      ! GNAMBU[+,-] = ( F2[-,-]  G2[-,+] )                                          !
      !                                                                             !
      ! WITH                                                                        ! 
      ! G1[+,-](i,j) = < C[i,up](z) * c[j,up] >                                     !
      ! G2[-,+](i,j) = < c[i,do](z) * C[j,do] >                                     !
      ! F1[+,+](i,j) = < C[i,up](z) * C[j,do] >                                     !
      ! F2[-,-](i,j) = < c[i,do](z) * c[j,up] > = TRANSPOSE(F1[+,+]) IF H IS REAL   !
      !-----------------------------------------------------------------------------!
 
      DO ipm=1,2
 mipm = 3-ipm

         GNAMBU%correl    (ipm,mipm)%fctn   = zero
         GNAMBU%correlstat(ipm,mipm)%rc%mat = zero

       ! NORMAL PART

       ! SPIN UP
         GNAMBU%correl    (ipm,mipm)%fctn  (bndsup(1):bndsup(2),bndsup(1):bndsup(2),:)  =  G(1)%correl    (ipm,mipm)%fctn   
         GNAMBU%correlstat(ipm,mipm)%rc%mat(bndsup(1):bndsup(2),bndsup(1):bndsup(2))    =  G(1)%correlstat(ipm,mipm)%rc%mat

       ! SPIN DOWN
         GNAMBU%correl    (ipm,mipm)%fctn  (bndsdo(1):bndsdo(2),bndsdo(1):bndsdo(2),:)  =  G(2)%correl    (mipm,ipm)%fctn   
         GNAMBU%correlstat(ipm,mipm)%rc%mat(bndsdo(1):bndsdo(2),bndsdo(1):bndsdo(2))    =  G(2)%correlstat(mipm,ipm)%rc%mat

       ! ANOMALOUS PART

       IF(PRESENT(GF))THEN
         ! UP/DOWN PART
         GNAMBU%correl    (ipm,mipm)%fctn  (bndsup(1):bndsup(2),bndsdo(1):bndsdo(2),:)  = GF(1)%correl    (ipm,ipm)%fctn   
         GNAMBU%correlstat(ipm,mipm)%rc%mat(bndsup(1):bndsup(2),bndsdo(1):bndsdo(2))    = GF(1)%correlstat(ipm,ipm)%rc%mat

       ! DOWN/UP PART

#ifdef _complex
         GNAMBU%correl    (ipm,mipm)%fctn  (bndsdo(1):bndsdo(2),bndsup(1):bndsup(2),:)  = GF(2)%correl    (mipm,mipm)%fctn  
         GNAMBU%correlstat(ipm,mipm)%rc%mat(bndsdo(1):bndsdo(2),bndsup(1):bndsup(2))    = GF(2)%correlstat(mipm,mipm)%rc%mat 
#else
         DO iw=1,GF(1)%Nw
           GNAMBU%correl  (ipm,mipm)%fctn  (bndsdo(1):bndsdo(2),bndsup(1):bndsup(2),iw) = TRANSPOSE(GF(1)%correl(ipm,ipm)%fctn(:,:,iw))
         ENDDO
         GNAMBU%correlstat(ipm,mipm)%rc%mat(bndsdo(1):bndsdo(2),bndsup(1):bndsup(2))    = TRANSPOSE(GF(1)%correlstat(ipm,ipm)%rc%mat)   
#endif

       ENDIF
      ENDDO

    ELSE

      !----------------------------------------------------------------------------------------------------------------------------------!
      ! WE JUST NEED TO PARTICLE-HOLE TRANSFORM THE MATRIX ELEMENTS IN THE LOWER-RIGHT CORNER OF GNAMBU THAT WERE NOT DIRECTLY COMPUTED  !
      !----------------------------------------------------------------------------------------------------------------------------------!

      CALL new_mask(MASK,GNAMBU%correl(2,1)%MM%MASK)

      !-------------------!
      ! NORMAL SPIN DOWN  !
      !-------------------!

      DO iorb=Nc+1,Nc*2
 DO jorb=Nc+1,Nc*2
 IF(.NOT.MASK%mat(iorb,jorb))THEN
       swapvec                                  = GNAMBU%correl(2,1)%fctn(iorb,jorb,:) 
       GNAMBU%correl(2,1)%fctn(iorb,jorb,:)     = GNAMBU%correl(1,2)%fctn(iorb,jorb,:) 
       GNAMBU%correl(1,2)%fctn(iorb,jorb,:)     = swapvec
       dswapv                                   = GNAMBU%correlstat(2,1)%rc%mat(iorb,jorb) 
       GNAMBU%correlstat(2,1)%rc%mat(iorb,jorb) = GNAMBU%correlstat(1,2)%rc%mat(iorb,jorb) 
       GNAMBU%correlstat(1,2)%rc%mat(iorb,jorb) = dswapv
      ENDIF
 ENDDO
 ENDDO

      !--------------------!
      ! ANOMALOUS UP-RIGHT !
      !--------------------!

      DO iorb=1,Nc
 DO jorb=Nc+1,Nc*2
 IF(.NOT.MASK%mat(iorb,jorb))THEN
       swapvec                                  = GNAMBU%correl(2,2)%fctn(iorb,jorb,:) 
       GNAMBU%correl(2,2)%fctn(iorb,jorb,:)     = GNAMBU%correl(1,1)%fctn(iorb,jorb,:) 
       GNAMBU%correl(1,1)%fctn(iorb,jorb,:)     = swapvec
       dswapv                                   = GNAMBU%correlstat(2,2)%rc%mat(iorb,jorb) 
       GNAMBU%correlstat(2,2)%rc%mat(iorb,jorb) = GNAMBU%correlstat(1,1)%rc%mat(iorb,jorb) 
       GNAMBU%correlstat(1,1)%rc%mat(iorb,jorb) = dswapv
      ENDIF
 ENDDO
 ENDDO

      CALL correl2vec(GNAMBU%correl(1,2))
 CALL correl2vec(GNAMBU%correl(2,1))
      CALL masked_matrix2vec(GNAMBU%correlstat(1,2))
 CALL masked_matrix2vec(GNAMBU%correlstat(2,1))
      CALL delete_mask(MASK)

    ENDIF

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

END MODULE 
