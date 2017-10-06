MODULE eigen_class

  USE rcmatrix_class
 
  IMPLICIT NONE

  REAL(DBL), PARAMETER, PRIVATE    :: zero=0.0_DBL,one=1.0_DBL
  LOGICAL,   PARAMETER, PRIVATE    :: F=.FALSE.,T=.TRUE.

  TYPE eigen_type
    INTEGER             :: rank      = 0                  ! NECESSARY TO IDENTIFY THE EIGENPAIR
    LOGICAL             :: converged = F                  ! TRUE IF CONVERGED
    REAL(DBL)           :: val       = zero               ! Hermitian matrix has real eigenvalues
    TYPE(rcvector_type) :: vec                            ! EIGENVECTOR
    REAL(DBL)           :: lanczos_vecp(1000),rdist 
    INTEGER             :: lanczos_iter
    INTEGER             :: dim_space
    INTEGER             :: dim_sector
  END TYPE 
 
  INTERFACE new_eigen
    MODULE PROCEDURE allocate_new_eigen
    MODULE PROCEDURE new_eigen_from_scratch
    MODULE PROCEDURE copy_eigen
  END INTERFACE
  
  TYPE eigenlist_type
   INTEGER                     ::  neigen = 0
   TYPE(eigen_type), POINTER   ::  eigen(:) => NULL()
  END TYPE 


  INTERFACE add_eigen
    MODULE PROCEDURE add_eigen_,add_eigen__
  END INTERFACE


CONTAINS



!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************


  SUBROUTINE new_eigenlist(list,neigen)
    TYPE(eigenlist_type), INTENT(INOUT) :: list
    INTEGER,              INTENT(IN)    :: neigen
    CALL delete_eigenlist(list)
    list%neigen = neigen
    IF(list%neigen>0) ALLOCATE(list%eigen(neigen))
  END SUBROUTINE 

  
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE remove_eigen(rank,list)
    TYPE(eigenlist_type)                :: list
    INTEGER,              INTENT(IN)    :: rank
    INTEGER                             :: ieigen,jeigen,truerank
    TYPE(eigenlist_type)                :: tmp

    truerank = rank_eigen_in_list(rank,list)

    IF(truerank==0) then
      write(*,*) 'rank : ', rank
      write(*,*) 'i2c  : ', i2c(rank-1)
      write(*,*) 'c2s  : ', c2s(i2c(1))
      write(*,*) "ERROR IN remove_eigen: E"//c2s(i2c(rank-1))//" ISNT IN LIST!"
      stop 'critical, termine'
    endif

    CALL copy_eigenlist(tmp,list)
    CALL  new_eigenlist(list,list%neigen-1)

    jeigen = 0
    DO ieigen=1,tmp%neigen
      IF(truerank/=ieigen)THEN
        jeigen = jeigen + 1
        CALL copy_eigen(list%eigen(jeigen),tmp%eigen(ieigen))
      ENDIF
    ENDDO

    CALL delete_eigenlist(tmp) 

  END SUBROUTINE
  
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

        !-----------------------------!

  SUBROUTINE add_eigen_(eigen,list)
    TYPE(eigenlist_type), INTENT(INOUT) :: list
    TYPE(eigen_type),     INTENT(IN)    :: eigen
    INTEGER                             :: ieigen
    TYPE(eigenlist_type)                :: tmp
    CALL copy_eigenlist(tmp,list)
    CALL  new_eigenlist(list,list%neigen+1)
    DO ieigen=1,tmp%neigen
      CALL copy_eigen(list%eigen(ieigen),tmp%eigen(ieigen))
    ENDDO
    CALL copy_eigen(list%eigen(list%neigen),eigen) 
    CALL delete_eigenlist(tmp) 
  END SUBROUTINE

        !-----------------------------!

  SUBROUTINE add_eigen__(n,VECP,VALP,list)
    TYPE(eigenlist_type), INTENT(INOUT) :: list
    INTEGER                             :: ieigen,n,i
    TYPE(eigenlist_type)                :: tmp
    REAL(8)                             :: VALP(:)
#ifdef _complex
    COMPLEX(DBL)                        :: VECP(:,:)
#else
    REAL(DBL)                           :: VECP(:,:)
#endif
    CALL delete_eigenlist(list)
    list%neigen = n
    IF(n==0) stop 'error add_eigen 0 dimension space'
    ALLOCATE(list%eigen(n))
    do i=1,n
      allocate(list%eigen(i)%vec%rc(size(VECP,1)))
      list%eigen(i)%vec%rc    =  VECP(:,i)
      list%eigen(i)%val       =  VALP(i)
      list%eigen(i)%vec%n     =  size(VECP(:,1))
      list%eigen(i)%converged = .true.
      list%eigen(i)%rank      =  i
      list%eigen(i)%dim_space =  size(VECP(:,1))
    enddo 
  END SUBROUTINE

        !-----------------------------!

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
  
  SUBROUTINE allocate_new_eigen(eigen,n) 
    ! ALLOCATE SINGLE EIGENSTATE
    TYPE(eigen_type), INTENT(INOUT) :: eigen
    INTEGER,          INTENT(IN)    :: n
    CALL delete_eigen(eigen)
    CALL new_rcvector(eigen%vec,n)
  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  subroutine eigen_allocate_vec(eigen,vec,N)
  implicit none
  TYPE(eigen_type),    INTENT(INOUT) :: eigen
  TYPE(rcvector_type), optional      :: vec
  integer,optional                   :: N
   if(present(vec))then
    CALL new_rcvector(eigen%vec,vec)
   else
    if(.not.present(N)) stop 'error allocate eigen vec but no dimension'
    call new_rcvector(eigen%vec,N)
   endif
  end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE new_eigen_from_scratch(eigen,val,vec,CONVERGED,RANK,no_vector) 

    !----------------------------------------!
    ! ALLOCATE & INITIALIZE SINGLE EIGENPAIR !
    !----------------------------------------!

    TYPE(eigen_type),    INTENT(INOUT) :: eigen
    REAL(DBL),           INTENT(IN)    :: val
    TYPE(rcvector_type), INTENT(IN)    :: vec
    LOGICAL, OPTIONAL,   INTENT(IN)    :: CONVERGED,no_vector
    INTEGER, OPTIONAL,   INTENT(IN)    :: RANK
    LOGICAL                            :: no_vec

    CALL delete_eigen(eigen)
    eigen%val = val

                   no_vec=.false.
    if(present(no_vector)) then
     if(no_vector) no_vec=.true.
    endif 
    if(.not.no_vec) CALL new_rcvector(eigen%vec,vec)

    IF(PRESENT(CONVERGED)) eigen%converged = CONVERGED
    IF(PRESENT(RANK))      eigen%rank      = RANK

  END SUBROUTINE


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE copy_eigen(EIGENOUT,EIGENIN) 

    TYPE(eigen_type), INTENT(INOUT) :: EIGENOUT
    TYPE(eigen_type), INTENT(IN)    :: EIGENIN

     CALL delete_eigen(EIGENOUT)
     CALL new_eigen_from_scratch(EIGENOUT,EIGENIN%val,EIGENIN%vec, &
                               & EIGENIN%converged,EIGENIN%rank,no_vector=.not.ASSOCIATED(EIGENIN%vec%rc))

     EIGENOUT%lanczos_vecp  = EIGENIN%lanczos_vecp
     EIGENOUT%rdist         = EIGENIN%rdist
     EIGENOUT%lanczos_iter  = EIGENIN%lanczos_iter
     EIGENOUT%dim_space     = EIGENIN%dim_space 
     EIGENOUT%dim_sector    = EIGENIN%dim_sector

  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE copy_eigenlist(listout,listin)
    TYPE(eigenlist_type), INTENT(INOUT) :: listout
    TYPE(eigenlist_type), INTENT(IN)    :: listin
    INTEGER                             :: ieigen
    CALL new_eigenlist(listout,listin%neigen)  
    DO ieigen=1,listout%neigen
      CALL copy_eigen(listout%eigen(ieigen),listin%eigen(ieigen))
    ENDDO
  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE delete_eigen(eigen)
    TYPE(eigen_type), INTENT(INOUT) :: eigen
    CALL delete_rcvector(eigen%vec)
    eigen%val       = zero
    eigen%rank      = 0
    eigen%converged = F
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE delete_eigenlist(list)
    TYPE(eigenlist_type), INTENT(INOUT) :: list
    INTEGER                             :: ieigen   
    DO ieigen=1,list%neigen
      CALL delete_eigen(list%eigen(ieigen))
    ENDDO
    IF(ASSOCIATED(list%eigen)) DEALLOCATE(list%eigen,STAT=istati)
    list%neigen = 0
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

         !-----------------!

  logical function is_eigen_in_window(E_,window)
  real(DBL)                :: E_,E0,E1
  REAL(DBL), INTENT(IN)    :: window(2)
   E0 = MINVAL(window); E1 = MAXVAL(window)
   is_eigen_in_window = E_>=E0.AND.E_<=E1
  end function

         !-----------------!

  SUBROUTINE filter_eigen(list,window) 
    TYPE(eigenlist_type)     :: list 
    REAL(DBL), INTENT(IN)    :: window(2)
    REAL(DBL)                :: E0,E1
    INTEGER                  :: ieigen,rank
    INTEGER                  :: rank2remove(list%neigen),neigen2remove

    E0 = MINVAL(window,1); E1 = MAXVAL(window,1)

    !--------------------------------------------------------! 
    ! EXPUNGE ALL EIGENSTATES IN list OUTSIDE [E0,E1] window !
    !--------------------------------------------------------!

    rank2remove   = 0 ; neigen2remove = 0

    if(list%neigen>size(list%eigen(:)))then
       write(*,*)  ' list%neigen  :  ' , list%neigen
       write(*,*)  ' size in list :  ' , size(list%eigen)
       write(*,*)  'error, list%neigen is bigger than the actual size of list%eigen , critical error , now stop'
    endif
 
   if(list%neigen>0)then 
    DO ieigen=1,list%neigen
      ! IF EIGENVALUE OUTSIDE window THEN REMOVE EIGENSTATE FROM list 
      IF(list%eigen(ieigen)%val<E0.OR.list%eigen(ieigen)%val>E1)THEN
       rank = list%eigen(ieigen)%rank
       WRITE(log_unit,'(a,f20.12,2(a,E9.3),a)') "# EIGENSTATE E"//c2s(i2c(rank-1))//" : ",list%eigen(ieigen)%val,&
       " : |E"//c2s(i2c(rank-1))//"-E0| = ",list%eigen(ieigen)%val-E0," > ",E1-E0," => EXPUNGE"
       neigen2remove = neigen2remove + 1
       rank2remove(neigen2remove) = rank
      ENDIF
    ENDDO
   endif

   if(neigen2remove>0)then
     DO ieigen=1,neigen2remove
      CALL remove_eigen(rank2remove(ieigen),list)
     ENDDO
   endif

  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  FUNCTION min_eigen(list) 
 
    ! FIND THE LOWEST EIGENVALUE IN list
 
    REAL(DBL)                        :: min_eigen
    TYPE(eigenlist_type), INTENT(IN) :: list
    INTEGER                          :: ieigen
    min_eigen = huge_
    DO ieigen=1,list%neigen
      IF(list%eigen(ieigen)%val<min_eigen) min_eigen = list%eigen(ieigen)%val
    ENDDO

  END FUNCTION 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  FUNCTION sum_boltz(list,beta,E0) 
    !---------------------------!
    ! SUM THE BOLTZMANN WEIGHTS !
    !---------------------------!
    REAL(DBL)                        :: sum_boltz
    REAL(DBL),            INTENT(IN) :: beta,E0
    TYPE(eigenlist_type), INTENT(IN) :: list
    INTEGER                          :: ieigen

    if(.not.FLAG_FULL_ED_GREEN)then
     sum_boltz = SUM((/(DEXPc(-beta*(list%eigen(ieigen)%val-E0)),ieigen=1,list%neigen)/))
    else
     sum_boltz = 0.d0
     do ieigen=1,list%neigen
      if(beta*(list%eigen(ieigen)%val-E0) <= dEmax0 )then
       sum_boltz = sum_boltz + DEXPc(-beta*(list%eigen(ieigen)%val-E0)) 
      endif
     enddo
    endif

  END FUNCTION 

  FUNCTION sum_boltz_times_E(list,beta,E0)
    !---------------------------!
    ! SUM THE BOLTZMANN WEIGHTS !
    !---------------------------!
    REAL(DBL)                        :: sum_boltz_times_E
    REAL(DBL),            INTENT(IN) :: beta,E0
    TYPE(eigenlist_type), INTENT(IN) :: list
    INTEGER                          :: ieigen

    if(.not.FLAG_FULL_ED_GREEN)then
     sum_boltz_times_E = SUM((/((list%eigen(ieigen)%val-E0)*DEXPc(-beta*(list%eigen(ieigen)%val-E0)),ieigen=1,list%neigen)/))
    else
     sum_boltz_times_E = 0.d0
     do ieigen=1,list%neigen
      if(beta*(list%eigen(ieigen)%val-E0) <= dEmax0 )then
       sum_boltz_times_E = sum_boltz_times_E + (list%eigen(ieigen)%val-E0)*DEXPc(-beta*(list%eigen(ieigen)%val-E0))
      endif
     enddo
    endif

  END FUNCTION

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE print_eigenlist(list,UNIT) 
    TYPE(eigenlist_type), INTENT(IN) :: list 
    INTEGER, OPTIONAL,    INTENT(IN) :: UNIT
    INTEGER                          :: ieigen,unit_
    unit_ = log_unit
    IF(PRESENT(UNIT)) unit_ = UNIT
    DO ieigen=1,list%neigen
      WRITE(unit_,'(a,f0.12)') "# E"//c2s(i2c(list%eigen(ieigen)%rank-1))//" = ",list%eigen(ieigen)%val
    ENDDO
  END SUBROUTINE 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE write_raw_eigenlist(list,UNIT) 
    TYPE(eigenlist_type), INTENT(IN) :: list 
    INTEGER,              INTENT(IN) :: UNIT
    INTEGER                          :: ieigen
    DO ieigen=1,list%neigen
      WRITE(UNIT,*) T ! there is another item 
      WRITE(UNIT,*) list%eigen(ieigen)%val
      WRITE(UNIT,*) list%eigen(ieigen)%converged
      WRITE(UNIT,*) list%eigen(ieigen)%rank
      CALL write_raw_rcvector(list%eigen(ieigen)%vec,UNIT)
    ENDDO
    WRITE(UNIT,*) F ! last item
  END SUBROUTINE

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  SUBROUTINE read_raw_eigenlist(list,UNIT) 
    TYPE(eigenlist_type), INTENT(INOUT) :: list 
    INTEGER,             INTENT(IN)     :: UNIT
    TYPE(eigen_type)                    :: eigen
    LOGICAL                             :: new_item
    DO
      READ(UNIT,*) new_item
      IF(.NOT.new_item) EXIT  
      READ(UNIT,*) eigen%val 
      READ(UNIT,*) eigen%converged 
      READ(UNIT,*) eigen%rank 
      CALL read_raw_rcvector(eigen%vec,UNIT)
      CALL add_eigen(eigen,list)
    ENDDO
  END SUBROUTINE 


!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  FUNCTION rank_eigen_in_list(rank,list) RESULT(eigenrank)
    INTEGER                          :: eigenrank
    INTEGER,              INTENT(IN) :: rank
    TYPE(eigenlist_type), INTENT(IN) :: list
    INTEGER                          :: ieigen
    eigenrank = find_rank(rank,(/(list%eigen(ieigen)%rank,ieigen=1,list%neigen)/))
    if(eigenrank==0) then
     write(*,*) 'DID NOT FIND RANK IN EIGENVALUE LIST'
     write(*,*) 'neigen in list   : ', list%neigen
     write(*,*) 'rank list        : ', (/(list%eigen(ieigen)%rank,ieigen=1,list%neigen)/)
     write(*,*) 'looking for rank : ', rank
     write(*,*) 'size list        : ', size(list%eigen)
     stop 'CRITICAL ERROR'
    endif
  END FUNCTION 

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

END MODULE 


