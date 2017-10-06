module random

  use linalg
  use genvar
  use common_def 

  real(8),parameter,    private :: rerror=1.d-3
  INTEGER,              PRIVATE :: seedsize
  CHARACTER(LEN=5),     private :: csize
  INTEGER, ALLOCATABLE, PRIVATE :: seed(:)
  CHARACTER(LEN=100),   PRIVATE :: SEEDFILEOUT
  CHARACTER(LEN=100),   PRIVATE :: fmtseed

  INTERFACE randomize_matrix
    module procedure randomize_matrix_r,randomize_matrix_c
  END INTERFACE

  INTERFACE randomize_rank2
    MODULE PROCEDURE randomize_rank2r,randomize_rank2c
  END INTERFACE

  INTERFACE randomize_rank1
    MODULE PROCEDURE randomize_rank1r,randomize_rank1c
  END INTERFACE

  INTERFACE get_random_vec_by_seed
    MODULE PROCEDURE get_random_vec_by_seed_c,get_random_vec_by_seed_r
  END INTERFACE 

  INTERFACE random_vec
    MODULE PROCEDURE random_vec_r,random_vec_c
  END INTERFACE

  INTERFACE randomize_mat
    MODULE PROCEDURE randomize_mat_r,randomize_mat_c
  END INTERFACE

contains


!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

  subroutine initialize_random_numbers(iseed,rank)
  implicit none
  integer :: iseed, rank
   call init_rantab
   call rand_init(iseed=iseed+10*rank,FILEIN='seed')
  end subroutine

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************


  SUBROUTINE write_seed(UNIT,ZESEED)
  implicit none
    INTEGER, OPTIONAL, INTENT(IN)    :: UNIT
    INTEGER, OPTIONAL, INTENT(INOUT) :: ZESEED(:)
    INTEGER                          :: unit_
    CALL RANDOM_SEED(GET=seed)
    IF(PRESENT(UNIT))THEN
      unit_ = UNIT
      CALL dump_message(TEXT="# seed value =",UNIT=unit_)
      IF(PRESENT(ZESEED))THEN
        IF(SIZE(ZESEED)/=seedsize) STOP "ERROR IN write_seed: INCONSISTENT SEED SIZES!"
        WRITE(unit_,fmtseed) ZESEED
      ELSE
        WRITE(unit_,fmtseed) seed
      ENDIF
      CALL flush(unit_)
    ELSE
      CALL open_safe(unit_,SEEDFILEOUT,"UNKNOWN","WRITE",get_unit=.true.)
      IF(PRESENT(ZESEED))THEN
        IF(SIZE(ZESEED)/=seedsize) STOP "ERROR IN write_seed: INCONSISTENT SEED SIZES!"
        WRITE(unit_,fmtseed) ZESEED 
      ELSE
        WRITE(unit_,fmtseed) seed
      ENDIF
      CALL close_safe(unit_)
    ENDIF
  END SUBROUTINE


!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

  SUBROUTINE randomize_rank2r(tab) 
  implicit none
    REAL(DBL),    INTENT(INOUT) :: tab(:,:)
    REAL(DBL),    ALLOCATABLE   :: ritab(:,:)
    INTEGER                     :: rmin,rmax,cmin,cmax
    rmin = LBOUND(tab,1); rmax = UBOUND(tab,1)
    cmin = LBOUND(tab,2); cmax = UBOUND(tab,2)
    ALLOCATE(ritab(rmin:rmax,cmin:cmax))
    ! REAL PART
    CALL RANDOM_NUMBER(ritab)
    tab = ritab
    DEALLOCATE(ritab)
  END SUBROUTINE

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

  SUBROUTINE randomize_rank2c(tab)
  implicit none
    COMPLEX(DBL), INTENT(INOUT)  ::    tab(:,:)
    REAL(DBL),    ALLOCATABLE    ::  ritab(:,:)
    INTEGER                      ::  rmin,rmax,cmin,cmax

    rmin = LBOUND(tab,1); rmax = UBOUND(tab,1)
    cmin = LBOUND(tab,2); cmax = UBOUND(tab,2)
    ALLOCATE(ritab(rmin:rmax,cmin:cmax))
    ! REAL PART
    CALL RANDOM_NUMBER(ritab)
    tab = ritab
    ! IMAGINARY PART
    CALL RANDOM_NUMBER(ritab)
    tab = tab + imi * ritab
    DEALLOCATE(ritab)

  END SUBROUTINE

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

  SUBROUTINE randomize_rank1r(vec) 
  implicit none
    REAL(DBL),    INTENT(INOUT) ::   vec(:)
    REAL(DBL),    ALLOCATABLE   :: rivec(:)
    INTEGER                     :: rmin,rmax,cmin,cmax
    rmin = LBOUND(vec,1); rmax = UBOUND(vec,1)
    ALLOCATE(rivec(rmin:rmax))
    ! REAL PART
    CALL RANDOM_NUMBER(rivec)
    vec = rivec
    DEALLOCATE(rivec)
  END SUBROUTINE

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************


  SUBROUTINE randomize_rank1c(vec)
  implicit none
    COMPLEX(DBL), INTENT(INOUT) :: vec(:)
    REAL(DBL)                   :: rivec(LBOUND(vec,1):UBOUND(vec,1))
    ! REAL PART
    CALL RANDOM_NUMBER(rivec)
    vec = rivec
    ! IMAGINARY PART
    CALL RANDOM_NUMBER(rivec)
    vec = vec + imi * rivec
  END SUBROUTINE


!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************


  SUBROUTINE GENRAN(seed_,randvec,ncpt)
  implicit none
    REAL(DBL), INTENT(INOUT) :: randvec(:)
    INTEGER,   INTENT(IN)    :: ncpt
    INTEGER,   INTENT(INOUT) :: seed_(:)
    ! GENERATES RANDOM VECTOR WITH seed_
    CALL RANDOM_SEED(SIZE=seedsize)
    IF(SIZE(seed_)/=seedsize) STOP "ERROR IN GENRAN: INCONSISTENT SEED SIZE!"
    CALL RANDOM_SEED(PUT=seed_)
    CALL RANDOM_NUMBER(randvec(1:ncpt))
    CALL RANDOM_SEED(GET=seed_)
    IF(.NOT.ALLOCATED(seed)) ALLOCATE(seed(seedsize))
    seed = seed_
    IF(iproc==1) CALL write_seed()
  END SUBROUTINE


!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

  SUBROUTINE build_seed_and_shot_in_genran(seed_,iseed)
  implicit none
    INTEGER   :: seedsize
    INTEGER   :: seed_(:),iseed,i
    CALL RANDOM_SEED(SIZE=seedsize)
    IF(SIZE(seed_)/=seedsize) STOP "ERROR IN build_seed_and_shot_in_genran: INCONSISTENT SEED SIZE!"
    do i=1,seedsize
      seed_(i)=iseed+i
    enddo
    CALL RANDOM_SEED(PUT=seed_)
  END SUBROUTINE

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************


  SUBROUTINE get_random_vec_by_seed_r(randvec,seedd)
  implicit none

    REAL(8),   INTENT(INOUT) :: randvec(:)
    INTEGER                  :: seedd,seedsize,i
    integer,allocatable      :: seed_(:)

    CALL RANDOM_SEED(SIZE=seedsize)
    allocate(seed_(seedsize))
    do i=1,seedsize
      seed_(i)=seedd+i
    enddo
    CALL RANDOM_SEED(PUT=seed_)
    CALL RANDOM_NUMBER(randvec(:))
    deallocate(seed_)

  END SUBROUTINE

       !--------------------!

  SUBROUTINE get_random_vec_by_seed_c(randvec,seedd)
  implicit none
    
    complex(8),INTENT(INOUT) :: randvec(:)
    real(8)                  :: v1(size(randvec)),v2(size(randvec))
    INTEGER                  :: seedd,seedsize,i
    integer,allocatable      :: seed_(:)

    CALL RANDOM_SEED(SIZE=seedsize)
    allocate(seed_(seedsize))
    do i=1,seedsize
      seed_(i)=seedd+i
    enddo
    CALL RANDOM_SEED(PUT=seed_)
    CALL RANDOM_NUMBER(v1(:))
    CALL RANDOM_NUMBER(v2(:))
    randvec=CMPLX(v1,v2,8)
    deallocate(seed_)

  END SUBROUTINE

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

  pure integer function ran_siz()
   ran_siz=seedsize
  end function

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

! Fonction      : Generateur de nombres pseudo-aleatoires
! References    : [1] Numerical Recipes, Fortran 90, p.1142.
! Notes         : Remplace le generateur de P. L'Ecuyer.
        FUNCTION GGUBFS(DSEED)
        IMPLICIT NONE
        INTEGER,PARAMETER               :: K4B=SELECTED_INT_KIND(9)
        INTEGER(K4B),INTENT(INOUT)      :: DSEED
        real(8)                          :: GGUBFS
        INTEGER(K4B),PARAMETER          :: IA=16807,IM=2147483647,IQ=127773,IR=2836
        real(8),SAVE                     :: AM
        INTEGER(K4B),SAVE               :: IX=-1,IY=-1,K
        IF(DSEED<=0.OR.IY<0) THEN
           AM=NEAREST(1.0,-1.0)/IM
           IY=IOR(IEOR(888889999,ABS(DSEED)),1)
           IX=IEOR(777755555,ABS(DSEED))
           DSEED=ABS(DSEED)+1
        ENDIF
        IX=IEOR(IX,ISHFT(IX,13))
        IX=IEOR(IX,ISHFT(IX,-17))
        IX=IEOR(IX,ISHFT(IX,5))
        K=IY/IQ
        IY=IA*(IY-K*IQ)-IR*K
        IF(IY<0) IY=IY+IM
        GGUBFS=AM*IOR(IAND(IM,IEOR(IX,IY)),1)
        END FUNCTION GGUBFS

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

subroutine choose_among_c(array,isite)
implicit none
integer               :: isite,i,j,k,l,siz
complex(8),intent(in) :: array(:)
complex(8)            :: array2(size(array))
real(8)                :: tot,subtot,ran

 siz=size(array)
 tot=real(sum(array))
 array2=real(array)/tot
 ran=drand1()

 subtot=0.d0
 do i=1,siz
  subtot=subtot+array2(i)
  if(subtot>ran)then
   isite=i
   return
  endif
 enddo

end subroutine

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

subroutine random_vec_r(vec)
real(8)    :: vec(:)
integer    :: i,j
do i=1,size(vec)
  vec(i)=dran_tab(i)
enddo
end subroutine

subroutine random_vec_c(vec)
complex(8) :: vec(:)
integer    :: i,j
do i=1,size(vec)
  vec(i)=dran_tabc(i)
enddo
end subroutine

function random_mat_c(n)
integer    :: n
complex(8) :: random_mat_c(n,n)
integer    :: i,j
do i=1,n
 do j=1,n
  random_mat_c(i,j)=dran_tabc(i+j*2)
 enddo
enddo
end function

subroutine randomize_mat_c(mat)
integer    :: n
complex(8) :: mat(:,:)
integer    :: i,j
do i=1,size(mat,1)
 do j=1,size(mat,2)
  mat(i,j)=dran_tabc(i+j*2)
 enddo
enddo
end subroutine

subroutine randomize_mat_r(mat)
integer    :: n
real(8)    :: mat(:,:)
integer    :: i,j
do i=1,size(mat,1)
 do j=1,size(mat,2)
  mat(i,j)=dran_tab(i+j*2)
 enddo
enddo
end subroutine


!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

subroutine uniform_choice_site(site,m)
implicit none
logical :: site(:)
integer :: i,j,k,l,m,N
 N=size(site)
 site=.false.
 j=0
 do
  k=iirand1(N)
  if(.not.site(k))then
   j=j+1
   site(k)=.true.
  endif
  if(j==m)exit
 enddo
end subroutine

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

function drand1()
implicit none
real(8) :: drand1
real(4) :: r
 call random_number(r)
 drand1=dble(r)
end function

 !------------!

function crand1()
implicit none
complex(8) :: crand1
real(8)     :: rtemp1,rtemp2
 rtemp1=(-1.d0+2.d0*drand1())
 rtemp2=(-1.d0+2.d0*drand1())
 crand1=CMPLX(rtemp1,rtemp2,kind=8)
end function

 !------------!

function iirand1(k)
implicit none
integer,intent(in) :: k
integer            :: iirand1
 iirand1=INT(drand1()*dble(k))+1
 if(iirand1>k) iirand1=k
 if(iirand1<1) iirand1=1
return
end function

 !------------!

function iirand1_interval(k1,k2,withoutzero)
implicit none
integer,intent(in) :: k1,k2
integer            :: iirand1_interval
logical            :: withoutzero
 10 iirand1_interval=INT(drand1()*dble(k2-k1+1))
 iirand1_interval=iirand1_interval+k1
 if(iirand1_interval>k2) then
  iirand1_interval=k2
  write(*,*) 'danger iirand1 out of interval'
 endif
 if(iirand1_interval<k1)then
   iirand1_interval=k1
   write(*,*) 'danger iirand1 out of interval'
 endif
 if(withoutzero.and.iirand1_interval==0) goto 10
return
end function


!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

subroutine rand_init(iseed,FILEIN,ZESEED)
implicit none
integer,optional                 :: iseed
CHARACTER(LEN=*),optional        :: FILEIN
INTEGER, OPTIONAL                :: ZESEED(:)
INTEGER                          :: unit_ 
integer                          :: i

    CALL RANDOM_SEED(SIZE=seedsize)
    IF(ALLOCATED(seed))THEN
      IF(SIZE(seed)/=seedsize)THEN
        DEALLOCATE(seed)
        ALLOCATE(seed(seedsize))
      ENDIF
    ELSE
      ALLOCATE(seed(seedsize))
    ENDIF
    WRITE(fmtseed,*) "(",seedsize,"(I0,X))"
    csize = c2s(i2c(seedsize))

   if(present(iseed))then
     seed(1)=iseed
     do i=2,seedsize
       seed(i)=seed(i-1)+1
     enddo
     CALL RANDOM_SEED(PUT=seed)
   endif

   if(present(FILEIN))then
    ! READ SEED
    CALL open_safe(unit_,TRIM(ADJUSTL(FILEIN)),"UNKNOWN","READ",get_unit=.true.)
    READ(unit_,*,END=12,ERR=12) seed
    CALL RANDOM_SEED(PUT=seed)
    12 continue
    CALL close_safe(unit_)
    ! OUTPUT FILE=INPUT FILE
    SEEDFILEOUT = TRIM(ADJUSTL(FILEIN))
   endif

    !
    IF(PRESENT(ZESEED))THEN
      IF(SIZE(ZESEED)/=seedsize) STOP "ERROR IN rand_init: INCONSISTENT SEED SIZES!"
      ZESEED = seed
    ENDIF

end subroutine

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

real(8) function dran_tab(jjj)
implicit none
integer             :: jj
integer, intent(in) :: jjj
  jj       = modi(jjj,size(ran_tab)-1)
  dran_tab = ran_tab(jj)
return
end function

 !======!

complex(8) function dran_tabc(jjj)
implicit none
integer             :: jj
integer, intent(in) :: jjj
  jj=modi(jjj,size(ran_tab-2))
  dran_tabc=CMPLX(-1.d0+2.d0*ran_tab(jj),-1.d0+2.d0*ran_tab(jj+1),kind=8)
return
end function

 !======!

subroutine init_rantab(iseed_)
implicit none
integer          :: i,j,k,l,m
integer,optional :: iseed_
integer          :: iseed

  if(present(iseed_))then
    iseed=iseed_
  else
    iseed=1234567
  endif

  write(*,*) 'initialize RANDOM number with seed: ', iseed
  call rand_init(iseed=iseed)
  if(messages2) write(*,*) ' my rank ' , rank
  do i=0,size(ran_tab)-1
   ran_tab(i)=floor_(drand1(),3)
  enddo
  if(messages2) write(*,*) 'ran_tab initiated',maxval(ran_tab),minval(ran_tab)
  if(messages2) write(*,*) 'my rank : ', rank
  if(messages2) write(*,*) 'iseed   : ', iseed
  call rand_init(iseed=iseed+rank)
  if(messages2) write(*,*) 'my rank :', rank
end subroutine

 !======!

subroutine test_random_number
implicit none
integer :: i
 call init_rantab
 do i=1,10
  write(*,*) 'i,rantab : ', i, ran_tab(i),dran_tab(i)
 enddo 

 do i=1,10
  write(*,*) 'i,drand1 : ', i, drand1()
 enddo
 stop 

end subroutine

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

function RANVECSPHERE()
implicit none
real(8)  :: chi1,chi2,chicarre
real(8)  :: RANVECSPHERE(3),r(3),norm
integer :: j

         20 CONTINUE
         chi1=1.d0-2.d0*drand1()
         chi2=1.d0-2.d0*drand1()
         chicarre = chi1**2 + chi2**2
         if(chicarre<1.d0) then
          r(1)= 2.d0*chi1*dsqrt(1.d0 - chicarre)
          r(2)= 2.d0*chi2*dsqrt(1.d0 - chicarre)
          r(3)= 1.d0 - 2.d0*chicarre
         else
          goto 20
         endif

         norm = norme(r)
         if(norm>1.d-5) then
          r=r/norm
         else
          goto 20
         endif

         RANVECSPHERE(:)=r(:)

         return

         12 norm=0.d0
         do j=1,3
          r(j)=-1.d0 + 2.d0*drand1()
          norm=norm+r(j)**2
         enddo
         if(norm>1.d0) goto 12
         r=r/sqrt(norm)

return
end function

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

function  RANVECCIRCLE()
implicit none
real(8) RANVECCIRCLE(2),dr3,x1,x2,dnorm
real(8) ddCos,ddSin,v2,v1,dr1,dr2
          18 continue
          dr1=drand1()
          dr2=drand1()
          v1=2.d0*dr1-1.d0
          v2=dr2
          dnorm=v1**2+v2**2
          if(dnorm.gt.1.d0) goto 18
          ddCos=(v1**2-v2**2)/dnorm
          ddSin=2.d0*v1*v2/dnorm
          RANVECCIRCLE(1)=ddCos
          RANVECCIRCLE(2)=ddSin
end function

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

function gaussian2()
implicit none
real(8) :: gaussian2,dr3,x1,x2
     19 continue
     x1=-1.d0+2.d0*drand1()
     x2=-1.d0+2.d0*drand1()
     if(x1**2+x2**2>1) goto 19
     gaussian2=DABS(DSQRT(-2.d0*DLOG(x1**2+x2**2)/(x1**2+x2**2))*x1/10.d0)
return
end function

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

function RAN(iseed)
implicit none
real(8)           :: RAN
integer,optional ::iseed
real(4)          :: r
  if(.not.(present(iseed)))then
    call random_number(r)
    RAN=real(r)
  else
    call rand_init(iseed=iseed)
  endif
end function

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

      real(8) FUNCTION gaussian()
      implicit none

  !-----------------------------------------------------------------------!
  ! Returns a normally distributed deviate with zero mean and unit        !
  ! variance. The routine uses the Box-Muller transformation of uniform   !
  ! deviates. For a more efficient implementation of this algorithm,      !
  ! see Press et al., Numerical Recipes, Sec. 7.2.                        !
  !-----------------------------------------------------------------------!

       real(8) :: R, X, Y, PGRAND
 10    X = 2.d0*drand1() - 1.d0
       Y = 2.d0*drand1() - 1.d0
       R = X**2 + Y**2
       IF (R>=1.d0 .or. R<1.d-20) GOTO 10
       gaussian = X*SQRT(-2.d0*LOG(R)/R)
 
     END function

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

subroutine randomvec (a,movfig)
implicit none
integer                          ::  j
real(8), dimension(3),intent(out) ::  a
integer, intent(in)              ::  movfig
 a=0 
 if(movfig==1)a(1) =  1
 if(movfig==2)a(2) =  1
 if(movfig==3)a(3) =  1
 if(movfig==4)a(1) = -1
 if(movfig==5)a(2) = -1
 if(movfig==6)a(3) = -1
return
end subroutine

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

function random_gaussian_vec(n)
implicit none
integer :: n,i
real(8),dimension(n) :: random_gaussian_vec
 do i=1,n
  random_gaussian_vec(i)=gaussian()
 enddo
return
end function

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

subroutine randomdir(hasard,a,movfig)
implicit none
integer                          ::  j
real(8),dimension(3),intent(out) ::  a
real(8),intent(in)               ::  hasard
integer ,intent(out)             ::  movfig
 a=0
 movfig=0
 if(hasard<1./6.)then
 a(1)= 1
 movfig=1
 endif
 if(hasard>=1./6..AND. hasard<2./6.)then
 a(2)= 1
 movfig=2
 endif
 if(hasard>=2./6..AND. hasard<3./6.)then
 a(3)= 1
 movfig=3
 endif
 if(hasard>=3./6..AND. hasard<4./6.)then
 a(1)=-1
 movfig=4
 endif
 if(hasard>=4./6..AND. hasard<5./6.)then
 a(2)=-1
 movfig=5
 endif
 if(hasard>=5./6..AND. hasard<=1.)then
 a(3)=-1
 movfig=6
 endif
 if(movfig==0) movfig=6

return
end subroutine

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

       !---------------!

subroutine randomize_vec(A,amp,flag,realfluc,kk2)
implicit none
integer          :: i,j,k,l,m,siz1,siz2
complex(8)       :: A(:),ctemp
real(8)           :: maxA
real(8),optional  :: amp
logical,optional :: flag,realfluc
integer,optional :: kk2
integer          :: kk

if(present(kk2))then
 kk=kk2
else
 kk=0
endif

maxA=1.d0
if(present(amp)) then
 maxA=maxA*amp
else
 maxA=maxA*rerror
endif
 
siz1=size(A(:))
do i=1,siz1
  if(.not.present(flag).and..not.present(kk2))then
   ctemp=crand1()*maxA
  else
   ctemp=CMPLX(-1.d0+2.d0*dran_tab(i+j*2+50*kk),-1.d0+2.d0*dran_tab(2*i+j+32*kk),kind=8)*maxA
  endif
  if(present(realfluc)) ctemp=real(ctemp)
  A(i)=A(i)+ctemp
enddo
end subroutine

       !---------------!

subroutine randomize_matrix_c(A,amp,flag,kk2)
implicit none
integer          :: i,j,k,l,m,siz1,siz2
complex(8)       :: A(:,:),ctemp
real(8)          :: maxA
real(8),optional :: amp
logical,optional :: flag
integer          :: kk
integer,optional :: kk2

if(present(kk2))then
 kk=kk2
else
 kk=0
endif

maxA=1.d0
if(present(amp)) then
 maxA=maxA*amp
else
 maxA=maxA*rerror
endif

siz1=size(A(:,1))
do i=1,siz1
 do j=i,siz1
  if(.not.present(flag).and..not.present(kk2))then
   ctemp=crand1()*maxA
  else
   ctemp=CMPLX(-1.d0+2.d0*dran_tab(i+j*2+kk),-1.d0+2.d0*dran_tab(2*i+j+42*kk),kind=8)*maxA
  endif
  if(j==i) A(i,i)=A(i,i)+real(ctemp)
  A(i,j)=A(i,j)+ctemp
  A(j,i)=A(j,i)+conjg(ctemp)
 enddo
enddo
end subroutine

       !---------------!

subroutine randomize_matrix_r(A,amp,flag,kk2)
implicit none
integer          :: i,j,k,l,m,siz1,siz2
real(8)           :: A(:,:),rtemp,maxA
real(8),optional  :: amp
logical,optional :: flag
integer          :: kk
integer,optional :: kk2

if(present(kk2))then
 kk=kk2
else
 kk=0
endif

maxA=1.d0
if(present(amp)) then
 maxA=maxA*amp
else
 maxA=maxA*rerror
endif
siz1=size(A(:,1))
do i=1,siz1
 do j=i,siz1
 if(.not.present(flag).and..not.present(kk2))then
  rtemp=(-1.d0+2.d0*drand1())*maxA
 else
  rtemp=(-1.d0+2.d0*dran_tab(i+j*2+kk))*maxA
 endif
  if(j==i) A(i,i)=A(i,i)+rtemp
  A(i,j)=A(i,j)+rtemp
  A(j,i)=A(j,i)+rtemp
 enddo
enddo
end subroutine

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

    real FUNCTION RAN_(idum2)
    IMPLICIT NONE
    INTEGER, PARAMETER :: K4B=selected_int_kind(9)
    INTEGER(K4B),optional :: idum2
    integer,save :: idum
    !  Minimal random number generator of Park and Miller
    ! combined with a Marsaglia shift sequence. Returns a uniform random
    ! deviate between 0.0 and 1.0 (exclusive of the endpoint values). This
    ! fully portable, scalar generator has the traditional (not Fortran 90)
    ! calling sequence with a random deviate as the returned function value:
    ! call with idum a negative integer to initialize; thereafter, do not
    ! alter idum except to reinitialize. The period of this generator is
    ! about 3.1 × 1018.
    INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
    REAL, SAVE :: am
    INTEGER(K4B), SAVE :: ix=-1,iy=-1,k

    if(present(idum2)) idum=idum2

    if (idum <= 0 .or. iy < 0) then    ! Initialize.
       am=nearest(1.0,-1.0)/IM
       iy=ior(ieor(888889999,abs(idum)),1)
       ix=ieor(777755555,abs(idum))
       idum=abs(idum)+1                ! Set idum positive.
    end if
    ix=ieor(ix,ishft(ix,13))           ! Marsaglia shift sequence with period 2^32-1.
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    k=iy/IQ                            ! Park-Miller sequence by Schrage s method, period 2^31-2.
    iy=IA*(iy-k*IQ)-IR*k
    if (iy < 0) iy=iy+IM
    ran_=am*ior(iand(IM,ieor(ix,iy)),1) ! Combine the two generators with masking to ensure nonzero value.
  END FUNCTION 

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

      SUBROUTINE RANDOM_SIMPLE (seed, rnd)
      REAL(8)    a, am, q, r, alo, hi, test
      REAL(8)    rnd, seed
      INTEGER    i, ira
      a = 16807.0d0
      am = 2147483647.0d0
      q = 127773.0d0
      r = 2836.0d0
         hi = int (seed / q)
         alo = seed - q * hi
         test = a * alo - r * hi
         IF (test.gt.0.) then
            seed = test
         ELSE
            seed = test + am
         ENDIF
         rnd = seed / am
      RETURN
      END SUBROUTINE

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

! This random number generator  uses a negative seed           
! It gives numbers in [0,1] interval. Iseed gets updated.               

      REAL(8) FUNCTION RAN32 (IDUM)
      IMPLICIT REAL (8)(A - H, O - Z)
      PARAMETER (IM1 = 2147483563, IM2 = 2147483399, AM = 1.0D0 / IM1,  &
      IMM1 = IM1 - 1, IA1 = 40014, IA2 = 40692, IQ1 = 53668, IQ2 =      &
      52774, IR1 = 12211, IR2 = 3791, NTAB = 32, NDIJV = 1 + IMM1 /     &
      NTAB, EPS = 1.2D-10, RNMX = 1.0D0 - EPS) 
      DIMENSION IJV (NTAB)
      SAVE IJV, IJY, IJDUM2
      DATA IJDUM2 / 123456789 /, IJV / NTAB * 0 /, IJY / 0 /
!                                                                       
      IF (IDUM.LE.0) THEN
         IDUM = MAX ( - IDUM, 1)
         IJDUM2 = IDUM
         DO 11 J = NTAB + 8, 1, - 1
            K = IDUM / IQ1
            IDUM = IA1 * (IDUM - K * IQ1) - K * IR1
            IF (IDUM.LT.0) IDUM = IDUM + IM1
            IF (J.LE.NTAB) IJV (J) = IDUM
   11    END DO
         IJY = IJV (1)
      ENDIF
      K = IDUM / IQ1
      IDUM = IA1 * (IDUM - K * IQ1) - K * IR1
      IF (IDUM.LT.0) IDUM = IDUM + IM1
      K = IJDUM2 / IQ2
      IJDUM2 = IA2 * (IJDUM2 - K * IQ2) - K * IR2
      IF (IJDUM2.LT.0) IJDUM2 = IJDUM2 + IM2
      J = 1 + IJY / NDIJV
      IJY = IJV (J) - IJDUM2
      IJV (J) = IDUM
      IF (IJY.LT.1) IJY = IJY + IMM1
      RAN32 = MIN (AM * IJY, RNMX)
      RETURN
      END FUNCTION RAN32

!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************
!********************************************

end module
