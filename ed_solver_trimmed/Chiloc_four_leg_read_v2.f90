!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!

module cdagger

 INTERFACE ASSIGNMENT (=)
    MODULE PROCEDURE cdagger_copy
  END INTERFACE

  TYPE eigen
   real(8),allocatable :: val(:)
  END TYPE

  TYPE cdagger_mat
   integer    :: k_p,l_p
   integer    :: k_m,l_m
   real(8),   allocatable :: c_p(:,:) , c_m(:,:)
  END TYPE

contains

         !-------------------------!

  subroutine cdagger_copy(cout,cin)
  implicit none
  Type(cdagger_mat),intent(inout) :: cout
  Type(cdagger_mat),intent(in)    :: cin
   call kill_cdagger(cout)
   cout%k_p=cin%k_p
   cout%l_p=cin%l_p
   cout%k_m=cin%k_m
   cout%l_m=cin%l_m
   call allocate_dagger(cout)
   cout%c_p=cin%c_p
   cout%c_m=cin%c_m
  end subroutine

         !-------------------------!

  subroutine kill_cdagger(cdagger)
  implicit none
  Type(cdagger_mat) :: cdagger
      if(allocated(cdagger%c_p)) deallocate(cdagger%c_p)
      if(allocated(cdagger%c_m)) deallocate(cdagger%c_m)
  end subroutine

         !-------------------------!

  subroutine allocate_dagger(cdagger)
  implicit none
  Type(cdagger_mat) :: cdagger
   call allocate_dagger_p(cdagger)
   call allocate_dagger_m(cdagger)
  end subroutine

         !-------------------------!

  subroutine allocate_dagger_p(cdagger)
  implicit none
  Type(cdagger_mat) :: cdagger
    if(allocated(cdagger%c_p)) deallocate(cdagger%c_p)
    allocate(cdagger%c_p(cdagger%k_p,cdagger%l_p))
  end subroutine

         !-------------------------!

  subroutine allocate_dagger_m(cdagger)
  implicit none
  Type(cdagger_mat) :: cdagger
    if(allocated(cdagger%c_m)) deallocate(cdagger%c_m)
    allocate(cdagger%c_m(cdagger%k_m,cdagger%l_m))
  end subroutine

         !-------------------------!

end module

!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!

program read_chi_loc
use cdagger
use mpirout  ! CW_MPI : link to the mpi module
implicit none

interface
 subroutine  chi_tilde_loc( &
 & op, w1,  w2,  w3, PHI_EPS, beta, Z, gs_E, nsites, nup, ndn,  &
 &    cp_i_E,              dim_E_i,      & 
 &    cp_pup_E,          dim_E_pup,    & 
 &    cp_pdn_E,          dim_E_pdn,    & 
 &    cp_mup_E,          dim_E_mup,    & 
 &    cp_mdn_E,          dim_E_mdn,    & 
 &    cp_p2dn_E,         dim_E_p2dn,   & 
 &    cp_m2dn_E,         dim_E_m2dn,   & 
 &    cp_puppdn_E,    dim_E_puppdn, & 
 &    cp_muppdn_E,    dim_E_muppdn, & 
 &    cp_pupmdn_E,    dim_E_pupmdn, & 
 &    cp_mupmdn_E,    dim_E_mupmdn, & 
 &    n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,&
 &    cp_i_cdup,& 
 &    cp_i_cddn,& 
 &    cp_pup_cddn,& 
 &    cp_pdn_cdup,& 
 &    cp_pdn_cddn,& 
 &    cp_mup_cdup,& 
 &    cp_mup_cddn,& 
 &    cp_mdn_cddn,& 
 &    cp_mdn_cdup,& 
 &    cp_muppdn_cdup,& 
 &    cp_pupmdn_cddn,& 
 &    cp_mupmdn_cdup,& 
 &    cp_mupmdn_cddn,& 
 &    cp_m2dn_cddn, & 
 &    pDNDN,pUPDN )
implicit none
integer :: op,nsites,nup,ndn
complex(8) :: w1,w2,w3
real(8) :: PHI_EPS,beta,Z,gs_E
complex(8) :: pDNDN,pUPDN
integer :: dim_E_i      
integer :: dim_E_pup    
integer :: dim_E_pdn    
integer :: dim_E_mup    
integer :: dim_E_mdn    
integer :: dim_E_p2dn   
integer :: dim_E_m2dn   
integer :: dim_E_puppdn 
integer :: dim_E_muppdn   
integer :: dim_E_pupmdn   
integer :: dim_E_mupmdn   

integer :: n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28
real(8) :: cp_i_E(n2)
real(8) :: cp_pup_E(n1)
real(8) :: cp_pdn_E(n3)
real(8) :: cp_mup_E(n12)
real(8) :: cp_mdn_E(n16)
real(8) :: cp_p2dn_E(n9)
real(8) :: cp_m2dn_E(n28)
real(8) :: cp_puppdn_E(n5)
real(8) :: cp_muppdn_E(n13)
real(8) :: cp_pupmdn_E(n17)
real(8) :: cp_mupmdn_E(n24)

real(8) :: cp_i_cdup(n2,n1) !cp_i_cdup(n1,n2) The matrix passed to c++ routine has to be transposed.       
real(8) :: cp_i_cddn(n4,n3) !cp_i_cddn(n3,n4)       
real(8) :: cp_pup_cddn(n6,n5) !cp_pup_cddn(n5,n6)     
real(8) :: cp_pdn_cdup(n8,n7) !cp_pdn_cdup(n7,n8)
real(8) :: cp_pdn_cddn(n10,n9) !cp_pdn_cddn(n9,n10)     
real(8) :: cp_mup_cdup(n12,n11) !cp_mup_cdup(n11,n12)     
real(8) :: cp_mup_cddn(n14,n13) !cp_mup_cddn(n13,n14)     
real(8) :: cp_mdn_cddn(n16,n15) !cp_mdn_cddn(n15,n16)     
real(8) :: cp_mdn_cdup(n18,n17) !cp_mdn_cdup(n17,n18)     
real(8) :: cp_muppdn_cdup(n20,n19) !cp_muppdn_cdup(n19,n20)  
real(8) :: cp_pupmdn_cddn(n22,n21) !cp_pupmdn_cddn(n21,n22)  
real(8) :: cp_mupmdn_cdup(n24,n23) !cp_mupmdn_cdup(n23,n24)  
real(8) :: cp_mupmdn_cddn(n26,n25) !cp_mupmdn_cddn(n25,n26)  
real(8) :: cp_m2dn_cddn(n28,n27) !cp_m2dn_cddn(n27,n28)   
 end subroutine
end interface

!CW_MPI
!complex(8)                          :: Chi_lattice
!END CW MPI
complex(8)                          :: w1,w2,w3
integer                             :: min_up,max_up,min_dn,max_dn,nsector
integer,allocatable                 :: nup(:),ndn(:)
integer                             :: neigen,k,l
TYPE(eigen),allocatable             :: eigenval(:,:)
TYPE(cdagger_mat),allocatable       :: cup_mat(:,:),cdn_mat(:,:)
integer                             :: i,j
real(8)                             :: temp
real(8)                             :: PHI_EPS,beta,ZZ
integer                             :: nsites,op

integer                             :: dim_E_i
integer                             :: dim_E_pup
integer                             :: dim_E_pdn
integer                             :: dim_E_mup
integer                             :: dim_E_mdn
integer                             :: dim_E_p2dn
integer                             :: dim_E_m2dn
integer                             :: dim_E_puppdn
integer                             :: dim_E_muppdn
integer                             :: dim_E_pupmdn
integer                             :: dim_E_mupmdn
real(8),allocatable                 :: cp_i_E(:)
real(8),allocatable                 :: cp_pup_E(:)
real(8),allocatable                 :: cp_pdn_E(:)
real(8),allocatable                 :: cp_mup_E(:)
real(8),allocatable                 :: cp_mdn_E(:)
real(8),allocatable                 :: cp_p2dn_E(:)
real(8),allocatable                 :: cp_m2dn_E(:)
real(8),allocatable                 :: cp_puppdn_E(:)
real(8),allocatable                 :: cp_muppdn_E(:)
real(8),allocatable                 :: cp_pupmdn_E(:)
real(8),allocatable                 :: cp_mupmdn_E(:)
integer                             :: n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28
real(8),allocatable                 :: cp_i_cdup(:,:)
real(8),allocatable                 :: cp_i_cddn(:,:)
real(8),allocatable                 :: cp_pup_cddn(:,:)
real(8),allocatable                 :: cp_pdn_cdup(:,:)
real(8),allocatable                 :: cp_pdn_cddn(:,:)
real(8),allocatable                 :: cp_mup_cdup(:,:)
real(8),allocatable                 :: cp_mup_cddn(:,:)
real(8),allocatable                 :: cp_mdn_cddn(:,:)
real(8),allocatable                 :: cp_mdn_cdup(:,:)
real(8),allocatable                 :: cp_muppdn_cdup(:,:)
real(8),allocatable                 :: cp_pupmdn_cddn(:,:)
real(8),allocatable                 :: cp_mupmdn_cdup(:,:)
real(8),allocatable                 :: cp_mupmdn_cddn(:,:)
real(8),allocatable                 :: cp_m2dn_cddn(:,:)
integer,allocatable                 :: neigentab(:,:)
complex(8)                          :: pDNDN,pUPDN

real(8)                             :: gs_E
real(8),allocatable                 :: minE_list(:)
real(8)                             :: Z_test

!CW_MPI
integer                             :: i_,j_
complex(8),allocatable              :: Chi_lattice(:)
complex(8),allocatable              :: frequ_(:,:) 
integer,allocatable                 :: op_(:) 
!END CW_MPI


     call initialize_MPI !CW_MPI : init mpi

     PHI_EPS=1.d-9

     open(unit=1414,file='chilattice_vertex',form='unformatted')
     open(unit=1415,form='formatted')
     read(1414) ZZ,beta,nsites,min_up,max_up,min_dn,max_dn,nsector

     write(1415,*) 'partition function : ', ZZ
     write(1415,*) 'beta               : ', beta
     write(1415,*) 'nsites             : ', nsites
     write(1415,*) 'nsector            : ', nsector

     allocate(nup(nsector),ndn(nsector))
     if(allocated(cup_mat)) deallocate(cup_mat)
     if(allocated(cdn_mat)) deallocate(cdn_mat)
     allocate(cup_mat(min_up:max_up,min_dn:max_dn))
     allocate(cdn_mat(min_up:max_up,min_dn:max_dn))
     allocate(eigenval(min_up:max_up,min_dn:max_dn))
     allocate(neigentab(min_up:max_up,min_dn:max_dn))

     write(1415,*) 'start reading' 
     allocate(minE_list(nsector))
     do i=1,nsector
       write(1415,*) '!########################################'
       write(1415,*) 'reading sector / nsector: ', i,nsector
       read(1414) nup(i),ndn(i)
       read(1414) neigen
       neigentab(nup(i),ndn(i))=neigen
       allocate(eigenval(nup(i),ndn(i))%val(neigen))
       read(1414) eigenval(nup(i),ndn(i))%val
       minE_list(i)=minval(eigenval(nup(i),ndn(i))%val)

       write(1415,*) 'nup,ndn   : ',nup(i),ndn(i) 
       write(1415,'(a,200f9.3)') 'eigenvals : ',eigenval(nup(i),ndn(i))%val
       write(1415,'(a,f9.3)') 'the lowest energy is ', minE_list(i)
       read(1414) k,l
       cup_mat(nup(i),ndn(i))%k_p=k
       cup_mat(nup(i),ndn(i))%l_p=l
       call allocate_dagger_p(cup_mat(nup(i),ndn(i)))
       read(1414) cup_mat(nup(i),ndn(i))%c_p

       read(1414) k,l
       cup_mat(nup(i),ndn(i))%k_m=k
       cup_mat(nup(i),ndn(i))%l_m=l
       call allocate_dagger_m(cup_mat(nup(i),ndn(i)))
       read(1414) cup_mat(nup(i),ndn(i))%c_m

       read(1414) k,l
       cdn_mat(nup(i),ndn(i))%k_p=k
       cdn_mat(nup(i),ndn(i))%l_p=l
       call allocate_dagger_p(cdn_mat(nup(i),ndn(i)))
       read(1414) cdn_mat(nup(i),ndn(i))%c_p

       read(1414) k,l
       cdn_mat(nup(i),ndn(i))%k_m=k
       cdn_mat(nup(i),ndn(i))%l_m=l
       call allocate_dagger_m(cdn_mat(nup(i),ndn(i)))
       read(1414) cdn_mat(nup(i),ndn(i))%c_m

     enddo
     close(1414)
     close(1415)

!-------------------------------------------------------------------------!
! Ground state energy.                     
     gs_E = minval(minE_list)
     write(*,*) 'Ground state energy=', gs_E
! Partition function. 
!     Z_test = 0.0
!     do i=1,nsector
!         neigen=neigentab(nup(i),ndn(i))
!         do j=1,neigen
!            Z_test=Z_test+exp(-beta*(eigenval(nup(i),ndn(i))%val(j)-gs_E))
!         enddo
!     enddo
!     write(*,*) 'partition function=', Z_test
!-------------------------------------------------------------------------!


    !---------------------------!
    !... for a fix frequency ...!
    !---------------------------! 

!     w1=(0.,0.)
!     w2=(0.,0.)
!     w3=(0.,0.)

!CW_MPI check how many frequ in omega_list
     open(unit=4848, file="omega_list", form="formatted")
     j_=0
     do 
      read(4848, *,end=871) op,w1,w2,w3
      j_=j_+1
      write(*,*) "in file op         = ", op
      write(*,*) "in file w1, w2, w3 = ", w1, w2, w3
     enddo
   871 continue

     write(*,*) 'THERE ARE [x] FREQUENCIES IN FILE : ', j_
     write(*,*) 'allocating Chilattice array'
     if(allocated(Chi_lattice)) deallocate(Chi_lattice)
     allocate(Chi_lattice(j_),frequ_(j_,3),op_(j_))
     Chi_lattice = 0.d0
     op_         = 0
     frequ_      = 0.d0
     rewind(4848)

     do i_=1,j_
       read(4848, *) op,w1,w2,w3
       write(*,*) "op         = ", op
       write(*,*) "w1, w2, w3 = ", w1, w2, w3
       op_(i_)      = op
       frequ_(i_,1) = w1
       frequ_(i_,2) = w2
       frequ_(i_,3) = w3
     enddo
     close(4848)

   do i_=rank+1,j_,size2 !CW_MPI : this is the parallelisation of the loop over frequencies
 
     op=op_(i_)
     w1=frequ_(i_,1)
     w2=frequ_(i_,2)
     w3=frequ_(i_,3)

     write(*,*) 'MPI RANK = ',rank
     write(*,*) 'MPI TOTAL NUMBER OF CPUS = ', size2
     write(*,*) 'I AM DOING FREQUENCY NUMBER = ', i_
!END CW_MPI

!========================================================================!
!========================================================================!
!========================================================================!
!========================================================================!
!========================================================================!

     open(unit=5486, form="formatted")
     write(5486,*) "! Print out for each E_i sector eigenvalue list"
     write(5486,*) "! and matrix elements."

     do i=1,nsector

!We know from above:
! nup(i),ndn(i)                           nup,ndn of this sector
! neigentab(nup(i),ndn(i))                the number of eigenstates of this sector
! eigenval(nup(i),ndn(i))%val(:)          the list of eigenvalues in sector nup(i),ndn(i)  
! cup_mat(nup(i),ndn(i))%c_p              the matrix <j | c_up^+ | i> 
! cup_mat(nup(i),ndn(i))%k_p and %l_p     the dimensions (k,l) of the matrix above, note that _p stands for plus (c^+) and _m for minus (c)
! cup_mat(nup(i),ndn(i))%c_m              the matrix <j | c_up | i> 
! cdn_mat(nup(i),ndn(i))%c_(p,m)          same thing for c_dn and c_dn^+

! We need to copy this values to the arrays below, by playing with the indices nup(i)+1 ....
! Some things to keep in mind :

! a) if nup(i)=max_up, nup(i)+1 is not defined (max_up is the largest nup sector)

! b) here we are defining 2D arrays for cp_i_cdup(:,:), BUT eventualy the routine of Andreas will receive it as a 1D vector
!    (2D arrays actually don't really exist, they are Fortran construction. 
!    the array received by Andreas's routines is the 2D array defined in Fortran where all the columns are sticked together, starting
!    with the first column, then comes the second column etcetera. 
!    Note : We should check in Andreas code how he defines his matrices from a vector
!           Andreas is using a notation by "line" I think, instead of columns
!           so we need actually to defined :  cp_i_cdup = transpose( cup_mat(nup(i),ndn(i))%c_p )
!           and defined the dimensions of the array cp_i_cdup accordingly with the couple (n1,n2) 

!#############################################
! assignment of number of eigenvalues ########
! of sectors needed ##########################
!#############################################
! dim_E_i
    dim_E_i= neigentab(nup(i), ndn(i))
! dim_E_pup
! if it enters a sector which is beyond the problem, give it a dim_E = -1, such
! that in chi_tilde_loc(), this dimension will not be iterated. 
    if (nup(i)+1 > max_up) then 
        dim_E_pup= -1 
    else 
        dim_E_pup= neigentab(nup(i)+1, ndn(i))
    end if
! dim_E_pdn
    if (ndn(i)+1 > max_dn) then
        dim_E_pdn= -1
    else 
        dim_E_pdn= neigentab(nup(i), ndn(i)+1)
    end if
! dim_E_mup
    if (nup(i)-1 < min_up) then
        dim_E_mup= -1
    else
        dim_E_mup= neigentab(nup(i)-1, ndn(i))
    end if
! dim_E_mdn
    if (ndn(i)-1 < min_dn) then
        dim_E_mdn= -1
    else
        dim_E_mdn= neigentab(nup(i), ndn(i)-1)
    end if
! dim_E_p2dn
    if (ndn(i)+2 > max_dn) then
        dim_E_p2dn= -1
    else 
        dim_E_p2dn= neigentab(nup(i), ndn(i)+2)
    end if
! dim_E_m2dn
    if (ndn(i)-2 < min_dn) then
        dim_E_m2dn= -1
    else 
        dim_E_m2dn= neigentab(nup(i), ndn(i)-2)
    end if
! dim_E_puppdn
    if (nup(i)+1 > max_up .or. ndn(i)+1 > max_dn) then
        dim_E_puppdn= -1
    else 
        dim_E_puppdn= neigentab(nup(i)+1, ndn(i)+1)
    end if
! dim_E_muppdn
    if (nup(i)-1 < min_up .or. ndn(i)+1 > max_dn) then
        dim_E_muppdn= -1
    else
        dim_E_muppdn= neigentab(nup(i)-1, ndn(i)+1)
    end if
! dim_E_pupmdn
    if (nup(i)+1 > max_up .or. ndn(i)-1 < min_dn) then
        dim_E_pupmdn= -1
    else
        dim_E_pupmdn= neigentab(nup(i)+1, ndn(i)-1)
    end if
! dim_E_mupmdn
    if (nup(i)-1 < min_up .or. ndn(i)-1 < min_dn) then
        dim_E_mupmdn= -1
    else
        dim_E_mupmdn= neigentab(nup(i)-1, ndn(i)-1)
    end if

!###############################################
! assign n1,n2,... for allocating eigen ########
! value arrays and c^dagger matrices ###########
!###############################################
! dim_E_i
    n2=dim_E_i; n4=dim_E_i; n11=dim_E_i; n15=dim_E_i
! dim_E_pup
    if (dim_E_pup == -1) then
        n1=1; n6=1; n21=1   ! if nup(i)+1>max_up, then give n1=1, such that the 
                            ! array/matrix can be initialized. but since dim_E_pup=-1, the iteration in 
    else                    ! chi_tilde_loc() will not be executed, and thus the array/matrix will not be used anyway. 
        n1=dim_E_pup; n6=dim_E_pup; n21=dim_E_pup
    end if
! dim_E_pdn
    if (dim_E_pdn == -1) then
        n3=1; n8=1; n10=1; n19=1
    else
        n3=dim_E_pdn; n8=dim_E_pdn; n10=dim_E_pdn; n19=dim_E_pdn
    end if
! dim_E_mup
    if (dim_E_mup == -1) then
        n12=1; n14=1; n25=1
    else 
        n12=dim_E_mup; n14=dim_E_mup; n25=dim_E_mup
    end if
! dim_E_mdn
    if (dim_E_mdn == -1) then
        n16=1; n18=1; n23=1; n27=1
    else
        n16=dim_E_mdn; n18=dim_E_mdn; n23=dim_E_mdn; n27=dim_E_mdn
    end if
! dim_E_p2dn
    if (dim_E_p2dn == -1) then
        n9=1
    else 
        n9=dim_E_p2dn
    end if
! dim_E_m2dn
    if (dim_E_m2dn == -1) then
        n28=1
    else 
        n28=dim_E_m2dn
    end if
! dim_E_puppdn
    if (dim_E_puppdn == -1) then
        n5=1; n7=1
    else 
        n5=dim_E_puppdn; n7=dim_E_puppdn
    end if
! dim_E_muppdn
    if (dim_E_muppdn == -1) then
        n13=1; n20=1
    else 
        n13=dim_E_muppdn; n20=dim_E_muppdn
    end if
! dim_E_pupmdn
    if (dim_E_pupmdn == -1) then
        n17=1; n22=1
    else
        n17=dim_E_pupmdn; n22=dim_E_pupmdn
    end if
! dim_E_mupmdn
    if (dim_E_mupmdn == -1) then
        n24=1; n26=1
    else 
        n24=dim_E_mupmdn; n26=dim_E_mupmdn
    end if


     if(allocated(cp_i_E)) deallocate(cp_i_E, cp_pup_E, cp_pdn_E, cp_mup_E, cp_mdn_E, cp_p2dn_E, cp_m2dn_E, cp_puppdn_E, &
& cp_muppdn_E, cp_pupmdn_E, cp_mupmdn_E, cp_i_cdup,cp_i_cddn,                        & 
& cp_pup_cddn, cp_pdn_cdup, cp_pdn_cddn, cp_mup_cdup, cp_mup_cddn, cp_mdn_cddn,      & 
& cp_mdn_cdup, cp_muppdn_cdup, cp_pupmdn_cddn, cp_mupmdn_cdup, cp_mupmdn_cddn, cp_m2dn_cddn)

!################################################
! allocate arrays of eigenvalues ################
! ###############################################
     allocate( cp_i_E(n2) )
     allocate( cp_pup_E(n1) )
     allocate( cp_pdn_E(n3) )
     allocate( cp_mup_E(n12) )
     allocate( cp_mdn_E(n16) )
     allocate( cp_p2dn_E(n9) )
     allocate( cp_m2dn_E(n28) )
     allocate( cp_puppdn_E(n5) )
     allocate( cp_muppdn_E(n13) )
     allocate( cp_pupmdn_E(n17) )
     allocate( cp_mupmdn_E(n24) )

!#######################################
! allocate of c^dagger matrices ########
!#######################################
     allocate( cp_i_cdup(n2,n1) )
     allocate( cp_i_cddn(n4,n3) )
     allocate( cp_pup_cddn(n6,n5) )
     allocate( cp_pdn_cdup(n8,n7) )
     allocate( cp_pdn_cddn(n10,n9) )
     allocate( cp_mup_cdup(n12,n11) )
     allocate( cp_mup_cddn(n14,n13) )
     allocate( cp_mdn_cddn(n16,n15) )
     allocate( cp_mdn_cdup(n18,n17) )
     allocate( cp_muppdn_cdup(n20,n19) )
     allocate( cp_pupmdn_cddn(n22,n21) )
     allocate( cp_mupmdn_cdup(n24,n23) )
     allocate( cp_mupmdn_cddn(n26,n25) )
     allocate( cp_m2dn_cddn(n28,n27) )

!########################################
! assignment of lists of eigenvalues ####
!########################################
! dim_E_i
     cp_i_E = eigenval(nup(i), ndn(i))%val
! dim_E_pup
     if (dim_E_pup == -1) then
         cp_pup_E(1:n1) = 0.0
     else 
         cp_pup_E = eigenval(nup(i)+1, ndn(i))%val
     end if
! dim_E_pdn
     if (dim_E_pdn == -1) then
         cp_pdn_E(1:n3) = 0.0
     else 
         cp_pdn_E = eigenval(nup(i), ndn(i)+1)%val
     end if
! dim_E_mup
     if (dim_E_mup == -1) then
         cp_mup_E(1:n12) =0.0
     else 
         cp_mup_E = eigenval(nup(i)-1, ndn(i))%val
     end if
! dim_E_mdn
     if (dim_E_mdn == -1) then
         cp_mdn_E(1:n16) = 0.0
     else 
         cp_mdn_E = eigenval(nup(i), ndn(i)-1)%val
     end if
! dim_E_p2dn
     if (dim_E_p2dn == -1) then
         cp_p2dn_E(1:n9) = 0.0
     else 
         cp_p2dn_E = eigenval(nup(i), ndn(i)+2)%val
     end if
! dim_E_m2dn
     if (dim_E_m2dn == -1) then
         cp_m2dn_E(1:n28) = 0.0
     else 
         cp_m2dn_E = eigenval(nup(i), ndn(i)-2)%val
     end if
! dim_E_puppdn
     if (dim_E_puppdn == -1) then
         cp_puppdn_E(1:n5) = 0.0
     else 
         cp_puppdn_E = eigenval(nup(i)+1, ndn(i)+1)%val
     end if
! dim_E_muppdn
     if (dim_E_muppdn == -1) then
         cp_muppdn_E(1:n13) = 0.0
     else 
         cp_muppdn_E = eigenval(nup(i)-1, ndn(i)+1)%val
     end if
! dim_E_pupmdn
     if (dim_E_pupmdn == -1) then
         cp_pupmdn_E(1:n17) = 0.0
     else 
         cp_pupmdn_E = eigenval(nup(i)+1, ndn(i)-1)%val
     end if
! dim_E_mupmdn
     if (dim_E_mupmdn == -1) then
         cp_mupmdn_E(1:n24) = 0.0
     else 
         cp_mupmdn_E = eigenval(nup(i)-1, ndn(i)-1)%val
     end if 

!########################################
! assignment of c^dagger matrices #######
!########################################
! cp_i_cdup
     if (dim_E_pup == -1) then
         cp_i_cdup(1:n2, 1:n1) = 0.0
     else
         cp_i_cdup(1:n2, 1:n1) = transpose(cup_mat(nup(i),ndn(i))%c_p)
     end if
! cp_i_cddn
     if (dim_E_pdn == -1) then
         cp_i_cddn(1:n4, 1:n3) = 0.0 
     else 
         cp_i_cddn(1:n4, 1:n3) = transpose(cdn_mat(nup(i),ndn(i))%c_p)
     end if 
! cp_pup_cddn
     if (dim_E_puppdn == -1 .or. dim_E_pup == -1) then
         cp_pup_cddn(1:n6, 1:n5) = 0.0
     else
         cp_pup_cddn(1:n6, 1:n5) = transpose(cdn_mat(nup(i)+1,ndn(i))%c_p)
     end if
! cp_pdn_cdup
     if (dim_E_puppdn==-1 .or. dim_E_pdn==-1) then
         cp_pdn_cdup(1:n8, 1:n7) = 0.0
     else 
         cp_pdn_cdup(1:n8, 1:n7) = transpose(cup_mat(nup(i),ndn(i)+1)%c_p)
     end if
! cp_pdn_cddn
     if (dim_E_p2dn==-1 .or. dim_E_pdn==-1) then
         cp_pdn_cddn(1:n10, 1:n9) = 0.0
     else 
         cp_pdn_cddn(1:n10, 1:n9) = transpose(cdn_mat(nup(i),ndn(i)+1)%c_p)
     end if
! cp_mup_cdup
     if (dim_E_mup==-1 ) then
         cp_mup_cdup(1:n12, 1:n11) = 0.0
     else
         cp_mup_cdup(1:n12, 1:n11) = transpose(cup_mat(nup(i)-1,ndn(i))%c_p)
     end if
! cp_mup_cddn
     if (dim_E_muppdn==-1 .or. dim_E_mup==-1) then
         cp_mup_cddn(1:n14, 1:n13) = 0.0
     else
         cp_mup_cddn(1:n14, 1:n13) = transpose(cdn_mat(nup(i)-1,ndn(i))%c_p)
     end if
! cp_mdn_cddn
     if (dim_E_mdn==-1) then
         cp_mdn_cddn(1:n16, 1:n15) = 0.0
     else
         cp_mdn_cddn(1:n16, 1:n15) = transpose(cdn_mat(nup(i),ndn(i)-1)%c_p)
     end if
! cp_mdn_cdup
     if (dim_E_pupmdn==-1 .or. dim_E_mdn==-1) then
         cp_mdn_cdup(1:n18, 1:n17) = 0.0
     else
         cp_mdn_cdup(1:n18, 1:n17) = transpose(cup_mat(nup(i),ndn(i)-1)%c_p)
     end if
! cp_muppdn_cdup
     if (dim_E_pdn==-1 .or. dim_E_muppdn==-1) then
         cp_muppdn_cdup(1:n20, 1:n19) = 0.0
     else
         cp_muppdn_cdup(1:n20, 1:n19) = transpose(cup_mat(nup(i)-1,ndn(i)+1)%c_p)
     end if
! cp_pupmdn_cddn
     if (dim_E_pup==-1 .or. dim_E_pupmdn==-1) then
         cp_pupmdn_cddn(1:n22, 1:n21) = 0.0
     else
         cp_pupmdn_cddn(1:n22, 1:n21) = transpose(cdn_mat(nup(i)+1,ndn(i)-1)%c_p)
     end if
! cp_mupmdn_cdup
     if (dim_E_mdn==-1 .or. dim_E_mupmdn==-1) then
         cp_mupmdn_cdup(1:n24, 1:n23) = 0.0
     else
         cp_mupmdn_cdup(1:n24, 1:n23) = transpose(cup_mat(nup(i)-1,ndn(i)-1)%c_p)
     end if
! cp_mupmdn_cddn
     if (dim_E_mup==-1 .or. dim_E_mupmdn==-1) then
         cp_mupmdn_cddn(1:n26, 1:n25) = 0.0
     else
         cp_mupmdn_cddn(1:n26, 1:n25) = transpose(cdn_mat(nup(i)-1,ndn(i)-1)%c_p)
     end if
! cp_m2dn_cddn
     if (dim_E_mdn==-1 .or. dim_E_m2dn==-1) then
         cp_m2dn_cddn(1:n28, 1:n27) = 0.0
     else
         cp_m2dn_cddn(1:n28, 1:n27) = transpose(cdn_mat(nup(i),ndn(i)-2)%c_p)
     end if

!############################################
! call chi_tilde_loc() ######################
!############################################

! ######### check scalar passing. ###################
!      write(*,*) 'w1, w2, w3 in fortran =', w1,w2,w3
!      write(*,*) 'PHI_EPS, beta, Z in fortran = ', PHI_EPS, beta, ZZ
!      write(*,*) 'sites in fortran = ', nsites
!      write(*,*) 'op in fortran  = ', op
!      write(*,*) "w1, w2, w3 in fortran = ", w1, w2, w3
!      write(*,*) "dim_E's in fortran = ", dim_E_i, dim_E_pup, dim_E_pdn, dim_E_mup, dim_E_mdn, dim_E_p2dn, dim_E_m2dn, &
!                                          dim_E_puppdn, dim_E_muppdn, dim_E_pupmdn, dim_E_mupmdn
!     write(5486,*) "########################################\n"
!     write(5486,"(3(a8, i3))") 'nsites=', nsites, 'nup=', nup(i), 'ndn=', ndn(i)

!     write(5486,*) "cp_i_E = "
!     write(5486,"(7(f8.4, tr1))"), cp_i_E(1:n2)
!
!     write(5486,*) "cp_pup_E = "
!     write(5486,"(7(f8.4, tr1))"), cp_pup_E(1:n1)
!
!     write(5486,*) "cp_pdn_E = "
!     write(5486,"(7(f8.4, tr1))"), cp_pdn_E(1:n3)
!
!     write(5486,*) "cp_mup_E = "
!     write(5486,"(7(f8.4, tr1))"), cp_mup_E(1:n12)
!
!     write(5486,*) "cp_mdn_E = "
!     write(5486,"(7(f8.4, tr1))"), cp_mdn_E(1:n16)
!
!     write(5486,*) "cp_p2dn_E = "
!     write(5486,"(7(f8.4, tr1))"), cp_p2dn_E(1:n9)
!
!     write(5486,*) "cp_m2dn_E = "
!     write(5486,"(7(f8.4, tr1))"), cp_m2dn_E(1:n28)
!
!     write(5486,*) "cp_puppdn_E = "
!     write(5486,"(7(f8.4, tr1))"), cp_puppdn_E(1:n5)
!
!     write(5486,*) "cp_muppdn_E = "
!     write(5486,"(7(f8.4, tr1))"), cp_muppdn_E(1:n13)
!
!     write(5486,*) "cp_pupmdn_E = "
!     write(5486,"(7(f8.4, tr1))"), cp_pupmdn_E(1:n17)
! 
!     write(5486,*) "cp_mupmdn_E = "
!     write(5486,"(7(f8.4, tr1))"), cp_mupmdn_E(1:n24)

!      write(5486, *) "cp_i_cdup = "
!      write(5486, "(7(f8.4, tr1))"), cp_i_cdup(1:n2,1:n1)
! 
!      write(5486, *) "cp_i_cddn = "
!      write(5486, "(7(f8.4, tr1))"), cp_i_cddn(1:n4,1:n3)
!
!      write(5486, *) "cp_pup_cddn = "
!      write(5486, "(7(f8.4, tr1))"), cp_pup_cddn(1:n6,1:n5)
!
!      write(5486, *) "cp_pdn_cdup = "
!      write(5486, "(7(f8.4, tr1))"), cp_pdn_cdup(1:n8,1:n7)
!
!      write(5486, *) "cp_pdn_cddn = "
!      write(5486, "(7(f8.4, tr1))"), cp_pdn_cddn(1:n10,1:n9)
!
!      write(5486, *) "cp_mup_cdup = "
!      write(5486, "(7(f8.4, tr1))"), cp_mup_cdup(1:n12,1:n11)
!
!      write(5486, *) "cp_mup_cddn = "
!      write(5486, "(7(f8.4, tr1))"), cp_mup_cddn(1:n14,1:n13)
!
!      write(5486, *) "cp_mdn_cddn = "
!      write(5486, "(7(f8.4, tr1))"), cp_mdn_cddn(1:n16,1:n15)
!      
!      write(5486, *) "cp_mdn_cdup = "
!      write(5486, "(7(f8.4, tr1))"), cp_mdn_cdup(1:n18,1:n17)
!
!      write(5486, *) "cp_muppdn_cdup = "
!      write(5486, "(7(f8.4, tr1))"), cp_muppdn_cdup(1:n20,1:n19)
!
!      write(5486, *) "cp_pupmdn_cddn = "
!      write(5486, "(7(f8.4, tr1))"), cp_pupmdn_cddn(1:n22,1:n21)
!
!      write(5486, *) "cp_mupmdn_cdup = "
!      write(5486, "(7(f8.4, tr1))"), cp_mupmdn_cdup(1:n24,1:n23)
!
!      write(5486, *) "cp_mupmdn_cddn = "
!      write(5486, "(7(f8.4, tr1))"), cp_mupmdn_cddn(1:n26,1:n25)
!
!      write(5486, *) "cp_m2dn_cddn = "
!      write(5486, "(7(f8.4, tr1))"), cp_m2dn_cddn(1:n28,1:n27)

     call chi_tilde_loc( &
 &     op, w1,  w2,  w3, PHI_EPS, beta, ZZ, gs_E, nsites, nup(i), ndn(i), & 
 &    cp_i_E,            dim_E_i,      &
 &    cp_pup_E,          dim_E_pup,    &
 &    cp_pdn_E,          dim_E_pdn,    &
 &    cp_mup_E,          dim_E_mup,    &
 &    cp_mdn_E,          dim_E_mdn,    &
 &    cp_p2dn_E,         dim_E_p2dn,   &
 &    cp_m2dn_E,         dim_E_m2dn,   &
 &    cp_puppdn_E,       dim_E_puppdn, &
 &    cp_muppdn_E,       dim_E_muppdn, &
 &    cp_pupmdn_E,       dim_E_pupmdn, &
 &    cp_mupmdn_E,       dim_E_mupmdn, &
 &    n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,&
 &    cp_i_cdup,       &
 &    cp_i_cddn,       &
 &    cp_pup_cddn,     &
 &    cp_pdn_cdup,     &
 &    cp_pdn_cddn,     &
 &    cp_mup_cdup,     &
 &    cp_mup_cddn,     &
 &    cp_mdn_cddn,     &
 &    cp_mdn_cdup,     &
 &    cp_muppdn_cdup,  &
 &    cp_pupmdn_cddn,  &
 &    cp_mupmdn_cdup,  &
 &    cp_mupmdn_cddn,  &
 &    cp_m2dn_cddn,    &
 &    pDNDN,pUPDN )
     
!      write(*,*) 'pDNDN, pUPDN = ', pDNDN, pUPDN
      Chi_lattice(i_) = Chi_lattice(i_) + pDNDN + pUPDN
!
     enddo
!     close(5486)

!========================================================================!
!========================================================================!
!========================================================================!
!========================================================================!
!========================================================================!

!CW_MPI 
    enddo

    call mpisum(Chi_lattice) !collect array Chi_lattice

     write(*,*) 'Chi is : ' , Chi_lattice

     if(rank==0)then
      !only rank=0 (master node) writes into file, otherwise it is a big mess
      ! all the nodes would write into the same file in random order
      write(*,*) 'Dumping Chi_lattice in file Chi_output'
      write(*,*) 'each line is the Chi_lat of the frequency contained in frequency input file'
      open(unit=8080,file='Chi_output')
      do i_=1,j_
       write(8080,*) Chi_lattice(i_)
      enddo
      close(8080)
     endif

     call finalize_MPI !CW_MPI : shut down mpi
!END CW_MPI

end program
