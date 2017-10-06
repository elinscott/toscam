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
implicit none

interface
 subroutine  chi_tilde_loc( &
 & op, w1,  w2,  w3, PHI_EPS, beta, Z, nsites, nup, ndn,  &
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
real(8) :: w1,w2,w3,PHI_EPS,beta,Z
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
real(8) :: cp_i_E(dim_E_i)
real(8) :: cp_pup_E(dim_E_pup)
real(8) :: cp_pdn_E(dim_E_pdn)
real(8) :: cp_mup_E(dim_E_mup)
real(8) :: cp_mdn_E(dim_E_mdn)
real(8) :: cp_p2dn_E(dim_E_p2dn)
real(8) :: cp_m2dn_E(dim_E_m2dn)
real(8) :: cp_puppdn_E(dim_E_puppdn)
real(8) :: cp_muppdn_E(dim_E_muppdn)
real(8) :: cp_pupmdn_E(dim_E_pupmdn)
real(8) :: cp_mupmdn_E(dim_E_mupmdn)
integer :: n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28
real(8) :: cp_i_cdup(n1,n2)       
real(8) :: cp_i_cddn(n3,n4)       
real(8) :: cp_pup_cddn(n5,n6)     
real(8) :: cp_pdn_cdup(n7,n8)     
real(8) :: cp_pdn_cddn(n9,n10)     
real(8) :: cp_mup_cdup(n11,n12)     
real(8) :: cp_mup_cddn(n13,n14)     
real(8) :: cp_mdn_cddn(n15,n16)     
real(8) :: cp_mdn_cdup(n17,n18)     
real(8) :: cp_muppdn_cdup(n19,n20)  
real(8) :: cp_pupmdn_cddn(n21,n22)  
real(8) :: cp_mupmdn_cdup(n23,n24)  
real(8) :: cp_mupmdn_cddn(n25,n26)  
real(8) :: cp_m2dn_cddn(n27,n28)   
 end subroutine
end interface

integer,parameter                   :: nw1=10,nw2=10,nw3=10
real(8),parameter                   :: w1max=3.14,w2max=3.14,w3max=3.14
complex(8)                          :: Chi_lattice(nw1,nw2,nw3)
real(8)                             :: w1,w2,w3
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

     PHI_EPS=1.d-9

     open(unit=1414,file='chilattice_vertex',form='unformatted')
     read(1414) ZZ,beta,nsites,min_up,max_up,min_dn,max_dn,nsector

     write(*,*) 'partition function : ', ZZ
     write(*,*) 'beta               : ', beta
     write(*,*) 'nsites             : ', nsites
     write(*,*) 'nsector            : ', nsector

     allocate(nup(nsector),ndn(nsector))
     if(allocated(cup_mat)) deallocate(cup_mat)
     if(allocated(cdn_mat)) deallocate(cdn_mat)
     allocate(cup_mat(min_up:max_up,min_dn:max_dn))
     allocate(cdn_mat(min_up:max_up,min_dn:max_dn))
     allocate(eigenval(min_up:max_up,min_dn:max_dn))
     allocate(neigentab(min_up:max_up,min_dn:max_dn))

     write(*,*) 'start reading' 
     do i=1,nsector
       write(*,*) 'reading sector / nsector: ', i,nsector
       read(1414) nup(i),ndn(i)
       read(1414) neigen
       neigentab(nup(i),ndn(i))=neigen
       allocate(eigenval(nup(i),ndn(i))%val(neigen))
       read(1414) eigenval(nup(i),ndn(i))%val

       write(*,*) 'nup,ndn   : ',nup(i),ndn(i) 
       write(*,'(a,200f9.3)') 'eigenvals : ',eigenval(nup(i),ndn(i))%val 
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


    !---------------------------!
    !... for a fix frequency ...!
    !---------------------------! 

     w1=0.
     w2=0.
     w3=0.
     Chi_lattice=0.d0
     op=1

!========================================================================!
!========================================================================!
!========================================================================!
!========================================================================!
!========================================================================!

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


!TO UPDATE BEGIN
 !    dim_E_i=      
 !    dim_E_pup=    
 !    dim_E_pdn=    
 !    dim_E_mup=    
 !    dim_E_mdn=    
 !    dim_E_p2dn=   
 !    dim_E_m2dn=   
 !    dim_E_puppdn= 
 !    dim_E_muppdn= 
 !    dim_E_pupmdn= 
 !    dim_E_mupmdn= 
 !    n1=
 !    n2=
 !    n3=
 !    n4=
 !    n5=
 !    n6=
 !    n7=
 !    n8=
 !    n9=
 !    n10=
 !    n11=
 !    n12=
 !    n13=
 !    n14=
 !    n15=
 !    n16=
 !    n17=
 !    n18=
 !    n19=
 !    n20=
 !    n21= 
 !    n22=
 !    n23=
 !    n24=
 !    n25=
 !    n26=
 !    n27=
 !    n28=
!TO UPDATE ENDS

     if(allocated(cp_i_E)) deallocate(cp_i_E, cp_pup_E, cp_pdn_E, cp_mup_E, cp_mdn_E, cp_p2dn_E, cp_m2dn_E, cp_puppdn_E, &
& cp_muppdn_E, cp_pupmdn_E, cp_mupmdn_E, cp_i_cdup,cp_i_cddn,                        & 
& cp_pup_cddn, cp_pdn_cdup, cp_pdn_cddn, cp_mup_cdup, cp_mup_cddn, cp_mdn_cddn,      & 
& cp_mdn_cdup, cp_muppdn_cdup, cp_pupmdn_cddn, cp_mupmdn_cdup, cp_mupmdn_cddn, cp_m2dn_cddn)

     allocate( cp_i_E(dim_E_i) )
     allocate( cp_pup_E(dim_E_pup) )
     allocate( cp_pdn_E(dim_E_pdn) )
     allocate( cp_mup_E(dim_E_mup) )
     allocate( cp_mdn_E(dim_E_mdn) )
     allocate( cp_p2dn_E(dim_E_p2dn) )
     allocate( cp_m2dn_E(dim_E_m2dn) )
     allocate( cp_puppdn_E(dim_E_puppdn) )
     allocate( cp_muppdn_E(dim_E_muppdn) )
     allocate( cp_pupmdn_E(dim_E_pupmdn) )
     allocate( cp_mupmdn_E(dim_E_mupmdn) )
     allocate( cp_i_cdup(n1,n2) )
     allocate( cp_i_cddn(n3,n4) )
     allocate( cp_pup_cddn(n5,n6) )
     allocate( cp_pdn_cdup(n7,n8) )
     allocate( cp_pdn_cddn(n9,n10) )
     allocate( cp_mup_cdup(n11,n12) )
     allocate( cp_mup_cddn(n13,n14) )
     allocate( cp_mdn_cddn(n15,n16) )
     allocate( cp_mdn_cdup(n17,n18) )
     allocate( cp_muppdn_cdup(n19,n20) )
     allocate( cp_pupmdn_cddn(n21,n22) )
     allocate( cp_mupmdn_cdup(n23,n24) )
     allocate( cp_mupmdn_cddn(n25,n26) )
     allocate( cp_m2dn_cddn(n27,n28) )

!TO UPDATE BEGINS
 !    cp_i_E=      
 !    cp_pup_E=    
 !    cp_pdn_E=    
 !    cp_mup_E=    
 !    cp_mdn_E=    
 !    cp_p2dn_E=   
 !    cp_m2dn_E=   
 !    cp_puppdn_E= 
 !    cp_muppdn_E= 
 !    cp_pupmdn_E= 
 !    cp_mupmdn_E= 

 !    cp_i_cdup=
 !    cp_i_cddn=
 !    cp_pup_cddn=
 !    cp_pdn_cdup=
 !    cp_pdn_cddn=
 !    cp_mup_cdup=
 !    cp_mup_cddn=
 !    cp_mdn_cddn=
 !    cp_mdn_cdup=
 !    cp_muppdn_cdup=
 !    cp_pupmdn_cddn=
 !    cp_mupmdn_cdup=
 !    cp_mupmdn_cddn=
 !    cp_m2dn_cddn=
!TO UPDATE ENDS

      call chi_tilde_loc( &
 &     op, w1,  w2,  w3, PHI_EPS, beta, ZZ, nsites, nup(i), ndn(i), & 
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
 &    cp_m2dn_cddn,&
 &    pDNDN,pUPDN )

      Chi_lattice(1,1,1)=Chi_lattice(1,1,1) + pDNDN + pUPDN

     enddo

!========================================================================!
!========================================================================!
!========================================================================!
!========================================================================!
!========================================================================!

     write(*,*) 'Chi is : ' , Chi_lattice(1,1,1)

end program
