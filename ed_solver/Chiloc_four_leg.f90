module vertex

  use correlations
  use namelistmod

PRIVATE

PUBLIC :: four_leg_vertex_matrices_routine

  INTERFACE ASSIGNMENT (=)
    MODULE PROCEDURE cdagger_copy
  END INTERFACE

  TYPE cdagger_mat
   integer    :: k_p,l_p
   integer    :: k_m,l_m
#ifdef _complex
   complex(8),allocatable :: c_p(:,:) , c_m(:,:)
#else
   real(8),   allocatable :: c_p(:,:) , c_m(:,:)
#endif 
  END TYPE

  logical,parameter :: verbose=.false.

contains

!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!

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

!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!

  subroutine four_leg_vertex_matrices_routine(AIM,GS)
  implicit none
  TYPE(eigensectorlist_type)                :: GS
  TYPE(cdagger_mat)                         :: cup(GS%nsector),cdn(GS%nsector)
  TYPE(cdagger_mat),allocatable             :: cup_mat(:,:),cdn_mat(:,:)
  TYPE(AIM_type)                            :: AIM
  integer                                   :: itot,i,j,k,nup,ndn,min_up,max_up,min_dn,max_dn
  real(8)                                   :: beta,ZZ

     CALL init_apply_C (AIM) 
     itot = AIM%bath%Nb+AIM%impurity%Nc
     ZZ   = partition(beta_ED,GS)

     call four_leg_vertex_matrices(AIM,GS,Cup_sector,apply_Cup,cup)
     call four_leg_vertex_matrices(AIM,GS,Cdo_sector,apply_Cdo,cdn)

     min_up=10000
 max_up=-10000
 min_dn=10000
 max_dn=-10000

      do i=1,GS%nsector
        nup=GS%es(i)%sector%updo%up%npart
        ndn=GS%es(i)%sector%updo%down%npart
        if(nup>max_up) max_up=nup
 if(ndn>max_dn) max_dn=ndn
 if(nup<min_up) min_up=nup
 if(ndn<min_dn) min_dn=ndn
      enddo
     write(*,*) 'nup ranging : ',min_up,max_up
     write(*,*) 'nup ranging : ',min_dn,max_dn
     if(allocated(cup_mat)) deallocate(cup_mat)
     if(allocated(cdn_mat)) deallocate(cdn_mat)
     allocate(cup_mat(min_up:max_up,min_dn:max_dn))
     allocate(cdn_mat(min_up:max_up,min_dn:max_dn))

     do i=1,GS%nsector
       nup=GS%es(i)%sector%updo%up%npart
       ndn=GS%es(i)%sector%updo%down%npart
       cup_mat(nup,ndn)=cup(i)
       cdn_mat(nup,ndn)=cdn(i)
     enddo

     open(unit=1414,file='chiloc_vertex',form='unformatted')
     write(1414) ZZ,beta_ED,itot,min_up,max_up,min_dn,max_dn,GS%nsector
     do i=1,GS%nsector
       nup=GS%es(i)%sector%updo%up%npart
       ndn=GS%es(i)%sector%updo%down%npart
       write(1414) nup,ndn
       write(1414) GS%es(i)%lowest%neigen
       write(1414) GS%es(i)%lowest%eigen(:)%val
       write(1414) cup(i)%k_p,cup(i)%l_p
       write(1414) cup(i)%c_p 
       write(1414) cup(i)%k_m,cup(i)%l_m
       write(1414) cup(i)%c_m
       write(1414) cdn(i)%k_p,cdn(i)%l_p
       write(1414) cdn(i)%c_p
       write(1414) cdn(i)%k_m,cdn(i)%l_m
       write(1414) cdn(i)%c_m
     enddo
     close(1414)

     write(*,*) 'cleaning'
     do i=min_up,max_up
     do j=min_dn,max_dn
      call kill_cdagger(cup_mat(i,j))
      call kill_cdagger(cdn_mat(i,j))
     enddo
     enddo
     deallocate(cup_mat)
     deallocate(cdn_mat)
     do i=1,GS%nsector
       call kill_cdagger(cup(i))
 call kill_cdagger(cdn(i))
     enddo

  return
  end subroutine

!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!

  SUBROUTINE four_leg_vertex_matrices(AIM,GS,Asector,applyA,cdagger)
  implicit none
    TYPE(eigensectorlist_type)                :: GS
    TYPE(cdagger_mat)                         :: cdagger(GS%nsector)
    TYPE(AIM_type)                            :: AIM
    TYPE(sector_type)                         :: Asec
    TYPE(eigensector_type)                    :: Apm_es(2)
    TYPE(eigensector_type), POINTER           :: es    => NULL()
    TYPE(eigen_type),       POINTER           :: eigen => NULL(), OPeigen => NULL()
    INTEGER                                   :: i1,j1,isec_back,ii,i,i_size,j_size
    INTEGER                                   :: ipm,jpm,kpm,iph,iorb,jorb,jjj1,jjj2,j
    LOGICAL                                   :: orb(1)
    INTEGER                                   :: isector,ieigen,uup,ddn,itot,i_,v(2),issz
    INTEGER                                   :: ktot,iiorb,jjorb,kk,jj,iorb_f,jorb_f

   !----------------------------------------------------!
    INTERFACE
    ! RETURNS SECTOR OF A,A^+|0> 
      SUBROUTINE Asector(Asec,pm,sector)
        use sector_class, only : sector_type
       TYPE(sector_type),      INTENT(INOUT) :: Asec
       TYPE(sector_type),      INTENT(IN)    :: sector
       CHARACTER(LEN=1),       INTENT(IN)    :: pm
      END SUBROUTINE
   !----------------------------------------------------!
    ! COMPUTES A,A^+|0>
      SUBROUTINE applyA(Aes,pm,MASK,es,rankr)
        use eigen_sector_class , only : eigensector_type
       TYPE(eigensector_type), INTENT(INOUT) :: Aes
       CHARACTER(LEN=1),       INTENT(IN)    :: pm
       TYPE(eigensector_type), INTENT(IN)    :: es
       INTEGER,                INTENT(IN)    :: rankr
       LOGICAL,                INTENT(IN)    :: MASK(:)
      END SUBROUTINE
    END INTERFACE
   !----------------------------------------------------!

    itot = AIM%bath%Nb+AIM%impurity%Nc

!=====================================================================!
!=====================================================================!
!=====================================================================!

    orb=.true.

    DO isector=1,GS%nsector
         es => GS%es(isector)
         uup  =  GS%es(isector)%sector%updo%up%npart
         ddn  =  GS%es(isector)%sector%updo%down%npart
         DO i=1,2
             CALL Asector(Asec,pm(i),es%sector)
             CALL new_eigensector(Apm_es(i),Asec)
         ENDDO

         if(verbose) write(*,*) '====================================='
         if(verbose) write(*,*) 'SECTOR , nup, ndn : ', isector,uup,ddn
         if(verbose) write(*,*) '====================================='

     !----------------------------------------------------!
         DO kpm=1,2

            if(kpm==1.and.verbose) write(*,*) '-------------C^+----------------'
            if(kpm==2.and.verbose) write(*,*) '--------------c-----------------'

            CALL delete_eigenlist(Apm_es(kpm)%lowest)
            if(verbose) write(*,*) 'starting loop over states in sector, there are [x] states : ', GS%es(isector)%lowest%neigen

            DO ieigen=1,GS%es(isector)%lowest%neigen
               CALL applyA(Apm_es(kpm),pm(kpm),orb,es,ieigen)
               if(verbose) write(*,*) 'applying A, neigen = ',Apm_es(kpm)%lowest%neigen
            ENDDO

            do ii=1,GS%nsector
              if(GS%es(ii)%sector%updo%up%npart==Apm_es(kpm)%sector%updo%up%npart.and. &
               & GS%es(ii)%sector%updo%down%npart==Apm_es(kpm)%sector%updo%down%npart)then
              exit
             endif
            enddo

            if(ii==GS%nsector+1)then
              if(verbose) write(*,*) '--sector does not exist--'
              if(verbose) write(*,*) 'original up and down : ', uup,ddn
              if(verbose) write(*,*) ' kpm = ', kpm
              if(kpm==1)then
                cdagger(isector)%k_p=1
                cdagger(isector)%l_p=1
                call allocate_dagger_p(cdagger(isector))
                cdagger(isector)%c_p(1,1)=0.
              else
                cdagger(isector)%k_m=1
                cdagger(isector)%l_m=1
                call allocate_dagger_m(cdagger(isector))
                cdagger(isector)%c_m(1,1)=0.
              endif
              goto 44
            endif

            i_size=GS%es(ii     )%lowest%neigen
            j_size=GS%es(isector)%lowest%neigen

            if(verbose) then
               write(*,*) 'applyA is in sector and max is         : ',ii,GS%nsector
               write(*,*) 'in target sector there are [x] states  : ',GS%es(ii)%lowest%neigen
               write(*,*) 'up,dn in original sector               : ',GS%es(isector)%sector%updo%up%npart,GS%es(isector)%sector%updo%down%npart
               write(*,*) 'up,dn in target sector                 : ',GS%es(ii)%sector%updo%up%npart,GS%es(ii)%sector%updo%down%npart
               write(*,*) 'there are [x] collected states from A  : ',Apm_es(kpm)%lowest%neigen
               write(*,*) 'destroy or creates                     : ',kpm,pm(kpm)
               write(*,*) 'target/initial sector sizes            : ',i_size,j_size
            endif

            if(Apm_es(kpm)%lowest%neigen/=GS%es(isector)%lowest%neigen)then
             write(*,*) 'obtained c|i> and initial sector size : ',Apm_es(kpm)%lowest%neigen,GS%es(isector)%lowest%neigen
             write(*,*) 'they should be the same, error'
             stop
            endif
 
            if(kpm==1)then
              cdagger(isector)%k_p=i_size 
              cdagger(isector)%l_p=j_size
              call allocate_dagger_p(cdagger(isector))
              do i=1,i_size
               do j=1,j_size
                cdagger(isector)%c_p(i,j)= scalprod( GS%es(ii)%lowest%eigen(i)%vec%rc , Apm_es(kpm)%lowest%eigen(j)%vec%rc)
               enddo
              enddo
              if(verbose) write(*,*) 'maxval(cdagger) : ', maxval(abs(cdagger(isector)%c_p))
            else
              cdagger(isector)%k_m=i_size 
              cdagger(isector)%l_m=j_size
              call allocate_dagger_m(cdagger(isector))
              do i=1,i_size
               do j=1,j_size
                cdagger(isector)%c_m(i,j)= scalprod( GS%es(ii)%lowest%eigen(i)%vec%rc , Apm_es(kpm)%lowest%eigen(j)%vec%rc)
               enddo
              enddo
              if(verbose) write(*,*) 'maxval(cdagger) : ', maxval(abs(cdagger(isector)%c_m))
            endif
           
   44    ENDDO
     !----------------------------------------------------!

    ENDDO
!=====================================================================!
!=====================================================================!
!=====================================================================!

    CALL delete_sector(Asec)
    DO iph=1,2
      CALL delete_eigensector(Apm_es(iph))
    ENDDO

return
end subroutine

!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!
!####################################################!

end module
