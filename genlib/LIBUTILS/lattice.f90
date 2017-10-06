module latticemod

 use unitcell_rout
 use latticerout
 use genvar
 use linalg 
 use lattice_symmetries, only : allocatesym

 IMPLICIT NONE

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
!********************************************
!********************************************
!********************************************

subroutine allocate_web(c,FLAG_OPEN_BC,FLAG_OPEN_BUT_PERIODIC_LINKED,FLAG_UNIT_CELL_INSIDE,FLAG_PLOT_3D,FLAG_LAT_T_COEF, &
                        & BORD,CROT,T1,T2,T3,routbuild,routbuild_,killsym,killallsym,killb,local_included, &
                        & dist_only,orbcut,cut_dist,rmaxdistr,k_space_only,reduceBZ,nx,ny,nz,nosplit)
implicit none

 type(web),intent(inout) :: c

 integer          :: N,i,j
 logical          :: notrans,killsym,FLAG_OPEN_BC,FLAG_OPEN_BUT_PERIODIC_LINKED,FLAG_UNIT_CELL_INSIDE,FLAG_PLOT_3D,FLAG_LAT_T_COEF,killallsym
 logical          :: killb,local_included,dist_only,k_space_only,reduceBZ,killallsymf
 real(8)          :: rmaxdistr,cut_dist
 integer          :: orbcut,nx,ny,nz
 real(8)          :: T1back(3),T2back(3),T3back(3),BORD(2),CROT(2),T1(3),T2(3),T3(3)
 logical,optional :: nosplit
 
interface
 subroutine routbuild_(in)
 use unitcell_rout
 implicit none
   type(unitcell),intent(inout) :: in
 end subroutine
 subroutine routbuild(in)
 use latticerout
 implicit none
   type(web),intent(inout) :: in
 end subroutine
end interface

 write(*,*) 'allocate web routine, first destroy lattice'
 call web_destructor(c)

 c%maxd=40
 c%maxvoisins=100
 c%maxconn=100
 c%FLAG_OPEN_BC=FLAG_OPEN_BC
 c%FLAG_OPEN_BUT_PERIODIC_LINKED=FLAG_OPEN_BUT_PERIODIC_LINKED
 c%FLAG_UNIT_CELL_INSIDE=FLAG_UNIT_CELL_INSIDE
 c%FLAG_PLOT_3D=FLAG_PLOT_3D
 notrans=c%FLAG_OPEN_BC.and..not.c%FLAG_OPEN_BUT_PERIODIC_LINKED

 write(*,*) 'define lattice quantities'
 c%teta(1)          =  BORD(1)     !periodic/anti-  boundary along T1
 c%teta(2)          =  BORD(2)     !periodic/anti-  boundary along T2
 c%symmetries       = .true.
 c%sym%centerrot    =  0.
 c%sym%centerrot(1) =  CROT(1)
 c%sym%centerrot(2) =  CROT(2)
 c%centre(1)        =  CROT(1)
 c%centre(2)        =  CROT(2)
 c%centre(3)        =  0.

 write(*,*) 'kill unitcell'
 call kill_unitcell(c%cell)
 write(*,*) 'build unitcell'
 call routbuild_(c%cell)
 if(maxval(abs(T1))<1.d-6) then
  c%T1 = [1.d0,0.d0,0.d0]
 else
  c%T1 = T1; 
 endif
 if(maxval(abs(T2))<1.d-6) then
  c%T2 = [0.d0,1.d0,0.d0]
 else
  c%T2 = T2;
 endif
 if(maxval(abs(T3))<1.d-6) then
  c%T3 = [0.d0,0.d0,1.d0]
 else
  c%T3 = T3;
 endif

 if(FLAG_LAT_T_COEF) then
   T1back=T1
   T2back=T2
   T3back=T3
   c%T1 = T1back(1)*c%cell%a+T1back(2)*c%cell%b+T1back(3)*c%cell%c
   c%T2 = T2back(1)*c%cell%a+T2back(2)*c%cell%b+T2back(3)*c%cell%c
   c%T3 = T3back(1)*c%cell%a+T3back(2)*c%cell%b+T3back(3)*c%cell%c
 endif
 if(FLAG_OPEN_BC)Then
  if(.not.FLAG_OPEN_BUT_PERIODIC_LINKED)then
    write(*,*) 'OPEN BOUNDARY CONDITION'
    c%open1=1; c%open2=1; c%open3=1
    c%T1(2)=0.
    c%T2(1)=0.
    c%T1(1)=c%T1(1)+0.5
    c%T2(2)=c%T2(2)+0.5
  endif
 endif
 
 write(*,*) 'build lattice'
 call routbuild(c)
 
 if(norme(T1)<1.d-4) T1=c%T1
 if(norme(T2)<1.d-4) T2=c%T2
 if(norme(T3)<1.d-4) then
  if(norme(c%T3)<1.d-4)then
   c%T3=[0.d0,0.d0,1.d0]
  endif
  T3=c%T3
 endif

 if(FLAG_OPEN_BC.and..not.FLAG_OPEN_BUT_PERIODIC_LINKED)then
    if(c%open1/=1.and.c%open2/=1.and.c%open3/=1)then
      write(*,*) 'building cluster with open boundary conditions'
      write(*,*) 'most probably something went wrong, check it out'
      write(*,*) 'if everything is fine, please remove this stoppoint'
      stop
    endif
 endif

 if(messages) write(*,*) 'cluster FIXED : boundary conditions / Center'

 if(.not.k_space_only)then 
   call findN(c)
   N=c%N
 else
   N=c%cell%q1
 endif

 if(N==0)then
   write(*,*) 'Found 0 site in lattice,error'
   write(*,*) ' N = ', N
   stop 'critical'
 endif

 if(.not.k_space_only)then
    killallsymf=killallsym.or.(c%FLAG_OPEN_BC.and..not.c%FLAG_OPEN_BUT_PERIODIC_LINKED)
    write(*,*) ' ---> allocate sym with killallsym : ', killallsymf
    call allocatesym(c%sym,c%N,c%cell%q1,notrans,killsym,killallsymf)
    write(*,*) 'danger, try to allocate large array : periodist, '
    allocate(c%periodist(c%N,c%N),c%periodist_open_bc(c%N,c%N),c%unitary_transform(2*c%N,2*c%N))
    write(*,*) '.....done......'
 else
    allocate(c%periodist(1,1),c%periodist_open_bc(1,1),c%unitary_transform(1,1))
 endif

 allocate(c%distances(c%N),c%distancesp(c%cell%q1,c%maxd))
 allocate(c%angles(c%N),c%latticemap(c%N),c%latticemapi(c%N))
 allocate(c%direc(c%N,c%N),c%distbord(c%N,4))
 allocate(c%bordure(c%N),c%whichbord(c%N))
 allocate(c%longvoisin(c%N,0:c%maxd,0:c%maxvoisins))
 allocate(c%nombrevoisin(c%N,0:c%maxd))
 allocate(c%ineigh(c%N,0:maxval(c%cell%nneigh)),c%phase(c%N,0:maxval(c%cell%nneigh)))
 allocate(c%nlink_(c%cell%q1,maxval(c%cell%nneigh)))
 allocate(c%maxdistsite(c%N,0:c%cell%q1))
 c%maxdistsite=0.d0
 allocate(c%centerCu(c%N,0:maxval(c%cell%nneigh)))
 allocate(c%phase4(c%N,maxval(c%cell%nneigh)),c%centreangle(c%N,maxval(c%cell%nneigh)),c%centredist(c%N))
 allocate(c%struct_link(c%N,c%cell%struct_nlink,2,c%cell%struct_nm))
 allocate(c%struct_pos(c%N,c%cell%struct_nm,3))
 allocate(c%cadran(c%N,maxval(c%cell%nneigh)))
 allocate(c%x(c%N,3),c%k_inv(c%N/c%cell%q1))
 allocate(c%site(0:c%N),c%site2(0:c%N),c%site3(0:c%N),c%xma(c%N,3))
 allocate(c%conv_unitcell_site_(c%cell%q1))
 if(allocated(c%major_site)) deallocate(c%major_site)
 allocate(c%major_site(c%N))

 if(messages) write(*,*) 'build unitcell'
 call unitcell_build(c%cell)
 if(messages) write(*,*) '---> done <---'

 if(.not.k_space_only)then
   if(messages) write(*,*) 'build cluster'
   if(.not.allocated(c%vecin))then
     call latticebuild(c)
   else
     call latticebuild(c,c%vecin)
   endif
   if(messages) write(*,*) 'build periodic structure'
   call periodist_init(c)
   call periodist_init_open_bc(c)
   if(messages) write(*,*) 'build max dist. between sites ...'
   call build_maxdistsite(c)
   if(messages) write(*,*) 'lattice built, now build neighbors .....'
   call buildneightable(c)
   if(messages) write(*,*) 'neigh table built ... now do reciprocal space ....'
   write(*,*) 'maxd_real is : ',c%maxd_real
   if(c%maxd_real==0) then
     write(*,*) 'maxd_real is 0!!! CRITICAL'
     c%maxd_real=1
    !stop 
   endif
 endif

 write(*,*) '---> start building reciprocal lattice <---'
 if(.not.k_space_only)then
  call reciproc(c%cell,c%T1,c%T2,c%T3,BORD,reduceBZ)
 else
  call reciproc(c%cell,c%T1,c%T2,c%T3,BORD,reduceBZ,nx,ny,nz)
 endif

 call copy_unitcellvar_to_cluster(c)

 if(messages) write(*,*) 'reciprocal space built'

 if(.not.k_space_only)then
  call fix_nlink_(c)
  call simple_test_lattice(c)
  call buildstruct(c)
  if(messages) write(*,*) 'build symtab'
  call buildsymtab(killb,dist_only,local_included,orbcut,cut_dist,rmaxdistr,c,c%sym,nosplit=nosplit)
  if(messages) write(*,*) 'symmetries done.....'
  if(messages) write(*,*) 'build center site'
  call build_center_cu(c)
  call build_usefull_arrays(c)

   if(c%N<100000) then
    call check_lattice_no_ambiguous_link(c)
    call dump_links_in_file(c)
   endif
 endif

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

end module
