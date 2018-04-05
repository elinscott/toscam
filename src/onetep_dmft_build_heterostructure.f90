!**************************************!
!**************************************!
!**************************************!
!**************************************!
!**************************************!
!**************************************!

module variables

use namelistmod
use StringManip, only : tostring
integer,parameter         :: nmaxentry=10
real(8),parameter         :: hartreeunits=27.2113839
real(8)                   :: separation,cutoff,same_atom_distance,planeshift
real(8)                   :: cell1shiftd1,cell1shiftd2,cell1shiftd3 
real(8)                   :: cell2shiftd1,cell2shiftd2,cell2shiftd3,dist_neigh,cutoff_rot,radius_atom
character(20)             :: filename1,filename2
logical                   :: remove_duplicate_unitcell,dmft_case,round_down_planes,symmetrize_structure,brokensym,check1
real(8)                   :: onlycorplanes
integer                   :: sandwitch,jj,jjj,nfrequ
real(8)                   :: open_boundary,fac,myaxis(3,3),cutoff_boundary_shifted
integer(4)                :: nplane1,nplane2,width
integer                   :: d3,d1,d2,planes_in_cell1,planes_in_cell2,rutile_dir
integer                   :: nbasis
integer                   :: ncor,firstplane,lastplane,obtainedplanes
real(8),allocatable       :: corU1(:),corU2(:),corJ1(:),corJ2(:)
character(20),allocatable :: atomsymbol(:)
character(40),allocatable :: atomConf(:)
character(20),allocatable :: corsymbol1(:),corsymbol2(:)
logical,allocatable       :: isitcor(:)
integer,allocatable       :: basistot(:),basisZ(:)
integer                   :: kkk,duplicatecellchoice
real(8)                   :: planer1,planer2
logical                   :: periodic_plane

contains


 subroutine set_var
  type(namelist_set) :: nm

   call namelist_init(nm,400,name_of_namelist='build_oxides')
   call putel_in_namelist(nm,duplicatecellchoice,'duplicatecellchoice',0,'Definition: if /=0, the code will choose the upper plane of the unitcell when there are duplicates at the interface belonging to different cells')
   call putel_in_namelist(nm,brokensym,'brokensym',.false.,'Definition: if true will consider a broken symmetry state, and the SELF ENERGY will be repeated with this pattern')
   call putel_in_namelist(nm,dmft_case,'dmft_case',.true.,'Definition: if true it will generate the case.dat file for dmft run') 
   call putel_in_namelist(nm,nfrequ,'nfrequ',156,'Definition: number of matsubara frequencies')
   call putel_in_namelist(nm,symmetrize_structure,'symmetrize_structure',.false.,'Definition: if true symmetrizes the sandwitch structure to minimize artefacts in dft')
   call putel_in_namelist(nm,filename1,'filename1','latio3.input','Definition: ')
   call putel_in_namelist(nm,filename2,'filename2','srtio3.input','Definition: ')
   call putel_in_namelist(nm,nplane1,'nplane1',2,'Definition: ')
   call putel_in_namelist(nm,nplane2,'nplane2',2,'Definition: ')
   call putel_in_namelist(nm,width,'width',2,'Definition: ')
   call putel_in_namelist(nm,remove_duplicate_unitcell,'remove_duplicate_unitcell',.true.,'Definition: ')
   call putel_in_namelist(nm,separation,'separation',2.85d0,'Definition: ')
   call putel_in_namelist(nm,d3,'d3',1,'Definition:heterostructure is stacked along d3')
   call putel_in_namelist(nm,d1,'d1',2,'Definition:')
   call putel_in_namelist(nm,d2,'d2',3,'Definition:')
   call putel_in_namelist(nm,cutoff,'cutoff',0.2d0,'Definition: same site distance cutoff')
   call putel_in_namelist(nm,cell1shiftd1,'cell1shiftd1',0.d0,'Definition:')   
   call putel_in_namelist(nm,cell1shiftd2,'cell1shiftd2',0.d0,'Definition:') 
   call putel_in_namelist(nm,cell1shiftd3,'cell1shiftd3',0.d0,'Definition:') 
   call putel_in_namelist(nm,cell2shiftd1,'cell2shiftd1',0.d0,'Definition:')  
   call putel_in_namelist(nm,planeshift,'planeshift',0.05d0,'Definition')
   call putel_in_namelist(nm,cell2shiftd2,'cell2shiftd2',0.d0,'Definition:') 
   call putel_in_namelist(nm,cell2shiftd3,'cell2shiftd3',0.d0,'Definition:')    
   call putel_in_namelist(nm,same_atom_distance,'same_atom_distance',0.5d0,'Definition:')
   call putel_in_namelist(nm,planes_in_cell1,'planes_in_cell1',2,'Definition:')
   call putel_in_namelist(nm,planes_in_cell2,'planes_in_cell2',2,'Definition:')
   call putel_in_namelist(nm,nbasis,'nbasis',4,'Definition: ')
   call putel_in_namelist(nm,ncor,'ncor',    1,'Definition: ')
   call putel_in_namelist(nm,dist_neigh,'dist_neigh',2.8d0,'Definition: max dist between nearest neighbor')
   call putel_in_namelist(nm,periodic_plane,'periodic_plane',.false.,'Definition: if true, when builds the planes, it will include atoms round the boundaries')
   call putel_in_namelist(nm,onlycorplanes,'onlycorplanes',-1.d0,'Definition: only keep atoms within Ti planes, this is the cutoff which includes the neighboring Oxygen, if negative, does not use this scheme')
   call putel_in_namelist(nm,round_down_planes,'round_down_planes',.false.,'Definition: if true the planes are rounded down to lowest integer, rather than rounded to nearest integer')
   call putel_in_namelist(nm,cutoff_rot,'cutoff_rot',0.3d0,'Definition: cutoff to consider atom as part of the plane perpendicular to apical direction for rotations')
   call putel_in_namelist(nm,rutile_dir,'rutile_dir',3,'Definition: direction of rutile axis, by default it is along the d3 direction, heterostructure direction')
   call putel_in_namelist(nm,cutoff_boundary_shifted,'cutoff_boundary_shifted',0.02d0,'Definition: will shift all atoms within this distance from the top boundary (cap), and put them back near the bottom boundary, if negative it does the opposite')

   allocate(corU1(ncor),corJ1(ncor),corU2(ncor),corJ2(ncor),corsymbol1(ncor),corsymbol2(ncor))
   allocate(atomsymbol(nbasis),basistot(nbasis),basisZ(nbasis),atomConf(nbasis)) 
   do i=1,nbasis
    call putel_in_namelist(nm,atomsymbol(i),'atomsymbol'//trim(adjustl(tostring(i))),' ','Definition: ')
    call putel_in_namelist(nm,basistot(i),'Ne'//trim(adjustl(tostring(i))),1,'Definition: ')
    call putel_in_namelist(nm,basisZ(i),'Ze'//trim(adjustl(tostring(i))),1,'Definition: ')
    call putel_in_namelist(nm,atomConf(i),'conf'//trim(adjustl(tostring(i))),' ','Definition: ')
   enddo
   do i=1,ncor   
    call putel_in_namelist(nm,corsymbol1(i),'corsymbolA'//trim(adjustl(tostring(i))),' ' ,'Definition: ')
    call putel_in_namelist(nm,corsymbol2(i),'corsymbolB'//trim(adjustl(tostring(i))),' ' ,'Definition: ')
    call putel_in_namelist(nm,corU1(i)     ,      'corUA'//trim(adjustl(tostring(i))),0.d0,'Definition: ')
    call putel_in_namelist(nm,corJ1(i)     ,      'corJA'//trim(adjustl(tostring(i))),0.d0,'Definition: ')
    call putel_in_namelist(nm,corU2(i)     ,      'corUB'//trim(adjustl(tostring(i))),0.d0,'Definition: ')
    call putel_in_namelist(nm,corJ2(i)     ,      'corJB'//trim(adjustl(tostring(i))),0.d0,'Definition: ')
   enddo
   call putel_in_namelist(nm,radius_atom,'radius_atom',2.7d0,'Definition: radius of atoms in onetep')
   call putel_in_namelist(nm,open_boundary,'open_boundary',0.0d0,'Definition: open boundary conditions')
   call putel_in_namelist(nm,sandwitch,'sandwitch',0,'Definition: if 1 creates a sandwitch structure (structure 2 in the middle), for embedding, if 2 creates a sandwitch structure with structure 1 in the center')
   call putel_in_namelist(nm,myaxis(:,1),'myaxis1',[1.d0,-1.d0,0.d0],'first axis for local rotations')
   call putel_in_namelist(nm,myaxis(:,2),'myaxis2',[1.d0, 1.d0,0.d0],'second axis for local rotations')
   call putel_in_namelist(nm,myaxis(:,3),'myaxis3',[0.d0, 0.d0,1.d0],'third axis for local rotations')
   call look_for_namelist_in_file(nm,'var.txt')
   call look_for_command_line_argument(nm)

   corU1=corU1/hartreeunits
   corU2=corU2/hartreeunits
   corJ1=corJ1/hartreeunits
   corJ2=corJ2/hartreeunits

   if(sandwitch==1)then
    if(nplane1<=1) then
      nplane1=2
      call system(" echo 'nplane1=2' >> var.txt ") 
    endif
    if(mod(nplane1,2)/=0) nplane1=nplane1+1
   elseif(sandwitch==2)then
    if(nplane2<=1) then
      nplane2=2
      call system(" echo 'nplane2=2' >> var.txt ")
    endif
    if(mod(nplane2,2)/=0) nplane2=nplane2+1
   endif 

   !----------------------------------------!
   !SANDWITCH 0:
   !  -|1-4 ||| 5-8 |-
   !SANDWITCH 1:
   ! -|3|4| ||| 5-8 (structure 2) ||| 1 | 2 |-
   !SANDWITCH 2:
   ! -|7|8| ||| 1-4 (structure 1) ||| 5 | 6 |-
   !----------------------------------------!

   if(sandwitch==0)then
     firstplane=1 
     lastplane=nplane2*planes_in_cell2 + nplane1*planes_in_cell1
   elseif(sandwitch==2)then
     firstplane=nplane2/2*planes_in_cell2 + nplane1*planes_in_cell1 + 1
      lastplane=firstplane-1
   elseif(sandwitch==1)then
     firstplane=nplane1/2*planes_in_cell1 + 1
     lastplane=firstplane-1
   elseif(sandwitch>2)then
     write(*,*) 'sandwitch > 2 not implemented'
     stop
   endif

   write(*,*) 'FIRST AND LAST PLANES ARE : ', firstplane,lastplane

   !Ne Ze
   !write(500,*) 'Sr  Sr  38 13 '//trim(adjustl(tostring(radius_atom)))  
   !write(500,*) 'Ti  Ti  22 10 '//trim(adjustl(tostring(radius_atom))) 
   !write(500,*) 'O    O   8  4 '//trim(adjustl(tostring(radius_atom))) 
   !write(500,*) 'La  La  57 20 '//trim(adjustl(tostring(radius_atom))) 
   !write(500,*) 'Al  Al  13  4 '//trim(adjustl(tostring(radius_atom)))  
 
 end subroutine

end module


PROGRAM oxides
 use variables
 IMPLICIT NONE
 integer                   :: natom
 integer                   :: i,j,k,ii
 real(8),allocatable       :: coordfull(:,:),coordfullcopy(:,:),coordhub(:,:),cell1(:,:),cell2(:,:),cell1copy(:,:),cell2copy(:,:)
 integer                   :: nhub,hubcounter,ncell1,ncell2,at
 character(3),allocatable  :: labelfullcopy(:),labelfull(:),labelcell1(:),labelcell2(:),labelcell1copy(:),labelcell2copy(:)
 integer,allocatable       :: planecell1(:),planecell2(:),planefullcopy(:,:),nentries(:),planefull(:,:),myatoms(:),hubequivalent(:)
 integer,allocatable       :: planeloc(:),planeloccopy(:),belongscell(:),belongscellcopy(:)
 integer,allocatable       :: symcell1(:),symcell2(:),symcell1copy(:),symcell2copy(:),symfullcopy(:),symfull(:)
 real(8)                   :: shift(3),trans1(3,3),trans2(3,3),tota(3),totb(3),totc(3),vecd3(3)
 logical                   :: duplicate
 real(8)                   :: shift1_sand(3),shift2_sand(3)
 character(6)              :: tempa

 call set_var


 vecd3=0.0 ; vecd3(d3)=1.d0
 
 ncell1    =  0
 ncell2    =  0 

 open(unit=100,file=trim(adjustl(filename1)))
 do 
  read(100,*,end=60)
  ncell1=ncell1+1
 enddo 
 60 continue
 ncell1=ncell1-1
 open(unit=100,file=trim(adjustl(filename2)))
 do
  read(100,*,end=61)
  ncell2=ncell2+1
 enddo
 61 continue
 ncell2=ncell2-1
 write(*,*) 'NCELL1 = ', ncell1
 write(*,*) 'NCELL2 = ', ncell2

 allocate(cell1(ncell1,3),cell2(ncell2,3))
 allocate(symcell1(ncell1),symcell2(ncell2))
 allocate(labelcell1(ncell1),labelcell2(ncell2))

 trans1=0. ; trans2=0.

 if(brokensym)then
  inquire(file=trim(adjustl(filename1))//".sym",exist=check1)
  if(.not.check1)then
   write(*,*) 'ERROR FILE SYM SHOULD BE PRESENT - using broken sym'
   stop
  endif 
  open(unit=100,file=trim(adjustl(filename1))//".sym")
   do i=1,ncell1
    read(100,*) tempa,symcell1(i)
    if(abs(symcell1(i))>1)then
     write(*,*) 'ERROR CELLSYM1 : now only implemented for +-1 for AF state'
     stop
    endif
   enddo
  close(100)
  inquire(file=trim(adjustl(filename2))//".sym",exist=check1)
  if(.not.check1)then
   write(*,*) 'ERROR FILE SYM FILE SHOULD BE PRESENT - using broken sym'
   stop
  endif 
  open(unit=100,file=trim(adjustl(filename2))//".sym")
   do i=1,ncell2
    read(100,*) tempa,symcell2(i)
    if(abs(symcell2(i))>1)then
     write(*,*) 'ERROR CELLSYM2 : now only implemented for +-1 for AF state'
     stop
    endif
   enddo
  close(100)
 endif

 open(unit=100,file=trim(adjustl(filename1)))
 read(100,*) trans1(1,1),trans1(2,2),trans1(3,3)
 do i=1,ncell1
  read(100,*) labelcell1(i),(cell1(i,j),j=1,3)
 enddo 
 close(100)
 open(unit=100,file=trim(adjustl(filename2)))
 read(100,*) trans2(1,1),trans2(2,2),trans2(3,3)
 do i=1,ncell2
  read(100,*) labelcell2(i),(cell2(i,j),j=1,3)
 enddo
 close(100)

 trans1=trans1*0.5291772084d0
 trans2=trans2*0.5291772084d0

 write(*,*) 'NCELL1 NCELL2 BEFORE CORRECTION : ', ncell1,ncell2
 write(*,*) 'TRANS1' , (trans1(i,i),i=1,3)
 write(*,*) 'TRANS2' , (trans2(i,i),i=1,3)

 if(.not.remove_duplicate_unitcell) goto 77

 allocate(cell1copy(ncell1,3),cell2copy(ncell2,3),symcell1copy(ncell1),symcell2copy(ncell2))
 allocate(labelcell1copy(ncell1),labelcell2copy(ncell2))
 ncell1=0
 ncell2=0
 cell1copy=cell1
 cell2copy=cell2
 symcell1copy=symcell1
 symcell2copy=symcell2
 labelcell1copy=labelcell1
 labelcell2copy=labelcell2

 write(*,*) 'CELL1 SCAN'
 do i=1,size(cell1copy,1)
  if(abs(cell1copy(i,d1)-trans1(d1,d1))>1.d-4.and.abs(cell1copy(i,d2)-trans1(d2,d2))>1.d-4)ncell1=ncell1+1
 enddo
 write(*,*) 'CELL2 SCAN'
 do i=1,size(cell2copy,1)
  if(abs(cell2copy(i,d1)-trans2(d1,d1))>1.d-4.and.abs(cell2copy(i,d2)-trans2(d2,d2))>1.d-4)ncell2=ncell2+1
 enddo
 write(*,*) 'CORRECTED NCELL1 NCELL2 : ', ncell1,ncell2
 deallocate(cell1,cell2,labelcell1,labelcell2,symcell1,symcell2)
 allocate(cell1(ncell1,3),cell2(ncell2,3))
 allocate(labelcell1(ncell1),labelcell2(ncell2),symcell1(ncell1),symcell2(ncell2))
 
 j=0
 do i=1,size(cell1copy,1)
  if(abs(cell1copy(i,d1)-trans1(d1,d1))>1.d-4.and.abs(cell1copy(i,d2)-trans1(d2,d2))>1.d-4)then
    j=j+1
    cell1(j,:)=cell1copy(i,:)
    labelcell1(j)=labelcell1copy(i)
    symcell1(j)=symcell1copy(i)
    cell1(j,d1)=cell1(j,d1)-cell1shiftd1*trans1(d1,d1)
    cell1(j,d2)=cell1(j,d2)-cell1shiftd2*trans1(d2,d2)
    cell1(j,d3)=cell1(j,d3)-cell1shiftd3*trans1(d3,d3)
    if(cell1(j,d1)>cutoff*trans1(d1,d1))cell1(j,d1)=cell1(j,d1)-trans1(d1,d1)
    if(cell1(j,d2)>cutoff*trans1(d2,d2))cell1(j,d2)=cell1(j,d2)-trans1(d2,d2)
  endif
 enddo
 j=0
 do i=1,size(cell2copy,1)
  if(abs(cell2copy(i,d1)-trans2(d1,d1))>1.d-4.and.abs(cell2copy(i,d2)-trans2(d2,d2))>1.d-4)then
    j=j+1
    cell2(j,:)=cell2copy(i,:)
    labelcell2(j)=labelcell2copy(i)
    symcell2(j)=symcell2copy(i)
    cell2(j,d1)=cell2(j,d1)-cell2shiftd1*trans2(d1,d1)
    cell2(j,d2)=cell2(j,d2)-cell2shiftd2*trans2(d2,d2)
    cell2(j,d3)=cell2(j,d3)-cell2shiftd3*trans2(d3,d3)
    if(cell2(j,d1)>cutoff*trans2(d1,d1))cell2(j,d1)=cell2(j,d1)-trans2(d1,d1)
    if(cell2(j,d2)>cutoff*trans2(d2,d2))cell2(j,d2)=cell2(j,d2)-trans2(d2,d2)
  endif
 enddo
 deallocate(cell1copy,cell2copy,labelcell1copy,labelcell2copy,symcell1copy,symcell2copy)

 77 continue

 allocate(planecell1(ncell1),planecell2(ncell2))
 do i=1,ncell1
  planer1=((cell1(i,d3)-minval(cell1(:,d3)))/trans1(d3,d3)-planeshift)*dble(planes_in_cell1)
  if(round_down_planes)then
   planecell1(i)= INT( planer1 )+1
  else
   planecell1(i)=NINT( planer1 )+1
  endif
  if(labelcell1(i)==trim(adjustl(corsymbol1(1)))) write(*,*) 'CELL1, PLANE N: ', i,labelcell1(i),planecell1(i),symcell1(i)
 enddo
 do i=1,ncell2 
  planer2=((cell2(i,d3)-minval(cell2(:,d3)))/trans2(d3,d3)-planeshift)*dble(planes_in_cell2)
  if(round_down_planes)then
   planecell2(i)= INT(planer2)+1
  else
   planecell2(i)=NINT(planer2)+1
  endif
  if(labelcell2(i)==trim(adjustl(corsymbol2(1)))) write(*,*) 'CELL2, PLANE N: ', i,labelcell2(i),planecell2(i),symcell2(i)
 enddo

 natom =  nplane1*width**2*ncell1  +  nplane2*width**2*ncell2

 tota=0.;totb=0.;totc=0.
 tota(d1)  = width  *trans1(d1,d1)
 totb(d2)  = width  *trans1(d2,d2)
 totc(d3)  = nplane1*trans1(d3,d3) + nplane2*trans2(d3,d3)  
 write(*,*) 'tota-totb CELL1 : ', tota(d1),totb(d2)
 tota(d1)  = width  *trans2(d1,d1)
 totb(d2)  = width  *trans2(d2,d2)
 write(*,*) 'tota-totb CELL2 : ', tota(d1),totb(d2)
 write(*,*) 'totc            : ', totc(d3)

 allocate(coordfull(natom,3),labelfull(natom),planefull(natom,nmaxentry),nentries(natom),symfull(natom),belongscell(natom),planeloc(natom))
 planefull=0;symfull=0;belongscell=0;planeloc=0

 write(*,*) 'PLEASE FOLLOW THIS PROCEDURE:'
 write(*,*) 'open cif file with VESTA, export as pdb'
 write(*,*) 'use cif2struct to obtain struct file'
 write(*,*) 'build filename.input, where the first line are the a,b,c numbers obtained from the struct file'
 write(*,*) 'then follows the atom label, coordinates, obtained from pdb file'

 write(*,*) 'CELL1'
 do i=1,ncell1
 write(*,*) labelcell1(i),(cell1(i,ii),ii=1,3),symcell1(i)
 enddo
 write(*,*) 'CELL2'
 do i=1,ncell2
 write(*,*) labelcell2(i),(cell2(i,ii),ii=1,3),symcell2(i)
 enddo
 

 !####################################################!
 !####################################################!
 !####################################################!

 hubcounter=0

 !####################################################!
 !####################################################!
 !####################################################!
 !WRITE STRUCTURE 1 COORDINATES

  shift=0.0; shift1_sand=0.0

  do k=1,nplane1
  if(sandwitch==1)then
   if(k<=nplane1/2)then
    shift1_sand = -dble(nplane1/2)*trans1(d3,:) + dble(nplane1)*trans1(d3,:) &
                & +(dble(nplane2)+separation)*trans2(d3,:)
   else
    shift1_sand = -dble(nplane1/2)*trans1(d3,:)
   endif
  elseif(sandwitch==2)then
    shift1_sand =   dble(nplane2/2)*trans2(d3,:)
  endif
  do i=1,width
  do j=1,width
  do at=1,ncell1
  hubcounter=hubcounter+1
  coordfull(hubcounter,:)=shift1_sand+shift + (dble(k-1))*trans1(d3,:) + &  
                                            & (dble(i-1))*trans1(d1,:) + &  
                                            & (dble(j-1))*trans1(d2,:) + cell1(at,:)
  labelfull(hubcounter)=trim(labelcell1(at))
  planefull(hubcounter,1)=(k-1) *planes_in_cell1 + planecell1(at)
  planeloc(hubcounter)=planecell1(at)
  belongscell(hubcounter)=1
  symfull(hubcounter)=symcell1(at)
  enddo
  enddo
  enddo
  enddo

 !####################################################!
 !####################################################!
 !####################################################!
 !WRITE STRUCTURE 2 COORDINATES

  if(sandwitch==1)then
    shift1_sand = -dble(nplane1/2)*trans1(d3,:)
  elseif(sandwitch==2)then
    shift1_sand =  dble(nplane2/2)*trans2(d3,:)
  else
    shift1_sand =  0.0
  endif

  shift = shift1_sand + dble(nplane1)*trans1(d3,:)

  shift2_sand=0.0

  shift(d3) = shift(d3) + separation * trans2(d3,d3)   
  totc(d3)  = dble(nplane1)*trans1(d3,d3) + (dble(nplane2)+separation)*trans2(d3,d3)


  do k=1,nplane2
  if(sandwitch==2.and.k>nplane2/2)then
   shift2_sand = - dble(nplane1)*trans1(d3,:) - (dble(nplane2)+separation)*trans2(d3,:)
  else
   shift2_sand = 0.0
  endif
  do i=1,width
  do j=1,width
  do at=1,ncell2
  hubcounter=hubcounter+1
  coordfull(hubcounter,:)=shift2_sand+shift+(dble(k-1))*trans2(d3,:)+ &
                                          & (dble(i-1))*trans2(d1,:)+ & 
                                          & (dble(j-1))*trans2(d2,:)+cell2(at,:)
  labelfull(hubcounter)=trim(labelcell2(at))
  symfull(hubcounter)=symcell2(at)
  planefull(hubcounter,1)=(k-1)*planes_in_cell2 + nplane1*planes_in_cell1 + planecell2(at)
  planeloc(hubcounter)=planecell2(at)
  belongscell(hubcounter)=2
  enddo
  enddo
  enddo
  enddo


 !####################################################!
 !####################################################!
 !####################################################!
 ! CLEAN STRUCTURE FROM DUPLICATE

  natom=0
  
  do i=1,size(coordfull,1)
   duplicate=.false.
   do j=1,i-1
    if(norm(coordfull(i,:)               -coordfull(j,:))<same_atom_distance .or. &
      &norm(coordfull(i,:)-totc(d3)*vecd3-coordfull(j,:))<same_atom_distance .or. &
      &norm(coordfull(i,:)+totc(d3)*vecd3-coordfull(j,:))<same_atom_distance &
      & )then
      duplicate=.true.
      exit
    endif
   enddo
   if(.not.duplicate) natom=natom+1
  enddo
  write(*,*) 'NUMBER OF ATOMS (CLEANED) : ', natom

  if(natom>size(coordfull,1)) then
   write(*,*) 'ERROR CLEANING'
   stop
  endif

  allocate(coordfullcopy(size(coordfull,1),3),labelfullcopy(size(labelfull)),planefullcopy(size(planefull,1),size(planefull,2)),planeloccopy(size(planeloc)))
  allocate(belongscellcopy(size(belongscell)))
  allocate(symfullcopy(size(symfull)))

  planefullcopy=0
  planeloccopy=0
  belongscellcopy=0

  coordfullcopy=coordfull
  labelfullcopy=labelfull
  planefullcopy=planefull
  planeloccopy=planeloc
  symfullcopy=symfull
  belongscellcopy=belongscell

  deallocate(coordfull,labelfull,planefull,symfull,belongscell,planeloc)
  allocate(coordfull(natom,3),labelfull(natom),planefull(natom,nmaxentry),symfull(natom),belongscell(natom),planeloc(natom))
  planefull=0
  planeloc=0 
  natom=0 
  do i=1,size(coordfullcopy,1)
   duplicate=.false.

   do j=1,i-1
    if(norm(coordfullcopy(i,:)               -coordfullcopy(j,:))<same_atom_distance .or. &
      &norm(coordfullcopy(i,:)-totc(d3)*vecd3-coordfullcopy(j,:))<same_atom_distance .or. &
      &norm(coordfullcopy(i,:)+totc(d3)*vecd3-coordfullcopy(j,:))<same_atom_distance &
      & )then
      duplicate=.true.
      if(symfullcopy(i)/=symfullcopy(j))then
        write(*,*) 'WARNING : duplicate but different SYMMETRY !!'
        write(*,*) 'They belong to cell : ', belongscellcopy(i),belongscellcopy(j)
        write(*,*) 'Labels of atoms     : ', labelfullcopy(i),labelfullcopy(j)
        write(*,*) 'Position of atoms, i: ', coordfullcopy(i,:)
        write(*,*) 'Position of atoms, j: ', coordfullcopy(j,:)
        write(*,*) 'duplicatecellchoice : ', duplicatecellchoice

        if(duplicatecellchoice/=0)then
           if(planeloccopy(i)<planeloccopy(j))then  
           !choosing i
           !------------->
           do jjj=1,natom
             if(norm(coordfull(jjj,:)               -coordfullcopy(i,:))<same_atom_distance .or. &
               &norm(coordfull(jjj,:)-totc(d3)*vecd3-coordfullcopy(i,:))<same_atom_distance .or. &
               &norm(coordfull(jjj,:)+totc(d3)*vecd3-coordfullcopy(i,:))<same_atom_distance &
               & )then
                   write(*,*) 'SWAPPING ELEMENT'
                   coordfull(jjj,:)=coordfullcopy(i,:)
                   planefull(jjj,:)=planefullcopy(i,:)
                   planeloc(jjj)   =planeloccopy(i)
                   labelfull(jjj  )=labelfullcopy(i)
                     symfull(jjj  )=symfullcopy(i)
                 belongscell(jjj)  =belongscell(i)
                  
                 if(cutoff_boundary_shifted>1.d-4)Then 
                    if(abs(coordfullcopy(i,d3)-totc(d3))<abs(cutoff_boundary_shifted))then
                       coordfull(jjj,:)=coordfullcopy(i,:)-totc(d3)*vecd3
                       planefull(jjj,1)=firstplane
                    endif
                 elseif(cutoff_boundary_shifted<-1.d-4)then
                    if(abs(coordfullcopy(i,d3))<abs(cutoff_boundary_shifted))then
                       coordfull(jjj,:)=coordfullcopy(i,:)+totc(d3)*vecd3
                       planefull(jjj,1)=lastplane
                    endif
                 endif

                 exit
             endif 
           enddo
           if(j==natom+1)Then
            write(*,*) 'ERROR site not found...'
            stop
           endif
          else
           !choosing j
           !-------------> 
          endif
        endif

      endif
      exit
    endif
   enddo

   if(.not.duplicate)then
      natom=natom+1
      jj=i
      coordfull(natom,:)=coordfullcopy(jj,:)
      planefull(natom,:)=planefullcopy(jj,:)
      planeloc(natom)   =planeloccopy(jj)
      labelfull(natom  )=labelfullcopy(jj)
        symfull(natom  )=symfullcopy(jj)
      belongscell(natom)=belongscell(jj)

     if(cutoff_boundary_shifted>1.d-4)then  
      if(abs(coordfullcopy(jj,d3)-totc(d3))<abs(cutoff_boundary_shifted))then
        coordfull(natom,:)=coordfullcopy(jj,:)-totc(d3)*vecd3
        write(*,*) 'WARNING      : ATOMS ON BOUNDARY - SHIFT THEM TO LOWER COORDINATES'
        write(*,*) 'PLANE WAS    : ', planefull(natom,1)
        planefull(natom,1)=firstplane
        write(*,*) 'PLANE IS NOW : ', planefull(natom,1)
      endif
     elseif(cutoff_boundary_shifted<-1.d-4)then
      if(abs(coordfullcopy(jj,d3))<abs(cutoff_boundary_shifted))then
        coordfull(natom,:)=coordfullcopy(jj,:)+totc(d3)*vecd3
        write(*,*) 'WARNING      : ATOMS ON BOUNDARY - SHIFT THEM TO LOWER COORDINATES'
        write(*,*) 'PLANE WAS    : ', planefull(natom,1)
        planefull(natom,1)=lastplane
        write(*,*) 'PLANE IS NOW : ', planefull(natom,1)
      endif
     endif
   endif

  enddo

  deallocate(coordfullcopy,labelfullcopy,planefullcopy,symfullcopy,belongscellcopy,planeloccopy)

  !--------------------------------------------------------------------!
   if(symmetrize_structure)then 
    do jj=1,natom
       if(abs(coordfull(jj,d3)-totc(d3))<abs(cutoff_boundary_shifted).and.coordfull(jj,d1)>0.0)then
         coordfull(jj,:)=coordfull(jj,:)-totc(d3)*vecd3
         write(*,*) 'WARNING      : ATOMS ON BOUNDARY - SHIFT THEM TO LOWER COORDINATES'
         write(*,*) 'PLANE WAS    : ', planefull(jj,1)
         planefull(jj,1)=firstplane

         if(onlycorplanes>0.) planefull(jj,1)=lastplane+20000
 
        write(*,*) 'PLANE IS NOW : ', planefull(jj,1)
       endif
       if(abs(coordfull(jj,d3)) <abs(cutoff_boundary_shifted).and.coordfull(jj,d1)<=0.0)then
         coordfull(jj,:)=coordfull(jj,:)+totc(d3)*vecd3
         write(*,*) 'WARNING      : ATOMS ON BOUNDARY - SHIFT THEM TO LOWER COORDINATES'
         write(*,*) 'PLANE WAS    : ', planefull(jj,1)
         planefull(jj,1)=lastplane

         if(onlycorplanes>0.) planefull(jj,1)=lastplane+20000

         write(*,*) 'PLANE IS NOW : ', planefull(jj,1)
       endif
    enddo
   endif

   nentries=1

   if(onlycorplanes>0.)then
      
    do jj=1,natom
     if(labelfull(jj)==trim(adjustl(corsymbol1(1))).or.labelfull(jj)==trim(adjustl(corsymbol2(1))))then
      !......... 
     else
      planefull(jj,:)=0
      nentries(jj)=0
     endif
    enddo

    do jj=1,natom
     if(labelfull(jj)==trim(adjustl(corsymbol1(1))).or.labelfull(jj)==trim(adjustl(corsymbol2(1))))then
       do jjj=1,natom
        if(jjj/=jj)then
        if(norm(coordfull(jj,:) -coordfull(jjj,:))< onlycorplanes )then
          nentries(jjj)=nentries(jjj)+1
          planefull(jjj,nentries(jjj))=planefull(jj,1)
        endif

         if(periodic_plane)then
          if(norm(coordfull(jj,:) -coordfull(jjj,:)-totc)< onlycorplanes )then
           nentries(jjj)=nentries(jjj)+1
           planefull(jjj,nentries(jjj))=planefull(jj,1)
          endif
          if(norm(coordfull(jj,:) -coordfull(jjj,:)+totc)< onlycorplanes )then
           nentries(jjj)=nentries(jjj)+1
           planefull(jjj,nentries(jjj))=planefull(jj,1)
          endif
         endif

        endif
       enddo
      endif
    enddo
   endif

  !--------------------------------------------------------------------!

  write(*,*) 'WRITING COORD FILE'

  open(unit=100,file='coord.xyz')
  write(100,*) natom
  do i=1,nbasis
    write(100,'(a)',advance='NO') trim(adjustl(atomsymbol(i)))//' '
  enddo
  write(100,*)
  do i=1,natom
     write(100,'(a,3f17.6)') trim(Adjustl(labelfull(i))),(coordfull(i,ii),ii=1,3)
  enddo
  close(100)

  open(unit=100,file='coord_zperio.xyz')
  write(100,*) natom*3
  do i=1,nbasis
    write(100,'(a)',advance='NO') trim(atomsymbol(i))//' '
  enddo
  write(100,*)
  do i=1,natom
     write(100,'(a,3f17.8)') trim(Adjustl(labelfull(i))),(coordfull(i,ii)                   ,ii=1,3)
  enddo
  do i=1,natom
     write(100,'(a,3f17.8)') trim(Adjustl(labelfull(i))),(coordfull(i,ii)-totc(d3)*vecd3(ii),ii=1,3)
  enddo
  do i=1,natom
     write(100,'(a,3f17.8)') trim(Adjustl(labelfull(i))),(coordfull(i,ii)+totc(d3)*vecd3(ii),ii=1,3)
  enddo
  close(100)

 write(*,*) 'BUIDING DMFT INPUTS'

 !####################################################!
 !####################################################!
 !####################################################!

 do i=1,natom
 write(*,*) 'ATOM / SYM : ', i, symfull(i)
  if(labelfull(i)==trim(adjustl(corsymbol1(1))).or.labelfull(i)==trim(adjustl(corsymbol2(1)))) then
      write(114,*) labelfull(i),coordfull(i,d3),planefull(i,1)
      if(nentries(i)/=1)then
        write(*,*) 'ERROR : a correlated atom can NOT belong to several planes'
        write(*,*) ' the planes are used to build mask_uniform in particular'
        write(*,*) 'atom and entries : ', labelfull(i),nentries(i)
        stop
      endif
  endif
 enddo

 hubcounter=0
 do i=1,natom
  write(102,'(a,100f10.3)') trim(labelfull(i)),coordfull(i,:)
 enddo

 hubcounter=0
 do i=1,natom
  do j=1,ncor
   if(labelfull(i)==trim(adjustl(corsymbol1(j))).or.labelfull(i)==trim(adjustl(corsymbol2(j))))then
    hubcounter=hubcounter+1
   endif
  enddo
 enddo

 write(103,*) 'THERE ARE [X] COR ATOMS : ', hubcounter
 allocate(myatoms(hubcounter),hubequivalent(hubcounter),isitcor(natom))
 isitcor=.false.
 hubcounter=0
 do i=1,natom
  do j=1,ncor
   if(labelfull(i)==trim(adjustl(corsymbol1(j))).or.labelfull(i)==trim(adjustl(corsymbol2(j))))then
    hubcounter=hubcounter+1
    myatoms(hubcounter)=i
    isitcor(i)=.true.
    exit
   endif
  enddo
 enddo
 do i=1,hubcounter
  write(103,*) myatoms(i),labelfull(myatoms(i)),coordfull(myatoms(i),:),planefull(myatoms(i),1)
 enddo

 open(unit=104,file='hub.xyz')
 write(104,*) hubcounter
 write(104,*) (trim(adjustl(corsymbol1(i)))//" "//trim(adjustl(corsymbol2(i))),i=1,ncor)
 do i=1,hubcounter
  write(104,'(a,3f17.8)') labelfull(myatoms(i)),coordfull(myatoms(i),:) 
 enddo
 close(104)

 if(brokensym)then
  open(unit=104,file='hubUP.xyz')
  j=0
  do i=1,hubcounter
   if(symfull(myatoms(i))>0) j=j+1
  enddo
  write(104,*) j
  write(104,*) (trim(adjustl(corsymbol1(i)))//" "//trim(adjustl(corsymbol2(i))),i=1,ncor)
  do i=1,hubcounter
   if(symfull(myatoms(i))>0) write(104,'(a,3f17.8)') labelfull(myatoms(i)),coordfull(myatoms(i),:)
  enddo
  close(104)
  open(unit=104,file='hubDO.xyz')
  j=0
  do i=1,hubcounter
   if(symfull(myatoms(i))<0) j=j+1
  enddo
  write(104,*) j
  write(104,*) (trim(adjustl(corsymbol1(i)))//" "//trim(adjustl(corsymbol2(i))),i=1,ncor)
  do i=1,hubcounter
   if(symfull(myatoms(i))<0) write(104,'(a,3f17.8)') labelfull(myatoms(i)),coordfull(myatoms(i),:)
  enddo
  close(104)

  open(unit=104,file='hubUPDO.xyz')
  write(104,*) hubcounter
  write(104,*) ("u"//trim(adjustl(corsymbol1(i)))//" u"//trim(adjustl(corsymbol2(i))),i=1,ncor),(" d"//trim(adjustl(corsymbol1(i)))//" d"//trim(adjustl(corsymbol2(i))),i=1,ncor)
  do i=1,hubcounter
   if(symfull(myatoms(i))<0) write(104,'(a,3f17.8)') "d"//trim(adjustl(labelfull(myatoms(i)))),coordfull(myatoms(i),:)
   if(symfull(myatoms(i))>0) write(104,'(a,3f17.8)') "u"//trim(adjustl(labelfull(myatoms(i)))),coordfull(myatoms(i),:)
  enddo
  close(104)

 endif


 write(*,*) 'BUILDING EQUIVALENT SELF ENERGIES PATTERN'

 hubequivalent=0

 do i=1,hubcounter
  hubequivalent(i)=i
  do j=1,i-1
   if( planefull(myatoms(i),1)==planefull(myatoms(j),1) )then
     hubequivalent(i)=j
     exit
   endif
  enddo
 enddo 

 write(*,*) 'BULDING MASK_UNIFORM MASKS'

 if(nplane1==0.or.nplane2==0)then
  if(brokensym)then
   do i=1,hubcounter
    if(symfull(myatoms(i))>0)then
      exit
    endif
   enddo
   if(i==hubcounter+1)then
    write(*,*) 'ERROR : symmetry - AF pattern, could not find candidate'
    stop
   endif
      hubequivalent=hubequivalent(i)
  else
      hubequivalent=hubequivalent(1)
  endif
 endif

 if(brokensym)then
  do i=1,hubcounter
    write(*,*) 'HUBBARD SITE / SYM : ', i,symfull(myatoms(i))
    hubequivalent(i) = symfull(myatoms(i))*hubequivalent(i)
  enddo
  do i=1,hubcounter
   if(hubequivalent(i)==-i)then
    !locate other candidate, which has + sign
    do j=1,hubcounter
     if(hubequivalent(j)==i)then
      exit
     endif
    enddo
    if(j==hubcounter+1)then
     write(*,*) 'ERROR SYMMETRY : DIDNT FIND A CANDIDATE FOR SYM(I)/=-I'
     stop
    endif
    !replace with this candidate
    where(hubequivalent== i) hubequivalent= j
    where(hubequivalent==-i) hubequivalent=-j
   endif
  enddo 
 endif

 open(unit=105,file='mask_uniform')
 open(unit=106,file='mask_uniform_coord')
 do i=1,hubcounter
  write(105,*)  i,hubequivalent(i)
  write(106,*)  coordfull(myatoms(i),:)    
  write(106,*)  coordfull(myatoms(abs(hubequivalent(i))),:)
 enddo
 close(105)
 close(106)

 open(unit=105,file='mask_dimer')
 do i=1,hubcounter
  write(105,*) i,'0 0 ' 
 enddo
 close(105)

 open(unit=105,file='mask_projections')
 write(105,*) 'T T T T T'
 close(105)

 open(unit=100,file='rotation.input')
 write(100,*) 'coord     !   positions of orginal cell'
 write(100,'(i5,a)')              natom , ' !   natom_small'
 write(100,*) trim(adjustl(corsymbol1(1))), ' !   ref'
 write(100,'(f7.3,a)')  cutoff_rot        , ' !   cutoff plane to determine which Oxygen are apical'
 vecd3=0.0
 vecd3(rutile_dir)=1.0
 write(100,'(3f7.3,a)') (vecd3(i),i=1,3), ' ! rutile axis - apical-cor is perpendicular to rutile'
 write(100,'(f7.3,a)') dist_neigh       , ' !   max_dist for bonds'
 write(100,*) 'temp      !   additional coordinate file'
 write(100,*) '1.0       !   flip  '
 write(100,*) '.true.    !   keep perp to rutile'
 write(100,*) '.true.    !   all planes include all 4-sites to define the 4-site plane'
 close(100)


 open(unit=100,file='axis_file')
 write(100,*) (myaxis(j,1),j=1,3) 
 write(100,*) (myaxis(j,2),j=1,3) 
 write(100,*) (myaxis(j,3),j=1,3)
 close(100)

 open(unit=100,file='coord') 
 write(100,'(3f17.8)') (tota(i),i=1,3) 
 write(100,'(3f17.8)') (totb(i),i=1,3)
 write(100,'(3f17.8)') (totc(i),i=1,3)
 do i=1,natom
   if(isitcor(i))then
    write(100,'(a,3f17.8)') trim(adjustl(corsymbol1(1))),(coordfull(i,j),j=1,3)
   else
    write(100,'(a,3f17.8)') trim(adjustl(labelfull(i))),(coordfull(i,j),j=1,3)
   endif
 enddo
 close(100)

 call system(" onetep.dmft.local.rotation ")

 open(unit=6121,file='stitch_coord')
 write(6121,*) nbasis
 write(6121,'(100a)')   (trim(adjustl(atomsymbol(i)))//" ",i=1,nbasis)
 write(6121,         *) (basisZ(i),i=1,nbasis)
 write(6121,'(3f17.8)') (tota(i),i=1,3)
 write(6121,'(3f17.8)') (totb(i),i=1,3)
 write(6121,'(3f17.8)') (totc(i),i=1,3)
 do i=1,natom
  write(6121,'(a,3f17.8)') trim(adjustl(labelfull(i))),(coordfull(i,j),j=1,3)
 enddo
 close(6121)

 if(nplane1==0.or.nplane2==0)then
 if(nplane1==0) fac=1.0/dble(nplane2)
 if(nplane2==0) fac=1.0/dble(nplane1)
 open(unit=6121,file='stitch_coord_shifted_left')
 write(6121,*) nbasis
 write(6121,'(100a)')   (trim(adjustl(atomsymbol(i)))//" ",i=1,nbasis)
 write(6121,         *) (basisZ(i),i=1,nbasis)
 write(6121,'(3f17.8)') (tota(i),i=1,3)
 write(6121,'(3f17.8)') (totb(i),i=1,3)
 write(6121,'(3f17.8)') (totc(i),i=1,3)
 do i=1,natom
  write(6121,'(a,3f17.8)') trim(adjustl(labelfull(i))),(coordfull(i,j)-totc(j)*fac,j=1,3) 
 enddo
 close(6121)
 open(unit=6121,file='stitch_coord_shifted_right')
 write(6121,*) nbasis
 write(6121,'(100a)')   (trim(adjustl(atomsymbol(i)))//" ",i=1,nbasis)
 write(6121,         *) (basisZ(i),i=1,nbasis)
 write(6121,'(3f17.8)') (tota(i),i=1,3)
 write(6121,'(3f17.8)') (totb(i),i=1,3)
 write(6121,'(3f17.8)') (totc(i),i=1,3)
 do i=1,natom
  write(6121,'(a,3f17.8)') trim(adjustl(labelfull(i))),(coordfull(i,j)+totc(j)*fac,j=1,3)
 enddo
 close(6121)
 endif

 call write_case_file

  i=nplane1*planes_in_cell1 + nplane2*planes_in_cell2

  if(sandwitch==0)then
  where(planefull==i+1) planefull=1
 elseif(sandwitch==1)then
  where(planefull==i+1) planefull=1
 elseif(sandwitch==2)then
  where(planefull==i+1) planefull=1
 elseif(sandwitch>2)then
  write(*,*) 'sandwitch > 2 not yet implemented'
  stop
 endif

 if(abs(onlycorplanes)>0.001)then
  if(minval(planefull)/=0)then
   write(*,*) 'WARNING : min values of planefull expected to be zero'
  endif
 endif

 k=1
 do i=1,natom
   j=minval(planefull)
   where(planefull==j) planefull=1000000-k
   if(minval(planefull)>1000000-50000)then
    exit
   endif
   k=k+1
 enddo

 write(*,*) 'THERE ARE [x] PLANES : ', k-1  ! zero removed, hence the -1
 obtainedplanes=k-1
 do i=1,k
  where(planefull==1000000-i) planefull=i-1
 enddo
 i=nplane1*planes_in_cell1 + nplane2*planes_in_cell2
 if(symmetrize_structure)i=i+1
 write(*,*) 'SHOULD BE : ', i
 if(i/=k-1)Then
   write(*,*) 'ERROR : number of planes do not match'
   stop
 endif

 open(unit=5050,file='group_of_atoms_atom_labels')
 do i=1,natom
  do j=1,nentries(i) 
   write(5050,'(a,i7)') trim(adjustl(labelfull(i))),planefull(i,j)
  enddo
 enddo
 close(5050)


 do kkk=1,obtainedplanes
  open(unit=5050,file='group_of_atoms_atom_labels_'//trim(adjustl(tostring(kkk)))//".xyz")

  jj=0
  do i=1,natom
   do j=1,nentries(i)
     if(planefull(i,j)==kkk) jj=jj+1 
   enddo
  enddo
  write(5050,*) jj 

  do i=1,nbasis
    write(5050,'(a)',advance='NO') trim(adjustl(atomsymbol(i)))//' '
  enddo
  write(5050,*)
  do i=1,natom
  do j=1,nentries(i)
     if(planefull(i,j)==kkk) write(5050,'(a,3f17.6)') trim(Adjustl(labelfull(i))),(coordfull(i,ii),ii=1,3)
  enddo
  enddo
  close(5050)
 enddo


 open(unit=5050,file='group_of_atoms_label')
 write(5050,*) obtainedplanes 
 do i=1,natom
  do j=1,nentries(i)
    write(5050,'(i10,i10)') i,planefull(i,j)
  enddo
 enddo
 close(5050)

 i=obtainedplanes
 write(*,*) 'CHECKING HOW MANY ATOMS IN EACH LAYER'
 call system(" for i in `seq 1 "//trim(adjustl(tostring(i))) //"` ; do grep $i group_of_atoms_atom_labels | wc -l ; done ; ")

 open(unit=7000,file='mask_u')
 open(unit=7001,file='mask_j')

 do j=1,ncor
  do i=1,natom
    if(labelfull(i)==trim(adjustl(corsymbol1(j))))then
      write(7000,*)    'UU='//trim(adjustl(tostring(corU1( j  )))) 
      write(7001,*) 'Jhund='//trim(adjustl(tostring(corJ1( j  ))))
    endif
  enddo
 enddo
 do j=1,ncor
  do i=1,natom
   if(trim(adjustl(corsymbol1(j)))/=trim(adjustl(corsymbol2(j))))then
    if(labelfull(i)==trim(adjustl(corsymbol2(j))))then
      write(7000,*)    'UU='//trim(adjustl(tostring(corU2( j  ))))
      write(7001,*) 'Jhund='//trim(adjustl(tostring(corJ2( j  ))))
    endif
   endif
  enddo
 enddo

 close(7000) 
 close(7001)

contains

 real(8) function norm(x)
 real(8) x(3)
  norm=dsqrt(sum(x**2))
 end function

 subroutine write_case_file

 open(unit=500,file='run.dat')  

 if(.not.dmft_case)then

 write(500,'(a)')'task                  : SINGLEPOINT'
 write(500,'(a)')'kernel_cutoff         : 4000.000000'
 write(500,'(a)')'cutoff_energy         : 850.000000  eV'
 write(500,'(a)')'fine_grid_scale       : 4.0'
 write(500,'(a)')'xc_functional         : PW92'
 write(500,'(a)')'do_properties         : T'
 write(500,'(a)')'homo_plot             : 8'
 write(500,'(a)')'lumo_plot             : 8'
 write(500,'(a)')'dos_smear             : 0.1 eV'
 write(500,'(a)')'ldos_smear            : 0.1 eV'
 write(500,'(a)')'num_eigenvalues       : 400'
 write(500,'(a)')'write_density_plot    : T'
 write(500,'(a)')'write_denskern        : T'
 write(500,'(a)')'write_tightbox_ngwfs  : T'
 write(500,'(a)')'read_denskern         : F'
 write(500,'(a)')'read_tightbox_ngwfs   : F'
 write(500,'(a)')'write_forces          : T'
 write(500,'(a)')'write_ngwf_plot       : F'
 write(500,'(a)')'write_ngwf_grad_plot  : T'
 write(500,'(a)')'cube_format           : T'
 write(500,'(a)')'kzero                 : 3.0'
 write(500,'(a)')'occmix                : 0.25'
 write(500,'(a)')'elec_cg_max           : 2'
 write(500,'(a)')'minit_lnv             : 20'
 write(500,'(a)')'maxit_lnv             : 8'
 write(500,'(a)')'maxit_pen             : 0'
 write(500,'(a)')'maxit_ngwf_cg         : 40'
 write(500,'(a)')'ngwf_threshold_orig   : 1.6e-6'
 write(500,'(a)')'elec_energy_tol       : 2e-7 Hartree'
 write(500,'(a)')'usespacefillingcurve  : F'
 write(500,'(a)')'initial_dens_realspace: T'
 write(500,'(a)')'maxit_palser_mano     : 400'
 write(500,'(a)')'output_detail         : VERBOSE'
 write(500,'(a)')'comms_group_size      : 1'
 write(500,'(a)')'lnv_check_trial_steps : T'
 write(500,'(a)')'paw                   : T'
 write(500,'(a)')'#devel_code            : NGWF_FD'
 
 else

 write(500,'(a)')'task                  : SINGLEPOINT'
 write(500,'(a)')'kernel_cutoff         : 4000.000000'
 write(500,'(a)')'cutoff_energy         : 850.000000  eV'
 write(500,'(a)')'fine_grid_scale       : 4.0'
 write(500,'(a)')'xc_functional         : PW92'
 write(500,'(a)')'do_properties         : F'
 write(500,'(a)')'homo_plot             : -1'
 write(500,'(a)')'lumo_plot             : -1'
 write(500,'(a)')'dos_smear             : 0.0 eV'
 write(500,'(a)')'ldos_smear            : 0.0 eV'
 write(500,'(a)')'num_eigenvalues       : 0'
 write(500,'(a)')'write_density_plot    : F'
 write(500,'(a)')'write_denskern        : T'
 write(500,'(a)')'write_tightbox_ngwfs  : T'
 write(500,'(a)')'read_denskern         : F'
 write(500,'(a)')'read_tightbox_ngwfs   : F'
 write(500,'(a)')'write_forces          : F'
 write(500,'(a)')'write_ngwf_plot       : F'
 write(500,'(a)')'write_ngwf_grad_plot  : F'
 write(500,'(a)')'cube_format           : F'
 write(500,'(a)')'kzero                 : 3.0'
 write(500,'(a)')'occmix                : 0.25'
 write(500,'(a)')'elec_cg_max           : 2'
 write(500,'(a)')'minit_lnv             : 2'
 write(500,'(a)')'maxit_lnv             : 8'
 write(500,'(a)')'maxit_pen             : 0'
 write(500,'(a)')'maxit_ngwf_cg         : 3'
 write(500,'(a)')'ngwf_threshold_orig   : 1.6e-6'
 write(500,'(a)')'elec_energy_tol       : 2e-7 Hartree'
 write(500,'(a)')'usespacefillingcurve  : F'
 write(500,'(a)')'initial_dens_realspace: T'
 write(500,'(a)')'output_detail         : VERBOSE'
 write(500,'(a)')'comms_group_size      : 1'
 write(500,'(a)')'lnv_check_trial_steps : T'
 write(500,'(a)')'paw                   : T'
 write(500,'(a)')'#devel_code            : NGWF_FD'

 write(500,*) 'spin_polarized        :  F'
 write(500,*) 'spin                  :  0'
 write(500,*) 'ngwf_analysis         :  F'
 write(500,*) 'optics_energy_max     :  0.0 eV'
 write(500,*) 'dmft_order_proj       :  0.0'
 write(500,*) 'dmft_switch_off_proj_order : F'
 write(500,*) 'dmft_nmu_loop         :  7'
 write(500,*) 'dmft_emin             : -1.0 Hartree'
 write(500,*) 'dmft_emax             :  1.0 Hartree'
 write(500,*) 'dmft_points           :  '//trim(adjustl(tostring(nfrequ)))
 write(500,*) 'dmft_smear            :  0.0005   Hartree'
 write(500,*) 'dmft_smear_shift      :  0.00     Hartree'
 write(500,*) 'dmft_smear_T          :  0.008    Hartree'
 write(500,*) 'dmft_smear_eta        :  0.008    Hartree'
 write(500,*) 'dmft_smear_w          :  0.045    Hartree'
 write(500,*) 'dmft_paramagnetic     :  T'
 write(500,*) 'dmft_cutoff_tail      :   1.00  Hartree'
 write(500,*) 'dmft_cutoff_small     :   0.005 Hartree'
 write(500,*) 'dmft_rotate_green     :  T'
 write(500,*) 'dmft_temp             :  0.00093 Hartree'
 write(500,*) 'dmft_optics           :  F'
 write(500,*) 'dmft_optics_i1        :  2'
 write(500,*) 'dmft_optics_i2        :  2'
 write(500,*) 'dmft_optics_window    :  0.22 Hartree'
 write(500,*) 'dmft_dos_min          : -0.23 Hartree'
 write(500,*) 'dmft_dos_max          :  0.23 Hartree'
 write(500,*) 'output_detail         :  VERBOSE'
 write(500,*) 'timings_level         :  0'
 write(500,*) 'print_qc              :  TRUE'
 write(500,*) 'kerfix                :  1'
 write(500,*) 'pen_param             :  20'
 write(500,*) 'maxit_palser_mano     : -1'
 write(500,*) 'maxit_kernel_fix      :  20'
 write(500,*) 'maxit_hotelling       :  0'
 write(500,*) 'delta_e_conv          :  T'
 write(500,*) 'write_converged_dkngwfs : F'
 
 endif

 write(500,*) '%block hubbard'

 do i=1,ncor
 do j=1,natom
 if( trim(adjustl(corsymbol1(i))) == trim(adjustl(labelfull(j)))) goto 55 
 enddo
 cycle
 55 continue
 write(500,*) trim(adjustl(corsymbol1(i))), ' 2  0.00 -10. 0.0 0.0'
 enddo

 do i=1,ncor
 do j=1,natom
 if( trim(adjustl(corsymbol2(i))) == trim(adjustl(labelfull(j)))) goto 65
 enddo
 cycle
 65 continue
 if(trim(adjustl(corsymbol2(i)))/=trim(adjustl(corsymbol1(i))))then
 write(500,*) trim(adjustl(corsymbol2(i))), ' 2  0.00 -10. 0.0 0.0'
 endif
 enddo

 write(500,*) '%endblock hubbard'

 write(500,*) 'hubbard_proj_mixing : 0.0'

 write(500,*)
 write(500,*)
 write(500,'(a)')'%block species'
 write(500,'(a)')'Ang'
 do i=1,nbasis
 do j=1,natom
 if( trim(adjustl(atomsymbol(i))) == trim(adjustl(labelfull(j)))) goto 56
 enddo
 cycle
 56 continue
  write(500,*) trim(adjustl(atomsymbol(i)))//' '//trim(adjustl(atomsymbol(i)))//' ',basistot(i),basisZ(i),trim(adjustl(tostring(radius_atom))) 
 enddo
 write(500,'(a)')'%endblock species'
 write(500,*)

 totc(d3)=totc(d3)+open_boundary

 write(500,'(a)')'%block lattice_cart'
 write(500,'(a)')'Ang'
 write(500,'(3f20.14)') (tota(i),i=1,3)
 write(500,'(3f20.14)') (totb(i),i=1,3)
 write(500,'(3f20.14)') (totc(i),i=1,3)
 write(500,'(a)')'%endblock lattice_cart'
 write(500,*)
 write(500,'(a)')'%block positions_abs'
 write(500,'(a)')'Ang'
 do i=1,natom
  write(500,'(a,3f17.8)') trim(adjustl(labelfull(i))),(coordfull(i,j),j=1,3)
 enddo
 write(500,'(a)')'%endblock positions_abs'
 write(500,*)
 write(500,'(a)')'%block species_pot'
 do i=1,nbasis  
 do j=1,natom
 if( trim(adjustl(atomsymbol(i))) == trim(adjustl(labelfull(j)))) goto 57
 enddo
 cycle
 57 continue
  write(500,*) trim(adjustl(atomsymbol(i)))//' "./pseudo/'//trim(adjustl(atomsymbol(i)))//'.LDA-PW-paw.abinit" '
 enddo

 !write(500,'(a,i0,a)')'Sr  "./pseudo/Sr.LDA-PW-paw.abinit"'
 !write(500,'(a,i0,a)')'Ti  "./pseudo/Ti.LDA-PW-paw.abinit"'
 !write(500,'(a,i0,a)')'O   "./pseudo/O.LDA-PW-paw.abinit "'
 !write(500,'(a,i0,a)')'La  "./pseudo/La.LDA-PW-paw.abinit"'
 !write(500,'(a,i0,a)')'Al  "./pseudo/Al.LDA-PW-paw.abinit"'

 write(500,'(a)')'%endblock species_pot'
 write(500,*)
 write(500,'(a)')'%block species_core_wf'
 do i=1,nbasis
 do j=1,natom
 if( trim(adjustl(atomsymbol(i))) == trim(adjustl(labelfull(j)))) goto 58
 enddo
 cycle
 58 continue
  write(500,*) trim(adjustl(atomsymbol(i)))//' "./pseudo/'//trim(adjustl(atomsymbol(i)))//'.LDA-PW-corewf.abinit" '
 enddo

 !write(500,'(a,i0,a)')'Sr  "./pseudo/Sr.LDA-PW-corewf.abinit"'
 !write(500,'(a,i0,a)')'Ti  "./pseudo/Ti.LDA-PW-corewf.abinit"'
 !write(500,'(a,i0,a)')'O   "./pseudo/O.LDA-PW-corewf.abinit "'
 !write(500,'(a,i0,a)')'La  "./pseudo/La.LDA-PW-corewf.abinit"'
 !write(500,'(a,i0,a)')'Al  "./pseudo/Al.LDA-PW-corewf.abinit"'
 write(500,'(a)')'%endblock species_core_wf'
 write(500,*)

 write(500,'(a)')'%block species_atomic_set'
 do i=1,nbasis
 do j=1,natom
 if( trim(adjustl(atomsymbol(i))) == trim(adjustl(labelfull(j)))) goto 59
 enddo
 cycle
 59 continue
  write(500,*) trim(adjustl(atomsymbol(i)))//' "SOLVE conf='//trim(adjustl(atomConf(i))//'"')
 enddo
 !write(500,'(a,i0,a)')'Sr  "SOLVE conf=4dX 4fX 5s1.2"'
 !write(500,'(a,i0,a)')'Ti  "SOLVE conf=4s1 3d2"'
 !write(500,'(a,i0,a)')'O   "SOLVE conf=2s2 2p4.6 3sX 3pX"'
 !write(500,'(a,i0,a)')'La  "SOLVE conf=4f0 5s2 5p6 5d1 5fX 5gX 6s1 6p0"'
 !write(500,'(a,i0,a)')'Al  "SOLVE conf=3s1.2"'
 write(500,'(a)')'%endblock species_atomic_set'
 write(500,*)
 close(500)

 end subroutine

END PROGRAM

!**************************************!
!**************************************!
!**************************************!
!**************************************!
!**************************************!
!**************************************!
