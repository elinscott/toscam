!  
!  E l e c t r o n i c   S t r u c t u r e   D a t a   F o r m a t
!  ---------------------------------------------------------------
!  
!                             E S D F
!                             =======
!  
!  Author: Chris J. Pickard (c)
!  Email : cp@min.uni-kiel.de
!  Place : Kiel, Germany
!  Date  : 5/6th August 1999
!  
!  Summary
!  -------
!  
!  This module is designed to simplify and enhance the input of data into
!  electronic structure codes (for example, CASTEP). It works from a
!  highly flexible input file. data is input in a "label <value>"
!  fashion, and input is independent of the ordering of the input
!  file. An important feature is the requirement that most inputs require
!  default settings to be supplied within the main program calling
!  ESDF. This means that rarely used variables will not clutter everyday
!  input files, and, even more usefully, "intelligence" may be built into
!  the main code as the defaults may be dependent of other set
!  variables. Block data may also be read in. Another important feature
!  is the ability to define "physical" values. This means that the input
!  files need not depend on the internal physical units used by the main
!  program.
!  
!  
!  History
!  -------
!  
!  ESDF has been written from scratch in F90, but is heavily based
!  (especially for the concept) on the FDF package developed by Alberto
!  Garcia and Jose Soler. It is not as "flexible" as FDF - there is no
!  provision for branching to further input files. This simplifies the
!  code, and I hope that it is still useful without this feature. Also,
!  the input and defaults are not dumped to a output file currently. I've
!  not found this a hindrance as of now.
!  
!  
!  Future
!  ------ 
!  
!  My intention is to make this release available to Alberto Garcia and
!  Jose Soler for their comments. It might be a good idea to use this as
!  a base for fully converting the FDF package to F90. Or it may remain
!  as a cut down version of FDF. I certainly hope that a package of the
!  FDF sort becomes widely used in the electronic structure community. My
!  experience has been very positive.
!  
!  
!  Usage
!  -----
!  
!  First, "Use esdf" wherever you wish to make use of its features. In
!  the main program call the initialisation routine: call
!  esdf_init('input.esdf'). "input.esdf" is the name of the input file -
!  it could be anything. This routine opens the input file, and reads
!  into a dynamically allocated storage array. The comments and blank
!  lines are stripped out. You are now ready to use the
!  esdf_functions. For example, if you want to read in the number of
!  atoms in your calculation, you would use: natom =
!  esdf_integer('NumberOfAtoms',1), where 'NumberOfAtoms' is the label to
!  search for in the input file, and '1' is the default. call esdf_close to
!  deallocate the data arrays. You may then open another input file using
!  esdf_init. It is not currently possible to open more that on input 
!  file at one time.
!  
!  
!  Syntax
!  ------
!  
!  The input file can contain comments. These are defined as anything to
!  the right of, and including, '#', ';', or '!'. It is straightforward
!  to modify the routine to accept further characters. Blank lines are
!  ignored -- use comments and blank lines to make you input file
!  readable.
!  
!  The "labels" are case insensitive (e.g. unitCell is equivalent to
!  UnItceLL) and punctuation insensitive (unit.cell is equivalent to
!  unit_cell is equivalent to unitcell). Punctuation characters are '.',
!  '_', and '-' at the moment. Again - use this feature to improve
!  readability.
!  
!  The following are equivalent ways of defining a physical quantity:
!  
!  "AgeOfUniverse = 24.d0 s" or "AgeOfUniverse : 24.d0 S" or
!  "AgeOfUniverse 24.d0 S"
!   
!  It would be read in by the main program in the following way:
!  
!  aou = esdf_physical('ageofuniverse',77.d0,ns)
!  
!  "aou" is the double precision variable, 77.d0 is the default number of
!  "ns" or nanoseconds. 24s will be converted automatically to its
!  equivalent number of nanoseconds.
!  
!  Block data should be placed in the input file as follows:
!  
!  %block cellvectors 
!  1.0 1.0 0.0 
!  0.0 1.0 1.0 
!  1.0 0.0 1.0 
!  %endblock cellvectors
!  
!  And it may be read:
!  
!    if (esdf_block('CellVectors',nlines))
!      if (nlines.ne.3) then (... break out here if the incorrect number
!  of lines)
!      do i=1,nlines
!        read(block_data(i),*) x,y,z
!      end do
!    end if
!  
!  
!  List of functions
!  -----------------
!  
!  Self explanatory:
!  
!  esdf_string(label,default)
!  esdf_integer(label,default)
!  esdf_single(label,default)
!  esdf_double(label,default)
!  esdf_physical(label,default,unit)
!  
!  A little more explanation:
!  
!  esdf_defined(label) is true if "label" found, false otherwise
!  
!  esdf_boolean(label,default) is true if "label yes/true/t (case/punct.insens)
!                              is false if"label no/false/f (case/punct.insens)
!  
!  The Help feature
!  ----------------
!  
!  The routine "esdf_help(helpword,searchword)" can be used to access the
!  information contained within the "esdf_key_mod" module.
!  
!  if "helpword" is "search" (case insensitive), then all variables whose
!  description contains "searchword" will be output.
!  
!  if "helpword" is "basic", "inter", "expert" or "dummy" the varibles of
!  that type will be displayed.
!  
!  if "helpword" is one of the valid labels, then a description of this
!  label will be output.
!  
!  
!  Finishing off
!  -------------
!  
!  Three routines, "esdf_dump", "esdf_warnout" and "esdf_close", can be
!  used to finish the use of ESDF. "esdf_dump" outputs a file "ESDF.esdf"
!  which could be used as an input file for further runs. "esdf_warnout"
!  outputs ESDF warnings to screen, and "esdf_close" deallocates the
!  allocated ESDF arrays.
!  
!  Contact the Author
!  ------------------
!  
!  This code is under development, and the author would be very happy to
!  receive comments by email. Any use in a commercial software package is 
!  forbidden without prior arrangement with the author (Chris J. Pickard).
!
!
!  Code cleanup by Nicholas Hine, 21 October 2011


Module esdf

  Use esdf_key
  Use constants, only: DP

  Implicit None

  ! Kind parameters 

  !qoh: Remove non standard integer kind
  !integer, private, parameter :: I4B = selected_int_kind(9)
  !qoh: Use ONETEP's standard DP from constants
  !integer, private, parameter :: DP  = Kind(1.d0)
  integer, private, parameter :: SP  = Kind(1.0)

  ! Set the length of the lines

  integer, public, parameter :: llength=80
  !cks 18/4/2000 the size of ndump has been changed to 20000 from its original value of 200
  !ndmh 08/09/2008 the size of ndump has been further increased to 200000
  integer, private, parameter ::  nphys = 70, ndump = 200000
  integer, private :: nrecords,nwarns,ndmp
  character(llength), private, dimension(:), allocatable :: llist,warns,dump 
  character(llength), private, dimension(:,:), allocatable :: tlist 

  ! The public block data array
  character(llength), public, dimension(:), allocatable :: block_data 

  ! Set the physical units database
  type phys_unit
     character(10) :: d,n ! d - dimension n - name
     real(kind=DP) :: u   ! u - unit
  end type phys_unit

  type(phys_unit), private, dimension(nphys) :: phy

  !
  !     We allow case variations in the units. This could be dangerous
  !     (meV --> MeV !!) in real life, but not in this restricted 
  !     field.
  !
  ! m - mass l - length t - time e - energy f - force p - pressure c- charge
  ! d - dipole mom - mom inert ef - efield s - surface tension
  !
  ! jd: 2010/06/08: Added a new category: surface tension, with two
  !                 units: ha/bohr**2 and J/m**2
  data phy(1)%d /'m'/;data phy(1)%n /'kg'/;data phy(1)%u /1.d0/
  data phy(2)%d /'m'/;data phy(2)%n /'g'/;data phy(2)%u /1.d-3/
  data phy(3)%d /'m'/;data phy(3)%n /'amu'/;data phy(3)%u /1.66054d-27/
  data phy(4)%d /'l'/;data phy(4)%n /'m'/;data phy(4)%u /1.d0/
  data phy(5)%d /'l'/;data phy(5)%n /'nm'/;data phy(5)%u /1.d-9/
  data phy(6)%d /'l'/;data phy(6)%n /'ang'/;data phy(6)%u /1.d-10/
  data phy(7)%d /'l'/;data phy(7)%n /'bohr'/;data phy(7)%u /0.529177d-10/
  data phy(8)%d /'t'/;data phy(8)%n /'aut'/;data phy(8)%u /2.41888468d-17/
  data phy(9)%d /'t'/;data phy(9)%n /'s'/;data phy(9)%u /1.d0/
  data phy(10)%d /'t'/;data phy(10)%n /'ns'/;data phy(10)%u /1.d-9/
  data phy(11)%d /'t'/;data phy(11)%n /'ps'/;data phy(11)%u /1.d-12/
  data phy(12)%d /'t'/;data phy(12)%n /'fs'/;data phy(12)%u /1.d-15/
  data phy(13)%d /'e'/;data phy(13)%n /'j'/;data phy(13)%u /1.d0/
  data phy(14)%d /'e'/;data phy(14)%n /'erg'/;data phy(14)%u /1.d-7/
  data phy(15)%d /'e'/;data phy(15)%n /'ev'/;data phy(15)%u /1.60219d-19/
  data phy(16)%d /'e'/;data phy(16)%n /'mev'/;data phy(16)%u /1.60219d-22/
  data phy(17)%d /'e'/;data phy(17)%n /'ry'/;data phy(17)%u /2.17991d-18/
  data phy(18)%d /'e'/;data phy(18)%n /'mry'/;data phy(18)%u /2.17991d-21/
  data phy(19)%d /'e'/;data phy(19)%n /'hartree'/;data phy(19)%u /4.35982d-18/
  data phy(20)%d /'e'/;data phy(20)%n /'kcal/mol'/;data phy(20)%u /6.94780d-21/
  data phy(21)%d /'e'/;data phy(21)%n /'mhartree'/;data phy(21)%u /4.35982d-21/
  data phy(22)%d /'e'/;data phy(22)%n /'kj/mol'/;data phy(22)%u /1.6606d-21/
  data phy(23)%d /'e'/;data phy(23)%n /'hz'/;data phy(23)%u /6.6262d-34/
  data phy(24)%d /'e'/;data phy(24)%n /'thz'/;data phy(24)%u /6.6262d-22/
  data phy(25)%d /'e'/;data phy(25)%n /'cm-1'/;data phy(25)%u /1.986d-23/
  data phy(26)%d /'e'/;data phy(26)%n /'cm^-1'/;data phy(26)%u /1.986d-23/
  data phy(27)%d /'e'/;data phy(27)%n /'cm**-1'/;data phy(27)%u /1.986d-23/
  data phy(28)%d /'f'/;data phy(28)%n /'N'/;data phy(28)%u /1.d0/
  data phy(29)%d /'f'/;data phy(29)%n /'ev/ang'/;data phy(29)%u /1.60219d-9/
  data phy(30)%d /'f'/;data phy(30)%n /'ry/bohr'/;data phy(30)%u /4.11943d-8/
  data phy(31)%d /'l'/;data phy(31)%n /'cm'/;data phy(31)%u /1.d-2/
  data phy(32)%d /'p'/;data phy(32)%n /'pa'/;data phy(32)%u /1.d0/
  data phy(33)%d /'p'/;data phy(33)%n /'mpa'/;data phy(33)%u /1.d6/
  data phy(34)%d /'p'/;data phy(34)%n /'gpa'/;data phy(34)%u /1.d9/
  data phy(35)%d /'p'/;data phy(35)%n /'atm'/;data phy(35)%u /1.01325d5/
  data phy(36)%d /'p'/;data phy(36)%n /'bar'/;data phy(36)%u /1.d5/
  data phy(37)%d /'p'/;data phy(37)%n /'mbar'/;data phy(37)%u /1.d11/
  data phy(38)%d /'p'/;data phy(38)%n /'ry/bohr**3'/;data phy(38)%u /1.47108d13/
  data phy(39)%d /'p'/;data phy(39)%n /'ev/ang**3'/;data phy(39)%u /1.60219d11/
  data phy(40)%d /'c'/;data phy(40)%n /'c'/;data phy(40)%u /1.d0/
  data phy(41)%d /'c'/;data phy(41)%n /'e'/;data phy(41)%u /1.602177d-19/
  data phy(42)%d /'d'/;data phy(42)%n /'C*m'/;data phy(42)%u /1.d0/
  data phy(43)%d /'d'/;data phy(43)%n /'D'/;data phy(43)%u /3.33564d-30/
  data phy(44)%d /'d'/;data phy(44)%n /'debye'/;data phy(44)%u /3.33564d-30/
  data phy(45)%d /'d'/;data phy(45)%n /'e*bohr'/;data phy(45)%u /8.47835d-30/
  data phy(46)%d /'d'/;data phy(46)%n /'e*ang'/;data phy(46)%u /1.602177d-29/
  data phy(47)%d /'mom'/;data phy(47)%n /'kg*m**2'/;data phy(47)%u /1.d0/
  data phy(48)%d /'mom'/;data phy(48)%n /'ry*fs**2'/;data phy(48)%u /2.1799d-48/
  data phy(49)%d /'ef'/;data phy(49)%n /'v/m'/;data phy(49)%u /1.d0/
  data phy(50)%d /'ef'/;data phy(50)%n /'v/nm'/;data phy(50)%u /1.d9/
  data phy(51)%d /'ef'/;data phy(51)%n /'v/ang'/;data phy(51)%u /1.d10/
  data phy(52)%d /'ef'/;data phy(52)%n /'v/bohr'/;data phy(52)%u /1.8897268d10/
  data phy(53)%d /'ef'/;data phy(53)%n /'ry/bohr/e'/;data phy(53)%u /2.5711273d11/
  data phy(54)%d /'ef'/;data phy(54)%n /'har/bohr/e'/;data phy(54)%u /5.1422546d11/
  data phy(55)%d /'e'/;data phy(55)%n /'k'/;data phy(55)%u /1.38066d-23/
  data phy(56)%d /'f'/;data phy(56)%n /'ha/bohr'/;data phy(56)%u /8.23886d-8/
  data phy(57)%d /'p'/;data phy(57)%n /'ha/bohr**3'/;data phy(57)%u /2.94216d13/
  data phy(58)%d /'1/l'/;data phy(58)%n /'1/ang'/;data phy(58)%u /1.d10/
  data phy(59)%d /'1/l'/;data phy(59)%n /'1/bohr'/;data phy(59)%u /1.889727d10/
  data phy(60)%d /'1/l'/;data phy(60)%n /'1/nm'/;data phy(60)%u /1.d9/
  data phy(61)%d /'s'/;data phy(61)%n /'ha/bohr**2'/;data phy(61)%u /1.556894d3/
  data phy(62)%d /'s'/;data phy(62)%n /'J/m**2'/;data phy(62)%u /1.d0/
  data phy(63)%d /'v'/;data phy(63)%n /'m/s'/;data phy(63)%u /1.d0/
  data phy(64)%d /'v'/;data phy(64)%n /'auv'/;data phy(64)%u /2.1876912633d6/
  data phy(65)%d /'v'/;data phy(65)%n /'ang/ns'/;data phy(65)%u /0.1d0/
  data phy(66)%d /'v'/;data phy(66)%n /'ang/ps'/;data phy(66)%u /1.d2/
  data phy(67)%d /'v'/;data phy(67)%n /'ang/fs'/;data phy(67)%u /1.d5/
  data phy(68)%d /'v'/;data phy(68)%n /'bohr/ns'/;data phy(68)%u /0.529177d-1/
  data phy(69)%d /'v'/;data phy(69)%n /'bohr/ps'/;data phy(69)%u /0.529177d2/
  data phy(70)%d /'v'/;data phy(70)%n /'bohr/fs'/;data phy(70)%u /0.529177d5/

contains

  subroutine esdf_init(filename,error)

    use utils, only: utils_alloc_check
    implicit none

    character(*), intent(in) :: filename
    integer, intent(out) :: error

    ! Local

    integer, parameter :: ncomm=3,ndiv=3
    integer :: unit,ierr,i,j,ic,nt,ndef
    character(llength) :: cjunk,ctemp,ctemp2
    character(1) :: comment(ncomm),divide(ndiv)
    logical :: inblock
    
    ! smmd : local variable to handle inclusion of external files
    integer :: ixf, xfile_num_loc,xfile_num_tot, xfile_unit(llength)
    character(llength) :: xfile_name(llength)

    ! Define comment characters
    data comment /'#',';','!'/
    data divide /' ','=',':'/

    ! Set error flag to success by default
    error = 0

    ! "reduce" the keyword list for comparison
    do i = 1,numkw
       ctemp = kw(i)%label
       ! qoh: Explicit character trunctation to avoid compiler warning
       ctemp2 = esdf_reduce(ctemp)
       kw(i)%label = ctemp2(1:30) 
    end do

    ! open the esdf file
    call esdf_file(unit,filename,ierr)
    cjunk = 'Unable to open main input file "'//trim(filename)//'"'

    if (ierr.eq.1) then
       write(6,*) 'ESDF WARNING: '//trim(cjunk)//' - using defaults'
       nrecords = 1
    else

       ! Count the number of records (excluding blank lines, commented lines, and included lines)
       nrecords = 0
       xfile_num_tot = 0

       do
          read(unit,'(a)',end=100,err=102,iostat=ierr) cjunk
          do j=1,ncomm
             ic=index(cjunk,comment(j))
             if (ic.gt.0) cjunk(ic:) = ' '
          end do

          ! smmd: Check for inclusion of external files via the "includefile" keyword
          ixf=index(esdf_reduce(cjunk),'includefile')
          if (len_trim(cjunk).gt.0 .and. ixf.eq.0) then
             nrecords = nrecords + 1
          end if

          ! smmd:  Count the number of records coming from external files
          if (ixf.gt.0) then
             xfile_num_loc=0
             ctemp = adjustl(cjunk)
             ixf = minval(index(ctemp,divide),mask=index(ctemp,divide)>0)
             ctemp = adjustl(ctemp(ixf:)) 
             do while(len_trim(ctemp).gt.0)
                ixf = minval(index(ctemp,divide),mask=index(ctemp,divide)>0)
                if (ixf.gt.1) then
                   xfile_num_tot=xfile_num_tot+1
                   xfile_num_loc=xfile_num_loc+1
                   xfile_name(xfile_num_tot) = adjustl(ctemp(:ixf-1))
                end if
                ctemp = adjustl(ctemp(ixf+1:))
             end do

             do ixf=xfile_num_tot-xfile_num_loc+1,xfile_num_tot
                call esdf_file(xfile_unit(ixf),xfile_name(ixf),ierr)
                cjunk = 'Unable to open secondary input file "'//trim(xfile_name(ixf))//'"'
                if (ierr.eq.1) then
                   write(6,*) 'ESDF WARNING: '//trim(cjunk)//' - using defaults'
                else
                   do
                      read(xfile_unit(ixf),'(a)',end=112,err=102,iostat=ierr) cjunk
                      do j=1,ncomm
                         ic=index(cjunk,comment(j))
                         if (ic.gt.0) cjunk(ic:) = ' '
                      end do
                      if (len_trim(cjunk).gt.0) then
                         nrecords = nrecords + 1
                      end if
                   end do
112                rewind(xfile_unit(ixf))
                end if
             end do
          end if 

       end do
100    rewind(unit)

    end if

    ! allocate the array to hold the records and tokens
    allocate(llist(nrecords),stat=ierr)
    call utils_alloc_check('esdf_init','llist',ierr)
    allocate(block_data(nrecords),stat=ierr)
    call utils_alloc_check('esdf_init','block_data',ierr)
    allocate(tlist(llength,nrecords),stat=ierr)
    call utils_alloc_check('esdf_init','tlist',ierr)
    allocate(warns(nrecords),stat=ierr)
    call utils_alloc_check('esdf_init','warns',ierr)
    allocate(dump(ndump),stat=ierr)
    call utils_alloc_check('esdf_init','dump',ierr)

    ! Set the number of warnings to zero
    nwarns = 0 ; warns = ' ' ; ndmp = 0 ; dump = ' '

    if (ierr.eq.1) then
       nrecords = 0
    else

       ! Read in the records
       nrecords = 0
       xfile_num_tot = 0
       do
          read(unit,'(a)',end=101,err=103,iostat=ierr) cjunk
          do j=1,ncomm
             ic=index(cjunk,comment(j))
             if (ic.gt.0) cjunk(ic:) = ' '
          end do

          ixf=index(esdf_reduce(cjunk),'includefile')

          if (len_trim(cjunk).gt.0 .and. ixf.eq.0) then
             nrecords=nrecords+1
             llist(nrecords) = adjustl(cjunk)
          end if

          ! smmd:  Read in the records coming from external files
          if (ixf.gt.0) then
             xfile_num_loc=0
             ctemp = adjustl(cjunk)
             ixf = minval(index(ctemp,divide),mask=index(ctemp,divide)>0)
             ctemp = adjustl(ctemp(ixf:)) 
             do while(len_trim(ctemp).gt.0)
                ixf = minval(index(ctemp,divide),mask=index(ctemp,divide)>0)
                if (ixf.gt.1) then
                   xfile_num_loc=xfile_num_loc+1
                   xfile_num_tot=xfile_num_tot+1
                end if
                ctemp = adjustl(ctemp(ixf+1:))
             end do
             do ixf=xfile_num_tot-xfile_num_loc+1,xfile_num_tot
                do
                   read(xfile_unit(ixf),'(a)',end=113,err=103,iostat=ierr) cjunk
                   do j=1,ncomm
                      ic=index(cjunk,comment(j))
                      if (ic.gt.0) cjunk(ic:) = ' '
                   end do
                   
                   if (len_trim(cjunk).gt.0) then
                      nrecords=nrecords+1
                      llist(nrecords) = adjustl(cjunk)
                   end if
                end do
113             close(xfile_unit(ixf))
             end do
          end if 

       end do
101    close(unit)

    end if

    ! Now read in the tokens from llist
    tlist = ' '

    do i=1,nrecords
       ctemp = llist(i)
       nt=0
       do while(len_trim(ctemp).gt.0)
          ic = minval(index(ctemp,divide),mask=index(ctemp,divide)>0)
          if (ic.gt.1) then
             nt=nt+1
             tlist(nt,i) = adjustl(ctemp(:ic-1))
          end if
          ctemp = adjustl(ctemp(ic+1:))
       end do
    end do

    ! Check if any of the "labels" in the input file are unrecognised
    inblock=.false.
    do i=1,nrecords
       ! Check if we are in a block
       if (esdf_reduce(tlist(1,i)).eq.'%block') then
          inblock = .true.
          ! Check if block label is recognised
! cks hack          if ((count(esdf_reduce(tlist(2,i)).eq.kw%label).eq.0)) then
!             ctemp='Label "'//trim(esdf_reduce(tlist(2,i)))//&
!                  &'" not in keyword list'
!             if (count(ctemp.eq.warns).eq.0) call esdf_warn(ctemp) 
!          end if
          ! Check if "label" is multiply defined in the input file
          ndef=0
          do j=1,nrecords
             if (esdf_reduce(tlist(2,i)).eq.esdf_reduce(tlist(2,j))) ndef=ndef+1
          end do
          ctemp='Label "'//trim(esdf_reduce(tlist(2,i)))//&
               &'" is multiply defined in the input file. '
          if ((ndef.gt.2).and.(count(ctemp.eq.warns).eq.0))&
               call esdf_warn(ctemp)
       end if
       ! Check it is in the list of keywords
!       if ((count(esdf_reduce(tlist(1,i)).eq.kw%label).eq.0)&
!            .and.(.not.inblock)) then
!          ctemp='Label "'//trim(esdf_reduce(tlist(1,i)))//&
!               &'" not in keyword list'
!          if (count(ctemp.eq.warns).eq.0) call esdf_warn(ctemp)
!       end if
       if (.not.inblock) then
          ! Check if "label" is multiply defined in the input file
          ndef=0
          do j=1,nrecords
             if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(tlist(1,j))) ndef=ndef+1
          end do
          ctemp='Label "'//trim(esdf_reduce(tlist(1,i)))//&
               &'" is multiply defined in the input file. '
          if ((ndef.gt.1).and.(count(ctemp.eq.warns).eq.0)) &
               call esdf_warn(ctemp)
       end if
       ! Check if we have left a block
       if (esdf_reduce(tlist(1,i)).eq.'%endblock') inblock= .false.

    end do

    return

    ! Error on read while counting records
102 write(6,'(3a,i6)') 'Error in esdf_init: counting records in file "', &
         trim(filename),'" failed with code ',ierr
    error = 1
    return

    ! Error on read while reading records
103 write(6,'(3a,i6)') 'Error in esdf_init: reading records in file "', &
         trim(filename),'" failed with code ',ierr
    error = 1
    return

  end subroutine esdf_init

  !  
  ! return the string attached to the "label"
  !
  function esdf_string(label,default)

    character(*), intent(in) :: label,default
    character(llength) :: esdf_string

    ! Local
    integer :: i
    character(llength) :: ctemp
    ! aam: For conversion of string to upper case
    integer   :: ipos      ! Loop variable over characters in string
    character :: letter    ! character in string

    ! Check "label" is defined
    call esdf_lblchk(label,'T')

    ! Set to default
    esdf_string = default

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          esdf_string = llist(i)(index(llist(i),trim(tlist(2,i))):)
          exit
       end if

    end do

    ! aam: Convert string to uppercase
    ! aam: Loop over characters in string
    do ipos=1,len(esdf_string)

       letter = esdf_string(ipos:ipos)
       if (letter >= 'a' .and. letter <= 'z') &
            esdf_string(ipos:ipos) = achar(iachar(letter)-32)

    end do

    ! Dump the string used 
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),':',trim(esdf_string)
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    return

101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_string'
    call esdf_warn(ctemp)
  end function esdf_string

  !  
  ! return the integer attached to the "label"
  !
  function esdf_integer(label,default)

    integer, intent(in) :: default
    character(*), intent(in) :: label
    integer :: esdf_integer

    ! Local
    integer :: i
    character(llength) :: ctemp

    ! Check "label" is defined
    call esdf_lblchk(label,'I')

    ! Set to default
    esdf_integer = default

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          read(tlist(2,i),*,err=100) esdf_integer
          exit
       end if

    end do

    ! Dump the value used 

    ndmp=ndmp+1
    write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),':',esdf_integer
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    return


100 ctemp = 'Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_integer'
    call esdf_die(ctemp)

101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_integer'
    call esdf_warn(ctemp)

  end function esdf_integer

  !  
  ! return the single precisioned value attached to the "label"
  !
  function esdf_single(label,default)

    real(SP), intent(in) :: default
    character(*), intent(in) :: label
    real(SP) :: esdf_single

    ! Local
    integer :: i
    character(llength) :: ctemp

    ! Check "label" is defined
    call esdf_lblchk(label,'S')

    ! Set to default
    esdf_single = default

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          read(tlist(2,i),*,err=100) esdf_single
          exit
       end if

    end do

    ! Dump the value used 
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),':',esdf_single
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    return

100 ctemp = 'Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_single'
    call esdf_die(ctemp)
101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_single'
    call esdf_warn(ctemp)

  end function esdf_single

  !  
  ! return the double precisioned value attached to the "label"
  !
  function esdf_double(label,default)

    real(kind=DP), intent(in) :: default
    character(*), intent(in) :: label
    real(kind=DP) :: esdf_double

    ! Local
    integer :: i
    character(llength) :: ctemp

    ! Check "label" is defined
    call esdf_lblchk(label,'D')
    
    ! Set to default
    esdf_double = default

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          read(tlist(2,i),*,err=100) esdf_double
          exit
       end if

    end do  
 
    ! Dump the value used 
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),':',esdf_double
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    return

100 esdf_double = default
    ctemp = 'Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_double'
    call esdf_die(ctemp)
101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_double'
    call esdf_warn(ctemp)

  end function esdf_double

  !  
  ! return the double precisioned physical value attached to the "label"
  ! Units converted to "dunit"
  !
  function esdf_physical(label,default,dunit)

    real(kind=DP), intent(in) :: default
    character(*), intent(in) :: label,dunit
    real(kind=DP) :: esdf_physical

    ! Local
    integer :: i
    character(llength) :: ctemp,iunit

    ! Check "label" is defined
    call esdf_lblchk(label,'P')

    ! Set to default
    esdf_physical = default

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          read(tlist(2,i),*,err=100) esdf_physical
          ! jd: Changed to '(a)', or otherwise parsing stopped at '/',
          !     so 'ha/bohr' would be interpreted as merely 'ha'.
          read(tlist(3,i),'(a)',err=100) iunit
          esdf_physical = esdf_convfac(iunit,dunit)*esdf_physical
          exit
       end if

    end do

    ! Dump the value used 
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),':',esdf_physical,' ',trim(dunit)
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    return

100 esdf_physical = default
    ctemp = 'Unable to parse "'//trim(esdf_reduce(label))//'" in esdf_physical'
    call esdf_die(ctemp)
101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_physical'
    call esdf_warn(ctemp)

  end function esdf_physical

  !
  ! Is the "label" defined in the input file
  !
  function esdf_defined(label)

    character(*), intent(in) :: label
    logical :: esdf_defined

    ! Local
    integer :: i
    character(llength) :: ctemp

    ! Check "label" is defined
    call esdf_lblchk(label,'E')

    ! Set to default
    esdf_defined = .false.

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          esdf_defined = .true.
          exit
       end if

    end do

    ! Dump the value used 
    if (esdf_defined) then

       ndmp=ndmp+1
       write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),':'
       if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    end if

    return

101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_defined'
    call esdf_warn(ctemp)

  end function esdf_defined

  !
  ! Is the "label" defined in the input file
  !
  function esdf_boolean(label,default)

    character(*), intent(in) :: label
    logical, intent(in) :: default
    logical :: esdf_boolean

    ! Local
    integer :: int_buffer(3)
    integer :: i
    character(llength) :: ctemp,positive(3),negative(3)

    data positive /'yes','true','t'/
    data negative /'no','false','f'/

    ! Check "label" is defined
    call esdf_lblchk(label,'L')

    ! Set to default
    esdf_boolean = default

    do i=1,nrecords

       ! Search in the first token for "label"
       ! The first instance is returned
       if (esdf_reduce(tlist(1,i)).eq.esdf_reduce(label)) then
          if (len_trim(tlist(2,i)).eq.0) then
             esdf_boolean = .true.
             exit
          end if

          int_buffer =index(positive,esdf_reduce(tlist(2,i)) )

          if ( Sum( int_buffer ) .gt.0 ) then

             esdf_boolean = .true.
             exit
          end if

          int_buffer =index(negative,esdf_reduce(tlist(2,i)))

          if (  Sum( int_buffer ) .gt.0 ) then

             esdf_boolean = .false.
             exit
          end if
          call esdf_die('Unable to parse boolean value')

       end if

    end do

    ! Dump the value used 
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101)  trim(esdf_reduce(label)),': ',esdf_boolean
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) ndmp=ndmp-1

    return

101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_boolean'
    call esdf_warn(ctemp)

  end function esdf_boolean

  function esdf_block(label,nlines)

    character(*), intent(in) :: label
    integer, intent(out) :: nlines
    logical :: esdf_block

    ! Local
    integer :: i,j
    character(llength) :: ctemp

    ! Check "label" is defined
    call esdf_lblchk(label,'B')

    ctemp ='Block "'//trim(esdf_reduce(label))//'" not closed correctly '
    
    esdf_block=.false.

    nlines = 0

    do i=1,nrecords
       if ((esdf_reduce(tlist(1,i)).eq.esdf_reduce('%block'))&
            .and.(esdf_reduce(tlist(2,i)).eq.esdf_reduce(label))) then
          esdf_block = .true.
          do while(esdf_reduce(tlist(1,i+nlines+1))&
               .Ne.esdf_reduce('%endblock'))
             nlines=nlines+1
             if (nlines+i.gt.nrecords) call esdf_die(ctemp)
             block_data(nlines)=llist(i+nlines)
          end do

          if (esdf_reduce(tlist(2,i+nlines+1)).Ne.esdf_reduce(label))&
               call esdf_die(ctemp)
          exit
       end if
    end do

    if (.not.esdf_block) return

    ! Dump the block 
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101) '%block ',trim(esdf_reduce(label)),': '
    if (count(dump(ndmp).eq.dump(1:ndmp-1)).gt.0) then
       ndmp=ndmp-1
       return
    end if
    do j=1,nlines
       ndmp=ndmp+1
       dump(ndmp)=block_data(j)
    end do
    ndmp=ndmp+1
    write(dump(ndmp),*,err=101) '%endblock ',trim(esdf_reduce(label)),': '

    return

101 ctemp = 'Unable to dump "'//trim(esdf_reduce(label))//'" in esdf_block'
    call esdf_warn(ctemp)
    
  end function esdf_block

  !
  ! Reduce the string to lower case and remove punctuation
  !
  function esdf_reduce(string)

    character(*), intent(in) :: string
    character(llength) :: esdf_reduce

    ! Local

    integer, parameter :: npunct=5
    integer :: iA,iZ,ishift,ic,i,ln
    character(llength) :: ctemp
    character(1) :: punct(npunct)
    logical :: keep_going

    ! Define the punctuation to be removed
    data punct /'.','_','-','"',''''/

    ! Initialise system dependent bounds in collating sequence
    iA = ichar('A');iZ = ichar('Z') 
    ishift = ichar('a')-iA 

    ! Initialise output
    ln = len(string)

    esdf_reduce(1:ln) = string(1:ln)
    if (ln < llength) esdf_reduce(ln+1:)=' '

    ! Drop all upper case characters to lower case
    do i=1,llength
       ic = ichar(esdf_reduce(i:i))
       if ((ic.ge.iA).and.(ic.le.iZ)) esdf_reduce(i:i) = char(ishift+ic) 
    end do

    ! Now remove punctuation
    do i=1,npunct

       keep_going=.true.

       do while(keep_going)
          ic = index(esdf_reduce,punct(i))
          if (ic.gt.0) then
             ctemp = esdf_reduce
             esdf_reduce(ic:)=ctemp(ic+1:)//' '
          else
             keep_going=.false.
          end If
       end do

    end do

    esdf_reduce = trim(adjustl(esdf_reduce))

  end function esdf_reduce

  !
  ! Find the conversion factor between physical units
  !
  function esdf_convfac(from,to)

    character(*), intent(in) :: from,to
    real(kind=DP) :: esdf_convfac

    ! Local
    integer :: i,ifrom,ito
    character(llength) :: ctemp

    ! Find the index numbers of the from and to units
    ifrom = 0 ; ito = 0
    do i=1,nphys
       if (esdf_reduce(from).eq.phy(i)%n) ifrom = i
       if (esdf_reduce(to).eq.phy(i)%n) ito = i
    end do

    ! Check that the units were recognised
    if (ifrom.eq.0) then
       ctemp = 'Units not recognised in input file : '//trim(esdf_reduce(from))
       call esdf_die(ctemp)
    end if

    if (ito.eq.0) then
       ctemp = 'Units not recognised in Program : '//trim(esdf_reduce(to))
       call esdf_die(ctemp)
    end if

    ! Check that from and to are of the same dimensions
    if (phy(ifrom)%d.Ne.phy(ito)%d) then
       ctemp = 'dimensions do not match : '//trim(esdf_reduce(from))&
            & //' vs '//trim(esdf_reduce(to))
       call esdf_die(ctemp)
    end if

    ! Set the conversion factor
    esdf_convfac = phy(ifrom)%u/phy(ito)%u

  end function esdf_convfac

  ! 
  ! Find an unused i/o unit
  !
  function esdf_unit(ierr)
    integer, intent(out) :: ierr
    integer :: esdf_unit
    ! Local 
    logical :: op
    ierr=0
    do esdf_unit=10,99
       inquire(unit=esdf_unit,opened=op,err=100)
       if (.not.op) return
    end do
    call esdf_warn('Unable to find a free i/o unit using esdf_unit')
    ierr = 1
    return
100 call esdf_die('Error opening files by esdf_unit')
  end function esdf_unit

  !
  ! open an old file
  !
  subroutine esdf_file(unit,filename,ierr)
    character(*), intent(in) :: filename
    integer, intent(out) :: unit,ierr
    logical :: ex
    unit = esdf_unit(ierr)
    if (ierr.gt.0) return
    inquire(file=trim(filename),exist=ex,err=100)
    if (.not.ex) goto 100
    open(unit=unit,file=trim(filename),form='formatted',status='old',err=100)
    return
100 ierr=1
    return
  end subroutine esdf_file

  ! open a new file

  subroutine esdf_newfile(unit,filename,ierr)
    character(*), intent(in) :: filename
    integer, intent(out) :: unit,ierr
    unit = esdf_unit(ierr)
    if (ierr.gt.0) return
    open(unit=unit,file=trim(filename),form='formatted',status='replace',err=100)
    return
100 ierr=1
    return
  end subroutine esdf_newfile

  !
  ! Check that the label is known, and used correctly
  !
  subroutine esdf_lblchk(string,typ)
    character(*), intent(in) :: string
    character(1), intent(in) :: typ
    ! Local
    character(llength) :: ctemp
    character(1) :: tp
    integer :: i
    ! Check if label is recognised
    i=count(esdf_reduce(string).eq.kw%label)
    ctemp = 'Label "'//trim(esdf_reduce(string))//'" not recognised in&
         & keyword list'
    if (i.eq.0) call esdf_die(ctemp)
    ctemp = 'Label "'//trim(esdf_reduce(string))//'" is multiply defined'
    if (i.gt.1) call esdf_die(ctemp)
    ctemp = 'Label "'//trim(esdf_reduce(string))//'" has been used with the wrong type'
    tp = ' '
    i=0
    do while(tp.eq.' ')
       i=i+1
       if (esdf_reduce(string).eq.kw(i)%label) tp=kw(i)%typ(1:1)
    end do
    if (typ.Ne.tp) call esdf_die(ctemp)
  end subroutine esdf_lblchk


  subroutine esdf_help(helpword,searchword)
  
    use comms, only: comms_abort !qoh

    Implicit None

    character(*), intent(inout) :: helpword,searchword

    ! Local
    integer  :: i,indx,indx2,ln
    character(20) :: ctyp,clev
    character(60) :: title,fm
    character(80) :: ctemp
    character(1)  :: cl

    helpword = esdf_reduce(helpword)
    searchword = esdf_reduce(searchword)

    if (esdf_reduce(helpword).eq.'search') then

       if (len_trim(searchword).lt.1) call esdf_die('help: "searchword" is empty')

       ! Search for useful keywords
       do i=1,numkw
          if ((index(kw(i)%label,trim(searchword)).gt.0).Or.&
               (index(kw(i)%dscrpt,trim(searchword)).gt.0)) then 
             indx=index(kw(i)%dscrpt,'!*')-1
             if (indx.eq.-1) &
                  call esdf_die('help: keyword description incorrectly formatted')
             title = kw(i)%dscrpt(1:indx)
             ln=len(trim(title))
             if (ln.gt.80) call esdf_die('help: keyword title too long')

             write(6,*) kw(i)%label,trim(title)

          end If
       end do
       call comms_abort !qoh


    end if

    ! All keywords, short description
    if ('all'.eq.helpword) then
       do i=1,numkw
          if (len_trim(kw(i)%label).gt.0) then 
             indx=index(kw(i)%dscrpt,'!*')-1
             if (indx.eq.-1) &
                  call esdf_die('help: keyword description incorrectly formatted')
             title = kw(i)%dscrpt(1:indx)
             ln=len(trim(title))
             if (ln.gt.80) call esdf_die('help: keyword title too long')

             write(6,*) kw(i)%label,trim(title)

          end If
       end do
       call comms_abort !qoh
    end If

    ! All specific levels of keywords
    if (any((/'basic ','inter ','expert','dummy '/).eq.helpword)) then

       cl = ' ' ! qoh: Initialise to prevent compiler warning
       select case(helpword)
       case('basic')  ; cl = 'B'
       case('inter')  ; cl = 'I'
       case('expert') ; cl = 'E'
       case('dummy')  ; cl = 'D'
       end select

       do i=1,numkw
          if (kw(i)%typ(3:3).eq.cl) then 
             indx=index(kw(i)%dscrpt,'!*')-1
             if (indx.eq.-1) &
                  call esdf_die('help: keyword description incorrectly formatted')
             title = kw(i)%dscrpt(1:indx)
             ln=len(trim(title))
             if (ln.gt.80) call esdf_die('help: keyword title too long')

             write(6,*) kw(i)%label,trim(title)

          end If
       end do
       call comms_abort !qoh
    end If

    ! More information about a specific keyword
    if (.not.any(kw%label.eq.helpword)) &
         call esdf_die('help: keyword not recognised')
    if (count(kw%label.eq.helpword).gt.1) &
         call esdf_die('help: keyword entry duplicated')
    do i=1,numkw
       if (kw(i)%label.eq.helpword) then 
          indx=index(kw(i)%dscrpt,'!*')+1
          if (indx.eq.1) &
               call esdf_die('help: keyword description incorrectly formatted')
          title = kw(i)%dscrpt(1:indx)
          ln=len(trim(title))
          if (ln.gt.80) &
               call esdf_die('help: keyword title too long')
          if (ln.le.9) write(fm,'("(",i2,"x,a",i1,")")') 40-ln/2,ln
          if (ln.gt.9) write(fm,'("(",i2,"x,a",i2,")")') 40-ln/2,ln
          write(6,fm) trim(title)
          write(6,*)
          select case(kw(i)%typ(1:1))
          case('I') ; ctyp ='integer'
          case('S') ; ctyp ='Single Precision'
          case('D') ; ctyp ='Double Precision'
          case('P') ; ctyp ='Physical'
          case('T') ; ctyp ='String'
          case('E') ; ctyp ='Defined'
          case('B') ; ctyp ='Block'
          case('L') ; ctyp ='Boolean'   
          end select
          select case(kw(i)%typ(3:3))
          case('B') ; clev ='Basic'
          case('I') ; clev ='Intermediate'
          case('E') ; clev ='Expert'
          case('D') ; clev ='Dummy'
          end select
          write(fm,'(a,i2,a)') '("type: ",a,',&
               78-(7+len_trim(clev))-(6+len_trim(ctyp)),'x," Level: ",a)'
          write(ctemp,fm) trim(ctyp),trim(clev)
          write(6,'(a)') trim(ctemp)
          write(6,*)
          indx=indx+1
          ln = len(trim(kw(i)%dscrpt))
          do while (indx.lt.ln)
             ctemp = kw(i)%dscrpt(indx:Min(indx+80,ln))
             indx2=index(ctemp,' ',back=.true.)
             indx=indx+len(ctemp(:indx2))
          end do

       end If
    end do

    call comms_abort

  end subroutine esdf_help

  ! 
  ! Stop execution due to an error cause by esdf
  !
  subroutine esdf_die(string)
    use comms, only: comms_abort
    character(*), intent(in) :: string
    write(6,'(a,a)') ' ESDF ERROR: ',trim(string)
    write(6,'(a)') ' Stopping now'
    call comms_abort    
  end subroutine esdf_die

  ! 
  ! Warning due to an error cause by esdf
  !
  subroutine esdf_warn(string)
    character(*), intent(in) :: string
    nwarns=nwarns+1
    warns(nwarns) = string
  end subroutine esdf_warn

  !
  ! Dump the warnings to screen
  !
  subroutine esdf_warnout
    integer :: i

    do i=1,nwarns
       write(6,*)'i=',i
       write(6,*)'warns(i)=',warns(i)
       write(6,*) 'ESDF WARNING: '//trim(warns(i))
    end do

  end subroutine esdf_warnout

  !
  ! Deallocate the data arrays --- call this before re-initialising
  !
  subroutine esdf_close

    use utils, only: utils_dealloc_check
    implicit none

    integer :: ierr

    deallocate(llist,stat=ierr)
    call utils_dealloc_check('esdf_close','llist',ierr)
    deallocate(block_data,stat=ierr)
    call utils_dealloc_check('esdf_close','block_data',ierr)
    deallocate(tlist,stat=ierr)
    call utils_dealloc_check('esdf_close','tlist',ierr)
    deallocate(warns,stat=ierr)
    call utils_dealloc_check('esdf_close','warns',ierr)
    deallocate(dump,stat=ierr)
    call utils_dealloc_check('esdf_close','dump',ierr)

  end subroutine esdf_close

  !
  ! Dump an input file which contains all set variables
  ! including defaults
  !
  subroutine esdf_dump(filename)

    ! Local
    integer :: unit_dump,ierr,i,j,indx
    character(*), intent(in) :: filename
    character(llength) :: cjunk

    ! open the ESDF.esdf file
    call esdf_newfile(unit_dump,filename,ierr)
    if (ierr.eq.1) call esdf_die('Unable to open main input file "ESDF.esdf"')

    indx = maxval(index(dump(1:ndmp),':'))

    do i=1,ndmp
       j=index(dump(i),':')
       if (j.gt.0) then
          cjunk = dump(i)(1:j-1)
          dump(i)=trim(cjunk)//repeat(' ',indx-j+1)//': '&
               & //trim(adjustl(dump(i)(j+1:)))//'#'
       end if
    end do

    indx = maxval(index(dump(1:ndmp),'#',back=.true.))

    do i=1,ndmp
       j=index(dump(i),'#',back=.true.)
       if (j.gt.0) then

          dump(i)=dump(i)(1:j-1)//repeat(' ',indx-j+1)//'#'
       end if
    end do

    do i=1,ndmp
       j=index(dump(i),':')
       if (j.gt.0) then
          cjunk = dump(i)(1:j-1)
          do j=1,numkw
             if (index(cjunk,trim(kw(j)%label)).gt.0) exit
          end do
          select case(kw(j)%typ(1:1))
          case('I') ; cjunk ='integer'
          case('S') ; cjunk ='Single Precision'
          case('D') ; cjunk ='Double Precision'
          case('P') ; cjunk ='Physical'
          case('T') ; cjunk ='String'
          case('E') ; cjunk ='Defined'
          case('B') ; cjunk ='Block'
          case('L') ; cjunk ='Boolean'   
          end select
          indx=index(kw(j)%dscrpt,'!*')
          dump(i)=trim(dump(i))//trim(kw(j)%dscrpt(1:indx-1))//' ('//trim(adjustl(cjunk))//')'
       end if
    end do

    do i=1,ndmp
       write(unit_dump,'(a)') adjustl(dump(i))
    end do

  end subroutine esdf_dump


  ! smmd : 12/07/2011
  ! Split an input line into pieces according to divide symbols
  subroutine esdf_line_divide(snum,sline,line)

    use constants,   only: stdout
    use comms,       only: pub_on_root

    ! local
    integer, intent(inout) :: snum
    character(llength), intent(in) :: line
    character(llength), intent(inout) :: sline(snum)
    character(1) :: divide(4)
    
    data divide /' ','=',':','	'/

    integer :: ixf, ntmp
    character(llength) :: cjunk

    ntmp = 0
    cjunk = adjustl(line)

    split_loop : do while(len_trim(cjunk).gt.0 .and. ntmp.lt.snum)
       ixf = minval(index(cjunk,divide),mask=index(cjunk,divide)>0)
       if (ixf .gt. 1) then
          ntmp = ntmp + 1
          sline(ntmp) = adjustl(cjunk(:ixf-1))
       end if 
       cjunk = adjustl(cjunk(ixf+1:))
    end do split_loop

    snum = ntmp

  end subroutine esdf_line_divide


  ! cks : 11/7/2001
  ! Dump to the standard output an input file which contains all set variables
  ! including defaults
  subroutine esdf_stdout_dump( )

    ! Local Variables
    integer :: i,j,indx, bpos
    character(llength) :: cjunk, keyword
    character(128)     :: keyword_buffer1, keyword_buffer2, keyword_buffer3
    logical :: inblock

    ! Find maximum index at which a ':' separator occurs
    indx = 0
    inblock = .false.
    do i=1,ndmp
       if (.not.inblock.and.(index(esdf_reduce(dump(i)),'%block')>0)) &
            inblock =.true.
       if (inblock.and.(index(esdf_reduce(dump(i)),'%endblock')>0)) &
            inblock =.false.
       if (.not.inblock) indx = max(indx,index(dump(i),':'))
    end do

    ! Re-format the dump lines so that the ':' characters line up
    inblock = .false.
    do i=1,ndmp
       if (inblock) then
          j = 0
       else
          j = index(dump(i),':')
       end if
       if (.not.inblock.and.(index(esdf_reduce(dump(i)),'%block')>0)) &
            inblock =.true.
       if (inblock.and.(index(esdf_reduce(dump(i)),'%endblock')>0)) then
          j = index(dump(i),':')
          inblock =.false.
       end if
       if (j.gt.0) then
          cjunk = dump(i)(1:j-1)
          keyword_buffer1 = trim(adjustl(dump(i)(j+1:)))
          dump(i) = trim(cjunk)//repeat(' ',indx-j+1)//': '// keyword_buffer1
       end if
    end do

    ! Find last incidence of '#' character on any dump line
    indx = maxval(index(dump(1:ndmp),'#',back=.true.))

    ! Re-format the dump lines so that the '#' characters line up
    do i=1,ndmp
       j=index(dump(i),'#',back=.true.)
       if (j.gt.0) then
          dump(i)=dump(i)(1:j-1)//repeat(' ',indx-j+1)//'#'
       end if
    end do

    ! Now add the keyword descriptions for each variable
    do i=1,ndmp
       if (inblock) then
          j = 0
       else
          j = index(dump(i),':')
       end if
       if (.not.inblock.and.(index(esdf_reduce(dump(i)),'%block')>0)) &
            inblock =.true.
       if (inblock.and.(index(esdf_reduce(dump(i)),'%endblock')>0)) then
          j = index(dump(i),':')
          inblock =.false.
       end if
       if (j.gt.0) then
          cjunk = dump(i)(1:j-1)
          ! qoh,jd: handle "%block" and "%endblock" manually
          bpos = index(cjunk,"%block")
          if ( bpos > 0 .and. len(cjunk) > bpos+6 ) cjunk = dump(i)(bpos+6:j-1)
          bpos = index(cjunk,"%endblock")
          if ( bpos > 0 .and. len(cjunk) > bpos+9 ) cjunk = dump(i)(bpos+9:j-1)
          keyword = trim(adjustl(cjunk))

          do j=1,numkw
             if (keyword == trim(adjustl(kw(j)%label))) exit
          ! jd: Originally was:
          !  if (index(cjunk,trim(kw(j)%label)).gt.0) exit
          !     ... but this gave false matches if one keyword was
          !         a substring of another
          end do

          select case(kw(j)%typ(1:1))
          case('I') ; cjunk ='Integer'
          case('S') ; cjunk ='Single Precision'
          case('D') ; cjunk ='Double Precision'
          case('P') ; cjunk ='Physical'
          case('T') ; cjunk ='String'
          case('E') ; cjunk ='Defined'
          case('B') ; cjunk ='Block'
          case('L') ; cjunk ='Boolean'   
          end select
          indx=index(kw(j)%dscrpt,'!*')

          ! cks: HACK in order to make it run on Franklin
          keyword_buffer1 =trim(dump(i))
          keyword_buffer2 =trim(kw(j)%dscrpt(1:indx-1))
          keyword_buffer3 =trim(adjustl(cjunk))
          dump(i) = trim(keyword_buffer1)//' '//trim(keyword_buffer2)//' ('//trim(keyword_buffer3)//')'
       end if
    end do

    do i=1,ndmp
       write(6,'(a)') adjustl(dump(i))
    end do

  end subroutine esdf_stdout_dump




end Module esdf
 
